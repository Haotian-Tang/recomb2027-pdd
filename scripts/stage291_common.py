#!/usr/bin/env python3
"""Stage 291 shared library: certified screening on estimated distances.

New file. Nothing pre-existing is modified. The screening geometry is taken
unchanged from the frozen Stage 44 routine `path_screen`; this module only

  (a) reads every bank family, including the ones with no ground truth,
  (b) re-derives the two ENDPOINT-LOCAL excesses of each edge from the frozen
      routine's own message dictionary, and cross-checks that their minimum
      reproduces the frozen witness exactly on every edge of every case,
  (c) projects the tree onto a retained reference subset (hull with degree-two
      suppression) and keeps, for each hull edge, the ordered list of original
      edges inside it, so that an exclusion certificate on the reduced panel is
      turned into exclusion certificates for ORIGINAL edges.

Hypotheses of the screen (Stage 44), restated and used verbatim here:
  H1  the reference topology and all edge lengths are known exactly;
  H2  every source lies strictly inside an ORIGINAL edge, at fractional
      position x with min(x, 1-x) >= rho, and carries weight p >= alpha;
  H3  the common mean pendant satisfies y >= pendant_floor > 0;
  H4  the supplied vector `upper` dominates the true model vector coordinatewise,
      upper_r >= d_r for every observed reference r.
Under H1-H4, if a source lies in original edge e then for every endpoint w of e
and every pair of observed leaves r,s whose path meets w and avoids the interior
of e,  (upper_r + upper_s - d(r,s))/2 - pendant_floor >= alpha * offset_w(e),
where offset_w(e) is the distance from w to the source point. An edge whose
measured excess falls below that bound therefore contains no source: this is an
EXCLUSION certificate. Retention is not evidence of presence, and alpha, rho and
pendant_floor are SUPPLIED, not estimated from the data.

Reduced panel. If references are removed, the hull edge J between surviving
branch vertices u,v consists of original edges e_1..e_m with lengths L_1..L_m
(ordered from u). A source inside e_k is at offset_u >= prefix_k + rho*L_k and
offset_v >= suffix_k + rho*L_k. Those per-original-edge bounds are used, which
reduces to the frozen threshold alpha*rho*L exactly when m = 1.
"""
import os
for _k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_k] = '1'
import hashlib
import json
import sys
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / 'scripts') not in sys.path:
    sys.path.insert(0, str(ROOT / 'scripts'))
BANK = ROOT / 'data/stage287_instances/bank'

from stage44_certificate import path_screen            # frozen screen, imported unchanged
from stage288_common import (serial, sha, bank_index, adjacency_lists, leaves_of,
                             bfs_from, predict, clade_removals, nearest_source_removal,
                             random_removal, run_pool, dependency_closure)

DUMMY_ALPHA, DUMMY_RHO, DUMMY_FLOOR = F(1, 10 ** 9), F(1, 4), F(1, 10 ** 9)


# --------------------------------------------------------------- bank loading

def load_any(entry):
    """Read one frozen bank case of ANY family, verifying its registered SHA256."""
    path = ROOT / entry['path']
    digest = sha(path)
    assert digest == entry['sha256'], f"bank case hash changed: {entry['name']}"
    raw = json.loads(path.read_text())
    sources = None
    if raw.get('generating_sources'):
        details = raw.get('source_details') or []
        sources = []
        for i, s in enumerate(raw['generating_sources']):
            p = s.get('p')
            if p is None and i < len(details):
                p = details[i].get('p')            # biotrees: weights live in source_details
            sources.append(dict(edge=int(s['edge']), p=(None if p is None else F(p)), x=F(s['x'])))
    pend = raw.get('generating_pendant', raw.get('pendant'))
    return dict(name=entry['name'], family=entry['family'], path=entry['path'], sha256=digest,
                edges=[(u, v, F(length)) for u, v, length in raw['edges']],
                observed={int(r): F(d) for r, d in raw['observed'].items()},
                cap=raw.get('cap'), pendant=(None if pend is None else F(pend)),
                sources=sources, provenance=raw.get('provenance', {}),
                reference_set_size=entry['reference_set_size'],
                K_generating=entry['K_generating'], tree_leaves=entry['tree_leaves'],
                distance_source=entry['distance_source'])


def model_vector(case, targets=None):
    """Exact PDD report of the published generating mixture, d_r = y + sum p_i dist(r,q_i)."""
    if case['sources'] is None or case['pendant'] is None:
        return None
    if any(s['p'] is None for s in case['sources']):
        return None
    return predict(case['edges'], case['sources'], case['pendant'],
                   sorted(targets if targets is not None else case['observed']))


# ------------------------------------------------- endpoint-local excesses

def two_smallest(values):
    first = second = None
    for v in values:
        if first is None or v < first:
            first, second = v, first
        elif second is None or v < second:
            second = v
    return first, second


def endpoint_raw(edges, upper):
    """Half-sums (upper_r + upper_s - d(r,s))/2 at both endpoints of every edge.

    Returns (raw, frozen_records). raw[j] = (raw_u, raw_v), each None when the
    endpoint has no eligible avoiding leaf pair. Excess at pendant floor f and
    error allowance eta is raw + eta - f, because adding eta to every coordinate
    of `upper` adds eta to every message and the floor enters additively.
    """
    retained, records, message = path_screen(edges, upper, DUMMY_FLOOR, DUMMY_ALPHA, DUMMY_RHO)
    adj = adjacency_lists(edges)
    raw = []
    for j, (u, v, _ell) in enumerate(edges):
        out = []
        for endpoint, away in ((u, v), (v, u)):
            a, b = two_smallest(message[endpoint, w][0] for w, _l in adj[endpoint] if w != away)
            out.append(None if b is None else (a + b) / 2)
        # Cross-check against the frozen witness, which is the minimum of the two.
        witness = records[j]['witness']
        mine = [x for x in out if x is not None]
        if witness is None:
            assert not mine, 'frozen screen reports no witness where one exists'
        else:
            assert min(mine) - DUMMY_FLOOR == witness['upper_excess'], 'endpoint split disagrees with frozen witness'
        raw.append(tuple(out))
    return raw, records, retained


# --------------------------------------------------------------- hull projection

def hull_projection(edges, keep):
    """Induced subtree on `keep` with degree-two suppression.

    Returns (hull_edges, members) where hull_edges[j] = (u, v, L) and
    members[j] = ordered list of original edge indices from u to v.
    """
    keep = set(keep)
    assert len(keep) >= 2
    adj = {}
    for j, (u, v, ell) in enumerate(edges):
        adj.setdefault(u, []).append((v, ell, j))
        adj.setdefault(v, []).append((u, ell, j))
    root = next(iter(keep))
    parent, order = {root: None}, [root]
    stack = [root]
    while stack:
        x = stack.pop()
        for y, _l, _j in adj[x]:
            if y not in parent:
                parent[y] = x
                order.append(y)
                stack.append(y)
    below = {}
    for x in reversed(order):
        below[x] = (1 if x in keep else 0) + sum(below[y] for y, _l, _j in adj[x] if parent.get(y) == x)
    total = len(keep)
    live = set()
    for j, (u, v, _ell) in enumerate(edges):
        child = u if parent.get(u) == v else v
        if 0 < below[child] < total:
            live.add(j)
    degree = {}
    for j in live:
        u, v, _ = edges[j]
        degree[u] = degree.get(u, 0) + 1
        degree[v] = degree.get(v, 0) + 1
    branch = {w for w, d in degree.items() if d != 2} | (keep & set(degree))
    hull_edges, members = [], []
    seen = set()
    for start in sorted(branch):
        for nxt, ell, j0 in adj[start]:
            if j0 not in live or j0 in seen:
                continue
            path = [j0]
            seen.add(j0)
            length = ell
            prev, cur = start, nxt
            while cur not in branch:
                step = next((y, l, j) for y, l, j in adj[cur] if j in live and j != path[-1])
                y, l, j = step
                path.append(j)
                seen.add(j)
                length += l
                prev, cur = cur, y
            hull_edges.append((start, cur, length))
            members.append(path)
    assert len(seen) == len(live)
    return hull_edges, members


def member_offsets(edges, hull_edges, members):
    """For each hull edge, per-original-edge (prefix, suffix, length) from endpoint u."""
    out = []
    for (u, _v, _L), path in zip(hull_edges, members):
        lengths = [edges[j][2] for j in path]
        node = u
        ordered = []
        for j, ell in zip(path, lengths):
            a, b, _ = edges[j]
            assert node in (a, b)
            ordered.append(ell)
            node = b if node == a else a
        prefixes, acc = [], F(0)
        for ell in ordered:
            prefixes.append(acc)
            acc += ell
        suffixes = [acc - prefixes[i] - ordered[i] for i in range(len(ordered))]
        out.append(list(zip(path, prefixes, suffixes, ordered)))
    return out


# --------------------------------------------------------------- the screen

def excluded_original(raw, offsets, eta, floor, alpha, rho):
    """Original edge indices certifiably excluded, given supplied H2/H3 bounds."""
    out = []
    for (raw_u, raw_v), items in zip(raw, offsets):
        eu = None if raw_u is None else raw_u + eta - floor
        ev = None if raw_v is None else raw_v + eta - floor
        for j, prefix, suffix, ell in items:
            bound_u = alpha * (prefix + rho * ell)
            bound_v = alpha * (suffix + rho * ell)
            if (eu is not None and eu < bound_u) or (ev is not None and ev < bound_v):
                out.append(j)
    return out


def identity_offsets(edges):
    return [[(j, F(0), F(0), ell)] for j, (_u, _v, ell) in enumerate(edges)]


# --------------------------------------------------------- panel-level screen

def screen_panel(case, keep):
    """Hull-project onto `keep`, run the frozen screen there, return raw excesses.

    Returns a dict with the hull, the per-original-edge offsets, the endpoint
    half-sums, and the set of original edges that the reduced panel can reach at
    all (edges off the hull are unconstrained by construction).
    """
    keep = sorted(keep)
    edges = case['edges']
    hull_edges, members = hull_projection(edges, keep)
    upper = {r: case['observed'][r] for r in keep}
    raw, records, frozen_keep = endpoint_raw(hull_edges, upper)
    offsets = member_offsets(edges, hull_edges, members)
    covered = sorted(j for path in members for j in path)
    return dict(hull_edges=hull_edges, members=members, offsets=offsets, raw=raw,
                covered=covered, n_hull=len(hull_edges), n_covered=len(covered))
