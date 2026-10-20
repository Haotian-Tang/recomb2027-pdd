#!/usr/bin/env python3
"""Stage 300 shared library: how does the published tool's OWN returned mixture fit?

THE COMPARISON.  Stage 298 measured, on real estimated distances, the L-infinity
residual of the TRUE generating configuration and the exact tolerance Theta_1 at
which a SINGLE source still explains the observations.  What was missing is the
answer a practitioner actually acts on: the mixture DecoDiPhy returns.  This
module evaluates an arbitrary configuration - truth, DecoDiPhy's answer, or the
best single source - at every observed reference leaf, in two norms:

    R      = max_r |pred_r - d_r|                     (L-infinity; Theta_1's norm)
    SSR    = sum_r (pred_r - d_r)^2                    (DecoDiPhy's own objective)

with pred_r = y + sum_i p_i * dist_T(r, q_i) and the frozen
dist(r, q) = min(d(r,u) + x L, d(r,v) + (1 - x) L).

R is computed in exact rational arithmetic (the configuration's floats are
converted by their exact binary value), so `R >= Theta_1` is an exact decision.
SSR is computed in float64 because the quantity it is compared against -
DecoDiPhy's reported loss - is itself float64; DecoDiPhy's objective is
`np.sum((d - (DA @ p + CAl @ x + y)) ** 2)` (decodiphy/optimize.py), i.e. the
sum of squared residuals of the same forward model.

FAIRNESS.  DecoDiPhy minimises SSR, not R.  The L2 analogue of Theta_1 is
therefore also provided: `best_single_source_l2` is the exact constrained
minimum of SSR over every original edge, every position x in [0, 1] and every
pendant y >= 0 with p = 1 - the same feasible set DecoDiPhy's own QP uses at
k = 1 (p >= 0, 0 <= t <= p, y >= 0, sum p = 1), searched exhaustively rather
than by hill climbing.

Nothing pre-existing is modified.  The encoding, the vertex-distance routine and
the bank accessors are imported unmodified from stage 288/289/298.
"""
import os
for _k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_k, '1')
import hashlib
import json
import sys
from collections import deque
from fractions import Fraction as F
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / 'scripts') not in sys.path:
    sys.path.insert(0, str(ROOT / 'scripts'))
import stage288_common as common                 # noqa: E402  frozen bank access
import stage289_decodiphy_bridge as bridge       # noqa: E402  frozen encoding
import stage298_common as s298                   # noqa: E402  frozen Theta_1 machinery

OUT = ROOT / 'results' / 'stage300'
S298 = ROOT / 'results' / 'stage298'
COLLECTIONS = ('krepp_ani', 'krepp_placement', 'travelers')
FAMILY_OF = {'krepp_ani': 'krepp_estimated_10k',
             'krepp_placement': 'krepp_placement_10k',
             'travelers': 'published_travelers_diarrhea'}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def serial(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {str(k): serial(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [serial(v) for v in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    return value


# ----------------------------------------------------------------- loading

def case_entries(collection):
    """Registered entries (name, path, sha256, ...) for one real collection."""
    if collection == 'krepp_placement':
        manifest = json.loads((S298 / 'stage298_placement_manifest.json').read_text())
        rows = []
        for c in manifest['cases']:
            rows.append(dict(name=c['name'], path=c['path'], sha256=c['sha256'],
                             family='krepp_placement_10k',
                             K_generating=c.get('K_generating'),
                             n_edges=c.get('n_edges')))
        return sorted(rows, key=lambda r: r['name'])
    family = FAMILY_OF[collection]
    rows = [dict(e) for e in common.bank_index() if e['family'] == family]
    return sorted(rows, key=lambda r: r['name'])


def load_case(entry):
    """Raw case JSON with its registered sha256 checked."""
    path = ROOT / entry['path']
    got = sha(path)
    assert got == entry['sha256'], 'case hash mismatch: %s' % path
    return json.loads(path.read_text())


def case_arrays(case):
    return s298.case_arrays(case)


def median_observed(dvals):
    return s298.median_observed(dvals)


def truth_configuration(case):
    """(sources, pendant) of the generating configuration, exact, or None."""
    srcs = case.get('generating_sources')
    if not srcs:
        return None
    pend = case.get('generating_pendant')
    if pend is None:
        pend = case.get('pendant')
    if pend is None:
        return None
    sources = [dict(edge=int(s['edge']), p=F(s['p']), x=F(s['x'])) for s in srcs]
    return sources, F(pend)


# ------------------------------------------------- exact forward evaluation

def _distance_maps(edges, vertices):
    return {v: bridge.vertex_distances(edges, v) for v in sorted(set(vertices))}


def predictions_exact(edges, leaves, sources, pendant, maps=None):
    """y + sum_i p_i * min(d(r,u) + x L, d(r,v) + (1 - x) L) at every leaf, exact."""
    need = [v for s in sources for v in edges[int(s['edge'])][:2]]
    maps = maps if maps is not None else _distance_maps(edges, need)
    pend = pendant if isinstance(pendant, F) else F(pendant)
    out = []
    for r in leaves:
        total = pend
        for s in sources:
            u, v, length = edges[int(s['edge'])]
            x = s['x'] if isinstance(s['x'], F) else F(s['x'])
            p = s['p'] if isinstance(s['p'], F) else F(s['p'])
            total += p * min(maps[u][r] + x * length, maps[v][r] + (1 - x) * length)
        out.append(total)
    return out


def predictions_interpolation(edges, leaves, sources, pendant, maps=None):
    """Independent route: dist(r, q) = (1 - x) d(r,u) + x d(r,v) (stage 296's form)."""
    need = [v for s in sources for v in edges[int(s['edge'])][:2]]
    maps = maps if maps is not None else _distance_maps(edges, need)
    pend = pendant if isinstance(pendant, F) else F(pendant)
    out = []
    for r in leaves:
        total = pend
        for s in sources:
            u, v, _ = edges[int(s['edge'])]
            x = s['x'] if isinstance(s['x'], F) else F(s['x'])
            p = s['p'] if isinstance(s['p'], F) else F(s['p'])
            total += p * ((1 - x) * maps[u][r] + x * maps[v][r])
        out.append(total)
    return out


def residuals_exact(edges, leaves, dvals, sources, pendant, route='min'):
    pred = (predictions_exact if route == 'min' else predictions_interpolation)(
        edges, leaves, sources, pendant)
    return [a - b for a, b in zip(pred, dvals)]


def linf(residuals):
    return max((abs(z) for z in residuals), default=F(0))


def ssr_float(residuals):
    return float(sum(float(z) ** 2 for z in residuals))


def fit_metrics(edges, leaves, dvals, sources, pendant):
    """Exact L-infinity residual and float64 SSR of one configuration."""
    res = residuals_exact(edges, leaves, dvals, sources, pendant)
    absr = [abs(z) for z in res]
    return dict(linf=linf(res), ssr=ssr_float(res),
                median_abs=sorted(absr)[len(absr) // 2] if absr else F(0),
                n_leaves=len(leaves))


# ------------------------------------- best single source in DecoDiPhy's norm

def _root_and_preorder(edges, leaves):
    """Root at a vertex of degree >= 3; DFS preorder with contiguous leaf slots."""
    adjacency = {}
    for index, (u, v, length) in enumerate(edges):
        adjacency.setdefault(u, []).append((v, float(length), index))
        adjacency.setdefault(v, []).append((u, float(length), index))
    root = None
    for v in sorted(adjacency):
        if len(adjacency[v]) >= 3:
            root = v
            break
    if root is None:
        root = min(adjacency)
    leafset = set(leaves)
    slot = {}
    lo = {}
    hi = {}
    parent = {root: None}
    order = []
    stack = [(root, iter(sorted(adjacency[root])), None, None)]
    lo[root] = 0
    counter = 0
    if root in leafset:
        slot[root] = counter
        counter += 1
    frames = [[root, list(sorted(adjacency[root])), 0, None, None]]
    edge_to_child = {}
    while frames:
        v, nbrs, pos, plen, pedge = frames[-1]
        if pos < len(nbrs):
            frames[-1][2] += 1
            u, length, eindex = nbrs[pos]
            if u == parent[v]:
                continue
            parent[u] = v
            edge_to_child[eindex] = u
            lo[u] = counter
            if u in leafset:
                slot[u] = counter
                counter += 1
            frames.append([u, [n for n in sorted(adjacency[u])], 0, length, eindex])
            order.append((u, v, length, eindex))
        else:
            hi[v] = counter
            frames.pop()
    return adjacency, root, parent, order, lo, hi, slot, edge_to_child


def best_single_source_l2(edges, leaves, dvals, return_profile=False):
    """Exact constrained minimum of SSR over every edge, x in [0,1], y >= 0, p = 1.

    Convex quadratic in (a = xL, y) per edge; the minimum is attained at the
    interior stationary point when feasible, otherwise on one of the three
    boundary lines a = 0, a = L, y = 0 (each a clamped one-dimensional
    quadratic).  All candidates are evaluated and the best kept.
    """
    adjacency, root, parent, order, lo, hi, slot, edge_to_child = _root_and_preorder(edges, leaves)
    n = len(leaves)
    d = np.empty(n, dtype=np.float64)
    for r in leaves:
        d[slot[r]] = float(dvals[leaves.index(r)]) if False else 0.0
    # aligned fill (leaves and dvals are aligned lists)
    for position, r in enumerate(leaves):
        d[slot[r]] = float(dvals[position])

    # distance from the root to every vertex, then the root column over leaves
    dist_root = {root: 0.0}
    queue = deque([root])
    while queue:
        v = queue.popleft()
        for u, length, _ in adjacency[v]:
            if u not in dist_root:
                dist_root[u] = dist_root[v] + length
                queue.append(u)
    col_root = np.empty(n, dtype=np.float64)
    for r in leaves:
        col_root[slot[r]] = dist_root[r]

    best = None
    profile = {}
    columns = {root: col_root}
    # iterative descent reusing the parent column
    children = {}
    for (u, v, length, eindex) in order:
        children.setdefault(v, []).append((u, length, eindex))
    stack = [root]
    live = {root: col_root}
    pending = {v: len(children.get(v, [])) for v in adjacency}
    while stack:
        w = stack.pop()
        col_w = live[w]
        kids = children.get(w, [])
        if not kids:
            live.pop(w, None)
            continue
        sum_h = float(np.sum(col_w) - np.sum(d))
        h = col_w - d
        sum_h2 = float(np.dot(h, h))
        for (c, length, eindex) in kids:
            a_lo, a_hi = lo[c], hi[c]
            n_in = a_hi - a_lo
            sum_in_h = float(np.sum(h[a_lo:a_hi])) if n_in else 0.0
            sig_h = sum_h - 2.0 * sum_in_h
            sig = float(n - 2 * n_in)
            L = length

            def value(a, y):
                return (sum_h2 + n * a * a + n * y * y + 2.0 * a * sig_h
                        + 2.0 * y * sum_h + 2.0 * a * y * sig)

            cands = []
            det = n * n - sig * sig
            if abs(det) > 1e-12 * max(1.0, n * n):
                a_star = (-n * sig_h + sig * sum_h) / det
                y_star = (-n * sum_h + sig * sig_h) / det
                if -1e-15 <= a_star <= L + 1e-15 and y_star >= -1e-15:
                    cands.append((min(max(a_star, 0.0), L), max(y_star, 0.0)))
            cands.append((min(max(-sig_h / n, 0.0), L), 0.0))          # y = 0
            cands.append((0.0, max(-sum_h / n, 0.0)))                   # a = 0
            cands.append((L, max(-(sum_h + L * sig) / n, 0.0)))         # a = L
            cands.extend([(0.0, 0.0), (L, 0.0)])
            va, vy = min(cands, key=lambda z: value(*z))
            vbest = value(va, vy)
            if return_profile:
                profile[eindex] = vbest
            if best is None or vbest < best['ssr']:
                # a is measured from the PARENT endpoint w; the bank measures x from
                # edges[e][0], which is the child on about half of the edges.
                x_from_w = (va / L if L > 0 else 0.0)
                u0 = edges[int(eindex)][0]
                best = dict(ssr=vbest, edge=int(eindex),
                            x=(x_from_w if u0 == w else 1.0 - x_from_w),
                            x_from_parent=x_from_w, parent=int(w), child=int(c),
                            a=va, y=vy, edge_length=L, n_in=int(n_in))
            col_c = col_w + length
            col_c[a_lo:a_hi] -= 2.0 * length
            live[c] = col_c
            stack.append(c)
        live.pop(w, None)
    if return_profile:
        return best, profile
    return best


# ------------------------------------------------------- stage 298 imports

def theta1_table():
    """Exact Theta_1 per case, imported from stage 298 (never recomputed here)."""
    table = {}
    analysis = json.loads((S298 / 'stage298_analysis.json').read_text())
    for row in analysis['small_cases']:
        ladder = {int(k): F(v) for k, v in row['theta'].items()}
        table[row['name']] = dict(theta1=ladder[1], theta_ladder=ladder,
                                  source='stage298_analysis.small_cases',
                                  own_median_observed=F(row['own_median_observed']),
                                  family=row['family'])
    td = json.loads((S298 / 'stage298_td_summary.json').read_text())
    for row in td['rows']:
        table[row['name']] = dict(theta1=F(row['theta1']), source='stage298_td_summary.rows',
                                  own_median_observed=F(row['median_observed']),
                                  family='published_travelers_diarrhea')
    return table


def truth_residual_table():
    """Exact L-infinity residual of the TRUE configuration, imported from stage 298."""
    table = {}
    res = json.loads((S298 / 'stage298_residuals.json').read_text())
    for row in res['rows']:
        if row.get('status') != 'ok':
            continue
        table[row['case']] = dict(residual_max=F(row['residual_max']),
                                  residual_median=F(row['residual_median']),
                                  n_shared=row['n_shared'], family=row['family'])
    return table
