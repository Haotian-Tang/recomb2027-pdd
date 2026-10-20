#!/usr/bin/env python3
"""Stage 298 shared library: Theta_k on REAL distance data.

PURPOSE.  Every incompleteness measurement in this project so far runs on
instances with synthetic planted sources.  Stage 296 built the instrument that
makes a real-data statement possible without ground truth: Theta_1, the exact
smallest tolerance at which a SINGLE source explains the observations, is a
rational number computable from the data alone.  This module computes it on
three real collections and places all of them on one axis.

Nothing pre-existing is modified.  Stage 296's exact rational LP
(`stage296_common.theta_support`), its support sweep and its independent
witness verifier are imported unmodified, so every number here is on the same
scale as stage 296's.

THE k = 1 CLOSED FORM (new here; needed only because the travellers'-diarrhoea
backbone has |E| = 21,147 and |R| = 10,575, where materialising the full
leaf-distance matrix the frozen sweep wants is 2.2e8 rationals).

For a single source at fraction x of edge e = (w, c) with length L, where w is
the parent endpoint and IN = observed leaves below c, OUT = the rest,

    dist(r, q) = d(r, w) + x L           (r in OUT)
               = d(r, c) + (1 - x) L
               = d(r, w) - x L           (r in IN),

so with h_r = d(r, w) - d_r and sigma_r = +1 on OUT, -1 on IN the deviation at
r is h_r + sigma_r x L + y.  Writing

    A_out = max_{OUT} h_r,  B_out = min_{OUT} h_r,
    A_in  = max_{IN}  h_r,  B_in  = min_{IN}  h_r,
    M(x)  = max(A_out + xL, A_in - xL),  m(x) = min(B_out + xL, B_in - xL),

the inner minimisation over y >= 0 of max_r |h_r + sigma_r x L + y| equals
(M - m)/2 when M + m <= 0 (the unconstrained optimum y* = -(M+m)/2 is feasible)
and M otherwise, i.e. exactly max(M(x), (M(x) - m(x))/2).  That is a maximum of
the six affine functions

    A_out + xL,  A_in - xL,  (A_out-B_out)/2,  (A_out-B_in)/2 + xL,
    (A_in-B_out)/2 - xL,  (A_in-B_in)/2,

which is minimised over x in [0, 1] exactly at an endpoint or at a crossing of
two of the six lines - finitely many exact rational candidates.  The four
aggregates for every edge at once come from one rerooting DP in O(|V|).

The closed form is NOT trusted: `theta1_all_edges` is validated against the
frozen `stage296_common.theta_support` by exact rational equality on complete
edge sweeps of small real cases and on random edges of the large ones, and the
attaining witness is always re-checked with the frozen
min(d(r,u)+xL, d(r,v)+(1-x)L) verifier.
"""
import os
for _k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_k] = '1'
import json
import sys
from collections import deque
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / 'scripts') not in sys.path:
    sys.path.insert(0, str(ROOT / 'scripts'))
import stage288_common as common            # noqa: E402  frozen bank access
import stage296_common as s296              # noqa: E402  frozen exact LP

OUT = ROOT / 'results' / 'stage298'
BANK = ROOT / 'data' / 'stage287_instances' / 'bank'
ZERO, ONE = F(0), F(1)


# ------------------------------------------------------------------- loading

def bank_cases(family):
    """(name, path) for every case of a bank family, sorted."""
    index = common.bank_index()
    return sorted((e['name'], e) for e in index if e['family'] == family)


def read_case(entry):
    """Load one bank case with its registered sha256 checked."""
    path = ROOT / entry['path']
    assert common.sha(path) == entry['sha256'], 'bank case hash mismatch: %s' % path
    return json.loads(path.read_text())


def case_arrays(case):
    """(edges as [(u, v, F)], sorted observed leaves, aligned observations)."""
    edges = [(int(u), int(v), F(L)) for u, v, L in case['edges']]
    leaves = sorted(int(r) for r in case['observed'])
    dvals = [F(case['observed'][str(r)]) if str(r) in case['observed']
             else F(case['observed'][r]) for r in leaves]
    return edges, leaves, dvals


def median_observed(dvals):
    """Stage 296's convention: sorted(values)[len // 2]."""
    s = sorted(dvals)
    return s[len(s) // 2]


# ------------------------------------------------- exact min of max of lines

def min_max_affine(lines, lo=ZERO, hi=ONE):
    """Exact min over [lo, hi] of max_i (slope_i * x + const_i).

    Convex piecewise linear: the minimiser is an endpoint or a crossing of two
    of the lines.  Returns (value, argmin) as exact rationals.
    """
    cands = [lo, hi]
    n = len(lines)
    for i in range(n):
        si, bi = lines[i]
        for j in range(i + 1, n):
            sj, bj = lines[j]
            if si == sj:
                continue
            x = (bj - bi) / (si - sj)
            if lo < x < hi:
                cands.append(x)
    best, arg = None, None
    for x in cands:
        v = max(s * x + b for s, b in lines)
        if best is None or v < best:
            best, arg = v, x
    return best, arg


# ------------------------------------------------------------- rerooting DP

def _root_tree(edges, nv):
    adj = [[] for _ in range(nv)]
    for i, (u, v, L) in enumerate(edges):
        adj[u].append((v, L, i))
        adj[v].append((u, L, i))
    present = [v for v in range(nv) if adj[v]]
    assert present, 'no vertices'
    root = present[0]
    for v in present:
        if len(adj[v]) > 1:
            root = v
            break
    parent = [-1] * nv
    pedge = [-1] * nv
    plen = [None] * nv
    order = [root]
    seen = [False] * nv
    seen[root] = True
    dq = deque([root])
    while dq:
        x = dq.popleft()
        for y, L, i in adj[x]:
            if not seen[y]:
                seen[y] = True
                parent[y], pedge[y], plen[y] = x, i, L
                order.append(y)
                dq.append(y)
    assert all(seen[v] for v in present), 'tree is disconnected'
    children = [[] for _ in range(nv)]
    for v in order[1:]:
        children[parent[v]].append(v)
    return adj, root, parent, pedge, plen, order, children


def _cmb_max(a, b):
    if a is None:
        return b
    if b is None:
        return a
    return a if a > b else b


def _cmb_min(a, b):
    if a is None:
        return b
    if b is None:
        return a
    return a if a < b else b


def edge_aggregates(edges, leaves, dvals, nv=None):
    """For every edge: (A_in, B_in, A_out, B_out, child, L) in the parent frame.

    A_* / B_* are max / min of h_r = d(r, w) - d_r over the two sides of the
    edge, w being the endpoint nearer the root.  `None` marks an empty side.
    One rerooting DP, O(|V|) exact rational operations.
    """
    if nv is None:
        nv = 1 + max(max(u, v) for u, v, _ in edges)
    obs = {r: d for r, d in zip(leaves, dvals)}
    adj, root, parent, pedge, plen, order, children = _root_tree(edges, nv)
    Aup = [None] * nv
    Bup = [None] * nv
    for v in reversed(order):
        a = b = None
        if v in obs:
            a = b = -obs[v]
        for c in children[v]:
            if Aup[c] is not None:
                a = _cmb_max(a, Aup[c] + plen[c])
                b = _cmb_min(b, Bup[c] + plen[c])
        Aup[v], Bup[v] = a, b
    Adn = [None] * nv
    Bdn = [None] * nv
    for w in order:
        ch = children[w]
        if not ch:
            continue
        vals = [(Aup[c] + plen[c] if Aup[c] is not None else None,
                 Bup[c] + plen[c] if Bup[c] is not None else None) for c in ch]
        m = len(ch)
        pre_a = [None] * (m + 1)
        pre_b = [None] * (m + 1)
        suf_a = [None] * (m + 1)
        suf_b = [None] * (m + 1)
        for i in range(m):
            pre_a[i + 1] = _cmb_max(pre_a[i], vals[i][0])
            pre_b[i + 1] = _cmb_min(pre_b[i], vals[i][1])
        for i in range(m - 1, -1, -1):
            suf_a[i] = _cmb_max(suf_a[i + 1], vals[i][0])
            suf_b[i] = _cmb_min(suf_b[i + 1], vals[i][1])
        base_a = _cmb_max(Adn[w], -obs[w] if w in obs else None)
        base_b = _cmb_min(Bdn[w], -obs[w] if w in obs else None)
        for i, c in enumerate(ch):
            a = _cmb_max(base_a, _cmb_max(pre_a[i], suf_a[i + 1]))
            b = _cmb_min(base_b, _cmb_min(pre_b[i], suf_b[i + 1]))
            Adn[c] = None if a is None else a + plen[c]
            Bdn[c] = None if b is None else b + plen[c]
    out = []
    for i, (u, v, L) in enumerate(edges):
        c = v if parent[v] == u and pedge[v] == i else u
        assert parent[c] != -1 and pedge[c] == i, 'edge %d not in the rooting' % i
        A_in = None if Aup[c] is None else Aup[c] + L
        B_in = None if Bup[c] is None else Bup[c] + L
        A_out = None if Adn[c] is None else Adn[c] - L
        B_out = None if Bdn[c] is None else Bdn[c] - L
        out.append((A_in, B_in, A_out, B_out, c, L))
    return out


def theta1_edge(A_in, B_in, A_out, B_out, L):
    """theta({e}) and the attaining x, measured from the parent endpoint."""
    lines = []
    if A_out is not None:
        lines.append((L, A_out))
    if A_in is not None:
        lines.append((-L, A_in))
    if A_out is not None and B_out is not None:
        lines.append((ZERO, (A_out - B_out) / 2))
    if A_out is not None and B_in is not None:
        lines.append((L, (A_out - B_in) / 2))
    if A_in is not None and B_out is not None:
        lines.append((-L, (A_in - B_out) / 2))
    if A_in is not None and B_in is not None:
        lines.append((ZERO, (A_in - B_in) / 2))
    if not lines:
        return None, None
    return min_max_affine(lines)


def theta1_all_edges(edges, leaves, dvals, nv=None):
    """Exact theta({e}) for every edge, plus Theta_1 and its argmin.

    Returns dict(theta=[...], x=[...], best=index, theta1=value).
    """
    aggs = edge_aggregates(edges, leaves, dvals, nv=nv)
    thetas, xs, best = [], [], None
    for i, (A_in, B_in, A_out, B_out, c, L) in enumerate(aggs):
        t, x = theta1_edge(A_in, B_in, A_out, B_out, L)
        thetas.append(t)
        # express x in the orientation of edges[i] = (u, v): x is the fraction
        # from u toward v, while the DP measures it from the parent endpoint.
        if x is not None and edges[i][0] != c:
            xs.append(x)
        elif x is not None:
            xs.append(ONE - x)
        else:
            xs.append(None)
        if t is not None and (best is None or t < thetas[best]):
            best = i
    return dict(theta=thetas, x=xs, best=best,
                theta1=(None if best is None else thetas[best]))


# ----------------------------------------------------- frozen cross-checking

def frozen_theta_support(edges, leaves, dvals, support, cutoff=None):
    """stage296_common.theta_support on one support, building its columns with
    the frozen BFS only from that support's endpoints (so it never needs the
    full |V| x |R| leaf-distance matrix)."""
    ends = sorted({x for e in support for x in edges[e][:2]})
    dist = common.bfs_from([(u, v, L) for u, v, L in edges], ends)
    cols = []
    for e in support:
        u, v, _L = edges[e]
        cols.append(tuple(dist[u][r] for r in leaves))
        cols.append(tuple(dist[v][r] for r in leaves))
    return s296.theta_support(cols, dvals, len(support), cutoff=cutoff)


def verify_witness_exact(case_edges, support, z, k, observed_map, leaves):
    """Frozen independent verifier: recompute with min(d(r,u)+xL, d(r,v)+(1-x)L)."""
    return s296.independent_deviation(case_edges, support, z, k, observed_map,
                                      leaves=leaves)


def serial(value):
    return common.serial(value)

