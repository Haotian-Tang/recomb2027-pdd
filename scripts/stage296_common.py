#!/usr/bin/env python3
"""Stage 296 shared library: the eps-tolerant source-budget problem.

MODEL.  Known reference tree T (positive rational edge lengths), observed
reference set R with exact rational observations d_r, k sources strictly
interior to k distinct ORIGINAL edges, weights p_i > 0 summing to one, common
mean pendant y > 0.  For a point q at fraction x of edge (u, v) and ANY vertex
r, dist(r, q) = (1 - x) d(r, u) + x d(r, v) exactly: the path from r to q
leaves the edge through whichever endpoint separates r from q, and both cases
give that interpolation.  With endpoint masses

    alpha_i = p_i (1 - x_i) > 0,     beta_i = p_i x_i > 0,

the model reads  d_r = y + sum_i (alpha_i D[u_i][r] + beta_i D[v_i][r])  and is
AFFINE in (alpha, beta, y) >= 0 under sum_i (alpha_i + beta_i) = 1.

EPS-TOLERANT PROBLEM.  Given a tolerance eps >= 0, the minimum number k of
strict interior sources on distinct original edges with
|pred_r - d_r| <= eps for every r in R.

THRESHOLD FORM (what this module computes).  For a FIXED support S put

    theta(S) = min over the CLOSED set {alpha, beta, y >= 0, sum = 1}
               of max_r |pred_r - d_r|                               (one LP),
    Theta_k  = min over all supports S of size k of theta(S).

  (i)  For every eps > theta(S) the STRICT system on S is feasible.  If z_c
       attains theta(S) and z_s is any strictly feasible box point with finite
       deviation M, then for L in (0, 1] the point (1-L) z_c + L z_s is
       strictly feasible and, the prediction being affine, has deviation at
       most (1-L) theta(S) + L M, which tends to theta(S) as L -> 0.
  (ii) Theta_k is non-increasing in k: a size-(k-1) solution is reproduced by
       any size-k support containing it with the extra endpoint masses zero,
       which the CLOSED system allows.
  (iii) Hence for eps > 0 the eps-tolerant minimum is m(eps) = min{k : Theta_k
       < eps}, and m(0) = min{k : some size-k support is STRICTLY feasible at
       exact equality}, decided here by a separate max-slack LP.  At the
       measure-zero coincidence eps == Theta_k > 0 the strict system may or may
       not be feasible; the registered convention is the strict '<' above (the
       conservative, larger-m choice), and strict_feasible_at() can settle such
       a tie exactly if one ever arises.

No part of the frozen structural localisation (leaf envelope h, weighted
deficit bound, positive-anchor filters, and therefore solve_partial,
stage292_all_minimisers, stage295_reduced_minimum) is used for eps > 0: all of
it is derived from exact equality and is unsound under a tolerance.  Supports
are enumerated exhaustively and every decision is an exact rational LP written
here.  The frozen machinery is called only in the eps = 0 validation, where it
is sound.

Nothing pre-existing is modified.
"""
import os
for _k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_k] = '1'
import itertools
import sys
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / 'scripts') not in sys.path:
    sys.path.insert(0, str(ROOT / 'scripts'))
import stage288_common as common            # noqa: E402  bank access, regimes, pool

ZERO, ONE = F(0), F(1)
BLAND_AFTER = 60


# --------------------------------------------------------------- exact LP core

def lp_minimize(n, cost, rows, rhs, rels, iteration_cap=20000):
    """Exact rational LP:  min cost.z  s.t.  row.z {'eq'|'ge'} rhs,  z >= 0.

    Dense two-phase tableau simplex over Fractions.  Dantzig pricing, switching
    to Bland's rule after BLAND_AFTER iterations so termination is guaranteed
    under degeneracy.  Returns dict(status, value, z, iterations).
    """
    m = len(rows)
    ncol = n
    surplus = {}
    for j, rel in enumerate(rels):
        if rel == 'ge':
            surplus[j] = ncol
            ncol += 1
    nstruct = ncol
    T, basis = [], []
    for j, row in enumerate(rows):
        r = list(row) + [ZERO] * (nstruct - n)
        if rels[j] == 'ge':
            r[surplus[j]] = -ONE
        b = rhs[j]
        if b < 0:
            r = [-x for x in r]
            b = -b
        T.append(r + [ONE if i == j else ZERO for i in range(m)] + [b])
        basis.append(nstruct + j)
    total = nstruct + m
    obj = [ZERO] * (total + 1)
    for j in range(m):
        Tj = T[j]
        for i in range(nstruct):
            obj[i] -= Tj[i]
        obj[total] -= Tj[total]
    T.append(obj)
    status, iters = _simplex(T, basis, nstruct, total, iteration_cap)
    if status != 'optimal':
        return dict(status='phase1_' + status, value=None, z=None, iterations=iters)
    if -T[-1][total] > 0:
        return dict(status='infeasible', value=None, z=None, iterations=iters)
    keep = list(range(m))
    for i in range(m):
        if basis[i] >= nstruct:
            col = next((j for j in range(nstruct) if T[i][j] != 0), None)
            if col is None:
                keep.remove(i)
            else:
                _pivot(T, basis, i, col, total)
    if len(keep) < m:
        T = [T[i] for i in keep] + [T[-1]]
        basis = [basis[i] for i in keep]
    c = list(cost) + [ZERO] * (total - n)
    obj = [c[i] for i in range(total)] + [ZERO]
    for i, bi in enumerate(basis):
        if c[bi] != 0:
            f, Ti = c[bi], T[i]
            for j in range(total + 1):
                if Ti[j] != 0:
                    obj[j] -= f * Ti[j]
    T[-1] = obj
    status, it2 = _simplex(T, basis, nstruct, total, iteration_cap)
    iters += it2
    if status != 'optimal':
        return dict(status='phase2_' + status, value=None, z=None, iterations=iters)
    z = [ZERO] * n
    for i, bi in enumerate(basis):
        if bi < n:
            z[bi] = T[i][total]
    return dict(status='optimal', value=-T[-1][total], z=z, iterations=iters)


def _pivot(T, basis, row, col, total):
    Tr = T[row]
    piv = Tr[col]
    if piv != 1:
        T[row] = Tr = [x / piv if x != 0 else ZERO for x in Tr]
    for i in range(len(T)):
        if i == row:
            continue
        Ti = T[i]
        f = Ti[col]
        if f == 0:
            continue
        T[i] = [a - f * b if b != 0 else a for a, b in zip(Ti, Tr)]
    basis[row] = col


def _simplex(T, basis, nstruct, total, iteration_cap):
    nrow = len(T) - 1
    for it in range(iteration_cap):
        obj = T[-1]
        col = -1
        if it < BLAND_AFTER:
            best = ZERO
            for j in range(nstruct):
                if obj[j] < best:
                    best, col = obj[j], j
        else:
            for j in range(nstruct):
                if obj[j] < 0:
                    col = j
                    break
        if col < 0:
            return 'optimal', it
        row, best_ratio = -1, None
        for i in range(nrow):
            a = T[i][col]
            if a > 0:
                ratio = T[i][total] / a
                if (best_ratio is None or ratio < best_ratio
                        or (ratio == best_ratio and basis[i] < basis[row])):
                    best_ratio, row = ratio, i
        if row < 0:
            return 'unbounded', it
        _pivot(T, basis, row, col, total)
    return 'iteration_cap', iteration_cap


# ------------------------------------------------------- instance preparation

def leaf_vectors(edges, leaves):
    """Dv[v] = tuple of exact d(v, r) aligned with `leaves`, for every vertex."""
    byleaf = common.bfs_from(edges, leaves)
    verts = sorted({x for u, v, _ in edges for x in (u, v)})
    return {v: tuple(byleaf[r][v] for r in leaves) for v in verts}


def support_cols(Dv, edges, support):
    """Column vectors (alpha_1, beta_1, ..., alpha_k, beta_k) of the affine map."""
    cols = []
    for e in support:
        u, v, _ = edges[e]
        cols.append(Dv[u])
        cols.append(Dv[v])
    return cols


# ------------------------------------------------- theta(S): min-max deviation

def theta_support(cols, dvals, k, cutoff=None, prune_on='ge', seed=(), max_rounds=400):
    """Exact min over the closed set of max_r |pred_r - d_r| for a fixed support.

    Row generation: the LP is solved on an active subset of the 2|R| tolerance
    rows.  That is a RELAXATION, so its value is an exact lower bound on the
    full value and the incumbent `cutoff` prunes as soon as the bound reaches
    it (prune_on='ge') or exceeds it (prune_on='gt').  Optimality is certified
    by recomputing every leaf deviation exactly.
    """
    n = 2 * k + 2
    cost = [ZERO] * n
    cost[n - 1] = ONE
    nleaves = len(dvals)
    active = [j for j in seed if 0 <= j < nleaves]
    if not active:
        active = sorted({0, nleaves // 2, nleaves - 1})
    rounds = 0
    while True:
        rounds += 1
        rows = [[ONE] * (2 * k) + [ZERO, ZERO]]
        rhs = [ONE]
        rels = ['eq']
        for j in active:
            cr = [col[j] for col in cols]
            rows.append([-x for x in cr] + [-ONE, ONE])
            rhs.append(-dvals[j])
            rels.append('ge')
            rows.append(cr + [ONE, ONE])
            rhs.append(dvals[j])
            rels.append('ge')
        out = lp_minimize(n, cost, rows, rhs, rels)
        if out['status'] != 'optimal':
            return dict(status=out['status'], value=None, z=None, rounds=rounds,
                        active=active, lower_bound=None)
        value, z = out['value'], out['z']
        if cutoff is not None and (value >= cutoff if prune_on == 'ge' else value > cutoff):
            return dict(status='pruned', value=None, z=z, rounds=rounds,
                        active=active, lower_bound=value)
        y = z[2 * k]
        worst, worst_j = None, None
        live = [(i, z[i]) for i in range(2 * k) if z[i] != 0]
        for j in range(nleaves):
            s = y
            for i, zi in live:
                s += zi * cols[i][j]
            dev = s - dvals[j]
            if dev < 0:
                dev = -dev
            if worst is None or dev > worst:
                worst, worst_j = dev, j
        if worst <= value:
            return dict(status='optimal', value=value, z=z, rounds=rounds,
                        active=active, lower_bound=value)
        if worst_j in active or rounds > max_rounds:
            return dict(status='row_generation_stalled', value=None, z=z,
                        rounds=rounds, active=active, lower_bound=value)
        active.append(worst_j)


def strict_feasible_at(cols, dvals, k, eps):
    """Exact: is the STRICT system (alpha, beta, y > 0) feasible at tolerance eps?

    Maximise the common slack s under alpha_i >= s, beta_i >= s, y >= s and the
    2|R| tolerance rows; the strict system is feasible iff the optimum is > 0.
    """
    n = 2 * k + 2
    cost = [ZERO] * n
    cost[n - 1] = -ONE
    rows = [[ONE] * (2 * k) + [ZERO, ZERO]]
    rhs = [ONE]
    rels = ['eq']
    for i in range(2 * k + 1):
        row = [ZERO] * n
        row[i] = ONE
        row[n - 1] = -ONE
        rows.append(row)
        rhs.append(ZERO)
        rels.append('ge')
    for j, d in enumerate(dvals):
        cr = [col[j] for col in cols]
        rows.append([-x for x in cr] + [-ONE, ZERO])
        rhs.append(-d - eps)
        rels.append('ge')
        rows.append(cr + [ONE, ZERO])
        rhs.append(d - eps)
        rels.append('ge')
    out = lp_minimize(n, cost, rows, rhs, rels)
    if out['status'] != 'optimal':
        return dict(status=out['status'], slack=None, z=None)
    return dict(status='optimal', slack=-out['value'], z=out['z'])


# --------------------------------------------- Theta_k: exhaustive support sweep

def theta_min_over_supports(Dv, edges, dvals, k, cutoff=None, hot_rows=8):
    """Theta_k, exhaustively over all C(|E|, k) supports, exact."""
    best, best_support, best_z = cutoff, None, None
    hot, visited, pruned = [], 0, 0
    for support in itertools.combinations(range(len(edges)), k):
        visited += 1
        out = theta_support(support_cols(Dv, edges, support), dvals, k,
                            cutoff=best, seed=hot[:hot_rows])
        if out['status'] == 'pruned':
            pruned += 1
            continue
        if out['status'] != 'optimal':
            return dict(status=out['status'], support=support, visited=visited)
        hot = list(dict.fromkeys(out['active'] + hot))[:hot_rows]
        if best is None or out['value'] < best:
            best, best_support, best_z = out['value'], support, out['z']
    return dict(status='complete', theta=best, support=best_support, z=best_z,
                visited=visited, pruned=pruned)


def first_strict_zero(Dv, edges, dvals, k, hot_rows=8):
    """First support of size k that is STRICTLY feasible at exact equality.

    theta(S) = 0 is necessary, so the sweep prunes on the exact bound
    theta(S) > 0 and the max-slack LP is run only on survivors.
    """
    hot, tested, zero = [], 0, 0
    for support in itertools.combinations(range(len(edges)), k):
        tested += 1
        cols = support_cols(Dv, edges, support)
        out = theta_support(cols, dvals, k, cutoff=ZERO, prune_on='gt', seed=hot[:hot_rows])
        if out['status'] == 'pruned':
            continue
        if out['status'] != 'optimal':
            return dict(status=out['status'], support=support, tested=tested)
        hot = list(dict.fromkeys(out['active'] + hot))[:hot_rows]
        if out['value'] != 0:
            continue
        zero += 1
        sf = strict_feasible_at(cols, dvals, k, ZERO)
        if sf['status'] != 'optimal':
            return dict(status='slack_' + sf['status'], support=support, tested=tested)
        if sf['slack'] > 0:
            return dict(status='found', support=support, z=sf['z'], slack=sf['slack'],
                        tested=tested, zero_theta_seen=zero)
    return dict(status='none', tested=tested, zero_theta_seen=zero)


# ------------------------------------------------------- independent checking

def witness_sources(support, z, k):
    """LP endpoint masses -> (edge, p, x) sources; p may be 0 and x may hit 0 or 1."""
    out = []
    for i, e in enumerate(support):
        alpha, beta = z[2 * i], z[2 * i + 1]
        p = alpha + beta
        out.append(dict(edge=int(e), p=p, x=(beta / p if p != 0 else ZERO),
                        alpha=alpha, beta=beta))
    return out, z[2 * k]


def independent_deviation(edges, support, z, k, observed, leaves=None):
    """Recompute the leaf equations on the ORIGINAL tree with the frozen
    min(d(r,u) + xL, d(r,v) + (1-x)L) form - an implementation independent of
    the interpolation the LP uses - and return the exact max |pred - d|."""
    sources, y = witness_sources(support, z, k)
    live = [dict(edge=s['edge'], p=s['p'], x=s['x']) for s in sources if s['p'] != 0]
    if not live:
        return None
    targets = sorted(observed) if leaves is None else list(leaves)
    pred = common.predict(edges, live, y, targets)
    return max(abs(pred[r] - observed[r]) for r in targets)
