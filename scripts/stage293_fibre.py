#!/usr/bin/env python3
"""Stage 293: exact rational analysis of a fixed-support fibre in the PDD model.

For a fixed support S the model is affine in z = (p_0..p_{k-1}, t_0..t_{k-1}, y)
with t_i = p_i x_i (stage289_support_feasibility.build_system).  The set of
explanations on S -- the FIBRE -- is

    Phi(S) = { z : A z = b,  y > 0,  p_i > 0,  t_i > 0,  p_i - t_i > 0 }.

This module measures how large Phi(S) is, in exact rational arithmetic.

Geometry used throughout
------------------------
Write Phi_closed for the same system with all strict inequalities relaxed to
">= 0".  If Phi(S) is non-empty then Phi_closed = closure(Phi(S)): for z in
Phi_closed and z0 in Phi(S) every point (1-th) z + th z0, th in (0,1], is
strictly feasible, so z is a limit of points of Phi(S).  Hence the supremum of
any continuous functional over Phi(S) equals its supremum over Phi_closed, which
is what the linear programs below return.  Extremes are suprema: they may be
approached but not attained inside the open set.

Phi_closed is BOUNDED: p_i >= 0 with sum_i p_i = 1 bounds p, 0 <= t_i <= p_i
bounds t, and 0 <= y <= min_r d_r bounds y.  No range is infinite.

Substituting p_i = t_i + s_i with s_i >= 0 puts the system in standard form
    max c'w  s.t.  M w = b,  w >= 0,     w = (t_0..t_{k-1}, s_0..s_{k-1}, y)
solved by a Fraction simplex with Bland's rule (finite, no cycling, no rounding).

The position range needs x_i = t_i / p_i, a linear-fractional objective; it is
handled by the Charnes-Cooper homogenisation (_scale_free_max), exact as well.
"""
from fractions import Fraction as F
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from stage289_support_feasibility import build_system, rref


# --------------------------------------------------------------- exact simplex

class LPError(Exception):
    pass


def _pivot(tab, basis, row, col):
    scale = tab[row][col]
    tab[row] = [v / scale for v in tab[row]]
    for j in range(len(tab)):
        if j != row and tab[j][col]:
            factor = tab[j][col]
            tab[j] = [a - factor * b for a, b in zip(tab[j], tab[row])]
    basis[row] = col


def _simplex(tab, basis, n_cols, allowed):
    """Phase-2 simplex on an explicit tableau; Bland's rule. Maximise."""
    guard = 0
    while True:
        guard += 1
        if guard > 20000:
            raise LPError('iteration_cap')
        col = None
        for j in range(n_cols):
            if j in allowed and tab[-1][j] < 0:
                col = j
                break
        if col is None:
            return
        row, best = None, None
        for i in range(len(tab) - 1):
            if tab[i][col] > 0:
                ratio = tab[i][-1] / tab[i][col]
                if best is None or ratio < best or (ratio == best and basis[i] < basis[row]):
                    row, best = i, ratio
        if row is None:
            raise LPError('unbounded')
        _pivot(tab, basis, row, col)


def lp_max(matrix, rhs, cost):
    """Exact max cost'w s.t. matrix w = rhs, w >= 0.  Returns (value, w)."""
    m, n = len(matrix), len(cost)
    rows = []
    for row, value in zip(matrix, rhs):
        row, value = list(row), F(value)
        if value < 0:
            row, value = [-v for v in row], -value
        rows.append(row + [value])
    # phase 1: one artificial variable per row
    tab = [rows[i][:n] + [F(1) if j == i else F(0) for j in range(m)] + [rows[i][n]]
           for i in range(m)]
    basis = [n + i for i in range(m)]
    objective = [F(0)] * (n + m + 1)
    for i in range(m):
        objective = [a + b for a, b in zip(objective, tab[i])]
    tab.append([-v for v in objective])
    _simplex(tab, basis, n + m, set(range(n)))
    if -tab[-1][-1] != 0:
        raise LPError('infeasible')
    # drive artificials out of the basis where possible
    for i in range(m):
        if basis[i] >= n:
            col = next((j for j in range(n) if tab[i][j]), None)
            if col is not None:
                _pivot(tab, basis, i, col)
    keep = [i for i in range(m) if basis[i] < n]
    tab = [tab[i] for i in keep] + [tab[-1]]
    basis = [basis[i] for i in keep]
    # phase 2
    reduced = [-c for c in cost] + [F(0)] * m + [F(0)]
    for i, b in enumerate(basis):
        if reduced[b]:
            factor = reduced[b]
            reduced = [a - factor * v for a, v in zip(reduced, tab[i])]
    tab[-1] = reduced
    _simplex(tab, basis, n, set(range(n)))
    w = [F(0)] * n
    for i, b in enumerate(basis):
        if b < n:
            w[b] = tab[i][-1]
    return tab[-1][-1], w


def _reduced_equalities(edges, observed, support):
    """RREF of the leaf equations: at most 2k+1 independent rows."""
    rows, width = build_system(edges, observed, support)
    reduced, pivots = rref(rows, width)
    kept = [row for row in reduced if any(row[:-1])]
    inconsistent = any(not any(row[:-1]) and row[-1] for row in reduced)
    return kept, pivots, width, inconsistent


def _standard_form(kept, width, k):
    """Columns for w = (t_0..t_{k-1}, s_0..s_{k-1}, y) with p_i = t_i + s_i."""
    matrix, rhs = [], []
    for row in kept:
        new = [row[i] + row[k + i] for i in range(k)] + [row[i] for i in range(k)] + [row[2 * k]]
        matrix.append(new)
        rhs.append(row[-1])
    return matrix, rhs


def _coord_cost(kind, index, k):
    """Objective in w-space for p_index / t_index / y."""
    cost = [F(0)] * (2 * k + 1)
    if kind == 'p':
        cost[index] = F(1)
        cost[k + index] = F(1)
    elif kind == 't':
        cost[index] = F(1)
    elif kind == 'y':
        cost[2 * k] = F(1)
    return cost


def _scale_free_max(matrix, rhs, k, index, sign):
    """Exact sup of sign * t_index / p_index by Charnes-Cooper homogenisation.

    w = what / sigma with sigma > 0; add p_index(what) = 1 and maximise t_index.
    """
    n = 2 * k + 1
    rows = [list(row) + [-value] for row, value in zip(matrix, rhs)]
    norm = [F(0)] * n + [F(0)]
    norm[index] = F(1)
    norm[k + index] = F(1)
    rows.append(norm)
    rhs2 = [F(0)] * len(matrix) + [F(1)]
    cost = [F(0)] * (n + 1)
    cost[index] = F(sign)
    value, _ = lp_max(rows, rhs2, cost)
    return value * sign


def fibre_report(edges, observed, support, truth=None, want_x=True):
    """Exact size of the fibre over `support` given `observed`."""
    support = list(support)
    k = len(support)
    if len(set(support)) != k:
        return dict(status='invalid', reason='repeated_edge')
    kept, pivots, width, inconsistent = _reduced_equalities(edges, observed, support)
    if inconsistent:
        return dict(status='infeasible', reason='inconsistent_equalities', rank=len(pivots))
    matrix, rhs = _standard_form(kept, width, k)
    out = dict(status='ok', rank=len(pivots), width=width,
               free_parameters=width - len(pivots), n_equations=len(kept))
    try:
        ranges = {}
        for kind, count in (('p', k), ('y', 1)):
            for i in range(count):
                hi, point_hi = lp_max(matrix, rhs, _coord_cost(kind, i, k))
                lo, point_lo = lp_max(matrix, rhs, [-v for v in _coord_cost(kind, i, k)])
                ranges[(kind, i)] = (-lo, hi, point_lo, point_hi)
    except LPError as exc:
        return dict(out, status='lp_failed', reason=str(exc))
    if any(hi < lo for lo, hi, _, _ in ranges.values()):
        return dict(out, status='lp_failed', reason='crossed_bounds')
    p_ranges = [(ranges[('p', i)][0], ranges[('p', i)][1]) for i in range(k)]
    y_lo, y_hi = ranges[('y', 0)][0], ranges[('y', 0)][1]
    out['p_bounds'] = [[str(lo), str(hi)] for lo, hi in p_ranges]
    out['y_bounds'] = [str(y_lo), str(y_hi)]
    out['p_sup_range'] = str(max((hi - lo for lo, hi in p_ranges), default=F(0)))
    out['y_sup_range'] = str(y_hi - y_lo)
    out['p_range_argmax_index'] = (max(range(k), key=lambda i: p_ranges[i][1] - p_ranges[i][0])
                                   if k else None)
    x_bounds = None
    if want_x:
        try:
            x_bounds = []
            for i in range(k):
                hi = _scale_free_max(matrix, rhs, k, i, 1)
                lo = _scale_free_max(matrix, rhs, k, i, -1)
                x_bounds.append((lo, hi))
        except LPError as exc:
            return dict(out, status='lp_failed', reason='x_' + str(exc))
        out['x_bounds'] = [[str(lo), str(hi)] for lo, hi in x_bounds]
        out['x_sup_range'] = str(max((hi - lo for lo, hi in x_bounds), default=F(0)))
    if truth is not None:
        tp = {s['edge']: F(s['p']) for s in truth}
        tx = {s['edge']: F(s['x']) for s in truth}
        if set(tp) == set(support):
            dev_p, dev_x, inside = [], [], True
            for i, edge in enumerate(support):
                lo, hi = p_ranges[i]
                dev_p.append(max(hi - tp[edge], tp[edge] - lo))
                if not (lo <= tp[edge] <= hi):
                    inside = False
                if x_bounds is not None:
                    xlo, xhi = x_bounds[i]
                    dev_x.append(max(xhi - tx[edge], tx[edge] - xlo))
                    if not (xlo <= tx[edge] <= xhi):
                        inside = False
            out['truth_inside_bounds'] = inside
            out['max_p_deviation_from_truth'] = str(max(dev_p))
            out['p_deviation_from_truth'] = [str(v) for v in dev_p]
            if dev_x:
                out['max_x_deviation_from_truth'] = str(max(dev_x))
                out['x_deviation_from_truth'] = [str(v) for v in dev_x]
    if k and out['p_range_argmax_index'] is not None:
        i = out['p_range_argmax_index']
        _, _, point_lo, point_hi = ranges[('p', i)]
        out['extreme_pair'] = [_to_sources(support, point, k) for point in (point_lo, point_hi)]
    return out


def _to_sources(support, w, k):
    """w = (t, s, y) -> {edge, p, x} with p = t + s; x undefined when p = 0."""
    out = []
    for i, edge in enumerate(support):
        t, s = w[i], w[k + i]
        p = t + s
        out.append(dict(edge=edge, p=str(p), x=(str(t / p) if p else None)))
    return dict(sources=out, pendant=str(w[2 * k]))
