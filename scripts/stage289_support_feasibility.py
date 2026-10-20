#!/usr/bin/env python3
"""Stage 289: exact fixed-support feasibility for the PDD model.

For a fixed set S of distinct original edges the model

    d_r = y + sum_i p_i * min(d(r,u_i) + x_i L_i, d(r,v_i) + (1 - x_i) L_i)

is affine in (p_i, t_i = p_i x_i, y): which endpoint of edge i a leaf r reaches
first is a property of the tree, not of x_i.  Deciding whether S admits an exact
explanation is therefore rational linear feasibility, solved here in Fractions.
A fixed-support decision is NOT a minimality claim; it only says whether this
particular support can reproduce the observations exactly.
"""
from fractions import Fraction as F
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from stage289_decodiphy_bridge import geometry, vertex_distances


def rref(matrix, width):
    """Exact reduced row echelon form; returns (rows, pivots)."""
    rows = [list(r) for r in matrix]
    pivots, rank = [], 0
    for column in range(width):
        pick = next((j for j in range(rank, len(rows)) if rows[j][column]), None)
        if pick is None:
            continue
        rows[rank], rows[pick] = rows[pick], rows[rank]
        scale = rows[rank][column]
        rows[rank] = [v / scale for v in rows[rank]]
        for j in range(len(rows)):
            if j != rank and rows[j][column]:
                factor = rows[j][column]
                rows[j] = [a - factor * b for a, b in zip(rows[j], rows[rank])]
        pivots.append(column)
        rank += 1
    return rows, pivots


def build_system(edges, observed, support):
    """Equality rows for variables (p_0..p_{k-1}, t_0..t_{k-1}, y)."""
    k = len(support)
    width = 2 * k + 1
    table = {}
    for index in support:
        for v in edges[index][:2]:
            if v not in table:
                table[v] = vertex_distances(edges, v)
    rows = []
    for r, value in sorted(observed.items()):
        row = [F(0)] * width + [F(value)]
        row[2 * k] = F(1)          # the common mean pendant y
        for i, index in enumerate(support):
            u, v, length = edges[index]
            du, dv = table[u][r], table[v][r]
            if dv - du == length:            # leaf lies beyond u
                row[i] += du
                row[k + i] += length
            elif du - dv == length:          # leaf lies beyond v
                row[i] += dv + length
                row[k + i] -= length
            else:
                raise ValueError('leaf is not separated by the edge')
        rows.append(row)
    rows.append([F(1)] * k + [F(0)] * (k + 1) + [F(1)])
    return rows, width


def feasible_support(edges, observed, support):
    """Exact verdict for one fixed support.

    Returns verdict feasible / infeasible / unresolved, and whether the support
    is identifiable (the observations pin down a unique (p, t, y)).
    """
    if len(set(support)) != len(support):
        return dict(verdict='invalid', reason='repeated_edge')
    rows, width = build_system(edges, observed, support)
    reduced, pivots = rref(rows, width)
    if any(not any(row[:-1]) and row[-1] for row in reduced):
        return dict(verdict='infeasible', reason='inconsistent_equalities',
                    rank=len(pivots), identifiable=len(pivots) == width)
    k = len(support)
    origin, basis = nullspace(reduced, pivots, width)
    functionals = strict_functionals(k, width)
    identifiable = not basis
    if identifiable:
        point = origin
        violations = sorted({name.rsplit('_', 1)[0] for name, g in functionals
                             if _dot(g, point) <= 0})
        if violations:
            return dict(verdict='infeasible', reason='+'.join(violations),
                        rank=len(pivots), identifiable=True)
    else:
        point = _resolve_underdetermined(origin, basis, functionals)
        if point is None:
            return dict(verdict='unresolved', reason='no_strict_point_found',
                        rank=len(pivots), identifiable=False, free_parameters=len(basis))
    p_values = point[:k]
    t_values = point[k:2 * k]
    return dict(verdict='feasible', rank=len(pivots), identifiable=identifiable,
                free_parameters=len(basis),
                sources=[dict(edge=index, p=str(pi), x=str(ti / pi))
                         for index, pi, ti in zip(support, p_values, t_values)],
                pendant=str(point[-1]))

def nullspace(reduced, pivots, width):
    """Exact particular solution and nullspace basis of the reduced system."""
    free = [j for j in range(width) if j not in pivots]
    origin = [F(0)] * width
    for row, column in zip(reduced, pivots):
        origin[column] = row[-1]
    basis = []
    for j in free:
        vector = [F(0)] * width
        vector[j] = F(1)
        for row, column in zip(reduced, pivots):
            vector[column] = -row[j]
        basis.append(vector)
    return origin, basis


def strict_functionals(k, width):
    """Rows g with the requirement g . z > 0: y, every p_i, every t_i, every p_i - t_i."""
    rows = []
    row = [F(0)] * width
    row[2 * k] = F(1)
    rows.append(('pendant', row))
    for i in range(k):
        row = [F(0)] * width
        row[i] = F(1)
        rows.append(('weight_%d' % i, row))
        row = [F(0)] * width
        row[k + i] = F(1)
        rows.append(('position_lower_%d' % i, row))
        row = [F(0)] * width
        row[i], row[k + i] = F(1), F(-1)
        rows.append(('position_upper_%d' % i, row))
    return rows


def _dot(a, b):
    return sum((x * y for x, y in zip(a, b)), F(0))


def _resolve_underdetermined(origin, basis, functionals):
    """Look for a strictly feasible point in the affine set; certify it exactly."""
    try:
        import numpy as np
        from scipy.optimize import linprog
    except ImportError:
        return None
    constant = [float(_dot(g, origin)) for _, g in functionals]
    gradient = [[float(_dot(g, vector)) for vector in basis] for _, g in functionals]
    m = len(basis)
    # maximise s subject to  gradient . lam + constant >= s,  s <= 1
    A = np.array([[-value for value in row] + [1.0] for row in gradient] + [[0.0] * m + [1.0]])
    b = np.array(constant + [1.0])
    cost = np.array([0.0] * m + [-1.0])
    answer = linprog(cost, A_ub=A, b_ub=b, bounds=[(None, None)] * (m + 1), method='highs')
    if not answer.success or -answer.fun <= 0:
        return None
    lam = [F(float(value)) for value in answer.x[:m]]
    point = [origin[i] + sum((z * vector[i] for z, vector in zip(lam, basis)), F(0))
             for i in range(len(origin))]
    if all(_dot(g, point) > 0 for _, g in functionals):
        return point
    return None
