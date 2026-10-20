#!/usr/bin/env python3
"""Stage43 exact equality-reduced local LP, with frozen Stage41 outer search.

Equality elimination, parallel-row reduction and 1D interval intersection are
classical engineering. This module makes no new asymptotic or novelty claim.
All arithmetic is fractions.Fraction; explicit zero faces and strict slack
are preserved. Earlier stage sources remain read-only.
"""
from fractions import Fraction as F
from itertools import combinations
from stage23_two_source_envelope import two_pass_envelope, inverse
from stage25_deficit_extension import tree
from stage27_general_fpt import Budget, ResourceLimit, square_solve, dot
from stage26_missing_leaf import observed_backbone, observe
from stage41_partial_fpt import mixed_cover, deficit_sets, realize, lift_mixed


def _count(budget, key, value=1):
    budget.counts[key] = budget.counts.get(key, 0) + value
    budget.check()


def _parallel_reduce(constraints):
    """Positive normalization retains the strongest parallel inequality."""
    unique = {}
    for a, rhs in constraints:
        a, rhs = tuple(map(F, a)), F(rhs)
        first = next((x for x in a if x), None)
        if first is None:
            if rhs < 0:
                return None
            continue
        scale = abs(first)
        row, bound = tuple(x / scale for x in a), rhs / scale
        if row not in unique or bound < unique[row]:
            unique[row] = bound
    return unique


def _affine_equalities(equalities, dimension):
    """Exact RREF: x = origin + sum_j parameter[j] * basis[j]."""
    matrix = [[F(x) for x in a] + [F(b)] for a, b in equalities]
    pivots, rank = [], 0
    for column in range(dimension):
        row = next((j for j in range(rank, len(matrix)) if matrix[j][column]), None)
        if row is None:
            continue
        matrix[rank], matrix[row] = matrix[row], matrix[rank]
        scale = matrix[rank][column]
        matrix[rank] = [x / scale for x in matrix[rank]]
        for j in range(len(matrix)):
            if j != rank and matrix[j][column]:
                scale = matrix[j][column]
                matrix[j] = [x - scale * y for x, y in zip(matrix[j], matrix[rank])]
        pivots.append(column)
        rank += 1
    if any(not any(row[:-1]) and row[-1] for row in matrix):
        return None
    free = [j for j in range(dimension) if j not in pivots]
    origin = [F(0)] * dimension
    for i, column in enumerate(pivots):
        origin[column] = matrix[i][-1]
    basis = []
    for column in free:
        vector = [F(0)] * dimension
        vector[column] = F(1)
        for i, pivot in enumerate(pivots):
            vector[pivot] = -matrix[i][column]
        basis.append(tuple(vector))
    return tuple(origin), basis, rank


def strict_vertex(constraints, dimension, budget):
    """Exact strict feasibility of the bounded Stage41 common-slack LP.

    A returned point satisfies the original inequalities and t > 0. Opposing
    tight inequalities are equalities, solved first. The injective affine
    parametrization maps the remaining bounded polytope bijectively; thus a
    positive maximum of affine t is attained at a reduced-polytope vertex.
    Vertex enumeration remains complete, including lower-dimensional faces.
    """
    _count(budget, 'local_lp_calls')
    rows = _parallel_reduce(constraints)
    if rows is None:
        _count(budget, 'constant_row_rejections')
        return None
    equalities, consumed = [], set()
    for row, bound in rows.items():
        if row in consumed:
            continue
        opposite = tuple(-x for x in row)
        if opposite in rows:
            if bound < -rows[opposite]:
                _count(budget, 'opposing_row_rejections')
                return None
            if bound == -rows[opposite]:
                equalities.append((row, bound))
                consumed.update((row, opposite))
    affine = _affine_equalities(equalities, dimension)
    if affine is None:
        _count(budget, 'inconsistent_equality_rejections')
        return None
    origin, basis, rank = affine
    _count(budget, 'eliminated_equality_rank', rank)
    reduced = _parallel_reduce([(tuple(dot(a, vector) for vector in basis), rhs - dot(a, origin))
                               for a, rhs in rows.items() if a not in consumed])
    if reduced is None:
        _count(budget, 'constant_row_rejections')
        return None
    reduced = list(reduced.items())
    def lift(parameters):
        point = tuple(origin[i] + sum((z * vector[i] for z, vector in zip(parameters, basis)), F(0))
                      for i in range(dimension))
        if point[-1] <= 0:
            return None
        assert all(dot(a, point) <= rhs for a, rhs in constraints)
        return point
    if not basis:
        _count(budget, 'zero_dimensional_lp_calls')
        return lift(())
    if len(basis) == 1:
        _count(budget, 'interval_lp_calls')
        lower = max((rhs / a[0] for a, rhs in reduced if a[0] < 0), default=None)
        upper = min((rhs / a[0] for a, rhs in reduced if a[0] > 0), default=None)
        assert lower is not None and upper is not None, 'Expected bounded common-slack polytope'
        if lower > upper:
            return None
        chosen = upper if basis[0][-1] >= 0 else lower
        return lift((chosen,))
    _count(budget, 'enumerated_lp_calls')
    for chosen in combinations(range(len(reduced)), len(basis)):
        budget.touch('vertex_bases')
        z = square_solve([reduced[i][0] for i in chosen], [reduced[i][1] for i in chosen])
        if z is not None and origin[-1] + sum((a * v[-1] for a, v in zip(z, basis)), F(0)) > 0 \
                and all(dot(a, z) <= rhs for a, rhs in reduced):
            return lift(z)
    return None


def solve_mixed(edges, observed, marked=(), cap=3, seconds=30, max_states=200000,
                max_bases=2000000, max_faces=500000):
    """All leaves of a reduced tree observed; marked internal vertex atoms legal."""
    if not isinstance(cap, int) or isinstance(cap, bool) or cap < 1:
        raise ValueError('Positive integer source budget required')
    budget = Budget(seconds, max_states, max_bases, max_faces)
    stats = {'deficit_sets': 0, 'cover_feasible_faces': 0, 'strict_faces': 0}
    best = None
    def report(reason=None):
        answer = {'status': 'complete', 'feasible': best is not None,
                  'statistics': dict(stats, **budget.counts)}
        if reason:
            answer['reason'] = reason
        if best is not None:
            answer.update(best, minimum_sources_within_cap=best['source_count'])
        return answer
    try:
        adj, lookup = tree(edges)
        leaves = {v for v in adj if len(adj[v]) == 1}
        marked = set(marked)
        if set(observed) != leaves or not marked <= set(adj) - leaves:
            raise ValueError('All backbone leaves and only internal vertex marks required')
        obs = {v: F(d) for v, d in observed.items()}
        h, _ = two_pass_envelope(adj, obs)
        if any(h[r] != obs[r] for r in leaves):
            return report('inconsistent_leaf_envelope')
        b, Y = inverse(adj, h)
        if Y <= 0:
            return report('nonpositive_baseline_pendant')
        positive = {v for v in adj if b[v] > 0}
        negative = {v for v in adj if b[v] < 0}
        if len(positive) > (2 * cap if cap <= 2 else 3 * cap - 2) or len(negative) > max(0, 4 * cap - 8):
            return report('baseline_mass_count_obstruction')
        eligible = {v for v in adj if 3 <= len(adj[v]) <= cap}
        for S in deficit_sets(adj, positive & eligible, eligible, max(0, cap - 2), budget):
            stats['deficit_sets'] += 1
            A = set(S)
            for v in S:
                A.update(adj[v])
            if negative - A:
                continue
            fixed = positive - A
            rows = {v: [F(0)] * len(S) for v in A}
            for j, v in enumerate(S):
                rows[v][j] = -sum(1 / length for length in adj[v].values()) / 2
                for u, length in adj[v].items():
                    rows[u][j] += 1 / (2 * length)
            ordered = sorted(A)
            # A mask explicitly chooses the positive affected coordinates;
            # the complement is an equality face, including boundary faces.
            for mask in range(1 << len(ordered)):
                budget.touch('mass_faces')
                P = fixed | {v for j, v in enumerate(ordered) if mask & (1 << j)}
                cover = mixed_cover(P, marked, lookup, edges, cap)
                if cover is None or best is not None and cover['count'] >= best['source_count']:
                    continue
                stats['cover_feasible_faces'] += 1
                d = len(S) + 1
                lower_t = (F(0),) * len(S) + (F(-1),)
                upper_t = (F(0),) * len(S) + (F(1),)
                constraints = [(lower_t, F(0)), (upper_t, F(1))]
                for j in range(len(S)):
                    row = [-F(i == j) for i in range(len(S))] + [F(1)]
                    constraints.append((tuple(row), F(0)))
                constraints.append((tuple(F(len(adj[v]) - 2, 2) for v in S) + (F(1),), Y))
                for v in ordered:
                    if v in P:
                        constraints.append((tuple(-a for a in rows[v]) + (F(1),), b[v]))
                    else:
                        constraints.extend([(tuple(rows[v]) + (F(0),), -b[v]),
                                            (tuple(-a for a in rows[v]) + (F(0),), b[v])])
                x = strict_vertex(constraints, d, budget)
                if x is None:
                    continue
                stats['strict_faces'] += 1
                delta = x[:-1]
                mass = {v: b[v] + (dot(rows[v], delta) if v in A else 0) for v in P}
                y = Y - sum(F(len(adj[v]) - 2, 2) * z for v, z in zip(S, delta))
                best = realize(edges, mass, cover, y, S, delta)
                if best['source_count'] == 1:
                    return report()
        return report()
    except ResourceLimit as exc:
        return {'status': 'unknown', 'reason': str(exc), 'incumbent': best,
                'feasibility_established': best is not None,
                'statistics': dict(stats, **budget.counts)}


def solve_partial(edges, observed, cap=3, **options):
    """At least two observed original leaves, free positive mean pendant."""
    if len(observed) < 2:
        raise ValueError('Prototype requires at least two observed leaves')
    observed = {v: F(d) for v, d in observed.items()}
    backbone = observed_backbone(edges, observed)
    marked = {x['vertex'] for x in backbone['saturation_violations']}
    answer = solve_mixed(backbone['edges'], observed, marked, cap, **options)
    answer['marked_vertices'] = sorted(marked)
    answer['backbone_edges'] = backbone['edges']
    mixed = answer if answer.get('feasible') else answer.get('incumbent')
    if mixed is not None:
        lifted = lift_mixed(backbone, mixed['atoms'], mixed['pendant'])
        assert observe(edges, lifted['sources'], lifted['pendant'], observed) == observed
        answer['original_witness'] = lifted
    return answer
