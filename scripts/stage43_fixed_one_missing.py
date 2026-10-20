#!/usr/bin/env python3
"""Exact strict PDD: prescribed pendant and exactly one missing original leaf.

New Stage43 extension. Retains the possibly degree-two projection gate, tests
the strict hidden-edge capacity, and permits a mandatory (not isolate-only)
gate singleton. Rational prototype, explicit resource UNKNOWNs. Earlier stages
are imported read-only. Generic LP and matching are established machinery.
"""
from fractions import Fraction as F
from itertools import combinations

from stage23_two_source_envelope import two_pass_envelope, inverse
from stage27_general_fpt import Budget, ResourceLimit, minimum_forest_cover, dot
from stage43_partial_fpt import strict_vertex


def tree(edges):
    """Tree geometry retaining degree-two vertices, unlike older reduced helpers."""
    adj, lookup = {}, {}
    for i, (u, v, length) in enumerate(edges):
        assert u != v and length > 0 and (min(u, v), max(u, v)) not in lookup
        adj.setdefault(u, {})[v] = F(length)
        adj.setdefault(v, {})[u] = F(length)
        lookup[min(u, v), max(u, v)] = i
    assert adj and len(edges) == len(adj) - 1
    seen = {min(adj)}
    queue = list(seen)
    for v in queue:
        for u in adj[v]:
            if u not in seen:
                seen.add(u)
                queue.append(u)
    assert len(seen) == len(adj)
    return adj, lookup


def deficit_sets(adj, anchors, eligible, weight_limit, budget):
    """Canonical anchored growth permitting one eligible zero-weight vertex."""
    for size in range(len(anchors) + 1):
        for seeds in combinations(sorted(anchors), size):
            budget.check()
            current = frozenset(seeds)
            weight = sum(len(adj[v]) - 2 for v in current)
            if weight > weight_limit:
                continue
            rejected = frozenset(anchors - current)
            frontier = frozenset(u for v in current for u in adj[v]
                                 if u in eligible and u not in current and u not in rejected)
            stack = [(current, rejected, frontier, weight)]
            while stack:
                budget.check()
                current, rejected, frontier, weight = stack.pop()
                budget.touch('growth_states')
                if not frontier:
                    yield tuple(sorted(current))
                    continue
                v = min(frontier)
                remaining = frontier - {v}
                stack.append((current, rejected | {v}, remaining, weight))
                next_weight = weight + len(adj[v]) - 2
                if next_weight <= weight_limit:
                    included = current | {v}
                    expanded = remaining | frozenset(u for u in adj[v]
                        if u in eligible and u not in included and u not in rejected)
                    stack.append((included, rejected, expanded, next_weight))


def forced_cover(positive, gate, dummy, lookup_augmented, edges_augmented, cap):
    """Attach a private leaf: its mandatory edge is the mandatory singleton."""
    if gate not in positive:
        return None
    selected = minimum_forest_cover(set(positive) | {dummy}, lookup_augmented,
                                    edges_augmented, cap)
    if selected is None:
        return None
    dummy_index = len(edges_augmented) - 1
    assert dummy_index in selected
    return {'edges': [e for e in selected if e != dummy_index], 'hidden': True,
            'count': len(selected)}


def solve_fixed_one_missing(edges, observed, pendant, cap=3, seconds=30,
                            max_states=200000, max_bases=2000000, max_faces=500000):
    if not isinstance(cap, int) or isinstance(cap, bool) or cap < 1:
        raise ValueError('Positive integer source budget required')
    y0 = F(pendant)
    if y0 <= 0:
        raise ValueError('Prescribed pendant must be positive')
    edges = [(u, v, F(length)) for u, v, length in edges]
    obs = {v: F(d) for v, d in observed.items()}
    budget = Budget(seconds, max_states, max_bases, max_faces)
    stats = {'deficit_sets': 0, 'cover_feasible_faces': 0, 'strict_faces': 0,
             'zero_weight_gate_sets': 0}
    best = None

    def report(reason=None):
        result = {'status': 'complete', 'feasible': best is not None,
                  'statistics': dict(stats, **budget.counts)}
        if reason:
            result['reason'] = reason
        if best is not None:
            result.update(best, minimum_sources_within_cap=best['source_count'])
        return result

    try:
        adj_T, _ = tree(edges)
        leaves = {v for v in adj_T if len(adj_T[v]) == 1}
        if not set(obs) < leaves or len(leaves - set(obs)) != 1:
            raise ValueError('Exactly one original leaf must be unobserved')
        if any(len(n) == 2 for n in adj_T.values()):
            raise ValueError('Original tree must be reduced')
        hidden = next(iter(leaves - set(obs)))
        gate, ell = next(iter(adj_T[hidden].items()))
        hidden_edge = next(i for i, (u, v, _) in enumerate(edges) if hidden in (u, v))
        if len(edges) == 1:
            r = next(iter(obs))
            distance = obs[r] - y0
            if 0 < distance < ell:
                x = distance / ell if edges[0][0] == r else 1 - distance / ell
                best = {'sources': [{'edge': 0, 'p': F(1), 'x': x}],
                        'pendant': y0, 'source_count': 1, 'hidden_source': True,
                        'deficit_vertices': [], 'deficits': []}
            return report('single_original_edge')
        original_index = [i for i in range(len(edges)) if i != hidden_edge]
        H = [edges[i] for i in original_index]
        adj, lookup = tree(H)
        assert set(obs) == {v for v in adj if len(adj[v]) == 1}
        dummy = min(adj_T) - 1
        H_augmented = H + [(gate, dummy, F(1))]
        lookup_augmented = dict(lookup)
        lookup_augmented[min(gate, dummy), max(gate, dummy)] = len(H)
        h, _ = two_pass_envelope(adj, obs)
        minus_upper, _ = two_pass_envelope(adj, {v: -d for v, d in obs.items()})
        upper = {v: -value for v, value in minus_upper.items()}
        if any(h[r] != obs[r] for r in obs):
            return report('inconsistent_leaf_envelope')
        b, Y = inverse(adj, h)
        if Y < y0:
            return report('projected_pendant_below_prescribed')
        positive = {v for v in adj if b[v] > 0}
        negative = {v for v in adj if b[v] < 0}
        if len(positive) > (2 if cap == 1 else 3 * cap - 1) or \
                len(negative) > (0 if cap == 1 else 4 * cap - 5):
            return report('baseline_mass_count_obstruction')
        eligible = set() if cap == 1 else {
            v for v in adj if 3 <= len(adj[v]) <= cap or (v == gate and len(adj[v]) == 2)}
        for S in deficit_sets(adj, positive & eligible, eligible, max(0, cap - 2), budget):
            stats['deficit_sets'] += 1
            stats['zero_weight_gate_sets'] += int(gate in S and len(adj[gate]) == 2)
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
            yrow = tuple(F(2 - len(adj[v]), 2) for v in S)
            crow = rows.get(gate, [F(0)] * len(S))
            for mask in range(1 << len(ordered)):
                budget.touch('mass_faces')
                P = fixed | {v for j, v in enumerate(ordered) if mask >> j & 1}
                if not P or len(P) > 2 * cap:
                    continue
                ordinary = minimum_forest_cover(P, lookup, H, cap)
                covers = []
                if ordinary is not None:
                    covers.append({'edges': ordinary, 'hidden': False, 'count': len(ordinary)})
                forced = forced_cover(P, gate, dummy, lookup_augmented, H_augmented, cap)
                if forced is not None:
                    covers.append(forced)
                for cover in covers:
                    if best is not None and cover['count'] >= best['source_count']:
                        continue
                    stats['cover_feasible_faces'] += 1
                    dimension = len(S) + 1
                    constraints = [((F(0),) * len(S) + (F(-1),), F(0)),
                                   ((F(0),) * len(S) + (F(1),), F(1))]
                    for j, v in enumerate(S):
                        unit = tuple(F(i == j) for i in range(len(S)))
                        constraints.append((tuple(-a for a in unit) + (F(1),), F(0)))
                        constraints.append((unit + (F(0),), upper[v] - h[v]))
                    for v in ordered:
                        if v in P:
                            constraints.append((tuple(-a for a in rows[v]) + (F(1),), b[v]))
                        else:
                            constraints.extend([(tuple(rows[v]) + (F(0),), -b[v]),
                                                (tuple(-a for a in rows[v]) + (F(0),), b[v])])
                    if cover['hidden']:
                        constraints.append((tuple(-a for a in yrow) + (F(1),), Y - y0))
                        constraints.append((tuple(a - ell * c for a, c in zip(yrow, crow)) +
                                            (F(1),), ell * b[gate] - Y + y0))
                    else:
                        constraints.extend([(yrow + (F(0),), y0 - Y),
                                            (tuple(-a for a in yrow) + (F(0),), Y - y0)])
                    point = strict_vertex(constraints, dimension, budget)
                    if point is None:
                        continue
                    stats['strict_faces'] += 1
                    delta = point[:-1]
                    mass = {v: b[v] + (dot(rows[v], delta) if v in A else 0) for v in P}
                    projected_y = Y + dot(yrow, delta)
                    residual = dict(mass)
                    selected = cover['edges']
                    degree = {v: sum(v in H[i][:2] for i in selected) for v in P}
                    sources = []
                    if cover['hidden']:
                        difference = projected_y - y0
                        p = mass[gate] if degree[gate] == 0 else (mass[gate] + difference / ell) / 2
                        distance = difference / p
                        x = distance / ell if edges[hidden_edge][0] == gate else 1 - distance / ell
                        sources.append({'edge': hidden_edge, 'p': p, 'x': x})
                        residual[gate] -= p
                    for i in selected:
                        u, v, _ = H[i]
                        a, c = residual[u] / degree[u], residual[v] / degree[v]
                        sources.append({'edge': original_index[i], 'p': a + c, 'x': c / (a + c)})
                    assert sum(s['p'] for s in sources) == 1
                    assert all(s['p'] > 0 and 0 < s['x'] < 1 for s in sources)
                    best = {'sources': sources, 'pendant': y0, 'source_count': len(sources),
                            'hidden_source': cover['hidden'], 'positive_mass': mass,
                            'deficit_vertices': list(S), 'deficits': list(delta),
                            'projected_pendant': projected_y, 'hidden_gate': gate}
                    if best['source_count'] == 1:
                        return report()
        return report()
    except ResourceLimit as exc:
        return {'status': 'unknown', 'reason': str(exc), 'incumbent': best,
                'feasibility_established': best is not None,
                'statistics': dict(stats, **budget.counts)}
