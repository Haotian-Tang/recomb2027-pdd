#!/usr/bin/env python3
"""Exact partial-leaf PDD via marked projection gates.

Stage41 extension; older writers are not invoked. The search implements the
Stage40 slack/zero-face proof with Stage27 anchored frontier enumeration.
Classical matching and exact active bases are tools, not novelty claims.
Resource stops return unknown. This prototype is for small-budget audits.
"""
from fractions import Fraction as F
from itertools import combinations
from stage23_two_source_envelope import two_pass_envelope, inverse
from stage25_deficit_extension import tree
from stage27_general_fpt import Budget, ResourceLimit, minimum_forest_cover, square_solve, dot
from stage26_missing_leaf import observed_backbone, observe


def mixed_cover(positive, marked, lookup, edges, cap):
    P = set(positive)
    if not P or len(P) > 2 * cap:
        return None
    paired = set()
    for u, v in combinations(sorted(P), 2):
        if (u, v) in lookup:
            paired.update((u, v))
    isolated = P - paired
    if not isolated <= marked or len(isolated) > cap:
        return None
    selected = minimum_forest_cover(paired, lookup, edges, cap - len(isolated)) if paired else []
    if selected is None:
        return None
    return {'edges': selected, 'vertices': sorted(isolated), 'count': len(selected) + len(isolated)}


def deficit_sets(adj, anchors, eligible, maxweight, budget):
    """Frozen Stage27 include/exclude logic, extracted without side effects."""
    for nr in range(min(len(anchors), maxweight) + 1):
        for seeds in combinations(sorted(anchors), nr):
            budget.check()
            current = frozenset(seeds)
            weight = sum(len(adj[v]) - 2 for v in current)
            if weight > maxweight:
                continue
            rejected = frozenset(anchors - current)
            frontier = frozenset(u for v in current for u in adj[v]
                                 if u in eligible and u not in current and u not in rejected)
            stack = [(current, rejected, frontier, weight)]
            while stack:
                budget.check()
                current, rejected, frontier, weight = stack.pop()
                if not frontier:
                    budget.touch('growth_states')
                    yield tuple(sorted(current))
                    continue
                v = min(frontier)
                remaining = frontier - {v}
                stack.append((current, rejected | {v}, remaining, weight))
                nw = weight + len(adj[v]) - 2
                if nw <= maxweight:
                    included = current | {v}
                    expanded = remaining | frozenset(u for u in adj[v]
                        if u in eligible and u not in included and u not in rejected)
                    stack.append((included, rejected, expanded, nw))


def strict_vertex(constraints, dimension, budget):
    """Find a rational feasible active-basis vertex with last coordinate >0."""
    for chosen in combinations(range(len(constraints)), dimension):
        budget.touch('vertex_bases')
        x = square_solve([constraints[i][0] for i in chosen],
                         [constraints[i][1] for i in chosen])
        if x is not None and x[-1] > 0 and all(dot(a, x) <= rhs for a, rhs in constraints):
            return x
    return None


def realize(edges, mass, cover, pendant, deficit_vertices, deficits):
    selected = cover['edges']
    degree = {v: sum(v in edges[i][:2] for i in selected) for v in mass}
    atoms = [{'kind': 'vertex', 'vertex': v, 'p': mass[v]} for v in cover['vertices']]
    for i in selected:
        u, v, _ = edges[i]
        a, b = mass[u] / degree[u], mass[v] / degree[v]
        atoms.append({'kind': 'edge', 'edge': i, 'p': a + b, 'x': b / (a + b)})
    assert sum(a['p'] for a in atoms) == 1 and pendant > 0
    return {'atoms': atoms, 'pendant': pendant, 'source_count': len(atoms),
            'positive_mass': mass, 'deficit_vertices': list(deficit_vertices), 'deficits': list(deficits)}


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


def lift_mixed(backbone, atoms, pendant):
    """Lift each distinct backbone atom to a distinct original open edge."""
    y, cost, sources = F(pendant), F(0), []
    E = backbone['original_edges']
    for atom in atoms:
        p = atom['p']
        if atom['kind'] == 'vertex':
            gate = atom['vertex']
            edge = backbone['gate_edges'][gate][0]
            u, v, length = E[edge]
            distance = min(length / 2, y / 2)
            x = distance / length if u == gate else 1 - distance / length
            cost += p * distance
        else:
            coordinate = atom['x'] * backbone['edges'][atom['edge']][2]
            for part in backbone['chains'][atom['edge']]:
                if part['start'] < coordinate < part['end']:
                    edge = part['edge']
                    x = (coordinate - part['start']) / part['length']
                    if E[edge][0] != part['u']:
                        x = 1 - x
                    break
                if coordinate == part['end']:
                    gate = part['v']
                    edge = backbone['gate_edges'][gate][0]
                    u, v, length = E[edge]
                    distance = min(length / 2, y / 2)
                    x = distance / length if u == gate else 1 - distance / length
                    cost += p * distance
                    break
            else:
                raise AssertionError('Strict backbone point not found in original chain')
        sources.append({'edge': edge, 'p': p, 'x': x})
    assert len({s['edge'] for s in sources}) == len(sources)
    assert y - cost >= y / 2 > 0
    return {'sources': sources, 'pendant': y - cost, 'source_count': len(sources)}


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
