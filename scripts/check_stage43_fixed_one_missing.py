#!/usr/bin/env python3
"""Independent original-edge support oracle and boundary controls for Stage43.

The oracle uses mass equalities on the ORIGINAL tree with prescribed pendant,
not gate projection or deficit geometry. Vertex averages test strict positivity.
These finite exact checks are implementation evidence, not a theorem proof.
"""
from datetime import datetime, timezone
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import json
import random
import time

from stage25_benchmark import geometry, reduce_equalities, unique_solution, caterpillar, balanced
from stage27_general_fpt import Budget
from stage42_benchmark import observations, verify, serial
from stage43_fixed_one_missing import solve_fixed_one_missing, forced_cover, deficit_sets, tree

ROOT = Path(__file__).resolve().parents[1]


def endpoint_feasible(endpoints, observed, distance, pendant):
    n = len(endpoints)
    rows = [[F(1)] * n + [F(1)]]
    rows += [[distance[r][v] for v in endpoints] + [d - pendant] for r, d in observed.items()]
    reduced = reduce_equalities(rows, n)
    if reduced is None:
        return None
    rank = len(reduced)
    total = [F(0)] * n
    accepted = 0
    for columns in combinations(range(n), rank):
        basic = unique_solution([[row[j] for j in columns] + [row[-1]] for row in reduced], rank)
        if basic is None or any(v < 0 for v in basic):
            continue
        accepted += 1
        for j, value in zip(columns, basic):
            total[j] += value
        if all(v > 0 for v in total):
            point = [v / accepted for v in total]
            assert all(sum(a * v for a, v in zip(row[:-1], point)) == row[-1] for row in rows)
            return point
    return None


def oracle(case):
    edges, obs, y0 = case['edges'], case['observed'], case['pendant']
    _, _, distance = geometry(edges)
    cache = {}
    for k in range(1, case['cap'] + 1):
        for selected in combinations(range(len(edges)), k):
            U = tuple(sorted({v for i in selected for v in edges[i][:2]}))
            if U not in cache:
                cache[U] = endpoint_feasible(U, obs, distance, y0)
            point = cache[U]
            if point is None:
                continue
            mass = dict(zip(U, point))
            degree = {v: sum(v in edges[i][:2] for i in selected) for v in U}
            sources = []
            for i in selected:
                u, v, _ = edges[i]
                a, b = mass[u] / degree[u], mass[v] / degree[v]
                sources.append({'edge': i, 'p': a + b, 'x': b / (a + b)})
            return {'minimum': k, 'witness': {'sources': sources, 'pendant': y0}}
    return {'minimum': None}


def fixtures():
    out = []
    rng = random.Random(431014)
    shapes = [('star3', [(3, i, F(1)) for i in range(3)]),
              ('star4', [(4, i, F(i + 2, 3)) for i in range(4)]),
              ('caterpillar4', caterpillar(4)), ('caterpillar5', caterpillar(5)),
              ('balanced6', balanced(6))]
    for name, edges in shapes:
        adj, leaves, _ = geometry(edges)
        for rep in range(10):
            k = 1 + rep % 3
            chosen = rng.sample(range(len(edges)), k)
            weights = [rng.randint(1, 5) for _ in chosen]
            sources = [{'edge': i, 'p': F(w, sum(weights)), 'x': F(rng.randint(1, 6), 7)}
                       for i, w in zip(chosen, weights)]
            y0 = F(1 + rep % 3, 3)
            omitted = leaves[rep % len(leaves)]
            observed = observations(edges, sources, y0, [r for r in leaves if r != omitted])
            for cap in sorted({max(1, k - 1), k, min(3, k + 1)}):
                out.append({'name': f'{name}_rep{rep}_cap{cap}', 'edges': edges,
                            'observed': observed, 'pendant': y0, 'cap': cap,
                            'generating_sources': sources, 'omitted': omitted})
    for degree in (3, 4):
        edges = [(degree, r, F(1)) for r in range(degree)]
        for value in (F(3, 2), F(2), F(5, 2), F(3), F(7, 2)):
            for cap in (1, 2, 3):
                expected = 1 if F(2) < value < F(3) else \
                    (2 if degree == 3 and value == 2 and cap >= 2 else None)
                out.append({'name': f'star{degree}_boundary_d{value}_cap{cap}', 'edges': edges,
                            'observed': {r: value for r in range(degree - 1)},
                            'pendant': F(1), 'cap': cap, 'expected': expected})
    for reverse in (False, True):
        edges = [(1, 0, F(2))] if reverse else [(0, 1, F(2))]
        for value in (F(1), F(2), F(3)):
            out.append({'name': f'single_edge_reverse{reverse}_d{value}', 'edges': edges,
                        'observed': {0: value}, 'pendant': F(1), 'cap': 1,
                        'expected': 1 if value == 2 else None})
    # Force the K=2 zero-weight deficit at the unsuppressed gate.
    edges = [(3, r, F(1)) for r in range(3)]
    for orientation in range(8):
        changed = [(v, u, ell) if orientation >> j & 1 else (u, v, ell)
                   for j, (u, v, ell) in enumerate(edges)]
        out.append({'name': f'zero_weight_gate_orientation{orientation}', 'edges': changed,
                    'observed': {0: F(2), 1: F(2)}, 'pendant': F(1), 'cap': 2,
                    'expected': 2, 'require_degree_two_deficit': True})
    out.append({'name': 'inconsistent_observed_envelope', 'edges': edges,
                'observed': {0: F(1), 1: F(5)}, 'pendant': F(1), 'cap': 3,
                'expected': None})
    return out


def cover_controls():
    checks = 0
    for H, gate in [([(3, 0, F(1)), (3, 1, F(1))], 3),
                    ([(4, 0, F(1)), (4, 1, F(1)), (4, 5, F(1)),
                      (5, 2, F(1)), (5, 3, F(1))], 4)]:
        adj, lookup = tree(H)
        dummy = min(adj) - 1
        augmented = H + [(gate, dummy, F(1))]
        extended = dict(lookup)
        extended[min(gate, dummy), max(gate, dummy)] = len(H)
        V = sorted(adj)
        for mask in range(1, 1 << len(V)):
            P = {v for j, v in enumerate(V) if mask >> j & 1}
            result = forced_cover(P, gate, dummy, extended, augmented, len(P))
            expected = None
            if gate in P:
                candidates = [i for i, (u, v, _) in enumerate(H) if u in P and v in P]
                for k in range(len(candidates) + 1):
                    if any(({v for i in Fset for v in H[i][:2]} | {gate}) == P
                           for Fset in combinations(candidates, k)):
                        expected = k + 1
                        break
            assert (None if result is None else result['count']) == expected
            checks += 1
    return checks


def growth_controls():
    H = [(4, 0, F(1)), (4, 1, F(1)), (4, 5, F(1)),
         (5, 6, F(1)), (6, 2, F(1)), (6, 3, F(1))]
    adj, _ = tree(H)
    checks = 0
    for cap in (2, 3, 4, 5):
        eligible = {v for v in adj if 2 <= len(adj[v]) <= cap}
        order = sorted(eligible)
        for mask in range(1 << len(order)):
            anchors = {v for j, v in enumerate(order) if mask >> j & 1}
            actual = list(deficit_sets(adj, anchors, eligible, cap - 2, Budget()))
            expected = set()
            for smask in range(1 << len(order)):
                S = {v for j, v in enumerate(order) if smask >> j & 1}
                if sum(len(adj[v]) - 2 for v in S) > cap - 2:
                    continue
                left, valid = set(S), True
                while left:
                    seed = min(left)
                    left.remove(seed)
                    component, queue = {seed}, [seed]
                    for v in queue:
                        for u in adj[v]:
                            if u in left:
                                left.remove(u)
                                component.add(u)
                                queue.append(u)
                    valid &= bool(component & anchors)
                if valid:
                    expected.add(tuple(sorted(S)))
            assert len(actual) == len(set(actual)) and set(actual) == expected
            checks += 1
    return checks


def main():
    cases = fixtures()
    protocol = {'created_utc': datetime.now(timezone.utc).isoformat(),
                'purpose': 'Finite exact correctness and boundary controls, not a speed benchmark',
                'case_seed': 431014, 'cases': serial(cases),
                'code_sha256': {name: sha256((ROOT / 'scripts' / name).read_bytes()).hexdigest()
                    for name in ['stage43_fixed_one_missing.py', 'check_stage43_fixed_one_missing.py',
                                 'stage43_partial_fpt.py']}}
    protocol_path = ROOT / 'results/stage43_fixed_one_protocol.json'
    result_path = ROOT / 'results/stage43_fixed_one_audit.json'
    if protocol_path.exists() or result_path.exists():
        raise SystemExit('Refusing to overwrite existing Stage43 fixed-one evidence')
    protocol_path.write_text(json.dumps(protocol, indent=2) + '\n')
    started = time.perf_counter()
    cover_checks, growth_checks = cover_controls(), growth_controls()
    records, witnesses = [], 0
    for index, case in enumerate(cases):
        raw = solve_fixed_one_missing(case['edges'], case['observed'], case['pendant'],
                                      cap=case['cap'], seconds=20)
        assert raw['status'] == 'complete', (case['name'], raw)
        reference = oracle(case)
        minimum = raw.get('minimum_sources_within_cap')
        assert minimum == reference['minimum'], (case['name'], serial(raw), serial(reference))
        if 'expected' in case:
            assert minimum == case['expected'], (case['name'], minimum, case['expected'])
        if minimum is not None:
            verify(case, raw)
            verify(case, reference['witness'])
            assert raw['pendant'] == reference['witness']['pendant'] == case['pendant']
            witnesses += 2
        if case.get('require_degree_two_deficit'):
            assert raw['deficit_vertices'] == [3] and raw['deficits'][0] > 0
        records.append({'name': case['name'], 'minimum': minimum,
                        'solver': serial(raw), 'oracle': serial(reference)})
        if (index + 1) % 25 == 0:
            print(f'checked {index + 1}/{len(cases)}', flush=True)
    stopped = solve_fixed_one_missing(cases[0]['edges'], cases[0]['observed'],
                                     cases[0]['pendant'], cap=3, max_states=0)
    assert stopped['status'] == 'unknown' and 'minimum_sources_within_cap' not in stopped
    result = {'status': 'PASS', 'created_utc': datetime.now(timezone.utc).isoformat(),
              'protocol_sha256': sha256(protocol_path.read_bytes()).hexdigest(),
              'cases': len(records), 'disagreements': 0, 'verified_witnesses': witnesses,
              'forced_cover_exhaustive_checks': cover_checks, 'growth_exhaustive_checks': growth_checks,
              'resource_limit_returns_unknown': True, 'elapsed_seconds': time.perf_counter() - started,
              'records': records, 'scope': 'Exact synthetic finite controls, not biological validation or proof.'}
    result_path.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'records'}, indent=2))


if __name__ == '__main__':
    main()
