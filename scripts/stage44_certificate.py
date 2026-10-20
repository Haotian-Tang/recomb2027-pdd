#!/usr/bin/env python3
"""Exact controls for linear-time, pendant-assisted path-excess screening.

Tree-metric identity and message passing are classical. This script checks a
conditional PDD corollary and compares it with the frozen Stage43 screen.
It imports the frozen geometry generator solely to use exactly the same data.
The new screening routine uses no distance matrix and no optimizer.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json
import random
from stage43_theory_robust_mass import cherry_tree, forward, bounds, tree

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/stage44_certificate.json'
PROTOCOL = ROOT / 'results/stage44_certificate_protocol.json'
PROTOCOL_HASH = '76f78ffb148785e65f17fa1b7a6958609091cd74d63ab789d1739b7a2392ea46'


def smallest(items, count=3):
    """Constant-size insertion; linear in input size, including high degrees."""
    answer = []
    for item in items:
        index = 0
        while index < len(answer) and answer[index][0] <= item[0]:
            index += 1
        answer.insert(index, item)
        if len(answer) > count:
            answer.pop()
    return answer


def path_screen(edges, upper, pendant_floor, alpha, rho):
    """Return safe candidate edges and endpoint-local best path certificates.

    Every original leaf must be in upper. Positive tree lengths; no degree-two
    vertices. Every allowed source has weight >= alpha and fractional distance
    >= rho from both endpoints. The allowed pendant is >= pendant_floor.
    The input dictionary is an upper bound on the unknown true leaf vector.
    """
    adjacency = {}
    for u, v, ell in edges:
        assert ell > 0
        adjacency.setdefault(u, {})[v] = ell
        adjacency.setdefault(v, {})[u] = ell
    leaves = {v for v in adjacency if len(adjacency[v]) == 1}
    assert set(upper) == leaves and len(edges) + 1 == len(adjacency)
    assert all(len(row) == 1 or len(row) >= 3 for row in adjacency.values())
    assert pendant_floor > 0 and alpha > 0 and 0 < rho <= F(1, 2)
    root = next((v for v in adjacency if v not in leaves), next(iter(adjacency)))
    parent, order = {root: None}, [root]
    for v in order:
        for u in adjacency[v]:
            if u != parent[v]:
                assert u not in parent
                parent[u] = v
                order.append(u)
    assert len(order) == len(adjacency)
    # message[v,u] = (min upper(r)-dist(v,r) toward u, attaining leaf r).
    message = {}
    for v in reversed(order):
        options = [(upper[v], v)] if v in leaves else []
        options += [message[v, u] for u in adjacency[v] if parent.get(u) == v]
        best = min(options, key=lambda x: x[0])
        if parent[v] is not None:
            p = parent[v]
            message[p, v] = (best[0] - adjacency[v][p], best[1])
    for v in order:
        # Include self leaf as a separate option, needed for the two-leaf tree.
        options = [(upper[v], None, v)] if v in leaves else []
        options += [(message[v, u][0], u, message[v, u][1])
                    for u in adjacency[v] if (v, u) in message]
        first_two = smallest(options, 2)
        for u in adjacency[v]:
            if parent.get(u) == v:
                best = next(item for item in first_two if item[1] != u)
                message[u, v] = (best[0] - adjacency[v][u], best[2])
    first_three = {v: smallest([(message[v, u][0], u, message[v, u][1])
                               for u in adjacency[v]], 3) for v in adjacency}
    records = []
    retained = []
    for j, (u, v, ell) in enumerate(edges):
        options = []
        for endpoint, away in ((u, v), (v, u)):
            choices = [item for item in first_three[endpoint] if item[1] != away]
            if len(choices) >= 2:
                excess = (choices[0][0] + choices[1][0]) / 2 - pendant_floor
                options.append({'upper_excess': excess, 'endpoint': endpoint,
                                'leaf_pair': [choices[0][2], choices[1][2]]})
        witness = min(options, key=lambda x: x['upper_excess']) if options else None
        threshold = alpha * rho * ell
        if witness is None or witness['upper_excess'] >= threshold:
            retained.append(j)
        records.append({'edge_index': j, 'threshold': threshold, 'witness': witness})
    return retained, records, message


def brute_pair_check(geometry, upper, floor, records, message):
    edges, adjacency, leaves, distance = geometry
    for (v, u), (score, witness_leaf) in message.items():
        branch = [r for r in leaves if distance[v][r] == adjacency[v][u] + distance[u][r]]
        assert score == min(upper[r] - distance[v][r] for r in branch)
        assert witness_leaf in branch and score == upper[witness_leaf] - distance[v][witness_leaf]
    for edge, record in zip(edges, records):
        u, v, _ = edge
        options = []
        for endpoint, away in ((u, v), (v, u)):
            for i, r in enumerate(leaves):
                for s in leaves[i+1:]:
                    # Path r--s contains endpoint and avoids the edge interior.
                    through_endpoint = distance[r][s] == distance[r][endpoint] + distance[endpoint][s]
                    r_wrong_side = distance[endpoint][r] == adjacency[endpoint][away] + distance[away][r]
                    s_wrong_side = distance[endpoint][s] == adjacency[endpoint][away] + distance[away][s]
                    if through_endpoint and not r_wrong_side and not s_wrong_side:
                        options.append((upper[r] + upper[s] - distance[r][s]) / 2 - floor)
        if options:
            assert record['witness']['upper_excess'] == min(options)
        else:
            assert record['witness'] is None


def old_screen(geometry, nominal, epsilon, floor, alpha, rho):
    _, _, _, gap, intervals = bounds(geometry, nominal, epsilon, floor)
    excluded = {v for v, (_, hi) in intervals.items() if hi < alpha * rho}
    kept = [j for j, (u, v, _) in enumerate(geometry[0]) if u not in excluded and v not in excluded]
    return kept, gap


def frozen_cases():
    rng = random.Random(430317)
    archived = json.loads((ROOT / 'results/stage43_theory_robust_mass.json').read_text())['mixed_cases']
    index = 0
    for number in (2, 3, 4, 6, 10):
        for unequal in (False, True):
            geometry = cherry_tree(number, unequal)
            edges, adjacency, leaves, _ = geometry
            for trial in range(20):
                count = 1 + trial % min(8, len(edges))
                chosen = rng.sample(range(len(edges)), count)
                raw = [rng.randrange(1, 12) for _ in chosen]
                total = sum(raw)
                sources = [(e, F(rng.randrange(1, 10), 10), F(w, total)) for e, w in zip(chosen, raw)]
                pendant = F(1 + trial % 4, 7)
                masses, potential = forward(geometry, sources, pendant)
                epsilon = F(1, (1000, 10**6, 10**12)[trial % 3])
                nominal = {r: potential[r] + epsilon * F(rng.randrange(-4, 5), 4) for r in leaves}
                floor = pendant if trial % 2 == 0 else pendant / 2
                alpha = min(p for _, _, p in sources)
                rho = min(min(x, 1-x) for _, x, _ in sources)
                kept, gap = old_screen(geometry, nominal, epsilon, floor, alpha, rho)
                old = archived[index]
                assert old['cherries'] == number and old['unequal_lengths'] == unequal and old['trial'] == trial
                assert old['source_count'] == count and old['epsilon'] == str(epsilon)
                assert old['gap_upper'] == str(gap) and old['candidate_edges'] == len(kept)
                assert old['pendant_floor_equals_truth'] == (floor == pendant)
                yield index, geometry, sources, masses, potential, epsilon, nominal, floor, alpha, rho, kept
                index += 1
    assert index == len(archived)


def compare_record(old, new, edges, truth):
    a, b = set(old), set(new)
    assert truth <= a and truth <= b
    relation = ('equal' if a == b else 'path_strictly_tighter' if b < a else
                'mass_strictly_tighter' if a < b else 'incomparable')
    return {'original_edges': len(edges), 'true_edges': sorted(truth), 'mass_retained': sorted(a),
            'path_retained': sorted(b), 'intersection_retained': sorted(a & b), 'relation': relation}


def summarize(records):
    counts = {}
    for row in records:
        counts[row['relation']] = counts.get(row['relation'], 0) + 1
    return {'cases': len(records), 'relations': counts,
            'total_original_edges': sum(row['original_edges'] for row in records),
            'total_mass_retained': sum(len(row['mass_retained']) for row in records),
            'total_path_retained': sum(len(row['path_retained']) for row in records),
            'total_intersection_retained': sum(len(row['intersection_retained']) for row in records),
            'all_true_edges_retained': True}


def main(write=True):
    assert hashlib.sha256(PROTOCOL.read_bytes()).hexdigest() == PROTOCOL_HASH
    protocol = json.loads(PROTOCOL.read_text())
    for name, expected in protocol['frozen_inputs'].items():
        assert hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == expected
    mixed, hetero = [], []
    for index, geometry, sources, masses, potential, epsilon, nominal, floor, alpha, rho, old in frozen_cases():
        upper = {r: nominal[r] + epsilon for r in geometry[2]}
        new, records, messages = path_screen(geometry[0], upper, floor, alpha, rho)
        brute_pair_check(geometry, upper, floor, records, messages)
        assert all(row['witness'] is None or row['witness']['upper_excess'] >= 0 for row in records)
        item = compare_record(old, new, geometry[0], {e for e, _, _ in sources})
        item.update({'case_index': index, 'alpha': alpha, 'rho': rho, 'pendant_floor': floor,
                     'epsilon': epsilon})
        mixed.append(item)
        upper_h = {r: potential[r] + epsilon * F(1 + i % 4, 4) for i, r in enumerate(geometry[2])}
        lower_h = {r: potential[r] - epsilon * F(1 + i % 5, 5) for i, r in enumerate(geometry[2])}
        assert all(lower_h[r] <= potential[r] <= upper_h[r] for r in geometry[2])
        kept_h, rec_h, msg_h = path_screen(geometry[0], upper_h, floor, alpha, rho)
        brute_pair_check(geometry, upper_h, floor, rec_h, msg_h)
        assert {e for e, _, _ in sources} <= set(kept_h)
        hetero.append({'case_index': index, 'retained': kept_h, 'original_edges': len(geometry[0])})
    one = []
    for number in (2, 3, 4, 6, 10):
        for unequal in (False, True):
            geometry = cherry_tree(number, unequal)
            for chosen in (0, len(geometry[0]) - 1):
                _, potential = forward(geometry, [(chosen, F(1, 2), F(1))], F(1))
                alpha, rho = F(1, 4), F(1, 4)
                threshold = alpha * rho * min(ell for _, _, ell in geometry[0])
                for ratio in (F(0), F(1, 4), F(1, 2), F(3, 4), F(1), F(5, 4), F(2)):
                    epsilon = ratio * threshold
                    nominal = {r: potential[r] for r in geometry[2]}
                    old, _ = old_screen(geometry, nominal, epsilon, F(1), alpha, rho)
                    new, records, messages = path_screen(geometry[0], {r: nominal[r] + epsilon for r in geometry[2]},
                                                         F(1), alpha, rho)
                    brute_pair_check(geometry, {r: nominal[r] + epsilon for r in geometry[2]}, F(1), records, messages)
                    if ratio < 1:
                        assert new == [chosen]
                    item = compare_record(old, new, geometry[0], {chosen})
                    item.update({'cherries': number, 'unequal_lengths': unequal, 'true_edge': chosen,
                                 'error_over_sharp_uniform_threshold': ratio, 'epsilon': epsilon})
                    one.append(item)
    star = tree([('c', f'a{i}', F(1)) for i in range(3)])
    _, target = forward(star, [(0, F(3, 8), F(1))], F(1))
    _, alternative = forward(star, [(0, F(1, 2), F(3, 4)), (1, F(1, 4), F(1, 4))], F(1))
    delta = {r: alternative[r] - target[r] for r in star[2]}
    epsilon = max(abs(value) for value in delta.values())
    assert epsilon == F(1, 16)
    at, records, _ = path_screen(star[0], {r: target[r] + epsilon for r in star[2]}, F(1), F(1, 4), F(1, 4))
    below, _, _ = path_screen(star[0], {r: target[r] + epsilon * F(999, 1000) for r in star[2]}, F(1), F(1, 4), F(1, 4))
    assert {0, 1} <= set(at) and below == [0]
    # Degenerate two-leaf tree: there is no eligible avoiding leaf pair.
    two = tree([('a', 'b', F(2))])
    _, pot_two = forward(two, [(0, F(1, 3), F(1))], F(1))
    kept_two, rec_two, msg_two = path_screen(two[0], pot_two, F(1), F(1, 4), F(1, 4))
    brute_pair_check(two, pot_two, F(1), rec_two, msg_two)
    assert kept_two == [0] and rec_two[0]['witness'] is None
    # An infeasible pendant floor is detectable without source-margin assumptions.
    _, bad_rec, _ = path_screen(star[0], {r: target[r] for r in star[2]}, F(2), F(1, 4), F(1, 4))
    assert any(row['witness'] and row['witness']['upper_excess'] < 0 for row in bad_rec)
    result = {'stage': 44, 'status': 'all_exact_assertions_passed', 'protocol_sha256': PROTOCOL_HASH,
              'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'claim_boundary': protocol['claim_boundary'], 'mixed_summary': summarize(mixed),
              'one_source_summary': summarize(one), 'heteroscedastic_cases': len(hetero),
              'heteroscedastic_total_retained': sum(len(x['retained']) for x in hetero),
              'strict_one_source_cases': sum(x['error_over_sharp_uniform_threshold'] < 1 for x in one),
              'sharpness': {'alpha': F(1, 4), 'rho': F(1, 4), 'minimum_edge_length': F(1),
                            'exact_sup_difference': epsilon, 'leaf_differences': delta,
                            'retained_at_equality': at, 'retained_just_below': below},
              'two_leaf_control': 'passed', 'negative_excess_infeasibility_control': 'passed',
              'mixed_cases': mixed, 'heteroscedastic_cases_details': hetero, 'one_source_cases': one}
    if write:
        assert not OUT.exists(), 'Refusing to overwrite results; use --check'
        OUT.write_text(json.dumps(result, indent=2, default=str) + '\n')
    print(json.dumps({key: result[key] for key in ['status', 'mixed_summary', 'one_source_summary',
                                                'heteroscedastic_cases', 'strict_one_source_cases', 'sharpness']},
                     indent=2, default=str))
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    main(write=not args.check)
