#!/usr/bin/env python3
"""Read-only independent replay of frozen inputs, returned points and summaries.

Imports no inference, benchmark, geometry or oracle implementation. Distances
are rebuilt by rooting the two components at each interior source point.
Does not rerun solvers or certify minimality of a returned feasible witness.
"""
from collections import Counter
from datetime import datetime, timezone
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
from statistics import median
import json

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return sha256(path.read_bytes()).hexdigest()


def adjacency(case):
    graph = {}
    for index, (u, v, length) in enumerate(case['edges']):
        length = F(length)
        assert length > 0 and u != v
        graph.setdefault(u, []).append((v, length, index))
        graph.setdefault(v, []).append((u, length, index))
    assert len(case['edges']) == len(graph) - 1
    assert all(len(n) in (1, 3) for n in graph.values())
    return graph


def check_witness(case, witness, cap):
    graph = adjacency(case)
    sources = witness['sources']
    assert 1 <= len(sources) <= cap
    assert len({s['edge'] for s in sources}) == len(sources)
    assert sum(F(s['p']) for s in sources) == 1
    y = F(witness['pendant'])
    assert y > 0
    if case['task'] == 'fixed_one_missing':
        assert y == F(case['pendant'])
    prediction = {int(r): y for r in case['observed']}
    for source in sources:
        p, x, e = F(source['p']), F(source['x']), source['edge']
        assert p > 0 and 0 < x < 1 and 0 <= e < len(case['edges'])
        u, v, length = case['edges'][e]
        length = F(length)
        distance = {u: x * length, v: (1 - x) * length}
        todo = [u, v]
        for a in todo:
            for b, ell, index in graph[a]:
                if index != e and b not in distance:
                    distance[b] = distance[a] + ell
                    todo.append(b)
        assert len(distance) == len(graph)
        for r in prediction:
            prediction[r] += p * distance[r]
    assert prediction == {int(r): F(d) for r, d in case['observed'].items()}


def fingerprint(case):
    return json.dumps([sorted((min(u, v), max(u, v), str(F(length)))
                             for u, v, length in case['edges']),
                       sorted((int(r), str(F(d))) for r, d in case['observed'].items()),
                       case['cap'], str(F(case['pendant'])) if 'pendant' in case else None])


def main():
    protocol_path = ROOT / 'results/stage44_confirmation_protocol.json'
    result_path = ROOT / 'results/stage44_confirmation.json'
    protocol, result = [json.loads(p.read_text()) for p in (protocol_path, result_path)]
    assert result['protocol_sha256'] == sha(protocol_path)
    assert result['summary']['status'] == 'PASS'
    for name, expected in protocol['script_hashes'].items():
        assert sha(ROOT / 'scripts' / name) == expected
    for name, expected in protocol['prior_protocol_hashes'].items():
        assert sha(ROOT / 'results' / name) == expected
    old = json.loads((ROOT / 'results/stage42_benchmark_protocol.json').read_text())['case_list']
    previous = {fingerprint(c) for c in old}
    fingerprints = set()
    cases = {c['name']: c for c in protocol['cases']}
    assert len(cases) == 48 and protocol['repetitions'] == 2
    for case in cases.values():
        graph = adjacency(case)
        leaves = {v for v, n in graph.items() if len(n) == 1}
        assert set(map(int, case['observed'])) == leaves - set(case['hidden_leaves'])
        if case['task'] == 'fixed_one_missing':
            assert len(case['hidden_leaves']) == 1
        check_witness(case, {'sources': case['generating_sources'],
                            'pendant': case['generating_pendant']}, case['generating_count'])
        key = fingerprint(case)
        assert key not in previous and key not in fingerprints
        fingerprints.add(key)
    witnesses, jobs = 0, 0
    stream_index = {(r['name'], r['repetition'], r['method']): r['result'] for r in
                    map(json.loads, (ROOT / 'results/stage44_confirmation_runs.jsonl').read_text().splitlines())}
    assert len(stream_index) == 256
    for record in result['records']:
        case = cases[record['name']]
        assert len(record['repetitions']) == 2
        answers = set()
        for rep in record['repetitions']:
            assert set(rep['runs']) == set(protocol['methods'][case['task']])
            for method, run in rep['runs'].items():
                jobs += 1
                assert stream_index[record['name'], rep['repetition'], method] == run
                assert run['status'] in ('complete', 'unknown')
                if run['status'] == 'unknown':
                    assert run['minimum'] is None
                    continue
                answers.add((run['feasible'], run['minimum']))
                if run['feasible']:
                    assert run['minimum'] == len(run['witness']['sources'])
                    check_witness(case, run['witness'], case['cap'])
                    witnesses += 1
                else:
                    assert run['minimum'] is None and 'witness' not in run
        assert len(answers) <= 1
    assert jobs == 256
    for task, methods in protocol['methods'].items():
        records = [r for r in result['records'] if r['task'] == task]
        stored = result['summary']['tasks'][task]
        for method in methods:
            runs = [rep['runs'][method] for r in records for rep in r['repetitions']]
            assert dict(Counter(r['status'] for r in runs)) == stored['methods'][method]['status_counts']
        for a, b in combinations(methods, 2):
            joint = [r for r in records if all(rep['runs'][m]['status'] == 'complete'
                     for rep in r['repetitions'] for m in (a, b))]
            data = {m: [median(rep['runs'][m]['child_cpu_seconds'] for rep in r['repetitions'])
                        for r in joint] for m in (a, b)}
            expected = stored['pairwise'][f'{a} vs {b}']
            assert len(joint) == expected['joint_complete_both_repetitions']
            for m in (a, b):
                assert (median(data[m]) if joint else None) == expected['median_case_cpu_seconds'][m]
            assert (median(x / y for x, y in zip(data[a], data[b])) if joint else None) == \
                expected['median_case_ratio_first_over_second']
    posthoc_protocol_path = ROOT / 'results/stage44_confirmation_posthoc_protocol.json'
    posthoc_path = ROOT / 'results/stage44_confirmation_posthoc.json'
    posthoc_protocol = json.loads(posthoc_protocol_path.read_text())
    posthoc = json.loads(posthoc_path.read_text())
    assert posthoc['protocol_sha256'] == sha(posthoc_protocol_path)
    assert posthoc_protocol['primary_results_sha256'] == sha(result_path)
    assert posthoc_protocol['primary_protocol_sha256'] == sha(protocol_path)
    for name, expected in posthoc_protocol['script_hashes'].items():
        assert sha(ROOT / 'scripts' / name) == expected
    posthoc_witnesses = 0
    for run in posthoc['runs'].values():
        assert run['status'] == 'complete' and run['minimum'] == 5
        check_witness(cases[posthoc['name']], run['witness'], 5)
        posthoc_witnesses += 1
    report = {'status': 'PASS', 'created_utc': datetime.now(timezone.utc).isoformat(),
              'protocol_sha256': sha(protocol_path), 'result_sha256': sha(result_path),
              'audit_script_sha256': sha(Path(__file__)), 'new_input_cases': len(cases),
              'generating_witnesses_checked': len(cases), 'returned_witnesses_replayed': witnesses,
              'jobs_checked': jobs, 'stream_records_agree': True,
              'separate_posthoc_witnesses_checked': posthoc_witnesses,
              'primary_results_unchanged_after_posthoc': True,
              'all_completed_answers_consistent': True, 'paired_medians_recomputed': True,
              'frozen_source_hashes_checked': len(protocol['script_hashes']),
              'scope': 'Independent read-only witness/summary replay; no solver rerun or independent proof of infeasible/minimal outcomes.'}
    (ROOT / 'results/stage44_confirmation_replay.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
