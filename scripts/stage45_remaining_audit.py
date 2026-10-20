#!/usr/bin/env python3
"""Read-only evidence audit; no solver imports, writes, or timing repetitions."""
import argparse
from collections import Counter
from itertools import combinations
import json
from pathlib import Path
from statistics import median

from stage44_confirmation_replay import check_witness, fingerprint, sha

ROOT = Path(__file__).resolve().parents[1]


def read(relative):
    return json.loads((ROOT / relative).read_text())


def audit():
    manifest_path = ROOT / 'MANIFEST.json'
    manifest_count = None
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
        for name, expected in manifest.items():
            assert sha(ROOT / name) == expected, name
        manifest_count = len(manifest)
    p = read('results/stage44_confirmation_protocol.json')
    primary = read('results/stage44_confirmation.json')
    postp = read('results/stage44_confirmation_posthoc_protocol.json')
    post = read('results/stage44_confirmation_posthoc.json')
    newp = read('results/stage45_remaining_protocol.json')
    new = read('results/stage45_remaining_checks.json')
    checked_source_names = set()
    for protocol in (p, postp, newp):
        for name, expected in protocol.get('script_hashes', protocol.get('source_hashes', {})).items():
            assert sha(ROOT / 'scripts' / name) == expected, name
            checked_source_names.add(name)
    for name, expected in p['prior_protocol_hashes'].items():
        assert sha(ROOT / 'results' / name) == expected
    for name, expected in newp['previous_evidence_sha256'].items():
        assert sha(ROOT / name) == expected
    assert primary['protocol_sha256'] == sha(ROOT / 'results/stage44_confirmation_protocol.json')
    assert post['protocol_sha256'] == sha(ROOT / 'results/stage44_confirmation_posthoc_protocol.json')
    assert postp['primary_results_sha256'] == sha(ROOT / 'results/stage44_confirmation.json')
    assert postp['primary_protocol_sha256'] == primary['protocol_sha256']
    assert new['protocol_sha256'] == sha(ROOT / 'results/stage45_remaining_protocol.json')
    cases = {case['name']: case for case in p['cases']}
    assert len(cases) == 48
    old = read('results/stage42_benchmark_protocol.json')['case_list']
    old_keys = {fingerprint(c) for c in old}
    new_keys = set()
    for case in cases.values():
        check_witness(case, {'sources': case['generating_sources'],
                            'pendant': case['generating_pendant']}, case['generating_count'])
        key = fingerprint(case)
        assert key not in old_keys and key not in new_keys
        new_keys.add(key)
    stream = {(r['name'], r['repetition'], r['method']): r['result'] for r in
              map(json.loads, (ROOT / 'results/stage44_confirmation_runs.jsonl').read_text().splitlines())}
    assert len(stream) == 256
    completed = {name: [] for name in cases}
    primary_witnesses = 0
    for record in primary['records']:
        for rep in record['repetitions']:
            for method, run in rep['runs'].items():
                assert stream[record['name'], rep['repetition'], method] == run
                assert run['status'] in ('complete', 'unknown')
                if run['status'] == 'complete':
                    completed[record['name']].append((method, run))
                    primary_witnesses += int(run['feasible'])
                else:
                    assert run['minimum'] is None
    assert primary_witnesses == 192
    for task, methods in p['methods'].items():
        records = [r for r in primary['records'] if r['task'] == task]
        stored = primary['summary']['tasks'][task]
        for method in methods:
            runs = [rep['runs'][method] for r in records for rep in r['repetitions']]
            assert dict(Counter(r['status'] for r in runs)) == stored['methods'][method]['status_counts']
        for a, b in combinations(methods, 2):
            joint = [r for r in records if all(rep['runs'][m]['status'] == 'complete'
                     for rep in r['repetitions'] for m in (a, b))]
            values = {m: [median(rep['runs'][m]['child_cpu_seconds'] for rep in r['repetitions'])
                          for r in joint] for m in (a, b)}
            expected = stored['pairwise'][f'{a} vs {b}']
            assert len(joint) == expected['joint_complete_both_repetitions']
            for m in (a, b):
                assert (median(values[m]) if joint else None) == expected['median_case_cpu_seconds'][m]
            assert (median(x / y for x, y in zip(values[a], values[b])) if joint else None) == \
                expected['median_case_ratio_first_over_second']
    assert all(run['status'] == 'complete' for run in post['runs'].values())
    completed[post['name']] += list(post['runs'].items())
    selected = []
    for name, runs in completed.items():
        methods = {m for m, run in runs}
        selected += [(name, method) for method in newp['required_methods'][cases[name]['task']]
                     if method not in methods]
    assert selected == [(job['name'], job['method']) for job in newp['jobs']]
    assert len(selected) == 5
    assert new['records'] == list(map(json.loads, (ROOT / 'results/stage45_remaining_runs.jsonl').read_text().splitlines()))
    assert len(new['records']) == 5
    for record in new['records']:
        assert record['run']['status'] == 'complete'
        assert record['exact_answer_agrees']
        completed[record['name']].append((record['method'], record['run']))
    witnesses, counts = 0, Counter()
    for name, runs in completed.items():
        assert set(newp['required_methods'][cases[name]['task']]) <= {m for m, run in runs}
        answers = {(run['feasible'], run['minimum']) for m, run in runs}
        assert len(answers) == 1
        feasible, minimum = next(iter(answers))
        counts['feasible' if feasible else 'infeasible'] += 1
        for method, run in runs:
            if run['feasible']:
                check_witness(cases[name], run['witness'], cases[name]['cap'])
                assert len(run['witness']['sources']) == run['minimum']
                witnesses += 1
            else:
                assert run['minimum'] is None and 'witness' not in run
    assert witnesses == 197 and counts == {'feasible': 40, 'infeasible': 8}
    certp = read('results/stage44_certificate_protocol.json')
    cert = read('results/stage44_certificate.json')
    assert cert['protocol_sha256'] == sha(ROOT / 'results/stage44_certificate_protocol.json')
    assert cert['script_sha256'] == sha(ROOT / 'scripts/stage44_certificate.py')
    for name, expected in certp['frozen_inputs'].items():
        assert sha(ROOT / name) == expected
    build = read('results/stage45_paper_build.json')
    for name, expected in build['source_sha256'].items():
        assert sha(ROOT / name) == expected
    for name, expected in build['pdf_sha256'].items():
        assert sha(ROOT / 'exports' / name) == expected
    return dict(status='PASS', manifest_files_checked=manifest_count,
                frozen_source_files_checked=len(checked_source_names),
                input_cases=48, generating_witnesses_replayed=48,
                primary_jobs_checked=256, primary_witnesses_replayed=192,
                separate_posthoc_witnesses_replayed=5, all_returned_witnesses_replayed=197,
                paired_original_and_structural_answer_agreements=48,
                answer_counts=dict(counts), missing_pair_selection_rederived=True,
                primary_paired_medians_recomputed=True,
                certificate_source_hashes_verified=True, manuscript_build_hashes_verified=True,
                files_written=0, solver_calls=0,
                scope='Exact witness/record/hash replay and cross-method answer agreement, not formal minimality proof, external review, biological validation or novelty clearance.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--audit', action='store_true', required=True,
                        help='Read and check frozen evidence without running solvers or writing files')
    parser.parse_args()
    print(json.dumps(audit(), indent=2))
