#!/usr/bin/env python3
"""Frozen, bounded correctness follow-up; no solver changes or timing claims.

Select every Stage44 case lacking completed answers from both the original
support oracle and the current structural method, after its separate posthoc
check. Run only the missing method, once, with a 30-second wall cap.
"""
import argparse
from datetime import datetime, timezone
from fractions import Fraction as F
import json
from pathlib import Path
import time

from stage42_benchmark import sha
from stage44_confirmation import isolated_run
from stage44_confirmation_replay import check_witness

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT / 'results/stage45_remaining_protocol.json'
RESULT = ROOT / 'results/stage45_remaining_checks.json'
STREAM = ROOT / 'results/stage45_remaining_runs.jsonl'
REPORT = ROOT / 'review/stage45_remaining_checks.md'
EVIDENCE = [
    'results/stage44_confirmation_protocol.json',
    'results/stage44_confirmation.json',
    'results/stage44_confirmation_runs.jsonl',
    'results/stage44_confirmation_posthoc_protocol.json',
    'results/stage44_confirmation_posthoc.json',
]
REQUIRED = {
    'free_pendant': ['original_support', 'stage43_equality_reduced'],
    'fixed_one_missing': ['original_fixed_oracle', 'stage43_fixed_one'],
}


def read(relative):
    return json.loads((ROOT / relative).read_text())


def selected_jobs():
    primary = read('results/stage44_confirmation.json')
    prior_posthoc = read('results/stage44_confirmation_posthoc.json')
    cases = {c['name']: c for c in read('results/stage44_confirmation_protocol.json')['cases']}
    jobs, already_paired = [], 0
    for record in primary['records']:
        completed = {m for rep in record['repetitions'] for m, r in rep['runs'].items()
                     if r['status'] == 'complete'}
        runs = [r for rep in record['repetitions'] for r in rep['runs'].values()
                if r['status'] == 'complete']
        if prior_posthoc['name'] == record['name']:
            completed |= {m for m, r in prior_posthoc['runs'].items() if r['status'] == 'complete'}
            runs += [r for r in prior_posthoc['runs'].values() if r['status'] == 'complete']
        answers = {(r['feasible'], r['minimum']) for r in runs}
        assert len(answers) == 1, 'Every retained case must have one consistent existing answer'
        missing = [m for m in REQUIRED[record['task']] if m not in completed]
        if not missing:
            already_paired += 1
        for method in missing:
            feasible, minimum = next(iter(answers))
            jobs.append(dict(name=record['name'], method=method, case=cases[record['name']],
                             previous_completed_methods=sorted(completed),
                             expected_answer=dict(feasible=feasible, minimum=minimum)))
    return jobs, already_paired


def prepare():
    assert not any(p.exists() for p in (PROTOCOL, RESULT, STREAM)), 'Never overwrite evidence'
    parent = read('results/stage44_confirmation_protocol.json')
    posthoc = read('results/stage44_confirmation_posthoc_protocol.json')
    hashes = {**parent['script_hashes'], **posthoc['script_hashes']}
    for name in (Path(__file__).name, 'stage44_confirmation_replay.py'):
        hashes[name] = sha(ROOT / 'scripts' / name)
    for name, expected in hashes.items():
        assert sha(ROOT / 'scripts' / name) == expected
    jobs, already_paired = selected_jobs()
    assert len(jobs) == 5 and already_paired == 43
    assert len({j['name'] for j in jobs}) == len(jobs)
    protocol = dict(
        stage=45, registered_utc=datetime.now(timezone.utc).isoformat(),
        purpose='Close remaining exact-answer agreement gaps; no performance or new-input evidence',
        selection='All and only Stage44 cases missing a completed original-support oracle or current structural answer after the Stage44 posthoc result. Selection is outcome-dependent and explicitly posthoc.',
        required_methods=REQUIRED, jobs=jobs, already_paired_cases=already_paired,
        per_job_hard_wall_seconds=30.0, per_job_address_space_bytes=1024**3,
        repetitions=1, maximum_solver_wall_budget_seconds=150.0,
        escalation='None. One attempt per missing method; UNKNOWN retained at 30 seconds.',
        execution_order='Primary registered case order, then required-method order; fresh sequential Linux fork per job',
        thread_environment={key: '1' for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS')},
        analysis='Compare completed exact feasible/minimum answers to every prior completed answer. Replay every returned feasible witness using the independent component-distance checker. A witness alone does not verify minimality or NO. Never replace primary timeouts, repetitions, timing summaries, or outputs; no aggregate speed claims.',
        source_hashes=hashes,
        previous_evidence_sha256={name: sha(ROOT / name) for name in EVIDENCE},
    )
    PROTOCOL.write_text(json.dumps(protocol, indent=2) + '\n')
    print(json.dumps(dict(protocol=str(PROTOCOL), sha256=sha(PROTOCOL),
                          jobs=[(j['name'], j['method']) for j in jobs]), indent=2))


def run():
    assert not RESULT.exists() and not STREAM.exists(), 'No reruns or evidence replacement'
    protocol = json.loads(PROTOCOL.read_text())
    for name, expected in protocol['source_hashes'].items():
        assert sha(ROOT / 'scripts' / name) == expected
    for name, expected in protocol['previous_evidence_sha256'].items():
        assert sha(ROOT / name) == expected
    jobs, already_paired = selected_jobs()
    assert jobs == protocol['jobs'] and already_paired == protocol['already_paired_cases']
    started, records = time.perf_counter(), []
    for job in jobs:
        case = dict(job['case'])
        case['edges'] = [(u, v, F(ell)) for u, v, ell in case['edges']]
        case['observed'] = {int(r): F(d) for r, d in case['observed'].items()}
        if 'pendant' in case:
            case['pendant'] = F(case['pendant'])
        answer = isolated_run(job['method'], case, protocol['per_job_hard_wall_seconds'],
                              protocol['per_job_address_space_bytes'])
        agrees = None
        replayed = False
        if answer['status'] == 'complete':
            agrees = {k: answer[k] for k in ('feasible', 'minimum')} == job['expected_answer']
            if answer['feasible']:
                check_witness(job['case'], answer['witness'], job['case']['cap'])
                assert len(answer['witness']['sources']) == answer['minimum']
                replayed = True
        record = dict(name=job['name'], method=job['method'], previous_answer=job['expected_answer'],
                      run=answer, exact_answer_agrees=agrees,
                      independent_component_distance_replay=replayed)
        records.append(record)
        with STREAM.open('a') as handle:
            handle.write(json.dumps(record) + '\n')
        print(json.dumps({k: v for k, v in record.items() if k != 'run'} |
                         dict(status=answer['status'], minimum=answer.get('minimum'))), flush=True)
    old_unchanged = all(sha(ROOT / name) == expected
                        for name, expected in protocol['previous_evidence_sha256'].items())
    assert old_unchanged
    completed = sum(r['run']['status'] == 'complete' for r in records)
    consistent = all(r['exact_answer_agrees'] is not False for r in records)
    no_errors = all(r['run']['status'] in ('complete', 'unknown') for r in records)
    summary = dict(status='PASS' if consistent and no_errors else 'FAIL',
                   completed_jobs=completed,
                   unknown_jobs=sum(r['run']['status'] == 'unknown' for r in records),
                   new_returned_witnesses_replayed=sum(r['independent_component_distance_replay'] for r in records),
                   all_completed_exact_answers_agree=consistent,
                   paired_original_oracle_and_current_structural_cases=already_paired + completed if consistent else None,
                   total_registered_cases=48, previous_evidence_unchanged=old_unchanged,
                   remaining_cases=[r['name'] for r in records if r['run']['status'] != 'complete'])
    result = dict(stage=45, created_utc=datetime.now(timezone.utc).isoformat(),
                  protocol_sha256=sha(PROTOCOL), records=records, summary=summary,
                  diagnostic_elapsed_seconds=time.perf_counter() - started,
                  scope='Posthoc exact-answer agreement only; shared exact arithmetic utilities limit full implementation independence. No formal proof, biological validation, novelty clearance, or performance evidence.')
    RESULT.write_text(json.dumps(result, indent=2) + '\n')
    lines = [
        '# Stage45: bounded completion of remaining exact-answer checks', '',
        'This is a correctness-only posthoc follow-up. Every Stage44 input, solver, result, timeout and timing table is unchanged. No new random inputs and no solver changes were introduced.', '',
        '## Selection and fixed budget', '',
        'The selection rule is exhaustive: after the separate Stage44 posthoc check, choose every input for which either the original support oracle or the current structural method has not completed. This selected five inputs; 43/48 already had both methods complete. The older structural implementation alone is not accepted as an independent original-support oracle.', '',
        'The protocol and all imported local source files were hashed before running these five jobs. Each missing method received one fresh single-threaded job, 30 seconds hard wall time and 1 GiB address space. No escalation or retries were allowed. The maximum solver wall allocation was 150 seconds. Returned witnesses were separately replayed with exact rational component distances.', '',
        '| Case | Missing method run | Existing answer | Follow-up status | Agreement |',
        '|---|---|---|---|---|',
    ]
    for row in records:
        prior = row['previous_answer']
        answer_text = f"minimum {prior['minimum']}" if prior['feasible'] else 'NO within registered cap'
        agreement = 'yes' if row['exact_answer_agrees'] else ('UNKNOWN' if row['exact_answer_agrees'] is None else 'DISAGREEMENT')
        lines.append(f"| `{row['name']}` | `{row['method']}` | {answer_text} | {row['run']['status']} | {agreement} |")
    lines += ['', '## Result and limits', '',
        f"Completed {completed}/5 follow-up jobs; {summary['unknown_jobs']} remain UNKNOWN. All completed exact answers agree: {consistent}. Newly returned feasible witnesses replayed: {summary['new_returned_witnesses_replayed']}. Both the original oracle and current structural method have now completed on {summary['paired_original_oracle_and_current_structural_cases']}/48 registered cases.", '',
        'Agreement on an exact minimum or NO is stronger than replaying a feasible witness alone: witness replay checks feasibility but not absence of a smaller solution. These comparisons use different original-support and structural enumeration routes, but share some exact arithmetic helpers; they are not formal verification or fully independent software implementations. No timing or performance conclusion is drawn from the longer selected runs.', '',
        f"Protocol SHA256: `{sha(PROTOCOL)}`.",
        f"Result SHA256: `{sha(RESULT)}`.",
        f"Previous evidence hashes remain unchanged: {old_unchanged}.", '',
        'Files: `results/stage45_remaining_protocol.json`, `results/stage45_remaining_checks.json`, `results/stage45_remaining_runs.jsonl`, `scripts/stage45_remaining_checks.py`.', '',
    ]
    REPORT.write_text('\n'.join(lines))
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('phase', choices=('prepare', 'run'))
    args = parser.parse_args()
    prepare() if args.phase == 'prepare' else run()
