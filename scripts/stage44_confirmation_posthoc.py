#!/usr/bin/env python3
"""Separate posthoc correctness-only attempt on the sole all-UNKNOWN case.

Never modifies the primary 2-second benchmark or uses these timings in ratios.
"""
import argparse
from datetime import datetime, timezone
from fractions import Fraction as F
import json
from pathlib import Path
from stage42_benchmark import sha
from stage44_confirmation import isolated_run
ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT / 'results/stage44_confirmation_posthoc_protocol.json'
RESULT = ROOT / 'results/stage44_confirmation_posthoc.json'
PRIMARY_P = ROOT / 'results/stage44_confirmation_protocol.json'
PRIMARY_R = ROOT / 'results/stage44_confirmation.json'
NAME = 'fixed_random_15_n12_cap5_g6_hidden0'


def prepare():
    assert not PROTOCOL.exists() and not RESULT.exists()
    parent, primary = json.loads(PRIMARY_P.read_text()), json.loads(PRIMARY_R.read_text())
    records = [r for r in primary['records'] if all(v['status'] == 'unknown'
        for rep in r['repetitions'] for v in rep['runs'].values())]
    assert len(records) == 1 and records[0]['name'] == NAME
    chosen = next(c for c in parent['cases'] if c['name'] == NAME)
    hashes = dict(parent['script_hashes'])
    hashes[Path(__file__).name] = sha(__file__)
    for name, expected in hashes.items():
        assert sha(ROOT / 'scripts' / name) == expected
    protocol = dict(stage=44, registered_utc=datetime.now(timezone.utc).isoformat(),
        purpose='Posthoc correctness-only follow-up on the sole case unknown to every primary method; no tuning and no primary benchmark alteration',
        selection='Chosen after primary results by the all-methods-UNKNOWN criterion; this is not held-out performance evidence',
        case=chosen, methods=['original_fixed_oracle', 'stage43_fixed_one'], repetitions=1,
        per_job_hard_wall_seconds=30.0, per_job_address_space_bytes=1024**3,
        analysis='Preserve primary 2 s UNKNOWN statuses and all primary completion/timing tables. Compare exact answers if both complete. Validate any witness on original tree with prescribed pendant. These extended-cap timings cannot enter primary benchmark ratios.',
        primary_protocol_sha256=sha(PRIMARY_P), primary_results_sha256=sha(PRIMARY_R),
        script_hashes=hashes)
    PROTOCOL.write_text(json.dumps(protocol, indent=2) + '\n')
    print(json.dumps(dict(protocol=str(PROTOCOL), sha256=sha(PROTOCOL), case=NAME), indent=2))


def run():
    assert not RESULT.exists()
    protocol = json.loads(PROTOCOL.read_text())
    assert sha(PRIMARY_P) == protocol['primary_protocol_sha256']
    assert sha(PRIMARY_R) == protocol['primary_results_sha256']
    for name, expected in protocol['script_hashes'].items():
        assert sha(ROOT / 'scripts' / name) == expected
    case = dict(protocol['case'])
    case['edges'] = [(u, v, F(ell)) for u, v, ell in case['edges']]
    case['observed'] = {int(r): F(d) for r, d in case['observed'].items()}
    case['pendant'] = F(case['pendant'])
    runs = {}
    for method in protocol['methods']:
        runs[method] = isolated_run(method, case, protocol['per_job_hard_wall_seconds'],
                                    protocol['per_job_address_space_bytes'])
        print(method, runs[method]['status'], runs[method].get('minimum'), flush=True)
    answers = {(v['feasible'], v['minimum']) for v in runs.values() if v['status'] == 'complete'}
    validation_failures = sum(v.get('independent_witness_verified') is False for v in runs.values())
    result = dict(stage=44, created_utc=datetime.now(timezone.utc).isoformat(),
        protocol_sha256=sha(PROTOCOL), primary_results_unchanged=sha(PRIMARY_R) == protocol['primary_results_sha256'],
        name=case['name'], runs=runs, completed_answers_agree=len(answers) <= 1,
        verified_witnesses=sum(v.get('independent_witness_verified', False) for v in runs.values()),
        status='PASS' if len(answers) <= 1 and not validation_failures and all(v['status'] != 'error' for v in runs.values()) else 'FAIL',
        scope='Posthoc exact correctness only; primary 2-second benchmark remains unchanged')
    RESULT.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'runs'}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('phase', choices=('prepare', 'run'))
    args = parser.parse_args()
    prepare() if args.phase == 'prepare' else run()
