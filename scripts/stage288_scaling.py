#!/usr/bin/env python3
"""Stage 288, steps 1-2: scaling of the frozen structural solver against the
frozen original-edge enumeration baseline, on the frozen stage 287 real-data
instance bank.

prepare  registers the grid, input hashes and code hashes before any solve.
run      executes every job in a fresh forked child with a wall and an
         address-space limit, retaining timeouts; witnesses are verified in the
         parent against the original tree, outside child timing.
analyze  completion rates, joint-completion medians, N- and K-dependence.

Nothing pre-existing is modified.
"""
import argparse
import json
import platform
import random
import sys
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from fractions import Fraction as F
from pathlib import Path
from statistics import median

import stage288_common as common
from stage288_common import ROOT, serial, sha, load_case, verify_witness, run_pool
from stage42_benchmark import baseline
from stage43_partial_fpt import solve_partial

PROTOCOL = ROOT / 'results/stage288_scaling_protocol.json'
RESULT = ROOT / 'results/stage288_scaling.json'
STREAM = ROOT / 'results/stage288_scaling_runs.jsonl'

SEED = 2880001
N_BINS = [(20, 29), (30, 49), (50, 79), (80, 129), (130, 249), (250, 361), (1000, 100000)]
K_VALUES = (2, 3, 5, 7, 10, 12)
PER_CELL = 6
METHODS = ('original_support', 'stage43_equality_reduced')
TIER1 = dict(repetitions=2, wall_seconds=60.0, address_space_bytes=4 * 1024 ** 3)
TIER2 = dict(repetitions=1, wall_seconds=600.0, address_space_bytes=8 * 1024 ** 3,
             applies_to_K_at_least=7, method='stage43_equality_reduced')
TIER3 = dict(repetitions=1, wall_seconds=600.0, address_space_bytes=16 * 1024 ** 3,
             caps=(1, 2, 3), family='exact_16k')
WORKERS = 20


def grid():
    """Deterministic stratified selection over (reference-set-size bin, K)."""
    index = [e for e in common.bank_index() if e['family'].startswith('exact')]
    cells = defaultdict(list)
    for entry in index:
        n = entry['reference_set_size']
        for low, high in N_BINS:
            if low <= n <= high:
                if entry['K_generating'] in K_VALUES:
                    cells[(low, high, entry['K_generating'])].append(entry)
                break
    selected = []
    for key in sorted(cells):
        members = sorted(cells[key], key=lambda e: (e['reference_set_size'], e['name']))
        take = min(PER_CELL, len(members))
        picks = [members[0]] if take == 1 else [
            members[round(j * (len(members) - 1) / (take - 1))] for j in range(take)]
        seen = set()
        for entry in picks:
            if entry['name'] in seen:
                continue
            seen.add(entry['name'])
            selected.append(dict(entry, cell=f'N{key[0]}-{key[1]}_K{key[2]}'))
    return selected


def prepare():
    assert not PROTOCOL.exists() and not RESULT.exists() and not STREAM.exists(), \
        'Never overwrite frozen evidence'
    cases = grid()
    counts = Counter(c['cell'] for c in cases)
    protocol = dict(
        stage=288, track='scaling', registered_utc=datetime.now(timezone.utc).isoformat(),
        seed=SEED,
        purpose=('Measure how the frozen structural solver and the frozen original-edge '
                 'enumeration baseline scale with reference-set size N and source budget K '
                 'on real reference trees; locate the regime where enumeration stops '
                 'completing while the structural method still finishes. Not a comparison '
                 'with the published DecoDiPhy software.'),
        selection=('All exactly-explainable frozen bank cases binned by reference-set size '
                   f'{N_BINS} and generating K {K_VALUES}; up to {PER_CELL} per cell, taken at '
                   'evenly spaced ranks after sorting the cell by (reference_set_size, name). '
                   'No solver was run before this protocol was written; no outcome-dependent '
                   'filtering.'),
        tiers=dict(
            tier1=dict(TIER1, methods=list(METHODS),
                       description='Primary paired comparison, both methods, two fresh runs.'),
            tier2=dict(TIER2, description=('Prespecified structural-only frontier probe at a '
                                           'longer wall limit for budgets where tier 1 is '
                                           'expected to censor both methods. Single run; not '
                                           'paired, never substituted into tier-1 statistics.')),
            tier3=dict(TIER3, methods=list(METHODS),
                       description=('Prespecified large-N probe on the single 15,941-leaf real '
                                    'case at reduced caps. Answers feasibility within the '
                                    'reduced cap, which is NOT the true minimum of that case.'))),
        execution=('Jobs run in forked children, at most %d concurrently, OMP/OPENBLAS/MKL '
                   'threads pinned to 1. Order is a single seeded shuffle registered here, so '
                   'repetitions and methods of one case are not systematically co-scheduled. '
                   'Child CPU time (ru_utime+ru_stime via wait4) is the primary cost measure; '
                   'process wall time is recorded alongside so co-scheduling inflation can be '
                   'audited.' % WORKERS),
        analysis=('Retain every job including timeouts and errors. Timeouts stay UNKNOWN and are '
                  'never imputed. Medians and paired ratios use only cases where both methods '
                  'completed both tier-1 repetitions. Witness verification recomputes the leaf '
                  'equations on the original tree in exact arithmetic in the parent process.'),
        limitations=[
            'Minimum strict source count, not recovery of the generating mixture.',
            'Two repetitions do not justify confidence intervals.',
            'The enumeration baseline is a transparent exact oracle, not DecoDiPhy.',
            'Real exactly-explainable mixtures exist at N=20-361 and at N=15,941 only; the '
            'N-fit has no real anchor between 361 and 15,941.',
            'Concurrent execution on a shared machine; CPU time is reported for that reason.'],
        n_bins=N_BINS, k_values=list(K_VALUES), per_cell=PER_CELL,
        cases=[{k: v for k, v in c.items()} for c in cases],
        cell_counts=dict(sorted(counts.items())),
        script_hashes={n: sha(ROOT / 'scripts' / n)
                       for n in common.dependency_closure('stage288_scaling.py')},
        bank_protocol_sha256=sha(ROOT / 'results/stage287_instance_bank_protocol.json'),
        machine=dict(platform=platform.platform(), python=sys.version, cpu_count=len(os.sched_getaffinity(0))
                     if hasattr(__import__('os'), 'sched_getaffinity') else None))
    PROTOCOL.write_text(json.dumps(serial(protocol), indent=2) + '\n')
    print(json.dumps(dict(protocol=str(PROTOCOL), sha256=sha(PROTOCOL), cases=len(cases),
                          cells=len(counts), cell_counts=dict(sorted(counts.items()))), indent=2))


def worker(spec):
    case, method, cap = spec['case'], spec['method'], spec['cap']
    start = time.perf_counter()
    if method == 'original_support':
        raw = baseline(dict(edges=case['edges'], observed=case['observed'], cap=cap))
        answer = dict(status=raw['status'], feasible=raw['feasible'], minimum=raw['minimum'],
                      statistics=raw['statistics'])
        if 'witness' in raw:
            answer['witness'] = raw['witness']
    else:
        raw = solve_partial(case['edges'], case['observed'], cap=cap, seconds=None,
                            max_states=None, max_bases=None, max_faces=None)
        answer = dict(status=raw['status'], feasible=raw.get('feasible'),
                      minimum=raw.get('minimum_sources_within_cap'),
                      statistics=raw['statistics'],
                      marked_vertices=raw.get('marked_vertices'),
                      backbone_edges=len(raw.get('backbone_edges', [])))
        if 'original_witness' in raw:
            answer['witness'] = raw['original_witness']
    answer['solver_wall_seconds'] = time.perf_counter() - start
    return answer


def make_post(cases_by_name):
    def post(job, answer):
        if 'witness' not in answer:
            return
        case = cases_by_name[job['meta']['name']]
        check = verify_witness(case['edges'], case['observed'], job['meta']['cap'],
                               answer['witness'])
        answer['independent_witness_verification'] = check
        answer['independent_witness_verified'] = check['verified']
        answer['witness_count_matches_minimum'] = (check.get('source_count') == answer.get('minimum'))
    return post


def run():
    assert not RESULT.exists(), 'Never overwrite measured runs'
    protocol = json.loads(PROTOCOL.read_text())
    for name, expected in protocol['script_hashes'].items():
        assert sha(ROOT / 'scripts' / name) == expected, f'Frozen source changed: {name}'
    assert sha(ROOT / 'results/stage287_instance_bank_protocol.json') == protocol['bank_protocol_sha256']
    cases = {}
    for entry in protocol['cases']:
        cases[entry['name']] = load_case(entry)
        cases[entry['name']]['cell'] = entry['cell']
    print(f'loaded {len(cases)} verified bank cases', flush=True)

    jobs = []
    for entry in protocol['cases']:
        case = cases[entry['name']]
        for repetition in range(1, TIER1['repetitions'] + 1):
            for method in METHODS:
                jobs.append(dict(
                    meta=dict(tier=1, name=case['name'], family=case['family'], cell=case['cell'],
                              N=case['reference_set_size'], K=case['K_generating'],
                              cap=case['cap'], method=method, repetition=repetition),
                    spec=dict(case=case, method=method, cap=case['cap'])))
    rng = random.Random(SEED)
    order = list(range(len(jobs)))
    rng.shuffle(order)
    jobs = [jobs[i] for i in order]
    print(f'tier 1: {len(jobs)} jobs, wall {TIER1["wall_seconds"]}s, {WORKERS} workers', flush=True)
    started = time.perf_counter()
    records = run_pool(jobs, worker, TIER1['wall_seconds'], TIER1['address_space_bytes'],
                       WORKERS, STREAM, post=make_post(cases), progress=40)
    print(f'tier 1 done in {(time.perf_counter() - started) / 60:.1f} min', flush=True)

    jobs2 = [dict(meta=dict(tier=2, name=c['name'], family=c['family'], cell=c['cell'],
                            N=c['reference_set_size'], K=c['K_generating'], cap=c['cap'],
                            method=TIER2['method'], repetition=1),
                  spec=dict(case=c, method=TIER2['method'], cap=c['cap']))
             for c in cases.values() if c['K_generating'] >= TIER2['applies_to_K_at_least']]
    print(f'tier 2: {len(jobs2)} jobs, wall {TIER2["wall_seconds"]}s', flush=True)
    started = time.perf_counter()
    records2 = run_pool(jobs2, worker, TIER2['wall_seconds'], TIER2['address_space_bytes'],
                        WORKERS, STREAM, post=make_post(cases), progress=10)
    print(f'tier 2 done in {(time.perf_counter() - started) / 60:.1f} min', flush=True)

    jobs3 = []
    for c in cases.values():
        if c['family'] != TIER3['family']:
            continue
        for cap in TIER3['caps']:
            for method in METHODS:
                jobs3.append(dict(meta=dict(tier=3, name=c['name'], family=c['family'],
                                            cell=c['cell'], N=c['reference_set_size'],
                                            K=c['K_generating'], cap=cap, method=method,
                                            repetition=1),
                                  spec=dict(case=c, method=method, cap=cap)))
    print(f'tier 3: {len(jobs3)} jobs, wall {TIER3["wall_seconds"]}s', flush=True)
    records3 = run_pool(jobs3, worker, TIER3['wall_seconds'], TIER3['address_space_bytes'],
                        min(WORKERS, 6), STREAM, post=make_post(cases), progress=2)

    out = dict(stage=288, track='scaling', created_utc=datetime.now(timezone.utc).isoformat(),
               protocol_sha256=sha(PROTOCOL), script_hashes=protocol['script_hashes'],
               records=records + records2 + records3)
    RESULT.write_text(json.dumps(serial(out), indent=2) + '\n')
    print('wrote', RESULT, 'records', len(out['records']))


import os  # noqa: E402  (used by prepare for affinity)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('phase', choices=('prepare', 'run'))
    args = parser.parse_args()
    prepare() if args.phase == 'prepare' else run()
