#!/usr/bin/env python3
"""Stage 288 addendum: a direct falsification test of the Lemma-3 projection.

Lemma 3 states that for a FIXED observed subset R, the minimum strict
original-edge source count on the original tree T equals the minimum mixed
edge/marked-vertex atom count on the reduced backbone (B, M). It does NOT say
that the answer is unchanged when R shrinks - dropping observations drops
equations, so the minimum can only fall or stay equal.

This script tests the lemma as stated, on real reference trees with real
missing references:
  E_original_enumeration  original-edge support enumeration on T, observed = R
  B_marks_retained        solve_partial on T, observed = R (reduced backbone)
Lemma 3 predicts E == B exactly, including infeasibility.

The panel is restricted to budgets and sizes where the enumeration oracle can
finish; that restriction is on cost only and is declared in the protocol.
"""
import argparse
import json
import platform
import random
import sys
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from statistics import median

import stage288_common as common
from stage288_common import ROOT, serial, sha, load_case, verify_witness, run_pool
from stage288_incompleteness import removal_draws, mark_report
from stage42_benchmark import baseline
from stage43_partial_fpt import solve_partial

PROTOCOL = ROOT / 'results/stage288_lemma3_protocol.json'
RESULT = ROOT / 'results/stage288_lemma3.json'
STREAM = ROOT / 'results/stage288_lemma3_runs.jsonl'

SEED = 2880005
MAX_N = 40
K_VALUES = (2, 3)
PER_K = 10
LEVELS = ('one', '0.25', '0.50')
WALL, ADDRESS_SPACE, WORKERS = 120.0, 4 * 1024 ** 3, 20


def panel():
    index = [e for e in common.bank_index()
             if e['family'].startswith('exact') and e['reference_set_size'] <= MAX_N
             and e['K_generating'] in K_VALUES]
    selected = []
    for k in K_VALUES:
        members = sorted([e for e in index if e['K_generating'] == k],
                         key=lambda e: (e['reference_set_size'], e['name']))
        take = min(PER_K, len(members))
        picks = [members[0]] if take == 1 else [
            members[round(j * (len(members) - 1) / (take - 1))] for j in range(take)]
        seen = set()
        for entry in picks:
            if entry['name'] not in seen:
                seen.add(entry['name'])
                selected.append(dict(entry, cell=f'K{k}'))
    return selected


def prepare():
    assert not PROTOCOL.exists(), 'Never overwrite frozen evidence'
    cases = panel()
    protocol = dict(
        stage=288, track='lemma3_projection_test',
        registered_utc=datetime.now(timezone.utc).isoformat(), seed=SEED,
        purpose=('Falsification test of the budget-preserving projection lemma on real '
                 'reference trees with real missing references: original-edge enumeration on '
                 'the original tree versus the reduced-backbone structural solver, on the same '
                 'observed subset. The lemma predicts exact agreement including infeasibility.'),
        variants={'E_original_enumeration': 'stage42_benchmark.baseline on the original tree with '
                                            'the retained observed subset and a free pendant.',
                  'B_marks_retained': 'stage43_partial_fpt.solve_partial on the same subset.'},
        selection=(f'Exactly-explainable bank cases with reference-set size <= {MAX_N} and '
                   f'K in {K_VALUES}, up to {PER_K} per K at evenly spaced ranks after sorting '
                   'by (reference_set_size, name). The size and budget restriction is a COST '
                   'restriction: the enumeration oracle visits C(|E|,k) supports and does not '
                   'finish beyond this range within any budget available here.'),
        removal_draws=(f'The stage 288 incompleteness removal draws at levels {LEVELS}, all '
                       'regimes and replicates, plus the full observed set as a control.'),
        budgets=dict(wall_seconds=WALL, address_space_bytes=ADDRESS_SPACE, workers=WORKERS,
                     repetitions=1),
        analysis=('Every job retained including timeouts. Agreement is assessed only where both '
                  'variants completed. Witnesses are verified against the original tree.'),
        limitations=['Small reference sets and budgets 2-3 only; no evidence about larger K.',
                     'A single run per job: the measured quantity is an exact integer.',
                     'Agreement supports the implementations and the lemma on these inputs; it '
                     'is not a proof.'],
        cases=cases,
        script_hashes={n: sha(ROOT / 'scripts' / n)
                       for n in common.dependency_closure('stage288_lemma3.py')},
        bank_protocol_sha256=sha(ROOT / 'results/stage287_instance_bank_protocol.json'),
        machine=dict(platform=platform.platform(), python=sys.version))
    PROTOCOL.write_text(json.dumps(serial(protocol), indent=2) + '\n')
    print(json.dumps(dict(protocol=str(PROTOCOL), sha256=sha(PROTOCOL), cases=len(cases)), indent=2))


def worker(spec):
    case, variant, cap, observed = spec['case'], spec['variant'], spec['cap'], spec['observed']
    start = time.perf_counter()
    if variant == 'E_original_enumeration':
        raw = baseline(dict(edges=case['edges'], observed=observed, cap=cap))
        answer = dict(status=raw['status'], feasible=raw['feasible'], minimum=raw['minimum'],
                      statistics=raw['statistics'])
        if 'witness' in raw:
            answer['witness'] = raw['witness']
    else:
        raw = solve_partial(case['edges'], observed, cap=cap, seconds=None, max_states=None,
                            max_bases=None, max_faces=None)
        answer = dict(status=raw['status'], feasible=raw.get('feasible'),
                      minimum=raw.get('minimum_sources_within_cap'), reason=raw.get('reason'),
                      marked_vertices=raw.get('marked_vertices'), statistics=raw['statistics'])
        if 'original_witness' in raw:
            answer['witness'] = raw['original_witness']
    answer['solver_wall_seconds'] = time.perf_counter() - start
    return answer


def run():
    assert not RESULT.exists(), 'Never overwrite measured runs'
    protocol = json.loads(PROTOCOL.read_text())
    for name, expected in protocol['script_hashes'].items():
        assert sha(ROOT / 'scripts' / name) == expected, f'Frozen source changed: {name}'
    cases = {e['name']: load_case(e) for e in protocol['cases']}
    jobs, observed_by_key = [], {}
    for entry in protocol['cases']:
        case = cases[entry['name']]
        base = dict(name=case['name'], family=case['family'], N=case['reference_set_size'],
                    K=case['K_generating'], cap=case['cap'])
        draws = [dict(level='none', regime='none', replicate=0, removed=[],
                      n_removed=0, n_retained=case['reference_set_size'])]
        for level in LEVELS:
            draws.extend(removal_draws(case, level))
        for draw in draws:
            retained = sorted(set(case['observed']) - set(draw['removed']))
            observed = {r: case['observed'][r] for r in retained}
            report, _ = mark_report(case['edges'], retained)
            key = f"{case['name']}|{draw['level']}|{draw['regime']}|{draw['replicate']}"
            observed_by_key[key] = observed
            meta = dict(base, level=draw['level'], regime=draw['regime'],
                        replicate=draw['replicate'], n_removed=draw['n_removed'],
                        n_retained=draw['n_retained'], n_marks=report['n_marks'],
                        observed_key=key)
            for variant in ('E_original_enumeration', 'B_marks_retained'):
                jobs.append(dict(meta=dict(meta, variant=variant),
                                 spec=dict(case=case, variant=variant, cap=case['cap'],
                                           observed=observed)))
    rng = random.Random(SEED)
    rng.shuffle(jobs)

    def post(job, answer):
        if 'witness' not in answer:
            return
        case = cases[job['meta']['name']]
        check = verify_witness(case['edges'], observed_by_key[job['meta']['observed_key']],
                               job['meta']['cap'], answer['witness'])
        answer['independent_witness_verification'] = check
        answer['independent_witness_verified'] = check['verified']

    print(f'{len(jobs)} jobs, wall {WALL}s, {WORKERS} workers', flush=True)
    records = run_pool(jobs, worker, WALL, ADDRESS_SPACE, WORKERS, STREAM, post=post, progress=50)

    index = {(r['name'], r['level'], r['regime'], r['replicate'], r['variant']): r for r in records}
    rows = []
    for key, rec in index.items():
        if rec['variant'] != 'B_marks_retained':
            continue
        other = index.get(key[:4] + ('E_original_enumeration',))
        rows.append(dict(name=key[0], N=rec['N'], K=rec['K'], level=key[1], regime=key[2],
                         replicate=key[3], n_removed=rec['n_removed'],
                         n_retained=rec['n_retained'], n_marks=rec['n_marks'],
                         B_status=rec['result']['status'], B_min=rec['result'].get('minimum'),
                         E_status=other['result']['status'] if other else None,
                         E_min=other['result'].get('minimum') if other else None))
    comparable = [r for r in rows if r['B_status'] == 'complete' and r['E_status'] == 'complete']
    summary = dict(
        jobs=len(records), draws=len(rows),
        status_by_variant={v: dict(Counter(r['result']['status'] for r in records
                                           if r['variant'] == v))
                           for v in sorted({r['variant'] for r in records})},
        comparable_draws=len(comparable),
        agree=sum(r['B_min'] == r['E_min'] for r in comparable),
        disagreements=[r for r in comparable if r['B_min'] != r['E_min']][:40],
        agreement_by_level={l: dict(comparable=sum(r['level'] == l for r in comparable),
                                    agree=sum(r['level'] == l and r['B_min'] == r['E_min']
                                              for r in comparable))
                            for l in sorted({r['level'] for r in rows})},
        agreement_by_regime={g: dict(comparable=sum(r['regime'] == g for r in comparable),
                                     agree=sum(r['regime'] == g and r['B_min'] == r['E_min']
                                               for r in comparable))
                             for g in sorted({r['regime'] for r in rows})},
        minimum_distribution=dict(Counter(str(r['B_min']) for r in comparable)),
        median_cpu_seconds={v: median([float(r['result']['child_cpu_seconds']) for r in records
                                       if r['variant'] == v and r['result']['status'] == 'complete']
                                      or [0]) or None
                            for v in sorted({r['variant'] for r in records})},
        witness_verified=sum(bool(r['result'].get('independent_witness_verified'))
                             for r in records if 'independent_witness_verified' in r['result']),
        witness_failed=sum(not r['result']['independent_witness_verified'] for r in records
                           if 'independent_witness_verified' in r['result']))
    RESULT.write_text(json.dumps(serial(dict(
        stage=288, track='lemma3_projection_test',
        created_utc=datetime.now(timezone.utc).isoformat(), protocol_sha256=sha(PROTOCOL),
        summary=summary, rows=rows, records=records)), indent=2) + '\n')
    print(json.dumps(serial(summary), indent=2)[:5000])


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('phase', choices=('prepare', 'run'))
    args = parser.parse_args()
    prepare() if args.phase == 'prepare' else run()
