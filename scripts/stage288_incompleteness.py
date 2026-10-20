#!/usr/bin/env python3
"""Stage 288, steps 3-4: reference incompleteness on real reference trees.

Four exact answers per (case, removal regime, removal level):
  A  minimum strict source count on the FULL observed set            (reference)
  B  free pendant, observed subset, Lemma-3 gate marks RETAINED      (theory: B == A)
  C  naive subsetting: induced backbone, degree-two vertices suppressed,
     marks DISCARDED (marked=())                                     (practitioner)
  D  prescribed pendant with arbitrary missing references            (W[1]-hard regime)
  D1 prescribed pendant, exactly one missing reference (frozen FPT solver)

Phases
  prepare  register grid, input hashes, code hashes
  census   solver-free combinatorial census of the Lemma-3 mark set over ALL
           exactly-explainable bank cases and every removal draw
  control  synthetic d-arm star positive control: the paper's mark-necessity
           example, which the measurement apparatus must detect
  run      the four-way solver comparison on the registered panel

Nothing pre-existing is modified.
"""
import argparse
import hashlib
import json
import platform
import sys
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from fractions import Fraction as F
from pathlib import Path

import stage288_common as common
from stage288_common import (ROOT, serial, sha, load_case, verify_witness, run_pool,
                             leaves_of, clade_removals, nearest_source_removal,
                             random_removal, predict)
from stage26_missing_leaf import observed_backbone
from stage41_partial_fpt import lift_mixed
from stage43_partial_fpt import solve_partial, solve_mixed
from stage43_fixed_one_missing import solve_fixed_one_missing
from check_stage43_fixed_one_missing import oracle

PROTOCOL = ROOT / 'results/stage288_incompleteness_protocol.json'
CENSUS = ROOT / 'results/stage288_mark_census.json'
CONTROL = ROOT / 'results/stage288_mark_control.json'
RESULT = ROOT / 'results/stage288_incompleteness.json'
STREAM = ROOT / 'results/stage288_incompleteness_runs.jsonl'

SEED = 2880003
LEVELS = ('one', '0.10', '0.25', '0.50', '0.75')
REGIMES = ('random', 'clade', 'nearest')
REPLICATES = {'random': 2, 'clade': 2, 'nearest': 1}
N_BINS = [(20, 29), (30, 49), (50, 79), (80, 129), (130, 249), (250, 361), (1000, 100000)]
K_VALUES = (2, 3, 5, 7, 10, 12)
PER_CELL = {2: 4, 3: 4, 5: 4, 7: 4, 10: 2, 12: 1}
BC_WALL, BC_AS = 30.0, 4 * 1024 ** 3
D_WALL, D_AS = 60.0, 4 * 1024 ** 3
D_LEVELS = ('one', '0.25', '0.50')
WORKERS = 20
CONTROL_ARMS = (3, 4, 5, 6, 7, 8, 9)


def removal_count(n, level):
    if level == 'one':
        return 1
    count = max(1, int(round(float(level) * n)))
    return count if n - count >= 3 else None


def removal_draws(case, level):
    """Every registered removal draw for one case and one level."""
    leaves = sorted(case['observed'])
    n = len(leaves)
    count = removal_count(n, level)
    if count is None:
        return []
    draws = []
    for rep in range(REPLICATES['random']):
        key = f"{SEED}|{case['name']}|{level}|random|{rep}"
        seed = int(hashlib.sha256(key.encode()).hexdigest()[:8], 16)
        draws.append(dict(regime='random', replicate=rep, removed=random_removal(leaves, count, seed),
                          detail={'seed': seed}))
    for rep, clade in enumerate(clade_removals(case['edges'], leaves, count, REPLICATES['clade'])):
        draws.append(dict(regime='clade', replicate=rep, removed=clade['removed'],
                          detail={'clade_vertex': clade['clade_vertex'],
                                  'complement': clade['complement']}))
    draws.append(dict(regime='nearest', replicate=0,
                      removed=nearest_source_removal(case['edges'], case['generating_sources'],
                                                     leaves, count),
                      detail={}))
    for draw in draws:
        draw['n_removed'] = len(draw['removed'])
        draw['n_retained'] = n - len(draw['removed'])
        draw['level'] = level
    return [d for d in draws if d['n_retained'] >= 3]


def mark_report(edges, retained):
    """Lemma-3 gate marks for an observed subset; pure combinatorics."""
    backbone = observed_backbone(edges, sorted(retained))
    violations = backbone['saturation_violations']
    return dict(n_marks=len(violations),
                mark_gate_degrees=sorted(v['original_degree'] for v in violations),
                mark_hull_degrees=sorted(v['hull_degree'] for v in violations),
                backbone_edges=len(backbone['edges']),
                hull_vertices=len(backbone['hull_vertices'])), backbone


def panel():
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
        take = min(PER_CELL[key[2]], len(members))
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
    assert not PROTOCOL.exists(), 'Never overwrite frozen evidence'
    cases = panel()
    protocol = dict(
        stage=288, track='reference_incompleteness',
        registered_utc=datetime.now(timezone.utc).isoformat(), seed=SEED,
        purpose=('Test on real reference trees whether the Lemma-3 budget-preserving '
                 'projection holds when references are missing, and whether discarding the '
                 'branch-gate marks - what a practitioner who simply subsets a reference tree '
                 'does - inflates the inferred minimum source count.'),
        variants={
            'A_full': 'solve_partial on the full observed leaf set; the reference answer.',
            'B_marks_retained': 'solve_partial on the retained subset; marks are the backbone '
                                'saturation violations, retained automatically.',
            'C_naive_subset': 'solve_mixed on the induced suppressed backbone with marked=(); '
                              'the naive practitioner. Witnesses are lifted to the original '
                              'tree and verified there.',
            'D_prescribed_pendant': 'exact original-edge enumeration oracle with the pendant '
                                    'prescribed at the generating value and an arbitrary '
                                    'missing-reference subset; the W[1]-hard regime.',
            'D1_fixed_one_missing': 'frozen FPT prescribed-pendant solver, applicable only when '
                                    'exactly one reference is removed.'},
        levels=list(LEVELS), regimes=list(REGIMES), replicates=REPLICATES,
        removal_rules=('one = exactly one reference removed; a decimal level removes '
                       'round(level*N) references and is skipped when fewer than three would '
                       'remain. random draws without replacement from a registered seed; clade '
                       'removes the whole-clade leaf set whose size is closest to the target, '
                       'taking the two closest distinct sets; nearest removes the references '
                       'with the smallest exact tree distance to any true source point.'),
        caps='cap = K_generating for every variant of a case.',
        budgets=dict(bc_wall_seconds=BC_WALL, bc_address_space_bytes=BC_AS,
                     d_wall_seconds=D_WALL, d_address_space_bytes=D_AS,
                     d_levels=list(D_LEVELS), workers=WORKERS),
        per_cell=PER_CELL,
        budget_calibration=('A solver-free structural census of the bank was run before this '
                            'protocol, and a 12-case timing pilot was run to size the wall '
                            'limits. No solver was modified. Per-cell counts are reduced at '
                            'K>=10 because the pilot showed the structural solver does not '
                            'complete at cap 10 within a 60 s wall limit; those cells are kept '
                            'so the censoring is reported rather than hidden.'),
        prediction=('Lemma 3 predicts B == A at every level and regime. The four-leaf-star '
                    'example predicts C > B only when the mark set is non-empty, which on a '
                    'tree of maximum internal degree three is impossible. The registered census '
                    'phase measures how often the mark set is non-empty on real trees.'),
        analysis=('Retain every job including timeouts and errors; never impute. Count '
                  'comparisons use only (case, draw) pairs where the compared variants both '
                  'completed. Report the full distribution of C - A, not only extremes.'),
        limitations=[
            'Exactly-explainable families only; noisy families have no exact minimum to compare.',
            'A cap of K_generating means an inflation beyond K is observed as infeasibility '
            'within cap, not as a numeric count.',
            'Single run per (case, draw, variant): the measured quantity is an exact integer, '
            'not a timing.',
            'No claim about recovery of the generating support.'],
        cases=cases,
        script_hashes={n: sha(ROOT / 'scripts' / n)
                       for n in common.dependency_closure('stage288_incompleteness.py')},
        bank_protocol_sha256=sha(ROOT / 'results/stage287_instance_bank_protocol.json'),
        machine=dict(platform=platform.platform(), python=sys.version))
    PROTOCOL.write_text(json.dumps(serial(protocol), indent=2) + '\n')
    print(json.dumps(dict(protocol=str(PROTOCOL), sha256=sha(PROTOCOL), panel_cases=len(cases),
                          cells=len(Counter(c['cell'] for c in cases))), indent=2))


# --------------------------------------------------------------- census

def census():
    """Solver-free: how often is the Lemma-3 mark set non-empty on real trees?"""
    assert not CENSUS.exists(), 'Never overwrite frozen evidence'
    index = [e for e in common.bank_index() if e['family'].startswith('exact')]
    started = time.perf_counter()
    draws, per_case, degree_hist = [], [], Counter()
    for position, entry in enumerate(index):
        case = load_case(entry)
        adjacency = common.adjacency_lists(case['edges'])
        internal = [len(adjacency[v]) for v in adjacency if len(adjacency[v]) > 1]
        degree_hist.update(internal)
        row = dict(name=case['name'], family=case['family'], N=case['reference_set_size'],
                   K=case['K_generating'], max_internal_degree=max(internal),
                   n_internal_degree_at_least_4=sum(1 for d in internal if d >= 4),
                   draws_with_marks=0, draws=0, max_marks=0)
        for level in LEVELS:
            for draw in removal_draws(case, level):
                report, _ = mark_report(case['edges'], set(case['observed']) - set(draw['removed']))
                record = dict(name=case['name'], family=case['family'], N=row['N'], K=row['K'],
                              level=level, regime=draw['regime'], replicate=draw['replicate'],
                              n_removed=draw['n_removed'], n_retained=draw['n_retained'], **report)
                draws.append(record)
                row['draws'] += 1
                row['draws_with_marks'] += int(report['n_marks'] > 0)
                row['max_marks'] = max(row['max_marks'], report['n_marks'])
        per_case.append(row)
        if (position + 1) % 100 == 0:
            print(f'  census {position + 1}/{len(index)} '
                  f'({time.perf_counter() - started:.0f}s)', flush=True)
    summary = dict(
        exact_cases=len(per_case), total_draws=len(draws),
        draws_with_non_empty_mark_set=sum(d['n_marks'] > 0 for d in draws),
        cases_with_any_marked_draw=sum(r['draws_with_marks'] > 0 for r in per_case),
        cases_with_an_internal_degree_at_least_4=sum(r['max_internal_degree'] >= 4 for r in per_case),
        internal_degree_histogram=dict(sorted(degree_hist.items())),
        by_regime={r: dict(draws=sum(d['regime'] == r for d in draws),
                           with_marks=sum(d['regime'] == r and d['n_marks'] > 0 for d in draws))
                   for r in REGIMES},
        by_level={l: dict(draws=sum(d['level'] == l for d in draws),
                          with_marks=sum(d['level'] == l and d['n_marks'] > 0 for d in draws))
                  for l in LEVELS},
        marked_gate_degree_histogram=dict(sorted(Counter(
            g for d in draws for g in d['mark_gate_degrees']).items())),
        note=('A gate mark needs hull degree >= 3 and original degree > hull degree, so a tree '
              'whose internal vertices all have degree 3 can never produce one; on such trees '
              'C and B are the same computation.'))
    CENSUS.write_text(json.dumps(serial(dict(
        stage=288, track='mark_census', created_utc=datetime.now(timezone.utc).isoformat(),
        protocol_sha256=sha(PROTOCOL), summary=summary, per_case=per_case, draws=draws)),
        indent=2) + '\n')
    print(json.dumps(summary, indent=2)[:4000])


# --------------------------------------------------------------- control

def star_case(arms, lengths=None):
    """Unit star with `arms`+1 arms; the omitted arm carries the single source."""
    lengths = lengths or [F(1)] * (arms + 1)
    centre = 1000
    edges = [(centre, i, lengths[i]) for i in range(arms + 1)]
    sources = [dict(edge=arms, p=F(1), x=F(1, 2))]
    pendant = F(1)
    leaves = list(range(arms + 1))
    observed = predict(edges, sources, pendant, leaves)
    return dict(name=f'star_control_arms{arms}' + ('' if lengths[0] == lengths[1] else '_unequal'),
                edges=edges, observed=observed, cap=arms + 1, generating_sources=sources,
                generating_pendant=pendant, omitted=arms)


def control():
    """The paper's mark-necessity example: the apparatus must detect inflation here."""
    assert not CONTROL.exists(), 'Never overwrite frozen evidence'
    rows = []
    cases = [star_case(a) for a in CONTROL_ARMS]
    cases.append(star_case(4, [F(1), F(3, 2), F(2), F(5, 2), F(1)]))
    for case in cases:
        retained = sorted(r for r in case['observed'] if r != case['omitted'])
        subset = {r: case['observed'][r] for r in retained}
        report, backbone = mark_report(case['edges'], retained)
        full = solve_partial(case['edges'], case['observed'], cap=case['cap'], seconds=None,
                             max_states=None, max_bases=None, max_faces=None)
        kept = solve_partial(case['edges'], subset, cap=case['cap'], seconds=None,
                             max_states=None, max_bases=None, max_faces=None)
        naive = solve_mixed(backbone['edges'], subset, marked=(), cap=case['cap'], seconds=None,
                            max_states=None, max_bases=None, max_faces=None)
        row = dict(name=case['name'], observed_arms=len(retained), total_arms=len(case['observed']),
                   cap=case['cap'], **report,
                   A_full=full.get('minimum_sources_within_cap'),
                   B_marks_retained=kept.get('minimum_sources_within_cap'),
                   C_naive_subset=naive.get('minimum_sources_within_cap'),
                   B_marked_vertices=kept.get('marked_vertices'))
        for label, answer, edges, obs in (
                ('A_full', full, case['edges'], case['observed']),
                ('B_marks_retained', kept, case['edges'], subset)):
            if 'original_witness' in answer:
                row[label + '_witness'] = verify_witness(edges, obs, case['cap'],
                                                         answer['original_witness'])['verified']
        if naive.get('feasible'):
            lifted = lift_mixed(backbone, naive['atoms'], naive['pendant'])
            row['C_naive_subset_witness'] = verify_witness(case['edges'], subset, case['cap'],
                                                           lifted)['verified']
        row['inflation_C_minus_A'] = (None if row['C_naive_subset'] is None or row['A_full'] is None
                                      else row['C_naive_subset'] - row['A_full'])
        rows.append(row)
        print(row, flush=True)
    CONTROL.write_text(json.dumps(serial(dict(
        stage=288, track='mark_control', created_utc=datetime.now(timezone.utc).isoformat(),
        protocol_sha256=sha(PROTOCOL),
        purpose=('Synthetic positive control for the step 3-4 apparatus: the manuscript '
                 'four-leaf-unit-star mark-necessity example, generalised to d observed arms '
                 'and to unequal arm lengths. Synthetic, not a real reference tree.'),
        rows=rows)), indent=2) + '\n')


# --------------------------------------------------------------- solver run

def worker(spec):
    case, variant, cap = spec['case'], spec['variant'], spec['cap']
    observed = spec['observed']
    start = time.perf_counter()
    if variant in ('A_full', 'B_marks_retained'):
        raw = solve_partial(case['edges'], observed, cap=cap, seconds=None, max_states=None,
                            max_bases=None, max_faces=None)
        answer = dict(status=raw['status'], feasible=raw.get('feasible'),
                      minimum=raw.get('minimum_sources_within_cap'),
                      reason=raw.get('reason'), marked_vertices=raw.get('marked_vertices'),
                      backbone_edges=len(raw.get('backbone_edges', [])),
                      statistics=raw['statistics'])
        if 'original_witness' in raw:
            answer['witness'] = raw['original_witness']
            answer['witness_frame'] = 'original'
    elif variant == 'C_naive_subset':
        backbone = observed_backbone(case['edges'], sorted(observed))
        raw = solve_mixed(backbone['edges'], observed, marked=(), cap=cap, seconds=None,
                          max_states=None, max_bases=None, max_faces=None)
        answer = dict(status=raw['status'], feasible=raw.get('feasible'),
                      minimum=raw.get('minimum_sources_within_cap'), reason=raw.get('reason'),
                      marked_vertices=[], backbone_edges=len(backbone['edges']),
                      discarded_marks=[v['vertex'] for v in backbone['saturation_violations']],
                      statistics=raw['statistics'])
        if raw.get('feasible'):
            answer['witness'] = lift_mixed(backbone, raw['atoms'], raw['pendant'])
            answer['witness_frame'] = 'original'
    elif variant == 'D_prescribed_pendant':
        raw = oracle(dict(edges=case['edges'], observed=observed, cap=cap,
                          pendant=case['generating_pendant']))
        answer = dict(status='complete', feasible=raw['minimum'] is not None,
                      minimum=raw['minimum'])
        if 'witness' in raw:
            answer['witness'] = raw['witness']
            answer['witness_frame'] = 'original'
    elif variant == 'D1_fixed_one_missing':
        raw = solve_fixed_one_missing(case['edges'], observed, case['generating_pendant'],
                                      cap=cap, seconds=None, max_states=None, max_bases=None,
                                      max_faces=None)
        answer = dict(status=raw['status'], feasible=raw.get('feasible'),
                      minimum=raw.get('minimum_sources_within_cap'), reason=raw.get('reason'),
                      statistics=raw['statistics'])
        if 'sources' in raw:
            answer['witness'] = {'sources': raw['sources'], 'pendant': raw['pendant']}
            answer['witness_frame'] = 'original'
    else:
        raise ValueError(variant)
    answer['solver_wall_seconds'] = time.perf_counter() - start
    return answer


def make_post(cases_by_name, observed_by_key):
    def post(job, answer):
        if 'witness' not in answer:
            return
        meta = job['meta']
        case = cases_by_name[meta['name']]
        observed = observed_by_key[meta['observed_key']]
        check = verify_witness(case['edges'], observed, meta['cap'], answer['witness'])
        answer['independent_witness_verification'] = check
        answer['independent_witness_verified'] = check['verified']
        answer['witness_count_matches_minimum'] = check.get('source_count') == answer.get('minimum')
        if meta['variant'].startswith('D'):
            answer['prescribed_pendant_respected'] = (
                F(answer['witness']['pendant']) == case['generating_pendant'])
    return post


def run():
    assert not RESULT.exists(), 'Never overwrite measured runs'
    protocol = json.loads(PROTOCOL.read_text())
    for name, expected in protocol['script_hashes'].items():
        assert sha(ROOT / 'scripts' / name) == expected, f'Frozen source changed: {name}'
    cases = {}
    for entry in protocol['cases']:
        cases[entry['name']] = load_case(entry)
        cases[entry['name']]['cell'] = entry['cell']
    print(f'loaded {len(cases)} verified bank cases', flush=True)

    observed_by_key, bc_jobs, d_jobs = {}, [], []
    for entry in protocol['cases']:
        case = cases[entry['name']]
        cap = case['cap']
        key = case['name'] + '|full'
        observed_by_key[key] = case['observed']
        base = dict(name=case['name'], family=case['family'], cell=case['cell'],
                    N=case['reference_set_size'], K=case['K_generating'], cap=cap)
        bc_jobs.append(dict(meta=dict(base, variant='A_full', level='none', regime='none',
                                      replicate=0, n_removed=0, n_retained=case['reference_set_size'],
                                      n_marks=0, observed_key=key),
                            spec=dict(case=case, variant='A_full', cap=cap,
                                      observed=case['observed'])))
        for level in LEVELS:
            for draw in removal_draws(case, level):
                retained = sorted(set(case['observed']) - set(draw['removed']))
                observed = {r: case['observed'][r] for r in retained}
                report, _ = mark_report(case['edges'], retained)
                key = f"{case['name']}|{level}|{draw['regime']}|{draw['replicate']}"
                observed_by_key[key] = observed
                meta = dict(base, level=level, regime=draw['regime'], replicate=draw['replicate'],
                            n_removed=draw['n_removed'], n_retained=draw['n_retained'],
                            n_marks=report['n_marks'],
                            mark_gate_degrees=report['mark_gate_degrees'],
                            observed_key=key)
                for variant in ('B_marks_retained', 'C_naive_subset'):
                    bc_jobs.append(dict(meta=dict(meta, variant=variant),
                                        spec=dict(case=case, variant=variant, cap=cap,
                                                  observed=observed)))
                if level in D_LEVELS and draw['replicate'] == 0:
                    d_jobs.append(dict(meta=dict(meta, variant='D_prescribed_pendant'),
                                       spec=dict(case=case, variant='D_prescribed_pendant',
                                                 cap=cap, observed=observed)))
                    if level == 'one':
                        d_jobs.append(dict(meta=dict(meta, variant='D1_fixed_one_missing'),
                                           spec=dict(case=case, variant='D1_fixed_one_missing',
                                                     cap=cap, observed=observed)))
    post = make_post(cases, observed_by_key)
    print(f'B/C jobs {len(bc_jobs)}, D jobs {len(d_jobs)}', flush=True)
    started = time.perf_counter()
    records = run_pool(bc_jobs, worker, BC_WALL, BC_AS, WORKERS, STREAM, post=post, progress=250)
    print(f'B/C done in {(time.perf_counter() - started) / 60:.1f} min', flush=True)
    started = time.perf_counter()
    records += run_pool(d_jobs, worker, D_WALL, D_AS, WORKERS, STREAM, post=post, progress=50)
    print(f'D done in {(time.perf_counter() - started) / 60:.1f} min', flush=True)
    RESULT.write_text(json.dumps(serial(dict(
        stage=288, track='reference_incompleteness',
        created_utc=datetime.now(timezone.utc).isoformat(), protocol_sha256=sha(PROTOCOL),
        script_hashes=protocol['script_hashes'], records=records)), indent=2) + '\n')
    print('wrote', RESULT, 'records', len(records))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('phase', choices=('prepare', 'census', 'control', 'run'))
    args = parser.parse_args()
    {'prepare': prepare, 'census': census, 'control': control, 'run': run}[args.phase]()
