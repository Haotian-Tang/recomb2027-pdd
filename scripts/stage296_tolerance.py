#!/usr/bin/env python3
"""Stage 296: the eps-tolerant source budget - noise as the controlled variable.

Stage 294 could not separate "the exact-model downward bias is an artefact of
the exact model" from "DecoDiPhy is too noisy an estimator of k to see it".
This stage removes the tool from the loop and makes the tolerance eps the
controlled variable, so the question becomes a property of the MODEL.

PHASES
  prepare    register the panel, the grid, every input and code hash
  cross      Theta_k on the full panel and on every stage-288 removal draw
  probe      Theta_k at the full panel only, for K = 5, 7, 10 (scale ceiling)
  validate   eps = 0 agreement against the FROZEN exact machinery
  calibrate  real estimation error placed in the same units as eps

The removal regimes, levels, seeds and draw plan are imported from
stage288_incompleteness (which imports stage288_common), never re-derived, so
the cells are comparable with the stage-288 exact results.

The frozen structural localisation is used ONLY in `validate`, at eps = 0,
where it is sound.  Every eps > 0 decision comes from the exhaustive support
enumeration in stage296_common, decided by an exact rational LP.

Nothing pre-existing is modified.
"""
import argparse
import json
import platform
import sys
import time
from collections import defaultdict
from datetime import datetime, timezone
from fractions import Fraction as F
from math import comb
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import stage288_common as common
import stage296_common as s296
import stage288_incompleteness as s288
from stage288_common import ROOT, serial, sha, load_case, verify_witness, run_pool
from stage43_partial_fpt import solve_partial

OUT = ROOT / 'results/stage296'
PROTOCOL = OUT / 'stage296_protocol.json'
CROSS_STREAM = OUT / 'stage296_cross_runs.jsonl'
PROBE_STREAM = OUT / 'stage296_probe_runs.jsonl'
VALIDATE_STREAM = OUT / 'stage296_validate_runs.jsonl'
CALIBRATION = OUT / 'stage296_calibration.json'

WORKERS = 32
CROSS_WALL, CROSS_AS = 3600.0, 6 * 1024 ** 3
PROBE_WALL, PROBE_AS = 10800.0, 8 * 1024 ** 3
VALIDATE_WALL, VALIDATE_AS = 120.0, 6 * 1024 ** 3

CROSS_CELLS = (
    ('exact_10k', 3, (20, 29), 5), ('exact_10k', 3, (30, 49), 5),
    ('exact_10k', 3, (50, 79), 3), ('exact_10k', 3, (80, 129), 2),
    ('exact_biotrees', 3, (20, 29), 5), ('exact_biotrees', 3, (30, 49), 5),
    ('exact_biotrees', 3, (80, 129), 2),
    ('exact_biotrees', 2, (20, 29), 5), ('exact_biotrees', 2, (30, 49), 5),
    ('exact_biotrees', 2, (80, 129), 3), ('exact_biotrees', 2, (130, 361), 3),
)
PROBE_CELLS = (
    ('exact_10k', 5, (20, 29), 2), ('exact_10k', 5, (30, 49), 2),
    ('exact_10k', 5, (50, 79), 1), ('exact_biotrees', 5, (20, 29), 2),
    ('exact_biotrees', 5, (30, 49), 1), ('exact_biotrees', 7, (20, 29), 2),
    ('exact_biotrees', 7, (30, 49), 2), ('exact_10k', 10, (20, 29), 2),
    ('exact_10k', 10, (30, 49), 2), ('exact_biotrees', 10, (20, 29), 2),
)
PROBE_SUPPORT_BUDGET = {3: 120000, 4: 400000}
EPS_LADDER_RELATIVE = ('0', '1e-4', '3.16e-4', '1e-3', '3.16e-3', '1e-2',
                       '3.16e-2', '1e-1', '3.16e-1', '1')


def select(cells):
    index = common.bank_index()
    out = []
    for family, K, (low, high), count in cells:
        members = sorted((e for e in index if e['family'] == family
                          and e['K_generating'] == K
                          and low <= e['reference_set_size'] <= high),
                         key=lambda e: (e['reference_set_size'], e['name']))
        if not members:
            continue
        take = min(count, len(members))
        picks = [members[0]] if take == 1 else [
            members[round(j * (len(members) - 1) / (take - 1))] for j in range(take)]
        seen = set()
        for entry in picks:
            if entry['name'] in seen:
                continue
            seen.add(entry['name'])
            out.append(dict(entry, cell=f'{family}_K{K}_N{low}-{high}'))
    return out


def leafsets(case):
    """The full panel plus every registered stage-288 removal draw."""
    leaves = sorted(case['observed'])
    out = [dict(level='full', regime='none', replicate=0, removed=[],
                n_removed=0, n_retained=len(leaves), detail={})]
    for level in s288.LEVELS:
        for draw in s288.removal_draws(case, level):
            out.append(dict(draw))
    return out


# ------------------------------------------------------------------ prepare

def prepare():
    OUT.mkdir(parents=True, exist_ok=True)
    assert not PROTOCOL.exists(), 'Never overwrite frozen evidence'
    cross, probe = select(CROSS_CELLS), select(PROBE_CELLS)
    for row in cross + probe:
        case = load_case(row)
        row['n_edges'] = len(case['edges'])
        row['min_generating_p'] = str(min(s['p'] for s in case['generating_sources']))
        row['median_observed'] = str(sorted(case['observed'].values())[len(case['observed']) // 2])
        row['n_leafsets'] = len(leafsets(case))
        row['supports_k2'] = comb(row['n_edges'], 2)
    for row in probe:
        row['probe_k_list'] = [1, 2] + [k for k in (3, 4)
                                        if comb(row['n_edges'], k) <= PROBE_SUPPORT_BUDGET[k]]
    protocol = dict(
        stage=296, track='eps_tolerant_source_budget',
        registered_utc=datetime.now(timezone.utc).isoformat(),
        purpose=('Settle the stage-294 ambiguity by removing the heuristic from the loop and '
                 'making noise the controlled variable: the eps-tolerant minimum source count '
                 'is a property of the model, computed exactly.'),
        definition=('Given observations d_r and a tolerance eps >= 0, the minimum number k of '
                    'strict interior sources on distinct original edges, with weights p_i > 0 '
                    'summing to one and a common mean pendant y > 0, such that '
                    '|y + sum_i p_i dist_T(r, q_i) - d_r| <= eps for every observed r.'),
        method=dict(
            threshold_form=('For a fixed support S, theta(S) = min over the CLOSED set '
                            '{alpha, beta, y >= 0, sum = 1} of max_r |pred_r - d_r|, one exact '
                            'rational LP.  Theta_k = min over all C(|E|, k) supports of size k. '
                            'Theta_k is non-increasing in k and, for every eps > theta(S), the '
                            'STRICT system on S is feasible (segment to any strictly feasible '
                            'box point; the prediction is affine).  Hence the eps-tolerant '
                            'minimum is m(eps) = min{k : Theta_k < eps} for eps > 0, and the '
                            'whole eps axis is determined by the exact rational vector '
                            '(Theta_1, Theta_2, ...) - no eps ladder is run, the ladder is a '
                            'reporting choice applied afterwards.'),
            m_at_zero=('m(0) is decided separately because the closed LP admits degenerate '
                       'configurations (zero endpoint mass, zero pendant): the smallest k for '
                       'which some support with theta(S) = 0 is STRICTLY feasible, decided by a '
                       'max-common-slack LP whose optimum must be > 0.'),
            top_of_ladder=('Theta_K = 0 for every leaf set, certified per job by verifying the '
                           'generating witness on the retained observations with the frozen '
                           'independent verifier, so sweeps stop at k = K - 1.'),
            lp=('Dense two-phase tableau simplex over Fractions written for this stage '
                '(stage296_common.lp_minimize), Dantzig pricing with a switch to Bland after 60 '
                'iterations.  Row generation over the 2|R| tolerance rows: the active-subset LP '
                'is a relaxation, so its value is an exact lower bound, which both certifies '
                'optimality when no leaf is violated and prunes exactly against the incumbent.'),
            forbidden=('The frozen structural localisation - leaf envelope h, weighted deficit '
                       'bound, positive-anchor filters, solve_partial, stage292_all_minimisers, '
                       'stage295_reduced_minimum - is derived from exact equality and is unsound '
                       'under a tolerance.  It is used only in the eps = 0 validation phase.'),
            verification=('Every reported witness is re-evaluated on the ORIGINAL tree with the '
                          'frozen min(d(r,u) + xL, d(r,v) + (1-x)L) form - an implementation '
                          'independent of the interpolation the LP uses - and the recomputed '
                          'max |pred - d| must equal theta exactly.')),
        draw_plan=dict(
            imported_from='stage288_incompleteness.removal_draws',
            seed=s288.SEED, levels=list(s288.LEVELS), regimes=list(s288.REGIMES),
            replicates=dict(s288.REPLICATES),
            note='Identical draws to stage 288, so the eps = 0 column reproduces that stage.'),
        eps_ladder_relative_to_median_observed=list(EPS_LADDER_RELATIVE),
        within_eps_rule=('Every incompleteness comparison is made against that case own '
                         'full-panel answer AT THE SAME eps, never against the eps = 0 answer, '
                         'so the noise-induced collapse and the incompleteness-induced collapse '
                         'are not confounded.'),
        scale_ceiling=dict(
            cross_arm='k <= K - 1 <= 2 exhaustively; K in {2, 3} only.',
            probe_arm=('k <= 2 always, k = 3 when C(|E|, 3) <= %d, k = 4 when C(|E|, 4) <= %d.'
                       % (PROBE_SUPPORT_BUDGET[3], PROBE_SUPPORT_BUDGET[4])),
            reason=('Exhaustive enumeration costs C(|E|, k) exact LPs at roughly 6 ms each; the '
                    'frozen localisation that would avoid it is unsound under a tolerance.')),
        budgets=dict(workers=WORKERS, cross_wall_seconds=CROSS_WALL,
                     cross_address_space_bytes=CROSS_AS, probe_wall_seconds=PROBE_WALL,
                     probe_address_space_bytes=PROBE_AS,
                     validate_wall_seconds=VALIDATE_WALL,
                     validate_address_space_bytes=VALIDATE_AS,
                     calibration='Timing pilot stage296_pilot.py, measurement only.'),
        analysis=('Retain every job including timeouts and errors; never impute.  A comparison '
                  'uses only (case, draw) pairs where both the draw and its own full-panel '
                  'baseline completed.'),
        limitations=[
            'Exactly explainable families only: a tolerance has no meaning without an exact '
            'reference answer to compare against.',
            'The cross arm covers K in {2, 3}; K >= 5 is observed only through Theta_1, Theta_2 '
            '(and Theta_3, Theta_4 where the support budget allows) at the full panel.',
            'eps is an L-infinity tolerance; a median residual is not the eps a dataset needs.',
            'At the measure-zero coincidence eps == Theta_k > 0 the convention m(eps) = '
            'min{k : Theta_k < eps} is used, which is the conservative (larger-m) choice.',
            'No claim about which support is recovered, and no claim about any estimator.'],
        cross_cases=cross, probe_cases=probe,
        script_hashes={n: sha(ROOT / 'scripts' / n)
                       for n in common.dependency_closure('stage296_tolerance.py')},
        bank_protocol_sha256=sha(ROOT / 'results/stage287_instance_bank_protocol.json'),
        stage288_protocol_sha256=sha(ROOT / 'results/stage288_incompleteness_protocol.json'),
        machine=dict(platform=platform.platform(), python=sys.version))
    PROTOCOL.write_text(json.dumps(serial(protocol), indent=2) + '\n')
    print(json.dumps(dict(protocol=str(PROTOCOL), sha256=sha(PROTOCOL),
                          cross_cases=len(cross), probe_cases=len(probe),
                          cross_jobs=sum(r['n_leafsets'] for r in cross),
                          probe_jobs=len(probe)), indent=2))


# -------------------------------------------------------------------- worker

def theta_worker(spec):
    case, keep, klist = spec['case'], spec['leaves'], spec['k_list']
    edges, obs = case['edges'], case['observed']
    dvals = [obs[r] for r in keep]
    K = len(case['generating_sources'])
    Dv = s296.leaf_vectors(edges, keep)
    answer = dict(status='complete', n_retained=len(keep), K=K, thetas={}, witnesses={},
                  sweep_seconds={}, visited={}, pruned={})
    generating = dict(sources=[dict(edge=s['edge'], p=str(s['p']), x=str(s['x']))
                               for s in case['generating_sources']],
                      pendant=str(case['generating_pendant']))
    answer['generating_witness'] = verify_witness(edges, {r: obs[r] for r in keep},
                                                  case['cap'], generating)
    zero_at = None
    for k in klist:
        if k >= K:
            break
        t0 = time.perf_counter()
        out = s296.theta_min_over_supports(Dv, edges, dvals, k)
        elapsed = time.perf_counter() - t0
        answer['sweep_seconds'][k] = round(elapsed, 3)
        if out['status'] != 'complete':
            answer.update(status='sweep_' + out['status'], failed_k=k)
            return answer
        theta = out['theta']
        answer['thetas'][k] = theta
        answer['visited'][k] = out['visited']
        answer['pruned'][k] = out['pruned']
        dev = s296.independent_deviation(edges, out['support'], out['z'], k, obs, leaves=keep)
        answer['witnesses'][k] = dict(
            support=list(out['support']), theta=theta, independent_deviation=dev,
            independent_agrees=(dev == theta),
            sources=[{kk: vv for kk, vv in s.items()}
                     for s in s296.witness_sources(out['support'], out['z'], k)[0]],
            pendant=out['z'][2 * k])
        if theta == 0:
            zero_at = k
            break
    answer['theta_K_is_zero_certified'] = bool(answer['generating_witness'].get('verified'))
    m0, m0_detail = None, None
    if zero_at is not None:
        for k in range(zero_at, K):
            t0 = time.perf_counter()
            fz = s296.first_strict_zero(Dv, edges, dvals, k)
            answer['sweep_seconds'][f'strict{k}'] = round(time.perf_counter() - t0, 3)
            if fz['status'] == 'found':
                m0, m0_detail = k, dict(support=list(fz['support']), slack=fz['slack'],
                                        tested=fz['tested'])
                break
            if fz['status'] != 'none':
                answer.update(status='strict_' + fz['status'], failed_k=k)
                return answer
    swept = [k for k in klist if k < K]
    highest_swept = max(swept) if swept else 0
    if m0 is None:
        if highest_swept >= K - 1:
            m0 = K
            m0_detail = dict(by='generating_witness',
                             note='every k < K was swept and admits no strict exact explanation')
        else:
            m0_detail = dict(by='undetermined',
                             note='sweeps stopped below K - 1; m(0) is only bracketed')
    answer['m_at_zero'] = m0
    answer['m_at_zero_bounds'] = [highest_swept + 1 if m0 is None else m0, K]
    answer['m_at_zero_detail'] = m0_detail
    return answer


def check_witness(job, answer):
    """Parent-side independent check, outside child timing."""
    flags = []
    for k, w in (answer.get('witnesses') or {}).items():
        if not w.get('independent_agrees'):
            flags.append(f'k{k}_witness_mismatch')
    if answer.get('status') == 'complete' and not answer.get('theta_K_is_zero_certified'):
        flags.append('generating_witness_unverified')
    answer['parent_flags'] = flags


# ----------------------------------------------------------------- run cross

def run_arm(which):
    protocol = json.loads(PROTOCOL.read_text())
    cases = protocol['cross_cases'] if which == 'cross' else protocol['probe_cases']
    stream = CROSS_STREAM if which == 'cross' else PROBE_STREAM
    assert not stream.exists(), 'Never overwrite frozen evidence'
    wall, space = ((CROSS_WALL, CROSS_AS) if which == 'cross' else (PROBE_WALL, PROBE_AS))
    jobs = []
    for row in cases:
        case = load_case(row)
        sets = leafsets(case) if which == 'cross' else leafsets(case)[:1]
        klist = ([1, 2] if which == 'cross' else row['probe_k_list'])
        for ls in sets:
            keep = sorted(set(case['observed']) - set(ls['removed']))
            jobs.append(dict(
                meta=dict(arm=which, name=row['name'], family=row['family'],
                          cell=row['cell'], K=row['K_generating'],
                          N=row['reference_set_size'], n_edges=len(case['edges']),
                          level=ls['level'], regime=ls['regime'], replicate=ls['replicate'],
                          n_removed=ls['n_removed'], n_retained=len(keep),
                          removed=ls['removed'], detail=ls['detail'],
                          min_generating_p=row['min_generating_p'],
                          median_observed=row['median_observed'], k_list=klist),
                spec=dict(case=case, leaves=keep, k_list=klist)))
    jobs.sort(key=lambda j: -j['meta']['n_edges'])
    print(f'{which}: {len(jobs)} jobs over {len(cases)} cases', flush=True)
    t0 = time.perf_counter()
    run_pool(jobs, theta_worker, wall, space, WORKERS, stream, post=check_witness,
             progress=25)
    print(f'{which} done in {time.perf_counter() - t0:.1f} s', flush=True)


# ------------------------------------------------------------------ validate

def frozen_worker(spec):
    case, keep = spec['case'], spec['leaves']
    observed = {r: case['observed'][r] for r in keep}
    raw = solve_partial(case['edges'], observed, cap=case['cap'], seconds=None,
                        max_states=None, max_bases=None, max_faces=None)
    return dict(status=raw.get('status', 'complete'),
                minimum=raw.get('minimum_sources_within_cap'), reason=raw.get('reason'))


def validate():
    protocol = json.loads(PROTOCOL.read_text())
    assert not VALIDATE_STREAM.exists(), 'Never overwrite frozen evidence'
    jobs = []
    for row in protocol['cross_cases']:
        case = load_case(row)
        for ls in leafsets(case):
            keep = sorted(set(case['observed']) - set(ls['removed']))
            jobs.append(dict(
                meta=dict(name=row['name'], K=row['K_generating'], level=ls['level'],
                          regime=ls['regime'], replicate=ls['replicate'],
                          n_retained=len(keep), removed=ls['removed']),
                spec=dict(case=case, leaves=keep)))
    print(f'validate: {len(jobs)} frozen runs', flush=True)
    run_pool(jobs, frozen_worker, VALIDATE_WALL, VALIDATE_AS, WORKERS, VALIDATE_STREAM,
             progress=100)


# ----------------------------------------------------------------- calibrate

def calibrate():
    """Real estimation error, in the same units as eps, on the panel geometry."""
    assert not CALIBRATION.exists(), 'Never overwrite frozen evidence'
    protocol = json.loads(PROTOCOL.read_text())
    index = {e['name']: e for e in common.bank_index()}
    rows = []
    for row in protocol['cross_cases'] + protocol['probe_cases']:
        case = load_case(row)
        exact = case['observed']
        median_observed = sorted(exact.values())[len(exact) // 2]
        stem = row['name'][:-len('_exact')]
        partners = []
        if row['family'] == 'exact_10k':
            partners = [(stem + '_krepp_estimated', 'krepp_ani')]
        else:
            base = stem[:-len('_noise0')]
            partners = [(f'{base}_noise0_published', 'published_noise0')]
            partners += [(f'{base}_noise{n}_l{l}_published', f'published_noise{n}_l{l}')
                         for n in (1, 2) for l in (100, 1000, 10000)]
        for partner_name, label in partners:
            entry = index.get(partner_name)
            if entry is None:
                continue
            raw = json.loads((ROOT / entry['path']).read_text())
            assert sha(ROOT / entry['path']) == entry['sha256']
            noisy = {int(r): F(d) for r, d in raw['observed'].items()}
            shared = sorted(set(noisy) & set(exact))
            if not shared:
                continue
            resid = sorted(abs(noisy[r] - exact[r]) for r in shared)
            rows.append(dict(
                case=row['name'], family=row['family'], K=row['K_generating'],
                N=row['reference_set_size'], partner=partner_name, label=label,
                n_shared=len(shared), n_exact_leaves=len(exact),
                median_observed=str(median_observed),
                residual_median=str(resid[len(resid) // 2]),
                residual_max=str(resid[-1]),
                residual_median_rel=float(resid[len(resid) // 2] / median_observed),
                residual_max_rel=float(resid[-1] / median_observed)))
    CALIBRATION.write_text(json.dumps(serial(dict(
        purpose=('The smallest tolerance at which the TRUE configuration still explains the '
                 'real data is exactly max_r |d_estimated(r) - d_exact(r)| on the shared leaf '
                 'set: the L-infinity residual.  It is reported with the median residual and '
                 'both relative to the median observed distance, so it can be placed on the '
                 'same axis as Theta_k.'),
        caveats=['krepp distances are real ANI estimates on the same tree with the same '
                 'generating sources; the exact partner is the model-implied distance vector.',
                 'published_biotrees is not a clean graded noise ladder: source weights are '
                 'unpublished and the bank assumes equal weights, so the noise0 residual '
                 'already exceeds the sequence-noise increment (stage 291).',
                 'krepp observed leaf sets can be a strict subset of the exact ones; residuals '
                 'use the shared leaves only and n_shared records the overlap.'],
        rows=rows)), indent=2) + '\n')
    print(json.dumps(dict(path=str(CALIBRATION), rows=len(rows)), indent=2))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('phase', choices=['prepare', 'cross', 'probe', 'validate', 'calibrate'])
    args = ap.parse_args()
    if args.phase == 'prepare':
        prepare()
    elif args.phase in ('cross', 'probe'):
        run_arm(args.phase)
    elif args.phase == 'validate':
        validate()
    else:
        calibrate()
