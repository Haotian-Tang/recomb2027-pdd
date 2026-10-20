#!/usr/bin/env python3
"""Stage297 validation of the combined reduced all-minimisers routine.

Three tracks, no track is allowed to patch a disagreement:

  A  REPLAY.  Every probe point recorded in results/stage292_validate_replay.jsonl
     (130 cases x 5 regimes x the probes of their bisections, plus e = 0; the
     removal order is regenerated with the same frozen `removal_order`) is
     recomputed with the Stage292 routine AND with all three Stage297 modes, and
     the complete minimiser SETS are compared element for element -- the stored
     records keep only counts, so the sets are recomputed on both sides here.
  B  FRESH.  The 500 instances of results/stage292_validate_fresh.jsonl are
     regenerated from `stage292_validate.random_instance` (same seeds, same
     generator) and compared against the independent exhaustive C(|E|,k)
     enumeration `stage292_validate.brute_all`, against the Stage292 routine and
     between reduced and full (2^|A|) face enumeration.
  C  HIGH K.  The full-panel baseline of all 676 Stage292 bridge cases
     (K in {2,3,5,7}) is recomputed with the Stage292 routine and with Stage297
     'combined', and both are compared with the value Stage292 recorded.  This
     is the only track that reaches K = 5 and K = 7, the regime the K = 10 run
     extends.  On a deterministic subset it also checks that with the reduction
     disabled the Stage297 face count equals the face count of the Stage292
     enumeration run with full_mask=True -- the correctness signal Stage295 used.
"""
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
sys.dont_write_bytecode = True

import stage297_common as common                                          # noqa: E402
from stage297_common import OUTDIR, encode, run_pool                      # noqa: E402
from stage290_erasure_core import (load_bank_case, support_of, CaseTables,  # noqa: E402
                                   removal_order, build_struct)
from stage292_all_minimisers import all_minimum_supports, backbone_patterns  # noqa: E402
from stage292_validate import random_instance, brute_all                  # noqa: E402
from stage297_reduced_all_minimisers import all_minimum_supports_reduced as reduced  # noqa: E402
from stage27_general_fpt import Budget                                    # noqa: E402
from stage26_missing_leaf import observed_backbone                        # noqa: E402
from stage297_reduced_all_minimisers import minimum_patterns              # noqa: E402

WORKERS = int(__import__('os').environ.get('STAGE297_WORKERS', '16'))
MEMORY = 6 << 30
WALL_A, WALL_B, WALL_C = 3600.0, 900.0, 5400.0


def _sets(res):
    return sorted(map(tuple, res['supports']))


def _arms(edges, observed, cap, want_full_mask=False):
    """The Stage292 routine and every Stage297 mode on one instance."""
    arms = {}
    t0 = time.perf_counter()
    arms['s292'] = all_minimum_supports(edges, observed, cap=cap)
    arms['s292_seconds'] = time.perf_counter() - t0
    for mode in ('combined', 'reduced_two_pass', 'frozen'):
        t0 = time.perf_counter()
        arms[mode] = reduced(edges, observed, cap=cap, mode=mode)
        arms[mode + '_seconds'] = time.perf_counter() - t0
    if want_full_mask:
        t0 = time.perf_counter()
        arms['combined_full_mask'] = reduced(edges, observed, cap=cap,
                                             mode='combined', full_mask=True)
        arms['combined_full_mask_seconds'] = time.perf_counter() - t0
    return arms


def _compare(arms, extra=None):
    ref = arms['s292']
    row = dict(extra or {})
    row['s292_status'] = ref['status']
    row['s292_minimum'] = ref['minimum']
    row['s292_n'] = ref['n_supports']
    row['s292_seconds'] = arms['s292_seconds']
    row['s292_frozen_minimum'] = ref.get('frozen_minimum')
    ref_set = _sets(ref)
    for mode in ('combined', 'reduced_two_pass', 'frozen', 'combined_full_mask'):
        if mode not in arms:
            continue
        res = arms[mode]
        row[mode + '_status'] = res['status']
        row[mode + '_minimum'] = res['minimum']
        row[mode + '_n'] = res['n_supports']
        row[mode + '_seconds'] = arms[mode + '_seconds']
        row['agree_minimum_' + mode] = bool(res['minimum'] == ref['minimum'])
        row['agree_supports_' + mode] = bool(_sets(res) == ref_set)
        row['agree_status_' + mode] = bool(res['status'] == ref['status'])
        row[mode + '_masks_enumerated'] = res['statistics']['masks_enumerated']
        row[mode + '_masks_full'] = res['statistics']['masks_full_enumeration']
        row[mode + '_faces_to_lp'] = res['statistics']['faces_to_lp']
        row[mode + '_oracle_calls'] = res['statistics']['oracle_calls']
    row['first_pass_minimum_reduced'] = arms['reduced_two_pass'].get('reduced_minimum')
    row['agree_minimum_pass'] = bool(
        arms['reduced_two_pass'].get('reduced_minimum') == ref.get('frozen_minimum'))
    if not all(v for k, v in row.items() if k.startswith('agree_')):
        row['disagreement'] = {k: encode(_sets(v)) for k, v in arms.items()
                               if isinstance(v, dict)}
        row['disagreement']['statuses'] = {k: v['status'] for k, v in arms.items()
                                           if isinstance(v, dict)}
    return row


# ------------------------------------------------------------------ track A

def replay_point(spec):
    rec = spec['record']
    entry = spec['entry']
    case = common.load_verified(entry)
    support = support_of(case)
    K = len(support)
    observed = case['observed']
    tables = CaseTables(case['edges'])
    st = build_struct(case['edges'], support, sorted(observed))
    orders = {}
    out = []
    for pr in spec['probes']:
        if pr['e'] == 0:
            keys = sorted(observed)
        else:
            key = (pr['regime'], pr['seed'])
            if key not in orders:
                orders[key] = removal_order(case, st, pr['regime'], pr['seed'],
                                            tables=tables)
            removed = set(orders[key][:pr['e']])
            keys = [r for r in sorted(observed) if r not in removed]
        sub = {r: observed[r] for r in keys}
        arms = _arms(case['edges'], sub, K)
        row = _compare(arms, dict(
            track='A', name=rec['name'], family=rec['family'], K=K,
            n_edges=len(case['edges']), R_size=len(observed),
            R_remaining=len(keys), regime=pr['regime'], seed=pr['seed'],
            e=pr['e'], recorded_minimum=pr.get('recorded_minimum'),
            recorded_n_supports=pr.get('recorded_n_supports')))
        if pr.get('recorded_minimum') is not None:
            row['agree_recorded_s292'] = bool(
                arms['s292']['minimum'] == pr['recorded_minimum'] and
                arms['s292']['n_supports'] == pr['recorded_n_supports'])
            row['agree_recorded_combined'] = bool(
                arms['combined']['minimum'] == pr['recorded_minimum'] and
                arms['combined']['n_supports'] == pr['recorded_n_supports'])
        out.append(row)
    return out


# ------------------------------------------------------------------ track B

def fresh_point(spec):
    inst = random_instance(spec['seed'])
    cap = inst['K']
    arms = _arms(inst['edges'], inst['observed'], cap, want_full_mask=True)
    t0 = time.perf_counter()
    bf = brute_all(inst['edges'], inst['observed'], cap)
    brute_seconds = time.perf_counter() - t0
    bf_set = sorted(map(tuple, bf['supports']))
    row = _compare(arms, dict(
        track='B', seed=spec['seed'], K=inst['K'], nleaf=inst['nleaf'],
        n_edges=len(inst['edges']), R_size=len(inst['R']),
        hidden=inst['nleaf'] - len(inst['R']), true_support=inst['support'],
        brute_minimum=bf['minimum'], brute_n=len(bf['supports']),
        brute_visited=bf['visited'], brute_seconds=brute_seconds))
    for mode in ('s292', 'combined', 'reduced_two_pass', 'frozen', 'combined_full_mask'):
        row['agree_brute_' + mode] = bool(_sets(arms[mode]) == bf_set and
                                          arms[mode]['minimum'] == bf['minimum'])
    if not all(v for k, v in row.items() if k.startswith('agree_')):
        row.setdefault('disagreement', {})['brute'] = encode(bf_set)
        row['disagreement']['instance'] = encode(inst)
    return [row]


# ------------------------------------------------------------------ track C

def face_count_check(edges, observed, cap, kstar):
    """Stage292's enumeration with full_mask=True vs Stage297's, face for face."""
    E = [(u, v, l) for u, v, l in edges]
    B = observed_backbone(E, {int(r): d for r, d in observed.items()})
    marked = {x['vertex'] for x in B['saturation_violations']}
    obs = {int(r): d for r, d in observed.items()}
    b292 = Budget(None, None, None, None)
    s292 = {}
    pats292, reason292 = backbone_patterns(B['edges'], obs, marked, cap, b292,
                                           accept=lambda c: c == kstar,
                                           full_mask=True, stats=s292)
    b297 = Budget(None, None, None, None)
    from stage297_reduced_all_minimisers import _new_stats
    s297 = _new_stats()
    pats297, minimum297, reason297 = minimum_patterns(B['edges'], obs, marked, cap,
                                                      b297, s297, full_mask=True,
                                                      kstar=kstar)
    b297r = Budget(None, None, None, None)
    s297r = _new_stats()
    pats297r, _, _ = minimum_patterns(B['edges'], obs, marked, cap, b297r, s297r,
                                      full_mask=False, kstar=kstar)
    return dict(
        faces_292_full=b292.counts.get('mass_faces'),
        faces_297_full=b297.counts.get('mass_faces'),
        faces_297_reduced=b297r.counts.get('mass_faces'),
        masks_297_full=s297['masks_enumerated'],
        masks_297_reduced=s297r['masks_enumerated'],
        skipped_forced=s297r['masks_skipped_forced'],
        deficit_sets_292=s292.get('deficit_sets'), deficit_sets_297=s297['deficit_sets'],
        patterns_292=len(pats292), patterns_297=len(pats297),
        patterns_297_reduced=len(pats297r),
        reason_292=reason292, reason_297=reason297,
        agree_face_count=bool(b292.counts.get('mass_faces') == b297.counts.get('mass_faces')),
        agree_pattern_count=bool(len(pats292) == len(pats297)),
        agree_reduced_patterns=bool(sorted(p['P'] for p in pats297) ==
                                    sorted(p['P'] for p in pats297r)))


def baseline_point(spec):
    entry = spec['entry']
    case = common.load_verified(entry)
    support = support_of(case)
    K = len(support)
    arms = _arms(case['edges'], case['observed'], K)
    row = _compare(arms, dict(
        track='C', name=entry['name'], family=entry['family'], K=K,
        n_edges=len(case['edges']), R_size=entry['reference_set_size'],
        true_support=support,
        recorded_minimum=spec['recorded_minimum'],
        recorded_n_supports=spec['recorded_n_supports'],
        recorded_unique_true=spec['recorded_unique_true']))
    row['agree_recorded_s292'] = bool(
        arms['s292']['minimum'] == spec['recorded_minimum'] and
        arms['s292']['n_supports'] == spec['recorded_n_supports'])
    row['agree_recorded_combined'] = bool(
        arms['combined']['minimum'] == spec['recorded_minimum'] and
        arms['combined']['n_supports'] == spec['recorded_n_supports'])
    if spec.get('face_check') and arms['s292']['minimum'] is not None:
        row['face_check'] = face_count_check(case['edges'], case['observed'], K,
                                             arms['s292']['minimum'])
        for key in ('agree_face_count', 'agree_pattern_count', 'agree_reduced_patterns'):
            row['agree_' + key] = bool(row['face_check'][key])
    if not all(v for k, v in row.items() if k.startswith('agree_')):
        row.setdefault('disagreement', {})['note'] = 'see arm fields'
    return [row]


# ---------------------------------------------------------------------------

def plan_A():
    index = {e['name']: e for e in common.bank_index()}
    rows = [json.loads(l) for l in (ROOT / 'results/stage292_validate_replay.jsonl').open()]
    by_case = {}
    for r in rows:
        if 'row_status' in r:
            continue
        by_case.setdefault(r['name'], []).append(
            dict(regime=r['regime'], seed=r['seed'], e=r['e'],
                 recorded_minimum=r.get('recorded_minimum'),
                 recorded_n_supports=r.get('recorded_n_supports')))
    jobs = []
    for name, probes in sorted(by_case.items()):
        jobs.append(dict(track='A', record=dict(name=name, family=index[name]['family']),
                         entry=index[name], probes=probes))
    jobs.sort(key=lambda j: -len(j['probes']))
    return jobs


def plan_B():
    return [dict(track='B', seed=s) for s in range(1, 501)]


def plan_C(face_every=10):
    index = {e['name']: e for e in common.bank_index()}
    rows = []
    for f in ('stage292_bridge_K2_3_5_7.jsonl', 'stage292_bridge_K7_resume.jsonl'):
        rows += [json.loads(l) for l in (ROOT / 'results' / f).open()]
    jobs = []
    for i, r in enumerate(sorted([r for r in rows if r.get('kind') == 'bridge'],
                                 key=lambda z: (z['K'], z['R_size'], z['name']))):
        jobs.append(dict(track='C', entry=index[r['name']],
                         recorded_minimum=r['baseline_minimum'],
                         recorded_n_supports=r['baseline_n_supports'],
                         recorded_unique_true=r['baseline_unique_true'],
                         face_check=bool(i % face_every == 0)))
    jobs.sort(key=lambda j: -j['entry']['reference_set_size'])
    return jobs


def main():
    tracks = sys.argv[1] if len(sys.argv) > 1 else 'ABC'
    if 'A' in tracks:
        run_pool(plan_A(), replay_point, WORKERS, WALL_A, MEMORY,
                 OUTDIR / 'stage297_validate_replay.jsonl', label='A', every=10)
    if 'B' in tracks:
        run_pool(plan_B(), fresh_point, WORKERS, WALL_B, MEMORY,
                 OUTDIR / 'stage297_validate_fresh.jsonl', label='B', every=50)
    if 'C' in tracks:
        run_pool(plan_C(), baseline_point, WORKERS, WALL_C, MEMORY,
                 OUTDIR / 'stage297_validate_baseline.jsonl', label='C', every=20)


if __name__ == '__main__':
    main()
