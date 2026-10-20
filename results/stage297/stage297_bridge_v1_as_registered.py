#!/usr/bin/env python3
"""Stage297 bridge: the unknown-support erasure threshold at K = 10 (and K = 12).

Stage292 measured t_unknown over the exact bank at K in {2,3,5,7} and stopped
there: the FROZEN minimum pass its all-minimisers routine calls returned
'unknown' after 3600 s on the smallest K = 10 case
(results/stage292_k10_ceiling.txt), and all eight K = 10 jobs of its stride-3
batch hit the per-case wall (results/stage292_bridge_K10_12.jsonl).  This script
repeats that measurement with the Stage297 combined routine, which pays the
reduced enumeration once for the minimum and the complete minimiser set
together.

Everything else is Stage292's, imported rather than restated:
  * the threshold search `stage292_bridge.find_threshold` -- exponential
    bracketing then bisection, both bracket endpoints evaluated;
  * the rank scanner `stage292_bridge.RankScanner` for t_known;
  * the removal orders `stage292_orders.order_for` with REGIMES
    (random 11, random 23, whole_clade 11, whole_clade 23,
     adversarial_nearest) and MIN_REMAINING = 3;
  * the uniqueness predicate: minimum == K, exactly one minimiser, and that
    minimiser is the true support.
A probe whose status is not decisive aborts that regime with
threshold_status='unknown_probe:<status>'; it is never read as a failure.
Every probe is retained with its counters so a resource ceiling can be
reported with numbers instead of being quietly absorbed into the claim.
"""
import json
import os
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
sys.dont_write_bytecode = True

import stage297_common as common                                          # noqa: E402
from stage297_common import OUTDIR, run_pool                              # noqa: E402
from stage290_erasure_core import support_of, build_struct, erasure_analysis  # noqa: E402
from stage292_bridge import RankScanner, find_threshold, ProbeUnknown, REGIMES, MIN_REMAINING  # noqa: E402
from stage292_orders import order_for                                     # noqa: E402
from stage292_all_minimisers import VertexDistances                       # noqa: E402
from stage297_reduced_all_minimisers import all_minimum_supports_reduced as reduced  # noqa: E402

WORKERS = int(os.environ.get('STAGE297_WORKERS', '8'))
MEMORY = 10 << 30


def _probe_counters(res):
    st = res['statistics']
    return dict(status=res['status'], minimum=res['minimum'],
                n_supports=res['n_supports'], reason=res.get('reason'),
                deficit_sets=st['deficit_sets'], faces_to_lp=st['faces_to_lp'],
                masks_enumerated=st['masks_enumerated'],
                masks_full_enumeration=st['masks_full_enumeration'],
                masks_skipped_forced=st['masks_skipped_forced'],
                feasible_patterns=st['feasible_patterns'],
                configurations=st['configurations'],
                oracle_calls=st['oracle_calls'], lift_products=st['lift_products'],
                supports_visited=st['supports_visited'],
                free_bits_max=st['free_bits_max'], A_bits_max=st['A_bits_max'],
                timing=res.get('timing'), budget_counts=res.get('budget_counts'))


def job(spec):
    entry = spec['entry']
    case = common.load_verified(entry)
    support = support_of(case)
    K = len(support)
    observed = case['observed']
    edges = case['edges']
    leaves = sorted(observed)
    st = build_struct(edges, support, leaves)
    dists = VertexDistances(edges)
    seconds = spec['probe_seconds']
    screen = spec.get('screen_rows', 48 if len(leaves) > 2000 else 0)
    out = dict(kind='bridge', stage=297, name=entry['name'], family=entry['family'],
               path=entry['path'], sha256=entry['sha256'], K=K,
               n_edges=len(edges), R_size=len(observed),
               original_leaves=case['original_leaves'], true_support=support,
               support_is_matching=st['matching'],
               n_pendant_support=len(st['pendant_support']),
               probe_seconds=seconds, screen_rows=screen,
               mode=spec.get('mode', 'combined'), regimes=[])

    if spec.get('want_known', True):
        t0 = time.perf_counter()
        known = erasure_analysis(st, want_witness=False)
        scanner = RankScanner(st)
        out['cogirth'] = known['delta']
        out['cogirth_tolerance'] = known['tolerance']
        out['s'] = known['s']
        out['rank_full'] = scanner.rank([])
        out['rank_identity'] = bool(out['rank_full'] == known['s'] + 1)
        out['known_seconds'] = time.perf_counter() - t0
    else:
        known, scanner = None, None

    t0 = time.perf_counter()
    base = reduced(edges, observed, cap=K, seconds=seconds, mode=out['mode'],
                   screen_rows=screen, distances=dists)
    out['baseline'] = _probe_counters(base)
    out['baseline_status'] = base['status']
    out['baseline_minimum'] = base['minimum']
    out['baseline_n_supports'] = base['n_supports']
    out['baseline_minimum_equals_K'] = bool(base['minimum'] == K)
    out['baseline_unique_true'] = bool(
        base['status'] == 'complete' and base['minimum'] == K and
        base['n_supports'] == 1 and list(base['supports'][0]) == support)
    out['baseline_seconds'] = time.perf_counter() - t0
    out['baseline_supports'] = base['supports'][:50]

    if spec.get('cross_check_modes'):
        t0 = time.perf_counter()
        other = reduced(edges, observed, cap=K, seconds=seconds,
                        mode='reduced_two_pass', screen_rows=screen, distances=dists)
        out['cross_check'] = dict(
            mode='reduced_two_pass', status=other['status'],
            minimum=other['minimum'], n_supports=other['n_supports'],
            reduced_minimum=other.get('reduced_minimum'),
            agree_minimum=bool(other['minimum'] == base['minimum']),
            agree_supports=bool(sorted(map(tuple, other['supports'])) ==
                                sorted(map(tuple, base['supports']))),
            seconds=time.perf_counter() - t0)

    if not spec.get('want_thresholds', True):
        return [out]

    emax = len(observed) - MIN_REMAINING
    for regime, seed in REGIMES:
        if spec.get('regimes') and regime not in spec['regimes']:
            continue
        order = order_for(edges, st, regime, seed,
                          sources=case['generating_sources'], dists=dists)
        t_known = None
        if scanner is not None:
            for e in range(1, emax + 1):
                if scanner.rank(order[:e]) != known['s'] + 1:
                    t_known = e - 1
                    break
        probes = {}

        def ok(e):
            if e in probes:
                if probes[e]['ok'] is None:
                    raise ProbeUnknown(probes[e]['status'])
                return probes[e]['ok']
            removed = set(order[:e])
            sub = {r: observed[r] for r in leaves if r not in removed}
            t1 = time.perf_counter()
            res = reduced(edges, sub, cap=K, seconds=seconds, mode=out['mode'],
                          stop_after=2, screen_rows=screen, distances=dists)
            if res['status'] == 'complete':
                good = bool(res['minimum'] == K and res['n_supports'] == 1 and
                            list(res['supports'][0]) == support)
            elif res['status'] == 'early_exit':
                good = False
            else:
                good = None
            probes[e] = dict(e=e, ok=good, seconds=time.perf_counter() - t1,
                             **_probe_counters(res))
            if good is None:
                raise ProbeUnknown(res['status'])
            return good

        t2 = time.perf_counter()
        try:
            t_unknown, status = find_threshold(ok, emax, probes)
        except ProbeUnknown as exc:
            t_unknown, status = None, 'unknown_probe:' + str(exc)
        out['regimes'].append(dict(
            regime=regime, seed=seed, t_known_sequence=t_known,
            t_unknown_sequence=t_unknown, emax=emax, threshold_status=status,
            n_probes=len(probes), seconds=time.perf_counter() - t2,
            fraction=(None if t_unknown is None else t_unknown / len(observed)),
            probe_table=sorted(probes.values(), key=lambda z: z['e'])))
    return [out]


def build_jobs(tag, Ks, per_size, probe_seconds, wall, want_thresholds=True,
               cross_check=False, families=('exact',), limit=None, want_known=True):
    entries = common.bank_cases(families=families, Ks=Ks)
    chosen = common.stratified(entries, per_size) if per_size else entries
    chosen.sort(key=lambda e: (e['K_generating'], e['reference_set_size'], e['name']))
    if limit:
        chosen = chosen[:limit]
    return [dict(tag=tag, entry=e, probe_seconds=probe_seconds, wall=wall,
                 memory=MEMORY, want_thresholds=want_thresholds,
                 want_known=want_known, cross_check_modes=cross_check)
            for e in chosen]


PLANS = {
    # name: (Ks, cases per distinct |R|, probe seconds, per-case wall, thresholds?)
    'pilot_k10_baseline': dict(Ks=(10,), per_size=1, probe_seconds=1200, wall=2400,
                               want_thresholds=False, cross_check=True),
    'pilot_k10_threshold': dict(Ks=(10,), per_size=1, probe_seconds=1200, wall=7200,
                                want_thresholds=True, cross_check=False, limit=6),
    'k12_baseline': dict(Ks=(12,), per_size=None, probe_seconds=3600, wall=4200,
                         want_thresholds=False, cross_check=False, want_known=False),
}


def main():
    name = sys.argv[1]
    plan = dict(PLANS[name])
    if len(sys.argv) > 2:
        plan['per_size'] = int(sys.argv[2])
    if len(sys.argv) > 3:
        plan['limit'] = int(sys.argv[3])
    jobs = build_jobs(name, **plan)
    run_pool(jobs, job, WORKERS, plan['wall'], MEMORY,
             OUTDIR / ('stage297_bridge_%s.jsonl' % name), label=name, every=1)


if __name__ == '__main__':
    main()
