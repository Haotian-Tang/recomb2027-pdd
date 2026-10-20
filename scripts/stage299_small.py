#!/usr/bin/env python3
"""Stage 299 step 2: Theta_1 under reference removal on the three small arms.

exact_partner / krepp_ani / krepp_placement, every registered leaf set.  Each
leaf set is decided TWICE and verified once:

  * frozen per-support LP (stage296_common.theta_support) on every one of the
    |E| edges, no cutoff - this is the number that is reported;
  * stage 298's k = 1 closed form on the same leaf set, compared edge by edge
    by exact rational equality.  Removal is what makes this comparison worth
    repeating: a thinned panel leaves many edges with retained references on
    one side only, a regime stage 298 never exercised, and the travellers' arm
    depends on the closed form being right there;
  * the attaining witness re-evaluated on the original tree with the frozen
    min(d(r,u) + xL, d(r,v) + (1-x)L) verifier.

Every run is retained, including failures.
"""
import json
import os
import sys
import time
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import stage288_common as common               # noqa: E402
import stage296_common as s296                 # noqa: E402
import stage298_common as s298                 # noqa: E402
import stage299_common as s299                 # noqa: E402

OUT = s299.OUT
PROTOCOL = OUT / 'stage299_protocol.json'


def worker(spec):
    edges, leaves, dvals, obs = spec['edges'], spec['leaves'], spec['dvals'], spec['obs']
    t0 = time.perf_counter()
    fz = s299.theta1_frozen_all_edges(edges, leaves, dvals)
    t_frozen = time.perf_counter() - t0
    if fz['status'] != 'ok':
        return dict(status=fz['status'], edge=fz.get('edge'), frozen_seconds=t_frozen)
    t0 = time.perf_counter()
    cf = s298.theta1_all_edges(edges, leaves, dvals)
    t_closed = time.perf_counter() - t0
    mism = [i for i in range(len(edges)) if cf['theta'][i] != fz['theta'][i]]
    theta1 = fz['theta1']
    att = s299.attaining(fz['theta'], theta1)
    srcs, y = s296.witness_sources((fz['best'],), fz['z'], 1)
    dev = s296.independent_deviation(edges, (fz['best'],), fz['z'], 1, obs, leaves=leaves)
    deg = s299.degrees(edges)
    med = s299.median_observed(dvals)
    return dict(
        status='complete', theta1=theta1, argmin=fz['best'],
        attaining=att, n_attaining=len(att),
        attaining_terminal=[i for i in att if s299.is_terminal(deg, edges[i])],
        argmin_is_terminal=s299.is_terminal(deg, edges[fz['best']]),
        argmin_x=srcs[0]['x'], argmin_p=srcs[0]['p'], pendant=y,
        closed_form_theta1=cf['theta1'], closed_form_argmin=cf['best'],
        closed_form_x=cf['x'][fz['best']],
        closed_form_x_agrees=(cf['x'][fz['best']] == srcs[0]['x']),
        closed_form_mismatches=len(mism), mismatch_examples=mism[:5],
        closed_form_theta1_agrees=(cf['theta1'] == theta1),
        independent_deviation=dev, witness_agrees=(dev is not None and dev == theta1),
        median_observed=med, theta1_over_median_observed=float(theta1 / med) if med else None,
        n_leaves=len(leaves), n_edges=len(edges),
        frozen_seconds=t_frozen, closed_seconds=t_closed,
        seconds_per_support=t_frozen / max(1, len(edges)))


def main():
    protocol = json.loads(PROTOCOL.read_text())
    stream = OUT / 'stage299_small_runs.jsonl'
    assert not stream.exists(), 'Never overwrite frozen evidence'
    arms = set(s299.SMALL_ARMS)
    cases = {}
    for row in protocol['cases']:
        if row['arm'] not in arms:
            continue
        cases[(row['arm'], row['name'])] = s299.prepare_case(
            s299.load_json(row['path'], row['sha256']), row['name'])
    print('loaded %d cases' % len(cases), flush=True)
    jobs = []
    for row in protocol['leafsets']:
        if row['arm'] not in arms:
            continue
        case = cases[(row['arm'], row['name'])]
        removed = set(int(r) for r in row['removed'])
        # the exact_partner arm is evaluated on the KREPP panel, leaf for leaf, so
        # the 0-3 references that exist only in the exact panel are never observed
        # here; that is what makes the exact and estimated arms paired.
        outside = (set(int(r) for r in (row['exact_only_references'] or []))
                   if row['arm'] == 'exact_partner' else set())
        leaves = [r for r in sorted(case['observed'])
                  if r not in removed and r not in outside]
        assert len(leaves) == row['n_retained'], (row['name'], row['level'], row['regime'],
                                                  len(leaves), row['n_retained'])
        dvals = [case['observed'][r] for r in leaves]
        obs = {r: d for r, d in zip(leaves, dvals)}
        meta = dict(arm=row['arm'], name=row['name'], stem=row['stem'], family=row['family'],
                    level=row['level'], regime=row['regime'], replicate=row['replicate'],
                    n_removed=row['n_removed'], n_retained=row['n_retained'],
                    removed_fraction=(row['n_removed'] / float(row['n_removed'] + row['n_retained'])),
                    K_generating=row['K_generating'], n_edges=row['n_edges'],
                    reference_set_size=row['reference_set_size'],
                    exact_partner=row['exact_partner'], detail=row['detail'],
                    removed=row['removed'])
        jobs.append(dict(meta=meta, spec=dict(edges=case['edges'], leaves=leaves,
                                              dvals=dvals, obs=obs)))
    jobs.sort(key=lambda j: (j['meta']['level'] != 'full', j['meta']['arm'], j['meta']['name']))
    print('%d jobs (%d full panels)' % (len(jobs), sum(1 for j in jobs
                                                       if j['meta']['level'] == 'full')), flush=True)
    t0 = time.time()
    res = common.run_pool(jobs, worker, protocol['budgets']['small_wall_seconds'],
                          protocol['budgets']['address_space_bytes'],
                          protocol['budgets']['workers'], stream, progress=1000)
    el = time.time() - t0
    ok = sum(1 for r in res if r['result'].get('status') == 'complete')
    mism = sum(r['result'].get('closed_form_mismatches', 0) for r in res)
    wit = sum(1 for r in res if r['result'].get('witness_agrees'))
    print('done in %.0f s: %d/%d complete, %d closed-form mismatches, %d witnesses agree'
          % (el, ok, len(res), mism, wit))
    summary = dict(stage=299, step='small_arms', n_jobs=len(jobs), complete=ok,
                   closed_form_mismatches=mism, witness_agree=wit,
                   elapsed_seconds=el,
                   failures=[dict(r, result={k: v for k, v in r['result'].items()
                                             if k in ('status', 'reason', 'exception', 'edge')})
                             for r in res if r['result'].get('status') != 'complete'][:200])
    (OUT / 'stage299_small_summary.json').write_text(json.dumps(common.serial(summary), indent=1))
    print('wrote', OUT / 'stage299_small_summary.json')


if __name__ == '__main__':
    main()

