#!/usr/bin/env python3
"""Stage 299 step 3: Theta_1 under reference removal on the 29 real stool samples.

|E| = 21,147 and |R| = 10,575.  The frozen per-support LP costs 0.512 core-
seconds per edge here (stage 298, measured), i.e. 3.0 core-hours for ONE leaf
set and 2,270 core-hours for the 754 registered leaf sets, so the sweep uses
stage 298's k = 1 closed form, which that stage validated against the frozen LP
on 653,807 edges of these same trees with 0 disagreements.  Stage 299 adds the
part stage 298 could not cover - thinned panels, where many edges hold retained
references on one side only - two ways:

  * step 1 (stage299_validate) compares the two routes edge by edge on hard
    thinnings of a real case and on a synthetic star with an empty side;
  * step 2 (stage299_small) compares them edge by edge on EVERY leaf set of the
    three small arms, which is 2.3 million double-decided edges including every
    removal regime and level;
  * and here, per leaf set, the attaining edge plus a seeded random sample of
    `TD_SAMPLE_EDGES` edges are re-decided with the frozen LP, and the attaining
    witness is re-verified with the frozen min(...) form.

Every run is retained, including failures.
"""
import hashlib
import json
import random
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
    cf = s298.theta1_all_edges(edges, leaves, dvals)
    t_closed = time.perf_counter() - t0
    theta1 = cf['theta1']
    if theta1 is None:
        return dict(status='closed_form_empty', closed_seconds=t_closed)
    att = s299.attaining(cf['theta'], theta1)
    best = cf['best']
    deg = s299.degrees(edges)
    # frozen LP spot checks: the attaining edge plus a seeded random sample
    rng = random.Random(spec['sample_seed'])
    sample = sorted({best} | set(rng.randrange(len(edges))
                                 for _ in range(s299.TD_SAMPLE_EDGES)))
    t0 = time.perf_counter()
    checks, mism = [], []
    z_best = None
    for i in sample:
        fz = s298.frozen_theta_support(edges, leaves, dvals, (i,))
        if fz['status'] != 'optimal':
            return dict(status='lp_' + fz['status'], edge=i, closed_seconds=t_closed)
        agrees = (fz['value'] == cf['theta'][i])
        checks.append(dict(edge=i, frozen=fz['value'], closed=cf['theta'][i], agrees=agrees))
        if not agrees:
            mism.append(i)
        if i == best:
            z_best = fz['z']
    t_frozen = time.perf_counter() - t0
    srcs, y = s296.witness_sources((best,), z_best, 1)
    dev = s296.independent_deviation(edges, (best,), z_best, 1, obs, leaves=leaves)
    med = s299.median_observed(dvals)
    return dict(
        status='complete', theta1=theta1, argmin=best,
        attaining=att[:200], n_attaining=len(att),
        attaining_terminal=[i for i in att[:200] if s299.is_terminal(deg, edges[i])],
        argmin_is_terminal=s299.is_terminal(deg, edges[best]),
        argmin_x=cf['x'][best], lp_argmin_x=srcs[0]['x'], pendant=y,
        lp_argmin_theta=next(ch['frozen'] for ch in checks if ch['edge'] == best),
        frozen_checks=checks, n_frozen_checks=len(checks),
        closed_form_mismatches=len(mism), mismatch_examples=mism[:5],
        frozen_theta1_agrees=(next(ch['frozen'] for ch in checks if ch['edge'] == best) == theta1),
        independent_deviation=dev, witness_agrees=(dev is not None and dev == theta1),
        median_observed=med, theta1_over_median_observed=float(theta1 / med) if med else None,
        n_leaves=len(leaves), n_edges=len(edges),
        closed_seconds=t_closed, frozen_seconds=t_frozen)


def main():
    protocol = json.loads(PROTOCOL.read_text())
    stream = OUT / 'stage299_td_runs.jsonl'
    assert not stream.exists(), 'Never overwrite frozen evidence'
    cases = {}
    for row in protocol['cases']:
        if row['arm'] != s299.TD_ARM:
            continue
        cases[row['name']] = s299.prepare_case(
            s299.load_json(row['path'], row['sha256']), row['name'])
    print('loaded %d travellers samples' % len(cases), flush=True)
    jobs = []
    for row in protocol['leafsets']:
        if row['arm'] != s299.TD_ARM:
            continue
        case = cases[row['name']]
        removed = set(int(r) for r in row['removed'])
        leaves = [r for r in sorted(case['observed']) if r not in removed]
        assert len(leaves) == row['n_retained']
        dvals = [case['observed'][r] for r in leaves]
        key = '%s|%s|%s|%s|%d' % (299, row['name'], row['level'], row['regime'],
                                  row['replicate'])
        meta = dict(arm=row['arm'], name=row['name'], stem=row['stem'], family=row['family'],
                    level=row['level'], regime=row['regime'], replicate=row['replicate'],
                    n_removed=row['n_removed'], n_retained=row['n_retained'],
                    removed_fraction=(row['n_removed'] / float(row['n_removed'] + row['n_retained'])),
                    n_edges=row['n_edges'], reference_set_size=row['reference_set_size'],
                    stage298_theta1=row['stage298_theta1'],
                    stage298_argmin=row['stage298_argmin'],
                    detail={k: v for k, v in row['detail'].items() if k != 'source'},
                    sample_seed=int(hashlib.sha256(key.encode()).hexdigest()[:8], 16))
        jobs.append(dict(meta=meta, spec=dict(edges=case['edges'], leaves=leaves, dvals=dvals,
                                              obs={r: d for r, d in zip(leaves, dvals)},
                                              sample_seed=meta['sample_seed'])))
    jobs.sort(key=lambda j: (j['meta']['level'] != 'full', j['meta']['name']))
    print('%d jobs (%d full panels)' % (len(jobs), sum(1 for j in jobs
                                                       if j['meta']['level'] == 'full')), flush=True)
    t0 = time.time()
    res = common.run_pool(jobs, worker, protocol['budgets']['td_wall_seconds'],
                          protocol['budgets']['address_space_bytes'],
                          protocol['budgets']['workers'], stream, progress=50)
    el = time.time() - t0
    ok = [r for r in res if r['result'].get('status') == 'complete']
    # cross-stage check: the full panel must reproduce stage 298 exactly
    cross = []
    for r in ok:
        if r['level'] != 'full':
            continue
        cross.append(dict(name=r['name'], stage299=str(r['result']['theta1']),
                          stage298=r['stage298_theta1'],
                          agrees=(F(r['result']['theta1']) == F(r['stage298_theta1'])),
                          argmin_agrees=(r['result']['argmin'] == r['stage298_argmin'])))
    mism = sum(r['result'].get('closed_form_mismatches', 0) for r in ok)
    checks = sum(r['result'].get('n_frozen_checks', 0) for r in ok)
    summary = dict(stage=299, step='travelers', n_jobs=len(jobs), complete=len(ok),
                   elapsed_seconds=el,
                   frozen_spot_checks=checks, closed_form_mismatches=mism,
                   witness_agree=sum(1 for r in ok if r['result'].get('witness_agrees')),
                   cross_stage298=dict(n=len(cross), agree=sum(c['agrees'] for c in cross),
                                       argmin_agree=sum(c['argmin_agrees'] for c in cross),
                                       rows=cross),
                   failures=[dict(name=r['name'], level=r['level'], regime=r['regime'],
                                  status=r['result'].get('status'),
                                  reason=r['result'].get('reason'))
                             for r in res if r['result'].get('status') != 'complete'])
    (OUT / 'stage299_td_summary.json').write_text(json.dumps(common.serial(summary), indent=1))
    print('done in %.0f s: %d/%d complete; %d frozen spot checks, %d mismatches; '
          'stage-298 agreement %d/%d (argmin %d/%d)'
          % (el, len(ok), len(res), checks, mism,
             summary['cross_stage298']['agree'], len(cross),
             summary['cross_stage298']['argmin_agree'], len(cross)))


if __name__ == '__main__':
    main()

