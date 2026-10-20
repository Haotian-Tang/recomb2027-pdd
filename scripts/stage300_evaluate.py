#!/usr/bin/env python3
"""Stage 300 step 2: evaluate every configuration's fit, exactly.

For one collection this writes two record files:

  stage300_fits_<collection>.jsonl       one row per DecoDiPhy run (harvested or
                                         new): R_dd = max_r |pred_r - d_r| as an
                                         exact rational, SSR recomputed from the
                                         returned configuration, the tool's own
                                         reported loss, and the exact decision
                                         R_dd >= Theta_1.
  stage300_baselines_<collection>.jsonl  one row per case: Theta_1 (imported from
                                         stage 298), the TRUE configuration's
                                         L-infinity residual (imported, and
                                         independently recomputed here) with its
                                         SSR, and Lambda_1 - the exhaustive
                                         minimum of DecoDiPhy's own objective
                                         over single sources.

Theta_1 and the truth residual are IMPORTED from stage 298 and never recomputed
differently; the recomputation that does happen is a verification and is
reported as such (`truth_linf_matches_stage298`).
"""
import argparse
import json
import sys
import time
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import stage300_common as s300                   # noqa: E402


def configuration(run):
    """(sources, pendant) exactly as the tool returned them, as exact rationals."""
    sources = [dict(edge=int(s['edge']), p=F(float(s['p'])), x=F(float(s['x'])))
               for s in run['sources']]
    return sources, F(float(run['pendant_total'] or 0.0))


def runs_for(collection):
    """Harvested runs plus any new stage-300 runs, tagged by origin."""
    rows = []
    harvest = s300.OUT / 'stage300_harvest.jsonl'
    if harvest.exists():
        with open(harvest) as handle:
            for line in handle:
                d = json.loads(line)
                if d['collection'] == collection:
                    d['origin'] = 'harvested'
                    rows.append(d)
    for new in sorted(s300.OUT.glob('stage300_runs_*.jsonl')):
        with open(new) as handle:
            for line_no, line in enumerate(handle, 1):
                d = json.loads(line)
                if d.get('collection') != collection:
                    continue
                result = d.get('result', {})
                rows.append(dict(kind='dd_new', stage=300, collection=collection,
                                 family=d.get('family'), name=d['name'],
                                 arm=d['method'], source_stage=300,
                                 source_file='results/stage300/' + new.name,
                                 source_line=line_no, tag=dict(),
                                 status=result.get('status'), k=result.get('k'),
                                 reported_loss=result.get('loss'),
                                 sources=result.get('sources'),
                                 pendant_total=result.get('pendant_total'),
                                 stop_reason=result.get('stop_reason'),
                                 ladder=result.get('ladder'),
                                 jplace_writer_crashed=result.get('jplace_writer_crashed'),
                                 writer_exception=result.get('writer_exception'),
                                 unmapped_anchors=result.get('unmapped_anchors'),
                                 exception=result.get('exception'),
                                 argv_tail=[a for a in (result.get('argv') or [])
                                            if not a.startswith('/')],
                                 solver_wall_seconds=result.get('solver_wall_seconds'),
                                 harness_wall_seconds=d.get('harness_wall_seconds'),
                                 maximum_rss_kib=d.get('maximum_rss_kib'),
                                 timed_out=d.get('timed_out'), origin='stage300'))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--collection', required=True)
    ap.add_argument('--limit', type=int, default=0)
    ap.add_argument('--skip-lambda1', action='store_true')
    args = ap.parse_args()

    collection = args.collection
    theta = s300.theta1_table()
    truth_table = s300.truth_residual_table()
    entries = s300.case_entries(collection)
    if args.limit:
        entries = entries[:args.limit]
    by_case = {}
    for row in runs_for(collection):
        by_case.setdefault(row['name'], []).append(row)

    fits_path = s300.OUT / ('stage300_fits_%s.jsonl' % collection)
    base_path = s300.OUT / ('stage300_baselines_%s.jsonl' % collection)
    fits = open(fits_path, 'w')
    bases = open(base_path, 'w')
    t_start = time.perf_counter()
    for index, entry in enumerate(entries, 1):
        name = entry['name']
        case = s300.load_case(entry)
        edges, leaves, dvals = s300.case_arrays(case)
        t1 = theta.get(name)
        theta1 = t1['theta1'] if t1 else None
        median_obs = t1['own_median_observed'] if t1 else s300.median_observed(dvals)
        t0 = time.perf_counter()

        # ---- baselines: truth and the best single source in DecoDiPhy's norm
        truth = s300.truth_configuration(case)
        truth_row = dict(status='no_truth')
        if truth is not None:
            metrics = s300.fit_metrics(edges, leaves, dvals, truth[0], truth[1])
            imported = truth_table.get(name)
            truth_row = dict(status='ok', k_true=len(truth[0]),
                             linf=metrics['linf'], ssr=metrics['ssr'],
                             median_abs=metrics['median_abs'],
                             imported_residual_max=(imported['residual_max']
                                                    if imported else None),
                             truth_linf_matches_stage298=(
                                 None if not imported
                                 else metrics['linf'] == imported['residual_max']))
        lam = None
        if not args.skip_lambda1:
            lam = s300.best_single_source_l2(edges, leaves, dvals)
            lam_sources = [dict(edge=lam['edge'], p=F(1), x=F(lam['x']))]
            lam_metrics = s300.fit_metrics(edges, leaves, dvals, lam_sources, F(lam['y']))
            lam = dict(lam, linf_at_l2_optimum=lam_metrics['linf'],
                       ssr_recomputed=lam_metrics['ssr'])
        bases.write(json.dumps(s300.serial(dict(
            kind='case_baseline', stage=300, collection=collection, name=name,
            family=entry['family'], path=entry['path'], sha256=entry['sha256'],
            n_edges=len(edges), n_leaves=len(leaves),
            theta1=theta1, theta1_source=(t1['source'] if t1 else None),
            theta_ladder=(t1.get('theta_ladder') if t1 else None),
            own_median_observed=median_obs, truth=truth_row, lambda1=lam,
            seconds=time.perf_counter() - t0))) + '\n')
        bases.flush()

        # ---- every DecoDiPhy run on this case
        for run in by_case.get(name, []):
            row = dict(kind='fit', stage=300, collection=collection, name=name,
                       arm=run['arm'], origin=run['origin'],
                       source_stage=run['source_stage'], source_file=run['source_file'],
                       source_line=run['source_line'], tag=run.get('tag'),
                       status=run['status'], k=run.get('k'),
                       reported_loss=run.get('reported_loss'),
                       stop_reason=run.get('stop_reason'),
                       jplace_writer_crashed=run.get('jplace_writer_crashed'),
                       unmapped_anchors=run.get('unmapped_anchors'),
                       exception=run.get('exception'),
                       solver_wall_seconds=run.get('solver_wall_seconds'),
                       theta1=theta1, n_leaves=len(leaves))
            if run['status'] == 'complete' and run.get('sources'):
                sources, pendant = configuration(run)
                metrics = s300.fit_metrics(edges, leaves, dvals, sources, pendant)
                sum_p = float(sum(float(s['p']) for s in run['sources']))
                row.update(evaluated=True, r_dd=metrics['linf'], ssr_dd=metrics['ssr'],
                           median_abs_dd=metrics['median_abs'], sum_p=sum_p,
                           pendant_total=run.get('pendant_total'),
                           distinct_edges=len({int(s['edge']) for s in run['sources']}),
                           support=sorted(int(s['edge']) for s in run['sources']),
                           min_p=min(float(s['p']) for s in run['sources']),
                           r_dd_float=float(metrics['linf']),
                           r_dd_ge_theta1=(None if theta1 is None
                                           else bool(metrics['linf'] >= theta1)),
                           r_dd_over_theta1=(None if theta1 is None or theta1 == 0
                                             else float(metrics['linf'] / theta1)),
                           loss_recompute_abs_diff=(
                               None if run.get('reported_loss') is None
                               else abs(metrics['ssr'] - float(run['reported_loss']))),
                           loss_recompute_rel_diff=(
                               None if not run.get('reported_loss')
                               else abs(metrics['ssr'] - float(run['reported_loss']))
                               / abs(float(run['reported_loss']))))
            else:
                row.update(evaluated=False)
            fits.write(json.dumps(s300.serial(row)) + '\n')
        fits.flush()
        if index % 10 == 0 or index == len(entries):
            print('[%s] %s %d/%d  %.1f s elapsed'
                  % (time.strftime('%H:%M:%S'), collection, index, len(entries),
                     time.perf_counter() - t_start), flush=True)
    fits.close()
    bases.close()
    print('wrote %s and %s' % (fits_path, base_path), flush=True)


if __name__ == '__main__':
    main()
