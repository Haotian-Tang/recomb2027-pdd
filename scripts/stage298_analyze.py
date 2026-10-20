#!/usr/bin/env python3
"""Stage 298 analysis: assemble the real-data tolerance axis."""
import json
import os
import sys
from collections import Counter
from fractions import Fraction as F

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
import stage288_common as common            # noqa: E402
import stage298_common as s298              # noqa: E402

OUT = s298.OUT


def q(vals, qs=(0.0, 0.25, 0.5, 0.75, 1.0)):
    """Linearly interpolated quantiles (the numpy.percentile convention).

    NOT the sorted(...)[n // 2] order statistic: for even n that picks the upper
    middle element rather than the midpoint, which is the defect a standing
    review finding raised against the stage-296 report's even-n subgroup
    medians.  The stage-296 DATA convention for `median_observed` - the axis
    denominator - is deliberately left untouched in stage298_common so the
    residual extension still reproduces stage 296's rows exactly; only the
    reporting quantiles here are interpolated.
    """
    if not vals:
        return None
    s = sorted(vals)
    out = {}
    for p in qs:
        pos = p * (len(s) - 1)
        lo = int(pos)
        hi = min(lo + 1, len(s) - 1)
        out[str(p)] = s[lo] + (pos - lo) * (s[hi] - s[lo])
    return out


def m_of_eps(thetas, eps, censored_from):
    """m(eps) = min{k : Theta_k < eps}, stage 296's strict convention."""
    for k in sorted(thetas):
        if thetas[k] < eps:
            return k, 'exact'
    if censored_from is None:
        return None, 'exhausted'
    return censored_from, 'lower_bound'


def main():
    protocol = json.loads((OUT / 'stage298_protocol.json').read_text())
    runs = [json.loads(l) for l in open(OUT / 'stage298_sweep_runs.jsonl')]
    censored = json.loads((OUT / 'stage298_sweep_censored.json').read_text())
    resid = {r['case']: r for r in json.loads((OUT / 'stage298_residuals.json').read_text())['rows']}
    manifest = json.loads((OUT / 'stage298_placement_manifest.json').read_text())
    td = json.loads((OUT / 'stage298_td_summary.json').read_text())['rows']
    meta = {s['name']: s for s in protocol['small_cases']}

    theta = {}
    verif = dict(k1_edges_compared=0, k1_closed_form_mismatches=0,
                 witnesses_checked=0, witnesses_agree=0, cells_complete=0,
                 cells_failed=0)
    for r in runs:
        a = r['result']
        if a['status'] != 'complete':
            verif['cells_failed'] += 1
            continue
        verif['cells_complete'] += 1
        theta.setdefault(r['name'], {})[r['k']] = F(a['theta'])
        if r['k'] == 1:
            verif['k1_edges_compared'] += r['n_edges']
            verif['k1_closed_form_mismatches'] += a['closed_form_mismatches']
        if a.get('witness_agrees') is not None:
            verif['witnesses_checked'] += 1
            verif['witnesses_agree'] += int(a['witness_agrees'])
    cens_k = {}
    for c in censored['censored']:
        cens_k.setdefault(c['name'], set()).add(c['k'])

    # ------------------------------------------------------------ small cases
    table = []
    for name, th in sorted(theta.items()):
        m = meta[name]
        rr = resid.get(name, {})
        med_own = F(rr['own_median_observed']) if rr.get('own_median_observed') else None
        cens = sorted(cens_k.get(name, []))
        first_cens = cens[0] if cens else None
        row = dict(name=name, family=m['family'], instance=m['instance'],
                   K_generating=m['K_generating'], N=m['reference_set_size'],
                   n_edges=m['n_edges'],
                   own_median_observed=str(med_own) if med_own else None,
                   theta={str(k): str(v) for k, v in sorted(th.items())},
                   theta_float={str(k): float(v) for k, v in sorted(th.items())},
                   theta_rel_own={str(k): float(v / med_own) for k, v in sorted(th.items())}
                   if med_own else None,
                   censored_from_k=first_cens)
        if rr.get('status') == 'ok':
            rmed, rmax = F(rr['residual_median']), F(rr['residual_max'])
            row.update(residual_median=str(rmed), residual_max=str(rmax),
                       residual_median_float=float(rmed), residual_max_float=float(rmax),
                       residual_max_rel_own=rr['residual_max_rel_own'],
                       residual_median_rel_own=rr['residual_median_rel_own'],
                       theta1_over_residual_max=float(th[1] / rmax) if rmax else None,
                       residual_max_over_theta1=float(rmax / th[1]),
                       residual_median_over_theta1=float(rmed / th[1]),
                       one_source_within_max_residual=bool(th[1] < rmax),
                       one_source_within_median_residual=bool(th[1] < rmed))
            mm, kind = m_of_eps(th, rmax, first_cens)
            row['m_at_max_residual'], row['m_at_max_residual_kind'] = mm, kind
            mm, kind = m_of_eps(th, rmed, first_cens)
            row['m_at_median_residual'], row['m_at_median_residual_kind'] = mm, kind
        table.append(row)

    by_family = {}
    for fam in sorted({r['family'] for r in table}):
        sel = [r for r in table if r['family'] == fam]
        ok = [r for r in sel if 'residual_max' in r]
        by_family[fam] = dict(
            n=len(sel),
            n_with_residual=len(ok),
            theta1_abs=q([r['theta_float']['1'] for r in sel]),
            theta1_rel_own=q([r['theta_rel_own']['1'] for r in sel]),
            theta2_over_theta1=q([r['theta_float']['2'] / r['theta_float']['1']
                                  for r in sel if '2' in r['theta_float']]),
            theta3_over_theta1=q([r['theta_float']['3'] / r['theta_float']['1']
                                  for r in sel if '3' in r['theta_float']]),
            n_theta3=sum(1 for r in sel if '3' in r['theta_float']),
            residual_max_rel_own=q([r['residual_max_rel_own'] for r in ok]),
            residual_median_rel_own=q([r['residual_median_rel_own'] for r in ok]),
            residual_max_over_theta1=q([r['residual_max_over_theta1'] for r in ok]),
            residual_median_over_theta1=q([r['residual_median_over_theta1'] for r in ok]),
            n_one_source_within_max_residual=sum(1 for r in ok
                                                 if r['one_source_within_max_residual']),
            n_one_source_within_median_residual=sum(1 for r in ok
                                                    if r['one_source_within_median_residual']),
            m_at_max_residual=dict(Counter(str(r['m_at_max_residual']) + ':' +
                                           r['m_at_max_residual_kind'] for r in ok)),
            m_at_median_residual=dict(Counter(str(r['m_at_median_residual']) + ':' +
                                              r['m_at_median_residual_kind'] for r in ok)),
            accuracy_factor_to_leave_collapsed=q(
                [r['residual_max_over_theta1'] for r in ok
                 if r['one_source_within_max_residual']]),
            by_K={str(k): dict(
                n=sum(1 for r in sel if r['K_generating'] == k),
                theta1_rel_own=q([r['theta_rel_own']['1'] for r in sel
                                  if r['K_generating'] == k]),
                n_one_source_within_max_residual=sum(
                    1 for r in ok if r['K_generating'] == k
                    and r['one_source_within_max_residual']))
                for k in sorted({r['K_generating'] for r in sel})})

    # ------------------------------------------------- travellers' diarrhoea
    td_rows = []
    for r in td:
        td_rows.append(dict(
            name=r['name'], sample=r['sample'], group=r['group'],
            n_edges=r['n_edges'], n_leaves=r['n_leaves'],
            median_observed=r['median_observed'],
            theta1=r['theta1'], theta1_float=float(F(r['theta1'])),
            theta1_rel=r['theta1_over_median_observed'],
            argmin_edge=r['argmin_edge'], argmin_endpoints=r['argmin_endpoints'],
            argmin_is_terminal=r['argmin_is_terminal'],
            argmin_terminal_leaf_labels=r['argmin_terminal_leaf_labels'],
            argmin_x=float(F(r['argmin_x'])) if r['argmin_x'] is not None else None,
            closed_form_mismatches=r['closed_form_mismatches'],
            edges_covered=r['edges_covered'],
            witness_agrees=r['verification'].get('witness_agrees'),
            frozen_sweep_seconds=r['frozen_sweep_seconds']))
    # transfer the measured krepp relative error onto the travellers' axis
    kr = [r for r in table if r['family'] == 'krepp_estimated_10k' and 'residual_max' in r]
    kr_max_rel = sorted(r['residual_max_rel_own'] for r in kr)
    kr_med_rel = sorted(r['residual_median_rel_own'] for r in kr)
    ref = dict(krepp_residual_max_rel_median=q(kr_max_rel, (0.5,))['0.5'],
               krepp_residual_median_rel_median=q(kr_med_rel, (0.5,))['0.5'])
    td_summary = dict(
        n=len(td_rows),
        theta1_abs=q([r['theta1_float'] for r in td_rows]),
        theta1_rel=q([r['theta1_rel'] for r in td_rows]),
        median_observed=q([float(F(r['median_observed'])) for r in td_rows]),
        n_argmin_terminal=sum(1 for r in td_rows if r['argmin_is_terminal']),
        distinct_argmin_edges=len({r['argmin_edge'] for r in td_rows}),
        closed_form_mismatches=sum(r['closed_form_mismatches'] for r in td_rows),
        edges_decided_twice=sum(r['edges_covered'] for r in td_rows),
        witnesses_agree=sum(1 for r in td_rows if r['witness_agrees']),
        reference_relative_error=ref,
        n_below_krepp_max_rel=sum(1 for r in td_rows
                                  if r['theta1_rel'] < ref['krepp_residual_max_rel_median']),
        n_below_krepp_median_rel=sum(1 for r in td_rows
                                     if r['theta1_rel'] < ref['krepp_residual_median_rel_median']),
        accuracy_factor_vs_krepp_max=q([ref['krepp_residual_max_rel_median'] / r['theta1_rel']
                                        for r in td_rows]),
        theta2_status='not_computed: see cost_of_theta2_on_travelers')

    report = dict(
        stage=298, track='real_data_tolerance_axis',
        quantile_convention=('reported quantiles are linearly interpolated '
                             '(numpy.percentile convention); the axis denominator '
                             'median_observed keeps stage 296\'s sorted(...)[n // 2] '
                             'data convention so the residual extension reproduces '
                             'stage 296\'s rows exactly'),
        counts=protocol['counts'],
        verification=dict(
            verif,
            placement_conversion=dict(
                n_with_placement_file=manifest['n_instances_with_placement_file'],
                converted=manifest['n_converted'], refused=manifest['n_refused'],
                refusal_reasons=dict(Counter(r['reason'] for r in manifest['refused'])),
                dropped_entries_are_exactly_the_queries=sum(
                    1 for c in manifest['cases']
                    if c['provenance']['dropped_entries_are_exactly_the_queries'])),
            residual_extension=json.loads(
                (OUT / 'stage298_residuals.json').read_text())['validation_against_stage296'],
            travelers_closed_form_mismatches=td_summary['closed_form_mismatches'],
            travelers_edges_decided_twice=td_summary['edges_decided_twice']),
        censoring=dict(
            n_cells=len(censored['censored']),
            by_k=dict(Counter(str(c['k']) for c in censored['censored'])),
            by_reason=dict(Counter(c['reason'] for c in censored['censored'])),
            support_cap=protocol['support_cap'],
            measured_seconds_per_support=q([r['result']['seconds_per_support']
                                            for r in runs
                                            if r['result']['status'] == 'complete']),
            cost_of_theta2_on_travelers=None),
        by_family=by_family, travelers=td_summary, travelers_rows=td_rows,
        small_cases=table)
    (OUT / 'stage298_analysis.json').write_text(
        json.dumps(common.serial(report), indent=1))
    print('wrote stage298_analysis.json')
    for fam, d in by_family.items():
        print('%-22s n=%3d  Theta_1/med median %.4f  resid_max/med %.4f  '
              'one-source-within-max-residual %d/%d'
              % (fam, d['n'], d['theta1_rel_own']['0.5'],
                 (d['residual_max_rel_own'] or {}).get('0.5', float('nan')),
                 d['n_one_source_within_max_residual'], d['n_with_residual']))
    print('travelers: Theta_1/median_observed quartiles',
          {k: round(v, 5) for k, v in td_summary['theta1_rel'].items()})
    print('travelers below krepp max-residual line: %d/%d; below median-residual line: %d/%d'
          % (td_summary['n_below_krepp_max_rel'], td_summary['n'],
             td_summary['n_below_krepp_median_rel'], td_summary['n']))


if __name__ == '__main__':
    main()

