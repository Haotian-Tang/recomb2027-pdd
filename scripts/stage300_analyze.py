#!/usr/bin/env python3
"""Stage 300 analysis: the three-way comparison, in both norms.

Reads the per-run fit records and the per-case baselines and writes
results/stage300/stage300_analysis.json.  Every comparison of an L-infinity
quantity against Theta_1 is an exact rational decision; quantiles are reported
as floats.  One run per (case, arm) is designated canonical - stage 289's first
repetition where it exists, then stage 294, then stage 297's default seed, then
this stage's own runs - and the spread across the remaining replicates is
reported instead of being averaged away.
"""
import json
import sys
from collections import Counter, defaultdict
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import stage300_common as s300                   # noqa: E402

STAGE_ORDER = {289: 0, 294: 1, 297: 2, 300: 3}


def quantiles(values):
    if not values:
        return None
    s = sorted(float(v) for v in values)
    def q(f):
        if len(s) == 1:
            return s[0]
        position = f * (len(s) - 1)
        low = int(position)
        high = min(low + 1, len(s) - 1)
        return s[low] + (s[high] - s[low]) * (position - low)
    return dict(n=len(s), min=s[0], q25=q(0.25), median=q(0.5), q75=q(0.75), max=s[-1])


def load(collection):
    fits = [json.loads(l) for l in open(s300.OUT / ('stage300_fits_%s.jsonl' % collection))]
    bases = [json.loads(l) for l in open(s300.OUT / ('stage300_baselines_%s.jsonl' % collection))]
    return fits, {b['name']: b for b in bases}


def canonical(runs):
    done = [r for r in runs if r.get('evaluated')]
    if not done:
        return None
    return sorted(done, key=lambda r: (STAGE_ORDER.get(r['source_stage'], 9),
                                       r.get('source_line') or 0))[0]


def analyse(collection):
    fits, bases = load(collection)
    by_arm = defaultdict(lambda: defaultdict(list))
    for row in fits:
        by_arm[row['arm']][row['name']].append(row)
    out = dict(collection=collection, n_cases=len(bases), arms={},
               availability={}, per_case={})
    for arm, cases in sorted(by_arm.items()):
        rows, failures, replicate_spread = [], Counter(), dict(k=Counter(), decision=Counter())
        for name, runs in sorted(cases.items()):
            base = bases.get(name)
            best = canonical(runs)
            for r in runs:
                if not r.get('evaluated'):
                    failures[(r.get('status') or 'unknown')] += 1
            if best is None:
                continue
            theta1 = F(base['theta1']) if base and base['theta1'] else None
            r_dd = F(best['r_dd'])
            lam = (base.get('lambda1') or {}) if base else {}
            truth = (base.get('truth') or {}) if base else {}
            entry = dict(
                name=name, arm=arm, k=best['k'], n_leaves=best['n_leaves'],
                source_stage=best['source_stage'], source_file=best['source_file'],
                r_dd=str(r_dd), r_dd_float=float(r_dd),
                theta1_float=(float(theta1) if theta1 is not None else None),
                r_dd_ge_theta1=(None if theta1 is None else bool(r_dd >= theta1)),
                r_dd_over_theta1=(None if not theta1 else float(r_dd / theta1)),
                ssr_dd=best['ssr_dd'], reported_loss=best['reported_loss'],
                loss_rel_diff=best.get('loss_recompute_rel_diff'),
                lambda1_ssr=lam.get('ssr'),
                ssr_over_lambda1=(None if not lam.get('ssr') else best['ssr_dd'] / lam['ssr']),
                dd_beats_one_source_in_ssr=(None if not lam.get('ssr')
                                            else best['ssr_dd'] < lam['ssr']),
                lambda1_linf_float=(float(F(lam['linf_at_l2_optimum']))
                                    if lam.get('linf_at_l2_optimum') else None),
                truth_linf=truth.get('linf'), truth_ssr=truth.get('ssr'),
                truth_k=truth.get('k_true'),
                truth_matches_stage298=truth.get('truth_linf_matches_stage298'),
                min_p=best.get('min_p'), sum_p=best.get('sum_p'),
                stop_reason=best.get('stop_reason'),
                writer_crashed=best.get('jplace_writer_crashed'),
                n_runs=len(runs), n_evaluated=sum(1 for r in runs if r.get('evaluated')))
            if truth.get('linf'):
                tl = F(truth['linf'])
                entry.update(
                    dd_beats_truth_in_linf=bool(r_dd < tl),
                    truth_ge_theta1=(None if theta1 is None else bool(tl >= theta1)),
                    truth_over_theta1=(None if not theta1 else float(tl / theta1)),
                    dd_beats_truth_in_ssr=(None if truth.get('ssr') is None
                                           else best['ssr_dd'] < truth['ssr']),
                    truth_beats_one_source_in_ssr=(
                        None if truth.get('ssr') is None or not lam.get('ssr')
                        else truth['ssr'] < lam['ssr']),
                    ordering_linf=''.join(sorted(
                        [('T', float(tl)), ('D', float(r_dd)),
                         ('1', float(theta1) if theta1 is not None else float('inf'))],
                        key=lambda z: z[1])[i][0] for i in range(3)))
            rows.append(entry)
            ks = {r['k'] for r in runs if r.get('evaluated')}
            decisions = {r.get('r_dd_ge_theta1') for r in runs if r.get('evaluated')}
            replicate_spread['k'][len(ks)] += 1
            replicate_spread['decision'][len(decisions)] += 1
        cell = dict(
            n_cases_with_answer=len(rows),
            availability=len(rows) / max(1, len(bases)),
            failures={str(k): v for k, v in failures.items()},
            k_returned=quantiles([r['k'] for r in rows]),
            n_r_dd_ge_theta1=sum(1 for r in rows if r['r_dd_ge_theta1']),
            n_r_dd_lt_theta1=sum(1 for r in rows if r['r_dd_ge_theta1'] is False),
            r_dd_over_theta1=quantiles([r['r_dd_over_theta1'] for r in rows
                                        if r['r_dd_over_theta1'] is not None]),
            ssr_over_lambda1=quantiles([r['ssr_over_lambda1'] for r in rows
                                        if r['ssr_over_lambda1'] is not None]),
            n_dd_beats_one_source_in_ssr=sum(1 for r in rows
                                             if r['dd_beats_one_source_in_ssr']),
            loss_rel_diff=quantiles([r['loss_rel_diff'] for r in rows
                                     if r['loss_rel_diff'] is not None]),
            replicates=dict(distinct_k=dict(replicate_spread['k']),
                            distinct_decision=dict(replicate_spread['decision'])),
            stop_reasons=dict(Counter(r['stop_reason'] for r in rows)),
            writer_crashed=sum(1 for r in rows if r['writer_crashed']))
        truth_rows = [r for r in rows if r.get('truth_linf')]
        if truth_rows:
            cell.update(
                n_with_truth=len(truth_rows),
                n_dd_beats_truth_in_linf=sum(1 for r in truth_rows
                                             if r['dd_beats_truth_in_linf']),
                n_dd_beats_truth_in_ssr=sum(1 for r in truth_rows
                                            if r['dd_beats_truth_in_ssr']),
                n_truth_ge_theta1=sum(1 for r in truth_rows if r['truth_ge_theta1']),
                n_truth_beats_one_source_in_ssr=sum(
                    1 for r in truth_rows if r['truth_beats_one_source_in_ssr']),
                truth_over_theta1=quantiles([r['truth_over_theta1'] for r in truth_rows
                                             if r['truth_over_theta1'] is not None]),
                ordering_linf=dict(Counter(r['ordering_linf'] for r in truth_rows)),
                truth_linf_matches_stage298=dict(Counter(
                    str(r['truth_matches_stage298']) for r in truth_rows)))
        out['arms'][arm] = cell
        out['per_case'][arm] = rows
    # baseline-only facts, arm independent
    base_rows = list(bases.values())
    out['baselines'] = dict(
        n=len(base_rows),
        theta1=quantiles([float(F(b['theta1'])) for b in base_rows if b['theta1']]),
        lambda1_ssr=quantiles([b['lambda1']['ssr'] for b in base_rows if b.get('lambda1')]),
        lambda1_ssr_recompute_max_rel=max(
            [abs(b['lambda1']['ssr_recomputed'] - b['lambda1']['ssr'])
             / max(1e-300, abs(b['lambda1']['ssr']))
             for b in base_rows if b.get('lambda1')], default=None),
        truth_linf_matches_stage298=dict(Counter(
            str((b.get('truth') or {}).get('truth_linf_matches_stage298'))
            for b in base_rows)))
    return out


def extra_arms(collection):
    """Summaries of the two fairness arms: the re-solved QP and the panel used."""
    out = {}
    resolved_path = s300.OUT / ('stage300_resolved_%s.jsonl' % collection)
    if resolved_path.exists():
        rows = [json.loads(l) for l in open(resolved_path)]
        cell = {}
        for arm in sorted({r['arm'] for r in rows}):
            sel = [r for r in rows if r['arm'] == arm]
            ok = [r for r in sel if r['resolved']['status'] == 'ok']
            rel = [abs(r['resolved']['ssr'] - r['reported_loss'])
                   / max(abs(r['reported_loss']), 1e-300)
                   for r in ok if r['reported_loss']]
            cell[arm] = dict(
                n=len(sel), statuses=dict(Counter(r['resolved']['status'] for r in sel)),
                n_ok=len(ok),
                n_linf_ge_theta1=sum(1 for r in ok if r.get('linf_ge_theta1')),
                linf_over_theta1=quantiles([r['linf_over_theta1'] for r in ok
                                            if r.get('linf_over_theta1') is not None]),
                resolved_ssr=quantiles([r['resolved']['ssr'] for r in ok]),
                resolved_vs_reported_loss_rel=quantiles(rel))
        out['resolved'] = cell
    prune_path = s300.OUT / ('stage300_prune_%s.jsonl' % collection)
    if prune_path.exists():
        rows = [json.loads(l) for l in open(prune_path)]
        dropped = [r for r in rows if r['n_kept'] < r['n_observed']]
        cell = dict(n_cases=len(rows), n_cases_with_references_dropped=len(dropped),
                    kept_fraction_when_dropped=quantiles([r['kept_fraction'] for r in dropped]),
                    n_kept_when_dropped=quantiles([r['n_kept'] for r in dropped]),
                    n_observed_when_dropped=quantiles([r['n_observed'] for r in dropped]),
                    dropped_cases=[dict(name=r['name'], n_observed=r['n_observed'],
                                        n_kept=r['n_kept'],
                                        n_with_distance=r['survive']['n_with_distance'])
                                   for r in dropped],
                    arms={})
        for arm in ('dd_default', 'dd_true_k'):
            sel = [r for r in rows if arm in r['arms']
                   and r['arms'][arm].get('r_dd_kept_over_theta1_kept') is not None]
            if not sel:
                continue
            cell['arms'][arm] = dict(
                n=len(sel),
                n_r_dd_kept_ge_theta1_kept=sum(
                    1 for r in sel if r['arms'][arm]['r_dd_kept_ge_theta1_kept']),
                ratio=quantiles([r['arms'][arm]['r_dd_kept_over_theta1_kept'] for r in sel]),
                dropped_only=dict(
                    n=len([r for r in sel if r['n_kept'] < r['n_observed']]),
                    n_ge=sum(1 for r in sel if r['n_kept'] < r['n_observed']
                             and r['arms'][arm]['r_dd_kept_ge_theta1_kept'])))
        out['fair_panel'] = cell
    return out


def main():
    collections = [c for c in s300.COLLECTIONS
                   if (s300.OUT / ('stage300_fits_%s.jsonl' % c)).exists()]
    report = dict(stage=300, step='analysis', collections=collections, by_collection={})
    for collection in collections:
        report['by_collection'][collection] = analyse(collection)
        report['by_collection'][collection].update(extra_arms(collection))
        cell = report['by_collection'][collection]
        print('==', collection, 'cases', cell['n_cases'])
        for arm, a in cell['arms'].items():
            print('  %-14s answer=%3d/%3d  R>=Theta1 %3d/%3d  ratio med %s  '
                  'SSR<Lambda1 %3d  k med %s'
                  % (arm, a['n_cases_with_answer'], cell['n_cases'],
                     a['n_r_dd_ge_theta1'], a['n_cases_with_answer'],
                     None if not a['r_dd_over_theta1'] else round(a['r_dd_over_theta1']['median'], 3),
                     a['n_dd_beats_one_source_in_ssr'],
                     None if not a['k_returned'] else a['k_returned']['median']))
    (s300.OUT / 'stage300_analysis.json').write_text(json.dumps(s300.serial(report), indent=1))
    print('wrote', s300.OUT / 'stage300_analysis.json')


if __name__ == '__main__':
    main()
