#!/usr/bin/env python3
"""Stage 308: who actually contributes to the K = 10 threshold correlations.

Two questions the published records answer but the paper did not state.

1. The reported Spearman correlations between panel size and the number of
   references erasable are computed over measurements, one per
   (case, regime, deletion order).  Under the random and whole-clade regimes
   many cases carry two deletion orders, so the measurement count and the case
   count differ and the three regimes do not rest on the same cases.

2. Whether the regime contrast survives giving every case equal weight.  This
   script restricts to the cases with a completed measurement in all three
   regimes, takes the within-case median threshold across deletion orders, and
   recomputes the correlation on that common cohort.

Nothing is re-solved: the script reads the committed bridge records and reuses
stage297_analyze's own de-duplication rule and rank-correlation function.  The
nominal row-level p-values are printed only to document what the naive test
would give; they do not account for repeated deletion orders within a case or
for cases sharing a parent phylogeny, and the paper does not report them.

Writes results/stage308/stage308_cohort_audit.json.
"""
import json
import os
import re
import statistics
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import stage297_analyze as s297                         # noqa: E402

K_TARGET = 10
REGIMES = ('adversarial_nearest', 'random', 'whole_clade')


def nominal_p(rho, n):
    """Two-sided p of the usual t approximation, for documentation only."""
    from scipy import stats
    import math
    if n < 3 or rho is None or abs(rho) >= 1:
        return None
    t = rho * math.sqrt((n - 2) / (1 - rho * rho))
    return float(2 * stats.t.sf(abs(t), n - 2))


def parent_tree(name):
    """Upstream phylogeny, i.e. the case name with its _k<K>_s<seed>... suffix cut."""
    m = re.search(r'_k\d+_s\d+', name)
    return name[:m.start()] if m else name


def main():
    paths = sorted(s297.OUTDIR.glob('stage297_bridge_*.jsonl'))
    rows = s297.bridge_rows(paths)
    good = [r for r in rows if r.get('kind') == 'bridge' and r.get('regimes')]

    # same identity rule as stage297_analyze.threshold_table
    seen = set()
    per_regime = {g: [] for g in REGIMES}
    for r in sorted(good, key=lambda z: z['name']):
        if not r['baseline_unique_true'] or r['K'] != K_TARGET:
            continue
        for g in r['regimes']:
            ident = (r['name'], g['regime'], g['seed'])
            if ident in seen:
                continue
            seen.add(ident)
            t = g['t_unknown_sequence']
            if t is None or g['regime'] not in per_regime:
                continue
            per_regime[g['regime']].append(
                dict(name=r['name'], tree=parent_tree(r['name']),
                     seed=g['seed'], R_size=r['R_size'], t=t))

    rowwise, by_case = {}, {}
    for g, recs in per_regime.items():
        rho = s297.spearman([x['R_size'] for x in recs], [x['t'] for x in recs])
        cases = sorted({x['name'] for x in recs})
        orders = {}
        for x in recs:
            orders.setdefault(x['name'], []).append(x)
        rowwise[g] = dict(
            measurements=len(recs), distinct_cases=len(cases),
            distinct_parent_trees=len({x['tree'] for x in recs}),
            cases_with_two_orders=sum(1 for v in orders.values() if len(v) > 1),
            spearman_rowwise=rho, nominal_p_rowwise=nominal_p(rho, len(recs)))
        by_case[g] = {n: dict(R_size=v[0]['R_size'],
                              t_median=statistics.median([x['t'] for x in v]),
                              n_orders=len(v))
                      for n, v in orders.items()}

    common = sorted(set.intersection(*(set(by_case[g]) for g in REGIMES)))
    cohort = {}
    for g in REGIMES:
        xs = [by_case[g][n]['R_size'] for n in common]
        ys = [by_case[g][n]['t_median'] for n in common]
        cohort[g] = dict(n_cases=len(common),
                         spearman_within_case_median=s297.spearman(xs, ys))

    payload = dict(stage=308, kind='k10_cohort_audit', K=K_TARGET,
                   rowwise=rowwise, common_cohort=cohort,
                   common_cases=common,
                   common_parent_trees=sorted({parent_tree(n) for n in common}))
    out_dir = os.path.join(HERE, os.pardir, 'results', 'stage308')
    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, 'stage308_cohort_audit.json'), 'w') as fh:
        json.dump(payload, fh, indent=1, sort_keys=True)

    print('%-20s %5s %6s %6s %8s %12s %12s' % ('regime', 'meas', 'cases',
                                               'trees', '2 orders', 'rho rowwise',
                                               'nominal p'))
    for g in REGIMES:
        b = rowwise[g]
        print('%-20s %5d %6d %6d %8d %12.4f %12.3g'
              % (g, b['measurements'], b['distinct_cases'],
                 b['distinct_parent_trees'], b['cases_with_two_orders'],
                 b['spearman_rowwise'], b['nominal_p_rowwise']))
    print('\ncommon cohort: %d cases measured in all three regimes, '
          '%d parent trees' % (len(common), len(payload['common_parent_trees'])))
    for g in REGIMES:
        print('  %-20s rho on within-case medians %.4f'
              % (g, cohort[g]['spearman_within_case_median']))


if __name__ == '__main__':
    main()
