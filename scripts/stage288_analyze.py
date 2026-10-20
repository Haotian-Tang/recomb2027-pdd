#!/usr/bin/env python3
"""Stage 288 analysis: scaling (steps 1-2) and reference incompleteness (3-4).

Reads only stage 288 result files plus the frozen protocols. Timeouts are never
imputed; every rate is reported over the jobs that actually ran.
"""
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median

import stage288_common as common
from stage288_common import ROOT, serial, sha

OUT = ROOT / 'results/stage288_analysis.json'


def ols(xs, ys):
    n = len(xs)
    if n < 3:
        return None
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    if sxx == 0:
        return None
    slope = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx
    intercept = my - slope * mx
    sst = sum((y - my) ** 2 for y in ys)
    ssr = sum((y - (slope * x + intercept)) ** 2 for x, y in zip(xs, ys))
    return dict(n=n, slope=slope, intercept=intercept,
                r_squared=(1 - ssr / sst) if sst > 0 else None)


def quantiles(values):
    if not values:
        return None
    s = sorted(values)
    def q(p):
        if len(s) == 1:
            return s[0]
        position = p * (len(s) - 1)
        low = int(math.floor(position))
        high = min(low + 1, len(s) - 1)
        return s[low] + (s[high] - s[low]) * (position - low)
    return dict(n=len(s), min=s[0], q25=q(0.25), median=q(0.5), q75=q(0.75), q90=q(0.9), max=s[-1])


# ------------------------------------------------------------------ scaling

def scaling_analysis():
    data = json.loads((ROOT / 'results/stage288_scaling.json').read_text())
    records = data['records']
    tier1 = [r for r in records if r['tier'] == 1]
    out = dict(jobs_total=len(records),
               jobs_by_tier=dict(Counter(str(r['tier']) for r in records)),
               status_by_tier_method={})
    for tier in sorted({r['tier'] for r in records}):
        for method in sorted({r['method'] for r in records if r['tier'] == tier}):
            subset = [r for r in records if r['tier'] == tier and r['method'] == method]
            out['status_by_tier_method'][f'tier{tier}:{method}'] = dict(
                Counter(r['result']['status'] for r in subset))

    # per-case, per-method tier-1 completion of BOTH repetitions
    by_case = defaultdict(dict)
    for r in tier1:
        by_case[r['name']].setdefault(r['method'], []).append(r)
    meta = {r['name']: r for r in tier1}
    complete, cpu = defaultdict(dict), defaultdict(dict)
    for name, methods in by_case.items():
        for method, runs in methods.items():
            done = [x for x in runs if x['result']['status'] == 'complete']
            complete[name][method] = len(done) == len(runs) and len(runs) == 2
            if complete[name][method]:
                cpu[name][method] = median(float(x['result']['child_cpu_seconds']) for x in runs)

    cells = defaultdict(lambda: dict(cases=0, enum=0, struct=0, N=[], K=None))
    for name in by_case:
        cell = meta[name]['cell']
        cells[cell]['cases'] += 1
        cells[cell]['N'].append(meta[name]['N'])
        cells[cell]['K'] = meta[name]['K']
        cells[cell]['enum'] += int(complete[name].get('original_support', False))
        cells[cell]['struct'] += int(complete[name].get('stage43_equality_reduced', False))
    cell_table = {}
    for cell in sorted(cells):
        v = cells[cell]
        cell_table[cell] = dict(cases=v['cases'], K=v['K'], N_min=min(v['N']), N_max=max(v['N']),
                                enumeration_complete=v['enum'], structural_complete=v['struct'],
                                enumeration_rate=v['enum'] / v['cases'],
                                structural_rate=v['struct'] / v['cases'])
    out['tier1_cells'] = cell_table
    out['tier1_frontier_cells_enumeration_zero_structural_positive'] = sorted(
        c for c, v in cell_table.items()
        if v['enumeration_rate'] == 0 and v['structural_rate'] > 0)
    out['tier1_cells_both_zero'] = sorted(
        c for c, v in cell_table.items() if v['enumeration_rate'] == 0 and v['structural_rate'] == 0)

    joint = [n for n in by_case
             if complete[n].get('original_support') and complete[n].get('stage43_equality_reduced')]
    ratios = [cpu[n]['original_support'] / cpu[n]['stage43_equality_reduced'] for n in joint]
    out['tier1_joint_completions'] = dict(
        cases=len(joint),
        enumeration_median_cpu_seconds=median(cpu[n]['original_support'] for n in joint) if joint else None,
        structural_median_cpu_seconds=median(cpu[n]['stage43_equality_reduced'] for n in joint) if joint else None,
        median_casewise_ratio_enumeration_over_structural=median(ratios) if ratios else None,
        ratio_quantiles=quantiles(ratios),
        structural_faster_cases=sum(r > 1 for r in ratios),
        enumeration_faster_cases=sum(r < 1 for r in ratios),
        ties=sum(r == 1 for r in ratios),
        by_K={str(k): dict(cases=sum(meta[n]['K'] == k for n in joint),
                           median_ratio=median([cpu[n]['original_support'] / cpu[n]['stage43_equality_reduced']
                                                for n in joint if meta[n]['K'] == k]) or None)
              for k in sorted({meta[n]['K'] for n in joint})})

    # answer agreement across every completed tier-1 run of the same case
    disagreements = []
    for name, methods in by_case.items():
        answers = {(x['result'].get('feasible'), x['result'].get('minimum'))
                   for runs in methods.values() for x in runs if x['result']['status'] == 'complete'}
        if len(answers) > 1:
            disagreements.append(dict(name=name, answers=sorted(map(str, answers))))
    out['tier1_completed_answer_disagreements'] = disagreements
    minima = {n: sorted({x['result']['minimum'] for runs in by_case[n].values()
                         for x in runs if x['result']['status'] == 'complete'})
              for n in by_case}
    out['tier1_minimum_equals_generating_K'] = dict(
        cases_with_any_completion=sum(1 for n in minima if minima[n]),
        minimum_equals_K=sum(1 for n in minima if minima[n] and minima[n] == [meta[n]['K']]),
        minimum_below_K=sorted(n for n in minima if minima[n] and minima[n][0] is not None
                               and minima[n][0] < meta[n]['K']))

    # N- and K-dependence over completions only
    ndep, kdep = {}, {}
    for method in ('original_support', 'stage43_equality_reduced'):
        for k in sorted({meta[n]['K'] for n in by_case}):
            points = [(meta[n]['N'], cpu[n][method]) for n in by_case
                      if complete[n].get(method) and meta[n]['K'] == k and cpu[n][method] > 0]
            fit = ols([math.log10(x) for x, _ in points], [math.log10(y) for _, y in points])
            ndep[f'{method}:K{k}'] = dict(points=len(points), fit=fit,
                                          N_range=[min(x for x, _ in points), max(x for x, _ in points)]
                                          if points else None)
        for cell_prefix in sorted({meta[n]['cell'].split('_K')[0] for n in by_case}):
            points = [(meta[n]['K'], cpu[n][method]) for n in by_case
                      if complete[n].get(method) and meta[n]['cell'].startswith(cell_prefix + '_K')
                      and cpu[n][method] > 0]
            fit = ols([x for x, _ in points], [math.log10(y) for _, y in points])
            kdep[f'{method}:{cell_prefix}'] = dict(
                points=len(points), fit=fit,
                K_values=sorted({x for x, _ in points}),
                factor_per_unit_K=(10 ** fit['slope']) if fit else None)
    out['N_dependence_log10cpu_vs_log10N'] = ndep
    out['K_dependence_log10cpu_vs_K'] = kdep

    # witnesses
    witness = [r for r in records if 'independent_witness_verified' in r['result']]
    out['witness_verification'] = dict(
        runs_returning_a_witness=len(witness),
        verified=sum(bool(r['result']['independent_witness_verified']) for r in witness),
        failed=sum(not r['result']['independent_witness_verified'] for r in witness),
        count_matches_minimum=sum(bool(r['result'].get('witness_count_matches_minimum')) for r in witness),
        failures=[dict(name=r['name'], tier=r['tier'], method=r['method'],
                       detail=r['result'].get('independent_witness_verification'))
                  for r in witness if not r['result']['independent_witness_verified']][:20],
        by_method=dict(Counter(r['method'] for r in witness)))
    out['runs_with_no_witness_though_feasible'] = [
        dict(name=r['name'], tier=r['tier'], method=r['method'])
        for r in records if r['result'].get('feasible') and 'witness' not in r['result']][:20]

    # tier 2 and tier 3
    tier2 = [r for r in records if r['tier'] == 2]
    out['tier2_structural_long_budget'] = dict(
        jobs=len(tier2), status=dict(Counter(r['result']['status'] for r in tier2)),
        by_K={str(k): dict(jobs=sum(r['K'] == k for r in tier2),
                           complete=sum(r['K'] == k and r['result']['status'] == 'complete' for r in tier2),
                           median_cpu_seconds=median([float(r['result']['child_cpu_seconds'])
                                                      for r in tier2 if r['K'] == k
                                                      and r['result']['status'] == 'complete'] or [0]) or None)
              for k in sorted({r['K'] for r in tier2})})
    out['tier3_large_N_reduced_cap'] = [
        dict(name=r['name'], N=r['N'], cap=r['cap'], method=r['method'],
             status=r['result']['status'], reason=r['result'].get('reason'),
             feasible=r['result'].get('feasible'), minimum=r['result'].get('minimum'),
             cpu_seconds=float(r['result']['child_cpu_seconds']),
             rss_mib=float(r['result']['maximum_rss_kib']) / 1024)
        for r in records if r['tier'] == 3]

    done1 = [r for r in tier1 if r['result']['status'] == 'complete'
             and float(r['result']['process_wall_seconds']) > 0.05]
    out['concurrency_audit_cpu_over_wall'] = quantiles(
        [float(r['result']['child_cpu_seconds']) / float(r['result']['process_wall_seconds'])
         for r in done1])
    return out


# --------------------------------------------------------- incompleteness

def incompleteness_analysis():
    data = json.loads((ROOT / 'results/stage288_incompleteness.json').read_text())
    census = json.loads((ROOT / 'results/stage288_mark_census.json').read_text())
    control = json.loads((ROOT / 'results/stage288_mark_control.json').read_text())
    records = data['records']
    out = dict(jobs_total=len(records),
               status_by_variant={v: dict(Counter(r['result']['status'] for r in records
                                                  if r['variant'] == v))
                                  for v in sorted({r['variant'] for r in records})},
               mark_census_summary=census['summary'],
               positive_control=[{k: v for k, v in row.items() if k != 'B_marked_vertices'}
                                 for row in control['rows']])

    index = {}
    for r in records:
        index[(r['name'], r.get('level', 'none'), r.get('regime', 'none'),
               r.get('replicate', 0), r['variant'])] = r
    reference = {r['name']: r for r in records if r['variant'] == 'A_full'}
    out['A_full'] = dict(
        cases=len(reference),
        complete=sum(r['result']['status'] == 'complete' for r in reference.values()),
        feasible=sum(bool(r['result'].get('feasible')) for r in reference.values()),
        minimum_equals_K=sum(r['result'].get('minimum') == r['K'] for r in reference.values()),
        status=dict(Counter(r['result']['status'] for r in reference.values())))

    comparisons, by_regime, by_level, inflation = [], defaultdict(Counter), defaultdict(Counter), []
    for key, rec in index.items():
        if rec['variant'] != 'B_marks_retained':
            continue
        name, level, regime, replicate, _ = key
        a = reference.get(name)
        b = rec
        cvar = index.get((name, level, regime, replicate, 'C_naive_subset'))
        row = dict(name=name, N=rec['N'], K=rec['K'], cell=rec['cell'], level=level, regime=regime,
                   replicate=replicate, n_removed=rec['n_removed'], n_retained=rec['n_retained'],
                   n_marks=rec.get('n_marks', 0),
                   A_status=a['result']['status'] if a else None,
                   A_min=a['result'].get('minimum') if a else None,
                   B_status=b['result']['status'], B_min=b['result'].get('minimum'),
                   C_status=cvar['result']['status'] if cvar else None,
                   C_min=cvar['result'].get('minimum') if cvar else None,
                   C_discarded_marks=len(cvar['result'].get('discarded_marks') or []) if cvar else None)
        row['AB_comparable'] = row['A_status'] == 'complete' and row['B_status'] == 'complete'
        row['BC_comparable'] = row['B_status'] == 'complete' and row['C_status'] == 'complete'
        row['AC_comparable'] = row['A_status'] == 'complete' and row['C_status'] == 'complete'
        row['B_equals_A'] = row['AB_comparable'] and row['A_min'] == row['B_min']
        row['C_equals_B'] = row['BC_comparable'] and row['B_min'] == row['C_min']
        row['C_minus_A'] = (row['C_min'] - row['A_min']
                            if row['AC_comparable'] and row['A_min'] is not None
                            and row['C_min'] is not None else None)
        row['C_infeasible_while_A_feasible'] = bool(
            row['AC_comparable'] and row['A_min'] is not None and row['C_min'] is None)
        comparisons.append(row)
        if row['AB_comparable']:
            by_regime[regime]['AB_comparable'] += 1
            by_regime[regime]['B_equals_A'] += int(row['B_equals_A'])
            by_level[level]['AB_comparable'] += 1
            by_level[level]['B_equals_A'] += int(row['B_equals_A'])
        if row['AC_comparable']:
            by_regime[regime]['AC_comparable'] += 1
            by_regime[regime]['C_equals_A'] += int(row['C_minus_A'] == 0)
            by_level[level]['AC_comparable'] += 1
            by_level[level]['C_equals_A'] += int(row['C_minus_A'] == 0)
            if row['C_minus_A'] is not None:
                inflation.append(row['C_minus_A'])

    ab = [r for r in comparisons if r['AB_comparable'] and r['A_min'] is not None
          and r['B_min'] is not None]

    def undercount_table(keyfn, keys):
        table = {}
        for key in keys:
            subset = [r for r in ab if keyfn(r) == key]
            deltas = [r['B_min'] - r['A_min'] for r in subset]
            table[str(key)] = dict(
                comparable=len(subset),
                B_equals_A=sum(d == 0 for d in deltas),
                B_below_A=sum(d < 0 for d in deltas),
                B_above_A=sum(d > 0 for d in deltas),
                fraction_B_below_A=(sum(d < 0 for d in deltas) / len(deltas)) if deltas else None,
                mean_B_minus_A=(sum(deltas) / len(deltas)) if deltas else None,
                worst_B_minus_A=min(deltas) if deltas else None,
                distribution=dict(Counter(deltas)))
        return table

    out['undercount_B_versus_A'] = dict(
        comparable=len(ab),
        distribution=dict(Counter(r['B_min'] - r['A_min'] for r in ab)),
        B_below_A=sum(r['B_min'] < r['A_min'] for r in ab),
        B_above_A=sum(r['B_min'] > r['A_min'] for r in ab),
        quantiles=quantiles([r['B_min'] - r['A_min'] for r in ab]),
        by_regime=undercount_table(lambda r: r['regime'], sorted({r['regime'] for r in ab})),
        by_level=undercount_table(lambda r: r['level'], sorted({r['level'] for r in ab})),
        by_K=undercount_table(lambda r: r['K'], sorted({r['K'] for r in ab})),
        by_N_bin=undercount_table(lambda r: r['cell'].split('_K')[0],
                                  sorted({r['cell'].split('_K')[0] for r in ab})),
        by_removed_fraction_decile=undercount_table(
            lambda r: round(r['n_removed'] / (r['n_removed'] + r['n_retained']), 1),
            sorted({round(r['n_removed'] / (r['n_removed'] + r['n_retained']), 1) for r in ab})))
    out['four_way'] = dict(
        draws=len(comparisons),
        AB_comparable=sum(r['AB_comparable'] for r in comparisons),
        B_equals_A=sum(r['B_equals_A'] for r in comparisons),
        B_differs_from_A=[r for r in comparisons if r['AB_comparable'] and not r['B_equals_A']][:40],
        BC_comparable=sum(r['BC_comparable'] for r in comparisons),
        C_equals_B=sum(r['C_equals_B'] for r in comparisons),
        C_differs_from_B=[r for r in comparisons if r['BC_comparable'] and not r['C_equals_B']][:40],
        AC_comparable=sum(r['AC_comparable'] for r in comparisons),
        inflation_distribution=dict(Counter(inflation)),
        inflation_quantiles=quantiles(inflation),
        C_infeasible_while_A_feasible=sum(r['C_infeasible_while_A_feasible'] for r in comparisons),
        draws_with_non_empty_mark_set=sum(r['n_marks'] > 0 for r in comparisons),
        by_regime={k: dict(v) for k, v in by_regime.items()},
        by_level={k: dict(v) for k, v in by_level.items()},
        by_K={str(k): dict(draws=sum(r['K'] == k for r in comparisons),
                           AB_comparable=sum(r['K'] == k and r['AB_comparable'] for r in comparisons),
                           B_equals_A=sum(r['K'] == k and r['B_equals_A'] for r in comparisons),
                           AC_comparable=sum(r['K'] == k and r['AC_comparable'] for r in comparisons))
              for k in sorted({r['K'] for r in comparisons})},
        by_N_bin={b: dict(draws=sum(r['cell'].split('_K')[0] == b for r in comparisons),
                          AB_comparable=sum(r['cell'].split('_K')[0] == b and r['AB_comparable']
                                            for r in comparisons),
                          B_equals_A=sum(r['cell'].split('_K')[0] == b and r['B_equals_A']
                                         for r in comparisons))
                  for b in sorted({r['cell'].split('_K')[0] for r in comparisons})})

    for variant in ('D_prescribed_pendant', 'D1_fixed_one_missing'):
        subset = [r for r in records if r['variant'] == variant]
        out[variant] = dict(
            jobs=len(subset), status=dict(Counter(r['result']['status'] for r in subset)),
            complete=sum(r['result']['status'] == 'complete' for r in subset),
            feasible=sum(bool(r['result'].get('feasible')) for r in subset),
            infeasible_within_cap=sum(r['result']['status'] == 'complete'
                                      and not r['result'].get('feasible') for r in subset),
            by_level={l: dict(jobs=sum(r['level'] == l for r in subset),
                              complete=sum(r['level'] == l and r['result']['status'] == 'complete'
                                           for r in subset))
                      for l in sorted({r['level'] for r in subset})},
            by_regime={g: dict(jobs=sum(r['regime'] == g for r in subset),
                               complete=sum(r['regime'] == g and r['result']['status'] == 'complete'
                                            for r in subset))
                       for g in sorted({r['regime'] for r in subset})},
            by_K={str(k): dict(jobs=sum(r['K'] == k for r in subset),
                               complete=sum(r['K'] == k and r['result']['status'] == 'complete'
                                            for r in subset))
                  for k in sorted({r['K'] for r in subset})},
            by_N_bin={b: dict(jobs=sum(r['cell'].split('_K')[0] == b for r in subset),
                              complete=sum(r['cell'].split('_K')[0] == b
                                           and r['result']['status'] == 'complete' for r in subset))
                      for b in sorted({r['cell'].split('_K')[0] for r in subset})},
            median_cpu_seconds_completed=median([float(r['result']['child_cpu_seconds'])
                                                 for r in subset
                                                 if r['result']['status'] == 'complete'] or [0]) or None)
    both = [(r, index.get((r['name'], r['level'], r['regime'], r['replicate'], 'D1_fixed_one_missing')))
            for r in records if r['variant'] == 'D_prescribed_pendant' and r['level'] == 'one']
    pairs = [(a, b) for a, b in both if b is not None]
    out['D_vs_D1_one_missing'] = dict(
        pairs=len(pairs),
        both_complete=sum(a['result']['status'] == 'complete' and b['result']['status'] == 'complete'
                          for a, b in pairs),
        agree=sum(a['result'].get('minimum') == b['result'].get('minimum') for a, b in pairs
                  if a['result']['status'] == 'complete' and b['result']['status'] == 'complete'),
        disagreements=[dict(name=a['name'], regime=a['regime'], D=a['result'].get('minimum'),
                            D1=b['result'].get('minimum')) for a, b in pairs
                       if a['result']['status'] == 'complete' and b['result']['status'] == 'complete'
                       and a['result'].get('minimum') != b['result'].get('minimum')][:20])

    witness = [r for r in records if 'independent_witness_verified' in r['result']]
    out['witness_verification'] = dict(
        runs_returning_a_witness=len(witness),
        verified=sum(bool(r['result']['independent_witness_verified']) for r in witness),
        failed=sum(not r['result']['independent_witness_verified'] for r in witness),
        by_variant=dict(Counter(r['variant'] for r in witness)),
        failures=[dict(name=r['name'], variant=r['variant'], level=r.get('level'),
                       regime=r.get('regime'),
                       detail=r['result'].get('independent_witness_verification'))
                  for r in witness if not r['result']['independent_witness_verified']][:20],
        prescribed_pendant_respected=sum(bool(r['result'].get('prescribed_pendant_respected'))
                                         for r in witness if r['variant'].startswith('D')))
    out['_comparison_rows'] = comparisons
    return out


def main():
    result = dict(stage=288, created_utc=__import__('datetime').datetime.now(
        __import__('datetime').timezone.utc).isoformat())
    result['scaling'] = scaling_analysis()
    result['incompleteness'] = incompleteness_analysis()
    result['protocol_hashes'] = {
        name: sha(ROOT / 'results' / name) for name in
        ['stage288_scaling_protocol.json', 'stage288_incompleteness_protocol.json']
        if (ROOT / 'results' / name).exists()}
    OUT.write_text(json.dumps(serial(result), indent=2) + '\n')
    printable = json.loads(json.dumps(serial(result)))
    printable['incompleteness'].pop('_comparison_rows', None)
    print(json.dumps(printable, indent=2)[:12000])
    print('WROTE', OUT)


if __name__ == '__main__':
    main()
