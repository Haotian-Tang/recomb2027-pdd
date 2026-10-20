#!/usr/bin/env python3
"""Stage 293 analysis: abundance and position error under reference loss.

W1 REUSE.  The mixture-level error is the Wasserstein-1 distance between the
recovered and the generating mixture under the tree metric, computed by calling
stage289_headtohead_analyze.wasserstein on a cost matrix built by
stage289_headtohead_analyze.point_distances -- the same two functions stage 289
used, imported, not reimplemented, so the numbers are directly comparable.

Reported quantities per draw
  support_exact_match            the returned minimum-size support is the true one
  max_p_abs_error / tv_p         |p - p_true|_inf and total variation, exact
  max_x_abs_error                |x - x_true|_inf, exact, only when support matches
  w1_to_truth                    tree-metric W1, stage 289 definition
  fibre free_parameters          dim of the set of explanations on a fixed support
  p_sup_range / x_sup_range      the fibre's l_inf diameter in p and in x, exact
  max_p_deviation_from_truth     largest p error still exactly explaining the data
"""
import os
for _k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_k] = '1'
import json
import math
import statistics as stats
import sys
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import stage288_common as common
import stage288_incompleteness as inc
from stage289_headtohead_analyze import wasserstein, point_distances, float_adjacency, distances_from

OUTDIR = ROOT / 'results/stage293'
LEVELS = ['none'] + list(inc.LEVELS)
REGIMES = ['random', 'clade', 'nearest']


def load_stream(name):
    path = OUTDIR / f'stage293_{name}_runs.jsonl'
    if not path.exists():
        return {}
    out = {}
    with open(path) as handle:
        for line in handle:
            row = json.loads(line)
            key = (row['name'], row['level'], row['regime'], row['replicate'])
            if name == 'decodiphy':
                key = key + (row['arm'],)
            out[key] = row
    return out


def quantiles(values):
    values = sorted(values)
    if not values:
        return None
    def q(f):
        if len(values) == 1:
            return values[0]
        position = f * (len(values) - 1)
        low = int(math.floor(position))
        high = min(low + 1, len(values) - 1)
        return values[low] + (values[high] - values[low]) * (position - low)
    return dict(n=len(values), min=values[0], q25=q(0.25), median=q(0.5),
                q75=q(0.75), q90=q(0.9), max=values[-1],
                mean=sum(values) / len(values))


def describe(values, digits=4):
    summary = quantiles(values)
    if summary is None:
        return None
    return {k: (round(v, digits) if isinstance(v, float) else v) for k, v in summary.items()}


def main():
    exact = load_stream('exact')
    fibre_true = load_stream('fibre_true')
    fibre_found = load_stream('fibre_found')
    dd = load_stream('decodiphy')
    print('records: exact=%d fibre_true=%d fibre_found=%d decodiphy=%d'
          % (len(exact), len(fibre_true), len(fibre_found), len(dd)), flush=True)

    cases = {}
    for entry in inc.panel():
        cases[entry['name']] = common.load_case(entry)

    rows = []
    for key, row in sorted(exact.items()):
        name, level, regime, replicate = key
        case = cases[name]
        result = row['result']
        truth = case['generating_sources']
        true_support = sorted(int(s['edge']) for s in truth)
        record = dict(name=name, family=row['family'], cell=row['cell'], level=level,
                      regime=regime, replicate=replicate,
                      n_leaves=row['n_leaves'], n_removed=row['n_removed'],
                      n_retained=row['n_retained'],
                      removed_fraction=row['removed_fraction'],
                      k_true=row['K_generating'], cap=row['cap'],
                      exact_status=result['status'], minimum=result.get('minimum'),
                      exact_wall=result.get('process_wall_seconds'),
                      support=result.get('support'),
                      support_exact_match=result.get('support_exact_match'),
                      n_true_edges_recovered=result.get('n_true_edges_recovered'),
                      witness_verified=result.get('independent_witness_verified'),
                      explains_removed=result.get('explains_removed_references'))
        if result.get('minimum') is not None:
            record['count_deficit'] = row['K_generating'] - result['minimum']
        for field in ('max_p_abs_error', 'max_x_abs_error', 'total_variation_p',
                      'pendant_abs_error'):
            if field in result:
                record[field] = float(F(result[field]))
        # ---- fibre over the true support
        ft = fibre_true.get(key)
        if ft:
            f = ft['result']
            record['fibre_true_status'] = f['status']
            record['true_support_feasible'] = (f['status'] == 'ok')
            if f['status'] == 'ok':
                record['fibre_true_free_parameters'] = f['free_parameters']
                record['fibre_true_p_sup_range'] = float(F(f['p_sup_range']))
                record['fibre_true_y_sup_range'] = float(F(f['y_sup_range']))
                if 'x_sup_range' in f:
                    record['fibre_true_x_sup_range'] = float(F(f['x_sup_range']))
                if 'max_p_deviation_from_truth' in f:
                    record['fibre_true_max_p_dev'] = float(F(f['max_p_deviation_from_truth']))
                if 'max_x_deviation_from_truth' in f:
                    record['fibre_true_max_x_dev'] = float(F(f['max_x_deviation_from_truth']))
                record['fibre_true_truth_inside'] = f.get('truth_inside_bounds')
        ff = fibre_found.get(key)
        if ff and ff['result']['status'] == 'ok':
            f = ff['result']
            record['fibre_found_free_parameters'] = f['free_parameters']
            record['fibre_found_p_sup_range'] = float(F(f['p_sup_range']))
            if 'x_sup_range' in f:
                record['fibre_found_x_sup_range'] = float(F(f['x_sup_range']))
        rows.append(record)

    # ---------------- W1, reusing the stage 289 functions verbatim
    print('computing W1 (stage289 wasserstein + point_distances)...', flush=True)
    by_case = defaultdict(list)
    for record in rows:
        by_case[record['name']].append(record)
    for name, records in by_case.items():
        case = cases[name]
        # exactly the stage 289 input convention: Fraction lengths, floated inside
        edges = case['edges']
        adjacency = float_adjacency(edges)
        truth = [dict(edge=s['edge'], p=float(F(s['p'])), x=float(F(s['x'])))
                 for s in case['generating_sources']]
        needed = {v for s in truth for v in edges[s['edge']][:2]}
        wanted = []
        for record in records:
            key = (record['name'], record['level'], record['regime'], record['replicate'])
            witness = exact[key]['result'].get('witness')
            if witness:
                sources = [dict(edge=int(s['edge']), p=float(F(s['p'])), x=float(F(s['x'])))
                           for s in witness['sources']]
                wanted.append((record, 'exact', sources))
                needed.update(v for s in sources for v in edges[s['edge']][:2])
            for arm in ('dd_default', 'dd_true_k'):
                row = dd.get(key + (arm,))
                if row and row['result'].get('status') == 'complete':
                    sources = [dict(edge=int(s['edge']), p=float(s['p']), x=float(s['x']))
                               for s in row['result']['sources']]
                    wanted.append((record, arm, sources))
                    needed.update(v for s in sources for v in edges[s['edge']][:2])
        table = {v: distances_from(adjacency, v) for v in needed}
        for record, arm, sources in wanted:
            cost = point_distances(edges, adjacency, table, sources, truth)
            value = wasserstein(cost, [s['p'] for s in sources], [s['p'] for s in truth])
            record[f'w1_{arm}'] = value
            record[f'placement_error_{arm}'] = sum(
                t['p'] * min(cost[i][j] for i in range(len(sources)))
                for j, t in enumerate(truth))
            if arm != 'exact':
                row = dd[(record['name'], record['level'], record['regime'],
                          record['replicate'], arm)]['result']
                record[f'k_{arm}'] = row['k']
                record[f'loss_{arm}'] = row['loss']
                record[f'support_match_{arm}'] = (
                    sorted({int(s['edge']) for s in row['sources']})
                    == sorted({int(s['edge']) for s in case['generating_sources']}))
                record[f'crash_{arm}'] = bool(row.get('jplace_writer_crashed'))
                # exact p error when the support matches and k == k_true
                if record[f'support_match_{arm}']:
                    got = {int(s['edge']): float(s['p']) for s in row['sources']}
                    tru = {int(s['edge']): float(F(s['p'])) for s in case['generating_sources']}
                    if len(got) == len(tru):
                        record[f'max_p_abs_error_{arm}'] = max(abs(got[e] - tru[e]) for e in tru)
    # DecoDiPhy runs that produced no usable answer
    for record in rows:
        key = (record['name'], record['level'], record['regime'], record['replicate'])
        for arm in ('dd_default', 'dd_true_k'):
            row = dd.get(key + (arm,))
            record[f'status_{arm}'] = row['result']['status'] if row else 'missing'
            if row:
                record[f'crash_{arm}'] = bool(row['result'].get('jplace_writer_crashed'))

    (OUTDIR / 'stage293_records.jsonl').write_text(
        '\n'.join(json.dumps(r) for r in rows) + '\n')
    print('wrote %d analysed records' % len(rows), flush=True)
    summary = summarise(rows)
    (OUTDIR / 'stage293_analysis.json').write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary['headline'], indent=1))


def summarise(rows):
    complete = [r for r in rows if r['exact_status'] == 'complete']
    matched = [r for r in complete if r.get('support_exact_match')]
    out = dict(n_records=len(rows), n_exact_complete=len(complete),
               n_support_match=len(matched))
    # verification
    verified = [r for r in complete if r.get('witness_verified') is not None]
    out['verification'] = dict(
        witnesses_checked=len(verified),
        all_witnesses_verified=all(r['witness_verified'] for r in verified),
        truth_inside_fibre_checked=sum(1 for r in rows if r.get('fibre_true_truth_inside') is not None),
        truth_inside_fibre_all=all(r.get('fibre_true_truth_inside') is not False for r in rows),
        witnesses_explaining_removed=sum(1 for r in complete if r.get('explains_removed') is True),
        witnesses_with_removed_checked=sum(1 for r in complete if r.get('explains_removed') is not None))
    # Invariant 1: the generating mixture explains every subset of the observations,
    # so the minimum over the retained set can never exceed K_true.
    over = [dict(name=r['name'], level=r['level'], regime=r['regime'],
                 minimum=r['minimum'], k_true=r['k_true'])
            for r in complete if r.get('minimum') is not None and r['minimum'] > r['k_true']]
    # Invariant 2: for the same reason the TRUE support must stay exactly feasible
    # on every retained set; an 'infeasible' verdict would be a bug, not a finding.
    infeasible = [dict(name=r['name'], level=r['level'], regime=r['regime'],
                       status=r['fibre_true_status'])
                  for r in rows if r.get('fibre_true_status') not in (None, 'ok')]
    out['invariants'] = dict(
        minimum_never_exceeds_k_true=not over, n_violations_minimum=len(over),
        violations_minimum=over[:10],
        true_support_always_feasible=not infeasible,
        n_true_support_not_ok=len(infeasible), true_support_not_ok=infeasible[:10],
        n_fibre_true_evaluated=sum(1 for r in rows if r.get('fibre_true_status') is not None))

    def cut(rows_in, **filters):
        out_rows = rows_in
        for field, value in filters.items():
            out_rows = [r for r in out_rows if r.get(field) == value]
        return out_rows

    # by level, by regime
    grid = {}
    for level in LEVELS:
        for regime in (['none'] if level == 'none' else REGIMES):
            cell = cut(complete, level=level, regime=regime)
            if not cell:
                continue
            cell_matched = [r for r in cell if r.get('support_exact_match')]
            entry = dict(
                n=len(cell), n_support_match=len(cell_matched),
                support_match_rate=len(cell_matched) / len(cell),
                median_removed_fraction=round(stats.median(
                    r['removed_fraction'] for r in cell), 4),
                count_deficit_rate=(sum(1 for r in cell if (r.get('count_deficit') or 0) > 0)
                                    / len(cell)),
                w1_exact=describe([r['w1_exact'] for r in cell if r.get('w1_exact') is not None]),
                w1_exact_support_match=describe(
                    [r['w1_exact'] for r in cell_matched if r.get('w1_exact') is not None]),
                max_p_abs_error_support_match=describe(
                    [r['max_p_abs_error'] for r in cell_matched if 'max_p_abs_error' in r]),
                max_x_abs_error_support_match=describe(
                    [r['max_x_abs_error'] for r in cell_matched if 'max_x_abs_error' in r]),
                tv_p_support_match=describe(
                    [r['total_variation_p'] for r in cell_matched if 'total_variation_p' in r]),
                true_support_feasible_rate=(
                    sum(1 for r in cell if r.get('true_support_feasible')) / len(cell)),
                fibre_true_positive_dim_rate=(
                    sum(1 for r in cell if (r.get('fibre_true_free_parameters') or 0) > 0)
                    / max(1, sum(1 for r in cell if r.get('fibre_true_free_parameters') is not None))),
                fibre_true_free_parameters=describe(
                    [r['fibre_true_free_parameters'] for r in cell
                     if r.get('fibre_true_free_parameters') is not None], 2),
                fibre_true_p_sup_range=describe(
                    [r['fibre_true_p_sup_range'] for r in cell
                     if r.get('fibre_true_p_sup_range') is not None]),
                fibre_true_x_sup_range=describe(
                    [r['fibre_true_x_sup_range'] for r in cell
                     if r.get('fibre_true_x_sup_range') is not None]),
                fibre_true_max_p_dev=describe(
                    [r['fibre_true_max_p_dev'] for r in cell
                     if r.get('fibre_true_max_p_dev') is not None]),
                w1_dd_default=describe([r['w1_dd_default'] for r in cell
                                        if r.get('w1_dd_default') is not None]),
                w1_dd_true_k=describe([r['w1_dd_true_k'] for r in cell
                                       if r.get('w1_dd_true_k') is not None]),
                dd_default_support_match_rate=(
                    sum(1 for r in cell if r.get('support_match_dd_default')) /
                    max(1, sum(1 for r in cell if r.get('status_dd_default') == 'complete'))),
                dd_true_k_support_match_rate=(
                    sum(1 for r in cell if r.get('support_match_dd_true_k')) /
                    max(1, sum(1 for r in cell if r.get('status_dd_true_k') == 'complete'))),
                dd_default_complete=sum(1 for r in cell if r.get('status_dd_default') == 'complete'),
                dd_true_k_complete=sum(1 for r in cell if r.get('status_dd_true_k') == 'complete'),
                dd_default_crash=sum(1 for r in cell if r.get('crash_dd_default')),
                dd_true_k_crash=sum(1 for r in cell if r.get('crash_dd_true_k')))
            grid[f'{level}|{regime}'] = entry
    out['grid'] = grid

    # positive-dimensional fibres: the diameters
    nonident = [r for r in rows if (r.get('fibre_true_free_parameters') or 0) > 0]
    out['non_identifiable'] = dict(
        n=len(nonident),
        of_total_with_feasible_true_support=sum(
            1 for r in rows if r.get('fibre_true_free_parameters') is not None),
        free_parameters=describe([r['fibre_true_free_parameters'] for r in nonident], 2),
        p_sup_range=describe([r['fibre_true_p_sup_range'] for r in nonident]),
        x_sup_range=describe([r['fibre_true_x_sup_range'] for r in nonident
                              if r.get('fibre_true_x_sup_range') is not None]),
        max_p_dev=describe([r['fibre_true_max_p_dev'] for r in nonident
                            if r.get('fibre_true_max_p_dev') is not None]),
        by_level={level: describe([r['fibre_true_p_sup_range'] for r in nonident
                                   if r['level'] == level])
                  for level in LEVELS},
        at_full_panel=describe([r['fibre_true_p_sup_range'] for r in nonident
                                if r['level'] == 'none']))

    # does a CORRECT support still give wrong abundances?
    out['correct_support_still_wrong'] = {}
    for level in LEVELS:
        cell = [r for r in rows if r['level'] == level and r.get('support_exact_match')
                and 'max_p_abs_error' in r]
        identifiable = [r for r in cell if r.get('fibre_true_free_parameters') == 0]
        loose = [r for r in cell if (r.get('fibre_true_free_parameters') or 0) > 0]
        out['correct_support_still_wrong'][level] = dict(
            n=len(cell),
            n_identifiable=len(identifiable), n_positive_dimensional=len(loose),
            p_error_identifiable=describe([r['max_p_abs_error'] for r in identifiable], 12),
            p_error_positive_dimensional=describe([r['max_p_abs_error'] for r in loose]),
            worst_case_p_error_positive_dimensional=describe(
                [r['fibre_true_max_p_dev'] for r in loose
                 if r.get('fibre_true_max_p_dev') is not None]))

    # exact vs heuristic degradation
    out['exact_vs_heuristic'] = {}
    for level in LEVELS:
        cell = [r for r in rows if r['level'] == level]
        both = [r for r in cell if r.get('w1_exact') is not None
                and r.get('w1_dd_default') is not None]
        out['exact_vs_heuristic'][level] = dict(
            n_both=len(both),
            w1_exact=describe([r['w1_exact'] for r in both]),
            w1_dd_default=describe([r['w1_dd_default'] for r in both]),
            w1_dd_true_k=describe([r['w1_dd_true_k'] for r in cell
                                   if r.get('w1_dd_true_k') is not None]),
            exact_better=sum(1 for r in both if r['w1_exact'] < r['w1_dd_default']),
            heuristic_better=sum(1 for r in both if r['w1_dd_default'] < r['w1_exact']),
            tied=sum(1 for r in both if r['w1_dd_default'] == r['w1_exact']))

    out['headline'] = dict(
        records=len(rows), exact_complete=len(complete),
        support_match_overall=len(matched),
        verification_all_witnesses_verified=out['verification']['all_witnesses_verified'],
        verification_truth_inside_fibre_all=out['verification']['truth_inside_fibre_all'],
        invariant_minimum_le_k_true=out['invariants']['minimum_never_exceeds_k_true'],
        invariant_true_support_feasible=out['invariants']['true_support_always_feasible'],
        non_identifiable_true_support=out['non_identifiable']['n'],
        non_identifiable_median_p_range=(out['non_identifiable']['p_sup_range'] or {}).get('median'))
    return out


if __name__ == '__main__':
    main()
