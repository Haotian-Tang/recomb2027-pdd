#!/usr/bin/env python3
"""Stage 289: build the comparison tables from the measured summary and per-case records."""
import json
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median

ROOT = Path(__file__).resolve().parents[1]
PERCASE = ROOT / 'results/stage289/stage289_headtohead_percase.json'
SUMMARY = ROOT / 'results/stage289/stage289_headtohead_summary.json'
TABLES = ROOT / 'results/stage289/stage289_tables.md'
DD = ('dd_default', 'dd_true_k', 'dd_resolution')
EXACT = 'exact_stage43_partial'


def fmt(value, digits=3):
    if value is None:
        return 'n/a'
    if isinstance(value, float):
        if value == 0:
            return '0'
        if abs(value) < 1e-3 or abs(value) >= 1e5:
            return f'{value:.{digits}g}'
        return f'{value:.{digits}f}'
    return str(value)


def table(rows, header):
    widths = [max(len(str(r[i])) for r in [header] + rows) for i in range(len(header))]
    out = ['| ' + ' | '.join(str(h).ljust(w) for h, w in zip(header, widths)) + ' |',
           '|' + '|'.join('-' * (w + 2) for w in widths) + '|']
    for row in rows:
        out.append('| ' + ' | '.join(str(v).ljust(w) for v, w in zip(row, widths)) + ' |')
    return '\n'.join(out)


def minimal_bound(case):
    """Upper bound on the true minimum: the generator used k_true sources."""
    exact = case['methods'][EXACT]
    if exact['status'] == 'complete' and exact.get('minimum') is not None:
        return exact['minimum'], 'exact'
    if case['exactly_explained']:
        return case['k_true'], 'generator'
    return None, None


def main():
    percase = json.loads(PERCASE.read_text())
    summary = json.loads(SUMMARY.read_text())
    lines = ['# Stage 289 comparison tables', '',
             f'Cases analysed: {summary["cases_analysed"]}; '
             f'cases with incomplete run records: {summary["cases_with_incomplete_run_records"]}.', '']

    lines += ['## 1. Completion and runtime (median child CPU seconds, 2 repetitions per case)', '']
    rows = []
    for family, block in summary['families'].items():
        for method in (EXACT,) + DD:
            entry = block['by_method'][method]
            rows.append([family, method, entry['runs'],
                         entry['status_counts'].get('complete', 0),
                         entry['status_counts'].get('unknown', 0),
                         entry['status_counts'].get('error', 0),
                         fmt(entry['median_cpu_seconds_completed']),
                         fmt(entry['total_cpu_seconds_all_runs'], 4)])
    lines += [table(rows, ['family', 'method', 'cases', 'complete', 'timeout', 'error',
                           'median CPU s', 'total CPU s']), '']

    lines += ['## 2. Source count returned, against the generating count', '']
    rows = []
    for family, block in summary['families'].items():
        for method in DD:
            entry = block['by_method'][method]
            total = entry['status_counts'].get('complete', 0)
            rows.append([family, method, total, entry['k_equals_k_true'],
                         entry['k_below_k_true'], entry['k_above_k_true'],
                         json.dumps(entry['stop_reasons'])])
    lines += [table(rows, ['family', 'arm', 'complete', 'k = k_true', 'k < k_true',
                           'k > k_true', 'stop reason']), '']

    lines += ['## 3. Support recovery and exact feasibility of the returned support', '']
    rows = []
    for family, block in summary['families'].items():
        for method in DD:
            entry = block['by_method'][method]
            total = entry['status_counts'].get('complete', 0)
            feas = entry['exact_feasibility']
            rows.append([family, method, total, entry['support_exact_match'],
                         fmt(entry['recovered_fraction_of_true_edges']),
                         feas.get('feasible', 0), feas.get('infeasible', 0),
                         feas.get('unresolved', 0) + feas.get('invalid', 0),
                         fmt(entry['median_w1_to_truth']), fmt(entry['median_loss'])])
        entry = block['by_method'][EXACT]
        rows.append([family, EXACT, entry['status_counts'].get('complete', 0),
                     entry['support_exact_match'], '-', entry['feasible'], entry['infeasible'],
                     '-', '-', '-'])
    lines += [table(rows, ['family', 'method', 'complete', 'support = truth',
                           'true edges recovered', 'support exactly feasible',
                           'support provably infeasible', 'other', 'median W1', 'median loss']), '']

    lines += ['## 4. By generating source count', '']
    for family, block in summary['families'].items():
        rows = []
        for k, entry in sorted(block['by_k'].items(), key=lambda kv: int(kv[0])):
            rows.append([k, entry['cases'], entry['exact_complete'], fmt(entry['exact_median_cpu']),
                         *[entry[f'{m}_support_exact_match'] for m in DD],
                         *[entry[f'{m}_exact_feasible'] for m in DD],
                         *[fmt(entry[f'{m}_median_cpu']) for m in DD]])
        lines += [f'### {family}', '',
                  table(rows, ['k_true', 'cases', 'exact complete', 'exact median CPU s',
                               *[f'{m} support=truth' for m in DD],
                               *[f'{m} feasible' for m in DD],
                               *[f'{m} median CPU s' for m in DD]]), '']

    lines += ['## 5. By reference-set size', '']
    for family, block in summary['families'].items():
        rows = []
        order = {'<=50': 0, '51-150': 1, '151-500': 2}
        for key, entry in sorted(block['by_reference_size'].items(), key=lambda kv: order.get(kv[0], 9)):
            rows.append([key, entry['cases'], entry['exact_complete'], fmt(entry['exact_median_cpu']),
                         *[entry[f'{m}_support_exact_match'] for m in DD],
                         *[fmt(entry[f'{m}_median_cpu']) for m in DD]])
        lines += [f'### {family}', '',
                  table(rows, ['leaves', 'cases', 'exact complete', 'exact median CPU s',
                               *[f'{m} support=truth' for m in DD],
                               *[f'{m} median CPU s' for m in DD]]), '']

    lines += ['## 6. Minimality on exactly explainable data', '',
              'An explanation is counted only when the returned support is certified exactly '
              'feasible. The reference minimum is the exact solver\'s answer where it completed, '
              'otherwise the generating count, which is an upper bound on the true minimum.', '']
    rows = []
    for family in ('exact_10k', 'exact_biotrees'):
        subset = [c for c in percase if c['family'] == family]
        if not subset:
            continue
        for method in DD:
            counts = Counter()
            for case in subset:
                record = case['methods'][method]
                if record['status'] != 'complete':
                    counts['not complete'] += 1
                    continue
                bound, _ = minimal_bound(case)
                if record.get('exact_feasibility') != 'feasible':
                    counts['no exact explanation returned'] += 1
                elif bound is None:
                    counts['no reference minimum'] += 1
                elif record['k'] > bound:
                    counts['non-minimal explanation'] += 1
                elif record['k'] < bound:
                    counts['below the reference minimum'] += 1
                else:
                    counts['minimal explanation'] += 1
            rows.append([family, method, len(subset),
                         counts['minimal explanation'], counts['non-minimal explanation'],
                         counts['below the reference minimum'],
                         counts['no exact explanation returned'], counts['not complete']])
    lines += [table(rows, ['family', 'arm', 'cases', 'exact and minimal', 'exact but non-minimal',
                           'below reference minimum', 'no exact explanation', 'incomplete']), '']

    lines += ['## 7. Identifiability of the returned support', '']
    rows = []
    for family, _ in summary['families'].items():
        subset = [c for c in percase if c['family'] == family]
        for method in DD:
            values = [c['methods'][method] for c in subset if c['methods'][method]['status'] == 'complete']
            rows.append([family, method, len(values),
                         sum(r.get('exact_feasibility') == 'feasible' and r.get('identifiable') is False
                             for r in values),
                         sum(bool(r.get('boundary_placements')) for r in values),
                         sum(bool(r.get('repeated_edge')) for r in values)])
    lines += [table(rows, ['family', 'arm', 'complete', 'feasible but non-identifiable',
                           'placement at a vertex', 'repeated edge']), '']

    TABLES.write_text('\n'.join(lines) + '\n')
    print('\n'.join(lines)[:3000])


if __name__ == '__main__':
    main()
