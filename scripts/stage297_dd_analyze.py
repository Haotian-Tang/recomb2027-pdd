#!/usr/bin/env python3
"""Stage297 analysis of the DecoDiPhy seed sweep.

Per case: the spread across seeds in the inferred k, in the returned support and
in W1 to truth, plus the fixed-seed control (arm default_1142, the stage 293
worker imported unmodified, against arm seeded_1142, the same seed through the
one-flag copy).  W1 is computed here with
`stage289_headtohead_analyze.wasserstein` and `point_distances`, the same
tree-metric transport distance stages 289/293/294 report.
"""
import json
import sys
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
sys.dont_write_bytecode = True

from fractions import Fraction as F                                        # noqa: E402
import stage297_common as common                                          # noqa: E402
from stage292_analysis import describe, spearman                          # noqa: E402
from stage289_headtohead_analyze import (float_adjacency, distances_from,  # noqa: E402
                                         point_distances, wasserstein)

OUTDIR = ROOT / 'results/stage297'
OUT = OUTDIR / 'stage297_dd_analysis.json'


def load_runs():
    rows = []
    for p in sorted(OUTDIR.glob('stage297_dd_seeds_*.jsonl')):
        rows.extend(json.loads(l) for l in p.open())
    return rows


def w1_for_case(entry, runs):
    """W1 to truth for every complete run of one case, in one distance table."""
    case = common.load_verified(entry)
    edges = case['edges']
    adjacency = float_adjacency(edges)
    truth = [dict(edge=int(s['edge']), p=float(F(s['p'])), x=float(F(s['x'])))
             for s in case['generating_sources']]
    needed = {v for s in truth for v in edges[s['edge']][:2]}
    for r in runs:
        for s in (r.get('sources') or []):
            needed.update(edges[int(s['edge'])][:2])
    table = {v: distances_from(adjacency, v) for v in needed}
    out = {}
    for r in runs:
        if r.get('status') != 'complete' or not r.get('sources'):
            out[r['arm']] = None
            continue
        sources = [dict(edge=int(s['edge']), p=float(F(s['p'])), x=float(F(s['x'])))
                   for s in r['sources']]
        cost = point_distances(edges, adjacency, table, sources, truth)
        out[r['arm']] = dict(
            w1_to_truth=wasserstein(cost, [s['p'] for s in sources],
                                    [s['p'] for s in truth]),
            weighted_placement_error=sum(
                t['p'] * min(cost[i][j] for i in range(len(sources)))
                for j, t in enumerate(truth)))
    return out


def jaccard(a, b):
    A, B = set(a), set(b)
    return len(A & B) / len(A | B) if (A or B) else 1.0


def per_case(rows, index):
    cases = defaultdict(list)
    for r in rows:
        if r.get('kind') != 'dd_seed':
            continue
        cases[r['name']].append(r)
    out = []
    for name, runs in sorted(cases.items()):
        entry = index[name]
        w1 = w1_for_case(entry, runs)
        by_arm = {r['arm']: r for r in runs}
        for arm, value in w1.items():
            if value:
                by_arm[arm].update(value)
        complete = [r for r in runs if r.get('status') == 'complete']
        seeded = [r for r in runs if r['arm'].startswith('seeded_')]
        seeded_complete = [r for r in seeded if r.get('status') == 'complete']
        ks = [r['k'] for r in seeded_complete]
        w1s = [r.get('w1_to_truth') for r in seeded_complete
               if r.get('w1_to_truth') is not None]
        sups = [tuple(sorted(r['support'])) for r in seeded_complete if r['support']]
        default = by_arm.get('default_1142')
        control = by_arm.get('seeded_1142')
        row = dict(
            name=name, family=entry['family'], K_generating=entry['K_generating'],
            R_size=entry['reference_set_size'], tree_leaves=entry['tree_leaves'],
            distance_source=entry['distance_source'],
            exactly_explained=entry['exactly_explained'],
            n_runs=len(runs), n_complete=len(complete),
            statuses=dict(Counter(r.get('status') for r in runs)),
            writer_crashed=sum(bool(r.get('jplace_writer_crashed')) for r in runs),
            seeds_complete=len(seeded_complete),
            k_values={r['arm']: r.get('k') for r in sorted(runs, key=lambda z: z['arm'])},
            k_distinct=len(set(ks)), k_min=min(ks) if ks else None,
            k_max=max(ks) if ks else None,
            k_range=(max(ks) - min(ks)) if ks else None,
            k_true=entry['K_generating'],
            support_distinct=len(set(sups)),
            support_sets=[list(s) for s in sorted(set(sups))],
            mean_pairwise_jaccard=(None if len(sups) < 2 else
                                   sum(jaccard(a, b) for a, b in combinations(sups, 2)) /
                                   len(list(combinations(sups, 2)))),
            w1_values={r['arm']: r.get('w1_to_truth')
                       for r in sorted(runs, key=lambda z: z['arm'])},
            w1_min=min(w1s) if w1s else None, w1_max=max(w1s) if w1s else None,
            w1_median=describe(w1s)['median'] if w1s else None,
            w1_range=(max(w1s) - min(w1s)) if w1s else None,
            w1_relative_range=((max(w1s) - min(w1s)) / describe(w1s)['median']
                               if w1s and describe(w1s)['median'] else None),
            loss_values={r['arm']: r.get('loss')
                         for r in sorted(runs, key=lambda z: z['arm'])},
            support_exact_match_any=any(r.get('support_exact_match')
                                        for r in seeded_complete),
            support_exact_match_all=(all(r.get('support_exact_match')
                                         for r in seeded_complete)
                                     if seeded_complete else None),
            n_true_edges_recovered_range=(
                (max(r['n_true_edges_recovered'] for r in seeded_complete) -
                 min(r['n_true_edges_recovered'] for r in seeded_complete))
                if seeded_complete else None),
            seconds=describe([r['solver_wall_seconds'] for r in runs
                              if r.get('solver_wall_seconds') is not None]))
        if default and control:
            both = (default.get('status') == 'complete' and
                    control.get('status') == 'complete')
            row['fixed_seed_control'] = dict(
                comparable=bool(both),
                k_equal=bool(both and default['k'] == control['k']),
                support_equal=bool(both and default['support'] == control['support']),
                loss_equal=bool(both and default['loss'] == control['loss']),
                w1_equal=bool(both and default.get('w1_to_truth') ==
                              control.get('w1_to_truth')),
                default_k=default.get('k'), control_k=control.get('k'),
                default_support=default.get('support'),
                control_support=control.get('support'))
        out.append(row)
    return out


def aggregate(cases):
    def block(rows, label):
        usable = [r for r in rows if r['seeds_complete'] >= 2]
        ctrl = [r['fixed_seed_control'] for r in rows
                if r.get('fixed_seed_control', {}).get('comparable')]
        return dict(
            label=label, cases=len(rows), cases_with_two_or_more_seeds=len(usable),
            availability=describe([r['n_complete'] / r['n_runs'] for r in rows]),
            k_identical_across_seeds=sum(1 for r in usable if r['k_distinct'] == 1),
            k_distinct_distribution=dict(Counter(r['k_distinct'] for r in usable)),
            k_range=describe([r['k_range'] for r in usable if r['k_range'] is not None]),
            support_identical_across_seeds=sum(1 for r in usable if r['support_distinct'] == 1),
            support_distinct_distribution=dict(Counter(r['support_distinct'] for r in usable)),
            mean_pairwise_jaccard=describe([r['mean_pairwise_jaccard'] for r in usable
                                            if r['mean_pairwise_jaccard'] is not None]),
            w1_median=describe([r['w1_median'] for r in usable
                                if r['w1_median'] is not None]),
            w1_range=describe([r['w1_range'] for r in usable
                               if r['w1_range'] is not None]),
            w1_relative_range=describe([r['w1_relative_range'] for r in usable
                                        if r['w1_relative_range'] is not None]),
            w1_range_zero=sum(1 for r in usable if r['w1_range'] == 0),
            support_exact_match_all=sum(1 for r in usable if r['support_exact_match_all']),
            support_exact_match_any=sum(1 for r in usable if r['support_exact_match_any']),
            fixed_seed_control=dict(
                comparable=len(ctrl),
                k_equal=sum(c['k_equal'] for c in ctrl),
                support_equal=sum(c['support_equal'] for c in ctrl),
                loss_equal=sum(c['loss_equal'] for c in ctrl),
                w1_equal=sum(c['w1_equal'] for c in ctrl)),
            writer_crashed_runs=sum(r['writer_crashed'] for r in rows),
            seconds=describe([r['seconds']['median'] for r in rows if r['seconds']]))
    out = dict(all=block(cases, 'all cases'))
    for family in sorted({r['family'] for r in cases}):
        out['family|' + family] = block([r for r in cases if r['family'] == family], family)
    for K in sorted({r['K_generating'] for r in cases}):
        out['K|%d' % K] = block([r for r in cases if r['K_generating'] == K], 'K=%d' % K)
    for family in sorted({r['family'] for r in cases}):
        for K in sorted({r['K_generating'] for r in cases if r['family'] == family}):
            sub = [r for r in cases if r['family'] == family and r['K_generating'] == K]
            out['cell|%s|K%d' % (family, K)] = block(sub, '%s K=%d' % (family, K))
    usable = [r for r in cases if r['seeds_complete'] >= 2]
    out['correlations'] = dict(
        k_range_vs_R=spearman([r['R_size'] for r in usable],
                              [r['k_range'] or 0 for r in usable]),
        support_distinct_vs_R=spearman([r['R_size'] for r in usable],
                                       [r['support_distinct'] for r in usable]),
        w1_relative_range_vs_R=spearman(
            [r['R_size'] for r in usable if r['w1_relative_range'] is not None],
            [r['w1_relative_range'] for r in usable if r['w1_relative_range'] is not None]))
    return out


def main():
    rows = load_runs()
    index = {e['name']: e for e in common.bank_index()}
    cases = per_case(rows, index)
    out = dict(runs=len(rows), cases=len(cases),
               run_status=dict(Counter(r.get('status') or r.get('row_status')
                                       for r in rows)),
               arms=dict(Counter(r.get('arm') for r in rows if r.get('kind') == 'dd_seed')),
               aggregate=aggregate(cases), per_case=cases)
    OUT.write_text(json.dumps(out, indent=1, default=str))
    print('written', OUT)
    print('runs', out['runs'], out['run_status'])
    a = out['aggregate']['all']
    for key in ('cases', 'cases_with_two_or_more_seeds', 'k_identical_across_seeds',
                'support_identical_across_seeds', 'k_distinct_distribution',
                'support_distinct_distribution', 'fixed_seed_control', 'w1_range_zero',
                'writer_crashed_runs'):
        print('  %-34s %s' % (key, a[key]))
    print('  w1_relative_range', a['w1_relative_range'])


if __name__ == '__main__':
    main()
