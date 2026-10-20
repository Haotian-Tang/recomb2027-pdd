#!/usr/bin/env python3
"""Stage 299 step 4: the exact-model comparator, from frozen stage-296 records.

Stage 296 already ran the SAME removal draws through the SAME instrument on 43
exactly explainable cases (exact_biotrees and exact_10k, K in {2, 3}): 1,075
draws plus 43 full panels, each with Theta_1 and its attaining witness.  This
script re-analyses those frozen records - no new runs - so the real-data arms
of stage 299 have an exact-model baseline measured on the identical statistic.

It is NOT the only comparator: the stage-299 exact_partner arm re-runs the
exact distances on the krepp trees and the krepp leaf sets, pairing exact and
estimated draw by draw.  Both are reported.
"""
import json
import sys
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import stage288_common as common               # noqa: E402
import stage299_common as s299                 # noqa: E402

OUT = s299.OUT
CROSS = ROOT / 'results' / 'stage296' / 'stage296_cross_runs.jsonl'


def med(values):
    s = sorted(values)
    return None if not s else (s[len(s) // 2] if len(s) % 2 else
                               (s[len(s) // 2 - 1] + s[len(s) // 2]) / 2)


def quart(values):
    s = sorted(values)
    if not s:
        return (None, None, None)
    return (s[len(s) // 4], med(s), s[(3 * len(s)) // 4])


def main():
    out_path = OUT / 'stage299_exact_baseline.json'
    assert not out_path.exists(), 'Never overwrite frozen evidence'
    rows = [json.loads(line) for line in CROSS.read_text().splitlines() if line.strip()]
    index = s299.bank_by_name()
    byname = defaultdict(dict)
    for r in rows:
        res = r['result']
        if res.get('status') != 'complete' or '1' not in res.get('thetas', {}):
            continue
        wit = res['witnesses']['1']
        key = (r['level'], r['regime'], r['replicate'])
        byname[r['name']][key] = dict(
            theta1=F(res['thetas']['1']), support=wit['support'][0],
            x=F(wit['sources'][0]['x']), n_retained=r['n_retained'],
            n_removed=r['n_removed'], K=r['K'], N=r['N'], family=r['family'],
            median_observed=F(r['median_observed']))
    cells = defaultdict(list)
    per_draw, rises = [], []
    trees = {}
    for name, draws in sorted(byname.items()):
        base = draws.get(('full', 'none', 0))
        if base is None:
            continue
        if name not in trees:
            entry = index[name]
            case = s299.prepare_case(s299.load_json(entry['path'], entry['sha256']), name)
            trees[name] = case['edges']
        edges = trees[name]
        deg = s299.degrees(edges)
        for (level, regime, rep), d in sorted(draws.items()):
            if level == 'full':
                continue
            fall = d['theta1'] / base['theta1'] if base['theta1'] else None
            if d['theta1'] > base['theta1']:
                rises.append(dict(name=name, level=level, regime=regime, replicate=rep,
                                  full=str(base['theta1']), draw=str(d['theta1'])))
            mig = s299.point_distance(edges, base['support'], base['x'], d['support'], d['x'])
            row = dict(name=name, family=d['family'], K=d['K'], N=d['N'],
                       level=level, regime=regime, replicate=rep,
                       n_removed=d['n_removed'], n_retained=d['n_retained'],
                       removed_fraction=d['n_removed'] / float(d['n_removed'] + d['n_retained']),
                       theta1=float(d['theta1']), theta1_full=float(base['theta1']),
                       fall_ratio=(float(fall) if fall is not None else None),
                       strict_fall=bool(d['theta1'] < base['theta1']),
                       migration=float(mig),
                       migration_rel=float(mig / base['median_observed']),
                       same_edge=bool(d['support'] == base['support']),
                       hops=s299.edge_hops(edges, base['support'], d['support']),
                       argmin_is_terminal=s299.is_terminal(deg, edges[d['support']]),
                       full_argmin_is_terminal=s299.is_terminal(deg, edges[base['support']]))
            per_draw.append(row)
            cells[(level, regime)].append(row)
    table = []
    for (level, regime), rs in sorted(cells.items()):
        fr = [r['fall_ratio'] for r in rs if r['fall_ratio'] is not None]
        q = quart(fr)
        table.append(dict(level=level, regime=regime, n=len(rs),
                          fall_ratio_q1=q[0], fall_ratio_median=q[1], fall_ratio_q3=q[2],
                          relative_drop_median=(None if q[1] is None else 1 - q[1]),
                          strict_fall_fraction=sum(r['strict_fall'] for r in rs) / float(len(rs)),
                          same_edge_fraction=sum(r['same_edge'] for r in rs) / float(len(rs)),
                          migration_rel_median=med([r['migration_rel'] for r in rs]),
                          hops_median=med([r['hops'] for r in rs]),
                          terminal_fraction=sum(r['argmin_is_terminal'] for r in rs) / float(len(rs))))
    pooled = {}
    for regime in sorted({r['regime'] for r in per_draw}):
        rs = [r for r in per_draw if r['regime'] == regime]
        q = quart([r['fall_ratio'] for r in rs])
        pooled[regime] = dict(n=len(rs), fall_ratio_median=q[1],
                              relative_drop_median=(None if q[1] is None else 1 - q[1]),
                              strict_fall_fraction=sum(r['strict_fall'] for r in rs) / float(len(rs)),
                              migration_rel_median=med([r['migration_rel'] for r in rs]),
                              same_edge_fraction=sum(r['same_edge'] for r in rs) / float(len(rs)))
    out = dict(stage=299, step='exact_baseline_from_stage296',
               source=str(CROSS.relative_to(ROOT)), source_sha256=common.sha(CROSS),
               n_records=len(rows), n_cases=len(byname), n_draws=len(per_draw),
               monotonicity=dict(rises=len(rises), examples=rises[:5]),
               by_level_regime=table, pooled_by_regime=pooled, per_draw=per_draw)
    out_path.write_text(json.dumps(common.serial(out), indent=1))
    print('cases %d, draws %d, rises %d' % (len(byname), len(per_draw), len(rises)))
    for regime, v in sorted(pooled.items()):
        print('  %-8s n=%4d median fall ratio %.4f  strict-fall %.3f  migration/med %.4f'
              % (regime, v['n'], v['fall_ratio_median'], v['strict_fall_fraction'],
                 v['migration_rel_median']))


if __name__ == '__main__':
    main()

