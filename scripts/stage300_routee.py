#!/usr/bin/env python3
"""Stage 300 verification route E: Lambda_1 against DecoDiPhy's own k = 1 answer.

At k = 1 DecoDiPhy hill-climbs over the same node set that `best_single_source_l2`
enumerates exhaustively, and optimises the same objective, so the loss in the
first rung of its ladder must be >= Lambda_1 up to floating point.  Equality on
most cases also checks the position/pendant conventions end to end: a sign or
orientation error would show up as a systematic gap.
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import stage300_common as s300                   # noqa: E402
import stage300_evaluate as s300eval             # noqa: E402


def main():
    report = dict(stage=300, step='route_e', by_collection={})
    for collection in s300.COLLECTIONS:
        base_path = s300.OUT / ('stage300_baselines_%s.jsonl' % collection)
        if not base_path.exists():
            continue
        lam = {}
        for line in open(base_path):
            row = json.loads(line)
            if row.get('lambda1'):
                lam[row['name']] = row['lambda1']['ssr']
        pairs = []
        for run in s300eval.runs_for(collection):
            if run.get('status') != 'complete' or not run.get('ladder'):
                continue
            first = [r for r in run['ladder'] if r.get('k') == 1]
            if not first or run['name'] not in lam:
                continue
            dd1 = float(first[0]['loss'])
            mine = lam[run['name']]
            pairs.append(dict(name=run['name'], arm=run['arm'],
                              source_stage=run['source_stage'],
                              dd_k1_loss=dd1, lambda1=mine,
                              gap=dd1 - mine,
                              rel=(dd1 - mine) / max(abs(mine), 1e-300)))
        if not pairs:
            continue
        rel = sorted(p['rel'] for p in pairs)
        report['by_collection'][collection] = dict(
            n=len(pairs),
            n_dd_below_lambda1=sum(1 for p in pairs if p['rel'] < -1e-9),
            n_equal_within_1e_9=sum(1 for p in pairs if abs(p['rel']) <= 1e-9),
            rel_min=rel[0], rel_median=rel[len(rel) // 2], rel_max=rel[-1],
            worst_below=sorted([p for p in pairs if p['rel'] < -1e-9],
                               key=lambda p: p['rel'])[:5],
            largest_above=sorted(pairs, key=lambda p: -p['rel'])[:5])
        cell = report['by_collection'][collection]
        print('%-18s n=%4d  dd<Lambda1: %d  equal: %d  rel median %.2e  max %.2e'
              % (collection, cell['n'], cell['n_dd_below_lambda1'],
                 cell['n_equal_within_1e_9'], cell['rel_median'], cell['rel_max']),
              flush=True)
    (s300.OUT / 'stage300_route_e.json').write_text(json.dumps(s300.serial(report), indent=1))
    print('wrote', s300.OUT / 'stage300_route_e.json')


if __name__ == '__main__':
    main()
