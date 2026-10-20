#!/usr/bin/env python3
"""Stage297 add-on: does the stage 297 harness reproduce stage 289's DecoDiPhy
answers at the documented default seed?

Stage294 reported that the support DecoDiPhy returns did not agree between its
own harness and stage 289's (575 compared, 567 agree on k, 471 on support) at
the same default seed.  The stage 297 sweep contains a seed-1142 arm run through
the stage 293 worker unmodified, so the same comparison can be made again for
the cases the two sets share.  Nothing is re-fitted: both sides are read off
their recorded runs.
"""
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
sys.dont_write_bytecode = True

OUT = ROOT / 'results/stage297/stage297_dd_vs_stage289.json'


def main():
    mine = {}
    for p in sorted((ROOT / 'results/stage297').glob('stage297_dd_seeds_*.jsonl')):
        for line in p.open():
            r = json.loads(line)
            if r.get('kind') == 'dd_seed' and r['arm'] in ('default_1142', 'seeded_1142'):
                mine.setdefault(r['name'], {})[r['arm']] = r
    theirs = {}
    for line in (ROOT / 'results/stage289/stage289_headtohead_runs.jsonl').open():
        r = json.loads(line)
        if r['method'] != 'dd_default':
            continue
        theirs.setdefault(r['name'], {})[r['repetition']] = r['result']
    rows = []
    for name in sorted(set(mine) & set(theirs)):
        a = mine[name].get('default_1142')
        reps = theirs[name]
        for rep, b in sorted(reps.items()):
            if a.get('status') != 'complete' or b.get('status') != 'complete':
                rows.append(dict(name=name, repetition=rep, comparable=False,
                                 mine_status=a.get('status'), theirs_status=b.get('status')))
                continue
            sup_a = sorted(int(s['edge']) for s in a['sources'])
            sup_b = sorted(int(s['edge']) for s in b['sources'])
            rows.append(dict(
                name=name, family=a['family'], K_generating=a['K_generating'],
                R_size=a['R_size'], repetition=rep, comparable=True,
                k_mine=a['k'], k_theirs=b['k'], k_equal=bool(a['k'] == b['k']),
                support_equal=bool(sup_a == sup_b),
                jaccard=(len(set(sup_a) & set(sup_b)) / len(set(sup_a) | set(sup_b))
                         if (sup_a or sup_b) else 1.0),
                loss_mine=a['loss'], loss_theirs=b['loss'],
                loss_equal=bool(a['loss'] == b['loss']),
                support_mine=sup_a, support_theirs=sup_b))
    good = [r for r in rows if r['comparable']]
    stage289_internal = []
    for name, reps in sorted(theirs.items()):
        if len(reps) < 2:
            continue
        b1, b2 = reps[1], reps[2]
        if b1.get('status') != 'complete' or b2.get('status') != 'complete':
            continue
        s1 = sorted(int(s['edge']) for s in b1['sources'])
        s2 = sorted(int(s['edge']) for s in b2['sources'])
        stage289_internal.append(dict(name=name, k_equal=bool(b1['k'] == b2['k']),
                                      support_equal=bool(s1 == s2)))
    out = dict(
        cases_in_both=len(set(mine) & set(theirs)),
        comparisons=len(rows), comparable=len(good),
        k_equal=sum(r['k_equal'] for r in good),
        support_equal=sum(r['support_equal'] for r in good),
        loss_equal=sum(r['loss_equal'] for r in good),
        jaccard_median=(sorted(r['jaccard'] for r in good)[len(good) // 2] if good else None),
        by_family={f: dict(
            n=sum(1 for r in good if r['family'] == f),
            k_equal=sum(r['k_equal'] for r in good if r['family'] == f),
            support_equal=sum(r['support_equal'] for r in good if r['family'] == f))
            for f in sorted({r['family'] for r in good})},
        disagreements=[{k: v for k, v in r.items() if k != 'jaccard'}
                       for r in good if not (r['k_equal'] and r['support_equal'])][:50],
        n_disagreements=sum(1 for r in good if not (r['k_equal'] and r['support_equal'])),
        stage289_internal_repetitions=dict(
            n=len(stage289_internal),
            k_equal=sum(r['k_equal'] for r in stage289_internal),
            support_equal=sum(r['support_equal'] for r in stage289_internal)),
        incomparable={'%s|%s' % k: v for k, v in
                      Counter((r.get('mine_status'), r.get('theirs_status'))
                              for r in rows if not r['comparable']).items()},
        note=('both sides are DecoDiPhy 1.5.0 at the documented default seed 1142, '
              'method dd_default, full panel; stage 297 ran it through '
              'stage293_abundance.dd_worker unmodified'))
    OUT.write_text(json.dumps(out, indent=1, default=str))
    print('written', OUT)
    for key in ('cases_in_both', 'comparable', 'k_equal', 'support_equal', 'loss_equal',
                'n_disagreements', 'by_family', 'stage289_internal_repetitions'):
        print('  %-32s %s' % (key, out[key]))


if __name__ == '__main__':
    main()
