#!/usr/bin/env python3
"""Stage 310: one case where support uniqueness fails before the count does.

Across the adversarial removal measurements the first loss of uniqueness is far
more often a competing support of the SAME size than a drop in the minimum
count.  This script exhibits a single case in full so the mechanism can be read
off rather than inferred: it recomputes the complete set of minimum-size strict
supports on the full panel and after each deletion in the case's own
nearest-to-source order, using the frozen all-minimisers routine.

Writes results/stage310/stage310_uniqueness_before_count.json.
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import stage300_common as s300                          # noqa: E402
from stage292_all_minimisers import all_minimum_supports  # noqa: E402

CASE = 'bio_bees-bosseret__bees-bosseret.tre_k2_s3_noise0_exact'
ORDER = [20, 13, 17]


def minsupports(edges, leaves, dvals, cap):
    observed = {r: d for r, d in zip(leaves, dvals)}
    res = all_minimum_supports(edges, observed, cap=cap)
    if not isinstance(res, dict) or res.get('status') != 'complete':
        raise RuntimeError('all-minimisers enumeration did not complete: %s' % res)
    sups = sorted(tuple(sorted(s)) for s in res['supports'])
    if (not sups or res.get('n_supports') != len(sups)
            or len(set(sups)) != len(sups)
            or any(len(s) != res.get('minimum') for s in sups)):
        raise RuntimeError('inconsistent completed minimum-support result: %s' % res)
    return res, sups


def main():
    import stage288_common as common
    entry = next((dict(e) for e in common.bank_index() if e['name'] == CASE), None)
    if entry is None:
        raise SystemExit('case not found: %s' % CASE)

    case = s300.load_case(entry)
    edges, leaves, dvals = s300.case_arrays(case)
    truth = s300.truth_configuration(case)
    k_true = len(truth[0])
    weights = [str(s['p']) for s in truth[0]]
    print('%s\n  %d observed references, K=%d, weights %s'
          % (CASE, len(leaves), k_true, ', '.join(weights)))

    steps, keep = [], list(leaves)
    for i in range(len(ORDER) + 1):
        if i:
            drop = ORDER[i - 1]
            if drop not in keep:
                raise SystemExit('reference %d is not in the panel' % drop)
            keep.remove(drop)
        idx = [leaves.index(r) for r in keep]
        res, sups = minsupports(edges, keep, [dvals[j] for j in idx], k_true)
        steps.append(dict(removed=ORDER[:i], n_references=len(keep),
                          minimum=res['minimum'],
                          n_supports=len(sups),
                          supports=[list(s) for s in sups],
                          status=res['status']))
        print('  removed %-12s refs %2d  minimum %s  supports %s'
              % (ORDER[:i] or 'none', len(keep),
                 steps[-1]['minimum'], steps[-1]['supports']))

    out_dir = os.path.join(HERE, os.pardir, 'results', 'stage310')
    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir,
                           'stage310_uniqueness_before_count.json'), 'w') as fh:
        json.dump(dict(stage=310, kind='uniqueness_before_count', case=CASE,
                       k_true=k_true, weights=weights,
                       deletion_order=ORDER, steps=steps), fh,
                  indent=1, sort_keys=True)


if __name__ == '__main__':
    main()
