#!/usr/bin/env python3
"""Stage 310: one case where support uniqueness fails before the count does.

Across the adversarial removal measurements the first loss of uniqueness is far
more often a competing support of the SAME size than a drop in the minimum
count.  This script exhibits a single case in full so the mechanism can be read
off rather than inferred: it recomputes the complete set of minimum-size strict
supports on the full panel and after each deletion in the case's own
nearest-to-source order, using the frozen all-minimisers routine.

It also distinguishes that order's safe prefix from the worst-case safe radius
by checking every single-reference deletion and one destroying pair.

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


def safe_prefix_and_radius(steps, singles, pair, true_support):
    """Certify preservation of the original minimum support, not just a count."""
    def preserves(row):
        return (row['status'] == 'complete'
                and row['minimum'] == len(true_support)
                and row['supports'] == [true_support])

    if not steps or not preserves(steps[0]):
        raise RuntimeError('full panel does not have the unique true support')
    failures = [i for i, row in enumerate(steps) if not preserves(row)]
    if not failures:
        raise RuntimeError('deletion order did not reach a loss of uniqueness')
    prefix = failures[0] - 1
    all_single_unique = bool(singles) and all(preserves(row) for row in singles)
    radius = (1 if all_single_unique and pair['status'] == 'complete'
              and not preserves(pair) else None)
    return prefix, all_single_unique, radius


def main():
    import stage288_common as common
    entry = next((dict(e) for e in common.bank_index() if e['name'] == CASE), None)
    if entry is None:
        raise SystemExit('case not found: %s' % CASE)

    case = s300.load_case(entry)
    edges, leaves, dvals = s300.case_arrays(case)
    truth = s300.truth_configuration(case)
    k_true = len(truth[0])
    true_support = sorted(s['edge'] for s in truth[0])
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

    # One order supplies an upper bound on the arbitrary-deletion guarantee.
    # Every singleton plus one destroying pair certifies radius one; checking
    # all pairs is unnecessary.  All calls must complete, via minsupports.
    singles = []
    for r in leaves:
        keep = [x for x in leaves if x != r]
        idx = [leaves.index(x) for x in keep]
        res, sups = minsupports(edges, keep, [dvals[j] for j in idx], k_true)
        singles.append(dict(removed=[r], n_references=len(keep),
                            n_supports=len(sups), minimum=res['minimum'],
                            supports=[list(s) for s in sups],
                            status=res['status']))
    pair = [ORDER[0], ORDER[2]]
    keep = [x for x in leaves if x not in pair]
    idx = [leaves.index(x) for x in keep]
    res, pair_sups = minsupports(edges, keep, [dvals[j] for j in idx], k_true)
    pair_record = dict(removed=pair, n_references=len(keep),
                       n_supports=len(pair_sups), minimum=res['minimum'],
                       supports=[list(s) for s in pair_sups],
                       status=res['status'])
    prefix, all_single_unique, radius = safe_prefix_and_radius(
        steps, singles, pair_record, true_support)
    print('  all %d single-reference deletions keep a unique support: %s'
          % (len(singles), all_single_unique))
    print('  deleting %s directly leaves %d minimum supports %s'
          % (pair, len(pair_sups), pair_record['supports']))
    print('  nearest-order safe prefix %d, worst-case safe erasure radius %s'
          % (prefix, radius))

    out_dir = os.path.join(HERE, os.pardir, 'results', 'stage310')
    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir,
                           'stage310_uniqueness_before_count.json'), 'w') as fh:
        json.dump(dict(stage=310, kind='uniqueness_before_count', case=CASE,
                       k_true=k_true, weights=weights,
                       deletion_order=ORDER, steps=steps,
                       single_deletions=singles, pair_deletion=pair_record,
                       all_single_deletions_unique=all_single_unique,
                       nearest_order_safe_prefix=prefix,
                       worst_case_safe_radius=radius), fh,
                  indent=1, sort_keys=True)


if __name__ == '__main__':
    main()
