#!/usr/bin/env python3
"""Reproduce the float early-return in stage301_common.choose_constraints.

The function scores constraints in float and keeps the tightest, recomputing
those EXACTLY.  That design is sound, but the float pass carries its own
rejection:

    con = envelope_constraint(view, s, eps, exact=False)
    if con['rho'] < 0:
        return None, dict(infeasible_at=s)

so a rho that is positive in exact arithmetic but negative after rounding is
rejected before the exact recomputation can run.  This script builds a strictly
feasible single-source instance on which that happens and checks the exact
recheck repair. It exits with an error if the witness fails or either
implementation wrongly rejects either tolerance. Requires numpy, through
stage301_common.constraint_mask.

Three-leaf star, centre 0, edge lengths 3/10, 3/5, 7/10.  One source on the
first edge at x = 1 - 1e-30 from the centre, p = 1, y = 1.  The source is
strictly interior, the weight and the pendant are strictly positive, so the
instance is inside the model.
"""
import os
import sys
from fractions import Fraction as F

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
import stage301_common as s301                          # noqa: E402

DELTA = F(3, 10 ** 31)


def instance():
    edges = [(0, 1, F(3, 10)), (0, 2, F(3, 5)), (0, 3, F(7, 10))]
    observed = {1: F(1) + DELTA, 2: F(19, 10) - DELTA, 3: F(2) - DELTA}
    return edges, observed


def witness_check(edges, observed):
    """Confirm the observations really are the model's values, exactly."""
    x = 1 - F(1, 10 ** 30)            # fraction of edge (0,1) from the centre
    lengths = {(u, v): length for u, v, length in edges}
    L1 = lengths[(0, 1)]
    y = F(1)
    if not (0 < x < 1 and y > 0 and all(length > 0 for length in lengths.values())):
        raise AssertionError('the witness must be strictly inside the model')
    # distance from the source to each leaf, by the original-edge formula
    d = {1: y + (1 - x) * L1,
         2: y + x * L1 + lengths[(0, 2)],
         3: y + x * L1 + lengths[(0, 3)]}
    return d == observed, d


def patched_choose(view, eps, n_edges, sample, n_keep, rng):
    """choose_constraints with the exact recheck before rejecting."""
    leaves = view['leaves']
    pick = list(leaves) if len(leaves) <= sample else rng.sample(list(leaves), sample)
    scored = []
    for s in pick:
        con = s301.envelope_constraint(view, s, eps, exact=False)
        if con['rho'] < 0:
            con = s301.envelope_constraint(view, s, eps, exact=True)
            if con['rho'] < 0:
                return None, dict(infeasible_at=s)
        size = int(s301.constraint_mask(view, con, n_edges).sum())
        scored.append((size, s))
    scored.sort()
    cons = []
    for s in [t for _sz, t in scored[:n_keep]]:
        con = s301.envelope_constraint(view, s, eps, exact=True)
        if con['rho'] < 0:
            return None, dict(infeasible_at=s)
        cons.append(con)
    return cons, dict(sampled=len(pick), kept=len(cons))


def main():
    import random
    edges, observed = instance()
    ok, d = witness_check(edges, observed)
    print('witness reproduces the observations exactly:', ok)
    if not ok:
        raise AssertionError('the witness does not reproduce the observations')
    view = s301.backbone_view(edges, observed)
    n_edges = len(edges)

    bad = good = 0
    for eps in (F(0), F(1, 10 ** 32)):
        ex = [s301.envelope_constraint(view, s, eps, exact=True)['rho']
              for s in view['leaves']]
        fl = [s301.envelope_constraint(view, s, eps, exact=False)['rho']
              for s in view['leaves']]
        if min(ex) != DELTA + eps:
            raise AssertionError('unexpected exact envelope margin')
        cur, _ = s301.choose_constraints(view, eps, n_edges, 10, 3,
                                         random.Random(0))
        fix, _ = patched_choose(view, eps, n_edges, 10, 3, random.Random(0))
        bad += int(cur is None)
        good += int(fix is None)
        print('eps=%-12s min exact rho=%s  min float rho=%.3g  current=%s  patched=%s'
              % (eps, min(ex), float(min(fl)),
                 'infeasible' if cur is None else 'ok',
                 'infeasible' if fix is None else 'ok'))
    print('wrong rejections: current %d/2, patched %d/2' % (bad, good))
    if bad or good:
        raise AssertionError('a verified strictly feasible instance was rejected')


if __name__ == '__main__':
    main()
