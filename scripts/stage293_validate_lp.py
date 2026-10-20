#!/usr/bin/env python3
"""Stage 293 validation: the exact simplex against an independent float LP.

Three checks, none of which trusts the other implementation's flags:
  1. random dense standard-form LPs: exact optimum vs scipy HiGHS optimum;
  2. exact optimal w satisfies M w = b and w >= 0 in Fractions;
  3. on real bank cases, the fibre's p/x bounds bracket the generating
     parameters, and the fibre's free_parameters agrees with the rank
     deficiency reported independently by stage289_support_feasibility.
"""
import os
for _k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_k] = '1'
import json
import random
import sys
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import numpy as np
from scipy.optimize import linprog
import stage288_common as common
import stage288_incompleteness as inc
from stage293_fibre import lp_max, LPError, fibre_report
from stage289_support_feasibility import feasible_support

OUT = ROOT / 'results/stage293/stage293_lp_validation.json'


def random_lp_checks(trials=300, seed=2930001):
    rng = random.Random(seed)
    agree, infeasible_both, mismatch = 0, 0, []
    for trial in range(trials):
        m, n = rng.randint(1, 4), rng.randint(2, 7)
        matrix = [[F(rng.randint(-4, 4)) for _ in range(n)] for _ in range(m)]
        w0 = [F(rng.randint(0, 5)) for _ in range(n)]          # a feasible point
        rhs = [sum((a * b for a, b in zip(row, w0)), F(0)) for row in matrix]
        cost = [F(rng.randint(-5, 5)) for _ in range(n)]
        try:
            value, w = lp_max(matrix, rhs, cost)
            exact_status = 'optimal'
        except LPError as exc:
            value, w, exact_status = None, None, str(exc)
        answer = linprog([-float(v) for v in cost],
                         A_eq=np.array([[float(v) for v in row] for row in matrix]),
                         b_eq=np.array([float(v) for v in rhs]),
                         bounds=(0, None), method='highs')
        if exact_status == 'optimal':
            residual = max(abs(sum((a * b for a, b in zip(row, w)), F(0)) - r)
                           for row, r in zip(matrix, rhs))
            assert residual == 0, 'exact solution violates the equalities'
            assert all(v >= 0 for v in w), 'exact solution is negative'
            if answer.status == 0 and abs(-answer.fun - float(value)) <= 1e-6 * max(1.0, abs(float(value))):
                agree += 1
            elif answer.status == 3:                 # float solver says unbounded
                mismatch.append(dict(trial=trial, exact=str(value), scipy='unbounded'))
            else:
                mismatch.append(dict(trial=trial, exact=str(value),
                                     scipy=(None if answer.fun is None else -answer.fun),
                                     scipy_status=int(answer.status)))
        elif exact_status == 'unbounded':
            if answer.status == 3:
                agree += 1
            else:
                mismatch.append(dict(trial=trial, exact='unbounded', scipy_status=int(answer.status)))
        else:
            infeasible_both += 1
    return dict(trials=trials, agreed=agree, exact_infeasible=infeasible_both,
                mismatches=mismatch[:20], n_mismatch=len(mismatch))


def bank_checks(limit=12):
    out = []
    for entry in inc.panel()[:limit]:
        if entry['reference_set_size'] > 200:
            continue
        case = common.load_case(entry)
        truth = case['generating_sources']
        support = sorted(s['edge'] for s in truth)
        report = fibre_report(case['edges'], case['observed'], support, truth=truth)
        verdict = feasible_support(case['edges'], case['observed'], support)
        row = dict(name=entry['name'], n=entry['reference_set_size'], k=len(truth),
                   status=report['status'], free_parameters=report.get('free_parameters'),
                   sf_identifiable=verdict.get('identifiable'),
                   sf_free_parameters=verdict.get('free_parameters'),
                   sf_verdict=verdict['verdict'],
                   truth_inside_bounds=report.get('truth_inside_bounds'),
                   p_sup_range=report.get('p_sup_range'),
                   x_sup_range=report.get('x_sup_range'))
        row['free_parameters_agree'] = (
            report.get('free_parameters') == (0 if verdict.get('identifiable')
                                              else verdict.get('free_parameters')))
        out.append(row)
    return out


if __name__ == '__main__':
    result = dict(random_lp=random_lp_checks(), bank=bank_checks())
    result['all_truth_inside'] = all(r['truth_inside_bounds'] is not False for r in result['bank'])
    result['all_free_parameters_agree'] = all(r['free_parameters_agree'] for r in result['bank'])
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(result, indent=1))
    print(json.dumps({k: v for k, v in result['random_lp'].items() if k != 'mismatches'}))
    print('truth inside bounds:', result['all_truth_inside'],
          '| free-parameter counts agree:', result['all_free_parameters_agree'])
    for row in result['bank']:
        print(row['name'][:42], 'k=%s' % row['k'], 'fp=%s' % row['free_parameters'],
              'sf_id=%s' % row['sf_identifiable'], 'p_range=%s' % row['p_sup_range'],
              'x_range=%s' % row['x_sup_range'])
