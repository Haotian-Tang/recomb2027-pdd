#!/usr/bin/env python3
"""Reproduction B - a measurement printed in the paper, on real estimated distances.

PAPER CLAIM REPRODUCED
  sections/experiments.tex, paragraph "Real distances":
  "The worst-leaf residual of the true configuration exceeds Theta_1 in 296
   of 296 ANI cases, median ratio 2.9".
  appendices/experiments_details.tex, "Fit comparison on real distances":
  "On the 296 ANI cases the true configuration's worst-leaf residual is at or
   above Theta_1 in 296 of 296".

WHAT IT DOES
  For every one of the 296 krepp ANI bank cases it recomputes, from the case
  files alone,

    Theta_1    the exact smallest tolerance at which ONE source explains the
               estimated observations - stage 298's exact O(|V|) closed form
               (stage298_common.theta1_all_edges), validated in stage 298
               against the frozen per-support LP on 653,807 edges of these
               trees with 0 disagreements;
    residual   max_r |d_estimated(r) - d_exact(r)| over the shared leaf set.
               The generating configuration reproduces d_exact exactly, so this
               is precisely the worst-leaf residual of the true configuration
               against the estimated observations - stage 298's frozen formula
               (stage298_residuals.residual_row).

  It reports how many of the 296 have residual >= Theta_1 and the quantiles of
  the ratio, and checks both against the registered record
  results/stage298/stage298_analysis.json, key
  by_family/krepp_estimated_10k/residual_max_over_theta1.  Quantiles follow
  that record's stated convention (linear interpolation, numpy.percentile
  "linear" method), evaluated here in exact rational arithmetic.
  Every intermediate value is an exact rational; only the printed summary and
  the comparison tolerance are floating point.
  Dependencies: Python >= 3.12 standard library only.

USAGE
  python3 repro/repro_b_theta1_vs_residual.py          # all 296, about 6 s
  python3 repro/repro_b_theta1_vs_residual.py --n 40
"""
import argparse
import json
import sys
import time
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))

import stage288_common as common            # noqa: E402
import stage298_common as s298              # noqa: E402
import stage298_residuals as s298r          # noqa: E402

QUANTS = ('0.0', '0.25', '0.5', '0.75', '1.0')


def quantile_linear(sorted_vals, q):
    """numpy.percentile 'linear' method, in exact rational arithmetic."""
    n = len(sorted_vals)
    pos = F(q) * (n - 1)
    lo = int(pos)
    hi = min(lo + 1, n - 1)
    frac = pos - lo
    return sorted_vals[lo] + frac * (sorted_vals[hi] - sorted_vals[lo])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, default=0, help='limit to the first n cases (0 = all 296)')
    args = ap.parse_args()

    index = {e['name']: e for e in common.bank_index()}
    names = sorted(n for n in index if index[n]['family'] == 'krepp_estimated_10k')
    if args.n:
        names = names[:args.n]
    print('krepp ANI cases: %d' % len(names))

    ratios, n_ge, rows = [], 0, []
    t0 = time.perf_counter()
    for i, name in enumerate(names, 1):
        noisy = s298.read_case(index[name])
        exact = s298.read_case(index[name.replace('_krepp_estimated', '_exact')])
        resid = s298r.residual_row(exact, noisy)
        edges, leaves, dvals = s298.case_arrays(noisy)
        cf = s298.theta1_all_edges(edges, leaves, dvals)
        theta1 = cf['theta1']
        ratio = F(resid['residual_max']) / theta1
        ratios.append(ratio)
        n_ge += (ratio >= 1)
        rows.append(dict(name=name, n_leaves=len(leaves), n_edges=len(edges),
                         theta1=str(theta1), residual_max=str(resid['residual_max']),
                         ratio=float(ratio)))
        if i % 100 == 0 or i == len(names):
            print('  %3d/%d  %.1fs' % (i, len(names), time.perf_counter() - t0), flush=True)

    ratios.sort()
    qs = {q: quantile_linear(ratios, q) for q in QUANTS}
    rec = json.loads((ROOT / 'results' / 'stage298' / 'stage298_analysis.json').read_text())
    fam = rec['by_family']['krepp_estimated_10k']
    rec_q, rec_n = fam['residual_max_over_theta1'], fam['n_with_residual']

    print('\n' + '=' * 74)
    print('cases with worst-leaf residual >= Theta_1 : %d of %d   (paper: 296 of 296)'
          % (n_ge, len(names)))
    print('median ratio residual / Theta_1          : %.4f           (paper: 2.9)'
          % float(qs['0.5']))
    print('\nquantiles of residual / Theta_1, recomputed vs registered record:')
    print('  %-6s %18s %18s %s' % ('q', 'recomputed', 'recorded', 'match'))
    all_match = True
    for q in QUANTS:
        a, b = float(qs[q]), rec_q[q]
        m = abs(a - b) <= 1e-9 * max(1.0, abs(b))
        all_match &= m
        print('  %-6s %18.10f %18.10f %s' % (q, a, b, 'yes' if m else 'NO'))
    full = len(names) == rec_n
    print('\ncase count matches record (%d == %d): %s' % (len(names), rec_n, full))
    (ROOT / 'repro_b_output.json').write_text(json.dumps(
        dict(n_cases=len(names), n_residual_ge_theta1=n_ge,
             quantiles={q: float(qs[q]) for q in QUANTS},
             recorded_quantiles=rec_q, rows=rows), indent=1))
    print('per-case table written to repro_b_output.json')
    ok = all_match and (n_ge == len(names)) and full
    print('RESULT: ' + ('PASS - reproduces the paper claim' if ok else 'FAIL'))
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
