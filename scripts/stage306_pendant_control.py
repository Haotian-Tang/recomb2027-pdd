#!/usr/bin/env python3
"""Stage 306: pendant-only control for the generating-configuration baseline.

Stage 300 scores the generating configuration with its true edges, positions,
weights AND its true common pendant y, while Theta_1 and Lambda_1 re-optimise
both the position and the pendant of their single source.  The comparison is
therefore not like for like.  This script grants the generating configuration
the one free parameter the baselines already have -- the common pendant -- and
nothing else, then repeats both comparisons.

With g_r = sum_j p_j dist_T(q_j, r) the geometric part of the generating
configuration and a_r = dhat_r - g_r, the pendant that minimises the worst-leaf
residual is the clipped midrange and the one that minimises the sum of squares
is the clipped mean:

    y*_inf = max(0, (min_r a_r + max_r a_r) / 2)
    y*_2   = max(0, mean_r a_r)

Both are closed forms; no search and no solver is involved.  y* = 0 is the
infimum of the strict y > 0 model rather than an attained minimum, so the
script reports how many optima sit on that boundary.

Residuals are exact rationals, as in stage 300; the sum of squares is float64,
as Lambda_1 is.  Writes results/stage306/stage306_pendant_control.json.
"""
import json
import os
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import stage300_common as s300                      # noqa: E402

COLLECTION = 'krepp_ani'


def median(xs):
    s = sorted(xs)
    n = len(s)
    if not n:
        return None
    return s[n // 2] if n % 2 else (s[n // 2 - 1] + s[n // 2]) / 2


def main():
    theta = s300.theta1_table()
    entries = s300.case_entries(COLLECTION)
    base = {}
    path = s300.OUT / ('stage300_baselines_%s.jsonl' % COLLECTION)
    with open(path) as fh:
        for line in fh:
            row = json.loads(line)
            base[row['name']] = row

    rows = []
    for entry in entries:
        name = entry['name']
        case = s300.load_case(entry)
        edges, leaves, dvals = s300.case_arrays(case)
        truth = s300.truth_configuration(case)
        if truth is None:
            continue
        sources, y0 = truth
        res = s300.residuals_exact(edges, leaves, dvals, sources, y0)
        linf_fixed = max(abs(z) for z in res)
        b = base.get(name, {})
        stored = b.get('truth', {}).get('linf')
        if stored is not None and F(stored) != linf_fixed:
            raise SystemExit('stage300 mismatch on %s' % name)

        # a_r = dhat_r - g_r = y0 - res_r
        a = [y0 - z for z in res]
        y_inf = max(F(0), (min(a) + max(a)) / 2)
        y_two = max(F(0), sum(a) / len(a))
        linf_free = max(abs(y_inf - z) for z in a)
        ssr_fixed = float(sum(float(z) ** 2 for z in res))
        ssr_free = float(sum((float(y_two) - float(z)) ** 2 for z in a))

        t1 = theta.get(name)
        theta1 = t1['theta1'] if t1 else None
        lam = b.get('lambda1') or {}
        lam_ssr = lam.get('ssr')

        rows.append(dict(
            name=name, k_true=len(sources), n_leaves=len(leaves),
            y_true=str(y0), y_inf=str(y_inf), y_two=str(y_two),
            y_inf_is_zero=(y_inf == 0), y_two_is_zero=(y_two == 0),
            theta1=(str(theta1) if theta1 is not None else None),
            linf_fixed=str(linf_fixed), linf_free=str(linf_free),
            ratio_fixed=(float(linf_fixed / theta1) if theta1 else None),
            ratio_free=(float(linf_free / theta1) if theta1 else None),
            above_theta1_fixed=(bool(linf_fixed > theta1) if theta1 else None),
            above_theta1_free=(bool(linf_free > theta1) if theta1 else None),
            ties_theta1_free=(bool(linf_free == theta1) if theta1 else None),
            lambda1_ssr=lam_ssr, ssr_fixed=ssr_fixed, ssr_free=ssr_free,
            above_lambda1_fixed=(bool(ssr_fixed > lam_ssr)
                                 if lam_ssr is not None else None),
            above_lambda1_free=(bool(ssr_free > lam_ssr)
                                if lam_ssr is not None else None),
        ))

    have_t = [r for r in rows if r['theta1'] is not None]
    have_l = [r for r in rows if r['lambda1_ssr'] is not None]
    beat_t = [r for r in have_t if not r['above_theta1_free']]
    beat_l = [r for r in have_l if not r['above_lambda1_free']]
    summary = dict(
        collection=COLLECTION, n_cases=len(rows),
        linf=dict(
            n=len(have_t),
            above_theta1_fixed=sum(1 for r in have_t if r['above_theta1_fixed']),
            above_theta1_free=sum(1 for r in have_t if r['above_theta1_free']),
            ties_free=sum(1 for r in have_t if r['ties_theta1_free']),
            median_ratio_fixed=median([r['ratio_fixed'] for r in have_t]),
            median_ratio_free=median([r['ratio_free'] for r in have_t]),
            n_y_zero=sum(1 for r in have_t if r['y_inf_is_zero']),
            n_beat=len(beat_t),
            n_beat_with_y_zero=sum(1 for r in beat_t if r['y_inf_is_zero']),
            beat_min_gap=min([r['ratio_free'] for r in beat_t], default=None),
        ),
        ssr=dict(
            n=len(have_l),
            above_lambda1_fixed=sum(1 for r in have_l if r['above_lambda1_fixed']),
            above_lambda1_free=sum(1 for r in have_l if r['above_lambda1_free']),
            n_y_zero=sum(1 for r in have_l if r['y_two_is_zero']),
            n_beat=len(beat_l),
            n_beat_with_y_zero=sum(1 for r in beat_l if r['y_two_is_zero']),
        ),
    )
    out_dir = s300.OUT.parent / 'stage306'
    out_dir.mkdir(parents=True, exist_ok=True)
    payload = dict(stage=306, kind='pendant_only_control',
                   summary=summary, cases=rows)
    with open(out_dir / 'stage306_pendant_control.json', 'w') as fh:
        json.dump(payload, fh, indent=1, sort_keys=True)

    s = summary['linf']
    print('worst leaf : above Theta_1  %d/%d fixed -> %d/%d refitting y only'
          % (s['above_theta1_fixed'], s['n'], s['above_theta1_free'], s['n']))
    print('             median ratio  %.4f -> %.4f   ties %d'
          % (s['median_ratio_fixed'], s['median_ratio_free'], s['ties_free']))
    print('             y* = 0 in %d of %d; %d beat Theta_1, %d of them at y* = 0'
          % (s['n_y_zero'], s['n'], s['n_beat'], s['n_beat_with_y_zero']))
    q = summary['ssr']
    print('sum squares: above Lambda_1 %d/%d fixed -> %d/%d refitting y only'
          % (q['above_lambda1_fixed'], q['n'], q['above_lambda1_free'], q['n']))
    print('             y* = 0 in %d of %d; %d beat Lambda_1, %d of them at y* = 0'
          % (q['n_y_zero'], q['n'], q['n_beat'], q['n_beat_with_y_zero']))


if __name__ == '__main__':
    main()
