#!/usr/bin/env python3
"""Stage 296 analysis: eps-threshold distributions, the within-eps
incompleteness table, and the position of real estimation error on that axis.

Everything is derived from the exact rational vectors (Theta_1, Theta_2, ...)
retained by stage296_tolerance.py.  No eps ladder was run: for eps > 0 the
eps-tolerant minimum is m(eps) = min{k : Theta_k < eps}, so a ladder is a
reporting choice applied here, and the exact fall region of a draw against its
own full-panel baseline is the union over k of the half-open intervals
(Theta_k^draw, Theta_k^full].
"""
import json
import sys
from collections import Counter, defaultdict
from decimal import Decimal
from fractions import Fraction as F
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from stage288_common import ROOT, serial
from scipy.stats import spearmanr

OUT = ROOT / 'results/stage296'
LADDER = ('0', '1e-4', '3.16e-4', '1e-3', '3.16e-3', '1e-2', '3.16e-2',
          '1e-1', '3.16e-1', '1')


def read(path):
    return [json.loads(line) for line in Path(path).read_text().splitlines() if line.strip()]


def thetas_of(record):
    """Exact Theta_k for every k, using Theta_k = 0 for k at or beyond the first zero."""
    res = record['result']
    out = {int(k): F(v) for k, v in res['thetas'].items()}
    return out


def m_of_eps(th, K, eps, m0):
    if eps == 0:
        return m0
    for k in range(1, K + 1):
        if th.get(k, F(0)) < eps:
            return k
    return K


def quantiles(values, qs=(0.0, 0.25, 0.5, 0.75, 1.0)):
    if not values:
        return {}
    s = sorted(values)
    return {str(q): float(s[min(len(s) - 1, int(q * (len(s) - 1) + 0.5))]) for q in qs}


def main():
    protocol = json.loads((OUT / 'stage296_protocol.json').read_text())
    cross = read(OUT / 'stage296_cross_runs.jsonl')
    probe = read(OUT / 'stage296_probe_runs.jsonl')
    validate = read(OUT / 'stage296_validate_runs.jsonl')
    calib = json.loads((OUT / 'stage296_calibration.json').read_text())
    report = dict(counts=dict(cross=len(cross), probe=len(probe), validate=len(validate)))

    ok = [r for r in cross if r['result'].get('status') == 'complete']
    report['counts']['cross_complete'] = len(ok)
    report['counts']['cross_failed'] = len(cross) - len(ok)
    report['failures'] = [dict(name=r['name'], level=r['level'], regime=r['regime'],
                               status=r['result'].get('status'),
                               reason=r['result'].get('reason'))
                          for r in cross if r['result'].get('status') != 'complete'][:40]
    flags = [dict(name=r['name'], level=r['level'], regime=r['regime'],
                  flags=r['result'].get('parent_flags'))
             for r in cross + probe if r['result'].get('parent_flags')]
    report['witness_flags'] = flags[:40]
    report['witness_flag_count'] = len(flags)

    # ------------------------------------------------- 1. eps = 0 validation
    key = lambda r: (r['name'], r['level'], r['regime'], r['replicate'])
    frozen = {key(r): r['result'] for r in validate}
    agree = disagree = missing = censored = 0
    disagreements = []
    for r in ok:
        f = frozen.get(key(r))
        if f is None:
            missing += 1
            continue
        if f.get('status') != 'complete' or f.get('minimum') is None:
            censored += 1
            continue
        mine = r['result'].get('m_at_zero')
        if mine is None:
            censored += 1
            continue
        if mine == f['minimum']:
            agree += 1
        else:
            disagree += 1
            disagreements.append(dict(name=r['name'], level=r['level'], regime=r['regime'],
                                      replicate=r['replicate'], mine=mine, frozen=f['minimum'],
                                      n_retained=r['n_retained']))
    report['eps0_validation'] = dict(agree=agree, disagree=disagree, missing=missing,
                                     censored=censored, disagreements=disagreements[:40])

    # --------------------------------------- 2. thresholds at the full panel
    full = {r['name']: r for r in ok if r['level'] == 'full'}
    rows = []
    for name, r in sorted(full.items()):
        th = thetas_of(r)
        K, m0 = r['K'], r['result']['m_at_zero']
        med = F(r['median_observed'])
        row = dict(name=name, family=r['family'], cell=r['cell'], K=K, N=r['N'],
                   n_edges=r['n_edges'], m_at_zero=m0,
                   min_generating_p=float(F(r['min_generating_p'])),
                   median_observed=float(med))
        for k in sorted(th):
            row[f'theta_{k}'] = float(th[k])
            row[f'theta_{k}_rel'] = float(th[k] / med)
        if m0 is not None and m0 >= 2:
            first = th.get(m0 - 1, F(0))
            row['first_drop_eps'] = float(first)
            row['first_drop_eps_rel'] = float(first / med)
            row['first_drop_from'] = m0
        rows.append(row)
    report['full_panel'] = rows

    have = [r for r in rows if 'first_drop_eps_rel' in r]
    report['threshold_distribution'] = dict(
        n=len(have),
        absolute=quantiles([r['first_drop_eps'] for r in have]),
        relative=quantiles([r['first_drop_eps_rel'] for r in have]),
        by_K={str(K): quantiles([r['first_drop_eps_rel'] for r in have if r['K'] == K])
              for K in sorted({r['K'] for r in have})},
        by_family={f: quantiles([r['first_drop_eps_rel'] for r in have if r['family'] == f])
                   for f in sorted({r['family'] for r in have})},
        by_N_bin={b: quantiles([r['first_drop_eps_rel'] for r in have
                                if b == ('20-49' if r['N'] < 50 else '50-129' if r['N'] < 130
                                         else '130+')])
                  for b in ('20-49', '50-129', '130+')})
    def corr(xs, ys):
        if len(xs) < 4:
            return None
        rho, p = spearmanr(xs, ys)
        return dict(spearman=round(float(rho), 4), p=float(p), n=len(xs))
    report['threshold_scaling'] = dict(
        vs_K=corr([r['K'] for r in have], [r['first_drop_eps_rel'] for r in have]),
        vs_N=corr([r['N'] for r in have], [r['first_drop_eps_rel'] for r in have]),
        vs_min_p=corr([r['min_generating_p'] for r in have],
                      [r['first_drop_eps_rel'] for r in have]),
        vs_N_within_K3=corr([r['N'] for r in have if r['K'] == 3],
                            [r['first_drop_eps_rel'] for r in have if r['K'] == 3]),
        vs_min_p_within_K3=corr([r['min_generating_p'] for r in have if r['K'] == 3],
                                [r['first_drop_eps_rel'] for r in have if r['K'] == 3]))

    # ----------------------------- 3. within-eps incompleteness (exact + ladder)
    byname = defaultdict(list)
    for r in ok:
        if r['level'] != 'full':
            byname[r['name']].append(r)
    real_eps = {}
    for row in calib['rows']:
        if row['label'] == 'krepp_ani':
            real_eps[row['case']] = dict(krepp_median=F(row['residual_median']),
                                         krepp_max=F(row['residual_max']))
    ladder_rows = []
    per_draw = []
    for name, draws in sorted(byname.items()):
        base = full.get(name)
        if base is None:
            continue
        thf, K = thetas_of(base), base['K']
        m0f, med = base['result']['m_at_zero'], F(base['median_observed'])
        for d in draws:
            thd, m0d = thetas_of(d), d['result']['m_at_zero']
            intervals = []
            for k in range(1, K + 1):
                lo, hi = thd.get(k, F(0)), thf.get(k, F(0))
                if lo < hi:
                    intervals.append((lo, hi))
            sup = max((hi for _, hi in intervals), default=F(0))
            per_draw.append(dict(
                name=name, K=K, N=base['N'], level=d['level'], regime=d['regime'],
                replicate=d['replicate'], n_retained=d['n_retained'],
                falls_somewhere=bool(intervals),
                fall_sup_rel=float(sup / med) if intervals else 0.0,
                m0_full=m0f, m0_draw=m0d,
                falls_at_zero=(m0d is not None and m0f is not None and m0d < m0f),
                rises_at_zero=(m0d is not None and m0f is not None and m0d > m0f)))
            steps = [(s, F(Decimal(s)) * med) for s in LADDER]
            steps += [(lbl, eps) for lbl, eps in real_eps.get(name, {}).items()]
            for lbl, eps in steps:
                mf = m_of_eps(thf, K, eps, m0f)
                md = m_of_eps(thd, K, eps, m0d)
                ladder_rows.append(dict(name=name, K=K, level=d['level'], regime=d['regime'],
                                        eps_rel=lbl, m_full=mf, m_draw=md,
                                        fall=md < mf, rise=md > mf))
    report['per_draw'] = per_draw

    table = []
    for s in list(LADDER) + ['krepp_median', 'krepp_max']:
        sel = [x for x in ladder_rows if x['eps_rel'] == s]
        live = [x for x in sel if x['m_full'] > 1]
        entry = dict(eps_rel=s, n_draws=len(sel),
                     mean_m_full=round(sum(x['m_full'] for x in sel) / max(1, len(sel)), 4),
                     mean_m_draw=round(sum(x['m_draw'] for x in sel) / max(1, len(sel)), 4),
                     fall_rate=round(sum(x['fall'] for x in sel) / max(1, len(sel)), 4),
                     rise_rate=round(sum(x['rise'] for x in sel) / max(1, len(sel)), 4),
                     n_live=len(live),
                     fall_rate_live=round(sum(x['fall'] for x in live) / max(1, len(live)), 4),
                     share_m_full_is_1=round(1 - len(live) / max(1, len(sel)), 4))
        for regime in ('random', 'clade', 'nearest'):
            g = [x for x in sel if x['regime'] == regime]
            gl = [x for x in live if x['regime'] == regime]
            entry[f'fall_{regime}'] = round(sum(x['fall'] for x in g) / max(1, len(g)), 4)
            entry[f'fall_live_{regime}'] = round(sum(x['fall'] for x in gl) / max(1, len(gl)), 4)
        base = entry['fall_random'] or None
        entry['ratio_clade_random'] = round(entry['fall_clade'] / base, 3) if base else None
        entry['ratio_nearest_random'] = round(entry['fall_nearest'] / base, 3) if base else None
        bl = entry['fall_live_random'] or None
        entry['ratio_live_clade_random'] = round(entry['fall_live_clade'] / bl, 3) if bl else None
        entry['ratio_live_nearest_random'] = round(entry['fall_live_nearest'] / bl, 3) if bl else None
        table.append(entry)
    report['within_eps_table'] = table

    by_level = []
    for level in ('one', '0.10', '0.25', '0.50', '0.75'):
        for s in LADDER:
            sel = [x for x in ladder_rows if x['eps_rel'] == s and x['level'] == level]
            if not sel:
                continue
            by_level.append(dict(level=level, eps_rel=s, n=len(sel),
                                 fall_rate=round(sum(x['fall'] for x in sel) / len(sel), 4)))
    report['within_eps_by_level'] = by_level
    report['fall_region'] = dict(
        n_draws=len(per_draw),
        share_falls_somewhere=round(sum(x['falls_somewhere'] for x in per_draw)
                                    / max(1, len(per_draw)), 4),
        share_falls_at_zero=round(sum(x['falls_at_zero'] for x in per_draw)
                                  / max(1, len(per_draw)), 4),
        share_rises_at_zero=round(sum(x['rises_at_zero'] for x in per_draw)
                                  / max(1, len(per_draw)), 4),
        sup_relative=quantiles([x['fall_sup_rel'] for x in per_draw if x['falls_somewhere']]),
        by_regime={rg: dict(
            n=sum(1 for x in per_draw if x['regime'] == rg),
            share_falls_somewhere=round(
                sum(x['falls_somewhere'] for x in per_draw if x['regime'] == rg)
                / max(1, sum(1 for x in per_draw if x['regime'] == rg)), 4),
            sup_relative=quantiles([x['fall_sup_rel'] for x in per_draw
                                    if x['regime'] == rg and x['falls_somewhere']]))
            for rg in ('random', 'clade', 'nearest')})

    # --------------------------------------------------- 4. reality calibration
    cal = defaultdict(dict)
    for row in calib['rows']:
        cal[row['case']][row['label']] = row
    placed = []
    for r in have:
        entry = dict(name=r['name'], K=r['K'], N=r['N'],
                     first_drop_eps=r['first_drop_eps'],
                     first_drop_eps_rel=r['first_drop_eps_rel'],
                     theta_1_rel=r.get('theta_1_rel'))
        for label, row in sorted(cal.get(r['name'], {}).items()):
            ratio_max = float(F(row['residual_max']) / F(str(F(r['first_drop_eps']))))  \
                if r['first_drop_eps'] > 0 else None
            entry[label] = dict(residual_median=float(F(row['residual_median'])),
                                residual_max=float(F(row['residual_max'])),
                                residual_median_rel=row['residual_median_rel'],
                                residual_max_rel=row['residual_max_rel'],
                                median_over_first_drop=(float(F(row['residual_median'])
                                                        / F(r['first_drop_eps']))
                                                        if r['first_drop_eps'] > 0 else None),
                                max_over_first_drop=ratio_max)
        placed.append(entry)
    for e in placed:
        base = full[e['name']]
        th, K, m0 = thetas_of(base), base['K'], base['result']['m_at_zero']
        for label in ('krepp_ani', 'published_noise1_l10000', 'published_noise2_l10000'):
            if label in e:
                e[label]['m_at_median_residual'] = m_of_eps(
                    th, K, F(str(e[label]['residual_median'])), m0)
                e[label]['m_at_max_residual'] = m_of_eps(
                    th, K, F(str(e[label]['residual_max'])), m0)
                e[label]['m_at_zero'] = m0
    report['calibration_placed'] = placed
    for label in ('krepp_ani', 'published_noise1_l10000', 'published_noise2_l10000'):
        sel = [e for e in placed if label in e and 'm_at_max_residual' in e[label]]
        if sel:
            report.setdefault('m_at_real_error', {})[label] = dict(
                n=len(sel),
                m_at_zero=dict(Counter(e[label]['m_at_zero'] for e in sel)),
                m_at_median_residual=dict(Counter(e[label]['m_at_median_residual'] for e in sel)),
                m_at_max_residual=dict(Counter(e[label]['m_at_max_residual'] for e in sel)))
    for label in ('krepp_ani', 'published_noise0', 'published_noise1_l100',
                  'published_noise1_l10000', 'published_noise2_l100',
                  'published_noise2_l10000'):
        vals = [e[label]['median_over_first_drop'] for e in placed
                if label in e and e[label]['median_over_first_drop'] is not None]
        maxs = [e[label]['max_over_first_drop'] for e in placed
                if label in e and e[label]['max_over_first_drop'] is not None]
        if vals:
            report.setdefault('calibration_summary', {})[label] = dict(
                n=len(vals), median_residual_over_first_drop=quantiles(vals),
                max_residual_over_first_drop=quantiles(maxs),
                share_median_above_first_drop=round(sum(v > 1 for v in vals) / len(vals), 4),
                share_max_above_first_drop=round(sum(v > 1 for v in maxs) / len(maxs), 4))

    # ------------------------------------------------------- 5. probe + ceiling
    probe_rows = []
    for r in probe:
        if r['result'].get('status') != 'complete':
            probe_rows.append(dict(name=r['name'], K=r['K'], N=r['N'],
                                   status=r['result'].get('status')))
            continue
        th = thetas_of(r)
        med = F(r['median_observed'])
        probe_rows.append(dict(name=r['name'], family=r['family'], K=r['K'], N=r['N'],
                               n_edges=r['n_edges'], k_swept=sorted(th),
                               m_at_zero=r['result'].get('m_at_zero'),
                               m_at_zero_bounds=r['result'].get('m_at_zero_bounds'),
                               **{f'theta_{k}_rel': float(th[k] / med) for k in sorted(th)},
                               seconds=r['result'].get('sweep_seconds')))
    report['probe'] = probe_rows
    visited = defaultdict(int)
    seconds = 0.0
    maxsupports = 0
    for r in cross + probe:
        res = r['result']
        for k, v in (res.get('visited') or {}).items():
            visited[int(k)] += v
            maxsupports = max(maxsupports, v)
        seconds += res.get('process_wall_seconds', 0.0)
    report['ceiling'] = dict(
        supports_enumerated={str(k): v for k, v in sorted(visited.items())},
        total_supports=sum(visited.values()),
        largest_single_enumeration=maxsupports,
        max_k_reached=max(visited) if visited else None,
        max_edges=max((r['n_edges'] for r in cross + probe), default=None),
        max_reference_set=max((r['n_retained'] for r in cross + probe), default=None),
        total_child_wall_seconds=round(seconds, 1))

    (OUT / 'stage296_analysis.json').write_text(json.dumps(serial(report), indent=2) + '\n')
    print(json.dumps(dict(
        cross=report['counts'], eps0=dict(agree=agree, disagree=disagree, censored=censored),
        witness_flags=report['witness_flag_count'],
        thresholds=report['threshold_distribution']['relative'],
        fall_region=report['fall_region']['share_falls_somewhere'],
        ceiling=report['ceiling']), indent=2))


if __name__ == '__main__':
    main()
