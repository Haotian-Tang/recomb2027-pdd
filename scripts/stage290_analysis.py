#!/usr/bin/env python3
"""Stage290 analysis: aggregate the survey and the bridge, test the structural
law for the cogirth, and deliver the known-vs-unknown verdict."""
import json
import math
import os
import statistics as stat
import sys
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt                                   # noqa: E402
import numpy as np                                                # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
sys.dont_write_bytecode = True
from stage290_erasure_core import encode                          # noqa: E402

SURVEY = ROOT / 'results/stage290_survey.jsonl'
BRIDGE = ROOT / 'results/stage290_bridge.jsonl'
OUT = ROOT / 'results/stage290_analysis.json'
FIG1 = ROOT / 'results/stage290_fig_tolerance.png'
FIG2 = ROOT / 'results/stage290_fig_bridge.png'


def read(path):
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def frac(num, den):
    return None if den == 0 else num / den


def spearman(x, y):
    n = len(x)
    if n < 3:
        return None
    def rank(v):
        order = sorted(range(n), key=lambda i: v[i])
        r = [0.0] * n
        i = 0
        while i < n:
            j = i
            while j + 1 < n and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    rx, ry = rank(x), rank(y)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = math.sqrt(sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry))
    return None if den == 0 else num / den


def survey_report(rows):
    ok = [r for r in rows if r['status'] == 'PASS']
    bank = [r for r in ok if r['kind'] == 'bank']
    panel = [r for r in ok if r['kind'] == 'panel']
    rep = dict(n_rows=len(rows), n_pass=len(ok), n_fail=len(rows) - len(ok),
               failures=[{k: v for k, v in r.items() if k in ('spec', 'error_type', 'error')}
                         for r in rows if r['status'] != 'PASS'],
               n_bank=len(bank), n_panel=len(panel))

    checks = dict(
        rank_identity_bank=sum(r['rank_identity'] for r in bank),
        rank_identity_panel=sum(r['rank_identity'] for r in panel),
        reduction_checked=sum('reduction_agrees' in r for r in bank),
        reduction_agrees=sum(r.get('reduction_agrees') is True for r in bank),
        witness_checked=sum(r.get('witness_checked') is True for r in bank + panel),
        brute_complete=sum(r.get('brute_status') == 'complete' for r in bank),
        brute_agrees=sum(r.get('brute_agrees') is True for r in bank),
        brute_limited=sum(r.get('brute_status') == 'limit' for r in bank),
        brute_subsets_tested=sum(r.get('brute_tested') or 0 for r in bank))
    checks['brute_agreement_rate'] = frac(checks['brute_agrees'], checks['brute_complete'])
    rep['checks'] = checks

    tol = [r['tolerance'] for r in bank]
    rep['bank_tolerance'] = dict(
        distribution=dict(sorted(Counter(tol).items())),
        zero_fraction=frac(sum(t == 0 for t in tol), len(tol)),
        max=max(tol) if tol else None,
        by_family={f: dict(sorted(Counter(r['tolerance'] for r in bank
                                          if r['family'] == f).items()))
                   for f in sorted({r['family'] for r in bank})},
        by_K={k: dict(sorted(Counter(r['tolerance'] for r in bank
                                     if r['K_generating'] == k).items()))
              for k in sorted({r['K_generating'] for r in bank})},
        by_tree={t: dict(sorted(Counter(r['tolerance'] for r in bank
                                        if r['name'].split('_k')[0] == t).items()))
                 for t in sorted({r['name'].split('_k')[0] for r in bank})})

    # structural law: cogirth == smallest leaf split of any support edge
    law = [r for r in bank if r['s'] == r['P_size']]
    rep['structural_law'] = dict(
        n_full_rank_cut=len(law),
        n_cases=len(bank),
        holds=sum(r['delta'] == r['min_small_side'] for r in law),
        violations=[dict(name=r['name'], delta=r['delta'],
                         min_small_side=r['min_small_side'], s=r['s'],
                         P=r['P_size']) for r in law
                    if r['delta'] != r['min_small_side']][:20],
        upper_bound_holds_all=sum(r['delta'] <= r['min_small_side'] for r in bank),
        deficient_cut_cases=sum(r['s'] < r['P_size'] for r in bank))
    rep['structural_law']['rate'] = frac(rep['structural_law']['holds'], len(law))

    rep['pendant'] = dict(
        with_pendant_support=sum(r['any_pendant_support'] for r in bank),
        tolerance_zero_given_pendant=sum(r['tolerance'] == 0 for r in bank
                                         if r['any_pendant_support']),
        tolerance_zero_given_no_pendant=sum(r['tolerance'] == 0 for r in bank
                                            if not r['any_pendant_support']),
        n_no_pendant=sum(not r['any_pendant_support'] for r in bank))

    rep['fibre'] = dict(
        matching_support=sum(r['matching'] for r in bank),
        dim_zero=sum(r['fibre_dim_endpoint'] == 0 for r in bank),
        dim_distribution=dict(sorted(Counter(r['fibre_dim_endpoint'] for r in bank).items())),
        dim_2H_distribution=dict(sorted(Counter(r['fibre_dim_2H'] for r in bank).items())))

    rep['correlations'] = dict(
        tolerance_vs_panel_size=spearman([r['R_size'] for r in bank], tol),
        tolerance_vs_K=spearman([r['K_generating'] for r in bank], tol),
        tolerance_vs_min_small_side=spearman([r['min_small_side'] for r in bank], tol),
        s_vs_2H=spearman([2 * r['H'] for r in bank], [r['s'] for r in bank]))

    pr = defaultdict(list)
    for r in panel:
        pr[(r['panel_leaves'], r['K'], r['mode'])].append(r['tolerance'])
    rep['panels'] = [dict(leaves=k[0], K=k[1], mode=k[2], n=len(v),
                          tolerances=sorted(v), median=stat.median(v))
                     for k, v in sorted(pr.items())]
    rep['panel_summary'] = dict(
        any_mode_zero=frac(sum(r['tolerance'] == 0 for r in panel if r['mode'] == 'any'),
                           sum(r['mode'] == 'any' for r in panel)),
        internal_mode_zero=frac(sum(r['tolerance'] == 0 for r in panel if r['mode'] == 'internal'),
                                sum(r['mode'] == 'internal' for r in panel)),
        internal_tolerance=dict(sorted(Counter(r['tolerance'] for r in panel
                                               if r['mode'] == 'internal').items())),
        any_tolerance=dict(sorted(Counter(r['tolerance'] for r in panel
                                          if r['mode'] == 'any').items())),
        max_leaves=max((r['panel_leaves'] for r in panel), default=None))
    return rep, bank, panel


def bridge_report(rows):
    ok = [r for r in rows if r['status'] == 'PASS']
    rep = dict(n_rows=len(rows), n_pass=len(ok), n_fail=len(rows) - len(ok),
               failures=[{k: v for k, v in r.items() if k in ('spec', 'error_type', 'error')}
                         for r in rows if r['status'] != 'PASS'])
    gt = [r for r in ok if r['baseline_ground_truth'] is not None]
    rep['validation'] = dict(
        ground_truth_cases=len(gt),
        ground_truth_agrees=sum(r['baseline_ground_truth'] for r in gt),
        ground_truth_rate=frac(sum(r['baseline_ground_truth'] for r in gt), len(gt)),
        frozen_solver_agrees=sum(r['frozen_agrees'] for r in ok),
        rank_identity=sum(r['rank_identity'] for r in ok),
        baseline_unique_true=sum(r['baseline_unique_true'] for r in ok),
        prefilter_margin_min=min((r['prefilter_margin'] for r in ok
                                  if r['prefilter_margin'] is not None), default=None),
        prefilter_survivor_rate=frac(sum(r['prefilter_survivors'] for r in ok),
                                     sum(r['prefilter_candidates'] for r in ok)))
    pairs = []
    for r in ok:
        for g in r['regimes']:
            pairs.append(dict(name=r['name'], family=r['family'], K=r['K'],
                              R=r['R_size'], regime=g['regime'], seed=g['seed'],
                              cogirth_tolerance=r['cogirth_tolerance'],
                              t_known=g['t_known_sequence'],
                              t_unknown=g['t_unknown_sequence'],
                              status=g['threshold_status'], emax=g['emax'],
                              fibre_dim=r['fibre_dim'], s=r['s'],
                              pendant=r['n_pendant_support']))
    rep['n_measurements'] = len(pairs)
    by = defaultdict(list)
    for p in pairs:
        by[p['regime']].append(p)
    out = {}
    for regime, ps in sorted(by.items()):
        known = [p['t_known'] for p in ps if p['t_known'] is not None]
        unk = [p for p in ps if p['t_unknown'] is not None]
        both = [p for p in ps if p['t_unknown'] is not None and p['t_known'] is not None]
        out[regime] = dict(
            n=len(ps),
            censored_unknown=sum(p['t_unknown'] is None for p in ps),
            censored_known=sum(p['t_known'] is None for p in ps),
            t_known_median=stat.median(known) if known else None,
            t_known_distribution=dict(sorted(Counter(known).items())),
            t_unknown_median=stat.median([p['t_unknown'] for p in unk]) if unk else None,
            t_unknown_distribution=dict(sorted(Counter(p['t_unknown'] for p in unk).items())),
            unknown_ge_cogirth=sum(p['t_unknown'] >= p['cogirth_tolerance'] for p in unk),
            unknown_lt_cogirth=sum(p['t_unknown'] < p['cogirth_tolerance'] for p in unk),
            unknown_ge_known=sum(p['t_unknown'] >= p['t_known'] for p in both),
            unknown_lt_known=sum(p['t_unknown'] < p['t_known'] for p in both),
            unknown_eq_known=sum(p['t_unknown'] == p['t_known'] for p in both),
            n_both=len(both),
            gap_median=stat.median([p['t_unknown'] - p['t_known'] for p in both]) if both else None,
            gap_max=max([p['t_unknown'] - p['t_known'] for p in both], default=None),
            spearman_known_unknown=spearman([p['t_known'] for p in both],
                                            [p['t_unknown'] for p in both]),
            spearman_cogirth_unknown=spearman([p['cogirth_tolerance'] for p in unk],
                                              [p['t_unknown'] for p in unk]),
            spearman_R_unknown=spearman([p['R'] for p in unk],
                                        [p['t_unknown'] for p in unk]))
    rep['by_regime'] = out
    allb = [p for p in pairs if p['t_unknown'] is not None and p['t_known'] is not None]
    allu = [p for p in pairs if p['t_unknown'] is not None]
    rep['overall'] = dict(
        n_both=len(allb), n_unknown_measured=len(allu),
        unknown_ge_known=sum(p['t_unknown'] >= p['t_known'] for p in allb),
        unknown_lt_known=sum(p['t_unknown'] < p['t_known'] for p in allb),
        unknown_ge_cogirth=sum(p['t_unknown'] >= p['cogirth_tolerance'] for p in allu),
        violations_unknown_lt_cogirth=[p for p in allu
                                       if p['t_unknown'] < p['cogirth_tolerance']][:20],
        violations_unknown_lt_known=[p for p in allb
                                     if p['t_unknown'] < p['t_known']][:20],
        spearman_known_unknown=spearman([p['t_known'] for p in allb],
                                        [p['t_unknown'] for p in allb]),
        spearman_cogirth_unknown=spearman([p['cogirth_tolerance'] for p in allu],
                                          [p['t_unknown'] for p in allu]))
    return rep, pairs


def figures(bank, panel, pairs):
    try:
        sys.path.insert(0, str(ROOT / 'scripts'))
        from stage290_figstyle import apply_figure_style
        apply_figure_style()
    except Exception:
        pass
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.0))
    ax = axes[0]
    fams = sorted({r['family'] for r in bank})
    width = 0.8 / max(1, len(fams))
    tolmax = max(r['tolerance'] for r in bank)
    xs = np.arange(tolmax + 2)
    for i, f in enumerate(fams):
        cnt = Counter(r['tolerance'] for r in bank if r['family'] == f)
        ax.bar(xs + i * width - 0.4, [cnt.get(int(x), 0) for x in xs], width,
               label='%s (n=%d)' % (f, sum(cnt.values())))
    ax.set_xlabel('erasure tolerance  (cogirth - 1)')
    ax.set_ylabel('cases')
    ax.set_title('A  Bank: 893 exactly explainable cases')
    ax.set_xticks(xs)
    ax.legend(fontsize=7)

    ax = axes[1]
    for mode, colour in (('any', 'tab:blue'), ('internal', 'tab:orange')):
        pts = [(r['panel_leaves'], r['tolerance']) for r in panel if r['mode'] == mode]
        if not pts:
            continue
        jx = np.array([p[0] for p in pts], dtype=float)
        jy = np.array([p[1] for p in pts], dtype=float)
        jy = jy + np.random.default_rng(0).uniform(-0.12, 0.12, size=jy.shape)
        ax.scatter(jx, jy, s=14, alpha=0.6, color=colour,
                   label='support on %s edges' % mode)
    ax.set_xscale('log')
    ax.set_xlabel('panel size (leaves)')
    ax.set_ylabel('erasure tolerance')
    ax.set_title('B  Frozen panels, random supports')
    ax.legend(fontsize=7)

    ax = axes[2]
    markers = {'random': 'o', 'whole_clade': 's', 'adversarial_nearest': '^'}
    for regime, mk in markers.items():
        ps = [p for p in pairs if p['regime'] == regime
              and p['t_unknown'] is not None and p['t_known'] is not None]
        if not ps:
            continue
        ax.scatter([p['t_known'] for p in ps], [p['t_unknown'] for p in ps],
                   s=18, alpha=0.6, marker=mk, label=regime)
    lim = [0, max([p['t_unknown'] for p in pairs if p['t_unknown'] is not None] +
                  [p['t_known'] for p in pairs if p['t_known'] is not None] + [1]) + 1]
    ax.plot(lim, lim, 'k--', lw=0.8, label='equality')
    ax.set_xlabel('known-support tolerance along the sequence')
    ax.set_ylabel('unknown-support threshold')
    ax.set_title('C  Bridge')
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(FIG1, dpi=200)

    fig2, ax2 = plt.subplots(1, 2, figsize=(9.5, 4.0))
    a = ax2[0]
    regimes = sorted({p['regime'] for p in pairs})
    data = [[p['t_unknown'] for p in pairs if p['regime'] == r and p['t_unknown'] is not None]
            for r in regimes]
    a.boxplot([d if d else [0] for d in data], tick_labels=regimes)
    a.set_ylabel('unknown-support threshold (references erased)')
    a.set_title('A  Threshold by removal regime')
    a.tick_params(axis='x', rotation=20)
    b = ax2[1]
    ps = [p for p in pairs if p['t_unknown'] is not None]
    b.scatter([p['R'] for p in ps], [p['t_unknown'] for p in ps], s=16, alpha=0.6)
    b.set_xlabel('panel size |R|')
    b.set_ylabel('unknown-support threshold')
    b.set_title('B  Threshold vs panel size')
    fig2.tight_layout()
    fig2.savefig(FIG2, dpi=200)
    return [str(FIG1.relative_to(ROOT)), str(FIG2.relative_to(ROOT))]


def main():
    srows = read(SURVEY)
    brows = read(BRIDGE) if BRIDGE.exists() else []
    srep, bank, panel = survey_report(srows)
    brep, pairs = bridge_report(brows) if brows else ({}, [])
    figs = figures(bank, panel, pairs) if os.environ.get('STAGE290_FIGS') else []
    out = dict(stage=290, survey=srep, bridge=brep, figures=figs)
    with OUT.open('w') as fh:
        json.dump(encode(out), fh, indent=2)
        fh.write('\n')
    print(json.dumps(encode(dict(
        survey_checks=srep['checks'],
        tolerance=srep['bank_tolerance']['distribution'],
        zero_fraction=srep['bank_tolerance']['zero_fraction'],
        structural_law=srep['structural_law']['rate'],
        panels=srep['panel_summary'],
        bridge_overall=brep.get('overall'),
        bridge_validation=brep.get('validation'))), indent=2))


if __name__ == '__main__':
    main()
