#!/usr/bin/env python3
"""Stage292 analysis: validation agreement, the extended threshold table and
the panel-size ladder."""
import json
import math
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/stage292_analysis.json'


def load(name):
    p = ROOT / 'results' / name
    if not p.exists():
        return []
    return [json.loads(l) for l in p.open()]


def ranks(xs):
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    out = [0.0] * len(xs)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and xs[order[j + 1]] == xs[order[i]]:
            j += 1
        r = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            out[order[k]] = r
        i = j + 1
    return out


def spearman(xs, ys):
    if len(xs) < 3:
        return None
    a, b = ranks(xs), ranks(ys)
    ma, mb = sum(a) / len(a), sum(b) / len(b)
    num = sum((x - ma) * (y - mb) for x, y in zip(a, b))
    da = math.sqrt(sum((x - ma) ** 2 for x in a))
    db = math.sqrt(sum((y - mb) ** 2 for y in b))
    if da == 0 or db == 0:
        return None
    return num / (da * db)


def q(xs, p):
    if not xs:
        return None
    s = sorted(xs)
    i = min(len(s) - 1, max(0, int(round(p * (len(s) - 1)))))
    return s[i]


def describe(xs):
    if not xs:
        return None
    return dict(n=len(xs), median=statistics.median(xs), mean=sum(xs) / len(xs),
                q25=q(xs, 0.25), q75=q(xs, 0.75), min=min(xs), max=max(xs))


def validation():
    A = [r for r in load('stage292_validate_replay.jsonl')]
    B = [r for r in load('stage292_validate_fresh.jsonl')]
    Aok = [r for r in A if 'row_status' not in r]
    Bok = [r for r in B if 'row_status' not in r]
    sp = [r['bf_seconds'] / max(r['new_seconds'], 1e-12) for r in Aok]
    out = dict(
        track_A=dict(
            probe_rows=len(Aok), failed_jobs=len(A) - len(Aok),
            cases=len({r['name'] for r in Aok}),
            agree_minimum=sum(r['agree_minimum'] for r in Aok),
            agree_support_sets=sum(r['agree_supports'] for r in Aok),
            compared_to_recorded=sum(1 for r in Aok if r.get('agree_recorded') is not None),
            agree_recorded=sum(1 for r in Aok if r.get('agree_recorded')),
            compared_to_exhaustive_exact=sum(1 for r in Aok
                                             if r.get('agree_exact_supports') is not None),
            agree_exhaustive_exact=sum(1 for r in Aok if r.get('agree_exact_supports')),
            disagreements=[{k: r[k] for k in
                            ('name', 'K', 'regime', 'seed', 'e', 'new_minimum',
                             'new_n', 'bf_minimum', 'bf_n', 'recorded_n_supports')}
                           for r in Aok if 'disagreement' in r],
            status_counts=dict(Counter(r['new_status'] for r in Aok)),
            speedup_vs_stage290_brute=describe(sp),
            new_seconds=describe([r['new_seconds'] for r in Aok]),
            brute_seconds=describe([r['bf_seconds'] for r in Aok]),
            n_supports_distribution=dict(Counter(r['new_n'] for r in Aok))),
        track_B=dict(
            cases=len(Bok), failed_jobs=len(B) - len(Bok),
            agree_minimum=sum(r['agree_minimum'] for r in Bok),
            agree_support_sets=sum(r['agree_supports'] for r in Bok),
            agree_reduced_vs_full_mask=sum(r['agree_full_mask'] for r in Bok),
            with_hidden_leaves=sum(1 for r in Bok if r['hidden'] > 0),
            K_distribution=dict(Counter(r['K'] for r in Bok)),
            brute_subsets_tested=sum(r['brute_visited'] for r in Bok),
            oracle_calls=sum(r['oracle_calls'] for r in Bok),
            multi_minimiser_cases=sum(1 for r in Bok if r['new_n'] > 1),
            disagreements=[{k: r[k] for k in ('seed', 'K', 'nleaf', 'R_size',
                                              'hidden', 'new_minimum', 'new_n',
                                              'brute_minimum', 'brute_n')}
                           for r in Bok if 'disagreement' in r]))
    return out


def bridge():
    rows = []
    for name in ('stage292_bridge_K2_3_5_7.jsonl', 'stage292_bridge_K7_resume.jsonl',
                 'stage292_bridge_K10_12.jsonl', 'stage292_bridge.jsonl'):
        rows.extend(load(name))
    good = [r for r in rows if r.get('kind') == 'bridge']
    failed = [r for r in rows if r.get('kind') != 'bridge']
    out = dict(cases=len(good), failed=len(failed),
               failed_detail=dict(Counter(r.get('status') for r in failed)),
               by_K=dict(Counter(r['K'] for r in good)),
               baseline_unique_true=sum(r['baseline_unique_true'] for r in good),
               baseline_minimum_equals_K=sum(r['baseline_minimum'] == r['K'] for r in good),
               rank_identity=sum(r['rank_identity'] for r in good),
               cogirth_tolerance_zero=sum(r['cogirth_tolerance'] == 0 for r in good))
    table = defaultdict(list)
    frac = defaultdict(list)
    sizes = defaultdict(list)
    status = Counter()
    bound_ok = bound_n = 0
    known_vs_unknown = defaultdict(lambda: ([], []))
    for r in good:
        if not r['baseline_unique_true']:
            continue
        for g in r['regimes']:
            status[(g['regime'], g['threshold_status'].split(':')[0])] += 1
            t = g['t_unknown_sequence']
            if t is None:
                continue
            key = (r['K'], g['regime'])
            table[key].append(t)
            frac[key].append(t / r['R_size'])
            sizes[key].append((r['R_size'], t))
            tk = g['t_known_sequence']
            if tk is not None:
                known_vs_unknown[g['regime']][0].append(tk)
                known_vs_unknown[g['regime']][1].append(t)
            if r['cogirth_tolerance'] is not None:
                bound_n += 1
                bound_ok += int(t >= r['cogirth_tolerance'])
    out['threshold_status'] = {'%s|%s' % k: v for k, v in sorted(status.items())}
    out['t_unknown'] = {'K%d|%s' % k: describe(v) for k, v in sorted(table.items())}
    out['fraction_erasable'] = {'K%d|%s' % k: describe(v) for k, v in sorted(frac.items())}
    out['spearman_t_vs_panel_size'] = {
        'K%d|%s' % k: spearman([a for a, _ in v], [b for _, b in v])
        for k, v in sorted(sizes.items()) if len(v) >= 3}
    byreg = defaultdict(list)
    for (K, reg), v in sizes.items():
        byreg[reg].extend(v)
    out['spearman_t_vs_panel_size_pooled'] = {
        reg: spearman([a for a, _ in v], [b for _, b in v]) for reg, v in byreg.items()}
    out['spearman_known_vs_unknown'] = {
        reg: spearman(a, b) for reg, (a, b) in known_vs_unknown.items() if len(a) >= 3}
    out['cogirth_bound_respected'] = dict(ok=bound_ok, n=bound_n)
    allt = defaultdict(list)
    for (K, reg), v in table.items():
        allt[reg].extend(v)
    out['t_unknown_all_K'] = {reg: describe(v) for reg, v in sorted(allt.items())}
    allf = defaultdict(list)
    for (K, reg), v in frac.items():
        allf[reg].extend(v)
    out['fraction_erasable_all_K'] = {reg: describe(v) for reg, v in sorted(allf.items())}
    out['seconds_per_case'] = describe([sum(g['seconds'] for g in r['regimes'])
                                        for r in good])
    return out


def panels():
    # every spec that failed in the first pass is re-queued in the re-run file,
    # so the first pass contributes only its successful rows
    rows = ([r for r in load('stage292_panels.jsonl') if r.get('kind') == 'panel']
            + load('stage292_panels_rerun.jsonl'))
    good = [r for r in rows if r.get('kind') == 'panel']
    failed = [r for r in rows if r.get('kind') != 'panel']
    out = dict(jobs=len(good), failed=len(failed),
               failed_detail=[dict(status=r.get('status'), spec=r.get('spec'))
                              for r in failed],
               baseline_unique_true=sum(r['baseline_unique_true'] for r in good),
               sizes=sorted({r['panel_leaves'] for r in good}))
    series = defaultdict(list)
    for r in good:
        if not r['baseline_unique_true']:
            continue
        panel = 'wolv2_backbone' if 'backbone_blen_conserved' in r['panel'] else 'wolv2_placement'
        for g in r['regimes']:
            if g['t_unknown'] is None:
                continue
            series[(panel, r['K'], g['regime'])].append(
                dict(N=r['panel_leaves'], t=g['t_unknown'], seed=r['seed'],
                     fraction=g['t_unknown'] / r['panel_leaves'],
                     status=g['threshold_status']))
    out['series'] = {'%s|K%d|%s' % k: sorted(v, key=lambda z: (z['N'], z['seed']))
                     for k, v in sorted(series.items())}
    out['spearman_t_vs_N'] = {}
    out['loglog_slope'] = {}
    for k, v in sorted(series.items()):
        N = [z['N'] for z in v]
        t = [z['t'] for z in v]
        out['spearman_t_vs_N']['%s|K%d|%s' % k] = spearman(N, t)
        pts = [(math.log(a), math.log(b)) for a, b in zip(N, t) if b > 0]
        if len(pts) >= 3:
            mx = sum(p[0] for p in pts) / len(pts)
            my = sum(p[1] for p in pts) / len(pts)
            den = sum((p[0] - mx) ** 2 for p in pts)
            out['loglog_slope']['%s|K%d|%s' % k] = (
                sum((p[0] - mx) * (p[1] - my) for p in pts) / den if den else None)
    byreg = defaultdict(list)
    for (panel, K, reg), v in series.items():
        byreg[reg].extend((z['N'], z['t']) for z in v)
    out['spearman_t_vs_N_pooled'] = {reg: spearman([a for a, _ in v], [b for _, b in v])
                                     for reg, v in byreg.items()}
    out['median_t_by_size'] = {}
    for reg in ('adversarial_nearest', 'random'):
        d = defaultdict(list)
        for (panel, K, r2), v in series.items():
            if r2 != reg:
                continue
            for z in v:
                d[z['N']].append(z['t'])
        out['median_t_by_size'][reg] = {str(n): describe(v) for n, v in sorted(d.items())}
    return out


def main():
    doc = dict(stage=292, validation=validation(), bridge=bridge(), panels=panels())
    OUT.write_text(json.dumps(doc, indent=2, default=str) + '\n')
    print(json.dumps(dict(validation=doc['validation']['track_A'],
                          trackB=doc['validation']['track_B']), indent=2,
                     default=str)[:2500])
    print('written', OUT)


if __name__ == '__main__':
    main()

