#!/usr/bin/env python3
"""Stage 301 step 5: analysis of every stage-301 record."""
import json
import statistics as st
import sys
from fractions import Fraction as F
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import stage301_common as s301   # noqa: E402

OUT = s301.OUT


def rows(name):
    p = OUT / name
    if not p.exists():
        return []
    return [json.loads(l) for l in p.read_text().splitlines() if l.strip()]


def q(xs, p):
    xs = sorted(xs)
    if not xs:
        return None
    i = min(len(xs) - 1, int(p * (len(xs) - 1) + 0.5))
    return xs[i]


def validation():
    recs = rows('stage301_validate_runs_main.jsonl') + rows('stage301_validate_runs_full.jsonl')
    A = [r for r in recs if r['arm'] == 'A']
    A2 = [r for r in recs if r['arm'] == 'A2']
    C = [r for r in recs if r['arm'] == 'C_count']
    B = [r for r in recs if r['arm'] == 'B_exhaustive']
    agree, disagree, pruned_agree, incomplete = 0, [], 0, 0
    prune_ratios, by_k = [], {}
    for r in A:
        a = r['result']
        if a.get('status') != 'complete':
            incomplete += 1
            continue
        got = a.get('theta')
        if got is not None and F(got) == F(r['recorded_theta']):
            agree += 1
            nc, ex = a.get('n_candidates'), a.get('exhaustive')
            if nc is not None and ex:
                prune_ratios.append(nc / ex)
                by_k.setdefault(r['k'], []).append(nc / ex)
                if nc < ex:
                    pruned_agree += 1
        else:
            disagree.append(dict(name=r['name'], k=r['k'], status=a.get('status'),
                                 recorded=r['recorded_theta'], got=str(got)))
    # Theta_k == 0 makes eps = Theta_k / 2 = 0, which is NOT strictly below
    # Theta_k; returning 0 there is correct, so those rows are excluded.
    A2nd = [r for r in A2 if F(r['recorded_theta']) != 0]
    viol2 = [dict(name=r['name'], k=r['k'], eps=r['eps'], got=str(r['result'].get('theta')))
             for r in A2nd if r['result'].get('theta') is not None
             and F(r['result']['theta']) <= F(r['eps'])]
    refuted2 = sum(1 for r in A2nd if r['result'].get('status') == 'refuted')
    cbad = [dict(name=r['name'], k=r['k'], visited=r['result']['visited'],
                 expect=r['result']['exhaustive'])
            for r in C if r['result']['visited'] != r['result']['exhaustive']]
    wit = [r for r in A if r['result'].get('independent_deviation') is not None]
    witbad = [r['name'] for r in wit
              if F(r['result']['independent_deviation']) != F(r['result']['theta'])]
    return dict(
        n_arm_A=len(A), n_complete=agree + len(disagree), n_agree=agree,
        n_disagree=len(disagree), disagreements=disagree[:20],
        n_incomplete=incomplete,
        n_agree_with_screen_active=pruned_agree,
        n_arm_A2=len(A2), n_arm_A2_nondegenerate=len(A2nd),
        n_A2_violations=len(viol2), A2_violations=viol2[:20],
        n_A2_refuted_outright=refuted2,
        n_A2_degenerate_zero_theta=len(A2) - len(A2nd),
        n_arm_C=len(C), n_C_violations=len(cbad), C_violations=cbad[:10],
        n_arm_B_exhaustive=len(B),
        n_witness_rechecked=len(wit), n_witness_mismatch=len(witbad),
        prune_ratio_median=(q(prune_ratios, 0.5) if prune_ratios else None),
        prune_ratio_min=(min(prune_ratios) if prune_ratios else None),
        n_vacuous=sum(1 for x in prune_ratios if x >= 1.0),
        n_measured=len(prune_ratios),
        prune_by_k={k: dict(n=len(v), median=q(v, 0.5), min=min(v),
                            vacuous=sum(1 for x in v if x >= 1.0))
                    for k, v in sorted(by_k.items())})


def reach():
    recs = rows('stage301_reach_runs.jsonl')
    out, best = [], {}
    sps = []
    for r in recs:
        a = r.get('result', {})
        if a.get('status') != 'complete':
            continue
        done = [g for g in a['rungs'] if g.get('status') == 'complete']
        if not done:
            continue
        top = max(g['k'] for g in done)
        row = dict(name=r['name'], K=r['K'], n_edges=r['n_edges'],
                   n_leaves=r['n_retained'], top_k=top,
                   stop=[g.get('status') for g in a['rungs'] if g.get('status') != 'complete'][:1],
                   eps_vacuous_pos=a['meta'].get('eps_vacuous_pos'),
                   median_edge_length=a.get('median_edge_length'),
                   rungs=[dict(k=g['k'], theta=g.get('theta'),
                               n_candidates=g.get('n_candidates'),
                               exhaustive=g.get('exhaustive'),
                               prune_ratio=g.get('prune_ratio'),
                               vacuous=g.get('vacuous'),
                               seconds=g.get('seconds'),
                               seconds_per_support=g.get('seconds_per_support'),
                               witness_agrees=g.get('witness_agrees'))
                          for g in done])
        for g in done:
            if g.get('seconds_per_support'):
                sps.append(g['seconds_per_support'])
            key = g['k']
            if key not in best or r['n_retained'] > best[key]['n_leaves']:
                best[key] = dict(n_leaves=r['n_retained'], n_edges=r['n_edges'],
                                 name=r['name'])
        out.append(row)
    nonvac = [g for row in out for g in row['rungs']
              if g['k'] >= 2 and not g['vacuous']]
    vac = [g for row in out for g in row['rungs'] if g['k'] >= 2 and g['vacuous']]
    return dict(n_cases=len(out),
                largest_leafset_by_k={k: v for k, v in sorted(best.items())},
                n_rungs_k_ge_2=len(nonvac) + len(vac),
                n_rungs_nonvacuous=len(nonvac), n_rungs_vacuous=len(vac),
                prune_ratio_nonvacuous_median=q([g['prune_ratio'] for g in nonvac], 0.5),
                prune_ratio_nonvacuous_min=(min(g['prune_ratio'] for g in nonvac)
                                            if nonvac else None),
                seconds_per_support_median=q(sps, 0.5),
                seconds_per_support_range=[min(sps), max(sps)] if sps else None,
                witness_mismatches=sum(1 for row in out for g in row['rungs']
                                       if g['witness_agrees'] is False),
                cases=out)


def td():
    recs = rows('stage301_td_runs.jsonl')
    tdr = [r for r in recs if r['arm'] == 'td' and r['result'].get('status') == 'complete']
    kr = [r for r in recs if r['arm'] == 'krepp' and r['result'].get('status') == 'complete']
    summary = json.loads((ROOT / 'results/stage298/stage298_td_summary.json').read_text())['rows']
    th1 = {row['name']: F(row['theta1']) for row in summary}
    med = {row['name']: F(row['median_observed']) for row in summary}
    tdout, soundness = [], []
    for r in tdr:
        a = r['result']
        name = r['name']
        b = {int(k): F(v['eps_refute']) for k, v in a['bounds'].items()}
        t1 = th1.get(name)
        if t1 is not None:
            soundness.append(dict(name=name, eps_R1=str(b[1]), theta1=str(t1),
                                  ok=bool(b[1] <= t1), ratio=float(b[1] / t1)))
        tdout.append(dict(name=name, n_edges=a['n_edges'], n_leaves=a['n_leaves'],
                          theta1=str(t1) if t1 else None,
                          median_observed=str(med.get(name)),
                          eps_refute={k: str(v) for k, v in b.items()},
                          eps_refute_over_median={k: float(v / med[name]) for k, v in b.items()}
                          if name in med else None,
                          theta1_over_median=(float(t1 / med[name]) if name in med else None),
                          max_leaf_envelope_slack=a['max_leaf_envelope_slack'],
                          seconds=a['seconds_view'] + a['seconds_bounds']))
    krout, kbad = [], []
    ratios = {}
    for r in kr:
        a = r['result']
        known = {int(k): F(v) for k, v in r['known_thetas'].items()}
        b = {int(k): F(v['eps_refute']) for k, v in a['bounds'].items()}
        for k, t in known.items():
            if k in b:
                ok = b[k] <= t
                if not ok:
                    kbad.append(dict(name=r['name'], k=k, eps_R=str(b[k]), theta=str(t)))
                if t > 0:
                    ratios.setdefault(k, []).append(float(b[k] / t))
        krout.append(dict(name=r['name'], n_leaves=a['n_leaves']))
    return dict(
        n_td=len(tdout), n_krepp=len(krout),
        td_soundness_checks=len(soundness),
        td_soundness_violations=sum(1 for s in soundness if not s['ok']),
        krepp_soundness_checks=sum(len(v) for v in ratios.values()),
        krepp_soundness_violations=len(kbad), krepp_violations=kbad[:10],
        tightness_eps_R_over_theta={k: dict(n=len(v), median=q(v, 0.5),
                                            p90=q(v, 0.9), max=max(v))
                                    for k, v in sorted(ratios.items())},
        td_eps_refute_median={k: q([F(r['eps_refute'][k]) for r in tdout
                                    if k in r['eps_refute']], 0.5)
                              for k in [str(i) for i in range(1, 9)]
                              if any(k in r['eps_refute'] for r in tdout)},
        td=tdout, td_soundness=soundness)


def grid():
    recs = rows('stage301_grid_runs_v2.jsonl')
    out = []
    for r in recs:
        a = r.get('result', {})
        if a.get('status') != 'complete':
            continue
        done = [g for g in a['rungs'] if g.get('status') == 'complete']
        out.append(dict(name=r['name'], K=r['K'], regime=r['regime'],
                        level=r['level'], replicate=r['replicate'],
                        n_full=r['n_full'], n_retained=r['n_retained'],
                        top_k=max([g['k'] for g in done], default=0),
                        thetas={g['k']: g.get('theta') for g in done},
                        vacuous={g['k']: g.get('vacuous') for g in done}))
    cross = {}
    for row in out:
        base = next((o for o in out if o['name'] == row['name']
                     and o['regime'] == 'none'), None)
        if base is None or row['regime'] == 'none':
            continue
        for k, t in row['thetas'].items():
            bt = base['thetas'].get(k)
            if bt is None or t is None:
                continue
            key = (row['regime'], row['level'], k)
            cross.setdefault(str(key), []).append(float(F(t) / F(bt)) if F(bt) else None)
    summ = {key: dict(n=len(v), median_ratio=q([x for x in v if x is not None], 0.5),
                      n_fall=sum(1 for x in v if x is not None and x < 1),
                      n_rise=sum(1 for x in v if x is not None and x > 1))
            for key, v in sorted(cross.items())}
    return dict(n_draws=len(out),
                top_k_reached=max([r['top_k'] for r in out], default=0),
                k_above_3=sum(1 for r in out if r['top_k'] >= 4),
                cross=summ, draws=out)


def main():
    payload = dict(stage=301, validation=validation(), reach=reach(),
                   real_data=td(), grid=grid())
    (OUT / 'stage301_analysis.json').write_text(
        json.dumps(payload, indent=1, default=str) + '\n')
    brief = {k: (v if not isinstance(v, dict) else
                 {kk: vv for kk, vv in v.items()
                  if kk not in ('cases', 'td', 'draws', 'td_soundness')})
             for k, v in payload.items()}
    print(json.dumps(brief, indent=1, default=str))


if __name__ == '__main__':
    main()
