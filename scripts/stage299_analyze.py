#!/usr/bin/env python3
"""Stage 299 step 5: analysis.

Reads the frozen stage-299 run streams and produces, per arm:

  * the verification ledger (double-decided edges, witness agreement, the
    monotonicity CHECK - removal deletes LP rows, so a strict rise is a bug -
    and the cross-stage agreement with stage 298 at the full panel);
  * the regime x level table of the Theta_1 fall, which is the real-data
    analogue of stage 288's count-drop table;
  * the dose-response: how much reference loss is needed before Theta_1 drops
    below the measured estimator residual, for two lines (worst leaf, median
    leaf) in two variants (recomputed on the retained leaves - primary - and
    fixed at the full panel);
  * attaining-edge migration, in exact tree distance and in hops.

All comparisons are made against the case's OWN full-panel Theta_1 in the SAME
arm.  Nothing is imputed; incomplete runs are counted and excluded.
"""
import csv
import json
import sys
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import stage288_common as common               # noqa: E402
import stage299_common as s299                 # noqa: E402

OUT = s299.OUT
LEVEL_ORDER = ['one', '0.10', '0.25', '0.50', '0.75']


def med(values):
    s = sorted(v for v in values if v is not None)
    if not s:
        return None
    return s[len(s) // 2] if len(s) % 2 else (s[len(s) // 2 - 1] + s[len(s) // 2]) / 2.0


def quart(values):
    s = sorted(v for v in values if v is not None)
    if not s:
        return (None, None, None)
    return (s[len(s) // 4], med(s), s[(3 * len(s)) // 4])


def frac(num, den):
    return None if not den else num / float(den)


def read_runs(path):
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def main():
    out_path = OUT / 'stage299_analysis.json'
    assert not out_path.exists(), 'Never overwrite frozen evidence'
    small_path = OUT / 'stage299_small_runs.jsonl'
    td_path = OUT / 'stage299_td_runs.jsonl'
    runs = read_runs(small_path) + read_runs(td_path)
    protocol = json.loads((OUT / 'stage299_protocol.json').read_text())
    case_index = {(c['arm'], c['name']): c for c in protocol['cases']}

    # ---- stage 298's full-panel Theta_1, for the cross-stage check
    s298_theta = {}
    for line in (s299.S298 / 'stage298_sweep_runs.jsonl').read_text().splitlines():
        if not line.strip():
            continue
        rec = json.loads(line)
        if rec['result'].get('k') == 1 and rec['result'].get('status') == 'complete':
            s298_theta[rec['name']] = F(rec['result']['theta'])
    for row in json.loads((s299.S298 / 'stage298_td_summary.json').read_text())['rows']:
        s298_theta[row['name']] = F(row['theta1'])

    cache = {}
    bank = s299.bank_by_name()

    def case_of(arm, name):
        if (arm, name) not in cache:
            row = case_index[(arm, name)]
            cache[(arm, name)] = s299.prepare_case(
                s299.load_json(row['path'], row['sha256']), name)
        return cache[(arm, name)]

    base, draws = {}, defaultdict(list)
    incomplete = []
    for r in runs:
        if r['result'].get('status') != 'complete':
            incomplete.append(dict(arm=r['arm'], name=r['name'], level=r['level'],
                                   regime=r['regime'], replicate=r['replicate'],
                                   status=r['result'].get('status'),
                                   reason=r['result'].get('reason')))
            continue
        key = (r['arm'], r['name'])
        if r['level'] == 'full':
            base[key] = r
        else:
            draws[key].append(r)

    # ---- verification ledger
    ver = dict(n_runs=len(runs), n_complete=len(runs) - len(incomplete),
               n_incomplete=len(incomplete), incomplete=incomplete[:50],
               double_decided_edges=0, closed_form_mismatches=0,
               closed_form_theta1_agree=0, closed_form_theta1_total=0,
               closed_form_x_agree=0, closed_form_x_total=0,
               witness_agree=0, witness_total=0,
               td_frozen_spot_checks=0, td_spot_mismatches=0)
    for r in runs:
        res = r['result']
        if res.get('status') != 'complete':
            continue
        if r['arm'] == s299.TD_ARM:
            ver['td_frozen_spot_checks'] += res.get('n_frozen_checks', 0)
            ver['td_spot_mismatches'] += res.get('closed_form_mismatches', 0)
            ver['closed_form_theta1_total'] += 1
            ver['closed_form_theta1_agree'] += int(bool(res.get('frozen_theta1_agrees')))
        else:
            ver['double_decided_edges'] += res.get('n_edges', 0)
            ver['closed_form_mismatches'] += res.get('closed_form_mismatches', 0)
            ver['closed_form_theta1_total'] += 1
            ver['closed_form_theta1_agree'] += int(bool(res.get('closed_form_theta1_agrees')))
            ver['closed_form_x_total'] += 1
            ver['closed_form_x_agree'] += int(bool(res.get('closed_form_x_agrees')))
        ver['witness_total'] += 1
        ver['witness_agree'] += int(bool(res.get('witness_agrees')))
    cross = []
    for (arm, name), r in sorted(base.items()):
        want = s298_theta.get(name)
        if want is None:
            continue
        cross.append(dict(arm=arm, name=name, agrees=(F(r['result']['theta1']) == want),
                          stage299=str(r['result']['theta1']), stage298=str(want)))
    ver['cross_stage298'] = dict(n=len(cross), agree=sum(c['agrees'] for c in cross),
                                 disagreements=[c for c in cross if not c['agrees']][:10])

    # ---- per-draw table
    per_draw, rises = [], []
    for (arm, name), rows in sorted(draws.items()):
        b = base.get((arm, name))
        if b is None:
            continue
        case = case_of(arm, name)
        edges = case['edges']
        deg = s299.degrees(edges)
        bres = b['result']
        theta_full = F(bres['theta1'])
        x_full = F(bres['argmin_x'] if arm != s299.TD_ARM else bres['argmin_x'])
        e_full = int(bres['argmin'])
        med_full = F(bres['median_observed'])
        # estimator residual, if this arm has an exact partner
        partner = case_index[(arm, name)].get('exact_partner')
        exact_obs = None
        if partner and arm in ('krepp_ani', 'krepp_placement'):
            pentry = bank[partner]
            exact_obs = s299.prepare_case(s299.load_json(pentry['path'], pentry['sha256']),
                                          partner)['observed']
        panel = sorted(case['observed'])
        if arm == 'exact_partner':
            outside = set(int(x) for x in (case_index[(arm, name)]['exact_only_references'] or []))
            panel = [r for r in panel if r not in outside]
        resid_full_max = resid_full_med = None
        if exact_obs is not None:
            devs = sorted(abs(case['observed'][r] - exact_obs[r]) for r in panel)
            resid_full_max, resid_full_med = devs[-1], devs[len(devs) // 2]
        for r in rows:
            res = r['result']
            theta = F(res['theta1'])
            if theta > theta_full:
                rises.append(dict(arm=arm, name=name, level=r['level'], regime=r['regime'],
                                  replicate=r['replicate'], full=str(theta_full),
                                  draw=str(theta)))
            # the travellers arm records no removed list (up to 7,931 ids per draw);
            # it is needed only for the residual lines, which that arm has not got
            removed = set(int(x) for x in (r.get('removed') or []))
            retained = ([x for x in panel if x not in removed]
                        if exact_obs is not None else None)
            e_draw = int(res['argmin'])
            x_draw = F(res['argmin_x'])
            mig = s299.point_distance(edges, e_full, x_full, e_draw, x_draw)
            row = dict(arm=arm, name=name, stem=r['stem'], family=r['family'],
                       K_generating=r.get('K_generating'), n_edges=r['n_edges'],
                       reference_set_size=r['reference_set_size'],
                       level=r['level'], regime=r['regime'], replicate=r['replicate'],
                       n_removed=r['n_removed'], n_retained=r['n_retained'],
                       removed_fraction=r['removed_fraction'],
                       theta1=float(theta), theta1_full=float(theta_full),
                       fall_ratio=(float(theta / theta_full) if theta_full else None),
                       strict_fall=bool(theta < theta_full),
                       theta1_over_median_observed=res.get('theta1_over_median_observed'),
                       migration=float(mig), migration_rel=float(mig / med_full),
                       hops=s299.edge_hops(edges, e_full, e_draw),
                       same_edge=bool(e_draw == e_full),
                       n_attaining=res.get('n_attaining'),
                       argmin_is_terminal=bool(res.get('argmin_is_terminal')),
                       full_argmin_is_terminal=bool(bres.get('argmin_is_terminal')))
            if exact_obs is not None:
                devs = sorted(abs(case['observed'][x] - exact_obs[x]) for x in retained)
                rmax, rmed = devs[-1], devs[len(devs) // 2]
                row.update(
                    residual_max_retained=float(rmax), residual_median_retained=float(rmed),
                    residual_max_full=float(resid_full_max),
                    residual_median_full=float(resid_full_med),
                    below_max_retained=bool(theta < rmax),
                    below_median_retained=bool(theta < rmed),
                    below_max_full=bool(theta < resid_full_max),
                    below_median_full=bool(theta < resid_full_med),
                    full_below_max_retained=bool(theta_full < resid_full_max),
                    full_below_median_retained=bool(theta_full < resid_full_med))
            per_draw.append(row)

    # ---- aggregate tables
    cells = defaultdict(list)
    for row in per_draw:
        cells[(row['arm'], row['level'], row['regime'])].append(row)
    table = []
    for (arm, level, regime), rs in sorted(cells.items(),
                                           key=lambda kv: (kv[0][0], LEVEL_ORDER.index(kv[0][1]),
                                                           kv[0][2])):
        q = quart([r['fall_ratio'] for r in rs])
        table.append(dict(arm=arm, level=level, regime=regime, n=len(rs),
                          n_cases=len({r['name'] for r in rs}),
                          removed_fraction_median=med([r['removed_fraction'] for r in rs]),
                          fall_ratio_q1=q[0], fall_ratio_median=q[1], fall_ratio_q3=q[2],
                          relative_drop_median=(None if q[1] is None else 1 - q[1]),
                          strict_fall_fraction=frac(sum(r['strict_fall'] for r in rs), len(rs)),
                          same_edge_fraction=frac(sum(r['same_edge'] for r in rs), len(rs)),
                          migration_rel_median=med([r['migration_rel'] for r in rs]),
                          hops_median=med([r['hops'] for r in rs]),
                          terminal_fraction=frac(sum(r['argmin_is_terminal'] for r in rs), len(rs)),
                          theta1_over_median_observed_median=med(
                              [r['theta1_over_median_observed'] for r in rs])))
    pooled = defaultdict(dict)
    for arm in sorted({r['arm'] for r in per_draw}):
        for regime in sorted({r['regime'] for r in per_draw if r['arm'] == arm}):
            rs = [r for r in per_draw if r['arm'] == arm and r['regime'] == regime]
            q = quart([r['fall_ratio'] for r in rs])
            pooled[arm][regime] = dict(
                n=len(rs), fall_ratio_median=q[1],
                relative_drop_median=(None if q[1] is None else 1 - q[1]),
                strict_fall_fraction=frac(sum(r['strict_fall'] for r in rs), len(rs)),
                migration_rel_median=med([r['migration_rel'] for r in rs]),
                same_edge_fraction=frac(sum(r['same_edge'] for r in rs), len(rs)),
                terminal_fraction=frac(sum(r['argmin_is_terminal'] for r in rs), len(rs)))
    # structured / random ratios, per arm and per level
    ratios = defaultdict(dict)
    for arm, by_regime in pooled.items():
        rnd = by_regime.get('random')
        if not rnd:
            continue
        for regime, v in by_regime.items():
            if regime == 'random':
                continue
            ratios[arm][regime] = dict(
                relative_drop_ratio=(None if not rnd['relative_drop_median'] else
                                     v['relative_drop_median'] / rnd['relative_drop_median']),
                strict_fall_ratio=(None if not rnd['strict_fall_fraction'] else
                                   v['strict_fall_fraction'] / rnd['strict_fall_fraction']),
                migration_ratio=(None if not rnd['migration_rel_median'] else
                                 v['migration_rel_median'] / rnd['migration_rel_median']))
    ratios_by_level = defaultdict(dict)
    for (arm, level, regime), rs in cells.items():
        base_rs = cells.get((arm, level, 'random'))
        if not base_rs or regime == 'random':
            continue
        b = 1 - quart([r['fall_ratio'] for r in base_rs])[1]
        v = 1 - quart([r['fall_ratio'] for r in rs])[1]
        ratios_by_level[arm].setdefault(regime, {})[level] = dict(
            relative_drop_ratio=(None if not b else v / b),
            strict_fall_ratio=frac(sum(r['strict_fall'] for r in rs) / float(len(rs)),
                                   sum(r['strict_fall'] for r in base_rs) / float(len(base_rs))
                                   if any(r['strict_fall'] for r in base_rs) else 0))

    # ---- dose-response
    dose = {}
    for arm in ('krepp_ani', 'krepp_placement'):
        rs_all = [r for r in per_draw if r['arm'] == arm and 'below_max_retained' in r]
        if not rs_all:
            continue
        dose[arm] = {}
        for variant, key, startkey in (
                ('worst_leaf_retained', 'below_max_retained', 'full_below_max_retained'),
                ('median_leaf_retained', 'below_median_retained', 'full_below_median_retained'),
                ('worst_leaf_full_panel', 'below_max_full', 'full_below_max_retained'),
                ('median_leaf_full_panel', 'below_median_full', 'full_below_median_retained')):
            already = sorted({r['name'] for r in rs_all if r[startkey]})
            live = [r for r in rs_all if r['name'] not in set(already)]
            per_level = {}
            for level in LEVEL_ORDER:
                for regime in sorted({r['regime'] for r in live}):
                    sub = [r for r in live if r['level'] == level and r['regime'] == regime]
                    if not sub:
                        continue
                    per_level.setdefault(regime, {})[level] = dict(
                        n=len(sub), n_below=sum(r[key] for r in sub),
                        fraction_below=frac(sum(r[key] for r in sub), len(sub)))
            firsts = defaultdict(list)
            for regime in sorted({r['regime'] for r in live}):
                for name in sorted({r['name'] for r in live}):
                    hit = None
                    for level in LEVEL_ORDER:
                        sub = [r for r in live if r['name'] == name and r['regime'] == regime
                               and r['level'] == level]
                        if sub and sum(r[key] for r in sub) * 2 >= len(sub):
                            hit = level
                            break
                    firsts[regime].append(dict(name=name, first_level=hit))
            dose[arm][variant] = dict(
                already_below_at_full=len(already),
                already_below_examples=already[:10],
                n_cases_in_dose_response=len({r['name'] for r in live}),
                by_regime_level=per_level,
                first_crossing_counts={
                    regime: {lvl: sum(1 for f in fs if f['first_level'] == lvl)
                             for lvl in LEVEL_ORDER + [None]}
                    for regime, fs in firsts.items()},
                never_crosses={regime: sum(1 for f in fs if f['first_level'] is None)
                               for regime, fs in firsts.items()})

    # ---- paired exact vs estimated, same leaf sets
    pairs = []
    byk = {(r['arm'], r['name'], r['level'], r['regime'], r['replicate']): r for r in per_draw}
    for r in per_draw:
        if r['arm'] != 'krepp_ani':
            continue
        partner = case_index[('krepp_ani', r['name'])]['exact_partner']
        e = byk.get(('exact_partner', partner, r['level'], r['regime'], r['replicate']))
        if e is None:
            continue
        pairs.append(dict(name=r['name'], level=r['level'], regime=r['regime'],
                          replicate=r['replicate'],
                          exact_fall=e['fall_ratio'], noisy_fall=r['fall_ratio'],
                          exact_strict=e['strict_fall'], noisy_strict=r['strict_fall'],
                          exact_migration_rel=e['migration_rel'],
                          noisy_migration_rel=r['migration_rel'],
                          exact_same_edge=e['same_edge'], noisy_same_edge=r['same_edge']))
    paired = dict(n=len(pairs),
                  both_fall=sum(1 for p in pairs if p['exact_strict'] and p['noisy_strict']),
                  exact_only=sum(1 for p in pairs if p['exact_strict'] and not p['noisy_strict']),
                  noisy_only=sum(1 for p in pairs if p['noisy_strict'] and not p['exact_strict']),
                  neither=sum(1 for p in pairs if not p['noisy_strict'] and not p['exact_strict']),
                  median_exact_fall=med([p['exact_fall'] for p in pairs]),
                  median_noisy_fall=med([p['noisy_fall'] for p in pairs]),
                  same_edge_exact=frac(sum(p['exact_same_edge'] for p in pairs), len(pairs)),
                  same_edge_noisy=frac(sum(p['noisy_same_edge'] for p in pairs), len(pairs)))

    # ---- full-panel census
    full_rows = []
    for (arm, name), r in sorted(base.items()):
        res = r['result']
        full_rows.append(dict(arm=arm, name=name, theta1=float(F(res['theta1'])),
                              theta1_over_median_observed=res.get('theta1_over_median_observed'),
                              argmin=res['argmin'], n_attaining=res.get('n_attaining'),
                              argmin_is_terminal=bool(res.get('argmin_is_terminal')),
                              n_leaves=res.get('n_leaves'), n_edges=res.get('n_edges')))
    full_census = {}
    for arm in sorted({r['arm'] for r in full_rows}):
        rs = [r for r in full_rows if r['arm'] == arm]
        full_census[arm] = dict(
            n=len(rs),
            theta1_over_median_observed_median=med([r['theta1_over_median_observed'] for r in rs]),
            terminal_fraction=frac(sum(r['argmin_is_terminal'] for r in rs), len(rs)),
            n_attaining_median=med([r['n_attaining'] for r in rs]),
            unique_attainer_fraction=frac(sum(1 for r in rs if r['n_attaining'] == 1), len(rs)))

    out = dict(
        stage=299, step='analysis',
        inputs=dict(small_runs_sha256=common.sha(small_path), td_runs_sha256=common.sha(td_path),
                    n_records=len(runs), protocol_digest=protocol['draw_plan']['digest']),
        verification=ver,
        monotonicity=dict(check='Theta_1 can only fall or stay; a strict rise is a bug',
                          n_draws=len(per_draw), rises=len(rises), examples=rises[:10]),
        full_panel=full_census, full_panel_rows=full_rows,
        by_arm_level_regime=table, pooled_by_arm_regime=pooled,
        structured_vs_random=ratios, structured_vs_random_by_level=ratios_by_level,
        dose_response=dose, paired_exact_vs_estimated=paired)
    out_path.write_text(json.dumps(common.serial(out), indent=1))

    fields = sorted({k for r in per_draw for k in r})
    with (OUT / 'stage299_per_draw.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for r in per_draw:
            writer.writerow(r)
    print('draws %d; rises %d; incomplete %d' % (len(per_draw), len(rises), len(incomplete)))
    print('cross-stage298 %d/%d; closed-form mismatches %d over %d double-decided edges; '
          'witnesses %d/%d' % (ver['cross_stage298']['agree'], ver['cross_stage298']['n'],
                               ver['closed_form_mismatches'], ver['double_decided_edges'],
                               ver['witness_agree'], ver['witness_total']))
    for arm in sorted(pooled):
        print(' %s:' % arm)
        for regime, v in sorted(pooled[arm].items()):
            print('   %-14s n=%5d fall %.4f  strict %.3f  migration/med %.4f  same-edge %.3f'
                  % (regime, v['n'], v['fall_ratio_median'] or -1, v['strict_fall_fraction'],
                     v['migration_rel_median'] or -1, v['same_edge_fraction']))
    print('wrote', out_path)


if __name__ == '__main__':
    main()

