#!/usr/bin/env python3
"""Stage297 analysis: validation agreement, the K = 10 threshold table and the
revised flatness correlations.

`describe`, `q`, `ranks` and `spearman` are IMPORTED from stage292_analysis so
the K = 10 numbers are computed exactly the way the K in {2,3,5,7} numbers were.
"""
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
sys.dont_write_bytecode = True

from stage292_analysis import describe, q, ranks, spearman                 # noqa: E402

OUTDIR = ROOT / 'results/stage297'
OUT = OUTDIR / 'stage297_analysis.json'


def load(path):
    p = Path(path)
    if not p.exists():
        return []
    return [json.loads(l) for l in p.open()]


def agreement_block(rows, modes, extra_keys=()):
    ok = [r for r in rows if 'row_status' not in r]
    bad = [r for r in rows if 'row_status' in r]
    out = dict(rows=len(ok), failed_jobs=len(bad),
               failed_detail=dict(Counter(r.get('row_status') for r in bad)),
               s292_status=dict(Counter(r['s292_status'] for r in ok)))
    for mode in modes:
        for field in ('minimum', 'supports', 'status'):
            key = 'agree_%s_%s' % (field, mode)
            vals = [r[key] for r in ok if key in r]
            out['%s|agree_%s' % (mode, field)] = dict(n=len(vals), ok=sum(vals))
        secs = [r[mode + '_seconds'] for r in ok if mode + '_seconds' in r]
        ref = [r['s292_seconds'] for r in ok if mode + '_seconds' in r]
        out['%s|seconds' % mode] = describe(secs)
        if secs:
            out['%s|speedup_vs_s292' % mode] = describe(
                [a / max(b, 1e-12) for a, b in zip(ref, secs)])
    mp = [r['agree_minimum_pass'] for r in ok if 'agree_minimum_pass' in r]
    out['reduced_minimum_pass_equals_frozen'] = dict(n=len(mp), ok=sum(mp))
    for key in extra_keys:
        vals = [r[key] for r in ok if key in r and r[key] is not None]
        out[key] = dict(n=len(vals), ok=sum(bool(v) for v in vals))
    out['disagreements'] = [
        {k: v for k, v in r.items()
         if k in ('track', 'name', 'seed', 'K', 'regime', 'e', 'R_remaining',
                  's292_minimum', 's292_n', 'combined_minimum', 'combined_n',
                  'recorded_minimum', 'recorded_n_supports', 'brute_minimum',
                  'brute_n') or k.startswith('agree_')}
        for r in ok if 'disagreement' in r]
    out['n_disagreements'] = len(out['disagreements'])
    return out


def validation():
    A = load(OUTDIR / 'stage297_validate_replay.jsonl')
    B = load(OUTDIR / 'stage297_validate_fresh.jsonl')
    C = load(OUTDIR / 'stage297_validate_baseline.jsonl')
    modes = ('combined', 'reduced_two_pass', 'frozen')
    out = dict(
        track_A=dict(agreement_block(A, modes,
                                     ('agree_recorded_s292', 'agree_recorded_combined')),
                     cases=len({r['name'] for r in A if 'row_status' not in r}),
                     K_distribution=dict(Counter(r['K'] for r in A if 'row_status' not in r))),
        track_B=dict(agreement_block(B, modes + ('combined_full_mask',),
                                     ('agree_brute_s292', 'agree_brute_combined',
                                      'agree_brute_reduced_two_pass',
                                      'agree_brute_frozen',
                                      'agree_brute_combined_full_mask')),
                     K_distribution=dict(Counter(r['K'] for r in B if 'row_status' not in r)),
                     brute_subsets_tested=sum(r.get('brute_visited', 0) for r in B
                                              if 'row_status' not in r)),
        track_C=dict(agreement_block(C, modes,
                                     ('agree_recorded_s292', 'agree_recorded_combined',
                                      'agree_agree_face_count',
                                      'agree_agree_pattern_count',
                                      'agree_agree_reduced_patterns')),
                     K_distribution=dict(Counter(r['K'] for r in C if 'row_status' not in r))))
    faces = [r['face_check'] for r in C if r.get('face_check')]
    if faces:
        ratio = [f['faces_297_full'] / max(f['faces_297_reduced'], 1) for f in faces]
        out['face_reduction'] = dict(
            checks=len(faces),
            full_equals_292=sum(f['agree_face_count'] for f in faces),
            pattern_count_equal=sum(f['agree_pattern_count'] for f in faces),
            reduced_patterns_equal_full=sum(f['agree_reduced_patterns'] for f in faces),
            faces_full=describe([f['faces_297_full'] for f in faces]),
            faces_reduced=describe([f['faces_297_reduced'] for f in faces]),
            reduction_factor=describe(ratio))
    return out


def bridge_rows(paths):
    rows = []
    for p in paths:
        rows.extend(load(p))
    return rows


def threshold_table(good, label):
    table, frac, sizes, status = (defaultdict(list), defaultdict(list),
                                  defaultdict(list), Counter())
    excluded = Counter()
    # one case can appear in more than one batch (a case that timed out in the
    # five-regime batch is re-measured in a regime-restricted one), so a
    # (case, regime, seed) measurement is taken once; repeats are counted and
    # checked for agreement instead of being pooled twice.
    seen, repeats, repeats_agree = {}, 0, 0
    for r in good:
        if not r['baseline_unique_true']:
            excluded[('baseline_not_unique_true', r['K'])] += 1
            continue
        for g in r['regimes']:
            ident = (r['name'], g['regime'], g['seed'])
            if ident in seen:
                repeats += 1
                repeats_agree += int(seen[ident] == g['t_unknown_sequence'])
                continue
            seen[ident] = g['t_unknown_sequence']
            status[(g['regime'], g['threshold_status'].split(':')[0])] += 1
            t = g['t_unknown_sequence']
            if t is None:
                continue
            key = (r['K'], g['regime'])
            table[key].append(t)
            frac[key].append(t / r['R_size'])
            sizes[key].append((r['R_size'], t))
    out = dict(label=label, cases=len(good),
               minimum_equals_K=sum(bool(r.get('baseline_minimum_equals_K',
                                               r['baseline_minimum'] == r['K']))
                                    for r in good),
               distinct_cases=len({r['name'] for r in good}),
               repeated_measurements=repeats,
               repeated_measurements_agree=repeats_agree,
               cases_unique_true=sum(r['baseline_unique_true'] for r in good),
               distinct_cases_unique_true=len({r['name'] for r in good
                                               if r['baseline_unique_true']}),
               excluded={'%s|K%d' % k: v for k, v in sorted(excluded.items())},
               threshold_status={'%s|%s' % k: v for k, v in sorted(status.items())},
               t_unknown={'K%d|%s' % k: describe(v) for k, v in sorted(table.items())},
               fraction_erasable={'K%d|%s' % k: describe(v) for k, v in sorted(frac.items())},
               spearman_t_vs_panel_size={'K%d|%s' % k: spearman([a for a, _ in v],
                                                               [b for _, b in v])
                                         for k, v in sorted(sizes.items()) if len(v) >= 3})
    byreg = defaultdict(list)
    for (K, reg), v in sizes.items():
        byreg[reg].extend(v)
    out['spearman_t_vs_panel_size_pooled'] = {
        reg: spearman([a for a, _ in v], [b for _, b in v])
        for reg, v in byreg.items() if len(v) >= 3}
    out['n_pairs_pooled'] = {reg: len(v) for reg, v in byreg.items()}
    allt, allf = defaultdict(list), defaultdict(list)
    for (K, reg), v in table.items():
        allt[reg].extend(v)
    for (K, reg), v in frac.items():
        allf[reg].extend(v)
    out['t_unknown_all_K'] = {reg: describe(v) for reg, v in sorted(allt.items())}
    out['fraction_erasable_all_K'] = {reg: describe(v) for reg, v in sorted(allf.items())}
    return out


def ceiling(good, failed):
    """What the measured counters say when a case or probe does not finish."""
    probe_status = Counter()
    probe_seconds = defaultdict(list)
    for r in good:
        bs = r.get('baseline', {})
        probe_status[('baseline', bs.get('status'))] += 1
        probe_seconds[('baseline', bs.get('status'))].append(r.get('baseline_seconds'))
        for g in r.get('regimes', []):
            for p in g['probe_table']:
                probe_status[(g['regime'], p['status'])] += 1
                probe_seconds[(g['regime'], p['status'])].append(p['seconds'])
    counters = defaultdict(list)
    for r in good:
        for g in r.get('regimes', []):
            for p in g['probe_table']:
                if p['status'] != 'complete':
                    for key in ('deficit_sets', 'faces_to_lp', 'masks_enumerated',
                                'masks_skipped_forced', 'oracle_calls', 'lift_products'):
                        counters[p['status'] + '|' + key].append(p[key])
    return dict(
        probe_status={'%s|%s' % k: v for k, v in sorted(probe_status.items(),
                                                        key=lambda z: str(z[0]))},
        probe_seconds={'%s|%s' % k: describe([x for x in v if x is not None])
                       for k, v in sorted(probe_seconds.items(), key=lambda z: str(z[0]))},
        incomplete_probe_counters={k: describe(v) for k, v in sorted(counters.items())},
        failed_jobs=len(failed),
        failed_detail=dict(Counter(r.get('row_status') or r.get('status') for r in failed)),
        failed_specs=[dict(status=r.get('row_status') or r.get('status'),
                           name=(r.get('spec') or {}).get('entry', {}).get('name'),
                           R_size=(r.get('spec') or {}).get('entry', {}).get('reference_set_size'),
                           wall_seconds=r.get('wall_seconds'))
                      for r in failed])


def k10():
    paths = sorted(OUTDIR.glob('stage297_bridge_*.jsonl'))
    rows = bridge_rows(paths)
    allgood = [r for r in rows if r.get('kind') == 'bridge']
    failed = [r for r in rows if r.get('kind') != 'bridge']
    # the pilot batch measured baselines only; baseline statistics are reported over
    # the distinct cases that carry a threshold search, with the pilot alongside
    good, pilot, seen_names = [], [], set()
    for r in sorted(allgood, key=lambda z: (not z.get('regimes'), z['name'])):
        if not r.get('regimes'):
            pilot.append(r)
            continue
        if r['name'] in seen_names:
            good.append(r)
            continue
        seen_names.add(r['name'])
        good.append(r)
    baseline_cases = {}
    for r in good:
        baseline_cases.setdefault(r['name'], r)
    good_unique = list(baseline_cases.values())
    out = dict(files=[p.name for p in paths], records=len(allgood),
               cases=len(good_unique), threshold_records=len(good),
               pilot_baseline_only=len(pilot), failed=len(failed),
               by_K=dict(Counter(r['K'] for r in good_unique)),
               by_R=dict(Counter(r['R_size'] for r in good_unique)),
               baseline_status=dict(Counter(r['baseline_status'] for r in good_unique)),
               baseline_minimum_equals_K=sum(r['baseline_minimum_equals_K'] for r in good_unique),
               baseline_unique_true=sum(r['baseline_unique_true'] for r in good_unique),
               baseline_minimum_distribution=dict(
                   Counter('K%d|min%s' % (r['K'], r['baseline_minimum']) for r in good_unique)),
               baseline_seconds=describe([r['baseline_seconds'] for r in good_unique]),
               baseline_supports_n=dict(Counter(r['baseline_n_supports'] for r in good_unique)))
    cross = [r['cross_check'] for r in allgood if r.get('cross_check')]
    if cross:
        out['mode_cross_check'] = dict(
            n=len(cross), agree_minimum=sum(c['agree_minimum'] for c in cross),
            agree_supports=sum(c['agree_supports'] for c in cross),
            combined_seconds=describe([r['baseline_seconds'] for r in allgood
                                       if r.get('cross_check')]),
            two_pass_seconds=describe([c['seconds'] for c in cross]))
    if good:
        out['thresholds_stage297'] = threshold_table(good, 'stage297 K>=10')
    out['pilot_baselines'] = dict(
        n=len(pilot),
        unique_true=sum(r['baseline_unique_true'] for r in pilot),
        minimum_equals_K=sum(r['baseline_minimum_equals_K'] for r in pilot),
        seconds=describe([r['baseline_seconds'] for r in pilot]))
    out['ceiling'] = ceiling(allgood, failed)
    return out


def stage292_reference():
    rows = []
    for f in ('stage292_bridge_K2_3_5_7.jsonl', 'stage292_bridge_K7_resume.jsonl'):
        rows.extend(load(ROOT / 'results' / f))
    good = [r for r in rows if r.get('kind') == 'bridge']
    return threshold_table(good, 'stage292 K in {2,3,5,7}')


def pooled(stage292_good, stage297_good):
    return threshold_table(stage292_good + stage297_good, 'pooled K in {2,3,5,7,10}')


def main():
    rows292 = []
    for f in ('stage292_bridge_K2_3_5_7.jsonl', 'stage292_bridge_K7_resume.jsonl'):
        rows292.extend(load(ROOT / 'results' / f))
    good292 = [r for r in rows292 if r.get('kind') == 'bridge']
    rows297 = bridge_rows(sorted(OUTDIR.glob('stage297_bridge_*.jsonl')))
    good297 = [r for r in rows297 if r.get('kind') == 'bridge' and r.get('regimes')]
    good297.sort(key=lambda z: z['name'])
    out = dict(validation=validation(), k10=k10(),
               stage292_reference=stage292_reference())
    if good297:
        out['pooled'] = pooled(good292, good297)
    OUT.write_text(json.dumps(out, indent=1, default=str))
    print('written', OUT)
    v = out['validation']
    for track in ('track_A', 'track_B', 'track_C'):
        b = v[track]
        print('%s rows=%d failed=%d disagreements=%d' %
              (track, b['rows'], b['failed_jobs'], b['n_disagreements']))
        for key in sorted(b):
            if '|agree_supports' in key or key.startswith('agree_') or \
                    key == 'reduced_minimum_pass_equals_frozen':
                print('   %-44s %s' % (key, b[key]))
    if 'face_reduction' in v:
        print('face_reduction', {k: v['face_reduction'][k] for k in
                                 ('checks', 'full_equals_292', 'pattern_count_equal',
                                  'reduced_patterns_equal_full')})
    k = out['k10']
    print('K>=10 cases=%d failed=%d unique_true=%d min==K=%d' %
          (k['cases'], k['failed'], k['baseline_unique_true'],
           k['baseline_minimum_equals_K']))
    print('baseline minima:', k['baseline_minimum_distribution'])
    print('baseline seconds:', k['baseline_seconds'])
    if 'mode_cross_check' in k:
        print('mode cross-check:', k['mode_cross_check'])
    if 'thresholds_stage297' in k:
        t = k['thresholds_stage297']
        print('t_unknown:', json.dumps(t['t_unknown'], default=str)[:800])
        print('spearman pooled:', t['spearman_t_vs_panel_size_pooled'])


if __name__ == '__main__':
    main()
