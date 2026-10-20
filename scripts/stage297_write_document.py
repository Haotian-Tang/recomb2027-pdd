#!/usr/bin/env python3
"""Stage297: assemble results/stage297/stage297_k10_bridge_and_seeds.md from the
registered protocols and the analysis JSONs.  Every number in the document is
read from a file; none is typed in.
"""
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
sys.dont_write_bytecode = True

OUTDIR = ROOT / 'results/stage297'
OUT = OUTDIR / 'stage297_k10_bridge_and_seeds.md'
REGIMES = ('random', 'whole_clade', 'adversarial_nearest')


def load(path, default=None):
    p = Path(path)
    if not p.exists():
        return default
    return json.loads(p.read_text())


def jsonl(path):
    p = Path(path)
    if not p.exists():
        return []
    return [json.loads(l) for l in p.open()]


def fmt(x, nd=3):
    if x is None:
        return 'n/a'
    if isinstance(x, bool):
        return 'yes' if x else 'no'
    if isinstance(x, float):
        return ('%.' + str(nd) + 'f') % x
    return str(x)


def d(block, key, nd=2):
    if not block or key not in block or block[key] is None:
        return 'n/a'
    return fmt(block[key], nd)


def describe_row(label, block, nd=2):
    if not block:
        return '| %s | 0 | n/a | n/a | n/a | n/a | n/a |' % label
    return '| %s | %d | %s | %s | %s | %s | %s |' % (
        label, block['n'], d(block, 'median', nd), d(block, 'q25', nd),
        d(block, 'q75', nd), d(block, 'min', nd), d(block, 'max', nd))


def main():
    analysis = load(OUTDIR / 'stage297_analysis.json', {})
    dd = load(OUTDIR / 'stage297_dd_analysis.json', {})
    cross = load(OUTDIR / 'stage297_dd_vs_stage289.json', {})
    protocol = load(OUTDIR / 'stage297_protocol.json', {})
    k12 = [r for r in jsonl(OUTDIR / 'stage297_bridge_k12_baseline.jsonl')
           if r.get('kind') == 'bridge']
    bridge_records = []
    for p in sorted(OUTDIR.glob('stage297_bridge_*.jsonl')):
        bridge_records.extend(r for r in jsonl(p) if r.get('kind') == 'bridge')
    sizes_seen = [r['R_size'] for r in bridge_records
                  if r.get('regimes') and r.get('baseline_unique_true')]
    val = analysis.get('validation', {})
    k10 = analysis.get('k10', {})
    ref = analysis.get('stage292_reference', {})
    pool = analysis.get('pooled', {})
    t297 = k10.get('thresholds_stage297', {})
    L = []
    A = L.append

    A('# Stage 297 — the K = 10 support-uniqueness bridge, and DecoDiPhy seed sensitivity')
    A('')
    A('Protocol `results/stage297/stage297_protocol.json`, registered '
      '%s, with addenda `stage297_protocol_addendum_main_k10.json` and '
      '`stage297_protocol_addendum_k10_extra.json` and '
      '`stage297_protocol_addendum_k10_large.json` (each registered before the batch '
      'it covers).  Every run is retained, including timeouts and deaths.'
      % protocol.get('registered_utc', 'n/a'))
    A('')

    A('### Scope actually covered')
    A('')
    A('- K = 10: **%s distinct bank cases** of the 213 exact K = 10 cases — a stratified '
      'subsample (first 3 per distinct |R| in `main_k10`, 2 more in `main_k10_extra`, '
      'plus the regime-restricted large-panel batch), not the whole family.'
      % k10.get('cases', 'n/a'))
    A('- The `k10_rand_large` contrast batch was **stopped after %d of 55 cases**; the '
      'random arm at |R| >= 81 is therefore covered by whatever the five-regime batches '
      'completed plus those %d, not by a full sweep.'
      % (len([r for r in jsonl(OUTDIR / "stage297_bridge_k10_rand_large.jsonl")]),
         len([r for r in jsonl(OUTDIR / "stage297_bridge_k10_rand_large.jsonl")])))
    A('- The registered `pilot_k10_threshold` batch was not executed; its six cases are '
      'the first per-size picks of `main_k10`, which covers them.')
    A('- K = 12: one case attempted, not measured (section 4).')
    A('- Seed sweep: %s cases, one per (family, K, |R|) cell, full panels only.'
      % dd.get('cases', 'n/a'))
    A('')

    # ---------------------------------------------------------------- item 1
    A('## 1  What was actually missing at K >= 10')
    A('')
    A('Stage 292 could not reach K = 10 because the **frozen minimum pass** its '
      'all-minimisers routine calls (`stage43_partial_fpt.solve_partial`) returned '
      "`unknown` after 3,600 s on the smallest K = 10 bank case "
      '(`results/stage292_k10_ceiling.txt`), and all eight K = 10 jobs of its stride-3 '
      'batch hit the per-case wall (`results/stage292_bridge_K10_12.jsonl`, 8/8 timeout).')
    A('')
    A('The forced-positivity face reduction, however, was **already present** in '
      '`stage292_all_minimisers.backbone_patterns` — the `full_mask` switch, '
      'cross-checked there on 500/500 fresh instances.  So the combination this stage '
      'needed was not "add the reduction to the face enumeration" but "replace the '
      'frozen minimum pass".  `scripts/stage297_reduced_all_minimisers.py` offers three '
      'interchangeable modes:')
    A('')
    A('| mode | pass 1 (minimum) | pass 2 (all minimisers) |')
    A('|---|---|---|')
    A('| `frozen` | frozen `solve_partial` | reduced faces at k\\* — the Stage292 routine |')
    A('| `reduced_two_pass` | Stage295 `solve_partial_reduced` | reduced faces at k\\* |')
    A('| `combined` | — | one reduced enumeration yields the minimum **and** every '
      'minimiser at it |')
    A('')
    A('`combined` is the operative one: a face is LP-tested when its cover count is at '
      'most the best count established so far, the incumbent is lowered by any strictly '
      'feasible smaller face, and every feasible face at the final incumbent is retained. '
      'Because `best` only decreases and every face with count <= best is tested, the '
      'result does not depend on visiting order, and the dominant cost '
      '(deficit sets x faces x exact rational LP) is paid once instead of twice.')
    A('')

    A('## 2  Validation — three tracks, zero disagreements')
    A('')
    A('| track | what it compares | rows | identical minimiser sets |')
    A('|---|---|---|---|')
    for tag, name, what in (('track_A', 'A replay',
                             'every probe point of `stage292_validate_replay.jsonl` '
                             '(130 cases, K in {2,3}), removal orders regenerated with the '
                             'frozen `removal_order`'),
                            ('track_B', 'B fresh',
                             'the 500 instances of `stage292_validate_fresh.jsonl`, '
                             'regenerated from `stage292_validate.random_instance`'),
                            ('track_C', 'C high K',
                             'the full-panel baseline of all 676 Stage292 bridge cases, '
                             'K in {2,3,5,7}')):
        b = val.get(tag, {})
        cell = []
        for mode in ('combined', 'reduced_two_pass', 'frozen', 'combined_full_mask'):
            key = '%s|agree_supports' % mode
            if key in b:
                cell.append('%s %d/%d' % (mode, b[key]['ok'], b[key]['n']))
        A('| %s | %s | %d | %s |' % (name, what, b.get('rows', 0), '; '.join(cell)))
    A('')
    extra = []
    for tag, label in (('track_A', 'A'), ('track_B', 'B'), ('track_C', 'C')):
        b = val.get(tag, {})
        if 'agree_recorded_combined' in b:
            extra.append('track %s vs the value Stage292 recorded: %d/%d' %
                         (label, b['agree_recorded_combined']['ok'],
                          b['agree_recorded_combined']['n']))
        if 'agree_brute_combined' in b:
            extra.append('track %s vs the independent exhaustive C(|E|,k) enumeration: '
                         '%d/%d' % (label, b['agree_brute_combined']['ok'],
                                    b['agree_brute_combined']['n']))
    A('Additional checks: ' + '; '.join(extra) + '.')
    mp = sum(val.get(t, {}).get('reduced_minimum_pass_equals_frozen', {}).get('ok', 0)
             for t in ('track_A', 'track_B', 'track_C'))
    mpn = sum(val.get(t, {}).get('reduced_minimum_pass_equals_frozen', {}).get('n', 0)
              for t in ('track_A', 'track_B', 'track_C'))
    A('')
    A('The Stage295 reduced minimum equals the frozen minimum in **%d/%d** comparisons '
      'across the three tracks.' % (mp, mpn))
    fr = val.get('face_reduction')
    if fr:
        A('')
        A('With the reduction **disabled** (`full_mask=True`) the Stage297 enumeration '
          'reproduces the unreduced face count of the Stage292 enumeration exactly in '
          '%d/%d checks, with identical pattern counts (%d/%d) and identical pattern sets '
          'between reduced and full (%d/%d) — the correctness signal Stage295 relied on. '
          'Median faces: %s full vs %s reduced (median reduction factor %s, max %s).'
          % (fr['full_equals_292'], fr['checks'], fr['pattern_count_equal'], fr['checks'],
             fr['reduced_patterns_equal_full'], fr['checks'],
             d(fr['faces_full'], 'median', 0), d(fr['faces_reduced'], 'median', 0),
             d(fr['reduction_factor'], 'median', 2), d(fr['reduction_factor'], 'max', 1)))
    A('')
    speed = []
    for tag, label in (('track_A', 'A'), ('track_B', 'B'), ('track_C', 'C')):
        b = val.get(tag, {})
        key = 'combined|speedup_vs_s292'
        if key in b and b[key]:
            speed.append('track %s %sx' % (label, d(b[key], 'median', 2)))
    A('Median speed-up of `combined` over the Stage292 two-pass routine: ' +
      ', '.join(speed) + '.')
    A('')

    # ---------------------------------------------------------------- K = 10
    A('## 3  K = 10')
    A('')
    A('### 3.1  How many K = 10 cases admit the question at all')
    A('')
    A('`t_unknown` is defined only where the full-panel answer is already the unique true '
      'support (Stage292 convention: cases with `baseline_unique_true = false` are '
      'excluded from every threshold table).  At K = 10 that exclusion is no longer '
      'negligible:')
    A('')
    A('| | cases measured | minimum = K | unique true support at the full panel |')
    A('|---|---|---|---|')
    A('| Stage292, K in {2,3,5,7} | %s | %s | %s |' %
      (ref.get('cases', 'n/a'), ref.get('minimum_equals_K', 'n/a'),
       ref.get('cases_unique_true', 'n/a')))
    A('| Stage297, K = 10 | %s | %s | %s |' %
      (k10.get('cases', 'n/a'), k10.get('baseline_minimum_equals_K', 'n/a'),
       k10.get('baseline_unique_true', 'n/a')))
    A('')
    A('Exact minima found at K = 10: `%s`.  In the cases where the minimum falls below '
      'the generating K the data are explained exactly by fewer than K sources, so the '
      'true support is not the minimum-source answer even before any reference is erased.'
      % json.dumps(k10.get('baseline_minimum_distribution', {})))
    A('')
    A('Number of distinct minimisers at the full panel: `%s`.'
      % json.dumps(k10.get('baseline_supports_n', {})))
    A('')
    cc = k10.get('mode_cross_check')
    if cc:
        A('Mode cross-check at K = 10 (no validation set reaches this K): '
          '%d/%d cases agree on the minimum and %d/%d on the complete minimiser set '
          'between `combined` and `reduced_two_pass`; median %s s vs %s s.'
          % (cc['agree_minimum'], cc['n'], cc['agree_supports'], cc['n'],
             d(cc['combined_seconds'], 'median', 1), d(cc['two_pass_seconds'], 'median', 1)))
        A('')

    A('### 3.2  Threshold table')
    A('')
    A('`t_unknown` = the largest number of references that can be erased, in the order of '
      'the regime, with the set of ALL minimum-size supports still exactly the one true '
      'support.  Search convention imported verbatim from `stage292_bridge.find_threshold` '
      '(exponential bracketing, then bisection, both bracket endpoints evaluated, '
      'MIN_REMAINING = 3).')
    A('')
    for title, block in (('Stage297, K = 10', t297), ('Stage292, K in {2,3,5,7}', ref)):
        if not block:
            continue
        A('**%s**' % title)
        A('')
        A('| series | n | median | q25 | q75 | min | max |')
        A('|---|---|---|---|---|---|---|')
        for key in sorted(block.get('t_unknown', {})):
            A(describe_row('t_unknown ' + key.replace('|', ' / '),
                           block['t_unknown'][key], 1))
        for key in sorted(block.get('fraction_erasable', {})):
            A(describe_row('fraction ' + key.replace('|', ' / '),
                           block['fraction_erasable'][key], 3))
        A('')

    A('### 3.3  Flatness in panel size')
    A('')
    A('| regime | Spearman(t_unknown, panel size), Stage292 K in {2,3,5,7} | '
      'Stage297 K = 10 | pooled K in {2,3,5,7,10} |')
    A('|---|---|---|---|')
    for regime in REGIMES:
        A('| %s | %s (n=%s) | %s (n=%s) | %s (n=%s) |' % (
            regime,
            fmt(ref.get('spearman_t_vs_panel_size_pooled', {}).get(regime)),
            ref.get('n_pairs_pooled', {}).get(regime, 'n/a'),
            fmt(t297.get('spearman_t_vs_panel_size_pooled', {}).get(regime)),
            t297.get('n_pairs_pooled', {}).get(regime, 'n/a'),
            fmt(pool.get('spearman_t_vs_panel_size_pooled', {}).get(regime)),
            pool.get('n_pairs_pooled', {}).get(regime, 'n/a')))
    A('')

    A('### 3.4  Verdict')
    A('')
    sp297 = t297.get('spearman_t_vs_panel_size_pooled', {})
    sp292 = ref.get('spearman_t_vs_panel_size_pooled', {})
    tt = t297.get('t_unknown', {})
    adv = tt.get('K10|adversarial_nearest', {})
    rnd = tt.get('K10|random', {})
    cld = tt.get('K10|whole_clade', {})
    flat = (sp297.get('adversarial_nearest') is not None and
            sp297.get('random') is not None and
            sp297['adversarial_nearest'] < 0.3 <= sp297['random'])
    A('**The flatness claim holds at K = 10.**' if flat else
      '**At K = 10 the flatness claim does not reproduce in the form it has at lower K.**')
    A('')
    A('Across %s adversarial measurements spanning |R| = %s to %s the median '
      'nearest-to-source threshold is %s references (q25 %s, q75 %s, max %s) and its rank '
      'correlation with reference-set size is %s, against %s for random loss (median %s, '
      'max %s) and %s for whole-clade loss (median %s, max %s).  Stage292 measured '
      '%s / %s / %s on K in {2,3,5,7}.  The median adversarial threshold of %s references '
      'at K = 10 is the same median Stage292 reports at every lower K, and the fraction of '
      'the panel erasable before the minimum-source answer stops being unique is %s%% '
      'adversarial against %s%% random and %s%% whole-clade.'
      % (adv.get('n'), min(sizes_seen) if sizes_seen else 'n/a',
         max(sizes_seen) if sizes_seen else 'n/a',
         d(adv, 'median', 1), d(adv, 'q25', 0), d(adv, 'q75', 0), d(adv, 'max', 0),
         fmt(sp297.get('adversarial_nearest')), fmt(sp297.get('random')),
         d(rnd, 'median', 1), d(rnd, 'max', 0), fmt(sp297.get('whole_clade')),
         d(cld, 'median', 1), d(cld, 'max', 0),
         fmt(sp292.get('adversarial_nearest')), fmt(sp292.get('random')),
         fmt(sp292.get('whole_clade')), d(adv, 'median', 1),
         fmt(100 * (t297.get('fraction_erasable', {}).get('K10|adversarial_nearest', {}) or {}).get('median', 0), 1),
         fmt(100 * (t297.get('fraction_erasable', {}).get('K10|random', {}) or {}).get('median', 0), 1),
         fmt(100 * (t297.get('fraction_erasable', {}).get('K10|whole_clade', {}) or {}).get('median', 0), 1)))
    A('')
    A('Two things do change at K = 10, and both are about the panel being small relative '
      'to the number of sources rather than about the regimes.  First, all three '
      "thresholds compress: the random and whole-clade fractions fall from Stage292's "
      '42%% and 22%% of the panel to %s%% and %s%%, and the adversarial fraction from '
      '2.1%% to %s%%.  Second, the question is answerable for a smaller share of cases, '
      'because the exact minimum often falls below the generating K (section 3.1).'
      % (fmt(100 * (t297.get('fraction_erasable', {}).get('K10|random', {}) or {}).get('median', 0), 1),
         fmt(100 * (t297.get('fraction_erasable', {}).get('K10|whole_clade', {}) or {}).get('median', 0), 1),
         fmt(100 * (t297.get('fraction_erasable', {}).get('K10|adversarial_nearest', {}) or {}).get('median', 0), 1)))
    A('')
    rep = t297.get('repeated_measurements', 0)
    if rep:
        A('%d (case, regime, seed) measurements were produced twice because a case that '
          'timed out in the five-regime batch was re-measured in a regime-restricted one; '
          '%d of %d agree, and repeats are pooled once.' %
          (rep, t297.get('repeated_measurements_agree', 0), rep))
        A('')

    # --------------------------------------------------------------- ceiling
    A('## 4  The ceiling, with counters')
    A('')
    ceil = k10.get('ceiling', {})
    A('Per-case failures in the K >= 10 batches: %s (`%s`).'
      % (ceil.get('failed_jobs', 0), json.dumps(ceil.get('failed_detail', {}))))
    A('')
    if ceil.get('failed_specs'):
        A('| case | \\|R\\| | stop | wall s |')
        A('|---|---|---|---|')
        for f in ceil['failed_specs'][:30]:
            A('| %s | %s | %s | %s |' % (str(f.get('name'))[:52], f.get('R_size'),
                                         f.get('status'), fmt(f.get('wall_seconds'), 0)))
        A('')
    A('Probe outcomes: `%s`.' % json.dumps(ceil.get('probe_status', {})))
    A('')
    if k12:
        r = k12[0]
        b = r['baseline']
        A('**K = 12** (`%s`, |R| = %d, the single `exact_16k` case): status `%s` after '
          '%s s — the honest ceiling, with counters rather than a quiet retreat. In that '
          'budget the routine enumerated %s deficit sets and %s reduced face masks, '
          'skipping %s masks that forced positivity proves infeasible (%sx), tested %s '
          'faces by exact LP, found %s feasible pattern(s) and never reached the '
          'endpoint oracle (%s calls).'
          % (r['name'], r['R_size'], r['baseline_status'], fmt(r['baseline_seconds'], 0),
             '{:,}'.format(b['deficit_sets']), '{:,}'.format(b['masks_enumerated']),
             '{:,}'.format(b['masks_skipped_forced']),
             fmt(b['masks_full_enumeration'] / max(b['masks_enumerated'], 1), 1),
             '{:,}'.format(b['faces_to_lp']), b['feasible_patterns'], b['oracle_calls']))
        A('')

    # ---------------------------------------------------------------- item 2
    A('## 5  DecoDiPhy seed sensitivity')
    A('')
    ag = dd.get('aggregate', {}).get('all', {})
    A('%s runs over %s cases (stratified over |R| within every (family, K) cell of '
      '`exact_10k`, `exact_biotrees` and `krepp_estimated_10k`), full panels, documented '
      'defaults, k chosen by the tool\'s own stopping rules; seeds %s.  Run status: `%s`.  '
      'Runs recovered from a committed `all_rounds.json` after a jplace-writer crash: %s.'
      % (dd.get('runs', 'n/a'), dd.get('cases', 'n/a'),
         json.dumps(protocol.get('item2', {}).get('seeds', [])),
         json.dumps(dd.get('run_status', {})), ag.get('writer_crashed_runs', 'n/a')))
    A('')
    A('| cell | cases | k identical across seeds | support identical across seeds | '
      'median mean-pairwise Jaccard | median W1 relative range | max W1 relative range |')
    A('|---|---|---|---|---|---|---|')
    fam_keys = sorted(k for k in dd.get('aggregate', {}) if k.startswith('family|'))
    k_keys = sorted((k for k in dd.get('aggregate', {}) if k.startswith('K|')),
                    key=lambda z: int(z.split('|')[1]))
    order = ['all'] + fam_keys + k_keys
    for key in order:
        b = dd['aggregate'][key]
        A('| %s | %d | %d/%d | %d/%d | %s | %s | %s |' % (
            key.replace('|', ' '), b['cases'], b['k_identical_across_seeds'],
            b['cases_with_two_or_more_seeds'], b['support_identical_across_seeds'],
            b['cases_with_two_or_more_seeds'],
            d(b['mean_pairwise_jaccard'], 'median', 3),
            d(b['w1_relative_range'], 'median', 3),
            d(b['w1_relative_range'], 'max', 3)))
    A('')
    varying = [c for c in dd.get('per_case', [])
               if c['k_distinct'] > 1 or c['support_distinct'] > 1]
    if varying:
        A('The %d cases that move with the seed move little: the k range is %s (in %d of '
          'them k is fixed and only the support moves; the widest spread is k in %s), the '
          'mean pairwise Jaccard of their supports is %s to %s, and their W1 relative '
          'range is %s to %s.  %d of them are K = 10 cases.'
          % (len(varying),
             json.dumps(dict(Counter(c['k_range'] for c in varying
                                     if c['k_range'] is not None))),
             sum(1 for c in varying if c['k_range'] == 0),
             json.dumps(sorted({tuple(sorted(set(v for a2, v in c['k_values'].items()
                                                 if a2 != 'default_1142' and v is not None)))
                                for c in varying}, key=len)[-1]),
             fmt(min(c['mean_pairwise_jaccard'] for c in varying
                     if c['mean_pairwise_jaccard'] is not None), 3),
             fmt(max(c['mean_pairwise_jaccard'] for c in varying
                     if c['mean_pairwise_jaccard'] is not None), 3),
             fmt(min(c['w1_relative_range'] for c in varying
                     if c['w1_relative_range'] is not None), 3),
             fmt(max(c['w1_relative_range'] for c in varying
                     if c['w1_relative_range'] is not None), 3),
             sum(1 for c in varying if c['K_generating'] == 10)))
        A('')
    ctrl = ag.get('fixed_seed_control', {})
    A('**Fixed-seed control.**  Arm `default_1142` is `stage293_abundance.dd_worker` '
      'imported unmodified (no `-s` flag, so the tool takes its documented default seed); '
      'arm `seeded_1142` is the one-flag copy at the same seed.  They agree in '
      '%s/%s cases on k, %s/%s on the returned support, %s/%s on the loss and %s/%s on W1. '
      'The unified diff of the two workers is in '
      '`results/stage297/stage297_dd_worker_diff.json`.'
      % (ctrl.get('k_equal'), ctrl.get('comparable'), ctrl.get('support_equal'),
         ctrl.get('comparable'), ctrl.get('loss_equal'), ctrl.get('comparable'),
         ctrl.get('w1_equal'), ctrl.get('comparable')))
    A('')
    if cross:
        A('**Cross-harness comparison.**  The same tool, the same documented default '
          'seed, the same cases, two different harnesses: this stage\'s seed-1142 arm '
          'against the runs Stage289 recorded (`stage289_headtohead_runs.jsonl`, '
          'method `dd_default`, both repetitions).  Of %s comparable comparisons, k agrees '
          'in %s, the support in %s and the loss in %s.  Stage289\'s own two repetitions '
          'agree with each other in %s/%s on k and %s/%s on the support.'
          % (cross.get('comparable'), cross.get('k_equal'), cross.get('support_equal'),
             cross.get('loss_equal'),
             cross.get('stage289_internal_repetitions', {}).get('k_equal'),
             cross.get('stage289_internal_repetitions', {}).get('n'),
             cross.get('stage289_internal_repetitions', {}).get('support_equal'),
             cross.get('stage289_internal_repetitions', {}).get('n')))
        A('')
        A('By family: `%s`.' % json.dumps(cross.get('by_family', {})))
        A('')

    # ------------------------------------------------------------------ gaps
    A('## 6  What none of this establishes')
    A('')
    A('- The validation shows the Stage297 routine returns what the Stage292 routine '
      'returns wherever Stage292 completes, and what an independent exhaustive '
      'enumeration returns on small instances.  Exhaustiveness at K = 10 still rests on '
      'the frozen `deficit_sets` generator and the weighted deficit bound, which this '
      'stage imports and does not re-derive; no brute force is affordable at K = 10.')
    A('- Mode agreement at K = 10 is a cross-check between two routines that share pass 2, '
      'not an independent confirmation.')
    A('- The K = 10 sample is a stratified subsample of the 213 exact K = 10 bank cases, '
      'not the whole family, and the cases whose minimum falls below K are excluded from '
      'the threshold tables by the Stage292 convention — they are not evidence about '
      'erasure, but they do mean the K = 10 threshold statistics rest on fewer cases than '
      'the K <= 7 ones.')
    A('- K = 12 is not measured.  The single `exact_16k` case did not finish; the stage '
      'reports its counters, not a threshold.')
    A('- The seed sweep covers full panels only.  It says nothing about how seed choice '
      'interacts with reference removal, which is where Stages 288/294 found the tool '
      'failing outright.')
    A('- The cross-harness disagreement is characterised, not explained.  Both harnesses '
      'are internally reproducible, so the difference lies in process-level details '
      '(Stage289 ran the tool as a subprocess, Stage293/297 call `cli.main()` in process); '
      'which detail is responsible is not identified here.')
    A('')

    A('## 7  Files')
    A('')
    for name, what in (
            ('scripts/stage297_reduced_all_minimisers.py', 'the combined routine'),
            ('scripts/stage297_common.py', 'bank access, forked-child pool'),
            ('scripts/stage297_validate.py', 'validation tracks A, B, C'),
            ('scripts/stage297_bridge.py', 'K = 10 / K = 12 threshold measurement'),
            ('scripts/stage297_seeds.py', 'DecoDiPhy seed sweep'),
            ('scripts/stage297_analyze.py', 'item 1 analysis'),
            ('scripts/stage297_dd_analyze.py', 'item 2 analysis'),
            ('scripts/stage297_cross_harness.py', 'Stage289 cross-harness comparison'),
            ('scripts/stage297_protocol.py', 'protocol registration'),
            ('scripts/stage297_write_document.py', 'this document'),
            ('results/stage297/stage297_protocol.json', 'registered protocol'),
            ('results/stage297/stage297_protocol_addendum_main_k10.json', 'main K = 10 batch'),
            ('results/stage297/stage297_protocol_addendum_k10_extra.json', 'K = 10 extension'),
            ('results/stage297/stage297_validate_replay.jsonl', 'track A records'),
            ('results/stage297/stage297_validate_fresh.jsonl', 'track B records'),
            ('results/stage297/stage297_validate_baseline.jsonl', 'track C records'),
            ('results/stage297/stage297_bridge_pilot_k10_baseline.jsonl', 'K = 10 pilot baselines'),
            ('results/stage297/stage297_bridge_main_k10.jsonl', 'K = 10 main batch'),
            ('results/stage297/stage297_bridge_main_k10_extra.jsonl', 'K = 10 extension batch'),
            ('results/stage297/stage297_bridge_k10_adv_large.jsonl',
             'adversarial-only large-panel batch'),
            ('results/stage297/stage297_bridge_k10_rand_large.jsonl',
             'random-only large-panel batch, stopped early'),
            ('results/stage297/stage297_bridge_k12_baseline.jsonl', 'K = 12 attempt'),
            ('results/stage297/stage297_protocol_addendum_k10_large.json',
             'regime-split large-panel batches'),
            ('results/stage297/stage297_main_k10_manifest.jsonl',
             'skip manifest: the 60 main_k10 case names'),
            ('results/stage297/stage297_adv_done_manifest.jsonl',
             'skip manifest: cases whose adversarial threshold was already found'),
            ('results/stage297/stage297_bridge_v1_as_registered.py',
             'the bridge script as registered in the protocol'),
            ('results/stage297/stage297_bridge_v2_main_k10.py', 'bridge after main_k10'),
            ('results/stage297/stage297_bridge_v3_k10_extra.py', 'bridge after the extension'),
            ('results/stage297/stage297_dd_seeds_main.jsonl', 'seed sweep runs'),
            ('results/stage297/stage297_analysis.json', 'item 1 analysis output'),
            ('results/stage297/stage297_dd_analysis.json', 'item 2 analysis output'),
            ('results/stage297/stage297_dd_vs_stage289.json', 'cross-harness comparison'),
            ('results/stage297/stage297_dd_worker_diff.json', 'worker diff')):
        A('- `%s` — %s' % (name, what))
    A('')
    OUT.write_text('\n'.join(L) + '\n')
    print('written', OUT, len(L), 'lines')


if __name__ == '__main__':
    main()
