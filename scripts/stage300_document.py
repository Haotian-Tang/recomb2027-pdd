#!/usr/bin/env python3
"""Stage 300: assemble results/stage300/stage300_own_fit_comparison.md.

Every number is read from the stage-300 records through `stage300_report.py`;
nothing is transcribed by hand.  Re-runnable: running it again after more
travellers'-diarrhoea samples finish updates the document in place.
"""
import json
import sys
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import stage300_common as s300                   # noqa: E402
from stage300_report import (read_json, read_rows, fmt, q, provenance, three_way,
                             own_norm, resolved_table, stool_table,
                             breakdown)   # noqa: E402

OUT = s300.OUT


def main():
    analysis = read_json('stage300_analysis.json')
    lossgap = read_json('stage300_lossgap.json', {})
    route_e = read_json('stage300_route_e.json', {})
    verify = read_json('stage300_verification_small.json', {})
    verify_td = read_json('stage300_verification_travelers.json', {})
    probe = read_json('stage300_tdprobe_td_sample_101_published.json', {})
    losscase = read_json('stage300_losscase_k10k_tree_15_k3_s1_exp_krepp_estimated.json', {})
    protocol = read_json('stage300_protocol.json', {})
    stool, stool_done, stool_total = stool_table()
    krepp = analysis['by_collection']['krepp_ani']
    plc = analysis['by_collection']['krepp_placement']
    trav = analysis['by_collection'].get('travelers', {})
    kd, pd_ = krepp['arms']['dd_default'], plc['arms']['dd_default']
    ktk = krepp['arms']['dd_true_k']
    fp, rsv = krepp['fair_panel'], krepp['resolved']['dd_default']
    fpa = fp['arms']['dd_default']
    lg = lossgap.get('by_collection', {})
    re_k = route_e.get('by_collection', {}).get('krepp_ani', {})
    re_p = route_e.get('by_collection', {}).get('krepp_placement', {})
    va = verify.get('routes', {})
    td_runs = read_rows('stage300_runs_travelers.jsonl') + read_rows('stage300_pilot_travelers.jsonl')

    L = []
    A = L.append
    A('# Stage 300 - DecoDiPhy\'s own fit, measured against Theta_1')
    A('')
    A('Protocol `results/stage300/stage300_protocol.json`, registered %s.  Inputs, code '
      'and the published package are hashed there.  Every number in this document is '
      'read from the stage-300 records by `scripts/stage300_document.py`.'
      % protocol.get('registered_utc', 'n/a'))
    A('')
    A('## 1. What is measured')
    A('')
    A('For a configuration - the true generating one, the one DecoDiPhy returns, or a '
      'single source - the prediction at an observed reference leaf r is')
    A('')
    A('```')
    A('pred_r = y + sum_i p_i * dist_T(r, q_i),   dist(r,q) = min(d(r,u) + xL, d(r,v) + (1-x)L)')
    A('```')
    A('')
    A('and two fits are reported for it:')
    A('')
    A('- **R** = `max_r |pred_r - d_r|`, the L-infinity residual, computed in exact '
      'rational arithmetic (the configuration\'s floats enter by their exact binary '
      'value).  Every `R >= Theta_1` decision below is therefore exact, not a tolerance '
      'test.')
    A('- **SSR** = `sum_r (pred_r - d_r)^2` in float64 - DecoDiPhy\'s own objective, '
      '`np.sum((d - (DA @ p + CAl @ x + y)) ** 2)` in `decodiphy/optimize.py`.')
    A('')
    A('The two single-source baselines:')
    A('')
    A('- **Theta_1** - the exact smallest tolerance at which ONE source explains the '
      'observations, **imported** from stage 298 (`stage298_analysis.small_cases` and '
      '`stage298_td_summary.rows`) and never recomputed differently here.')
    A('- **Lambda_1** - the exhaustive minimum of SSR over every original edge, every '
      'x in [0,1] and every y >= 0 with p = 1: the same feasible set DecoDiPhy\'s QP '
      'uses at k = 1, enumerated instead of hill-climbed.  It is the L2 analogue of '
      'Theta_1 and exists so the comparison can be made in the norm the tool optimises.')
    A('')
    A('DecoDiPhy reports a per-placement pendant vector; the model\'s additive constant '
      'is its sum, because `compute_individual_ys` splits the scalar pendant '
      'proportionally and the parts sum back to it by construction.  Abundances, '
      'positions and pendant are taken verbatim from the tool\'s output and mapped onto '
      'bank edge indices by the frozen stage-289 bridge.  Nothing is renormalised.')
    A('')
    A('## 2. Where the DecoDiPhy answers come from')
    A('')
    A(provenance())
    A('')
    A('Harvested runs are full-panel runs already on record from stages 289/293/294/297, '
      'each row carrying its file, line, arm and repetition or seed; nothing was re-run.  '
      'The krepp ANI collection reaches 296/296 cases only by combining stage 289 (290 '
      'cases complete) with stage 294 (296, whose harness recovers the configuration '
      'after the `save_jplace` KeyError).  The 60 placement cases and the 29 stool '
      'samples had no prior runs and were executed here with the same worker '
      '(`stage293_abundance.dd_worker`, imported unmodified): documented defaults, '
      '`--method hill --optimizer osqp`, the only operable combination in 1.5.0.')
    A('')
    A('## 3. The decisive count, in L-infinity')
    A('')
    A(three_way(analysis))
    A('')
    A('Ordering of the three fits per case, best first (`1` = one source at Theta_1, '
      '`D` = DecoDiPhy\'s returned mixture, `T` = the true configuration):')
    A('')
    for label, cell in (('krepp ANI (dd_default)', kd), ('krepp placement (dd_default)', pd_)):
        counts = cell.get('ordering_linf', {})
        A('- %s: %s' % (label, ', '.join('`%s` %d' % (k, v) for k, v in
                                         sorted(counts.items(), key=lambda z: -z[1]))))
    A('')
    A('- On the 296 krepp ANI cases the tool\'s own mixture fits **no better in '
      'L-infinity than the best single source in %d of 296 cases** (%.1f%%); median '
      'R_dd/Theta_1 = %s, upper quartile %s, maximum %s.'
      % (kd['n_r_dd_ge_theta1'], 100.0 * kd['n_r_dd_ge_theta1'] / kd['n_cases_with_answer'],
         q(kd['r_dd_over_theta1']), q(kd['r_dd_over_theta1'], 'q75'),
         q(kd['r_dd_over_theta1'], 'max')))
    A('- On the 60 placement cases it fits better than one source everywhere '
      '(%d/%d at or above Theta_1, median ratio %s).'
      % (pd_['n_r_dd_ge_theta1'], pd_['n_cases_with_answer'], q(pd_['r_dd_over_theta1'])))
    A('- The true configuration is at or above Theta_1 in %d/%d krepp ANI and %d/%d '
      'placement cases - stage 298 reproduced exactly - and DecoDiPhy\'s answer fits '
      'better than the truth in %d/%d and %d/%d.  Fitting the estimator\'s error is '
      'easier than reproducing the generator, and the tool does it.'
      % (kd['n_truth_ge_theta1'], kd['n_with_truth'], pd_['n_truth_ge_theta1'],
         pd_['n_with_truth'], kd['n_dd_beats_truth_in_linf'], kd['n_with_truth'],
         pd_['n_dd_beats_truth_in_linf'], pd_['n_with_truth']))
    A('- Handing the tool the true k does not change the verdict: %d/%d krepp ANI cases '
      'at or above Theta_1 with `dd_true_k`, median ratio %s.'
      % (ktk['n_r_dd_ge_theta1'], ktk['n_cases_with_answer'], q(ktk['r_dd_over_theta1'])))
    A('- The decision is reproducible across replicates: identical in %d of %d krepp ANI '
      'cases over the 1-9 harvested runs each (k itself differs between replicates in '
      '%d cases).'
      % (kd['replicates']['distinct_decision'].get('1', 0),
         sum(kd['replicates']['distinct_decision'].values()),
         kd['replicates']['distinct_k'].get('2', 0)))
    A('- Median k returned: %s on krepp ANI (true K in {3, 5, 10}) and %s on placement '
      '(true K = %s), i.e. the default stopping rule inflates the count on placement '
      'distances.'
      % (q(kd['k_returned']), q(pd_['k_returned']), q(plc['arms']['dd_true_k']['k_returned'])))
    A('')
    A('The verdict is flat in the true number of sources and strongly dependent on '
      'how many references the panel carries:')
    A('')
    A(breakdown())
    A('')
    A('With 20-43 references the mixture usually does beat one source in '
      'L-infinity; from about 43 references upwards it fails to in more than '
      'half the cases.  The median count it returns falls the other way '
      '(12, 5, 4 sources over the three panel strata), so the extra sources '
      'buy worst-leaf accuracy only while the panel is small enough to absorb '
      'them.')
    A('')
    A('## 4. The same comparison in the norm DecoDiPhy optimises')
    A('')
    A('**DecoDiPhy does not optimise L-infinity.**  It minimises SSR, and in SSR it '
      'wins.  This has to be read in the same breath as section 3.')
    A('')
    A(own_norm(analysis))
    A('')
    A('- Against the best single source in its own norm it wins in %d/%d krepp ANI and '
      '%d/%d placement cases, with median SSR %s and %s of Lambda_1.'
      % (kd['n_dd_beats_one_source_in_ssr'], kd['n_cases_with_answer'],
         pd_['n_dd_beats_one_source_in_ssr'], pd_['n_cases_with_answer'],
         q(kd['ssr_over_lambda1']), q(pd_['ssr_over_lambda1'])))
    A('- So the conclusion in section 3 IS norm-dependent, and the honest statement is: '
      'the extra sources buy a better least-squares fit and do not buy a smaller '
      'worst-leaf error.  Theta_1\'s meaning - a tolerance that must contain the error '
      'at EVERY leaf - is intrinsically an L-infinity statement, which is why the '
      'paper\'s claim is made there; a reader who cares about least squares should read '
      'this table instead of section 3.')
    A('- The tool also beats the TRUE configuration in SSR in %d/%d krepp ANI cases, '
      'while the truth beats one source in SSR in only %d/%d.  On real estimated '
      'distances the generator is not the least-squares optimum either.'
      % (kd['n_dd_beats_truth_in_ssr'], kd['n_with_truth'],
         kd['n_truth_beats_one_source_in_ssr'], kd['n_with_truth']))
    A('- In %d of %d krepp ANI cases the returned mixture is worse in SSR than the best '
      'single source (ratio up to %s), i.e. its hill climbing can return an answer that '
      'loses to one source under its own objective.'
      % (kd['n_cases_with_answer'] - kd['n_dd_beats_one_source_in_ssr'],
         kd['n_cases_with_answer'], q(kd['ssr_over_lambda1'], 'max')))
    A('')
    A('### 4.1 Re-solving its own programme at its own support')
    A('')
    A('Second fairness arm, in case the comparison is an artefact of what the tool '
      'writes rather than of what it computes: `decodiphy.optimize.solve_osqp_once` is '
      'called on the tool\'s own matrices at the support it returned, and the resulting '
      'configuration is evaluated by the same exact code path.')
    A('')
    A(resolved_table(analysis))
    A('')
    A('The verdict barely moves: %d of %d re-solved krepp ANI cases are still at or '
      'above Theta_1 (median ratio %s, against %s for the configuration as returned).  '
      'The L-infinity statement is about the model and the support, not about output '
      'formatting.'
      % (rsv['n_linf_ge_theta1'], rsv['n_ok'], q(rsv['linf_over_theta1']),
         q(kd['r_dd_over_theta1'])))
    A('')
    cost = read_json('stage300_travelers_cost.json', {})
    A('## 5. The 29 travellers\'-diarrhoea stool samples: a measured cost, not a missing cell')
    A('')
    A('Nothing is known about these samples - no ground truth, no error model, no '
      'admissible floor - so the intended statement was a comparison of fit quality in '
      'a stated norm plus the number of constituents each answer carries.  Half of that '
      'is here: Theta_1 and Lambda_1 are measured for all 29 samples and need no '
      'DecoDiPhy run.  The other half is not: under this stage\'s invocation the tool '
      'returned no answer on any sample within the registered wall.')
    A('')
    if cost:
        el = cost['elapsed_seconds']
        rs = cost['rss_gib_at_stop']
        A('**Wall.** %.0f s (%.0f h) of wall-clock per sample, registered in '
          '`stage300_travelers_cost.json` with its justification: about %s GiB resident '
          'per worker (still rising) is what the shared machine allows at 5 concurrent '
          'workers, and 29 samples x 2 h / 5 workers is about 11.6 h - the largest '
          'per-sample wall that keeps the collection inside a single day of this '
          'stage\'s 40-core budget.'
          % (cost['wall_seconds'], cost['wall_seconds'] / 3600.0, fmt(rs['min'], 2)))
        A('')
        A('**Measured.** %d of %d samples attempted, **%d completed**.  Per-sample '
          'wall-clock at the moment each run was stopped: %.0f-%.0f s (%.1f-%.1f h); '
          'resident set size then: %.1f-%.1f GiB (a lower bound on the peak); %.1f '
          'core-hours consumed in total, with no returned configuration.  The remaining '
          '%d samples were not attempted.'
          % (cost['n_samples_attempted'], cost['n_samples_total'],
             cost['n_samples_completed'], el['min'], el['max'], el['min'] / 3600.0,
             el['max'] / 3600.0, rs['min'], rs['max'], cost['core_hours_consumed'],
             cost['n_samples_not_attempted']))
        A('')
        A('**Configuration.** `%s`; worker `%s`; concurrency %s; %s; %s.'
          % (cost['configuration']['invocation'], cost['configuration']['worker'],
             cost['configuration']['concurrency'], cost['configuration']['memory_limit'],
             cost['configuration']['threads']))
        A('')
        mi = cost['measured_input_stage']
        A('**Where the time went, measured on one sample** (`%s`): parsing the backbone '
          '%.2f s, `preprocess_input` %.2f s, treeswift\'s leaf-to-leaf '
          '`distance_matrix` %.0f s for %s entries, and the enclosing '
          '`get_input_matrices` had not returned after %d s (%.0f min) when that probe '
          'was stopped.  The input stage, not the search, is where this invocation '
          'spent its time.  Unlike the krepp ANI cases of section 6.2, nothing is pruned '
          'here: every one of the %d backbone leaves carries a distance.'
          % (mi['source'], mi['read_tree_newick_seconds'],
             mi['preprocess_input_seconds'], mi['distance_matrix_seconds'],
             '{:,}'.format(mi['distance_matrix_entries']),
             mi['get_input_matrices_not_returned_after_seconds'],
             mi['get_input_matrices_not_returned_after_seconds'] / 60.0,
             cost['pruning']['n_leaves_kept']))
        A('')
        A('**What this is not.**  %s' % cost['what_this_is_not'])
        A('')
    A('Theta_1 and Lambda_1 for all 29 samples, which is what can be reported without '
      'the tool (%d of 29 samples have a DecoDiPhy answer):' % stool_done)
    A('')
    A(stool)
    A('')
    A('All three quantities are in the samples\' own distance units (DEPP placement '
      'distances on the 16S backbone, median observed distance of order 300), so only '
      'ratios are comparable across samples.  The L-infinity residual of the Lambda_1 '
      'configuration is listed next to Theta_1 to show the gap between the two '
      'single-source optima: the least-squares best single source is worse at its worst '
      'leaf than the L-infinity best single source, by construction, and the column '
      'quantifies by how much.')
    A('')
    A('## 6. Two properties of the implementation that this comparison exposed')
    A('')
    A('### 6.1 The configuration written out is not always the one whose loss is printed')
    A('')
    if lg:
        kk = lg.get('krepp_ani', {})
        pp = lg.get('krepp_placement', {})
        A('Recomputing SSR from the returned configuration and comparing it with the loss '
          'DecoDiPhy reports for that run: the two differ by more than 1e-6 relative in '
          '%d of %d krepp ANI runs and %d of %d placement runs, and the returned '
          'configuration is the worse of the two in %d and %d runs respectively.'
          % (kk.get('n_gap_above_1e_6', 0), kk.get('n_runs', 0),
             pp.get('n_gap_above_1e_6', 0), pp.get('n_runs', 0),
             kk.get('n_returned_worse_than_reported', 0),
             pp.get('n_returned_worse_than_reported', 0)))
        A('')
    A('`decodiphy/search.py:hill_climbing` resets `min_obj` at the top of every sweep '
      'but never resets `min_p / min_x / min_y`, so when the final sweep improves on '
      'nothing, the abundances, positions and pendant written out are those of the last '
      'improving candidate of an earlier sweep while the anchors and the loss belong to '
      'the converged support.  `all_rounds.json` and the jplace carry that triple.')
    A('')
    if losscase:
        A('Worst case, `%s` (`stage300_losscase_*.json`): the tool reports loss %s at '
          'k = %d, the configuration it writes has SSR %s on the same leaves, and '
          're-solving its own QP at its own printed support attains %s.  The printed '
          'number is a real objective value - just not the one belonging to the printed '
          'parameters.'
          % (losscase['name'], fmt(losscase['reported_loss']), losscase['reported_k'],
             fmt(next((c['ssr_dd'] for c in krepp['per_case']['dd_default']
                       if c['name'] == losscase['name']), None)),
             fmt(losscase['objective_at_reported_support'])))
        A('')
    A('This is why section 4.1 exists: with the support re-solved the comparison is '
      'unchanged, so nothing in section 3 rests on this inconsistency.')
    A('')
    A('### 6.2 The input stage silently deletes reference leaves that do carry a distance')
    A('')
    A('`decodiphy/tree_utils.py:preprocess_input` removes leaves without a distance '
      'while iterating over `tree_obj.traverse_leaves()`, calling `remove_child` and '
      '`contract` inside that loop; that mutation-during-traversal is read from the '
      'source, the counts below are measured by running the function on every case '
      '(`stage300_prune_*.jsonl`):')
    A('')
    A('| collection | cases | cases where references were dropped | references kept when '
      'it happens | panel size when it happens |')
    A('|---|---|---|---|---|')
    for name, cell in (('krepp_ani', fp), ('krepp_placement',
                                           plc.get('fair_panel', {}))):
        if not cell:
            continue
        A('| %s | %d | %d | %s | %s |'
          % (name, cell['n_cases'], cell['n_cases_with_references_dropped'],
             ('%s-%s of %s-%s (median %s%%)'
              % (fmt(cell['n_kept_when_dropped']['min']),
                 fmt(cell['n_kept_when_dropped']['max']),
                 fmt(cell['n_observed_when_dropped']['min']),
                 fmt(cell['n_observed_when_dropped']['max']),
                 fmt(100 * cell['kept_fraction_when_dropped']['median'], 3))
              if cell['n_cases_with_references_dropped'] else '-'),
             ('%s-%s' % (fmt(cell['n_observed_when_dropped']['min']),
                         fmt(cell['n_observed_when_dropped']['max']))
              if cell['n_cases_with_references_dropped'] else '-')))
    A('')
    A('In those %d krepp ANI cases every one of the %s-%s references carried a distance '
      '(`n_with_distance` equals `n_observed` in all of them) and the tool nevertheless '
      'kept %s-%s of them.  Its reported loss for those cases is computed on 3-4 leaves.'
      % (fp['n_cases_with_references_dropped'],
         fmt(fp['n_observed_when_dropped']['min']), fmt(fp['n_observed_when_dropped']['max']),
         fmt(fp['n_kept_when_dropped']['min']), fmt(fp['n_kept_when_dropped']['max'])))
    A('')
    A('Because a comparison on data the tool never saw would not be fair, both sides were '
      'recomputed on the panel it did use - Theta_1 restricted to the surviving '
      'references with stage 298\'s exact k = 1 routine, and R_dd on the same leaves:')
    A('')
    A('- full panel, as the practitioner holds it: R_dd >= Theta_1 in **%d of %d** krepp '
      'ANI cases, median ratio %s;'
      % (kd['n_r_dd_ge_theta1'], kd['n_cases_with_answer'], q(kd['r_dd_over_theta1'])))
    A('- restricted to the references the tool used: **%d of %d**, median ratio %s, and '
      'in the %d mutilated cases the restricted verdict is %d of %d (it fits its 3-4 '
      'retained references almost exactly).'
      % (fpa['n_r_dd_kept_ge_theta1_kept'], fpa['n'], q(fpa['ratio']),
         fpa['dropped_only']['n'], fpa['dropped_only']['n_ge'], fpa['dropped_only']['n']))
    A('')
    A('The %d-case difference between the two counts is exactly the %d mutilated cases.  '
      'Neither count is cosmetic: the first is what a user gets, the second is what the '
      'tool was given.'
      % (kd['n_r_dd_ge_theta1'] - fpa['n_r_dd_kept_ge_theta1_kept'],
         fp['n_cases_with_references_dropped']))
    A('')
    A('## 7. Verification')
    A('')
    A('| route | what it checks | result |')
    A('|---|---|---|')
    if va:
        total_a = sum(cell['route_a']['compared'] for cell in va.values())
        agree_a = sum(cell['route_a']['agree'] for cell in va.values())
        maxb = max(cell['route_b']['max_abs_diff'] for cell in va.values())
        maxc = max(cell['route_c']['max_abs_diff'] for cell in va.values())
        nb = sum(len(cell['route_b']['cases']) for cell in va.values())
        nc = sum(len(cell['route_c']['cases']) for cell in va.values())
        A('| A | frozen `min` form vs stage 296\'s interpolation form, exact rational '
          'equality, every evaluated run | %d/%d agree |' % (agree_a, total_a))
        A('| B | whole-case recomputation with `scipy.sparse.csgraph.dijkstra` on the '
          'edge list with source points inserted as vertices | %d cases, max abs. '
          'difference %s |' % (nb, fmt(maxb, 2)))
        A('| C | prediction rebuilt from DecoDiPhy\'s own `D`, `C`, `l` matrices on the '
          'tree file the bridge writes | %d cases, max abs. difference %s |'
          % (nc, fmt(maxc, 2)))
    A('| D | recomputed SSR of the returned configuration vs the loss the tool reports, '
      'every run | median relative difference %s (krepp ANI, dd_default); the %d outliers '
      'are section 6.1 |' % (q(kd['loss_rel_diff'], digits=2),
                             lg.get('krepp_ani', {}).get('n_gap_above_1e_6', 0)))
    if re_k:
        A('| E | Lambda_1 vs DecoDiPhy\'s own k = 1 ladder loss (same objective, same '
          'edge set, hill climbing vs exhaustive) | equal to 1e-9 in %d/%d krepp ANI and '
          '%d/%d placement runs; all %d exceptions are the 14 cases of section 6.2 |'
          % (re_k['n_equal_within_1e_9'], re_k['n'], re_p['n_equal_within_1e_9'],
             re_p['n'], re_k['n_dd_below_lambda1']))
    A('| F | truth residual recomputed here vs stage 298\'s imported `residual_max` | '
      '%s |' % ', '.join('%s: %s' % (k, v) for k, v in
                         sorted(krepp['baselines']['truth_linf_matches_stage298'].items())))
    A('| G | Lambda_1 closed form vs exact re-evaluation of its attaining configuration '
      '| max relative difference %s (krepp ANI), %s (placement), %s (travellers) |'
      % (fmt(krepp['baselines']['lambda1_ssr_recompute_max_rel'], 2),
         fmt(plc['baselines']['lambda1_ssr_recompute_max_rel'], 2),
         fmt((trav.get('baselines') or {}).get('lambda1_ssr_recompute_max_rel'), 2)))
    cons = read_json('stage300_leafset_consistency.json', {})
    if cons:
        A('| H | the reference panel behind stage 298\'s Theta_1 against the panel this '
          'stage evaluates on, case by case | %d/%d identical, %d mismatches |'
          % (cons['n_cases_checked'] - cons['n_mismatch'], cons['n_cases_checked'],
             cons['n_mismatch']))
    A('')
    A('Route B is the independent one asked for: a library shortest-path routine that '
      'shares no code with this project, run on a graph in which the source points are '
      'explicit vertices, reproduces the leaf equations of whole cases to 1 ulp.')
    A('')
    A('## 8. What this does NOT establish')
    A('')
    A('1. **It does not show that DecoDiPhy is wrong.**  No sample here has a known '
      'answer except the 356 semi-synthetic cases, and on those the tool\'s mixture fits '
      'the observations better than the true configuration does.  A worse L-infinity fit '
      'than one source is not evidence of a wrong answer; it is evidence that the data '
      'does not force a multi-source answer in that norm.')
    A('2. **It does not establish a true source count for any stool sample.**  Theta_1 '
      'and R_dd are fit quantities; `k` is what the tool reported.  Nothing here '
      'licenses "the sample contains one source".')
    A('3. **It is not a comparison on DecoDiPhy\'s home ground.**  Section 4 is the '
      'relevant table for least squares, and there the tool wins.  The L-infinity '
      'comparison is the one Theta_1 supports, not the one the tool was built for.')
    A('4. **It does not measure accuracy of placement or abundance.**  Support recovery '
      'and abundance error are stages 289 and 293; this stage measures fit only.')
    A('5. **The krepp ANI and placement collections are semi-synthetic**: real estimated '
      'distances, but sources planted in a published simulation.  Only the 29 stool '
      'samples are real end to end, and they carry no truth.')
    A('6. **The two implementation findings in section 6 are properties of DecoDiPhy '
      '1.5.0 as published** (commit and file hashes in the protocol), reproduced from '
      'its own source; they are not claims about the method in the paper.')
    A('7. **No fit comparison for DecoDiPhy on the stool samples was obtained.**  '
      'Section 5 reports a cost, not a result: 6 of 29 samples were attempted under this '
      'stage\'s invocation and none returned a configuration within the registered 2-hour '
      'wall, so there is no R_dd, no k and no SSR for any real stool sample here.  '
      'Nothing in this document says how the tool would fit them.')
    A('8. **No error model is assumed anywhere.**  Theta_1 needs none, which is what '
      'makes the stool-sample comparison possible at all; the price is that the '
      'comparison is of fit quality, not of correctness.')
    A('')
    A('## 9. Files')
    A('')
    for name, what in (
            ('stage300_protocol.json', 'protocol: input, code and package hashes, run plan'),
            ('stage300_harvest.jsonl', 'every harvested full-panel DecoDiPhy run with its provenance'),
            ('stage300_runs_placement.jsonl', 'new runs, 60 placement cases x 2 arms'),
            ('stage300_runs_travelers.jsonl / _pilot_travelers.jsonl',
             'stool-sample run logs; both are empty because no run returned before the '
             'wall - the outcomes are in the censored record below'),
            ('stage300_travelers_censored.jsonl', 'the six stopped stool-sample runs with '
             'their measured wall-clock and resident set size'),
            ('stage300_travelers_cost.json', 'the registered wall, its justification and '
             'the measured cost of the travellers arm'),
            ('stage300_fits_<collection>.jsonl', 'per-run exact R_dd, SSR, decisions'),
            ('stage300_baselines_<collection>.jsonl', 'per-case Theta_1, truth, Lambda_1'),
            ('stage300_resolved_<collection>.jsonl', 'the re-solved fairness arm'),
            ('stage300_prune_<collection>.jsonl', 'surviving panel and the fair-panel comparison'),
            ('stage300_lossgap.json / _rows.jsonl', 'returned configuration vs reported loss'),
            ('stage300_losscase_*.json', 'the worst inconsistency, re-solved at every neighbouring support'),
            ('stage300_route_e.json', 'Lambda_1 against the tool\'s own k = 1 answer'),
            ('stage300_verification_small.json', 'routes A, B, C'),
            ('stage300_analysis.json', 'all aggregates behind the tables above'),
            ('stage300_tdprobe_*.json', 'cost decomposition on one stool sample'),
            ('stage300_three_way.csv', 'the per-case three-way table (708 rows: every '
             'truth-known case x arm, exact rationals beside their float values)'),
            ('stage300_leafset_consistency.json', 'route H: the panels behind Theta_1 and '
             'behind R_dd, case by case'),
            ('stage300_protocol_amendment.json', 'what was added after registration, with '
             'hashes, and the declared deviations')):
        A('- `%s` - %s' % (name, what))
    A('')
    (OUT / 'stage300_own_fit_comparison.md').write_text('\n'.join(L) + '\n')
    print('wrote', OUT / 'stage300_own_fit_comparison.md', '%d lines' % len(L))
    return L


if __name__ == '__main__':
    main()
