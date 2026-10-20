#!/usr/bin/env python3
"""Stage 298: write stage298_real_data_axis.md and the compact figure data."""
import json
import os
import sys
from fractions import Fraction as F

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
import stage288_common as common            # noqa: E402
import stage298_common as s298              # noqa: E402

OUT = s298.OUT
DOC = OUT / 'stage298_real_data_axis.md'
FIGDATA = OUT / 'stage298_figure_data.json'


def fm(x, n=4):
    if x is None:
        return 'n/a'
    return ('%.' + str(n) + 'g') % x


def qline(d, n=4):
    if not d:
        return 'n/a'
    return ' / '.join(fm(d[k], n) for k in ('0.0', '0.25', '0.5', '0.75', '1.0'))


def main():
    A = json.loads((OUT / 'stage298_analysis.json').read_text())
    DEG = json.loads((OUT / 'stage298_witness_degeneracy.json').read_text())
    P = json.loads((OUT / 'stage298_protocol.json').read_text())
    M = json.loads((OUT / 'stage298_placement_manifest.json').read_text())
    K2 = json.loads((OUT / 'stage298_td_k2_cost.json').read_text())
    td = A['travelers']
    tdr = A['travelers_rows']
    ka = A['by_family']['krepp_estimated_10k']
    kp = A['by_family']['krepp_placement_10k']
    L = []
    w = L.append
    w('# Stage 298 - the real-data tolerance axis')
    w('')
    w('`Theta_1` is the exact smallest tolerance at which a SINGLE source explains the')
    w('observations: the minimum over every original edge, and over the closed mass set')
    w('`{alpha, beta, y >= 0, alpha + beta = 1}`, of the largest absolute residual at an')
    w('observed reference leaf. It is a rational number computed from the distance data')
    w('alone - no error model, no ground truth, no admissible floor. Stage 296 proved')
    w('`m(eps) = min{k : Theta_k < eps}` for `eps > 0`, so the vector')
    w('`(Theta_1, Theta_2, ...)` fixes the whole tolerance axis.')
    w('')
    w('**The question this answers for a practitioner: is my sample distinguishable from')
    w('a single source at all, given how accurate my distances are?** If the L-infinity')
    w('error of the pipeline that produced the distances exceeds `Theta_1`, then no')
    w('method - exact or heuristic - can justify reporting a mixture: the data admit a')
    w('one-source explanation within their own noise.')
    w('')
    w('All three collections below are REAL measured distances. Every `Theta_k` is')
    w('decided by frozen stage-296 code; every attaining witness is re-verified with the')
    w('frozen `min(d(r,u) + xL, d(r,v) + (1-x)L)` form on the original tree.')
    w('')
    w('## 1. The 29 travellers\'-diarrhoea stool samples (real 16S, no ground truth)')
    w('')
    w('Backbone: |E| = %d original edges, |R| = %d observed reference leaves, one'
      % (tdr[0]['n_edges'], tdr[0]['n_leaves']))
    w('tree shared by all 29 samples; distances are DEPP placement estimates published')
    w('with the DecoDiPhy benchmark. `Theta_1` was computed exhaustively over every edge')
    w('TWICE - by the frozen per-support LP and by the stage-298 closed form - giving')
    w('%s independent agreements and %d disagreements.'
      % ('{:,}'.format(td['edges_decided_twice']), td['closed_form_mismatches']))
    w('')
    w('| sample | group | Theta_1 (obs. units) | Theta_1 / median obs. | attaining edge | terminal | x (closed form) |')
    w('|---|---|---|---|---|---|---|')
    for r in sorted(tdr, key=lambda z: z['theta1_rel']):
        w('| %s | %s | %s | %s | %d (%d,%d) | %s | %s |'
          % (r['sample'], r['group'], fm(r['theta1_float'], 5), fm(r['theta1_rel'], 3),
             r['argmin_edge'], r['argmin_endpoints'][0], r['argmin_endpoints'][1],
             'yes' if r['argmin_is_terminal'] else 'no', fm(r['argmin_x'], 3)))
    w('')
    w('Distribution over the 29 samples (min / Q1 / median / Q3 / max):')
    w('')
    w('- `Theta_1` in observation units: %s' % qline(td['theta1_abs'], 5))
    w('- `Theta_1` / median observed distance: %s' % qline(td['theta1_rel'], 3))
    w('- median observed distance: %s' % qline(td['median_observed'], 5))
    w('- attaining edge terminal (one endpoint a leaf) in %d of %d samples; %d distinct'
      % (td['n_argmin_terminal'], td['n'], td['distinct_argmin_edges']))
    w('  attaining edges across the 29 samples.')
    w('- attaining witness independently verified in %d of %d samples.'
      % (td['witnesses_agree'], td['n']))
    w('- the attaining position need not be unique: `theta(S)` is a minimum over the')
    w('  CLOSED mass set, so it can be attained on a face.  The `x` column is the closed')
    w('  form\'s attaining position; the frozen LP may return a different point of the')
    w('  same optimal face, and in %d of %d samples the two differ while giving the'
      % (DEG['travelers']['x_differs'], td['n']))
    w('  identical exact `Theta_1`.  No travellers\' witness places the source at an edge')
    w('  endpoint; %d of %d small-case k=1 witnesses do, and no witness anywhere carries'
      % (DEG['k1']['boundary_x'], DEG['k1']['n']))
    w('  zero mass.  This is expected and is why `m(eps)` uses a STRICT inequality: stage')
    w('  296 proved the strict system (`alpha, beta, y > 0`, hence `0 < x < 1`) is')
    w('  feasible for every `eps > theta(S)`, so a degenerate attainer at the threshold')
    w('  never enters an `m(eps)` verdict.')
    w('')
    w('### Theta_2 on this backbone')
    w('')
    w('Not affordable, and the cost is measured rather than estimated: %d random'
      % K2['n_probe_supports'])
    w('two-edge supports of sample `%s` under the frozen LP cost' % K2['sample'])
    w('%.3f s per support (min %.3f, max %.3f), and the exhaustive sweep is'
      % (K2['seconds_per_support']['median'], K2['seconds_per_support']['min'],
         K2['seconds_per_support']['max']))
    w('C(%d, 2) = %s supports, i.e. %s core-hours for ONE sample.'
      % (K2['n_edges'], '{:,}'.format(K2['exhaustive_support_count']),
         '{:,.0f}'.format(K2['measured_core_hours_for_exhaustive_theta2'])))
    w('Stage 296 showed the frozen structural localisation is unsound for `eps > 0`, so')
    w('the sweep cannot be pruned structurally, and incumbent pruning does not help')
    w('because the dominant cost is building each support\'s columns with the frozen BFS,')
    w('which happens before the LP decides anything. `Theta_2 <= Theta_1` always, so the')
    w('one-source statement below does not depend on it.')
    w('')
    w('## 2. All 296 real krepp ANI cases')
    w('')
    w('`krepp_estimated_10k`: real ANI distances estimated from reads, N = %d-%d'
      % (min(r['N'] for r in A['small_cases'] if r['family'] == 'krepp_estimated_10k'),
         max(r['N'] for r in A['small_cases'] if r['family'] == 'krepp_estimated_10k')))
    w('reference leaves, K_true in {3, 5, 10}. Each case has an exact partner in the')
    w('bank - the model-implied distance vector for the same tree and the same')
    w('generating taxa - so the real estimation error is measurable in the same units as')
    w('`eps`: the L-infinity residual `max_r |d_est(r) - d_exact(r)|` on the shared leaf')
    w('set, which is stage 296\'s definition, applied here to all 296 cases (stage 296')
    w('measured 24 of them; the extension reproduces all 24 exactly).')
    w('')
    w('| quantity (min / Q1 / median / Q3 / max) | krepp ANI (n=%d) | krepp placement (n=%d) |'
      % (ka['n'], kp['n']))
    w('|---|---|---|')
    for lab, key, nd in (('Theta_1 / median observed', 'theta1_rel_own', 3),
                         ('Theta_2 / Theta_1', 'theta2_over_theta1', 3),
                         ('Theta_3 / Theta_1', 'theta3_over_theta1', 3),
                         ('residual (worst leaf) / median observed', 'residual_max_rel_own', 3),
                         ('residual (median leaf) / median observed', 'residual_median_rel_own', 3),
                         ('worst-leaf residual / Theta_1', 'residual_max_over_theta1', 3),
                         ('median-leaf residual / Theta_1', 'residual_median_over_theta1', 3)):
        w('| %s | %s | %s |' % (lab, qline(ka[key], nd), qline(kp[key], nd)))
    w('')
    w('**The headline count.** The real worst-leaf residual already exceeds `Theta_1` in')
    w('%d of %d krepp ANI cases (%.1f%%) and in %d of %d placement cases (%.1f%%): at a'
      % (ka['n_one_source_within_max_residual'], ka['n_with_residual'],
         100.0 * ka['n_one_source_within_max_residual'] / ka['n_with_residual'],
         kp['n_one_source_within_max_residual'], kp['n_with_residual'],
         100.0 * kp['n_one_source_within_max_residual'] / kp['n_with_residual']))
    w('tolerance large enough to contain the measured error at EVERY leaf, the model\'s')
    w('own minimum is one source. Using the median-leaf residual instead - a tolerance')
    w('that contains only half the leaves\' error - the count is %d of %d (krepp ANI) and'
      % (ka['n_one_source_within_median_residual'], ka['n_with_residual']))
    w('%d of %d (placement).' % (kp['n_one_source_within_median_residual'], kp['n_with_residual']))
    w('')
    w('`m(eps)` at `eps` = worst-leaf residual, krepp ANI: %s'
      % json.dumps(ka['m_at_max_residual'], sort_keys=True))
    w('')
    w('`m(eps)` at `eps` = median-leaf residual, krepp ANI: %s'
      % json.dumps(ka['m_at_median_residual'], sort_keys=True))
    w('')
    w('(`k:exact` means `Theta_k < eps <= Theta_{k-1}` with both computed; `k:lower_bound`')
    w('means the sweep was censored at that k, so `m(eps) >= k`.)')
    w('')
    w('**How much accuracy would move a sample off the collapsed side.** For the cases')
    w('whose worst-leaf residual exceeds `Theta_1`, the ratio residual / `Theta_1` is the')
    w('factor by which the pipeline\'s L-infinity distance error would have to shrink')
    w('before a mixture becomes justifiable at all: %s (min / Q1 / median / Q3 / max)'
      % qline(ka['accuracy_factor_to_leave_collapsed'], 3))
    w('for krepp ANI and %s for placement.' % qline(kp['accuracy_factor_to_leave_collapsed'], 3))
    w('')
    w('By K_true, krepp ANI:')
    w('')
    w('| K_true | n | Theta_1 / median obs. (median) | one source within worst-leaf residual |')
    w('|---|---|---|---|')
    for k, d in sorted(ka['by_K'].items(), key=lambda z: int(z[0])):
        w('| %s | %d | %s | %d |' % (k, d['n'], fm(d['theta1_rel_own']['0.5'], 3),
                                     d['n_one_source_within_max_residual']))
    w('')
    w('## 3. The second real distance estimator: krepp placement distances')
    w('')
    w('Outcome: **converted, %d of %d instances, 0 refused.**'
      % (M['n_converted'], M['n_instances_with_placement_file']))
    w('')
    w('The placement-distance tables sit in the Demixing_Data instance directories as')
    w('`distances_from_placement.txt` for %d of the 656 instances (trees 15 and 17 only).'
      % M['n_instances_with_placement_file'])
    w('Conversion used the frozen stage-287 pipeline unchanged - `build_exact_case` for')
    w('the normal form, `attach_published_observed` for the observation vector, and')
    w('`verify_case` for conformance - plus strict positivity of every observed value.')
    w('All %d converted cases carry byte-identical edges, leaf labels and generating'
      % M['n_converted'])
    w('sources to the bank\'s exact case for the same instance, so their observations are')
    w('comparable leaf by leaf with the exact partner and with the ANI estimates.')
    w('')
    w('One documented difference from the stage-287 krepp branch: that branch required')
    w('`set(published labels) == set(reference leaf labels)`. The placement files report')
    w('a distance for EVERY leaf of the tree, including the K query taxa the model turns')
    w('into interior sources and removes from the reference tree, so the published label')
    w('set is a strict superset. The rule applied here is containment - every reference')
    w('leaf must carry a positive placement distance - with the query-taxon entries')
    w('dropped, which the model requires because observations live at reference leaves')
    w('only. The dropped entries are exactly the K query taxa in %d of %d cases.'
      % (A['verification']['placement_conversion']['dropped_entries_are_exactly_the_queries'],
         M['n_converted']))
    w('')
    w('`truemash` remains unusable for the reason the gap audit records (every ANI tree')
    w('it produces carries negative branch lengths, so no PDD-conformant reference tree')
    w('exists); no reconstructed or repaired tree was substituted and no simulated noise')
    w('was used as a stand-in.')
    w('')
    w('## 4. Where the travellers\' samples sit')
    w('')
    w('The travellers\' samples have no exact partner, so their estimator error cannot be')
    w('measured from the data. Transferring the measured krepp relative error as a')
    w('REFERENCE SCALE (a different estimator on a different backbone - illustrative,')
    w('not a measurement of DEPP accuracy): the median krepp worst-leaf residual is')
    w('%s x the median observed distance and the median-leaf residual %s x. Against'
      % (fm(td['reference_relative_error']['krepp_residual_max_rel_median'], 3),
         fm(td['reference_relative_error']['krepp_residual_median_rel_median'], 3)))
    w('those two lines, %d of %d and %d of %d travellers\' samples respectively have'
      % (td['n_below_krepp_max_rel'], td['n'], td['n_below_krepp_median_rel'], td['n']))
    w('`Theta_1` below the line - i.e. sit on the collapsed side, where the data admit a')
    w('one-source explanation. Moving a sample off that line needs the L-infinity')
    w('distance error reduced by a factor of %s (min / Q1 / median / Q3 / max) relative'
      % qline(td['accuracy_factor_vs_krepp_max'], 3))
    w('to the krepp worst-leaf scale.')
    w('')
    w('## 5. Verification and censoring')
    w('')
    v = A['verification']
    w('- `Theta_k` cells completed: %d; failed: %d.' % (v['cells_complete'], v['cells_failed']))
    w('- k = 1 was run support by support with no cutoff on all %d small real cases:'
      % (ka['n'] + kp['n']))
    w('  %s edges compared between the frozen LP and the closed form, %d disagreements.'
      % ('{:,}'.format(v['k1_edges_compared']), v['k1_closed_form_mismatches']))
    w('- travellers: %s edges decided twice, %d disagreements.'
      % ('{:,}'.format(v['travelers_edges_decided_twice']),
         v['travelers_closed_form_mismatches']))
    w('- attaining witnesses re-verified with the frozen min(...) form: %d of %d agree.'
      % (v['witnesses_agree'], v['witnesses_checked']))
    w('- residual extension validated against stage 296: %d of %d rows agree exactly.'
      % (v['residual_extension']['agree'], v['residual_extension']['n_compared']))
    w('- censored cells: %d (%s), support cap %s; measured seconds per support %s.'
      % (A['censoring']['n_cells'], json.dumps(A['censoring']['by_k'], sort_keys=True),
         json.dumps(A['censoring']['support_cap'], sort_keys=True),
         qline(A['censoring']['measured_seconds_per_support'], 3)))
    w('')
    w('## 6. What this does NOT establish')
    w('')
    for s in [
        'It does not give a source count for any real sample. The travellers\' samples '
        'have no ground truth at all, and for the krepp and placement cases the '
        'generating taxa are known but the estimated distances are not exactly '
        'explainable by any configuration, so the minimum at tolerance zero is not the '
        'true count either.',
        'It does not measure the accuracy of DEPP placement distances. The travellers\' '
        'axis is reported in the sample\'s own units; the krepp relative residual used '
        'in section 4 is a reference scale from a different estimator on a different '
        'backbone, not a measurement of this pipeline\'s error.',
        'It does not show that the samples ARE single-source. Theta_1 below the assay '
        'error means the data cannot JUSTIFY a mixture, not that no mixture is present; '
        'the one-source explanation is admissible, not identified.',
        'It says nothing about reference incompleteness on real data. Theta_1 is '
        'computed at the full reference panel only; the removal regimes of stages 288-294 '
        'were not run here.',
        'Theta_2 and Theta_3 are not available everywhere: %d cells are censored by the '
        'registered support cap, and no Theta_2 exists for any travellers\' sample. '
        'Where a cell is censored, m(eps) is reported as a lower bound, never '
        'interpolated.' % A['censoring']['n_cells'],
        'It does not compare against DecoDiPhy or any heuristic. The axis is a property '
        'of the model and the data; no method was run on these samples here.',
        'The krepp residual is measured against the model-implied distances for the '
        'known generating taxa on the positive-length reference tree. It is the error '
        'of the estimator relative to that configuration, and it inherits the bank\'s '
        'pairing choice (the authors\' ANI branch-length trees carry negative lengths '
        'and are not PDD-conformant).',
    ]:
        w('- ' + s)
    w('')
    w('## 7. Figure')
    w('')
    w('![Real-data tolerance axis](stage298_real_data_axis.png)')
    w('')
    w('**(a)** Each point is one real case: its exact `Theta_1` against the measured')
    w('worst-leaf residual of the estimator that produced its distances, both divided by')
    w('that case\'s median observed distance. Points above the dashed identity sit on the')
    w('collapsed side - a tolerance containing the real error at every leaf already')
    w('admits a one-source explanation. **(b)** The 29 travellers\' samples have no')
    w('measurable residual, so their `Theta_1` is shown against the krepp relative-error')
    w('scale (dotted: median-leaf, dashed: worst-leaf) as a reference from a different')
    w('estimator, not a measurement of DEPP accuracy.')
    w('')
    w('## 8. Files')
    w('')
    for f, d in [
        ('results/stage298/stage298_protocol.json', 'protocol, registered before any timed run'),
        ('results/stage298/stage298_placement_manifest.json', 'placement conversion, per-instance outcome and hashes'),
        ('results/stage298/stage298_sweep_runs.jsonl', 'every Theta_k run on the 356 small real cases'),
        ('results/stage298/stage298_sweep_censored.json', 'every censored (case, k) cell with its support count'),
        ('results/stage298/stage298_td_runs.jsonl', 'every travellers\' edge-chunk run'),
        ('results/stage298/stage298_td_summary.json', 'per-sample travellers\' result and witness verification'),
        ('results/stage298/stage298_td_k2_cost.json', 'measured cost of Theta_2 on the travellers\' backbone'),
        ('results/stage298/stage298_residuals.json', 'real estimation residuals, all 356 cases'),
        ('results/stage298/stage298_analysis.json', 'assembled axis'),
        ('results/stage298/stage298_witness_degeneracy.json',
         'degeneracy census of the attaining witnesses'),
        ('results/stage298/stage298_protocol_amendment.json',
         'post-registration code hashes'),
        ('results/stage298/stage298_real_data_axis.png', 'the figure'),
        ('data/stage298_placement/cases/', '60 converted placement cases'),
    ]:
        w('- `%s` - %s' % (f, d))
    DOC.write_text('\n'.join(L) + '\n')
    print('wrote', DOC)

    fig = dict(
        krepp=[dict(theta1_rel=r['theta_rel_own']['1'],
                    resid_max_rel=r['residual_max_rel_own'],
                    resid_med_rel=r['residual_median_rel_own'],
                    K=r['K_generating'], N=r['N'])
               for r in A['small_cases'] if r['family'] == 'krepp_estimated_10k'
               and 'residual_max' in r],
        placement=[dict(theta1_rel=r['theta_rel_own']['1'],
                        resid_max_rel=r['residual_max_rel_own'],
                        resid_med_rel=r['residual_median_rel_own'],
                        K=r['K_generating'], N=r['N'])
                   for r in A['small_cases'] if r['family'] == 'krepp_placement_10k'
                   and 'residual_max' in r],
        travelers=[dict(theta1_rel=r['theta1_rel'], sample=r['sample'],
                        terminal=r['argmin_is_terminal']) for r in tdr],
        reference=td['reference_relative_error'],
        counts=dict(krepp_collapsed=ka['n_one_source_within_max_residual'],
                    krepp_n=ka['n_with_residual'],
                    placement_collapsed=kp['n_one_source_within_max_residual'],
                    placement_n=kp['n_with_residual'],
                    td_below_max=td['n_below_krepp_max_rel'], td_n=td['n']))
    FIGDATA.write_text(json.dumps(fig))
    print('wrote', FIGDATA)


if __name__ == '__main__':
    main()

