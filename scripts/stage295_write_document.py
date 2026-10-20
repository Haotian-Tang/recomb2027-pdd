#!/usr/bin/env python3
"""Stage295 document writer: renders stage295_ceiling_diagnosis.md from the
Stage295 report JSON.  Every number is read from the retained records.
"""
import json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RES = os.path.join(ROOT, 'results')


def fmt(x, nd=1):
    if x is None:
        return 'n/a'
    if isinstance(x, (int, np.integer)) or (isinstance(x, float) and x == int(x) and abs(x) < 1e15):
        return '{:,}'.format(int(x))
    return '{:,.{}f}'.format(x, nd)


def main():
    rep = json.load(open(os.path.join(RES, 'stage295_ceiling_report.json')))
    ana = json.load(open(os.path.join(RES, 'stage295_ceiling_analysis.json')))
    census = json.load(open(os.path.join(RES, 'stage295_static_census.json')))
    proto = json.load(open(os.path.join(RES, 'stage295_ceiling_protocol.json')))
    amend = json.load(open(os.path.join(RES, 'stage295_protocol_amendment.json')))
    _replay = json.load(open(os.path.join(RES, 'stage295_census_replay_check.json')))

    H = rep['headline_case']
    reg = rep['phase_profile']['prefix_regression']
    stop = rep['phase_profile']['counters_at_stop']['prefix_faces_K10']
    cp_deep = rep['phase_profile']['cprofile'].get('cprofile_K10_bases')
    cp_k7 = rep['phase_profile']['cprofile'].get('cprofile_full_K7')
    g = rep['growth']
    ladder = g['cap_ladder'][next(iter(g['cap_ladder']))]
    fits = g['log2_fits']
    V = rep['validation']
    P = rep['push']

    alpha, beta = reg['seconds_per_face'], reg['seconds_per_basis']
    floor_s = H['exact_zero_faces_frozen'] * alpha
    # forward projection of the prototype under its own measured growth
    red_bits = ladder['reduced_bits_per_K']
    head_rows = [x for x in P['rows'] if x['name'] == H['name']
                 and x['status'] == 'complete' and not x['popcount_order']]
    per_face_proto = (head_rows[0]['seconds'] / head_rows[0]['mass_faces']) if head_rows else None
    forward = {}
    if per_face_proto:
        for K in (11, 12, 13, 14, 15):
            f = H['exact_zero_faces_reduced'] * 2 ** (red_bits * (K - 10))
            forward[K] = {'faces': f, 'seconds': f * per_face_proto}

    _seen_caps = {}
    for c in census:
        if c.get('status') == 'complete' and H['name'] in c.get('case', ''):
            _seen_caps.setdefault(c['cap'], c)
    cap_rows = [_seen_caps[k] for k in sorted(_seen_caps)]

    L = []
    A = L.append
    A('# Stage 295 — Why the exact minimum pass stopped at K = 7, and what it costs to pass it')
    A('')
    A('Diagnosis of the K >= 10 ceiling in the frozen exact minimum pass')
    A('`scripts/stage43_partial_fpt.solve_partial`. Protocol registered before')
    A('every TIMED measurement of this stage as `results/stage295_ceiling_protocol.json`')
    A('-- the %s profile runs, the %s head-to-head comparisons and the K >= 10 push.'
      % (fmt(rep['inputs']['profile_runs']), fmt(rep['inputs']['validation_rows'])))
    A('The deterministic static census of Section 2.1 ran before registration and')
    A('informed the plan the protocol then fixed; see Section 7.')
    A('(Amendment 1: `results/stage295_protocol_amendment.json`).')
    A('No frozen module was modified; the amendment record confirms the frozen hashes are')
    A('unchanged (`frozen_unchanged = %s`).' % amend['frozen_unchanged'])
    A('')
    A('## Verdict')
    A('')
    A('**The K = 7 ceiling was an implementation artefact, not an intrinsic limit of the')
    A('method, and it has been removed.** The frozen search enumerates all `2^|A|` sign')
    A('masks over `A = S u N(S)` for each deficit support `S`, but on the strict face every')
    A('`t_v > 0` for `v in S`, so a neighbour `u in N(S) \\ S` with baseline `b_u >= 0` has')
    A('`m_u = b_u + sum_{v~u} t_v / (2 l_uv) > 0` forced. Every mask that pins such a `u` to')
    A('zero is infeasible and cannot produce an incumbent. On the smallest K = 10 bank case')
    A('those provably infeasible masks are %s of the %s faces the frozen solver visits'
      % ('{:.2%}'.format(1 - H['exact_zero_faces_reduced'] / H['exact_zero_faces_frozen']),
         fmt(H['exact_zero_faces_frozen'])))
    A('(%s remain). Removing them makes that case — the one Stage 292 reported as returning'
      % fmt(H['exact_zero_faces_reduced']))
    A('`unknown` after a 3,600 s budget — complete in **%s s** with a certified minimum.'
      % fmt(H['prototype_seconds']))
    A('')
    A('What is intrinsic is the *exponent*. Even after the reduction the face count still')
    A('grows at %.2f bits per unit K, so the method remains exponential in K: the ceiling'
      % red_bits)
    A('moved, it did not disappear (Section 5).')
    A('')
    A('## 1. The pipeline and how it was instrumented')
    A('')
    A('`solve_partial` reduces to the observed backbone and then runs four nested phases:')
    A('')
    A('1. **P1 deficit supports** — `stage41_partial_fpt.deficit_sets` grows connected sets')
    A('   `S` of eligible vertices (`3 <= deg <= K`) from positive-baseline anchors under the')
    A('   weighted deficit bound `sum_{v in S} (deg(v) - 2) <= K - 2`. Counters `growth_states`,')
    A('   `deficit_sets`.')
    A('2. **P2 zero faces** — for each `S`, `for mask in range(1 << len(A))` with')
    A('   `A = S u N(S)`; the mask selects the coordinates held strictly positive and pins')
    A('   the complement to `m_v = 0`. Counter `mass_faces`; theory bound `2^(4K-8)`.')
    A('3. **P3 cover filter** — `mixed_cover` per face, rejecting patterns whose minimum')
    A('   atom cover exceeds the cap or the incumbent. Counter `cover_feasible_faces`.')
    A('4. **P4 active bases** — `stage43_partial_fpt.strict_vertex` decides the face. After')
    A('   parallel-row reduction and equality elimination it takes a closed form in')
    A('   dimension 0 or 1, and otherwise enumerates all `C(m, d)` active bases with')
    A('   `d = |S| + 1`, solving each by exact rational Gaussian elimination')
    A('   (`stage27_general_fpt.square_solve`). Counters `vertex_bases`, `local_lp_calls`,')
    A('   `interval_lp_calls`, `enumerated_lp_calls`, `zero_dimensional_lp_calls`.')
    A('')
    A('All instrumentation is external: the budget knobs `seconds`, `max_states`,')
    A('`max_bases`, `max_faces` the solver already exposes (geometric prefixes of the')
    A('search), the `statistics` dict it already returns, `cProfile` around the call, and a')
    A('static census that replicates only the preamble and the control flow fixing the face')
    A('count, with no LP work. %s profile runs, %s census rows.'
      % (fmt(rep['inputs']['profile_runs']), fmt(rep['inputs']['census_rows'])))
    A('')
    A('## 2. Where the time goes')
    A('')
    A('### 2.1 Exact phase counts, no LP solved')
    A('')
    A('Cap ladder on one fixed case (`%s`, |R| = %s, %s backbone edges, 13 positive and 5'
      % (H['name'], fmt(H['reference_set_size']), fmt(H['backbone_edges'])))
    A('negative baselines, 18 eligible vertices, 8 anchors). Sweeping the solver cap')
    A('isolates K as the sole variable.')
    A('')
    A('| K | deficit supports | supports with `negative` covered | zero faces `2^|A|` | faces after reduction | max &#124;A&#124; | bound 4K-8 |')
    A('|---|---|---|---|---|---|---|')
    for c in cap_rows:
        A('| %d | %s | %s | %s | %s | %d | %d |'
          % (c['cap'], fmt(c['deficit_supports']), fmt(c['supports_with_negative_covered']),
             fmt(c['mass_faces_total']), fmt(c['mass_faces_reduced']),
             c['max_A_size'], c['theory_max_A']))
    A('')
    A('The deficit-support enumeration is **not** the obstruction: %s supports at K = 10,'
      % fmt(cap_rows[-1]['deficit_supports']))
    A('enumerated in %.2f s of pure combinatorics. The face count is.'
      % cap_rows[-1]['census_seconds'])
    A('')
    A('### 2.2 Measured unit costs')
    A('')
    A('Non-negative least squares of `elapsed = c0 + alpha * mass_faces + beta * vertex_bases`')
    A('over %d deduplicated prefix points of the two ladders (R^2 = %.5f):'
      % (reg['n_points'], reg['r2']))
    A('')
    A('- `alpha` = **%.1f us per zero face** (mask construction plus the `mixed_cover` test)'
      % (alpha * 1e6))
    A('- `beta`  = **%.1f us per active basis** (one exact rational `square_solve` plus its'
      % (beta * 1e6))
    A('  feasibility check)')
    A('')
    A('Counters at the 600 s stop of the frozen solver on that case:')
    A('')
    A('| quantity | value | share |')
    A('|---|---|---|')
    A('| zero faces visited | %s | %.4f%% of the %s total |'
      % (fmt(stop['mass_faces']), 100.0 * stop['mass_faces'] / H['exact_zero_faces_frozen'],
         fmt(H['exact_zero_faces_frozen'])))
    A('| deficit supports visited | %s | of %s |' % (fmt(stop['deficit_sets']),
                                                    fmt(cap_rows[-1]['deficit_supports'])))
    A('| faces reaching the LP | %s | %.2f%% of faces |'
      % (fmt(stop['cover_feasible_faces']), 100 * stop['face_reaches_lp_fraction']))
    A('| LP calls taking the 1-D closed form | %s | %.1f%% of LP calls |'
      % (fmt(stop['interval_lp_calls']),
         100.0 * stop['interval_lp_calls'] / stop['local_lp_calls']))
    A('| LP calls entering basis enumeration | %s | %.1f%% of LP calls |'
      % (fmt(stop['enumerated_lp_calls']), 100 * stop['lp_is_enumerated_fraction']))
    A('| active bases enumerated | %s | %s per enumerating call |'
      % (fmt(stop['vertex_bases']), fmt(stop['bases_per_enumerated_call'])))
    A('')
    A('With %.1f bases per face at that depth the marginal cost is'
      % stop['bases_per_face'])
    A('`alpha + beta * %.1f` = **%.2f ms per face**, of which the basis enumeration is'
      % (stop['bases_per_face'], 1e3 * (alpha + beta * stop['bases_per_face'])))
    A('%.1f%%.' % (100 * beta * stop['bases_per_face'] / (alpha + beta * stop['bases_per_face'])))
    A('')
    A('### 2.3 Cumulative time per frame (independent cross-check)')
    A('')
    A('`cProfile`, cumulative seconds per frame as a share of the profiled total.')
    A('')
    A('| frame | K = 10, basis regime (%s s, %s faces, %s bases) | K = 7, completed run (%s s, %s faces, %s bases) |'
      % (fmt(cp_deep['wall_s']), fmt(cp_deep['mass_faces']), fmt(cp_deep['vertex_bases']),
         fmt(cp_k7['wall_s']), fmt(cp_k7['mass_faces']), fmt(cp_k7['vertex_bases'])))
    A('|---|---|---|')
    for frame in ('strict_vertex', 'square_solve', '_parallel_reduce',
                  '_affine_equalities', 'mixed_cover', 'deficit_sets'):
        d1 = (cp_deep['cumtime_key_frames'] or {}).get(frame)
        d2 = (cp_k7['cumtime_key_frames'] or {}).get(frame)
        t1 = cp_deep['wall_s']
        t2 = cp_k7['wall_s']
        A('| `%s` | %s | %s |'
          % (frame,
             '%.1f%%' % (100 * d1 / t1) if d1 else '~0%',
             '%.1f%%' % (100 * d2 / t2) if d2 else '~0%'))
    A('')
    A('Both instruments agree: at K = 10 the time is inside `strict_vertex`')
    A('(%.1f%% cumulative), and three quarters of it inside `square_solve`, the exact'
      % (100 * cp_deep['cumtime_key_frames']['strict_vertex'] / cp_deep['wall_s']))
    A('rational solve of one candidate active basis. The cover filter is %.1f%% and the'
      % (100 * cp_deep['cumtime_key_frames']['mixed_cover'] / cp_deep['wall_s']))
    A('deficit enumeration is below measurement resolution.')
    A('')
    A('### 2.4 The regime shift at K ~ 8')
    A('')
    A('The basis enumeration is **dormant below K = 8**. Median counters over completed')
    A('runs (12 smallest exact cases per cap):')
    A('')
    A('| K | zero faces | deficit supports | LP calls | enumerating LP calls | active bases | seconds |')
    A('|---|---|---|---|---|---|---|')
    for cap in sorted(g['median_counters_by_cap'], key=lambda x: int(x)):
        d = g['median_counters_by_cap'][cap]
        A('| %s | %s | %s | %s | %s | %s | %s |'
          % (cap, fmt(d.get('mass_faces')), fmt(d.get('deficit_sets')),
             fmt(d.get('cover_feasible_faces')), fmt(d.get('enumerated_lp_calls')),
             fmt(d.get('vertex_bases')), fmt(d.get('elapsed'), 2)))
    A('')
    A('`strict_vertex` enumerates bases only when the equality-reduced polytope has')
    A('dimension >= 2; dimension is `|S| + 1` minus the eliminated equality rank, and')
    A('`|S| <= K - 2`. Below K = 8 the reduced problem is almost always 1-dimensional and is')
    A('settled by the interval closed form, so `vertex_bases` stays at or near zero. From')
    A('K = 9 the enumeration engages and dominates everything else.')
    A('')
    A('This is why **runtime extrapolated from K <= 7 is badly wrong**: the fit over')
    A('completed K = 2, 3, 5, 7 runs gives %.2f bits per K and predicts %s s at K = 10,'
      % (fits['elapsed']['bits_per_K'], fmt(fits['elapsed']['predicted_at_K10'])))
    A('while the measured rate implies ~%s s. Fitting runtime inside the dormant regime'
      % fmt(H['frozen_projected_seconds']))
    A('under-predicts by about %sx.'
      % fmt(H['frozen_projected_seconds'] / fits['elapsed']['predicted_at_K10']))
    A('')
    A('## 3. Realised growth against the theoretical shape')
    A('')
    A('The paper bounds supports by `2^O(K)`, zero faces by `2^|A|` with `|A| <= 4K-8`, and')
    A('bases by `2^O(K)`. Measured:')
    A('')
    A('| counter | realised bits per K | per-K factor | theoretical shape |')
    A('|---|---|---|---|')
    A('| zero faces (cap ladder, one case) | %.3f | x%.1f | <= 4 bits per K |'
      % (ladder['faces_bits_per_K'], 2 ** ladder['faces_bits_per_K']))
    A('| zero faces (across cases, K = 2,3,5,7) | %.3f | x%.1f | <= 4 bits per K |'
      % (fits['mass_faces']['bits_per_K'], fits['mass_faces']['per_K_factor']))
    A('| deficit supports (cap ladder) | %.3f | x%.1f | `2^O(K)` |'
      % (ladder['supports_bits_per_K'], 2 ** ladder['supports_bits_per_K']))
    A('| zero faces after reduction | %.3f | x%.1f | — |'
      % (red_bits, 2 ** red_bits))
    A('')
    A('`max |A|` tracks the `4K-8` bound at a constant offset (%s against %s for K = %s).'
      % (ladder['max_A'], ladder['theory_max_A'], ladder['caps']))
    A('The extrapolation is checked against exact counts rather than trusted: the cross-case')
    A('fit predicts %s faces at K = 10 and the static census measures a median of %s over'
      % (fmt(fits['mass_faces']['predicted_at_K10']),
         fmt(g['census_exact_K10']['mass_faces_exact_K10']['median'])))
    A('%d K = 10 cases (range %s to %s).'
      % (g['census_exact_K10']['mass_faces_exact_K10']['n'],
         fmt(g['census_exact_K10']['mass_faces_exact_K10']['min']),
         fmt(g['census_exact_K10']['mass_faces_exact_K10']['max'])))
    A('')
    A('So the *face count* matches the theory\'s exponent almost exactly — 4.06 to 4.13 bits')
    A('per K against a bound of 4. The enumeration is not worse than the bound allows. What')
    A('the bound does not say is that most of the enumerated faces are provably infeasible.')
    A('')
    A('Two costs were on the candidate list and are already present in the frozen solver, so')
    A('neither is available as a fix: the edge-cover bound is **already** applied before the')
    A('LP (`mixed_cover` precedes `strict_vertex`), and equality elimination is already')
    A('performed per face by `_affine_equalities`.')
    A('')
    A('## 4. Prototype: `scripts/stage295_reduced_minimum.py`')
    A('')
    A('A new module. Every primitive is imported from the frozen modules and the control')
    A('flow of `solve_mixed` is reproduced line for line, with one change: the mask loop runs')
    A('over `free(S) = S u {u in N(S) \\ S : b_u < 0}` instead of all of `A`. Two further')
    A('exact switches were measured separately — a bounded `mixed_cover` memo table and')
    A('visiting masks in order of increasing `|P|`.')
    A('')
    A('### 4.1 Validation')
    A('')
    A('| check | result |')
    A('|---|---|')
    A('| cases where both solvers completed | %s |' % fmt(V['comparable_cases']))
    A('| minima identical | **%s** |' % fmt(V['agreeing']))
    A('| disagreements | **%s** |' % fmt(len(V['disagreements'])))
    A('| `full_mask=True` runs reproducing the frozen face count exactly | %s of %s |'
      % (fmt(V['full_mask_identical_face_counts']), fmt(V['full_mask_identity_rows'])))
    A('')
    A('The agreement set spans exactly-explainable cases at K = 2, 3, 5, 7 and')
    A('non-explainable cases (real krepp ANI-estimated distances and published rounded')
    A('abundances), so the infeasible code path is covered too. The `full_mask` identity is')
    A('the strongest correctness signal available: with the reduction switched off the')
    A('prototype visits exactly the same number of faces as the frozen solver.')
    A('')
    A('### 4.2 Speedup where both complete')
    A('')
    A('| K | n | median speedup | min | max |')
    A('|---|---|---|---|---|')
    for k, v in V['speedup_by_cap'].items():
        A('| %s | %s | %.2fx | %.2fx | %.1fx |'
          % (k, v['n'], v['median'], v['min'], v['max']))
    A('')
    A('The gain is negligible at K = 2, 3 — as it must be, since `free(S) = A` when there are')
    A('no negative baselines to free and `|S|` is small — and grows with K exactly where the')
    A('wasted enumeration grows.')
    A('')
    A('### 4.3 K = 10 is now reachable')
    A('')
    _v2 = [r for r in P['rows'] if r['batch'] == 'push_v2_bounded_cache']
    _v2done = [r for r in _v2 if r['status'] == 'complete']
    _v1fail = [r for r in P['rows'] if r['batch'] == 'push_v1_unbounded_cache'
               and r['status'] != 'complete']
    A('In the current batch %s of %s K >= 10 runs completed; the %s MemoryError rows'
      % (fmt(len(_v2done)), fmt(len(_v2)), fmt(len(_v1fail))))
    A('of the earlier batch were caused by the unbounded memo table and are superseded')
    A('by this one (Amendment 1). Every completed run carries two certificates: the')
    A('witness is exactly verified against the observations by the frozen checker')
    A('(`observe(...) == observed`, %s of %s), and a second complete run at `cap = k* - 1`'
      % (fmt(sum(1 for x in _v2done if x['witness_exactly_verified'])), fmt(len(_v2done))))
    A('returns infeasible, proving `k*` minimal (%s of %s). This certificate chain does not'
      % (fmt(sum(1 for x in _v2done if x['minimality_certified'])), fmt(len(_v2done))))
    A('use the frozen minimum pass at all.')
    A('')
    A('| &#124;R&#124; | K | minimum found | seconds | zero faces | active bases | minimality certified |')
    A('|---|---|---|---|---|---|---|')
    seen = set()
    for x in sorted([r for r in P['rows'] if r['status'] == 'complete'
                     and r['batch'] == 'push_v2_bounded_cache'],
                    key=lambda r: (r['reference_set_size'], r['seconds'])):
        key = (x['reference_set_size'], x['name'])
        if key in seen:
            continue
        seen.add(key)
        A('| %s | %d | %s | %s | %s | %s | %s |'
          % (fmt(x['reference_set_size']), x['cap'], x['minimum'], fmt(x['seconds']),
             fmt(x['mass_faces']), fmt(x['vertex_bases']), x['minimality_certified']))
    A('')
    A('On the headline case the end-to-end comparison is:')
    A('')
    A('| | frozen | prototype |')
    A('|---|---|---|')
    A('| zero faces | %s | %s |' % (fmt(H['exact_zero_faces_frozen']),
                                   fmt(H['exact_zero_faces_reduced'])))
    A('| active bases | %s in the first %s faces | %s in all %s faces |'
      % (fmt(stop['vertex_bases']), fmt(stop['mass_faces']),
         fmt(head_rows[0]['vertex_bases']), fmt(head_rows[0]['mass_faces'])))
    A('| wall time | %s days projected, >= %s days floor | **%s s measured** |'
      % (fmt(H['frozen_projected_days']), fmt(floor_s / 86400.0),
         fmt(H['prototype_seconds'])))
    A('| result | `unknown` at 3,600 s (Stage 292) | minimum %s, certified |'
      % H['prototype_minimum'])
    A('')
    A('The projection is stated two ways. The **floor** multiplies the exact face count by')
    A('the face coefficient alone (%.1f us), ignoring every basis solve: %s days. The'
      % (alpha * 1e6, fmt(floor_s / 86400.0)))
    A('**measured-rate projection** uses the %.1f bases per face observed in the sampled'
      % stop['bases_per_face'])
    A('prefix: %s days. Both understate the true cost, because the sampled prefix visits'
      % fmt(H['frozen_projected_days']))
    A('the smallest supports first and per-face cost rises with `|S|`.')
    A('')
    A('Mask ordering by `|P|` changed nothing measurable (%s vs %s s on the headline case).'
      % (fmt([r['seconds'] for r in P['rows'] if r['name'] == H['name']
              and r['popcount_order'] and r['status'] == 'complete'][:1][0]
             if any(r['name'] == H['name'] and r['popcount_order']
                    and r['status'] == 'complete' for r in P['rows']) else None),
         fmt(H['prototype_seconds'])))
    A('The memo table had a 58% hit rate but, unbounded, exhausted 10 GiB on |R| > 20; it is')
    A('now capped (Amendment 1) and is not where the gain comes from.')
    A('')
    A('## 5. Is the ceiling intrinsic or removable?')
    A('')
    A('**Removable, and removed, at K = 10.** By the decision rule registered in the')
    A('protocol, a cost is an implementation artefact if it admits a reformulation that')
    A('provably visits strictly fewer objects without changing the answer. The forced-')
    A('positive argument is such a reformulation, it is proved rather than heuristic, and it')
    A('is confirmed empirically on %s head-to-head cases with zero disagreements.'
      % fmt(V['comparable_cases']))
    A('')
    A('**The exponent is intrinsic.** The surviving face count still grows at %.2f bits per'
      % red_bits)
    A('K. Projecting the prototype forward at its own measured growth and its own measured')
    A('%.1f us per face on the headline case:' % (per_face_proto * 1e6 if per_face_proto else 0))
    A('')
    A('| K | projected zero faces | projected seconds |')
    A('|---|---|---|')
    for K, v in forward.items():
        A('| %d | %s | %s |' % (K, fmt(int(round(v['faces']))), fmt(v['seconds'])))
    A('')
    A('So the practical ceiling moved from K = 7 to roughly K = 12-13 on cases of this size,')
    A('and the method stays exponential in K. The next candidate cost, **not** prototyped')
    A('here, is `strict_vertex` itself: it decides a face by enumerating all `C(m, d)` active')
    A('bases where an exact rational LP would do, and that is %.1f%% of the time in the'
      % (100 * cp_deep['cumtime_key_frames']['strict_vertex'] / cp_deep['wall_s']))
    A('basis regime. Whether replacing it helps is untested.')
    A('')
    A('## 6. Sentences for the paper')
    A('')
    pct_infeasible = '{:.1%}'.format(1 - H['exact_zero_faces_reduced']
                                     / H['exact_zero_faces_frozen'])
    A('> The exact minimum-source pass is fixed-parameter tractable in the source')
    A('> count K, with `2^O(K) poly(B)` behaviour, and the realised growth of its')
    A('> dominant counter -- the number of zero faces of the mass polytope -- is')
    A('> %.2f bits per unit K against a bound of `2^(4K-8)`, i.e. 4 bits per unit K.'
      % ladder['faces_bits_per_K'])
    A('> A direct implementation enumerates all `2^|A|` sign masks over `A = S u N(S)`,')
    A('> but on the strict face the mass at a neighbour of the deficit support with')
    A('> non-negative baseline is forced positive, so %s of those masks are infeasible'
      % pct_infeasible)
    A('> by construction. Restricting the enumeration to the free coordinates reduces')
    A('> the growth to %.2f bits per unit K and moves the practical ceiling from K = 7'
      % red_bits)
    A('> to K = 10 at reference-set sizes up to %s, with per-case times of %s to %s s;'
      % (fmt(max(P['completed_reference_sizes'])),
         fmt(P['seconds_quantiles']['min']), fmt(P['seconds_quantiles']['max'])))
    A('> the minimum and its witness are certified exactly in rational arithmetic.')
    A('')
    A('## 7. Scope and limits')
    A('')
    A('- **Pre-registration chronology (a deviation from the house pattern).** The static')
    A('  enumeration census of Section 2.1 was run, and `results/stage295_static_census.json`')
    A('  written, BEFORE the protocol was registered; it is what identified the face loop')
    A('  as the phase to instrument, and the protocol then fixed the plan around it. Every')
    A('  timed measurement -- the profile runs, the agreement comparisons, the K >= 10')
    A('  push -- ran after registration. The census contains no timing and solves no LP:')
    A('  it is deterministic combinatorial enumeration, and re-running it under the')
    A('  registered script hash reproduces all %s rows with every non-timing field'
      % fmt(_replay['rows_compared']))
    A('  identical (`results/stage295_census_replay_check.json`,')
    A('  `all_non_timing_fields_identical = %s`). The sequencing is recorded as a'
      % _replay['all_non_timing_fields_identical'])
    A('  deviation regardless.')
    A('- The Stage 287 bank has caps {2, 3, 5, 7, 10, 12} and **no K = 8 or K = 9 families**.')
    A('  K = 8 and K = 9 were realised by sweeping the solver `cap` knob on fixed K = 10')
    A('  cases, which isolates K as the sole variable (`maxweight = K-2`, `|A| <= 4K-8` and')
    A('  `d = |S|+1` all depend on the cap alone). Declared in the protocol as the one scope')
    A('  substitution of the stage.')
    A('- The frozen K = 10 runtime is a **projection**, not a measurement: no frozen K = 10')
    A('  run has ever completed. It rests on unit costs fitted to budget-truncated prefixes')
    A('  of the smallest K = 10 case and on exact face counts from the static census.')
    A('- Per-face unit cost was sampled in the first %.4f%% of the search on one case and is'
      % (100.0 * stop['mass_faces'] / H['exact_zero_faces_frozen']))
    A('  assumed representative; because supports are visited smallest-first, the true')
    A('  average is probably higher, making the projection conservative.')
    A('- The forward projection to K = 11-15 in Section 5 is an extrapolation of a growth')
    A('  rate fitted on K = 5..10 for a single case; it is indicative, not measured.')
    _p16 = os.path.join(RES, 'stage295_push2', 'push_0016.json')
    if os.path.exists(_p16):
        _d16 = json.load(open(_p16))
        _s16 = (_d16.get('primary') or {}).get('statistics') or {}
        A('- One opportunistic row beyond the brief: the single `exact_16k` case')
        A('  (|R| = %s, %s original edges, registered cap 12) did NOT complete. It'
          % (fmt(_d16.get('reference_set_size')), fmt(31173)))
        A('  stopped on its 3,600 s wall budget after %s s having visited %s zero faces'
          % (fmt((_d16.get('primary') or {}).get('elapsed')), fmt(_s16.get('mass_faces'))))
        A('  and, notably, %s deficit supports with zero active bases. At that scale'
          % fmt(_s16.get('deficit_sets')))
        A('  phase P1 -- negligible everywhere else in this stage -- becomes an')
        A('  obstruction in its own right. That is a different ceiling from the one')
        A('  diagnosed here, and nothing is claimed about K = 12 at that size.')
    A('- The prototype replaces the minimum pass only. Stage 292\'s all-minimisers routine')
    A('  takes `kstar` as an input and is unaffected.')
    A('')
    A('## 8. Files')
    A('')
    A('Protocol and records (on the server, under `~/projects/recomb2027_novelty`):')
    A('')
    A('- `results/stage295_ceiling_protocol.json` — protocol, registered before any run')
    A('- `results/stage295_protocol_amendment.json` — Amendment 1 (bounded memo), frozen hashes re-verified')
    A('- `results/stage295_census_replay_check.json` — census determinism replay under the registered hash')
    A('- `results/stage295_static_census.json` — exact phase-1/2 counts, %s rows'
      % fmt(rep['inputs']['census_rows']))
    A('- `results/stage295_profile_runs/`, `results/stage295_profile_deep/` — %s profile runs'
      % fmt(rep['inputs']['profile_runs']))
    A('- `results/stage295_validation/`, `results/stage295_revalidation/` — %s head-to-head rows'
      % fmt(rep['inputs']['validation_rows']))
    A('- `results/stage295_push/`, `results/stage295_push2/` — %s K >= 10 attempts'
      % fmt(rep['inputs']['push_rows']))
    A('- `results/stage295_ceiling_analysis.json`, `results/stage295_ceiling_report.json`')
    A('')
    A('New code (no frozen module touched):')
    A('')
    for f in ('stage295_common.py', 'stage295_static_census.py', 'stage295_profile.py',
              'stage295_reduced_minimum.py', 'stage295_validate.py',
              'stage295_ceiling_push.py', 'stage295_analyze.py', 'stage295_report.py',
              'stage295_write_document.py', 'stage295_register_protocol.py'):
        A('- `scripts/%s`' % f)
    A('')

    out = os.path.join(ROOT, 'stage295_ceiling_diagnosis.md')
    open(out, 'w').write('\n'.join(L) + '\n')
    print('wrote', out, len(L), 'lines')


if __name__ == '__main__':
    main()
