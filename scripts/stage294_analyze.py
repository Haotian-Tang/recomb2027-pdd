#!/usr/bin/env python3
"""Stage 294 analysis: the CHANGE reference removal induces on a noisy inference.

Two steps.

  metrics  per run, against the bank's published truth: the returned support,
           how much of the generating support it recovers, the tree-metric
           Wasserstein-1 distance to the generating mixture and the
           abundance-weighted placement error.  W1 and the cost matrix are
           stage289_headtohead_analyze.wasserstein / point_distances, imported
           unchanged, so the numbers sit on the stage 289 scale.

  summary  the count-change distributions.  The reference point is NOT the
           published truth -- DecoDiPhy already mis-estimates k on almost every
           full-panel krepp case -- but the SAME arm on the SAME case at level
           "none".  A draw contributes only when both it and its own full-panel
           baseline returned an answer; nothing is imputed.

Also computed here: the exact reference numbers recomputed from the registered
stage 288 analysis for side-by-side comparison, the paired exact_matched arm on
the same draws, and how close whole-clade removal is to the taxonomic regime.
"""
import os
for _key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
             'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_key] = '1'
import argparse
import collections
import json
import multiprocessing as mp
import statistics
import sys
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import stage294_common as s294                                    # noqa: E402
from stage289_headtohead_analyze import (float_adjacency, distances_from,      # noqa: E402
                                         point_distances, wasserstein)

OUTDIR = ROOT / 'results/stage294'
PROTOCOL = OUTDIR / 'stage294_protocol.json'
DD_STREAM = OUTDIR / 'stage294_dd_runs.jsonl'
EXACT_STREAM = OUTDIR / 'stage294_exact_matched_runs.jsonl'
METRICS = OUTDIR / 'stage294_dd_metrics.jsonl'
SUMMARY = OUTDIR / 'stage294_analysis.json'
STAGE288_ANALYSIS = ROOT / 'results/stage288_analysis.json'

KEY = ('name', 'level', 'regime', 'replicate')
LEVELS = ('none', 'one', '0.10', '0.25', '0.50', '0.75')
REGIMES = ('random', 'clade', 'nearest', 'taxon')


def key_of(row):
    return tuple(row[k] for k in KEY)


def quantiles(values):
    if not values:
        return None
    values = sorted(values)
    def q(p):
        return values[min(len(values) - 1, int(p * len(values)))]
    return dict(n=len(values), min=values[0], q10=q(0.10), q25=q(0.25),
                median=statistics.median(values), q75=q(0.75), q90=q(0.90), max=values[-1],
                mean=sum(values) / len(values))


# ----------------------------------------------------------------- metrics

def case_metrics(payload):
    """Truth-referenced metrics for every DecoDiPhy run of one case."""
    entry, rows = payload
    case = s294.load_noisy_case(entry)
    edges = case['edges']
    adjacency = float_adjacency(edges)
    truth = [dict(edge=s['edge'], p=float(s['p']), x=float(s['x']))
             for s in case['generating_sources']]
    true_support = {s['edge'] for s in truth}
    needed = {v for s in truth for v in edges[s['edge']][:2]}
    for row in rows:
        for s in row['result'].get('sources') or []:
            needed.update(edges[s['edge']][:2])
    table = {v: distances_from(adjacency, v) for v in needed}
    out = []
    for row in rows:
        result = row['result']
        record = {k: row[k] for k in ('name', 'arm', 'level', 'regime', 'replicate',
                                      'K_generating', 'n_leaves', 'n_removed', 'n_retained',
                                      'removed_fraction', 'reference_set_size')}
        record.update(status=result.get('status'),
                      jplace_writer_crashed=bool(result.get('jplace_writer_crashed')),
                      stop_reason=result.get('stop_reason'),
                      loss=result.get('loss'),
                      wall_seconds=result.get('process_wall_seconds'),
                      exception=(result.get('exception') or result.get('reason')),
                      detail_rank=(row.get('detail') or {}).get('rank'),
                      detail_order=(row.get('detail') or {}).get('order'))
        if result.get('status') == 'complete':
            sources = result['sources']
            support = [s['edge'] for s in sources]
            cost = point_distances(edges, adjacency, table, sources, truth)
            weights = [float(s['p']) for s in sources]
            w1 = None
            if weights and abs(sum(weights) - 1.0) < 1e-6:
                w1 = wasserstein(cost, weights, [s['p'] for s in truth])
            record.update(
                k=result['k'], support=sorted(support),
                distinct_support=len(set(support)),
                repeated_edge=len(set(support)) != len(support),
                support_exact_match=set(support) == true_support,
                n_true_edges_recovered=len(set(support) & true_support),
                true_support_fraction=len(set(support) & true_support) / len(true_support),
                abundance_mass_sum=sum(weights),
                w1_to_truth=w1,
                weighted_placement_error=sum(
                    t['p'] * min(cost[i][j] for i in range(len(sources)))
                    for j, t in enumerate(truth)),
                boundary_placements=sum(s['x'] <= 1e-9 or s['x'] >= 1 - 1e-9 for s in sources))
        out.append(record)
    return out


def run_metrics(workers):
    assert not METRICS.exists(), 'Never overwrite computed metrics'
    entries = s294.bank_entries(s294.NOISY_FAMILY)
    rows = collections.defaultdict(list)
    with DD_STREAM.open() as handle:
        for line in handle:
            row = json.loads(line)
            rows[row['name']].append(row)
    payloads = [(entries[name], rows[name]) for name in sorted(rows)]
    done = 0
    with METRICS.open('a') as out:
        with mp.get_context('fork').Pool(workers) as pool:
            for records in pool.imap_unordered(case_metrics, payloads, chunksize=1):
                for record in records:
                    out.write(json.dumps(record) + '\n')
                out.flush()
                done += 1
                if done % 50 == 0:
                    print(f'{done}/{len(payloads)} cases', flush=True)
    print('metrics written', METRICS)


# ----------------------------------------------------------------- summary

def paired_changes(rows, value_key, complete):
    """delta = value(draw) - value(full panel), same arm, same case."""
    base = {}
    for row in rows:
        if row['level'] == 'none' and complete(row):
            base[(row['name'], row['arm'])] = row[value_key]
    deltas = []
    for row in rows:
        if row['level'] == 'none':
            continue
        anchor = base.get((row['name'], row['arm']))
        if anchor is None or not complete(row):
            continue
        deltas.append(dict(name=row['name'], arm=row['arm'], level=row['level'],
                           regime=row['regime'], replicate=row['replicate'],
                           K=row.get('K_generating'), n_removed=row['n_removed'],
                           removed_fraction=row['removed_fraction'],
                           baseline=anchor, value=row[value_key],
                           delta=row[value_key] - anchor))
    return deltas, base


def change_block(deltas):
    values = [d['delta'] for d in deltas]
    if not values:
        return dict(comparable=0)
    return dict(comparable=len(values),
                fell=sum(1 for v in values if v < 0),
                unchanged=sum(1 for v in values if v == 0),
                rose=sum(1 for v in values if v > 0),
                fraction_fell=sum(1 for v in values if v < 0) / len(values),
                fraction_rose=sum(1 for v in values if v > 0) / len(values),
                mean_delta=sum(values) / len(values),
                distribution=dict(collections.Counter(values)),
                quantiles=quantiles(values))


def by_cell(deltas):
    out = dict(overall=change_block(deltas), by_level={}, by_regime={}, by_level_regime={},
               by_K={})
    for level in LEVELS[1:]:
        subset = [d for d in deltas if d['level'] == level]
        out['by_level'][level] = change_block(subset)
    for regime in REGIMES:
        subset = [d for d in deltas if d['regime'] == regime]
        out['by_regime'][regime] = change_block(subset)
        for level in LEVELS[1:]:
            cell = [d for d in subset if d['level'] == level]
            out['by_level_regime'][f'{regime}|{level}'] = change_block(cell)
    for k in sorted({d['K'] for d in deltas if d['K'] is not None}):
        out['by_K'][str(k)] = change_block([d for d in deltas if d['K'] == k])
    random_fraction = out['by_regime'].get('random', {}).get('fraction_fell')
    out['structured_over_random_ratio'] = {
        regime: (out['by_regime'][regime]['fraction_fell'] / random_fraction
                 if random_fraction else None)
        for regime in REGIMES if regime != 'random' and out['by_regime'].get(regime, {}).get('comparable')}
    return out


def availability(rows):
    out = {}
    for arm in sorted({r['arm'] for r in rows}):
        block = dict(by_level={}, by_regime={}, by_level_regime={})
        for label, selector in (('by_level', 'level'), ('by_regime', 'regime')):
            for value in sorted({r[selector] for r in rows}):
                subset = [r for r in rows if r['arm'] == arm and r[selector] == value]
                if not subset:
                    continue
                block[label][value] = dict(
                    runs=len(subset),
                    complete=sum(1 for r in subset if r['status'] == 'complete'),
                    fraction_complete=sum(1 for r in subset if r['status'] == 'complete') / len(subset),
                    recovered_after_writer_crash=sum(
                        1 for r in subset if r['status'] == 'complete' and r['jplace_writer_crashed']),
                    statuses=dict(collections.Counter(r['status'] for r in subset)))
        for regime in REGIMES:
            for level in LEVELS[1:]:
                subset = [r for r in rows if r['arm'] == arm and r['regime'] == regime
                          and r['level'] == level]
                if subset:
                    block['by_level_regime'][f'{regime}|{level}'] = dict(
                        runs=len(subset),
                        complete=sum(1 for r in subset if r['status'] == 'complete'),
                        fraction_complete=sum(1 for r in subset if r['status'] == 'complete') / len(subset))
        out[arm] = block
    return out


def accuracy_block(rows):
    """Truth-referenced accuracy, for the record; not the primary quantity."""
    out = {}
    for arm in sorted({r['arm'] for r in rows}):
        block = {}
        for level in LEVELS:
            subset = [r for r in rows if r['arm'] == arm and r['level'] == level
                      and r['status'] == 'complete']
            if not subset:
                continue
            w1 = [r['w1_to_truth'] for r in subset if r.get('w1_to_truth') is not None]
            block[level] = dict(
                runs=len(subset),
                support_exact_match=sum(1 for r in subset if r['support_exact_match']),
                mean_true_support_fraction=sum(r['true_support_fraction'] for r in subset) / len(subset),
                w1_to_truth=quantiles(w1),
                weighted_placement_error=quantiles([r['weighted_placement_error'] for r in subset]),
                k=quantiles([r['k'] for r in subset]),
                k_equals_k_true=sum(1 for r in subset if r['k'] == r['K_generating']))
        out[arm] = block
    return out


def taxon_proxy():
    """How close is whole-clade removal to dropping whole taxonomic groups?"""
    ranks = s294.load_ranks()
    plan = s294.draw_plan(ranks)
    by_case_level = collections.defaultdict(dict)
    cases = {}
    for noisy, _exact, meta in plan:
        cases[meta['name']] = noisy
        by_case_level[(meta['name'], meta['level'])][(meta['regime'], meta['replicate'])] = meta
    jaccard = collections.defaultdict(list)
    mono = collections.Counter()
    rank_used = collections.Counter()
    achieved = collections.defaultdict(list)
    bipart_cache = {}
    for (name, level), draws in by_case_level.items():
        if level == 'none':
            continue
        taxa = [m for (r, _), m in draws.items() if r == 'taxon']
        if not taxa:
            continue
        case = cases[name]
        leaves = sorted(case['observed'])
        if name not in bipart_cache:
            bipart_cache[name] = bipartitions(case['edges'], leaves)
        parts = bipart_cache[name]
        for meta in taxa:
            rank_used[meta['detail']['rank']] += 1
            achieved[level].append(meta['n_removed'] / max(1, meta['detail']['target_count']))
            removed = frozenset(meta['removed'])
            mono['taxon_draw_is_a_clade' if removed in parts else 'taxon_draw_is_not_a_clade'] += 1
            groups = s294.taxon_groups(ranks, case['leaf_labels'], leaves, meta['detail']['rank'])
            for group in meta['detail']['groups']:
                mono['group_monophyletic' if frozenset(groups[group]) in parts
                     else 'group_not_monophyletic'] += 1
            for other_regime in ('clade', 'nearest', 'random'):
                best = 0.0
                for (regime, _rep), other in draws.items():
                    if regime != other_regime:
                        continue
                    a, b = set(meta['removed']), set(other['removed'])
                    best = max(best, len(a & b) / len(a | b))
                jaccard[f'{other_regime}|{level}'].append(best)
                jaccard[other_regime].append(best)
    return dict(rank_used=dict(rank_used), monophyly=dict(mono),
                achieved_over_target_by_level={k: quantiles(v) for k, v in achieved.items()},
                best_jaccard_to_other_regime={k: quantiles(v) for k, v in jaccard.items()})


def bipartitions(edges, leaves):
    """Every leaf-set induced by removing one edge, as frozensets."""
    adjacency = collections.defaultdict(list)
    for index, (u, v, _) in enumerate(edges):
        adjacency[u].append((v, index))
        adjacency[v].append((u, index))
    allset = set(leaves)
    parts = set()
    for index, (u, _v, _length) in enumerate(edges):
        seen = {u}
        stack = [u]
        while stack:
            x = stack.pop()
            for y, eidx in adjacency[x]:
                if eidx == index or y in seen:
                    continue
                seen.add(y)
                stack.append(y)
        side = frozenset(seen & allset)
        parts.add(side)
        parts.add(frozenset(allset - side))
    return parts


def stage289_agreement(rows):
    """Reproducibility: at level "none" this stage repeats the stage 289 full-panel run.

    Same tool, same arm, same case, same documented defaults, same unremoved
    panel -- so k and the returned support must agree exactly.  Disagreement
    would mean the stage 294 bridge differs from the stage 289 one; agreement is
    recorded rather than assumed.
    """
    path = ROOT / 'results/stage289/stage289_headtohead_runs.jsonl'
    if not path.exists():
        return None
    earlier = {}
    with path.open() as handle:
        for line in handle:
            row = json.loads(line)
            if (row['family'] != s294.NOISY_FAMILY
                    or row['method'] not in ('dd_default', 'dd_true_k')
                    or row['repetition'] != 1
                    or row['result'].get('status') != 'complete'):
                continue
            earlier[(row['name'], row['method'])] = dict(
                k=row['result']['k'], loss=row['result']['loss'],
                support=sorted(s['edge'] for s in row['result']['sources']))
    out = collections.Counter()
    mismatches = []
    for row in rows:
        if row['level'] != 'none' or row['status'] != 'complete':
            continue
        before = earlier.get((row['name'], row['arm']))
        if before is None:
            out['no_stage289_counterpart'] += 1
            continue
        same_k = before['k'] == row['k']
        same_support = before['support'] == row['support']
        out['compared'] += 1
        out['k_agrees'] += same_k
        out['support_agrees'] += same_support
        if not (same_k and same_support) and len(mismatches) < 20:
            mismatches.append(dict(name=row['name'], arm=row['arm'], stage289=before,
                                   stage294=dict(k=row['k'], support=row['support'],
                                                 loss=row['loss'])))
    return dict(counts=dict(out), mismatches=mismatches)


def summarise():
    protocol = json.loads(PROTOCOL.read_text())
    rows = [json.loads(line) for line in METRICS.open()]
    for row in rows:
        row['K_generating'] = row['K_generating']
    complete = lambda r: r['status'] == 'complete' and r.get('k') is not None
    dd = {}
    for arm in ('dd_default', 'dd_true_k'):
        subset = [r for r in rows if r['arm'] == arm]
        deltas, base = paired_changes(subset, 'k', complete)
        dd[arm] = dict(cells=by_cell(deltas),
                       baseline_cases_with_an_answer=len(base),
                       baseline_k=quantiles([v for v in base.values()]))
    exact = None
    if EXACT_STREAM.exists():
        erows = []
        for line in EXACT_STREAM.open():
            row = json.loads(line)
            erows.append(dict(name=row['name'], arm='exact_matched', level=row['level'],
                              regime=row['regime'], replicate=row['replicate'],
                              K_generating=row['K_generating'], n_removed=row['n_removed'],
                              n_retained=row['n_retained'],
                              removed_fraction=row['removed_fraction'],
                              status=row['result']['status'],
                              minimum=row['result'].get('minimum'),
                              support_exact_match=row['result'].get('support_exact_match'),
                              jplace_writer_crashed=False))
        ecomplete = lambda r: r['status'] == 'complete' and r.get('minimum') is not None
        deltas, base = paired_changes(erows, 'minimum', ecomplete)
        kby = {r['name']: r['K_generating'] for r in erows}
        exact = dict(cells=by_cell(deltas), runs=len(erows),
                     availability=availability(erows),
                     baseline_cases_with_an_answer=len(base),
                     baseline_minimum_equals_K=sum(1 for (name, _arm), value in base.items()
                                                   if kby.get(name) == value),
                     baseline_minimum=quantiles(list(base.values())))
        # direction agreement on draws where both the exact and the noisy arm answered
        dd_delta = {(d['name'], d['level'], d['regime'], d['replicate']): d['delta']
                    for d in paired_changes([r for r in rows if r['arm'] == 'dd_default'],
                                            'k', complete)[0]}
        table = collections.Counter()
        for d in deltas:
            key = (d['name'], d['level'], d['regime'], d['replicate'])
            if key in dd_delta:
                table[(('exact_fell' if d['delta'] < 0 else
                        'exact_same' if d['delta'] == 0 else 'exact_rose'),
                       ('noisy_fell' if dd_delta[key] < 0 else
                        'noisy_same' if dd_delta[key] == 0 else 'noisy_rose'))] += 1
        exact['paired_direction_table'] = {f'{a}|{b}': n for (a, b), n in sorted(table.items())}
    stage288 = json.loads(STAGE288_ANALYSIS.read_text())['incompleteness']['undercount_B_versus_A']
    summary = dict(
        stage=294, protocol_sha256=s294.sha(PROTOCOL),
        plan_sha256=protocol['plan_sha256'],
        runs=dict(dd=sum(1 for _ in DD_STREAM.open()),
                  exact_matched=(sum(1 for _ in EXACT_STREAM.open())
                                 if EXACT_STREAM.exists() else 0)),
        availability=availability(rows),
        count_change=dd,
        accuracy_against_truth=accuracy_block(rows),
        exact_matched=exact,
        stage288_exact_reference=stage288,
        taxon_regime=taxon_proxy(),
        stage289_full_panel_agreement=stage289_agreement(rows))
    SUMMARY.write_text(json.dumps(summary, indent=1) + '\n')
    print(json.dumps(dict(
        dd_default_overall=dd['dd_default']['cells']['overall'],
        dd_default_by_regime={k: v['fraction_fell']
                              for k, v in dd['dd_default']['cells']['by_regime'].items()},
        exact_matched_overall=(exact or {}).get('cells', {}).get('overall'),
        stage288_fraction_B_below_A=stage288['B_below_A'] / stage288['comparable']), indent=1))




# ------------------------------------------------- matched exact/noisy view

def matched_comparison():
    """The exact and the noisy answer on the SAME draw, on one common draw set.

    Three things make the headline by-regime numbers of the two arms not
    directly comparable, and all three are removed here:
      * the exact arm runs only at K in {3, 5}, the noisy arm at K in {3, 5, 10};
      * DecoDiPhy returns no answer at all on a large and REGIME-DEPENDENT share
        of the draws, so a rate computed over its completed runs is conditioned
        on survival; here a failure is kept as its own outcome, 'no_answer';
      * a draw counts only when the arm's own full-panel baseline exists.
    """
    rows = [json.loads(line) for line in METRICS.open() if '"dd_default"' in line]
    dd = {}
    base_dd = {}
    for row in rows:
        if row['arm'] != 'dd_default':
            continue
        if row['level'] == 'none':
            if row['status'] == 'complete':
                base_dd[row['name']] = row['k']
            continue
        dd[key_of(row)] = row
    exact, base_exact = {}, {}
    for line in EXACT_STREAM.open():
        row = json.loads(line)
        value = row['result'].get('minimum')
        if row['level'] == 'none':
            if row['result']['status'] == 'complete' and value is not None:
                base_exact[row['name']] = value
            continue
        exact[key_of(row)] = row
    def direction(delta):
        return 'fell' if delta < 0 else 'rose' if delta > 0 else 'same'
    cells = collections.defaultdict(lambda: collections.Counter())
    agreement = collections.defaultdict(lambda: collections.Counter())
    deltas = collections.defaultdict(list)
    for key, erow in exact.items():
        name = key[0]
        if name not in base_exact or name not in base_dd:
            continue
        if erow['result']['status'] != 'complete' or erow['result'].get('minimum') is None:
            continue
        drow = dd.get(key)
        if drow is None:
            continue
        level, regime = erow['level'], erow['regime']
        e_delta = erow['result']['minimum'] - base_exact[name]
        cells[(regime, level)]['draws'] += 1
        cells[(regime, level)]['exact_' + direction(e_delta)] += 1
        deltas[('exact', regime, level)].append(e_delta)
        if drow['status'] != 'complete' or drow.get('k') is None:
            cells[(regime, level)]['noisy_no_answer'] += 1
            agreement[(regime, level)][f'{direction(e_delta)}|no_answer'] += 1
            continue
        n_delta = drow['k'] - base_dd[name]
        cells[(regime, level)]['noisy_' + direction(n_delta)] += 1
        deltas[('noisy', regime, level)].append(n_delta)
        agreement[(regime, level)][f'{direction(e_delta)}|{direction(n_delta)}'] += 1
    def block(counter):
        total = counter['draws']
        answered = total - counter['noisy_no_answer']
        out = dict(draws=total, **{k: v for k, v in counter.items() if k != 'draws'})
        out['exact_fraction_fell'] = counter['exact_fell'] / total if total else None
        out['noisy_fraction_fell_of_answered'] = (counter['noisy_fell'] / answered
                                                  if answered else None)
        out['noisy_fraction_fell_of_all_draws'] = counter['noisy_fell'] / total if total else None
        out['noisy_fraction_no_answer'] = counter['noisy_no_answer'] / total if total else None
        return out
    by_regime_level = {f'{r}|{l}': block(v) for (r, l), v in sorted(cells.items())}
    by_regime, by_level, overall = {}, {}, collections.Counter()
    for (regime, level), counter in cells.items():
        for scope, store in ((regime, by_regime), (level, by_level)):
            acc = store.setdefault(scope, collections.Counter())
            acc.update(counter)
        overall.update(counter)
    agg = collections.Counter()
    for counter in agreement.values():
        agg.update(counter)
    deltas_summary = {f'{arm}|{regime}|{level}': quantiles(values)
                      for (arm, regime, level), values in sorted(deltas.items())}
    random_exact = by_regime['random']['exact_fell'] / by_regime['random']['draws']
    random_noisy_all = by_regime['random']['noisy_fell'] / by_regime['random']['draws']
    ratios = {}
    for regime in REGIMES:
        if regime == 'random' or regime not in by_regime:
            continue
        counter = by_regime[regime]
        ratios[regime] = dict(
            exact=(counter['exact_fell'] / counter['draws']) / random_exact if random_exact else None,
            noisy_over_all_draws=((counter['noisy_fell'] / counter['draws']) / random_noisy_all
                                  if random_noisy_all else None))
    return dict(
        definition=('common draw set: K in {3, 5}, both arms have a full-panel baseline, the '
                    'exact arm completed; a DecoDiPhy failure is retained as noisy_no_answer'),
        overall=block(overall), by_regime={k: block(v) for k, v in sorted(by_regime.items())},
        by_level={k: block(v) for k, v in sorted(by_level.items())},
        by_regime_level=by_regime_level,
        structured_over_random_ratio=ratios,
        direction_agreement_overall=dict(agg),
        direction_agreement_by_regime={f'{r}|{l}': dict(v) for (r, l), v in sorted(agreement.items())},
        delta_quantiles=deltas_summary)


def nontrivial_monophyly():
    """Monophyly of removed taxonomic groups, excluding trivially monophyletic singletons."""
    ranks = s294.load_ranks()
    plan = s294.draw_plan(ranks)
    counts = collections.Counter()
    sizes = collections.Counter()
    cache = {}
    for noisy, _exact, meta in plan:
        if meta['regime'] != 'taxon':
            continue
        leaves = sorted(noisy['observed'])
        if meta['name'] not in cache:
            cache[meta['name']] = bipartitions(noisy['edges'], leaves)
        parts = cache[meta['name']]
        groups = s294.taxon_groups(ranks, noisy['leaf_labels'], leaves, meta['detail']['rank'])
        for name in meta['detail']['groups']:
            members = groups[name]
            sizes[len(members)] += 1
            if len(members) < 2:
                counts['singleton_trivially_monophyletic'] += 1
                continue
            counts['size_at_least_two'] += 1
            counts['size_at_least_two_monophyletic' if frozenset(members) in parts
                   else 'size_at_least_two_not_monophyletic'] += 1
    return dict(counts=dict(counts), group_size_histogram=dict(sorted(sizes.items())),
                fraction_monophyletic_among_groups_of_size_at_least_two=(
                    counts['size_at_least_two_monophyletic'] / counts['size_at_least_two']
                    if counts['size_at_least_two'] else None))


def write_matched():
    path = OUTDIR / 'stage294_matched_comparison.json'
    assert not path.exists(), 'Never overwrite a computed comparison'
    payload = dict(stage=294, matched=matched_comparison(),
                   taxon_group_monophyly=nontrivial_monophyly())
    path.write_text(json.dumps(payload, indent=1) + '\n')
    print(json.dumps(dict(overall=payload['matched']['overall'],
                          ratios=payload['matched']['structured_over_random_ratio'],
                          monophyly=payload['taxon_group_monophyly']['counts']), indent=1))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('command', choices=('metrics', 'summary', 'matched'))
    parser.add_argument('--workers', type=int, default=30)
    args = parser.parse_args()
    if args.command == 'matched':
        write_matched()
    elif args.command == 'metrics':
        run_metrics(args.workers)
    else:
        summarise()
