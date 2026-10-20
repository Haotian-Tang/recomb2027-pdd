#!/usr/bin/env python3
"""Stage 289 analysis: what DecoDiPhy and the exact solver return on the same cases.

Per case and method the run stream stores the returned placements; this script
derives, against the bank's published ground truth:
  * the number of sources returned,
  * whether the returned support is the generating support,
  * whether the returned support admits an exact explanation at all
    (exact rational decision, independent of either solver),
  * the Wasserstein-1 distance between the returned mixture and the generating
    mixture under the tree metric, in branch-length units,
  * runtime, as median child CPU seconds over the repetitions.
Timeouts are retained and never imputed.
"""
import json
import sys
from collections import Counter, defaultdict
from fractions import Fraction as F
from pathlib import Path
from statistics import median

import numpy as np
from scipy.optimize import linprog

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from stage289_decodiphy_bridge import load_case
from stage289_support_feasibility import feasible_support

PROTOCOL = ROOT / 'results/stage289/stage289_headtohead_protocol.json'
STREAM = ROOT / 'results/stage289/stage289_headtohead_runs.jsonl'
PERCASE = ROOT / 'results/stage289/stage289_headtohead_percase.json'
SUMMARY = ROOT / 'results/stage289/stage289_headtohead_summary.json'
DD = ('dd_default', 'dd_true_k', 'dd_resolution')


def float_adjacency(edges):
    adjacency = {}
    for u, v, length in edges:
        value = float(length)
        adjacency.setdefault(u, []).append((v, value))
        adjacency.setdefault(v, []).append((u, value))
    return adjacency


def distances_from(adjacency, start):
    out = {start: 0.0}
    stack = [start]
    while stack:
        v = stack.pop()
        for u, length in adjacency[v]:
            if u not in out:
                out[u] = out[v] + length
                stack.append(u)
    return out


def point_distances(edges, adjacency, table, left, right):
    """Tree distance matrix between two lists of {edge, x} placement points."""
    matrix = []
    for a in left:
        ua, va, la = edges[a['edge']]
        la, xa = float(la), float(a['x'])
        row = []
        for b in right:
            ub, vb, lb = edges[b['edge']]
            lb, xb = float(lb), float(b['x'])
            if a['edge'] == b['edge']:
                row.append(abs(xa - xb) * la)
                continue
            row.append(min(table[p][q] + sa + sb
                           for p, sa in ((ua, xa * la), (va, (1 - xa) * la))
                           for q, sb in ((ub, xb * lb), (vb, (1 - xb) * lb))))
        matrix.append(row)
    return matrix


def wasserstein(cost, weights_a, weights_b):
    """Balanced optimal transport between two discrete mixtures."""
    n, m = len(weights_a), len(weights_b)
    cost = np.asarray(cost, dtype=float).reshape(-1)
    rows = []
    for i in range(n):
        row = np.zeros((n, m))
        row[i, :] = 1
        rows.append(row.reshape(-1))
    for j in range(m):
        row = np.zeros((n, m))
        row[:, j] = 1
        rows.append(row.reshape(-1))
    equalities = np.array(rows)
    target = np.concatenate([np.asarray(weights_a, float), np.asarray(weights_b, float)])
    answer = linprog(cost, A_eq=equalities, b_eq=target, bounds=(0, None), method='highs')
    return float(answer.fun) if answer.success else None


def analyse_case(entry, runs):
    case = load_case(ROOT / entry['path'])
    edges = case['edges']
    adjacency = float_adjacency(edges)
    truth = [dict(edge=s['edge'], p=float(F(s['p'])), x=float(F(s['x'])))
             for s in case['generating_sources']]
    needed = {v for s in truth for v in edges[s['edge']][:2]}
    for run in runs.values():
        for value in run:
            for s in value.get('sources', []) or []:
                needed.update(edges[s['edge']][:2])
            witness = value.get('witness')
            if witness:
                for s in witness['sources']:
                    needed.update(edges[s['edge']][:2])
    table = {v: distances_from(adjacency, v) for v in needed}
    true_support = {s['edge'] for s in truth}
    out = dict(name=entry['name'], family=entry['family'], leaves=entry['reference_set_size'],
               n_edges=len(edges), k_true=entry['K_generating'], cap=entry['cap'],
               exactly_explained=entry['exactly_explained'], methods={})
    for method, values in runs.items():
        statuses = [v['status'] for v in values]
        record = dict(status=statuses[0], statuses=statuses,
                      cpu_seconds=median(v['child_cpu_seconds'] for v in values),
                      wall_seconds=median(v['process_wall_seconds'] for v in values),
                      max_rss_mib=max(v['maximum_rss_kib'] for v in values) / 1024)
        value = values[0]
        if method == 'exact_stage43_partial':
            record.update(feasible=value.get('feasible'), minimum=value.get('minimum'),
                          reason=value.get('reason'))
            witness = value.get('witness')
            if witness:
                sources = [dict(edge=s['edge'], p=float(F(s['p'])), x=float(F(s['x'])))
                           for s in witness['sources']]
                record['support'] = sorted(s['edge'] for s in sources)
                record['support_exact_match'] = set(record['support']) == true_support
                record['n_true_edges_recovered'] = len(set(record['support']) & true_support)
                cost = point_distances(edges, adjacency, table, sources, truth)
                record['w1_to_truth'] = wasserstein(cost, [s['p'] for s in sources],
                                                    [s['p'] for s in truth])
                record['weighted_placement_error'] = sum(
                    t['p'] * min(cost[i][j] for i in range(len(sources)))
                    for j, t in enumerate(truth))
                record['pendant'] = float(F(witness['pendant']))
            record['answers_agree_across_repetitions'] = len(
                {(v.get('feasible'), v.get('minimum'), v['status']) for v in values}) == 1
        elif value['status'] == 'complete':
            sources = value['sources']
            support = [s['edge'] for s in sources]
            record.update(k=value['k'], loss=value['loss'], stop_reason=value['stop_reason'],
                          support=sorted(support),
                          repeated_edge=len(set(support)) != len(support),
                          support_exact_match=set(support) == true_support,
                          n_true_edges_recovered=len(set(support) & true_support),
                          pendant=value['pendant_total'],
                          boundary_placements=sum(s['x'] <= 1e-9 or s['x'] >= 1 - 1e-9 for s in sources))
            cost = point_distances(edges, adjacency, table, sources, truth)
            record['w1_to_truth'] = wasserstein(cost, [s['p'] for s in sources],
                                                [s['p'] for s in truth])
            record['weighted_placement_error'] = sum(
                t['p'] * min(cost[i][j] for i in range(len(sources))) for j, t in enumerate(truth))
            verdict = feasible_support(edges, case['observed'], sorted(set(support)))
            record['exact_feasibility'] = verdict['verdict']
            record['exact_feasibility_reason'] = verdict.get('reason')
            record['identifiable'] = verdict.get('identifiable')
            record['free_parameters'] = verdict.get('free_parameters')
            record['k_values_agree_across_repetitions'] = len({v.get('k') for v in values}) == 1
        out['methods'][method] = record
    raw_pendant = case.get('generating_pendant')
    truth_pendant = float(F(raw_pendant)) if raw_pendant is not None else None
    out['generating_pendant'] = truth_pendant
    for method, record in out['methods'].items():
        if record.get('pendant') is not None and truth_pendant:
            record['pendant_relative_error'] = abs(record['pendant'] - truth_pendant) / truth_pendant
    return out


def _analyse(payload):
    return analyse_case(*payload)


def bucket(leaves):
    return '<=50' if leaves <= 50 else '51-150' if leaves <= 150 else '151-500'


def summarise(percase):
    summary = {}
    for family in sorted({c['family'] for c in percase}):
        subset = [c for c in percase if c['family'] == family]
        block = dict(cases=len(subset), by_method={}, by_k={}, by_reference_size={})
        for method in ('exact_stage43_partial',) + DD:
            runs = [c['methods'][method] for c in subset if method in c['methods']]
            complete = [r for r in runs if r['status'] == 'complete']
            entry = dict(runs=len(runs), status_counts=dict(Counter(r['status'] for r in runs)),
                         median_cpu_seconds_completed=median([r['cpu_seconds'] for r in complete]) if complete else None,
                         total_cpu_seconds_all_runs=sum(r['cpu_seconds'] for r in runs))
            if method == 'exact_stage43_partial':
                entry.update(
                    feasible=sum(bool(r.get('feasible')) for r in complete),
                    infeasible=sum(r.get('feasible') is False for r in complete),
                    reasons=dict(Counter(r.get('reason') for r in complete if r.get('feasible') is False)),
                    minimum_equals_k_true=sum(r.get('minimum') == c['k_true']
                                              for r, c in zip(runs, subset) if r['status'] == 'complete'),
                    support_exact_match=sum(bool(r.get('support_exact_match')) for r in complete),
                    repetition_disagreements=sum(not r.get('answers_agree_across_repetitions', True) for r in runs))
            else:
                entry.update(
                    k_equals_k_true=sum(r.get('k') == c['k_true'] for r, c in zip(runs, subset)
                                        if r['status'] == 'complete'),
                    k_below_k_true=sum(r.get('k', 0) < c['k_true'] for r, c in zip(runs, subset)
                                       if r['status'] == 'complete'),
                    k_above_k_true=sum(r.get('k', 0) > c['k_true'] for r, c in zip(runs, subset)
                                       if r['status'] == 'complete'),
                    support_exact_match=sum(bool(r.get('support_exact_match')) for r in complete),
                    exact_feasibility=dict(Counter(r.get('exact_feasibility') for r in complete)),
                    stop_reasons=dict(Counter(r.get('stop_reason') for r in complete)),
                    boundary_placement_cases=sum(bool(r.get('boundary_placements')) for r in complete),
                    repeated_edge_cases=sum(bool(r.get('repeated_edge')) for r in complete),
                    median_w1_to_truth=median([r['w1_to_truth'] for r in complete
                                               if r.get('w1_to_truth') is not None]) if complete else None,
                    median_loss=median([r['loss'] for r in complete]) if complete else None,
                    recovered_fraction_of_true_edges=(
                        sum(r.get('n_true_edges_recovered', 0) for r in complete)
                        / sum(c['k_true'] for c, r in zip(subset, runs) if r['status'] == 'complete')
                        if complete else None))
            block['by_method'][method] = entry
        for key, selector in (('by_k', lambda c: c['k_true']),
                              ('by_reference_size', lambda c: bucket(c['leaves']))):
            for value in sorted({selector(c) for c in subset}, key=str):
                group = [c for c in subset if selector(c) == value]
                block[key][str(value)] = dict(
                    cases=len(group),
                    exact_complete=sum(g['methods']['exact_stage43_partial']['status'] == 'complete' for g in group),
                    exact_median_cpu=median([g['methods']['exact_stage43_partial']['cpu_seconds'] for g in group]),
                    **{f'{m}_support_exact_match': sum(bool(g['methods'][m].get('support_exact_match')) for g in group)
                       for m in DD},
                    **{f'{m}_exact_feasible': sum(g['methods'][m].get('exact_feasibility') == 'feasible' for g in group)
                       for m in DD},
                    **{f'{m}_median_cpu': median([g['methods'][m]['cpu_seconds'] for g in group]) for m in DD},
                    **{f'{m}_median_w1': median([g['methods'][m]['w1_to_truth'] for g in group
                                                 if g['methods'][m].get('w1_to_truth') is not None])
                       if any(g['methods'][m].get('w1_to_truth') is not None for g in group) else None
                       for m in DD})
        summary[family] = block
    return summary


def main():
    protocol = json.loads(PROTOCOL.read_text())
    index = {entry['name']: entry for entry in protocol['case_index']}
    runs = defaultdict(lambda: defaultdict(list))
    with STREAM.open() as stream:
        for line in stream:
            record = json.loads(line)
            runs[record['name']][record['method']].append(record['result'])
    expected = set(('exact_stage43_partial',) + DD)
    percase, incomplete, pending = [], [], []
    for name, by_method in runs.items():
        if set(by_method) != expected or any(len(v) != protocol['repetitions']
                                             for v in by_method.values()):
            incomplete.append(name)
            continue
        pending.append((index[name], dict(by_method)))
    import multiprocessing as mp
    with mp.get_context('fork').Pool(12) as pool:
        percase = pool.map(_analyse, pending, chunksize=4)
    percase.sort(key=lambda c: (c['family'], c['name']))
    PERCASE.write_text(json.dumps(percase, indent=1) + '\n')
    summary = dict(stage=289, cases_analysed=len(percase),
                   cases_with_incomplete_run_records=len(incomplete),
                   incomplete_case_names=sorted(incomplete)[:20],
                   protocol_sha256=__import__('hashlib').sha256(PROTOCOL.read_bytes()).hexdigest(),
                   runs_sha256=__import__('hashlib').sha256(STREAM.read_bytes()).hexdigest(),
                   families=summarise(percase))
    SUMMARY.write_text(json.dumps(summary, indent=1) + '\n')
    print(json.dumps(summary['families'], indent=1)[:4000])


if __name__ == '__main__':
    main()
