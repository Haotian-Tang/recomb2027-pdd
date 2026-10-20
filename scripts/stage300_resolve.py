#!/usr/bin/env python3
"""Stage 300 fairness arm: re-solve DecoDiPhy's OWN quadratic programme at the
support it returned, with its own optimiser, and evaluate THAT configuration.

Why this exists.  For many runs the recomputed SSR of the configuration written
to `all_rounds.json` (and from there to the jplace) is far above the loss the
tool reports for it.  `decodiphy/search.py:hill_climbing` resets `min_obj` at
the top of every sweep but never resets `min_p / min_x / min_y`, so when the
final sweep improves nothing the parameters written out are those of the last
improving candidate of an earlier sweep while the anchors and the loss are those
of the converged support.  The returned triple can therefore be internally
inconsistent.

To keep the comparison fair, the published tool is given the benefit of its own
optimiser: `decodiphy.optimize.solve_osqp_once` is called on the tool's own
matrices at the tool's own returned anchors, and the resulting configuration is
evaluated with the same exact code path as everything else.  Both numbers are
reported - the configuration a user reads, and the configuration its optimiser
actually attains on that support.
"""
import argparse
import json
import sys
import tempfile
import time
from fractions import Fraction as F
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
sys.path.insert(0, str(ROOT / 'data' / 'stage289_decodiphy' / 'repo'))
import stage300_common as s300                   # noqa: E402
import stage300_evaluate as s300eval             # noqa: E402
import stage300_analyze as s300an                # noqa: E402


def resolve_case(case, edges, leaves, dvals, run):
    from treeswift import read_tree_newick
    from decodiphy.tree_utils import preprocess_input, compute_individual_ys
    from decodiphy.matrix import get_input_matrices
    from decodiphy.optimize import solve_osqp_once
    from stage289_decodiphy_bridge import write_inputs, to_bank_sources
    bank = dict(edges=edges, observed={r: v for r, v in zip(leaves, dvals)},
                generating_sources=case.get('generating_sources') or [])
    with tempfile.TemporaryDirectory(prefix='stage300_resolve_') as work:
        mapping = write_inputs(bank, work)
        text = (Path(work) / 'tree.nwk').read_text().strip()
        distances = {}
        for line in (Path(work) / 'distances.txt').read_text().strip().split('\n'):
            a, b = line.split()
            distances[a] = float(b)
        tree_obj = read_tree_newick(text)
        preprocess_input(tree_obj, distances)
        d, l, C, D, index_to_node, node_to_index, index_to_leaf, branch_lengths, \
            terminal_lengths = get_input_matrices(tree_obj, distances)
    labels = [s['label'] for s in run['sources']]
    if any(label not in node_to_index for label in labels):
        return dict(status='anchor_absent_from_pruned_tree')
    anchors = np.array([node_to_index[label] for label in labels], dtype=int)
    if len(set(anchors.tolist())) != len(anchors):
        return dict(status='repeated_anchor', k=len(anchors))
    obj, p, x, y, _ = solve_osqp_once(d, l.reshape(-1), C.T, D.T, anchors,
                                      index_to_node, len(anchors))
    if p is None:
        return dict(status='osqp_failed')
    ys = compute_individual_ys([index_to_node[i] for i in anchors], x, y,
                               terminal_lengths, branch_lengths)
    srcs = to_bank_sources(bank, [index_to_node[i] for i in anchors], p, x, mapping)
    unmapped = [s for s in srcs if s is None]
    sources = [dict(edge=int(s['edge']), p=F(float(s['p'])), x=F(float(s['x'])))
               for s in srcs]
    metrics = s300.fit_metrics(edges, leaves, dvals, sources, F(float(sum(ys))))
    return dict(status='ok', k=len(anchors), objective=float(obj),
                ssr=metrics['ssr'], linf=metrics['linf'],
                pendant_total=float(sum(ys)), sum_p=float(np.sum(p)),
                min_p=float(np.min(p)), x=[float(v) for v in x],
                p=[float(v) for v in p],
                support=sorted(int(s['edge']) for s in srcs))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--collection', required=True)
    ap.add_argument('--arms', default='dd_default,dd_true_k')
    ap.add_argument('--limit', type=int, default=0)
    ap.add_argument('--out', default='')
    args = ap.parse_args()
    arms = tuple(args.arms.split(','))
    entries = {e['name']: e for e in s300.case_entries(args.collection)}
    by_case = {}
    for run in s300eval.runs_for(args.collection):
        if run['status'] == 'complete' and run.get('sources'):
            by_case.setdefault((run['name'], run['arm']), []).append(run)
    theta = s300.theta1_table()
    out_path = s300.OUT / (args.out or ('stage300_resolved_%s.jsonl' % args.collection))
    handle = open(out_path, 'w')
    names = sorted({n for (n, a) in by_case})
    if args.limit:
        names = names[:args.limit]
    t0 = time.perf_counter()
    for position, name in enumerate(names, 1):
        case = s300.load_case(entries[name])
        edges, leaves, dvals = s300.case_arrays(case)
        for arm in arms:
            runs = by_case.get((name, arm))
            if not runs:
                continue
            evaluated = [dict(r, evaluated=True) for r in runs]
            best = sorted(evaluated, key=lambda r: (
                s300an.STAGE_ORDER.get(r['source_stage'], 9), r.get('source_line') or 0))[0]
            try:
                answer = resolve_case(case, edges, leaves, dvals, best)
            except BaseException as exc:
                answer = dict(status='error', exception=repr(exc)[:300])
            t1 = theta.get(name)
            row = dict(kind='resolved', stage=300, collection=args.collection, name=name,
                       arm=arm, source_stage=best['source_stage'],
                       source_line=best['source_line'],
                       reported_loss=best.get('reported_loss'),
                       returned_k=best.get('k'),
                       theta1=(t1['theta1'] if t1 else None), resolved=answer)
            if answer.get('status') == 'ok' and t1:
                row['linf_ge_theta1'] = bool(F(answer['linf']) >= t1['theta1'])
                row['linf_over_theta1'] = float(F(answer['linf']) / t1['theta1'])
            handle.write(json.dumps(s300.serial(row)) + '\n')
        if position % 25 == 0:
            print('[%s] %s %d/%d %.0f s' % (time.strftime('%H:%M:%S'), args.collection,
                                            position, len(names),
                                            time.perf_counter() - t0), flush=True)
            handle.flush()
    handle.close()
    print('wrote', out_path, flush=True)


if __name__ == '__main__':
    main()
