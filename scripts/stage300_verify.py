#!/usr/bin/env python3
"""Stage 300 verification: three independent routes to the same leaf equations.

ROUTE A (every evaluated run).  The frozen min form
    dist(r, q) = min(d(r,u) + xL, d(r,v) + (1-x)L)
and stage 296's interpolation form
    dist(r, q) = (1 - x) d(r,u) + x d(r,v)
must agree exactly in a tree: whichever endpoint the path to r leaves the edge
by, the other candidate exceeds it by twice the detour.  Any disagreement is a
bug in the geometry, so both are computed in exact rational arithmetic and
compared for equality, not closeness.

ROUTE B (whole cases, independent implementation).  The k source points are
INSERTED into the edge list as new vertices that split their edges, and the
distance from every observed leaf to every source point is obtained from
`scipy.sparse.csgraph.dijkstra` on the resulting weighted graph - a library
shortest-path routine that shares no code with this project.  The model
prediction is then reassembled from those distances and compared, as float64,
with the exact prediction.

ROUTE C (DecoDiPhy's own matrices).  For runs whose anchors survive the tool's
own pruning unchanged, the prediction is recomputed as
`D[:, anchors] @ p + (C[:, anchors] * l[anchors]) @ (p * x) + y` with D, C, l
built by `decodiphy.matrix.get_input_matrices` from the very tree file the
bridge writes - i.e. using the published tool's own geometry, not ours.
"""
import argparse
import json
import sys
import time
from fractions import Fraction as F
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
sys.path.insert(0, str(ROOT / 'data' / 'stage289_decodiphy' / 'repo'))
import stage300_common as s300                   # noqa: E402
import stage300_evaluate as s300eval             # noqa: E402


def route_b(edges, leaves, sources, pendant):
    """Predictions via scipy Dijkstra on the edge list with source points inserted."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import dijkstra
    vertices = sorted({v for e in edges for v in e[:2]})
    index = {v: i for i, v in enumerate(vertices)}
    rows, cols, data = [], [], []
    next_index = len(vertices)
    split = {}
    for j, s in enumerate(sources):
        u, v, length = edges[int(s['edge'])]
        x = float(s['x'])
        q = next_index
        next_index += 1
        split.setdefault(int(s['edge']), []).append((q, x))
    for e, (u, v, length) in enumerate(edges):
        if e in split:
            continue
        rows += [index[u], index[v]]
        cols += [index[v], index[u]]
        data += [float(length), float(length)]
    for e, points in split.items():
        u, v, length = edges[e]
        points = sorted(points, key=lambda z: z[1])
        previous, previous_x = index[u], 0.0
        for q, x in points:
            rows += [previous, q]
            cols += [q, previous]
            data += [(x - previous_x) * float(length)] * 2
            previous, previous_x = q, x
        rows += [previous, index[v]]
        cols += [index[v], previous]
        data += [(1.0 - previous_x) * float(length)] * 2
    n = next_index
    graph = coo_matrix((data, (rows, cols)), shape=(n, n)).tocsr()
    source_nodes = [q for e in split for (q, _) in split[e]]
    dist = dijkstra(graph, directed=False, indices=source_nodes)
    node_of = {}
    counter = 0
    for e in split:
        for (q, x) in split[e]:
            node_of[q] = counter
            counter += 1
    pred = np.full(len(leaves), float(pendant), dtype=np.float64)
    for j, s in enumerate(sources):
        # the j-th source's inserted node, matched by (edge, x)
        candidates = [q for (q, x) in split[int(s['edge'])]
                      if abs(x - float(s['x'])) < 1e-15]
        q = candidates[0]
        rowvec = dist[node_of[q]]
        pred += float(s['p']) * np.array([rowvec[index[r]] for r in leaves])
    return pred


def route_c(case, run, edges, leaves):
    """Predictions from DecoDiPhy's own matrices; None when its pruning moved an anchor."""
    import tempfile
    from treeswift import read_tree_newick
    from decodiphy.tree_utils import preprocess_input
    from decodiphy.matrix import get_input_matrices
    from stage289_decodiphy_bridge import write_inputs
    bank = dict(edges=edges, observed={r: F(case['observed'][str(r)]) for r in leaves},
                generating_sources=case.get('generating_sources') or [])
    with tempfile.TemporaryDirectory(prefix='stage300_verify_') as work:
        mapping = write_inputs(bank, work)
        text = (Path(work) / 'tree.nwk').read_text().strip()
        distances = {}
        for line in (Path(work) / 'distances.txt').read_text().strip().split('\n'):
            a, b = line.split()
            distances[a] = float(b)
        original = read_tree_newick(text)
        tree_obj = read_tree_newick(text)
        preprocess_input(tree_obj, distances)
        d, l, C, D, index_to_node, node_to_index, index_to_leaf, _, _ = \
            get_input_matrices(tree_obj, distances)
    labels = [s['label'] for s in run['sources']]
    if any(label not in node_to_index for label in labels):
        return None, 'anchor_absent_from_pruned_tree'
    pruned_nodes = tree_obj.label_to_node(selection='all')
    original_nodes = original.label_to_node(selection='all')
    for label in labels:
        if pruned_nodes[label].parent.label != original_nodes[label].parent.label:
            return None, 'anchor_reparented_by_pruning'
    anchors = [node_to_index[label] for label in labels]
    p = np.array([float(s['p']) for s in run['sources']])
    x = np.array([float(s['x_decodiphy']) for s in run['sources']])
    y = float(run['pendant_total'] or 0.0)
    DA = D[:, anchors]
    CAl = C[:, anchors] * l[anchors].reshape(-1)
    pred = DA @ p + CAl @ (p * x) + y
    leaf_order = {label: i for i, label in enumerate(index_to_leaf)}
    leafset = {v for v in {a for e in edges for a in e[:2]}
               if sum(1 for e in edges if v in e[:2]) == 1}
    want = []
    for r in leaves:
        label = ('L%d' if r in leafset else 'N%d') % r
        want.append(pred[leaf_order[label]])
    return np.array(want), 'ok'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--collections', default='krepp_ani,krepp_placement,travelers')
    ap.add_argument('--route-a-limit', type=int, default=0,
                    help='0 = every evaluated run')
    ap.add_argument('--route-b-cases', type=int, default=4)
    ap.add_argument('--route-c-cases', type=int, default=4)
    ap.add_argument('--out', default='stage300_verification.json')
    args = ap.parse_args()

    report = dict(stage=300, step='verification', routes={}, started=time.time())
    for collection in args.collections.split(','):
        entries = {e['name']: e for e in s300.case_entries(collection)}
        runs = {}
        for run in s300eval.runs_for(collection):
            if run['status'] == 'complete' and run.get('sources'):
                runs.setdefault(run['name'], []).append(run)
        cell = dict(route_a=dict(compared=0, agree=0, disagree=[]),
                    route_b=dict(cases=[], max_abs_diff=0.0),
                    route_c=dict(cases=[], skipped={}, max_abs_diff=0.0))
        names = sorted(runs)
        for position, name in enumerate(names):
            case = s300.load_case(entries[name])
            edges, leaves, dvals = s300.case_arrays(case)
            for run in runs[name]:
                sources, pendant = s300eval.configuration(run)
                a1 = s300.predictions_exact(edges, leaves, sources, pendant)
                a2 = s300.predictions_interpolation(edges, leaves, sources, pendant)
                cell['route_a']['compared'] += 1
                if a1 == a2:
                    cell['route_a']['agree'] += 1
                else:
                    cell['route_a']['disagree'].append(
                        dict(name=name, arm=run['arm'], source_stage=run['source_stage']))
                if args.route_a_limit and cell['route_a']['compared'] >= args.route_a_limit:
                    break
            if position < args.route_b_cases:
                run = runs[name][0]
                sources, pendant = s300eval.configuration(run)
                exact = np.array([float(z) for z in s300.predictions_exact(
                    edges, leaves, sources, pendant)])
                viab = route_b(edges, leaves, run['sources'], float(run['pendant_total'] or 0.0))
                diff = float(np.max(np.abs(exact - viab)))
                cell['route_b']['cases'].append(dict(name=name, arm=run['arm'], k=run['k'],
                                                     max_abs_diff=diff,
                                                     scale=float(np.max(np.abs(exact)))))
                cell['route_b']['max_abs_diff'] = max(cell['route_b']['max_abs_diff'], diff)
            if position < args.route_c_cases:
                run = runs[name][0]
                sources, pendant = s300eval.configuration(run)
                exact = np.array([float(z) for z in s300.predictions_exact(
                    edges, leaves, sources, pendant)])
                try:
                    viac, why = route_c(case, run, edges, leaves)
                except BaseException as exc:
                    viac, why = None, repr(exc)[:200]
                if viac is None:
                    cell['route_c']['skipped'][why] = cell['route_c']['skipped'].get(why, 0) + 1
                else:
                    diff = float(np.max(np.abs(exact - viac)))
                    cell['route_c']['cases'].append(dict(name=name, arm=run['arm'],
                                                         k=run['k'], max_abs_diff=diff))
                    cell['route_c']['max_abs_diff'] = max(cell['route_c']['max_abs_diff'], diff)
            if args.route_a_limit and cell['route_a']['compared'] >= args.route_a_limit:
                break
        report['routes'][collection] = cell
        print(collection, 'route A %d/%d agree' % (cell['route_a']['agree'],
                                                   cell['route_a']['compared']),
              'route B max diff %.3e' % cell['route_b']['max_abs_diff'],
              'route C cases %d max diff %.3e skipped %s'
              % (len(cell['route_c']['cases']), cell['route_c']['max_abs_diff'],
                 cell['route_c']['skipped']), flush=True)
    report['elapsed_seconds'] = time.time() - report.pop('started')
    (s300.OUT / args.out).write_text(json.dumps(s300.serial(report), indent=1))
    print('wrote', s300.OUT / args.out)


if __name__ == '__main__':
    main()
