#!/usr/bin/env python3
"""Stage 300 diagnostic: where does an impossible reported loss come from?

One krepp ANI case (k10k_tree_15_k3_s1_exp) is reported by DecoDiPhy with
loss 3.6e-20 at k = 2 - a near-exact fit to 51 real ANI distances, which the
exact model says is impossible (its Theta_2 is strictly positive).  Either the
evaluation in this stage is wrong, or the number the tool prints does not belong
to the support it prints.  This script decides between the two by re-solving the
tool's own quadratic programme, with its own routine, at the printed support and
at every single-anchor neighbour of it, and reporting the attained objectives.
"""
import json
import sys
import tempfile
from fractions import Fraction as F
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
sys.path.insert(0, str(ROOT / 'data' / 'stage289_decodiphy' / 'repo'))
import stage300_common as s300                   # noqa: E402
import stage300_evaluate as s300eval             # noqa: E402

NAME = sys.argv[1] if len(sys.argv) > 1 else 'k10k_tree_15_k3_s1_exp_krepp_estimated'
ARM = sys.argv[2] if len(sys.argv) > 2 else 'dd_default'


def main():
    from treeswift import read_tree_newick
    from decodiphy.tree_utils import preprocess_input
    from decodiphy.matrix import get_input_matrices
    from decodiphy.optimize import solve_osqp_once
    from stage289_decodiphy_bridge import write_inputs
    entry = {e['name']: e for e in s300.case_entries('krepp_ani')}[NAME]
    case = s300.load_case(entry)
    edges, leaves, dvals = s300.case_arrays(case)
    runs = [r for r in s300eval.runs_for('krepp_ani')
            if r['name'] == NAME and r['arm'] == ARM and r['status'] == 'complete']
    run = runs[0]
    bank = dict(edges=edges, observed={r: v for r, v in zip(leaves, dvals)},
                generating_sources=case.get('generating_sources') or [])
    with tempfile.TemporaryDirectory(prefix='stage300_losscase_') as work:
        write_inputs(bank, work)
        text = (Path(work) / 'tree.nwk').read_text().strip()
        distances = {}
        for line in (Path(work) / 'distances.txt').read_text().strip().split('\n'):
            a, b = line.split()
            distances[a] = float(b)
        tree_obj = read_tree_newick(text)
        preprocess_input(tree_obj, distances)
        d, l, C, D, index_to_node, node_to_index, index_to_leaf, _, _ = \
            get_input_matrices(tree_obj, distances)
    labels = [s['label'] for s in run['sources']]
    anchors = np.array([node_to_index[a] for a in labels], dtype=int)
    k = len(anchors)
    at_support = solve_osqp_once(d, l.reshape(-1), C.T, D.T, anchors, index_to_node, k)[0]
    neighbours = []
    n_nodes = len(index_to_node)
    for position in range(k):
        for candidate in range(n_nodes):
            if candidate in anchors.tolist():
                continue
            trial = anchors.copy()
            trial[position] = candidate
            obj = solve_osqp_once(d, l.reshape(-1), C.T, D.T, trial, index_to_node, k)[0]
            neighbours.append(dict(position=position, candidate=index_to_node[candidate],
                                   objective=float(obj)))
    neighbours.sort(key=lambda z: z['objective'])
    report = dict(stage=300, step='loss_case', name=NAME, arm=ARM,
                  reported_loss=run['reported_loss'], reported_anchors=labels,
                  reported_k=run['k'],
                  objective_at_reported_support=float(at_support),
                  best_neighbour=neighbours[0], worst_neighbour=neighbours[-1],
                  n_neighbours=len(neighbours),
                  min_neighbour_objective=neighbours[0]['objective'],
                  ladder=run.get('ladder'),
                  conclusion=('the printed loss is attained at no support within one '
                              'anchor swap of the printed support'
                              if min(neighbours[0]['objective'], float(at_support))
                              > 100 * (run['reported_loss'] or 0) + 1e-12
                              else 'the printed loss is attained at a neighbouring support'))
    (s300.OUT / ('stage300_losscase_%s.json' % NAME)).write_text(
        json.dumps(s300.serial(report), indent=1))
    print(json.dumps({k2: v for k2, v in report.items() if k2 != 'ladder'}, indent=1))
    print('ladder:', json.dumps(run.get('ladder'))[:600])


if __name__ == '__main__':
    main()
