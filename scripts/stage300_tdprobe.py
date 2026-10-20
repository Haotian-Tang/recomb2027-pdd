#!/usr/bin/env python3
"""Stage 300 probe: where does DecoDiPhy's time go on a 10,575-leaf backbone?

The travellers'-diarrhoea samples are two orders of magnitude larger than the
krepp cases, and the tool's wall time there is dominated by its input stage, not
its search.  This times the stages separately on one sample so the cost can be
reported rather than guessed: Newick parsing, `preprocess_input`, treeswift's
leaf-to-leaf `distance_matrix`, and the whole of `get_input_matrices`.
"""
import json
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
sys.path.insert(0, str(ROOT / 'data' / 'stage289_decodiphy' / 'repo'))
import stage300_common as s300                   # noqa: E402

NAME = sys.argv[1] if len(sys.argv) > 1 else 'td_sample_101_published'


def main():
    from treeswift import read_tree_newick
    from decodiphy.tree_utils import preprocess_input
    from decodiphy.matrix import get_input_matrices
    from stage289_decodiphy_bridge import write_inputs
    entry = {e['name']: e for e in s300.case_entries('travelers')}[NAME]
    case = s300.load_case(entry)
    edges, leaves, dvals = s300.case_arrays(case)
    timings = {}
    t0 = time.perf_counter()
    bank = dict(edges=edges, observed={r: v for r, v in zip(leaves, dvals)},
                generating_sources=[])
    work = Path(tempfile.mkdtemp(prefix='stage300_tdprobe_'))
    write_inputs(bank, work)
    timings['write_inputs'] = time.perf_counter() - t0
    text = (work / 'tree.nwk').read_text().strip()
    distances = {}
    for line in (work / 'distances.txt').read_text().strip().split('\n'):
        a, b = line.split()
        distances[a] = float(b)
    t0 = time.perf_counter()
    tree_obj = read_tree_newick(text)
    timings['read_tree_newick'] = time.perf_counter() - t0
    t0 = time.perf_counter()
    pruned = preprocess_input(tree_obj, distances)
    timings['preprocess_input'] = time.perf_counter() - t0
    n_after = sum(1 for _ in tree_obj.traverse_leaves())
    print(json.dumps(dict(stage='partial', pruned=bool(pruned), n_leaves_after=n_after,
                          timings=timings)), flush=True)
    t0 = time.perf_counter()
    dm = tree_obj.distance_matrix(leaf_labels=True)
    timings['distance_matrix'] = time.perf_counter() - t0
    timings['distance_matrix_entries'] = sum(len(v) for v in dm.values())
    del dm
    print(json.dumps(dict(stage='distance_matrix', timings=timings)), flush=True)
    t0 = time.perf_counter()
    out = get_input_matrices(tree_obj, distances)
    timings['get_input_matrices_total'] = time.perf_counter() - t0
    timings['n_nodes'] = len(out[4])
    timings['n_leaves'] = len(out[6])
    report = dict(stage=300, step='td_probe', name=NAME, pruned=bool(pruned),
                  n_leaves_after_prune=n_after, timings=timings)
    (s300.OUT / ('stage300_tdprobe_%s.json' % NAME)).write_text(json.dumps(report, indent=1))
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    main()
