#!/usr/bin/env python3
"""Stage 300, the fair-panel arm: which references did DecoDiPhy actually use?

DecoDiPhy 1.5.0 removes reference leaves that carry no distance, and does so
while iterating over the leaves it is removing (`decodiphy/tree_utils.py`:
`preprocess_input` loops over `tree_obj.traverse_leaves()` and calls
`parent.remove_child(l)` / `parent.contract()` inside the loop).  On the two
estimated-distance collections, where the estimator leaves a few references
without a value, the traversal is corrupted and the cascade deletes references
that DO carry a distance.  This script measures the surviving panel for every
case, and - because a comparison on data the tool never saw would not be fair -
recomputes BOTH sides of the comparison on the panel it did see:

    Theta_1 restricted to the surviving references (stage 298's exact k = 1
    routine, imported unmodified), and
    R_dd restricted to the same references.

Nothing here changes the full-panel numbers; it adds the restricted ones so both
can be reported.
"""
import argparse
import json
import sys
import tempfile
import time
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
sys.path.insert(0, str(ROOT / 'data' / 'stage289_decodiphy' / 'repo'))
import stage300_common as s300                   # noqa: E402
import stage300_evaluate as s300eval             # noqa: E402
import stage300_analyze as s300an                # noqa: E402
import stage298_common as s298                   # noqa: E402


def surviving_leaves(edges, leaves, dvals):
    """Labels and bank ids of the references DecoDiPhy keeps after preprocess_input."""
    from treeswift import read_tree_newick
    from decodiphy.tree_utils import preprocess_input
    from stage289_decodiphy_bridge import write_inputs
    bank = dict(edges=edges, observed={r: v for r, v in zip(leaves, dvals)},
                generating_sources=[])
    with tempfile.TemporaryDirectory(prefix='stage300_prune_') as work:
        write_inputs(bank, work)
        text = (Path(work) / 'tree.nwk').read_text().strip()
        distances = {}
        for line in (Path(work) / 'distances.txt').read_text().strip().split('\n'):
            a, b = line.split()
            distances[a] = float(b)
    tree_before = read_tree_newick(text)
    n_before = sum(1 for _ in tree_before.traverse_leaves())
    tree_obj = read_tree_newick(text)
    pruned = preprocess_input(tree_obj, distances)
    kept = [l.label for l in tree_obj.traverse_leaves()]
    kept_ids = sorted(int(label[1:]) for label in kept if label.startswith('L'))
    return dict(pruned_flag=bool(pruned), n_tree_leaves_before=n_before,
                n_kept=len(kept), kept_labels_non_leaf=[l for l in kept
                                                        if not l.startswith('L')],
                kept_ids=kept_ids, n_with_distance=len(distances))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--collection', required=True)
    ap.add_argument('--arms', default='dd_default,dd_true_k')
    ap.add_argument('--limit', type=int, default=0)
    args = ap.parse_args()
    entries = {e['name']: e for e in s300.case_entries(args.collection)}
    by_case = {}
    for run in s300eval.runs_for(args.collection):
        if run['status'] == 'complete' and run.get('sources'):
            by_case.setdefault((run['name'], run['arm']), []).append(run)
    theta_full = s300.theta1_table()
    path = s300.OUT / ('stage300_prune_%s.jsonl' % args.collection)
    handle = open(path, 'w')
    names = sorted(entries)
    if args.limit:
        names = names[:args.limit]
    t0 = time.perf_counter()
    for position, name in enumerate(names, 1):
        case = s300.load_case(entries[name])
        edges, leaves, dvals = s300.case_arrays(case)
        survive = surviving_leaves(edges, leaves, dvals)
        kept = [r for r in leaves if r in set(survive['kept_ids'])]
        kept_d = [dvals[leaves.index(r)] for r in kept]
        restricted = None
        if kept and len(kept) < len(leaves):
            res = s298.theta1_all_edges(edges, kept, kept_d)
            restricted = dict(theta1=res['theta1'], argmin=res['best'])
        row = dict(kind='prune', stage=300, collection=args.collection, name=name,
                   n_observed=len(leaves), n_kept=len(kept),
                   kept_fraction=len(kept) / max(1, len(leaves)),
                   survive=survive,
                   theta1_full=(theta_full[name]['theta1'] if name in theta_full else None),
                   theta1_restricted=(restricted['theta1'] if restricted else None),
                   arms={})
        for arm in args.arms.split(','):
            runs = by_case.get((name, arm))
            if not runs:
                continue
            best = sorted(runs, key=lambda r: (s300an.STAGE_ORDER.get(r['source_stage'], 9),
                                               r.get('source_line') or 0))[0]
            sources, pendant = s300eval.configuration(best)
            full = s300.fit_metrics(edges, leaves, dvals, sources, pendant)
            cell = dict(k=best['k'], reported_loss=best.get('reported_loss'),
                        r_dd_full=full['linf'], ssr_full=full['ssr'])
            if kept:
                part = s300.fit_metrics(edges, kept, kept_d, sources, pendant)
                cell.update(r_dd_kept=part['linf'], ssr_kept=part['ssr'],
                            n_kept=len(kept))
                reference = (restricted['theta1'] if restricted
                             else (theta_full[name]['theta1'] if name in theta_full else None))
                if reference:
                    cell['r_dd_kept_ge_theta1_kept'] = bool(part['linf'] >= reference)
                    cell['r_dd_kept_over_theta1_kept'] = float(part['linf'] / reference)
            row['arms'][arm] = cell
        handle.write(json.dumps(s300.serial(row)) + '\n')
        if position % 25 == 0:
            handle.flush()
            print('[%s] %s %d/%d %.0f s' % (time.strftime('%H:%M:%S'), args.collection,
                                            position, len(names),
                                            time.perf_counter() - t0), flush=True)
    handle.close()
    print('wrote', path, flush=True)


if __name__ == '__main__':
    main()
