#!/usr/bin/env python3
"""Stage 289: reproduction check of the DecoDiPhy example shipped with the tool.

Compares a locally produced all_rounds.json / output.jplace against the copies
committed in the DecoDiPhy repository (produced by the tool's author).
Writes a JSON report; no pre-existing file is modified.
"""
import json
import sys
from pathlib import Path


def load(path):
    return json.loads(Path(path).read_text())


def round_key(r):
    return (r['k'], r['rounds'])


def compare_rounds(published, produced):
    out = []
    for a, b in zip(published, produced):
        anchors_a, anchors_b = list(a['anchors']), list(b['anchors'])
        pa, pb = list(a['p']), list(b['p'])
        xa, xb = list(a['x']), list(b['x'])
        ya = a['y'] if isinstance(a['y'], list) else [a['y']]
        yb = b['y'] if isinstance(b['y'], list) else [b['y']]
        same_set = set(anchors_a) == set(anchors_b)
        same_order = anchors_a == anchors_b
        la, lb = a['loss'], b['loss']
        entry = dict(k_published=a['k'], k_produced=b['k'],
                     rounds_published=a['rounds'], rounds_produced=b['rounds'],
                     anchor_set_identical=same_set, anchor_order_identical=same_order,
                     loss_published=la, loss_produced=lb)
        if isinstance(la, (int, float)) and isinstance(lb, (int, float)):
            entry['loss_abs_diff'] = abs(la - lb)
            entry['loss_rel_diff'] = abs(la - lb) / abs(la) if la else None
        if same_order and len(pa) == len(pb):
            entry['max_abs_p_diff'] = max(abs(u - v) for u, v in zip(pa, pb))
            entry['max_abs_x_diff'] = max(abs(u - v) for u, v in zip(xa, xb))
            entry['max_abs_y_diff'] = max(abs(u - v) for u, v in zip(ya, yb))
        out.append(entry)
    return out


def compare_jplace(published, produced):
    def placements(doc):
        return {q['n']: q['p'][0] for q in doc['placements']}
    a, b = placements(published), placements(produced)
    shared = sorted(set(a) & set(b))
    rows = []
    for q in shared:
        pa, pb = a[q], b[q]
        rows.append(dict(query=q, edge_published=pa[0], edge_produced=pb[0],
                         edge_identical=pa[0] == pb[0],
                         abundance_abs_diff=abs(pa[1] - pb[1]),
                         x_abs_diff=abs(pa[2] - pb[2]),
                         y_abs_diff=abs(pa[3] - pb[3])))
    return dict(n_published=len(a), n_produced=len(b), n_shared=len(shared),
                only_published=sorted(set(a) - set(b)), only_produced=sorted(set(b) - set(a)),
                all_edges_identical=all(r['edge_identical'] for r in rows),
                max_abundance_abs_diff=max((r['abundance_abs_diff'] for r in rows), default=None),
                max_x_abs_diff=max((r['x_abs_diff'] for r in rows), default=None),
                max_y_abs_diff=max((r['y_abs_diff'] for r in rows), default=None),
                per_query=rows)


def main():
    pub_dir, loc_dir, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
    report = dict(published_dir=str(pub_dir), produced_dir=str(loc_dir))
    report['rounds'] = compare_rounds(load(pub_dir / 'all_rounds.json'),
                                      load(loc_dir / 'all_rounds.json'))
    report['n_rounds_published'] = len(load(pub_dir / 'all_rounds.json'))
    report['n_rounds_produced'] = len(load(loc_dir / 'all_rounds.json'))
    report['jplace'] = compare_jplace(load(pub_dir / 'output.jplace'), load(loc_dir / 'output.jplace'))
    pub_tree = (pub_dir / 'labeled_tree.trees')
    loc_tree = (loc_dir / 'labeled_tree.trees')
    if pub_tree.exists() and loc_tree.exists():
        report['labeled_tree_identical'] = pub_tree.read_text().strip() == loc_tree.read_text().strip()
    report['final_anchor_set_identical'] = report['rounds'][-1]['anchor_set_identical']
    out.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({k: v for k, v in report.items() if k not in ('rounds',)}, indent=2)[:3000])
    print('--- last 3 rounds ---')
    print(json.dumps(report['rounds'][-3:], indent=2))


if __name__ == '__main__':
    main()
