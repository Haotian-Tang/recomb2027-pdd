#!/usr/bin/env python3
"""Stage 307: how far apart are the two branch-length calibrations of one ANI tree?

The instance bank pairs the benchmark's estimated ANI distance vectors with the
positive-length backbone tree_X.tre rather than with the ANI-calibrated
tree_X-ani_blen.tre the benchmark ships, because the latter carries negative
branch lengths and is not a metric tree.  This script measures what that choice
changes, so the manuscript can state it rather than leave it in the build log.

No source count is recomputed and the negative-length tree is never solved; the
only question asked is whether the two calibrations induce the same geometry.
Writes results/stage307/stage307_tree_pairing.json.

The upstream Newick inputs are not bundled in git or the bulk archive. See
README's tree-pairing check for a hash-verified download of just these files.
"""
import json
import os
import re
import statistics
import sys
from collections import defaultdict

TREES = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir,
                     'data', 'stage287_instances', 'demixing_10k', 'trees')


def parse_newick(text):
    """Return (adjacency, leaf_names).  Nodes are ints; edges carry lengths.

    Square-bracket comments (the `[&R]` rooting marker these files carry) are
    stripped first; otherwise the marker parses as a leaf label.
    """
    text = re.sub(r'\[[^\]]*\]', '', text).strip().rstrip(';')
    adj = defaultdict(list)
    leaves = {}
    nxt = [0]

    def new_node():
        nxt[0] += 1
        return nxt[0] - 1

    def parse(i, parent):
        if text[i] == '(':
            me = new_node()
            i += 1
            while True:
                i = parse(i, me)
                if text[i] == ',':
                    i += 1
                    continue
                if text[i] == ')':
                    i += 1
                    break
                raise ValueError('bad newick at %d' % i)
        else:
            me = None
        j = i
        while j < len(text) and text[j] not in '(),:;':
            j += 1
        label = text[i:j]
        i = j
        if me is None:
            me = new_node()
            leaves[me] = label
        length = 0.0
        if i < len(text) and text[i] == ':':
            j = i + 1
            while j < len(text) and (text[j].isdigit() or text[j] in '.-+eE'):
                j += 1
            length = float(text[i + 1:j])
            i = j
        if parent is not None:
            adj[parent].append((me, length))
            adj[me].append((parent, length))
        return i

    parse(0, None)
    return adj, leaves


def all_pairs(adj, leaves):
    """Path-length sums between every pair of leaves, by label."""
    out = {}
    ids = sorted(leaves)
    for s in ids:
        dist = {s: 0.0}
        stack = [s]
        while stack:
            u = stack.pop()
            for v, w in adj[u]:
                if v not in dist:
                    dist[v] = dist[u] + w
                    stack.append(v)
        for t in ids:
            if leaves[s] < leaves[t]:
                out[(leaves[s], leaves[t])] = dist[t]
    return out


def edge_set(adj):
    seen = set()
    for u, nbrs in adj.items():
        for v, w in nbrs:
            key = (min(u, v), max(u, v))
            if key not in seen:
                seen.add(key)
                yield w


def topology_splits(adj, leaves):
    """Unrooted edge bipartitions, ignoring lengths and degree-two roots.

    A set removes duplicate splits induced by subdivisions or a degree-two
    Newick root. Zero-length edges are retained as combinatorial edges.
    """
    names = sorted(leaves.values())
    if len(set(names)) != len(names):
        raise ValueError('duplicate leaf labels')
    bits = {name: 1 << i for i, name in enumerate(names)}
    full = (1 << len(names)) - 1
    root = next(iter(adj))
    parent, order = {root: None}, [root]
    for u in order:
        for v, _ in adj[u]:
            if v not in parent:
                parent[v] = u
                order.append(v)
    masks = {}
    for u in reversed(order):
        mask = bits[leaves[u]] if u in leaves else 0
        for v, _ in adj[u]:
            if parent.get(v) == u:
                mask |= masks[v]
        masks[u] = mask
    return {min(masks[v], full ^ masks[v]) for v in order[1:]}


def main():
    stems = sys.argv[1:] or ['tree_10']
    missing = [stem + suffix + '.tre' for stem in stems
               for suffix in ('', '-ani_blen')
               if not os.path.isfile(os.path.join(TREES, stem + suffix + '.tre'))]
    if missing:
        raise SystemExit('missing upstream Newick inputs: %s; see the tree-pairing '
                         'download instructions in README.md' % ', '.join(missing))
    rows = []
    for stem in stems:
        a, la = parse_newick(open(os.path.join(TREES, stem + '.tre')).read())
        b, lb = parse_newick(open(os.path.join(TREES, stem + '-ani_blen.tre')).read())
        names_a, names_b = sorted(la.values()), sorted(lb.values())
        same_labels = names_a == names_b
        if not same_labels:
            raise ValueError('leaf label mismatch for %s' % stem)
        same_topology = topology_splits(a, la) == topology_splits(b, lb)
        if not same_topology:
            raise ValueError('unrooted topology mismatch for %s' % stem)
        wb = list(edge_set(b))
        pa, pb = all_pairs(a, la), all_pairs(b, lb)
        shared = sorted(set(pa) & set(pb))
        ratios = [pb[k] / pa[k] for k in shared if pa[k] > 0]
        reldiff = [abs(pb[k] - pa[k]) / pa[k] for k in shared if pa[k] > 0]
        rows.append(dict(
            stem=stem, n_leaves=len(la), same_leaf_labels=same_labels,
            same_topology=same_topology,
            n_pairs=len(shared), n_edges_ani=len(wb),
            n_negative_edges_ani=sum(1 for w in wb if w < 0),
            n_zero_edges_ani=sum(1 for w in wb if w == 0),
            median_path_ratio=statistics.median(ratios),
            median_abs_rel_diff=statistics.median(reldiff),
        ))
        r = rows[-1]
        print('%-10s leaves %3d  labels match %s  pairs %5d | ANI tree: %2d negative, '
              '%d zero of %d edges' % (stem, r['n_leaves'], r['same_leaf_labels'],
                                       r['n_pairs'], r['n_negative_edges_ani'],
                                       r['n_zero_edges_ani'], r['n_edges_ani']))
        print('%-10s median ANI/original path-sum ratio %.4f, median |relative '
              'difference| %.2f%%' % ('', r['median_path_ratio'],
                                      100 * r['median_abs_rel_diff']))
    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir,
                           'results', 'stage307')
    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, 'stage307_tree_pairing.json'), 'w') as fh:
        json.dump(dict(stage=307, kind='tree_pairing_geometry', trees=rows),
                  fh, indent=1, sort_keys=True)


if __name__ == '__main__':
    main()
