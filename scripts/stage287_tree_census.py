#!/usr/bin/env python3
"""Stage 287 (part 2) -- census of every fetched reference tree.

Counts leaves from the Newick itself (not from commas), records exact-rational
edge-length pathologies that bear on the PDD model's positivity requirement,
and writes results/stage287_tree_census.json.
"""
import hashlib, json, os, sys
from datetime import datetime, timezone
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
from stage286_newick import parse_newick
DEST = os.path.join(ROOT, 'data', 'stage287_instances')

def sha256_file(p):
    h = hashlib.sha256()
    with open(p, 'rb') as fh:
        for b in iter(lambda: fh.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()

rows = {}
for group in ['demixing_10k/trees', 'biotrees/trees', 'sixteen_k',
              'travelers_diarrhea', 'krepp_shared']:
    gdir = os.path.join(DEST, group)
    if not os.path.isdir(gdir):
        continue
    for fn in sorted(os.listdir(gdir)):
        p = os.path.join(gdir, fn)
        if not os.path.isfile(p) or not fn.endswith(('.tre', '.nwk', '.rerooted', '.trees')):
            continue
        try:
            with open(p) as fh:
                txt = fh.read().strip().split('\n')[0]
            if not txt:
                rows[group + '/' + fn] = {'error': 'empty file', 'bytes': os.path.getsize(p)}
                continue
            t, diag = parse_newick(txt)
            lens = [t.elen[i] for i in range(t.n) if i != t.root and t.elen[i] is not None]
            labs = [t.label[i] for i in range(t.n) if not t.children[i]]
            rows[group + '/' + fn] = {
                'bytes': os.path.getsize(p), 'sha256': sha256_file(p),
                'n_leaves': diag['n_leaves'], 'n_vertices': diag['n_vertices'],
                'n_negative_edges': diag['n_negative_edges'],
                'n_zero_edges': diag['n_zero_edges'],
                'n_missing_edges': diag['n_missing_edges'],
                'n_unifurcations': diag['n_unifurcations'],
                'root_degree': diag['root_degree'],
                'n_duplicate_leaf_labels': diag['n_duplicate_leaf_labels'],
                'n_unlabelled_leaves': diag['n_unlabelled_leaves'],
                'min_edge_length': str(min(lens)) if lens else None,
                'max_edge_length': str(max(lens)) if lens else None,
                'pdd_positivity_ok': (diag['n_negative_edges'] == 0
                                      and diag['n_zero_edges'] == 0
                                      and diag['n_missing_edges'] == 0),
            }
        except Exception as exc:                            # noqa: BLE001
            rows[group + '/' + fn] = {'error': '%s: %s' % (type(exc).__name__, exc),
                                      'bytes': os.path.getsize(p)}
out = {'stage': 287, 'part': 'tree_census',
       'registered_utc': datetime.now(timezone.utc).isoformat(),
       'n_trees': len(rows), 'trees': rows}
with open(os.path.join(ROOT, 'results', 'stage287_tree_census.json'), 'w') as fh:
    json.dump(out, fh, indent=1)
ok = [(v['n_leaves'], k) for k, v in rows.items() if 'n_leaves' in v]
ok.sort(reverse=True)
print('trees parsed: %d (errors %d)' % (len(ok), sum(1 for v in rows.values() if 'error' in v)))
print('LARGEST 15 by leaf count:')
for n, k in ok[:15]:
    v = rows[k]
    print('  %7d leaves  neg=%-5d zero=%-4d dup=%-3d  %s' % (n, v['n_negative_edges'], v['n_zero_edges'], v['n_duplicate_leaf_labels'], k))
print('errors:', {k: v['error'] for k, v in rows.items() if 'error' in v})
