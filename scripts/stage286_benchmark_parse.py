#!/usr/bin/env python3
"""Stage 286 -- parse the parent work's public benchmark (read-only inventory).

Establishes, for the DecoDiPhy / Demixing public data already on disk under
``data/stage29_distance_target`` (plus the related earlier public inputs in
``data/stage08_distance_public`` and ``data/stage15_public_pairs``):

  * the real column schema of every aggregate distance file, taken from the
    authors' own reader code rather than guessed;
  * how many distinct query instances exist, keyed by (tree, K, replicate,
    abundance distribution);
  * each instance's reference-set size, counted as table rows and, where the
    tree is on disk, cross-checked against a proper Newick leaf count;
  * the true source count and true abundances;
  * which distance estimators exist for each instance, and in particular which
    instances carry an exact/true distance vector rather than only an estimate.

Writes ``results/stage286_benchmark_parse.json``. Reads only; creates nothing
else.
"""
import collections
import hashlib
import json
import os
import re
import sys
from datetime import datetime, timezone
from fractions import Fraction

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
from stage286_newick import parse_newick  # noqa: E402

S29 = os.path.join(ROOT, 'data', 'stage29_distance_target')
S08 = os.path.join(ROOT, 'data', 'stage08_distance_public')
S15 = os.path.join(ROOT, 'data', 'stage15_public_pairs')
OUT = os.path.join(ROOT, 'results', 'stage286_benchmark_parse.json')

# --------------------------------------------------------------------------
# Column schemas, each with the evidence that establishes it.
# --------------------------------------------------------------------------
SCHEMAS = {
    'per_instance/distances.txt': {
        'columns': ['REF_ID', 'TRUE_DISTANCE'],
        'separator': 'whitespace',
        'header': False,
        'evidence': ('Two columns observed in every file; values carry ~17 '
                     'significant digits consistent with a float64 repr of a '
                     'tree-derived distance. DecoDiPhy cli.py reads a bare '
                     '(label, distance) table in exactly this form '
                     '(decodiphy/cli.py: distances = {l[0]: float(l[1])}).'),
    },
    'per_instance/estimated_distances_pr-h14_dth4.txt': {
        'columns': ['REF_ID', 'ESTIMATED_DISTANCE', 'PORTION_MAPPED',
                    'TRUE_DISTANCE'],
        'separator': 'tab',
        'header': True,
        'evidence': ('Literal header row in the file; confirmed by '
                     'find_sol_10k.py, which reads this file and takes '
                     'est = float(d[1]), portion = float(d[2]) and '
                     'true = float(d[3]).'),
    },
    'Demixing__10k_dataset__dist.txt': {
        'columns': ['tree', 'k', 's', 'dist', 'REF_ID', 'ESTIMATED_DISTANCE',
                    'PORTION_MAPPED', 'TRUE_DISTANCE'],
        'separator': 'tab',
        'header': False,
        'evidence': ('Four key columns prepended by commands.sh-style awk '
                     '(-v var1=tree var2=k var3=s var4=dist); the remaining '
                     'four are the per-instance estimated-distance table, '
                     'whose header is spelled out verbatim in '
                     'krepp_distances.txt.'),
    },
    'Demixing__10k_dataset__krepp_distances.txt': {
        'columns': ['tree', 'k', 's', 'dist', 'REF_ID', 'ESTIMATED_DISTANCE',
                    'PORTION_MAPPED', 'TRUE_DISTANCE'],
        'separator': 'tab',
        'header': 'one header row per instance group',
        'evidence': 'Header row repeated inside the file for every group.',
    },
    'Demixing__10k_dataset__krepp_placement_distances.txt': {
        'columns': ['tree', 'k', 's', 'dist', 'REF_ID', 'PLACEMENT_DISTANCE'],
        'separator': 'tab',
        'header': False,
        'evidence': ('Six columns; row count per group equals the FULL tree '
                     'leaf count (queries included), unlike the other files '
                     'whose row count is full-minus-K.'),
    },
    'Demixing__10k_dataset__truemash.txt': {
        'columns': ['tree', 'k', 's', 'dist', 'REF_ID', 'TRUE_MASH_DISTANCE'],
        'separator': 'tab',
        'header': False,
        'evidence': ('Six columns; equals the per-instance '
                     'true_mash_distances.txt table with the four key columns '
                     'prepended. find_sol_10k.py writes that per-instance file '
                     'as (label, true distance) computed on the ANI '
                     'branch-length tree.'),
    },
    'Demixing__10k_dataset__all_metrics_est_mash_corrected.txt': {
        'columns': ['tree', 'k', 's', 'dist', 'jaccard', 'unifrac', 'wunifrac',
                    'est_k', 'wasser', 'bc'],
        'separator': 'tab',
        'header': False,
        'evidence': ("Named verbatim by the authors' plotting script "
                     "draw_10k.R: names(est_mash_corrected_data) <- c('tree', "
                     "'k', 's', 'dist', 'jaccard', 'unifrac', 'wunifrac', "
                     "'est_k', 'wasser', 'bc')."),
    },
    'per_instance/read_counts_table-h14_dth.txt': {
        'columns': ['#FeatureID', 'ESTIMATED'],
        'separator': 'tab',
        'header': True,
        'evidence': ('Literal header "#FeatureID\\tESTIMATED"; values are a '
                     'read-mapping-derived relative abundance ESTIMATE over '
                     'reference genomes, not the ground-truth mixture.'),
    },
    'per_instance/queries.txt': {
        'columns': ['QUERY_LABEL', 'TRUE_ABUNDANCE'],
        'separator': 'whitespace',
        'header': False,
        'evidence': ('find_sol_10k.py reads queries.txt and takes '
                     'true_p = [float(q[1]) for q in queries] with q[0] the '
                     'query taxon label; k = len(queries).'),
    },
}


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for blk in iter(lambda: fh.read(1 << 20), b''):
            h.update(blk)
    return h.hexdigest()


def read_rows(path, sep=None):
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line.strip():
                continue
            rows.append(line.split('\t') if sep == '\t' else line.split())
    return rows


# --------------------------------------------------------------------------
# 1. local reference trees
# --------------------------------------------------------------------------
def parse_local_trees():
    out = {}
    for fn in sorted(os.listdir(S29)):
        if not fn.endswith('__original.tre'):
            continue
        tid = fn.split('__')[0]
        path = os.path.join(S29, fn)
        with open(path) as fh:
            txt = fh.read().strip().split('\n')[0]
        t, diag = parse_newick(txt)
        lv = [i for i in range(t.n) if not t.children[i]]
        lens = [t.elen[i] for i in range(t.n) if i != t.root and t.elen[i] is not None]
        out[tid] = {
            'file': os.path.relpath(path, ROOT),
            'sha256': sha256_file(path),
            'newick_diagnostics': diag,
            'n_leaves_counted_from_newick': len(lv),
            'n_commas_in_file': txt.count(','),
            'min_edge_length': str(min(lens)) if lens else None,
            'max_edge_length': str(max(lens)) if lens else None,
            'has_negative_edges': diag['n_negative_edges'] > 0,
            'leaf_labels_sample': sorted(t.label[i] for i in lv)[:5],
        }
    return out


# --------------------------------------------------------------------------
# 2. aggregate files
# --------------------------------------------------------------------------
AGG = {
    'dist': 'Demixing__10k_dataset__dist.txt',
    'truemash': 'Demixing__10k_dataset__truemash.txt',
    'krepp': 'Demixing__10k_dataset__krepp_distances.txt',
    'krepp_placement': 'Demixing__10k_dataset__krepp_placement_distances.txt',
    'all_metrics_est_mash_corrected':
        'Demixing__10k_dataset__all_metrics_est_mash_corrected.txt',
}


def parse_aggregates():
    agg = {}
    groups = {}
    for tag, fn in AGG.items():
        path = os.path.join(S29, fn)
        rows = read_rows(path, sep='\t')
        g = collections.OrderedDict()
        ncols = collections.Counter()
        nheader = 0
        short = []
        for r in rows:
            ncols[len(r)] += 1
            if len(r) > 4 and r[4] == 'REF_ID':
                nheader += 1
                continue
            key = (r[0], r[1], r[2], r[3])
            g.setdefault(key, []).append(r[4:])
            if tag == 'dist' and len(r) == 7:
                short.append(key)
        agg[tag] = {
            'file': fn,
            'sha256': sha256_file(path),
            'n_lines': len(rows),
            'n_columns_histogram': dict(ncols),
            'n_embedded_header_rows': nheader,
            'n_instance_groups': len(g),
            'tree_ids': sorted({k[0] for k in g}, key=int),
            'k_values': sorted({int(k[1]) for k in g}),
            'replicates': sorted({int(k[2]) for k in g}),
            'distributions': sorted({k[3] for k in g}),
            'rows_per_group_min': min(len(v) for v in g.values()),
            'rows_per_group_max': max(len(v) for v in g.values()),
            'n_groups_with_short_rows': len(set(short)),
            'short_row_groups': sorted(set(short))[:20],
        }
        groups[tag] = g
    return agg, groups


# --------------------------------------------------------------------------
# 3. instance inventory from the GitHub manifests
# --------------------------------------------------------------------------
def manifest_inventory():
    inv = {}
    dd = json.load(open(os.path.join(S29, 'data_10k_root.json')))
    dm = json.load(open(os.path.join(S29, 'Demixing_Data.tree.json')))

    inst = collections.defaultdict(set)
    for e in dm['tree']:
        p = e['path'].split('/')
        if p[0] == '10k_dataset' and len(p) == 3 and e['type'] == 'blob' and \
                re.match(r'tree_\d+_k\d+_s\d+_(uniform|exp)$', p[1]):
            inst[p[1]].add(p[2])
    rows = {}
    for name, files in inst.items():
        m = re.match(r'tree_(\d+)_k(\d+)_s(\d+)_(\w+)$', name)
        rows[name] = {
            'tree': m.group(1), 'K': int(m.group(2)),
            'replicate': int(m.group(3)), 'distribution': m.group(4),
            'files': sorted(files),
        }
    inv['demixing_10k_instances'] = rows

    trees10k = {}
    for e in dm['tree']:
        p = e['path']
        if p.startswith('10k_dataset/tree_') and p.endswith(('.tre', '.rerooted')) \
                and p.count('/') == 1 and e['type'] == 'blob':
            trees10k[p.split('/')[-1]] = e['size']
    inv['demixing_10k_tree_blobs'] = trees10k
    inv['decodiphy_data_10k_tree_blobs'] = {
        e['path']: e['size'] for e in dd['tree']
        if e['type'] == 'blob' and e['path'].endswith(('.tre', '.rerooted'))}
    inv['decodiphy_data_10k_instance_dirs'] = sorted(
        e['path'] for e in dd['tree'] if e['type'] == 'tree'
        and re.match(r'tree_\d+_k\d+_s\d+_(uniform|exp)$', e['path']))

    bio = collections.defaultdict(set)
    for e in dm['tree']:
        p = e['path'].split('/')
        if p[0] == 'biotrees' and len(p) == 4 and e['type'] == 'blob':
            bio[p[1] + '/' + p[2]].add(p[3])
    brows = {}
    for name, files in bio.items():
        ds, inm = name.split('/')
        m = re.search(r'^(.*)_k(\d+)_s(\d+)_noise(\d)(?:_l(\d+))?$', inm)
        if not m:
            brows[name] = {'dataset': ds, 'unparsed': True, 'files': sorted(files)}
            continue
        brows[name] = {
            'dataset': ds, 'tree_file': m.group(1), 'K': int(m.group(2)),
            'seed': int(m.group(3)), 'noise_type': int(m.group(4)),
            'sequence_length': int(m.group(5)) if m.group(5) else None,
            'files': sorted(files),
        }
    inv['biotrees_instances'] = brows
    inv['biotrees_tree_blobs'] = {
        e['path'].split('/')[-1]: e['size'] for e in dm['tree']
        if e['type'] == 'blob' and e['path'].count('/') == 1
        and e['path'].startswith('biotrees/') and e['path'].endswith(('.tre', '.nwk'))}
    inv['sixteen_k'] = {
        'backbone': [(e['path'], e['size']) for e in dm['tree']
                     if e['path'].startswith('16k/') and e['path'].endswith('.nwk')],
        'query_instances': sorted(e['path'] for e in dm['tree']
                                 if e['path'].startswith('16k/simulated/')
                                 and e['path'].endswith('query.txt')),
    }
    inv['travelers_diarrhea'] = sorted(
        (e['path'], e['size']) for e in dm['tree']
        if e['path'].startswith('travelers_diarrhea/') and e['type'] == 'blob'
        and e['path'].count('/') == 1)
    return inv


# --------------------------------------------------------------------------
# 4. local per-instance files
# --------------------------------------------------------------------------
def parse_local_instances():
    out = {}
    for tid in ['tree_4', 'tree_7', 'tree_15', 'tree_17', 'tree_1322']:
        rec = {}
        dpath = os.path.join(S29, tid + '__distances.txt')
        rpath = os.path.join(S29, tid + '__read_counts.txt')
        if os.path.exists(dpath):
            rows = read_rows(dpath)
            rec['distances'] = {
                'file': os.path.basename(dpath), 'sha256': sha256_file(dpath),
                'n_rows': len(rows), 'n_columns': sorted({len(r) for r in rows}),
                'interpretation': 'TRUE distances (tree-derived), float-rounded',
            }
        if os.path.exists(rpath):
            rows = read_rows(rpath)
            body = [r for r in rows if not r[0].startswith('#')]
            tot = sum(Fraction(r[1]) for r in body)
            rec['read_counts'] = {
                'file': os.path.basename(rpath), 'sha256': sha256_file(rpath),
                'n_rows': len(body), 'sum_of_values': str(tot),
                'sum_as_float': float(tot),
                'interpretation': ('krepp read-mapping ESTIMATE of relative '
                                   'abundance over reference genomes; NOT the '
                                   'ground-truth mixture'),
            }
        out[tid + '_k3_s1_uniform'] = rec
    # stage08 copies carry queries.txt => ground truth for 5 instances
    s08 = {}
    if os.path.isdir(S08):
        for d in sorted(os.listdir(S08)):
            dd = os.path.join(S08, d)
            if not os.path.isdir(dd):
                continue
            rec = {'files': sorted(os.listdir(dd))}
            qp = os.path.join(dd, 'queries.txt')
            if os.path.exists(qp):
                q = read_rows(qp)
                rec['queries'] = [[r[0], r[1]] for r in q]
                rec['K_true'] = len(q)
                rec['true_abundance_sum'] = str(sum(Fraction(r[1]) for r in q))
            tp = os.path.join(dd, 'true_mash_distances.txt')
            if os.path.exists(tp):
                rec['n_true_mash_rows'] = len(read_rows(tp))
            ep = os.path.join(dd, 'estimated_distances_pr-h14_dth4.txt')
            if os.path.exists(ep):
                rows = read_rows(ep, sep='\t')
                rec['estimated_header'] = rows[0]
                rec['n_estimated_rows'] = len(rows) - 1
            ap = os.path.join(dd, 'ani.tre')
            if os.path.exists(ap):
                with open(ap) as fh:
                    t, diag = parse_newick(fh.read().strip().split('\n')[0])
                rec['ani_tree_diagnostics'] = diag
                rec['ani_tree_sha256'] = sha256_file(ap)
            s08[d] = rec
    return out, s08


def main():
    res = {
        'stage': 286,
        'registered_utc': datetime.now(timezone.utc).isoformat(),
        'purpose': ('Read-only inventory and schema determination for the '
                    'parent work public benchmark; no instance is built and no '
                    'solver is run here.'),
        'sources': json.load(open(os.path.join(S29, 'downloads.json'))),
        'file_schemas': SCHEMAS,
    }
    res['local_reference_trees'] = parse_local_trees()
    agg, groups = parse_aggregates()
    res['aggregate_files'] = agg

    # per-instance reference-set sizes and estimator availability
    inv = manifest_inventory()
    per = {}
    for tag in ['dist', 'truemash', 'krepp', 'krepp_placement']:
        for key, rows in groups[tag].items():
            name = 'tree_%s_k%s_s%s_%s' % (key[0], key[1], key[2], key[3])
            e = per.setdefault(name, {'tree': key[0], 'K': int(key[1]),
                                      'replicate': int(key[2]),
                                      'distribution': key[3],
                                      'estimators': {}})
            e['estimators'][tag] = len(rows)
    for key, rows in groups['all_metrics_est_mash_corrected'].items():
        name = 'tree_%s_k%s_s%s_%s' % (key[0], key[1], key[2], key[3])
        if name in per:
            per[name]['published_est_k'] = rows[0][3].strip()
            per[name]['published_metrics'] = {
                'jaccard': rows[0][0].strip(), 'unifrac': rows[0][1].strip(),
                'wunifrac': rows[0][2].strip(), 'wasser': rows[0][4].strip(),
                'bc': rows[0][5].strip()}
    for name, rec in inv['demixing_10k_instances'].items():
        e = per.setdefault(name, {'tree': rec['tree'], 'K': rec['K'],
                                  'replicate': rec['replicate'],
                                  'distribution': rec['distribution'],
                                  'estimators': {}})
        e['repo_files'] = rec['files']
        e['has_ground_truth_queries'] = 'queries.txt' in rec['files']
        e['has_true_distance_table'] = 'distances.txt' in rec['files']
        e['has_krepp_estimated_table'] = \
            'estimated_distances_pr-h14_dth4.txt' in rec['files']
        e['has_true_mash_table'] = 'true_mash_distances.txt' in rec['files']
    # reference-set size and implied full-tree leaf count
    for name, e in per.items():
        n = e['estimators'].get('dist') or e['estimators'].get('krepp') \
            or e['estimators'].get('truemash')
        e['reference_set_size_from_tables'] = n
        e['implied_full_tree_leaves'] = (n + e['K']) if n else None
        pl = e['estimators'].get('krepp_placement')
        e['full_tree_leaves_from_placement_table'] = pl
    res['per_instance_10k'] = per
    res['inventory'] = inv

    # cross-check implied full-tree leaf counts against on-disk Newick counts
    checks = []
    for tid, rec in res['local_reference_trees'].items():
        nl = rec['n_leaves_counted_from_newick']
        implied = sorted({e['implied_full_tree_leaves'] for e in per.values()
                          if e['tree'] == tid.split('_')[1]
                          and e['implied_full_tree_leaves']})
        checks.append({'tree': tid, 'newick_leaf_count': nl,
                       'implied_full_tree_leaves_from_tables': implied,
                       'consistent': implied == [nl] if implied else None})
    res['leaf_count_cross_check'] = checks

    loc, s08 = parse_local_instances()
    res['local_per_instance_files'] = loc
    res['stage08_public_instances'] = s08
    res['stage15_public_pairs'] = json.load(open(os.path.join(S15, 'source_manifest.json')))

    sizes = sorted({e['reference_set_size_from_tables'] for e in per.values()
                    if e['reference_set_size_from_tables']})
    res['summary'] = {
        'n_query_instances_with_any_distance_table': len(per),
        'n_query_instances_in_repo_with_ground_truth': sum(
            1 for e in per.values() if e.get('has_ground_truth_queries')),
        'n_with_true_distance_table': sum(
            1 for e in per.values() if e.get('has_true_distance_table')),
        'n_with_krepp_estimated_table': sum(
            1 for e in per.values() if e.get('has_krepp_estimated_table')),
        'n_with_true_mash_table': sum(
            1 for e in per.values() if e.get('has_true_mash_table')),
        'reference_set_sizes_observed': sizes,
        'reference_set_size_min': sizes[0], 'reference_set_size_max': sizes[-1],
        'K_values': sorted({e['K'] for e in per.values()}),
        'n_biotrees_instances': len(inv['biotrees_instances']),
        'n_16k_query_instances': len(inv['sixteen_k']['query_instances']),
    }
    with open(OUT, 'w') as fh:
        json.dump(res, fh, indent=1, sort_keys=False)
    print(json.dumps(res['summary'], indent=1))
    print('leaf cross-check:', json.dumps(res['leaf_count_cross_check'], indent=1))
    print('wrote', OUT)


if __name__ == '__main__':
    main()
