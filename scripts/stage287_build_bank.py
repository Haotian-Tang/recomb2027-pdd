#!/usr/bin/env python3
"""Stage 287 (part 3) -- build and freeze the exact PDD instance bank.

Converts the public mixture instances fetched in part 1 into the project's
exact rational case format (the format of
``results/stage44_confirmation_protocol.json``) and freezes the result as

    data/stage287_instances/bank/cases/<family>/<case>.json
    results/stage287_instance_bank_protocol.json

Two parallel families are produced wherever the data allow:

  * an EXACT family, whose observed vector is recomputed from the real tree and
    the published generating mixture in ``Fraction`` arithmetic, and
  * an ESTIMATED/PUBLISHED family, whose observed vector is the published
    decimal table read with ``Fraction(str)``; these are not exactly
    explainable and are labelled as such.

No solver is run. Every exact-family case is validated by requiring that the
stated generating mixture reproduce every observed leaf distance exactly, with
the full-tree and pruned-tree computations as independent code paths.
"""
import collections
import hashlib
import json
import os
import re
import sys
import time
from datetime import datetime, timezone
from fractions import Fraction

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
from stage286_newick import parse_newick                      # noqa: E402
from stage287_pdd_cases import (build_exact_case, CaseError,   # noqa: E402
                                undirected, dists_from, leaf_labels)

DEST = os.path.join(ROOT, 'data', 'stage287_instances')
BANK = os.path.join(DEST, 'bank')
CASES = os.path.join(BANK, 'cases')
PROTOCOL = os.path.join(ROOT, 'results', 'stage287_instance_bank_protocol.json')
SIXTEEN_K_MAX = 20            # frozen rule: materialise 16k cases with K <= 20

_TREE_CACHE = {}


def load_tree(path):
    if path not in _TREE_CACHE:
        with open(path) as fh:
            txt = fh.read().strip().split('\n')[0]
        _TREE_CACHE[path] = parse_newick(txt)[0]
    return _TREE_CACHE[path]


def sha256_bytes(b):
    return hashlib.sha256(b).hexdigest()


def read_two_col(path):
    out = {}
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) >= 2 and not f[0].startswith('#'):
                out[f[0]] = f[1]
    return out


_NUMRE = re.compile(r'^[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?$')


def is_num(v):
    return bool(v) and bool(_NUMRE.match(v.strip()))


def read_estimated(path):
    """REF_ID / ESTIMATED_DISTANCE / PORTION_MAPPED / TRUE_DISTANCE."""
    est, por, tru = {}, {}, {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) < 2 or f[0] == 'REF_ID':
                continue
            est[f[0]] = f[1]
            if len(f) > 2:
                por[f[0]] = f[2]
            if len(f) > 3:
                tru[f[0]] = f[3]
    return est, por, tru


def exact_normalised(weights):
    """Exact rational weights summing to 1 from decimal/float inputs."""
    fr = {k: Fraction(v) if isinstance(v, str) else Fraction(v)
          for k, v in weights.items()}
    tot = sum(fr.values())
    if tot == 0:
        raise CaseError('zero_weight_total')
    return {k: v / tot for k, v in fr.items()}, tot


def write_case(family, case, provenance):
    d = os.path.join(CASES, family)
    os.makedirs(d, exist_ok=True)
    rec = dict(case)
    rec['provenance'] = provenance
    blob = json.dumps(rec, indent=1, sort_keys=False).encode()
    path = os.path.join(d, case['name'] + '.json')
    with open(path, 'wb') as fh:
        fh.write(blob)
    return {'family': family, 'name': case['name'],
            'path': os.path.relpath(path, ROOT), 'bytes': len(blob),
            'sha256': sha256_bytes(blob),
            'reference_set_size': case['observed_leaves'],
            'K_generating': case.get('generating_count'),
            'cap': case['cap'],
            'tree_leaves': provenance.get('full_tree_leaves'),
            'distance_source': provenance['distance_source'],
            'exactly_explained': provenance['exactly_explained']}


# ---------------------------------------------------------------------------
# 10k dataset
# ---------------------------------------------------------------------------
def build_10k(index, quarantine, census):
    idir = os.path.join(DEST, 'demixing_10k', 'instances')
    tdir = os.path.join(DEST, 'demixing_10k', 'trees')
    names = sorted(os.listdir(idir))
    for nm in names:
        m = re.match(r'tree_(\d+)_k(\d+)_s(\d+)_(uniform|exp)$', nm)
        if not m:
            continue
        tid, K, rep, dist = m.group(1), int(m.group(2)), int(m.group(3)), m.group(4)
        tpath = os.path.join(tdir, 'tree_%s.tre' % tid)
        qpath = os.path.join(idir, nm, 'queries.txt')
        dpath = os.path.join(idir, nm, 'distances.txt')
        epath = os.path.join(idir, nm, 'estimated_distances_pr-h14_dth4.txt')
        if not (os.path.exists(tpath) and os.path.exists(qpath)):
            quarantine.append({'instance': nm, 'family': 'exact_10k',
                               'reason': 'missing_input_file'})
            continue
        cen = census['trees'].get('demixing_10k/trees/tree_%s.tre' % tid, {})
        tree = load_tree(tpath)
        raw = {k: v for k, v in read_two_col(qpath).items() if is_num(v)}
        try:
            qn, tot = exact_normalised(raw)
            case = build_exact_case(tree, qn, cap=K, name='k10k_%s_exact' % nm)
        except CaseError as exc:
            quarantine.append({
                'instance': nm, 'family': 'exact_10k', 'reason': str(exc),
                'tree': 'tree_%s.tre' % tid,
                'tree_negative_edges': cen.get('n_negative_edges'),
                'tree_zero_edges': cen.get('n_zero_edges')})
            continue
        prov = {
            'collection': 'DecoDiPhy/Demixing 10k_dataset',
            'repo': 'shayesteh99/Demixing_Data',
            'commit': '5d3a1141158994af3fd6e0712fe6aaa14b52903f',
            'instance': nm, 'tree_file': 'tree_%s.tre' % tid,
            'tree_sha256': cen.get('sha256'),
            'full_tree_leaves': cen.get('n_leaves'),
            'K_true': K, 'replicate': rep, 'abundance_distribution': dist,
            'ground_truth': 'queries.txt (query taxon, true abundance)',
            'weights_renormalised_exactly': str(tot) != '1',
            'distance_source': 'exact_tree_derived',
            'exactly_explained': True,
        }
        index.append(write_case('exact_10k', case, prov))
        lbl2id = {v: k for k, v in case['leaf_labels'].items()}

        if os.path.exists(dpath):
            pub = {l: v for l, v in read_two_col(dpath).items() if is_num(v)}
            if set(pub) == set(lbl2id):
                c2 = dict(case)
                c2['name'] = 'k10k_%s_published_true' % nm
                c2['observed'] = {lbl2id[l]: str(Fraction(v))
                                  for l, v in pub.items()}
                c2.pop('generating_pendant', None)
                p2 = dict(prov)
                p2.update(distance_source='published_true_rounded',
                          exactly_explained=False,
                          note=('Authors published true distances as float '
                                'reprs; they differ from the exact rational '
                                'values at the 1e-16 level.'),
                          source_file='distances.txt')
                index.append(write_case('published_true_10k', c2, p2))
            else:
                quarantine.append({'instance': nm,
                                   'family': 'published_true_10k',
                                   'reason': 'leaf_set_mismatch'})
        if os.path.exists(epath):
            est, por, _t = read_estimated(epath)
            if set(est) == set(lbl2id):
                good = {l: v for l, v in est.items() if is_num(v)}
                missing = sorted(int(lbl2id[l]) for l in est
                                 if not is_num(est[l]))
                c3 = dict(case)
                c3['name'] = 'k10k_%s_krepp_estimated' % nm
                c3['observed'] = {lbl2id[l]: str(Fraction(v))
                                  for l, v in good.items()}
                c3['observed_leaves'] = len(c3['observed'])
                c3['hidden_leaves'] = missing
                c3['portion_mapped'] = {lbl2id[l]: por.get(l)
                                        for l in good}
                c3.pop('generating_pendant', None)
                p3 = dict(prov)
                p3.update(distance_source='krepp_estimated_ani',
                          exactly_explained=False,
                          n_references_without_estimate=len(missing),
                          source_file='estimated_distances_pr-h14_dth4.txt',
                          note=('ANI distances estimated from reads by krepp. '
                                'The authors pair these with the ANI '
                                'branch-length tree tree_%s-ani_blen.tre, '
                                'which carries negative branch lengths and is '
                                'therefore not PDD-conformant; the bank pairs '
                                'them with the positive-length tree instead.'
                                % tid))
                index.append(write_case('krepp_estimated_10k', c3, p3))
            else:
                quarantine.append({'instance': nm,
                                   'family': 'krepp_estimated_10k',
                                   'reason': 'leaf_set_mismatch'})


# ---------------------------------------------------------------------------
# biotrees
# ---------------------------------------------------------------------------
def seeded_weights(seed, K):
    """Reproduce find_sol.py's abundance draw for noise type 0."""
    rs = np.random.RandomState(seed)
    p = rs.rand(K)
    p = p / p.sum()
    while p.min() < 0.01 / K:
        p = rs.rand(K)
        p = p / p.sum()
    return p


def build_biotrees(index, quarantine, census, recon):
    idir = os.path.join(DEST, 'biotrees', 'instances')
    tdir = os.path.join(DEST, 'biotrees', 'trees')
    if not os.path.isdir(idir):
        return
    for ds in sorted(os.listdir(idir)):
        for inst in sorted(os.listdir(os.path.join(idir, ds))):
            m = re.match(r'^(.*)_k(\d+)_s(\d+)_noise(\d)(?:_l(\d+))?$', inst)
            if not m:
                continue
            tf, K, seed = m.group(1), int(m.group(2)), int(m.group(3))
            ntype = int(m.group(4))
            slen = int(m.group(5)) if m.group(5) else None
            dpath = os.path.join(idir, ds, inst, 'distances.txt')
            tpath = os.path.join(tdir, tf)
            if not (os.path.exists(dpath) and os.path.exists(tpath)):
                quarantine.append({'instance': '%s/%s' % (ds, inst),
                                   'family': 'biotrees',
                                   'reason': 'missing_input_file'})
                continue
            cen = census['trees'].get('biotrees/trees/' + tf, {})
            tree = load_tree(tpath)
            lab = leaf_labels(tree)
            allq = set(lab.values())
            pub = {l: v for l, v in read_two_col(dpath).items() if is_num(v)}
            Q = sorted(allq - set(pub))
            if len(Q) != K:
                quarantine.append({'instance': '%s/%s' % (ds, inst),
                                   'family': 'biotrees',
                                   'reason': 'support_size_%d_expected_%d'
                                             % (len(Q), K)})
                continue
            prov_base = {
                'collection': 'Demixing_Data biotrees',
                'repo': 'shayesteh99/Demixing_Data',
                'commit': '5d3a1141158994af3fd6e0712fe6aaa14b52903f',
                'instance': '%s/%s' % (ds, inst), 'dataset': ds,
                'tree_file': tf, 'tree_sha256': cen.get('sha256'),
                'full_tree_leaves': cen.get('n_leaves'),
                'K_true': K, 'seed': seed, 'noise_type': ntype,
                'sequence_length': slen,
                'ground_truth': ('support recovered as leaves(full tree) minus '
                                 'rows of distances.txt; weights reproduced '
                                 'from the generator seed'),
            }
            # ---- exact family, noise type 0 only -------------------------
            if ntype == 0:
                try:
                    p = seeded_weights(seed, K)
                    adj = undirected(tree)
                    lv = {v: l for v, l in lab.items()}
                    lab2v = {l: v for v, l in lv.items()}
                    R = sorted(pub)
                    M = np.zeros((len(R), K))
                    for j, q in enumerate(Q):
                        dd = dists_from(adj, lab2v[q])
                        for i, r in enumerate(R):
                            M[i, j] = float(dd[lab2v[r]])
                    yv = np.array([float(pub[r]) for r in R])
                    w, _res, _rk, _sv = np.linalg.lstsq(M, yv, rcond=None)
                    resid = float(np.max(np.abs(M @ w - yv)))
                    scale = float(np.max(np.abs(yv)))
                    order_w = np.argsort(w)
                    order_p = np.argsort(p)
                    dev = float(np.max(np.abs(np.sort(w) - np.sort(p))))
                    gap = float(np.min(np.diff(np.sort(p)))) if K > 1 else 1.0
                    ok = (resid <= 1e-9 * max(scale, 1e-12) and dev <= 1e-7
                          and gap > 10 * max(dev, 1e-15))
                    recon.append({'instance': '%s/%s' % (ds, inst),
                                  'max_residual': resid,
                                  'max_weight_deviation': dev,
                                  'min_weight_gap': gap, 'accepted': bool(ok)})
                    if not ok:
                        raise CaseError('weight_reconstruction_unverified')
                    assign = {}
                    for a, b in zip(order_w, order_p):
                        assign[Q[a]] = float(p[b])
                    qn, tot = exact_normalised(assign)
                    case = build_exact_case(tree, qn, cap=K,
                                            name='bio_%s__%s_exact' % (ds, inst))
                except CaseError as exc:
                    quarantine.append({
                        'instance': '%s/%s' % (ds, inst),
                        'family': 'exact_biotrees', 'reason': str(exc),
                        'tree': tf,
                        'tree_negative_edges': cen.get('n_negative_edges'),
                        'tree_zero_edges': cen.get('n_zero_edges')})
                    case = None
                if case is not None:
                    prov = dict(prov_base)
                    prov.update(distance_source='exact_tree_derived',
                                exactly_explained=True,
                                weights_renormalised_exactly=str(tot) != '1',
                                weight_reconstruction_max_residual=resid,
                                weight_reconstruction_max_deviation=dev)
                    index.append(write_case('exact_biotrees', case, prov))
            # ---- published family (all noise types) -----------------------
            try:
                pmap = {q: Fraction(1, K) for q in Q}
                shell = build_exact_case(tree, pmap, cap=K,
                                         name='shell')
            except CaseError as exc:
                quarantine.append({'instance': '%s/%s' % (ds, inst),
                                   'family': 'published_biotrees',
                                   'reason': 'structure_' + str(exc),
                                   'tree': tf})
                continue
            lbl2id = {v: k for k, v in shell['leaf_labels'].items()}
            c = dict(shell)
            tag = 'noise%d' % ntype + ('_l%d' % slen if slen else '')
            c['name'] = 'bio_%s__%s_published' % (ds, inst)
            c['observed'] = {lbl2id[l]: str(Fraction(v)) for l, v in pub.items()}
            c['observed_leaves'] = len(c['observed'])
            c.pop('generating_pendant', None)
            c['generating_sources'] = [
                {'edge': s['edge'], 'x': s['x'], 'p': None,
                 'query_label': s['query_label']}
                for s in shell['source_details']]
            prov = dict(prov_base)
            prov.update(
                distance_source=('published_true_rounded' if ntype == 0
                                 else 'published_noisy_' + tag),
                exactly_explained=False,
                source_file='distances.txt', noise_tag=tag,
                note=('Source support is published; source weights are not '
                      'published for this family. x positions are those of the '
                      'real pruned attachment points.'))
            index.append(write_case('published_biotrees', c, prov))


# ---------------------------------------------------------------------------
# 16k (WoL v2 backbone, real HMP-derived compositions)
# ---------------------------------------------------------------------------
def build_16k(index, quarantine, census, skipped):
    base = os.path.join(DEST, 'sixteen_k')
    tpath = os.path.join(base, 'backbone-WoLv2.nwk')
    if not os.path.exists(tpath):
        return
    cen = census['trees'].get('sixteen_k/backbone-WoLv2.nwk', {})
    tree = load_tree(tpath)
    sim = os.path.join(base, 'simulated')
    for panel in sorted(os.listdir(sim)):
        for s in sorted(os.listdir(os.path.join(sim, panel))):
            qp = os.path.join(sim, panel, s, 'query.txt')
            raw = read_two_col(qp)
            K = len(raw)
            if K > SIXTEEN_K_MAX:
                skipped.append({'instance': '16k/%s/%s' % (panel, s), 'K': K,
                                'reason': 'K_above_materialisation_cap_%d'
                                          % SIXTEEN_K_MAX})
                continue
            try:
                qn, tot = exact_normalised({k: v for k, v in raw.items()
                                            if is_num(v)})
                case = build_exact_case(tree, qn, cap=K,
                                        name='k16k_%s_%s_exact' % (panel, s))
            except CaseError as exc:
                quarantine.append({'instance': '16k/%s/%s' % (panel, s),
                                   'family': 'exact_16k', 'reason': str(exc),
                                   'K': K})
                continue
            prov = {
                'collection': 'Demixing_Data 16k (WoL v2 backbone)',
                'repo': 'shayesteh99/Demixing_Data',
                'commit': '5d3a1141158994af3fd6e0712fe6aaa14b52903f',
                'instance': '16k/simulated/%s/%s' % (panel, s),
                'tree_file': 'backbone-WoLv2.nwk',
                'tree_sha256': cen.get('sha256'),
                'full_tree_leaves': cen.get('n_leaves'),
                'K_true': K, 'panel': panel, 'sample': s,
                'ground_truth': 'query.txt (genome, relative abundance)',
                'weights_renormalised_exactly': str(tot) != '1',
                'distance_source': 'exact_tree_derived',
                'exactly_explained': True,
            }
            index.append(write_case('exact_16k', case, prov))



# ---------------------------------------------------------------------------
# travellers' diarrhoea: real stool-sample mixtures on a real 16S backbone
# ---------------------------------------------------------------------------
def build_travelers(index, quarantine, census):
    base = os.path.join(DEST, 'travelers_diarrhea')
    tpath = os.path.join(base, 'tree_backbone.nwk')
    dpath = os.path.join(base, 'mixture_distances-all_references.tsv')
    if not (os.path.exists(tpath) and os.path.exists(dpath)):
        return
    cen = census['trees'].get('travelers_diarrhea/tree_backbone.nwk', {})
    tree = load_tree(tpath)
    lab = leaf_labels(tree)
    if cen.get('n_negative_edges') or cen.get('n_zero_edges'):
        quarantine.append({'instance': 'travelers_diarrhea',
                           'family': 'published_travelers_diarrhea',
                           'reason': 'non_positive_edge_length'})
        return
    adj = undirected(tree)
    leaves = {v for v in adj if len(adj[v]) == 1}
    R = sorted(leaves, key=lambda v: lab[v])
    vid = {r: i for i, r in enumerate(R)}
    nxt = len(R)
    for v in sorted(adj):
        if v not in vid:
            vid[v] = nxt
            nxt += 1
    out_edges = []
    done = set()
    for u in sorted(adj):
        for v, L in sorted(adj[u].items()):
            key = (min(u, v), max(u, v))
            if key in done:
                continue
            done.add(key)
            out_edges.append([vid[u], vid[v], str(L)])
    lab2id = {lab[r]: str(vid[r]) for r in R}
    with open(dpath) as fh:
        header = fh.readline().rstrip('\n').split('\t')
        samples = header[1:]
        cols = {s: {} for s in samples}
        for line in fh:
            f = line.rstrip('\n').split('\t')
            for s, v in zip(samples, f[1:]):
                if is_num(v):
                    cols[s][f[0]] = v
    meta = {}
    mp = os.path.join(base, 'metadata.tsv')
    if os.path.exists(mp):
        with open(mp) as fh:
            for line in fh:
                f = line.rstrip('\n').split('\t')
                if len(f) >= 3 and not f[0].startswith('#'):
                    meta[f[0]] = {'sample': f[1], 'group': f[2]}
    for s in samples:
        vals = cols[s]
        if set(vals) - set(lab2id):
            quarantine.append({'instance': 'travelers_diarrhea/%s' % s,
                               'family': 'published_travelers_diarrhea',
                               'reason': 'row_label_not_a_leaf'})
            continue
        case = {
            'name': 'td_sample_%s_published' % s,
            'task': 'free_pendant',
            'edges': out_edges,
            'observed': {lab2id[l]: str(Fraction(v)) for l, v in vals.items()},
            'cap': None,
            'hidden_leaves': sorted(int(lab2id[l]) for l in lab2id
                                    if l not in vals),
            'pendant': None,
            'generating_sources': None,
            'generating_count': None,
            'original_leaves': len(R),
            'observed_leaves': len(vals),
            'vertices': len(adj),
            'leaf_labels': {lab2id[l]: l for l in lab2id},
        }
        prov = {
            'collection': "Demixing_Data travellers' diarrhoea (real 16S data)",
            'repo': 'shayesteh99/Demixing_Data',
            'commit': '5d3a1141158994af3fd6e0712fe6aaa14b52903f',
            'instance': 'travelers_diarrhea/sample_%s' % s,
            'tree_file': 'tree_backbone.nwk',
            'tree_sha256': cen.get('sha256'),
            'full_tree_leaves': cen.get('n_leaves'),
            'K_true': None, 'sample_metadata': meta.get(s),
            'ground_truth': ('none: this is a real stool sample, the true '
                             'mixture composition is unknown'),
            'distance_source': 'depp_placement_estimated_16S',
            'exactly_explained': False,
            'source_file': 'mixture_distances-all_references.tsv',
            'note': ('Real biological data. Observed distances are those the '
                     'authors computed from DEPP placements of the sample '
                     'reads against the backbone; the budget (cap) is not '
                     'determined by the data and must be declared downstream.'),
        }
        index.append(write_case('published_travelers_diarrhea', case, prov))


def main():
    t0 = time.time()
    os.makedirs(CASES, exist_ok=True)
    census = json.load(open(os.path.join(ROOT, 'results',
                                         'stage287_tree_census.json')))
    index, quarantine, recon, skipped = [], [], [], []
    build_10k(index, quarantine, census)
    print('10k done: %d cases, %.0fs' % (len(index), time.time() - t0), flush=True)
    build_16k(index, quarantine, census, skipped)
    print('16k done: %d cases, %.0fs' % (len(index), time.time() - t0), flush=True)
    build_travelers(index, quarantine, census)
    print('travelers done: %d cases, %.0fs' % (len(index), time.time() - t0), flush=True)
    build_biotrees(index, quarantine, census, recon)
    print('biotrees done: %d cases, %.0fs' % (len(index), time.time() - t0), flush=True)

    # reference panels: large real trees carried as panels, no mixture
    panels = []
    for key, v in census['trees'].items():
        if 'n_leaves' not in v or v['n_leaves'] < 1000:
            continue
        panels.append({'file': 'data/stage287_instances/' + key,
                       'leaves': v['n_leaves'], 'vertices': v['n_vertices'],
                       'negative_edges': v['n_negative_edges'],
                       'zero_edges': v['n_zero_edges'],
                       'duplicate_leaf_labels': v['n_duplicate_leaf_labels'],
                       'pdd_positivity_ok': v['pdd_positivity_ok'],
                       'sha256': v['sha256']})
    panels.sort(key=lambda r: -r['leaves'])

    fam = collections.Counter(r['family'] for r in index)
    by_size = collections.defaultdict(collections.Counter)
    for r in index:
        n = r['reference_set_size']
        b = ('<=50' if n <= 50 else '51-150' if n <= 150 else
             '151-500' if n <= 500 else '501-2000' if n <= 2000 else '>2000')
        by_size[r['family']][b] += 1
    by_k = collections.defaultdict(collections.Counter)
    for r in index:
        by_k[r['family']][r['K_generating']] += 1
    qreason = collections.Counter('%s:%s' % (q['family'], q['reason'].split(':')[0])
                                 for q in quarantine)
    out = {
        'stage': 287,
        'registered_utc': datetime.now(timezone.utc).isoformat(),
        'purpose': ('Frozen real-data PDD instance bank for the Phase 2 '
                    'experiments. Built before any solver is run; no solver '
                    'call occurs in this stage.'),
        'case_format': {
            'reference': 'results/stage44_confirmation_protocol.json',
            'edges': '[u, v, length-string] with integer vertex ids and exact '
                     'rational length strings; leaves are ids 0..N-1',
            'observed': 'map from leaf id (string) to exact rational distance '
                        'string',
            'cap': 'source budget K',
            'hidden_leaves': 'empty in the bank; reference-incompleteness '
                             'variants are derived downstream',
            'pendant': 'the generating common mean pendant, supplied so that a '
                       'prescribed-pendant task needs no recomputation',
            'generating_sources': '[{edge: index into edges, p: rational '
                                  'string, x: rational string measured from '
                                  'the first endpoint of that edge}]',
            'generating_pendant': 'common mean pendant y = sum_i p_i y_i',
            'extra_keys': ['leaf_labels', 'source_details', 'portion_mapped',
                           'provenance'],
        },
        'conversion': {
            'model': ('d_r = y + sum_i p_i dist_T(r, q_i) with '
                      'y = sum_i p_i y_i; derived from the parent generator '
                      'd_r = sum_i p_i dist_Tfull(q_i, r) by pruning the K '
                      'query leaves.'),
            'arithmetic': 'Fraction(str) throughout; no float round trip',
            'validation': ('every exact-family case is accepted only if the '
                           'full-tree computation and the pruned-tree '
                           'edge-parameterised computation agree exactly for '
                           'every observed leaf'),
            'weight_normalisation': ('published decimal weights are converted '
                                     'to exact rationals and divided by their '
                                     'exact sum so that sum p = 1 holds '
                                     'exactly; the published decimals sum to 1 '
                                     'only up to float rounding'),
            'sixteen_k_materialisation_cap': SIXTEEN_K_MAX,
        },
        'repositories': json.load(
            open(os.path.join(DEST, 'download_manifest.json')))['repositories'],
        'download_manifest': 'data/stage287_instances/download_manifest.json',
        'tree_census': 'results/stage287_tree_census.json',
        'counts': {
            'n_cases': len(index),
            'by_family': dict(fam),
            'by_family_and_reference_set_size': {k: dict(v) for k, v in by_size.items()},
            'by_family_and_generating_K': {k: dict(v) for k, v in by_k.items()},
            'n_quarantined': len(quarantine),
            'quarantine_reasons': dict(qreason),
            'n_skipped_by_cap': len(skipped),
        },
        'reference_set_size_extremes': {
            'min': min((r['reference_set_size'] for r in index), default=None),
            'max': max((r['reference_set_size'] for r in index), default=None),
        },
        'weight_reconstruction': {
            'n_attempted': len(recon),
            'n_accepted': sum(1 for r in recon if r['accepted']),
            'max_residual_over_accepted': max(
                [r['max_residual'] for r in recon if r['accepted']] or [None]),
            'max_deviation_over_accepted': max(
                [r['max_weight_deviation'] for r in recon if r['accepted']] or [None]),
            'rejected': [r for r in recon if not r['accepted']][:50],
        },
        'reference_panels': panels,
        'cases': index,
        'quarantine': quarantine,
        'skipped_by_cap': skipped[:50],
    }
    with open(PROTOCOL, 'w') as fh:
        json.dump(out, fh, indent=1)
    with open(os.path.join(BANK, 'index.json'), 'w') as fh:
        json.dump({'cases': index}, fh, indent=1)
    print(json.dumps(out['counts'], indent=1))
    print('extremes', out['reference_set_size_extremes'])
    print('reconstruction', {k: v for k, v in out['weight_reconstruction'].items()
                             if k != 'rejected'})
    print('wrote', PROTOCOL, '%.0fs' % (time.time() - t0))


if __name__ == '__main__':
    main()
