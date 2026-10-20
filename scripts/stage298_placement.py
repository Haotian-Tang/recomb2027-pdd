#!/usr/bin/env python3
"""Stage 298 step 1: convert the krepp PLACEMENT-distance instances into
bank-conformant PDD cases - or refuse them with a machine-readable reason.

The gap audit records that the 60 krepp placement-distance instances in
bo1929/shared.krepp were fetched (they sit in the Demixing_Data instance
directories as distances_from_placement.txt) but never converted to bank
cases.  This script applies exactly the stage-287 pipeline: the frozen
`build_exact_case` normal-form conversion, the frozen
`attach_published_observed` observation attachment, and the frozen
`stage287_verify_bank.verify_case` conformance check, plus an explicit
positivity check on every observed value.  Nothing pre-existing is modified:
cases are written under data/stage298_placement/.

ONE DOCUMENTED DIFFERENCE from the stage-287 krepp_estimated_10k branch.  That
branch required set(published labels) == set(reference leaf labels).  The
placement files report a distance for EVERY leaf of the tree, including the K
query taxa that the model turns into interior sources and removes from the
reference tree, so the published label set is a strict SUPERSET of the
reference leaf set.  The rule applied here is therefore containment -
every reference leaf must carry a placement distance - with the K
query-taxon entries dropped, which is what the model requires (observations
live at reference leaves only).  The count of dropped entries is recorded per
case.
"""
import hashlib
import json
import os
import re
import sys
import time
from datetime import datetime, timezone
from fractions import Fraction

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
from stage287_pdd_cases import (build_exact_case, CaseError,          # noqa: E402
                                attach_published_observed)
from stage287_build_bank import (load_tree, read_two_col, is_num,     # noqa: E402
                                 exact_normalised, sha256_bytes)
from stage287_verify_bank import verify_case                          # noqa: E402

SRC = os.path.join(ROOT, 'data', 'stage287_instances', 'demixing_10k')
DEST = os.path.join(ROOT, 'data', 'stage298_placement')
OUT = os.path.join(ROOT, 'results', 'stage298')
MANIFEST = os.path.join(OUT, 'stage298_placement_manifest.json')
FILE = 'distances_from_placement.txt'
# the registered bank schema (stage287_instance_bank_protocol.json)
REQUIRED = ['name', 'task', 'edges', 'observed', 'cap', 'hidden_leaves',
            'pendant', 'generating_sources', 'generating_count', 'vertices',
            'original_leaves', 'observed_leaves']


def census_tree(tid):
    prot = json.load(open(os.path.join(ROOT, 'results',
                                       'stage287_tree_census.json')))
    key = 'demixing_10k/trees/tree_%s.tre' % tid
    return prot['trees'].get(key, {})


def main():
    assert not os.path.exists(MANIFEST), 'Never overwrite frozen evidence'
    os.makedirs(os.path.join(DEST, 'cases'), exist_ok=True)
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    names = sorted(d for d in os.listdir(os.path.join(SRC, 'instances'))
                   if os.path.exists(os.path.join(SRC, 'instances', d, FILE)))
    converted, refused, source_files = [], [], []
    for nm in names:
        m = re.match(r'tree_(\d+)_k(\d+)_s(\d+)_(uniform|exp)$', nm)
        if not m:
            refused.append({'instance': nm, 'reason': 'unparsable_instance_name'})
            continue
        tid, K, rep, dist = m.group(1), int(m.group(2)), int(m.group(3)), m.group(4)
        idir = os.path.join(SRC, 'instances', nm)
        tpath = os.path.join(SRC, 'trees', 'tree_%s.tre' % tid)
        ppath = os.path.join(idir, FILE)
        qpath = os.path.join(idir, 'queries.txt')
        if not (os.path.exists(tpath) and os.path.exists(qpath)):
            refused.append({'instance': nm, 'reason': 'missing_input_file'})
            continue
        source_files.append({'path': os.path.relpath(ppath, ROOT),
                             'sha256': sha256_bytes(open(ppath, 'rb').read())})
        cen = census_tree(tid)
        tree = load_tree(tpath)
        raw = {k: v for k, v in read_two_col(qpath).items() if is_num(v)}
        try:
            qn, tot = exact_normalised(raw)
            case = build_exact_case(tree, qn, cap=K, name='plc_%s_exact' % nm)
        except CaseError as exc:
            refused.append({'instance': nm, 'reason': 'exact_case:%s' % exc})
            continue
        lbl2id = {v: k for k, v in case['leaf_labels'].items()}
        pub_raw = read_two_col(ppath)
        non_numeric = sorted(l for l, v in pub_raw.items() if not is_num(v))
        pub = {l: v for l, v in pub_raw.items() if is_num(v)}
        missing_refs = sorted(l for l in lbl2id if l not in pub)
        if missing_refs:
            refused.append({'instance': nm, 'reason': 'reference_leaf_without_placement_distance',
                            'n_missing': len(missing_refs),
                            'examples': missing_refs[:5]})
            continue
        extra = sorted(l for l in pub if l not in lbl2id)
        non_positive = sorted(l for l in lbl2id if Fraction(pub[l]) <= 0)
        if non_positive:
            refused.append({'instance': nm, 'reason': 'non_positive_observed_distance',
                            'n_non_positive': len(non_positive),
                            'examples': non_positive[:5]})
            continue
        c = attach_published_observed(case, {l: pub[l] for l in lbl2id}, lbl2id)
        c['name'] = 'plc_%s_krepp_placement' % nm
        c['observed_leaves'] = len(c['observed'])
        c['hidden_leaves'] = []
        missing_schema = [k for k in REQUIRED if k not in c]
        if missing_schema:
            refused.append({'instance': nm, 'reason': 'schema_keys_missing:%s'
                            % ','.join(missing_schema)})
            continue
        problems, nbad, dev, skip = verify_case(c)
        if problems:
            refused.append({'instance': nm, 'reason': 'verify_case:%s'
                            % ','.join(sorted(set(problems)))})
            continue
        if skip != 'published_observations':
            refused.append({'instance': nm, 'reason': 'unexpected_verify_mode:%s' % skip})
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
            'distance_source': 'krepp_placement_estimated',
            'exactly_explained': False,
            'source_file': FILE,
            'source_file_sha256': source_files[-1]['sha256'],
            'n_published_entries': len(pub_raw),
            'n_non_numeric_entries': len(non_numeric),
            'n_entries_dropped_as_query_taxa': len(extra),
            'dropped_entries_are_exactly_the_queries': sorted(extra) == sorted(raw),
            'containment_rule': ('placement files cover every leaf of the tree; the rule '
                                 'applied is that every reference leaf carries a distance, '
                                 'with the K query-taxon entries dropped'),
            'note': ('Distances to the mixture sample estimated by krepp placement '
                     '(distances_from_placement.txt), paired with the positive-length '
                     'reference tree, exactly as krepp_estimated_10k is.'),
            'observed_max_deviation_from_true_configuration': dev,
        }
        blob = json.dumps(c, sort_keys=True).encode()
        rel = os.path.join('data', 'stage298_placement', 'cases', c['name'] + '.json')
        with open(os.path.join(ROOT, rel), 'wb') as fh:
            fh.write(blob)
        converted.append({
            'name': c['name'], 'family': 'krepp_placement_10k', 'path': rel,
            'sha256': sha256_bytes(blob), 'instance': nm,
            'K_generating': K, 'reference_set_size': len(c['observed']),
            'n_edges': len(c['edges']), 'exact_partner_instance': nm,
            'provenance': prov})
    manifest = {
        'stage': 298, 'step': 'placement_conversion',
        'registered_utc': datetime.now(timezone.utc).isoformat(),
        'source_file_name': FILE,
        'n_instances_with_placement_file': len(names),
        'n_converted': len(converted), 'n_refused': len(refused),
        'rule': ('frozen build_exact_case normal form + frozen '
                 'attach_published_observed + frozen verify_case, with '
                 'containment of the reference leaf set in the published label '
                 'set and strict positivity of every observed value'),
        'cases': converted, 'refused': refused,
        'source_files': source_files,
        'elapsed_seconds': round(time.time() - t0, 2),
    }
    with open(MANIFEST, 'w') as fh:
        json.dump(manifest, fh, indent=1, sort_keys=True)
    print('converted %d, refused %d, %.1fs' % (len(converted), len(refused),
                                               time.time() - t0))
    for r in refused[:10]:
        print('  refused', r)
    if converted:
        ks = sorted({c['K_generating'] for c in converted})
        print('  K', ks, 'N range', min(c['reference_set_size'] for c in converted),
              max(c['reference_set_size'] for c in converted),
              '|E| range', min(c['n_edges'] for c in converted),
              max(c['n_edges'] for c in converted))
        print('  dropped-entries-are-queries:',
              sum(1 for c in converted
                  if c['provenance']['dropped_entries_are_exactly_the_queries']),
              '/', len(converted))


if __name__ == '__main__':
    main()

