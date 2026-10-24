#!/usr/bin/env python3
"""Audit the two shared stool clade masks from frozen records; no solver runs.

Inputs: results/stage299/stage299_td_runs.jsonl AND stage299_protocol.json
(both in the bulk archive), bulk_manifest.json, the bank index, and all 29
hash-registered travellers' case JSONs. The protocol contains explicit removed
IDs that the compact run log omits. Runs are joined to that plan by their full
identity; equal removal counts alone are not evidence of equal masks.

Run from the repository root after fetching the documented bulk archive:
    python3 -B scripts/stage311_shared_mask_audit.py
Or use --bulk-root PATH for an existing copy of those two frozen bulk files.
--output PATH redirects the result for isolated audits. No inputs are changed.
"""
import argparse
import hashlib
import json
import statistics
import sys
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LEVEL = '0.75'
PROTOCOL = 'results/stage299/stage299_protocol.json'
RUNS = 'results/stage299/stage299_td_runs.jsonl'
INDEX = 'data/stage287_instances/bank/index.json'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as fh:
        for block in iter(lambda: fh.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def object_sha(value):
    return hashlib.sha256(json.dumps(value, separators=(',', ':'),
                                     sort_keys=True).encode()).hexdigest()


def identity(row):
    return (row['arm'], row['name'], row['level'], row['regime'],
            row['replicate'])


def unique_rows(rows, description):
    out = {}
    for row in rows:
        key = identity(row)
        require(key not in out, 'duplicate %s identity: %r' % (description, key))
        out[key] = row
    return out


def checked_bulk(root, relative, manifest):
    path = root / relative
    require(path.is_file(), 'missing bulk input %s; fetch the README bulk archive '
            'or pass --bulk-root' % path)
    entry = manifest[relative]
    require(path.stat().st_size == entry['bytes'], 'size mismatch: ' + relative)
    require(sha256(path) == entry['sha256'], 'SHA-256 mismatch: ' + relative)
    return path


def audit(bulk_root):
    manifest_path = ROOT / 'bulk_manifest.json'
    manifest = {e['path']: e for e in json.loads(manifest_path.read_text())['members']}
    protocol_path = checked_bulk(bulk_root, PROTOCOL, manifest)
    runs_path = checked_bulk(bulk_root, RUNS, manifest)
    protocol = json.loads(protocol_path.read_text())
    with runs_path.open() as fh:
        runs = unique_rows((json.loads(line) for line in fh if line.strip()), 'run')
    plan = unique_rows((r for r in protocol['leafsets'] if r['arm'] == 'travelers'),
                       'registered plan')
    require(set(runs) == set(plan), 'run identities differ from registered stool plan')

    index_path = ROOT / INDEX
    bank = {c['name']: c for c in json.loads(index_path.read_text())['cases']}
    cases = [c for c in protocol['cases'] if c['arm'] == 'travelers']
    names = {c['name'] for c in cases}
    require(len(cases) == len(names) == 29, 'expected exactly 29 registered stool cases')
    require({k[1] for k in plan} == names, 'plan sample set differs from cases')
    common = None
    case_hashes = {}
    for c in cases:
        name = c['name']
        require(name in bank, 'case missing from bank index: ' + name)
        require((c['path'], c['sha256']) == (bank[name]['path'], bank[name]['sha256']),
                'protocol/index case mismatch: ' + name)
        path = ROOT / c['path']
        require(sha256(path) == c['sha256'], 'case SHA-256 mismatch: ' + name)
        raw = json.loads(path.read_text())
        require(raw['name'] == name, 'case name mismatch')
        panel = tuple(sorted(int(r) for r in raw['observed']))
        require(len(panel) == len(raw['observed']) == 10575, 'unexpected observed panel')
        tree_id = raw['provenance']['tree_sha256']
        geometry = (raw['edges'], panel, raw['leaf_labels'], tree_id)
        if common is None:
            common = geometry
        else:
            require(geometry == common, 'backbone, metric, reference labels or panel differ: ' + name)
        require(c['reference_set_size'] == len(panel) and c['n_edges'] == len(raw['edges']),
                'case metadata size mismatch: ' + name)
        case_hashes[c['path']] = c['sha256']
    edges, panel, leaf_labels, tree_id = common
    panel_set = set(panel)
    require(len(edges) == 21147, 'unexpected backbone edge count')
    degree = defaultdict(int)
    for u, v, length in edges:
        require(u != v and F(length) > 0, 'invalid backbone edge')
        degree[u] += 1
        degree[v] += 1
    require({v for v, deg in degree.items() if deg == 1} == panel_set,
            'observed panel is not the complete backbone leaf set')

    for key, row in runs.items():
        p = plan[key]
        require(row['arm'] == 'travelers', 'non-stool run in stool log')
        require(row['result'].get('status') == 'complete', 'incomplete run: %r' % (key,))
        removed = p['removed']
        require(removed == sorted(set(removed)) and set(removed) <= panel_set,
                'registered mask is not a unique sorted panel subset: %r' % (key,))
        require(len(removed) == p['n_removed'] and
                p['n_retained'] == len(panel) - len(removed), 'plan count mismatch')
        for field in ('n_removed', 'n_retained', 'n_edges', 'reference_set_size',
                      'stage298_theta1', 'stage298_argmin'):
            require(row[field] == p[field], 'run/plan %s mismatch: %r' % (field, key))
        detail = {k: v for k, v in p['detail'].items() if k != 'source'}
        require(row['detail'] == detail, 'run/plan detail mismatch: %r' % (key,))
        require(row['result']['n_leaves'] == p['n_retained'] and
                row['result']['n_edges'] == len(edges), 'result geometry counts mismatch')
        require(F(row['result']['theta1']) >= 0, 'negative theta: %r' % (key,))

    base = {}
    groups = {g: {} for g in ('clade', 'random')}
    for name in sorted(names):
        baseline = runs.get(('travelers', name, 'full', 'none', 0))
        require(baseline is not None, 'missing full-panel baseline: ' + name)
        base[name] = F(baseline['result']['theta1'])
        require(base[name] > 0, 'nonpositive full-panel theta: ' + name)
        require(base[name] == F(baseline['stage298_theta1']), 'full baseline disagrees with Stage298')
        for g in groups:
            keys = [k for k in runs if k[1] == name and k[2] == LEVEL and k[3] == g]
            require(len(keys) == 2 and {k[4] for k in keys} == {0, 1},
                    'expected exactly two %s draws for %s' % (g, name))
            groups[g][name] = [runs[k] for k in sorted(keys)]

    def loss(row):
        value = 1 - F(row['result']['theta1']) / base[row['name']]
        require(0 <= value <= 1, 'invalid removal loss')
        return value

    masks = defaultdict(list)
    for rs in groups['clade'].values():
        for row in rs:
            masks[tuple(plan[identity(row)]['removed'])].append(row)
    require(len(masks) == 2, 'expected exactly two actual clade deletion sets')
    ordered = sorted(masks, key=lambda mask: (len(mask), mask))
    smaller, larger = map(set, ordered)
    require(smaller < larger, 'the two explicit deletion sets are not strictly nested')
    per_set = []
    for mask in ordered:
        rs = masks[mask]
        require(len(rs) == 29 and {r['name'] for r in rs} == names,
                'a clade mask is not represented once in each sample')
        require(len({(r['detail']['clade_vertex'], r['detail']['complement']) for r in rs}) == 1,
                'same mask has inconsistent clade provenance')
        value = statistics.median([loss(r) for r in rs])
        per_set.append(dict(n_removed=len(mask), n_results=len(rs),
                            fraction_of_panel=float(F(len(mask), len(panel))),
                            median_loss=float(value), median_loss_exact=str(value),
                            removed_ids_sha256=object_sha(mask),
                            clade_detail=rs[0]['detail'],
                            replicates=sorted({r['replicate'] for r in rs})))
    sample_losses = []
    for name in sorted(names):
        cs = [loss(r) for r in groups['clade'][name]]
        rs = [loss(r) for r in groups['random'][name]]
        sample_losses.append(dict(name=name, clade_losses_exact=[str(x) for x in cs],
                                  random_losses_exact=[str(x) for x in rs],
                                  both_clade_exceed_both_random=min(cs) > max(rs)))
    clade_losses = [loss(r) for rs in groups['clade'].values() for r in rs]
    random_losses = [loss(r) for rs in groups['random'].values() for r in rs]
    return dict(stage=311, kind='shared_mask_audit', level=LEVEL, panel_size=len(panel),
                n_samples=len(names), n_clade_results=len(clade_losses),
                n_clade_deletion_sets=len(masks), deletion_sets=per_set,
                samples_where_both_clade_exceed_both_random=sum(
                    x['both_clade_exceed_both_random'] for x in sample_losses),
                min_clade_loss=float(min(clade_losses)), max_random_loss=float(max(random_losses)),
                min_clade_loss_exact=str(min(clade_losses)),
                max_random_loss_exact=str(max(random_losses)), per_sample=sample_losses,
                mask_nesting=dict(strict_subset=True, smaller_size=len(smaller),
                                  larger_size=len(larger), intersection_size=len(smaller & larger),
                                  larger_minus_smaller=len(larger - smaller)),
                backbone=dict(n_edges=len(edges), n_observed=len(panel), tree_source_sha256=tree_id,
                              exact_ordered_edges_sha256=object_sha(edges),
                              observed_ids_sha256=object_sha(panel),
                              leaf_labels_sha256=object_sha(leaf_labels),
                              identical_across_samples=True),
                provenance=dict(inputs_sha256={
                    PROTOCOL: manifest[PROTOCOL]['sha256'], RUNS: manifest[RUNS]['sha256'],
                    'bulk_manifest.json': sha256(manifest_path), INDEX: sha256(index_path),
                    **case_hashes}, matched_run_plan_identities=len(runs),
                    mask_identity_encoding='SHA-256 of compact JSON of sorted original reference IDs',
                    method='Explicit registered masks joined to complete frozen runs by '
                           '(arm,name,level,regime,replicate), with matching counts and detail; '
                           'case files prove common exact tree, reference labels and panel.',
                    limitation='The compact run log omits removed IDs; identity is established '
                               'through its frozen registered plan and the producer that consumes it. '
                               'No solver or removal-plan generation is rerun.'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bulk-root', type=Path, default=ROOT)
    parser.add_argument('--output', type=Path, default=ROOT / 'results/stage311/stage311_shared_mask_audit.json')
    args = parser.parse_args()
    try:
        out = audit(args.bulk_root)
    except (ValueError, KeyError, OSError, TypeError, json.JSONDecodeError) as exc:
        print('audit failed: %s' % exc, file=sys.stderr)
        return 1
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=1, sort_keys=True) + '\n')
    print('%d clade results = %d samples x %d deletion sets on one %d-leaf backbone'
          % (out['n_clade_results'], out['n_samples'], out['n_clade_deletion_sets'], out['panel_size']))
    for s in out['deletion_sets']:
        print('  set of %5d references (%.3f%% of panel), %2d results, median Theta_1 loss %.1f%%'
              % (s['n_removed'], 100 * s['fraction_of_panel'], s['n_results'], 100 * s['median_loss']))
    print('  both clade losses exceed both random losses in %d of %d samples'
          % (out['samples_where_both_clade_exceed_both_random'], out['n_samples']))
    print('  smallest clade loss %.1f%%, largest random loss %.1f%%'
          % (100 * out['min_clade_loss'], 100 * out['max_random_loss']))
    print('  exact registered deletion sets are strictly nested: True')
    return 0


if __name__ == '__main__':
    sys.exit(main())
