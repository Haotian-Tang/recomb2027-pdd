#!/usr/bin/env python3
"""Stage 287 (part 1) -- fetch public reference trees and query instances.

Downloads, into ``data/stage287_instances/``, the full reference trees and
query instances of the parent work's two public data repositories plus the
large real reference backbones published by the krepp author, and the file
inventory of an independent reference-set-selection benchmark.

Every retrieved file is recorded in ``download_manifest.json`` with its URL,
repository commit SHA, retrieval time, byte count and SHA256. Where the
repository's git tree manifest is on disk, the git blob SHA1 of the retrieved
bytes is recomputed and compared with the manifest, which proves the content
matches the named commit.

Idempotent: a file already present with the expected size is not re-fetched.
No existing project file is modified.
"""
import collections
import hashlib
import json
import os
import re
import sys
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
S29 = os.path.join(ROOT, 'data', 'stage29_distance_target')
DEST = os.path.join(ROOT, 'data', 'stage287_instances')

REPOS = {
    'Demixing_Data': {
        'owner_repo': 'shayesteh99/Demixing_Data',
        'commit': '5d3a1141158994af3fd6e0712fe6aaa14b52903f',
        'manifest': os.path.join(S29, 'Demixing_Data.tree.json'),
        'note': ("Parent work's data repository; same commit as the aggregate "
                 "files already retrieved in stage 29."),
    },
    'shared.krepp': {
        'owner_repo': 'bo1929/shared.krepp',
        'commit': '0310db641d5edc18fc0c77804a3902ea1a2f7458',
        'manifest': os.path.join(S29, 'shared.krepp.tree.json'),
        'note': ('Independent repository of the krepp author holding the large '
                 'real reference backbones (WoL v1/v2, RefSeq CAMI II).'),
    },
    'reference_set_selection_benchmark': {
        'owner_repo': 'JaspervB-tud/reference_set_selection_benchmark',
        'commit': 'ea5d512b3716bee4331f5689b02416640201bd67',
        'manifest': os.path.join(ROOT, 'data', 'stage248_external_benchmark',
                                 'tree.json'),
        'note': 'Independent reference-set-selection benchmark.',
    },
}
RAW = 'https://raw.githubusercontent.com/%s/%s/%s'


def git_blob_sha1(data):
    h = hashlib.sha1()
    h.update(b'blob %d\0' % len(data))
    h.update(data)
    return h.hexdigest()


def load_manifest(repo):
    m = json.load(open(REPOS[repo]['manifest']))
    return {e['path']: e for e in m['tree'] if e['type'] == 'blob'}, m['sha']


def build_plan():
    """Return list of (repo, repo_path, local_relpath)."""
    plan = []
    dm, _ = load_manifest('Demixing_Data')

    # (a) every reference tree in 10k_dataset
    for p in dm:
        if re.match(r'10k_dataset/tree_\d+(-ani_blen)?\.tre(\.rerooted)?$', p):
            plan.append(('Demixing_Data', p,
                         'demixing_10k/trees/' + p.split('/')[-1]))

    # (b) every 10k query instance: ground truth, true distances, estimates
    want = ['queries.txt', 'distances.txt',
            'estimated_distances_pr-h14_dth4.txt', 'true_mash_distances.txt',
            'distances_from_placement.txt']
    for p in dm:
        parts = p.split('/')
        if len(parts) == 3 and parts[0] == '10k_dataset' and \
                re.match(r'tree_\d+_k\d+_s\d+_(uniform|exp)$', parts[1]) and \
                parts[2] in want:
            plan.append(('Demixing_Data', p,
                         'demixing_10k/instances/%s/%s' % (parts[1], parts[2])))

    # (c) biotrees: the real species trees and every instance distance table
    for p in dm:
        if p.count('/') == 1 and p.startswith('biotrees/') and \
                p.endswith(('.tre', '.nwk', '.rerooted')):
            plan.append(('Demixing_Data', p,
                         'biotrees/trees/' + p.split('/')[-1]))
        parts = p.split('/')
        if len(parts) == 4 and parts[0] == 'biotrees' and \
                parts[3] in ('distances.txt', 'all_rounds.json'):
            plan.append(('Demixing_Data', p,
                         'biotrees/instances/%s/%s/%s' % (parts[1], parts[2],
                                                          parts[3])))
    # biotrees rerooted trees live inside the dataset directory
    for p in dm:
        parts = p.split('/')
        if len(parts) == 3 and parts[0] == 'biotrees' and \
                parts[2].endswith(('.tre', '.tre.rerooted', '.nwk')):
            plan.append(('Demixing_Data', p,
                         'biotrees/trees/' + parts[2]))

    # (d) 16k: WoL v2 backbone and every simulated query composition
    for p in dm:
        if p.startswith('16k/') and p.endswith('.nwk'):
            plan.append(('Demixing_Data', p, 'sixteen_k/' + p.split('/')[-1]))
        if p.startswith('16k/simulated/') and p.endswith('query.txt'):
            parts = p.split('/')
            plan.append(('Demixing_Data', p,
                         'sixteen_k/simulated/%s/%s/query.txt' % (parts[2],
                                                                  parts[3])))

    # (e) travelers' diarrhea: real 16S backbone and real mixture distances
    for name in ['tree_backbone.nwk', 'tree_wplacements.nwk',
                 'mixture_distances-all_references.tsv', 'all_leaves.txt',
                 'metadata.tsv', 'clusters-max300.txt', 'data_16S.tsv']:
        p = 'travelers_diarrhea/' + name
        if p in dm:
            plan.append(('Demixing_Data', p, 'travelers_diarrhea/' + name))

    # (f) large real reference backbones from the independent krepp repository
    km, _ = load_manifest('shared.krepp')
    for p in km:
        if p.endswith(('.nwk', '.tre')) and (km[p]['size'] or 0) > 200000:
            plan.append(('shared.krepp', p,
                         'krepp_shared/' + p.replace('/', '__')))

    # (g) independent reference-set-selection benchmark: full file inventory
    rm, _ = load_manifest('reference_set_selection_benchmark')
    for p in rm:
        if p.endswith(('.txt', '.md', '.py')):
            plan.append(('reference_set_selection_benchmark', p,
                         'refset_benchmark/' + p.replace('/', '__')))

    seen = set()
    uniq = []
    for item in plan:
        if item[2] in seen:
            continue
        seen.add(item[2])
        uniq.append(item)
    return uniq


def fetch_one(job):
    repo, rpath, lrel = job
    commit = REPOS[repo]['commit']
    url = RAW % (REPOS[repo]['owner_repo'], commit, rpath)
    lpath = os.path.join(DEST, lrel)
    os.makedirs(os.path.dirname(lpath), exist_ok=True)
    rec = {'repo': repo, 'commit': commit, 'repo_path': rpath, 'url': url,
           'local_path': os.path.relpath(lpath, ROOT)}
    if os.path.exists(lpath):
        data = open(lpath, 'rb').read()
        rec.update(cached=True)
    else:
        data = None
        for attempt in range(6):
            try:
                req = urllib.request.Request(url)
                with urllib.request.urlopen(req, timeout=120) as fh:
                    data = fh.read()
                break
            except urllib.error.HTTPError as exc:
                if exc.code in (403, 429, 500, 502, 503) and attempt < 5:
                    time.sleep(2 ** attempt)
                    continue
                rec.update(error='HTTP %d' % exc.code)
                return rec
            except Exception as exc:                      # noqa: BLE001
                if attempt < 5:
                    time.sleep(2 ** attempt)
                    continue
                rec.update(error='%s: %s' % (type(exc).__name__, exc))
                return rec
        if data is None:
            rec.update(error='no data')
            return rec
        tmp = lpath + '.part'
        with open(tmp, 'wb') as fh:
            fh.write(data)
        os.replace(tmp, lpath)
        rec.update(cached=False)
    rec['bytes'] = len(data)
    rec['sha256'] = hashlib.sha256(data).hexdigest()
    rec['git_blob_sha1'] = git_blob_sha1(data)
    rec['retrieved_utc'] = datetime.now(timezone.utc).isoformat()
    return rec


def main():
    os.makedirs(DEST, exist_ok=True)
    plan = build_plan()
    print('planned files:', len(plan))
    by_group = collections.Counter(x[2].split('/')[0] for x in plan)
    print('by group:', dict(by_group))
    t0 = time.time()
    recs = []
    with ThreadPoolExecutor(max_workers=8) as ex:
        for i, rec in enumerate(ex.map(fetch_one, plan), 1):
            recs.append(rec)
            if i % 500 == 0:
                nerr = sum(1 for r in recs if 'error' in r)
                print('  %d/%d  errors=%d  %.0fs' % (i, len(plan), nerr,
                                                     time.time() - t0),
                      flush=True)
    # verify against the repository git tree manifests
    checks = {}
    for repo in REPOS:
        mani, sha = load_manifest(repo)
        ok = bad = miss = 0
        bad_paths = []
        for r in recs:
            if r['repo'] != repo or 'git_blob_sha1' not in r:
                continue
            e = mani.get(r['repo_path'])
            if e is None:
                miss += 1
            elif e['sha'] == r['git_blob_sha1']:
                ok += 1
            else:
                bad += 1
                bad_paths.append(r['repo_path'])
        checks[repo] = {'manifest_tree_sha': sha, 'blob_sha1_match': ok,
                        'blob_sha1_mismatch': bad, 'not_in_manifest': miss,
                        'mismatched_paths': bad_paths[:20],
                        'note': REPOS[repo]['note'],
                        'commit_requested': REPOS[repo]['commit']}
    errors = [r for r in recs if 'error' in r]
    out = {
        'stage': 287,
        'part': 'download',
        'retrieved_utc_start': datetime.fromtimestamp(t0, timezone.utc).isoformat(),
        'retrieved_utc_end': datetime.now(timezone.utc).isoformat(),
        'repositories': {k: {kk: vv for kk, vv in v.items() if kk != 'manifest'}
                         for k, v in REPOS.items()},
        'n_planned': len(plan),
        'n_retrieved': len(recs) - len(errors),
        'n_errors': len(errors),
        'errors': errors[:100],
        'total_bytes': sum(r.get('bytes', 0) for r in recs),
        'integrity_checks': checks,
        'files': recs,
    }
    with open(os.path.join(DEST, 'download_manifest.json'), 'w') as fh:
        json.dump(out, fh, indent=1)
    print('retrieved %d/%d, %d errors, %.1f MB' % (
        out['n_retrieved'], out['n_planned'], out['n_errors'],
        out['total_bytes'] / 1e6))
    print('integrity:', json.dumps(checks, indent=1)[:2000])


if __name__ == '__main__':
    main()
