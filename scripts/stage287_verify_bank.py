#!/usr/bin/env python3
"""Stage 287 (part 4) -- read-back verification of the frozen instance bank.

Reads every case file written by ``stage287_build_bank.py`` from disk, checks
its SHA256 against the protocol, and for every exact-family case recomputes

    y + sum_i p_i * min(d_T(r,u_i) + x_i L_i, d_T(r,v_i) + (1-x_i) L_i)

from the stored edge list and stored generating mixture alone, requiring exact
equality with the stored observed vector. Also checks structural invariants:
positive edge lengths, x strictly interior, distinct source edges, weights
summing to 1, observed keys being leaf ids.

Writes results/stage287_instance_bank_verification.json.
"""
import collections
import hashlib
import json
import os
import sys
from datetime import datetime, timezone
from fractions import Fraction

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROTOCOL = os.path.join(ROOT, 'results', 'stage287_instance_bank_protocol.json')


def dists_from(adj, src):
    d = {src: Fraction(0)}
    stack = [src]
    while stack:
        u = stack.pop()
        du = d[u]
        for v, L in adj[u].items():
            if v not in d:
                d[v] = du + L
                stack.append(v)
    return d


def verify_case(case):
    problems = []
    edges = [(int(u), int(v), Fraction(L)) for u, v, L in case['edges']]
    adj = collections.defaultdict(dict)
    for u, v, L in edges:
        adj[u][v] = L
        adj[v][u] = L
        if L <= 0:
            problems.append('non_positive_edge')
    deg = {v: len(adj[v]) for v in adj}
    leaves = {v for v in adj if deg[v] == 1}
    obs = {int(k): Fraction(v) for k, v in case['observed'].items()}
    hidden = set(case.get('hidden_leaves') or [])
    if not set(obs) <= leaves:
        problems.append('observed_key_not_a_leaf')
    if hidden & set(obs):
        problems.append('hidden_leaf_also_observed')
    exact_mode = 'generating_pendant' in case
    srcs = case.get('generating_sources') or []
    if (not srcs or any(s.get('p') is None for s in srcs)
            or case.get('pendant') is None):
        return problems, None, None, 'no_generating_weights'
    p = [Fraction(s['p']) for s in srcs]
    if sum(p) != 1:
        problems.append('weights_not_summing_to_one')
    if any(q <= 0 for q in p):
        problems.append('non_positive_weight')
    eidx = [int(s['edge']) for s in srcs]
    if len(set(eidx)) != len(eidx):
        problems.append('sources_share_an_edge')
    y = Fraction(case['generating_pendant'] if exact_mode else case['pendant'])
    if y <= 0:
        problems.append('non_positive_pendant')
    need = sorted({edges[i][0] for i in eidx} | {edges[i][1] for i in eidx})
    dmap = {ep: dists_from(adj, ep) for ep in need}
    pred = {r: y for r in obs}
    for s, pi in zip(srcs, p):
        u, v, L = edges[int(s['edge'])]
        x = Fraction(s['x'])
        if not (0 < x < 1):
            problems.append('x_not_strictly_interior')
        for r in pred:
            a = dmap[u][r] + x * L
            b = dmap[v][r] + (1 - x) * L
            pred[r] += pi * (a if a <= b else b)
    bad = [r for r in obs if pred[r] != obs[r]]
    dev = max((abs(float(pred[r] - obs[r])) for r in obs), default=0.0)
    if not exact_mode:
        return problems, None, dev, 'published_observations'
    return problems, len(bad), dev, None


def main():
    prot = json.load(open(PROTOCOL))
    stats = collections.Counter()
    fam_ok = collections.Counter()
    fam_fail = collections.Counter()
    failures = []
    hash_bad = []
    nmax = collections.defaultdict(int)
    nmin = collections.defaultdict(lambda: 10 ** 9)
    devs = collections.defaultdict(list)
    for rec in prot['cases']:
        path = os.path.join(ROOT, rec['path'])
        blob = open(path, 'rb').read()
        if hashlib.sha256(blob).hexdigest() != rec['sha256']:
            hash_bad.append(rec['path'])
        case = json.loads(blob)
        fam = rec['family']
        n = len(case['observed'])
        nmax[fam] = max(nmax[fam], n)
        nmin[fam] = min(nmin[fam], n)
        problems, nbad, dev, skip = verify_case(case)
        if dev is not None:
            devs[fam].append(dev)
        if skip:
            stats[skip] += 1
            if problems:
                fam_fail[fam] += 1
                failures.append({'case': rec['path'], 'problems': problems})
            else:
                fam_ok[fam] += 1
            continue
        if problems or nbad:
            fam_fail[fam] += 1
            failures.append({'case': rec['path'], 'problems': problems,
                             'n_mismatched_observations': nbad})
        else:
            fam_ok[fam] += 1
            stats['exact_verified'] += 1
    out = {
        'stage': 287, 'part': 'verification',
        'verified_utc': datetime.now(timezone.utc).isoformat(),
        'n_cases_read_back': len(prot['cases']),
        'n_sha256_mismatches': len(hash_bad),
        'sha256_mismatches': hash_bad[:20],
        'n_exact_cases_recomputed_and_matching': stats['exact_verified'],
        'n_structure_only_cases': stats['no_generating_weights'],
        'n_published_observation_cases': stats['published_observations'],
        'max_abs_model_minus_observed_by_family': {
            f: max(v) for f, v in devs.items()},
        'median_abs_model_minus_observed_by_family': {
            f: sorted(v)[len(v) // 2] for f, v in devs.items()},
        'by_family_passing': dict(fam_ok),
        'by_family_failing': dict(fam_fail),
        'failures': failures[:50],
        'reference_set_size_by_family': {
            f: {'min': nmin[f], 'max': nmax[f]} for f in nmax},
    }
    with open(os.path.join(ROOT, 'results',
                           'stage287_instance_bank_verification.json'), 'w') as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps({k: v for k, v in out.items() if k != 'failures'}, indent=1))
    if failures:
        print('FAILURES:', json.dumps(failures[:5], indent=1))


if __name__ == '__main__':
    main()
