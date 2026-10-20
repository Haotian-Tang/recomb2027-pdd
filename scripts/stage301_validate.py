#!/usr/bin/env python3
"""Stage 301 step 1: validate the transported localisation.

Arm A  replication   - every stage-296 cross/probe run with a recorded Theta_k
                       is recomputed with the mass screen at eps = Theta_k.
                       The localised value must equal the recorded one EXACTLY.
Arm A2 one-sidedness - the same cases at eps = Theta_k / 2, where the routine
                       must never return a value <= eps (it may only certify
                       Theta_k > eps or refute).
Arm B  independent   - fresh exhaustive sweep vs localised on small cases over a
                       ladder of eps, both computed here.
Arm C  visit count   - with the screen disabled the enumerator must visit
                       exactly C(|E|, k) supports.

Every surviving witness is re-verified with the frozen independent verifier.
"""
import json
import os
import sys
import time
from fractions import Fraction as F
from itertools import combinations
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import stage288_common as common            # noqa: E402
import stage296_common as s296              # noqa: E402
import stage301_common as s301              # noqa: E402

OUT = s301.OUT
PROTOCOL = OUT / 'stage301_protocol.json'
WORKERS = 40
ADDRESS = 8 * 1024 ** 3
K2_EDGES = 160
K2_BIG = 40


def case_for(row, index):
    entry = next(e for e in index if e['name'] == row['name'])
    case = common.load_case(entry)
    removed = set(row['removed'])
    obs = {r: d for r, d in case['observed'].items() if r not in removed}
    assert len(obs) == row['n_retained'], (len(obs), row['n_retained'])
    return case, obs


def worker(spec):
    edges, obs, k, eps = spec['edges'], spec['obs'], spec['k'], spec['eps']
    mode = spec['mode']
    leaves = sorted(obs)
    dvals = tuple(obs[r] for r in leaves)
    t0 = time.perf_counter()
    if mode == 'exhaustive':
        Dv = s296.leaf_vectors(edges, leaves)
        out = s296.theta_min_over_supports(Dv, edges, dvals, k)
        out['seconds'] = time.perf_counter() - t0
        out['mode'] = mode
        if out.get('support') is not None:
            out['independent_deviation'] = s296.independent_deviation(
                edges, out['support'], out['z'], k, obs, leaves)
        return out
    if mode == 'count':
        view = s301.backbone_view(edges, obs)
        it, n, stats = s301.candidate_plan(view, eps, k, len(edges), hard_cap=None,
                                           disabled=True)
        visited = sum(1 for _ in it)
        return dict(status='complete', mode=mode, visited=visited,
                    exhaustive=comb(len(edges), k), screen=stats,
                    seconds=time.perf_counter() - t0)
    Dv = s296.leaf_vectors(edges, leaves)
    out = s301.theta_localised(edges, obs, k, eps, Dv=Dv, dvals=dvals, leaves=leaves)
    out['seconds'] = time.perf_counter() - t0
    out['mode'] = mode
    out['exhaustive'] = comb(len(edges), k)
    if out.get('support') is not None:
        out['independent_deviation'] = s296.independent_deviation(
            edges, out['support'], out['z'], k, obs, leaves)
    return out


def build_jobs(limit_a, limit_b):
    index = common.bank_index()
    rows = [json.loads(l) for l in
            (ROOT / 'results/stage296/stage296_cross_runs.jsonl').read_text().splitlines()]
    rows += [json.loads(l) for l in
             (ROOT / 'results/stage296/stage296_probe_runs.jsonl').read_text().splitlines()]
    jobs = []
    seen = 0
    # All k = 1 rows; k >= 2 rows only where the fallback exhaustive sweep is
    # affordable (n_edges <= K2_EDGES) plus a fixed random sample of larger ones,
    # so that a vacuous screen cannot turn one validation into a 1800 s timeout.
    import random as _r
    big = [i for i, row in enumerate(rows)
           if '2' in (row['result'].get('thetas') or {}) and row['n_edges'] > K2_EDGES]
    keep_big = set(_r.Random(301).sample(big, min(K2_BIG, len(big))))
    for i, row in enumerate(rows):
        th = row['result'].get('thetas') or {}
        if not th:
            continue
        if limit_a is not None and seen >= limit_a:
            break
        seen += 1
        case, obs = case_for(row, index)
        for ks, val in sorted(th.items()):
            k, theta = int(ks), F(val)
            if k >= 2 and row['n_edges'] > K2_EDGES and i not in keep_big:
                continue
            for arm, eps in (('A', theta), ('A2', theta / 2)):
                jobs.append(dict(
                    meta=dict(arm=arm, name=row['name'], family=row['family'],
                              K=row['K'], n_edges=row['n_edges'],
                              n_retained=row['n_retained'], regime=row['regime'],
                              level=row['level'], replicate=row['replicate'],
                              k=k, eps=str(eps), recorded_theta=str(theta)),
                    spec=dict(edges=case['edges'], obs=obs, k=k, eps=eps,
                              mode='localised')))
    # Arm B/C on small cases: fresh exhaustive and the disabled-screen count.
    small = [e for e in index if e['family'] == 'exact_biotrees'
             and e['reference_set_size'] <= 24]
    small = sorted(small, key=lambda e: e['reference_set_size'])[:limit_b]
    for entry in small:
        case = common.load_case(entry)
        obs = case['observed']
        for k in (1, 2):
            jobs.append(dict(meta=dict(arm='B_exhaustive', name=entry['name'],
                                       family=entry['family'], k=k,
                                       n_edges=len(case['edges']),
                                       n_retained=len(obs)),
                             spec=dict(edges=case['edges'], obs=obs, k=k,
                                       eps=None, mode='exhaustive')))
            jobs.append(dict(meta=dict(arm='C_count', name=entry['name'],
                                       family=entry['family'], k=k,
                                       n_edges=len(case['edges']),
                                       n_retained=len(obs)),
                             spec=dict(edges=case['edges'], obs=obs, k=k,
                                       eps=F(0), mode='count')))
    return jobs


def main():
    assert PROTOCOL.exists(), 'register the protocol first'
    tag = sys.argv[1]
    limit_a = None if sys.argv[2] == 'all' else int(sys.argv[2])
    limit_b = int(sys.argv[3]) if len(sys.argv) > 3 else 12
    stream = OUT / ('stage301_validate_runs_%s.jsonl' % tag)
    jobs = build_jobs(limit_a, limit_b)
    print('jobs', len(jobs), flush=True)
    results = common.run_pool(jobs, worker, 1800, ADDRESS, WORKERS, str(stream))
    print('done', len(results), flush=True)


if __name__ == '__main__':
    main()
