#!/usr/bin/env python3
"""Stage292 extended bridge: the unknown-support erasure threshold at K >= 5.

Stage290 could only measure t_unknown at K = 2 (|R| <= 130) and K = 3
(|R| <= 45) because it enumerated all minimum-size supports by brute force at
cost C(|E|, K).  This script repeats exactly the same measurement with the
Stage292 structural all-minimisers routine, over the WHOLE exact bank:
K in {2, 3, 5, 7, 10, 12} at every reference-set size the bank provides.

Definitions and conventions are Stage290's, unchanged:
  t_known    largest e with rank(leaf-report matrix after erasing the first e
             references of the order) still s + 1, found by an exhaustive scan
  t_unknown  largest e for which the set of ALL minimum-size supports is
             exactly the single true support, found by exponential bracketing
             then bisection with both bracket endpoints evaluated
  regimes    random (seeds 11, 23), whole_clade (seeds 11, 23),
             adversarial_nearest (deterministic)
  MIN_REMAINING = 3 references are never erased.
A probe whose status is not decisive aborts that regime with
threshold_status='unknown_probe'; it is never silently read as a failure.
"""
import json
import multiprocessing as mp
import os
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
sys.dont_write_bytecode = True

from stage290_erasure_core import (load_bank_case, support_of, build_struct,  # noqa: E402
                                   erasure_analysis, distances_mod, rank_mod,
                                   PRIME)
from stage292_all_minimisers import all_minimum_supports, VertexDistances  # noqa: E402
from stage292_harness import encode, isolated_run                        # noqa: E402
from stage292_orders import order_for                                    # noqa: E402

WORKERS = int(os.environ.get('STAGE292_WORKERS', '28'))
OUT = ROOT / (os.environ.get('STAGE292_OUT') or
              ('results/stage292_bridge%s.jsonl' %
               (('_K' + '_'.join(sys.argv[1:])) if len(sys.argv) > 1 else '')))
SKIP = os.environ.get('STAGE292_SKIP', '')
MEMORY = 10 << 30
MIN_REMAINING = 3
REGIMES = (('random', 11), ('random', 23), ('whole_clade', 11),
           ('whole_clade', 23), ('adversarial_nearest', 0))
PROBE_SECONDS = {2: 120, 3: 120, 5: 300, 7: 600, 10: 240, 12: 240}
STRIDE = int(os.environ.get('STAGE292_STRIDE', '1'))


class ProbeUnknown(Exception):
    pass


class RankScanner(object):
    """Stage290's rank-after-prefix-erasure scanner (same construction)."""

    def __init__(self, st, p=PRIME):
        self.p = p
        Plist, Rlist = sorted(st['P']), sorted(st['R'])
        dmod, index = distances_mod(st, Plist, p)
        nR, nP = len(Rlist), len(Plist)
        M = np.empty((nR + 1, nP + 1), dtype=np.int64)
        cols = np.array([index[r] for r in Rlist], dtype=np.int64)
        for j in range(nP):
            M[:nR, j] = dmod[j, cols]
        M[:nR, nP] = 1
        M[nR, :nP] = 1
        M[nR, nP] = 0
        self.M, self.nR = M, nR
        self.row_of = {r: i for i, r in enumerate(Rlist)}

    def rank(self, removed):
        drop = {self.row_of[r] for r in removed}
        keep = [i for i in range(self.nR) if i not in drop] + [self.nR]
        return rank_mod(self.M[keep], self.p)


def find_threshold(ok, emax, probes):
    if emax < 1:
        return None, 'degenerate'
    if not ok(1):
        return 0, 'found'
    lo, hi = 1, 2
    while hi <= emax and ok(hi):
        lo, hi = hi, hi * 2
    if hi > emax:
        if ok(emax):
            return None, 'censored'
        hi = emax
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if ok(mid):
            lo = mid
        else:
            hi = mid
    return lo, 'found'


def job(spec):
    case = load_bank_case(ROOT / spec['path'])
    support = support_of(case)
    K = len(support)
    observed = case['observed']
    edges = case['edges']
    leaves = sorted(observed)
    st = build_struct(edges, support, leaves)
    dists = VertexDistances(edges)
    seconds = PROBE_SECONDS.get(K, 900)
    screen = 48 if len(leaves) > 2000 else 0
    known = erasure_analysis(st, want_witness=False)
    scanner = RankScanner(st)
    rank_full = scanner.rank([])
    t0 = time.perf_counter()
    base = all_minimum_supports(edges, observed, cap=K, seconds=seconds,
                                screen_rows=screen, distances=dists)
    out = dict(kind='bridge', name=spec['name'], family=spec['family'],
               path=spec['path'], sha256=spec['sha256'], K=K,
               n_edges=len(edges), R_size=len(observed),
               original_leaves=case['original_leaves'], true_support=support,
               support_is_matching=st['matching'],
               n_pendant_support=len(st['pendant_support']),
               cogirth=known['delta'], cogirth_tolerance=known['tolerance'],
               s=known['s'], rank_full=rank_full,
               rank_identity=bool(rank_full == known['s'] + 1),
               baseline_status=base['status'], baseline_minimum=base['minimum'],
               baseline_n_supports=base['n_supports'],
               baseline_unique_true=bool(base['status'] == 'complete' and
                                         base['minimum'] == K and
                                         base['n_supports'] == 1 and
                                         list(base['supports'][0]) == support),
               baseline_seconds=time.perf_counter() - t0,
               frozen_minimum=base.get('frozen_minimum'), regimes=[])
    emax = len(observed) - MIN_REMAINING
    for regime, seed in REGIMES:
        order = order_for(edges, st, regime, seed,
                          sources=case['generating_sources'], dists=dists)
        t_known = None
        for e in range(1, emax + 1):
            if scanner.rank(order[:e]) != known['s'] + 1:
                t_known = e - 1
                break
        probes = {}

        def ok(e):
            if e in probes:
                if probes[e][0] is None:
                    raise ProbeUnknown(probes[e][3])
                return probes[e][0]
            removed = set(order[:e])
            sub = {r: observed[r] for r in leaves if r not in removed}
            t1 = time.perf_counter()
            res = all_minimum_supports(edges, sub, cap=K, seconds=seconds,
                                       stop_after=2, screen_rows=screen,
                                       distances=dists)
            if res['status'] == 'complete':
                good = bool(res['minimum'] == K and res['n_supports'] == 1 and
                            list(res['supports'][0]) == support)
            elif res['status'] == 'early_exit':
                good = False
            else:
                good = None
            probes[e] = (good, res['minimum'], res['n_supports'], res['status'],
                         time.perf_counter() - t1)
            if good is None:
                raise ProbeUnknown(res['status'])
            return good

        t2 = time.perf_counter()
        try:
            t_unknown, status = find_threshold(ok, emax, probes)
        except ProbeUnknown as exc:
            t_unknown, status = None, 'unknown_probe:' + str(exc)
        out['regimes'].append(dict(
            regime=regime, seed=seed, t_known_sequence=t_known,
            t_unknown_sequence=t_unknown, emax=emax, threshold_status=status,
            n_probes=len(probes), seconds=time.perf_counter() - t2,
            fraction=(None if t_unknown is None else t_unknown / len(observed)),
            probe_table=sorted([dict(e=e, ok=v[0], minimum=v[1], n_supports=v[2],
                                     status=v[3], seconds=v[4])
                                for e, v in probes.items()],
                               key=lambda z: z['e'])))
    return [out]


def _run(spec):
    res = isolated_run(job, spec, spec['wall'], MEMORY)
    res['spec'] = spec
    return res


def main():
    index = json.loads((ROOT / 'data/stage287_instances/bank/index.json').read_text())['cases']
    want = {int(a) for a in sys.argv[1:]} or None
    jobs = []
    for x in index:
        if not x['family'].startswith('exact'):
            continue
        K = x['K_generating']
        if want is not None and K not in want:
            continue
        wall = 1800.0 if K <= 5 else (5400.0 if K <= 7 else 1200.0)
        jobs.append(dict(wall=wall, **x))
    jobs.sort(key=lambda j: (j['K_generating'], j['reference_set_size']))
    if STRIDE > 1:
        jobs = jobs[::STRIDE]
    if SKIP:
        done = set()
        for part in SKIP.split(','):
            p = ROOT / part
            if p.exists():
                for line in p.open():
                    rec = json.loads(line)
                    if rec.get('kind') == 'bridge':
                        done.add(rec['name'])
        jobs = [j for j in jobs if j['name'] not in done]
    print('%d bridge jobs (stride %d)' % (len(jobs), STRIDE), flush=True)
    with OUT.open('x') as fh:
        with mp.Pool(WORKERS, maxtasksperchild=1) as pool:
            done = 0
            for res in pool.imap_unordered(_run, jobs, chunksize=1):
                if res['status'] == 'ok':
                    for row in res['value']:
                        fh.write(json.dumps(encode(row)) + '\n')
                else:
                    fh.write(json.dumps(encode(res)) + '\n')
                fh.flush()
                done += 1
                if done % 10 == 0 or done == len(jobs):
                    print('%d/%d' % (done, len(jobs)), flush=True)
    print('written', OUT, flush=True)


if __name__ == '__main__':
    main()

