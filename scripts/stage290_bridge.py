#!/usr/bin/env python3
"""Stage290 bridge: from known-support erasure tolerance to unknown-support
identifiability.

With the support KNOWN the reference-erasure tolerance is cogirth - 1: the
number of reference reports that may be erased with the fixed-support fibre
dimension unchanged, worst case over every erasure set.  With the support
UNKNOWN - the manuscript's actual problem - the analogue is how many references
may be erased before the exact minimum-source answer stops being unique, that
is before either more than one support of minimum size explains the remaining
observations or the minimum size drops below the generating count.

Per case and removal regime this measures
  cogirth_tolerance   delta - 1, worst case over ALL erasure sets (Stage196 DP)
  t_known_sequence    largest e with rank(R minus first e) still s + 1, found by
                      an exhaustive scan over e
  t_unknown_sequence  largest e with exactly one minimum-size support, of the
                      generating size, found by exponential bracketing and
                      bisection with the bracket endpoints re-verified
and asks whether the known-support quantity predicts the unknown-support one.
Every reported support is confirmed by the exact rational Stage25/42 endpoint
polytope oracle.
"""
import json
import multiprocessing as mp
import os
import sys
import time
import traceback
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
sys.dont_write_bytecode = True

from stage43_partial_fpt import solve_partial                              # noqa: E402
from stage290_erasure_core import (load_bank_case, support_of, build_struct,   # noqa: E402
                                   erasure_analysis, min_supports, CaseTables,
                                   removal_order, encode, distances_mod,
                                   rank_mod, PRIME)

WORKERS = int(os.environ.get('STAGE290_WORKERS', '16'))
OUT = ROOT / 'results/stage290_bridge.jsonl'
REGIMES = (('random', 11), ('random', 23), ('whole_clade', 11),
           ('whole_clade', 23), ('adversarial_nearest', 0))
GROUND_TRUTH_EDGES = 90          # exhaustive exact enumeration below this size
MIN_REMAINING = 3                # never erase below this many references


class RankScanner(object):
    """Rank of the fixed-support leaf-report matrix after erasing a prefix."""

    def __init__(self, st, p=PRIME):
        self.p = p
        self.Plist = sorted(st['P'])
        self.Rlist = sorted(st['R'])
        dmod, index = distances_mod(st, self.Plist, p)
        nR, nP = len(self.Rlist), len(self.Plist)
        M = np.empty((nR + 1, nP + 1), dtype=np.int64)
        cols = np.array([index[r] for r in self.Rlist], dtype=np.int64)
        for j in range(nP):
            M[:nR, j] = dmod[j, cols]
        M[:nR, nP] = 1
        M[nR, :nP] = 1
        M[nR, nP] = 0
        self.M = M
        self.row_of = {r: i for i, r in enumerate(self.Rlist)}
        self.nR = nR

    def rank(self, removed):
        drop = {self.row_of[r] for r in removed}
        keep = [i for i in range(self.nR) if i not in drop] + [self.nR]
        return rank_mod(self.M[keep], self.p)


def unknown_ok(tables, observed, keys, K, true_support):
    sub = {r: observed[r] for r in keys}
    ms = min_supports(tables, sub, K, rng_seed=7)
    ok = (ms['status'] == 'complete' and ms['minimum'] == K and
          ms['n_supports'] == 1 and ms['supports'][0] == true_support)
    return bool(ok), ms


def find_threshold(tables, observed, order, K, true_support, emax, probes):
    """First e with a non-unique minimum-size answer; exponential + bisection.

    Bracketing assumes the failure set is upward closed in e along this
    sequence; both bracket endpoints are evaluated explicitly, so the reported
    threshold t always satisfies `ok at t` and `not ok at t+1`.
    """
    keys_all = sorted(observed)

    def ok(e):
        if e in probes:
            return probes[e][0]
        removed = set(order[:e])
        res = unknown_ok(tables, observed, [r for r in keys_all if r not in removed],
                         K, true_support)
        probes[e] = (res[0], res[1]['minimum'], res[1]['n_supports'], res[1]['status'])
        return res[0]

    if emax < 1:
        return None, 'degenerate'
    if not ok(1):
        return 0, 'found'
    lo, hi = 1, 2
    while hi <= emax and ok(hi):
        lo, hi = hi, hi * 2
    if hi > emax:
        if ok(emax):
            return None, 'censored'            # still unique at the largest e
        hi = emax
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if ok(mid):
            lo = mid
        else:
            hi = mid
    return lo, 'found'


def run_case(spec):
    t0 = time.perf_counter()
    try:
        case = load_bank_case(ROOT / spec['path'])
        support = support_of(case)
        K = len(support)
        observed = case['observed']
        st = build_struct(case['edges'], support, sorted(observed))
        tables = CaseTables(case['edges'])
        known = erasure_analysis(st, want_witness=True)
        scanner = RankScanner(st)
        rank_full = scanner.rank([])
        gt = len(case['edges']) <= GROUND_TRUTH_EDGES

        base = min_supports(tables, observed, K, rng_seed=7)
        base_ok = (base['status'] == 'complete' and base['minimum'] == K and
                   base['n_supports'] == 1 and base['supports'][0] == support)
        gt_agree = None
        if gt:
            exact = min_supports(tables, observed, K, rng_seed=7, exhaustive_exact=True)
            gt_agree = bool(exact['minimum'] == base['minimum'] and
                            sorted(map(tuple, exact['supports'] or [])) ==
                            sorted(map(tuple, base['supports'] or [])))
        frozen = solve_partial(case['edges'], observed, cap=K, seconds=None,
                               max_states=None, max_bases=None, max_faces=None)
        out = dict(kind='bridge', name=spec['name'], family=spec['family'],
                   path=spec['path'], sha256=spec['sha256'], K=K,
                   n_edges=len(case['edges']), R_size=len(observed),
                   original_leaves=case['original_leaves'], true_support=support,
                   support_is_matching=st['matching'],
                   n_pendant_support=len(st['pendant_support']),
                   cogirth=known['delta'], cogirth_tolerance=known['tolerance'],
                   s=known['s'], rank_full=rank_full,
                   rank_identity=bool(rank_full == known['s'] + 1),
                   fibre_dim=known['fibre_dim_endpoint'],
                   baseline_minimum=base['minimum'],
                   baseline_n_supports=base['n_supports'],
                   baseline_unique_true=bool(base_ok),
                   baseline_ground_truth=gt_agree,
                   prefilter_survivors=base['statistics']['prefilter_survivors'],
                   prefilter_candidates=base['statistics']['prefilter_candidates'],
                   prefilter_margin=base.get('prefilter_margin'),
                   frozen_solver_minimum=frozen.get('minimum_sources_within_cap'),
                   frozen_agrees=bool(frozen.get('minimum_sources_within_cap') ==
                                      base['minimum']),
                   regimes=[])
        emax = len(observed) - MIN_REMAINING
        for regime, seed in REGIMES:
            order = removal_order(case, st, regime, seed, tables=tables)
            # known-support side: exhaustive scan, rank is cheap
            t_known = None
            for e in range(1, emax + 1):
                if scanner.rank(order[:e]) != known['s'] + 1:
                    t_known = e - 1
                    break
            probes = {}
            t_unknown, tstatus = find_threshold(tables, observed, order, K,
                                                support, emax, probes)
            out['regimes'].append(dict(
                regime=regime, seed=seed, t_known_sequence=t_known,
                t_unknown_sequence=t_unknown, emax=emax,
                threshold_status=tstatus, probes=len(probes),
                probe_table=sorted([dict(e=e, ok=v[0], minimum=v[1],
                                         n_supports=v[2], status=v[3])
                                    for e, v in probes.items()],
                                   key=lambda z: z['e']),
                order_head=order[:20]))
        out['status'] = 'PASS'
    except Exception as exc:
        out = dict(kind='bridge', status='FAIL', spec=spec,
                   error_type=type(exc).__name__, error=str(exc),
                   traceback=traceback.format_exc())
    out['wall_seconds'] = time.perf_counter() - t0
    return out


def jobs():
    index = json.loads((ROOT / 'data/stage287_instances/bank/index.json').read_text())['cases']
    out = []
    for x in index:
        if x['family'] == 'exact_biotrees' and x['K_generating'] == 2 and x['reference_set_size'] <= 130:
            out.append(x)
        elif x['family'] == 'exact_10k' and x['K_generating'] == 3 and x['reference_set_size'] <= 45:
            out.append(x)
        elif x['family'] == 'exact_biotrees' and x['K_generating'] == 3 and x['reference_set_size'] <= 45:
            out.append(x)
    return sorted(out, key=lambda z: (z['K_generating'], z['reference_set_size']))


def main():
    allj = jobs()
    print('bridge jobs: %d' % len(allj), flush=True)
    with OUT.open('x') as fh:
        with mp.Pool(WORKERS, maxtasksperchild=4) as pool:
            done = 0
            for res in pool.imap_unordered(run_case, allj, chunksize=1):
                fh.write(json.dumps(encode(res)) + '\n')
                fh.flush()
                done += 1
                print('%d/%d %s %s %.0fs' % (done, len(allj), res.get('name'),
                                             res['status'], res['wall_seconds']),
                      flush=True)
    print('written', OUT)


if __name__ == '__main__':
    main()
