#!/usr/bin/env python3
"""Stage 301 step 3: measured reach and measured pruning of the localisation.

Theta_k is computed by a LADDER: Theta_1 exhaustively over the |E| single-edge
supports, then Theta_k screened at eps = Theta_{k-1}, which is a valid upper
bound because Theta_k is non-increasing in k (stage 296).  Every surviving
support is decided by the frozen stage-296 LP, so the ladder is exact; the only
thing the screen changes is how many supports are visited.

Recorded per rung: candidates, C(|E|, k), the pruning ratio, seconds per
support, whether the screen was VACUOUS (candidates == C(|E|, k)), and the
exact Theta_k.  A rung is abandoned when the candidate count exceeds the
registered cap or the wall expires, and the abandoned rung IS the measured
ceiling for that case.
"""
import json
import sys
import time
from fractions import Fraction as F
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import stage288_common as common            # noqa: E402
import stage296_common as s296              # noqa: E402
import stage301_common as s301              # noqa: E402

OUT = s301.OUT
WORKERS = 40
ADDRESS = 10 * 1024 ** 3
CAP = 60000


def ladder(edges, obs, kmax, wall, cap=CAP):
    leaves = sorted(obs)
    dvals = tuple(obs[r] for r in leaves)
    Dv = s296.leaf_vectors(edges, leaves)
    view = s301.backbone_view(edges, obs)
    pos, neg = s301.critical_eps(view)
    rungs = []
    t_start = time.perf_counter()
    theta = None
    for k in range(1, kmax + 1):
        total = comb(len(edges), k)
        t0 = time.perf_counter()
        if k == 1:
            out = s296.theta_min_over_supports(Dv, edges, dvals, 1)
            n_cand, vacuous, visited = total, True, out['visited']
            screen = dict(vacuous=True, bootstrap=True)
        else:
            it, n_cand, screen = s301.candidate_plan(view, theta, k, len(edges),
                                                     hard_cap=cap)
            if it is None or (n_cand is not None and n_cand > cap):
                rungs.append(dict(k=k, status='ceiling_candidates', exhaustive=total,
                                  n_candidates=n_cand, screen=screen,
                                  seconds=time.perf_counter() - t0))
                break
            if n_cand == 0:
                rungs.append(dict(k=k, status='refuted_at_previous_theta',
                                  exhaustive=total, n_candidates=0, screen=screen,
                                  seconds=time.perf_counter() - t0))
                break
            best, best_S, best_z, hot, visited, pruned = theta, None, None, [], 0, 0
            for S in it:
                visited += 1
                r = s296.theta_support(s296.support_cols(Dv, edges, S), dvals, k,
                                       cutoff=best, seed=hot[:8])
                if r['status'] == 'pruned':
                    pruned += 1
                    continue
                if r['status'] != 'optimal':
                    return rungs, dict(status='lp_' + r['status'], k=k)
                hot = list(dict.fromkeys(r['active'] + hot))[:8]
                if best is None or r['value'] < best:
                    best, best_S, best_z = r['value'], S, r['z']
            out = dict(status='complete', theta=best, support=best_S, z=best_z,
                       visited=visited, pruned=pruned)
            vacuous = bool(screen.get('vacuous'))
        el = time.perf_counter() - t0
        theta = out['theta']
        row = dict(k=k, status='complete', theta=theta, support=out.get('support'),
                   exhaustive=total, n_candidates=n_cand, visited=out['visited'],
                   prune_ratio=(float(n_cand) / total if total else None),
                   vacuous=vacuous, screen=screen, seconds=el,
                   seconds_per_support=(el / max(1, out['visited'])))
        if out.get('support') is not None and out.get('z') is not None:
            dev = s296.independent_deviation(edges, out['support'], out['z'], k,
                                             obs, leaves)
            row['independent_deviation'] = dev
            row['witness_agrees'] = (dev == theta)
        rungs.append(row)
        if time.perf_counter() - t_start > wall:
            rungs.append(dict(k=k + 1, status='ceiling_wall'))
            break
    return rungs, dict(status='complete',
                       eps_vacuous_pos=(max(pos.values()) if pos else F(0)),
                       n_pos_total=len(pos), n_neg_total=len(neg))


def worker(spec):
    t0 = time.perf_counter()
    rungs, meta = ladder(spec['edges'], spec['obs'], spec['kmax'], spec['wall'])
    lengths = sorted(l for _u, _v, l in spec['edges'])
    return dict(status='complete', rungs=rungs, meta=meta,
                n_edges=len(spec['edges']), n_leaves=len(spec['obs']),
                median_edge_length=lengths[len(lengths) // 2],
                min_edge_length=lengths[0],
                seconds=time.perf_counter() - t0)


def main():
    assert (OUT / 'stage301_protocol.json').exists()
    per_cell = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    kmax = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    wall = float(sys.argv[3]) if len(sys.argv) > 3 else 1200.0
    global CAP
    if len(sys.argv) > 4:
        CAP = int(sys.argv[4])
    index = common.bank_index()
    exact = [e for e in index if e['family'] in ('exact_10k', 'exact_biotrees')]
    cells = {}
    for e in exact:
        n = e['reference_set_size']
        band = ('20-29' if n < 30 else '30-49' if n < 50 else '50-99' if n < 100
                else '100-199' if n < 200 else '200-361')
        cells.setdefault((e['K_generating'], band), []).append(e)
    jobs = []
    for (K, band), group in sorted(cells.items()):
        group = sorted(group, key=lambda e: e['name'])[:per_cell]
        for entry in group:
            case = common.load_case(entry)
            jobs.append(dict(meta=dict(arm='reach', name=entry['name'],
                                       family=entry['family'], K=K, band=band,
                                       n_edges=len(case['edges']),
                                       n_retained=len(case['observed'])),
                             spec=dict(edges=case['edges'], obs=case['observed'],
                                       kmax=min(kmax, K), wall=wall)))
    print('jobs', len(jobs), flush=True)
    common.run_pool(jobs, worker, wall + 600, ADDRESS, WORKERS,
                    str(OUT / 'stage301_reach_runs.jsonl'))
    print('done', flush=True)


if __name__ == '__main__':
    main()
