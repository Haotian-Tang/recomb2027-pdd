#!/usr/bin/env python3
"""Stage292 panel ladder: is the adversarial erasure threshold flat in panel size?

Planted (NOT published) mixtures on the frozen Stage287 reference panels.
Planting rule, fully deterministic given (panel, size, K, seed):
  * the panel Newick is parsed by the frozen Stage286 parser, the root
    artefact is suppressed exactly as Stage290 does, and every degree-two
    vertex is suppressed so the tree is reduced;
  * for a ladder size L < |panel|, L leaves are drawn uniformly without
    replacement with random.Random(10**6 + 1000*size_rank + seed) and the
    induced subtree (frozen `observed_backbone`) is taken and renumbered;
  * K edges are drawn uniformly without replacement from ALL edges with
    random.Random(seed); weights are w_i / sum(w) with w_i uniform on 1..9;
    positions are x_i = a_i / 12 with a_i uniform on 1..11;
  * the pendant is y = min_r (sum_i p_i dist(r, q_i)) / 10, so every
    observation is exact rational and the instance is exactly explainable by
    construction with the planted support.
These are planted instances used to vary panel size at fixed generative
process; they are not published mixtures and carry no biological claim.
"""
import json
import multiprocessing as mp
import os
import random
import sys
import time
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
sys.dont_write_bytecode = True

from stage26_missing_leaf import observed_backbone                      # noqa: E402
from stage290_erasure_core import build_struct                          # noqa: E402
from stage290_survey import panel_edges                                 # noqa: E402
from stage292_all_minimisers import all_minimum_supports, VertexDistances  # noqa: E402
from stage292_harness import encode, isolated_run, sha                  # noqa: E402
from stage292_orders import order_for                                   # noqa: E402

WORKERS = int(os.environ.get('STAGE292_WORKERS', '26'))
OUT = ROOT / 'results/stage292_panels.jsonl'
MEMORY = 10 << 30
MIN_REMAINING = 3
REGIMES = (('adversarial_nearest', 0), ('random', 11))
BIG_PANEL = 'data/stage287_instances/krepp_shared/results__expt-emp__emp-pp-WoLv2__placement_tree_extended.nwk'
WOL2 = 'data/stage287_instances/krepp_shared/misc_data__backbone_blen_conserved-WoLv2.nwk'
LADDER = (64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 47151)
WOL_LADDER = (100, 1000, 15953)
KS = (3, 5)
SEEDS = (101, 202, 303)


def suppress_degree_two(edges):
    adj = {}
    for u, v, L in edges:
        L = F(L)
        adj.setdefault(u, {})[v] = L
        adj.setdefault(v, {})[u] = L
    for v in list(adj):
        if v in adj and len(adj[v]) == 2:
            (a, La), (b, Lb) = list(adj[v].items())
            del adj[a][v], adj[b][v], adj[v]
            adj[a][b] = adj[b][a] = La + Lb
    out, seen = [], set()
    for u in sorted(adj):
        for v, L in adj[u].items():
            key = (min(u, v), max(u, v))
            if key not in seen:
                seen.add(key)
                out.append((u, v, L))
    return out


def renumber(edges):
    adj = {}
    for u, v, _ in edges:
        adj.setdefault(u, []).append(v)
        adj.setdefault(v, []).append(u)
    leaves = sorted(v for v in adj if len(adj[v]) == 1)
    rest = sorted(v for v in adj if len(adj[v]) != 1)
    index = {v: i for i, v in enumerate(leaves)}
    index.update({v: i + len(leaves) for i, v in enumerate(rest)})
    return ([(index[u], index[v], F(L)) for u, v, L in edges],
            list(range(len(leaves))))


def build_panel(path, size, size_rank, seed):
    edges, leaves, _ = panel_edges(ROOT / path)
    edges = suppress_degree_two(edges)
    edges, leaves = renumber(edges)
    if size is not None and size < len(leaves):
        rng = random.Random(10 ** 6 + 1000 * size_rank + seed)
        keep = sorted(rng.sample(leaves, size))
        B = observed_backbone([(u, v, F(L)) for u, v, L in edges], keep)
        edges, leaves = renumber(B['edges'])
    return edges, leaves


def plant(edges, leaves, K, seed):
    rng = random.Random(seed)
    dists = VertexDistances(edges)
    chosen = sorted(rng.sample(range(len(edges)), K))
    raw = [rng.randint(1, 9) for _ in chosen]
    total = sum(raw)
    sources = [{'edge': i, 'p': F(w, total), 'x': F(rng.randint(1, 11), 12)}
               for i, w in zip(chosen, raw)]
    base = {r: F(0) for r in leaves}
    for s in sources:
        u, v, L = edges[s['edge']]
        x, p = s['x'], s['p']
        cu, cv = dists.column(u), dists.column(v)
        for r in leaves:
            base[r] += p * min(cu[r] + x * L, cv[r] + (1 - x) * L)
    y = min(base.values()) / 10
    assert y > 0
    return sources, y, {r: y + base[r] for r in leaves}, dists


def find_threshold(ok, emax, probes):
    """Stage290 convention: exponential bracketing, then bisection, both
    bracket endpoints explicitly evaluated."""
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
    t0 = time.perf_counter()
    edges, leaves = build_panel(spec['panel'], spec['size'], spec['size_rank'],
                                spec['seed'])
    sources, y, observed, dists = plant(edges, leaves, spec['K'], spec['seed'])
    support = sorted(s['edge'] for s in sources)
    K = spec['K']
    screen = 48 if len(leaves) > 2000 else 0
    build_seconds = time.perf_counter() - t0
    base = all_minimum_supports(edges, observed, cap=K, screen_rows=screen,
                                distances=dists)
    out = dict(kind='panel', panel=spec['panel'], size=spec['size'],
               K=K, seed=spec['seed'], panel_leaves=len(leaves),
               n_edges=len(edges), true_support=support,
               pendant=y, build_seconds=build_seconds,
               baseline_status=base['status'], baseline_minimum=base['minimum'],
               baseline_n_supports=base['n_supports'],
               baseline_unique_true=bool(base['status'] == 'complete' and
                                         base['minimum'] == K and
                                         base['n_supports'] == 1 and
                                         list(base['supports'][0]) == support),
               baseline_seconds=None, regimes=[])
    st = build_struct(edges, support, leaves)
    emax = len(leaves) - MIN_REMAINING
    for regime, seed in REGIMES:
        order = order_for(edges, st, regime, seed, sources=sources, dists=dists)
        probes = {}

        def ok(e):
            if e in probes:
                return probes[e][0]
            removed = set(order[:e])
            sub = {r: observed[r] for r in leaves if r not in removed}
            t1 = time.perf_counter()
            res = all_minimum_supports(edges, sub, cap=K, stop_after=2,
                                       screen_rows=screen, distances=dists)
            good = bool(res['status'] in ('complete',) and res['minimum'] == K and
                        res['n_supports'] == 1 and
                        list(res['supports'][0]) == support)
            probes[e] = (good, res['minimum'], res['n_supports'], res['status'],
                         time.perf_counter() - t1)
            return good

        t2 = time.perf_counter()
        t_unknown, status = find_threshold(ok, emax, probes)
        out['regimes'].append(dict(
            regime=regime, seed=seed, t_unknown=t_unknown, emax=emax,
            threshold_status=status, n_probes=len(probes),
            fraction=(None if t_unknown is None else t_unknown / len(leaves)),
            seconds=time.perf_counter() - t2,
            probe_table=sorted([dict(e=e, ok=v[0], minimum=v[1], n_supports=v[2],
                                     status=v[3], seconds=v[4])
                                for e, v in probes.items()], key=lambda z: z['e'])))
    return [out]


def main():
    jobs = []
    for rank, size in enumerate(LADDER):
        for K in KS:
            for seed in SEEDS:
                jobs.append(dict(panel=BIG_PANEL, size=size, size_rank=rank,
                                 K=K, seed=seed))
    for rank, size in enumerate(WOL_LADDER):
        for K in KS:
            for seed in SEEDS:
                jobs.append(dict(panel=WOL2, size=size, size_rank=100 + rank,
                                 K=K, seed=seed))
    jobs.sort(key=lambda j: j['size'])
    print('%d panel jobs' % len(jobs), flush=True)
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
                print('%d/%d %s' % (done, len(jobs), res['status']), flush=True)
    print('written', OUT, flush=True)


def _run(spec):
    wall = 600.0 if spec['size'] <= 2048 else 7200.0
    res = isolated_run(job, spec, wall, MEMORY)
    res['spec'] = spec
    return res


if __name__ == '__main__':
    main()

