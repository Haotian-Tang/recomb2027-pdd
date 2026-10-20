#!/usr/bin/env python3
"""Stage 301 step 6b: Theta_inf by exact column generation, certified.

Theta_inf = min over ALL mass distributions on the vertex set (sum = 1, m >= 0)
and y >= 0 of max_r |pred_r - d_r|.  Every size-k configuration is such a mass,
so Theta_inf <= Theta_k for every k, and with stage 298's Theta_1 this brackets
Theta_2, Theta_3, ... from both sides.

Column generation.  The restricted master keeps all 2|R| tolerance rows and a
small vertex set W; HiGHS proposes in floating point.  Pricing is exact and
costs O(|V|): for a signed leaf vector s the potential

    Phi_s(v) = sum_r s_r d(v, r)      obeys   Phi_s(c) = Phi_s(v) + l (S - 2 A_c)

down every edge (v, c) of length l, with S = sum_r s_r and A_c the s-mass of
c's subtree.  Two tree passes give Phi_s at every vertex.

Certification.  For ANY rational s with sum|s_r| <= 1 and sum s_r <= 0, weak
duality gives the exact rational bound

    L(s) = sum_r d_r s_r - max_v Phi_s(v)   <=   Theta_inf   <=   Theta_k.

The float solver only ever proposes s; L(s) is recomputed in Fractions.
"""
import json
import sys
import time
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import stage288_common as common            # noqa: E402
import stage301_common as s301              # noqa: E402

OUT = s301.OUT
WORKERS = 12
ADDRESS = 28 * 1024 ** 3
DEN = 10 ** 12


def rooted(edges):
    adj = {}
    for u, v, l in edges:
        adj.setdefault(u, {})[v] = l
        adj.setdefault(v, {})[u] = l
    root = min(adj)
    parent, order = {root: None}, [root]
    for v in order:
        for u in adj[v]:
            if u not in parent:
                parent[u] = v
                order.append(u)
    return adj, parent, order, root


def potential(adj, parent, order, root, s, zero):
    """Phi_s(v) = sum_r s_r d(v, r) at every vertex, two passes, any numeric type."""
    A = {v: s.get(v, zero) for v in adj}
    for v in reversed(order[1:]):
        A[parent[v]] += A[v]
    S = A[root]
    B = {v: zero for v in adj}
    for v in reversed(order[1:]):
        B[parent[v]] += B[v] + A[v] * adj[v][parent[v]]
    phi = {root: B[root]}
    for v in order[1:]:
        p = parent[v]
        phi[v] = phi[p] + adj[v][p] * (S - 2 * A[v])
    return phi


def repair(sv, leaves):
    s = {r: F(round(x * DEN), DEN) for r, x in zip(leaves, sv) if round(x * DEN) != 0}
    pos = sum(x for x in s.values() if x > 0)
    neg = -sum(x for x in s.values() if x < 0)
    if pos == 0 and neg == 0:
        return {}
    if pos > neg:
        scale = neg / pos
        s = {r: (x * scale if x > 0 else x) for r, x in s.items()}
        pos = neg
    total = pos + neg
    if total > 0:
        s = {r: x / total for r, x in s.items()}
    return s


def solve(edges, obs, rounds, add_per_round, seed, time_budget):
    import numpy as np
    from scipy.optimize import linprog
    leaves = sorted(obs)
    dq = np.array([float(obs[r]) for r in leaves])
    adjf = {}
    for u, v, l in edges:
        adjf.setdefault(u, {})[v] = float(l)
        adjf.setdefault(v, {})[u] = float(l)
    adjF, parent, order, root = rooted([(u, v, F(l)) for u, v, l in edges])
    verts = sorted(adjf)
    # seed columns: the two ends of every edge attaining a small single-edge fit
    import random
    rng = random.Random(seed)
    W = sorted(set(rng.sample(verts, min(len(verts), 24))))
    cols = {}
    t0 = time.perf_counter()
    best, best_meta = F(0), None
    history = []
    for it in range(rounds):
        for w in W:
            if w not in cols:
                d = s301._float_dist(adjf, w)
                cols[w] = np.array([d[r] for r in leaves])
        M = np.stack([cols[w] for w in W], axis=1)      # |R| x |W|
        nW = len(W)
        n = nW + 2
        A = np.zeros((2 * len(leaves), n))
        A[:len(leaves), :nW] = -M
        A[:len(leaves), nW] = -1.0
        A[:len(leaves), nW + 1] = -1.0
        A[len(leaves):, :nW] = M
        A[len(leaves):, nW] = 1.0
        A[len(leaves):, nW + 1] = -1.0
        b = np.concatenate([-dq, dq])
        cost = np.zeros(n)
        cost[nW + 1] = 1.0
        Aeq = np.zeros((1, n))
        Aeq[0, :nW] = 1.0
        res = linprog(cost, A_ub=A, b_ub=b, A_eq=Aeq, b_eq=[1.0],
                      bounds=[(0, None)] * n, method='highs')
        if not res.success:
            break
        marg = res.ineqlin.marginals
        sv = (-marg[:len(leaves)]) - (-marg[len(leaves):])
        s = repair(list(sv), leaves)
        if s:
            phi = potential(adjF, parent, order, root, s, F(0))
            L = sum(F(obs[r]) * x for r, x in s.items()) - max(phi.values())
            if L > best:
                best, best_meta = L, dict(iteration=it, n_cols=nW,
                                          float_master=float(res.fun),
                                          n_live=len(s))
        history.append(dict(iteration=it, n_cols=nW, master=float(res.fun),
                            bound=str(best)))
        # price every vertex exactly in O(|V|): add the largest Phi_s
        if not s:
            break
        phif = potential(adjf, parent, order, root,
                         {r: float(x) for r, x in s.items()}, 0.0)
        cand = sorted((v for v in verts if v not in set(W)),
                      key=lambda v: -phif[v])[:add_per_round]
        if not cand:
            break
        # keep the master small: basic columns plus the new candidates
        keep = {w for w, mv in zip(W, res.x[:nW]) if mv > 1e-12}
        W = sorted(keep | set(cand) | set(W[:24]))
        if time.perf_counter() - t0 > time_budget:
            break
    return best, best_meta, history, time.perf_counter() - t0


def worker(spec):
    case = json.loads(Path(spec['path']).read_text())
    edges = [(u, v, F(l)) for u, v, l in case['edges']]
    obs = {int(r): F(d) for r, d in case['observed'].items()}
    L, meta, hist, el = solve(edges, obs, spec['rounds'], spec['add'],
                              spec['seed'], spec['budget'])
    return dict(status='complete', bound=L, meta=meta, history=hist[-6:],
                n_edges=len(edges), n_leaves=len(obs), seconds=el)


def main():
    assert (OUT / 'stage301_protocol.json').exists()
    rounds = int(sys.argv[1]) if len(sys.argv) > 1 else 14
    add = int(sys.argv[2]) if len(sys.argv) > 2 else 12
    budget = float(sys.argv[3]) if len(sys.argv) > 3 else 1200.0
    nkrepp = int(sys.argv[4]) if len(sys.argv) > 4 else 120
    proto = json.loads((OUT / 'stage301_protocol.json').read_text())
    jobs = [dict(meta=dict(arm='td', name=Path(r['path']).stem),
                 spec=dict(path=str(ROOT / r['path']), rounds=rounds, add=add,
                           seed=301, budget=budget))
            for r in proto['td_cases']]
    sweeps = [json.loads(l) for l in
              (ROOT / 'results/stage298/stage298_sweep_runs.jsonl').read_text().splitlines()]
    known = {}
    for r in sweeps:
        if r['result'].get('theta') is not None:
            known.setdefault(r['name'], {})[r['k']] = r['result']['theta']
    index = common.bank_index()
    manifest = json.loads((ROOT / 'results/stage298/stage298_placement_manifest.json').read_text())
    pl = {i['name']: i['path'] for i in manifest['cases']}
    n = 0
    for name, th in sorted(known.items()):
        if n >= nkrepp:
            break
        entry = next((e for e in index if e['name'] == name), None)
        path = (ROOT / entry['path']) if entry else (ROOT / pl[name] if name in pl else None)
        if path is None or not Path(path).exists():
            continue
        n += 1
        jobs.append(dict(meta=dict(arm='krepp', name=name, known_thetas=th),
                         spec=dict(path=str(path), rounds=rounds, add=add,
                                   seed=301, budget=budget)))
    print('jobs', len(jobs), flush=True)
    common.run_pool(jobs, worker, budget + 900, ADDRESS, WORKERS,
                    str(OUT / 'stage301_thetainf_runs.jsonl'))
    print('done', flush=True)


if __name__ == '__main__':
    main()
