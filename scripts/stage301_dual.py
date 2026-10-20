#!/usr/bin/env python3
"""Stage 301 step 6: a certified lower bound on Theta_k for EVERY k.

Relax the support entirely.  For any mass m >= 0 on the vertex set with
sum m = 1 and any y >= 0 the predictions are y + sum_v m_v d(v, r), and every
size-k configuration is such a mass (endpoint masses of its k edges).  So

    Theta_inf := min_{m, y} max_r |pred_r - d_r|   satisfies  Theta_inf <= Theta_k

for every k.  Theta_inf is one LP with |V| + 2 columns and 2|R| rows.  Its LP
DUAL gives a certificate that needs no optimality: for ANY rational vector s on
the observed leaves with

        sum_r |s_r| <= 1      and      sum_r s_r <= 0,

weak duality gives the exact rational bound

    L(s) = sum_r d_r s_r  -  max_{v in V} sum_r d(v, r) s_r   <=  Theta_inf.

(The two-leaf case s = (+1/2, -1/2) reduces to (d_r - d_s - d(r,s))/2.)

A good s is found in FLOATING POINT with HiGHS on a subset of leaves; the bound
it yields is then recomputed in EXACT rational arithmetic, streaming one frozen
BFS per leaf, so the float solver only proposes and never certifies.
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
WORKERS = 20
ADDRESS = 24 * 1024 ** 3
DEN = 10 ** 9


def repair(svals):
    """Rational s with sum|s| <= 1 and sum s <= 0, from any float vector."""
    s = [F(round(x * DEN), DEN) for x in svals]
    pos = sum(x for x in s if x > 0)
    neg = -sum(x for x in s if x < 0)
    if pos == 0 and neg == 0:
        return s
    if pos > neg:
        scale = neg / pos
        s = [x * scale if x > 0 else x for x in s]
        pos = neg
    total = pos + neg
    if total > 0:
        s = [x / total for x in s]
    return s


def certify(edges, leaves, dvals, Q, s):
    """Exact L(s); one frozen BFS per leaf with s_r != 0, streamed."""
    acc = {}
    live = [(q, sv) for q, sv in zip(Q, s) if sv != 0]
    for q, sv in live:
        dist = common.bfs_from(edges, [q])[q]
        if not acc:
            acc = {v: sv * d for v, d in dist.items()}
        else:
            for v, d in dist.items():
                acc[v] += sv * d
    if not acc:
        return F(0), 0
    obj = sum(dvals[leaves.index(q)] * sv for q, sv in live)
    return obj - max(acc.values()), len(live)


def propose(edges, obs, Q, seed_vertices=None):
    """Float HiGHS dual on the Q-restricted problem; mass free over all V."""
    import numpy as np
    from scipy.optimize import linprog
    verts = sorted({x for u, v, _l in edges for x in (u, v)})
    idx = {v: i for i, v in enumerate(verts)}
    nV = len(verts)
    M = np.zeros((len(Q), nV))
    adj = {}
    for u, v, l in edges:
        adj.setdefault(u, {})[v] = float(l)
        adj.setdefault(v, {})[u] = float(l)
    for j, q in enumerate(Q):
        d = s301._float_dist(adj, q)
        for v, x in d.items():
            M[j, idx[v]] = x
    dq = np.array([float(obs[q]) for q in Q])
    n = nV + 2                       # m, y, tau
    A = np.zeros((2 * len(Q), n))
    b = np.zeros(2 * len(Q))
    A[:len(Q), :nV] = -M
    A[:len(Q), nV] = -1.0
    A[:len(Q), nV + 1] = -1.0
    b[:len(Q)] = -dq
    A[len(Q):, :nV] = M
    A[len(Q):, nV] = 1.0
    A[len(Q):, nV + 1] = -1.0
    b[len(Q):] = dq
    cost = np.zeros(n)
    cost[nV + 1] = 1.0
    Aeq = np.zeros((1, n))
    Aeq[0, :nV] = 1.0
    res = linprog(cost, A_ub=A, b_ub=b, A_eq=Aeq, b_eq=[1.0],
                  bounds=[(0, None)] * n, method='highs')
    if not res.success:
        return None, dict(status=res.message[:120])
    marg = res.ineqlin.marginals          # <= 0
    u = -marg[:len(Q)]
    w = -marg[len(Q):]
    return (u - w), dict(status='optimal', float_value=float(res.fun))


def pick_leaves(obs, nq, seed):
    import random
    leaves = sorted(obs)
    if len(leaves) <= nq:
        return leaves
    order = sorted(leaves, key=lambda r: obs[r])
    picked = set(order[:nq // 4] + order[-nq // 4:])
    rng = random.Random(seed)
    while len(picked) < nq:
        picked.add(rng.choice(leaves))
    return sorted(picked)


def worker(spec):
    t0 = time.perf_counter()
    case = json.loads(Path(spec['path']).read_text())
    edges = [(u, v, F(l)) for u, v, l in case['edges']]
    obs = {int(r): F(d) for r, d in case['observed'].items()}
    leaves = sorted(obs)
    dvals = [obs[r] for r in leaves]
    best, best_meta = F(0), None
    for nq in spec['nq_ladder']:
        Q = pick_leaves(obs, nq, spec['seed'])
        sv, meta = propose(edges, obs, Q)
        if sv is None:
            continue
        s = repair(list(sv))
        t1 = time.perf_counter()
        bound, nlive = certify(edges, leaves, dvals, Q, s)
        if bound > best:
            best, best_meta = bound, dict(nq=len(Q), n_live=nlive,
                                          float_value=meta.get('float_value'),
                                          certify_seconds=time.perf_counter() - t1)
    return dict(status='complete', bound=best, meta=best_meta,
                n_edges=len(edges), n_leaves=len(leaves),
                seconds=time.perf_counter() - t0)


def main():
    assert (OUT / 'stage301_protocol.json').exists()
    nq = [int(x) for x in (sys.argv[1] if len(sys.argv) > 1 else '40,120').split(',')]
    proto = json.loads((ROOT / 'results/stage301/stage301_protocol.json').read_text())
    jobs = []
    for row in proto['td_cases']:
        jobs.append(dict(meta=dict(arm='td', name=Path(row['path']).stem),
                         spec=dict(path=str(ROOT / row['path']), nq_ladder=nq,
                                   seed=301)))
    sweeps = [json.loads(l) for l in
              (ROOT / 'results/stage298/stage298_sweep_runs.jsonl').read_text().splitlines()]
    known = {}
    for r in sweeps:
        if r['result'].get('theta') is not None:
            known.setdefault(r['name'], {})[r['k']] = r['result']['theta']
    index = common.bank_index()
    manifest = json.loads((ROOT / 'results/stage298/stage298_placement_manifest.json').read_text())
    pl = {i['name']: i['path'] for i in manifest['cases']}
    n_added = 0
    for name, th in sorted(known.items()):
        if n_added >= int(sys.argv[2] if len(sys.argv) > 2 else 120):
            break
        entry = next((e for e in index if e['name'] == name), None)
        path = (ROOT / entry['path']) if entry else (ROOT / pl[name] if name in pl else None)
        if path is None or not Path(path).exists():
            continue
        n_added += 1
        jobs.append(dict(meta=dict(arm='krepp', name=name, known_thetas=th),
                         spec=dict(path=str(path), nq_ladder=nq, seed=301)))
    print('jobs', len(jobs), flush=True)
    common.run_pool(jobs, worker, 3600, ADDRESS, WORKERS,
                    str(OUT / 'stage301_dual_runs.jsonl'))
    print('done', flush=True)


if __name__ == '__main__':
    main()
