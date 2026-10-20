#!/usr/bin/env python3
"""Stage290 survey: reference-erasure tolerance and the fibre-dimension law
across the frozen Stage287 instance bank and the frozen large panels.

For every job the survey records
  s          minimum P-to-R vertex cut  (= rank of the leaf-report matrix - 1)
  delta      cogirth of the leaf-report matroid (Stage196 budgeted recurrence)
  tolerance  delta - 1, the number of references that may be erased with the
             fixed-support fibre dimension unchanged
and three independent checks: the modular rank identity rank = s + 1 under two
primes, the unreduced dynamic program, and exhaustive erasure of observed-leaf
subsets (brute force) wherever that is affordable.
"""
import json
import multiprocessing as mp
import os
import random
import sys
import time
import traceback
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
sys.dont_write_bytecode = True

from stage286_newick import parse_newick                                  # noqa: E402
from stage290_erasure_core import (load_bank_case, support_of, build_struct,  # noqa: E402
                                   erasure_analysis, erasure_analysis_full,
                                   brute_force_cogirth, fixed_support_rank,
                                   encode, PRIME, PRIME2, root_order)

WORKERS = int(os.environ.get('STAGE290_WORKERS', '16'))
OUT = ROOT / 'results/stage290_survey.jsonl'
BRUTE_LIMIT = 150000
FULL_DP_MAX_VERTICES = 3000
PANEL_K = (3, 5, 10, 12)
PANEL_REPLICATES = 5


def split_profile(st):
    """Leaf-split sizes of the support edges, and support spread."""
    adj, edges, R = st['adj'], st['edges'], st['R']
    r0 = min(adj)
    parent, order = root_order(adj, r0)
    cnt = dict.fromkeys(order, 0)
    for v in reversed(order):
        cnt[v] += int(v in R)
        if parent[v] is not None:
            cnt[parent[v]] += cnt[v]
    total = len(R)
    smalls = []
    for e in st['support']:
        u, v, _ = edges[e]
        child = v if parent.get(v) == u else u
        a = cnt[child]
        smalls.append(min(a, total - a))
    return dict(support_small_side=sorted(smalls), min_small_side=min(smalls),
                n_pendant_support=len(st['pendant_support']),
                any_pendant_support=bool(st['pendant_support']))


def analyse(st, want_full_dp, want_brute):
    out = erasure_analysis(st, want_witness=True)
    out.update(split_profile(st))
    rank1 = fixed_support_rank(st, PRIME)
    rank2 = fixed_support_rank(st, PRIME2)
    out.update(rank_p1=rank1, rank_p2=rank2,
               rank_identity=bool(rank1 == rank2 == out['s'] + 1))
    if want_full_dp:
        full = erasure_analysis_full(st)
        out.update(full_dp_s=full['s'], full_dp_delta=full['delta'],
                   full_dp_nodes=full['n_nodes'],
                   reduction_agrees=bool(full['s'] == out['s'] and
                                         full['delta'] == out['delta']))
    if want_brute:
        bf = brute_force_cogirth(st, limit=BRUTE_LIMIT)
        out.update(brute_status=bf['status'], brute_cogirth=bf.get('cogirth'),
                   brute_tested=bf['tested'], brute_full_rank=bf['full_rank'],
                   brute_agrees=(None if bf['status'] != 'complete'
                                 else bool(bf['cogirth'] == out['delta'])))
    return out


def panel_edges(path):
    tree, diag = parse_newick(Path(path).read_text())
    adj = {}
    for i in range(tree.n):
        p = tree.parent[i]
        if p is None:
            continue
        L = tree.elen[i]
        adj.setdefault(i, []).append((p, L))
        adj.setdefault(p, []).append((i, L))
    # suppress an unlabelled degree-1 or degree-2 root artefact
    root = tree.root
    if len(adj.get(root, [])) == 1 and not tree.label[root]:
        (w, L), = adj[root]
        adj[w] = [(x, y) for x, y in adj[w] if x != root]
        del adj[root]
    elif len(adj.get(root, [])) == 2 and not tree.label[root]:
        (a, La), (b, Lb) = adj[root]
        adj[a] = [(x, y) for x, y in adj[a] if x != root] + [(b, La + Lb)]
        adj[b] = [(x, y) for x, y in adj[b] if x != root] + [(a, La + Lb)]
        del adj[root]
    verts = sorted(adj)
    index = {v: i for i, v in enumerate(verts)}
    edges, seen = [], set()
    for v in verts:
        for w, L in adj[v]:
            key = (min(v, w), max(v, w))
            if key in seen:
                continue
            seen.add(key)
            assert L is not None and L > 0, 'panel must be positivity conformant'
            edges.append((index[v], index[w], L))
    leaves = [index[v] for v in verts if len(adj[v]) == 1]
    return edges, leaves, diag


def job_bank(spec):
    case = load_bank_case(ROOT / spec['path'])
    support = support_of(case)
    st = build_struct(case['edges'], support, sorted(case['observed']))
    nv = len(st['adj'])
    res = analyse(st, nv <= FULL_DP_MAX_VERTICES, True)
    res.update(kind='bank', family=spec['family'], name=spec['name'],
               path=spec['path'], sha256=spec['sha256'],
               K_generating=spec['K_generating'], tree_leaves=spec['tree_leaves'],
               original_leaves=case['original_leaves'], vertices=nv,
               n_edges=len(case['edges']),
               hidden=len(case.get('hidden_leaves') or []))
    return res


def job_panel(spec):
    edges, leaves, diag = panel_edges(ROOT / spec['file'])
    rng = random.Random(spec['seed'])
    leafset = set(leaves)
    if spec['mode'] == 'internal':
        pool = [i for i, (u, v, _) in enumerate(edges)
                if u not in leafset and v not in leafset]
    else:
        pool = list(range(len(edges)))
    support = sorted(rng.sample(pool, spec['K']))
    st = build_struct(edges, support, leaves)
    res = analyse(st, False, False)
    res.update(kind='panel', file=spec['file'], mode=spec['mode'], K=spec['K'],
               seed=spec['seed'], panel_leaves=len(leaves), vertices=len(st['adj']),
               n_edges=len(edges), newick_leaves=diag['n_leaves'])
    return res


def run(spec):
    t0 = time.perf_counter()
    try:
        res = job_bank(spec) if spec['kind'] == 'bank' else job_panel(spec)
        res['status'] = 'PASS'
    except Exception as exc:                                  # retained, never dropped
        res = dict(kind=spec['kind'], status='FAIL', spec=spec,
                   error_type=type(exc).__name__, error=str(exc),
                   traceback=traceback.format_exc())
    res['wall_seconds'] = time.perf_counter() - t0
    return res


def jobs():
    index = json.loads((ROOT / 'data/stage287_instances/bank/index.json').read_text())['cases']
    out = [dict(kind='bank', **x) for x in index if x['family'].startswith('exact_')]
    protocol = json.loads((ROOT / 'results/stage287_instance_bank_protocol.json').read_text())
    panels = [p for p in protocol['reference_panels'] if p.get('pdd_positivity_ok')]
    panels = sorted({p['leaves']: p for p in sorted(panels, key=lambda z: z['file'])}.values(),
                    key=lambda z: -z['leaves'])
    for p in panels:
        for K in PANEL_K:
            for mode in ('any', 'internal'):
                for rep in range(PANEL_REPLICATES):
                    out.append(dict(kind='panel', file=p['file'], K=K, mode=mode,
                                    seed=1000 * K + 7 * rep + (1 if mode == 'internal' else 0),
                                    leaves=p['leaves']))
    return out


def main():
    allj = jobs()
    nbank = sum(j['kind'] == 'bank' for j in allj)
    print('jobs: %d bank + %d panel' % (nbank, len(allj) - nbank), flush=True)
    # big panels last and single-file to keep peak memory bounded
    allj.sort(key=lambda j: (j['kind'] == 'panel', j.get('leaves', 0)))
    with OUT.open('x') as fh:
        with mp.Pool(WORKERS, maxtasksperchild=8) as pool:
            done = 0
            for res in pool.imap_unordered(run, allj, chunksize=1):
                fh.write(json.dumps(encode(res)) + '\n')
                fh.flush()
                done += 1
                if done % 100 == 0:
                    print('%d/%d' % (done, len(allj)), flush=True)
    print('written', OUT)


if __name__ == '__main__':
    main()
