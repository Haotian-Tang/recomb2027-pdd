#!/usr/bin/env python3
"""Stage292 validation of the all-minimisers routine.

Track A  replay: every probe point recorded in results/stage290_bridge.jsonl
         (130 cases x 5 regimes x the probes of its bisection, plus e = 0) is
         recomputed with the new structural routine AND re-run with the
         Stage290 brute force, and the three answers are compared.
Track B  fresh small instances generated here, where the complete C(|E|,k)
         enumeration against the frozen Stage25 endpoint oracle is affordable,
         compared support-set for support-set; the reduced zero-face
         enumeration is additionally compared against the full 2^|A| one.

Nothing is patched to agree: every disagreement is written out in full.
"""
import json
import multiprocessing as mp
import os
import random
import sys
import time
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
sys.dont_write_bytecode = True

from stage25_benchmark import geometry, fixed_support_polytope          # noqa: E402
from stage26_missing_leaf import observe                                # noqa: E402
from stage290_erasure_core import (load_bank_case, support_of, CaseTables,  # noqa: E402
                                   min_supports, removal_order, build_struct)
from stage292_all_minimisers import all_minimum_supports                # noqa: E402
from stage292_harness import encode, isolated_run                       # noqa: E402

WORKERS = int(os.environ.get('STAGE292_WORKERS', '28'))
OUT_A = ROOT / 'results/stage292_validate_replay.jsonl'
OUT_B = ROOT / 'results/stage292_validate_fresh.jsonl'
MEMORY = 6 << 30
WALL_A = 5400.0
WALL_B = 900.0
EXACT_BRUTE_EDGES = 70        # exhaustive_exact brute force below this size


# ---------------------------------------------------------------------------
# track A
# ---------------------------------------------------------------------------

def replay_case(spec):
    rec = spec['record']
    case = load_bank_case(ROOT / rec['path'])
    support = support_of(case)
    K = len(support)
    observed = case['observed']
    tables = CaseTables(case['edges'])
    st = build_struct(case['edges'], support, sorted(observed))
    exact_ok = len(case['edges']) <= EXACT_BRUTE_EDGES
    probes = [dict(regime='baseline', seed=None, e=0)]
    for reg in rec['regimes']:
        for pr in reg['probe_table']:
            probes.append(dict(regime=reg['regime'], seed=reg['seed'], e=pr['e'],
                               recorded_minimum=pr['minimum'],
                               recorded_n_supports=pr['n_supports'],
                               recorded_ok=pr['ok']))
    orders = {}
    out = []
    for pr in probes:
        key = (pr['regime'], pr['seed'])
        if pr['e'] == 0:
            keys = sorted(observed)
        else:
            if key not in orders:
                orders[key] = removal_order(case, st, pr['regime'], pr['seed'],
                                            tables=tables)
            removed = set(orders[key][:pr['e']])
            keys = [r for r in sorted(observed) if r not in removed]
        sub = {r: observed[r] for r in keys}
        t0 = time.perf_counter()
        new = all_minimum_supports(case['edges'], sub, cap=K)
        t1 = time.perf_counter()
        bf = min_supports(tables, sub, K, rng_seed=7)
        t2 = time.perf_counter()
        ex = min_supports(tables, sub, K, rng_seed=7, exhaustive_exact=True) \
            if exact_ok else None
        t3 = time.perf_counter()
        new_set = sorted(map(tuple, new['supports']))
        bf_set = sorted(map(tuple, bf['supports'] or []))
        row = dict(name=rec['name'], family=rec['family'], K=K,
                   n_edges=len(case['edges']), R_size=len(observed),
                   regime=pr['regime'], seed=pr['seed'], e=pr['e'],
                   R_remaining=len(keys),
                   new_status=new['status'], new_minimum=new['minimum'],
                   new_n=new['n_supports'],
                   bf_status=bf['status'], bf_minimum=bf['minimum'],
                   bf_n=bf['n_supports'],
                   agree_minimum=bool(new['minimum'] == bf['minimum']),
                   agree_supports=bool(new_set == bf_set),
                   recorded_minimum=pr.get('recorded_minimum'),
                   recorded_n_supports=pr.get('recorded_n_supports'),
                   new_seconds=t1 - t0, bf_seconds=t2 - t1,
                   exact_seconds=(t3 - t2) if exact_ok else None,
                   frozen_minimum=new.get('frozen_minimum'),
                   oracle_calls=new['statistics']['oracle_calls'],
                   bf_supports_visited=bf['statistics']['supports_visited'],
                   configurations=new['statistics']['configurations'])
        if ex is not None:
            ex_set = sorted(map(tuple, ex['supports'] or []))
            row.update(exact_minimum=ex['minimum'], exact_n=ex['n_supports'],
                       agree_exact_minimum=bool(new['minimum'] == ex['minimum']),
                       agree_exact_supports=bool(new_set == ex_set))
        if pr.get('recorded_minimum') is not None:
            row['agree_recorded'] = bool(
                new['minimum'] == pr['recorded_minimum'] and
                new['n_supports'] == pr['recorded_n_supports'])
        if not (row['agree_minimum'] and row['agree_supports']) or \
                row.get('agree_recorded') is False or \
                row.get('agree_exact_supports') is False:
            row['disagreement'] = dict(new=new_set, brute=bf_set,
                                       true_support=support,
                                       new_full=encode(new))
        out.append(row)
    return out


# ---------------------------------------------------------------------------
# track B
# ---------------------------------------------------------------------------

def random_instance(seed):
    rng = random.Random(seed)
    nleaf = rng.randint(7, 12)
    edges = [[0, nleaf, F(rng.randint(1, 30), rng.choice([1, 2, 3, 5, 7]))],
             [1, nleaf, F(rng.randint(1, 30), rng.choice([1, 2, 3, 5, 7]))],
             [2, nleaf, F(rng.randint(1, 30), rng.choice([1, 2, 3, 5, 7]))]]
    nextv = nleaf + 1
    internal = [nleaf]
    for leaf in range(3, nleaf):
        if rng.random() < 0.78:
            j = rng.randrange(len(edges))
            u, v, L = edges[j]
            a = L * F(rng.randint(1, 9), 10)
            w = nextv
            nextv += 1
            edges[j] = [u, w, a]
            edges.append([w, v, L - a])
            edges.append([leaf, w, F(rng.randint(1, 30), rng.choice([1, 2, 3, 5, 7]))])
            internal.append(w)
        else:
            w = rng.choice(internal)
            edges.append([leaf, w, F(rng.randint(1, 30), rng.choice([1, 2, 3, 5, 7]))])
    edges = [tuple(e) for e in edges]
    K = rng.choice([1, 2, 2, 3, 3, 4])
    chosen = rng.sample(range(len(edges)), K)
    raw = [rng.randint(1, 9) for _ in chosen]
    total = sum(raw)
    sources = [{'edge': i, 'p': F(w, total), 'x': F(rng.randint(1, 11), 12)}
               for i, w in zip(chosen, raw)]
    y = F(rng.randint(1, 40), rng.choice([1, 2, 4]))
    leaves = sorted(v for v in range(nleaf))
    nobs = rng.randint(3, nleaf)
    R = sorted(rng.sample(leaves, nobs))
    observed = observe(edges, sources, y, R)
    return dict(edges=edges, observed=observed, sources=sources, pendant=y,
                support=sorted(chosen), K=K, nleaf=nleaf, R=R)


def brute_all(edges, observed, cap):
    """Independent exhaustive enumeration: every k-subset, exact oracle."""
    _, leaves, dist = geometry(edges)
    Rlist = sorted(observed)
    stats = dict(equality_rows_processed=0, vertex_bases_attempted=0)
    cache = {}
    visited = 0
    for k in range(1, cap + 1):
        found = []
        for sup in combinations(range(len(edges)), k):
            visited += 1
            ep = tuple(sorted({v for i in sup for v in edges[i][:2]}))
            if ep not in cache:
                cache[ep] = fixed_support_polytope(ep, Rlist, dist, observed,
                                                   stats) is not None
            if cache[ep]:
                found.append(list(sup))
        if found:
            return dict(minimum=k, supports=found, visited=visited)
    return dict(minimum=None, supports=[], visited=visited)


def fresh_case(spec):
    inst = random_instance(spec['seed'])
    cap = inst['K']
    t0 = time.perf_counter()
    new = all_minimum_supports(inst['edges'], inst['observed'], cap=cap)
    t1 = time.perf_counter()
    full = all_minimum_supports(inst['edges'], inst['observed'], cap=cap,
                                full_mask=True)
    t2 = time.perf_counter()
    bf = brute_all(inst['edges'], inst['observed'], cap)
    t3 = time.perf_counter()
    new_set = sorted(map(tuple, new['supports']))
    full_set = sorted(map(tuple, full['supports']))
    bf_set = sorted(map(tuple, bf['supports']))
    row = dict(seed=spec['seed'], K=inst['K'], nleaf=inst['nleaf'],
               n_edges=len(inst['edges']), R_size=len(inst['R']),
               hidden=inst['nleaf'] - len(inst['R']),
               true_support=inst['support'],
               new_status=new['status'], new_minimum=new['minimum'],
               new_n=new['n_supports'], brute_minimum=bf['minimum'],
               brute_n=len(bf['supports']), brute_visited=bf['visited'],
               agree_minimum=bool(new['minimum'] == bf['minimum']),
               agree_supports=bool(new_set == bf_set),
               agree_full_mask=bool(new_set == full_set and
                                    new['minimum'] == full['minimum']),
               oracle_calls=new['statistics']['oracle_calls'],
               configurations=new['statistics']['configurations'],
               new_seconds=t1 - t0, full_seconds=t2 - t1, brute_seconds=t3 - t2)
    if not (row['agree_minimum'] and row['agree_supports'] and row['agree_full_mask']):
        row['disagreement'] = dict(new=new_set, full=full_set, brute=bf_set,
                                   instance=encode(inst), new_full=encode(new))
    return [row]


# ---------------------------------------------------------------------------

def run(spec):
    fn = replay_case if spec['track'] == 'A' else fresh_case
    wall = WALL_A if spec['track'] == 'A' else WALL_B
    res = isolated_run(fn, spec, wall, MEMORY)
    res['spec'] = {k: v for k, v in spec.items() if k != 'record'}
    if spec['track'] == 'A':
        res['spec']['name'] = spec['record']['name']
    return res


def main():
    track = sys.argv[1] if len(sys.argv) > 1 else 'AB'
    if 'A' in track:
        records = [json.loads(l) for l in (ROOT / 'results/stage290_bridge.jsonl').open()]
        jobs = [dict(track='A', record=r) for r in records if r['status'] == 'PASS']
        jobs.sort(key=lambda j: -j['record']['n_edges'])
        _execute(jobs, OUT_A)
    if 'B' in track:
        jobs = [dict(track='B', seed=s) for s in range(1, 501)]
        _execute(jobs, OUT_B)


def _execute(jobs, out):
    print('%d jobs -> %s' % (len(jobs), out), flush=True)
    with out.open('x') as fh:
        with mp.Pool(WORKERS, maxtasksperchild=4) as pool:
            done = 0
            for res in pool.imap_unordered(run, jobs, chunksize=1):
                if res['status'] == 'ok':
                    for row in res['value']:
                        fh.write(json.dumps(encode(row)) + '\n')
                else:
                    fh.write(json.dumps(encode(dict(row_status=res['status'],
                                                    **res))) + '\n')
                fh.flush()
                done += 1
                if done % 20 == 0 or done == len(jobs):
                    print('%d/%d' % (done, len(jobs)), flush=True)
    print('written', out, flush=True)


if __name__ == '__main__':
    main()

