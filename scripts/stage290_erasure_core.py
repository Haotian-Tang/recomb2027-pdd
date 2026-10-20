#!/usr/bin/env python3
"""Stage290 core: reference-erasure tolerance at scale, and the bridge from
known-support erasure tolerance to unknown-support identifiability.

Known support.  For a fixed source support F on a tree T with support-edge
endpoint set P and observed reference-leaf set R, the fixed-support leaf-report
matrix [ d(r,v) | 1 ; 1..1 | 0 ] has rank s + 1, where s is the minimum number
of vertices whose deletion separates P from R (equivalently, by Menger, the
maximum number of vertex-disjoint P-to-R paths).  The fibre dimension in
endpoint-mass coordinates is |P| - s, which is 2H - s when the support is a
matching (H = |F|).  The reference-erasure tolerance is cogirth - 1 with

    cogirth = min { #R-leaves still joined to P after deleting <= s-1 vertices }.

This module recomputes the cogirth with the Stage196 recurrence (states D, K00,
K01, K11 over a deletion-budget axis), in a form that scales to 10^5-vertex
trees, and validates it against the frozen Stage196 implementation and against
brute force.

Unknown support.  The analogue asks how many references may be erased before
more than one support of minimum size explains the remaining observations.
This module enumerates ALL minimum-size supports with the Stage25/42 endpoint
polytope oracle, optionally behind a batched singular-value prefilter whose
survivors are always confirmed in exact rational arithmetic.

No pre-existing project file is modified by this module.
"""
import os
for _k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_k] = '1'
import hashlib
import itertools
import json
import random
import sys
from fractions import Fraction as F
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
sys.dont_write_bytecode = True

from stage25_benchmark import geometry, fixed_support_polytope   # noqa: E402

PRIME = 2147483647            # 2**31 - 1, prime; products fit in int64
PRIME2 = 2147483629           # second prime for the modular rank cross-check
BIG = 1 << 40                 # DP sentinel for +infinity (costs are integers)
STATE_D, STATE_K00, STATE_K01, STATE_K11 = 0, 1, 2, 3
ACCEPTED_ROOT_STATES = (STATE_D, STATE_K00, STATE_K11)
# (left_state, right_state, target_state): transcribed from Stage196 merge_state
MERGE = ((0, 0, 0), (0, 1, 0), (0, 3, 0),
         (1, 0, 1), (2, 0, 2), (3, 0, 3),
         (1, 1, 1), (2, 2, 2), (2, 3, 3), (3, 2, 3), (3, 3, 3))


def sha256_path(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def encode(x):
    if isinstance(x, F):
        return str(x)
    if isinstance(x, dict):
        return {str(k): encode(v) for k, v in x.items()}
    if isinstance(x, (list, tuple, set)):
        return [encode(v) for v in x]
    if isinstance(x, np.integer):
        return int(x)
    if isinstance(x, np.floating):
        return float(x)
    return x


# --------------------------------------------------------------------------
# case loading
# --------------------------------------------------------------------------

def load_bank_case(path):
    """Stage287 bank case -> exact arithmetic (Stage44 loader convention)."""
    case = json.loads(Path(path).read_text())
    case['edges'] = [(u, v, F(length)) for u, v, length in case['edges']]
    case['observed'] = {int(r): F(value) for r, value in case['observed'].items()}
    if case.get('pendant') is not None:
        case['pendant'] = F(case['pendant'])
    for key in ('generating_sources', 'source_details'):
        for s in case.get(key) or []:
            s['p'] = F(s['p'])
            s['x'] = F(s['x'])
    return case


def support_of(case):
    return sorted({s['edge'] for s in case['generating_sources']})


# --------------------------------------------------------------------------
# structure
# --------------------------------------------------------------------------

def build_struct(edges, support_edges, observed_ids):
    adj = {}
    for i, (u, v, L) in enumerate(edges):
        assert u != v and L > 0
        adj.setdefault(u, []).append((v, i))
        adj.setdefault(v, []).append((u, i))
    assert len(edges) == len(adj) - 1, 'edge count must be |V| - 1'
    P = set()
    for e in support_edges:
        P.add(edges[e][0])
        P.add(edges[e][1])
    leaves = {v for v in adj if len(adj[v]) == 1}
    R = set(observed_ids)
    assert R <= set(adj), 'observed ids must be vertices'
    return dict(adj=adj, edges=edges, support=list(support_edges), P=P, R=R,
                leaves=leaves, H=len(support_edges),
                matching=len(P) == 2 * len(support_edges),
                pendant_support=sorted(e for e in support_edges
                                       if edges[e][0] in leaves or edges[e][1] in leaves))


def root_order(adj, r0):
    parent = {r0: None}
    order = [r0]
    for u in order:
        for w, _ in adj[u]:
            if w not in parent:
                parent[w] = u
                order.append(w)
    assert len(order) == len(adj), 'input must be connected'
    return parent, order


def reduced_dp_tree(st):
    """Collapse every maximal P-free hanging branch into one weighted leaf.

    Inside a P-free branch only two options are non-dominated: keep the branch
    (budget 0, charge its observed-leaf count) or cut its root (budget 1,
    charge 0); any deeper cut spends at least as much budget and retains at
    least as many leaves.  A weighted synthetic leaf has exactly those two
    options, so the reduction is exact.
    """
    adj, P, R = st['adj'], st['P'], st['R']
    r0 = min(P)
    parent, order = root_order(adj, r0)
    cntP = dict.fromkeys(order, 0)
    cntR = dict.fromkeys(order, 0)
    for v in reversed(order):
        cntP[v] += int(v in P)
        cntR[v] += int(v in R)
        p = parent[v]
        if p is not None:
            cntP[p] += cntP[v]
            cntR[p] += cntR[v]
    skel = [v for v in order if cntP[v] > 0]
    inskel = set(skel)
    index = {v: i for i, v in enumerate(skel)}
    children = [[] for _ in skel]
    isP = [v in P for v in skel]
    weight = [int(v in R) for v in skel]
    real = list(skel)
    for v in skel[1:]:
        children[index[parent[v]]].append(index[v])
    n_skel = len(skel)
    for v in skel:
        for w, _ in adj[v]:
            if w == parent[v] or w in inskel or cntR[w] == 0:
                continue
            j = len(children)
            children.append([])
            isP.append(False)
            weight.append(cntR[w])
            real.append(w)
            children[index[v]].append(j)
    return dict(children=children, isP=isP, weight=weight, real=real,
                order=list(range(len(children))), n_skeleton=n_skel,
                n_nodes=len(children), root=0)


def dp_profile(tree, Kmax, want_witness=True):
    """Stage196 recurrence on the (reduced) tree.

    f[b] = minimum number of observed leaves still joined to P when exactly b
    vertices are deleted (BIG if no deletion set of that size exists).
    """
    children, isP, weight = tree['children'], tree['isP'], tree['weight']
    K = Kmax
    tab = [None] * tree['n_nodes']
    choice = [None] * tree['n_nodes'] if want_witness else None
    for v in reversed(tree['order']):
        val = np.full((4, K + 1), BIG, dtype=np.int64)
        if K >= 1:
            val[STATE_D, 1] = 0
        if isP[v]:
            val[STATE_K11, 0] = weight[v]
        else:
            val[STATE_K00, 0] = 0
            val[STATE_K01, 0] = weight[v]
        steps = []
        for ch in children[v]:
            right = tab[ch]
            new = np.full((4, K + 1), BIG, dtype=np.int64)
            arg = np.full((4, K + 1, 3), -1, dtype=np.int32) if want_witness else None
            for s1, s2, tg in MERGE:
                A, B = val[s1], right[s2]
                for b1 in range(K + 1):
                    a = int(A[b1])
                    if a >= BIG:
                        continue
                    cand = B[:K + 1 - b1] + a
                    tgt = new[tg, b1:]
                    if want_witness:
                        better = np.nonzero(cand < tgt)[0]
                        if better.size:
                            tgt[better] = cand[better]
                            arg[tg, b1 + better, 0] = s1
                            arg[tg, b1 + better, 1] = s2
                            arg[tg, b1 + better, 2] = b1
                    else:
                        np.minimum(tgt, cand, out=tgt)
            val = new
            if want_witness:
                steps.append(arg)
            tab[ch] = None
        tab[v] = val
        if want_witness:
            choice[v] = steps
    final = tab[tree['root']]
    f = np.min(final[list(ACCEPTED_ROOT_STATES), :], axis=0)
    best_state = {b: min((int(final[s, b]), s) for s in ACCEPTED_ROOT_STATES)[1]
                  for b in range(K + 1)}
    return f, best_state, choice


def witness_separator(tree, choice, budget, state):
    """Replay the argmin trace; returns the deleted ORIGINAL vertices."""
    out = set()
    stack = [(tree['root'], len(tree['children'][tree['root']]) - 1, budget, state)]
    while stack:
        v, step, b, s = stack.pop()
        if step < 0:
            if s == STATE_D:
                out.add(tree['real'][v])
            continue
        s1, s2, b1 = (int(x) for x in choice[v][step][s, b])
        assert s1 >= 0, 'traceback entered an infeasible state'
        stack.append((v, step - 1, b1, s1))
        ch = tree['children'][v][step]
        stack.append((ch, len(tree['children'][ch]) - 1, b - b1, s2))
    return out


def components_objective(st, C):
    """Independent recomputation on the FULL tree (Stage196 components check)."""
    adj, P, R = st['adj'], st['P'], st['R']
    remaining = set(adj) - set(C)
    cost, D = 0, []
    while remaining:
        start = next(iter(remaining))
        remaining.discard(start)
        todo, comp = [start], [start]
        while todo:
            u = todo.pop()
            for w, _ in adj[u]:
                if w in remaining:
                    remaining.discard(w)
                    todo.append(w)
                    comp.append(w)
        cs = set(comp)
        if cs & P:
            leaves = sorted(cs & R)
            D.extend(leaves)
            cost += len(leaves)
    return cost, sorted(D)


def erasure_analysis(st, want_witness=True):
    """s (min P-to-R vertex cut), cogirth delta, erasure tolerance delta - 1."""
    tree = reduced_dp_tree(st)
    Kmax = len(st['P'])
    f, best_state, choice = dp_profile(tree, Kmax, want_witness=want_witness)
    zeros = np.nonzero(f == 0)[0]
    assert zeros.size, 'deleting all of P must disconnect P from R'
    s = int(zeros[0])
    out = dict(s=s, H=st['H'], P_size=len(st['P']), matching=st['matching'],
               R_size=len(st['R']),
               profile=[int(x) if x < BIG else None for x in f],
               fibre_dim_endpoint=len(st['P']) - s,
               fibre_dim_2H=2 * st['H'] - s,
               reduced_nodes=tree['n_nodes'], skeleton_nodes=tree['n_skeleton'])
    if s == 0:
        out.update(delta=0, tolerance=-1, witness_C=[], witness_D=[],
                   witness_checked=False)
        return out
    b_star = int(np.argmin(f[:s]))
    delta = int(f[b_star])
    out.update(delta=delta, tolerance=delta - 1)
    if want_witness:
        C = witness_separator(tree, choice, b_star, best_state[b_star])
        cost, D = components_objective(st, C)
        out.update(witness_C=sorted(C), witness_D=D, witness_cost=cost,
                   witness_checked=bool(cost == delta and len(C) == b_star))
    else:
        out.update(witness_checked=False)
    return out


# --------------------------------------------------------------------------
# exact / modular rank of the fixed-support leaf-report matrix
# --------------------------------------------------------------------------

def frac_mod(q, p):
    return (q.numerator % p) * pow(q.denominator % p, p - 2, p) % p


def distances_mod(st, targets, p):
    """BFS distance (mod p) from each target vertex to every vertex."""
    adj, edges = st['adj'], st['edges']
    elen = [frac_mod(L, p) for _, _, L in edges]
    verts = sorted(adj)
    index = {v: i for i, v in enumerate(verts)}
    out = np.zeros((len(targets), len(verts)), dtype=np.int64)
    for t, src in enumerate(targets):
        d = {src: 0}
        order = [src]
        for u in order:
            du = d[u]
            for w, ei in adj[u]:
                if w not in d:
                    d[w] = (du + elen[ei]) % p
                    order.append(w)
        row = out[t]
        for v, i in index.items():
            row[i] = d[v]
    return out, index


def rank_mod(M, p):
    """Rank of an integer matrix over GF(p) (p < 2**31)."""
    A = (M % p).astype(np.int64, copy=True)
    m, k = A.shape
    row = 0
    for col in range(k):
        if row >= m:
            break
        nz = np.nonzero(A[row:, col])[0]
        if nz.size == 0:
            continue
        piv = row + int(nz[0])
        if piv != row:
            A[[row, piv]] = A[[piv, row]]
        inv = pow(int(A[row, col]), p - 2, p)
        A[row] = (A[row] * inv) % p
        colvals = A[row + 1:, col]
        hit = np.nonzero(colvals)[0]
        if hit.size:
            A[row + 1:][hit] = (A[row + 1:][hit] -
                                colvals[hit, None] * A[row][None, :]) % p
        row += 1
    return row


def fixed_support_rank(st, p=PRIME):
    """rank of [ d(r,v) | 1 ; 1..1 | 0 ] over GF(p).  The law predicts s + 1."""
    Plist, Rlist = sorted(st['P']), sorted(st['R'])
    dmod, index = distances_mod(st, Plist, p)
    M = np.empty((len(Rlist) + 1, len(Plist) + 1), dtype=np.int64)
    cols = np.array([index[r] for r in Rlist], dtype=np.int64)
    for j in range(len(Plist)):
        M[:len(Rlist), j] = dmod[j, cols]
    M[:len(Rlist), len(Plist)] = 1
    M[len(Rlist), :len(Plist)] = 1
    M[len(Rlist), len(Plist)] = 0
    return rank_mod(M, p)


def full_dp_tree(st):
    """Unreduced DP tree over every vertex (cross-check for the reduction)."""
    adj, P, R = st['adj'], st['P'], st['R']
    r0 = min(P)
    parent, order = root_order(adj, r0)
    index = {v: i for i, v in enumerate(order)}
    children = [[] for _ in order]
    for v in order[1:]:
        children[index[parent[v]]].append(index[v])
    return dict(children=children, isP=[v in P for v in order],
                weight=[int(v in R) for v in order], real=list(order),
                order=list(range(len(order))), n_skeleton=len(order),
                n_nodes=len(order), root=0)


def erasure_analysis_full(st):
    """Same quantities as erasure_analysis but without the branch collapse."""
    tree = full_dp_tree(st)
    f, _, _ = dp_profile(tree, len(st['P']), want_witness=False)
    zeros = np.nonzero(f == 0)[0]
    s = int(zeros[0])
    delta = int(np.min(f[:s])) if s > 0 else 0
    return dict(s=s, delta=delta, tolerance=delta - 1, n_nodes=tree['n_nodes'])


def brute_force_cogirth(st, limit=40000, p=PRIME):
    """Independent cogirth by exhaustive erasure of observed-leaf subsets.

    The cogirth is the least |X| with rank(rows R minus X) < rank(rows R); this
    enumerates X by increasing size using modular rank, which shares no code
    path with the budgeted dynamic program.
    """
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
    full = rank_mod(M, p)
    tested = 0
    for size in range(1, nR + 1):
        total = 1
        for j in range(size):
            total = total * (nR - j) // (j + 1)
        if tested + total > limit:
            return dict(status='limit', cogirth=None, tested=tested,
                        full_rank=full, next_size=size)
        for X in itertools.combinations(range(nR), size):
            tested += 1
            keep = [i for i in range(nR) if i not in X] + [nR]
            if rank_mod(M[keep], p) < full:
                return dict(status='complete', cogirth=size, tested=tested,
                            full_rank=full, witness=[Rlist[i] for i in X])
    return dict(status='complete', cogirth=None, tested=tested, full_rank=full)


# --------------------------------------------------------------------------
# unknown support: enumeration of ALL minimum-size supports
# --------------------------------------------------------------------------

def endpoints_of(edges, support):
    return tuple(sorted({v for i in support for v in edges[i][:2]}))


class CaseTables(object):
    """Per-case precomputation shared by every erasure level of that case."""

    def __init__(self, edges):
        self.edges = edges
        _, leaves, dist = geometry(edges)
        self.dist = dist
        self.leaves = leaves
        verts = sorted({v for u, w, _ in edges for v in (u, w)})
        self.vindex = {v: i for i, v in enumerate(verts)}
        self.lindex = {r: i for i, r in enumerate(leaves)}
        D = np.empty((len(leaves), len(verts)), dtype=np.float64)
        for i, r in enumerate(leaves):
            row = dist[r]
            for v, j in self.vindex.items():
                D[i, j] = float(row[v])
        self.D = D


def min_supports(tables, observed, cap, rng_seed=0, extra_rows=4, tol=1e-7,
                 exhaustive_exact=False, max_candidates=4_000_000):
    """All minimum-size supports admitting a strictly positive witness.

    With `exhaustive_exact` every candidate goes to the exact rational oracle
    (ground truth).  Otherwise a batched singular-value prefilter preselects
    candidates; consistency over Q forces the augmented matrix to be column
    rank deficient, so the prefilter can only err by a floating-point false
    negative, whose margin is reported.  Every survivor is confirmed exactly.
    """
    edges = tables.edges
    Rlist = sorted(observed)
    rowsel = np.array([tables.lindex[r] for r in Rlist], dtype=np.int64)
    yv = np.array([float(observed[r]) for r in Rlist], dtype=np.float64)
    rng = random.Random(rng_seed)
    nE = len(edges)
    stats = dict(supports_visited=0, endpoint_sets_solved=0, endpoint_cache_hits=0,
                 equality_rows_processed=0, vertex_bases_attempted=0,
                 prefilter_candidates=0, prefilter_survivors=0)
    margin = None
    cache = {}
    for k in range(1, cap + 1):
        total = 1
        for j in range(k):
            total = total * (nE - j) // (j + 1)
        if total > max_candidates:
            return dict(status='candidate_limit', minimum=None, supports=None,
                        k_aborted=k, candidates=int(total), statistics=stats)
        supports = list(itertools.combinations(range(nE), k))
        stats['supports_visited'] += len(supports)
        eps = [endpoints_of(edges, s) for s in supports]
        uniq = sorted(set(eps))
        if exhaustive_exact or len(Rlist) < 2 * k + 6:
            keep_set = set(uniq)
        else:
            keep_set = set()
            by_q = {}
            for ep in uniq:
                by_q.setdefault(len(ep), []).append(ep)
            for q, group in sorted(by_q.items()):
                nrows = min(len(Rlist), max(2 * q + 4, q + 10))
                if nrows < q + 2:
                    keep_set.update(group)
                    continue
                pick = np.array(sorted(rng.sample(range(len(Rlist)), nrows)))
                # Centring over the sampled rows removes the all-ones pendant
                # column; column scaling removes the severe ill-conditioning of
                # phylogenetic distances.  The test is whether the centred
                # observation vector lies in the column space of the centred
                # endpoint-distance columns, which it must if the support is
                # feasible, so the prefilter has no false negatives beyond
                # floating-point noise and its margin is reported.
                Dsel = tables.D[rowsel[pick], :]
                Dc = Dsel - Dsel.mean(axis=0, keepdims=True)
                scale = np.linalg.norm(Dc, axis=0)
                scale[scale == 0] = 1.0
                Dn = Dc / scale
                yc = yv[pick] - yv[pick].mean()
                ny = np.linalg.norm(yc)
                yn = yc / (ny if ny > 0 else 1.0)
                for start_i in range(0, len(group), 20000):
                    part = group[start_i:start_i + 20000]
                    C = np.empty((len(part), nrows, q), dtype=np.float64)
                    for i, ep in enumerate(part):
                        for j, v in enumerate(ep):
                            C[i, :, j] = Dn[:, tables.vindex[v]]
                    U, S, _ = np.linalg.svd(C, full_matrices=False)
                    live = S > (1e-12 * np.maximum(S[:, :1], 1e-300))
                    proj = np.einsum('nrj,r->nj', U, yn) * live
                    resid = np.sqrt(np.maximum(0.0, 1.0 - np.sum(proj ** 2, axis=1)))
                    sel = resid < tol
                    for i, ep in enumerate(part):
                        if sel[i]:
                            keep_set.add(ep)
                    if np.any(~sel):
                        w = float(np.min(resid[~sel]))
                        margin = w if margin is None else min(margin, w)
            stats['prefilter_candidates'] += len(uniq)
            stats['prefilter_survivors'] += len(keep_set)
        found = []
        for sup, ep in zip(supports, eps):
            if ep not in keep_set:
                continue
            if ep not in cache:
                stats['endpoint_sets_solved'] += 1
                cache[ep] = fixed_support_polytope(ep, Rlist, tables.dist,
                                                   observed, stats)
            else:
                stats['endpoint_cache_hits'] += 1
            if cache[ep] is not None:
                found.append(list(sup))
        if found:
            return dict(status='complete', minimum=k, supports=found,
                        n_supports=len(found), statistics=stats,
                        prefilter_margin=margin)
    return dict(status='complete', minimum=None, supports=[], n_supports=0,
                statistics=stats, prefilter_margin=margin)


# --------------------------------------------------------------------------
# removal regimes
# --------------------------------------------------------------------------

def source_point_distance(tables, case, r):
    best = None
    for s in case['generating_sources']:
        u, v, L = case['edges'][s['edge']]
        d = min(tables.dist[r][u] + s['x'] * L,
                tables.dist[r][v] + (1 - s['x']) * L)
        best = d if best is None or d < best else best
    return best


def removal_order(case, st, regime, seed, tables=None):
    """A nested removal sequence over the observed reference leaves."""
    Rlist = sorted(st['R'])
    rng = random.Random(seed)
    if regime == 'random':
        order = Rlist[:]
        rng.shuffle(order)
        return order
    if regime == 'adversarial_nearest':
        if tables is None:
            tables = CaseTables(case['edges'])
        return sorted(Rlist, key=lambda r: (source_point_distance(tables, case, r), r))
    if regime == 'whole_clade':
        adj = st['adj']
        Rset = set(Rlist)
        r0 = min(adj)
        parent, order = root_order(adj, r0)
        sub = {}
        for v in reversed(order):
            acc = [v] if v in Rset else []
            for w, _ in adj[v]:
                if parent.get(w) == v:
                    acc.extend(sub[w])
            sub[v] = acc
        clades = [sorted(sub[v]) for v in order
                  if v != r0 and 1 <= len(sub[v]) <= max(1, len(Rlist) // 2)]
        rng.shuffle(clades)
        seen, out = set(), []
        for cl in clades:
            for r in cl:
                if r not in seen:
                    seen.add(r)
                    out.append(r)
        for r in Rlist:
            if r not in seen:
                seen.add(r)
                out.append(r)
        return out
    raise ValueError('unknown regime ' + regime)
