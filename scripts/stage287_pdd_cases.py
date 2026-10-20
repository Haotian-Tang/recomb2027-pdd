#!/usr/bin/env python3
"""Stage 287 -- conversion of public mixture instances into the project's exact
PDD case format.

Generative model of the parent work (DecoDiPhy ``find_sol.py`` /
``find_sol_10k.py``): a mixture is simulated by choosing K leaves of a FULL
reference phylogeny as query taxa, drawing abundances p, and reporting, for
every remaining (reference) leaf r,

    d_r = sum_i p_i * dist_{T_full}(q_i, r).

Pruning the query leaves turns each q_i into a point source located where its
neighbour a_i used to sit, at parameter x_i inside an edge of the pruned
reference tree T, with pendant y_i = length of the edge (q_i, a_i); the common
mean pendant is y = sum_i p_i * y_i. Hence

    d_r = y + sum_i p_i * dist_T(r, q_i),

which is the project's PDD model. This module performs that conversion in exact
``Fraction`` arithmetic and verifies it: the full-tree computation and the
pruned-tree edge-parameterised computation are independent code paths and must
agree exactly.
"""
from fractions import Fraction
import collections

ORIGINAL_LEAF = 'leaf'


def undirected(tree):
    """Adjacency dict {v: {u: Fraction}} from a parsed Newick tree.

    The Newick root is suppressed when it has degree 2, so the result is the
    unrooted tree the PDD model works with.
    """
    adj = collections.defaultdict(dict)
    for i in range(tree.n):
        if i == tree.root:
            continue
        p = tree.parent[i]
        L = tree.elen[i]
        if L is None:
            raise ValueError('missing edge length at vertex %d' % i)
        adj[i][p] = L
        adj[p][i] = L
    for i in range(tree.n):
        adj.setdefault(i, {})
    r = tree.root
    if len(adj[r]) == 2:
        (a, la), (b, lb) = list(adj[r].items())
        del adj[a][r]
        del adj[b][r]
        adj[a][b] = la + lb
        adj[b][a] = la + lb
        del adj[r]
    return dict(adj)


def leaf_labels(tree):
    return {i: tree.label[i] for i in range(tree.n) if not tree.children[i]}


def dists_from(adj, src):
    """Exact distances from ``src`` to every vertex (iterative DFS)."""
    d = {src: Fraction(0)}
    stack = [src]
    while stack:
        u = stack.pop()
        du = d[u]
        for v, L in adj[u].items():
            if v not in d:
                d[v] = du + L
                stack.append(v)
    return d


def _reduce(adj, removed, keep_leaves):
    """Remove ``removed`` vertices, then iteratively drop dangling vertices."""
    radj = {v: dict(e) for v, e in adj.items() if v not in removed}
    for v in radj:
        for u in list(radj[v]):
            if u in removed:
                del radj[v][u]
    changed = True
    while changed:
        changed = False
        for v in list(radj):
            if v in keep_leaves:
                continue
            if len(radj[v]) <= 1:
                for u in list(radj[v]):
                    del radj[u][v]
                del radj[v]
                changed = True
    return radj


def contract(radj, keep_leaves):
    """Contract degree-2 chains.

    Returns (edges, offset) where ``edges`` is a list of
    (u, v, length, [(chain_vertex, distance_from_u), ...]) over the kept
    vertices, and ``offset`` maps each contracted degree-2 vertex to
    (edge_index, distance_from_u).
    """
    kept = {v for v in radj if v in keep_leaves or len(radj[v]) != 2}
    edges = []
    offset = {}
    seen = set()
    for u in sorted(kept):
        for first in sorted(radj[u]):
            prev, cur = u, first
            length = radj[u][first]
            chain = []
            while cur not in kept:
                chain.append((cur, length))
                nxt = [w for w in radj[cur] if w != prev]
                if len(nxt) != 1:
                    raise ValueError('chain walk failed at %d' % cur)
                prev, cur = cur, nxt[0]
                length = length + radj[prev][cur]
            key = (min(u, cur), max(u, cur))
            if key in seen:
                continue
            seen.add(key)
            idx = len(edges)
            edges.append((u, cur, length, chain))
            for cv, doff in chain:
                offset[cv] = (idx, doff)
    return edges, offset


class CaseError(Exception):
    pass


def build_exact_case(tree, query_p, cap=None, name=None):
    """Convert (full tree, {query label: Fraction p}) into an exact PDD case.

    Raises :class:`CaseError` with a machine-readable reason when the instance
    cannot be put in the model's normal form (sources strictly inside distinct
    original edges of the reference tree, all lengths positive).
    """
    lab = leaf_labels(tree)
    by_label = {}
    for v, l in lab.items():
        by_label.setdefault(l, []).append(v)
    dup = sorted(l for l, vs in by_label.items() if len(vs) > 1)
    if dup:
        raise CaseError('duplicate_leaf_labels:%d' % len(dup))
    label_to_v = {l: vs[0] for l, vs in by_label.items()}
    missing = [l for l in query_p if l not in label_to_v]
    if missing:
        raise CaseError('query_not_in_tree:%d' % len(missing))

    adj = undirected(tree)
    allv = set(adj)
    leaves = {v for v in allv if len(adj[v]) == 1}
    Q = [label_to_v[l] for l in sorted(query_p)]
    qlab = {label_to_v[l]: l for l in query_p}
    if any(q not in leaves for q in Q):
        raise CaseError('query_not_a_leaf')
    R = sorted(leaves - set(Q), key=lambda v: lab[v])
    if not R:
        raise CaseError('no_reference_leaves')

    # --- path A: observed distances from the FULL tree -------------------
    p = {q: Fraction(query_p[qlab[q]]) for q in Q}
    if any(v <= 0 for v in p.values()):
        raise CaseError('non_positive_weight')
    if sum(p.values()) != 1:
        raise CaseError('weights_do_not_sum_to_one')
    obs_full = {r: Fraction(0) for r in R}
    pend = {}
    for q in Q:
        dq = dists_from(adj, q)
        for r in R:
            if r not in dq:
                raise CaseError('disconnected_tree')
            obs_full[r] += p[q] * dq[r]
        (a_q, y_q), = adj[q].items()
        pend[q] = (a_q, y_q)
    y = sum(p[q] * pend[q][1] for q in Q)

    # --- reference tree ---------------------------------------------------
    radj = _reduce(adj, set(Q), set(R))
    if set(R) - set(radj):
        raise CaseError('reference_leaf_lost_in_reduction')
    edges, offset = contract(radj, set(R))
    anchors = {}
    for q in Q:
        a_q = pend[q][0]
        if a_q not in offset:
            if a_q in radj:
                raise CaseError('source_at_vertex_not_inside_edge')
            raise CaseError('source_vertex_removed')
        anchors[q] = offset[a_q]
    if len({e for e, _ in anchors.values()}) != len(Q):
        raise CaseError('sources_share_an_edge')

    for (u, v, L, _ch) in edges:
        if L <= 0:
            raise CaseError('non_positive_edge_length')
    for q, (ei, doff) in anchors.items():
        L = edges[ei][2]
        if not (0 < doff < L):
            raise CaseError('source_not_strictly_interior')
    if y <= 0:
        raise CaseError('non_positive_pendant')

    # --- integer ids ------------------------------------------------------
    vid = {}
    for i, r in enumerate(R):
        vid[r] = i
    nxt = len(R)
    for v in sorted(radj):
        if v not in vid:
            vid[v] = nxt
            nxt += 1
    out_edges = [[vid[u], vid[v], str(L)] for (u, v, L, _c) in edges]

    # --- path B: verify through the pruned-tree edge parameterisation -----
    endpoints = sorted({edges[ei][0] for ei, _ in anchors.values()} |
                       {edges[ei][1] for ei, _ in anchors.values()})
    dmap = {}
    tadj = collections.defaultdict(dict)
    for (u, v, L, _c) in edges:
        tadj[u][v] = L
        tadj[v][u] = L
    for ep in endpoints:
        dmap[ep] = dists_from(dict(tadj), ep)
    obs_edge = {r: Fraction(y) for r in R}
    srcs = []
    for q in Q:
        ei, doff = anchors[q]
        u, v, L, _c = edges[ei]
        x = doff / L
        srcs.append({'edge': ei, 'p': str(p[q]), 'x': str(x),
                     'query_label': qlab[q], 'pendant': str(pend[q][1])})
        for r in R:
            a = dmap[u][r] + x * L
            b = dmap[v][r] + (1 - x) * L
            obs_edge[r] += p[q] * (a if a <= b else b)
    mismatch = [lab[r] for r in R if obs_edge[r] != obs_full[r]]
    if mismatch:
        raise CaseError('exact_verification_failed:%d' % len(mismatch))

    case = {
        'name': name,
        'task': 'free_pendant',
        'edges': out_edges,
        'observed': {str(vid[r]): str(obs_full[r]) for r in R},
        'cap': int(cap if cap is not None else len(Q)),
        'hidden_leaves': [],
        'pendant': str(y),
        'generating_sources': [{'edge': s['edge'], 'p': s['p'], 'x': s['x']}
                               for s in srcs],
        'generating_pendant': str(y),
        'generating_count': len(Q),
        'original_leaves': len(R),
        'observed_leaves': len(R),
        'vertices': len(radj),
        'leaf_labels': {str(vid[r]): lab[r] for r in R},
        'source_details': srcs,
    }
    return case


def attach_published_observed(case, values, label_to_id):
    """Return a copy of ``case`` whose observed vector is a published table."""
    out = dict(case)
    obs = {}
    for lbl, val in values.items():
        vid = label_to_id.get(lbl)
        if vid is None:
            continue
        obs[vid] = str(Fraction(val))
    out['observed'] = obs
    out['observed_leaves'] = len(obs)
    out.pop('generating_pendant', None)
    return out
