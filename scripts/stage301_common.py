#!/usr/bin/env python3
"""Stage 301 shared library: transporting the frozen structural localisation
to the eps-tolerant source-budget problem.

SETTING.  Known reference tree T (positive rational lengths, reduced), observed
reference leaves R with exact rational d_r, k sources interior to k distinct
ORIGINAL edges, weights p_i > 0 summing to one, common mean pendant y > 0.
Stage 296 defines, for a fixed support S,

    theta(S) = min over the CLOSED set {alpha, beta, y >= 0, sum = 1}
               of max_r |pred_r - d_r|,      Theta_k = min_{|S| = k} theta(S),

and establishes m(eps) = min{k : Theta_k < eps} for eps > 0.  Stage 296 also
records that the frozen localisation (leaf envelope h, weighted deficit bound,
positive-anchor filters) is unsound for eps > 0 because the envelope bounds the
potential pointwise only under exact equality.

WHAT THIS MODULE ADDS.  Each frozen ingredient is transported to a TOLERANT
form that is sound for every eps >= 0.  Write the potential of a candidate
configuration as f(v) = y + sum_w m(w) d(v, w) (m >= 0, sum m = 1, supp(m)
inside the atom vertices of S) and phi = f - h with h(v) = max_r (d_r - d(v,r)).

  F-A  (Lipschitz)   |f(u) - f(v)| <= d(u, v), because sum_w m(w) = 1.
  F-B  (lower env.)  theta(S) <= eps  =>  phi >= -eps everywhere.
  F-C  (exposure)    v is EXPOSED if some branch at v carries no mass; then
                     f(v) = f(r) - d(v, r) for any observed r in that branch,
                     hence phi(v) <= eps.  Observed leaves also have phi <= eps.
  F-D  (deficit)     Sigma_eps := {v : phi(v) > eps} contains only surrounded
                     internal vertices, so 3 <= deg(v) <= k on Sigma_eps and
                     sum_{v in Sigma_eps} (deg(v) - 2) <= k - 2.  PURELY
                     STRUCTURAL: no observation enters, so it transports as is.
  F-E  (mass band)   for v outside N[Sigma_eps],
                     |m(v) - b(v)| <= tau_eps(v) := sum_{u~v} min(eps/l_uv, 1),
                     where b = inverse(h) is the baseline mass.  The cap at 1
                     uses F-A for f and the 1-Lipschitzness of h.
  F-F  (anchors)     every connected component of T[Sigma_eps] contains a
                     vertex with b(v) > 0 -- the UNRELAXED positive set.

Consequences used as screens, all necessary conditions for theta(S) <= eps:
  S1  NEG(eps) := {v : b(v) < -tau_eps(v)}  is contained in N[Sigma_eps];
  S2  POS(eps) := {v : b(v) >  tau_eps(v)}  satisfies POS \\ N[Sigma_eps]
      contained in supp(m), so the atoms of S must cover it and
      |POS(eps) \\ N[Sigma_eps]| <= 2k;
  S3  h(r) - d_r <= 2 eps at every observed leaf r;
  S4  Y >= -eps (|R| - 1) for the baseline pendant Y.

SECOND, INDEPENDENT TOLERANT LOCALISATION (the envelope in quantitative form).
For observed leaves r, s put Delta = d_r - d_s and let proj(v) be the point
where v meets the r-s path.  Then f(r) - f(s) = d(r,s) - 2 Psi with
Psi = sum_w m(w) d(s, proj(w)) >= 0, so theta(S) <= eps forces

    Psi <= rho(r, s, eps) := (d(r,s) - Delta)/2 + eps                    (S5)

and in particular SOME atom vertex of S has d(s, proj(v)) <= rho.  Choosing
r = argmax_r (d_r - d(r,s)) makes rho = (d_s - h(s))/2 + eps <= eps, so every
candidate support must contain an edge whose atom meets the subtree cut off at
path-distance rho from s.  This screen does not degrade as eps grows past the
edge lengths, which is exactly where S1/S2 go vacuous.

VALIDATION CONTRACT.  Nothing here decides a theta: every surviving support is
handed to the FROZEN stage296_common.theta_support, and every attaining witness
is re-checked with the FROZEN min(d(r,u)+xL, d(r,v)+(1-x)L) verifier.  With
screening disabled the candidate enumerator must visit exactly C(|E|, k)
supports.  Nothing pre-existing is modified.
"""
import os
for _k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_k] = '1'
import itertools
import sys
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / 'scripts') not in sys.path:
    sys.path.insert(0, str(ROOT / 'scripts'))

import stage288_common as common                      # noqa: E402  bank, bfs, pool
import stage296_common as s296                        # noqa: E402  frozen tolerant LP
from stage23_two_source_envelope import two_pass_envelope, inverse   # noqa: E402
from stage25_deficit_extension import tree            # noqa: E402
from stage26_missing_leaf import observed_backbone    # noqa: E402
from stage41_partial_fpt import deficit_sets          # noqa: E402
from stage27_general_fpt import Budget                # noqa: E402

OUT = ROOT / 'results/stage301'
ZERO, ONE = F(0), F(1)


# --------------------------------------------------------------- backbone view

def backbone_view(edges, observed):
    """Backbone of the observed leaves plus every quantity the screens need.

    `atoms[i]` is the frozenset of BACKBONE vertices at which original edge i
    can place mass: the two ends of its backbone edge when i lies on a chain,
    otherwise the gate vertex (or the ends of the backbone edge the gate is
    interior to).  Mass of any source on edge i lives in atoms[i], which is all
    the screens use about i.
    """
    edges = [(u, v, F(l)) for u, v, l in edges]
    obs = {int(r): F(d) for r, d in observed.items()}
    bb = observed_backbone(edges, sorted(obs))
    B = bb['edges']
    adj, _lookup = tree(B)
    leaves = sorted(v for v in adj if len(adj[v]) == 1)
    assert set(leaves) == set(obs), 'backbone leaves must be the observed set'
    h, _ = two_pass_envelope(adj, obs)
    b, Y = inverse(adj, h)
    o2b, gate, interior = bb['original_to_backbone'], bb['gate'], bb['interior_vertex']
    atoms = []
    for i, (u, v, _l) in enumerate(edges):
        if i in o2b:
            j = o2b[i]['backbone_edge']
            atoms.append(frozenset((B[j][0], B[j][1])))
        else:
            g = gate[u]
            if g in interior:
                j = interior[g]['edge']
                atoms.append(frozenset((B[j][0], B[j][1])))
            else:
                atoms.append(frozenset((g,)))
    return dict(edges=edges, backbone=B, adj=adj, leaves=leaves, obs=obs,
                h=h, b=b, Y=Y, atoms=atoms,
                deg={v: len(adj[v]) for v in adj},
                envelope_slack={r: h[r] - obs[r] for r in leaves})


def tau(view, v, eps):
    """Per-vertex tolerance band of F-E: sum_{u~v} min(eps / l_uv, 1)."""
    if eps == 0:
        return ZERO
    return sum(min(F(eps) / l, ONE) for l in view['adj'][v].values())


def sign_sets(view, eps):
    """POS(eps) and NEG(eps) of S1/S2, plus the unrelaxed anchor set of F-F."""
    b = view['b']
    pos, neg, anchors = set(), set(), set()
    for v in view['adj']:
        t = tau(view, v, eps)
        if b[v] > t:
            pos.add(v)
        elif b[v] < -t:
            neg.add(v)
        if b[v] > 0:
            anchors.add(v)
    return pos, neg, anchors


def leaf_envelope_test(view, eps):
    """S3: h(r) - d_r <= 2 eps at every observed leaf (necessary for any k)."""
    two = 2 * F(eps)
    return all(s <= two for s in view['envelope_slack'].values())


def pendant_test(view, eps):
    """S4: Y >= -eps (|R| - 1)."""
    return view['Y'] >= -F(eps) * (len(view['leaves']) - 1)


# ------------------------------------------------- mass screen: required atoms

def required_vertex_sets(view, eps, k, budget_seconds=600):
    """Every (Sigma, Required) the mass screen admits for support size k.

    Sigma runs over the frozen `deficit_sets` enumeration (weighted deficit
    bound, positive anchors, eligible degrees), which F-D and F-F justify
    verbatim under a tolerance.  Required = POS(eps) \\ N[Sigma] must be covered
    by the atoms of the support and has at most 2k elements.
    """
    pos, neg, anchors = sign_sets_sharp(view, eps)
    adj, deg = view['adj'], view['deg']
    eligible = {v for v in adj if 3 <= deg[v] <= k}
    out = []
    budget = Budget(budget_seconds, 10 ** 7, 10 ** 9, 10 ** 9)
    for S in deficit_sets(adj, anchors & eligible, eligible, max(0, k - 2), budget):
        A = set(S)
        for v in S:
            A.update(adj[v])
        if neg - A:
            continue
        # POS \ Sigma (not POS \ N[Sigma]) is forced into supp(m): for v not in
        # Sigma every Sigma-neighbour has phi > eps, which only raises m(v).
        required = frozenset(pos - set(S))
        if len(required) > 2 * k:
            continue
        out.append((tuple(S), required))
    seen, uniq = set(), []
    for S, req in out:
        if req not in seen:
            seen.add(req)
            uniq.append((S, req))
    return uniq, dict(n_pos=len(pos), n_neg=len(neg), n_anchor=len(anchors),
                      n_sigma=len(out), n_required_sets=len(uniq))


def supports_covering(view, required, k, n_edges):
    """All k-subsets of original edges whose atoms cover `required`.

    Exact enumeration: pick a set cover of `required` by at most k edges, then
    fill the remaining slots with every other edge.  Yields each support once.
    """
    atoms = view['atoms']
    req = sorted(required)
    if not req:
        yield from itertools.combinations(range(n_edges), k)
        return
    by_vertex = {w: [i for i in range(n_edges) if w in atoms[i]] for w in req}
    if any(not by_vertex[w] for w in req):
        return
    emitted = set()

    def cores(idx, chosen):
        """Minimal-ish covers: assign the first uncovered required vertex."""
        covered = set()
        for i in chosen:
            covered |= atoms[i]
        rest = [w for w in req if w not in covered]
        if not rest:
            yield tuple(sorted(chosen))
            return
        if len(chosen) == k:
            return
        for i in by_vertex[rest[0]]:
            if i in chosen:
                continue
            yield from cores(idx + 1, chosen + [i])

    for core in set(cores(0, [])):
        slots = k - len(core)
        if slots == 0:
            if core not in emitted:
                emitted.add(core)
                yield core
            continue
        pool = [i for i in range(n_edges) if i not in core]
        for extra in itertools.combinations(pool, slots):
            S = tuple(sorted(core + extra))
            if S not in emitted:
                emitted.add(S)
                yield S


def mass_screen_candidates(view, eps, k, n_edges, cap=None):
    """Union over admissible Sigma of the supports the mass screen allows."""
    sets, stats = required_vertex_sets(view, eps, k)
    if not sets:
        return [], dict(stats, admissible=False, n_candidates=0)
    seen = []
    pool = set()
    for _S, req in sets:
        for S in supports_covering(view, req, k, n_edges):
            pool.add(S)
            if cap is not None and len(pool) > cap:
                return None, dict(stats, admissible=True, n_candidates=len(pool),
                                  overflow=True)
    seen = sorted(pool)
    return seen, dict(stats, admissible=True, n_candidates=len(seen), overflow=False)


# ---------------------------------------------------- exact theta on a support

def leaf_distance_cache(edges, leaves):
    """Lazy exact d(v, .) over the observed leaves, built with the frozen BFS."""
    store = {}

    def get(v):
        if v not in store:
            dist = common.bfs_from(edges, [v])[v]
            store[v] = tuple(dist[r] for r in leaves)
        return store[v]
    return get, store


def theta_of_support(getcol, edges, support, dvals, cutoff=None, seed=()):
    """Frozen stage-296 LP on one support; columns from the lazy cache."""
    cols = []
    for e in support:
        u, v, _l = edges[e]
        cols.append(getcol(u))
        cols.append(getcol(v))
    return s296.theta_support(cols, dvals, len(support), cutoff=cutoff, seed=list(seed))


def theta_min_over(candidates, getcol, edges, dvals, k, cutoff=None, hot_rows=8):
    """Minimum of theta over an explicit candidate list, exact, with pruning."""
    best, best_S, best_z = cutoff, None, None
    hot, visited, pruned = [], 0, 0
    for S in candidates:
        visited += 1
        out = theta_of_support(getcol, edges, S, dvals, cutoff=best, seed=hot[:hot_rows])
        if out['status'] == 'pruned':
            pruned += 1
            continue
        if out['status'] != 'optimal':
            return dict(status=out['status'], support=S, visited=visited)
        hot = list(dict.fromkeys(out['active'] + hot))[:hot_rows]
        if best is None or out['value'] < best:
            best, best_S, best_z = out['value'], S, out['z']
    return dict(status='complete', theta=best, support=best_S, z=best_z,
                visited=visited, pruned=pruned)


# ------------------------------------------- envelope screen (S5), quantitative

def _float_dist(adj, root):
    """Float single-source distances; used only to CHOOSE constraint leaves."""
    d = {root: 0.0}
    stack = [root]
    while stack:
        v = stack.pop()
        dv = d[v]
        for u, l in adj[v].items():
            if u not in d:
                d[u] = dv + float(l)
                stack.append(u)
    return d


def _exact_dist(adj, root):
    d = {root: ZERO}
    stack = [root]
    while stack:
        v = stack.pop()
        dv = d[v]
        for u, l in adj[v].items():
            if u not in d:
                d[u] = dv + l
                stack.append(u)
    return d


def envelope_constraint(view, s, eps, exact=True):
    """One S5 constraint anchored at observed leaf s.

    Returns rho, the attaining partner r, and proj[v] = d(s, proj_{rs}(v)),
    the Gromov product (d(r,s) + d(s,v) - d(r,v)) / 2.  A support may place
    mass only in atoms meeting {v : proj[v] <= rho}: otherwise every atom of
    the support carries Psi > rho, contradicting S5.
    """
    adj, obs, leaves = view['adj'], view['obs'], view['leaves']
    dist = (_exact_dist if exact else _float_dist)(adj, s)
    best_r, best_val = None, None
    for r in leaves:
        if r == s:
            continue
        val = obs[r] - dist[r]
        if best_val is None or val > best_val:
            best_val, best_r = val, r
    dr = (_exact_dist if exact else _float_dist)(adj, best_r)
    drs = dist[best_r]
    rho = (drs - (obs[best_r] - obs[s])) / 2 + (F(eps) if exact else float(eps))
    proj = {v: (drs + dist[v] - dr[v]) / 2 for v in adj}
    return dict(s=s, r=best_r, rho=rho, proj=proj, d_rs=drs)


def constraint_mask(view, con, n_edges):
    """Boolean over ORIGINAL edges: may this edge be in a support at all."""
    import numpy as np
    proj, rho, atoms = con['proj'], con['rho'], view['atoms']
    allowed = {v for v in proj if proj[v] <= rho}
    mask = np.zeros(n_edges, dtype=bool)
    for i in range(n_edges):
        if atoms[i] & allowed:
            mask[i] = True
    return mask


def choose_constraints(view, eps, n_edges, sample, n_keep, rng):
    """Sample observed leaves, score each constraint by its kernel size (float),
    and return the n_keep tightest, recomputed EXACTLY."""
    leaves = view['leaves']
    pick = list(leaves) if len(leaves) <= sample else rng.sample(list(leaves), sample)
    scored = []
    for s in pick:
        con = envelope_constraint(view, s, eps, exact=False)
        if con['rho'] < 0:
            return None, dict(infeasible_at=s)
        size = int(constraint_mask(view, con, n_edges).sum())
        scored.append((size, s))
    scored.sort()
    keep = [s for _sz, s in scored[:n_keep]]
    cons = []
    for s in keep:
        con = envelope_constraint(view, s, eps, exact=True)
        if con['rho'] < 0:
            return None, dict(infeasible_at=s)
        cons.append(con)
    return cons, dict(sampled=len(pick), kept=len(cons),
                      float_kernel_sizes=[sz for sz, _s in scored[:n_keep]])


def envelope_survivors_k2(view, cons, n_edges):
    """Every pair of original edges surviving all S5 constraints, exactly."""
    import numpy as np
    masks = np.stack([constraint_mask(view, con, n_edges) for con in cons])
    sizes = masks.sum(axis=1)
    order = np.argsort(sizes)
    masks = masks[order]
    base = np.flatnonzero(masks[0])
    if base.size == 0:
        return None, dict(kernel_sizes=sizes[order].tolist(), empty_kernel=True)
    notm = ~masks
    packed = np.packbits(notm.T.copy(), axis=1)          # (n_edges, ceil(C/8))
    surv = []
    allv = np.arange(n_edges)
    for e in base:
        other = allv[allv != e]
        bad = np.bitwise_and(packed[other], packed[e]).any(axis=1)
        good = other[~bad]
        for f in good:
            surv.append((int(e), int(f)) if e < f else (int(f), int(e)))
    surv = sorted(set(surv))
    return surv, dict(kernel_sizes=sizes[order].tolist(), base_kernel=int(base.size),
                      n_survivors=len(surv))


def envelope_kernel_cover(view, cons, n_edges, k):
    """Greedy pairwise-disjoint kernels: more than k of them refutes size k."""
    import numpy as np
    masks = [constraint_mask(view, con, n_edges) for con in cons]
    sizes = [int(m.sum()) for m in masks]
    order = sorted(range(len(masks)), key=lambda i: sizes[i])
    chosen, union = [], np.zeros(n_edges, dtype=bool)
    for i in order:
        if not (masks[i] & union).any():
            chosen.append(i)
            union |= masks[i]
    return chosen, [sizes[i] for i in chosen]


def candidate_plan(view, eps, k, n_edges, hard_cap=None, disabled=False):
    """Plan the support enumeration under the mass screen.

    Returns (iterator, n_candidates_or_None, stats).  `n_candidates` is None
    when the screen is VACUOUS (some admissible Sigma leaves nothing required),
    in which case the iterator is the full C(|E|, k) sweep and the stage reports
    the screen as vacuous rather than as a success.
    """
    if disabled:
        total = 1
        for i in range(k):
            total = total * (n_edges - i) // (i + 1)
        return itertools.combinations(range(n_edges), k), total, \
            dict(disabled=True, vacuous=True)
    sets, stats = required_vertex_sets(view, eps, k)
    if not sets:
        return iter(()), 0, dict(stats, admissible=False, vacuous=False)
    if any(len(req) == 0 for _S, req in sets):
        total = 1
        for i in range(k):
            total = total * (n_edges - i) // (i + 1)
        return itertools.combinations(range(n_edges), k), total, \
            dict(stats, admissible=True, vacuous=True)
    pool = set()
    for _S, req in sets:
        for S in supports_covering(view, req, k, n_edges):
            pool.add(S)
            if hard_cap is not None and len(pool) > hard_cap:
                return None, None, dict(stats, admissible=True, vacuous=False,
                                        overflow=True)
    out = sorted(pool)
    return iter(out), len(out), dict(stats, admissible=True, vacuous=False,
                                     overflow=False)


def theta_localised(edges, observed, k, eps, hard_cap=2 * 10 ** 6, Dv=None, dvals=None,
                    leaves=None, cutoff=None):
    """Theta_k under the mass screen, every survivor decided by the frozen LP."""
    view = backbone_view(edges, observed)
    if not leaf_envelope_test(view, eps) or not pendant_test(view, eps):
        return dict(status='refuted', reason='leaf_or_pendant_test', theta=None,
                    visited=0, n_candidates=0, screen=dict(vacuous=False))
    it, n_cand, stats = candidate_plan(view, eps, k, len(edges), hard_cap=hard_cap)
    if it is None:
        return dict(status='overflow', theta=None, visited=0,
                    n_candidates=None, screen=stats)
    if n_cand == 0:
        return dict(status='refuted', reason='mass_screen', theta=None,
                    visited=0, n_candidates=0, screen=stats)
    if leaves is None:
        leaves = sorted(observed)
        dvals = tuple(F(observed[r]) for r in leaves)
    if Dv is None:
        Dv = s296.leaf_vectors(edges, leaves)
    best, best_S, best_z, hot, visited, pruned = cutoff, None, None, [], 0, 0
    for S in it:
        visited += 1
        out = s296.theta_support(s296.support_cols(Dv, edges, S), dvals, k,
                                 cutoff=best, seed=hot[:8])
        if out['status'] == 'pruned':
            pruned += 1
            continue
        if out['status'] != 'optimal':
            return dict(status='lp_' + out['status'], support=S, visited=visited,
                        screen=stats, theta=None)
        hot = list(dict.fromkeys(out['active'] + hot))[:8]
        if best is None or out['value'] < best:
            best, best_S, best_z = out['value'], S, out['z']
    decided = best is not None and best <= eps
    return dict(status='complete', theta=best, support=best_S, z=best_z,
                visited=visited, pruned=pruned, n_candidates=n_cand,
                screen=stats, is_exact_Theta_k=bool(decided),
                certifies_greater_than_eps=bool(not decided))


# ------------------------------------- sharpened tolerance bands and Theta_k >=

def bands(view, v):
    """Per-incident-edge caps of the sharpened F-E bound.

    phi(u) - phi(v) = (f(u) - f(v)) - (h(u) - h(v)) with |f(u) - f(v)| <= l and
    |phi| <= eps off Sigma, so the increment lies in
        [max(-2 eps, -l - dh), min(2 eps, l - dh)],  dh = h(u) - h(v),
    giving   m(v) - b(v) <= sum_u min(eps/l, (l - dh)/(2l)) =: tau_up(v)
    and      m(v) - b(v) >= -sum_u min(eps/l, (l + dh)/(2l)) =: -tau_dn(v).
    Both caps lie in [0, 1], so this dominates the plain min(eps/l, 1) band.
    Their eps-free sums are 1 - b(v) and deg(v) - 1 + b(v), which reproduce the
    trivial 0 <= m(v) <= 1 in the limit.
    """
    h, adj = view['h'], view['adj']
    out = []
    for u, l in adj[v].items():
        dh = h[u] - h[v]
        out.append((l, (l - dh) / (2 * l), (l + dh) / (2 * l)))
    return out


def tau_up(view, v, eps):
    return sum(min(F(eps) / l, cu) for l, cu, _cd in bands(view, v))


def tau_dn(view, v, eps):
    return sum(min(F(eps) / l, cd) for l, _cu, cd in bands(view, v))


def _critical(items, target):
    """Smallest eps >= 0 with sum_u min(eps/l_u, c_u) >= target, or None.

    Each term is piecewise linear in eps with a single breakpoint at c_u l_u.
    """
    if target <= 0:
        return ZERO
    if sum(c for _l, c in items) < target:
        return None
    breaks = sorted(set([ZERO] + [c * l for l, c in items]))
    for a, bnd in zip(breaks, breaks[1:] + [None]):
        slope = sum(1 / l for l, c in items if c * l > a)
        const = sum(c for l, c in items if c * l <= a)
        if slope == 0:
            continue
        x = (target - const) / slope
        if x >= a and (bnd is None or x <= bnd):
            return x
    return breaks[-1]


def critical_eps(view):
    """Exact eps at which each vertex leaves POS / leaves NEG.

    v is in POS(eps) iff eps < pos[v]; v is in NEG(eps) iff eps < neg[v].
    """
    pos, neg = {}, {}
    for v in view['adj']:
        bv = view['b'][v]
        bd = bands(view, v)
        if bv > 0:
            x = _critical([(l, cd) for l, _cu, cd in bd], bv)
            if x is not None:
                pos[v] = x
        elif bv < 0:
            x = _critical([(l, cu) for l, cu, _cd in bd], -bv)
            if x is not None:
                neg[v] = x
    return pos, neg


def refutation_threshold(view, k):
    """Largest eps the frozen count obstruction still refutes, exactly.

    S1/S2 give NEG(eps) inside N[Sigma] and POS(eps) inside supp(m) union
    Sigma, so |NEG| <= max(0, 4k - 8) and |POS| <= 3k - 2 (2k for k <= 2) are
    necessary for theta(S) <= eps at some size-k support.  Both sets shrink as
    eps grows, so the refuted region is an interval [0, eps_R) and
    Theta_k >= eps_R.
    """
    pos, neg = critical_eps(view)
    pos_cap = (2 * k) if k <= 2 else (3 * k - 2)
    neg_cap = max(0, 4 * k - 8)
    pv = sorted(pos.values(), reverse=True)
    nv = sorted(neg.values(), reverse=True)
    e_pos = pv[pos_cap] if len(pv) > pos_cap else ZERO
    e_neg = nv[neg_cap] if len(nv) > neg_cap else ZERO
    return max(e_pos, e_neg), dict(n_pos_total=len(pv), n_neg_total=len(nv),
                                   pos_cap=pos_cap, neg_cap=neg_cap,
                                   eps_from_pos=e_pos, eps_from_neg=e_neg)


def sign_sets_sharp(view, eps):
    """POS(eps) / NEG(eps) under the sharpened bands, plus unrelaxed anchors."""
    pos, neg, anchors = set(), set(), set()
    for v in view['adj']:
        bv = view['b'][v]
        if bv > 0:
            anchors.add(v)
            if bv > tau_dn(view, v, eps):
                pos.add(v)
        elif bv < 0 and bv + tau_up(view, v, eps) < 0:
            neg.add(v)
    return pos, neg, anchors
