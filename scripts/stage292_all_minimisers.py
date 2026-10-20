#!/usr/bin/env python3
r"""Stage292 -- complete exact enumeration of ALL minimum-size strict PDD source
supports, replacing the C(|E|, K) brute force used by Stage290.

MODEL.  Reduced tree T with positive rational edge lengths; observed reference
set R (a subset of the leaves of T) with exact rational observations d_r; k
sources strictly interior to k distinct edges with weights p_i > 0 summing to
one; a common mean pendant y > 0; and

    d_r = y + sum_i p_i * dist_T(r, q_i),  r in R.

MASS FORM.  For q at fraction x of edge (u, v) and ANY vertex r,
dist(r, q) = (1 - x) d(r, u) + x d(r, v).  A configuration is therefore
equivalent to a vertex mass vector m >= 0 with sum_v m_v = 1 and
d_r = y + sum_v m_v d(r, v).  A set F of k edges is a support of the
configuration iff m is strictly positive exactly on the endpoint set
U = V(F); conversely any m that is strictly positive exactly on V(F) splits
over F with every x_i strictly interior (give each v in U a positive share to
each incident edge of F).  Hence

    F is a minimum-size support  <=>  V(F) is feasible and |F| is minimum,

and feasibility of an endpoint set is a single exact rational linear
feasibility problem -- the frozen Stage25 oracle `fixed_support_polytope`
(all endpoint masses and the pendant strictly positive).  The only
combinatorial question is WHICH endpoint sets to test.  Brute force tests all
C(|E|, k) of them.  This module tests O(K)-local candidates only.

WHAT IS ENUMERATED, AND WHY IT IS EXHAUSTIVE
--------------------------------------------
(1) BACKBONE PROJECTION (frozen `observed_backbone`).  Unobserved leaves are
    pruned.  For a vertex w outside the observed hull with gate g (its nearest
    hull vertex) and off-distance a(w), d(r, w) = d(r, g) + a(w) for every
    r in R, so mass at w acts as mass at g plus a constant a(w)*m_w added to
    the pendant.  Hull vertices of hull-degree two are interior points of
    backbone chains, and mass there splits exactly onto the chain endpoints.
    So (m, y) on T is feasible iff its push-forward (m^, y^) on the backbone B
    is feasible, with y^ = y + sum_w m_w a(w) >= y > 0.  B is reduced and its
    leaf set is exactly R, so the all-leaves theory applies on B.

(2) DEFICIT SUPPORTS.  On B let h be the leaf envelope h_v = max_r (d_r -
    d(v, r)); every feasible potential f satisfies f >= h and f = h + t with
    t >= 0 and t = 0 at leaves.  t_v > 0 iff EVERY branch at v carries mass,
    i.e. iff every branch at v contains a source point.  Each source point
    lies in exactly one branch of v, so deg(v) <= k; and taking Y = the
    minimal subtree of B spanning the k source points, deg_Y(v) = deg_B(v)
    for every v with t_v > 0 and #leaves(Y) <= k, whence

        sum_{v: t_v > 0} (deg_B(v) - 2) <= #leaves(Y) - 2 <= k - 2.

    This is the project's weighted deficit bound.  The admissible supports S
    are therefore enumerated by the FROZEN generator `deficit_sets`
    (imported, not re-implemented), which grows connected sets of eligible
    vertices (3 <= deg <= cap) from positive-baseline anchors under that
    weight bound.  Exhaustiveness of that generator is inherited from the
    frozen solver and is cross-checked here on every call by comparing the
    resulting minimum with the frozen `solve_partial`.

(3) ZERO FACES.  With S fixed the mass is affine in t: m_v(t) = b_v +
    (A t)_v, non-zero only on A = S u N(S).  A positive pattern is a choice of
    which coordinates of A vanish.  Vertices outside A keep their baseline
    sign, so negative baselines outside A are infeasible.  Within A, a vertex
    u in N(S) \ S has (A t)_u = sum_{v in S, v~u} t_v / (2 l_uv) > 0 because
    all t_v > 0 strictly on the strict face, so if b_u >= 0 then m_u > 0 is
    FORCED.  Only S itself and the strictly negative baselines of N(S) are
    free, which reduces the face enumeration from 2^|A| to 2^|free|; the full
    2^|A| enumeration is available via full_mask=True and is checked against
    the reduced one in the validation harness.

(4) STRICT FEASIBILITY of a face is decided exactly by the frozen Stage43
    common-slack vertex routine `strict_vertex` on the same constraint system
    the frozen solver builds.

(5) MINIMUM COVERS.  A feasible pattern P is realised by atoms: an EDGE atom
    on a backbone edge with both endpoints in P, or a VERTEX atom at a marked
    backbone vertex (a gate with unobserved material behind it).  The minimum
    number of atoms covering P is |P| - nu(G[P]) with nu the maximum matching
    of the forest induced on P; this module enumerates ALL covers of that
    size, not one.

(6) LIFTS.  Every atom is lifted back to the original tree.  An edge atom on
    backbone edge j sits at an unknown coordinate in (0, L_j), so its possible
    original edges are all original edges of chain j, plus all omitted edges
    behind any chain-interior hull vertex of j (the coordinate may land
    exactly on such a vertex).  A vertex atom at gate g lifts to any omitted
    edge behind g.  These candidate lists are pairwise disjoint across atoms,
    so the lifts of a configuration are their Cartesian product.  The product
    is a SUPERSET of the lifts; every element is then decided exactly by the
    Stage25 oracle, so no false positive survives and no true support is lost.

HONEST SCOPE.  Exhaustiveness rests on (a) the weighted deficit bound and the
positive-anchor growth of the frozen `deficit_sets`, which this module imports
rather than re-derives, and (b) the exactness of the frozen Stage25 endpoint
oracle.  Resource stops are reported as `status != 'complete'` and never as a
uniqueness verdict.  No pre-existing project file is modified.
"""
import os
for _k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_k] = '1'
import itertools
import sys
from collections import deque
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
sys.dont_write_bytecode = True

from stage23_two_source_envelope import two_pass_envelope, inverse       # noqa: E402
from stage25_deficit_extension import tree                               # noqa: E402
from stage27_general_fpt import Budget, ResourceLimit, dot               # noqa: E402
from stage26_missing_leaf import observed_backbone                       # noqa: E402
from stage41_partial_fpt import deficit_sets, mixed_cover                # noqa: E402
from stage43_partial_fpt import strict_vertex, solve_partial             # noqa: E402
from stage25_benchmark import fixed_support_polytope                     # noqa: E402


# ---------------------------------------------------------------------------
# distances from a handful of vertices (never the full all-pairs matrix)
# ---------------------------------------------------------------------------

class VertexDistances(object):
    """d(v, .) for a few vertices v of a possibly very large tree, cached."""

    def __init__(self, edges, max_columns=24):
        self.adj = {}
        for u, v, L in edges:
            L = F(L)
            assert L > 0 and u != v
            self.adj.setdefault(u, []).append((v, L))
            self.adj.setdefault(v, []).append((u, L))
        self.cache = {}
        self.order = deque()
        self.max_columns = max_columns
        self.bfs_calls = 0

    def column(self, v):
        col = self.cache.get(v)
        if col is None:
            self.bfs_calls += 1
            # bounded memoisation: on a 10^5-vertex tree one exact rational
            # distance column costs tens of MB, so old columns are evicted.
            while len(self.order) >= self.max_columns:
                self.cache.pop(self.order.popleft(), None)
            self.order.append(v)
            col = {v: F(0)}
            order = [v]
            for a in order:
                da = col[a]
                for b, L in self.adj[a]:
                    if b not in col:
                        col[b] = da + L
                        order.append(b)
            self.cache[v] = col
        return col

    def rows(self, endpoints, Rlist):
        """{r: {v: d(r, v)}} for r in Rlist, v in endpoints."""
        cols = [self.column(v) for v in endpoints]
        return {r: {v: col[r] for v, col in zip(endpoints, cols)} for r in Rlist}


# ---------------------------------------------------------------------------
# covers of a positive pattern by edge atoms and marked-vertex atoms
# ---------------------------------------------------------------------------

def _incidence(P, edgesB, lookupB):
    inc = {v: [] for v in P}
    for u, v in itertools.combinations(sorted(P), 2):
        j = lookupB.get((u, v))
        if j is not None:
            inc[u].append((j, v))
            inc[v].append((j, u))
    return inc


def _forest_matching(P, inc):
    """Maximum matching of the forest induced on P, by leaf removal."""
    active = set(P)
    degree = {v: len(inc[v]) for v in P}
    queue = deque(sorted(v for v in P if degree[v] == 1))
    nu = 0
    while queue:
        u = queue.popleft()
        if u not in active or degree[u] != 1:
            continue
        w = next((o for _, o in inc[u] if o in active), None)
        if w is None:
            continue
        nu += 1
        active.discard(u)
        active.discard(w)
        for a in (u, w):
            for _, o in inc[a]:
                if o in active:
                    degree[o] -= 1
                    if degree[o] == 1:
                        queue.append(o)
    return nu


def cover_count(P, inc, marked):
    """Minimum number of atoms covering P, or None if P is unrealisable."""
    for v in P:
        if not inc[v] and v not in marked:
            return None
    return len(P) - _forest_matching(P, inc)


def all_min_covers(P, inc, marked, target, limit=200000):
    """Every cover of P by exactly `target` atoms ('e', j) / ('v', v)."""
    Pl = sorted(P)
    found = set()
    state = {'visits': 0}

    def rec(start, covered, atoms):
        state['visits'] += 1
        if state['visits'] > limit:
            raise ResourceLimit('cover_enumeration_limit')
        j = start
        while j < len(Pl) and Pl[j] in covered:
            j += 1
        if j == len(Pl):
            if len(atoms) == target:
                found.add(frozenset(atoms))
            return
        remaining = sum(1 for v in Pl[j:] if v not in covered)
        if len(atoms) + (remaining + 1) // 2 > target:
            return
        v = Pl[j]
        for e, other in inc[v]:
            rec(j + 1, covered | {v, other}, atoms + (('e', e),))
        if v in marked:
            rec(j + 1, covered | {v}, atoms + (('v', v),))

    rec(0, frozenset(), ())
    return found



# ---------------------------------------------------------------------------
# feasible positive patterns of the backbone mass vector
# ---------------------------------------------------------------------------

def backbone_patterns(edgesB, observed, marked, cap, budget, accept,
                      full_mask=False, stats=None):
    """All strictly feasible positive patterns P whose cover count `accept`s.

    Returns (patterns, reason).  `reason` is None when the enumeration ran;
    otherwise it names the baseline obstruction that makes every pattern
    infeasible within `cap`.
    """
    stats = stats if stats is not None else {}
    adj, lookupB = tree(edgesB)
    leaves = {v for v in adj if len(adj[v]) == 1}
    marked = set(marked)
    obs = {v: F(d) for v, d in observed.items()}
    if set(obs) != leaves:
        raise ValueError('backbone leaves must be exactly the observed set')
    h, _ = two_pass_envelope(adj, obs)
    if any(h[r] != obs[r] for r in leaves):
        return [], 'inconsistent_leaf_envelope'
    b, Y = inverse(adj, h)
    if Y <= 0:
        return [], 'nonpositive_baseline_pendant'
    positive = {v for v in adj if b[v] > 0}
    negative = {v for v in adj if b[v] < 0}
    if len(positive) > (2 * cap if cap <= 2 else 3 * cap - 2) \
            or len(negative) > max(0, 4 * cap - 8):
        return [], 'baseline_mass_count_obstruction'
    eligible = {v for v in adj if 3 <= len(adj[v]) <= cap}
    out = []
    for S in deficit_sets(adj, positive & eligible, eligible, max(0, cap - 2), budget):
        stats['deficit_sets'] = stats.get('deficit_sets', 0) + 1
        A = set(S)
        for v in S:
            A.update(adj[v])
        if negative - A:
            continue
        fixed = positive - A
        rows = {v: [F(0)] * len(S) for v in A}
        for j, v in enumerate(S):
            rows[v][j] = -sum(1 / length for length in adj[v].values()) / 2
            for u, length in adj[v].items():
                rows[u][j] += 1 / (2 * length)
        ordered = sorted(A)
        Sset = set(S)
        if full_mask:
            free = list(ordered)
        else:
            # m_u > 0 is forced for u in N(S) \ S with b_u >= 0 (see docstring)
            free = [v for v in ordered if v in Sset or b[v] < 0]
        forced = [v for v in ordered if v not in set(free)]
        base = fixed | set(forced)
        for mask in range(1 << len(free)):
            budget.touch('mass_faces')
            P = base | {v for j, v in enumerate(free) if mask & (1 << j)}
            if not P or len(P) > 2 * cap:
                continue
            inc = _incidence(P, edgesB, lookupB)
            count = cover_count(P, inc, marked)
            if count is None or count > cap or not accept(count):
                continue
            stats['faces_to_lp'] = stats.get('faces_to_lp', 0) + 1
            d = len(S) + 1
            constraints = [((F(0),) * len(S) + (F(-1),), F(0)),
                           ((F(0),) * len(S) + (F(1),), F(1))]
            for j in range(len(S)):
                constraints.append((tuple([-F(i == j) for i in range(len(S))] + [F(1)]), F(0)))
            constraints.append((tuple(F(len(adj[v]) - 2, 2) for v in S) + (F(1),), Y))
            for v in ordered:
                if v in P:
                    constraints.append((tuple(-a for a in rows[v]) + (F(1),), b[v]))
                else:
                    constraints.extend([(tuple(rows[v]) + (F(0),), -b[v]),
                                        (tuple(-a for a in rows[v]) + (F(0),), b[v])])
            x = strict_vertex(constraints, d, budget)
            if x is None:
                continue
            delta = x[:-1]
            mass = {v: b[v] + (dot(rows[v], delta) if v in A else F(0)) for v in P}
            y = Y - sum((F(len(adj[v]) - 2, 2) * z for v, z in zip(S, delta)), F(0))
            assert all(m > 0 for m in mass.values()) and y > 0
            assert sum(mass.values()) == 1
            stats['feasible_patterns'] = stats.get('feasible_patterns', 0) + 1
            out.append(dict(S=tuple(S), P=frozenset(P), inc=inc, count=count,
                            mass=mass, pendant=y))
    return out, None


# ---------------------------------------------------------------------------
# lifting backbone atoms to original edges
# ---------------------------------------------------------------------------

def lift_tables(B):
    """original-edge candidate lists for every backbone atom."""
    behind = {}
    for i in B['omitted_edges']:
        u = B['original_edges'][i][0]
        behind.setdefault(B['gate'][u], []).append(i)
    edge_atom = []
    for j, chain in enumerate(B['chains']):
        cand = [part['edge'] for part in chain]
        for part in chain[:-1]:
            cand.extend(behind.get(part['v'], ()))
        edge_atom.append(sorted(set(cand)))
    vertex_atom = {g: sorted(behind.get(g, ())) for g in B['gate_edges']}
    return edge_atom, vertex_atom


# ---------------------------------------------------------------------------
# main entry point
# ---------------------------------------------------------------------------

def all_minimum_supports(edges, observed, cap, seconds=None, max_states=None,
                         max_bases=None, max_faces=None, max_oracle=400000,
                         stop_after=None, full_mask=False, kstar=None,
                         screen_rows=0, distances=None, check_frozen=True):
    """Complete set of minimum-size strict source supports on the ORIGINAL tree.

    stop_after   stop once this many distinct supports are confirmed (the
                 verdict is then a lower bound; status becomes 'early_exit')
    kstar        skip the frozen minimum pass and use this value
    screen_rows  if > 0, pre-screen each candidate on that many observation
                 rows (an exact necessary condition) before the full oracle
    """
    E = [(u, v, F(L)) for u, v, L in edges]
    observed = {int(r): F(d) for r, d in observed.items()}
    Rlist = sorted(observed)
    stats = dict(deficit_sets=0, faces_to_lp=0, feasible_patterns=0,
                 configurations=0, lift_products=0, oracle_calls=0,
                 oracle_cached=0, screen_rejects=0, bfs_calls=0,
                 supports_visited=0, endpoint_sets_solved=0,
                 endpoint_cache_hits=0, equality_rows_processed=0,
                 vertex_bases_attempted=0)
    budget = Budget(seconds, max_states, max_bases, max_faces)
    out = dict(status='complete', minimum=None, supports=[], n_supports=0,
               statistics=stats, reason=None)
    try:
        B = observed_backbone(E, Rlist)
        marked = {x['vertex'] for x in B['saturation_violations']}
        if kstar is None:
            frozen = solve_partial(E, observed, cap=cap, seconds=seconds,
                                   max_states=max_states, max_bases=max_bases,
                                   max_faces=max_faces)
            if frozen.get('status') != 'complete':
                out.update(status='frozen_' + str(frozen.get('status')),
                           reason=frozen.get('reason'))
                return out
            kstar = frozen.get('minimum_sources_within_cap')
            out['frozen_minimum'] = kstar
            out['frozen_reason'] = frozen.get('reason')
        if kstar is None:
            out['reason'] = out.get('frozen_reason') or 'infeasible_within_cap'
            return out
        pats, reason = backbone_patterns(B['edges'], observed, marked, cap, budget,
                                         accept=lambda c: c == kstar,
                                         full_mask=full_mask, stats=stats)
        if reason is not None or not pats:
            out.update(minimum=None, reason=reason or 'no_pattern_at_frozen_minimum',
                       status='disagreement' if reason is None else 'complete')
            return out
        edge_atom, vertex_atom = lift_tables(B)
        configs = set()
        for pat in pats:
            configs |= all_min_covers(pat['P'], pat['inc'], marked, kstar)
        stats['configurations'] = len(configs)
        if distances is None:
            distances = VertexDistances(E)
        cache = {}
        supports = set()
        screen = Rlist[:screen_rows] if screen_rows else None
        if screen is not None and len(screen) >= len(Rlist):
            screen = None
        for cfg in sorted(configs, key=lambda z: sorted(z)):
            lists = []
            for kind, obj in sorted(cfg):
                lists.append(edge_atom[obj] if kind == 'e' else vertex_atom[obj])
            size = 1
            for L in lists:
                size *= len(L)
            stats['lift_products'] += size
            if stats['oracle_calls'] + size > max_oracle:
                out['status'] = 'oracle_limit'
                break
            for combo in itertools.product(*lists):
                budget.check()
                sup = tuple(sorted(combo))
                if len(set(sup)) != len(sup) or sup in supports:
                    continue
                stats['supports_visited'] += 1
                ep = tuple(sorted({v for i in sup for v in E[i][:2]}))
                good = cache.get(ep)
                if good is None:
                    stats['oracle_calls'] += 1
                    rows = distances.rows(ep, Rlist)
                    if screen is not None:
                        sub = {r: rows[r] for r in screen}
                        if fixed_support_polytope(ep, screen, sub, observed, stats) is None:
                            stats['screen_rejects'] += 1
                            cache[ep] = False
                            good = False
                    if good is None:
                        good = fixed_support_polytope(ep, Rlist, rows, observed,
                                                      stats) is not None
                        cache[ep] = good
                else:
                    stats['oracle_cached'] += 1
                if good:
                    supports.add(sup)
                    if stop_after is not None and len(supports) >= stop_after:
                        out['status'] = 'early_exit'
                        break
            if out['status'] != 'complete':
                break
        stats['bfs_calls'] = distances.bfs_calls
        out.update(minimum=kstar, supports=sorted(list(s) for s in supports),
                   n_supports=len(supports))
        if out['status'] == 'complete' and not supports:
            out['status'] = 'disagreement'
            out['reason'] = 'no_lift_confirmed_at_frozen_minimum'
    except ResourceLimit as exc:
        out.update(status='resource_limit', reason=str(exc))
    return out


def unique_true(result, true_support, K):
    """Stage290's `ok` predicate, recomputed from an all-minimisers result."""
    return bool(result['status'] == 'complete' and result['minimum'] == K and
                result['n_supports'] == 1 and
                list(result['supports'][0]) == list(true_support))

