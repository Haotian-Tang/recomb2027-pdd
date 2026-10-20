#!/usr/bin/env python3
r"""Stage295 prototype: the frozen Stage43 minimum pass with a provably
equivalent, strictly smaller zero-face enumeration.

NOTHING FROZEN IS MODIFIED.  Every primitive is imported from the frozen
modules -- `observed_backbone`, `observe`, `tree`, `two_pass_envelope`,
`inverse`, `deficit_sets`, `mixed_cover`, `realize`, `lift_mixed`, `Budget`,
`ResourceLimit`, `dot` and Stage43's `strict_vertex`.  The control flow of
`stage43_partial_fpt.solve_mixed` is reproduced line for line except for the
sign-mask loop.

WHY THE SMALLER ENUMERATION IS EQUIVALENT
-----------------------------------------
With the deficit support S fixed, the frozen search enumerates all 2^|A| sign
patterns over A = S u N(S), where the mask selects the set P of coordinates
constrained strictly positive and the complement is pinned to m_v = 0 exactly.

On every face the frozen constraint block already forces t_v >= s > 0 for all
v in S (rows ``[-e_j, 1] . (t, s) <= 0``).  For u in N(S) \ S the mass row
carries only the positive incidence terms, so

    m_u = b_u + sum_{v in S, v ~ u} t_v / (2 l_uv)  >  b_u .

Hence if b_u >= 0 then m_u > 0 strictly on the whole strict face, and any mask
that places u outside P -- i.e. imposes m_u = 0 -- is infeasible.  Those masks
cannot produce an incumbent, so skipping them cannot change the reported
minimum.  Only S itself and the strictly negative baselines of N(S) \ S are
genuinely free:

    free(S) = S u {u in N(S) \ S : b_u < 0},     2^|free| <= 2^|A| .

This is the same reduction Stage292 states and cross-checks against the full
2^|A| enumeration; here it is applied to the MINIMUM pass rather than to the
all-minimisers pass.  It does not change the worst-case exponent (if every
neighbour has a negative baseline then free(S) = A); it removes provably
infeasible work, and the size of that saving is an empirical question measured
in Stage295.

Two further switches, both exact, are measured separately:

  memo_cover        `mixed_cover` depends only on (P, marked, lookup, edges,
                    cap); within one case the last four are fixed, so its
                    value is a pure function of frozenset(P) and is cached.
  popcount_order    visiting masks in order of increasing |P| tends to produce
                    a small incumbent earlier, which strengthens the frozen
                    ``cover['count'] >= best['source_count']`` cut.  Visiting
                    order cannot change a complete search's minimum.

`full_mask=True` restores the frozen 2^|A| enumeration so the two can be
compared inside one harness.
"""
from fractions import Fraction as F
from itertools import combinations

from stage23_two_source_envelope import two_pass_envelope, inverse
from stage25_deficit_extension import tree
from stage26_missing_leaf import observed_backbone, observe
from stage27_general_fpt import Budget, ResourceLimit, dot
from stage41_partial_fpt import mixed_cover, deficit_sets, realize, lift_mixed
from stage43_partial_fpt import strict_vertex


def _masks(n, popcount_order):
    if not popcount_order:
        return range(1 << n)
    def gen():
        for size in range(n + 1):
            for chosen in combinations(range(n), size):
                m = 0
                for j in chosen:
                    m |= 1 << j
                yield m
    return gen()


def solve_mixed_reduced(edges, observed, marked=(), cap=3, seconds=30,
                        max_states=200000, max_bases=2000000, max_faces=500000,
                        full_mask=False, memo_cover=True, popcount_order=False,
                        max_cache=200000):
    """All leaves of a reduced tree observed; marked internal vertex atoms legal."""
    if not isinstance(cap, int) or isinstance(cap, bool) or cap < 1:
        raise ValueError('Positive integer source budget required')
    budget = Budget(seconds, max_states, max_bases, max_faces)
    stats = {'deficit_sets': 0, 'cover_feasible_faces': 0, 'strict_faces': 0,
             'skipped_forced_faces': 0, 'cover_memo_hits': 0, 'cover_calls': 0,
             'cover_cache_full': 0}
    best = None
    cache = {}

    def cover_of(P):
        if not memo_cover:
            stats['cover_calls'] += 1
            return mixed_cover(P, marked, lookup, edges, cap)
        key = frozenset(P)
        if key in cache:
            stats['cover_memo_hits'] += 1
            return cache[key]
        stats['cover_calls'] += 1
        value = mixed_cover(P, marked, lookup, edges, cap)
        # A bounded cache: beyond max_cache entries nothing further is
        # stored, so memory stays flat.  The cache is a pure function
        # table, so capping it cannot change any answer.
        if len(cache) < max_cache:
            cache[key] = value
        else:
            stats['cover_cache_full'] = 1
        return value

    def report(reason=None):
        answer = {'status': 'complete', 'feasible': best is not None,
                  'statistics': dict(stats, **budget.counts)}
        if reason:
            answer['reason'] = reason
        if best is not None:
            answer.update(best, minimum_sources_within_cap=best['source_count'])
        return answer
    try:
        adj, lookup = tree(edges)
        leaves = {v for v in adj if len(adj[v]) == 1}
        marked = set(marked)
        if set(observed) != leaves or not marked <= set(adj) - leaves:
            raise ValueError('All backbone leaves and only internal vertex marks required')
        obs = {v: F(d) for v, d in observed.items()}
        h, _ = two_pass_envelope(adj, obs)
        if any(h[r] != obs[r] for r in leaves):
            return report('inconsistent_leaf_envelope')
        b, Y = inverse(adj, h)
        if Y <= 0:
            return report('nonpositive_baseline_pendant')
        positive = {v for v in adj if b[v] > 0}
        negative = {v for v in adj if b[v] < 0}
        if len(positive) > (2 * cap if cap <= 2 else 3 * cap - 2) \
                or len(negative) > max(0, 4 * cap - 8):
            return report('baseline_mass_count_obstruction')
        eligible = {v for v in adj if 3 <= len(adj[v]) <= cap}
        for S in deficit_sets(adj, positive & eligible, eligible, max(0, cap - 2), budget):
            stats['deficit_sets'] += 1
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
                free = [v for v in ordered if v in Sset or b[v] < 0]
            forced = [v for v in ordered if v not in set(free)]
            stats['skipped_forced_faces'] += (1 << len(ordered)) - (1 << len(free))
            base = fixed | set(forced)
            for mask in _masks(len(free), popcount_order):
                budget.touch('mass_faces')
                P = base | {v for j, v in enumerate(free) if mask & (1 << j)}
                cover = cover_of(P)
                if cover is None or best is not None and cover['count'] >= best['source_count']:
                    continue
                stats['cover_feasible_faces'] += 1
                d = len(S) + 1
                lower_t = (F(0),) * len(S) + (F(-1),)
                upper_t = (F(0),) * len(S) + (F(1),)
                constraints = [(lower_t, F(0)), (upper_t, F(1))]
                for j in range(len(S)):
                    row = [-F(i == j) for i in range(len(S))] + [F(1)]
                    constraints.append((tuple(row), F(0)))
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
                stats['strict_faces'] += 1
                delta = x[:-1]
                mass = {v: b[v] + (dot(rows[v], delta) if v in A else 0) for v in P}
                y = Y - sum(F(len(adj[v]) - 2, 2) * z for v, z in zip(S, delta))
                best = realize(edges, mass, cover, y, S, delta)
                if best['source_count'] == 1:
                    return report()
        return report()
    except ResourceLimit as exc:
        return {'status': 'unknown', 'reason': str(exc), 'incumbent': best,
                'feasibility_established': best is not None,
                'statistics': dict(stats, **budget.counts)}


def solve_partial_reduced(edges, observed, cap=3, **options):
    """At least two observed original leaves, free positive mean pendant."""
    if len(observed) < 2:
        raise ValueError('Prototype requires at least two observed leaves')
    observed = {v: F(d) for v, d in observed.items()}
    backbone = observed_backbone(edges, observed)
    marked = {x['vertex'] for x in backbone['saturation_violations']}
    answer = solve_mixed_reduced(backbone['edges'], observed, marked, cap, **options)
    answer['marked_vertices'] = sorted(marked)
    answer['backbone_edges'] = backbone['edges']
    mixed = answer if answer.get('feasible') else answer.get('incumbent')
    if mixed is not None:
        lifted = lift_mixed(backbone, mixed['atoms'], mixed['pendant'])
        assert observe(edges, lifted['sources'], lifted['pendant'], observed) == observed
        answer['original_witness'] = lifted
    return answer
