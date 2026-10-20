#!/usr/bin/env python3
r"""Stage297 -- ALL minimum-size strict PDD source supports with the Stage295
forced-positivity reduction applied to BOTH passes.

WHAT THIS ADDS TO STAGE292
--------------------------
`stage292_all_minimisers.all_minimum_supports` decides uniqueness of the
minimum-source explanation in two passes:

  pass 1  the FROZEN `stage43_partial_fpt.solve_partial` returns the minimum
          source count k* within the cap;
  pass 2  `backbone_patterns(accept = (count == k*))` enumerates every
          strictly feasible zero face whose cover count is k*, and the covers
          and lifts of those faces give every minimum-size support.

Pass 2 already restricts the sign-mask loop to the free coordinates
free(S) = S u {u in N(S)\S : b_u < 0} (Stage292 docstring item 3, cross-checked
there against the full 2^|A| enumeration on 500 fresh instances).  Pass 1 does
not: the frozen solver enumerates all 2^|A| masks, and that is where the K >= 10
ceiling sat (results/stage292_k10_ceiling.txt: `solve_partial` unknown after
3600 s on the smallest K = 10 bank case).  Stage295 removed exactly that
obstruction for the minimum question
(`stage295_reduced_minimum.solve_partial_reduced`).

This module supplies the combination, in three interchangeable modes:

  mode='frozen'            pass 1 = frozen `solve_partial`.  Reproduces
                           Stage292's routine exactly; the reference arm.
  mode='reduced_two_pass'  pass 1 = Stage295 `solve_partial_reduced`, pass 2
                           unchanged.  Same two-pass structure, reduced minimum.
  mode='combined'          ONE pass.  The deficit sets and the reduced face
                           masks are enumerated once; a face is LP-tested when
                           its cover count is <= the best count established so
                           far, the incumbent count is lowered whenever a
                           strictly feasible face beats it, and every feasible
                           face achieving the final incumbent is retained.
                           Minimum and complete minimiser set come out of the
                           same enumeration, so the dominant cost
                           (deficit sets x faces x exact rational LP) is paid
                           once instead of twice.

WHY 'combined' RETURNS THE SAME ANSWER
--------------------------------------
Both passes range over the same objects: the deficit supports S produced by the
frozen `deficit_sets` generator and, for each, the positive patterns P reachable
from the free coordinates.  Let c(P) be the cover count of a strictly feasible
pattern.  The minimum source count within the cap is min_P c(P) over strictly
feasible P (Stage292 docstring items 2-5), and the minimum-size supports are the
lifts of the covers of the P attaining it.  The single pass computes both:

  * a face with c(P) > best cannot attain the minimum, so skipping its LP
    changes neither the minimum nor the retained set;
  * a face with c(P) < best is LP-tested; if strictly feasible it becomes the
    new incumbent and the retained set is reset to it;
  * a face with c(P) == best is LP-tested and retained if strictly feasible.

Since every face with c(P) <= best is LP-tested and `best` only decreases, at
the end `best` = min over feasible faces and the retained set is exactly the
feasible faces attaining it, independent of visiting order.  The forced-
positivity argument for free(S) is Stage295's, verbatim: for u in N(S)\S,
m_u = b_u + sum_{v in S, v~u} t_v/(2 l_uv) > b_u, and t_v > 0 strictly on the
strict face, so b_u >= 0 forces m_u > 0 and every mask that pins m_u = 0 is
infeasible.  `full_mask=True` restores the full 2^|A| enumeration so the two can
be compared inside one harness.

NOTHING FROZEN AND NOTHING PRE-EXISTING IS MODIFIED.  `_incidence`,
`cover_count`, `all_min_covers`, `lift_tables` and `VertexDistances` are
imported from Stage292; `deficit_sets`, `strict_vertex`, `observed_backbone`,
`two_pass_envelope`, `inverse`, `tree`, `fixed_support_polytope`, `Budget`,
`ResourceLimit` and `dot` from the frozen modules; `solve_partial_reduced` from
Stage295.  Resource stops are reported as `status != 'complete'` and are never
read as a uniqueness verdict.
"""
import os
for _k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_k] = '1'
import itertools
import sys
import time
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
sys.dont_write_bytecode = True

from stage23_two_source_envelope import two_pass_envelope, inverse        # noqa: E402
from stage25_deficit_extension import tree                                # noqa: E402
from stage27_general_fpt import Budget, ResourceLimit, dot                # noqa: E402
from stage26_missing_leaf import observed_backbone                        # noqa: E402
from stage41_partial_fpt import deficit_sets                              # noqa: E402
from stage43_partial_fpt import strict_vertex, solve_partial              # noqa: E402
from stage25_benchmark import fixed_support_polytope                      # noqa: E402
from stage292_all_minimisers import (VertexDistances, _incidence,         # noqa: E402
                                     cover_count, all_min_covers, lift_tables)
from stage295_reduced_minimum import solve_partial_reduced                # noqa: E402

MODES = ('combined', 'reduced_two_pass', 'frozen')


def _new_stats():
    return dict(deficit_sets=0, faces_to_lp=0, feasible_patterns=0,
                configurations=0, lift_products=0, oracle_calls=0,
                oracle_cached=0, screen_rejects=0, bfs_calls=0,
                supports_visited=0, endpoint_sets_solved=0,
                endpoint_cache_hits=0, equality_rows_processed=0,
                vertex_bases_attempted=0, masks_enumerated=0,
                masks_full_enumeration=0, masks_skipped_forced=0,
                faces_visited=0, incumbent_drops=0, patterns_discarded=0,
                free_bits_max=0, A_bits_max=0)


def minimum_patterns(edgesB, observed, marked, cap, budget, stats,
                     full_mask=False, kstar=None):
    """Minimum cover count over strictly feasible faces AND every face at it.

    kstar is not None  -> only faces whose cover count equals kstar are
                          LP-tested (Stage292's pass 2, unchanged).
    kstar is None      -> single-pass mode: the incumbent starts empty and is
                          lowered by any strictly feasible smaller face.

    Returns (patterns, minimum, reason).
    """
    adj, lookupB = tree(edgesB)
    leaves = {v for v in adj if len(adj[v]) == 1}
    marked = set(marked)
    obs = {v: F(d) for v, d in observed.items()}
    if set(obs) != leaves:
        raise ValueError('backbone leaves must be exactly the observed set')
    h, _ = two_pass_envelope(adj, obs)
    if any(h[r] != obs[r] for r in leaves):
        return [], None, 'inconsistent_leaf_envelope'
    b, Y = inverse(adj, h)
    if Y <= 0:
        return [], None, 'nonpositive_baseline_pendant'
    positive = {v for v in adj if b[v] > 0}
    negative = {v for v in adj if b[v] < 0}
    if len(positive) > (2 * cap if cap <= 2 else 3 * cap - 2) \
            or len(negative) > max(0, 4 * cap - 8):
        return [], None, 'baseline_mass_count_obstruction'
    eligible = {v for v in adj if 3 <= len(adj[v]) <= cap}
    best = kstar
    out = []
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
            # m_u > 0 forced for u in N(S)\S with b_u >= 0 (see module docstring)
            free = [v for v in ordered if v in Sset or b[v] < 0]
        forced = [v for v in ordered if v not in set(free)]
        stats['masks_full_enumeration'] += 1 << len(ordered)
        stats['masks_enumerated'] += 1 << len(free)
        stats['masks_skipped_forced'] += (1 << len(ordered)) - (1 << len(free))
        stats['free_bits_max'] = max(stats['free_bits_max'], len(free))
        stats['A_bits_max'] = max(stats['A_bits_max'], len(ordered))
        base = fixed | set(forced)
        for mask in range(1 << len(free)):
            budget.touch('mass_faces')
            stats['faces_visited'] += 1
            P = base | {v for j, v in enumerate(free) if mask & (1 << j)}
            if not P or len(P) > 2 * cap:
                continue
            inc = _incidence(P, edgesB, lookupB)
            count = cover_count(P, inc, marked)
            if count is None or count > cap:
                continue
            if best is not None and count > best:
                continue
            if kstar is not None and count != kstar:
                continue
            stats['faces_to_lp'] += 1
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
            stats['feasible_patterns'] += 1
            record = dict(S=tuple(S), P=frozenset(P), inc=inc, count=count,
                          mass=mass, pendant=y)
            if best is None or count < best:
                if out:
                    stats['incumbent_drops'] += 1
                    stats['patterns_discarded'] += len(out)
                best = count
                out = [record]
            else:
                out.append(record)
    return out, best, None


def all_minimum_supports_reduced(edges, observed, cap, seconds=None,
                                 max_states=None, max_bases=None,
                                 max_faces=None, max_oracle=400000,
                                 stop_after=None, full_mask=False,
                                 mode='combined', kstar=None, screen_rows=0,
                                 distances=None, minimum_knobs=None):
    """Complete set of minimum-size strict source supports on the ORIGINAL tree.

    Signature and return shape follow `stage292_all_minimisers.
    all_minimum_supports`; `mode` selects which pass computes the minimum.
    """
    if mode not in MODES:
        raise ValueError('mode must be one of %r' % (MODES,))
    E = [(u, v, F(L)) for u, v, L in edges]
    observed = {int(r): F(d) for r, d in observed.items()}
    Rlist = sorted(observed)
    stats = _new_stats()
    budget = Budget(seconds, max_states, max_bases, max_faces)
    out = dict(status='complete', minimum=None, supports=[], n_supports=0,
               statistics=stats, reason=None, mode=mode, full_mask=bool(full_mask))
    timing = {}
    try:
        B = observed_backbone(E, Rlist)
        marked = {x['vertex'] for x in B['saturation_violations']}
        knobs = dict(seconds=seconds, max_states=max_states,
                     max_bases=max_bases, max_faces=max_faces)
        if minimum_knobs:
            knobs.update(minimum_knobs)
        if kstar is None and mode != 'combined':
            t0 = time.perf_counter()
            if mode == 'frozen':
                first = solve_partial(E, observed, cap=cap, **knobs)
            else:
                first = solve_partial_reduced(E, observed, cap=cap,
                                              full_mask=full_mask, **knobs)
            timing['minimum_pass_seconds'] = time.perf_counter() - t0
            out['minimum_pass_status'] = first.get('status')
            out['minimum_pass_statistics'] = first.get('statistics')
            if first.get('status') != 'complete':
                out.update(status=('frozen_' if mode == 'frozen' else 'reduced_')
                           + str(first.get('status')), reason=first.get('reason'))
                out['timing'] = timing
                return out
            kstar = first.get('minimum_sources_within_cap')
            out['first_pass_minimum'] = kstar
            out['first_pass_reason'] = first.get('reason')
            if mode == 'frozen':
                out['frozen_minimum'] = kstar
            else:
                out['reduced_minimum'] = kstar
            if kstar is None:
                out['reason'] = out.get('first_pass_reason') or 'infeasible_within_cap'
                out['timing'] = timing
                return out
        t0 = time.perf_counter()
        pats, minimum, reason = minimum_patterns(B['edges'], observed, marked, cap,
                                                 budget, stats,
                                                 full_mask=full_mask, kstar=kstar)
        timing['pattern_pass_seconds'] = time.perf_counter() - t0
        if reason is not None or not pats:
            out.update(minimum=None, reason=reason or
                       ('no_pattern_at_frozen_minimum' if kstar is not None
                        else 'infeasible_within_cap'),
                       status='disagreement' if (reason is None and kstar is not None)
                       else 'complete')
            out['timing'] = timing
            return out
        edge_atom, vertex_atom = lift_tables(B)
        configs = set()
        for pat in pats:
            configs |= all_min_covers(pat['P'], pat['inc'], marked, minimum)
        stats['configurations'] = len(configs)
        if distances is None:
            distances = VertexDistances(E)
        cache = {}
        supports = set()
        screen = Rlist[:screen_rows] if screen_rows else None
        if screen is not None and len(screen) >= len(Rlist):
            screen = None
        t0 = time.perf_counter()
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
        timing['lift_pass_seconds'] = time.perf_counter() - t0
        stats['bfs_calls'] = distances.bfs_calls
        out.update(minimum=minimum, supports=sorted(list(s) for s in supports),
                   n_supports=len(supports))
        if out['status'] == 'complete' and not supports:
            out['status'] = 'disagreement'
            out['reason'] = 'no_lift_confirmed_at_minimum'
    except ResourceLimit as exc:
        out.update(status='resource_limit', reason=str(exc))
    out['timing'] = timing
    out['budget_counts'] = dict(budget.counts)
    return out


def unique_true(result, true_support, K):
    """Stage290's `ok` predicate, recomputed from an all-minimisers result."""
    return bool(result['status'] == 'complete' and result['minimum'] == K and
                result['n_supports'] == 1 and
                list(result['supports'][0]) == list(true_support))
