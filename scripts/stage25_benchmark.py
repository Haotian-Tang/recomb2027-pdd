#!/usr/bin/env python3
"""Same-task exact <=3 benchmark: structural solver vs support enumeration.

The independent baseline uses rational equality reduction and vertices of each
fixed-support nonnegative endpoint polytope. It does not use deficit localization.
Its wall-time limit yields unknown, never infeasible. Main writes Stage25 only.
"""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
from datetime import datetime, timezone
import argparse, hashlib, json, platform, sys, time
from stage24_three_source import solve_three

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def geometry(edges):
    adj = {}
    for u, v, length in edges:
        assert u != v and length > 0
        adj.setdefault(u, {})[v] = F(length)
        adj.setdefault(v, {})[u] = F(length)
    leaves = sorted(v for v in adj if len(adj[v]) == 1)
    dist = {}
    for leaf in leaves:
        d = {leaf: F(0)}
        todo = [leaf]
        for v in todo:
            for w, length in adj[v].items():
                if w not in d:
                    d[w] = d[v] + length
                    todo.append(w)
        assert len(d) == len(adj)
        dist[leaf] = d
    return adj, leaves, dist


def reduce_equalities(rows, variables, stats=None):
    """Streaming exact RREF; None certifies a contradictory equality row."""
    basis = {}
    for row in rows:
        if stats is not None:
            stats['equality_rows_processed'] += 1
        row = list(row)
        for pivot in sorted(basis):
            if row[pivot]:
                a = row[pivot]
                row = [x - a * y for x, y in zip(row, basis[pivot])]
        pivot = next((i for i in range(variables) if row[i]), None)
        if pivot is None:
            if row[-1]:
                return None
            continue
        a = row[pivot]
        row = [x / a for x in row]
        for old in basis:
            if basis[old][pivot]:
                a = basis[old][pivot]
                basis[old] = [x - a * y for x, y in zip(basis[old], row)]
        basis[pivot] = row
    return [basis[p] for p in sorted(basis)]


def unique_solution(rows, variables):
    reduced = reduce_equalities(rows, variables)
    if reduced is None or len(reduced) != variables:
        return None
    return [row[-1] for row in reduced]


def fixed_support_polytope(endpoints, leaves, dist, observed, stats=None):
    """All endpoint masses and pendant strictly positive, or None.

    The nonnegative equality polytope is bounded: masses sum to one and pendant
    is bounded by every finite leaf observation. Strict feasibility holds iff
    each coordinate is positive at some vertex; their average is a witness.
    """
    q = len(endpoints) + 1
    def rows():
        yield [F(1)] * len(endpoints) + [F(0), F(1)]
        for r in leaves:
            yield [dist[r][v] for v in endpoints] + [F(1), observed[r]]
    reduced = reduce_equalities(rows(), q, stats)
    if reduced is None:
        return None
    if len(reduced) == q:
        x = [row[-1] for row in reduced]
        return x if all(v > 0 for v in x) else None
    rank = len(reduced)
    aggregate = [F(0)] * q
    accepted = 0
    for columns in combinations(range(q), rank):
        if stats is not None:
            stats['vertex_bases_attempted'] += 1
        equations = [[row[j] for j in columns] + [row[-1]] for row in reduced]
        basic = unique_solution(equations, rank)
        if basic is None or any(v < 0 for v in basic):
            continue
        accepted += 1
        for j, value in zip(columns, basic):
            aggregate[j] += value
        if all(v > 0 for v in aggregate):
            x = [v / accepted for v in aggregate]
            assert all(sum(a*b for a, b in zip(row[:-1], x)) == row[-1]
                       for row in reduced)
            return x
    return None


def split_witness(edges, support, endpoints, x):
    mass = dict(zip(endpoints, x[:-1]))
    degree = {v: sum(v in edges[i][:2] for i in support) for v in endpoints}
    sources = []
    for i in support:
        u, v, _ = edges[i]
        a, b = mass[u] / degree[u], mass[v] / degree[v]
        sources.append({'edge': i, 'p': a+b, 'x': b/(a+b)})
    return {'feasible': True, 'minimum_sources_within_cap': len(support),
            'support': list(support), 'sources': sources, 'pendant': x[-1]}


def enumerate_three(edges, observed, leaves, dist, time_limit_seconds=10.0):
    start = time.perf_counter()
    stats = dict(supports_visited=0, endpoint_sets_solved=0,
                 endpoint_cache_hits=0, equality_rows_processed=0,
                 vertex_bases_attempted=0)
    cache = {}
    for k in range(1, 4):
        for support in combinations(range(len(edges)), k):
            if time.perf_counter() - start > time_limit_seconds:
                return {'status': 'time_limit', 'completed': False,
                        'last_source_count': k, 'statistics': stats,
                        'seconds': time.perf_counter()-start,
                        'scope': 'Unknown; not an infeasibility certificate.'}
            stats['supports_visited'] += 1
            U = tuple(sorted({v for i in support for v in edges[i][:2]}))
            if U in cache:
                stats['endpoint_cache_hits'] += 1
                x = cache[U]
            else:
                stats['endpoint_sets_solved'] += 1
                x = fixed_support_polytope(U, leaves, dist, observed, stats)
                cache[U] = x
            if x is not None:
                return {'status': 'complete', 'completed': True,
                        'solution': split_witness(edges, support, U, x),
                        'statistics': stats, 'seconds': time.perf_counter()-start}
    return {'status': 'complete', 'completed': True,
            'solution': {'feasible': False}, 'statistics': stats,
            'seconds': time.perf_counter()-start}


def verify_witness(edges, observed, dist, solution):
    if not solution['feasible']:
        return
    sources = solution['sources']
    assert 1 <= len(sources) <= 3
    assert len({s['edge'] for s in sources}) == len(sources)
    assert sorted(s['edge'] for s in sources) == solution['support']
    assert len(sources) == solution['minimum_sources_within_cap']
    assert sum(s['p'] for s in sources) == 1 and solution['pendant'] > 0
    assert all(s['p'] > 0 and 0 < s['x'] < 1 for s in sources)
    for r, value in observed.items():
        predicted = solution['pendant']
        for s in sources:
            u, v, _ = edges[s['edge']]
            predicted += s['p']*((1-s['x'])*dist[r][u]+s['x']*dist[r][v])
        assert predicted == value


def generated_observations(edges, dist, sources, pendant):
    return {r: pendant+sum(s['p']*((1-s['x'])*d[edges[s['edge']][0]]+
                                  s['x']*d[edges[s['edge']][1]]) for s in sources)
            for r, d in dist.items()}


def caterpillar(n):
    assert n >= 3
    if n == 3:
        pairs = [(i, 3) for i in range(3)]
    else:
        internal = list(range(n, 2*n-2))
        pairs = [(0, internal[0]), (1, internal[0])]
        pairs += [(i+2, v) for i, v in enumerate(internal[1:-1])]
        pairs += [(n-2, internal[-1]), (n-1, internal[-1])]
        pairs += list(zip(internal, internal[1:]))
    return [(u, v, F(2+i % 7, 1+i % 3)) for i, (u, v) in enumerate(pairs)]


def balanced(n):
    edges = [(0, 3, F(1)), (1, 3, F(2)), (2, 3, F(3))]
    depth = {0: 1, 1: 1, 2: 1}
    next_id = 4
    while len(depth) < n:
        leaf = min(depth, key=lambda v: (depth[v], v))
        i = next(i for i, (u, v, _) in enumerate(edges) if leaf in (u, v))
        u, v, old_length = edges.pop(i)
        parent = v if u == leaf else u
        new_internal, new_leaf = next_id, next_id+1
        next_id += 2
        edges.extend([(parent, new_internal, old_length),
                      (leaf, new_internal, F(1+leaf % 5, 2)),
                      (new_leaf, new_internal, F(1+new_leaf % 7, 3))])
        d = depth[leaf]+1
        depth[leaf] = depth[new_leaf] = d
    return edges


def compact_solution(sol):
    keys = ['feasible', 'minimum_sources_within_cap', 'support', 'sources',
            'pendant', 'deficit_vertex', 'deficit', 'candidate_potentials_examined', 'reason']
    return {k: sol[k] for k in keys if k in sol}


def cases():
    out = []
    for shape, n, k in [('caterpillar', 3, 3), ('caterpillar', 6, 1),
                         ('balanced', 8, 2), ('caterpillar', 8, 3),
                         ('balanced', 12, 3), ('caterpillar', 16, 3)]:
        edges = caterpillar(n) if shape == 'caterpillar' else balanced(n)
        _, _, dist = geometry(edges)
        m = len(edges)
        selected = [m-1] if k == 1 else ([m//3, m-1] if k == 2 else ([0, 1, 2] if m == 3 else [m//3, 2*m//3, m-1]))
        weights = [F(i+1, k*(k+1)//2) for i in range(k)]
        sources = [dict(edge=e, p=w, x=F(2+i, 7)) for i, (e, w) in enumerate(zip(selected, weights))]
        y = F(2, 5)
        obs = generated_observations(edges, dist, sources, y)
        out.append(dict(name=f'{shape}_{n}_generated_{k}', kind='synthetic',
                        edges=edges, observed=obs, generating_sources=sources, generating_pendant=y))
    # Three-arm supports contain one hidden degree of freedom: exercise vertex
    # averaging, both with a strict witness and at the y=0 closure boundary.
    for value, label in [(F(2), 'strict_three'), (F(1), 'closed_boundary_only')]:
        out.append(dict(name=f'unit_star_3_{label}', kind='synthetic_affine_polytope_control',
                        edges=[(i, 3, F(1)) for i in range(3)], observed={i: value for i in range(3)}))
    # Equal observations admit a closed centre atom, but no <=3 strict support.
    for n in (4, 8, 12):
        out.append(dict(name=f'unit_star_{n}_requires_{n}', kind='synthetic_strict_obstruction',
                        edges=[(i, n, F(1)) for i in range(n)], observed={i: F(2) for i in range(n)}))
    saved = json.loads((ROOT/'results/stage24_public_generating.json').read_text())
    for c in saved['cases']:
        out.append(dict(name=c['case'], kind='frozen_public_simulator_generating_geometry',
                        edges=[(u, v, F(l)) for u, v, l in c['edges']],
                        observed={int(r): F(d) for r, d in c['observed'].items()},
                        source_commit=saved['source_commit'], branch_transform=c['positive_ANI_transform'],
                        prior_scope=saved['scope']))
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--baseline-seconds', type=float, default=10.0)
    parser.add_argument('--only', default=None)
    parser.add_argument('--output', default=str(ROOT/'results/stage25_benchmark.json'))
    args = parser.parse_args()
    records = []
    for c in cases():
        if args.only and args.only not in c['name']:
            continue
        edges, observed = c['edges'], c['observed']
        begun = time.perf_counter()
        adj, leaves, dist = geometry(edges)
        preprocessing = time.perf_counter()-begun
        begun = time.perf_counter()
        fast = solve_three(edges, observed)
        fast_seconds = time.perf_counter()-begun
        verify_witness(edges, observed, dist, fast)
        slow = enumerate_three(edges, observed, leaves, dist, args.baseline_seconds)
        if slow['completed']:
            verify_witness(edges, observed, dist, slow['solution'])
            assert slow['solution']['feasible'] == fast['feasible']
            if fast['feasible']:
                assert slow['solution']['minimum_sources_within_cap'] == fast['minimum_sources_within_cap']
        records.append(c | {'leaf_count': len(leaves), 'vertex_count': len(adj),
                            'edge_count': len(edges), 'baseline_geometry_preprocessing_seconds': preprocessing,
                            'structural_solver_seconds': fast_seconds, 'structural_solution': compact_solution(fast),
                            'baseline': slow})
        print(json.dumps({'case': c['name'], 'leaves': len(leaves), 'fast_seconds': fast_seconds,
                          'baseline_seconds': slow['seconds'], 'baseline_status': slow['status'],
                          'minimum': fast.get('minimum_sources_within_cap'), 'supports': slow['statistics']['supports_visited']}), flush=True)
    output = dict(created_utc=datetime.now(timezone.utc).isoformat(), script_sha256=sha(Path(__file__)),
                  dependencies={str(p.relative_to(ROOT)): sha(p) for p in [ROOT/'scripts/stage24_three_source.py', ROOT/'scripts/stage23_two_source_envelope.py', ROOT/'results/stage24_public_generating.json']},
                  python=sys.version, platform=platform.platform(), baseline_limit_seconds=args.baseline_seconds,
                  baseline_order='increasing source count then lexicographic input-edge indices; no true support supplied',
                  timing_scope='Single-run wall times; structural timer includes its own adjacency/envelope preprocessing. Enumeration timer excludes separately recorded leaf-distance preprocessing. Witness checks excluded from both. No runtime ratio for censored baseline.',
                  baseline_method='Streaming rational equality reduction; fixed-support compact nonnegative-polytope vertex averaging; no deficit theorem, no floating LP, caches identical endpoint sets.',
                  claim_scope='Pilot implementations of the same exact <=3 minimum problem, not a comparison to optimized DecoDiPhy or a biological benchmark. Public cases use frozen rational generating geometry, not estimated observations.',
                  cases=records, completed_comparisons=sum(c['baseline']['completed'] for c in records),
                  censored_comparisons=sum(not c['baseline']['completed'] for c in records), all_completed_comparisons_agree=True)
    Path(args.output).write_text(json.dumps(output, indent=2, default=str)+'\n')


if __name__ == '__main__':
    main()
