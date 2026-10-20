#!/usr/bin/env python3
"""Preregistered, bounded same-task exact PDD prototype comparison.

Original-support baseline uses no projection/envelope/deficit localization.
Frozen Stage41 solver is imported read-only. No result is a recovery claim.
"""
import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
import argparse
from collections import Counter
from datetime import datetime, timezone
from fractions import Fraction as F
import hashlib
from itertools import combinations
import json
from pathlib import Path
import platform
import random
import resource
import signal
import sys
import tempfile
import time
import traceback
from stage25_benchmark import balanced, caterpillar, geometry, fixed_support_polytope, split_witness
from stage41_partial_fpt import solve_partial

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT / 'results/stage42_benchmark_protocol.json'
RESULT = ROOT / 'results/stage42_benchmark.json'
DEPENDENCIES = ['stage42_benchmark.py', 'stage41_partial_fpt.py', 'stage27_general_fpt.py',
                'stage25_benchmark.py', 'stage26_missing_leaf.py', 'stage23_two_source_envelope.py',
                'stage25_deficit_extension.py', 'stage24_three_source.py']


def serial(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {str(k): serial(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [serial(x) for x in value]
    return value


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def independent_geometry(edges):
    adjacency = {}
    for u, v, length in edges:
        adjacency.setdefault(u, []).append((v, length))
        adjacency.setdefault(v, []).append((u, length))
    distance = {}
    for origin in adjacency:
        found = {origin: F(0)}
        queue = [origin]
        for vertex in queue:
            for neighbor, length in adjacency[vertex]:
                if neighbor not in found:
                    found[neighbor] = found[vertex] + length
                    queue.append(neighbor)
        assert len(found) == len(adjacency)
        distance[origin] = found
    return adjacency, distance


def observations(edges, sources, pendant, observed):
    _, distance = independent_geometry(edges)
    return {r: pendant + sum(s['p'] * min(
        distance[r][edges[s['edge']][0]] + s['x'] * edges[s['edge']][2],
        distance[r][edges[s['edge']][1]] + (1 - s['x']) * edges[s['edge']][2])
        for s in sources) for r in observed}


def verify(case, witness):
    """Check original-tree point distances, independently of both solvers."""
    sources = witness['sources']
    y = F(witness['pendant'])
    sources = [dict(edge=s['edge'], p=F(s['p']), x=F(s['x'])) for s in sources]
    assert 1 <= len(sources) <= case['cap']
    assert len({s['edge'] for s in sources}) == len(sources)
    assert sum(s['p'] for s in sources) == 1 and y > 0
    assert all(s['p'] > 0 and 0 < s['x'] < 1 for s in sources)
    assert observations(case['edges'], sources, y, case['observed']) == case['observed']
    return len(sources)


def fixtures():
    out = []
    rng = random.Random(421014)
    for shape_index, (shape, build) in enumerate([('caterpillar', caterpillar), ('balanced', balanced)]):
        for size_index, n in enumerate((6, 10, 16)):
            edges = build(n)
            adjacency, _ = independent_geometry(edges)
            leaves = sorted(v for v in adjacency if len(adjacency[v]) == 1)
            for k in (2, 3, 4):
                selected_leaves = [leaves[(j * (n - 1)) // (k - 1)] for j in range(k)]
                selected = [next(i for i, e in enumerate(edges) if v in e[:2]) for v in selected_leaves]
                weights = [F(j + 1, k * (k + 1) // 2) for j in range(k)]
                sources = [dict(edge=e, p=p, x=F(2 + j, 7)) for j, (e, p) in enumerate(zip(selected, weights))]
                fraction = (F(0), F(1, 4), F(1, 2))[(shape_index + size_index + k) % 3]
                hidden_count = min(n - k, int(n * fraction))
                available = [r for r in leaves if r not in selected_leaves]
                rng.shuffle(available)
                hidden = set(available[:hidden_count])
                observed = [r for r in leaves if r not in hidden]
                out.append(dict(name=f'{shape}_n{n}_k{k}_missing{hidden_count}', group='binary_generated',
                    edges=edges, observed=observations(edges, sources, F(2, 5), observed), cap=k,
                    generating_sources=sources, generating_pendant=F(2, 5),
                    requested_missing_fraction=fraction, expected_minimum=None))
    # Replace two redundant small generated cases before preregistration with
    # an audited minimum-four fixture requiring two nonzero deficit vertices.
    out = [c for c in out if not c['name'].startswith(('caterpillar_n6_k4_', 'balanced_n6_k4_'))]
    alternating = [(100, 0, F(1)), (100, 1, F(1)), (101, 2, F(1)),
                   (102, 3, F(1)), (102, 4, F(1)), (100, 101, F(1)), (101, 102, F(1))]
    sources = [dict(edge=i, p=F(1, 4), x=F(1, 2)) for i in (0, 1, 3, 4)]
    values = observations(alternating, sources, F(1), range(5))
    out.append(dict(name='alternating4_full', group='essential_two_deficits', edges=alternating,
                    observed=values, cap=4, expected_minimum=4, expected_feasible=True,
                    known_generating_deficit_vertices=[100, 102], generating_sources=sources,
                    generating_pendant=F(1)))
    changed = alternating[:5] + [(100, 2000, F(1, 2)), (2000, 101, F(1, 2)),
                                (2000, 2001, F(1)), alternating[6]]
    out.append(dict(name='alternating4_one_missing', group='essential_two_deficits', edges=changed,
                    observed=values, cap=4, expected_minimum=4, expected_feasible=True,
                    known_generating_deficit_vertices=[100, 102], generating_sources=sources,
                    generating_pendant=F(1)))
    for degree, cap in ((3, 2), (3, 3), (4, 3), (4, 4), (6, 3), (6, 4)):
        edges = [(100, i, F(1)) for i in range(degree)]
        out.append(dict(name=f'star{degree}_strict_cap{cap}', group='known_strict_budget_boundary',
                        edges=edges, observed={i: F(2) for i in range(degree)}, cap=cap,
                        expected_minimum=degree, expected_feasible=degree <= cap))
    for degree in (3, 4):
        out.append(dict(name=f'star{degree}_closure_only', group='known_closure_boundary',
                        edges=[(100, i, F(1)) for i in range(degree)],
                        observed={i: F(1) for i in range(degree)}, cap=degree,
                        expected_minimum=None, expected_feasible=False))
    for hubs, cap in ((2, 2), (2, 3), (3, 3), (4, 4)):
        edges, observed, hidden_edges = [], [], []
        for j in range(hubs):
            hub = 100 + 2 * j
            for t in (0, 1):
                leaf = 10 * j + t
                observed.append(leaf)
                edges.append((hub, leaf, F(2 + j + t, 2)))
            hidden_edges.append(len(edges))
            edges.append((hub, 10 * j + 2, F(3, 2)))
            if j < hubs - 1:
                middle = hub + 1
                edges.extend([(hub, middle, F(1)), (middle, hub + 2, F(4, 3)),
                              (middle, 10 * j + 3, F(1))])
                observed.append(10 * j + 3)
        sources = [dict(edge=e, p=F(j + 1, hubs * (hubs + 1) // 2), x=F(1, 2))
                   for j, e in enumerate(hidden_edges)]
        out.append(dict(name=f'marked_hubs{hubs}_cap{cap}', group='multifurcating_marked_gates',
                        edges=edges, observed=observations(edges, sources, F(1), observed), cap=cap,
                        generating_sources=sources, generating_pendant=F(1), expected_minimum=None))
    assert len(out) == 30
    for case in out:
        adjacency, _ = independent_geometry(case['edges'])
        case['vertices'] = len(adjacency)
        case['original_leaves'] = sum(len(n) == 1 for n in adjacency.values())
        case['observed_leaves'] = len(case['observed'])
        case['missing_fraction'] = F(case['original_leaves'] - len(case['observed']), case['original_leaves'])
        assert all(len(n) != 2 for n in adjacency.values())
    return out


def baseline(case):
    """Increasing original-edge supports; rational endpoint LPs and caching."""
    edges, obs = case['edges'], case['observed']
    _, _, distance = geometry(edges)
    stats = dict(supports_visited=0, endpoint_sets_solved=0, endpoint_cache_hits=0,
                 equality_rows_processed=0, vertex_bases_attempted=0)
    cache = {}
    for k in range(1, case['cap'] + 1):
        for support in combinations(range(len(edges)), k):
            stats['supports_visited'] += 1
            endpoints = tuple(sorted({v for i in support for v in edges[i][:2]}))
            if endpoints not in cache:
                stats['endpoint_sets_solved'] += 1
                cache[endpoints] = fixed_support_polytope(endpoints, sorted(obs), distance, obs, stats)
            else:
                stats['endpoint_cache_hits'] += 1
            if cache[endpoints] is not None:
                witness = split_witness(edges, support, endpoints, cache[endpoints])
                return dict(status='complete', feasible=True, minimum=k, witness=witness, statistics=stats)
    return dict(status='complete', feasible=False, minimum=None, statistics=stats)


def worker(method, case):
    start = time.perf_counter()
    if method == 'original_support':
        result = baseline(case)
    else:
        raw = solve_partial(case['edges'], case['observed'], cap=case['cap'], seconds=None,
                            max_states=None, max_bases=None, max_faces=None)
        result = dict(status=raw['status'], feasible=raw.get('feasible'),
                      minimum=raw.get('minimum_sources_within_cap'), statistics=raw['statistics'],
                      marked_vertices=raw['marked_vertices'], backbone_edges=len(raw['backbone_edges']),
                      deficit_vertices=raw.get('deficit_vertices'), deficits=raw.get('deficits'))
        if 'original_witness' in raw:
            result['witness'] = raw['original_witness']
    result['solver_wall_seconds'] = time.perf_counter() - start
    return result


def isolated_run(method, case, limit, memory_bytes):
    """Linux fork + wait4 records actual child CPU, including killed jobs."""
    with tempfile.TemporaryDirectory(prefix='stage42_bench_') as temporary:
        result_path = Path(temporary) / 'result.json'
        start = time.perf_counter()
        pid = os.fork()
        if pid == 0:
            try:
                resource.setrlimit(resource.RLIMIT_AS, (memory_bytes, memory_bytes))
                answer = worker(method, case)
            except BaseException as exc:
                answer = dict(status='error', exception=repr(exc), traceback=traceback.format_exc())
            try:
                result_path.write_text(json.dumps(serial(answer)))
            finally:
                os._exit(0)
        timed_out = False
        while True:
            waited, status, usage = os.wait4(pid, os.WNOHANG)
            if waited:
                break
            if time.perf_counter() - start >= limit:
                timed_out = True
                os.kill(pid, signal.SIGKILL)
                _, status, usage = os.wait4(pid, 0)
                break
            time.sleep(0.01)
        elapsed = time.perf_counter() - start
        if timed_out:
            answer = dict(status='unknown', reason='hard_wall_timeout', minimum=None)
        elif result_path.exists():
            answer = json.loads(result_path.read_text())
        else:
            answer = dict(status='error', reason='worker_exited_without_result')
        answer.update(process_wall_seconds=elapsed, child_cpu_seconds=usage.ru_utime + usage.ru_stime,
                      maximum_rss_kib=usage.ru_maxrss, wait_status=status)
        if 'witness' in answer:
            answer['independent_witness_verified'] = verify(case, answer['witness']) == answer['minimum']
            assert answer['independent_witness_verified']
        return answer


def prepare():
    assert not PROTOCOL.exists(), 'Protocol is immutable after preregistration; do not silently overwrite'
    protocol = dict(stage=42, registered_utc=datetime.now(timezone.utc).isoformat(),
        description='Fixed deterministic exact minimum-source tasks; prototype limits, not biological recovery',
        per_job_hard_wall_seconds=5.0, per_job_address_space_bytes=1024**3,
        process_isolation='fresh Linux fork; preloaded imports shared equally; sequential jobs; wait4 CPU/RSS',
        thread_environment={key: os.environ[key] for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS')},
        method_order='alternate original_support then marked_deficit by even/odd case index',
        analysis='All cases retained; unknown is never infeasible. Compare minima only where both complete. No timeout imputation in speed ratios.',
        scope_caveats=['synthetic small exact tasks', 'one timed repetition, no uncertainty intervals',
            'baseline is exact endpoint support enumeration, not DecoDiPhy',
            'frozen partial prototype contains sorting and O(KN) witness-check overhead',
            'imports and case generation excluded; independent witness verification after timing',
            'memory limit is address space, not peak RSS; both methods receive identical limit'],
        case_list=fixtures(), script_hashes={n: sha(ROOT / 'scripts' / n) for n in DEPENDENCIES},
        machine=dict(platform=platform.platform(), python=sys.version, processor=platform.processor(),
                     cpu_count=os.cpu_count()))
    PROTOCOL.write_text(json.dumps(serial(protocol), indent=2) + '\n')
    print(json.dumps(dict(protocol=str(PROTOCOL), sha256=sha(PROTOCOL), cases=len(protocol['case_list'])), indent=2))


def run():
    assert not RESULT.exists(), 'Do not overwrite a measured run'
    protocol = json.loads(PROTOCOL.read_text())
    for n, expected in protocol['script_hashes'].items():
        assert sha(ROOT / 'scripts' / n) == expected, f'Preregistered code changed: {n}'
    records = []
    start = time.perf_counter()
    for index, encoded_case in enumerate(protocol['case_list']):
        case = dict(encoded_case)
        case['edges'] = [(u, v, F(length)) for u, v, length in case['edges']]
        case['observed'] = {int(r): F(value) for r, value in case['observed'].items()}
        methods = ['original_support', 'marked_deficit']
        if index % 2:
            methods.reverse()
        record = dict(name=case['name'], group=case['group'], cap=case['cap'],
                      vertices=case['vertices'], edges=len(case['edges']), original_leaves=case['original_leaves'],
                      observed_leaves=case['observed_leaves'], missing_fraction=case['missing_fraction'], runs={})
        for method in methods:
            record['runs'][method] = isolated_run(method, case, protocol['per_job_hard_wall_seconds'],
                                                  protocol['per_job_address_space_bytes'])
        a, b = [record['runs'][method] for method in ('original_support', 'marked_deficit')]
        record['both_complete'] = a['status'] == b['status'] == 'complete'
        if record['both_complete']:
            record['exact_agreement'] = (a['feasible'], a['minimum']) == (b['feasible'], b['minimum'])
            assert record['exact_agreement'], (case['name'], a, b)
        for method, measured in record['runs'].items():
            if measured['status'] == 'complete' and 'expected_feasible' in case:
                assert measured['feasible'] == case['expected_feasible']
                if measured['feasible']:
                    assert measured['minimum'] == case['expected_minimum']
        records.append(record)
        print(f"{index + 1:02d}/30 {case['name']}: " + ' | '.join(
            f"{method}={value['status']} min={value.get('minimum')} wall={value['process_wall_seconds']:.3f}s"
            for method, value in record['runs'].items()), flush=True)
    summary = dict(cases=len(records), paired_complete=sum(r['both_complete'] for r in records),
        paired_disagreements=sum(r.get('exact_agreement') is False for r in records),
        methods={method: dict(status_counts=dict(Counter(r['runs'][method]['status'] for r in records)),
            certified_minimum_counts=dict(Counter(str(r['runs'][method]['minimum']) for r in records
                if r['runs'][method]['status'] == 'complete')),
            verified_witnesses=sum(r['runs'][method].get('independent_witness_verified', False) for r in records),
            total_cpu_seconds=sum(r['runs'][method]['child_cpu_seconds'] for r in records),
            total_process_wall_seconds=sum(r['runs'][method]['process_wall_seconds'] for r in records))
            for method in ('original_support', 'marked_deficit')},
        total_wall_seconds=time.perf_counter() - start)
    result = dict(stage=42, created_utc=datetime.now(timezone.utc).isoformat(), protocol_sha256=sha(PROTOCOL),
                  script_hashes=protocol['script_hashes'], records=records, summary=summary)
    RESULT.write_text(json.dumps(serial(result), indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('phase', choices=['prepare', 'run'])
    args = parser.parse_args()
    prepare() if args.phase == 'prepare' else run()
