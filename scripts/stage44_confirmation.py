#!/usr/bin/env python3
"""Preregistered new-instance confirmation, with no solver changes or tuning.

Prepare writes all exact inputs and transitive code hashes before any solve.
Run retains every case/repetition and every timeout; analysis never imputes timeouts.
"""
import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
import argparse
import ast
from collections import Counter
from datetime import datetime, timezone
from fractions import Fraction as F
from itertools import combinations
import json
from pathlib import Path
import platform
import random
import resource
import signal
from statistics import median
import sys
import tempfile
import time
import traceback
from stage42_benchmark import baseline, independent_geometry, observations, serial, sha, verify
from stage41_partial_fpt import solve_partial as old_solve
from stage43_partial_fpt import solve_partial as new_solve
from stage43_fixed_one_missing import solve_fixed_one_missing
from check_stage43_fixed_one_missing import oracle

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT / 'results/stage44_confirmation_protocol.json'
RESULT = ROOT / 'results/stage44_confirmation.json'
STREAM = ROOT / 'results/stage44_confirmation_runs.jsonl'
FREE = ('original_support', 'stage41_active_basis', 'stage43_equality_reduced')
FIXED = ('original_fixed_oracle', 'stage43_fixed_one')
SEED = 4410142026


def random_binary(n, rng):
    """Uniform edge insertion from a labelled three-leaf star; random orientation."""
    edges = [(n, i) for i in range(3)]
    next_internal = n + 1
    for leaf in range(3, n):
        at = rng.randrange(len(edges))
        u, v = edges.pop(at)
        edges.extend([(u, next_internal), (next_internal, v), (next_internal, leaf)])
        next_internal += 1
    rng.shuffle(edges)
    result = []
    for j, (u, v) in enumerate(edges):
        if rng.randrange(2):
            u, v = v, u
        result.append((u, v, F(rng.randint(2, 19), rng.randint(2, 9))))
    assert len({e[2] for e in result}) > 1
    return result


def source_sample(edges, count, mode, rng, forced=None, forbidden=None):
    adjacency, _ = independent_geometry(edges)
    terminal = [i for i, (u, v, _) in enumerate(edges) if len(adjacency[u]) == 1 or len(adjacency[v]) == 1]
    internal = [i for i in range(len(edges)) if i not in terminal]
    eligible = [i for i in range(len(edges)) if i != forbidden]
    chosen = [] if forced is None else [forced]
    if mode == 'mixed':
        for pool in (internal, terminal):
            valid = [i for i in pool if i in eligible and i not in chosen]
            if valid and len(chosen) < count:
                chosen.append(rng.choice(valid))
    elif mode == 'terminal':
        valid = [i for i in terminal if i in eligible and i not in chosen]
        chosen.extend(rng.sample(valid, min(count - len(chosen), len(valid))))
    valid = [i for i in eligible if i not in chosen]
    chosen.extend(rng.sample(valid, count - len(chosen)))
    rng.shuffle(chosen)
    weights = [rng.randint(1, 9) for _ in chosen]
    return [dict(edge=e, p=F(w, sum(weights)), x=F(rng.randint(1, 10), 11))
            for e, w in zip(chosen, weights)]


def fixtures():
    cases = []
    rng = random.Random(SEED)
    for i in range(32):
        n = (6, 8, 10, 12)[(i // 4) % 4]
        tightened = i >= 16
        cap = 2 + i % 4
        generating_count = cap + int(tightened)
        edges = random_binary(n, rng)
        adjacency, _ = independent_geometry(edges)
        leaves = sorted(v for v in adjacency if len(adjacency[v]) == 1)
        mode = 'terminal' if i % 3 == 0 else 'mixed'
        sources = source_sample(edges, generating_count, mode, rng)
        hidden_count = (0, 1, n // 3, n // 2)[(i // 4 + i) % 4]
        hidden = rng.sample(leaves, hidden_count)
        pendant = F(rng.randint(1, 7), rng.randint(2, 9))
        observed = observations(edges, sources, pendant, [r for r in leaves if r not in hidden])
        cases.append(dict(name=f'free_random_{i:02d}_n{n}_cap{cap}_g{generating_count}_missing{hidden_count}',
            task='free_pendant', group='budget_tightened' if tightened else 'budget_at_generator',
            edges=edges, observed=observed, cap=cap, hidden_leaves=hidden,
            generating_sources=sources, generating_pendant=pendant, source_selection_mode=mode))
    for i in range(16):
        n = (6, 8, 10, 12)[i // 4]
        cap = 2 + i % 4
        tightened = i >= 8
        generating_count = cap + int(tightened)
        edges = random_binary(n, rng)
        adjacency, _ = independent_geometry(edges)
        leaves = sorted(v for v in adjacency if len(adjacency[v]) == 1)
        hidden = rng.choice(leaves)
        hidden_edge = next(j for j, e in enumerate(edges) if hidden in e[:2])
        hidden_present = i % 2 == 0
        sources = source_sample(edges, generating_count, 'mixed', rng,
            forced=hidden_edge if hidden_present else None,
            forbidden=None if hidden_present else hidden_edge)
        pendant = F(rng.randint(1, 7), rng.randint(2, 9))
        observed = observations(edges, sources, pendant, [r for r in leaves if r != hidden])
        cases.append(dict(name=f'fixed_random_{i:02d}_n{n}_cap{cap}_g{generating_count}_hidden{int(hidden_present)}',
            task='fixed_one_missing', group='budget_tightened' if tightened else 'budget_at_generator',
            edges=edges, observed=observed, cap=cap, pendant=pendant, hidden_leaves=[hidden],
            generating_sources=sources, generating_pendant=pendant, source_selection_mode='mixed',
            hidden_source_in_generator=hidden_present))
    fingerprints = set()
    for case in cases:
        adjacency, _ = independent_geometry(case['edges'])
        assert all(len(n) in (1, 3) for n in adjacency.values())
        assert len(case['edges']) == len(adjacency) - 1
        case.update(vertices=len(adjacency), original_leaves=sum(len(n) == 1 for n in adjacency.values()),
                    observed_leaves=len(case['observed']), generating_count=len(case['generating_sources']))
        case['source_edge_types'] = ['terminal' if any(len(adjacency[v]) == 1 for v in case['edges'][s['edge']][:2])
                                     else 'internal' for s in case['generating_sources']]
        case['cherries'] = sum(sum(len(adjacency[u]) == 1 for u, _ in neighbors) == 2
                              for neighbors in adjacency.values() if len(neighbors) == 3)
        case['missing_fraction'] = F(len(case['hidden_leaves']), case['original_leaves'])
        check = dict(case, cap=case['generating_count'])
        assert verify(check, {'sources': case['generating_sources'], 'pendant': case['generating_pendant']}) == case['generating_count']
        key = json.dumps(serial([case['edges'], case['observed'], case['cap'], case.get('pendant')]), sort_keys=True)
        assert key not in fingerprints
        fingerprints.add(key)
    assert len(cases) == 48
    return cases


def dependency_closure():
    pending, result = [Path(__file__).name], set()
    while pending:
        name = pending.pop()
        if name in result:
            continue
        result.add(name)
        parsed = ast.parse((ROOT / 'scripts' / name).read_text())
        for node in ast.walk(parsed):
            names = [node.module] if isinstance(node, ast.ImportFrom) else \
                    [item.name for item in node.names] if isinstance(node, ast.Import) else []
            for module in names:
                if module:
                    candidate = module.split('.')[0] + '.py'
                    if (ROOT / 'scripts' / candidate).exists():
                        pending.append(candidate)
    return sorted(result)


def prepare():
    assert not PROTOCOL.exists() and not RESULT.exists() and not STREAM.exists(), 'Never overwrite frozen evidence'
    cases = fixtures()  # Arithmetic input checks only; no algorithm or oracle is called.
    protocol = dict(stage=44, registered_utc=datetime.now(timezone.utc).isoformat(), seed=SEED,
        purpose='New-instance exact algorithm confirmation; neither biological validation nor novelty evidence',
        generator='Uniform random original-edge insertion from a labelled three-leaf star; every length/orientation/source choice seeded; no output-dependent filtering',
        selection='32 free and 16 fixed-one-missing exact generated cases; no solver calls before protocol freeze; no solver tuning from outcomes',
        repetitions=2, per_job_hard_wall_seconds=2.0, per_job_address_space_bytes=1024**3,
        methods={'free_pendant': FREE, 'fixed_one_missing': FIXED},
        execution_order='Repetition outermost; cases in registered order; rotate method list by (case_index + repetition) modulo method count; reverse if floor(case_index/method_count)+repetition is odd',
        process_isolation='Fresh Linux fork, preloaded imports, sequential jobs, child CPU/RSS via wait4; wall cap includes serialization and all child work',
        timing='Imports and fixture construction excluded; existing solver internal verification included; independent original-tree witness verification outside child timing',
        analysis='Retain all 256 jobs. Timeouts remain UNKNOWN. Check stable exact answers across all completed methods/repetitions. Pairwise median timings use per-case median of two repetitions only when both methods complete both. No timeout imputation or paired speed ratios using censored cases.',
        limitations=['Small synthetic exact instances; minimum-source task, not recovery of generating composition',
            'Two repetitions do not justify confidence intervals or population-level speed claims',
            'Budget below generating source count need not be infeasible; true minimum determined only by completed exact solvers',
            'Original-support oracle is a transparent exact enumerator, not published DecoDiPhy software',
            'All methods are frozen before these new cases are solved; standard equality reduction remains engineering'],
        thread_environment={key: os.environ[key] for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS')},
        cases=cases, script_hashes={n: sha(ROOT / 'scripts' / n) for n in dependency_closure()},
        prior_protocol_hashes={name: sha(ROOT / 'results' / name) for name in
            ['stage43_benchmark_protocol.json', 'stage43_fixed_one_protocol.json', 'stage43_fixed_one_extra_protocol.json']},
        machine=dict(platform=platform.platform(), python=sys.version, cpu_count=os.cpu_count()))
    PROTOCOL.write_text(json.dumps(serial(protocol), indent=2) + '\n')
    print(json.dumps(dict(protocol=str(PROTOCOL), sha256=sha(PROTOCOL), cases=len(cases), jobs=256,
                         frozen_source_files=len(protocol['script_hashes'])), indent=2))


def worker(method, case):
    start = time.perf_counter()
    if method == 'original_support':
        answer = baseline(case)
    elif method == 'original_fixed_oracle':
        raw = oracle(case)
        answer = dict(status='complete', feasible=raw['minimum'] is not None, **raw)
    else:
        options = dict(cap=case['cap'], seconds=None, max_states=None, max_bases=None, max_faces=None)
        if method == 'stage43_fixed_one':
            raw = solve_fixed_one_missing(case['edges'], case['observed'], case['pendant'], **options)
        else:
            solve = old_solve if method == 'stage41_active_basis' else new_solve
            raw = solve(case['edges'], case['observed'], **options)
        answer = dict(status=raw['status'], feasible=raw.get('feasible'),
            minimum=raw.get('minimum_sources_within_cap'), statistics=raw['statistics'],
            deficit_vertices=raw.get('deficit_vertices'), deficits=raw.get('deficits'))
        if method == 'stage43_fixed_one' and 'sources' in raw:
            answer['witness'] = {'sources': raw['sources'], 'pendant': raw['pendant']}
        elif 'original_witness' in raw:
            answer['witness'] = raw['original_witness']
    answer['solver_wall_seconds'] = time.perf_counter() - start
    return answer


def isolated_run(method, case, limit, memory_bytes):
    with tempfile.TemporaryDirectory(prefix='stage44_confirmation_') as temporary:
        path = Path(temporary) / 'result.json'
        start = time.perf_counter()
        pid = os.fork()
        if pid == 0:
            try:
                resource.setrlimit(resource.RLIMIT_AS, (memory_bytes, memory_bytes))
                answer = worker(method, case)
            except BaseException as exc:
                answer = dict(status='error', exception=repr(exc), traceback=traceback.format_exc())
            try:
                path.write_text(json.dumps(serial(answer)))
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
            time.sleep(0.005)
        elapsed = time.perf_counter() - start
        if timed_out:
            answer = dict(status='unknown', reason='hard_wall_timeout', minimum=None)
        elif path.exists():
            answer = json.loads(path.read_text())
        else:
            answer = dict(status='error', reason='worker_exited_without_result')
        answer.update(process_wall_seconds=elapsed, child_cpu_seconds=usage.ru_utime + usage.ru_stime,
                      maximum_rss_kib=usage.ru_maxrss, wait_status=status)
        if 'witness' in answer:
            try:
                answer['independent_witness_verified'] = verify(case, answer['witness']) == answer['minimum']
                if case['task'] == 'fixed_one_missing':
                    answer['prescribed_pendant_verified'] = F(answer['witness']['pendant']) == case['pendant']
                    assert answer['prescribed_pendant_verified']
                assert answer['independent_witness_verified']
            except BaseException as exc:
                answer['independent_witness_verified'] = False
                answer['validation_error'] = repr(exc)
        return answer


def summarize(protocol, records):
    summary = dict(cases=len(records), jobs=sum(len(rep['runs']) for r in records for rep in r['repetitions']), tasks={})
    for task, methods in protocol['methods'].items():
        subset = [r for r in records if r['task'] == task]
        task_result = dict(cases=len(subset), methods={}, pairwise={}, completed_answer_disagreements=0,
                           witness_validation_failures=0)
        for r in subset:
            completed = [v for rep in r['repetitions'] for v in rep['runs'].values() if v['status'] == 'complete']
            r['all_completed_answers_agree'] = len({(v['feasible'], v['minimum']) for v in completed}) <= 1
            task_result['completed_answer_disagreements'] += int(not r['all_completed_answers_agree'])
        for method in methods:
            runs = [rep['runs'][method] for r in subset for rep in r['repetitions']]
            both = [r for r in subset if all(rep['runs'][method]['status'] == 'complete' for rep in r['repetitions'])]
            failures = sum(v.get('independent_witness_verified') is False for v in runs)
            task_result['witness_validation_failures'] += failures
            task_result['methods'][method] = dict(status_counts=dict(Counter(v['status'] for v in runs)),
                cases_complete_both_repetitions=len(both),
                feasible_cases_both_repetitions=sum(r['repetitions'][0]['runs'][method]['feasible'] for r in both),
                minimum_counts_both_repetitions=dict(Counter(str(r['repetitions'][0]['runs'][method]['minimum']) for r in both)),
                verified_witnesses=sum(v.get('independent_witness_verified', False) for v in runs),
                cpu_seconds_all_jobs_including_timeouts=sum(v['child_cpu_seconds'] for v in runs))
        for a, b in combinations(methods, 2):
            joint = [r for r in subset if all(rep['runs'][m]['status'] == 'complete' for rep in r['repetitions'] for m in (a, b))]
            values = {m: [median(rep['runs'][m]['child_cpu_seconds'] for rep in r['repetitions']) for r in joint] for m in (a, b)}
            task_result['pairwise'][f'{a} vs {b}'] = dict(joint_complete_both_repetitions=len(joint),
                median_case_cpu_seconds={m: median(values[m]) if joint else None for m in (a, b)},
                median_case_ratio_first_over_second=median(x / y for x, y in zip(values[a], values[b])) if joint else None,
                first_lower_cpu_cases=sum(x < y for x, y in zip(values[a], values[b])),
                second_lower_cpu_cases=sum(y < x for x, y in zip(values[a], values[b])))
        summary['tasks'][task] = task_result
    summary['status'] = 'PASS' if all(not x['completed_answer_disagreements'] and not x['witness_validation_failures']
        and all(not m['status_counts'].get('error') for m in x['methods'].values()) for x in summary['tasks'].values()) else 'FAIL'
    return summary


def run():
    assert not RESULT.exists() and not STREAM.exists(), 'Never overwrite measured runs'
    protocol = json.loads(PROTOCOL.read_text())
    for name, expected in protocol['script_hashes'].items():
        assert sha(ROOT / 'scripts' / name) == expected, f'Frozen source changed: {name}'
    for name, expected in protocol['prior_protocol_hashes'].items():
        assert sha(ROOT / 'results' / name) == expected
    records = [dict(name=c['name'], task=c['task'], group=c['group'], cap=c['cap'],
                    original_leaves=c['original_leaves'], missing_fraction=c['missing_fraction'],
                    generating_count=c['generating_count'], repetitions=[]) for c in protocol['cases']]
    started = time.perf_counter()
    for repetition in range(protocol['repetitions']):
        for index, encoded in enumerate(protocol['cases']):
            case = dict(encoded)
            case['edges'] = [(u, v, F(length)) for u, v, length in case['edges']]
            case['observed'] = {int(r): F(value) for r, value in case['observed'].items()}
            if 'pendant' in case:
                case['pendant'] = F(case['pendant'])
            methods = list(protocol['methods'][case['task']])
            shift = (index + repetition) % len(methods)
            methods = methods[shift:] + methods[:shift]
            if (index // len(methods) + repetition) % 2:
                methods.reverse()
            rep = dict(repetition=repetition + 1, execution_order=methods, runs={})
            for method in methods:
                value = isolated_run(method, case, protocol['per_job_hard_wall_seconds'], protocol['per_job_address_space_bytes'])
                rep['runs'][method] = value
                with STREAM.open('a') as stream:
                    stream.write(json.dumps(dict(name=case['name'], repetition=repetition + 1, method=method, result=value)) + '\n')
            records[index]['repetitions'].append(rep)
            print(f'rep {repetition + 1} {index + 1:02d}/48 {case["name"]}: ' + ' | '.join(
                f'{m}={rep["runs"][m]["status"]} min={rep["runs"][m].get("minimum")} cpu={rep["runs"][m]["child_cpu_seconds"]:.5f}'
                for m in methods), flush=True)
    summary = summarize(protocol, records)
    summary['total_wall_seconds'] = time.perf_counter() - started
    result = dict(stage=44, created_utc=datetime.now(timezone.utc).isoformat(), protocol_sha256=sha(PROTOCOL),
                  script_hashes=protocol['script_hashes'], records=records, summary=summary)
    RESULT.write_text(json.dumps(serial(result), indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('phase', choices=('prepare', 'run'))
    args = parser.parse_args()
    prepare() if args.phase == 'prepare' else run()
