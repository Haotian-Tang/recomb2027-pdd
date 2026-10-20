#!/usr/bin/env python3
"""Frozen three-way exact benchmark of the unchanged Stage42 30 cases.

New solver optimization was selected from the Stage42 reported bottleneck.
Correctness checks precede preregistration, but timed measurements do not.
"""
import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
import argparse
from collections import Counter
from datetime import datetime, timezone
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import platform
import resource
import signal
import sys
import tempfile
import time
import traceback
from stage42_benchmark import baseline, serial, sha, verify
from stage41_partial_fpt import solve_partial as old_solve
from stage43_partial_fpt import solve_partial as new_solve
ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT / 'results/stage43_benchmark_protocol.json'
RESULT = ROOT / 'results/stage43_benchmark.json'
PARENT_PROTOCOL = ROOT / 'results/stage42_benchmark_protocol.json'
METHODS = ('original_support', 'stage41_active_basis', 'stage43_equality_reduced')


def worker(method, case):
    start = time.perf_counter()
    if method == 'original_support':
        result = baseline(case)
    else:
        solve = old_solve if method == 'stage41_active_basis' else new_solve
        raw = solve(case['edges'], case['observed'], cap=case['cap'], seconds=None,
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
    with tempfile.TemporaryDirectory(prefix='stage43_bench_') as temporary:
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
            time.sleep(0.01)
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
            answer['independent_witness_verified'] = verify(case, answer['witness']) == answer['minimum']
            assert answer['independent_witness_verified']
        return answer


def prepare():
    assert not PROTOCOL.exists(), 'Do not overwrite a registered protocol'
    parent = json.loads(PARENT_PROTOCOL.read_text())
    dependencies = list(parent['script_hashes']) + ['stage43_partial_fpt.py', 'stage43_benchmark.py',
                                                    'stage43_solver_validation.py']
    for name, expected in parent['script_hashes'].items():
        assert sha(ROOT / 'scripts' / name) == expected
    protocol = dict(stage=43, registered_utc=datetime.now(timezone.utc).isoformat(),
        parent_protocol_sha256=sha(PARENT_PROTOCOL),
        description='Same frozen Stage42 exact minimum-source tasks; classical equality-reduction engineering',
        selection='All 30 Stage42 cases copied verbatim; optimization chosen from Stage42, no Stage43 timed tuning',
        per_job_hard_wall_seconds=5.0, per_job_address_space_bytes=1024**3,
        process_isolation='fresh Linux fork; all methods preloaded; sequential; wait4 CPU/RSS',
        methods=METHODS, method_order='rotate three methods by case_index mod 3; reverse when floor(case_index/3) is odd',
        thread_environment={key: os.environ[key] for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS')},
        analysis='All cases retained. Unknown is not infeasible. Pairwise minima only where both complete; no timeout-imputed speed ratios.',
        scope_caveats=parent['scope_caveats'] + [
            'Single timed trial per method/case; no confidence intervals or universal speed claims',
            'This follow-up suite is already known; improvement is not out-of-sample evidence',
            'Equality elimination and row reduction are classical LP tools, not a novelty claim'],
        case_list=parent['case_list'], script_hashes={n: sha(ROOT/'scripts'/n) for n in dependencies},
        pre_timing_validation_sha256=sha(ROOT/'results/stage43_solver_validation.json'),
        machine=dict(platform=platform.platform(), python=sys.version, processor=platform.processor(), cpu_count=os.cpu_count()))
    PROTOCOL.write_text(json.dumps(serial(protocol), indent=2)+'\n')
    print(json.dumps(dict(protocol=str(PROTOCOL), sha256=sha(PROTOCOL), cases=len(protocol['case_list'])), indent=2))


def run():
    assert not RESULT.exists(), 'Do not overwrite a measured run'
    protocol = json.loads(PROTOCOL.read_text())
    assert sha(PARENT_PROTOCOL) == protocol['parent_protocol_sha256']
    assert sha(ROOT/'results/stage43_solver_validation.json') == protocol['pre_timing_validation_sha256']
    for name, expected in protocol['script_hashes'].items():
        assert sha(ROOT/'scripts'/name) == expected, f'Preregistered code changed: {name}'
    records, started = [], time.perf_counter()
    for index, encoded in enumerate(protocol['case_list']):
        case = dict(encoded)
        case['edges'] = [(u, v, F(length)) for u, v, length in case['edges']]
        case['observed'] = {int(r): F(value) for r, value in case['observed'].items()}
        shift = index % 3
        methods = list(METHODS[shift:] + METHODS[:shift])
        if index // 3 % 2:
            methods.reverse()
        record = dict(name=case['name'], group=case['group'], cap=case['cap'],
                      vertices=case['vertices'], edges=len(case['edges']), original_leaves=case['original_leaves'],
                      observed_leaves=case['observed_leaves'], missing_fraction=case['missing_fraction'],
                      execution_order=methods, runs={})
        for method in methods:
            record['runs'][method] = isolated_run(method, case, protocol['per_job_hard_wall_seconds'],
                                                  protocol['per_job_address_space_bytes'])
        completed = [value for value in record['runs'].values() if value['status'] == 'complete']
        record['completed_methods_agree'] = len({(x['feasible'], x['minimum']) for x in completed}) <= 1
        assert record['completed_methods_agree'], (case['name'], record)
        for measured in completed:
            if 'expected_feasible' in case:
                assert measured['feasible'] == case['expected_feasible']
                if measured['feasible']:
                    assert measured['minimum'] == case['expected_minimum']
        records.append(record)
        print(f'{index+1:02d}/30 {case["name"]}: ' + ' | '.join(
            f'{m}={record["runs"][m]["status"]} min={record["runs"][m].get("minimum")} cpu={record["runs"][m]["child_cpu_seconds"]:.5f}s'
            for m in METHODS), flush=True)
    summary = dict(cases=len(records), disagreements=sum(not r['completed_methods_agree'] for r in records),
        pairwise_complete={f'{a} vs {b}':sum(r['runs'][a]['status']==r['runs'][b]['status']=='complete' for r in records)
            for i,a in enumerate(METHODS) for b in METHODS[i+1:]},
        methods={m:dict(status_counts=dict(Counter(r['runs'][m]['status'] for r in records)),
            verified_witnesses=sum(r['runs'][m].get('independent_witness_verified', False) for r in records),
            total_cpu_seconds=sum(r['runs'][m]['child_cpu_seconds'] for r in records),
            total_vertex_bases=sum(r['runs'][m].get('statistics',{}).get('vertex_bases',0) for r in records)) for m in METHODS},
        total_wall_seconds=time.perf_counter()-started)
    result = dict(stage=43, created_utc=datetime.now(timezone.utc).isoformat(), protocol_sha256=sha(PROTOCOL),
                  script_hashes=protocol['script_hashes'], records=records, summary=summary)
    RESULT.write_text(json.dumps(serial(result), indent=2)+'\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('phase', choices=['prepare', 'run'])
    args = parser.parse_args()
    prepare() if args.phase == 'prepare' else run()
