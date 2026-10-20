#!/usr/bin/env python3
"""Stage 289 exact-cost probe: how far beyond the matched 120 s budget the exact
solver needs on the source counts where it timed out.  One repetition, a longer
budget, and the solver's own internal limit so that its statistics are returned.
"""
import os
for _key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_key] = '1'
import json
import multiprocessing as mp
import resource
import signal
import sys
import tempfile
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from stage289_decodiphy_bridge import load_case
from stage43_partial_fpt import solve_partial
from stage289_headtohead import sha, serial

PROTOCOL = ROOT / 'results/stage289/stage289_headtohead_protocol.json'
OUT = ROOT / 'results/stage289/stage289_exact_cost_probe.json'
LIMIT = 1800
INTERNAL = 1750
MEMORY = 8 * 2 ** 30


def one(entry):
    case = load_case(ROOT / entry['path'])
    with tempfile.TemporaryDirectory(prefix='stage289_probe_') as temporary:
        path = Path(temporary) / 'result.json'
        start = time.perf_counter()
        pid = os.fork()
        if pid == 0:
            try:
                resource.setrlimit(resource.RLIMIT_AS, (MEMORY, MEMORY))
                began = time.perf_counter()
                raw = solve_partial(case['edges'], case['observed'], cap=entry['cap'],
                                    seconds=INTERNAL, max_states=None, max_bases=None,
                                    max_faces=None)
                answer = dict(status=raw['status'], feasible=raw.get('feasible'),
                              minimum=raw.get('minimum_sources_within_cap'),
                              reason=raw.get('reason'), statistics=raw.get('statistics'),
                              feasibility_established=raw.get('feasibility_established'),
                              solver_wall_seconds=time.perf_counter() - began)
            except BaseException as exc:
                answer = dict(status='error', exception=repr(exc)[:300],
                              traceback=traceback.format_exc()[-1500:])
            try:
                path.write_text(json.dumps(serial(answer)))
            finally:
                os._exit(0)
        timed_out = False
        while True:
            waited, status, usage = os.wait4(pid, os.WNOHANG)
            if waited:
                break
            if time.perf_counter() - start >= LIMIT:
                timed_out = True
                os.kill(pid, signal.SIGKILL)
                _, status, usage = os.wait4(pid, 0)
                break
            time.sleep(0.01)
        if timed_out:
            answer = dict(status='unknown', reason='hard_wall_timeout')
        elif path.exists():
            answer = json.loads(path.read_text())
        else:
            answer = dict(status='error', reason='worker_exited_without_result')
        answer.update(name=entry['name'], family=entry['family'], cap=entry['cap'],
                      k_true=entry['K_generating'], leaves=entry['reference_set_size'],
                      process_wall_seconds=time.perf_counter() - start,
                      child_cpu_seconds=usage.ru_utime + usage.ru_stime,
                      maximum_rss_kib=usage.ru_maxrss)
        return answer


def main():
    assert not OUT.exists(), 'Never overwrite measured runs'
    protocol = json.loads(PROTOCOL.read_text())
    index = protocol['case_index']
    chosen = []
    for family in ('exact_10k', 'exact_biotrees'):
        for k in (7, 10):
            group = sorted([e for e in index if e['family'] == family and e['K_generating'] == k],
                           key=lambda e: (e['reference_set_size'], e['name']))
            if not group:
                continue
            step = max(1, len(group) // 6)
            chosen.extend(group[::step][:6])
    with mp.get_context('fork').Pool(12) as pool:
        records = pool.map(one, chosen, chunksize=1)
    OUT.write_text(json.dumps(dict(
        stage=289, created_utc=datetime.now(timezone.utc).isoformat(),
        purpose='Cost of the exact solver beyond the 120 s matched budget.',
        protocol_sha256=sha(PROTOCOL), per_job_hard_wall_seconds=LIMIT,
        solver_internal_seconds=INTERNAL, repetitions=1, parallel_workers=12,
        cases=len(records), records=records), indent=1) + '\n')
    for r in records:
        print(r['family'], r['name'][:46], 'k', r['k_true'], 'leaves', r['leaves'],
              r['status'], r.get('reason'), 'min', r.get('minimum'),
              'cpu %.1f' % r['child_cpu_seconds'], flush=True)


if __name__ == '__main__':
    main()
