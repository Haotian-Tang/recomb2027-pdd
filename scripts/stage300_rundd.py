#!/usr/bin/env python3
"""Stage 300: new DecoDiPhy runs on the real collections that have none yet.

Harvest first, run second.  Stages 289/293/294/297 already hold full-panel
DecoDiPhy answers for every krepp ANI case; the 60 converted krepp PLACEMENT
cases (stage 298) and the 29 travellers'-diarrhoea stool samples have none, so
they are executed here with the SAME worker the harvested runs used:
`stage293_abundance.dd_worker`, imported unmodified.  Only
`--method hill --optimizer osqp` is operable in DecoDiPhy 1.5.0, which is what
its documented defaults select.

Every run is retained, including hard-wall timeouts, memory kills and the
`save_jplace` KeyError the worker recovers from.  One forked child per job,
RLIMIT_AS inside, hard wall outside, OMP/OPENBLAS/MKL_NUM_THREADS=1.
"""
import os
for _k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_k] = '1'
import argparse
import json
import resource
import signal
import sys
import tempfile
import time
import traceback
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
sys.path.insert(0, str(ROOT / 'data' / 'stage289_decodiphy' / 'repo'))
import stage300_common as s300                   # noqa: E402
from stage293_abundance import dd_worker         # noqa: E402  frozen DD invocation


def build_jobs(collection, methods):
    jobs = []
    for entry in s300.case_entries(collection):
        case = s300.load_case(entry)
        k_true = len(case.get('generating_sources') or [])
        for method in methods:
            if method == 'dd_true_k' and not k_true:
                continue
            jobs.append(dict(collection=collection, name=entry['name'],
                             path=entry['path'], sha256=entry['sha256'],
                             family=entry['family'], method=method, k_true=k_true,
                             n_edges=len(case['edges']), n_leaves=len(case['observed'])))
    return jobs


def _child(spec, path, memory_bytes):
    try:
        resource.setrlimit(resource.RLIMIT_AS, (memory_bytes, memory_bytes))
        case = s300.load_case(spec)
        edges = [(int(u), int(v), F(L)) for u, v, L in case['edges']]
        observed = {int(r): F(v) for r, v in case['observed'].items()}
        sources = case.get('generating_sources') or []
        answer = dd_worker(dict(edges=edges, observed=observed,
                                generating_sources=sources,
                                method=spec['method'], k_true=spec['k_true']))
    except BaseException as exc:
        answer = dict(status='error', exception=repr(exc)[:500],
                      traceback=traceback.format_exc()[-2000:])
    try:
        Path(path).write_text(json.dumps(s300.serial(answer)))
    finally:
        os._exit(0)


def run_pool(jobs, out_path, workers, wall, memory_bytes, log_every=1):
    """Keep `workers` forked children alive; retain every outcome."""
    pending = list(jobs)
    running = {}
    done = 0
    handle = open(out_path, 'a')
    temp = tempfile.TemporaryDirectory(prefix='stage300_pool_')
    while pending or running:
        while pending and len(running) < workers:
            spec = pending.pop(0)
            path = Path(temp.name) / ('%d_%s.json' % (len(running) + done, os.getpid()))
            start = time.perf_counter()
            pid = os.fork()
            if pid == 0:
                handle.close()
                _child(spec, str(path), memory_bytes)
            running[pid] = dict(spec=spec, start=start, path=path)
        time.sleep(0.05)
        for pid in list(running):
            info = running[pid]
            waited, status, usage = os.wait4(pid, os.WNOHANG)
            timed_out = False
            if not waited:
                if time.perf_counter() - info['start'] < wall:
                    continue
                timed_out = True
                os.kill(pid, signal.SIGKILL)
                waited, status, usage = os.wait4(pid, 0)
            elapsed = time.perf_counter() - info['start']
            if timed_out:
                answer = dict(status='unknown', reason='hard_wall_timeout')
            elif Path(info['path']).exists():
                answer = json.loads(Path(info['path']).read_text())
            else:
                answer = dict(status='error', reason='child_wrote_no_result',
                              wait_status=status)
            row = dict(kind='dd_run', stage=300, **{k: v for k, v in info['spec'].items()})
            row.update(result=answer, harness_wall_seconds=elapsed,
                       wait_status=status, maximum_rss_kib=usage.ru_maxrss,
                       child_cpu_seconds=usage.ru_utime + usage.ru_stime,
                       timed_out=timed_out)
            handle.write(json.dumps(s300.serial(row)) + '\n')
            handle.flush()
            running.pop(pid)
            done += 1
            if done % log_every == 0:
                print('[%s] %d/%d done  last=%s/%s status=%s k=%s %.1fs rss=%.1fGB'
                      % (time.strftime('%H:%M:%S'), done, len(jobs),
                         info['spec']['name'], info['spec']['method'],
                         answer.get('status'), answer.get('k'), elapsed,
                         usage.ru_maxrss / 1048576), flush=True)
    handle.close()
    temp.cleanup()
    return done


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--collection', required=True)
    ap.add_argument('--methods', default='dd_default')
    ap.add_argument('--workers', type=int, default=8)
    ap.add_argument('--wall', type=float, default=1800.0)
    ap.add_argument('--memory-gb', type=float, default=8.0)
    ap.add_argument('--limit', type=int, default=0)
    ap.add_argument('--only', default='')
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    methods = tuple(m for m in args.methods.split(',') if m)
    jobs = build_jobs(args.collection, methods)
    if args.only:
        wanted = set(args.only.split(','))
        jobs = [j for j in jobs if j['name'] in wanted]
    if args.limit:
        jobs = jobs[:args.limit]
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    print('collection=%s methods=%s jobs=%d workers=%d wall=%.0fs mem=%.0fGB out=%s'
          % (args.collection, methods, len(jobs), args.workers, args.wall,
             args.memory_gb, out_path), flush=True)
    t0 = time.perf_counter()
    done = run_pool(jobs, out_path, args.workers, args.wall,
                    int(args.memory_gb * 1024 ** 3))
    print('finished %d jobs in %.1f s' % (done, time.perf_counter() - t0), flush=True)


if __name__ == '__main__':
    main()
