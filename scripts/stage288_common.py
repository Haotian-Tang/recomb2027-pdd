#!/usr/bin/env python3
"""Stage 288 shared library.

New file; nothing pre-existing is read-write. Provides frozen-bank access,
reference-removal regimes, an independent exact witness verifier, and a
forking job pool that keeps the stage 44 measurement conventions (fresh child
process per job, RLIMIT_AS, hard wall cap, every run retained) while allowing
bounded concurrency.
"""
import os
for _k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_k] = '1'
import hashlib
import json
import resource
import signal
import sys
import tempfile
import time
import traceback
from collections import deque
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / 'scripts') not in sys.path:
    sys.path.insert(0, str(ROOT / 'scripts'))
BANK = ROOT / 'data/stage287_instances/bank'


def serial(value):
    """Exactly the stage42_benchmark convention: Fractions become strings."""
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {str(k): serial(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [serial(x) for x in value]
    return value


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def bank_index():
    return json.loads((BANK / 'index.json').read_text())['cases']


def load_case(entry):
    """Read one frozen bank case, verifying its registered SHA256."""
    path = ROOT / entry['path']
    digest = sha(path)
    assert digest == entry['sha256'], f"bank case hash changed: {entry['name']}"
    raw = json.loads(path.read_text())
    case = dict(
        name=entry['name'], family=entry['family'], path=entry['path'], sha256=digest,
        edges=[(u, v, F(length)) for u, v, length in raw['edges']],
        observed={int(r): F(d) for r, d in raw['observed'].items()},
        cap=raw['cap'],
        generating_sources=[dict(edge=s['edge'], p=F(s['p']), x=F(s['x']))
                            for s in raw['generating_sources']],
        generating_pendant=F(raw['generating_pendant']),
        reference_set_size=entry['reference_set_size'],
        K_generating=entry['K_generating'])
    return case


# ---------------------------------------------------------------- tree helpers

def adjacency_lists(edges):
    adj = {}
    for u, v, length in edges:
        adj.setdefault(u, []).append((v, length))
        adj.setdefault(v, []).append((u, length))
    return adj


def leaves_of(edges):
    adj = adjacency_lists(edges)
    return sorted(v for v in adj if len(adj[v]) == 1)


def bfs_from(edges, roots):
    """Exact distances from each root to every vertex; independent of solvers."""
    adj = adjacency_lists(edges)
    out = {}
    for start in roots:
        dist = {start: F(0)}
        queue = deque([start])
        while queue:
            x = queue.popleft()
            for y, length in adj[x]:
                if y not in dist:
                    dist[y] = dist[x] + length
                    queue.append(y)
        out[start] = dist
    return out


def predict(edges, sources, pendant, targets):
    """d_r = y + sum_i p_i * min(d(r,u)+x_i L, d(r,v)+(1-x_i) L), exact."""
    endpoints = sorted({v for s in sources for v in edges[s['edge']][:2]})
    dist = bfs_from(edges, endpoints)
    out = {}
    for r in targets:
        total = F(pendant)
        for s in sources:
            u, v, length = edges[s['edge']]
            x = F(s['x'])
            total += F(s['p']) * min(dist[u][r] + x * length,
                                     dist[v][r] + (1 - x) * length)
        out[r] = total
    return out


def verify_witness(edges, observed, cap, witness):
    """Recompute the leaf equations on the original tree; trust no solver flag."""
    try:
        sources = [dict(edge=int(s['edge']), p=F(s['p']), x=F(s['x']))
                   for s in witness['sources']]
        y = F(witness['pendant'])
        checks = {
            'edge_indices_valid': all(0 <= s['edge'] < len(edges) for s in sources),
            'count_within_cap': 1 <= len(sources) <= cap,
            'distinct_edges': len({s['edge'] for s in sources}) == len(sources),
            'weights_sum_to_one': sum(s['p'] for s in sources) == 1,
            'weights_positive': all(s['p'] > 0 for s in sources),
            'strictly_interior': all(0 < s['x'] < 1 for s in sources),
            'pendant_positive': y > 0,
        }
        if not checks['edge_indices_valid']:
            return {'verified': False, 'checks': checks, 'reason': 'edge_index_out_of_range'}
        pred = predict(edges, sources, y, sorted(observed))
        deviation = max((abs(pred[r] - observed[r]) for r in observed), default=F(0))
        checks['reproduces_observations'] = deviation == 0
        return {'verified': all(checks.values()), 'checks': checks,
                'source_count': len(sources), 'max_abs_deviation': str(deviation)}
    except BaseException as exc:
        return {'verified': False, 'reason': repr(exc)[:300]}


# ------------------------------------------------------- removal regimes

def _rooted(edges):
    adj = {}
    for u, v, length in edges:
        adj.setdefault(u, []).append(v)
        adj.setdefault(v, []).append(u)
    root = next((v for v in sorted(adj) if len(adj[v]) > 1), sorted(adj)[0])
    parent = {root: None}
    order = [root]
    stack = [root]
    while stack:
        x = stack.pop()
        for y in sorted(adj[x]):
            if y not in parent:
                parent[y] = x
                order.append(y)
                stack.append(y)
    counts = {}
    for x in reversed(order):
        children = [y for y in adj[x] if parent.get(y) == x]
        counts[x] = 1 if not children and len(adj[x]) == 1 else sum(counts[y] for y in children)
    return adj, parent, order, counts, root


def _subtree_leaves(adj, parent, v):
    out, stack = [], [v]
    while stack:
        x = stack.pop()
        children = [y for y in adj[x] if parent.get(y) == x]
        if not children:
            out.append(x)
        else:
            stack.extend(children)
    return sorted(out)


def clade_removals(edges, leaves, target, replicates):
    """Whole-clade removals whose size is closest to `target`, deterministic."""
    adj, parent, order, counts, root = _rooted(edges)
    total = len(leaves)
    leafset = set(leaves)
    candidates = []
    for v in order:
        if v == root:
            continue
        size = counts[v]
        for complement in (False, True):
            drop = total - size if complement else size
            if 1 <= drop <= total - 3:
                candidates.append((abs(drop - target), drop, v, complement))
    candidates.sort()
    chosen, seen = [], set()
    for _, drop, v, complement in candidates:
        below = _subtree_leaves(adj, parent, v)
        removed = sorted(leafset - set(below)) if complement else below
        key = tuple(removed)
        if key in seen:
            continue
        seen.add(key)
        chosen.append({'removed': removed, 'clade_vertex': v, 'complement': complement})
        if len(chosen) == replicates:
            break
    return chosen


def nearest_source_removal(edges, sources, leaves, count):
    """Adversarial: drop the references closest to the true source points."""
    endpoints = sorted({v for s in sources for v in edges[s['edge']][:2]})
    dist = bfs_from(edges, endpoints)
    score = {}
    for r in leaves:
        best = None
        for s in sources:
            u, v, length = edges[s['edge']]
            x = F(s['x'])
            value = min(dist[u][r] + x * length, dist[v][r] + (1 - x) * length)
            best = value if best is None else min(best, value)
        score[r] = best
    return sorted(sorted(leaves, key=lambda r: (score[r], r))[:count])


def random_removal(leaves, count, seed):
    import random as _random
    rng = _random.Random(seed)
    return sorted(rng.sample(list(leaves), count))


# ------------------------------------------------------------- job pool

def _child(worker, spec, path, address_space_bytes):
    try:
        resource.setrlimit(resource.RLIMIT_AS, (address_space_bytes, address_space_bytes))
        answer = worker(spec)
    except BaseException as exc:
        answer = {'status': 'error', 'exception': repr(exc)[:400],
                  'traceback': traceback.format_exc()[-2000:]}
    try:
        path.write_text(json.dumps(serial(answer)))
    finally:
        os._exit(0)


def run_pool(jobs, worker, wall_seconds, address_space_bytes, workers,
             stream_path, post=None, progress=None):
    """Run `jobs` as forked children, at most `workers` at a time.

    jobs      - list of dicts with keys 'meta' (JSON-serialisable record) and
                'spec' (arbitrary live objects handed to `worker`).
    worker    - called in the child; returns the answer dict.
    post      - optional post(job, answer) run in the PARENT, outside child
                timing, for independent verification.
    Every job yields exactly one record, including timeouts and errors.
    """
    results = [None] * len(jobs)
    stream = open(stream_path, 'a')
    temporary = tempfile.TemporaryDirectory(prefix='stage288_pool_')
    pending = list(range(len(jobs)))
    live = {}
    done = 0
    started_at = time.perf_counter()
    try:
        while pending or live:
            while pending and len(live) < workers:
                index = pending.pop(0)
                path = Path(temporary.name) / f'{index}.json'
                sys.stdout.flush()
                sys.stderr.flush()
                stream.flush()
                begin = time.perf_counter()
                pid = os.fork()
                if pid == 0:
                    stream.close()
                    _child(worker, jobs[index]['spec'], path, address_space_bytes)
                live[pid] = (index, path, begin)
            finished = []
            for pid, (index, path, begin) in list(live.items()):
                waited, status, usage = os.wait4(pid, os.WNOHANG)
                timed_out = False
                if not waited:
                    if time.perf_counter() - begin < wall_seconds:
                        continue
                    timed_out = True
                    os.kill(pid, signal.SIGKILL)
                    _, status, usage = os.wait4(pid, 0)
                elapsed = time.perf_counter() - begin
                if timed_out:
                    answer = {'status': 'unknown', 'reason': 'hard_wall_timeout', 'minimum': None}
                elif path.exists():
                    answer = json.loads(path.read_text())
                else:
                    answer = {'status': 'error', 'reason': 'worker_exited_without_result'}
                answer.update(process_wall_seconds=elapsed,
                              child_cpu_seconds=usage.ru_utime + usage.ru_stime,
                              maximum_rss_kib=usage.ru_maxrss, wait_status=status)
                if post is not None:
                    try:
                        post(jobs[index], answer)
                    except BaseException as exc:
                        answer['post_error'] = repr(exc)[:300]
                record = dict(jobs[index]['meta'], result=serial(answer))
                results[index] = record
                stream.write(json.dumps(record) + '\n')
                finished.append(pid)
                done += 1
                if progress and done % progress == 0:
                    rate = done / max(1e-9, time.perf_counter() - started_at)
                    print(f'  {done}/{len(jobs)} done, {rate:.2f} jobs/s, '
                          f'{(len(jobs) - done) / max(rate, 1e-9) / 60:.1f} min left', flush=True)
            for pid in finished:
                live.pop(pid)
            if not finished:
                time.sleep(0.005)
    finally:
        stream.close()
        temporary.cleanup()
    return results


def dependency_closure(entry_name):
    """Transitive hashes of project scripts reachable by import from `entry_name`."""
    import ast
    pending, found = [entry_name], set()
    while pending:
        name = pending.pop()
        if name in found:
            continue
        path = ROOT / 'scripts' / name
        if not path.exists():
            continue
        found.add(name)
        for node in ast.walk(ast.parse(path.read_text())):
            names = ([node.module] if isinstance(node, ast.ImportFrom)
                     else [a.name for a in node.names] if isinstance(node, ast.Import) else [])
            for module in names:
                if module:
                    pending.append(module.split('.')[0] + '.py')
    return sorted(found)
