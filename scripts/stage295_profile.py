#!/usr/bin/env python3
"""Stage295 external profiler for the frozen Stage43 minimum pass.

The frozen modules are imported and NEVER modified.  All instrumentation is
external and uses only (a) the budget knobs `seconds`, `max_states`,
`max_bases`, `max_faces` that `solve_partial` already exposes, (b) the
`statistics` dict it already returns, and (c) `cProfile` around the call.

Three measurement modes:

  solve     one call with a given knob setting; records wall time, status and
            every counter at stop.  A geometric ladder over `max_faces` or
            `max_bases` turns this into a cumulative work/time curve whose
            slope separates per-face overhead from per-basis LP cost.
  cprofile  the same call under cProfile; records per-function tottime so the
            phase shares can be cross-checked against the ladder regression.

Each run executes in a forked child with RLIMIT_AS and a SIGALRM backstop, so
a run that exhausts memory or ignores its own budget cannot take the harness
down.  Every run is retained, including timeouts and failures.
"""
import json, os, resource, signal, sys, time
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
ROOT = os.path.dirname(HERE)

COUNTERS = ('deficit_sets', 'cover_feasible_faces', 'strict_faces',
            'growth_states', 'vertex_bases', 'mass_faces',
            'local_lp_calls', 'interval_lp_calls', 'enumerated_lp_calls',
            'zero_dimensional_lp_calls', 'eliminated_equality_rank',
            'constant_row_rejections', 'opposing_row_rejections',
            'inconsistent_equality_rejections')


from stage295_common import load_case as _load


def load_case(relpath):
    edges, obs, entry, raw = _load(relpath)
    return edges, obs


def summarize(answer):
    """JSON-safe projection of a solver answer: no Fraction leaks."""
    stats = answer.get('statistics', {}) or {}
    out = {'status': answer.get('status'), 'reason': answer.get('reason'),
           'feasible': bool(answer.get('feasible')),
           'minimum': answer.get('minimum_sources_within_cap'),
           'statistics': {k: stats.get(k, 0) for k in COUNTERS}}
    inc = answer.get('incumbent')
    if inc is not None:
        out['incumbent_source_count'] = inc.get('source_count')
    return out


class _Alarm(Exception):
    pass


def _child(row, path):
    knobs = dict(row.get('knobs') or {})
    payload = {'row': row}
    try:
        resource.setrlimit(resource.RLIMIT_AS,
                           (int(row.get('mem_gb', 6) * 2**30),) * 2)
    except (ValueError, OSError) as exc:
        payload['rlimit_error'] = str(exc)

    def _fire(signum, frame):
        raise _Alarm('wall_backstop')
    signal.signal(signal.SIGALRM, _fire)
    signal.alarm(int(row.get('wall', 900)))
    import stage43_partial_fpt as frozen
    edges, obs = load_case(row['case'])
    t0 = time.perf_counter()
    try:
        if row.get('mode') == 'cprofile':
            import cProfile, pstats
            pr = cProfile.Profile()
            pr.enable()
            answer = frozen.solve_partial(edges, obs, row['cap'], **knobs)
            pr.disable()
            payload['elapsed'] = time.perf_counter() - t0
            payload.update(summarize(answer))
            st = pstats.Stats(pr)
            funcs = []
            for func, entry in st.stats.items():
                cc, nc, tt, ct = entry[0], entry[1], entry[2], entry[3]
                funcs.append({'file': os.path.basename(func[0]), 'line': func[1],
                              'name': func[2], 'ncalls': nc, 'primcalls': cc,
                              'tottime': tt, 'cumtime': ct})
            funcs.sort(key=lambda r: -r['tottime'])
            payload['profile'] = funcs[:45]
            payload['profile_total_tt'] = st.total_tt
        else:
            answer = frozen.solve_partial(edges, obs, row['cap'], **knobs)
            payload['elapsed'] = time.perf_counter() - t0
            payload.update(summarize(answer))
        payload['ok'] = True
    except _Alarm as exc:
        payload.update(ok=False, error='wall_backstop', message=str(exc),
                       elapsed=time.perf_counter() - t0)
    except MemoryError as exc:
        payload.update(ok=False, error='MemoryError', message=str(exc)[:300],
                       elapsed=time.perf_counter() - t0)
    except BaseException as exc:
        payload.update(ok=False, error=type(exc).__name__,
                       message=str(exc)[:600], elapsed=time.perf_counter() - t0)
    signal.alarm(0)
    with open(path, 'w') as fh:
        json.dump(payload, fh, default=str)


def run_plan(plan, outdir, workers=8):
    os.makedirs(outdir, exist_ok=True)
    pending, live, done = list(enumerate(plan)), {}, []
    while pending or live:
        while pending and len(live) < workers:
            i, row = pending.pop(0)
            path = os.path.join(outdir, 'run_%04d.json' % i)
            pid = os.fork()
            if pid == 0:
                try:
                    _child(row, path)
                finally:
                    os._exit(0)
            live[pid] = (i, row, path, time.perf_counter())
        pid, status = os.wait()
        if pid not in live:
            continue
        i, row, path, t0 = live.pop(pid)
        wall = time.perf_counter() - t0
        try:
            payload = json.load(open(path))
        except Exception as exc:
            payload = {'row': row, 'ok': False, 'error': 'no_payload',
                       'message': '%s exit=%s' % (exc, status)}
        payload['harness_wall'] = wall
        payload['exit_status'] = status
        payload['index'] = i
        done.append(payload)
        print(json.dumps({k: payload.get(k) for k in
                          ('index', 'ok', 'status', 'reason', 'error', 'elapsed')} |
                         {'tag': row.get('tag'), 'cap': row.get('cap'),
                          'name': row.get('name'),
                          'mass_faces': payload.get('statistics', {}).get('mass_faces'),
                          'vertex_bases': payload.get('statistics', {}).get('vertex_bases')}),
              flush=True)
    done.sort(key=lambda r: r['index'])
    return done


if __name__ == '__main__':
    plan = json.load(open(sys.argv[1]))
    outdir = sys.argv[2]
    workers = int(sys.argv[3]) if len(sys.argv) > 3 else 8
    rows = run_plan(plan, outdir, workers)
    json.dump(rows, open(os.path.join(outdir, 'all_runs.json'), 'w'), indent=1)
    print('runs:', len(rows))
