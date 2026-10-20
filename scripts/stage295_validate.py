#!/usr/bin/env python3
"""Stage295 validation and head-to-head: frozen minimum pass vs the Stage295
reduced-face prototype.

Per case the child runs the frozen `stage43_partial_fpt.solve_partial` and
each requested prototype variant on identical inputs, and records the minimum
each returns together with its counters and wall time.  The agreement rule is
fixed in the protocol: on every case where both complete, the minimum source
count must be identical.  A disagreement is reported, never patched.
"""
import json, os, resource, signal, sys, time
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
ROOT = os.path.dirname(HERE)
from stage295_common import load_case
from stage295_profile import COUNTERS, summarize


class _Alarm(Exception):
    pass


VARIANTS = {
    'frozen':            None,
    'reduced':           dict(full_mask=False, memo_cover=False, popcount_order=False),
    'reduced_memo':      dict(full_mask=False, memo_cover=True,  popcount_order=False),
    'reduced_memo_pop':  dict(full_mask=False, memo_cover=True,  popcount_order=True),
    'full_mask_memo':    dict(full_mask=True,  memo_cover=True,  popcount_order=False),
}


def _one(solver, edges, obs, cap, knobs, variant):
    import stage43_partial_fpt as frozen
    import stage295_reduced_minimum as proto
    t0 = time.perf_counter()
    try:
        if solver == 'frozen':
            ans = frozen.solve_partial(edges, obs, cap, **knobs)
        else:
            ans = proto.solve_partial_reduced(edges, obs, cap, **dict(knobs, **VARIANTS[variant]))
        out = summarize(ans)
        st = ans.get('statistics', {}) or {}
        for extra in ('skipped_forced_faces', 'cover_memo_hits', 'cover_calls'):
            out['statistics'][extra] = st.get(extra, 0)
        out['elapsed'] = time.perf_counter() - t0
        out['ok'] = True
    except BaseException as exc:
        out = {'ok': False, 'error': type(exc).__name__, 'message': str(exc)[:400],
               'elapsed': time.perf_counter() - t0}
    return out


def _child(row, path):
    payload = {'row': row}
    try:
        resource.setrlimit(resource.RLIMIT_AS, (int(row.get('mem_gb', 6) * 2**30),) * 2)
    except (ValueError, OSError) as exc:
        payload['rlimit_error'] = str(exc)

    def _fire(signum, frame):
        raise _Alarm('wall_backstop')
    signal.signal(signal.SIGALRM, _fire)
    signal.alarm(int(row.get('wall', 1800)))
    try:
        edges, obs, entry, raw = load_case(row['case'])
        payload['reference_set_size'] = len(obs)
        payload['generating_K'] = raw.get('cap')
        knobs = dict(row.get('knobs') or {})
        results = {}
        for variant in row['variants']:
            results[variant] = _one('frozen' if variant == 'frozen' else 'proto',
                                    edges, obs, row['cap'], knobs, variant)
        payload['results'] = results
        mins = {v: r.get('minimum') for v, r in results.items()
                if r.get('ok') and r.get('status') == 'complete'}
        payload['minima'] = mins
        payload['agree'] = (len(set(mins.values())) <= 1)
        payload['ok'] = True
    except _Alarm as exc:
        payload.update(ok=False, error='wall_backstop', message=str(exc))
    except BaseException as exc:
        payload.update(ok=False, error=type(exc).__name__, message=str(exc)[:600])
    signal.alarm(0)
    with open(path, 'w') as fh:
        json.dump(payload, fh, default=str)


def run_plan(plan, outdir, workers=8):
    os.makedirs(outdir, exist_ok=True)
    pending, live, done = list(enumerate(plan)), {}, []
    while pending or live:
        while pending and len(live) < workers:
            i, row = pending.pop(0)
            path = os.path.join(outdir, 'val_%04d.json' % i)
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
        try:
            payload = json.load(open(path))
        except Exception as exc:
            payload = {'row': row, 'ok': False, 'error': 'no_payload',
                       'message': '%s exit=%s' % (exc, status)}
        payload['harness_wall'] = time.perf_counter() - t0
        payload['index'] = i
        done.append(payload)
        print(json.dumps({'index': i, 'cap': row.get('cap'), 'name': row.get('name'),
                          'agree': payload.get('agree'), 'minima': payload.get('minima'),
                          'err': payload.get('error'),
                          'secs': {v: round(r.get('elapsed', 0), 2)
                                   for v, r in (payload.get('results') or {}).items()}}),
              flush=True)
    done.sort(key=lambda r: r['index'])
    return done


if __name__ == '__main__':
    plan = json.load(open(sys.argv[1]))
    outdir = sys.argv[2]
    workers = int(sys.argv[3]) if len(sys.argv) > 3 else 8
    rows = run_plan(plan, outdir, workers)
    json.dump(rows, open(os.path.join(outdir, 'all_validation.json'), 'w'), indent=1)
    both = [r for r in rows if r.get('ok') and len(r.get('minima') or {}) > 1]
    bad = [r for r in both if not r.get('agree')]
    print('validation rows:', len(rows), 'comparable:', len(both), 'DISAGREEMENTS:', len(bad))
    for r in bad:
        print('DISAGREE', r['row']['name'], r['row']['cap'], r['minima'])
