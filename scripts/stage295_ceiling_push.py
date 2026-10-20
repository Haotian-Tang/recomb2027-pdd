#!/usr/bin/env python3
"""Stage295 ceiling push: run the reduced-face prototype at K >= 10 and certify
the answer without the frozen minimum pass.

For each case the child does two things:

  (1) solve at cap = K.  `solve_partial_reduced` inherits the frozen assertion
      `observe(edges, lifted['sources'], lifted['pendant'], observed) ==
      observed`, so any returned witness is exactly verified against the
      observations -- feasibility at k* is certified by construction.
  (2) if a minimum k* > 1 is returned, solve again at cap = k* - 1.  A complete
      infeasible answer there proves no explanation with fewer than k* sources
      exists, so k* is the exact minimum.  This certificate chain uses the
      prototype and the frozen witness checker only; it does not depend on the
      frozen minimum pass, which cannot finish at this K.
"""
import json, os, resource, signal, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
ROOT = os.path.dirname(HERE)
from stage295_common import load_case
from stage295_profile import summarize, COUNTERS


class _Alarm(Exception):
    pass


def _child(row, path):
    payload = {'row': row}
    try:
        resource.setrlimit(resource.RLIMIT_AS, (int(row.get('mem_gb', 10) * 2**30),) * 2)
    except (ValueError, OSError) as exc:
        payload['rlimit_error'] = str(exc)

    def _fire(signum, frame):
        raise _Alarm('wall_backstop')
    signal.signal(signal.SIGALRM, _fire)
    signal.alarm(int(row.get('wall', 4200)))
    try:
        import stage295_reduced_minimum as proto
        edges, obs, entry, raw = load_case(row['case'])
        payload['reference_set_size'] = len(obs)
        payload['generating_K'] = raw.get('cap')
        payload['n_generating_sources'] = len(raw.get('generating_sources') or [])
        knobs = dict(row.get('knobs') or {})
        opts = dict(full_mask=False, memo_cover=True,
                    popcount_order=row.get('popcount_order', False))
        t0 = time.perf_counter()
        ans = proto.solve_partial_reduced(edges, obs, row['cap'], **dict(knobs, **opts))
        primary = summarize(ans)
        st = ans.get('statistics', {}) or {}
        for extra in ('skipped_forced_faces', 'cover_memo_hits', 'cover_calls'):
            primary['statistics'][extra] = st.get(extra, 0)
        primary['elapsed'] = time.perf_counter() - t0
        primary['witness_exactly_verified'] = 'original_witness' in ans
        payload['primary'] = primary
        kstar = primary.get('minimum')
        if primary.get('status') == 'complete' and kstar and kstar > 1:
            t1 = time.perf_counter()
            low = proto.solve_partial_reduced(edges, obs, kstar - 1, **dict(knobs, **opts))
            payload['minimality_check'] = {
                'cap': kstar - 1, 'status': low.get('status'),
                'reason': low.get('reason'), 'feasible': bool(low.get('feasible')),
                'elapsed': time.perf_counter() - t1,
                'certifies_minimality': (low.get('status') == 'complete'
                                         and not low.get('feasible'))}
        payload['ok'] = True
    except _Alarm as exc:
        payload.update(ok=False, error='wall_backstop', message=str(exc))
    except BaseException as exc:
        payload.update(ok=False, error=type(exc).__name__, message=str(exc)[:600])
    signal.alarm(0)
    with open(path, 'w') as fh:
        json.dump(payload, fh, default=str)


def main():
    plan = json.load(open(sys.argv[1]))
    outdir = sys.argv[2]
    workers = int(sys.argv[3]) if len(sys.argv) > 3 else 6
    os.makedirs(outdir, exist_ok=True)
    pending, live, done = list(enumerate(plan)), {}, []
    while pending or live:
        while pending and len(live) < workers:
            i, row = pending.pop(0)
            p = os.path.join(outdir, 'push_%04d.json' % i)
            pid = os.fork()
            if pid == 0:
                try:
                    _child(row, p)
                finally:
                    os._exit(0)
            live[pid] = (i, row, p, time.perf_counter())
        pid, status = os.wait()
        if pid not in live:
            continue
        i, row, p, t0 = live.pop(pid)
        try:
            payload = json.load(open(p))
        except Exception as exc:
            payload = {'row': row, 'ok': False, 'error': 'no_payload',
                       'message': '%s exit=%s' % (exc, status)}
        payload['index'] = i
        payload['harness_wall'] = time.perf_counter() - t0
        done.append(payload)
        pr = payload.get('primary') or {}
        print(json.dumps({'index': i, 'name': row.get('name'), 'cap': row.get('cap'),
                          'pop': row.get('popcount_order'),
                          'status': pr.get('status'), 'min': pr.get('minimum'),
                          'secs': round(pr.get('elapsed', 0), 1),
                          'faces': (pr.get('statistics') or {}).get('mass_faces'),
                          'bases': (pr.get('statistics') or {}).get('vertex_bases'),
                          'minimality': (payload.get('minimality_check') or {}).get('certifies_minimality'),
                          'err': payload.get('error')}), flush=True)
    done.sort(key=lambda r: r['index'])
    json.dump(done, open(os.path.join(outdir, 'all_push.json'), 'w'), indent=1)
    ok = [r for r in done if (r.get('primary') or {}).get('status') == 'complete']
    print('push rows:', len(done), 'completed:', len(ok))


if __name__ == '__main__':
    main()
