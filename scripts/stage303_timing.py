#!/usr/bin/env python3
"""Stage303: the paired timing comparison against DecoDiPhy 1.5.0.

The main text quotes a median speed ratio. This script is the record behind it.
It reads the frozen head-to-head runs and, for every case where both arms
completed, forms the per-case ratio of DecoDiPhy's time to the exact method's
time, taking the median over repetitions within each arm first.

Three timing fields are reported because they differ by more than a factor of
one and a half, and two DecoDiPhy arms because the tool searching for its own k
does strictly more work than the exact method solving a fixed K:

  process_wall_seconds  end-to-end, what a user waits
  solver_wall_seconds   the solve alone
  child_cpu_seconds     CPU consumed by the forked child

  dd_true_k             the oracle supplies k; the like-for-like comparison
  dd_default            the tool chooses its own k

The exact arm does not complete on every run. Cases where it timed out are
absent from the ratio, so every number here is favourable to the exact method
and is an upper bound on what a user would see. Hardware was not recorded for
this batch (Appendix H), so these are ratios on one machine, not absolute times.

Run from the repository root:  python3 -B scripts/stage303_timing.py
Writes results/stage303/stage303_timing.json. No third-party imports.
"""
import hashlib
import json
import statistics
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUNS = ROOT / 'results' / 'stage289' / 'stage289_headtohead_runs.jsonl'
OUT = ROOT / 'results' / 'stage303' / 'stage303_timing.json'

EXACT = 'exact_stage43_partial'
FIELDS = ['process_wall_seconds', 'solver_wall_seconds', 'child_cpu_seconds']
ARMS = ['dd_true_k', 'dd_default', 'dd_resolution']
DONE = {'ok', 'complete', 'completed', 'success'}


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def describe(v):
    v = sorted(v)
    q = statistics.quantiles(v, n=4) if len(v) >= 4 else [None, None, None]
    return dict(n=len(v), median=statistics.median(v), q25=q[0], q75=q[2],
                min=v[0], max=v[-1])


def main():
    if not RUNS.exists():
        print('missing input (fetch the bulk archive, see README):',
              RUNS.relative_to(ROOT))
        return 1
    rows = [json.loads(line) for line in RUNS.open()]

    status = defaultdict(lambda: defaultdict(int))
    for r in rows:
        status[r['method']][str(r['result'].get('status'))] += 1

    out = {
        'stage': 303,
        'what': 'paired per-case wall-time ratio, DecoDiPhy 1.5.0 over the exact method',
        'pairing': ('per case name, median over repetitions within each arm, then the '
                    'ratio; only cases where both arms completed'),
        'caveats': [
            'cases where the exact arm did not complete are absent, so every ratio '
            'is favourable to the exact method',
            'hardware was not recorded for this batch, so these are ratios on one '
            'machine and not absolute times',
            'DecoDiPhy was run outside the per-read pipeline it was designed for, '
            'on one seed family and only its operable code paths',
        ],
        'input': {str(RUNS.relative_to(ROOT)): sha256(RUNS)},
        'runs': len(rows),
        'distinct_cases': len({r['name'] for r in rows}),
        'status_by_method': {m: dict(sorted(d.items())) for m, d in sorted(status.items())},
        'ratios': {},
    }

    for field in FIELDS:
        per = defaultdict(lambda: defaultdict(list))
        for r in rows:
            if str(r['result'].get('status')) not in DONE:
                continue
            v = r['result'].get(field)
            if v is None or v <= 0:
                continue
            per[r['name']][r['method']].append(v)
        for arm in ARMS:
            ratios, faster = [], 0
            for name, m in per.items():
                if EXACT not in m or arm not in m:
                    continue
                e = statistics.median(m[EXACT])
                d = statistics.median(m[arm])
                if e <= 0:
                    continue
                ratios.append(d / e)
                faster += int(d > e)
            if ratios:
                rec = describe(ratios)
                rec['exact_faster_in'] = faster
                out['ratios']['%s|%s' % (field, arm)] = rec

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=1) + '\n')
    print('wrote', OUT.relative_to(ROOT))
    for k, v in out['ratios'].items():
        print('  %-38s n=%4d  median=%7.2f  q25=%6.2f  q75=%7.2f  exact faster %d'
              % (k, v['n'], v['median'], v['q25'], v['q75'], v['exact_faster_in']))
    return 0


if __name__ == '__main__':
    sys.exit(main())
