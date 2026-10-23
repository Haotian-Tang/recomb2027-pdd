#!/usr/bin/env python3
"""Stage305: do the unfinished threshold searches explain the regime gap?

A search that does not finish might be a hard one, so reporting medians over
completed searches alone could manufacture the gap between nearest-to-sources
and random deletion. This script bounds the median of the complete grid without
assuming anything about why a search stopped.

Along a fixed deletion order the generating support, once it stops being the
unique minimum-size support, cannot recover that property by deleting more
observations: deleting an observation deletes constraints, a competing
explanation stays feasible and a minimum that has already fallen cannot rise.
So every recorded probe constrains an unfinished search from one side --- the
largest deletion count still verified unique is a lower bound on its threshold,
the smallest failing count less one is an upper bound --- and a cell with no
probe at all is left free. The median is monotone in each entry, so filling
every free cell with its admissible extreme bounds the median of the complete
grid deterministically.

These are NOT confidence intervals. They say nothing about cases whose
generating support was not the unique minimum on the full panel, and nothing
about the rank correlations.

Run from the repository root:  python3 -B scripts/stage305_missing_thresholds.py
Writes results/stage305/stage305_missing_thresholds.json. No third-party imports.
"""
import glob
import hashlib
import json
import statistics
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results' / 'stage305' / 'stage305_missing_thresholds.json'

INPUTS = ['results/stage292_bridge_K2_3_5_7.jsonl',
          'results/stage292_bridge_K7_resume.jsonl']
GLOB = 'results/stage297/stage297_bridge_*.jsonl'

# one adversarial order per case, two seeds for each unstructured regime
SEEDS_PER_CASE = {'adversarial_nearest': 1, 'random': 2, 'whole_clade': 2}
UNBOUNDED = 10 ** 9


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    paths = [ROOT / p for p in INPUTS] + sorted(ROOT.glob(GLOB))
    missing = [str(p.relative_to(ROOT)) for p in paths if not p.exists()]
    if missing:
        print('missing inputs (fetch the bulk archive, see README):', *missing, sep='\n  ')
        return 1

    seen = set()
    cells = defaultdict(dict)
    status = defaultdict(lambda: defaultdict(int))
    baseline_cases = set()
    for p in paths:
        with p.open() as fh:
            for line in fh:
                r = json.loads(line)
                if r.get('kind') != 'bridge' or not r.get('baseline_unique_true'):
                    continue
                baseline_cases.add(r['name'])
                for g in r.get('regimes', []):
                    ident = (r['name'], g['regime'], g['seed'])
                    if ident in seen:
                        continue
                    seen.add(ident)
                    reg = g['regime']
                    status[reg][str(g.get('threshold_status')).split(':')[0]] += 1
                    t = g.get('t_unknown_sequence')
                    if t is not None:
                        cells[reg][ident] = (t, t)
                        continue
                    probes = g.get('probe_table') or []
                    good = [q['e'] for q in probes if q.get('ok') is True]
                    bad = [q['e'] for q in probes if q.get('ok') is False]
                    cells[reg][ident] = (max(good) if good else 0,
                                         (min(bad) - 1) if bad else UNBOUNDED)

    out = {
        'stage': 305,
        'what': 'deterministic bounds on the median deletion threshold when some searches did not finish',
        'monotonicity': ('along a fixed deletion order the property is lost once and for all, so a '
                         'completed probe bounds an unfinished search from one side'),
        'caveats': ['not confidence intervals',
                    'restricted to cases whose generating support is the unique minimum on the full panel',
                    'says nothing about the rank correlations'],
        'inputs': {str(p.relative_to(ROOT)): sha256(p) for p in paths},
        'baseline_unique_cases': len(baseline_cases),
        'by_regime': {},
    }

    for reg in sorted(cells):
        grid = len(baseline_cases) * SEEDS_PER_CASE.get(reg, 1)
        vals = list(cells[reg].values())
        exact = [lo for lo, hi in vals if lo == hi]
        pad = max(0, grid - len(vals))
        lows = [lo for lo, _ in vals] + [0] * pad
        highs = [hi for _, hi in vals] + [UNBOUNDED] * pad
        out['by_regime'][reg] = dict(
            balanced_grid=grid, completed=len(exact), recorded_cells=len(vals),
            observed_median=statistics.median(exact) if exact else None,
            bounded_median_low=statistics.median(lows),
            bounded_median_high=statistics.median(highs),
            threshold_status=dict(sorted(status[reg].items())))

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=1) + '\n')
    print('wrote', OUT.relative_to(ROOT))
    print('  %d baseline-unique cases' % out['baseline_unique_cases'])
    for reg, d in out['by_regime'].items():
        print('  %-20s completed %5d / grid %5d  observed median %-6s  bounded [%s, %s]'
              % (reg, d['completed'], d['balanced_grid'], d['observed_median'],
                 d['bounded_median_low'], d['bounded_median_high']))
    return 0


if __name__ == '__main__':
    sys.exit(main())
