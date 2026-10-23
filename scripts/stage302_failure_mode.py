#!/usr/bin/env python3
"""Stage302: what the breaking erasure actually does.

The erasure threshold recorded by `find_threshold` (scripts/stage290_bridge.py)
is the LAST deletion count at which the true support is still the unique
minimum-size support: the routine's own contract is `ok at t` and `not ok at
t+1`.  The success predicate requires all three of

    minimum == K,  n_supports == 1,  supports[0] == true support,

so the first failure at e = t + 1 need not be a tie at the same size: the
minimum can instead have fallen below K, in which case the true support is no
longer minimum at all.  This script classifies the recorded probe at e = t + 1
for every case that contributes a threshold to the pooled panel of Figure 2(c),
so the two events can be reported separately.  The same pass recomputes the
rank correlation of the threshold with panel size, which the paper reports as
the contrast between the adversarial and the random regime.

Inputs are the records stage297_analyze.py pools, filtered the same way
(`baseline_unique_true`) and de-duplicated the same way: one case can appear in
more than one batch, because a case that timed out in the five-regime batch is
re-measured in a regime-restricted one, so a (case, regime, seed) measurement is
counted once.  Repeats are counted and checked for agreement rather than pooled
twice.  De-duplicating by plotted coordinates would be wrong, since distinct
cases can share a [panel size, threshold] pair, and de-duplicating by case alone
would be wrong, since different seeds are planned distinct measurements.  `results/stage292_bridge_K2_3_5_7.jsonl` ships in the
bulk archive rather than the repository tree; the script reports which inputs
it found and refuses to write a record if any are missing.

Run from the repository root:  python3 -B scripts/stage302_failure_mode.py
Writes results/stage302/stage302_failure_mode.json.  No third-party imports.
"""
import hashlib
import json
import math
import statistics
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results' / 'stage302' / 'stage302_failure_mode.json'

INPUTS = ['results/stage292_bridge_K2_3_5_7.jsonl',
          'results/stage292_bridge_K7_resume.jsonl']
GLOB = 'results/stage297/stage297_bridge_*.jsonl'


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def spearman(xs, ys):
    """Rank correlation with midranks for ties, and the two-sided t statistic."""
    n = len(xs)
    if n < 3:
        return None, None

    def ranks(v):
        order = sorted(range(n), key=lambda i: v[i])
        out = [0.0] * n
        i = 0
        while i < n:
            j = i
            while j + 1 < n and v[order[j + 1]] == v[order[i]]:
                j += 1
            for t in range(i, j + 1):
                out[order[t]] = (i + j) / 2.0 + 1.0
            i = j + 1
        return out

    a, b = ranks(xs), ranks(ys)
    ma, mb = sum(a) / n, sum(b) / n
    num = sum((p - ma) * (q - mb) for p, q in zip(a, b))
    den = math.sqrt(sum((p - ma) ** 2 for p in a) * sum((q - mb) ** 2 for q in b))
    if den == 0:
        return None, None
    rho = num / den
    if abs(rho) >= 1.0:
        return rho, None
    return rho, rho * math.sqrt((n - 2) / (1 - rho * rho))


def main():
    paths = [ROOT / p for p in INPUTS] + sorted(ROOT.glob(GLOB))
    missing = [str(p.relative_to(ROOT)) for p in paths if not p.exists()]
    if missing:
        print('missing inputs (fetch the bulk archive, see README):', *missing, sep='\n  ')
        return 1

    rows = []
    for p in paths:
        with p.open() as fh:
            for line in fh:
                r = json.loads(line)
                if r.get('kind') == 'bridge' and r.get('baseline_unique_true'):
                    rows.append(r)

    mode = defaultdict(lambda: defaultdict(int))
    thresholds = defaultdict(list)
    panel = defaultdict(list)
    cases = defaultdict(set)
    seen, repeats, repeats_agree = {}, 0, 0
    for r in rows:
        K = r['K']
        for g in r.get('regimes', []):
            ident = (r['name'], g['regime'], g['seed'])
            t = g.get('t_unknown_sequence')
            if ident in seen:
                repeats += 1
                repeats_agree += int(seen[ident] == t)
                continue
            seen[ident] = t
            if t is None:
                continue
            reg = g['regime']
            cases[reg].add(r['name'])
            thresholds[reg].append(t)
            panel[reg].append((r['R_size'], t))
            probe = {q['e']: q for q in g.get('probe_table', [])}.get(t + 1)
            if probe is None:
                mode[reg]['no_probe_recorded'] += 1
                continue
            minimum, n_supports = probe.get('minimum'), probe.get('n_supports')
            if minimum is None:
                mode[reg]['probe_incomplete'] += 1
            elif minimum < K:
                mode[reg]['minimum_fell_below_K'] += 1
            elif minimum == K and (n_supports or 0) > 1:
                mode[reg]['competing_support_same_size'] += 1
            else:
                mode[reg]['unique_but_not_true_support'] += 1

    out = {
        'stage': 302,
        'what': 'classification of the first erasure that breaks uniqueness, at e = t_unknown + 1',
        'threshold_convention': ('t_unknown is the last deletion count at which the true support is '
                                 'still the unique minimum; the first failure is at t_unknown + 1'),
        'filter': 'baseline_unique_true, same as stage297_analyze.py pooled',
        'inputs': {str(p.relative_to(ROOT)): sha256(p) for p in paths},
        'bridge_records_used': len(rows),
        'measurements_after_dedup': len(seen),
        'repeated_measurements_skipped': repeats,
        'repeated_measurements_agreeing': repeats_agree,
        'by_regime': {},
    }
    for reg in sorted(thresholds):
        v = thresholds[reg]
        sizes = [a for a, _ in panel[reg]]
        rho, tstat = spearman(sizes, [b for _, b in panel[reg]])
        out['by_regime'][reg] = {
            'n_thresholds': len(v),
            'distinct_cases': len(cases[reg]),
            'median_t_unknown': statistics.median(v),
            'min_t_unknown': min(v),
            'max_t_unknown': max(v),
            'panel_size_min': min(sizes),
            'panel_size_max': max(sizes),
            'spearman_t_vs_panel_size': rho,
            'spearman_t_statistic': tstat,
            'first_failure_mode': dict(sorted(mode[reg].items())),
        }

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=1) + '\n')
    print('wrote', OUT.relative_to(ROOT))
    print(f"  {repeats} repeated (case, regime, seed) measurements skipped, "
          f"{repeats_agree} of them agreeing with the first reading")
    for reg, d in out['by_regime'].items():
        print(f"  {reg:22s} n={d['n_thresholds']:5d}  median={d['median_t_unknown']}"
              f"  cases={d['distinct_cases']:5d}  panel {d['panel_size_min']}-{d['panel_size_max']}"
              f"  rho={d['spearman_t_vs_panel_size']:.4f}  {d['first_failure_mode']}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
