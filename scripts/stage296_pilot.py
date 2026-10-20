#!/usr/bin/env python3
"""Stage 296 pilot: time the exhaustive eps-tolerant sweeps so the registered
panel can be sized.  Measurement only; no protocol, no retained verdicts."""
import sys, time, json
from fractions import Fraction as F
sys.path.insert(0, 'scripts')
import stage288_common as common, stage296_common as s296

idx = common.bank_index()
by = {}
for e in idx:
    if e['family'].startswith('exact'):
        by.setdefault((e['family'], e['K_generating']), []).append(e)
for v in by.values():
    v.sort(key=lambda e: (e['reference_set_size'], e['name']))

def pick(family, K, target):
    v = by[(family, K)]
    return min(v, key=lambda e: (abs(e['reference_set_size'] - target), e['name']))

rows = []
for fam, K, target, ks in [('exact_10k', 3, 27, (1, 2)), ('exact_10k', 3, 40, (2,)),
                           ('exact_10k', 3, 60, (2,)), ('exact_10k', 3, 80, (2,)),
                           ('exact_10k', 3, 110, (2,)), ('exact_biotrees', 3, 101, (2,)),
                           ('exact_biotrees', 2, 101, (1,)), ('exact_biotrees', 2, 361, (1,)),
                           ('exact_10k', 5, 27, (2, 3)), ('exact_10k', 10, 27, (2,))]:
    entry = pick(fam, K, target)
    case = common.load_case(entry)
    leaves = sorted(case['observed']); dvals = [case['observed'][r] for r in leaves]
    t0 = time.perf_counter(); Dv = s296.leaf_vectors(case['edges'], leaves); tprep = time.perf_counter() - t0
    row = dict(name=entry['name'], family=fam, K=K, N=len(leaves), E=len(case['edges']),
               prep_seconds=round(tprep, 2))
    for k in ks:
        t0 = time.perf_counter()
        out = s296.theta_min_over_supports(Dv, case['edges'], dvals, k)
        t = time.perf_counter() - t0
        row[f'k{k}_seconds'] = round(t, 1)
        row[f'k{k}_theta'] = float(out['theta']) if out.get('theta') is not None else None
        row[f'k{k}_visited'] = out.get('visited'); row[f'k{k}_pruned'] = out.get('pruned')
        row[f'k{k}_ms_per_support'] = round(1000 * t / max(1, out.get('visited', 1)), 2)
    rows.append(row); print(json.dumps(row), flush=True)

# k = 4 reach probe on the smallest K=5 case, with a wall budget
entry = pick('exact_10k', 5, 20)
case = common.load_case(entry)
leaves = sorted(case['observed']); dvals = [case['observed'][r] for r in leaves]
Dv = s296.leaf_vectors(case['edges'], leaves)
import itertools
E = len(case['edges'])
t0 = time.perf_counter(); best = None; visited = 0; hot = []
for support in itertools.combinations(range(E), 4):
    visited += 1
    out = s296.theta_support(s296.support_cols(Dv, case['edges'], support), dvals, 4,
                             cutoff=best, seed=hot[:8])
    if out['status'] == 'optimal':
        hot = list(dict.fromkeys(out['active'] + hot))[:8]
        if best is None or out['value'] < best:
            best = out['value']
    if time.perf_counter() - t0 > 900:
        break
el = time.perf_counter() - t0
from math import comb
print(json.dumps(dict(k4_probe=entry['name'], N=len(leaves), E=E, total=comb(E, 4),
                      visited=visited, seconds=round(el, 1),
                      ms_per_support=round(1000 * el / visited, 2),
                      projected_hours=round(comb(E, 4) * el / visited / 3600, 2))), flush=True)
