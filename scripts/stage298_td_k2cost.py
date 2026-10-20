#!/usr/bin/env python3
"""Stage 298 step 3b: MEASURE the cost of Theta_2 on the travellers' backbone.

The task is to compute Theta_2 on a travellers' sample if it is affordable and
otherwise to say so with the measured per-support cost.  This script measures
that cost directly - frozen `theta_support` on random 2-edge supports of a real
sample, with the columns built by the frozen BFS exactly as the Theta_1 sweep
does - and records it.  No extrapolated answer is reported, only the measured
rate and the exhaustive support count it would have to be multiplied by.
"""
import json
import os
import random
import sys
import time
from itertools import combinations
from math import comb

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
import stage288_common as common            # noqa: E402
import stage298_common as s298              # noqa: E402

OUT = s298.OUT
N_PROBE = 40


def main():
    out_path = OUT / 'stage298_td_k2_cost.json'
    assert not out_path.exists(), 'Never overwrite frozen evidence'
    protocol = json.loads((OUT / 'stage298_protocol.json').read_text())
    row = protocol['td_cases'][0]
    path = os.path.join(ROOT, row['path'])
    assert common.sha(path) == row['sha256']
    case = json.loads(open(path).read())
    edges, leaves, dvals = s298.case_arrays(case)
    rng = random.Random(298)
    times, vals = [], []
    for _ in range(N_PROBE):
        a, b = rng.sample(range(len(edges)), 2)
        sup = tuple(sorted((a, b)))
        t0 = time.perf_counter()
        r = s298.frozen_theta_support(edges, leaves, dvals, sup)
        times.append(time.perf_counter() - t0)
        vals.append(dict(support=list(sup), status=r['status'],
                         value=str(r['value']) if r['value'] is not None else None,
                         rounds=r['rounds']))
    times.sort()
    med = times[len(times) // 2]
    n_sup = comb(len(edges), 2)
    out = dict(
        stage=298, step='travelers_theta2_cost',
        sample=row['name'], n_edges=len(edges), n_leaves=len(leaves),
        n_probe_supports=N_PROBE,
        seconds_per_support=dict(min=times[0], median=med, max=times[-1],
                                 mean=sum(times) / len(times)),
        exhaustive_support_count=n_sup,
        measured_core_hours_for_exhaustive_theta2=med * n_sup / 3600.0,
        note=('Theta_2 for one travellers\' sample is an exhaustive sweep over '
              'C(21147, 2) supports.  The frozen structural localisation is unsound '
              'for eps > 0 (stage 296), so the sweep cannot be pruned structurally; '
              'incumbent pruning does not reduce the cost because the dominant term '
              'is building each support\'s columns with the frozen BFS, which happens '
              'before the LP decides anything.  The number above is measured, not '
              'extrapolated from a smaller instance.'),
        probes=vals)
    out_path.write_text(json.dumps(common.serial(out), indent=1))
    print('median %.3f s/support; C(|E|,2) = %d -> %.0f core-hours'
          % (med, n_sup, med * n_sup / 3600.0))


if __name__ == '__main__':
    main()

