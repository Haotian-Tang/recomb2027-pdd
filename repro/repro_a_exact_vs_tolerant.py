#!/usr/bin/env python3
"""Reproduction A - the exact solver check (cheap).

PAPER CLAIM REPRODUCED
  sections/experiments.tex:  "At eps = 0 the tolerant program returns the frozen
  exact solver's minimum on all 1,118 leaf sets checked against it."
  appendices/experiments_details.tex: "At eps = 0 the tolerant program and the
  frozen exact minimum agree on 1,118 of 1,118 leaf sets, 0 disagreements."

WHAT IT DOES
  For a deterministic sample of the 1,118 registered leaf sets it re-runs BOTH
  registered code paths - the frozen structural solver
  (stage43_partial_fpt.solve_partial, via stage296_tolerance.frozen_worker) and
  the eps-tolerant rational program at eps = 0
  (stage296_tolerance.theta_worker, m(0) branch) - and checks

    1. the two agree with each other on every leaf set;
    2. each reproduces the value recorded in the frozen evidence streams
       results/stage296/stage296_validate_runs.jsonl  (frozen minimum)
       results/stage296/stage296_cross_runs.jsonl     (tolerant m(0)).

  Leaf sets are taken in increasing cost order (n_retained, then name, level,
  regime, replicate), so --n selects the cheapest n of the 1,118 and the
  selection is reproducible.  Dependencies: Python >= 3.12 standard library only.

USAGE
  python3 repro/repro_a_exact_vs_tolerant.py            # 120 leaf sets
  python3 repro/repro_a_exact_vs_tolerant.py --n 50
"""
import argparse
import json
import os
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))

import stage288_common as common            # noqa: E402
import stage296_tolerance as s296t          # noqa: E402

S296 = ROOT / 'results' / 'stage296'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, default=120,
                    help='number of leaf sets to re-run (cheapest first)')
    args = ap.parse_args()

    protocol = json.loads((S296 / 'stage296_protocol.json').read_text())
    by_name = {r['name']: r for r in protocol['cross_cases']}

    def stream(path):
        return [json.loads(l) for l in path.read_text().splitlines() if l.strip()]

    key = lambda r: (r['name'], r['level'], r['regime'], r['replicate'])
    frozen_rec = {key(r): r for r in stream(S296 / 'stage296_validate_runs.jsonl')}
    cross_rec = {key(r): r for r in stream(S296 / 'stage296_cross_runs.jsonl')}
    shared = sorted(set(frozen_rec) & set(cross_rec),
                    key=lambda k: (frozen_rec[k]['n_retained'],) + k)
    print('registered leaf sets: frozen %d, tolerant %d, joined %d'
          % (len(frozen_rec), len(cross_rec), len(shared)))
    sel = shared[:args.n]
    print('re-running the %d cheapest\n' % len(sel))

    hdr = ('%-44s %5s %3s %9s %9s %9s %9s %7s' %
           ('case', '|R|', 'K', 'frozen', 'rec.froz', 'm(0)', 'rec.m(0)', 'sec'))
    print(hdr); print('-' * len(hdr))

    n_agree = n_match_frozen = n_match_tolerant = 0
    disagreements = []
    t_start = time.perf_counter()
    for k in sel:
        row = frozen_rec[k]
        case = common.load_case(by_name[row['name']])
        keep = sorted(set(case['observed']) - set(row['removed']))
        t0 = time.perf_counter()
        frz = s296t.frozen_worker(dict(case=case, leaves=keep))
        tol = s296t.theta_worker(dict(case=case, leaves=keep,
                                      k_list=list(range(1, len(case['generating_sources']) + 1))))
        el = time.perf_counter() - t0
        a, b = frz.get('minimum'), tol.get('m_at_zero')
        rec_a = row['result'].get('minimum')
        rec_b = cross_rec[k]['result'].get('m_at_zero')
        n_agree += (a == b)
        n_match_frozen += (a == rec_a)
        n_match_tolerant += (b == rec_b)
        if not (a == b == rec_a == rec_b):
            disagreements.append(dict(key=k, frozen=a, tolerant=b,
                                      recorded_frozen=rec_a, recorded_tolerant=rec_b))
        print('%-44s %5d %3d %9s %9s %9s %9s %7.2f'
              % (row['name'][:44], row['n_retained'], row['K'], a, rec_a, b, rec_b, el))

    n = len(sel)
    print('\n' + '=' * 72)
    print('tolerant m(0) == frozen minimum          : %d of %d' % (n_agree, n))
    print('frozen minimum == recorded frozen stream : %d of %d' % (n_match_frozen, n))
    print('tolerant m(0) == recorded cross stream   : %d of %d' % (n_match_tolerant, n))
    print('disagreements                            : %d' % len(disagreements))
    print('wall seconds                             : %.1f' % (time.perf_counter() - t_start))
    if disagreements:
        print(json.dumps(disagreements, indent=1))
    ok = n_agree == n_match_frozen == n_match_tolerant == n and not disagreements
    print('RESULT: ' + ('PASS - reproduces the paper claim on this sample'
                        if ok else 'FAIL'))
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
