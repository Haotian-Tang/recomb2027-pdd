#!/usr/bin/env python3
"""Stage 291 addendum: pendant-floor sensitivity of the certificate count.

The exclusion rule is  raw + eta - floor < alpha * offset, so with alpha*offset
small the supplied pendant floor is the dominant term. This sweeps the floor on
the 29 travellers' diarrhoea samples (no ground truth) and, as a labelled
control, on exact_10k and exact_biotrees where the true pendant y is known.
Registers results/stage291_floor_protocol.json before measuring.
"""
import json, sys, statistics, collections
from datetime import datetime, timezone
from fractions import Fraction as F
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from stage291_common import (bank_index, load_any, model_vector, screen_panel, sha, serial,
                             run_pool, dependency_closure)

PROTO = ROOT / 'results/stage291_floor_protocol.json'
STREAM = ROOT / 'results/stage291_floor_runs.jsonl'
RATIOS = (F(1, 100), F(1, 20), F(1, 10), F(1, 4), F(1, 2), F(3, 4), F(1), F(3, 2))
ETA = [F(0), F(1, 1000), F(1, 100), F(1, 10), F(1), F(10)]
FAMILIES = ('published_travelers_diarrhea', 'exact_10k', 'exact_biotrees')


def worker(spec):
    case = load_any(spec['entry'])
    keep = sorted(case['observed'])
    pan = screen_panel(case, keep)
    n = len(case['edges'])
    truth = model_vector(case)
    true_edges = sorted({s['edge'] for s in case['sources']}) if case['sources'] else []
    cap = case['cap'] or 10
    alpha, rho = F(1, 2 * cap), F(1, 10)
    reference = min(case['observed'].values())
    rows = []
    for ratio in RATIOS:
        floor = ratio * reference
        if floor <= 0:
            continue
        exact = {}
        for (raw_u, raw_v), items in zip(pan['raw'], pan['offsets']):
            for j, pre, suf, ell in items:
                best = None
                if raw_u is not None:
                    best = alpha * (pre + rho * ell) + floor - raw_u
                if raw_v is not None:
                    cand = alpha * (suf + rho * ell) + floor - raw_v
                    best = cand if best is None else max(best, cand)
                if best is not None:
                    exact[j] = best
        rows.append(dict(ratio=str(ratio), floor=float(floor),
                         floor_over_true_pendant=(float(floor / case['pendant'])
                                                  if case['pendant'] else None),
                         floor_admissible=(bool(case['pendant'] is not None
                                                and floor <= case['pendant'])),
                         excluded={str(e): sum(1 for v in exact.values() if v > e) for e in ETA},
                         true_excluded={str(e): sum(1 for j in true_edges
                                                    if j in exact and exact[j] > e) for e in ETA}))
    return dict(status='ok', family=case['family'], name=case['name'], n_edges=n,
                refs=len(keep), K=case['K_generating'], n_true=len(true_edges),
                median_observed=float(statistics.median(sorted(case['observed'].values()))),
                min_observed=float(reference),
                true_pendant=(float(case['pendant']) if case['pendant'] else None),
                rows=rows)


def main():
    if not PROTO.exists():
        PROTO.write_text(json.dumps(dict(
            stage=291, addendum='pendant_floor_sensitivity',
            created_utc=datetime.now(timezone.utc).isoformat(),
            parent_protocol='results/stage291_protocol_v2.json',
            parent_sha256=sha(ROOT / 'results/stage291_protocol_v2.json'),
            purpose=('Quantify how much of the certified exclusion count is produced by the '
                     'supplied pendant floor rather than by the data, on the real samples and '
                     'on labelled controls where admissibility of the floor is decidable.'),
            frozen_inputs={name: sha(ROOT / name) for name in
                           ('data/stage287_instances/bank/index.json',
                            'scripts/stage44_certificate.py', 'scripts/stage291_common.py',
                            'scripts/stage291_floor_sweep.py')},
            dependency_closure={n: sha(ROOT / 'scripts' / n)
                                for n in dependency_closure('stage291_floor_sweep.py')},
            families=list(FAMILIES), floor_ratios=[str(r) for r in RATIOS],
            floor_rule='floor = ratio * min_r observed_r', alpha='1/(2 cap)', rho='1/10',
            eta_grid=[str(e) for e in ETA],
            claim_boundary=('No ground truth exists for published_travelers_diarrhea. Only '
                            'certificate counts under explicitly supplied bounds are reported; '
                            'no recovery, accuracy or biological identification claim is made.')),
            indent=2) + '\n')
        print('registered', PROTO)
    assert not STREAM.exists(), 'refusing to append to an existing stream'
    entries = [e for e in bank_index() if e['family'] in FAMILIES]
    jobs = [dict(meta=dict(family=e['family'], name=e['name'], sha256=e['sha256']),
                 spec=dict(entry=e)) for e in entries]
    print('cases', len(jobs), flush=True)
    run_pool(jobs, worker, 3600.0, 4 * 1024 ** 3, 24, STREAM, progress=200)
    print('done ->', STREAM)


if __name__ == '__main__':
    main()
