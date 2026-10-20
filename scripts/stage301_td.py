#!/usr/bin/env python3
"""Stage 301 step 2: certified lower bounds on Theta_k for real data.

The frozen count obstruction, transported, gives a support-free certificate

    Theta_k >= eps_R(k) = max( (4k-7)-th largest critical eps of NEG(.),
                               (3k-1)-th largest critical eps of POS(.) ),

computed from one envelope pass.  Combined with the exactly known Theta_1 of
stage 298 this brackets every Theta_k and therefore bounds m(eps) from below:
for eps <= eps_R(k) the data do not admit k sources.

Arms
  td        29 travellers'-diarrhoea stool samples (no ground truth).
  krepp     real ANI and real placement cases, where stage 298 already knows
            Theta_1..Theta_3 exactly, so eps_R(k) <= Theta_k is a SOUNDNESS
            CHECK on real estimated distances and the ratio measures tightness.
"""
import json
import os
import sys
import time
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import stage288_common as common            # noqa: E402
import stage298_common as s298              # noqa: E402
import stage301_common as s301              # noqa: E402

OUT = s301.OUT
KMAX = 8
WORKERS = 40
ADDRESS = 10 * 1024 ** 3


def worker(spec):
    t0 = time.perf_counter()
    case = json.loads(Path(spec['path']).read_text())
    edges = [(u, v, F(l)) for u, v, l in case['edges']]
    obs = {int(r): F(d) for r, d in case['observed'].items()}
    view = s301.backbone_view(edges, obs)
    t_view = time.perf_counter() - t0
    t1 = time.perf_counter()
    pos, neg = s301.critical_eps(view)
    bounds = {}
    for k in range(1, KMAX + 1):
        eps_r, meta = s301.refutation_threshold(view, k)
        bounds[k] = dict(eps_refute=eps_r, **{a: b for a, b in meta.items()})
    pv = sorted(pos.values(), reverse=True)
    return dict(status='complete', n_edges=len(edges), n_leaves=len(obs),
                n_backbone_edges=len(view['backbone']),
                Y=view['Y'], sum_b=sum(view['b'].values()),
                max_leaf_envelope_slack=max(view['envelope_slack'].values()),
                n_pos_total=len(pos), n_neg_total=len(neg),
                eps_vacuous_pos=(pv[0] if pv else F(0)),
                bounds=bounds, seconds_view=t_view,
                seconds_bounds=time.perf_counter() - t1)


def main():
    assert (OUT / 'stage301_protocol.json').exists()
    proto = json.loads((OUT / 'stage301_protocol.json').read_text())
    jobs = []
    for row in proto['td_cases']:
        path = ROOT / row['path']
        assert common.sha(path) == row['sha256'], 'input hash mismatch'
        jobs.append(dict(meta=dict(arm='td', name=Path(row['path']).stem,
                                   path=row['path']),
                         spec=dict(path=str(path))))
    index = common.bank_index()
    sweeps = [json.loads(l) for l in
              (ROOT / 'results/stage298/stage298_sweep_runs.jsonl').read_text().splitlines()]
    known = {}
    for r in sweeps:
        if r['result'].get('theta') is not None:
            known.setdefault(r['name'], {})[r['k']] = r['result']['theta']
    manifest = json.loads((ROOT / 'results/stage298/stage298_placement_manifest.json').read_text())
    pl = {}
    for item in (manifest if isinstance(manifest, list) else manifest.get('cases', [])):
        pl[item['name']] = item['path']
    seen = set()
    for name, th in sorted(known.items()):
        if name in seen:
            continue
        seen.add(name)
        entry = next((e for e in index if e['name'] == name), None)
        path = ROOT / entry['path'] if entry else (ROOT / pl[name] if name in pl else None)
        if path is None or not Path(path).exists():
            continue
        jobs.append(dict(meta=dict(arm='krepp', name=name, path=str(path),
                                   known_thetas=th),
                         spec=dict(path=str(path))))
    print('jobs', len(jobs), flush=True)
    stream = OUT / 'stage301_td_runs.jsonl'
    common.run_pool(jobs, worker, 1800, ADDRESS, WORKERS, str(stream))
    print('done', flush=True)


if __name__ == '__main__':
    main()
