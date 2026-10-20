#!/usr/bin/env python3
"""Stage 300 protocol registration: inputs, code and plan, hashed before the runs.

House pattern: every input case, every imported stage-298 result, every harvest
source and every code file this stage executes is hashed here, before the timed
work, so that any later change is visible.  Two run groups were launched before
this file was written - the 60 krepp-placement cases and the first
travellers'-diarrhoea pilot sample - from the identical `stage300_rundd.py` and
`stage300_common.py` whose hashes are registered below; that is recorded in
`launched_before_registration` rather than hidden.
"""
import json
import platform
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import stage300_common as s300                   # noqa: E402

VENV = ROOT / 'data' / 'stage289_decodiphy' / 'venv' / 'bin' / 'python'
REPO = ROOT / 'data' / 'stage289_decodiphy' / 'repo' / 'decodiphy'
CODE = ['scripts/stage300_common.py', 'scripts/stage300_rundd.py',
        'scripts/stage300_harvest.py', 'scripts/stage300_evaluate.py',
        'scripts/stage300_verify.py', 'scripts/stage300_protocol.py',
        'scripts/stage288_common.py', 'scripts/stage289_decodiphy_bridge.py',
        'scripts/stage293_abundance.py', 'scripts/stage296_common.py',
        'scripts/stage298_common.py']
IMPORTED = ['results/stage298/stage298_analysis.json',
            'results/stage298/stage298_td_summary.json',
            'results/stage298/stage298_residuals.json',
            'results/stage298/stage298_placement_manifest.json']
HARVEST = ['results/stage289/stage289_headtohead_runs.jsonl',
           'results/stage293/stage293_decodiphy_runs.jsonl',
           'results/stage294/stage294_dd_runs.jsonl',
           'results/stage297/stage297_dd_seeds_main.jsonl']


def hashes(paths):
    return {p: s300.sha(ROOT / p) for p in paths}


def main():
    inputs = {}
    for collection in s300.COLLECTIONS:
        entries = s300.case_entries(collection)
        inputs[collection] = dict(
            n_cases=len(entries),
            cases=[dict(name=e['name'], path=e['path'], sha256=e['sha256'])
                   for e in entries])
    versions = subprocess.run(
        [str(VENV), '-c', 'import json, sys, numpy, scipy, osqp, decodiphy; '
         'from importlib.metadata import version; '
         'print(json.dumps(dict(python=sys.version.split()[0], numpy=numpy.__version__, '
         'scipy=scipy.__version__, osqp=version("osqp"), '
         'treeswift=version("treeswift"), decodiphy=decodiphy.__version__)))'],
        capture_output=True, text=True)
    protocol = dict(
        stage=300, step='own_fit_comparison',
        registered_utc=time.strftime('%Y-%m-%dT%H:%M:%S+00:00', time.gmtime()),
        question=('On real estimated distances, how does the mixture DecoDiPhy actually '
                  'returns fit the data, measured in the same units as the exact '
                  'single-source tolerance Theta_1?'),
        measurements=dict(
            r_dd=('R_dd = max_r |y + sum_i p_i dist_T(r, q_i) - d_r| over the observed '
                  'reference leaves, with the configuration exactly as the tool returned '
                  'it (bank edge index, p, x from the frozen stage-289 bridge; y = the '
                  'scalar pendant, which DecoDiPhy splits across placements only for its '
                  'jplace output and whose parts sum back to it by construction - '
                  'decodiphy/tree_utils.compute_individual_ys).  Exact rational.'),
            theta1='imported from stage 298; never recomputed differently',
            truth_linf=('imported from stage 298 (residual_max) and independently '
                        'recomputed here as a verification'),
            ssr=('sum_r (pred_r - d_r)^2, float64 - DecoDiPhy\'s own objective '
                 '(decodiphy/optimize.py: np.sum((d - (DA @ p + CAl @ x + y)) ** 2))'),
            lambda1=('exhaustive minimum of SSR over every original edge, x in [0,1], '
                     'y >= 0, p = 1: the L2 analogue of Theta_1, i.e. the best single '
                     'source in the norm the tool optimises')),
        norm_caveat=('DecoDiPhy minimises SSR, not the L-infinity residual.  Both norms '
                     'are therefore reported for every configuration, and the single-'
                     'source baseline is reported in both (Theta_1 and Lambda_1).'),
        collections=inputs,
        imported_stage298=hashes(IMPORTED),
        harvest_sources=hashes(HARVEST),
        code=hashes(CODE),
        decodiphy_package={('data/stage289_decodiphy/repo/decodiphy/' + p.name):
                           s300.sha(p) for p in sorted(REPO.glob('*.py'))},
        environment=dict(node=platform.platform(), venv=str(VENV),
                         versions=(json.loads(versions.stdout) if versions.returncode == 0
                                   else versions.stderr[-500:])),
        run_plan=dict(
            harvested=('every full-panel DecoDiPhy run recorded by stages 289, 293, 294 '
                       'and 297 on the krepp ANI family, with provenance (file, line, '
                       'arm, repetition or seed); nothing re-run'),
            new=dict(krepp_placement='60 cases x {dd_default, dd_true_k}, wall 1800 s, '
                                     'RLIMIT_AS 8 GB, 20 workers',
                     travelers='29 samples x {dd_default}, wall 14400 s, RLIMIT_AS 60 GB, '
                               'workers sized from the pilot\'s peak RSS'),
            retention='every run retained, including timeouts, memory kills and crashes',
            threads='OMP/OPENBLAS/MKL_NUM_THREADS=1 in every child',
            cores_budget=40),
        launched_before_registration=[
            dict(group='krepp_placement dd_default+dd_true_k', jobs=120,
                 out='results/stage300/stage300_runs_placement.jsonl'),
            dict(group='travelers pilot (1 sample, dd_default)', jobs=1,
                 out='results/stage300/stage300_pilot_travelers.jsonl'),
            dict(group='smoke test (2 placement cases)', jobs=4,
                 out='results/stage300/stage300_smoke.jsonl')],
        verification_plan=dict(
            route_a='exact rational equality of the min form and the interpolation form '
                    'on every evaluated run',
            route_b='whole-case recomputation with scipy.sparse.csgraph.dijkstra on the '
                    'edge list with the source points inserted as vertices',
            route_c='prediction rebuilt from DecoDiPhy\'s own D, C, l matrices on the '
                    'tree file the bridge writes',
            route_d='recomputed SSR of the returned configuration against the loss the '
                    'tool itself reported, for every run',
            route_e='Lambda_1 against the tool\'s own k = 1 ladder loss, which it '
                    'computes by hill climbing over the same edge set'))
    path = s300.OUT / 'stage300_protocol.json'
    path.write_text(json.dumps(protocol, indent=1))
    print('registered', path, 'cases:',
          {k: v['n_cases'] for k, v in inputs.items()})
    print('versions', protocol['environment']['versions'])


if __name__ == '__main__':
    main()
