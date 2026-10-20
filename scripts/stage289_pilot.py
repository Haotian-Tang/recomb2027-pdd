#!/usr/bin/env python3
"""Stage 289 pilot: time every arm on a handful of bank cases (no protocol yet)."""
import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
            'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import contextlib, io, json, sys, time
from fractions import Fraction as F
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from stage289_decodiphy_bridge import load_case, write_inputs, to_bank_sources
from stage289_support_feasibility import feasible_support
from stage43_partial_fpt import solve_partial

import decodiphy.cli as ddcli


def run_decodiphy(argv, outdir):
    saved = sys.argv
    sys.argv = ['decodiphy'] + argv
    start = time.perf_counter()
    try:
        with contextlib.redirect_stdout(io.StringIO()) as buf:
            ddcli.main()
    finally:
        sys.argv = saved
    rounds = json.loads((Path(outdir) / 'all_rounds.json').read_text())
    return rounds, time.perf_counter() - start, buf.getvalue()


def pilot(case_path, work):
    case = load_case(case_path)
    k_true = len(case['generating_sources'])
    work = Path(work)
    mapping = write_inputs(case, work)
    tree, dist = str(work / 'tree.nwk'), str(work / 'distances.txt')
    report = dict(name=case['name'], leaves=len(case['observed']), edges=len(case['edges']),
                  k_true=k_true, cap=case['cap'])
    arms = {
        'dd_default': [],
        'dd_true_k': ['-k', str(k_true), '-f', '1'],
        'dd_resolution': ['--min_p', '1e-12', '-k', str(k_true + 3), '-f', '0'],
    }
    for name, extra in arms.items():
        outdir = work / name
        outdir.mkdir(exist_ok=True)
        try:
            rounds, seconds, _ = run_decodiphy(['-t', tree, '-d', dist, '-o', str(outdir)] + extra, outdir)
            final = rounds[-1]
            sources = to_bank_sources(case, final['anchors'], final['p'], final['x'], mapping)
            verdict = feasible_support(case['edges'], case['observed'], [s['edge'] for s in sources])
            report[name] = dict(seconds=round(seconds, 3), k=final['k'], loss=final['loss'],
                                verdict=verdict['verdict'])
        except Exception as exc:
            report[name] = dict(error=repr(exc)[:200])
    start = time.perf_counter()
    answer = solve_partial(case['edges'], case['observed'], cap=case['cap'],
                           seconds=float(os.environ.get('EXACT_SECONDS', '120')),
                           max_states=None, max_bases=None, max_faces=None)
    report['exact'] = dict(seconds=round(time.perf_counter() - start, 3), status=answer['status'],
                           feasible=answer.get('feasible'),
                           minimum=answer.get('minimum_sources_within_cap'),
                           reason=answer.get('reason'))
    return report


if __name__ == '__main__':
    out = []
    for path in sys.argv[1:]:
        work = Path('/tmp/stage289_pilot') / Path(path).stem
        work.mkdir(parents=True, exist_ok=True)
        out.append(pilot(path, work))
        print(json.dumps(out[-1]), flush=True)
