#!/usr/bin/env python3
"""Stage297 item 2: DecoDiPhy 1.5.0 seed sensitivity at the full panel.

Stage289 ran the tool only at its documented default seed (1142) and recorded
seed sensitivity as unmeasured.  Stage294 then found that the SUPPORT the tool
returns is not reproducible across harness process details (its
stage289_full_panel_agreement counts: 575 compared, 567 agree on k, 471 agree on
support), which makes the gap sharper: part of what stages 289/293/294 attribute
to the heuristic may be run-to-run variation rather than method error.

This stage measures it directly.  Each (case, arm) pair is one forked child, so
every run is a separate process exactly as in stage 293.

  arm 'default_1142'   `stage293_abundance.dd_worker` IMPORTED UNMODIFIED, no -s
                       flag, so DecoDiPhy takes its documented default seed
                       (cli.py:39, default=1142).  This is byte-identical to the
                       invocation of stages 289/293/294.
  arm 'seeded_<s>'     `dd_worker_seeded`, a verbatim copy of that worker whose
                       ONLY change is `-s <seed>` appended to the argument list.
                       The unified diff between the two function bodies is
                       written into the protocol, and arm 'seeded_1142' repeats
                       the default seed through the copy: any difference between
                       'default_1142' and 'seeded_1142' is process-level
                       nondeterminism at a fixed seed, not seed sensitivity.

Only `-s` ever differs; `--method hill --optimizer osqp` and every other
documented default are untouched, and k is always chosen by the tool's own
stopping rules (method 'dd_default'), which is the quantity whose spread is in
question.  W1 to truth is NOT computed here: the workers record the tool's own
placements and the analysis step recomputes W1 with
`stage289_headtohead_analyze.wasserstein` on them.
"""
import contextlib
import difflib
import inspect
import io
import json
import os
import shutil
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
sys.dont_write_bytecode = True

import stage297_common as common                                          # noqa: E402
from stage297_common import OUTDIR, run_pool                              # noqa: E402
from stage293_abundance import dd_worker                                  # noqa: E402

WORKERS = int(os.environ.get('STAGE297_DD_WORKERS', '8'))
DD_WALL = float(os.environ.get('STAGE297_DD_WALL', '900'))
DD_AS = 3 * 1024 ** 3
SEEDS = (1142, 1, 7, 101, 2024, 31337)
DEFAULT_SEED = 1142


def dd_worker_seeded(spec):
    """Verbatim copy of stage293_abundance.dd_worker; ONLY `-s <seed>` added."""
    import decodiphy.cli as ddcli
    from stage289_decodiphy_bridge import write_inputs, to_bank_sources
    case = dict(edges=spec['edges'], observed=spec['observed'],
                generating_sources=spec['generating_sources'])
    with tempfile.TemporaryDirectory(prefix='stage297_dd_') as work:
        mapping = write_inputs(case, work)
        outdir = Path(work) / spec['method']
        outdir.mkdir(parents=True, exist_ok=True)
        base = ['-t', str(Path(work) / 'tree.nwk'), '-d', str(Path(work) / 'distances.txt'),
                '-o', str(outdir), '-s', str(spec['seed'])]
        if spec['method'] == 'dd_default':
            argv = base
        elif spec['method'] == 'dd_true_k':
            argv = base + ['-k', str(spec['k_true']), '-f', '1']
        else:
            raise ValueError(spec['method'])
        saved, buf = sys.argv, io.StringIO()
        sys.argv = ['decodiphy'] + argv
        start = time.perf_counter()
        crash = None
        try:
            with contextlib.redirect_stdout(buf):
                ddcli.main()
        except BaseException as exc:
            crash = repr(exc)[:300]
        finally:
            sys.argv = saved
        seconds = time.perf_counter() - start
        log = buf.getvalue()
        rounds_path = outdir / 'all_rounds.json'
        if not rounds_path.exists():
            raise RuntimeError('decodiphy produced no all_rounds.json; crash=%s' % crash)
        rounds = json.loads(rounds_path.read_text())
        final = rounds[-1]
        unmapped = [a for a in final['anchors'] if a not in mapping['edge_of_label']]
        if unmapped:
            shutil.rmtree(outdir, ignore_errors=True)
            return dict(status='unmapped_anchor', argv=argv, unmapped_anchors=unmapped,
                        anchors=list(final['anchors']), k=len(final['anchors']),
                        loss=final['loss'], jplace_writer_crashed=crash is not None,
                        writer_exception=crash, solver_wall_seconds=seconds)
        sources = to_bank_sources(case, final['anchors'], final['p'], final['x'], mapping)
        y = final['y']
        stop = ('low_abundance' if 'Stopped by low abundance' in log else
                'near_zero_residue' if 'Stopped by near-zero residue' in log else 'k_exhausted')
        shutil.rmtree(outdir, ignore_errors=True)
    return dict(status='complete', argv=argv, k=len(sources), loss=final['loss'],
                sources=sources, stop_reason=stop,
                jplace_writer_crashed=crash is not None, writer_exception=crash,
                pendant_total=float(sum(y)) if isinstance(y, list) else float(y),
                ladder=[dict(k=r['k'], rounds=r['rounds'], loss=r['loss']) for r in rounds],
                solver_wall_seconds=seconds)


def worker_diff():
    """Unified diff of the two worker bodies, for the protocol."""
    a = inspect.getsource(dd_worker).splitlines()
    b = inspect.getsource(dd_worker_seeded).splitlines()
    return [l for l in difflib.unified_diff(a, b, 'stage293.dd_worker',
                                            'stage297.dd_worker_seeded', lineterm='', n=1)]


def job(spec):
    entry = spec['entry']
    case = common.load_verified(entry)
    payload = dict(edges=case['edges'], observed=case['observed'],
                   generating_sources=case['generating_sources'],
                   method='dd_default', k_true=entry['K_generating'],
                   seed=spec['seed'])
    t0 = time.perf_counter()
    if spec['arm'] == 'default_1142':
        answer = dd_worker({k: v for k, v in payload.items() if k != 'seed'})
    else:
        answer = dd_worker_seeded(payload)
    answer['harness_wall_seconds'] = time.perf_counter() - t0
    truth = sorted({int(s['edge']) for s in case['generating_sources']})
    support = (sorted(int(s['edge']) for s in answer['sources'])
               if answer.get('status') == 'complete' else None)
    row = dict(kind='dd_seed', stage=297, name=entry['name'], family=entry['family'],
               path=entry['path'], sha256=entry['sha256'],
               K_generating=entry['K_generating'], R_size=entry['reference_set_size'],
               tree_leaves=entry['tree_leaves'], n_edges=len(case['edges']),
               distance_source=entry['distance_source'],
               exactly_explained=entry['exactly_explained'],
               arm=spec['arm'], seed=spec['seed'], true_support=truth,
               status=answer.get('status'), k=answer.get('k'),
               loss=answer.get('loss'), stop_reason=answer.get('stop_reason'),
               support=support,
               support_exact_match=(None if support is None else set(support) == set(truth)),
               n_true_edges_recovered=(None if support is None
                                       else len(set(support) & set(truth))),
               repeated_edge=(None if support is None
                              else len(set(support)) != len(support)),
               pendant_total=answer.get('pendant_total'),
               jplace_writer_crashed=answer.get('jplace_writer_crashed'),
               writer_exception=answer.get('writer_exception'),
               solver_wall_seconds=answer.get('solver_wall_seconds'),
               harness_wall_seconds=answer['harness_wall_seconds'],
               argv_tail=[a for a in answer.get('argv', []) if not a.startswith('/')],
               sources=answer.get('sources'), ladder=answer.get('ladder'),
               unmapped_anchors=answer.get('unmapped_anchors'))
    return [row]


def plan(per_cell, families=('exact_10k', 'exact_biotrees', 'krepp_estimated_10k'),
         seeds=SEEDS):
    """`per_cell` cases for every (family, K) cell, stratified over |R|."""
    jobs = []
    chosen = []
    for family in families:
        entries = [e for e in common.bank_index() if e['family'] == family]
        for K in sorted({e['K_generating'] for e in entries}):
            cell = [e for e in entries if e['K_generating'] == K]
            chosen.extend(common.stratified(cell, per_cell))
    chosen.sort(key=lambda e: (e['family'], e['K_generating'], e['reference_set_size']))
    for entry in chosen:
        jobs.append(dict(entry=entry, arm='default_1142', seed=DEFAULT_SEED))
        for s in seeds:
            jobs.append(dict(entry=entry, arm='seeded_%d' % s, seed=s))
    return jobs


def main():
    per_cell = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    tag = sys.argv[2] if len(sys.argv) > 2 else 'main'
    jobs = plan(per_cell)
    (OUTDIR).mkdir(parents=True, exist_ok=True)
    (OUTDIR / 'stage297_dd_worker_diff.json').write_text(
        json.dumps(dict(diff=worker_diff(), seeds=list(SEEDS)), indent=1))
    run_pool(jobs, job, WORKERS, DD_WALL, DD_AS,
             OUTDIR / ('stage297_dd_seeds_%s.jsonl' % tag), label='dd_seeds', every=25)


if __name__ == '__main__':
    main()
