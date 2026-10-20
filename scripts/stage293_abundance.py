#!/usr/bin/env python3
"""Stage 293: what reference incompleteness costs the ABUNDANCES and POSITIONS.

Stages 288/290 measured counts and support identity under reference loss.  This
stage measures the parameters a practitioner reports: the weights p_i and the
interior positions x_i, plus the mixture-level error and, when the data do not
pin the parameters down, the exact size of the set of equally valid answers.

Case panel, removal regimes, removal levels, seeds and draws are IMPORTED from
stage288_incompleteness (panel, removal_draws) so that every draw here is the
same draw stage 288 measured the count collapse on.  The tree is always the full
reference tree; reference loss removes OBSERVATIONS, never edges -- the stage 288
"B_marks_retained" convention.

Arms
  exact        solve_partial on the retained observations: minimum source count
               and one minimum-size witness (the answer a practitioner gets).
  fibre_true   exact fibre over the TRUE support given the retained data: is the
               correct support still feasible, is it still identifiable, and how
               far apart can two explanations on it be (stage293_fibre).
  fibre_found  the same analysis over the support the solver returned, when that
               differs from the true support.
  dd_default   DecoDiPhy 1.5.0, documented defaults, k chosen by its own rules.
  dd_true_k    DecoDiPhy 1.5.0 with -k K_true -f 1, the k asymmetry made explicit.

Nothing pre-existing is written.  Every run is retained, including timeouts.
"""
import os
for _key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
             'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_key] = '1'
import argparse
import contextlib
import hashlib
import io
import json
import platform
import shutil
import sys
import tempfile
import time
from datetime import datetime, timezone
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import stage288_common as common
import stage288_incompleteness as inc
from stage288_common import sha, serial, run_pool, load_case, verify_witness, predict
from stage293_fibre import fibre_report

OUTDIR = ROOT / 'results/stage293'
PROTOCOL = OUTDIR / 'stage293_protocol.json'
STREAMS = dict(exact=OUTDIR / 'stage293_exact_runs.jsonl',
               fibre_true=OUTDIR / 'stage293_fibre_true_runs.jsonl',
               fibre_found=OUTDIR / 'stage293_fibre_found_runs.jsonl',
               decodiphy=OUTDIR / 'stage293_decodiphy_runs.jsonl')

WORKERS = 20
EXACT_WALL, EXACT_AS = 60.0, 3 * 1024 ** 3
FIBRE_WALL, FIBRE_AS = 300.0, 3 * 1024 ** 3
DD_WALL, DD_AS = 120.0, 3 * 1024 ** 3
DD_METHODS = ('dd_default', 'dd_true_k')
VENV = ROOT / 'data/stage289_decodiphy/venv'
FROZEN = ('stage43_partial_fpt.py', 'stage41_partial_fpt.py', 'stage27_general_fpt.py',
          'stage26_missing_leaf.py', 'stage25_deficit_extension.py',
          'stage23_two_source_envelope.py', 'stage42_benchmark.py',
          'stage288_common.py', 'stage288_incompleteness.py',
          'stage289_support_feasibility.py', 'stage289_decodiphy_bridge.py',
          'stage289_headtohead_analyze.py')
NEW = ('stage293_fibre.py', 'stage293_abundance.py', 'stage293_validate_lp.py',
       'stage293_analyze.py')


# ------------------------------------------------------------------ job plan

def draw_plan():
    """Every (case, level, regime, replicate) this stage measures, plus the
    full-panel baseline.  Deterministic; identical across arms by construction."""
    plan = []
    for entry in inc.panel():
        case = load_case(entry)
        leaves = sorted(case['observed'])
        base = dict(name=entry['name'], family=entry['family'], cell=entry['cell'],
                    reference_set_size=entry['reference_set_size'],
                    K_generating=entry['K_generating'], cap=case['cap'],
                    n_leaves=len(leaves))
        plan.append((case, dict(base, level='none', regime='none', replicate=0,
                                removed=[], n_removed=0, n_retained=len(leaves),
                                removed_fraction=0.0)))
        for level in inc.LEVELS:
            for draw in inc.removal_draws(case, level):
                plan.append((case, dict(base, level=level, regime=draw['regime'],
                                        replicate=draw['replicate'],
                                        removed=draw['removed'],
                                        n_removed=draw['n_removed'],
                                        n_retained=draw['n_retained'],
                                        removed_fraction=draw['n_removed'] / len(leaves),
                                        detail=draw['detail'])))
    return plan


def retained_observations(case, removed):
    drop = set(removed)
    return {r: d for r, d in case['observed'].items() if r not in drop}


def plan_digest(plan):
    """Hash of the whole draw list so every arm can be shown to share it."""
    payload = [[m['name'], m['level'], m['regime'], m['replicate'], list(m['removed'])]
               for _, m in plan]
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()


# -------------------------------------------------------------------- workers

def exact_worker(spec):
    from stage43_partial_fpt import solve_partial
    start = time.perf_counter()
    raw = solve_partial(spec['edges'], spec['observed'], cap=spec['cap'], seconds=None,
                        max_states=None, max_bases=None, max_faces=None)
    answer = dict(status=raw['status'], feasible=raw.get('feasible'),
                  minimum=raw.get('minimum_sources_within_cap'), reason=raw.get('reason'),
                  marked_vertices=raw.get('marked_vertices'),
                  statistics=raw.get('statistics'),
                  solver_wall_seconds=time.perf_counter() - start)
    if 'original_witness' in raw:
        answer['witness'] = raw['original_witness']
        answer['witness_frame'] = 'original'
    return answer


def fibre_worker(spec):
    start = time.perf_counter()
    answer = fibre_report(spec['edges'], spec['observed'], spec['support'],
                          truth=spec.get('truth'), want_x=True)
    answer['solver_wall_seconds'] = time.perf_counter() - start
    answer['support'] = list(spec['support'])
    return answer


def dd_worker(spec):
    """One DecoDiPhy run on the retained distances; documented defaults only."""
    import decodiphy.cli as ddcli
    from stage289_decodiphy_bridge import write_inputs, to_bank_sources
    case = dict(edges=spec['edges'], observed=spec['observed'],
                generating_sources=spec['generating_sources'])
    with tempfile.TemporaryDirectory(prefix='stage293_dd_') as work:
        mapping = write_inputs(case, work)
        outdir = Path(work) / spec['method']
        outdir.mkdir(parents=True, exist_ok=True)
        base = ['-t', str(Path(work) / 'tree.nwk'), '-d', str(Path(work) / 'distances.txt'),
                '-o', str(outdir)]
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
            # DecoDiPhy 1.5.0 writes all_rounds.json (cli.py:379) BEFORE the jplace
            # writer (cli.py:384), and save_jplace looks its anchors up in the PRUNED
            # tree.  When post_process_output maps a placement back onto an original
            # edge that pruning contracted away, the writer raises KeyError although
            # the fit itself finished and is already on disk.  Recover the tool's own
            # committed placements; nothing is re-fitted and nothing is tuned.
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
            # an anchor with no parent edge in the original tree (e.g. the root):
            # retained as an unusable answer rather than silently dropped
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


# -------------------------------------------------- independent verification

def make_exact_post(by_key):
    """Recompute the leaf equations on the ORIGINAL tree; trust no solver flag."""
    def post(job, answer):
        meta = job['meta']
        case, retained = by_key[(meta['name'], meta['level'], meta['regime'], meta['replicate'])]
        witness = answer.get('witness')
        if not witness:
            return
        check = verify_witness(case['edges'], retained, meta['cap'], witness)
        answer['independent_witness_verification'] = check
        answer['independent_witness_verified'] = check['verified']
        answer['witness_count_matches_minimum'] = check.get('source_count') == answer.get('minimum')
        support = sorted({int(s['edge']) for s in witness['sources']})
        truth = {int(s['edge']) for s in case['generating_sources']}
        answer['support'] = support
        answer['support_exact_match'] = set(support) == truth
        answer['n_true_edges_recovered'] = len(set(support) & truth)
        # does the witness also reproduce the references that were REMOVED?
        removed = sorted(set(case['observed']) - set(retained))
        if removed:
            sources = [dict(edge=int(s['edge']), p=F(s['p']), x=F(s['x']))
                       for s in witness['sources']]
            pred = predict(case['edges'], sources, F(witness['pendant']), removed)
            deviation = max(abs(pred[r] - case['observed'][r]) for r in removed)
            answer['max_abs_deviation_on_removed_references'] = str(deviation)
            answer['explains_removed_references'] = deviation == 0
        if answer['support_exact_match']:
            got = {int(s['edge']): (F(s['p']), F(s['x'])) for s in witness['sources']}
            true = {int(s['edge']): (F(s['p']), F(s['x'])) for s in case['generating_sources']}
            answer['p_abs_error'] = {str(e): str(abs(got[e][0] - true[e][0])) for e in true}
            answer['x_abs_error'] = {str(e): str(abs(got[e][1] - true[e][1])) for e in true}
            answer['max_p_abs_error'] = str(max(abs(got[e][0] - true[e][0]) for e in true))
            answer['max_x_abs_error'] = str(max(abs(got[e][1] - true[e][1]) for e in true))
            answer['total_variation_p'] = str(sum(abs(got[e][0] - true[e][0]) for e in true) / 2)
            answer['pendant_abs_error'] = str(abs(F(witness['pendant'])
                                                  - case['generating_pendant']))
    return post


# ------------------------------------------------------------------ protocol

def prepare():
    assert not PROTOCOL.exists(), 'Never overwrite a registered protocol'
    OUTDIR.mkdir(parents=True, exist_ok=True)
    plan = draw_plan()
    cases = sorted({m['name'] for _, m in plan})
    by_level = {}
    for _, m in plan:
        by_level.setdefault(m['level'], {}).setdefault(m['regime'], 0)
        by_level[m['level']][m['regime']] += 1
    protocol = dict(
        stage=293, track='abundance_and_position_cost_of_reference_incompleteness',
        registered_utc=datetime.now(timezone.utc).isoformat(),
        question=('When references go missing, how wrong are the inferred abundances p and '
                  'interior positions x -- separately from the count collapse of stage 288?'),
        design=dict(
            panel='imported stage288_incompleteness.panel()',
            draws='imported stage288_incompleteness.removal_draws(); seeds unchanged',
            seed=inc.SEED, levels=list(inc.LEVELS), regimes=list(inc.REGIMES),
            replicates=dict(inc.REPLICATES),
            removal_convention=('observations are removed; the reference tree keeps every edge '
                                '(stage 288 variant B)'),
            baseline='level "none" is the full panel, the no-loss reference point',
            arms=dict(
                exact='stage43_partial.solve_partial on the retained observations',
                fibre_true=('exact fibre over the true support given the retained data; '
                            'LP bounds in Fractions, suprema over the closure'),
                fibre_found='the same over the support the exact solver returned, when different',
                dd_default='DecoDiPhy 1.5.0 defaults, k by its own stopping rules',
                dd_true_k='DecoDiPhy 1.5.0 with -k K_true -f 1'),
            wasserstein=('reusing stage289_headtohead_analyze.wasserstein together with '
                         'point_distances: balanced optimal transport between the recovered and '
                         'generating mixtures under the tree metric, computed in the analysis '
                         'step so the stage 289 numbers are directly comparable'),
            limits=dict(workers=WORKERS, exact_wall_seconds=EXACT_WALL,
                        fibre_wall_seconds=FIBRE_WALL, decodiphy_wall_seconds=DD_WALL,
                        address_space_bytes=EXACT_AS),
            decodiphy_tuning='none; documented defaults, only the k control named in the arm',
            decodiphy_missing_references=(
                'DecoDiPhy handles missing references itself: preprocess_input prunes leaves '
                'absent from the distance file and contracts their parents, and '
                'post_process_output maps placements back to the original tree. Its jplace '
                'writer (cli.py:384) then looks those anchors up in the PRUNED tree and can '
                'raise KeyError; all_rounds.json is already written at cli.py:379, so this '
                'stage reads the placements the tool itself committed and records '
                'jplace_writer_crashed. No refit, no tuning, no fallback method.')),
        n_cases=len(cases), n_jobs_per_arm=len(plan), draws_by_level_and_regime=by_level,
        plan_sha256=plan_digest(plan),
        predictions=[
            'H1 the abundance error grows with the removed fraction and is worse under the '
            'structured regimes, mirroring the count collapse of stage 288.',
            'H2 a retained, CORRECT support can still carry large abundance error, because the '
            'fibre over the true support gains free parameters as references are removed; this '
            'is a separate effect from the count collapse.',
            'H3 positive-dimensional fibres are not rare and their diameters are large enough '
            'that a single reported p is not a well-defined quantity.',
            'H4 no prediction is registered for which of the exact method and the heuristic '
            'degrades more gracefully.'],
        inputs=dict(bank_index=str(sha(ROOT / 'data/stage287_instances/bank/index.json')),
                    bank_protocol=str(sha(ROOT / 'results/stage287_instance_bank_protocol.json')),
                    stage288_protocol=str(sha(ROOT / 'results/stage288_incompleteness_protocol.json'))),
        code_hashes={name: sha(ROOT / 'scripts' / name)
                     for name in sorted(set(FROZEN + NEW)) if (ROOT / 'scripts' / name).exists()},
        dependency_closure=common.dependency_closure('stage293_abundance.py'),
        environment=dict(python=platform.python_version(), platform=platform.platform(),
                         project_venv=str(ROOT / '.venv/bin/python'),
                         decodiphy_venv=str(VENV / 'bin/python'),
                         decodiphy_version='1.5.0'))
    PROTOCOL.write_text(json.dumps(protocol, indent=1))
    print(json.dumps(dict(cases=len(cases), jobs_per_arm=len(plan),
                          plan_sha256=protocol['plan_sha256'],
                          by_level={k: sum(v.values()) for k, v in by_level.items()}), indent=1))


def _check_protocol(plan):
    registered = json.loads(PROTOCOL.read_text())
    assert registered['plan_sha256'] == plan_digest(plan), 'draw plan differs from the protocol'
    return registered


# ----------------------------------------------------------------------- runs

def run_exact():
    stream = STREAMS['exact']
    assert not stream.exists(), 'Never overwrite measured runs'
    plan = draw_plan()
    _check_protocol(plan)
    jobs, by_key = [], {}
    for case, meta in plan:
        retained = retained_observations(case, meta['removed'])
        key = (meta['name'], meta['level'], meta['regime'], meta['replicate'])
        by_key[key] = (case, retained)
        jobs.append(dict(meta=dict(meta, arm='exact'),
                         spec=dict(edges=case['edges'], observed=retained, cap=meta['cap'])))
    print(f'exact arm: {len(jobs)} jobs', flush=True)
    run_pool(jobs, exact_worker, EXACT_WALL, EXACT_AS, WORKERS, stream,
             post=make_exact_post(by_key), progress=100)


def run_fibre_true():
    stream = STREAMS['fibre_true']
    assert not stream.exists(), 'Never overwrite measured runs'
    plan = draw_plan()
    _check_protocol(plan)
    jobs = []
    for case, meta in plan:
        retained = retained_observations(case, meta['removed'])
        support = sorted(int(s['edge']) for s in case['generating_sources'])
        jobs.append(dict(meta=dict(meta, arm='fibre_true', support=support),
                         spec=dict(edges=case['edges'], observed=retained, support=support,
                                   truth=case['generating_sources'])))
    print(f'fibre_true arm: {len(jobs)} jobs', flush=True)
    run_pool(jobs, fibre_worker, FIBRE_WALL, FIBRE_AS, WORKERS, stream, progress=100)


def run_fibre_found():
    """Fibre over the support the exact solver actually returned, when it differs."""
    stream = STREAMS['fibre_found']
    assert not stream.exists(), 'Never overwrite measured runs'
    plan = draw_plan()
    _check_protocol(plan)
    found = {}
    with open(STREAMS['exact']) as handle:
        for line in handle:
            row = json.loads(line)
            support = row['result'].get('support')
            if support and not row['result'].get('support_exact_match'):
                found[(row['name'], row['level'], row['regime'], row['replicate'])] = support
    jobs = []
    for case, meta in plan:
        key = (meta['name'], meta['level'], meta['regime'], meta['replicate'])
        if key not in found:
            continue
        retained = retained_observations(case, meta['removed'])
        jobs.append(dict(meta=dict(meta, arm='fibre_found', support=found[key]),
                         spec=dict(edges=case['edges'], observed=retained,
                                   support=found[key], truth=None)))
    print(f'fibre_found arm: {len(jobs)} jobs', flush=True)
    run_pool(jobs, fibre_worker, FIBRE_WALL, FIBRE_AS, WORKERS, stream, progress=100)


def run_decodiphy():
    stream = STREAMS['decodiphy']
    assert not stream.exists(), 'Never overwrite measured runs'
    plan = draw_plan()
    _check_protocol(plan)
    jobs = []
    for case, meta in plan:
        retained = retained_observations(case, meta['removed'])
        for method in DD_METHODS:
            jobs.append(dict(meta=dict(meta, arm=method),
                             spec=dict(edges=case['edges'], observed=retained, method=method,
                                       k_true=meta['K_generating'],
                                       generating_sources=case['generating_sources'])))
    print(f'decodiphy arms: {len(jobs)} jobs', flush=True)
    run_pool(jobs, dd_worker, DD_WALL, DD_AS, WORKERS, stream, progress=200)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('command', choices=['prepare', 'run-exact', 'run-fibre-true',
                                            'run-fibre-found', 'run-decodiphy'])
    name = parser.parse_args().command
    {'prepare': prepare, 'run-exact': run_exact, 'run-fibre-true': run_fibre_true,
     'run-fibre-found': run_fibre_found, 'run-decodiphy': run_decodiphy}[name]()
