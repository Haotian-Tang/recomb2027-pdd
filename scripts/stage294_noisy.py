#!/usr/bin/env python3
"""Stage 294: does the reference-incompleteness count collapse survive ESTIMATION ERROR?

Every statement stages 288/290/292/293 make about reference incompleteness is
measured on exactly explainable instances in exact rational arithmetic.  On real
krepp ANI distances the exact solver correctly answers "no exact explanation"
and stops (stage 289), so the question must be put to the tool that does run on
noisy data, against published ground truth.

Arms
  dd_default      DecoDiPhy 1.5.0, documented defaults, k chosen by its own
                  stopping rules.  THE ARM UNDER TEST: the inferred source count
                  is the quantity reference removal is supposed to depress.
  dd_true_k       DecoDiPhy 1.5.0 with -k K_true -f 1, so support and abundance
                  accuracy can be read separately from k selection.
  exact_matched   the frozen exact solver on the PAIRED exactly explainable
                  instance -- same tree, same sources, same retained reference
                  set, exact tree-derived distances instead of ANI estimates.
                  Restricted to K in {3, 5}: stage 292 established that the
                  frozen minimum pass does not complete at K = 10 even at a
                  3,600 s budget, so K = 10 is not attempted rather than burned
                  as timeouts.

Both DecoDiPhy workers and the exact worker are IMPORTED from stage293_abundance
unchanged, so the tool invocation is byte-identical to the exact-data stage and
the two sets of runs are directly comparable.

The baseline that matters is NOT the truth.  DecoDiPhy already mis-estimates k
on 273 of 290 full-panel krepp cases (stage 289), so this stage measures the
CHANGE induced by reference removal relative to the same tool's own answer on
the same case at the full panel: level "none" is that reference point and is
run for every case.

Nothing pre-existing is written.  Every run is retained, including failures.
"""
import os
for _key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
             'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_key] = '1'
import argparse
import json
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import stage288_common as common                       # noqa: E402
import stage288_incompleteness as inc                  # noqa: E402
import stage294_common as s294                         # noqa: E402
from stage288_common import sha, run_pool              # noqa: E402
from stage293_abundance import dd_worker, exact_worker, make_exact_post   # noqa: E402

OUTDIR = ROOT / 'results/stage294'
PROTOCOL = OUTDIR / 'stage294_protocol.json'
STREAMS = dict(dd=OUTDIR / 'stage294_dd_runs.jsonl',
               exact_matched=OUTDIR / 'stage294_exact_matched_runs.jsonl')
PILOT = OUTDIR / 'stage294_pilot.json'

WORKERS = 30                      # <= 32 cores, one thread per child
DD_WALL, DD_AS = 120.0, 3 * 1024 ** 3        # stage 293 limits
EXACT_WALL, EXACT_AS = 30.0, 4 * 1024 ** 3   # stage 288 variant-B limits
DD_METHODS = ('dd_default', 'dd_true_k')
EXACT_K = (3, 5)

VENV = ROOT / 'data/stage289_decodiphy/venv'
FROZEN = ('stage43_partial_fpt.py', 'stage41_partial_fpt.py', 'stage27_general_fpt.py',
          'stage26_missing_leaf.py', 'stage25_deficit_extension.py',
          'stage23_two_source_envelope.py', 'stage42_benchmark.py',
          'stage288_common.py', 'stage288_incompleteness.py',
          'stage289_decodiphy_bridge.py', 'stage289_headtohead_analyze.py',
          'stage293_abundance.py', 'stage293_fibre.py')
NEW = ('stage294_common.py', 'stage294_noisy.py', 'stage294_taxonomy.py',
       'stage294_analyze.py')
TAXONOMY_MANIFEST = s294.TAXDIR / 'stage294_taxonomy_manifest.json'


# ----------------------------------------------------------------- job build

def dd_jobs(plan):
    jobs = []
    for noisy, _exact, meta in plan:
        retained = s294.retained_observations(noisy, meta['removed'])
        for method in DD_METHODS:
            jobs.append(dict(meta=dict(meta, arm=method),
                             spec=dict(edges=noisy['edges'], observed=retained, method=method,
                                       k_true=meta['K_generating'],
                                       generating_sources=noisy['generating_sources'])))
    return jobs


def exact_jobs(plan):
    """The paired exact instance on the SAME retained reference set.

    The noisy panel is a subset of the exact panel (krepp returns no ANI
    estimate for 0-3 references per case); the exact arm is restricted to the
    noisy panel first so that the two arms see the same references.
    """
    jobs, by_key = [], {}
    for noisy, exact, meta in plan:
        if meta['K_generating'] not in EXACT_K:
            continue
        keep = set(noisy['observed']) - set(meta['removed'])
        retained = {r: d for r, d in exact['observed'].items() if r in keep}
        key = (meta['name'], meta['level'], meta['regime'], meta['replicate'])
        by_key[key] = (exact, retained)
        jobs.append(dict(meta=dict(meta, arm='exact_matched'),
                         spec=dict(edges=exact['edges'], observed=retained, cap=meta['cap'])))
    return jobs, by_key


# ------------------------------------------------------------------ protocol

def prepare():
    assert not PROTOCOL.exists(), 'Never overwrite a registered protocol'
    assert TAXONOMY_MANIFEST.exists(), 'Register the taxonomy manifest first'
    OUTDIR.mkdir(parents=True, exist_ok=True)
    ranks = s294.load_ranks()
    plan = s294.draw_plan(ranks)
    manifest = json.loads(TAXONOMY_MANIFEST.read_text())
    by_cell = {}
    for _, _, meta in plan:
        by_cell.setdefault(meta['level'], {}).setdefault(meta['regime'], 0)
        by_cell[meta['level']][meta['regime']] += 1
    n_exact = sum(1 for _, _, m in plan if m['K_generating'] in EXACT_K)
    protocol = dict(
        stage=294, track='noisy_reference_incompleteness',
        registered_utc=datetime.now(timezone.utc).isoformat(),
        question=('Does the downward bias that reference incompleteness puts on the inferred '
                  'constituent count survive real estimation error, and does structured loss '
                  'still hurt more than random loss, when the inference is done by the only '
                  'tool that runs on noisy data?'),
        design=dict(
            panel=('the complete krepp_estimated_10k family, 296 cases, real krepp ANI '
                   'distances estimated from reads; no subsampling'),
            paired_exact=('each case is paired 1:1 by instance stem with the exact_10k case '
                          'built from the same tree, the same generating sources and the same '
                          'abundances; verified identical edges and generating_sources in all '
                          '296 pairs, identical observed leaf sets in 135 (the rest differ only '
                          'by the 1-3 references for which krepp returned no ANI estimate, and '
                          'the noisy panel is a subset of the exact panel in all 296)'),
            draws=('levels, regimes, replicate counts, the removal-count rule and the three '
                   'synthetic regimes are imported from stage288_incompleteness / '
                   'stage288_common; the random seed key uses the PAIRED EXACT case name so a '
                   'random draw is byte-identical to the stage 288 draw wherever the two panels '
                   'agree'),
            seed=inc.SEED, levels=list(inc.LEVELS),
            regimes=list(inc.REGIMES) + ['taxon'],
            replicates=dict(dict(inc.REPLICATES), taxon=len(s294.TAXON_ORDERS)),
            taxon_regime=dict(
                rule=('drop WHOLE taxonomic groups at the finest operable rank (genus, else '
                      'species), taking groups largest-first (replicate 0) or smallest-first '
                      '(replicate 1) while neither overshooting the target removal count nor '
                      'leaving fewer than three references'),
                source=manifest,
                granularity_note=('whole-group granularity means the achieved removal count can '
                                  'fall short of the nominal level; both are recorded and the '
                                  'analysis bins on the achieved fraction as well as the level'),
                panel_breadth=('these benchmark panels are taxonomically narrow -- a median of '
                               'one phylum, one class, one order and one family per case -- so '
                               'no rank above genus is operable; this is a property of the '
                               'published 10k benchmark, reported rather than worked around')),
            removal_convention=('observations are removed; the reference tree keeps every edge. '
                                'DecoDiPhy prunes the leaves absent from its distance file '
                                'itself (stage 293 convention, stage 288 variant B)'),
            baseline=('level "none" is the full panel. The measured quantity is the CHANGE '
                      'relative to the same tool, same arm, same case at level "none" -- not '
                      'the deviation from the published truth, which is already large at the '
                      'full panel'),
            arms=dict(
                dd_default='DecoDiPhy 1.5.0 documented defaults, k by its own stopping rules',
                dd_true_k='DecoDiPhy 1.5.0 with -k K_true -f 1',
                exact_matched=('frozen stage43_partial.solve_partial on the paired exact '
                               'instance restricted to the same retained references; '
                               'K in %s only, because stage 292 established the frozen minimum '
                               'pass does not complete at K = 10 within 3,600 s' % (EXACT_K,))),
            workers_imported_verbatim=('dd_worker, exact_worker and make_exact_post are imported '
                                       'from stage293_abundance without modification'),
            determinism=('DecoDiPhy 1.5.0 at seed 1142 is deterministic: the stage 289 krepp '
                         'runs agree in k, loss and support across both repetitions in 575 of '
                         '575 cases where both completed, so one run per draw is recorded'),
            wasserstein=('stage289_headtohead_analyze.wasserstein with point_distances, applied '
                         'in the analysis step, so the numbers are comparable to stage 289'),
            limits=dict(workers=WORKERS, dd_wall_seconds=DD_WALL, dd_address_space_bytes=DD_AS,
                        exact_wall_seconds=EXACT_WALL, exact_address_space_bytes=EXACT_AS),
            decodiphy_tuning='none; documented defaults, only the k control named in the arm',
            failure_policy=('every run is retained. DecoDiPhy 1.5.0 writes all_rounds.json '
                            '(cli.py:379) before its jplace writer (cli.py:384), and that writer '
                            'raises on anchors that pruning contracted away; a run whose fit '
                            'finished and whose all_rounds.json exists is recorded as a result '
                            'with jplace_writer_crashed set, exactly as in stage 293. A run with '
                            'no all_rounds.json is recorded as an error with its exception')),
        n_cases=len({m['name'] for _, _, m in plan}),
        n_draws=len(plan), n_dd_jobs=2 * len(plan), n_exact_matched_jobs=n_exact,
        draws_by_level_and_regime=by_cell,
        plan_sha256=s294.plan_digest(plan),
        predictions=[
            'H1 the inferred count falls as references are removed, relative to the same tool '
            'on the same case at the full panel, at a rate comparable to the 22.9 per cent of '
            'draws / 59.1 per cent at 75 per cent loss stage 288 measured exactly.',
            'H2 structured loss (clade, nearest, taxon) depresses the count more than random '
            'loss, as it does on exact data (about threefold in draw frequency).',
            'H3 no prediction is registered for whether estimation error amplifies or masks the '
            'effect; a null result is reported as plainly as a positive one, and would mean the '
            'thesis is a statement about the exact model and must be scoped to it.'],
        inputs=dict(
            bank_index=sha(ROOT / 'data/stage287_instances/bank/index.json'),
            bank_protocol=sha(ROOT / 'results/stage287_instance_bank_protocol.json'),
            stage288_protocol=sha(ROOT / 'results/stage288_incompleteness_protocol.json'),
            stage289_headtohead_protocol=sha(
                ROOT / 'results/stage289/stage289_headtohead_protocol.json'),
            stage293_protocol=sha(ROOT / 'results/stage293/stage293_protocol.json'),
            taxonomy_manifest=sha(TAXONOMY_MANIFEST)),
        code_hashes={name: sha(ROOT / 'scripts' / name)
                     for name in sorted(set(FROZEN + NEW)) if (ROOT / 'scripts' / name).exists()},
        dependency_closure=common.dependency_closure('stage294_noisy.py'),
        environment=dict(python=platform.python_version(), platform=platform.platform(),
                         project_venv=str(ROOT / '.venv/bin/python'),
                         decodiphy_venv=str(VENV / 'bin/python'),
                         decodiphy_version='1.5.0'))
    PROTOCOL.write_text(json.dumps(protocol, indent=1) + '\n')
    print(json.dumps(dict(cases=protocol['n_cases'], draws=protocol['n_draws'],
                          dd_jobs=protocol['n_dd_jobs'],
                          exact_jobs=protocol['n_exact_matched_jobs'],
                          plan_sha256=protocol['plan_sha256']), indent=1))


def _check(plan):
    registered = json.loads(PROTOCOL.read_text())
    assert registered['plan_sha256'] == s294.plan_digest(plan), 'draw plan differs from protocol'
    amendment = OUTDIR / 'stage294_protocol_amendment.json'
    amended = {}
    if amendment.exists():
        amended = {Path(k).name: v.get('sha256_after_amendment')
                   for k, v in json.loads(amendment.read_text())['files'].items()}
    for name, digest in registered['code_hashes'].items():
        if name in ('stage294_analyze.py',):
            continue
        current = sha(ROOT / 'scripts' / name)
        if current == digest:
            continue
        assert current == amended.get(name), f'source changed since protocol: {name}'
    return registered


# --------------------------------------------------------------------- runs

def pilot():
    """A small timing pilot: never overwritten runs, written to its own file."""
    assert not PILOT.exists(), 'Never overwrite a registered pilot'
    plan = s294.draw_plan()
    picked = []
    seen = set()
    for noisy, exact, meta in plan:
        key = (meta['K_generating'], meta['level'], meta['regime'])
        if key in seen or meta['level'] not in ('none', '0.50', '0.75'):
            continue
        seen.add(key)
        picked.append((noisy, exact, meta))
    picked = picked[:24]
    out = []
    for noisy, exact, meta in picked:
        retained = s294.retained_observations(noisy, meta['removed'])
        for method in DD_METHODS:
            start = time.perf_counter()
            try:
                answer = dd_worker(dict(edges=noisy['edges'], observed=retained, method=method,
                                        k_true=meta['K_generating'],
                                        generating_sources=noisy['generating_sources']))
                status = answer['status']
                extra = dict(k=answer.get('k'), crashed=answer.get('jplace_writer_crashed'))
            except BaseException as exc:
                status, extra = 'error', dict(exception=repr(exc)[:200])
            out.append(dict(arm=method, name=meta['name'], level=meta['level'],
                            regime=meta['regime'], K=meta['K_generating'],
                            n_retained=meta['n_retained'], status=status,
                            seconds=round(time.perf_counter() - start, 3), **extra))
        if meta['K_generating'] in EXACT_K:
            keep = set(noisy['observed']) - set(meta['removed'])
            robs = {r: d for r, d in exact['observed'].items() if r in keep}
            start = time.perf_counter()
            answer = exact_worker(dict(edges=exact['edges'], observed=robs, cap=meta['cap']))
            out.append(dict(arm='exact_matched', name=meta['name'], level=meta['level'],
                            regime=meta['regime'], K=meta['K_generating'],
                            n_retained=meta['n_retained'], status=answer['status'],
                            minimum=answer.get('minimum'),
                            seconds=round(time.perf_counter() - start, 3)))
        print(json.dumps(out[-1]), flush=True)
    PILOT.write_text(json.dumps(out, indent=1) + '\n')


def run_dd():
    stream = STREAMS['dd']
    assert not stream.exists(), 'Never overwrite measured runs'
    plan = s294.draw_plan()
    _check(plan)
    jobs = dd_jobs(plan)
    print(f'decodiphy arms: {len(jobs)} jobs on {WORKERS} workers', flush=True)
    run_pool(jobs, dd_worker, DD_WALL, DD_AS, WORKERS, stream, progress=250)


def run_exact_matched():
    stream = STREAMS['exact_matched']
    assert not stream.exists(), 'Never overwrite measured runs'
    plan = s294.draw_plan()
    _check(plan)
    jobs, by_key = exact_jobs(plan)
    print(f'exact_matched arm: {len(jobs)} jobs on {WORKERS} workers', flush=True)
    run_pool(jobs, exact_worker, EXACT_WALL, EXACT_AS, WORKERS, stream,
             post=make_exact_post(by_key), progress=250)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('command', choices=('prepare', 'pilot', 'run-dd', 'run-exact-matched'))
    name = parser.parse_args().command
    {'prepare': prepare, 'pilot': pilot, 'run-dd': run_dd,
     'run-exact-matched': run_exact_matched}[name]()
