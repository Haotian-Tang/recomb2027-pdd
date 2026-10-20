#!/usr/bin/env python3
"""Stage297 protocol registration.  Run BEFORE any timed batch; refuses to
overwrite.  An addendum file is written for every batch whose specification is
fixed later (the main K = 10 sweep is sized from the pilot).
"""
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
sys.dont_write_bytecode = True

import stage297_common as common                                          # noqa: E402
from stage297_common import OUTDIR, sha                                   # noqa: E402

PROTOCOL = OUTDIR / 'stage297_protocol.json'

NEW_CODE = ('stage297_reduced_all_minimisers.py', 'stage297_common.py',
            'stage297_validate.py', 'stage297_bridge.py', 'stage297_seeds.py',
            'stage297_protocol.py')
IMPORTED_CODE = ('stage292_all_minimisers.py', 'stage292_bridge.py',
                 'stage292_orders.py', 'stage292_harness.py',
                 'stage292_validate.py', 'stage292_analysis.py',
                 'stage295_reduced_minimum.py', 'stage295_common.py',
                 'stage290_erasure_core.py', 'stage293_abundance.py',
                 'stage293_fibre.py', 'stage289_decodiphy_bridge.py',
                 'stage289_headtohead_analyze.py', 'stage288_common.py',
                 'stage288_incompleteness.py', 'stage43_partial_fpt.py',
                 'stage41_partial_fpt.py', 'stage27_general_fpt.py',
                 'stage26_missing_leaf.py', 'stage25_deficit_extension.py',
                 'stage25_benchmark.py', 'stage23_two_source_envelope.py')
INPUTS = ('data/stage287_instances/bank/index.json',
          'results/stage292_validate_replay.jsonl',
          'results/stage292_validate_fresh.jsonl',
          'results/stage292_bridge_K2_3_5_7.jsonl',
          'results/stage292_bridge_K7_resume.jsonl',
          'results/stage292_bridge_K10_12.jsonl',
          'results/stage292_k10_ceiling.txt',
          'results/stage292_analysis.json')


def job_digest(names):
    import hashlib
    return hashlib.sha256(json.dumps(sorted(names)).encode()).hexdigest()


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    assert not PROTOCOL.exists(), 'Never overwrite a registered protocol'
    import stage297_validate as validate
    import stage297_bridge as bridge
    plans = dict(
        validate_A=[j['entry']['name'] for j in validate.plan_A()],
        validate_B=['fresh_seed_%d' % j['seed'] for j in validate.plan_B()],
        validate_C=[j['entry']['name'] for j in validate.plan_C()],
        bridge_pilot_k10_baseline=[j['entry']['name'] for j in
                                   bridge.build_jobs('pilot_k10_baseline',
                                                     **bridge.PLANS['pilot_k10_baseline'])],
        bridge_pilot_k10_threshold=[j['entry']['name'] for j in
                                    bridge.build_jobs('pilot_k10_threshold',
                                                      **bridge.PLANS['pilot_k10_threshold'])],
        bridge_k12_baseline=[j['entry']['name'] for j in
                             bridge.build_jobs('k12_baseline',
                                               **bridge.PLANS['k12_baseline'])],
    )
    dd_plan = None
    try:
        import stage297_seeds as seeds
        dd_plan = [j['entry']['name'] for j in seeds.plan(1)]
    except Exception as exc:                 # decodiphy venv only
        dd_plan = ['IMPORT_DEFERRED:%r' % exc]
    protocol = dict(
        stage=297,
        registered_utc=datetime.now(timezone.utc).isoformat(),
        status='REGISTERED_BEFORE_ANY_TIMED_RUN_OF_THE_BATCHES_IT_COVERS',
        budget=dict(cores_allowed=24, validate_workers=16, bridge_workers=8,
                    dd_workers=8, note='at most 24 of the 62 server cores at any time'),
        item1=dict(
            question=('does the unknown-support erasure threshold stay flat in '
                      'reference-set size at K = 10, the regime Stage292 could not '
                      'reach because the FROZEN minimum pass its all-minimisers '
                      'routine calls did not complete'),
            new_module='scripts/stage297_reduced_all_minimisers.py',
            modes=dict(
                frozen='pass 1 = frozen stage43_partial_fpt.solve_partial (Stage292 reference arm)',
                reduced_two_pass='pass 1 = stage295_reduced_minimum.solve_partial_reduced',
                combined=('one pass: reduced face enumeration computes the minimum and '
                          'the complete minimiser set together')),
            already_present_in_stage292=('the forced-positivity face reduction is ALREADY '
                                         'in stage292_all_minimisers.backbone_patterns '
                                         '(full_mask switch, 500/500 agreement on its fresh '
                                         'track); what Stage295 adds is the reduced MINIMUM '
                                         'pass, which is where the K >= 10 ceiling sat'),
            validation=dict(
                track_A=('replay of every probe point in stage292_validate_replay.jsonl '
                         '(130 cases, 5339 points, K in {2,3}); removal orders regenerated '
                         'with the frozen removal_order; complete minimiser SETS compared '
                         'between the Stage292 routine and all three Stage297 modes'),
                track_B=('the 500 fresh instances of stage292_validate_fresh.jsonl '
                         'regenerated from stage292_validate.random_instance, compared '
                         'against the independent exhaustive C(|E|,k) enumeration '
                         'stage292_validate.brute_all and between reduced and full 2^|A| '
                         'face enumeration'),
                track_C=('full-panel baseline of all 676 Stage292 bridge cases '
                         '(K in {2,3,5,7}): Stage292 routine vs Stage297 combined vs the '
                         'value Stage292 recorded, plus, on every 10th case, equality of '
                         'the Stage297 full-mask face count with the face count of the '
                         'Stage292 enumeration run with full_mask=True'),
                stop_rule='any disagreement is reported, never patched'),
            measurement=dict(
                regimes='stage292_bridge.REGIMES, imported verbatim',
                min_remaining='stage292_bridge.MIN_REMAINING = 3, imported verbatim',
                threshold_search=('stage292_bridge.find_threshold, imported verbatim: '
                                  'exponential bracketing then bisection with both '
                                  'bracket endpoints evaluated'),
                predicate=('minimum == K, exactly one minimiser, and that minimiser is '
                           'the true support (stage290 ok)'),
                indecisive_probe=("threshold_status='unknown_probe:<status>'; never read "
                                  'as a failure'),
                retained='every probe with its counters, including timeouts and deaths'),
            batches=dict(
                pilot_k10_baseline=('one case per distinct |R| at K = 10 (20 cases, |R| 20-353), '
                                    'full-panel baseline only, probe budget 1200 s, wall 2400 s, '
                                    'plus a reduced_two_pass cross-check per case'),
                pilot_k10_threshold='the 6 smallest K = 10 cases, full threshold search',
                k12_baseline='the single exact_16k K = 12 case (|R| = 15941), baseline only',
                main_k10='SPECIFICATION FIXED BY ADDENDUM once the pilot has calibrated cost')),
        item2=dict(
            question=('how much of the structural error attributed to DecoDiPhy is '
                      'seed-dependent: spread in inferred k, in the returned support and '
                      'in W1 to truth across seeds, per case'),
            tool='DecoDiPhy 1.5.0 at data/stage289_decodiphy/{repo,venv}',
            invocation=('method dd_default, documented defaults, k by the tool own '
                        'stopping rules; --method hill --optimizer osqp are the defaults '
                        'and are untouched; only -s differs between arms'),
            arms=dict(
                default_1142=('stage293_abundance.dd_worker imported UNMODIFIED, no -s '
                              'flag, so the tool takes its documented default seed 1142 '
                              '(cli.py:39) -- byte-identical to stages 289/293/294'),
                seeded=('stage297_seeds.dd_worker_seeded, a verbatim copy whose only '
                        'change is "-s <seed>"; the unified diff is recorded in '
                        'results/stage297/stage297_dd_worker_diff.json'),
                seeded_1142=('the default seed repeated through the copy: any difference '
                             'from default_1142 is process-level nondeterminism at a '
                             'fixed seed, not seed sensitivity')),
            seeds=[1142, 1, 7, 101, 2024, 31337],
            panel='full panel only; no reference removal in this item',
            cases=('stratified over |R| within every (family, K) cell of exact_10k, '
                   'exact_biotrees and krepp_estimated_10k'),
            prior=dict(
                stage289='single seed 1142 only; seed sensitivity recorded as unmeasured',
                stage294=('stage289_full_panel_agreement counts: compared 575, k_agrees 567, '
                          'support_agrees 471 -- the returned support is not reproducible '
                          'across harness process details at a fixed seed')),
            recovered_outputs=('DecoDiPhy 1.5.0 can raise in its jplace writer after '
                               'all_rounds.json is already committed; stage 293 recovers the '
                               'tool own placements from that file and flags '
                               'jplace_writer_crashed.  Those runs count as results and are '
                               'reported separately'),
            w1=('stage289_headtohead_analyze.wasserstein with point_distances, computed in '
                'the analysis step from the placements the workers record')),
        plans={k: dict(n=len(v), digest=job_digest(v)) for k, v in plans.items()},
        plan_names=plans,
        dd_plan=dict(n=len(dd_plan), digest=job_digest(dd_plan), names=dd_plan),
        code_sha256={name: sha('scripts/' + name) for name in NEW_CODE},
        imported_code_sha256={name: sha('scripts/' + name) for name in IMPORTED_CODE},
        input_sha256={name: sha(name) for name in INPUTS},
        environment=dict(
            python=sys.version.split()[0],
            project_venv=subprocess.run(
                ['.venv/bin/python', '-c', 'import numpy,sys;print(sys.version.split()[0],numpy.__version__)'],
                cwd=str(ROOT), capture_output=True, text=True).stdout.strip(),
            decodiphy_venv=subprocess.run(
                ['data/stage289_decodiphy/venv/bin/python', '-c',
                 'import numpy,scipy,decodiphy,sys;print(sys.version.split()[0],numpy.__version__,scipy.__version__,getattr(decodiphy,"__version__","1.5.0"))'],
                cwd=str(ROOT), capture_output=True, text=True).stdout.strip(),
            threads='OMP/OPENBLAS/MKL/VECLIB/NUMEXPR_NUM_THREADS=1 in every worker'),
        safety=('new files only, all under results/stage297 and scripts/stage297_*; '
                'nothing pre-existing is modified, renamed or deleted'),
    )
    PROTOCOL.write_text(json.dumps(protocol, indent=1, sort_keys=False))
    print('registered', PROTOCOL)
    print('plans:', {k: v['n'] for k, v in protocol['plans'].items()},
          'dd', protocol['dd_plan']['n'])


if __name__ == '__main__':
    main()
