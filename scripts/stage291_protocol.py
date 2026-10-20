#!/usr/bin/env python3
"""Register the Stage 291 protocol BEFORE any measurement is taken.

Writes results/stage291_protocol.json and refuses to overwrite it.
"""
import hashlib, json, sys
from datetime import datetime, timezone
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from stage291_common import sha, dependency_closure
import stage291_screen_run as runner

OUT = ROOT / 'results/stage291_protocol_v2.json'
SUPERSEDES = ('results/stage291_protocol.json registered the same design with an eta grid that '
              'stopped at 1e-1, below the ANI error range of the krepp family. No measurement '
              'was recorded under it: the only run was an 8-case dry run left in place as '
              'results/stage291_full_pilot_runs.jsonl. v2 extends the grid to 3e-1 and 1, adds '
              'tau deciles and the two-sided exclusion count, and is the protocol every '
              'reported number is produced under.')

FROZEN = [
    'data/stage287_instances/bank/index.json',
    'results/stage287_instance_bank_protocol.json',
    'scripts/stage44_certificate.py',
    'scripts/stage43_theory_robust_mass.py',
    'scripts/stage288_common.py',
    'scripts/stage291_common.py',
    'scripts/stage291_screen_run.py',
]

TRACKS = {
    'full': {'families': ['exact_10k', 'exact_biotrees', 'published_true_10k',
                          'published_biotrees', 'krepp_estimated_10k', 'exact_16k'],
             'panels': 'full observed panel only'},
    'removal': {'families': ['exact_10k', 'exact_biotrees', 'published_true_10k',
                             'published_biotrees', 'krepp_estimated_10k'],
                'panels': 'full panel plus Stage 288 removal regimes'},
    'real': {'families': ['published_travelers_diarrhea'],
             'panels': 'full observed panel only; no ground truth exists'},
}

CLAIM_BOUNDARY = {
    'what_a_retained_edge_means': ('Nothing. Retention by the screen is the absence of an '
        'exclusion certificate, not evidence that a source is present in that edge.'),
    'what_an_excluded_edge_means': ('Under H1-H4 below, no source of the generating mixture '
        'lies strictly inside that original edge. The certificate is conditional on the '
        'supplied bounds; it is not a statistical test and carries no error rate.'),
    'H1': 'The reference topology and every edge length are known exactly.',
    'H2': ('Every source lies strictly inside an original edge at fractional position x with '
           'min(x,1-x) >= rho and carries weight p >= alpha. alpha and rho are SUPPLIED.'),
    'H3': 'The common mean pendant satisfies y >= floor > 0. floor is SUPPLIED.',
    'H4': ('The screened vector dominates the true model vector coordinatewise: '
           'observed_r + eta >= d_r for every observed reference r. eta is SUPPLIED.'),
    'oracle_setting': ('alpha, rho and floor taken from the published generating mixture. This '
        'is the most favourable admissible choice and is reported as an upper bound on the '
        'screening power attainable with correctly supplied bounds, not as a usable procedure.'),
    'eta_plus': ('The smallest eta satisfying H4 for a given case and panel, computed from the '
        'published truth. Certificates at eta >= eta_plus are sound; the behaviour at smaller '
        'eta measures what an under-stated error bound costs.'),
    'reduced_panels': ('With references removed the screen runs on the hull of the surviving '
        'panel. For a hull edge made of original edges e_1..e_m the per-original-edge offset '
        'bounds prefix_k + rho*L_k and suffix_k + rho*L_k are used; this reduces to the frozen '
        'Stage 44 threshold alpha*rho*L exactly when m = 1. Original edges off the hull carry '
        'no certificate at any eta and are counted as not excluded.'),
    'no_recovery_claim': ('No track estimates a source support, a source count or a weight. '
        'Nothing here is a recovery-accuracy measurement.'),
    'travelers_diarrhea': ('No ground truth exists for these 29 real samples. Only the number '
        'of certifiably excluded reference edges under explicitly stated supplied bounds is '
        'reported. No accuracy, recovery or biological identification claim is made or implied.'),
    'stage130': ('The Stage 130 metric-uncertainty support certificate was inspected. '
        'results/stage130_prior_metric_uncertainty_candidate.json is an analytic candidate tied '
        'to the Stage 125/128 constructed worlds and the only covector implementation in the '
        'repository is scripts/stage126_structure_joint_validation.py:covector, which is '
        'likewise world-specific. No general implementation exists, so nothing from Stage 130 '
        'is reused here. Its zero-sum covector family also excludes the Stage 44 pair covector, '
        'whose coefficients sum to one, which is exactly why the pair screen needs a pendant '
        'floor while the Stage 130 constants do not.'),
}


def main():
    assert not OUT.exists(), 'protocol already registered; refusing to overwrite'
    record = dict(
        stage=291, track='noise_and_completeness', protocol_version=2,
        supersedes=SUPERSEDES,
        created_utc=datetime.now(timezone.utc).isoformat(),
        purpose=('Measure certified edge-exclusion power on estimated (not exactly explainable) '
                 'distance vectors, check the certificates against published ground truth, and '
                 'join the error axis to the Stage 288 reference-incompleteness axis.'),
        frozen_inputs={name: sha(ROOT / name) for name in FROZEN},
        dependency_closure={name: sha(ROOT / 'scripts' / name)
                            for name in dependency_closure('stage291_screen_run.py')},
        tracks=TRACKS,
        eta_grid=[str(e) for e in runner.ETA_GRID],
        removal={'levels': [str(x) for x in runner.LEVELS], 'regimes': list(runner.REGIMES),
                 'replicates': runner.REPLICATES, 'seed': runner.SEED,
                 'definitions': 'scripts/stage288_common.py, unchanged'},
        hypothesis_settings={
            'oracle': 'alpha=min p_i, rho=min min(x_i,1-x_i), floor=true common mean pendant',
            'assumed_mid': 'alpha=1/(2 cap), rho=1/10, floor=y/2',
            'assumed_weak': 'alpha=1/(4 cap), rho=1/100, floor=y/10',
            'data_scaled': 'alpha=1/(2 cap), rho=1/10, floor=(min observed)/2',
        },
        limits={'wall_seconds_per_case': runner.WALL, 'address_space_bytes': runner.ADDRESS_SPACE,
                'workers': 28, 'note': 'every run retained, including timeouts and errors'},
        claim_boundary=CLAIM_BOUNDARY)
    OUT.write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps({'written': str(OUT), 'sha256': sha(OUT),
                      'frozen_inputs': len(record['frozen_inputs']),
                      'closure': len(record['dependency_closure'])}, indent=2))


if __name__ == '__main__':
    main()
