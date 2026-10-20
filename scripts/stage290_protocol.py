#!/usr/bin/env python3
"""Stage290 protocol registration: input and code hashes, fixed before any
result is produced.  Writes results/stage290_protocol.json exclusively.
"""
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
sys.dont_write_bytecode = True
from stage290_erasure_core import sha256_path, encode          # noqa: E402

OUT = ROOT / 'results/stage290_protocol.json'
CODE = ['scripts/stage290_erasure_core.py', 'scripts/stage290_survey.py',
        'scripts/stage290_bridge.py', 'scripts/stage290_protocol.py',
        'scripts/stage290_validate_core.py', 'scripts/stage290_analysis.py']
REUSED = ['scripts/stage196_public_erasure_validation.py',
          'scripts/stage25_benchmark.py', 'scripts/stage42_benchmark.py',
          'scripts/stage43_partial_fpt.py', 'scripts/stage286_newick.py',
          'data/stage196_public_erasure_cases.json',
          'results/stage196_public_erasure_validation.json',
          'results/stage196_public_erasure_protocol.json',
          'data/stage287_instances/bank/index.json',
          'results/stage287_instance_bank_protocol.json']


def main():
    index = json.loads((ROOT / 'data/stage287_instances/bank/index.json').read_text())['cases']
    exact = [x for x in index if x['family'].startswith('exact_')]
    bad = [x['path'] for x in exact if sha256_path(ROOT / x['path']) != x['sha256']]
    bank_protocol = json.loads((ROOT / 'results/stage287_instance_bank_protocol.json').read_text())
    panels = [p for p in bank_protocol['reference_panels'] if p.get('pdd_positivity_ok')]
    panels = sorted({p['leaves']: p for p in sorted(panels, key=lambda z: z['file'])}.values(),
                    key=lambda z: -z['leaves'])
    badp = [p['file'] for p in panels if sha256_path(ROOT / p['file']) != p['sha256']]
    protocol = dict(
        stage=290,
        registered_utc=datetime.now(timezone.utc).isoformat(),
        status='REGISTERED_BEFORE_EXECUTION',
        title='Reference-erasure tolerance at scale and the unknown-support '
              'identifiability bridge',
        question=[
            'Does the Thread-F erasure tolerance (cogirth - 1) of the leaf-report '
            'matroid stay at zero across real reference trees, or was the single '
            '1322-leaf measurement an artefact of one tree?',
            'Does the fibre-dimension law |P| - s hold on real trees at scale?',
            'With the support UNKNOWN, how many references can be erased before '
            'more than one minimum-size support explains the observations, and '
            'does the known-support theory predict that threshold?'],
        code_sha256={f: sha256_path(ROOT / f) for f in CODE if (ROOT / f).exists()},
        reused_inputs_sha256={f: sha256_path(ROOT / f) for f in REUSED},
        bank=dict(n_exact_cases=len(exact), hash_mismatches=bad,
                  families=sorted({x['family'] for x in exact})),
        panels=[dict(file=p['file'], leaves=p['leaves'], sha256=p['sha256']) for p in panels],
        panel_hash_mismatches=badp,
        plan=dict(
            survey='Every exactly explainable bank case (support = generating '
                   'support, R = observed references) plus random supports of '
                   'K in {3,5,10,12} on every positivity-conformant frozen panel, '
                   '5 replicates per (panel, K, mode) with mode in {any edge, '
                   'internal edge only}.',
            checks=['modular rank identity rank = s + 1 under two primes',
                    'unreduced dynamic program vs branch-collapsed dynamic program',
                    'exhaustive erasure of observed-leaf subsets (brute force) '
                    'wherever at most 40000 subsets are needed',
                    'independent component recomputation of every DP witness'],
            bridge='exact_biotrees K=2 with |R| <= 130, exact_10k K=3 and '
                   'exact_biotrees K=3 with |R| <= 45; removal regimes random '
                   '(2 seeds), whole_clade (2 seeds), adversarial_nearest; up to '
                   '14 erasures, stopping 3 steps past the first failure.',
            bridge_checks=['every minimum-size support confirmed by the exact '
                           'rational Stage25/42 endpoint polytope oracle',
                           'minimum size cross-checked against the frozen '
                           'solve_partial solver at the full panel',
                           'exhaustive exact enumeration (no prefilter) as ground '
                           'truth on every case with at most 90 edges']),
        limits=dict(workers=16, cores_available=62,
                    brute_force_subset_limit=150000,
                    full_dp_max_vertices=3000,
                    threshold_search='exponential bracketing then bisection, both endpoints evaluated',
                    min_remaining_references=3,
                    prefilter_relative_tolerance=1e-9,
                    candidate_limit=4000000),
        determinism=dict(removal_seeds=dict(random=[11, 23], whole_clade=[11, 23],
                                            adversarial_nearest=[0]),
                         panel_seed_rule='1000*K + 7*replicate + (1 if internal)',
                         prefilter_row_seed=7),
        caveats=[
            'The 1322-leaf public tree of Thread F had negative ANI branch lengths '
            'floored upstream at 1e-4 and its fitted reports are synthetic exact-'
            'model values; the same holds for every exactly explainable bank case. '
            'Statements here are about real tree geometry with synthetic exact '
            'observations, not about biological recovery from measured distances.',
            'Panels carry no mixtures, so panel supports are drawn at random; '
            'because s and the cogirth depend only on topology and support and not '
            'on the edge metric, this is a statement about real tree topology.',
            'The erasure tolerance and the fibre dimension are metric invariant; '
            'the unknown-support threshold is not, since it depends on the observed '
            'values.'],
        environment=dict(python=platform.python_version(),
                         platform=platform.platform(), executable=sys.executable),
        execution_authorized=True)
    with OUT.open('x') as fh:
        json.dump(encode(protocol), fh, indent=2)
        fh.write('\n')
    print(json.dumps(dict(path=str(OUT.relative_to(ROOT)), sha256=sha256_path(OUT),
                          n_exact_cases=len(exact), bank_mismatches=len(bad),
                          panels=len(panels), panel_mismatches=len(badp)), indent=2))


if __name__ == '__main__':
    main()
