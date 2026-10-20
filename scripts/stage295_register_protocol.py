#!/usr/bin/env python3
"""Stage295 protocol registration: diagnose the K>=10 ceiling of the frozen
minimum pass `stage43_partial_fpt.solve_partial`.

Registered BEFORE any measurement run, following the house pattern: input and
code hashes, the full measurement plan, the resource budget, and the decision
rules that will be applied to the result.
"""
import json, os, subprocess, sys, time
from datetime import datetime, timezone

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
ROOT = os.path.dirname(HERE)
from stage295_common import sha

FROZEN = ['scripts/stage43_partial_fpt.py', 'scripts/stage41_partial_fpt.py',
          'scripts/stage27_general_fpt.py', 'scripts/stage26_missing_leaf.py',
          'scripts/stage25_deficit_extension.py',
          'scripts/stage23_two_source_envelope.py',
          'scripts/stage292_all_minimisers.py']
NEW = ['scripts/stage295_common.py', 'scripts/stage295_static_census.py',
       'scripts/stage295_profile.py']
INPUTS = ['data/stage287_instances/bank/index.json',
          'stage295_tmp/census_plan.json', 'stage295_tmp/profile_plan.json']


def hashes(paths):
    out = {}
    for p in paths:
        full = os.path.join(ROOT, p)
        out[p] = sha(full) if os.path.exists(full) else None
    return out


def main():
    plan = json.load(open(os.path.join(ROOT, 'stage295_tmp/profile_plan.json')))
    census = json.load(open(os.path.join(ROOT, 'stage295_tmp/census_plan.json')))
    protocol = {
        'stage': 295,
        'title': 'Diagnosis of the K>=10 ceiling in the frozen exact minimum pass',
        'registered_utc': datetime.now(timezone.utc).isoformat(),
        'question': (
            'Stage292 established that K>=10 is unreachable because the frozen '
            'minimum pass stage43_partial_fpt.solve_partial returns status '
            '"unknown" after a 3600 s budget even on the smallest K=10 bank '
            'case (|R|=20, 37 backbone edges). This stage determines WHICH '
            'phase of that search explodes, whether the realised growth in K '
            'matches the paper\'s 2^O(K) poly(B) shape, and whether the '
            'dominant cost is removable without touching the frozen solver.'),
        'pipeline_phases': {
            'P1_deficit_supports': 'stage41_partial_fpt.deficit_sets; counters '
                                   'growth_states, deficit_sets',
            'P2_zero_faces': 'for mask in range(1 << len(A)) with A = S u N(S); '
                             'counter mass_faces; theory bound 2^(4K-8)',
            'P3_cover_filter': 'stage41_partial_fpt.mixed_cover per face; '
                               'counter cover_feasible_faces',
            'P4_active_bases': 'stage43_partial_fpt.strict_vertex basis '
                               'enumeration C(m, d), d = |S|+1; counters '
                               'vertex_bases, local_lp_calls, interval_lp_calls, '
                               'enumerated_lp_calls, zero_dimensional_lp_calls, '
                               'eliminated_equality_rank'},
        'instrumentation': (
            'Entirely external. The frozen modules are imported, never '
            'modified: (a) the budget knobs seconds/max_states/max_bases/'
            'max_faces that solve_partial already exposes are used to take '
            'geometric work prefixes, so the cumulative time-vs-work curve '
            'separates per-face overhead from per-basis LP cost; (b) the '
            'statistics dict solve_partial already returns supplies every '
            'counter at stop; (c) cProfile around the call supplies '
            'function-level tottime shares as an independent cross-check; '
            '(d) a static census replicates only the preamble and the '
            'control flow that fixes the face count, with no LP work, giving '
            'exact phase-1 and phase-2 counts.'),
        'measurement_plan': {
            'static_census_rows': len(census),
            'profile_rows': len(plan),
            'tags': sorted({r['tag'] for r in plan}),
            'k_ladder_note': (
                'The Stage287 bank contains caps {2,3,5,7,10,12} only; there '
                'are no K=8 or K=9 bank families. K=8 and K=9 are therefore '
                'realised by sweeping the solver cap knob on a fixed K=10 '
                'case, which isolates cap as the sole variable (maxweight = '
                'cap-2, |A| <= 4*cap-8, d = |S|+1 all depend on cap alone). '
                'This is declared as the one scope substitution of the stage.'),
        },
        'resource_budget': {'workers': 16, 'cores_available': 62,
                            'per_run_mem_gb': 6,
                            'per_run_rlimit': 'RLIMIT_AS enforced in a forked child',
                            'per_run_wall_backstop_s': 'SIGALRM above the solver seconds budget',
                            'retention': 'every run retained, including timeouts and failures'},
        'decision_rules': {
            'intrinsic_vs_artefact': (
                'The ceiling is called INTRINSIC if the realised per-K growth '
                'factor of the dominant counter is within the shape the paper '
                'bound allows (face count <= 2^(4K-8), supports and bases '
                '2^O(K)) AND no reformulation of the dominant phase reduces '
                'the realised exponent. It is called an IMPLEMENTATION '
                'ARTEFACT if the dominant cost admits a reformulation that '
                'provably visits strictly fewer objects without changing the '
                'answer.'),
            'prototype_gate': (
                'A prototype is built only if the profile identifies a single '
                'dominant and clearly removable cost. Any prototype lives in '
                'a NEW module and must agree with the frozen solver on the '
                'minimum source count on EVERY case where both complete. A '
                'single disagreement is a stop-and-report, not something to '
                'patch around.'),
            'reporting': (
                'Counts and time shares are reported, not impressions. Any '
                'number resting on a truncated or budget-limited run is '
                'labelled as such.'),
        },
        'code_sha256_frozen': hashes(FROZEN),
        'code_sha256_new': hashes(NEW),
        'input_sha256': hashes(INPUTS),
        'python': sys.version.split()[0],
    }
    out = os.path.join(ROOT, 'results/stage295_ceiling_protocol.json')
    json.dump(protocol, open(out, 'w'), indent=1)
    print('registered', out)
    print('frozen hashes:', json.dumps(protocol['code_sha256_frozen'], indent=1)[:400])


if __name__ == '__main__':
    main()
