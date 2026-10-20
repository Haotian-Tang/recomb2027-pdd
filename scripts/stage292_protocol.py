#!/usr/bin/env python3
"""Register the Stage 292 protocol BEFORE any measurement is recorded.

Writes results/stage292_protocol.json and refuses to overwrite it.
"""
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
sys.dont_write_bytecode = True
from stage292_harness import sha                                        # noqa: E402

OUT = ROOT / 'results/stage292_protocol.json'

NEW_CODE = ['scripts/stage292_all_minimisers.py', 'scripts/stage292_harness.py',
            'scripts/stage292_orders.py', 'scripts/stage292_validate.py',
            'scripts/stage292_bridge.py', 'scripts/stage292_panels.py',
            'scripts/stage292_pilot.py']

FROZEN = ['scripts/stage23_two_source_envelope.py',
          'scripts/stage25_deficit_extension.py',
          'scripts/stage25_benchmark.py',
          'scripts/stage26_missing_leaf.py',
          'scripts/stage27_general_fpt.py',
          'scripts/stage41_partial_fpt.py',
          'scripts/stage42_benchmark.py',
          'scripts/stage43_partial_fpt.py',
          'scripts/stage286_newick.py',
          'scripts/stage290_erasure_core.py',
          'scripts/stage290_bridge.py',
          'scripts/stage290_survey.py',
          'results/stage290_bridge.jsonl',
          'results/stage290_protocol.json',
          'data/stage287_instances/bank/index.json',
          'results/stage287_instance_bank_protocol.json']

PANELS = ['data/stage287_instances/krepp_shared/results__expt-emp__emp-pp-WoLv2__placement_tree_extended.nwk',
          'data/stage287_instances/krepp_shared/misc_data__backbone_blen_conserved-WoLv2.nwk']


def main():
    if OUT.exists():
        raise SystemExit('protocol already registered: ' + str(OUT))
    index = json.loads((ROOT / 'data/stage287_instances/bank/index.json').read_text())
    mismatch = []
    exact = [x for x in index['cases'] if x['family'].startswith('exact')]
    for x in exact:
        if sha(x['path']) != x['sha256']:
            mismatch.append(x['name'])
    doc = {
        'stage': 292,
        'registered_utc': datetime.now(timezone.utc).isoformat(),
        'status': 'REGISTERED_BEFORE_EXECUTION',
        'title': 'Exact all-minimisers enumeration and the unknown-support erasure '
                 'threshold at K >= 5 and at panel scale',
        'question': [
            'Can the complete set of minimum-size strict source supports be decided '
            'exactly without the C(|E|,K) brute force Stage290 used, by enumerating '
            'only the O(K)-local deficit supports, zero faces, minimum covers and '
            'lifts?',
            'Does the resulting routine agree with brute force wherever brute force '
            'is affordable?',
            'Is the Stage290 finding -- that the adversarial (nearest-to-source) '
            'erasure threshold does not grow with panel size while the random one '
            'does -- still true at K = 5 and above, and on reference panels spanning '
            'more than two orders of magnitude in size?'],
        'new_code_sha256': {p: sha(p) for p in NEW_CODE},
        'reused_inputs_sha256': {p: sha(p) for p in FROZEN},
        'panel_sha256': {p: sha(p) for p in PANELS},
        'bank': {'n_exact_cases': len(exact), 'hash_mismatches': mismatch,
                 'families': sorted({x['family'] for x in exact})},
        'method': {
            'all_minimisers': (
                'Backbone projection (frozen observed_backbone); leaf envelope and '
                'baseline mass (frozen two_pass_envelope / inverse); admissible '
                'deficit supports from the FROZEN generator deficit_sets under the '
                'weighted deficit bound sum(deg-2) <= K-2 and 3 <= deg <= K; zero '
                'faces of the affine mass map over A = S u N(S), reduced by the '
                'forced-positivity of N(S)\\S vertices with non-negative baseline; '
                'exact strict feasibility by the frozen Stage43 common-slack vertex '
                'routine; ALL minimum atom covers of each feasible pattern; and all '
                'original-edge lifts of every atom, each decided by the frozen '
                'Stage25 exact rational endpoint oracle fixed_support_polytope.'),
            'exhaustiveness_rests_on': [
                'the weighted deficit bound and positive-anchor growth implemented by '
                'the frozen deficit_sets generator, imported and not re-derived here',
                'the exactness of the frozen Stage25 endpoint polytope oracle',
                'the lift candidate lists being supersets of the true lifts'],
            'minimum_cross_check': ('every call also runs the frozen solve_partial and '
                                    'records its minimum; a mismatch is reported as a '
                                    'disagreement and never silently repaired'),
            'threshold_convention': (
                "Stage290's: t_unknown is the largest e along a fixed removal order "
                'for which the set of all minimum-size supports is exactly the single '
                'true support; found by exponential bracketing then bisection with '
                'both bracket endpoints explicitly evaluated; MIN_REMAINING = 3 '
                'references are never erased; a probe with an indecisive status '
                'aborts the regime as unknown_probe rather than counting as failure.'),
            'regimes': ['random seed 11', 'random seed 23', 'whole_clade seed 11',
                        'whole_clade seed 23', 'adversarial_nearest deterministic'],
        },
        'tracks': {
            'validate_A': ('every probe point recorded in results/stage290_bridge.jsonl '
                           '(130 cases x 5 regimes x the probes of each bisection) '
                           'plus e = 0, recomputed with the new routine and re-run '
                           'with the Stage290 brute force; cases with <= 70 edges also '
                           'get the prefilter-free exhaustive brute force'),
            'validate_B': ('500 freshly generated small instances (7-12 leaves, K in '
                           '1..4, random hidden-leaf subsets) compared support-set for '
                           'support-set against the complete C(|E|,k) enumeration, and '
                           'the reduced zero-face enumeration compared against the full '
                           '2^|A| one'),
            'bridge': ('all 893 exactly-explainable bank cases, K in {2,3,5,7,10,12}, '
                       'five regimes each'),
            'panels': ('planted mixtures on the frozen 47,151-leaf WoL v2 placement '
                       'panel sub-sampled to 64...47151 leaves and on the 15,953-leaf '
                       'WoL v2 backbone sub-sampled to 100/1000/15953, K in {3,5}, '
                       'three seeds, adversarial and random regimes'),
        },
        'claim_boundary': {
            'planted_instances': ('the panel-ladder mixtures are PLANTED by the rule in '
                                  'scripts/stage292_panels.py, not published mixtures; '
                                  'they vary panel size at a fixed generative process '
                                  'and carry no biological claim'),
            'resource_stops': ('a non-complete routine status is recorded as such and '
                               'never read as a uniqueness verdict'),
            'what_t_unknown_is_not': ('t_unknown is a property of one removal ORDER, not '
                                      'a worst case over all erasure sets; the worst case '
                                      'over sets is the known-support cogirth quantity'),
        },
        'resources': {'workers': 28, 'memory_bytes_per_child': 10 << 30,
                      'per_case_wall_seconds': {'K<=5': 1800, 'K<=7': 5400,
                                                'K>=10': 10800,
                                                'panels<=2048': 600,
                                                'panels>2048': 7200}},
        'safety': 'New files only under stage 292; no pre-existing file is modified.',
    }
    OUT.write_text(json.dumps(doc, indent=2) + '\n')
    print(json.dumps({k: v for k, v in doc.items()
                      if k not in ('new_code_sha256', 'reused_inputs_sha256',
                                   'panel_sha256')}, indent=2)[:2500])
    print('bank hash mismatches:', mismatch)
    print('written', OUT)


if __name__ == '__main__':
    main()

