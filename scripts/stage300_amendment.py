#!/usr/bin/env python3
"""Stage 300 protocol amendment: what was added after registration, and why.

The registered protocol covers the inputs, the harvest, the exact evaluation and
the verification plan.  Three things were added while the stage ran, each in
response to something the measurements showed; they are hashed here rather than
folded silently into the original file.
"""
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import stage300_common as s300                   # noqa: E402

ADDED = {
    'scripts/stage300_resolve.py':
        'fairness arm: re-solve DecoDiPhy\'s own QP (solve_osqp_once) at the support it '
        'returned, because the recomputed SSR of the returned configuration disagreed '
        'with the loss the tool reports for it',
    'scripts/stage300_lossgap.py':
        'quantifies that disagreement over every evaluated run and records the '
        'endpoint-substitution trigger condition',
    'scripts/stage300_losscase.py':
        'the worst single case, re-solved at its printed support and at every '
        'single-anchor neighbour of it',
    'scripts/stage300_prune.py':
        'fair-panel arm: measures the reference panel DecoDiPhy keeps after '
        'preprocess_input and recomputes BOTH Theta_1 (stage 298\'s exact k=1 routine) '
        'and R_dd on that panel',
    'scripts/stage300_routee.py':
        'verification route E: Lambda_1 against the tool\'s own k = 1 ladder loss',
    'scripts/stage300_tdprobe.py':
        'cost decomposition of the tool\'s input stage on one 10,575-leaf stool sample',
    'scripts/stage300_analyze.py': 'aggregation of all records into stage300_analysis.json',
    'scripts/stage300_report.py': 'table builders read by the document generator',
    'scripts/stage300_document.py': 'assembles stage300_own_fit_comparison.md',
    'scripts/stage300_table.py': 'per-case three-way CSV',
    'scripts/stage300_censor.py':
        'registers the 7,200 s per-sample wall for the travellers collection and '
        'records the six censored runs with the wall-clock and resident set size '
        'measured at the moment they were stopped',
    'scripts/stage300_amendment.py': 'this file',
}


def main():
    amendment = dict(
        stage=300, step='protocol_amendment',
        registered_utc=time.strftime('%Y-%m-%dT%H:%M:%S+00:00', time.gmtime()),
        relation_to_protocol=('additions to results/stage300/stage300_protocol.json; '
                              'nothing registered there was changed, and the code hashes '
                              'registered there are unchanged on disk'),
        code_added={path: dict(sha256=s300.sha(ROOT / path), why=why)
                    for path, why in sorted(ADDED.items())},
        code_registered_unchanged={path: s300.sha(ROOT / path) for path in
                                   json.loads((s300.OUT / 'stage300_protocol.json')
                                              .read_text())['code']},
        run_configuration_actually_used=dict(
            krepp_placement='60 cases x {dd_default, dd_true_k}, 20 workers, wall 1800 s, '
                            'RLIMIT_AS 8 GB; 120/120 runs complete',
            travelers_pilot='1 sample (td_sample_101_published), 1 worker, wall 14400 s, '
                            'RLIMIT_AS 60 GB',
            travelers_batch='the other 28 samples, 5 workers, wall 21600 s, '
                            'RLIMIT_AS 30 GB',
            cores='at most 26 of 62 in flight at any time (20 placement workers, then '
                  '5 travellers workers + 1 pilot + 1 probe), within the 40-core budget',
            krepp_ani='no new runs: harvested from stages 289/293/294/297'),
        travelers_wall=dict(
            seconds=7200.0,
            registered_in='results/stage300/stage300_travelers_cost.json',
            why='the registered plan left the travellers wall to be sized from the '
                'pilot; the pilot did not return, so a wall was set from the budget '
                'instead: 29 samples x 2 h / 5 memory-safe concurrent workers is '
                'about 11.6 h, the largest per-sample wall that fits a single day '
                'of this stage\'s 40-core allocation',
            outcome='6 of 29 samples attempted, 0 completed, 12.9 core-hours '
                    'consumed, 9.9-10.8 GiB resident per run at the stop'),
        deviations_from_the_registered_plan=[
            'The original plan said the travellers workers would be "sized from the '
            'pilot\'s peak RSS"; the pilot had not finished when the batch was launched, '
            'so 5 workers x 30 GB RLIMIT_AS was chosen from its RSS trajectory '
            '(5.9 -> 7.5 GB over 50 minutes) and the measured machine headroom.',
            'The analysis, report and table scripts were written after registration and '
            'are hashed in this amendment, as is the house pattern.',
            'The travellers collection is reported as a measured cost, not as a fit '
            'comparison: 6 of 29 samples attempted, 0 completed within the 7,200 s '
            'wall, 23 not attempted.  Theta_1 (imported from stage 298) and Lambda_1 '
            '(computed here) are reported for all 29.'])
    path = s300.OUT / 'stage300_protocol_amendment.json'
    path.write_text(json.dumps(amendment, indent=1))
    print('wrote', path)
    print(json.dumps({k: v['sha256'][:12] for k, v in amendment['code_added'].items()},
                     indent=1))


if __name__ == '__main__':
    main()
