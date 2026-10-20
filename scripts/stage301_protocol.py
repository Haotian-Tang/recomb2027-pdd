#!/usr/bin/env python3
"""Stage 301 step 0: register the protocol (inputs + code hashes) before timing.

Nothing is overwritten; the protocol is written once and every later step
asserts the hashes still match.
"""
import json
import os
import sys
import hashlib
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import stage288_common as common            # noqa: E402
import stage301_common as s301              # noqa: E402

OUT = s301.OUT
OUT.mkdir(parents=True, exist_ok=True)
PROTOCOL = OUT / 'stage301_protocol.json'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    assert not PROTOCOL.exists(), 'Never overwrite a registered protocol'
    s298_proto = json.loads((ROOT / 'results/stage298/stage298_protocol.json').read_text())
    cross = [json.loads(l) for l in
             (ROOT / 'results/stage296/stage296_cross_runs.jsonl').read_text().splitlines()]
    probe = [json.loads(l) for l in
             (ROOT / 'results/stage296/stage296_probe_runs.jsonl').read_text().splitlines()]
    code = {}
    for name in ['stage301_common.py', 'stage301_protocol.py', 'stage301_validate.py',
                 'stage301_prune.py', 'stage301_td.py', 'stage301_grid.py',
                 'stage301_analyze.py',
                 'stage296_common.py', 'stage298_common.py', 'stage288_common.py',
                 'stage43_partial_fpt.py', 'stage41_partial_fpt.py',
                 'stage23_two_source_envelope.py', 'stage26_missing_leaf.py']:
        p = ROOT / 'scripts' / name
        code[name] = sha(p) if p.exists() else None
    payload = {
        'stage': 301,
        'created_utc': datetime.now(timezone.utc).isoformat(),
        'question': 'Does the frozen structural localisation transport to the '
                    'eps-tolerant problem, and if so how far does it move the '
                    'reach of Theta_k?',
        'decision_rules': {
            'validation': 'A localised Theta_k must equal the exhaustive Theta_k '
                          'EXACTLY (rational equality) on every case where both '
                          'run.  Any disagreement stops the stage and is reported.',
            'screen_contract': 'A screen may only be used after a written proof '
                               'that it is a necessary condition for theta(S) <= '
                               'eps, and after the validation arm passes.',
            'reading': 'If min theta over survivors <= eps the localised value IS '
                       'Theta_k; otherwise the run certifies Theta_k > eps only.',
            'vacuity': 'Pruning is reported as supports visited / C(|E|, k); a '
                       'ratio of 1.0 is reported as vacuous, not as a success.',
        },
        'screens': {
            'S1': 'NEG(eps) = {v : b(v) < -tau} must lie in N[Sigma]',
            'S2': 'POS(eps) = {v : b(v) >  tau} minus N[Sigma] must be covered '
                  'by the atoms of the support; at most 2k of them',
            'S3': 'h(r) - d_r <= 2 eps at every observed leaf',
            'S4': 'Y >= -eps (|R| - 1)',
            'S5': 'mass-weighted mean projection from leaf s onto the s-r path '
                  'is at most rho = (d(r,s) - (d_r - d_s))/2 + eps',
            'tau': 'tau_eps(v) = sum_{u~v} min(eps / l_uv, 1)',
            'anchors': 'UNRELAXED {v : b(v) > 0}; relaxing the anchor set would '
                       'be unsound (the proof gives b(v*) > 0 strictly)',
        },
        'validation_sources': {
            'stage296_cross_runs': {'n': len(cross),
                                    'sha256': sha(ROOT / 'results/stage296/stage296_cross_runs.jsonl')},
            'stage296_probe_runs': {'n': len(probe),
                                    'sha256': sha(ROOT / 'results/stage296/stage296_probe_runs.jsonl')},
            'stage298_td_summary': {'sha256': sha(ROOT / 'results/stage298/stage298_td_summary.json')},
        },
        'td_cases': s298_proto['td_cases'],
        'code_sha256': code,
        'resources': {'max_workers': 40, 'address_space_bytes': 8 * 1024 ** 3,
                      'note': 'at most 40 of the 62 cores, per the stage brief'},
    }
    PROTOCOL.write_text(json.dumps(payload, indent=1) + '\n')
    print(json.dumps({k: v for k, v in payload.items() if k != 'td_cases'}, indent=1)[:1500])


if __name__ == '__main__':
    main()
