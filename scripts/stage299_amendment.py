#!/usr/bin/env python3
"""Stage 299 amendment: every code change made after the protocol was registered.

The protocol hashed the stage-299 scripts at registration time.  Three of them
changed afterwards.  None of them decides a theta - every number in this stage
is produced by frozen stage-296/298 code - but the record must say so
explicitly, with the final hashes.
"""
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import stage288_common as common               # noqa: E402
import stage299_common as s299                 # noqa: E402

OUT = s299.OUT


def main():
    path = OUT / 'stage299_protocol_amendment.json'
    assert not path.exists(), 'Never overwrite frozen evidence'
    protocol = json.loads((OUT / 'stage299_protocol.json').read_text())
    final = {}
    for name in sorted(protocol['code_hashes']):
        p = ROOT / 'scripts' / name
        final[name] = common.sha(p) if p.exists() else 'ABSENT'
    final['stage299_amendment.py'] = 'self'
    changed = sorted(n for n in final
                     if final[n] != protocol['code_hashes'].get(n) and final[n] != 'self')
    out = dict(
        stage=299, step='protocol_amendment',
        registered_utc=datetime.now(timezone.utc).isoformat(),
        purpose='Record the post-registration code changes and the final hashes.',
        changes=[
            dict(script='stage299_validate.py', versions=3,
                 what=('Two successive defects in the VALIDATOR\'S OWN reference computation for '
                       'a point pair lying on the SAME edge: v1 inserted a single subdivision '
                       'vertex for both points, v2 kept the two handles in a dict keyed by edge '
                       'index, which collides.  Both made the reference read 0.  '
                       '`stage299_common.point_distance` and `edge_hops`, the code under test, '
                       'were never changed; v3 agrees with the subdivision reference on 126/126 '
                       'pairs including every same-edge pair.'),
                 evidence_retained=['stage299_validation.json (v1, failed)',
                                    'stage299_validation_v2.json (v2, failed)',
                                    'stage299_validation_v3.json (v3, passed)'],
                 affects_measurements=False),
            dict(script='stage299_small.py', versions=2,
                 what=('The exact_partner arm is evaluated on the KREPP panel, leaf for leaf, so '
                       'the 0-3 references present only in the exact panel (median 1, max 3) '
                       'must never be observed in that arm.  The first version did not subtract '
                       'them and aborted on its own leaf-count assertion before writing any '
                       'record; the registered leaf sets are unchanged.'),
                 evidence_retained=['no partial stream was written'],
                 affects_measurements=False),
            dict(script='stage299_analyze.py', versions=3,
                 what=('Two robustness fixes in ANALYSIS code only: the travellers arm records '
                       'neither a removed-reference list (up to 7,931 ids per draw) nor '
                       'K_generating, so both accesses were made optional.  No measurement is '
                       'produced by this script; it reads the frozen streams.'),
                 evidence_retained=['stage299_analysis.json is the only output'],
                 affects_measurements=False)],
        derived_outputs=dict(
            stage299_regime_ratios_json=('structured/random ratios of the median fall in Theta_1, '
                                         'computed deterministically from the registered per-cell '
                                         'table in stage299_analysis.json; ratios are suppressed '
                                         'where the random median fall is zero'),
            stage299_real_incompleteness_md='the stage report',
            stage299_real_incompleteness_png='the stage figure'),
        protocol_code_hashes=protocol['code_hashes'],
        final_code_hashes=final,
        changed_since_registration=changed,
        output_hashes={name: common.sha(OUT / name) for name in sorted(
            p.name for p in OUT.iterdir() if p.is_file())})
    path.write_text(json.dumps(common.serial(out), indent=1))
    print('wrote', path)
    print('changed since registration:', changed)


if __name__ == '__main__':
    main()

