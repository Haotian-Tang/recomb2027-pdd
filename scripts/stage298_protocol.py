#!/usr/bin/env python3
"""Stage 298 protocol registration.  Run BEFORE any timed measurement."""
import json
import os
import sys
from datetime import datetime, timezone
from math import comb

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
import stage288_common as common            # noqa: E402
import stage298_common as s298              # noqa: E402

OUT = s298.OUT
PROTOCOL = OUT / 'stage298_protocol.json'

FROZEN = ['stage288_common.py', 'stage296_common.py', 'stage287_pdd_cases.py',
          'stage287_build_bank.py', 'stage287_verify_bank.py', 'stage286_newick.py',
          'stage42_benchmark.py', 'stage43_partial_fpt.py']
NEW = ['stage298_common.py', 'stage298_placement.py', 'stage298_sweep.py',
       'stage298_td.py', 'stage298_residuals.py', 'stage298_analyze.py']

WORKERS = 40
SUPPORT_CAP = {'1': 10 ** 9, '2': 10 ** 9, '3': 250000}
WALL = {'1': 1800.0, '2': 3600.0, '3': 3600.0, 'td': 2400.0}
ADDRESS_SPACE = 8 * 1024 ** 3
TD_CHUNK = 500


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    assert not PROTOCOL.exists(), 'Never overwrite frozen evidence'
    index = {e['name']: e for e in common.bank_index()}
    manifest = json.loads((OUT / 'stage298_placement_manifest.json').read_text())

    small = []
    for name, e in sorted(index.items()):
        if e['family'] != 'krepp_estimated_10k':
            continue
        case = json.loads(open(os.path.join(ROOT, e['path'])).read())
        small.append(dict(name=name, family=e['family'], path=e['path'],
                          sha256=e['sha256'],
                          K_generating=len(case['generating_sources'] or []),
                          reference_set_size=len(case['observed']),
                          n_edges=len(case['edges']),
                          instance=case['provenance']['instance'],
                          distance_source=case['provenance']['distance_source']))
    for c in manifest['cases']:
        small.append(dict(name=c['name'], family=c['family'], path=c['path'],
                          sha256=c['sha256'], K_generating=c['K_generating'],
                          reference_set_size=c['reference_set_size'],
                          n_edges=c['n_edges'], instance=c['instance'],
                          distance_source=c['provenance']['distance_source']))
    td = []
    for name, e in sorted(index.items()):
        if e['family'] != 'published_travelers_diarrhea':
            continue
        case = json.loads(open(os.path.join(ROOT, e['path'])).read())
        meta = case['provenance'].get('sample_metadata') or {}
        td.append(dict(name=name, family=e['family'], path=e['path'],
                       sha256=e['sha256'], sample=meta.get('sample'),
                       group=meta.get('group'),
                       reference_set_size=len(case['observed']),
                       n_edges=len(case['edges']),
                       distance_source=case['provenance']['distance_source']))

    protocol = dict(
        stage=298, track='real_data_tolerance_axis',
        registered_utc=datetime.now(timezone.utc).isoformat(),
        purpose=(
            'Every incompleteness measurement in this project so far runs on instances '
            'with synthetic planted sources.  Stage 296 built the instrument that makes a '
            'real-data statement possible without ground truth: Theta_1, the exact '
            'smallest tolerance at which a SINGLE source explains the observations, is a '
            'rational number computable from the data alone.  This stage computes it on '
            'three real collections and places them on one axis, in observation units and '
            'relative to each sample\'s median observed distance.'),
        claims_tested=[
            'Theta_1 is computable exactly on a real 16S backbone with |E| = 21,147 and '
            '|R| = 10,575 (29 travellers\' diarrhoea stool samples, no ground truth).',
            'Across all 296 real krepp ANI cases - not the 15-case subsample of stage 296 - '
            'the measured estimation residual against the exact partner is compared with '
            'Theta_1, Theta_2, Theta_3 to give the implied m(eps).',
            'A second real distance estimator (krepp placement distances) is placed on the '
            'same axis, or refused with a machine-readable reason.'],
        decision_rules=[
            'Theta_k is decided ONLY by frozen stage-296 code: theta_support, '
            'theta_min_over_supports, independent_deviation, witness_sources.',
            'The stage-298 k=1 closed form is accepted only where it agrees with the '
            'frozen LP by exact rational equality; every edge of every case is compared.',
            'A (case, k) cell is RUN when C(|E|, k) <= support_cap[k], CENSORED otherwise; '
            'censored cells are reported with their support count and the measured '
            'seconds per support, never estimated away.',
            'k = 2 and k = 3 are swept with cutoff = Theta_{k-1}; if no support is '
            'strictly below, Theta_k = Theta_{k-1} exactly and is recorded as such.',
            'm(eps) = min{k : Theta_k < eps} for eps > 0, stage 296\'s registered '
            'convention (strict inequality).',
            'No source count is claimed for any real sample: no ground truth exists for '
            'the travellers\' samples, and for the krepp cases the generating taxa are '
            'known but the estimated distances are not exactly explainable.',
            'Stage 291\'s inadmissible-floor construction is not used anywhere.'],
        workers=WORKERS, support_cap=SUPPORT_CAP, wall_seconds=WALL,
        address_space_bytes=ADDRESS_SPACE, td_chunk_edges=TD_CHUNK,
        small_cases=small, td_cases=td,
        counts=dict(krepp_ani=sum(1 for s in small if s['family'] == 'krepp_estimated_10k'),
                    krepp_placement=sum(1 for s in small if s['family'] == 'krepp_placement_10k'),
                    travelers=len(td),
                    cells_k1=sum(1 for s in small if comb(s['n_edges'], 1) <= SUPPORT_CAP['1']),
                    cells_k2=sum(1 for s in small if comb(s['n_edges'], 2) <= SUPPORT_CAP['2']),
                    cells_k3=sum(1 for s in small if comb(s['n_edges'], 3) <= SUPPORT_CAP['3']),
                    td_chunks=sum(-(-s['n_edges'] // TD_CHUNK) for s in td)),
        code_hashes={f: common.sha(os.path.join(ROOT, 'scripts', f))
                     for f in FROZEN + NEW
                     if os.path.exists(os.path.join(ROOT, 'scripts', f))},
        input_hashes=dict(
            placement_manifest=common.sha(OUT / 'stage298_placement_manifest.json'),
            stage296_calibration=common.sha(os.path.join(
                ROOT, 'results', 'stage296', 'stage296_calibration.json')),
            stage296_protocol=common.sha(os.path.join(
                ROOT, 'results', 'stage296', 'stage296_protocol.json')),
            bank_protocol=common.sha(os.path.join(
                ROOT, 'results', 'stage287_instance_bank_protocol.json'))),
        imported_without_recomputation=[
            'results/stage296/stage296_calibration.json: the 24 krepp ANI residual rows '
            'stage 296 measured (15 of which it placed against Theta_k).  Stage 298 '
            'extends the same formula to all 356 real cases and requires exact rational '
            'equality on those 24 rows as the validation of the extension.'],
    )
    PROTOCOL.write_text(json.dumps(common.serial(protocol), indent=1, sort_keys=True))
    print('registered', PROTOCOL)
    print(json.dumps(protocol['counts'], indent=1))
    miss = [f for f in FROZEN + NEW if f not in protocol['code_hashes']]
    print('missing code hashes:', miss)


if __name__ == '__main__':
    main()

