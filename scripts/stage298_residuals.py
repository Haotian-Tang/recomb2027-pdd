#!/usr/bin/env python3
"""Stage 298 step 4: the real estimation error, in the same units as eps.

Stage 296 defined it and computed it for 24 of the 296 krepp ANI cases:
the L-infinity residual max_r |d_estimated(r) - d_exact(r)| on the shared leaf
set, with the median residual and both relative to the exact partner's median
observed distance.  This script applies that definition - the same formula,
the same shared-leaf rule, the same sorted(...)[len // 2] median convention -
to all 296 krepp cases and all 60 placement cases, and VALIDATES the extension
by requiring exact rational equality with every one of stage 296's 24 rows.

It also checks that each placement case converted in step 1 carries the same
reference tree and leaf labelling as the bank's exact case for the same
instance, so the two observation vectors are comparable leaf by leaf.
"""
import json
import os
import sys
from fractions import Fraction as F

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
import stage288_common as common            # noqa: E402
import stage298_common as s298              # noqa: E402

OUT = s298.OUT
S296_CALIB = ROOT + '/results/stage296/stage296_calibration.json'


def residual_row(exact_case, noisy_case):
    """Stage 296's calibration formula, unchanged."""
    exact = {int(r): F(d) for r, d in exact_case['observed'].items()}
    median_observed = sorted(exact.values())[len(exact) // 2]
    noisy = {int(r): F(d) for r, d in noisy_case['observed'].items()}
    shared = sorted(set(noisy) & set(exact))
    if not shared:
        return None
    resid = sorted(abs(noisy[r] - exact[r]) for r in shared)
    own = sorted(noisy.values())
    own_median = own[len(own) // 2]
    return dict(
        n_shared=len(shared), n_exact_leaves=len(exact),
        median_observed=median_observed,
        own_median_observed=own_median,
        residual_median=resid[len(resid) // 2],
        residual_max=resid[-1],
        residual_median_rel=float(resid[len(resid) // 2] / median_observed),
        residual_max_rel=float(resid[-1] / median_observed),
        residual_median_rel_own=float(resid[len(resid) // 2] / own_median),
        residual_max_rel_own=float(resid[-1] / own_median))


def main():
    out_path = OUT / 'stage298_residuals.json'
    assert not out_path.exists(), 'Never overwrite frozen evidence'
    index = {e['name']: e for e in common.bank_index()}
    manifest = json.loads((OUT / 'stage298_placement_manifest.json').read_text())
    placement = {c['instance']: c for c in manifest['cases']}

    def load(entry):
        path = os.path.join(ROOT, entry['path'])
        assert common.sha(path) == entry['sha256'], 'hash mismatch %s' % entry['path']
        return json.loads(open(path).read())

    rows, tree_checks = [], []
    for name, entry in sorted(index.items()):
        if entry['family'] != 'krepp_estimated_10k':
            continue
        inst = json.loads(open(os.path.join(ROOT, entry['path'])).read())['provenance']['instance'] \
            if 'provenance' in json.loads(open(os.path.join(ROOT, entry['path'])).read()) else None
        stem = name[:-len('_krepp_estimated')]
        ex = index.get(stem + '_exact')
        if ex is None:
            rows.append(dict(case=name, family=entry['family'], label='krepp_ani',
                             status='no_exact_partner'))
            continue
        r = residual_row(load(ex), load(entry))
        rows.append(dict(case=name, exact_partner=stem + '_exact',
                         family='krepp_estimated_10k', label='krepp_ani',
                         instance=inst, status='ok', **r))
    for inst, c in sorted(placement.items()):
        ex = index.get('k10k_%s_exact' % inst)
        if ex is None:
            rows.append(dict(case=c['name'], family='krepp_placement_10k',
                             label='krepp_placement', status='no_exact_partner'))
            continue
        exact_case = load(ex)
        path = os.path.join(ROOT, c['path'])
        assert common.sha(path) == c['sha256'], 'hash mismatch %s' % c['path']
        pcase = json.loads(open(path).read())
        same_edges = exact_case['edges'] == pcase['edges']
        same_labels = exact_case['leaf_labels'] == pcase['leaf_labels']
        same_sources = exact_case['generating_sources'] == pcase['generating_sources']
        tree_checks.append(dict(instance=inst, same_edges=same_edges,
                                same_leaf_labels=same_labels,
                                same_generating_sources=same_sources))
        r = residual_row(exact_case, pcase)
        rows.append(dict(case=c['name'], exact_partner='k10k_%s_exact' % inst,
                         family='krepp_placement_10k', label='krepp_placement',
                         instance=inst, status='ok', **r))

    # --- validation against stage 296's 24 rows, exact equality
    calib = json.loads(open(S296_CALIB).read())
    mine = {r['case']: r for r in rows if r.get('status') == 'ok'}
    agree, disagree = 0, []
    for old in calib['rows']:
        if old['label'] != 'krepp_ani':
            continue
        new = mine.get(old['partner'])
        if new is None:
            disagree.append(dict(case=old['partner'], reason='absent_from_extension'))
            continue
        same = (F(old['residual_median']) == F(new['residual_median'])
                and F(old['residual_max']) == F(new['residual_max'])
                and F(old['median_observed']) == F(new['median_observed'])
                and old['n_shared'] == new['n_shared'])
        if same:
            agree += 1
        else:
            disagree.append(dict(case=old['partner'],
                                 old=dict(median=old['residual_median'],
                                          max=old['residual_max'],
                                          n_shared=old['n_shared']),
                                 new=dict(median=str(new['residual_median']),
                                          max=str(new['residual_max']),
                                          n_shared=new['n_shared'])))
    out = dict(
        stage=298, step='residuals',
        definition=calib['purpose'],
        note=('residual_*_rel keeps stage 296\'s denominator (the exact partner\'s '
              'median observed distance); residual_*_rel_own uses the case\'s own '
              'median observed distance, which is the only denominator available '
              'for the travellers\' samples, so the two collections can share an axis.'),
        n_rows=len(rows),
        validation_against_stage296=dict(n_compared=agree + len(disagree),
                                        agree=agree, disagree=disagree),
        tree_identity_checks=tree_checks,
        rows=rows)
    out_path.write_text(json.dumps(common.serial(out), indent=1))
    print('rows %d; stage-296 agreement %d/%d; tree identity: edges %d/%d, labels %d/%d'
          % (len(rows), agree, agree + len(disagree),
             sum(t['same_edges'] for t in tree_checks), len(tree_checks),
             sum(t['same_leaf_labels'] for t in tree_checks), len(tree_checks)))
    if disagree:
        print('DISAGREEMENTS', json.dumps(disagree[:3], indent=1))


if __name__ == '__main__':
    main()

