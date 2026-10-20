#!/usr/bin/env python3
"""Stage 299 step 0: register the protocol BEFORE any timed run.

Builds the complete draw plan for all four arms, verifies that its random
draws reproduce stage288_incompleteness.removal_draws exactly wherever the
panels coincide, hashes every input case file and every script in the import
closure, and writes results/stage299/stage299_protocol.json.
"""
import json
import sys
from datetime import datetime, timezone
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import stage288_common as common               # noqa: E402
import stage288_incompleteness as inc          # noqa: E402
import stage294_common as s294                 # noqa: E402
import stage299_common as s299                 # noqa: E402

OUT = s299.OUT
PROTOCOL = OUT / 'stage299_protocol.json'
WORKERS = 40
SMALL_WALL, TD_WALL = 600.0, 1800.0
AS_BYTES = 8 * 1024 ** 3


def build_plan():
    index = s299.bank_by_name()
    manifest = json.loads((s299.S298 / 'stage298_placement_manifest.json').read_text())
    td_summary = json.loads((s299.S298 / 'stage298_td_summary.json').read_text())
    p298 = json.loads((s299.S298 / 'stage298_protocol.json').read_text())
    ranks = s294.load_ranks()
    rows, cases, checks = [], [], dict(random_import_checked=0, random_import_agree=0,
                                       sources_identical=0, sources_compared=0)

    def add(arm, name, path, sha256, meta, draws):
        for d in draws:
            rows.append(dict(arm=arm, name=name, level=d['level'], regime=d['regime'],
                             replicate=d['replicate'], removed=d['removed'],
                             raw_removed=d['raw_removed'], n_removed=d['n_removed'],
                             n_retained=d['n_retained'], detail=d['detail'], **meta))
        cases.append(dict(arm=arm, name=name, path=path, sha256=sha256,
                          n_leafsets=len(draws), **meta))

    # ---- krepp ANI + its paired exact partner, identical leaf sets
    for name, entry in sorted(index.items()):
        if entry['family'] != 'krepp_estimated_10k':
            continue
        stem = name[:-len('_krepp_estimated')]
        ex = index[stem + '_exact']
        noisy_raw = s299.load_json(entry['path'], entry['sha256'])
        exact_raw = s299.load_json(ex['path'], ex['sha256'])
        noisy = s299.prepare_case(noisy_raw, name)
        exact = s299.prepare_case(exact_raw, ex['name'])
        checks['sources_compared'] += 1
        checks['sources_identical'] += int(noisy['generating_sources'] == exact['generating_sources'])
        draws = s299.small_draws(noisy, ex['name'], ranks)
        # import check: identical panels must reproduce the frozen random draws
        if sorted(noisy['observed']) == sorted(exact['observed']):
            frozen_case = common.load_case(ex)
            for level in s299.LEVELS:
                frozen = [d for d in inc.removal_draws(frozen_case, level)
                          if d['regime'] == 'random']
                mine = [d for d in draws if d['level'] == level and d['regime'] == 'random']
                for a, b in zip(frozen, mine):
                    checks['random_import_checked'] += 1
                    checks['random_import_agree'] += int(
                        [int(r) for r in a['removed']] == [int(r) for r in b['removed']])
        meta = dict(stem=stem, family=entry['family'], K_generating=entry['K_generating'],
                    reference_set_size=len(noisy['observed']), n_edges=len(noisy['edges']),
                    median_observed=str(s299.median_observed(list(noisy['observed'].values()))),
                    exact_partner=ex['name'],
                    exact_only_references=sorted(set(exact['observed']) - set(noisy['observed'])))
        add('krepp_ani', name, entry['path'], entry['sha256'], meta, draws)
        add('exact_partner', ex['name'], ex['path'], ex['sha256'],
            dict(meta, paired_with=name), draws)

    # ---- krepp placement cases (stage 298 step 1 conversions)
    for row in sorted(manifest['cases'], key=lambda r: r['instance']):
        case = s299.prepare_case(s299.load_json(row['path'], row['sha256']), row['name'])
        ex = index['k10k_%s_exact' % row['instance']]
        draws = s299.small_draws(case, ex['name'], ranks)
        meta = dict(stem='k10k_%s' % row['instance'], family='krepp_placement_10k',
                    K_generating=row['K_generating'],
                    reference_set_size=len(case['observed']), n_edges=len(case['edges']),
                    median_observed=str(s299.median_observed(list(case['observed'].values()))),
                    exact_partner=ex['name'], exact_only_references=[])
        add('krepp_placement', row['name'], row['path'], row['sha256'], meta, draws)

    # ---- travellers' diarrhoea: nearest regime anchored on stage 298's witness
    td_rows = {r['name']: r for r in td_summary['rows']}
    for row in sorted(p298['td_cases'], key=lambda r: r['name']):
        case = s299.prepare_case(s299.load_json(row['path'], row['sha256']), row['name'])
        s = td_rows[row['name']]
        full_source = dict(edge=int(s['argmin_edge']), p=F(1), x=F(str(s['argmin_x'])))
        draws = s299.td_draws(case, full_source)
        meta = dict(stem=row['sample'], family=row['family'], K_generating=None,
                    reference_set_size=len(case['observed']), n_edges=len(case['edges']),
                    median_observed=str(s299.median_observed(list(case['observed'].values()))),
                    exact_partner=None, exact_only_references=[],
                    stage298_theta1=str(s['theta1']), stage298_argmin=int(s['argmin_edge']),
                    stage298_argmin_x=str(s['argmin_x']))
        add('travelers', row['name'], row['path'], row['sha256'], meta, draws)
    return rows, cases, checks


ONE_STR = '1'


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    assert not PROTOCOL.exists(), 'Never overwrite frozen evidence'
    rows, cases, checks = build_plan()
    assert checks['random_import_agree'] == checks['random_import_checked'], \
        'random draw import does not reproduce stage 288: %s' % checks
    assert checks['sources_identical'] == checks['sources_compared'], checks
    by_arm = {}
    for r in rows:
        by_arm.setdefault(r['arm'], []).append(r)
    code = {}
    for entry in ('stage299_common.py', 'stage299_protocol.py', 'stage299_small.py',
                  'stage299_td.py', 'stage299_validate.py', 'stage299_analyze.py',
                  'stage298_common.py', 'stage296_common.py', 'stage288_common.py',
                  'stage288_incompleteness.py', 'stage294_common.py'):
        path = ROOT / 'scripts' / entry
        code[entry] = common.sha(path) if path.exists() else 'NOT_YET_WRITTEN'
    protocol = dict(
        stage=299, track='real_data_reference_incompleteness',
        registered_utc=datetime.now(timezone.utc).isoformat(),
        purpose=('Close the loop the paper needs: run the stage-288 reference-removal '
                 'regimes through stage 298\'s ground-truth-free instrument (Theta_1) on '
                 'real estimated distances, and ask whether the regime ordering that holds '
                 'in the exact model survives on real data.'),
        quantity=('Theta_1 = min over every original edge e, and over the CLOSED mass set '
                  '{alpha, beta, y >= 0, alpha + beta = 1}, of max_r |pred_r - d_r| with '
                  'pred_r = y + (1-x) d(r,u) + x d(r,v) on e.  It is the exact smallest '
                  'tolerance at which a SINGLE source explains the retained observations, '
                  'needs no error model, no ground truth and no admissible pendant floor, '
                  'and is a rational number.'),
        arms=dict(
            exact_partner=('the 296 exact_10k partners of the krepp cases, evaluated on the '
                           'krepp panel and its identical leaf sets - the paired exact-model '
                           'baseline'),
            krepp_ani='the 296 krepp ANI cases (real k-mer distance estimates)',
            krepp_placement='the 60 krepp placement cases converted in stage 298 step 1',
            travelers=('the 29 travellers\'-diarrhoea stool samples, DEPP placement distances '
                       'on a 21,147-edge / 10,575-leaf 16S backbone, no ground truth')),
        draw_plan=dict(
            imported_from='stage288_incompleteness.removal_draws',
            seed=inc.SEED, levels=list(s299.LEVELS), regimes=list(inc.REGIMES) + ['taxon'],
            replicates=dict(inc.REPLICATES), taxon_orders=list(s294.TAXON_ORDERS),
            taxon_ranks=list(s294.TAXON_RANKS),
            seed_key_convention=('stage 294\'s: the random seed key uses the PAIRED EXACT case '
                                 'name, so the exact, ANI and placement arms share leaf sets '
                                 'leaf for leaf'),
            clade_bookkeeping=('clade_removals returns whole-clade TREE leaves; the effective '
                               'removal is that set intersected with the arm\'s observed panel, '
                               'and every count and level uses the effective one.  Both sets '
                               'are recorded per draw.'),
            travelers_nearest=('no ground truth exists for a stool sample, so the nearest-to-'
                               'source regime is anchored on the stage-298 full-panel Theta_1 '
                               'attaining point through the same frozen nearest_source_removal '
                               'routine, and is reported as `nearest_theta1`, never pooled with '
                               'the truth-based `nearest` of the other arms.  No taxonomic '
                               'regime: the 16S backbone carries no accession labels.'),
            import_verification=checks,
            digest=s299.draw_digest(rows),
            n_leafsets=len(rows),
            n_leafsets_by_arm={k: len(v) for k, v in sorted(by_arm.items())}),
        decision=dict(
            small_arms=('Theta_1 decided by the FROZEN per-support LP '
                        '(stage296_common.theta_support) on every one of the |E| edges with no '
                        'cutoff, and independently by the stage-298 k=1 closed form; exact '
                        'rational equality required edge by edge; the attaining witness '
                        're-verified with the frozen min(d(r,u)+xL, d(r,v)+(1-x)L) form.'),
            travelers_arm=('Theta_1 decided by the stage-298 closed form over all 21,147 edges '
                           '- licensed by stage 298, which decided 653,807 edges of these same '
                           'trees twice with 0 disagreements - with the attaining edge plus a '
                           'seeded random sample of %d edges per leaf set re-decided by the '
                           'frozen per-support LP, and the witness re-verified with the frozen '
                           'min(...) form.  The frozen LP alone would cost 0.512 core-seconds '
                           'per edge, i.e. 3.0 core-hours per leaf set and 2,270 core-hours '
                           'for the 754 leaf sets.' % s299.TD_SAMPLE_EDGES),
            ties=('the attaining set (every edge whose theta equals Theta_1) is recorded; the '
                  'canonical attainer is the smallest such edge index, stage 298\'s '
                  'convention; migration is reported canonically AND as the minimum distance '
                  'between the two attaining sets.')),
        claims_tested=[
            'C1 MONOTONICITY CHECK, not a finding: removing references deletes LP rows, so '
            'Theta_1 can only fall or stay.  Any strict rise is a stop-and-report bug.',
            'C2 CROSS-STAGE: the full-panel Theta_1 of every case must equal stage 298\'s '
            'recorded value exactly (356 small cases, 29 travellers\' samples).',
            'C3 REGIME ORDERING (the measurement): per (arm, level, regime), the median fall '
            'ratio Theta_1(draw)/Theta_1(full) and the fraction of draws that fall strictly.  '
            'The exact-model comparator is the exact_partner arm on identical leaf sets, plus '
            'stage 296\'s 1,075 frozen exact draws re-analysed without new runs.  Reported '
            'with equal prominence whether or not the ordering survives.',
            'C4 DOSE-RESPONSE: the smallest removal level at which Theta_1 drops below the '
            'measured estimator residual.  Two lines - worst-leaf and median-leaf - each in '
            'two variants: recomputed on the RETAINED leaves (primary; it is what a '
            'practitioner at the reduced panel actually faces) and fixed at the full panel '
            '(secondary; the literal reading of a horizontal line).  Cases already below the '
            'line at the full panel are excluded and counted, never recorded as a threshold '
            'of zero.',
            'C5 ATTAINING-EDGE MIGRATION: exact tree distance between the full-panel attaining '
            'point and the post-removal attaining point, normalised by the case\'s own median '
            'observed distance, plus topological hop count and the terminal-edge fraction.'],
        not_established=[
            'No source count is claimed for any real sample.',
            'Theta_2 and above are not computed here; stage 298 measured Theta_2 on the '
            'travellers\' backbone at 39,594 core-hours per sample.',
            'The travellers\' arm has no estimator residual, so it carries no dose-response.',
            'The `nearest_theta1` regime is data-derived, not truth-derived, and is not '
            'comparable draw-by-draw with the truth-based `nearest`.'],
        budgets=dict(workers=WORKERS, small_wall_seconds=SMALL_WALL, td_wall_seconds=TD_WALL,
                     address_space_bytes=AS_BYTES,
                     calibration=('sized from stage 298\'s measured costs, not guessed: '
                                  '0.0051 s median per support on the small cases (356 full '
                                  'panels = 216 core-seconds) and 3.7 s per full closed-form '
                                  'sweep / 0.512 core-seconds per frozen LP edge on the '
                                  'travellers\' backbone.')),
        analysis=('Retain every run including timeouts and errors; never impute.  Every '
                  'comparison is made against the case\'s OWN full-panel Theta_1 in the SAME '
                  'arm, so noise-induced and incompleteness-induced collapse are not '
                  'confounded.  Normalisation is the case\'s own median observed distance '
                  '(stage 298\'s convention), the only denominator the travellers\' samples '
                  'admit.'),
        cases=cases, leafsets=rows, code_hashes=code)
    PROTOCOL.write_text(json.dumps(common.serial(protocol), indent=1))
    print('registered %s' % PROTOCOL)
    print('leaf sets %d: %s' % (len(rows), {k: len(v) for k, v in sorted(by_arm.items())}))
    print('cases %d; import checks %s' % (len(cases), checks))
    print('digest %s' % protocol['draw_plan']['digest'])


if __name__ == '__main__':
    main()

