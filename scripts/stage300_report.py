#!/usr/bin/env python3
"""Stage 300 report: assembles stage300_own_fit_comparison.md from the records.

Every number in the document is read from the stage-300 JSON/JSONL records at
run time; none is transcribed by hand.
"""
import json
import sys
from collections import Counter
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import stage300_common as s300                   # noqa: E402

OUT = s300.OUT


def read_json(name, default=None):
    path = OUT / name
    return json.loads(path.read_text()) if path.exists() else default


def read_rows(name):
    path = OUT / name
    if not path.exists():
        return []
    return [json.loads(l) for l in open(path)]


def fmt(value, digits=3):
    if value is None:
        return 'n/a'
    if isinstance(value, float):
        if value != 0 and (abs(value) < 1e-3 or abs(value) >= 1e5):
            return '%.*e' % (digits - 1, value)
        return ('%.*f' % (digits, value)).rstrip('0').rstrip('.')
    return str(value)


def q(cell, key='median', digits=3):
    if not cell:
        return 'n/a'
    return fmt(cell.get(key), digits)


def provenance():
    lines = ['| collection | cases | DecoDiPhy runs harvested | from | runs newly executed | runs not usable |',
             '|---|---|---|---|---|---|']
    harvest = read_rows('stage300_harvest.jsonl')
    new = []
    for name in ('stage300_runs_placement.jsonl', 'stage300_runs_travelers.jsonl',
                 'stage300_pilot_travelers.jsonl',
                 'stage300_travelers_censored.jsonl'):
        new.extend(read_rows(name))
    for collection in s300.COLLECTIONS:
        cases = len(s300.case_entries(collection))
        h = [r for r in harvest if r['collection'] == collection]
        n = [r for r in new if r.get('collection') == collection]
        stages = sorted({r['source_stage'] for r in h})
        bad = sum(1 for r in h if r['status'] != 'complete' or not r['sources'])
        bad += sum(1 for r in n if (r.get('result') or {}).get('status') != 'complete')
        lines.append('| %s | %d | %d | %s | %d | %d |'
                     % (collection, cases, len(h),
                        ', '.join('stage %d' % s for s in stages) or '-', len(n), bad))
    return '\n'.join(lines)


def three_way(analysis):
    lines = ['| collection | arm | cases with an answer | k returned (median) | '
             'R_dd/Theta_1 median | R_dd >= Theta_1 | truth/Theta_1 median | '
             'truth >= Theta_1 | DecoDiPhy beats truth (L-inf) |', '|---|' + '---|' * 8]
    for collection in ('krepp_ani', 'krepp_placement'):
        cell = analysis['by_collection'].get(collection)
        if not cell:
            continue
        for arm in ('dd_default', 'dd_true_k', 'dd_resolution'):
            a = cell['arms'].get(arm)
            if not a:
                continue
            lines.append('| %s | %s | %d/%d | %s | %s | %d/%d | %s | %d/%d | %d/%d |'
                         % (collection, arm, a['n_cases_with_answer'], cell['n_cases'],
                            q(a['k_returned']), q(a['r_dd_over_theta1']),
                            a['n_r_dd_ge_theta1'], a['n_cases_with_answer'],
                            q(a.get('truth_over_theta1')),
                            a.get('n_truth_ge_theta1', 0), a.get('n_with_truth', 0),
                            a.get('n_dd_beats_truth_in_linf', 0), a.get('n_with_truth', 0)))
    return '\n'.join(lines)


def own_norm(analysis):
    lines = ['| collection | arm | SSR(DecoDiPhy)/Lambda_1 median | beats one source in SSR | '
             'beats truth in SSR | truth beats one source in SSR | '
             'recomputed SSR vs its reported loss (median rel.) |', '|---|' + '---|' * 6]
    for collection in ('krepp_ani', 'krepp_placement', 'travelers'):
        cell = analysis['by_collection'].get(collection)
        if not cell:
            continue
        for arm in ('dd_default', 'dd_true_k'):
            a = cell['arms'].get(arm)
            if not a:
                continue
            lines.append('| %s | %s | %s | %d/%d | %s | %s | %s |'
                         % (collection, arm, q(a['ssr_over_lambda1']),
                            a['n_dd_beats_one_source_in_ssr'], a['n_cases_with_answer'],
                            ('%d/%d' % (a['n_dd_beats_truth_in_ssr'], a['n_with_truth'])
                             if 'n_with_truth' in a else 'no truth'),
                            ('%d/%d' % (a['n_truth_beats_one_source_in_ssr'], a['n_with_truth'])
                             if 'n_with_truth' in a else 'no truth'),
                            q(a['loss_rel_diff'], digits=2)))
    return '\n'.join(lines)


def resolved_table(analysis):
    lines = ['| collection | arm | cases re-solved | R(re-solved) >= Theta_1 | '
             'ratio median | returned-configuration verdict |', '|---|' + '---|' * 5]
    for collection in ('krepp_ani', 'krepp_placement'):
        cell = analysis['by_collection'].get(collection, {})
        res = (cell.get('resolved') or {})
        for arm in ('dd_default', 'dd_true_k'):
            a = res.get(arm)
            if not a:
                continue
            ret = cell['arms'].get(arm, {})
            lines.append('| %s | %s | %d/%d | %d/%d | %s | %d/%d |'
                         % (collection, arm, a['n_ok'], a['n'],
                            a['n_linf_ge_theta1'], a['n_ok'],
                            q(a['linf_over_theta1']),
                            ret.get('n_r_dd_ge_theta1', 0),
                            ret.get('n_cases_with_answer', 0)))
    return '\n'.join(lines)


def breakdown():
    """R_dd >= Theta_1 broken down by true K and by reference-panel size."""
    import csv as _csv
    path = OUT / 'stage300_three_way.csv'
    if not path.exists():
        return ''
    rows = [r for r in _csv.DictReader(open(path)) if r['arm'] == 'dd_default']
    lines = ['| collection | stratum | cases | R_dd >= Theta_1 | R_dd/Theta_1 median | '
             'k returned (median) |', '|---|' + '---|' * 5]
    for collection in ('krepp_ani', 'krepp_placement'):
        sel = [r for r in rows if r['collection'] == collection and r['K_true']]
        if not sel:
            continue
        strata = []
        for K in sorted({int(r['K_true']) for r in sel}):
            strata.append(('true K = %d' % K, [r for r in sel if int(r['K_true']) == K]))
        ordered = sorted(sel, key=lambda r: int(r['n_leaves']))
        third = len(ordered) // 3
        for label, group in (('panel, lower third', ordered[:third]),
                             ('panel, middle third', ordered[third:2 * third]),
                             ('panel, upper third', ordered[2 * third:])):
            sizes = [int(r['n_leaves']) for r in group]
            strata.append(('%s (%d-%d references)' % (label, min(sizes), max(sizes)), group))
        for label, group in strata:
            ratios = sorted(float(r['r_dd_over_theta1']) for r in group)
            ks = sorted(int(r['k_returned']) for r in group)
            ge = sum(1 for r in group if r['r_dd_ge_theta1'] == 'True')
            lines.append('| %s | %s | %d | %d (%.0f%%) | %s | %s |'
                         % (collection, label, len(group), ge, 100.0 * ge / len(group),
                            fmt(ratios[len(ratios) // 2]), ks[len(ks) // 2]))
    return '\n'.join(lines)

def stool_table():
    fits = read_rows('stage300_fits_travelers.jsonl')
    bases = {b['name']: b for b in read_rows('stage300_baselines_travelers.jsonl')}
    td298 = json.loads((ROOT / 'results' / 'stage298' /
                        'stage298_td_summary.json').read_text())
    groups = {r['name']: (r.get('sample'), r.get('group')) for r in td298['rows']}
    censored = {r['name']: r for r in read_rows('stage300_travelers_censored.jsonl')}
    lines = ['| sample | group | references | Theta_1 | Lambda_1 (SSR) | '
             'L-inf residual of the Lambda_1 configuration | DecoDiPhy answer |',
             '|---|' + '---|' * 6]
    done = 0
    for name in sorted(bases, key=lambda n: groups.get(n, (n, ''))[0] or n):
        base = bases[name]
        rows = [r for r in fits if r['name'] == name and r.get('evaluated')]
        sample, group = groups.get(name, (name, ''))
        theta1 = float(F(base['theta1']))
        lam = (base.get('lambda1') or {})
        lam_linf = (float(F(lam['linf_at_l2_optimum']))
                    if lam.get('linf_at_l2_optimum') else None)
        if rows:
            done += 1
            row = rows[0]
            answer = ('k = %d, R_dd = %s (%s x Theta_1), SSR = %s'
                      % (row['k'], fmt(row['r_dd_float']),
                         fmt(row['r_dd_float'] / theta1), fmt(row['ssr_dd'])))
        elif name in censored:
            cell = censored[name]['result']
            answer = ('censored: stopped at %.0f s (%.1f h) with %.1f GiB resident, '
                      'no answer returned'
                      % (cell['elapsed_seconds'], cell['elapsed_seconds'] / 3600.0,
                         cell['rss_kib_at_stop'] / 1048576.0))
        else:
            answer = 'not attempted (see section 5)'
        lines.append('| %s | %s | %d | %s | %s | %s | %s |'
                     % (sample, group, base['n_leaves'], fmt(theta1),
                        fmt(lam.get('ssr')), fmt(lam_linf), answer))
    return '\n'.join(lines), done, len(bases)


def main():
    analysis = read_json('stage300_analysis.json')
    lossgap = read_json('stage300_lossgap.json', {})
    route_e = read_json('stage300_route_e.json', {})
    verify = read_json('stage300_verification_small.json', {})
    verify_td = read_json('stage300_verification_travelers.json', {})
    probe = read_json('stage300_tdprobe_td_sample_101_published.json', {})
    losscase = read_json('stage300_losscase_k10k_tree_15_k3_s1_exp_krepp_estimated.json', {})
    protocol = read_json('stage300_protocol.json', {})
    stool, stool_done, stool_total = stool_table()
    table = dict(provenance=provenance(), three_way=three_way(analysis),
                 own_norm=own_norm(analysis), resolved=resolved_table(analysis),
                 stool=stool)
    payload = dict(analysis=analysis, lossgap=lossgap, route_e=route_e, verify=verify,
                   verify_td=verify_td, probe=probe, losscase=losscase,
                   protocol=protocol, tables=table,
                   stool_done=stool_done, stool_total=stool_total)
    (OUT / 'stage300_report_payload.json').write_text(json.dumps(s300.serial(payload), indent=1))
    print(json.dumps(dict(stool_done=stool_done, stool_total=stool_total), indent=1))
    for key, value in table.items():
        print('\n### ' + key + '\n' + value)


if __name__ == '__main__':
    main()
