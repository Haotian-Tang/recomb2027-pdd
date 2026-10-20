#!/usr/bin/env python3
"""Stage295 report: assemble the ceiling diagnosis from the raw Stage295 records.

Writes results/stage295_ceiling_report.json and stage295_ceiling_diagnosis.md.
Every number in the document is computed here from the retained run records;
nothing is transcribed by hand.
"""
import json, os, glob, sys
from collections import defaultdict
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RES = os.path.join(ROOT, 'results')


def jload(p, default=None):
    try:
        return json.load(open(p))
    except Exception:
        return default


def runs_in(subdir, prefix):
    out = []
    for p in sorted(glob.glob(os.path.join(RES, subdir, prefix + '*.json'))):
        if os.path.basename(p).startswith('all_'):
            continue
        d = jload(p)
        if d:
            out.append(d)
    return out


def q(xs, p):
    return float(np.percentile(xs, p)) if xs else None


def main():
    census = jload(os.path.join(RES, 'stage295_static_census.json'), [])
    ana = jload(os.path.join(RES, 'stage295_ceiling_analysis.json'), {})
    val = runs_in('stage295_validation', 'val_') + runs_in('stage295_revalidation', 'val_')
    push_v1 = runs_in('stage295_push', 'push_')
    push_v2 = runs_in('stage295_push2', 'push_')
    for r in push_v1:
        r['batch'] = 'push_v1_unbounded_cache'
    for r in push_v2:
        r['batch'] = 'push_v2_bounded_cache'
    push = push_v1 + push_v2
    rep = {'inputs': {'census_rows': len(census), 'validation_rows': len(val),
                      'push_rows': len(push),
                      'profile_runs': ana.get('n_profile_runs')}}

    # ---- validation -----------------------------------------------------
    agree, disagree, speed = 0, [], defaultdict(list)
    identity_ok, identity_rows = 0, 0
    for r in val:
        if not r.get('ok'):
            continue
        res = r.get('results') or {}
        mins = r.get('minima') or {}
        if len(mins) > 1:
            if r.get('agree'):
                agree += 1
            else:
                disagree.append({'name': r['row']['name'], 'cap': r['row']['cap'],
                                 'minima': mins})
        tag = r['row'].get('tag', '')
        f = res.get('frozen') or {}
        if tag.startswith('identity_') and f.get('ok'):
            fm = res.get('full_mask_memo') or {}
            if fm.get('ok'):
                identity_rows += 1
                if (f.get('statistics', {}).get('mass_faces')
                        == fm.get('statistics', {}).get('mass_faces')):
                    identity_ok += 1
        p = res.get('reduced_memo') or {}
        if f.get('ok') and p.get('ok') and p.get('elapsed', 0) > 0 \
                and f.get('status') == 'complete':
            speed[r['row']['cap']].append(f['elapsed'] / p['elapsed'])
    rep['validation'] = {
        'comparable_cases': agree + len(disagree),
        'agreeing': agree, 'disagreements': disagree,
        'full_mask_identity_rows': identity_rows,
        'full_mask_identical_face_counts': identity_ok,
        'speedup_by_cap': {str(k): {'n': len(v), 'median': q(v, 50),
                                    'min': min(v), 'max': max(v)}
                           for k, v in sorted(speed.items())}}

    # ---- ceiling push ---------------------------------------------------
    rows = []
    for r in push:
        pr = r.get('primary') or {}
        mc = r.get('minimality_check') or {}
        rows.append({'batch': r.get('batch'), 'name': r['row']['name'], 'cap': r['row']['cap'],
                     'reference_set_size': r.get('reference_set_size'),
                     'popcount_order': r['row'].get('popcount_order'),
                     'status': pr.get('status'), 'minimum': pr.get('minimum'),
                     'seconds': pr.get('elapsed'),
                     'mass_faces': (pr.get('statistics') or {}).get('mass_faces'),
                     'vertex_bases': (pr.get('statistics') or {}).get('vertex_bases'),
                     'skipped_forced_faces': (pr.get('statistics') or {}).get('skipped_forced_faces'),
                     'witness_exactly_verified': pr.get('witness_exactly_verified'),
                     'minimality_certified': mc.get('certifies_minimality'),
                     'error': r.get('error')})
    done = [x for x in rows if x['status'] == 'complete']
    rep['push'] = {
        'rows': rows, 'n_rows': len(rows), 'n_complete': len(done),
        'n_minimality_certified': sum(1 for x in done if x['minimality_certified']),
        'n_witness_verified': sum(1 for x in done if x['witness_exactly_verified']),
        'completed_reference_sizes': sorted({x['reference_set_size'] for x in done}),
        'seconds_quantiles': {'min': min((x['seconds'] for x in done), default=None),
                              'median': q([x['seconds'] for x in done], 50),
                              'max': max((x['seconds'] for x in done), default=None)},
        'minima_observed': sorted({x['minimum'] for x in done if x['minimum']}),
        'failures': [x for x in rows if x['status'] != 'complete'],
        'failures_current_batch': [x for x in rows if x['status'] != 'complete'
                                   and x['batch'] == 'push_v2_bounded_cache'],
        'note': ('push_v1 ran with the unbounded cover cache that Amendment 1 '
                 'fixed; its MemoryError rows are superseded by push_v2 and are '
                 'retained only for the record.')}

    # ---- headline case: frozen projection vs prototype measurement -------
    HEAD = 'k10k_tree_16_k10_s5_exp_exact'
    cen = next((c for c in census if HEAD in c.get('case', '') and c.get('cap') == 10), None)
    proj = ana.get('projection_K10') or {}
    head_proto = [x for x in done if x['name'] == HEAD and not x['popcount_order']
                  and x['batch'] == 'push_v2_bounded_cache'] or \
                 [x for x in done if x['name'] == HEAD and not x['popcount_order']]
    rep['headline_case'] = {
        'name': HEAD,
        'reference_set_size': cen.get('reference_set_size') if cen else None,
        'backbone_edges': cen.get('backbone_edges') if cen else None,
        'exact_zero_faces_frozen': cen.get('mass_faces_total') if cen else None,
        'exact_zero_faces_reduced': cen.get('mass_faces_reduced') if cen else None,
        'face_reduction_factor': cen.get('reduction_factor') if cen else None,
        'frozen_projected_seconds': proj.get('projected_frozen_seconds'),
        'frozen_projected_days': proj.get('projected_frozen_days'),
        'prototype_seconds': head_proto[0]['seconds'] if head_proto else None,
        'prototype_minimum': head_proto[0]['minimum'] if head_proto else None,
        'prototype_minimality_certified': head_proto[0]['minimality_certified'] if head_proto else None,
        'end_to_end_speedup': (proj.get('projected_frozen_seconds') / head_proto[0]['seconds']
                               if head_proto and head_proto[0]['seconds'] and
                               proj.get('projected_frozen_seconds') else None)}
    rep['phase_profile'] = {'prefix_regression': ana.get('phase_cost_K10'),
                            'counters_at_stop': ana.get('lp_breakdown'),
                            'cprofile': {k: {'wall_s': v['wall_s'],
                                             'mass_faces': v['mass_faces'],
                                             'vertex_bases': v['vertex_bases'],
                                             'cumtime_key_frames': v['cumtime_key_frames'],
                                             'tottime_shares': v['shares']}
                                         for k, v in (ana.get('cprofile') or {}).items() if v}}
    rep['growth'] = ana.get('growth')
    json.dump(rep, open(os.path.join(RES, 'stage295_ceiling_report.json'), 'w'),
              indent=1, default=float)
    print(json.dumps({k: rep[k] for k in ('inputs', 'validation', 'headline_case')},
                     indent=1, default=float)[:2600])
    print('\nPUSH COMPLETE ROWS')
    for x in sorted(done, key=lambda x: (x['reference_set_size'], x['name'])):
        print('  |R|=%-6d cap=%-3d min=%-3s %7.1fs faces=%-10d bases=%-9d minimality=%s  %s'
              % (x['reference_set_size'], x['cap'], x['minimum'], x['seconds'],
                 x['mass_faces'], x['vertex_bases'], x['minimality_certified'], x['name'][:42]))
    print('\nPUSH FAILURES')
    for x in rep['push']['failures']:
        print('   batch=', x['batch'])
        print('  |R|=%-6d cap=%-3d %s %s' % (x['reference_set_size'] or -1, x['cap'],
                                             x['error'], x['name'][:50]))


if __name__ == '__main__':
    main()
