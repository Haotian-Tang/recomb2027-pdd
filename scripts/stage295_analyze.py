#!/usr/bin/env python3
"""Stage295 analysis: where the frozen minimum pass spends its time, how its
counters grow in K, and what the reduced-face prototype buys.

Reads only Stage295 raw records.  Produces results/stage295_ceiling_analysis.json
and a plain-text report on stdout.
"""
import json, os, sys, glob
from collections import defaultdict
import numpy as np
from scipy.optimize import nnls

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RES = os.path.join(ROOT, 'results')


def load_runs(pattern):
    rows = []
    for p in sorted(glob.glob(pattern)):
        if p.endswith('all_runs.json') or p.endswith('all_validation.json'):
            continue
        try:
            rows.append(json.load(open(p)))
        except Exception:
            pass
    return rows


def st(run, key, default=0):
    return (run.get('statistics') or {}).get(key, default)


# ---------------------------------------------------------------- phase cost
def phase_regression(runs, tag):
    """Fit elapsed = c0 + alpha*mass_faces + beta*vertex_bases on prefix runs.

    The two prefix ladders (max_faces and max_bases) move work along different
    directions in the (faces, bases) plane, which makes the two marginal costs
    separately identifiable.  Non-negative least squares; saturated duplicates
    are collapsed so repeated time-limited points cannot dominate the fit.
    """
    pts = {}
    for r in runs:
        if not r.get('ok') or r['row'].get('tag') not in tag:
            continue
        f, b, t = st(r, 'mass_faces'), st(r, 'vertex_bases'), r.get('elapsed')
        if not t or f == 0:
            continue
        key = (round(f, -max(0, len(str(f)) - 3)), round(b, -max(0, len(str(b)) - 3)))
        if key not in pts or t < pts[key][2]:
            pts[key] = (f, b, t)
    pts = sorted(pts.values())
    if len(pts) < 3:
        return None
    A = np.array([[1.0, f, b] for f, b, t in pts])
    y = np.array([t for f, b, t in pts])
    coef, _ = nnls(A, y)
    pred = A @ coef
    ss_res = float(((y - pred) ** 2).sum())
    ss_tot = float(((y - y.mean()) ** 2).sum())
    return {'n_points': len(pts), 'const_s': coef[0],
            'seconds_per_face': coef[1], 'seconds_per_basis': coef[2],
            'r2': 1 - ss_res / ss_tot if ss_tot > 0 else None,
            'points': [{'mass_faces': f, 'vertex_bases': b, 'elapsed_s': t}
                       for f, b, t in pts]}


def lp_breakdown(runs, tag):
    """Counter shares at stop for the deepest run of a tag."""
    sel = [r for r in runs if r.get('ok') and r['row'].get('tag') == tag]
    if not sel:
        return None
    r = max(sel, key=lambda r: st(r, 'mass_faces'))
    faces = st(r, 'mass_faces')
    lp = st(r, 'local_lp_calls')
    enum = st(r, 'enumerated_lp_calls')
    return {'tag': tag, 'elapsed_s': r.get('elapsed'), 'stop_reason': r.get('reason'),
            'mass_faces': faces, 'deficit_sets': st(r, 'deficit_sets'),
            'growth_states': st(r, 'growth_states'),
            'cover_feasible_faces': st(r, 'cover_feasible_faces'),
            'local_lp_calls': lp,
            'zero_dimensional_lp_calls': st(r, 'zero_dimensional_lp_calls'),
            'interval_lp_calls': st(r, 'interval_lp_calls'),
            'enumerated_lp_calls': enum,
            'vertex_bases': st(r, 'vertex_bases'),
            'strict_faces': st(r, 'strict_faces'),
            'face_reaches_lp_fraction': lp / faces if faces else None,
            'lp_is_enumerated_fraction': enum / lp if lp else None,
            'bases_per_enumerated_call': st(r, 'vertex_bases') / enum if enum else None,
            'bases_per_face': st(r, 'vertex_bases') / faces if faces else None}


FILE_PHASE = {
    'fractions.py': 'exact_rational_arithmetic',
    'abc': 'exact_rational_arithmetic',
    '<frozen abc>': 'exact_rational_arithmetic',
    'numbers.py': 'exact_rational_arithmetic',
}
NAME_PHASE = {
    ('stage41_partial_fpt.py', 'deficit_sets'): 'P1_deficit_supports',
    ('stage43_partial_fpt.py', 'solve_mixed'): 'P2_face_loop',
    ('stage41_partial_fpt.py', 'mixed_cover'): 'P3_cover_filter',
    ('stage41_partial_fpt.py', 'minimum_forest_cover'): 'P3_cover_filter',
    ('stage27_general_fpt.py', 'minimum_forest_cover'): 'P3_cover_filter',
    ('stage43_partial_fpt.py', '_parallel_reduce'): 'P4a_lp_row_reduction',
    ('stage43_partial_fpt.py', '_affine_equalities'): 'P4b_equality_elimination',
    ('stage43_partial_fpt.py', 'strict_vertex'): 'P4c_basis_enumeration',
    ('stage27_general_fpt.py', 'square_solve'): 'P4c_basis_enumeration',
    ('stage27_general_fpt.py', 'dot'): 'P4d_feasibility_checks',
    ('stage43_partial_fpt.py', 'lift'): 'P4d_feasibility_checks',
}


def classify(f):
    key = (f['file'], f['name'])
    if key in NAME_PHASE:
        return NAME_PHASE[key]
    if f['file'] in FILE_PHASE:
        return FILE_PHASE[f['file']]
    if f['file'] in ('~', '<string>') or f['name'].startswith('<'):
        return 'interpreter_builtins'
    return 'other_' + f['file']


def profile_shares(runs, tag):
    sel = [r for r in runs if r.get('ok') and r['row'].get('tag') == tag and r.get('profile')]
    if not sel:
        return None
    r = sel[0]
    total = r.get('profile_total_tt') or sum(f['tottime'] for f in r['profile'])
    buckets = defaultdict(float)
    cum = {}
    for f in r['profile']:
        buckets[classify(f)] += f['tottime']
        cum[(f['file'], f['name'])] = f['cumtime']
    return {'tag': tag, 'profiled_total_tt_s': total,
            'wall_s': r.get('elapsed'), 'mass_faces': st(r, 'mass_faces'),
            'vertex_bases': st(r, 'vertex_bases'),
            'shares': {k: v / total for k, v in sorted(buckets.items(), key=lambda kv: -kv[1])},
            'seconds': {k: v for k, v in sorted(buckets.items(), key=lambda kv: -kv[1])},
            'cumtime_key_frames': {k: cum.get(v) for k, v in {
                'strict_vertex': ('stage43_partial_fpt.py', 'strict_vertex'),
                '_parallel_reduce': ('stage43_partial_fpt.py', '_parallel_reduce'),
                '_affine_equalities': ('stage43_partial_fpt.py', '_affine_equalities'),
                'square_solve': ('stage27_general_fpt.py', 'square_solve'),
                'mixed_cover': ('stage41_partial_fpt.py', 'mixed_cover'),
                'deficit_sets': ('stage41_partial_fpt.py', 'deficit_sets'),
                'solve_mixed': ('stage43_partial_fpt.py', 'solve_mixed')}.items()},
            'top_functions': [{'name': f['name'], 'file': f['file'], 'ncalls': f['ncalls'],
                               'tottime': f['tottime'], 'cumtime': f['cumtime']}
                              for f in r['profile'][:15]]}


def growth(runs, census):
    """Realised growth in K of each counter, against the theoretical shape."""
    per_cap = defaultdict(lambda: defaultdict(list))
    for r in runs:
        tag = r['row'].get('tag', '')
        if not (r.get('ok') and tag.startswith('full_K') and r.get('status') == 'complete'):
            continue
        cap = r['row']['cap']
        for key in ('mass_faces', 'deficit_sets', 'growth_states', 'vertex_bases',
                    'cover_feasible_faces', 'local_lp_calls', 'enumerated_lp_calls'):
            per_cap[cap][key].append(st(r, key))
        per_cap[cap]['elapsed'].append(r.get('elapsed'))
    med = {cap: {k: float(np.median(v)) for k, v in d.items()} for cap, d in per_cap.items()}
    caps = sorted(med)
    fits = {}
    for key in ('mass_faces', 'deficit_sets', 'vertex_bases', 'elapsed'):
        xs = [c for c in caps if med[c].get(key, 0) > 0]
        if len(xs) < 3:
            continue
        ys = [np.log2(med[c][key]) for c in xs]
        slope, intercept = np.polyfit(xs, ys, 1)
        fits[key] = {'bits_per_K': float(slope), 'intercept_bits': float(intercept),
                     'caps_used': xs,
                     'predicted_at_K10': float(2 ** (slope * 10 + intercept)),
                     'per_K_factor': float(2 ** slope)}
    # exact K=10 face count from the static census, as a check on the fit
    cen = {}
    for row in census:
        if row.get('status') == 'complete' and row['cap'] == 10:
            cen.setdefault('mass_faces_exact_K10', []).append(row['mass_faces_total'])
            cen.setdefault('mass_faces_reduced_K10', []).append(row['mass_faces_reduced'])
    for k, v in list(cen.items()):
        cen[k] = {'median': float(np.median(v)), 'min': min(v), 'max': max(v), 'n': len(v)}
    # cap-ladder growth on one fixed case (isolates cap as the sole variable)
    ladder = defaultdict(dict)
    for row in census:
        if row.get('status') == 'complete':
            ladder[row['case']][row['cap']] = row
    cap_ladder = {}
    for case, byc in ladder.items():
        if len(byc) >= 5:
            caps2 = sorted(byc)
            cap_ladder[os.path.basename(case)] = {
                'caps': caps2,
                'mass_faces': [byc[k]['mass_faces_total'] for k in caps2],
                'mass_faces_reduced': [byc[k]['mass_faces_reduced'] for k in caps2],
                'deficit_supports': [byc[k]['deficit_supports'] for k in caps2],
                'max_A': [byc[k]['max_A_size'] for k in caps2],
                'theory_max_A': [byc[k]['theory_max_A'] for k in caps2],
                'faces_bits_per_K': float(np.polyfit(caps2,
                    [np.log2(byc[k]['mass_faces_total']) for k in caps2], 1)[0]),
                'reduced_bits_per_K': float(np.polyfit(caps2,
                    [np.log2(byc[k]['mass_faces_reduced']) for k in caps2], 1)[0]),
                'supports_bits_per_K': float(np.polyfit(caps2,
                    [np.log2(byc[k]['deficit_supports']) for k in caps2], 1)[0])}
    return {'median_counters_by_cap': med, 'log2_fits': fits,
            'census_exact_K10': cen, 'cap_ladder': cap_ladder}


def main():
    runs = load_runs(os.path.join(RES, 'stage295_profile_runs', 'run_*.json'))
    runs += load_runs(os.path.join(RES, 'stage295_profile_deep', 'run_*.json'))
    census = json.load(open(os.path.join(RES, 'stage295_static_census.json')))
    out = {
        'n_profile_runs': len(runs),
        'phase_cost_K10': phase_regression(runs, ('prefix_faces_K10', 'prefix_bases_K10')),
        'phase_cost_K7case_cap10': phase_regression(runs, ('prefix_faces_K7case_cap10',)),
        'lp_breakdown': {t: lp_breakdown(runs, t) for t in
                         ('prefix_faces_K10', 'prefix_bases_K10',
                          'prefix_faces_K7case_cap10', 'full_K10', 'full_K7')},
        'cprofile': {t: profile_shares(runs, t) for t in
                     ('cprofile_K10_deep', 'cprofile_K10_bases', 'cprofile_K10', 'cprofile_cap9', 'cprofile_cap8',
                      'cprofile_full_K7', 'cprofile_full_K5')},
        'growth': growth(runs, census),
    }
    # projected frozen runtime at K=10 from measured marginal costs
    reg = out['phase_cost_K10']
    bd = out['lp_breakdown']['prefix_faces_K10']
    if reg and bd:
        cen = out['growth']['census_exact_K10']
        faces = cen['mass_faces_exact_K10']['median']
        red = cen['mass_faces_reduced_K10']['median']
        bpf = bd['bases_per_face'] or 0
        per_face = reg['seconds_per_face'] + reg['seconds_per_basis'] * bpf
        out['projection_K10'] = {
            'exact_faces_median': faces, 'reduced_faces_median': red,
            'measured_bases_per_face': bpf,
            'seconds_per_face_total': per_face,
            'projected_frozen_seconds': faces * per_face,
            'projected_frozen_days': faces * per_face / 86400.0,
            'projected_reduced_seconds_same_unit_cost': red * per_face,
            'face_reduction_factor': faces / red if red else None,
            'note': ('Unit cost measured on budget-truncated prefix runs of the '
                     'smallest K=10 case; the projection assumes the sampled '
                     'per-face cost is representative of the whole search.')}
    json.dump(out, open(os.path.join(RES, 'stage295_ceiling_analysis.json'), 'w'),
              indent=1, default=float)
    print(json.dumps({k: out[k] for k in ('n_profile_runs', 'phase_cost_K10',
                                          'projection_K10')}, indent=1, default=float)[:3000])
    print('\nLP BREAKDOWN')
    for t, v in out['lp_breakdown'].items():
        if v:
            print(' ', t, json.dumps(v, default=float))
    print('\nGROWTH FITS')
    print(json.dumps(out['growth']['log2_fits'], indent=1, default=float))
    print('\nCAP LADDER')
    print(json.dumps(out['growth']['cap_ladder'], indent=1, default=float))
    print('\nCPROFILE SHARES')
    for t, v in out['cprofile'].items():
        if v:
            print(' ', t, json.dumps(v['shares'], default=float), 'wall', v['wall_s'])


if __name__ == '__main__':
    main()
