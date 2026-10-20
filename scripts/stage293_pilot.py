#!/usr/bin/env python3
"""Stage 293 pilot: exercise every arm once before the production run."""
import os, sys, json, time
for _k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[_k] = '1'
from pathlib import Path
from fractions import Fraction as F
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import stage288_common as common
import stage288_incompleteness as inc
import stage293_abundance as ab

which = sys.argv[1] if len(sys.argv) > 1 else 'all'
entries = [e for e in inc.panel() if e['K_generating'] == 5 and e['reference_set_size'] < 80]
entry = entries[0]
case = common.load_case(entry)
draws = inc.removal_draws(case, '0.50')
draw = [d for d in draws if d['regime'] == 'nearest'][0]
retained = ab.retained_observations(case, draw['removed'])
print('case', entry['name'], 'N', entry['reference_set_size'], 'K', entry['K_generating'],
      'retained', len(retained), 'removed', draw['n_removed'], flush=True)
truth_support = sorted(int(s['edge']) for s in case['generating_sources'])

if which in ('all', 'exact'):
    t = time.perf_counter()
    a = ab.exact_worker(dict(edges=case['edges'], observed=retained, cap=case['cap']))
    print('EXACT status=%s min=%s support=%s  %.2fs' % (
        a['status'], a.get('minimum'),
        sorted({int(s['edge']) for s in a['witness']['sources']}) if 'witness' in a else None,
        time.perf_counter()-t), flush=True)
    print('   true support', truth_support, flush=True)

if which in ('all', 'fibre'):
    t = time.perf_counter()
    f = ab.fibre_worker(dict(edges=case['edges'], observed=retained, support=truth_support,
                             truth=case['generating_sources']))
    print('FIBRE_TRUE status=%s free=%s p_sup_range=%s x_sup_range=%s maxdev_p=%s  %.2fs' % (
        f['status'], f.get('free_parameters'),
        (float(F(f['p_sup_range'])) if 'p_sup_range' in f else None),
        (float(F(f['x_sup_range'])) if 'x_sup_range' in f else None),
        (float(F(f['max_p_deviation_from_truth'])) if 'max_p_deviation_from_truth' in f else None),
        time.perf_counter()-t), flush=True)
    print('   truth_inside_bounds =', f.get('truth_inside_bounds'), flush=True)

if which in ('all', 'dd'):
    for method in ('dd_default', 'dd_true_k'):
        t = time.perf_counter()
        try:
            d = ab.dd_worker(dict(edges=case['edges'], observed=retained, method=method,
                                  k_true=entry['K_generating'],
                                  generating_sources=case['generating_sources']))
            print('%s k=%s loss=%.3e stop=%s support=%s  %.2fs' % (
                method, d['k'], d['loss'], d['stop_reason'],
                sorted({s['edge'] for s in d['sources']}), time.perf_counter()-t), flush=True)
        except BaseException as exc:
            print(method, 'FAILED', repr(exc)[:400], flush=True)
