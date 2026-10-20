#!/usr/bin/env python3
"""Stage292 pilot timing probe (no measurement is recorded from this script)."""
import sys, json, time
sys.path.insert(0, 'scripts'); sys.dont_write_bytecode = True
from pathlib import Path
from stage290_erasure_core import load_bank_case, support_of, build_struct, removal_order, CaseTables
from stage292_all_minimisers import all_minimum_supports

idx = json.loads(Path('data/stage287_instances/bank/index.json').read_text())['cases']
for K, fam in ((5, 'exact_10k'), (5, 'exact_biotrees'), (7, 'exact_biotrees'),
               (10, 'exact_10k'), (10, 'exact_biotrees'), (12, 'exact_16k')):
    rows = [x for x in idx if x['family'] == fam and x['K_generating'] == K]
    rows.sort(key=lambda z: z['reference_set_size'])
    picks = sorted({0, len(rows) // 2, len(rows) - 1})
    for pick in picks:
        spec = rows[pick]; case = load_bank_case(Path(spec['path']))
        sup = support_of(case); obs = case['observed']
        t0 = time.perf_counter()
        r = all_minimum_supports(case['edges'], obs, cap=K, seconds=300)
        t1 = time.perf_counter()
        print('%-16s K=%2d R=%5d E=%5d e=0: %-14s min=%s n=%s %.2fs cfg=%s orc=%s true_ok=%s' % (
            fam, K, len(obs), len(case['edges']), r['status'], r['minimum'], r['n_supports'],
            t1 - t0, r['statistics']['configurations'], r['statistics']['oracle_calls'],
            r['n_supports'] == 1 and r['supports'] and r['supports'][0] == sup), flush=True)
        if r['status'] != 'complete':
            continue
        tables = CaseTables(case['edges']); st = build_struct(case['edges'], sup, sorted(obs))
        order = removal_order(case, st, 'adversarial_nearest', 0, tables=tables)
        for e in (1, 2, 4, 8):
            if e > len(obs) - 3:
                break
            rem = set(order[:e]); sub = {k: v for k, v in obs.items() if k not in rem}
            t0 = time.perf_counter()
            rr = all_minimum_supports(case['edges'], sub, cap=K, seconds=300)
            t1 = time.perf_counter()
            print('      adv e=%d: %-14s min=%s n=%s %.2fs cfg=%s orc=%s' % (
                e, rr['status'], rr['minimum'], rr['n_supports'], t1 - t0,
                rr['statistics']['configurations'], rr['statistics']['oracle_calls']), flush=True)

