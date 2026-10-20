#!/usr/bin/env python3
"""Stage292: characterise the K = 10 ceiling on the smallest K = 10 bank cases.

Separates the cost of the frozen minimum pass (solve_partial) from the cost of
the all-minimisers pass, at e = 0 only. Diagnostic; no threshold is measured.
"""
import json, sys, time
sys.path.insert(0, 'scripts'); sys.dont_write_bytecode = True
from pathlib import Path
from fractions import Fraction as F
from stage43_partial_fpt import solve_partial
from stage290_erasure_core import load_bank_case, support_of
from stage292_all_minimisers import all_minimum_supports

idx = json.loads(Path('data/stage287_instances/bank/index.json').read_text())['cases']
rows = [x for x in idx if x['K_generating'] in (10, 12)]
rows.sort(key=lambda z: z['reference_set_size'])
for spec in rows[:4] + rows[-1:]:
    case = load_bank_case(Path(spec['path'])); sup = support_of(case)
    K = len(sup); obs = case['observed']
    E = [(u, v, F(L)) for u, v, L in case['edges']]
    t0 = time.perf_counter()
    fr = solve_partial(E, obs, cap=K, seconds=3600, max_states=None,
                       max_bases=None, max_faces=None)
    t1 = time.perf_counter()
    print('%-44s K=%2d R=%5d frozen_solve_partial %-10s min=%s %.1fs' % (
        spec['name'][:44], K, len(obs), fr.get('status'),
        fr.get('minimum_sources_within_cap'), t1 - t0), flush=True)
    t2 = time.perf_counter()
    r = all_minimum_supports(case['edges'], obs, cap=K, seconds=3600,
                             kstar=fr.get('minimum_sources_within_cap'))
    t3 = time.perf_counter()
    print('    all_minimisers %-14s min=%s n=%s %.1fs  patterns=%s cfg=%s oracle=%s' % (
        r['status'], r['minimum'], r['n_supports'], t3 - t2,
        r['statistics']['feasible_patterns'], r['statistics']['configurations'],
        r['statistics']['oracle_calls']), flush=True)

