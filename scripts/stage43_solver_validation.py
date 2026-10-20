#!/usr/bin/env python3
"""Exact adversarial differential checks; no timing comparison or fixture selection."""
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import random
from stage27_general_fpt import Budget, dot
from stage41_partial_fpt import strict_vertex as original
import stage43_partial_fpt as revised
from stage42_benchmark import fixtures, verify
ROOT = Path(__file__).resolve().parents[1]


def compare(rows, dim):
    before = original(rows, dim, Budget())
    after = revised.strict_vertex(rows, dim, Budget())
    assert (before is None) == (after is None), (rows, before, after)
    if after is not None:
        assert after[-1] > 0 and all(dot(a, after) <= rhs for a, rhs in rows)


def run():
    rng = random.Random(431015)
    generic = 0
    for dimension in (1, 2, 3):
        for trial in range(20):
            rows = []
            for j in range(dimension):
                row = tuple(F(i == j) for i in range(dimension))
                rows.extend([(row, F(1)), (tuple(-x for x in row), F(0))])
            for _ in range(3):
                rows.append((tuple(F(rng.randrange(-2, 3)) for _ in range(dimension)), F(rng.randrange(-1, 4), 2)))
            if trial % 2 == 0:
                a = tuple(F(rng.randrange(-2, 3)) for _ in range(dimension))
                b = F(rng.randrange(-1, 2), 2)
                rows.extend([(a, b), (tuple(-2*x for x in a), -2*b)])
            if trial % 5 == 0:
                rows.extend([((F(0),) * dimension, F(0)), (rows[0][0], F(2))])
            compare(rows, dimension)
            generic += 1
    actual_faces = 0
    real_solver = revised.strict_vertex
    def differential(rows, dimension, budget):
        nonlocal actual_faces
        answer = real_solver(rows, dimension, budget)
        reference = original(rows, dimension, Budget())
        assert (answer is None) == (reference is None), (rows, answer, reference)
        actual_faces += 1
        return answer
    revised.strict_vertex = differential
    witnesses = 0
    answers = []
    try:
        for case in fixtures():
            answer = revised.solve_partial(case['edges'], case['observed'], cap=case['cap'],
                                           seconds=None, max_states=None, max_faces=None, max_bases=None)
            assert answer['status'] == 'complete'
            if answer['feasible']:
                assert verify(case, answer['original_witness']) == answer['minimum_sources_within_cap']
                witnesses += 1
            if 'expected_feasible' in case:
                assert answer['feasible'] == case['expected_feasible']
                if answer['feasible']:
                    assert answer['minimum_sources_within_cap'] == case['expected_minimum']
            answers.append(dict(name=case['name'], feasible=answer['feasible'],
                                minimum=answer.get('minimum_sources_within_cap')))
    finally:
        revised.strict_vertex = real_solver
    out = dict(status='PASS', generic_bounded_lp_differential_cases=generic,
               encountered_local_lp_differential_cases=actual_faces, fixture_cases=len(answers),
               independently_verified_witnesses=witnesses, answers=answers,
               solver_sha256=hashlib.sha256((ROOT/'scripts/stage43_partial_fpt.py').read_bytes()).hexdigest(),
               limitation='Exact LP equivalence and witness checks, not a timing measurement or empirical novelty claim')
    path = ROOT/'results/stage43_solver_validation.json'
    path.write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k != 'answers'}, indent=2))


if __name__ == '__main__':
    run()
