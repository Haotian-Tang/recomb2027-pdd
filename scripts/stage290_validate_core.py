#!/usr/bin/env python3
"""Stage290 step 1: verify the reused Thread-F machinery before trusting it.

Three independent comparisons on the seven frozen Stage196 fixtures:
  (a) Stage290 scaled DP   vs  the saved Stage196 result records;
  (b) Stage290 scaled DP   vs  the frozen Stage196 solve_dp run live;
  (c) Stage290 scaled DP   vs  the frozen Stage196 exhaustive separator
      enumeration (brute force over every vertex subset of size < s).
Also checks the Menger/rank identity rank(leaf-report matrix) = s + 1 on the
same fixtures with unit edge lengths and with random positive rational lengths.
"""
import json
import random
import sys
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
sys.dont_write_bytecode = True

import stage196_public_erasure_validation as S196       # noqa: E402
from stage290_erasure_core import (build_struct, erasure_analysis,              # noqa: E402
                                   fixed_support_rank, encode, PRIME, PRIME2)


class NoBudget(object):
    def check(self):
        pass

    def charge(self, key, maximum=None):
        pass

    counts = __import__('collections').Counter()


def to_edges(incidence, lengths=None):
    labels = list(dict.fromkeys(v for e in incidence for v in e))
    index = {v: i for i, v in enumerate(labels)}
    out = []
    for i, (a, b) in enumerate(incidence):
        L = F(1) if lengths is None else lengths[i]
        out.append((index[a], index[b], L))
    return out, index, labels


def main():
    cases = json.loads((ROOT / 'data/stage196_public_erasure_cases.json').read_text())['cases']
    saved = {r['id']: r for r in json.loads(
        (ROOT / 'results/stage196_public_erasure_validation.json').read_text())['records']}
    rng = random.Random(29012026)
    rows = []
    for case in cases:
        cid = case['id']
        incidence = case['graph_inputs']['incidence']
        support = case['graph_inputs']['support']
        edges, index, labels = to_edges(incidence)
        st = build_struct(edges, support, sorted(v for v in range(len(labels))
                                                 if sum(1 for e in edges if v in e[:2]) == 1))
        mine = erasure_analysis(st, want_witness=True)

        # (a) saved Stage196 record
        rec = saved[cid]['record']
        saved_s = rec['rank_certificate']['s']
        saved_tol = rec['thresholds']['information_erasure_tolerance']

        # (b) frozen Stage196 DP, run live on the same input
        budget = S196.Budget()
        model = S196.topology(incidence, support, budget)
        rank = S196.verify_rank_certificate(model, case['rank_certificate'], budget)
        raw = S196.solve_dp(model, rank['s'], budget)
        frozen_delta = raw['delta'] if raw['status'] == 'DP_CANDIDATE' else None

        # (c) frozen exhaustive separator enumeration
        brute = None
        if case['checker_only']['exhaustive_separator_check']:
            enum = S196.enumerate_separators(model, rank['s'], S196.Budget())
            brute = enum['minimum'] if enum['status'].startswith('COMPLETE') else None

        # rank identity, unit lengths and random positive rational lengths
        rank_unit = fixed_support_rank(st, PRIME)
        rand_edges = [(u, v, F(rng.randint(1, 97), rng.randint(1, 31)))
                      for u, v, _ in edges]
        st2 = build_struct(rand_edges, support, sorted(st['R']))
        rank_rand = fixed_support_rank(st2, PRIME)
        rank_rand2 = fixed_support_rank(st2, PRIME2)
        mine2 = erasure_analysis(st2, want_witness=False)

        rows.append(dict(
            id=cid, H=st['H'], P=len(st['P']), R=len(st['R']),
            stage290_s=mine['s'], stage196_saved_s=saved_s,
            stage290_delta=mine['delta'], stage196_live_delta=frozen_delta,
            brute_force_delta=brute,
            stage290_tolerance=mine['tolerance'], stage196_saved_tolerance=saved_tol,
            witness_checked=mine['witness_checked'],
            witness_C=mine['witness_C'], witness_D=mine['witness_D'],
            rank_unit=rank_unit, rank_random_lengths=rank_rand,
            rank_random_second_prime=rank_rand2,
            s_random_lengths=mine2['s'],
            delta_random_lengths=mine2['delta'],
            rank_identity_unit=(rank_unit == mine['s'] + 1),
            rank_identity_random=(rank_rand == mine2['s'] + 1 == rank_rand2),
            metric_invariance=(mine2['s'] == mine['s'] and mine2['delta'] == mine['delta']),
            agree_saved=(mine['s'] == saved_s and mine['tolerance'] == saved_tol),
            agree_live=(frozen_delta == mine['delta']),
            agree_brute=(brute is None or brute == mine['delta'])))
    summary = dict(
        stage=290, check='reuse_verification_against_frozen_stage196',
        n_cases=len(rows),
        all_agree_saved=all(r['agree_saved'] for r in rows),
        all_agree_live=all(r['agree_live'] for r in rows),
        all_agree_brute=all(r['agree_brute'] for r in rows),
        brute_forced=sum(r['brute_force_delta'] is not None for r in rows),
        all_witnesses_checked=all(r['witness_checked'] for r in rows),
        all_rank_identity_unit=all(r['rank_identity_unit'] for r in rows),
        all_rank_identity_random=all(r['rank_identity_random'] for r in rows),
        all_metric_invariant=all(r['metric_invariance'] for r in rows),
        rows=rows)
    out = ROOT / 'results/stage290_reuse_verification.json'
    with out.open('x') as fh:
        json.dump(encode(summary), fh, indent=2)
        fh.write('\n')
    print(json.dumps({k: v for k, v in summary.items() if k != 'rows'}, indent=2))


if __name__ == '__main__':
    main()
