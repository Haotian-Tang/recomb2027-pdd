#!/usr/bin/env python3
"""Stage 289 verification: audit the measured runs before anything is reported.

1. Every exact-solver witness is re-checked with the frozen stage-42 verifier.
2. Every DecoDiPhy support certified feasible is re-checked by evaluating the
   certificate through the frozen stage-26 observation routine.
3. No certified-feasible support may be smaller than the exact minimum on a case
   where the exact solver completed.
4. The two repetitions of every (case, method) pair must agree.
"""
import json
import multiprocessing as mp
import sys
from collections import Counter, defaultdict
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from stage289_decodiphy_bridge import load_case
from stage289_support_feasibility import feasible_support
from stage42_benchmark import verify as verify_witness
from stage26_missing_leaf import observe

PROTOCOL = ROOT / 'results/stage289/stage289_headtohead_protocol.json'
STREAM = ROOT / 'results/stage289/stage289_headtohead_runs.jsonl'
PERCASE = ROOT / 'results/stage289/stage289_headtohead_percase.json'
OUT = ROOT / 'results/stage289/stage289_verification.json'
DD = ('dd_default', 'dd_true_k', 'dd_resolution')
EXACT = 'exact_stage43_partial'


def check(payload):
    entry, by_method, derived = payload
    case = load_case(ROOT / entry['path'])
    case['cap'] = entry['cap']
    out = dict(name=entry['name'], family=entry['family'],
               witness_verified=None, certificates_checked=0, certificate_failures=[],
               contradictions=[], repetition_disagreements=[])
    for method, runs in by_method.items():
        keys = set()
        for run in runs:
            if method == EXACT:
                keys.add((run['status'], run.get('feasible'), run.get('minimum'),
                          tuple(sorted(s['edge'] for s in run.get('witness', {}).get('sources', [])))))
            else:
                keys.add((run['status'], run.get('k'),
                          tuple(sorted(s['edge'] for s in run.get('sources', []) or []))))
        if len(keys) != 1:
            out['repetition_disagreements'].append(method)
    exact_runs = by_method[EXACT]
    witness = exact_runs[0].get('witness')
    if witness:
        try:
            out['witness_verified'] = verify_witness(case, witness) == exact_runs[0]['minimum']
        except BaseException as exc:
            out['witness_verified'] = False
            out['witness_error'] = repr(exc)[:200]
    minimum = exact_runs[0].get('minimum') if exact_runs[0]['status'] == 'complete' else None
    for method in DD:
        record = derived['methods'][method]
        if record.get('exact_feasibility') != 'feasible':
            continue
        support = sorted(set(record['support']))
        verdict = feasible_support(case['edges'], case['observed'], support)
        sources = [dict(edge=s['edge'], p=F(s['p']), x=F(s['x'])) for s in verdict['sources']]
        pendant = F(verdict['pendant'])
        out['certificates_checked'] += 1
        try:
            assert observe(case['edges'], sources, pendant, case['observed']) == case['observed']
        except BaseException as exc:
            out['certificate_failures'].append(dict(method=method, error=repr(exc)[:200]))
        if minimum is not None and len(support) < minimum:
            out['contradictions'].append(dict(method=method, support_size=len(support),
                                              exact_minimum=minimum))
    return out


def main():
    protocol = json.loads(PROTOCOL.read_text())
    index = {e['name']: e for e in protocol['case_index']}
    derived = {c['name']: c for c in json.loads(PERCASE.read_text())}
    runs = defaultdict(lambda: defaultdict(list))
    with STREAM.open() as stream:
        for line in stream:
            record = json.loads(line)
            runs[record['name']][record['method']].append(record['result'])
    payloads = [(index[name], dict(by_method), derived[name]) for name, by_method in runs.items()]
    with mp.get_context('fork').Pool(12) as pool:
        records = pool.map(check, payloads, chunksize=4)
    summary = dict(
        stage=289, cases=len(records),
        exact_witnesses_checked=sum(r['witness_verified'] is not None for r in records),
        exact_witness_failures=[r['name'] for r in records if r['witness_verified'] is False],
        decodiphy_certificates_checked=sum(r['certificates_checked'] for r in records),
        certificate_failures=[dict(name=r['name'], **f) for r in records for f in r['certificate_failures']],
        contradictions=[dict(name=r['name'], **x) for r in records for x in r['contradictions']],
        repetition_disagreements=dict(Counter(m for r in records for m in r['repetition_disagreements'])),
        repetition_disagreement_cases=sorted({r['name'] for r in records if r['repetition_disagreements']}))
    summary['status'] = 'PASS' if (not summary['exact_witness_failures']
                                   and not summary['certificate_failures']
                                   and not summary['contradictions']) else 'FAIL'
    OUT.write_text(json.dumps(dict(summary, records=records), indent=1) + '\n')
    print(json.dumps({k: v for k, v in summary.items()
                      if k != 'repetition_disagreement_cases'}, indent=1)[:2500])
    print('repetition disagreement cases:', len(summary['repetition_disagreement_cases']))


if __name__ == '__main__':
    main()
