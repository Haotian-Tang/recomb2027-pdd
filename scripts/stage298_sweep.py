#!/usr/bin/env python3
"""Stage 298 step 2: exact Theta_1, Theta_2, Theta_3 on the small real cases.

Collections: krepp_estimated_10k (296 real ANI cases from the stage-287 bank)
and krepp_placement_10k (60 real placement-distance cases converted in step 1).

Everything that decides a number is frozen stage-296 code:
`theta_support` (exact rational two-phase simplex with row generation),
`theta_min_over_supports` (exhaustive sweep with exact incumbent pruning) and
`independent_deviation` (recomputes the leaf equations on the original tree
with the min(d(r,u)+xL, d(r,v)+(1-x)L) form).

k = 1 is run support by support with NO cutoff, so every edge gets its exact
theta({e}); those values are compared against stage 298's k=1 closed form by
exact rational equality, which is what licenses the closed form on the
travellers'-diarrhoea backbone where the frozen sweep's leaf-distance matrix
would be 2.2e8 rationals.

k = 2 and k = 3 use the frozen sweep with cutoff = Theta_{k-1}: Theta_k is
non-increasing in k, so any support whose exact lower bound reaches
Theta_{k-1} cannot improve it.  When nothing is strictly below, Theta_k =
Theta_{k-1} exactly and the witness is the (k-1) witness padded with a
zero-mass edge, which the closed system allows.

A (case, k) cell is RUN when C(|E|, k) <= the registered support cap and
CENSORED otherwise; censored cells are recorded with their support count and
the measured cost per support, never estimated away.
"""
import json
import os
import sys
import time
from fractions import Fraction as F
from itertools import combinations
from math import comb

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
import stage288_common as common            # noqa: E402
import stage296_common as s296              # noqa: E402
import stage298_common as s298              # noqa: E402

OUT = s298.OUT
PROTOCOL = OUT / 'stage298_protocol.json'


def worker(spec):
    case, k, cutoff = spec['case'], spec['k'], spec['cutoff']
    edges, leaves, dvals = s298.case_arrays(case)
    obs = {int(r): F(d) for r, d in case['observed'].items()}
    t0 = time.perf_counter()
    Dv = s296.leaf_vectors(edges, leaves)
    t_prep = time.perf_counter() - t0
    out = {'status': 'error', 'k': k}
    if k == 1:
        t0 = time.perf_counter()
        thetas, best = [], None
        for i in range(len(edges)):
            r = s296.theta_support(s296.support_cols(Dv, edges, (i,)), dvals, 1)
            if r['status'] != 'optimal':
                return dict(status='lp_' + r['status'], edge=i, k=1)
            thetas.append(r['value'])
            if best is None or r['value'] < thetas[best]:
                best, bz = i, r['z']
        t_sweep = time.perf_counter() - t0
        cf = s298.theta1_all_edges(edges, leaves, dvals)
        mism = [i for i in range(len(edges)) if cf['theta'][i] != thetas[i]]
        out = dict(status='complete', k=1, theta=thetas[best], support=[best],
                   z=bz, visited=len(edges), pruned=0,
                   strictly_below_previous=True,
                   closed_form_theta1=cf['theta1'],
                   closed_form_argmin=cf['best'],
                   closed_form_mismatches=len(mism),
                   closed_form_mismatch_examples=mism[:5],
                   per_edge_theta=[str(t) for t in thetas],
                   seconds_per_support=t_sweep / max(1, len(edges)),
                   prep_seconds=t_prep, sweep_seconds=t_sweep)
    else:
        t0 = time.perf_counter()
        r = s296.theta_min_over_supports(Dv, edges, dvals, k,
                                         cutoff=F(cutoff) if cutoff is not None else None)
        t_sweep = time.perf_counter() - t0
        if r['status'] != 'complete':
            return dict(status='sweep_' + r['status'], k=k,
                        visited=r.get('visited'), sweep_seconds=t_sweep)
        out = dict(status='complete', k=k, theta=r['theta'],
                   support=(list(r['support']) if r['support'] is not None else None),
                   z=r['z'], visited=r['visited'], pruned=r['pruned'],
                   strictly_below_previous=r['support'] is not None,
                   seconds_per_support=t_sweep / max(1, r['visited']),
                   prep_seconds=t_prep, sweep_seconds=t_sweep)
    # independent witness verification with the frozen min(...) form
    if out.get('support') is not None and out.get('z') is not None:
        dev = s296.independent_deviation(edges, tuple(out['support']), out['z'],
                                         k, obs, leaves=leaves)
        out['independent_deviation'] = dev
        out['witness_agrees'] = (dev is not None and dev == out['theta'])
        srcs, y = s296.witness_sources(tuple(out['support']), out['z'], k)
        out['witness'] = dict(
            pendant=y, sources=[dict(edge=s['edge'], p=s['p'], x=s['x']) for s in srcs])
    out['median_observed'] = s298.median_observed(dvals)
    out['n_edges'] = len(edges)
    out['n_leaves'] = len(leaves)
    return out


def cells(protocol, k):
    jobs = []
    for row in protocol['small_cases']:
        n = comb(row['n_edges'], k)
        if n > protocol['support_cap'][str(k)]:
            continue
        jobs.append(row)
    return jobs


def main():
    protocol = json.loads(PROTOCOL.read_text())
    stream = OUT / 'stage298_sweep_runs.jsonl'
    assert not stream.exists(), 'Never overwrite frozen evidence'
    cases = {}
    for row in protocol['small_cases']:
        path = os.path.join(ROOT, row['path'])
        assert common.sha(path) == row['sha256'], 'input hash mismatch %s' % row['path']
        cases[row['name']] = json.loads(open(path).read())
    theta_prev = {}
    censored = []
    for k in (1, 2, 3):
        rows = cells(protocol, k)
        skipped = [r for r in protocol['small_cases'] if r not in rows]
        for r in skipped:
            censored.append(dict(name=r['name'], family=r['family'], k=k,
                                 n_edges=r['n_edges'],
                                 n_supports=comb(r['n_edges'], k),
                                 reason='support_cap'))
        jobs = []
        for row in rows:
            prev = theta_prev.get(row['name'])
            if k > 1 and prev is None:
                censored.append(dict(name=row['name'], family=row['family'], k=k,
                                     n_edges=row['n_edges'],
                                     n_supports=comb(row['n_edges'], k),
                                     reason='previous_k_unavailable'))
                continue
            jobs.append(dict(
                meta=dict(name=row['name'], family=row['family'], k=k,
                          K_generating=row['K_generating'],
                          reference_set_size=row['reference_set_size'],
                          n_edges=row['n_edges'],
                          n_supports=comb(row['n_edges'], k),
                          cutoff=(str(prev) if prev is not None else None)),
                spec=dict(case=cases[row['name']], k=k,
                          cutoff=(str(prev) if prev is not None else None))))
        print('k=%d: %d runs (%d censored so far)' % (k, len(jobs), len(censored)), flush=True)
        t0 = time.time()
        res = common.run_pool(jobs, worker, protocol['wall_seconds'][str(k)],
                              protocol['address_space_bytes'], protocol['workers'],
                              stream, progress=25)
        print('  k=%d done in %.0f s' % (k, time.time() - t0), flush=True)
        for rec in res:
            a = rec['result']
            if a.get('status') == 'complete' and a.get('theta') is not None:
                theta_prev[rec['name']] = F(a['theta'])
            else:
                theta_prev.pop(rec['name'], None)
                censored.append(dict(name=rec['name'], family=rec['family'], k=k,
                                     n_edges=rec['n_edges'],
                                     n_supports=rec['n_supports'],
                                     reason=a.get('status', 'unknown'),
                                     detail=a.get('reason')))
    (OUT / 'stage298_sweep_censored.json').write_text(
        json.dumps(common.serial(dict(censored=censored, n=len(censored))), indent=1))
    print('censored cells:', len(censored))


if __name__ == '__main__':
    main()

