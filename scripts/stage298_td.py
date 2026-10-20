#!/usr/bin/env python3
"""Stage 298 step 3: Theta_1 on the 29 real travellers'-diarrhoea stool samples.

The backbone is |E| = 21,147 edges over |R| = 10,575 observed reference leaves.
The frozen sweep `theta_min_over_supports` wants the whole |V| x |R| leaf
distance matrix (2.2e8 rationals), so the sweep is run here edge by edge, each
support's two columns built with the FROZEN BFS `stage288_common.bfs_from`
called only from that edge's two endpoints, and each support's theta decided by
the FROZEN `stage296_common.theta_support`.  Nothing about the decision is new
code; only the loop that feeds it is.

In parallel the stage-298 k = 1 closed form is evaluated for every edge of
every sample and compared to the frozen LP by exact rational equality, so all
29 x 21,147 = 613,263 edges are decided twice by independent implementations.

The attaining witness is re-verified with the frozen
min(d(r,u)+xL, d(r,v)+(1-x)L) verifier on the original tree.
"""
import json
import os
import sys
import time
from fractions import Fraction as F

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
import stage288_common as common            # noqa: E402
import stage296_common as s296              # noqa: E402
import stage298_common as s298              # noqa: E402

OUT = s298.OUT
PROTOCOL = OUT / 'stage298_protocol.json'


def worker(spec):
    """Exact theta({e}) for a contiguous chunk of edges, frozen LP per edge."""
    edges, leaves, dvals = spec['edges'], spec['leaves'], spec['dvals']
    lo, hi, closed = spec['lo'], spec['hi'], spec['closed']
    t0 = time.perf_counter()
    thetas, mism, best = [], [], None
    for i in range(lo, hi):
        ends = sorted({edges[i][0], edges[i][1]})
        dist = common.bfs_from(edges, ends)
        u, v, _L = edges[i]
        cols = (tuple(dist[u][r] for r in leaves), tuple(dist[v][r] for r in leaves))
        r = s296.theta_support(cols, dvals, 1)
        if r['status'] != 'optimal':
            return dict(status='lp_' + r['status'], edge=i)
        thetas.append(r['value'])
        if closed[i - lo] != r['value']:
            mism.append(i)
        if best is None or r['value'] < thetas[best - lo]:
            best, bz = i, r['z']
    el = time.perf_counter() - t0
    return dict(status='complete', lo=lo, hi=hi, n=hi - lo,
                chunk_min=thetas[best - lo], chunk_argmin=best, z=bz,
                closed_form_mismatches=len(mism), mismatch_examples=mism[:5],
                seconds_per_support=el / max(1, hi - lo),
                theta_floats=[float(t) for t in thetas])


def main():
    protocol = json.loads(PROTOCOL.read_text())
    stream = OUT / 'stage298_td_runs.jsonl'
    summary_path = OUT / 'stage298_td_summary.json'
    assert not stream.exists(), 'Never overwrite frozen evidence'
    chunk = protocol['td_chunk_edges']
    out_rows = []
    for row in protocol['td_cases']:
        path = os.path.join(ROOT, row['path'])
        assert common.sha(path) == row['sha256'], 'input hash mismatch %s' % row['path']
        case = json.loads(open(path).read())
        edges, leaves, dvals = s298.case_arrays(case)
        obs = {int(r): F(d) for r, d in case['observed'].items()}
        med = s298.median_observed(dvals)
        t0 = time.perf_counter()
        cf = s298.theta1_all_edges(edges, leaves, dvals)
        t_cf = time.perf_counter() - t0
        jobs = []
        for lo in range(0, len(edges), chunk):
            hi = min(lo + chunk, len(edges))
            jobs.append(dict(
                meta=dict(name=row['name'], lo=lo, hi=hi),
                spec=dict(edges=edges, leaves=leaves, dvals=dvals, lo=lo, hi=hi,
                          closed=cf['theta'][lo:hi])))
        print('%s: closed form %.1f s, Theta_1(closed) = %s, %d chunks'
              % (row['name'], t_cf, float(cf['theta1']), len(jobs)), flush=True)
        t0 = time.time()
        res = common.run_pool(jobs, worker, protocol['wall_seconds']['td'],
                              protocol['address_space_bytes'], protocol['workers'],
                              stream, progress=10)
        el = time.time() - t0
        mins = [(F(r['result']['chunk_min']), r['result']['chunk_argmin'])
                for r in res if r['result'].get('status') == 'complete']
        incomplete = [r['meta'] if 'meta' in r else dict(lo=r['lo'], hi=r['hi'])
                      for r in res if r['result'].get('status') != 'complete']
        mismatches = sum(r['result'].get('closed_form_mismatches', 0) for r in res
                         if r['result'].get('status') == 'complete')
        covered = sum(r['result']['n'] for r in res if r['result'].get('status') == 'complete')
        frozen_theta1, frozen_arg = (min(mins) if mins else (None, None))
        # frozen witness at the attaining edge, verified independently
        ver = {}
        if frozen_arg is not None:
            fz = s298.frozen_theta_support(edges, leaves, dvals, (frozen_arg,))
            dev = s298.verify_witness_exact(edges, (frozen_arg,), fz['z'], 1, obs, leaves)
            srcs, y = s296.witness_sources((frozen_arg,), fz['z'], 1)
            ver = dict(frozen_lp_value=fz['value'],
                       independent_deviation=dev,
                       witness_agrees=(dev == fz['value'] == frozen_theta1),
                       witness=dict(pendant=y, sources=[dict(edge=s['edge'], p=s['p'],
                                                             x=s['x']) for s in srcs]))
        u, v, L = edges[frozen_arg] if frozen_arg is not None else (None, None, None)
        deg = {}
        for a, b, _l in edges:
            deg[a] = deg.get(a, 0) + 1
            deg[b] = deg.get(b, 0) + 1
        terminal_ends = [x for x in (u, v) if deg.get(x) == 1]
        out_rows.append(dict(
            name=row['name'], sample=row['sample'], group=row['group'],
            n_edges=len(edges), n_leaves=len(leaves),
            median_observed=med,
            theta1_closed_form=cf['theta1'], theta1_closed_argmin=cf['best'],
            theta1=frozen_theta1, argmin_edge=frozen_arg,
            theta1_over_median_observed=(float(frozen_theta1 / med)
                                         if frozen_theta1 is not None else None),
            argmin_endpoints=[u, v], argmin_edge_length=L,
            argmin_x=(cf['x'][frozen_arg] if frozen_arg is not None else None),
            argmin_is_terminal=bool(terminal_ends),
            argmin_terminal_leaf_labels=[case['leaf_labels'].get(str(x))
                                         for x in terminal_ends],
            closed_form_mismatches=mismatches, edges_covered=covered,
            incomplete_chunks=incomplete,
            closed_form_seconds=t_cf, frozen_sweep_seconds=el,
            verification=ver))
        print('  frozen Theta_1 = %s at edge %s; closed-form mismatches %d; %.0f s'
              % (float(frozen_theta1) if frozen_theta1 is not None else None,
                 frozen_arg, mismatches, el), flush=True)
        summary_path.write_text(json.dumps(common.serial(dict(
            stage=298, step='travelers_diarrhea', rows=out_rows)), indent=1))
    print('wrote', summary_path)


if __name__ == '__main__':
    main()

