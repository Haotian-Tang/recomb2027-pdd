#!/usr/bin/env python3
"""Stage 301 step 4: the eps x incompleteness cross grid ABOVE K = 3.

Stage 296's cross grid stops at K in {2, 3} because its exhaustive sweep caps
at k = 4 around |R| ~ 25-29.  Here the same grid is run with the transported
localisation driving the Theta ladder, on exact-model cases with K in {5, 7,
10}, under the stage-288/299 removal regimes.  Each draw reports the exact
Theta vector it reached and the rung at which it stopped, so the ceiling the
localisation achieves is measured rather than projected.
"""
import json
import sys
import time
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import stage288_common as common            # noqa: E402
import stage301_common as s301              # noqa: E402
import stage301_prune as reach              # noqa: E402

OUT = s301.OUT
WORKERS = 40
ADDRESS = 10 * 1024 ** 3
LEVELS = (F(1, 4), F(1, 2), F(3, 4))


def worker(spec):
    t0 = time.perf_counter()
    rungs, meta = reach.ladder(spec['edges'], spec['obs'], spec['kmax'], spec['wall'])
    return dict(status='complete', rungs=rungs, meta=meta,
                n_retained=len(spec['obs']), seconds=time.perf_counter() - t0)


def draws_for(case, replicates, rng_seed):
    """Removal draws in the stage-288/299 regimes: random, clade, nearest."""
    import random
    edges, obs = case['edges'], case['observed']
    leaves = sorted(obs)
    out = [dict(regime='none', level='full', replicate=0, removed=[])]
    rng = random.Random(rng_seed)
    for level in LEVELS:
        count = int(len(leaves) * level)
        if count < 1 or len(leaves) - count < 3:
            continue
        for rep in range(replicates):
            seed = rng.randrange(2 ** 32)
            removed = common.random_removal(leaves, count, seed)
            out.append(dict(regime='random', level=str(level), replicate=rep,
                            removed=sorted(removed), seed=seed))
        try:
            cl = common.clade_removals(edges, leaves, count, replicates)
        except Exception:
            cl = []
        for rep, item in enumerate(cl):
            out.append(dict(regime='clade', level=str(level), replicate=rep,
                            removed=sorted(item['removed']),
                            clade_vertex=item['clade_vertex'],
                            complement=item['complement']))
        try:
            nr = common.nearest_source_removal(edges, case['generating_sources'],
                                               leaves, count)
            out.append(dict(regime='nearest', level=str(level), replicate=0,
                            removed=sorted(nr)))
        except Exception:
            pass
    return out


def main():
    assert (OUT / 'stage301_protocol.json').exists()
    per_cell = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    replicates = int(sys.argv[2]) if len(sys.argv) > 2 else 2
    kmax = int(sys.argv[3]) if len(sys.argv) > 3 else 10
    wall = float(sys.argv[4]) if len(sys.argv) > 4 else 900.0
    nmax = int(sys.argv[5]) if len(sys.argv) > 5 else 60
    tag = sys.argv[6] if len(sys.argv) > 6 else 'v1'
    index = common.bank_index()
    pool = [e for e in index if e['family'] in ('exact_10k', 'exact_biotrees')
            and e['K_generating'] >= 5 and e['reference_set_size'] <= nmax]
    cells = {}
    for e in pool:
        cells.setdefault(e['K_generating'], []).append(e)
    jobs = []
    for K, group in sorted(cells.items()):
        for entry in sorted(group, key=lambda e: e['name'])[:per_cell]:
            case = common.load_case(entry)
            for draw in draws_for(case, replicates, hash(entry['name']) % (2 ** 31)):
                obs = {r: d for r, d in case['observed'].items()
                       if r not in set(draw['removed'])}
                if len(obs) < 4:
                    continue
                jobs.append(dict(
                    meta=dict(arm='grid', name=entry['name'], family=entry['family'],
                              K=K, n_edges=len(case['edges']),
                              n_full=len(case['observed']), n_retained=len(obs),
                              regime=draw['regime'], level=draw['level'],
                              replicate=draw['replicate'],
                              n_removed=len(draw['removed'])),
                    spec=dict(edges=case['edges'], obs=obs,
                              kmax=min(kmax, K), wall=wall)))
    print('jobs', len(jobs), flush=True)
    common.run_pool(jobs, worker, wall + 600, ADDRESS, WORKERS,
                    str(OUT / ('stage301_grid_runs_%s.jsonl' % tag)))
    print('done', flush=True)


if __name__ == '__main__':
    main()
