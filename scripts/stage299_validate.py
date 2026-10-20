#!/usr/bin/env python3
"""Stage 299 step 1: validate the only new code that decides anything.

Two things in stage 299 are not frozen: the attaining-point geometry
(`point_distance`, `edge_hops`) and the use of stage 298's k = 1 closed form on
REDUCED leaf sets, where one side of an edge can hold no retained reference at
all.  Both are checked here against independent computations before any
measurement run.

GEOMETRY.  point_distance(e1, x1, e2, x2) is checked by SUBDIVISION: the two
points are inserted as real vertices (splitting their edges), and the distance
is read off the frozen BFS `stage288_common.bfs_from` on the subdivided tree.
edge_hops is checked the same way on a unit-length copy.

ONE-SIDED PANELS.  theta1_all_edges is compared edge by edge against the frozen
`stage296_common.theta_support` on leaf sets that have been thinned hard enough
that many edges carry retained references on one side only - including a
hand-built star where the entire retained panel sits on one side of the
attaining edge.
"""
import json
import random
import sys
import time
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import stage288_common as common               # noqa: E402
import stage296_common as s296                 # noqa: E402
import stage298_common as s298                 # noqa: E402
import stage299_common as s299                 # noqa: E402

OUT = s299.OUT


def subdivide(edges, e1, x1, e2, x2):
    """Insert the two points as vertices; return (edges', q1, q2)."""
    nxt = 1 + max(max(u, v) for u, v, _ in edges)
    out = []
    q1 = q2 = None            # separate handles: a dict keyed by edge index
    for i, (u, v, L) in enumerate(edges):   # collides when e1 == e2
        if i == e1 and i == e2:
            lo, hi = (F(x1), F(x2)) if F(x1) <= F(x2) else (F(x2), F(x1))
            w1, w2 = nxt, nxt + 1
            nxt += 2
            out.append((u, w1, lo * L))
            out.append((w1, w2, (hi - lo) * L))
            out.append((w2, v, (1 - hi) * L))
            q1, q2 = (w1, w2) if F(x1) <= F(x2) else (w2, w1)
        elif i == e1 or i == e2:
            x = x1 if i == e1 else x2
            w = nxt
            nxt += 1
            if i == e1:
                q1 = w
            else:
                q2 = w
            out.append((u, w, F(x) * L))
            out.append((w, v, (1 - F(x)) * L))
        else:
            out.append((u, v, L))
    return out, q1, q2


def check_geometry(edges, pairs, rng):
    rows = []
    for e1, e2 in pairs:
        x1 = F(rng.randint(1, 9), 10)
        x2 = F(rng.randint(1, 9), 10)
        mine = s299.point_distance(edges, e1, x1, e2, x2)
        sub, q1, q2 = subdivide(edges, e1, x1, e2, x2)
        ref = common.bfs_from(sub, [q1])[q1][q2]
        unit = [(u, v, F(1)) for u, v, _L in sub]
        hops_ref = common.bfs_from(unit, [q1])[q1][q2]
        # a subdivided path of h hops between the two inserted vertices crosses
        # h - 2 whole edges plus the two halves, so the edge-hop count is h - 2
        hops_mine = s299.edge_hops(edges, e1, e2)
        # the subdivided path between the inserted vertices crosses the two
        # halves plus the whole edges between them: h - 2 hops when the two
        # points sit on different edges, h - 1 when they share one
        hops_ref_edges = int(hops_ref) - (1 if e1 == e2 else 2)
        rows.append(dict(e1=e1, e2=e2, x1=str(x1), x2=str(x2),
                         distance_agrees=(mine == ref),
                         mine=str(mine), reference=str(ref),
                         hops_mine=hops_mine, hops_reference=hops_ref_edges,
                         hops_agree=(hops_mine == hops_ref_edges)))
    return rows


def check_closed_form(case_name, edges, leaves, dvals, label):
    """Every edge decided twice on one (possibly thinned) leaf set."""
    t0 = time.perf_counter()
    cf = s298.theta1_all_edges(edges, leaves, dvals)
    fz = s299.theta1_frozen_all_edges(edges, leaves, dvals)
    assert fz['status'] == 'ok', fz
    mism = [i for i in range(len(edges)) if cf['theta'][i] != fz['theta'][i]]
    one_sided = 0
    aggs = s298.edge_aggregates(edges, leaves, dvals)
    for A_in, B_in, A_out, B_out, _c, _L in aggs:
        if A_in is None or A_out is None:
            one_sided += 1
    return dict(case=case_name, label=label, n_edges=len(edges), n_leaves=len(leaves),
                one_sided_edges=one_sided, mismatches=len(mism),
                mismatch_examples=mism[:5],
                theta1_closed=str(cf['theta1']), theta1_frozen=str(fz['theta1']),
                theta1_agrees=(cf['theta1'] == fz['theta1']),
                seconds=time.perf_counter() - t0)


def star_case():
    """A hand-built star where a hard thinning leaves one side empty."""
    edges = [(0, 1, F(1)), (0, 2, F(2)), (0, 3, F(3)), (3, 4, F(1)), (3, 5, F(5))]
    obs = {1: F('3/2'), 2: F('5/2'), 4: F('9/2'), 5: F('17/2')}
    return edges, obs


def main():
    # v1 is retained: its three geometry disagreements were a defect in the
    # reference computation above (one inserted vertex for a same-edge pair),
    # not in point_distance/edge_hops, which v2 confirms.
    out_path = OUT / 'stage299_validation_v3.json'
    assert not out_path.exists(), 'Never overwrite frozen evidence'
    rng = random.Random(2990001)
    protocol = json.loads((OUT / 'stage299_protocol.json').read_text())
    geometry, closed = [], []

    # --- geometry on three real trees of very different size
    picks = []
    for arm in ('krepp_ani', 'krepp_placement', 'travelers'):
        row = sorted([c for c in protocol['cases'] if c['arm'] == arm],
                     key=lambda c: c['n_edges'])[-1]
        picks.append(row)
    for row in picks:
        case = s299.prepare_case(s299.load_json(row['path'], row['sha256']), row['name'])
        edges = case['edges']
        m = len(edges)
        pairs = [(rng.randrange(m), rng.randrange(m)) for _ in range(40)]
        pairs += [(0, 0), (0, 1)]
        t0 = time.perf_counter()
        rows = check_geometry(edges, pairs, rng)
        for r in rows:
            r['case'] = row['name']
            r['n_edges'] = m
        geometry.extend(rows)
        print('%s: %d geometry pairs, %.1f s' % (row['name'], len(rows), time.perf_counter() - t0),
              flush=True)

    # --- closed form vs frozen LP on thinned leaf sets of a small real case
    row = sorted([c for c in protocol['cases'] if c['arm'] == 'krepp_ani'],
                 key=lambda c: c['n_edges'])[len([c for c in protocol['cases']
                                                  if c['arm'] == 'krepp_ani']) // 2]
    case = s299.prepare_case(s299.load_json(row['path'], row['sha256']), row['name'])
    leaves_all = sorted(case['observed'])
    for keep in (len(leaves_all), max(3, len(leaves_all) // 4), 4, 3):
        sub = sorted(rng.sample(leaves_all, keep)) if keep < len(leaves_all) else leaves_all
        edges, leaves, dvals = s299.case_arrays(case, sub)
        closed.append(check_closed_form(row['name'], edges, leaves, dvals,
                                        'thinned_to_%d' % keep))
        print('  thinned to %d: mismatches %d, one-sided edges %d'
              % (keep, closed[-1]['mismatches'], closed[-1]['one_sided_edges']), flush=True)

    # --- the synthetic star, retained panel entirely on one side
    edges, obs = star_case()
    for sub in ([1, 2, 4, 5], [4, 5], [1, 2], [1, 4]):
        leaves = sorted(sub)
        dvals = [obs[r] for r in leaves]
        closed.append(check_closed_form('synthetic_star', edges, leaves, dvals,
                                        'leaves_' + '_'.join(str(r) for r in leaves)))
    out = dict(stage=299, step='validate',
               geometry=dict(n=len(geometry),
                             distance_agree=sum(r['distance_agrees'] for r in geometry),
                             hops_agree=sum(r['hops_agree'] for r in geometry),
                             rows=geometry),
               closed_form=dict(n=len(closed),
                                total_mismatches=sum(r['mismatches'] for r in closed),
                                theta1_agree=sum(r['theta1_agrees'] for r in closed),
                                rows=closed))
    out_path.write_text(json.dumps(common.serial(out), indent=1))
    print('geometry %d/%d distances, %d/%d hops; closed form %d/%d theta1, %d mismatches'
          % (out['geometry']['distance_agree'], len(geometry),
             out['geometry']['hops_agree'], len(geometry),
             out['closed_form']['theta1_agree'], len(closed),
             out['closed_form']['total_mismatches']))
    assert out['geometry']['distance_agree'] == len(geometry)
    assert out['geometry']['hops_agree'] == len(geometry)
    assert out['closed_form']['total_mismatches'] == 0
    assert out['closed_form']['theta1_agree'] == len(closed)
    print('wrote', out_path)


if __name__ == '__main__':
    main()

