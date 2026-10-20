#!/usr/bin/env python3
"""Stage 299 shared library: Theta_1 under reference removal, on REAL distances.

THE LOOP THIS CLOSES.  Every reference-incompleteness measurement in this
project (stages 288-295) runs on instances with synthetic planted sources, and
needs either ground truth or an admissible pendant floor.  Stage 298 built an
instrument that needs neither - Theta_1, the exact smallest tolerance at which
a SINGLE source explains the observations - but ran it only at the full
reference panel and said so.  This stage runs the stage-288 removal regimes
through that instrument on four collections: the 296 exact partners (paired
baseline), the 296 krepp ANI cases, the 60 krepp placement cases and the 29
travellers'-diarrhoea stool samples.

NOTHING HERE DECIDES A NUMBER.  Every theta is decided by frozen code:
  stage296_common.theta_support         exact rational two-phase simplex,
  stage296_common.independent_deviation frozen min(d(r,u)+xL, d(r,v)+(1-x)L),
  stage298_common.theta1_all_edges      the k=1 closed form, which stage 298
                                        validated against the frozen LP on
                                        653,807 edges with 0 disagreements,
  stage288_common.bfs_from / predict / random_removal / clade_removals /
  nearest_source_removal, stage294_common.taxon_removal.
This module only builds leaf sets, routes them, and measures geometry
(attaining-edge migration) that no frozen routine reports.

DRAW KEYING.  stage288_incompleteness.removal_draws seeds its random draws on
`case['name']`.  The noisy and exact partners have different names, so stage
294 fixed the convention: the seed key uses the PAIRED EXACT case name, which
makes the exact, ANI and placement arms share leaf sets leaf for leaf.  Stage
299 follows stage 294 and VERIFIES the import: on every case whose krepp panel
equals its exact panel, the random draws built here must equal
stage288_incompleteness.removal_draws(exact_case, level) element by element.

ONE BOOKKEEPING FIX.  clade_removals returns the whole-clade set of TREE
leaves, which for a panel that is a strict subset of the tree's leaves can
contain unobserved ones.  Removal is defined on the observed panel, so every
draw here records both the raw clade leaf set and the EFFECTIVE removal
(raw intersected with the arm's panel); all counts and levels use the
effective one.  The retained set - the only thing the LP sees - is identical
either way.
"""
import os
for _k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_k] = '1'
import hashlib
import json
import sys
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / 'scripts') not in sys.path:
    sys.path.insert(0, str(ROOT / 'scripts'))
import stage288_common as common               # noqa: E402  frozen
import stage288_incompleteness as inc          # noqa: E402  frozen draw plan
import stage294_common as s294                 # noqa: E402  frozen taxon regime
import stage296_common as s296                 # noqa: E402  frozen exact LP
import stage298_common as s298                 # noqa: E402  frozen closed form

OUT = ROOT / 'results' / 'stage299'
S298 = ROOT / 'results' / 'stage298'
ZERO, ONE = F(0), F(1)

LEVELS = tuple(inc.LEVELS)                     # ('one','0.10','0.25','0.50','0.75')
SMALL_ARMS = ('exact_partner', 'krepp_ani', 'krepp_placement')
TD_ARM = 'travelers'
TD_SAMPLE_EDGES = 20                           # frozen-LP spot checks per TD leaf set


# --------------------------------------------------------------- bank access

def bank_by_name():
    return {e['name']: e for e in common.bank_index()}


def load_json(relpath, sha256=None):
    path = ROOT / relpath
    if sha256 is not None:
        assert common.sha(path) == sha256, 'input hash mismatch: %s' % relpath
    return json.loads(path.read_text())


def prepare_case(raw, name):
    """Normalise one raw bank case to the frozen stage-288/294 representation.

    The frozen removal regimes index edges with Fraction lengths, key observed
    and leaf_labels by INT leaf id and carry generating_sources with Fraction p
    and x.  Reproducing stage 288 random draws requires the same int keys it
    used (sorted int order differs from sorted string order), which the
    registered import check verifies.
    """
    return dict(
        name=name,
        edges=[(int(u), int(v), F(L)) for u, v, L in raw['edges']],
        observed={int(r): F(d) for r, d in raw['observed'].items()},
        cap=raw.get('cap'),
        generating_sources=[dict(edge=int(s['edge']), p=F(s['p']), x=F(s['x']))
                            for s in (raw.get('generating_sources') or [])],
        leaf_labels={int(k): v for k, v in (raw.get('leaf_labels') or {}).items()})


def case_arrays(case, leaves=None):
    """(edges, leaves, dvals) with `leaves` optionally restricted."""
    edges = [(int(u), int(v), F(L)) for u, v, L in case['edges']]
    obs = {int(r): F(d) for r, d in case['observed'].items()}
    keys = sorted(obs) if leaves is None else sorted(int(r) for r in leaves)
    return edges, keys, [obs[r] for r in keys]


def median_observed(dvals):
    s = sorted(dvals)
    return s[len(s) // 2]


# ------------------------------------------------------------- removal draws

def _random_draw(exact_name, level, rep, leaves, count):
    key = '%s|%s|%s|random|%d' % (inc.SEED, exact_name, level, rep)
    seed = int(hashlib.sha256(key.encode()).hexdigest()[:8], 16)
    return common.random_removal(leaves, count, seed), seed


def small_draws(case, exact_name, ranks=None):
    """Full panel plus every removal draw for one small case, all levels.

    `case` supplies the tree and the panel; `sources` the generating sources
    for the nearest regime; `ranks`/`labels` the stage-294 taxonomic regime.
    """
    edges = case['edges']
    sources = case['generating_sources']
    labels = case['leaf_labels'] or None
    leaves = sorted(case['observed'])
    panel = set(leaves)
    n = len(leaves)
    out = [dict(level='full', regime='none', replicate=0, removed=[], raw_removed=[],
                n_removed=0, n_retained=n, detail={})]
    for level in LEVELS:
        count = inc.removal_count(n, level)
        if count is None:
            continue
        draws = []
        for rep in range(inc.REPLICATES['random']):
            removed, seed = _random_draw(exact_name, level, rep, leaves, count)
            draws.append(dict(regime='random', replicate=rep, removed=removed,
                              detail={'seed': seed, 'seed_key_case': exact_name}))
        for rep, clade in enumerate(common.clade_removals(edges, leaves, count,
                                                          inc.REPLICATES['clade'])):
            draws.append(dict(regime='clade', replicate=rep, removed=clade['removed'],
                              detail={'clade_vertex': clade['clade_vertex'],
                                      'complement': clade['complement']}))
        draws.append(dict(regime='nearest', replicate=0,
                          removed=common.nearest_source_removal(edges, sources, leaves, count),
                          detail={'sources': len(sources)}))
        if ranks is not None and labels is not None:
            for rep, order in enumerate(s294.TAXON_ORDERS):
                drawn = s294.taxon_removal(ranks, labels, leaves, count, order)
                if drawn is None:
                    continue
                draws.append(dict(regime='taxon', replicate=rep, removed=drawn['removed'],
                                  detail={'rank': drawn['rank'], 'groups': drawn['groups'],
                                          'order': order, 'target_count': count,
                                          'n_groups_at_rank': drawn['n_groups_at_rank']}))
        for d in draws:
            raw = sorted(int(r) for r in d['removed'])
            eff = sorted(r for r in raw if r in panel)
            d['raw_removed'] = raw
            d['removed'] = eff
            d['n_removed'] = len(eff)
            d['n_retained'] = n - len(eff)
            d['level'] = level
        out.extend([d for d in draws if d['n_retained'] >= 3 and d['n_removed'] >= 1])
    return out


def td_draws(case, full_source):
    """Travellers' samples: no ground truth, so `nearest` is data-derived.

    The nearest-to-source regime needs source points.  There are none for a
    real stool sample, so the regime is run against the stage-298 full-panel
    Theta_1 attaining point - the data's own best single-source explanation -
    through the SAME frozen nearest_source_removal routine.  It is reported
    under the distinct name `nearest_theta1` and never pooled with the
    truth-based `nearest` of the other three arms.
    """
    edges = case['edges']
    leaves = sorted(case['observed'])
    panel = set(leaves)
    n = len(leaves)
    out = [dict(level='full', regime='none', replicate=0, removed=[], raw_removed=[],
                n_removed=0, n_retained=n, detail={})]
    for level in LEVELS:
        count = inc.removal_count(n, level)
        if count is None:
            continue
        draws = []
        for rep in range(inc.REPLICATES['random']):
            removed, seed = _random_draw(case['name'], level, rep, leaves, count)
            draws.append(dict(regime='random', replicate=rep, removed=removed,
                              detail={'seed': seed, 'seed_key_case': case['name']}))
        for rep, clade in enumerate(common.clade_removals(edges, leaves, count,
                                                          inc.REPLICATES['clade'])):
            draws.append(dict(regime='clade', replicate=rep, removed=clade['removed'],
                              detail={'clade_vertex': clade['clade_vertex'],
                                      'complement': clade['complement']}))
        draws.append(dict(regime='nearest_theta1', replicate=0,
                          removed=common.nearest_source_removal(edges, [full_source],
                                                                leaves, count),
                          detail={'source': dict(edge=int(full_source['edge']),
                                                 x=str(full_source['x'])),
                                  'origin': 'stage298_full_panel'}))
        for d in draws:
            raw = sorted(int(r) for r in d['removed'])
            eff = sorted(r for r in raw if r in panel)
            d['raw_removed'] = raw
            d['removed'] = eff
            d['n_removed'] = len(eff)
            d['n_retained'] = n - len(eff)
            d['level'] = level
        out.extend([d for d in draws if d['n_retained'] >= 3 and d['n_removed'] >= 1])
    return out


def draw_digest(plan_rows):
    """sha256 over the full registered draw plan, so the run cannot drift."""
    h = hashlib.sha256()
    for row in plan_rows:
        h.update(('%s|%s|%s|%s|%d|%s\n' % (
            row['arm'], row['name'], row['level'], row['regime'], row['replicate'],
            ','.join(str(r) for r in row['removed']))).encode())
    return h.hexdigest()


# ------------------------------------------------------------- Theta_1 routes

def theta1_frozen_all_edges(edges, leaves, dvals):
    """Theta_1 by the frozen per-support LP on EVERY edge, no cutoff."""
    Dv = s296.leaf_vectors(edges, leaves)
    thetas, best, bz = [], None, None
    for i in range(len(edges)):
        r = s296.theta_support(s296.support_cols(Dv, edges, (i,)), dvals, 1)
        if r['status'] != 'optimal':
            return dict(status='lp_' + r['status'], edge=i)
        thetas.append(r['value'])
        if best is None or r['value'] < thetas[best]:
            best, bz = i, r['z']
    return dict(status='ok', theta=thetas, best=best, z=bz, theta1=thetas[best])


def attaining(thetas, theta1):
    return [i for i, t in enumerate(thetas) if t is not None and t == theta1]


def degrees(edges):
    deg = {}
    for u, v, _L in edges:
        deg[u] = deg.get(u, 0) + 1
        deg[v] = deg.get(v, 0) + 1
    return deg


def is_terminal(deg, edge):
    u, v, _L = edge
    return bool(deg.get(u) == 1 or deg.get(v) == 1)


# --------------------------------------------------- geometry of the attainer

def point_distance(edges, e1, x1, e2, x2, dist=None):
    """Exact tree distance between two interior points (e, x), x from edges[e][0].

    On a tree the path from a point inside e1 to a point inside e2 leaves e1
    through one of its endpoints and enters e2 through one of its endpoints,
    so the distance is the minimum over the four endpoint pairs of
    offset1 + d(a, b) + offset2.  Same edge is handled directly.
    `dist` may supply a prebuilt common.bfs_from table keyed by e1's endpoints.
    """
    u1, v1, L1 = edges[e1]
    u2, v2, L2 = edges[e2]
    x1, x2 = F(x1), F(x2)
    if e1 == e2:
        return abs(x1 - x2) * L1
    if dist is None:
        dist = common.bfs_from(edges, sorted({u1, v1}))
    off1 = {u1: x1 * L1, v1: (ONE - x1) * L1}
    off2 = {u2: x2 * L2, v2: (ONE - x2) * L2}
    best = None
    for a in (u1, v1):
        for b in (u2, v2):
            value = off1[a] + dist[a][b] + off2[b]
            if best is None or value < best:
                best = value
    return best


def edge_hops(edges, e1, e2, hop=None):
    """Topological hops between two edges: 0 if identical or adjacent."""
    if e1 == e2:
        return 0
    u1, v1, _ = edges[e1]
    u2, v2, _ = edges[e2]
    if hop is None:
        unit = [(u, v, ONE) for u, v, _L in edges]
        hop = common.bfs_from(unit, sorted({u1, v1}))
    best = None
    for a in (u1, v1):
        for b in (u2, v2):
            value = hop[a][b]
            if best is None or value < best:
                best = value
    return int(best)


def serial(value):
    return common.serial(value)

