#!/usr/bin/env python3
"""Stage292 removal orders.

Identical in definition to Stage290's `removal_order`, but the
adversarial (nearest-to-source) order is computed from a handful of BFS
distance columns instead of the dense |R| x |V| float matrix, so it also
works on 10^4-leaf panels.  Equality with the Stage290 order is checked in
the validation harness on every case where both are affordable.
"""
import sys
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
sys.dont_write_bytecode = True

from stage290_erasure_core import removal_order                         # noqa: E402
from stage292_all_minimisers import VertexDistances                     # noqa: E402


def nearest_source_distance(edges, sources, Rlist, dists):
    best = {}
    for s in sources:
        u, v, L = edges[s['edge']]
        x = F(s['x'])
        L = F(L)
        cu, cv = dists.column(u), dists.column(v)
        for r in Rlist:
            d = min(cu[r] + x * L, cv[r] + (1 - x) * L)
            if r not in best or d < best[r]:
                best[r] = d
    return best


def order_for(edges, st, regime, seed, sources=None, dists=None):
    if regime != 'adversarial_nearest':
        return removal_order(None, st, regime, seed, tables=None)
    Rlist = sorted(st['R'])
    if dists is None:
        dists = VertexDistances(edges)
    best = nearest_source_distance(edges, sources, Rlist, dists)
    return sorted(Rlist, key=lambda r: (best[r], r))

