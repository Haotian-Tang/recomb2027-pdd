#!/usr/bin/env python3
"""Stage295 static enumeration census of the frozen Stage43 outer search.

Replicates ONLY the preamble and the control flow of
`stage43_partial_fpt.solve_mixed` that decides HOW MANY zero faces the search
must visit, without performing any linear-programming work.  The frozen
modules are imported, never modified.

For a case with source cap K the frozen search visits, in order:

  (1) deficit supports S yielded by the frozen generator `deficit_sets`
      (counter `growth_states`);
  (2) for each S surviving the `negative - A` test, all 2^|A| sign masks over
      A = S u N(S) (counter `mass_faces`);
  (3) for each mask surviving `mixed_cover`, an active-basis enumeration of
      size at most C(m, d) with d = |S| + 1 (counter `vertex_bases`).

Phases (1) and (2) are pure combinatorics on the backbone and are counted here
exactly.  Phase (3) needs the LP and is measured separately.
"""
import json, os, sys, time
from fractions import Fraction as F
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from stage26_missing_leaf import observed_backbone
from stage25_deficit_extension import tree
from stage23_two_source_envelope import two_pass_envelope, inverse
from stage27_general_fpt import Budget, ResourceLimit
from stage41_partial_fpt import deficit_sets, mixed_cover

ROOT = os.path.dirname(HERE)


from stage295_common import load_case as _load


def load_case(relpath):
    edges, obs, entry, raw = _load(relpath)
    return raw, edges, obs


def preamble(edges, observed, cap):
    """Everything solve_partial computes before the deficit loop."""
    bb = observed_backbone(edges, observed)
    marked = {x['vertex'] for x in bb['saturation_violations']}
    adj, lookup = tree(bb['edges'])
    h, _ = two_pass_envelope(adj, observed)
    leaves = {v for v in adj if len(adj[v]) == 1}
    if any(h[r] != observed[r] for r in leaves):
        return None, 'inconsistent_leaf_envelope'
    b, Y = inverse(adj, h)
    if Y <= 0:
        return None, 'nonpositive_baseline_pendant'
    positive = {v for v in adj if b[v] > 0}
    negative = {v for v in adj if b[v] < 0}
    if len(positive) > (2 * cap if cap <= 2 else 3 * cap - 2) \
            or len(negative) > max(0, 4 * cap - 8):
        return None, 'baseline_mass_count_obstruction'
    eligible = {v for v in adj if 3 <= len(adj[v]) <= cap}
    return dict(bb=bb, marked=marked, adj=adj, lookup=lookup, b=b, Y=Y,
                positive=positive, negative=negative, eligible=eligible,
                backbone_edges=bb['edges']), None


def census(relpath, cap, max_states=5_000_000, seconds=900, face_cap=10**15):
    """Exact counts for phases (1) and (2); no LP is solved."""
    d, edges, obs = load_case(relpath)
    t0 = time.perf_counter()
    pre, stop = preamble(edges, obs, cap)
    out = {'case': relpath, 'cap': cap, 'reference_set_size': len(obs),
           'tree_edges': len(edges)}
    if pre is None:
        out.update(status='preamble_stop', reason=stop,
                   preamble_seconds=time.perf_counter() - t0)
        return out
    adj, negative, b = pre['adj'], pre['negative'], pre['b']
    out.update(backbone_edges=len(pre['backbone_edges']),
               positive=len(pre['positive']), negative=len(negative),
               eligible=len(pre['eligible']),
               anchors=len(pre['positive'] & pre['eligible']),
               marked=len(pre['marked']),
               max_degree=max(len(adj[v]) for v in adj))
    budget = Budget(seconds, max_states, None, None)
    supports = 0
    survived = 0
    faces = 0
    faces_reduced = 0
    free_hist = {}
    size_hist = {}
    A_hist = {}
    biggest = None
    truncated = None
    try:
        for S in deficit_sets(adj, pre['positive'] & pre['eligible'],
                              pre['eligible'], max(0, cap - 2), budget):
            supports += 1
            A = set(S)
            for v in S:
                A.update(adj[v])
            size_hist[len(S)] = size_hist.get(len(S), 0) + 1
            if negative - A:
                continue
            survived += 1
            A_hist[len(A)] = A_hist.get(len(A), 0) + 1
            faces += 1 << len(A)
            # Stage292 reduction: m_u > 0 is forced for u in N(S) \\ S with
            # b_u >= 0, so only S and the strictly negative baselines are free.
            Sset = set(S)
            free = sum(1 for v in A if v in Sset or b[v] < 0)
            faces_reduced += 1 << free
            free_hist[str(free)] = free_hist.get(str(free), 0) + 1
            if biggest is None or len(A) > biggest[1]:
                biggest = (len(S), len(A))
            if faces > face_cap:
                truncated = 'face_cap'
                break
    except ResourceLimit as exc:
        truncated = str(exc)
    out.update(status='truncated' if truncated else 'complete',
               truncation=truncated,
               deficit_supports=supports, supports_with_negative_covered=survived,
               mass_faces_total=faces,
               mass_faces_reduced=faces_reduced,
               free_size_histogram=dict(sorted(free_hist.items(), key=lambda kv: int(kv[0]))),
               reduction_factor=(faces / faces_reduced) if faces_reduced else None,
               support_size_histogram={str(k): v for k, v in sorted(size_hist.items())},
               A_size_histogram={str(k): v for k, v in sorted(A_hist.items())},
               max_support_size=max(size_hist) if size_hist else 0,
               max_A_size=max(A_hist) if A_hist else 0,
               theory_max_A=max(0, 4 * cap - 8),
               theory_face_bound=1 << max(0, 4 * cap - 8),
               census_seconds=time.perf_counter() - t0)
    return out


if __name__ == '__main__':
    plan = json.load(open(sys.argv[1]))
    rows = []
    for item in plan:
        row = census(item['case'], item['cap'],
                     seconds=item.get('seconds', 900),
                     face_cap=item.get('face_cap', 10**15))
        rows.append(row)
        print(json.dumps(row), flush=True)
    json.dump(rows, open(sys.argv[2], 'w'), indent=1)
