"""Independent exact-arithmetic check of the two-source support identifiability corollary.

Claim (Supplementary Appendix E, Corollary): on a reduced positive-rational tree with every leaf
observed and exact distances, the occupied original edges of an at-most-two-source explanation are
determined by the observations, with no assumption that the source weights are equal or known.

The check draws random reduced trees, places 1 or 2 strict interior sources with UNEQUAL random
rational weights, forms the exact leaf distances, and verifies the four steps of the proof:
  1. the lower Lipschitz envelope h equals the true vertex potential f at every vertex (S = {});
  2. the tree-distance inverse recovers the mean pendant y exactly;
  3. the true support covers the positive-mass set P;
  4. the <=2-edge cover whose endpoint union equals P is unique and equals the true support.

A control subdivides one edge to create a degree-two internal vertex and confirms that the reduced
hypothesis is load-bearing: f ceases to equal h on a substantial fraction of draws.

All arithmetic is exact (fractions.Fraction). Run: python verify_two_source_corollary.py
"""
from fractions import Fraction as F
import collections, itertools, random


def rand_reduced_tree(nleaf, rng):
    """Random tree whose internal vertices all have degree >= 3, with positive rational lengths."""
    while True:
        n_int = rng.randint(1, max(1, nleaf - 2))
        verts, leaves = list(range(n_int)), list(range(n_int, n_int + nleaf))
        edges = [(v, rng.randrange(v)) for v in range(1, n_int)]
        edges += [(L, rng.randrange(n_int)) for L in leaves]
        deg = collections.Counter()
        for a, b in edges:
            deg[a] += 1
            deg[b] += 1
        if all(deg[v] >= 3 for v in verts) and len(edges) == n_int - 1 + nleaf:
            return verts, leaves, [(a, b, F(rng.randint(1, 9), rng.randint(1, 5))) for a, b in edges]


def all_pairs(nv, edges):
    adj = collections.defaultdict(list)
    for a, b, l in edges:
        adj[a].append((b, l))
        adj[b].append((a, l))
    D = [[None] * nv for _ in range(nv)]
    for s in range(nv):
        D[s][s] = F(0)
        stack = [s]
        while stack:
            u = stack.pop()
            for w, l in adj[u]:
                if D[s][w] is None:
                    D[s][w] = D[s][u] + l
                    stack.append(w)
    return D, adj


def main(n_trials=400, n_control=300, seed=20261020):
    rng = random.Random(seed)
    verified = 0
    failures = []
    for _ in range(n_trials):
        verts, leaves, edges = rand_reduced_tree(rng.randint(4, 9), rng)
        nv = len(verts) + len(leaves)
        D, adj = all_pairs(nv, edges)
        k = rng.choice([1, 2])
        eidx = rng.sample(range(len(edges)), k)
        ps = [F(rng.randint(1, 20)) for _ in range(k)]
        tot = sum(ps)
        ps = [p / tot for p in ps]
        if k == 2 and ps[0] == ps[1]:
            ps = [F(1, 3), F(2, 3)]                      # force unequal weights
        xs = [F(rng.randint(1, 6), 7) for _ in range(k)]
        y = F(rng.randint(1, 5), rng.randint(1, 4))

        def potential(v):
            return y + sum(p * min(D[v][edges[e][0]] + x * edges[e][2],
                                   D[v][edges[e][1]] + (1 - x) * edges[e][2])
                           for e, p, x in zip(eidx, ps, xs))

        dobs = {r: potential(r) for r in leaves}
        h = [max(dobs[r] - D[v][r] for r in leaves) for v in range(nv)]
        f = [potential(v) for v in range(nv)]
        if h != f:
            failures.append("f != h")
            continue
        deg = collections.Counter()
        for a, b, l in edges:
            deg[a] += 1
            deg[b] += 1
        W = sum(l for _, _, l in edges)
        tau = [2 - deg[v] for v in range(nv)]
        Lh = [sum((h[v] - h[u]) / l for u, l in adj[v]) for v in range(nv)]
        m = [(F(tau[v]) - Lh[v]) / 2 for v in range(nv)]
        if (sum(F(tau[v]) * h[v] for v in range(nv)) - W) / 2 != y:
            failures.append("y not recovered")
            continue
        P = {v for v in range(nv) if m[v] > 0}
        true_support = {frozenset((edges[e][0], edges[e][1])) for e in eidx}
        exact_covers = []
        for size in (1, 2):
            for comb in itertools.combinations(range(len(edges)), size):
                union = set()
                for ci in comb:
                    union |= {edges[ci][0], edges[ci][1]}
                if union == P:
                    exact_covers.append({frozenset((edges[ci][0], edges[ci][1])) for ci in comb})
        if len(exact_covers) != 1 or exact_covers[0] != true_support:
            failures.append("cover not unique")
            continue
        verified += 1
    print(f"reduced trees, k in {{1,2}}, unequal rational weights: {verified}/{n_trials} verified")
    print(f"failures: {collections.Counter(failures) or 'none'}")

    rng = random.Random(777)
    broke = 0
    for _ in range(n_control):
        verts, leaves, edges = rand_reduced_tree(rng.randint(4, 7), rng)
        nv = len(verts) + len(leaves)
        a, b, l = edges[rng.randrange(len(edges))]
        mid = nv
        nv += 1
        edges = [e for e in edges if not (e[0] == a and e[1] == b)]
        edges += [(a, mid, l / 2), (mid, b, l / 2)]        # degree-two internal vertex
        D, _ = all_pairs(nv, edges)
        eidx = rng.sample(range(len(edges)), 2)
        ps, xs = [F(1, 3), F(2, 3)], [F(rng.randint(1, 6), 7) for _ in range(2)]
        y = F(rng.randint(1, 5), 2)

        def potential(v):
            return y + sum(p * min(D[v][edges[e][0]] + x * edges[e][2],
                                   D[v][edges[e][1]] + (1 - x) * edges[e][2])
                           for e, p, x in zip(eidx, ps, xs))

        dobs = {r: potential(r) for r in leaves}
        if [max(dobs[r] - D[v][r] for r in leaves) for v in range(nv)] != [potential(v) for v in range(nv)]:
            broke += 1
    print(f"control with a degree-two internal vertex: f != h in {broke}/{n_control} draws")
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
