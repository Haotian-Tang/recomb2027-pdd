#!/usr/bin/env python3
"""Exact controls for pendant-assisted noisy PDD endpoint-mass bounds.

Independent tree distances, source interpolation, envelope, and inverse mass
calculations. This is a proof-control program, not a global interval optimizer.
No previous project solver or numerical optimization package is imported.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json
import random

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/stage43_theory_robust_mass.json"


def tree(edges):
    adjacency = {}
    for u, v, length in edges:
        assert length > 0
        adjacency.setdefault(u, {})[v] = length
        adjacency.setdefault(v, {})[u] = length
    assert len(edges) + 1 == len(adjacency)
    assert all(len(row) == 1 or len(row) >= 3 for row in adjacency.values())
    leaves = [v for v, row in adjacency.items() if len(row) == 1]
    distances = {}
    for start in adjacency:
        row = {start: F(0)}
        queue = [start]
        for u in queue:
            for v, length in adjacency[u].items():
                if v not in row:
                    row[v] = row[u] + length
                    queue.append(v)
        assert len(row) == len(adjacency)
        distances[start] = row
    return edges, adjacency, leaves, distances


def cherry_tree(number, unequal=False):
    pairs = []
    if number == 2:
        pairs.append(("c0", "c1"))
    else:
        pairs += [(f"v{i}", f"v{i+1}") for i in range(number - 3)]
        pairs += [("c0", "v0")]
        pairs += [(f"c{i}", f"v{i-1}") for i in range(1, number - 1)]
        pairs += [(f"c{number-1}", f"v{number-3}")]
    pairs += [(f"c{i}", f"{arm}{i}") for i in range(number) for arm in ("a", "b")]
    edges = [(u, v, F(1 + i % 5, 1 + i % 3) if unequal else F(1))
             for i, (u, v) in enumerate(pairs)]
    return tree(edges)


def forward(geometry, sources, pendant):
    edges, adjacency, leaves, distances = geometry
    assert sum(p for _, _, p in sources) == 1
    masses = {v: F(0) for v in adjacency}
    for edge_index, fraction, weight in sources:
        assert 0 < fraction < 1 and weight > 0
        u, v, _ = edges[edge_index]
        masses[u] += weight * (1 - fraction)
        masses[v] += weight * fraction
    assert len(set(e for e, _, _ in sources)) == len(sources)
    potential = {v: pendant + sum((masses[u] * distances[v][u] for u in adjacency), F(0))
                 for v in adjacency}
    return masses, potential


def baseline(geometry, observations):
    edges, adjacency, leaves, distances = geometry
    h = {v: max(observations[r] - distances[v][r] for r in leaves) for v in adjacency}
    masses = {v: (F(2 - len(row)) - sum(((h[v] - h[u]) / length
              for u, length in row.items()), F(0))) / 2 for v, row in adjacency.items()}
    y = (sum((F(2 - len(adjacency[v])) * h[v] for v in adjacency), F(0))
         - sum((length for _, _, length in edges), F(0))) / 2
    return h, masses, y


def bounds(geometry, observations, epsilon, pendant_floor):
    _, adjacency, leaves, _ = geometry
    h, b, y = baseline(geometry, observations)
    gap = y - pendant_floor + (len(leaves) - 1) * epsilon
    intervals = {}
    for v, row in adjacency.items():
        conductance = sum((1 / length for length in row.values()), F(0))
        neighbor_coefficient = max([1 / (length * (len(adjacency[u]) - 2))
                                   for u, length in row.items() if len(adjacency[u]) >= 3]
                                  + [F(0)])
        lower = b[v] - epsilon * conductance
        if len(row) >= 3:
            lower -= gap * conductance / (len(row) - 2)
        upper = b[v] + epsilon * conductance + gap * neighbor_coefficient
        intervals[v] = (lower, upper)
    return h, b, y, gap, intervals


def main(write=True):
    rng = random.Random(430317)
    records = []
    strict_count = 0
    for number in (2, 3, 4, 6, 10):
        for unequal in (False, True):
            geometry = cherry_tree(number, unequal)
            edges, adjacency, leaves, _ = geometry
            for trial in range(20):
                count = 1 + trial % min(8, len(edges))
                chosen = rng.sample(range(len(edges)), count)
                raw = [rng.randrange(1, 12) for _ in chosen]
                total = sum(raw)
                sources = [(edge, F(rng.randrange(1, 10), 10), F(w, total))
                           for edge, w in zip(chosen, raw)]
                pendant = F(1 + trial % 4, 7)
                masses, potential = forward(geometry, sources, pendant)
                epsilon = F(1, (1000, 10**6, 10**12)[trial % 3])
                nominal = {r: potential[r] + epsilon * F(rng.randrange(-4, 5), 4)
                           for r in leaves}
                # A lower bound on the actual pendant is sufficient.
                floor = pendant if trial % 2 == 0 else pendant / 2
                h, b, y, gap, intervals = bounds(geometry, nominal, epsilon, floor)
                true_h, _, true_y = baseline(geometry, {r: potential[r] for r in leaves})
                deficit = {v: potential[v] - true_h[v] for v in adjacency}
                assert all(value >= 0 for value in deficit.values())
                assert all(deficit[r] == 0 for r in leaves)
                deficit_budget = sum((F(len(adjacency[v]) - 2) * deficit[v]
                                      for v in adjacency if len(adjacency[v]) >= 3), F(0))
                assert deficit_budget == 2 * (true_y - pendant)
                assert gap >= 0 and deficit_budget <= 2 * gap
                assert all(abs(true_h[v] - h[v]) <= epsilon for v in adjacency)
                assert abs(true_y - y) <= (len(leaves) - 1) * epsilon
                assert all(lo <= masses[v] <= hi for v, (lo, hi) in intervals.items())
                alpha = min(p for _, _, p in sources)
                rho = min(min(x, 1 - x) for _, x, _ in sources)
                excluded = {v for v, (_, upper) in intervals.items() if upper < alpha * rho}
                candidate_edges = [j for j, (u, v, _) in enumerate(edges)
                                   if u not in excluded and v not in excluded]
                assert set(chosen) <= set(candidate_edges)
                strict_count += len(adjacency)
                records.append({"cherries": number, "unequal_lengths": unequal, "trial": trial,
                                "source_count": count, "epsilon": epsilon,
                                "pendant_floor_equals_truth": floor == pendant,
                                "gap_upper": gap, "true_deficit_budget": deficit_budget,
                                "excluded_vertices": len(excluded), "candidate_edges": len(candidate_edges)})

    one_source = []
    for number in (2, 3, 4, 6, 10):
        for unequal in (False, True):
            geometry = cherry_tree(number, unequal)
            edges, adjacency, leaves, _ = geometry
            for chosen in (0, len(edges) - 1):
                masses, potential = forward(geometry, [(chosen, F(1, 2), F(1))], F(1))
                alpha, rho = F(1, 4), F(1, 4)
                max_conductance = max(sum((1 / ell for ell in row.values()), F(0))
                                      for row in adjacency.values())
                separation = alpha * rho / (len(leaves) * max_conductance)
                epsilon = separation / 2
                _, _, _, gap, intervals = bounds(geometry, {r: potential[r] for r in leaves}, epsilon, F(1))
                excluded = {v for v, (_, upper) in intervals.items() if upper < alpha * rho}
                candidate_edges = [j for j, (u, v, _) in enumerate(edges)
                                   if u not in excluded and v not in excluded]
                assert candidate_edges == [chosen]
                one_source.append({"cherries": number, "unequal_lengths": unequal,
                                   "true_edge": chosen, "alpha": alpha, "rho": rho,
                                   "proved_separation": separation, "tested_error": epsilon,
                                   "remaining_edges": candidate_edges})

    # Unknown pendant: equal exact observations with 1 and 3 well-inside sources.
    star = tree([(f"a{i}", "c", F(1)) for i in range(3)])
    _, one_pot = forward(star, [(0, F(11, 12), F(1))], F(13, 12))
    dense = [(0, F(1, 2), F(1, 3)), (1, F(3, 4), F(1, 3)), (2, F(3, 4), F(1, 3))]
    _, dense_pot = forward(star, dense, F(1))
    assert all(one_pot[r] == dense_pot[r] for r in star[2])

    # Known pendant: 2-source and 3-source generators with the SAME endpoint mass.
    quartet = tree([("a", "u", F(1)), ("b", "u", F(1)), ("u", "v", F(1)),
                    ("v", "c", F(1)), ("v", "d", F(1))])
    two_m, two_pot = forward(quartet, [(0, F(1, 2), F(1, 2)), (3, F(1, 2), F(1, 2))], F(1))
    three_sources = [(0, F(1, 3), F(3, 8)), (2, F(1, 2), F(1, 4)), (3, F(2, 3), F(3, 8))]
    three_m, three_pot = forward(quartet, three_sources, F(1))
    assert two_m == three_m and two_pot == three_pot
    assert min(p for _, _, p in three_sources) == F(1, 4)
    assert min(min(x, 1 - x) for _, x, _ in three_sources) == F(1, 3)

    result = {"stage": 43, "status": "all_exact_assertions_passed",
              "scope": "Proof controls for conditional mass bounds; no interval optimizer or external validation",
              "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "random_seed": 430317, "mixed_case_count": len(records),
              "vertex_intervals_checked": strict_count,
              "one_source_separation_count": len(one_source),
              "mixed_cases": records, "one_source_cases": one_source,
              "counterexamples": {
                  "free_pendant_one_vs_three": {"equal_leaf_observations": True,
                                                 "one_pendant": F(13, 12), "three_pendant": F(1),
                                                 "three_min_abundance": F(1, 3), "three_min_fraction": F(1, 4)},
                  "fixed_pendant_two_vs_three": {"equal_full_potential_and_mass": True,
                                                  "pendant": F(1), "three_min_abundance": F(1, 4),
                                                  "three_min_fraction": F(1, 3)}}}
    if write:
        if OUT.exists():
            raise FileExistsError(f"Refusing to overwrite {OUT}; use --check")
        OUT.write_text(json.dumps(result, indent=2, default=str) + "\n")
    print(json.dumps({key: result[key] for key in ("status", "mixed_case_count", "vertex_intervals_checked",
                                                  "one_source_separation_count")}))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    main(write=not args.check)
