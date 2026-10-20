#!/usr/bin/env python3
"""Exact positive colored-SUM reduction on a strictly binary tree.

The extra calibration arm absorbs the small positive binary-core deformation.
This is a constructive candidate proof implementation, not a priority claim.
No previous-stage inference, construction, LP, or tree utility is imported.

For k groups, c=1/[2(k+B)], tau=c/16, q=1/2-tau, Delta=1/2-tau.
The root has unit sentinel and calibration edges and a positive binary fan.
Every active arm has length one and one unit hidden edge at each reciprocal
gate 1/(1+a). Total fan length is epsilon=c/[100*k*(Amax+2)].
The source budget is k+1, prescribed pendant one. Observations are
d_z=3+k*c-2*tau, d_cal=d_z-2*q, and d_i=d_z-2*c.

If h_i is the root-to-port distance and B_ij is the shared root path length,
one source at active gate s_i satisfies (diag(s)+B)p=c*1+h/2. Solve that
system exactly, put p_cal=1-sum(p), x_cal=q/p_cal, and choose the common
hidden-edge position to match the absolute sentinel distance. The same
construction works with fractional group means and at most two atoms in
one group, supplying a k+2-source promise witness without solving SUM.
"""
from fractions import Fraction as F
from itertools import product


def serial(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {str(k): serial(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [serial(v) for v in value]
    return value


def solve_linear(matrix, rhs):
    """Independent exact rational Gauss-Jordan solve of a square system."""
    n = len(rhs)
    a = [[F(x) for x in row] + [F(b)] for row, b in zip(matrix, rhs)]
    assert len(a) == n and all(len(row) == n + 1 for row in a)
    for column in range(n):
        pivot = next((r for r in range(column, n) if a[r][column]), None)
        if pivot is None:
            raise ValueError('Singular square system')
        a[column], a[pivot] = a[pivot], a[column]
        scale = a[column][column]
        a[column] = [x / scale for x in a[column]]
        for row in range(n):
            if row != column and a[row][column]:
                scale = a[row][column]
                a[row] = [x - scale * y for x, y in zip(a[row], a[column])]
    return [row[-1] for row in a]


def generate_case(groups, target, name='binary_colored_sum', shape='caterpillar'):
    """Construct JSON-ready data; no SUM search or inference is performed."""
    if not isinstance(target, int) or isinstance(target, bool) or target <= 0:
        raise ValueError('Positive integer target required')
    if len(groups) < 2:
        raise ValueError('At least two groups required')
    normalized = []
    for group in groups:
        if not group or any(not isinstance(a, int) or isinstance(a, bool) or a <= 0
                            for a in group):
            raise ValueError('Nonempty positive integer groups required')
        normalized.append(sorted(set(group)))
    if shape not in ('caterpillar', 'balanced'):
        raise ValueError('shape must be caterpillar or balanced')
    k = len(normalized)
    amax = max(map(max, normalized))
    c = F(1, 2 * (k + target))
    tau, q, delta = c / 16, F(1, 2) - c / 16, F(1, 2) - c / 16
    epsilon = c / (100 * k * (amax + 2))
    edge_length = epsilon / (k - 1)
    edges = [('root', 'sentinel', F(1)), ('root', 'calibration', F(1))]
    ports = [None] * k
    paths = [None] * k
    core_count = 0

    def fan(parent, indices, path):
        nonlocal core_count
        if len(indices) == 1:
            ports[indices[0]], paths[indices[0]] = parent, path
            return
        vertex = f'core{core_count}'
        core_count += 1
        edge_id = len(edges)
        edges.append((parent, vertex, edge_length))
        split = 1 if shape == 'caterpillar' else len(indices) // 2
        fan(vertex, indices[:split], path + [edge_id])
        fan(vertex, indices[split:], path + [edge_id])

    fan('root', list(range(k)), [])
    assert core_count == k - 1
    h = [len(path) * edge_length for path in paths]
    overlap = [[len(set(left) & set(right)) * edge_length for right in paths]
               for left in paths]
    options = []
    for i, group in enumerate(normalized):
        previous, position = ports[i], F(0)
        arm_options = []
        for value in sorted(group, reverse=True):
            coordinate = F(1, 1 + value)
            gate, hidden = f'g{i}_{value}', f'z{i}_{value}'
            edges.append((previous, gate, coordinate - position))
            hidden_edge = len(edges)
            edges.append((gate, hidden, F(1)))
            arm_options.append(dict(value=value, gate=gate, position=coordinate,
                                    hidden_edge=hidden_edge, hidden_length=F(1)))
            previous, position = gate, coordinate
        edges.append((previous, f'r{i}', 1 - position))
        options.append(arm_options)
    d_sentinel = 3 + k * c - 2 * tau
    observed = {f'r{i}': d_sentinel - 2 * c for i in range(k)}
    observed['sentinel'] = d_sentinel
    observed['calibration'] = d_sentinel - 2 * q
    case = dict(name=name, edges=edges, observed=observed, pendant=F(1), cap=k + 1,
                reduction=dict(groups=normalized, target=target, k=k, amax=amax,
                               c=c, tau=tau, q=q, delta=delta, epsilon=epsilon,
                               fan_edge_length=edge_length, shape=shape, ports=ports,
                               root_distances=h, overlap=overlap, options=options,
                               calibration_edge=1, original_edges=len(edges),
                               observed_leaves=k + 2, missing_leaves=sum(map(len, normalized))))
    validate_tree(case)
    return serial(case)


def validate_tree(case):
    adj = {}
    for u, v, length in case['edges']:
        assert u != v and F(length) > 0
        assert v not in adj.setdefault(u, {})
        adj[u][v] = F(length)
        adj.setdefault(v, {})[u] = F(length)
    assert len(case['edges']) == len(adj) - 1
    assert all(len(neighbors) in (1, 3) for neighbors in adj.values())
    reached, queue = {next(iter(adj))}, [next(iter(adj))]
    for u in queue:
        for v in adj[u]:
            if v not in reached:
                reached.add(v)
                queue.append(v)
    assert reached == set(adj)
    assert all(len(adj[v]) == 1 for v in case['observed'])
    return adj


def sum_solutions(groups, target):
    return [list(values) for values in product(*groups) if sum(values) == target]


def _witness_from_choices(case, choices):
    info = case['reduction']
    k, c = info['k'], F(info['c'])
    h = list(map(F, info['root_distances']))
    overlap = [[F(x) for x in row] for row in info['overlap']]
    means = [sum(F(weight) * value for value, weight in atoms) for atoms in choices]
    assert len(means) == k and sum(means) == info['target']
    effective = [1 / (1 + mean) for mean in means]
    matrix = [[overlap[i][j] + (effective[i] if i == j else 0)
               for j in range(k)] for i in range(k)]
    p = solve_linear(matrix, [c + value / 2 for value in h])
    total = sum(p)
    calibration_mass = 1 - total
    calibration_position = F(info['q']) / calibration_mass
    hidden_moment = F(info['delta']) + k * c - sum(
        mass * (height + coordinate) for mass, height, coordinate in zip(p, h, effective))
    theta = hidden_moment / total
    assert all(mass > 0 for mass in p)
    assert calibration_mass > 0 and 0 < calibration_position < 1 and 0 < theta < 1
    sources = []
    for i, atoms in enumerate(choices):
        assert sum(F(weight) for _, weight in atoms) == 1
        for value, weight in atoms:
            weight = F(weight)
            assert weight > 0 and value in info['groups'][i]
            option = next(o for o in info['options'][i] if o['value'] == value)
            mass = p[i] * weight * (1 + value) / (1 + means[i])
            sources.append(dict(edge=option['hidden_edge'], p=mass, x=theta))
    sources.append(dict(edge=info['calibration_edge'], p=calibration_mass,
                        x=calibration_position))
    record = dict(sources=sources, convex_choices=choices, means=means,
                  source_count=len(sources), pendant=case['pendant'], active_masses=p,
                  calibration_mass=calibration_mass, calibration_position=calibration_position,
                  hidden_moment=hidden_moment, common_hidden_position=theta,
                  ideal_active_masses=[c * (1 + mean) for mean in means],
                  active_mass_l1_error=sum(abs(mass - c * (1 + mean))
                                           for mass, mean in zip(p, means)))
    replay(case, sources)
    return serial(record)


def witness_for_solution(case, values):
    """Return a strict k+1-source witness for the specified discrete solution."""
    info = case['reduction']
    if len(values) != info['k'] or sum(values) != info['target']:
        raise ValueError('Not a solution with one value per group')
    if any(value not in group for value, group in zip(values, info['groups'])):
        raise ValueError('Value absent from its group')
    return _witness_from_choices(case, [[(value, F(1))] for value in values])['sources']


def promised_witness(case):
    """Polynomial strict k+2 witness when B is in the continuous sum interval."""
    info = case['reduction']
    groups, target = info['groups'], info['target']
    low, high = sum(map(min, groups)), sum(map(max, groups))
    if not low <= target <= high:
        raise ValueError('Target outside the continuous sum interval')
    remaining = F(target - low)
    choices = []
    for group in groups:
        lo, hi = min(group), max(group)
        amount = min(remaining, F(hi - lo))
        remaining -= amount
        mean = lo + amount
        if mean == lo:
            atoms = [(lo, F(1))]
        elif mean == hi:
            atoms = [(hi, F(1))]
        else:
            weight = (mean - lo) / (hi - lo)
            atoms = [(lo, 1 - weight), (hi, weight)]
        choices.append(atoms)
    assert remaining == 0
    result = _witness_from_choices(case, choices)
    assert result['source_count'] <= case['cap'] + 1
    return result


def replay(case, sources):
    """Direct original-tree rational distance replay and strict binary checks."""
    adj = validate_tree(case)
    assert len({source['edge'] for source in sources}) == len(sources)
    assert sum(F(source['p']) for source in sources) == 1
    assert all(F(source['p']) > 0 and 0 < F(source['x']) < 1 for source in sources)
    predictions = {}
    for ref in case['observed']:
        distances, queue = {ref: F(0)}, [ref]
        for vertex in queue:
            for neighbor, length in adj[vertex].items():
                if neighbor not in distances:
                    distances[neighbor] = distances[vertex] + length
                    queue.append(neighbor)
        prediction = F(case['pendant'])
        for source in sources:
            u, v, length = case['edges'][source['edge']]
            length, x, mass = F(length), F(source['x']), F(source['p'])
            prediction += mass * min(distances[u] + x * length,
                                     distances[v] + (1 - x) * length)
        predictions[ref] = prediction
    assert predictions == {ref: F(value) for ref, value in case['observed'].items()}
    return serial(predictions)
