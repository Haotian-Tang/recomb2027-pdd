#!/usr/bin/env python3
"""Stage 289 bridge: stage-287 bank case <-> DecoDiPhy 1.5.0 command-line inputs.

The bank stores an unrooted edge list with exact rational lengths.  DecoDiPhy
reads a Newick tree and a two-column distance file in double precision, and
reports placements as (node label, x) where the placement lies at fraction x of
the parent edge of that node, measured from the node towards its parent.  This
module writes those inputs, records the label/edge dictionary, and maps the
output back onto bank edge indices.  No bank file is modified.
"""
from fractions import Fraction as F
import json
from pathlib import Path


def geometry(edges):
    adjacency = {}
    for index, (u, v, length) in enumerate(edges):
        adjacency.setdefault(u, {})[v] = (F(length), index)
        adjacency.setdefault(v, {})[u] = (F(length), index)
    return adjacency


def load_case(path):
    case = json.loads(Path(path).read_text())
    case['edges'] = [(u, v, F(length)) for u, v, length in case['edges']]
    case['observed'] = {int(r): F(value) for r, value in case['observed'].items()}
    return case


def choose_root(adjacency):
    """An internal vertex of degree >= 3: every bank edge is then a parent edge."""
    best = None
    for v, nb in adjacency.items():
        if len(nb) >= 3 and (best is None or v < best):
            best = v
    if best is None:
        raise ValueError('No vertex of degree >= 3')
    return best


def label_of(v, leaves):
    return ('L%d' if v in leaves else 'N%d') % v


def rooted_newick(edges, root):
    """Newick with unique non-numeric labels; returns (text, parent, label maps)."""
    adjacency = geometry(edges)
    leaves = {v for v in adjacency if len(adjacency[v]) == 1}
    parent, order, seen = {root: None}, [root], {root}
    stack = [root]
    while stack:
        v = stack.pop()
        for u in sorted(adjacency[v]):
            if u not in seen:
                seen.add(u)
                parent[u] = v
                order.append(u)
                stack.append(u)
    children = {v: [] for v in adjacency}
    for v, p in parent.items():
        if p is not None:
            children[p].append(v)
    text = {}
    for v in reversed(order):
        label = label_of(v, leaves)
        inner = label if not children[v] else (
            '(' + ','.join(text[c] for c in children[v]) + ')' + label)
        if parent[v] is None:
            text[v] = inner + ';'
        else:
            length = adjacency[v][parent[v]][0]
            text[v] = '%s:%r' % (inner, float(length))
    label_to_vertex = {label_of(v, leaves): v for v in adjacency}
    edge_of_label = {label_of(v, leaves): adjacency[v][parent[v]][1]
                     for v in adjacency if parent[v] is not None}
    return text[root], parent, label_to_vertex, edge_of_label, leaves


def write_inputs(case, directory):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    adjacency = geometry(case['edges'])
    root = choose_root(adjacency)
    newick, parent, label_to_vertex, edge_of_label, leaves = rooted_newick(case['edges'], root)
    (directory / 'tree.nwk').write_text(newick + '\n')
    lines = ['%s\t%r' % (label_of(r, leaves), float(value))
             for r, value in sorted(case['observed'].items())]
    (directory / 'distances.txt').write_text('\n'.join(lines) + '\n')
    mapping = dict(root=root, edge_of_label=edge_of_label,
                   label_to_vertex=label_to_vertex,
                   parent={str(k): v for k, v in parent.items()})
    (directory / 'mapping.json').write_text(json.dumps(mapping))
    return mapping


def to_bank_sources(case, anchors, p_values, x_values, mapping):
    """DecoDiPhy (label, x) -> bank {edge, p, x} with x measured from edges[e][0]."""
    sources = []
    for label, p, x in zip(anchors, p_values, x_values):
        index = mapping['edge_of_label'][label]
        u, v, _ = case['edges'][index]
        vertex = mapping['label_to_vertex'][label]
        x_bank = float(x) if vertex == u else 1.0 - float(x)
        sources.append(dict(edge=index, p=float(p), x=x_bank, label=label, x_decodiphy=float(x)))
    return sources


def vertex_distances(edges, start):
    adjacency = geometry(edges)
    distance = {start: F(0)}
    stack = [start]
    while stack:
        v = stack.pop()
        for u, (length, _) in adjacency[v].items():
            if u not in distance:
                distance[u] = distance[v] + length
                stack.append(u)
    return distance


def point_distance_table(edges):
    """All-pairs vertex distances as Fractions (small trees only)."""
    return {v: vertex_distances(edges, v) for v in geometry(edges)}


def _frac(value):
    """Exact Fraction: Fraction, rational string, int, or binary-exact float."""
    return value if isinstance(value, F) else F(value)


def source_point_distance(edges, table, a, b):
    """Tree distance between two points, each given as {edge, x} with x from edges[e][0]."""
    ua, va, la = edges[a['edge']]
    ub, vb, lb = edges[b['edge']]
    xa, xb = _frac(a['x']), _frac(b['x'])
    if a['edge'] == b['edge']:
        return abs(xa - xb) * la
    return min(table[p][q] + sa + sb
               for p, sa in ((ua, xa * la), (va, (1 - xa) * la))
               for q, sb in ((ub, xb * lb), (vb, (1 - xb) * lb)))


def residual_vector(edges, sources, pendant, observed):
    """Model minus observed at every observed leaf, in floating point."""
    table = {}
    for s in sources:
        for v in edges[s['edge']][:2]:
            if v not in table:
                table[v] = vertex_distances(edges, v)
    out = {}
    for r, value in observed.items():
        total = _frac(pendant)
        for s in sources:
            u, v, length = edges[s['edge']]
            x = _frac(s['x'])
            total += _frac(s['p']) * min(
                table[u][r] + x * length, table[v][r] + (1 - x) * length)
        out[r] = float(total) - float(value)
    return out
