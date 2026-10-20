"""Exact-arithmetic Newick reader for PDD instance construction (stage 286/287).

Parses Newick strings into an undirected, integer-indexed tree whose edge
lengths are exact ``fractions.Fraction`` values obtained with
``Fraction(str)`` -- never through ``float`` -- so that every downstream
distance is an exact rational.

The reader is deliberately strict about what it accepts and reports
pathologies (negative or missing edge lengths, duplicate labels, unifurcations)
instead of silently repairing them.
"""
from fractions import Fraction
import re

_NUM = re.compile(r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?$')


def _frac(token):
    """Exact rational from a Newick length token (no float round trip)."""
    token = token.strip()
    if token == '':
        return None
    if not _NUM.match(token):
        raise ValueError('non-numeric edge length %r' % token)
    if 'e' in token or 'E' in token:
        mant, _, expo = token.lower().partition('e')
        return Fraction(mant) * Fraction(10) ** int(expo)
    return Fraction(token)


class Tree(object):
    """Undirected tree with integer vertex ids and exact rational edges."""

    __slots__ = ('parent', 'children', 'label', 'elen', 'root', 'n')

    def __init__(self):
        self.parent = []
        self.children = []
        self.label = []
        self.elen = []
        self.root = 0
        self.n = 0

    def _new(self, parent, label=None, elen=None):
        i = self.n
        self.n += 1
        self.parent.append(parent)
        self.children.append([])
        self.label.append(label)
        self.elen.append(elen)
        if parent is not None:
            self.children[parent].append(i)
        return i

    def leaves(self):
        return [i for i in range(self.n) if not self.children[i] and i != self.root] + \
               ([self.root] if not self.children[self.root] else [])

    def leaf_labels(self):
        return [self.label[i] for i in self.leaves()]


def parse_newick(text):
    """Parse one Newick tree. Returns (Tree, diagnostics dict)."""
    s = text.strip()
    if s.startswith('[') and ']' in s:          # leading [&R] / [&U] comment
        s = s[s.index(']') + 1:].strip()
    if not s.endswith(';'):
        s = s + ';'
    t = Tree()
    root = t._new(None)
    t.root = root
    cur = root
    i = 0
    token = []
    pending_label = True   # True when the next token belongs to the current node
    while i < len(s):
        ch = s[i]
        if ch == '(':
            cur = t._new(cur)
            i += 1
        elif ch == ',':
            cur = t._new(t.parent[cur])
            i += 1
        elif ch == ')':
            cur = t.parent[cur]
            i += 1
        elif ch == ';':
            break
        elif ch == ':':
            j = i + 1
            while j < len(s) and s[j] not in '(),;':
                j += 1
            t.elen[cur] = _frac(s[i + 1:j])
            i = j
        elif ch in ' \t\n\r':
            i += 1
        else:
            j = i
            while j < len(s) and s[j] not in '(),;:':
                j += 1
            lab = s[i:j].strip().strip("'\"")
            if lab:
                t.label[cur] = lab
            i = j
    # diagnostics
    diag = {}
    leaves = [i for i in range(t.n) if not t.children[i]]
    diag['n_vertices'] = t.n
    diag['n_leaves'] = len(leaves)
    diag['n_internal'] = t.n - len(leaves)
    diag['root_degree'] = len(t.children[root])
    diag['n_negative_edges'] = sum(1 for i in range(t.n)
                                  if i != root and t.elen[i] is not None and t.elen[i] < 0)
    diag['n_zero_edges'] = sum(1 for i in range(t.n)
                               if i != root and t.elen[i] == 0)
    diag['n_missing_edges'] = sum(1 for i in range(t.n)
                                  if i != root and t.elen[i] is None)
    diag['n_unifurcations'] = sum(1 for i in range(t.n) if len(t.children[i]) == 1)
    labs = [t.label[i] for i in leaves]
    diag['n_unlabelled_leaves'] = sum(1 for l in labs if not l)
    diag['n_duplicate_leaf_labels'] = len(labs) - len(set(labs))
    diag['n_labelled_internal'] = sum(1 for i in range(t.n)
                                      if t.children[i] and t.label[i])
    return t, diag


def leaf_label_set(text):
    t, d = parse_newick(text)
    return set(t.label[i] for i in range(t.n) if not t.children[i]), d
