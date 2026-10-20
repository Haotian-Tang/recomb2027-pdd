#!/usr/bin/env python3
"""Stage295 shared loading: frozen Stage287 bank cases with hash verification.

Follows the Stage288 loader convention exactly (integer leaf keys, Fraction
lengths and observations, registered SHA256 verified on every read).
"""
import hashlib, json, os
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
INDEX = os.path.join(ROOT, 'data/stage287_instances/bank/index.json')


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def bank_index():
    return json.load(open(INDEX))['cases']


_BY_PATH = None


def load_case(relpath):
    """Read one frozen bank case by its registered relative path."""
    global _BY_PATH
    if _BY_PATH is None:
        _BY_PATH = {e['path']: e for e in bank_index()}
    entry = _BY_PATH[relpath]
    path = os.path.join(ROOT, relpath)
    digest = sha(path)
    assert digest == entry['sha256'], 'bank case hash changed: %s' % entry['name']
    raw = json.loads(open(path).read())
    edges = [(u, v, F(length)) for u, v, length in raw['edges']]
    observed = {int(r): F(d) for r, d in raw['observed'].items()}
    return edges, observed, entry, raw
