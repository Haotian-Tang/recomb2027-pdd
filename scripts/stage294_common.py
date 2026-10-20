#!/usr/bin/env python3
"""Stage 294 shared library: the stage 288 incompleteness design on NOISY data.

Stage 288/290/293 measured reference incompleteness on exactly explainable
instances in exact rational arithmetic.  This stage asks the same question of
the krepp_estimated_10k family, whose observed distances are real ANI estimates
from reads (median |model - observed| = 0.134) and which therefore has no exact
explanation at all.  The tool that does run on such data is DecoDiPhy 1.5.0.

Design rules IMPORTED, not re-derived:
  * levels, regimes, replicate counts, the removal-count rule and the three
    removal regimes come from stage288_incompleteness / stage288_common;
  * the random-draw seed key uses the PAIRED EXACT case name, so a random draw
    here is byte-identical to the stage 288 draw whenever the two panels agree
    (135 of 296 pairs; the rest differ only by the 1-3 references for which
    krepp returned no ANI estimate).

Added here: a fourth, taxonomic removal regime built from the Web of Life NCBI
rank table, which drops whole taxonomic groups at a rank rather than a clade.

Nothing pre-existing is read-write; every path written is new and under 294.
"""
import os
for _key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
             'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_key] = '1'
import collections
import hashlib
import json
import lzma
import re
import sys
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / 'scripts') not in sys.path:
    sys.path.insert(0, str(ROOT / 'scripts'))
import stage288_common as common                      # noqa: E402
import stage288_incompleteness as inc                 # noqa: E402
from stage288_common import sha, serial               # noqa: E402

BANK = ROOT / 'data/stage287_instances/bank'
TAXDIR = ROOT / 'data/stage294_taxonomy'
RANKS_FILE = TAXDIR / 'raw/taxonomy_ncbi_ranks.tsv.xz'
NOISY_FAMILY = 'krepp_estimated_10k'
EXACT_FAMILY = 'exact_10k'
TAXON_RANKS = ('genus', 'species')      # finest operable rank, in this order
TAXON_ORDERS = ('large_first', 'small_first')
UNASSIGNED = '<unassigned>'


# ------------------------------------------------------------------ bank I/O

def stem(name):
    """Strip the family suffix: the two families share one instance stem."""
    return re.sub(r'_(exact|krepp_estimated)$', '', name)


def bank_entries(family):
    return {e['name']: e for e in common.bank_index() if e['family'] == family}


def load_noisy_case(entry):
    """Read one krepp_estimated_10k case.

    The noisy schema carries no 'generating_pendant' key -- under ANI estimates
    there is no pendant that reproduces the data -- but it does carry 'pendant',
    the generating mean pendant of the paired exact instance.  It is recorded as
    published truth and never used as a solver input.
    """
    path = ROOT / entry['path']
    digest = sha(path)
    assert digest == entry['sha256'], f"bank case hash changed: {entry['name']}"
    raw = json.loads(path.read_text())
    assert 'generating_pendant' not in raw, 'noisy schema changed'
    return dict(
        name=entry['name'], family=entry['family'], path=entry['path'], sha256=digest,
        edges=[(u, v, F(length)) for u, v, length in raw['edges']],
        observed={int(r): F(d) for r, d in raw['observed'].items()},
        cap=raw['cap'],
        generating_sources=[dict(edge=s['edge'], p=F(s['p']), x=F(s['x']))
                            for s in raw['generating_sources']],
        published_pendant=F(raw['pendant']),
        leaf_labels={int(k): v for k, v in raw['leaf_labels'].items()},
        provenance=raw['provenance'],
        reference_set_size=entry['reference_set_size'],
        K_generating=entry['K_generating'])


def load_exact_case(entry):
    """The paired exactly explainable instance: same tree, same sources."""
    return common.load_case(entry)


def pair_panel():
    """Every krepp case with its paired exact case; the full family, no subsample."""
    noisy = bank_entries(NOISY_FAMILY)
    exact = bank_entries(EXACT_FAMILY)
    by_stem = {stem(n): e for n, e in exact.items()}
    out = []
    for name in sorted(noisy):
        entry = noisy[name]
        twin = by_stem[stem(name)]
        out.append(dict(noisy=entry, exact=twin, stem=stem(name)))
    assert len(out) == len(noisy)
    return out


# --------------------------------------------------------------- taxonomy

def load_ranks():
    """genome accession -> {rank: taxon name}, from the frozen WoL NCBI table."""
    table = {}
    with lzma.open(RANKS_FILE, 'rt') as handle:
        header = handle.readline().rstrip('\n').split('\t')
        for line in handle:
            fields = line.rstrip('\n').split('\t')
            table[fields[0]] = dict(zip(header[1:], fields[1:]))
    return table


def taxon_groups(ranks, labels, leaves, rank):
    groups = collections.defaultdict(list)
    for leaf in leaves:
        name = ranks[labels[leaf]][rank].strip() or UNASSIGNED
        groups[name].append(leaf)
    return {k: sorted(v) for k, v in groups.items()}


def taxon_removal(ranks, labels, leaves, count, order):
    """Drop WHOLE taxonomic groups at the finest operable rank.

    Groups are taken in a deterministic order -- largest first, or smallest
    first -- and a group is taken only if it does not overshoot the target count
    and leaves at least three references, exactly the stage 288 retention floor.
    Whole-group granularity means the achieved count can fall short of the
    target; both are recorded and the analysis bins on the achieved fraction as
    well as the nominal level.
    """
    n = len(leaves)
    for rank in TAXON_RANKS:
        groups = taxon_groups(ranks, labels, leaves, rank)
        if len(groups) < 2:
            continue
        if order == 'large_first':
            keys = sorted(groups, key=lambda k: (-len(groups[k]), k))
        else:
            keys = sorted(groups, key=lambda k: (len(groups[k]), k))
        chosen, removed = [], []
        for key in keys:
            size = len(groups[key])
            if len(removed) + size <= count and n - len(removed) - size >= 3:
                chosen.append(key)
                removed.extend(groups[key])
        if not chosen:
            continue
        return dict(removed=sorted(removed), rank=rank, groups=chosen,
                    n_groups_at_rank=len(groups), target=count, order=order)
    return None


def is_monophyletic(edges, leaves_all, subset):
    """Is `subset` exactly the leaf set on one side of a single edge?"""
    subset = set(subset)
    adjacency = collections.defaultdict(list)
    for index, (u, v, _) in enumerate(edges):
        adjacency[u].append((v, index))
        adjacency[v].append((u, index))
    allset = set(leaves_all)
    for index, (u, v, _) in enumerate(edges):
        seen = {u}
        stack = [u]
        while stack:
            x = stack.pop()
            for y, eidx in adjacency[x]:
                if eidx == index or y in seen:
                    continue
                seen.add(y)
                stack.append(y)
        side = seen & allset
        if side == subset or (allset - side) == subset:
            return True
    return False


# -------------------------------------------------------------- draw plan

def removal_draws(noisy_case, exact_case, level, ranks):
    """Every registered removal draw for one noisy case and one level.

    random / clade / nearest follow stage288_incompleteness.removal_draws with
    the seed key taken from the PAIRED EXACT case name; taxon is new here.
    """
    leaves = sorted(noisy_case['observed'])
    n = len(leaves)
    count = inc.removal_count(n, level)
    if count is None:
        return []
    draws = []
    for rep in range(inc.REPLICATES['random']):
        key = f"{inc.SEED}|{exact_case['name']}|{level}|random|{rep}"
        seed = int(hashlib.sha256(key.encode()).hexdigest()[:8], 16)
        draws.append(dict(regime='random', replicate=rep,
                          removed=common.random_removal(leaves, count, seed),
                          detail={'seed': seed, 'seed_key_case': exact_case['name']}))
    for rep, clade in enumerate(common.clade_removals(noisy_case['edges'], leaves, count,
                                                      inc.REPLICATES['clade'])):
        draws.append(dict(regime='clade', replicate=rep, removed=clade['removed'],
                          detail={'clade_vertex': clade['clade_vertex'],
                                  'complement': clade['complement']}))
    draws.append(dict(regime='nearest', replicate=0,
                      removed=common.nearest_source_removal(
                          noisy_case['edges'], noisy_case['generating_sources'], leaves, count),
                      detail={}))
    for rep, order in enumerate(TAXON_ORDERS):
        drawn = taxon_removal(ranks, noisy_case['leaf_labels'], leaves, count, order)
        if drawn is None:
            continue
        draws.append(dict(regime='taxon', replicate=rep, removed=drawn['removed'],
                          detail={'rank': drawn['rank'], 'groups': drawn['groups'],
                                  'order': order, 'target_count': count,
                                  'n_groups_at_rank': drawn['n_groups_at_rank']}))
    for draw in draws:
        draw['n_removed'] = len(draw['removed'])
        draw['n_retained'] = n - len(draw['removed'])
        draw['level'] = level
    return [d for d in draws if d['n_retained'] >= 3 and d['n_removed'] >= 1]


def draw_plan(ranks=None):
    """(noisy case, exact case, meta) for every draw, plus the full-panel baseline."""
    ranks = load_ranks() if ranks is None else ranks
    plan = []
    for pair in pair_panel():
        noisy = load_noisy_case(pair['noisy'])
        exact = load_exact_case(pair['exact'])
        leaves = sorted(noisy['observed'])
        exact_only = sorted(set(exact['observed']) - set(noisy['observed']))
        base = dict(name=noisy['name'], stem=pair['stem'], exact_name=exact['name'],
                    family=noisy['family'], reference_set_size=noisy['reference_set_size'],
                    K_generating=noisy['K_generating'], cap=noisy['cap'],
                    n_leaves=len(leaves),
                    exact_panel_extra_references=exact_only,
                    panels_identical=not exact_only)
        plan.append((noisy, exact, dict(base, level='none', regime='none', replicate=0,
                                        removed=[], n_removed=0, n_retained=len(leaves),
                                        removed_fraction=0.0, detail={})))
        for level in inc.LEVELS:
            for draw in removal_draws(noisy, exact, level, ranks):
                plan.append((noisy, exact, dict(
                    base, level=level, regime=draw['regime'], replicate=draw['replicate'],
                    removed=draw['removed'], n_removed=draw['n_removed'],
                    n_retained=draw['n_retained'],
                    removed_fraction=draw['n_removed'] / len(leaves),
                    detail=draw['detail'])))
    return plan


def retained_observations(case, removed):
    drop = set(removed)
    return {r: d for r, d in case['observed'].items() if r not in drop}


def plan_digest(plan):
    payload = [[m['name'], m['level'], m['regime'], m['replicate'], list(m['removed'])]
               for _, _, m in plan]
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()
