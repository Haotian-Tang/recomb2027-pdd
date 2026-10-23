#!/usr/bin/env python3
"""Stage304: attaining-edge SETS, not the canonical representative.

The closed-form one-source routine returns a per-edge optimum; the recorded
`argmin` is whichever edge the scan reached first among the optima. Asking
whether reference loss moved the best position is therefore two questions:

  did the canonical representative change,          argmin(R) != argmin(R')
  was the move forced,                              A(R) cap A(R') = empty

with A(R) = {e : theta_R({e}) = Theta_1(R)}. They differ: ties are common.
Stage299 stores `attaining` truncated to the first 200 entries alongside the
untruncated count `n_attaining`, so for a draw with more than 200 optima the
intersection cannot always be decided from the stream alone. For each unresolved
comparison, this script recomputes the draw's full attaining set with the frozen
closed form and intersects it with the full panel's complete saved set. The
results and per-case counts are recorded in `resolved_undecidable`; a truncated
full-panel set cannot be resolved by this procedure.

Truncation would in principle make the terminal-edge count a lower bound, since
a terminal optimum beyond entry 200 would be invisible. It does not here: every
truncated list already contains a terminal witness within its saved prefix, so
no record is left undecided and the count is exact. The cross-classification
that establishes this is `terminal_classification`.

Run from the repository root:  python3 -B scripts/stage304_attaining_sets.py
Writes results/stage304/stage304_attaining_sets.json. No third-party imports.
"""
import hashlib
import json
import sys
from fractions import Fraction
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUNS = ROOT / 'results' / 'stage299' / 'stage299_td_runs.jsonl'
OUT = ROOT / 'results' / 'stage304' / 'stage304_attaining_sets.json'

# Resolve each undecidable comparison by recomputing all draw-edge optima with
# the frozen closed form, then intersecting the full attaining set with the
# complete saved full-panel set. An empty intersection means the move was forced.
# Requires the stage299 protocol and case files; no LP is used.
RESOLVE = True


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()



def resolve_cell(rows, full, regime, level):
    """Resolve truncated draw lists via full closed-form attaining sets."""
    import sys as _sys
    _sys.path.insert(0, str(ROOT / 'scripts'))
    import stage298_common as s298
    import stage299_common as s299

    proto_path = ROOT / 'results' / 'stage299' / 'stage299_protocol.json'
    if not proto_path.exists():
        return None
    proto = json.loads(proto_path.read_text())
    cases = {}
    for row in proto['cases']:
        if row['arm'] != s299.TD_ARM:
            continue
        cases[row['name']] = s299.prepare_case(
            s299.load_json(row['path'], row['sha256']), row['name'])
    sets = {}
    for row in proto['leafsets']:
        if row['arm'] != s299.TD_ARM:
            continue
        sets[(row['name'], str(row['level']), row['regime'])] = row

    per_case, disjoint = [], 0
    for r in rows:
        if r['regime'] != regime or str(r['level']) != str(level):
            continue
        f = full.get(r['name'])
        if f is None:
            continue
        a_full = set(f.get('attaining') or [])
        a_new = set(r['result'].get('attaining') or [])
        trunc = ((f.get('n_attaining') or 0) > len(a_full) or
                 (r['result'].get('n_attaining') or 0) > len(a_new))
        if (a_full & a_new) or not trunc:
            continue
        if (f.get('n_attaining') or 0) > len(a_full):
            per_case.append(dict(name=r['name'], decided=False,
                                 why='full-panel list truncated too'))
            continue
        row = sets.get((r['name'], str(level), regime))
        case = cases.get(r['name'])
        if row is None or case is None:
            per_case.append(dict(name=r['name'], decided=False, why='leaf set not found'))
            continue
        removed = set(int(x) for x in row['removed'])
        leaves = [q for q in sorted(case['observed']) if q not in removed]
        dvals = [case['observed'][q] for q in leaves]
        cf = s298.theta1_all_edges(case['edges'], leaves, dvals)
        a_draw_full = set(s299.attaining(cf['theta'], cf['theta1']))
        hits = sorted(a_full & a_draw_full)
        per_case.append(dict(name=r['name'], decided=True, old_optima=len(a_full),
                             still_optimal=len(hits), disjoint=not hits))
        disjoint += int(not hits)
    return dict(checked=sum(1 for p in per_case if p['decided']), disjoint=disjoint,
                intersecting=sum(1 for p in per_case if p['decided'] and not p['disjoint']),
                per_case=per_case)


def main():
    if not RUNS.exists():
        print('missing input (fetch the bulk archive, see README):',
              RUNS.relative_to(ROOT))
        return 1
    rows = [json.loads(line) for line in RUNS.open()]
    full = {r['name']: r['result'] for r in rows
            if r['regime'] == 'none' and str(r['level']) == 'full'}

    n_att = [r['result']['n_attaining'] for r in rows
             if r['result'].get('n_attaining') is not None]
    truncated = sum(1 for r in rows
                    if isinstance(r['result'].get('attaining'), list)
                    and (r['result'].get('n_attaining') or 0) > len(r['result']['attaining']))

    out = {
        'stage': 304,
        'what': 'attaining-edge sets against the canonical representative',
        'definitions': {
            'representative_changed': "argmin(R) != argmin(R')",
            'move_forced': "the attaining sets are disjoint",
            'undecidable': 'at least one attaining list is truncated at 200 entries',
        },
        'input': {str(RUNS.relative_to(ROOT)): sha256(RUNS)},
        'leaf_sets': len(rows),
        'full_panels': len(full),
        'with_tied_optima': sum(1 for x in n_att if x > 1),
        'max_tied_optima': max(n_att),
        'attaining_lists_truncated': truncated,
        'representative_is_terminal': sum(1 for r in rows
                                          if r['result'].get('argmin_is_terminal')),
        'attaining_set_contains_terminal': sum(1 for r in rows
                                               if r['result'].get('attaining_terminal')),
        'by_regime_level': {},
    }

    tc = defaultdict(int)
    for r in rows:
        res = r['result']
        att = res.get('attaining') or []
        trunc = (res.get('n_attaining') or 0) > len(att)
        is_full = (r['regime'] == 'none' and str(r['level']) == 'full')
        tc['%s|%s|%s' % ('full' if is_full else 'draw',
                         'truncated' if trunc else 'complete',
                         'terminal' if res.get('attaining_terminal') else 'none')] += 1
    out['terminal_classification'] = dict(sorted(tc.items()))
    out['terminal_undecided'] = (tc.get('draw|truncated|none', 0)
                                 + tc.get('full|truncated|none', 0))

    cells = defaultdict(lambda: dict(n=0, representative_changed=0, intersecting=0,
                                     disjoint=0, undecidable=0))
    for r in rows:
        if r['regime'] == 'none':
            continue
        f = full.get(r['name'])
        if f is None:
            continue
        key = '%s|%s' % (r['regime'], r['level'])
        c = cells[key]
        c['n'] += 1
        a_full = set(f.get('attaining') or [])
        a_new = set(r['result'].get('attaining') or [])
        trunc = ((f.get('n_attaining') or 0) > len(a_full) or
                 (r['result'].get('n_attaining') or 0) > len(a_new))
        if f.get('argmin') != r['result'].get('argmin'):
            c['representative_changed'] += 1
        if a_full & a_new:
            c['intersecting'] += 1
        elif trunc:
            c['undecidable'] += 1
        else:
            c['disjoint'] += 1
    out['resolved_undecidable'] = {}
    for key in sorted(cells):
        c = cells[key]
        if RESOLVE and c['undecidable']:
            regime, level = key.split('|')
            res = resolve_cell(rows, full, regime, level)
            if res is not None:
                out['resolved_undecidable'][key] = res
                c['disjoint_after_resolution'] = c['disjoint'] + res['disjoint']
        out['by_regime_level'][key] = c

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=1) + '\n')
    print('wrote', OUT.relative_to(ROOT))
    print('  %d leaf sets, %d with tied optima (max %d), %d lists truncated'
          % (out['leaf_sets'], out['with_tied_optima'], out['max_tied_optima'],
             out['attaining_lists_truncated']))
    print('  representative terminal in %d, attaining set contains a terminal edge in %d '
          '(%d undecided)' % (out['representative_is_terminal'],
                              out['attaining_set_contains_terminal'],
                              out['terminal_undecided']))
    print('  terminal classification:', out['terminal_classification'])
    for key, c in out['by_regime_level'].items():
        print('  %-24s n=%3d  rep changed %3d  intersect %3d  disjoint %3d%s  undecidable %3d'
              % (key, c['n'], c['representative_changed'], c['intersecting'], c['disjoint'],
                 ('->%d' % c['disjoint_after_resolution']) if 'disjoint_after_resolution' in c else '',
                 c['undecidable']))
    return 0


if __name__ == '__main__':
    sys.exit(main())
