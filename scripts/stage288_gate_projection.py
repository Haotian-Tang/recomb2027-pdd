#!/usr/bin/env python3
"""Stage 288 addendum: does a true source actually project onto a MARKED gate?

A Lemma-3 count gap needs more than a non-empty mark set: some source must
project onto a marked branch gate, because only then can a mixed explanation
require a vertex atom that the mark-free reduced problem cannot express. This
is pure combinatorics on the frozen bank - no solver is called. It uses exactly
the removal draws registered in the stage 288 incompleteness protocol.
"""
import json
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from fractions import Fraction as F

import stage288_common as common
from stage288_common import ROOT, serial, sha, load_case
from stage288_incompleteness import removal_draws, LEVELS, PROTOCOL
from stage26_missing_leaf import observed_backbone

OUT = ROOT / 'results/stage288_gate_projection.json'


def classify(case, retained):
    """For each generating source, locate its gate in the observed hull."""
    backbone = observed_backbone(case['edges'], sorted(retained))
    hull = set(backbone['hull_vertices'])
    marked = {v['vertex']: v for v in backbone['saturation_violations']}
    interior = backbone['interior_vertex']
    on_backbone = backbone['original_to_backbone']
    rows = []
    for source in case['generating_sources']:
        index = source['edge']
        u, v, length = case['edges'][index]
        if index in on_backbone:
            part = on_backbone[index]
            x = source['x'] if u == part['u'] else 1 - source['x']
            coordinate = part['start'] + x * length
            at_gate = None
            for other, info in interior.items():
                if info['edge'] == part['backbone_edge'] and info['coordinate'] == coordinate:
                    at_gate = other
                    break
            rows.append(dict(position='backbone_edge_interior' if at_gate is None
                             else 'suppressed_degree_two_gate',
                             gate=at_gate, marked=False, gate_original_degree=None))
        else:
            gate = backbone['gate'][u]
            assert backbone['gate'][v] == gate
            rows.append(dict(position='off_hull', gate=gate, marked=gate in marked,
                             gate_original_degree=marked[gate]['original_degree'] if gate in marked
                             else (len(common.adjacency_lists(case['edges'])[gate])
                                   if gate in hull else None),
                             gate_in_hull=gate in hull,
                             gate_is_observed_leaf=gate in retained))
    return rows, len(marked)


def main():
    index = [e for e in common.bank_index() if e['family'].startswith('exact')]
    started = time.perf_counter()
    draws, per_case = [], []
    for position, entry in enumerate(index):
        case = load_case(entry)
        row = dict(name=case['name'], family=case['family'], N=case['reference_set_size'],
                   K=case['K_generating'], draws=0, draws_with_a_source_at_a_marked_gate=0,
                   max_sources_at_marked_gates=0)
        for level in LEVELS:
            for draw in removal_draws(case, level):
                retained = set(case['observed']) - set(draw['removed'])
                rows, n_marks = classify(case, retained)
                at_marked = sum(r['marked'] for r in rows)
                record = dict(name=case['name'], N=row['N'], K=row['K'], level=level,
                              regime=draw['regime'], replicate=draw['replicate'],
                              n_removed=draw['n_removed'], n_retained=draw['n_retained'],
                              n_marks=n_marks, sources_at_marked_gates=at_marked,
                              marked_gate_degrees=sorted(r['gate_original_degree'] for r in rows
                                                         if r['marked']),
                              positions=dict(Counter(r['position'] for r in rows)))
                draws.append(record)
                row['draws'] += 1
                row['draws_with_a_source_at_a_marked_gate'] += int(at_marked > 0)
                row['max_sources_at_marked_gates'] = max(row['max_sources_at_marked_gates'], at_marked)
        per_case.append(row)
        if (position + 1) % 200 == 0:
            print(f'  {position + 1}/{len(index)} ({time.perf_counter() - started:.0f}s)', flush=True)
    summary = dict(
        exact_cases=len(per_case), total_draws=len(draws),
        draws_with_non_empty_mark_set=sum(d['n_marks'] > 0 for d in draws),
        draws_with_a_source_at_a_marked_gate=sum(d['sources_at_marked_gates'] > 0 for d in draws),
        cases_with_a_source_at_a_marked_gate=sum(r['draws_with_a_source_at_a_marked_gate'] > 0
                                                 for r in per_case),
        source_position_totals=dict(sum((Counter(d['positions']) for d in draws), Counter())),
        marked_gate_degrees_hit=dict(sorted(Counter(
            g for d in draws for g in d['marked_gate_degrees']).items())),
        by_regime={g: dict(draws=sum(d['regime'] == g for d in draws),
                           with_source_at_marked_gate=sum(d['regime'] == g
                                                          and d['sources_at_marked_gates'] > 0
                                                          for d in draws))
                   for g in sorted({d['regime'] for d in draws})},
        by_level={l: dict(draws=sum(d['level'] == l for d in draws),
                          with_source_at_marked_gate=sum(d['level'] == l
                                                         and d['sources_at_marked_gates'] > 0
                                                         for d in draws))
                  for l in LEVELS},
        marked_case_draws=[d for d in draws if d['n_marks'] > 0])
    OUT.write_text(json.dumps(serial(dict(
        stage=288, track='gate_projection', created_utc=datetime.now(timezone.utc).isoformat(),
        protocol_sha256=sha(PROTOCOL), summary=summary, per_case=per_case)), indent=2) + '\n')
    printable = dict(summary)
    printable['marked_case_draws'] = printable['marked_case_draws'][:6]
    print(json.dumps(serial(printable), indent=2)[:6000])


if __name__ == '__main__':
    main()
