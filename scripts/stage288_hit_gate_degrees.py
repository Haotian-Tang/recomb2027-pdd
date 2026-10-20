#!/usr/bin/env python3
"""Stage 288 addendum: hull degree of the marked gates that a true source hits.

The manuscript four-leaf-star argument gives a count gap of 1 -> d where d is
the number of observed arms at the gate, i.e. its HULL degree. This records
that number for every real removal draw in which a generating source projects
onto a marked gate, so the predicted inflation on the real tree is quantified
even where the solver cannot confirm it within budget.
"""
import json
from collections import Counter
from datetime import datetime, timezone

import stage288_common as common
from stage288_common import ROOT, serial, sha, load_case
from stage288_incompleteness import removal_draws, LEVELS, PROTOCOL
from stage26_missing_leaf import observed_backbone

OUT = ROOT / 'results/stage288_hit_gate_degrees.json'


def main():
    index = [e for e in common.bank_index() if e['family'].startswith('exact')]
    census = json.loads((ROOT / 'results/stage288_gate_projection.json').read_text())
    wanted = {r['name'] for r in census['per_case'] if r['draws_with_a_source_at_a_marked_gate'] > 0}
    print('cases to re-examine:', sorted(wanted), flush=True)
    rows = []
    for entry in index:
        if entry['name'] not in wanted:
            continue
        case = load_case(entry)
        for level in LEVELS:
            for draw in removal_draws(case, level):
                retained = sorted(set(case['observed']) - set(draw['removed']))
                backbone = observed_backbone(case['edges'], retained)
                marked = {v['vertex']: v for v in backbone['saturation_violations']}
                if not marked:
                    continue
                hits = []
                for source in case['generating_sources']:
                    u, v, _ = case['edges'][source['edge']]
                    if source['edge'] in backbone['original_to_backbone']:
                        continue
                    gate = backbone['gate'][u]
                    if gate in marked:
                        hits.append(dict(edge=source['edge'], p=str(source['p']), gate=gate,
                                         hull_degree=marked[gate]['hull_degree'],
                                         original_degree=marked[gate]['original_degree'],
                                         predicted_gap_one_to=marked[gate]['hull_degree']))
                rows.append(dict(name=case['name'], N=case['reference_set_size'],
                                 K=case['K_generating'], level=level, regime=draw['regime'],
                                 replicate=draw['replicate'], n_removed=draw['n_removed'],
                                 n_retained=draw['n_retained'], n_marks=len(marked),
                                 sources_at_marked_gates=len(hits), hits=hits))
                print(rows[-1]['level'], rows[-1]['regime'], rows[-1]['replicate'],
                      'marks', len(marked), 'hits', len(hits),
                      [h['hull_degree'] for h in hits], flush=True)
    summary = dict(
        draws_examined=len(rows),
        draws_with_a_hit=sum(r['sources_at_marked_gates'] > 0 for r in rows),
        hit_hull_degree_histogram=dict(sorted(Counter(
            h['hull_degree'] for r in rows for h in r['hits']).items())),
        hit_original_degree_histogram=dict(sorted(Counter(
            h['original_degree'] for r in rows for h in r['hits']).items())),
        max_predicted_gap=max([h['hull_degree'] for r in rows for h in r['hits']], default=None),
        by_regime={g: sum(r['regime'] == g and r['sources_at_marked_gates'] > 0 for r in rows)
                   for g in ('random', 'clade', 'nearest')})
    OUT.write_text(json.dumps(serial(dict(
        stage=288, track='hit_gate_degrees', created_utc=datetime.now(timezone.utc).isoformat(),
        protocol_sha256=sha(PROTOCOL), summary=summary, rows=rows)), indent=2) + '\n')
    print(json.dumps(serial(summary), indent=2))


if __name__ == '__main__':
    main()
