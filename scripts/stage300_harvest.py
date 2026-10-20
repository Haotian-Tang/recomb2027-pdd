#!/usr/bin/env python3
"""Stage 300 step 1: harvest existing full-panel DecoDiPhy answers.

Stages 289, 293, 294 and 297 recorded thousands of DecoDiPhy runs, each with the
configuration the tool returned (bank edge index, abundance p, position x,
scalar pendant).  Re-running them would change nothing and cost hours, so the
full-panel runs are harvested here, with their provenance (file, line, arm,
repetition or seed) attached to every row.  Nothing is recomputed in this step.

Full panel means: stage 289 (every run is at the full reference panel), and
`n_removed == 0` for the removal-draw files of stages 293/294.  Stage 297's
seed study is full-panel by construction.  Runs whose status is not `complete`
are harvested too - a crash is a result.
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import stage300_common as s300                   # noqa: E402

RESULTS = ROOT / 'results'
SOURCES = [
    dict(stage=289, path='stage289/stage289_headtohead_runs.jsonl', style='stage289'),
    dict(stage=293, path='stage293/stage293_decodiphy_runs.jsonl', style='draws'),
    dict(stage=294, path='stage294/stage294_dd_runs.jsonl', style='draws'),
    dict(stage=297, path='stage297/stage297_dd_seeds_main.jsonl', style='stage297'),
]
FAMILY_TO_COLLECTION = {v: k for k, v in s300.FAMILY_OF.items()}


def harvest():
    rows = []
    for src in SOURCES:
        path = RESULTS / src['path']
        digest = s300.sha(path)
        with open(path) as handle:
            for line_no, line in enumerate(handle, 1):
                d = json.loads(line)
                family = d.get('family')
                collection = FAMILY_TO_COLLECTION.get(family)
                if collection is None:
                    continue                      # exact families: not this stage's object
                if src['style'] == 'stage289':
                    arm, full, tag = d['method'], True, dict(repetition=d.get('repetition'))
                elif src['style'] == 'draws':
                    arm, full = d.get('arm'), d.get('n_removed', 0) == 0
                    tag = dict(replicate=d.get('replicate'), regime=d.get('regime'),
                               level=d.get('level'))
                else:
                    arm, full, tag = d.get('arm'), True, dict(seed=d.get('seed'))
                if not full or arm is None or str(arm).startswith('exact'):
                    continue
                result = d.get('result', d)
                rows.append(dict(
                    kind='dd_harvest', stage=300, collection=collection, family=family,
                    name=d['name'], arm=arm, source_stage=src['stage'],
                    source_file=src['path'], source_sha256=digest, source_line=line_no,
                    tag=tag, status=result.get('status'), k=result.get('k'),
                    reported_loss=result.get('loss'),
                    sources=result.get('sources'),
                    pendant_total=result.get('pendant_total'),
                    stop_reason=result.get('stop_reason'),
                    ladder=result.get('ladder'),
                    jplace_writer_crashed=result.get('jplace_writer_crashed'),
                    writer_exception=result.get('writer_exception'),
                    unmapped_anchors=result.get('unmapped_anchors'),
                    exception=result.get('exception'),
                    argv_tail=[a for a in (result.get('argv') or []) if not a.startswith('/')],
                    solver_wall_seconds=result.get('solver_wall_seconds')))
    return rows


def main():
    out = s300.OUT
    out.mkdir(parents=True, exist_ok=True)
    rows = harvest()
    path = out / 'stage300_harvest.jsonl'
    with open(path, 'w') as handle:
        for row in rows:
            handle.write(json.dumps(s300.serial(row)) + '\n')
    summary = {}
    for row in rows:
        key = '%s|%s|%d' % (row['collection'], row['arm'], row['source_stage'])
        cell = summary.setdefault(key, dict(runs=0, complete=0, cases=set(),
                                            complete_cases=set(), crashed_writer=0))
        cell['runs'] += 1
        cell['cases'].add(row['name'])
        if row['status'] == 'complete' and row['sources'] is not None:
            cell['complete'] += 1
            cell['complete_cases'].add(row['name'])
        if row.get('jplace_writer_crashed'):
            cell['crashed_writer'] += 1
    for key in sorted(summary):
        cell = summary[key]
        print('%-40s runs=%4d complete=%4d cases=%3d complete_cases=%3d writer_crash=%3d'
              % (key, cell['runs'], cell['complete'], len(cell['cases']),
                 len(cell['complete_cases']), cell['crashed_writer']))
    print('harvested %d rows -> %s' % (len(rows), path))


if __name__ == '__main__':
    main()
