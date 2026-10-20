#!/usr/bin/env python3
"""Stage 300 diagnostic: is the configuration returned the one whose loss is reported?

DecoDiPhy writes `all_rounds.json` and then, AFTER hill climbing has finished,
substitutes anchors whose position sits at an edge endpoint (cli.py: the
`max(opt_x) > 1 - 1e-3` and `min(opt_x) < 1e-4` blocks).  The substitution is
accepted when the new position is interior - `if new_x[i] < 1 - 1e-3` - without
requiring the objective not to worsen, and the final record copies the loss of
the previous round (`"loss": all_rounds[-1]["loss"]`).  The configuration a user
reads can therefore differ from the configuration whose loss is printed.

This script measures how often that matters: the recomputed SSR of the returned
configuration against the reported loss, for every evaluated run, with the
trigger condition (a returned position within 1e-4 of 0 or within 1e-3 of 1)
recorded alongside.  With `--rerun NAME` the case is re-run with its output
directory kept, so the ladder's own anchors can be compared with the final
record's anchors directly.
"""
import argparse
import json
import sys
import tempfile
from collections import Counter
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
sys.path.insert(0, str(ROOT / 'data' / 'stage289_decodiphy' / 'repo'))
import stage300_common as s300                   # noqa: E402
import stage300_evaluate as s300eval             # noqa: E402

TRIGGER_LOW, TRIGGER_HIGH = 1e-4, 1e-3


def classify(collection, arms=('dd_default', 'dd_true_k')):
    fits = {}
    for line in open(s300.OUT / ('stage300_fits_%s.jsonl' % collection)):
        row = json.loads(line)
        if row.get('evaluated') and row['arm'] in arms:
            fits.setdefault((row['name'], row['arm']), []).append(row)
    harvest = {}
    for row in s300eval.runs_for(collection):
        if row.get('status') == 'complete' and row.get('sources'):
            harvest.setdefault((row['name'], row['arm'], row['source_line'],
                                row['source_stage']), []).append(row)
    out = []
    for key, rows in sorted(fits.items()):
        for row in rows:
            run = harvest.get((row['name'], row['arm'], row['source_line'],
                               row['source_stage']))
            positions = []
            if run:
                positions = [float(s.get('x_decodiphy', s['x'])) for s in run[0]['sources']]
            reported = row.get('reported_loss')
            ssr = row['ssr_dd']
            rel = (None if reported is None
                   else abs(ssr - reported) / max(abs(reported), 1e-300))
            out.append(dict(collection=collection, name=row['name'], arm=row['arm'],
                            source_stage=row['source_stage'],
                            source_line=row['source_line'], k=row['k'],
                            reported_loss=reported, ssr_returned=ssr,
                            rel_diff=rel, stop_reason=row.get('stop_reason'),
                            ssr_worse_than_reported=(None if reported is None
                                                     else ssr > reported),
                            endpoint_trigger=bool(positions and
                                                  (min(positions) < TRIGGER_LOW or
                                                   max(positions) > 1 - TRIGGER_HIGH)),
                            min_x_dd=(min(positions) if positions else None),
                            max_x_dd=(max(positions) if positions else None)))
    return out


def rerun(collection, name, method='dd_default'):
    """Re-run one case keeping all_rounds.json, to exhibit the substitution."""
    import contextlib
    import io
    import decodiphy.cli as ddcli
    from stage289_decodiphy_bridge import write_inputs, to_bank_sources
    entry = {e['name']: e for e in s300.case_entries(collection)}[name]
    case = s300.load_case(entry)
    edges, leaves, dvals = s300.case_arrays(case)
    bank = dict(edges=edges, observed={r: v for r, v in zip(leaves, dvals)},
                generating_sources=case.get('generating_sources') or [])
    work = Path(tempfile.mkdtemp(prefix='stage300_rerun_'))
    mapping = write_inputs(bank, work)
    outdir = work / method
    outdir.mkdir(parents=True, exist_ok=True)
    argv = ['-t', str(work / 'tree.nwk'), '-d', str(work / 'distances.txt'),
            '-o', str(outdir)]
    if method == 'dd_true_k':
        argv += ['-k', str(len(bank['generating_sources'])), '-f', '1']
    saved, buf = sys.argv, io.StringIO()
    sys.argv = ['decodiphy'] + argv
    crash = None
    try:
        with contextlib.redirect_stdout(buf):
            ddcli.main()
    except BaseException as exc:
        crash = repr(exc)[:300]
    finally:
        sys.argv = saved
    rounds = json.loads((outdir / 'all_rounds.json').read_text())
    final = rounds[-1]
    ladder_last = [r for r in rounds[:-1]][-1] if len(rounds) > 1 else None
    report = dict(collection=collection, name=name, method=method, crash=crash,
                  n_rounds=len(rounds),
                  final=dict(k=final['k'], anchors=list(final['anchors']),
                             loss=final['loss'], x=list(final['x']),
                             p=list(final['p'])),
                  ladder_last=(None if ladder_last is None else
                               dict(k=ladder_last['k'], anchors=list(ladder_last['anchors']),
                                    loss=ladder_last['loss'], x=list(ladder_last['x']),
                                    p=list(ladder_last['p']))),
                  anchors_changed_after_hill_climbing=(
                      None if ladder_last is None else
                      sorted(final['anchors']) != sorted(ladder_last['anchors'])))
    for tag, entry_ in (('final', final), ('ladder_last', ladder_last)):
        if entry_ is None:
            continue
        unmapped = [a for a in entry_['anchors'] if a not in mapping['edge_of_label']]
        if unmapped:
            report[tag + '_ssr'] = dict(status='unmapped', unmapped=unmapped)
            continue
        srcs = to_bank_sources(bank, entry_['anchors'], entry_['p'], entry_['x'], mapping)
        y = entry_['y']
        pendant = F(float(sum(y)) if isinstance(y, list) else float(y))
        sources = [dict(edge=int(s['edge']), p=F(float(s['p'])), x=F(float(s['x'])))
                   for s in srcs]
        metrics = s300.fit_metrics(edges, leaves, dvals, sources, pendant)
        report[tag + '_ssr'] = dict(status='ok', ssr=metrics['ssr'],
                                    linf=str(metrics['linf']),
                                    linf_float=float(metrics['linf']))
    return report


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--collections', default='krepp_ani,krepp_placement')
    ap.add_argument('--rerun', default='')
    ap.add_argument('--out', default='stage300_lossgap.json')
    args = ap.parse_args()
    report = dict(stage=300, step='loss_gap', by_collection={}, reruns=[])
    rows_all = []
    for collection in args.collections.split(','):
        rows = classify(collection)
        rows_all.extend(rows)
        gap = [r for r in rows if r['rel_diff'] is not None and r['rel_diff'] > 1e-6]
        report['by_collection'][collection] = dict(
            n_runs=len(rows), n_gap_above_1e_6=len(gap),
            n_gap_and_endpoint_trigger=sum(1 for r in gap if r['endpoint_trigger']),
            n_gap_without_trigger=sum(1 for r in gap if not r['endpoint_trigger']),
            n_returned_worse_than_reported=sum(1 for r in rows
                                               if r['ssr_worse_than_reported']),
            stop_reasons_in_gap=dict(Counter(r['stop_reason'] for r in gap)),
            worst=sorted(gap, key=lambda r: -r['rel_diff'])[:5])
    (s300.OUT / 'stage300_lossgap_rows.jsonl').write_text(
        '\n'.join(json.dumps(s300.serial(r)) for r in rows_all) + '\n')
    if args.rerun:
        collection, name, method = args.rerun.split(':')
        report['reruns'].append(rerun(collection, name, method))
    (s300.OUT / args.out).write_text(json.dumps(s300.serial(report), indent=1))
    for collection, cell in report['by_collection'].items():
        print(collection, {k: v for k, v in cell.items() if k != 'worst'})
    for r in report['reruns']:
        print('rerun', r['name'], 'anchors_changed', r['anchors_changed_after_hill_climbing'],
              'final_ssr', r.get('final_ssr'), 'ladder_last_ssr', r.get('ladder_last_ssr'),
              'reported_final_loss', r['final']['loss'])
    print('wrote', s300.OUT / args.out)


if __name__ == '__main__':
    main()
