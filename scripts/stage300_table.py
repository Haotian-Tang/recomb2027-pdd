#!/usr/bin/env python3
"""Stage 300: per-case three-way table as CSV (one row per case and arm).

The pooled counts live in stage300_analysis.json; this is the per-case form the
paper's supplement needs, with the exact rationals alongside their float values.
"""
import csv
import json
import sys
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import stage300_common as s300                   # noqa: E402

FIELDS = ['collection', 'name', 'arm', 'n_leaves', 'n_edges', 'K_true', 'k_returned',
          'theta1', 'theta1_float', 'r_dd', 'r_dd_float', 'r_dd_over_theta1',
          'r_dd_ge_theta1', 'truth_linf', 'truth_linf_float', 'truth_over_theta1',
          'truth_ge_theta1', 'ordering_linf', 'lambda1_ssr', 'ssr_dd', 'truth_ssr',
          'ssr_over_lambda1', 'dd_beats_one_source_in_ssr', 'dd_beats_truth_in_ssr',
          'reported_loss', 'loss_rel_diff', 'stop_reason', 'n_kept_by_tool',
          'theta1_restricted', 'r_dd_kept', 'r_dd_kept_over_theta1_kept',
          'resolved_status', 'resolved_ssr', 'resolved_linf', 'resolved_linf_over_theta1',
          'source_stage', 'source_file', 'source_line']


def main():
    analysis = json.loads((s300.OUT / 'stage300_analysis.json').read_text())
    rows = []
    for collection in ('krepp_ani', 'krepp_placement', 'travelers'):
        cell = analysis['by_collection'].get(collection)
        if not cell:
            continue
        bases = {}
        base_path = s300.OUT / ('stage300_baselines_%s.jsonl' % collection)
        if base_path.exists():
            for line in open(base_path):
                row = json.loads(line)
                bases[row['name']] = row
        prune = {}
        prune_path = s300.OUT / ('stage300_prune_%s.jsonl' % collection)
        if prune_path.exists():
            for line in open(prune_path):
                row = json.loads(line)
                prune[row['name']] = row
        resolved = {}
        res_path = s300.OUT / ('stage300_resolved_%s.jsonl' % collection)
        if res_path.exists():
            for line in open(res_path):
                row = json.loads(line)
                resolved[(row['name'], row['arm'])] = row
        for arm, cases in cell['per_case'].items():
            if arm not in ('dd_default', 'dd_true_k'):
                continue
            for case in cases:
                base = bases.get(case['name'], {})
                pr = prune.get(case['name'], {})
                pra = (pr.get('arms') or {}).get(arm, {})
                rs = resolved.get((case['name'], arm), {})
                rsv = rs.get('resolved', {})
                truth = base.get('truth') or {}
                rows.append(dict(
                    collection=collection, name=case['name'], arm=arm,
                    n_leaves=case['n_leaves'], n_edges=base.get('n_edges'),
                    K_true=truth.get('k_true'), k_returned=case['k'],
                    theta1=base.get('theta1'), theta1_float=case.get('theta1_float'),
                    r_dd=case['r_dd'], r_dd_float=case['r_dd_float'],
                    r_dd_over_theta1=case.get('r_dd_over_theta1'),
                    r_dd_ge_theta1=case.get('r_dd_ge_theta1'),
                    truth_linf=case.get('truth_linf'),
                    truth_linf_float=(float(F(case['truth_linf']))
                                      if case.get('truth_linf') else None),
                    truth_over_theta1=case.get('truth_over_theta1'),
                    truth_ge_theta1=case.get('truth_ge_theta1'),
                    ordering_linf=case.get('ordering_linf'),
                    lambda1_ssr=case.get('lambda1_ssr'), ssr_dd=case.get('ssr_dd'),
                    truth_ssr=case.get('truth_ssr'),
                    ssr_over_lambda1=case.get('ssr_over_lambda1'),
                    dd_beats_one_source_in_ssr=case.get('dd_beats_one_source_in_ssr'),
                    dd_beats_truth_in_ssr=case.get('dd_beats_truth_in_ssr'),
                    reported_loss=case.get('reported_loss'),
                    loss_rel_diff=case.get('loss_rel_diff'),
                    stop_reason=case.get('stop_reason'),
                    n_kept_by_tool=pr.get('n_kept'),
                    theta1_restricted=pr.get('theta1_restricted'),
                    r_dd_kept=pra.get('r_dd_kept'),
                    r_dd_kept_over_theta1_kept=pra.get('r_dd_kept_over_theta1_kept'),
                    resolved_status=rsv.get('status'), resolved_ssr=rsv.get('ssr'),
                    resolved_linf=rsv.get('linf'),
                    resolved_linf_over_theta1=rs.get('linf_over_theta1'),
                    source_stage=case['source_stage'], source_file=case['source_file'],
                    source_line=None))
    path = s300.OUT / 'stage300_three_way.csv'
    with open(path, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, extrasaction='ignore')
        writer.writeheader()
        for row in sorted(rows, key=lambda r: (r['collection'], r['name'], r['arm'])):
            writer.writerow(row)
    print('wrote %s with %d rows' % (path, len(rows)))


if __name__ == '__main__':
    main()
