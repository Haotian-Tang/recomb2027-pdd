#!/usr/bin/env python3
"""Stage 291 joint table: certified exclusion power by source count, estimation
error allowance and reference completeness. Reads the stage 291 run streams only."""
import json, collections, statistics, sys
from fractions import Fraction as F
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from stage291_common import sha

ETA = json.loads((ROOT / 'results/stage291_protocol_v2.json').read_text())['eta_grid']
ETAF = [float(F(e)) for e in ETA]
OUT = ROOT / 'results/stage291_joint_table.json'
med = lambda a: float(statistics.median(a)) if a else None


def curve_by_K():
    acc = collections.defaultdict(lambda: [[] for _ in ETA])
    for line in open(ROOT / 'results/stage291_full_runs.jsonl'):
        d = json.loads(line); res = d['result']
        if res.get('status') != 'ok':
            continue
        p = res['panels'][0]; o = p['hypotheses'].get('oracle')
        if o is None:
            continue
        n = p['n_edges']
        for i, k in enumerate(o['grid_excluded']):
            acc[(res['family'], res['K'])][i].append(k / n)
    return {f'{a}|K{b}': dict(n=len(v[0]), median=[med(x) for x in v],
                              mean=[(sum(x) / len(x) if x else None) for x in v],
                              frac_any=[(sum(1 for y in x if y > 0) / len(x) if x else None) for x in v])
            for (a, b), v in acc.items()}


def joint_table():
    cells = collections.defaultdict(lambda: [[] for _ in ETA])
    cover = collections.defaultdict(list)
    sound = collections.defaultdict(list)
    for line in open(ROOT / 'results/stage291_removal_runs.jsonl'):
        d = json.loads(line); res = d['result']
        if res.get('status') != 'ok':
            continue
        fam, K = res['family'], res['K']
        for panel in res['panels']:
            if panel.get('status') != 'ok':
                continue
            o = panel['hypotheses'].get('oracle')
            if o is None:
                continue
            n = panel['n_edges']
            level = '0' if panel['regime'] == 'none' else panel['level']
            key = (fam, K, panel['regime'], level)
            for i, k in enumerate(o['grid_excluded']):
                cells[key][i].append(k / n)
            cover[key].append(panel['n_covered'] / n)
            if 'excluded_at_eta_plus' in o:
                sound[key].append(o['excluded_at_eta_plus'] / n)
    out = {}
    for key, v in cells.items():
        out['|'.join(map(str, key))] = dict(
            n=len(v[0]), coverage_median=med(cover[key]),
            power_at_eta_plus_median=med(sound[key]),
            frac_any_at_eta_plus=(sum(1 for x in sound[key] if x > 0) / len(sound[key])
                                  if sound[key] else None),
            median=[med(x) for x in v],
            frac_any=[(sum(1 for y in x if y > 0) / len(x) if x else None) for x in v])
    return out


def main():
    result = dict(stage=291, eta_grid=ETA, eta_float=ETAF,
                  protocol_sha256=sha(ROOT / 'results/stage291_protocol_v2.json'),
                  curve_by_family_K=curve_by_K(), joint=joint_table())
    OUT.write_text(json.dumps(result, indent=2) + '\n')
    print('wrote', OUT)


if __name__ == '__main__':
    main()
