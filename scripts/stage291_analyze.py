#!/usr/bin/env python3
"""Stage 291 analysis. Reads only stage291_* run streams; writes stage291_analysis.json."""
import json, math, statistics, sys, collections
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from stage291_common import sha

PROTOCOL = json.loads((ROOT / 'results/stage291_protocol_v2.json').read_text())
ETA = [F(s) for s in PROTOCOL['eta_grid']]
ETAF = [float(e) for e in ETA]
OUT = ROOT / 'results/stage291_analysis.json'


def med(xs):
    return float(statistics.median(xs)) if xs else None


def quant(xs, f):
    if not xs:
        return None
    s = sorted(xs)
    return float(s[min(len(s) - 1, int(f * (len(s) - 1) + 0.5))])


def bounds_valid(h, oracle):
    """Are the supplied H2/H3 bounds admissible for this instance?

    alpha must not exceed the smallest true source weight, rho the smallest true
    interior margin, and the pendant floor the true common mean pendant. The
    oracle record carries all three true values, so validity is decidable here.
    """
    if oracle is None:
        return None
    return bool(h['alpha_f'] <= oracle['alpha_f'] and h['rho_f'] <= oracle['rho_f']
                and h['floor_f'] <= oracle['floor_f'])


def stream(path):
    with open(path) as fh:
        for line in fh:
            yield json.loads(line)


def collapse_eta(grid, n, frac):
    """Smallest grid eta whose excluded fraction is below `frac` (None if never)."""
    for e, k in zip(ETAF, grid):
        if n and k / n < frac:
            return e
    return None


# ----------------------------------------------------------------- step 1
def step1():
    per_family = collections.defaultdict(lambda: collections.defaultdict(list))
    biotrees = collections.defaultdict(lambda: collections.defaultdict(list))
    paired = collections.defaultdict(list)
    base = {}
    violations, consistency = [], collections.defaultdict(list)
    curves = collections.defaultdict(lambda: [[] for _ in ETA])
    ncases = collections.Counter()
    for row in stream(ROOT / 'results/stage291_full_runs.jsonl'):
        res = row['result']
        if res.get('status') != 'ok':
            violations.append(dict(name=row['name'], reason='non_ok_case', status=res.get('status')))
            continue
        fam = res['family']
        ncases[fam] += 1
        panel = res['panels'][0]
        n = panel['n_edges']
        acc = per_family[fam]
        if panel['eta_plus_f'] is not None:
            acc['eta_plus'].append(panel['eta_plus_f'])
            acc['eta_abs'].append(panel['eta_abs_f'])
            acc['eta_plus_rel'].append(panel['eta_plus_f'] / res['median_observed'])
        oracle = panel['hypotheses'].get('oracle')
        for hname, h in panel['hypotheses'].items():
            key = f'{hname}'
            valid = bounds_valid(h, oracle)
            acc[f'bounds_valid_{key}'].append(1 if valid else 0)
            if 'excluded_at_eta_plus' in h:
                acc[f'power_sound_{key}'].append(h['excluded_at_eta_plus'] / n)
                acc[f'true_excl_sound_{key}'].append(h['true_excluded_at_eta_plus'])
                acc[f'any_true_excl_{key}'].append(1 if h['true_excluded_at_eta_plus'] else 0)
                if valid:
                    acc[f'power_sound_valid_{key}'].append(h['excluded_at_eta_plus'] / n)
                if h['true_excluded_at_eta_plus'] > 0:
                    violations.append(dict(name=row['name'], family=fam, hypothesis=hname,
                                           reason='true_source_edge_excluded_at_sound_eta',
                                           bounds_valid=valid,
                                           count=h['true_excluded_at_eta_plus']))
            if h['tau_true_max'] is not None and panel['eta_plus_f'] is not None:
                acc[f'tau_true_max_{key}'].append(h['tau_true_max'])
                if panel['eta_plus_f'] > 0:
                    acc[f'margin_ratio_{key}'].append(h['tau_true_max'] / panel['eta_plus_f'])
            if h['n_witnessed'] == n and h['tau']:
                consistency[(fam, hname)].append(
                    (h['tau']['min'], panel['eta_plus_f']))
            for i, k in enumerate(h['grid_excluded']):
                curves[(fam, hname)][i].append(k / n)
            acc[f'collapse01_{key}'].append(collapse_eta(h['grid_excluded'], n, 0.01))
            acc[f'collapse10_{key}'].append(collapse_eta(h['grid_excluded'], n, 0.10))
            acc[f'collapse50_{key}'].append(collapse_eta(h['grid_excluded'], n, 0.50))
        if fam == 'published_biotrees':
            tag = res['noise_tag'] or 'noise0'
            key = (tag, res['seqlen'])
            biotrees[key]['eta_plus'].append(panel['eta_plus_f'])
            o = panel['hypotheses'].get('oracle')
            if o:
                biotrees[key]['power_sound'].append(o['excluded_at_eta_plus'] / n)
            stem = row['name'].split('_noise')[0]
            if tag == 'noise0':
                base[stem] = panel['eta_plus_f']
            else:
                paired[key].append((stem, panel['eta_plus_f']))
    summary = {}
    for fam, acc in per_family.items():
        s = {'cases': ncases[fam]}
        for k, v in acc.items():
            vv = [x for x in v if x is not None]
            if not vv:
                s[k] = {'n': 0, 'all_none': True}
                continue
            s[k] = dict(n=len(vv), median=med(vv), q25=quant(vv, .25), q75=quant(vv, .75),
                        min=float(min(vv)), max=float(max(vv)),
                        n_none=len(v) - len(vv))
        summary[fam] = s
    bio = {}
    for key, acc in biotrees.items():
        inc = []
        for stem, val in paired.get(key, []):
            if stem in base:
                inc.append(val - base[stem])
        bio[f'{key[0]}|{key[1]}'] = dict(
            cases=len(acc['eta_plus']), eta_plus_median=med(acc['eta_plus']),
            eta_plus_q25=quant(acc['eta_plus'], .25), eta_plus_q75=quant(acc['eta_plus'], .75),
            power_sound_oracle_median=med(acc.get('power_sound', [])),
            paired_increment_over_noise0_median=med(inc), paired_n=len(inc))
    cons = {}
    for (fam, hname), pairs in consistency.items():
        lows = [a for a, b in pairs]
        bad = [1 for a, b in pairs if b is not None and a > b + 1e-18]
        cons[f'{fam}|{hname}'] = dict(n=len(pairs), certified_lower_bound_median=med(lows),
                                      certified_lower_bound_q90=quant(lows, .90),
                                      exceeds_true_eta_plus=len(bad))
    curve_out = {f'{a}|{b}': dict(eta=[str(e) for e in ETA],
                                  median_excluded_fraction=[med(v) for v in rows])
                 for (a, b), rows in curves.items()}
    return dict(summary=summary, biotrees_by_noise=bio, consistency=cons,
                curves=curve_out, violations=violations[:200],
                n_violations=len(violations),
                n_violations_valid_bounds=sum(1 for v in violations
                                              if v.get('bounds_valid') is True),
                violations_by_hypothesis=dict(collections.Counter(
                    v.get('hypothesis') for v in violations)))


# ----------------------------------------------------------------- step 2
def step2():
    cells = collections.defaultdict(lambda: collections.defaultdict(list))
    grid_cells = collections.defaultdict(lambda: [[] for _ in ETA])
    viol = []
    for row in stream(ROOT / 'results/stage291_removal_runs.jsonl'):
        res = row['result']
        if res.get('status') != 'ok':
            viol.append(dict(name=row['name'], reason='non_ok_case', status=res.get('status')))
            continue
        fam = res['family']
        for panel in res['panels']:
            if panel.get('status') != 'ok':
                continue
            n = panel['n_edges']
            level = panel['level'] if panel['regime'] != 'none' else '0'
            key = (fam, panel['regime'], level)
            cells[key]['coverage'].append(panel['n_covered'] / n)
            cells[key]['kept_frac'].append(panel['kept'] / res['refs'])
            cells[key]['true_covered_frac'].append(
                len(panel['true_covered']) / max(1, len(panel['true_edges'])))
            if panel['eta_plus_f'] is not None:
                cells[key]['eta_plus'].append(panel['eta_plus_f'])
            oracle = panel['hypotheses'].get('oracle')
            for hname, h in panel['hypotheses'].items():
                valid = bounds_valid(h, oracle)
                if 'excluded_at_eta_plus' in h:
                    cells[key][f'power_sound_{hname}'].append(h['excluded_at_eta_plus'] / n)
                    cells[key][f'true_excl_{hname}'].append(h['true_excluded_at_eta_plus'])
                    cells[key][f'bounds_valid_{hname}'].append(1 if valid else 0)
                    if h['true_excluded_at_eta_plus'] > 0:
                        viol.append(dict(name=row['name'], regime=panel['regime'],
                                         level=level, hypothesis=hname, bounds_valid=valid,
                                         reason='true_source_edge_excluded_at_sound_eta'))
                if hname == 'oracle':
                    for i, k in enumerate(h['grid_excluded']):
                        grid_cells[key][i].append(k / n)
    out = {}
    for key, acc in cells.items():
        s = {}
        for k, v in acc.items():
            vv = [x for x in v if x is not None]
            s[k] = dict(n=len(vv), median=med(vv), q25=quant(vv, .25), q75=quant(vv, .75)) if vv else None
        out['|'.join(map(str, key))] = s
    grids = {'|'.join(map(str, key)): dict(eta=[str(e) for e in ETA],
                                           median_excluded_fraction=[med(v) for v in rows],
                                           n=[len(v) for v in rows])
             for key, rows in grid_cells.items()}
    return dict(cells=out, oracle_grids=grids, violations=viol[:200], n_violations=len(viol),
                n_violations_valid_bounds=sum(1 for v in viol if v.get('bounds_valid') is True),
                violations_by_hypothesis=dict(collections.Counter(
                    v.get('hypothesis') for v in viol)))


# ----------------------------------------------------------------- real data
def real():
    rows = []
    for row in stream(ROOT / 'results/stage291_real_runs.jsonl'):
        res = row['result']
        if res.get('status') != 'ok':
            rows.append(dict(name=row['name'], status=res.get('status')))
            continue
        panel = res['panels'][0]
        h = panel['hypotheses']['data_scaled']
        n = panel['n_edges']
        rows.append(dict(name=res['name'], refs=res['refs'], n_edges=n,
                         n_witnessed=h['n_witnessed'],
                         median_observed=res['median_observed'], median_edge=res['median_edge'],
                         alpha=h['alpha'], rho=h['rho'], floor=h['floor_f'],
                         tau=h['tau'],
                         excluded_fraction={str(e): k / n for e, k in zip(ETA, h['grid_excluded'])},
                         certified_error_lower_bound=(h['tau']['min']
                                                      if h['n_witnessed'] == n else None)))
    ok = [r for r in rows if 'tau' in r]
    agg = dict(samples=len(rows), complete_witness=sum(1 for r in ok if r['certified_error_lower_bound'] is not None),
               certified_error_lower_bound_median=med([r['certified_error_lower_bound'] for r in ok
                                                       if r['certified_error_lower_bound'] is not None]),
               excluded_fraction_median_by_eta={str(e): med([r['excluded_fraction'][str(e)] for r in ok])
                                                for e in ETA})
    return dict(samples=rows, aggregate=agg)


def main():
    result = dict(stage=291, protocol_sha256=sha(ROOT / 'results/stage291_protocol_v2.json'),
                  inputs={name: sha(ROOT / 'results' / name) for name in
                          ('stage291_full_runs.jsonl', 'stage291_removal_runs.jsonl',
                           'stage291_real_runs.jsonl')},
                  step1=step1(), step2=step2(), real=real())

    OUT.write_text(json.dumps(result, indent=2) + '\n')
    print('wrote', OUT)
    for part in ('step1', 'step2'):
        print(part, 'violations', result[part]['n_violations'],
              'of which with valid supplied bounds', result[part]['n_violations_valid_bounds'],
              result[part]['violations_by_hypothesis'])


if __name__ == '__main__':
    main()
