#!/usr/bin/env python3
"""Stage 291 runner: certified exclusion power under estimation error and
reference incompleteness. Writes only stage291_* files.

Tracks
  full     every bank case with a published generating mixture, full panel
  removal  the same cases under the Stage 288 removal regimes
  real     the 29 travellers' diarrhoea samples (no ground truth; certificates only)

Measurement. For a fixed panel and fixed supplied bounds (alpha, rho, floor),
an original edge j is certifiably excluded at error allowance eta if and only if
eta < tau_j, where
    tau_j = max over endpoints w of  alpha * offset_w(j) + floor - raw_w ,
raw_w being the endpoint half-sum returned by the frozen Stage 44 screen. tau_j
is therefore the largest error allowance at which the exclusion certificate for
edge j still holds, and the whole eta-curve follows from the tau vector.
"""
import os
for _k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_k] = '1'
import argparse, binascii, hashlib, json, math, statistics, sys, time
from datetime import datetime, timezone
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from stage291_common import (ROOT as CROOT, bank_index, load_any, model_vector, screen_panel,
                             leaves_of, serial, sha, run_pool, clade_removals,
                             nearest_source_removal, random_removal)

PROTOCOL_DEFAULT = 'results/stage291_protocol_v2.json'
SEED = 2910007
NEG = float('-inf')

ETA_GRID = ([F(0), F(3, 10), F(1)] + [F(1, 10 ** k) for k in (16, 14, 12, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1)]
            + [F(3, 10 ** k) for k in (7, 6, 5, 4, 3, 2)])
ETA_GRID = sorted(set(ETA_GRID))
LEVELS = (F(1, 10), F(1, 4), F(1, 2), F(3, 4))
REGIMES = ('random', 'clade', 'nearest')
REPLICATES = {'random': 2, 'clade': 2, 'nearest': 1}
WALL, ADDRESS_SPACE = 3600.0, 4 * 1024 ** 3


def hypotheses(case):
    """Supplied H2/H3 bounds. 'oracle' uses the published truth and is the most
    favourable admissible choice; the others are stated assumptions."""
    cap = case['cap'] or 10
    out = {}
    src = case['sources']
    y = case['pendant']
    if src and all(s['p'] is not None for s in src) and y is not None:
        out['oracle'] = dict(alpha=min(s['p'] for s in src),
                             rho=min(min(s['x'], 1 - s['x']) for s in src), floor=y)
    if y is not None:
        out['assumed_mid'] = dict(alpha=F(1, 2 * cap), rho=F(1, 10), floor=y / 2)
        out['assumed_weak'] = dict(alpha=F(1, 4 * cap), rho=F(1, 100), floor=y / 10)
    out['data_scaled'] = dict(alpha=F(1, 2 * cap), rho=F(1, 10),
                              floor=min(case['observed'].values()) / 2)
    return out


def tau_vector(pan, alpha, rho, floor):
    """tau_j for every original edge that the panel can witness at all.

    An edge whose two hull endpoints both lack an eligible avoiding leaf pair
    carries no certificate at any eta and is simply absent from the result."""
    exact = {}
    for (raw_u, raw_v), items in zip(pan['raw'], pan['offsets']):
        for j, pre, suf, ell in items:
            best = None
            if raw_u is not None:
                best = alpha * (pre + rho * ell) + floor - raw_u
            if raw_v is not None:
                cand = alpha * (suf + rho * ell) + floor - raw_v
                best = cand if best is None else max(best, cand)
            if best is not None:
                exact[j] = best
    return exact


def quantiles(values):
    if not values:
        return {}
    s = sorted(values)
    q = lambda f: s[min(len(s) - 1, int(f * (len(s) - 1) + 0.5))]
    return dict(min=s[0], q05=q(0.05), q10=q(0.10), q25=q(0.25), median=q(0.5),
                q75=q(0.75), q90=q(0.90), q95=q(0.95), max=s[-1])


def panel_record(case, keep, truth_vector):
    n_edges = len(case['edges'])
    pan = screen_panel(case, keep)
    covered = set(pan['covered'])
    eta_plus = eta_abs = None
    if truth_vector is not None:
        diffs = [truth_vector[r] - case['observed'][r] for r in keep]
        eta_plus = max([F(0)] + diffs)
        eta_abs = max(abs(d) for d in diffs)
    true_edges = sorted({s['edge'] for s in case['sources']}) if case['sources'] else []
    out = dict(kept=len(keep), n_edges=n_edges, n_hull=pan['n_hull'], n_covered=len(covered),
               true_edges=true_edges, true_covered=[j for j in true_edges if j in covered],
               eta_plus=(None if eta_plus is None else str(eta_plus)),
               eta_abs=(None if eta_abs is None else str(eta_abs)),
               eta_plus_f=(None if eta_plus is None else float(eta_plus)),
               eta_abs_f=(None if eta_abs is None else float(eta_abs)), hypotheses={})
    for name, h in hypotheses(case).items():
        exact = tau_vector(pan, h['alpha'], h['rho'], h['floor'])
        witnessed = sorted(exact)
        grid = [sum(1 for v in exact.values() if v > e) for e in ETA_GRID]
        tau_true = {str(j): (float(exact[j]) if j in exact else None) for j in true_edges}
        present = [exact[j] for j in true_edges if j in exact]
        rec = dict(alpha=str(h['alpha']), rho=str(h['rho']), floor=str(h['floor']),
                   alpha_f=float(h['alpha']), rho_f=float(h['rho']), floor_f=float(h['floor']),
                   n_witnessed=len(witnessed), tau=quantiles([float(exact[j]) for j in witnessed]),
                   grid_excluded=grid, tau_true=tau_true,
                   tau_true_max=(float(max(present)) if present else None))
        if eta_plus is not None:
            rec['excluded_at_eta_plus'] = sum(1 for v in exact.values() if v > eta_plus)
            rec['excluded_at_eta_abs'] = sum(1 for v in exact.values() if v > eta_abs)
            rec['true_excluded_at_eta_plus'] = sum(1 for j in true_edges
                                                   if j in exact and exact[j] > eta_plus)
            rec['floor_valid'] = bool(case['pendant'] is not None and h['floor'] <= case['pendant'])
        out['hypotheses'][name] = rec
    return out


def panels_for(case, track):
    keep_all = sorted(case['observed'])
    yield dict(regime='none', level='0', replicate=0), keep_all
    if track != 'removal':
        return
    leaves = keep_all
    total = len(leaves)
    for level in LEVELS:
        count = int(level * total)
        if count < 1 or total - count < 3:
            continue
        for rep in range(REPLICATES['random']):
            tag = binascii.crc32(case['name'].encode()) % 9973
            drop = set(random_removal(leaves, count, SEED + 101 * rep + tag))
            yield dict(regime='random', level=str(level), replicate=rep), sorted(set(leaves) - drop)
        for rep, item in enumerate(clade_removals(case['edges'], leaves, count, REPLICATES['clade'])):
            yield (dict(regime='clade', level=str(level), replicate=rep,
                        clade_vertex=item['clade_vertex'], complement=item['complement']),
                   sorted(set(leaves) - set(item['removed'])))
        if case['sources']:
            drop = set(nearest_source_removal(case['edges'], case['sources'], leaves, count))
            yield dict(regime='nearest', level=str(level), replicate=0), sorted(set(leaves) - drop)


def worker(spec):
    case = load_any(spec['entry'])
    truth = model_vector(case)
    rows = []
    for meta, keep in panels_for(case, spec['track']):
        if len(keep) < 2:
            rows.append(dict(meta, status='skipped_panel_too_small'))
            continue
        rec = panel_record(case, keep, truth)
        rows.append(dict(meta, status='ok', **rec))
    return dict(status='ok', family=case['family'], name=case['name'],
                refs=len(case['observed']), K=case['K_generating'],
                cap=case['cap'], tree_leaves=case['tree_leaves'],
                noise_tag=case['provenance'].get('noise_tag'),
                seqlen=case['provenance'].get('sequence_length'),
                median_observed=float(statistics.median(sorted(case['observed'].values()))),
                median_edge=float(statistics.median(sorted(e[2] for e in case['edges']))),
                min_edge=float(min(e[2] for e in case['edges'])),
                truth_available=truth is not None, panels=rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--track', required=True, choices=('full', 'removal', 'real'))
    parser.add_argument('--families', default='')
    parser.add_argument('--workers', type=int, default=28)
    parser.add_argument('--limit', type=int, default=0)
    parser.add_argument('--tag', default='')
    parser.add_argument('--protocol', default=PROTOCOL_DEFAULT)
    args = parser.parse_args()
    protocol = json.loads((ROOT / args.protocol).read_text())
    for name, expected in protocol['frozen_inputs'].items():
        assert sha(ROOT / name) == expected, f'input changed since protocol: {name}'
    families = ([f for f in args.families.split(',') if f] or
                protocol['tracks'][args.track]['families'])
    entries = [e for e in bank_index() if e['family'] in families]
    entries.sort(key=lambda e: (e['family'], e['name']))
    if args.limit:
        entries = entries[:args.limit]
    stream = ROOT / f'results/stage291_{args.track}{args.tag}_runs.jsonl'
    assert not stream.exists(), f'refusing to append to existing stream {stream}'
    jobs = [dict(meta=dict(family=e['family'], name=e['name'], sha256=e['sha256'],
                           track=args.track), spec=dict(entry=e, track=args.track))
            for e in entries]
    print(f'track={args.track} cases={len(jobs)} workers={args.workers}', flush=True)
    began = time.perf_counter()
    run_pool(jobs, worker, WALL, ADDRESS_SPACE, args.workers, stream,
             progress=max(1, len(jobs) // 40))
    print(f'done in {time.perf_counter() - began:.1f}s -> {stream}', flush=True)


if __name__ == '__main__':
    main()
