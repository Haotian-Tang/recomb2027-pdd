#!/usr/bin/env python3
"""Stage 312: what a tolerance diagnostic certifies on an existing candidate.

The practitioner criterion in the discussion is positive: for a verified strict
K-source fit with worst-leaf residual r_fit on the same tree and panel as
Theta_1, every budget r_fit < eps < Theta_1 has 2 <= m(eps) <= K.  The
measurements reported elsewhere in the paper are mostly about what fails, so
this script asks the opposite question on the published ANI benchmark: how
often does an already-available candidate supply such an interval, and what
does one worked case look like?

A heuristic's float output is not a strict configuration: the weights need not
sum to exactly one, a position may sit at an endpoint and a pendant may be
marginally negative.  Each candidate is therefore mapped to a strict witness by
one uniform rule, with delta = 1e-12 and s_d the case's registered full-panel
median observation,

    p_i' = p_i / sum_j p_j,
    x_i' = (1 - delta) x_i + delta / 2,
    y'   = max(0, y) + delta * s_d,

and then rescored in exact rational arithmetic.  Nothing is re-fitted and no
candidate is reselected: the canonical run per (case, arm) is stage 300's.

Reads the registered cases, Stage 298/300 baselines, and the bulk files
stage300_harvest.jsonl and stage300_fits_krepp_ani.jsonl. Fetch the README's bulk
archive first, or pass --bulk-root PATH to use an existing copy. All strict
witnesses and their candidate identities are retained; the example additionally
retains every residual and a three-reference lower certificate for Theta_1.
Writes results/stage312/stage312_positive_intervals.json; --output redirects it.
"""
import argparse
import json
import os
import statistics
import sys
from fractions import Fraction as F
from pathlib import Path

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import stage300_common as s300                          # noqa: E402
import stage300_analyze as s300an                        # noqa: E402

COLLECTION = 'krepp_ani'
ARM = 'dd_default'
DELTA = F(1, 10 ** 12)
EXAMPLE = 'k10k_tree_9_k5_s1_uniform_krepp_estimated'
EXAMPLE_EPS = F(3, 100)
ROOT = Path(HERE).parent
HARVEST = 'results/stage300/stage300_harvest.jsonl'
FITS = 'results/stage300/stage300_fits_krepp_ani.jsonl'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def identity(run):
    return tuple(run[k] for k in ('name', 'arm', 'source_stage',
                                  'source_file', 'source_line'))


def registered_candidates(bulk_root):
    """Join Stage 300's evaluated canonical choices to their saved parameters."""
    manifest = {r['path']: r for r in
                json.loads((ROOT / 'bulk_manifest.json').read_text())['members']}
    inputs = {}
    for relative in (HARVEST, FITS):
        path = bulk_root / relative
        require(path.is_file(), 'missing bulk input %s; fetch the README bulk '
                'archive or pass --bulk-root' % path)
        require(path.stat().st_size == manifest[relative]['bytes']
                and s300.sha(path) == manifest[relative]['sha256'],
                'bulk hash/size mismatch: ' + relative)
        inputs[relative] = s300.sha(path)
    with (bulk_root / HARVEST).open() as handle:
        raw = [r for r in map(json.loads, handle) if r['collection'] == COLLECTION]
    # The two bulk inputs may live separately from the public checkout.
    # Normalize public Stage 300 runs as stage300_evaluate.runs_for does;
    # selection itself remains stage300_analyze.canonical on frozen fits.
    own_runs = sorted((ROOT / 'results/stage300').glob('stage300_runs_*.jsonl'))
    for public in own_runs:
        relative = public.relative_to(ROOT).as_posix()
        inputs[relative] = s300.sha(public)
        with public.open() as handle:
            for line_no, line in enumerate(handle, 1):
                row = json.loads(line)
                if row.get('collection') != COLLECTION:
                    continue
                result = row.get('result', {})
                raw.append(dict(name=row['name'], arm=row['method'],
                                source_stage=300, source_file=relative,
                                source_line=line_no, status=result.get('status'),
                                sources=result.get('sources'),
                                pendant_total=result.get('pendant_total')))
    candidates = {}
    for run in raw:
        if run.get('arm') == ARM:
            key = identity(run)
            require(key not in candidates, 'duplicate raw candidate: %r' % (key,))
            candidates[key] = run
    groups = {}
    with (bulk_root / FITS).open() as handle:
        for line in handle:
            row = json.loads(line)
            if row.get('collection') == COLLECTION and row.get('arm') == ARM:
                groups.setdefault(row['name'], []).append(row)
    selected = {}
    for name, fits in groups.items():
        chosen = s300an.canonical(fits)
        require(chosen is not None, 'no evaluated canonical candidate: ' + name)
        run = candidates.get(identity(chosen))
        require(run is not None and run.get('status') == 'complete'
                and run.get('sources'), 'canonical parameters unavailable: ' + name)
        selected[name] = (run, chosen)
    return selected, inputs


def strictify(run, scale, n_edges):
    """Apply the uniform map only when the original candidate meets its domain.

    Weights must be positive, original positions must lie in [0, 1], edges
    must be distinct and valid, and the registered scale must be positive.
    An invalid candidate receives a reason; out-of-range positions are never
    repaired merely because contraction would bring them inside the edge.
    """
    try:
        src = run['sources']
        weights = [F(s['p']) for s in src]
        positions = [F(s['x']) for s in src]
        indices = [F(s['edge']) for s in src]
        pendant = F(run['pendant_total'])
        scale = F(scale)
    except (ValueError, TypeError, KeyError, ZeroDivisionError, OverflowError):
        return None, 'invalid_numeric_input'
    if not src:
        return None, 'empty_sources'
    if scale <= 0:
        return None, 'nonpositive_scale'
    if any(p <= 0 for p in weights):
        return None, 'nonpositive_weight'
    if any(not 0 <= position <= 1 for position in positions):
        return None, 'position_out_of_range'
    if any(e.denominator != 1 or not 0 <= e < n_edges for e in indices):
        return None, 'invalid_edge'
    tot = sum(weights)
    sources, seen = [], set()
    for edge, p, position in zip(indices, weights, positions):
        e = int(edge)
        if e in seen:                      # distinct original edges only
            return None, 'repeated_edge'
        seen.add(e)
        x = (1 - DELTA) * position + DELTA / 2
        sources.append(dict(edge=e, p=p / tot, x=x))
    y = max(F(0), pendant) + DELTA * scale
    require(all(s['p'] > 0 and 0 < s['x'] < 1 for s in sources)
            and sum(s['p'] for s in sources) == 1 and y > 0,
            'strictification postcondition failed')
    return (sources, y), 'ok'


def example_lower_certificate(edges, leaves, dvals):
    """Four-point necessity: report perturbations close this gap by at most 2 eps."""
    refs = [0, 24, 36]
    observed = dict(zip(leaves, dvals))
    maps = s300._distance_maps(edges, refs)
    sums = [observed[r] + maps[refs[(i + 1) % 3]][refs[(i + 2) % 3]]
            for i, r in enumerate(refs)]
    ordered = sorted(sums)
    bound = (ordered[-1] - ordered[-2]) / 2
    return dict(references=refs, observed_plus_opposite_distance=list(map(str, sums)),
                lower_bound=str(bound),
                reason='For one source plus a pendant, the largest two four-point '
                       'sums are equal. Changing each report by at most eps can '
                       'close the observed gap by at most 2*eps.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bulk-root', type=Path, default=ROOT)
    parser.add_argument('--output', type=Path,
                        default=ROOT / 'results/stage312/stage312_positive_intervals.json')
    args = parser.parse_args()
    theta = s300.theta1_table()
    entries = {e['name']: e for e in s300.case_entries(COLLECTION)}
    selected, inputs = registered_candidates(args.bulk_root)
    require(len(entries) == 296 and set(selected) == set(entries),
            'expected canonical candidates for all 296 registered ANI cases')
    baseline_path = ROOT / 'results/stage300/stage300_baselines_krepp_ani.jsonl'
    with baseline_path.open() as handle:
        baselines = {r['name']: r for r in map(json.loads, handle)}
    for relative in ('data/stage287_instances/bank/index.json',
                     'results/stage298/stage298_analysis.json',
                     'results/stage298/stage298_td_summary.json',
                     baseline_path.relative_to(ROOT).as_posix()):
        inputs[relative] = s300.sha(ROOT / relative)

    rows, example = [], None
    for name in sorted(entries):
        run, selected_fit = selected[name]
        t1 = theta.get(name)
        require(t1 is not None and t1['theta1'] > 0, 'missing/nonpositive Theta_1: ' + name)
        entry, base = entries[name], baselines[name]
        require(s300.sha(ROOT / entry['path']) == entry['sha256'], 'case hash mismatch: ' + name)
        case = s300.load_case(entries[name])
        edges, leaves, dvals = s300.case_arrays(case)
        require(base['path'] == entry['path'] and base['sha256'] == entry['sha256']
                and base['n_leaves'] == len(leaves) and base['n_edges'] == len(edges)
                and F(base['theta1']) == t1['theta1'] == F(selected_fit['theta1'])
                and t1['own_median_observed'] == s300.median_observed(dvals)
                == F(base['own_median_observed']), 'baseline/panel mismatch: ' + name)
        origin = {k: run.get(k) for k in ('arm', 'source_stage', 'source_file',
                                         'source_line', 'source_sha256')}
        made, why = strictify(run, t1['own_median_observed'], len(edges))
        if made is None:
            rows.append(dict(name=name, status=why, canonical_run=origin,
                             case_path=entry['path'], case_sha256=entry['sha256']))
            continue
        sources, y = made
        require(sum(s['p'] for s in sources) == 1, 'strict weights do not sum to one')
        residuals = s300.residuals_exact(edges, leaves, dvals, sources, y)
        r_fit = max(map(abs, residuals))
        theta1 = t1['theta1']
        rec = dict(name=name, status='ok', k_fit=len(sources),
                   n_leaves=len(leaves), n_edges=len(edges),
                   r_fit=str(r_fit), theta1=str(theta1),
                   separates=bool(r_fit < theta1),
                   relative_width=(float((theta1 - r_fit) / theta1)
                                   if r_fit < theta1 else None),
                   stored_r_dd_over_theta1=run.get('r_dd_over_theta1'),
                   canonical_run=origin, case_path=entry['path'], case_sha256=entry['sha256'],
                   scale=str(t1['own_median_observed']),
                   witness=s300.serial(dict(sources=sources, pendant=y)))
        rows.append(rec)
        if name == EXAMPLE:
            example = dict(rec, n_edges=len(edges), eps=str(EXAMPLE_EPS),
                           lower=2, upper=len(sources),
                           interval_holds=bool(r_fit < EXAMPLE_EPS < theta1),
                           generating_k=(len(s300.truth_configuration(case)[0])
                                         if s300.truth_configuration(case) else None),
                           residuals=[dict(reference=r, residual=str(z))
                                      for r, z in zip(leaves, residuals)],
                           one_source_lower_certificate=example_lower_certificate(
                               edges, leaves, dvals))
            require(F(example['one_source_lower_certificate']['lower_bound']) == theta1,
                    'example lower certificate does not match Theta_1')

    ok = [r for r in rows if r['status'] == 'ok']
    sep = [r for r in ok if r['separates']]
    summary = dict(collection=COLLECTION, arm=ARM, n_cases=len(rows),
                   n_scored=len(ok), n_separating=len(sep),
                   median_relative_width=(statistics.median(
                       sorted(r['relative_width'] for r in sep)) if sep else None))
    require(example is not None and example['interval_holds'],
            'worked example does not certify its stated interval')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('w') as fh:
        json.dump(dict(stage=312, kind='positive_tolerance_intervals',
                       summary=summary, example=example, cases=rows,
                       strictification=dict(delta=str(DELTA),
                           float_conversion='exact binary value via Fraction',
                           scale='registered sorted full-panel observations[len // 2]'),
                       provenance=dict(inputs_sha256=inputs,
                           canonical_rule='stage300_analyze.canonical on evaluated frozen fits; '
                                          'join saved parameters by case, arm, source stage/file/line',
                           theta2_used=False, refitting=False)), fh,
                  indent=1, sort_keys=True)

    if example:
        print('%s\n  %d references, %d edges, strict %d-source witness'
              % (example['name'], example['n_leaves'], example['n_edges'],
                 example['k_fit']))
        print('  r_fit   = %s = %.17g' % (example['r_fit'], float(F(example['r_fit']))))
        print('  Theta_1 = %s = %.17g' % (example['theta1'], float(F(example['theta1']))))
        print('  at eps = %s:  %s  ->  %d <= m(eps) <= %d   (generating K = %s)'
              % (example['eps'], example['interval_holds'],
                 example['lower'], example['upper'], example['generating_k']))
    print('separating candidates: %d of %d scored (%d cases), median relative '
          'width %.4f' % (summary['n_separating'], summary['n_scored'],
                          summary['n_cases'], summary['median_relative_width']))


if __name__ == '__main__':
    main()
