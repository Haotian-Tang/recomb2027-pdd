#!/usr/bin/env python3
"""Stage 300: register the travellers'-diarrhoea wall and record the censored runs.

The 29 stool samples are two orders of magnitude larger than the krepp cases
(10,575 references, 21,147 edges against 20-129 and 39-255).  Under this
stage's invocation - DecoDiPhy's documented defaults, `-t tree -d distances`,
one forked child per sample - no sample returned an answer within hours, so a
wall is registered here and the runs that exceeded it are recorded as censored
with the cost they consumed.  The numbers below are measured, not estimated:
elapsed time and resident set size were read from the process table at the
moment each run was stopped.

WALL.  7,200 s (2 h) of wall-clock per sample.  Justification: 5 concurrent
workers is what the measured footprint (about 10 GiB resident per worker, still
rising) allows on this shared machine, and 29 samples x 2 h / 5 workers is
about 11.6 h - the largest per-sample wall that keeps the whole collection
inside a single day of the stage's 40-core budget.  Every run in flight had
reached or passed that wall when it was stopped.

This is a statement about the cost of THIS invocation on THIS machine, not
about the tool: the parent publication ran DecoDiPhy on this very collection.
DecoDiPhy's README documents `--warm_start [OUTDIR]/all_rounds.json` precisely
because its search "can be time extensive", i.e. the intended workflow for
expensive inputs is to resume a search across runs rather than to drive one run
to the low-abundance stop, which is what this stage did.
"""
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import stage300_common as s300                   # noqa: E402

WALL_SECONDS = 7200.0
STOPPED_UTC = '2026-10-19T11:02:40+00:00'

# measured with `ps -eo pid,ppid,etimes,rss` at 2026-10-19T11:02:29Z, immediately
# before the processes were killed; rss is the resident set size at that instant,
# which is a lower bound on the peak.
OBSERVED = [
    dict(name='td_sample_101_published', group='pilot', pid=1571104,
         elapsed_seconds=9922, rss_kib=11367548, workers_in_flight=1,
         registered_wall_at_launch=14400.0),
    dict(name='td_sample_10_published', group='batch', pid=1573299,
         elapsed_seconds=7286, rss_kib=10450144, workers_in_flight=5,
         registered_wall_at_launch=21600.0),
    dict(name='td_sample_125_published', group='batch', pid=1573300,
         elapsed_seconds=7286, rss_kib=10469684, workers_in_flight=5,
         registered_wall_at_launch=21600.0),
    dict(name='td_sample_147_published', group='batch', pid=1573301,
         elapsed_seconds=7286, rss_kib=10466116, workers_in_flight=5,
         registered_wall_at_launch=21600.0),
    dict(name='td_sample_156_published', group='batch', pid=1573302,
         elapsed_seconds=7286, rss_kib=10432704, workers_in_flight=5,
         registered_wall_at_launch=21600.0),
    dict(name='td_sample_160_published', group='batch', pid=1573303,
         elapsed_seconds=7286, rss_kib=10408532, workers_in_flight=5,
         registered_wall_at_launch=21600.0),
]
PAIRING = ('the five batch children are paired with the first five samples of the '
           '28-sample queue in queue order (scripts/stage300_rundd.py builds jobs in '
           'sorted case order and the pool launches them in that order); the pilot '
           'child ran the single sample its command line names')


def main():
    rows = []
    for entry in OBSERVED:
        rows.append(dict(
            kind='dd_run_censored', stage=300, collection='travelers',
            name=entry['name'], method='dd_default',
            result=dict(status='unknown', reason='stopped_at_stage300_wall',
                        wall_seconds_registered=WALL_SECONDS,
                        elapsed_seconds=entry['elapsed_seconds'],
                        rss_kib_at_stop=entry['rss_kib'],
                        returned_anything=False),
            harness_wall_seconds=entry['elapsed_seconds'],
            maximum_rss_kib=entry['rss_kib'], timed_out=True,
            workers_in_flight=entry['workers_in_flight'],
            registered_wall_at_launch=entry['registered_wall_at_launch'],
            stopped_utc=STOPPED_UTC, pid=entry['pid'], pairing=PAIRING))
    path = s300.OUT / 'stage300_travelers_censored.jsonl'
    with open(path, 'w') as handle:
        for row in rows:
            handle.write(json.dumps(s300.serial(row)) + '\n')
    elapsed = [r['harness_wall_seconds'] for r in rows]
    rss = [r['maximum_rss_kib'] / 1048576.0 for r in rows]
    summary = dict(
        stage=300, step='travelers_cost', registered_utc=time.strftime(
            '%Y-%m-%dT%H:%M:%S+00:00', time.gmtime()),
        wall_seconds=WALL_SECONDS, wall_justification=__doc__.split('WALL.')[1].split(
            'This is a statement')[0].strip(),
        n_samples_total=29, n_samples_attempted=len(rows), n_samples_completed=0,
        n_samples_not_attempted=29 - len(rows),
        elapsed_seconds=dict(min=min(elapsed), max=max(elapsed),
                             total=sum(elapsed)),
        core_hours_consumed=sum(elapsed) / 3600.0,
        rss_gib_at_stop=dict(min=min(rss), max=max(rss)),
        configuration=dict(
            invocation='decodiphy -t tree.nwk -d distances.txt -o OUTDIR (documented '
                       'defaults: --method hill, --quick 1, --radius 2, --seed 1142, '
                       '--fix_k 0, --optimizer osqp, --min_p 0.01)',
            worker='stage293_abundance.dd_worker, imported unmodified',
            concurrency='1 (pilot) and 5 (batch), chosen from the measured resident set '
                        'size on a shared 62-core machine with other users\' jobs running',
            memory_limit='RLIMIT_AS 60 GB (pilot) and 30 GB (batch); neither was reached',
            threads='OMP/OPENBLAS/MKL_NUM_THREADS=1'),
        measured_input_stage=dict(
            source='results/stage300/stage300_tdprobe_td_sample_101_published.json',
            read_tree_newick_seconds=0.0994, preprocess_input_seconds=1.793,
            distance_matrix_seconds=40.243, distance_matrix_entries=111820050,
            get_input_matrices_not_returned_after_seconds=3777,
            note='treeswift\'s leaf-to-leaf distance matrix alone materialises 111.8 '
                 'million entries; the enclosing get_input_matrices had not returned '
                 'after 63 minutes when the probe was stopped'),
        pruning=dict(preprocess_input_pruned=False, n_leaves_kept=10575,
                     note='unlike the krepp ANI cases, every backbone leaf carries a '
                          'distance here, so the input stage keeps the full panel'),
        what_this_is_not=(
            'This is the cost of this stage\'s invocation on this machine, not a '
            'property of DecoDiPhy: the parent publication ran it on this same '
            'collection.  Its README documents --warm_start because the search "can be '
            'time extensive", i.e. for inputs of this size the intended workflow resumes '
            'a search across runs instead of driving one run to the low-abundance stop. '
            'No conclusion about the tool\'s behaviour on these samples is drawn here.'),
        what_is_reported_instead=(
            'Theta_1 and Lambda_1 for all 29 samples, which need no DecoDiPhy run: '
            'Theta_1 is imported from stage 298 and Lambda_1 was computed here in 1.1 s '
            'per sample.'))
    (s300.OUT / 'stage300_travelers_cost.json').write_text(json.dumps(summary, indent=1))
    print('wrote', path, 'and', s300.OUT / 'stage300_travelers_cost.json')
    print('attempted %d of 29 samples, completed %d, %.1f core-hours consumed, '
          'elapsed %d-%d s, rss %.1f-%.1f GiB'
          % (summary['n_samples_attempted'], summary['n_samples_completed'],
             summary['core_hours_consumed'], min(elapsed), max(elapsed),
             min(rss), max(rss)))


if __name__ == '__main__':
    main()
