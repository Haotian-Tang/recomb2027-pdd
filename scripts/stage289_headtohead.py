#!/usr/bin/env python3
"""Stage 289: head-to-head of DecoDiPhy 1.5.0 against the frozen exact solver.

prepare  registers the protocol (cases, input hashes, code hashes, environment,
         matched limits) before any solve.
run      executes every (case, method, repetition) in an isolated forked child
         with an address-space limit and a hard wall limit, retaining every run
         including timeouts.

The two tools answer different questions.  DecoDiPhy minimises a squared
distance residual for a number of placements it either takes as input or selects
by its own stopping rules; the exact solver returns the minimum number of
sources that reproduces the observations exactly, or reports that none exists
within the budget.  Three DecoDiPhy arms make the asymmetry explicit; all use
the tool's documented defaults except for the k control named in the arm.
"""
import os
for _key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
             'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_key] = '1'
import argparse
import contextlib
import hashlib
import io
import json
import multiprocessing as mp
import platform
import resource
import shutil
import signal
import subprocess
import sys
import tempfile
import time
import traceback
from datetime import datetime, timezone
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from stage289_decodiphy_bridge import load_case, write_inputs, to_bank_sources
from stage43_partial_fpt import solve_partial

PROTOCOL = ROOT / 'results/stage289/stage289_headtohead_protocol.json'
STREAM = ROOT / 'results/stage289/stage289_headtohead_runs.jsonl'
BANK_INDEX = ROOT / 'data/stage287_instances/bank/index.json'
BANK_PROTOCOL = ROOT / 'results/stage287_instance_bank_protocol.json'
VENV = ROOT / 'data/stage289_decodiphy/venv'
REPO = ROOT / 'data/stage289_decodiphy/repo'
FAMILIES = ('exact_10k', 'exact_biotrees', 'krepp_estimated_10k')
FROZEN = ('stage43_partial_fpt.py', 'stage41_partial_fpt.py', 'stage27_general_fpt.py',
          'stage26_missing_leaf.py', 'stage25_deficit_extension.py', 'stage23_two_source_envelope.py',
          'stage42_benchmark.py')
NEW = ('stage289_headtohead.py', 'stage289_decodiphy_bridge.py', 'stage289_support_feasibility.py',
       'stage289_repro_check.py')
METHODS = ('exact_stage43_partial', 'dd_default', 'dd_true_k', 'dd_resolution')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def serial(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {str(k): serial(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [serial(x) for x in value]
    return value


def decodiphy_argv(method, case, tree, dist, outdir):
    k_true = len(case['generating_sources'])
    base = ['-t', tree, '-d', dist, '-o', outdir]
    if method == 'dd_default':
        return base
    if method == 'dd_true_k':
        return base + ['-k', str(k_true), '-f', '1']
    if method == 'dd_resolution':
        return base + ['--min_p', '1e-12', '-k', str(k_true + 3), '-f', '0']
    raise ValueError(method)


def run_decodiphy(method, case, work):
    import decodiphy.cli as ddcli
    outdir = Path(work) / method
    outdir.mkdir(parents=True, exist_ok=True)
    argv = decodiphy_argv(method, case, str(Path(work) / 'tree.nwk'),
                          str(Path(work) / 'distances.txt'), str(outdir))
    saved = sys.argv
    sys.argv = ['decodiphy'] + argv
    buf = io.StringIO()
    start = time.perf_counter()
    try:
        with contextlib.redirect_stdout(buf):
            ddcli.main()
    finally:
        sys.argv = saved
    seconds = time.perf_counter() - start
    log = buf.getvalue()
    rounds = json.loads((outdir / 'all_rounds.json').read_text())
    final = rounds[-1]
    mapping = json.loads((Path(work) / 'mapping.json').read_text())
    sources = to_bank_sources(case, final['anchors'], final['p'], final['x'], mapping)
    ladder = [dict(k=r['k'], rounds=r['rounds'], loss=r['loss'], anchors=list(r['anchors']))
              for r in rounds]
    stop = ('low_abundance' if 'Stopped by low abundance' in log else
            'near_zero_residue' if 'Stopped by near-zero residue' in log else 'k_exhausted')
    y = final['y']
    answer = dict(status='complete', argv=argv, k=len(sources), loss=final['loss'],
                  sources=sources, pendant_total=float(sum(y)) if isinstance(y, list) else float(y),
                  pendant_per_placement=list(y) if isinstance(y, list) else [float(y)],
                  stop_reason=stop, ladder=ladder, solver_wall_seconds=seconds)
    shutil.rmtree(outdir, ignore_errors=True)
    return answer


def run_exact(case):
    start = time.perf_counter()
    raw = solve_partial(case['edges'], case['observed'], cap=case['cap'],
                        seconds=None, max_states=None, max_bases=None, max_faces=None)
    answer = dict(status=raw['status'], feasible=raw.get('feasible'),
                  minimum=raw.get('minimum_sources_within_cap'),
                  reason=raw.get('reason'), statistics=raw.get('statistics'),
                  solver_wall_seconds=time.perf_counter() - start)
    if 'original_witness' in raw:
        answer['witness'] = raw['original_witness']
    return answer


def worker(method, case, work):
    if method == 'exact_stage43_partial':
        return run_exact(case)
    return run_decodiphy(method, case, work)


def isolated_run(method, case, work, limit, memory_bytes):
    """One job in a forked child: address-space limit inside, wall limit outside."""
    with tempfile.TemporaryDirectory(prefix='stage289_') as temporary:
        path = Path(temporary) / 'result.json'
        start = time.perf_counter()
        pid = os.fork()
        if pid == 0:
            try:
                resource.setrlimit(resource.RLIMIT_AS, (memory_bytes, memory_bytes))
                answer = worker(method, case, work)
            except BaseException as exc:
                answer = dict(status='error', exception=repr(exc)[:500],
                              traceback=traceback.format_exc()[-2000:])
            try:
                path.write_text(json.dumps(serial(answer)))
            finally:
                os._exit(0)
        timed_out = False
        while True:
            waited, status, usage = os.wait4(pid, os.WNOHANG)
            if waited:
                break
            if time.perf_counter() - start >= limit:
                timed_out = True
                os.kill(pid, signal.SIGKILL)
                _, status, usage = os.wait4(pid, 0)
                break
            time.sleep(0.005)
        elapsed = time.perf_counter() - start
        if timed_out:
            answer = dict(status='unknown', reason='hard_wall_timeout', minimum=None)
        elif path.exists():
            answer = json.loads(path.read_text())
        else:
            answer = dict(status='error', reason='worker_exited_without_result')
        answer.update(process_wall_seconds=elapsed,
                      child_cpu_seconds=usage.ru_utime + usage.ru_stime,
                      maximum_rss_kib=usage.ru_maxrss, wait_status=status)
        return answer


_SETTINGS = {}


def run_one_case(payload):
    index, entry = payload
    case = load_case(ROOT / entry['path'])
    case['cap'] = entry['cap']
    records = []
    with tempfile.TemporaryDirectory(prefix='stage289_case_') as work:
        write_inputs(case, work)
        for repetition in range(_SETTINGS['repetitions']):
            methods = list(METHODS)
            shift = (index + repetition) % len(methods)
            methods = methods[shift:] + methods[:shift]
            if (index // len(methods) + repetition) % 2:
                methods.reverse()
            for method in methods:
                value = isolated_run(method, case, work, _SETTINGS['limit'], _SETTINGS['memory'])
                records.append(dict(name=case['name'], family=entry['family'],
                                    repetition=repetition + 1, method=method,
                                    execution_order=methods, result=value))
    return index, records


def initializer(settings):
    _SETTINGS.update(settings)
    import decodiphy.cli  # imported in the pool worker so no child pays import time


def environment_record():
    freeze = subprocess.run([str(VENV / 'bin/pip'), 'freeze'], capture_output=True, text=True).stdout
    def git(*args):
        return subprocess.run(['git', '-C', str(REPO)] + list(args),
                              capture_output=True, text=True).stdout.strip()
    sources = sorted((REPO / 'decodiphy').glob('*.py'))
    return dict(
        tool='DecoDiPhy', version='1.5.0',
        repository='https://github.com/shayesteh99/DecoDiPhy',
        commit=git('rev-parse', 'HEAD'), describe=git('describe', '--tags'),
        commit_date=git('log', '-1', '--format=%ci'),
        install='python3 -m venv; pip install . from the cloned repository',
        interpreter=subprocess.run([str(VENV / 'bin/python'), '-V'],
                                   capture_output=True, text=True).stdout.strip(),
        pip_freeze=freeze.strip().split('\n'),
        package_file_sha256={p.name: sha(p) for p in sources},
        host=dict(platform=platform.platform(), processor=platform.processor()))


def select_cases():
    index = json.loads(BANK_INDEX.read_text())['cases']
    chosen = [entry for entry in index if entry['family'] in FAMILIES]
    chosen.sort(key=lambda e: (e['family'], e['name']))
    return chosen


def prepare(limit, memory_gib, repetitions, workers):
    assert not PROTOCOL.exists(), 'Never overwrite a registered protocol'
    PROTOCOL.parent.mkdir(parents=True, exist_ok=True)
    entries = select_cases()
    cases = []
    for entry in entries:
        path = ROOT / entry['path']
        digest = sha(path)
        assert digest == entry['sha256'], f'Bank case changed: {entry["name"]}'
        payload = json.loads(path.read_text())
        payload.update(family=entry['family'], bank_path=entry['path'], bank_sha256=digest,
                       reference_set_size=entry['reference_set_size'],
                       K_generating=entry['K_generating'],
                       exactly_explained=entry['exactly_explained'],
                       distance_source=entry['distance_source'])
        cases.append(payload)
    protocol = dict(
        stage=289, registered_utc=datetime.now(timezone.utc).isoformat(),
        purpose='Head-to-head comparison of the published DecoDiPhy heuristic with the frozen '
                'exact structural solver on identical stage-287 bank instances.',
        decodiphy=environment_record(),
        exact_solver=dict(module='stage43_partial_fpt', entry='solve_partial',
                          options='cap=case cap; seconds/max_states/max_bases/max_faces all None, '
                                  'as in the stage-44 confirmation protocol'),
        methods={
            'exact_stage43_partial': 'minimum number of sources within the case budget, exact rational arithmetic',
            'dd_default': 'documented defaults: hill climbing, quick neighbourhood radius 2, osqp, '
                          'seed 1142, min_p 0.01, k selected by DecoDiPhy stopping rules',
            'dd_true_k': 'documented defaults plus -k K_true -f 1, i.e. exactly the generating number of placements',
            'dd_resolution': 'documented defaults plus --min_p 1e-12 -k K_true+3 -f 0, i.e. the '
                             'tool\'s own near-zero-residue stop with the abundance stop effectively disabled',
        },
        decodiphy_defaults_used=dict(method='hill', quick='1', radius=2, seed=1142, fix_k='0',
                                     optimizer='osqp', min_p=0.01),
        asymmetry_note='DecoDiPhy takes the number of placements as an input or selects it by stopping '
                       'rules; the exact solver returns the minimum. The three arms cover its own '
                       'selection, the oracle value of k, and its highest-resolution setting.',
        input_encoding='Exact rational edge lengths and observed distances are written to the Newick '
                       'tree and the distance file as shortest round-trip doubles; DecoDiPhy is a '
                       'floating-point tool. The exact solver reads the unrounded rationals.',
        per_job_hard_wall_seconds=limit, per_job_address_space_bytes=memory_gib * 2 ** 30,
        repetitions=repetitions, parallel_workers=workers,
        bank_index_sha256=sha(BANK_INDEX), bank_protocol_sha256=sha(BANK_PROTOCOL),
        script_hashes={name: sha(ROOT / 'scripts' / name) for name in FROZEN + NEW},
        case_count=len(cases),
        case_index=[dict(name=c['name'], family=c['family'], path=c['bank_path'],
                         sha256=c['bank_sha256'], reference_set_size=c['reference_set_size'],
                         K_generating=c['K_generating'], cap=c['cap'],
                         exactly_explained=c['exactly_explained']) for c in cases])
    PROTOCOL.write_text(json.dumps(protocol, indent=1) + '\n')
    print('registered', PROTOCOL, 'cases', len(cases))


def run():
    assert not STREAM.exists(), 'Never overwrite measured runs'
    protocol = json.loads(PROTOCOL.read_text())
    for name, expected in protocol['script_hashes'].items():
        if name == 'stage289_headtohead.py':
            continue
        assert sha(ROOT / 'scripts' / name) == expected, f'Source changed: {name}'
    for entry in protocol['case_index']:
        assert sha(ROOT / entry['path']) == entry['sha256'], f'Bank case changed: {entry["name"]}'
    settings = dict(limit=protocol['per_job_hard_wall_seconds'],
                    memory=protocol['per_job_address_space_bytes'],
                    repetitions=protocol['repetitions'])
    payloads = list(enumerate(protocol['case_index']))
    started = time.perf_counter()
    done = 0
    with STREAM.open('a') as stream:
        with mp.get_context('fork').Pool(protocol['parallel_workers'], initializer,
                                         (settings,)) as pool:
            for index, records in pool.imap_unordered(run_one_case, payloads, chunksize=1):
                for record in records:
                    stream.write(json.dumps(record) + '\n')
                stream.flush()
                done += 1
                if done % 10 == 0 or done == len(payloads):
                    print(f'{done}/{len(payloads)} cases, {time.perf_counter() - started:.0f}s',
                          flush=True)
    print('finished', done, 'cases in', round(time.perf_counter() - started, 1), 's')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('phase', choices=('prepare', 'run'))
    parser.add_argument('--seconds', type=int, default=300)
    parser.add_argument('--memory-gib', type=int, default=8)
    parser.add_argument('--repetitions', type=int, default=2)
    parser.add_argument('--workers', type=int, default=12)
    args = parser.parse_args()
    if args.phase == 'prepare':
        prepare(args.seconds, args.memory_gib, args.repetitions, args.workers)
    else:
        run()
