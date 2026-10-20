#!/usr/bin/env python3
"""Stage297 shared harness: bank access with hash verification, forked job
execution and result streaming.

The forked-child pattern, the JSON encoder and the SHA helper are IMPORTED from
stage292_harness so that this stage runs jobs exactly the way stages 292-295
did (one child per job, RLIMIT_AS, hard wall, every run retained including
timeouts and deaths).  No pre-existing file is written.
"""
import json
import multiprocessing as mp
import os
import sys
from pathlib import Path

for _k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
           'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_k] = '1'

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
sys.dont_write_bytecode = True

from stage292_harness import sha, encode, isolated_run                    # noqa: E402

OUTDIR = ROOT / 'results/stage297'
BANK_INDEX = ROOT / 'data/stage287_instances/bank/index.json'


def bank_index():
    return json.loads(BANK_INDEX.read_text())['cases']


def bank_cases(families=None, Ks=None):
    out = []
    for entry in bank_index():
        if families is not None and not any(entry['family'].startswith(f) for f in families):
            continue
        if Ks is not None and entry['K_generating'] not in Ks:
            continue
        out.append(entry)
    return out


def load_verified(entry):
    """Load one frozen bank case after checking its registered SHA256."""
    from stage290_erasure_core import load_bank_case
    path = ROOT / entry['path']
    digest = sha(entry['path'])
    assert digest == entry['sha256'], 'bank case hash changed: %s' % entry['name']
    return load_bank_case(path)


def stratified(entries, per_size, key='reference_set_size', seed=0):
    """`per_size` cases for every distinct value of `key`, deterministically."""
    import random
    buckets = {}
    for e in entries:
        buckets.setdefault(e[key], []).append(e)
    chosen = []
    for value in sorted(buckets):
        group = sorted(buckets[value], key=lambda z: z['name'])
        rng = random.Random('%s|%s|%s' % (seed, key, value))
        rng.shuffle(group)
        chosen.extend(group[:per_size])
    return chosen


def _call(payload):
    """Module-level trampoline: run one spec in a forked child (picklable)."""
    function, spec, wall, memory = payload
    res = isolated_run(function, spec, spec.get('wall', wall),
                       spec.get('memory', memory))
    res['spec'] = {k: v for k, v in spec.items() if k != 'record'}
    if 'record' in spec:
        res['spec']['record_name'] = spec['record'].get('name')
    return res


def run_pool(jobs, function, workers, wall, memory, out_path, label='', every=10):
    """Forked-child pool; one JSON line per job, every failure retained."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    payloads = [(function, j, wall, memory) for j in jobs]
    print('%s: %d jobs -> %s' % (label, len(jobs), out_path), flush=True)
    with out_path.open('x') as fh:
        with mp.Pool(workers, maxtasksperchild=1) as pool:
            done = 0
            for res in pool.imap_unordered(_call, payloads, chunksize=1):
                if res['status'] == 'ok':
                    rows = res['value']
                    for row in (rows if isinstance(rows, list) else [rows]):
                        fh.write(json.dumps(encode(row)) + '\n')
                else:
                    fh.write(json.dumps(encode(dict(row_status=res['status'], **res))) + '\n')
                fh.flush()
                done += 1
                if done % every == 0 or done == len(jobs):
                    print('  %s %d/%d' % (label, done, len(jobs)), flush=True)
    print('written %s' % out_path, flush=True)
