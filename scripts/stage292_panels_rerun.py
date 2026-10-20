#!/usr/bin/env python3
"""Stage292 panel re-run: only the panel specs that failed in
results/stage292_panels.jsonl (memory or wall), under the bounded distance
cache. Writes a separate file; the original records are left untouched."""
import json
import multiprocessing as mp
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
sys.dont_write_bytecode = True

import stage292_panels as P                                             # noqa: E402
from stage292_harness import encode, isolated_run                       # noqa: E402

OUT = ROOT / 'results/stage292_panels_rerun.jsonl'
WORKERS = int(os.environ.get('STAGE292_WORKERS', '8'))


def _run(spec):
    wall = 1800.0 if spec['size'] <= 2048 else 12600.0
    res = isolated_run(P.job, spec, wall, 12 << 30)
    res['spec'] = spec
    return res


def main():
    rows = [json.loads(l) for l in (ROOT / 'results/stage292_panels.jsonl').open()]
    jobs = [r['spec'] for r in rows if r.get('kind') != 'panel']
    print('%d failed panel specs to redo' % len(jobs), flush=True)
    with OUT.open('x') as fh:
        with mp.Pool(min(WORKERS, max(1, len(jobs))), maxtasksperchild=1) as pool:
            for res in pool.imap_unordered(_run, jobs, chunksize=1):
                payload = res['value'][0] if res['status'] == 'ok' else res
                fh.write(json.dumps(encode(payload)) + '\n')
                fh.flush()
                print(res['status'], res['spec']['size'], res['spec']['K'],
                      res['spec']['seed'], flush=True)
    print('written', OUT, flush=True)


if __name__ == '__main__':
    main()

