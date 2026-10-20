#!/usr/bin/env python3
"""Stage 289 diagnostic: characterise the DecoDiPhy failures seen in the run.

Re-runs the affected cases outside the registered protocol and records whether
the search result survived the crash (all_rounds.json present) and what the
command-line tool does end to end.  Diagnostic only; not a protocol measurement.
"""
import json, subprocess, sys
from collections import Counter
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from stage289_decodiphy_bridge import load_case, write_inputs
DD = str(ROOT / 'data/stage289_decodiphy/venv/bin/decodiphy')
OUT = ROOT / 'results/stage289/stage289_failure_diagnostic.json'

def main():
    protocol = json.loads((ROOT / 'results/stage289/stage289_headtohead_protocol.json').read_text())
    index = {e['name']: e for e in protocol['case_index']}
    failures = {}
    for line in (ROOT / 'results/stage289/stage289_headtohead_runs.jsonl').open():
        r = json.loads(line)
        if r['result']['status'] == 'error':
            failures.setdefault(r['name'], {})[r['method']] = r['result'].get('exception')
    records = []
    for name, methods in sorted(failures.items()):
        entry = index[name]
        case = load_case(ROOT / entry['path'])
        work = Path('/tmp/stage289_diag') / name
        work.mkdir(parents=True, exist_ok=True)
        write_inputs(case, work)
        k_true = len(case['generating_sources'])
        for method, exception in sorted(methods.items()):
            extra = {'dd_default': [], 'dd_true_k': ['-k', str(k_true), '-f', '1'],
                     'dd_resolution': ['--min_p', '1e-12', '-k', str(k_true + 3), '-f', '0']}[method]
            out = work / method
            out.mkdir(exist_ok=True)
            proc = subprocess.run([DD, '-t', str(work / 'tree.nwk'), '-d', str(work / 'distances.txt'),
                                   '-o', str(out)] + extra, capture_output=True, text=True, timeout=900)
            rounds = out / 'all_rounds.json'
            partial = sorted(out.glob('all_rounds_*.json'))
            payload = dict(name=name, family=entry['family'], leaves=entry['reference_set_size'],
                           k_true=k_true, method=method, protocol_exception=exception,
                           exit_code=proc.returncode,
                           stderr_tail=proc.stderr.strip().split('\n')[-1] if proc.stderr.strip() else '',
                           all_rounds_written=rounds.exists(),
                           intermediate_rounds_written=[p.name for p in partial],
                           output_jplace_written=(out / 'output.jplace').exists())
            if rounds.exists():
                data = json.loads(rounds.read_text())
                payload['final_k'] = data[-1]['k']
                payload['final_loss'] = data[-1]['loss']
                payload['final_anchors'] = data[-1]['anchors']
            records.append(payload)
            print(json.dumps({k: v for k, v in payload.items() if k != 'final_anchors'}), flush=True)
    OUT.write_text(json.dumps(dict(stage=289, purpose='Diagnostic for the DecoDiPhy failures in the '
                                   'registered run; re-run outside the protocol.',
                                   failures=len(records), records=records), indent=1) + '\n')
    print(Counter(r['protocol_exception'][:40] for r in records))

if __name__ == '__main__':
    main()
