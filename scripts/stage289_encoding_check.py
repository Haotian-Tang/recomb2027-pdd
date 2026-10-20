#!/usr/bin/env python3
"""Stage 289 encoding check: the bank case, written as DecoDiPhy input, is
exactly explainable at the generating support in DecoDiPhy's own parameterisation."""
import json, subprocess, sys
from fractions import Fraction as F
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from stage289_decodiphy_bridge import load_case, write_inputs, residual_vector

ROOT = Path(__file__).resolve().parents[1]
DD = ROOT / 'data/stage289_decodiphy/venv/bin/decodiphy'

def main():
    case_path, workdir = Path(sys.argv[1]), Path(sys.argv[2])
    case = load_case(case_path)
    workdir.mkdir(parents=True, exist_ok=True)
    mapping = write_inputs(case, workdir)
    label_of_edge = {v: k for k, v in mapping['edge_of_label'].items()}
    truth = case['generating_sources']
    anchors = [label_of_edge[s['edge']] for s in truth]
    (workdir / 'anchors.txt').write_text('\n'.join(anchors) + '\n')
    res = residual_vector(case['edges'], [dict(edge=s['edge'], p=F(s['p']), x=F(s['x'])) for s in truth],
                          F(case['generating_pendant']), case['observed'])
    print('case', case['name'], 'leaves', len(case['observed']), 'edges', len(case['edges']), 'K', len(truth))
    print('max |model-observed| at generating truth (float):', max(abs(v) for v in res.values()))
    out = subprocess.run([str(DD), '-t', str(workdir / 'tree.nwk'), '-d', str(workdir / 'distances.txt'),
                          '--score', str(workdir / 'anchors.txt'), '-o', str(workdir)],
                         capture_output=True, text=True, timeout=600)
    print('--- decodiphy --score ---')
    print(out.stdout[-1500:])
    if out.returncode:
        print('STDERR', out.stderr[-2000:])
    print('true p:', [float(F(s['p'])) for s in truth])
    print('true x (from edges[e][0]):', [float(F(s['x'])) for s in truth])
    print('anchor vertex is edges[e][0]?:', [mapping['label_to_vertex'][a] == case['edges'][s['edge']][0]
                                             for a, s in zip(anchors, truth)])
    print('true pendant y:', float(F(case['generating_pendant'])))

if __name__ == '__main__':
    main()
