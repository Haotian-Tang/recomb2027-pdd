#!/usr/bin/env python3
"""Stage 293 diagnostic: where DecoDiPhy's returned anchor labels land."""
import os, sys, json, io, contextlib, tempfile
for _k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[_k] = '1'
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import stage288_common as common
import stage288_incompleteness as inc
import stage293_abundance as ab
from stage289_decodiphy_bridge import write_inputs

entries = [e for e in inc.panel() if e['K_generating'] == 5 and e['reference_set_size'] < 80]
case = common.load_case(entries[0])
draw = [d for d in inc.removal_draws(case, '0.50') if d['regime'] == 'nearest'][0]
retained = ab.retained_observations(case, draw['removed'])
shim = dict(edges=case['edges'], observed=retained,
            generating_sources=case['generating_sources'])
work = tempfile.mkdtemp(prefix='stage293_diag_')
mapping = write_inputs(shim, work)
print('root vertex:', mapping['root'])
print('n labels with a parent edge:', len(mapping['edge_of_label']))
inv = {v: k for k, v in mapping['label_to_vertex'].items()}
print('root label:', inv.get(mapping['root']))
print("is 'N41' a known label?", 'N41' in mapping['label_to_vertex'],
      "| has parent edge?", 'N41' in mapping['edge_of_label'])
outdir = Path(work) / 'dd'; outdir.mkdir()
import decodiphy.cli as ddcli
sys.argv = ['decodiphy', '-t', str(Path(work)/'tree.nwk'), '-d', str(Path(work)/'distances.txt'),
            '-o', str(outdir)]
buf = io.StringIO()
with contextlib.redirect_stdout(buf):
    ddcli.main()
rounds = json.loads((outdir / 'all_rounds.json').read_text())
final = rounds[-1]
print('anchors:', final['anchors'])
print('x      :', [round(float(v), 4) for v in final['x']])
print('p      :', [round(float(v), 4) for v in final['p']])
unknown = [a for a in final['anchors'] if a not in mapping['edge_of_label']]
print('anchors without a parent edge:', unknown)
for a in unknown:
    print('   ', a, '-> vertex', mapping['label_to_vertex'].get(a),
          '| is root:', mapping['label_to_vertex'].get(a) == mapping['root'],
          '| parent:', mapping['parent'].get(str(mapping['label_to_vertex'].get(a))))
print('log tail:', buf.getvalue()[-400:].replace(chr(10), ' | '))
