#!/usr/bin/env python3
"""Prepare a local, inspectable research source package; does not publish it."""
from pathlib import Path
from datetime import datetime, timezone
import ast
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
import zipfile

ROOT = Path(__file__).resolve().parents[1]


def main():
    seeds = ['stage42_benchmark', 'stage41_partial_fpt', 'stage41_fixed_pendant_hardness_check',
             'stage41_fixed_pendant_adversarial', 'stage41_boundary', 'check_stage41_partial',
             'check_stage41_partial_independent', 'stage42_benchmark_plot', 'stage42_build_paper',
             'stage42_package', 'stage16_binary']
    seen = set()
    todo = seeds[:]
    while todo:
        name = todo.pop()
        if name in seen:
            continue
        seen.add(name)
        for node in ast.walk(ast.parse((ROOT / 'scripts' / f'{name}.py').read_text())):
            imports = ([node.module.split('.')[0]] if isinstance(node, ast.ImportFrom) and node.module
                       else [alias.name.split('.')[0] for alias in node.names] if isinstance(node, ast.Import) else [])
            todo.extend(dep for dep in imports if (ROOT / 'scripts' / f'{dep}.py').exists())
    files = {f'scripts/{name}.py': ROOT / 'scripts' / f'{name}.py' for name in seen}
    for pattern in ('manuscript/stage42*.tex', 'literature/stage42*.bib', 'results/stage42_benchmark*.json',
                    'results/stage41*.json'):
        for path in ROOT.glob(pattern):
            files[str(path.relative_to(ROOT))] = path
    for name in ('stage42_paper.pdf', 'stage42_main.pdf', 'stage42_supplement.pdf', 'stage42_benchmark.pdf'):
        files[f'exports/{name}'] = ROOT / 'exports' / name
    files['README.md'] = ROOT / 'manuscript/stage42_REPRODUCIBILITY.md'
    files['tools/tectonic-0.17.0/provenance.json'] = ROOT / 'tools/tectonic-0.17.0/provenance.json'
    manifest = {name: hashlib.sha256(path.read_bytes()).hexdigest() for name, path in sorted(files.items())}
    archive = ROOT / 'exports/stage42_review_package.zip'
    with zipfile.ZipFile(archive, 'w', compression=zipfile.ZIP_DEFLATED) as output:
        for name, path in sorted(files.items()):
            output.write(path, 'pdd_research_draft/' + name)
        output.writestr('pdd_research_draft/MANIFEST.json', json.dumps(manifest, indent=2) + '\n')
    with tempfile.TemporaryDirectory(prefix='stage42-package-') as temp:
        with zipfile.ZipFile(archive) as packed:
            assert packed.testzip() is None
            packed.extractall(temp)
        location = Path(temp) / 'pdd_research_draft'
        for name, digest in manifest.items():
            assert hashlib.sha256((location / name).read_bytes()).hexdigest() == digest
        # Isolated import and example with site-packages disabled: tests packaging,
        # not another timing repetition of the frozen benchmark.
        code = """import sys
sys.path.insert(0,'scripts')
from fractions import Fraction as F
from stage41_partial_fpt import solve_partial
import stage42_benchmark
answer=solve_partial([(10,i,F(1)) for i in range(4)],{i:F(5,2) for i in range(3)},cap=1)
assert answer['status']=='complete' and answer['minimum_sources_within_cap']==1
print('stdlib-only packaged import and documented example PASS')
"""
        checked = subprocess.run([sys.executable, '-S', '-c', code], cwd=location,
                                 capture_output=True, text=True, timeout=30)
        assert checked.returncode == 0, checked.stderr
    report = {'stage': 42, 'created_utc': datetime.now(timezone.utc).isoformat(), 'passed': True,
              'archive': str(archive.relative_to(ROOT)), 'bytes': archive.stat().st_size,
              'sha256': hashlib.sha256(archive.read_bytes()).hexdigest(), 'files': len(files),
              'scripts': len(seen), 'contents_sha256': manifest, 'isolated_example': checked.stdout.strip(),
              'scope': 'Local source/document review package only. No repository created, published or submitted.'}
    (ROOT / 'results/stage42_package.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({k: v for k, v in report.items() if k != 'contents_sha256'}, indent=2))


if __name__ == '__main__':
    main()
