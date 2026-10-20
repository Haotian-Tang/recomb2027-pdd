#!/usr/bin/env python3
"""Create and inspect a local Stage43 source/review package. No publishing."""
from pathlib import Path
from datetime import datetime, timezone
import ast
import hashlib
import json
import subprocess
import sys
import tempfile
import zipfile

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    seeds = ['stage43_partial_fpt', 'stage43_solver_validation', 'stage43_benchmark',
             'stage43_fixed_one_missing', 'check_stage43_fixed_one_missing',
             'check_stage43_fixed_one_extra', 'stage43_theory_robust_mass',
             'stage43_build_paper', 'stage43_package', 'stage16_binary',
             'stage41_fixed_pendant_hardness_check', 'stage41_fixed_pendant_adversarial',
             'stage41_boundary', 'check_stage41_partial', 'check_stage41_partial_independent']
    seen, todo = set(), list(seeds)
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
    for pattern in ['manuscript/stage43*', 'literature/stage43*', 'review/stage43*',
                    'results/stage43*.json', 'logs/stage43_fixed_one_development.txt',
                    'results/stage42_benchmark*.json', 'results/stage41*.json']:
        for path in ROOT.glob(pattern):
            if path.name not in {'stage43_package.json', 'stage43_integrity.json'}:
                files[str(path.relative_to(ROOT))] = path
    for name in ['manuscript/stage42_hardness_appendix.tex', 'literature/stage42_verified_references.bib',
                 'tools/tectonic-0.17.0/provenance.json']:
        files[name] = ROOT / name
    for name in ['stage43_paper.pdf', 'stage43_main.pdf', 'stage43_supplement.pdf']:
        files[f'exports/{name}'] = ROOT / 'exports' / name
    files['README.md'] = ROOT / 'manuscript/stage43_REPRODUCIBILITY.md'
    manifest = {name: sha(path) for name, path in sorted(files.items())}
    archive = ROOT / 'exports/stage43_review_package.zip'
    with zipfile.ZipFile(archive, 'w', compression=zipfile.ZIP_DEFLATED) as packed:
        for name, path in sorted(files.items()):
            packed.write(path, 'pdd_stage43/' + name)
        packed.writestr('pdd_stage43/MANIFEST.json', json.dumps(manifest, indent=2) + '\n')
    with tempfile.TemporaryDirectory(prefix='stage43-package-') as temp:
        with zipfile.ZipFile(archive) as packed:
            assert packed.testzip() is None
            packed.extractall(temp)
        location = Path(temp) / 'pdd_stage43'
        for name, digest in manifest.items():
            assert sha(location / name) == digest
        code = """import sys
sys.path.insert(0,'scripts')
from fractions import Fraction as F
from stage43_fixed_one_missing import solve_fixed_one_missing
import stage43_benchmark
import check_stage43_fixed_one_extra
answer=solve_fixed_one_missing([(3,i,F(1)) for i in range(3)],{0:F(2),1:F(2)},F(1),cap=2)
assert answer['status']=='complete' and answer['minimum_sources_within_cap']==2
assert answer['deficit_vertices']==[3]
print('PASS stdlib-only imports and documented degree-two-gate example')
"""
        checked = subprocess.run([sys.executable, '-S', '-c', code], cwd=location,
                                 capture_output=True, text=True, timeout=30)
        assert checked.returncode == 0, checked.stderr
    report = {'stage': 43, 'created_utc': datetime.now(timezone.utc).isoformat(), 'passed': True,
              'archive': str(archive.relative_to(ROOT)), 'bytes': archive.stat().st_size,
              'sha256': sha(archive), 'files': len(files), 'scripts': len(seen),
              'contents_sha256': manifest, 'isolated_example': checked.stdout.strip(),
              'scope': 'Local review package only. No publication, external contact or submission.'}
    (ROOT / 'results/stage43_package.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({k: v for k, v in report.items() if k != 'contents_sha256'}, indent=2))


if __name__ == '__main__':
    main()
