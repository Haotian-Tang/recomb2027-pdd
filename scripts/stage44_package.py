#!/usr/bin/env python3
"""Create and smoke-check a local Stage44 review package; never publishes."""
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
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def script_closure(seeds):
    seen, todo = set(), list(seeds)
    while todo:
        name = todo.pop()
        if name in seen:
            continue
        path = ROOT / 'scripts' / f'{name}.py'
        assert path.is_file(), f'Missing script dependency: {path}'
        seen.add(name)
        for node in ast.walk(ast.parse(path.read_text())):
            imports = ([node.module.split('.')[0]] if isinstance(node, ast.ImportFrom) and node.module
                       else [alias.name.split('.')[0] for alias in node.names] if isinstance(node, ast.Import) else [])
            todo.extend(dep for dep in imports if (ROOT / 'scripts' / f'{dep}.py').is_file())
    return seen


def main():
    seeds = [p.stem for p in (ROOT / 'scripts').glob('stage44*.py')]
    # Retain historical reproducibility tools and the explicitly dynamic
    # Stage16 hardness helper, which an AST-only import walk cannot discover.
    seeds += ['stage43_partial_fpt', 'stage43_solver_validation', 'stage43_benchmark',
              'stage43_fixed_one_missing', 'check_stage43_fixed_one_missing',
              'check_stage43_fixed_one_extra', 'stage43_theory_robust_mass',
              'stage43_build_paper', 'stage16_binary',
              'stage41_fixed_pendant_hardness_check', 'stage41_fixed_pendant_adversarial',
              'stage41_boundary', 'check_stage41_partial', 'check_stage41_partial_independent']
    scripts = script_closure(seeds)
    files = {f'scripts/{name}.py': ROOT / 'scripts' / f'{name}.py' for name in scripts}
    excluded = {'stage44_package.json', 'stage44_integrity.json',
                'stage43_package.json', 'stage43_integrity.json'}
    patterns = ['manuscript/stage44*', 'literature/stage44*', 'review/stage44*',
                'results/stage44*', 'logs/stage44*',
                'manuscript/stage43*', 'literature/stage43*', 'review/stage43*',
                'results/stage43*.json', 'logs/stage43_fixed_one_development.txt',
                'results/stage42_benchmark*.json', 'results/stage41*.json']
    for pattern in patterns:
        for path in ROOT.glob(pattern):
            if path.is_file() and path.name not in excluded:
                files[str(path.relative_to(ROOT))] = path
    for name in ['manuscript/stage42_hardness_appendix.tex',
                 'literature/stage42_verified_references.bib',
                 'tools/tectonic-0.17.0/provenance.json']:
        files[name] = ROOT / name
    for name in ['stage44_paper.pdf', 'stage44_main.pdf', 'stage44_supplement.pdf']:
        files[f'exports/{name}'] = ROOT / 'exports' / name
    files['README.md'] = ROOT / 'manuscript/stage44_REPRODUCIBILITY.md'
    assert all(path.is_file() for path in files.values())
    assert not any(name.endswith('/tectonic') for name in files)
    assert 'results/stage44_package.json' not in files and 'results/stage44_integrity.json' not in files

    # Packaging must not silently include PDFs stale relative to the audited build.
    build = json.loads((ROOT / 'results/stage44_paper_build.json').read_text())
    for name, expected in build['source_sha256'].items():
        assert sha(ROOT / name) == expected, f'Manuscript changed since build: {name}'
        assert name in files
    for name, expected in build['pdf_sha256'].items():
        assert sha(ROOT / 'exports' / name) == expected
    for protocol_name in ['stage44_confirmation_protocol.json', 'stage44_confirmation_posthoc_protocol.json']:
        protocol = json.loads((ROOT / 'results' / protocol_name).read_text())
        for name, expected in protocol['script_hashes'].items():
            assert sha(ROOT / 'scripts' / name) == expected
            assert f'scripts/{name}' in files
    cert = json.loads((ROOT / 'results/stage44_certificate_protocol.json').read_text())
    for name, expected in cert['frozen_inputs'].items():
        assert sha(ROOT / name) == expected and name in files

    manifest = {name: sha(path) for name, path in sorted(files.items())}
    archive = ROOT / 'exports/stage44_review_package.zip'
    with zipfile.ZipFile(archive, 'w', compression=zipfile.ZIP_DEFLATED) as packed:
        for name, path in sorted(files.items()):
            packed.write(path, 'pdd_stage44/' + name)
        packed.writestr('pdd_stage44/MANIFEST.json', json.dumps(manifest, indent=2) + '\n')
    smoke_reports = []
    with tempfile.TemporaryDirectory(prefix='stage44-package-', dir='/tmp') as temporary:
        with zipfile.ZipFile(archive) as packed:
            assert packed.testzip() is None
            packed.extractall(temporary)
        location = Path(temporary) / 'pdd_stage44'
        for name, expected in manifest.items():
            assert sha(location / name) == expected
        code = """import sys
sys.path.insert(0, 'scripts')
from fractions import Fraction as F
from stage44_certificate import path_screen
from stage43_fixed_one_missing import solve_fixed_one_missing
edges = [(3, i, F(1)) for i in range(3)]
upper = {0: F(13,8)+F(1,32), 1: F(19,8)+F(1,32), 2: F(19,8)+F(1,32)}
retained, witnesses, _ = path_screen(edges, upper, F(1), F(1,4), F(1,4))
assert retained == [0]
answer = solve_fixed_one_missing(edges, {0:F(2), 1:F(2)}, F(1), cap=2)
assert answer['status'] == 'complete' and answer['minimum_sources_within_cap'] == 2
assert answer['deficit_vertices'] == [3]
print('PASS: standard-library-only certificate example and fixed-one degree-two-gate example')
"""
        commands = [('documented_examples', [sys.executable, '-S', '-c', code]),
                    ('independent_evidence_replay', [sys.executable, '-S', 'scripts/stage44_confirmation_replay.py'])]
        for name, command in commands:
            checked = subprocess.run(command, cwd=location, capture_output=True, text=True, timeout=30)
            assert checked.returncode == 0, (name, checked.stdout, checked.stderr)
            smoke_reports.append({'name': name, 'passed': True, 'stdout': checked.stdout.strip(),
                                  'stderr': checked.stderr.strip()})
        replay = json.loads((location / 'results/stage44_confirmation_replay.json').read_text())
        assert replay['jobs_checked'] == 256 and replay['returned_witnesses_replayed'] == 192
        assert replay['separate_posthoc_witnesses_checked'] == 2
        # Replay is permitted to rewrite only its own audit summary in the
        # isolated extraction. No archived scientific input/result is changed.
        for name, expected in manifest.items():
            if name != 'results/stage44_confirmation_replay.json':
                assert sha(location / name) == expected
    # Detect any concurrent alteration while taking the package snapshot.
    for name, expected in manifest.items():
        assert sha(files[name]) == expected, f'Packaging input changed concurrently: {name}'
    report = dict(stage=44, created_utc=datetime.now(timezone.utc).isoformat(), passed=True,
        archive=str(archive.relative_to(ROOT)), bytes=archive.stat().st_size, sha256=sha(archive),
        files=len(files), scripts=len(scripts), zip_entries=len(files) + 1,
        manifest_verified=True, build_source_and_pdf_hashes_verified=True,
        frozen_protocol_dependencies_verified=True, smoke_checks=smoke_reports,
        isolated_python_no_site_packages=True, no_timing_benchmark_reruns=True,
        missing_required_paths=[], contents_sha256=manifest,
        scope='Local review package only; no publication, external contact or submission. Compiler provenance is included, executable is excluded.')
    (ROOT / 'results/stage44_package.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({k: v for k, v in report.items() if k not in ('contents_sha256', 'smoke_checks')}, indent=2))
    for smoke in smoke_reports:
        print(smoke['name'], 'PASS')


if __name__ == '__main__':
    main()
