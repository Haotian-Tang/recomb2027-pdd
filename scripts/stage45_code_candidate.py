#!/usr/bin/env python3
"""Create a local unpublished scientific-code candidate; no publishing calls."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import zipfile

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    target = ROOT / 'release/stage45_public_candidate'
    archive = ROOT / 'exports/stage45_code_candidate.zip'
    report_path = ROOT / 'results/stage45_code_candidate.json'
    assert not any(p.exists() for p in (target, archive, report_path)), 'No snapshot overwrite'
    old = json.loads((ROOT / 'results/stage44_package.json').read_text())
    omitted_scripts = {'scripts/stage43_build_paper.py', 'scripts/stage44_build_paper.py',
                       'scripts/stage44_package.py'}
    files = {name: ROOT / name for name in old['contents_sha256']
             if name.startswith('scripts/') and name not in omitted_scripts}
    # Keep scientific historical evidence, excluding document/package bookkeeping.
    for name in old['contents_sha256']:
        if name.startswith('results/') and not any(word in name for word in
                ('paper_build', 'document_integrity', 'package.json')):
            files[name] = ROOT / name
    for name in ['scripts/stage45_remaining_checks.py', 'scripts/stage45_remaining_audit.py',
                 'scripts/stage45_build_paper.py', 'results/stage45_remaining_protocol.json',
                 'results/stage45_remaining_checks.json', 'results/stage45_remaining_runs.jsonl',
                 'results/stage45_remaining_integrity.json', 'results/stage45_paper_build.json',
                 'tools/tectonic-0.17.0/provenance.json']:
        files[name] = ROOT / name
    build = json.loads((ROOT / 'results/stage45_paper_build.json').read_text())
    for name, expected in build['source_sha256'].items():
        assert sha(ROOT / name) == expected
        files[name] = ROOT / name
    for name, expected in build['pdf_sha256'].items():
        assert sha(ROOT / 'exports' / name) == expected
        files[f'exports/{name}'] = ROOT / 'exports' / name
    files['README.md'] = ROOT / 'manuscript/stage45_CODE_README.md'
    assert all(path.is_file() for path in files.values())
    assert not any(name.startswith(('review/', 'logs/')) for name in files)
    assert not any(word in name.lower() for name in files for word in
                   ('expert_brief', 'contact', 'outreach', 'email', 'supplementary materials'))
    assert not any(name.endswith('/tectonic') for name in files)
    assert 'LICENSE' not in files
    for name in files:
        if name != 'README.md' and name in old['contents_sha256']:
            assert sha(files[name]) == old['contents_sha256'][name], name
    status = dict(state='NOT YET PUBLISHED', repository_url=None, authors=None, license=None,
                  pending_user_decisions=['Public repository URL and visibility',
                                          'Authorship and affiliation metadata',
                                          'License and release permission'],
                  notice='Local candidate only. No license is granted by this snapshot. No upload, contact or submission has occurred.')
    source_manifest = {name: sha(path) for name, path in sorted(files.items())}
    target.parent.mkdir(exist_ok=True)
    smoke = []
    with tempfile.TemporaryDirectory(prefix='stage45-code-build-', dir='/tmp') as staging:
        staged = Path(staging) / 'pdd_stage45_code'
        for name, path in files.items():
            destination = staged / name
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(path, destination)
        (staged / 'PUBLICATION_STATUS.json').write_text(json.dumps(status, indent=2) + '\n')
        manifest = dict(source_manifest)
        manifest['PUBLICATION_STATUS.json'] = sha(staged / 'PUBLICATION_STATUS.json')
        (staged / 'MANIFEST.json').write_text(json.dumps(dict(sorted(manifest.items())), indent=2) + '\n')
        temporary_archive = Path(staging) / 'stage45_code_candidate.zip'
        with zipfile.ZipFile(temporary_archive, 'w', compression=zipfile.ZIP_DEFLATED) as packed:
            for path in sorted(staged.rglob('*')):
                if path.is_file():
                    packed.write(path, 'pdd_stage45_code/' + str(path.relative_to(staged)))
        with tempfile.TemporaryDirectory(prefix='stage45-code-smoke-', dir='/tmp') as extraction:
            with zipfile.ZipFile(temporary_archive) as packed:
                assert packed.testzip() is None
                packed.extractall(extraction)
            extracted = Path(extraction) / 'pdd_stage45_code'
            example = """import sys
sys.path.insert(0, 'scripts')
from fractions import Fraction as F
from stage43_fixed_one_missing import solve_fixed_one_missing
from stage44_certificate import path_screen
edges = [(3, i, F(1)) for i in range(3)]
answer = solve_fixed_one_missing(edges, {0:F(2), 1:F(2)}, F(1), cap=2)
assert answer['status'] == 'complete' and answer['minimum_sources_within_cap'] == 2
upper = {0:F(13,8)+F(1,32), 1:F(19,8)+F(1,32), 2:F(19,8)+F(1,32)}
retained, certificates, diagnostics = path_screen(edges, upper, F(1), F(1,4), F(1,4))
assert retained == [0]
print('PASS: exact fixed-pendant and conditional path-screen examples')
"""
            commands = [('read_only_audit', [sys.executable, '-S', 'scripts/stage45_remaining_audit.py', '--audit']),
                        ('documented_examples', [sys.executable, '-S', '-c', example])]
            for name, command in commands:
                run = subprocess.run(command, cwd=extracted, capture_output=True, text=True, timeout=30)
                assert run.returncode == 0, (name, run.stdout, run.stderr)
                smoke.append(dict(name=name, passed=True, stdout=run.stdout.strip(), stderr=run.stderr.strip()))
            audit = json.loads(smoke[0]['stdout'])
            assert audit['all_returned_witnesses_replayed'] == 197
            assert audit['paired_original_and_structural_answer_agreements'] == 48
            assert audit['files_written'] == audit['solver_calls'] == 0
            for name, expected in manifest.items():
                assert sha(extracted / name) == expected, name
        for name, expected in source_manifest.items():
            assert sha(files[name]) == expected, f'Concurrent source change: {name}'
        shutil.copytree(staged, target)
        shutil.copyfile(temporary_archive, archive)
    report = dict(stage=45, created_utc=datetime.now(timezone.utc).isoformat(), status='PASS',
                  state='NOT YET PUBLISHED', candidate=str(target.relative_to(ROOT)),
                  archive=str(archive.relative_to(ROOT)), archive_sha256=sha(archive),
                  archive_bytes=archive.stat().st_size, files=len(manifest)+1,
                  scripts=sum(name.startswith('scripts/') for name in manifest),
                  uncompressed_bytes=sum(path.stat().st_size for path in target.rglob('*') if path.is_file()),
                  manifest_sha256=sha(target / 'MANIFEST.json'), contents_sha256=manifest,
                  smoke_checks=smoke, python_no_site_packages=True, network_used=False,
                  frozen_stage44_manifest_hashes_preserved=True, no_timing_reruns=True,
                  excluded=['review and contact/email materials', 'expert brief',
                            'user-provided predecessor supplementary PDF', 'old review PDFs',
                            'compiler binary and TeX cache', 'outer outreach packagers'],
                  pending_metadata=status['pending_user_decisions'], license_granted=False,
                  publication_or_external_contact=False)
    report_path.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({k: v for k, v in report.items() if k not in ('contents_sha256', 'smoke_checks')}, indent=2))
    for check in smoke:
        print(check['name'], 'PASS')


if __name__ == '__main__':
    main()
