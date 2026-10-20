#!/usr/bin/env python3
"""Build and inspect the self-contained Stage42 LaTeX review manuscript.

The official Tectonic release is installed locally and its digest recorded.
PyMuPDF is used only for document checks, page separation and preview images.
No scientific result writer is invoked by this build.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import os
import re
import subprocess
import sys

try:
    import pymupdf as pdf
except ImportError:
    sys.path.insert(0, '/tmp/h2170_pdfdeps')
    import pymupdf as pdf

ROOT = Path(__file__).resolve().parents[1]
EXPORT = ROOT / 'exports'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    EXPORT.mkdir(exist_ok=True)
    (ROOT / 'logs').mkdir(exist_ok=True)
    compiler = ROOT / 'tools/tectonic-0.17.0/tectonic'
    env = dict(os.environ, XDG_CACHE_HOME='/tmp/recomb_stage42_tex_cache')
    command = [str(compiler), '--untrusted', '--keep-logs', '--keep-intermediates',
               '--outdir', str(EXPORT), 'stage42_paper.tex']
    built = subprocess.run(command, cwd=ROOT / 'manuscript', env=env,
                           text=True, capture_output=True, timeout=180)
    (ROOT / 'logs/stage42_paper_build.log').write_text(built.stdout + built.stderr)
    if built.returncode:
        raise RuntimeError(built.stdout + built.stderr)
    log = (EXPORT / 'stage42_paper.log').read_text()
    failures = [line for line in log.splitlines() if re.search(
        r'Overfull|Missing character|undefined references|Citation .* undefined|Reference .* undefined|^!', line)]
    assert not failures, failures
    document = pdf.open(EXPORT / 'stage42_paper.pdf')
    texts = [page.get_text() for page in document]
    assert all(abs(page.rect.width - 612) < .1 and abs(page.rect.height - 792) < .1 for page in document)
    assert not any('\ufffd' in t or '??' in t for t in texts)
    reference_page = next(i for i, text in enumerate(texts) if text.startswith('References'))
    appendix_page = next(i for i, text in enumerate(texts) if re.match(r'A\s+The hardness construction', text))
    body_pages = reference_page - 1  # The first page is separate title/abstract.
    assert body_pages <= 10, body_pages
    clipped = []
    for number, page in enumerate(document, 1):
        for block in page.get_text('dict')['blocks']:
            for line in block.get('lines', []):
                for span in line['spans']:
                    x0, y0, x1, y1 = span['bbox']
                    if x0 < 35 or x1 > page.rect.width - 35 or y0 < 15 or y1 > page.rect.height - 15:
                        clipped.append({'page': number, 'text': span['text'], 'bbox': span['bbox']})
    assert not clipped, clipped
    for name, first, last in [('stage42_main', 0, appendix_page - 1),
                               ('stage42_supplement', appendix_page, len(document) - 1)]:
        part = pdf.open()
        part.insert_pdf(document, from_page=first, to_page=last)
        part.set_metadata({'title': name.replace('_', ' '), 'subject': 'Internal research draft'})
        part.save(EXPORT / f'{name}.pdf', garbage=4, deflate=True)
        part.close()
    preview_dir = EXPORT / 'stage42_preview'
    preview_dir.mkdir(exist_ok=True)
    # Selected full-resolution pages include title, dense math, experiment and appendix.
    indices = sorted({0, 2, min(reference_page - 1, 7), appendix_page,
                      min(appendix_page + 4, len(document) - 1)})
    for i in indices:
        document[i].get_pixmap(matrix=pdf.Matrix(1.4, 1.4)).save(preview_dir / f'page_{i+1:02d}.png')
    # Vector contact sheet for a final all-page layout inspection.
    cols, tile_w, tile_h = 4, 204, 264
    rows = (len(document) + cols - 1) // cols
    contact = pdf.open()
    page = contact.new_page(width=cols * tile_w, height=rows * tile_h)
    for i in range(len(document)):
        left, top = (i % cols) * tile_w, (i // cols) * tile_h
        page.show_pdf_page(pdf.Rect(left, top, left + tile_w, top + tile_h), document, i)
    page.get_pixmap(matrix=pdf.Matrix(1.7, 1.7)).save(preview_dir / 'contact_sheet.png')
    contact.close()
    sources = [ROOT / 'manuscript/stage42_paper.tex',
               ROOT / 'manuscript/stage42_hardness_appendix.tex',
               ROOT / 'manuscript/stage42_experiment_section.tex',
               ROOT / 'literature/stage42_verified_references.bib',
               EXPORT / 'stage42_benchmark.pdf']
    report = {'stage': 42, 'created_utc': datetime.now(timezone.utc).isoformat(), 'passed': True,
              'compiler': 'Tectonic 0.17.0, verified official musl release',
              'compiler_sha256': sha(compiler), 'pymupdf_version': pdf.VersionBind,
              'paper_size_points': [612, 792], 'font_size_points': 10, 'margins_inches': 1,
              'total_pages': len(document), 'separate_title_abstract_pages': 1,
              'main_body_pages': body_pages, 'bibliography_pages': appendix_page - reference_page,
              'supplement_pages': len(document) - appendix_page,
              'unresolved_references': False, 'missing_glyph_warnings': False,
              'overfull_boxes': False, 'clipped_spans': clipped,
              'source_sha256': {str(p.relative_to(ROOT)): sha(p) for p in sources},
              'pdf_sha256': {n: sha(EXPORT / n) for n in
                  ('stage42_paper.pdf', 'stage42_main.pdf', 'stage42_supplement.pdf')},
              'scope': 'Build and layout checks, not proof certification or submission; author metadata remains a draft placeholder.'}
    (ROOT / 'results/stage42_paper_build.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({k: v for k, v in report.items() if not k.endswith('sha256')}, indent=2))


if __name__ == '__main__':
    main()
