# Exact phylogenetic distance deconvolution

![reproductions](https://github.com/Haotian-Tang/recomb2027-pdd/actions/workflows/reproductions.yml/badge.svg)

Solvers, instance bank and measurement records for

> **Report the tolerance, not the count: exact phylogenetic distance deconvolution
> with incomplete references.** Haotian Tang, Mohamed bin Zayed University of
> Artificial Intelligence.
> Preprint: [arXiv:2610.30776](https://arxiv.org/abs/2610.30776) \[q-bio.PE\],
> [doi:10.48550/arXiv.2610.30776](https://doi.org/10.48550/arXiv.2610.30776).
> The supplementary appendices are an ancillary file on the arXiv record.

## The problem

A reference tree `T` with known, positive, rational edge lengths is given. `k`
sources sit strictly inside `k` distinct original edges, carry weights
`p_i > 0` summing to one, and share a common mean pendant `y > 0`. At an
observed reference leaf `r` the model predicts

```
d_r = y + sum_i p_i * dist_T(r, q_i)
```

Given `T` and the observations `{d_r}` on a reference subset `R`, what is the
smallest number of sources that explains them? This repository holds the exact
rational algorithms that answer that question, the frozen instance bank they
were run on, and the saved records behind every measurement in the paper.

Three things are measured, and the repository is organised around them:

1. **The exact minimum source count** - an exact rational algorithm, with a
   complexity classification (fixed-parameter tractable with a free mean
   pendant, W[1]-hard with a prescribed pendant and arbitrary missing
   references).
2. **The tolerance `Theta_k`** - for a support `S` of `k` edges,
   `theta(S)` is the exact minimum over the closed mass set of
   `max_r |pred_r - d_r|`, one rational linear program; `Theta_k` is its
   minimum over all size-`k` supports, and the minimum count within a tolerance
   is `m(eps) = min{k : Theta_k < eps}`.
3. **What incomplete reference sets cost**, in the exact model and on real
   estimated distances.

All arithmetic is exact: `fractions.Fraction` throughout, no floating point in
any decision. **The solvers, the tolerant program, the `Theta_1` closed form
and both reproductions below need nothing beyond the Python 3.12 standard
library.**

## Layout

```
scripts/      170 modules: the frozen solvers, their helpers, every
              measurement and analysis driver named in the table below, and
              one standalone proof check. Flat and importable by bare module
              name - each script resolves the repository root as the parent
              of scripts/.
data/         the frozen instance bank (the 1,218 cases the paper's claims
              rest on), the placement-distance cases, and the Web of Life
              taxonomy tables used by the taxonomic removal regime.
results/      the saved records: protocols, per-run streams and the analysis
              JSONs that the paper's numbers are read from.
figures/      the paper's six figures, the plotting code that exists, and
              under figures/data/ the payload one of them is redrawn from.
repro/        two end-to-end reproductions (see below).
provenance/   the manuscript's own claim-to-record maps, carried over as is.
```

## Install

```bash
git clone https://github.com/Haotian-Tang/recomb2027-pdd.git
cd recomb2027-pdd
python3 --version          # 3.12 or newer
```

Nothing else is needed for the reproductions below. The measurements were run
on CPython 3.12.3, x86_64 Linux, glibc 2.39.

For the bank-construction and plotting scripts only:

```bash
python3 -m venv .venv && . .venv/bin/activate
pip install -r requirements.txt     # numpy 2.5.4, scipy 1.18.1, matplotlib 3.11.2
```

`requirements.txt` says which script needs which package. No solver depends on
any of them.

## Reproductions

All three run from a clean checkout with no third-party packages installed.
Measured on one core of an x86_64 Linux host (62 cores, 94 GiB, no GPU use).

### A. The exact solver check - 6.0 s at --n 20, 4 min 40 s at the default 120

```bash
python3 -B repro/repro_a_exact_vs_tolerant.py --n 20     # 6.0 s
python3 -B repro/repro_a_exact_vs_tolerant.py            # 120 leaf sets, 4 min 40 s
```

Reproduces: *"At eps = 0 the tolerant program returns the frozen exact solver's
minimum on all 1,118 leaf sets checked against it"* (results section) and
*"the tolerant program and the frozen exact minimum agree on 1,118 of 1,118
leaf sets, 0 disagreements"* (supplement).

It re-runs **both** registered code paths on a deterministic, cheapest-first
sample of the 1,118 registered leaf sets - the frozen structural solver
(`stage43_partial_fpt.solve_partial`) and the eps-tolerant rational program at
`eps = 0` (`stage296_tolerance.theta_worker`, `m(0)` branch) - and checks that
they agree with each other *and* each reproduces the value in the frozen
evidence streams `results/stage296/stage296_validate_runs.jsonl` and
`results/stage296/stage296_cross_runs.jsonl`.

### B. A paper number on real estimated distances - 5.6 s

```bash
python3 -B repro/repro_b_theta1_vs_residual.py           # all 296 cases, 5.6 s
```

Reproduces: *"The worst-leaf residual of the true configuration exceeds
`Theta_1` in 296 of 296 ANI cases, median ratio 2.9"* (results section).

For each of the 296 krepp ANI cases it recomputes `Theta_1` from the case file
with stage 298's exact closed form and the worst-leaf residual with stage 298's
frozen formula, then checks the count and all five quantiles of the ratio
against the registered record
`results/stage298/stage298_analysis.json`.

### C. The two-source corollary check - 0.44 s

```bash
python3 -B scripts/verify_two_source_corollary.py
```

An independent exact-arithmetic check of the corollary added to supplementary
Appendix E, *Two-source support identifiability*: on a reduced tree with every
leaf observed and exact distances, the occupied original edges of an
at-most-two-source explanation are determined by the observations, with no
assumption that the weights are equal or known. The script draws random reduced
trees, places one or two strict interior sources with **unequal** random
rational weights, and verifies the four steps of the proof in turn; a control
subdivides one edge to confirm that the reduced hypothesis is load-bearing.
It prints

```
reduced trees, k in {1,2}, unequal rational weights: 400/400 verified
failures: none
control with a degree-two internal vertex: f != h in 96/300 draws
```

which are the two figures the supplement quotes. Standard library only.

## Paper claim to script to record

Every claim in the paper that cites a number, with the script that produced it and the record
it was read from: [docs/CLAIMS.md](docs/CLAIMS.md).

## Figures

`figures/` holds the six figures of the submission as vectors, plus the
plotting code that exists. Main-text and supplement figures are numbered
independently, as the two documents have separate `.aux` files.

* `model_projection.pdf` (**main-text Figure 1**, model schematic) - produced by
  `figures/support/make_model_figure.py`. Run
  `python3 -B figures/support/make_model_figure.py`; it overwrites
  `figures/model_projection.pdf` and `.png`.
* `incompleteness.pdf` (**main-text Figure 2**, "What incomplete references
  cost") - produced by `figures/make_fig_incompleteness.py`, which redraws all
  three panels from `figures/data/fig_incompleteness_data.json` and calls no
  solver. Run it from the repository root with
  `python3 figures/make_fig_incompleteness.py`; it needs numpy and matplotlib
  and overwrites `figures/incompleteness.pdf` and `.png`. Because matplotlib
  stamps a creation date into the PDF, a re-run shows those two tracked files
  as modified even when the content is unchanged.
* `performance.pdf` (**supplement Figure 3**, frozen confirmation) - produced by
  `figures/support/make_figures.py`, which reads only the archived records
  under `figures/support/figure_archive/` and calls no solver. Run it with
  `python3 -B figures/support/make_figures.py`; it overwrites
  `figures/performance.pdf` and `.png` and writes an ignored audit JSON under
  `figures/support/`.
* `tolerance_axis.pdf` (**main-text Figure 3**), `real_incompleteness.pdf`
  (**supplement Figure 1**) and `real_data_axis.pdf` (**supplement Figure 2**) -
  **the plotting code for these three was written ad hoc and was not retained
  in the project tree.** The figures are included as vectors, the rasters they
  were wrapped from are in `results/stage296/`, `results/stage299/` and
  `results/stage298/`, and the numeric payloads behind every plotted value are
  `results/stage296/stage296_payload.json`,
  `results/stage299/stage299_regime_ratios.json` and
  `results/stage298/stage298_figure_data.json`. A reader can check any plotted
  number against those payloads; redrawing the panels pixel for pixel would
  need the plotting code rewritten.


## Bulk data archive

The repository carries the analysis records the paper's numbers are read from.
The per-run streams behind them - every individual solver run, with its status,
timing and witness - are 446.7 MB uncompressed and are published as a separate
archive rather than tracked in git.

```
recomb2027_pdd_bulk_v1.tar.gz     62.07 MB compressed, 446.69 MB unpacked
                                  4,498 files
sha256  e1646b3a44a3ec7513127b7a84d0e8f6532a60bc6cb44fa5a2668234e1b9a492
```

`bulk_manifest.json` in this repository lists every member with its byte count
and SHA-256. To fetch and place it:

```bash
# 1. download the v1.0.0 release asset into the repository root
curl --fail --location --output recomb2027_pdd_bulk_v1.tar.gz \
  https://github.com/Haotian-Tang/recomb2027-pdd/releases/download/v1.0.0/recomb2027_pdd_bulk_v1.tar.gz
# 2. check it
sha256sum -c <<<"e1646b3a44a3ec7513127b7a84d0e8f6532a60bc6cb44fa5a2668234e1b9a492  recomb2027_pdd_bulk_v1.tar.gz"
# 3. unpack from the repository root - member paths are already root-relative,
#    so the files land beside the git-tracked ones under results/ and data/
tar -xzf recomb2027_pdd_bulk_v1.tar.gz -C .
# 4. verify every member against the manifest
python3 - <<'PY'
import hashlib, json, pathlib
m = json.load(open('bulk_manifest.json'))
bad = []
for e in m['members']:
    p = pathlib.Path(e['path'])
    if not p.exists():
        bad.append((e['path'], 'missing')); continue
    h = hashlib.sha256(p.read_bytes()).hexdigest()
    if h != e['sha256']:
        bad.append((e['path'], 'sha256'))
print(f"{len(m['members']) - len(bad)} of {len(m['members'])} members verified")
for b in bad[:20]: print(' ', *b)
PY
```

The archive contains, by path:

* `results/stage28*…stage301/*_runs.jsonl` and the other per-run streams - the
  individual solver runs behind the aggregates in this repository.
* `results/stage299/stage299_protocol.json` (71.7 MB) - the registered leaf-set
  plan of the real-data removal measurement, 24,226 leaf sets written out in
  full.
* `results/stage299/stage299_per_draw.csv` (7.4 MB) - one row per draw.
* `data/stage287_instances/bank/cases/published_biotrees/` (4,172 cases) and
  `.../published_true_10k/` (296 cases) - the published-distance bank families.
  No claim in the paper rests on them; they are here so
  `data/stage287_instances/bank/index.json`, which hashes all 5,686 registered
  cases, verifies completely.

Nothing in the archive is needed by either reproduction, or by any number in
the claim table that is not marked **(bulk)**.

## What is deliberately excluded, and why

**The parent authors' implementation.** The project tree carries a clone of
<https://github.com/shayesteh99/DecoDiPhy> and its virtual environment (546 MB)
under `data/stage289_decodiphy/`. Republishing someone else's code inside this
repository would be wrong, so it is excluded entirely. The comparison was run
against commit `a26aa618aaf50eaba80d12bb3e3d71a201d8617e` (2026-07-10), MIT
licensed; install it from its own repository as `requirements.txt` describes.
The two implementation properties the paper reports were measured against that
commit and have since been fixed upstream.

**The upstream benchmark downloads** (136 MB): `demixing_10k`, `krepp_shared`,
`biotrees`, `sixteen_k`, `travelers_diarrhea` and `refset_benchmark` under
`data/stage287_instances/`. These are third-party published data, not ours to
redistribute. `data/stage287_instances/download_manifest.json` records every
source URL and hash, and `scripts/stage287_fetch_instances.py` re-fetches them.
The bank *derived* from them - the instances every measurement actually ran on -
is published here in full.

**Abandoned research threads.** The project has about 300 numbered stages;
stages 42-45 and 286-301 back this paper and stages 1-41 and 46-285 do not.
Their scripts and the 2.1 GB of records they produced are excluded. Notably
`results/stage291_removal_runs.jsonl` alone is 421 MB and belongs to a
screening thread no paper claim uses.

**Virtual environments, build outputs and scratch directories** - `.venv/`,
`__pycache__/`, `stage288_tmp/`, `stage295_tmp/`, `scratch/`, the vendored
`tectonic` and `minimap2` binaries under `tools/`, and the manuscript build
directories. Covered by `.gitignore`.

**The manuscript sources.** This repository is the code and data release; the
paper is the paper.

## Sizes, measured

| Part | Files | Bytes | MB |
|---|---|---|---|
| `scripts/` | 170 | 1,793,137 | 1.79 |
| `data/` (bank subset, placement cases, taxonomy) | 1,285 | 93,275,403 | 93.28 |
| `results/` (analysis records) | 547 | 57,695,549 | 57.70 |
| `figures/` | 30 | 2,291,627 | 2.29 |
| `repro/`, `provenance/`, `README.md` and the other root files | 11 | 1,327,861 | 1.33 |
| **git-tracked total** | **2,043** | **156,383,577** | **156.38** |
| bulk archive (release asset) | 4,498 | 62,072,358 compressed | 62.07 compressed / 446.69 unpacked |

The byte column is exact for this initial release and sums to the total;
the MB column is each row rounded independently. All 2,039 files from the
prepared public folder are retained. Four missing historical confirmation
source files were restored byte-for-byte against their frozen protocol
hashes, bringing the release to 2,043 files. Documentation and figure-output
paths were corrected; saved scientific inputs, results, and figures were
preserved. The historical confirmation driver refuses to overwrite its
archived protocol and results.

Largest tracked file: `results/stage288_incompleteness.json`, 13.40 MB - below
GitHub's 50 MB warning and 100 MB refusal. The 29 stool-sample case files are
54.33 MB of the bank subset; they are kept in git because they are the only
panel on which the paper's real-data headline can be re-measured.

## Provenance

`provenance/experiment_provenance.json` and `provenance/claim_evidence_map.json`
are the manuscript's own machine-readable maps from claim to protocol, result
and raw-stream hash, carried over unchanged. They were written on 17 October
2026 and therefore cover the frozen-confirmation lineage; the table above is
current as of the submitted manuscript and supersedes them where they differ.

Every protocol JSON under `results/` records the script hashes frozen before
its runs, the registered seeds and budgets, and the machine. Every per-run
stream retains all runs, including failures and timeouts.

## Citing this work

Cite the paper for the results and this repository for the code and records.
`CITATION.cff` in the root carries both, so GitHub's "Cite this repository"
button produces them directly.

```bibtex
@article{tang2026tolerance,
  title   = {Report the tolerance, not the count: exact phylogenetic distance
             deconvolution with incomplete references},
  author  = {Tang, Haotian},
  journal = {arXiv},
  year    = {2026},
  eprint  = {2610.30776},
  archivePrefix = {arXiv},
  primaryClass  = {q-bio.PE},
  doi     = {10.48550/arXiv.2610.30776},
  url     = {https://arxiv.org/abs/2610.30776}
}
```

## Licence

MIT, see `LICENSE`. Copyright (c) 2026 Haotian Tang.
