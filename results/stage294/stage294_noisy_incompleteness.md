# Stage 294 — Does the reference-incompleteness count collapse survive estimation error?

Registered 2026-10-18T13:43:29.013490+00:00. Protocol `results/stage294/stage294_protocol.json`, draw-plan digest `879fce98bba2023c…`.

## The question, and why it needed a different tool

Stages 288, 290, 292 and 293 measure reference incompleteness on exactly explainable instances in exact rational arithmetic. On real krepp ANI distances the exact solver correctly answers "no exact explanation" and stops, so the same question has to be put to the tool that does run on noisy data — DecoDiPhy 1.5.0, documented defaults — against the benchmark's published ground truth.

The baseline is **not** the truth. DecoDiPhy already mis-estimates the source count on 283 of 296 full-panel krepp cases (stage 289 found the same on 273 of 290), so the quantity isolated here is the **change** induced by reference removal relative to the same tool, same arm, same case, at the full panel. Level `none` is that reference point and was run for every case.

## What was run

* Panel: the complete `krepp_estimated_10k` family, 296 cases, real krepp ANI estimates from reads, N 20–128, K ∈ {3, 5, 10}. No subsampling.
* Draws: 10656 — level `none` plus five removal levels × four regimes, with levels, regimes, replicate counts, the removal-count rule and the three synthetic regimes imported from `stage288_incompleteness` / `stage288_common`; the random-draw seed key uses the paired exact case name, so a random draw is byte-identical to the stage 288 draw wherever the two panels agree (135 of 296 pairs).
* DecoDiPhy arms: `dd_default` (k by its own stopping rules — the arm under test) and `dd_true_k` (`-k K_true -f 1`). 21312 runs over both arms, of which 16399 returned an answer and 4913 failed; the arm under test, `dd_default`, accounts for 8417 answers and 2239 failures; 1343 of the `dd_default` answers were recovered from the tool's own committed `all_rounds.json` after its jplace writer raised, exactly as stage 293 did. Every run is retained.
* Paired exact arm: `exact_matched`, 7200 runs of the frozen `stage43_partial.solve_partial` on the paired exact instance — same tree, same generating sources, same retained reference set, exact distances instead of ANI estimates. K ∈ {3, 5} only: stage 292 established the frozen minimum pass does not complete at K = 10 within 3,600 s, so K = 10 was not attempted rather than burned as timeouts. 200 baseline runs and all 7200 runs completed.
* Both DecoDiPhy workers and the exact worker are imported from `stage293_abundance` unmodified, so the tool invocation is byte-identical to the exact-data stage.

## Answer

**1. On the paired exact instances the stage 288 phenomenon reproduces on this panel.** Of 7000 comparable draws the exact minimum falls in 1979 (28.3%) and rises in 0 — one-sided, as the model requires. Stage 288 measured 450/1968 = 22.9% on its own panel.

**2. On real estimated distances the effect does not carry over as a count collapse.** Over 8121 draws where DecoDiPhy answered both the draw and its own full-panel baseline, the inferred count falls in 36.1% of draws and **rises in 20.9%**, with Δk spanning -18 to +16 and a mean of -0.549. The one-sidedness that makes the exact statement a bias is absent.

**3. The regime ordering inverts.** On the identical draws, random removal is the regime under which the noisy count falls *most* often, and the structured regimes least:

| structured / random ratio in the fall rate | clade | nearest | taxon |
|---|---|---|---|
| exact arm, same draws | 2.39 | 2.69 | 1.28 |
| noisy arm (`dd_default`), same draws | 0.54 | 0.80 | 0.73 |

Stage 288 measured 3.15 (clade) and 2.88 (nearest) over random on exact data; the exact arm here reproduces that ordering, the noisy arm reverses it.

**4. What structured removal reliably produces on noisy data is not a lower count but no answer at all.** DecoDiPhy availability falls from 100 % at the full panel to 48.2% at 75 % loss, and the failures are regime dependent: random 10.0% failed, clade 31.2% failed, nearest 12.2% failed, taxon 28.2% failed. Any rate computed over completed runs alone is conditioned on that survival, which is why the matched comparison below reports both denominators.

**5. At the level of a single draw the noisy answer does not track the exact one.** On the 7000 common draws, when the exact minimum falls DecoDiPhy reproduces the fall in 580 of 1979 such draws (29.3%), gives no answer in 814 (41.1%), goes the opposite way in 329 (16.6%) and reports no change in 256 (12.9%).
  Conversely, where the exact minimum is unchanged DecoDiPhy reports a fall in 1452 of 5021 draws (28.9%).

**Scope this forces on the thesis.** "An incomplete reference set biases the inferred constituent count downward, and which references are missing matters far more than how many" is, on this evidence, a statement about the exact model and the exact minimum. It reproduces exactly, draw for draw, on exact distances. On estimated distances, with the only tool that runs there, the bias is not one-sided, the regime ordering reverses, and the dominant consequence of structured loss is tool failure. Both halves are measured on the same draws, so this is not a panel artefact.

## Count change by level — `dd_default`, all 296 cases (K ∈ {3, 5, 10})

| level | draws compared | fell | unchanged | rose | mean Δk | q10 | q90 | availability |
|---|---|---|---|---|---|---|---|---|
| one | 2011 | 17.6% | 68.4% | 14.0% | -0.014 | -1 | 1 | 97.1% |
| 0.10 | 1923 | 31.5% | 47.4% | 21.2% | -0.222 | -2 | 1 | 92.8% |
| 0.25 | 1728 | 41.1% | 33.5% | 25.4% | -0.371 | -3 | 3 | 83.4% |
| 0.50 | 1460 | 49.5% | 27.3% | 23.3% | -0.912 | -5 | 2 | 70.5% |
| 0.75 | 999 | 54.1% | 23.2% | 22.7% | -2.030 | -9 | 2 | 48.2% |

## Count change by regime — `dd_default`, all 296 cases

| regime | draws compared | fell | unchanged | rose | mean Δk | worst Δk | availability |
|---|---|---|---|---|---|---|---|
| random | 2663 | 41.4% | 40.4% | 18.3% | -0.811 | -17 | 90.0% |
| clade | 2035 | 29.5% | 47.0% | 23.5% | -0.257 | -18 | 68.8% |
| nearest | 1299 | 33.3% | 35.7% | 30.9% | -0.246 | -16 | 87.8% |
| taxon | 2124 | 37.5% | 47.1% | 15.4% | -0.684 | -15 | 71.8% |

`dd_true_k` is fixed at K_true by construction and its Δk is 0 in all 7682 comparable draws; it carries the support and abundance accuracy, not the count.

## Exact versus noisy on the same draws

Common draw set: K in {3, 5}, both arms have a full-panel baseline, the exact arm completed; a DecoDiPhy failure is retained as noisy_no_answer. 7000 draws.

| level | draws | exact: fell | noisy: fell (all draws) | noisy: rose | noisy: no answer | noisy: fell (of answered) |
|---|---|---|---|---|---|---|
| one | 1400 | 0.9% | 18.1% | 13.1% | 2.5% | 18.6% |
| 0.10 | 1400 | 5.6% | 30.1% | 19.4% | 6.1% | 32.0% |
| 0.25 | 1400 | 17.8% | 35.8% | 20.4% | 16.1% | 42.7% |
| 0.50 | 1400 | 46.7% | 34.4% | 16.9% | 29.1% | 48.4% |
| 0.75 | 1400 | 70.4% | 26.8% | 10.3% | 51.1% | 54.7% |

| regime | draws | exact: fell | noisy: fell (all draws) | noisy: rose | noisy: no answer | noisy: fell (of answered) |
|---|---|---|---|---|---|---|
| random | 2000 | 16.4% | 38.0% | 15.7% | 9.6% | 42.0% |
| clade | 2000 | 39.4% | 20.6% | 16.1% | 30.8% | 29.8% |
| nearest | 1000 | 44.3% | 30.5% | 27.8% | 11.1% | 34.3% |
| taxon | 2000 | 21.0% | 27.8% | 10.4% | 27.5% | 38.3% |

Direction agreement over the same draws:

| exact \ noisy | fell | unchanged | rose | no answer |
|---|---|---|---|---|
| fell | 580 | 256 | 329 | 814 |
| unchanged | 1452 | 2122 | 793 | 654 |

The exact arm never raises the minimum, so there is no `rose` row.

## The exact arm in its own right

At the full panel the exact minimum equals K_generating in 198 of 200 cases.

| level | draws | fell | mean Δ(minimum) | worst |
|---|---|---|---|---|
| one | 1400 | 0.9% | -0.009 | -1 |
| 0.10 | 1400 | 5.6% | -0.059 | -2 |
| 0.25 | 1400 | 17.8% | -0.224 | -2 |
| 0.50 | 1400 | 46.7% | -0.684 | -3 |
| 0.75 | 1400 | 70.4% | -1.346 | -4 |

| regime | draws | fell | mean Δ(minimum) | worst |
|---|---|---|---|---|
| random | 2000 | 16.4% | -0.208 | -3 |
| clade | 2000 | 39.4% | -0.734 | -4 |
| nearest | 1000 | 44.3% | -0.752 | -4 |
| taxon | 2000 | 21.0% | -0.306 | -4 |

## Reproducibility of the full-panel baseline

Level `none` repeats the stage 289 full-panel run: same tool, same arm, same case, same defaults, same unremoved panel. Of 575 comparable runs the inferred count agrees in 567 and the returned support in 471. Where the support differs the objective values agree to a median absolute difference of 3.4e-10 on a loss of order 7e-2 (4.8e-9 relative), i.e. the tool is landing on numerically equivalent optima with different supports; k, the quantity under test, is the stable part. This is itself worth reporting: on noisy data the support DecoDiPhy returns is not reproducible across harness-level process details, while its count mostly is. 13 runs have no stage 289 counterpart because that run errored there.

## Accuracy against the published truth (secondary)

| arm | level | runs | support exactly recovered | mean fraction of true support edges | median W1 to truth | median weighted placement error |
|---|---|---|---|---|---|---|
| dd_default | none | 296 | 0 | 0.248 | 0.0844 | 0.0466 |
| dd_default | one | 2011 | 0 | 0.230 | 0.0857 | 0.0492 |
| dd_default | 0.10 | 1923 | 0 | 0.201 | 0.0885 | 0.0539 |
| dd_default | 0.25 | 1728 | 0 | 0.186 | 0.0939 | 0.0579 |
| dd_default | 0.50 | 1460 | 0 | 0.131 | 0.0983 | 0.0663 |
| dd_default | 0.75 | 999 | 0 | 0.097 | 0.1042 | 0.0761 |
| dd_true_k | none | 292 | 0 | 0.185 | 0.0772 | 0.0408 |
| dd_true_k | one | 1976 | 0 | 0.165 | 0.0798 | 0.0435 |
| dd_true_k | 0.10 | 1854 | 0 | 0.152 | 0.0830 | 0.0475 |
| dd_true_k | 0.25 | 1644 | 0 | 0.143 | 0.0846 | 0.0528 |
| dd_true_k | 0.50 | 1382 | 0 | 0.122 | 0.0889 | 0.0581 |
| dd_true_k | 0.75 | 834 | 0 | 0.091 | 0.1008 | 0.0691 |

W1 is `stage289_headtohead_analyze.wasserstein` with `point_distances`, imported unchanged, so the numbers sit on the stage 289 scale (balanced optimal transport under the tree metric, branch-length units). The generating support is never recovered exactly, at any level, in either arm; the mixture error degrades monotonically but mildly with reference loss, roughly 23–31 % from the full panel to 75 % loss, while the count change above is far larger and two-sided.

## The taxonomic removal regime

It exists, and it is built from real taxonomy rather than a proxy.

* Source: Web of Life (WoL), https://github.com/biocore/wol, head commit `ff92d2251c0ea72299a2ead1aecfd748a9acd698` (2023-02-16T20:28:37Z); table `taxonomy_ncbi_ranks.tsv.xz` with 10575 genomes, sha256 `1e443141312cab3e…`. Three further tables (NCBI and GTDB lineages, genome metadata) were fetched and checksummed for provenance but are not used by the regime. Full manifest: `data/stage294_taxonomy/stage294_taxonomy_manifest.json`.
* Match: all 441 distinct accessions used by the noisy panel are in the table, 16631/16631 observations — a complete match, no approximation.
* Rule: drop **whole taxonomic groups** at the finest operable rank, largest-first (replicate 0) or smallest-first (replicate 1), without overshooting the target removal count and never leaving fewer than three references. Rank used: genus in 1920 draws, species in 1040.
* Panel breadth is the binding constraint and is a property of the published benchmark: species spans ≥ 2 groups in all 296 cases, genus in 192, family and order in 134 each, phylum and class in 28 each, kingdom in none. Coarser ranks are left unused by design — a regime defined only where the rank is operable would change the case set between cells and break the draw-for-draw pairing. 162 of 296 panels span a single family.
* Whole-group granularity means the achieved removal count can fall short of the nominal level: the median achieved/target ratio is 1.00 at every level, the mean falls from 0.89 at level 0.10 to 0.73 at 0.75. Both counts are recorded per draw.

**How close is whole-clade removal as a proxy?** Not very, at the level of a draw. The taxonomic draw is itself a clade of the reference tree in 832 of 2960 draws, and the best Jaccard overlap between a taxonomic draw and any whole-clade draw at the same level has median 0.15 (mean 0.27) -- closer than to the nearest-to-source draws (0.04) but no closer than to random draws (0.13). The individual removed groups are mostly monophyletic -- 2713 of 4341 groups of size at least two (62.5 per cent) -- but 16703 of the removed groups are single genomes, so a taxonomic draw is typically a scattered union of small groups rather than one subtree. On the measured behaviour the taxonomic regime sits between random and clade: its exact-arm fall rate (21.0 per cent) is nearer random (16.4 per cent) than clade (39.4 per cent).

## What this does not establish

* It does not show that reference incompleteness is harmless on noisy data. It shows that the *inferred count* does not move downward in a one-sided way, and that reference loss instead produces large two-sided count error and outright tool failure.
* It does not separate "the effect is an artefact of the exact model" from "DecoDiPhy is too noisy an estimator of k to see it". Both are consistent with these data. The second is strongly suggested — its count is wrong on 283 of 296 full-panel cases before any reference is removed, and its support is not even reproducible across harness process ancestry at equal loss — but no experiment here settles it, and no second noisy estimator exists in this release to triangulate with (stage 289: only `--method hill --optimizer osqp` is operable in DecoDiPhy 1.5.0).
* It does not measure the exact minimum on noisy data. There is none: the exact solver answers "no exact explanation" on krepp ANI distances in about 9 ms. Any claim of the form "the true minimum under noise falls" remains unmeasured by any method available here.
* The exact arm covers K ∈ {3, 5} only, 200 of the 296 cases, because the frozen minimum pass does not complete at K = 10. The exact/noisy comparison therefore says nothing about K = 10; the noisy arm alone covers it, and its by-K block shows the same two-sidedness there.
* DecoDiPhy availability is regime dependent, so the by-regime noisy rates are not a clean randomised comparison. The tables give both denominators, but no imputation can recover what the tool did not answer, and the draws it fails on are plausibly the harder ones.
* The taxonomic regime is genus-or-species level only on this benchmark; nothing here tests removal of whole families, orders or phyla, which is the scale at which real reference databases are most obviously incomplete. A panel whose references span several families would be needed.
* Single run per draw. DecoDiPhy at seed 1142 was verified deterministic within stage 289 (575/575 identical across repetitions), so this costs nothing in k, but the support instability measured against stage 289 shows the returned support is a weaker statistic than that determinism suggests.
* No statistical test is reported. The quantities are proportions over a registered draw plan, and the draws within a case are not independent.
* Nothing here speaks to abundance or position error under noise beyond the W1 and placement-error columns, which are truth-referenced, not change-referenced; stage 293 covers those on exact data.

## Files

* `results/stage294/stage294_protocol.json` — registered protocol, input and code hashes, draw-plan digest
* `results/stage294/stage294_protocol_amendment.json` — the one post-registration correction (text only, no measured value, no removal draw and no run input changed): the taxonomy manifest asserted that no rank above genus is operable on this benchmark, which its own breadth table contradicts, and the string was replaced by the measured per-rank counts; the file records both checksums and the superseded text
* `results/stage294/stage294_pilot.json` — timing pilot used to size the wall limits
* `results/stage294/stage294_dd_runs.jsonl` — 21312 DecoDiPhy runs, failures retained
* `results/stage294/stage294_exact_matched_runs.jsonl` — 7200 paired exact runs
* `results/stage294/stage294_dd_metrics.jsonl` — per-run truth-referenced metrics (support, W1, placement error)
* `results/stage294/stage294_analysis.json` — count-change distributions, availability, accuracy, taxon diagnostics
* `results/stage294/stage294_matched_comparison.json` — the common-draw exact/noisy comparison
* `data/stage294_taxonomy/` — fetched WoL tables with URLs, checksums and the match report
* `scripts/stage294_common.py` — bank loader, paired-case resolution, draw plan, taxonomic regime
* `scripts/stage294_noisy.py` — protocol registration and the three run arms
* `scripts/stage294_taxonomy.py` — taxonomy fetch provenance and match report
* `scripts/stage294_analyze.py` — metrics, summary and matched comparison
* `scripts/stage294_report.py` — this document

