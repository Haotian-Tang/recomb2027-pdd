# Exact phylogenetic distance deconvolution

Solvers, instance bank and measurement records for

> **Report the tolerance, not the count: exact phylogenetic distance deconvolution
> with incomplete references.** Haotian Tang, Mohamed bin Zayed University of
> Artificial Intelligence. Manuscript for RECOMB 2027.

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

Every claim in the measurement section (`sections/experiments.tex`) and in the
measurements appendix (`appendices/experiments_details.tex`) that cites a
number. Scripts are under `scripts/`, records under `results/` and `data/`.
Records marked **(bulk)** are in the data archive, not in this repository - see
"Bulk data archive" below.

### Instances

| Paper claim | Script | Record |
|---|---|---|
| 893 exactly explainable profiles over 23 published reference phylogenies, 20-361 observed leaves, `K` in {2,3,5,7,10,12}, plus one 15,941-leaf backbone | `stage287_fetch_instances.py` -> `stage287_pdd_cases.py` -> `stage287_build_bank.py`; census `stage287_tree_census.py`; hashes re-checked by `stage287_verify_bank.py` | `data/stage287_instances/bank/index.json`; `results/stage287_instance_bank_protocol.json`; `results/stage287_tree_census.json`; `results/stage287_instance_bank_verification.json` |
| 296 ANI and 60 placement profiles from the parent benchmark, each paired with the exact profile of the same tree and taxa | `stage286_benchmark_parse.py`, `stage287_build_bank.py`, `stage298_placement.py` | `results/stage286_benchmark_parse.json`; bank families `krepp_estimated_10k`, `exact_10k`; `data/stage298_placement/`; `results/stage298/stage298_placement_manifest.json` |
| 29 travellers'-diarrhoea stool samples on a 10,575-leaf, 21,147-edge 16S backbone with no ground truth | `stage287_build_bank.py` | bank family `published_travelers_diarrhea`; `results/stage298/stage298_analysis.json` `/travelers_rows` |
| Published mixtures exist only to 361 references; every larger panel carries planted sources under a documented deterministic rule | `stage287_build_bank.py`, `stage292_panels.py` | `results/stage287_instance_bank_protocol.json`; `results/stage292_panels.jsonl` |

### The tolerant program

| Paper claim | Script | Record |
|---|---|---|
| At `eps = 0` the tolerant program and the frozen exact minimum agree on 1,118 of 1,118 leaf sets, 0 disagreements | `stage296_tolerance.py` phases `validate` and `cross`; frozen path `stage43_partial_fpt.solve_partial` | `results/stage296/stage296_validate_runs.jsonl`; `results/stage296/stage296_cross_runs.jsonl`; `results/stage296/stage296_analysis.json` `/counts/validate` — **reproduction A** |
| Supports enumerated exhaustively, 6,611,481 of them | `stage296_common.py` (`theta_min_over_supports`), `stage296_tolerance.py` | `results/stage296/stage296_analysis.json` `/ceiling/total_supports` |
| Every reported witness re-evaluated through the independent `min(d(r,u)+xL, d(r,v)+(1-x)L)` form, 0 mismatches | `stage296_common.py` (`independent_deviation`) | `results/stage296/stage296_analysis.json`; `results/stage296/stage296_cross_runs.jsonl` per-row `independent_agrees` |

### One-sidedness

| Paper claim | Script | Record |
|---|---|---|
| Nothing rose in 23,545 removal draws on the real collections | `stage299_small.py`, `stage299_td.py`, `stage299_analyze.py` | `results/stage299/stage299_analysis.json` `/monotonicity/n_draws`; raw `results/stage299/stage299_small_runs.jsonl` **(bulk)** |
| Nothing rose in 1,075 exact draws at each of ten tolerances | `stage296_tolerance.py` phase `cross`, `stage296_analyze.py` | `results/stage296/stage296_analysis.json` `/within_eps_table[*]/n_draws`; `results/stage296/stage296_cross_runs.jsonl` |
| Nothing rose in 502 paired comparisons on a ladder reaching `k = 8` (411 falls, 0 rises) | `stage301_grid.py`, `stage301_analyze.py` | `results/stage301/stage301_analysis.json`; `results/stage301/stage301_tolerant_localisation.md`; raw `results/stage301/stage301_grid_runs_v2.jsonl` **(bulk)** |
| The exact minimum fell in 450 of 1,968 comparable draws on 101 real reference trees, by up to six sources (Table, exact-model removal) | `stage288_incompleteness.py`, `stage288_analyze.py` | `results/stage288_analysis.json` `/incompleteness/undercount_B_versus_A/{B_below_A, comparable}`, `/incompleteness/A_full/cases` |
| Fall rate by level: 0.2%, 4.9%, 16.2%, 37.1%, 59.1% (1/414, 20/406, 63/389, 139/375, 227/384) | `stage288_analyze.py` | `results/stage288_analysis.json` `/incompleteness/undercount_B_versus_A/by_level/{one,0.10,0.25,0.50,0.75}` |
| By regime: whole clade 32.5% (252/776), nearest to sources 29.7% (115/387), uniformly at random 10.3% (83/805) | `stage288_incompleteness.py`, `stage288_analyze.py` | `results/stage288_analysis.json` `/incompleteness/undercount_B_versus_A/by_regime/{clade,nearest,random}` |
| Those ratios hold between 2.0 and 4.0 at every non-vacuous tolerance (Figure 3b) | `stage296_analyze.py` | `results/stage296/stage296_analysis.json` `/within_eps_table`; `results/stage296/stage296_payload.json` |
| Whole-taxon removal is indistinguishable from random on ANI (1.1x at 75% loss) and below random on placement | `stage294_taxonomy.py`, `stage299_analyze.py` | `results/stage299/stage299_regime_ratios.json`; `results/stage299/stage299_analysis.json` `/by_arm_level_regime` |
| Its granularity undershoots the nominal level, a median 0.733 where 0.75 was asked (0.386 on placement trees) | `stage299_analyze.py` | `results/stage299/stage299_analysis.json` `/by_arm_level_regime[*]/removed_fraction_median` |
| Panel breadth confines the taxonomic regime to genus (1,920 draws) and species (1,040) | `stage294_taxonomy.py` | `results/stage294/stage294_analysis.json`; `data/stage294_taxonomy/` |
| Four regimes at five levels; draw plans imported from the first incompleteness protocol, checked element for element on all 1,350 random draws of the 135 cases whose estimated panel equals its exact panel | `stage288_common.py` (`clade_removals`, `nearest_source_removal`, `random_removal`), `stage299_common.py` | `results/stage288_incompleteness_protocol.json`; `results/stage288_incompleteness.json` `/removal_draws`; `results/stage299/stage299_validation.json` |
| Figure 2, "What incomplete references cost" - (a) median fall in `Theta_1` against the share of references removed, exact arm, by regime; (b) the same draws on real ANI estimates of the same trees and the identical leaf sets; (c) references erasable before the minimum-size explanation stops being unique, against panel size, pooled over `K` = 2-10. Caption counts: 592 draws per regime and level, 296 for nearest to sources; (c) pools 723 adversarial and 1,336 random thresholds | `figures/make_fig_incompleteness.py` (redraws it from the payload; needs numpy and matplotlib) | payload `figures/data/fig_incompleteness_data.json`, derived from `results/stage299/stage299_regime_ratios.json` (`fall` block) and `results/stage297/stage297_analysis.json` (`threshold_vs_panel` block); output `figures/incompleteness.{pdf,png}` |

### Uniqueness

| Paper claim | Script | Record |
|---|---|---|
| Erasing a median of two references nearest the sources leaves a second support of minimum size, at every `K` measured | `stage290_erasure_core.py`, `stage292_all_minimisers.py`, `stage297_reduced_all_minimisers.py` | `results/stage292_analysis.json`; `results/stage297/stage297_analysis.json` `/k10/thresholds_stage297/t_unknown`; `results/stage292_support_uniqueness.md` |
| The all-minimisers routine agrees with brute force on 5,839 of 5,839 comparisons; of those, 5,339 replay the recorded minimiser sets and the remaining 500 are fresh instances checked against an independent binomial enumeration | `stage297_validate.py` | `results/stage297/stage297_analysis.json` `/validation/track_A/rows` (5,339 replay) and `/validation/track_B/rows` (500 fresh, `brute_subsets_tested` 148,233); raw `results/stage297/stage297_validate_replay.jsonl` **(bulk)** |
| One of them a case whose complete minimiser set contains 4,443 supports | `stage292_validate.py` | `results/stage292_analysis.json` `/validation/track_A/n_supports_distribution/4443` |
| Adversarial threshold rank correlation with panel size 0.203 at `p = 0.15` over 52 cases spanning 20 to 353 references | `stage297_analyze.py` | `results/stage297/stage297_analysis.json` `/k10/thresholds_stage297/spearman_t_vs_panel_size/K10\|adversarial_nearest` |
| Random loss tolerates a median of 7 at `rho = 0.710`, `p = 1.9e-9`; whole clades 8 at `rho = 0.819`, `p = 2.4e-15` | `stage297_analyze.py` | same key, `K10\|random` and `K10\|whole_clade`; medians at `/k10/thresholds_stage297/t_unknown` |
| On planted panels to 47,151 leaves the random threshold grows linearly in panel size | `stage292_panels.py`, `stage292_analysis.py` | `results/stage292_analysis.json` `/panels`; `results/stage292_panels.jsonl` |
| At `K = 10`, 67 of 99 have minimum = `K` and 55 a unique true support, against 663 of 676 at `K <= 7` | `stage297_bridge.py`, `stage297_analyze.py`; `stage292_bridge.py` | `results/stage297/stage297_analysis.json` `/k10/{cases, baseline_minimum_equals_K, baseline_unique_true}`; `results/stage292_analysis.json` `/bridge/{cases, baseline_unique_true}` |
| Erasable fraction 2.1% pooled over `K <= 7` and 1.4% at `K = 10`, against 42.2% and 9.7% at random | `stage297_analyze.py` | `results/stage297/stage297_analysis.json` `/k10/thresholds_stage297/fraction_erasable` |

### The tolerance

| Paper claim | Script | Record |
|---|---|---|
| With no reference removed the minimum first drops below the exact count at a median tolerance 0.0147 of the case's median observed distance, over 43 cases | `stage296_tolerance.py`, `stage296_analyze.py` | `results/stage296/stage296_analysis.json` `/threshold_distribution`, `/full_panel[*]/first_drop_eps_rel` |
| It tracks the smallest constituent (`rho = 0.665`, `p = 1.1e-6`), not the panel size (`rho = 0.187`, `p = 0.23`) | `stage296_analyze.py` | `results/stage296/stage296_analysis.json` `/threshold_scaling/vs_min_p/spearman`, `/threshold_scaling/vs_N/spearman` |
| Real ANI distances deviate by 0.34 of a median observed distance at the worst leaf, 7.3 to 14.5 times that tolerance in 15 of 15 cases | `stage296_tolerance.py` phase `calibrate` | `results/stage296/stage296_calibration.json`; `results/stage296/stage296_analysis.json` |
| At a tolerance holding that error the model's own minimum is one source in all 15, against three at tolerance zero in 14 of them | `stage296_analyze.py` | `results/stage296/stage296_analysis.json`; `results/stage296/stage296_tolerance_titration.md` |
| Figure 3 (tolerance axis) | plotting code not retained - see "Figures" | `results/stage296/stage296_payload.json` (every plotted value); raster `results/stage296/stage296_tolerance_axis.png`; vector `figures/tolerance_axis.pdf` |

### Real distances

| Paper claim | Script | Record |
|---|---|---|
| The worst-leaf residual of the true configuration exceeds `Theta_1` in 296 of 296 ANI cases, median ratio 2.9 | `stage298_sweep.py`, `stage298_residuals.py`, `stage298_analyze.py`; closed form in `stage298_common.py` | `results/stage298/stage298_analysis.json` `/by_family/krepp_estimated_10k/residual_max_over_theta1`; `results/stage298/stage298_residuals.json` — **reproduction B** |
| and in 50 of 60 placement cases | `stage298_placement.py`, `stage298_analyze.py` | `results/stage298/stage298_analysis.json` `/by_family` placement block |
| At 75% loss the median fall in `Theta_1` is 0.34 at random against 0.63 whole clade and 0.61 nearest (1.9x and 1.8x), where the paired exact arm on identical leaf sets gives 2.9x and 2.6x (Figure 2) | `stage299_small.py`, `stage299_exact_baseline.py`, `stage299_analyze.py` | `results/stage299/stage299_analysis.json` `/by_arm_level_regime`; `results/stage299/stage299_regime_ratios.json`; per-draw table `results/stage299/stage299_per_draw.csv` **(bulk)** |
| On the 10,575-reference backbone, discarding 7,931 references at random costs 1.8% of `Theta_1` and comparably many as whole clades costs 90% (89.6%) | `stage299_td.py`, `stage299_analyze.py` | `results/stage299/stage299_analysis.json` `/by_arm_level_regime` travellers rows; `results/stage299/stage299_real_incompleteness.md`; raw `results/stage299/stage299_td_runs.jsonl` **(bulk)** |
| The attaining edge is never terminal on a real stool sample, in none of 29 full panels and none of 725 draws | `stage298_td.py`, `stage299_td.py` | `results/stage298/stage298_analysis.json` `/travelers/n_argmin_terminal`, `/travelers_rows[*]/argmin_is_terminal`; `results/stage299/stage299_analysis.json` |
| Nearest-to-attainer removal has already moved the attaining edge in 76% of stool draws at 25% loss; whole-clade removal of three quarters moves it in every draw, a median 22 edges away | `stage299_td.py`, `stage299_analyze.py` | `results/stage299/stage299_analysis.json` `/by_arm_level_regime[*]/same_edge_fraction` and hop medians; `results/stage299/stage299_real_incompleteness.md` |
| `Theta_2` on the stool backbone costs a measured 0.638 s per support over `C(21147,2)` supports, 39,594 core-hours for one sample | `stage298_td_k2cost.py` | `results/stage298/stage298_td_k2_cost.json` |
| A certified two-sided bracket yields a median `m(eps) >= 9`, but at `eps` about 1.75e-4 of a median observed distance | `stage301_common.py`, `stage301_prune.py`, `stage301_analyze.py` | `results/stage301/stage301_analysis.json`; `results/stage301/stage301_tolerant_localisation.md` |
| Dropping the support-size restriction fits these data essentially exactly, at a median 0.0019 of `Theta_1` | `stage301_thetainf.py` | `results/stage301/stage301_analysis.json`; raw `results/stage301/stage301_dual_runs.jsonl` **(bulk)** |
| Supplement Figures 1 and 2 (real incompleteness; real data axis) | plotting code not retained - see "Figures" | `results/stage299/stage299_regime_ratios.json`, `results/stage298/stage298_figure_data.json`; rasters `results/stage299/stage299_real_incompleteness.png`, `results/stage298/stage298_real_data_axis.png`; vectors `figures/real_incompleteness.pdf`, `figures/real_data_axis.pdf` |

### The published heuristic

The comparison runs DecoDiPhy 1.5.0 at pinned commit
`a26aa618aaf50eaba80d12bb3e3d71a201d8617e` (2026-07-10). That code is **not**
vendored here; see "What is excluded".

| Paper claim | Script | Record |
|---|---|---|
| On 675 exactly explainable cases the exact method is a median 45.6 times cheaper than DecoDiPhy 1.5.0 | `stage289_headtohead.py`, `stage289_exact_cost_probe.py`, `stage289_headtohead_analyze.py` | `results/stage289/stage289_decodiphy_comparison.md` (the 45.6x line); `results/stage289/stage289_headtohead_percase.json`; `results/stage289/stage289_exact_cost_probe.json`; raw `results/stage289/stage289_headtohead_runs.jsonl` **(bulk)** |
| Whole-support recovery table: `K` = 2, 3, 5, 7, 10 at a 120 s wall - exact 100/99.1/98.6/90.8/2.8%, own `k` 99.2/96.4/80.9/63.0/18.3%, true `k` 99.2/97.7/94.5/93.3/61.0%, resolution arm 99.2/96.8/90.5/90.8/55.4% | `stage289_headtohead.py`, `stage289_headtohead_analyze.py` | `results/stage289/stage289_decodiphy_comparison.md` (recovery table); `results/stage289/stage289_headtohead_percase.json`; `results/stage289/stage289_headtohead_summary.json` |
| At `K = 10` the heuristic given the true `k` still places 94% of individual sources on the correct edge | `stage293_abundance.py`, `stage293_analyze.py` | `results/stage293/stage293_analysis.json`; `results/stage293/stage293_abundance_cost.md` |
| Inferred `k` identical across six seeds in 125 of 132 cases; fixed-seed control agrees 132 of 132 | `stage297_seeds.py`, `stage297_dd_analyze.py`, `stage297_cross_harness.py` | `results/stage297/stage297_dd_analysis.json`; raw `results/stage297/stage297_dd_seeds_main.jsonl` **(bulk)** |
| The returned mixture fits better than the generating configuration in 282 of 296, in the worst-leaf norm and in its own sum of squares | `stage300_evaluate.py`, `stage300_analyze.py`; baselines in `stage300_common.py` | `results/stage300/stage300_analysis.json` `/by_collection/krepp_ani/arms/dd_default/{n_dd_beats_truth_in_linf, n_dd_beats_truth_in_ssr}` |
| In the sum of squares it beats the best single source in 280 of 296 at a median 0.40 of the single-source optimum, while the true configuration beats one source in only 2 of 296 | `stage300_common.py` (`best_single_source_l2`), `stage300_analyze.py` | `results/stage300/stage300_analysis.json` `.../n_dd_beats_one_source_in_ssr`, `.../ssr_over_lambda1/median` |
| Its mixture is no better than the best single source at the worst leaf in 121 of 296, flat in `K` (38/44/41% at `K` = 3, 5, 10) and panel-size dependent (14 of 98, 55 of 98, 52 of 100) | `stage300_analyze.py` | `results/stage300/stage300_analysis.json` `.../n_r_dd_ge_theta1`, `.../by_K`, `.../by_panel_tertile` |
| Re-solving the tool's own quadratic programme at its own support moves that count to 116 of 290; restricting both sides to the references it used moves it to 107 of 296 | `stage300_resolve.py`, `stage300_routee.py`, `stage300_analyze.py` | `results/stage300/stage300_analysis.json` `/by_collection/krepp_ani/resolved/dd_default/n_linf_ge_theta1`, `.../ordering_linf` |
| The configuration written out is not always the one whose loss is printed, the two differing by more than 1e-6 relative in 496 of 1,738 runs | `stage300_harvest.py`, `stage300_verify.py`, `stage300_analyze.py` | `results/stage300/stage300_analysis.json`; raw `results/stage300/stage300_harvest.jsonl` **(bulk)**, `results/stage300/stage300_baselines_krepp_ani.jsonl` |
| The input stage deletes reference leaves that do carry a distance, keeping 3 or 4 of 100-107 references in 14 of 296 cases | `stage300_prune.py` | `results/stage300/stage300_prune_krepp_ani.jsonl` (fields `n_observed`, `n_kept`) |
| On the 60 placement cases the mixture fits better than one source everywhere | `stage300_evaluate.py`, `stage300_analyze.py` | `results/stage300/stage300_analysis.json` `/by_collection/krepp_placement` |
| 6 of 29 stool samples attempted under documented defaults, 0 returned an answer within 7,200 s | `stage300_rundd.py`, `stage300_censor.py` | `results/stage300/stage300_analysis.json`; `results/stage300/stage300_report_payload.json` |

### Supplementary corollary, Appendix E

| Paper claim | Script | Record |
|---|---|---|
| Corollary (Two-source support identifiability): on a reduced tree with every leaf observed and exact distances, the occupied original edges of an at-most-two-source explanation are determined, with no assumption that the weights are equal or known. "Checked on 400 random reduced trees with unequal rational weights in exact arithmetic, all in agreement" | `scripts/verify_two_source_corollary.py` | self-contained, no stored record: prints `400/400 verified` and `failures: none` in 0.44 s on the standard library alone |
| "The reduced hypothesis is doing work: after subdividing one edge, `f` ceased to equal `h` in 96 of 300 draws" | `scripts/verify_two_source_corollary.py`, degree-two control | prints `control with a degree-two internal vertex: f != h in 96/300 draws` |

### Verification ledgers (supplement)

| Paper claim | Script | Record |
|---|---|---|
| 24,226 of 24,226 leaf sets completed with no failure and no timeout | `stage299_analyze.py` | `results/stage299/stage299_analysis.json` `/inputs/n_records`, `/verification/{n_runs, n_complete}` |
| 2,638,368 edges decided twice, frozen per-support program and closed form, 0 disagreements | `stage299_validate.py`, `stage299_small.py` | `results/stage299/stage299_analysis.json` `/verification/double_decided_edges` |
| 15,831 further spot checks on the stool backbone, 0 disagreements; 385 of 385 full-panel values reproduced | `stage299_td.py`, `stage299_analyze.py` | `results/stage299/stage299_analysis.json` `/verification/td_frozen_spot_checks` |
| On the full-panel axis, 653,807 edges decided twice with 0 disagreements and 859 of 859 attaining witnesses re-verified | `stage298_td.py`, `stage298_sweep.py`, `stage298_analyze.py` | `results/stage298/stage298_analysis.json` `/verification/{witnesses_checked, witnesses_agree}`; `results/stage298/stage298_td_summary.json`; raw `results/stage298/stage298_td_runs.jsonl` **(bulk)** |
| 4,656 of 4,656 returned witnesses verified on the original tree; all 563 prescribed-pendant witnesses reproduced the prescribed pendant | `stage288_common.py` (`verify_witness`), `stage288_analyze.py`, `stage43_fixed_one_missing.py` | `results/stage288_analysis.json` |
| The two implementations agree on `Theta_1` in 23,472 of 23,472 leaf sets and on the attaining coordinate in 11,930 | `stage299_validate.py`, `stage299_analyze.py` | `results/stage299/stage299_analysis.json` `/verification`; `results/stage298/stage298_witness_degeneracy.json` |

### Frozen confirmation of the exact implementations (supplement appendix)

| Paper claim | Script | Record |
|---|---|---|
| 48 registered inputs (32 free-`y`, 16 prescribed-`y` with one missing leaf), seed 4410142026, 16 transitive source hashes frozen before solving; all 256 primary runs retained | `stage44_confirmation.py` | `results/stage44_confirmation_protocol.json`; `results/stage44_confirmation.json`; `results/stage44_confirmation_integrity.json` |
| Free `y`: structural 32/32 vs enumeration 30/32; median casewise ratio 2.19 on 30 joint completions (26 vs 4 wins). Prescribed `y`: 13/16 vs 14/16, median ratio 1.48 on 12 joint completions | `stage44_confirmation_analyze.py` | `results/stage44_confirmation.json` |
| All 192 primary returned witnesses passed separate original-tree distance checks | `stage44_confirmation_replay.py` | `results/stage44_confirmation_replay.json` |
| Two posthoc 30 s follow-ups supplied the missing comparisons, yielding agreement on all 48 inputs | `stage44_confirmation_posthoc.py`, `stage45_remaining_checks.py` | `results/stage44_confirmation_posthoc{,_protocol}.json`; `results/stage45_remaining_{checks,protocol,integrity}.json` |
| 12 binary reduction inputs and two integer rescalings matched the predicted source thresholds; all 14 oracle calls completed (7 yes, 7 no) | `stage50_confirmation.py`, `stage50_binary_reduction.py`, `stage50_closed_screened.py`, `stage50_closed_oracle.py` (all match the protocol's frozen source hashes) | `results/stage50_confirmation_protocol.json`; `results/stage50_confirmation.json` |
| Development suites (not used for any reported timing) | `stage42_benchmark.py`, `stage43_benchmark.py`, `stage43_solver_validation.py`, `stage43_fixed_one_missing.py` | `results/stage42_benchmark{,_protocol}.json`; `results/stage43_benchmark{,_protocol}.json`; `results/stage43_solver_validation.json`; `results/stage43_fixed_one{,_extra}_audit.json` |
| Supplement Figure 3 (frozen confirmation) | `figures/support/make_figures.py` | `figures/performance.pdf`; archive under `figures/support/figure_archive/` |

### Solver lineage

The solvers themselves, which every row above runs through:

| Component | Module |
|---|---|
| Exact minimum source count with a free positive mean pendant and at least two observed original leaves | `scripts/stage43_partial_fpt.py` (`solve_partial`) |
| Prescribed pendant with one missing leaf | `scripts/stage43_fixed_one_missing.py` |
| Robust-mass theory helpers | `scripts/stage43_theory_robust_mass.py` |
| Reduction certificate / closed support oracle | `scripts/stage44_certificate.py` |
| Bank access, removal regimes, independent witness verifier, forking job pool | `scripts/stage288_common.py` |
| Complete minimiser sets at a given count | `scripts/stage292_all_minimisers.py`, `scripts/stage297_reduced_all_minimisers.py` |
| Forced-positivity face reduction for the minimum pass | `scripts/stage295_reduced_minimum.py` |
| `eps`-tolerant rational LP, row generation, exhaustive support enumeration | `scripts/stage296_common.py` |
| Exact `O(\|V\|)` closed form for `Theta_1` | `scripts/stage298_common.py` (`edge_aggregates`, `theta1_edge`, `theta1_all_edges`) |
| Transported screens and the certified `Theta_k` bracket | `scripts/stage301_common.py` |

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

## Licence

MIT, see `LICENSE`. Copyright (c) 2026 Haotian Tang.
