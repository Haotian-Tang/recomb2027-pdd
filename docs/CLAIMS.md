# Claim-to-record index

This page maps numerical claims in the paper to their producing scripts and
archived records. All file paths below are relative to the repository root.
See the [README](../README.md) for reproduction commands and the
[bulk archive instructions](../README.md#bulk-data-archive) for records marked **(bulk)**.

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
| Figure 2, "What incomplete references cost" - (a) median fall in `Theta_1` against the share of references removed, exact arm, by regime; (b) the same draws on real ANI estimates of the same trees and the identical leaf sets; (c) references erasable before the minimum-size explanation stops being unique, against panel size, pooled over `K` = 2-10. Caption counts: 592 draws per regime and level, 296 for nearest to sources; (c) pools 713 adversarial and 1,330 random thresholds | `figures/make_fig_incompleteness.py` (redraws it from the payload; needs numpy and matplotlib) | payload `figures/data/fig_incompleteness_data.json`, derived from `results/stage299/stage299_regime_ratios.json` (`fall` block) and `results/stage297/stage297_analysis.json` (`threshold_vs_panel` block); output `figures/incompleteness.{pdf,png}` |

### Uniqueness

| Paper claim | Script | Record |
|---|---|---|
| A median of two references nearest the sources can be erased before the generating support stops being the unique smallest explanation, at every `K` measured | `stage290_erasure_core.py`, `stage292_all_minimisers.py`, `stage297_reduced_all_minimisers.py` | `results/stage292_analysis.json`; `results/stage297/stage297_analysis.json` `/k10/thresholds_stage297/t_unknown`; `results/stage302/stage302_failure_mode.json` `/by_regime/*/median_t_unknown`; `results/stage292_support_uniqueness.md` |
| The erasure that breaks it leaves a same-size competing support in 671 of 713 adversarial measurements and lowers the minimum below `K` in the other 42 | `stage302_failure_mode.py` | `results/stage302/stage302_failure_mode.json` `/by_regime/adversarial_nearest/first_failure_mode` |
| The all-minimisers routine agrees with brute force on 5,839 of 5,839 comparisons; of those, 5,339 replay the recorded minimiser sets and the remaining 500 are fresh instances checked against an independent binomial enumeration | `stage297_validate.py` | `results/stage297/stage297_analysis.json` `/validation/track_A/rows` (5,339 replay) and `/validation/track_B/rows` (500 fresh, `brute_subsets_tested` 148,233); raw `results/stage297/stage297_validate_replay.jsonl` **(bulk)** |
| One of them a case whose complete minimiser set contains 4,443 supports | `stage292_validate.py` | `results/stage292_analysis.json` `/validation/track_A/n_supports_distribution/4443` |
| Pooled over `K`, the adversarial threshold's rank correlation with panel size is 0.10 over 713 measurements spanning 20 to 361 references, against 0.82 over 1,330 random ones | `stage302_failure_mode.py` | `results/stage302/stage302_failure_mode.json` `/by_regime/{adversarial_nearest,random}/{spearman_t_vs_panel_size, n_thresholds, panel_size_min, panel_size_max}`; the same values are recomputable from `figures/data/fig_incompleteness_data.json` `/threshold_vs_panel` |
| Pooled over `K`, random loss tolerates a median of 32.5 references and whole clades 18 | `stage302_failure_mode.py` | `results/stage302/stage302_failure_mode.json` `/by_regime/{random,whole_clade}/median_t_unknown` |
| On planted panels to 47,151 leaves the random threshold grows linearly in panel size | `stage292_panels.py`, `stage292_analysis.py` | `results/stage292_analysis.json` `/panels`; `results/stage292_panels.jsonl` |
| Unfinished threshold searches cannot explain the regime gap: bounding every unrecorded cell of a balanced grid over the 718 baseline-unique cases by its admissible extreme, the median of the complete grid lies in [2, 2] for nearest-to-sources, [32, 36] for random and [17, 19] for whole clades | `stage305_missing_thresholds.py` | `results/stage305/stage305_missing_thresholds.json` `/by_regime` |
| Supplement: at `K = 10`, 67 of 99 have minimum = `K` and 55 a unique true support, against 663 of 676 at `K <= 7` | `stage297_bridge.py`, `stage297_analyze.py`; `stage292_bridge.py` | `results/stage297/stage297_analysis.json` `/k10/{cases, baseline_minimum_equals_K, baseline_unique_true}`; `results/stage292_analysis.json` `/bridge/{cases, baseline_unique_true}` |
| Erasable fraction 2.1% pooled over `K <= 7` and 1.4% at `K = 10`, against 42.2% and 9.7% at random | `stage297_analyze.py` | `results/stage297/stage297_analysis.json` `/k10/thresholds_stage297/fraction_erasable` |

### The tolerance

| Paper claim | Script | Record |
|---|---|---|
| With no reference removed the minimum first drops below the exact count at a median tolerance 0.0147 of the case's median observed distance, over 43 cases | `stage296_tolerance.py`, `stage296_analyze.py` | `results/stage296/stage296_analysis.json` `/threshold_distribution`, `/full_panel[*]/first_drop_eps_rel` |
| It tracks the smallest constituent far more closely than the panel size (`rho = 0.665` against `0.187`) | `stage296_analyze.py` | `results/stage296/stage296_analysis.json` `/threshold_scaling/vs_min_p/spearman`, `/threshold_scaling/vs_N/spearman` |
| Real ANI distances deviate by 0.34 of a median observed distance at the worst leaf, 7.3 to 14.5 times that tolerance in 15 of 15 cases | `stage296_tolerance.py` phase `calibrate` | `results/stage296/stage296_calibration.json`; `results/stage296/stage296_analysis.json` |
| At a tolerance holding that error the model's own minimum is one source in all 15, against three at tolerance zero in 14 of them | `stage296_analyze.py` | `results/stage296/stage296_analysis.json`; `results/stage296/stage296_tolerance_titration.md` |
| Figure 3 (tolerance axis) | `figures/make_fig_tolerance_axis.py` | `results/stage296/stage296_payload.json` (every plotted value, read directly by the script); vector `figures/tolerance_axis.pdf` |

### Real distances

| Paper claim | Script | Record |
|---|---|---|
| The worst-leaf residual of the true configuration exceeds `Theta_1` in 296 of 296 ANI cases, median ratio 2.9 | `stage298_sweep.py`, `stage298_residuals.py`, `stage298_analyze.py`; closed form in `stage298_common.py` | `results/stage298/stage298_analysis.json` `/by_family/krepp_estimated_10k/residual_max_over_theta1`; `results/stage298/stage298_residuals.json` — **reproduction B** |
| Not an artefact of the worst-leaf norm: under the sum-of-squares analogue `Lambda_1` the generating configuration fails to beat the best single source in 294 of 296 | `stage300_evaluate.py`, `stage300_analyze.py` | `results/stage300/stage300_report_payload.json` `/analysis/by_collection/krepp_ani/arms/dd_default/n_truth_beats_one_source_in_ssr` (= 2 of 296) |
| and in 50 of 60 placement cases | `stage298_placement.py`, `stage298_analyze.py` | `results/stage298/stage298_analysis.json` `/by_family` placement block |
| At 75% loss the median fall in `Theta_1` is 0.34 at random against 0.63 whole clade and 0.61 nearest (1.9x and 1.8x), where the paired exact arm on identical leaf sets gives 2.9x and 2.6x (Figure 2) | `stage299_small.py`, `stage299_exact_baseline.py`, `stage299_analyze.py` | `results/stage299/stage299_analysis.json` `/by_arm_level_regime`; `results/stage299/stage299_regime_ratios.json`; per-draw table `results/stage299/stage299_per_draw.csv` **(bulk)** |
| On the 10,575-reference backbone, discarding 7,931 references at random costs 1.8% of `Theta_1` and comparably many as whole clades costs 90% (89.6%) | `stage299_td.py`, `stage299_analyze.py` | `results/stage299/stage299_analysis.json` `/by_arm_level_regime` travellers rows; `results/stage299/stage299_real_incompleteness.md`; raw `results/stage299/stage299_td_runs.jsonl` **(bulk)** |
| No terminal edge attains `Theta_1` on any of the 29 full panels, whose attaining sets are complete; 92 of the 725 draws carry one among their tied optima, never as the canonical representative | `stage298_td.py`, `stage299_td.py`, `stage304_attaining_sets.py` | `results/stage304/stage304_attaining_sets.json` `/representative_is_terminal` (0), `/attaining_set_contains_terminal` (92), `/terminal_undecided` (0) with the `/terminal_classification` breakdown 633 + 52 + 40 + 29 that makes the count exact |
| Nearest-to-attainer removal changes the canonical attaining edge in 22 of 29 stool draws at 25% loss, ten of which still share an optimal edge with the full panel, so the move is forced in 12; whole-clade removal of three quarters forces it in all 58, a median 22 edges away | `stage299_td.py`, `stage304_attaining_sets.py` | `results/stage304/stage304_attaining_sets.json` `/by_regime_level/{nearest_theta1\|0.25, clade\|0.75}`; the nine truncated comparisons were resolved by recomputing both attaining sets in full, 2 of 9 disjoint |
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
| Against DecoDiPhy 1.5.0 given the true `k`, on the 675 exactly generated profiles where the exact method returned an explanation it is a median 36 times cheaper end to end and 54 times on solver time alone; it timed out on 217 of the 892 such profiles. On the 282 noisy profiles it rejects exact explainability a median 93 times faster than the heuristic fits them, a different task | `stage303_timing.py` | `results/stage303/stage303_timing.json` `/ratios/*/by_stratum`; raw `results/stage289/stage289_headtohead_runs.jsonl` **(bulk)** |
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

