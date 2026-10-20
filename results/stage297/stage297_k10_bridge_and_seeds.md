# Stage 297 — the K = 10 support-uniqueness bridge, and DecoDiPhy seed sensitivity

Protocol `results/stage297/stage297_protocol.json`, registered 2026-10-18T16:33:11.025030+00:00, with addenda `stage297_protocol_addendum_main_k10.json` and `stage297_protocol_addendum_k10_extra.json` and `stage297_protocol_addendum_k10_large.json` (each registered before the batch it covers).  Every run is retained, including timeouts and deaths.

### Scope actually covered

- K = 10: **99 distinct bank cases** of the 213 exact K = 10 cases — a stratified subsample (first 3 per distinct |R| in `main_k10`, 2 more in `main_k10_extra`, plus the regime-restricted large-panel batch), not the whole family.
- The `k10_rand_large` contrast batch was **stopped after 10 of 55 cases**; the random arm at |R| >= 81 is therefore covered by whatever the five-regime batches completed plus those 10, not by a full sweep.
- The registered `pilot_k10_threshold` batch was not executed; its six cases are the first per-size picks of `main_k10`, which covers them.
- K = 12: one case attempted, not measured (section 4).
- Seed sweep: 132 cases, one per (family, K, |R|) cell, full panels only.

## 1  What was actually missing at K >= 10

Stage 292 could not reach K = 10 because the **frozen minimum pass** its all-minimisers routine calls (`stage43_partial_fpt.solve_partial`) returned `unknown` after 3,600 s on the smallest K = 10 bank case (`results/stage292_k10_ceiling.txt`), and all eight K = 10 jobs of its stride-3 batch hit the per-case wall (`results/stage292_bridge_K10_12.jsonl`, 8/8 timeout).

The forced-positivity face reduction, however, was **already present** in `stage292_all_minimisers.backbone_patterns` — the `full_mask` switch, cross-checked there on 500/500 fresh instances.  So the combination this stage needed was not "add the reduction to the face enumeration" but "replace the frozen minimum pass".  `scripts/stage297_reduced_all_minimisers.py` offers three interchangeable modes:

| mode | pass 1 (minimum) | pass 2 (all minimisers) |
|---|---|---|
| `frozen` | frozen `solve_partial` | reduced faces at k\* — the Stage292 routine |
| `reduced_two_pass` | Stage295 `solve_partial_reduced` | reduced faces at k\* |
| `combined` | — | one reduced enumeration yields the minimum **and** every minimiser at it |

`combined` is the operative one: a face is LP-tested when its cover count is at most the best count established so far, the incumbent is lowered by any strictly feasible smaller face, and every feasible face at the final incumbent is retained. Because `best` only decreases and every face with count <= best is tested, the result does not depend on visiting order, and the dominant cost (deficit sets x faces x exact rational LP) is paid once instead of twice.

## 2  Validation — three tracks, zero disagreements

| track | what it compares | rows | identical minimiser sets |
|---|---|---|---|
| A replay | every probe point of `stage292_validate_replay.jsonl` (130 cases, K in {2,3}), removal orders regenerated with the frozen `removal_order` | 5339 | combined 5339/5339; reduced_two_pass 5339/5339; frozen 5339/5339 |
| B fresh | the 500 instances of `stage292_validate_fresh.jsonl`, regenerated from `stage292_validate.random_instance` | 500 | combined 500/500; reduced_two_pass 500/500; frozen 500/500; combined_full_mask 500/500 |
| C high K | the full-panel baseline of all 676 Stage292 bridge cases, K in {2,3,5,7} | 676 | combined 676/676; reduced_two_pass 676/676; frozen 676/676 |

Additional checks: track A vs the value Stage292 recorded: 5209/5209; track B vs the independent exhaustive C(|E|,k) enumeration: 500/500; track C vs the value Stage292 recorded: 676/676.

The Stage295 reduced minimum equals the frozen minimum in **6515/6515** comparisons across the three tracks.

With the reduction **disabled** (`full_mask=True`) the Stage297 enumeration reproduces the unreduced face count of the Stage292 enumeration exactly in 68/68 checks, with identical pattern counts (68/68) and identical pattern sets between reduced and full (68/68) — the correctness signal Stage295 relied on. Median faces: 168 full vs 72 reduced (median reduction factor 2.00, max 91.8).

Median speed-up of `combined` over the Stage292 two-pass routine: track A 1.44x, track B 2.27x, track C 2.22x.

## 3  K = 10

### 3.1  How many K = 10 cases admit the question at all

`t_unknown` is defined only where the full-panel answer is already the unique true support (Stage292 convention: cases with `baseline_unique_true = false` are excluded from every threshold table).  At K = 10 that exclusion is no longer negligible:

| | cases measured | minimum = K | unique true support at the full panel |
|---|---|---|---|
| Stage292, K in {2,3,5,7} | 676 | 674 | 663 |
| Stage297, K = 10 | 99 | 67 | 55 |

Exact minima found at K = 10: `{"K10|min10": 67, "K10|minNone": 12, "K10|min9": 16, "K10|min8": 3, "K10|min7": 1}`.  In the cases where the minimum falls below the generating K the data are explained exactly by fewer than K sources, so the true support is not the minimum-source answer even before any reference is erased.

Number of distinct minimisers at the full panel: `{"1": 72, "0": 12, "3": 6, "5": 4, "7": 4, "9": 1}`.

Mode cross-check at K = 10 (no validation set reaches this K): 15/16 cases agree on the minimum and 15/16 on the complete minimiser set between `combined` and `reduced_two_pass`; median 141.3 s vs 287.8 s.

### 3.2  Threshold table

`t_unknown` = the largest number of references that can be erased, in the order of the regime, with the set of ALL minimum-size supports still exactly the one true support.  Search convention imported verbatim from `stage292_bridge.find_threshold` (exponential bracketing, then bisection, both bracket endpoints evaluated, MIN_REMAINING = 3).

**Stage297, K = 10**

| series | n | median | q25 | q75 | min | max |
|---|---|---|---|---|---|---|
| t_unknown K10 / adversarial_nearest | 52 | 2.0 | 1 | 3 | 0 | 9 |
| t_unknown K10 / random | 54 | 7.0 | 3 | 26 | 0 | 85 |
| t_unknown K10 / whole_clade | 59 | 8 | 4 | 19 | 0 | 117 |
| fraction K10 / adversarial_nearest | 52 | 0.014 | 0.004 | 0.028 | 0.000 | 0.125 |
| fraction K10 / random | 54 | 0.097 | 0.038 | 0.215 | 0.000 | 0.423 |
| fraction K10 / whole_clade | 59 | 0.100 | 0.054 | 0.196 | 0.000 | 0.331 |

**Stage292, K in {2,3,5,7}**

| series | n | median | q25 | q75 | min | max |
|---|---|---|---|---|---|---|
| t_unknown K2 / adversarial_nearest | 120 | 2.0 | 1 | 4 | 0 | 115 |
| t_unknown K2 / random | 240 | 74.0 | 30 | 122 | 1 | 337 |
| t_unknown K2 / whole_clade | 240 | 37.5 | 17 | 74 | 0 | 328 |
| t_unknown K3 / adversarial_nearest | 218 | 2.0 | 1 | 4 | 0 | 25 |
| t_unknown K3 / random | 436 | 33.0 | 19 | 70 | 0 | 300 |
| t_unknown K3 / whole_clade | 436 | 20.0 | 8 | 42 | 0 | 262 |
| t_unknown K5 / adversarial_nearest | 217 | 2 | 1 | 4 | 0 | 27 |
| t_unknown K5 / random | 423 | 25 | 10 | 52 | 0 | 271 |
| t_unknown K5 / whole_clade | 434 | 12.0 | 5 | 31 | 0 | 212 |
| t_unknown K7 / adversarial_nearest | 106 | 2.0 | 1 | 4 | 0 | 9 |
| t_unknown K7 / random | 177 | 32 | 10 | 59 | 0 | 264 |
| t_unknown K7 / whole_clade | 210 | 16.0 | 5 | 37 | 0 | 186 |
| fraction K2 / adversarial_nearest | 120 | 0.020 | 0.008 | 0.043 | 0.000 | 0.676 |
| fraction K2 / random | 240 | 0.647 | 0.427 | 0.767 | 0.022 | 0.979 |
| fraction K2 / whole_clade | 240 | 0.375 | 0.208 | 0.536 | 0.000 | 0.946 |
| fraction K3 / adversarial_nearest | 218 | 0.023 | 0.010 | 0.052 | 0.000 | 0.519 |
| fraction K3 / random | 436 | 0.488 | 0.310 | 0.660 | 0.000 | 0.941 |
| fraction K3 / whole_clade | 436 | 0.290 | 0.137 | 0.483 | 0.000 | 0.954 |
| fraction K5 / adversarial_nearest | 217 | 0.023 | 0.010 | 0.054 | 0.000 | 0.372 |
| fraction K5 / random | 423 | 0.353 | 0.207 | 0.500 | 0.000 | 0.849 |
| fraction K5 / whole_clade | 434 | 0.188 | 0.089 | 0.290 | 0.000 | 0.728 |
| fraction K7 / adversarial_nearest | 106 | 0.014 | 0.006 | 0.031 | 0.000 | 0.220 |
| fraction K7 / random | 177 | 0.262 | 0.137 | 0.393 | 0.000 | 0.742 |
| fraction K7 / whole_clade | 210 | 0.138 | 0.059 | 0.204 | 0.000 | 0.522 |

### 3.3  Flatness in panel size

| regime | Spearman(t_unknown, panel size), Stage292 K in {2,3,5,7} | Stage297 K = 10 | pooled K in {2,3,5,7,10} |
|---|---|---|---|
| random | 0.832 (n=1276) | 0.710 (n=54) | 0.820 (n=1330) |
| whole_clade | 0.647 (n=1320) | 0.819 (n=59) | 0.651 (n=1379) |
| adversarial_nearest | 0.095 (n=661) | 0.203 (n=52) | 0.099 (n=713) |

### 3.4  Verdict

**The flatness claim holds at K = 10.**

Across 52 adversarial measurements spanning |R| = 20 to 353 the median nearest-to-source threshold is 2.0 references (q25 1, q75 3, max 9) and its rank correlation with reference-set size is 0.203, against 0.710 for random loss (median 7.0, max 85) and 0.819 for whole-clade loss (median 8, max 117).  Stage292 measured 0.095 / 0.832 / 0.647 on K in {2,3,5,7}.  The median adversarial threshold of 2.0 references at K = 10 is the same median Stage292 reports at every lower K, and the fraction of the panel erasable before the minimum-source answer stops being unique is 1.4% adversarial against 9.7% random and 10.0% whole-clade.

Two things do change at K = 10, and both are about the panel being small relative to the number of sources rather than about the regimes.  First, all three thresholds compress: the random and whole-clade fractions fall from Stage292's 42% and 22% of the panel to 9.7% and 10.0%, and the adversarial fraction from 2.1% to 1.4%.  Second, the question is answerable for a smaller share of cases, because the exact minimum often falls below the generating K (section 3.1).

18 (case, regime, seed) measurements were produced twice because a case that timed out in the five-regime batch was re-measured in a regime-restricted one; 18 of 18 agree, and repeats are pooled once.

## 4  The ceiling, with counters

Per-case failures in the K >= 10 batches: 41 (`{"timeout": 41}`).

| case | \|R\| | stop | wall s |
|---|---|---|---|
| bio_tilapia-Ciezarek__tilapia-Ciezarek.tre.rerooted_ | 81 | timeout | 7200 |
| bio_tilapia-Ciezarek__tilapia-Ciezarek.tre.rerooted_ | 81 | timeout | 7200 |
| bio_1kp__1kp-concat-fig2.tre_k10_s10_noise0_exact | 93 | timeout | 7200 |
| bio_1kp__1kp-concat-fig2.tre_k10_s7_noise0_exact | 93 | timeout | 7200 |
| bio_pancrustacean-Bernot__pancrustacean-Bernot.tre.r | 95 | timeout | 7200 |
| bio_pancrustacean-Bernot__pancrustacean-Bernot.tre.r | 95 | timeout | 7200 |
| k10k_tree_15_k10_s1_exp_exact | 101 | timeout | 7200 |
| k10k_tree_15_k10_s2_exp_exact | 101 | timeout | 7200 |
| k10k_tree_15_k10_s3_exp_exact | 101 | timeout | 7200 |
| bio_beetles-Johnson__beetles-Johnson.tre_k10_s10_noi | 117 | timeout | 7200 |
| bio_beetles-Johnson__beetles-Johnson.tre_k10_s8_nois | 117 | timeout | 7200 |
| k10k_tree_11_k10_s3_uniform_exact | 122 | timeout | 7200 |
| k10k_tree_11_k10_s4_uniform_exact | 122 | timeout | 7200 |
| k10k_tree_11_k10_s5_uniform_exact | 122 | timeout | 7200 |
| bio_fish-Troyer__fish-Troyer.tre.rerooted_k10_s7_noi | 142 | timeout | 7200 |
| bio_fish-Troyer__fish-Troyer.tre.rerooted_k10_s9_noi | 142 | timeout | 7200 |
| bio_hemipteroid-johnson__hemipteroid-johnson.tre.rer | 183 | timeout | 7200 |
| bio_hemipteroid-johnson__hemipteroid-johnson.tre.rer | 183 | timeout | 7200 |
| bio_hemipteroid-johnson__hemipteroid-johnson.tre.rer | 183 | timeout | 7200 |
| bio_fish_hughes__fish_hughes_1105examl.tre.rerooted_ | 295 | timeout | 7200 |
| bio_fish_hughes__fish_hughes_1105examl.tre.rerooted_ | 295 | timeout | 7200 |
| bio_fish_hughes__fish_hughes_1105examl.tre.rerooted_ | 295 | timeout | 7200 |
| bio_birds-stiller__birds-stiller-castlepro.tre_k10_s | 353 | timeout | 7200 |
| k10k_tree_8_k10_s1_exp_exact | 43 | timeout | 7200 |
| bio_tilapia-Ciezarek__tilapia-Ciezarek.tre.rerooted_ | 81 | timeout | 7200 |
| bio_1kp__1kp-concat-fig2.tre_k10_s6_noise0_exact | 93 | timeout | 7200 |
| bio_pancrustacean-Bernot__pancrustacean-Bernot.tre.r | 95 | timeout | 7200 |
| k10k_tree_15_k10_s5_uniform_exact | 101 | timeout | 7200 |
| bio_beetles-Johnson__beetles-Johnson.tre_k10_s2_nois | 117 | timeout | 7200 |
| bio_beetles-Johnson__beetles-Johnson.tre_k10_s6_nois | 117 | timeout | 7200 |

Probe outcomes: `{"adversarial_nearest|complete": 153, "adversarial_nearest|early_exit": 110, "adversarial_nearest|resource_limit": 13, "baseline|complete": 121, "baseline|resource_limit": 16, "random|complete": 408, "random|early_exit": 151, "random|oracle_limit": 7, "random|resource_limit": 12, "whole_clade|complete": 349, "whole_clade|early_exit": 168, "whole_clade|resource_limit": 5}`.

**K = 12** (`k16k_n100_SRS014343_s1_exact`, |R| = 15941, the single `exact_16k` case): status `resource_limit` after 3600 s — the honest ceiling, with counters rather than a quiet retreat. In that budget the routine enumerated 28,238,250 deficit sets and 1,086,980,096 reduced face masks, skipping 337,885,396,992 masks that forced positivity proves infeasible (311.8x), tested 17 faces by exact LP, found 1 feasible pattern(s) and never reached the endpoint oracle (0 calls).

## 5  DecoDiPhy seed sensitivity

924 runs over 132 cases (stratified over |R| within every (family, K) cell of `exact_10k`, `exact_biotrees` and `krepp_estimated_10k`), full panels, documented defaults, k chosen by the tool's own stopping rules; seeds [1142, 1, 7, 101, 2024, 31337].  Run status: `{"complete": 924}`.  Runs recovered from a committed `all_rounds.json` after a jplace-writer crash: 0.

| cell | cases | k identical across seeds | support identical across seeds | median mean-pairwise Jaccard | median W1 relative range | max W1 relative range |
|---|---|---|---|---|---|---|
| all | 132 | 125/132 | 123/132 | 1.000 | 0.000 | 1.287 |
| family exact_10k | 27 | 25/27 | 25/27 | 1.000 | 0.000 | 0.925 |
| family exact_biotrees | 60 | 58/60 | 57/60 | 1.000 | 0.000 | 1.287 |
| family krepp_estimated_10k | 45 | 42/45 | 41/45 | 1.000 | 0.000 | 0.001 |
| K 2 | 12 | 12/12 | 12/12 | 1.000 | 0.000 | 1.287 |
| K 3 | 34 | 32/34 | 32/34 | 1.000 | 0.000 | 0.925 |
| K 5 | 37 | 36/37 | 36/37 | 1.000 | 0.000 | 0.732 |
| K 7 | 12 | 12/12 | 12/12 | 1.000 | 0.000 | 0.570 |
| K 10 | 37 | 33/37 | 31/37 | 1.000 | 0.000 | 1.188 |

The 9 cases that move with the seed move little: the k range is {"1": 6, "0": 2, "3": 1} (in 2 of them k is fixed and only the support moves; the widest spread is k in [9, 10, 12]), the mean pairwise Jaccard of their supports is 0.759 to 0.981, and their W1 relative range is 0.000 to 1.188.  6 of them are K = 10 cases.

**Fixed-seed control.**  Arm `default_1142` is `stage293_abundance.dd_worker` imported unmodified (no `-s` flag, so the tool takes its documented default seed); arm `seeded_1142` is the one-flag copy at the same seed.  They agree in 132/132 cases on k, 132/132 on the returned support, 132/132 on the loss and 132/132 on W1. The unified diff of the two workers is in `results/stage297/stage297_dd_worker_diff.json`.

**Cross-harness comparison.**  The same tool, the same documented default seed, the same cases, two different harnesses: this stage's seed-1142 arm against the runs Stage289 recorded (`stage289_headtohead_runs.jsonl`, method `dd_default`, both repetitions).  Of 264 comparable comparisons, k agrees in 256, the support in 246 and the loss in 210.  Stage289's own two repetitions agree with each other in 1182/1182 on k and 1182/1182 on the support.

By family: `{"exact_10k": {"n": 54, "k_equal": 52, "support_equal": 50}, "exact_biotrees": {"n": 120, "k_equal": 120, "support_equal": 112}, "krepp_estimated_10k": {"n": 90, "k_equal": 84, "support_equal": 84}}`.

## 6  What none of this establishes

- The validation shows the Stage297 routine returns what the Stage292 routine returns wherever Stage292 completes, and what an independent exhaustive enumeration returns on small instances.  Exhaustiveness at K = 10 still rests on the frozen `deficit_sets` generator and the weighted deficit bound, which this stage imports and does not re-derive; no brute force is affordable at K = 10.
- Mode agreement at K = 10 is a cross-check between two routines that share pass 2, not an independent confirmation.
- The K = 10 sample is a stratified subsample of the 213 exact K = 10 bank cases, not the whole family, and the cases whose minimum falls below K are excluded from the threshold tables by the Stage292 convention — they are not evidence about erasure, but they do mean the K = 10 threshold statistics rest on fewer cases than the K <= 7 ones.
- K = 12 is not measured.  The single `exact_16k` case did not finish; the stage reports its counters, not a threshold.
- The seed sweep covers full panels only.  It says nothing about how seed choice interacts with reference removal, which is where Stages 288/294 found the tool failing outright.
- The cross-harness disagreement is characterised, not explained.  Both harnesses are internally reproducible, so the difference lies in process-level details (Stage289 ran the tool as a subprocess, Stage293/297 call `cli.main()` in process); which detail is responsible is not identified here.

## 7  Files

- `scripts/stage297_reduced_all_minimisers.py` — the combined routine
- `scripts/stage297_common.py` — bank access, forked-child pool
- `scripts/stage297_validate.py` — validation tracks A, B, C
- `scripts/stage297_bridge.py` — K = 10 / K = 12 threshold measurement
- `scripts/stage297_seeds.py` — DecoDiPhy seed sweep
- `scripts/stage297_analyze.py` — item 1 analysis
- `scripts/stage297_dd_analyze.py` — item 2 analysis
- `scripts/stage297_cross_harness.py` — Stage289 cross-harness comparison
- `scripts/stage297_protocol.py` — protocol registration
- `scripts/stage297_write_document.py` — this document
- `results/stage297/stage297_protocol.json` — registered protocol
- `results/stage297/stage297_protocol_addendum_main_k10.json` — main K = 10 batch
- `results/stage297/stage297_protocol_addendum_k10_extra.json` — K = 10 extension
- `results/stage297/stage297_validate_replay.jsonl` — track A records
- `results/stage297/stage297_validate_fresh.jsonl` — track B records
- `results/stage297/stage297_validate_baseline.jsonl` — track C records
- `results/stage297/stage297_bridge_pilot_k10_baseline.jsonl` — K = 10 pilot baselines
- `results/stage297/stage297_bridge_main_k10.jsonl` — K = 10 main batch
- `results/stage297/stage297_bridge_main_k10_extra.jsonl` — K = 10 extension batch
- `results/stage297/stage297_bridge_k10_adv_large.jsonl` — adversarial-only large-panel batch
- `results/stage297/stage297_bridge_k10_rand_large.jsonl` — random-only large-panel batch, stopped early
- `results/stage297/stage297_bridge_k12_baseline.jsonl` — K = 12 attempt
- `results/stage297/stage297_protocol_addendum_k10_large.json` — regime-split large-panel batches
- `results/stage297/stage297_main_k10_manifest.jsonl` — skip manifest: the 60 main_k10 case names
- `results/stage297/stage297_adv_done_manifest.jsonl` — skip manifest: cases whose adversarial threshold was already found
- `results/stage297/stage297_bridge_v1_as_registered.py` — the bridge script as registered in the protocol
- `results/stage297/stage297_bridge_v2_main_k10.py` — bridge after main_k10
- `results/stage297/stage297_bridge_v3_k10_extra.py` — bridge after the extension
- `results/stage297/stage297_dd_seeds_main.jsonl` — seed sweep runs
- `results/stage297/stage297_analysis.json` — item 1 analysis output
- `results/stage297/stage297_dd_analysis.json` — item 2 analysis output
- `results/stage297/stage297_dd_vs_stage289.json` — cross-harness comparison
- `results/stage297/stage297_dd_worker_diff.json` — worker diff

