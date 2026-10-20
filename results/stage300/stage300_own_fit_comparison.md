# Stage 300 - DecoDiPhy's own fit, measured against Theta_1

Protocol `results/stage300/stage300_protocol.json`, registered 2026-10-19T08:21:54+00:00.  Inputs, code and the published package are hashed there.  Every number in this document is read from the stage-300 records by `scripts/stage300_document.py`.

## 1. What is measured

For a configuration - the true generating one, the one DecoDiPhy returns, or a single source - the prediction at an observed reference leaf r is

```
pred_r = y + sum_i p_i * dist_T(r, q_i),   dist(r,q) = min(d(r,u) + xL, d(r,v) + (1-x)L)
```

and two fits are reported for it:

- **R** = `max_r |pred_r - d_r|`, the L-infinity residual, computed in exact rational arithmetic (the configuration's floats enter by their exact binary value).  Every `R >= Theta_1` decision below is therefore exact, not a tolerance test.
- **SSR** = `sum_r (pred_r - d_r)^2` in float64 - DecoDiPhy's own objective, `np.sum((d - (DA @ p + CAl @ x + y)) ** 2)` in `decodiphy/optimize.py`.

The two single-source baselines:

- **Theta_1** - the exact smallest tolerance at which ONE source explains the observations, **imported** from stage 298 (`stage298_analysis.small_cases` and `stage298_td_summary.rows`) and never recomputed differently here.
- **Lambda_1** - the exhaustive minimum of SSR over every original edge, every x in [0,1] and every y >= 0 with p = 1: the same feasible set DecoDiPhy's QP uses at k = 1, enumerated instead of hill-climbed.  It is the L2 analogue of Theta_1 and exists so the comparison can be made in the norm the tool optimises.

DecoDiPhy reports a per-placement pendant vector; the model's additive constant is its sum, because `compute_individual_ys` splits the scalar pendant proportionally and the parts sum back to it by construction.  Abundances, positions and pendant are taken verbatim from the tool's output and mapped onto bank edge indices by the frozen stage-289 bridge.  Nothing is renormalised.

## 2. Where the DecoDiPhy answers come from

| collection | cases | DecoDiPhy runs harvested | from | runs newly executed | runs not usable |
|---|---|---|---|---|---|
| krepp_ani | 296 | 2683 | stage 289, stage 294, stage 297 | 0 | 52 |
| krepp_placement | 60 | 0 | - | 120 | 0 |
| travelers | 29 | 0 | - | 6 | 6 |

Harvested runs are full-panel runs already on record from stages 289/293/294/297, each row carrying its file, line, arm and repetition or seed; nothing was re-run.  The krepp ANI collection reaches 296/296 cases only by combining stage 289 (290 cases complete) with stage 294 (296, whose harness recovers the configuration after the `save_jplace` KeyError).  The 60 placement cases and the 29 stool samples had no prior runs and were executed here with the same worker (`stage293_abundance.dd_worker`, imported unmodified): documented defaults, `--method hill --optimizer osqp`, the only operable combination in 1.5.0.

## 3. The decisive count, in L-infinity

| collection | arm | cases with an answer | k returned (median) | R_dd/Theta_1 median | R_dd >= Theta_1 | truth/Theta_1 median | truth >= Theta_1 | DecoDiPhy beats truth (L-inf) |
|---|---|---|---|---|---|---|---|---|
| krepp_ani | dd_default | 296/296 | 7 | 0.905 | 121/296 | 2.865 | 296/296 | 282/296 |
| krepp_ani | dd_true_k | 292/296 | 5 | 0.916 | 125/292 | 2.842 | 292/292 | 282/292 |
| krepp_ani | dd_resolution | 289/296 | 6 | 0.881 | 115/289 | 2.888 | 289/289 | 275/289 |
| krepp_placement | dd_default | 60/60 | 20 | 0.439 | 0/60 | 1.49 | 50/60 | 60/60 |
| krepp_placement | dd_true_k | 60/60 | 5 | 0.589 | 3/60 | 1.49 | 50/60 | 60/60 |

Ordering of the three fits per case, best first (`1` = one source at Theta_1, `D` = DecoDiPhy's returned mixture, `T` = the true configuration):

- krepp ANI (dd_default): `D1T` 175, `1DT` 107, `1TD` 14
- krepp placement (dd_default): `D1T` 50, `DT1` 10

- On the 296 krepp ANI cases the tool's own mixture fits **no better in L-infinity than the best single source in 121 of 296 cases** (40.9%); median R_dd/Theta_1 = 0.905, upper quartile 1.103, maximum 4.129.
- On the 60 placement cases it fits better than one source everywhere (0/60 at or above Theta_1, median ratio 0.439).
- The true configuration is at or above Theta_1 in 296/296 krepp ANI and 50/60 placement cases - stage 298 reproduced exactly - and DecoDiPhy's answer fits better than the truth in 282/296 and 60/60.  Fitting the estimator's error is easier than reproducing the generator, and the tool does it.
- Handing the tool the true k does not change the verdict: 125/292 krepp ANI cases at or above Theta_1 with `dd_true_k`, median ratio 0.916.
- The decision is reproducible across replicates: identical in 296 of 296 krepp ANI cases over the 1-9 harvested runs each (k itself differs between replicates in 8 cases).
- Median k returned: 7 on krepp ANI (true K in {3, 5, 10}) and 20 on placement (true K = 5), i.e. the default stopping rule inflates the count on placement distances.

The verdict is flat in the true number of sources and strongly dependent on how many references the panel carries:

| collection | stratum | cases | R_dd >= Theta_1 | R_dd/Theta_1 median | k returned (median) |
|---|---|---|---|---|---|
| krepp_ani | true K = 3 | 100 | 38 (38%) | 0.919 | 7 |
| krepp_ani | true K = 5 | 100 | 44 (44%) | 0.868 | 7 |
| krepp_ani | true K = 10 | 96 | 39 (41%) | 0.865 | 6 |
| krepp_ani | panel, lower third (20-43 references) | 98 | 14 (14%) | 0.442 | 12 |
| krepp_ani | panel, middle third (43-53 references) | 98 | 55 (56%) | 1.021 | 5 |
| krepp_ani | panel, upper third (55-128 references) | 100 | 52 (52%) | 1.032 | 4 |
| krepp_placement | true K = 3 | 20 | 0 (0%) | 0.463 | 22 |
| krepp_placement | true K = 5 | 20 | 0 (0%) | 0.438 | 20 |
| krepp_placement | true K = 10 | 20 | 0 (0%) | 0.445 | 20 |
| krepp_placement | panel, lower third (24-29 references) | 20 | 0 (0%) | 0.279 | 17 |
| krepp_placement | panel, middle third (31-101 references) | 20 | 0 (0%) | 0.465 | 21 |
| krepp_placement | panel, upper third (106-108 references) | 20 | 0 (0%) | 0.627 | 26 |

With 20-43 references the mixture usually does beat one source in L-infinity; from about 43 references upwards it fails to in more than half the cases.  The median count it returns falls the other way (12, 5, 4 sources over the three panel strata), so the extra sources buy worst-leaf accuracy only while the panel is small enough to absorb them.

## 4. The same comparison in the norm DecoDiPhy optimises

**DecoDiPhy does not optimise L-infinity.**  It minimises SSR, and in SSR it wins.  This has to be read in the same breath as section 3.

| collection | arm | SSR(DecoDiPhy)/Lambda_1 median | beats one source in SSR | beats truth in SSR | truth beats one source in SSR | recomputed SSR vs its reported loss (median rel.) |
|---|---|---|---|---|---|---|
| krepp_ani | dd_default | 0.4 | 280/296 | 282/296 | 2/296 | 4.7e-06 |
| krepp_ani | dd_true_k | 0.456 | 281/292 | 282/292 | 2/292 | 1.5e-15 |
| krepp_placement | dd_default | 0.097 | 60/60 | 60/60 | 8/60 | 0.01 |
| krepp_placement | dd_true_k | 0.215 | 60/60 | 60/60 | 8/60 | 1.1e-15 |

- Against the best single source in its own norm it wins in 280/296 krepp ANI and 60/60 placement cases, with median SSR 0.4 and 0.097 of Lambda_1.
- So the conclusion in section 3 IS norm-dependent, and the honest statement is: the extra sources buy a better least-squares fit and do not buy a smaller worst-leaf error.  Theta_1's meaning - a tolerance that must contain the error at EVERY leaf - is intrinsically an L-infinity statement, which is why the paper's claim is made there; a reader who cares about least squares should read this table instead of section 3.
- The tool also beats the TRUE configuration in SSR in 282/296 krepp ANI cases, while the truth beats one source in SSR in only 2/296.  On real estimated distances the generator is not the least-squares optimum either.
- In 16 of 296 krepp ANI cases the returned mixture is worse in SSR than the best single source (ratio up to 13.555), i.e. its hill climbing can return an answer that loses to one source under its own objective.

### 4.1 Re-solving its own programme at its own support

Second fairness arm, in case the comparison is an artefact of what the tool writes rather than of what it computes: `decodiphy.optimize.solve_osqp_once` is called on the tool's own matrices at the support it returned, and the resulting configuration is evaluated by the same exact code path.

| collection | arm | cases re-solved | R(re-solved) >= Theta_1 | ratio median | returned-configuration verdict |
|---|---|---|---|---|---|
| krepp_ani | dd_default | 290/296 | 116/290 | 0.883 | 121/296 |
| krepp_ani | dd_true_k | 286/292 | 120/286 | 0.911 | 125/292 |
| krepp_placement | dd_default | 60/60 | 0/60 | 0.439 | 0/60 |
| krepp_placement | dd_true_k | 60/60 | 3/60 | 0.589 | 3/60 |

The verdict barely moves: 116 of 290 re-solved krepp ANI cases are still at or above Theta_1 (median ratio 0.883, against 0.905 for the configuration as returned).  The L-infinity statement is about the model and the support, not about output formatting.

## 5. The 29 travellers'-diarrhoea stool samples: a measured cost, not a missing cell

Nothing is known about these samples - no ground truth, no error model, no admissible floor - so the intended statement was a comparison of fit quality in a stated norm plus the number of constituents each answer carries.  Half of that is here: Theta_1 and Lambda_1 are measured for all 29 samples and need no DecoDiPhy run.  The other half is not: under this stage's invocation the tool returned no answer on any sample within the registered wall.

**Wall.** 7200 s (2 h) of wall-clock per sample, registered in `stage300_travelers_cost.json` with its justification: about 9.93 GiB resident per worker (still rising) is what the shared machine allows at 5 concurrent workers, and 29 samples x 2 h / 5 workers is about 11.6 h - the largest per-sample wall that keeps the collection inside a single day of this stage's 40-core budget.

**Measured.** 6 of 29 samples attempted, **0 completed**.  Per-sample wall-clock at the moment each run was stopped: 7286-9922 s (2.0-2.8 h); resident set size then: 9.9-10.8 GiB (a lower bound on the peak); 12.9 core-hours consumed in total, with no returned configuration.  The remaining 23 samples were not attempted.

**Configuration.** `decodiphy -t tree.nwk -d distances.txt -o OUTDIR (documented defaults: --method hill, --quick 1, --radius 2, --seed 1142, --fix_k 0, --optimizer osqp, --min_p 0.01)`; worker `stage293_abundance.dd_worker, imported unmodified`; concurrency 1 (pilot) and 5 (batch), chosen from the measured resident set size on a shared 62-core machine with other users' jobs running; RLIMIT_AS 60 GB (pilot) and 30 GB (batch); neither was reached; OMP/OPENBLAS/MKL_NUM_THREADS=1.

**Where the time went, measured on one sample** (`results/stage300/stage300_tdprobe_td_sample_101_published.json`): parsing the backbone 0.10 s, `preprocess_input` 1.79 s, treeswift's leaf-to-leaf `distance_matrix` 40 s for 111,820,050 entries, and the enclosing `get_input_matrices` had not returned after 3777 s (63 min) when that probe was stopped.  The input stage, not the search, is where this invocation spent its time.  Unlike the krepp ANI cases of section 6.2, nothing is pruned here: every one of the 10575 backbone leaves carries a distance.

**What this is not.**  This is the cost of this stage's invocation on this machine, not a property of DecoDiPhy: the parent publication ran it on this same collection.  Its README documents --warm_start because the search "can be time extensive", i.e. for inputs of this size the intended workflow resumes a search across runs instead of driving one run to the low-abundance stop. No conclusion about the tool's behaviour on these samples is drawn here.

Theta_1 and Lambda_1 for all 29 samples, which is what can be reported without the tool (0 of 29 samples have a DecoDiPhy answer):

| sample | group | references | Theta_1 | Lambda_1 (SSR) | L-inf residual of the Lambda_1 configuration | DecoDiPhy answer |
|---|---|---|---|---|---|---|
| S10 | TD | 10575 | 31.635 | 8.57e+05 | 47.049 | censored: stopped at 7286 s (2.0 h) with 10.0 GiB resident, no answer returned |
| S101 | TD | 10575 | 20.848 | 99399.26 | 27.164 | censored: stopped at 9922 s (2.8 h) with 10.8 GiB resident, no answer returned |
| S125 | TD | 10575 | 20.511 | 46623.039 | 38.087 | censored: stopped at 7286 s (2.0 h) with 10.0 GiB resident, no answer returned |
| S147 | TD | 10575 | 14.343 | 57747.355 | 20.787 | censored: stopped at 7286 s (2.0 h) with 10.0 GiB resident, no answer returned |
| S156 | TD | 10575 | 55.527 | 6.51e+06 | 87.585 | censored: stopped at 7286 s (2.0 h) with 9.9 GiB resident, no answer returned |
| S160 | TD | 10575 | 22.99 | 8.51e+05 | 25.705 | censored: stopped at 7286 s (2.0 h) with 9.9 GiB resident, no answer returned |
| S3 | TD | 10575 | 24.34 | 1.26e+05 | 41.448 | not attempted (see section 5) |
| S45 | TD | 10575 | 33.718 | 5.39e+05 | 37.426 | not attempted (see section 5) |
| S50070 | TD | 10575 | 26.637 | 3.53e+05 | 38.762 | not attempted (see section 5) |
| S50076 | TD | 10575 | 29.053 | 1.21e+06 | 32.447 | not attempted (see section 5) |
| S50395 | TD | 10575 | 38.804 | 1.53e+06 | 48.826 | not attempted (see section 5) |
| S510 | HT | 10575 | 32.377 | 5.55e+05 | 42.122 | not attempted (see section 5) |
| S521 | HT | 10575 | 56.943 | 2.76e+06 | 100.235 | not attempted (see section 5) |
| S523 | HT | 10575 | 41.698 | 8.04e+05 | 65.954 | not attempted (see section 5) |
| S528 | HT | 10575 | 41.041 | 1.40e+06 | 67.69 | not attempted (see section 5) |
| S538 | HT | 10575 | 66.433 | 3.25e+06 | 96.943 | not attempted (see section 5) |
| S540 | HT | 10575 | 32.091 | 7.66e+05 | 43.847 | not attempted (see section 5) |
| S560 | HT | 10575 | 21.4 | 85042.935 | 39.366 | not attempted (see section 5) |
| S6128 | TD | 10575 | 12.11 | 24533.915 | 17.983 | not attempted (see section 5) |
| S6163 | TD | 10575 | 20.538 | 4.92e+05 | 21.239 | not attempted (see section 5) |
| S6165 | TD | 10575 | 40.58 | 4.02e+06 | 65.002 | not attempted (see section 5) |
| S6168 | TD | 10575 | 36.595 | 1.90e+06 | 44.221 | not attempted (see section 5) |
| S678 | TD | 10575 | 48.503 | 3.81e+06 | 58.778 | not attempted (see section 5) |
| S715 | TD | 10575 | 24.047 | 4.04e+05 | 25.34 | not attempted (see section 5) |
| S76 | TD | 10575 | 54.741 | 5.27e+06 | 84.837 | not attempted (see section 5) |
| S78 | TD | 10575 | 36.596 | 2.22e+06 | 46.482 | not attempted (see section 5) |
| S80129 | TD | 10575 | 36.597 | 1.19e+06 | 53.562 | not attempted (see section 5) |
| S80142 | TD | 10575 | 28.715 | 5.86e+05 | 34.701 | not attempted (see section 5) |
| S80152 | TD | 10575 | 44.567 | 3.61e+06 | 76.405 | not attempted (see section 5) |

All three quantities are in the samples' own distance units (DEPP placement distances on the 16S backbone, median observed distance of order 300), so only ratios are comparable across samples.  The L-infinity residual of the Lambda_1 configuration is listed next to Theta_1 to show the gap between the two single-source optima: the least-squares best single source is worse at its worst leaf than the L-infinity best single source, by construction, and the column quantifies by how much.

## 6. Two properties of the implementation that this comparison exposed

### 6.1 The configuration written out is not always the one whose loss is printed

Recomputing SSR from the returned configuration and comparing it with the loss DecoDiPhy reports for that run: the two differ by more than 1e-6 relative in 496 of 1738 krepp ANI runs and 59 of 120 placement runs, and the returned configuration is the worse of the two in 1451 and 89 runs respectively.

`decodiphy/search.py:hill_climbing` resets `min_obj` at the top of every sweep but never resets `min_p / min_x / min_y`, so when the final sweep improves on nothing, the abundances, positions and pendant written out are those of the last improving candidate of an earlier sweep while the anchors and the loss belong to the converged support.  `all_rounds.json` and the jplace carry that triple.

Worst case, `k10k_tree_15_k3_s1_exp_krepp_estimated` (`stage300_losscase_*.json`): the tool reports loss 3.60e-20 at k = 2, the configuration it writes has SSR 0.31 on the same leaves, and re-solving its own QP at its own printed support attains 8.13e-25.  The printed number is a real objective value - just not the one belonging to the printed parameters.

This is why section 4.1 exists: with the support re-solved the comparison is unchanged, so nothing in section 3 rests on this inconsistency.

### 6.2 The input stage silently deletes reference leaves that do carry a distance

`decodiphy/tree_utils.py:preprocess_input` removes leaves without a distance while iterating over `tree_obj.traverse_leaves()`, calling `remove_child` and `contract` inside that loop; that mutation-during-traversal is read from the source, the counts below are measured by running the function on every case (`stage300_prune_*.jsonl`):

| collection | cases | cases where references were dropped | references kept when it happens | panel size when it happens |
|---|---|---|---|---|
| krepp_ani | 296 | 14 | 3-4 of 100-107 (median 3.738%) | 100-107 |
| krepp_placement | 60 | 0 | - | - |

In those 14 krepp ANI cases every one of the 100-107 references carried a distance (`n_with_distance` equals `n_observed` in all of them) and the tool nevertheless kept 3-4 of them.  Its reported loss for those cases is computed on 3-4 leaves.

Because a comparison on data the tool never saw would not be fair, both sides were recomputed on the panel it did use - Theta_1 restricted to the surviving references with stage 298's exact k = 1 routine, and R_dd on the same leaves:

- full panel, as the practitioner holds it: R_dd >= Theta_1 in **121 of 296** krepp ANI cases, median ratio 0.905;
- restricted to the references the tool used: **107 of 296**, median ratio 0.839, and in the 14 mutilated cases the restricted verdict is 0 of 14 (it fits its 3-4 retained references almost exactly).

The 14-case difference between the two counts is exactly the 14 mutilated cases.  Neither count is cosmetic: the first is what a user gets, the second is what the tool was given.

## 7. Verification

| route | what it checks | result |
|---|---|---|
| A | frozen `min` form vs stage 296's interpolation form, exact rational equality, every evaluated run | 2751/2751 agree |
| B | whole-case recomputation with `scipy.sparse.csgraph.dijkstra` on the edge list with source points inserted as vertices | 16 cases, max abs. difference 8.3e-17 |
| C | prediction rebuilt from DecoDiPhy's own `D`, `C`, `l` matrices on the tree file the bridge writes | 16 cases, max abs. difference 8.3e-17 |
| D | recomputed SSR of the returned configuration vs the loss the tool reports, every run | median relative difference 4.7e-06 (krepp ANI, dd_default); the 496 outliers are section 6.1 |
| E | Lambda_1 vs DecoDiPhy's own k = 1 ladder loss (same objective, same edge set, hill climbing vs exhaustive) | equal to 1e-9 in 2451/2631 krepp ANI and 120/120 placement runs; all 114 exceptions are the 14 cases of section 6.2 |
| F | truth residual recomputed here vs stage 298's imported `residual_max` | True: 296 |
| G | Lambda_1 closed form vs exact re-evaluation of its attaining configuration | max relative difference 2.7e-13 (krepp ANI), 2.7e-13 (placement), 3.3e-12 (travellers) |
| H | the reference panel behind stage 298's Theta_1 against the panel this stage evaluates on, case by case | 385/385 identical, 0 mismatches |

Route B is the independent one asked for: a library shortest-path routine that shares no code with this project, run on a graph in which the source points are explicit vertices, reproduces the leaf equations of whole cases to 1 ulp.

## 8. What this does NOT establish

1. **It does not show that DecoDiPhy is wrong.**  No sample here has a known answer except the 356 semi-synthetic cases, and on those the tool's mixture fits the observations better than the true configuration does.  A worse L-infinity fit than one source is not evidence of a wrong answer; it is evidence that the data does not force a multi-source answer in that norm.
2. **It does not establish a true source count for any stool sample.**  Theta_1 and R_dd are fit quantities; `k` is what the tool reported.  Nothing here licenses "the sample contains one source".
3. **It is not a comparison on DecoDiPhy's home ground.**  Section 4 is the relevant table for least squares, and there the tool wins.  The L-infinity comparison is the one Theta_1 supports, not the one the tool was built for.
4. **It does not measure accuracy of placement or abundance.**  Support recovery and abundance error are stages 289 and 293; this stage measures fit only.
5. **The krepp ANI and placement collections are semi-synthetic**: real estimated distances, but sources planted in a published simulation.  Only the 29 stool samples are real end to end, and they carry no truth.
6. **The two implementation findings in section 6 are properties of DecoDiPhy 1.5.0 as published** (commit and file hashes in the protocol), reproduced from its own source; they are not claims about the method in the paper.
7. **No fit comparison for DecoDiPhy on the stool samples was obtained.**  Section 5 reports a cost, not a result: 6 of 29 samples were attempted under this stage's invocation and none returned a configuration within the registered 2-hour wall, so there is no R_dd, no k and no SSR for any real stool sample here.  Nothing in this document says how the tool would fit them.
8. **No error model is assumed anywhere.**  Theta_1 needs none, which is what makes the stool-sample comparison possible at all; the price is that the comparison is of fit quality, not of correctness.

## 9. Files

- `stage300_protocol.json` - protocol: input, code and package hashes, run plan
- `stage300_harvest.jsonl` - every harvested full-panel DecoDiPhy run with its provenance
- `stage300_runs_placement.jsonl` - new runs, 60 placement cases x 2 arms
- `stage300_runs_travelers.jsonl / _pilot_travelers.jsonl` - stool-sample run logs; both are empty because no run returned before the wall - the outcomes are in the censored record below
- `stage300_travelers_censored.jsonl` - the six stopped stool-sample runs with their measured wall-clock and resident set size
- `stage300_travelers_cost.json` - the registered wall, its justification and the measured cost of the travellers arm
- `stage300_fits_<collection>.jsonl` - per-run exact R_dd, SSR, decisions
- `stage300_baselines_<collection>.jsonl` - per-case Theta_1, truth, Lambda_1
- `stage300_resolved_<collection>.jsonl` - the re-solved fairness arm
- `stage300_prune_<collection>.jsonl` - surviving panel and the fair-panel comparison
- `stage300_lossgap.json / _rows.jsonl` - returned configuration vs reported loss
- `stage300_losscase_*.json` - the worst inconsistency, re-solved at every neighbouring support
- `stage300_route_e.json` - Lambda_1 against the tool's own k = 1 answer
- `stage300_verification_small.json` - routes A, B, C
- `stage300_analysis.json` - all aggregates behind the tables above
- `stage300_tdprobe_*.json` - cost decomposition on one stool sample
- `stage300_three_way.csv` - the per-case three-way table (708 rows: every truth-known case x arm, exact rationals beside their float values)
- `stage300_leafset_consistency.json` - route H: the panels behind Theta_1 and behind R_dd, case by case
- `stage300_protocol_amendment.json` - what was added after registration, with hashes, and the declared deviations

