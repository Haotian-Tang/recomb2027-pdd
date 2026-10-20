# Stage 293 — what missing references cost the abundances and the positions

Stages 288 and 290 measured what reference incompleteness does to the *number* of sources and to the *identity* of the support. This stage measures the two quantities a practitioner actually reports: the weights `p_i` and the interior positions `x_i`, together with the mixture-level error and, where the retained data do not pin the parameters down, the exact size of the set of equally valid answers.

Protocol `results/stage293/stage293_protocol.json`, registered 2026-10-18T06:38:37Z before any run, draw-plan digest `e8bef62351be6c27`.

## 1. Design, and what it inherits

The case panel, the three removal regimes, the five removal levels, the seeds and therefore every individual draw are **imported** from `stage288_incompleteness` (`panel`, `removal_draws`); the registered digest above is over the resulting draw list. Each draw here is the same draw stage 288 measured the count collapse on, so the abundance cost and the count collapse can be read against each other case by case. As in stage 288, reference loss removes *observations*; the reference tree keeps every edge.

This was checked rather than assumed: the 2,525 removal draws measured here are an exact set match, on `(case, level, regime, replicate)`, to the 2,525 draws in stage 288's `B_marks_retained` arm — 0 draws present in one and not the other.

101 cases from the frozen exact families, 2626 draws (including the 101 full-panel baselines), four measured arms:

- **exact** — `solve_partial` on the retained observations: the minimum source count and one minimum-size witness.
- **fibre_true** — the exact fibre over the *true* support given the retained data: is the correct support still feasible, is it still identifiable, and how far apart can two explanations on it be.
- **fibre_found** — the same analysis over the support the solver returned, when that differs from the true support.
- **dd_default / dd_true_k** — DecoDiPhy 1.5.0 (`hill` + `osqp`), documented defaults; `dd_true_k` additionally gets `-k K_true -f 1`.

### The Wasserstein-1 quantity is stage 289's, reused

The mixture-level error is the Wasserstein-1 distance between the recovered and the generating mixture under the tree metric, in branch-length units. It is computed by importing and calling **`stage289_headtohead_analyze.wasserstein`** on a cost matrix built by **`stage289_headtohead_analyze.point_distances`** — the same two functions stage 289 used, with the same input convention (edge lengths as `Fraction`, floated inside). Nothing was reimplemented, so the numbers are directly comparable to stage 289's. Balanced optimal transport compares mixtures with different numbers of atoms natively, which is what makes the count asymmetry tractable.

### The fibre measurement

For a fixed support the model is affine in `(p_i, t_i = p_i x_i, y)`, so the set of explanations is a polytope. `stage293_fibre` reduces the leaf equations to row echelon form in `Fraction`s and runs an exact rational simplex (Bland's rule) to get the supremum and infimum of every coordinate. Reported as the fibre *diameter* is the l-infinity diameter, `max_i (sup p_i - inf p_i)` and likewise for `x`; because the strictly-feasible set is dense in its closure whenever it is non-empty, these suprema are exactly the suprema over valid explanations, approached but not necessarily attained. Position ranges need `x_i = t_i/p_i`, a linear-fractional objective, handled by Charnes-Cooper homogenisation. The fibre is bounded (`sum p_i = 1`, `p_i >= 0`, `0 <= y <= min_r d_r`), so no range is infinite.

## 2. What actually ran, and what did not

| removal level | draws | exact: witness returned | dd_default: completed | dd_true_k: completed |
|---|---|---|---|---|
| full panel | 101 | 85 (84%) | 100 (99%) | 100 (99%) |
| one reference | 505 | 424 (84%) | 498 (99%) | 497 (98%) |
| 10% | 505 | 414 (82%) | 463 (92%) | 458 (91%) |
| 25% | 505 | 404 (80%) | 402 (80%) | 391 (77%) |
| 50% | 505 | 391 (77%) | 311 (62%) | 304 (60%) |
| 75% | 505 | 432 (86%) | 211 (42%) | 197 (39%) |

The exact arm returned a witness on 2150 of 2626 draws; the 476 misses are hard 60 s wall timeouts and are concentrated entirely in the large-K cases (K=7: 156, K=10: 294, K=12: 26), matching the K >= 10 unavailability stage 289 already reported. Availability of the exact arm is therefore flat in the removal level — the misses are a property of the case, not of the loss. DecoDiPhy is the opposite: it answers on 99% of full panels and 42% of 75%-loss panels. That collapse is analysed in section 7.

## 3. Verification

- Every one of the 2150 returned witnesses was re-checked by recomputing the leaf equations on the original tree in exact rational arithmetic (`stage288_common.verify_witness`, independent of any solver flag): 2150/2150 reproduce the retained observations exactly, 0 failures.
- The true support was exactly feasible on all 2626 retained sets, as it must be since the generating mixture explains every subset of the observations. An infeasible verdict would have been a bug in the fibre code, not a finding; none occurred.
- The minimum never exceeded `K_true` on any draw (True), the second consequence of the same fact.
- The generating parameters lay inside the computed fibre bounds on every draw where the bounds were computed, a direct check that the LP brackets the truth.
- The exact simplex was validated against `scipy` HiGHS on 300 random dense LPs (300/300 agreement, 0 mismatches) and returns a p-range of exactly 0 on every identifiable fibre; free-parameter counts agree with the independent `stage289_support_feasibility.feasible_support`. See `results/stage293/stage293_lp_validation.json`.

- **Generalisation check.** Of the 2065 witnesses that explain the retained references exactly, only 771 (37.3%) also reproduce the *removed* references exactly. An exact fit on an incomplete panel is therefore weak evidence about the full panel; among draws that also recovered the correct support the figure is 753/992.

## 4. (a) Abundance and position error when the support is still correct

Restricted to draws where the returned minimum-size support **is** the generating support, with the error measured against the known generating mixture in exact rational arithmetic. `max |p - p_true|` is the l-infinity abundance error, `max |x - x_true|` the l-infinity position error.

| level | regime | support correct / answered | rate | p err med | p err p90 | p err max | x err med | x err p90 | x err max |
|---|---|---|---|---|---|---|---|---|---|
| full panel | - | 83/85 | 98% | 0.0000 | 0.0041 | 0.2390 | 0.0000 | 0.0147 | 0.9400 |
| one reference | random | 165/170 | 97% | 0.0000 | 0.0031 | 0.2390 | 0.0000 | 0.0881 | 0.9400 |
| one reference | whole-clade | 166/170 | 98% | 0.0000 | 0.0082 | 0.2390 | 0.0000 | 0.1775 | 0.9400 |
| one reference | nearest-to-source | 78/84 | 93% | 0.0000 | 0.0000 | 0.1375 | 0.1050 | 0.4481 | 0.5068 |
| 10% | random | 159/167 | 95% | 0.0000 | 0.0054 | 0.2390 | 0.0000 | 0.3064 | 0.9400 |
| 10% | whole-clade | 103/165 | 62% | 0.0000 | 0.0127 | 0.1375 | 0.0000 | 0.0568 | 0.4547 |
| 10% | nearest-to-source | 16/82 | 20% | 0.0000 | 0.0000 | 0.0590 | 0.1753 | 0.4907 | 0.4986 |
| 25% | random | 133/166 | 80% | 0.0000 | 0.0000 | 0.2390 | 0.0000 | 0.4748 | 0.9400 |
| 25% | whole-clade | 53/158 | 34% | 0.0000 | 0.0000 | 0.1375 | 0.0000 | 0.0000 | 0.3240 |
| 25% | nearest-to-source | 7/80 | 9% | 0.0000 | 0.0000 | 0.0000 | 0.3154 | 0.4986 | 0.4986 |
| 50% | random | 74/163 | 45% | 0.0000 | 0.0000 | 0.0140 | 0.1482 | 0.4313 | 0.5719 |
| 50% | whole-clade | 10/152 | 7% | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| 50% | nearest-to-source | 2/76 | 3% | 0.0000 | 0.0000 | 0.0000 | 0.3832 | 0.4375 | 0.4510 |
| 75% | random | 25/172 | 15% | 0.0000 | 0.0000 | 0.0408 | 0.2042 | 0.4986 | 0.6975 |
| 75% | whole-clade | 1/173 | 1% | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| 75% | nearest-to-source | 0/87 | 0% | - | - | - | - | - | - |

The median abundance error is exactly zero in all 15 of 15 populated cells, and so is the median position error in 9 of them. The difference is in *which* cells the position error becomes non-zero: one reference/nearest-to-source 0.105, 10%/nearest-to-source 0.175, 25%/nearest-to-source 0.315, 50%/random 0.148, 50%/nearest-to-source 0.383, 75%/random 0.204. Position error appears under nearest-to-source removal at every level, including the removal of a *single* reference, and under random removal only once half the panel is gone; abundance error never becomes the typical case at all. The tails behave differently again: the largest position error over these cells is 0.940 against 0.239 for abundances.

The read is that when the exact method returns the right support it returns the right weights, and what reference loss takes first is the position — specifically, the offset of a source along its own edge, which is the quantity the nearest leaves resolve. The high-loss cells of this table are also strongly survivorship-selected (2 to 25 surviving draws under the structured regimes), so their small errors are not evidence that heavy loss is benign; section 6 is where the honest uncertainty is.

Sharper: among the 757 support-correct draws whose fibre is 0-dimensional — the data pin the parameters down — the largest abundance error observed anywhere is 0.000e+00. Where the correct support is identifiable at all, the exact method recovers the generating abundances exactly, at every removal level and under every regime. All remaining abundance error lives in the positive-dimensional fibres of section 6.

## 5. (b) Mixture-level error (Wasserstein-1, tree metric)

W1 between the recovered and generating mixture, over **all** answered draws, split by whether the support was recovered. Branch-length units.

| level | regime | n correct | W1 med | W1 p90 | n wrong | W1 med | W1 p90 | W1 max |
|---|---|---|---|---|---|---|---|---|
| full panel | - | 83 | 0.0000 | 0.0002 | 2 | 0.0060 | 0.0096 | 0.0105 |
| one reference | random | 165 | 0.0000 | 0.0005 | 5 | 0.0102 | 0.0105 | 0.0105 |
| one reference | whole-clade | 166 | 0.0000 | 0.0020 | 4 | 0.0060 | 0.0105 | 0.0105 |
| one reference | nearest-to-source | 78 | 0.0002 | 0.0070 | 6 | 0.0083 | 0.0213 | 0.0296 |
| 10% | random | 159 | 0.0000 | 0.0054 | 8 | 0.0106 | 0.0199 | 0.0372 |
| 10% | whole-clade | 103 | 0.0000 | 0.0008 | 62 | 0.0128 | 0.0792 | 0.1487 |
| 10% | nearest-to-source | 16 | 0.0031 | 0.0193 | 66 | 0.0182 | 0.0531 | 0.1510 |
| 25% | random | 133 | 0.0000 | 0.0149 | 33 | 0.0089 | 0.0404 | 0.1105 |
| 25% | whole-clade | 53 | 0.0000 | 0.0000 | 105 | 0.0167 | 0.2503 | 0.4705 |
| 25% | nearest-to-source | 7 | 0.0063 | 0.0346 | 73 | 0.0241 | 0.1045 | 0.3417 |
| 50% | random | 74 | 0.0013 | 0.0187 | 89 | 0.0126 | 0.0502 | 0.0911 |
| 50% | whole-clade | 10 | 0.0000 | 0.0000 | 142 | 0.0339 | 0.2557 | 0.6646 |
| 50% | nearest-to-source | 2 | 0.0219 | 0.0243 | 74 | 0.0504 | 0.1716 | 0.4594 |
| 75% | random | 25 | 0.0007 | 0.0319 | 147 | 0.0256 | 0.0791 | 0.4130 |
| 75% | whole-clade | 1 | 0.0000 | 0.0000 | 172 | 0.0623 | 0.3041 | 0.5982 |
| 75% | nearest-to-source | 0 | - | - | 87 | 0.0572 | 0.2602 | 0.6390 |

## 6. (c) When the fibre is positive-dimensional: how far apart two valid answers can be

Over the true support, 1827 of 2626 draws (69.6%) leave the parameters underdetermined: the retained observations are consistent with a continuum of explanations on the correct support. Median 3 free parameters, up to 21.

| level | positive-dim / evaluated | rate | p diam med | p diam p90 | p diam max | x diam med | x diam p90 | worst p error med | p90 | max |
|---|---|---|---|---|---|---|---|---|---|---|
| full panel | 21/101 | 21% | 0.2271 | 0.4413 | 0.5520 | 0.9755 | 1.0000 | 0.1790 | 0.2948 | 0.5322 |
| one reference | 169/505 | 33% | 0.0451 | 0.3835 | 0.6356 | 1.0000 | 1.0000 | 0.0406 | 0.2948 | 0.5322 |
| 10% | 288/505 | 57% | 0.0000 | 0.3741 | 0.9880 | 1.0000 | 1.0000 | 0.0000 | 0.2800 | 0.9797 |
| 25% | 379/505 | 75% | 0.0000 | 0.5055 | 1.0000 | 1.0000 | 1.0000 | 0.0000 | 0.3861 | 0.9797 |
| 50% | 473/505 | 94% | 0.1293 | 0.7842 | 1.0000 | 1.0000 | 1.0000 | 0.0898 | 0.6234 | 0.9897 |
| 75% | 497/505 | 98% | 0.4207 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 0.3050 | 0.9684 | 0.9984 |

"worst p error" is the largest deviation from the generating abundance that still explains the retained observations exactly — the honest error bar on a reported abundance when the answer is a range.

**This is not an artefact of reference loss.** At the full panel 21 of 101 cases (21%) are already parameter-non-identifiable on the true support, with a median abundance diameter of 0.227 and a maximum of 0.552. Stage 290 reported the same phenomenon on the full bank (70/893 cases); this panel is enriched for large K, which is where non-identifiability concentrates, so the rate is higher. Reference loss widens an uncertainty that is present from the start.

Within the draws that recovered the correct support, 318 of 1075 (29.6%) still have a positive-dimensional fibre, and in 230 of those 318 the position of at least one source is entirely unconstrained (diameter 1.0, the whole edge) while the abundances remain pinned. "Correct support" and "determined parameters" are different properties, and the second fails much more often than the first.

On the 1075 draws where the solver returned a support other than the true one, 990 (92%) of those returned supports are themselves positive-dimensional, with abundance diameter up to 0.951. A collapsed answer is not only the wrong support, it is frequently also an underdetermined one.

## 7. Does the abundance error track the count collapse, or is it separate?

| K_true - minimum | draws | support correct | W1 med | W1 p90 | W1 max | p err med (correct supp.) | x err med (correct supp.) |
|---|---|---|---|---|---|---|---|
| 0 | 1611 | 67% | 0.0003 | 0.0423 | 0.6646 | 0.0000 | 0.0000 |
| 1 | 270 | 0% | 0.0306 | 0.1711 | 0.6057 | - | - |
| 2 | 134 | 0% | 0.0523 | 0.2560 | 0.5269 | - | - |
| 3 | 52 | 0% | 0.0560 | 0.2061 | 0.6390 | - | - |
| 4+ | 83 | 0% | 0.0491 | 0.2689 | 0.5982 | - | - |

**Both, in different coordinates.** At the mixture level the error is almost entirely the count collapse: median W1 is 2.56e-04 when the count survives and 4.34e-02 when it does not, a factor of 170. Losing a source is simply a bigger transport cost than misplacing one. A count deficit also forces the support to be wrong, so the two cannot be separated by conditioning on the support.

But the *abundance* error is not the count collapse in disguise, and it is not a separate effect either in the sense of a smooth degradation with lost references. Conditional on the right support and an identifiable fibre it is exactly zero at every level; conditional on the right support and a positive-dimensional fibre, a reported abundance is not a well-defined quantity at all, and the correct object is the interval of section 6. The thing reference loss degrades smoothly and separately is the **position**: the x error grows with the removed fraction (section 4) and is worst under nearest-to-source removal, which is exactly the regime that strips the leaves whose distances would have resolved a source's offset along its own edge.

## 8. Exact method vs the published heuristic under the same reference loss

Two separate comparisons are needed, because the methods fail in different ways.

### 8.1 Availability

| level | draws | exact: answered | dd_default: answered | dd hard failures | dd writer crashes (recovered) |
|---|---|---|---|---|---|
| full panel | 101 | 84% | 99% | 1 | 0 |
| one reference | 505 | 84% | 99% | 7 | 46 |
| 10% | 505 | 82% | 92% | 42 | 138 |
| 25% | 505 | 80% | 80% | 103 | 198 |
| 50% | 505 | 77% | 62% | 194 | 178 |
| 75% | 505 | 86% | 42% | 294 | 157 |

DecoDiPhy 1.5.0 handles missing references by its own documented route: `preprocess_input` prunes leaves absent from the distance file and contracts their parents, and `post_process_output` maps placements back to the original tree. Two distinct failure modes appear once that path is exercised, neither of them present at the full panel:

1. **Hard failures** — 1320 of 5252 runs (both arms) produce no output at all, raised inside DecoDiPhy on its own defaults: `AttributeError` 878, `cannot be empty` 214, `k too small` 177, `MemoryError` 26, `IndexError` 22, `UnboundLocalError` 3. These are concentrated at high loss, not at large N: 2 failures over the 202 full-panel runs against 602 over the 1010 runs at 75% loss. Only 26 of the failures are `MemoryError` against this stage's 3 GiB address-space cap, so the cap is not the explanation.
2. **Writer crashes** — a further 1400 runs fit successfully and then raise `KeyError` in `save_jplace`, which looks the back-mapped anchors up in the *pruned* tree. Because `all_rounds.json` is written at `cli.py:379` and the writer runs at `cli.py:384`, the placements the tool itself committed are intact; this stage reads them and records `jplace_writer_crashed=true`. Nothing was refitted and no default was changed. A user running the CLI normally would see a traceback and no `.jplace`.

### 8.2 Accuracy, conditional on both methods answering

| level | n both answered | exact W1 med | p90 | dd_default W1 med | p90 | dd_true_k W1 med | p90 | exact strictly better |
|---|---|---|---|---|---|---|---|---|
| full panel | 85 | 0.0000 | 0.0007 | 0.0000 | 0.0079 | 0.0000 | 0.0019 | 83/85 |
| one reference | 422 | 0.0000 | 0.0039 | 0.0000 | 0.0100 | 0.0000 | 0.0056 | 404/422 |
| 10% | 384 | 0.0001 | 0.0280 | 0.0021 | 0.0462 | 0.0005 | 0.0425 | 362/384 |
| 25% | 323 | 0.0050 | 0.0540 | 0.0100 | 0.1091 | 0.0078 | 0.1019 | 296/323 |
| 50% | 239 | 0.0191 | 0.0967 | 0.0278 | 0.1726 | 0.0277 | 0.1788 | 210/239 |
| 75% | 176 | 0.0336 | 0.1727 | 0.0510 | 0.2314 | 0.0612 | 0.2438 | 146/176 |

The exact method is better on the large majority of paired draws at every level, but the *gap closes* as references are removed: at the full panel the p90 ratio dd_default:exact is about 11x, at 75% loss about 1.3x. Read correctly this is not the heuristic degrading gracefully — it is the exact method having further to fall, because it starts at zero error. In absolute terms both degrade, and the heuristic is worse at every level.

The k asymmetry is handled as in stage 289 and is not an incidental detail: the exact method returns the *minimum* count, which the count collapse drives downward, while DecoDiPhy either picks k by its own stopping rules (`dd_default`) or is handed the true k (`dd_true_k`). Giving the heuristic the true k does not rescue it — its median W1 at 75% loss is no better than with its own k, because with the references gone there is no information to place the extra sources with. W1 compares mixtures of different sizes natively, so no arm is penalised for returning the wrong number.

Support recovery, among runs that answered: exact 1075/2150 (50%), dd_default 569/1985 (29%), dd_true_k 772/1947 (40%).

## 9. What this does **not** establish

- Nothing here is about noisy distances. Every case is exactly explainable by construction, and the whole measurement is exact rational feasibility. Real distance estimates carry error (stage 291 measured 1.34e-1 ANI error against a certified-screening tolerance of ~1e-3), and under noise "the true support is still feasible" is not even the right question. The fibre diameters are therefore a lower bound on practical uncertainty, not a model of it.
- The abundance errors in section 4 are conditioned on the exact method returning the true support, which it does on 50% of answered draws overall and on 3% at 75% loss. The high-loss cells of that table rest on few draws (26 support-correct draws at 75%) and are survivorship-selected toward easy cases; they are not evidence that heavy reference loss is harmless.
- The exact arm did not answer on 476 of 2,626 draws (hard 60 s timeouts, all at K >= 7, mostly K >= 10). Every comparison in section 8.2 is conditional on both methods answering, and that conditioning is not random: it removes the large-K cases from both arms. The availability table in 8.1 is the unconditional statement.
- The fibre diameter is an l-infinity diameter over coordinates, computed per coordinate. It is the exact width of the set of valid answers in each parameter separately; it is not a volume, not a posterior, and carries no probability statement about where in the interval the truth lies.
- Positions are reported as x along the generating edge. When the returned support differs from the true support there is no coordinate correspondence at all, so no position error is defined for those draws; they appear only in the W1 and availability tables.
- The DecoDiPhy failure counts are for version 1.5.0 with `hill` + `osqp`, the only operable configuration in this install, on inputs produced by this stage's bridge. They are reproducible failures of the published implementation on incomplete panels, not a claim about the method as published, and not a claim that the underlying optimisation is unsound.
- The panel is 101 of the 893 exactly-explainable bank cases, selected by stage 288's N x K grid, and is enriched for large K relative to the bank. The non-identifiability rates here are therefore higher than the bank-wide rates stage 290 reported, and the two should not be quoted interchangeably.
- Three removal regimes with 2 + 2 + 1 replicates per level do not characterise the worst case. Stage 290 showed the adversarial threshold is the one that matters and that it does not grow with panel size; nothing here contradicts or extends that.
- No claim is made about which support a practitioner would choose when several minimum-size supports exist. The exact solver returns one witness; enumerating all minimum-size supports is the stage 290 bridge computation, which remains limited to K = 2 and 3.

## 10. Files

- `results/stage293/stage293_protocol.json` — registered protocol: panel, draws, limits, code and input hashes
- `results/stage293/stage293_exact_runs.jsonl` — 2,626 exact-solver runs with independent witness verification
- `results/stage293/stage293_fibre_true_runs.jsonl` — 2,626 exact fibre analyses over the true support
- `results/stage293/stage293_fibre_found_runs.jsonl` — 1,075 fibre analyses over the returned support
- `results/stage293/stage293_decodiphy_runs.jsonl` — 5,252 DecoDiPhy runs (both arms), failures retained
- `results/stage293/stage293_records.jsonl` — per-draw joined records including W1
- `results/stage293/stage293_analysis.json` — summary grids
- `results/stage293/stage293_lp_validation.json` — exact-simplex validation against scipy HiGHS
- `results/stage293/stage293_exact_runs.aborted_attempt1.jsonl` — 313 runs from a first launch stopped after 3 min for a too-short job clock; retained, superseded, not used in any number above

Scripts: `stage293_fibre.py` (exact rational fibre LP), `stage293_abundance.py` (protocol + four arms), `stage293_analyze.py` (joining and W1 reuse), `stage293_validate_lp.py`, `stage293_report.py` (this document).
