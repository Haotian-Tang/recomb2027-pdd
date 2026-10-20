# Stage 296 — the ε-tolerant source budget: noise as the controlled variable

*Exactly explainable instances from the stage-287 bank; the published heuristic is not in the
loop anywhere in this stage.*

## 0. What this settles

Stage 294 showed that the exact-model downward count bias does not appear on real ANI distances
and could not separate two explanations: the effect is an artefact of the exact model, or
DecoDiPhy is too noisy an estimator of the source count to see it. Making the tolerance ε the
controlled variable and removing the estimator entirely settles it.

**The bias is a property of the model, not of exactness.** A tolerance alone — no reference
removed — already lowers the model's own minimum source count below the truth, at a median
ε of 0.0147 × the median observed distance
(0.0023 in observation units) over 43 cases. Reference
incompleteness and estimation noise are therefore the same mechanism acting on the same
quantity, and the exact arm of stage 288 is the ε → 0 limit of this stage, reproduced here
exactly (0.202 of 1075 draws fall, 0 rise).

**But the effect is unobservable at current distance-estimation accuracy.** The real ANI error
of the krepp family, measured on the same trees against their exact partners, is
0.15 × the median observed distance at the median leaf and 0.34 × at the
worst leaf — 7.3× to
14.5× above the tolerance at which
the count has already dropped, in 100 %
of cases. At a tolerance wide enough to contain the real error at every leaf, the *model's own*
ε-tolerant minimum is **1 source in 15 of 15 cases**, against a true count of 3 (2 in one case).
There is nothing left for any estimator to recover: the question "how many sources" has no
informative answer at that accuracy, independently of the method used to ask it.

**The honest claim for the paper** is therefore sharper than stage 294's: *in the exact model,
and up to tolerances roughly an order of magnitude below current ANI accuracy, reference
incompleteness removes constituents one-sidedly and which references are missing matters far
more than how many; at the accuracy real ANI estimates actually deliver, the model's own answer
has already collapsed to a single source, so the effect is real but not currently observable.*
Distance estimation would need to be about **7× (median leaf) to 15× (worst leaf) more
accurate** before the effect becomes measurable at all.

A structural consequence worth stating separately: under the exact definition the ε-tolerant
minimum **can never rise** when references are removed, at any ε ≥ 0. Removing an observation
deletes LP rows, so Θ_k can only fall. Observed: 0 rises in 1075 draws at every
tolerance tested. DecoDiPhy's two-sided movement in stage 294 (36.1 % down, 20.9 % up) is
therefore not the model's behaviour at any noise level; it is estimator variance on top of a
quantity that has already collapsed.

![Left: where the first-drop tolerance sits against real and simulated distance error. Right: within-ε fall rate by removal regime, with the share of cases whose baseline has already collapsed.]({{artifact:3afb0167-3b5d-41e2-bd3b-2b86e32e5c10}})

## 1. Definition and method

**Problem.** Given observations d_r on a known reference tree T and a tolerance ε ≥ 0, the
ε-tolerant source budget is the minimum number k of strict interior sources on k distinct
original edges, with weights p_i > 0 summing to one and a common mean pendant y > 0, such that
|y + Σ_i p_i·dist_T(r, q_i) − d_r| ≤ ε for every observed r.

**Affine form.** For a point q at fraction x of edge (u, v) and any vertex r,
dist(r, q) = (1−x)·d(r, u) + x·d(r, v) exactly. With endpoint masses α_i = p_i(1−x_i) and
β_i = p_i·x_i the model is affine in (α, β, y) ≥ 0 under Σ_i(α_i + β_i) = 1, so for a fixed
support the tolerance gives 2|R| inequality rows and feasibility is exact rational linear
feasibility.

**No ε ladder was run.** For a fixed support S let θ(S) be the minimum over the *closed* set
{α, β, y ≥ 0, Σ = 1} of max_r |pred_r − d_r| — one exact LP — and Θ_k = min over all C(|E|, k)
supports of size k. Two facts make the whole ε axis a function of the exact rational vector
(Θ_1, Θ_2, …):

1. *For every ε > θ(S) the strict system on S is feasible.* If z_c attains θ(S) and z_s is any
   strictly feasible box point with finite deviation M, then (1−λ)z_c + λz_s is strictly
   feasible for every λ ∈ (0, 1] and, the prediction being affine, has deviation at most
   (1−λ)θ(S) + λM → θ(S).
2. *Θ_k is non-increasing in k*, because the closed system lets a size-k support reproduce any
   size-(k−1) solution with the extra endpoint masses set to zero.

Hence m(ε) = min{k : Θ_k < ε} for ε > 0. m(0) is decided separately, because the closed LP
admits degenerate configurations: it is the smallest k for which some support with θ(S) = 0 is
*strictly* feasible, decided by a max-common-slack LP whose optimum must be > 0. Θ_K = 0 holds
for every leaf set and is certified per job by verifying the generating witness on the retained
observations, so sweeps stop at k = K − 1. The ladder in the tables below is a reporting choice
applied to the exact thresholds afterwards.

**The LP.** A dense two-phase tableau simplex over `fractions.Fraction`, written for this stage
(`stage296_common.lp_minimize`): Dantzig pricing with a switch to Bland's rule after 60
iterations so termination is guaranteed under degeneracy. Row generation over the 2|R| tolerance
rows — the active-subset LP is a relaxation, so its optimum is an exact lower bound, which both
certifies optimality when no leaf is violated and prunes exactly against the incumbent. Supports
are enumerated exhaustively; no structural shortcut is used.

**The design constraint was respected.** None of the frozen structural localisation (leaf
envelope h, weighted deficit bound, positive-anchor filters, and therefore `solve_partial`,
`stage292_all_minimisers`, `stage295_reduced_minimum`) is used for ε > 0 — all of it is derived
from exact equality and is unsound under a tolerance. The frozen machinery appears only in the
ε = 0 validation, where it is sound.

**Independent verification.** Every reported witness is re-evaluated on the original tree with
the frozen `min(d(r,u) + xL, d(r,v) + (1−x)L)` form — an implementation independent of the
interpolation the LP uses — and the recomputed max |pred − d| must equal θ exactly.
Mismatches: **0** over all 6,611,481 enumerated supports and
1118 + 18 jobs.

**Design.** Removal regimes, levels, seeds and the draw plan are imported from
`stage288_incompleteness.removal_draws`, not re-derived, so the cells are comparable with stage
288: levels {1 reference, 10 %, 25 %, 50 %, 75 %} × regimes {random ×2, whole clade ×2,
nearest-to-sources ×1}. Cross arm: 43 exactly explainable cases with K ∈ {2, 3},
|R| from 27 to 361, each at the full panel and at
all 25 registered draws — 1118 leaf sets. Probe arm: 18 cases with
K ∈ {5, 7, 10} at the full panel only.

## 2. ε = 0 validation against the frozen machinery

| check | result |
|---|---|
| leaf sets compared (full panel + every draw) | 1118 |
| m(0) from the new LP equals `solve_partial`'s minimum | **1118 / 1118** |
| disagreements | **0** |
| frozen runs censored or missing | 0 |
| witness recomputation mismatches | 0 |

Every one of the 1118 cross-arm leaf sets was also solved by the frozen
`solve_partial` (all 1118 completed; total 27 s). The new ε-tolerant LP
agrees with the frozen exact minimum on every leaf set, including the 1
case whose full-panel minimum is below its generating K. No disagreement arose, so no
stop-and-report condition was triggered.

## 3. Noise alone: the ε at which the count first drops

Threshold = Θ_{m(0)−1}, the exact smallest tolerance at which the ε-tolerant minimum falls
below the exact minimum, for the full panel (no references removed).

| group | cases | min | Q1 | median | Q3 | max |
|---|---|---|---|---|---|---|
| all cross-arm cases, relative to median observed distance | 43 | 1.6e-05 | 0.0066 | **0.0147** | 0.0273 | 0.246 |
| all cross-arm cases, absolute (observation units) | 43 | 3.9e-06 | 0.001 | **0.0023** | 0.00592 | 0.104 |
| K = 2 (relative) | 16 | 0.000853 | 0.00839 | **0.0208** | 0.0295 | 0.246 |
| K = 3 (relative) | 27 | 1.6e-05 | 0.00594 | **0.0135** | 0.0267 | 0.1 |
| exact_10k (relative) | 15 | 1.6e-05 | 0.0171 | **0.0227** | 0.04 | 0.1 |
| exact_biotrees (relative) | 28 | 0.000363 | 0.00594 | **0.0119** | 0.0234 | 0.246 |
| &#124;R&#124; 20-49 (relative) | 30 | 1.6e-05 | 0.00594 | **0.0135** | 0.0208 | 0.143 |
| &#124;R&#124; 50-129 (relative) | 10 | 0.00448 | 0.0126 | **0.024** | 0.04 | 0.1 |
| &#124;R&#124; 130+ (relative) | 3 | 0.0066 | 0.0295 | **0.0295** | 0.246 | 0.246 |

**Scaling.**

| against | Spearman ρ | p | n |
|---|---|---|---|
| minimum abundance in the generating mixture | **+0.665** | 1.1e-06 | 43 |
| minimum abundance, within K = 3 only | +0.643 | 0.0003 | 27 |
| reference-set size &#124;R&#124; | +0.187 | 0.23 | 43 |
| reference-set size, within K = 3 only | +0.180 | 0.37 | 27 |
| source count K (cross arm, K ∈ {2, 3}) | −0.174 | 0.26 | 43 |
| source count K (adding the probe cases where m(0) is determined) | **−0.365** | 0.012 | 47 |

The threshold is governed by the **smallest constituent**, not by the size of the reference
panel: ρ = +0.665 against the minimum generating abundance
(p = 1e-06), against ρ = +0.187 (n.s.) for |R|. A
mixture with a 1 % constituent loses it at a tolerance an order of magnitude smaller than one
whose smallest constituent carries 36 %. Adding more references does not protect the count.

Against the source count, the cross arm alone (K ∈ {2, 3}) is not conclusive, but the probe
arm resolves it: for the four K = 5 cases where the exhaustive k ≤ 4 sweep establishes
m(0) = 5, the first drop occurs at a *median* of 1.4 × 10⁻³ relative — an order of magnitude
below the K = 2 (2.1 × 10⁻²) and K = 3 (1.3 × 10⁻²) medians, ρ = −0.365 (p = 0.012, n = 47).
Losing *one* constituent gets easier as K grows, while collapsing all the way to two sources
gets harder (median Θ₂ rises from 1.2 × 10⁻² at K = 3 to 4.3 × 10⁻² at K = 10).

## 4. Noise against incompleteness, compared within ε

Every draw is compared against **its own full-panel answer at the same ε**. Rates are over all
1075 draws; "fall rate where baseline > 1" restricts to draws whose full-panel
answer has not yet collapsed to a single source, which is the only subset where incompleteness
can still be detected.

| ε / median distance | draws | mean m(ε) full panel | mean m(ε) draw | fall rate | rise rate | baseline already 1 | fall rate where baseline > 1 | random | clade | nearest | clade/random | nearest/random |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 (exact) | 1075 | 2.605 | 2.355 | 0.202 | 0.000 | 0.000 | 0.202 | 0.081 | 0.279 | 0.288 | 3.43 | 3.54 |
| 1e-4 | 1075 | 2.581 | 2.333 | 0.205 | 0.000 | 0.000 | 0.205 | 0.088 | 0.279 | 0.288 | 3.16 | 3.26 |
| 3.16e-4 | 1075 | 2.581 | 2.322 | 0.216 | 0.000 | 0.000 | 0.216 | 0.105 | 0.284 | 0.302 | 2.71 | 2.89 |
| 1e-3 | 1075 | 2.535 | 2.282 | 0.211 | 0.000 | 0.023 | 0.216 | 0.102 | 0.274 | 0.302 | 2.68 | 2.96 |
| 3.16e-3 | 1075 | 2.512 | 2.221 | 0.245 | 0.000 | 0.023 | 0.251 | 0.133 | 0.295 | 0.367 | 2.23 | 2.77 |
| 1e-2 | 1075 | 2.233 | 1.945 | 0.252 | 0.000 | 0.116 | 0.285 | 0.144 | 0.288 | 0.395 | 2.00 | 2.74 |
| 3.16e-2 | 1075 | 1.721 | 1.532 | 0.175 | 0.000 | 0.372 | 0.279 | 0.091 | 0.200 | 0.293 | 2.21 | 3.23 |
| 1e-1 | 1075 | 1.349 | 1.230 | 0.119 | 0.000 | 0.651 | 0.341 | 0.056 | 0.130 | 0.223 | 2.33 | 4.00 |
| 3.16e-1 | 1075 | 1.000 | 1.000 | 0.000 | 0.000 | 1.000 | 0.000 | 0.000 | 0.000 | 0.000 | — | — |
| 1 | 1075 | 1.000 | 1.000 | 0.000 | 0.000 | 1.000 | 0.000 | 0.000 | 0.000 | 0.000 | — | — |
| real ANI, median leaf* | 375 | 1.200 | 1.128 | 0.072 | 0.000 | 0.800 | 0.360 | 0.027 | 0.087 | 0.133 | 3.25 | 4.99 |
| real ANI, worst leaf* | 375 | 1.000 | 1.000 | 0.000 | 0.000 | 1.000 | 0.000 | 0.000 | 0.000 | 0.000 | — | — |

\* The two reality rows use, per case, that case's own measured ANI residual as ε (15 of the 43
cases have a krepp partner, giving 375 draws).

| removal level | draws | ε=0 | ε=1e-3 | ε=1e-2 | ε=3.16e-2 | ε=1e-1 | ε=3.16e-1 | falls at some ε |
|---|---|---|---|---|---|---|---|---|
| 1 reference | 215 | 0.000 | 0.005 | 0.019 | 0.019 | 0.009 | 0.000 | 0.149 |
| 10 % | 215 | 0.042 | 0.056 | 0.074 | 0.051 | 0.042 | 0.000 | 0.340 |
| 25 % | 215 | 0.112 | 0.116 | 0.153 | 0.112 | 0.093 | 0.000 | 0.563 |
| 50 % | 215 | 0.298 | 0.307 | 0.372 | 0.256 | 0.181 | 0.000 | 0.856 |
| 75 % | 215 | 0.558 | 0.572 | 0.642 | 0.437 | 0.270 | 0.000 | 0.986 |

Reading:

* **The regime ordering survives all the way up.** Structured loss stays 2.0–3.4× (whole clade)
  and 2.7–4.0× (nearest to sources) worse than random removal at every tolerance where the
  comparison is non-vacuous — including at the real-ANI median-leaf tolerance (3.2× and 5.0×).
  At ε = 0 the ratios are 3.43 and
  3.54, reproducing stage 288's 3.15 / 2.88 on an
  independently drawn panel.
* **The incompleteness signal does not blur into the noise-induced collapse; it is extinguished
  by it.** Conditional on the baseline still being informative, the fall rate *rises* with ε
  (0.202 at ε = 0 → 0.341 at
  ε = 10⁻¹). What kills the measurement is vacuity of the baseline: the share of draws whose
  full-panel answer is already 1 goes 0.116 → 
  0.372 → 0.651 → 1.000
  across ε = 10⁻², 10^−1.5, 10⁻¹, 10^−0.5.
* **The signal dies between ε = 10⁻¹ and ε = 10^−0.5** relative to the median observed distance
  — i.e. at a tolerance of 10–30 % of a typical distance, every case's full-panel answer is a
  single source and no incompleteness effect of any kind can be observed.
* **Nothing ever rises.** rise rate is exactly 0 at every ε, as the monotonicity of Θ_k requires.
* **A tolerance makes one new thing visible**: single-reference removal, which never changes the
  exact answer (fall rate 0.000 at ε = 0), does change it for 1.9 % of draws at ε = 10⁻² and for
  0.149 of draws at some tolerance.
* Across all 1075 draws, 57.9 % have a non-empty
  fall region (ε-interval union where the draw's answer is below its own baseline);
  20.2 % fall already at ε = 0. By regime:
  random 50.0 %,
  clade 51.2 %,
  nearest 87.0 %.

## 5. Where real estimation error sits on that axis

For each case the real error is measured against its own exact partner on the shared leaves:
ε_real = max_r |d_estimated(r) − d_exact(r)| (worst leaf) is exactly the smallest tolerance at
which the true configuration still explains the estimated data; the median-leaf residual is
reported alongside because it is the quantity stage 294 quoted (0.134 for the krepp family).

| source | median-leaf residual ÷ median observed distance | worst-leaf residual ÷ median observed distance |
|---|---|---|
| real ANI (krepp), 15 panel cases | 0.152 (IQR 0.101–0.265) | 0.338 (IQR 0.251–0.566) |
| **first drop below the true count** | **0.0147** (IQR 0.0066–0.0273) | — |

| noise source | cases | median-leaf residual ÷ first-drop ε | worst-leaf residual ÷ first-drop ε | cases whose median residual exceeds the first drop |
|---|---|---|---|---|
| real ANI (krepp), measured against the exact partner | 15 | 7.34 (4.29–12.3) | 14.5 (10.9–31.6) | 100 % |
| published noise-free partner (round-trip control) | 28 | 7.8e-15 (3.07e-15–1.52e-14) | 3.31e-14 (1.77e-14–6.84e-14) | 0 % |
| simulated noise 1, 100 bp | 28 | 6.34 (3.14–17.3) | 23.5 (14.8–32.3) | 93 % |
| simulated noise 1, 10 kb | 28 | 6.22 (3.13–17.4) | 23.2 (14.9–31.9) | 93 % |
| simulated noise 2, 100 bp | 28 | 0.249 (0.0713–0.478) | 1.25 (0.4–1.66) | 18 % |
| simulated noise 2, 10 kb | 28 | 0.00249 (0.000713–0.00478) | 0.0125 (0.004–0.0166) | 0 % |

**The model's own answer at real accuracy** (full panel, no references removed, 15 exact_10k
cases with a krepp partner):

| tolerance | m(ε) = 1 | m(ε) = 2 | m(ε) = 3 |
|---|---|---|---|
| ε = 0 (exact) | 0 | 1 | 14 |
| ε = that case's median-leaf ANI residual | 12 | 3 | 0 |
| ε = that case's worst-leaf ANI residual | **15** | 0 | 0 |

Real ANI error sits **far above** the ε at which the signal dies, as stage 291's screening
collapse at η ≈ 10⁻³ suggested. The simulated ladder brackets the threshold and shows the
regime is reachable in principle: simulated noise 2 at 10 kb sits a factor ~400 *below* the
first drop (0 % of cases above it) and leaves the count intact; simulated noise 2 at 100 bp
straddles it (18 % above); simulated noise 1 and real ANI are far above it. The published
noise-free partner reproduces the exact observations to 10⁻¹⁴ of the threshold, which validates
the pairing used for this calibration.

**How much accuracy would be needed.** To bring the worst-leaf error below the first-drop
tolerance, ANI-derived distances would have to be
14.5× more accurate at the median
case (10.9–31.6×
over the interquartile range); using the median leaf instead,
7.3×
(4.3–12.3×).
That is roughly one order of magnitude.

## 6. The scale ceiling that was actually reached

Exhaustive enumeration costs C(|E|, k) exact LPs; the frozen localisation that would avoid it is
unsound under a tolerance, so there is no shortcut available here.

| quantity | value |
|---|---|
| supports enumerated, each decided by an exact rational LP | **6,611,481** |
| by size | k=1: 137,578 · k=2: 3,957,806 · k=3: 581,616 · k=4: 1,934,481 |
| largest single exhaustive sweep | 341,055 supports |
| largest k reached exhaustively | **4** |
| largest reference set | 361 leaves (719 edges) |
| total child CPU | 14.4 h over 32 workers |
| measured rate | 5.4–8.9 ms per support (&#124;R&#124; = 25 → 129) |

The practical frontier at ≈1 h of single-core time per leaf set:
k = 4 at |R| ≈ 25–29 (|E| ≤ 55), k = 3 at |R| ≈ 46 (|E| ≤ 89), k = 2 at |R| ≈ 110 (|E| ≈ 213),
k = 1 at |R| = 361 (|E| = 719). The cross arm therefore used K ∈ {2, 3}, where k ≤ K − 1 ≤ 2
makes the full ε axis computable on every leaf set; for K ∈ {5, 7, 10} only the full panel was
swept, and m(0) is determined only where the sweep reached K − 1 (4 of 18 probe cases).

| K | &#124;R&#124; | edges | k swept exhaustively | m(0) | Θ₁ | Θ₂ | Θ₃ | Θ₄ |
|---|---|---|---|---|---|---|---|---|
| 5 | 25 | 47 | ≤4 | 5 | 0.0424 | 0.0337 | 0.0102 | 9.96e-05 |
| 5 | 25 | 47 | ≤4 | 5 | 0.123 | 0.0176 | 0.0163 | 0.00178 |
| 5 | 25 | 47 | ≤4 | 5 | 0.193 | 0.0118 | 0.00505 | 0.0015 |
| 5 | 29 | 55 | ≤4 | 5 | 0.0501 | 0.0178 | 0.0105 | 0.00133 |
| 5 | 31 | 59 | ≤3 | [4, 5] | 0.0711 | 0.0216 | 0.0193 | — |
| 5 | 32 | 61 | ≤3 | [4, 5] | 0.0974 | 0.0382 | 0.00965 | — |
| 5 | 48 | 93 | ≤2 | [3, 5] | 0.149 | 0.0204 | — | — |
| 5 | 51 | 99 | ≤2 | [3, 5] | 0.052 | 0.025 | — | — |
| 7 | 23 | 43 | ≤4 | [5, 7] | 0.0996 | 0.0523 | 0.0147 | 0.0137 |
| 7 | 29 | 55 | ≤4 | [5, 7] | 0.119 | 0.0578 | 0.0215 | 0.00967 |
| 7 | 41 | 79 | ≤3 | [4, 7] | 0.0209 | 0.00724 | 0.0055 | — |
| 7 | 41 | 79 | ≤3 | [4, 7] | 0.158 | 0.016 | 0.0104 | — |
| 10 | 20 | 37 | ≤4 | [5, 10] | 0.0361 | 0.0302 | 0.0236 | 0.00865 |
| 10 | 20 | 37 | ≤4 | [5, 10] | 0.089 | 0.0731 | 0.0156 | 0.00978 |
| 10 | 26 | 49 | ≤4 | [5, 10] | 0.0948 | 0.0441 | 0.0267 | 0.0139 |
| 10 | 27 | 51 | ≤4 | [5, 10] | 0.0652 | 0.0301 | 0.0254 | 0.0106 |
| 10 | 40 | 77 | ≤3 | [4, 10] | 0.0969 | 0.0566 | 0.0198 | — |
| 10 | 46 | 89 | ≤3 | [4, 10] | 0.0553 | 0.041 | 0.0294 | — |

## 7. What this does **not** establish

1. **Nothing about DecoDiPhy or any other estimator.** This stage measures the model's own
   ε-tolerant minimum. That the heuristic's inferred count moves two-sidedly (stage 294) is not
   explained here beyond showing that two-sided movement is impossible for the exact quantity,
   so it must come from the estimator.
2. **No K > 3 cross grid.** Reference incompleteness crossed with ε was measured only for
   K ∈ {2, 3}. The probe arm gives Θ₁…Θ₄ at the full panel for K ∈ {5, 7, 10} but no draws,
   so the regime ordering is untested above K = 3.
3. **m(0) is bracketed, not determined, for 14 of the 18 probe cases**, because the sweep
   stopped below K − 1. Their Θ₃ and Θ₄ are exact, but their first-drop threshold is not
   observed — only bounded above by the largest swept Θ.
4. **ε is an L∞ tolerance.** A median residual is not the ε a dataset needs; both are reported
   and they differ by a factor of ~2 on real ANI data.
5. **Exactly explainable instances only.** A tolerance threshold is defined relative to an exact
   reference answer, so the noisy bank families cannot carry this measurement; they enter only
   as the calibration partners of their exact twins.
6. **The calibration is a residual, not an error model.** ε_real is the deviation between the
   estimated distances and the exact partner's model-implied distances on the same geometry. It
   does not decompose into estimator bias, tree error, or the bank's equal-weight assumption,
   and for the biotrees families stage 291's caveat on unpublished source weights still stands
   (here the noise-free round-trip is clean, 10⁻¹⁴, but that only validates the data pairing).
7. **No claim about which support is recovered.** Only the minimum count is measured; position
   recovery under a tolerance is untouched (stage 293 did that at ε = 0).
8. **The convention at ε exactly equal to a threshold.** m(ε) = min{k : Θ_k < ε} uses strict
   inequality, the conservative (larger-m) choice; at the measure-zero coincidence ε = Θ_k the
   strict system may or may not be feasible. No tabulated cell depends on it.
9. **No sampling claim.** The panel is a registered deterministic selection from the stage-287
   bank (evenly spaced within family × K × |R| cells), not a random sample of trees, so the
   quartiles describe this panel.

## 8. Files

Under `results/stage296/` on the project host:

| file | contents |
|---|---|
| `stage296_protocol.json` | registered before any timed run: panel, grid, definitions, budgets, input and code hashes |
| `stage296_cross_runs.jsonl` | 1118 cross-arm records (exact Θ vectors, witnesses, timings), every run retained |
| `stage296_probe_runs.jsonl` | 18 ceiling-probe records |
| `stage296_validate_runs.jsonl` | 1118 frozen `solve_partial` runs for the ε = 0 validation |
| `stage296_calibration.json` | residuals of every estimated partner against its exact twin |
| `stage296_analysis.json` | all tables above, per-case and per-draw |
| `scripts/stage296_common.py` | exact rational LP, row generation, exhaustive sweeps, independent verifier |
| `scripts/stage296_tolerance.py` | protocol, cross, probe, validate, calibrate phases |
| `scripts/stage296_analyze.py` | analysis |
| `scripts/stage296_pilot.py` | timing pilot used to size the panel (measurement only) |

## Appendix — per-case thresholds (cross arm, full panel)

All values relative to that case's median observed distance. "first drop" = Θ_{m(0)−1}.

| case | family | K | &#124;R&#124; | min p | m(0) | Θ₁ / median d | Θ₂ / median d | first drop | real ANI median-leaf | real ANI worst-leaf |
|---|---|---|---|---|---|---|---|---|---|---|
| k10k_tree_16_k3_s1_exp | 10k | 3 | 27 | 0.157 | 3 | 0.0637 | 0.0464 | **0.0464** | 0.0995 | 0.229 |
| k10k_tree_16_k3_s2_exp | 10k | 3 | 27 | 0.019 | 3 | 0.181 | 0.00663 | **0.00663** | 0.102 | 0.279 |
| k10k_tree_16_k3_s3_exp | 10k | 3 | 27 | 0.145 | 3 | 0.0565 | 0.0185 | **0.0185** | 0.0515 | 0.226 |
| k10k_tree_16_k3_s4_uniform | 10k | 3 | 27 | 0.220 | 3 | 0.0291 | 0.0195 | **0.0195** | 0.143 | 0.202 |
| k10k_tree_16_k3_s5_uniform | 10k | 3 | 27 | 0.159 | 3 | 0.0394 | 0.00861 | **0.00861** | 0.106 | 0.272 |
| k10k_tree_17_k3_s1_exp | 10k | 3 | 31 | 0.157 | 2 | 0.014 | 0 | **0.014** | 0.161 | 0.218 |
| k10k_tree_17_k3_s4_uniform | 10k | 3 | 31 | 0.220 | 3 | 0.0373 | 0.0273 | **0.0273** | 0.152 | 0.285 |
| k10k_tree_14_k3_s3_exp | 10k | 3 | 34 | 0.145 | 3 | 0.0984 | 0.0267 | **0.0267** | 0.258 | 0.338 |
| k10k_tree_9_k3_s2_exp | 10k | 3 | 47 | 0.019 | 3 | 0.123 | 1.6e-05 | **1.6e-05** | 0.272 | 0.756 |
| k10k_tree_9_k3_s5_uniform | 10k | 3 | 47 | 0.159 | 3 | 0.125 | 0.0171 | **0.0171** | 0.3 | 0.551 |
| k10k_tree_8_k3_s1_exp | 10k | 3 | 50 | 0.157 | 3 | 0.181 | 0.1 | **0.1** | 0.429 | 0.64 |
| k10k_tree_10_k3_s1_exp | 10k | 3 | 58 | 0.157 | 3 | 0.104 | 0.04 | **0.04** | 0.476 | 0.58 |
| k10k_tree_12_k3_s5_uniform | 10k | 3 | 58 | 0.159 | 3 | 0.0943 | 0.0409 | **0.0409** | 0.226 | 0.444 |
| k10k_tree_15_k3_s1_exp | 10k | 3 | 108 | 0.157 | 3 | 0.0567 | 0.0227 | **0.0227** | 0.0963 | 0.417 |
| k10k_tree_11_k3_s5_uniform | 10k | 3 | 129 | 0.159 | 3 | 0.0835 | 0.0261 | **0.0261** | 0.0588 | 0.591 |
| bio_bees-bosseret__bees-bosseret.tre_k2_s10_no | biotrees | 2 | 28 | 0.026 | 2 | 0.0208 | — | **0.0208** | — | — |
| bio_bees-bosseret__bees-bosseret.tre_k2_s2_noi | biotrees | 2 | 28 | 0.056 | 2 | 0.0421 | — | **0.0421** | — | — |
| bio_bees-bosseret__bees-bosseret.tre_k2_s4_noi | biotrees | 2 | 28 | 0.361 | 2 | 0.0321 | — | **0.0321** | — | — |
| bio_bees-bosseret__bees-bosseret.tre_k2_s7_noi | biotrees | 2 | 28 | 0.089 | 2 | 0.00659 | — | **0.00659** | — | — |
| bio_bees-bosseret__bees-bosseret.tre_k2_s9_noi | biotrees | 2 | 28 | 0.020 | 2 | 0.00604 | — | **0.00604** | — | — |
| bio_mammals-song__mammals-song-castlespro.tre_ | biotrees | 2 | 34 | 0.026 | 2 | 0.018 | — | **0.018** | — | — |
| bio_mammals-song__mammals-song-castlespro.tre_ | biotrees | 2 | 34 | 0.203 | 2 | 0.143 | — | **0.143** | — | — |
| bio_birds-jarvis__birds-jarvis.tre.rerooted_k2 | biotrees | 2 | 46 | 0.026 | 2 | 0.00839 | — | **0.00839** | — | — |
| bio_birds-jarvis__birds-jarvis.tre.rerooted_k2 | biotrees | 2 | 46 | 0.361 | 2 | 0.0147 | — | **0.0147** | — | — |
| bio_birds-jarvis__birds-jarvis.tre.rerooted_k2 | biotrees | 2 | 46 | 0.020 | 2 | 0.000853 | — | **0.000853** | — | — |
| bio_tilapia-Ciezarek__tilapia-Ciezarek.tre.rer | biotrees | 2 | 89 | 0.026 | 2 | 0.024 | — | **0.024** | — | — |
| bio_pancrustacean-Bernot__pancrustacean-Bernot | biotrees | 2 | 103 | 0.026 | 2 | 0.0234 | — | **0.0234** | — | — |
| bio_beetles-Johnson__beetles-Johnson.tre_k2_s9 | biotrees | 2 | 125 | 0.020 | 2 | 0.0126 | — | **0.0126** | — | — |
| bio_fish-Troyer__fish-Troyer.tre.rerooted_k2_s | biotrees | 2 | 150 | 0.026 | 2 | 0.0295 | — | **0.0295** | — | — |
| bio_mammals-foley__mammals-foley-concat-heutra | biotrees | 2 | 239 | 0.361 | 2 | 0.246 | — | **0.246** | — | — |
| bio_birds-stiller__birds-stiller-castlepro.tre | biotrees | 2 | 361 | 0.020 | 2 | 0.0066 | — | **0.0066** | — | — |
| bio_bees-bosseret__bees-bosseret.tre_k3_s10_no | biotrees | 3 | 27 | 0.015 | 3 | 0.0482 | 0.0051 | **0.0051** | — | — |
| bio_bees-bosseret__bees-bosseret.tre_k3_s2_noi | biotrees | 3 | 27 | 0.026 | 3 | 0.201 | 0.0119 | **0.0119** | — | — |
| bio_bees-bosseret__bees-bosseret.tre_k3_s4_noi | biotrees | 3 | 27 | 0.220 | 3 | 0.0782 | 0.0135 | **0.0135** | — | — |
| bio_bees-bosseret__bees-bosseret.tre_k3_s7_noi | biotrees | 3 | 27 | 0.059 | 3 | 0.259 | 0.00373 | **0.00373** | — | — |
| bio_bees-bosseret__bees-bosseret.tre_k3_s9_noi | biotrees | 3 | 27 | 0.010 | 3 | 0.226 | 0.00258 | **0.00258** | — | — |
| bio_mammals-song__mammals-song-castlespro.tre_ | biotrees | 3 | 33 | 0.015 | 3 | 0.112 | 0.00594 | **0.00594** | — | — |
| bio_mammals-song__mammals-song-castlespro.tre_ | biotrees | 3 | 33 | 0.159 | 3 | 0.128 | 0.0288 | **0.0288** | — | — |
| bio_birds-jarvis__birds-jarvis.tre.rerooted_k3 | biotrees | 3 | 45 | 0.015 | 3 | 0.0148 | 0.00428 | **0.00428** | — | — |
| bio_birds-jarvis__birds-jarvis.tre.rerooted_k3 | biotrees | 3 | 45 | 0.220 | 3 | 0.162 | 0.00785 | **0.00785** | — | — |
| bio_birds-jarvis__birds-jarvis.tre.rerooted_k3 | biotrees | 3 | 45 | 0.010 | 3 | 0.027 | 0.000363 | **0.000363** | — | — |
| bio_tilapia-Ciezarek__tilapia-Ciezarek.tre.rer | biotrees | 3 | 88 | 0.015 | 3 | 0.134 | 0.0109 | **0.0109** | — | — |
| bio_beetles-Johnson__beetles-Johnson.tre_k3_s9 | biotrees | 3 | 124 | 0.010 | 3 | 0.218 | 0.00448 | **0.00448** | — | — |