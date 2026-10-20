# Stage 292 — Deciding uniqueness of the minimum-source explanation, and the erasure threshold at K ≥ 5 and at panel scale

**Question this stage was opened to settle.** Stage 290 could measure the unknown-support
erasure threshold only at K = 2 (|R| ≤ 130) and K = 3 (|R| ≤ 45), because it decided
uniqueness by enumerating every edge subset — C(|E|, K) exact feasibility tests. The
manuscript's most quotable sentence, *the adversarial (nearest-to-source) threshold does not
improve with panel size*, was therefore established only on 130 of 893 bank cases, at two
values of K, on trees of at most 130 references.

**Outcome.** A new exact all-minimisers routine replaces the brute force, agrees with it on
every one of 5,839 head-to-head comparisons, and extends the measurement to 676 bank cases at
K ∈ {2, 3, 5, 7} and to planted mixtures on reference panels from 64 to 47,151 leaves. The
adversarial threshold has a pooled Spearman correlation of **0.096** with reference-set size on
the bank and **0.045** with panel size on the ladder, against **0.83** and **0.98** for uniform
random loss, whose log–log slope in panel size is **1.0**. The flatness claim is now
established, not provisional, over 2.9 orders of magnitude in panel size and for K up to 7.

![Panel (a): threshold versus panel size on planted mixtures, log–log, nearest-to-source versus uniform random removal, K = 3 and K = 5, 70 jobs over sizes 64–47,151. Panel (b): median and interquartile range of the threshold on the 663 frozen-bank cases whose true support is the unique minimiser at the full panel, by K and removal regime. Higher is more robust. The dotted guides in (a) are slope 1 and slope 0, not fits.]({{artifact:bc4195ce-758e-47f2-bc20-8f2d7bd48dbe}})

---

## 1. The all-minimisers routine

`scripts/stage292_all_minimisers.py`, entry point
`all_minimum_supports(edges, observed, cap)`. Exact rational arithmetic throughout; every
frozen module is imported, never modified or re-implemented.

### 1.1 The reduction that makes the problem combinatorial

For a source at fraction *x* of edge (u, v) and **any** vertex r,
dist(r, q) = (1 − x)·d(r, u) + x·d(r, v). A configuration is therefore equivalent to a vertex
mass vector m ≥ 0 with Σ m = 1 and d_r = y + Σ_v m_v d(r, v). A set F of k edges is a support
of that configuration exactly when m is strictly positive precisely on the endpoint set
U = V(F); conversely, any m strictly positive precisely on V(F) splits over F with every x_i
strictly interior. So

> F is a support ⟺ V(F) is feasible,

and feasibility of an endpoint set is one exact rational linear-feasibility problem — the
frozen Stage 25 oracle `fixed_support_polytope`, which demands all endpoint masses **and** the
pendant strictly positive. The entire combinatorial content is *which endpoint sets to test*.
Brute force tests all C(|E|, k). This routine tests an O(K)-local candidate list.

### 1.2 What is enumerated, in order

1. **Backbone projection** (frozen `observed_backbone`). For w outside the observed hull with
   gate g and off-distance a(w), d(r, w) = d(r, g) + a(w) for every r ∈ R, so hidden mass acts
   as mass at g plus a constant added to the pendant; mass at a hull vertex of hull-degree two
   splits exactly onto its chain endpoints. Hence (m, y) on T is feasible iff the push-forward
   (m̂, ŷ) is feasible on the backbone B, with ŷ = y + Σ_w m_w a(w) ≥ y > 0. B is reduced and
   its leaf set is exactly R, so the all-leaves theory applies.
2. **Deficit supports.** On B, every feasible potential is f = h + t with h the leaf envelope
   and t ≥ 0 vanishing at leaves. t_v > 0 iff every branch at v carries mass, i.e. iff every
   branch at v contains a source point. Each source point lies in exactly one branch of v, so
   deg(v) ≤ k; and with Y the minimal subtree of B spanning the k source points,
   deg_Y(v) = deg_B(v) at every deficit vertex and #leaves(Y) ≤ k, so
   **Σ_{t_v > 0} (deg_B(v) − 2) ≤ #leaves(Y) − 2 ≤ k − 2**. Admissible supports S are produced
   by the frozen generator `deficit_sets` under exactly that weight bound with 3 ≤ deg ≤ cap.
3. **Zero faces.** With S fixed the mass is affine in t and differs from the baseline only on
   A = S ∪ N(S). A positive pattern is a choice of which coordinates of A vanish. Vertices
   outside A keep their baseline sign. Inside A, for u ∈ N(S) \ S the increment
   Σ_{v ∈ S, v~u} t_v / (2 l_uv) is strictly positive on the strict face, so **m_u > 0 is forced
   whenever b_u ≥ 0**; only S itself and the strictly negative baselines of N(S) are free. This
   cuts the face enumeration from 2^|A| to 2^|free| and is a logical deduction, not an
   approximation — it is checked against the full 2^|A| enumeration in validation track B.
4. **Strict feasibility** of each face is decided exactly by the frozen Stage 43 common-slack
   vertex routine `strict_vertex` on the same constraint system the frozen solver builds.
5. **Minimum covers.** A feasible pattern P is realised by atoms: an edge atom on a backbone
   edge with both endpoints in P, or a vertex atom at a *marked* backbone vertex (a gate with
   unobserved material behind it). The minimum number of atoms is |P| − ν(G[P]) with ν the
   maximum matching of the forest induced on P; the routine enumerates **all** covers of that
   size by depth-first search with a ⌈uncovered/2⌉ bound, not one.
6. **Lifts.** An edge atom on backbone edge j sits at an unknown coordinate in (0, L_j), so its
   possible original edges are all original edges of chain j plus all omitted edges behind any
   chain-interior hull vertex of j; a vertex atom at gate g lifts to any omitted edge behind g.
   These candidate lists are pairwise disjoint across atoms, so the lifts of one configuration
   are their Cartesian product. The product is a **superset**; every element is then decided by
   the exact oracle, so no false positive survives and no true support is lost.

The minimum size itself is taken from the frozen `solve_partial` on every call, and the
structural enumeration is required to reproduce it; a mismatch is reported as
`status='disagreement'`, never repaired. Across 676 bank cases the two agreed 676/676.

### 1.3 Exhaustiveness — what it rests on, and what it does not

The argument is complete *given* three inputs:

- the weighted deficit bound and the positive-anchor growth implemented by the frozen
  `deficit_sets` generator, which this module **imports rather than re-derives** — the bound is
  restated above from the model equation, but the generator's completeness is inherited from
  stages 27/40/41 and is not re-proved here;
- the exactness of the frozen Stage 25 endpoint oracle;
- the lift candidate lists being supersets of the true lifts, which follows from the projection
  in step 1.

Given those, every original support of minimum size pushes forward to a backbone configuration
of exactly k\* atoms whose positive pattern is one of the enumerated feasible patterns, whose
atom set is one of the enumerated minimum covers, and whose edges lie in the enumerated
candidate lists. A resource stop is reported as a non-`complete` status and is **never** read as
a uniqueness verdict; in the bridge runs an indecisive probe aborts that regime as
`unknown_probe` rather than counting as a failure to be unique.

Two conditions are worth naming explicitly. First, the routine answers "all supports of size
**≤ cap**"; with `cap = K` it cannot distinguish "no explanation with ≤ K sources" from
"an explanation with more". Second, the lift product is capped (`max_oracle`); the cap was never
reached in any run reported here, but if it were, the status would say `oracle_limit`.

## 2. Validation against brute force — 5,839 comparisons, zero disagreements

Registered before execution in `results/stage292_protocol_v2.json`.

### Track A — replay of every Stage 290 probe point

Every probe recorded in `results/stage290_bridge.jsonl` (130 cases × 5 regimes × the probes of
each bisection), plus e = 0: **5,339 measurements**, each recomputed with the new routine and
re-run with the Stage 290 brute force.

| comparison | agreement |
|---|---|
| minimum source count, new vs Stage 290 brute force | **5,339 / 5,339** |
| complete minimiser **set**, new vs Stage 290 brute force | **5,339 / 5,339** |
| (minimum, count) vs the values recorded by Stage 290 | **5,209 / 5,209** |
| complete set vs prefilter-free exhaustive enumeration (|E| ≤ 70) | **2,421 / 2,421** |
| routine status | 5,339 `complete`, 0 resource stops |
| disagreements | **0** |

The agreement includes the hard cases: the largest minimiser set enumerated exactly contained
**4,443** supports, and the routine returned it identically to brute force. Median work per
probe: **3 exact oracle calls** against **20,706** brute-force edge subsets — a median of
3,980 subsets per oracle call. Median wall time 0.030 s against 0.43 s, a median speed-up of
**10.2×** (IQR 4.2–35.4×, maximum 1,486×). The speed-up is modest here *because these are
precisely the cases brute force could afford*; at K = 5 with |R| = 51 the brute force exceeds
its own 4 M-candidate limit and returns nothing at all.

### Track B — fresh instances against complete enumeration

500 freshly generated instances (7–12 leaves, K ∈ 1..4, random hidden-leaf subsets; 430 with
hidden leaves, 181 with more than one minimiser), each compared against the complete
C(|E|, k) enumeration decided by the same exact oracle:

| comparison | agreement |
|---|---|
| minimum source count | **500 / 500** |
| complete minimiser set | **500 / 500** |
| reduced zero-face enumeration vs full 2^\|A\| enumeration | **500 / 500** |

148,233 edge subsets were tested by brute force against 3,286 exact oracle calls by the routine.

After the one later code change (a bound on the distance-column memoisation, §6), 30 further
fresh instances were re-checked: 30/30 on all three comparisons.

## 3. Extended coverage in K

`scripts/stage292_bridge.py`, same definitions and conventions as Stage 290: t_unknown is the
largest e along a fixed removal order for which the set of **all** minimum-size supports is
exactly the single true support, found by **exponential bracketing then bisection with both
bracket endpoints explicitly evaluated**; MIN_REMAINING = 3 references are never erased;
regimes are random (seeds 11, 23), whole-clade (seeds 11, 23) and adversarial nearest-to-source
(deterministic).

| K | bank cases measured | of available | status |
|---|---|---|---|
| 2 | 120 | 120 | complete, now at full reference-set size (Stage 290 capped at 130) |
| 3 | 220 | 220 | complete, now at full reference-set size (Stage 290 capped at 45) |
| 5 | 220 | 220 | **complete — no Stage 290 measurement existed** |
| 7 | 116 | 119 | 3 cases hit the per-case wall; the batch otherwise ran to completion |
| 10 | 0 | 213 | **does not run**: 8/8 attempted cases hit a 1,200 s wall |
| 12 | 0 | 1 | not attempted after the K = 10 result |

676 cases in total, against 130 in Stage 290. 3,257 of 3,315 regime measurements resolved;
58 aborted honestly as `unknown_probe` (50 of them in the random regime, where the bracket runs
furthest).

**The honest ceiling at K ≥ 10 is in the frozen minimum pass, not in the new routine.** On the
*smallest* K = 10 bank case (|R| = 20, 37 edges) the frozen `solve_partial` returned
`status='unknown'` after its full 3,600 s budget (`results/stage292_k10_ceiling.txt`). The
all-minimisers extension cannot be faster than the minimum it is asked to enumerate at, so
K ≥ 10 is out of reach until the frozen deficit/zero-face search is improved — not because the
uniqueness question became harder.

## 4. Panel-size ladder on planted mixtures

`scripts/stage292_panels.py`. Planting rule, deterministic given (panel, size, K, seed) and
documented in the script header: the frozen Newick is parsed by the Stage 286 parser, the root
artefact is suppressed as Stage 290 does and every degree-two vertex is suppressed; for ladder
size L the leaves are drawn with `random.Random(10**6 + 1000·size_rank + seed)` and the induced
subtree is taken and renumbered; K edges are drawn uniformly from all edges with
`random.Random(seed)`; weights are w_i/Σw with w_i uniform on 1..9; positions are x_i = a_i/12
with a_i uniform on 1..11; the pendant is y = min_r (Σ_i p_i dist(r, q_i))/10. **These are
planted instances, not published mixtures**: they vary panel size at a fixed generative process
and carry no biological claim. Seeds 101, 202, 303. Panels: the 47,151-leaf WoL v2 placement
tree sub-sampled to 64…47,151, and the 15,953-leaf WoL v2 backbone sub-sampled to
100/1,000/15,953. 70 of 78 jobs returned; the true support is the unique minimiser at the full
panel in 70/70.

Median references erasable, pooling K = 3 and K = 5 and both panels:

| panel leaves | 64 | 128 | 256 | 512 | 1,024 | 2,048 | 4,096 | 8,192 | 15,953 | 16,384 | 47,151 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| nearest to source | 5 | 4 | 6 | 3 | 3 | 6 | 2 | 2.5 | 4 | 10 | 7 |
| uniformly at random | 24 | 59 | 131 | 297 | 793 | 748 | 2,863 | 3,957 | 9,665 | 13,756 | 15,629 |

Pooled Spearman of threshold against panel size: **0.045** adversarial, **0.981** random.
Log–log slopes: adversarial −0.11, −0.07, +0.02, +0.05 across the four panel × K series;
random +0.97, +1.06, +1.02, +0.96. As a *fraction* of the panel the adversarial threshold
collapses from 7.8 % at 64 leaves to **0.01 %** (median 7 of 47,151 references) at the top of
the ladder, while the random fraction stays between 33 % and 84 % throughout.

## 5. Revised threshold table (frozen bank, 663 cases whose true support is the unique minimiser at the full panel)

Median t_unknown (IQR), and Spearman correlation with reference-set size:

| K | nearest to source | whole clade | uniformly at random |
|---|---|---|---|
| 2 | **2** (1–4), ρ = −0.04 | 37.5 (17–74), ρ = 0.65 | 74 (30–122), ρ = 0.85 |
| 3 | **2** (1–4), ρ = −0.01 | 20 (8–42), ρ = 0.60 | 33 (19–70), ρ = 0.82 |
| 5 | **2** (1–4), ρ = 0.20 | 12 (5–31), ρ = 0.67 | 25 (10–52), ρ = 0.86 |
| 7 | **2** (1–4), ρ = 0.12 | 16 (5–37), ρ = 0.73 | 32 (10–59), ρ = 0.78 |
| pooled | ρ = **0.095** | ρ = 0.647 | ρ = 0.832 |

As a fraction of the panel, pooled over K: **2.1 %** adversarial, 22.0 % whole-clade, 42.2 %
random (Stage 290 reported 4 %, 36 %, 52 % on its K = 2/3 subset).

Known-support comparison, carried over unchanged from Stage 290's design:

- the cogirth bound is never violated — t_unknown ≥ cogirth − 1 in **3,257 / 3,257** measurements —
  and remains nearly vacuous;
- the fibre-dimension rank identity rank = s + 1 held in **676 / 676** cases;
- the sequence-specific known-support quantity t_known predicts t_unknown well under random
  (ρ = 0.74) and whole-clade (ρ = 0.90) loss but **not** under adversarial loss (ρ = 0.15).

13 of 676 bank cases are already non-identifiable at the full panel (the true support is not the
unique minimiser); in 2 of those the minimum is strictly below K. They are excluded from the
threshold table and recorded in the raw files.

## 6. What is now established, and what remains provisional

**Established (was provisional after Stage 290).**

1. Deciding whether the minimum-source explanation is *unique* does not require brute force. The
   structural routine returns the complete minimiser set exactly and agreed with brute force on
   5,839/5,839 comparisons with zero disagreements, including a case with 4,443 minimisers.
2. The adversarial erasure threshold does not improve with panel size. Established on 676 bank
   cases at K ∈ {2, 3, 5, 7} (ρ = 0.095) and on a planted ladder spanning 64 → 47,151 leaves
   (ρ = 0.045, log–log slope ≈ 0), against ρ = 0.83 / 0.98 and slope 1.0 for random loss.
3. The threshold is flat in K as well: median 2 references at every K measured.
4. The contrast between the three regimes is not an artefact of K = 2/3 or of |R| ≤ 130.

**Still provisional or open.**

- **K ≥ 10 is unmeasured.** 8/8 attempted cases hit the wall and the smallest K = 10 case defeats
  the frozen minimum pass at 3,600 s. Nothing in this stage supports extrapolating the flatness
  claim beyond K = 7.
- **K = 7 coverage is 116/119**, not complete; the 3 missing cases are among the largest.
- **The ladder is planted, not published.** The bank supplies published mixtures only up to 361
  references; the 10³–10⁴-leaf evidence comes from the planting rule above. A reader who rejects
  the planting rule is left with the bank result (ρ = 0.095 up to 361 references).
- **Exhaustiveness is inherited, not re-proved.** The deficit-support generator is the frozen
  one; this stage restates its bound and tests its consequences against brute force on 5,839
  measurements, which is strong evidence but not a proof of the generator.
- **58 of 3,315 regime measurements are `unknown_probe`** and contribute nothing either way.
- Two panel jobs (1,024/K=5/seed 101 and 2,048/K=5/seed 202) timed out twice and are missing from
  the ladder; six further jobs were still running when the batch was closed and are recorded as
  such.

## 7. Files, protocol and reproduction

All paths are under `~/projects/recomb2027_novelty` on `ssh:blip`. No pre-existing file was
modified.

**Protocol** (registered before the measurements it covers):
`results/stage292_protocol.json` → superseded, before any measurement, by
`results/stage292_protocol_v2.json`; `results/stage292_protocol_v2_addendum.json` covers the
K = 7 resume batch and the K ≥ 10 batch; `results/stage292_protocol_v2_addendum2.json` covers the
panel re-run and records the one code change made after measurements had begun (a bound on the
exact-distance memoisation table, which cannot change a returned value, verified by re-running
30 track-B instances).

**Code (new):** `scripts/stage292_all_minimisers.py` (the routine),
`scripts/stage292_harness.py`, `scripts/stage292_orders.py`, `scripts/stage292_validate.py`,
`scripts/stage292_bridge.py`, `scripts/stage292_panels.py`,
`scripts/stage292_panels_rerun.py`, `scripts/stage292_analysis.py`,
`scripts/stage292_k10_ceiling.py`, `scripts/stage292_pilot.py`.

**Raw records:** `results/stage292_validate_replay.jsonl` (5,339 rows),
`results/stage292_validate_fresh.jsonl` (500), `results/stage292_bridge_K2_3_5_7.jsonl` (599
cases), `results/stage292_bridge_K7_resume.jsonl` (80 jobs, 77 cases + 3 timeouts), `results/stage292_bridge_K10_12.jsonl` (8
timeouts, retained), `results/stage292_panels.jsonl`, `results/stage292_panels_rerun.jsonl`,
`results/stage292_k10_ceiling.txt`.

**Analysis:** `results/stage292_analysis.json`.

**Removal orders** are Stage 290's. The adversarial order is recomputed from a handful of BFS
distance columns instead of the dense |R| × |V| float matrix so that it also works on 10⁴-leaf
panels; it was verified identical to Stage 290's order on the cases where both are affordable,
and the whole bridge reproduces Stage 290's t_unknown and t_known values exactly on the spot
checks run before the batch.

## 8. Declared reductions in scope

- K = 10 and K = 12 were attempted on a size-stratified 1-in-3 subsample (72 of 214 cases) and
  **produced no measurements**: the 8 cases reached before the batch was stopped all hit the
  1,200 s per-case wall. No K ≥ 10 number is reported.
- K = 7 was measured on **116 of 119** cases (3 wall timeouts). The missing cases
  are among the largest reference sets, so the K = 7 Spearman against |R| is the weakest entry in
  the table.
- The panel ladder returned **70 of 78** jobs; 2 specs timed out twice and 6 were still running
  when the batch was closed. Coverage at both ends of the ladder (64 and 47,151 leaves, both K)
  is complete.
- Panel-ladder mixtures are **planted**, with the rule and seeds given in §4; they are not
  published mixtures.
- The exhaustive, prefilter-free brute force in track A was run only on cases with |E| ≤ 70
  (2,421 of 5,339 probe points); the remaining 2,918 were compared against Stage 290's
  singular-value-prefiltered brute force, which is the same enumeration with a float pre-screen
  whose survivors are confirmed exactly.

