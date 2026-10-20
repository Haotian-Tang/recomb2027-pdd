# Stage 301 - transporting the structural localisation to the tolerant problem

## Verdict

The conjecture holds. Every ingredient of the frozen localisation (`stage43_partial_fpt.solve_mixed`) survives the move from exact equality to a tolerance, but **not** by substituting `d -> d - eps`. A uniform shift of the observations leaves the baseline masses unchanged and moves only the baseline pendant (`inverse(h + c) = (b, Y + c)`), so the substitution the stage brief proposed is by itself inert. What transports is the DECOMPOSITION behind the localisation, with the sign sets replaced by threshold sets.

Let `f(v) = y + sum_w m(w) d(v,w)` be the potential of a candidate configuration (`m >= 0`, `sum m = 1`, mass only at the atom vertices of the support), `h(v) = max_r (d_r - d(v,r))` the frozen leaf envelope, `b = inverse(h)` the baseline mass and `phi = f - h`.

## Per-ingredient soundness

| # | Frozen ingredient | Survives? | Tolerant form |
|---|---|---|---|
| I1 | leaf-envelope consistency `h(r) == d_r`, else reject | relaxed | `h(r) - d_r <= 2 eps`; the equality test is an upper-bound test and dies, the inequality is necessary |
| I2 | baseline pendant `Y > 0`, else reject | relaxed, vacuous | `Y >= -eps (|R| - 1)`; the coefficient sum `sum_v |2 - deg(v)| = 2|R| - 2` makes this empty at any useful eps |
| I3 | mass-count obstruction `|positive| <= 3k-2`, `|negative| <= 4k-8` | YES, with threshold sets | same counts applied to `POS(eps)` / `NEG(eps)` below |
| I4 | weighted deficit bound `sum_{v in Sigma} (deg(v)-2) <= k-2`, `3 <= deg(v) <= k` | YES, UNCHANGED | purely structural: it constrains the mass configuration and never touches the observations |
| I5 | `negative` must lie inside `N[Sigma]` | YES, with threshold sets | `NEG(eps)` inside `N[Sigma_eps]` |
| I6 | positive-anchor filter: every component of `T[Sigma]` holds a vertex with `b(v) > 0` | YES, UNCHANGED and NOT relaxable | the proof returns `b(v*) > 0` strictly; widening the anchor set to `b(v) > -tau` would be unsound |
| I7 | atom cover of the positive-mass set by at most k atoms | YES, with threshold sets | the atoms of the support must cover `POS(eps) \ Sigma_eps` |
| I8 | the local strict-slack LP in `(delta, t)` | NO | under a tolerance `phi` is no longer supported on `Sigma`, so the dimension reduction from `|V|` to `|Sigma| + 1` is lost; every surviving support is handed to the frozen stage-296 LP instead |

### The four facts the transport rests on

1. `|f(u) - f(v)| <= d(u,v)` because the masses sum to one.
2. `theta(S) <= eps` implies `phi >= -eps` everywhere (`f(v) >= f(r) - d(v,r) >= d_r - eps - d(v,r)`).
3. A vertex with a mass-free branch is EXPOSED and then `f(v) = f(r) - d(v,r)` exactly for an observed `r` in that branch, so `phi(v) <= eps`. Observed leaves satisfy `phi(r) <= eps` directly. Hence `Sigma_eps := {v : phi(v) > eps}` contains only surrounded internal vertices; the counting argument on the Steiner tree of the mass then gives I4 verbatim, because it never mentions `d`.
4. Off `Sigma_eps` the increment `phi(u) - phi(v) = (f(u) - f(v)) - (h(u) - h(v))` is squeezed by both the tolerance and the Lipschitz bound, giving the one-sided bands

       m(v) - b(v) <=  tau_up(v) = sum_{u~v} min(eps/l, (l - dh)/(2l))
       m(v) - b(v) >= -tau_dn(v) = -sum_{u~v} min(eps/l, (l + dh)/(2l))

   with `dh = h(u) - h(v)`, hence the threshold sets

       POS(eps) = {v : b(v) >  tau_dn(v)}   (forced into supp(m))
       NEG(eps) = {v : b(v) + tau_up(v) < 0}   (infeasible outside N[Sigma])

   The eps-free limits of the two band sums are `1 - b(v)` and `deg(v) - 1 + b(v)`, which reproduce exactly `0 <= m(v) <= 1`. That is the boundary of the method, quantified in the next section.

### What does NOT survive, with the reason

* **The equality test `h(r) == d_r`.** It is an upper bound on the envelope and has no tolerant analogue beyond `h(r) - d_r <= 2 eps`.
* **`Y > 0`.** The relaxation is `Y >= -eps(|R| - 1)`; on the travellers' backbone (`|R| = 10,575`) that is empty for every eps of interest. Reported as vacuous, not as a filter.
* **The reduced local LP (I8).** The frozen solver parametrises `f = h + delta` with `delta` supported on `Sigma`. Under a tolerance `phi` is nonzero at every vertex (bounded by eps, but nonzero), so the reduction is false. This is the one ingredient that genuinely dies, and it is why stage 301 screens supports but still pays the frozen LP per survivor.

## Validation (exact rational equality, non-negotiable)

| Arm | What it checks | Result |
|---|---|---|
| A | localised `Theta_k` at `eps = Theta_k` vs the stage-296 exhaustive value | **2284 / 2284 exact agreement, 0 disagreements** |
| A | of those, runs where the screen actually pruned | 1183 |
| A2 | at `eps = Theta_k / 2` the routine must never return a value `<= eps` | 2034 non-degenerate runs, **0 violations** (273 rows excluded: `Theta_k = 0` makes `eps = 0` not strictly below) |
| A2 | of those, refuted outright with no LP at all | 317 |
| C | with the screen DISABLED the enumerator visits exactly `C(|E|, k)` | 40 cases, **0 violations** |
| - | attaining witnesses re-verified with the frozen `min(d(r,u)+xL, d(r,v)+(1-x)L)` form | 2284 checked, **0 mismatches** |

Real-data soundness, where stage 298 already knows `Theta_k` exactly:

* travellers `eps_R(1) <= Theta_1`: 29 / 29, **0 violations**
* real ANI and real placement `eps_R(k) <= Theta_k`: 978 checks, **0 violations**

## Measured pruning and the vacuity boundary

On the Theta ladder (`Theta_1` exhaustive, then `Theta_k` screened at `eps = Theta_{k-1}`), over 146 exact-model cases and 185 rungs with k >= 2:

At `eps = Theta_k` across the 2284 validation runs the screen was non-vacuous on 577 of 703 k = 2 runs (median ratio 0.0230, best 3.35e-05) and on 10 of 10 k = 4 runs (median 0.0051, best 5.61e-06); at k = 1 it is vacuous on 969 of 1555.

* non-vacuous on **63 of 185** rungs; median candidates / `C(|E|,k)` = **0.0225** (about 44x), best **2.18e-05**
* vacuous on **122** rungs - reported as vacuous, not as a success
* seconds per support: median **0.0070 s**, range 0.0024-0.0792 s (stage 296 measured 0.0054-0.0089 s, i.e. the same scale)

**The boundary, and it is structural.** `b(v)` lies in `[1 - deg(v), 1]` and `tau` is capped by `deg(v)`, so once eps exceeds the incident edge lengths `POS(eps)` and `NEG(eps)` are provably empty and the screen is vacuous by construction. The tolerances at which `Theta_k` is attained on real estimated distances sit above that threshold, which is why the localisation does not unlock `Theta_2` on the travellers' samples.

## Measured reach

| k | stage 296 frontier (|R|) | stage 301 (|R| reached) |
|---|---|---|
| 1 | ~361 | **361** |
| 2 | ~110 | **361** |
| 3 | ~46 | **100** |
| 4 | ~27 | **43** |
| 5 | not reached | **41** |
| 6 | not reached | **23** |
| 7 | not reached | **23** |

Budget per case: single core, ladder wall 1500 s, candidate cap 60,000 - SMALLER than the ~1 h single core that defined the stage-296 frontier. k = 5, 6, 7 were not reachable at all in stage 296.

## eps x incompleteness above K = 3

144 draws on exact-model cases with K >= 5 under the stage-288/299 removal regimes; the ladder reached k = 8, and 105 draws reached k >= 4 - stage 296's grid stops at K in {2, 3}.

Across every cell: **502** paired comparisons, **411** falls, **0** rises. Stage 296 proved Theta_k cannot rise when references are removed (removal deletes LP rows); this is the first confirmation above k = 3.

Cells with k >= 3, the territory stage 296 could not enter:

| regime | level | k | n | median Theta_k(draw) / Theta_k(full) | falls | rises |
|---|---|---|---|---|---|---|
| clade | 1/2 | 3 | 18 | 0.196 | 17 | 0 |
| clade | 1/2 | 4 | 8 | 0.000 | 8 | 0 |
| clade | 1/2 | 5 | 7 | 0.000 | 3 | 0 |
| clade | 1/2 | 6 | 1 | 0.000 | 1 | 0 |
| clade | 1/4 | 3 | 18 | 0.628 | 14 | 0 |
| clade | 1/4 | 4 | 7 | 0.942 | 5 | 0 |
| clade | 1/4 | 5 | 7 | 0.000 | 3 | 0 |
| clade | 1/4 | 6 | 2 | 0.000 | 2 | 0 |
| clade | 3/4 | 3 | 18 | 0.000 | 18 | 0 |
| clade | 3/4 | 4 | 7 | 0.000 | 7 | 0 |
| clade | 3/4 | 5 | 3 | 0.000 | 3 | 0 |
| nearest | 1/2 | 3 | 9 | 0.000 | 9 | 0 |
| nearest | 1/2 | 4 | 3 | 0.000 | 3 | 0 |
| nearest | 1/2 | 5 | 2 | 0.000 | 2 | 0 |
| nearest | 1/4 | 3 | 9 | 0.928 | 7 | 0 |
| nearest | 1/4 | 4 | 3 | 0.633 | 3 | 0 |
| nearest | 1/4 | 5 | 1 | 0.469 | 1 | 0 |
| nearest | 3/4 | 3 | 9 | 0.000 | 9 | 0 |
| nearest | 3/4 | 4 | 2 | 0.000 | 2 | 0 |
| random | 1/2 | 3 | 18 | 0.621 | 15 | 0 |
| random | 1/2 | 4 | 8 | 0.535 | 7 | 0 |
| random | 1/2 | 5 | 5 | 0.000 | 3 | 0 |
| random | 1/2 | 6 | 1 | 0.000 | 1 | 0 |
| random | 1/4 | 3 | 18 | 1.000 | 8 | 0 |
| random | 1/4 | 4 | 8 | 1.000 | 4 | 0 |
| random | 1/4 | 5 | 7 | 0.791 | 3 | 0 |
| random | 1/4 | 6 | 1 | 0.000 | 1 | 0 |
| random | 3/4 | 3 | 18 | 0.000 | 18 | 0 |
| random | 3/4 | 4 | 7 | 0.000 | 7 | 0 |
| random | 3/4 | 5 | 3 | 0.000 | 3 | 0 |

The regime ordering of stages 288/299 survives above k = 3: at 75 % loss and k = 1 the median ratio is 0.066 for whole-clade and 0.144 for nearest-reference against 0.654 for random, and structured loss keeps collapsing Theta_k to exactly 0 at every k it reaches.

## Theta_2 on the 29 travellers'-diarrhoea samples

**Not delivered exactly, and the reason is measured.** What stage 301 delivers instead is a certified two-sided bracket for every k, from one envelope pass and no support enumeration:

    eps_R(k)  <=  Theta_k  <=  Theta_1   (stage 298)

with `eps_R(k) = max( (4k-7)-th largest critical eps of NEG, (3k-1)-th largest critical eps of POS )`, the transported count obstruction. Medians over the 29 samples, in the samples' own units:

| k | median `eps_R(k)` | / `Theta_1` | / median observed distance | reading |
|---|---|---|---|---|
| 1 | 0.5910 | 0.0175 | 1.99e-03 | for eps <= `eps_R(1)` no 1-source explanation exists, so m(eps) >= 2 |
| 2 | 0.5294 | 0.0159 | 1.77e-03 | for eps <= `eps_R(2)` no 2-source explanation exists, so m(eps) >= 3 |
| 3 | 0.2787 | 0.0075 | 9.16e-04 | for eps <= `eps_R(3)` no 3-source explanation exists, so m(eps) >= 4 |
| 4 | 0.1633 | 0.0049 | 5.38e-04 | for eps <= `eps_R(4)` no 4-source explanation exists, so m(eps) >= 5 |
| 5 | 0.1189 | 0.0032 | 4.03e-04 | for eps <= `eps_R(5)` no 5-source explanation exists, so m(eps) >= 6 |
| 6 | 0.0911 | 0.0025 | 2.85e-04 | for eps <= `eps_R(6)` no 6-source explanation exists, so m(eps) >= 7 |
| 7 | 0.0719 | 0.0018 | 2.27e-04 | for eps <= `eps_R(7)` no 7-source explanation exists, so m(eps) >= 8 |
| 8 | 0.0543 | 0.0016 | 1.75e-04 | for eps <= `eps_R(8)` no 8-source explanation exists, so m(eps) >= 9 |

`Theta_1` itself is 0.0399-0.2052 of the median observed distance (median 0.1058) on these samples, so the bracket is wide: the certificate reaches a median 0.0175 of `Theta_1`.

On the real collections where `Theta_k` IS known exactly, the same certificate is much tighter - median `eps_R(k) / Theta_k` of 0.336, 0.413, 0.332 for k = 1, 2, 3 over 978 checks - so the looseness is a large-panel effect, not a defect of the bound.

### Why the support-free relaxation cannot close the gap

Dropping the support-size restriction entirely gives `Theta_inf = min over all mass distributions`, a valid lower bound on every `Theta_k`. Column generation with exact O(|V|) pricing drives the master to a median of **0.0019 of `Theta_1`** (range 0.0002-0.0107) on the 29 samples: an unrestricted mixture fits these data essentially exactly, so `Theta_inf` is worthless as a lower bound on `Theta_2` and the count obstruction (`eps_R(2)`, median 0.0159 of `Theta_1`) is the stronger device by an order of magnitude. The support-size restriction is doing all of the work.

## What this does NOT establish

1. It does not give `Theta_2` or `Theta_3` on the 29 travellers' samples. Only the bracket above is certified.
2. The certified minimum counts hold at tolerances of order 1e-4 to 2e-3 of the median observed distance, while stage 298 measured real estimator residuals at 0.15-0.34 of it. The counts are therefore certified but sit far below any plausible error floor; they are NOT a claim about how many organisms are in a stool sample.
3. The travellers' samples have no ground truth and no measurable estimator residual, so nothing here is validated against a true source count.
4. The screen is a NECESSARY condition only. When the minimum over survivors exceeds eps the run certifies `Theta_k > eps` and nothing more.
5. The vacuity boundary is measured on these trees, not proved to be tight; a sharper band than `tau_up / tau_dn` may exist.
6. I8 is lost, so the per-survivor cost is still a full stage-296 LP; the speed-up is entirely in how many supports are visited.
7. The k >= 2 validation sample is bounded by tree size (see deviations in the protocol): all k = 1 rows plus k = 2 rows with |E| <= 160 and a fixed random sample of 40 larger ones.

## Record hygiene and declared reductions

* `stage301_validate_runs_full.jsonl` is a PARTIAL first validation run (1,063 records) that was stopped because vacuous-screen k = 2 rows fall back to a full `C(|E|,2)` sweep and the breadth would have taken ~16 h. It is retained; `stage301_validate_runs_main.jsonl` is the bounded re-run and both are read by the analysis.
* `stage301_grid_runs.jsonl` (v1) is retained but SUPERSEDED: `clade_removals` returns dicts and the first driver treated them as leaf lists, so the whole-clade draws removed nothing (`n_retained == n_full`). `stage301_grid_runs_v2.jsonl` is the corrected run and is the one analysed.
* `stage301_dual_runs.jsonl` is the first Theta_inf attempt with a row-restricted master; it is sound but gave a zero bound on the travellers because a 40-150 leaf subset is fitted exactly by free mass. `stage301_thetainf_runs.jsonl` is the column-generation version.
* Validation breadth: all k = 1 rows of stage 296 plus k = 2 rows with |E| <= 160 and a fixed random sample of 40 larger ones (seed 301).
* Reach and grid use `per_cell = 6` and `per_cell = 3` bank cases per cell, 2 replicates per removal level, |R| <= 60 for the grid.
* Everything ran single-core per job with at most 40 concurrent workers.

## Files

* `scripts/stage301_common.py` - transported screens, exact critical-eps machinery, candidate planner
* `scripts/stage301_protocol.py`, `scripts/stage301_protocol_amendment.py`
* `scripts/stage301_validate.py` - arms A, A2, B, C
* `scripts/stage301_prune.py` - the Theta ladder, reach and pruning
* `scripts/stage301_grid.py` - eps x incompleteness above K = 3
* `scripts/stage301_td.py` - certified brackets on real data
* `scripts/stage301_thetainf.py` - Theta_inf column generation
* `scripts/stage301_analyze.py`, `scripts/stage301_report.py`
* `results/stage301/` - protocol, amendment, every run record, `stage301_analysis.json`

