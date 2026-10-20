# Stage 289 — Head-to-head with the published DecoDiPhy implementation

Phase 2, track "Comparison". All artefacts live on `ssh:blip` under
`~/projects/recomb2027_novelty`, stage number 289. No pre-existing file was modified.

---

## 1. The tool that was run

| field | value |
|---|---|
| tool | DecoDiPhy 1.5.0 |
| repository | `https://github.com/shayesteh99/DecoDiPhy` |
| commit | `a26aa618aaf50eaba80d12bb3e3d71a201d8617e` |
| `git describe --tags` | `v1.5.0-11-ga26aa61` (11 commits after the v1.5.0 tag). `git diff v1.5.0..HEAD` touches five files: the README, the Colab tutorial, an R plotting script, a deleted PDF, and one bug fix in `decodiphy/cli.py` — three call sites of `compute_individual_ys` in the anchor post-processing now pass node labels instead of indices. Every other module, including `search.py`, `optimize.py` and `matrix.py`, is unchanged from the tag. HEAD was used, i.e. the latest published code including that fix. |
| commit date | 2026-07-10 |
| install | `python3 -m venv` + `pip install .` from the clone, in a dedicated environment (`data/stage289_decodiphy/venv`); the project's shared `.venv` was not touched |
| interpreter | Python 3.12.3 (the README suggests conda with Python 3.11; no conda and no 3.11 interpreter exists on this host, so 3.12.3 was used) |
| resolved dependencies | `numpy==2.5.4`, `scipy==1.18.1`, `cvxpy==1.9.3`, `osqp==1.1.3`, `treeswift==1.1.51` |
| host | Linux 6.8.0, x86_64, single-threaded (`OMP/OPENBLAS/MKL/VECLIB/NUMEXPR_NUM_THREADS=1`) |

The clone was checked against two copies of the same sources archived earlier in this project by a
different route (`data/stage09_prior/v1.5.0__decodiphy__*.py`, `data/stage07_distance_pilot/decodiphy__*.py`):
`search.py`, `optimize.py` and `matrix.py` are byte-identical (SHA-256 match) across all three, and the
tagged `v1.5.0:decodiphy/search.py` has the same digest, so the clone is the published code and not a
local variant.

Two properties of the released version constrain what can be compared, and were not changed:

- `search.py` v1.5.0 raises `NotImplementedError` for `--method exhaustive`, `closest`,
  `closest_iterative` and `reverse_hill`, with the message that they require a legacy cvxpy solver that
  is no longer available. Only the default `--method hill` is operable.
- `optimize.py` never imports cvxpy, so `--optimizer cvxpy` and the `--score` diagnostic raise
  `NameError: name 'cp' is not defined`. Only the default `--optimizer osqp` is operable.

Everything below therefore uses the tool's documented defaults: `--method hill`, `--quick 1`,
`--radius 2`, `--seed 1142`, `--optimizer osqp`, `--min_p 0.01` (the default when the input is a
distance file). Nothing was tuned in either direction.

---

## 2. Reproduction check

**What was reproduced.** The example distributed with the tool. The repository ships
`Example/krepp_multiplacement.jplace` (70 MB of krepp read placements) together with the author's own
output in `Example/output/` — the committed `output.jplace` records the invocation
`decodiphy -j Example/krepp_multiplacement.jplace -o Example/output` run on the author's machine. The
same command was re-run here into a fresh directory and compared field by field
(`scripts/stage289_repro_check.py`, report `results/stage289/stage289_repro_jplace.json`).

**What "reproduces" means here, and how close it was.**

| quantity | agreement |
|---|---|
| number of search rounds | 13 published, 13 produced (k = 1…12, then the final record) |
| anchor labels, every round | identical, and in identical order, in all 13 rounds |
| final placement count | 11 published, 11 produced |
| final placement edges (jplace `edge_num`) | identical for all 11 |
| final abundances | max absolute difference 6.5 × 10⁻¹³ |
| final relative positions `x` | max absolute difference 4.4 × 10⁻¹¹ |
| final terminal lengths `y` | max absolute difference 1.8 × 10⁻¹³ |
| labelled backbone tree | byte-identical after whitespace strip |
| loss, final round | relative difference 3.5 × 10⁻⁵ |

The one visible difference is confined to the k = 12 round, which both runs discard: its hill-climbing
optimum differs by 4.5 × 10⁻² in one `x` coordinate and 3.5 × 10⁻⁵ relatively in loss. Both runs stop
at the same place for the same reason ("Stopped by low abundance!") and emit the same 11 placements.
This is an installation that reproduces the parent work's own published example to floating-point
tolerance, with identical discrete output.

---

## 3. Comparison protocol

Registered before any solve as `results/stage289/stage289_headtohead_protocol.json`, carrying the
environment record above, the per-case SHA-256 of every input, the SHA-256 of every frozen solver
source and of the four new stage-289 sources, and the matched limits. The frozen bank is referenced by
hash (`bank_index_sha256 bf53a419…`, `bank_protocol_sha256 1aebf33c…`) and every case file's digest was
re-checked at run time.

**Instances.** All 1,188 bank cases of the three families the task names:

| family | cases | reference leaves | generating K | distances |
|---|---|---|---|---|
| `exact_10k` | 296 | 20–129 | 3, 5, 10 | exact, tree-derived |
| `exact_biotrees` | 596 | 20–361 | 2, 3, 5, 7, 10 | exact, tree-derived |
| `krepp_estimated_10k` | 296 | 20–128 | 3, 5, 10 | real krepp ANI estimates |

**Both tools see the same numbers.** The bank stores exact rational edge lengths and observed
distances. DecoDiPhy receives them as a Newick tree and a two-column distance file written with
shortest round-trip doubles; the exact solver reads the unrounded rationals. The encoding was validated
end to end: on a case where DecoDiPhy's returned support equals the generating support, the exact
rational solve on that support reproduces DecoDiPhy's reported abundances and positions to ~10⁻⁶, and
recovers the bank's generating `p`, `x` and `y` exactly.

**Handling the k asymmetry.** DecoDiPhy takes the number of placements as an input or chooses it by its
own stopping rules; the exact solver returns the minimum. Three DecoDiPhy arms make that explicit, all
on documented defaults except for the named k control:

| arm | command | what it answers |
|---|---|---|
| `dd_default` | `decodiphy -t T -d D -o OUT` | what a user gets with no k knowledge at all |
| `dd_true_k` | `… -k K_true -f 1` | what it gets when an oracle supplies the generating k |
| `dd_resolution` | `… --min_p 1e-12 -k K_true+3 -f 0` | its own near-zero-residue stop with the abundance stop effectively disabled, and room to exceed the true k |

The exact arm is `stage43_partial_fpt.solve_partial(edges, observed, cap=case cap)` with
`seconds/max_states/max_bases/max_faces = None`, i.e. byte-identical options to the stage-44
confirmation protocol.

**Execution.** The stage-44 house pattern: every job in its own forked child with
`RLIMIT_AS = 8 GiB` and an external hard wall of **120 s**, identical for both tools; two fresh
repetitions per (case, method) in a rotated, alternating order; every run retained including timeouts;
results streamed to `results/stage289/stage289_headtohead_runs.jsonl` (9,504 records). 12 worker
processes. Input files are written once per case outside the timed region, and `decodiphy.cli` is
imported in the pool worker before forking, so neither tool is charged file I/O or import time.

**Audit before interpretation** (`results/stage289/stage289_verification.json`, status **PASS**):

- 675 exact-solver witnesses re-verified with the frozen stage-42 verifier — 0 failures;
- 2,295 DecoDiPhy supports certified feasible were re-checked by evaluating the certificate through
  the frozen stage-26 observation routine — 0 failures;
- 0 contradictions (no certified-feasible support smaller than a completed exact minimum);
- the two repetitions agreed on status, source count and support for every (case, method) pair.

**One derived quantity needs defining.** For a fixed set of distinct edges, the PDD model is affine in
`(p_i, t_i = p_i x_i, y)`, so whether that support admits an *exact* explanation is rational linear
feasibility. `scripts/stage289_support_feasibility.py` decides it in `Fraction` arithmetic: it reduces
the equality system exactly and either returns the unique `(p, t, y)` and checks `p > 0`, `0 < t < p`,
`y > 0`, or — when the support is not identifiable from the observations — searches the solution
subspace for a strictly feasible point and certifies it exactly. This turns "DecoDiPhy's loss is
10⁻¹⁷" into a yes/no question: *is the returned support one that could possibly be right?*

---

## 4. Results

### 4.1 Completion and cost

| family              | method                | cases | complete | timeout | error | median CPU s | total CPU s |
|---------------------|-----------------------|-------|----------|---------|-------|--------------|-------------|
| exact_10k           | exact_stage43_partial | 296   | 200      | 96      | 0     | 0.023        | 11536.1182  |
| exact_10k           | dd_default            | 296   | 296      | 0       | 0     | 1.554        | 757.6926    |
| exact_10k           | dd_true_k             | 296   | 296      | 0       | 0     | 1.563        | 858.7838    |
| exact_10k           | dd_resolution         | 296   | 296      | 0       | 0     | 1.568        | 822.6662    |
| exact_biotrees      | exact_stage43_partial | 596   | 475      | 121     | 0     | 0.050        | 15993.9890  |
| exact_biotrees      | dd_default            | 596   | 596      | 0       | 0     | 3.138        | 2888.2689   |
| exact_biotrees      | dd_true_k             | 596   | 596      | 0       | 0     | 3.224        | 2951.9050   |
| exact_biotrees      | dd_resolution         | 596   | 596      | 0       | 0     | 3.230        | 2977.0703   |
| krepp_estimated_10k | exact_stage43_partial | 296   | 293      | 3       | 0     | 0.009        | 363.3036    |
| krepp_estimated_10k | dd_default            | 296   | 290      | 0       | 6     | 3.849        | 1551.3806   |
| krepp_estimated_10k | dd_true_k             | 296   | 285      | 0       | 11    | 1.634        | 809.2068    |
| krepp_estimated_10k | dd_resolution         | 296   | 289      | 0       | 7     | 2.424        | 986.0584    |

On the 675 exactly-explainable cases where the exact solver completed, the median per-case ratio of
DecoDiPhy CPU time to exact-solver CPU time is **45.6×** (`dd_true_k`; 45.8× for `dd_default`,
45.7× for `dd_resolution`). The exact method is not the slow one in that regime — it is roughly two
orders of magnitude faster.

Where it stops is sharp and depends on the source count, not on the reference-set size:

| generating K | cases | exact solver completes within 120 s | median exact CPU s (all runs) |
|---|---|---|---|
| 2 | 120 | 100 % | 0.024 |
| 3 | 220 | 100 % | 0.020 |
| 5 | 220 | 100 % | 0.058 |
| 7 | 119 | 91.6 % | 3.10 |
| 10 | 213 | 2.8 % | 120 (budget) |

A separate probe at a 1,800 s wall and a 1,750 s internal budget
(`results/stage289/stage289_exact_cost_probe.json`, 18 cases, 1 repetition) shows the 120 s budget is
not the binding constraint: all 12 sampled K = 10 cases exhaust 1,750 s (reference sets of 20 to 231
leaves alike), while all 6 sampled K = 7 cases finish in 0.1–11.3 s. At exhaustion 5 of the 12 K = 10
runs had an incumbent feasible solution but no proof of minimality. The wall is the source count.

### 4.2 Accuracy on exactly explainable data

| family              | method                | complete | support = truth | true edges recovered | support exactly feasible | support provably infeasible | other | median W1 | median loss |
|---------------------|-----------------------|----------|-----------------|----------------------|--------------------------|-----------------------------|-------|-----------|-------------|
| exact_10k           | dd_default            | 296      | 190             | 0.874                | 210                      | 83                          | 3     | 2.01e-05  | 3.56e-17    |
| exact_10k           | dd_true_k             | 296      | 228             | 0.943                | 256                      | 29                          | 11    | 1.66e-07  | 2.65e-17    |
| exact_10k           | dd_resolution         | 296      | 210             | 0.941                | 244                      | 42                          | 10    | 4.81e-08  | 3.41e-17    |
| exact_10k           | exact_stage43_partial | 200      | 196             | -                    | 200                      | 0                           | -     | -         | -           |
| exact_biotrees      | dd_default            | 596      | 433             | 0.917                | 446                      | 147                         | 3     | 1.54e-08  | 6.95e-17    |
| exact_biotrees      | dd_true_k             | 596      | 555             | 0.983                | 569                      | 19                          | 8     | 1.06e-08  | 5.32e-17    |
| exact_biotrees      | dd_resolution         | 596      | 547             | 0.986                | 570                      | 17                          | 9     | 1.06e-08  | 5.54e-17    |
| exact_biotrees      | exact_stage43_partial | 475      | 473             | -                    | 475                      | 0                           | -     | -         | -           |
| krepp_estimated_10k | dd_default            | 290      | 0               | 0.238                | 0                        | 284                         | 6     | 0.084     | 0.006       |
| krepp_estimated_10k | dd_true_k             | 285      | 0               | 0.200                | 0                        | 285                         | 0     | 0.076     | 0.006       |
| krepp_estimated_10k | dd_resolution         | 289      | 0               | 0.230                | 0                        | 287                         | 2     | 0.081     | 0.006       |
| krepp_estimated_10k | exact_stage43_partial | 293      | 0               | -                    | 0                        | 293                         | -     | -         | -           |

Pooling the two exact families by generating source count (1,188 → 892 exactly-explainable cases):

| K | cases | exact solver: support = truth | `dd_default` | `dd_true_k` | `dd_resolution` |
|---|---|---|---|---|---|
| 2 | 120 | 100 % | 99.2 % | 99.2 % | 99.2 % |
| 3 | 220 | 99.1 % | 96.4 % | 97.7 % | 96.8 % |
| 5 | 220 | 98.6 % | 80.9 % | 94.5 % | 90.5 % |
| 7 | 119 | 90.8 % | 63.0 % | 93.3 % | 90.8 % |
| 10 | 213 | 2.8 % | 18.3 % | 61.0 % | 55.4 % |

The exact-solver column is end-to-end: it counts a case only if the solver completed within 120 s
*and* its witness is the generating support, so the K = 7 and K = 10 entries are dominated by
timeouts, not by wrong answers. On completed cases the witness is the generating support in 196/200
(`exact_10k`) and 473/475 (`exact_biotrees`); of the 6 exceptions, 4 are alternative explanations of
the same minimal cardinality and 2 have a minimum strictly below the generating count.

The same picture at the level of individual edges rather than whole supports:

| K | `dd_default` | `dd_true_k` | `dd_resolution` |
|---|---|---|---|
| 2 | 99.6 % | 99.6 % | 99.6 % |
| 3 | 98.9 % | 99.2 % | 99.2 % |
| 5 | 95.8 % | 98.9 % | 98.8 % |
| 7 | 92.3 % | 98.7 % | 98.9 % |
| 10 | 82.7 % | 94.1 % | 94.3 % |

This is the central nuance: at K = 10 with the true k supplied, DecoDiPhy places **94 % of the
individual sources on the correct edge** but gets **the whole support right only 61 %** of the time.
Its errors are near-misses on one or two components, not wholesale failures. The Wasserstein-1
distance between the returned mixture and the generating mixture, under the tree metric, is
8.0 × 10⁻⁹ substitutions/site on the 783 `dd_true_k` cases that recover the support exactly and
3.6 × 10⁻³ on the 109 that do not — a median error three to four orders of magnitude below a typical
branch length even when the support is wrong.

### 4.3 Minimality

An explanation is counted only when the returned support is certified exactly feasible. The reference minimum is the exact solver's answer where it completed, otherwise the generating count, which is an upper bound on the true minimum.

| family         | arm           | cases | exact and minimal | exact but non-minimal | below reference minimum | no exact explanation | incomplete |
|----------------|---------------|-------|-------------------|-----------------------|-------------------------|----------------------|------------|
| exact_10k      | dd_default    | 296   | 193               | 4                     | 13                      | 86                   | 0          |
| exact_10k      | dd_true_k     | 296   | 254               | 2                     | 0                       | 40                   | 0          |
| exact_10k      | dd_resolution | 296   | 216               | 4                     | 24                      | 52                   | 0          |
| exact_biotrees | dd_default    | 596   | 439               | 2                     | 5                       | 150                  | 0          |
| exact_biotrees | dd_true_k     | 596   | 569               | 0                     | 0                       | 27                   | 0          |
| exact_biotrees | dd_resolution | 596   | 556               | 6                     | 8                       | 26                   | 0          |

Three distinct failure modes, all certified in exact arithmetic:

1. **No exact explanation returned** (86 / 40 / 52 cases in `exact_10k`; 150 / 27 / 26 in
   `exact_biotrees`). The returned support provably cannot reproduce the observations exactly, even
   though the data are exactly explainable and the reported loss is at machine precision — the median
   `dd_default` loss on `exact_10k` is 3.6 × 10⁻¹⁷. A floating-point residual of 10⁻¹⁷ is not evidence
   that the structure is right.
2. **Exact but non-minimal** (6 `dd_default`, 2 `dd_true_k`, 10 `dd_resolution` across both exact
   families). The returned support is a genuine exact explanation using more sources than necessary.
   Rare.
3. **Below the reference minimum** (33 distinct cases, all with K = 10, all cases where the exact
   solver timed out). Here DecoDiPhy returned a *certified exact* explanation with **fewer** sources
   than were used to generate the instance — so for these published benchmark instances the generating
   count is not the minimum. The exact solver completed on none of them, so this is a lower-bound
   statement, not a disagreement: no contradiction with a completed exact answer occurs anywhere in
   the run.

### 4.4 Identifiability

| family              | arm           | complete | feasible but non-identifiable | placement at a vertex | repeated edge |
|---------------------|---------------|----------|-------------------------------|-----------------------|---------------|
| exact_10k           | dd_default    | 296      | 42                            | 17                    | 0             |
| exact_10k           | dd_true_k     | 296      | 75                            | 19                    | 0             |
| exact_10k           | dd_resolution | 296      | 62                            | 18                    | 0             |
| exact_biotrees      | dd_default    | 596      | 46                            | 20                    | 0             |
| exact_biotrees      | dd_true_k     | 596      | 69                            | 25                    | 0             |
| exact_biotrees      | dd_resolution | 596      | 75                            | 27                    | 0             |
| krepp_estimated_10k | dd_default    | 290      | 0                             | 258                   | 0             |
| krepp_estimated_10k | dd_true_k     | 285      | 0                             | 223                   | 0             |
| krepp_estimated_10k | dd_resolution | 289      | 0                             | 256                   | 0             |

On 42–75 exactly-explainable cases per arm the returned support is exactly feasible but **not
identifiable**: the observations do not pin down a unique `(p, x, y)` on that support, so a continuum
of mixtures explains the data equally well. DecoDiPhy reports one point from that continuum with no
indication that it is not unique; the tool's own `check_identifiability` is only consulted on the
`--score` path, which is inoperable in this release. 17–27 cases per arm place a source exactly at a
vertex (`x` within 10⁻⁹ of 0 or 1), which is outside the PDD model's strict interior requirement. No
arm ever returned two placements on one edge.

### 4.5 Real krepp-estimated distances

Nothing changes qualitatively with the source count here, so the summary is flat:

| arm | completed cases | support = truth | individual true edges recovered | median W1 to truth | median CPU s |
|---|---|---|---|---|---|
| `dd_default` | 290 / 296 | **0** | 20–29 % | 0.079–0.087 | 2.95–5.50 |
| `dd_true_k` | 285 / 296 | **0** | 14–22 % | 0.073–0.079 | 0.75–5.61 |
| `dd_resolution` | 289 / 296 | **0** | 22–25 % | 0.079–0.083 | 1.86–4.08 |
| exact solver | 293 / 296 | — | — | — | 0.009 |

On all 293 cases where it completed, the exact solver returns *no exact explanation exists within the
budget*, in a median 9 ms: 237 rejected by a non-positive baseline pendant, 47 by the leaf-envelope
consistency test, 4 by the baseline mass-count obstruction, and 5 after exhausting the search. That is
the correct answer — these distances are ANI estimates and are not exactly explainable — and it is the
only thing the exact method has to say here. DecoDiPhy returns a least-squares mixture whose placements
are, on this data, mostly wrong at edge resolution: not one of the 864 completed DecoDiPhy runs
recovered the whole generating support of any instance.

Two further observations on this family. First, `dd_default` and `dd_resolution` choose a source count
different from the truth on 273 and 272 of 290/289 cases respectively (`dd_default`: 122 below, 151
above), stopping on low abundance in 276 of 290 cases — on noisy data the abundance stop, not the
residual stop, decides the answer. Second, 161 of the 296 krepp cases have 1–3 reference leaves for
which krepp returned no distance estimate; DecoDiPhy prunes those leaves internally and the exact
solver handles the partial observation natively through its observed backbone.

### 4.6 Robustness

DecoDiPhy exits with an unhandled exception on **11 of the 1,188 cases (0.9 %)**, all in
`krepp_estimated_10k` — the only family with unobserved reference leaves
(`results/stage289/stage289_failure_diagnostic.json`, 24 re-runs through the command-line entry point,
21 reproducing):

- **9 cases, `KeyError` in the jplace writer.** When the distance file does not cover every leaf,
  `preprocess_input` prunes the tree; the final placements are then mapped back onto the *unpruned*
  tree by `post_process_output`, while `save_jplace` builds its label dictionary from the *pruned*
  tree. Direct inspection of one case: at the moment `save_jplace` runs the tree has 248 of its
  original 252 nodes, and the returned anchor `N129` is one of the four that were removed.
  `all_rounds.json` has already been written at that point, so the search result survives; the
  primary `output.jplace` does not.
- **2 cases, `ValueError: 'a' cannot be empty unless no samples are taken`** inside the hill-climbing
  restart (`np.random.choice` on an empty candidate pool) in the `-f 1` arm. Here no final
  `all_rounds.json` is written and the run produces nothing.
- For one case (3 arms) the `KeyError` did not reproduce when the same inputs were re-run through the
  command-line entry point rather than in-process. I did not establish why.

---

## 5. Where the heuristic is adequate, and where it is not

**Stated plainly: for two or three sources the heuristic is adequate, and the exact method adds
nothing a practitioner would notice.** At K = 2 and K = 3 both tools recover the generating support on
96–100 % of cases, and the absolute runtime difference (≈ 0.02 s versus ≈ 1 s) is irrelevant at that
scale. If real samples are two- or three-source mixtures, the published tool is the right tool.

**Between four and seven sources the exact method strictly dominates, and this is the practical
argument.** At K = 5–7 it is simultaneously more accurate (98.6 % and 90.8 % whole-support recovery
end-to-end, against 63–94 % for the heuristic) and about 45× cheaper per case. The dominance is not a
tie-break on speed; the heuristic's errors here are real structural errors that an exact answer does
not make.

**At ten sources the exact method is simply unavailable** — 2.8 % completion at 120 s and 0 % at
1,750 s — and the heuristic is the only option. Given the true k it still places 94 % of individual
sources correctly, so it degrades gracefully rather than failing; but it returns the whole generating
support only 61 % of the time, and on 40 % of cases the support it returns provably cannot reproduce
the data exactly.

**The weakest part of the heuristic is its own choice of k, and that is exactly what the exact
formulation supplies.** With the oracle k, whole-support recovery at K = 5 is 94.5 %; with the tool's
default stopping rules it is 80.9 %, and at K = 10 the gap is 61.0 % against 18.3 %. The minimum-source
problem is a principled k-selection rule, and wherever it can be solved it selects better than
`min_p = 0.01` does.

**On real estimated distances the exact method is a screening certificate, not an estimator.** It
answers "these distances admit no exact explanation" in 9 ms and says nothing further; the heuristic
answers with a mixture whose edge-level accuracy on this data is 14–29 %. These are not competing
answers to the same question, and a paper should not present them as such.

**The honest summary for a reviewer.** The exact method's practical argument is real but bounded: it
is a better algorithm than the published heuristic for mixtures of up to roughly seven sources, both
in accuracy and in cost, on 892 exactly-explainable instances from the parent work's own benchmark and
from an independent phylogenetic benchmark. Beyond that it stops running, and on noisy real distances
it changes role from estimator to certifier. Any claim that the exact method replaces DecoDiPhy would
be false.

---

## 6. What this comparison does not establish

1. **It is not DecoDiPhy in its intended pipeline.** The tool is built around per-read krepp
   placements (`-j file.jplace`); here it received a distance vector (`-t/-d`), because that is what
   the PDD model and the frozen exact solver consume and because the bank stores distances. For the
   exact families the distance vector is a lossless encoding of the model, so nothing is lost; for
   `krepp_estimated_10k` the per-read placements may carry information the averaged ANI distances do
   not, and the parent work states that individual krepp placements consolidate more accurately. The
   krepp accuracy numbers here are therefore **not** a measurement of DecoDiPhy's accuracy in the
   pipeline it was designed for.
2. **Only one random seed.** Both repetitions used the documented default seed (1142); DecoDiPhy
   reseeds internally, so the repetitions measure timing variance, not search variance. Hill climbing
   starts from a random anchor set, so some of the structural errors reported here may be seed-specific.
   Sensitivity to the seed was not measured.
3. **Only the operable code paths.** `--method exhaustive/closest/closest_iterative/reverse_hill` and
   `--optimizer cvxpy` raise in v1.5.0, so the comparison covers hill climbing with OSQP only.
4. **It does not reproduce the parent paper's experiments.** Only the example shipped with the tool
   was reproduced. No figure or table of the published work was recomputed.
5. **No large reference sets.** The largest tree here has 361 leaves. The bank's 15,941-leaf
   `exact_16k` case and the 10,575-leaf travellers'-diarrhoea samples were not run in this track.
6. **Feasibility verdicts are about the exact rational data.** DecoDiPhy solves a double-precision
   version of the same system. A support called "provably infeasible" is infeasible for the exact
   instance; the effect of input rounding on that verdict was not quantified, though the observed
   structural errors (wrong edge, not wrong digit) are far larger than double rounding.
7. **Runtimes are indicative, not a benchmark.** One machine, single-threaded, 12 concurrent jobs,
   one wall budget. Ratios within a case are meaningful; absolute seconds are not portable.
8. **"Minimum" means minimum within the case budget.** The exact solver searches up to `cap = K_true`.
   Since the data are explainable with `K_true` sources, that minimum is the global minimum for the
   exact families; for `krepp_estimated_10k` the verdict is "no exact explanation with at most
   `K_true` sources", not "no exact explanation at all".

---

## 7. Files

All paths relative to `~/projects/recomb2027_novelty` on `ssh:blip`.

| file | contents |
|---|---|
| `data/stage289_decodiphy/repo` | the DecoDiPhy clone at commit `a26aa618…` |
| `data/stage289_decodiphy/venv` | the dedicated environment |
| `data/stage289_decodiphy/repro/jplace` | re-run of the shipped example |
| `results/stage289/stage289_headtohead_protocol.json` | registered protocol: environment record, limits, case hashes, code hashes |
| `results/stage289/stage289_headtohead_runs.jsonl` | 9,504 raw run records, timeouts and errors retained |
| `results/stage289/stage289_headtohead_percase.json` | derived per-case metrics |
| `results/stage289/stage289_headtohead_summary.json` | aggregated summary |
| `results/stage289/stage289_tables.md` | the full comparison tables, including the by-reference-size breakdown not reproduced above |
| `results/stage289/stage289_verification.json` | the audit (status PASS) |
| `results/stage289/stage289_repro_jplace.json` | the reproduction comparison |
| `results/stage289/stage289_exact_cost_probe.json` | the 1,800 s cost probe |
| `results/stage289/stage289_failure_diagnostic.json` | the DecoDiPhy failure diagnostic |
| `scripts/stage289_*.py` | 7 new scripts: bridge, feasibility, harness, analysis, report, probe, verification, diagnostic |

No pre-existing file was modified, renamed or deleted.

---

## 8. Deviations from the brief

- **Case coverage: none.** All 296 `exact_10k`, all 596 `exact_biotrees` and all 296
  `krepp_estimated_10k` cases were run, with two repetitions of all four methods (9,504 jobs). The
  bank's other families (`published_true_10k`, `published_biotrees`, `exact_16k`,
  `published_travelers_diarrhea`) are outside this track's brief and were not run.
- **The exact-cost probe is a subsample**: 18 cases (6 at K = 7, 12 at K = 10), one repetition, chosen
  by stratifying on reference-set size. Its conclusion ("K = 10 does not complete at 1,750 s") is
  provisional on those 18 cases; the 120 s result covers all 213 K = 10 cases.
- **Python 3.12.3 instead of the README's suggested 3.11** — no conda and no 3.11 interpreter on the
  host. All pinned dependency ranges in `pyproject.toml` are satisfied.
- **Per-job wall budget 120 s**, chosen from a pilot and justified post hoc by the 1,800 s probe;
  identical for both tools.
- **One seed** (the documented default), as noted in §6.2.
