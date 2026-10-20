# Stage 287 — Frozen real-data PDD instance bank

**Phase 1, track "Instances".** Built on `blip` under
`~/projects/recomb2027_novelty`. No solver was run in this track. All outputs
are new files with stage numbers 286–287; nothing pre-existing was modified.

---

## 1. What was established about the parent work's public benchmark (Step 1)

Script `scripts/stage286_benchmark_parse.py` → `results/stage286_benchmark_parse.json`
(read-only; 2.4 MB).

### 1.1 Column schemas — taken from the authors' own reader code, not guessed

| File | Columns | Evidence |
|---|---|---|
| per-instance `queries.txt` | `QUERY_LABEL`, `TRUE_ABUNDANCE` | `find_sol_10k.py`: `true_p = [float(q[1]) for q in queries]`, `k = len(queries)` |
| per-instance `distances.txt` | `REF_ID`, `TRUE_DISTANCE` | two columns in every file; `decodiphy/cli.py` reads `{l[0]: float(l[1])}` |
| per-instance `estimated_distances_pr-h14_dth4.txt` | `REF_ID`, `ESTIMATED_DISTANCE`, `PORTION_MAPPED`, `TRUE_DISTANCE` | literal header row; `find_sol_10k.py` uses `d[1]`, `d[2]`, `d[3]` for est / portion / true |
| per-instance `read_counts_table-h14_dth.txt` | `#FeatureID`, `ESTIMATED` | literal header; this is a krepp read-mapping **estimate** of abundance, *not* ground truth |
| `dist.txt` | `tree`, `k`, `s`, `dist` + the four estimated-distance columns | four key columns prepended by the `awk -v var1..var4` idiom in `commands.sh`; the four payload columns are spelled out verbatim in `krepp_distances.txt` |
| `krepp_distances.txt` | same 8, with one header row per instance group | header repeated inside the file |
| `krepp_placement_distances.txt` | `tree`, `k`, `s`, `dist`, `REF_ID`, `PLACEMENT_DISTANCE` | 6 columns; rows per group equal the **full** tree leaf count, not full − K |
| `truemash.txt` | `tree`, `k`, `s`, `dist`, `REF_ID`, `TRUE_MASH_DISTANCE` | equals per-instance `true_mash_distances.txt` with keys prepended |
| `all_metrics_est_mash_corrected.txt` | `tree`, `k`, `s`, `dist`, `jaccard`, `unifrac`, `wunifrac`, `est_k`, `wasser`, `bc` | named verbatim by the authors' `draw_10k.R` |

### 1.2 Instance inventory

* **656 distinct query instances** in `10k_dataset`, keyed by
  (reference tree, K, replicate, abundance distribution), over **22 reference
  trees**, with **K ∈ {3, 5, 10}**, 5 replicates and 2 abundance distributions
  (`uniform`, `exp`).
* **Reference-set sizes 20 – 129** (39 distinct values), counted as table rows
  and cross-checked against proper Newick leaf counts.
* Ground truth: all 656 publish `queries.txt` (query taxon + true abundance),
  i.e. the true source count **and** the true abundances.
* Estimator availability: 656 carry a true distance table (`distances.txt`);
  656 carry krepp ANI estimates (`estimated_distances_pr-h14_dth4.txt`,
  with `PORTION_MAPPED`); 208 carry `true_mash_distances.txt`; 60 carry
  krepp placement distances; 140 carry the published accuracy metrics
  (including the authors' estimated source count `est_k`).
* Aggregate files: `dist.txt` 35,622 lines / 656 groups / 22 trees;
  `truemash.txt` 10,116 / 208 / 7; `krepp_distances.txt` 4,050 / 60 / 2 (trees
  15 and 17 only); `krepp_placement_distances.txt` 4,350 / 60 / 2;
  `all_metrics_est_mash_corrected.txt` 140 rows / 5 trees.
  **263 of the 656 `dist.txt` groups contain rows with only 7 fields**, i.e. the
  `TRUE_DISTANCE` column is absent for part of the table.

### 1.3 The generative model, confirmed from the data

`find_sol.py` / `find_sol_10k.py` simulate a mixture by selecting K **leaves of
the full reference phylogeny** as queries, drawing abundances p, and reporting
`d_r = Σ_i p_i · dist_Tfull(q_i, r)` for the remaining leaves. Pruning the K
query leaves turns each q_i into a point source at the old position of its
neighbour, with pendant y_i equal to the query's own pendant edge, so that
`d_r = y + Σ_i p_i dist_T(r, q_i)` with `y = Σ_i p_i y_i` — the project's PDD
model. The prediction *full-tree leaf count = reference-set size + K* holds for
every on-disk tree (tree_1322: 30, tree_17: 34, tree_7: 50, tree_4: 89,
tree_15: 111).

### 1.4 A model-conformance problem in the parent data

`tree_1322.tre` (5), `tree_4.tre` (7) and `tree_7.tre` (6) carry **negative
branch lengths**, as do **all** `-ani_blen` trees. The PDD model requires
positive edge lengths, so these instances cannot enter the exact family
(see §3.3). The earlier public inputs in `data/stage08_distance_public`
(5 instances with `queries.txt`, `ani.tre`, krepp estimates and true-mash
distances) and `data/stage15_public_pairs` (2 `queries.txt` files) were parsed
and are subsumed by the full retrieval below.

---

## 2. What was fetched (Step 2)

Script `scripts/stage287_fetch_instances.py` →
`data/stage287_instances/download_manifest.json`.

**7,897 files, 91.8 MB, 0 errors.** Every file's **git blob SHA1 was recomputed
from the retrieved bytes and matched the repository's git tree manifest at the
named commit** (7,843 / 20 / 34 matches, 0 mismatches), which proves the content
is exactly that commit's. URL, commit, retrieval time, byte count and SHA256 are
recorded per file. No host was blocked.

| Repository | Commit | What was taken |
|---|---|---|
| `shayesteh99/Demixing_Data` | `5d3a1141158994af3fd6e0712fe6aaa14b52903f` | all 10k reference trees; `queries.txt`, `distances.txt`, `estimated_distances_pr-h14_dth4.txt`, `true_mash_distances.txt`, `distances_from_placement.txt` for all 656 instances; all 19 `biotrees` species trees and 4,172 instance `distances.txt`; WoL v2 backbone + 630 `16k` query compositions; the travellers' diarrhoea 16S backbone and real mixture distances |
| `bo1929/shared.krepp` | `0310db641d5edc18fc0c77804a3902ea1a2f7458` | 20 large real reference backbones (WoL v1/v2, RefSeq CAMI II) — an **independent** repository (krepp's author) |
| `JaspervB-tud/reference_set_selection_benchmark` | `ea5d512b3716bee4331f5689b02416640201bd67` | complete `.txt` / `.md` / `.py` inventory (34 files) |

### 2.1 The independent reference-set-selection benchmark contributes no PDD instances

The repository holds **no phylogenetic tree and no distance data of any kind**:
its 58 blobs are 18 `.txt` (NCBI accession lists for the family / genus / order
experiments), 13 `.py`, 7 `.sh`, 6 `.ipynb`, 9 `.pdf` figures, 3 `.md` and a
licence. Its reference sets are accession lists and its raw reads live on
Zenodo. It is conceptually on point for reference incompleteness — it measures
how reference-set *selection* changes profiling accuracy — but it is a
read-mapping benchmark, not a tree-based one, so it cannot supply a reference
tree, an observed distance vector or a ground-truth mixture. **It is recorded
for provenance and contributes zero cases.** The second independent data source
in the bank is therefore `bo1929/shared.krepp` (reference panels) together with
the 19 independently published species trees redistributed in `biotrees`.

### 2.2 Maximum real reference-set size obtained

Script `scripts/stage287_tree_census.py` → `results/stage287_tree_census.json`;
102 trees parsed, leaves counted from the Newick.

| Real tree | Leaves | Strictly positive edges? |
|---|---|---|
| `shared.krepp` RefSeq CAMI II, 380k with duplicates | **379,490** | no — 130,479 negative, 202,255 zero |
| `shared.krepp` RefSeq CAMI II dup, placement-extended | 371,557 | no — 178,910 zero |
| `shared.krepp` WoL v2 placement-extended (EMP / HMI) | **47,151** | **yes** |
| `shared.krepp` WoL v1 placement-extended | 31,177 | yes |
| WoL v2 backbone (`16k/backbone-WoLv2.nwk`, and 4 `shared.krepp` variants) | **15,953** | yes |
| travellers' diarrhoea 16S backbone | 10,575 | yes |
| `biotrees` prokaryotes WoL (CASTLES-Pro) | 10,575 | no — 132 zero |
| `biotrees` bees-Piskulich | 4,591 | yes (no mixtures published) |
| `biotrees` mammals-Upham | 4,099 | no — 1 zero |

* **Largest real reference tree of any kind: 379,490 leaves.**
* **Largest real tree that satisfies the PDD positivity requirement: 47,151
  leaves**; the largest clean reference *backbone* is WoL v2 at **15,953**.
* **Largest real tree that comes with published mixtures *and* ground truth:
  the WoL v2 backbone, 15,953 leaves** (the `16k` family, reference set 15,941).
* **Largest real tree with published *per-instance* ground truth at small K:
  `birds-stiller-castlepro.tre`, 363 leaves** (reference set 361).

This is 1.5–3 orders of magnitude beyond the 6–12 leaves used in the submitted
paper. Two `shared.krepp` trees could not be parsed because their edge lengths
carry jplace-style `{n}` annotations, and `tree_8-ani_blen.tre.rerooted` is
empty in the repository itself; all three are duplicates of trees that did
parse.

---

## 3. The frozen bank (Step 3)

Built by `scripts/stage287_build_bank.py` (library
`scripts/stage287_pdd_cases.py`, parser `scripts/stage286_newick.py`) and
checked by `scripts/stage287_verify_bank.py`.

* `results/stage287_instance_bank_protocol.json` — schema, provenance,
  conversion rules, per-case SHA256, counts, quarantine, panels (3.2 MB)
* `data/stage287_instances/bank/cases/<family>/<case>.json` — **5,686 case
  files** (bank 349 MB including the retrieved inputs)
* `data/stage287_instances/bank/index.json` — case index
* `results/stage287_instance_bank_verification.json` — read-back verification

### 3.1 Case schema

Exactly the format of `results/stage44_confirmation_protocol.json`:
`edges` = list of `[u, v, "length"]` with integer vertex ids (leaves are
`0..N-1`) and exact rational length strings; `observed` = map from leaf id
(string) to an exact rational distance string; `cap` = source budget;
`hidden_leaves` = list of leaf ids; `pendant` = the generating common mean
pendant, supplied so a prescribed-pendant task needs no recomputation;
`generating_sources` = `[{edge: index into edges, p: rational, x: rational
measured from the first endpoint of that edge}]`; `generating_pendant` = y.
Added keys: `leaf_labels` (leaf id → taxon label), `source_details`,
`portion_mapped` (krepp family), `provenance` (repo, commit, instance, tree
file + SHA256, K, seed, distance source, exact-explainability flag).
`hidden_leaves` is empty in the exact families: reference-incompleteness
variants are derived downstream, not frozen here.

All lengths and weights go through `Fraction(str)`; no float round trip.
Published decimal weights are converted to exact rationals and divided by their
exact sum, because the published decimals sum to 1 only up to float rounding.

### 3.2 Families and counts

| Family | Cases | Reference-set size (min–max) | Generating K | Observed distances | Exactly explainable |
|---|---|---|---|---|---|
| `exact_10k` | 296 | 20 – 129 | 3 / 5 / 10 (100/100/96) | recomputed exactly from the tree + published mixture | **yes** |
| `published_true_10k` | 296 | 20 – 129 | 3 / 5 / 10 | the authors' published true distances | no (float-rounded) |
| `krepp_estimated_10k` | 296 | 20 – 128 | 3 / 5 / 10 | krepp ANI estimates from reads | no |
| `exact_biotrees` | 596 | 20 – 361 | 2/3/5/7/10 (120/120/120/119/117) | recomputed exactly | **yes** |
| `published_biotrees` | 4,172 | 20 – 361 | 2/3/5/7/10 | published `noise0` and `noise1/noise2` × `l∈{100,1000,10000}` | no |
| `exact_16k` | 1 | 15,941 | 12 | recomputed exactly | **yes** |
| `published_travelers_diarrhea` | 29 | 10,575 | unknown (real sample) | DEPP-placement 16S distances | no |
| **total** | **5,686** | **20 – 15,941** | | | 893 exact |

By reference-set size across all families: ≤50 → 1,656; 51–150 → 2,400;
151–500 → 1,600; >2,000 → 30 (total 5,686).

### 3.3 Ground truth carried by each instance

* **10k (296 exact):** published query taxa **and** published true abundances
  (`queries.txt`); source positions x and the common mean pendant y are derived
  exactly from the real tree.
* **biotrees (596 exact):** the *support* is published only implicitly — it was
  recovered as `leaves(full tree) \ rows(distances.txt)` — and the weights are
  not published at all. They were **reproduced from the generator seed**:
  `find_sol.py` calls `random.seed(s); np.random.seed(s)` and, for noise type 0,
  draws `true_p = np.random.rand(K)` (normalised, redrawn while
  `min(p) < 0.01/K`) as the first NumPy draw. The reproduced weights were
  matched to the recovered support by least squares against the published
  distances: **596 of 596 accepted**, maximum residual 6.7 × 10⁻¹⁵, maximum
  deviation between fitted and seed-reproduced weights 1.8 × 10⁻¹⁴, with the
  minimum weight gap always ≫ the deviation so the assignment is unambiguous.
* **`published_biotrees`:** support and x positions are real and recorded;
  `p` is `null`, because the weights are not recoverable for the noisy variants
  from published data alone.
* **16k (1 exact):** real HMP-derived composition over the WoL v2 backbone
  (`query.txt`), abundances summing to exactly 1.
* **travellers' diarrhoea (29):** **no ground truth** — these are 29 real stool
  samples (metadata group `TD` / control is recorded per case); `cap` is `null`
  and must be declared downstream.

### 3.4 Validation outcome

Every exact-family case is accepted only if two independent code paths agree
exactly: the full-tree computation `Σ_i p_i dist_Tfull(q_i, r)` and the
pruned-tree edge parameterisation
`y + Σ_i p_i min(d_T(r,u_i) + x_i L_i, d_T(r,v_i) + (1−x_i) L_i)`.
`scripts/stage287_verify_bank.py` then **re-read all 5,686 case files from
disk**, re-checked their SHA256 against the protocol, and recomputed the second
path from the stored edge list and stored mixture alone:

* **0 SHA256 mismatches**, **0 failing cases** in any family.
* **893 / 893 exact cases reproduce every observed leaf distance exactly**
  (maximum |model − observed| = 0 in `exact_10k`, `exact_biotrees`, `exact_16k`).
* `published_true_10k`: maximum |model − observed| = 1.75 × 10⁻¹⁶,
  median 7.9 × 10⁻¹⁷ — pure float-repr rounding, as expected.
* `krepp_estimated_10k`: maximum 0.292, median 0.134 — genuinely estimated and
  not exactly explainable, which is the point of the later noise work.
* Structural invariants checked on every case: positive edge lengths, x strictly
  interior, distinct source edges, Σp = 1, y > 0, observed keys are leaf ids,
  hidden and observed leaves disjoint.

### 3.5 Quarantined and not built — 1,088 quarantined, 616 skipped

| Reason | Count | Meaning |
|---|---|---|
| `exact_10k:non_positive_edge_length` | 326 | instances on the 12 of 22 `10k` trees that carry negative or zero branch lengths; **no `-ani_blen` tree is positivity-conformant** |
| `exact_10k:source_at_vertex_not_inside_edge` | 20 | after pruning, a source sits at a vertex, not strictly inside an edge |
| `exact_10k:source_not_strictly_interior` | 8 | x = 0 or 1 |
| `exact_10k:sources_share_an_edge` | 6 | two sources land in the same original edge |
| `exact_16k:source_vertex_removed` | 13 | the real HMP compositions are not generated with the separation constraint, so pruning removes whole clades; **13 of the 14 small-K instances degenerate, leaving 1** |
| `biotrees:missing_input_file` | 715 | instance directories that publish `all_rounds.json` but no `distances.txt` |
| skipped by the frozen K ≤ 20 materialisation cap (`16k`) | 616 | see §4 |

Instance families deliberately **not** built, and why:

1. **`true_mash` distances (208 instances).** These are distances computed *on
   the ANI branch-length tree*, and every `-ani_blen` tree carries negative
   branch lengths, so there is no positivity-conformant tree to pair them with.
   The tables are fetched and frozen; no case was created.
2. **`biotrees` `plants-onekp` (1,178 leaves), `turtle-thomas` (593),
   `frogs-feng` (309).** These three datasets publish `pruned_tree.trees` but no
   `distances.txt` for any of their 150 instances, and `pruned_tree.trees` has
   all labels replaced by integers, so neither the observed vector nor the
   support can be recovered. 900 of the 5,072 `biotrees` instance directories
   lack `distances.txt`; 715 of them were reached locally and quarantined, the
   other 185 publish nothing we fetched. This is why the largest
   ground-truth-bearing `biotrees` tree in the bank is 363 leaves rather than
   1,178.
3. **`krepp_placement_distances` (60 instances) and `read_counts` tables.** Not
   converted; the placement tables are indexed over the full tree including the
   query taxa and the read-count tables are abundance estimates, not distances.
4. **Large reference panels.** 27 real trees of ≥ 1,000 leaves (17 of them
   positivity-conformant, up to 47,151 leaves) are frozen as **panels** with
   leaf counts, pathology counts and SHA256 — trees only, no mixture. Generating
   mixtures on them is a Phase 2 protocol decision with its own seeds, so no
   synthetic mixture was invented here.

---

## 4. Deviations from the brief

1. **`16k` materialisation cap.** Only 14 of the 630 `16k` instances have
   K ≤ 20; the frozen rule `SIXTEEN_K_MAX = 20` materialises those and **skips
   616**, because a case on the 15,953-leaf backbone costs ~2 MB and exact
   rational distance computation scales with K × N (K reaches 1,095). The
   skipped instances are listed with their K and are reconstructible from the
   frozen inputs. Combined with the 13 degenerate supports, the `16k` family
   contains **one** case.
2. **`biotrees` weights are reconstructed, not published.** See §3.3. The
   reconstruction is verified to 10⁻¹⁴ on all 596 instances, but it is a
   reconstruction from the generator seed, not a published ground truth.
3. **`krepp_estimated_10k` is paired with the positive-length tree
   `tree_X.tre`, not the `-ani_blen` tree the authors used**, because every ANI
   tree violates positivity. The substitution is recorded in each case's
   provenance note.
4. **`published_true_10k` / `published_biotrees` / `krepp_estimated_10k` /
   `published_travelers_diarrhea` are not exactly explainable by construction**
   and are labelled as such; only the three `exact_*` families carry
   `generating_pendant` and pass exact verification.
5. **References without a krepp estimate are moved to `hidden_leaves`.** Some
   `estimated_distances_pr-h14_dth4.txt` rows have an empty
   `ESTIMATED_DISTANCE`; those leaves are dropped from `observed` and listed as
   hidden, which is why `krepp_estimated_10k` reaches only 128 observed leaves
   where the exact family reaches 129.
6. **The independent benchmark contributes no cases** (§2.1).
7. **No solver was run**, as required. The least-squares step used to match
   reconstructed weights to the recovered support is instance reconstruction,
   not a source-budget solve.

---

## 5. Files written (all new)

| Path | Contents |
|---|---|
| `scripts/stage286_newick.py` | exact-arithmetic Newick reader + pathology diagnostics |
| `scripts/stage286_benchmark_parse.py` | Step 1 inventory and schema determination |
| `scripts/stage287_fetch_instances.py` | Step 2 retrieval with per-file provenance and git-blob verification |
| `scripts/stage287_tree_census.py` | leaf counts and edge pathologies for every fetched tree |
| `scripts/stage287_pdd_cases.py` | full-tree → PDD case conversion in exact arithmetic |
| `scripts/stage287_build_bank.py` | Step 3 bank builder |
| `scripts/stage287_verify_bank.py` | read-back verification of the frozen bank |
| `results/stage286_benchmark_parse.json` | benchmark inventory, schemas, leaf-count cross-check |
| `results/stage287_tree_census.json` | census of 102 trees |
| `results/stage287_instance_bank_protocol.json` | **the frozen protocol**: schema, provenance, per-case hashes, counts, quarantine, panels |
| `results/stage287_instance_bank_verification.json` | verification outcome |
| `data/stage287_instances/` | 7,897 retrieved files + `download_manifest.json` + `bank/` (5,686 cases, `index.json`) |
| `logs_stage287_fetch.txt`, `logs_stage287_census.txt`, `logs_stage287_bank.txt`, `logs_stage287_verify.txt` | run logs |

## 6. What this determines for Phase 2

* **Scaling in N with published ground truth** is testable on a clean grid:
  N ∈ [20, 361] with K ∈ {2, 3, 5, 7, 10} from 892 exact cases on 21 real trees,
  plus one exact case at N = 15,941.
* **Scaling to large N** must use the frozen reference panels (up to 47,151
  leaves positivity-conformant, 15,953 for the WoL v2 backbone) with
  Phase-2-generated seeded mixtures — no real published small-K mixture exists
  above N = 361.
* **Reference incompleteness** has a real, non-synthetic instance of itself
  already in the bank: the krepp estimated tables are missing estimates for some
  references, and those leaves are already recorded as `hidden_leaves`.
* **Noise work** has three graded, labelled levels on the same structures:
  float-rounding (10⁻¹⁶), simulated sequence noise (`noise1`/`noise2` ×
  l ∈ {100, 1000, 10000}) and real ANI estimation error (median 0.134).
* **Biological data** is present: 29 real stool-sample mixtures against a real
  10,575-leaf 16S backbone, with no ground truth — so claims there must be
  about certificates and screening, not about recovery accuracy.
