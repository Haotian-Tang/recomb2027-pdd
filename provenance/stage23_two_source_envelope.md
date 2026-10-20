# Stage23: exact recovery of the support with at most two strict sources

2026-10-13. Candidate specialization/strengthening found while reading the predecessor's supplement. This is an internal AI derivation with exact finite controls, not external review or publication-priority clearance. The predecessor already proves equal-weight two-source edge identifiability in SB.3.6 Claim5, PDF pp21–23, lines1136–1155. The present statement allows arbitrary positive weights but keeps **every competing explanation restricted to at most two consolidated strict sources**.

## 1. Precise statement and classical ingredients

Let T be a finite positive-length tree with no degree-two vertices and with every original leaf observed. The two-leaf single-edge case is included. Let a strict explanation use at most two distinct open edges, with positive abundances summing to one, attachment points q_i, and free positive mean pendant y. Define its potential at every original vertex by

    f(v)=Σ_i p_i dist_T(v,q_i)+y,       d_r=f(r).

**Proposition.** The observations determine every original-vertex potential, its equivalent endpoint mass vector and its mean pendant uniquely among all such at-most-two-source explanations. Specifically,

    f(v)=h(v):=max_(r observed leaf) [d_r−dist_T(v,r)].

The distinct strict placement-edge support is also unique. This does not, in general, identify individual positions or abundances on adjacent edges, or any individual pendant lengths. It does not exclude explanations with three or more sources.

The envelope formula is classical: it is the lower McShane–Whitney 1-Lipschitz extension, not a new extension construction. Petrakis, *McShane–Whitney extensions in constructive analysis*, LMCS16(1:18),2020, introduction pp18:2–18:3, explicitly gives the supremum formula and its lower bound on every Lipschitz extension. [Primary full text](https://lmcs.episciences.org/6105/pdf). The tree inverse, a two-pass tree maximum-distance transform, and uniqueness of a forest perfect matching are also established tools. The candidate additional observation here is saturation of that lower envelope at **all original vertices** for the stated small source class, followed by its consequences in this PDD model. This combination still needs a dedicated broader prior comparison.

## 2. Empty-direction proof of envelope saturation

For any vertex v and any observed leaf r, the triangle inequality gives

    dist_T(r,q_i)−dist_T(r,v) ≤ dist_T(v,q_i).

Multiplying by nonnegative p_i and using Σp_i=1 yields d_r−dist_T(r,v)≤f(v), hence h(v)≤f(v).

If v is an original leaf, the term r=v attains equality. If v is internal, it has degree at least three. Removing v produces at least three components, whereas at most two components contain the two attachment points. Choose a component containing neither attachment. It contains an original leaf r, because T is finite and the branch cannot terminate at an unobserved original internal vertex. The path from that r to each attachment passes through v. Thus

    dist_T(r,q_i)=dist_T(r,v)+dist_T(v,q_i)

for both sources, so d_r−dist_T(r,v)=f(v). This proves equality at every original vertex. The proof works for unequal weights and does not require y to be known or to match a competitor's mean in advance. The assertion is about original vertices; the envelope can differ from the affine interpolation of these values inside a source edge.

Every competing explanation with at most two sources satisfies the same argument, so all of them have the same f. For general k, an internal vertex can have a source in every direction, and the equality argument can fail.

## 3. Endpoint inversion and support uniqueness

Let L be the reciprocal-length Laplacian and τ_v=2−deg(v). The existing tree inverse identifies

    m=(τ−Lf)/2,
    y=[τ^T f−Σ_e ell_e]/2.

These are therefore fixed by d within the at-most-two model. A strict source on uv contributes positive mass to both endpoints. The set P={v:m_v>0} is exactly the union of the endpoints of its used edges. Any competing strict used edge must have both endpoints in P, and the used edges must cover P.

- |P|=2: one edge is necessary and is unique; two distinct strict edges cannot have only two distinct endpoints in a tree.
- |P|=3: the two used edges are adjacent. The induced forest on P must be the two-edge path, so both edges are forced.
- |P|=4: two used edges must be disjoint and form a perfect matching of the induced forest. A forest has at most one perfect matching: two different ones would have a cycle in their symmetric difference.

This also handles competing supports that overlap in one edge, without relying on an exhaustive drawing of query quartets. A genuine two-distinct-edge explanation cannot have a one-edge strict competitor. Consequently the number one versus two is identifiable **within this restricted class**, after same-edge aggregation. It is not a claim that an unrestricted generating sample has that count.

For nonadjacent used edges, their disjoint endpoint masses also determine each edge's abundance and internal position. For adjacent edges u–v–w, write the private endpoint masses as a=m_u>0 and b=m_w>0, with m_v=1−a−b>0. Any first-edge abundance p∈(a,1−b) splits the middle mass positively and gives valid strict positions. There is a continuous family. If both abundances are specified as1/2, that splitting and the positions are fixed when feasible. In all cases, individual terminal lengths remain visible only through their weighted mean.

This recovers the predecessor's Claim5 support conclusion under its assumptions and explains why its adjacent-edge Claim3 is compatible with that conclusion. The strengthening to unknown weights is our inference, not the wording of Claim5.

## 4. Complete exact test for at most two sources

For arbitrary input d on all leaves of a positive reduced tree:

1. Compute h(v)=max_r[d_r−dist(v,r)].
2. Reject if h(r)≠d_r at any observed leaf. Otherwise an inconsistent observation could be silently raised by the envelope.
3. Apply the inverse to obtain m,y. Reject if any m_v<0 or y≤0.
4. Reject if more than four vertices have positive mass. Inspect edges with two positive endpoints and test whether one or two cover all positive vertices. Reject if none does.
5. Split each positive endpoint mass among its selected incident edges, obtaining strict positions and positive weights. Return the resulting explanation.

Necessity follows from §§2–3. Sufficiency follows from the inverse identity and positive endpoint splitting; Σm=1 follows from Στ=2 and the Laplacian row sums. The reconstructed potential agrees with every input leaf. The cover, if it exists, is unique. A rejection means no explanation with at most two strict sources; it does not decide unrestricted strict feasibility.

**Complexity.** The envelope has the standard two-pass recurrence on a rooted tree. The downward pass maximizes over leaves in the subtree after subtracting traversed edge lengths. The upward pass combines the parent's outside value, the local observation when present, and the best sibling value. Keeping the largest two child contributions permits all exclusions in linear work, including vertices of high degree. Inversion uses one traversal of the edges. When |P|≤4, the induced forest has at most three edges, so the residual cover test is constant-sized after a linear scan. Thus this gives O(|V|) rational arithmetic operations and O(|V|) rational storage. This is an arithmetic count, not an O(|V|) bit-time claim: sums, maxima, divisions by input edge lengths and a bounded number of source ratios have polynomial rational bit length in the encoded input. Their exact arithmetic therefore gives polynomial bit complexity.

The implemented solver does not build a full distance matrix. For a returned explanation, it additionally checks its exact potential by at most four endpoint-to-all-vertices tree scans. The finite test harness separately builds all-pairs distances as an independent small-instance reference, so the harness runtime is not a large-scale performance benchmark.

## 5. Exact controls and necessary scope

`scripts/stage23_two_source_envelope.py` imports only Python's standard library. Its140 positive controls span five trees: a single edge, an unequal three-star, a quartet, a five-leaf caterpillar, and a high-ratio four-star. Every one- and two-edge support is used; the pairs have three explicit weight/position settings including equal weights, unequal weights, and near-boundary positions. The pendant is10^-9. All controls verify original path distances, the direct all-leaf maximum, the two-pass maximum, inversion, unique strict edge cover, and the reconstructed source witness. They passed in approximately0.139s. This finite set accompanies the proof and is not itself a proof of universal identifiability.

Additional controls deliberately limit the conclusion:

**Adjacent supports, nonunique parameters.** On a unit three-star, sources on armsA/B with radial positions(1/2,1/2), weights(1/2,1/2), and pendant1 give (2,2,5/2). The same edges with radial positions(3/4,3/8) and weights(1/3,2/3) give exactly the same observations and pendant. The edge support is fixed while parameters differ.

**Three sources can have different supports.** On the unit quartet with edges UA,UB,UV,VC,VD, use UA/UB sources of weight1/6 at fraction3/4 from U, and a UV source of weight2/3 at fraction15/16 from U, with pendant11/8. This gives all four observations equal to3. Its mirror uses VC,VD,UV with the same weight multiset and pendant. All positions are strict. The two internal potential pairs are(9/4,2) and(2,9/4); their common leaf envelope gives(2,2), so saturation fails at one internal vertex in each world. This root-supplied control is independently verified in the script and is consistent with the predecessor's symmetric same-k counterexample theme.

**Missing leaf.** If only A and B are observed in the unit three-star, one midpoint source on the unobserved C arm with pendant1/2 and the equal A/B midpoint pair with pendant1 both give(2,2). Observing C distinguishes them. The all-leaf requirement is material.

**Degree-two vertices.** On a four-edge unit path with only its endpoints as references, two equally weighted midpoint sources on the first/fourth edges and on the second/third edges both give(3,3) with pendant1. Their strict edge supports differ. The reference is not reduced; the empty-direction argument is unavailable at its degree-two vertices.

**Strict versus closed.** Equal d=(2,2,2) on the unit three-star is rejected by the at-most-two strict solver. Its envelope places all mass at the center; no positive-endpoint cover exists. A closed center source fits, and a strict three-arm explanation exists. This prevents replacing the cover test by nonnegative inverse masses alone.

**Inconsistent d.** The input(100,1,1) on that star is rejected because the envelope changes an observed leaf value. This checks the necessary forward consistency gate.

Script SHA256: `527930cd95e0ab9c684bb04bea5b2e63172794c2ac5000754449e718a316ab1a`.

Result `results/stage23_two_source_envelope.json` SHA256: `8379bd17735f7a64e331d07747f8a4009634f941a3701ac2d57208c59edeb1e3`.

## 6. Consequence for the project, not a novelty conclusion

The new reading corrects an overly negative general impression of PDD identifiability: the predecessor already has a positive two-source theorem. The argument here offers a concise, conditional exact characterization of the at-most-two class and a standard-tree-DP implementation. It neither removes the general variable-K hardness nor settles the general strict n−1 conjecture. It is not a noisy distance estimator, a robust biological recovery procedure, or evidence of speed/accuracy superiority over DecoDiPhy. Before treating it as an independent contribution, its saturation/reconstruction combination needs comparison with earlier tree placement, tree-metric mixture and sparse-support recovery results; the generic Lipschitz formula is explicitly excluded from novelty.
