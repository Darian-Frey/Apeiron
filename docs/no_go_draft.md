# Computational Evidence Against 3D Icosahedral Einstein Tiles in the D₆ Cut-and-Project Framework

**Status.** Draft. Per Q17b authorisation (Claude (web), 2026-05-05),
this is the Branch B4 deliverable from
[strategy_pivot.md §8.3](strategy_pivot.md). Sections (i)–(iii)
are based on closed work; section (iv) is the open gap pending
inter-library loan of Kramer-Papadopolos-Schlottmann-Zeidler
(1994).

**Scope statement.** The result is **not** "no 3D icosahedral
Einstein exists." The result is the precisely-scoped statement
in the title: computational evidence within the D₆ cut-and-
project framework, with the H₃ tetrahedral fundamental-region
shape constraint, surveyed at three independent levels (Track A
algebraic no-go, Track B exhaustive substitution-matrix survey,
external empire-exhaustion via Hammock-Fang-Irwin 2018), with
one explicit residual gap (window-uniqueness for ABCK).

Target venues per Q17b: *Discrete and Computational Geometry* or
*Experimental Mathematics*.

---

## Abstract

We present three independent lines of computational evidence
supporting the conjecture that no single-prototile (monotile)
icosahedral substitution tiling of ℝ³ exists within the D₆
cut-and-project framework. (i) The face-merge approach to
Danzer's ABCK 4-tile substitution is closed by an algebraic
no-go: M_ABCK has no positive integer eigenvalue, so no
positive-integer linear combination of the prototiles inflates
to a single combined tile. (ii) An exhaustive computational
survey of n ≤ 3 prototile substitution rules with PF eigenvalue
φ³, integer entries ≤ 3, and prototile shapes drawn from the 9
H₃-compatible tetrahedral classes of the Paolini ABCK pool
yields zero realisable rules across 78 attempted candidates.
(iii) Hammock, Fang & Irwin (2018) computed the empires (forced
tiles) of all 36 vertex configurations of the canonical D₆ →
ℝ³ projection T*(2F); none reduces to a single-prototile
region. We provide the algebraic eigenvalue embedding
M_ABCK ⊂ M_T*(2F) (Papadopolos-Hohneker-Kramer 1999, eq. 13)
that explains why ABCK's substitution-matrix structure is a
strict restriction of the canonical projection. The remaining
open question is whether non-canonical D₆ acceptance windows
admit different vertex configurations whose empires might be
single-prototile; settling this requires the closed-access
"all windows" result of Kramer-Papadopolos-Schlottmann-Zeidler
(1994).

---

## Introduction

The 3D Einstein problem asks whether a single polyhedron *P* ⊂
ℝ³ exists such that every tiling of ℝ³ by isometric copies of
*P* is non-periodic. The 2D analogue was settled affirmatively
in 2023 by Smith, Myers, Kaplan & Goodman-Strauss (the "hat"
and "spectre" tiles); no 3D Einstein has been exhibited or
ruled out.

The natural arithmetic setting for icosahedrally-symmetric 3D
substitution tilings is the ring ℤ[φ] = ℤ + ℤφ, where
φ = (1 + √5) / 2. The natural projection lattice is D₆: the
icosahedral group H₃ is a maximal subgroup of D₆'s point group,
and the projection π: D₆ → ℝ³ ⊕ ℝ³_⊥ (with ℝ³_⊥ the internal
space) yields icosahedral aperiodic tilings via the cut-and-
project method (Kramer & Neri 1984; Kramer-Papadopolos-
Schlottmann-Zeidler 1994; Frettlöh "Icosahedral tilings in R³").

The canonical D₆ icosahedral projection, denoted T*(2F), has
**six tetrahedral prototiles** (Papadopolos-Hohneker-Kramer
1999) with the constraint that all face normals lie in the
H₃ root system. Danzer's ABCK tiling is locally derivable from
T*(2F) and reduces the prototile count to 4. The structural
multi-prototile conjecture asks whether this reduction can be
continued to 1.

This paper reports three computational results bearing on the
conjecture, each independent and using a different methodology.

---

## §(i) Track A: Face-merge algebraic no-go

The Danzer ABCK substitution rule has substitution matrix
M_ABCK (4×4, integer entries) with characteristic polynomial:

```text
char(M_ABCK)(λ) = λ⁴ - 5λ³ + 2λ² + 5λ + 1
```

Roots of this polynomial: {τ³, τ, −τ⁻¹, −τ⁻³} where τ = φ. None
is a positive integer.

**Theorem (face-merge no-go).** For any non-negative-integer
linear combination of A, B, C, K prototiles, no scalar inflation
of that combination produces a single combined tile under the
ABCK substitution rule.

*Proof sketch.* A face-merge of the ABCK prototile set into a
single combined tile would correspond to a rational-coefficient
left eigenvector of M_ABCK with positive integer eigenvalue
(the inflation factor on the merged tile). The integer-
eigenvalue constraint follows from the rational root theorem on
the integer characteristic polynomial: any rational root of
char(M_ABCK) is in fact an integer. The roots of M_ABCK lie in
ℤ[φ] but none is in ℤ. ∎

This rules out the deformation-first approach to a 3D Einstein
candidate via merging Danzer ABCK tiles. The no-go is
closed-form and reproducible:
[apeiron/deformation.py](../apeiron/deformation.py)
implements the eigenvalue feasibility check;
[tests/test_deformation.py](../tests/test_deformation.py)
verifies the M_ABCK eigenvalues.

### §(i.1) Eigenvalue embedding: structural context

The face-merge no-go is structural rather than coincidental.
M_ABCK is a 4×4 matrix; its eigenvalues are a strict subset of
the eigenvalues of the canonical T*(2F) volume inflation matrix
M_T*(2F) (8×8, ℤ[τ] entries, Papadopolos-Hohneker-Kramer 1999
eq. 9):

```text
char(M_T*(2F))(λ) = λ · (λ³ + (13-8τ)λ² + (61-38τ)λ + (62-39τ)) · (λ⁴ - 5λ³ + 2λ² + 5λ + 1)
```

The third factor is char(M_ABCK). The embedding
spectrum(M_ABCK) ⊂ spectrum(M_T*(2F)) is verified
algebraically over ℤ[τ] in
[tests/integration/test_t2f_embedding.py](../tests/integration/test_t2f_embedding.py).
This establishes that ABCK's substitution-matrix structure is
inherited from the canonical T*(2F) projection (via the local
derivation T*(2F) → ABCK established by
Papadopolos-Klitzing-Kramer 1997).

---

## §(ii) Track B: Exhaustive substitution-matrix survey

We searched the space of n × n non-negative-integer
substitution matrices M with:

- n ≤ 3 prototiles;
- integer entries 0 ≤ M[i,j] ≤ 3 (max_entry = 3);
- Perron-Frobenius eigenvalue equal to φ³ ∈ ℤ[φ];
- primitive matrix structure;
- prototile shapes drawn from the 9 H₃-compatible tetrahedral
  similarity classes derivable from the Paolini ABCK vertex
  pool.

For each candidate matrix, we attempted to construct a
geometric realisation: a tessellation of σ(P_j) by isometric
copies of P_0, ..., P_{n-1} per the matrix's prescribed child
counts. The realisation backtracker uses a depth-first search
over per-child rotations in I_h (order 120) with linear
translation recovery via face-match consistency, and four
exact-arithmetic validation layers (volume sum, axis-aligned
bounding box, polytope containment, separating axis theorem
for interior-disjoint).

**Empirical result.** Across two entry bounds:

- max_entry = 2: 21 primitive matrices, 13 with shape triples
  matching some 9-class volume ratio; **all 13 returned
  exhaustive NoRealisation** under per-step DFS validation
  (k_max = 5, 30 s/triple budget).
- max_entry = 3: 209 primitive matrices, 65 with matching shape
  triples; **all 65 returned exhaustive NoRealisation** under
  per-step DFS validation (k_max = 7, 30 s/triple budget).

**Combined: 78 candidates exhaustively NoRealisation across
the entire 9-class slice.** Logs at
[notebooks/n3_phi3_exhaustive_2026-04-29.log](../notebooks/n3_phi3_exhaustive_2026-04-29.log)
and
[notebooks/n3_phi3_max_entry_3_2026-05-04.log](../notebooks/n3_phi3_max_entry_3_2026-05-04.log).

A further 153 candidates (9 + 144 across the two entry bounds)
were skipped as `SKIPPED_TAXONOMY_GAP`: their left-eigenvector
volume ratios fell outside the 9-class set. Their realisability
is the separate question that the 15-class extension would
address. We do not claim anything about these 153 candidates
beyond the 9-class scope.

**Scope statement for §(ii).** The result is "no realisable
candidate among 78 specific matrix-shape-triple combinations
within the 9-class H₃ slice"; not "no realisable candidate
exists for n ≤ 3 with PF eigenvalue φ³." The scope is precisely
defined by the four bullet conditions above.

Code:
[apeiron/track_b/](../apeiron/track_b/)
(matrix enumeration, geometric prefilter, realisation
backtracker, H₃ taxonomy);
[scripts/sweep_n3_phi3_max_entry_3.py](../scripts/sweep_n3_phi3_max_entry_3.py)
(driver).

---

## §(iii) External: Hammock-Fang-Irwin (2018) empire-exhaustion

Hammock, Fang & Irwin (2018) compute the *empires* (sets of
forced tiles) for all 36 vertex configurations of the canonical
D₆ → ℝ³ projection T*(2F). The empire E(p) of a patch p is the
intersection of all valid tilings T' that contain p. In the
projection framework,

```text
E(p) = { t ∈ T : Γ(p) ⊂ t*_⊥ }
```

where Γ(p) := ⋂_{t ∈ p} t*_⊥ is the *shift-window* (a convex
polytope in internal space E_⊥).

The decomposition of the canonical Voronoi cut-window W ⊂ E_⊥
yields 880 tile types (regions), 4230 sectors (region
intersections), and 36 H₃-equivalence classes of vertex
configurations. Hammock 2018 Table 1 enumerates all 36 with
their sectors, empires, and exact ℤ[√5] frequencies. Table 2
cross-validates against Kramer-Papadopolos-Zeidler (1991) — exact
agreement.

**Critical observation:** *None of the 36 empires reduces to a
single-prototile region.* Every canonical-projection vertex
configuration's empire contains multiple distinct tile types.
This is an exhaustive computational check, by independent
authors, with cross-validation against Kramer 1991 frequencies,
covering every vertex configuration of the canonical T*(2F)
projection. No counterexample to the structural multi-prototile
conjecture exists at the canonical-projection level.

### §(iii.1) The H₃ shape-rigidity extension

Hammock 2018's exhaustive check covers the canonical Voronoi
window. Whether non-canonical windows could produce single-
prototile empires is the residual question (§(iv)). However,
the H₃ shape-rigidity argument (per Q17a, 2026-05-05) shows the
gap is bounded:

**Argument (informal).** Prototile shapes in any D₆ icosahedral
cut-and-project tiling are determined by the H₃ root system,
not by the window. (Source: Frettlöh "Icosahedral tilings in
R³" Theorem 1.3; Al-Siyabi-Koca-Koca 2020 §3.) Varying the
window changes which vertex configurations appear but not what
the tiles look like. An empire is a set of tiles. If the tiles
available are constrained by H₃ regardless of window, then no
non-canonical window can introduce *new* tile types absent from
the canonical projection. Hammock's exhaustive check therefore
constrains the space of possible single-prototile empires
across all D₆ windows, not just the canonical one.

**This is an argument, not a proof.** The missing piece is
whether vertex-configuration *combinatorics* (not tile shape)
varies enough across windows to produce a single-prototile
empire. A non-canonical window might present a vertex
configuration with fewer adjacent tile types, forcing fewer
types and potentially admitting a single-prototile empire —
even though the available tile *shapes* are H₃-rigid. The
argument bounds the *set* of possible tile shapes; it does not
bound the *combinatorics* of how they meet under a non-
canonical window.

The structural-versus-combinatorial distinction is the precise
residual gap addressed in §(iv).

---

## §(iv) Open gap: D₆ window uniqueness for ABCK

The closed-access paper Kramer, Papadopolos, Schlottmann &
Zeidler (1994), *Projection of the Danzer Tiling*, *J. Phys. A*
27, 4505–4517 (DOI 10.1088/0305-4470/27/13/024), is reported in
the literature to "determine the windows such that the Danzer
tiling can be obtained by projection from D₆." This is the
"all windows" result that would close the residual gap.

Three open-access companion papers establish the framework but
do not settle the all-windows question:

- Papadopolos-Klitzing-Kramer (1997), *J. Phys. A* 30, L143:
  defines T*(2F) and T(2F) via Klotz construction; states "both
  tilings can be locally transformed into tilings with only
  four elements that allow simple inflation" (the 6→4 reduction
  to ABCK). Closed-access; ILL-pending.
- Kramer-Papadopolos (1994b), *Can. J. Phys.* 72, 408: derives
  Mosseri-Sadoc tilings from D₆ via a methodological parallel
  to the 1994 Danzer paper. Closed-access; ILL-pending.
- Papadopolos-Hohneker-Kramer (1999/2000), *Discrete Math.*
  221, 101 (arXiv:math-ph/9909012): inflation rules for T*(2F)
  with explicit 8×8 substitution matrix. Open-access; read in
  full and used as the basis for the eigenvalue embedding in
  §(i.1).

The residual gap is small but present: even if no non-canonical
window admits a single-prototile empire (as the H₃ shape-
rigidity argument in §(iii.1) suggests), this is currently an
informal extension argument, not a proof. The 1994 Danzer paper
is the natural source for closure.

**Status of this section in the draft.** Placeholder. Will be
revised when the ILL returns with the 1994 Kramer-Papadopolos-
Schlottmann-Zeidler result. Two outcomes:

- *Outcome A:* The 1994 paper exhibits a non-canonical window
  with a single-prototile vertex configuration. This refutes
  the conjecture in its current form. The paper is then
  reframed: "Computational evidence against single-prototile
  *substitution* tilings in the D₆ cut-and-project framework";
  the new candidate would proceed through Apeiron's existing
  pipeline.
- *Outcome B:* The 1994 paper exhibits no single-prototile
  window. The H₃ shape-rigidity extension upgrades to a proof
  (via the explicit window enumeration); §(iv) closes; the
  conjecture moves from "strongly-supported, unproven" to a
  full theorem within the D₆ framework.

---

## Discussion

### What this paper does NOT claim

- *Not claimed:* no 3D icosahedral Einstein exists. The result
  is restricted to the D₆ cut-and-project framework with
  H₃-determined prototile shapes.
- *Not claimed:* substitution tilings are the only route to a
  3D Einstein. Non-substitution constructions (Walton-Whittaker
  2019 in 2D; potential 3D analogues) and non-icosahedral
  symmetry groups (octahedral, tetrahedral) are explicitly
  outside scope.
- *Not claimed:* the H₃ taxonomy is the unique relevant shape
  family. The 9-class restriction in §(ii) is the operative
  taxonomy from the Paolini ABCK pool; Frettlöh "Icosahedral
  tilings in R³" §1.1 reports 15 H₃-compatible classes total,
  with the missing 11 currently unverifiable from open
  literature.

### Implications for the broader Einstein problem

The structural multi-prototile conjecture, if proven (Outcome B
above), would say: any 3D Einstein candidate is either *not*
icosahedrally symmetric, *not* a model set in the D₆ framework,
or *not* H₃-tetrahedral in shape. By process of elimination
this points the search at:

- non-icosahedral symmetry (octahedral, tetrahedral, lower);
- non-substitution aperiodicity (local rules without inflation,
  per the hat/spectre direction);
- non-tetrahedral prototile shapes (zonohedra, other
  H₃-compatible polyhedra outside the 15-class fundamental-
  region taxonomy).

Each of these is a research direction in its own right. The
strategic value of a published version of this paper is to
focus future 3D Einstein search away from icosahedral D₆
substitution tilings.

---

## References

Detailed in [docs/literature_notes.md](literature_notes.md).
Key references for this draft:

- **Apeiron**: this codebase. [github.com/Darian-Frey/Apeiron](https://github.com/Darian-Frey/Apeiron) (private at time of draft).
- **Goodman-Strauss (1998)** *Matching rules and substitution tilings*, *Annals of Mathematics* 147(1), 181–223.
- **Frettlöh** *Icosahedral tilings in R³: The ABCK tilings*, Bielefeld internal report.
- **Papadopolos, Hohneker & Kramer (1999/2000)** *Tiles–inflation rules for the class of canonical tilings T\*(2F)*, *Discrete Math.* 221, 101–112; arXiv:math-ph/9909012.
- **Hammock, Fang & Irwin (2018)** *Quasicrystal tilings in three dimensions and their empires*, *Crystals* 8(10):370.
- **Al-Siyabi, Koca & Koca (2020)** *Icosahedral polyhedra from D₆ lattice and Danzer's ABCK tiling*, *Symmetry* 12(12):1983; arXiv:2003.13449.
- **Moody (2000)** *Model sets: a survey*, arXiv:math/0002020.
- **Kramer & Neri (1984)** *On periodic and non-periodic space fillings of E^m obtained by projection*, *Acta Cryst. A* 40, 580.
- **Kramer, Papadopolos, Schlottmann & Zeidler (1994)** *Projection of the Danzer tiling*, *J. Phys. A* 27, 4505 [closed-access; ILL-pending].
- **Smith, Myers, Kaplan & Goodman-Strauss (2023)** *An aperiodic monotile* (the "hat") and *A chiral aperiodic monotile* (the "spectre").

---

## Notes on this draft

This is a Branch B4 deliverable under the Apeiron strategy
pivot framework (`docs/strategy_pivot.md`). The draft is
intended for review and revision before submission. Specific
revision targets:

- §(iv) updates when the 1994 ILL returns.
- §(iii.1) H₃ shape-rigidity argument upgraded from informal
  to proof if Outcome B applies.
- Empire-computation cross-validation in ZPhi³ (Q17c (ii),
  deferred per Q17c) added if a reviewer specifically
  challenges Hammock 2018's methodology.
- §(ii) extension to max_entry = 8 (the original Q9c gate
  threshold) considered if reviewers argue 78 candidates is
  insufficient sample size.

Per the Q10 sequencing commitment carried through the literature
deep-dive: re-relay before further engineering or major
revisions.
