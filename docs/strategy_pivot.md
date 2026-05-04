# Strategy pivot — post-(c.1) closure, pre-(d.2) Baake-Grimm read

**Status.** Decision-forcing document, not exhaustive. Written per
Q15c (Claude (web), 2026-05-04). Its job is to make the next
strategic decision legible, not to answer it. Commit before reading
Baake-Grimm Vol. 1 Ch. 7 so the reading has a target.

**Date.** 2026-05-04.

---

## §1 Current state of evidence

**Track A — Danzer face-merge formalism (closed).** Q7a closed via a
closed-form no-go: M_ABCK has eigenvalues {φ³, φ, −1/φ, −1/φ³}, none
of which are positive integers. By the integer-eigenvalue feasibility
upper bound (rational-root-theorem on the integer characteristic
polynomial), no face-merge of the ABCK tile set produces a single-
prototile substitution. Phase 1.5 + 1.6 infrastructure landed; no
candidate found.

**Track B — n≤3 PF=φ³ icosahedral substitution (closed).**
Exhaustive computational survey across two entry bounds:

- max_entry=2: 21 primitive matrices, 13 with shape triples in the
  9-class H₃ taxonomy → all 13 returned exhaustive NoRealisation.
- max_entry=3: 209 primitive matrices, 65 with shape triples in the
  same taxonomy → all 65 returned exhaustive NoRealisation.

Combined: **78 candidates exhaustively NoRealisation across the
9-class H₃ slice; 0 Realised.** Plus 153 SKIPPED_TAXONOMY_GAP
candidates whose left-eigenvector ratios fall outside the 9-class
volume set — separate question requiring the 15-class extension.

The headline result is published-ready in scope:

> No icosahedral substitution monotile with PF eigenvalue φ³ exists
> among n ≤ 3 prototile rules with integer entries ≤ 3 and prototile
> shapes drawn from the 9 H₃-compatible tetrahedral classes
> derivable from the Paolini ABCK vertex pool. This is the first
> systematic computational survey of this family.

(See `docs/literature_review.md` headline result + sweep logs at
`notebooks/n3_phi3_*_2026-*.log`.)

**Cut-and-project literature gap (open).** Q14b authorised reading
Kramer-Neri 1984 + Kramer-Papadopolos-Schlottmann-Zeidler 1994.
Both firmly closed-access (agent-verified Unpaywall +
OpenAlex queries). Three open-access substitutes read instead:
Al-Siyabi-Koca-Koca 2020, remaining Frettlöh "Icosahedral tilings
in R³" pp.8-9, Tilings Encyclopedia (no dedicated ABCK page).

Substitutes confirm:

- ABCK is a D₆ projection (Frettlöh Theorem 1.3); parameterised by
  integer pair (m₁, m₂) with m₁+m₂ even (Al-Siyabi eq. 11).
- All published F-type icosahedral tilings (ABCK,
  Socolar-Steinhardt, T^(2F)) are multi-prototile.
- Al-Siyabi explicitly does NOT use a window-CPS — they use direct
  projection of the (m₁, m₂) lattice subset.

Substitutes leave open:

- **Q14b(c)** window-deformation literature for monotile candidate
  generation;
- **Q14b(d)** window as continuous parameter vs discrete choice;
- **Q14b(e)** window → prototile-shape map.

The (c)/(d)/(e) gap is exactly the infrastructure-design-relevant
slice. Without those answers, building a "window-search" CSP
proceeds on incomplete foundations.

---

## §2 The structural conjecture and its consequences

**Structural multi-prototile conjecture (working hypothesis,
unverified).**

> Every D₆ CPS tiling of ℝ³ with icosahedral symmetry requires more
> than one prototile type. If true, icosahedral CPS cannot contain
> a 3D Einstein, and Apeiron's search must either (a) look outside
> icosahedral symmetry, (b) look outside the CPS framework, or (c)
> abandon the substitution-tiling approach entirely.

**Falsifier.** A single D₆ CPS tiling of ℝ³ with icosahedral
symmetry whose prototile set has cardinality 1. Constructive
counterexample (an explicit window producing a single-prototile
tiling) refutes immediately. A theorem proving that the framework
admits a monotile (without exhibiting one) also refutes.

**Status.** Inference based on absence — the published F-type and
P-type tiling families are uniformly multi-prototile, AND the
substitute reads suggest (without proving) that windows pick
vertex classes rather than shapes (Frettlöh's three windows
produce the same {A, B, C, K} shape family). But this absence may
reflect "no one has tried" or "the framework structurally
precludes it"; the substitutes do not distinguish.

**Three branches of consequence:**

- **Branch A — conjecture false.** A D₆ CPS monotile exists. The
  search continues within icosahedral CPS, with broadened scope
  (alternative windows, larger search space). Apeiron's existing
  Track A and Track B infrastructure remains relevant; a new
  "window-search" CSP becomes the next infrastructure piece.
- **Branch B — conjecture true.** Icosahedral CPS is the wrong
  framework. Apeiron pivots away from D₆ projection. Existing
  infrastructure for ZPhi arithmetic, icosahedral symmetry,
  corona BFS, recognisability remains but is no longer load-
  bearing for the central search.
- **Branch C — conjecture undecidable in reasonable time.** No
  proof or counterexample within the project's timescale.
  Pragmatic pivot regardless: treat the conjecture as if true and
  proceed to Branch B options.

The Q15b ruling sequence — (d.3) write this document → (d.2)
Baake-Grimm Ch. 7 read → (d.1) acquire originals if needed — is
designed to discriminate between A, B, C.

---

## §3 Branch A search strategy

If Baake-Grimm reveals a counterexample or refutes the conjecture
otherwise, the icosahedral CPS framework remains viable and the
next infrastructure phase is window-search. Concrete progression:

**§3.1 Close the 15-class taxonomy gap.** The 153 SKIPPED_TAXONOMY_GAP
candidates from the n≤3 sweeps have eigenvector ratios that may or
may not be realisable; the missing 11 H₃ classes are needed to
attempt them. Path: acquire ABCK-Book (Danzer-Sonneborn-van
Ophuysen 1993) coordinates via institutional library or interlibrary
loan, cross-validate against Frettlöh's notation (warned to contain
errors), encode the 6 missing classes that complete the 9-pool to
15, re-run the sweeps. **Cost: 1-2 weeks of acquisition + encoding +
cross-validation; ~8 hours of compute for the re-sweep.**

**§3.2 Window-search CSP infrastructure.** If Branch A is correct,
the right search is over D₆ acceptance domains rather than
substitution matrices. Required new infrastructure:

- D₆ lattice arithmetic (additive overhead on top of existing ℤ[φ]³
  + icosahedral group infrastructure).
- Acceptance-domain (window) representation in the internal space
  ℝ³_⊥.
- Projection-tile machinery: given a lattice + window, generate the
  tile set; verify single-prototile property.
- Search loop: enumerate or optimise window shapes; test each for
  monotile property + aperiodicity.

**Cost: ~weeks of infrastructure engineering, dependent on whether
windows are continuous (Q14b(d) unanswered) or discrete (also
unanswered).**

**§3.3 Beyond icosahedral CPS within Branch A.** If neither §3.1 nor
§3.2 produces a candidate, broaden to other root lattices — B₆, E₆,
or higher-dimensional projection schemes (D₇, E₇). The required
expansion: ZPhi arithmetic generalises to other quadratic-irrational
rings as needed. **Cost: substantial; effectively a new project.**

---

## §4 Branch B / C pivot options

If Baake-Grimm confirms the structural conjecture (or leaves it
unresolved per Branch C), the icosahedral CPS route is closed and
Apeiron's central search must pivot. Four concrete alternatives,
each with cost/benefit:

**§4.1 (B1) Non-icosahedral symmetry groups.**
Octahedral (order 48), tetrahedral (order 24), or other crystallo-
graphic 3D point groups. Lower symmetry but potentially richer
monotile structure: lower-order groups have fewer constraints on
tile shapes. The hat / spectre 2D monotiles use lower-order
symmetry (cyclic, dihedral) than Penrose's 5-fold. By analogy, a
3D monotile may live in a lower-symmetry group than icosahedral.

- **Pro:** Direct continuation of substitution-tiling formalism.
  Existing pillar-1/2/3 machinery applies. Smaller search space
  per group, but more groups to scan.
- **Con:** Loses the algebraic elegance of icosahedral / ℤ[φ]
  structure. ZPhi arithmetic may not be the right algebra; ℤ[√2]
  for octahedral, etc.

**§4.2 (B2) Non-substitution aperiodicity (local rules without
inflation).** Socolar-Taylor 2011, Walton-Whittaker 2019 direction.
Aperiodicity enforced by local matching rules on tile shape +
decoration, without an underlying substitution.

- **Pro:** The hat / spectre proofs (2023) use a non-substitution
  case-analysis approach. Computational verification scales (the
  hat's proof was computer-assisted). 3D analogues are unexplored.
- **Con:** Pillar-2 (recognisability) and pillar-3 (inflation
  argument) become moot; pillar-4 (no non-hierarchical) becomes
  the central question and is harder. The four-pillar architecture
  in CLAUDE.md §6.3 would need restructuring.

**§4.3 (B3) Decorated tiles.** Walton-Whittaker direction
explicitly. Tile shape + edge / face / vertex decorations enforce
aperiodicity via matching rules.

- **Pro:** Walton-Whittaker shows this works in 2D; 3D analogue is
  natural. Existing matching-rule infrastructure (Goodman-Strauss
  framework) directly applies.
- **Con:** Changes the definition of "tile" (CLAUDE.md §2.2).
  Pillar 4 must show every tiling respecting the decorations is
  hierarchical — strictly stronger than every shape-tiling being
  hierarchical. Apeiron's stated target (CLAUDE.md §2.1) is
  shape-only Strong Einstein; this pivot is a research-direction
  change, not just an infrastructure change.

**§4.4 (B4) Reframe as no-go theorem programme.** If the structural
conjecture is true, *prove it*. Apeiron's contribution becomes the
no-go theorem rather than the candidate.

- **Pro:** Publishable negative result. The 78-candidate
  exhaustive NoRealisation is already partial evidence;
  formalising the structural argument would be the natural
  research contribution. Aligns with Goodman-Strauss's
  precedent (matching-rule sufficiency theorems are negative-
  flavoured by establishing what's necessary).
- **Con:** Different mathematics — algebraic / topological,
  possibly Dehn-invariant or face-pairing arguments — than the
  current verification-and-search architecture. The
  infrastructure built so far is not directly reusable. Apeiron's
  identity changes from "find a 3D Einstein" to "prove no
  icosahedral CPS Einstein exists."

---

## §5 Decision gate

After Baake-Grimm Vol. 1 Ch. 7 is read, the decision is one of:

**Gate A — Branch A activates** if Baake-Grimm contains:

- A counterexample (a D₆ CPS monotile, possibly mentioned in
  passing or as a citation).
- A theorem stating that D₆ CPS admits monotile windows for some
  parameter range.
- An explicit treatment of window-deformation as a candidate-
  generation strategy.

In any of these cases, §3 progression starts. Re-relay
authorising window-search infrastructure (and possibly the
15-class extension first) becomes the natural next step.

**Gate B — Branch B activates** if Baake-Grimm contains:

- A theorem stating that all D₆ CPS icosahedral tilings have ≥ 2
  prototile types (proves the structural conjecture).
- A general result on minimum-prototile-count for CPS tilings
  that implies multi-prototile-ness for icosahedral D₆.
- An explicit no-go remark on icosahedral monotile candidates.

Re-relay authorising one of (B1)-(B4) under §4 becomes the next
step. The choice among B1-B4 should consider the existing
infrastructure (which favours B2/B3 over B1/B4, since pillar-1/2/3
machinery survives) and the user's strategic preference (B4 has
the highest publication-impact-per-time ratio if true).

**Gate C — Branch C activates** if Baake-Grimm contains:

- Discussion of the conjecture-relevant questions but no decisive
  resolution.
- A reference to Kramer-Papadopolos 1994 as containing the
  resolution (in which case (d.1) acquire-the-originals becomes
  active).
- Silence on the structural question entirely (in which case
  (d.1) is the next step regardless).

Re-relay authorising (d.1) acquisition or a pragmatic Branch B
pivot becomes the next step.

---

## §6 What this document is and is not

**This document is** a decision-forcing artefact. Its purpose is
to enumerate the strategic options and the criteria for selecting
among them, before the next round of reading produces evidence
that should resolve the choice.

**This document is not** a research roadmap, an infrastructure
plan, or a commitment to any specific pivot. The recommended next
step is reading Baake-Grimm Ch. 7 with the structural conjecture
active; the infrastructure consequences of the resulting decision
are deferred to the gate-resolution relay.

**Sequencing commitment carried forward.** No infrastructure work
lands before the gate decision. The pattern Apeiron has been
pushing back against — "build infrastructure, find no survivors,
build more infrastructure" — would re-emerge if any of §3 or §4
were started without the gate first.
