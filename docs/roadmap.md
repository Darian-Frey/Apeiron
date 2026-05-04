# Roadmap

Research plan for Apeiron. Expands CLAUDE.md §6 and §8 into a
phase-by-phase schedule. Treat phase boundaries as commitments to
reach a reviewable checkpoint, not as hard deadlines. For
day-to-day operational state see [STATUS.md](../STATUS.md).

## Phase 0 — Foundations

**Complete.** All six §5.1 modules shipped with full acceptance
criteria, oracle tests, and integration tests:

- [x] `zphi.py` — exact ℤ[φ] arithmetic.
- [x] `symmetry.py` — icosahedral group I (60 rotations) under ×2
      storage convention.
- [x] `polyhedron.py` — convex hull with coplanar-face merging,
      validation predicates, isometry action, RTH integration test.
- [x] `substitution.py` — `SubstitutionRule`, `substitution_matrix`,
      `is_primitive`, `perron_frobenius_in_zphi` (P3 + Fibonacci
      oracles).
- [x] `corona.py` — face-to-face placement, completeness predicates
      (`incidence_defect`, `has_interior_overlap`), `corona_1` BFS
      engine, `corona_2` second-shell extension. Cube and rhombic
      dodecahedron acceptance fixtures.
- [x] `hierarchy.py` — `Supertile` lazy-recursive form, pillar 2
      (`is_recognisable` with Ammann–Beenker oracle), pillar 3
      (`inflation_argument` with P3 oracle), pillar 4 framework
      (`FourthPillarArgument` protocol; concrete implementations
      live per-candidate).
- [x] `viz.py` and `notebooks/` — interactive HTML rendering for
      the RTH, the Danzer prototiles, and the cube's first corona.

## Phase 1 — Track A (deformation-first)

**In progress.** The Danzer ABCK 4-tile set is the baseline
candidate; a successful Track A would identify a minimal-deformation
merge of the four tiles into one that retains the substitution
property.

- [x] Encode Danzer's 4 prototiles with exact ℤ[φ]³ vertices
      (`candidates/danzer/{A,B,C,K}.json` from Frettlöh Table 1).
- [x] Pillar 1 verified: Frettlöh's substitution matrix recovered
      exactly, primitive, PF eigenvalue = φ³ = ZPhi(1, 2).
- [x] Pillar-4 stub (`candidates/danzer/fourth_pillar.py`) cites the
      published proof (Danzer 1989 Thm 2; Goodman-Strauss 1998;
      Frettlöh Thm 1.2) and exercises the protocol's import path
      end-to-end.
- [x] **27B-β** done in `b3d8e55` + `164e6dd` — Paolini POV-Ray
      source machine-extracted into
      `candidates/danzer/paolini_dissection.json` with all 25
      children verified by volume conservation, multiplicity match,
      and I_h membership; promoted to canonical `danzer_rule`
      fixture per Claude (web)'s Q4a 2026-04-29 ruling. (Frettlöh's
      live URLs were dead; Paolini's POV-Ray is the canonical
      machine-readable source.)
- [x] **27D** done in `5bfa3a5` + `875ba02` — full four-pillar
      pipeline returns `InflationArgument(pf=φ³, radius=2)` on the
      canonical real-geometry rule; pillar 2 succeeds uniformly
      across all four prototiles at level 1 and at level 2 for
      σ(A) using the Goodman-Strauss `position_signature`.
      Empirical recognisability radii: σ(A)→2, σ(B)→6, σ(C)→1,
      σ(K)→1, σ²(A)→4. Track A's baseline confirmed: pillars 1 + 2
      + 3 hold for ABCK with real geometry.
- [x] **Phase 1.5** done in commits `27f84f1`, `6732fbc`, `4406900`,
      per Claude (web)'s Q6 ruling 2026-04-29 (face-merge formalism,
      brute-force enumeration with cheap filters, scaffold-time
      stub).
      `apeiron/deformation.py` provides `face_adjacent_pairs`,
      `merge_two_tiles`, `prioritise_merge_candidate`, and
      `scaffold_merge_candidate`. The face-merge enumeration over
      σ(A,B,C,K) finds 24 face-adjacent pairs, ZERO of which pass
      Q6b's combinatorial-feasibility filter — no σ-count vector
      matches any 2-tile composition with the right structure.
- [x] **Phase 1.5 no-go (2-tile face-merges)** — face-merging two
      ABCK tetrahedra into a single tile cannot yield a 1-tile
      substitution rule with PF in ℤ[φ]; the σ-count vector
      mismatch is decisive. Concrete: σ(A)+σ(B) = (0, 5, 3, 10),
      no factorisation through (1, 1, 0, 0). Tested in
      `tests/test_deformation.py::TestPrioritiseMergeCandidateOnDanzer::test_passes_all_aggregate_is_zero_for_danzer`.
- [x] **Phase 1.6** done in `42b3d9a` — closed-form no-go on the
      entire face-merge track for ABCK at every composition size.
      Per Claude (web)'s Q7a/b 2026-04-29 ruling: a face-merge
      produces a valid 1-prototile substitution iff the substitution
      matrix has a positive-integer eigenvalue (with a non-negative-
      integer eigenvector). M_ABCK's eigenvalues are
      {φ³, φ, -1/φ, -1/φ³}; none is a positive integer. Therefore
      no face-merge of any number of children yields a 1-prototile
      substitution rule with PF in ℤ[φ].
      Implementation: ``feasibility_upper_bound(rule)`` does the
      closed-form check via the rational root theorem on the integer
      characteristic polynomial; ``enumerate_face_merge_compositions``
      short-circuits to empty for any rule that fails the check.
      Tested at k_max up to 8 on the canonical Paolini Danzer rule.
- [x] **Track A no-go theorem (face-merge formalism)** — recorded.
      The face-merge deformation search over Danzer's ABCK is closed
      with no survivors. Either (i) Track B (substitution-first), (ii)
      vertex perturbation around a hypothetical near-solution, or
      (iii) a different Track A formalism entirely (continuous
      parametric families per Q6a's secondary path).

## Phase 2 — Track B (substitution-first)

**In progress** as of `d994384` (2026-04-29). Track A's face-merge
formalism is closed by the Phase 1.6 no-go; Track B is the next
direction per Claude (web)'s Q7c ruling. Sub-package layout
(`apeiron/track_b/`):

- [x] `matrix_search.py` — `enumerate_primitive_matrices(n,
      pf_target, *, max_entry=10)` yields primitive non-negative
      integer substitution matrices with target PF eigenvalue
      computed exactly in ℤ[φ]. Deduplicated under simultaneous
      row/column permutation. n=2 algebraic survey at max_entry=5:
        - PF=φ:   1 candidate (Fibonacci [[0,1],[1,1]]).
        - PF=φ²:  1 candidate (Penrose P3 [[1,1],[1,2]]).
        - PF=φ³:  5 candidates — concrete 2-tile algebraic
          alternatives to Danzer's 4-tile rule. Each has trace=4,
          det=−1; geometric realisation pending.
- [x] `geometric_prefilter.py` ✅ in `a84593d`. Filter 1 (vertex-
      class, ZPhi PF eigenvector) implemented for n=2; filters 2
      (dihedral-angle) and 3 (Euler) are documented stubs returning
      None. All 5 PF=φ³ candidates pass filter 1.
- [x] `realisation.py` API ✅ in `f97d24c` per Claude (web) Q8
      ruling 2026-04-29. Result types
      (`Realised | NoRealisation | Inconclusive`),
      `ChildPlacement`, `SearchProgress`, public
      `realise(matrix, pf_target, *, max_search_seconds=300)`
      function. Fibonacci oracle returns Realised via a manual
      witness; non-Fibonacci inputs return Inconclusive with
      ``fraction_searched=0`` and a reason string referencing the
      pending CSP.
- [x] **Realisation CSP body** ✅ in commits `72d222f`, `589c4ca`,
      `d49c3d7`, `6949bec`, `c0ef1d9` per Claude (web)'s Q8a/c
      ruling. Three building blocks plus three validation layers
      now compose end-to-end:

      Building blocks:
        * `translation_offset_from_face_match` — per-edge ZPhi
          offset solver (6 vertex permutations).
        * `propagate_translations_along_tree` — tree-DFS chain of
          per-edge offsets.
        * `FaceMatchEdge` dataclass.

      Search loop (`_search_realisation_for_parent`):
        * Iterates DFS-discovery-ordered tree topologies.
        * Per edge: 4×4 face-pair sequence.
        * Per non-root child: rotation pool (60 proper).
        * Child 0's rotation fixed to identity, translation to
          origin (Q8 meta-1).

      Validation layers (in order, with early bail):
        1. Volume sum: pure ℤ[φ] check that
           ``Σ vol(child) = pf_target × vol(parent)``. Bails entire
           search if mismatched (saves all rotation iterations).
        2. AABB: every placed vertex inside the inflated parent's
           axis-aligned bounding box (cheap pre-filter).
        3. Polytope containment: every placed vertex inside the
           inflated parent polytope (centroid-relative, per face).
        4. SAT (Separating Axis Theorem): pairwise interior-disjoint
           check on all C(k,2) pairs. Tangent-at-boundary counts as
           disjoint per face-to-face semantics.

      With all four layers passing, "Realised" is mathematically
      equivalent to a face-to-face dissection (by inclusion-
      exclusion: vol-sum + non-overlap + containment ⇒ coverage,
      since any non-covered region would be open in the parent's
      interior with measure 0, which is empty).

      Remaining limitation:
        * **k > 3**: outer iteration ``(k-1)! × 16^(k-1) × 60^(k-1)``
          becomes intractable past k=3 without fail-first ordering
          per Q8 meta-3. The CSP body is a complete decision
          procedure for k ≤ 3 today.
- [x] **DFS-style backtracker with face-match pruning** ✅ in
      `551510a`. Replaces the triple-nested
      ``trees × face-pair-seqs × rotation-seqs`` iteration with a
      recursive DFS that explores ``(face_pair, rotation)`` triples
      per edge and prunes inconsistent face-match choices at the
      source via ``translation_offset_from_face_match``. All 5
      n=2 PF=φ³ candidates now answered definitively (in 0.00s
      with synthetic identical-volume tets — volume mismatch
      rejects every candidate by design).
- [ ] **Shape derivation from eigenvector** to close the loop on
      ``realise(matrix, pf_target)`` without requiring caller-supplied
      prototile_shapes. Per Q8e's "fix canonical shape from volume
      ratio first": derive icosahedral-compatible tetrahedra whose
      volume ratio matches the left eigenvector of M. Currently
      realise() rejects wrong shapes definitively; with shape
      derivation it could attempt the *right* shapes too.
- [ ] **Fail-first rotation ordering optimisation** for k > 3:
      sort rotation pool by per-rotation constraint score (number of
      face-match violations when placed first). Currently DFS
      explores rotations in fixed order; sorting could give 5–10×
      speedup for large k. Optional; the existing DFS handles
      k=2–3 in seconds.
- [ ] Integration test on the n=2 PF=φ³ survey: run
      `realise` on each of the 5 candidates with appropriately-shaped
      prototiles. Most or all are expected to return NoRealisation;
      any survivor is a major Track B finding.

## Phase 3 — Recognisability and beyond

The recognisability machinery (pillar 2) and the inflation
argument (pillar 3) are already implemented and oracle-tested.
What remains here is candidate-specific:

- [ ] Compute the recognisability radius for any surviving Phase 1
      / Phase 2 candidate.
- [ ] Confirm the inflation argument terminates with a witness
      rather than a failure.

## Phase 4 — The fourth pillar

The pillar-4 framework is complete; concrete implementations are
written per candidate, and only after pillars 1–3 are established
for that candidate.

- [ ] Per-candidate `fourth_pillar.py` implementations (case
      analysis ruling out non-hierarchical tilings).
- [ ] If the fourth pillar fails for every candidate that passed
      pillars 1–3: write up the no-go theorem.

## Non-goals

- A product. A library for third parties.
- Performance before correctness.
- A 2D reimplementation. Prior 2D work is reference-only; see
  CLAUDE.md §5.3.
