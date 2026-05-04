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
- [ ] `geometric_prefilter.py` — three filters before realisation
      attempt per Q7c: vertex-class consistency, dihedral-angle
      commensurability with icosahedral angles, Euler-characteristic
      consistency. Eliminates algebraically-realisable but
      geometrically-infeasible matrices before the expensive CSP.
- [ ] `realisation.py` — vertex-placement CSP given a filtered
      matrix candidate. Find vertex positions such that the
      substitution children fit face-to-face in the tile's local
      frame. Uses `polyhedron.py` and `corona.py` primitives.
- [ ] Integration test: 1D Fibonacci known to realise; random
      non-realisable matrix must fail filter 1 or 2.

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
