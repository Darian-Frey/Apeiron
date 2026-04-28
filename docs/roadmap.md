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
- [ ] **27B-β** — transcribe real geometric dissection (translation
      + rotation per child of σ(X) for X ∈ {A,B,C,K}) from
      Frettlöh Figure 2 / Tilings Encyclopedia interactive view
      into `candidates/danzer/dissection_notes.md` and update the
      `SubstitutionRule`.
- [ ] **27D** — full four-pillar pipeline run on the Danzer
      baseline using the real geometric dissection.
- [ ] Systematic search for minimal-deformation merges of the four
      tiles into one. Interval arithmetic over ℤ[φ] for the
      parameterised tile family. Each candidate gets its own
      `candidates/danzer_merge_*/` directory with prototile JSON,
      substitution rule, and `fourth_pillar.py` slot.
- [ ] No-go theorem write-up if every parametric merge fails the
      pipeline.

## Phase 2 — Track B (substitution-first)

**Not started.** Pursue if Track A reaches a no-go conclusion;
the no-go result then informs which algebraic constraints Track B
should impose.

- [ ] Algebraic search over ℤ[φ]-linear σ with eigenvalue φ² on
      small alphabets (n = 1, 2, 3).
- [ ] Primitivity filter on substitution matrices.
- [ ] Geometric realisation search for surviving σ candidates.

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
