# Roadmap

Research plan for Apeiron. This document expands CLAUDE.md §6 and §8 into a
phase-by-phase schedule. Treat phase boundaries as commitments to reach a
reviewable checkpoint, not as hard deadlines.

## Phase 0 — Foundations (current)

- [ ] Fix coplanar-face merging in `polyhedron.py` (CLAUDE.md §5.1).
- [ ] Rhombic triacontahedron and rhombic dodecahedron test fixtures.
- [ ] Scaffold `substitution.py`: `SubstitutionRule`, `substitution_matrix`,
      `is_primitive`.
- [ ] Scaffold `corona.py`: `CoronaConfig`, `corona_1`, `corona_2`.
- [ ] First-corona completion on the unit cube passes as a sanity check.

## Phase 1 — Track A (deformation-first)

- [ ] Encode Danzer's 4-tile set with exact ℤ[φ]³ vertices.
- [ ] Verify its substitution rule in Apeiron.
- [ ] Systematic search for minimal-deformation merges of the 4 tiles into
      1, with interval arithmetic over ℤ[φ].

## Phase 2 — Track B (substitution-first)

- [ ] Algebraic search over ℤ[φ]-linear σ with eigenvalue φ² on small
      alphabets (n = 1, 2, 3).
- [ ] Primitivity filter on substitution matrices.
- [ ] Geometric realisation search for surviving σ candidates.

## Phase 3 — Recognisability

- [ ] Border-forcing detector on a candidate supertile.
- [ ] Formal statement of the recognisability theorem for a surviving
      candidate, with the local neighbourhood radius computed.

## Phase 4 — The fourth pillar

- [ ] Case analysis ruling out non-hierarchical tilings.
- [ ] If this pillar fails for every candidate: write up the no-go theorem.

## Non-goals

- A product. A library for third parties. Performance before correctness.
- A 2D reimplementation. Prior 2D work is reference-only; see CLAUDE.md §5.3.
