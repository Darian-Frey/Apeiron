# Integration tests

## Rule

A test lives in `tests/integration/` **iff it imports from two or more core
modules** under `apeiron/`. Tests that exercise only one core module live
in the top-level `tests/` directory as `test_<module>.py`.

This rule makes placement mechanical: count the `from apeiron.X import …`
lines. Two or more distinct `X`'s ⇒ integration.

Unit tests that happen to import shared fixtures defined under
`apeiron.util` are **not** integration tests — `util` is plumbing, not a
core verification module.

## Current tenants

- `test_rhombic_triacontahedron.py` — RTH fixture loaded from
  `candidates/rhombic_triacontahedron.json` via `apeiron.util`,
  validated against `apeiron.polyhedron`'s convex-hull machinery,
  checked for icosahedral symmetry under the 60 elements of `I` from
  `apeiron.symmetry`. First integration test; established the
  convention.
- `test_rhombic_dodecahedron.py` — RD fixture loaded similarly; the
  acceptance oracle for `apeiron.corona.corona_1` with per-vertex
  expected counts (3-valent vs 4-valent vertex types meet 4 vs 6
  tiles respectively). Spans `polyhedron`, `symmetry`, `corona`,
  `util`.
- `test_danzer_abck.py` — Track A's first candidate: Frettlöh's
  Danzer ABCK 4-tile aperiodic set. Loads all four prototiles,
  verifies tetrahedron shape, builds the `SubstitutionRule`
  (Frettlöh's matrix with placeholder dissection geometry per
  sub-commit 27B-α; real geometry from Frettlöh Figure 2
  transcription deferred to 27B-β), runs pillars 1 and 3 of the
  four-pillar verifier on it via
  `apeiron.hierarchy.inflation_argument`, exercises the
  `FourthPillarArgument` protocol via the
  `candidates/danzer/fourth_pillar.py` stub. Spans `polyhedron`,
  `substitution`, `hierarchy`, `corona`, `util`.

## Planned

- `test_corona_2_acceptance.py` — currently lives inline in
  `tests/test_corona.py::TestCoronaTwoCube`; technically spans
  `corona`, `polyhedron`, and `symmetry`. Could be promoted to
  integration if the convention is read strictly. Defer until
  `corona` gets multi-prototile support, at which point the test
  surface grows beyond what fits cleanly in a single-module test
  file.
- `test_substitution_realisation.py` — applies a substitution rule
  from `substitution` to a concrete polyhedron from `polyhedron`
  and checks that the inflated tile dissects into the specified
  multiset. Waits on Track A's 27B-β (real Danzer geometric
  dissection) — that commit will produce the first realisation
  fixture.
- `test_four_pillars_end_to_end.py` — pillar 1 → pillar 2 → pillar
  3 → pillar 4 chain on a real candidate via real (not synthetic)
  recognisability data. Track A's sub-commit 27D is the natural
  home for this once 27B-β has landed.
