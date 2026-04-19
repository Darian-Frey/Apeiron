# Integration tests

## Rule

A test lives in `tests/integration/` **iff it imports from two or more core
modules** under `apeiron/`. Tests that exercise only one core module live
in the top-level `tests/` directory as `test_<module>.py`.

This rule makes placement mechanical: count the `from apeiron.X import …`
lines. Two or more distinct `X`'s ⇒ integration.

## Expected tenants

Listed in approximate order of arrival. Each will be added when the core
modules it spans are stable.

- `test_rhombic_triacontahedron.py` — constructs the RTH fixture via
  `polyhedron.py` and asserts that the icosahedral group I (from
  `symmetry.py`) permutes its 30 faces transitively. First integration
  test; sets the convention.
- `test_substitution_realisation.py` — applies a substitution rule from
  `substitution.py` to a concrete polyhedron from `polyhedron.py` and
  checks that the inflated tile dissects into the specified multiset of
  isometric copies.
- `test_corona_on_periodic_tiles.py` — the cube and rhombic dodecahedron
  corona oracles. Touches `polyhedron.py`, `symmetry.py`, and `corona.py`.
- `test_hierarchy_pillars.py` — the four-pillar chain: substitution →
  recognisability → aperiodicity. Spans `substitution.py`, `corona.py`,
  and `hierarchy.py`.

Unit tests that happen to import shared fixtures defined under
`apeiron.util` are **not** integration tests — `util` is plumbing, not a
core verification module.
