# Apeiron

Research framework for the **3D Einstein problem**: the search for a single
polyhedron *P* ⊂ ℝ³ such that every tiling of ℝ³ by isometric copies of *P*
is non-periodic.

Exact arithmetic over ℤ[φ], icosahedral symmetry, substitution tilings,
four-pillar proof apparatus. See [CLAUDE.md](./CLAUDE.md) for the full
problem statement, algebraic framework, repository conventions, and
current priorities; [STATUS.md](./STATUS.md) for the rolling
operational state (current focus, commit log, open relays);
[docs/decisions.md](./docs/decisions.md) for the chronological log of
non-obvious design decisions.

## Status

The six-module pipeline (CLAUDE.md §5.1) is complete:
[`zphi`](apeiron/zphi.py) → [`symmetry`](apeiron/symmetry.py) →
[`polyhedron`](apeiron/polyhedron.py) →
[`substitution`](apeiron/substitution.py) →
[`corona`](apeiron/corona.py) → [`hierarchy`](apeiron/hierarchy.py).
The four-pillar proof structure (substitution exists / recognisability /
aperiodicity-from-recognisability / no-non-hierarchical-tiling) is
encoded; pillar 4 is candidate-specific by design (concrete
implementations live at `candidates/<name>/fourth_pillar.py`).

**Track A** (CLAUDE.md §6.1, deformation-first) is currently in
progress on Danzer's published 4-tile aperiodic baseline. Pillar 1
verification holds against Frettlöh's substitution matrix
(PF eigenvalue = φ³). Pending: real geometric dissection (sub-commit
27B-β) followed by full four-pillar pipeline run (27D).

## Layout

- [`apeiron/`](apeiron/) — core package: ℤ[φ] exact arithmetic,
  icosahedral group, polyhedra with hull + coplanar-merge,
  substitution rules with primitivity / Perron–Frobenius,
  corona BFS engine, hierarchy / four-pillar apparatus,
  shared utilities (`util.py`), and a `viz.py` for plotly-based
  3D rendering.
- [`tests/`](tests/) — pytest suite, one file per core module
  plus [`tests/integration/`](tests/integration/) for cross-module
  acceptance tests (rule: a test is integration ⇔ it imports from
  ≥ 2 core modules).
- [`candidates/`](candidates/) — candidate tile definitions: JSON
  vertex lists per the loader's author-facing schema, plus
  per-candidate `fourth_pillar.py` implementations of the
  pillar-4 protocol.
- [`docs/`](docs/) — algebra notes, research roadmap, annotated
  literature, decision log.
- [`notebooks/`](notebooks/) — runnable scripts that produce
  interactive HTML visualisations: the rhombic triacontahedron
  ([`rth_demo.py`](notebooks/rth_demo.py)), the Danzer ABCK tiles
  ([`danzer_demo.py`](notebooks/danzer_demo.py)), the cube's
  26-neighbour first corona
  ([`cube_corona_demo.py`](notebooks/cube_corona_demo.py)). HTML
  outputs are gitignored; regenerate with
  `.venv/bin/python notebooks/<demo>.py`.

## Development

```sh
uv venv
uv pip install -e ".[dev]"            # core + test + lint
uv pip install -e ".[dev,viz]"        # add plotly / pyvista for notebooks/
.venv/bin/python -m pytest             # full suite, ~23 s
ruff check .
mypy
```

The build backend is `hatchling`; environment / dependency management
is `uv` (per the 2026-04-19 design decision in
[docs/decisions.md](docs/decisions.md) and STATUS.md). System Python
without a venv works for everything except the corona acceptance
fixtures, which depend on `scipy.spatial.ConvexHull` (declared in the
core dependencies).

## Discipline

This is a correctness-first codebase. The verification pipeline is
exact in ℤ[φ] from end to end — no float arithmetic anywhere except
visualisation (in `viz.py`) and the convex-hull *combinatorial*
oracle (`scipy.ConvexHull` returns triangle index triples; the float
coordinates it computed are immediately discarded; the resulting
combinatorics get validated exactly via three ZPhi predicates in
[`polyhedron.py`](apeiron/polyhedron.py)). See
[docs/decisions.md](docs/decisions.md) and the CLAUDE.md §7
testing/style sections.
