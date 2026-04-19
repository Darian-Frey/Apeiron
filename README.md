# Apeiron

Research framework for the **3D Einstein problem**: the search for a single
polyhedron *P* ⊂ ℝ³ such that every tiling of ℝ³ by isometric copies of *P*
is non-periodic.

Exact arithmetic over ℤ[φ], icosahedral symmetry, substitution tilings.

See [CLAUDE.md](./CLAUDE.md) for the full problem statement, algebraic
framework, repository layout, conventions, and current priorities.

## Layout

- `apeiron/` — core package (ℤ[φ] arithmetic, symmetry, polyhedra,
  substitution, corona BFS, hierarchy).
- `tests/` — pytest suite, one file per module.
- `candidates/` — candidate tile definitions.
- `docs/` — algebra notes, research roadmap, annotated literature.
- `notebooks/` — exploration and visualisation.

## Development

```sh
pip install -e ".[dev]"
pytest
ruff check .
mypy
```
