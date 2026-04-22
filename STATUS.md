# STATUS

Living status document for the Apeiron build-out. Updated whenever a
commit lands, a decision is reached, a relay to Claude (web) is sent or
returned, or an attempt fails in a way worth remembering.

- **Current date:** 2026-04-20.
- **Scope of this file:** *what is happening right now*, *what has been
  tried and abandoned*, *what questions are open*. Durable guidance
  (invariants, conventions, research strategy) lives in
  [CLAUDE.md](CLAUDE.md); this file is operational.

---

## Current focus

Module 5 of the CLAUDE.md §5.1 sequence: `corona.py`. The algebraic
substitution layer is done (`SubstitutionRule`, `substitution_matrix`,
`is_primitive`, `perron_frobenius_in_zphi`) and oracle-tested against
Penrose P3 (PF = φ²) and Fibonacci (PF = φ). Next: the `CoronaConfig`
dataclass, `corona_1`, `corona_2`, and canonical hashing for BFS
deduplication.

---

## Environment

- **Python:** `/home/azathoth/miniconda3/bin/python3` (CPython 3.13.12),
  managed via `uv`.
- **Venv:** `./.venv/` (created by `uv venv` on 2026-04-20).
  `.gitignore`'d.
- **Install command:** `uv pip install -e '.[dev]'`. Currently
  installed: apeiron (editable), numpy 2.4.4, scipy 1.17.1,
  sympy 1.14.0, networkx 3.6.1, pytest 9.0.3, pytest-cov 7.1.0,
  mypy 1.20.1, ruff 0.15.11, coverage 7.13.5, plus transitive deps.
- **Test runner:** `.venv/bin/python3 -m pytest` (or
  `source .venv/bin/activate` + `pytest`). Commits 1–8 were written
  against system pytest 7.4.4; the suite was re-verified under venv
  pytest 9.0.3 on 2026-04-20 before commit C started. 173/173 pass
  in 0.29 s under the new runner.

---

## Commits

Reverse-chronological. Authoritative log is `git log`; this list is for
quick orientation.

Reverse-chronological. Authoritative log is `git log`; this list is for
quick orientation.

- **15** (pending) — `feat(substitution)`: SubstitutionRule,
  substitution_matrix, is_primitive, perron_frobenius_in_zphi.
- **14** `3822783` — `docs`: decisions log + post-E STATUS update.
- **13** `fba7ae7` — `test(integration)`: rhombic triacontahedron
  fixture and I-transitivity.
- **12** `dd36f96` — `feat(util)`: reject interior input vertices by
  default in `load_candidate`.
- **11** `2847b34` — `feat(util)`: `load_candidate` JSON loader with
  `scale_denom` schema.
- **10** `7b07848` — `feat(polyhedron)`: hull-then-merge constructor
  via scipy oracle, exact ZPhi validation, coplanar face merge.
- **9** `1c1e22a` — `docs`: add STATUS.md project status file.
- **8** `6c9e6e0` — `feat(polyhedron)`: Polyhedron dataclass, canonical
  form, isometry action.
- **7** `2848146` — `feat(polyhedron)`: validation predicates
  (orientation, interiority, Euler).
- **6** `1d0d0fa` — `feat(symmetry)`: icosahedral group I over ℤ[φ]
  with ×2 storage.
- **5** `eb8c537` — `docs`: CLAUDE.md §3 — matrix/vertex entries in
  (½)ℤ[φ] under ×2 storage.
- **4** `10524d9` — `docs(tests)`: integration test convention.
- **3** `ca66913` — `docs`: CLAUDE.md §5 acceptance criteria rewrite.
- **2** `a472a5e` — `feat(zphi)`: exact ℤ[φ] arithmetic.
- **1** `a678272` — `chore`: scaffold repo layout.

Test totals (pre-commit-15 working tree): 277 passing in 1.89 s under
venv pytest 9.0.3.

---

## Module sequence (CLAUDE.md §5.1)

- **1 `zphi.py`** — done in `a472a5e`.
- **2 `symmetry.py`** — done in `1d0d0fa`.
- **3 `polyhedron.py`** — done via five sub-commits A–E, details below.
- **4 `substitution.py`** — done in commit 15 (pending). Penrose P3
  and Fibonacci acceptance oracles both pass.
- **5 `corona.py`** — current focus.
- **6 `hierarchy.py`** — not started.

## Implementation plan — `polyhedron.py`

Five-sub-commit breakdown, decided 2026-04-19:

- **A** — validation predicates (`exact_orientation`, `is_in_hull`,
  `check_euler`). **Done** in `2848146`.
- **B** — `Polyhedron` dataclass, canonical form, isometry action.
  **Done** in `6c9e6e0`.
- **C** — hull-then-merge constructor: scipy oracle, exact validation,
  coplanar face merge. **Done** in `7b07848`.
- **D** — `load_candidate(path) -> Polyhedron` in `util.py` with JSON
  schema carrying `scale_denom`. **Done** in `2847b34`; strict-mode
  flag added in `dd36f96` (rejects interior input vertices by default
  to catch vertex-set spec errors at load time).
- **E** — integration test: rhombic triacontahedron fixture and
  I-transitivity on faces. **Done** in `fba7ae7`, 14 tests across
  four classes (combinatorics, rhombus shape, icosahedral symmetry,
  face transitivity, ring correctness).

---

## Upstream relays with Claude (web)

Claude (web) authors `CLAUDE.md` and is the architectural upstream;
Shane relays Q&A. See the collab-relay memory in
`~/.claude/projects/-home-azathoth-Apeiron/memory/`.

- **2026-04-19 — Is the repo greenfield?** (§5 of CLAUDE.md claimed
  partially-built code; disk was empty.) **Resolved** — §5 was
  forward-looking; rewritten as acceptance criteria in `ca66913`.
- **2026-04-19 — Build backend** (hatchling vs uv vs poetry).
  **Resolved** — hatchling as PEP-517 backend, uv for env/dep
  management.
- **2026-04-19 — Rhombic triacontahedron vertex convention.**
  **Resolved** — 8 cube vertices, 12 dodecahedral cyclic perms,
  12 icosahedral cyclic perms; 32 vertices / 60 edges / 30 rhombi.
- **2026-04-19 — Where does the RTH I-transitivity test live?**
  **Resolved** — `tests/integration/` per the ≥-2-core-imports rule,
  codified in `tests/integration/README.md` (`10524d9`).
- **2026-04-19 — Icosahedral group matrices: ℤ[φ] or (½)ℤ[φ]?**
  **Resolved** — (½)ℤ[φ]; denominator globally bounded by 2 via the
  I-invariant lattice L. Stored as 2·g numerators under implicit
  denom 2. CLAUDE.md §3 rewritten in `eb8c537`.
- **2026-04-20 — Generator choice for `symmetry.py`.** Non-blocking,
  flagged for later relay if Claude (web) has a preferred canonical
  triple. Current choice: ROT_5 about (0, 1, φ), ROT_3 about
  (1, 1, 1), ROT_2 about the x-axis.
- **2026-04-20 — Convex-hull strategy: scipy or from-scratch?**
  **Resolved** — scipy as a *combinatorial* oracle (index triples
  only, floats discarded at the boundary), followed by three exact
  ZPhi validation predicates: orientation, interiority, Euler.
  Fallback to exact hull required if validation fires; currently a
  stub (`NotImplementedError`).
- **2026-04-20 — On-disk vertex convention for candidate files.**
  **Resolved** — ℤ[φ] author-form on disk (vertices as
  `[[a, b], [a, b], [a, b]]` triples), with `scale_denom` (default 1;
  only 2 otherwise permitted) handled at load time. Loader
  `load_candidate(path) -> Polyhedron` lives in `util.py`, not
  `polyhedron.py` — the latter stays pure geometry.
- **2026-04-20 — RTH vertex coordinates.** First draft (cube ×φ /
  dodec ×φ / icos ×φ all sharing cyclic-perm axes at z = ±φ²) gave a
  20-vertex convex hull, not 32 — the dodec vertices were strictly
  interior to edges between icos vertices. **Resolved** — corrected
  spec uses `(0, ±φ, ±1/φ)` for non-cube dodec positions (different
  z-levels from the icos set), landing all 32 intended vertices as
  hull extreme points. Full rationale in `docs/decisions.md`.
- **2026-04-20 — Strict-mode flag on the candidate loader.** Should
  `load_candidate` reject inputs that scipy drops as interior?
  **Resolved** — yes, default-on (`allow_interior_inputs=False`).
  The RTH-first-draft failure is exactly the silent-wrong-polyhedron
  case that needs loud surfacing at load time.
- **2026-04-20 — Acceptance oracle for `substitution.py`.** Which 2D
  Penrose-style rule to use as the PF / primitivity oracle?
  **Resolved** — Penrose P3 (thick/thin rhombus, PF = φ²) as primary
  for structural similarity to the 3D target (rhombic faces, φ-edge
  ratios); Fibonacci (PF = φ) as secondary diagnostic specifically to
  catch any code that hard-codes φ² as "the" canonical PF.

---

## Attempts that failed and what we did instead

### scipy install via `/usr/bin/python3 -m pip`

**Tried** (2026-04-20): `python3 -m pip install --user 'scipy>=1.11'`
against `/usr/bin/python3`.

**Failed** with `No module named pip`. The system Python doesn't have
pip reachable; the apt-packaged `python3-scipy` exists (1.11.4) but
installing it would need sudo and would be a system-level change.

**Doing instead:** set up a project-local uv venv
(`uv venv && uv pip install -e '.[dev]'`), added `scipy>=1.11` to
`pyproject.toml`, re-synced — scipy 1.17.1 is now installed in the
venv. All test runs switched to `.venv/bin/python3 -m pytest`.

### `_` discard-variable in `test_symmetry.py`

**Tried** (2026-04-19): `with pytest.raises(IndexError): _ = v[3]`
to catch an expected IndexError.

**Failed** lint (ruff: unused-variable on `_`).

**Doing instead:** bare expression statement `with pytest.raises(...): v[3]`.

### `_raw_tetrahedron()[0]` in direct `Polyhedron(...)` construction tests

**Tried** (2026-04-20, commit B development): use the raw (non-
canonical) tetrahedron vertex list to test face-level rejection
branches in `Polyhedron.__post_init__`.

**Failed** because the vertex-canonical-order check fires *before*
the face-level checks, so none of the five face-rejection tests
hit their intended branch — they all raised the wrong error first.

**Doing instead:** added a `_canonical_tetra_vertices()` helper
returning vertices pre-sorted by the canonical key; the face-rejection
tests use this helper so they bypass the vertex-order gate and hit
the face-level gates as intended.

### Trusting scipy's simplex vertex order to be outward-oriented

**Tried** (2026-04-20, commit C development): use the simplex index
triples returned by `scipy.spatial.ConvexHull` directly — under the
assumption that `(b - a) × (c - a)` applied to each simplex gives the
outward normal consistently.

**Failed** on the stored-×2 unit cube: `from_vertices` produced 12
triangular faces instead of 6 quadrilaterals. Diagnosis — scipy
returns simplices with mixed orientations (some outward-facing, some
inward-facing) even within a single hull face. scipy's `hull.equations`
attribute is the authoritative outward normal per simplex, but it's
float metadata that we deliberately don't consume.

**Doing instead:** added `_orient_simplices_outward`, a pure-Z[φ]
reorientation pass that runs immediately after scipy and flips any
simplex whose exact normal points toward the hull centroid. Uses the
unscaled centroid direction `sum_v - n_count * a` (positive-scalar
multiple of `centroid - a`) to avoid division, and raises
`ValueError` if the centroid happens to be coplanar with a face
(which is impossible for a valid convex hull). After reorientation,
`_same_oriented_plane` correctly groups same-face simplices and the
coplanar-merge produces one quadrilateral per cube face.

### Claiming matrices of I have entries in ℤ[φ] (Cartesian)

**Tried** (2026-04-19, initial §3 of CLAUDE.md): state that every
g ∈ I is a SL(3, ℤ[φ]) matrix in the standard Cartesian basis.

**Failed** mathematically — the 5-fold rotation about (0, 1, φ) has
entries in (½)ℤ[φ]; the half-factors from cos(72°) = (φ−1)/2 and the
Rodrigues denominator do not cancel. `det(S) = 1` verified independently.

**Doing instead:** storage convention is ×2 scaling everywhere
(CLAUDE.md §3.2 and §3.3 as of `eb8c537`). `Rotation.matrix` holds
2·g with pure ℤ[φ] entries; `apply`/`compose` discharge the implicit
denominator 2 via `_halve_zphi`, which asserts exactness and
fails loudly on any parity violation. The L-invariant lattice (e₁ =
(1, φ, 0), e₂ = (φ, 0, 1), e₃ = (0, 1, φ); det V = −2φ²) is the
reason the denominator never grows under composition.

---

## Open questions / things to watch

- **Exact-arithmetic hull fallback.** Commit C ships a
  `NotImplementedError` stub for the case where scipy's validation
  fails. We don't hit this path for the standard fixtures
  (tetrahedron, cube, RTH) but the stub needs replacing before any
  candidate tile with near-degenerate geometry is attempted.
- **CLAUDE.md §7.1 phrasing.** Says "No SciPy for core arithmetic".
  Adding scipy as a combinatorial oracle is compatible with the
  *intent* (scipy is not used for ℤ[φ] arithmetic) but the letter of
  the rule reads more broadly. Proposed resolution deferred until after
  commit C: relay to Claude (web) for a one-sentence §7.1 clarification
  at the same time as any other CLAUDE.md edits we accumulate.
- **Canonical face-cycle direction.** Current canonicalisation rotates
  (preserving direction), which assumes the caller has provided
  outward-oriented faces. scipy's output is documented as outward-
  normal consistent; if this ever turns out to be wrong we'll see it
  immediately via `is_in_hull` rejecting all inputs.
