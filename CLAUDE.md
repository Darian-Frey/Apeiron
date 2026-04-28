# CLAUDE.md — Apeiron

**Handoff document for Claude Code.** Read this in full before touching any code in this repository. It defines the problem, the algebraic setting, the architecture, the conventions, and — critically — what is currently broken and what to prioritise next.

---

## 1. Project identity

**Apeiron** is a research framework for the **3D Einstein problem**: the search for a single polyhedron *P* ⊂ ℝ³ such that *every* tiling of ℝ³ by isometric copies of *P* is non-periodic.

This is not an application. It is not a product. It is a mathematical research codebase whose purpose is to:

1. Provide **exact-arithmetic primitives** for reasoning about polyhedra with vertices in ℤ[φ]³.
2. Enumerate and verify **corona configurations** (a tile surrounded by its neighbours) under icosahedral symmetry.
3. Specify and test **substitution rules** (inflation maps) over ℤ[φ]-modules.
4. Serve as the computational backbone for a multi-year research programme aimed at either constructing a 3D einstein or proving a no-go theorem.

The author (Shane Hartley, GitHub: `Darian-Frey`) is treating this as a 2–4 year programme. Code quality matters; speed does not.

---

## 2. Problem statement (precise)

### 2.1 Target definitions

We work with three increasingly strong notions of aperiodicity. Every function, test, and proof obligation in this codebase must be tagged with which definition it is establishing.

- **Weak aperiodic:** there exists at least one non-periodic tiling by *P*.
- **Strong aperiodic (Einstein):** every tiling by *P* is non-periodic.
- **Strict Einstein:** strong, using only orientation-preserving isometries of ℝ³ (SO(3) plus translations — no reflections, no screw displacements).

**The target of Apeiron is the Strong Einstein definition**, with the Strict Einstein
version as a stretch goal.

Track A (Danzer deformation, §6.1) starts from a baseline that is Strong but not Strict:
the published ABCK substitution uses orientation-reversing isometries (I_h, order 120),
confirmed by both Koca et al. (arXiv 2003.13449, eqs. 19–31) and Paolini's POV-Ray
source (`danzer.dmf.unicatt.it`). Any single-tile Einstein candidate derived by merging
ABCK tiles inherits this I_h structure unless the deformation specifically eliminates all
orientation-reversing children. Strict Einstein status for a Track A candidate is
achievable only if such elimination is demonstrated; it is not guaranteed by the merge
alone and should be treated as a secondary check, not an assumption.

Track B (substitution-first, §6.2) begins with no constraint on I vs I_h and may produce
candidates of either type. The I-only (Strict) constraint should be imposed explicitly
when searching for Track B substitution rules if Strict Einstein is the goal.

The Schmitt-Conway-Danzer biprism is excluded from consideration not because it uses
I_h but because it satisfies only the Weak definition (every tiling is non-periodic, but
some tilings exist only via screw symmetry rather than being forced aperiodic).

### 2.2 What counts as a tile

A **tile** is a closed, topologically-ball polyhedron *P* ⊂ ℝ³ with:

- Finitely many flat faces, each a simple polygon.
- All vertices in ℤ[φ]³ (see §3 for why this is required).
- Interior connected; ∂*P* a topological 2-sphere.

A **tiling** by *P* is a covering of ℝ³ by isometric copies of *P* with pairwise interior-disjoint interiors, face-to-face unless explicitly noted otherwise.

### 2.3 Why this is hard in 3D (for context)

- No tractable cyclotomic analogue exists in 3D. ℤ[ζ₈] for Ammann-Beenker and ℤ[ζ₅] for Penrose both exploit 2D-specific structure; the closest 3D analogue is ℤ[φ] with icosahedral symmetry, which is *weaker* than cyclotomic.
- Corona combinatorics explode: a planar tile has O(10¹–10²) distinct local configurations; a 3D tile routinely has O(10⁴–10⁶) even for modest face counts.
- The 2023 hat/spectre proof leaned on a computer-assisted case analysis of local configurations. Scaling that approach to 3D is an architectural question, not just a compute question.

---

## 3. Algebraic framework

### 3.1 The coordinate ring: ℤ[φ]

Let φ = (1 + √5) / 2, the golden ratio. We work over the ring ℤ[φ] = ℤ + ℤφ.

**Why ℤ[φ] and not ℤ[ζ₅] or a 6D cut-and-project lattice?**

- Icosahedral symmetry (the highest-symmetry rotation group of ℝ³ other than SO(3) itself) has eigenvalues in ℚ(√5).
- The natural inflation factor for icosahedral substitution tilings is **φ²**, which lives in ℤ[φ] with integer coefficients.
- ℤ[φ] is a PID (in fact the ring of integers of ℚ(√5)), making exact arithmetic and equality testing trivial.
- Cut-and-project from ℤ⁶ into ℝ³ via the canonical icosahedral slope lands vertices in ℤ[φ]³. Using ℤ[φ]³ directly bypasses the lift-and-project machinery.

### 3.2 Representation

All coordinates are triples (a, b) with a, b ∈ ℤ, representing a + bφ. See `zphi.py` for the canonical implementation.

Tile vertices formally live in **(½)ℤ[φ]³**, not in ℤ[φ]³. The factor of 2 is the minimal denominator that makes the icosahedral rotation group I act integrally in the standard Cartesian basis (see §3.3). It is globally bounded: no finer denominator is ever needed, at any stage of the pipeline.

**Storage convention: ×2 scaling at the data layer.** Every vertex is stored as 2·v ∈ ℤ[φ]³, so the storage type is pure `Mat1[ZPhi]` with no denominator field. A vertex that appears in mathematical text as (0, 1, φ) is stored as (0, 2, 2φ). The factor of 2 surfaces only in docstrings and in human-readable coordinate accessors; plane equations, face-merge coplanarity tests, edge equality, and all corona incidence arithmetic happen in pure ℤ[φ] with no denominators anywhere in the verification pipeline.

**Do not use floats for tile geometry under any circumstances.** Floats are acceptable only in `viz/` for rendering and in numerical sanity checks — never in the verification pipeline.

### 3.3 Symmetry

The relevant group is the **icosahedral rotation group I** (order 60), generated by rotations of orders 2, 3, and 5. The full icosahedral group I_h (order 120) includes reflections; for Strict Einstein work, only I is used.

**I does not act by SL(3, ℤ[φ]) matrices in the standard Cartesian basis.** The 5-fold rotation about (0, 1, φ) has entries in (½)ℤ[φ] (the (φ−1)/2 and 1/2 factors from cos(72°) and the Rodrigues denominator do not cancel). There are two ways to recover exact arithmetic — a change of basis, or a controlled denominator — and Apeiron uses the second.

**The denominator is globally bounded by 2.** Let L ⊂ ℝ³ be the ℤ[φ]-lattice spanned by e₁ = (1, φ, 0), e₂ = (φ, 0, 1), e₃ = (0, 1, φ). L contains all 12 icosahedron vertices (cyclic permutations of (0, ±1, ±φ)) and is I-invariant, so I ⊂ SL(3, ℤ[φ]) when expressed in the L-basis. Change-of-basis back to Cartesian uses the matrix V with columns e₁, e₂, e₃; det V = −2φ². Its inverse has entries in (1/det V)·ℤ[φ] = (2−φ)/2 · ℤ[φ] = (½)ℤ[φ] (since (2−φ) + (φ−1) = 1 generates ℤ over ℤ[φ]). Therefore every g ∈ I satisfies [g]_Cartesian ∈ (½) M₃(ℤ[φ]), uniformly in g.

Consequently, composition within I does not grow the denominator: a product of two (½)ℤ[φ] matrices lies a priori in (¼)ℤ[φ], but must equal some h ∈ I ⊂ (½)ℤ[φ], so the numerator of the product is always even and the division by 2 is exact.

**Storage.** `Rotation` is a `dataclass(frozen=True)` with a single field `matrix: Mat3[ZPhi]` — the numerator of the matrix under a fixed, class-level implicit denominator of 2. The denominator is not a parameter of the type and is never exposed at the API boundary. The 60 matrices of I are precomputed and verified once at import time (`|I| = 60`, closure under composition, inverse closure, norm preservation).

**Action on positions.** Positions stored as 2v (§3.2) compose with rotations stored as 2·[g]_Cartesian via `(2v)' = [g]_num · (2v) / 2`, where the division is always exact for the same parity reason as composition above. In practice this is the one place the factor of 2 appears in the arithmetic, and it is an exact integer-division check — the test that it is always exact is itself a correctness invariant of the pipeline.

### 3.4 Substitution

A **substitution rule** σ on *P* is a linear map σ: ℝ³ → ℝ³ with eigenvalue λ > 1 such that σ(*P*) admits a dissection into isometric copies of *P*. For icosahedral systems:

- λ = φ² is the canonical choice.
- σ acts on ℤ[φ]³ as a ℤ[φ]-linear map; its matrix in the standard basis has entries in ℤ[φ].
- The **substitution matrix** *M* counts how many copies of each prototile type (here, orientation class of *P*) appear in σ(*P*). For a monotile with *n* orientation classes, *M* is an *n* × *n* non-negative integer matrix.
- **Primitivity** of *M* (some power has all-positive entries) implies a unique Perron-Frobenius eigenvalue, which is the necessary algebraic witness for aperiodicity of the substitution hierarchy.

---

## 4. Repository layout

```
Apeiron/
├── CLAUDE.md                     # this file
├── README.md                     # human-facing project description
├── pyproject.toml                # build/dependency config (uv or poetry)
├── apeiron/                      # core package
│   ├── __init__.py
│   ├── zphi.py                   # ℤ[φ] exact arithmetic
│   ├── symmetry.py               # icosahedral group I (and I_h) matrices
│   ├── polyhedron.py             # Polyhedron class, face/edge/vertex ops
│   ├── substitution.py           # Substitution rules and matrices
│   ├── corona.py                 # Corona BFS engine (verification core)
│   ├── hierarchy.py              # Supertile construction and recognisability
│   └── util.py                   # Hashing, canonical forms, small helpers
├── tests/                        # pytest suite, one file per module
│   ├── test_zphi.py
│   ├── test_symmetry.py
│   ├── test_polyhedron.py
│   ├── test_substitution.py
│   ├── test_corona.py
│   └── test_hierarchy.py
├── candidates/                   # candidate tile definitions (JSON or Python)
│   └── README.md
├── notebooks/                    # exploration and visualisation
└── docs/
    ├── algebra.md                # expanded notes on ℤ[φ], icosahedral action
    ├── roadmap.md                # research phases and milestones
    └── literature.md             # annotated bibliography
```

The seven core modules (`zphi`, `symmetry`, `polyhedron`, `substitution`, `corona`, `hierarchy`, `util`) are intentionally kept flat. Do not introduce subpackages without discussion.

---

## 5. Implementation order and known pitfalls

This repo is greenfield. The seven-module layout in §4 is implemented in a specific order, and each module has defined acceptance criteria. A module "ships" when all criteria are met and its tests pass — not when the code compiles.

### 5.1 Implementation order and acceptance criteria

Implement strictly in this order. Each module depends on the previous ones being stable.

**1. `zphi.py`** — exact arithmetic over ℤ[φ].

Acceptance:

- `ZPhi(a: int, b: int)` dataclass representing a + bφ, frozen and hashable.
- Operators: `+`, `-`, `*`, unary `-`, equality. Division is explicitly not supported.
- `.conjugate()` method (Galois conjugation: a + bφ ↦ a + b(1−φ) = (a+b) − bφ).
- `.norm()` method for debugging (integer-valued).
- Tests: ring axioms (associativity, distributivity), φ² = 1 + φ, conjugation involution, closure under arithmetic, hash stability.
- `mypy --strict` passes. 100% line coverage.

**2. `symmetry.py`** — icosahedral rotation group I.

Acceptance:

- 60 rotation matrices, each a 3×3 matrix over ZPhi, constructed from generators (rotations of orders 2, 3, 5).
- Group closure verified in a test (the set is closed under composition).
- Identity present; every element has its inverse in the set.
- Action on ZPhi³ preserves the standard quadratic form (norm-preserving).
- Tests: |I| = 60, closure under composition, inverse closure, no floating-point entries anywhere.

**3. `polyhedron.py`** — polyhedra with ℤ[φ]³ vertices. **This is where the known pitfall lives — see §5.2.**

Acceptance:

- `Polyhedron` dataclass: frozen vertex list, frozen face list (each face an ordered vertex-index cycle).
- Construction from a vertex set via convex hull, followed by a coplanar-face-merge pass (see §5.2).
- Face-merge uses exact ZPhi plane equality, never a float tolerance.
- Isometry action: `apply(g: Rotation) → Polyhedron`.
- Canonical form for hash/equality (vertex ordering normalised).
- Tests: rhombic triacontahedron fixture with 32 vertices, 60 edges, 30 rhombic faces, Euler χ = 2; rhombic dodecahedron fixture as a secondary check.

**4. `substitution.py`** — substitution rules and their matrices.

Acceptance:

- `SubstitutionRule` dataclass: linear map (ZPhi-valued 3×3 matrix) plus a dissection specification mapping the inflated prototile to a multiset of positioned copies.
- `substitution_matrix(rule) → np.ndarray` returning the n×n non-negative integer matrix.
- `is_primitive(M)` check (some power of M has all-positive entries).
- Perron-Frobenius eigenvalue extraction (as an algebraic number in ℚ(√5) where possible).
- Tests: a known 2D Penrose-style rule as an oracle (primitive, PF eigenvalue φ²).

**5. `corona.py`** — corona BFS engine.

Acceptance:

- `CoronaConfig` dataclass: central tile plus a set of positioned neighbour tiles.
- `corona_1(P)`: enumerate all complete first-shell configurations modulo I.
- `corona_2(config)`: extend to second shell.
- Canonical hashing of configurations for BFS deduplication.
- Tests: cube has a unique `corona_1` (itself, face-to-face); rhombic dodecahedron has a unique `corona_1`. Both are sanity-check oracles from periodic tilings — not proof obligations for the Einstein problem.

**6. `hierarchy.py`** — supertile construction and recognisability.

Acceptance:

- Supertile construction from a substitution rule.
- Recognisability test: given a tile, determine its supertile membership from a bounded local neighbourhood.
- Every function tagged with which of the four pillars (§6.3) it establishes.
- Tests: oracle checks against known 2D hierarchies where recognisability is established in the literature.

### 5.2 Predicted pitfall: coplanar face merging in `polyhedron.py`

When constructing a Polyhedron from a vertex set via convex hull, standard hull algorithms emit **triangulated** faces. For tiles whose natural faces are non-triangular — which is almost all icosahedral-compatible tiles (pentagons, rhombi, decagons) — adjacent coplanar triangles that should be *one face* arrive as separate faces.

**Why this is predicted, not observed:** prior (un-persisted) implementations of similar code ran into this exact failure mode. Nothing in the current repo has been attempted yet.

**Downstream effect if unfixed:** the CoronaEngine BFS in `corona.py` cannot run correctly — face-to-face adjacency is ill-defined when "faces" are triangulation artefacts. Two tiles sharing a pentagonal face appear as sharing three triangular sub-faces, breaking adjacency counts, corona completion tests, and the eventual supertile recognisability checks.

**Required approach:** run a coplanar-face-merge pass immediately after hull construction:

1. Group triangles whose supporting planes agree exactly (exact ZPhi plane equality, not a float tolerance — non-negotiable).
2. Within each group, compute the union polygon by edge-matching — triangles sharing an edge with the same pair of vertices collapse that edge.
3. Recover the ordered vertex cycle of the merged face.
4. Verify the merged face is a simple polygon (no self-intersections, single connected boundary).

The rhombic triacontahedron fixture is the acceptance test: 30 rhombic faces, not 60 triangulated sub-faces.

### 5.3 Prior work that informs this

There is prior 2D work on Ammann-Beenker (ℤ[ζ₈]) and an SD-14 candidate tile from a separate project. **That code lives elsewhere and is not part of Apeiron.** Reference it only for algorithmic patterns, not for import. The 3D framework is a clean re-architecture.

---

## 6. Research strategy (two parallel tracks)

Both tracks should be pursued simultaneously. A candidate from either unblocks the full verification pipeline.

### 6.1 Track A — Deformation-first

Start from a known 3D aperiodic tile *set* (Danzer's 4-tile set is canonical) and attempt to:

- Identify minimal-deformation paths that merge multiple tiles into one.
- Check whether the merged tile retains the substitution property.
- Use interval arithmetic over ℤ[φ] to systematically explore parameterised tile families.

**Risk:** may be provably impossible — the 4-tile structure of Danzer may be essential.

### 6.2 Track B — Substitution-first

Define an abstract substitution rule σ with eigenvalue φ² on an *n*-tile alphabet (start with *n* = 1), and:

- Search for ℤ[φ]-linear σ with primitive substitution matrix.
- For each candidate σ, attempt to construct a geometric realisation — an actual polyhedron *P* whose inflation σ(*P*) dissects into isometric copies of *P*.

**Risk:** the algebraic candidate σ may have no geometric realisation.

### 6.3 The four pillars of a proof

Any eventual proof that a candidate *P* is a Strong Einstein decomposes into four pillars:

1. **Substitution exists and is primitive.** Standard linear algebra over ℤ[φ].
2. **Recognisability (border forcing).** Every tile's supertile membership is determined by a bounded local neighbourhood. This is historically the hard pillar.
3. **Aperiodicity from recognisability.** Standard argument: a period vector *v* implies σⁿ(*v*) is also a period, contradiction by finite local complexity.
4. **No non-hierarchical tiling exists.** The genuinely Einstein step — showing no tiling by *P* escapes the substitution hierarchy. In 2D this is where the hat proof does its heavy case analysis.

Code in `hierarchy.py` should be organised explicitly around these four pillars, with each tagged with which pillar it is establishing.

---

## 7. Development conventions

### 7.1 Language and dependencies

- **Python 3.11+** (we use structural pattern matching and the latest typing features).
- Core dependencies: `numpy`, `sympy` (for occasional symbolic checks), `networkx` (for corona adjacency graphs), `pytest`. No SciPy for core arithmetic — we roll exact ℤ[φ] ourselves.
- Visualisation (in `notebooks/` and an optional `viz/` module): `pyvista` or `plotly` for 3D, `matplotlib` for 2D projections. **Never** in the core verification path.

### 7.2 Type discipline

- Full type hints on every public function. `mypy --strict` should pass.
- The core ℤ[φ] element is a `ZPhi` dataclass, not a tuple. Never accept raw tuples at API boundaries.
- `Polyhedron`, `Face`, `Edge`, `Vertex` are all dataclasses with `frozen=True` where possible for hashability.

### 7.3 Exactness discipline

- **No floats in the verification pipeline.** Ever. Equality of tile positions, planes, face normals, and edge directions must be decidable exactly in ℤ[φ].
- Any function that takes or returns a float outside `viz/` is a bug.
- Numerical sanity checks (e.g. "does this match what `scipy.spatial.ConvexHull` would return in floats?") are fine as *additional* tests, never as the primary check.

### 7.4 Testing

- `pytest` with `pytest-cov`; aim for ≥ 90% line coverage on all core modules.
- Every bug fix ships with a regression test.
- Each module has a corresponding `test_*.py`. Cross-module integration tests live in `tests/integration/`.
- For the hull merge fix specifically: use a rhombic triacontahedron and a rhombic dodecahedron as fixtures. Both have known face structures.

### 7.5 Style

- `ruff` for linting and formatting with default config plus line length 100.
- Docstrings: numpy style.
- Commits: conventional-commits (`fix:`, `feat:`, `refactor:`, `test:`, `docs:`).
- No `TODO` comments without a GitHub issue number.

### 7.6 Performance

This is explicitly a correctness-first codebase. Do not optimise before the verification pipeline runs end-to-end on a trivial test tile. When optimisation does come, it will likely be:

- Canonical-form hashing of tile positions (for BFS deduplication).
- Precomputed symmetry orbits.
- Parallel corona exploration via `multiprocessing`.

Not sooner.

---

## 8. Immediate priorities (in order)

1. **Fix the coplanar face merging bug in `polyhedron.py`.** This unblocks everything. See §5.1.
2. **Write the rhombic triacontahedron test fixture.** It's the canonical icosahedral-symmetric polyhedron and will anchor every downstream test.
3. **Scaffold `substitution.py`** — define the `SubstitutionRule` dataclass and the `substitution_matrix` / `is_primitive` functions. No geometric realisation yet, just the algebraic layer.
4. **Scaffold `corona.py`** — define `CoronaConfig`, `corona_1`, `corona_2` computation. A first-corona completion test on a trivial periodic tile (cube) should pass as a sanity check.
5. **Write `docs/roadmap.md`** — expand §6 into a proper phase-by-phase research plan with time estimates.

Do **not** attempt to implement the recognisability test in `hierarchy.py` until steps 1–4 are done and stable.

---

## 9. Literature references

These should be in `docs/literature.md` with full annotations. Minimum reading for anyone contributing:

- Smith, D., Myers, J. S., Kaplan, C. S., & Goodman-Strauss, C. (2023). *An aperiodic monotile.* — the hat paper.
- Smith, D., Myers, J. S., Kaplan, C. S., & Goodman-Strauss, C. (2023). *A chiral aperiodic monotile.* — the spectre paper.
- Goodman-Strauss, C. (1998). *Matching rules and substitution tilings.* Annals of Mathematics, 147(1), 181–223. — the theoretical backbone for the four-pillar proof structure.
- Danzer, L. (1989). *Three-dimensional analogs of the planar Penrose tilings and quasicrystals.* Discrete Mathematics, 76, 1–7. — the 4-tile 3D aperiodic set.
- Socolar, J. E. S., & Taylor, J. M. (2011). *An aperiodic hexagonal tile.* Journal of Combinatorial Theory A, 118, 2207–2231. — cautionary example on non-local matching rules.
- Schmitt, P. (1988). *An aperiodic prototile in space.* Unpublished manuscript; the SCD biprism construction.

---

## 10. Communication contract for AI assistants

If you are Claude Code (or any AI assistant) working on this repo:

- **Never introduce floating-point arithmetic into the core verification path.** If a problem seems to require it, stop and raise the issue explicitly.
- **Never claim a proof step is done without the test that verifies it.** This codebase is a research artefact; false positives are worse than slow progress.
- **Flag every epistemic uncertainty.** If you're implementing a substitution rule whose primitivity you haven't verified, say so in the docstring and the commit message.
- **Prefer small, reviewable commits.** The author reviews every change. A 500-line diff is a rejected PR.
- **When in doubt about mathematical content, stop and ask.** The algebra is load-bearing; a subtle error in ℤ[φ] arithmetic or icosahedral group structure will poison every result downstream.

This is open-ended research. The goal is not to finish quickly. The goal is to be correct.
