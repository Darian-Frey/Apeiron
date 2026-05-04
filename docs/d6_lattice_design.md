# D₆ lattice arithmetic — design document

**Status.** Design document, not code. Authorised by Claude (web)
Q16b (2026-05-04, archived at
[relays/q16_2026-05-04.md](relays/q16_2026-05-04.md)) as the
concrete action during the P1 wait (institutional ILL of
Kramer-Papadopolos-Schlottmann-Zeidler 1994). Engineering follows
this design; design proceeds during the wait.

**Scope.** Specify the data structures, operations, and interfaces
needed to implement Route 1 — a computational test of the refined
structural multi-prototile conjecture (`strategy_pivot.md` §2,
sharpened in `literature_review.md` §(e)) within the
Moody-bounded D₆ CPS framework.

**Out of scope.** Code. Implementation. Routes 2 and 3 (per Q16c,
those are wrong-scope for Route 1 work). The bridging from this
design to existing Apeiron primitives is described but not
realised.

---

## §1. The mathematical setting

### 1.1 The D₆ root lattice

D₆ is the rank-6 root lattice with simple roots α_i (i = 1, ..., 6)
and the Cartan matrix:

```text
      [ 2 -1  0  0  0  0]
      [-1  2 -1  0  0  0]
C  =  [ 0 -1  2 -1  0  0]
      [ 0  0 -1  2 -1 -1]
      [ 0  0  0 -1  2  0]
      [ 0  0  0 -1  0  2]
```

A D₆ point can be expressed in two equivalent ways:

1. As `λ = Σᵢ nᵢ αᵢ` where `nᵢ ∈ ℤ`. (root-coordinate basis)
2. As `λ = Σᵢ mᵢ lᵢ` where `lᵢ` is the orthonormal Euclidean basis
   for ℝ⁶ and `mᵢ ∈ ℤ` with the constraint `Σmᵢ = even`.
   (l-basis, per Al-Siyabi-Koca-Koca eq. (1))

The l-basis is more convenient for the projection because the
projection maps `lᵢ` to explicit ℝ³ + ℝ³_⊥ vectors per
Al-Siyabi-Koca-Koca eq. (5).

### 1.2 The icosahedral subgroup

The point group of D₆ has order 2⁵ · 6! = 23040. The icosahedral
group H₃ (order 120 including reflections; order 60 for proper
rotations I) is a maximal subgroup. Generators of H₃ in terms of
D₆ reflections (Al-Siyabi-Koca-Koca p.3):

```text
R_1 = r_1 r_5,   R_2 = r_2 r_4,   R_3 = r_3 r_6
```

where rᵢ is the simple reflection associated with αᵢ.

### 1.3 The projection map

ℝ⁶ decomposes as ℝ³ ⊕ ℝ³_⊥ (the *physical* and *internal* spaces).
The decomposition is invariant under H₃ but inverts the algebraic
conjugation σ = -τ⁻¹.

Per Al-Siyabi-Koca-Koca eq. (5), the projection of l_i to
(ℝ³, ℝ³_⊥) is given by an explicit 6×6 matrix whose entries are
in ℤ[τ] (with τ = (1+√5)/2 = ZPhi(0,1) in Apeiron's storage).
The first three rows project to physical space; the last three
to internal space.

### 1.4 Cut-and-project framework

A *window* W ⊂ ℝ³_⊥ is a compact subset satisfying Moody's W1
condition (W = closure(int(W))). The model set is

```text
Λ(W) = {π_∥(λ) : λ ∈ D₆, π_⊥(λ) ∈ W}
```

where π_∥, π_⊥ are projections to physical and internal space
respectively.

For ABCK specifically (Frettlöh §1.3, Theorem 1.3), the windows
for the three vertex classes are:

- Class I: a regular dodecahedron with five-fold-star decoration
  (long edge length 4τ/√(τ+2)).
- Class II: a regular dodecahedron of edge length 2.
- Class III: the great dodecahedron (long edge length 2τ).

The full ABCK model set's window is the disjoint union of these
three.

### 1.5 The structural multi-prototile conjecture (refined)

Per `literature_review.md` §(e):

> For every D₆ cut-and-project tiling of ℝ³ with icosahedral
> symmetry where the prototile-shape structure is determined by
> the H₃ tetrahedral fundamental region (the ABCK / icosian-
> model-set family per Moody §4.1), the tiling requires more
> than one prototile type.

Route 1 tests this conjecture by enumerating *all* D₆ vertex-
class arrangements compatible with the H₃ fundamental region
and checking each for single-tile-type tessellation.

---

## §2. Data structures

### 2.1 `D6Point`

```text
@dataclass(frozen=True, slots=True)
class D6Point:
    """Point in the D₆ root lattice, l-basis storage.

    Stores the integer 6-tuple (m_1, ..., m_6) with the
    constraint sum(m) is even. Constraint enforced at
    construction.
    """
    m: tuple[int, int, int, int, int, int]

    def __post_init__(self) -> None:
        if sum(self.m) % 2 != 0:
            raise ValueError(
                "D₆ l-basis constraint violated: "
                f"Σm_i = {sum(self.m)} (must be even)"
            )
```

**Operations:**
- `__add__(self, other: D6Point) -> D6Point` — pointwise.
- `__sub__(self, other: D6Point) -> D6Point` — pointwise.
- `__neg__(self) -> D6Point` — negation.
- `__hash__` from frozen dataclass.

**Conversion:**
- `to_root_basis(self) -> tuple[int, ...]` — returns
  (n_1, ..., n_6) such that λ = Σ n_i α_i. Inverse of
  `from_root_basis`. Note that not every (m_1, ..., m_6) with
  Σm even has a clean root-basis representation; the conversion
  uses the Cartan-matrix inverse and may return rationals (rare;
  only for D₆ weight-lattice points outside the root lattice).

### 2.2 `D6Window`

```text
class D6Window(Protocol):
    """Window in the internal space ℝ³_⊥.

    Represented as a convex polytope (or union thereof for the
    multi-window case like ABCK's three vertex-class windows).
    """

    def contains(self, point: Vec3_perp) -> bool: ...
    def vertices(self) -> tuple[Vec3_perp, ...]: ...
    def faces(self) -> tuple[tuple[int, ...], ...]: ...
    # Decomposition into connected components, for multi-prototile
    # tilings:
    def components(self) -> tuple[D6Window, ...]: ...
```

`Vec3_perp` is structurally identical to existing `Vec3` (ZPhi³),
distinguished only by its semantic role as living in internal
space. Per CLAUDE.md §3 storage convention, ×2-stored.

**Concrete implementations needed:**

1. `ConvexPolyWindow(vertices: tuple[Vec3_perp, ...], faces: ...)`
   — single convex polytope window. `contains` via halfspace
   intersection (existing `_point_in_convex_polyhedron` from
   `apeiron.track_b.realisation` should adapt directly with the
   ZPhi³ → Vec3_perp type rename).
2. `UnionWindow(components: tuple[ConvexPolyWindow, ...])` —
   disjoint union; `contains` is logical-or over components.
   ABCK's window is `UnionWindow([class_I, class_II, class_III])`.

### 2.3 `D6Lattice`

```text
@dataclass(frozen=True, slots=True)
class D6Lattice:
    """The D₆ root lattice with its icosahedral subgroup."""

    # Projection matrices: 6×6 matrices over ZPhi.
    # First 3 rows project to physical ℝ³;
    # last 3 rows project to internal ℝ³_⊥.
    projection: Mat6_ZPhi

    # The 60 (or 120) elements of I (or I_h) acting on ℝ⁶.
    icosahedral_subgroup: tuple[Mat6_int, ...]
```

The `projection` matrix is the explicit ZPhi-valued 6×6 from
Al-Siyabi-Koca-Koca eq. (5). It is fixed (one canonical D₆
embedding); not a parameter.

The `icosahedral_subgroup` is the 60 proper rotations of I (or
120 of I_h) realised as integer-valued 6×6 matrices acting on the
l-basis. These compose with `projection` to produce the
icosahedral action on (ℝ³, ℝ³_⊥).

---

## §3. Operations

### 3.1 Projection: D6Point → (Vec3_phys, Vec3_perp)

```text
def project(p: D6Point, lattice: D6Lattice) -> tuple[Vec3, Vec3_perp]:
    """Apply the canonical D₆ → ℝ³ + ℝ³_⊥ projection.

    Returns (physical, internal) where:
        physical = π_∥(p) ∈ ℝ³
        internal = π_⊥(p) ∈ ℝ³_⊥
    """
    full = lattice.projection @ p.as_vec6()  # 6-vector over ZPhi
    return (
        Vec3(full[0], full[1], full[2]),       # physical
        Vec3_perp(full[3], full[4], full[5]),  # internal
    )
```

**Storage convention.** D6Point stores raw integers. Projection
introduces ZPhi entries via the projection matrix. Output Vec3 +
Vec3_perp use the existing ×2 storage convention (CLAUDE.md §3.2).
The projection matrix is pre-scaled by 2 to make the output
storage form direct.

### 3.2 Membership in the model set

```text
def is_in_model_set(p: D6Point, window: D6Window,
                    lattice: D6Lattice) -> bool:
    """Is p's projection in Λ(W)?"""
    _, internal = project(p, lattice)
    return window.contains(internal)
```

### 3.3 Nearest D₆ point to a ℝ³ position

This is the **dual** problem: given a position v ∈ ℝ³ (in
Apeiron's ZPhi³ representation), find the D₆ point whose
physical projection equals v (if any).

```text
def nearest_d6_point(v: Vec3, lattice: D6Lattice) -> D6Point | None:
    """Returns the D₆ point p with π_∥(p) = v, or None.

    The map π_∥: D₆ → L (the projected lattice) is many-to-one
    in general, but for icosahedral D₆ in the standard embedding
    it's one-to-one (per Al-Siyabi-Koca-Koca §2).

    Strategy: solve the linear system
        projection_phys @ m = v
    over ℤ. Since `projection_phys` has ZPhi entries, this
    requires lifting to a 6D ZPhi-linear system; solution is
    integer iff v ∈ L.
    """
```

Note: existence depends on whether v is in the projected lattice
L = π_∥(D₆). Not every ZPhi³ point is in L.

### 3.4 Orbit under I_h

```text
def orbit_under_ih(p: D6Point, lattice: D6Lattice) -> set[D6Point]:
    """The 120 (or fewer if p is fixed) images of p under I_h."""
    return {D6Point(g @ p.as_vec6_int())
            for g in lattice.icosahedral_subgroup}
```

### 3.5 Tile-type partitioning (the load-bearing operation for Route 1)

Given a model set Λ(W) with a multi-component window
W = W_1 ∪ ... ∪ W_k, every model-set point lies in exactly one
W_i (by disjointness). The point's "vertex class" is the index i.

For ABCK: 3 vertex classes (I, II, III), corresponding to the
three windows from Frettlöh §1.3.

The **prototile structure** of the tiling is determined by the
combinatorial arrangement of vertex classes around each tile, plus
the H₃ fundamental-region cell structure.

```text
def tile_types_from_window(window: D6Window,
                            lattice: D6Lattice,
                            radius: ZPhi) -> tuple[Polyhedron, ...]:
    """Compute the prototile set of Λ(W) restricted to a ball.

    Strategy:
    1. Enumerate all D₆ points whose physical projection is within
       `radius` of the origin.
    2. Filter to those whose internal projection is in W.
    3. Build the H₃-fundamental-region tessellation around each
       lattice point.
    4. Group resulting tiles by similarity class (existing
       `apeiron.track_b.taxonomy.is_h3_compatible` machinery).
    5. Return one representative per class.
    """
```

This is the function whose output cardinality the conjecture
claims is always > 1 for icosahedral D₆ CPS.

---

## §4. Interface to existing Apeiron primitives

### 4.1 ZPhi arithmetic

`apeiron.zphi.ZPhi` is reused unchanged. The 6×6 projection matrix
is `Mat6_ZPhi`, a new type analogous to existing `Mat3` from
`apeiron.symmetry`.

### 4.2 Vec3, Mat3, Polyhedron

`Vec3` is reused for both physical and internal-space points. A
phantom-type `Vec3_perp` aliasing `Vec3` is recommended only if
mypy strictness is desired; otherwise use Vec3 directly with a
runtime-only `space: Literal["phys", "perp"]` tag for debug logs.

`Polyhedron` is reused unchanged. The tile-type partitioning in
§3.5 returns `tuple[Polyhedron, ...]`.

### 4.3 Icosahedral group

`apeiron.symmetry.ICOSAHEDRAL` is reused for the 60-element
proper-rotation list in **physical** space. The new
`icosahedral_subgroup` field of `D6Lattice` is the **6D** list,
acting on the l-basis. The two are conjugate via the projection
matrix:

```text
g_phys = projection @ g_6d @ projection_inverse
```

Verified once at module load (analogous to how
`apeiron.symmetry` verifies |I| = 60 + closure + inverse closure
on import).

### 4.4 Polyhedron containment (existing primitive)

Apeiron already has exact-ZPhi containment via
`apeiron.track_b.realisation._point_in_convex_polyhedron` and
`_convex_polyhedra_interior_disjoint`. Both extend straightforwardly
to internal-space windows: rename the input type but keep the
ZPhi³ arithmetic.

### 4.5 The H₃ tetrahedron taxonomy (existing primitive)

`apeiron.track_b.taxonomy.build_h3_tetrahedra` returns the
9 (or 15, post-extension) similarity classes of H₃-compatible
tetrahedra. The §3.5 tile-type partitioning uses this for
classification of computed tiles into similarity classes.

---

## §5. Module layout

Proposed (post-design, pre-implementation):

```text
apeiron/track_d/                       # new track for cut-and-project
    __init__.py
    d6_lattice.py                      # D6Point, D6Lattice, projection
    windows.py                         # D6Window, ConvexPolyWindow, UnionWindow
    cps_search.py                      # tile_types_from_window, Route 1 driver
```

`track_d` is the natural name (after Track A: deformation, Track B:
substitution-search, Track C — never named — was the earlier
decoration pivot, now superseded). The track_b infrastructure that
gets reused (ZPhi arithmetic, polyhedron primitives) does not need
to move.

---

## §6. Route 1 driver pseudocode

```text
def test_route_1(
    lattice: D6Lattice = canonical_d6_lattice(),
    window_search: WindowSearch = enumerate_h3_compatible_windows(),
    radius: ZPhi = phi_squared,
) -> RouteOneResult:
    """Test the structural multi-prototile conjecture exhaustively.

    For each candidate window in `window_search`:
      1. Compute the tile types via tile_types_from_window.
      2. If exactly 1 type: COUNTEREXAMPLE (refutes conjecture).
      3. If ≥ 2 types: NoRealisation for the monotile question.

    `enumerate_h3_compatible_windows` is the load-bearing piece.
    Per the refined conjecture, the search is over windows whose
    prototile-shape structure is determined by the H₃ tetrahedral
    fundamental region. Concretely: windows whose decomposition
    into open subsets is consistent with the icosian model-set
    family per Moody §4.1.

    For ABCK specifically, the canonical window is a UnionWindow of
    three components (the three vertex-class windows). Variants
    explored:
      - Drop a vertex class (e.g., use only I and II).
      - Replace a class window with a different polytope of the
        same I_h-orbit.
      - Use a single connected window (the monotile candidate).

    Termination: finite combinatorial enumeration over the search
    space (which is bounded per Moody §4.1's framework
    exhaustion). The search space's exact cardinality is
    determined by the equivalence under I_h symmetry — likely
    ≤ thousands of windows after symmetry reduction.
    """
```

---

## §7. Open questions (for the post-P1 re-relay)

The following depend on Kramer-Papadopolos-Schlottmann-Zeidler
(1994) content; the design above is internally consistent but
some choices may need revisiting:

1. **Q14b(d) — Window as continuous parameter vs discrete choice?**
   The design above treats windows as an enumerable set
   (UnionWindow with ConvexPolyWindow components). If
   Kramer-Papadopolos treats windows as continuous parameters
   (e.g., parameterised by a vertex deformation), the
   `WindowSearch` interface would need to support continuous
   optimisation rather than enumeration.

2. **Q14b(e) — Window → prototile shape map.** The design assumes
   that prototile shapes emerge from H₃-fundamental-region
   combinatorics applied to a given window. If Kramer-Papadopolos
   states a more direct window → tile relation (e.g., "the tile is
   the projection of the window's pre-image under the canonical
   inflation"), §3.5 would need restructuring.

3. **Choice of canonical D₆ embedding.** The design above adopts
   the Al-Siyabi-Koca-Koca embedding (eq. 5 of arXiv:2003.13449).
   Kramer-Papadopolos may use a different embedding;
   compatibility with their notation is a check to do.

4. **Acceptance domain for the canonical ABCK model set.** Per
   Frettlöh §1.3, three explicit windows. Kramer-Papadopolos
   likely contains the canonical primary source for these windows
   and may define them differently.

If the P1 read clarifies any of (1)-(4), the corresponding §2-§6
sections will be revised before implementation.

---

## §8. What this document is and is not

**This document is** a pre-implementation design specification for
Route 1. It can be reviewed, revised, and committed before any
code lands. It standardises the data structures and operations so
that implementation is mechanical engineering rather than
research-during-coding.

**This document is not** an implementation, a binding
architectural decision, or a commitment that Route 1 will land.
Per the strict relay-then-act sequencing established in
strategy_pivot.md and the Q15-Q16 chain, the implementation
authorisation is a separate decision made after P1 (Kramer-
Papadopolos read) returns.

**Sequencing commitment carried forward.** No `apeiron/track_d/`
module is created before the post-P1 re-relay authorising
implementation. This document is the design artifact; engineering
follows separately.
