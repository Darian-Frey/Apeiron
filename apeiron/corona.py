"""Corona BFS engine — the verification core (CLAUDE.md §5.1).

A *corona* around a central tile ``P`` is the Moore neighbourhood of
``P`` under a face-to-face tiling: every other tile that shares a face,
an edge, or a vertex with ``P``, positioned so no interior-overlap
occurs and the full solid angle around every edge and vertex of ``P``
is covered (zero angular defect).

The acceptance criteria for this module (CLAUDE.md §5.1 item 5) are:

- ``CoronaConfig`` dataclass: central tile plus a set of positioned
  neighbour tiles.
- ``corona_1(P)``: enumerate all complete first-shell configurations
  modulo I.
- ``corona_2(config)``: extend to second shell.
- Canonical hashing of configurations for BFS deduplication.
- Tests: cube has a unique ``corona_1``; rhombic dodecahedron has a
  unique ``corona_1``.

**Scope of this commit.** This is sub-commit A of a four-part
breakdown:

- **A (this commit)** — data model: ``PlacedTile``, ``CoronaConfig``,
  canonical ordering of neighbours, equality & hashing.
- **B** — geometric primitives: face-to-face placement enumeration,
  interior-overlap predicate, exact-ℤ[φ] ``angular_defect``.
- **C** — ``corona_1(P)`` BFS engine with constraint propagation;
  cube acceptance test.
- **D** — ``corona_2(config)`` and the rhombic-dodecahedron
  acceptance test.

Splitting follows the pattern used for ``polyhedron.py`` (commits A–E);
the rationale is that each sub-commit ships a reviewable unit of ~200
lines with its own tests, rather than a single ~800-line drop.

**Face-to-face vs complete Moore neighbourhood.** Internally the BFS
(sub-commit D) works face-to-face with constraint propagation —
placing a neighbour can force additional placements to fill the edge-
and vertex-wedges where that neighbour now participates, and any
placement creating an interior overlap or an unfillable wedge is
pruned. Externally ``corona_1`` returns the complete closed
neighbourhood: for a cube, all 26 neighbours (6 face + 12 edge + 8
vertex) in a 3×3×3 block minus the centre. This internal/external
distinction is the one thing most likely to confuse maintainers a
year out — flagging it here so it's captured at the module level.

**Two orthogonal predicates for corona completeness.** Per Claude
(web) 2026-04-22, corona completeness at the BFS level requires two
independent invariants, both checked separately:

1. ``incidence_defect(config, feature, expected) == 0`` for every
   edge and vertex of the central — combinatorial closure. Returns
   ``expected - actual_count``; zero means the expected number of
   tiles meet at the feature. For a face-transitive tile with
   dihedral angle dividing 2π, the expected edge-count is
   ``2π / dihedral`` and the expected vertex-count follows from the
   tile's solid angle.
2. ``has_interior_overlap(config) == False`` — geometric non-overlap.
   Returns ``True`` if any pair of placed tiles (including central
   with a neighbour, or two neighbours) has interior-overlapping
   volumes.

Both invariants are necessary: a configuration can satisfy (1) with
incidence counts correct but still overlap geometrically (e.g. two
tiles stacked on the same face), and it can satisfy (2) without
closure (too few tiles, leaving gaps). The BFS prunes on both.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass

from apeiron.polyhedron import Polyhedron, exact_orientation
from apeiron.symmetry import ICOSAHEDRAL, Mat3, Rotation, Vec3
from apeiron.zphi import ZPhi

__all__ = [
    "CoronaConfig",
    "Edge",
    "PlacedTile",
    "Vertex",
    "face_to_face_placements",
    "find_rotation",
    "has_interior_overlap",
    "incidence_defect",
]


_ZPHI_ZERO = ZPhi(0, 0)


@dataclass(frozen=True, slots=True)
class PlacedTile:
    """A placed copy of the monotile in a corona: translation + rotation.

    ``translation`` is in the ×2 storage form of CLAUDE.md §3.2.
    ``rotation`` is an element of the icosahedral group I, stored as
    ``2·g`` under the denominator-2 convention of CLAUDE.md §3.3.

    Semantically this is the corona analogue of
    ``apeiron.substitution.PositionedTile`` — both carry a translation
    and a rotation — but the types intentionally differ: ``PlacedTile``
    belongs to a corona and will eventually grow adjacency metadata
    (which face / edges / vertices of the central it participates in;
    which wedge of its own boundary is filled by which other neighbour)
    that does not belong on ``PositionedTile``. Per Claude (web)'s
    2026-04-20 guidance the fields are duplicated rather than factored
    through a shared ``Placement`` type; if a third corona-like
    consumer arises later, refactor then.

    Equality and hash are induced by the frozen dataclass.
    """

    translation: Vec3
    rotation: Rotation


def _placed_tile_key(pt: PlacedTile) -> tuple[int, ...]:
    """Total order on ``PlacedTile`` for deterministic neighbour ordering.

    Sorts first by translation (lex on raw ZPhi coefficient pairs
    component-by-component), then by rotation matrix entries in row-
    major order under the same raw-coefficient lex convention. Not the
    real-value ordering — this is a deterministic surrogate, which is
    all ``CoronaConfig``'s canonical-order invariant requires.
    """
    t = pt.translation
    m = pt.rotation.matrix
    return (
        t.x.a, t.x.b, t.y.a, t.y.b, t.z.a, t.z.b,
        m.row0.x.a, m.row0.x.b, m.row0.y.a, m.row0.y.b, m.row0.z.a, m.row0.z.b,
        m.row1.x.a, m.row1.x.b, m.row1.y.a, m.row1.y.b, m.row1.z.a, m.row1.z.b,
        m.row2.x.a, m.row2.x.b, m.row2.y.a, m.row2.y.b, m.row2.z.a, m.row2.z.b,
    )


@dataclass(frozen=True, slots=True)
class CoronaConfig:
    """A corona configuration: central tile and its placed neighbours.

    Fields
    ------
    central : Polyhedron
        The tile being coronated. All neighbours are placed copies of
        ``central`` under icosahedral rotations and ×2-stored
        translations. (The monotile case is the research target;
        multi-prototile support would require extending ``PlacedTile``
        with a prototile index, which has not been done yet.)
    neighbours : tuple of PlacedTile
        Placed copies of the central surrounding it, in canonical
        order (sorted by ``_placed_tile_key``). Pairwise distinct.

    Invariants (enforced by ``__post_init__``)
    -----------------------------------------
    - ``neighbours`` is sorted ascending by ``_placed_tile_key``.
    - All neighbours are pairwise distinct.

    Construction
    ------------
    Prefer ``CoronaConfig.from_neighbours``, which sorts and
    deduplicates the input iterable. Direct construction is
    permissible when the caller guarantees canonical form (e.g., when
    producing a ``CoronaConfig`` from another, already-canonical one
    via a pure-rotation action).

    Equality and hashing
    --------------------
    Dataclass equality is pointwise: ``c1 == c2`` iff they share the
    same central and the same neighbour tuple. This is *not*
    equivalence up to the action of I on the whole configuration —
    that is a separate orbit-level notion that will be added when the
    BFS engine (sub-commit C) needs it for deduplication.
    """

    central: Polyhedron
    neighbours: tuple[PlacedTile, ...]

    def __post_init__(self) -> None:
        keys = [_placed_tile_key(n) for n in self.neighbours]
        if keys != sorted(keys):
            raise ValueError(
                "CoronaConfig.neighbours must be in canonical order; "
                "use CoronaConfig.from_neighbours to construct from "
                "an unordered iterable."
            )
        if len(set(self.neighbours)) != len(self.neighbours):
            raise ValueError("CoronaConfig has duplicate neighbours.")

    @classmethod
    def from_neighbours(
        cls,
        central: Polyhedron,
        neighbours: Iterable[PlacedTile],
    ) -> CoronaConfig:
        """Construct a canonical ``CoronaConfig`` from any iterable.

        Deduplicates by value, sorts by ``_placed_tile_key``. Duplicate
        neighbours in the input are silently dropped — BFS will
        routinely produce the same placement via multiple search
        paths, and raising on duplicates would push dedup responsibility
        to every caller.
        """
        ordered = sorted(set(neighbours), key=_placed_tile_key)
        return cls(central=central, neighbours=tuple(ordered))


# -- geometric primitives (sub-commit B) ------------------------------


def find_rotation(target: Mat3) -> Rotation | None:
    """Return the ``Rotation`` in ``ICOSAHEDRAL`` whose matrix equals
    ``target``, or ``None`` if the target matrix is not in ``I``.

    ``target`` is a ``Mat3`` of ZPhi entries in the ×2 storage
    convention of CLAUDE.md §3.3 — i.e. the numerator of ``2·g`` for
    the sought rotation ``g``. This helper exists so that the eventual
    optimisation (a precomputed index over orbit representatives of
    the tile's faces under I, per the 2026-04-20 Claude (web) relay)
    is a one-function swap with no call-site changes. The current
    implementation is an O(60) linear scan, which is more than fast
    enough for BFS-era usage and wrong to optimise before the
    pipeline runs end-to-end.
    """
    for g in ICOSAHEDRAL:
        if g.matrix == target:
            return g
    return None


def _reversed_face_cycle(cycle: tuple[Vec3, ...]) -> tuple[Vec3, ...]:
    """Return the cycle seen from the opposite side.

    For a face cycle ``(v_0, v_1, ..., v_{n-1})`` traversed CCW from
    outside the owning polyhedron, the face-to-face neighbour's
    matching face has the same vertex positions but traversed CCW
    from outside the *neighbour* — which is the reversal
    ``(v_0, v_{n-1}, v_{n-2}, ..., v_1)``. Starting vertex preserved
    so callers can canonicalise directly.
    """
    if len(cycle) < 2:
        return cycle
    return (cycle[0],) + tuple(reversed(cycle[1:]))


def _find_face_isometry(
    source_pts: tuple[Vec3, ...],
    target_pts: tuple[Vec3, ...],
) -> tuple[Rotation, Vec3] | None:
    """Find ``(g, t)`` with ``g ∈ ICOSAHEDRAL`` such that
    ``g.apply(source_pts[i]) + t == target_pts[i]`` for every ``i``,
    or ``None`` if no such isometry exists.

    Requires ``len(source_pts) == len(target_pts) >= 3`` and the first
    three source points non-collinear. For each ``g ∈ I`` the
    algorithm checks that ``g`` rotates the first two edge-vectors
    ``a_1 - a_0`` and ``a_2 - a_0`` to the matching target edges; if
    both hold, ``t`` is solved from ``a_0 ↦ b_0`` and the remaining
    correspondences (for ``n > 3`` faces) are verified.

    Inputs and outputs are in the ×2 storage convention of CLAUDE.md
    §3.2: differences of ×2-stored ``Vec3`` values are themselves
    ×2-stored, so ``Rotation.apply`` consumes and produces the right
    thing with no conversion.
    """
    n = len(source_pts)
    if n < 3 or n != len(target_pts):
        return None
    a0, a1, a2 = source_pts[0], source_pts[1], source_pts[2]
    b0, b1, b2 = target_pts[0], target_pts[1], target_pts[2]
    u0 = a1 - a0
    u1 = a2 - a0
    v0 = b1 - b0
    v1 = b2 - b0
    for g in ICOSAHEDRAL:
        if g.apply(u0) != v0:
            continue
        if g.apply(u1) != v1:
            continue
        t = b0 - g.apply(a0)
        # Verify remaining correspondences (triangular faces have
        # none; polygonal faces have one check per vertex beyond the
        # first three).
        ok = True
        for a_i, b_i in zip(source_pts[3:], target_pts[3:], strict=True):
            if g.apply(a_i) + t != b_i:
                ok = False
                break
        if ok:
            return g, t
    return None


def face_to_face_placements(
    central: Polyhedron,
    central_face_index: int,
) -> list[PlacedTile]:
    """Enumerate all face-to-face placements of a copy of ``central``
    against ``central.faces[central_face_index]``.

    A face-to-face placement ``(t, g)`` is one where some face of
    ``g·central + t`` coincides with the central face at
    ``central_face_index`` *with opposite orientation* — the two
    polyhedra are on opposite sides of their shared face plane, which
    is the standard face-to-face tiling condition.

    Algorithm
    ---------
    Let the central face be ``f_c`` with vertex cycle
    ``(v_0, v_1, ..., v_{n-1})`` (CCW from outside central). The
    face-to-face target cycle is the reversal
    ``(v_0, v_{n-1}, ..., v_1)``. For each face ``f_j`` of ``central``
    with the same vertex count and each cyclic rotation ``s`` of
    ``f_j``'s cycle, compute the unique rigid isometry mapping the
    rotated ``f_j`` cycle to the reversed target, via
    ``_find_face_isometry``. Keep only those isometries whose rotation
    lies in ``ICOSAHEDRAL`` (all others are rejected by that helper).
    Deduplicate across cyclic rotations.

    Complexity for a polyhedron with ``F`` faces of size ``n`` is
    ``O(F · n · |I|)`` = ``O(60 F n)``. For the rhombic
    triacontahedron (``F = 30``, ``n = 4``) that's ``7200`` rotation
    tries per central face — sub-millisecond in Python.

    Returns a list of distinct ``PlacedTile`` in canonical order (by
    ``_placed_tile_key``). The caller decides which placements to
    admit into a corona — this function only determines which
    ``(t, g)`` are geometrically valid as face-to-face with the given
    central face.
    """
    if not 0 <= central_face_index < len(central.faces):
        raise IndexError(
            f"central_face_index {central_face_index} out of range "
            f"[0, {len(central.faces)})."
        )
    c_face = central.faces[central_face_index]
    c_pts = tuple(central.vertices[i] for i in c_face)
    target_cycle = _reversed_face_cycle(c_pts)
    n = len(target_cycle)

    placements: set[PlacedTile] = set()
    for face in central.faces:
        if len(face) != n:
            continue
        source_base = tuple(central.vertices[i] for i in face)
        for shift in range(n):
            rotated_source = tuple(source_base[(shift + i) % n] for i in range(n))
            isometry = _find_face_isometry(rotated_source, target_cycle)
            if isometry is None:
                continue
            g, t = isometry
            placements.add(PlacedTile(translation=t, rotation=g))
    return sorted(placements, key=_placed_tile_key)


# -- feature types and completeness predicates (sub-commit C) ---------


@dataclass(frozen=True, slots=True)
class Vertex:
    """A vertex of a central polyhedron, referenced by its index.

    The index is into ``CoronaConfig.central.vertices``. Vertex
    features are one of the two feature types (the other being
    ``Edge``) at which corona completeness is checked via
    ``incidence_defect``.
    """

    index: int

    def __post_init__(self) -> None:
        if self.index < 0:
            raise ValueError(f"Vertex index must be non-negative; got {self.index}.")


@dataclass(frozen=True, slots=True)
class Edge:
    """An edge of a central polyhedron, as a sorted vertex-index pair.

    The invariant ``lo < hi`` canonicalises the representation: the
    two vertices of an undirected edge always appear in ascending
    index order, so equality and hashing don't depend on the order
    the caller specifies.
    """

    lo: int
    hi: int

    def __post_init__(self) -> None:
        if self.lo < 0 or self.hi < 0:
            raise ValueError(
                f"Edge indices must be non-negative; got ({self.lo}, {self.hi})."
            )
        if self.lo >= self.hi:
            raise ValueError(
                f"Edge requires lo < hi; got ({self.lo}, {self.hi})."
            )


def _edges_of(P: Polyhedron) -> frozenset[tuple[int, int]]:
    """Undirected edges of ``P`` as sorted ``(lo, hi)`` index pairs."""
    edges: set[tuple[int, int]] = set()
    for face in P.faces:
        n = len(face)
        for idx in range(n):
            u, v = face[idx], face[(idx + 1) % n]
            edges.add((min(u, v), max(u, v)))
    return frozenset(edges)


def _placed_vertices(
    P: Polyhedron, placement: PlacedTile,
) -> list[Vec3]:
    """Vertices of ``P`` transformed by ``placement``.

    Each local vertex ``v`` of ``P`` becomes
    ``placement.rotation.apply(v) + placement.translation``, in the
    ×2 storage form shared across the pipeline (CLAUDE.md §3.2).
    """
    return [
        placement.rotation.apply(v) + placement.translation
        for v in P.vertices
    ]


def incidence_defect(
    config: CoronaConfig,
    feature: Vertex | Edge,
    *,
    expected: int,
) -> int:
    """Combinatorial closure defect at a feature of the central.

    Returns ``expected - actual_count``, where ``actual_count`` is the
    number of placed tiles in ``config`` (including the central
    itself) that include this feature in their placed boundary —
    a vertex coincidence for a ``Vertex`` feature, a vertex-pair
    coincidence *and* an edge relation in the placed tile's local
    frame for an ``Edge`` feature.

    A **complete** corona (per Claude (web) 2026-04-22) requires
    ``incidence_defect(config, f, expected=...) == 0`` at every edge
    and vertex of the central *together with*
    ``not has_interior_overlap(config)``. The two predicates are
    orthogonal and both necessary — see the module docstring.

    The ``expected`` count is tiling-specific: it depends on how the
    tile is intended to be used (dihedral angle divides 2π → expected
    edge count is the quotient), not only on the tile itself. For
    ergonomics the count is passed explicitly rather than derived
    from local geometry, since exact-angle reasoning is not in
    general reducible to integer arithmetic in a form that is both
    clean and correct.

    Raises ``ValueError`` if the feature doesn't belong to the
    central polyhedron (out-of-range vertex index, or an edge that
    isn't in the central's edge set).
    """
    central = config.central
    n_vertices = len(central.vertices)
    if isinstance(feature, Vertex):
        if not 0 <= feature.index < n_vertices:
            raise ValueError(
                f"Vertex index {feature.index} is outside "
                f"[0, {n_vertices})."
            )
        actual = _count_vertex_incidences(config, feature.index)
        return expected - actual
    if isinstance(feature, Edge):
        if feature.hi >= n_vertices:
            raise ValueError(
                f"Edge index {feature.hi} is outside [0, {n_vertices})."
            )
        if (feature.lo, feature.hi) not in _edges_of(central):
            raise ValueError(
                f"Edge ({feature.lo}, {feature.hi}) is not an edge of "
                "the central polyhedron."
            )
        actual = _count_edge_incidences(config, feature.lo, feature.hi)
        return expected - actual
    raise TypeError(
        f"feature must be Vertex or Edge; got {type(feature).__name__}."
    )


def _count_vertex_incidences(
    config: CoronaConfig, vertex_index: int,
) -> int:
    """Number of placed tiles whose placed vertex set includes the
    central's vertex at ``vertex_index``. The central counts as 1.
    """
    P = config.central
    target = P.vertices[vertex_index]
    count = 1   # the central itself
    for placement in config.neighbours:
        if target in _placed_vertices(P, placement):
            count += 1
    return count


def _count_edge_incidences(
    config: CoronaConfig, lo: int, hi: int,
) -> int:
    """Number of placed tiles whose placed boundary includes the
    central's edge ``(lo, hi)``. The central counts as 1.

    A placement is counted iff both central-vertex positions appear
    among its placed vertices *and* the corresponding local-frame
    indices form an edge of the prototype — merely having the
    positions coincide is not enough, since two non-adjacent vertices
    of the prototype can happen to land at the edge endpoints after
    rotation + translation (e.g. across a face diagonal).
    """
    P = config.central
    target_lo = P.vertices[lo]
    target_hi = P.vertices[hi]
    local_edges = _edges_of(P)
    count = 1   # the central itself (caller already validated the edge)
    for placement in config.neighbours:
        placed = _placed_vertices(P, placement)
        try:
            i = placed.index(target_lo)
            j = placed.index(target_hi)
        except ValueError:
            continue
        if (min(i, j), max(i, j)) in local_edges:
            count += 1
    return count


def has_interior_overlap(config: CoronaConfig) -> bool:
    """Return ``True`` iff any pair of placed tiles in ``config`` has
    interior-overlapping volumes.

    Uses face-plane separating-axis reasoning: for each pair of
    placements (central with a neighbour, or two neighbours), the
    pair does *not* interior-overlap iff some face of one polytope
    strictly separates the other — every vertex of the other lies on
    the non-interior side of that face's supporting plane or exactly
    on the plane. The check is run for every face of each polytope;
    if no separator is found, we conservatively report overlap.

    **Completeness.** Face-plane SAT is sufficient for corona
    configurations where non-face-sharing neighbours meet at edges or
    vertices — for the cube, rhombic dodecahedron, and RTH fixtures,
    the face planes containing the shared feature always separate.
    Pathological pairs touching only via skew edges would require a
    full edge-edge SAT sweep; promote to that when a failing
    configuration ever surfaces, not before.

    **Exactness.** Every comparison goes through ``exact_orientation``
    (polyhedron.py) in ℤ[φ]. No floats anywhere.
    """
    central = config.central
    central_placement = PlacedTile(
        translation=Vec3(_ZPHI_ZERO, _ZPHI_ZERO, _ZPHI_ZERO),
        rotation=Rotation.identity(),
    )
    all_placements: list[PlacedTile] = [central_placement, *config.neighbours]
    # Precompute placed vertices once per placement.
    placed_vertex_lists: list[list[Vec3]] = [
        _placed_vertices(central, p) for p in all_placements
    ]

    n = len(all_placements)
    for i in range(n):
        for j in range(i + 1, n):
            if _placements_interior_overlap(
                central,
                placed_vertex_lists[i],
                placed_vertex_lists[j],
            ):
                return True
    return False


def _placements_interior_overlap(
    P: Polyhedron,
    verts_a: Sequence[Vec3],
    verts_b: Sequence[Vec3],
) -> bool:
    """Do two placed copies of ``P`` have interior-overlapping volumes?

    Tries each face of ``P`` as a separating plane under both
    placements; if any separates, the pair does not interior-overlap.
    Conservative return if no separator is found via face planes
    alone.
    """
    if _any_face_separates(P, verts_a, verts_b):
        return False
    if _any_face_separates(P, verts_b, verts_a):
        return False
    return True


def _any_face_separates(
    P: Polyhedron,
    verts_own: Sequence[Vec3],
    verts_other: Sequence[Vec3],
) -> bool:
    """``True`` iff some face of ``P`` (placed at ``verts_own``)
    strictly separates ``verts_other`` from ``verts_own``'s interior.

    "Strictly separates" means every ``verts_other`` vertex lies on
    the non-interior side of the face's supporting plane or exactly
    on the plane itself. Equality of signs is decided by
    ``exact_orientation``.
    """
    for face_indices in P.faces:
        if len(face_indices) < 3:
            continue
        a = verts_own[face_indices[0]]
        b = verts_own[face_indices[1]]
        c = verts_own[face_indices[2]]
        interior_sign = _face_interior_sign(
            verts_own, face_indices, a, b, c,
        )
        if interior_sign == 0:
            continue   # degenerate; try another face
        strictly_separates = True
        for v in verts_other:
            orient = exact_orientation(a, b, c, v)
            if _sign_of(orient) == interior_sign:
                strictly_separates = False
                break
        if strictly_separates:
            return True
    return False


def _face_interior_sign(
    verts_own: Sequence[Vec3],
    face_indices: tuple[int, ...],
    a: Vec3,
    b: Vec3,
    c: Vec3,
) -> int:
    """Sign of a non-face vertex's orientation against ``(a, b, c)``.

    For a convex polytope, all vertices not on the given face lie on
    the same side; the first non-zero orientation pins the interior
    sign. Returns ``0`` if every non-face vertex is coplanar with the
    face (degenerate polytope; caller handles).
    """
    face_set = set(face_indices)
    for idx, v in enumerate(verts_own):
        if idx in face_set:
            continue
        s = _sign_of(exact_orientation(a, b, c, v))
        if s != 0:
            return s
    return 0


def _sign_of(z: ZPhi) -> int:
    """Integer sign of a ZPhi via comparison to zero. +1 / 0 / -1."""
    if z > _ZPHI_ZERO:
        return 1
    if z < _ZPHI_ZERO:
        return -1
    return 0
