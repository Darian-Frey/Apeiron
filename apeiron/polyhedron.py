"""Polyhedra with Z[phi]^3 vertices and exact validation predicates.

Per CLAUDE.md §5.1 and the 2026-04-20 hull-construction decision, the
construction of a Polyhedron from a vertex set proceeds in two phases:

1. A combinatorial oracle (typically ``scipy.spatial.ConvexHull``)
   returns a triangulation of the convex hull as triangle index
   triples into the input vertex list. The oracle's float coordinates
   are discarded at the boundary — only the discrete combinatorics
   cross into the verification pipeline. See the scipy-as-oracle
   feedback memory.
2. An exact validation pass in pure Z[phi] confirms the oracle's
   output via three predicates: orientation, interiority, and Euler's
   characteristic. If any predicate fires, the candidate falls back to
   an exact-arithmetic hull; that fallback path is a required escape
   hatch, not an optimisation target.

This module provides the three predicates. The ``Polyhedron`` dataclass,
its isometry action, and the hull-then-merge constructor follow in
subsequent commits.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

from apeiron.symmetry import Rotation, Vec3
from apeiron.zphi import ZPhi

__all__ = [
    "Polyhedron",
    "check_euler",
    "exact_orientation",
    "is_in_hull",
]


_ZPHI_ZERO = ZPhi(0, 0)


def exact_orientation(a: Vec3, b: Vec3, c: Vec3, d: Vec3) -> ZPhi:
    """Signed volume of the tetrahedron (a, b, c, d), as an exact ZPhi.

    Computes the scalar triple product ``(b - a) · ((c - a) × (d - a))``.
    The result is 6× the signed volume, but only the sign and the
    zero-test are consumed by orientation predicates, so the ×6 factor
    is left in place — never divided out, never rationalised.

    Sign interpretation (right-hand rule for the oriented triangle
    (a, b, c)):

    * ``> 0`` — d is on the positive side of the plane (a, b, c).
    * ``< 0`` — d is on the negative side.
    * ``== 0`` — d is coplanar with (a, b, c).

    The predicate is agnostic to any scaling convention on the input
    vectors (including the ×2 storage of CLAUDE.md §3.2): sign and
    zero are invariant under uniform scaling.
    """
    ba = b - a
    ca = c - a
    da = d - a
    cross_x = ca.y * da.z - ca.z * da.y
    cross_y = ca.z * da.x - ca.x * da.z
    cross_z = ca.x * da.y - ca.y * da.x
    return ba.x * cross_x + ba.y * cross_y + ba.z * cross_z


def is_in_hull(
    point: Vec3,
    hull_vertices: Sequence[Vec3],
    faces: Sequence[tuple[int, int, int]],
) -> bool:
    """Decide whether ``point`` lies in the convex hull.

    The hull is specified by ``hull_vertices`` and a triangulation
    ``faces``, where each face is a triple of indices into
    ``hull_vertices``. The function is agnostic to face winding order:
    for each face it derives the inward sign from a reference hull
    vertex not on that face — the first such vertex whose orientation
    against the face is non-zero.

    A point is in the hull iff, for every face, its orientation is
    either zero (on the face's supporting plane) or the same strict
    sign as the inward reference. Strict opposite sign on any face
    means the point is strictly outside and the function returns
    ``False``.

    Raises ``ValueError`` if the hull has fewer than four vertices, or
    if some face admits no off-plane hull vertex (indicating a
    degenerate, lower-dimensional hull).
    """
    hv = list(hull_vertices)
    n = len(hv)
    if n < 4:
        raise ValueError(
            f"Convex hull in R^3 requires at least 4 vertices; got {n}."
        )
    for (i, j, k) in faces:
        a, b, c = hv[i], hv[j], hv[k]
        ref_sign = _inward_reference_orientation(a, b, c, hv, (i, j, k))
        point_orient = exact_orientation(a, b, c, point)
        # Strict opposite sign means strictly outside this face.
        # Product of ZPhi is negative iff the signs strictly disagree.
        if point_orient * ref_sign < _ZPHI_ZERO:
            return False
    return True


def check_euler(
    num_vertices: int,
    faces: Sequence[tuple[int, int, int]],
) -> bool:
    """Verify the Euler characteristic ``V - E + F == 2`` for a convex
    polytope (topological sphere).

    The edge count is derived from the triangulation: each triangle
    contributes three undirected edges, deduplicated by sorted index
    pair. Faces are assumed to be index triples into a vertex list of
    size ``num_vertices``.
    """
    edges: set[tuple[int, int]] = set()
    for (i, j, k) in faces:
        for (u, v) in ((i, j), (j, k), (i, k)):
            edges.add((min(u, v), max(u, v)))
    return num_vertices - len(edges) + len(faces) == 2


def _inward_reference_orientation(
    a: Vec3,
    b: Vec3,
    c: Vec3,
    hull_vertices: Sequence[Vec3],
    skip: tuple[int, int, int],
) -> ZPhi:
    """Return the orientation of the first off-plane hull vertex
    against the face (a, b, c), excluding the face's own indices.

    Used by ``is_in_hull`` to pin the inward sign convention without
    trusting any particular face winding order. Raises ``ValueError``
    if every non-face hull vertex is coplanar with (a, b, c) —
    indicating a degenerate hull.
    """
    skip_set = set(skip)
    for r, vr in enumerate(hull_vertices):
        if r in skip_set:
            continue
        s = exact_orientation(a, b, c, vr)
        if s != _ZPHI_ZERO:
            return s
    raise ValueError(
        "No off-plane hull vertex for face {}; hull is degenerate "
        "(not full-dimensional).".format(skip)
    )


# -- canonical form helpers -------------------------------------------


_Vec3Key = tuple[int, int, int, int, int, int]


def _vec3_key(v: Vec3) -> _Vec3Key:
    """Total order on Vec3 by raw (a, b) coefficient tuples.

    This is a deterministic lexicographic key, not the real-value
    ordering — either works for canonicalisation, and raw coefficients
    are cheaper and avoid invoking the ZPhi sign algorithm.
    """
    return (v.x.a, v.x.b, v.y.a, v.y.b, v.z.a, v.z.b)


def _canonical_cycle(cycle: Sequence[int]) -> tuple[int, ...]:
    """Rotate ``cycle`` to start at its minimum index, preserving direction.

    A face cycle is a *directed* sequence — reversing it inverts the
    outward-normal orientation of the face. So canonicalisation only
    rotates, never reflects.
    """
    c = tuple(cycle)
    n = len(c)
    start = min(range(n), key=lambda i: c[i])
    return c[start:] + c[:start]


# -- Polyhedron -------------------------------------------------------


@dataclass(frozen=True, slots=True)
class Polyhedron:
    """A convex polyhedron with Z[phi]^3 vertices (×2 stored per §3.2).

    Invariants (enforced by ``__post_init__``):

    - ``vertices`` is a tuple of ``Vec3`` in canonical order — sorted
      ascending by ``_vec3_key``, pairwise distinct, with at least 4
      entries.
    - Each face is a tuple of vertex indices forming a closed ordered
      cycle (a directed polygon boundary). Triangular and polygonal
      faces are both representable; the merged-pentagonal-face case of
      the rhombic triacontahedron motivates polygonal support
      (CLAUDE.md §5.2).
    - Each face cycle is in canonical rotation (starts at its smallest
      index; direction preserved).
    - ``faces`` itself is sorted lexicographically.

    Construct via ``Polyhedron.from_raw`` unless you are deliberately
    providing already-canonical data.
    """

    vertices: tuple[Vec3, ...]
    faces: tuple[tuple[int, ...], ...]

    def __post_init__(self) -> None:
        if len(self.vertices) < 4:
            raise ValueError(
                f"Polyhedron requires at least 4 vertices; got {len(self.vertices)}."
            )
        if len(set(self.vertices)) != len(self.vertices):
            raise ValueError("Polyhedron vertices must be pairwise distinct.")
        keys = [_vec3_key(v) for v in self.vertices]
        if keys != sorted(keys):
            raise ValueError(
                "Polyhedron.vertices must be in canonical order; use "
                "Polyhedron.from_raw to build from an unsorted vertex list."
            )
        n = len(self.vertices)
        seen_faces: set[tuple[int, ...]] = set()
        for face in self.faces:
            if len(face) < 3:
                raise ValueError(f"Face has fewer than 3 vertices: {face}.")
            if len(set(face)) != len(face):
                raise ValueError(f"Face has repeated vertex indices: {face}.")
            for idx in face:
                if not 0 <= idx < n:
                    raise ValueError(
                        f"Face index {idx} out of range [0, {n}) in {face}."
                    )
            if _canonical_cycle(face) != face:
                raise ValueError(
                    f"Face {face} is not in canonical cycle order; use "
                    "Polyhedron.from_raw."
                )
            if face in seen_faces:
                raise ValueError(f"Duplicate face: {face}.")
            seen_faces.add(face)
        if tuple(self.faces) != tuple(sorted(self.faces)):
            raise ValueError(
                "Polyhedron.faces must be sorted; use Polyhedron.from_raw."
            )

    @classmethod
    def from_raw(
        cls,
        vertices: Sequence[Vec3],
        faces: Sequence[Sequence[int]],
    ) -> Polyhedron:
        """Construct a canonical Polyhedron from possibly-unordered input.

        Sorts the vertex list by ``_vec3_key``, remaps face indices
        through the sort permutation, rotates each face cycle to start
        at its smallest index (direction preserved), and sorts the face
        list lexicographically. The result satisfies every invariant
        checked by ``__post_init__``.

        Face winding is the caller's responsibility — this method does
        not reorient faces, because reversing a cycle inverts the
        outward normal and changes the polyhedron's "inside" half-space.
        """
        order = sorted(range(len(vertices)), key=lambda i: _vec3_key(vertices[i]))
        sorted_vertices = tuple(vertices[i] for i in order)
        remap = {old: new for new, old in enumerate(order)}
        remapped = [tuple(remap[i] for i in face) for face in faces]
        canonical = tuple(sorted(_canonical_cycle(f) for f in remapped))
        return cls(vertices=sorted_vertices, faces=canonical)

    def apply(self, g: Rotation) -> Polyhedron:
        """Return the polyhedron rotated by ``g``.

        Rotation generally permutes canonical vertex order, so the
        result is re-canonicalised via ``from_raw``. Face combinatorics
        are preserved (the rotation is an orientation-preserving
        isometry).
        """
        rotated = [g.apply(v) for v in self.vertices]
        return Polyhedron.from_raw(rotated, self.faces)
