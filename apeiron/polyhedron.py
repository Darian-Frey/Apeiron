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

from apeiron.symmetry import Vec3
from apeiron.zphi import ZPhi

__all__ = [
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
