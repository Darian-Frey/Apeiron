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

import math
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

# Float approximation of phi, used solely for the scipy combinatorial-
# oracle call (see the scipy-oracle feedback memory). Never consumed
# back into the Z[phi] verification pipeline: scipy returns index
# triples, the original Vec3 values drive every downstream
# verification, and this float is discarded at the oracle boundary.
_PHI_FLOAT = (1.0 + math.sqrt(5.0)) / 2.0


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

    @classmethod
    def from_vertices(cls, vertices: Sequence[Vec3]) -> Polyhedron:
        """Build a canonical Polyhedron from a vertex set via hull + merge.

        Phase 1 — *combinatorial oracle*: ``scipy.spatial.ConvexHull``
        returns a triangulation of the hull as triangle index triples
        into ``vertices``. The oracle's float coordinates are consumed
        only to drive the oracle itself; they never cross into the
        Z[phi] pipeline (see the scipy-oracle feedback memory).

        Phase 2 — *exact validation in Z[phi]*. Three checks:

        - Every input vertex lies in the hull (``is_in_hull``). Catches
          the worst scipy failure mode — a truly-outside vertex dropped
          as interior.
        - Every scipy-classified-interior vertex is *strictly* interior
          (no face passes through it). Catches near-degenerate cases
          where scipy misclassifies a true hull vertex as interior.
        - V - E + F == 2 on the triangulation (``check_euler``).

        Phase 3 — *coplanar face merge* (CLAUDE.md §5.2). Co-oriented
        coplanar triangles are grouped by exact Z[phi] plane equality
        and their oriented boundary edges are assembled into the
        enclosing polygon cycle. Rhombic-triacontahedron-style faces
        (2 triangles → 1 rhombus) and pentagonal faces (3 triangles
        → 1 pentagon) are handled uniformly.

        Phase 4 — feed the merged face list through ``from_raw`` for
        canonicalisation.

        Raises ``ValueError`` if any phase-2 validation fires, or
        ``NotImplementedError`` if scipy's QhullError fires —  the
        exact-arithmetic hull fallback is a planned escape hatch but
        not yet implemented. For well-separated vertex sets
        (tetrahedra, cubes, the icosahedral fixtures) this path is
        not expected to be hit.
        """
        vs = list(vertices)
        if len(vs) < 4:
            raise ValueError(
                f"Polyhedron.from_vertices requires at least 4 vertices; "
                f"got {len(vs)}."
            )

        # Local imports: keep scipy off the module-level dependency path
        # so ``import apeiron.polyhedron`` works even in environments
        # that only need the predicates or the Polyhedron dataclass.
        import numpy as np
        from scipy.spatial import ConvexHull, QhullError

        coords = np.array([_vec3_to_float_triple(v) for v in vs], dtype=float)
        try:
            hull = ConvexHull(coords)
        except QhullError as e:
            raise NotImplementedError(
                "scipy.spatial.ConvexHull failed; the exact-arithmetic "
                "hull fallback is declared in CLAUDE.md §5.1 but not yet "
                "implemented. This vertex set needs the fallback."
            ) from e

        simplices_orig: list[tuple[int, int, int]] = [
            (int(s[0]), int(s[1]), int(s[2])) for s in hull.simplices
        ]
        hull_orig_indices = sorted({v for s in simplices_orig for v in s})
        index_map = {orig: new for new, orig in enumerate(hull_orig_indices)}
        hull_vs = [vs[i] for i in hull_orig_indices]
        hull_simplices = [
            tuple(index_map[i] for i in s) for s in simplices_orig
        ]

        # scipy does not return simplex vertex triples with a
        # consistent outward orientation. Reorient each in exact Z[phi]
        # before validation and coplanar-merge, both of which depend
        # on adjacent same-face simplices having matching normals.
        hull_simplices = _orient_simplices_outward(hull_vs, hull_simplices)

        _validate_scipy_output(
            all_vertices=vs,
            hull_orig_indices=set(hull_orig_indices),
            hull_vs=hull_vs,
            hull_simplices=hull_simplices,
        )

        merged_faces = _merge_coplanar_faces(hull_vs, hull_simplices)
        return cls.from_raw(hull_vs, merged_faces)


# -- hull-oracle helpers ------------------------------------------------


def _zphi_to_float(z: ZPhi) -> float:
    """Float embedding of z, used solely at the scipy oracle boundary."""
    return z.a + z.b * _PHI_FLOAT


def _vec3_to_float_triple(v: Vec3) -> tuple[float, float, float]:
    """Three-float representation of v, used solely for scipy."""
    return (_zphi_to_float(v.x), _zphi_to_float(v.y), _zphi_to_float(v.z))


def _triangle_plane(a: Vec3, b: Vec3, c: Vec3) -> tuple[Vec3, ZPhi]:
    """Return ``(n, d)`` for the plane through (a, b, c).

    ``n`` is the unnormalised cross product ``(b - a) × (c - a)`` (for
    outward-oriented triangles it points outward); ``d = n · a`` is the
    plane offset so any point ``p`` on the plane satisfies
    ``n.x * p.x + n.y * p.y + n.z * p.z == d``.
    """
    ba = b - a
    ca = c - a
    nx = ba.y * ca.z - ba.z * ca.y
    ny = ba.z * ca.x - ba.x * ca.z
    nz = ba.x * ca.y - ba.y * ca.x
    normal = Vec3(nx, ny, nz)
    offset = nx * a.x + ny * a.y + nz * a.z
    return normal, offset


def _orient_simplices_outward(
    hull_vs: Sequence[Vec3],
    simplices: Sequence[tuple[int, int, int]],
) -> list[tuple[int, int, int]]:
    """Reorient each simplex so its exact-Z[phi] normal points outward.

    ``scipy.spatial.ConvexHull`` does not guarantee a consistent outward
    vertex ordering across simplices — each simplex's outward equation
    lives in ``hull.equations`` as float metadata we deliberately don't
    consume (see the scipy-oracle feedback memory). Without reorienting,
    two coplanar simplices on the same face can have antiparallel
    ``_triangle_plane`` normals, which ``_same_oriented_plane`` correctly
    rejects — and the coplanar-face-merge fails to fire.

    Orientation is fixed in pure Z[phi] against the (unscaled) centroid
    direction. For a simplex with first vertex ``a`` and exact normal
    ``n``, we compare ``n`` to the vector ``sum(hull_vs) - n_count * a``,
    which is a positive-scalar multiple of ``centroid - a``. For an
    outward-pointing ``n`` this dot product is negative (the centroid
    lies on the inward side of every face of a convex hull). A positive
    dot product means ``n`` points inward, and we swap the last two
    simplex indices to flip the orientation.
    """
    n_count = len(hull_vs)
    sum_x = _ZPHI_ZERO
    sum_y = _ZPHI_ZERO
    sum_z = _ZPHI_ZERO
    for v in hull_vs:
        sum_x = sum_x + v.x
        sum_y = sum_y + v.y
        sum_z = sum_z + v.z
    oriented: list[tuple[int, int, int]] = []
    for (i, j, k) in simplices:
        a = hull_vs[i]
        normal, _ = _triangle_plane(hull_vs[i], hull_vs[j], hull_vs[k])
        dx = sum_x - n_count * a.x
        dy = sum_y - n_count * a.y
        dz = sum_z - n_count * a.z
        orient = normal.x * dx + normal.y * dy + normal.z * dz
        if orient == _ZPHI_ZERO:
            raise ValueError(
                f"Hull centroid is coplanar with simplex ({i}, {j}, {k}); "
                "the hull is degenerate (not full-dimensional)."
            )
        if orient > _ZPHI_ZERO:
            oriented.append((i, k, j))   # flip — normal was inward
        else:
            oriented.append((i, j, k))   # keep — normal is outward
    return oriented


def _same_oriented_plane(
    n1: Vec3, d1: ZPhi, n2: Vec3, d2: ZPhi,
) -> bool:
    """True iff ``(n1, d1)`` and ``(n2, d2)`` describe the same oriented
    plane — parallel normals pointing in the same direction, with the
    offsets consistently scaled.

    Tests performed, each in exact Z[phi]:

    1. ``n1 × n2 == 0`` (normals parallel, possibly antiparallel).
    2. ``n1 · n2 > 0`` (same direction; rejects antiparallel case).
    3. For the first non-zero component of ``n1``, the offsets match
       under the implied scalar: ``d1 * n2.c == d2 * n1.c``.
    """
    cross_x = n1.y * n2.z - n1.z * n2.y
    cross_y = n1.z * n2.x - n1.x * n2.z
    cross_z = n1.x * n2.y - n1.y * n2.x
    if cross_x != _ZPHI_ZERO or cross_y != _ZPHI_ZERO or cross_z != _ZPHI_ZERO:
        return False
    dot = n1.x * n2.x + n1.y * n2.y + n1.z * n2.z
    if not (dot > _ZPHI_ZERO):
        return False
    for c1, c2 in ((n1.x, n2.x), (n1.y, n2.y), (n1.z, n2.z)):
        if c1 != _ZPHI_ZERO:
            return d1 * c2 == d2 * c1
    # n1 is the zero vector — degenerate (collinear) triangle.
    return False


def _group_coplanar_triangles(
    hull_vs: Sequence[Vec3],
    triangles: Sequence[tuple[int, int, int]],
) -> list[list[int]]:
    """Partition triangle indices into equivalence classes of triangles
    sharing the same oriented supporting plane. O(n²) pairwise plane
    comparisons via union-find.
    """
    n = len(triangles)
    planes = [
        _triangle_plane(hull_vs[i], hull_vs[j], hull_vs[k])
        for (i, j, k) in triangles
    ]
    parent = list(range(n))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x: int, y: int) -> None:
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[rx] = ry

    for i in range(n):
        for j in range(i + 1, n):
            if _same_oriented_plane(*planes[i], *planes[j]):
                union(i, j)

    groups: dict[int, list[int]] = {}
    for i in range(n):
        groups.setdefault(find(i), []).append(i)
    return list(groups.values())


def _reconstruct_polygon_cycle(
    triangles_group: Sequence[tuple[int, int, int]],
) -> tuple[int, ...]:
    """Recover the ordered vertex cycle of the polygon tiled by a group
    of co-oriented coplanar triangles.

    For an outward-oriented triangulation of a convex polygon, every
    interior edge appears once as ``(u, v)`` and once as ``(v, u)``
    across adjacent triangles; those oriented edges cancel. The
    remaining oriented edges are the polygon boundary, each oriented
    in the polygon's CCW direction. Walking the boundary from its
    smallest-index vertex produces the canonical cycle.
    """
    oriented_edges: set[tuple[int, int]] = set()
    for (a, b, c) in triangles_group:
        for edge in ((a, b), (b, c), (c, a)):
            if edge in oriented_edges:
                raise ValueError(
                    f"Duplicate oriented edge {edge} in triangle group — "
                    "non-manifold coplanar configuration."
                )
            oriented_edges.add(edge)
    boundary: dict[int, int] = {
        u: v for (u, v) in oriented_edges if (v, u) not in oriented_edges
    }
    if not boundary:
        raise ValueError(
            "Coplanar triangle group has no boundary; the triangulation "
            "is fully closed or the triangles overlap incorrectly."
        )
    start = min(boundary)
    cycle = [start]
    current = start
    while True:
        if current not in boundary:
            raise ValueError(
                f"Boundary walk broke at vertex {current}; polygon is "
                "disconnected."
            )
        nxt = boundary[current]
        if nxt == start:
            break
        cycle.append(nxt)
        current = nxt
        if len(cycle) > len(boundary):
            raise ValueError(
                "Boundary walk failed to close within one revolution — "
                "polygon has a disconnected component."
            )
    if len(cycle) != len(boundary):
        raise ValueError(
            "Boundary cycle does not cover every boundary edge; polygon "
            "has multiple connected components."
        )
    return tuple(cycle)


def _merge_coplanar_faces(
    hull_vs: Sequence[Vec3],
    triangles: Sequence[tuple[int, int, int]],
) -> list[tuple[int, ...]]:
    """Merge co-oriented coplanar triangles into polygonal faces.

    Single-triangle groups pass through unchanged; multi-triangle groups
    are reassembled into their enclosing polygon cycle via
    ``_reconstruct_polygon_cycle``.
    """
    groups = _group_coplanar_triangles(hull_vs, triangles)
    merged: list[tuple[int, ...]] = []
    for group in groups:
        tris = [triangles[i] for i in group]
        if len(tris) == 1:
            merged.append(tris[0])
        else:
            merged.append(_reconstruct_polygon_cycle(tris))
    return merged


def _is_strictly_interior(
    point: Vec3,
    hull_vs: Sequence[Vec3],
    faces: Sequence[tuple[int, int, int]],
) -> bool:
    """True iff ``point`` has non-zero exact orientation against every
    face — i.e. lies in the open interior of the hull, with no face's
    supporting plane passing through it.
    """
    for (i, j, k) in faces:
        if exact_orientation(hull_vs[i], hull_vs[j], hull_vs[k], point) == _ZPHI_ZERO:
            return False
    return True


def _validate_scipy_output(
    all_vertices: Sequence[Vec3],
    hull_orig_indices: set[int],
    hull_vs: Sequence[Vec3],
    hull_simplices: Sequence[tuple[int, int, int]],
) -> None:
    """Exact Z[phi] validation of scipy's hull triangulation.

    Three checks, each producing a descriptive ``ValueError`` on
    failure. See ``Polyhedron.from_vertices`` for the full contract.
    """
    if not check_euler(len(hull_vs), hull_simplices):
        raise ValueError(
            f"scipy hull triangulation violates V - E + F = 2 "
            f"(V={len(hull_vs)}, F={len(hull_simplices)}); "
            "non-manifold output."
        )
    for orig_idx, v in enumerate(all_vertices):
        if not is_in_hull(v, hull_vs, hull_simplices):
            on_hull = "on hull" if orig_idx in hull_orig_indices else "interior"
            raise ValueError(
                f"scipy classified input vertex {orig_idx} as {on_hull}, "
                "but exact Z[phi] orientation places it strictly outside "
                "the hull."
            )
        if orig_idx in hull_orig_indices:
            continue
        # scipy declared this one interior; verify it is *strictly*
        # interior, not on a face boundary.
        if not _is_strictly_interior(v, hull_vs, hull_simplices):
            raise ValueError(
                f"scipy classified input vertex {orig_idx} as interior, "
                "but exact orientation places it on a face boundary — "
                "scipy missed a hull vertex."
            )
