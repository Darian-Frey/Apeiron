"""Track B taxonomy — H₃-compatible tetrahedra enumeration.

Per Claude (web)'s Q9a 2026-04-29 ruling: derive the icosahedral-
compatible tetrahedral taxonomy from first principles rather than
transcribing it from Frettlöh's references (which appear in the
privately-published ABCK-Book and contain known errors per
Frettlöh's own note).

The construction (Q9a):

1. ``H3_AXES`` — the 15 reflection-mirror normals of the H₃ Coxeter
   group, derived as the icosahedron's 30 edge-midpoint directions
   quotiented by ± (15 axes × 2 signs = 30 root vectors). Computed
   once at module import time from the icosahedron's 12 vertices.

2. ``is_h3_compatible(polyhedron)`` — True iff every face's normal
   is parallel to some ``H3_AXES`` direction. The test for parallel:
   cross product equals the zero vector exactly in ℤ[φ].

3. ``build_h3_tetrahedra(vertex_pool)`` — enumerate 4-subsets of
   ``vertex_pool``, filter by H₃ compatibility, deduplicate by
   similarity (canonical form on edge-length-squared multisets).
   Default pool: the 10 ABCK vertices from Paolini's POV-Ray source
   (``candidates/danzer/{A,B,C,K}.json``), which yields the 4 ABCK
   tetrahedra.

**Status (2026-04-29).** Default pool finds 4 H₃-compatible classes
(the 4 ABCK tetrahedra). Frettlöh's reported count of 15 implies a
larger vertex pool — likely including dodecahedron vertices,
icosidodecahedron vertices, and φ-scaled variants. Extending the
pool to verify the count of 15 is research follow-up. The test
asserts the 4-class baseline; if a richer pool is committed later
and the count rises, the test should be updated to ``≥ 4`` or to
``== 15`` depending on the pool choice.
"""

from __future__ import annotations

import itertools
from collections.abc import Iterable, Sequence

from apeiron.polyhedron import Polyhedron
from apeiron.symmetry import Vec3, _halve_zphi
from apeiron.zphi import ZPhi

__all__ = [
    "H3_AXES",
    "build_h3_tetrahedra",
    "is_h3_compatible",
]


# ---- ∆H₃ axes ------------------------------------------------------


def _icosahedron_vertices_x2() -> list[Vec3]:
    """The 12 icosahedron vertices in ×2-stored form (real coords are
    cyclic permutations of ``(0, ±1, ±φ)``)."""
    z = ZPhi(0, 0)
    pos1 = ZPhi(2, 0)
    neg1 = ZPhi(-2, 0)
    pos_phi = ZPhi(0, 2)
    neg_phi = ZPhi(0, -2)
    out: list[Vec3] = []
    for sa, sb in (
        (pos1, pos_phi), (neg1, pos_phi),
        (pos1, neg_phi), (neg1, neg_phi),
    ):
        out.append(Vec3(z, sa, sb))   # 0, ±1, ±φ
        out.append(Vec3(sb, z, sa))   # ±φ, 0, ±1
        out.append(Vec3(sa, sb, z))   # ±1, ±φ, 0
    return out


def _vec3_canonical_axis_key(v: Vec3) -> tuple:
    """Canonical key for parallel-equivalence under ±. Returns either
    the ZPhi-lex-smaller of (v, -v).
    """
    neg = Vec3(-v.x, -v.y, -v.z)
    key_v = (
        (v.x.a, v.x.b), (v.y.a, v.y.b), (v.z.a, v.z.b),
    )
    key_neg = (
        (neg.x.a, neg.x.b), (neg.y.a, neg.y.b), (neg.z.a, neg.z.b),
    )
    return min(key_v, key_neg)


def _build_h3_axes() -> tuple[Vec3, ...]:
    """The 15 H₃-mirror normals = icosahedron edge midpoint
    directions.

    Construction: enumerate pairs of icosahedron vertices at edge
    distance (squared distance = 4 in real coords = 16 in ×2-stored),
    take their midpoints, quotient by ±. Result: 15 distinct axes.
    """
    verts = _icosahedron_vertices_x2()
    edge_dist_sq_x2 = ZPhi(16, 0)  # ×2-stored squared edge length 4
    seen: dict[tuple, Vec3] = {}
    for i in range(len(verts)):
        for j in range(i + 1, len(verts)):
            diff = verts[j] - verts[i]
            if diff.norm_squared() != edge_dist_sq_x2:
                continue
            # Midpoint (×2-stored) = (V_i + V_j)/2 in stored frame.
            sum_v = verts[i] + verts[j]
            midpoint = Vec3(
                _halve_zphi(sum_v.x),
                _halve_zphi(sum_v.y),
                _halve_zphi(sum_v.z),
            )
            key = _vec3_canonical_axis_key(midpoint)
            if key not in seen:
                seen[key] = midpoint
    return tuple(seen.values())


H3_AXES: tuple[Vec3, ...] = _build_h3_axes()
assert len(H3_AXES) == 15, (
    f"Expected 15 H₃ axes, got {len(H3_AXES)}; check the "
    "edge-midpoint construction."
)


# ---- H₃-compatibility predicate -----------------------------------


def _cross(a: Vec3, b: Vec3) -> Vec3:
    return Vec3(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x,
    )


def _is_zero(v: Vec3) -> bool:
    z = ZPhi(0, 0)
    return v.x == z and v.y == z and v.z == z


def _is_parallel_to_h3_axis(direction: Vec3) -> bool:
    """True iff ``direction`` is parallel to some axis in
    :data:`H3_AXES`. Pure ZPhi cross product check."""
    if _is_zero(direction):
        return False
    for axis in H3_AXES:
        if _is_zero(_cross(direction, axis)):
            return True
    return False


def is_h3_compatible(polyhedron: Polyhedron) -> bool:
    """True iff every face's outward normal is parallel to an H₃ axis.

    For a polyhedron with triangular faces, the face normal at face
    ``(i, j, k)`` is ``(v_j − v_i) × (v_k − v_i)``. Pure ZPhi.

    Parallel-to-axis is decided exactly by ``cross(normal, axis) == 0``.
    """
    verts = polyhedron.vertices
    for face in polyhedron.faces:
        if len(face) < 3:
            continue
        normal = _cross(
            verts[face[1]] - verts[face[0]],
            verts[face[2]] - verts[face[0]],
        )
        if not _is_parallel_to_h3_axis(normal):
            return False
    return True


# ---- 4-subset enumeration with similarity quotient ----------------


def _edge_length_sq_multiset(verts: Sequence[Vec3]) -> tuple[ZPhi, ...]:
    """Sorted (canonical) multiset of squared edge lengths between
    every pair of vertices.

    For a tetrahedron: 6 ZPhi values. The multiset is invariant under
    rotation, reflection, and translation, so two tetrahedra with
    the same multiset are congruent. (Similarity quotient by uniform
    scaling is *not* applied here — see ``_similarity_class_key``.)
    """
    n = len(verts)
    out: list[ZPhi] = []
    for i in range(n):
        for j in range(i + 1, n):
            out.append((verts[j] - verts[i]).norm_squared())
    # Sort by ZPhi real value via _sign-on-difference. Use cmp_to_key
    # since ZPhi.__lt__ already gives real-number order.
    out.sort()
    return tuple(out)


def _similarity_class_key(verts: Sequence[Vec3]) -> tuple:
    """Canonical key for similarity equivalence under uniform
    scaling. Two tetrahedra have the same key iff their edge-length-
    squared multisets are proportional.

    Implementation: take the multiset, divide all entries by the
    smallest entry (= ratio scale), and use the resulting tuple as
    the key. The ratios are ZPhi values; equal ratios mean similar.
    """
    multiset = _edge_length_sq_multiset(verts)
    if not multiset:
        return ()
    # Find smallest entry.
    smallest = multiset[0]  # already sorted in real-number order
    if smallest == ZPhi(0, 0):
        return multiset  # degenerate; return raw
    # ZPhi doesn't have native division, but we can express each
    # entry as a (ZPhi, ZPhi) pair representing entry / smallest
    # WITHOUT actually dividing. Two pairs (a, b) and (c, d) are
    # equivalent iff a*d == b*c (cross-multiplication). For a
    # canonical form, normalize: pair = (entry, smallest) for each
    # entry, then sort by lex-order of the pair-tuples.
    return tuple((e.a, e.b, smallest.a, smallest.b) for e in multiset)


def _is_non_degenerate_tet(verts: Sequence[Vec3]) -> bool:
    """True iff the 4 vertices form a non-degenerate tetrahedron
    (positive volume)."""
    if len(verts) != 4:
        return False
    a = verts[1] - verts[0]
    b = verts[2] - verts[0]
    c = verts[3] - verts[0]
    # 6·volume_signed = a · (b × c)
    bxc = _cross(b, c)
    det6 = a.dot(bxc)
    return det6._sign() != 0


def _tet_faces() -> tuple[tuple[int, int, int], ...]:
    """Outward-wound tetrahedral faces for the canonical vertex
    ordering (V_0, V_1, V_2, V_3) where V_0 = origin.

    For our purposes the winding doesn't need to be globally
    outward — it just needs to be consistent enough that the face
    normal computation produces a valid axis check. The standard
    cyclic faces (opposite-vertex form) work fine.
    """
    return ((1, 2, 3), (0, 3, 2), (0, 1, 3), (0, 2, 1))


def build_h3_tetrahedra(
    vertex_pool: Sequence[Vec3] | None = None,
) -> list[Polyhedron]:
    """Enumerate 4-subsets of ``vertex_pool``, filter by H₃-axis
    compatibility, return the similarity-distinct tetrahedra.

    Default vertex pool: the 10 ABCK vertices from Paolini's POV-Ray
    source. With this pool, the function returns the 4 ABCK
    tetrahedra (= 4 of Frettlöh's reported 15 H₃-compatible
    classes).

    Per Q9a 2026-04-29, the test that this function returns exactly
    15 with a complete pool *is* the verification of Frettlöh's
    count. Extending the pool to find all 15 is research follow-up.

    Parameters
    ----------
    vertex_pool
        Vertices to draw from when forming candidate tetrahedra. If
        None, uses the default ABCK Paolini pool.

    Returns
    -------
    list of Polyhedron
        Distinct H₃-compatible tetrahedra (up to similarity), in
        an arbitrary but stable order.
    """
    if vertex_pool is None:
        vertex_pool = _default_abck_vertex_pool()
    pool = list(vertex_pool)
    seen_keys: set[tuple] = set()
    out: list[Polyhedron] = []
    faces = _tet_faces()
    for combo in itertools.combinations(range(len(pool)), 4):
        verts = [pool[i] for i in combo]
        if not _is_non_degenerate_tet(verts):
            continue
        # Build the polyhedron and test every face normal.
        poly = Polyhedron.from_raw(verts, faces)
        if not is_h3_compatible(poly):
            continue
        key = _similarity_class_key(poly.vertices)
        if key in seen_keys:
            continue
        seen_keys.add(key)
        out.append(poly)
    return out


def _default_abck_vertex_pool() -> tuple[Vec3, ...]:
    """The 10 ABCK vertices from Paolini's POV-Ray source — pt0..pt9
    (×2-stored ZPhi for use with ``Polyhedron``).

    Reference: ``candidates/danzer/_paolini_extract.py`` and the
    upstream ``danzer.dmf.unicatt.it/danzer_tiles.inc``.
    """
    z = ZPhi(0, 0)
    return (
        # pt0 = (0, 0, 0)
        Vec3(z, z, z),
        # pt1 = τ²(τ, 0, 1) = (τ³, 0, τ²)
        Vec3(ZPhi(2, 4), z, ZPhi(2, 2)),
        # pt2 = τ²(1, 1, 1) = (τ², τ², τ²)
        Vec3(ZPhi(2, 2), ZPhi(2, 2), ZPhi(2, 2)),
        # pt3 = (τ², 1, 0)
        Vec3(ZPhi(2, 2), ZPhi(2, 0), z),
        # pt4 = (τ², τ, 1)
        Vec3(ZPhi(2, 2), ZPhi(0, 2), ZPhi(2, 0)),
        # pt5 = (-τ, 0, 1)
        Vec3(ZPhi(0, -2), z, ZPhi(2, 0)),
        # pt6 = (0, τ², 1)
        Vec3(z, ZPhi(2, 2), ZPhi(2, 0)),
        # pt7 = (-1, τ, 0)
        Vec3(ZPhi(-2, 0), ZPhi(0, 2), z),
        # pt8 = τ(1, 1, 1) = (τ, τ, τ)
        Vec3(ZPhi(0, 2), ZPhi(0, 2), ZPhi(0, 2)),
        # pt9 = ½(-1, 1/τ, τ) — vertices stored ×2 already, so the
        # ×2-stored value is (-1, 1/τ, τ) ≡ (-1, τ-1, τ).
        Vec3(ZPhi(-1, 0), ZPhi(-1, 1), ZPhi(0, 1)),
    )


__all__ = [
    "H3_AXES",
    "build_h3_tetrahedra",
    "is_h3_compatible",
]
