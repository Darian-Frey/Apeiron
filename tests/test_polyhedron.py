"""Tests for apeiron.polyhedron.

Scope for this commit (A of the polyhedron implementation): the three
validation predicates used to check a convex-hull oracle's output in
exact Z[phi].

- ``exact_orientation``: signed volume / scalar triple product.
- ``is_in_hull``: point-in-polytope decision via per-face sign checks.
- ``check_euler``: V − E + F == 2 for a topological sphere.

Downstream commits add the Polyhedron dataclass, the hull-then-merge
constructor, and the RTH integration test (the canonical icosahedral
fixture) per CLAUDE.md §5.1 acceptance criteria.
"""

from __future__ import annotations

import pytest

from apeiron.polyhedron import (
    Polyhedron,
    _merge_coplanar_faces,
    _orient_simplices_outward,
    _reconstruct_polygon_cycle,
    _same_oriented_plane,
    _triangle_plane,
    check_euler,
    exact_orientation,
    is_in_hull,
)
from apeiron.symmetry import ICOSAHEDRAL, ROT_2, ROT_3, ROT_5, Rotation, Vec3
from apeiron.zphi import ZPhi


def _v(x: tuple[int, int], y: tuple[int, int], z: tuple[int, int]) -> Vec3:
    """Compact Vec3 constructor from three (a, b) ZPhi coefficient pairs."""
    return Vec3(ZPhi(*x), ZPhi(*y), ZPhi(*z))


# -- exact_orientation -------------------------------------------------


class TestExactOrientation:
    def test_standard_tetrahedron_positive(self) -> None:
        # Vertices of a right-handed tetrahedron scaled to stored-×2 coords.
        # (c - a) × (d - a) = (0, 2, 0) × (0, 0, 2) = (4, 0, 0).
        # (b - a) . (4, 0, 0) = 2 * 4 = 8.
        a = _v((0, 0), (0, 0), (0, 0))
        b = _v((2, 0), (0, 0), (0, 0))
        c = _v((0, 0), (2, 0), (0, 0))
        d = _v((0, 0), (0, 0), (2, 0))
        assert exact_orientation(a, b, c, d) == ZPhi(8, 0)

    def test_coplanar_is_zero(self) -> None:
        a = _v((0, 0), (0, 0), (0, 0))
        b = _v((2, 0), (0, 0), (0, 0))
        c = _v((0, 0), (2, 0), (0, 0))
        d = _v((2, 0), (2, 0), (0, 0))   # still on the plane z = 0
        assert exact_orientation(a, b, c, d) == ZPhi(0, 0)

    def test_opposite_side_is_negative(self) -> None:
        a = _v((0, 0), (0, 0), (0, 0))
        b = _v((2, 0), (0, 0), (0, 0))
        c = _v((0, 0), (2, 0), (0, 0))
        d = _v((0, 0), (0, 0), (-2, 0))   # z < 0
        assert exact_orientation(a, b, c, d) == ZPhi(-8, 0)

    def test_swap_reverses_sign(self) -> None:
        # Swapping any two vertices of the oriented triangle flips sign.
        a = _v((0, 0), (0, 0), (0, 0))
        b = _v((2, 0), (0, 0), (0, 0))
        c = _v((0, 0), (2, 0), (0, 0))
        d = _v((0, 0), (0, 0), (2, 0))
        base = exact_orientation(a, b, c, d)
        assert exact_orientation(b, a, c, d) == -base
        assert exact_orientation(a, c, b, d) == -base

    def test_scale_invariance_of_sign(self) -> None:
        # Uniform scaling by a positive integer preserves sign.
        a = _v((0, 0), (0, 0), (0, 0))
        b = _v((2, 0), (0, 0), (0, 0))
        c = _v((0, 0), (2, 0), (0, 0))
        d = _v((0, 0), (0, 0), (2, 0))
        a2 = _v((0, 0), (0, 0), (0, 0))
        b2 = _v((6, 0), (0, 0), (0, 0))
        c2 = _v((0, 0), (6, 0), (0, 0))
        d2 = _v((0, 0), (0, 0), (6, 0))
        orig = exact_orientation(a, b, c, d)
        scaled = exact_orientation(a2, b2, c2, d2)
        assert orig > ZPhi(0, 0)
        assert scaled > ZPhi(0, 0)

    def test_phi_cubed_volume(self) -> None:
        # Tetrahedron with b, c, d along the axes at distance phi.
        # volume = phi * phi * phi = phi^3 = 2*phi + 1 = ZPhi(1, 2).
        a = _v((0, 0), (0, 0), (0, 0))
        b = _v((0, 1), (0, 0), (0, 0))
        c = _v((0, 0), (0, 1), (0, 0))
        d = _v((0, 0), (0, 0), (0, 1))
        assert exact_orientation(a, b, c, d) == ZPhi(1, 2)

    def test_matches_icosahedron_vertex_triangle(self) -> None:
        # Three adjacent icosahedron vertices stored in ×2 form:
        #   V1 = (0, 1, phi)   stored (0, 2, 2*phi)
        #   V2 = (0, -1, phi)  stored (0, -2, 2*phi)
        #   V3 = (1, phi, 0)   stored (2, 2*phi, 0)
        # The apex at the origin gives a non-zero signed volume; we
        # only require a definite sign here.
        v0 = _v((0, 0), (0, 0), (0, 0))
        v1 = _v((0, 0), (2, 0), (0, 2))
        v2 = _v((0, 0), (-2, 0), (0, 2))
        v3 = _v((2, 0), (0, 2), (0, 0))
        vol = exact_orientation(v0, v1, v2, v3)
        assert vol != ZPhi(0, 0)


# -- is_in_hull --------------------------------------------------------


def _unit_cube() -> tuple[list[Vec3], list[tuple[int, int, int]]]:
    """Stored-×2 unit cube with vertices at (0, 2)×(0, 2)×(0, 2).

    Indexing: (i, j, k) ↦ 4*i + 2*j + k. Triangulation: 2 triangles
    per face, 12 triangles total. Winding is intentionally mixed so
    the ``is_in_hull`` sign-agnostic behaviour is exercised.
    """
    vs: list[Vec3] = []
    for i in (0, 1):
        for j in (0, 1):
            for k in (0, 1):
                vs.append(_v((2 * i, 0), (2 * j, 0), (2 * k, 0)))

    def idx(i: int, j: int, k: int) -> int:
        return 4 * i + 2 * j + k

    faces = [
        # z = 0
        (idx(0, 0, 0), idx(1, 0, 0), idx(0, 1, 0)),
        (idx(1, 0, 0), idx(1, 1, 0), idx(0, 1, 0)),
        # z = 1
        (idx(0, 0, 1), idx(0, 1, 1), idx(1, 0, 1)),
        (idx(1, 0, 1), idx(0, 1, 1), idx(1, 1, 1)),
        # y = 0
        (idx(0, 0, 0), idx(0, 0, 1), idx(1, 0, 0)),
        (idx(1, 0, 0), idx(0, 0, 1), idx(1, 0, 1)),
        # y = 1
        (idx(0, 1, 0), idx(1, 1, 0), idx(0, 1, 1)),
        (idx(1, 1, 0), idx(1, 1, 1), idx(0, 1, 1)),
        # x = 0
        (idx(0, 0, 0), idx(0, 1, 0), idx(0, 0, 1)),
        (idx(0, 1, 0), idx(0, 1, 1), idx(0, 0, 1)),
        # x = 1
        (idx(1, 0, 0), idx(1, 0, 1), idx(1, 1, 0)),
        (idx(1, 1, 0), idx(1, 0, 1), idx(1, 1, 1)),
    ]
    return vs, faces


class TestIsInHull:
    def test_every_vertex_is_in_hull(self) -> None:
        vs, faces = _unit_cube()
        for v in vs:
            assert is_in_hull(v, vs, faces)

    def test_center_is_in_hull(self) -> None:
        # Real (0.5, 0.5, 0.5) → stored (1, 1, 1).
        vs, faces = _unit_cube()
        center = _v((1, 0), (1, 0), (1, 0))
        assert is_in_hull(center, vs, faces)

    def test_face_centre_is_in_hull(self) -> None:
        # Real (0.5, 0.5, 0) → stored (1, 1, 0).
        vs, faces = _unit_cube()
        face_centre = _v((1, 0), (1, 0), (0, 0))
        assert is_in_hull(face_centre, vs, faces)

    @pytest.mark.parametrize(
        "x,y,z",
        [
            ((4, 0), (1, 0), (1, 0)),    # past x = 2
            ((-1, 0), (1, 0), (1, 0)),   # past x = 0
            ((1, 0), (3, 0), (1, 0)),    # past y = 2
            ((1, 0), (1, 0), (-1, 0)),   # past z = 0
        ],
    )
    def test_outside_points_rejected(
        self,
        x: tuple[int, int],
        y: tuple[int, int],
        z: tuple[int, int],
    ) -> None:
        vs, faces = _unit_cube()
        assert not is_in_hull(_v(x, y, z), vs, faces)

    def test_sign_agnostic_winding(self) -> None:
        # Reverse every face's winding; is_in_hull must still decide
        # interior/exterior correctly, because the inward sign is
        # re-derived per face from a reference hull vertex.
        vs, faces = _unit_cube()
        flipped = [(k, j, i) for (i, j, k) in faces]
        inside = _v((1, 0), (1, 0), (1, 0))
        outside = _v((4, 0), (1, 0), (1, 0))
        assert is_in_hull(inside, vs, flipped)
        assert not is_in_hull(outside, vs, flipped)

    def test_rejects_fewer_than_four_vertices(self) -> None:
        vs = [
            _v((0, 0), (0, 0), (0, 0)),
            _v((2, 0), (0, 0), (0, 0)),
            _v((0, 0), (2, 0), (0, 0)),
        ]
        with pytest.raises(ValueError):
            is_in_hull(_v((1, 0), (1, 0), (0, 0)), vs, [(0, 1, 2)])

    def test_degenerate_coplanar_hull_raises(self) -> None:
        # All four vertices on the plane z = 0; the face (0, 1, 2) has
        # no off-plane reference vertex.
        vs = [
            _v((0, 0), (0, 0), (0, 0)),
            _v((2, 0), (0, 0), (0, 0)),
            _v((0, 0), (2, 0), (0, 0)),
            _v((2, 0), (2, 0), (0, 0)),
        ]
        with pytest.raises(ValueError):
            is_in_hull(_v((1, 0), (1, 0), (0, 0)), vs, [(0, 1, 2)])


# -- check_euler -------------------------------------------------------


class TestCheckEuler:
    def test_tetrahedron(self) -> None:
        # V = 4, E = 6, F = 4; χ = 2.
        faces = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]
        assert check_euler(4, faces)

    def test_triangulated_cube(self) -> None:
        # V = 8, E = 18 (12 cube edges + 6 face diagonals), F = 12; χ = 2.
        _, faces = _unit_cube()
        assert check_euler(8, faces)

    def test_icosahedron_combinatorics(self) -> None:
        # The regular icosahedron has V = 12, E = 30, F = 20.
        # We just need any index-triple combinatorics that satisfy
        # V - E + F = 2 with V = 12 and F = 20.
        faces = [
            (0, 1, 2), (0, 2, 3), (0, 3, 4), (0, 4, 5), (0, 5, 1),
            (6, 7, 8), (6, 8, 9), (6, 9, 10), (6, 10, 11), (6, 11, 7),
            (1, 7, 2), (2, 8, 3), (3, 9, 4), (4, 10, 5), (5, 11, 1),
            (7, 2, 8), (8, 3, 9), (9, 4, 10), (10, 5, 11), (11, 1, 7),
        ]
        assert check_euler(12, faces)

    def test_rejects_wrong_relation(self) -> None:
        # A single triangle with three vertices: V = 3, E = 3, F = 1.
        # χ = 3 - 3 + 1 = 1, not 2.
        assert not check_euler(3, [(0, 1, 2)])

    def test_duplicated_face_deduplicates_edges_correctly(self) -> None:
        # Each edge in the set is counted once regardless of how many
        # faces share it. This is what makes the V - E + F relation
        # hold for the triangulated-sphere combinatorics.
        faces = [(0, 1, 2), (0, 1, 2)]  # degenerate: same face twice
        # V = 3, E = 3, F = 2; χ = 3 - 3 + 2 = 2. Numerically χ = 2 but
        # the combinatorics aren't a valid triangulation. check_euler
        # answers the relation as-asked, not the validity of the
        # complex; that's a separate check if ever needed.
        assert check_euler(3, faces)


# -- Polyhedron fixtures -----------------------------------------------


def _raw_tetrahedron() -> tuple[list[Vec3], list[tuple[int, int, int]]]:
    """Right-angled tetrahedron at the origin, stored ×2.

    Real vertices: (0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1); stored
    as Vec3 with components doubled. Face cycles are oriented with
    outward-pointing normals (each cycle is CCW as viewed from
    outside the tetrahedron).
    """
    vs = [
        _v((0, 0), (0, 0), (0, 0)),   # v0 origin
        _v((2, 0), (0, 0), (0, 0)),   # v1 +x
        _v((0, 0), (2, 0), (0, 0)),   # v2 +y
        _v((0, 0), (0, 0), (2, 0)),   # v3 +z
    ]
    faces = [
        (1, 2, 3),   # opposite v0; normal +x+y+z
        (0, 2, 3),   # opposite v1; normal -x
        (0, 3, 2),   # opposite v2; normal -y
        (0, 2, 1),   # opposite v3; normal -z
    ]
    return vs, faces


def _canonical_tetra_vertices() -> tuple[Vec3, ...]:
    """Canonical-order vertex tuple for the tetrahedron fixture.

    Sorted by the raw-coefficient lex key used inside ``Polyhedron``,
    suitable for direct ``Polyhedron(...)`` construction in tests that
    need to bypass ``from_raw``.
    """
    return (
        _v((0, 0), (0, 0), (0, 0)),   # origin
        _v((0, 0), (0, 0), (2, 0)),   # +z
        _v((0, 0), (2, 0), (0, 0)),   # +y
        _v((2, 0), (0, 0), (0, 0)),   # +x
    )


# -- Polyhedron construction and canonical form ------------------------


class TestPolyhedronConstruction:
    def test_from_raw_sorts_vertices(self) -> None:
        # Reverse the raw vertex order; from_raw must sort back.
        vs, fs = _raw_tetrahedron()
        reversed_vs = list(reversed(vs))
        # Remap faces to match the reversed indexing.
        n = len(vs)
        reversed_fs = [tuple(n - 1 - i for i in f) for f in fs]
        P = Polyhedron.from_raw(reversed_vs, reversed_fs)
        Q = Polyhedron.from_raw(vs, fs)
        assert P == Q
        assert hash(P) == hash(Q)

    def test_from_raw_remaps_face_indices(self) -> None:
        vs, fs = _raw_tetrahedron()
        P = Polyhedron.from_raw(vs, fs)
        # Every face index stays within range after remap.
        for face in P.faces:
            for idx in face:
                assert 0 <= idx < len(P.vertices)

    def test_vertices_in_canonical_order(self) -> None:
        vs, fs = _raw_tetrahedron()
        P = Polyhedron.from_raw(vs, fs)
        # Raw coefficients in ascending lexicographic order.
        keys = [(v.x.a, v.x.b, v.y.a, v.y.b, v.z.a, v.z.b) for v in P.vertices]
        assert keys == sorted(keys)

    def test_face_cycles_start_at_min_index(self) -> None:
        vs, fs = _raw_tetrahedron()
        P = Polyhedron.from_raw(vs, fs)
        for face in P.faces:
            assert face[0] == min(face)

    def test_faces_sorted(self) -> None:
        vs, fs = _raw_tetrahedron()
        P = Polyhedron.from_raw(vs, fs)
        assert list(P.faces) == sorted(P.faces)

    def test_rejects_fewer_than_four_vertices(self) -> None:
        with pytest.raises(ValueError, match="at least 4 vertices"):
            Polyhedron(
                vertices=(
                    _v((0, 0), (0, 0), (0, 0)),
                    _v((2, 0), (0, 0), (0, 0)),
                    _v((0, 0), (2, 0), (0, 0)),
                ),
                faces=((0, 1, 2),),
            )

    def test_rejects_duplicate_vertices(self) -> None:
        with pytest.raises(ValueError, match="pairwise distinct"):
            Polyhedron(
                vertices=(
                    _v((0, 0), (0, 0), (0, 0)),
                    _v((0, 0), (0, 0), (0, 0)),
                    _v((0, 0), (2, 0), (0, 0)),
                    _v((2, 0), (0, 0), (0, 0)),
                ),
                faces=((0, 1, 2),),
            )

    def test_rejects_unsorted_vertices(self) -> None:
        # Largest vertex first — not canonical.
        with pytest.raises(ValueError, match="canonical order"):
            Polyhedron(
                vertices=(
                    _v((2, 0), (0, 0), (0, 0)),
                    _v((0, 0), (0, 0), (0, 0)),
                    _v((0, 0), (0, 0), (2, 0)),
                    _v((0, 0), (2, 0), (0, 0)),
                ),
                faces=((0, 1, 2),),
            )

    def test_rejects_degenerate_face(self) -> None:
        vs = _canonical_tetra_vertices()
        with pytest.raises(ValueError, match="fewer than 3"):
            Polyhedron(vertices=vs, faces=((0, 1),))

    def test_rejects_repeated_face_indices(self) -> None:
        vs = _canonical_tetra_vertices()
        with pytest.raises(ValueError, match="repeated vertex indices"):
            Polyhedron(vertices=vs, faces=((0, 1, 1),))

    def test_rejects_out_of_range_face_index(self) -> None:
        vs = _canonical_tetra_vertices()
        with pytest.raises(ValueError, match="out of range"):
            Polyhedron(vertices=vs, faces=((0, 1, 99),))

    def test_rejects_non_canonical_cycle(self) -> None:
        # (1, 2, 0) — valid cycle rotation but doesn't start at min.
        vs, fs = _raw_tetrahedron()
        canonical = Polyhedron.from_raw(vs, fs)
        bad_face = tuple(reversed(canonical.faces[0]))  # also not canonical
        with pytest.raises(ValueError, match="canonical cycle order"):
            Polyhedron(vertices=canonical.vertices, faces=(bad_face,))

    def test_rejects_duplicate_faces(self) -> None:
        vs = _canonical_tetra_vertices()
        with pytest.raises(ValueError, match="Duplicate face"):
            Polyhedron(vertices=vs, faces=((0, 1, 2), (0, 1, 2)))

    def test_rejects_unsorted_faces(self) -> None:
        vs = _canonical_tetra_vertices()
        with pytest.raises(ValueError, match="must be sorted"):
            Polyhedron(vertices=vs, faces=((1, 2, 3), (0, 1, 2)))


# -- Polyhedron equality and hashability ------------------------------


class TestPolyhedronEquality:
    def test_from_raw_is_canonical_regardless_of_vertex_order(self) -> None:
        vs, fs = _raw_tetrahedron()
        A = Polyhedron.from_raw(vs, fs)
        # Scramble vertex order differently; remap faces accordingly.
        perm = [2, 0, 3, 1]
        inv = {old: new for new, old in enumerate(perm)}
        scrambled_vs = [vs[p] for p in perm]
        scrambled_fs = [tuple(inv[i] for i in f) for f in fs]
        B = Polyhedron.from_raw(scrambled_vs, scrambled_fs)
        assert A == B
        assert hash(A) == hash(B)

    def test_different_polyhedra_distinct(self) -> None:
        vs, fs = _raw_tetrahedron()
        T = Polyhedron.from_raw(vs, fs)
        # Translate one vertex; different polyhedron.
        moved_vs = list(vs)
        moved_vs[3] = _v((0, 0), (0, 0), (4, 0))  # +z moved further
        moved = Polyhedron.from_raw(moved_vs, fs)
        assert T != moved

    def test_hashable_in_set(self) -> None:
        vs, fs = _raw_tetrahedron()
        T = Polyhedron.from_raw(vs, fs)
        assert {T, T} == {T}


# -- Polyhedron isometry action ---------------------------------------


class TestPolyhedronApply:
    def test_identity_fixes_polyhedron(self) -> None:
        vs, fs = _raw_tetrahedron()
        T = Polyhedron.from_raw(vs, fs)
        assert T.apply(Rotation.identity()) == T

    def test_apply_preserves_vertex_and_face_counts(self) -> None:
        vs, fs = _raw_tetrahedron()
        T = Polyhedron.from_raw(vs, fs)
        for g in (ROT_2, ROT_3, ROT_5):
            R = T.apply(g)
            assert len(R.vertices) == len(T.vertices)
            assert len(R.faces) == len(T.faces)

    def test_apply_then_inverse_returns_original(self) -> None:
        vs, fs = _raw_tetrahedron()
        T = Polyhedron.from_raw(vs, fs)
        # Sample a subset of I for speed; the full group is exercised
        # in the RTH integration test.
        for g in ICOSAHEDRAL[:12]:
            assert T.apply(g).apply(g.inverse()) == T

    def test_apply_composes_as_hg_on_positions(self) -> None:
        """Applying g then h equals applying h ∘ g to positions.

        ``Rotation.compose(self, other)`` returns ``self ∘ other``
        (function composition). So ``T.apply(g).apply(h)`` applies g
        then h to each position, which is ``h ∘ g`` applied once.
        """
        vs, fs = _raw_tetrahedron()
        T = Polyhedron.from_raw(vs, fs)
        g, h = ROT_2, ROT_3
        assert T.apply(g).apply(h) == T.apply(h.compose(g))

    def test_apply_generically_permutes_vertices(self) -> None:
        """A generic rotation should permute the canonical vertex
        ordering (the rotated polyhedron's canonical vertex list is
        not the same tuple as the original's, element-wise)."""
        vs, fs = _raw_tetrahedron()
        T = Polyhedron.from_raw(vs, fs)
        R = T.apply(ROT_5)
        # Tetrahedron is not symmetric under ROT_5, so positions change:
        assert R.vertices != T.vertices


# -- _same_oriented_plane ---------------------------------------------


class TestSameOrientedPlane:
    def test_identical_plane(self) -> None:
        n = Vec3(ZPhi(0, 0), ZPhi(0, 0), ZPhi(4, 0))
        d = ZPhi(0, 0)
        assert _same_oriented_plane(n, d, n, d)

    def test_positive_scalar_multiple(self) -> None:
        # Both planes are z = 0 with same-direction +z normal.
        n1 = Vec3(ZPhi(0, 0), ZPhi(0, 0), ZPhi(2, 0))
        d1 = ZPhi(0, 0)
        n2 = Vec3(ZPhi(0, 0), ZPhi(0, 0), ZPhi(6, 0))  # 3x n1
        d2 = ZPhi(0, 0)
        assert _same_oriented_plane(n1, d1, n2, d2)

    def test_antiparallel_normals_rejected(self) -> None:
        n1 = Vec3(ZPhi(0, 0), ZPhi(0, 0), ZPhi(4, 0))
        n2 = Vec3(ZPhi(0, 0), ZPhi(0, 0), ZPhi(-4, 0))
        d1 = d2 = ZPhi(0, 0)
        assert not _same_oriented_plane(n1, d1, n2, d2)

    def test_parallel_different_offset_rejected(self) -> None:
        # z = 0 vs z = 2, both +z normals: different planes.
        n = Vec3(ZPhi(0, 0), ZPhi(0, 0), ZPhi(4, 0))
        d1 = ZPhi(0, 0)
        d2 = ZPhi(8, 0)
        assert not _same_oriented_plane(n, d1, n, d2)

    def test_nonparallel_rejected(self) -> None:
        n1 = Vec3(ZPhi(0, 0), ZPhi(0, 0), ZPhi(4, 0))    # +z
        n2 = Vec3(ZPhi(4, 0), ZPhi(0, 0), ZPhi(0, 0))    # +x
        d = ZPhi(0, 0)
        assert not _same_oriented_plane(n1, d, n2, d)


# -- _reconstruct_polygon_cycle --------------------------------------


class TestReconstructPolygon:
    def test_single_triangle_passthrough(self) -> None:
        # A one-triangle "group" has no interior edges; its boundary
        # is the triangle itself.
        cycle = _reconstruct_polygon_cycle([(0, 1, 2)])
        assert cycle == (0, 1, 2)

    def test_two_triangles_to_quad(self) -> None:
        # Triangles (0, 1, 2) and (0, 2, 3) share edge (0, 2) ↔ (2, 0);
        # interior edge cancels; boundary is the quad 0-1-2-3.
        cycle = _reconstruct_polygon_cycle([(0, 1, 2), (0, 2, 3)])
        assert cycle == (0, 1, 2, 3)

    def test_three_triangles_to_pentagon(self) -> None:
        # Fan-triangulated pentagon from vertex 0:
        # (0, 1, 2), (0, 2, 3), (0, 3, 4).
        # Interior edges (0, 2) and (0, 3) cancel; boundary is
        # 0-1-2-3-4.
        cycle = _reconstruct_polygon_cycle([(0, 1, 2), (0, 2, 3), (0, 3, 4)])
        assert cycle == (0, 1, 2, 3, 4)

    def test_duplicate_oriented_edge_raises(self) -> None:
        # Same directed edge (0, 1) in two triangles — non-manifold.
        with pytest.raises(ValueError, match="Duplicate oriented edge"):
            _reconstruct_polygon_cycle([(0, 1, 2), (0, 1, 3)])


# -- _orient_simplices_outward ---------------------------------------


class TestOrientSimplicesOutward:
    def test_fixes_inward_normal(self) -> None:
        # Flat square base at z=0 of a pyramid; the two simplices on
        # it should reorient consistently to outward (-z for z=0).
        # Construct 5 vertices: 4 base corners + 1 apex.
        hull_vs = [
            _v((0, 0), (0, 0), (0, 0)),    # 0
            _v((2, 0), (0, 0), (0, 0)),    # 1
            _v((0, 0), (2, 0), (0, 0)),    # 2
            _v((2, 0), (2, 0), (0, 0)),    # 3
            _v((1, 0), (1, 0), (2, 0)),    # 4 apex
        ]
        # Two base simplices, one inward-oriented, one outward.
        simplices = [
            (0, 1, 3),   # (b-a)×(c-a) gives +z — INWARD for z=0 face
            (0, 3, 2),   # gives -z — OUTWARD
        ]
        oriented = _orient_simplices_outward(hull_vs, simplices)
        # Both should have outward (-z) normals after re-orientation:
        for (i, j, k) in oriented:
            n, _ = _triangle_plane(hull_vs[i], hull_vs[j], hull_vs[k])
            # For base simplices, outward normal has negative z.
            assert n.z < ZPhi(0, 0)

    def test_keeps_already_outward(self) -> None:
        hull_vs = [
            _v((0, 0), (0, 0), (0, 0)),
            _v((2, 0), (0, 0), (0, 0)),
            _v((0, 0), (2, 0), (0, 0)),
            _v((1, 0), (1, 0), (2, 0)),
        ]
        # Apex simplex (0, 1, 2) has +z normal; this is INWARD for the
        # base of the tet, since the apex is above z=0. Reorient keeps
        # outward, which means swapping.
        oriented = _orient_simplices_outward(hull_vs, [(0, 1, 2)])
        # The resulting triangle's outward normal should be -z.
        (i, j, k) = oriented[0]
        n, _ = _triangle_plane(hull_vs[i], hull_vs[j], hull_vs[k])
        assert n.z < ZPhi(0, 0)


# -- Polyhedron.from_vertices ----------------------------------------


def _stored_cube_vertices() -> list[Vec3]:
    """Stored-×2 unit cube vertex list (8 vertices)."""
    return [
        _v((2 * i, 0), (2 * j, 0), (2 * k, 0))
        for i in (0, 1) for j in (0, 1) for k in (0, 1)
    ]


def _stored_rhombic_dodec_vertices() -> list[Vec3]:
    """Stored-×2 rhombic dodecahedron: 14 vertices.

    Real vertices: 8 cube corners at (±1, ±1, ±1) (degree-3), plus 6
    face centres at (±2, 0, 0) and cyclic permutations (degree-4).
    Stored form multiplies every component by 2.
    """
    vs: list[Vec3] = []
    for x in (-1, 1):
        for y in (-1, 1):
            for z in (-1, 1):
                vs.append(_v((2 * x, 0), (2 * y, 0), (2 * z, 0)))
    for sign in (-1, 1):
        for pos in range(3):
            coords = [0, 0, 0]
            coords[pos] = 4 * sign   # real ±2 → stored ±4
            vs.append(_v((coords[0], 0), (coords[1], 0), (coords[2], 0)))
    return vs


class TestPolyhedronFromVertices:
    def test_tetrahedron(self) -> None:
        vs = list(_canonical_tetra_vertices())
        P = Polyhedron.from_vertices(vs)
        assert len(P.vertices) == 4
        assert len(P.faces) == 4
        for face in P.faces:
            assert len(face) == 3   # all triangular

    def test_cube_merges_to_six_quadrilaterals(self) -> None:
        P = Polyhedron.from_vertices(_stored_cube_vertices())
        assert len(P.vertices) == 8
        assert len(P.faces) == 6
        for face in P.faces:
            assert len(face) == 4   # every face is a quadrilateral

    def test_cube_euler_relation(self) -> None:
        # Polygonal Euler: V - E + F = 2, with edge count derived from
        # the polygonal face list.
        P = Polyhedron.from_vertices(_stored_cube_vertices())
        edges: set[tuple[int, int]] = set()
        for face in P.faces:
            n = len(face)
            for idx in range(n):
                u, v = face[idx], face[(idx + 1) % n]
                edges.add((min(u, v), max(u, v)))
        # Cube: V = 8, E = 12, F = 6.
        assert len(P.vertices) == 8
        assert len(edges) == 12
        assert len(P.faces) == 6
        assert len(P.vertices) - len(edges) + len(P.faces) == 2

    def test_rhombic_dodecahedron_counts(self) -> None:
        # RD: 14 vertices, 24 edges, 12 rhombic faces; V - E + F = 2.
        P = Polyhedron.from_vertices(_stored_rhombic_dodec_vertices())
        assert len(P.vertices) == 14
        assert len(P.faces) == 12
        for face in P.faces:
            assert len(face) == 4   # rhombic
        edges: set[tuple[int, int]] = set()
        for face in P.faces:
            n = len(face)
            for idx in range(n):
                u, v = face[idx], face[(idx + 1) % n]
                edges.add((min(u, v), max(u, v)))
        assert len(edges) == 24
        assert len(P.vertices) - len(edges) + len(P.faces) == 2

    def test_filters_interior_vertex(self) -> None:
        cube = _stored_cube_vertices()
        interior = _v((1, 0), (1, 0), (1, 0))   # real (0.5, 0.5, 0.5)
        P = Polyhedron.from_vertices(cube + [interior])
        assert len(P.vertices) == 8
        assert interior not in P.vertices

    def test_rejects_fewer_than_four_vertices(self) -> None:
        vs = [
            _v((0, 0), (0, 0), (0, 0)),
            _v((2, 0), (0, 0), (0, 0)),
            _v((0, 0), (2, 0), (0, 0)),
        ]
        with pytest.raises(ValueError, match="at least 4"):
            Polyhedron.from_vertices(vs)

    def test_hull_result_is_canonical(self) -> None:
        # Shuffled input order should not affect the output.
        cube = _stored_cube_vertices()
        shuffled = [cube[i] for i in (7, 2, 4, 0, 3, 5, 1, 6)]
        A = Polyhedron.from_vertices(cube)
        B = Polyhedron.from_vertices(shuffled)
        assert A == B

    def test_merge_round_trip_via_from_vertices(self) -> None:
        # End-to-end: every cube face has 4 distinct vertex indices,
        # and every pair of adjacent indices in a face cycle shares
        # an edge in the polytope (consecutive face cycle ⇒ edge).
        P = Polyhedron.from_vertices(_stored_cube_vertices())
        for face in P.faces:
            assert len(set(face)) == len(face)
            for idx in range(len(face)):
                u, v = face[idx], face[(idx + 1) % len(face)]
                assert u != v

    def test_rotation_action_commutes_with_hull(self) -> None:
        # from_vertices(rotated) == from_vertices(original).apply(g)
        # for g in the icosahedral group. Rotation preserves which
        # points are on the hull.
        cube = _stored_cube_vertices()
        P = Polyhedron.from_vertices(cube)
        g = ROT_3
        rotated_inputs = [g.apply(v) for v in cube]
        Q = Polyhedron.from_vertices(rotated_inputs)
        assert Q == P.apply(g)


# -- _merge_coplanar_faces end-to-end --------------------------------


class TestMergeCoplanarFaces:
    def test_trivial_single_triangle(self) -> None:
        # Four tetrahedron vertices; one triangle → stays one triangle.
        hull_vs = list(_canonical_tetra_vertices())
        merged = _merge_coplanar_faces(hull_vs, [(0, 1, 2)])
        assert merged == [(0, 1, 2)]

    def test_two_coplanar_triangles_merge(self) -> None:
        # Square base of a pyramid triangulated as two triangles;
        # they share a diagonal and merge.
        hull_vs = [
            _v((0, 0), (0, 0), (0, 0)),
            _v((2, 0), (0, 0), (0, 0)),
            _v((2, 0), (2, 0), (0, 0)),
            _v((0, 0), (2, 0), (0, 0)),
            _v((1, 0), (1, 0), (2, 0)),   # apex; off-plane for merge test
        ]
        # Both base triangles outward-oriented (-z):
        tris = [(0, 3, 2), (0, 2, 1)]
        merged = _merge_coplanar_faces(hull_vs, tris)
        # One merged quadrilateral face.
        assert len(merged) == 1
        assert len(merged[0]) == 4
        assert set(merged[0]) == {0, 1, 2, 3}
