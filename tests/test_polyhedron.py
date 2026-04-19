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
