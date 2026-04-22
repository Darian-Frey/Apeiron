"""Tests for apeiron.corona.

Scope for this commit (sub-commit A of the corona.py plan): the data
model — ``PlacedTile``, ``CoronaConfig``, canonical neighbour
ordering, equality and hashing. Geometric primitives (face-to-face
placement, interior-overlap, angular-defect) and the BFS engine
itself come in sub-commits B through D.

Sanity fixtures (cube, rhombic dodecahedron — CLAUDE.md §8) are the
acceptance oracles for ``corona_1`` in sub-commits C and D, not this
one; here they appear only as handy ``Polyhedron`` instances for the
data-model tests.
"""

from __future__ import annotations

import pytest

from apeiron.corona import (
    CoronaConfig,
    Edge,
    PlacedTile,
    Vertex,
    _edges_of,
    face_to_face_placements,
    find_rotation,
    has_interior_overlap,
    incidence_defect,
)
from apeiron.polyhedron import Polyhedron
from apeiron.symmetry import ICOSAHEDRAL, ROT_2, ROT_3, ROT_5, Mat3, Rotation, Vec3
from apeiron.zphi import ZPhi


def _v(x: tuple[int, int], y: tuple[int, int], z: tuple[int, int]) -> Vec3:
    return Vec3(ZPhi(*x), ZPhi(*y), ZPhi(*z))


def _unit_cube() -> Polyhedron:
    """Stored-×2 unit cube with corners at (0, 2)×(0, 2)×(0, 2)."""
    corners = [
        _v((2 * i, 0), (2 * j, 0), (2 * k, 0))
        for i in (0, 1) for j in (0, 1) for k in (0, 1)
    ]
    return Polyhedron.from_vertices(corners)


# -- PlacedTile --------------------------------------------------------


class TestPlacedTile:
    def test_construction(self) -> None:
        t = _v((4, 0), (0, 0), (0, 0))
        g = ROT_3
        pt = PlacedTile(translation=t, rotation=g)
        assert pt.translation == t
        assert pt.rotation == g

    def test_equality(self) -> None:
        t = _v((4, 0), (0, 0), (0, 0))
        pt1 = PlacedTile(translation=t, rotation=Rotation.identity())
        pt2 = PlacedTile(translation=t, rotation=Rotation.identity())
        assert pt1 == pt2
        assert hash(pt1) == hash(pt2)

    def test_distinct_on_translation(self) -> None:
        pt1 = PlacedTile(_v((4, 0), (0, 0), (0, 0)), Rotation.identity())
        pt2 = PlacedTile(_v((-4, 0), (0, 0), (0, 0)), Rotation.identity())
        assert pt1 != pt2

    def test_distinct_on_rotation(self) -> None:
        t = _v((4, 0), (0, 0), (0, 0))
        pt1 = PlacedTile(t, Rotation.identity())
        pt2 = PlacedTile(t, ROT_2)
        assert pt1 != pt2

    def test_hashable_in_set(self) -> None:
        pt = PlacedTile(_v((4, 0), (0, 0), (0, 0)), Rotation.identity())
        assert {pt, pt} == {pt}


# -- CoronaConfig from_neighbours (happy paths) -----------------------


class TestCoronaConfigFromNeighbours:
    def test_basic_construction(self) -> None:
        cube = _unit_cube()
        pt = PlacedTile(_v((4, 0), (0, 0), (0, 0)), Rotation.identity())
        c = CoronaConfig.from_neighbours(cube, [pt])
        assert c.central == cube
        assert c.neighbours == (pt,)

    def test_deduplicates_identical_neighbours(self) -> None:
        cube = _unit_cube()
        pt = PlacedTile(_v((4, 0), (0, 0), (0, 0)), Rotation.identity())
        # Triplicate input; from_neighbours must collapse to one.
        c = CoronaConfig.from_neighbours(cube, [pt, pt, pt])
        assert c.neighbours == (pt,)

    def test_sorts_by_canonical_key(self) -> None:
        cube = _unit_cube()
        pt_neg = PlacedTile(_v((-4, 0), (0, 0), (0, 0)), Rotation.identity())
        pt_pos = PlacedTile(_v((4, 0), (0, 0), (0, 0)), Rotation.identity())
        c = CoronaConfig.from_neighbours(cube, [pt_pos, pt_neg])
        # pt_neg has smaller translation.x.a, so it comes first.
        assert c.neighbours == (pt_neg, pt_pos)

    def test_sort_breaks_ties_on_rotation(self) -> None:
        # Same translation, different rotations: tie-break uses the
        # rotation matrix coefficients.
        cube = _unit_cube()
        t = _v((0, 0), (0, 0), (0, 0))
        pt_id = PlacedTile(t, Rotation.identity())
        pt_rot = PlacedTile(t, ROT_2)
        c1 = CoronaConfig.from_neighbours(cube, [pt_id, pt_rot])
        c2 = CoronaConfig.from_neighbours(cube, [pt_rot, pt_id])
        # Whichever rotation has the smaller matrix key comes first;
        # from_neighbours is order-invariant on input.
        assert c1 == c2
        assert c1.neighbours[0] in (pt_id, pt_rot)

    def test_empty_neighbour_iterable(self) -> None:
        cube = _unit_cube()
        c = CoronaConfig.from_neighbours(cube, [])
        assert c.neighbours == ()


# -- CoronaConfig direct construction and validation ------------------


class TestCoronaConfigValidation:
    def test_canonical_direct_construction_allowed(self) -> None:
        cube = _unit_cube()
        pt_a = PlacedTile(_v((-4, 0), (0, 0), (0, 0)), Rotation.identity())
        pt_b = PlacedTile(_v((4, 0), (0, 0), (0, 0)), Rotation.identity())
        # Already sorted by translation.x.a.
        c = CoronaConfig(central=cube, neighbours=(pt_a, pt_b))
        assert c.neighbours == (pt_a, pt_b)

    def test_rejects_unsorted_neighbours(self) -> None:
        cube = _unit_cube()
        pt_a = PlacedTile(_v((-4, 0), (0, 0), (0, 0)), Rotation.identity())
        pt_b = PlacedTile(_v((4, 0), (0, 0), (0, 0)), Rotation.identity())
        with pytest.raises(ValueError, match="canonical order"):
            CoronaConfig(central=cube, neighbours=(pt_b, pt_a))

    def test_rejects_duplicate_neighbours(self) -> None:
        cube = _unit_cube()
        pt = PlacedTile(_v((4, 0), (0, 0), (0, 0)), Rotation.identity())
        with pytest.raises(ValueError, match="duplicate"):
            CoronaConfig(central=cube, neighbours=(pt, pt))


# -- CoronaConfig equality and hashing --------------------------------


class TestCoronaConfigEquality:
    def test_equal_after_order_normalisation(self) -> None:
        cube = _unit_cube()
        pt_a = PlacedTile(_v((-4, 0), (0, 0), (0, 0)), Rotation.identity())
        pt_b = PlacedTile(_v((4, 0), (0, 0), (0, 0)), Rotation.identity())
        c1 = CoronaConfig.from_neighbours(cube, [pt_a, pt_b])
        c2 = CoronaConfig.from_neighbours(cube, [pt_b, pt_a])
        assert c1 == c2
        assert hash(c1) == hash(c2)

    def test_distinct_neighbour_sets_are_unequal(self) -> None:
        cube = _unit_cube()
        pt_a = PlacedTile(_v((-4, 0), (0, 0), (0, 0)), Rotation.identity())
        pt_b = PlacedTile(_v((4, 0), (0, 0), (0, 0)), Rotation.identity())
        c1 = CoronaConfig.from_neighbours(cube, [pt_a])
        c2 = CoronaConfig.from_neighbours(cube, [pt_b])
        assert c1 != c2

    def test_distinct_centrals_are_unequal(self) -> None:
        cube_a = _unit_cube()
        # Build a distinct polyhedron by shifting the cube along +x.
        shifted_corners = [
            _v((2 * i + 4, 0), (2 * j, 0), (2 * k, 0))
            for i in (0, 1) for j in (0, 1) for k in (0, 1)
        ]
        cube_b = Polyhedron.from_vertices(shifted_corners)
        assert cube_a != cube_b   # precondition
        pt = PlacedTile(_v((4, 0), (0, 0), (0, 0)), Rotation.identity())
        c1 = CoronaConfig.from_neighbours(cube_a, [pt])
        c2 = CoronaConfig.from_neighbours(cube_b, [pt])
        assert c1 != c2

    def test_hashable_in_set(self) -> None:
        cube = _unit_cube()
        pt = PlacedTile(_v((4, 0), (0, 0), (0, 0)), Rotation.identity())
        c = CoronaConfig.from_neighbours(cube, [pt])
        assert {c, c} == {c}

    def test_equal_configs_survive_through_frozen_set(self) -> None:
        cube = _unit_cube()
        pt_a = PlacedTile(_v((-4, 0), (0, 0), (0, 0)), Rotation.identity())
        pt_b = PlacedTile(_v((4, 0), (0, 0), (0, 0)), Rotation.identity())
        c1 = CoronaConfig.from_neighbours(cube, [pt_a, pt_b])
        c2 = CoronaConfig.from_neighbours(cube, [pt_b, pt_a])
        assert len({c1, c2}) == 1


# -- Smoke check: the canonical-ordering surface ------------------------


class TestCanonicalOrderingAcrossIGroup:
    """The canonical key is a total order on PlacedTile; any sequence
    of PlacedTiles drawn from the icosahedral-group action produces a
    deterministic neighbour tuple.

    This is not yet a BFS-dedup test (that's sub-commit C). It's a
    precondition: verify that sorting is stable and consistent across
    the 60 elements of I acting on a single base PlacedTile.
    """

    def test_all_rotations_produce_distinct_keys_at_fixed_translation(self) -> None:
        cube = _unit_cube()
        t = _v((4, 0), (0, 0), (0, 0))
        placements = [PlacedTile(t, g) for g in ICOSAHEDRAL]
        # 60 distinct rotations → 60 distinct PlacedTiles at this
        # translation, hence 60 distinct canonical keys.
        c = CoronaConfig.from_neighbours(cube, placements)
        assert len(c.neighbours) == 60

    def test_sorting_is_deterministic_across_permutations(self) -> None:
        cube = _unit_cube()
        t = _v((4, 0), (0, 0), (0, 0))
        placements = [PlacedTile(t, g) for g in ICOSAHEDRAL]
        c1 = CoronaConfig.from_neighbours(cube, placements)
        c2 = CoronaConfig.from_neighbours(cube, reversed(placements))
        assert c1.neighbours == c2.neighbours


# -- find_rotation -----------------------------------------------------


class TestFindRotation:
    def test_finds_identity(self) -> None:
        g = find_rotation(Rotation.identity().matrix)
        assert g == Rotation.identity()

    def test_finds_each_generator(self) -> None:
        for gen in (ROT_2, ROT_3, ROT_5):
            found = find_rotation(gen.matrix)
            assert found == gen

    def test_finds_every_element_of_i(self) -> None:
        # Round-trip: each rotation in ICOSAHEDRAL is found by its own
        # matrix.
        for g in ICOSAHEDRAL:
            assert find_rotation(g.matrix) == g

    def test_returns_none_for_non_i_matrix(self) -> None:
        # The identity's zero matrix is not a rotation — not in I.
        zero_mat = Mat3.from_ints(((0, 0),) * 9)
        assert find_rotation(zero_mat) is None

    def test_returns_none_for_non_rotation_matrix(self) -> None:
        # 3x diagonal matrix (not a rotation — determinant 27, not 1).
        # 3*I in real coords → 6*I in ×2 storage.
        scaled = Mat3.from_ints((
            (6, 0), (0, 0), (0, 0),
            (0, 0), (6, 0), (0, 0),
            (0, 0), (0, 0), (6, 0),
        ))
        assert find_rotation(scaled) is None


# -- face_to_face_placements ------------------------------------------


class TestFaceToFacePlacements:
    def test_rejects_out_of_range_face_index(self) -> None:
        cube = _unit_cube()
        with pytest.raises(IndexError, match="out of range"):
            face_to_face_placements(cube, 999)
        with pytest.raises(IndexError, match="out of range"):
            face_to_face_placements(cube, -1)

    def test_cube_each_face_has_placements(self) -> None:
        cube = _unit_cube()
        for idx in range(len(cube.faces)):
            placements = face_to_face_placements(cube, idx)
            assert len(placements) >= 1, f"face {idx} has no placements"

    def test_cube_placement_count_is_12_per_face(self) -> None:
        # For the axis-aligned cube with icosahedral-compatible choice
        # of x-axis as a 2-fold axis of I, |I ∩ O| = 12: 1 identity +
        # 3 face-axis 2-folds + 8 body-diagonal 3-folds. Each of these
        # 12 elements yields one face-to-face placement per central
        # face. The 6-face × 4-cyclic-shift = 24 iteration counts each
        # placement twice (face stabiliser of order 2 in I ∩ O).
        cube = _unit_cube()
        for idx in range(len(cube.faces)):
            assert len(face_to_face_placements(cube, idx)) == 12

    def test_cube_placements_are_distinct(self) -> None:
        cube = _unit_cube()
        placements = face_to_face_placements(cube, 0)
        assert len(set(placements)) == len(placements)

    def test_cube_placements_in_canonical_order(self) -> None:
        # face_to_face_placements returns a list sorted by
        # _placed_tile_key.
        from apeiron.corona import _placed_tile_key
        cube = _unit_cube()
        placements = face_to_face_placements(cube, 0)
        keys = [_placed_tile_key(p) for p in placements]
        assert keys == sorted(keys)

    def test_placement_face_coincides_with_target(self) -> None:
        # Each returned placement must actually produce a copy of the
        # central with some face coinciding with the reversed central
        # face.
        cube = _unit_cube()
        central_face = cube.faces[0]
        central_vertices = tuple(cube.vertices[i] for i in central_face)
        # Target cycle: reversed central face.
        target = {central_vertices[0]} | set(central_vertices[1:])
        for placement in face_to_face_placements(cube, 0):
            # Compute placed vertices.
            placed = [
                placement.rotation.apply(v) + placement.translation
                for v in cube.vertices
            ]
            # Some face of the placed cube has the same vertex set as
            # the central face.
            placed_set = set(placed)
            assert target.issubset(placed_set), (
                "placement does not contain the target face's vertices"
            )

    def test_rhombic_triacontahedron_face_to_face_counts(self) -> None:
        # I-transitivity on RTH faces means every central face gets
        # the same number of face-to-face placements.
        from apeiron.util import load_candidate
        from pathlib import Path
        rth_path = Path(__file__).resolve().parent.parent / "candidates" / "rhombic_triacontahedron.json"
        if not rth_path.exists():
            pytest.skip(f"RTH fixture not found at {rth_path}")
        rth = load_candidate(rth_path)
        counts = [len(face_to_face_placements(rth, i)) for i in range(len(rth.faces))]
        # All faces are in a single I-orbit (transitivity), so all
        # counts are equal.
        assert len(set(counts)) == 1, f"placement counts vary across faces: {counts}"
        # And the count is positive — the RTH does admit face-to-face
        # placements.
        assert counts[0] >= 1


# -- Vertex and Edge feature types ------------------------------------


class TestVertexAndEdge:
    def test_vertex_construction(self) -> None:
        assert Vertex(index=0).index == 0
        assert Vertex(index=31).index == 31

    def test_vertex_rejects_negative(self) -> None:
        with pytest.raises(ValueError, match="non-negative"):
            Vertex(index=-1)

    def test_vertex_equality_and_hash(self) -> None:
        assert Vertex(3) == Vertex(3)
        assert hash(Vertex(3)) == hash(Vertex(3))
        assert Vertex(3) != Vertex(4)

    def test_edge_construction(self) -> None:
        e = Edge(lo=0, hi=5)
        assert e.lo == 0
        assert e.hi == 5

    def test_edge_rejects_unordered(self) -> None:
        with pytest.raises(ValueError, match="lo < hi"):
            Edge(lo=5, hi=0)

    def test_edge_rejects_degenerate(self) -> None:
        with pytest.raises(ValueError, match="lo < hi"):
            Edge(lo=3, hi=3)

    def test_edge_rejects_negative(self) -> None:
        with pytest.raises(ValueError, match="non-negative"):
            Edge(lo=-1, hi=2)

    def test_edge_equality_and_hash(self) -> None:
        assert Edge(1, 4) == Edge(1, 4)
        assert hash(Edge(1, 4)) == hash(Edge(1, 4))
        assert Edge(1, 4) != Edge(2, 4)


# -- _edges_of --------------------------------------------------------


class TestEdgesOf:
    def test_cube_has_12_edges(self) -> None:
        assert len(_edges_of(_unit_cube())) == 12

    def test_tetrahedron_has_6_edges(self) -> None:
        tetra_verts = [
            _v((0, 0), (0, 0), (0, 0)),
            _v((2, 0), (0, 0), (0, 0)),
            _v((0, 0), (2, 0), (0, 0)),
            _v((0, 0), (0, 0), (2, 0)),
        ]
        tetra = Polyhedron.from_vertices(tetra_verts)
        assert len(_edges_of(tetra)) == 6

    def test_edges_are_sorted_pairs(self) -> None:
        for (lo, hi) in _edges_of(_unit_cube()):
            assert lo < hi


# -- incidence_defect --------------------------------------------------


def _face_sharer_at_plus_x() -> PlacedTile:
    """Cube face-sharer: neighbour through the +x face of the
    unit-cube fixture (central vertices at stored (0..2)^3), i.e.
    identity rotation translated by stored (+2, 0, 0) — real +1 shift
    by one cube edge — so the neighbour occupies stored (2..4)^[x]
    × (0..2)^[y, z] and its -x face coincides with the central's +x
    face.
    """
    return PlacedTile(
        translation=_v((2, 0), (0, 0), (0, 0)),
        rotation=Rotation.identity(),
    )


def _edge_sharer_at_plus_x_plus_y() -> PlacedTile:
    """Cube edge-sharer: neighbour shifted by stored (+2, +2, 0) —
    real (+1, +1, 0) — sharing only the edge from stored (2, 2, 0) to
    (2, 2, 2) with the central.
    """
    return PlacedTile(
        translation=_v((2, 0), (2, 0), (0, 0)),
        rotation=Rotation.identity(),
    )


def _vertex_sharer_at_plus_xyz() -> PlacedTile:
    """Cube vertex-sharer: neighbour shifted by stored (+2, +2, +2) —
    real (+1, +1, +1) — sharing only the vertex (2, 2, 2) with the
    central.
    """
    return PlacedTile(
        translation=_v((2, 0), (2, 0), (2, 0)),
        rotation=Rotation.identity(),
    )


class TestIncidenceDefect:
    def test_empty_config_vertex_defect(self) -> None:
        cube = _unit_cube()
        empty = CoronaConfig.from_neighbours(cube, [])
        # Only the central contributes; defect = expected - 1.
        assert incidence_defect(empty, Vertex(0), expected=8) == 7

    def test_empty_config_edge_defect(self) -> None:
        cube = _unit_cube()
        empty = CoronaConfig.from_neighbours(cube, [])
        edges = list(_edges_of(cube))
        lo, hi = edges[0]
        assert incidence_defect(empty, Edge(lo, hi), expected=4) == 3

    def test_face_sharer_reduces_shared_face_defects(self) -> None:
        # The face-sharer across +x touches 4 central vertices and
        # 4 central edges. Defects at those features drop by 1
        # relative to the empty config.
        cube = _unit_cube()
        empty = CoronaConfig.from_neighbours(cube, [])
        with_sharer = CoronaConfig.from_neighbours(
            cube, [_face_sharer_at_plus_x()]
        )
        # Vertices of the central's +x face: (2, *, *) — stored
        # coordinates x = 2.
        shared_vertex_indices = [
            i for i, v in enumerate(cube.vertices) if v.x == ZPhi(2, 0)
        ]
        assert len(shared_vertex_indices) == 4
        for idx in shared_vertex_indices:
            assert incidence_defect(with_sharer, Vertex(idx), expected=8) == (
                incidence_defect(empty, Vertex(idx), expected=8) - 1
            )

    def test_non_shared_vertex_defect_unchanged(self) -> None:
        cube = _unit_cube()
        empty = CoronaConfig.from_neighbours(cube, [])
        with_sharer = CoronaConfig.from_neighbours(
            cube, [_face_sharer_at_plus_x()]
        )
        # Vertices on the central's -x face (x = 0) aren't touched.
        for i, v in enumerate(cube.vertices):
            if v.x == ZPhi(0, 0):
                assert incidence_defect(
                    with_sharer, Vertex(i), expected=8,
                ) == incidence_defect(empty, Vertex(i), expected=8)

    def test_edge_sharer_contributes_to_shared_edge_only(self) -> None:
        # The edge-diagonal neighbour at (+4, +4, 0) shares only the
        # edge (2, 2, 0)-(2, 2, 2) with the central — a single
        # central edge, two central vertices.
        cube = _unit_cube()
        with_edge_sharer = CoronaConfig.from_neighbours(
            cube, [_edge_sharer_at_plus_x_plus_y()]
        )
        # The shared-edge endpoints are the two central vertices with
        # x = 2 and y = 2.
        shared = [
            i for i, v in enumerate(cube.vertices)
            if v.x == ZPhi(2, 0) and v.y == ZPhi(2, 0)
        ]
        assert len(shared) == 2
        lo, hi = min(shared), max(shared)
        # That edge must be in the cube's edge set for this test to
        # be meaningful.
        assert (lo, hi) in _edges_of(cube)
        # Defect at the edge is 4 - 2 (central + edge-sharer) = 2.
        assert incidence_defect(
            with_edge_sharer, Edge(lo, hi), expected=4,
        ) == 2

    def test_rejects_invalid_vertex_index(self) -> None:
        cube = _unit_cube()
        empty = CoronaConfig.from_neighbours(cube, [])
        with pytest.raises(ValueError, match="outside"):
            incidence_defect(empty, Vertex(999), expected=8)

    def test_rejects_non_edge_pair(self) -> None:
        cube = _unit_cube()
        empty = CoronaConfig.from_neighbours(cube, [])
        # Pick two cube vertices that aren't connected by an edge —
        # the face-diagonal pair (0, 0, 0) and (2, 2, 0). These are
        # indices 0 and 6 in the canonical-order cube vertex list
        # (sorted by (x, y, z)).
        with pytest.raises(ValueError, match="not an edge"):
            incidence_defect(empty, Edge(0, 7), expected=4)

    def test_rejects_non_feature_type(self) -> None:
        cube = _unit_cube()
        empty = CoronaConfig.from_neighbours(cube, [])
        with pytest.raises(TypeError, match="Vertex or Edge"):
            incidence_defect(empty, "not a feature", expected=8)   # type: ignore[arg-type]


# -- has_interior_overlap ---------------------------------------------


class TestHasInteriorOverlap:
    def test_empty_config_no_overlap(self) -> None:
        cube = _unit_cube()
        empty = CoronaConfig.from_neighbours(cube, [])
        assert not has_interior_overlap(empty)

    def test_face_sharer_no_overlap(self) -> None:
        cube = _unit_cube()
        c = CoronaConfig.from_neighbours(cube, [_face_sharer_at_plus_x()])
        assert not has_interior_overlap(c)

    def test_edge_sharer_no_overlap(self) -> None:
        cube = _unit_cube()
        c = CoronaConfig.from_neighbours(
            cube, [_edge_sharer_at_plus_x_plus_y()]
        )
        assert not has_interior_overlap(c)

    def test_vertex_sharer_no_overlap(self) -> None:
        cube = _unit_cube()
        c = CoronaConfig.from_neighbours(
            cube, [_vertex_sharer_at_plus_xyz()]
        )
        assert not has_interior_overlap(c)

    def test_coincident_neighbour_overlaps(self) -> None:
        # Neighbour placed exactly on top of central: full volume
        # overlap.
        cube = _unit_cube()
        coincident = PlacedTile(
            translation=_v((0, 0), (0, 0), (0, 0)),
            rotation=Rotation.identity(),
        )
        c = CoronaConfig.from_neighbours(cube, [coincident])
        assert has_interior_overlap(c)

    def test_strictly_nested_neighbour_overlaps(self) -> None:
        # Neighbour shifted by a small amount (stored +1 in x,
        # real +0.5) — its interior overlaps the central's interior
        # along a slab.
        cube = _unit_cube()
        shifted = PlacedTile(
            translation=_v((1, 0), (0, 0), (0, 0)),
            rotation=Rotation.identity(),
        )
        c = CoronaConfig.from_neighbours(cube, [shifted])
        assert has_interior_overlap(c)

    def test_disjoint_neighbour_no_overlap(self) -> None:
        # Far-away neighbour: no chance of overlap.
        cube = _unit_cube()
        far = PlacedTile(
            translation=_v((20, 0), (20, 0), (20, 0)),
            rotation=Rotation.identity(),
        )
        c = CoronaConfig.from_neighbours(cube, [far])
        assert not has_interior_overlap(c)

    def test_two_non_overlapping_neighbours(self) -> None:
        # Two face-sharers on opposite faces; neither overlaps the
        # central nor each other.
        cube = _unit_cube()
        plus_x = _face_sharer_at_plus_x()
        minus_x = PlacedTile(
            translation=_v((-2, 0), (0, 0), (0, 0)),
            rotation=Rotation.identity(),
        )
        c = CoronaConfig.from_neighbours(cube, [plus_x, minus_x])
        assert not has_interior_overlap(c)

    def test_two_overlapping_neighbours(self) -> None:
        # Two neighbours at the same face-sharing position (+x face)
        # but different rotations: the ROT_3 rotation of the cube
        # about the (1, 1, 1) body diagonal is still a cube at the
        # same bounding volume, so the two placed tiles occupy
        # identical real regions and interior-overlap.
        cube = _unit_cube()
        plus_x_a = _face_sharer_at_plus_x()
        plus_x_b = PlacedTile(
            translation=_v((2, 0), (0, 0), (0, 0)),
            rotation=ROT_3,
        )
        c = CoronaConfig.from_neighbours(cube, [plus_x_a, plus_x_b])
        assert has_interior_overlap(c)
