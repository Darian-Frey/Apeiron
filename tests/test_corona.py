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

from apeiron.corona import CoronaConfig, PlacedTile
from apeiron.polyhedron import Polyhedron
from apeiron.symmetry import ICOSAHEDRAL, ROT_2, ROT_3, Rotation, Vec3
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
