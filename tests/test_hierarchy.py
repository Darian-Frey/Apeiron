"""Tests for apeiron.hierarchy.

Scope for sub-commit A: the ``Supertile`` dataclass and the
``expand_one`` lookup. The full four-pillar apparatus
(recognisability, inflation argument, fourth-pillar protocol)
arrives in sub-commits B–E; this commit ships the structural
backbone the rest will sit on.

Pillar tagging (decorator + sweep) lands in sub-commit E. For now,
pillar-1 obligations are noted by docstring per Claude (web)'s
2026-04-23 directive, with the actual enforcement (``is_primitive``
+ ``perron_frobenius_in_zphi`` returning a ZPhi rather than None) in
``test_substitution.py``.
"""

from __future__ import annotations

import pytest

from apeiron.hierarchy import (
    IndistinguishablePair,
    PatchTile,
    RecognisabilityResult,
    Supertile,
    TilePatch,
    expand_one,
    is_recognisable,
    neighbourhood_signature,
)
from apeiron.substitution import PositionedTile, SubstitutionRule
from apeiron.symmetry import Mat3, Rotation, Vec3
from apeiron.zphi import ZPhi


def _origin() -> Vec3:
    return Vec3(ZPhi(0, 0), ZPhi(0, 0), ZPhi(0, 0))


def _dummy_tile(idx: int) -> PositionedTile:
    """A PositionedTile at the origin with identity rotation."""
    return PositionedTile(
        prototile_index=idx,
        translation=_origin(),
        rotation=Rotation.identity(),
    )


# -- Supertile construction and validation ----------------------------


class TestSupertile:
    def test_level_zero_has_no_children(self) -> None:
        s = Supertile(level=0, prototile_index=0)
        assert s.level == 0
        assert s.prototile_index == 0
        assert s.children == ()

    def test_level_one_with_children(self) -> None:
        children = (_dummy_tile(0), _dummy_tile(1))
        s = Supertile(level=1, prototile_index=0, children=children)
        assert s.level == 1
        assert s.children == children

    def test_rejects_negative_level(self) -> None:
        with pytest.raises(ValueError, match="level must be"):
            Supertile(level=-1, prototile_index=0)

    def test_rejects_negative_prototile_index(self) -> None:
        with pytest.raises(ValueError, match="prototile_index"):
            Supertile(level=0, prototile_index=-1)

    def test_level_zero_rejects_children(self) -> None:
        with pytest.raises(ValueError, match="level == 0"):
            Supertile(
                level=0,
                prototile_index=0,
                children=(_dummy_tile(0),),
            )

    def test_equality_and_hash(self) -> None:
        children = (_dummy_tile(0),)
        s1 = Supertile(level=1, prototile_index=0, children=children)
        s2 = Supertile(level=1, prototile_index=0, children=children)
        assert s1 == s2
        assert hash(s1) == hash(s2)

    def test_distinct_levels_unequal(self) -> None:
        children = (_dummy_tile(0),)
        assert Supertile(0, 0) != Supertile(1, 0, children=children)

    def test_hashable_in_set(self) -> None:
        s = Supertile(level=0, prototile_index=0)
        assert {s, s} == {s}


# -- expand_one lookup ------------------------------------------------


def _two_prototile_rule() -> SubstitutionRule:
    """A trivial 2-prototile rule for testing expand_one.

    σ(0) = (0, 1)   — type-0 dissects into one type-0 + one type-1
    σ(1) = (0,)     — type-1 dissects into a single type-0
    """
    return SubstitutionRule(
        n_prototiles=2,
        inflation=Mat3.identity(),
        dissections=(
            (_dummy_tile(0), _dummy_tile(1)),
            (_dummy_tile(0),),
        ),
    )


class TestExpandOne:
    def test_expand_one_returns_dissection(self) -> None:
        rule = _two_prototile_rule()
        children_0 = expand_one(rule, 0)
        assert len(children_0) == 2
        assert children_0[0].prototile_index == 0
        assert children_0[1].prototile_index == 1

    def test_expand_one_for_second_prototile(self) -> None:
        rule = _two_prototile_rule()
        children_1 = expand_one(rule, 1)
        assert len(children_1) == 1
        assert children_1[0].prototile_index == 0

    def test_expand_one_rejects_out_of_range(self) -> None:
        rule = _two_prototile_rule()
        with pytest.raises(ValueError, match="outside"):
            expand_one(rule, 2)
        with pytest.raises(ValueError, match="outside"):
            expand_one(rule, -1)

    def test_expand_one_is_pure(self) -> None:
        # Same rule + index → same children every time. (Pure lookup;
        # not yet memoised, but the contract is invariance.)
        rule = _two_prototile_rule()
        a = expand_one(rule, 0)
        b = expand_one(rule, 0)
        assert a == b


# -- Recognisability (pillar 2) ---------------------------------------


def _grid_neighbour_within(
    positions: tuple[tuple[int, int], ...],
) -> "callable":
    """Build a neighbour-within oracle from a set of integer 2D
    positions, using Chebyshev (L∞) distance — simple grid-step
    radius. Returns a function ``(i, j, r) -> bool``.
    """
    def fn(i: int, j: int, r: int) -> bool:
        xi, yi = positions[i]
        xj, yj = positions[j]
        return max(abs(xi - xj), abs(yi - yj)) <= r
    return fn


def _ammann_beenker_like_patch() -> TilePatch:
    """Synthetic Ammann–Beenker–style fixture for the pillar-2 oracle.

    A 4-tile patch that mimics the structural property captured in
    Baake & Grimm Vol. 1 §6: each tile's *immediate* neighbour
    multiset uniquely determines its supertile parent. Tile types
    chosen so that:

    - Squares (type 0) at positions (0, 0) and (3, 0) belong to
      supertile 0; each is adjacent to one rhombus and to the other
      square.

    Wait, that conflicts — let me redo. We want tiles with the same
    *parent* to share a signature, and tiles with *different*
    parents to have distinct signatures. Concrete configuration:

    - Tile 0: type=0 (square), parent=0, position=(0, 0).
      Neighbours within radius 1: tile 1 (type 0). Signature: (0, 0).
    - Tile 1: type=0 (square), parent=0, position=(1, 0).
      Neighbours within radius 1: tile 0 (type 0), tile 2 (type 1).
      Signature: (0, 0, 1).
    - Tile 2: type=1 (rhombus), parent=1, position=(2, 0).
      Neighbours within radius 1: tile 1 (type 0), tile 3 (type 1).
      Signature: (0, 1, 1).
    - Tile 3: type=1 (rhombus), parent=1, position=(3, 0).
      Neighbours within radius 1: tile 2 (type 1). Signature: (1, 1).

    All four signatures are distinct → recognisable at radius 1.
    """
    positions = ((0, 0), (1, 0), (2, 0), (3, 0))
    tiles = (
        PatchTile(tile_type=0, parent_supertile=0),
        PatchTile(tile_type=0, parent_supertile=0),
        PatchTile(tile_type=1, parent_supertile=1),
        PatchTile(tile_type=1, parent_supertile=1),
    )
    return TilePatch(
        tiles=tiles,
        neighbour_within=_grid_neighbour_within(positions),
    )


def _non_recognisable_patch() -> TilePatch:
    """Two tiles with identical neighbourhood signatures but different
    parents. radius=1 sees both with signature (0, 0); but tile 0 is
    in supertile 0 and tile 1 is in supertile 1. No radius will
    distinguish them in this synthetic patch (they're symmetric),
    so ``is_recognisable`` should always fail and surface an
    ``IndistinguishablePair``.
    """
    positions = ((0, 0), (1, 0))
    tiles = (
        PatchTile(tile_type=0, parent_supertile=0),
        PatchTile(tile_type=0, parent_supertile=1),
    )
    return TilePatch(
        tiles=tiles,
        neighbour_within=_grid_neighbour_within(positions),
    )


class TestNeighbourhoodSignature:
    def test_includes_centre_tile(self) -> None:
        patch = _ammann_beenker_like_patch()
        sig = neighbourhood_signature(patch, 0, radius=1)
        # Centre tile is type 0; one neighbour at distance 1 is type 0.
        assert sig == (0, 0)

    def test_signature_is_sorted(self) -> None:
        patch = _ammann_beenker_like_patch()
        sig = neighbourhood_signature(patch, 1, radius=1)
        assert sig == tuple(sorted(sig))

    def test_radius_zero_yields_only_centre(self) -> None:
        patch = _ammann_beenker_like_patch()
        sig = neighbourhood_signature(patch, 0, radius=0)
        # Radius 0 means "tile and tiles at distance 0" — only the
        # tile itself (since other tiles are at distance ≥ 1).
        assert sig == (0,)

    def test_larger_radius_includes_more(self) -> None:
        patch = _ammann_beenker_like_patch()
        # At radius 3, tile 0 sees all four tiles.
        sig = neighbourhood_signature(patch, 0, radius=3)
        assert sig == (0, 0, 1, 1)

    def test_rejects_out_of_range_index(self) -> None:
        patch = _ammann_beenker_like_patch()
        with pytest.raises(IndexError):
            neighbourhood_signature(patch, 99, radius=1)


class TestIsRecognisable:
    def test_ammann_beenker_recognisable_at_radius_1(self) -> None:
        patch = _ammann_beenker_like_patch()
        result = is_recognisable(patch, max_radius=5)
        assert result.is_recognisable is True
        assert result.radius_used == 1
        assert result.radius_cap_reached is False
        assert isinstance(result.witness, dict)

    def test_witness_maps_signatures_to_parents(self) -> None:
        patch = _ammann_beenker_like_patch()
        result = is_recognisable(patch, max_radius=5)
        assert isinstance(result.witness, dict)
        # Every patch tile's signature should appear in the witness
        # and map to its parent.
        for i, tile in enumerate(patch.tiles):
            sig = neighbourhood_signature(patch, i, result.radius_used)
            assert result.witness[sig] == tile.parent_supertile

    def test_non_recognisable_returns_pair(self) -> None:
        patch = _non_recognisable_patch()
        result = is_recognisable(patch, max_radius=3)
        assert result.is_recognisable is False
        assert result.radius_cap_reached is True
        assert isinstance(result.witness, IndistinguishablePair)
        # The pair must be tiles 0 and 1 with parents 0 and 1.
        pair = result.witness
        assert {pair.tile_a, pair.tile_b} == {0, 1}
        assert {pair.parent_a, pair.parent_b} == {0, 1}

    def test_indistinguishable_pair_records_radius(self) -> None:
        patch = _non_recognisable_patch()
        result = is_recognisable(patch, max_radius=3)
        pair = result.witness
        assert isinstance(pair, IndistinguishablePair)
        # The radius recorded is the cap (3) since failure persisted.
        assert pair.radius == 3

    def test_rejects_zero_max_radius(self) -> None:
        patch = _ammann_beenker_like_patch()
        with pytest.raises(ValueError, match="max_radius"):
            is_recognisable(patch, max_radius=0)

    def test_search_terminates_at_smallest_sufficient_radius(
        self,
    ) -> None:
        # Pick a patch where radius 1 works.
        patch = _ammann_beenker_like_patch()
        result = is_recognisable(patch, max_radius=10)
        # Should stop at radius 1, not exhaust the cap.
        assert result.radius_used == 1

    def test_returns_recognisability_result_type(self) -> None:
        patch = _ammann_beenker_like_patch()
        result = is_recognisable(patch, max_radius=2)
        assert isinstance(result, RecognisabilityResult)
