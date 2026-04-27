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

from apeiron.hierarchy import Supertile, expand_one
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
