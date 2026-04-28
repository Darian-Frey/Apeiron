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

from apeiron.corona import CoronaConfig
from apeiron.hierarchy import (
    FourthPillarArgument,
    HierarchicalCounterexample,
    HierarchicalWitness,
    IndistinguishablePair,
    InflationArgument,
    InflationFailure,
    PatchTile,
    PlacedSubtile,
    RecognisabilityResult,
    Supertile,
    TilePatch,
    expand_one,
    expand_supertile,
    expand_supertile_with_parents,
    inflation_argument,
    is_recognisable,
    make_euclidean_squared_oracle,
    neighbourhood_signature,
    patch_from_supertile,
)
from apeiron.polyhedron import Polyhedron
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


# -- Inflation argument (pillar 3) ------------------------------------


def _penrose_p3_rule() -> SubstitutionRule:
    """Penrose P3 thick/thin rhombus rule used as the pillar-3 oracle.

    Substitution matrix [[2, 1], [1, 1]]: thick → 2 thick + 1 thin;
    thin → 1 thick + 1 thin. Characteristic polynomial x² − 3x + 1;
    PF eigenvalue (3 + √5)/2 = 1 + φ = φ² ≈ 2.618. Primitive (M
    itself is all-positive). Recognisable in 2D at radius 1 in the
    standard fixture (Baake & Grimm Vol. 1 §4); for the algebraic
    pillar-3 test we supply a synthetic ``RecognisabilityResult``
    directly rather than building a P3 patch — pillar 3 only needs
    the *outcome* of pillar 2, not the geometric machinery that
    produced it.
    """
    return SubstitutionRule(
        n_prototiles=2,
        inflation=Mat3.identity(),
        dissections=(
            (_dummy_tile(0), _dummy_tile(0), _dummy_tile(1)),
            (_dummy_tile(0), _dummy_tile(1)),
        ),
    )


def _trivial_one_prototile_rule() -> SubstitutionRule:
    """A 1-prototile rule with σ(0) = (0,) — substitution matrix
    [[1]], PF eigenvalue 1. Primitive but degenerate; the inflation
    argument's λ > 1 condition fails.
    """
    return SubstitutionRule(
        n_prototiles=1,
        inflation=Mat3.identity(),
        dissections=((_dummy_tile(0),),),
    )


def _non_primitive_rule() -> SubstitutionRule:
    """A reducible rule: σ(0) = (0,), σ(1) = (1,). Substitution
    matrix is the 2×2 identity; not primitive.
    """
    return SubstitutionRule(
        n_prototiles=2,
        inflation=Mat3.identity(),
        dissections=((_dummy_tile(0),), (_dummy_tile(1),)),
    )


def _success_recognisability(radius: int = 1) -> RecognisabilityResult:
    """Synthetic successful recognisability witness.

    Pillar 3's input is the *outcome* of pillar 2 — a structurally
    valid ``RecognisabilityResult`` with ``is_recognisable=True``.
    The witness map's contents don't affect pillar 3's reasoning;
    only the boolean and the radius do.
    """
    return RecognisabilityResult(
        is_recognisable=True,
        radius_used=radius,
        radius_cap_reached=False,
        witness={(0,): 0, (1,): 1},
    )


def _failure_recognisability() -> RecognisabilityResult:
    """Synthetic failing recognisability witness, for the
    not-recognisable branch of pillar 3.
    """
    return RecognisabilityResult(
        is_recognisable=False,
        radius_used=5,
        radius_cap_reached=True,
        witness=IndistinguishablePair(
            tile_a=0, tile_b=1, radius=5, parent_a=0, parent_b=1,
        ),
    )


class TestInflationArgumentP3:
    """Penrose P3 as the primary pillar-3 oracle."""

    def test_p3_yields_inflation_argument(self) -> None:
        result = inflation_argument(
            _penrose_p3_rule(), _success_recognisability(),
        )
        assert isinstance(result, InflationArgument)

    def test_witness_carries_phi_squared_eigenvalue(self) -> None:
        result = inflation_argument(
            _penrose_p3_rule(), _success_recognisability(),
        )
        assert isinstance(result, InflationArgument)
        # P3's PF eigenvalue is φ² = 1 + φ = ZPhi(1, 1).
        assert result.pf_eigenvalue == ZPhi(1, 1)

    def test_witness_records_recognisability_radius(self) -> None:
        result = inflation_argument(
            _penrose_p3_rule(), _success_recognisability(radius=2),
        )
        assert isinstance(result, InflationArgument)
        assert result.recognisability_radius == 2


class TestInflationArgumentFailures:
    def test_non_primitive_rule_fails(self) -> None:
        result = inflation_argument(
            _non_primitive_rule(), _success_recognisability(),
        )
        assert isinstance(result, InflationFailure)
        assert result.reason == "not primitive"

    def test_pf_eq_one_rule_fails(self) -> None:
        # Trivial 1×1 rule with PF = 1.
        result = inflation_argument(
            _trivial_one_prototile_rule(), _success_recognisability(),
        )
        assert isinstance(result, InflationFailure)
        assert result.reason == "pf <= 1"

    def test_non_recognisable_input_fails(self) -> None:
        result = inflation_argument(
            _penrose_p3_rule(), _failure_recognisability(),
        )
        assert isinstance(result, InflationFailure)
        assert result.reason == "not recognisable"

    def test_failure_carries_detail_string(self) -> None:
        result = inflation_argument(
            _non_primitive_rule(), _success_recognisability(),
        )
        assert isinstance(result, InflationFailure)
        assert result.detail
        assert "pillar 1" in result.detail.lower() or "primitive" in result.detail.lower()


class TestInflationArgumentStructure:
    def test_argument_is_frozen(self) -> None:
        arg = InflationArgument(
            pf_eigenvalue=ZPhi(1, 1), recognisability_radius=1,
        )
        with pytest.raises(Exception):
            arg.recognisability_radius = 2  # type: ignore[misc]

    def test_argument_is_hashable(self) -> None:
        a = InflationArgument(
            pf_eigenvalue=ZPhi(1, 1), recognisability_radius=1,
        )
        b = InflationArgument(
            pf_eigenvalue=ZPhi(1, 1), recognisability_radius=1,
        )
        assert hash(a) == hash(b)
        assert {a, b} == {a}

    def test_failure_is_frozen(self) -> None:
        f = InflationFailure(reason="x", detail="y")
        with pytest.raises(Exception):
            f.reason = "z"  # type: ignore[misc]


# -- Fourth pillar framework (pillar 4) -------------------------------


def _v(x: tuple[int, int], y: tuple[int, int], z: tuple[int, int]) -> Vec3:
    return Vec3(ZPhi(*x), ZPhi(*y), ZPhi(*z))


def _trivial_corona_config() -> CoronaConfig:
    """Smallest valid CoronaConfig: a tetrahedron central with no
    neighbours. Used as a stand-in CoronaConfig for protocol-level
    tests; the geometric content doesn't affect pillar-4 *framework*
    tests, since this commit ships no concrete pillar-4 logic.
    """
    tetra_verts = [
        _v((0, 0), (0, 0), (0, 0)),
        _v((2, 0), (0, 0), (0, 0)),
        _v((0, 0), (2, 0), (0, 0)),
        _v((0, 0), (0, 0), (2, 0)),
    ]
    tetra = Polyhedron.from_vertices(tetra_verts)
    return CoronaConfig.from_neighbours(tetra, [])


class TestHierarchicalWitness:
    def test_construction(self) -> None:
        w = HierarchicalWitness(
            supertile_level=2,
            parent_prototile_index=0,
            rationale="case 1: configuration matches branch A of σ²(P_0)",
        )
        assert w.supertile_level == 2
        assert w.parent_prototile_index == 0
        assert "branch" in w.rationale

    def test_is_frozen(self) -> None:
        w = HierarchicalWitness(
            supertile_level=1, parent_prototile_index=0, rationale="x",
        )
        with pytest.raises(Exception):
            w.supertile_level = 2  # type: ignore[misc]

    def test_is_hashable(self) -> None:
        a = HierarchicalWitness(
            supertile_level=1, parent_prototile_index=0, rationale="x",
        )
        b = HierarchicalWitness(
            supertile_level=1, parent_prototile_index=0, rationale="x",
        )
        assert hash(a) == hash(b)
        assert {a, b} == {a}


class TestHierarchicalCounterexample:
    def test_construction(self) -> None:
        c = HierarchicalCounterexample(
            rationale="no σⁿ embedding for n ≤ 5: edge mismatch on face 3",
            search_bound=5,
        )
        assert c.search_bound == 5
        assert "edge mismatch" in c.rationale

    def test_is_frozen(self) -> None:
        c = HierarchicalCounterexample(rationale="x", search_bound=3)
        with pytest.raises(Exception):
            c.search_bound = 5  # type: ignore[misc]


class TestFourthPillarArgumentProtocol:
    """Protocol-satisfaction tests.

    The framework only ships the protocol; concrete implementations
    live per-candidate at ``candidates/<name>/fourth_pillar.py``.
    These tests verify the contract: a class implementing the two
    methods with the right signatures satisfies the protocol via
    ``isinstance`` (because the protocol is ``@runtime_checkable``).
    """

    def test_minimal_implementation_satisfies_protocol(self) -> None:
        class _Always(FourthPillarArgument):
            def local_configurations(self) -> frozenset[CoronaConfig]:
                return frozenset()

            def verify_hierarchical(
                self, _config: CoronaConfig,
            ) -> HierarchicalWitness | HierarchicalCounterexample:
                return HierarchicalWitness(
                    supertile_level=1,
                    parent_prototile_index=0,
                    rationale="trivial",
                )

        impl = _Always()
        assert isinstance(impl, FourthPillarArgument)

    def test_missing_method_fails_protocol_check(self) -> None:
        class _Incomplete:
            def local_configurations(self) -> frozenset[CoronaConfig]:
                return frozenset()
            # No verify_hierarchical defined.

        impl = _Incomplete()
        assert not isinstance(impl, FourthPillarArgument)

    def test_concrete_implementation_round_trip(self) -> None:
        """A concrete implementation can be queried; the witness
        type returned matches the protocol.
        """
        config = _trivial_corona_config()

        class _OneConfigYesPillar4(FourthPillarArgument):
            def local_configurations(self) -> frozenset[CoronaConfig]:
                return frozenset({config})

            def verify_hierarchical(
                self, c: CoronaConfig,
            ) -> HierarchicalWitness | HierarchicalCounterexample:
                if c == config:
                    return HierarchicalWitness(
                        supertile_level=1,
                        parent_prototile_index=0,
                        rationale="single-config trivial witness",
                    )
                return HierarchicalCounterexample(
                    rationale="config not in the local set",
                    search_bound=1,
                )

        impl = _OneConfigYesPillar4()
        assert isinstance(impl, FourthPillarArgument)
        configs = impl.local_configurations()
        assert len(configs) == 1
        verdict = impl.verify_hierarchical(next(iter(configs)))
        assert isinstance(verdict, HierarchicalWitness)

    def test_concrete_implementation_can_return_counterexample(self) -> None:
        class _AlwaysCounterexample(FourthPillarArgument):
            def local_configurations(self) -> frozenset[CoronaConfig]:
                return frozenset()

            def verify_hierarchical(
                self, _c: CoronaConfig,
            ) -> HierarchicalWitness | HierarchicalCounterexample:
                return HierarchicalCounterexample(
                    rationale="placeholder negative",
                    search_bound=2,
                )

        impl = _AlwaysCounterexample()
        verdict = impl.verify_hierarchical(_trivial_corona_config())
        assert isinstance(verdict, HierarchicalCounterexample)


# -- expand_supertile (level-N supertile flattening) ------------------


def _phi_inflation_mat() -> Mat3:
    """φ·I, the canonical Danzer linear inflation."""
    phi = ZPhi(0, 1)
    z = ZPhi(0, 0)
    return Mat3(
        Vec3(phi, z, z),
        Vec3(z, phi, z),
        Vec3(z, z, phi),
    )


def _placeholder_dissection_for(matrix_column: list[int]) -> tuple[
    PositionedTile, ...,
]:
    """Build a tuple of PositionedTile (origin / identity) whose
    prototile_index multiset matches ``matrix_column``. Used for
    expand_supertile tests where the geometric content of children
    doesn't matter — only the type counts.
    """
    children: list[PositionedTile] = []
    for type_idx, count in enumerate(matrix_column):
        for _ in range(count):
            children.append(
                PositionedTile(
                    prototile_index=type_idx,
                    translation=_origin(),
                    rotation=Rotation.identity(),
                )
            )
    return tuple(children)


def _danzer_placeholder_rule() -> SubstitutionRule:
    """Frettlöh's Danzer 4×4 substitution matrix as a SubstitutionRule
    with placeholder dissection geometry. Identical in spirit to the
    helper in tests/integration/test_danzer_abck.py; duplicated here
    rather than imported so this test file stays unit-scoped (it
    only needs the *combinatorial* structure of σ for
    expand_supertile coverage).
    """
    matrix = [
        [0, 0, 1, 0],
        [3, 2, 0, 1],
        [2, 1, 2, 0],
        [6, 4, 2, 1],
    ]
    n = 4
    dissections = tuple(
        _placeholder_dissection_for([matrix[t][col] for t in range(n)])
        for col in range(n)
    )
    return SubstitutionRule(
        n_prototiles=n,
        inflation=_phi_inflation_mat(),
        dissections=dissections,
    )


class TestExpandSupertile:
    def test_level_zero_returns_single_leaf(self) -> None:
        rule = _two_prototile_rule()
        leaves = expand_supertile(rule, prototile_index=0, level=0)
        assert len(leaves) == 1
        assert isinstance(leaves[0], PlacedSubtile)
        assert leaves[0].prototile_index == 0
        assert leaves[0].rotation == Rotation.identity()

    def test_level_one_matches_expand_one(self) -> None:
        # At level 1 the leaf count + types should agree with
        # expand_one's direct dissection (only the position/rotation
        # accumulation in expand_supertile is novel).
        rule = _two_prototile_rule()
        leaves = expand_supertile(rule, prototile_index=0, level=1)
        children = expand_one(rule, 0)
        assert len(leaves) == len(children)
        leaf_types = sorted(leaf.prototile_index for leaf in leaves)
        child_types = sorted(c.prototile_index for c in children)
        assert leaf_types == child_types

    def test_rejects_negative_level(self) -> None:
        rule = _two_prototile_rule()
        with pytest.raises(ValueError, match="level must be"):
            expand_supertile(rule, 0, -1)

    def test_rejects_out_of_range_prototile_index(self) -> None:
        rule = _two_prototile_rule()
        with pytest.raises(ValueError, match="outside"):
            expand_supertile(rule, 99, 1)

    def test_count_matches_substitution_matrix_power(self) -> None:
        # For Danzer's 4×4 matrix, leaf counts at level N equal the
        # column sums of M^N. Cross-checks expand_supertile against
        # an independent sanity check (numpy matrix power).
        import numpy as np
        rule = _danzer_placeholder_rule()
        M = np.array([
            [0, 0, 1, 0],
            [3, 2, 0, 1],
            [2, 1, 2, 0],
            [6, 4, 2, 1],
        ])
        for level in (1, 2, 3):
            for col in range(4):
                expected = int(np.linalg.matrix_power(M, level)[:, col].sum())
                leaves = expand_supertile(rule, col, level)
                assert len(leaves) == expected, (
                    f"level={level}, col={col}: expected {expected} "
                    f"leaves, got {len(leaves)}"
                )

    def test_leaf_type_distribution_matches_matrix_column(self) -> None:
        # Per-prototile counts at level N match column N of M^level.
        import numpy as np
        rule = _danzer_placeholder_rule()
        M = np.array([
            [0, 0, 1, 0],
            [3, 2, 0, 1],
            [2, 1, 2, 0],
            [6, 4, 2, 1],
        ])
        from collections import Counter
        for level in (1, 2):
            for col in range(4):
                leaves = expand_supertile(rule, col, level)
                actual = Counter(leaf.prototile_index for leaf in leaves)
                expected_col = np.linalg.matrix_power(M, level)[:, col]
                for type_idx in range(4):
                    assert actual[type_idx] == int(expected_col[type_idx]), (
                        f"level={level}, col={col}, type={type_idx}: "
                        f"expected {expected_col[type_idx]}, got "
                        f"{actual[type_idx]}"
                    )

    def test_returns_tuple_not_list(self) -> None:
        # Frozen output: callers can hash / share without copying.
        rule = _two_prototile_rule()
        leaves = expand_supertile(rule, 0, 1)
        assert isinstance(leaves, tuple)


class TestExpandSupertileWithParents:
    """Verify the level-1 parent tagging used for pillar-2 patch
    construction: each leaf in σⁿ(P) carries the index ``i ∈ [0, k)``
    of the level-1 child of ``P`` from which it descends.
    """

    def test_level_zero_single_leaf_parent_zero(self) -> None:
        # Degenerate: no level-1 dissection has been applied.
        rule = _two_prototile_rule()
        tagged = expand_supertile_with_parents(rule, 0, 0)
        assert len(tagged) == 1
        leaf, parent_id = tagged[0]
        assert leaf.prototile_index == 0
        assert parent_id == 0

    def test_level_one_parent_ids_are_dissection_positions(self) -> None:
        # At level 1 each leaf is the i-th child, so parent_id == i.
        rule = _two_prototile_rule()
        tagged = expand_supertile_with_parents(rule, 0, 1)
        children = expand_one(rule, 0)
        assert len(tagged) == len(children)
        for i, (leaf, parent_id) in enumerate(tagged):
            assert parent_id == i
            assert leaf.prototile_index == children[i].prototile_index

    def test_leaf_component_matches_expand_supertile(self) -> None:
        # The leaf-only projection must be equal element-wise (same
        # geometry, same iteration order) to expand_supertile's output.
        rule = _danzer_placeholder_rule()
        for level in (0, 1, 2, 3):
            plain = expand_supertile(rule, 0, level)
            tagged = expand_supertile_with_parents(rule, 0, level)
            assert len(plain) == len(tagged)
            for a, (b, _parent) in zip(plain, tagged, strict=True):
                assert a == b

    def test_rejects_negative_level(self) -> None:
        rule = _two_prototile_rule()
        with pytest.raises(ValueError, match="level must be"):
            expand_supertile_with_parents(rule, 0, -1)

    def test_rejects_out_of_range_prototile_index(self) -> None:
        rule = _two_prototile_rule()
        with pytest.raises(ValueError, match="outside"):
            expand_supertile_with_parents(rule, 99, 1)

    def test_parent_ids_partition_leaves_correctly(self) -> None:
        # For each parent id i, the leaves tagged i are exactly the
        # level-(N-1) expansion of the i-th level-1 child of P. So the
        # count of leaves with parent_id=i equals the column-sum of
        # M^(N-1) at column = type of that child.
        import numpy as np
        from collections import Counter
        rule = _danzer_placeholder_rule()
        M = np.array([
            [0, 0, 1, 0],
            [3, 2, 0, 1],
            [2, 1, 2, 0],
            [6, 4, 2, 1],
        ])
        for level in (2, 3):
            for col in range(4):
                tagged = expand_supertile_with_parents(rule, col, level)
                groups = Counter(parent_id for _leaf, parent_id in tagged)
                children = rule.dissections[col]
                M_pow = np.linalg.matrix_power(M, level - 1)
                for i, child in enumerate(children):
                    expected = int(M_pow[:, child.prototile_index].sum())
                    assert groups[i] == expected, (
                        f"level={level}, col={col}, parent_id={i}: "
                        f"expected {expected}, got {groups[i]}"
                    )
                # All parent_ids are within [0, k).
                for parent_id in groups:
                    assert 0 <= parent_id < len(children)

    def test_returns_tuple_of_tuples(self) -> None:
        # Frozen output, hashable shape; the inner pairs are tuples too.
        rule = _two_prototile_rule()
        tagged = expand_supertile_with_parents(rule, 0, 1)
        assert isinstance(tagged, tuple)
        for entry in tagged:
            assert isinstance(entry, tuple)
            assert len(entry) == 2
            assert isinstance(entry[0], PlacedSubtile)
            assert isinstance(entry[1], int)


# -- TilePatch from level-N expansion (pillar-2 bridge) ---------------


class TestEuclideanSquaredOracle:
    """Smoke tests for the squared-Euclidean neighbour oracle factory.
    Distances are computed in ZPhi (no float fallback)."""

    def test_self_neighbour_at_radius_zero(self) -> None:
        positions = (
            Vec3(ZPhi(0, 0), ZPhi(0, 0), ZPhi(0, 0)),
            Vec3(ZPhi(2, 0), ZPhi(0, 0), ZPhi(0, 0)),
        )
        oracle = make_euclidean_squared_oracle(positions, ZPhi(1, 0))
        assert oracle(0, 0, 0) is True
        assert oracle(1, 1, 0) is True

    def test_distinct_points_outside_radius_zero(self) -> None:
        positions = (
            Vec3(ZPhi(0, 0), ZPhi(0, 0), ZPhi(0, 0)),
            Vec3(ZPhi(2, 0), ZPhi(0, 0), ZPhi(0, 0)),
        )
        # |pos_0 - pos_1|² = 4; threshold = 0² · 1 = 0 → not within.
        oracle = make_euclidean_squared_oracle(positions, ZPhi(1, 0))
        assert oracle(0, 1, 0) is False

    def test_threshold_inclusive(self) -> None:
        # |pos_0 - pos_1|² = 4 = 1² · 4, so r=1 with step=4 should hit.
        positions = (
            Vec3(ZPhi(0, 0), ZPhi(0, 0), ZPhi(0, 0)),
            Vec3(ZPhi(2, 0), ZPhi(0, 0), ZPhi(0, 0)),
        )
        oracle = make_euclidean_squared_oracle(positions, ZPhi(4, 0))
        assert oracle(0, 1, 1) is True

    def test_radius_grows_quadratically(self) -> None:
        # |pos_0 - pos_1|² = 16; threshold = r² · 1 = r².
        # r=3 → threshold 9 < 16, no; r=4 → threshold 16 = 16, yes.
        positions = (
            Vec3(ZPhi(0, 0), ZPhi(0, 0), ZPhi(0, 0)),
            Vec3(ZPhi(4, 0), ZPhi(0, 0), ZPhi(0, 0)),
        )
        oracle = make_euclidean_squared_oracle(positions, ZPhi(1, 0))
        assert oracle(0, 1, 3) is False
        assert oracle(0, 1, 4) is True

    def test_phi_step_irrational_threshold(self) -> None:
        # Step = phi, so threshold at r=2 is 4*phi.
        # |pos_0 - pos_1|² = 4 (rational). 4 < 4*phi (since phi > 1)? Yes.
        # So r=2 should include them.
        positions = (
            Vec3(ZPhi(0, 0), ZPhi(0, 0), ZPhi(0, 0)),
            Vec3(ZPhi(2, 0), ZPhi(0, 0), ZPhi(0, 0)),
        )
        oracle = make_euclidean_squared_oracle(positions, ZPhi(0, 1))
        assert oracle(0, 1, 2) is True

    def test_rejects_negative_radius(self) -> None:
        positions = (Vec3(ZPhi(0, 0), ZPhi(0, 0), ZPhi(0, 0)),)
        oracle = make_euclidean_squared_oracle(positions, ZPhi(1, 0))
        with pytest.raises(ValueError, match="radius must be"):
            oracle(0, 0, -1)

    def test_rejects_out_of_range_index(self) -> None:
        positions = (Vec3(ZPhi(0, 0), ZPhi(0, 0), ZPhi(0, 0)),)
        oracle = make_euclidean_squared_oracle(positions, ZPhi(1, 0))
        with pytest.raises(IndexError):
            oracle(0, 5, 1)
        with pytest.raises(IndexError):
            oracle(-1, 0, 1)

    def test_position_snapshot_is_immutable(self) -> None:
        # Mutating the source list must not change oracle behaviour.
        positions: list[Vec3] = [
            Vec3(ZPhi(0, 0), ZPhi(0, 0), ZPhi(0, 0)),
            Vec3(ZPhi(2, 0), ZPhi(0, 0), ZPhi(0, 0)),
        ]
        oracle = make_euclidean_squared_oracle(positions, ZPhi(4, 0))
        positions[1] = Vec3(ZPhi(100, 0), ZPhi(0, 0), ZPhi(0, 0))
        # Oracle still uses the original snapshot: r=1 still hits.
        assert oracle(0, 1, 1) is True


class TestPatchFromSupertile:
    """Verify the bridge between expand_supertile_with_parents and
    pillar-2's TilePatch contract.
    """

    def test_level_zero_single_tile_patch(self) -> None:
        rule = _two_prototile_rule()
        patch = patch_from_supertile(
            rule, 0, level=0, radius_step_squared=ZPhi(1, 0),
        )
        assert isinstance(patch, TilePatch)
        assert len(patch.tiles) == 1
        tile = patch.tiles[0]
        assert tile.tile_type == 0
        assert tile.parent_supertile == 0

    def test_level_one_tile_count_and_types(self) -> None:
        rule = _two_prototile_rule()
        patch = patch_from_supertile(
            rule, 0, level=1, radius_step_squared=ZPhi(1, 0),
        )
        children = expand_one(rule, 0)
        assert len(patch.tiles) == len(children)
        for i, (tile, child) in enumerate(
            zip(patch.tiles, children, strict=True)
        ):
            assert tile.tile_type == child.prototile_index
            assert tile.parent_supertile == i

    def test_parent_supertile_matches_tagging(self) -> None:
        # Every tile's parent_supertile equals the level-1 ancestor
        # index returned by expand_supertile_with_parents.
        rule = _danzer_placeholder_rule()
        patch = patch_from_supertile(
            rule, 0, level=2, radius_step_squared=ZPhi(1, 0),
        )
        tagged = expand_supertile_with_parents(rule, 0, 2)
        assert len(patch.tiles) == len(tagged)
        for tile, (_leaf, parent_id) in zip(
            patch.tiles, tagged, strict=True,
        ):
            assert tile.parent_supertile == parent_id

    def test_neighbour_oracle_is_callable(self) -> None:
        rule = _two_prototile_rule()
        patch = patch_from_supertile(
            rule, 0, level=1, radius_step_squared=ZPhi(1, 0),
        )
        # Self-distance is zero; r=0 should hit.
        assert patch.neighbour_within(0, 0, 0) is True

    def test_oracle_uses_leaf_positions(self) -> None:
        # Pick a level-1 patch; verify the oracle's notion of
        # "within radius r" actually reflects leaf positions.
        rule = _danzer_placeholder_rule()
        # All placeholder children are at the origin / identity, so all
        # leaf positions coincide and any radius ≥ 0 should connect them.
        patch = patch_from_supertile(
            rule, 0, level=1, radius_step_squared=ZPhi(1, 0),
        )
        n = len(patch.tiles)
        for i in range(n):
            for j in range(n):
                assert patch.neighbour_within(i, j, 0) is True


# -- Pillar 2 / 3 / 4 tag coverage -----------------------------------


class TestPillarTagging:
    """Sweep: every pillar-establishing function or protocol in
    ``hierarchy.py`` must carry the ``@pillar(n)`` decorator's tag.
    Together with ``TestPillar1Tagging`` in ``test_substitution.py``,
    this verifies the four-pillar coverage is complete and stable
    against future refactors that might silently strip a decorator.
    """

    def test_is_recognisable_is_pillar_2(self) -> None:
        assert is_recognisable._pillar == 2

    def test_inflation_argument_is_pillar_3(self) -> None:
        assert inflation_argument._pillar == 3

    def test_fourth_pillar_argument_is_pillar_4(self) -> None:
        assert FourthPillarArgument._pillar == 4

    def test_all_four_pillars_have_at_least_one_tagged_target(self) -> None:
        # Cross-check against the pillar-1 functions in
        # apeiron.substitution. This test exists to enforce that no
        # pillar is unimplemented or silently un-tagged at the
        # module-graph level.
        from apeiron.substitution import is_primitive, perron_frobenius_in_zphi
        pillar_targets: dict[int, list[object]] = {
            1: [is_primitive, perron_frobenius_in_zphi],
            2: [is_recognisable],
            3: [inflation_argument],
            4: [FourthPillarArgument],
        }
        for pillar_n, targets in pillar_targets.items():
            for t in targets:
                assert getattr(t, "_pillar", None) == pillar_n, (
                    f"target {t} expected _pillar = {pillar_n}, "
                    f"got {getattr(t, '_pillar', None)}"
                )
