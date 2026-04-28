"""Danzer ABCK 4-tile baseline.

Track A (CLAUDE.md §6.1) starts here: encode Danzer's published
4-tile aperiodic set, run it through the pipeline, then perturb
toward a single-tile candidate via deformation search.

Sub-commit 27A landed prototile shape verification. **27B-α (this
file's additions)** builds the Danzer ``SubstitutionRule`` with
Frettlöh's matrix and placeholder geometric dissections — children
positioned at origin with identity rotation. The dissection's
*combinatorial* content (how many of each prototile per σ(X)) is
correct and matches Frettlöh exactly; the *geometric* content
(translation + rotation per child) is deferred to 27B-β, which
will transcribe Frettlöh Figure 2 / the Tilings Encyclopedia
interactive view into a ``dissection_notes.md`` sidecar before
encoding the full geometry.

This is sufficient for pillar-1 verification (which only consults
the substitution matrix derived from prototile_index counts), and
is what Claude (web) and the user agreed to as the matrix-only
acceptance step.

The fourth-pillar stub lands in 27C; the full pipeline run lands
in 27D.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from apeiron.polyhedron import Polyhedron
from apeiron.substitution import (
    PositionedTile,
    SubstitutionRule,
    is_primitive,
    perron_frobenius_in_zphi,
    substitution_matrix,
)
from apeiron.symmetry import Mat3, Rotation, Vec3
from apeiron.util import load_candidate
from apeiron.zphi import ZPhi

_DANZER_DIR = Path("candidates/danzer")


@pytest.fixture(scope="module")
def danzer_tiles() -> dict[str, Polyhedron]:
    """Load all four Danzer prototiles once per module."""
    return {
        letter: load_candidate(_DANZER_DIR / f"{letter}.json")
        for letter in ("A", "B", "C", "K")
    }


# -- structural shape acceptance --------------------------------------


class TestDanzerProtileShapes:
    @pytest.mark.parametrize("letter", ["A", "B", "C", "K"])
    def test_is_tetrahedron(self, danzer_tiles, letter) -> None:
        tile = danzer_tiles[letter]
        assert len(tile.vertices) == 4
        assert len(tile.faces) == 4

    @pytest.mark.parametrize("letter", ["A", "B", "C", "K"])
    def test_every_face_is_triangle(self, danzer_tiles, letter) -> None:
        tile = danzer_tiles[letter]
        for face in tile.faces:
            assert len(face) == 3

    @pytest.mark.parametrize("letter", ["A", "B", "C", "K"])
    def test_six_edges(self, danzer_tiles, letter) -> None:
        tile = danzer_tiles[letter]
        edges: set[tuple[int, int]] = set()
        for face in tile.faces:
            n = len(face)
            for i in range(n):
                u, v = face[i], face[(i + 1) % n]
                edges.add((min(u, v), max(u, v)))
        assert len(edges) == 6

    @pytest.mark.parametrize("letter", ["A", "B", "C", "K"])
    def test_euler_characteristic_two(self, danzer_tiles, letter) -> None:
        tile = danzer_tiles[letter]
        edges: set[tuple[int, int]] = set()
        for face in tile.faces:
            n = len(face)
            for i in range(n):
                u, v = face[i], face[(i + 1) % n]
                edges.add((min(u, v), max(u, v)))
        assert (
            len(tile.vertices) - len(edges) + len(tile.faces) == 2
        )


# -- prototiles are pairwise distinct ---------------------------------


class TestDanzerProtileDistinctness:
    def test_all_four_are_distinct_polyhedra(self, danzer_tiles) -> None:
        polyhedra = [danzer_tiles[letter] for letter in "ABCK"]
        assert len(set(polyhedra)) == 4


# -- K's half-integer storage -----------------------------------------


class TestDanzerKHalfInteger:
    """Tile K is the first fixture exercising scale_denom=2 in the
    wild. Frettlöh notes vertex 4 as class IV with half-integer
    coordinates relative to the icosahedral basis; the loader stores
    everything ×2 so the polyhedron is in pure ℤ[φ]³ post-load.
    """

    def test_k_loads_via_scale_denom_two(self, danzer_tiles) -> None:
        # No exception during fixture build is the assertion here;
        # this test just confirms the load completed.
        assert "K" in danzer_tiles

    def test_k_vertices_are_integer_zphi_after_loading(
        self, danzer_tiles,
    ) -> None:
        # After scale_denom=2 doubling on input then ÷2 in the loader,
        # every stored coordinate is an integer ZPhi (a, b ∈ ℤ).
        from apeiron.zphi import ZPhi
        K = danzer_tiles["K"]
        for v in K.vertices:
            for comp in (v.x, v.y, v.z):
                assert isinstance(comp, ZPhi)
                assert isinstance(comp.a, int)
                assert isinstance(comp.b, int)


# -- prototile_index metadata is internally consistent ----------------


class TestDanzerProtileIndices:
    """Each Danzer JSON declares ``prototile_index`` (0=A, 1=B, 2=C,
    3=K). Verify the on-disk indices match the convention used by
    the substitution matrix in ``danzer/substitution.json``.

    The loader currently silently ignores extra top-level keys, so
    this test reads ``prototile_index`` directly from the JSON.
    """

    def test_prototile_indices_match_letter_order(self) -> None:
        expected = {"A": 0, "B": 1, "C": 2, "K": 3}
        for letter, idx in expected.items():
            with (_DANZER_DIR / f"{letter}.json").open() as f:
                payload = json.load(f)
            assert payload["prototile_index"] == idx


# -- Sub-commit 27B-α: SubstitutionRule + pillar-1 verification -------


def _danzer_substitution_metadata() -> dict:
    """Read ``candidates/danzer/substitution.json`` once."""
    with (_DANZER_DIR / "substitution.json").open() as f:
        return json.load(f)


def _zero_vec3() -> Vec3:
    """Origin in ×2 storage (the dummy translation for placeholder
    children)."""
    return Vec3(ZPhi(0, 0), ZPhi(0, 0), ZPhi(0, 0))


def _phi_inflation_matrix() -> Mat3:
    """The Danzer linear inflation σ = φ·I.

    Per Claude (web)'s 2026-04-23 correction, the *linear* inflation
    factor is φ (= ZPhi(0, 1)); volume scales as the cube and that
    cube is the PF eigenvalue of the substitution matrix
    (φ³ = ZPhi(1, 2)). The two are not the same quantity; the
    inflation matrix carries the linear factor.
    """
    phi = ZPhi(0, 1)
    z = ZPhi(0, 0)
    return Mat3(
        Vec3(phi, z, z),
        Vec3(z, phi, z),
        Vec3(z, z, phi),
    )


def _placeholder_dissection(
    matrix_column: list[int],
) -> tuple[PositionedTile, ...]:
    """Build a tuple of ``PositionedTile`` whose ``prototile_index``
    multiset matches ``matrix_column``.

    All children share the origin / identity-rotation placement.
    Geometric content is deliberately stubbed; this commit's job is
    to verify the substitution-matrix shape (pillar 1 only).
    Sub-commit 27B-β replaces these with real translations and
    rotations transcribed from Frettlöh Figure 2.
    """
    children: list[PositionedTile] = []
    for type_idx, count in enumerate(matrix_column):
        for _ in range(count):
            children.append(
                PositionedTile(
                    prototile_index=type_idx,
                    translation=_zero_vec3(),
                    rotation=Rotation.identity(),
                )
            )
    return tuple(children)


def _build_danzer_rule_with_placeholder_geometry() -> SubstitutionRule:
    """Construct a Danzer ``SubstitutionRule`` with Frettlöh's matrix
    and placeholder-only dissection geometry (origin + identity for
    every child).

    The substitution matrix is column-conventioned: column ``i`` of
    ``M`` is the count vector of σ(prototile_i)'s children by type.
    Frettlöh's matrix:

    ::

        M  =  A  B  C  K
           A [[0, 0, 1, 0],
           B  [3, 2, 0, 1],
           C  [2, 1, 2, 0],
           K  [6, 4, 2, 1]]
    """
    metadata = _danzer_substitution_metadata()
    matrix = metadata["substitution_matrix"]
    n = len(matrix)
    dissections = tuple(
        _placeholder_dissection([matrix[type_idx][col] for type_idx in range(n)])
        for col in range(n)
    )
    return SubstitutionRule(
        n_prototiles=n,
        inflation=_phi_inflation_matrix(),
        dissections=dissections,
    )


@pytest.fixture(scope="module")
def danzer_rule() -> SubstitutionRule:
    """The Danzer ``SubstitutionRule`` (matrix-only; placeholder
    dissection geometry).

    Module-scoped so every pillar-1 test reuses a single
    construction.
    """
    return _build_danzer_rule_with_placeholder_geometry()


class TestDanzerSubstitutionMatrix:
    """Pillar 1 — the substitution-matrix recovery and primitivity
    checks against Frettlöh's published matrix.
    """

    def test_matrix_shape(self, danzer_rule) -> None:
        m = substitution_matrix(danzer_rule)
        assert m.shape == (4, 4)

    def test_matrix_matches_frettloh(self, danzer_rule) -> None:
        # Reconstructed matrix must equal the on-disk reference
        # (transcribed from Frettlöh's report).
        expected = _danzer_substitution_metadata()["substitution_matrix"]
        m = substitution_matrix(danzer_rule)
        for j in range(4):
            for i in range(4):
                assert m[j, i] == expected[j][i], (
                    f"M[{j},{i}] = {m[j, i]} but Frettlöh says "
                    f"{expected[j][i]}"
                )

    def test_matrix_column_sums_match_child_counts(
        self, danzer_rule,
    ) -> None:
        # Column j of M is the multiset count vector of σ(prototile_j)'s
        # children — its sum is the total number of children. Useful
        # cross-check that the placeholder dissection construction
        # didn't drop any children.
        m = substitution_matrix(danzer_rule)
        # σ(A) = 11 children (3+2+6), σ(B) = 7, σ(C) = 5, σ(K) = 2.
        assert int(m[:, 0].sum()) == 11
        assert int(m[:, 1].sum()) == 7
        assert int(m[:, 2].sum()) == 5
        assert int(m[:, 3].sum()) == 2

    def test_matrix_is_primitive(self, danzer_rule) -> None:
        # Pillar 1 acceptance: M is primitive (some power has all-
        # positive entries). Wielandt bound for a 4×4 matrix is
        # (4-1)² + 1 = 10.
        m = substitution_matrix(danzer_rule)
        assert is_primitive(m)


class TestDanzerPerronFrobenius:
    """Pillar 1 — Perron–Frobenius eigenvalue extraction. Frettlöh
    states the PF eigenvalue is φ³ = ZPhi(1, 2); the other
    eigenvalues are τ, −τ⁻¹, −τ⁻³ (all in ℤ[φ]).
    """

    def test_pf_eigenvalue_is_phi_cubed(self, danzer_rule) -> None:
        m = substitution_matrix(danzer_rule)
        pf = perron_frobenius_in_zphi(m)
        assert pf is not None, (
            "PF eigenvalue not recovered in ℤ[φ]; expected φ³ = "
            "ZPhi(1, 2)"
        )
        assert pf == ZPhi(1, 2), f"expected ZPhi(1, 2), got {pf}"

    def test_pf_eigenvalue_matches_metadata(
        self, danzer_rule,
    ) -> None:
        # The on-disk substitution metadata declares the expected PF.
        m = substitution_matrix(danzer_rule)
        pf = perron_frobenius_in_zphi(m)
        meta_pf = _danzer_substitution_metadata()["expected_pf_eigenvalue"]
        assert pf is not None
        assert pf == ZPhi(meta_pf[0], meta_pf[1])


class TestDanzerInflationMatrix:
    """Pillar-1 metadata cross-check: the inflation matrix on the
    SubstitutionRule encodes the *linear* factor (φ), not the
    volume-scaling factor (φ³ = PF eigenvalue).
    """

    def test_inflation_is_phi_times_identity(
        self, danzer_rule,
    ) -> None:
        # Row 0 column 0: ZPhi(0, 1) = φ. Off-diagonal: 0.
        m = danzer_rule.inflation
        phi = ZPhi(0, 1)
        zero = ZPhi(0, 0)
        assert m.row0.x == phi
        assert m.row0.y == zero
        assert m.row0.z == zero
        assert m.row1.x == zero
        assert m.row1.y == phi
        assert m.row1.z == zero
        assert m.row2.x == zero
        assert m.row2.y == zero
        assert m.row2.z == phi

    def test_inflation_factor_in_metadata_is_phi(self) -> None:
        # The on-disk metadata records [0, 1] = ZPhi(0, 1) = φ as
        # the linear inflation factor. Guards against a future
        # confusion-with-PF regression.
        meta = _danzer_substitution_metadata()
        assert meta["inflation_linear_factor"] == [0, 1]


# -- Sub-commit 27C: fourth-pillar stub -------------------------------


class TestDanzerFourthPillarStub:
    """The ``candidates/danzer/fourth_pillar.py`` stub satisfies the
    ``FourthPillarArgument`` protocol, raises ``NotImplementedError``
    on its methods, and cites the published proof.

    First end-to-end exercise of the
    ``apeiron.hierarchy.FourthPillarArgument`` protocol's import +
    instantiation path. Future deformation-search candidates will
    inherit the ``candidates/<name>/fourth_pillar.py`` convention
    that this module establishes.
    """

    def test_module_imports(self) -> None:
        from candidates.danzer import fourth_pillar
        assert hasattr(fourth_pillar, "DanzerABCKFourthPillar")

    def test_class_instantiates(self) -> None:
        from candidates.danzer.fourth_pillar import DanzerABCKFourthPillar
        impl = DanzerABCKFourthPillar()
        assert impl is not None

    def test_satisfies_protocol(self) -> None:
        from apeiron.hierarchy import FourthPillarArgument
        from candidates.danzer.fourth_pillar import DanzerABCKFourthPillar
        impl = DanzerABCKFourthPillar()
        assert isinstance(impl, FourthPillarArgument)

    def test_local_configurations_raises_not_implemented(self) -> None:
        from candidates.danzer.fourth_pillar import DanzerABCKFourthPillar
        impl = DanzerABCKFourthPillar()
        with pytest.raises(NotImplementedError):
            impl.local_configurations()

    def test_verify_hierarchical_raises_not_implemented(self) -> None:
        from candidates.danzer.fourth_pillar import DanzerABCKFourthPillar
        from apeiron.corona import CoronaConfig
        # Build any valid CoronaConfig; the stub should refuse
        # regardless of input.
        from apeiron.polyhedron import Polyhedron
        any_tetra_verts = [
            Vec3(ZPhi(0, 0), ZPhi(0, 0), ZPhi(0, 0)),
            Vec3(ZPhi(2, 0), ZPhi(0, 0), ZPhi(0, 0)),
            Vec3(ZPhi(0, 0), ZPhi(2, 0), ZPhi(0, 0)),
            Vec3(ZPhi(0, 0), ZPhi(0, 0), ZPhi(2, 0)),
        ]
        any_polytope = Polyhedron.from_vertices(any_tetra_verts)
        any_config = CoronaConfig.from_neighbours(any_polytope, [])
        impl = DanzerABCKFourthPillar()
        with pytest.raises(NotImplementedError):
            impl.verify_hierarchical(any_config)

    def test_error_message_contains_citation(self) -> None:
        from candidates.danzer.fourth_pillar import DanzerABCKFourthPillar
        impl = DanzerABCKFourthPillar()
        try:
            impl.local_configurations()
        except NotImplementedError as exc:
            msg = str(exc)
            # All three published sources should be referenced.
            assert "Danzer 1989" in msg
            assert "Goodman-Strauss 1998" in msg
            assert "Frettlöh" in msg or "Frettloh" in msg
        else:
            raise AssertionError("expected NotImplementedError")


# -- Pillars 1 + 3 end-to-end through hierarchy.inflation_argument ----


class TestDanzerPillarsOneAndThree:
    """End-to-end exercise of ``hierarchy.inflation_argument`` on the
    Danzer rule using a synthetic ``RecognisabilityResult`` for the
    pillar-2 input.

    This complements the unit-level pillar-3 tests in
    ``test_hierarchy.py`` (which use a hand-built P3 rule) by
    confirming the same predicate also produces a valid
    ``InflationArgument`` on Track A's first real candidate. Until
    27B-β lands the geometric dissection and real recognisability
    can be computed on a Danzer patch, the pillar-2 input here is
    a synthetic ``is_recognisable=True`` witness — same pattern as
    ``test_hierarchy.py``'s P3 tests.
    """

    def _success_recognisability(self, radius: int = 1):
        from apeiron.hierarchy import RecognisabilityResult
        return RecognisabilityResult(
            is_recognisable=True,
            radius_used=radius,
            radius_cap_reached=False,
            witness={(0,): 0, (1,): 1, (2,): 2, (3,): 3},
        )

    def test_inflation_argument_succeeds_on_danzer(
        self, danzer_rule,
    ) -> None:
        from apeiron.hierarchy import (
            InflationArgument,
            inflation_argument,
        )
        result = inflation_argument(
            danzer_rule, self._success_recognisability(),
        )
        assert isinstance(result, InflationArgument)

    def test_inflation_argument_pf_is_phi_cubed(
        self, danzer_rule,
    ) -> None:
        from apeiron.hierarchy import (
            InflationArgument,
            inflation_argument,
        )
        result = inflation_argument(
            danzer_rule, self._success_recognisability(),
        )
        assert isinstance(result, InflationArgument)
        assert result.pf_eigenvalue == ZPhi(1, 2)

    def test_inflation_argument_records_recognisability_radius(
        self, danzer_rule,
    ) -> None:
        from apeiron.hierarchy import (
            InflationArgument,
            inflation_argument,
        )
        result = inflation_argument(
            danzer_rule, self._success_recognisability(radius=3),
        )
        assert isinstance(result, InflationArgument)
        assert result.recognisability_radius == 3

    def test_inflation_argument_fails_on_failing_recognisability(
        self, danzer_rule,
    ) -> None:
        # If pillar 2 fails (synthetic is_recognisable=False),
        # pillar 3 must report "not recognisable" rather than
        # producing a witness.
        from apeiron.hierarchy import (
            IndistinguishablePair,
            InflationFailure,
            RecognisabilityResult,
            inflation_argument,
        )
        bad_result = RecognisabilityResult(
            is_recognisable=False,
            radius_used=5,
            radius_cap_reached=True,
            witness=IndistinguishablePair(
                tile_a=0, tile_b=1, radius=5, parent_a=0, parent_b=1,
            ),
        )
        result = inflation_argument(danzer_rule, bad_result)
        assert isinstance(result, InflationFailure)
        assert result.reason == "not recognisable"
