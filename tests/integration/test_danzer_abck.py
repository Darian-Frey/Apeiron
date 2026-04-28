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


# -- Q5b: stabiliser of each ABCK prototile in I_h --------------------


class TestDanzerProtileStabiliser:
    """Per Claude (web) Q5b ruling 2026-04-29: the stabiliser of each
    ABCK prototile in I_h is trivial.

    The argument from Frettlöh: the four vertices of each tetrahedron
    belong to distinct CPS vertex classes (I/II/III/IV), so no
    non-identity element of I_h can permute the vertices of a single
    tile while fixing the tile as a set. Therefore Stab_{I_h}(A) =
    Stab_{I_h}(B) = Stab_{I_h}(C) = Stab_{I_h}(K) = {identity}.

    This is the hypothesis under which ``position_signature`` (Q5a /
    Goodman-Strauss atlas form) elides the centre-stabiliser quotient
    for ABCK. If the assertion ever fails, ``position_signature`` for
    that prototile would need the explicit modulo-stabiliser pass.
    """

    @pytest.mark.parametrize("letter", ["A", "B", "C", "K"])
    def test_only_identity_fixes_tile(self, danzer_tiles, letter) -> None:
        from apeiron.symmetry import ICOSAHEDRAL, ImproperRotation, Rotation

        tile = danzer_tiles[letter]
        original = frozenset(tile.vertices)

        identity = Rotation.identity()
        fixed_by: list = []
        # 60 proper + 60 improper = 120 elements of I_h.
        for g in ICOSAHEDRAL:
            transformed = frozenset(g.apply(v) for v in tile.vertices)
            if transformed == original:
                fixed_by.append(g)
            transformed_imp = frozenset(
                ImproperRotation(g).apply(v) for v in tile.vertices
            )
            if transformed_imp == original:
                fixed_by.append(ImproperRotation(g))

        assert len(fixed_by) == 1, (
            f"Stab_{{I_h}}({letter}) has {len(fixed_by)} elements, "
            f"expected 1 (trivial). Non-identity stabilisers would "
            f"break the position_signature elision; see Q5b."
        )
        assert fixed_by[0] == identity, (
            f"Stab_{{I_h}}({letter}) contains a non-identity element."
        )


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


def _build_danzer_rule_from_paolini_dissection() -> SubstitutionRule:
    """Construct the canonical Danzer ``SubstitutionRule`` with real
    per-child geometry.

    **Source**: Paolini's POV-Ray substitution implementation, SVN
    trunk at
    ``https://svn.dmf.unicatt.it/svn/projects/animations/danzer/trunk/``,
    machine-extracted by ``candidates/danzer/_paolini_extract.py``
    into ``candidates/danzer/paolini_dissection.json``. Each of the
    25 children carries a translation in ×2-stored ℤ[φ]³ and a
    rotation in I (orientation-preserving, ``Rotation``) or in
    I_h \\ I (orientation-reversing, ``ImproperRotation``).

    **Why Paolini, not Koca** (arXiv 2003.13449, eqs. 19–31): per
    Claude (web) Q4a ruling 2026-04-29, Paolini is the canonical
    source. Both sources produce valid, isometrically-equivalent
    σ(K) dissections, but Koca's text leaves B's pose in σ(K)
    implicit (Fig. 8 only — no algebraic equation), and the natural
    "B at identity at origin" guess places B's vertex 1 outside
    τK. A canonical source that requires figure-reading for at
    least one child is not canonical. Paolini is complete (all 25
    children algebraically specified), and our extractor verified
    the result via three independent checks (volume conservation,
    multiplicity match, I_h membership; commits ``b3d8e55`` and
    ``b7174b1``).

    **Cross-validations against Koca (where Koca gives equations)**:

    - σ(K) K-child: Koca eq. (19) ``g_K = [[0,-1,0],[0,0,-1],[1,0,0]]``,
      ``t_K = ½(τ, τ², 0)``. Paolini's ``trK_2`` places the K-child at
      a different absolute position within τK but isometrically
      equivalent (pairwise distances among placed vertices match);
      not a per-child position match, but consistent on isometry-
      class. Frame difference, not bug.
    - σ(B), σ(C), σ(A): Koca's matrices appear in eqs. (20)–(31) but
      are not exhaustively cross-checked here; the Paolini volume
      conservation (sum of children's vols = τ³·vol(parent)) and
      multiplicity match against Frettlöh's matrix are the
      independent witnesses for these.

    **Children-only-from-Paolini (figure-dependent in Koca)**:

    - σ(K) B-child: only Paolini gives an explicit pose
      (``trK_1 = scale<-1,-1,-1>; translate(pt2)``, an
      ``ImproperRotation(identity)`` at translation
      ``(τ², τ², τ²)``). Koca's text describes the placement as
      "the origin coincides with one of the vertices of B" with
      the rest in Fig. 8 — not algebraically reproducible from text.
    """
    json_path = _DANZER_DIR / "paolini_dissection.json"
    if not json_path.exists():
        raise FileNotFoundError(
            f"{json_path} missing — run "
            "candidates/danzer/_paolini_extract.py to regenerate."
        )
    from apeiron.symmetry import ICOSAHEDRAL, ImproperRotation

    data = json.loads(json_path.read_text())
    type_to_idx = {"A": 0, "B": 1, "C": 2, "K": 3}

    by_parent: dict[str, list] = {"A": [], "B": [], "C": [], "K": []}
    for record in data["children"]:
        by_parent[record["parent"]].append(record)

    dissections: list[tuple[PositionedTile, ...]] = []
    for parent_letter in "ABCK":
        children: list[PositionedTile] = []
        for record in sorted(by_parent[parent_letter], key=lambda r: r["child_index"]):
            t_x2 = record["translation_x2"]
            translation = Vec3(
                ZPhi(*t_x2[0]), ZPhi(*t_x2[1]), ZPhi(*t_x2[2]),
            )
            proper = ICOSAHEDRAL[record["icosahedral_index"]]
            rotation = proper if record["is_proper"] else ImproperRotation(proper)
            children.append(PositionedTile(
                prototile_index=type_to_idx[record["child_type"]],
                translation=translation,
                rotation=rotation,
            ))
        dissections.append(tuple(children))

    return SubstitutionRule(
        n_prototiles=4,
        inflation=_phi_inflation_matrix(),
        dissections=tuple(dissections),
    )


@pytest.fixture(scope="module")
def danzer_rule() -> SubstitutionRule:
    """The canonical Danzer ``SubstitutionRule`` (Paolini real geometry).

    Per Claude (web) Q4a ruling 2026-04-29; see
    ``_build_danzer_rule_from_paolini_dissection`` for source citation
    and Koca cross-validations. Module-scoped.
    """
    return _build_danzer_rule_from_paolini_dissection()


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


class TestDanzerPatchBridge:
    """Exercise the algebraic→combinatorial pillar-2 bridge on Track
    A's first candidate: ``patch_from_supertile`` on the canonical
    Paolini real-geometry Danzer rule produces a ``TilePatch`` of the
    expected shape.

    With real geometry, leaves at level ≥ 1 land at distinct positions
    in ℤ[φ]³ (×2-stored), so the squared-Euclidean oracle separates
    them. ``is_recognisable`` with the multiset / shell signatures
    still returns False on this rule (the signatures are too coarse;
    see Q5 in ``candidates/danzer/dissection_notes.md``); the
    refined ``position_signature`` (Goodman-Strauss atlas form) is
    where pillar 2 actually decides.
    """

    def test_level_one_patch_matches_dissection_size(
        self, danzer_rule,
    ) -> None:
        from apeiron.hierarchy import patch_from_supertile
        # σ(A) per Frettlöh column 0: 0 A + 3 B + 2 C + 6 K = 11.
        patch = patch_from_supertile(
            danzer_rule, prototile_index=0, level=1,
            radius_step_squared=ZPhi(1, 0),
        )
        assert len(patch.tiles) == 11

    def test_level_one_parent_ids_are_unique(
        self, danzer_rule,
    ) -> None:
        # At level 1 each child IS its own level-1 supertile, so
        # parent_supertile values are 0..k-1 with no repeats.
        from apeiron.hierarchy import patch_from_supertile
        patch = patch_from_supertile(
            danzer_rule, prototile_index=0, level=1,
            radius_step_squared=ZPhi(1, 0),
        )
        parent_ids = [t.parent_supertile for t in patch.tiles]
        assert sorted(parent_ids) == list(range(11))

    def test_level_two_count_matches_matrix_power(
        self, danzer_rule,
    ) -> None:
        # σ²(A) leaf count = sum over level-1 children c of |σ(c)|.
        # Frettlöh dissection sizes: |σ(A)|=11, |σ(B)|=7, |σ(C)|=5,
        # |σ(K)|=2. σ(A) has 0 A + 3 B + 2 C + 6 K, so
        # σ²(A) has 0·11 + 3·7 + 2·5 + 6·2 = 0 + 21 + 10 + 12 = 43.
        from apeiron.hierarchy import patch_from_supertile
        patch = patch_from_supertile(
            danzer_rule, prototile_index=0, level=2,
            radius_step_squared=ZPhi(1, 0),
        )
        assert len(patch.tiles) == 43

    def test_level_two_parent_ids_partition_into_eleven_groups(
        self, danzer_rule,
    ) -> None:
        # The 43 level-2 leaves descend from the 11 level-1 children;
        # parent_supertile partitions them into 11 groups.
        from collections import Counter
        from apeiron.hierarchy import patch_from_supertile
        patch = patch_from_supertile(
            danzer_rule, prototile_index=0, level=2,
            radius_step_squared=ZPhi(1, 0),
        )
        groups = Counter(t.parent_supertile for t in patch.tiles)
        assert sorted(groups.keys()) == list(range(11))
        # Group sizes are |σ(child_type)| for each level-1 child.
        # σ(A) at column 0 has 3 B + 2 C + 6 K, so groups should
        # contain 3 of size 7 (B), 2 of size 5 (C), 6 of size 2 (K).
        size_counts = Counter(groups.values())
        assert size_counts == Counter({7: 3, 5: 2, 2: 6})

    def test_recognisable_fails_with_multiset_signature(
        self, danzer_rule,
    ) -> None:
        # The default multiset signature is provably too coarse for
        # ABCK (per Q5 in dissection_notes.md). is_recognisable must
        # report IndistinguishablePair until position_signature lands.
        from apeiron.hierarchy import (
            IndistinguishablePair,
            is_recognisable,
            patch_from_supertile,
        )
        patch = patch_from_supertile(
            danzer_rule, prototile_index=0, level=1,
            radius_step_squared=ZPhi(4, 0),
        )
        result = is_recognisable(patch, max_radius=3)
        assert result.is_recognisable is False
        assert isinstance(result.witness, IndistinguishablePair)


class TestPaoliniDissectionExtraction:
    """Smoke + sanity checks for the Paolini POV-Ray transcription
    (commit b3d8e55). The extractor itself is in
    ``candidates/danzer/_paolini_extract.py``; the JSON it produces is
    the canonical Paolini encoding under Q4a if Claude (web) selects
    that source.

    Volume conservation: for each σ(X), the sum of children's volumes
    equals τ³·vol(X). The child counts match Frettlöh's matrix.
    Every child rotation is recoverable as either a member of
    ``ICOSAHEDRAL`` (proper) or ``ImproperRotation(g)`` for some
    ``g ∈ ICOSAHEDRAL``.
    """

    def test_json_has_25_records_with_correct_counts(self) -> None:
        json_path = _DANZER_DIR / "paolini_dissection.json"
        if not json_path.exists():
            pytest.skip("paolini_dissection.json not yet generated")
        data = json.loads(json_path.read_text())
        records = data["children"]
        assert len(records) == 25
        # Per-parent counts.
        from collections import Counter
        parent_counts = Counter(r["parent"] for r in records)
        assert parent_counts == Counter({"A": 11, "B": 7, "C": 5, "K": 2})

    def test_per_parent_child_type_distribution_matches_frettloh(self) -> None:
        json_path = _DANZER_DIR / "paolini_dissection.json"
        if not json_path.exists():
            pytest.skip("paolini_dissection.json not yet generated")
        data = json.loads(json_path.read_text())
        from collections import Counter
        # Frettlöh column convention: σ(X) → counts.
        expected = {
            "A": Counter({"A": 0, "B": 3, "C": 2, "K": 6}),
            "B": Counter({"A": 0, "B": 2, "C": 1, "K": 4}),
            "C": Counter({"A": 1, "B": 0, "C": 2, "K": 2}),
            "K": Counter({"A": 0, "B": 1, "C": 0, "K": 1}),
        }
        actual = {
            parent: Counter(
                r["child_type"] for r in data["children"] if r["parent"] == parent
            )
            for parent in "ABCK"
        }
        for parent in "ABCK":
            for child_type in "ABCK":
                assert actual[parent].get(child_type, 0) == expected[parent].get(child_type, 0), (
                    f"σ({parent}) → {child_type}: expected "
                    f"{expected[parent].get(child_type, 0)}, "
                    f"got {actual[parent].get(child_type, 0)}"
                )

    def test_rotations_resolve_in_icosahedral(self) -> None:
        # Every record has icosahedral_index ∈ [0, 60).
        json_path = _DANZER_DIR / "paolini_dissection.json"
        if not json_path.exists():
            pytest.skip("paolini_dissection.json not yet generated")
        data = json.loads(json_path.read_text())
        for r in data["children"]:
            assert 0 <= r["icosahedral_index"] < 60, (
                f"{r['parent']} child {r['child_index']}: "
                f"icosahedral_index {r['icosahedral_index']} out of range"
            )
            assert isinstance(r["is_proper"], bool)

    def test_translations_are_zphi_pairs(self) -> None:
        # Every translation is three (a, b) integer pairs.
        json_path = _DANZER_DIR / "paolini_dissection.json"
        if not json_path.exists():
            pytest.skip("paolini_dissection.json not yet generated")
        data = json.loads(json_path.read_text())
        for r in data["children"]:
            t = r["translation_x2"]
            assert len(t) == 3
            for component in t:
                assert len(component) == 2
                a, b = component
                assert isinstance(a, int) and isinstance(b, int)

    def test_volume_conservation(self) -> None:
        # For each σ(X), sum of children's prototile volumes (in real
        # coords) equals τ³·vol(X). Pure floats — this is a oracle
        # check on Frettlöh's matrix multiplicities + child types,
        # not exact arithmetic.
        import math
        PHI = (1.0 + math.sqrt(5.0)) / 2.0

        def real_z(a: int, b: int) -> float:
            return a + b * PHI

        prototile_vertices = {
            "A": [(0,0,0), (real_z(1,2), 0, real_z(1,1)),
                  (real_z(1,1),)*3, (real_z(1,1), 1, 0)],
            "B": [(0,0,0), (real_z(1,2), 0, real_z(1,1)),
                  (real_z(1,1),)*3, (real_z(1,1), real_z(0,1), 1)],
            "C": [(0,0,0), (real_z(0,-1), 0, 1),
                  (real_z(1,1),)*3, (0, real_z(1,1), 1)],
            "K": [(0,0,0), (-1, real_z(0,1), 0),
                  (real_z(0,1),)*3, (-0.5, -0.5+0.5*PHI, 0.5*PHI)],
        }

        def vol(verts: list[tuple[float, float, float]]) -> float:
            v0, v1, v2, v3 = verts
            a = tuple(v1[i]-v0[i] for i in range(3))
            b = tuple(v2[i]-v0[i] for i in range(3))
            c = tuple(v3[i]-v0[i] for i in range(3))
            det = (
                a[0]*(b[1]*c[2] - b[2]*c[1])
              - a[1]*(b[0]*c[2] - b[2]*c[0])
              + a[2]*(b[0]*c[1] - b[1]*c[0])
            )
            return abs(det) / 6.0

        volumes = {k: vol(v) for k, v in prototile_vertices.items()}
        tau3 = PHI ** 3
        counts = {
            "A": {"A": 0, "B": 3, "C": 2, "K": 6},
            "B": {"A": 0, "B": 2, "C": 1, "K": 4},
            "C": {"A": 1, "B": 0, "C": 2, "K": 2},
            "K": {"A": 0, "B": 1, "C": 0, "K": 1},
        }
        for parent in "ABCK":
            expected = tau3 * volumes[parent]
            actual = sum(
                n * volumes[ct] for ct, n in counts[parent].items()
            )
            assert abs(actual - expected) < 1e-9, (
                f"σ({parent}): vol mismatch — expected {expected}, "
                f"got {actual}, ratio {actual/expected}"
            )


class TestDanzerFullPipeline:
    """27D — full end-to-end four-pillar pipeline on the canonical
    Paolini Danzer rule. Pillars 1 + 2 (position_signature) + 3 in
    one chain.

    Per Claude (web) Q5a 2026-04-29: position_signature is the
    Goodman-Strauss atlas form, the canonical pillar-2 signature
    for 3D substitution tilings. Earlier multiset / shell signatures
    are too coarse for ABCK and produce IndistinguishablePair
    counterexamples; position_signature uses exact ℤ[φ] relative-
    position-and-rotation tuples and should resolve them.
    """

    def test_level_one_leaves_are_geometrically_distinct(
        self, danzer_rule,
    ) -> None:
        # σ(A) = 11 children — at level 1 they land at distinct
        # ℤ[φ]³ positions (×2-stored). Some siblings may share a
        # translation if face-glued at the same vertex, so we don't
        # require 11 unique positions; we DO require > 1 (placeholder
        # rule had all leaves at origin, giving exactly 1).
        from apeiron.hierarchy import expand_supertile_with_parents
        tagged = expand_supertile_with_parents(danzer_rule, 0, 1)
        positions = {leaf.translation for leaf, _parent in tagged}
        assert len(positions) > 1

    @pytest.mark.parametrize(
        ("prototile_index", "letter", "expected_max_radius_used"),
        [
            (0, "A", 5),  # σ(A) at radius 2; 5 is a comfortable cap
            (1, "B", 8),  # σ(B) needs radius 6 — see radius-cap note
            (2, "C", 5),  # σ(C) at radius 1
            (3, "K", 5),  # σ(K) at radius 1
        ],
    )
    def test_pillar_2_succeeds_with_position_signature_level_1(
        self, danzer_rule, prototile_index, letter, expected_max_radius_used,
    ) -> None:
        # The 27D milestone, exercised across all four prototiles:
        # with position_signature on the canonical Paolini real-
        # geometry rule, pillar 2 succeeds for every σ(X), X ∈ {A,B,C,K}.
        #
        # Empirical recognisability radii (radius_step_squared=ZPhi(4,0),
        # i.e., r=k corresponds to k unit edge lengths of context):
        #   σ(A) at radius 2  (11 leaves)
        #   σ(B) at radius 6  ( 7 leaves) — needs the most context
        #   σ(C) at radius 1  ( 5 leaves)
        #   σ(K) at radius 1  ( 2 leaves)
        # The σ(B) requirement of 6 reflects that B's level-1 dissection
        # has multiple K-children that are I_h-related until distant
        # context disambiguates. Theoretical bound for an n=4 alphabet
        # with inflation φ is ~φ^8 ≈ 47; observed ≤ 6 is well within.
        from apeiron.hierarchy import (
            is_recognisable, patch_from_supertile, position_signature,
        )
        patch = patch_from_supertile(
            danzer_rule, prototile_index=prototile_index, level=1,
            radius_step_squared=ZPhi(4, 0),
        )
        recog = is_recognisable(
            patch, max_radius=expected_max_radius_used,
            signature_fn=position_signature,
        )
        assert recog.is_recognisable is True, (
            f"Pillar 2 failed at level=1 for σ({letter}) with "
            f"max_radius={expected_max_radius_used}; witness={recog.witness}"
        )

    def test_pillar_2_succeeds_at_level_2_for_A(
        self, danzer_rule,
    ) -> None:
        # σ²(A) has 43 leaves. Pillar 2 must succeed at level 2 too —
        # going up the supertile hierarchy must not break recognisability.
        # Empirically: radius_used = 4 with step²=ZPhi(4, 0).
        from apeiron.hierarchy import (
            is_recognisable, patch_from_supertile, position_signature,
        )
        patch = patch_from_supertile(
            danzer_rule, prototile_index=0, level=2,
            radius_step_squared=ZPhi(4, 0),
        )
        recog = is_recognisable(
            patch, max_radius=5, signature_fn=position_signature,
        )
        assert recog.is_recognisable is True, (
            f"Pillar 2 failed at level=2 for σ²(A); witness={recog.witness}"
        )

    def test_full_pipeline_emits_inflation_argument(
        self, danzer_rule,
    ) -> None:
        # Pillars 1 + 2 + 3 in one chain on the canonical rule. The
        # resulting InflationArgument is the Track A baseline witness
        # for sub-commit 27D: the rule is primitive, recognisable, and
        # therefore aperiodic by the inflation argument.
        from apeiron.hierarchy import (
            InflationArgument, inflation_argument, is_recognisable,
            patch_from_supertile, position_signature,
        )
        patch = patch_from_supertile(
            danzer_rule, prototile_index=0, level=1,
            radius_step_squared=ZPhi(4, 0),
        )
        recog = is_recognisable(
            patch, max_radius=5, signature_fn=position_signature,
        )
        result = inflation_argument(danzer_rule, recog)
        assert isinstance(result, InflationArgument)
        assert result.pf_eigenvalue == ZPhi(1, 2)  # φ³
        assert result.recognisability_radius == recog.radius_used

    def test_aperiodicity_witness_one_call_baseline(
        self, danzer_rule,
    ) -> None:
        # The deformation-search building block: aperiodicity_witness
        # runs pillars 1 + 2 + 3 across all four prototiles in one
        # call. The Track A baseline returns InflationArgument with
        # PF = φ³ and recognisability_radius = max across {σ(A)→2,
        # σ(B)→6, σ(C)→1, σ(K)→1} = 6.
        from apeiron.hierarchy import (
            InflationArgument, aperiodicity_witness,
        )
        result = aperiodicity_witness(
            danzer_rule, level=1, radius_step_squared=ZPhi(4, 0),
            max_radius=8,
        )
        assert isinstance(result, InflationArgument)
        assert result.pf_eigenvalue == ZPhi(1, 2)
        assert result.recognisability_radius == 6
