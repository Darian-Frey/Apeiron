"""Tests for apeiron.deformation — Track A Phase 1.5 face-merge
deformation-search infrastructure."""

from __future__ import annotations

import json
import shutil
import tempfile
from pathlib import Path

import pytest

from apeiron.deformation import (
    FaceAdjacentPair,
    MergeFiltering,
    enumerate_face_merge_compositions,
    face_adjacent_pairs,
    feasibility_upper_bound,
    merge_two_tiles,
    placed_vertices,
    prioritise_merge_candidate,
    scaffold_merge_candidate,
)
from apeiron.substitution import PositionedTile, SubstitutionRule
from apeiron.symmetry import ICOSAHEDRAL, ImproperRotation, Mat3, Rotation, Vec3
from apeiron.util import load_candidate
from apeiron.zphi import ZPhi


@pytest.fixture
def tmp_candidates_dir() -> Path:
    """Temp directory for scaffold tests; cleaned up after each test."""
    tmp = Path(tempfile.mkdtemp())
    yield tmp
    shutil.rmtree(tmp)


# -- scaffold_merge_candidate (Q6c) -----------------------------------


class TestScaffoldMergeCandidate:
    """Per Claude (web) Q6c ruling 2026-04-29: every face-merge
    candidate gets an auto-generated ``fourth_pillar.py`` stub at
    scaffold time, regardless of whether it later survives pillars
    1+2+3.
    """

    def test_creates_directory(self, tmp_candidates_dir: Path) -> None:
        path = scaffold_merge_candidate(
            name="merge_ab_test",
            parent_letters=("A", "B"),
            shared_face_description="face 2 of A coincides with face 0 of B",
            candidates_dir=tmp_candidates_dir,
        )
        assert path.exists()
        assert path.is_dir()
        assert path == tmp_candidates_dir / "merge_ab_test"

    def test_scaffolds_fourth_pillar_and_init(
        self, tmp_candidates_dir: Path,
    ) -> None:
        path = scaffold_merge_candidate(
            name="merge_ab_test",
            parent_letters=("A", "B"),
            shared_face_description="shared face details",
            candidates_dir=tmp_candidates_dir,
        )
        assert (path / "fourth_pillar.py").exists()
        assert (path / "__init__.py").exists()

    def test_fourth_pillar_module_imports(
        self, tmp_candidates_dir: Path,
    ) -> None:
        path = scaffold_merge_candidate(
            name="merge_ab_test",
            parent_letters=("A", "B"),
            shared_face_description="shared face",
            candidates_dir=tmp_candidates_dir,
        )
        # Dynamically import the generated module.
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            "test_scaffold_module", path / "fourth_pillar.py",
        )
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        # The PascalCase'd class name must exist.
        cls = mod.MergeAbTestFourthPillar
        assert cls is not None

    def test_generated_class_satisfies_protocol(
        self, tmp_candidates_dir: Path,
    ) -> None:
        from apeiron.hierarchy import FourthPillarArgument

        path = scaffold_merge_candidate(
            name="merge_ab_test",
            parent_letters=("A", "B"),
            shared_face_description="shared face",
            candidates_dir=tmp_candidates_dir,
        )
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            "test_scaffold_module", path / "fourth_pillar.py",
        )
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        instance = mod.MergeAbTestFourthPillar()
        assert isinstance(instance, FourthPillarArgument)

    def test_methods_raise_not_implemented_with_context(
        self, tmp_candidates_dir: Path,
    ) -> None:
        from apeiron.corona import CoronaConfig
        path = scaffold_merge_candidate(
            name="merge_test_ctx",
            parent_letters=("A", "B"),
            shared_face_description="distinctive description string",
            candidates_dir=tmp_candidates_dir,
        )
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            "test_scaffold_ctx", path / "fourth_pillar.py",
        )
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        instance = mod.MergeTestCtxFourthPillar()
        # The error must include the context (parent letters + face desc).
        with pytest.raises(NotImplementedError, match="distinctive description"):
            instance.local_configurations()
        with pytest.raises(NotImplementedError, match="distinctive description"):
            instance.verify_hierarchical(_config=None)  # type: ignore[arg-type]

    def test_rejects_existing_directory(
        self, tmp_candidates_dir: Path,
    ) -> None:
        # First call creates; second call refuses to overwrite.
        scaffold_merge_candidate(
            name="merge_ab_test",
            parent_letters=("A", "B"),
            shared_face_description="x",
            candidates_dir=tmp_candidates_dir,
        )
        with pytest.raises(FileExistsError):
            scaffold_merge_candidate(
                name="merge_ab_test",
                parent_letters=("A", "B"),
                shared_face_description="x",
                candidates_dir=tmp_candidates_dir,
            )

    def test_rejects_invalid_name(
        self, tmp_candidates_dir: Path,
    ) -> None:
        for bad_name in ("", "../escape", ".hidden", "with/slash"):
            with pytest.raises(ValueError, match="directory-safe"):
                scaffold_merge_candidate(
                    name=bad_name,
                    parent_letters=("A",),
                    shared_face_description="x",
                    candidates_dir=tmp_candidates_dir,
                )

    def test_pillar_4_approach_appears_in_docstring(
        self, tmp_candidates_dir: Path,
    ) -> None:
        # Custom approach text is captured in the generated module.
        custom_approach = (
            "Inherit from Penrose P3's matching-rule restriction and "
            "verify the merged tile's edge labels are I_h-equivariant."
        )
        path = scaffold_merge_candidate(
            name="merge_with_custom_approach",
            parent_letters=("A", "K"),
            shared_face_description="shared face",
            pillar_4_approach=custom_approach,
            candidates_dir=tmp_candidates_dir,
        )
        text = (path / "fourth_pillar.py").read_text()
        assert custom_approach in text

    def test_parent_letters_recorded(
        self, tmp_candidates_dir: Path,
    ) -> None:
        path = scaffold_merge_candidate(
            name="merge_abk_test",
            parent_letters=("A", "B", "K"),
            shared_face_description="x",
            candidates_dir=tmp_candidates_dir,
        )
        text = (path / "fourth_pillar.py").read_text()
        assert "A+B+K" in text


# -- Q6a: face_adjacent_pairs + merge_two_tiles -----------------------


_DANZER = Path("candidates/danzer")


def _phi_inflation() -> Mat3:
    phi = ZPhi(0, 1); z = ZPhi(0, 0)
    return Mat3(Vec3(phi, z, z), Vec3(z, phi, z), Vec3(z, z, phi))


@pytest.fixture(scope="module")
def paolini_rule() -> SubstitutionRule:
    """Build the canonical Paolini-real-geometry Danzer rule from
    ``candidates/danzer/paolini_dissection.json``."""
    data = json.loads((_DANZER / "paolini_dissection.json").read_text())
    type_to_idx = {"A": 0, "B": 1, "C": 2, "K": 3}
    by_parent: dict[str, list] = {"A": [], "B": [], "C": [], "K": []}
    for r in data["children"]:
        by_parent[r["parent"]].append(r)
    dissections = []
    for parent in "ABCK":
        children = []
        for r in sorted(by_parent[parent], key=lambda r: r["child_index"]):
            t = r["translation_x2"]
            translation = Vec3(ZPhi(*t[0]), ZPhi(*t[1]), ZPhi(*t[2]))
            proper = ICOSAHEDRAL[r["icosahedral_index"]]
            rotation = proper if r["is_proper"] else ImproperRotation(proper)
            children.append(PositionedTile(
                prototile_index=type_to_idx[r["child_type"]],
                translation=translation, rotation=rotation,
            ))
        dissections.append(tuple(children))
    return SubstitutionRule(
        n_prototiles=4, inflation=_phi_inflation(),
        dissections=tuple(dissections),
    )


@pytest.fixture(scope="module")
def abck_prototiles():
    """Load A, B, C, K Polyhedra in prototile-index order (0=A, 1=B, 2=C, 3=K)."""
    return [load_candidate(_DANZER / f"{x}.json") for x in "ABCK"]


class TestPlacedVertices:
    """Smoke tests for the placed_vertices helper."""

    def test_identity_pose_returns_local_vertices(self, abck_prototiles) -> None:
        # Identity rotation + zero translation = unchanged vertices.
        proto = abck_prototiles[0]
        zero = Vec3(ZPhi(0, 0), ZPhi(0, 0), ZPhi(0, 0))
        result = placed_vertices(proto, (zero, Rotation.identity()))
        assert result == proto.vertices

    def test_translation_only_shifts_all(self, abck_prototiles) -> None:
        proto = abck_prototiles[0]
        offset = Vec3(ZPhi(2, 0), ZPhi(4, 0), ZPhi(0, 2))
        result = placed_vertices(proto, (offset, Rotation.identity()))
        for v_local, v_placed in zip(proto.vertices, result, strict=True):
            assert v_placed == v_local + offset


class TestFaceAdjacentPair:
    """Validation of the FaceAdjacentPair dataclass."""

    def test_rejects_i_ge_j(self) -> None:
        with pytest.raises(ValueError, match="i < j"):
            FaceAdjacentPair(
                parent_index=0, i=2, j=1,
                shared_vertices=frozenset([
                    Vec3(ZPhi(a, 0), ZPhi(0, 0), ZPhi(0, 0))
                    for a in (0, 1, 2)
                ]),
            )

    def test_rejects_wrong_shared_count(self) -> None:
        # Tetrahedral face = 3 vertices; anything else rejected.
        bad = frozenset([Vec3(ZPhi(0, 0), ZPhi(0, 0), ZPhi(0, 0))])
        with pytest.raises(ValueError, match="exactly 3"):
            FaceAdjacentPair(
                parent_index=0, i=0, j=1, shared_vertices=bad,
            )


class TestFaceAdjacentPairsOnDanzer:
    """Per Claude (web) Q6a 2026-04-29: enumerate face-adjacent pairs
    for each parent's σ in the canonical Paolini Danzer rule. The
    counts are an empirical baseline — total ≤ O(50) per Q6b.
    """

    def test_sigma_K_has_exactly_one_pair(
        self, paolini_rule, abck_prototiles,
    ) -> None:
        # σ(K) = B + K, glued on a triangular face. Exactly one
        # face-adjacent pair (the dissection itself).
        pairs = face_adjacent_pairs(paolini_rule, 3, abck_prototiles)
        assert len(pairs) == 1
        pair = pairs[0]
        assert pair.parent_index == 3
        assert (pair.i, pair.j) == (0, 1)
        assert len(pair.shared_vertices) == 3

    def test_sigma_C_pair_count(
        self, paolini_rule, abck_prototiles,
    ) -> None:
        pairs = face_adjacent_pairs(paolini_rule, 2, abck_prototiles)
        assert len(pairs) == 4

    def test_sigma_B_pair_count(
        self, paolini_rule, abck_prototiles,
    ) -> None:
        pairs = face_adjacent_pairs(paolini_rule, 1, abck_prototiles)
        assert len(pairs) == 8

    def test_sigma_A_pair_count(
        self, paolini_rule, abck_prototiles,
    ) -> None:
        pairs = face_adjacent_pairs(paolini_rule, 0, abck_prototiles)
        assert len(pairs) == 11

    def test_total_pairs_within_q6b_estimate(
        self, paolini_rule, abck_prototiles,
    ) -> None:
        # Claude (web) Q6b 2026-04-29: "Call it O(50) merge candidates
        # before symmetry reduction." Empirically 11+8+4+1 = 24, well
        # within bounds and also small enough for brute-force pillars
        # 1+2+3 verification per candidate.
        total = sum(
            len(face_adjacent_pairs(paolini_rule, p, abck_prototiles))
            for p in range(4)
        )
        assert total == 24

    def test_pair_indices_in_range(
        self, paolini_rule, abck_prototiles,
    ) -> None:
        for parent_idx in range(4):
            n_children = len(paolini_rule.dissections[parent_idx])
            for pair in face_adjacent_pairs(
                paolini_rule, parent_idx, abck_prototiles,
            ):
                assert 0 <= pair.i < pair.j < n_children

    def test_rejects_out_of_range_parent(
        self, paolini_rule, abck_prototiles,
    ) -> None:
        with pytest.raises(ValueError, match="parent_index"):
            face_adjacent_pairs(paolini_rule, 99, abck_prototiles)


class TestMergeTwoTiles:
    """Per Claude (web) Q6a 2026-04-29: tetrahedra glued on a shared
    triangular face produce a 5-vertex / 9-edge / 6-face polyhedron
    (Euler χ = 2). Tests confirm the combinatorics on each of the
    24 face-adjacent pairs in Danzer.
    """

    def _check_merge_combinatorics(self, merged) -> None:
        assert len(merged.vertices) == 5
        assert len(merged.faces) == 6
        edges: set[tuple[int, int]] = set()
        for face in merged.faces:
            n = len(face)
            for k in range(n):
                edges.add(tuple(sorted([face[k], face[(k + 1) % n]])))
        assert len(edges) == 9
        assert (
            len(merged.vertices) - len(edges) + len(merged.faces) == 2
        )
        for face in merged.faces:
            assert len(face) == 3

    def test_sigma_K_merge(self, paolini_rule, abck_prototiles) -> None:
        # σ(K)'s only pair: B + K.
        pairs = face_adjacent_pairs(paolini_rule, 3, abck_prototiles)
        pair = pairs[0]
        ci = paolini_rule.dissections[3][pair.i]
        cj = paolini_rule.dissections[3][pair.j]
        merged = merge_two_tiles(
            abck_prototiles[ci.prototile_index],
            (ci.translation, ci.rotation),
            abck_prototiles[cj.prototile_index],
            (cj.translation, cj.rotation),
            pair.shared_vertices,
        )
        self._check_merge_combinatorics(merged)

    def test_every_face_adjacent_pair_merges_cleanly(
        self, paolini_rule, abck_prototiles,
    ) -> None:
        # All 24 pairs across σ(A,B,C,K) must produce a well-formed
        # 5-vertex / 6-face merged polyhedron — sanity check on the
        # whole face-adjacency enumeration.
        for parent_idx in range(4):
            for pair in face_adjacent_pairs(
                paolini_rule, parent_idx, abck_prototiles,
            ):
                ci = paolini_rule.dissections[parent_idx][pair.i]
                cj = paolini_rule.dissections[parent_idx][pair.j]
                merged = merge_two_tiles(
                    abck_prototiles[ci.prototile_index],
                    (ci.translation, ci.rotation),
                    abck_prototiles[cj.prototile_index],
                    (cj.translation, cj.rotation),
                    pair.shared_vertices,
                )
                self._check_merge_combinatorics(merged)

    def test_rejects_wrong_shared_count(self, abck_prototiles) -> None:
        zero = Vec3(ZPhi(0, 0), ZPhi(0, 0), ZPhi(0, 0))
        with pytest.raises(ValueError, match="3 shared"):
            merge_two_tiles(
                abck_prototiles[0],
                (zero, Rotation.identity()),
                abck_prototiles[1],
                (zero, Rotation.identity()),
                frozenset([zero]),  # only 1 vertex
            )


# -- Q6b: prioritise_merge_candidate ----------------------------------


class TestMergeFiltering:
    """Validation of the MergeFiltering dataclass + passes_all logic."""

    def test_passes_all_requires_all_filters(self) -> None:
        # Filter 1 fails → passes_all False.
        f = MergeFiltering(
            combinatorial_feasibility=False,
            combinatorial_k=None,
            face_count_strictly_decreases=True,
            vertex_class_compatible=True,
            detail="x",
        )
        assert f.passes_all is False

    def test_passes_all_when_all_true(self) -> None:
        f = MergeFiltering(
            combinatorial_feasibility=True,
            combinatorial_k=4,
            face_count_strictly_decreases=True,
            vertex_class_compatible=True,
            detail="x",
        )
        assert f.passes_all is True

    def test_unimplemented_filter_3_treated_as_pass(self) -> None:
        # vertex_class_compatible=None means "not checked" — gate
        # accepts. Once the filter lands, None → True/False.
        f = MergeFiltering(
            combinatorial_feasibility=True,
            combinatorial_k=2,
            face_count_strictly_decreases=True,
            vertex_class_compatible=None,
            detail="x",
        )
        assert f.passes_all is True


class TestPrioritiseMergeCandidateOnDanzer:
    """Per Claude (web) Q6b 2026-04-29: cheap filters applied before
    pillar-2 BFS. For ABCK 2-tile face-merges, filter 1 (combinatorial
    feasibility) is provably incompatible — no σ-count vector matches
    any 2-tile composition with the right structure. The filter
    implementation here confirms this and constitutes no-go evidence
    per the roadmap's Phase 1.5.
    """

    def test_no_face_merge_passes_combinatorial_feasibility(
        self, paolini_rule, abck_prototiles,
    ) -> None:
        # Across all 24 face-adjacent pairs in Danzer, no merge has
        # σ-count vector = k·composition for any integer k. This is
        # the key no-go finding.
        survivors = []
        for parent_idx in range(4):
            for pair in face_adjacent_pairs(
                paolini_rule, parent_idx, abck_prototiles,
            ):
                result = prioritise_merge_candidate(
                    paolini_rule, pair, abck_prototiles,
                )
                if result.combinatorial_feasibility:
                    survivors.append((parent_idx, pair))
        assert survivors == [], (
            f"Expected no Danzer face-merge to be combinatorially "
            f"feasible (no-go evidence per Q6b 2026-04-29); got "
            f"{len(survivors)} survivors: {survivors}"
        )

    def test_face_count_filter_passes_for_all_tetrahedral_merges(
        self, paolini_rule, abck_prototiles,
    ) -> None:
        # Tetrahedral merges always satisfy F_merged = F_a + F_b - 2
        # < F_a + F_b. Filter 2 is essentially a sanity check for
        # tetrahedral candidates; it's the prepared check for future
        # polyhedral candidates.
        for parent_idx in range(4):
            for pair in face_adjacent_pairs(
                paolini_rule, parent_idx, abck_prototiles,
            ):
                result = prioritise_merge_candidate(
                    paolini_rule, pair, abck_prototiles,
                )
                assert result.face_count_strictly_decreases is True

    def test_filter_1_detail_records_count_vectors(
        self, paolini_rule, abck_prototiles,
    ) -> None:
        # The MergeFiltering.detail string captures composition + σ_count
        # for each candidate — useful for debugging the no-go in a
        # future investigation.
        pair = face_adjacent_pairs(paolini_rule, 3, abck_prototiles)[0]
        result = prioritise_merge_candidate(
            paolini_rule, pair, abck_prototiles,
        )
        assert "composition=" in result.detail
        assert "σ_count=" in result.detail
        assert result.combinatorial_feasibility is False

    def test_passes_all_aggregate_is_zero_for_danzer(
        self, paolini_rule, abck_prototiles,
    ) -> None:
        # The single bottom-line aggregate: 0 of 24 Danzer face-merge
        # candidates pass all filters. The full deformation search
        # over face-merges of ABCK terminates here without proceeding
        # to pillar-1+2+3 verification.
        passes = 0
        total = 0
        for parent_idx in range(4):
            for pair in face_adjacent_pairs(
                paolini_rule, parent_idx, abck_prototiles,
            ):
                total += 1
                result = prioritise_merge_candidate(
                    paolini_rule, pair, abck_prototiles,
                )
                if result.passes_all:
                    passes += 1
        assert total == 24
        assert passes == 0


# -- Phase 1.6: feasibility_upper_bound + k-tuple enumeration ---------


def _trivial_pf2_rule() -> SubstitutionRule:
    """1-prototile rule σ(0) = 2 copies of 0 → M = [[2]], PF = 2.
    Used as the smoke-test counterpart for Phase 1.6 — its
    feasibility_upper_bound returns True and every composition
    survives the k-tuple enumeration.
    """
    z = ZPhi(0, 0)
    origin = Vec3(z, z, z)
    g = ICOSAHEDRAL[0]
    child = PositionedTile(prototile_index=0, translation=origin, rotation=g)
    return SubstitutionRule(
        n_prototiles=1, inflation=Mat3.identity(),
        dissections=((child, child),),
    )


def _phi_squared_rule() -> SubstitutionRule:
    """2-prototile rule with substitution matrix [[2, 1], [1, 1]],
    eigenvalues φ² and 1/φ². Neither is a positive integer, so
    feasibility_upper_bound returns False.
    """
    z = ZPhi(0, 0)
    origin = Vec3(z, z, z)
    g = ICOSAHEDRAL[0]
    return SubstitutionRule(
        n_prototiles=2, inflation=Mat3.identity(),
        dissections=(
            (PositionedTile(0, origin, g),
             PositionedTile(0, origin, g),
             PositionedTile(1, origin, g)),
            (PositionedTile(0, origin, g),
             PositionedTile(1, origin, g)),
        ),
    )


class TestFeasibilityUpperBound:
    """Per Claude (web) Q7a/b 2026-04-29 ruling: closed-form check
    for positive-integer eigenvalue, short-circuiting the entire
    face-merge enumeration when no such eigenvalue exists.
    """

    def test_danzer_returns_false(self, paolini_rule) -> None:
        # M_ABCK eigenvalues = {φ³, φ, -1/φ, -1/φ³} — none is a
        # positive integer.
        assert feasibility_upper_bound(paolini_rule) is False

    def test_phi_squared_rule_returns_false(self) -> None:
        # M = [[2,1],[1,1]] has eigenvalues φ² and 1/φ² — also no
        # positive integer.
        assert feasibility_upper_bound(_phi_squared_rule()) is False

    def test_trivial_pf2_rule_returns_true(self) -> None:
        # M = [[2]] has eigenvalue 2 — positive integer.
        assert feasibility_upper_bound(_trivial_pf2_rule()) is True


class TestEnumerateFaceMergeCompositions:
    """Phase 1.6 k-tuple enumeration: brute-force confirmation of
    the closed-form no-go.
    """

    def test_danzer_yields_no_survivors_at_k_max_4(
        self, paolini_rule,
    ) -> None:
        # The Phase 1.6 closing result: every composition fails
        # the eigenvector check on M_ABCK. (Short-circuited via
        # feasibility_upper_bound; result independent of k_max.)
        survivors = enumerate_face_merge_compositions(
            paolini_rule, k_max=4,
        )
        assert survivors == ()

    def test_danzer_yields_no_survivors_at_k_max_8(
        self, paolini_rule,
    ) -> None:
        # The closed-form argument is k-independent; confirm by
        # raising k_max well beyond ABCK's alphabet.
        survivors = enumerate_face_merge_compositions(
            paolini_rule, k_max=8,
        )
        assert survivors == ()

    def test_trivial_pf2_yields_every_composition(self) -> None:
        # M = [[2]]: every (k,) composition for k ∈ {1, 2, 3, 4}
        # satisfies σ(k) = 2·k.
        survivors = enumerate_face_merge_compositions(
            _trivial_pf2_rule(), k_max=4,
        )
        assert survivors == ((1,), (2,), (3,), (4,))

    def test_rejects_zero_k_max(self, paolini_rule) -> None:
        with pytest.raises(ValueError, match="k_max must be"):
            enumerate_face_merge_compositions(paolini_rule, k_max=0)


class TestPhase1_6_NoGo:
    """The single integration test that records the Phase 1.6 closing
    result for the project's research record.
    """

    def test_face_merge_track_closed_for_abck(self, paolini_rule) -> None:
        # The closed-form Track A face-merge no-go, per Claude (web)
        # Q7a 2026-04-29: M_ABCK has no positive-integer eigenvalue,
        # therefore no face-merge composition of any size yields a
        # 1-prototile substitution rule with PF in ℤ[φ]. This test
        # is the citable artifact.
        assert feasibility_upper_bound(paolini_rule) is False
        # Higher-k enumerations all return empty.
        for k in (2, 3, 4, 5, 6, 8):
            assert enumerate_face_merge_compositions(
                paolini_rule, k_max=k,
            ) == ()
