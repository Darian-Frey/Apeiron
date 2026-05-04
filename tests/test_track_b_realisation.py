"""Tests for apeiron.track_b.realisation — the realisation CSP API
+ Fibonacci oracle witness.

Per Claude (web)'s Q8 ruling 2026-04-29: the realise() function
returns one of Realised | NoRealisation | Inconclusive. Current
status: API + Fibonacci oracle as a manual witness; algorithmic CSP
for general inputs is a follow-up commit.
"""

from __future__ import annotations

import numpy as np
import pytest

from apeiron.polyhedron import Polyhedron
from apeiron.symmetry import Rotation, Vec3
from apeiron.track_b import (
    ChildPlacement,
    FaceMatchEdge,
    Inconclusive,
    NoRealisation,
    Realised,
    SearchProgress,
    enumerate_primitive_matrices,
    propagate_translations_along_tree,
    realise,
    translation_offset_from_face_match,
)
from apeiron.zphi import ZPhi


class TestResultTypes:
    """Validation of the Realised / NoRealisation / Inconclusive
    dataclasses and the ChildPlacement / SearchProgress helpers.
    """

    def test_realised_is_frozen(self) -> None:
        z = ZPhi(0, 0)
        v = Vec3(z, z, z)
        cp = ChildPlacement(0, v, Rotation.identity())
        proto = Polyhedron.from_raw(
            [
                Vec3(z, z, z),
                Vec3(ZPhi(2, 0), z, z),
                Vec3(z, ZPhi(2, 0), z),
                Vec3(z, z, ZPhi(2, 0)),
            ],
            [(0, 1, 2), (0, 3, 1), (0, 2, 3), (1, 3, 2)],
        )
        r = Realised(
            prototile_shapes=(proto,),
            children_per_parent=((cp,),),
            fraction_searched=ZPhi(1, 0),
        )
        with pytest.raises(Exception):
            r.fraction_searched = ZPhi(0, 0)  # type: ignore[misc]

    def test_norealisation_default_fraction_is_one(self) -> None:
        # NoRealisation is "provably no realisation in the bounded
        # space" so fraction_searched must be 1.
        nr = NoRealisation(reason="empty space")
        assert nr.fraction_searched == ZPhi(1, 0)

    def test_inconclusive_carries_partial_state(self) -> None:
        z = ZPhi(0, 0)
        progress = SearchProgress(
            fraction_searched=ZPhi(1, 0) - ZPhi(1, 0),  # demo
            realised_partial=(),
        )
        inc = Inconclusive(
            fraction_searched=z,
            partial=progress,
            reason="budget exhausted",
        )
        assert inc.partial is progress

    def test_child_placement_carries_pose(self) -> None:
        z = ZPhi(0, 0)
        cp = ChildPlacement(
            prototile_index=2,
            translation=Vec3(z, z, z),
            rotation=Rotation.identity(),
        )
        assert cp.prototile_index == 2
        assert cp.rotation == Rotation.identity()


class TestFibonacciOracle:
    """Q8f Fibonacci oracle: M=[[0,1],[1,1]], PF=φ. The realisation
    is constructible by hand (boxes packed along x). realise() must
    return Realised with a complete witness.
    """

    def test_fibonacci_returns_realised(self) -> None:
        m = np.array([[0, 1], [1, 1]])
        result = realise(m, ZPhi(0, 1))
        assert isinstance(result, Realised)

    def test_fibonacci_has_two_prototile_shapes(self) -> None:
        m = np.array([[0, 1], [1, 1]])
        result = realise(m, ZPhi(0, 1))
        assert isinstance(result, Realised)
        assert len(result.prototile_shapes) == 2

    def test_fibonacci_child_count_matches_matrix(self) -> None:
        # Column 0 = σ(P_0) = 0·P_0 + 1·P_1 → 1 child.
        # Column 1 = σ(P_1) = 1·P_0 + 1·P_1 → 2 children.
        m = np.array([[0, 1], [1, 1]])
        result = realise(m, ZPhi(0, 1))
        assert isinstance(result, Realised)
        assert len(result.children_per_parent[0]) == 1
        assert len(result.children_per_parent[1]) == 2

    def test_fibonacci_child_types_match_matrix(self) -> None:
        m = np.array([[0, 1], [1, 1]])
        result = realise(m, ZPhi(0, 1))
        assert isinstance(result, Realised)
        # σ(P_0)'s only child is P_1 (since M[1,0] = 1).
        assert result.children_per_parent[0][0].prototile_index == 1
        # σ(P_1)'s children: one of each type. Multiset check.
        types_in_p1 = sorted(
            cp.prototile_index for cp in result.children_per_parent[1]
        )
        assert types_in_p1 == [0, 1]

    def test_fibonacci_fraction_searched_is_one(self) -> None:
        # Manual witness covers the entire (trivial) search space.
        m = np.array([[0, 1], [1, 1]])
        result = realise(m, ZPhi(0, 1))
        assert isinstance(result, Realised)
        assert result.fraction_searched == ZPhi(1, 0)

    def test_fibonacci_p1_box_is_phi_long(self) -> None:
        # P_1 is the "long" Fibonacci tile; its x-extent is φ.
        # Vertices are ×2-stored, so max x-coordinate should be 2φ
        # (= ZPhi(0, 2) in stored form).
        m = np.array([[0, 1], [1, 1]])
        result = realise(m, ZPhi(0, 1))
        assert isinstance(result, Realised)
        p1 = result.prototile_shapes[1]
        max_x = max(v.x for v in p1.vertices)
        assert max_x == ZPhi(0, 2)


class TestPlaceholderForGeneralInputs:
    """All non-Fibonacci inputs currently return Inconclusive with
    fraction_searched = 0. This is the expected status until the
    algorithmic CSP lands.
    """

    def test_phi_cubed_candidate_returns_inconclusive(self) -> None:
        # The 5 PF=φ³ candidates from the Track B survey all currently
        # return Inconclusive — that's the expected placeholder
        # behavior, not a result.
        m = np.array([[0, 1], [1, 4]])
        result = realise(m, ZPhi(1, 2))
        assert isinstance(result, Inconclusive)
        assert result.fraction_searched == ZPhi(0, 0)

    def test_inconclusive_reason_mentions_prototile_shapes(self) -> None:
        # When called without prototile_shapes, the realise() call
        # returns Inconclusive whose reason directs the caller to
        # supply shapes per Q8b.
        m = np.array([[0, 1], [1, 4]])
        result = realise(m, ZPhi(1, 2))
        assert isinstance(result, Inconclusive)
        assert "prototile_shapes" in result.reason

    def test_penrose_p3_returns_inconclusive(self) -> None:
        # Q8f's stress-test oracle (P3, PF=φ²) currently returns
        # Inconclusive — to be flipped to Realised when the CSP lands.
        m = np.array([[1, 1], [1, 2]])
        result = realise(m, ZPhi(1, 1))
        assert isinstance(result, Inconclusive)

    def test_all_5_phi_cubed_candidates_inconclusive(self) -> None:
        # Sweep across the Track B PF=φ³ survey; every result is
        # Inconclusive at this commit.
        count = 0
        for matrix in enumerate_primitive_matrices(
            2, ZPhi(1, 2), max_entry=5,
        ):
            result = realise(matrix, ZPhi(1, 2))
            assert isinstance(result, Inconclusive)
            count += 1
        assert count == 5


class TestTranslationOffsetFromFaceMatch:
    """Per Q8a 2026-04-29: per-edge constraint solver. Given fixed
    rotations and a triangular face match between two prototiles,
    recover ``t_b − t_a`` exactly in ZPhi.
    """

    def _unit_tetrahedron_vertices(self) -> tuple[Vec3, Vec3, Vec3, Vec3]:
        # ×2-stored unit tetrahedron with one vertex at origin and the
        # other three on the positive coordinate axes at (1, 0, 0),
        # (0, 1, 0), (0, 0, 1). All in pure ℤ.
        z = ZPhi(0, 0)
        two = ZPhi(2, 0)
        return (
            Vec3(z, z, z),
            Vec3(two, z, z),
            Vec3(z, two, z),
            Vec3(z, z, two),
        )

    def test_identity_gluing_yields_zero_offset(self) -> None:
        # Two identical tetrahedra placed at identity / identity-rot,
        # sharing the face opposite vertex 0. Trivial offset (0, 0, 0).
        verts = self._unit_tetrahedron_vertices()
        offset = translation_offset_from_face_match(
            rotation_a=Rotation.identity(),
            rotation_b=Rotation.identity(),
            face_indices_a=(1, 2, 3),
            face_indices_b=(1, 2, 3),
            vertices_a=verts, vertices_b=verts,
        )
        assert offset == Vec3(ZPhi(0, 0), ZPhi(0, 0), ZPhi(0, 0))

    def test_incongruent_faces_return_none(self) -> None:
        # Tetrahedron face (0, 1, 2) is the right triangle in z=0
        # plane (legs 2 along x, y). Face (0, 1, 3) is a different
        # right triangle. Not congruent — no permutation matches.
        verts = self._unit_tetrahedron_vertices()
        offset = translation_offset_from_face_match(
            rotation_a=Rotation.identity(),
            rotation_b=Rotation.identity(),
            face_indices_a=(0, 1, 2),
            face_indices_b=(0, 1, 3),
            vertices_a=verts, vertices_b=verts,
        )
        assert offset is None

    def test_recovers_paolini_sigma_k_offset(self) -> None:
        # Real-world cross-check: σ(K)'s B-child + K-child face-
        # adjacency in the canonical Paolini Danzer rule. The helper
        # recovers exactly cj.translation − ci.translation.
        import json
        from pathlib import Path
        from apeiron.deformation import face_adjacent_pairs, placed_vertices
        from apeiron.symmetry import ICOSAHEDRAL, ImproperRotation, Mat3
        from apeiron.substitution import PositionedTile, SubstitutionRule
        from apeiron.util import load_candidate

        DANZER = Path("candidates/danzer")
        prototiles = [load_candidate(DANZER / f"{x}.json") for x in "ABCK"]
        z = ZPhi(0, 0)
        phi = ZPhi(0, 1)
        inflation = Mat3(
            Vec3(phi, z, z), Vec3(z, phi, z), Vec3(z, z, phi),
        )
        data = json.loads((DANZER / "paolini_dissection.json").read_text())
        type_to_idx = {"A": 0, "B": 1, "C": 2, "K": 3}
        by_parent: dict[str, list] = {p: [] for p in "ABCK"}
        for r in data["children"]:
            by_parent[r["parent"]].append(r)
        dissections = []
        for parent in "ABCK":
            children = []
            for r in sorted(
                by_parent[parent], key=lambda r: r["child_index"],
            ):
                t = r["translation_x2"]
                translation = Vec3(
                    ZPhi(*t[0]), ZPhi(*t[1]), ZPhi(*t[2]),
                )
                proper = ICOSAHEDRAL[r["icosahedral_index"]]
                rotation = (
                    proper if r["is_proper"] else ImproperRotation(proper)
                )
                children.append(PositionedTile(
                    prototile_index=type_to_idx[r["child_type"]],
                    translation=translation, rotation=rotation,
                ))
            dissections.append(tuple(children))
        rule = SubstitutionRule(
            n_prototiles=4, inflation=inflation,
            dissections=tuple(dissections),
        )

        # σ(K)'s only face-adjacent pair: i=0 (B-child, ImproperRotation),
        # j=1 (K-child, proper Rotation).
        pair = face_adjacent_pairs(rule, 3, prototiles)[0]
        ci = rule.dissections[3][pair.i]
        cj = rule.dissections[3][pair.j]
        proto_a = prototiles[ci.prototile_index]
        proto_b = prototiles[cj.prototile_index]
        placed_a = placed_vertices(proto_a, (ci.translation, ci.rotation))
        placed_b = placed_vertices(proto_b, (cj.translation, cj.rotation))
        face_a_idx = tuple(
            i for i, v in enumerate(placed_a) if v in pair.shared_vertices
        )
        face_b_idx = tuple(
            j for j, v in enumerate(placed_b) if v in pair.shared_vertices
        )
        assert len(face_a_idx) == 3 and len(face_b_idx) == 3

        offset = translation_offset_from_face_match(
            ci.rotation, cj.rotation,
            face_a_idx, face_b_idx,
            proto_a.vertices, proto_b.vertices,
        )
        expected = cj.translation - ci.translation
        assert offset == expected

    def test_rejects_non_triangular_face(self) -> None:
        verts = self._unit_tetrahedron_vertices()
        with pytest.raises(ValueError, match="exactly 3"):
            translation_offset_from_face_match(
                rotation_a=Rotation.identity(),
                rotation_b=Rotation.identity(),
                face_indices_a=(0, 1, 2, 3),  # type: ignore[arg-type]
                face_indices_b=(1, 2, 3),
                vertices_a=verts, vertices_b=verts,
            )


class TestFaceMatchEdge:
    def test_rejects_self_loop(self) -> None:
        with pytest.raises(ValueError, match="differ from"):
            FaceMatchEdge(
                parent=0, new=0,
                face_indices_parent=(0, 1, 2),
                face_indices_new=(0, 1, 2),
            )

    def test_rejects_non_triangular_parent_face(self) -> None:
        with pytest.raises(ValueError, match="3 vertex indices on parent"):
            FaceMatchEdge(
                parent=0, new=1,
                face_indices_parent=(0, 1, 2, 3),  # type: ignore[arg-type]
                face_indices_new=(0, 1, 2),
            )

    def test_rejects_non_triangular_new_face(self) -> None:
        with pytest.raises(ValueError, match="3 vertex indices on new"):
            FaceMatchEdge(
                parent=0, new=1,
                face_indices_parent=(0, 1, 2),
                face_indices_new=(0, 1),  # type: ignore[arg-type]
            )


class TestPropagateTranslationsAlongTree:
    """Per Q8a 2026-04-29: tree-DFS translation propagation. Composes
    translation_offset_from_face_match over a fixed tree topology.
    """

    def _unit_tetrahedron(self) -> tuple[Vec3, ...]:
        z = ZPhi(0, 0)
        two = ZPhi(2, 0)
        return (
            Vec3(z, z, z),
            Vec3(two, z, z),
            Vec3(z, two, z),
            Vec3(z, z, two),
        )

    def test_two_child_identity_gluing(self) -> None:
        # Two identical tetrahedra glued at face (1, 2, 3) ↔ (1, 2, 3),
        # both at identity rotation. Translations: child 0 at origin,
        # child 1 also at origin (trivial overlap, OK for unit test).
        verts = self._unit_tetrahedron()
        edges = [FaceMatchEdge(
            parent=0, new=1,
            face_indices_parent=(1, 2, 3),
            face_indices_new=(1, 2, 3),
        )]
        result = propagate_translations_along_tree(
            rotations=[Rotation.identity(), Rotation.identity()],
            edges=edges,
            prototile_indices=[0, 0],
            prototile_vertices=[verts],
        )
        assert result is not None
        assert result == (
            Vec3(ZPhi(0, 0), ZPhi(0, 0), ZPhi(0, 0)),
            Vec3(ZPhi(0, 0), ZPhi(0, 0), ZPhi(0, 0)),
        )

    def test_child_zero_at_origin(self) -> None:
        # Convention: child 0 always at (0, 0, 0).
        verts = self._unit_tetrahedron()
        edges = [FaceMatchEdge(
            parent=0, new=1,
            face_indices_parent=(1, 2, 3),
            face_indices_new=(1, 2, 3),
        )]
        result = propagate_translations_along_tree(
            rotations=[Rotation.identity(), Rotation.identity()],
            edges=edges,
            prototile_indices=[0, 0],
            prototile_vertices=[verts],
        )
        assert result is not None
        assert result[0] == Vec3(ZPhi(0, 0), ZPhi(0, 0), ZPhi(0, 0))

    def test_incompatible_face_returns_none(self) -> None:
        # Incongruent face pair → propagate returns None.
        verts = self._unit_tetrahedron()
        edges = [FaceMatchEdge(
            parent=0, new=1,
            face_indices_parent=(0, 1, 2),
            face_indices_new=(0, 1, 3),
        )]
        result = propagate_translations_along_tree(
            rotations=[Rotation.identity(), Rotation.identity()],
            edges=edges,
            prototile_indices=[0, 0],
            prototile_vertices=[verts],
        )
        assert result is None

    def test_real_paolini_sigma_k(self) -> None:
        # Cross-check on Paolini's σ(K) dissection. Tree: B-child → K-child.
        # Propagated translations match cj.translation − ci.translation.
        import json
        from pathlib import Path
        from apeiron.deformation import face_adjacent_pairs, placed_vertices
        from apeiron.symmetry import ICOSAHEDRAL, ImproperRotation, Mat3
        from apeiron.substitution import PositionedTile, SubstitutionRule
        from apeiron.util import load_candidate

        DANZER = Path("candidates/danzer")
        prototiles = [load_candidate(DANZER / f"{x}.json") for x in "ABCK"]
        z = ZPhi(0, 0)
        phi = ZPhi(0, 1)
        inflation = Mat3(
            Vec3(phi, z, z), Vec3(z, phi, z), Vec3(z, z, phi),
        )
        data = json.loads((DANZER / "paolini_dissection.json").read_text())
        type_to_idx = {"A": 0, "B": 1, "C": 2, "K": 3}
        by_parent: dict[str, list] = {p: [] for p in "ABCK"}
        for r in data["children"]:
            by_parent[r["parent"]].append(r)
        dissections = []
        for parent in "ABCK":
            children = []
            for r in sorted(
                by_parent[parent], key=lambda r: r["child_index"],
            ):
                t = r["translation_x2"]
                translation = Vec3(
                    ZPhi(*t[0]), ZPhi(*t[1]), ZPhi(*t[2]),
                )
                proper = ICOSAHEDRAL[r["icosahedral_index"]]
                rotation = (
                    proper if r["is_proper"] else ImproperRotation(proper)
                )
                children.append(PositionedTile(
                    prototile_index=type_to_idx[r["child_type"]],
                    translation=translation, rotation=rotation,
                ))
            dissections.append(tuple(children))
        rule = SubstitutionRule(
            n_prototiles=4, inflation=inflation,
            dissections=tuple(dissections),
        )

        pair = face_adjacent_pairs(rule, 3, prototiles)[0]
        ci = rule.dissections[3][pair.i]
        cj = rule.dissections[3][pair.j]
        placed_a = placed_vertices(
            prototiles[ci.prototile_index], (ci.translation, ci.rotation),
        )
        placed_b = placed_vertices(
            prototiles[cj.prototile_index], (cj.translation, cj.rotation),
        )
        face_a = tuple(
            i for i, v in enumerate(placed_a) if v in pair.shared_vertices
        )
        face_b = tuple(
            j for j, v in enumerate(placed_b) if v in pair.shared_vertices
        )

        edges = [FaceMatchEdge(
            parent=0, new=1,
            face_indices_parent=face_a,
            face_indices_new=face_b,
        )]
        result = propagate_translations_along_tree(
            rotations=[ci.rotation, cj.rotation],
            edges=edges,
            prototile_indices=[ci.prototile_index, cj.prototile_index],
            prototile_vertices=[tuple(p.vertices) for p in prototiles],
        )
        assert result is not None
        assert result[0] == Vec3(ZPhi(0, 0), ZPhi(0, 0), ZPhi(0, 0))
        # Relative offset should equal cj.translation − ci.translation.
        assert result[1] == cj.translation - ci.translation

    def test_rejects_length_mismatch(self) -> None:
        verts = self._unit_tetrahedron()
        with pytest.raises(ValueError, match="prototile_indices length"):
            propagate_translations_along_tree(
                rotations=[Rotation.identity(), Rotation.identity()],
                edges=[],
                prototile_indices=[0],  # length 1 vs rotations length 2
                prototile_vertices=[verts],
            )

    def test_rejects_wrong_edge_count(self) -> None:
        verts = self._unit_tetrahedron()
        # 3 children but only 1 edge — tree needs 2 edges.
        with pytest.raises(ValueError, match="needs"):
            propagate_translations_along_tree(
                rotations=[
                    Rotation.identity(),
                    Rotation.identity(),
                    Rotation.identity(),
                ],
                edges=[FaceMatchEdge(
                    parent=0, new=1,
                    face_indices_parent=(1, 2, 3),
                    face_indices_new=(1, 2, 3),
                )],
                prototile_indices=[0, 0, 0],
                prototile_vertices=[verts],
            )

    def test_rejects_disordered_edges(self) -> None:
        # Edges must be in DFS order — first edge's parent must be 0.
        verts = self._unit_tetrahedron()
        with pytest.raises(ValueError, match="not yet"):
            propagate_translations_along_tree(
                rotations=[Rotation.identity()] * 3,
                edges=[
                    # Edge with parent=2 before 2 has been placed.
                    FaceMatchEdge(
                        parent=2, new=1,
                        face_indices_parent=(1, 2, 3),
                        face_indices_new=(1, 2, 3),
                    ),
                    FaceMatchEdge(
                        parent=0, new=2,
                        face_indices_parent=(1, 2, 3),
                        face_indices_new=(1, 2, 3),
                    ),
                ],
                prototile_indices=[0, 0, 0],
                prototile_vertices=[verts],
            )


class TestSearchLoopIntegration:
    """End-to-end: realise() now runs the search loop for non-Fibonacci
    inputs when prototile_shapes is provided. **Validation is bounding-
    box only at this iteration**, so Realised is a necessary-condition
    witness, not sufficient — overlap detection is the next refinement.
    """

    def _unit_tetrahedron(self) -> Polyhedron:
        z = ZPhi(0, 0)
        two = ZPhi(2, 0)
        return Polyhedron.from_raw(
            [
                Vec3(z, z, z),
                Vec3(two, z, z),
                Vec3(z, two, z),
                Vec3(z, z, two),
            ],
            [(0, 1, 2), (0, 3, 1), (0, 2, 3), (1, 3, 2)],
        )

    def test_no_prototile_shapes_returns_inconclusive(self) -> None:
        m = np.array([[0, 1], [1, 4]])
        result = realise(m, ZPhi(1, 2))
        assert isinstance(result, Inconclusive)
        assert "prototile_shapes" in result.reason

    def test_volume_mismatch_rejected_before_search(self) -> None:
        # M=[[0,1],[1,4]] with identical-volume tetrahedra: σ(P_0)
        # column = [0, 1] → 1 P_1 child. Children's vol sum = vol(P_1)
        # = vol(P_0); required pf · vol(P_0) = (1+2φ) · vol(P_0).
        # 1 ≠ 1+2φ → volume check rejects σ(P_0) immediately, before
        # the k=5 σ(P_1) is even considered.
        m = np.array([[0, 1], [1, 4]])
        proto = self._unit_tetrahedron()
        result = realise(
            m, ZPhi(1, 2),
            prototile_shapes=(proto, proto),
            max_search_seconds=10,
        )
        assert isinstance(result, NoRealisation)

    def test_penrose_p3_with_identical_protos_rejected_by_volume(
        self,
    ) -> None:
        # M=[[1,1],[1,2]] (Penrose P3): σ(P_0) sum-of-vols = 2·V if
        # P_0 and P_1 are identical-volume; required pf·V where
        # pf=φ². 2 ≠ φ² so the volume-sum check rejects.
        m = np.array([[1, 1], [1, 2]])
        proto = self._unit_tetrahedron()
        result = realise(
            m, ZPhi(1, 1),
            prototile_shapes=(proto, proto),
            max_search_seconds=10,
        )
        assert isinstance(result, NoRealisation)

    def test_trivial_pf2_volume_passes_overlap_rejects(self) -> None:
        # M=[[2]] (PF=2) with identical-shape children: volume-sum
        # passes (2V = pf·V), but pairwise SAT overlap-detection
        # rejects every rotation assignment (two identical tetrahedra
        # at the origin always overlap). Result: NoRealisation —
        # confirming SAT correctly catches the trivial overlap that
        # the previous bounding-box-only iteration missed.
        m = np.array([[2]])
        proto = self._unit_tetrahedron()
        result = realise(
            m, ZPhi(2, 0),
            prototile_shapes=(proto,),
            max_search_seconds=10,
        )
        assert isinstance(result, NoRealisation)


class TestSatOverlap:
    """Per-pair separating-axis-theorem tests on tetrahedral
    prototiles."""

    def _unit_tet_verts(self) -> tuple[Vec3, ...]:
        z = ZPhi(0, 0)
        two = ZPhi(2, 0)
        return (
            Vec3(z, z, z),
            Vec3(two, z, z),
            Vec3(z, two, z),
            Vec3(z, z, two),
        )

    def _unit_tet_faces(self) -> tuple[tuple[int, int, int], ...]:
        return ((0, 1, 2), (0, 3, 1), (0, 2, 3), (1, 3, 2))

    def test_identical_tets_at_origin_overlap(self) -> None:
        from apeiron.track_b.realisation import (
            _convex_polyhedra_interior_disjoint,
        )
        v = self._unit_tet_verts()
        f = self._unit_tet_faces()
        assert _convex_polyhedra_interior_disjoint(v, f, v, f) is False

    def test_translated_apart_are_disjoint(self) -> None:
        from apeiron.track_b.realisation import (
            _convex_polyhedra_interior_disjoint,
        )
        v_a = self._unit_tet_verts()
        z = ZPhi(0, 0)
        offset = Vec3(ZPhi(10, 0), z, z)
        v_b = tuple(p + offset for p in v_a)
        f = self._unit_tet_faces()
        assert _convex_polyhedra_interior_disjoint(v_a, f, v_b, f) is True

    def test_face_tangent_is_disjoint(self) -> None:
        # Two unit tetrahedra glued at the face (1, 2, 3) (face
        # opposite the origin vertex). After reflecting B through
        # this plane, their interiors are disjoint and they touch
        # only on the shared face. SAT should report disjoint.
        from apeiron.track_b.realisation import (
            _convex_polyhedra_interior_disjoint,
        )
        v_a = self._unit_tet_verts()
        # Reflect B's vertex 0 (= origin) through the (1, 2, 3) face
        # plane x + y + z = 2. Reflection of (0,0,0) is (4/3, 4/3, 4/3)
        # in real (= ZPhi(4, 0) / 3, not in ZPhi). Skip to a simpler
        # gluing: place B's vertex 0 outside A, sharing face (1, 2, 3).
        # Concretely, B is the mirror image of A across the plane
        # through (1, 2, 3). For a regular-corner tet (1,2,3 at
        # x+y+z=2 plane), B's vertex 0 lands at (4/3, 4/3, 4/3) — not
        # ZPhi.  Instead build a tangent case via SAT directly: A and
        # B both contain a shared face but their interiors are on
        # opposite sides. Construct B as A's reflection across z=0:
        z = ZPhi(0, 0)
        two = ZPhi(2, 0)
        v_b = (
            Vec3(z, z, z),
            Vec3(two, z, z),
            Vec3(z, two, z),
            Vec3(z, z, ZPhi(-2, 0)),  # mirror vertex 3 below z=0
        )
        f = self._unit_tet_faces()
        # A and B share the face (0, 1, 2) (the z=0 triangle); their
        # interiors are on opposite sides (z > 0 for A, z < 0 for B).
        assert _convex_polyhedra_interior_disjoint(
            v_a, f, v_b, f,
        ) is True


class TestInputValidation:
    def test_rejects_non_square_matrix(self) -> None:
        m = np.array([[1, 0, 0], [0, 1, 0]])
        with pytest.raises(ValueError, match="square"):
            realise(m, ZPhi(0, 1))

    def test_rejects_mismatched_prototile_count(self) -> None:
        m = np.array([[0, 1], [1, 1]])
        z = ZPhi(0, 0)
        proto = Polyhedron.from_raw(
            [
                Vec3(z, z, z),
                Vec3(ZPhi(2, 0), z, z),
                Vec3(z, ZPhi(2, 0), z),
                Vec3(z, z, ZPhi(2, 0)),
            ],
            [(0, 1, 2), (0, 3, 1), (0, 2, 3), (1, 3, 2)],
        )
        with pytest.raises(ValueError, match="prototile_shapes length"):
            realise(m, ZPhi(0, 1), prototile_shapes=(proto,))
