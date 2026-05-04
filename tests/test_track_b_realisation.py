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
    Inconclusive,
    NoRealisation,
    Realised,
    SearchProgress,
    enumerate_primitive_matrices,
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

    def test_inconclusive_reason_mentions_csp_pending(self) -> None:
        m = np.array([[0, 1], [1, 4]])
        result = realise(m, ZPhi(1, 2))
        assert isinstance(result, Inconclusive)
        assert "CSP" in result.reason

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
