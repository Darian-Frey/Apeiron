"""Tests for apeiron.track_b.taxonomy — H₃-compatible tetrahedral
classification per Claude (web)'s Q9a 2026-04-29 ruling.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from apeiron.symmetry import Vec3
from apeiron.track_b.taxonomy import (
    H3_AXES,
    build_h3_tetrahedra,
    is_h3_compatible,
)
from apeiron.track_b.taxonomy import _similarity_class_key
from apeiron.util import load_candidate
from apeiron.zphi import ZPhi


_DANZER = Path("candidates/danzer")


class TestH3Axes:
    """The 15 H₃-mirror normals = icosahedron edge midpoint
    directions, derived from first principles (per Q9a 2026-04-29).
    """

    def test_count_is_fifteen(self) -> None:
        # 30 edge midpoints quotiented by ± gives 15 axes.
        assert len(H3_AXES) == 15

    def test_axes_are_zphi(self) -> None:
        for axis in H3_AXES:
            assert isinstance(axis, Vec3)
            for comp in (axis.x, axis.y, axis.z):
                assert isinstance(comp, ZPhi)
                assert isinstance(comp.a, int)
                assert isinstance(comp.b, int)

    def test_axes_pairwise_distinct(self) -> None:
        # No duplicate axes after the ± quotient.
        assert len(set(H3_AXES)) == 15


class TestIsH3Compatible:
    """Per-face normal direction check on a Polyhedron."""

    @pytest.mark.parametrize("letter", ["A", "B", "C", "K"])
    def test_abck_prototiles_are_h3_compatible(self, letter) -> None:
        # All 4 ABCK prototiles have face normals parallel to H₃ axes.
        # This is the foundational test — Frettlöh's classification
        # places ABCK in the H₃-compatible family by construction.
        proto = load_candidate(_DANZER / f"{letter}.json")
        assert is_h3_compatible(proto)


class TestBuildH3Tetrahedra:
    """Enumerate H₃-compatible tetrahedra from a vertex pool, modulo
    similarity.

    Status (2026-04-29): default ABCK pool (10 vertices from Paolini)
    yields 9 similarity classes. Frettlöh's reported count of 15
    awaits a richer vertex pool — research follow-up. The 4 ABCK
    shapes are confirmed among the 9.
    """

    def test_default_pool_yields_nine_classes(self) -> None:
        # Empirically, the 10-vertex Paolini pool yields 9 H₃-
        # compatible similarity classes (4 ABCK + 5 additional
        # tetrahedra emerging from other 4-subsets).
        tets = build_h3_tetrahedra()
        assert len(tets) == 9

    def test_all_four_abck_shapes_are_present(self) -> None:
        # Every ABCK prototile must appear (up to similarity) among
        # the enumerated classes. Without this, the enumerator is
        # broken.
        tets = build_h3_tetrahedra()
        tet_keys = {_similarity_class_key(p.vertices) for p in tets}
        for letter in "ABCK":
            proto = load_candidate(_DANZER / f"{letter}.json")
            assert _similarity_class_key(proto.vertices) in tet_keys, (
                f"ABCK prototile {letter} missing from enumerated H₃ "
                "tetrahedra."
            )

    def test_each_returned_polyhedron_is_h3_compatible(self) -> None:
        # The enumeration must self-validate.
        for tet in build_h3_tetrahedra():
            assert is_h3_compatible(tet)
            assert len(tet.vertices) == 4
            assert len(tet.faces) == 4

    def test_classes_are_pairwise_dissimilar(self) -> None:
        # No two classes share a similarity key.
        tets = build_h3_tetrahedra()
        keys = [_similarity_class_key(p.vertices) for p in tets]
        assert len(set(keys)) == len(keys)
