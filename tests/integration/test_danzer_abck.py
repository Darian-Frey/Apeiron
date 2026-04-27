"""Danzer ABCK 4-tile baseline.

Track A (CLAUDE.md §6.1) starts here: encode Danzer's published
4-tile aperiodic set, run it through the pipeline, then perturb
toward a single-tile candidate via deformation search.

Sub-commit 27A: load and structurally verify all four prototiles
(A, B, C, K). Each is a tetrahedron with 4 vertices, 6 edges, 4
triangular faces, Euler characteristic 2. K's vertices include a
class-IV half-integer position (per Frettlöh's note) so its
``scale_denom`` is 2 — this is the first wild-input exercise of the
loader's half-integer path.

The substitution rule itself (matrix + dissection) lands in 27B.
The fourth-pillar stub lands in 27C. The full pipeline run on the
4-tile set lands in 27D.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from apeiron.polyhedron import Polyhedron
from apeiron.util import load_candidate

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
        import json
        expected = {"A": 0, "B": 1, "C": 2, "K": 3}
        for letter, idx in expected.items():
            with (_DANZER_DIR / f"{letter}.json").open() as f:
                payload = json.load(f)
            assert payload["prototile_index"] == idx
