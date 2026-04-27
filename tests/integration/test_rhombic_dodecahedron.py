"""Rhombic dodecahedron acceptance test for ``corona_1``.

Per CLAUDE.md §5.1 acceptance criterion 5: "rhombic dodecahedron has
a unique ``corona_1``." Together with the cube test, this is the
second sanity-check oracle for the BFS engine — the RD exercises a
case the cube doesn't:

- Two distinct vertex types in the tiling. The RD has 8 three-valent
  (cube-corner) vertices where 4 tiles meet and 6 four-valent
  (octahedral) vertices where 6 tiles meet. The ``corona_1`` API's
  per-feature ``expected_vertex_count`` callable is exercised here
  for the first time.
- Three tiles per edge instead of four (the RD dihedral is 120°).

Cross-module per the ``≥-2-core-imports`` rule
(``tests/integration/README.md``): touches ``corona``, ``util``, and
``polyhedron``.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path

import pytest

from apeiron.corona import (
    Edge,
    Vertex,
    corona_1,
    has_interior_overlap,
    incidence_defect,
)
from apeiron.polyhedron import Polyhedron
from apeiron.util import load_candidate

_RD_PATH = Path("candidates/rhombic_dodecahedron.json")


def _vertex_degree_map(P: Polyhedron) -> dict[int, int]:
    """Number of faces incident to each vertex.

    For the RD, returns 3 at cube-corner vertices and 4 at octahedral
    vertices.
    """
    degree: Counter[int] = Counter()
    for face in P.faces:
        for idx in face:
            degree[idx] += 1
    return dict(degree)


@pytest.fixture(scope="module")
def rd_corona() -> tuple[Polyhedron, tuple]:
    """Load the rhombic dodecahedron and compute its first corona.

    Module-scoped so the ~6 s ``corona_1`` invocation runs once for
    the whole test class. Expected counts: 3 tiles per edge (uniform);
    4 tiles per 3-valent vertex, 6 per 4-valent vertex.
    """
    rd = load_candidate(_RD_PATH)
    degree = _vertex_degree_map(rd)

    def expected_vertex(v: Vertex) -> int:
        return 4 if degree[v.index] == 3 else 6

    configs = corona_1(
        rd,
        expected_edge_count=3,
        expected_vertex_count=expected_vertex,
    )
    return rd, configs


# -- structural acceptance --------------------------------------------


class TestRhombicDodecahedronCorona:
    def test_returns_exactly_one_canonical_configuration(
        self, rd_corona,
    ) -> None:
        _, configs = rd_corona
        assert len(configs) == 1

    def test_unique_corona_has_18_neighbours(self, rd_corona) -> None:
        _, configs = rd_corona
        # 12 face-sharers (one per RD face) + 6 vertex-sharers at the
        # 4-valent (octahedral) positions = 18.
        assert len(configs[0].neighbours) == 18

    def test_no_interior_overlap(self, rd_corona) -> None:
        _, configs = rd_corona
        assert not has_interior_overlap(configs[0])

    def test_every_edge_has_incidence_count_three(self, rd_corona) -> None:
        _, configs = rd_corona
        cfg = configs[0]
        # Build the edge set from the canonicalised central's faces.
        edges: set[tuple[int, int]] = set()
        for face in cfg.central.faces:
            n = len(face)
            for idx in range(n):
                u, v = face[idx], face[(idx + 1) % n]
                edges.add((min(u, v), max(u, v)))
        assert len(edges) == 24
        for (lo, hi) in edges:
            assert incidence_defect(cfg, Edge(lo, hi), expected=3) == 0

    def test_three_valent_vertices_have_count_four(self, rd_corona) -> None:
        _, configs = rd_corona
        cfg = configs[0]
        degree = _vertex_degree_map(cfg.central)
        three_valent_count = 0
        for i in range(len(cfg.central.vertices)):
            if degree[i] == 3:
                three_valent_count += 1
                assert incidence_defect(cfg, Vertex(i), expected=4) == 0
        assert three_valent_count == 8

    def test_four_valent_vertices_have_count_six(self, rd_corona) -> None:
        _, configs = rd_corona
        cfg = configs[0]
        degree = _vertex_degree_map(cfg.central)
        four_valent_count = 0
        for i in range(len(cfg.central.vertices)):
            if degree[i] == 4:
                four_valent_count += 1
                assert incidence_defect(cfg, Vertex(i), expected=6) == 0
        assert four_valent_count == 6
