"""Smoke tests for apeiron.viz.

The viz module's output is visual and the real QA is eyeballing the
generated HTML. These tests guard against *structural* regressions:
the render functions produce plotly ``Figure`` objects, the
``Mesh3d`` / ``Scatter3d`` traces contain the expected number of
vertices and face triangles, and nothing raises on the standard
fixtures (cube, RTH, I-orbit of the RTH).

``plotly`` is in the ``[viz]`` optional extras. These tests
``importorskip`` so a plain ``dev`` environment without plotly can
still run the rest of the suite.
"""

from __future__ import annotations

import pytest

go = pytest.importorskip("plotly.graph_objects")

from apeiron.polyhedron import Polyhedron   # noqa: E402
from apeiron.symmetry import ICOSAHEDRAL, ROT_2, Vec3   # noqa: E402
from apeiron.viz import render_polyhedra, render_polyhedron, render_rotations   # noqa: E402
from apeiron.zphi import ZPhi   # noqa: E402


def _v(x: tuple[int, int], y: tuple[int, int], z: tuple[int, int]) -> Vec3:
    return Vec3(ZPhi(*x), ZPhi(*y), ZPhi(*z))


def _unit_cube() -> Polyhedron:
    corners = [
        _v((2 * i, 0), (2 * j, 0), (2 * k, 0))
        for i in (0, 1) for j in (0, 1) for k in (0, 1)
    ]
    return Polyhedron.from_vertices(corners)


# -- render_polyhedron ------------------------------------------------


class TestRenderPolyhedron:
    def test_returns_plotly_figure(self) -> None:
        P = _unit_cube()
        fig = render_polyhedron(P)
        assert isinstance(fig, go.Figure)

    def test_mesh_vertex_count_matches(self) -> None:
        P = _unit_cube()
        fig = render_polyhedron(P)
        mesh = next(t for t in fig.data if isinstance(t, go.Mesh3d))
        assert len(mesh.x) == len(P.vertices) == 8

    def test_mesh_triangle_count_matches_fan_triangulation(self) -> None:
        # Cube: 6 quadrilateral faces → fan-triangulated to 2 triangles
        # each → 12 triangles total.
        P = _unit_cube()
        fig = render_polyhedron(P)
        mesh = next(t for t in fig.data if isinstance(t, go.Mesh3d))
        assert len(mesh.i) == 12
        assert len(mesh.j) == 12
        assert len(mesh.k) == 12

    def test_edge_trace_present(self) -> None:
        P = _unit_cube()
        fig = render_polyhedron(P)
        scatter_traces = [t for t in fig.data if isinstance(t, go.Scatter3d)]
        # At least one scatter trace for edges (and optionally one for
        # vertices when show_vertices=True, the default).
        assert len(scatter_traces) >= 1

    def test_vertex_trace_omitted_when_disabled(self) -> None:
        P = _unit_cube()
        fig = render_polyhedron(P, show_vertices=False)
        # Only the edges scatter trace; no vertex markers.
        scatter_traces = [t for t in fig.data if isinstance(t, go.Scatter3d)]
        assert len(scatter_traces) == 1


# -- render_polyhedra -------------------------------------------------


class TestRenderPolyhedra:
    def test_multiple_polyhedra(self) -> None:
        P = _unit_cube()
        fig = render_polyhedra([P, P])
        meshes = [t for t in fig.data if isinstance(t, go.Mesh3d)]
        assert len(meshes) == 2

    def test_per_polyhedron_translations(self) -> None:
        P = _unit_cube()
        fig = render_polyhedra(
            [P, P],
            translations=[(0.0, 0.0, 0.0), (5.0, 0.0, 0.0)],
        )
        meshes = [t for t in fig.data if isinstance(t, go.Mesh3d)]
        # Second mesh is translated +5 in x.
        assert max(meshes[1].x) - max(meshes[0].x) == pytest.approx(5.0)

    def test_rejects_mismatched_color_count(self) -> None:
        P = _unit_cube()
        with pytest.raises(ValueError, match="polyhedra count"):
            render_polyhedra([P, P], colors=["red"])

    def test_rejects_mismatched_translation_count(self) -> None:
        P = _unit_cube()
        with pytest.raises(ValueError, match="polyhedra count"):
            render_polyhedra([P, P], translations=[(0.0, 0.0, 0.0)])


# -- render_rotations -------------------------------------------------


class TestRenderRotations:
    def test_all_sixty_icosahedral_rotations(self) -> None:
        P = _unit_cube()
        fig = render_rotations(P, rotations=list(ICOSAHEDRAL))
        meshes = [t for t in fig.data if isinstance(t, go.Mesh3d)]
        assert len(meshes) == 60

    def test_small_subset(self) -> None:
        P = _unit_cube()
        fig = render_rotations(P, rotations=[ROT_2])
        meshes = [t for t in fig.data if isinstance(t, go.Mesh3d)]
        assert len(meshes) == 1
