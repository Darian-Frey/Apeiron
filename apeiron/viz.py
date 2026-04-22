"""3D visualisation of polyhedra, coronas, and symmetry actions.

This module is **outside** the exact-arithmetic verification pipeline
(CLAUDE.md §7.1, §7.3). Floats appear freely here — at the viz
boundary, ×2-stored ``Vec3`` values are divided by 2 and converted to
Python floats for plotly, and the rendered output is consumed only by
humans and HTML files. No function in this module is ever imported by
the core verification pipeline. The inverse rule — ``viz`` imports
core but never the reverse — is what keeps the exact arithmetic
uncontaminated.

Public API
----------

``render_polyhedron(P)`` — a single ``Polyhedron`` as a plotly
``Figure`` with translucent faces (fan-triangulated for polygonal
faces), edge lines, and vertex markers.

``render_polyhedra(Ps)`` — multiple polyhedra in one scene with
distinct colours. Useful for visualising coronas once ``corona.py``
produces them.

``render_rotations(P, rotations, spacing)`` — a grid of translated
copies of ``P``, each under one rotation. Sanity-check for I-symmetry
by eyeball.

Dependencies come from the ``[viz]`` optional extras in
``pyproject.toml``. Install with ``uv pip install -e '.[viz]'``.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import TYPE_CHECKING

from apeiron.polyhedron import Polyhedron
from apeiron.symmetry import Rotation, Vec3
from apeiron.zphi import ZPhi

if TYPE_CHECKING:
    import plotly.graph_objects as go

__all__ = [
    "render_polyhedra",
    "render_polyhedron",
    "render_rotations",
]


_PHI_FLOAT: float = (1.0 + math.sqrt(5.0)) / 2.0


def _zphi_to_real_float(z: ZPhi) -> float:
    """Convert a ×2-stored ZPhi component to its real-valued float.

    Storage is ``2·real``, so the float embedding is
    ``(a + b * phi) / 2``. The ``/ 2`` is the one point at which the
    ×2 convention surfaces on the viz side; in the core pipeline it
    stays implicit (CLAUDE.md §3.2).
    """
    return (z.a + z.b * _PHI_FLOAT) / 2.0


def _vec3_to_xyz(v: Vec3) -> tuple[float, float, float]:
    """Convert a ×2-stored Vec3 to a (float, float, float) triple."""
    return (
        _zphi_to_real_float(v.x),
        _zphi_to_real_float(v.y),
        _zphi_to_real_float(v.z),
    )


def _fan_triangulate(
    face: tuple[int, ...],
) -> list[tuple[int, int, int]]:
    """Fan-triangulate a polygonal face index cycle from its first vertex.

    Returns a list of triangle index triples. Each triangle is
    ``(face[0], face[k], face[k+1])`` for ``k`` in ``1 .. n-2``. Works
    for any convex face; non-convex faces aren't produced by the
    core pipeline (convex hull + coplanar merge gives only convex
    polygonal faces).
    """
    return [(face[0], face[k], face[k + 1]) for k in range(1, len(face) - 1)]


def _collect_edges(faces: Sequence[tuple[int, ...]]) -> set[tuple[int, int]]:
    """Undirected edge set as ``(min, max)`` index pairs."""
    edges: set[tuple[int, int]] = set()
    for face in faces:
        n = len(face)
        for idx in range(n):
            u, v = face[idx], face[(idx + 1) % n]
            edges.add((min(u, v), max(u, v)))
    return edges


def render_polyhedron(
    P: Polyhedron,
    *,
    color: str = "lightblue",
    opacity: float = 0.5,
    edge_color: str = "black",
    edge_width: float = 3.0,
    show_vertices: bool = True,
    vertex_color: str = "crimson",
    vertex_size: float = 4.0,
    title: str | None = None,
) -> "go.Figure":
    """Render a single ``Polyhedron`` as a plotly ``Figure``.

    The figure contains three traces: a translucent ``Mesh3d`` for the
    faces (fan-triangulated), a ``Scatter3d`` line trace for the
    edges, and an optional ``Scatter3d`` marker trace for the
    vertices.

    The aspect-ratio mode is ``"data"`` so the visualised polyhedron
    preserves its true proportions — crucial for eyeballing
    φ-related edge ratios on icosahedral-symmetric tiles like the
    rhombic triacontahedron.
    """
    import plotly.graph_objects as go

    xyz = [_vec3_to_xyz(v) for v in P.vertices]
    xs = [p[0] for p in xyz]
    ys = [p[1] for p in xyz]
    zs = [p[2] for p in xyz]

    # Face mesh (fan-triangulated).
    i_idx: list[int] = []
    j_idx: list[int] = []
    k_idx: list[int] = []
    for face in P.faces:
        for (a, b, c) in _fan_triangulate(face):
            i_idx.append(a)
            j_idx.append(b)
            k_idx.append(c)
    mesh = go.Mesh3d(
        x=xs,
        y=ys,
        z=zs,
        i=i_idx,
        j=j_idx,
        k=k_idx,
        color=color,
        opacity=opacity,
        flatshading=True,
        name="faces",
        showlegend=False,
    )

    # Edge lines. Scatter3d uses None-separated segments.
    edge_x: list[float | None] = []
    edge_y: list[float | None] = []
    edge_z: list[float | None] = []
    for (u, v) in _collect_edges(P.faces):
        edge_x.extend([xs[u], xs[v], None])
        edge_y.extend([ys[u], ys[v], None])
        edge_z.extend([zs[u], zs[v], None])
    edges = go.Scatter3d(
        x=edge_x,
        y=edge_y,
        z=edge_z,
        mode="lines",
        line={"color": edge_color, "width": edge_width},
        name="edges",
        showlegend=False,
        hoverinfo="skip",
    )

    traces: list = [mesh, edges]
    if show_vertices:
        verts = go.Scatter3d(
            x=xs,
            y=ys,
            z=zs,
            mode="markers",
            marker={"color": vertex_color, "size": vertex_size},
            name="vertices",
            showlegend=False,
            hovertemplate="v[%{text}]: (%{x:.3f}, %{y:.3f}, %{z:.3f})<extra></extra>",
            text=[str(i) for i in range(len(xyz))],
        )
        traces.append(verts)

    fig = go.Figure(data=traces)
    fig.update_layout(
        title=title,
        scene={
            "aspectmode": "data",
            "xaxis_title": "x",
            "yaxis_title": "y",
            "zaxis_title": "z",
        },
        margin={"l": 0, "r": 0, "t": 32 if title else 0, "b": 0},
    )
    return fig


def render_polyhedra(
    polyhedra: Sequence[Polyhedron],
    *,
    colors: Sequence[str] | None = None,
    opacity: float = 0.5,
    translations: Sequence[tuple[float, float, float]] | None = None,
    show_edges: bool = True,
    edge_color: str = "black",
    edge_width: float = 2.0,
    show_vertices: bool = False,
    title: str | None = None,
) -> "go.Figure":
    """Render several polyhedra in one scene.

    ``translations`` shifts each polyhedron in real (float) coordinates
    before rendering — useful for side-by-side comparison or for
    spacing out a set of rotated copies. If ``None``, polyhedra are
    rendered in place (which is what you want for a corona: every
    tile is already placed at its correct position in the overlay).

    ``colors`` is per-polyhedron; if ``None`` a small rotating palette
    is used.
    """
    import plotly.graph_objects as go

    palette = [
        "lightblue",
        "lightsalmon",
        "lightgreen",
        "plum",
        "khaki",
        "lightsteelblue",
        "wheat",
    ]
    if colors is None:
        colors = [palette[i % len(palette)] for i in range(len(polyhedra))]
    if translations is None:
        translations = [(0.0, 0.0, 0.0)] * len(polyhedra)
    if len(colors) != len(polyhedra) or len(translations) != len(polyhedra):
        raise ValueError(
            "colors and translations, when provided, must match the "
            "polyhedra count."
        )

    traces: list = []
    for P, color, (tx, ty, tz) in zip(polyhedra, colors, translations, strict=True):
        xyz = [_vec3_to_xyz(v) for v in P.vertices]
        xs = [p[0] + tx for p in xyz]
        ys = [p[1] + ty for p in xyz]
        zs = [p[2] + tz for p in xyz]
        i_idx: list[int] = []
        j_idx: list[int] = []
        k_idx: list[int] = []
        for face in P.faces:
            for (a, b, c) in _fan_triangulate(face):
                i_idx.append(a)
                j_idx.append(b)
                k_idx.append(c)
        traces.append(
            go.Mesh3d(
                x=xs, y=ys, z=zs,
                i=i_idx, j=j_idx, k=k_idx,
                color=color,
                opacity=opacity,
                flatshading=True,
                showlegend=False,
            )
        )
        if show_edges:
            edge_x: list[float | None] = []
            edge_y: list[float | None] = []
            edge_z: list[float | None] = []
            for (u, v) in _collect_edges(P.faces):
                edge_x.extend([xs[u], xs[v], None])
                edge_y.extend([ys[u], ys[v], None])
                edge_z.extend([zs[u], zs[v], None])
            traces.append(
                go.Scatter3d(
                    x=edge_x, y=edge_y, z=edge_z,
                    mode="lines",
                    line={"color": edge_color, "width": edge_width},
                    showlegend=False,
                    hoverinfo="skip",
                )
            )
        if show_vertices:
            traces.append(
                go.Scatter3d(
                    x=xs, y=ys, z=zs,
                    mode="markers",
                    marker={"color": edge_color, "size": 3.0},
                    showlegend=False,
                    hoverinfo="skip",
                )
            )

    fig = go.Figure(data=traces)
    fig.update_layout(
        title=title,
        scene={"aspectmode": "data"},
        margin={"l": 0, "r": 0, "t": 32 if title else 0, "b": 0},
    )
    return fig


def render_rotations(
    P: Polyhedron,
    rotations: Sequence[Rotation],
    *,
    spacing: float | None = None,
    per_row: int = 6,
    **polyhedra_kwargs: object,
) -> "go.Figure":
    """Lay out ``P`` rotated by each element of ``rotations`` in a grid.

    Grid layout: ``per_row`` per row in x, further rows marching in y,
    further slabs in z. ``spacing`` is the real-coordinate step
    between grid cells; defaults to a bounding-box-diameter-aware
    value derived from ``P`` itself so cells never overlap.

    Useful for a quick eyeball check that every element of the
    icosahedral group produces a distinct oriented copy of ``P`` —
    pass ``ICOSAHEDRAL`` from ``apeiron.symmetry`` as the rotation
    list.
    """
    # Compute a reasonable spacing from P's real-coord bounding-box
    # diameter if none given.
    if spacing is None:
        spacing = _real_bounding_box_diameter(P) * 1.2

    rotated: list[Polyhedron] = [P.apply(g) for g in rotations]
    translations: list[tuple[float, float, float]] = []
    for idx in range(len(rotated)):
        row = idx // per_row
        col = idx % per_row
        translations.append((col * spacing, -row * spacing, 0.0))

    return render_polyhedra(
        rotated,
        translations=translations,
        **polyhedra_kwargs,  # type: ignore[arg-type]
    )


def _real_bounding_box_diameter(P: Polyhedron) -> float:
    """Diameter of the real-coord bounding box of ``P``.

    Used by ``render_rotations`` to pick a default grid spacing; a
    float-only computation, never consumed by core code.
    """
    xyz = [_vec3_to_xyz(v) for v in P.vertices]
    xs = [p[0] for p in xyz]
    ys = [p[1] for p in xyz]
    zs = [p[2] for p in xyz]
    dx = max(xs) - min(xs)
    dy = max(ys) - min(ys)
    dz = max(zs) - min(zs)
    return max(dx, dy, dz)
