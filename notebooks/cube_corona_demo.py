"""Render the cube's complete first corona to interactive HTML.

Exercises the ``corona`` ↔ ``viz`` integration end-to-end: builds a
unit cube, runs ``corona_1`` to find its unique 26-neighbour Moore
configuration, then renders the central + every neighbour as a
single plotly scene.

Run from the repo root::

    .venv/bin/python notebooks/cube_corona_demo.py

The corona_1 computation takes ~3 s; the render is ~1 s. Output:
``notebooks/cube_corona.html``. The 27 unit cubes (central + 26)
form a 3×3×3 block; click-and-drag rotates, scroll zooms.
"""

from __future__ import annotations

import subprocess
import sys
import time
from pathlib import Path

from apeiron.corona import corona_1
from apeiron.polyhedron import Polyhedron
from apeiron.symmetry import Vec3
from apeiron.viz import _vec3_to_xyz, render_polyhedra
from apeiron.zphi import ZPhi

REPO_ROOT = Path(__file__).resolve().parent.parent
OUTPUT_HTML = REPO_ROOT / "notebooks" / "cube_corona.html"


def _v(*components: int) -> Vec3:
    """Build a Vec3 from six ints — pairs (a, b) per coordinate."""
    return Vec3(
        ZPhi(components[0], components[1]),
        ZPhi(components[2], components[3]),
        ZPhi(components[4], components[5]),
    )


def main() -> int:
    cube_verts = [
        _v(2 * i, 0, 2 * j, 0, 2 * k, 0)
        for i in (0, 1)
        for j in (0, 1)
        for k in (0, 1)
    ]
    cube = Polyhedron.from_vertices(cube_verts)

    print("computing corona_1(cube)...")
    t0 = time.time()
    coronas = corona_1(cube, expected_edge_count=4, expected_vertex_count=8)
    print(f"  done in {time.time() - t0:.2f}s; found {len(coronas)} canonical")
    assert len(coronas) == 1
    corona = coronas[0]
    print(f"  {len(corona.neighbours)} neighbours")

    # Build the per-tile polyhedron + float translation.
    polyhedra: list[Polyhedron] = [corona.central]
    translations: list[tuple[float, float, float]] = [(0.0, 0.0, 0.0)]
    for placement in corona.neighbours:
        polyhedra.append(corona.central.apply(placement.rotation))
        translations.append(_vec3_to_xyz(placement.translation))

    # Distinct colours so the central is visually identifiable.
    colors = ["crimson"] + ["lightblue"] * len(corona.neighbours)

    fig = render_polyhedra(
        polyhedra,
        translations=translations,
        colors=colors,
        opacity=0.4,
        show_edges=True,
        edge_color="navy",
        edge_width=2.0,
        show_vertices=False,
        title=(
            "Cube first corona — 1 central (crimson) + 26 neighbours "
            "(lightblue): 6 face + 12 edge + 8 vertex sharers, "
            "forming a 3×3×3 block"
        ),
    )
    fig.write_html(OUTPUT_HTML, include_plotlyjs="cdn")
    print(f"wrote {OUTPUT_HTML.relative_to(REPO_ROOT)}")
    _try_open(OUTPUT_HTML)
    return 0


def _try_open(path: Path) -> None:
    for opener in ("xdg-open", "open"):
        try:
            subprocess.run(
                [opener, str(path)],
                check=False,
                capture_output=True,
                timeout=2,
            )
            return
        except (FileNotFoundError, subprocess.TimeoutExpired):
            continue


if __name__ == "__main__":
    sys.exit(main())
