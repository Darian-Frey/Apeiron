"""Render the rhombic triacontahedron to an interactive HTML file.

Quick-demo script for the ``apeiron.viz`` module: loads the RTH
candidate, renders it as a plotly ``Figure``, writes the output to
``notebooks/rth.html`` for browser viewing.

Run from the repo root::

    .venv/bin/python notebooks/rth_demo.py

Opens the HTML in the default browser (if ``xdg-open`` / ``open`` is
available) or just prints the path. Click-and-drag rotates; scroll
zooms; vertices show their index + (x, y, z) on hover.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

from apeiron.symmetry import ICOSAHEDRAL
from apeiron.util import load_candidate
from apeiron.viz import render_polyhedron, render_rotations


REPO_ROOT = Path(__file__).resolve().parent.parent
CANDIDATE_PATH = REPO_ROOT / "candidates" / "rhombic_triacontahedron.json"
SINGLE_HTML = REPO_ROOT / "notebooks" / "rth.html"
ORBIT_HTML = REPO_ROOT / "notebooks" / "rth_icosahedral_orbit.html"


def main() -> int:
    if not CANDIDATE_PATH.exists():
        print(f"error: candidate not found at {CANDIDATE_PATH}", file=sys.stderr)
        return 1

    rth = load_candidate(CANDIDATE_PATH)
    print(
        f"Loaded RTH: {len(rth.vertices)} vertices, {len(rth.faces)} faces."
    )

    # Single RTH.
    fig = render_polyhedron(
        rth,
        color="lightblue",
        opacity=0.55,
        edge_color="navy",
        edge_width=3.0,
        vertex_color="crimson",
        vertex_size=4.0,
        title="Rhombic triacontahedron (32 vertices, 30 rhombic faces)",
    )
    fig.write_html(SINGLE_HTML, include_plotlyjs="cdn")
    print(f"wrote {SINGLE_HTML.relative_to(REPO_ROOT)}")

    # All 60 icosahedral rotations, laid out in a grid. Visual
    # sanity-check that every g in I produces a genuinely rotated
    # copy of the RTH (no two cells in the grid look identical).
    orbit_fig = render_rotations(
        rth,
        rotations=list(ICOSAHEDRAL),
        per_row=10,
        opacity=0.6,
        show_vertices=False,
        title="Rhombic triacontahedron under all 60 elements of I",
    )
    orbit_fig.write_html(ORBIT_HTML, include_plotlyjs="cdn")
    print(f"wrote {ORBIT_HTML.relative_to(REPO_ROOT)}")

    # Try to open the single-RTH HTML in the browser.
    _try_open(SINGLE_HTML)
    return 0


def _try_open(path: Path) -> None:
    """Best-effort cross-platform open-in-browser. Silent on failure —
    the file is written either way and the caller can open manually.
    """
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
