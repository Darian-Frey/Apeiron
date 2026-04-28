"""Render the four Danzer ABCK prototiles to interactive HTML.

Quick-demo script for the ``apeiron.viz`` module on Track A's
baseline candidates: loads each prototile from
``candidates/danzer/{A,B,C,K}.json``, renders them side-by-side as
a plotly scene, writes ``notebooks/danzer.html`` for browser
viewing.

Run from the repo root::

    .venv/bin/python notebooks/danzer_demo.py

The four tetrahedra are translated apart in the rendered scene so
they don't all visually pile up at the origin (each prototile, in
its canonical Frettlöh-Table-1 form, has one vertex at the
origin). Click-and-drag rotates; scroll zooms; hover on a vertex
to see its (x, y, z) coordinates.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

from apeiron.util import load_candidate
from apeiron.viz import render_polyhedra

REPO_ROOT = Path(__file__).resolve().parent.parent
DANZER_DIR = REPO_ROOT / "candidates" / "danzer"
OUTPUT_HTML = REPO_ROOT / "notebooks" / "danzer.html"


def main() -> int:
    letters = ("A", "B", "C", "K")
    tiles = []
    for letter in letters:
        path = DANZER_DIR / f"{letter}.json"
        if not path.exists():
            print(f"error: {path} not found", file=sys.stderr)
            return 1
        tiles.append(load_candidate(path))
    print(
        "Loaded Danzer ABCK: "
        + ", ".join(
            f"{letter} (V={len(t.vertices)} F={len(t.faces)})"
            for letter, t in zip(letters, tiles, strict=True)
        )
    )

    # Spread the four tetrahedra along x by the per-tile diameter.
    # Each Danzer tile fits in a roughly φ²-side bounding box, so a
    # spacing of 4.0 keeps them visually separate.
    spacing = 4.0
    translations = [(i * spacing, 0.0, 0.0) for i in range(len(tiles))]
    palette = ["lightblue", "lightsalmon", "lightgreen", "khaki"]

    fig = render_polyhedra(
        tiles,
        translations=translations,
        colors=palette,
        opacity=0.6,
        show_edges=True,
        edge_color="navy",
        edge_width=2.5,
        show_vertices=True,
        title=(
            "Danzer ABCK prototiles (Frettlöh Table 1) — "
            "A blue · B salmon · C green · K khaki"
        ),
    )
    fig.write_html(OUTPUT_HTML, include_plotlyjs="cdn")
    print(f"wrote {OUTPUT_HTML.relative_to(REPO_ROOT)}")
    _try_open(OUTPUT_HTML)
    return 0


def _try_open(path: Path) -> None:
    """Best-effort cross-platform open-in-browser; silent on failure."""
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
