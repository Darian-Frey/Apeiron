"""Shared helpers: candidate-file loaders, canonical forms, small utilities.

Keep this module small. If a helper grows past a few functions, promote it
to its own module rather than letting ``util`` bloat.
"""

from __future__ import annotations

import json
from pathlib import Path

from apeiron.polyhedron import Polyhedron
from apeiron.symmetry import Vec3
from apeiron.zphi import ZPhi

__all__ = ["load_candidate"]


_VALID_SCALE_DENOMS: tuple[int, ...] = (1, 2)


def load_candidate(
    path: Path | str,
    *,
    allow_interior_inputs: bool = False,
) -> Polyhedron:
    """Load a candidate tile JSON file and return its ``Polyhedron``.

    **Author-facing JSON schema.** Vertex coordinates are written as
    integer-ZPhi pairs in natural mathematical form; the loader handles
    the ×2 storage convention (CLAUDE.md §3.2, §3.3) internally.

    ::

        {
          "name": "rhombic_triacontahedron",
          "scale_denom": 1,
          "vertices": [
            [[a, b], [a, b], [a, b]],
            ...
          ]
        }

    - ``name`` (str, required): non-empty tile identifier.
    - ``scale_denom`` (int, optional, default 1): denominator applied to
      every real vertex coordinate. Must be 1 or 2. The ``2`` case
      arises when a tile's vertices are naturally in (½)ℤ[φ]³ — e.g. a
      substitution-rule image whose scaling is not integer-ZPhi. Any
      other value raises ``ValueError``.
    - ``vertices`` (list, required, length ≥ 4): each vertex is a
      3-element list; each element is a 2-int pair ``[a, b]``
      representing ``a + b*phi``.

    **Loader semantics.** Every real vertex ``v`` becomes ``2*v /
    scale_denom`` in ×2 storage form. Since ``scale_denom ∈ {1, 2}``
    and the input components are integers, the division is always
    exact. The resulting ``Vec3`` list feeds ``Polyhedron.from_vertices``
    for hull construction + exact ZPhi validation + coplanar-face-merge.

    Any additional top-level keys in the JSON are ignored silently, so
    richer metadata (provenance, intended symmetry class, research
    track, aperiodicity claim) can appear in candidate files and be
    consumed by future tooling without breaking this loader.

    Parameters
    ----------
    path : Path or str
        File path to the candidate JSON.
    allow_interior_inputs : bool, keyword-only, default False
        If ``False`` (default), raise ``ValueError`` when any input
        vertex ends up classified as interior to the convex hull of
        the others. This catches silent coordinate-spec errors — the
        RTH-era bug where four intended vertices turned out to be
        collinear with two others, dropping them as non-extreme and
        producing a 20-vertex polytope from a 32-vertex input. For
        vertex-set specs this is always a spec error, so default-on
        is correct. Set ``True`` for legitimate cases where the input
        deliberately includes interior points (e.g., an
        over-generated substitution-rule output with pending pruning).

    Raises
    ------
    FileNotFoundError
        If ``path`` does not exist.
    ValueError
        If the JSON schema is malformed, ``scale_denom`` is outside
        ``{1, 2}``, the vertex set does not form a valid convex
        polytope (propagated from ``Polyhedron.from_vertices``), or
        ``allow_interior_inputs`` is ``False`` and any input vertex
        is classified as interior.
    """
    p = Path(path)
    with p.open("r", encoding="utf-8") as f:
        raw = json.load(f)

    if not isinstance(raw, dict):
        raise ValueError(
            f"{p}: top-level JSON must be an object; got {type(raw).__name__}."
        )

    name = raw.get("name")
    if not isinstance(name, str) or not name:
        raise ValueError(
            f"{p}: missing or empty required field 'name' (non-empty str)."
        )

    if "vertices" not in raw:
        raise ValueError(f"{p}: missing required field 'vertices' (list).")
    raw_vertices = raw["vertices"]
    if not isinstance(raw_vertices, list):
        raise ValueError(
            f"{p}: 'vertices' must be a list; got {type(raw_vertices).__name__}."
        )

    scale_denom = raw.get("scale_denom", 1)
    if isinstance(scale_denom, bool) or not isinstance(scale_denom, int):
        raise ValueError(
            f"{p}: scale_denom must be an int; got {type(scale_denom).__name__}."
        )
    if scale_denom not in _VALID_SCALE_DENOMS:
        raise ValueError(
            f"{p}: scale_denom must be in {_VALID_SCALE_DENOMS}; got "
            f"{scale_denom}."
        )

    vertices = [
        _parse_vertex(raw_v, scale_denom, p, i)
        for i, raw_v in enumerate(raw_vertices)
    ]
    polyhedron = Polyhedron.from_vertices(vertices)

    if not allow_interior_inputs:
        hull_set = set(polyhedron.vertices)
        interior_idx = [i for i, v in enumerate(vertices) if v not in hull_set]
        if interior_idx:
            raise ValueError(
                f"{p}: {len(interior_idx)} of {len(vertices)} input vertices "
                f"were classified as interior to the convex hull of the others "
                f"(indices: {interior_idx}). For a vertex-set spec this is "
                "always a coordinate error — intended hull vertices should be "
                "extreme. Pass allow_interior_inputs=True if the interior "
                "points are deliberate (e.g. an over-generated substitution-"
                "rule output to be pruned)."
            )

    return polyhedron


def _parse_vertex(
    raw_v: object,
    scale_denom: int,
    path: Path,
    index: int,
) -> Vec3:
    """Parse one vertex: three ``[a, b]`` pairs → Vec3 in ×2 storage.

    The transformation is ``stored = 2 * author / scale_denom``, which
    for ``scale_denom ∈ {1, 2}`` and integer author coefficients is an
    exact integer operation. ``scale_denom = 2`` cancels the ×2
    multiplier and passes the author coefficients through unchanged
    (the case where the natural (½)ℤ[φ] form was written directly
    as integer ZPhi); ``scale_denom = 1`` doubles every coefficient.
    """
    if not isinstance(raw_v, list) or len(raw_v) != 3:
        raise ValueError(
            f"{path}: vertex[{index}] must be a 3-element list; got {raw_v!r}."
        )
    components: list[ZPhi] = []
    for j, raw_c in enumerate(raw_v):
        if (
            not isinstance(raw_c, list)
            or len(raw_c) != 2
            or not all(
                isinstance(x, int) and not isinstance(x, bool) for x in raw_c
            )
        ):
            raise ValueError(
                f"{path}: vertex[{index}] component {j} must be a "
                f"2-element list of ints; got {raw_c!r}."
            )
        a, b = raw_c
        stored_a = (2 * a) // scale_denom
        stored_b = (2 * b) // scale_denom
        components.append(ZPhi(stored_a, stored_b))
    return Vec3(components[0], components[1], components[2])
