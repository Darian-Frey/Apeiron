"""Deformation-search infrastructure for Track A (CLAUDE.md §6.1).

Per Claude (web) Q6 ruling 2026-04-29:

- **Q6a** Face-merge formalism: enumerate pairs of children in a
  σ-dissection that share a face; the merged tile is the polyhedron
  formed by removing the shared face. Vertex perturbation and
  parametric families are deferred secondary tools.
- **Q6b** Brute-force enumeration with cheap prioritisation filters
  applied before pillar-2 BFS (combinatorial PF feasibility, face-
  count sanity, vertex-class compatibility).
- **Q6c** Auto-generate ``candidates/<name>/fourth_pillar.py`` stub
  for every face-merge candidate at scaffold time, regardless of
  pillar status. The directory tree records every candidate
  *considered*, not just every candidate *surviving*.

This module is the entry point for Phase 1.5 of the roadmap. It does
not run any pillar checks itself — it produces candidates and helpers
that the existing ``hierarchy.aperiodicity_witness`` consumes.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from apeiron.polyhedron import Polyhedron
from apeiron.substitution import SubstitutionRule
from apeiron.symmetry import ImproperRotation, Rotation, Vec3

__all__ = [
    "FaceAdjacentPair",
    "MergeFiltering",
    "face_adjacent_pairs",
    "merge_two_tiles",
    "placed_vertices",
    "prioritise_merge_candidate",
    "scaffold_merge_candidate",
]


def placed_vertices(
    prototile: Polyhedron,
    pose: tuple[Vec3, Rotation | ImproperRotation],
) -> tuple[Vec3, ...]:
    """Apply ``(translation, rotation)`` to ``prototile.vertices``.

    For each vertex ``v`` in the prototile's local frame, returns
    ``rotation.apply(v) + translation``. Rotation is the
    Apeiron-storage form (×2-stored numerator); translation is the
    placed origin in ×2-stored ℤ[φ]³.
    """
    translation, rotation = pose
    return tuple(rotation.apply(v) + translation for v in prototile.vertices)


@dataclass(frozen=True, slots=True)
class FaceAdjacentPair:
    """A pair of face-adjacent children in some σ(parent).

    Per Claude (web) Q6a 2026-04-29: two children of σ(parent) share
    a face iff their placed vertex sets share *exactly three* points.
    For ABCK tetrahedral prototiles this is the only face-adjacency
    notion; for higher-face-count prototiles the framework would
    extend to "share exactly k points where k is the polygon size".

    ``i`` and ``j`` (with ``i < j``) index into
    ``rule.dissections[parent_index]``. ``shared_vertices`` is the
    frozenset of three placed Vec3s that lie on the shared face.
    """

    parent_index: int
    i: int
    j: int
    shared_vertices: frozenset[Vec3]

    def __post_init__(self) -> None:
        if self.i >= self.j:
            raise ValueError(
                f"FaceAdjacentPair requires i < j; got i={self.i}, j={self.j}."
            )
        if len(self.shared_vertices) != 3:
            raise ValueError(
                "Tetrahedral face-adjacency requires exactly 3 "
                f"shared vertices; got {len(self.shared_vertices)}."
            )


def face_adjacent_pairs(
    rule: SubstitutionRule,
    parent_index: int,
    prototile_polyhedra: Sequence[Polyhedron],
) -> tuple[FaceAdjacentPair, ...]:
    """Enumerate all pairs of children in ``σ(rule.prototile_{parent_index})``
    that share a triangular face.

    Per Claude (web) Q6a: brute-force pairwise comparison of placed
    vertex sets. For tetrahedra, two children share a face iff
    ``|placed(c_i) ∩ placed(c_j)| == 3``. The candidate set per
    parent is bounded by ``C(|σ(parent)|, 2)`` — for ABCK this is
    at most ``C(11, 2) = 55``, dropping to a small number after the
    "exactly 3 shared points" filter.

    Parameters
    ----------
    rule
        SubstitutionRule whose dissection of prototile_{parent_index}
        is being examined.
    parent_index
        Which of the rule's prototiles to enumerate.
    prototile_polyhedra
        Per-prototile-index Polyhedron, used to compute placed vertex
        sets via ``placed_vertices``. ``prototile_polyhedra[t]`` is
        the polyhedron for prototile-type ``t``.
    """
    if not 0 <= parent_index < rule.n_prototiles:
        raise ValueError(
            f"parent_index {parent_index} outside [0, {rule.n_prototiles})."
        )
    children = rule.dissections[parent_index]
    placed: list[frozenset[Vec3]] = []
    for child in children:
        proto = prototile_polyhedra[child.prototile_index]
        placed.append(frozenset(
            placed_vertices(proto, (child.translation, child.rotation))
        ))
    pairs: list[FaceAdjacentPair] = []
    for i in range(len(children)):
        for j in range(i + 1, len(children)):
            shared = placed[i] & placed[j]
            if len(shared) == 3:
                pairs.append(FaceAdjacentPair(
                    parent_index=parent_index,
                    i=i, j=j,
                    shared_vertices=shared,
                ))
    return tuple(pairs)


def merge_two_tiles(
    prototile_a: Polyhedron,
    pose_a: tuple[Vec3, Rotation | ImproperRotation],
    prototile_b: Polyhedron,
    pose_b: tuple[Vec3, Rotation | ImproperRotation],
    shared_vertices: frozenset[Vec3],
) -> Polyhedron:
    """Merge two face-adjacent tetrahedra along their shared face.

    Per Claude (web) Q6a: the merged tile is the polyhedron formed
    by removing the shared interior face. For tetrahedra glued on
    a triangular face, the result has ``4 + 4 − (3 shared) = 5``
    vertices and ``4 + 4 − 2 = 6`` triangular faces (the two faces
    incident to the shared triangle are both removed).

    Face winding: prototiles' faces are oriented outward in their
    local frame. A ``Rotation`` placement preserves orientation;
    an ``ImproperRotation`` (det = -1) reverses it, so faces from a
    tile placed via ``ImproperRotation`` have their cycles reversed
    here to keep outward-pointing normals on the merged polyhedron.

    Parameters
    ----------
    prototile_a, prototile_b
        The two prototiles being glued (their local-frame Polyhedron
        objects).
    pose_a, pose_b
        Each tile's (translation, rotation) placement in the parent's
        σ-frame.
    shared_vertices
        The three Vec3s lying on the shared face — must match exactly
        three vertices of each placed tile.

    Returns
    -------
    Polyhedron
        Canonical merged polyhedron, suitable for downstream pillar-1
        / pillar-2 / pillar-4 work.
    """
    if len(shared_vertices) != 3:
        raise ValueError(
            f"Tetrahedral merge requires 3 shared vertices; "
            f"got {len(shared_vertices)}."
        )

    placed_a = list(placed_vertices(prototile_a, pose_a))
    placed_b = list(placed_vertices(prototile_b, pose_b))

    placed_a_set = set(placed_a)
    placed_b_set = set(placed_b)
    if not shared_vertices <= placed_a_set:
        raise ValueError(
            "shared_vertices not a subset of placed prototile_a's vertices."
        )
    if not shared_vertices <= placed_b_set:
        raise ValueError(
            "shared_vertices not a subset of placed prototile_b's vertices."
        )

    # Winding: reverse face cycles for tiles placed via ImproperRotation
    # so the merged polyhedron's faces all point outward.
    def _maybe_reverse(face: tuple[int, ...], pose) -> tuple[int, ...]:
        _t, g = pose
        if isinstance(g, ImproperRotation):
            return tuple(reversed(face))
        return face

    # For each tile, drop the face whose vertex set equals shared_vertices.
    def _faces_minus_shared(
        prototile: Polyhedron, placed: list[Vec3], pose,
    ) -> list[list[Vec3]]:
        out: list[list[Vec3]] = []
        for face_idx in prototile.faces:
            face_v = {placed[i] for i in face_idx}
            if face_v == shared_vertices:
                continue
            face_idx = _maybe_reverse(face_idx, pose)
            out.append([placed[i] for i in face_idx])
        return out

    merged_face_v_lists = (
        _faces_minus_shared(prototile_a, placed_a, pose_a)
        + _faces_minus_shared(prototile_b, placed_b, pose_b)
    )

    merged_vertex_set = placed_a_set | placed_b_set
    merged_vertex_list = list(merged_vertex_set)
    vertex_to_idx = {v: i for i, v in enumerate(merged_vertex_list)}
    face_index_lists = [
        [vertex_to_idx[v] for v in face_verts]
        for face_verts in merged_face_v_lists
    ]
    return Polyhedron.from_raw(merged_vertex_list, face_index_lists)


# -- Q6b prioritisation filter ---------------------------------------


@dataclass(frozen=True, slots=True)
class MergeFiltering:
    """Outcome of running Q6b prioritisation filters on a merge candidate.

    Per Claude (web) Q6b 2026-04-29: cheap checks applied *before* the
    expensive pillar-2 BFS. A candidate that fails any filter is
    eliminated from the deformation search; a candidate that passes
    all filters proceeds to pillar-1+2+3 verification via
    ``hierarchy.aperiodicity_witness``.

    Three filters:

    - ``combinatorial_feasibility``: is there an integer ``k`` such
      that the merged tile's σ-decomposition equals ``k`` copies of
      itself? For a 2-tile face-merge this requires the parent's σ
      count vector to be ``k·(component multiset)`` exactly. For
      ABCK 2-tile merges this filter ALWAYS fails (no σ-count vector
      matches any 2-tile composition); the failure constitutes
      no-go evidence per the roadmap's Phase 1.5.
    - ``face_count_strictly_decreases``: the merged polyhedron has
      strictly fewer faces than the sum of its components' face
      counts. For tetrahedral face-merge this is automatic
      (``F_merged = F_a + F_b - 2`` since the shared face is removed
      from both). Always True for valid tetrahedral merges; included
      for symmetry with future polyhedral candidates.
    - ``vertex_class_compatible``: Frettlöh's vertex classification
      (I/II/III/IV) is consistent across the merged tile.
      *Not implemented yet* — requires per-prototile vertex-class
      data in the candidate JSON, which ABCK encodes only as
      docstring text. Filter returns ``None`` (= "unchecked") for
      now; future deformation-search candidates may add structured
      vertex-class data and enable the check.

    ``passes_all`` is True iff every implemented filter (1 and 2)
    returned True, AND filter 3 (if implemented) returned True. Used
    as the gating condition for proceeding to pillar-1+2+3 work.
    """

    combinatorial_feasibility: bool
    combinatorial_k: int | None
    face_count_strictly_decreases: bool
    vertex_class_compatible: bool | None
    detail: str

    @property
    def passes_all(self) -> bool:
        if not self.combinatorial_feasibility:
            return False
        if not self.face_count_strictly_decreases:
            return False
        if self.vertex_class_compatible is False:
            return False
        return True


def prioritise_merge_candidate(
    rule: SubstitutionRule,
    pair: FaceAdjacentPair,
    prototile_polyhedra: Sequence[Polyhedron],
) -> MergeFiltering:
    """Apply Q6b's three prioritisation filters to a face-merge
    candidate. Cheap O(1) checks; do not call ``aperiodicity_witness``.

    Returns a ``MergeFiltering`` summarising each filter's outcome
    plus a human-readable ``detail`` string. The deformation search
    consumes the result via ``MergeFiltering.passes_all``.
    """
    children = rule.dissections[pair.parent_index]
    type_i = children[pair.i].prototile_index
    type_j = children[pair.j].prototile_index
    proto_a = prototile_polyhedra[type_i]
    proto_b = prototile_polyhedra[type_j]

    # --- Filter 1: combinatorial feasibility ---
    # Merged tile P_AB has composition {type_i: 1, type_j: 1} (or
    # {type_i: 2} if type_i == type_j, though for face-merges of
    # distinct tetrahedra we expect i ≠ j).
    composition: dict[int, int] = {}
    composition[type_i] = composition.get(type_i, 0) + 1
    composition[type_j] = composition.get(type_j, 0) + 1

    # σ(P_AB) count vector = σ(type_i) + σ(type_j).
    sigma_count: list[int] = [0] * rule.n_prototiles
    for t in (type_i, type_j):
        for child in rule.dissections[t]:
            sigma_count[child.prototile_index] += 1

    # Look for integer k such that sigma_count == k · composition_vec.
    composition_vec = [composition.get(t, 0) for t in range(rule.n_prototiles)]
    feasible = False
    k_found: int | None = None
    # k is determined by any non-zero entry in composition; check
    # consistency across all entries.
    nonzero_idx = next(
        (idx for idx, c in enumerate(composition_vec) if c > 0), None,
    )
    if nonzero_idx is not None and composition_vec[nonzero_idx] > 0:
        if sigma_count[nonzero_idx] % composition_vec[nonzero_idx] == 0:
            k_candidate = (
                sigma_count[nonzero_idx] // composition_vec[nonzero_idx]
            )
            if k_candidate >= 1 and all(
                sigma_count[t] == k_candidate * composition_vec[t]
                for t in range(rule.n_prototiles)
            ):
                feasible = True
                k_found = k_candidate

    # --- Filter 2: face-count strictly decreases ---
    # Tetrahedral merge: F_merged = F_a + F_b - 2 < F_a + F_b.
    f_total = len(proto_a.faces) + len(proto_b.faces)
    f_merged = f_total - 2  # shared face removed from both
    face_count_ok = f_merged < f_total

    # --- Filter 3: vertex-class compatibility ---
    # Not implemented; requires structured vertex-class data per
    # prototile (Frettlöh I/II/III/IV) which ABCK encodes only as
    # docstring text in the JSON. Future candidates may add a
    # vertex_classes field; for now this filter returns None.
    vertex_class_ok: bool | None = None

    detail_parts = [
        f"composition={composition_vec}",
        f"σ_count={sigma_count}",
    ]
    if feasible:
        detail_parts.append(f"k={k_found}")
    else:
        detail_parts.append("not k·composition for any integer k")
    detail_parts.append(f"f_merged={f_merged} < f_total={f_total}")

    return MergeFiltering(
        combinatorial_feasibility=feasible,
        combinatorial_k=k_found,
        face_count_strictly_decreases=face_count_ok,
        vertex_class_compatible=vertex_class_ok,
        detail="; ".join(detail_parts),
    )


# -- Q6c scaffold ----------------------------------------------------


_FOURTH_PILLAR_TEMPLATE = '''"""Fourth-pillar stub for the {display_name} merge candidate.

**Status: scaffold-time stub.** Generated by
``apeiron.deformation.scaffold_merge_candidate`` per Claude (web)'s
Q6c 2026-04-29 ruling. Two reasons for an auto-generated stub:

1. *Exercise the ``FourthPillarArgument`` protocol's import path*
   for every candidate — a missing module would surface only at
   pillar-4 time, potentially weeks after scaffold.
2. *Capture the geometric context* while it's fresh. The candidate
   was discovered as a face-merge of {parent_summary}; the proof
   approach for pillar 4 is recorded below.

**Pillar-4 approach (proposed, unverified).** {pillar_4_approach}

If this candidate survives pillars 1, 2, 3 — verified by
``apeiron.hierarchy.aperiodicity_witness`` — pillar 4 graduates from
this stub to a concrete implementation. Until then, both methods
raise ``NotImplementedError`` with the candidate's discovery
context.
"""

from __future__ import annotations

from apeiron.corona import CoronaConfig
from apeiron.hierarchy import (
    FourthPillarArgument,
    HierarchicalCounterexample,
    HierarchicalWitness,
)


class {class_name}(FourthPillarArgument):
    """Stub ``FourthPillarArgument`` for the {display_name} candidate."""

    _CONTEXT = (
        "Face-merge candidate from {parent_summary}; "
        "shared face: {shared_face_description}."
    )

    def local_configurations(self) -> frozenset[CoronaConfig]:
        raise NotImplementedError(
            "Pillar-4 local-configuration enumeration not yet "
            f"encoded for {{self._CONTEXT}}"
        )

    def verify_hierarchical(
        self, _config: CoronaConfig,
    ) -> HierarchicalWitness | HierarchicalCounterexample:
        raise NotImplementedError(
            "Pillar-4 hierarchical verification not yet "
            f"encoded for {{self._CONTEXT}}"
        )
'''


def _to_pascal_case(name: str) -> str:
    """Convert a snake_case or hyphen-case identifier to PascalCase
    for class naming. ``"danzer_merge_ab"`` -> ``"DanzerMergeAb"``.
    """
    parts = name.replace("-", "_").split("_")
    return "".join(p.capitalize() for p in parts if p)


def scaffold_merge_candidate(
    name: str,
    parent_letters: tuple[str, ...],
    shared_face_description: str,
    *,
    pillar_4_approach: str = (
        "Pillar 4 would require showing every face-to-face tiling by "
        "the merged tile arises from the ABCK hierarchy. Likely "
        "approach: extend Frettlöh Theorem 1.2's blue-edge mirror "
        "rule to account for the merged geometry."
    ),
    candidates_dir: Path = Path("candidates"),
) -> Path:
    """Create ``candidates/<name>/`` with an auto-generated
    ``fourth_pillar.py`` stub.

    Parameters
    ----------
    name
        Directory-safe name for the candidate, e.g.
        ``"danzer_merge_ab_face_3"``. Used both as the directory
        name and (PascalCase'd) as the ``FourthPillarArgument``
        subclass name.
    parent_letters
        The ABCK letters of the two prototiles being merged, e.g.
        ``("A", "B")``. Recorded in the stub's docstring for
        traceability.
    shared_face_description
        Human-readable description of the shared face, e.g.
        ``"face 2 of A coincides with face 0 of B at translation
        ½(τ², τ², τ²)"``. Appears in the stub's ``_CONTEXT`` string.
    pillar_4_approach
        Proposed approach to pillar 4 for this candidate. Defaults
        to a Frettlöh-style blue-edge mirror generalisation; pass
        a candidate-specific approach when one is known.
    candidates_dir
        Where to create the candidate directory. Defaults to
        ``candidates/`` at the project root.

    Returns
    -------
    Path
        The created candidate directory. The caller can then write
        the merged-tile JSON, substitution rule, etc., into it.

    Notes
    -----
    Per Q6c, called once per face-merge candidate *discovered* —
    not per candidate that survives pillar checks. This is the right
    research record to keep: every candidate considered has a
    directory, every directory has a fourth-pillar stub, every stub
    documents the discovery context. Survivors graduate to manual
    pillar-4 work.

    Raises
    ------
    FileExistsError
        If ``candidates_dir / name`` already exists. To rebuild a
        candidate, the caller deletes the directory first — a
        deliberate choice to prevent accidental overwriting of
        partial pillar-4 work.
    """
    if not name or "/" in name or name.startswith("."):
        raise ValueError(
            f"name must be a directory-safe identifier; got {name!r}"
        )
    candidate_dir = candidates_dir / name
    if candidate_dir.exists():
        raise FileExistsError(
            f"{candidate_dir} already exists. Delete it first to rebuild."
        )
    candidate_dir.mkdir(parents=True, exist_ok=False)

    parent_summary = "+".join(parent_letters)
    display_name = name
    class_name = _to_pascal_case(name) + "FourthPillar"

    fourth_pillar_text = _FOURTH_PILLAR_TEMPLATE.format(
        display_name=display_name,
        parent_summary=parent_summary,
        pillar_4_approach=pillar_4_approach,
        shared_face_description=shared_face_description,
        class_name=class_name,
    )
    (candidate_dir / "fourth_pillar.py").write_text(fourth_pillar_text)
    (candidate_dir / "__init__.py").write_text(
        f'"""Merge candidate {display_name}: {parent_summary}.\n\n'
        f"Auto-scaffolded by ``apeiron.deformation.scaffold_merge_candidate``.\n"
        f'"""\n'
    )
    return candidate_dir
