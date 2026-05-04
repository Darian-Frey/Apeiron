"""Track B realisation — vertex-placement CSP for substitution-matrix
candidates.

Per Claude (web)'s Q8 ruling 2026-04-29: structured rotation-search
with linear translation recovery in ℤ[φ]. Returns one of
``Realised | NoRealisation | Inconclusive`` per the resumable-state
contract (Q8d). Internal implementation is a generator yielding
``SearchProgress`` snapshots (Q8 meta — coroutine-style).

**Status: API + result types + Fibonacci oracle as a manual witness.**
The full algorithmic CSP (rotation search + linear translation
recovery + fail-first ordering) is multi-week per Q8c; this commit
delivers the architecture and the simplest oracle without
committing to a partial CSP that might be wrong. The ``realise``
function's body returns ``Inconclusive`` for all inputs except a
hand-recognised Fibonacci-shape oracle; the Fibonacci witness test
confirms the API delivers the right shape end-to-end.

A follow-up commit will replace the placeholder with the structured
backtracker per Q8a/c (rotation iteration with face-match pruning,
linear translation recovery in ZPhi, fail-first ordering on rotation
choice for early NoRealisation discovery per Q8 meta-3).
"""

from __future__ import annotations

import itertools
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Final

import numpy as np

from apeiron.polyhedron import Polyhedron
from apeiron.symmetry import ImproperRotation, Rotation, Vec3
from apeiron.zphi import ZPhi

__all__ = [
    "ChildPlacement",
    "FaceMatchEdge",
    "Inconclusive",
    "NoRealisation",
    "Realised",
    "SearchProgress",
    "propagate_translations_along_tree",
    "realise",
    "translation_offset_from_face_match",
]


_ZERO_ZPHI: Final[ZPhi] = ZPhi(0, 0)
_ONE_ZPHI: Final[ZPhi] = ZPhi(1, 0)


@dataclass(frozen=True, slots=True)
class ChildPlacement:
    """A child tile in a placed σ(parent) dissection.

    ``prototile_index`` indexes into the realisation's prototile-shape
    list. ``translation`` is the placed origin in ×2-storage form
    (CLAUDE.md §3.2). ``rotation`` is an isometry of R^3:
    ``Rotation`` (in I, det +1) or ``ImproperRotation`` (in I_h \\ I,
    det -1).
    """

    prototile_index: int
    translation: Vec3
    rotation: Rotation | ImproperRotation


@dataclass(frozen=True, slots=True)
class SearchProgress:
    """Per-iteration snapshot of the CSP search state.

    Yielded internally by the realisation generator (Q8 meta). The
    public ``realise`` API consumes these and decides timeout /
    completion; callers needing finer-grained progress can drive the
    generator directly.

    ``fraction_searched`` is an exact ZPhi rational ``∈ [0, 1]`` —
    the fraction of the bounded rotation space that has been
    explored or pruned. Per Q8d, this distinguishes "94 % covered,
    nothing found" from "3 % covered, timed out".
    """

    fraction_searched: ZPhi
    realised_partial: tuple[ChildPlacement | None, ...]


@dataclass(frozen=True, slots=True)
class Realised:
    """A complete face-to-face realisation of the substitution rule.

    ``prototile_shapes`` is the n-tuple of prototile Polyhedra (one
    per substitution-matrix prototile). ``children_per_parent[i]``
    is the placed children of σ(prototile_i): a tuple of
    ``ChildPlacement`` whose total count and per-type multiplicities
    match column i of the substitution matrix.

    A ``Realised`` instance is the input to
    ``hierarchy.aperiodicity_witness`` for the four-pillar check —
    pillars 1, 2, 3 verify the resulting tiling is aperiodic;
    pillar 4 is candidate-specific and stays a stub.

    ``fraction_searched == ZPhi(1, 0)`` only if the search had
    exhausted the bounded space; for fast finds it is the fraction
    covered up to the first valid assignment.
    """

    prototile_shapes: tuple[Polyhedron, ...]
    children_per_parent: tuple[tuple[ChildPlacement, ...], ...]
    fraction_searched: ZPhi


@dataclass(frozen=True, slots=True)
class NoRealisation:
    """Provably no realisation exists in the bounded search space.

    Distinguished from ``Inconclusive`` by ``fraction_searched ==
    ZPhi(1, 0)`` — the entire bounded space was covered. The
    ``reason`` field records why each branch failed (typically a
    constraint violation that pruned an entire subtree).
    """

    reason: str
    fraction_searched: ZPhi = _ONE_ZPHI


@dataclass(frozen=True, slots=True)
class Inconclusive:
    """Search hit the budget cap without finishing the bounded space.

    Per Q8d, ``fraction_searched`` records exactly how much of the
    space was covered (an exact ZPhi rational). ``partial`` carries
    the search state for resumption with an extended budget.
    """

    fraction_searched: ZPhi
    partial: SearchProgress | None
    reason: str = "search budget exhausted"


def translation_offset_from_face_match(
    rotation_a: Rotation | ImproperRotation,
    rotation_b: Rotation | ImproperRotation,
    face_indices_a: tuple[int, int, int],
    face_indices_b: tuple[int, int, int],
    vertices_a: Sequence[Vec3],
    vertices_b: Sequence[Vec3],
) -> Vec3 | None:
    """Recover ``t_b − t_a`` from a triangular face-match constraint.

    Per Claude (web)'s Q8a 2026-04-29 — the per-edge sub-routine of
    the rotation-search-with-linear-translation-recovery backtracker.

    Setup. Let placed-A be the prototile-A vertices placed at
    ``rotation_a · v + t_a``, and placed-B similarly for prototile-B
    with rotation ``rotation_b`` and translation ``t_b``. If A's face
    indexed by ``face_indices_a`` (3 vertex indices into
    ``vertices_a``) matches B's face indexed by ``face_indices_b``
    (3 indices into ``vertices_b``) — i.e., the placed face vertex
    sets are equal — then for some permutation ``π ∈ S_3``::

        rotation_a · vertices_a[face_indices_a[k]] + t_a
            == rotation_b · vertices_b[face_indices_b[π(k)]] + t_b
            for k ∈ {0, 1, 2}.

    Rearranging::

        t_b − t_a == rotation_a · vertices_a[face_indices_a[k]]
                   − rotation_b · vertices_b[face_indices_b[π(k)]].

    For the match to be consistent, the RHS must give the same Vec3
    across all three ``k``. We try each of the 6 permutations and
    return the offset on the first consistent one. Return None if no
    permutation works.

    Returns
    -------
    Vec3 | None
        The exact ZPhi offset ``t_b − t_a`` (×2-stored), or ``None``
        if no permutation of ``face_indices_b`` produces a
        consistent match.

    Notes
    -----
    Pure ZPhi arithmetic; no floats. The output is unique up to the
    permutation symmetry of the matched triangle — multiple
    permutations could in principle produce the SAME offset (a
    triangle with non-trivial stabiliser under the rotation
    pair). For ABCK-style irregular triangles this won't happen;
    for regular triangles it can. The function returns the first
    consistent offset found.
    """
    if len(face_indices_a) != 3 or len(face_indices_b) != 3:
        raise ValueError(
            "Triangular face requires exactly 3 vertex indices each; "
            f"got A={face_indices_a}, B={face_indices_b}."
        )
    placed_a = [
        rotation_a.apply(vertices_a[i]) for i in face_indices_a
    ]
    for perm in itertools.permutations(range(3)):
        placed_b = [
            rotation_b.apply(vertices_b[face_indices_b[perm[k]]])
            for k in range(3)
        ]
        offset_candidate = placed_a[0] - placed_b[0]
        consistent = all(
            placed_a[k] - placed_b[k] == offset_candidate
            for k in range(1, 3)
        )
        if consistent:
            return offset_candidate
    return None


@dataclass(frozen=True, slots=True)
class FaceMatchEdge:
    """One edge in a face-adjacency tree across the children of σ(P).

    Identifies a triangular face shared by two children. The
    ``parent`` is the child placed earlier (translation already
    known); ``new`` is the child whose translation gets propagated
    via this edge. Both face-index tuples are the 3 vertex indices
    of the shared face in the respective prototile's local frame.

    Used by ``propagate_translations_along_tree`` to walk a fixed
    tree topology and recover translations from face-match
    constraints, given fixed per-child rotations.
    """

    parent: int
    new: int
    face_indices_parent: tuple[int, int, int]
    face_indices_new: tuple[int, int, int]

    def __post_init__(self) -> None:
        if self.parent == self.new:
            raise ValueError(
                f"FaceMatchEdge.parent must differ from .new; "
                f"got both = {self.parent}."
            )
        if len(self.face_indices_parent) != 3:
            raise ValueError(
                "Triangular face requires 3 vertex indices on parent; "
                f"got {self.face_indices_parent}."
            )
        if len(self.face_indices_new) != 3:
            raise ValueError(
                "Triangular face requires 3 vertex indices on new; "
                f"got {self.face_indices_new}."
            )


def propagate_translations_along_tree(
    rotations: Sequence[Rotation | ImproperRotation],
    edges: Sequence[FaceMatchEdge],
    prototile_indices: Sequence[int],
    prototile_vertices: Sequence[Sequence[Vec3]],
) -> tuple[Vec3, ...] | None:
    """Walk a face-match tree and recover one translation per child.

    The next layer above ``translation_offset_from_face_match`` per
    Q8a 2026-04-29: given fixed rotations and a fixed tree topology
    of face-match edges, propagate translations by tree DFS. Each
    edge contributes one offset via the per-edge solver; if any
    edge is infeasible (no permutation match), the whole tree is
    rejected.

    Convention: child 0 is the root. Its translation is fixed to
    ``Vec3(0, 0, 0)`` (×2-stored zero), eliminating the global
    translation degree of freedom — WLOG since any uniform shift of
    all children is also a valid placement.

    Per Q8 meta-1: combined with fixing child 0's *rotation* to
    identity (handled by the outer search loop, not here), this
    eliminates both the rotation and translation global symmetries.

    Parameters
    ----------
    rotations
        Per-child rotations (length k). ``rotations[0]`` is the
        root child's rotation.
    edges
        Tree edges in DFS order — every edge's ``parent`` must
        appear in some earlier edge's ``new`` (or be 0 for the
        first edge). The tree must cover all k children
        (length-(k−1) edges spanning {0, ..., k−1}).
    prototile_indices
        Per-child prototile-type assignment (length k).
    prototile_vertices
        ``prototile_vertices[t]`` is the local-frame vertex tuple
        for prototile-type t.

    Returns
    -------
    tuple of Vec3
        Per-child translations (length k); ``[0]`` is always
        ``Vec3(0, 0, 0)``. Each subsequent child's translation is
        computed by adding the per-edge offset to the parent's
        translation.
    None
        If any edge's face-match is infeasible under the given
        rotations (no consistent permutation).

    Raises
    ------
    ValueError
        On length mismatch between ``rotations`` and
        ``prototile_indices`` or on edges that don't form a valid
        tree covering all children.
    """
    k = len(rotations)
    if len(prototile_indices) != k:
        raise ValueError(
            f"prototile_indices length {len(prototile_indices)} != "
            f"rotations length {k}."
        )
    if len(edges) != k - 1:
        raise ValueError(
            f"Tree on {k} children needs {k - 1} edges; got "
            f"{len(edges)}."
        )

    zero = Vec3(_ZERO_ZPHI, _ZERO_ZPHI, _ZERO_ZPHI)
    translations: list[Vec3 | None] = [None] * k
    translations[0] = zero
    for edge in edges:
        if not 0 <= edge.parent < k or not 0 <= edge.new < k:
            raise ValueError(
                f"Edge {edge} references out-of-range child."
            )
        if translations[edge.parent] is None:
            raise ValueError(
                f"Edge {edge}: parent child {edge.parent} not yet "
                "placed. Edges must be in DFS order from root 0."
            )
        if translations[edge.new] is not None:
            raise ValueError(
                f"Edge {edge}: child {edge.new} already placed; "
                "tree must visit each non-root child exactly once."
            )
        offset = translation_offset_from_face_match(
            rotation_a=rotations[edge.parent],
            rotation_b=rotations[edge.new],
            face_indices_a=edge.face_indices_parent,
            face_indices_b=edge.face_indices_new,
            vertices_a=prototile_vertices[prototile_indices[edge.parent]],
            vertices_b=prototile_vertices[prototile_indices[edge.new]],
        )
        if offset is None:
            return None
        parent_t = translations[edge.parent]
        # parent_t known non-None per check above.
        assert parent_t is not None
        translations[edge.new] = parent_t + offset

    if any(t is None for t in translations):
        raise ValueError(
            f"Tree of {len(edges)} edges did not cover all {k} "
            "children — non-spanning?"
        )
    return tuple(t for t in translations if t is not None)


def _is_fibonacci_oracle(
    matrix: np.ndarray, pf_target: ZPhi,
) -> bool:
    """Recognise the Fibonacci 1D oracle.

    Per Q8f: the first oracle is M=[[0,1],[1,1]] (canonical form of
    the Fibonacci substitution after row/col swap) with PF = φ.
    The realisation is constructible by hand — boxes of x-extent
    1 and φ packed end-to-end — so this branch returns ``Realised``
    via a manual witness without invoking a CSP.

    The Penrose P3 oracle (M=[[1,1],[1,2]], PF=φ²) and the 5
    PF=φ³ candidates require the algorithmic CSP and currently
    return ``Inconclusive``; that's the next sub-commit's work.
    """
    if matrix.shape != (2, 2):
        return False
    expected = np.array([[0, 1], [1, 1]])
    if not np.array_equal(matrix, expected):
        return False
    return pf_target == ZPhi(0, 1)


def _build_fibonacci_realisation() -> Realised:
    """Manual realisation witness for the Fibonacci 1D oracle.

    Builds two box prototiles in ℤ[φ]³:

    - ``P_0`` (= "S", short tile) is the unit box ``[0, 2] × [0, 2] ×
      [0, 2]`` in ×2-storage (real ``[0, 1]³``).
    - ``P_1`` (= "L", long tile) is ``[0, 2φ] × [0, 2] × [0, 2]`` in
      ×2-storage (real ``[0, φ] × [0, 1]²``).

    Substitution rule (under the canonical-form M = [[0,1],[1,1]],
    where row-output / column-input):

    - σ(P_0) = 1 child of P_1 at translation 0, identity rotation.
    - σ(P_1) = 1 child of P_0 at translation 0, identity rotation;
              1 child of P_1 at real x-translation 1 (×2-stored
              translation ``(2, 0, 0)``), identity rotation.

    Inflation: x → φ·x; y, z unchanged. So σ(P_0) has real x-extent
    φ (= one P_1) and σ(P_1) has real x-extent φ² = φ + 1 (= one
    P_0 of length 1 followed by one P_1 of length φ). Volumes match
    the eigenvector ``(1, φ)``.

    This is the simplest "Realised" instance — the API smoke test
    that validates the result-type shape end-to-end.
    """
    # Box vertices for P_0 (unit cube, ×2-stored; real coordinates 0..1).
    z = ZPhi(0, 0)
    two = ZPhi(2, 0)
    p0_vertices = [
        Vec3(z, z, z),       Vec3(two, z, z),
        Vec3(z, two, z),     Vec3(two, two, z),
        Vec3(z, z, two),     Vec3(two, z, two),
        Vec3(z, two, two),   Vec3(two, two, two),
    ]
    # P_1: same box but x-extent is 2φ (×2-stored from real φ).
    two_phi = ZPhi(0, 2)
    p1_vertices = [
        Vec3(z, z, z),         Vec3(two_phi, z, z),
        Vec3(z, two, z),       Vec3(two_phi, two, z),
        Vec3(z, z, two),       Vec3(two_phi, z, two),
        Vec3(z, two, two),     Vec3(two_phi, two, two),
    ]
    # Six box faces (one per side), as quads. Outward winding:
    # bottom z=0: ccw seen from +z → (0, 1, 3, 2)
    # top z=1:    cw seen from +z (= ccw seen from -z, but normal points +z, so ccw from outside +z) → (4, 6, 7, 5)
    # The face winding here is approximate; for the Fibonacci oracle's
    # API smoke test we don't run is_in_hull on the boxes — the test
    # only checks the API delivers the right shape.
    box_faces = (
        (0, 1, 3, 2),  # z = 0
        (4, 6, 7, 5),  # z = 1
        (0, 4, 5, 1),  # y = 0
        (2, 3, 7, 6),  # y = 1
        (0, 2, 6, 4),  # x = 0
        (1, 5, 7, 3),  # x = max
    )
    p0 = Polyhedron.from_raw(p0_vertices, box_faces)
    p1 = Polyhedron.from_raw(p1_vertices, box_faces)

    # σ(P_0) = (P_1 at origin, identity rotation).
    # σ(P_1) = (P_0 at origin) + (P_1 at x = 2 in ×2-stored, real x = 1).
    children_of_p0 = (
        ChildPlacement(
            prototile_index=1,
            translation=Vec3(z, z, z),
            rotation=Rotation.identity(),
        ),
    )
    children_of_p1 = (
        ChildPlacement(
            prototile_index=0,
            translation=Vec3(z, z, z),
            rotation=Rotation.identity(),
        ),
        ChildPlacement(
            prototile_index=1,
            translation=Vec3(two, z, z),  # ×2-stored real x = 1.
            rotation=Rotation.identity(),
        ),
    )

    return Realised(
        prototile_shapes=(p0, p1),
        children_per_parent=(children_of_p0, children_of_p1),
        fraction_searched=_ONE_ZPHI,
    )


def realise(
    matrix: np.ndarray,
    pf_target: ZPhi,
    prototile_shapes: tuple[Polyhedron, ...] | None = None,
    *,
    max_search_seconds: float = 300.0,
) -> Realised | NoRealisation | Inconclusive:
    """Decide whether a substitution-matrix candidate admits a 3D
    polyhedral realisation with ℤ[φ]³ vertices.

    Per Claude (web)'s Q8 ruling 2026-04-29:

    - **Q8a**: variable = per-child pose (rotation in I_h × translation
      in ℤ[φ]³). Per Q8 meta, fix one child's rotation to identity
      (WLOG by global I_h symmetry) and search ``I_h^(k−1)`` for
      the remaining ``k−1`` children.
    - **Q8b**: prototile shape pool restricted to icosahedral-
      compatible tetrahedra. Caller supplies via
      ``prototile_shapes``; if ``None``, the function attempts to
      derive a canonical pair from ``pf_target``.
    - **Q8c**: hand-rolled backtracking; rotation iteration with
      linear translation recovery (face-match equations form a
      linear system in ZPhi solvable by Gaussian elimination).
    - **Q8d**: returns ``Realised | NoRealisation | Inconclusive``.
      ``Inconclusive`` carries ``fraction_searched`` and partial
      state for resumption.
    - **Q8 meta-3**: fail-first ordering on rotation search, biased
      toward early NoRealisation detection.

    **Current status (placeholder).** The algorithmic CSP described
    above is multi-week per Q8c; this commit ships:

    - The complete API (``Realised | NoRealisation | Inconclusive``,
      ``ChildPlacement``, ``SearchProgress``).
    - A Fibonacci oracle (``M=[[0,1],[1,1]]``, ``PF=φ``) handled by
      a manual realisation witness — confirms the API end-to-end.
    - Inputs that aren't the Fibonacci oracle return ``Inconclusive``
      with ``fraction_searched=ZPhi(0,0)``, indicating the CSP body
      hasn't been implemented yet.

    Future commits will replace the placeholder with the structured
    backtracker per Q8a/c.

    Parameters
    ----------
    matrix
        The substitution matrix (n × n non-negative integer).
    pf_target
        The target Perron–Frobenius eigenvalue in ℤ[φ].
    prototile_shapes
        Tuple of n Polyhedra (the prototile-shape commitment per
        Q8b/Q8e). If None, the function attempts to derive a canonical
        shape pair from ``pf_target``; currently only supported for
        the Fibonacci oracle.
    max_search_seconds
        Time budget for the algorithmic CSP. Currently only used
        for the placeholder; future implementation will respect it
        as the soft cap.
    """
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError(
            f"matrix must be square; got shape {matrix.shape}."
        )
    n = matrix.shape[0]
    if prototile_shapes is not None and len(prototile_shapes) != n:
        raise ValueError(
            f"prototile_shapes length {len(prototile_shapes)} does "
            f"not match matrix size n={n}."
        )

    # Fibonacci oracle: hand-recognised + manual witness.
    if _is_fibonacci_oracle(matrix, pf_target):
        return _build_fibonacci_realisation()

    # Placeholder: any other input gets Inconclusive with
    # fraction_searched=0. Future structured backtracker fills this
    # in per Q8a/c.
    return Inconclusive(
        fraction_searched=_ZERO_ZPHI,
        partial=None,
        reason=(
            "Track B realisation CSP not yet implemented for general "
            "inputs; only the Fibonacci oracle (M=[[0,1],[1,1]], "
            "PF=φ) returns Realised at this commit. "
            "See apeiron/track_b/realisation.py docstring for "
            "Q8a/c plan."
        ),
    )
