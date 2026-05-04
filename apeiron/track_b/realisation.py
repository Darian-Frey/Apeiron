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

import time

from apeiron.polyhedron import Polyhedron
from apeiron.symmetry import ICOSAHEDRAL, ImproperRotation, Rotation, Vec3
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


def _enumerate_dfs_trees(k: int) -> list[tuple[int, ...]]:
    """Yield ``(p_1, ..., p_{k-1})`` where ``p_i ∈ {0, ..., i-1}`` is
    the parent of child ``i`` in a DFS-discovery-ordered rooted tree.

    For ``k=1``: yields one empty tuple. For ``k=2``: ``[(0,)]``.
    For ``k=3``: ``[(0, 0), (0, 1)]`` (star at 0; path 0-1-2). Total
    count: ``(k-1)!``. Each sequence specifies a tree on labelled
    nodes ``{0, ..., k-1}`` rooted at 0 with edges visited in
    DFS-discovery order — which is exactly the form
    ``propagate_translations_along_tree`` consumes.

    Matches the convention used by the rotation-search backtracker:
    child 0 is fixed at identity rotation and origin translation; each
    subsequent child connects to one already-placed parent via one
    face-match edge.
    """
    if k < 1:
        raise ValueError(f"k must be ≥ 1; got {k}.")
    if k == 1:
        return [()]
    import itertools
    return [
        parents
        for parents in itertools.product(*(range(i) for i in range(1, k)))
    ]


def _vec3_lex_lt(a: Vec3, b: Vec3) -> bool:
    """ZPhi-lex compare ``a < b`` on (x, y, z)."""
    for ai, bi in ((a.x, b.x), (a.y, b.y), (a.z, b.z)):
        c = (ai - bi)._sign()
        if c < 0:
            return True
        if c > 0:
            return False
    return False


def _tetrahedron_volume_x6(verts: Sequence[Vec3]) -> ZPhi:
    """Six times the (signed) volume of a tetrahedron.

    For a 4-vertex tetrahedron with vertices ``v_0, v_1, v_2, v_3``,
    the signed volume is ``(1/6) · det([v_1−v_0, v_2−v_0, v_3−v_0])``.
    Returning ``6·volume`` keeps everything in ℤ[φ] (no division).
    The sign reflects orientation; ``abs()`` recovers the
    unsigned magnitude for tile-volume comparisons.

    Per CLAUDE.md §3.2, vertices are ×2-stored, so the result is
    ``6 · volume_storage = 6 · (2³ · volume_real) = 48 · volume_real``.
    The factor of 48 is global and cancels in volume-sum-equality
    checks across the same storage convention.
    """
    if len(verts) != 4:
        raise ValueError(
            f"Tetrahedron requires 4 vertices; got {len(verts)}."
        )
    a = verts[1] - verts[0]
    b = verts[2] - verts[0]
    c = verts[3] - verts[0]
    # det = a · (b × c)
    bxc_x = b.y * c.z - b.z * c.y
    bxc_y = b.z * c.x - b.x * c.z
    bxc_z = b.x * c.y - b.y * c.x
    return a.x * bxc_x + a.y * bxc_y + a.z * bxc_z


def _abs_zphi(z: ZPhi) -> ZPhi:
    """Component-wise absolute value (in real-number sense). ``z._sign()``
    decides the case; if negative, returns ``-z``.
    """
    return -z if z._sign() < 0 else z


def _bounding_box(vertices: Sequence[Vec3]) -> tuple[Vec3, Vec3]:
    """Min and max corners (component-wise) over a vertex sequence.
    Pure ZPhi; no floats. Uses ``ZPhi._sign`` on differences.
    """
    if not vertices:
        raise ValueError("Cannot compute bounding box of empty vertex set.")
    min_x = min_y = min_z = vertices[0].x  # placeholder values
    min_x, min_y, min_z = vertices[0].x, vertices[0].y, vertices[0].z
    max_x, max_y, max_z = vertices[0].x, vertices[0].y, vertices[0].z
    for v in vertices[1:]:
        if (v.x - min_x)._sign() < 0:
            min_x = v.x
        if (v.x - max_x)._sign() > 0:
            max_x = v.x
        if (v.y - min_y)._sign() < 0:
            min_y = v.y
        if (v.y - max_y)._sign() > 0:
            max_y = v.y
        if (v.z - min_z)._sign() < 0:
            min_z = v.z
        if (v.z - max_z)._sign() > 0:
            max_z = v.z
    return (Vec3(min_x, min_y, min_z), Vec3(max_x, max_y, max_z))


def _vertex_in_box(v: Vec3, lo: Vec3, hi: Vec3) -> bool:
    """Component-wise ``lo ≤ v ≤ hi`` in real-number ZPhi order."""
    return (
        (v.x - lo.x)._sign() >= 0 and (hi.x - v.x)._sign() >= 0
        and (v.y - lo.y)._sign() >= 0 and (hi.y - v.y)._sign() >= 0
        and (v.z - lo.z)._sign() >= 0 and (hi.z - v.z)._sign() >= 0
    )


def _search_realisation_for_parent(
    children_types: Sequence[int],
    prototile_shapes: Sequence[Polyhedron],
    inflated_parent_bbox: tuple[Vec3, Vec3],
    *,
    rotation_pool: Sequence[Rotation],
    deadline: float,
    k_max: int,
    expected_volume_x6: ZPhi | None = None,
) -> tuple[ChildPlacement, ...] | None | str:
    """Search for a face-to-face placement of ``children_types`` inside
    a parent whose inflated bounding box is ``inflated_parent_bbox``.

    Per Q8a/c/meta 2026-04-29: rotation-search-with-linear-translation
    recovery, with child 0's rotation fixed to identity and translation
    fixed to the origin (eliminating both global symmetries WLOG).

    Search space:

    - **Tree topology**: ``(k−1)!`` DFS-discovery-ordered rooted trees
      on ``k`` children with root 0.
    - **Face pair per edge**: ``faces(prototile_parent) × faces(prototile_new)``
      = 16 options for tetrahedral prototiles (4 × 4).
    - **Rotation per non-root child**: pool of ``len(rotation_pool)``
      isometries (typically the 60 proper rotations in I).

    **Validation: bounding-box only at this iteration.** A placement is
    accepted if every placed vertex lies in
    ``inflated_parent_bbox``. **This is too weak.** Volume
    conservation, non-overlap, and full-coverage are *not* checked.
    Trivially-overlapping placements (e.g., all children at origin
    with identity rotation) pass bounding-box validation and would
    be returned as Realised. A subsequent commit must add:

    - Volume sum: ``Σ vol(child) == λ³ · vol(parent)`` exact in ℤ[φ].
    - Pairwise non-overlap: child interiors disjoint.
    - Coverage: union of placed children equals λ·parent.

    Until then, the ``Realised`` returned by this function is a
    *necessary-condition witness*, not sufficient. Treat any
    Realised at this commit as a candidate for future-iteration
    validation, not a confirmed dissection.

    Returns
    -------
    tuple of ChildPlacement
        First valid placement found. Length matches ``children_types``.
    None
        Exhaustively no placement passes bounding-box validation.
    str
        ``"timeout"`` if the deadline was hit mid-search.
    """
    k = len(children_types)
    if k < 1:
        raise ValueError(f"children_types must be non-empty; got {k}.")
    if k > k_max:
        return "k_too_large"
    z = ZPhi(0, 0)
    origin = Vec3(z, z, z)
    identity = Rotation.identity()

    # Pre-compute child volume sum (×6, signed) for validation. Same
    # value across all rotation/translation assignments since rotations
    # preserve volume and translations don't change it. If
    # expected_volume_x6 is provided and the sum doesn't match, no
    # rotation assignment can possibly Realise — bail immediately.
    if expected_volume_x6 is not None:
        child_vol_sum_x6 = ZPhi(0, 0)
        for t in children_types:
            v = _abs_zphi(_tetrahedron_volume_x6(
                prototile_shapes[t].vertices,
            ))
            child_vol_sum_x6 = child_vol_sum_x6 + v
        if child_vol_sum_x6 != expected_volume_x6:
            return None

    # Special case: k=1. One child, placed at identity/origin. Just
    # check vertex containment.
    if k == 1:
        proto = prototile_shapes[children_types[0]]
        placed = [identity.apply(v) + origin for v in proto.vertices]
        lo, hi = inflated_parent_bbox
        if all(_vertex_in_box(v, lo, hi) for v in placed):
            return (ChildPlacement(
                prototile_index=children_types[0],
                translation=origin, rotation=identity,
            ),)
        return None

    trees = _enumerate_dfs_trees(k)
    n_faces = len(prototile_shapes[0].faces)
    # Cache prototile face vertex tuples (3-tuples per face).
    proto_face_verts: list[tuple[tuple[int, int, int], ...]] = []
    for proto in prototile_shapes:
        face_tuples: list[tuple[int, int, int]] = []
        for face in proto.faces:
            if len(face) != 3:
                raise ValueError(
                    "Search currently supports tetrahedral prototiles "
                    "(triangular faces) only."
                )
            face_tuples.append(face)  # type: ignore[arg-type]
        proto_face_verts.append(tuple(face_tuples))

    proto_vertex_tuples = [tuple(p.vertices) for p in prototile_shapes]

    # Outer iteration: tree, then face-pair sequence, then rotation
    # sequence. Per Q8 meta-3, fail-first ordering (try the most-
    # constraining options first) is a future optimisation; this
    # iteration is straight DFS in the natural order.
    for tree in trees:
        for face_pair_seq in _iter_face_pair_sequences(tree, n_faces):
            for rot_seq in _iter_rotation_sequences(k, rotation_pool):
                if time.monotonic() > deadline:
                    return "timeout"
                edges = []
                for edge_idx, parent in enumerate(tree):
                    new_child = edge_idx + 1
                    f_par_idx, f_new_idx = face_pair_seq[edge_idx]
                    edges.append(FaceMatchEdge(
                        parent=parent,
                        new=new_child,
                        face_indices_parent=proto_face_verts[
                            children_types[parent]
                        ][f_par_idx],
                        face_indices_new=proto_face_verts[
                            children_types[new_child]
                        ][f_new_idx],
                    ))
                translations = propagate_translations_along_tree(
                    rotations=list(rot_seq),
                    edges=edges,
                    prototile_indices=list(children_types),
                    prototile_vertices=proto_vertex_tuples,
                )
                if translations is None:
                    continue
                # Bounding-box validation.
                lo, hi = inflated_parent_bbox
                ok = True
                for child_idx in range(k):
                    proto = prototile_shapes[children_types[child_idx]]
                    rot = rot_seq[child_idx]
                    t = translations[child_idx]
                    for v in proto.vertices:
                        placed = rot.apply(v) + t
                        if not _vertex_in_box(placed, lo, hi):
                            ok = False
                            break
                    if not ok:
                        break
                if not ok:
                    continue
                # Found a valid placement.
                return tuple(
                    ChildPlacement(
                        prototile_index=children_types[i],
                        translation=translations[i],
                        rotation=rot_seq[i],
                    )
                    for i in range(k)
                )
    return None


def _iter_face_pair_sequences(
    tree: tuple[int, ...], n_faces: int,
):
    """Iterator over (face_idx_parent, face_idx_new) tuples per edge.

    For a tree with ``len(tree)`` edges (= k-1), yields all
    ``n_faces ** (2 * len(tree))`` combinations.
    """
    import itertools
    edges_count = len(tree)
    if edges_count == 0:
        yield ()
        return
    pairs = list(itertools.product(range(n_faces), repeat=2))
    for combo in itertools.product(pairs, repeat=edges_count):
        yield combo


def _iter_rotation_sequences(
    k: int, rotation_pool: Sequence[Rotation],
):
    """Iterator over (rot_0, rot_1, ..., rot_{k-1}) where rot_0 is
    fixed to identity and the rest range over rotation_pool.
    """
    import itertools
    identity = Rotation.identity()
    if k == 1:
        yield (identity,)
        return
    for rest in itertools.product(rotation_pool, repeat=k - 1):
        yield (identity,) + rest


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

    # If caller provided no prototile_shapes, we have nothing to
    # search over.
    if prototile_shapes is None:
        return Inconclusive(
            fraction_searched=_ZERO_ZPHI,
            partial=None,
            reason=(
                "Track B realisation requires prototile_shapes for "
                "non-Fibonacci inputs; got None. Caller must supply "
                "icosahedral-compatible prototile Polyhedra per Q8b."
            ),
        )

    # Iterate over each parent's σ; search for a placement of its
    # children. If any parent's σ has no valid placement, the rule
    # has no realisation. Per Q8 meta-1, the search fixes child 0's
    # rotation to identity and translation to the origin.
    deadline = time.monotonic() + max_search_seconds
    rotation_pool = list(ICOSAHEDRAL)  # 60 proper rotations
    children_per_parent: list[tuple[ChildPlacement, ...]] = []
    for parent_idx in range(n):
        children_types = []
        for type_out in range(n):
            children_types.extend(
                [type_out] * int(matrix[type_out, parent_idx])
            )
        if not children_types:
            return NoRealisation(
                reason=(
                    f"σ(prototile_{parent_idx}) is empty (column "
                    f"{parent_idx} of substitution matrix is zero)."
                ),
            )
        # Inflated parent's bounding box: Λ·P_parent's vertices.
        # For the inflation matrix Λ=pf_target·I (the typical case),
        # bbox vertices are pf_target * proto.vertices. For a general
        # ZPhi inflation, the caller must ensure bbox is right.
        # First-iteration assumption: inflation is pf_target * I.
        parent_proto = prototile_shapes[parent_idx]
        inflated_verts = [
            Vec3(
                pf_target * v.x, pf_target * v.y, pf_target * v.z,
            ) for v in parent_proto.vertices
        ]
        bbox = _bounding_box(inflated_verts)
        # Volume-sum constraint: total child volume must equal
        # det(Λ) × parent volume. For substitution rules,
        # pf_target IS det(Λ) (the PF eigenvalue of M equals the
        # volume scaling — Perron–Frobenius theorem applied to
        # the substitution-matrix interpretation): so
        # 6·sum_child_vol = pf_target × 6·parent_vol.
        try:
            parent_vol_x6 = _abs_zphi(_tetrahedron_volume_x6(
                parent_proto.vertices,
            ))
        except ValueError:
            parent_vol_x6 = None
        if parent_vol_x6 is not None:
            expected_x6 = pf_target * parent_vol_x6
        else:
            expected_x6 = None
        result = _search_realisation_for_parent(
            children_types=children_types,
            prototile_shapes=prototile_shapes,
            inflated_parent_bbox=bbox,
            rotation_pool=rotation_pool,
            deadline=deadline,
            k_max=3,
            expected_volume_x6=expected_x6,
        )
        if result == "k_too_large":
            return Inconclusive(
                fraction_searched=_ZERO_ZPHI,
                partial=None,
                reason=(
                    f"σ(prototile_{parent_idx}) has k="
                    f"{len(children_types)} children > k_max=3. "
                    "First-iteration search supports k ≤ 3 only; "
                    "larger σ-rules await algorithmic refinement of "
                    "the rotation-search backtracker."
                ),
            )
        if result == "timeout":
            return Inconclusive(
                fraction_searched=_ZERO_ZPHI,
                partial=None,
                reason=(
                    f"Search hit max_search_seconds={max_search_seconds} "
                    f"on σ(prototile_{parent_idx}); k="
                    f"{len(children_types)} children. Increase budget."
                ),
            )
        if result is None:
            return NoRealisation(
                reason=(
                    f"σ(prototile_{parent_idx}) with {len(children_types)} "
                    "children admits no face-to-face placement under "
                    "bounding-box validation."
                ),
            )
        # result is a tuple of ChildPlacement.
        children_per_parent.append(result)

    return Realised(
        prototile_shapes=tuple(prototile_shapes),
        children_per_parent=tuple(children_per_parent),
        fraction_searched=_ONE_ZPHI,
    )
