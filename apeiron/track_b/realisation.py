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

from dataclasses import dataclass
from typing import Final

import numpy as np

from apeiron.polyhedron import Polyhedron
from apeiron.symmetry import ImproperRotation, Rotation, Vec3
from apeiron.zphi import ZPhi

__all__ = [
    "ChildPlacement",
    "Inconclusive",
    "NoRealisation",
    "Realised",
    "SearchProgress",
    "realise",
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
