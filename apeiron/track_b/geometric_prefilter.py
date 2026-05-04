"""Track B geometric pre-filter — three cheap necessary conditions
that an algebraic substitution matrix candidate must satisfy *before*
attempting geometric realisation.

Per Claude (web)'s Q7c 2026-04-29 ruling: filter 1 is fully
implementable at the matrix level; filters 2 and 3 require either
realisation hypothesis or extended algebraic machinery and are
provided as documented stubs returning ``None`` until those
machineries land.

Filter 1 — **Vertex-class consistency**. A 3D icosahedral substitution
tiling with PF eigenvalue ``λ`` has prototile *volume frequencies*
proportional to the right PF eigenvector ``v`` of the substitution
matrix. For the tiling to be expressible exactly in ℤ[φ]³, ``v``
must have ZPhi-valued entries (after suitable normalisation) and
all entries must be strictly positive (the Perron–Frobenius theorem
guarantees positivity for primitive matrices, but the ZPhi-ness is
the matrix-specific check).

Filter 2 — **Dihedral-angle commensurability** (stub). For 3D
icosahedral compatibility, dihedral angles must be in
``{π/2, π/3, 2π/3, π/5, 2π/5, 3π/5, 4π/5}`` and their supplements.
This requires the PF eigenvector's component ratios (= prototile
volume ratios) to factor as cubes of ZPhi-valued edge length ratios
that produce icosahedral angles. The full check needs algebraic
machinery for ZPhi cube roots and angle commensurability that hasn't
been built yet; the stub returns ``None`` to indicate "not checked".

Filter 3 — **Euler characteristic consistency** (stub). For a
prototile that is a topological ball, ``V − E + F = 2``. For
tetrahedra (V=4, E=6, F=4), this is automatic. For general convex
polyhedra the F count must be consistent with the substitution
matrix's implied face structure. Without a candidate realisation
the F count is undetermined; stub returns ``None``.

The composite ``prefilter`` returns a ``PrefilterResult`` with one
field per filter and a ``passes_all`` accessor that treats
unimplemented filters (``None``) as "not eliminating". Only
candidates passing every implemented filter proceed to the
realisation CSP.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final

import numpy as np

from apeiron.zphi import ZPhi

__all__ = [
    "PrefilterResult",
    "dihedral_angle_compatible",
    "euler_characteristic_consistent",
    "pf_eigenvector_in_zphi",
    "prefilter",
    "vertex_class_consistent",
]


_ZERO_ZPHI: Final[ZPhi] = ZPhi(0, 0)


def pf_eigenvector_in_zphi(
    matrix: np.ndarray,
    pf_target: ZPhi,
) -> tuple[ZPhi, ...] | None:
    """Compute the **right** PF eigenvector of ``matrix`` over ℤ[φ],
    or return ``None`` if it isn't ZPhi-expressible.

    **Right vs left eigenvector — semantics warning.** For a
    substitution matrix ``M`` with column convention
    ``M[j, i] = #type-j in σ(prototile_i)``:

    - The **right** eigenvector ``v`` (``M·v = λ·v``) gives the
      *frequency* vector — long-range relative occurrences of each
      prototile in tilings.
    - The **left** eigenvector ``w`` (``M^T·w = λ·w``) gives the
      *volume* vector — vol(P_i) ∝ w_i, since
      ``vol(σ(P_i)) = Σ_j M[j, i]·vol(P_j) = (M^T·w)[i] = λ·w_i``.

    For symmetric ``M`` the two coincide; for asymmetric they
    differ. Callers needing the volume vector should pass
    ``M.T`` (the transpose) to recover the left eigenvector via
    this same function.

    For an n×n integer matrix with eigenvalue ``pf_target ∈ ZPhi``,
    the kernel of ``matrix − pf_target · I`` is a ZPhi-module. We
    construct one non-zero vector in this kernel by Gaussian
    elimination over ℚ(√5), tracking entries in ZPhi exactly.

    For ``n=2``, closed form: if ``matrix[0,1] != 0``, the kernel
    direction is ``(matrix[0,1], pf_target − matrix[0,0])``. If
    that's zero (degenerate), fall through to the general kernel
    computation. (For primitive matrices the eigenspace is
    one-dimensional and one of the two row equations is non-
    degenerate, so this works.)

    For ``n ≥ 3``, currently raises ``NotImplementedError`` —
    Track B's first-pass focus is n=2; extending to higher n is
    a follow-up task.

    Returns
    -------
    tuple of ZPhi
        A non-zero eigenvector over ZPhi. Component sign convention:
        if the natural kernel direction has all components in ZPhi
        with at least one positive (real-valued), returns it;
        otherwise tries the negation. Caller should not rely on a
        canonical scaling.
    None
        If no ZPhi-valued non-zero kernel vector exists. (For
        integer matrices with ZPhi eigenvalue this is rare but
        possible if the natural construction yields all-zero
        candidate components.)
    """
    n = matrix.shape[0]
    if matrix.shape != (n, n):
        raise ValueError(f"Expected square matrix; got {matrix.shape}.")
    if n == 2:
        return _pf_eigenvector_in_zphi_n2(matrix, pf_target)
    if n == 3:
        return _pf_eigenvector_in_zphi_n3(matrix, pf_target)
    raise NotImplementedError(
        "pf_eigenvector_in_zphi supports n ∈ {2, 3}; "
        f"got n={n}. Extend via Gaussian elimination over ℤ[φ] for "
        "higher n as a follow-up."
    )


def _pf_eigenvector_in_zphi_n2(
    matrix: np.ndarray, pf_target: ZPhi,
) -> tuple[ZPhi, ...] | None:
    a = ZPhi(int(matrix[0, 0]), 0)
    b = ZPhi(int(matrix[0, 1]), 0)
    c = ZPhi(int(matrix[1, 0]), 0)
    d = ZPhi(int(matrix[1, 1]), 0)
    # Try the row-0 form first.
    v0 = b
    v1 = pf_target - a
    if v0 != _ZERO_ZPHI or v1 != _ZERO_ZPHI:
        if c * v0 + (d - pf_target) * v1 == _ZERO_ZPHI:
            return (v0, v1)
    # Fall through to row-1 form.
    v0 = pf_target - d
    v1 = c
    if v0 != _ZERO_ZPHI or v1 != _ZERO_ZPHI:
        if (a - pf_target) * v0 + b * v1 == _ZERO_ZPHI:
            return (v0, v1)
    return None


def _pf_eigenvector_in_zphi_n3(
    matrix: np.ndarray, pf_target: ZPhi,
) -> tuple[ZPhi, ...] | None:
    """For a 3×3 matrix M with PF eigenvalue λ ∈ ℤ[φ], the kernel
    of (M − λI) is the 1D PF eigenspace (for primitive M). A non-zero
    eigenvector is the cross product of any two linearly-independent
    rows of (M − λI) — the kernel direction, perpendicular to both
    rows in the row-space sense.

    Try the cross product of rows (0, 1); if zero, try (0, 2); if
    zero, try (1, 2). At least one cross product is non-zero since
    the kernel is 1D (= rank 2 matrix).
    """
    from apeiron.symmetry import Vec3
    # Build (M - λI) as 3 ZPhi rows.
    a_minus_λI = []
    for i in range(3):
        row = []
        for j in range(3):
            entry = ZPhi(int(matrix[i, j]), 0)
            if i == j:
                entry = entry - pf_target
            row.append(entry)
        a_minus_λI.append(Vec3(row[0], row[1], row[2]))

    def _cross(u, v):
        return Vec3(
            u.y * v.z - u.z * v.y,
            u.z * v.x - u.x * v.z,
            u.x * v.y - u.y * v.x,
        )

    def _is_zero(v):
        return v.x == _ZERO_ZPHI and v.y == _ZERO_ZPHI and v.z == _ZERO_ZPHI

    # Try each pair of rows.
    for i, j in [(0, 1), (0, 2), (1, 2)]:
        cand = _cross(a_minus_λI[i], a_minus_λI[j])
        if _is_zero(cand):
            continue
        # Verify: for the candidate to be a kernel vector, every row's
        # dot with cand must be zero.
        ok = True
        for k in range(3):
            if a_minus_λI[k].dot(cand) != _ZERO_ZPHI:
                ok = False
                break
        if ok:
            return (cand.x, cand.y, cand.z)
    return None


def vertex_class_consistent(
    matrix: np.ndarray,
    pf_target: ZPhi,
) -> bool:
    """Filter 1: PF right eigenvector is ZPhi-valued with all
    components positive in the real embedding.

    By the Perron–Frobenius theorem, a primitive non-negative matrix
    has a unique strictly positive eigenvector for its PF eigenvalue;
    we additionally require ZPhi-ness for icosahedral compatibility
    (vertices lie in ℤ[φ]³).

    For ``n=2`` only at present (per ``pf_eigenvector_in_zphi``).
    """
    eigenvector = pf_eigenvector_in_zphi(matrix, pf_target)
    if eigenvector is None:
        return False
    # All components must be strictly positive in real embedding.
    return all(component > _ZERO_ZPHI for component in eigenvector)


def dihedral_angle_compatible(
    matrix: np.ndarray,
    pf_target: ZPhi,
) -> bool | None:
    """Filter 2 (stub): dihedral angles consistent with icosahedral.

    Per Claude (web) Q7c 2026-04-29, the check is: the PF eigenvector's
    component ratios (= prototile volume ratios) factor as cubes of
    ZPhi edge length ratios that produce icosahedral-compatible
    dihedral angles ``{π/2, π/3, 2π/3, π/5, 2π/5, 3π/5, 4π/5}`` and
    supplements.

    The check needs:

    1. ZPhi cube roots — the volume-to-edge-length conversion. ℤ[φ]
       is a Dedekind domain but cube roots aren't generally in it;
       a candidate's volume ratio must factor through a cube
       structure compatible with icosahedral angles.
    2. Dihedral-angle parametrisation in ZPhi — angles via cosine
       (or via the icosahedral-rotation-eigenvalue framework already
       in ``apeiron.symmetry``).
    3. Compatibility check.

    This is a substantial bit of algebraic machinery to build. Until
    it lands the filter returns ``None`` to indicate "not checked",
    and ``PrefilterResult.passes_all`` treats ``None`` as
    "not eliminating".
    """
    return None


def euler_characteristic_consistent(
    matrix: np.ndarray,
    pf_target: ZPhi,
) -> bool | None:
    """Filter 3 (stub): substitution face structure is consistent
    with V − E + F = 2 for prototiles being topological balls.

    For tetrahedral candidates (V=4, E=6, F=4) this is automatic.
    For general convex polyhedra the F count must be consistent with
    the substitution rule's face-adjacency structure inside σ(P_i),
    which is undetermined without a realisation.

    Returns ``None`` until the realisation step provides face counts
    that the matrix-level check can validate.
    """
    return None


@dataclass(frozen=True, slots=True)
class PrefilterResult:
    """Outcome of running all three Q7c geometric prefilters.

    - ``vertex_class``: filter 1 (necessary, fully implemented).
    - ``dihedral_angle``: filter 2 (stub; ``None`` = not checked).
    - ``euler``: filter 3 (stub; ``None`` = not checked).
    - ``eigenvector``: the ZPhi PF right eigenvector (if recovered)
      for downstream callers (e.g., the realisation CSP needs the
      volume ratio to seed vertex placements).

    ``passes_all`` returns True iff every implemented filter
    returned True; ``None`` (unchecked) is treated as
    "not eliminating".
    """

    vertex_class: bool
    dihedral_angle: bool | None
    euler: bool | None
    eigenvector: tuple[ZPhi, ...] | None

    @property
    def passes_all(self) -> bool:
        if self.vertex_class is not True:
            return False
        if self.dihedral_angle is False:
            return False
        if self.euler is False:
            return False
        return True


def prefilter(
    matrix: np.ndarray,
    pf_target: ZPhi,
) -> PrefilterResult:
    """Apply all three Track B geometric prefilters to a candidate
    substitution matrix.

    Cheap; runs in O(n²) ZPhi operations for n=2. Caller proceeds to
    the realisation CSP iff ``result.passes_all is True``.
    """
    eigenvector = pf_eigenvector_in_zphi(matrix, pf_target)
    vertex_class = bool(
        eigenvector is not None
        and all(component > _ZERO_ZPHI for component in eigenvector)
    )
    return PrefilterResult(
        vertex_class=vertex_class,
        dihedral_angle=dihedral_angle_compatible(matrix, pf_target),
        euler=euler_characteristic_consistent(matrix, pf_target),
        eigenvector=eigenvector,
    )
