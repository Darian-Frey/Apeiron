"""Substitution rules and substitution matrices over Z[phi]-modules.

A substitution rule σ: R^3 → R^3 is a linear inflation with eigenvalue
λ > 1 (canonically φ² for icosahedral systems per CLAUDE.md §3.4)
such that σ(P) admits a dissection into isometric copies of P for each
prototile P.

This module provides the algebraic layer: construction of a rule, the
non-negative integer substitution matrix that counts prototile-type
multiplicities, primitivity detection (Wielandt bound), and Perron–
Frobenius eigenvalue extraction — both numerically and, where the
eigenvalue lands in Z[phi], exactly.

The algebraic functions operate on ``numpy.ndarray`` matrices directly
so they can be tested against 2D oracle rules (Penrose P3, Fibonacci)
without needing to stand up 3D geometric dissections. ``SubstitutionRule``
carries the full geometric structure for when the downstream
``hierarchy.py`` and ``corona.py`` modules need it.

**Storage note on the inflation matrix.** Unlike ``Rotation`` (which
stores ``2·g`` to absorb the implicit denominator of 2 from CLAUDE.md
§3.3), the inflation is stored as its raw Z[phi] entries. The
canonical inflation ``φ² · I`` has integer-Z[phi] entries directly, so
no scaling convention is needed at the matrix level. When an inflation
is applied to a ×2-stored position ``2v`` the result is ``2·(σv)``,
which is the correct ×2-stored transformed position — no halving step.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Final

import numpy as np

from apeiron.symmetry import ImproperRotation, Mat3, Rotation, Vec3
from apeiron.util import pillar
from apeiron.zphi import ZPhi

__all__ = [
    "PositionedTile",
    "SubstitutionRule",
    "is_primitive",
    "perron_frobenius_eigenvalue",
    "perron_frobenius_in_zphi",
    "substitution_matrix",
]


_PHI_FLOAT: Final[float] = (1.0 + math.sqrt(5.0)) / 2.0


@dataclass(frozen=True, slots=True)
class PositionedTile:
    """A positioned copy of a prototile, used in dissection specs.

    ``prototile_index`` selects which prototile in the alphabet this
    tile is a copy of. ``translation`` is the position of the copy's
    origin in ×2-storage form (CLAUDE.md §3.2). ``rotation`` is an
    isometry of R^3: either a ``Rotation`` (element of I, |I|=60,
    orientation-preserving) or an ``ImproperRotation`` (element of
    I_h \\ I, orientation-reversing). Track A's Danzer ABCK rule
    requires the I_h type per CLAUDE.md §2.1 (Claude (web) ruling
    2026-04-28).
    """

    prototile_index: int
    translation: Vec3
    rotation: Rotation | ImproperRotation


@dataclass(frozen=True, slots=True)
class SubstitutionRule:
    """A substitution rule on an alphabet of ``n_prototiles`` tiles.

    Fields
    ------
    n_prototiles : int
        Alphabet size, ``≥ 1``.
    inflation : Mat3
        The Z[phi]-linear inflation σ as a raw ``Mat3[ZPhi]`` — not
        under the ×2 storage convention of ``Rotation`` (see module
        docstring). For icosahedral substitution tilings the canonical
        choice is ``φ² · I``.
    dissections : tuple of tuples of PositionedTile
        ``dissections[i]`` is the multiset of positioned copies that
        make up ``σ(prototile_i)``. Length must equal ``n_prototiles``.
        Every child's ``prototile_index`` must be in
        ``range(n_prototiles)``.
    """

    n_prototiles: int
    inflation: Mat3
    dissections: tuple[tuple[PositionedTile, ...], ...]

    def __post_init__(self) -> None:
        if self.n_prototiles < 1:
            raise ValueError(
                f"n_prototiles must be ≥ 1; got {self.n_prototiles}."
            )
        if len(self.dissections) != self.n_prototiles:
            raise ValueError(
                f"dissections has length {len(self.dissections)} but "
                f"n_prototiles is {self.n_prototiles}."
            )
        for i, dissection in enumerate(self.dissections):
            for j, child in enumerate(dissection):
                if not 0 <= child.prototile_index < self.n_prototiles:
                    raise ValueError(
                        f"dissections[{i}][{j}].prototile_index = "
                        f"{child.prototile_index} is outside "
                        f"[0, {self.n_prototiles})."
                    )


def substitution_matrix(rule: SubstitutionRule) -> np.ndarray:
    """Return the ``n × n`` non-negative integer substitution matrix.

    Convention: ``M[j, i]`` is the number of type-``j`` copies in
    ``σ(prototile_i)``. Equivalently, the ``i``-th column is the
    vector of child-type counts produced by one inflation of
    prototile ``i``.

    The choice of column convention makes ``M`` act on the left on
    type-count column vectors: if ``v`` is the vector of tile counts
    in a pattern, ``M @ v`` is the vector of tile counts after one
    substitution. Primitivity and the Perron–Frobenius eigenvalue are
    transpose-invariant, so this convention is immaterial for those
    queries.
    """
    n = rule.n_prototiles
    M = np.zeros((n, n), dtype=int)
    for i, dissection in enumerate(rule.dissections):
        for child in dissection:
            M[child.prototile_index, i] += 1
    return M


@pillar(1)
def is_primitive(M: np.ndarray, *, max_power: int | None = None) -> bool:
    """Pillar 1: substitution exists and is primitive.

    Return ``True`` iff some power of ``M`` has all-positive entries.

    Uses Wielandt's bound: for an ``n × n`` non-negative matrix, if
    primitive, then ``M ** k`` has all-positive entries for some
    ``k ≤ (n − 1)² + 1``. Default ``max_power`` uses this bound;
    override for pathological cases or for exploratory use.

    The input must be a square two-dimensional array with non-negative
    integer or float entries. No exactness check is performed here
    (numpy integer arithmetic is exact; float rounding in the product
    could, in principle, turn a strictly-positive entry into zero, but
    the matrices we care about — substitution matrices — are integer-
    valued and numpy preserves them).
    """
    if M.ndim != 2 or M.shape[0] != M.shape[1]:
        raise ValueError(
            f"is_primitive requires a square matrix; got shape {M.shape}."
        )
    n = M.shape[0]
    if max_power is None:
        max_power = (n - 1) ** 2 + 1
    if max_power < 1:
        raise ValueError(f"max_power must be ≥ 1; got {max_power}.")
    power = M.copy()
    for _ in range(max_power):
        if np.all(power > 0):
            return True
        power = power @ M
    return False


def perron_frobenius_eigenvalue(M: np.ndarray) -> float:
    """Return the Perron–Frobenius eigenvalue of ``M`` as a float.

    For a primitive non-negative matrix, the PF eigenvalue is unique,
    simple, real, strictly positive, and strictly dominant. This
    function returns the largest real positive eigenvalue; it raises
    ``ValueError`` if ``M`` has no positive real eigenvalue (which
    indicates the matrix is not in the class where PF is meaningful).

    The return type is ``float``; for cases where the PF eigenvalue
    lies in Z[phi] exactly, use ``perron_frobenius_in_zphi`` to
    recover the exact representation.
    """
    eigenvalues = np.linalg.eigvals(M)
    real_positive = [
        float(e.real)
        for e in eigenvalues
        if abs(e.imag) < 1e-9 and e.real > 0
    ]
    if not real_positive:
        raise ValueError(
            "Matrix has no positive real eigenvalue; Perron–Frobenius "
            "eigenvalue is undefined for this input."
        )
    return max(real_positive)


@pillar(1)
def perron_frobenius_in_zphi(
    M: np.ndarray,
    *,
    max_search: int = 20,
    tolerance: float = 1e-7,
) -> ZPhi | None:
    """Pillar 1: substitution exists and is primitive.

    Return the PF eigenvalue as a ``ZPhi`` if it lies exactly in
    ``Z[phi]``; otherwise ``None``.

    Strategy: compute the PF eigenvalue numerically, enumerate candidate
    ``ZPhi(a, b)`` pairs close to it, and for each numerically-matching
    candidate plug it into the integer characteristic polynomial of
    ``M`` (computed exactly via ``sympy``) evaluated in Z[phi]. A
    candidate is confirmed iff ``char_poly(candidate) = 0`` exactly.

    The exact verification step is what distinguishes genuine Z[phi]
    eigenvalues from irrational PFs that happen to land near the
    Z[phi] lattice (which, since φ is irrational, will happen for
    arbitrarily precise but non-exact approximations). For non-Z[phi]
    PFs the char-poly residual is non-zero in Z[phi] and the candidate
    is rejected.

    ``max_search`` is a per-axis bound on the integer coefficients
    ``(a, b)`` tried; ``tolerance`` is the float matching tolerance for
    the initial numerical filter.
    """
    pf = perron_frobenius_eigenvalue(M)
    char_coeffs = _char_poly_int_coeffs(M)
    for a in range(-max_search, max_search + 1):
        for b in range(-max_search, max_search + 1):
            if abs(a + b * _PHI_FLOAT - pf) >= tolerance:
                continue
            candidate = ZPhi(a, b)
            if _eval_poly_at_zphi(char_coeffs, candidate) == ZPhi(0, 0):
                return candidate
    return None


# -- helpers ----------------------------------------------------------


def _char_poly_int_coeffs(M: np.ndarray) -> list[int]:
    """Return integer coefficients of the characteristic polynomial of
    ``M``, highest degree first.

    Uses ``sympy`` for exact integer arithmetic; assumes ``M`` has
    integer (or near-integer) entries. Output: list of ``n + 1``
    Python ints.
    """
    import sympy as sp

    msp = sp.Matrix(M.tolist())
    x = sp.Symbol("x")
    poly = msp.charpoly(x)
    return [int(c) for c in poly.all_coeffs()]


def _eval_poly_at_zphi(coeffs: Sequence[int], z: ZPhi) -> ZPhi:
    """Horner evaluation of a polynomial (coeffs highest-degree first)
    at a ``ZPhi`` point. Returns a ``ZPhi`` — the polynomial value.
    """
    result = ZPhi(0, 0)
    for c in coeffs:
        result = result * z + c
    return result
