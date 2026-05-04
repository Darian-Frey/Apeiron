"""Track B matrix search — enumerate primitive non-negative integer
substitution matrices with target PF eigenvalue in ℤ[φ].

Per Claude (web)'s Q7c 2026-04-29 ruling: the algebraic step of
Track B. For an n×n substitution matrix M to be a candidate,
char_poly(M) evaluated at the target ZPhi PF eigenvalue must be
exactly zero. The exhaustive search is over all M ∈ ([0, max_entry])^(n×n)
that pass the primitivity check; we deduplicate by simultaneous
row/column permutation (relabelling prototiles).

Complexity scales as ``(max_entry+1)^(n²) × n!`` for the canonical-
form check. For n=2, max_entry=10 this is ~14641 × 2 → trivially
fast. For n=3 it's ~10^9 × 6 — exhaustive enumeration is borderline
without pruning. n=4 needs structural pruning (the Danzer matrix is
the obvious sanity-check oracle for n=4).
"""

from __future__ import annotations

import itertools
from collections.abc import Iterator

import numpy as np

from apeiron.substitution import (
    _char_poly_int_coeffs,
    _eval_poly_at_zphi,
    is_primitive,
)
from apeiron.zphi import ZPhi

__all__ = ["enumerate_primitive_matrices"]


_ZERO_ZPHI = ZPhi(0, 0)


def _canonical_under_permutation(matrix: np.ndarray) -> tuple[int, ...]:
    """Return the lex-min flattened form of ``P · M · Pᵀ`` over all
    permutation matrices ``P``.

    Two substitution matrices that differ only by simultaneous
    row/column permutation describe the same substitution rule with
    relabelled prototiles, so they're equivalent. The lex-min
    representative serves as the canonical form for deduplication.
    """
    n = matrix.shape[0]
    best: tuple[int, ...] | None = None
    for perm in itertools.permutations(range(n)):
        permuted = matrix[list(perm), :][:, list(perm)]
        flat = tuple(int(x) for x in permuted.flatten())
        if best is None or flat < best:
            best = flat
    assert best is not None
    return best


def enumerate_primitive_matrices(
    n: int,
    pf_target: ZPhi,
    *,
    max_entry: int = 10,
) -> Iterator[np.ndarray]:
    """Yield primitive ``n × n`` non-negative integer matrices with
    Perron–Frobenius eigenvalue exactly ``pf_target`` in ℤ[φ], up to
    simultaneous row/column permutation.

    Each yielded matrix satisfies:

    1. All entries in ``[0, max_entry]``.
    2. ``is_primitive(M)`` per ``apeiron.substitution`` (Wielandt-bounded
       check that some power has all-positive entries).
    3. ``char_poly_M(pf_target) == ZPhi(0, 0)`` exactly. Note this is
       a *necessary but not sufficient* check that ``pf_target`` is
       the *Perron–Frobenius* eigenvalue specifically — it's
       sufficient that ``pf_target`` is an eigenvalue, but for
       primitive matrices the PF eigenvalue is the unique strictly
       largest real eigenvalue (Perron–Frobenius theorem). Caller
       can re-verify if needed via ``perron_frobenius_in_zphi``.
    4. The canonical form under simultaneous row/column permutation
       (i.e., ``M`` is lex-min in its ``S_n``-orbit).

    Parameters
    ----------
    n : int
        Alphabet size (number of prototile types). Start with ``n=2``;
        ``n=3`` and ``n=4`` are tractable but slower.
    pf_target : ZPhi
        Target PF eigenvalue. Common choices: ``ZPhi(0, 1) = φ``,
        ``ZPhi(1, 1) = φ²``, ``ZPhi(1, 2) = φ³``.
    max_entry : int, keyword-only
        Inclusive upper bound for each matrix entry. Default 10.

    Yields
    ------
    numpy.ndarray of shape (n, n), dtype int
        Primitive substitution matrices with the requested PF
        eigenvalue, in canonical form. Order is the natural
        lex-by-flat-tuple iteration order of the search.

    Raises
    ------
    ValueError
        If ``n < 1`` or ``max_entry < 1``.
    """
    if n < 1:
        raise ValueError(f"n must be ≥ 1; got {n}.")
    if max_entry < 1:
        raise ValueError(f"max_entry must be ≥ 1; got {max_entry}.")

    seen_canonical: set[tuple[int, ...]] = set()
    for entries in itertools.product(range(max_entry + 1), repeat=n * n):
        matrix = np.array(entries, dtype=int).reshape(n, n)

        if not is_primitive(matrix):
            continue

        coeffs = _char_poly_int_coeffs(matrix)
        if _eval_poly_at_zphi(coeffs, pf_target) != _ZERO_ZPHI:
            continue

        canonical = _canonical_under_permutation(matrix)
        if canonical in seen_canonical:
            continue
        seen_canonical.add(canonical)
        yield matrix
