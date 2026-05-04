"""Tests for apeiron.track_b.matrix_search — primitive substitution
matrix enumeration with target PF eigenvalue in ℤ[φ].

Per Claude (web)'s Q7c 2026-04-29 ruling: Track B's algebraic step.
"""

from __future__ import annotations

import numpy as np
import pytest

from apeiron.substitution import is_primitive, perron_frobenius_in_zphi
from apeiron.track_b import enumerate_primitive_matrices
from apeiron.track_b.matrix_search import _canonical_under_permutation
from apeiron.zphi import ZPhi


class TestEnumeratePrimitiveMatrices:
    """Algebraic survey: enumerate primitive matrices with target PF
    eigenvalues for n=2 and n=3.
    """

    def test_n2_pf_phi_yields_fibonacci(self) -> None:
        # n=2, PF=φ: exactly the Fibonacci substitution [[0,1],[1,1]].
        results = list(enumerate_primitive_matrices(
            2, ZPhi(0, 1), max_entry=5,
        ))
        assert len(results) == 1
        np.testing.assert_array_equal(results[0], [[0, 1], [1, 1]])

    def test_n2_pf_phi_squared_yields_penrose_p3(self) -> None:
        # n=2, PF=φ²: exactly the Penrose P3 thick/thin matrix
        # (in canonical form [[1,1],[1,2]] = relabelled [[2,1],[1,1]]).
        results = list(enumerate_primitive_matrices(
            2, ZPhi(1, 1), max_entry=5,
        ))
        assert len(results) == 1
        np.testing.assert_array_equal(results[0], [[1, 1], [1, 2]])

    def test_n2_pf_phi_cubed_yields_five_candidates(self) -> None:
        # n=2, PF=φ³: 5 distinct primitive matrices with this PF
        # eigenvalue. None of these is ABCK (which is n=4); they're
        # the 2-tile algebraic candidates per Track B.
        results = list(enumerate_primitive_matrices(
            2, ZPhi(1, 2), max_entry=5,
        ))
        assert len(results) == 5
        # Each must satisfy trace=4, det=-1 (the closed-form PF=φ³
        # condition for n=2).
        for matrix in results:
            assert matrix.trace() == 4
            # Use exact 2×2 det formula to avoid float rounding.
            a, b = int(matrix[0, 0]), int(matrix[0, 1])
            c, d = int(matrix[1, 0]), int(matrix[1, 1])
            assert a * d - b * c == -1

    def test_yielded_matrices_are_primitive(self) -> None:
        for matrix in enumerate_primitive_matrices(
            2, ZPhi(1, 2), max_entry=5,
        ):
            assert is_primitive(matrix)

    def test_yielded_matrices_have_target_pf(self) -> None:
        # PF eigenvalue recovered exactly via perron_frobenius_in_zphi.
        for matrix in enumerate_primitive_matrices(
            2, ZPhi(1, 1), max_entry=5,
        ):
            assert perron_frobenius_in_zphi(matrix) == ZPhi(1, 1)

    def test_no_duplicates_under_permutation(self) -> None:
        # No two yielded matrices share a canonical form under
        # simultaneous row/column permutation.
        canonicals = set()
        for matrix in enumerate_primitive_matrices(
            2, ZPhi(1, 2), max_entry=5,
        ):
            canonical = _canonical_under_permutation(matrix)
            assert canonical not in canonicals
            canonicals.add(canonical)

    def test_rejects_invalid_n(self) -> None:
        with pytest.raises(ValueError, match="n must be"):
            list(enumerate_primitive_matrices(0, ZPhi(0, 1), max_entry=5))

    def test_rejects_invalid_max_entry(self) -> None:
        with pytest.raises(ValueError, match="max_entry"):
            list(enumerate_primitive_matrices(2, ZPhi(0, 1), max_entry=0))

    def test_n2_pf_phi_with_small_max_entry_still_finds_fibonacci(
        self,
    ) -> None:
        # Fibonacci's largest entry is 1; max_entry=2 must still find it.
        results = list(enumerate_primitive_matrices(
            2, ZPhi(0, 1), max_entry=2,
        ))
        assert len(results) == 1
        np.testing.assert_array_equal(results[0], [[0, 1], [1, 1]])


class TestCanonicalUnderPermutation:
    """Internal helper: lex-min representative under S_n action."""

    def test_identity_is_canonical(self) -> None:
        # The lex-min permutation of any matrix is its own canonical
        # form iff no relabelling produces a lex-smaller flatten.
        matrix = np.array([[0, 1], [1, 1]])
        canonical = _canonical_under_permutation(matrix)
        assert canonical == (0, 1, 1, 1)

    def test_transposed_pairs_share_canonical_form(self) -> None:
        # Two matrices related by simultaneous row/column swap have
        # the same canonical form.
        a = np.array([[2, 1], [1, 1]])
        b = np.array([[1, 1], [1, 2]])  # rows + cols swapped
        assert _canonical_under_permutation(a) == _canonical_under_permutation(b)

    def test_distinct_matrices_have_distinct_canonical(self) -> None:
        a = np.array([[0, 1], [1, 1]])
        b = np.array([[1, 1], [1, 2]])
        assert _canonical_under_permutation(a) != _canonical_under_permutation(b)
