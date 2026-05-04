"""Tests for apeiron.track_b.geometric_prefilter — three Q7c
necessary conditions before realisation CSP.
"""

from __future__ import annotations

import numpy as np
import pytest

from apeiron.track_b import (
    PrefilterResult,
    enumerate_primitive_matrices,
    pf_eigenvector_in_zphi,
    prefilter,
    vertex_class_consistent,
)
from apeiron.track_b.geometric_prefilter import (
    dihedral_angle_compatible,
    euler_characteristic_consistent,
)
from apeiron.zphi import ZPhi


class TestPfEigenvectorInZPhi:
    """Filter 1 helper: recover the right PF eigenvector exactly in ZPhi."""

    def test_fibonacci_eigenvector(self) -> None:
        # M = [[0, 1], [1, 1]], PF = φ. Kernel direction (b, λ-a)
        # = (1, φ-0) = (1, φ).
        M = np.array([[0, 1], [1, 1]])
        result = pf_eigenvector_in_zphi(M, ZPhi(0, 1))
        assert result == (ZPhi(1, 0), ZPhi(0, 1))

    def test_penrose_p3_canonical_eigenvector(self) -> None:
        # M = [[1, 1], [1, 2]] (Penrose P3 in canonical form),
        # PF = φ². Kernel direction (b, λ-a) = (1, (1+φ)-1) = (1, φ).
        M = np.array([[1, 1], [1, 2]])
        result = pf_eigenvector_in_zphi(M, ZPhi(1, 1))
        assert result == (ZPhi(1, 0), ZPhi(0, 1))

    def test_phi_cubed_first_candidate(self) -> None:
        # M = [[0, 1], [1, 4]], PF = φ³ = ZPhi(1, 2).
        # Kernel direction (b, λ-a) = (1, 1+2φ - 0) = (1, 1+2φ).
        M = np.array([[0, 1], [1, 4]])
        result = pf_eigenvector_in_zphi(M, ZPhi(1, 2))
        assert result == (ZPhi(1, 0), ZPhi(1, 2))

    def test_eigenvector_satisfies_eigenvalue_equation(self) -> None:
        # For any primitive 2x2 matrix M with PF eigenvalue p,
        # the recovered v satisfies M @ v == p · v in ZPhi.
        for M in enumerate_primitive_matrices(2, ZPhi(1, 2), max_entry=5):
            v = pf_eigenvector_in_zphi(M, ZPhi(1, 2))
            assert v is not None
            # M @ v_real component-wise:
            mv0 = ZPhi(int(M[0, 0]), 0) * v[0] + ZPhi(int(M[0, 1]), 0) * v[1]
            mv1 = ZPhi(int(M[1, 0]), 0) * v[0] + ZPhi(int(M[1, 1]), 0) * v[1]
            # Compare to p · v:
            pv0 = ZPhi(1, 2) * v[0]
            pv1 = ZPhi(1, 2) * v[1]
            assert mv0 == pv0
            assert mv1 == pv1

    def test_n3_via_cross_product(self) -> None:
        # n=3 extension uses cross product of two rows of (M − λI).
        # Verify on a known n=3 case from the algebraic survey.
        M = np.array([[0, 0, 1], [1, 2, 2], [1, 2, 2]])
        # PF = φ³ = ZPhi(1, 2). The eigenvector should satisfy
        # (M − λI) v = 0.
        result = pf_eigenvector_in_zphi(M, ZPhi(1, 2))
        assert result is not None
        # Verify M @ v == λ v.
        for i in range(3):
            mv_i = ZPhi(0, 0)
            for j in range(3):
                mv_i = mv_i + ZPhi(int(M[i, j]), 0) * result[j]
            assert mv_i == ZPhi(1, 2) * result[i]

    def test_right_vs_left_eigenvector_for_asymmetric_M(self) -> None:
        # Bug-prevention: pf_eigenvector_in_zphi returns the RIGHT
        # eigenvector (frequency vector). For asymmetric M the LEFT
        # eigenvector (volume vector) is different. Callers needing
        # volumes pass M.T explicitly.
        M = np.array([[0, 0, 1], [1, 2, 2], [1, 2, 2]])
        right = pf_eigenvector_in_zphi(M, ZPhi(1, 2))
        left = pf_eigenvector_in_zphi(M.T, ZPhi(1, 2))
        assert right is not None
        assert left is not None
        # For this asymmetric M they must NOT be parallel as vectors;
        # if both were the same up to scale, the test is trivially
        # passed by symmetric Ms only.
        # Cross-multiply to detect proportionality: a/c == b/d iff
        # a*d == b*c.
        a, b, c = right
        x, y, z = left
        # Two vectors are parallel iff every cross-product component
        # is zero.
        assert not (
            a * y == x * b and a * z == x * c and b * z == y * c
        )

    def test_higher_n_not_yet_supported(self) -> None:
        # n=4+ still raises NotImplementedError; cross-product
        # shortcut works only for n=3.
        M = np.eye(4, dtype=int)
        with pytest.raises(NotImplementedError, match="n ∈"):
            pf_eigenvector_in_zphi(M, ZPhi(1, 0))

    def test_rejects_non_square(self) -> None:
        M = np.zeros((2, 3), dtype=int)
        with pytest.raises(ValueError, match="square"):
            pf_eigenvector_in_zphi(M, ZPhi(1, 0))


class TestVertexClassConsistent:
    """Filter 1: PF eigenvector in ZPhi with all positive components."""

    def test_fibonacci_passes(self) -> None:
        M = np.array([[0, 1], [1, 1]])
        assert vertex_class_consistent(M, ZPhi(0, 1)) is True

    def test_all_phi_cubed_candidates_pass(self) -> None:
        # Empirically all 5 PF=φ³ n=2 candidates have ZPhi-positive
        # eigenvectors per the earlier survey.
        for M in enumerate_primitive_matrices(2, ZPhi(1, 2), max_entry=5):
            assert vertex_class_consistent(M, ZPhi(1, 2)) is True


class TestStubFilters:
    """Filter 2 (dihedral) and filter 3 (Euler) currently return None
    (= not implemented). The composite passes_all treats None as
    'not eliminating'.
    """

    def test_dihedral_angle_filter_returns_none(self) -> None:
        M = np.array([[0, 1], [1, 1]])
        assert dihedral_angle_compatible(M, ZPhi(0, 1)) is None

    def test_euler_filter_returns_none(self) -> None:
        M = np.array([[0, 1], [1, 1]])
        assert euler_characteristic_consistent(M, ZPhi(0, 1)) is None


class TestPrefilterResult:
    """passes_all logic across (T/F/None) combinations of the three
    filters."""

    def test_passes_all_true_when_filter_1_true_and_others_none(
        self,
    ) -> None:
        r = PrefilterResult(
            vertex_class=True,
            dihedral_angle=None,
            euler=None,
            eigenvector=(ZPhi(1, 0), ZPhi(0, 1)),
        )
        assert r.passes_all is True

    def test_passes_all_false_when_filter_1_false(self) -> None:
        r = PrefilterResult(
            vertex_class=False,
            dihedral_angle=None,
            euler=None,
            eigenvector=None,
        )
        assert r.passes_all is False

    def test_passes_all_false_when_explicit_filter_failure(self) -> None:
        # If filter 2 or 3 is implemented and returns False, gate fails.
        r = PrefilterResult(
            vertex_class=True,
            dihedral_angle=False,
            euler=None,
            eigenvector=(ZPhi(1, 0), ZPhi(0, 1)),
        )
        assert r.passes_all is False
        r = PrefilterResult(
            vertex_class=True,
            dihedral_angle=None,
            euler=False,
            eigenvector=(ZPhi(1, 0), ZPhi(0, 1)),
        )
        assert r.passes_all is False


class TestPrefilterOnTrackBSurvey:
    """Composite check: run prefilter on every n=2 candidate from
    the matrix_search survey.
    """

    def test_all_5_phi_cubed_candidates_pass(self) -> None:
        # The Track B leads at n=2, PF=φ³: all 5 pass filter 1 and
        # the stubs don't eliminate any. Each is a candidate for
        # geometric realisation.
        survivors: list[np.ndarray] = []
        for M in enumerate_primitive_matrices(2, ZPhi(1, 2), max_entry=5):
            result = prefilter(M, ZPhi(1, 2))
            if result.passes_all:
                survivors.append(M)
        assert len(survivors) == 5

    def test_prefilter_returns_eigenvector_for_callers(self) -> None:
        # The realisation CSP needs the eigenvector to seed volume
        # ratios; prefilter exposes it on success.
        M = np.array([[0, 1], [1, 4]])
        result = prefilter(M, ZPhi(1, 2))
        assert result.eigenvector == (ZPhi(1, 0), ZPhi(1, 2))

    def test_n3_pf_phi_cubed_survey_at_max_entry_2(self) -> None:
        # Q9b 2026-04-29: n=3 algebraic survey at PF=φ³, max_entry=2
        # yields 21 primitive matrices. ALL pass filter 1 (positive
        # ZPhi³ eigenvector, vertex-class consistency). This is the
        # empirical fertility check — far from the "zero survivors"
        # gate from Q9c that would trigger a literature deep-dive.
        candidates = list(enumerate_primitive_matrices(
            3, ZPhi(1, 2), max_entry=2,
        ))
        assert len(candidates) == 21
        pass_count = sum(
            1 for M in candidates
            if prefilter(M, ZPhi(1, 2)).passes_all
        )
        assert pass_count == 21
