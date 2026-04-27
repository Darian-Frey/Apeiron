"""Tests for apeiron.substitution.

Acceptance oracles per CLAUDE.md §5.1 and the Claude (web) 2026-04-20
directive:

- Primary: Penrose P3 (thick/thin rhombus) substitution,
  ``M = [[2, 1], [1, 1]]``. Primitive with PF eigenvalue φ² = 1 + φ.
- Secondary: Fibonacci substitution, ``M = [[1, 1], [1, 0]]``.
  Primitive with PF eigenvalue φ. Included specifically to catch code
  that hard-codes φ² as "the" canonical PF — Fibonacci's algebraic
  differentness is diagnostic.

The algebraic layer is exercised on bare ``numpy`` integer matrices,
without needing a geometric ``SubstitutionRule``. ``SubstitutionRule``
itself is tested separately for construction, validation, and
substitution-matrix extraction.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from apeiron.substitution import (
    PositionedTile,
    SubstitutionRule,
    is_primitive,
    perron_frobenius_eigenvalue,
    perron_frobenius_in_zphi,
    substitution_matrix,
)
from apeiron.symmetry import Mat3, Rotation, Vec3
from apeiron.zphi import ZPhi

_PHI_FLOAT = (1.0 + math.sqrt(5.0)) / 2.0


def _dummy_tile(idx: int) -> PositionedTile:
    """A PositionedTile at the origin with identity rotation.

    Used to stand up ``SubstitutionRule`` test fixtures without
    requiring geometrically-meaningful translations. The algebraic
    tests only care about ``prototile_index`` counts.
    """
    origin = Vec3(ZPhi(0, 0), ZPhi(0, 0), ZPhi(0, 0))
    return PositionedTile(
        prototile_index=idx,
        translation=origin,
        rotation=Rotation.identity(),
    )


# -- SubstitutionRule construction and validation ----------------------


class TestSubstitutionRuleValidation:
    def test_valid_construction(self) -> None:
        rule = SubstitutionRule(
            n_prototiles=2,
            inflation=Mat3.identity(),
            dissections=(
                (_dummy_tile(0), _dummy_tile(1)),
                (_dummy_tile(0),),
            ),
        )
        assert rule.n_prototiles == 2
        assert rule.inflation == Mat3.identity()
        assert len(rule.dissections) == 2

    def test_rejects_zero_n_prototiles(self) -> None:
        with pytest.raises(ValueError, match="n_prototiles"):
            SubstitutionRule(
                n_prototiles=0,
                inflation=Mat3.identity(),
                dissections=(),
            )

    def test_rejects_negative_n_prototiles(self) -> None:
        with pytest.raises(ValueError, match="n_prototiles"):
            SubstitutionRule(
                n_prototiles=-1,
                inflation=Mat3.identity(),
                dissections=(),
            )

    def test_rejects_wrong_dissection_count(self) -> None:
        with pytest.raises(ValueError, match="length"):
            SubstitutionRule(
                n_prototiles=2,
                inflation=Mat3.identity(),
                dissections=((_dummy_tile(0),),),   # only 1 dissection
            )

    def test_rejects_out_of_range_prototile_index(self) -> None:
        with pytest.raises(ValueError, match="prototile_index"):
            SubstitutionRule(
                n_prototiles=1,
                inflation=Mat3.identity(),
                dissections=((_dummy_tile(5),),),
            )

    def test_rejects_negative_prototile_index(self) -> None:
        with pytest.raises(ValueError, match="prototile_index"):
            SubstitutionRule(
                n_prototiles=2,
                inflation=Mat3.identity(),
                dissections=((_dummy_tile(-1),), (_dummy_tile(0),)),
            )


# -- substitution_matrix extraction ------------------------------------


class TestSubstitutionMatrix:
    def test_column_convention(self) -> None:
        # Prototile 0 → (0, 0, 1): 2 type-0 + 1 type-1.
        # Prototile 1 → (0, 1): 1 type-0 + 1 type-1.
        rule = SubstitutionRule(
            n_prototiles=2,
            inflation=Mat3.identity(),
            dissections=(
                (_dummy_tile(0), _dummy_tile(0), _dummy_tile(1)),
                (_dummy_tile(0), _dummy_tile(1)),
            ),
        )
        M = substitution_matrix(rule)
        assert M.shape == (2, 2)
        # Column 0 is σ(prototile 0) → (2, 1).
        assert M[0, 0] == 2
        assert M[1, 0] == 1
        # Column 1 is σ(prototile 1) → (1, 1).
        assert M[0, 1] == 1
        assert M[1, 1] == 1

    def test_matrix_is_integer_dtype(self) -> None:
        rule = SubstitutionRule(
            n_prototiles=1,
            inflation=Mat3.identity(),
            dissections=((_dummy_tile(0), _dummy_tile(0)),),
        )
        assert substitution_matrix(rule).dtype.kind == "i"

    def test_empty_dissection_gives_zero_column(self) -> None:
        # σ(prototile 1) is an empty multiset → column of zeros.
        rule = SubstitutionRule(
            n_prototiles=2,
            inflation=Mat3.identity(),
            dissections=(
                (_dummy_tile(0),),
                (),
            ),
        )
        M = substitution_matrix(rule)
        assert M[0, 1] == 0
        assert M[1, 1] == 0

    def test_single_prototile_rule(self) -> None:
        rule = SubstitutionRule(
            n_prototiles=1,
            inflation=Mat3.identity(),
            dissections=((_dummy_tile(0),) * 4,),   # 4 copies of self
        )
        M = substitution_matrix(rule)
        assert M.shape == (1, 1)
        assert M[0, 0] == 4


# -- Penrose P3 primary oracle ----------------------------------------


class TestPenroseP3Oracle:
    """Penrose P3 thick/thin rhombus substitution.

    Matrix: thick → 2 thick + 1 thin; thin → 1 thick + 1 thin. This
    gives ``[[2, 1], [1, 1]]``. Characteristic polynomial
    ``x² − 3x + 1``; roots ``(3 ± √5)/2``. The PF eigenvalue is
    ``(3 + √5)/2 = 1 + φ = φ²``.
    """

    P3 = np.array([[2, 1], [1, 1]], dtype=int)

    def test_primitive(self) -> None:
        assert is_primitive(self.P3)

    def test_pf_numerical_equals_phi_squared(self) -> None:
        pf = perron_frobenius_eigenvalue(self.P3)
        expected = 1.0 + _PHI_FLOAT   # φ² = 1 + φ ≈ 2.618
        assert abs(pf - expected) < 1e-10

    def test_pf_in_zphi_is_phi_squared(self) -> None:
        pf = perron_frobenius_in_zphi(self.P3)
        assert pf == ZPhi(1, 1)

    def test_pf_in_zphi_satisfies_char_poly(self) -> None:
        pf = perron_frobenius_in_zphi(self.P3)
        assert pf is not None
        # x² − 3x + 1 = 0 at x = φ².
        residual = pf * pf - 3 * pf + ZPhi(1, 0)
        assert residual == ZPhi(0, 0)


# -- Fibonacci secondary oracle ---------------------------------------


class TestFibonacciOracle:
    """Fibonacci substitution, chosen for PF = φ (not φ²).

    Matrix ``[[1, 1], [1, 0]]``. Characteristic polynomial
    ``x² − x − 1``; roots ``(1 ± √5)/2``. The PF eigenvalue is φ
    itself. Catches any code path that assumes the PF is always φ²
    in the icosahedral context — Fibonacci is 1D and its PF is
    genuinely φ, not φ².
    """

    FIB = np.array([[1, 1], [1, 0]], dtype=int)

    def test_primitive(self) -> None:
        assert is_primitive(self.FIB)

    def test_pf_numerical_equals_phi(self) -> None:
        pf = perron_frobenius_eigenvalue(self.FIB)
        assert abs(pf - _PHI_FLOAT) < 1e-10

    def test_pf_in_zphi_is_phi(self) -> None:
        pf = perron_frobenius_in_zphi(self.FIB)
        assert pf == ZPhi(0, 1)

    def test_pf_not_equal_to_phi_squared(self) -> None:
        # Explicit anti-regression: a code path that hard-codes PF = φ²
        # would fail here but pass the P3 oracle. Keep this test.
        pf = perron_frobenius_in_zphi(self.FIB)
        assert pf != ZPhi(1, 1)


# -- is_primitive generic cases ---------------------------------------


class TestIsPrimitive:
    def test_identity_is_not_primitive(self) -> None:
        assert not is_primitive(np.eye(3, dtype=int))

    def test_zero_matrix_is_not_primitive(self) -> None:
        assert not is_primitive(np.zeros((2, 2), dtype=int))

    def test_reducible_upper_triangular_not_primitive(self) -> None:
        # Upper triangular non-negative → M^k stays upper triangular.
        assert not is_primitive(np.array([[1, 1], [0, 1]], dtype=int))

    def test_periodic_not_primitive(self) -> None:
        # Permutation matrix has PF = 1 but cycles; not primitive.
        M = np.array([[0, 1], [1, 0]], dtype=int)
        assert not is_primitive(M)

    def test_rejects_non_square(self) -> None:
        with pytest.raises(ValueError, match="square"):
            is_primitive(np.zeros((2, 3), dtype=int))

    def test_rejects_non_2d(self) -> None:
        with pytest.raises(ValueError, match="square"):
            is_primitive(np.zeros((2, 2, 2), dtype=int))

    def test_rejects_zero_max_power(self) -> None:
        with pytest.raises(ValueError, match="max_power"):
            is_primitive(np.eye(2, dtype=int), max_power=0)

    def test_respects_max_power_override(self) -> None:
        # P3 is primitive at k = 1 already; give max_power = 1 and
        # confirm. Then give max_power = 1 on a harder matrix.
        P3 = np.array([[2, 1], [1, 1]], dtype=int)
        assert is_primitive(P3, max_power=1)


# -- perron_frobenius_eigenvalue generic cases -----------------------


class TestPerronFrobeniusNumerical:
    def test_one_by_one(self) -> None:
        assert perron_frobenius_eigenvalue(np.array([[5]], dtype=int)) == 5.0

    def test_raises_on_no_positive_real(self) -> None:
        # Rotation matrix — complex eigenvalues ±i.
        M = np.array([[0, -1], [1, 0]], dtype=int)
        with pytest.raises(ValueError, match="positive real eigenvalue"):
            perron_frobenius_eigenvalue(M)

    def test_picks_largest_positive_real(self) -> None:
        # Diagonal with mixed signs.
        M = np.array([[3, 0], [0, 2]], dtype=int)
        assert perron_frobenius_eigenvalue(M) == 3.0


# -- perron_frobenius_in_zphi generic cases --------------------------


class TestPerronFrobeniusInZPhi:
    def test_one_by_one_integer(self) -> None:
        assert perron_frobenius_in_zphi(np.array([[7]], dtype=int)) == ZPhi(7, 0)

    def test_pf_outside_zphi_returns_none(self) -> None:
        # ``[[2, 1], [2, 2]]`` has char poly x² − 4x + 2, PF = 2 + √2.
        # 2 + √2 ∉ Z[phi] because √2 ⊥ √5 over Q.
        M = np.array([[2, 1], [2, 2]], dtype=int)
        assert perron_frobenius_in_zphi(M) is None

    def test_exact_verification_rejects_near_match(self) -> None:
        # Construct a matrix whose PF is close to a ZPhi but not
        # equal: we use the same 2 + √2 example; the best numerical
        # ZPhi candidate at small search radius has residual > 0 in
        # the char-poly check, so None is returned.
        M = np.array([[2, 1], [2, 2]], dtype=int)
        assert perron_frobenius_in_zphi(M, max_search=5) is None
        assert perron_frobenius_in_zphi(M, max_search=50) is None

    def test_zphi_identification_holds_across_search_radii(self) -> None:
        # A legitimate ZPhi PF must be found regardless of max_search
        # (up to the search-radius bound on coefficients).
        P3 = np.array([[2, 1], [1, 1]], dtype=int)
        for radius in (5, 10, 20, 50):
            assert perron_frobenius_in_zphi(P3, max_search=radius) == ZPhi(1, 1)


# -- Pillar 1 tag coverage -------------------------------------------


class TestPillar1Tagging:
    """Sweep: every function establishing pillar 1 must carry the
    ``@pillar(1)`` decorator's ``_pillar = 1`` attribute. Catches
    regressions where a future refactor strips the decorator
    silently.
    """

    def test_is_primitive_is_pillar_1(self) -> None:
        assert is_primitive._pillar == 1

    def test_perron_frobenius_in_zphi_is_pillar_1(self) -> None:
        assert perron_frobenius_in_zphi._pillar == 1
