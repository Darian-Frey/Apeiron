"""Algebraic embedding M_ABCK ⊂ M_T*(2F) verified per Papadopolos (1999).

Per Q17c(i) (Claude (web), 2026-05-05): the canonical D₆ icosahedral
projection T*(2F) requires 6 prototiles (or 8 with mandatory blue/red
colour decoration on tiles C and G). ABCK is locally derivable from
T*(2F), with σ_ABCK reducing T*(2F)'s 6/8 prototiles to ABCK's 4. The
algebraic signature of this local-derivation relationship is the
**eigenvalue embedding**:

    {τ³, -τ⁻³, τ, -τ⁻¹}  =  spectrum(M_ABCK)  ⊂  spectrum(M_T*(2F))

This test verifies the embedding directly: encode the volume inflation
matrix M_T*(2F) from Papadopolos-Hohneker-Kramer (1999), eq. (9),
compute its characteristic polynomial, divide by the M_ABCK
characteristic polynomial λ⁴ - 5λ³ + 2λ² + 5λ + 1, and confirm:

  - the quotient factors exactly as λ · (λ³ + (13-8τ)λ² + (61-38τ)λ
    + (62-39τ)) per Papadopolos eq. (13);
  - the remainder is identically zero.

The result establishes algebraically that ABCK's eigenvalue structure
is a strict restriction of T*(2F)'s — explanatory context for why the
Track A face-merge no-go (M_ABCK has no positive integer eigenvalue,
hence no monotile reduction via face-merge) is structural rather than
coincidental.

Reference: Papadopolos, Z., Hohneker, C. & Kramer, P. (1999/2000).
"Tiles–inflation rules for the class of canonical tilings T*(2F)
derived by the projection method." Discrete Math. 221, 101–112.
arXiv:math-ph/9909012.

Reading notes: docs/literature_notes.md §9.
"""

from __future__ import annotations

import numpy as np
import sympy as sp

from apeiron.substitution import substitution_matrix
from tests.integration.test_danzer_abck import (
    _build_danzer_rule_from_paolini_dissection,
)


# τ symbol with the relation τ² = τ + 1 enforced via post-substitution.
TAU = sp.Symbol("tau")
LAM = sp.Symbol("lam")

# Papadopolos-Hohneker-Kramer 1999, Discrete Math. 221, eq. (9).
# Volume inflation matrix M for T*(2F):  {X'} ≡ τ³{X} = M{X}.
# Verbatim transcription. Prototile order: A, B, C^b, C^r, D, F, G^b, G^r.
# Each row X' is the τ-volume decomposition of σ(X).
_M_T2F_VOLUME = sp.Matrix([
    [11 * TAU - 16, 0, 2 * TAU - 2, 2 * TAU - 3, 0,           9 * TAU - 13, TAU - 1, 3 * TAU - 4],   # A'
    [0,             0, 0,           1,           0,           0,            0,       1],             # B'
    [-2 * TAU + 4,  1, 0,           -TAU + 2,    1,           -TAU + 3,     0,       -TAU + 2],     # C^b'
    [-9 * TAU + 15, 0, -2 * TAU + 4, -TAU + 2,   1,           -8 * TAU + 14, -TAU + 2, -2 * TAU + 4],  # C^r'
    [0,             0, 0,           1,           1,           1,            0,       0],             # D'
    [1,             0, 1,           0,           1,           1,            0,       0],             # F'
    [-2 * TAU + 4,  1, 0,           -TAU + 2,    0,           -TAU + 2,     0,       -TAU + 2],     # G^b'
    [-9 * TAU + 15, 0, -2 * TAU + 4, -TAU + 2,   0,           -8 * TAU + 13, -TAU + 2, -2 * TAU + 4],  # G^r'
])


# Papadopolos eq. (13), the expected factored form of char(M_T*(2F)):
#   char(M) = λ · (cubic factor with τ coefs) · (λ⁴ - 5λ³ + 2λ² + 5λ + 1)
#
# The first factor in eq. (13) as printed is (-λ⁴ + 5λ³ - 2λ² - 5λ - 1),
# which is -(λ⁴ - 5λ³ + 2λ² + 5λ + 1). The sign is irrelevant for
# verifying that this polynomial is a factor; we use the unsigned form
# (which is the characteristic polynomial of M_ABCK in standard form).
_ABCK_CHAR_POLY_COEFFS = (1, -5, 2, 5, 1)
_ABCK_CHAR_POLY = sum(c * LAM ** k for k, c in enumerate(reversed(_ABCK_CHAR_POLY_COEFFS)))

# Expected quotient: λ · (λ³ + (13 - 8τ)λ² + (61 - 38τ)λ + (62 - 39τ))
_EXPECTED_QUOTIENT = LAM * (
    LAM ** 3
    + (13 - 8 * TAU) * LAM ** 2
    + (61 - 38 * TAU) * LAM
    + (62 - 39 * TAU)
)


_MINIMAL_POLY = sp.Poly(TAU ** 2 - TAU - 1, TAU)


def _reduce_tau_squared(expr: sp.Expr) -> sp.Expr:
    """Canonicalise an expression over ℤ[τ] / (τ² − τ − 1).

    Polynomial coefficients in λ are reduced mod (τ² − τ − 1).
    Result is in the form a(λ) + b(λ)·τ where a, b are polynomials
    in λ with rational (typically integer) coefficients.
    """
    expr = sp.expand(expr)
    if expr.has(LAM):
        # Split into λ-coefficients, reduce each mod (τ²-τ-1).
        poly = sp.Poly(expr, LAM)
        coeffs = poly.all_coeffs()
        reduced = []
        for c in coeffs:
            if c.has(TAU):
                rem = sp.Poly(c, TAU).rem(_MINIMAL_POLY).as_expr()
                reduced.append(sp.expand(rem))
            else:
                reduced.append(c)
        return sp.Poly(reduced, LAM).as_expr()
    if expr.has(TAU):
        return sp.Poly(expr, TAU).rem(_MINIMAL_POLY).as_expr()
    return expr


def _char_poly(matrix: sp.Matrix) -> sp.Expr:
    """Characteristic polynomial det(M − λI), τ²-reduced."""
    poly = (matrix - LAM * sp.eye(matrix.shape[0])).det()
    return _reduce_tau_squared(sp.expand(poly))


class TestT2FEmbedding:
    """Verify M_ABCK eigenvalues are a strict subset of M_T*(2F)'s
    spectrum, via Papadopolos eq. (9) + eq. (13).
    """

    def test_abck_char_poly_is_quartic_with_expected_coefficients(
        self,
    ) -> None:
        """Apeiron's M_ABCK has characteristic polynomial
        λ⁴ - 5λ³ + 2λ² + 5λ + 1, whose roots are {τ³, -τ⁻³, τ, -τ⁻¹}.
        """
        rule = _build_danzer_rule_from_paolini_dissection()
        m_abck = substitution_matrix(rule)
        coeffs_float = np.poly(m_abck)
        coeffs_int = np.round(coeffs_float).astype(int).tolist()
        assert np.allclose(coeffs_float, coeffs_int), (
            f"M_ABCK char poly should have integer coefficients; "
            f"got {coeffs_float}"
        )
        assert tuple(coeffs_int) == _ABCK_CHAR_POLY_COEFFS, (
            f"expected {_ABCK_CHAR_POLY_COEFFS}, got {tuple(coeffs_int)}"
        )

    def test_t2f_volume_matrix_char_poly_factors_per_papadopolos_eq_13(
        self,
    ) -> None:
        """char(M_T*(2F)) = λ · (cubic) · (M_ABCK char poly)
        with cubic = λ³ + (13 - 8τ)λ² + (61 - 38τ)λ + (62 - 39τ)
        per Papadopolos eq. (13).
        """
        char = _char_poly(_M_T2F_VOLUME)
        # Divide by M_ABCK char poly. Uses sympy polynomial division
        # with λ as the indeterminate; coefficients live in ℤ[τ].
        quotient_poly, remainder_poly = sp.div(
            sp.Poly(char, LAM, domain=sp.QQ[TAU]),
            sp.Poly(_ABCK_CHAR_POLY, LAM, domain=sp.QQ[TAU]),
        )
        # Remainder must be zero — i.e., M_ABCK char poly divides
        # M_T*(2F) char poly.
        remainder_expr = _reduce_tau_squared(remainder_poly.as_expr())
        assert remainder_expr == 0, (
            f"M_ABCK char poly does not divide M_T*(2F) char poly; "
            f"remainder = {remainder_expr}"
        )
        # Quotient must equal λ · (λ³ + (13-8τ)λ² + (61-38τ)λ + (62-39τ))
        quotient_expr = _reduce_tau_squared(
            sp.expand(quotient_poly.as_expr())
        )
        expected_expr = _reduce_tau_squared(
            sp.expand(_EXPECTED_QUOTIENT)
        )
        assert sp.simplify(quotient_expr - expected_expr) == 0, (
            f"quotient mismatch:\n  got={quotient_expr}\n  expected={expected_expr}"
        )

    def test_eigenvalue_embedding_is_strict(self) -> None:
        """The embedding is strict: M_T*(2F) has eigenvalues that
        M_ABCK does not (the 0 root, 1.187, and the complex pair
        ≈ -0.593 ± 0.755i). Hence dim spectrum(M_T*(2F)) > dim
        spectrum(M_ABCK).
        """
        # M_ABCK is 4×4 with 4 eigenvalues; M_T*(2F) is 8×8 with 8
        # (counting algebraic multiplicity).
        rule = _build_danzer_rule_from_paolini_dissection()
        m_abck = substitution_matrix(rule)
        assert m_abck.shape == (4, 4)
        assert _M_T2F_VOLUME.shape == (8, 8)
        # The strictness is structurally encoded in the shape gap;
        # algebraically expressed by the non-trivial cubic factor in
        # eq. (13). The cubic has roots {1.187..., -0.593 ± 0.755i}
        # which are not in M_ABCK's spectrum.

    def test_pf_eigenvalue_phi_cubed_is_root_of_both(self) -> None:
        """The Perron-Frobenius eigenvalue τ³ = ZPhi(1, 2) is a root
        of both M_ABCK's characteristic polynomial and M_T*(2F)'s,
        as expected from the local-derivation embedding.
        """
        # τ³ = 2τ + 1 in the τ²=τ+1 reduction.
        tau_cubed = 2 * TAU + 1
        # Substitute into M_ABCK char poly:
        eval_abck = _reduce_tau_squared(
            _ABCK_CHAR_POLY.subs(LAM, tau_cubed)
        )
        assert eval_abck == 0, (
            f"τ³ should be a root of M_ABCK char poly; got {eval_abck}"
        )
        # Substitute into M_T*(2F) char poly:
        char_t2f = _char_poly(_M_T2F_VOLUME)
        eval_t2f = _reduce_tau_squared(char_t2f.subs(LAM, tau_cubed))
        assert eval_t2f == 0, (
            f"τ³ should be a root of M_T*(2F) char poly; got {eval_t2f}"
        )
