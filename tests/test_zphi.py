"""Tests for apeiron.zphi — exact arithmetic over Z[phi].

Coverage targets (CLAUDE.md §5, §7.4): construction and type discipline,
arithmetic operations, the phi**2 = phi + 1 identity, integer powers and
the Fibonacci identity, Galois conjugation, field norm and trace, total
ordering via exact integer arithmetic, hashing, and immutability.

Numerical sanity checks against float approximations are used only as
supplementary cross-checks, never as the primary assertion (CLAUDE.md §7.4).
"""

from __future__ import annotations

import math

import pytest

from apeiron.zphi import ONE, PHI, ZERO, ZPhi

# A float approximation of phi, used ONLY for numerical sanity checks in
# the tests below. It must never appear in a production code path.
_PHI_FLOAT = (1.0 + math.sqrt(5.0)) / 2.0


def _as_float(z: ZPhi) -> float:
    """Float embedding of z, for numerical sanity checks only."""
    return z.a + z.b * _PHI_FLOAT


# -- construction and type discipline -----------------------------------


class TestConstruction:
    def test_basic_fields(self) -> None:
        z = ZPhi(3, -2)
        assert z.a == 3
        assert z.b == -2

    def test_from_int(self) -> None:
        assert ZPhi.from_int(7) == ZPhi(7, 0)
        assert ZPhi.from_int(0) == ZERO

    def test_module_constants(self) -> None:
        assert ZERO == ZPhi(0, 0)
        assert ONE == ZPhi(1, 0)
        assert PHI == ZPhi(0, 1)

    @pytest.mark.parametrize("bad", [1.0, 0.5, "1", None, [1]])
    def test_rejects_non_int_a(self, bad: object) -> None:
        with pytest.raises(TypeError):
            ZPhi(bad, 0)  # type: ignore[arg-type]

    @pytest.mark.parametrize("bad", [1.0, 0.5, "1", None, [1]])
    def test_rejects_non_int_b(self, bad: object) -> None:
        with pytest.raises(TypeError):
            ZPhi(0, bad)  # type: ignore[arg-type]

    def test_rejects_bool(self) -> None:
        # bool is an int subclass in Python; we reject it explicitly.
        with pytest.raises(TypeError):
            ZPhi(True, 0)  # type: ignore[arg-type]
        with pytest.raises(TypeError):
            ZPhi(0, False)  # type: ignore[arg-type]

    def test_from_int_rejects_bool(self) -> None:
        with pytest.raises(TypeError):
            ZPhi.from_int(True)  # type: ignore[arg-type]


# -- immutability -------------------------------------------------------


class TestImmutability:
    def test_frozen(self) -> None:
        z = ZPhi(1, 2)
        with pytest.raises(Exception):  # FrozenInstanceError
            z.a = 5  # type: ignore[misc]

    def test_hashable_in_set(self) -> None:
        s = {ZPhi(1, 2), ZPhi(1, 2), ZPhi(3, 4)}
        assert len(s) == 2

    def test_equal_objects_equal_hashes(self) -> None:
        assert hash(ZPhi(5, -7)) == hash(ZPhi(5, -7))


# -- equality -----------------------------------------------------------


class TestEquality:
    def test_equal(self) -> None:
        assert ZPhi(2, 3) == ZPhi(2, 3)

    def test_not_equal_components(self) -> None:
        assert ZPhi(2, 3) != ZPhi(2, 4)
        assert ZPhi(2, 3) != ZPhi(3, 3)

    def test_int_not_equal_zphi(self) -> None:
        # Strict: ZPhi(5, 0) is not Python int 5.
        assert ZPhi(5, 0) != 5
        assert 5 != ZPhi(5, 0)


# -- addition and subtraction -------------------------------------------


class TestAdditive:
    def test_add(self) -> None:
        assert ZPhi(1, 2) + ZPhi(3, 4) == ZPhi(4, 6)

    def test_add_int_right(self) -> None:
        assert ZPhi(1, 2) + 5 == ZPhi(6, 2)

    def test_add_int_left(self) -> None:
        assert 5 + ZPhi(1, 2) == ZPhi(6, 2)

    def test_sub(self) -> None:
        assert ZPhi(5, 7) - ZPhi(2, 3) == ZPhi(3, 4)

    def test_sub_int_right(self) -> None:
        assert ZPhi(5, 7) - 2 == ZPhi(3, 7)

    def test_sub_int_left(self) -> None:
        assert 5 - ZPhi(1, 2) == ZPhi(4, -2)

    def test_neg(self) -> None:
        assert -ZPhi(3, -4) == ZPhi(-3, 4)

    def test_pos_identity(self) -> None:
        z = ZPhi(3, -4)
        assert +z == z

    def test_additive_identity(self) -> None:
        z = ZPhi(7, -3)
        assert z + ZERO == z
        assert ZERO + z == z

    def test_additive_inverse(self) -> None:
        z = ZPhi(7, -3)
        assert z + (-z) == ZERO


# -- multiplication -----------------------------------------------------


class TestMultiplication:
    def test_phi_squared_identity(self) -> None:
        """The load-bearing identity: phi**2 = phi + 1."""
        assert PHI * PHI == PHI + ONE
        assert PHI * PHI == ZPhi(1, 1)

    def test_general_product(self) -> None:
        # (1 + 2*phi)(3 + 4*phi)
        #   = 3 + 4*phi + 6*phi + 8*phi**2
        #   = 3 + 10*phi + 8*(phi + 1)
        #   = 11 + 18*phi
        assert ZPhi(1, 2) * ZPhi(3, 4) == ZPhi(11, 18)

    def test_mul_int_right(self) -> None:
        assert ZPhi(3, -5) * 4 == ZPhi(12, -20)

    def test_mul_int_left(self) -> None:
        assert 4 * ZPhi(3, -5) == ZPhi(12, -20)

    def test_multiplicative_identity(self) -> None:
        z = ZPhi(7, -3)
        assert z * ONE == z
        assert ONE * z == z

    def test_commutativity(self) -> None:
        x, y = ZPhi(2, 3), ZPhi(-1, 5)
        assert x * y == y * x

    def test_associativity(self) -> None:
        x, y, z = ZPhi(2, 3), ZPhi(-1, 5), ZPhi(4, -2)
        assert (x * y) * z == x * (y * z)

    def test_distributivity(self) -> None:
        x, y, z = ZPhi(2, 3), ZPhi(-1, 5), ZPhi(4, -2)
        assert x * (y + z) == x * y + x * z


# -- powers -------------------------------------------------------------


class TestPower:
    def test_power_zero(self) -> None:
        assert ZPhi(5, -3) ** 0 == ONE

    def test_power_one(self) -> None:
        z = ZPhi(5, -3)
        assert z**1 == z

    def test_phi_power_fibonacci(self) -> None:
        # phi**n = F_{n-1} + F_n * phi, with F_0 = 0, F_1 = 1.
        fib = [0, 1]
        for _ in range(2, 20):
            fib.append(fib[-1] + fib[-2])
        for n in range(1, 20):
            assert PHI**n == ZPhi(fib[n - 1], fib[n])

    def test_negative_power_rejected(self) -> None:
        with pytest.raises(ValueError):
            PHI ** (-1)

    def test_power_rejects_non_int(self) -> None:
        with pytest.raises(TypeError):
            PHI ** 2.0  # type: ignore[operator]


# -- Galois conjugate, norm, trace --------------------------------------


class TestNumberField:
    def test_conjugate_of_phi(self) -> None:
        # conj(phi) = 1 - phi
        assert PHI.conjugate() == ONE - PHI
        assert PHI.conjugate() == ZPhi(1, -1)

    def test_conjugate_involution(self) -> None:
        z = ZPhi(7, -3)
        assert z.conjugate().conjugate() == z

    def test_conjugate_is_homomorphism(self) -> None:
        x, y = ZPhi(2, 3), ZPhi(-1, 5)
        assert (x + y).conjugate() == x.conjugate() + y.conjugate()
        assert (x * y).conjugate() == x.conjugate() * y.conjugate()

    def test_norm_is_product_with_conjugate(self) -> None:
        # N(z) = z * conj(z), which lies in Z (b-component is zero).
        for z in [ZPhi(2, 3), ZPhi(-1, 5), PHI, ONE, ZPhi(4, -7)]:
            product = z * z.conjugate()
            assert product.b == 0
            assert product.a == z.norm()

    def test_norm_multiplicative(self) -> None:
        x, y = ZPhi(2, 3), ZPhi(-1, 5)
        assert (x * y).norm() == x.norm() * y.norm()

    def test_norm_of_phi(self) -> None:
        # N(phi) = 0 + 1*0 - 1 = -1; phi is a unit.
        assert PHI.norm() == -1

    def test_trace(self) -> None:
        assert ZPhi(3, -2).trace() == 2 * 3 + (-2)
        assert PHI.trace() == 1
        assert ONE.trace() == 2

    def test_trace_additive(self) -> None:
        x, y = ZPhi(2, 3), ZPhi(-1, 5)
        assert (x + y).trace() == x.trace() + y.trace()


# -- ordering -----------------------------------------------------------


class TestOrdering:
    def test_sign_zero(self) -> None:
        assert ZERO._sign() == 0

    def test_sign_b_zero(self) -> None:
        assert ZPhi(5, 0)._sign() == 1
        assert ZPhi(-5, 0)._sign() == -1

    def test_sign_a_zero(self) -> None:
        assert PHI._sign() == 1
        assert ZPhi(0, -3)._sign() == -1

    def test_sign_same_signs(self) -> None:
        assert ZPhi(2, 3)._sign() == 1
        assert ZPhi(-2, -3)._sign() == -1

    @pytest.mark.parametrize(
        "a,b,expected",
        [
            (2, -1, 1),   # 2 - phi ≈ 0.382 > 0
            (1, -1, -1),  # 1 - phi ≈ -0.618 < 0
            (-2, 1, -1),  # -2 + phi ≈ -0.382 < 0
            (-1, 1, 1),   # -1 + phi ≈ 0.618 > 0
            (3, -2, -1),  # 3 - 2*phi ≈ -0.236 < 0
            (4, -2, 1),   # 4 - 2*phi ≈ 0.764 > 0
        ],
    )
    def test_sign_opposite_signs(self, a: int, b: int, expected: int) -> None:
        z = ZPhi(a, b)
        assert z._sign() == expected
        # Numerical sanity check (supplementary, per CLAUDE.md §7.4).
        f = _as_float(z)
        assert (f > 0) == (expected > 0)
        assert (f < 0) == (expected < 0)

    def test_lt(self) -> None:
        assert ONE < PHI
        assert PHI < PHI * PHI
        assert ZPhi(1, -1) < ZERO

    def test_le(self) -> None:
        assert ONE <= ONE
        assert ONE <= PHI

    def test_gt(self) -> None:
        assert PHI > ONE
        assert PHI * PHI > PHI

    def test_ge(self) -> None:
        assert PHI >= PHI
        assert PHI >= ONE

    def test_compare_with_int(self) -> None:
        assert PHI > 1
        assert PHI < 2
        assert ZERO <= 0
        assert 0 < PHI
        assert 2 > PHI

    def test_ordering_rejects_non_numeric(self) -> None:
        with pytest.raises(TypeError):
            PHI < "phi"  # type: ignore[operator]


# -- display ------------------------------------------------------------


class TestDisplay:
    def test_repr_roundtrip(self) -> None:
        z = ZPhi(3, -2)
        assert repr(z) == "ZPhi(3, -2)"
        # eval-based roundtrip is not claimed, but repr is informative.

    @pytest.mark.parametrize(
        "z,expected",
        [
            (ZPhi(0, 0), "0"),
            (ZPhi(5, 0), "5"),
            (ZPhi(0, 1), "phi"),
            (ZPhi(0, -1), "-phi"),
            (ZPhi(0, 3), "3*phi"),
            (ZPhi(2, 1), "2 + phi"),
            (ZPhi(2, -1), "2 - phi"),
            (ZPhi(2, 3), "2 + 3*phi"),
            (ZPhi(2, -3), "2 - 3*phi"),
        ],
    )
    def test_str_forms(self, z: ZPhi, expected: str) -> None:
        assert str(z) == expected


# -- mixed-type arithmetic rejection -----------------------------------


class TestMixedTypeGuards:
    def test_add_float_rejected(self) -> None:
        with pytest.raises(TypeError):
            ZPhi(1, 0) + 0.5  # type: ignore[operator]

    def test_sub_float_rejected(self) -> None:
        with pytest.raises(TypeError):
            ZPhi(1, 0) - 0.5  # type: ignore[operator]

    def test_mul_float_rejected(self) -> None:
        with pytest.raises(TypeError):
            ZPhi(1, 0) * 0.5  # type: ignore[operator]

    def test_radd_float_rejected(self) -> None:
        with pytest.raises(TypeError):
            0.5 + ZPhi(1, 0)  # type: ignore[operator]

    def test_rsub_float_rejected(self) -> None:
        with pytest.raises(TypeError):
            0.5 - ZPhi(1, 0)  # type: ignore[operator]
