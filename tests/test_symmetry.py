"""Tests for apeiron.symmetry — icosahedral rotation group I over Z[phi].

Coverage (CLAUDE.md §5.1 acceptance criteria for symmetry.py):

- Vec3 and Mat3 primitives: arithmetic, matmul, transpose, hash.
- The ×2 storage convention and its halve-assertion contract (§3.3).
- Generator orders: ROT_5^5 = ROT_3^3 = ROT_2^2 = identity.
- |I| = 60, identity present, closure under composition, inverse-closure.
- Norm-preservation under the group action (storage-coord form).
- No floating-point entries anywhere in the group.
"""

from __future__ import annotations

import pytest

from apeiron.symmetry import (
    ICOSAHEDRAL,
    ROT_2,
    ROT_3,
    ROT_5,
    ImproperRotation,
    Mat3,
    Rotation,
    Vec3,
    _halve_zphi,
    determinant,
    is_proper,
)
from apeiron.zphi import ZPhi

_Z = ZPhi(0, 0)
_ONE = ZPhi(1, 0)
_TWO = ZPhi(2, 0)
_PHI = ZPhi(0, 1)


# -- Vec3 ---------------------------------------------------------------


class TestVec3:
    def test_construction_and_indexing(self) -> None:
        v = Vec3(ZPhi(1, 0), ZPhi(2, 0), ZPhi(3, 0))
        assert v[0] == ZPhi(1, 0)
        assert v[1] == ZPhi(2, 0)
        assert v[2] == ZPhi(3, 0)

    def test_index_out_of_range(self) -> None:
        v = Vec3(_Z, _Z, _Z)
        with pytest.raises(IndexError):
            v[3]
        with pytest.raises(IndexError):
            v[-1]

    def test_add(self) -> None:
        u = Vec3(ZPhi(1, 0), ZPhi(2, 0), ZPhi(3, 0))
        w = Vec3(ZPhi(4, 0), ZPhi(5, 0), ZPhi(6, 0))
        assert u + w == Vec3(ZPhi(5, 0), ZPhi(7, 0), ZPhi(9, 0))

    def test_sub(self) -> None:
        u = Vec3(ZPhi(5, 0), ZPhi(7, 0), ZPhi(9, 0))
        w = Vec3(ZPhi(1, 0), ZPhi(2, 0), ZPhi(3, 0))
        assert u - w == Vec3(ZPhi(4, 0), ZPhi(5, 0), ZPhi(6, 0))

    def test_neg(self) -> None:
        v = Vec3(ZPhi(1, -1), ZPhi(0, 2), ZPhi(-3, 0))
        assert -v == Vec3(ZPhi(-1, 1), ZPhi(0, -2), ZPhi(3, 0))

    def test_dot(self) -> None:
        u = Vec3(ZPhi(1, 0), ZPhi(2, 0), ZPhi(3, 0))
        w = Vec3(ZPhi(4, 0), ZPhi(-1, 0), ZPhi(2, 0))
        # 1*4 + 2*(-1) + 3*2 = 4 - 2 + 6 = 8
        assert u.dot(w) == ZPhi(8, 0)

    def test_dot_with_phi(self) -> None:
        u = Vec3(_Z, _ONE, _PHI)   # stored form of the axis (0, 1, phi)/2
        # u . u = 1 + phi^2 = 1 + (phi + 1) = 2 + phi
        assert u.dot(u) == ZPhi(2, 1)

    def test_norm_squared(self) -> None:
        v = Vec3(ZPhi(1, 0), ZPhi(1, 0), ZPhi(1, 0))
        assert v.norm_squared() == ZPhi(3, 0)

    def test_equality_and_hash(self) -> None:
        u = Vec3(ZPhi(1, 2), ZPhi(3, 4), ZPhi(5, 6))
        w = Vec3(ZPhi(1, 2), ZPhi(3, 4), ZPhi(5, 6))
        assert u == w
        assert hash(u) == hash(w)
        assert {u, w} == {u}


# -- Mat3 ---------------------------------------------------------------


class TestMat3:
    def test_identity(self) -> None:
        I = Mat3.identity()
        assert I.row(0) == Vec3(_ONE, _Z, _Z)
        assert I.row(1) == Vec3(_Z, _ONE, _Z)
        assert I.row(2) == Vec3(_Z, _Z, _ONE)

    def test_from_ints_length_check(self) -> None:
        with pytest.raises(ValueError):
            Mat3.from_ints(((0, 0),) * 8)  # wrong length

    def test_from_ints_roundtrip(self) -> None:
        M = Mat3.from_ints((
            (1, 0), (2, 0), (3, 0),
            (4, 0), (5, 0), (6, 0),
            (7, 0), (8, 0), (9, 0),
        ))
        assert M.row0 == Vec3(ZPhi(1, 0), ZPhi(2, 0), ZPhi(3, 0))
        assert M.row2 == Vec3(ZPhi(7, 0), ZPhi(8, 0), ZPhi(9, 0))

    def test_row_index_out_of_range(self) -> None:
        I = Mat3.identity()
        with pytest.raises(IndexError):
            I.row(3)

    def test_col(self) -> None:
        M = Mat3.from_ints((
            (1, 0), (2, 0), (3, 0),
            (4, 0), (5, 0), (6, 0),
            (7, 0), (8, 0), (9, 0),
        ))
        assert M.col(0) == Vec3(ZPhi(1, 0), ZPhi(4, 0), ZPhi(7, 0))
        assert M.col(2) == Vec3(ZPhi(3, 0), ZPhi(6, 0), ZPhi(9, 0))

    def test_transpose(self) -> None:
        M = Mat3.from_ints((
            (1, 0), (2, 0), (3, 0),
            (4, 0), (5, 0), (6, 0),
            (7, 0), (8, 0), (9, 0),
        ))
        Mt = M.transpose()
        assert Mt.row0 == Vec3(ZPhi(1, 0), ZPhi(4, 0), ZPhi(7, 0))
        assert Mt.row2 == Vec3(ZPhi(3, 0), ZPhi(6, 0), ZPhi(9, 0))

    def test_transpose_involution(self) -> None:
        M = Mat3.from_ints((
            (1, 2), (3, -1), (0, 4),
            (5, 0), (-2, 3), (7, -1),
            (0, 0), (1, 1), (2, -2),
        ))
        assert M.transpose().transpose() == M

    def test_matmul_identity_vector(self) -> None:
        v = Vec3(ZPhi(2, 3), ZPhi(-1, 0), ZPhi(0, 5))
        assert Mat3.identity() @ v == v

    def test_matmul_identity_matrix(self) -> None:
        M = Mat3.from_ints((
            (1, 0), (2, 0), (3, 0),
            (4, 0), (5, 0), (6, 0),
            (7, 0), (8, 0), (9, 0),
        ))
        assert Mat3.identity() @ M == M
        assert M @ Mat3.identity() == M

    def test_matmul_vector(self) -> None:
        # Cyclic permutation (x, y, z) -> (y, z, x).
        T = Mat3.from_ints((
            (0, 0), (1, 0), (0, 0),
            (0, 0), (0, 0), (1, 0),
            (1, 0), (0, 0), (0, 0),
        ))
        v = Vec3(ZPhi(3, 0), ZPhi(5, 0), ZPhi(7, 0))
        assert T @ v == Vec3(ZPhi(5, 0), ZPhi(7, 0), ZPhi(3, 0))

    def test_matmul_matrix_associativity(self) -> None:
        A = Mat3.from_ints((
            (1, 0), (2, 0), (0, 1),
            (0, 0), (1, 0), (2, 0),
            (0, 1), (0, 0), (1, 0),
        ))
        B = Mat3.from_ints((
            (1, 1), (0, 0), (2, 0),
            (0, 0), (1, 1), (0, 0),
            (2, 0), (0, 0), (1, 1),
        ))
        C = Mat3.from_ints((
            (0, 0), (1, 0), (0, 0),
            (0, 1), (0, 0), (2, 0),
            (1, 0), (0, 0), (0, 1),
        ))
        assert (A @ B) @ C == A @ (B @ C)

    def test_equality_and_hash(self) -> None:
        M = Mat3.from_ints(((1, 2),) * 9)
        N = Mat3.from_ints(((1, 2),) * 9)
        assert M == N
        assert hash(M) == hash(N)


# -- _halve_zphi assertion ----------------------------------------------


class TestHalveAssertion:
    def test_halve_even_components(self) -> None:
        assert _halve_zphi(ZPhi(4, 6)) == ZPhi(2, 3)
        assert _halve_zphi(ZPhi(-2, 0)) == ZPhi(-1, 0)
        assert _halve_zphi(ZPhi(0, 0)) == ZPhi(0, 0)

    def test_halve_odd_a_raises(self) -> None:
        with pytest.raises(AssertionError):
            _halve_zphi(ZPhi(1, 0))
        with pytest.raises(AssertionError):
            _halve_zphi(ZPhi(3, 4))

    def test_halve_odd_b_raises(self) -> None:
        with pytest.raises(AssertionError):
            _halve_zphi(ZPhi(0, 1))
        with pytest.raises(AssertionError):
            _halve_zphi(ZPhi(4, 3))

    def test_halve_negative_odd_raises(self) -> None:
        with pytest.raises(AssertionError):
            _halve_zphi(ZPhi(-3, 0))
        with pytest.raises(AssertionError):
            _halve_zphi(ZPhi(0, -5))

    def test_vec3_halve_propagates(self) -> None:
        bad = Vec3(ZPhi(2, 0), ZPhi(1, 0), ZPhi(4, 0))  # middle component odd
        with pytest.raises(AssertionError):
            bad._halve()

    def test_mat3_halve_propagates(self) -> None:
        bad = Mat3.from_ints((
            (2, 0), (4, 0), (6, 0),
            (8, 0), (1, 0), (2, 0),  # (1, 0) is odd
            (0, 2), (0, 4), (0, 6),
        ))
        with pytest.raises(AssertionError):
            bad._halve()


# -- Rotation identity and generator properties ------------------------


class TestRotationIdentity:
    def test_identity_storage_form(self) -> None:
        # Identity real form is I_3; stored as 2·I_3 = diag(2, 2, 2).
        I = Rotation.identity()
        assert I.matrix == Mat3(
            Vec3(_TWO, _Z, _Z),
            Vec3(_Z, _TWO, _Z),
            Vec3(_Z, _Z, _TWO),
        )

    def test_identity_fixes_vectors(self) -> None:
        I = Rotation.identity()
        v = Vec3(ZPhi(3, 2), ZPhi(-1, 4), ZPhi(0, 5))
        assert I.apply(v) == v

    def test_identity_is_composition_unit(self) -> None:
        I = Rotation.identity()
        assert I.compose(ROT_5) == ROT_5
        assert ROT_5.compose(I) == ROT_5


class TestGeneratorOrders:
    def _power(self, g: Rotation, n: int) -> Rotation:
        result = Rotation.identity()
        for _ in range(n):
            result = result.compose(g)
        return result

    def test_rot_5_has_order_5(self) -> None:
        assert self._power(ROT_5, 5) == Rotation.identity()
        # And smaller powers are non-identity:
        for k in range(1, 5):
            assert self._power(ROT_5, k) != Rotation.identity()

    def test_rot_3_has_order_3(self) -> None:
        assert self._power(ROT_3, 3) == Rotation.identity()
        for k in range(1, 3):
            assert self._power(ROT_3, k) != Rotation.identity()

    def test_rot_2_has_order_2(self) -> None:
        assert self._power(ROT_2, 2) == Rotation.identity()
        assert ROT_2 != Rotation.identity()

    def test_rot_5_fixes_its_axis(self) -> None:
        # Axis (0, 1, phi) in ×2 storage is (0, 2, 2*phi).
        axis = Vec3(_Z, _TWO, ZPhi(0, 2))
        assert ROT_5.apply(axis) == axis

    def test_rot_3_permutes_cyclically(self) -> None:
        # T: (x, y, z) -> (y, z, x). In ×2 storage the stored coords
        # shift the same way.
        v = Vec3(ZPhi(2, 0), ZPhi(4, 0), ZPhi(6, 0))
        assert ROT_3.apply(v) == Vec3(ZPhi(4, 0), ZPhi(6, 0), ZPhi(2, 0))

    def test_rot_2_flips_y_and_z(self) -> None:
        v = Vec3(ZPhi(4, 2), ZPhi(6, 0), ZPhi(2, 4))
        assert ROT_2.apply(v) == Vec3(ZPhi(4, 2), ZPhi(-6, 0), ZPhi(-2, -4))


# -- The icosahedral group ---------------------------------------------


class TestIcosahedralGroup:
    def test_order_60(self) -> None:
        assert len(ICOSAHEDRAL) == 60

    def test_identity_present(self) -> None:
        assert Rotation.identity() in ICOSAHEDRAL

    def test_generators_present(self) -> None:
        assert ROT_5 in ICOSAHEDRAL
        assert ROT_3 in ICOSAHEDRAL
        assert ROT_2 in ICOSAHEDRAL

    def test_all_distinct(self) -> None:
        assert len(set(ICOSAHEDRAL)) == 60

    def test_closure_under_composition(self) -> None:
        group = set(ICOSAHEDRAL)
        for g in ICOSAHEDRAL:
            for h in ICOSAHEDRAL:
                assert g.compose(h) in group

    def test_inverse_closure(self) -> None:
        # For a rotation, inverse = transpose of the underlying real
        # matrix. In ×2 storage this is exactly the transpose of the
        # stored Mat3.
        group = set(ICOSAHEDRAL)
        for g in ICOSAHEDRAL:
            assert g.inverse() in group

    def test_inverse_is_group_inverse(self) -> None:
        I = Rotation.identity()
        for g in ICOSAHEDRAL:
            assert g.compose(g.inverse()) == I
            assert g.inverse().compose(g) == I

    def test_five_fold_axes_count(self) -> None:
        # I has exactly 24 rotations of order 5 (6 axes × 4 non-identity
        # powers each).
        def order(g: Rotation) -> int:
            I = Rotation.identity()
            if g == I:
                return 1
            h = g
            for k in range(2, 61):
                h = h.compose(g)
                if h == I:
                    return k
            raise AssertionError("Group element has order > 60 (impossible).")

        orders = [order(g) for g in ICOSAHEDRAL]
        assert orders.count(1) == 1    # identity
        assert orders.count(2) == 15   # 15 edge-midpoint 2-fold axes
        assert orders.count(3) == 20   # 10 face-pair 3-fold axes × 2 nontrivial powers
        assert orders.count(5) == 24   # 6 vertex-pair 5-fold axes × 4 nontrivial powers
        assert sum(orders) and len(orders) == 60


# -- Norm preservation --------------------------------------------------


class TestNormPreservation:
    def test_preserves_stored_norm_squared(self) -> None:
        # In ×2 storage, |v_stored|^2 = 4*|v_real|^2. Since rotations
        # preserve |v_real|^2, they preserve |v_stored|^2 too.
        panel = [
            Vec3(_TWO, _Z, _Z),                   # real (1, 0, 0)
            Vec3(_Z, _TWO, ZPhi(0, 2)),           # real (0, 1, phi) — vertex
            Vec3(_TWO, ZPhi(0, 2), _Z),           # real (1, phi, 0) — vertex
            Vec3(ZPhi(0, 2), _Z, _TWO),           # real (phi, 0, 1) — vertex
            Vec3(ZPhi(-2, 2), ZPhi(4, 0), _TWO),  # generic ×2-stored vector
            Vec3(_TWO, _TWO, _TWO),               # real (1, 1, 1)
        ]
        for g in ICOSAHEDRAL:
            for v in panel:
                gv = g.apply(v)
                assert gv.norm_squared() == v.norm_squared()


# -- ImproperRotation (I_h) --------------------------------------------


class TestImproperRotation:
    """Orientation-reversing isometries of R^3 — the 60 elements of
    I_h \\ I needed by Track A's Danzer ABCK substitution.
    """

    def test_apply_inverts_then_rotates(self) -> None:
        # Improper(r).apply(v) == r.apply(-v) == -r.apply(v).
        g = ICOSAHEDRAL[5]
        imp = ImproperRotation(g)
        v = Vec3(_TWO, _Z, ZPhi(0, 2))
        assert imp.apply(v) == g.apply(-v)
        assert imp.apply(v) == -g.apply(v)

    def test_apply_to_identity_is_negation(self) -> None:
        # Improper(identity) is point inversion: v ↦ -v.
        imp = ImproperRotation(Rotation.identity())
        v = Vec3(ZPhi(2, 4), ZPhi(-6, 0), ZPhi(0, -2))
        assert imp.apply(v) == -v

    def test_compose_improper_with_proper_is_improper(self) -> None:
        imp = ImproperRotation(ROT_5)
        result = imp.compose(ROT_3)
        assert isinstance(result, ImproperRotation)
        assert result.rotation == ROT_5.compose(ROT_3)

    def test_compose_improper_with_improper_is_proper(self) -> None:
        # Two parity flips cancel; Improper ∘ Improper ∈ I.
        imp_a = ImproperRotation(ROT_5)
        imp_b = ImproperRotation(ROT_3)
        result = imp_a.compose(imp_b)
        assert isinstance(result, Rotation)
        assert result == ROT_5.compose(ROT_3)

    def test_inverse_is_improper(self) -> None:
        imp = ImproperRotation(ROT_5)
        inv = imp.inverse()
        assert isinstance(inv, ImproperRotation)
        # inv ∘ imp == identity (a Rotation in I).
        round_trip = inv.compose(imp)
        assert isinstance(round_trip, Rotation)
        assert round_trip == Rotation.identity()

    def test_matmul_dispatches_like_compose(self) -> None:
        imp = ImproperRotation(ROT_5)
        assert imp @ ROT_3 == imp.compose(ROT_3)
        assert imp @ imp == imp.compose(imp)

    def test_apply_preserves_norm(self) -> None:
        # Reflections and rotations both preserve |v|^2; the I_h
        # action is norm-preserving end-to-end.
        v = Vec3(ZPhi(2, 0), ZPhi(0, 2), ZPhi(-2, 4))
        for g in ICOSAHEDRAL:
            imp = ImproperRotation(g)
            assert imp.apply(v).norm_squared() == v.norm_squared()

    def test_is_frozen_and_hashable(self) -> None:
        a = ImproperRotation(ROT_5)
        b = ImproperRotation(ROT_5)
        assert hash(a) == hash(b)
        assert {a, b} == {a}
        with pytest.raises(Exception):
            a.rotation = ROT_3  # type: ignore[misc]


class TestIsProperAndDeterminant:
    def test_is_proper_true_for_rotation(self) -> None:
        for g in ICOSAHEDRAL:
            assert is_proper(g) is True

    def test_is_proper_false_for_improper(self) -> None:
        for g in ICOSAHEDRAL:
            assert is_proper(ImproperRotation(g)) is False

    def test_determinant_plus_one_for_rotation(self) -> None:
        for g in ICOSAHEDRAL:
            assert determinant(g) == 1

    def test_determinant_minus_one_for_improper(self) -> None:
        for g in ICOSAHEDRAL:
            assert determinant(ImproperRotation(g)) == -1

    def test_rejects_other_types(self) -> None:
        with pytest.raises(TypeError, match="Rotation or ImproperRotation"):
            is_proper("not a rotation")  # type: ignore[arg-type]

    def test_determinant_returns_python_int(self) -> None:
        # Exact integer-valued, not a numpy scalar or float.
        d = determinant(ROT_5)
        assert isinstance(d, int) and not isinstance(d, bool)


# -- No-floats-anywhere invariant --------------------------------------


class TestNoFloats:
    def test_every_entry_is_zphi_with_int_components(self) -> None:
        # Every ZPhi has strict-int a and b components (enforced by
        # ZPhi's __post_init__). Verify the property across the group.
        for g in ICOSAHEDRAL:
            for i in range(3):
                row = g.matrix.row(i)
                for comp in (row.x, row.y, row.z):
                    assert isinstance(comp, ZPhi)
                    assert isinstance(comp.a, int) and not isinstance(comp.a, bool)
                    assert isinstance(comp.b, int) and not isinstance(comp.b, bool)
