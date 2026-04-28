"""Icosahedral rotation group I over Z[phi].

I is the order-60 orientation-preserving symmetry group of the regular
icosahedron. See CLAUDE.md §3.3.

**Storage convention.** The natural Cartesian representation of I has
entries in (1/2)Z[phi], not Z[phi] — the 5-fold rotation about the axis
(0, 1, phi) has half-integer entries even after cancelling the sqrt(5)
factors. Following CLAUDE.md §3.3, every rotation is stored as 2*g, i.e.
with a uniform implicit denominator of 2. The denominator is globally
bounded (it never grows under composition), so arithmetic is performed
in pure Z[phi] and the factor of 2 is discharged at a single point via
``_halve_zphi`` — a runtime assertion that both integer components of
the numerator are even before dividing.

**Positions.** Position vectors (stored as 2*v per CLAUDE.md §3.2) act
the same way: ``Rotation.apply`` computes ``self.matrix @ v`` — which
a priori lies in (1/4)Z[phi] — and halves component-wise. The halve
asserts exactness; if it fires, a bug upstream has violated the
I-invariance denominator bound.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from typing import Final, overload

from apeiron.zphi import ZPhi

__all__ = [
    "ICOSAHEDRAL",
    "ImproperRotation",
    "Mat3",
    "ROT_2",
    "ROT_3",
    "ROT_5",
    "Rotation",
    "Vec3",
    "determinant",
    "is_proper",
]


_Z = ZPhi(0, 0)
_ONE = ZPhi(1, 0)
_TWO = ZPhi(2, 0)


def _halve_zphi(z: ZPhi) -> ZPhi:
    """Divide ``z`` by 2, asserting exactness in Z[phi].

    The single place the ×2 storage convention is discharged. Raises
    ``AssertionError`` if either integer component is odd — which would
    mean a caller has produced a value outside the image of the ×2
    storage map, i.e. a bug upstream has violated the I-invariance
    denominator bound of CLAUDE.md §3.3.

    This assertion is a correctness invariant of the pipeline. Do not
    wrap in ``try``/``except``, downgrade to a warning, or add fallback
    behaviour.
    """
    if z.a % 2 != 0 or z.b % 2 != 0:
        raise AssertionError(
            f"Halving {z!r} would leave Z[phi]: both a={z.a} and b={z.b} "
            "must be even. This violates the denominator-2 invariant "
            "(CLAUDE.md §3.3); some upstream operation has produced a "
            "value outside the image of the ×2 storage map."
        )
    return ZPhi(z.a // 2, z.b // 2)


@dataclass(frozen=True, slots=True)
class Vec3:
    """A 3-vector with ZPhi components.

    Used both as a position (stored as 2*v per CLAUDE.md §3.2) and as a
    row of a ``Mat3``. The class is agnostic to the storage convention;
    the ×2 factor is discharged by ``Rotation.apply``, not here.
    """

    x: ZPhi
    y: ZPhi
    z: ZPhi

    def __getitem__(self, i: int) -> ZPhi:
        if i == 0:
            return self.x
        if i == 1:
            return self.y
        if i == 2:
            return self.z
        raise IndexError(f"Vec3 index out of range: {i}")

    def __add__(self, other: Vec3) -> Vec3:
        if not isinstance(other, Vec3):
            return NotImplemented  # type: ignore[return-value]
        return Vec3(self.x + other.x, self.y + other.y, self.z + other.z)

    def __sub__(self, other: Vec3) -> Vec3:
        if not isinstance(other, Vec3):
            return NotImplemented  # type: ignore[return-value]
        return Vec3(self.x - other.x, self.y - other.y, self.z - other.z)

    def __neg__(self) -> Vec3:
        return Vec3(-self.x, -self.y, -self.z)

    def dot(self, other: Vec3) -> ZPhi:
        return self.x * other.x + self.y * other.y + self.z * other.z

    def norm_squared(self) -> ZPhi:
        """Sum of squared components.

        For a position stored as 2*v, this is 4*|v|^2 in real
        coordinates. Rotations preserve this quantity exactly, so
        equality of ``norm_squared`` before and after ``Rotation.apply``
        is the storage-coord form of the norm-preservation invariant.
        """
        return self.dot(self)

    def _halve(self) -> Vec3:
        return Vec3(_halve_zphi(self.x), _halve_zphi(self.y), _halve_zphi(self.z))


@dataclass(frozen=True, slots=True)
class Mat3:
    """A 3×3 matrix with ZPhi entries, stored row-major as three ``Vec3`` rows."""

    row0: Vec3
    row1: Vec3
    row2: Vec3

    @classmethod
    def from_ints(cls, entries: tuple[tuple[int, int], ...]) -> Mat3:
        """Build a Mat3 from a 9-tuple of (a, b) pairs, row-major.

        Each ``(a, b)`` is interpreted as ``ZPhi(a, b) = a + b*phi``.
        """
        if len(entries) != 9:
            raise ValueError(f"Expected 9 (a, b) pairs, got {len(entries)}")
        v = [ZPhi(a, b) for (a, b) in entries]
        return cls(
            Vec3(v[0], v[1], v[2]),
            Vec3(v[3], v[4], v[5]),
            Vec3(v[6], v[7], v[8]),
        )

    @classmethod
    def identity(cls) -> Mat3:
        return cls(Vec3(_ONE, _Z, _Z), Vec3(_Z, _ONE, _Z), Vec3(_Z, _Z, _ONE))

    def row(self, i: int) -> Vec3:
        if i == 0:
            return self.row0
        if i == 1:
            return self.row1
        if i == 2:
            return self.row2
        raise IndexError(f"Mat3 row index out of range: {i}")

    def col(self, j: int) -> Vec3:
        return Vec3(self.row0[j], self.row1[j], self.row2[j])

    def transpose(self) -> Mat3:
        return Mat3(self.col(0), self.col(1), self.col(2))

    @overload
    def __matmul__(self, other: Vec3) -> Vec3: ...
    @overload
    def __matmul__(self, other: Mat3) -> Mat3: ...
    def __matmul__(self, other: Vec3 | Mat3) -> Vec3 | Mat3:
        if isinstance(other, Vec3):
            return Vec3(
                self.row0.dot(other),
                self.row1.dot(other),
                self.row2.dot(other),
            )
        if isinstance(other, Mat3):
            c0, c1, c2 = other.col(0), other.col(1), other.col(2)
            return Mat3(
                Vec3(self.row0.dot(c0), self.row0.dot(c1), self.row0.dot(c2)),
                Vec3(self.row1.dot(c0), self.row1.dot(c1), self.row1.dot(c2)),
                Vec3(self.row2.dot(c0), self.row2.dot(c1), self.row2.dot(c2)),
            )
        return NotImplemented  # type: ignore[return-value]

    def _halve(self) -> Mat3:
        return Mat3(self.row0._halve(), self.row1._halve(), self.row2._halve())


@dataclass(frozen=True, slots=True)
class Rotation:
    """An element of I, stored as 2·g (CLAUDE.md §3.3).

    The ``matrix`` field is the Z[phi] numerator of the rotation under
    a class-level implicit denominator of 2. The denominator is not a
    parameter of the type and is never exposed at the API boundary.
    """

    matrix: Mat3

    @classmethod
    def identity(cls) -> Rotation:
        # Identity in real coords is I_3; stored as 2·I_3.
        return cls(Mat3(Vec3(_TWO, _Z, _Z), Vec3(_Z, _TWO, _Z), Vec3(_Z, _Z, _TWO)))

    def apply(self, v: Vec3) -> Vec3:
        """Apply this rotation to a position vector in ×2 storage form.

        Input ``v`` is ``2*v_real``; output is ``2*(g·v_real)``. The
        product ``self.matrix @ v`` equals ``4*(g·v_real)``, which we
        halve to recover the ×2 form. The halve is exact by CLAUDE.md
        §3.3 and fires loudly on any parity violation.
        """
        return (self.matrix @ v)._halve()

    def compose(
        self, other: Rotation | ImproperRotation,
    ) -> Rotation | ImproperRotation:
        """Compose two isometries: ``self ∘ other``.

        ``Rotation ∘ Rotation = Rotation``: ``self.matrix @ other.matrix``
        equals ``4*(g ∘ h)``; halving returns ``2*(g ∘ h)``, the
        storage form. The denominator is preserved because ``g ∘ h ∈ I``
        and I fits in (1/2)Z[phi] uniformly.

        ``Rotation ∘ Improper = Improper``: ``r ∘ (s ∘ -I) = (r ∘ s) ∘ -I``,
        i.e., ``Improper(r ∘ s)``. Point inversion commutes with every
        rotation.
        """
        if isinstance(other, Rotation):
            return Rotation((self.matrix @ other.matrix)._halve())
        if isinstance(other, ImproperRotation):
            composed_proper = Rotation(
                (self.matrix @ other.rotation.matrix)._halve()
            )
            return ImproperRotation(composed_proper)
        return NotImplemented  # type: ignore[return-value]

    def __matmul__(
        self, other: Rotation | ImproperRotation,
    ) -> Rotation | ImproperRotation:
        if not isinstance(other, (Rotation, ImproperRotation)):
            return NotImplemented  # type: ignore[return-value]
        return self.compose(other)

    def inverse(self) -> Rotation:
        """Inverse of this rotation. For g ∈ SO(3), g^{-1} = g^T.

        In storage form ``2·g``, the inverse storage is ``2·g^T``,
        which is just the transpose of the stored matrix.
        """
        return Rotation(self.matrix.transpose())


@dataclass(frozen=True, slots=True)
class ImproperRotation:
    """An orientation-reversing isometry of R^3: a rotation in I composed
    with point inversion ``v ↦ -v`` (CLAUDE.md §2.1, Claude (web) Q2
    ruling 2026-04-28).

    Members of I_h \\ I (the 60 orientation-reversing elements of the
    full icosahedral group I_h, order 120). Determinant is ``-1``;
    members of ``Rotation`` (= I) have determinant ``+1``. Distinct
    types deliberately: orientation-preserving and orientation-reversing
    isometries are different mathematical objects, and conflating them
    into a single type would force every downstream consumer to check a
    flag to decide which kind it has — at which point the wrapper has
    been reinvented with worse structure.

    The action is ``v ↦ self.rotation.apply(-v)``, equivalently
    ``v ↦ -self.rotation.apply(v)``. Since point inversion ``-I``
    commutes with every ``g ∈ I``, the factorisation ``Improper(r) =
    r ∘ (-I) = (-I) ∘ r`` is unambiguous and ``self.rotation`` is the
    unique proper component.

    Used by Track A's Danzer ABCK substitution: Koca et al.
    (arXiv 2003.13449, eqs. 19–31) and Paolini's POV-Ray source both
    place several of the 25 children using orientation-reversing
    isometries — encoded here as ``ImproperRotation`` instances rather
    than as ``Rotation`` plus a parity flag.
    """

    rotation: Rotation

    def apply(self, v: Vec3) -> Vec3:
        """Apply to a position vector in ×2 storage form.

        ``Improper(r).apply(v) = r.apply(-v) = -r.apply(v)``. Both
        forms agree because point inversion commutes with rotation;
        the negation form is one fewer Mat3 multiplication so it's
        what we use.
        """
        return -self.rotation.apply(v)

    def compose(
        self, other: Rotation | ImproperRotation,
    ) -> Rotation | ImproperRotation:
        """Compose: ``self ∘ other``.

        ``Improper(r) ∘ Rotation(s) = Improper(r ∘ s)``;
        ``Improper(r) ∘ Improper(s) = Rotation(r ∘ s)``
        (the two parity flips cancel).
        """
        if isinstance(other, Rotation):
            return ImproperRotation(self.rotation.compose(other))
        if isinstance(other, ImproperRotation):
            return self.rotation.compose(other.rotation)
        return NotImplemented  # type: ignore[return-value]

    def __matmul__(
        self, other: Rotation | ImproperRotation,
    ) -> Rotation | ImproperRotation:
        if not isinstance(other, (Rotation, ImproperRotation)):
            return NotImplemented  # type: ignore[return-value]
        return self.compose(other)

    def inverse(self) -> ImproperRotation:
        """Inverse: ``Improper(r)^{-1} = Improper(r^{-1})``.

        ``(r ∘ -I)^{-1} = (-I)^{-1} ∘ r^{-1} = -I ∘ r^{-1}``, and
        ``-I`` commutes with ``r^{-1}``, so ``= r^{-1} ∘ -I =
        Improper(r^{-1})``. The proper-component inverse is the
        familiar transpose.
        """
        return ImproperRotation(self.rotation.inverse())


def is_proper(g: Rotation | ImproperRotation) -> bool:
    """Return ``True`` iff ``g`` is orientation-preserving (in I).

    The type-discrimination check used downstream of any code that
    accepts ``Rotation | ImproperRotation``. ``Rotation`` ↦ ``True``;
    ``ImproperRotation`` ↦ ``False``. No float, no determinant
    computation — the type is the answer.
    """
    if isinstance(g, Rotation):
        return True
    if isinstance(g, ImproperRotation):
        return False
    raise TypeError(
        f"is_proper expects Rotation or ImproperRotation; "
        f"got {type(g).__name__}."
    )


def determinant(g: Rotation | ImproperRotation) -> int:
    """Return the determinant of ``g``: ``+1`` for ``Rotation``,
    ``-1`` for ``ImproperRotation``.

    Exact integer-valued; no float fallback. Equivalent to
    ``+1 if is_proper(g) else -1``.
    """
    return 1 if is_proper(g) else -1


# -- generators ---------------------------------------------------------
#
# Three generators of I, each stored as 2·g per CLAUDE.md §3.3.
# BFS closure of any pair produces all 60 elements; the third is
# included for redundancy and for direct order-2 / order-3 / order-5
# tests.

# Order-5 rotation about the axis (0, 1, phi).
#   Real form (derived via Rodrigues):
#     S = (1/2) * [[phi-1, -phi,  1   ],
#                  [phi,    1,    phi-1],
#                  [-1,     phi-1, phi ]]
# Storage form 2·S has entries directly in Z[phi]:
ROT_5: Final[Rotation] = Rotation(Mat3.from_ints((
    (-1, 1),  (0, -1),  (1, 0),
    (0, 1),   (1, 0),   (-1, 1),
    (-1, 0),  (-1, 1),  (0, 1),
)))

# Order-3 rotation about the axis (1, 1, 1) (a dodecahedron vertex /
# icosahedron face-centre direction). Real form is the cyclic
# permutation (x, y, z) -> (y, z, x):
#     T = [[0, 1, 0], [0, 0, 1], [1, 0, 0]]
# Storage form 2·T:
ROT_3: Final[Rotation] = Rotation(Mat3.from_ints((
    (0, 0),  (2, 0),  (0, 0),
    (0, 0),  (0, 0),  (2, 0),
    (2, 0),  (0, 0),  (0, 0),
)))

# Order-2 rotation about the x-axis. The x-axis passes through the
# icosahedron edge midpoints (±phi, 0, 0) — midpoint of the edge
# between (phi, 0, 1) and (phi, 0, -1). Real form:
#     R = diag(1, -1, -1)
# Storage form 2·R = diag(2, -2, -2):
ROT_2: Final[Rotation] = Rotation(Mat3.from_ints((
    (2, 0),   (0, 0),   (0, 0),
    (0, 0),   (-2, 0),  (0, 0),
    (0, 0),   (0, 0),   (-2, 0),
)))


def _close_group(generators: Iterable[Rotation]) -> tuple[Rotation, ...]:
    """BFS closure under composition, starting from the identity.

    Iteratively composes each frontier element with every generator;
    adds unseen products to the queue. Terminates when the queue is
    empty. Every halve-assertion during the closure is a genuine bug,
    not a group-closure artefact — see CLAUDE.md §3.3.
    """
    gens = tuple(generators)
    identity = Rotation.identity()
    seen: set[Rotation] = {identity}
    order: list[Rotation] = [identity]
    head = 0
    while head < len(order):
        g = order[head]
        head += 1
        for h in gens:
            gh = g.compose(h)
            if gh not in seen:
                seen.add(gh)
                order.append(gh)
    return tuple(order)


ICOSAHEDRAL: Final[tuple[Rotation, ...]] = _close_group([ROT_5, ROT_3, ROT_2])
