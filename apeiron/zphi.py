"""Exact arithmetic over Z[phi] = Z + Z*phi, phi = (1 + sqrt(5)) / 2.

The canonical element is a frozen ``ZPhi`` dataclass holding integer
coefficients ``(a, b)`` representing ``a + b*phi``. Addition, subtraction,
negation, multiplication, integer powers, Galois conjugation, field norm,
field trace, equality, hashing, and total ordering are all exact integer
operations using the relation ``phi**2 = phi + 1``.

This module has no float dependency. Any function in this module that
accepts or returns a float is a bug (CLAUDE.md §7.3).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final

__all__ = ["ZPhi"]


def _require_int(name: str, value: object) -> int:
    """Reject anything that is not a strict Python ``int``.

    ``bool`` is a subclass of ``int`` in Python; we exclude it explicitly so
    that ``ZPhi(True, False)`` raises rather than silently becoming
    ``ZPhi(1, 0)``. Rejecting floats here is the first line of defence
    against float contamination of the verification pipeline.
    """
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(
            f"ZPhi.{name} must be a Python int; got {type(value).__name__}"
        )
    return value


@dataclass(frozen=True, slots=True, order=False)
class ZPhi:
    """Element ``a + b*phi`` of Z[phi], with ``a, b`` in Z.

    Parameters
    ----------
    a : int
        The rational-integer component.
    b : int
        The phi component.

    Notes
    -----
    The pair ``(a, b)`` is itself canonical: ``{1, phi}`` is a free
    Z-basis of Z[phi], so there is no reduction step. Equality and
    hashing are induced by the dataclass.
    """

    a: int
    b: int

    def __post_init__(self) -> None:
        _require_int("a", self.a)
        _require_int("b", self.b)

    # ---- constructors --------------------------------------------------

    @classmethod
    def from_int(cls, n: int) -> ZPhi:
        """Embed ``n`` in Z[phi] as ``n + 0*phi``."""
        _require_int("n", n)
        return cls(n, 0)

    # ---- arithmetic ----------------------------------------------------

    def __add__(self, other: ZPhi | int) -> ZPhi:
        if isinstance(other, ZPhi):
            return ZPhi(self.a + other.a, self.b + other.b)
        if isinstance(other, int) and not isinstance(other, bool):
            return ZPhi(self.a + other, self.b)
        return NotImplemented  # type: ignore[return-value]

    def __radd__(self, other: int) -> ZPhi:
        return self.__add__(other)

    def __sub__(self, other: ZPhi | int) -> ZPhi:
        if isinstance(other, ZPhi):
            return ZPhi(self.a - other.a, self.b - other.b)
        if isinstance(other, int) and not isinstance(other, bool):
            return ZPhi(self.a - other, self.b)
        return NotImplemented  # type: ignore[return-value]

    def __rsub__(self, other: int) -> ZPhi:
        if isinstance(other, int) and not isinstance(other, bool):
            return ZPhi(other - self.a, -self.b)
        return NotImplemented  # type: ignore[return-value]

    def __neg__(self) -> ZPhi:
        return ZPhi(-self.a, -self.b)

    def __pos__(self) -> ZPhi:
        return self

    def __mul__(self, other: ZPhi | int) -> ZPhi:
        # (a + b*phi)(c + d*phi)
        #   = a*c + (a*d + b*c)*phi + b*d*phi**2
        #   = a*c + (a*d + b*c)*phi + b*d*(phi + 1)
        #   = (a*c + b*d) + (a*d + b*c + b*d)*phi
        if isinstance(other, ZPhi):
            a, b = self.a, self.b
            c, d = other.a, other.b
            return ZPhi(a * c + b * d, a * d + b * c + b * d)
        if isinstance(other, int) and not isinstance(other, bool):
            return ZPhi(self.a * other, self.b * other)
        return NotImplemented  # type: ignore[return-value]

    def __rmul__(self, other: int) -> ZPhi:
        return self.__mul__(other)

    def __pow__(self, exponent: int) -> ZPhi:
        _require_int("exponent", exponent)
        if exponent < 0:
            raise ValueError(
                "ZPhi supports only non-negative integer powers; negative "
                "exponents would leave Z[phi] in general."
            )
        result = ZPhi(1, 0)
        base = self
        e = exponent
        while e > 0:
            if e & 1:
                result = result * base
            e >>= 1
            if e > 0:
                base = base * base
        return result

    # ---- number-field operations --------------------------------------

    def conjugate(self) -> ZPhi:
        """Galois conjugate: send ``phi`` to ``1 - phi``.

        ``conj(a + b*phi) = a + b*(1 - phi) = (a + b) - b*phi``.
        """
        return ZPhi(self.a + self.b, -self.b)

    def norm(self) -> int:
        """Field norm N(a + b*phi) = a**2 + a*b - b**2.

        This is ``(a + b*phi) * conj(a + b*phi)`` as an element of Z.
        """
        return self.a * self.a + self.a * self.b - self.b * self.b

    def trace(self) -> int:
        """Field trace ``tr(a + b*phi) = 2*a + b``."""
        return 2 * self.a + self.b

    # ---- ordering -----------------------------------------------------

    def _sign(self) -> int:
        """Return the sign of ``self`` as an integer in {-1, 0, +1}.

        Uses only integer arithmetic. The real embedding of ``a + b*phi``
        is ``(2*a + b + b*sqrt(5)) / 2``. With ``s = 2*a + b``, the sign
        of ``s + b*sqrt(5)`` is decided by comparing ``s**2`` to ``5 * b**2``
        when ``s`` and ``b`` have opposite signs; otherwise it is the sign
        of ``s`` (when ``b == 0``) or of ``b`` (when ``s == 0``).

        Since ``sqrt(5)`` is irrational, ``s + b*sqrt(5) == 0`` forces
        ``s == b == 0``, i.e., ``a == b == 0``.
        """
        a, b = self.a, self.b
        if a == 0 and b == 0:
            return 0
        s = 2 * a + b
        if b == 0:
            return 1 if s > 0 else -1
        if s == 0:
            return 1 if b > 0 else -1
        # Both s and b are non-zero.
        if (s > 0) == (b > 0):
            # Same sign: s + b*sqrt(5) shares that sign.
            return 1 if s > 0 else -1
        # Opposite signs: compare |s| to |b|*sqrt(5), i.e., s**2 vs 5*b**2.
        lhs = s * s
        rhs = 5 * b * b
        # lhs == rhs is impossible for integers (sqrt(5) irrational).
        if lhs > rhs:
            return 1 if s > 0 else -1
        return 1 if b > 0 else -1

    def __lt__(self, other: ZPhi | int) -> bool:
        return self._compare(other) < 0

    def __le__(self, other: ZPhi | int) -> bool:
        return self._compare(other) <= 0

    def __gt__(self, other: ZPhi | int) -> bool:
        return self._compare(other) > 0

    def __ge__(self, other: ZPhi | int) -> bool:
        return self._compare(other) >= 0

    def _compare(self, other: ZPhi | int) -> int:
        if isinstance(other, int) and not isinstance(other, bool):
            other = ZPhi(other, 0)
        if not isinstance(other, ZPhi):
            raise TypeError(
                f"ZPhi ordering requires ZPhi or int; got {type(other).__name__}"
            )
        return (self - other)._sign()

    # ---- display ------------------------------------------------------

    def __repr__(self) -> str:
        return f"ZPhi({self.a}, {self.b})"

    def __str__(self) -> str:
        if self.b == 0:
            return str(self.a)
        if self.a == 0:
            if self.b == 1:
                return "phi"
            if self.b == -1:
                return "-phi"
            return f"{self.b}*phi"
        sign = "+" if self.b > 0 else "-"
        mag = abs(self.b)
        tail = "phi" if mag == 1 else f"{mag}*phi"
        return f"{self.a} {sign} {tail}"


ZERO: Final[ZPhi] = ZPhi(0, 0)
ONE: Final[ZPhi] = ZPhi(1, 0)
PHI: Final[ZPhi] = ZPhi(0, 1)
