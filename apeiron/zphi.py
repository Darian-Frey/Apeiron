"""Exact arithmetic over Z[phi] = Z + Z*phi, phi = (1 + sqrt(5)) / 2.

The canonical element is a frozen ``ZPhi`` dataclass holding integer
coefficients ``(a, b)`` representing ``a + b*phi``. Addition, multiplication,
negation, equality, and hashing are all exact integer operations using the
relation ``phi^2 = phi + 1``.

No floats appear in this module's public API. See CLAUDE.md §3 and §7.3.
"""
