"""Track B — substitution-first search.

Per CLAUDE.md §6.2 and Claude (web)'s Q7c 2026-04-29 ruling: search
for primitive non-negative integer substitution matrices with target
Perron–Frobenius eigenvalue in ℤ[φ], then attempt geometric
realisation. Independent of Track A's face-merge no-go (commits
`4406900`, `42b3d9a`); the algebraic search has a much larger and
better-defined space than the geometric merge.

Sub-package layout (per Q7c, evolving as Track B builds out):

- ``matrix_search`` — primitive matrix enumerator (this commit).
- ``geometric_prefilter`` — three filters before realisation
  (vertex-class, dihedral-angle, Euler-characteristic). Pending.
- ``realisation`` — vertex-placement CSP given a filtered matrix
  candidate. Pending.

Track B is a sub-package, not a flat module, because it's a research
direction with multiple distinct components (algebraic search,
geometric pre-filter, CSP solver) that each warrant their own file.
The overall flat-module convention from CLAUDE.md §4 still holds for
``apeiron``'s seven core modules; ``track_b`` is a self-contained
exception authorised by Claude (web)'s Q7c ruling.
"""

from apeiron.track_b.geometric_prefilter import (
    PrefilterResult,
    pf_eigenvector_in_zphi,
    prefilter,
    vertex_class_consistent,
)
from apeiron.track_b.matrix_search import enumerate_primitive_matrices
from apeiron.track_b.realisation import (
    ChildPlacement,
    Inconclusive,
    NoRealisation,
    Realised,
    SearchProgress,
    realise,
)

__all__ = [
    "ChildPlacement",
    "Inconclusive",
    "NoRealisation",
    "PrefilterResult",
    "Realised",
    "SearchProgress",
    "enumerate_primitive_matrices",
    "pf_eigenvector_in_zphi",
    "prefilter",
    "realise",
    "vertex_class_consistent",
]
