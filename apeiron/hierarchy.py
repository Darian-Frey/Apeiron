"""Supertile construction, recognisability, and the four-pillar proof
apparatus (CLAUDE.md §6.3).

A *supertile* of order ``n`` is the inflated prototile ``σⁿ(P)``,
dissected into base prototiles. The dissection is recursive: a level-
``n`` supertile is the level-1 dissection of ``P`` with each child
itself a level-``n-1`` supertile (recursively, down to level 0 which
is just a positioned prototile copy). Lazy expansion is the canonical
form — we never materialise a flat level-``n`` decomposition because
the count grows like ``λⁿ`` (where ``λ = φ²`` for icosahedral systems)
and the level-1 children are the only structure ever reused.

Four pillars (per Claude (web) 2026-04-23 design relay):

1. **Substitution exists and is primitive.** Already covered by
   ``substitution.py`` (``is_primitive`` + ``perron_frobenius_in_zphi``).
   This module re-tags the relevant tests with ``@pillar(1)`` (added
   in sub-commit E) so the four-pillar coverage is greppable.
2. **Recognisability (border forcing).** Every tile's supertile
   membership is determined by a bounded local neighbourhood. The hard
   pillar; lands in sub-commit B.
3. **Aperiodicity from recognisability.** Standard inflation argument
   that any period vector ``v`` would imply ``σⁿ(v)`` is also a period,
   contradicting finite local complexity. Lands in sub-commit C.
4. **No non-hierarchical tiling.** The genuinely-Einstein step. Per
   Claude (web), this is *framework-only* in this module — concrete
   implementations of the ``FourthPillarArgument`` protocol live per-
   candidate in ``candidates/<name>/fourth_pillar.py``, since the
   case-analysis is candidate-specific (40 pages in the hat / spectre
   proofs). Lands in sub-commit D.

Sub-commit E adds the ``@pillar(n)`` decorator (in ``util.py``) and a
sweep that tags every relevant function in this module and the
relevant tests in ``test_substitution.py`` and elsewhere.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from apeiron.substitution import PositionedTile, SubstitutionRule

__all__ = [
    "Supertile",
    "expand_one",
]


@dataclass(frozen=True, slots=True)
class Supertile:
    """A level-``level`` supertile of prototile ``prototile_index``.

    Recursive dataclass. ``children`` is the *level-1* dissection of
    ``σ(P_prototile_index)`` — a tuple of ``PositionedTile`` instances,
    each carrying a (sub-)prototile index, a translation, and a
    rotation. To obtain the level-``n`` decomposition of this
    supertile, recursively apply ``σ`` to each child via
    ``expand_one``; do not pre-materialise flat level-``n`` lists,
    because the count grows like ``λⁿ`` (CLAUDE.md §6.3, Claude (web)
    2026-04-23 directive).

    A ``level == 0`` supertile is a single prototile copy with no
    children — the recursion base case.

    Pillar 1 (substitution exists and is primitive) is established
    upstream by ``substitution.is_primitive`` /
    ``substitution.perron_frobenius_in_zphi`` on the rule's
    substitution matrix; ``Supertile`` itself is purely structural.
    """

    level: int
    prototile_index: int
    children: tuple[PositionedTile, ...] = field(default_factory=tuple)

    def __post_init__(self) -> None:
        if self.level < 0:
            raise ValueError(
                f"Supertile.level must be ≥ 0; got {self.level}."
            )
        if self.prototile_index < 0:
            raise ValueError(
                f"Supertile.prototile_index must be ≥ 0; "
                f"got {self.prototile_index}."
            )
        if self.level == 0 and self.children:
            raise ValueError(
                f"Supertile.level == 0 must have empty children; "
                f"got {len(self.children)} children."
            )


def expand_one(
    rule: SubstitutionRule,
    prototile_index: int,
) -> tuple[PositionedTile, ...]:
    """Return the level-1 dissection of ``σ(prototile_index)``.

    Looks up the dissection from ``rule.dissections``. This is the
    one operation worth memoising (it is invariant per
    ``(rule, prototile_index)``); higher-level expansion is the
    repeated application of ``expand_one`` in a recursion or a
    breadth-first walk, never a flat materialisation.

    Raises ``ValueError`` if ``prototile_index`` is outside
    ``[0, rule.n_prototiles)``.
    """
    if not 0 <= prototile_index < rule.n_prototiles:
        raise ValueError(
            f"prototile_index {prototile_index} outside "
            f"[0, {rule.n_prototiles})."
        )
    return rule.dissections[prototile_index]
