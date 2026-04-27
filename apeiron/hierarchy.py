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

from collections.abc import Callable, Mapping
from dataclasses import dataclass, field

from apeiron.substitution import PositionedTile, SubstitutionRule

__all__ = [
    "IndistinguishablePair",
    "PatchTile",
    "RecognisabilityResult",
    "Supertile",
    "TilePatch",
    "expand_one",
    "is_recognisable",
    "neighbourhood_signature",
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


# -- Recognisability (pillar 2) — sub-commit B ------------------------


@dataclass(frozen=True, slots=True)
class IndistinguishablePair:
    """A counter-example to recognisability at a given radius.

    Two tiles whose radius-``r`` neighbourhood signatures agree but
    whose supertile parents differ — i.e., the local pattern at this
    radius is ambiguous about which supertile each tile belongs to.
    The actionable diagnostic for a researcher: either the radius
    needs to grow, or the tile definition admits genuinely
    non-recognisable configurations.

    Returned in the ``witness`` field of a failing
    ``RecognisabilityResult``.
    """

    tile_a: int
    tile_b: int
    radius: int
    parent_a: int
    parent_b: int


@dataclass(frozen=True, slots=True)
class RecognisabilityResult:
    """Outcome of an ``is_recognisable`` query.

    ``is_recognisable``: ``True`` iff every tile's radius-``radius_used``
    neighbourhood uniquely determines its supertile parent.

    ``radius_used``: the radius at which the conclusion was reached. On
    success, the smallest sufficient radius. On failure, the
    ``max_radius`` cap.

    ``radius_cap_reached``: ``True`` iff the search exhausted
    ``max_radius`` without conclusion. Distinguishes "provably
    non-recognisable at radius ``radius_used``" (in which case the
    failure is a witnessable property at that radius) from "unknown
    beyond cap, may or may not be recognisable at higher radii."

    ``witness``: on success, the map ``signature → parent_supertile``
    that establishes the recognisability function. On failure, an
    ``IndistinguishablePair`` exhibiting the ambiguity.
    """

    is_recognisable: bool
    radius_used: int
    radius_cap_reached: bool
    witness: Mapping[tuple, int] | IndistinguishablePair


@dataclass(frozen=True, slots=True)
class PatchTile:
    """A single tile in a finite patch of a tiling.

    ``tile_type``: the prototile-type identifier (e.g., 0 = square,
    1 = rhombus for Ammann–Beenker; arbitrary integers).
    ``parent_supertile``: which supertile contains this tile in the
    next-level decomposition. Recognisability tests whether this
    parent is determinable from the tile's bounded neighbourhood
    alone.
    """

    tile_type: int
    parent_supertile: int


@dataclass(frozen=True, slots=True)
class TilePatch:
    """A finite patch of tiles with an adjacency-radius oracle.

    ``tiles`` is the ordered tuple; tile indices are positions in
    this tuple. ``neighbour_within(i, j, radius)`` returns ``True``
    iff tiles ``i`` and ``j`` are within graph-distance ``radius`` —
    typically a step-count along a corona-adjacency graph, but the
    abstraction is agnostic to whether the underlying notion of
    radius is geometric (Euclidean), combinatorial (graph), or
    something else. Decoupling here lets ``is_recognisable`` work on
    Ammann–Beenker (silver-ratio inflation, not in ℤ[φ]) without
    having to embed silver-ratio coordinates in the project's
    golden-ratio number system.
    """

    tiles: tuple[PatchTile, ...]
    neighbour_within: Callable[[int, int, int], bool]


def neighbourhood_signature(
    patch: TilePatch,
    tile_index: int,
    radius: int,
) -> tuple[int, ...]:
    """Hashable summary of the tile-type multiset within ``radius``
    of ``tile_index`` (including the centre tile itself).

    Sorted-tuple of tile types makes the result canonical regardless
    of patch ordering. The signature is the equivalence class used
    by ``is_recognisable`` to group tiles for parent-uniqueness
    checking.

    A more refined signature (distinguishing neighbours by relative
    position, rotation, etc.) is reasonable for tiles where the
    multiset alone is insufficient; that's an extension, not the
    base case. For the Ammann–Beenker oracle the multiset signature
    suffices (Baake & Grimm Vol. 1 §6).
    """
    n = len(patch.tiles)
    if not 0 <= tile_index < n:
        raise IndexError(
            f"tile_index {tile_index} outside [0, {n})."
        )
    types: list[int] = [patch.tiles[tile_index].tile_type]
    for j in range(n):
        if j == tile_index:
            continue
        if patch.neighbour_within(tile_index, j, radius):
            types.append(patch.tiles[j].tile_type)
    return tuple(sorted(types))


def is_recognisable(
    patch: TilePatch,
    *,
    max_radius: int = 5,
) -> RecognisabilityResult:
    """Decide whether the supertile parent of every tile in ``patch``
    is uniquely determined by some bounded radius of local
    neighbourhood.

    Iterative search: try ``radius = 1, 2, ..., max_radius``. At each
    radius, group tiles by ``neighbourhood_signature``; if every
    group's tiles share a single ``parent_supertile`` value, return
    success with the witness map. If any group has tiles with
    distinct parents at this radius, record an ``IndistinguishablePair``
    and continue to the next radius. If the loop exhausts
    ``max_radius`` without success, return failure with the latest
    counter-example and ``radius_cap_reached=True``.

    Per Claude (web) 2026-04-23: the predicate cannot decide
    "no radius works" in general (it's undecidable from the
    substitution rule alone for some classes of rules). What it
    *can* decide is "radius ``r`` suffices?", and the failure mode
    surfaces an actionable distinguishing pair the researcher uses
    to either tighten the tile definition or widen the search.
    Default cap of 5 is practical; theoretical bounds are
    impractically loose (~λ^(2n) for an ``n``-letter alphabet with
    inflation factor ``λ``).

    Raises ``ValueError`` if ``max_radius < 1``.
    """
    if max_radius < 1:
        raise ValueError(
            f"max_radius must be ≥ 1; got {max_radius}."
        )
    last_counterexample: IndistinguishablePair | None = None
    for r in range(1, max_radius + 1):
        groups: dict[tuple[int, ...], list[int]] = {}
        for i in range(len(patch.tiles)):
            sig = neighbourhood_signature(patch, i, r)
            groups.setdefault(sig, []).append(i)
        counter: IndistinguishablePair | None = None
        for indices in groups.values():
            if len(indices) < 2:
                continue
            head_parent = patch.tiles[indices[0]].parent_supertile
            for j in indices[1:]:
                if patch.tiles[j].parent_supertile != head_parent:
                    counter = IndistinguishablePair(
                        tile_a=indices[0],
                        tile_b=j,
                        radius=r,
                        parent_a=head_parent,
                        parent_b=patch.tiles[j].parent_supertile,
                    )
                    break
            if counter is not None:
                break
        if counter is None:
            witness: dict[tuple[int, ...], int] = {
                sig: patch.tiles[indices[0]].parent_supertile
                for sig, indices in groups.items()
            }
            return RecognisabilityResult(
                is_recognisable=True,
                radius_used=r,
                radius_cap_reached=False,
                witness=witness,
            )
        last_counterexample = counter
    assert last_counterexample is not None
    return RecognisabilityResult(
        is_recognisable=False,
        radius_used=max_radius,
        radius_cap_reached=True,
        witness=last_counterexample,
    )
