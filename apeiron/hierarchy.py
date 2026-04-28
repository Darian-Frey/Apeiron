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

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import Protocol, runtime_checkable

from apeiron.corona import CoronaConfig
from apeiron.substitution import (
    PositionedTile,
    SubstitutionRule,
    is_primitive,
    perron_frobenius_in_zphi,
    substitution_matrix,
)
from apeiron.symmetry import Rotation, Vec3
from apeiron.util import pillar
from apeiron.zphi import ZPhi

__all__ = [
    "FourthPillarArgument",
    "HierarchicalCounterexample",
    "HierarchicalWitness",
    "IndistinguishablePair",
    "InflationArgument",
    "InflationFailure",
    "PatchTile",
    "PlacedSubtile",
    "RecognisabilityResult",
    "Supertile",
    "TilePatch",
    "expand_one",
    "expand_supertile",
    "expand_supertile_with_parents",
    "inflation_argument",
    "is_recognisable",
    "make_euclidean_squared_oracle",
    "neighbourhood_signature",
    "patch_from_supertile",
    "shell_neighbourhood_signature",
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


@dataclass(frozen=True, slots=True)
class PlacedSubtile:
    """A leaf-level positioned tile in a level-``N`` supertile expansion.

    Distinct from ``substitution.PositionedTile`` (which describes a
    child relative to its parent in a *single* substitution step) and
    from ``corona.PlacedTile`` (which describes a tile's position in
    a *corona* configuration). ``PlacedSubtile`` is the leaf in a
    flattened level-``N`` expansion: its ``translation`` and
    ``rotation`` are the accumulated isometry from the level-``N``
    parent's frame down to the leaf's local frame.

    Used by ``expand_supertile``; consumed by callers building a
    ``TilePatch`` for ``is_recognisable``.
    """

    prototile_index: int
    translation: Vec3
    rotation: Rotation


def expand_supertile(
    rule: SubstitutionRule,
    prototile_index: int,
    level: int,
) -> tuple[PlacedSubtile, ...]:
    """Recursively flatten ``σ^level(prototile_{prototile_index})``
    into a tuple of ``PlacedSubtile`` leaves.

    Recursion semantics, for a substitution rule with linear
    inflation matrix ``M`` and dissections ``rule.dissections``:

    - ``level == 0`` → one leaf at origin / identity rotation, of
      type ``prototile_index``. The recursion base case.
    - ``level == n`` → for each leaf at ``(T, G)`` of type ``t`` in
      the level-``(n-1)`` expansion, scale ``T`` by ``M`` (the
      single-step inflation, applied once per recursion level) and,
      for every child ``(t_c, g_c)`` of ``rule.dissections[t]``,
      yield a new leaf at ``(M·T + G·t_c, G ∘ g_c)``.

    The position update has two components: ``M·T`` accounts for the
    fact that going from level ``n-1`` to level ``n`` inflates the
    whole previous-level configuration by ``M``; ``G·t_c`` is the
    child's offset from its parent leaf, rotated by the parent's
    accumulated orientation ``G`` so the offset aligns with the
    parent's pose in the level-``n`` frame.

    Returns a tuple sized to ``column_sum(M, prototile_index)^level``
    in the limit (Perron–Frobenius) — for the cube-corona-sized
    test rules this is fine, but for high-level expansions on rules
    with large eigenvalues the count grows exponentially. Cap the
    ``level`` argument accordingly when calling.

    Raises ``ValueError`` if ``prototile_index`` is out of range or
    if ``level < 0``.
    """
    if not 0 <= prototile_index < rule.n_prototiles:
        raise ValueError(
            f"prototile_index {prototile_index} outside "
            f"[0, {rule.n_prototiles})."
        )
    if level < 0:
        raise ValueError(f"level must be ≥ 0; got {level}.")

    origin = Vec3(ZPhi(0, 0), ZPhi(0, 0), ZPhi(0, 0))
    identity = Rotation.identity()

    if level == 0:
        return (
            PlacedSubtile(
                prototile_index=prototile_index,
                translation=origin,
                rotation=identity,
            ),
        )

    leaves: tuple[PlacedSubtile, ...] = expand_supertile(
        rule, prototile_index, level - 1,
    )
    M = rule.inflation
    next_leaves: list[PlacedSubtile] = []
    for leaf in leaves:
        scaled_t = M @ leaf.translation
        for child in rule.dissections[leaf.prototile_index]:
            new_t = scaled_t + leaf.rotation.apply(child.translation)
            new_g = leaf.rotation.compose(child.rotation)
            next_leaves.append(
                PlacedSubtile(
                    prototile_index=child.prototile_index,
                    translation=new_t,
                    rotation=new_g,
                )
            )
    return tuple(next_leaves)


def expand_supertile_with_parents(
    rule: SubstitutionRule,
    prototile_index: int,
    level: int,
) -> tuple[tuple[PlacedSubtile, int], ...]:
    """Like ``expand_supertile``, but each leaf is paired with its
    level-1 parent index — the position ``0..k-1`` of the level-1
    child of ``P_{prototile_index}`` from which this leaf descends.

    The level-1 parent is what pillar-2 recognisability tests as the
    ``parent_supertile`` field of a ``PatchTile``: σⁿ(P) decomposes
    as ``⋃_i σⁿ⁻¹(c_i)`` where ``{c_i}`` is the level-1 dissection
    of ``P``, and the recognisability question is whether each leaf's
    membership in some ``σⁿ⁻¹(c_i)`` is determinable from its bounded
    local neighbourhood. The integer paired with each leaf in the
    output is exactly that membership index ``i``.

    Recursion semantics: at level 1, the ``i``-th child of
    ``rule.dissections[prototile_index]`` is tagged with parent index
    ``i``. At higher levels, the tag is preserved through subsequent
    expansion steps — a level-``n`` leaf inherits its level-``(n-1)``
    predecessor's tag, since the recursion ``σⁿ = σ(σⁿ⁻¹)`` applies
    ``σ`` to each predecessor and the new children inherit the
    predecessor's level-1 ancestor.

    For ``level == 0`` the dissection has not been applied; the single
    leaf is conventionally tagged ``0``. Callers should treat
    ``level == 0`` as a degenerate case for which "parent supertile"
    is not meaningful.

    The leaf component of each pair is identical (in geometry and in
    iteration order) to the corresponding entry in
    ``expand_supertile(rule, prototile_index, level)``; the only
    addition is the parent-index tag.

    Raises ``ValueError`` if ``prototile_index`` is out of range or
    if ``level < 0``.
    """
    if not 0 <= prototile_index < rule.n_prototiles:
        raise ValueError(
            f"prototile_index {prototile_index} outside "
            f"[0, {rule.n_prototiles})."
        )
    if level < 0:
        raise ValueError(f"level must be ≥ 0; got {level}.")

    origin = Vec3(ZPhi(0, 0), ZPhi(0, 0), ZPhi(0, 0))
    identity = Rotation.identity()

    if level == 0:
        return (
            (
                PlacedSubtile(
                    prototile_index=prototile_index,
                    translation=origin,
                    rotation=identity,
                ),
                0,
            ),
        )

    children = rule.dissections[prototile_index]
    leaves: list[tuple[PlacedSubtile, int]] = [
        (
            PlacedSubtile(
                prototile_index=c.prototile_index,
                translation=c.translation,
                rotation=c.rotation,
            ),
            i,
        )
        for i, c in enumerate(children)
    ]
    M = rule.inflation
    for _ in range(2, level + 1):
        next_leaves: list[tuple[PlacedSubtile, int]] = []
        for leaf, parent_id in leaves:
            scaled_t = M @ leaf.translation
            for child in rule.dissections[leaf.prototile_index]:
                new_t = scaled_t + leaf.rotation.apply(child.translation)
                new_g = leaf.rotation.compose(child.rotation)
                next_leaves.append(
                    (
                        PlacedSubtile(
                            prototile_index=child.prototile_index,
                            translation=new_t,
                            rotation=new_g,
                        ),
                        parent_id,
                    )
                )
        leaves = next_leaves
    return tuple(leaves)


def make_euclidean_squared_oracle(
    positions: Sequence[Vec3],
    radius_step_squared: ZPhi,
) -> Callable[[int, int, int], bool]:
    """Build a ``TilePatch.neighbour_within`` oracle from squared
    Euclidean distance.

    ``oracle(i, j, r)`` returns ``True`` iff
    ``|pos_i - pos_j|² ≤ r² · radius_step_squared``.

    Squared distances let the comparison stay in ℤ[φ] (no sqrt, no
    floats). The caller supplies ``radius_step_squared`` — the squared
    distance corresponding to "radius 1" — based on the prototile's
    characteristic edge length. Picking it equal to ``edge_length²``
    (or its ×2-storage scaling) gives the natural "r face-shares away"
    semantics for tiles that share a single edge per radius step.

    Positions are read once and snapshotted into a tuple so that later
    mutations to a caller-side list don't change the oracle's behaviour.
    """
    pos_tuple = tuple(positions)
    n = len(pos_tuple)

    def oracle(i: int, j: int, r: int) -> bool:
        if not 0 <= i < n:
            raise IndexError(f"i={i} outside [0, {n}).")
        if not 0 <= j < n:
            raise IndexError(f"j={j} outside [0, {n}).")
        if r < 0:
            raise ValueError(f"radius must be ≥ 0; got {r}.")
        diff = pos_tuple[i] - pos_tuple[j]
        sq_dist = diff.norm_squared()
        threshold = radius_step_squared * (r * r)
        return sq_dist <= threshold

    return oracle


def patch_from_supertile(
    rule: SubstitutionRule,
    prototile_index: int,
    level: int,
    *,
    radius_step_squared: ZPhi,
) -> TilePatch:
    """Build a ``TilePatch`` from ``σⁿ(P_{prototile_index})`` for
    pillar-2 recognisability testing.

    Each leaf of the level-``n`` expansion becomes a ``PatchTile`` with:

    - ``tile_type``: the leaf's prototile index.
    - ``parent_supertile``: the level-1 ancestor index (from
      ``expand_supertile_with_parents``) — i.e., which of the ``k``
      level-1 children of ``P`` this leaf descends from.

    The neighbour oracle is the squared-Euclidean version produced by
    ``make_euclidean_squared_oracle``. ``radius_step_squared`` is the
    squared distance corresponding to "radius 1" and must be chosen
    by the caller based on the prototile's edge length (in ×2-storage
    units, since leaf positions are ×2-stored per CLAUDE.md §3.2).

    Used in pillar-2 testing as: build a level-``n`` patch, call
    ``is_recognisable`` on it, and check whether the result's
    ``radius_used`` is bounded uniformly across ``n``.

    Raises ``ValueError`` if ``level < 0`` or if the inputs to
    ``expand_supertile_with_parents`` are out of range.
    """
    tagged = expand_supertile_with_parents(rule, prototile_index, level)
    tiles: tuple[PatchTile, ...] = tuple(
        PatchTile(
            tile_type=leaf.prototile_index,
            parent_supertile=parent_id,
        )
        for leaf, parent_id in tagged
    )
    positions: tuple[Vec3, ...] = tuple(
        leaf.translation for leaf, _ in tagged
    )
    oracle = make_euclidean_squared_oracle(positions, radius_step_squared)
    return TilePatch(tiles=tiles, neighbour_within=oracle)


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


def shell_neighbourhood_signature(
    patch: TilePatch,
    tile_index: int,
    radius: int,
) -> tuple[tuple[int, ...], ...]:
    """Refined signature: the per-distance multiset *sequence* up to
    ``radius`` of the centre tile.

    Returns ``(shell_0, shell_1, ..., shell_radius)``, where
    ``shell_k`` is the sorted tuple of tile types at distance exactly
    ``k`` from ``tile_index`` under the patch's neighbour oracle.
    Distance ``k`` is ``"within radius k but not within radius k-1"``;
    shell 0 is the centre tile alone.

    Strictly more discriminating than ``neighbourhood_signature``: the
    multiset signature at radius ``r`` is the disjoint union of
    ``shell_0..shell_r``, while the shell signature preserves the
    distance bucket each tile fell into. The Fibonacci tiling is the
    canonical case where the multiset coalesces tiles that the shell
    sequence separates.

    Caller-supplied as the ``signature_fn`` argument to
    ``is_recognisable`` for tilings where the multiset alone is
    insufficient — see the function's docstring.

    Assumes the oracle is *monotone in radius* — i.e., if tile ``j``
    is within radius ``r``, it is within radius ``r + 1``. The
    Euclidean oracle from ``make_euclidean_squared_oracle`` satisfies
    this; combinatorial / BFS-based oracles must satisfy it for the
    shell decomposition to be well-defined.

    Raises ``IndexError`` for out-of-range ``tile_index`` and
    ``ValueError`` for negative ``radius``.
    """
    n = len(patch.tiles)
    if not 0 <= tile_index < n:
        raise IndexError(
            f"tile_index {tile_index} outside [0, {n})."
        )
    if radius < 0:
        raise ValueError(f"radius must be ≥ 0; got {radius}.")
    shells: list[tuple[int, ...]] = []
    prev_ball: set[int] = set()
    for r in range(radius + 1):
        ball: set[int] = {
            j for j in range(n)
            if patch.neighbour_within(tile_index, j, r)
        }
        shell_indices = ball - prev_ball
        shell_types = sorted(
            patch.tiles[j].tile_type for j in shell_indices
        )
        shells.append(tuple(shell_types))
        prev_ball = ball
    return tuple(shells)


@pillar(2)
def is_recognisable(
    patch: TilePatch,
    *,
    max_radius: int = 5,
    signature_fn: Callable[[TilePatch, int, int], tuple] = (
        neighbourhood_signature
    ),
) -> RecognisabilityResult:
    """Pillar 2: recognisability / border forcing.

    Decide whether the supertile parent of every tile in ``patch``
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

    The ``signature_fn`` keyword swaps the per-radius signature
    function used to group tiles. Default is ``neighbourhood_signature``
    (multiset of types within the radius — sufficient for tilings
    with high local complexity). For tilings where the multiset
    coalesces distinct configurations, pass
    ``shell_neighbourhood_signature`` to refine by per-distance
    shells; the Fibonacci tiling is the canonical example.

    Raises ``ValueError`` if ``max_radius < 1``.
    """
    if max_radius < 1:
        raise ValueError(
            f"max_radius must be ≥ 1; got {max_radius}."
        )
    last_counterexample: IndistinguishablePair | None = None
    for r in range(1, max_radius + 1):
        groups: dict[tuple, list[int]] = {}
        for i in range(len(patch.tiles)):
            sig = signature_fn(patch, i, r)
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
            witness: dict[tuple, int] = {
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


# -- Inflation argument (pillar 3) — sub-commit C ---------------------


@dataclass(frozen=True, slots=True)
class InflationArgument:
    """Pillar-3 witness: aperiodicity by the inflation argument.

    The argument: given a primitive substitution rule with
    Perron–Frobenius eigenvalue ``λ > 1`` and a recognisability
    radius ``r`` for that rule, suppose for contradiction that some
    tiling ``T`` admits a non-zero period vector ``v``. By
    recognisability (pillar 2), every tile's supertile parent in
    ``T`` is determined by its radius-``r`` neighbourhood; so a
    period of ``T`` lifts to a period of the supertile tiling
    ``σ⁻¹(T)`` of the same length ``|v|``. But ``σ⁻¹`` shrinks the
    tiling by factor ``1/λ``, so the same period in supertile
    coordinates has length ``|v|/λ`` in metric units. Iterating
    yields periods of length ``|v|/λⁿ`` for every ``n``, all of
    which are valid periods of the original tiling — contradicting
    finite local complexity (which guarantees a positive lower
    bound on period length). Therefore no period exists.

    The witness records the two pillar inputs that combined to
    establish pillar 3:

    - ``pf_eigenvalue``: the Perron–Frobenius eigenvalue ``λ``,
      strictly greater than 1, recovered exactly from the rule's
      substitution matrix.
    - ``recognisability_radius``: the radius at which the
      recognisability witness was established (pillar 2).

    Per CLAUDE.md §6.3, this argument is "standard" — the
    interesting work is verifying its inputs (pillars 1 + 2) and
    the genuinely-Einstein step (pillar 4, framework only in this
    module).
    """

    pf_eigenvalue: ZPhi
    recognisability_radius: int


@dataclass(frozen=True, slots=True)
class InflationFailure:
    """The inflation argument cannot be invoked: at least one of its
    three preconditions fails.

    ``reason`` identifies which precondition failed:

    - ``"not primitive"`` — the substitution matrix has no
      strictly-positive power; pillar 1 is not satisfied.
    - ``"pf not in z[phi]"`` — the Perron–Frobenius eigenvalue is
      not exactly representable in ℤ[φ]. The rule may still be
      aperiodic, but the inflation argument needs an exact ZPhi
      eigenvalue for the pillar-3 witness; rules outside ℚ(√5)
      need a different eigenvalue ring or the float fallback (not
      yet implemented).
    - ``"pf <= 1"`` — the eigenvalue does not exceed 1, so the
      shrinking step in the inflation argument doesn't apply.
    - ``"not recognisable"`` — pillar 2 was not established at the
      given radius, so the period-lift step has no foundation.

    ``detail`` is a human-readable elaboration tied to the specific
    failure path.
    """

    reason: str
    detail: str


@pillar(3)
def inflation_argument(
    rule: SubstitutionRule,
    recognisability: RecognisabilityResult,
) -> InflationArgument | InflationFailure:
    """Pillar 3: aperiodicity from recognisability.

    Combine pillar 1 (primitivity) and pillar 2 (recognisability)
    to produce a pillar-3 witness, or identify which precondition
    failed.

    Three checks, in order:

    1. The rule's substitution matrix is primitive (some power has
       all-positive entries) — pillar 1.
    2. The Perron–Frobenius eigenvalue is recoverable in ℤ[φ] and
       is strictly greater than 1.
    3. The supplied ``recognisability`` result has
       ``is_recognisable == True`` — pillar 2.

    On all three passing, returns an ``InflationArgument`` carrying
    the eigenvalue and recognisability radius. On any failing,
    returns an ``InflationFailure`` identifying the first failed
    check.

    The structure deliberately separates concerns: the inflation
    argument *itself* is just the implication
    ``primitive ∧ recognisable ∧ λ > 1 ⟹ aperiodic``; the
    interesting machinery is in pillars 1, 2, and (separately)
    pillar 4. This function's job is to verify the implication's
    antecedent and witness it; it does not re-run primitivity or
    recognisability checks beyond the antecedent verification.
    """
    matrix = substitution_matrix(rule)
    if not is_primitive(matrix):
        return InflationFailure(
            reason="not primitive",
            detail=(
                f"Substitution matrix of shape {matrix.shape} has no "
                "strictly-positive power within the Wielandt bound; "
                "pillar 1 is not established."
            ),
        )
    pf = perron_frobenius_in_zphi(matrix)
    if pf is None:
        return InflationFailure(
            reason="pf not in z[phi]",
            detail=(
                "Perron–Frobenius eigenvalue is not exactly "
                "representable in ℤ[φ]; pillar 3 in this module "
                "requires a ZPhi eigenvalue."
            ),
        )
    if not (pf > ZPhi(1, 0)):
        return InflationFailure(
            reason="pf <= 1",
            detail=(
                f"Perron–Frobenius eigenvalue {pf} is not strictly "
                "greater than 1; the inflation step's shrinking "
                "factor is not less than 1."
            ),
        )
    if not recognisability.is_recognisable:
        return InflationFailure(
            reason="not recognisable",
            detail=(
                "Recognisability witness failed at "
                f"radius {recognisability.radius_used}; pillar 2 is "
                "not established, so the period-lift step has no "
                "foundation."
            ),
        )
    return InflationArgument(
        pf_eigenvalue=pf,
        recognisability_radius=recognisability.radius_used,
    )


# -- Fourth pillar (framework only) — sub-commit D --------------------
#
# Pillar 4 — "no non-hierarchical tiling exists" — is the genuinely-
# Einstein step (CLAUDE.md §6.3). Per Claude (web) 2026-04-23, this
# pillar is *not generic*: it requires tile-specific case analysis
# about which local configurations can occur, and the proof is
# non-decidable from the substitution rule alone (the hat / spectre
# proofs each spend ~40 pages on it).
#
# This module ships only the *framework* — the protocol every
# candidate's concrete fourth-pillar implementation must satisfy,
# plus the witness / counterexample dataclasses. Concrete
# implementations of ``FourthPillarArgument`` live per-candidate at
#
#     candidates/<name>/fourth_pillar.py
#
# and are written *only* once a candidate has passed pillars 1–3.
# Premature pillar-4 work on a candidate that fails primitivity,
# recognisability, or the inflation argument is wasted effort.


@dataclass(frozen=True, slots=True)
class HierarchicalWitness:
    """Pillar-4 success witness for a particular ``CoronaConfig``.

    The presence of this object asserts that the configuration arises
    from the substitution hierarchy — i.e., it's the local view of
    some tile in some level-``n`` supertile decomposition for some
    ``n``. The fields document the witness:

    - ``supertile_level``: the smallest ``n`` for which this
      configuration appears as a sub-pattern of the level-``n``
      supertile of some prototile.
    - ``parent_prototile_index``: the prototile whose level-``n``
      supertile contains this configuration.
    - ``rationale``: human-readable explanation of why this
      configuration is hierarchical, supplied by the candidate-
      specific implementation. Per CLAUDE.md §10, the rationale
      should never be hand-waving — it should describe the
      case-analysis branch this configuration falls into.

    The candidate's ``FourthPillarArgument.verify_hierarchical``
    method returns one of these on success. The framework itself
    does not construct ``HierarchicalWitness`` instances; concrete
    candidates do.
    """

    supertile_level: int
    parent_prototile_index: int
    rationale: str


@dataclass(frozen=True, slots=True)
class HierarchicalCounterexample:
    """Pillar-4 failure witness for a particular ``CoronaConfig``.

    A configuration that cannot be embedded in any level-``n``
    supertile decomposition. If a candidate exhibits even one such
    configuration, pillar 4 fails — the tile admits a tiling that
    doesn't arise from the substitution hierarchy, and the
    aperiodicity argument doesn't go through.

    - ``rationale``: human-readable explanation of why no embedding
      exists. Typically describes which face/edge/vertex constraint
      the configuration violates relative to every prototile's
      level-``n`` decomposition for ``n`` up to the search bound.
    - ``search_bound``: the largest level ``n`` checked before
      concluding non-embedding. ``HierarchicalCounterexample`` is
      *bounded-search* evidence; a configuration that survives every
      level up to ``search_bound`` is treated as a counterexample,
      but in principle a sufficiently-deep level might admit it.
      Concrete implementations should pick ``search_bound`` large
      enough that the conclusion is robust for the candidate at
      hand.
    """

    rationale: str
    search_bound: int


@pillar(4)
@runtime_checkable
class FourthPillarArgument(Protocol):
    """Pillar 4: no non-hierarchical tiling exists.

    Protocol for the candidate-specific pillar-4 argument.

    A concrete implementation lives at
    ``candidates/<name>/fourth_pillar.py`` and supplies two
    methods:

    - ``local_configurations() -> frozenset[CoronaConfig]``: every
      local configuration the candidate's tilings can exhibit. For
      a candidate that has passed pillars 1–3 (primitive substitution
      with recognisability radius ``r``), the relevant local
      configurations are typically the ``corona_r``-style
      neighbourhoods of every prototile under every rotation of the
      icosahedral group I. The concrete implementation is
      responsible for enumerating *all* of them; pillar 4 is the
      claim that the hierarchical embedding holds for every one.
    - ``verify_hierarchical(config: CoronaConfig) -> HierarchicalWitness
      | HierarchicalCounterexample``: for a given configuration,
      either return a hierarchical witness (success) or a
      counterexample (failure). The candidate establishes pillar 4
      iff every configuration in ``local_configurations()`` returns
      a ``HierarchicalWitness``.

    Per Claude (web), this protocol is the entire framework
    contribution of ``hierarchy.py`` to pillar 4. There is no
    generic algorithmic check; the work happens per-candidate in
    ``candidates/<name>/fourth_pillar.py``. The presence of a
    ``FourthPillarArgument`` implementation in that directory is
    itself the entry point to pillar-4 verification for the
    candidate.
    """

    def local_configurations(self) -> frozenset[CoronaConfig]:
        """Every local configuration the candidate's tilings exhibit."""
        ...

    def verify_hierarchical(
        self, config: CoronaConfig,
    ) -> HierarchicalWitness | HierarchicalCounterexample:
        """Decide hierarchical embedding for one configuration."""
        ...
