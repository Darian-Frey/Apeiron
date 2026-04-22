"""Corona BFS engine — the verification core (CLAUDE.md §5.1).

A *corona* around a central tile ``P`` is the Moore neighbourhood of
``P`` under a face-to-face tiling: every other tile that shares a face,
an edge, or a vertex with ``P``, positioned so no interior-overlap
occurs and the full solid angle around every edge and vertex of ``P``
is covered (zero angular defect).

The acceptance criteria for this module (CLAUDE.md §5.1 item 5) are:

- ``CoronaConfig`` dataclass: central tile plus a set of positioned
  neighbour tiles.
- ``corona_1(P)``: enumerate all complete first-shell configurations
  modulo I.
- ``corona_2(config)``: extend to second shell.
- Canonical hashing of configurations for BFS deduplication.
- Tests: cube has a unique ``corona_1``; rhombic dodecahedron has a
  unique ``corona_1``.

**Scope of this commit.** This is sub-commit A of a four-part
breakdown:

- **A (this commit)** — data model: ``PlacedTile``, ``CoronaConfig``,
  canonical ordering of neighbours, equality & hashing.
- **B** — geometric primitives: face-to-face placement enumeration,
  interior-overlap predicate, exact-ℤ[φ] ``angular_defect``.
- **C** — ``corona_1(P)`` BFS engine with constraint propagation;
  cube acceptance test.
- **D** — ``corona_2(config)`` and the rhombic-dodecahedron
  acceptance test.

Splitting follows the pattern used for ``polyhedron.py`` (commits A–E);
the rationale is that each sub-commit ships a reviewable unit of ~200
lines with its own tests, rather than a single ~800-line drop.

**Face-to-face vs complete Moore neighbourhood.** Internally the BFS
(sub-commit C) works face-to-face with constraint propagation —
placing a neighbour can force additional placements to fill the edge-
and vertex-wedges where that neighbour now participates, and any
placement creating an interior overlap or an unfillable wedge is
pruned. Externally ``corona_1`` returns the complete closed
neighbourhood: for a cube, all 26 neighbours (6 face + 12 edge + 8
vertex) in a 3×3×3 block minus the centre. This internal/external
distinction is the one thing most likely to confuse maintainers a
year out — flagging it here so it's captured at the module level.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass

from apeiron.polyhedron import Polyhedron
from apeiron.symmetry import Rotation, Vec3

__all__ = [
    "CoronaConfig",
    "PlacedTile",
]


@dataclass(frozen=True, slots=True)
class PlacedTile:
    """A placed copy of the monotile in a corona: translation + rotation.

    ``translation`` is in the ×2 storage form of CLAUDE.md §3.2.
    ``rotation`` is an element of the icosahedral group I, stored as
    ``2·g`` under the denominator-2 convention of CLAUDE.md §3.3.

    Semantically this is the corona analogue of
    ``apeiron.substitution.PositionedTile`` — both carry a translation
    and a rotation — but the types intentionally differ: ``PlacedTile``
    belongs to a corona and will eventually grow adjacency metadata
    (which face / edges / vertices of the central it participates in;
    which wedge of its own boundary is filled by which other neighbour)
    that does not belong on ``PositionedTile``. Per Claude (web)'s
    2026-04-20 guidance the fields are duplicated rather than factored
    through a shared ``Placement`` type; if a third corona-like
    consumer arises later, refactor then.

    Equality and hash are induced by the frozen dataclass.
    """

    translation: Vec3
    rotation: Rotation


def _placed_tile_key(pt: PlacedTile) -> tuple[int, ...]:
    """Total order on ``PlacedTile`` for deterministic neighbour ordering.

    Sorts first by translation (lex on raw ZPhi coefficient pairs
    component-by-component), then by rotation matrix entries in row-
    major order under the same raw-coefficient lex convention. Not the
    real-value ordering — this is a deterministic surrogate, which is
    all ``CoronaConfig``'s canonical-order invariant requires.
    """
    t = pt.translation
    m = pt.rotation.matrix
    return (
        t.x.a, t.x.b, t.y.a, t.y.b, t.z.a, t.z.b,
        m.row0.x.a, m.row0.x.b, m.row0.y.a, m.row0.y.b, m.row0.z.a, m.row0.z.b,
        m.row1.x.a, m.row1.x.b, m.row1.y.a, m.row1.y.b, m.row1.z.a, m.row1.z.b,
        m.row2.x.a, m.row2.x.b, m.row2.y.a, m.row2.y.b, m.row2.z.a, m.row2.z.b,
    )


@dataclass(frozen=True, slots=True)
class CoronaConfig:
    """A corona configuration: central tile and its placed neighbours.

    Fields
    ------
    central : Polyhedron
        The tile being coronated. All neighbours are placed copies of
        ``central`` under icosahedral rotations and ×2-stored
        translations. (The monotile case is the research target;
        multi-prototile support would require extending ``PlacedTile``
        with a prototile index, which has not been done yet.)
    neighbours : tuple of PlacedTile
        Placed copies of the central surrounding it, in canonical
        order (sorted by ``_placed_tile_key``). Pairwise distinct.

    Invariants (enforced by ``__post_init__``)
    -----------------------------------------
    - ``neighbours`` is sorted ascending by ``_placed_tile_key``.
    - All neighbours are pairwise distinct.

    Construction
    ------------
    Prefer ``CoronaConfig.from_neighbours``, which sorts and
    deduplicates the input iterable. Direct construction is
    permissible when the caller guarantees canonical form (e.g., when
    producing a ``CoronaConfig`` from another, already-canonical one
    via a pure-rotation action).

    Equality and hashing
    --------------------
    Dataclass equality is pointwise: ``c1 == c2`` iff they share the
    same central and the same neighbour tuple. This is *not*
    equivalence up to the action of I on the whole configuration —
    that is a separate orbit-level notion that will be added when the
    BFS engine (sub-commit C) needs it for deduplication.
    """

    central: Polyhedron
    neighbours: tuple[PlacedTile, ...]

    def __post_init__(self) -> None:
        keys = [_placed_tile_key(n) for n in self.neighbours]
        if keys != sorted(keys):
            raise ValueError(
                "CoronaConfig.neighbours must be in canonical order; "
                "use CoronaConfig.from_neighbours to construct from "
                "an unordered iterable."
            )
        if len(set(self.neighbours)) != len(self.neighbours):
            raise ValueError("CoronaConfig has duplicate neighbours.")

    @classmethod
    def from_neighbours(
        cls,
        central: Polyhedron,
        neighbours: Iterable[PlacedTile],
    ) -> CoronaConfig:
        """Construct a canonical ``CoronaConfig`` from any iterable.

        Deduplicates by value, sorts by ``_placed_tile_key``. Duplicate
        neighbours in the input are silently dropped — BFS will
        routinely produce the same placement via multiple search
        paths, and raising on duplicates would push dedup responsibility
        to every caller.
        """
        ordered = sorted(set(neighbours), key=_placed_tile_key)
        return cls(central=central, neighbours=tuple(ordered))
