# Danzer ABCK Substitution Dissections

**Status: template — to be filled by 27B-β.**

This file is the working artifact for sub-commit 27B-β: the
real geometric dissection of σ(X) for X ∈ {A, B, C, K} in
Frettlöh's ABCK substitution. Sub-commit 27B-α landed the matrix-
only `SubstitutionRule` with placeholder dissection geometry;
27B-β replaces those placeholders by transcribing each child tile's
translation and rotation from Frettlöh Figure 2 / the Tilings
Encyclopedia interactive view at
[tilings.math.uni-bielefeld.de](https://tilings.math.uni-bielefeld.de).

## Source priority

1. **Frettlöh, Bielefeld internal report**, Figure 2 — primary.
2. **Tilings Encyclopedia interactive view** — interactive
   cross-check; useful for confirming sign conventions and
   rotation orientations that may be ambiguous from the static
   figure.
3. **Senechal, *Quasicrystals and Geometry*, Chapter 7** —
   secondary cross-check on individual tile orientations.

## Form to record

For each prototile X ∈ {A, B, C, K}, σ(X) decomposes into a
specific multiset of children. The matrix is fixed (Frettlöh,
recovered in 27B-α):

```text
σ(A) → 0 A + 3 B + 2 C + 6 K   (11 children total)
σ(B) → 0 A + 2 B + 1 C + 4 K   (7  children total)
σ(C) → 1 A + 0 B + 2 C + 2 K   (5  children total)
σ(K) → 0 A + 1 B + 0 C + 1 K   (2  children total)
```

For each child of σ(X), record:

| field | type | meaning |
|---|---|---|
| `prototile_index` | int 0–3 | Which of A=0, B=1, C=2, K=3. |
| `translation` | `[[a, b], [a, b], [a, b]]` | Child's origin position in the σ(X)-frame, ×2-stored ZPhi. |
| `rotation` | `Rotation` matrix from `ICOSAHEDRAL` | Child's orientation; one of the 60 elements of I, identified by some canonical label or matrix transcription. |

Convention: child placements are expressed in the σ(X) frame
*before* applying the inflation. After populating, the loader will
combine these with the inflation matrix φ·I from the `SubstitutionRule`.

## σ(A) — 11 children

```text
TODO — transcribe from Frettlöh Figure 2 panel for σ(A).

Children, in any order:
- [ ]  child 1: prototile_index=?, translation=[[?,?], [?,?], [?,?]], rotation=?
- [ ]  child 2: ...
- [ ]  ...
- [ ] child 11: ...
```

## σ(B) — 7 children

```text
TODO — transcribe from Frettlöh Figure 2 panel for σ(B).

- [ ]  child 1: ...
- [ ]  ...
- [ ]  child 7: ...
```

## σ(C) — 5 children

```text
TODO — transcribe from Frettlöh Figure 2 panel for σ(C).

- [ ]  child 1: ...
- [ ]  ...
- [ ]  child 5: ...
```

## σ(K) — 2 children

```text
TODO — transcribe from Frettlöh Figure 2 panel for σ(K).

- [ ]  child 1: ...
- [ ]  child 2: ...
```

## Cross-checks once populated

- Per-tile child count sums to the column sum of the substitution
  matrix (already verified by `test_matrix_column_sums_match_child_counts`).
- After encoding, every child's rotation must be in
  `apeiron.symmetry.ICOSAHEDRAL` (use `find_rotation` from
  `apeiron.corona`); raise loudly if any rotation isn't.
- Children of σ(X) tile φ·X without overlap and without gaps; the
  validation pass (`has_interior_overlap`-style) on the assembled
  level-1 patch should confirm.
- After 27B-β lands, sub-commit 27D runs the full pipeline against
  the real geometry: pillar 2 (recognisability) on a level-N
  Danzer patch, then pillar 3 (inflation argument), then the
  pillar-4 stub.

## Outstanding decisions (for the relay session)

- **Rotation labelling.** The 60 elements of I are stored as a
  `tuple[Rotation, ...]` indexed 0–59 in `apeiron.symmetry`. Should
  this file record a child's rotation as the index (0–59), as the
  matrix entries directly, or as a generator-product expression
  (e.g., `ROT_5 ∘ ROT_3`)? The index is concise but binds to the
  current order of the BFS-closure construction in `_close_group`;
  if that ever changes, every recorded index breaks. Matrix entries
  are robust but verbose. Generator-product expressions are
  middle-ground but require a parser.
- **σ-frame coordinate origin.** Frettlöh's figure may use the
  "blue vertex" of each prototile as the σ-frame origin; our JSON
  has the prototile's origin at vertex 1 of Table 1 (which is the
  class-I vertex for A, B, C and not in the same place for K).
  Confirm the σ-frame origin convention before transcribing.
