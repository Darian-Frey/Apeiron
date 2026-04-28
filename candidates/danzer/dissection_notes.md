# Danzer ABCK Substitution Dissections

**Status: Paolini transcription complete (2026-04-28, commit `b3d8e55`).
Q4a–Q4c partially answered locally; final ruling pending from Claude
(web). See "Status update" below.**

## Status update (2026-04-28, post-`b7174b1`)

The Paolini path has been fully executed locally:

- ✅ Symbolic+numerical evaluator built
  (`candidates/danzer/_paolini_extract.py`). Mirrors POV-Ray's
  `transform { ... }` composition: `scale<-1,-1,-1>`, `translate`,
  `Axis_Rotate_Trans` (Rodrigues), `Reorient_Trans`, and
  `AllineaTriangoliTrans` (with `AVOIDINSTABILITY=0` per the Q4c
  proposal — the 40°-randrot kludge is a numerical hack and is
  disabled in our exact-arithmetic reading).
- ✅ All 25 Paolini transforms evaluated. Every output rotation lands
  in I_h (proper or improper); every translation snaps to ZPhi (×2
  storage) within tolerance 1e-9. Saved as
  `candidates/danzer/paolini_dissection.json`.
- ✅ Volume conservation verified: sum of children's volumes = τ³·vol(parent)
  for all four σ-rules, to machine precision. Tested in
  `tests/integration/test_danzer_abck.py::TestPaoliniDissectionExtraction`.
- ✅ Per-parent child counts match Frettlöh's matrix exactly (11/7/5/2).
- ✅ Proper/improper distribution: 12 proper, 13 improper. Confirms
  Track A operates in I_h, not I.

**Q4a still requires Claude (web) confirmation that Paolini is the
canonical source.** If Claude (web) prefers Koca, the existing JSON
serves as Q4b cross-validation: per-child position-match under a
single global isometry between the two sources. We have not yet
determined the global isometry; that's a downstream task once Q4a is
locked.

**Q4c implicitly answered:** the `AVOIDINSTABILITY=1` flag is
disabled in our extractor. Every transform produces a clean I_h
isometry without the kludge — confirming the kludge is purely a POV-Ray
floating-point hygiene measure, not a mathematical correction.

This file is the working artifact for sub-commit 27B-β: the
real geometric dissection of σ(X) for X ∈ {A, B, C, K} in
Danzer's ABCK substitution.

## Sources located (general-purpose agent, 2026-04-28)

The original transcription target — Frettlöh's interactive viewer at
`www.math.uni-bielefeld.de/baake/frettloe/abck/inflabck.html` (his
ref [14]) — is **dead** (404). Two live, complete sources were
located instead:

1. **Koca, Koca, Al-Ajmi (arXiv 2003.13449, *Symmetry* 12(12) 1983)** —
   "Icosahedral Polyhedra from D₆ Lattice and Danzer's ABCK Tiling".
   Sections "Construction of τK / τB / τC / τA" (pp. 11–15) and
   Table 3 (pp. 15–17) give explicit child rotations and translations
   in matrix form, all in ½ℤ[τ] in the Cartesian frame — directly
   compatible with our denominator-2 storage convention.
2. **Maurizio Paolini's POV-Ray source** at
   `http://danzer.dmf.unicatt.it/`, SVN trunk
   `https://svn.dmf.unicatt.it/svn/projects/animations/danzer/trunk/`.
   Files `danzer_tiles.inc` (prototile vertex coordinates),
   `danzer_transformations.inc` (25 named transforms `tr{X}_{n}`,
   one per child), and `danzer_recursion.inc` (dispatch table giving
   each child's tile-type id and orientation flag with a leading
   minus sign for orientation-reversing children).

Local mirrors saved at `/tmp/d6_abck.txt`, `/tmp/frettloh_ikosa.txt`,
`/tmp/danzer_pov/*.inc`.

## Frame verification (Apeiron's perspective)

✅ **Vertex coordinates: VERIFIED.** Paolini's `pt0..pt9` in
`danzer_tiles.inc` match Frettlöh Table 1 exactly, which already
matches our `candidates/danzer/{A,B,C,K}.json`:

| vertex | Paolini | Frettlöh / our JSON |
| --- | --- | --- |
| pt0 | `<0, 0, 0>` | A.V0 = B.V0 = C.V0 = K.V0 = (0, 0, 0) |
| pt1 | `τ²<τ, 0, 1> = <τ³, 0, τ²>` | A.V1 = B.V1 |
| pt2 | `τ²<1, 1, 1>` | A.V2 = B.V2 = C.V2 |
| pt3 | `<τ², 1, 0>` | A.V3 |
| pt4 | `<τ², τ, 1>` | B.V3 |
| pt5 | `<-τ, 0, 1>` | C.V1 |
| pt6 | `<0, τ², 1>` | C.V3 |
| pt7 | `<-1, τ, 0>` | K.V1 |
| pt8 | `τ<1, 1, 1>` | K.V2 |
| pt9 | `½<-1, 1/τ, τ>` | K.V3 |

✅ **Koca's rotation entries: VERIFIED in I.** Eq. (19)'s g_K =
`[[0,-1,0], [0,0,-1], [1,0,0]]` has det +1 and, when stored as
`2·g_K`, is a member of our `apeiron.symmetry.ICOSAHEDRAL` set.
Koca's matrices in ½ℤ[τ] map directly to our denominator-2 storage.

✅ **Koca's translations: VERIFIED in ℤ[φ] under ×2 storage.** t_K =
½(τ, τ², 0) becomes `(ZPhi(0,1), ZPhi(1,1), ZPhi(0,0))` after the
×2 multiplication — pure ℤ[φ], no fractional component.

✅ **Substitution uses I_h, not just I.** Both sources confirm:
Paolini's `scale<-1,-1,-1>` (point inversion) appears in `trK_1` and
several other transforms; Koca's rotation matrices in eqs.
(20)–(31) include det = -1 entries. CLAUDE.md §2.1 has been
patched (commit `ffe72b5`) to record this finding.

## Critical disagreement: Koca and Paolini are not the same dissection

❌ **Koca's σ(K) and Paolini's σ(K) place the K-child at different
absolute positions within τK.** They are isometrically equivalent
(pairwise distances among K-child's transformed vertices match),
but the choice of which 3 of K's 4 vertices land on τK's corners
differs.

- **Koca eq. (19)+ text:** g_K rotates K so that K's face opposite
  D₂ lines up with B's face opposite D₀ inside τK; t_K = ½(τ, τ², 0).
  Resulting K-child has all 4 vertices interior to τK (none coincide
  with τK's corners).
- **Paolini `trK_2` macro:** `AllineaTriangoliTrans(pt7, pt9, pt8,
  τ·pt0, τ·pt9, τ·pt7)` maps K's V1 → τV0, V2 → τV1, V3 → τV3.
  So 3 of K-child's 4 vertices land at τK corners; only V0 is
  interior.
- **B-child poses also differ.** Koca's text says "origin coincides
  with one of the vertices of B" (Fig. 8), without giving an
  explicit rotation. Paolini's `trK_1` is `scale<-1,-1,-1> +
  translate pt2`, an `ImproperRotation(identity)` followed by
  `(τ², τ², τ²)`. Naive "B at identity at origin" assumption
  (without a rotation) puts B's vertex 1 = (τ³, 0, τ²) outside τK,
  so identity is wrong; Paolini's improper-rotation form is
  needed to keep B inside.

This is a real choice, not a bug: the same combinatorial
substitution rule (matrix [[0,0,1,0],[3,2,0,1],[2,1,2,0],[6,4,2,1]])
admits multiple valid geometric dissections. Both Koca's and
Paolini's encodings produce face-to-face dissections of τK with
the right tile counts, but the resulting σ-rule has different
recognisability behaviour in pillar 2.

Apeiron must pick one. **This is Q4 for Claude (web)** — see
relay below.

## Form to record (once Q4 is decided)

For each child of σ(X), record:

| field | type | meaning |
| --- | --- | --- |
| `prototile_index` | int 0–3 | A=0, B=1, C=2, K=3. |
| `translation` | `[[a, b], [a, b], [a, b]]` | Child's origin position in the σ(X)-frame, ×2-stored ZPhi. |
| `rotation` | `Rotation` or `ImproperRotation` | Child's isometry; orientation-preserving in I (60 elts) or orientation-reversing in I_h \\ I (60 elts). |

Convention: child placements are expressed in the σ(X) frame
*before* the inflation matrix is applied. The substitution rule's
inflation `φ·I` scales translations during multi-level expansion;
each child's translation here is its position in the level-1 frame.

## σ(A) — 11 children

```text
TODO — pending Q4 (Koca vs Paolini choice).
Once chosen, transcribe per-child (rotation, translation) from
the agreed source.
```

## σ(B) — 7 children

```text
TODO — pending Q4.
```

## σ(C) — 5 children

```text
TODO — pending Q4.
```

## σ(K) — 2 children — preliminary (Paolini path)

```text
PENDING Q4 — given here as a sample of the encoding, NOT committed.

Both sources agree this dissection contains:
- 1 B-child (orientation-reversing, parity = -1)
- 1 K-child (orientation-preserving for Koca; orientation-preserving for Paolini)

Paolini's encoding (machine-readable from danzer_transformations.inc):

  trK_1 → applied to tile B:
    operation:    scale<-1,-1,-1>; translate pt2 - pt0 = pt2
    Apeiron form: rotation = ImproperRotation(Rotation.identity())
                  translation_real = (τ², τ², τ²)
                  translation_×2   = (2τ², 2τ², 2τ²)
                                   = (ZPhi(2,2), ZPhi(2,2), ZPhi(2,2))

  trK_2 → applied to tile K:
    operation:    AllineaTriangoliTrans(pt7, pt9, pt8, τ·pt0, τ·pt9, τ·pt7)
    Apeiron form: TBD — requires symbolic evaluation of POV-Ray
                  Reorient_Trans macro composition; see Q4.

Koca's encoding (eq. 19):

  K-child:
    g_K_real         = [[0,-1,0], [0,0,-1], [1,0,0]]   det = +1, ∈ I
    g_K_×2_stored    = [[0,-2,0], [0,0,-2], [2,0,0]]   ∈ ICOSAHEDRAL
    t_K_real         = ½(τ, τ², 0)
    t_K_×2_stored    = (τ, τ², 0)
                     = (ZPhi(0,1), ZPhi(1,1), ZPhi(0,0))

  B-child: pose left to Fig. 8, not algebraically stated.
```

## Cross-checks once populated

- Per-tile child count sums to the column sum of the substitution
  matrix (already verified by `test_matrix_column_sums_match_child_counts`).
- After encoding, every child's rotation must be in
  `apeiron.symmetry.ICOSAHEDRAL` (proper) or constructible as
  `ImproperRotation(g)` with `g ∈ ICOSAHEDRAL`.
- Children of σ(X) tile φ·X without overlap and without gaps; a
  validation pass on the assembled level-1 patch should confirm
  (face-to-face adjacency, no interior overlap).
- After encoding lands, sub-commit 27D runs the full pipeline:
  patch_from_supertile → is_recognisable (with shell signature) →
  inflation_argument → pillar-4 stub.

## Q4 to relay to Claude (web)

> **Question for Claude (web): Which dissection encoding for Track A
> — Koca or Paolini?**
>
> Frame-verification findings (locally, against `candidates/danzer/*.json`):
>
> - Frettlöh Table 1 vertex coordinates ↔ our JSONs: VERIFIED.
> - Koca's matrices in ½ℤ[τ] ↔ our `Rotation` storage convention:
>   compatible (entries map cleanly to ZPhi after ×2).
> - Koca eq. (19)'s g_K ∈ Apeiron's `ICOSAHEDRAL` set: VERIFIED.
> - I_h finding (improper rotations in the rule): VERIFIED in both
>   sources; CLAUDE.md §2.1 patched.
>
> **The problem:** Koca and Paolini place σ(K)'s K-child at
> isometrically-equivalent but different absolute positions within
> τK. Same combinatorial substitution matrix, same total dissection,
> but a different choice of which 3 K-child vertices coincide with
> τK's corners. Both are valid face-to-face dissections; pillar-2
> recognisability behaviour may differ.
>
> Three sub-questions:
>
> **Q4a (source choice):** Encode Koca or Paolini? Trade-off:
>
> - **Koca**: explicit (g, t) matrices for K, B, C, A children given
>   in eqs. (19)–(31). Direct ZPhi transcription. But: B's pose in
>   σ(K) is left to Fig. 8 (no equation), and the pose interpretation
>   I tried (B at identity at origin) puts B's vertex 1 outside τK.
>   So either Koca's text is incomplete or implicit; figure-reading
>   is required for at least one child.
> - **Paolini**: complete machine-readable POV-Ray source for all
>   25 children. But: many transforms use the
>   `AllineaTriangoliTrans` macro, which composes
>   `Reorient_Trans(N1, N2)` (a vector-alignment rotation) plus
>   `Axis_Rotate_Trans(rotAB, angle)` (a residual rotation). To
>   extract `(g, t)` for each child, I'd need to symbolically
>   evaluate these macros over ℚ(√5), which is non-trivial. Some
>   transforms also include a `randrot` "instability avoidance"
>   40°-rotation kludge that's clearly a numerical hack, not a
>   clean isometry — would need to disable.
>
> Recommendation: Paolini for completeness, but pre-process with a
> short symbolic-evaluation Python script (using SymPy in ℚ(√5)) to
> extract `(g, t)` for each `tr*_n` macro. Concur, override, or
> propose a third path?
>
> **Q4b (validation strategy):** Once one source is chosen,
> validation via the OTHER source means symbolically computing the
> isometry between Koca's and Paolini's K-child placements (they're
> isometrically equivalent, so a global isometry exists). If they
> disagree on the dissection topology (which children share which
> faces), then one of them has a bug — this is the test. Should the
> cross-source check be:
>
> - (a) per-child position-matching after applying a single global
>   isometry, OR
> - (b) topological face-adjacency check that's invariant under
>   global isometry?
>
> **Q4c (Paolini's `randrot` kludge):** Paolini's macro
> `AllineaTriangoliTrans` includes an `#ifdef AVOIDINSTABILITY` branch
> that interposes a 40°-rotation around the second normal. This is
> floating-point numerical hack. The flag is set globally
> (`#declare AVOIDINSTABILITY=1` in `danzer_macros.inc`). For exact
> ZPhi transcription, the kludge should be DISABLED — the macro
> reduces to the cleaner unstable-but-exact form. Confirm? If so,
> the symbolic evaluator I'd write would explicitly set
> `AVOIDINSTABILITY=0`.

## Outstanding decisions (post-Q4)

- **Rotation labelling in JSON.** Once child rotations are recorded,
  store as ICOSAHEDRAL index (compact but breaks if BFS order
  changes), as raw matrix entries (verbose but robust), or as a
  generator-product expression. Defer until Q4 lands; the choice
  is downstream.
- **σ-frame coordinate origin.** Confirmed: σ(X)'s frame has its
  origin at X's V0 (= (0,0,0) per Frettlöh Table 1, per
  `candidates/danzer/*.json`'s vertex 0). This is also Paolini's
  `pt0` and Koca's `D₀`. No frame disagreement here.
