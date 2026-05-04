# Apeiron literature deep-dive — paper manifest

Triggered after exhaustive computational search produced zero realisable
substitution monotile candidates. Each entry: title, authors, year, venue,
source URL, local filename, one-line summary of relevance.

---

## 1. Goodman-Strauss, *Matching Rules and Substitution Tilings* (1998)

- **Authors:** Chaim Goodman-Strauss
- **Year:** 1998
- **Venue:** Annals of Mathematics, 147(1), 181–223
- **DOI:** 10.2307/120988
- **Source URL:** https://strauss.hosted.uark.edu/papers/MRandST.pdf
  (author's University of Arkansas faculty page; open-access PDF, 45 pages,
  matches journal pagination)
- **Local file:** `goodman_strauss_1998_matching_rules.pdf`
- **Why on the read list:** The theoretical backbone of the four-pillar proof
  structure (CLAUDE.md §6.3). Establishes that any substitution tiling
  satisfying mild conditions admits local matching rules forcing the
  hierarchy — i.e. recognisability + hierarchical forcing in a single
  framework. Highest-priority read for both Track A and Track B.

---

## 2. Frettlöh, *Substitution tilings with statistical circular symmetry* (2008, EJC)

- **Authors:** Dirk Frettlöh
- **Year:** 2008
- **Venue:** European Journal of Combinatorics, 29(8), 1881–1893
- **arXiv:** math/0704.2521
- **Source URL:** https://arxiv.org/pdf/0704.2521 (arXiv preprint, 13 pages)
- **Local file:** `frettloh_2008_substitution_tilings_circular_symmetry.pdf`
- **Why on the read list:** Comprehensive treatment of substitution tilings
  whose tile-orientation distribution becomes circularly symmetric in the
  limit. Relevant to assessing whether icosahedrally-symmetric 3D
  substitutions can produce statistically isotropic hulls — a structural
  property that bears on the choice of inflation factor (φ² vs φ³) and on
  Pillar 4 (no non-hierarchical tilings).

---

## 3. Frettlöh, *About substitution tilings with statistical circular symmetry* (2008, conference)

- **Authors:** Dirk Frettlöh
- **Year:** 2008
- **Venue:** Conference proceedings (companion exposition; cites the EJC
  paper as [6] for the main theorem)
- **arXiv:** 0803.2172 (math.SP, 14 Mar 2008)
- **Source URL:** https://arxiv.org/pdf/0803.2172 (7 pages)
- **Local file:** `frettloh_2008_statistical_circular_symmetry.pdf`
- **Why on the read list:** Shorter expository companion to entry 2. Adds
  diffraction-spectrum and dynamical-systems context: tilings with
  statistical circular symmetry have circular-symmetric diffraction with
  trivial pure-point part. Useful framing for what the diffraction signature
  of a 3D icosahedral monotile would look like.

**Note on the user's "direct product variation" prompt:** the DPV
construction is due to Frank & Robinson, not Frettlöh; Frettlöh's published
work on tilings produced by direct-product–like constructions is the
statistical-circular-symmetry pair above (the broader 2008 EJC paper
includes examples in this family). The dedicated *Fibonacci direct product
variation tilings* paper (arXiv:2203.07743, J. Math. Phys. 2022, Baake et
al.) is the canonical DPV reference if needed later — not downloaded here
because Frettlöh is not an author and the construction is 1D-restricted.

---

## 4. Frettlöh & Harriss, *Parallelogram tilings, worms, and finite orientations* (2013)

- **Authors:** Dirk Frettlöh, Edmund O. Harriss
- **Year:** 2013 (preprint Feb 2012)
- **Venue:** Discrete & Computational Geometry, 49(3), 531–539
- **DOI:** 10.1007/s00454-012-9478-5
- **arXiv:** 1202.4686
- **Source URL:** https://arxiv.org/pdf/1202.4686 (9 pages including refs)
- **Local file:** `frettloh_harriss_2013_parallelogram_tilings.pdf`
- **Why on the read list:** Proves that any plane parallelogram tiling with
  finite protoset has tiles in only finitely many orientations
  (Theorem 2.5). The proof (Crossing Lemma + Travel Lemma via "worms")
  generalises in spirit to higher dimensions, with explicit remarks on the
  R^d (parallelotope) extension on p.7. Directly relevant to
  finite-vs-infinite orientation classes for an icosahedrally symmetric 3D
  monotile, and to the n in the n×n substitution matrix.

---

## 5. Walton & Whittaker, *An aperiodic tile with edge-to-edge orientational matching rules* (2019)

- **Authors:** James J. Walton, Michael F. Whittaker
- **Year:** 2019 submission (v2 Oct 2021); published in Journal of the
  Institute of Mathematics of Jussieu
- **arXiv:** 1907.10139 (math.MG)
- **Source URL:** https://arxiv.org/pdf/1907.10139 (27 pages)
- **Local file:** `walton_whittaker_2019_orientational_matching_rules.pdf`
- **Why on the read list:** Constructs a hexagonal monotile whose
  aperiodicity is enforced by genuinely nearest-neighbour (edge-to-edge)
  orientational matching rules — improving on Socolar–Taylor, where the
  rules are non-nearest-neighbour. Directly addresses the local vs nonlocal
  enforcement question. Hull is uniquely ergodic with pure point dynamical
  spectrum and regular model set structure.

**Attribution note:** the original prompt requested "Socolar 2019" on this
theme. No such Socolar paper exists; Socolar's work on this question is the
2010/2011 Socolar–Taylor papers (arXiv:1003.4279, arXiv:1009.1419). The
2019 paper that most directly continues the Socolar–Taylor programme on
local-vs-nonlocal enforcement is this Walton–Whittaker paper, which
explicitly cites Socolar–Taylor as the immediate predecessor and improves
on its locality. Substituted accordingly. If a Socolar-authored paper is
specifically required, recommend the 2011 *J. Combin. Theory A* hexagonal
tile paper (arXiv:1003.4279) as the canonical reference for that theme.

---

## 6. Kramer & Neri, *On periodic and non-periodic space fillings of E^m obtained by projection* (1984) — **NOT DOWNLOADED**

- **Authors:** Peter Kramer, Reinhardt Neri
- **Year:** 1984
- **Venue:** Acta Crystallographica Section A, 40(5), 580–587
- **DOI:** 10.1107/S0108767384001203
- **Erratum:** Acta Cryst. A41 (1985), DOI 10.1107/S0108767385001350
- **Source URL (paywalled):** https://doi.org/10.1107/S0108767384001203
  (IUCr abstract page: https://journals.iucr.org/paper?S0108767384001203;
  publisher PDF requires subscription; Wiley mirror likewise paywalled).
- **Local file:** *not downloaded* (closed-access; no open-access copy
  exists per Unpaywall and OpenAlex queries on 2026-05-04).
- **Why on the read list (Q11c / Q14b):** seminal paper introducing the
  D6 → R^3 cut-and-project method for icosahedral quasilattices — the
  upstream construction underlying Frettlöh's "Icosahedral tilings in R³"
  Theorem 1.3 result that ABCK is a model set from D6 via CPS. Q11c
  (algebraic CPS structure) and Q14b (window characterisation) require
  the original window definition from §3 of this paper, which subsequent
  papers cite without re-deriving.

---

## 7. Kramer, Papadopolos, Schlottmann & Zeidler, *Projection of the Danzer tiling* (1994) — **NOT DOWNLOADED**

- **Authors:** Peter Kramer, Zorka Papadopolos, Martin Schlottmann,
  Dieter Zeidler
- **Year:** 1994
- **Venue:** Journal of Physics A: Mathematical and General, 27(13),
  4505–4517
- **DOI:** **10.1088/0305-4470/27/13/024** (correction: the user-provided
  DOI 10.1088/0305-4470/27/13/021 actually resolves to Wong–Tomonaga,
  *Shell-model matrix elements in a neutron-proton quasi-spin formalism* —
  unrelated. The correct article ID for the Danzer projection paper is
  /024, verified via OpenAlex exact-title search.)
- **Source URL (paywalled):** https://doi.org/10.1088/0305-4470/27/13/024
  (IOPscience landing page; PDF requires subscription).
- **Local file:** *not downloaded* (closed-access; no open-access copy
  exists per Unpaywall and OpenAlex queries on 2026-05-04. No arXiv
  preprint exists — Papadopolos's earliest arXiv submission is 1998
  (math-ph/9909012, *Tiles–Inflation Rules for the Class of Canonical
  Tilings T*(2F)*); Schlottmann never uploaded to arXiv. Kramer's
  Tübingen publications page lists the paper but has no PDF link.)
- **Why on the read list (Q11c / Q14b):** the explicit construction of
  the Danzer ABCK tiling as a D6 → R^3 cut-and-project model set, with
  windows fully specified in §3–4. This is the paper that bridges
  Kramer–Neri 1984 (general projection method) and the modern statement
  in Frettlöh "Icosahedral tilings in R³" Theorem 1.3 that triggered the
  Q9c gate. Q14b (window characterisation in the internal space) cannot
  be answered from secondary sources alone.

### Recommended open-access substitutes

Both papers are subscription-only with no preprint, no faculty-page PDF,
no arXiv version, and no Wayback capture of the PDF (Wayback only
captured the IUCr 301/401 redirect page, never the actual PDF). The
following open-access papers cover the same material and should be read
in their place pending acquisition through institutional access:

1. **Al-Siyabi, Koca & Koca (2020),** *Icosahedral Polyhedra from D6
   lattice and Danzer's ABCK tiling*, arXiv:2003.13449 — explicitly
   reconstructs the ABCK tiling as a D6 projection with full window
   specifications (eqs. 19–31), and is already cited in CLAUDE.md §2.1.
   Open-access PDF: https://arxiv.org/pdf/2003.13449.
2. **Frettlöh, Harriss & Gähler (eds.), Tilings Encyclopedia entry
   "Danzer's ABCK tilings",** open-access at
   https://tilings.math.uni-bielefeld.de/substitution/danzers-abck-tiling/
   — consolidates the Kramer-Papadopolos-Schlottmann-Zeidler results in
   modern notation.
3. **Frettlöh,** *Icosahedral tilings in R³: The ABCK tilings*,
   https://www.math.uni-bielefeld.de/~frettloe/papers/ikosa.pdf —
   contains Theorem 1.3 (model-set status of ABCK) and references
   Kramer et al 1994 for the projection details. Open-access.

### Acquisition options for the originals

- Institutional library subscription to Acta Crystallographica A and
  J. Phys. A.
- Inter-library loan via any university library.
- Direct request to Z. Papadopolos (Universität Tübingen, retired) or
  the Kramer estate. Kramer (1933–2024) maintained an active
  publications page until ca. 2010; Papadopolos's contributions are
  archived at researchain.net.

Sci-hub and pirate sources were **not** used and should not be used.

---

## 8. Frettlöh, *Icosahedral tilings in R³: The ABCK tilings* (Bielefeld internal report)

- **Authors:** Dirk Frettlöh
- **Year:** undated (post-2008, pre-2013 by reference list)
- **Venue:** Bielefeld internal report
- **Source URL:** https://www.math.uni-bielefeld.de/~frettloe/papers/ikosa.pdf
  (author's faculty page; open-access PDF, 9 pages including refs)
- **Local file:** `frettloh_icosahedral_ABCK_D6_CPS.pdf`
- **Why on the read list (Q14b):** Theorem 1.3 — ABCK tilings are model
  sets from D₆ via CPS, with three explicit windows (regular dodecahedron,
  great dodecahedron, custom 5-fold-star polyhedron). Most direct
  open-access substitute for §7.4 of Baake-Grimm Vol. 1. Read in full
  during Q14b post-(c.1) literature phase; notes at
  `docs/literature_notes.md` §6 + §7.2.

---

## 9. Moody, *Model sets: a survey* (2000)

- **Authors:** Robert V. Moody
- **Year:** 2000
- **Venue:** Survey article (chapter in *From Quasicrystals to More
  Complex Systems*, Centre de Physique des Houches, ed. F. Axel,
  F. Dénoyer, J.-P. Gazeau)
- **arXiv:** math/0002020
- **Source URL:** https://arxiv.org/pdf/math/0002020 (28 pages)
- **Local file:** `moody_2000_model_sets_survey.pdf`
- **Why on the read list (Q15 substitute for Baake-Grimm Ch. 7):** "The
  canonical model-set survey" cited by Baake-Grimm Vol. 1; the most
  expository open-access treatment of cut-and-project schemes and model
  sets in the mathematical-quasicrystals tradition. Strongest substitute
  for §7.2 (Cut and project schemes and model sets).

---

## 10. Baake, Damanik & Grimm, *What is... aperiodic order?* (2015)

- **Authors:** Michael Baake, David Damanik, Uwe Grimm
- **Year:** 2015
- **Venue:** Notices of the American Mathematical Society
- **arXiv:** 1512.05104
- **Source URL:** https://arxiv.org/pdf/1512.05104 (5 pages)
- **Local file:** `baake_damanik_grimm_2015_what_is_aperiodic_order.pdf`
- **Why on the read list (Q15 substitute):** Light Notices-style
  introduction by two of the three Baake-Grimm Vol. 1 authors; surveys
  the same conceptual ground at a much higher level. Useful for
  cross-referencing language and framing only.

---

## 11. Richard & Strungaru (Baake-cite), *A short guide to pure point diffraction in cut-and-project sets*

- **Authors:** Christoph Richard, Nicolae Strungaru
- **Source URL:** open-access PDF (32 pages)
- **Local file:** `baake_short_guide_pure_point_diffraction_CPS.pdf`
- **Why on the read list (Q15 substitute):** Diffraction-focused; covers
  the spectral side of CPS in detail. Less directly relevant to the
  structural-multi-prototile question but worth keeping for the
  diffraction angle.

---

## Note on the Baake-Grimm Vol. 1 Ch. 7 acquisition

The target chapter (Baake & Grimm 2013, *Aperiodic Order, Vol. 1*, Ch. 7
"Tilings, Coverings and Their Diffraction") is **closed-access**.

Cambridge Core hosts no open sample chapter. Google Books preview was
inaccessible at the time of the agent's fetch. Author homepages do not
host draft chapters.

TOC confirmed via the CUP frontmatter PDF
(https://assets.cambridge.org/97805218/69911/frontmatter/9780521869911_frontmatter.pdf,
16 pages, open-access):

- 7.1 Silver mean chain via projection (p. 251)
- 7.2 Cut and project schemes and model sets (p. 263)
- 7.3 Cyclotomic model sets (p. 272)
- **7.4 Icosahedral model sets and beyond (p. 283)** — directly relevant
- 7.5 Alternative constructions (p. 289)

**Acquisition options:**

- UK: Cambridge UL, Oxford Bodleian, Open University (Grimm's home),
  Bielefeld (Baake's home).
- US: any R1 mathematics library (MIT, Stanford, Berkeley, UIUC,
  Princeton).
- Inter-library loan via any university library is the standard path.

**Fallback if institutional access not pursued:** treat entries §6 (Frettlöh
ABCK), §7 (Al-Siyabi-Koca-Koca 2020), §9 (Moody 2000) as the working
substitute for §7.4 + §7.2.

---

## 12. Hayashi et al. (2011), *Notes on Vertex Atlas of Danzer Tiling*

- **Authors:** Hiroko Hayashi, Yuu Kawachi, Kazushi Komatsu, Aya
  Konda, Miho Kurozoe, Fumihiko Nakano, Naomi Odawara, Rika Onda,
  Akinobu Sugio, Masatetsu Yamauchi
- **Year:** 2011
- **Venue:** Nihonkai Math. J. 22 (2011), 49–58
- **Source URL:** Local file (origin unspecified)
- **Local file:** `1339694050.pdf`
- **Topic:** **Different Danzer tiling** — the 2D Nischke-Danzer
  7-fold construction with 6 prototiles (three 2D triangles and
  their reflections). Theorem 1.1: 39 vertex atlases with arrows,
  29 without. Discussion of up-down generation and comparison
  with Penrose 5-fold.
- **Why on the read list (post-Q16 user contribution):** Added by
  user 2026-05-04 as a candidate substitute for
  Kramer-Papadopolos-Schlottmann-Zeidler 1994. **Off-topic:** does
  not address 3D Danzer ABCK, D₆ cut-and-project, window
  deformation, or the structural multi-prototile conjecture. The
  shared keyword "Danzer" refers to a different tiling from the
  same author. Documented here for completeness but not used in
  the synthesis.

---

## 13. Papadopolos, Hohneker & Kramer (1999/2000), *Tiles–inflation rules for the class of canonical tilings T\*(2F)*

- **Authors:** Z. Papadopolos, C. Hohneker, P. Kramer
- **Year:** 1999 (preprint), 2000 (published)
- **Venue:** Discrete Mathematics 221 (2000) 101–112
- **arXiv:** math-ph/9909012
- **Source URL:** https://arxiv.org/pdf/math-ph/9909012 (24 pages)
- **Local file:** `papadopolos_hohneker_kramer_1999_t2f_inflation.pdf`
- **Why on the read list (Q15+ supplementary, post-Path A):** the
  open-access successor to closed-access Kramer-Papadopolos-
  Schlottmann-Zeidler 1994. Identifies T*(2F) as the canonical D₆
  icosahedral projection (parent of ABCK), with **6 prototile
  tetrahedra** (or 8 with mandatory blue/red colour decoration on
  C and G). Provides verbatim 6×6 projection matrix B (eq. 1) and
  inflation matrix I_D6 (eq. 2) over ℤ[τ]. Read in full; notes at
  `docs/literature_notes.md` §9. Strengthens the structural
  multi-prototile conjecture by establishing the canonical
  projection requires ≥ 6 prototiles, with ABCK being a local
  reduction to 4.

---

## 14. Hammock, Fang & Irwin (2018), *Quasicrystal Tilings in Three Dimensions and Their Empires*

- **Authors:** Dugan Hammock, Fang Fang, Klee Irwin
- **Year:** 2018
- **Venue:** Crystals 8(10):370
- **DOI:** 10.3390/cryst8100370
- **Source URL:** https://www.mdpi.com/2073-4352/8/10/370 (CC-BY,
  16 pages)
- **Local file:** `hammock_fang_irwin_2018_quasicrystal_empires.pdf`
- **Why on the read list (Q15++ supplementary, 2026-05-05):**
  computes acceptance-domain sector decomposition for T^(D6) =
  T*(2F): **880 tile types** (regions of W), **4230 sectors**,
  **36 vertex configurations** (matching Papadopolos eq. 7
  exactly). Table 1 enumerates all 36 with sectors, **empires
  (forced tiles)**, and exact ℤ[√5] frequencies. Critical
  finding for Apeiron: **none of the 36 empires reduces to a
  single-prototile region** — strongest empirical evidence yet
  for the structural multi-prototile conjecture. §4 (eqs. 16-24)
  gives the empire computation methodology directly translatable
  to ZPhi³ implementation. Read in full; notes at
  `docs/literature_notes.md` §10.

---

## 15. Papadopolos, Klitzing & Kramer (1997), *Quasiperiodic icosahedral tilings from the six-dimensional bcc lattice* — ILP PRIORITY 1

- **Authors:** Z. Papadopolos, R. Klitzing, P. Kramer
- **Year:** 1997
- **Venue:** J. Phys. A: Math. Gen. 30 (1997) L143–L147
- **DOI:** 10.1088/0305-4470/30/5/004
- **Source URL:** closed-access (IOP). Length ~5 pp (Letter).
- **Local file:** *not downloaded* (closed-access; ILL pending).
- **Why on the read list (Q16++ supplementary, 2026-05-05):** per
  Claude (web), this Letter "defines BOTH T\*(2F) and T(2F) via
  Klotz construction from D₆", determines acceptance domains for
  both classes, and crucially establishes that **"both tilings
  can be locally transformed into tilings with only four elements
  that allow simple inflation"**. This is where the 6-tile →
  4-tile reduction (T\*(2F) → ABCK) is published; reading it
  clarifies the local-derivation mechanism that ABCK uses to
  reduce from the canonical projection's 6 prototiles to 4.
  Prerequisite to asking whether further 4 → 1 reduction is
  possible (the central question for Track A / Route 1).

---

## 16. Kramer & Papadopolos (1994b), *The Mosseri–Sadoc tilings derived from the root lattice D₆* — ILP PRIORITY 2

- **Authors:** P. Kramer, Z. Papadopolos
- **Year:** 1994
- **Venue:** Canadian Journal of Physics 72(7-8):408–414
- **DOI:** 10.1139/p94-057
- **Source URL:** closed-access (NRC Press). Length ~7 pp.
- **Local file:** *not downloaded* (closed-access; ILL pending).
- **Why on the read list (Q16++ supplementary, 2026-05-05):**
  methodological parallel to the closed-access 1994 Danzer
  projection paper (same authors, same year, same "determine all
  windows" objective, different tiling). Reading this gives the
  method for window enumeration in the D₆ family. If the method
  is clear from this paper, the Danzer window enumeration may be
  reproducible computationally without the 1994 Danzer paper.
  Add to the same ILP request as the 1994 Danzer paper and the
  1997 bcc-lattice Letter; all three should be obtainable from a
  single library visit.

---

## ILP request — three closed-access D₆ projection papers

For institutional inter-library loan via any UK or US R1
mathematics/physics library:

1. **Kramer, P., Papadopolos, Z., Schlottmann, M. & Zeidler, D.**
   "Projection of the Danzer tiling." *J. Phys. A: Math. Gen.*
   27 (1994) 4505–4517. DOI: 10.1088/0305-4470/27/13/024.
2. **Papadopolos, Z., Klitzing, R. & Kramer, P.** "Quasiperiodic
   icosahedral tilings from the six-dimensional bcc lattice."
   *J. Phys. A: Math. Gen.* 30 (1997) L143–L147. DOI:
   10.1088/0305-4470/30/5/004.
3. **Kramer, P. & Papadopolos, Z.** "The Mosseri–Sadoc tilings
   derived from the root lattice D₆." *Canadian Journal of
   Physics* 72 (1994) 408–414. DOI: 10.1139/p94-057.

The first contains the "all windows" Danzer result that would
clear Gate C for the structural multi-prototile conjecture. The
second and third provide methodological context (T\*(2F) → ABCK
reduction mechanism; Mosseri–Sadoc parallel "all windows"
computation).

---

## Note on duplicates added 2026-05-04

User contributed five files; three are duplicates of papers
already in this manifest, one is a duplicate of the Frettlöh
ikosa paper, and one is the Hayashi 2D paper (entry §12 above).
Verified by md5:

- `2003.13449v3.pdf`,
  `Icosahedral_Polyhedra_from_D6_lattice_and_Danzers.pdf`,
  `al_siyabi_koca_koca_2020_icosahedral_polyhedra_d6.pdf` — all
  md5 `0e0af0f0a08d2b0bf014ab719a4c0b79`. Same arXiv content
  (v3 ≡ v1 byte-identically; the apparent v1 we initially
  downloaded is in fact the v3).
- `ikosa.pdf`, `frettloh_icosahedral_ABCK_D6_CPS.pdf` — both md5
  `e214d8fd3b06dff6038b45e16d485bf4`. Same Frettlöh paper.
- `symmetry-12-01983-v2.pdf` — published *Symmetry* 12 (12) 1983
  version of arXiv:2003.13449. 15 pages (vs arXiv v1's 19 pages,
  due to journal two-column formatting). Subject metadata
  identical to the arXiv version's abstract. Treat as the same
  paper for synthesis purposes; no re-read required.

---

## Source-discovery notes

- **Goodman-Strauss:** open-access PDF on Arkansas faculty page;
  `https://strauss.hosted.uark.edu/papers/MRandST.pdf` resolved cleanly.
  No need for Wayback or JSTOR.
- **Frettlöh papers:** arXiv preprints found via search of his Bielefeld
  research page and arXiv author search.
- **Frettlöh & Harriss:** arXiv:1202.4686 located via search "Frettloh
  Harriss parallelogram tilings 2013".
- **Walton & Whittaker:** found via the literal search for the
  "local-vs-nonlocal enforcement" theme on Socolar's arXiv list (which
  showed no 2019 paper) and the 2019 arXiv listing for the relevant
  keywords.
- **Kramer & Neri 1984; Kramer, Papadopolos, Schlottmann & Zeidler 1994:**
  closed-access. Searched IUCr (Cloudflare-protected, subscription-gated),
  IOPscience (subscription-gated), Wiley Online Library (subscription),
  arXiv (no preprint exists for either paper), Wayback Machine (no PDF
  capture; only redirect/401 captures of IUCr), Unpaywall API
  (`is_oa=false` for both), OpenAlex API (`oa_status=closed` for both,
  `any_repository_has_fulltext=false`), Google Scholar (no [PDF] link),
  ResearchGate (publication-page access blocked), Semantic Scholar
  (`openAccessPdf.status=CLOSED`), Z. Papadopolos's Tübingen and
  ResearchGate pages, Kramer's Tübingen publications listing
  (`homepages.uni-tuebingen.de/peter.kramer/Sx1.SSx1.html` — paper is
  listed in bibliography but no PDF link). All paths exhausted.

All successful downloads verified with `file <path>` (PDF document) and
`pdfinfo` (page counts match expected venue pagination). No HTML error
pages, no sci-hub or pirate sources used.
