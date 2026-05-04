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
