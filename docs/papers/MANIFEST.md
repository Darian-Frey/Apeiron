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

All downloads verified with `file <path>` (PDF document) and `pdfinfo`
(page counts match expected venue pagination). No HTML error pages, no
sci-hub or pirate sources used.
