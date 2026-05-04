# Literature review — Q10 Output 2 (synthesis)

**Status.** Synthesis of the four-paper read recorded in
[literature_notes.md](literature_notes.md). Per the Q10 ruling
(2026-04-29), this is Output 2 of the post-Q9c-gate literature
deep-dive. After this commit, re-relay to Claude (web) with the
synthesis result before any new code lands.

**Authoritative claim level.** Section (a) reports what the
literature says; section (b) interprets it against Apeiron's
zero-Realised result; section (c) recommends the next move. The
recommendation is a judgment call — the relay-back is what
authorises any pivot.

---

## Headline result (post-Q14a, 2026-05-04)

**Precisely-scoped negative result.** Per Q14a's ruling, the
empirical finding is now stated with explicit scope:

> No icosahedral substitution monotile with PF eigenvalue φ³
> exists among n ≤ 3 prototile rules with integer entries ≤ 3 and
> prototile shapes drawn from the 9 H₃-compatible tetrahedral
> classes derivable from the Paolini ABCK vertex pool. This is the
> first systematic computational survey of this family.

Evidence: 78 candidates (13 at max_entry=2, 65 at max_entry=3) ran
exhaustive realisation against every volume-compatible shape
triple in the 9-class pool; every one returned NoRealisation
within the per-step DFS pruning + 30s/triple budget at k_max=7.
Plus 9 + 144 = 153 SKIPPED_TAXONOMY_GAP candidates whose left-
eigenvector ratios fall outside the 9-class volume set; their
realisability is the *separate* question that the 15-class
extension would address. Logs at
[notebooks/n3_phi3_max_entry_3_2026-05-04.log](../notebooks/n3_phi3_max_entry_3_2026-05-04.log)
and the prior `notebooks/n3_phi3_exhaustive_2026-04-29.log`.

**What this scoping does NOT claim.** Not "no 3D icosahedral
monotile exists." Not "the icosahedral substitution route is
exhausted." Not "the H₃ taxonomy is the wrong shape family." The
result holds within the named slice; broader claims require either
(a) the 15-class extension, (b) larger max_entry, or (c) a
genuinely different shape family (e.g., the cut-and-project pivot
under investigation per Q14b).

---

## (a) What the literature says about 3D icosahedral monotile substitution

**Direct results: none.** None of the four papers in the read
list addresses the 3D icosahedral monotile substitution problem.
Each addresses a related but distinct setting:

| Paper                        | Setting             | Result type                    |
|------------------------------|---------------------|--------------------------------|
| Goodman-Strauss 1998         | E^d, d > 1, *any* substitution tiling | Sufficiency: substitution → matching rules     |
| Frettlöh-Harriss 2013        | ℝ², parallelograms only          | Finite-orientations theorem                  |
| Frettlöh 2008 EJC + companion| ℝ², infinite-orientation substitutions | Equidistribution and circular symmetry       |
| Walton-Whittaker 2019        | ℝ², hexagonal monotile           | Edge-to-edge orientational matching example  |

Goodman-Strauss covers d > 1 and is the only one whose theorem
formally applies to ℝ³. The other three are 2D-specific.

**Indirect results bearing on Apeiron's setting:**

1. **Goodman-Strauss's main theorem is sufficiency-only.** It
   asserts that *if* a substitution tiling satisfies sibling-edge-
   to-edge with hereditary edges, *then* it can be enforced by
   matching rules on a finite set of marked prototiles. It says
   nothing about which substitutions exist, nor whether a single
   unmarked prototile can carry the substitution.

2. **The construction in Goodman-Strauss does not preserve
   monotile-ness.** Even if the input substitution has n = 1
   prototiles, the output decomposes into vertex tiles, edge
   tiles, small tiles, and big tiles (§3.1). The theorem is about
   *tile sets*, not monotiles.

3. **The "mild condition" is stricter in 3D than 2D.** Lemma 1.3
   (and the author's note on p.10) shows hereditary-vertex does
   not imply hereditary-edge in d > 2. An example is anticipated
   but not exhibited as of 1998.

4. **The locking mechanism (vertex wires) requires shared sibling
   vertices** (§4, p.37). The author leaves open whether
   substitutions exist where sibling positions cannot be fixed at
   all — i.e., whether there are substitutions that resist
   matching-rule enforcement entirely.

5. **Frettlöh-Harriss's parallelogram theorem does not generalise
   to 3D for free.** The author explicitly notes (p.7) that
   extending Lemma 2.4 to ℝ^d, d ≥ 2, by parallelotopes "would
   require further objects... pseudo-planes (and pseudo-
   hyperplanes) made up of linked parallelotopes sharing a
   parallel edge in common." The 3D case is sketched, not proved.
   And tetrahedra (Apeiron's H₃-compatible candidates) are NOT
   parallelotopes; even the proven 2D theorem does not apply.

6. **Frettlöh 2008 EJC + companion classify pinwheel-like
   substitutions.** Proposition 3.4 (EJC) gives the test: a
   primitive substitution σ in ℝ² is pinwheel-like iff some
   σ^n(T_i) contains two same-type tiles with α(T) − α(T') ∉ πℚ.
   Apeiron's icosahedral substitutions live in the order-60
   icosahedral group I (or order-120 I_h); all rotation-angle
   differences are rational multiples of π. **Therefore Apeiron's
   substitutions are NOT pinwheel-like;** the Frettlöh 2008
   machinery does not apply.

7. **Walton-Whittaker 2019 is an existence proof for a 2D marked
   monotile** (hexagonal, with edge decorations enforcing
   aperiodicity by orientational matching rules). It does not
   prove an *unmarked* 2D monotile exists, nor does it address 3D.
   Its relevance is meta: it demonstrates that decoration is a
   valid route to aperiodicity when shape alone does not enforce
   it. (The hat/spectre 2023 papers, not in the read list, do
   prove an unmarked 2D monotile exists by shape alone — but for
   3D no such result exists.)

**Bottom line of (a):** the literature has no direct theorem
about 3D icosahedral monotile substitution; the closest results
are sufficiency-only (Goodman-Strauss) or 2D-specific (Frettlöh,
Walton-Whittaker, Frettlöh-Harriss).

---

## (b) Is any known obstruction consistent with the exhaustive NoRealisation result?

**Apeiron's exhaustive-sweep result.** 170 (cand, triple) pairs
across 13 candidates from the n=3 PF=φ³ algebraic survey
(max_entry=2). Result: 0 Realised, 13 (all triples) NoRealisation,
0 Inconclusive, 2 Skipped (k > k_max=5). Source:
[notebooks/n3_phi3_exhaustive_2026-04-29.log](../notebooks/n3_phi3_exhaustive_2026-04-29.log).

**No obstruction in the read literature explains this.** None of
the four papers contains a theorem that would predict zero
realisable n=3 PF=φ³ icosahedral monotile substitutions. This is
significant: if the result is genuine, it is not a corollary of
known mathematics in the read sample.

**Three structural possibilities are consistent with what we know:**

### (b.1) The result is genuinely novel

The zero-Realised across 170 shape triples could be the first
recorded computational evidence for non-existence of
unmarked 3D icosahedral monotile substitutions in this slice of
the search space. This would be a meaningful negative result —
not a no-go theorem yet, but a substantive empirical input.

**Confidence in this reading: low to moderate.** The slice of
the search space is narrow (max_entry=2, n=3, PF=φ³). The
zero-Realised might evaporate with broader search bounds. Until
broader search confirms or refutes, this is hypothesis, not
result.

**Path to confirmation:** extend max_entry to 3, then 4, ...
The Q10 relay deferred this for compute reasons. If the result
is structural, the same outcome should hold at higher
max_entry. If structural, this is a publishable negative
result.

### (b.2) The result is an artefact of pipeline gaps

Three documented gaps in Apeiron's search that could in principle
mask realisable candidates:

- **Taxonomy gap.** 9 H₃ tetrahedron similarity classes from the
  Paolini ABCK pool, vs Frettlöh's reported 15. 6 candidate
  matrices in the n=3 survey have left-eigenvector volume ratios
  that do not match any of the 9 known classes — the search never
  attempts realisation for them.
- **k_max gap.** k_max=5 in the realisation backtracker. 2 of 13
  candidates at max_entry=2 require k > 5 (sum of matrix entries
  exceeds 5). They are skipped, not refuted.
- **Realisation pipeline gap.** Apeiron's `realise()` checks
  matchings of given prototile shapes against a given matrix. It
  does not yet *derive* shape candidates from the eigenvector
  ratio. The current pipeline assumes the answer to "what shapes"
  is in a finite list it knows about (the Paolini-derived 9-class
  taxonomy), which is the same list as gap 1.

**Confidence in this reading: moderate.** Each gap is documented
and addressable. The Q10c discussion noted that closing all
three is bounded compute, not new mathematics.

**Path to confirmation:** close the gaps and re-run. If after
closure some candidates Realise, the zero-Realised was an
artefact. If after closure the result remains zero, (b.1) is
strengthened.

### (b.3) The result is signal that we are searching the wrong space

Walton-Whittaker 2019 is a marked monotile. The hat/spectre
papers (2023, not in read list but cited in CLAUDE.md §9) are
unmarked monotiles. The hat-spectre route works in 2D. In 3D no
analogous unmarked monotile is known.

If 3D unmarked icosahedral monotile substitutions genuinely do
not exist, but 3D *marked* monotiles can carry the substitution
via decoration (à la Walton-Whittaker generalised to 3D), then
the entire Track B search has been on the wrong target. The right
target would be: marked icosahedral monotiles.

**Confidence in this reading: low.** The 2D unmarked monotile
exists (hat). If a 2D shape-alone monotile exists, no a priori
obstruction prevents a 3D one. But the 2D case took until 2023 to
exhibit; 3D could plausibly require decoration, and the absence
of any 3D unmarked monotile in the literature for 50+ years is at
least suggestive.

**Path to confirmation:** generalise Apeiron's `realise()` to
search for matched (decorated) candidates. This requires
modelling decoration as part of the prototile shape. Significant
infrastructure work.

### Bottom line of (b)

The exhaustive NoRealisation is **not explained** by the read
literature. (b.1), (b.2), (b.3) are all consistent with what we
know. The natural next step is to discriminate between them,
because each implies a different next move.

---

## (c) What new candidate families or formalisms emerge?

The Q10c spec called for one of three deliverable types: new
candidate matrices, infrastructure changes for a different
formalism, or a `strategy_pivot.md` documenting a known
obstruction. The literature read produced none of the three
directly, but it does generate three candidate next-moves keyed
to (b.1)-(b.3) above.

### (c.1) Close the documented Q9c gaps before pivoting

This corresponds to (b.2). Three discrete tasks:

1. **Extend the H₃ tetrahedron taxonomy to 15 classes.** Frettlöh's
   reported count is 15; Apeiron has 9 from the Paolini ABCK pool.
   The missing 6 require a richer vertex pool. Action: grep the
   Frettlöh papers (esp. the Bielefeld Tilings Encyclopedia) for
   the explicit 15-class list and import it.

2. **Lift k_max from 5 to 7 (or 8).** The 2 SKIPped n=3 candidates
   need k = 6 or 7. Compute cost: per-candidate budget grows roughly
   exponentially in k_max (since the DFS-tree branching factor is
   bounded by orbit size × matching-shape count). Probably affordable
   for the 2 specific candidates; not affordable as a default.

3. **Run n=3 at max_entry=3.** Survey produces 209 candidates
   (28.6 s of enumeration). Feasibility: ≈ 8 hours of realisation
   search at the current rate of 30 s/triple budget × ~17 triples per
   candidate avg × 209 candidates. Acceptable but not cheap.

**Cost.** Bounded engineering, no new mathematics. Output is a
fresh exhaustive-sweep result at the closed gaps. If still zero
Realised, (b.1) hypothesis strengthens; if some Realise, (b.2)
hypothesis is the answer.

**Risk.** None to the project; some compute time.

### (c.2) Pivot to marked icosahedral monotile candidates

This corresponds to (b.3). The infrastructure change is
substantial:

- Add a **decoration model** to `Polyhedron`: a finite set of
  edge / face / vertex decorations from a discrete alphabet (e.g.
  charges as in Walton-Whittaker, or color-matchings as in classic
  Wang/Berger constructions).
- Extend `realise()` to **search over decoration assignments** as
  well as shape matchings. The combinatorial cost is multiplicative
  in the alphabet size.
- Generalise `face_match_consistency` and the per-step DFS pruning
  to handle decoration constraints.

**Cost.** Significant infrastructure work, weeks-of-development
scope. Branching factor in the search increases.

**Risk.** Moderate. Decoration is a strict superset of unmarked,
so if the unmarked search has artefactually missed candidates,
decoration will too. But decoration opens a search space we know
contains existing solutions (Walton-Whittaker analogues).

### (c.3) Audit Apeiron's pipeline against Goodman-Strauss

This is a `strategy_pivot.md`-flavoured deliverable but smaller
in scope than (c.2). Two specific audits derived from the
Goodman-Strauss read:

- **Sibling-edge-to-edge audit.** Verify Apeiron's
  Frettlöh-derived Danzer ABCK rule satisfies sibling-edge-to-edge
  (not just sibling-vertex-to-vertex). If not, the four-pillar
  pipeline is operating on a substitution outside Goodman-Strauss's
  enforcement scope, and Track A's pillar 2 (recognisability)
  may be on shakier ground than assumed.
- **Recognisability radius vs κ.** Goodman-Strauss's κ can be
  arbitrarily large (Sadun's pinwheel). Apeiron's `is_recognisable`
  takes a fixed radius (currently 2). If our search's radius bound
  is below the true κ for some candidate, we are implicitly
  rejecting realisable rules. Action: characterise the
  relationship between our radius and Goodman-Strauss's κ;
  potentially make the radius adaptive.

**Cost.** Bounded; ≈ 1-2 days of code review and test additions.
No new infrastructure.

**Risk.** Low. Either an audit passes (small good news) or
identifies a real gap (which then informs further work).

### Recommended next move

The three options are not mutually exclusive. The recommended
sequence is:

1. **(c.3) audit first** — bounded scope, immediate. Either
   confirms the pipeline is sound or identifies a fixable gap.
   Outcome informs whether (b.1)/(b.2) is more likely.
2. **Then (c.1) gap-closure** — bounded compute. Output: an
   exhaustive sweep at max_entry=3 with k_max=7 and 15-class
   taxonomy. This is the empirical resolution between (b.1) and
   (b.2).
3. **Only if both above leave the result intact: (c.2) pivot to
   marked candidates** — this is the major pivot. The relay-back
   to Claude (web) should authorise it explicitly, not slip into
   it by default.

**Rationale for ordering.** (c.3) is cheapest and most likely to
inform; (c.1) is the empirical-completeness step before any
strategy pivot is justified; (c.2) is the most expensive and
requires architectural buy-in. Acting on (c.2) before (c.3) and
(c.1) risks pivoting on incomplete data — exactly the failure
mode Claude (web) flagged at the end of the Q10 ruling ("the
pattern of this project has been: build infrastructure, find no
survivors, build more infrastructure. The literature read breaks
that loop only if the findings are synthesised into an explicit
direction-change before the next build phase starts.").

The synthesis above is the explicit direction. The relay-back is
the authorisation.

---

## What this synthesis does NOT claim

- **Not claimed:** that no 3D icosahedral monotile substitution
  exists. The result is a slice of the search space, not a proof.
- **Not claimed:** that the Goodman-Strauss theorem implies an
  obstruction. The notes flag two relevant signals (sibling-edge-
  to-edge stricter in 3D; construction does not preserve monotile-
  ness) but neither is a derived obstruction.
- **Not claimed:** that decoration is necessary. Walton-Whittaker
  is a sufficiency demonstration in 2D, not a 3D necessity proof.
- **Not claimed:** that the Frettlöh circular-symmetry framework
  is irrelevant. It is irrelevant for *this question* (icosahedral
  finite-orientation case), not for substitution tilings in
  general.

## What is still unread (deferred per scope-limit)

- Goodman-Strauss 1998 Appendix (pp. 39–45): worked example.
  Deferred unless the relay-back flags need.
- Walton-Whittaker §6.3+: spectral / model-set details. Deferred.
- Frettlöh's "direct product variation" papers: the relay
  identified these as candidates but the agent's MANIFEST notes
  the construction is Frank-Robinson, not Frettlöh. The
  Baake-Frank-Grimm 2022 *Fibonacci direct product variation
  tilings* paper (arXiv 2203.07743) is the canonical DPV
  reference; not yet downloaded. Deferred unless the synthesis
  recommends it.
- Hat/spectre 2023 papers (Smith-Myers-Kaplan-Goodman-Strauss).
  In CLAUDE.md §9 minimum reading list. Deferred from the Q10
  cycle; relevant if (c.2) decoration-pivot becomes active.
- Senechal 1995 *Quasicrystals and Geometry* Chapter 7
  (cross-check on Frettlöh ABCK transcription). Already used in
  Track A; not directly relevant to the Q10 synthesis question.

---

## Re-relay action

Per the Q10 sequencing commitment, the next concrete step is to
relay this synthesis to Claude (web) with the recommended next
move (c.3 → c.1 → c.2). Authorisation for (c.1) gap-closure work
is the natural ask; authorisation for (c.2) decoration-pivot
should come only after (c.3) and (c.1) results are in.

---

## (d) Cut-and-project substitute synthesis (Q14b, post-(c.1) closure, 2026-05-04)

**Context.** Q14b authorised a Kramer-Neri 1984 +
Kramer-Papadopolos-Schlottmann-Zeidler 1994 read for the
cut-and-project literature phase. Both originals are firmly
closed-access (agent confirmed via Unpaywall + OpenAlex; no arXiv
preprint, no faculty PDF, no Wayback capture). Per Option C, three
substitute reads conducted: Al-Siyabi-Koca-Koca 2020 (full),
remaining pages of Frettlöh "Icosahedral tilings in R³" (pp.8–9),
and Tilings Encyclopedia ABCK page (does not exist; Frettlöh's
paper is the canonical Encyclopedia reference).

Notes recorded at [literature_notes.md §7](literature_notes.md).

### Findings against Q14b's five questions

**(a) Cut-and-project framework on icosahedral monotile
candidates.** *No published monotile candidate from D₆ projection
exists.* The published F-type icosahedral tiling family — ABCK,
Socolar-Steinhardt, T^(2F), and their mld relatives — is uniformly
multi-prototile (4 prototiles for ABCK and Socolar-Steinhardt).
The published P-type family (Ammann rhombohedra, two prototiles)
is also multi-prototile. The literature has not produced a
single-prototile icosahedral D₆-projection candidate; whether
this is because no one looked or because the framework
structurally precludes it cannot be settled from the substitutes
alone.

**(b) Windows yielding F-type vs other families.** *Partially
answerable.* F-type windows arise as projections of D₆ vertex-
class orbits — Frettlöh §1.3 gives three explicit windows
(regular dodecahedron edge length 2, great dodecahedron long edge
2τ, and a custom polyhedron with 5-fold stars on each pentagonal
face). The three are merged when the vertex classes are taken as
a single model set (Theorem 1.3). What windows produce P-type
(Ammann), B-type, or non-standard families is not in the
substitutes; Kramer-Papadopolos likely catalogues this.

**(c) Window-deformation literature for monotile candidate
generation.** *Not addressed in any substitute.* Al-Siyabi-Koca-
Koca explicitly use a *non-CPS* approach (direct projection of a
(m₁, m₂)-parameterised D₆ subset). Frettlöh §1.3 fixes the windows
and does not discuss deformation. No substitute touches the
"vary the window to find a monotile" question. This is the most
significant gap: window-deformation as a search method may be
unprecedented in the literature, OR it may be in Kramer-Neri /
Kramer-Papadopolos but unverifiable here.

**(d) Window as continuous parameter vs discrete choice.**
*Unanswered from substitutes.* The (m₁, m₂) integer-pair
parameter in Al-Siyabi-Koca-Koca is *discrete* (with the
constraint m₁ + m₂ even), but it is not a window — it is a
lattice-subset selector. Whether the *window* is a continuous
parameter (allowing infinitesimal deformations) or a discrete
choice (selecting from a finite alphabet of polytope shapes)
determines whether window-search is continuous optimisation vs
combinatorial enumeration. The substitutes do not settle this.

**(e) Relationship between windows and prototile shape.**
*Unanswered from substitutes.* In the substitute reading,
prototile shapes are derived from the H₃ root system (Frettlöh
§1.1) and from D₆ weight orbits (Al-Siyabi-Koca-Koca Table 1) —
not from window geometry. Whether the window directly determines
the prototile shape (in which case window-search is the right
search-space) or merely controls the tiling's vertex set (in
which case the prototile shape is determined separately by H₃
structure) is ambiguous. The substitutes lean toward the latter:
Frettlöh's three windows produce the same tiling-with-shapes
{A, B, C, K}, suggesting the windows pick vertex classes, not
shapes. If true, *window-search would NOT generate new
prototile candidates* — the framework is rigid in shape, flexible
only in vertex-class selection. But this is an inference from
absence of contrary evidence in the substitutes, not a proven
claim.

### Bottom line of (d)

The substitute reading **definitively confirms** the following:

1. ABCK is a D₆ projection (Theorem 1.3, Frettlöh).
2. Inflation by τ corresponds to a Fibonacci action on (m₁, m₂)
   in the lattice-subset framework (Al-Siyabi eq. 16).
3. F-type and P-type icosahedral tiling families are
   uniformly multi-prototile in the published literature.
4. Window-deformation as a candidate-generation strategy has *no
   visible precedent* in the open-access literature on F-type
   icosahedral tilings.

The substitute reading **leaves open** Q14b(c), (d), (e) — the
specific questions Claude (web) flagged as needed for
infrastructure-design decisions. Without those answers, building
a "window search" CSP (the c.2-replacement infrastructure
recommended at Q11c) would proceed on incomplete foundations.

### Three options for closing the cut-and-project gap

**(d.1) Acquire Kramer-Neri 1984 + Kramer-Papadopolos 1994 via
proper channels** (institutional library / interlibrary loan /
direct author contact). Direct path; bounded delay (days to
weeks). Recommended if any infrastructure work depends on
Q14b(c)/(d)/(e).

**(d.2) Read closed-access alternatives in the same family.**
Candidates discovered in the substitutes' references but not yet
inspected: Kramer & Andrle 2004 *J. Phys. A* 37 (Danzer tiles
from wavelet POV); Kasner & Böttger 1993 *IJMP-B* 7 (lattice
dynamics of F-type icosahedral quasicrystal); Roth 1993
*J. Phys. A* 26 (face-centred icosahedral local derivability);
Baake & Grimm 2013 Vol. 1 Chapter 7 (already in CLAUDE.md §7.1
as a working reference; chapter covers model-set / CPS
framework). All likely closed-access except Baake-Grimm, which
is institutional-library-only.

**(d.3) Re-read the substitutes with a different lens.** The
absence of monotile candidates in the published F-type / P-type
families is itself a meaningful empirical observation, possibly
sufficient for a strategy_pivot.md without further reading.
Hypothesis to test: "the icosahedral cut-and-project framework
inherently produces multi-prototile tilings; the search for a 3D
icosahedral monotile must therefore look outside this framework."
If true, the (c.2-replacement) cut-and-project pivot is the
*wrong* pivot, and either the original (c.2) decoration pivot or
a third direction (e.g., non-icosahedral 3D substitution) becomes
the right next move.

### Recommended next step

Re-relay to Claude (web) with the synthesis above. The relay
authorisation is for a decision among (d.1), (d.2), (d.3), or a
combination. Acting on any one without authorisation is the
failure mode the Q10 sequencing commitment was designed to
prevent.

---

## (e) Moody (2000) read result — Gate decision (post-Path A, 2026-05-04)

**Context.** Path A (substitute reading) authorised by user
2026-05-04 in lieu of acquiring closed-access Baake-Grimm Vol. 1
Ch. 7. Moody (2000) read in full as the strongest open-access
substitute for §7.2 of Baake-Grimm. Notes at
[literature_notes.md §8](literature_notes.md).

### Findings against the three Q15b/Q14b gate questions

**(1) Theorem ruling out icosahedral CPS monotiles?** *No.* Moody
never raises the question. Model sets are point sets; the
connection from model sets to tilings (and prototile counts) is
not Moody's subject. The §4.2 p-adic Robinson example
illustrates *one direction* — for tilings whose vertex sets are
model sets, multi-prototile-ness corresponds to window
decomposition into one open subset per tile type — but this is
not a theorem, just an illustration.

**(2) Window-deformation as candidate-generation strategy?**
*Partially.* Translation of the window is continuous and
parameterised by G (internal space group); but it produces only
the LI class of a generic model set, not structurally new
tilings. *Shape-deformation* of the window is **not** discussed
as a deformation type. Moody's framework treats windows as fixed
once chosen, not as continuous parameters of a search space.

**(3) Q14b(d) and (e)?** *Partial.*

- **Q14b(d) — windows continuous vs discrete.** The window is a
  *set* in Moody's framework, not a parameter from a
  parameterised family. Translation is continuous; shape
  variation is not formalised. **Window-shape parameter space is
  not in Moody.**
- **Q14b(e) — window → prototile-shape map.** Sharper answer
  via §4.2: in tilings whose vertex sets are model sets, each
  prototile type corresponds to a distinct open subset of internal
  space. **For ABCK the relation is more complex:** Frettlöh's 3
  vertex-class windows produce 4 prototiles via the H₃
  fundamental-region structure. Window decomposition determines
  *vertex classes*; tile shapes emerge from *combinatorial
  arrangement of vertex classes* under H₃.

### Refined structural multi-prototile conjecture

The original conjecture is sharpened by Moody's framework into:

> **Refined conjecture (post-Moody).** For every D₆ cut-and-
> project tiling of ℝ³ with icosahedral symmetry where the
> prototile-shape structure is determined by the H₃ tetrahedral
> fundamental region (the ABCK / icosian-model-set family per
> Moody §4.1), the tiling requires more than one prototile type.

Three escape routes that would refute it:

1. A D₆ vertex-class arrangement that yields a single tile type
   (counterexample to the conjecture as stated).
2. A non-H₃-fundamental-region prototile-shape structure on a D₆
   model set (would refute the *scope* by exhibiting a different
   shape family).
3. The D₆ weight lattice or an intermediate D₆ lattice yielding
   a monotile (would refute by exhibiting a sibling case).

Moody §4.1 explicitly notes (verbatim, p.9): "the only other
relevant lattices in 6-space are the D_6 weight lattice and the
lattices lying between the root and weight lattices." So the
conjecture's scope is well-defined: D₆ root + weight +
intermediates exhaust the icosahedral CPS framework per Moody.

### Gate resolution

Per [strategy_pivot.md §5](strategy_pivot.md):

- **Gate A** (Branch A, conjecture false) requires Baake-Grimm to
  contain a counterexample, a window-deformation theorem, or an
  explicit "monotile windows admitted" remark. Moody contains
  none of these. **Gate A NOT activated by this read.**
- **Gate B** (Branch B, conjecture true) requires Baake-Grimm to
  contain a no-go theorem, a minimum-prototile-count theorem, or
  an explicit no-go remark. Moody contains *partial inferential
  support* (the window-decomposition correspondence + "no
  monotile candidate appears anywhere in the published F-type
  family") but **no theorem proving the conjecture.** Gate B
  partially activated; not fully resolved.
- **Gate C** (Branch C, undecidable) requires Baake-Grimm to
  contain discussion without resolution, OR silence on the
  structural question, OR a reference to Kramer-Papadopolos as
  containing the resolution. Moody is **largely silent on the
  structural question**. Gate C is the closest match for what
  Moody returned: he refines the framework and gives partial
  evidence consistent with the conjecture, but does not settle
  it.

### Verdict: Gate C activated

The Moody read does **not** resolve the structural conjecture.
What it does:

1. **Sharpens** the conjecture to a more falsifiable form (the
   refined version above) by clarifying the window → prototile-
   shape relationship.
2. **Bounds** the framework's scope (D₆ root + weight +
   intermediates exhaust icosahedral CPS per Moody §4.1).
3. **Confirms** that no shape-deformation literature exists in
   the open-access substitutes; if it exists, it's in
   Kramer-Papadopolos 1994 (closed-access).

### Recommended next step (post-Moody)

Per Q15b's authorised sequence, Gate C activates the **(d.1)
acquire-the-originals** option. The recommended next move:

- **Pursue institutional access** to Kramer-Papadopolos-
  Schlottmann-Zeidler (1994) *Projection of the Danzer Tiling*,
  *J. Phys. A* 27, 4505–4517. The DOI is
  10.1088/0305-4470/27/13/024 (per the agent's correction of an
  earlier typo). Inter-library loan via any UK or US R1 maths
  library. Allow days to weeks.
- **In parallel**: pursue *one* of the three escape routes
  computationally, since they are concrete and Apeiron's existing
  D₆ infrastructure could be adapted. Specifically Route 1 (test
  whether any vertex-class arrangement under H₃ yields a single
  tile type) is a finite combinatorial check, ~weeks of
  engineering on top of existing Apeiron primitives.

The user's authorisation is the gate for either of these. Until
authorised, no infrastructure work or paper-acquisition activity
lands. The conjecture stands as "open, refined, and bounded."
