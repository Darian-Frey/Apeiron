# Literature notes — running notes during the Q10 read

**Purpose.** Per the Q10 ruling (Output 1), this file collects raw
notes during the literature deep-dive triggered by the exhaustive
n=3 PF=φ³ NoRealisation result.

**Format.** Each entry: (paper, §section, claim, relevance to Apeiron).
**No synthesis here.** Synthesis goes in `literature_review.md`
(Output 2) only after the reading is complete.

**Stop condition.** Goodman-Strauss 1998 and the two Frettlöh targets
(direct product variation; Parallelogram tilings, worms and finite
orientations 2013) fully read. Time-box: one week of reading.

**Reading priority** (per Q10b):

1. **Goodman-Strauss, C. (1998).** *Matching rules and substitution
   tilings.* Annals of Mathematics, 147(1), 181–223. — IN PROGRESS.
2. **Frettlöh, D.** *Substitution tilings with statistical circular
   symmetry* (2008, EJC + 2008 companion).
3. **Frettlöh, D., & Harriss, E. (2013).** *Parallelogram tilings,
   worms and finite orientations.*
4. **Walton & Whittaker (2019)** [substituted for "Socolar 2019",
   which does not exist; see MANIFEST].

---

## 1. Goodman-Strauss (1998), *Matching rules and substitution tilings*

**Note on section numbering.** The Q10 relay referenced "§3 + §6"
and "Theorem 6.1." The paper has only §1–§5 + Appendix; the
relay's labels are mislabeled. The substantive equivalents are
§3 (constructing matched prototiles) and §4 (proof of enforcement);
the **Main Theorem** is stated on p.2, no internal numbering.

### §0 / abstract / introduction (pp.1–3)

**Main Theorem (p.2).** Every substitution tiling of E^d, d > 1,
can be enforced with finite matching rules, subject to a mild
condition: the tiles are required to admit a set of *hereditary
edges* such that the substitution tiling is *sibling-edge-to-edge*.

> Relevance: This is the headline result. It is a **sufficiency**
> claim: substitution → matching rules. It says nothing about which
> substitutions exist. The active question (does the theorem imply
> an obstruction to 3D icosahedral monotile substitution?) is
> answered structurally by reading the construction itself in §3.

**Sibling-edge-to-edge condition (p.3).** "the substitution breaks
the edges of the parent into edges of the children and that
siblings' edges coincide if they overlap at all." Stated as
"mild" but not always satisfied — author acknowledges (p.10)
example may not be hard to find for d > 2.

> Relevance: 3D-specific obstruction signal. In d=2 hereditary
> vertices ⇔ hereditary edges (Lemma 1.3); in d > 2 the converse
> is "probably not true." Sibling-edge-to-edge is more restrictive
> in 3D than in 2D.

**Sharp at d = 1 (p.3).** "no 1-dimensional substitution tiling can
be enforced by matching rules, and our construction completely
breaks down when d = 1." Theorem requires d > 1 strictly.

**Two general aperiodic-tile methods cited (p.3).** (1) matching
rules forcing a substitution tiling, as in this theorem; (2)
matching rules forcing a quasiperiodic tiling (cut-and-project).
Le T.T.Q. has a similar theorem for certain quasiperiodic classes.

**Other forced aperiodic tilings (p.3).** "Exactly two other
classes of forced aperiodic tilings are known. The first was found
by P. Schmitt and altered by J.H. Conway and L. Danzer [SCD]. These
tilings are of E^3 and are not isotropic. A new class emerged in
1996: J. Kari tiled the plane with tiles reminiscent of Wang's
constructions of the early 1960s."

> Relevance: as of 1998, only three known classes — substitution
> tilings (this theorem), quasiperiodic, and SCD/Kari. The
> substitution route is the only one we are pursuing for 3D.

### §1 Definitions

**§1.1 Matching rules (p.4).** Restrictions on which tiles may be
adjacent. "The simplest form of these restrictions is simply to
require the tiles fit together — this is sufficient for our needs."

**§1.2 Inflation σ (p.5).** Inflation = expanding linear map; for
each child C_i, for every n ∈ ℕ there exists C_i^(n) ∈ G such that
σ^n C_i = C_i^(n) σ. (This is the recursion-compatibility
condition.)

**§1.2 Substitution-as-self-similar (p.6).** "Note the requirement
that an inflated copy of a prototile is congruent to the union of
its offspring is very strong indeed, and at first glance seems
not to be obeyed in some well known examples, such as the Penrose
rhombs, or fractal examples such as Thurston-Kenyon. However,
these examples can all be recomposed into a tiling with
substitution satisfying our requirements."

> Relevance: the framework demands strict self-similarity of
> substitution. Apeiron's substitution.py / hierarchy.py operate
> in the same setting.

**§1.2 Connected hierarchy (pp.7–8).** A tiling has *connected
hierarchy* iff every two points lie in a common supertile. Three
properties: (i) every bounded configuration is congruent to a
configuration in a connected-hierarchy tiling; (ii) in any
G-invariant probability measure, connected-hierarchy tilings have
measure 1; (iii) covered by infinite-level supertiles.

> Relevance: this is the measure-theoretic backbone of the
> "almost every tiling is hierarchical" argument. Apeiron's
> Pillar 4 (no non-hierarchical tilings) is the *strict*
> version — measure 1 is not measure-everything.

**§1.4 Hereditary edges (pp.9–10).** A set of edges is *hereditary*
if for every prototile A and every edge e of A, σ(e) is exactly
tiled by edges of children of A, and every k-facet of σ(e) is
tiled by k-facets of children. *Hereditary vertices* defined
analogously.

**§1.4 Sibling edge-to-edge / sibling vertex-to-vertex (p.10).**
Two definitions; for d=2, edge ⇒ vertex (Lemma 1.3). For d > 2,
vertex-to-vertex does **not** imply edge-to-edge ("probably not
true when d > 2").

> Relevance: **Direct 3D obstruction signal.** Apeiron's Track A
> (Danzer ABCK + Paolini real geometry) needs to be checked for
> sibling-edge-to-edge separately from sibling-vertex-to-vertex.
> The Frettlöh ABCK rule is sibling-vertex-to-vertex; whether it
> is sibling-edge-to-edge in 3D is non-obvious and unverified
> here. ACTION ITEM for review of `apeiron/substitution.py` and
> the Paolini dissection: confirm the Danzer rule is sibling-
> edge-to-edge, not just sibling-vertex-to-vertex.

**§1.4 Mild conditions take 1 (p.10).** "for the time being it is
unknown if the conditions are always satisfied" — i.e., the
hereditary-edge / sibling-edge-to-edge condition is not known
to be necessary, only sufficient (and an example failing it for
d > 2 is anticipated but not yet exhibited as of 1998).

**§1.5 Endo-, Meso-, Epi-vertices (pp.11–12).** Endovertex:
incident to interior edges of σ^n(A) for all n ≥ κ(v); also
incident to E(A^-) in σ(A^-). Mesovertex: incident to E(A^-) in
σ(A^-) but no σ^k(v) is incident to interior edges. Epivertex:
not incident to E(A^-) in σ(A^-).

> Relevance: terminology used in Apeiron's hierarchy.py for
> position_signature and the recognisability framework. No
> structural surprise.

### §2 Selecting structures (skeletons, keys, vertex hulls)

Dense technical construction. Notes elided — the substantive
content is the parameter κ.

**§2.1.3 Resolution κ (pp.14–15).** "κ is the resolution at which
we view the hierarchical structure." For each B ∈ S choose a
natural κ(B) such that there is a connected collection of edges
in σ^κ(B)(B) covering all relevant edges/sites/endovertices.

**§2.1.3 κ in practice (p.15).** "κ is often rather low; in fact,
many well known examples have κ = 1 and very few have κ > 2." The
Conway-Radin pinwheel and Robinson triangles have κ = 1; the
L-tiling requires κ = 2. **One can construct examples requiring
arbitrarily large κ, such as Sadun's generalized pinwheel
tilings.**

> Relevance: κ is a free parameter of the construction; for
> Apeiron's purposes, the recognisability radius corresponds
> roughly to κ. If we are searching at radius 2 and miss
> candidates at radius 5 or higher, we may be missing realisable
> rules. **ACTION ITEM:** check whether Apeiron's recognisability
> radius bound corresponds to κ, and whether we have implicitly
> set it too low.

### §3 Creating tiles and markings

**§3.1 Three flavors of prototile (p.27).** The construction
produces three types of marked prototile: **vertex tiles** (disks
or sectors of d-balls, marked at vertices), **edge tiles**
(δ-thickened images of edges), **small tiles** (the original
prototile minus vertex and edge regions, marked with packets
encoding the supertile address).

> Relevance: **CRITICAL FINDING for the active question.** Even
> if the original substitution is on a single prototile (n = 1),
> the construction expands it into a *finite collection of
> marked prototiles* T'. The theorem proves substitution tilings
> are enforceable by matching rules on a *set* of marked tiles,
> not on a single shape. The construction does NOT preserve
> monotile-ness.

**§3.5 Matching rules M (pp.31–32).** Three explicit rules:

1. **Terminal matching rule:** terminator marked R = X_κ...X_2 X_1
   incident to a vertex tile with terminal vertex marking v''
   such that v ∈ σ(X_2), v'' aligned correctly.
2. **Edge matching rule:** end of edge tile [e, Q_1, ..., Q_k]
   meets vertex tile with vertex marking [v, Q_1, ..., Q_k] with
   matching orientations and labels.
3. **Tiling rule:** tiles cover the plane with disjoint interiors.

Plus implicit: well-formed-packet rules (§2.4), compatibility
rules for vertex-marking stacks (§3.3), big-tile rule (§3.7).

**§3.6 Labelling (p.32).** "Note every supertile can only be
marked in finitely many ways, since each supertile can only be
marked in finitely many places and there are only finitely many
possible markings."

### §4 The proof of the Theorem

**§4 Proof structure (pp.33–37).** Inductive: any almost-well-
formed supertile of level n lies in an almost-well-formed
supertile of level n+1 (Proposition 4.1). Plus Lemma 4.2:
almost-well-formed of level n inside almost-well-formed of level
> n + κ ⇒ well-formed of level n.

**§4 Crucial mechanism (p.37).** "We used the condition in the
statement of the theorem very strongly: we needed a mechanism —
**vertex wires** — to fix the position of the vertices of σ(A) and
thus keep edges from propagating beyond their intended borders.
The vertex wires required vertices to be hereditary. Secondly, to
fix the position of sibling supertiles relative to some initial
supertile, we need some point in each we can say they share —
they are forced to share vertices."

> Relevance: the construction is built on **vertex wires**, which
> require hereditary vertices. The "we need some point in each
> we can say they share" is the locking mechanism. In a 3D
> icosahedral monotile substitution, whether such shared vertices
> exist between siblings is a geometric question; this is where
> the d > 2 sibling-edge-to-edge gap from §1.4 becomes load-
> bearing.

**§4 Open question raised by author (p.37).** "Other mechanisms
can be devised, exploiting other conditions one could make.
However it is not clear if an example exists in which the
relative positions of sibling tiles cannot be fixed at all."

> Relevance: **Open question still relevant in 1998 and unresolved
> in this paper:** are there substitution tilings whose sibling
> positions are genuinely not lockable? If yes, the matching-rule
> enforcement fails for them.

**§4 Lemma 4.3 (aperiodicity, p.37).** Every matching rule tiling
constructed in the proof is aperiodic. Proof: every well-formed
supertile is identifiable; an infinite cyclic isometry fixing the
tiling would have to fix some supertile larger than its translation
distance — contradiction.

**§4 Corollary 4.4 (p.37).** There are infinitely many aperiodic,
hierarchical matching rule tilings.

> Relevance: this is the main aperiodicity payoff. Note the result
> is about *matching rule tilings*, i.e. the *constructed* tilings
> in (M, T'), not the original substitution tilings (T, σ, S).
> The original tilings inherit aperiodicity from the matching
> rule tilings via the one-to-one correspondence (mod measure
> zero, p.37).

### §5 Worked example (pp.38+)

A 2D variation of one of Danzer's tilings (T = {X, Y}, S = {A, B,
C, D}, σ^6 illustrated). κ = 1. Note: this is Danzer's *2D*
tiling, not the 3D ABCK. Apeiron's Danzer encoding is the 3D
ABCK, a different object.

> Relevance: confirms the construction is illustrated only in 2D
> in this paper; the 3D case is theoretically covered by the
> Main Theorem (d > 1) but never worked out as an explicit
> example. Reading further into §5 (Appendix) only useful if we
> need a concrete decoration template; the algebraic obstruction
> question doesn't depend on it.

### Bottom-line answers to the active question (provisional, pre-synthesis)

**Active question:** "Does Theorem 6.1 (or its corollaries) imply
any obstruction to 3D icosahedral monotile substitution, or does
it only establish sufficiency conditions?"

**Provisional answer (NOT to be treated as synthesis until
literature_review.md):**

1. The theorem is **sufficiency-only**: substitution → matching
   rules. It does not assert any substitution exists for any given
   tile alphabet.
2. The construction **does not preserve monotile-ness**: starting
   from n=1, it produces a finite *set* of marked prototiles
   (vertex tiles, edge tiles, small tiles, big tiles). The theorem
   says nothing about whether the original tile *shape* can carry
   the substitution as a single unmarked monotile.
3. The "mild condition" (sibling edge-to-edge with hereditary
   edges) is **stricter in 3D than 2D** (Lemma 1.3 plus author's
   note p.10). This is a real, paper-flagged 3D-specific
   obstruction: in 3D, vertex-to-vertex does not automatically
   give edge-to-edge.
4. The construction's locking mechanism (**vertex wires**, p.37)
   requires hereditary vertices and shared sibling vertices. The
   author leaves open whether substitutions exist where sibling
   positions cannot be fixed at all (p.37).

These do **not** explain Apeiron's exhaustive zero-Realised
result directly, because Apeiron is checking whether *any*
substitution rule on n=3 prototiles can be geometrically realised
— a question prior to and independent of matching-rule
enforcement. The Goodman-Strauss theorem says: if you have a
realisation, the matching rules exist. Apeiron's question: do you
have a realisation? Goodman-Strauss is silent on that.

**Suggested follow-ups beyond the immediate read:**

- Verify Apeiron's Frettlöh-derived Danzer rule satisfies sibling-
  edge-to-edge (not just vertex-to-vertex). If not, the four-
  pillar pipeline may be working with a substitution outside
  Goodman-Strauss's enforcement scope.
- Check whether Apeiron's recognisability radius corresponds to
  Goodman-Strauss's κ; if our search bound caps κ implicitly at
  a value ≤ Sadun's pinwheel example, we may be missing high-κ
  realisable rules.
- Read §4.1.2 of the Frettlöh papers next for orientation-class
  constraints; see whether a 3D version of the parallelogram-
  tiling argument bounds the n in n×n substitution matrices.

---

## 2. Frettlöh & Harriss (2013), *Parallelogram tilings, worms and finite orientations*

### Setting (§1, pp.1–2)

**Object.** Locally finite vertex-to-vertex tilings of ℝ² by
parallelograms with finite protoset.

**Definition 1.1 (worm).** Given a tile T and edge e of T, the
unique adjacent tile T_1 sharing e has a parallel edge e_1 ≠ e;
T_2 shares e_1; and so on. The biinfinite sequence W = ...T_{-2},
T_{-1}, T, T_1, T_2,... is the *worm given by T and e*.

> Relevance: a worm is the parallelogram-tiling analogue of a
> de Bruijn line in canonical projection tilings. The whole
> structural argument relies on parallel-edge propagation — a
> property that requires the prototiles to be parallelograms /
> parallelotopes. Tetrahedra (Apeiron's Track A/B candidates) are
> NOT parallelotopes, so the theorem does NOT directly apply.

### §2 Three lemmas

**Lemma 2.1 (Cone Lemma, p.3).** A worm cannot bend "too much": it
is contained in the complement of two open cones C_1, C_2 of half-
angle α (the minimum interior angle of any prototile).

**Lemma 2.2 (No-loop, p.3).** A worm has no loop.

**Lemma 2.3 (Crossing, p.3).** Two worms cross at most once. Worms
in the same family (parallel defining edges) never cross.

**Lemma 2.4 (Travel, p.4).** Any two tiles can be connected by a
finite sequence of tiles contained in at most k = ⌈2π/α⌉ different
worms.

### §2 Main theorem

**Theorem 2.5 (p.6).** *Each prototile occurs in a finite number of
orientations in any parallelogram tiling with finite protoset.*
Bound: N := (1/(2m-1)) · ((2m)^{k+1} + 2m − 2), where m = number
of prototiles, k = ⌈2π/α⌉.

> Relevance: the bound N is exponential in k, polynomial in m. For
> α small (very thin parallelograms) k blows up, but for any
> *finite* m and α > 0 the orientation count is finite.

### §3 Remarks (p.7)

**Sharpness (Figure 6).** Theorem 2.5 fails if infinite prototiles
are allowed: a square tile can occur in infinitely many
orientations if we allow infinitely many distinct rectangle
prototiles around it.

**3D extension (p.7, last paragraph).** "Equivalent results to
those proved in this paper should apply to tilings of ℝ^d (d ≥ 2)
by parallelotopes. In order to generalise Lemma 2.4, however, we
need to enclose some tile T with a finite collection of tiles.
This would require further objects... pseudo-planes (and pseudo-
hyperplanes) made up of linked parallelotopes sharing a parallel
edge in common."

> Relevance: **Author explicitly notes the 3D extension is
> incomplete.** Lemma 2.4's proof ("worms encircle T to bound the
> turn count") doesn't transfer directly to ℝ³ — you cannot encircle
> a tile with a 1-parameter family of worms in 3D; you need
> 2-parameter pseudo-hyperplanes. The 15-class H₃ tetrahedron
> taxonomy completion (Apeiron's Q10b secondary gap) is **not**
> resolved by this paper, even in spirit. The taxonomy is for
> tetrahedra under icosahedral symmetry, not parallelotopes.
> **Bottom line: this paper does NOT inform the 15-class
> completion.** Apeiron's deferred follow-up (extend H₃ taxonomy
> to all 15) needs a different reference — likely Frettlöh's older
> ABCK papers or the Tilings Encyclopedia entry for Danzer.

---

## 3. Frettlöh (2008, EJC), *Substitution tilings with statistical circular symmetry*

### Setting (§§2–3, pp.2–4)

**Substitution σ in ℝ²** with prototiles T_1,...,T_m. Substitution
matrix S_σ = (|Φ_{ij}|). Self-similar: λT_i = ⋃_{T ∈ σ(T_i)} T.

**Definition 3.1.** A substitution tiling is *pinwheel-like* if
infinitely many distinct values α(T) occur (i.e., tiles appear in
infinitely many orientations).

**Definition 3.2 (statistical circular symmetry).** Tiles are
ordered by inclusion in supertiles; the orientation distribution
α(T_j) is uniformly distributed in [0, 2π).

### §3 Proposition 3.4 (p.5)

**Statement.** σ primitive in ℝ² with prototiles T_1,...,T_m. Then
σ is pinwheel-like ⇔ there exist n, i such that σ^n(T_i) contains
two tiles T, T' of the same type with **α(T) − α(T') ∉ πℚ**.

> Relevance: **direct test for whether Apeiron's icosahedral
> substitution candidates have finite or infinite orientation
> classes.** Apeiron's candidate substitution rotations live in the
> icosahedral group I (order 60) or I_h (order 120). For any two
> elements g, h ∈ I, the rotation angles satisfy α(g) − α(h) ∈
> {kπ/30 : k ∈ ℤ}, since I has rotation orders 1, 2, 3, 5 (eigen-
> angles 0, π, ±2π/3, ±2π/5, ±4π/5 in ℝ², and analogues from
> 3D rotation projection). All these are rational multiples of π.
> **Therefore: by Proposition 3.4, no icosahedral substitution
> tiling is pinwheel-like — orientations are finite.** The
> Frettlöh 2008 paper does NOT identify an obstruction in our
> setting; it identifies an obstruction for the *complementary*
> case (pinwheel-like, infinite orientations).

### §6 Main theorems

**Theorem 6.1 (p.8).** Pinwheel-like ⇒ statistical circular
symmetry.

**Theorem 6.3 (p.10).** Statistical circular symmetry ⇒ uniform
patch frequencies.

> Relevance: These are properties of the *complementary* case
> (infinite orientations). Apeiron's icosahedral case is the
> finite-orientation case; uniform patch frequencies hold there
> too via the standard Perron-Frobenius / minimality argument
> (Theorem 4.1 of the companion paper §4 below). No new constraint
> for our setting.

### §7 Remarks (p.11): PV substitution factors

"if the substitution factor is a PV number, then the tiling is
FLC under a certain (mild) condition." PV = Pisot-Vijayaraghavan.

> Relevance: φ² (Apeiron's canonical inflation) is a PV number
> (the smallest quadratic PV ≠ φ; note φ itself is also PV).
> So Apeiron's substitutions are FLC under the cited mild
> condition. No surprise; matches CLAUDE.md §3.1 motivation.

### Synthesis of Frettlöh 2008 EJC for the active question

> **Conclusion: this paper does NOT address Apeiron's setting.** It
> classifies substitution tilings whose orientations are infinite
> (pinwheel-like) and proves they have circular-symmetric
> diffraction. Icosahedral monotile substitutions have finite
> orientation classes (the icosahedral group is order 60), so
> Proposition 3.4's condition fails. The paper's machinery is
> orthogonal to our zero-Realised result.

---

## 4. Frettlöh (2008 companion), *About substitution tilings with statistical circular symmetry*

### Highlights (pp.1–6)

- Definition 1.1 (substitution tiling): every finite subset of 𝒯
  is congruent to a subset of some supertile σ^k(T_i).
- Theorem 2.2 (p.2): primitive substitution tiling with prototile
  in infinitely many orientations ⇒ statistical circular symmetry.
- Theorem 2.3 (p.3): in primitive substitution tilings, every
  prototile occurs in each of its orientations with the same
  frequency. **The frequency vector is the normed Perron-
  Frobenius eigenvector of the substitution matrix M_σ.**
- Theorems 3.1, 3.3 (pp.5): autocorrelation and diffraction are
  circular-symmetric in the statistical-circular-symmetry case.
- Theorem 4.1 (p.6): for primitive substitution tilings, ℝ²-action
  uniquely ergodic if finite orientations; E(2)-action uniquely
  ergodic if infinite orientations.

> Relevance: **Theorem 2.3's eigenvector statement is the basis for
> Apeiron's volume/frequency vector computation in
> `pf_eigenvector_in_zphi`.** The right eigenvector of the
> substitution matrix M_σ gives prototile *frequencies* in any
> primitive substitution tiling. (Apeiron's earlier bug: confusing
> right eigenvector with left eigenvector for *volume* matching.
> The companion paper does not explicitly disambiguate; the volume
> identity volume_i = (left eigenvector)_i comes from the linear-
> algebra identity λ · (left eigenvector M) = left eigenvector,
> NOT from this paper. Bug fix `a82715a` is in agreement with the
> companion paper insofar as the frequency vector is right-PF;
> the volume vector is left-PF, which is implicit but not stated
> here.)
>
> Otherwise: companion is purely a 2D-circular-symmetry paper;
> not directly relevant to the 3D zero-Realised question.

---

## 5. Walton & Whittaker (2019), *An aperiodic tile with edge-to-edge orientational matching rules*

### Setting (§1, pp.1–4)

**Tile.** A single hexagonal shape (Figure 1) with edge
decorations: black lines (R1) and charge symbols ±, oriented
clockwise/anticlockwise (R2).

**Matching rules.** Two tiles t_1, t_2 may meet along edge e iff:

- **R1:** black-line decorations are continuous across e.
- **R2:** if both charges at e are oriented clockwise, they must
  be opposite in sign.

**Theorem 1.1 (p.3).** Valid tilings exist; every valid tiling is
non-periodic.

> Relevance: This is a **2D marked monotile** — same shape but
> with edge decorations forcing aperiodicity. The decoration is
> non-trivial (uses orientation-dependent charges, not just colour
> matching). It does *not* prove a shape-only monotile exists;
> the rules R1, R2 cannot be replaced by tile geometry alone.

### §2 Aperiodicity argument (pp.4–8)

The R1 rule forces R1-triangles. Lemma 2.2 shows R1-edges have
unbounded length (no finite max R1-triangle). Theorem 2.4: any
valid tiling is non-periodic.

> Relevance: A "shape forces hierarchy" argument with charge
> propagation (R2) doing the locking. The locking mechanism is
> not the substitution itself but the charge-flip property along
> edge propagation. Mechanically distinct from
> Goodman-Strauss's vertex-wire mechanism (substitution-derived
> labels propagating along skeletons).

### §5 Substitution (p.18)

**Proposition 5.4.** Almost every valid tiling has a unique
locally-defined level-n supertiling for each n. Defects: P_∞
tilings (3-fold rotational symmetry, central infinite cycle),
n-cycle tilings (defect at level n), fault-line tilings.

> Relevance: confirms the matching-rule tiling is recognisable
> via substitution, with measure-zero defects. Same structural
> picture as Goodman-Strauss's enforced substitution tilings —
> connected hierarchy has measure 1.

### §6 Dynamical and spectral properties (pp.20–23)

- Corollary 6.3 (p.21): hull (Ω, ℝ²) is uniquely ergodic.
- §6.2: hull's maximal equicontinuous factor (MEF) is the 2D
  dyadic solenoid; pure point dynamical spectrum (modulo a
  charge-flip 2-fold cover); regular model set structure.

> Relevance: this gives an example of a marked monotile whose
> dynamics are essentially crystallographic (pure point spectrum,
> regular model set), distinguishing it from Penrose
> ((1+ε+ε²)-tilings) and Socolar-Taylor.

### Synthesis of Walton-Whittaker 2019 for the active question

> **Conclusion: this paper is meta-relevant, not directly
> applicable.** The active question is about 3D *unmarked*
> monotile substitutions. Walton-Whittaker is a 2D *marked*
> monotile. Its existence is a reminder that decoration can carry
> aperiodicity when shape alone does not. If Apeiron's exhaustive
> n=3 NoRealisation result is structural (no unmarked icosahedral
> 3D monotile substitution can be realised), the next-broader
> search would be marked monotiles — adding decorations to the
> tile and re-running the search. This is a strategy-pivot signal
> rather than an obstruction proof.

---

## Provisional cross-paper synthesis (PRE-Output-2; not final)

Per the Q10 ruling, the synthesis goes in `literature_review.md`
(Output 2). The notes below are PROVISIONAL and not authoritative.

**On the active question (does the literature imply an obstruction
to 3D icosahedral monotile substitution?):**

1. **Goodman-Strauss 1998:** sufficiency-only. No obstruction.
   Two relevant signals — sibling-edge-to-edge condition is
   stricter in 3D (Lemma 1.3); construction does not preserve
   monotile-ness (§3.1 produces multiple marked tile types from a
   single input prototile).
2. **Frettlöh-Harriss 2013:** does not apply (parallelogram-only;
   tetrahedra not covered). Author flags 3D parallelotope
   extension as incomplete.
3. **Frettlöh 2008 EJC + companion:** does not apply (icosahedral
   tilings are not pinwheel-like; orientations are finite).
4. **Walton-Whittaker 2019:** strategy-pivot signal. A *marked*
   monotile carries 2D aperiodicity when shape-only does not.
   The 3D analogue would be a marked icosahedral monotile,
   broadening the search beyond unmarked candidates.

**No direct obstruction to 3D icosahedral monotile substitution
appears in the read literature.** The zero-Realised result must be
explained by either:

- A finer obstruction in the 3D geometric realisation step
  (something Apeiron-specific to the H₃ taxonomy / volume / face-
  matching / tree-DFS pipeline);
- Insufficient search (max_entry too low, taxonomy incomplete,
  k_max too low — the documented Q9c gaps);
- A genuine non-existence result that the literature has not
  recorded (i.e., a novel result Apeiron is on track to discover);
- Or: the search space we are exploring is the wrong one — the
  pivot suggested by Walton-Whittaker (markings) or Goodman-
  Strauss (don't expect monotile-ness to survive) is the right
  next move.

**Action items emerging from the read** (for `literature_review.md`):

- Verify Apeiron's Frettlöh ABCK rule satisfies sibling-edge-to-
  edge (not just sibling-vertex-to-vertex).
- Audit Apeiron's recognisability radius vs Goodman-Strauss's κ.
- Decide whether the Q10c Output 3 should be (a) extending the
  search at higher max_entry / k_max, or (b) widening from
  unmarked to marked monotile candidates (cf. Walton-Whittaker),
  or (c) a strategy_pivot.md acknowledging that the unmarked
  monotile route may be a dead end and that marked candidates
  are the natural extension.

---

## Reading log

| Paper                                  | Status            | Pages read |
|----------------------------------------|-------------------|------------|
| Goodman-Strauss 1998                   | Read (§1–§5)      | 38 / 45    |
| Frettlöh-Harriss 2013                  | Read (full)       | 9 / 9      |
| Frettlöh 2008 EJC                      | Read (full)       | 13 / 13    |
| Frettlöh 2008 companion                | Read (full)       | 7 / 7      |
| Walton-Whittaker 2019                  | Read (§1–§6.2)    | 23 / 27    |

Goodman-Strauss Appendix (pp.39–45), Walton-Whittaker §6.3+
(pp.24–27): not read; deferred unless the synthesis flags a need.
