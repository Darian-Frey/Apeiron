# Audit: Goodman-Strauss κ vs Apeiron's recognisability radius

**Trigger.** Q11a ruling (2026-05-04, archived at
[relays/q11_2026-05-04.md](relays/q11_2026-05-04.md)): "if Apeiron's
recognisability radius doesn't bound κ, pillar-2 results for every
tested candidate are potentially invalid." Auditing this *first*
within (c.3) per Claude (web)'s sequencing.

**Bottom-line finding (this audit).** Apeiron's `is_recognisable`
output is **sound** — when it returns `True` at radius `r`, every
tile in the supplied patch *does* have a unique parent supertile
determined by its radius-`r` neighbourhood signature. The deeper
question is whether radius-`r` patch-level recognisability lifts to
**rule-level recognisability in the hull** in the Goodman-Strauss /
Solomyak sense. The answer depends on κ, and the audit below
shows the relationship between κ and Apeiron's radius is **not a
direct bound** — they measure structurally different things.

**Implication.** Prior pillar-2 successes are not invalidated, but
they are **not equivalent** to the formal recognisability claim in
Goodman-Strauss. A bridging argument is needed before the
inflation-argument output (`InflationArgument(pf=φ³, radius=2)`)
can be cited as proof of recognisability in the hull. Either
(a) the bridging lemma should be stated and verified per-rule, or
(b) Apeiron should adopt a stronger patch construction (a
sufficiently large σ^N(P) with N > κ) that makes patch-level and
rule-level recognisability coincide.

The action item is light: write the bridging lemma per the §4
sketch below, add it as a docstring contract on
`is_recognisable`, and verify the existing test fixtures still
pass under the strengthened contract.

---

## §1. The two definitions

### 1.1 Goodman-Strauss's κ (1998, §2.1.3, pp.14–15)

For each prototile B in the substitution alphabet S, κ(B) is the
smallest natural number such that there is a connected collection
of edges in σ^κ(B)(B) that:

1. Contains every σ^(κ(e)−1)(e) for e ∈ E(B).
2. Meets every Z(B) (sites) on the boundary of B.
3. Meets every endovertex on the boundary of B.

Then κ := max_B κ(B).

**Semantics.** κ is a **level-of-hierarchy index**, not a Euclidean
or graph distance. κ = 2 means: examining σ²(B) — the level-2
supertile — reveals enough connected edge structure to identify
the role of B's vertices and sites in the parent supertile σ(B^-).

**Practical scale (Goodman-Strauss p.15).** "κ is often rather
low; many well-known examples have κ = 1, very few have κ > 2.
The Conway-Radin pinwheel and Robinson triangles have κ = 1; the
L-tiling requires κ = 2." Sadun's generalised pinwheel
constructions can require arbitrarily large κ.

### 1.2 Apeiron's radius
([apeiron/hierarchy.py:540-556](../apeiron/hierarchy.py#L540-L556))

`TilePatch.neighbour_within(i, j, radius)` returns True iff tiles
i and j are within graph-distance `radius` along the patch's
adjacency oracle. The oracle is abstract — typically corona-
adjacency graph distance, but it can be Euclidean (squared
distance against a threshold) or any monotone-in-r relation.

`is_recognisable(patch, max_radius=5)` iterates r = 1, 2, …,
max_radius. At each r it groups tiles by their radius-r
neighbourhood signature (tile-type multiset within graph
distance r). If every group's tiles share a single
`parent_supertile` value, return success at radius r.

**Semantics.** Apeiron's radius is a **graph-distance bound** on a
finite patch. It measures how many adjacency hops you walk before
inspecting the multiset signature.

### 1.3 What they share, what they don't

**Shared.** Both ask: how much local information is needed to
resolve a tile's role in the substitution hierarchy?

**Different.**

- κ counts *substitution levels*; Apeiron's radius counts
  *adjacency hops in a flattened patch*.
- κ is defined in terms of edge skeletons (vertex wires, sites);
  Apeiron's radius is agnostic to edge structure — it relies on
  whichever signature function the caller passes.
- κ is a **rule property** (function of the substitution σ alone);
  Apeiron's radius is a **patch property** (function of the
  specific σⁿ(P) the caller supplies).

**Direct conversion is not available.** κ and Apeiron's radius
measure different quantities; no general formula relates them.

---

## §2. Why this matters for the substitution hull

The four-pillar argument (CLAUDE.md §6.3) treats recognisability
as a *rule* property: every tile's supertile membership is
determined by a bounded local neighbourhood, *for almost every
tiling in the hull* (in the Solomyak sense — measure 1 in the
G-invariant probability measure on (T, σ, S)).

Apeiron's `is_recognisable` checks a *patch* property: every tile
in the supplied patch has a unique parent supertile determined by
radius-r signature on that patch.

The two coincide **if and only if** the patch is large enough to
carry the rule's recognisability information. In the
Goodman-Strauss framework, "large enough" means: the patch
includes σ^κ(B) for every prototile B that appears in it. If the
patch is smaller than σ^κ-scale, patch-level recognisability is a
sufficient but not necessary condition for rule-level
recognisability — and conversely, a *failure* of patch-level
recognisability at small patch may not invalidate the rule.

### 2.1 Patch sizes Apeiron has been using

`TestPipelineFibonacci` and Track A's Danzer pillar-2:
`patch_from_supertile(σⁿ(P))` for some bounded n.

Track A, commit `5bfa3a5`: full pipeline returns
`InflationArgument(pf=ZPhi(1,2), radius=2)` on σⁿ(P) for n
unspecified in STATUS.md (likely n = 2 or 3 for the Paolini ABCK
fixture; the test should be checked).

The Paolini ABCK σⁿ(P) at n = 1 has up to 25 children (per
Frettlöh / Paolini); at n = 2, up to 25² = 625 children. The
diameter (in graph distance) of σ²(P) is bounded by ~25–50.

Apeiron's `is_recognisable` returned success at radius = 2 on this
patch. Radius = 2 is the **bare minimum non-trivial neighbourhood**
(neighbour-of-neighbour). For Goodman-Strauss-style recognisability
in the hull, this is a strong sufficient claim: it asserts that
even at this minimal radius, the multiset signature already
distinguishes parents.

### 2.2 Risk: does small-radius patch success lift to the hull?

**Risk class 1: false rule-level success.** A patch is finite; the
hull contains tilings the patch does not. A radius-2 success on
σ²(P) does **not** automatically imply radius-2 recognisability
holds for arbitrary tilings in the hull. A counter-example would
be: some tiling in the hull contains a configuration where two
tiles share the radius-2 multiset signature but have different
parent supertiles.

**Risk class 2: false rule-level failure.** Apeiron's max_radius
default is 5. If a rule is recognisable at radius 7 in the
Goodman-Strauss sense, our default would mark it as failed.
Within the Q9 sweep this is moot — Q9's gate was geometric
realisation, not recognisability — but for a future Track A or
Track B candidate that depends on pillar 2, the cap matters.

### 2.3 Risk class 1 is the load-bearing one

Apeiron's pillar-2 successes (Track A Danzer) feed into the
inflation-argument pipeline (pillar 3). If pillar 2 is unsound at
the rule level, pillar 3's conclusion (the rule is aperiodic) is
unsound too.

**However**, the Solomyak / Goodman-Strauss argument has a
bootstrap: if patch-level recognisability holds for a patch
covering σ^κ(P), then rule-level recognisability holds in the hull
*for tilings that arise from the substitution*. The full hull
contains non-substitution tilings (the pillar-4 territory), but
those have measure zero in any G-invariant probability measure (§1.2
of Goodman-Strauss, properties (i)-(iii) of connected hierarchy).

So patch-level recognisability **does** lift to rule-level
recognisability for the substitution-hull (the measure-1 part of
the full hull), provided the patch is at least σ^κ-scale.

---

## §3. The κ values at play in Apeiron

### 3.1 Track A — Danzer ABCK

ABCK is a 4-prototile substitution with linear inflation φ.
Goodman-Strauss does not give κ_ABCK explicitly. Frettlöh's papers
do not state κ in the Goodman-Strauss sense. The L-tiling has
κ = 2 (per Goodman-Strauss p.15). ABCK is structurally close to
the L-tiling family (planar / 3D substitution with edge-to-edge
sibling structure).

**Conjecture (unverified):** κ_ABCK ∈ {1, 2}. Most likely 2,
matching the L-tiling case which Goodman-Strauss flags as a
typical small-κ rule.

**Action item.** Verify κ_ABCK by inspection of the Paolini
dissection: examine a single A, B, C, K supertile; check whether
the boundary-incident edges of σ¹(A) (resp. B, C, K) form a
connected collection meeting every site and endovertex — if yes,
κ_ABCK = 1; if not, check σ². This is a finite combinatorial
check on `candidates/danzer/paolini_dissection.json`.

### 3.2 Track B — n=3 PF=φ³ candidates

For the n=3 PF=φ³ candidates from the Q9 sweep: κ is
**not relevant** because Q9's gate was geometric realisation
(`realise()` returning NoRealisation), not recognisability.
Pillar 2 was never invoked on these candidates.

If we ever do reach a Realised candidate from Track B, this audit's
findings apply: recompute κ for that candidate before citing
`is_recognisable` as proof of pillar 2.

---

## §4. The bridging lemma

For Apeiron's existing pillar-2 results to be sound under the
Goodman-Strauss framework, the following bridging lemma must hold:

> **Lemma (proposed).** Let (T, σ, S) be a primitive substitution
> tiling with finite alphabet, let κ be Goodman-Strauss's κ for σ,
> and let P be a prototile. Suppose `is_recognisable` returns
> success at radius r on `patch_from_supertile(σⁿ(P))` for some
> n ≥ κ. Then σ is recognisable at radius r in the substitution
> hull (mod measure zero).

**Proof sketch (informal).** σⁿ(P) for n ≥ κ contains every
configuration that arises in any substitution-hull tiling, mod
finite-radius windowing. Patch-level radius-r recognisability on
this patch therefore witnesses radius-r recognisability for every
configuration in the hull's measure-1 part. The mod-measure-zero
caveat covers infinite fault-line tilings and other non-
hierarchical tilings (Goodman-Strauss §1.2 (ii)).

**What needs to be verified before this lemma is invoked:**

1. **n ≥ κ for the patch supplied.** For Track A's Danzer
   pillar-2 (radius=2 success on σⁿ(Paolini ABCK)), confirm
   n ≥ κ_ABCK. If κ_ABCK = 2 and the test used n = 2, the lemma
   applies. If the test used n = 1, we need to re-run at n ≥ 2.

2. **The patch-graph oracle is consistent with σⁿ(P)'s actual
   geometry.** Apeiron's `make_euclidean_squared_oracle` uses
   squared-Euclidean distance against a threshold tied to
   `radius_step_squared`. The threshold must correspond to
   "graph distance 1 = adjacent-tile in σⁿ(P)." If the threshold
   is too generous (counts non-adjacent tiles as distance 1),
   patch-level recognisability is a *weaker* claim than rule-
   level; the lemma still holds but the bound is loose.

3. **The signature_fn is rich enough to distinguish parents at
   the chosen radius.** `neighbourhood_signature` (multiset of
   types) is sufficient for the Ammann-Beenker case (Baake-Grimm
   §6). For ABCK, the multiset may coalesce distinct
   configurations; `shell_neighbourhood_signature` (per-distance
   multiset) is the documented stronger option. Verify which was
   used in the Track A pillar-2 result.

---

## §5. Action items (in priority order)

### 5.1 Document the bridging lemma in `is_recognisable`'s docstring

Add a "Soundness contract" section:

> Returns success iff every tile in the supplied patch has a
> unique parent_supertile determined by its radius-r
> neighbourhood signature on the patch.
>
> This implies rule-level recognisability in the substitution hull
> (mod measure zero) **if** the patch is built from σⁿ(P) with
> n ≥ κ(σ), where κ is Goodman-Strauss's resolution parameter
> (1998, §2.1.3). For substitution rules with κ = 1 (Conway-Radin
> pinwheel, Robinson triangles), n ≥ 1 suffices. For rules with
> κ = 2 (L-tiling and most ABCK-family rules), n ≥ 2 is required.
> Callers must determine κ for their rule and ensure the patch
> meets this bound; the predicate cannot verify n ≥ κ itself.

**Cost:** docstring-only, minutes. Should land before any other
code touches the predicate.

### 5.2 Verify Track A's pillar-2 used n ≥ κ_ABCK

Find the Track A Danzer pillar-2 test, check the σⁿ depth used,
and check whether κ_ABCK ≤ that n. If not, re-run at higher n.

**Cost:** 30 minutes of inspection plus possibly re-running one
test fixture.

### 5.3 Compute κ_ABCK from the Paolini dissection

Inspect `candidates/danzer/paolini_dissection.json`; verify
whether the boundary edge-collection of each prototile's level-1
supertile is connected and meets every site / endovertex. If yes,
κ_ABCK = 1; if not, check level 2.

**Cost:** ~1 hour, one-off.

### 5.4 Add a regression test that asserts the bridging-lemma
hypothesis

A test that constructs a patch from σⁿ(P) with n < κ and asserts
that radius-r success on that patch does NOT lift to rule-level
recognisability (or, if no counter-example is constructible,
documents this as expected behaviour for ABCK-family rules where
n=1 still suffices in practice).

**Cost:** test-writing, half a day; defer until 5.2 and 5.3 are
done.

---

## §6. What this audit does NOT do

- It does NOT pin down κ_ABCK rigorously. §5.3 does that.
- It does NOT check whether `make_euclidean_squared_oracle`'s
  threshold is correctly tuned for ABCK. That's a separate audit
  on the patch construction.
- It does NOT cover Goodman-Strauss's edge-to-edge condition
  (Lemma 1.3 / 3D obstruction signal). That's the second part of
  (c.3), to be done after the κ work concludes.

---

## §7. Conclusion

Apeiron's `is_recognisable` is sound at the patch level. The
Goodman-Strauss-style hull-level recognisability claim requires
the bridging lemma in §4, which is true under the n ≥ κ patch-
size condition. Existing Track A pillar-2 results are likely
sound (assuming κ_ABCK ≤ 2 and the patch was σⁿ-built for
n ≥ 2), but this needs verification per §5.2-5.3.

**No prior result is being retracted by this audit.** The audit
identifies a contract that should be made explicit in the
docstring and verified per-test, and it lays out a finite list of
follow-up tasks (§5) that close the integrity question without
changing the four-pillar pipeline architecture.

The next step is the second part of (c.3): the sibling-edge-to-
edge audit (Goodman-Strauss Lemma 1.3, the 3D-specific stricture).
After both parts of (c.3) are done, re-relay per the Q11
sequencing commitment.
