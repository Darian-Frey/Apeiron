# Audit: sibling-edge-to-edge condition for Track A's Danzer rule

**Trigger.** Q11a ruling (2026-05-04), part 2 of (c.3): verify
Apeiron's Frettlöh-derived Danzer ABCK rule satisfies Goodman-
Strauss's *sibling-edge-to-edge* condition (1998, §1.4) — the
strictly-stronger-in-d>2 version of sibling-vertex-to-vertex.

**Bottom-line finding (this audit).** The sibling-edge-to-edge
condition has **never been verified for Track A's Danzer rule**.
The Paolini-derived dissection encodes the geometry but no
existing test in Apeiron checks it. The condition is a finite
combinatorial property of
`candidates/danzer/paolini_dissection.json` — verifiable in
bounded compute, not requiring new mathematics.

**Status.** Per the Q10 sequencing commitment, the verification
*test* is not run in this audit. This document specifies the
test, identifies the inputs, and defers execution to post-relay.
Re-relay should report: "audit specifies the test; not yet run;
authorisation to run is the natural ask."

---

## §1. The condition (Goodman-Strauss 1998, §1.4)

### 1.1 Hereditary edges

A set of edges E for a substitution tiling (T, σ, S) is
**hereditary** if for every prototile A and every edge e ∈ E(A):

- σ(e) ⊂ σ(A) is exactly tiled by edges of children A^+;
- every k-facet (0 ≤ k ≤ d-2) of σ(e) is tiled by k-facets of
  children of A.

In d = 3, edges are 2-facets (= 2D faces of the tetrahedral
prototiles). So "edge" in §1.4 = "face" in 3D usage.

### 1.2 Sibling-edge-to-edge

For each pair of sibling children B, C ∈ A^+ (children of A in
σ(A)): if any k-facet (k < d-1) of any edge e of B is incident
to C in σ(A), then this facet is **exactly coincident** to some
k-facet of some edge f of C.

In d = 3 (Apeiron's setting):

- **k = d - 1 = 2**: edges (faces) themselves. Two siblings B, C
  meeting along a 2D region — that region must be a complete face
  of each.
- **k = 1**: edge-curves (1D). Two siblings meeting along a 1D
  segment — that segment must be a complete edge of each face's
  boundary.
- **k = 0**: vertices. Two siblings meeting at a 0D point — that
  point must be a complete vertex of each.

### 1.3 The 3D stricture (Lemma 1.3 + p.10 note)

**In d = 2:** hereditary edges ⇔ hereditary vertices. Sibling-
edge-to-edge ⇔ sibling-vertex-to-vertex.

**In d ≥ 3:** hereditary edges ⇒ hereditary vertices, but
**not the converse**. Goodman-Strauss p.10: "the converse is
probably not true when d > 2: there are likely to be tilings
that are sibling vertex-to-vertex that are not sibling edge-to-
edge."

**Practical implication.** A 3D substitution that is sibling-
vertex-to-vertex (siblings only meet at shared vertices in clean
ways) may *fail* sibling-edge-to-edge (siblings meet along edge
or face regions in messy ways — partial overlaps, not
coincidences).

---

## §2. Why this matters for Track A

### 2.1 The Paolini dissection encodes the candidate geometry

`candidates/danzer/paolini_dissection.json` lists 25 children
across 4 parent types:

| Parent | A | B | C | K | Total children                 |
|--------|---|---|---|---|--------------------------------|
| A      | 0 | 3 | 2 | 6 | 11                             |
| B      | 0 | 2 | 1 | 4 | 7                              |
| C      | 1 | 0 | 2 | 2 | 5                              |
| K      | 0 | 1 | 0 | 1 | 2                              |

Each child is recorded by `(parent, child_index, child_type,
is_proper, icosahedral_index, translation_x2)`. Together with the
ABCK prototile JSON files this fully specifies the geometry of
σ(A), σ(B), σ(C), σ(K).

### 2.2 What the existing tests verify

Per STATUS.md and `tests/integration/test_danzer_abck.py`:

- **Volume conservation**: σ(A) volume = sum of children's
  volumes. (Pillar-1 algebraic check.)
- **Multiplicity match**: child counts match the substitution
  matrix M_ABCK. (Pillar-1.)
- **I_h membership**: each child's icosahedral_index ∈ [0, 60)
  with is_proper flag. (Symmetry-group sanity.)
- **Pillar-2 recognisability**: `is_recognisable` returns success
  at radius=2 on σⁿ(P). (Subject to the §4 bridging-lemma caveat
  in [audit_kappa_vs_radius.md](audit_kappa_vs_radius.md).)
- **Pillar-3 inflation argument**: the four-function pipeline
  composes; PF eigenvalue extracted as φ³.

**None of these check sibling-face-to-face.** Volume and
multiplicity are insensitive to whether siblings overlap at
their shared faces in clean coincidences or messy partial
intersections.

### 2.3 Why not catching this earlier matters

If sibling-face-to-face fails for Track A's Danzer rule, then:

- The rule falls outside Goodman-Strauss's enforcement scope.
  The matching-rule construction in §3 cannot be applied
  directly.
- Apeiron's Track A pillar-4 path (`fourth_pillar.py`) cites
  Goodman-Strauss 1998 as one source of the no-non-hierarchical
  argument. If the rule fails the precondition, this citation is
  void.
- The ABCK substitution may still be geometrically valid — it
  tiles ℝ³ and admits the substitution — but the *enforcement*
  by matching rules requires extra work (e.g., refining edge
  sets to make them hereditary).

If sibling-face-to-face holds, the rule is **inside** Goodman-
Strauss's framework and the existing pillar-4 citation stands.

### 2.4 What the literature says about ABCK and face-to-face

**Frettlöh (Tilings Encyclopedia / Theorem 1.2).** Every face-to-
face tiling by A, B, C, K with the blue-edge mirror matching rule
belongs to the substitution hierarchy. This is a *tiling-level*
face-to-face statement: the tiling, viewed as a covering of ℝ³
by isometric copies of A, B, C, K, is face-to-face.

**Tiling-level face-to-face does NOT imply sibling-face-to-face
within σ(A).** The tiling's face-to-face property is about how
the four prototile types meet across the entire tiling. Sibling-
face-to-face is about how children within σ(A) meet each other,
*before* the substitution propagates them into the tiling.

The two are related but not identical: tiling-level face-to-face
is necessary for sibling-face-to-face (otherwise σ(A) wouldn't be
face-to-face inside σ²(A)), but not sufficient. A substitution
could be tiling-face-to-face yet have sibling overlaps at the
substitution-level that get smoothed only after one or more
inflation steps.

### 2.5 The expected outcome

ABCK is a well-studied 3D substitution. The Paolini POV-Ray
rendering shows clean face-to-face meetings within σ(A), σ(B),
σ(C), σ(K). The expected outcome of the audit is:
**sibling-face-to-face holds.** But this is an expectation, not a
proof. The audit proposes to verify computationally.

---

## §3. The verification procedure

### 3.1 Inputs

- `candidates/danzer/paolini_dissection.json` (25 children).
- `candidates/danzer/{A,B,C,K}.json` (prototile geometry).
- `apeiron.zphi.ZPhi`, `apeiron.symmetry.ICOSAHEDRAL`,
  `apeiron.symmetry.ImproperRotation` for placing children.

### 3.2 The check

For each parent P ∈ {A, B, C, K}:

1. Build `placed_children: list[Polyhedron]` by applying each
   child's `(rotation, translation)` to its prototile geometry.
2. For each pair (i, j) of sibling children with i < j:
   1. Compute the intersection `placed_children[i] ∩
      placed_children[j]` exactly in ZPhi³.
   2. Classify the intersection by dimension: empty, point (0D),
      segment (1D), polygon (2D), polytope (3D).
   3. **3D intersection** ⇒ siblings overlap in interior.
      Substitution invalid (volume conservation would have
      caught this; sanity check).
   4. **2D intersection** ⇒ siblings share a face region. Check:
      is the intersection an *entire face* of both child[i] and
      child[j]? If yes, sibling-face-to-face holds for this
      pair. If the intersection is a strict sub-region of either
      child's face, sibling-face-to-face **fails**.
   5. **1D intersection** ⇒ siblings share an edge segment.
      Check: is the intersection an entire edge of both
      child[i] and child[j]? Similar test as 2D.
   6. **0D intersection** ⇒ siblings share a vertex. Trivially
      sibling-vertex-to-vertex.
   7. **Empty intersection** ⇒ no shared boundary. No constraint.

3. For sibling-face-to-face: every 2D intersection in (2.iv)
   above must be a full face of both children.

### 3.3 Implementation cost

- Exact-arithmetic polyhedron intersection in ZPhi: not
  implemented in Apeiron currently. Closest existing primitive
  is `_convex_polyhedra_interior_disjoint` (SAT-based, returns
  bool only).
- Face equality: implementable via canonical form on the merged
  face polygon (existing `Polyhedron.from_raw` machinery should
  handle this).
- Total: ~2 days of engineering for a self-contained
  `sibling_edge_to_edge_audit(rule, prototile_geom)` predicate.

### 3.4 Defer per Q10 sequencing

Per Q11a → Q10's pre-code-relay commitment: the implementation
should not land before re-relay. This audit's deliverable is the
specification above; execution awaits authorisation.

---

## §4. What to put in the relay-back

**Two findings from (c.3):**

1. **κ vs Apeiron's radius.** They measure different things.
   Apeiron's predicate is patch-level sound. Hull-level
   recognisability needs a bridging lemma (n ≥ κ patch-size
   condition) that should be added to `is_recognisable`'s
   docstring contract. Track A's pillar-2 result at radius=2 is
   likely sound (assuming κ_ABCK ≤ 2 and patch n ≥ κ_ABCK), but
   needs verification. Concrete action items at
   [audit_kappa_vs_radius.md §5](audit_kappa_vs_radius.md).

2. **Sibling-edge-to-edge.** Never verified. Existing tests
   check tiling-level properties (volume, multiplicity, I_h
   membership) but not sibling-level face-to-face within σ(A).
   The Paolini dissection looks face-to-face by visual
   inspection, but rigorous verification requires a new
   `sibling_edge_to_edge_audit` predicate (~2 days
   engineering).

**Re-relay sub-questions:**

- **Q12a.** Authorise running the sibling-face-to-face
  verification on Track A's Danzer rule? If yes, this is the
  first concrete (c.3) follow-up code task.
- **Q12b.** Authorise the κ-audit follow-ups? Specifically
  (audit_kappa_vs_radius.md §5.1) updating
  `is_recognisable`'s docstring with the bridging-lemma
  contract, and (§5.2) verifying Track A's pillar-2 patch was
  σⁿ for n ≥ κ_ABCK.
- **Q12c.** Given the audits don't change the (c.1) → (c.2)
  next steps, should we proceed with (c.1) gap-closure
  (taxonomy, k_max, max_entry=3) in parallel with the (c.3)
  follow-ups? Or serialise: finish (c.3) first, then (c.1),
  then the cut-and-project literature read?
