# Decisions

This file is a chronological log of *non-obvious* design choices — the
"we tried X, it failed because Y, so we do Z" record that doesn't
belong in [CLAUDE.md](../CLAUDE.md) (durable guidance) or
[STATUS.md](../STATUS.md) (operational state). If you're reviewing code
and wondering "why is this particular coordinate / representation /
convention chosen?", this is where the answer should live.

Entries accumulate. Never delete an entry — if a decision is
superseded, add a new entry that supersedes the old one and cross-link.

## Format

Each entry is a level-2 heading with ISO date, a one-line topic, and
the SHA of the commit that carries the decision into the code. Body:
**Context**, **Decision**, **Rationale**, and, where relevant,
**Alternatives considered**.

---

## 2026-04-20 — RTH vertex coordinate choice (fba7ae7)

**Context.** The rhombic triacontahedron is the canonical icosahedral
Catalan solid: 32 vertices (20 three-valent + 12 five-valent),
60 edges, 30 rhombic faces. It's the first integration fixture for the
`polyhedron.py` stack (commit E of the §5.1 breakdown), and its
coordinates need to live natively in ℤ[φ]³.

**Decision.** The `candidates/rhombic_triacontahedron.json` vertex set
is:

- 8 cube-type dodec vertices: `(±1, ±1, ±1)`
- 12 non-cube dodec vertices: cyclic permutations of `(0, ±φ, ±1/φ)`
- 12 icosahedral vertices: cyclic permutations of `(0, ±1, ±φ)`

All 32 coordinates are integer ZPhi (since `1/φ = φ − 1 ∈ ℤ[φ]`), so
`scale_denom = 1` in the loader schema.

**Rationale.** The first attempted spec — cube `(±φ, ±φ, ±φ)`, dodec
`(0, ±1, ±φ²)`, icos `(0, ±φ, ±φ²)` — produces **only a 20-vertex
convex hull**, not 32. Cause: the dodec and icos vertex sets share the
cyclic-perm axes `(x = 0, z = ±φ²)`, `(y = 0, x = ±φ²)`, and
`(z = 0, y = ±φ²)`. On each such axis, four points are collinear — two
dodec (at `y ∈ {±1}`) strictly between two icos (at `y ∈ {±φ}`). The
dodec points are convex combinations:

    (0, 1, φ²) = ((φ+1)/(2φ)) · (0, +φ, φ²) + ((φ-1)/(2φ)) · (0, -φ, φ²)

Both coefficients are in `(0, 1)`, so `(0, 1, φ²)` is strictly interior
to the edge between the two icos vertices. All 12 dodec non-cube
vertices fall into this pattern across the three cyclic perms × four
sign combinations, collapsing the "RTH" to a 20-vertex polytope (cube
+ icos only).

The corrected spec above uses `(0, ±φ, ±1/φ)` for the non-cube dodec
positions (the "large" slot holds `φ`, not `φ²`; the "small" slot
holds `1/φ = φ − 1`, not `1`). The dodec and icos z-levels in each
cyclic perm now differ: dodec `(0, ±φ, ±1/φ)` has z ∈ `{±1/φ}`, icos
`(0, ±1, ±φ)` has z ∈ `{±φ}`. No shared lines, no collinear
quadruples, every one of the 32 vertices is an extreme point of the
hull.

The scipy combinatorial oracle plus exact ZPhi validation caught this
on the first failing spec. The validation in `Polyhedron.from_vertices`
correctly rejected the dodec vertex as "classified interior but on a
face boundary", and hand-verification confirmed the collinearity.

**Alternatives considered.**

- *Scale the wrong spec differently.* Multiplying all coords by a
  constant doesn't fix collinearity — it's a projective property
  preserved under scaling.
- *Use non-standard cyclic-perm conventions.* Swapping dodec and icos
  naming doesn't help (symmetric in the two sets); permuting to odd
  permutations gives the same 12-point set as the even/cyclic ones.
- *From-scratch derivation as face centroids of the
  icosidodecahedron.* Correct but expensive; the dual-construction
  face centroids turn out to coincide with the icos + dodec vertex
  positions only at specific scales, which is exactly what the
  corrected spec encodes directly.

Cross-reference: this event drove commit
`dd36f96 feat(util): reject interior input vertices by default` —
`load_candidate` now raises by default when scipy drops any input
vertex, so a future spec error of this class surfaces at load time
with the specific interior indices, not later as a corona-BFS
anomaly.
