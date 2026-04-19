# Candidates

Candidate tile definitions for the 3D Einstein search.

Each candidate lives in its own file (JSON for pure data, Python for
parameterised families) and declares:

- Vertex set in Z[phi]^3 (exact; no floats).
- Intended symmetry class (icosahedral, tetrahedral, trivial, …).
- Research track (Track A — deformation-first — or Track B —
  substitution-first; see CLAUDE.md §6).
- Aperiodicity claim level being targeted (weak / strong / strict;
  CLAUDE.md §2.1).
- Provenance: where the candidate came from (literature reference,
  deformation of a known set, algebraic search output, …).

Nothing here is load-bearing until it has a passing corona and substitution
test in `tests/`.
