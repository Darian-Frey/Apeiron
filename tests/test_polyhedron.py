"""Tests for apeiron.polyhedron.

Critical fixtures (CLAUDE.md §7.4):

- Rhombic triacontahedron — 30 rhombic faces, icosahedral symmetry; the
  canonical target for the coplanar-face merging pass.
- Rhombic dodecahedron — 12 rhombic faces; secondary fixture for the merge
  pass.

Both have independently-known face structures, so the post-merge face count
and per-face vertex cycles are decidable ground truth.
"""
