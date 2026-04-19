"""Polyhedra with vertices in Z[phi]^3.

Defines ``Vertex``, ``Edge``, ``Face``, and ``Polyhedron`` dataclasses with
exact Z[phi] coordinates. Supports isometry application (via the matrices in
``symmetry``) and face/edge adjacency queries.

KNOWN BLOCKER (CLAUDE.md §5.1): construction from a vertex set via convex
hull must include a coplanar-face merging pass. Adjacent triangles whose
supporting planes are exactly equal in Z[phi] must be merged into a single
polygonal face before any downstream corona or adjacency logic runs.
Equality of supporting planes must be decided exactly in Z[phi] — float
tolerances are not acceptable.
"""
