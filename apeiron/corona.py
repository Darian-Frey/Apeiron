"""Corona BFS engine — the verification core.

A ``CoronaConfig`` is a tile surrounded by its immediate neighbours. The
first corona of P is the set of tiles sharing at least a face with P; the
k-th corona is defined inductively. The engine enumerates corona
configurations up to icosahedral symmetry and tests completion / forcing
properties.

Blocked on coplanar-face merging in ``polyhedron`` (CLAUDE.md §5.1): until
face adjacency is well-defined, corona completion tests return spurious
results.
"""
