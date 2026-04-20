"""Rhombic triacontahedron integration test.

First integration test under the ``≥ 2 core imports ⇒ integration``
rule (see [tests/integration/README.md](../integration/README.md)).
Spans ``apeiron.polyhedron``, ``apeiron.symmetry``, and
``apeiron.util``.

Acceptance criteria (CLAUDE.md §5.1 for ``polyhedron.py`` and
§3.3 for ``symmetry.py``):

- 32 vertices, 60 edges, 30 faces, Euler characteristic 2.
- Every face is a rhombus (4 vertices, 4 equal edges, parallelogram).
- The icosahedral rotation group I acts on the RTH as a symmetry:
  ``RTH.apply(g) == RTH`` for every g in ICOSAHEDRAL.
- I acts transitively on the 30 faces: the orbit of any one face
  under I covers the full face set.

The RTH fixture lives at ``candidates/rhombic_triacontahedron.json``
and is the first entry in the candidate library. Vertex-coordinate
decision: see ``docs/decisions.md``.
"""

from __future__ import annotations

from pathlib import Path

from apeiron.polyhedron import Polyhedron, _canonical_cycle
from apeiron.symmetry import ICOSAHEDRAL, Vec3
from apeiron.util import load_candidate
from apeiron.zphi import ZPhi

# Repo-relative path to the RTH fixture. tests are run from the repo
# root so relative-to-CWD works; if a future change moves that, switch
# to Path(__file__).resolve() based construction.
_RTH_PATH = Path("candidates/rhombic_triacontahedron.json")


def _rth() -> Polyhedron:
    """Load the RTH fixture. Factored out for reuse across tests."""
    return load_candidate(_RTH_PATH)


# -- combinatorics -----------------------------------------------------


class TestRhombicTriacontahedronCombinatorics:
    def test_thirtytwo_vertices(self) -> None:
        assert len(_rth().vertices) == 32

    def test_thirty_faces(self) -> None:
        assert len(_rth().faces) == 30

    def test_sixty_edges(self) -> None:
        P = _rth()
        edges: set[tuple[int, int]] = set()
        for face in P.faces:
            n = len(face)
            for idx in range(n):
                u, v = face[idx], face[(idx + 1) % n]
                edges.add((min(u, v), max(u, v)))
        assert len(edges) == 60

    def test_euler_characteristic_two(self) -> None:
        P = _rth()
        edges: set[tuple[int, int]] = set()
        for face in P.faces:
            n = len(face)
            for idx in range(n):
                u, v = face[idx], face[(idx + 1) % n]
                edges.add((min(u, v), max(u, v)))
        assert len(P.vertices) - len(edges) + len(P.faces) == 2


# -- rhombus shape -----------------------------------------------------


class TestRhombicFaces:
    def test_every_face_is_quadrilateral(self) -> None:
        for face in _rth().faces:
            assert len(face) == 4

    def test_every_face_has_equal_edge_lengths(self) -> None:
        P = _rth()
        for face in P.faces:
            v0, v1, v2, v3 = (P.vertices[i] for i in face)
            e0 = (v1 - v0).norm_squared()
            e1 = (v2 - v1).norm_squared()
            e2 = (v3 - v2).norm_squared()
            e3 = (v0 - v3).norm_squared()
            assert e0 == e1 == e2 == e3, (
                f"Face {face} has unequal edge lengths: {e0}, {e1}, {e2}, {e3}"
            )

    def test_every_face_is_parallelogram(self) -> None:
        # In a parallelogram the diagonals bisect each other, so
        # v0 + v2 == v1 + v3 (all component-wise in Z[phi]).
        P = _rth()
        for face in P.faces:
            v0, v1, v2, v3 = (P.vertices[i] for i in face)
            lhs = Vec3(v0.x + v2.x, v0.y + v2.y, v0.z + v2.z)
            rhs = Vec3(v1.x + v3.x, v1.y + v3.y, v1.z + v3.z)
            assert lhs == rhs, (
                f"Face {face} is not a parallelogram: "
                f"v0+v2={lhs} v1+v3={rhs}"
            )


# -- icosahedral symmetry ---------------------------------------------


class TestIcosahedralSymmetry:
    def test_every_rotation_preserves_polyhedron(self) -> None:
        P = _rth()
        for g in ICOSAHEDRAL:
            assert P.apply(g) == P

    def test_rotation_permutes_vertices_within_the_set(self) -> None:
        # Every g in I maps each RTH vertex to another RTH vertex.
        P = _rth()
        vertex_set = set(P.vertices)
        for g in ICOSAHEDRAL:
            for v in P.vertices:
                assert g.apply(v) in vertex_set


# -- face transitivity ------------------------------------------------


def _face_orbit_under_I(P: Polyhedron, face: tuple[int, ...]) -> set[tuple[int, ...]]:
    """Compute the orbit of ``face`` under the icosahedral group.

    Builds the permutation induced by each g in I on P's vertex
    indices (since I is a symmetry, g.apply(v_i) is another P-vertex,
    with a unique index). Applies each permutation to ``face`` and
    canonicalises the resulting cycle.
    """
    vertex_to_index: dict[Vec3, int] = {v: i for i, v in enumerate(P.vertices)}
    orbit: set[tuple[int, ...]] = set()
    for g in ICOSAHEDRAL:
        perm = [vertex_to_index[g.apply(v)] for v in P.vertices]
        rotated = tuple(perm[i] for i in face)
        orbit.add(_canonical_cycle(rotated))
    return orbit


class TestFaceTransitivity:
    def test_icosahedral_group_acts_transitively_on_faces(self) -> None:
        # I has order 60, RTH has 30 faces; each face's stabilizer has
        # order 2. The orbit of any face under I therefore has size 30
        # and covers the full face set.
        P = _rth()
        face_set = set(P.faces)
        orbit = _face_orbit_under_I(P, P.faces[0])
        assert orbit == face_set

    def test_orbit_size_is_thirty(self) -> None:
        P = _rth()
        orbit = _face_orbit_under_I(P, P.faces[0])
        assert len(orbit) == 30

    def test_every_face_has_stabilizer_of_order_two(self) -> None:
        # Orbit-stabilizer: |orbit| * |stabilizer| = |I| = 60.
        # With |orbit| = 30 (shown above), |stabilizer| = 2. Verify by
        # counting the g's that fix the face's canonical cycle.
        P = _rth()
        vertex_to_index: dict[Vec3, int] = {v: i for i, v in enumerate(P.vertices)}
        F0 = P.faces[0]
        stabilizer_count = 0
        for g in ICOSAHEDRAL:
            perm = [vertex_to_index[g.apply(v)] for v in P.vertices]
            rotated = tuple(perm[i] for i in F0)
            if _canonical_cycle(rotated) == F0:
                stabilizer_count += 1
        assert stabilizer_count == 2


# -- ring correctness --------------------------------------------------


class TestRingCorrectness:
    def test_no_floats_in_vertices_or_faces(self) -> None:
        # Every vertex component is a ZPhi with strict-int a, b;
        # every face index is a Python int. Ensures the loader's
        # scale-denom handling produced clean Z[phi] storage.
        P = _rth()
        for v in P.vertices:
            for comp in (v.x, v.y, v.z):
                assert isinstance(comp, ZPhi)
                assert isinstance(comp.a, int) and not isinstance(comp.a, bool)
                assert isinstance(comp.b, int) and not isinstance(comp.b, bool)
        for face in P.faces:
            for idx in face:
                assert isinstance(idx, int)

    def test_icosahedral_rotation_norm_preservation_across_vertices(self) -> None:
        # Every RTH vertex has the same squared-norm equivalence class
        # within its orbit under I (rotations preserve length exactly).
        P = _rth()
        for v in P.vertices:
            base = v.norm_squared()
            for g in ICOSAHEDRAL[:10]:   # subsample for speed
                assert g.apply(v).norm_squared() == base

