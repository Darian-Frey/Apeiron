"""Sibling-face-to-face audit on Track A's Paolini-derived Danzer rule.

Per Q12a authorisation (Claude (web), 2026-05-04). Verifies the
Goodman-Strauss §1.4 condition: for every pair of sibling tetrahedra
(B, C) within σ(A) (or σ(B), σ(C), σ(K)) that meet along a 2-facet
(face), the meeting region is *exactly coincident* to a complete face
of both B and C. The check has two passes:

(1) Vertex-coincidence pass. For each sibling pair, find shared
    vertex coordinates. Then for every face-pair (f_a from sibling i,
    f_b from sibling j), test whether they share a supporting plane.
    If yes, test whether the face-vertex sets coincide. A
    coplanar-but-different face-pair is the sibling-face-to-face
    *failure* signal.

(2) Interior-disjoint pass. Reuses Apeiron's existing exact-ZPhi
    SAT primitive ``_convex_polyhedra_interior_disjoint`` to verify
    no sibling pair has overlapping interior. (Volume conservation
    already implies this; SAT is the geometric verification.)

Per Claude (web)'s instruction: on failure, do **not** modify the
Paolini dissection — flag the offenders and document.

Run from repository root: ``python -m scripts.audit_sibling_face_to_face``.
"""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

from apeiron.polyhedron import Polyhedron
from apeiron.symmetry import ICOSAHEDRAL, ImproperRotation, Rotation, Vec3
from apeiron.track_b.realisation import _convex_polyhedra_interior_disjoint
from apeiron.util import load_candidate
from apeiron.zphi import ZPhi

_DANZER_DIR = Path("candidates/danzer")


def _place(
    proto: Polyhedron,
    rotation: Rotation | ImproperRotation,
    translation_x2: Vec3,
) -> Polyhedron:
    """Apply rotation + translation (in ×2 storage) to a prototile.

    Returns a fresh ``Polyhedron`` with same face-index structure but
    transformed vertex coordinates. Uses
    ``Polyhedron.apply(rotation)`` for the rotation (which preserves
    face structure) then re-wraps with the translated vertices.
    """
    rotated = proto.apply(rotation)
    new_vertices = tuple(v + translation_x2 for v in rotated.vertices)
    return Polyhedron(
        vertices=new_vertices,
        faces=rotated.faces,
    )


def _cross(a: Vec3, b: Vec3) -> Vec3:
    return Vec3(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x,
    )


def _is_zero(v: Vec3) -> bool:
    Z = ZPhi(0, 0)
    return v.x == Z and v.y == Z and v.z == Z


def _planes_coincide(
    face_a_verts: tuple[Vec3, Vec3, Vec3],
    face_b_verts: tuple[Vec3, Vec3, Vec3],
) -> bool:
    """Exact ZPhi test: do these two triangles support the same plane?

    Two triangles support the same plane iff their normals are
    parallel (cross product vanishes) AND a vertex of one lies on
    the plane of the other.
    """
    a0, a1, a2 = face_a_verts
    b0, b1, b2 = face_b_verts
    n_a = _cross(a1 - a0, a2 - a0)
    n_b = _cross(b1 - b0, b2 - b0)
    if not _is_zero(_cross(n_a, n_b)):
        return False
    # Normals are parallel; check b0 lies on plane of (a0, a1, a2).
    return (b0 - a0).dot(n_a) == ZPhi(0, 0)


def _vertex_set(verts: tuple[Vec3, Vec3, Vec3]) -> frozenset[Vec3]:
    return frozenset(verts)


def _inflate_polyhedron_phi(p: Polyhedron) -> Polyhedron:
    """Return τP — multiply each vertex by φ (linear inflation).

    Vertices stored in ×2 form; multiplying by φ = ZPhi(0, 1)
    preserves the storage convention because 2·(φ·v_real) = φ·(2·v_real).
    """
    phi = ZPhi(0, 1)
    new_vertices = tuple(
        Vec3(v.x * phi, v.y * phi, v.z * phi) for v in p.vertices
    )
    return Polyhedron(vertices=new_vertices, faces=p.faces)


def _coplanar_triangles_overlap_in_2d(
    face_a_verts: tuple[Vec3, Vec3, Vec3],
    face_b_verts: tuple[Vec3, Vec3, Vec3],
) -> bool:
    """SAT in 3D restricted to candidate axes lying in the shared plane.

    Caller must have already verified the two triangles are coplanar.
    Returns True iff the triangles overlap with positive 2D area in
    their shared plane (touching at a single edge or vertex returns
    False — that's a 1D / 0D meeting, not 2D overlap).

    Candidate separating axes: for each of the 3 edges of each
    triangle (6 total), take ``edge × plane_normal``. If any
    candidate axis separates the two triangles' projections (with
    strict inequality), they don't overlap. Else they do.
    """
    a0, a1, a2 = face_a_verts
    n = _cross(a1 - a0, a2 - a0)  # plane normal (any scaling)
    edges = [
        (a1 - a0), (a2 - a1), (a0 - a2),
        (face_b_verts[1] - face_b_verts[0]),
        (face_b_verts[2] - face_b_verts[1]),
        (face_b_verts[0] - face_b_verts[2]),
    ]
    for e in edges:
        axis = _cross(e, n)
        if _is_zero(axis):
            continue
        proj_a = [v.dot(axis) for v in face_a_verts]
        proj_b = [v.dot(axis) for v in face_b_verts]
        max_a, min_a = max(proj_a), min(proj_a)
        max_b, min_b = max(proj_b), min(proj_b)
        # Strict separation: closure-disjoint intervals.
        if max_a <= min_b or max_b <= min_a:
            return False
    return True


def _is_supertile_boundary_face(
    face_verts: tuple[Vec3, Vec3, Vec3],
    tau_P: Polyhedron,
) -> bool:
    """Does this face's supporting plane coincide with one of τP's
    boundary face-planes?

    True ⇒ the face lies on σ(P)'s outer boundary (expected for many
    children's faces; sibling co-tiling along τP boundary is **not**
    a sibling-face-to-face violation — it's just supertile-boundary
    decomposition).
    """
    for fp in tau_P.faces:
        if len(fp) != 3:
            continue
        boundary_verts = tuple(tau_P.vertices[k] for k in fp)
        if _planes_coincide(face_verts, boundary_verts):
            return True
    return False


def audit_parent(
    parent_letter: str,
    children_data: list,
    prototypes: dict[str, Polyhedron],
) -> dict:
    """Run the audit for one parent's σ(P).

    Returns a dict with keys:
      - 'placed': list of placed child ``Polyhedron``s
      - 'pairs_n': number of (i, j) pairs checked
      - 'interior_disjoint_failures': list of (i, j) where SAT
        returned not-disjoint
      - 'face_matches': count of face-pairs where same plane AND
        same vertex set (the expected face-to-face coincidences,
        not on σ(P)'s boundary)
      - 'boundary_co_tilings': count of face-pairs where the shared
        plane is a τP boundary face-plane (expected; not a violation)
      - 'internal_face_failures': list of dicts describing
        coplanar-but-different face pairs in the σ(P) interior —
        the actual sibling-face-to-face violations to investigate
    """
    placed: list[Polyhedron] = []
    child_types: list[str] = []
    for ch in sorted(children_data, key=lambda r: r["child_index"]):
        proto = prototypes[ch["child_type"]]
        rot_proper = ICOSAHEDRAL[ch["icosahedral_index"]]
        rot = rot_proper if ch["is_proper"] else ImproperRotation(rot_proper)
        t_x2 = ch["translation_x2"]
        translation = Vec3(
            ZPhi(*t_x2[0]), ZPhi(*t_x2[1]), ZPhi(*t_x2[2])
        )
        placed.append(_place(proto, rot, translation))
        child_types.append(ch["child_type"])

    tau_P = _inflate_polyhedron_phi(prototypes[parent_letter])

    interior_disjoint_failures: list[tuple[int, int]] = []
    internal_face_failures: list[dict] = []
    face_matches = 0
    boundary_co_tilings = 0

    for i in range(len(placed)):
        for j in range(i + 1, len(placed)):
            P_i = placed[i]
            P_j = placed[j]

            if not _convex_polyhedra_interior_disjoint(
                P_i.vertices, P_i.faces,
                P_j.vertices, P_j.faces,
            ):
                interior_disjoint_failures.append((i, j))

            for fi_idx, fi in enumerate(P_i.faces):
                if len(fi) != 3:
                    continue
                fi_verts = tuple(P_i.vertices[k] for k in fi)
                for fj_idx, fj in enumerate(P_j.faces):
                    if len(fj) != 3:
                        continue
                    fj_verts = tuple(P_j.vertices[k] for k in fj)
                    if not _planes_coincide(fi_verts, fj_verts):
                        continue
                    on_boundary = _is_supertile_boundary_face(fi_verts, tau_P)
                    if _vertex_set(fi_verts) == _vertex_set(fj_verts):
                        if on_boundary:
                            boundary_co_tilings += 1
                        else:
                            face_matches += 1
                    else:
                        if on_boundary:
                            boundary_co_tilings += 1
                            continue
                        # Internal coplanar face-pair with different
                        # vertex sets. Only a sibling-face-to-face
                        # violation if the triangles actually overlap
                        # in 2D — otherwise they're just two sibling
                        # faces in the same plane that don't meet.
                        if _coplanar_triangles_overlap_in_2d(
                            fi_verts, fj_verts,
                        ):
                            internal_face_failures.append({
                                "child_i": i,
                                "child_i_type": child_types[i],
                                "child_j": j,
                                "child_j_type": child_types[j],
                                "face_i": fi,
                                "face_j": fj,
                                "verts_i": fi_verts,
                                "verts_j": fj_verts,
                            })

    return {
        "placed": placed,
        "child_types": child_types,
        "pairs_n": len(placed) * (len(placed) - 1) // 2,
        "interior_disjoint_failures": interior_disjoint_failures,
        "internal_face_failures": internal_face_failures,
        "face_matches": face_matches,
        "boundary_co_tilings": boundary_co_tilings,
    }


def main() -> None:
    prototypes = {
        L: load_candidate(_DANZER_DIR / f"{L}.json") for L in "ABCK"
    }
    diss_data = json.loads(
        (_DANZER_DIR / "paolini_dissection.json").read_text()
    )
    by_parent: dict[str, list] = defaultdict(list)
    for ch in diss_data["children"]:
        by_parent[ch["parent"]].append(ch)

    print("=" * 70)
    print("Sibling-face-to-face audit — Paolini Danzer ABCK rule")
    print("=" * 70)

    total_failures = 0
    total_face_matches = 0
    total_interior_failures = 0

    for parent_letter in "ABCK":
        result = audit_parent(
            parent_letter, by_parent[parent_letter], prototypes,
        )
        n_children = len(result["placed"])
        print(
            f"\nσ({parent_letter}): {n_children} children, "
            f"{result['pairs_n']} sibling pairs"
        )
        print(
            f"  interior-disjoint pass: "
            f"{len(result['interior_disjoint_failures'])} failures"
        )
        if result["interior_disjoint_failures"]:
            for (i, j) in result["interior_disjoint_failures"]:
                print(
                    f"    [FAIL] children {i} ({result['child_types'][i]}) "
                    f"and {j} ({result['child_types'][j]}) overlap"
                )
        print(
            f"  internal face matches (true sibling face-to-face): "
            f"{result['face_matches']}"
        )
        print(
            f"  τP-boundary co-tilings (expected, not violations): "
            f"{result['boundary_co_tilings']}"
        )
        print(
            f"  internal face failures (true sibling-face-to-face "
            f"violations): {len(result['internal_face_failures'])}"
        )
        for fail in result["internal_face_failures"]:
            print(
                f"    [FAIL] child {fail['child_i']} "
                f"({fail['child_i_type']}) face {fail['face_i']} vs "
                f"child {fail['child_j']} ({fail['child_j_type']}) "
                f"face {fail['face_j']}: coplanar, internal, vertices differ"
            )

        total_failures += len(result["internal_face_failures"])
        total_interior_failures += len(result["interior_disjoint_failures"])
        total_face_matches += result["face_matches"]

    print()
    print("=" * 70)
    print(
        f"Summary: {total_face_matches} face-to-face coincidences, "
        f"{total_failures} face violations, "
        f"{total_interior_failures} interior overlaps"
    )
    if total_failures == 0 and total_interior_failures == 0:
        print(
            "VERDICT: sibling-face-to-face HOLDS for all sibling pairs "
            "in σ(A), σ(B), σ(C), σ(K)."
        )
    else:
        print(
            "VERDICT: sibling-face-to-face FAILS — see flagged pairs above."
        )
    print("=" * 70)


if __name__ == "__main__":
    main()
