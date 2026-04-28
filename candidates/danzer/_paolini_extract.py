"""Extract per-child (rotation, translation) pairs from Paolini's
Danzer ABCK POV-Ray source.

Local mirror at /tmp/danzer_pov/danzer_transformations.inc (25 transforms,
one per child of σ(A), σ(B), σ(C), σ(K)). This script transcribes each
transform as a Python expression mirroring the POV-Ray source, evaluates
it on the appropriate prototile's vertices using float Rodrigues, and
recovers (R_×2, t_×2) — the matrix and translation in Apeiron's
denominator-2 storage convention.

This is a *combinatorial oracle* in the CLAUDE.md memory sense: floats
are used to reach the rotation matrix; the result is then snapped to
ZPhi and verified to be in I_h (i.e., proper rotation in
``apeiron.symmetry.ICOSAHEDRAL`` or improper of the form -2g for some
g ∈ ICOSAHEDRAL). Any transform that fails to land in I_h is a bug —
either in the transcription or in the snap tolerance.

Output: ``candidates/danzer/paolini_dissection.json`` with per-child
records ``{"parent": "A"|"B"|"C"|"K", "child_index": int, "child_type":
"A"|"B"|"C"|"K", "is_proper": bool, "icosahedral_index": int,
"translation_x2": [[a,b], [a,b], [a,b]]}``.

NOT a runtime dependency — meant to be run once to produce the JSON,
then committed.

Status: pending Q4 ruling from Claude (web). If Q4a chooses Paolini,
this script's output becomes the canonical encoding. If Q4a chooses
Koca, this serves as the cross-validation (per-child position-match
under a global isometry, per Q4b).
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

PHI = (1.0 + math.sqrt(5.0)) / 2.0
TOL = 1e-9


def real(a: int, b: int) -> float:
    return a + b * PHI


# --- Paolini's pt0..pt9 in real coordinates ---

pt = {
    0: (real(0, 0),  real(0, 0),  real(0, 0)),
    1: (real(1, 2),  real(0, 0),  real(1, 1)),    # τ³, 0, τ²
    2: (real(1, 1),  real(1, 1),  real(1, 1)),    # τ²(1,1,1)
    3: (real(1, 1),  real(1, 0),  real(0, 0)),    # τ², 1, 0
    4: (real(1, 1),  real(0, 1),  real(1, 0)),    # τ², τ, 1
    5: (real(0, -1), real(0, 0),  real(1, 0)),    # -τ, 0, 1
    6: (real(0, 0),  real(1, 1),  real(1, 0)),    # 0, τ², 1
    7: (real(-1, 0), real(0, 1),  real(0, 0)),    # -1, τ, 0
    8: (real(0, 1),  real(0, 1),  real(0, 1)),    # τ(1,1,1)
    9: (-0.5,        -0.5 + 0.5*PHI,  0.5*PHI),   # ½(-1, 1/τ, τ); 1/τ = τ-1
}

# Per-prototile vertex aliases (Paolini convention).
prototile_vertices = {
    "A": [pt[0], pt[1], pt[2], pt[3]],   # pt[1A,2A,3A,4A] = pt[0,1,2,3]
    "B": [pt[0], pt[1], pt[2], pt[4]],   # pt[1B,2B,3B,4B] = pt[0,1,2,4]
    "C": [pt[0], pt[5], pt[2], pt[6]],   # pt[1C,2C,3C,4C] = pt[0,5,2,6]
    "K": [pt[0], pt[7], pt[8], pt[9]],   # pt[1K,2K,3K,4K] = pt[0,7,8,9]
}

prototile_index = {"A": 0, "B": 1, "C": 2, "K": 3}


# --- Vector / matrix primitives ---


def vsub(u, v): return tuple(u[i]-v[i] for i in range(3))
def vadd(u, v): return tuple(u[i]+v[i] for i in range(3))
def vscale(s, v): return tuple(s*v[i] for i in range(3))
def vneg(v): return tuple(-v[i] for i in range(3))
def vdot(u, v): return sum(u[i]*v[i] for i in range(3))
def vcross(u, v):
    return (u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0])
def vnorm(v): return math.sqrt(vdot(v, v))


def matvec(M, v):
    return tuple(sum(M[i][j]*v[j] for j in range(3)) for i in range(3))


def matmat(A, B):
    return tuple(tuple(sum(A[i][k]*B[k][j] for k in range(3)) for j in range(3)) for i in range(3))


I3 = ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))


def axis_rotate(axis, deg):
    """Rodrigues rotation matrix; axis need not be unit."""
    n = vnorm(axis)
    if n < 1e-14:
        return I3
    k = (axis[0]/n, axis[1]/n, axis[2]/n)
    th = math.radians(deg)
    c = math.cos(th); s = math.sin(th)
    K = ((0.0, -k[2], k[1]), (k[2], 0.0, -k[0]), (-k[1], k[0], 0.0))
    KK = matmat(K, K)
    return tuple(tuple(I3[i][j] + s*K[i][j] + (1-c)*KK[i][j] for j in range(3)) for i in range(3))


def reorient(v_from, v_to):
    """Rotation that maps direction(v_from) to direction(v_to)."""
    nf = vnorm(v_from); nt = vnorm(v_to)
    if nf < 1e-14 or nt < 1e-14:
        return I3
    f = (v_from[0]/nf, v_from[1]/nf, v_from[2]/nf)
    t = (v_to[0]/nt, v_to[1]/nt, v_to[2]/nt)
    cos_th = max(-1.0, min(1.0, vdot(f, t)))
    if cos_th > 1 - 1e-14:
        return I3
    if cos_th < -1 + 1e-14:
        if abs(f[0]) < 0.9: perp = vcross(f, (1.0, 0.0, 0.0))
        else: perp = vcross(f, (0.0, 1.0, 0.0))
        return axis_rotate(perp, 180)
    axis = vcross(f, t)
    return axis_rotate(axis, math.degrees(math.acos(cos_th)))


def alllinea_triangoli_trans(A, B, C, X, Y, Z):
    """Paolini's AllineaTriangoliTrans macro with AVOIDINSTABILITY=0
    (the 40°-randrot kludge disabled — it's a numerical hack, not a
    clean isometry).

    Computes the (R, t) such that: R(A) + t = X, R(B) + t = Y,
    R(C) + t = Z. R is the unique proper isometry mapping triangle
    ABC to triangle XYZ (same orientation).
    """
    N1 = vcross(vsub(B, A), vsub(C, A))
    N2 = vcross(vsub(Y, X), vsub(Z, X))
    R1 = reorient(N1, N2)
    rotAB = matvec(R1, vsub(B, A))
    R2 = reorient(rotAB, vsub(Y, X))
    R = matmat(R2, R1)
    t = vsub(X, matvec(R, A))
    return R, t


# --- Composable transforms (per POV-Ray transform { ... } block) ---


def affine_compose(T1, T2):
    """Compose two affine transforms (R1, t1), (R2, t2) — applied as
    v ↦ R2(R1·v + t1) + t2 = (R2·R1)·v + (R2·t1 + t2). Matches POV-Ray's
    left-to-right transform-block ordering."""
    R1, t1 = T1
    R2, t2 = T2
    return matmat(R2, R1), vadd(matvec(R2, t1), t2)


def trf_identity():
    return I3, (0.0, 0.0, 0.0)


def trf_scale_neg1():
    """scale<-1,-1,-1> — point inversion."""
    return ((-1.0, 0.0, 0.0), (0.0, -1.0, 0.0), (0.0, 0.0, -1.0)), (0.0, 0.0, 0.0)


def trf_translate(v):
    return I3, v


def trf_axis_rotate(axis, deg):
    return axis_rotate(axis, deg), (0.0, 0.0, 0.0)


def trf_alllinea(A, B, C, X, Y, Z):
    return alllinea_triangoli_trans(A, B, C, X, Y, Z)


# Composite tileK transforms (defined in danzer_tiles.inc):
# tile8Kcenter = pt9; tileKface1 = pt9; tileKblue = pt9-pt8; tileKface2 = pt9-pt7.

def trf_tileK90():
    return _seq([
        trf_translate(vneg(pt[9])),
        trf_axis_rotate(pt[9], 180),
        trf_scale_neg1(),
        trf_translate(pt[9]),
    ])

def trf_tileK180():
    return _seq([
        trf_translate(vneg(pt[9])),
        trf_axis_rotate(vsub(pt[9], pt[8]), 180),
        trf_translate(pt[9]),
    ])

def trf_tileK270():
    return _seq([
        trf_translate(vneg(pt[9])),
        trf_axis_rotate(vsub(pt[9], pt[7]), 180),
        trf_scale_neg1(),
        trf_translate(pt[9]),
    ])

def trf_tileKneg():
    return _seq([
        trf_translate(vneg(pt[9])),
        trf_scale_neg1(),
        trf_translate(pt[9]),
    ])


def _seq(steps):
    """Compose a sequence of transforms (left-to-right order)."""
    out = trf_identity()
    for step in steps:
        out = affine_compose(out, step)
    return out


# --- The 25 Paolini transforms ---


def TR_A_1():
    return _seq([
        trf_scale_neg1(),
        trf_alllinea(vneg(pt[1]), vneg(pt[0]), vneg(pt[2]),
                     vscale(PHI, pt[2]), vscale(PHI, pt[1]), vscale(PHI, pt[4])),
    ])

def TR_A_2():
    pt1, pt0, pt2 = pt[1], pt[0], pt[2]
    return _seq([
        trf_axis_rotate(vcross(vsub(pt1, pt0), vsub(pt2, pt0)), 180),
        trf_translate(vadd(pt0, vscale(PHI, vsub(pt1, pt0)))),
        trf_axis_rotate(vsub(pt1, pt0), 144),  # 2*72
    ])

def TR_A_3():
    pt1B, pt4B, pt2B = pt[0], pt[4], pt[1]   # B's vertex aliases
    pt2A, pt4A, pt1A = pt[1], pt[3], pt[0]
    inv_tau = PHI - 1   # 1/τ = τ - 1
    return _seq([
        trf_scale_neg1(),
        trf_alllinea(vneg(pt1B), vneg(pt4B), vneg(pt2B),
                     vscale(PHI, pt2A),
                     vadd(pt4A, vscale(inv_tau, pt2A)),
                     vadd(pt1A, vscale(inv_tau, pt2A))),
    ])

def TR_A_4():
    return _seq([trf_scale_neg1(), trf_translate(vsub(pt[2], pt[0]))])

def TR_A_5():
    pt3C, pt4C, pt1C = pt[2], pt[6], pt[0]
    pt1A, pt2A, pt4A = pt[0], pt[1], pt[3]
    inv_tau = PHI - 1
    return trf_alllinea(pt3C, pt4C, pt1C,
                        vscale(PHI, pt1A),
                        vadd(pt1A, vscale(inv_tau, pt2A)),
                        vscale(PHI, pt4A))

def TR_A_6():
    return _seq([
        trf_tileK270(),
        trf_axis_rotate(vsub(pt[8], pt[0]), 120),  # pt3K-pt1K = pt8-pt0
        trf_translate(vscale(PHI, vsub(pt[8], pt[0]))),
    ])

def TR_A_7():
    return _seq([
        trf_tileK180(),
        trf_axis_rotate(vsub(pt[8], pt[0]), 120),
        trf_translate(vscale(PHI, vsub(pt[8], pt[0]))),
    ])

def TR_A_8():
    return _seq([
        trf_tileK90(),
        trf_tileKneg(),
        trf_axis_rotate(vsub(pt[8], pt[0]), 120),
        trf_translate(vscale(PHI, vsub(pt[8], pt[0]))),
    ])

def TR_A_9():
    return _seq([
        trf_tileKneg(),
        trf_axis_rotate(vsub(pt[8], pt[0]), 120),
        trf_translate(vscale(PHI, vsub(pt[8], pt[0]))),
    ])

def TR_A_10():
    pt2K, pt3K, pt1K = pt[7], pt[8], pt[0]
    return _seq([
        trf_translate(vneg(pt2K)),
        trf_axis_rotate(vsub(pt3K, pt2K), 216),  # 3*72
        trf_translate(pt2K),
        trf_tileK90(),
        trf_tileKneg(),
        trf_axis_rotate(vsub(pt3K, pt1K), 120),
        trf_translate(vscale(PHI, vsub(pt3K, pt1K))),
    ])

def TR_A_11():
    pt2K, pt3K, pt1K = pt[7], pt[8], pt[0]
    return _seq([
        trf_translate(vneg(pt2K)),
        trf_axis_rotate(vsub(pt3K, pt2K), -216),
        trf_translate(pt2K),
        trf_tileKneg(),
        trf_axis_rotate(vsub(pt3K, pt1K), 120),
        trf_translate(vscale(PHI, vsub(pt3K, pt1K))),
    ])

def TR_B_1(): return TR_A_1()
def TR_B_2(): return TR_A_2()

def TR_B_3():
    return _seq([trf_scale_neg1(), trf_translate(vsub(pt[2], pt[0]))])

def TR_B_4(): return TR_A_6()
def TR_B_5(): return TR_A_7()
def TR_B_6(): return TR_A_8()
def TR_B_7(): return TR_A_9()


def TR_C_1():
    return _seq([
        trf_scale_neg1(),
        trf_alllinea(vneg(pt[0]), vneg(pt[1]), vneg(pt[2]),
                     vscale(PHI, pt[6]), vscale(PHI, pt[5]), vscale(PHI, pt[0])),
    ])

def TR_C_2():
    pt1C, pt2C, pt3C, pt4C = pt[0], pt[5], pt[2], pt[6]
    inv_tau = PHI - 1
    return _seq([
        trf_scale_neg1(),
        trf_alllinea(vneg(pt1C), vneg(pt2C), vneg(pt3C),
                     vscale(PHI, pt4C),
                     vadd(pt4C, vscale(inv_tau, pt3C)),
                     vscale(PHI, pt1C)),
    ])

def TR_C_3():
    pt3C, pt1C, pt4C, pt4A = pt[2], pt[0], pt[6], pt[3]
    rot1 = trf_axis_rotate(vsub(pt3C, pt1C), -60)
    # locnormal = pt4C - pt1C - dot(pt4A-pt1C, pt3C-pt1C)/dot(pt3C-pt1C,pt3C-pt1C)·(pt3C-pt1C)
    diff = vsub(pt3C, pt1C)
    proj_scalar = vdot(vsub(pt4A, pt1C), diff) / vdot(diff, diff)
    locnormal = vsub(vsub(pt4C, pt1C), vscale(proj_scalar, diff))
    return _seq([
        rot1,
        trf_axis_rotate(locnormal, 180),
        trf_translate(diff),
    ])

def TR_C_4():
    return trf_translate(vscale(PHI, vsub(pt[8], pt[0])))   # τ·(pt3K-pt1K)

def TR_C_5():
    return _seq([
        trf_tileK180(),
        trf_tileKneg(),
        trf_translate(vscale(PHI, vsub(pt[8], pt[0]))),
    ])


def TR_K_1():
    return _seq([trf_scale_neg1(), trf_translate(vsub(pt[2], pt[0]))])

def TR_K_2():
    return trf_alllinea(pt[7], pt[9], pt[8],
                        vscale(PHI, pt[0]), vscale(PHI, pt[9]), vscale(PHI, pt[7]))


# --- Dispatch table from danzer_recursion.inc inline comments ---
#   parent prototile -> list of (transform_func, child_prototile, is_proper).
#   is_proper inferred from the inline "tileX (negative)" comment;
#   absence of "(negative)" means proper.

DISPATCH = {
    "A": [
        (TR_A_1,  "B", False),
        (TR_A_2,  "B", True),
        (TR_A_3,  "B", False),
        (TR_A_4,  "C", False),
        (TR_A_5,  "C", True),
        (TR_A_6,  "K", False),
        (TR_A_7,  "K", True),
        (TR_A_8,  "K", True),
        (TR_A_9,  "K", False),
        (TR_A_10, "K", True),
        (TR_A_11, "K", True),
    ],
    "B": [
        (TR_B_1,  "B", False),
        (TR_B_2,  "B", True),
        (TR_B_3,  "C", False),
        (TR_B_4,  "K", False),
        (TR_B_5,  "K", True),
        (TR_B_6,  "K", True),
        (TR_B_7,  "K", False),
    ],
    "C": [
        (TR_C_1,  "A", False),
        (TR_C_2,  "C", False),
        (TR_C_3,  "C", True),
        (TR_C_4,  "K", True),
        (TR_C_5,  "K", False),
    ],
    "K": [
        (TR_K_1,  "B", None),  # determined below
        (TR_K_2,  "K", None),
    ],
}


# --- ZPhi snap ---


def zphi_snap(x: float) -> tuple[int, int]:
    best = None; best_err = float("inf")
    for b in range(-300, 301):
        a = round(x - b * PHI)
        err = abs((a + b * PHI) - x)
        if err < best_err:
            best_err, best = err, (a, b)
    if best_err > TOL:
        raise ValueError(f"x={x} doesn't snap to ZPhi: err={best_err:.3e}")
    return best


def stored_snap(x): return zphi_snap(2 * x)


# --- Main loop ---

def _main():
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from apeiron.symmetry import ICOSAHEDRAL, Mat3, Rotation, Vec3
    from apeiron.zphi import ZPhi
    import numpy as np

    icos_set = list(ICOSAHEDRAL)
    icos_lookup = {g: i for i, g in enumerate(icos_set)}

    def to_rotation(R_st):
        m = Mat3(
            Vec3(ZPhi(*R_st[0][0]), ZPhi(*R_st[0][1]), ZPhi(*R_st[0][2])),
            Vec3(ZPhi(*R_st[1][0]), ZPhi(*R_st[1][1]), ZPhi(*R_st[1][2])),
            Vec3(ZPhi(*R_st[2][0]), ZPhi(*R_st[2][1]), ZPhi(*R_st[2][2])),
        )
        return Rotation(m)

    records: list[dict] = []
    for parent, dispatch in DISPATCH.items():
        verts = prototile_vertices[parent]
        for child_idx, (trf_fn, child_type, is_proper_hint) in enumerate(dispatch):
            R_real, t_real = trf_fn()
            placed = [vadd(matvec(R_real, v), t_real) for v in verts]
            # Recover (R, t) by least-squares from placed vertices.
            P = np.array([vsub(verts[i], verts[0]) for i in (1, 2, 3)]).T
            Q = np.array([vsub(placed[i], placed[0]) for i in (1, 2, 3)]).T
            R_recovered = (Q @ np.linalg.inv(P))
            t_recovered = tuple(placed[0][i] - sum(R_recovered[i][j]*verts[0][j] for j in range(3)) for i in range(3))
            R_st = tuple(tuple(stored_snap(R_recovered[i][j]) for j in range(3)) for i in range(3))
            t_st = tuple(stored_snap(c) for c in t_recovered)
            # Determinant: in stored form, det(2R) = 8·det(R).
            det_real = float(np.linalg.det(R_recovered))
            is_proper = det_real > 0
            if abs(abs(det_real) - 1.0) > 1e-9:
                raise AssertionError(
                    f"{parent} child {child_idx} ({trf_fn.__name__}): "
                    f"determinant {det_real} is not ±1 (not an isometry)"
                )
            # Look up: if proper, R_st should match some ICOSAHEDRAL element.
            #          if improper, -R_st (entry-wise negation) should match.
            if is_proper:
                cand = to_rotation(R_st)
                lookup = icos_lookup.get(cand)
            else:
                R_st_neg = tuple(tuple((-R_st[i][j][0], -R_st[i][j][1]) for j in range(3)) for i in range(3))
                cand = to_rotation(R_st_neg)
                lookup = icos_lookup.get(cand)
            if lookup is None:
                raise AssertionError(
                    f"{parent} child {child_idx} ({trf_fn.__name__}): "
                    f"R not in I_h. R_real:\n{R_recovered}"
                )
            records.append({
                "parent": parent,
                "child_index": child_idx,
                "child_type": child_type,
                "is_proper": is_proper,
                "icosahedral_index": lookup,
                "translation_x2": [list(c) for c in t_st],
                "_transform_name": trf_fn.__name__,
                "_paolini_negative_hint": (
                    "negative" if is_proper_hint is False else
                    "positive" if is_proper_hint is True else "?"
                ),
            })
            print(
                f"  {parent} child {child_idx} ({trf_fn.__name__:8s}) -> "
                f"{child_type} | {'proper' if is_proper else 'improper'} | "
                f"icos[{lookup:2d}] | t={t_st}"
            )

    out_path = Path(__file__).resolve().parent / "paolini_dissection.json"
    out_path.write_text(json.dumps({
        "source": "Paolini POV-Ray (svn.dmf.unicatt.it/projects/animations/danzer/trunk)",
        "extractor": "candidates/danzer/_paolini_extract.py",
        "convention": "translation_x2 stores 2·t in ZPhi (denominator-2 convention, CLAUDE.md §3.2). icosahedral_index references apeiron.symmetry.ICOSAHEDRAL[i]; if is_proper is false, the rotation is ImproperRotation(ICOSAHEDRAL[icosahedral_index]).",
        "children": records,
    }, indent=2) + "\n")
    print(f"\nWrote {out_path} ({len(records)} children).")
    print(f"  σ(A): {sum(1 for r in records if r['parent']=='A')} children (expected 11)")
    print(f"  σ(B): {sum(1 for r in records if r['parent']=='B')} children (expected 7)")
    print(f"  σ(C): {sum(1 for r in records if r['parent']=='C')} children (expected 5)")
    print(f"  σ(K): {sum(1 for r in records if r['parent']=='K')} children (expected 2)")


if __name__ == "__main__":
    _main()
