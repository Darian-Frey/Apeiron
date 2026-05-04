"""Exhaustive n=3 PF=φ³ sweep at max_entry=3 (Q13c authorisation).

Per Claude (web) Q13c (2026-05-04): run n=3 max_entry=3 sweep with
the existing 9-class H₃ taxonomy (do not wait on the 15-class
extension — that's a separate follow-up). Document candidates
whose left-eigenvector ratios don't match any 9-class shape triple
as ``SKIPPED_TAXONOMY_GAP``. Lift k_max from 5 to 7 to reach the
SKIPPED-large-k cases.

Inputs:
- max_entry=3 → 209 primitive matrices at PF=φ³ (per Q9b).
- 9 H₃-compatible tetrahedral similarity classes from Paolini ABCK
  pool (per Q9a).
- k_max=7 (lifted from 5 per Q11a / Q12 follow-up).
- 30s/triple budget per realisation call.

Outputs:
- Per-candidate summary with verdict in:
    Realised, NoRealisation (all-triples or some-triples),
    Inconclusive, SKIPPED_LARGE_K, SKIPPED_TAXONOMY_GAP.

Usage:
    .venv/bin/python -m scripts.sweep_n3_phi3_max_entry_3 \\
        2>&1 | tee notebooks/n3_phi3_max_entry_3_$(date +%Y-%m-%d).log
"""

from __future__ import annotations

import os
import time
from typing import Iterator

import numpy as np

from apeiron.polyhedron import Polyhedron
from apeiron.track_b import (
    Inconclusive,
    Realised,
    enumerate_primitive_matrices,
    pf_eigenvector_in_zphi,
    realise,
)
from apeiron.track_b.taxonomy import build_h3_tetrahedra
from apeiron.zphi import ZPhi


PF_TARGET = ZPhi(1, 2)  # φ³
MAX_ENTRY = int(os.environ.get("APEIRON_MAX_ENTRY", "3"))
N_PROTOTILES = 3
K_MAX = int(os.environ.get("APEIRON_K_MAX", "7"))
PER_TRIPLE_BUDGET_SECONDS = float(
    os.environ.get("APEIRON_BUDGET", "30.0")
)


def _tetrahedron_volume(p: Polyhedron) -> ZPhi:
    """Six times the (signed) volume of a tetrahedron with ZPhi
    vertices. Returns absolute value."""
    v0, v1, v2, v3 = p.vertices
    a = v1 - v0
    b = v2 - v0
    c = v3 - v0
    cross = type(a)(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x,
    )
    vol_x6 = a.x * cross.x + a.y * cross.y + a.z * cross.z  # wrong
    # Recompute properly: 6V = |a · (b × c)|
    bc = type(a)(
        b.y * c.z - b.z * c.y,
        b.z * c.x - b.x * c.z,
        b.x * c.y - b.y * c.x,
    )
    vol_x6 = a.x * bc.x + a.y * bc.y + a.z * bc.z
    if vol_x6._sign() < 0:
        vol_x6 = -vol_x6
    return vol_x6


def _matching_shape_triples(
    eigvec: tuple[ZPhi, ZPhi, ZPhi],
    h3_pool: list[Polyhedron],
) -> Iterator[tuple[int, int, int]]:
    """Yield (i, j, k) — indices into h3_pool whose volume ratio
    matches the left eigenvector ratio.

    Match: volumes (v_i, v_j, v_k) proportional to
    (eigvec_0, eigvec_1, eigvec_2) up to a ZPhi scalar.
    Tested via cross-multiply: v_i · eigvec[1] = v_j · eigvec[0]
    and v_j · eigvec[2] = v_k · eigvec[1].
    """
    e0, e1, e2 = eigvec
    pool_vols = [_tetrahedron_volume(p) for p in h3_pool]
    n_pool = len(h3_pool)
    for i in range(n_pool):
        for j in range(n_pool):
            if i == j:
                continue
            if pool_vols[i] * e1 != pool_vols[j] * e0:
                continue
            for k in range(n_pool):
                if k in (i, j):
                    continue
                if pool_vols[j] * e2 != pool_vols[k] * e1:
                    continue
                yield (i, j, k)


def _candidate_summary(matrix: np.ndarray) -> str:
    """One-line representation of a 3×3 matrix."""
    rows = ["[" + " ".join(str(x) for x in r) + "]" for r in matrix]
    return ", ".join(rows)


def main() -> None:
    print(f"Exhaustive n=3 PF=φ³ sweep at max_entry={MAX_ENTRY}")
    print(f"k_max={K_MAX}, {PER_TRIPLE_BUDGET_SECONDS}s/triple budget")
    print(f"Per Q13c: 9-class H₃ taxonomy; unmatched marked SKIPPED_TAXONOMY_GAP")
    print("=" * 78)

    h3_pool = build_h3_tetrahedra()
    print(f"H₃-pool size: {len(h3_pool)}")

    candidates = list(enumerate_primitive_matrices(
        N_PROTOTILES, PF_TARGET, max_entry=MAX_ENTRY,
    ))
    print(f"Primitive candidates at max_entry={MAX_ENTRY}: {len(candidates)}")
    print("=" * 78)

    counts = {
        "Realised": 0,
        "NoRealisation_all_triples": 0,
        "NoRealisation_some_triples": 0,
        "Mixed": 0,
        "Inconclusive": 0,
        "SKIPPED_LARGE_K": 0,
        "SKIPPED_TAXONOMY_GAP": 0,
        "FILTER1_REJECT": 0,
    }
    realised_witnesses: list[tuple[int, np.ndarray, tuple[int, int, int]]] = []
    start = time.monotonic()

    for cand_idx, M in enumerate(candidates):
        cand_label = f"cand {cand_idx + 1}/{len(candidates)}"
        eigvec = pf_eigenvector_in_zphi(M.T, PF_TARGET)
        if eigvec is None:
            counts["FILTER1_REJECT"] += 1
            print(f"  [{cand_label}] FILTER1 reject (no positive ZPhi eigvec for M.T)")
            continue
        triples = list(_matching_shape_triples(eigvec, h3_pool))
        if not triples:
            counts["SKIPPED_TAXONOMY_GAP"] += 1
            print(
                f"  [{cand_label}] SKIPPED_TAXONOMY_GAP — left "
                f"eigvec ratio not in 9-class volume set; "
                f"{_candidate_summary(M)}"
            )
            continue

        # Run all matching triples for this candidate.
        per_triple = []
        any_realised = False
        any_inconclusive = False
        any_skipped_k = False
        for triple in triples:
            shapes = tuple(h3_pool[t] for t in triple)
            r = realise(
                M, PF_TARGET, prototile_shapes=shapes,
                max_search_seconds=PER_TRIPLE_BUDGET_SECONDS,
                k_max=K_MAX,
            )
            per_triple.append((triple, r))
            if isinstance(r, Realised):
                any_realised = True
                realised_witnesses.append((cand_idx, M, triple))
            elif isinstance(r, Inconclusive):
                if "k_max" in r.reason:
                    any_skipped_k = True
                else:
                    any_inconclusive = True

        # Classify the candidate.
        if any_realised:
            counts["Realised"] += 1
            verdict = "REALISED"
        elif any_inconclusive:
            counts["Inconclusive"] += 1
            verdict = "Inconclusive"
        elif any_skipped_k and not any_inconclusive:
            counts["SKIPPED_LARGE_K"] += 1
            verdict = f"SKIPPED_LARGE_K (k>{K_MAX})"
        else:
            # All triples returned NoRealisation.
            counts["NoRealisation_all_triples"] += 1
            verdict = f"NoRealisation (all {len(per_triple)} triples)"
        print(
            f"  [{cand_label}] {len(triples):3d} triples → {verdict}"
        )

    elapsed = time.monotonic() - start
    print("=" * 78)
    print(f"\nPer-candidate summary (elapsed: {elapsed:.1f}s):")
    for k, v in counts.items():
        print(f"  {k:32s}: {v:5d}")
    print(f"  total candidates              : {len(candidates):5d}")

    if realised_witnesses:
        print()
        print("REALISED WITNESSES:")
        for cand_idx, M, triple in realised_witnesses:
            print(f"  candidate {cand_idx}: {_candidate_summary(M)} — triple {triple}")
    else:
        print()
        print("No Realised witnesses across the entire sweep.")


if __name__ == "__main__":
    main()
