"""Algebraic viz of Track B's 5 candidate substitution matrices at
n=2, PF=φ³ (= ZPhi(1, 2)).

These candidates have no geometric realisation yet — that's the
output of the ``apeiron/track_b/realisation.py`` CSP, pending
Q8 ruling from Claude (web). What we *do* have for each candidate:

- The 2×2 matrix entries.
- The right Perron–Frobenius eigenvector in pure ℤ[φ] (= prototile
  *volume ratio* in the eventual realisation, per Q7c filter 1).
- Primitivity dynamics (smallest power M^k with all entries
  strictly positive).
- The substitution graph — arrows P_i → m·P_j with multiplicities.

This notebook renders all five candidates as a single interactive
HTML dashboard for at-a-glance comparison.

Run from the repo root::

    .venv/bin/python notebooks/track_b_phi_cubed_candidates.py
"""

from __future__ import annotations

import math
import subprocess
import sys
from pathlib import Path

import numpy as np

from apeiron.substitution import is_primitive
from apeiron.track_b import enumerate_primitive_matrices, prefilter
from apeiron.zphi import ZPhi

REPO_ROOT = Path(__file__).resolve().parent.parent
OUTPUT_HTML = REPO_ROOT / "notebooks" / "track_b_phi_cubed.html"

PHI_CUBED = ZPhi(1, 2)
PHI_FLOAT = (1.0 + math.sqrt(5.0)) / 2.0


def _zphi_to_float(z: ZPhi) -> float:
    return z.a + z.b * PHI_FLOAT


def _smallest_primitive_power(matrix: np.ndarray, cap: int = 20) -> int | None:
    """Smallest k ≥ 1 such that M^k is strictly positive entry-wise."""
    p = matrix.copy()
    for k in range(1, cap + 1):
        if np.all(p > 0):
            return k
        p = p @ matrix
    return None


def _matrix_to_html_table(matrix: np.ndarray) -> str:
    rows = []
    for i in range(matrix.shape[0]):
        cells = "".join(
            f'<td style="padding: 4px 12px; border: 1px solid #ccc; '
            f'text-align: center;">{int(matrix[i, j])}</td>'
            for j in range(matrix.shape[1])
        )
        rows.append(f"<tr>{cells}</tr>")
    return (
        '<table style="border-collapse: collapse; '
        'font-family: monospace; font-size: 14px;">'
        f'{"".join(rows)}</table>'
    )


def _zphi_eigenvector_to_str(eigenvector: tuple[ZPhi, ...]) -> str:
    parts = []
    for z in eigenvector:
        if z.b == 0:
            parts.append(f"{z.a}")
        elif z.a == 0:
            if z.b == 1:
                parts.append("φ")
            elif z.b == -1:
                parts.append("-φ")
            else:
                parts.append(f"{z.b}φ")
        else:
            sign = "+" if z.b > 0 else "-"
            mag = abs(z.b)
            tail = "φ" if mag == 1 else f"{mag}φ"
            parts.append(f"{z.a} {sign} {tail}")
    return f"({', '.join(parts)})"


def _build_candidate_panel(
    idx: int, matrix: np.ndarray, eigenvector: tuple[ZPhi, ...],
) -> str:
    """One HTML panel summarising a single candidate."""
    a, b = int(matrix[0, 0]), int(matrix[0, 1])
    c, d = int(matrix[1, 0]), int(matrix[1, 1])
    trace = a + d
    det = a * d - b * c
    prim_k = _smallest_primitive_power(matrix)
    eigenvector_str = _zphi_eigenvector_to_str(eigenvector)
    v0_real = _zphi_to_float(eigenvector[0])
    v1_real = _zphi_to_float(eigenvector[1])
    ratio = v1_real / v0_real
    # σ(P_0) gets `a` copies of P_0 and `c` copies of P_1.
    # σ(P_1) gets `b` copies of P_0 and `d` copies of P_1.
    sigma_lines = []
    for col, label in enumerate(("P_0", "P_1")):
        children: list[str] = []
        for row, child_label in enumerate(("P_0", "P_1")):
            n = int(matrix[row, col])
            if n == 0:
                continue
            if n == 1:
                children.append(child_label)
            else:
                children.append(f"{n}·{child_label}")
        sigma_lines.append(
            f"σ({label}) = " + (" + ".join(children) if children else "(empty)")
        )
    return f"""
<div style="border: 1px solid #888; padding: 14px; margin: 10px;
            border-radius: 8px; background: #fafafa; min-width: 320px;">
    <h3 style="margin: 0 0 8px 0;">Candidate {idx} — M = [[{a}, {b}], [{c}, {d}]]</h3>
    <div style="display: flex; gap: 24px; align-items: flex-start;">
        <div>
            <div style="font-weight: bold; margin-bottom: 4px;">Matrix M</div>
            {_matrix_to_html_table(matrix)}
        </div>
        <div style="font-family: monospace; font-size: 13px;
                    line-height: 1.6;">
            <div><b>trace</b> = {trace}</div>
            <div><b>det</b> = {det}</div>
            <div><b>PF eigenvalue</b> = φ³ = 1 + 2φ ≈ {PHI_FLOAT**3:.4f}</div>
            <div><b>primitivity power</b> k = {prim_k}</div>
            <div><b>PF eigenvector (ℤ[φ])</b> = {eigenvector_str}</div>
            <div><b>volume ratio v_1/v_0</b> ≈ {ratio:.4f}</div>
        </div>
        <div style="font-family: monospace; font-size: 13px;
                    line-height: 1.6;">
            <div style="font-weight: bold;">Substitution rule</div>
            <div>{sigma_lines[0]}</div>
            <div>{sigma_lines[1]}</div>
        </div>
    </div>
</div>
"""


def _build_volume_ratio_chart(
    candidates: list[tuple[np.ndarray, tuple[ZPhi, ...]]],
) -> str:
    """Plotly bar chart of v_0 vs v_1 (normalised so v_0 + v_1 = 1) for
    each candidate, showing prototile volume fractions."""
    import plotly.graph_objects as go

    labels = [f"M_{i + 1}" for i in range(len(candidates))]
    v0 = []
    v1 = []
    for _, ev in candidates:
        a = _zphi_to_float(ev[0])
        b = _zphi_to_float(ev[1])
        total = a + b
        v0.append(a / total)
        v1.append(b / total)

    fig = go.Figure(data=[
        go.Bar(name="P_0 vol fraction", x=labels, y=v0,
               marker_color="#4a90e2"),
        go.Bar(name="P_1 vol fraction", x=labels, y=v1,
               marker_color="#e8743b"),
    ])
    fig.update_layout(
        barmode="stack",
        title="Prototile volume fractions per candidate (normalised)",
        yaxis=dict(title="fraction of total volume"),
        height=400,
        font=dict(family="monospace", size=12),
    )
    return fig.to_html(full_html=False, include_plotlyjs=False)


def _build_substitution_graph_chart(
    candidates: list[tuple[np.ndarray, tuple[ZPhi, ...]]],
) -> str:
    """Plotly node-link diagram of each candidate's substitution graph
    (P_0 ↔ P_1 with multiplicity edges)."""
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    n_cand = len(candidates)
    fig = make_subplots(
        rows=1, cols=n_cand,
        subplot_titles=[f"M_{i + 1}" for i in range(n_cand)],
    )

    for col_idx, (matrix, _) in enumerate(candidates, start=1):
        # Layout: P_0 at (0, 0), P_1 at (1, 0).
        node_x = [0.0, 1.0]
        node_y = [0.0, 0.0]
        node_text = ["P_0", "P_1"]
        # Self-loops on P_0 (count = M[0,0]) and P_1 (M[1,1]),
        # cross-edges P_0 → P_1 (M[1,0]) and P_1 → P_0 (M[0,1]).
        edges_x: list[float | None] = []
        edges_y: list[float | None] = []
        annotations_text = []
        # Cross edges.
        m_01 = int(matrix[0, 1])  # P_1 → P_0 (rows = output)
        m_10 = int(matrix[1, 0])  # P_0 → P_1
        if m_10 > 0:
            edges_x.extend([0.05, 0.95, None])
            edges_y.extend([0.05, 0.05, None])
            annotations_text.append((0.5, 0.12, f"{m_10}"))
        if m_01 > 0:
            edges_x.extend([0.95, 0.05, None])
            edges_y.extend([-0.05, -0.05, None])
            annotations_text.append((0.5, -0.12, f"{m_01}"))
        # Self-loops as text annotations only (drawing actual loops in
        # plotly is fiddly; the multiplicity number is the important
        # info).
        m_00 = int(matrix[0, 0])
        m_11 = int(matrix[1, 1])
        if m_00 > 0:
            annotations_text.append((-0.18, 0.0, f"⟲{m_00}"))
        if m_11 > 0:
            annotations_text.append((1.18, 0.0, f"⟲{m_11}"))

        fig.add_trace(go.Scatter(
            x=edges_x, y=edges_y, mode="lines",
            line=dict(color="#888", width=2),
            hoverinfo="none", showlegend=False,
        ), row=1, col=col_idx)
        fig.add_trace(go.Scatter(
            x=node_x, y=node_y, mode="markers+text",
            marker=dict(size=40, color=["#4a90e2", "#e8743b"]),
            text=node_text, textposition="middle center",
            textfont=dict(color="white", size=14, family="monospace"),
            hoverinfo="text", showlegend=False,
        ), row=1, col=col_idx)
        for x, y, text in annotations_text:
            fig.add_annotation(
                x=x, y=y, text=text, showarrow=False,
                xref=f"x{col_idx if col_idx > 1 else ''}",
                yref=f"y{col_idx if col_idx > 1 else ''}",
                font=dict(size=12, family="monospace"),
            )

    fig.update_xaxes(showticklabels=False, range=[-0.4, 1.4],
                     showgrid=False, zeroline=False)
    fig.update_yaxes(showticklabels=False, range=[-0.3, 0.3],
                     showgrid=False, zeroline=False)
    fig.update_layout(
        title="Substitution graphs (edges labeled with multiplicities)",
        height=300,
        font=dict(family="monospace"),
    )
    return fig.to_html(full_html=False, include_plotlyjs=False)


def main() -> int:
    print("Enumerating Track B candidates at n=2, PF=φ³, max_entry=5…")
    candidates: list[tuple[np.ndarray, tuple[ZPhi, ...]]] = []
    for matrix in enumerate_primitive_matrices(2, PHI_CUBED, max_entry=5):
        result = prefilter(matrix, PHI_CUBED)
        if not result.passes_all:
            continue
        assert result.eigenvector is not None
        candidates.append((matrix, result.eigenvector))

    print(f"Found {len(candidates)} candidates passing all prefilters.")
    for i, (matrix, ev) in enumerate(candidates, start=1):
        print(f"  M_{i} = {matrix.tolist()}  eigenvector = {_zphi_eigenvector_to_str(ev)}")

    # Sanity: every candidate primitive.
    for matrix, _ in candidates:
        assert is_primitive(matrix)

    panels = "\n".join(
        _build_candidate_panel(i + 1, m, ev)
        for i, (m, ev) in enumerate(candidates)
    )
    chart_volumes = _build_volume_ratio_chart(candidates)
    chart_graphs = _build_substitution_graph_chart(candidates)

    html = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>Track B candidates — n=2, PF=φ³</title>
<script src="https://cdn.plot.ly/plotly-2.32.0.min.js"></script>
<style>
    body {{
        font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
        max-width: 1400px; margin: 24px auto; padding: 0 24px;
        line-height: 1.6;
    }}
    h1 {{ margin-bottom: 4px; }}
    .subtitle {{ color: #666; margin-bottom: 20px; }}
    .panel-row {{ display: flex; flex-wrap: wrap; }}
</style>
</head>
<body>
    <h1>Track B candidates — n = 2, PF = φ³ = ZPhi(1, 2)</h1>
    <p class="subtitle">
        Algebraic survey of primitive 2×2 substitution matrices with
        Perron–Frobenius eigenvalue φ³ at <code>max_entry=5</code>,
        deduplicated under simultaneous row/column permutation.
        Per <code>apeiron.track_b.matrix_search</code> (commit
        <code>d994384</code>) and
        <code>apeiron.track_b.geometric_prefilter</code>
        (commit <code>a84593d</code>) — Q7c implementation.
        All {len(candidates)} candidates pass filter 1
        (vertex-class consistency, ZPhi positive eigenvector);
        filters 2 + 3 are documented stubs and don't eliminate.
    </p>
    <p>
        <b>What this is.</b> Each candidate is an abstract 2-prototile
        substitution rule. None has a known 3D polyhedral
        realisation yet — that's the output of the
        <code>realisation.py</code> CSP, pending Claude (web)'s
        Q8 ruling.
    </p>
    <div class="panel-row">
        {panels}
    </div>
    <h2>Volume ratios (from PF eigenvectors)</h2>
    <p>
        The Perron–Frobenius eigenvector of each substitution matrix
        gives the *volume ratio* of the two prototiles in any
        eventual realisation. Bars below are normalised so each
        candidate's bars sum to 1 (i.e., each shows the fraction of
        total tile volume occupied by P_0 vs P_1).
    </p>
    {chart_volumes}
    <h2>Substitution graphs</h2>
    <p>
        Edge from P_i to P_j labelled <code>k</code> means
        σ(P_j) contains <code>k</code> copies of P_i. Self-loops
        (⟲<code>k</code>) mean σ(P_i) contains <code>k</code> copies
        of P_i itself.
    </p>
    {chart_graphs}
    <h2>Next step</h2>
    <p>
        Each candidate proceeds to <code>realisation.py</code>
        (pending Q8 ruling) — a vertex-placement CSP that decides
        whether the matrix admits a 3D polyhedral realisation with
        ℤ[φ]³ vertices. Output: <code>Realised | NoRealisation |
        Inconclusive</code>. Survivors then run
        <code>aperiodicity_witness</code> for the four-pillar check.
    </p>
</body>
</html>
"""
    OUTPUT_HTML.write_text(html)
    print(f"\nWrote {OUTPUT_HTML} ({OUTPUT_HTML.stat().st_size:,} bytes)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
