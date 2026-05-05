# arXiv submission v1 — cover note and metadata

**Status.** Submission package for the v1 preprint of
[no_go_draft.md](no_go_draft.md). Per Q19b authorisation
(Claude (web), 2026-05-05): post to arXiv now, before the
Kramer-Papadopolos-Schlottmann-Zeidler (1994) inter-library loan
returns. Precedence is time-stamped at posting; the §(iv)
placeholder is honest, scoped, and standard for computational
mathematics preprints with a known pending update.

**v2 commitment.** §(iv) will be updated upon receipt of
Kramer-Papadopolos-Schlottmann-Zeidler (1994) via interlibrary
loan; v2 will either close or refine the scope of the result per
the two outcomes documented in §(iv) of the v1 draft.

---

## Title

Computational Evidence Against 3D Icosahedral Einstein Tiles
Among H₃-Fundamental-Region Substitution Candidates

## arXiv categories

- **Primary:** math.CO (Combinatorics)
- **Cross-list:** math.MG (Metric Geometry)

Per Q19b: the result is mathematical, not physics, and the
target audience (Goodman-Strauss, Frettlöh, the aperiodic
tilings community) reads math.CO. Physics readers interested in
quasicrystals will find it via the D₆ / icosahedral keywords
regardless. **Do not** post to math-ph as primary.

## MSC 2020 classification

- **Primary:** 52C23 — Quasicrystals and aperiodic tilings.
- **Secondary:** 52C20 — Tilings in 2 dimensions; 52B10 —
  Three-dimensional polytopes; 11H56 — Multiplicative theory of
  algebraic numbers (for the eigenvalue embedding).

## Cover note (for arXiv "Comments" field)

> v1, 2026-05-05. §(iv) addresses an open question pending
> inter-library loan of Kramer-Papadopolos-Schlottmann-Zeidler
> (1994) "Projection of the Danzer Tiling," J. Phys. A 27,
> 4505–4517. v2 will resolve §(iv) per the two outcomes
> documented therein. Source repository: github.com/Darian-Frey/Apeiron
> (private until v2 finalised).

## Authors

- Shane Hartley (corresponding) — independent researcher.
- (Acknowledgement of Claude (web) and Claude Code as research
  collaborators; their contributions are listed in the
  acknowledgements section of the draft, not as co-authors per
  standard arXiv conventions.)

## Submission checklist

- [ ] PDF generated from `no_go_draft.md` via pandoc + LaTeX
      (template: any of arXiv's standard math article templates;
      `arxiv-style.cls` if minimal customisation desired).
- [ ] All ASCII typography rendered correctly: math-script
      symbols (Σᵢ, ⋂, τ³), Unicode subscripts (E_⊥, t\*\_⊥), and
      the H₃ + ℤ[φ] notation.
- [ ] Figures: none in the v1 draft (text-only). If a figure
      becomes useful for v2 (e.g., a sector decomposition
      diagram from Hammock 2018), it should be cited from the
      original paper rather than reproduced.
- [ ] Citation style: tightened per Q19a-3 — Papadopolos 2000
      (journal year); arXiv preprint number supplementary, not
      primary.
- [ ] References cross-checked against
      [docs/literature_notes.md](literature_notes.md) and
      [docs/papers/MANIFEST.md](papers/MANIFEST.md).
- [ ] Test-cite for the eigenvalue embedding:
      [tests/integration/test_t2f_embedding.py](../tests/integration/test_t2f_embedding.py)
      (Apeiron repository). The 4 tests pass on commit `9779e7a`
      and following.
- [ ] Sweep-log cite:
      [notebooks/n3_phi3_max_entry_3_2026-05-04.log](../notebooks/n3_phi3_max_entry_3_2026-05-04.log).

## Why arXiv first, journal submission later

Per Q19b: the contribution is complete and citable today
independent of the 1994 ILL outcome. Holding the work off arXiv
for a months-uncertain ILL is the wrong trade-off. arXiv
posting is time-stamped; a clearly-scoped "this is what
remains" §(iv) is a feature, not a bug. *Experimental
Mathematics* in particular publishes results with explicit open
questions regularly; the journal submission can follow the v1
preprint or wait for v2 depending on whether the §(iv)
resolution timeline aligns with the submission cycle.

## Next steps for the user

1. Generate PDF from `no_go_draft.md` (the user's call on
   tooling — pandoc + LaTeX is the simplest path).
2. Register/log in to arXiv with the user's institutional or
   independent account.
3. Submit under math.CO with math.MG cross-list, MSC 52C23
   primary.
4. Paste the cover note above into the "Comments" field.
5. Allow arXiv's automated review window (typically 1–3 days).
6. Receive arXiv ID; log it in `docs/no_go_draft.md`'s "Notes
   on this draft" section as the v1 reference.

The user retains full control over the submission timing and
account ownership. Apeiron's role is to keep the draft
publication-ready and the supporting materials (notes, code,
tests, sweep logs) cite-stable.
