# Literature

Annotated bibliography. Minimum reading for anyone contributing to Apeiron.
Expand each entry below with a paragraph of notes as the reference is
consulted — the point is to make the relevance to the 3D Einstein problem
explicit, not to restate the abstract.

## Primary references

- **Smith, D., Myers, J. S., Kaplan, C. S., & Goodman-Strauss, C. (2023).**
  *An aperiodic monotile.*
  — The hat paper. 2D. Establishes the computer-assisted local-configuration
  case analysis pattern that any 3D result will have to scale up.

- **Smith, D., Myers, J. S., Kaplan, C. S., & Goodman-Strauss, C. (2023).**
  *A chiral aperiodic monotile.*
  — The spectre paper. Strict Einstein in 2D via chirality; direct
  inspiration for the Strict Einstein stretch goal in 3D.

- **Goodman-Strauss, C. (1998).** *Matching rules and substitution tilings.*
  Annals of Mathematics, 147(1), 181–223.
  — Theoretical backbone for the four-pillar proof structure (CLAUDE.md
  §6.3).

- **Danzer, L. (1989).** *Three-dimensional analogs of the planar Penrose
  tilings and quasicrystals.* Discrete Mathematics, 76, 1–7.
  — The canonical 4-tile 3D aperiodic set. Starting point for Track A.

- **Socolar, J. E. S., & Taylor, J. M. (2011).** *An aperiodic hexagonal
  tile.* Journal of Combinatorial Theory A, 118, 2207–2231.
  — Cautionary example: non-local matching rules. Informs which proof
  obligations must be local.

- **Schmitt, P. (1988).** *An aperiodic prototile in space.* Unpublished
  manuscript.
  — The SCD biprism. Weak aperiodic only (via screw symmetry); the reason
  the Weak definition is insufficient for Apeiron's target.

## Working references for Track A

- **Frettlöh, D.** *Substitution tilings* (Bielefeld internal
  report). Online: tilings.math.uni-bielefeld.de.
  — Primary source for the Danzer ABCK encoding. Table 1 gives
  exact ℤ[φ]³ vertex coordinates for all four prototiles (with
  K's vertex 4 noted as half-integer relative to MF, requiring
  `scale_denom=2` storage). The substitution matrix is given
  explicitly; the dissection geometry is depicted in Figure 2 and
  the Tilings Encyclopedia interactive view, not listed
  algebraically. Frettlöh's Theorem 1.2 is the pillar-4 statement
  for Danzer (every face-to-face tiling by A,B,C,K with the
  blue-edge mirror matching rule belongs to the substitution
  hierarchy).

- **Senechal, M. (1995).** *Quasicrystals and Geometry.* Cambridge
  University Press.
  — Chapter 7: careful published transcription of the Danzer
  ABCK tiles into modern notation; useful as a secondary
  cross-check on Frettlöh Table 1.

- **Socolar, J. E. S., & Steinhardt, P. J. (1986).** *Quasicrystals.
  II. Unit-cell configurations.* Physical Review B, 34(2), 617–647.
  — Icosahedral tiles in Cartesian ℤ[φ]³ with explicit vertex lists.
  Cross-check for inflation factor (the Danzer ABCK substitution
  has linear inflation φ, volume scaling φ³ — these are not the
  same quantity).

- **Baake, M., & Grimm, U. (2013).** *Aperiodic Order, Volume 1: A
  Mathematical Invitation.* Cambridge University Press.
  — Chapter 6 covers the Ammann–Beenker tiling and its
  recognisability radius; used as the pillar-2 oracle fixture.
  The book's recognisability framework matches the
  ``hierarchy.is_recognisable`` predicate's contract.
