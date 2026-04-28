"""Fourth-pillar argument for the Danzer ABCK aperiodic baseline.

**Status: stub.** The published proof of Danzer ABCK aperiodicity
exists; this module *does not* encode it. What it ships is a
concrete class satisfying ``apeiron.hierarchy.FourthPillarArgument``
whose two methods raise ``NotImplementedError``. Two reasons for
the stub-only form, per Claude (web)'s 2026-04-23 directive:

1. *Exercise the protocol's import path.* Without any concrete
   implementation in the wild, the runtime-checkable typing on
   ``FourthPillarArgument`` has never run end-to-end.
2. *Establish the directory convention* ``candidates/<name>/`` as
   the parent of every candidate's fourth-pillar implementation.
   Future deformation-search candidates derived from Danzer
   (``candidates/danzer_merge_AB/``, ``candidates/danzer_merge_ABC/``,
   ...) will inherit this directory structure and write their own
   concrete ``fourth_pillar.py`` files alongside.

**The published argument** (encoding deferred — research-grade work
distinct from the deformation-search baseline):

- Danzer 1989, *Three-dimensional analogs of the planar Penrose
  tilings and quasicrystals*, Discrete Mathematics 76(1–3), 1–7.
  Theorem 2 establishes ABCK aperiodicity.
- Goodman-Strauss 1998, *Matching rules and substitution tilings*,
  Annals of Mathematics 147(1), 181–223. Provides the matching-rule
  recognisability framework underlying pillars 2 and 3 of our
  four-pillar decomposition.
- Frettlöh, *Substitution tilings* (Bielefeld internal report),
  Theorem 1.2: every face-to-face tiling by A, B, C, K satisfying
  the blue-edge mirror matching rule is an ABCK tiling — i.e., a
  member of the substitution hierarchy. This is the pillar-4
  statement proper for Danzer (no tiling escapes the hierarchy).

Encoding any of these as executable verification logic is the
genuinely-Einstein step (CLAUDE.md §6.3 pillar 4) and is left to a
future commit if and when a Danzer-derived candidate needs the
formal verifier to *run* rather than just *exist as a stub*.
"""

from __future__ import annotations

from apeiron.corona import CoronaConfig
from apeiron.hierarchy import (
    FourthPillarArgument,
    HierarchicalCounterexample,
    HierarchicalWitness,
)


class DanzerABCKFourthPillar(FourthPillarArgument):
    """Concrete-but-stub ``FourthPillarArgument`` for the Danzer ABCK
    set. Methods raise ``NotImplementedError`` with a citation to
    the published proof.
    """

    _CITATION = (
        "Danzer 1989, Theorem 2; "
        "Goodman-Strauss 1998 (matching rules); "
        "Frettlöh, Bielefeld report, Theorem 1.2."
    )

    def local_configurations(self) -> frozenset[CoronaConfig]:
        """Every local configuration the Danzer ABCK tiling can
        exhibit — the input set for ``verify_hierarchical``.

        Stub: encoding the configuration enumeration requires the
        full geometric dissection (which the matrix-only baseline
        sub-commit 27B-α deliberately deferred to 27B-β) and the
        published proof's case analysis, which is out of scope.
        """
        raise NotImplementedError(
            "Danzer pillar-4 local-configuration enumeration is "
            f"not yet encoded; published proof: {self._CITATION}"
        )

    def verify_hierarchical(
        self, _config: CoronaConfig,
    ) -> HierarchicalWitness | HierarchicalCounterexample:
        """Decide hierarchical embedding for one ``CoronaConfig``.

        Stub: encoding the per-configuration verification requires
        Frettlöh Theorem 1.2's case analysis of the blue-edge mirror
        matching rule and is out of scope for the deformation-search
        baseline.
        """
        raise NotImplementedError(
            "Danzer pillar-4 hierarchical verification is not yet "
            f"encoded; published proof: {self._CITATION}"
        )
