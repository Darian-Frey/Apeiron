"""Supertile construction and recognisability (border forcing).

Organised around the four pillars of a Strong Einstein proof (CLAUDE.md
§6.3):

1. Substitution exists and is primitive.
2. Recognisability: every tile's supertile membership is determined by a
   bounded local neighbourhood.
3. Aperiodicity from recognisability (standard inflation argument).
4. No non-hierarchical tiling exists — the genuinely Einstein step.

Every public function must be tagged in its docstring with which pillar it
establishes.

Blocked on ``substitution`` and ``corona`` (CLAUDE.md §5.2).
"""
