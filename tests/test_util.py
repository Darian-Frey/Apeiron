"""Tests for apeiron.util — candidate-file loading.

Scope: the ``load_candidate`` JSON loader. Each test constructs a
minimal JSON document via ``tmp_path``, invokes the loader, and
asserts either the resulting ``Polyhedron`` shape or a well-formed
``ValueError`` on malformed input.

Schema decisions (CLAUDE.md §3.2, §3.3; Claude (web) answer
2026-04-20) are the source of the acceptance criteria here.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

from apeiron.polyhedron import Polyhedron
from apeiron.symmetry import Vec3
from apeiron.util import load_candidate, pillar
from apeiron.zphi import ZPhi


def _write_json(tmp_path: Path, payload: Any, name: str = "candidate.json") -> Path:
    """Write ``payload`` to ``tmp_path / name`` and return the full path."""
    p = tmp_path / name
    p.write_text(json.dumps(payload), encoding="utf-8")
    return p


def _cube_payload(scale_denom: int = 1) -> dict[str, Any]:
    """Unit-cube vertex list under the loader's author-facing schema.

    Real coordinates: each axis component is 0 or 1. Under
    ``scale_denom=1`` the stored ×2 representation is 0 or 2. Under
    ``scale_denom=2`` the author writes the already-doubled form
    (0 or 2) so that stored = 2 * (0 or 2) / 2 = 0 or 2 — same result.
    """
    if scale_denom == 1:
        coords = (0, 1)
    elif scale_denom == 2:
        coords = (0, 2)
    else:
        raise ValueError(f"unsupported scale_denom in test helper: {scale_denom}")
    return {
        "name": "unit_cube",
        "scale_denom": scale_denom,
        "vertices": [
            [[i, 0], [j, 0], [k, 0]]
            for i in coords for j in coords for k in coords
        ],
    }


# -- happy-path loads ---------------------------------------------------


class TestLoadCandidateHappyPath:
    def test_cube_default_scale_denom(self, tmp_path: Path) -> None:
        payload = _cube_payload(scale_denom=1)
        p = _write_json(tmp_path, payload)
        P = load_candidate(p)
        assert isinstance(P, Polyhedron)
        assert len(P.vertices) == 8
        assert len(P.faces) == 6
        for face in P.faces:
            assert len(face) == 4   # quads

    def test_cube_omitted_scale_denom_defaults_to_one(self, tmp_path: Path) -> None:
        payload = _cube_payload(scale_denom=1)
        del payload["scale_denom"]   # exercise default
        p = _write_json(tmp_path, payload)
        P = load_candidate(p)
        assert len(P.vertices) == 8
        assert len(P.faces) == 6

    def test_cube_scale_denom_two(self, tmp_path: Path) -> None:
        # Same cube but author-coords already doubled; scale_denom=2
        # divides the ×2 multiplier back out, producing identical
        # stored coordinates.
        payload = _cube_payload(scale_denom=2)
        p = _write_json(tmp_path, payload)
        P = load_candidate(p)
        assert len(P.vertices) == 8
        assert len(P.faces) == 6

    def test_scale_denom_one_and_two_agree_on_same_tile(self, tmp_path: Path) -> None:
        A = load_candidate(_write_json(tmp_path, _cube_payload(1), "a.json"))
        B = load_candidate(_write_json(tmp_path, _cube_payload(2), "b.json"))
        assert A == B

    def test_accepts_str_and_path(self, tmp_path: Path) -> None:
        p = _write_json(tmp_path, _cube_payload(1))
        P_str = load_candidate(str(p))
        P_path = load_candidate(p)
        assert P_str == P_path

    def test_ignores_extra_top_level_keys(self, tmp_path: Path) -> None:
        payload = _cube_payload(1)
        payload["provenance"] = "made for testing"
        payload["symmetry"] = "cubic"
        payload["aperiodicity_claim"] = None
        p = _write_json(tmp_path, payload)
        # Extra keys must not cause failure; loader should ignore them.
        P = load_candidate(p)
        assert len(P.vertices) == 8

    def test_stored_coordinates_are_correct_zphi(self, tmp_path: Path) -> None:
        # Non-degenerate tetrahedron with one phi-valued coordinate.
        # Real (0, 0, phi) → stored (0, 0, 2*phi) = Vec3(_, _, ZPhi(0, 2)).
        payload = {
            "name": "phi_tet",
            "scale_denom": 1,
            "vertices": [
                [[0, 0], [0, 0], [0, 0]],
                [[1, 0], [0, 0], [0, 0]],
                [[0, 0], [1, 0], [0, 0]],
                [[0, 0], [0, 0], [0, 1]],   # real (0, 0, phi)
            ],
        }
        p = _write_json(tmp_path, payload)
        P = load_candidate(p)
        phi_vertex = Vec3(ZPhi(0, 0), ZPhi(0, 0), ZPhi(0, 2))
        assert phi_vertex in P.vertices


# -- scale_denom validation --------------------------------------------


class TestScaleDenomValidation:
    @pytest.mark.parametrize("bad", [0, -1, 3, 4, 100])
    def test_rejects_invalid_integer_scale_denom(
        self,
        tmp_path: Path,
        bad: int,
    ) -> None:
        payload = _cube_payload(1)
        payload["scale_denom"] = bad
        p = _write_json(tmp_path, payload)
        with pytest.raises(ValueError, match="scale_denom"):
            load_candidate(p)

    @pytest.mark.parametrize("bad", [1.0, "1", None, [1]])
    def test_rejects_non_int_scale_denom(
        self,
        tmp_path: Path,
        bad: object,
    ) -> None:
        payload = _cube_payload(1)
        payload["scale_denom"] = bad   # type: ignore[assignment]
        p = _write_json(tmp_path, payload)
        with pytest.raises(ValueError, match="scale_denom must be an int"):
            load_candidate(p)

    def test_rejects_bool_scale_denom(self, tmp_path: Path) -> None:
        # bool is a subclass of int in Python; we reject it explicitly
        # to match ZPhi's type discipline.
        payload = _cube_payload(1)
        payload["scale_denom"] = True   # type: ignore[assignment]
        p = _write_json(tmp_path, payload)
        with pytest.raises(ValueError, match="scale_denom must be an int"):
            load_candidate(p)


# -- schema rejections -------------------------------------------------


class TestSchemaRejections:
    def test_file_not_found(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError):
            load_candidate(tmp_path / "does_not_exist.json")

    def test_top_level_must_be_object(self, tmp_path: Path) -> None:
        p = _write_json(tmp_path, ["not", "an", "object"])
        with pytest.raises(ValueError, match="top-level JSON must be an object"):
            load_candidate(p)

    def test_missing_name(self, tmp_path: Path) -> None:
        payload = _cube_payload(1)
        del payload["name"]
        p = _write_json(tmp_path, payload)
        with pytest.raises(ValueError, match="'name'"):
            load_candidate(p)

    def test_empty_name_rejected(self, tmp_path: Path) -> None:
        payload = _cube_payload(1)
        payload["name"] = ""
        p = _write_json(tmp_path, payload)
        with pytest.raises(ValueError, match="'name'"):
            load_candidate(p)

    def test_non_string_name_rejected(self, tmp_path: Path) -> None:
        payload = _cube_payload(1)
        payload["name"] = 42   # type: ignore[assignment]
        p = _write_json(tmp_path, payload)
        with pytest.raises(ValueError, match="'name'"):
            load_candidate(p)

    def test_missing_vertices(self, tmp_path: Path) -> None:
        payload = _cube_payload(1)
        del payload["vertices"]
        p = _write_json(tmp_path, payload)
        with pytest.raises(ValueError, match="'vertices'"):
            load_candidate(p)

    def test_vertices_not_a_list(self, tmp_path: Path) -> None:
        payload = _cube_payload(1)
        payload["vertices"] = "not a list"   # type: ignore[assignment]
        p = _write_json(tmp_path, payload)
        with pytest.raises(ValueError, match="'vertices' must be a list"):
            load_candidate(p)

    def test_vertex_wrong_length(self, tmp_path: Path) -> None:
        payload = _cube_payload(1)
        payload["vertices"][0] = [[0, 0], [0, 0]]   # only 2 components
        p = _write_json(tmp_path, payload)
        with pytest.raises(ValueError, match="vertex\\[0\\]"):
            load_candidate(p)

    def test_vertex_component_wrong_length(self, tmp_path: Path) -> None:
        payload = _cube_payload(1)
        payload["vertices"][0] = [[0, 0, 0], [0, 0], [0, 0]]   # 3-elem
        p = _write_json(tmp_path, payload)
        with pytest.raises(ValueError, match="component 0"):
            load_candidate(p)

    def test_vertex_component_with_float(self, tmp_path: Path) -> None:
        payload = _cube_payload(1)
        payload["vertices"][0] = [[0.0, 0], [0, 0], [0, 0]]   # float
        p = _write_json(tmp_path, payload)
        with pytest.raises(ValueError, match="component 0"):
            load_candidate(p)

    def test_vertex_component_with_bool(self, tmp_path: Path) -> None:
        payload = _cube_payload(1)
        payload["vertices"][0] = [[True, 0], [0, 0], [0, 0]]   # bool
        p = _write_json(tmp_path, payload)
        with pytest.raises(ValueError, match="component 0"):
            load_candidate(p)

    def test_fewer_than_four_vertices(self, tmp_path: Path) -> None:
        payload = {
            "name": "degenerate",
            "scale_denom": 1,
            "vertices": [
                [[0, 0], [0, 0], [0, 0]],
                [[1, 0], [0, 0], [0, 0]],
                [[0, 0], [1, 0], [0, 0]],
            ],
        }
        p = _write_json(tmp_path, payload)
        with pytest.raises(ValueError, match="at least 4"):
            load_candidate(p)


# -- allow_interior_inputs flag ----------------------------------------


class TestAllowInteriorInputs:
    def _cube_plus_interior_payload(self) -> dict[str, Any]:
        """Unit cube plus the centroid (0.5, 0.5, 0.5). Under
        ``scale_denom=1`` that last point is genuinely interior to the
        hull of the cube corners (and lies on no face since the hull
        centroid is strictly inside).

        We bump the outer coordinates to 2 and add an interior point
        at (1, 1, 1) — real coord 1, stored 2, strictly inside the
        cube from 0..2.
        """
        corners = [
            [[i, 0], [j, 0], [k, 0]]
            for i in (0, 2) for j in (0, 2) for k in (0, 2)
        ]
        interior = [[[1, 0], [1, 0], [1, 0]]]
        return {
            "name": "cube_plus_interior",
            "scale_denom": 1,
            "vertices": corners + interior,
        }

    def test_default_rejects_interior_input(self, tmp_path: Path) -> None:
        payload = self._cube_plus_interior_payload()
        p = _write_json(tmp_path, payload)
        with pytest.raises(ValueError, match="interior to the convex hull"):
            load_candidate(p)

    def test_allow_flag_accepts_interior_input(self, tmp_path: Path) -> None:
        payload = self._cube_plus_interior_payload()
        p = _write_json(tmp_path, payload)
        P = load_candidate(p, allow_interior_inputs=True)
        # Interior point is silently filtered out; hull is the cube.
        assert len(P.vertices) == 8
        assert len(P.faces) == 6

    def test_error_message_includes_interior_index(self, tmp_path: Path) -> None:
        payload = self._cube_plus_interior_payload()
        p = _write_json(tmp_path, payload)
        with pytest.raises(ValueError, match=r"indices: \[8\]"):
            load_candidate(p)

    def test_error_message_counts_multiple_interiors(self, tmp_path: Path) -> None:
        payload = self._cube_plus_interior_payload()
        # Append a second interior point — duplicate of the first is
        # fine for this test, since both input-list occurrences appear
        # in interior_idx when they're not extreme.
        payload["vertices"].append([[1, 0], [1, 0], [1, 0]])
        p = _write_json(tmp_path, payload)
        with pytest.raises(ValueError, match=r"2 of 10 input vertices"):
            load_candidate(p)


# -- pillar decorator --------------------------------------------------


class TestPillarDecorator:
    def test_attaches_pillar_attribute_to_function(self) -> None:
        @pillar(2)
        def f() -> int:
            return 42

        assert f._pillar == 2
        assert f() == 42  # behaviour preserved

    def test_attaches_pillar_attribute_to_class(self) -> None:
        @pillar(4)
        class C:
            pass

        assert C._pillar == 4

    def test_pillar_zero_is_rejected(self) -> None:
        with pytest.raises(ValueError, match=r"\{1, 2, 3, 4\}"):
            pillar(0)

    def test_pillar_five_is_rejected(self) -> None:
        with pytest.raises(ValueError, match=r"\{1, 2, 3, 4\}"):
            pillar(5)

    def test_decorator_returns_target_unchanged(self) -> None:
        # Identity except for the attribute attachment.
        def f() -> str:
            return "hello"

        decorated = pillar(1)(f)
        assert decorated is f
        assert decorated() == "hello"

    def test_no_interior_no_error(self, tmp_path: Path) -> None:
        # Plain cube, no interior points.
        payload = _cube_payload(1)
        p = _write_json(tmp_path, payload)
        # Default strict mode, no error.
        P = load_candidate(p)
        assert len(P.vertices) == 8

    def test_allow_flag_keyword_only(self, tmp_path: Path) -> None:
        # Positional-argument form must not accept the flag, preventing
        # accidental mis-use where a caller passes a Path-like str and
        # then a truthy second positional.
        payload = _cube_payload(1)
        p = _write_json(tmp_path, payload)
        with pytest.raises(TypeError):
            load_candidate(p, True)   # type: ignore[misc]
