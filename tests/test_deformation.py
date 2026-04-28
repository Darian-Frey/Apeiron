"""Tests for apeiron.deformation — Track A Phase 1.5 face-merge
deformation-search infrastructure."""

from __future__ import annotations

import shutil
import tempfile
from pathlib import Path

import pytest

from apeiron.deformation import scaffold_merge_candidate


@pytest.fixture
def tmp_candidates_dir() -> Path:
    """Temp directory for scaffold tests; cleaned up after each test."""
    tmp = Path(tempfile.mkdtemp())
    yield tmp
    shutil.rmtree(tmp)


# -- scaffold_merge_candidate (Q6c) -----------------------------------


class TestScaffoldMergeCandidate:
    """Per Claude (web) Q6c ruling 2026-04-29: every face-merge
    candidate gets an auto-generated ``fourth_pillar.py`` stub at
    scaffold time, regardless of whether it later survives pillars
    1+2+3.
    """

    def test_creates_directory(self, tmp_candidates_dir: Path) -> None:
        path = scaffold_merge_candidate(
            name="merge_ab_test",
            parent_letters=("A", "B"),
            shared_face_description="face 2 of A coincides with face 0 of B",
            candidates_dir=tmp_candidates_dir,
        )
        assert path.exists()
        assert path.is_dir()
        assert path == tmp_candidates_dir / "merge_ab_test"

    def test_scaffolds_fourth_pillar_and_init(
        self, tmp_candidates_dir: Path,
    ) -> None:
        path = scaffold_merge_candidate(
            name="merge_ab_test",
            parent_letters=("A", "B"),
            shared_face_description="shared face details",
            candidates_dir=tmp_candidates_dir,
        )
        assert (path / "fourth_pillar.py").exists()
        assert (path / "__init__.py").exists()

    def test_fourth_pillar_module_imports(
        self, tmp_candidates_dir: Path,
    ) -> None:
        path = scaffold_merge_candidate(
            name="merge_ab_test",
            parent_letters=("A", "B"),
            shared_face_description="shared face",
            candidates_dir=tmp_candidates_dir,
        )
        # Dynamically import the generated module.
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            "test_scaffold_module", path / "fourth_pillar.py",
        )
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        # The PascalCase'd class name must exist.
        cls = mod.MergeAbTestFourthPillar
        assert cls is not None

    def test_generated_class_satisfies_protocol(
        self, tmp_candidates_dir: Path,
    ) -> None:
        from apeiron.hierarchy import FourthPillarArgument

        path = scaffold_merge_candidate(
            name="merge_ab_test",
            parent_letters=("A", "B"),
            shared_face_description="shared face",
            candidates_dir=tmp_candidates_dir,
        )
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            "test_scaffold_module", path / "fourth_pillar.py",
        )
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        instance = mod.MergeAbTestFourthPillar()
        assert isinstance(instance, FourthPillarArgument)

    def test_methods_raise_not_implemented_with_context(
        self, tmp_candidates_dir: Path,
    ) -> None:
        from apeiron.corona import CoronaConfig
        path = scaffold_merge_candidate(
            name="merge_test_ctx",
            parent_letters=("A", "B"),
            shared_face_description="distinctive description string",
            candidates_dir=tmp_candidates_dir,
        )
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            "test_scaffold_ctx", path / "fourth_pillar.py",
        )
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        instance = mod.MergeTestCtxFourthPillar()
        # The error must include the context (parent letters + face desc).
        with pytest.raises(NotImplementedError, match="distinctive description"):
            instance.local_configurations()
        with pytest.raises(NotImplementedError, match="distinctive description"):
            instance.verify_hierarchical(_config=None)  # type: ignore[arg-type]

    def test_rejects_existing_directory(
        self, tmp_candidates_dir: Path,
    ) -> None:
        # First call creates; second call refuses to overwrite.
        scaffold_merge_candidate(
            name="merge_ab_test",
            parent_letters=("A", "B"),
            shared_face_description="x",
            candidates_dir=tmp_candidates_dir,
        )
        with pytest.raises(FileExistsError):
            scaffold_merge_candidate(
                name="merge_ab_test",
                parent_letters=("A", "B"),
                shared_face_description="x",
                candidates_dir=tmp_candidates_dir,
            )

    def test_rejects_invalid_name(
        self, tmp_candidates_dir: Path,
    ) -> None:
        for bad_name in ("", "../escape", ".hidden", "with/slash"):
            with pytest.raises(ValueError, match="directory-safe"):
                scaffold_merge_candidate(
                    name=bad_name,
                    parent_letters=("A",),
                    shared_face_description="x",
                    candidates_dir=tmp_candidates_dir,
                )

    def test_pillar_4_approach_appears_in_docstring(
        self, tmp_candidates_dir: Path,
    ) -> None:
        # Custom approach text is captured in the generated module.
        custom_approach = (
            "Inherit from Penrose P3's matching-rule restriction and "
            "verify the merged tile's edge labels are I_h-equivariant."
        )
        path = scaffold_merge_candidate(
            name="merge_with_custom_approach",
            parent_letters=("A", "K"),
            shared_face_description="shared face",
            pillar_4_approach=custom_approach,
            candidates_dir=tmp_candidates_dir,
        )
        text = (path / "fourth_pillar.py").read_text()
        assert custom_approach in text

    def test_parent_letters_recorded(
        self, tmp_candidates_dir: Path,
    ) -> None:
        path = scaffold_merge_candidate(
            name="merge_abk_test",
            parent_letters=("A", "B", "K"),
            shared_face_description="x",
            candidates_dir=tmp_candidates_dir,
        )
        text = (path / "fourth_pillar.py").read_text()
        assert "A+B+K" in text
