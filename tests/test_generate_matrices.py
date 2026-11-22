"""Unit tests for classifyFunctions.generateMatrices and the lookup of VAF/CCF values
(audit findings 15 and 24: errors must not be swallowed).
MatrixGenerator is replaced by a recorder; no genome needed.
Run with: pytest tests/test_generate_matrices.py -v
"""

import sys

import pytest

from SigProfilerClusters import classifyFunctions
from SigProfilerClusters.classifyFunctions import (
    _lookup_by_key,
    _subclass_has_mutations,
    generateMatrices,
)

HEADER = "project\tsamples\tID\n"


def _subclass(tmp_path, name, lines):
    folder = tmp_path / "subclasses" / name
    folder.mkdir(parents=True)
    (folder / f"{name}.txt").write_text(HEADER + "".join(lines))


def _run(tmp_path, subclasses):
    generateMatrices("GRCh37", str(tmp_path / "log.txt"), str(tmp_path) + "/", "", subclasses)


class TestGenerateMatrices:
    def test_generates_the_matrices_of_subclasses_with_mutations(self, tmp_path, monkeypatch):
        calls = []
        monkeypatch.setattr(
            classifyFunctions.matGen,
            "SigProfilerMatrixGeneratorFunc",
            lambda *a, **k: calls.append(a[0]),
        )
        _subclass(tmp_path, "class1a", ["P\tS1\t.\n"])
        _subclass(tmp_path, "class1b", [])  # header only
        _run(tmp_path, ["class1a", "class1b"])
        assert calls == ["class1a"]

    def test_a_failure_is_raised(self, tmp_path, monkeypatch):
        def fail(*a, **k):
            raise RuntimeError("matrix generation failed")

        monkeypatch.setattr(classifyFunctions.matGen, "SigProfilerMatrixGeneratorFunc", fail)
        _subclass(tmp_path, "class1a", ["P\tS1\t.\n"])
        with pytest.raises(RuntimeError, match="matrix generation failed"):
            _run(tmp_path, ["class1a"])

    def test_standard_output_is_restored_after_success_and_failure(self, tmp_path, monkeypatch):
        original = sys.stdout
        monkeypatch.setattr(
            classifyFunctions.matGen, "SigProfilerMatrixGeneratorFunc", lambda *a, **k: None
        )
        _subclass(tmp_path, "class1a", ["P\tS1\t.\n"])
        _run(tmp_path, ["class1a"])
        assert sys.stdout is original

        def fail(*a, **k):
            raise RuntimeError("x")

        monkeypatch.setattr(classifyFunctions.matGen, "SigProfilerMatrixGeneratorFunc", fail)
        with pytest.raises(RuntimeError):
            _run(tmp_path, ["class1a"])
        assert sys.stdout is original

    def test_missing_subclass_folder_is_skipped(self, tmp_path, monkeypatch):
        calls = []
        monkeypatch.setattr(
            classifyFunctions.matGen,
            "SigProfilerMatrixGeneratorFunc",
            lambda *a, **k: calls.append(a[0]),
        )
        _run(tmp_path, ["class1a"])
        assert calls == []


class TestSubclassHasMutations:
    def test_header_only_is_empty(self, tmp_path):
        _subclass(tmp_path, "c", [])
        assert _subclass_has_mutations(str(tmp_path / "subclasses" / "c")) is False

    def test_a_line_besides_the_header_counts(self, tmp_path):
        _subclass(tmp_path, "c", ["P\tS1\t.\n"])
        assert _subclass_has_mutations(str(tmp_path / "subclasses" / "c")) is True

    def test_blank_lines_do_not_count(self, tmp_path):
        _subclass(tmp_path, "c", ["\n", "\n"])
        assert _subclass_has_mutations(str(tmp_path / "subclasses" / "c")) is False


class TestLookupByKey:
    LINE = ["P", "S1", ".", "G", "SNP", "1", "100", "100", "C", "T"]

    def test_key_without_prefix(self):
        assert _lookup_by_key({"S1": {"1:100:C:T": 0.4}}, "S1", self.LINE, "VAF") == 0.4

    def test_key_with_chr_prefix(self):
        assert _lookup_by_key({"S1": {"chr1:100:C:T": 0.5}}, "S1", self.LINE, "VAF") == 0.5

    def test_missing_mutation_names_sample_and_mutation(self):
        with pytest.raises(KeyError, match="1:100:C:T of sample S1"):
            _lookup_by_key({"S1": {}}, "S1", self.LINE, "VAF")

    def test_missing_sample_is_a_key_error(self):
        with pytest.raises(KeyError, match="sample S1"):
            _lookup_by_key({}, "S1", self.LINE, "CCF")
