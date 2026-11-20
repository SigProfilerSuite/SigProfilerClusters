"""Unit tests for the path and line helpers in SigProfilerClusters.convertToVCF.
No genome reference required.
Run with: pytest tests/test_convert_to_vcf.py -v
"""

import os

import pytest

from SigProfilerClusters.convertToVCF import (
    checkPaths,
    generateAllPaths,
    pullLineInfo,
)


# ─── checkPaths ──────────────────────────────────────────────


class TestCheckPaths:
    def test_creates_missing_directory(self, tmp_path):
        checkPaths(str(tmp_path), "SBS")
        assert (tmp_path / "SBS").is_dir()

    def test_recreates_existing_directory_empty(self, tmp_path):
        stale = tmp_path / "SBS"
        stale.mkdir()
        (stale / "old.txt").write_text("stale")
        checkPaths(str(tmp_path), "SBS")
        assert (tmp_path / "SBS").is_dir()
        assert os.listdir(str(tmp_path / "SBS")) == []

    def test_does_not_touch_sibling_directories(self, tmp_path):
        (tmp_path / "other").mkdir()
        (tmp_path / "other" / "keep.txt").write_text("keep")
        checkPaths(str(tmp_path), "SBS")
        assert (tmp_path / "other" / "keep.txt").exists()


# ─── generateAllPaths ────────────────────────────────────────


class TestGenerateAllPaths:
    def test_snv_layout(self, tmp_path):
        generateAllPaths(str(tmp_path), "96")
        clustered = tmp_path / "output" / "clustered"
        for name in ("DBS", "MBS", "omikli", "kataegis", "other"):
            assert (clustered / name).is_dir()
        assert (tmp_path / "output" / "nonClustered" / "SBS").is_dir()

    def test_indel_layout(self, tmp_path):
        generateAllPaths(str(tmp_path), "ID")
        assert (tmp_path / "output" / "clustered" / "ID").is_dir()
        assert (tmp_path / "output" / "nonClustered" / "ID").is_dir()

    def test_indel_layout_has_no_snv_folders(self, tmp_path):
        generateAllPaths(str(tmp_path), "ID")
        assert not (tmp_path / "output" / "clustered" / "kataegis").exists()

    def test_second_call_clears_previous_output(self, tmp_path):
        generateAllPaths(str(tmp_path), "96")
        leftover = tmp_path / "output" / "clustered" / "kataegis" / "old.vcf"
        leftover.write_text("stale")
        generateAllPaths(str(tmp_path), "96")
        assert not leftover.exists()


# ─── pullLineInfo ────────────────────────────────────────────

# Columns are whitespace separated: index 1 sample, 5 chrom, 6 position,
# 8 ref, 9 alt, 12 and 14 the two VAF-related values.
LINE = "x\tsample1\tx\tx\tx\t7\t12345\tx\tG\tT\tx\tx\t0.4\tx\t0.6\n"


class TestPullLineInfo:
    def test_with_vaf(self):
        assert pullLineInfo(LINE) == ("sample1", "7", "12345", "G", "T", "0.6", "0.4")

    def test_without_vaf(self):
        assert pullLineInfo(LINE, collectVAF=False) == (
            "sample1",
            "7",
            "12345",
            "G",
            "T",
        )

    def test_trailing_newline_is_ignored(self):
        assert pullLineInfo(LINE.rstrip("\n"), collectVAF=False) == pullLineInfo(
            LINE, collectVAF=False
        )

    def test_space_separated_columns(self):
        space_line = " ".join(LINE.strip().split("\t")) + "\n"
        assert pullLineInfo(space_line, collectVAF=False)[0] == "sample1"

    def test_short_line_without_vaf_columns_raises(self):
        short = "x\ts\tx\tx\tx\t1\t100\tx\tA\tC\n"
        with pytest.raises(IndexError):
            pullLineInfo(short, collectVAF=True)
        assert pullLineInfo(short, collectVAF=False) == ("s", "1", "100", "A", "C")

    def test_too_few_columns_raises(self):
        with pytest.raises(IndexError):
            pullLineInfo("only\ttwo\n", collectVAF=False)
