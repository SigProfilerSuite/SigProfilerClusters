"""Unit tests for convertToVCF.convertFiles (clustered/non-clustered text -> per-sample VCF).
Inputs are small dummy tables built with tmp_path. No genome reference required.
Run with: pytest tests/test_convert_files.py -v
"""

import os

import pytest

from SigProfilerClusters.convertToVCF import convertFiles, generateAllPaths

PROJECT = "proj"
SNV_CLASSES = {
    "class1a": "DBS",
    "class1b": "MBS",
    "class1c": "omikli",
    "class2": "kataegis",
    "class3": "other",
}
HEADER = "header\n"


def _row(sample, chrom, pos, ref, alt, vaf="0.5", group="1"):
    """Build a whitespace-separated row with the columns pullLineInfo reads:
    1 sample, 5 chrom, 6 pos, 8 ref, 9 alt, 12 group, 14 vaf."""
    cols = ["x"] * 15
    cols[1] = sample
    cols[5] = chrom
    cols[6] = pos
    cols[8] = ref
    cols[9] = alt
    cols[12] = group
    cols[14] = vaf
    return "\t".join(cols) + "\n"


def _write(path, rows):
    os.makedirs(os.path.dirname(str(path)), exist_ok=True)
    with open(str(path), "w") as f:
        f.write(HEADER)
        f.writelines(rows)


def _read(path):
    with open(str(path)) as f:
        return [line.rstrip("\n").split("\t") for line in f]


def _base(tmp_path, corrected=True):
    return tmp_path / "output" / ("vcf_files_corrected" if corrected else "vcf_files")


def _setup_id(tmp_path, clustered_rows, non_clustered_rows, corrected=True):
    base = _base(tmp_path, corrected)
    _write(base / f"{PROJECT}_clustered" / "INDEL" / f"{PROJECT}_clustered.txt", clustered_rows)
    _write(
        base / f"{PROJECT}_nonClustered" / "INDEL" / f"{PROJECT}_nonClustered.txt",
        non_clustered_rows,
    )
    generateAllPaths(str(tmp_path), "ID")


def _setup_snv(tmp_path, rows_by_class, non_clustered_rows, corrected=True):
    base = _base(tmp_path, corrected)
    for cls in SNV_CLASSES:
        _write(
            base / f"{PROJECT}_clustered" / "subclasses" / cls / f"{PROJECT}_clustered_{cls}.txt",
            rows_by_class.get(cls, []),
        )
    _write(base / f"{PROJECT}_nonClustered" / "SNV" / f"{PROJECT}_nonClustered.txt", non_clustered_rows)
    generateAllPaths(str(tmp_path), "96")


# ─── ID (indel) context ──────────────────────────────────────


class TestConvertFilesIndel:
    def test_clustered_and_non_clustered_vcfs_are_written(self, tmp_path):
        _setup_id(
            tmp_path,
            [_row("S1", "1", "100", "A", "AT")],
            [_row("S1", "2", "200", "GC", "G")],
        )
        convertFiles(str(tmp_path), "ID", PROJECT, True)

        clustered = _read(tmp_path / "output" / "clustered" / "ID" / "S1.vcf")
        assert clustered[0] == ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
        assert clustered[1] == ["1", "100", "S1", "A", "AT", ".", ".", "ClusteredAnalysis"]

        non_clustered = _read(tmp_path / "output" / "nonClustered" / "ID" / "S1.vcf")
        assert non_clustered[1] == ["2", "200", "S1", "GC", "G", ".", ".", "NonClustered"]

    def test_one_file_per_sample(self, tmp_path):
        rows = [
            _row("S1", "1", "100", "A", "AT"),
            _row("S1", "1", "150", "C", "CA"),
            _row("S2", "3", "300", "T", "TG"),
        ]
        _setup_id(tmp_path, rows, [_row("S1", "1", "1", "A", "AC")])
        convertFiles(str(tmp_path), "ID", PROJECT, True)

        out_dir = tmp_path / "output" / "clustered" / "ID"
        assert sorted(os.listdir(str(out_dir))) == ["S1.vcf", "S2.vcf"]
        assert len(_read(out_dir / "S1.vcf")) == 3  # header + two records
        assert len(_read(out_dir / "S2.vcf")) == 2  # header + one record

    def test_uncorrected_inputs_are_read_from_vcf_files(self, tmp_path):
        _setup_id(
            tmp_path,
            [_row("S1", "5", "55", "A", "AT")],
            [_row("S1", "6", "66", "G", "GA")],
            corrected=False,
        )
        convertFiles(str(tmp_path), "ID", PROJECT, False)
        assert _read(tmp_path / "output" / "clustered" / "ID" / "S1.vcf")[1][0] == "5"

    def test_missing_input_file_raises(self, tmp_path):
        generateAllPaths(str(tmp_path), "ID")
        with pytest.raises(FileNotFoundError):
            convertFiles(str(tmp_path), "ID", PROJECT, True)

    @pytest.mark.xfail(
        reason="empty clustered ID table: 'if len(allLines) == 0: pass' does not "
        "stop execution, so allLines[0] raises IndexError",
        raises=IndexError,
        strict=True,
    )
    def test_empty_clustered_table(self, tmp_path):
        _setup_id(tmp_path, [], [_row("S1", "1", "1", "A", "AT")])
        convertFiles(str(tmp_path), "ID", PROJECT, True)


# ─── SBS (SNV) context ───────────────────────────────────────


class TestConvertFilesSnv:
    def test_subclass_outputs_go_to_their_folders(self, tmp_path):
        rows = {
            "class1a": [_row("S1", "1", "10", "A", "C", vaf="0.3", group="7")],
            "class2": [_row("S1", "2", "20", "G", "T", vaf="0.6", group="8")],
        }
        _setup_snv(tmp_path, rows, [_row("S1", "9", "90", "C", "A")])
        convertFiles(str(tmp_path), "96", PROJECT, True)

        dbs = _read(tmp_path / "output" / "clustered" / "DBS" / "S1.vcf")
        assert dbs[1][:5] == ["1", "10", "S1", "A", "C"]
        kataegis = _read(tmp_path / "output" / "clustered" / "kataegis" / "S1.vcf")
        assert kataegis[1][:5] == ["2", "20", "S1", "G", "T"]
        assert not (tmp_path / "output" / "clustered" / "omikli" / "S1.vcf").exists()

    def test_info_field_carries_vaf_and_group(self, tmp_path):
        rows = {"class3": [_row("S1", "1", "10", "A", "C", vaf="0.25", group="4")]}
        _setup_snv(tmp_path, rows, [_row("S1", "9", "90", "C", "A")])
        convertFiles(str(tmp_path), "96", PROJECT, True)

        info = _read(tmp_path / "output" / "clustered" / "other" / "S1.vcf")[1][7]
        assert info.startswith("ClusteredAnalysis;VAF=0.25;groupNumber")
        assert info.endswith("4")

    def test_non_clustered_goes_to_sbs_folder(self, tmp_path):
        _setup_snv(tmp_path, {}, [_row("S1", "9", "90", "C", "A"), _row("S2", "8", "80", "T", "G")])
        convertFiles(str(tmp_path), "96", PROJECT, True)

        out_dir = tmp_path / "output" / "nonClustered" / "SBS"
        assert sorted(os.listdir(str(out_dir))) == ["S1.vcf", "S2.vcf"]
        assert _read(out_dir / "S2.vcf")[1] == ["8", "80", "S2", "T", "G", ".", ".", "NonClustered"]

    def test_empty_subclass_tables_are_skipped(self, tmp_path):
        _setup_snv(tmp_path, {}, [_row("S1", "9", "90", "C", "A")])
        convertFiles(str(tmp_path), "96", PROJECT, True)
        for folder in SNV_CLASSES.values():
            assert os.listdir(str(tmp_path / "output" / "clustered" / folder)) == []

    def test_multiple_samples_in_one_subclass(self, tmp_path):
        rows = {
            "class1c": [
                _row("S1", "1", "10", "A", "C"),
                _row("S2", "2", "20", "G", "T"),
                _row("S2", "2", "21", "G", "A"),
            ]
        }
        _setup_snv(tmp_path, rows, [_row("S1", "9", "90", "C", "A")])
        convertFiles(str(tmp_path), "96", PROJECT, True)

        out_dir = tmp_path / "output" / "clustered" / "omikli"
        assert sorted(os.listdir(str(out_dir))) == ["S1.vcf", "S2.vcf"]
        assert len(_read(out_dir / "S2.vcf")) == 3
