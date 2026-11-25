"""Audit finding 27: the samples that the analysis skips are collected with the reason and
reported together (screen, log and a table), instead of one warning that only reaches the
.err file or no message at all.
Run with: pytest tests/test_skipped_samples.py -v
"""

import pytest

import test_calculate_sample_imds as csi
from SigProfilerClusters import hotspot


class TestCalculateSampleImdsReportsTheSkippedSamples:
    def test_nothing_is_skipped_when_every_sample_has_data(self, tmp_path):
        csi._build(tmp_path)
        assert csi._run(tmp_path)["skipped_samples"] == {}

    def test_sample_without_an_original_file_has_its_reason(self, tmp_path):
        csi._build(tmp_path, sample="S1")
        csi._build(tmp_path, sample="S2", orig=False)
        skipped = csi._run(tmp_path, folders=("S1", "S2"))["skipped_samples"]
        assert list(skipped) == ["S2"]
        assert "no nearby mutations" in skipped["S2"]

    def test_the_same_with_the_density_correction(self, tmp_path, monkeypatch):
        csi._build(tmp_path, sample="S1", orig=False)
        csi._regions(monkeypatch, [])
        skipped = csi._run(tmp_path, correction=True)["skipped_samples"]
        assert "no nearby mutations" in skipped["S1"]

    def test_sample_without_distances_has_its_reason(self, tmp_path):
        csi._build(tmp_path)
        (tmp_path / "orig" / "S1_96_intradistance.txt").write_text("0\tS1\t1\t100\tC\tT\t0\n")
        skipped = csi._run(tmp_path)["skipped_samples"]
        assert "not enough distances" in skipped["S1"]

    def test_sample_without_simulation_files_has_its_reason(self, tmp_path):
        csi._build(tmp_path)
        for path in (tmp_path / "sims" / "S1").iterdir():
            path.unlink()
        assert "not enough distances" in csi._run(tmp_path)["skipped_samples"]["S1"]

    def test_hidden_entries_are_not_reported(self, tmp_path):
        csi._build(tmp_path)
        assert csi._run(tmp_path, folders=(".DS_Store",))["skipped_samples"] == {}


class TestReport:
    def test_prints_the_count_and_each_sample_sorted(self, tmp_path, capsys):
        hotspot.report_skipped_samples(
            {"B": "reason b", "A": "reason a"}, str(tmp_path / "skipped.txt"), str(tmp_path / "run.out")
        )
        out = capsys.readouterr().out
        assert "Skipped samples (2)" in out
        assert out.index("A: reason a") < out.index("B: reason b")

    def test_appends_to_the_log_without_deleting_it(self, tmp_path):
        log = tmp_path / "run.out"
        log.write_text("existing\n")
        hotspot.report_skipped_samples({"A": "reason a"}, str(tmp_path / "s.txt"), str(log))
        text = log.read_text()
        assert text.startswith("existing\n")
        assert "A: reason a" in text

    def test_writes_a_table_with_sample_and_reason(self, tmp_path):
        table = tmp_path / "skipped.txt"
        hotspot.report_skipped_samples(
            {"B": "reason b", "A": "reason a"}, str(table), str(tmp_path / "run.out")
        )
        assert table.read_text().splitlines() == ["sample\treason", "A\treason a", "B\treason b"]

    def test_without_skipped_samples_nothing_is_printed_or_written(self, tmp_path, capsys):
        log = tmp_path / "run.out"
        log.write_text("existing\n")
        hotspot.report_skipped_samples({}, str(tmp_path / "skipped.txt"), str(log))
        assert capsys.readouterr().out == ""
        assert log.read_text() == "existing\n"
        assert not (tmp_path / "skipped.txt").exists()

    def test_a_table_of_an_earlier_run_is_removed(self, tmp_path):
        table = tmp_path / "skipped.txt"
        table.write_text("old")
        hotspot.report_skipped_samples({}, str(table), str(tmp_path / "run.out"))
        assert not table.exists()
