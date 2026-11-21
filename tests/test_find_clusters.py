"""Unit tests for classifyFunctions.findClustersOfClusters (subclassification of clustered
SBS that carry a VAF/CCF value). Same synthetic events as test_find_clusters_novaf.py,
each with a similar VAF; some events mix very different VAFs and are declassified.
Inputs are built with tmp_path; the matrix generation step needs a genome and is stubbed.
Run with: pytest tests/test_find_clusters.py -v
"""

import os
import pickle

import pytest

from clustered_table_helpers import (
    GENOME,
    HEADER_VAF,
    IMD_CUTOFF,
    PROJECT,
    SUBCLASSES,
    entries as _entries,
    no_matrices as _no_matrices,
    populated as _populated,
    positions as _positions,
    row as _row,
    setup,
    subclass_path as _subclass_path,
)
from SigProfilerClusters import classifyFunctions
from SigProfilerClusters.classifyFunctions import findClustersOfClusters

VAF_COL = 14  # VAF/CCF column of the subclass files
GROUP_COL = 12  # group number column of the subclass files


@pytest.fixture(autouse=True)
def stub_matrix_generation(monkeypatch):
    monkeypatch.setattr(classifyFunctions, "generateMatrices", _no_matrices)


def _vrow(pos, vaf, **kwargs):
    return _row(pos, vaf=vaf, **kwargs)


# findClustersOfClusters never processes the last event of the table (see the xfail test
# below). Every test table therefore ends with this event, which is never reported.
SENTINEL = [_vrow(900000, 0.4), _vrow(900001, 0.4)]


def _run(tmp_path, rows, correction=False, ccf=False, imds=None, sentinel=True):
    rows = list(rows) + (SENTINEL if sentinel else [])
    clustered = setup(
        tmp_path,
        rows,
        HEADER_VAF,
        f"{PROJECT}_clustered_vaf.txt",
        imds,
        correction,
    )
    findClustersOfClusters(
        PROJECT,
        False,
        str(tmp_path) + "/",
        1000000,
        {GENOME: {"1": 0, "2": 1000000}},
        {sample: [] for sample in (imds or {"S1": 0})},
        None,
        GENOME,
        1,
        {},
        correction,
        ccf,
    )
    return clustered


# ─── event classes with similar VAFs ─────────────────────────


class TestEventClassesWithSimilarVaf:
    def test_adjacent_pair_is_dbs(self, tmp_path):
        clustered = _run(tmp_path, [_vrow(1000, 0.40), _vrow(1001, 0.42)])
        assert _populated(clustered) == ["class1", "class1a"]
        assert _positions(clustered, "class1a") == [1000, 1001]

    def test_adjacent_run_of_three_is_mbs(self, tmp_path):
        rows = [_vrow(1000, 0.40), _vrow(1001, 0.41), _vrow(1002, 0.42)]
        clustered = _run(tmp_path, rows)
        assert _populated(clustered) == ["class1", "class1b"]
        assert _positions(clustered, "class1b") == [1000, 1001, 1002]

    def test_two_mutations_apart_are_omikli(self, tmp_path):
        clustered = _run(tmp_path, [_vrow(1000, 0.40), _vrow(1050, 0.42)])
        assert _populated(clustered) == ["class1", "class1c"]

    def test_three_mutations_within_cutoff_are_omikli(self, tmp_path):
        rows = [_vrow(1000, 0.40), _vrow(1050, 0.42), _vrow(1100, 0.41)]
        clustered = _run(tmp_path, rows)
        assert _populated(clustered) == ["class1", "class1c"]
        assert _positions(clustered, "class1c") == [1000, 1050, 1100]

    def test_four_mutations_within_cutoff_are_class2(self, tmp_path):
        rows = [_vrow(1000 + 50 * i, vaf) for i, vaf in enumerate([0.40, 0.42, 0.44, 0.43])]
        clustered = _run(tmp_path, rows)
        assert "class2" in _populated(clustered)
        assert "class1" not in _populated(clustered)
        assert _positions(clustered, "class2") == [1000, 1050, 1100, 1150]

    def test_distance_equal_to_imd_cutoff_is_clustered(self, tmp_path):
        rows = [_vrow(1000, 0.40), _vrow(1000 + IMD_CUTOFF, 0.42)]
        assert _populated(_run(tmp_path, rows)) == ["class1", "class1c"]

    def test_vaf_is_kept_in_the_output_rows(self, tmp_path):
        clustered = _run(tmp_path, [_vrow(1000, 0.40), _vrow(1001, 0.42)])
        assert [cols[VAF_COL] for cols in _entries(clustered, "class1a")] == ["0.4", "0.42"]

    def test_class_label_is_the_last_column(self, tmp_path):
        clustered = _run(tmp_path, [_vrow(1000, 0.40), _vrow(1001, 0.42)])
        assert {cols[-1] for cols in _entries(clustered, "class1")} == {"ClassI"}
        assert {cols[-1] for cols in _entries(clustered, "class1a")} == {"ClassIA"}


class TestClass2Processivity:
    def _class2(self, tmp_path, refs):
        rows = [_vrow(1000 + 50 * i, 0.40 + 0.01 * i, ref=ref) for i, ref in enumerate(refs)]
        return _run(tmp_path, rows)

    def test_pyrimidine_event_goes_to_2y(self, tmp_path):
        clustered = self._class2(tmp_path, "CCCC")
        assert _positions(clustered, "class2Y") == [1000, 1050, 1100, 1150]
        assert _entries(clustered, "class2K") == []

    def test_keto_event_goes_to_2k(self, tmp_path):
        assert len(_entries(self._class2(tmp_path, "TGTG"), "class2K")) == 4

    def test_strong_event_goes_to_2s(self, tmp_path):
        assert len(_entries(self._class2(tmp_path, "CGCG"), "class2S")) == 4

    def test_mixed_event_goes_to_2n(self, tmp_path):
        assert len(_entries(self._class2(tmp_path, "CAGT"), "class2N")) == 4


# ─── events declassified by VAF ──────────────────────────────


class TestVafDeclassification:
    def test_dbs_with_different_vafs_is_not_an_event(self, tmp_path):
        clustered = _run(tmp_path, [_vrow(1000, 0.2), _vrow(1001, 0.6)])
        assert _populated(clustered) == ["class3"]
        assert _positions(clustered, "class3") == [1000, 1001]

    def test_declassified_rows_carry_class_and_reason(self, tmp_path):
        clustered = _run(tmp_path, [_vrow(1000, 0.2), _vrow(1001, 0.6)])
        assert {tuple(cols[-2:]) for cols in _entries(clustered, "class3")} == {("ClassIII", "vaf")}
        header = _subclass_path(clustered, "class3").read_text().splitlines()[0]
        assert header.endswith("reason")

    def test_vaf_difference_equal_to_cut_is_tolerated(self, tmp_path):
        clustered = _run(tmp_path, [_vrow(1000, 0.5), _vrow(1001, 0.4)])
        assert _populated(clustered) == ["class1", "class1a"]

    def test_vaf_difference_just_above_cut_declassifies(self, tmp_path):
        clustered = _run(tmp_path, [_vrow(1000, 0.50), _vrow(1001, 0.39)])
        assert _populated(clustered) == ["class3"]

    def test_larger_cut_with_ccf_keeps_the_event(self, tmp_path):
        rows = [_vrow(1000, 0.4), _vrow(1001, 0.6)]
        assert _populated(_run(tmp_path / "vaf", rows)) == ["class3"]
        assert _populated(_run(tmp_path / "ccf", rows, ccf=True)) == ["class1", "class1a"]

    def test_ccf_cut_still_declassifies_large_differences(self, tmp_path):
        clustered = _run(tmp_path, [_vrow(1000, 0.2), _vrow(1001, 0.6)], ccf=True)
        assert _populated(clustered) == ["class3"]

    def test_event_with_two_vaf_populations_is_split_into_events(self, tmp_path):
        rows = [_vrow(1000 + 50 * i, vaf) for i, vaf in enumerate([0.3, 0.3, 0.7, 0.7])]
        clustered = _run(tmp_path, rows)
        assert _populated(clustered) == ["class1", "class1c"]
        assert _positions(clustered, "class1c") == [1000, 1050, 1100, 1150]
        groups = [cols[GROUP_COL] for cols in _entries(clustered, "class1c")]
        assert groups[0] == groups[1] and groups[2] == groups[3] and groups[0] != groups[2]

    def test_mutation_with_another_vaf_is_left_out_of_the_event(self, tmp_path):
        rows = [_vrow(1000, 0.3), _vrow(1050, 0.7), _vrow(1100, 0.3)]
        clustered = _run(tmp_path, rows)
        assert _positions(clustered, "class1c") == [1000, 1100]
        assert _positions(clustered, "class3") == [1050]

    def test_distance_above_imd_cutoff_declassifies_even_with_equal_vaf(self, tmp_path):
        clustered = _run(tmp_path, [_vrow(1000, 0.4), _vrow(1200, 0.4)])
        assert _populated(clustered) == ["class3"]
        assert _positions(clustered, "class3") == [1000, 1200]


# ─── grouping of mutations into events ───────────────────────


class TestEventGrouping:
    def test_distant_mutations_form_separate_events(self, tmp_path):
        rows = [_vrow(1000, 0.4), _vrow(1001, 0.4), _vrow(500000, 0.4), _vrow(500001, 0.4)]
        clustered = _run(tmp_path, rows)
        assert _positions(clustered, "class1a") == [1000, 1001, 500000, 500001]

    def test_group_numbers_increase_with_each_event(self, tmp_path):
        rows = [_vrow(1000, 0.4), _vrow(1001, 0.4), _vrow(500000, 0.4), _vrow(500001, 0.4)]
        groups = [cols[GROUP_COL] for cols in _entries(_run(tmp_path, rows), "class1a")]
        assert groups == ["1", "1", "2", "2"]

    def test_different_chromosomes_are_separate_events(self, tmp_path):
        rows = [
            _vrow(1000, 0.4, chrom="1"),
            _vrow(1001, 0.4, chrom="1"),
            _vrow(1002, 0.4, chrom="2"),
            _vrow(1003, 0.4, chrom="2"),
        ]
        assert len(_entries(_run(tmp_path, rows), "class1a")) == 4

    def test_different_samples_are_separate_events(self, tmp_path):
        rows = [
            _vrow(1000, 0.4, sample="S1"),
            _vrow(1001, 0.4, sample="S1"),
            _vrow(1000, 0.4, sample="S2"),
            _vrow(1001, 0.4, sample="S2"),
        ]
        clustered = _run(tmp_path, rows, imds={"S1": IMD_CUTOFF, "S2": IMD_CUTOFF})
        assert [cols[1] for cols in _entries(clustered, "class1a")] == ["S1", "S1", "S2", "S2"]

    def test_per_sample_imd_cutoff_is_used(self, tmp_path):
        rows = [
            _vrow(1000, 0.4, sample="S1"),
            _vrow(1050, 0.4, sample="S1"),
            _vrow(1000, 0.4, sample="S2"),
            _vrow(1050, 0.4, sample="S2"),
        ]
        clustered = _run(tmp_path, rows, imds={"S1": 100, "S2": 10})
        assert [cols[1] for cols in _entries(clustered, "class1c")] == ["S1", "S1"]
        assert [cols[1] for cols in _entries(clustered, "class3")] == ["S2", "S2"]


# ─── outputs and options ─────────────────────────────────────


class TestOutputs:
    def test_every_subclass_file_has_a_header(self, tmp_path):
        clustered = _run(tmp_path, [_vrow(1000, 0.4), _vrow(1001, 0.4)])
        for name in SUBCLASSES:
            first = _subclass_path(clustered, name).read_text().splitlines()[0]
            assert first.startswith("project\tsamples\tID")
            assert "VAF/CCF" in first

    def test_clusters_of_clusters_imd_table_has_the_class2_events(self, tmp_path):
        rows = [_vrow(1000 + 50 * i, 0.4) for i in range(4)]
        rows += [_vrow(500000 + 50 * i, 0.4) for i in range(4)]
        clustered = _run(tmp_path, rows)
        table = clustered / f"{PROJECT}_clusters_of_clusters_imd.txt"
        lines = [line.split("\t") for line in table.read_text().splitlines() if line.strip()]
        assert lines[0][0] == "clust_group" and lines[0][-1] == "VAF/CCF"
        assert len(lines) == 1 + 8

    def test_previous_subclass_output_is_replaced(self, tmp_path):
        clustered = _run(tmp_path, [_vrow(1000, 0.4), _vrow(1001, 0.4)])
        stale = _subclass_path(clustered, "class1a").parent / "old.txt"
        stale.write_text("stale")
        findClustersOfClusters(
            PROJECT, False, str(tmp_path) + "/", 1000000, {GENOME: {"1": 0}},
            {"S1": []}, None, GENOME, 1, {}, False, False,
        )
        assert not stale.exists()

    def test_original_counts_pickle_counts_every_mutation(self, tmp_path):
        _run(tmp_path, [_vrow(1000, 0.4), _vrow(1001, 0.4), _vrow(1002, 0.4)])
        path = tmp_path / "output" / "simulations" / "data" / "originalCounts.pickle"
        with open(str(path), "rb") as f:
            # three mutations plus the two of the trailing sentinel event
            assert pickle.load(f) == {"S1": 5}

    def test_clustered_stats_pickle_is_written(self, tmp_path):
        clustered = _run(tmp_path, [_vrow(1000, 0.4), _vrow(1001, 0.4)])
        with open(str(clustered / "clusteredStats.pickle"), "rb") as f:
            assert set(pickle.load(f)) == {"I", "II", "III", "Ia", "Ib", "Ic"}

    def test_corrected_paths_are_used_with_correction(self, tmp_path):
        rows = [_vrow(1000, 0.4), _vrow(1001, 0.42)]
        clustered = _run(tmp_path, rows, correction=True)
        assert "vcf_files_corrected" in str(clustered)
        assert _positions(clustered, "class1a") == [1000, 1001]
        counts = tmp_path / "output" / "simulations" / "data" / "originalCounts_corrected.pickle"
        assert os.path.exists(str(counts))

    def test_missing_input_table_raises(self, tmp_path):
        (tmp_path / "output" / "simulations" / "data").mkdir(parents=True)
        with open(str(tmp_path / "output" / "simulations" / "data" / "imds.pickle"), "wb") as f:
            pickle.dump({}, f)
        with pytest.raises(FileNotFoundError):
            findClustersOfClusters(
                PROJECT, False, str(tmp_path) + "/", 1000000, {GENOME: {"1": 0}},
                {}, None, GENOME, 1, {}, False, False,
            )


# ─── last event of the table ─────────────────────────────────


class TestLastEvent:
    @pytest.mark.xfail(
        reason="findClustersOfClusters has no block that closes the last event after the "
        "read loop (findClustersOfClusters_noVAF has one), so the last event of the table "
        "is silently dropped",
        strict=True,
    )
    def test_last_event_of_the_table_is_classified(self, tmp_path):
        rows = [_vrow(1000, 0.4), _vrow(1001, 0.4)]
        clustered = _run(tmp_path, rows, sentinel=False)
        assert _positions(clustered, "class1a") == [1000, 1001]
