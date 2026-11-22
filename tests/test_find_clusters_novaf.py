"""Unit tests for classifyFunctions.findClustersOfClusters_noVAF (subclassification of
clustered SBS into DBS, MBS, omikli and kataegis-like classes).
Inputs are small synthetic clustered tables built with tmp_path. The matrix generation
step at the end of the function needs a genome, so it is replaced by a no-op here.
Run with: pytest tests/test_find_clusters_novaf.py -v
"""

import os
import pickle

import pytest

from clustered_table_helpers import (
    GENOME,
    HEADER,
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
from SigProfilerClusters.classifyFunctions import findClustersOfClusters_noVAF


@pytest.fixture(autouse=True)
def stub_matrix_generation(monkeypatch):
    monkeypatch.setattr(classifyFunctions, "generateMatrices", _no_matrices)


def _setup(tmp_path, rows, imds=None, correction=False, chrom_based=False):
    return setup(
        tmp_path, rows, HEADER, f"SNV/{PROJECT}_clustered.txt", imds, correction, chrom_based
    )


def _run(tmp_path, rows, correction=False, chrom_based=False, imds=None):
    clustered = _setup(tmp_path, rows, imds, correction, chrom_based)
    findClustersOfClusters_noVAF(
        PROJECT,
        chrom_based,
        str(tmp_path) + "/",
        1000000,
        {GENOME: {"1": 0, "2": 1000000}},
        {sample: [] for sample in (imds or {"S1": 0})},
        None,
        GENOME,
        1,
        {},
        correction,
    )
    return clustered


# ─── event classes ───────────────────────────────────────────


class TestEventClasses:
    def test_adjacent_pair_is_dbs(self, tmp_path):
        clustered = _run(tmp_path, [_row(1000), _row(1001)])
        assert _populated(clustered) == ["class1", "class1a"]
        assert _positions(clustered, "class1a") == [1000, 1001]

    def test_adjacent_run_of_three_is_mbs(self, tmp_path):
        clustered = _run(tmp_path, [_row(1000), _row(1001), _row(1002)])
        assert _populated(clustered) == ["class1", "class1b"]
        assert _positions(clustered, "class1b") == [1000, 1001, 1002]

    def test_two_mutations_apart_are_omikli(self, tmp_path):
        clustered = _run(tmp_path, [_row(1000), _row(1050)])
        assert _populated(clustered) == ["class1", "class1c"]

    def test_three_mutations_within_cutoff_are_omikli(self, tmp_path):
        clustered = _run(tmp_path, [_row(1000), _row(1050), _row(1100)])
        assert _populated(clustered) == ["class1", "class1c"]
        assert _positions(clustered, "class1c") == [1000, 1050, 1100]

    def test_four_mutations_within_cutoff_are_class2(self, tmp_path):
        rows = [_row(1000), _row(1050), _row(1100), _row(1150)]
        clustered = _run(tmp_path, rows)
        assert "class2" in _populated(clustered)
        assert "class1" not in _populated(clustered)
        assert _positions(clustered, "class2") == [1000, 1050, 1100, 1150]

    def test_distance_above_imd_cutoff_is_class3_and_not_written(self, tmp_path):
        # 200 > IMD_CUTOFF: the event is class3 and no subclass file receives it
        clustered = _run(tmp_path, [_row(1000), _row(1200)])
        assert _populated(clustered) == []

    def test_distance_equal_to_imd_cutoff_is_clustered(self, tmp_path):
        clustered = _run(tmp_path, [_row(1000), _row(1000 + IMD_CUTOFF)])
        assert _populated(clustered) == ["class1", "class1c"]

    def test_class_label_is_the_last_column(self, tmp_path):
        clustered = _run(tmp_path, [_row(1000), _row(1001)])
        assert {cols[-1] for cols in _entries(clustered, "class1")} == {"ClassI"}
        assert {cols[-1] for cols in _entries(clustered, "class1a")} == {"ClassIA"}


class TestKataegisWithOneGap:
    """audit finding 14: four or more mutations with a single distance above 1 are class II"""

    @pytest.mark.parametrize(
        "positions",
        [
            [1000, 1001, 1002, 1004],  # IMDs [1, 1, 2]
            [1000, 1001, 1003, 1004],  # [1, 2, 1]
            [1000, 1001, 1002, 1052],  # [1, 1, 50]
            [1000, 1002, 1003, 1004, 1005],  # [2, 1, 1, 1]
        ],
    )
    def test_four_or_more_with_one_gap_is_class2(self, tmp_path, positions):
        rows = [_row(p) for p in positions] + [_row(500000), _row(500001)]
        clustered = _run(tmp_path, rows)
        assert _positions(clustered, "class2") == positions
        assert "class1c" not in _populated(clustered)

    def test_event_at_the_end_of_the_table_is_class2(self, tmp_path):
        rows = [_row(1000), _row(1001)] + [_row(p) for p in (500000, 500001, 500002, 500004)]
        clustered = _run(tmp_path, rows)
        assert _positions(clustered, "class2") == [500000, 500001, 500002, 500004]

    def test_three_mutations_with_one_gap_stay_omikli(self, tmp_path):
        rows = [_row(1000), _row(1001), _row(1003), _row(500000), _row(500001)]
        clustered = _run(tmp_path, rows)
        assert _positions(clustered, "class1c") == [1000, 1001, 1003]
        assert _positions(clustered, "class2") == []

    def test_four_adjacent_mutations_stay_mbs(self, tmp_path):
        rows = [_row(p) for p in (1000, 1001, 1002, 1003)] + [_row(500000), _row(500001)]
        clustered = _run(tmp_path, rows)
        assert _positions(clustered, "class1b") == [1000, 1001, 1002, 1003]
        assert _positions(clustered, "class2") == []


class TestHasClass2Pattern:
    @pytest.mark.parametrize(
        "gaps, size, expected",
        [(0, 2, False), (0, 5, False), (1, 2, False), (1, 3, False), (1, 4, True),
         (1, 9, True), (2, 3, True), (2, 4, True)],
    )
    def test_pattern(self, gaps, size, expected):
        assert classifyFunctions.has_class2_pattern(gaps, size) is expected


class TestSingleMutationEvent:
    def test_single_mutation_event_in_the_middle_is_dropped(self, tmp_path):
        rows = [_row(1000), _row(500000), _row(900000), _row(900001)]
        clustered = _run(tmp_path, rows)
        assert _positions(clustered, "class1") == [900000, 900001]

    def test_single_mutation_event_at_the_end_is_dropped(self, tmp_path):
        # audit finding 3: it used to be written as class1/class1b (MBS)
        rows = [_row(1000), _row(1001), _row(500000)]
        clustered = _run(tmp_path, rows)
        assert _positions(clustered, "class1") == [1000, 1001]


# ─── class 2 processivity subclasses ─────────────────────────


class TestClass2Processivity:
    def _class2(self, tmp_path, refs):
        rows = [_row(1000 + 50 * i, ref=ref) for i, ref in enumerate(refs)]
        return _run(tmp_path, rows)

    def test_pyrimidine_event_goes_to_2y(self, tmp_path):
        clustered = self._class2(tmp_path, "CCCC")
        assert _positions(clustered, "class2Y") == [1000, 1050, 1100, 1150]
        assert _entries(clustered, "class2K") == []

    def test_keto_event_goes_to_2k(self, tmp_path):
        clustered = self._class2(tmp_path, "TGTG")
        assert len(_entries(clustered, "class2K")) == 4

    def test_strong_event_goes_to_2s(self, tmp_path):
        clustered = self._class2(tmp_path, "CGCG")
        assert len(_entries(clustered, "class2S")) == 4

    def test_mixed_event_goes_to_2n(self, tmp_path):
        clustered = self._class2(tmp_path, "CAGT")
        assert len(_entries(clustered, "class2N")) == 4


# ─── grouping of mutations into events ───────────────────────


class TestEventGrouping:
    def test_distant_mutations_form_separate_events(self, tmp_path):
        rows = [_row(1000), _row(1001), _row(500000), _row(500001)]
        clustered = _run(tmp_path, rows)
        assert _positions(clustered, "class1a") == [1000, 1001, 500000, 500001]

    def test_different_chromosomes_are_separate_events(self, tmp_path):
        rows = [_row(1000, chrom="1"), _row(1001, chrom="1"), _row(1002, chrom="2"), _row(1003, chrom="2")]
        clustered = _run(tmp_path, rows)
        assert len(_entries(clustered, "class1a")) == 4

    def test_different_samples_are_separate_events(self, tmp_path):
        rows = [
            _row(1000, sample="S1"),
            _row(1001, sample="S1"),
            _row(1000, sample="S2"),
            _row(1001, sample="S2"),
        ]
        clustered = _run(tmp_path, rows, imds={"S1": IMD_CUTOFF, "S2": IMD_CUTOFF})
        samples = [cols[1] for cols in _entries(clustered, "class1a")]
        assert samples == ["S1", "S1", "S2", "S2"]

    def test_per_sample_imd_cutoff_is_used(self, tmp_path):
        # a 50 bp distance is clustered for S1 (cutoff 100) but not for S2 (cutoff 10)
        rows = [_row(1000, sample="S1"), _row(1050, sample="S1"), _row(1000, sample="S2"), _row(1050, sample="S2")]
        clustered = _run(tmp_path, rows, imds={"S1": 100, "S2": 10})
        assert [cols[1] for cols in _entries(clustered, "class1c")] == ["S1", "S1"]

    def test_event_written_to_class2_gets_a_group_number(self, tmp_path):
        rows = [_row(1000 + 50 * i) for i in range(4)] + [_row(500000 + 50 * i) for i in range(4)]
        clustered = _run(tmp_path, rows)
        imd_file = clustered / f"{PROJECT}_clusters_of_clusters_imd.txt"
        lines = [l.split("\t") for l in imd_file.read_text().splitlines() if l.strip()]
        assert lines[0][0] == "clust_group"
        assert {row[0] for row in lines[1:]} == {"1", "2"}


# ─── outputs and options ─────────────────────────────────────


class TestOutputs:
    def test_every_subclass_file_has_a_header(self, tmp_path):
        clustered = _run(tmp_path, [_row(1000), _row(1001)])
        for name in SUBCLASSES:
            first = _subclass_path(clustered, name).read_text().splitlines()[0]
            assert first.startswith("project\tsamples\tID")

    def test_class3_header_has_reason_column(self, tmp_path):
        clustered = _run(tmp_path, [_row(1000), _row(1001)])
        header = _subclass_path(clustered, "class3").read_text().splitlines()[0]
        assert header.endswith("reason")

    def test_previous_subclass_output_is_replaced(self, tmp_path):
        clustered = _setup(tmp_path, [_row(1000), _row(1001)])
        stale = _subclass_path(clustered, "class1a")
        stale.parent.mkdir(parents=True)
        (stale.parent / "old.txt").write_text("stale")
        findClustersOfClusters_noVAF(
            PROJECT, False, str(tmp_path) + "/", 1000000, {GENOME: {"1": 0}},
            {"S1": []}, None, GENOME, 1, {}, False,
        )
        assert not (stale.parent / "old.txt").exists()

    def test_original_counts_pickle_counts_mutations_per_sample(self, tmp_path):
        _run(tmp_path, [_row(1000), _row(1001), _row(1002)])
        path = tmp_path / "output" / "simulations" / "data" / "originalCounts.pickle"
        with open(str(path), "rb") as f:
            assert pickle.load(f) == {"S1": 3}

    def test_clustered_stats_pickle_is_written(self, tmp_path):
        clustered = _run(tmp_path, [_row(1000), _row(1001)])
        with open(str(clustered / "clusteredStats.pickle"), "rb") as f:
            stats = pickle.load(f)
        assert set(stats) == {"I", "II", "III", "Ia", "Ib", "Ic"}

    def test_corrected_paths_are_used_with_correction(self, tmp_path):
        clustered = _run(tmp_path, [_row(1000), _row(1001)], correction=True)
        assert "vcf_files_corrected" in str(clustered)
        assert _positions(clustered, "class1a") == [1000, 1001]
        counts = tmp_path / "output" / "simulations" / "data" / "originalCounts_corrected.pickle"
        assert os.path.exists(str(counts))

    def test_chrom_based_uses_chrom_suffix_folders(self, tmp_path):
        clustered = _run(tmp_path, [_row(1000), _row(1001)], chrom_based=True)
        assert _positions(clustered, "class1a", chrom_based=True) == [1000, 1001]
        assert os.path.exists(str(clustered / "clusteredStats_chrom.pickle"))

    def test_missing_input_table_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            findClustersOfClusters_noVAF(
                PROJECT, False, str(tmp_path) + "/", 1000000, {GENOME: {"1": 0}},
                {}, None, GENOME, 1, {}, False,
            )
