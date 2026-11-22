"""Unit tests for the genome-wide density correction (correction=True) inside
classifyFunctions.findClustersOfClusters and findClustersOfClusters_noVAF.
Every test runs on both functions with the same synthetic table. Dense regions are given
as `regions` (sample -> sorted cumulative positions) and `imds_corrected` (sample ->
{region edge: corrected IMD cutoff}); the global sample cutoff is 100 bp.
The matrix generation step needs a genome and is stubbed.
Tests marked xfail record behavior that looks wrong and is planned for a later fix.
Run with: pytest tests/test_find_clusters_correction.py -v
"""

import pytest

from clustered_table_helpers import (
    GENOME,
    HEADER,
    HEADER_VAF,
    PROJECT,
    SUBCLASSES,
    entries,
    no_matrices,
    row,
    setup,
)
from SigProfilerClusters import classifyFunctions
from SigProfilerClusters.classifyFunctions import (
    findClustersOfClusters,
    findClustersOfClusters_noVAF,
)

GLOBAL_CUTOFF = 100
EDGE = 1_000_000  # right edge of the dense region, in cumulative genome coordinates
WINDOW = 1_000_000
CHROM_LENGTHS = {GENOME: {"1": 0, "2": 5_000_000}}
VARIANTS = ["novaf", "vaf"]


@pytest.fixture(autouse=True)
def stub_matrix_generation(monkeypatch):
    monkeypatch.setattr(classifyFunctions, "generateMatrices", no_matrices)


def _mutation(variant, pos, imd=1, sample="S1", chrom="1"):
    """imd is the distance to the previous mutation, as stored in the clustered table."""
    return row(pos, sample=sample, chrom=chrom, imd=imd, vaf=0.4 if variant == "vaf" else None)


def _run(
    tmp_path,
    variant,
    rows,
    regions=None,
    imds_corrected=None,
    correction=True,
    samples=("S1",),
):
    """Run one function. Both end with an adjacent pair on chromosome 2 so that the last
    event of the table (which findClustersOfClusters ignores) is never one under test."""
    rows = list(rows) + [_mutation(variant, 1000, chrom="2"), _mutation(variant, 1001, chrom="2")]
    imds = {sample: GLOBAL_CUTOFF for sample in samples}
    regions = {sample: [] for sample in samples} if regions is None else regions
    if variant == "vaf":
        header, table = HEADER_VAF, f"{PROJECT}_clustered_vaf.txt"
    else:
        header, table = HEADER, f"SNV/{PROJECT}_clustered.txt"
    clustered = setup(tmp_path, rows, header, table, imds, correction, False, imds_corrected)
    args = (
        PROJECT,
        False,
        str(tmp_path) + "/",
        WINDOW,
        CHROM_LENGTHS,
        regions,
        None,
        GENOME,
        1,
        {},
        correction,
    )
    if variant == "vaf":
        findClustersOfClusters(*args, False)
    else:
        findClustersOfClusters_noVAF(*args)
    return clustered


def _positions(clustered, name):
    """Positions on chromosome 1 reported in one subclass file."""
    return [int(cols[6]) for cols in entries(clustered, name) if cols[5] == "1"]


def _classified(clustered):
    """Subclass files that received chromosome 1 mutations."""
    return sorted(name for name in SUBCLASSES if _positions(clustered, name))


DENSE = {"S1": [EDGE]}


def _dense(cutoff):
    return {"S1": {EDGE: cutoff}}


# ─── mutations grouped by the corrected cutoff ───────────────


@pytest.mark.parametrize("variant", VARIANTS)
class TestCorrectedGrouping:
    def test_distant_mutations_in_a_dense_region_form_an_event(self, tmp_path, variant):
        rows = [_mutation(variant, 1000), _mutation(variant, 31000, imd=30000)]
        clustered = _run(tmp_path, variant, rows, DENSE, _dense(50000))
        assert _classified(clustered) == ["class1", "class1c"]
        assert _positions(clustered, "class1c") == [1000, 31000]

    def test_the_same_mutations_are_not_an_event_without_correction(self, tmp_path, variant):
        rows = [_mutation(variant, 1000), _mutation(variant, 31000, imd=30000)]
        clustered = _run(tmp_path, variant, rows, DENSE, _dense(50000), correction=False)
        assert _classified(clustered) == []

    def test_distance_above_the_corrected_cutoff_is_not_grouped(self, tmp_path, variant):
        rows = [_mutation(variant, 1000), _mutation(variant, 61000, imd=60000)]
        clustered = _run(tmp_path, variant, rows, DENSE, _dense(50000))
        assert _classified(clustered) == []

    def test_region_without_a_corrected_cutoff_uses_the_global_one(self, tmp_path, variant):
        rows = [_mutation(variant, 1000), _mutation(variant, 31000, imd=30000)]
        clustered = _run(tmp_path, variant, rows, DENSE, {"S1": {}})
        assert _classified(clustered) == []

    def test_corrected_cutoff_is_applied_per_sample(self, tmp_path, variant):
        rows = [
            _mutation(variant, 1000, sample="S1"),
            _mutation(variant, 31000, imd=30000, sample="S1"),
            _mutation(variant, 1000, sample="S2"),
            _mutation(variant, 31000, imd=30000, sample="S2"),
        ]
        regions = {"S1": [EDGE], "S2": []}
        clustered = _run(
            tmp_path, variant, rows, regions, {"S1": {EDGE: 50000}, "S2": {}}, samples=("S1", "S2")
        )
        samples = [cols[1] for cols in entries(clustered, "class1c") if cols[5] == "1"]
        assert samples == ["S1", "S1"]

    def test_original_counts_use_the_corrected_suffix(self, tmp_path, variant):
        rows = [_mutation(variant, 1000), _mutation(variant, 31000, imd=30000)]
        _run(tmp_path, variant, rows, DENSE, _dense(50000))
        data = tmp_path / "output" / "simulations" / "data"
        assert (data / "originalCounts_corrected.pickle").exists()
        assert not (data / "originalCounts.pickle").exists()


# ─── classification of grouped events under correction ───────


@pytest.mark.parametrize("variant", VARIANTS)
class TestCorrectedClassification:
    def test_pair_within_the_corrected_cutoff_is_omikli(self, tmp_path, variant):
        rows = [_mutation(variant, 1000), _mutation(variant, 4000, imd=3000)]
        clustered = _run(tmp_path, variant, rows, DENSE, _dense(5000))
        assert _classified(clustered) == ["class1", "class1c"]

    def test_adjacent_pair_is_dbs_in_a_dense_region(self, tmp_path, variant):
        rows = [_mutation(variant, 1000), _mutation(variant, 1001)]
        clustered = _run(tmp_path, variant, rows, DENSE, _dense(5000))
        assert _classified(clustered) == ["class1", "class1a"]

    def test_pair_within_the_global_cutoff_is_omikli(self, tmp_path, variant):
        rows = [_mutation(variant, 1000), _mutation(variant, 1050, imd=50)]
        clustered = _run(tmp_path, variant, rows, DENSE, _dense(5000))
        assert _classified(clustered) == ["class1", "class1c"]

    def test_four_mutations_in_a_dense_region_are_class2(self, tmp_path, variant):
        rows = [_mutation(variant, 1000)] + [
            _mutation(variant, 1000 + 3000 * i, imd=3000) for i in range(1, 4)
        ]
        clustered = _run(tmp_path, variant, rows, DENSE, _dense(5000))
        assert "class2" in _classified(clustered)
        assert _positions(clustered, "class2") == [1000, 4000, 7000, 10000]


# ─── behavior that looks wrong ───────────────────────────────


@pytest.mark.parametrize("variant", VARIANTS)
class TestSuspectedProblems:
    @pytest.mark.xfail(
        reason="with correction=True a distance above the global cutoff is not counted as "
        "failed when the sample has no dense region (distancesFailed requires regions), "
        "so the pair becomes omikli; correction=False declassifies it",
        strict=True,
    )
    def test_distance_above_the_global_cutoff_without_regions_is_not_an_event(
        self, tmp_path, variant
    ):
        rows = [_mutation(variant, 1000), _mutation(variant, 1200, imd=200)]
        clustered = _run(tmp_path, variant, rows)
        assert _positions(clustered, "class1c") == []

    @pytest.mark.xfail(
        reason="a distance above both the global and the corrected cutoff is still "
        "classified as omikli: distancesFailed passes the per-sample imds_corrected "
        "dictionary to cutoffCatch, which looks the sample up inside it",
        strict=True,
    )
    def test_distance_above_the_corrected_cutoff_is_not_an_event(self, tmp_path, variant):
        rows = [_mutation(variant, 1000), _mutation(variant, 9000, imd=8000)]
        clustered = _run(tmp_path, variant, rows, DENSE, _dense(5000))
        assert _positions(clustered, "class1c") == []

    @pytest.mark.xfail(
        reason="to verify: cutoffCatch compares the cumulative genomic position with "
        "windowSize, so the corrected cutoff is only applied in the first window of the genome",
        strict=True,
    )
    def test_corrected_cutoff_applies_beyond_the_first_window(self, tmp_path, variant):
        rows = [
            _mutation(variant, 2_000_000),
            _mutation(variant, 2_030_000, imd=30000),
        ]
        clustered = _run(
            tmp_path, variant, rows, {"S1": [3_000_000]}, {"S1": {3_000_000: 50000}}
        )
        assert _positions(clustered, "class1c") == [2_000_000, 2_030_000]
