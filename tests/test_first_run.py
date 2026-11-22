"""Unit tests for hotspot.first_run: the per-sample IMD cutoff and the split of the original
mutations into clustered and non-clustered tables.
Inputs are simulated inter-mutation distances: the background is exponential (mean of the
order of a genome divided by the mutation count, so short distances are rare) and the
original sample mixes short clustered distances with background-like ones.
No genome reference required.
Run with: pytest tests/test_first_run.py -v
"""

import numpy as np
import pytest

from SigProfilerClusters.hotspot import first_run

PROJECT = "proj"
GENOME = "GRCh37"
SIM_COUNT = 20
EDGE = 1_000_000


def _background(rng, n):
    return [max(1, int(x)) for x in rng.exponential(7.5e6, size=n)]


@pytest.fixture(scope="module")
def sims():
    rng = np.random.default_rng(5)
    return [_background(rng, 400) for _ in range(SIM_COUNT)]


def _orig(n_clustered, n_background, seed=11):
    rng = np.random.default_rng(seed)
    short = [int(x) for x in rng.integers(1, 100, size=n_clustered)]
    return short + _background(rng, n_background)


def _lines(distances, sample="S1"):
    """Rows as read from an intradistance file: IMD, sample, chr, position, ref, alt, plot."""
    return [
        [str(d), sample, "1", str(1000 + 10 * i), "C", "T", "1"]
        for i, d in enumerate(distances)
    ]


def _run(tmp_path, sims, distances, vaf=False, sample="S1", **kwargs):
    clust, non_clust = tmp_path / "clustered", tmp_path / "nonClustered"
    clust.mkdir(exist_ok=True)
    non_clust.mkdir(exist_ok=True)
    result = first_run(
        sims,
        _lines(distances, sample),
        distances,
        str(clust) + "/",
        str(non_clust) + "/",
        sample,
        True,
        SIM_COUNT,
        PROJECT,
        GENOME,
        vaf,
        **kwargs,
    )
    return result


def _table(tmp_path, kind):
    name = {"clustered": f"{PROJECT}_clustered.txt", "nonClustered": f"{PROJECT}_nonClustered.txt"}[kind]
    path = tmp_path / kind / name
    return [line.split("\t") for line in path.read_text().splitlines()] if path.exists() else None


# result tuple positions
Q_VAL, N_CLUSTERED, CUTOFF = 2, 4, 7


# ─── cutoff and split of a clustered sample ──────────────────


class TestClusteredSample:
    def test_returns_thirteen_values(self, tmp_path, sims):
        assert len(_run(tmp_path, sims, _orig(300, 100))) == 13

    def test_short_distances_are_clustered_and_background_is_not(self, tmp_path, sims):
        result = _run(tmp_path, sims, _orig(300, 100))
        assert result[N_CLUSTERED] == 300
        assert len(_table(tmp_path, "clustered")) == 300
        assert len(_table(tmp_path, "nonClustered")) == 100

    def test_excess_over_the_background_is_significant(self, tmp_path, sims):
        assert _run(tmp_path, sims, _orig(300, 100))[Q_VAL] < 0.01

    def test_cutoff_is_capped_at_10000(self, tmp_path, sims):
        assert _run(tmp_path, sims, _orig(300, 100))[CUTOFF] == 10000

    def test_clustered_rows_are_within_the_cutoff_and_the_rest_above_it(self, tmp_path, sims):
        cutoff = _run(tmp_path, sims, _orig(300, 100))[CUTOFF]
        assert all(int(row[-1]) <= cutoff for row in _table(tmp_path, "clustered"))
        assert all(int(row[-1]) > cutoff for row in _table(tmp_path, "nonClustered"))

    def test_every_mutation_goes_to_exactly_one_table(self, tmp_path, sims):
        _run(tmp_path, sims, _orig(300, 100))
        assert len(_table(tmp_path, "clustered")) + len(_table(tmp_path, "nonClustered")) == 400

    def test_row_layout(self, tmp_path, sims):
        _run(tmp_path, sims, _orig(300, 100))
        row = _table(tmp_path, "clustered")[0]
        assert row[:5] == [PROJECT, "S1", ".", GENOME, "SNP"]
        assert row[5:10] == ["1", "1000", "1000", "C", "T"]
        assert row[10] == "SOMATIC"
        assert len(row) == 13

    def test_vaf_column_is_added_for_subclassification(self, tmp_path, sims):
        _run(tmp_path, sims, _orig(300, 100), vaf=True)
        assert len(_table(tmp_path, "clustered")[0]) == 14

    def test_tables_are_appended_for_each_sample(self, tmp_path, sims):
        _run(tmp_path, sims, _orig(300, 100), sample="S1")
        _run(tmp_path, sims, _orig(300, 100, seed=12), sample="S2")
        samples = [row[1] for row in _table(tmp_path, "clustered")]
        assert samples.count("S1") == 300 and samples.count("S2") == 300


class TestSampleWithoutClustering:
    def test_background_like_sample_has_no_clustered_mutations(self, tmp_path, sims):
        result = _run(tmp_path, sims, _orig(0, 400, seed=2))
        assert result[N_CLUSTERED] == 0
        assert _table(tmp_path, "clustered") == []
        assert len(_table(tmp_path, "nonClustered")) == 400

    def test_cutoff_does_not_reach_any_observed_distance(self, tmp_path, sims):
        # the scan starts at the first IMD bin with observed mutations (audit finding 1), so
        # for a background-like sample the cutoff stays below its shortest distance
        distances = _orig(0, 400, seed=2)
        assert _run(tmp_path, sims, distances)[CUTOFF] <= min(distances)


# ─── correctionData mode ─────────────────────────────────────


class TestCorrectionDataMode:
    def test_no_tables_are_written(self, tmp_path, sims):
        _run(tmp_path, sims, _orig(300, 100), correctionData=True)
        assert _table(tmp_path, "clustered") is None
        assert _table(tmp_path, "nonClustered") is None

    def test_clustered_count_is_still_returned(self, tmp_path, sims):
        assert _run(tmp_path, sims, _orig(300, 100), correctionData=True)[N_CLUSTERED] == 300


# ─── corrected IMD of dense regions ──────────────────────────


class TestRegionCorrectedCutoff:
    chrom_lengths = {GENOME: {"1": 0}}
    far = 20000  # above the 10000 cap of the sample cutoff: only a regional cutoff can cluster it

    def _distances(self):
        # background without any short distance, so that only the planted ones can be clustered
        rng = np.random.default_rng(3)
        background = [d for d in _background(rng, 150) if d > 50000][:100]
        return background + [self.far] * 20

    def _run(self, tmp_path, sims, regions, imds_corrected):
        return _run(
            tmp_path,
            sims,
            self._distances(),
            correction=True,
            regions=regions,
            imds_corrected=imds_corrected,
            windowSize=1_000_000,
            chromLengths=self.chrom_lengths,
        )

    def test_without_correction_distant_mutations_are_not_clustered(self, tmp_path, sims):
        assert _run(tmp_path, sims, self._distances())[N_CLUSTERED] == 0

    def test_region_cutoff_clusters_mutations_above_the_sample_cutoff(self, tmp_path, sims):
        result = self._run(tmp_path, sims, [EDGE], {"S1": {EDGE: 30000}})
        assert result[N_CLUSTERED] == 20
        assert {int(row[-1]) for row in _table(tmp_path, "clustered")} == {self.far}

    def test_region_cutoff_below_the_distance_does_not_cluster(self, tmp_path, sims):
        assert self._run(tmp_path, sims, [EDGE], {"S1": {EDGE: 15000}})[N_CLUSTERED] == 0

    def test_without_regions_the_sample_cutoff_is_used(self, tmp_path, sims):
        assert self._run(tmp_path, sims, [], {"S1": {}})[N_CLUSTERED] == 0

    def test_region_without_a_corrected_cutoff_uses_the_sample_cutoff(self, tmp_path, sims):
        assert self._run(tmp_path, sims, [EDGE], {"S1": {}})[N_CLUSTERED] == 0
