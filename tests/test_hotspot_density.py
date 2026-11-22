"""Unit tests for hotspot.densityCorrection (search for genomic regions denser than the
simulations) and hotspot.getResults (merge of the per-process results).
Inputs are simulated mutation positions. No genome reference required.
Tests marked xfail record behavior that looks wrong and is planned for a later fix.
Run with: pytest tests/test_hotspot_density.py -v
"""

import warnings

import pytest

from SigProfilerClusters import hotspot
from SigProfilerClusters.hotspot import densityCorrection, getResults

WINDOW = 1000


def _simulations(n=12, windows=10, per_window=10):
    """Every simulation has the same number of mutations in each of the first windows."""
    positions = [w * WINDOW + 500 for w in range(windows)] * per_window
    return {str(i): list(positions) for i in range(n)}


def _correction(muts, sims):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return densityCorrection(muts, sims, WINDOW)


# ─── densityCorrection ───────────────────────────────────────


class TestDensityCorrection:
    def test_returns_a_list(self):
        muts = [w * WINDOW + 500 for w in range(10)] * 10
        assert isinstance(_correction(muts, _simulations()), list)

    def test_sample_as_dense_as_the_simulations_has_no_regions(self):
        muts = [w * WINDOW + 500 for w in range(10)] * 10
        assert _correction(muts, _simulations()) == []

    def test_sample_without_mutations_has_no_regions(self):
        assert _correction([], _simulations()) == []

    def test_simulations_without_mutations_have_no_regions(self):
        assert _correction([3 * WINDOW + 500] * 200, {str(i): [] for i in range(12)}) == []

    def test_window_far_denser_than_the_simulations_is_a_region(self):
        muts = [3 * WINDOW + 500] * 500 + [w * WINDOW + 500 for w in range(10) if w != 3] * 10
        assert len(_correction(muts, _simulations())) == 1

    def test_regions_depend_on_the_original_mutations(self):
        first = _correction([3 * WINDOW + 500] * 5000, _simulations())
        second = _correction([9 * WINDOW + 500] * 5000, _simulations())
        assert first != second

    def test_fewer_than_ten_simulations(self):
        assert isinstance(_correction([3 * WINDOW + 500] * 200, _simulations(n=5)), list)

    def test_region_is_the_right_edge_of_the_dense_window(self):
        muts = [3 * WINDOW + 500] * 500 + [w * WINDOW + 500 for w in range(10) if w != 3] * 10
        assert _correction(muts, _simulations()) == [4 * WINDOW]

    def test_dense_first_window_is_found(self):
        muts = [500] * 500 + [w * WINDOW + 500 for w in range(1, 10)] * 10
        assert _correction(muts, _simulations()) == [WINDOW]

    def test_dense_window_far_along_the_genome_is_found(self):
        muts = [50 * WINDOW + 500] * 500
        assert _correction(muts, _simulations()) == [51 * WINDOW]

    def test_uses_all_simulations_and_is_deterministic(self):
        muts = [3 * WINDOW + 500] * 500 + [w * WINDOW + 500 for w in range(10) if w != 3] * 10
        sims = _simulations(n=40)
        assert _correction(muts, sims) == _correction(muts, dict(reversed(sims.items())))


# ─── getResults ──────────────────────────────────────────────

RESULT_NAMES = [
    "imds",
    "y2s",
    "bincenters2s",
    "q_values",
    "interval_lines",
    "orig_mutations_samps",
    "lower_CIs",
    "upper_CIs",
    "avg_bin_counts_samp",
    "avg_simCounts",
    "std_simCounts",
    "upper_CIs_refined",
    "lower_CIs_refined",
    "imds_corrected",
    "y2s_corrected",
    "bincenters2s_corrected",
    "q_values_corrected",
    "interval_lines_corrected",
    "orig_mutations_samps_corrected",
    "lower_CIs_corrected",
    "upper_CIs_corrected",
    "avg_bin_counts_samp_corrected",
    "regions",  # a list that getResults does not keep
    "regionsSamps",
    "avg_simCounts_corrected",
    "std_simCounts_corrected",
    "upper_CIs_refined_corrected",
    "lower_CIs_refined_corrected",
]
STORED = [name for name in RESULT_NAMES if name != "regions"]


@pytest.fixture
def finals(monkeypatch):
    """The module-level dictionaries that getResults updates (created by hotSpotAnalysis)."""
    store = {}
    for name in STORED:
        store[name] = {}
        monkeypatch.setattr(hotspot, name + "Final", store[name], raising=False)
    return store


def _result(sample, value):
    return tuple([10, 20] if name == "regions" else {sample: (name, value)} for name in RESULT_NAMES)


class TestGetResults:
    def test_each_result_updates_its_own_dictionary(self, finals):
        getResults(_result("S1", 1))
        for name in STORED:
            assert finals[name] == {"S1": (name, 1)}

    def test_results_of_several_samples_accumulate(self, finals):
        getResults(_result("S1", 1))
        getResults(_result("S2", 2))
        for name in STORED:
            assert set(finals[name]) == {"S1", "S2"}

    def test_a_sample_reported_twice_is_overwritten(self, finals):
        getResults(_result("S1", 1))
        getResults(_result("S1", 2))
        assert finals["imds"] == {"S1": ("imds", 2)}

    def test_wrong_number_of_values_raises(self, finals):
        with pytest.raises(ValueError):
            getResults(_result("S1", 1)[:-1])
