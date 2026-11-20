"""Unit tests for classifyFunctions.processivitySubclassification and hotspot.refineIMD.
Inputs are small dummy events and simulated distance lists. No genome reference required.
Run with: pytest tests/test_classify_and_imd.py -v
"""

import io

import numpy as np
import pytest

from SigProfilerClusters.classifyFunctions import processivitySubclassification
from SigProfilerClusters.hotspot import refineIMD


def _event(refs):
    """One mutation per reference base; the reference base is column 8."""
    return [["S1", "1", str(100 + i), "x", "x", "x", "x", "x", ref, "x"] for i, ref in enumerate(refs)]


def _classify(refs):
    outs = {name: io.StringIO() for name in ("Y", "K", "S", "N")}
    processivitySubclassification(_event(refs), outs["Y"], outs["K"], outs["S"], outs["N"])
    return {name: out.getvalue() for name, out in outs.items()}


def _only(result, name):
    """Assert that exactly one output received the event and return its text."""
    assert result[name] != ""
    assert all(text == "" for key, text in result.items() if key != name)
    return result[name]


# ─── processivitySubclassification ───────────────────────────


class TestProcessivitySubclassification:
    def test_single_base_event_is_class2y(self):
        _only(_classify("CCCCCC"), "Y")

    def test_pyrimidine_alternation_is_class2y(self):
        _only(_classify("CTCTCT"), "Y")

    def test_long_run_of_one_base_is_class2y_even_below_purity_cutoff(self):
        # purity is 5/6 (< 0.9), but a run of 5 C exceeds 75% of the event length
        _only(_classify("CCCCCT"), "Y")

    def test_keto_alternation_is_class2k(self):
        _only(_classify("TGTGTG"), "K")

    def test_strong_alternation_is_class2s(self):
        _only(_classify("CGCGCG"), "S")

    def test_mixed_event_is_class2n(self):
        _only(_classify("CAGT"), "N")

    def test_output_lines_are_tab_separated_copies_of_the_event(self):
        text = _only(_classify("CCCC"), "Y")
        lines = text.splitlines()
        assert len(lines) == 4
        assert lines[0].split("\t")[8] == "C"
        assert lines[0].split("\t")[2] == "100"
        assert lines[3].split("\t")[2] == "103"

    def test_each_event_goes_to_exactly_one_output(self):
        for refs in ("CCCCCC", "CTCTCT", "TGTGTG", "CGCGCG", "CAGT", "ACGTACGT"):
            result = _classify(refs)
            assert sum(1 for text in result.values() if text) == 1, refs

    def test_all_lines_of_the_event_are_written(self):
        result = _classify("TGTGTG")
        assert len(result["K"].splitlines()) == 6

    def test_empty_event_raises(self):
        with pytest.raises(ZeroDivisionError):
            processivitySubclassification([], io.StringIO(), io.StringIO(), io.StringIO(), io.StringIO())


# ─── refineIMD ───────────────────────────────────────────────


def _sims(seed=0, n_sims=20):
    """Background simulations with a variable number of distances spread over [0, 1000)."""
    rng = np.random.default_rng(seed)
    return [list(rng.integers(0, 1000, size=int(rng.integers(30, 50)))) for _ in range(n_sims)]


def _refine(distances_orig, sims, left=0, right=1000, CI=2):
    empty = np.array([0.0])
    return refineIMD(sims, distances_orig, empty, empty, 0, left, right, CI, [0.0], [0.0])


class TestRefineIMD:
    def test_returns_six_values(self):
        rng = np.random.default_rng(1)
        result = _refine(list(rng.integers(0, 400, size=400)), _sims())
        assert len(result) == 6

    def test_strongly_clustered_sample_keeps_a_high_cutoff(self):
        rng = np.random.default_rng(1)
        clustered = list(rng.integers(0, 400, size=400))
        middle, q_value, avg_sim, std_sim, upper_ci, lower_ci = _refine(clustered, _sims())
        assert 0 < middle <= 1000
        assert q_value < 0.01

    def test_sample_without_excess_collapses_to_the_left_bound(self):
        rng = np.random.default_rng(2)
        background_like = list(rng.integers(0, 1000, size=40))
        middle = _refine(background_like, _sims())[0]
        assert middle == 0

    def test_cutoff_stays_inside_the_search_interval(self):
        rng = np.random.default_rng(3)
        clustered = list(rng.integers(0, 400, size=400))
        middle = _refine(clustered, _sims(), left=100, right=900)[0]
        assert 100 <= middle <= 900

    def test_more_clustering_never_gives_a_smaller_cutoff(self):
        rng = np.random.default_rng(4)
        weak = list(rng.integers(0, 1000, size=60))
        strong = list(rng.integers(0, 200, size=600))
        assert _refine(strong, _sims())[0] >= _refine(weak, _sims())[0]

    def test_average_of_simulations_is_reported(self):
        rng = np.random.default_rng(1)
        sims = _sims()
        result = _refine(list(rng.integers(0, 400, size=400)), sims)
        assert 30 <= result[2] <= 50  # average count in the first bin
