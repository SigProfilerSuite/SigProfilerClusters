"""Audit finding 115: when no simulated sample has two mutations on the same chromosome there
are no simulated distances at all; the analysis said so only through a ValueError of
multiprocessing ("Number of processes must be at least 1").
Run with: pytest tests/test_no_simulated_distances.py -v
"""

import inspect

import pytest

from SigProfilerClusters import hotspot
from SigProfilerClusters.exceptions import MissingSimulationsError


class TestCheckSimulatedDistances:
    def test_an_empty_folder_is_a_clear_error(self):
        with pytest.raises(MissingSimulationsError, match="No simulated sample"):
            hotspot.check_simulated_distances([], "/sims/proj_intradistance/")

    def test_hidden_entries_do_not_count(self):
        with pytest.raises(MissingSimulationsError):
            hotspot.check_simulated_distances(
                [".DS_Store", ".DS_Store_intradistance.txt"], "/sims/"
            )

    def test_the_message_says_what_to_do_and_where_it_looked(self):
        with pytest.raises(MissingSimulationsError) as error:
            hotspot.check_simulated_distances([], "/sims/proj_intradistance/")
        message = str(error.value)
        assert "/sims/proj_intradistance/" in message
        assert "more mutations" in message

    def test_samples_with_distances_pass(self):
        assert hotspot.check_simulated_distances(["S1", "S2"], "/sims/") is None


class TestCalledBeforeThePool:
    def test_hotspot_analysis_checks_before_creating_the_pool(self):
        source = inspect.getsource(hotspot.hotSpotAnalysis)
        assert "check_simulated_distances(" in source
        assert source.index("check_simulated_distances(") < source.index("mp.Pool(")
