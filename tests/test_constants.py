"""The scientific thresholds live in SigProfilerClusters.constants (audit finding 22).
The values equal the literals they replaced (nothing changes for the user), the literals are gone
from the code, and the modules read the constants where they use them.
Run with: pytest tests/test_constants.py -v
"""

import re
from pathlib import Path

import numpy as np
import pytest

import test_find_clusters_novaf as novaf
import test_first_run as first_run_tests
from SigProfilerClusters import classifyFunctions, constants
from SigProfilerClusters.hotspot import densityCorrectionOriginal

PACKAGE = Path(__file__).resolve().parent.parent / "SigProfilerClusters"


class TestValues:
    @pytest.mark.parametrize(
        "name, value",
        [
            ("IMD_CUTOFF_CAP", 10000),
            ("CLUSTERED_FRACTION", 0.90),
            ("CLUSTERED_FRACTION_CORRECTED", 0.95),
            ("Q_VALUE_THRESHOLD", 0.01),
            ("EVENT_GROUPING_DISTANCE", 10001),
            ("VAF_MAX_DIFFERENCE", 0.1),
            ("CCF_MAX_DIFFERENCE", 0.25),
            ("PROCESSIVE_RUN_FRACTION", 0.75),
            ("PROCESSIVE_PURITY", 0.9),
            ("DENSITY_FOLD_CHANGE", 1.25),
            ("DENSITY_MIN_MUTATIONS", 50),
            ("DENSITY_PEAK_THRESHOLD", 9.0),
            ("EVENT_PROBABILITY_FLOOR", 1e-05),
        ],
    )
    def test_value_is_the_historical_literal(self, name, value):
        assert getattr(constants, name) == value


def _code(module):
    """Source lines of a module without comments."""
    return [l for l in (PACKAGE / module).read_text().splitlines() if not l.lstrip().startswith("#")]


class TestLiteralsAreGone:
    @pytest.mark.parametrize(
        "module, pattern",
        [
            ("classifyFunctions.py", r"\b10001\b"),
            ("classifyFunctions.py", r"\* 0\.75"),
            ("classifyFunctions.py", r"purity(CG|TA) > 0\.9\b"),
            ("classifyFunctions.py", r"vaf_cut = 0\.(1|25)\b"),
            ("hotspot.py", r"distance_cut > 10000"),
            ("hotspot.py", r"current_q > 0\.01"),
            ("hotspot.py", r"percent_clustered = 0\.9"),
            ("hotspot.py", r"x > 1\.25"),
            ("hotspot.py", r"threshold = 9\.0"),
            ("hotspot.py", r"keys\(\)\), 10\)"),
            ("SigProfilerClusters.py", r"probs < 1e-05"),
        ],
    )
    def test_literal_is_not_in_the_code(self, module, pattern):
        assert not [l for l in _code(module) if re.search(pattern, l)]


class TestConstantsAreReadWhereTheyAreUsed:
    def test_cutoff_cap(self, tmp_path, monkeypatch):
        rng = np.random.default_rng(5)
        sims = [first_run_tests._background(rng, 400) for _ in range(first_run_tests.SIM_COUNT)]
        monkeypatch.setattr(constants, "IMD_CUTOFF_CAP", 100)
        result = first_run_tests._run(tmp_path, sims, first_run_tests._orig(300, 100))
        assert result[first_run_tests.CUTOFF] == 100

    def test_event_grouping_distance(self, tmp_path, monkeypatch):
        monkeypatch.setattr(classifyFunctions, "generateMatrices", novaf._no_matrices)
        rows = [novaf._row(1000), novaf._row(1050)]
        together = novaf._run(tmp_path / "default", rows)
        assert "class1c" in novaf._populated(together)  # one event of two mutations: omikli
        monkeypatch.setattr(constants, "EVENT_GROUPING_DISTANCE", 5)
        apart = novaf._run(tmp_path / "small", rows)
        assert "class1c" not in novaf._populated(apart)  # two separate events

    def test_density_minimum_mutations(self, monkeypatch):
        bins = np.arange(0, 101, 10)
        sims = {i: [b + 5 for b in range(0, 100, 10)] * 10 for i in range(3)}
        original = [5] * 100 + [b + 5 for b in range(10, 100, 10)] * 5
        assert len(densityCorrectionOriginal(original, sims, bins)) == 2
        monkeypatch.setattr(constants, "DENSITY_MIN_MUTATIONS", 10**6)
        assert densityCorrectionOriginal(original, sims, bins) == [6000000000]
