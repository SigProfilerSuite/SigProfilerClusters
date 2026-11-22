"""Golden tests of the sample IMD cutoff and of the clustered/non-clustered split
(audit finding #18). The inputs are frozen: simulated background distances and the
observed sample are drawn with the standard library `random` module and fixed seeds, so the
expected cutoffs and mutation identities are exact and reproducible.
The planted clusters of audit finding #1 (clusters with no mutations in the leading IMD bins, which
the cutoff scan used to miss) are checked across several backgrounds.
No genome reference required.
Run with: pytest tests/test_golden_cutoff.py -v
"""

import random

import pytest

from SigProfilerClusters.hotspot import first_run

N_SIMS = 20
SIM_SIZE = 400
GENOME_SCALE = 7.5e6  # mean background distance
CAP = 10000  # the global cutoff never exceeds it


def _background(rng, n):
    return [max(1, int(rng.expovariate(1 / GENOME_SCALE))) for _ in range(n)]


def _simulations(seed):
    rng = random.Random(seed)
    return [_background(rng, SIM_SIZE) for _ in range(N_SIMS)]


def _run(tmp_path, sims, observed, name="run"):
    clustered, non_clustered = tmp_path / name / "c", tmp_path / name / "n"
    clustered.mkdir(parents=True)
    non_clustered.mkdir(parents=True)
    rows = [[str(d), "S1", "1", str(1000 + 10 * i), "C", "T", "1"] for i, d in enumerate(observed)]
    result = first_run(
        sims, rows, observed, str(clustered) + "/", str(non_clustered) + "/",
        "S1", True, N_SIMS, "proj", "GRCh37", False,
    )
    lines = (clustered / "proj_clustered.txt").read_text().splitlines()
    return int(result[7]), sorted(int(line.split("\t")[-1]) for line in lines)


# planted samples -------------------------------------------------------------

def null_sample():
    return _background(random.Random(1), 400)


def adjacent_cluster():
    return [i % 99 + 1 for i in range(300)] + _background(random.Random(2), 100)


def s2_like():
    """40 clustered mutations at 300-495 bp, no adjacent pair, plus background."""
    return [300 + 5 * i for i in range(40)] + _background(random.Random(3), 260)


def s4_like():
    """18 clustered mutations at 100-202 bp, no adjacent pair, plus background."""
    return [100 + 6 * i for i in range(18)] + _background(random.Random(4), 282)


class TestFrozenCutoffs:
    sims = staticmethod(lambda: _simulations(2026))

    def test_null_sample_has_no_clustered_mutations(self, tmp_path):
        cutoff, clustered = _run(tmp_path, self.sims(), null_sample())
        assert (cutoff, clustered) == (4096, [])

    def test_cluster_with_adjacent_mutations_is_called_exactly(self, tmp_path):
        cutoff, clustered = _run(tmp_path, self.sims(), adjacent_cluster())
        assert cutoff == CAP
        assert clustered == sorted(i % 99 + 1 for i in range(300))

    def test_planted_clusters_are_called_with_this_background(self, tmp_path):
        cutoff, clustered = _run(tmp_path, self.sims(), s2_like(), name="s2")
        assert cutoff == CAP
        assert clustered == [300 + 5 * i for i in range(40)]
        _, clustered = _run(tmp_path, self.sims(), s4_like(), name="s4")
        assert clustered == [100 + 6 * i for i in range(18)]

    def test_same_input_gives_the_same_result(self, tmp_path):
        first = _run(tmp_path, self.sims(), adjacent_cluster(), name="first")
        second = _run(tmp_path, self.sims(), adjacent_cluster(), name="second")
        assert first == second


# backgrounds in which the old scan stopped in a leading bin and called nothing
EARLY_STOP_SEEDS = [2, 6, 12, 15, 18, 22, 27, 28]


class TestPlantedClustersWithLeadingSimulatedDistances:
    """With these backgrounds a simulation has a short distance in a bin where the sample has
    none; the scan used to stop there (cutoff 8) and call nothing (audit samples S2 and S4)."""

    sims = staticmethod(lambda: _simulations(6))

    def test_s2_like_cluster_is_detected(self, tmp_path):
        _, clustered = _run(tmp_path, self.sims(), s2_like())
        assert clustered == [300 + 5 * i for i in range(40)]

    def test_s4_like_cluster_is_detected(self, tmp_path):
        _, clustered = _run(tmp_path, self.sims(), s4_like())
        assert clustered == [100 + 6 * i for i in range(18)]

    @pytest.mark.parametrize("seed", EARLY_STOP_SEEDS)
    def test_clusters_are_detected_with_every_background(self, tmp_path, seed):
        _, s2 = _run(tmp_path, _simulations(seed), s2_like(), name="s2")
        _, s4 = _run(tmp_path, _simulations(seed), s4_like(), name="s4")
        assert s2 == [300 + 5 * i for i in range(40)]
        assert s4 == [100 + 6 * i for i in range(18)]

    def test_null_sample_stays_negative_with_this_background(self, tmp_path):
        _, clustered = _run(tmp_path, self.sims(), null_sample())
        assert clustered == []
