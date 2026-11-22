"""Unit tests for hotspot.calculateSampleIMDs: per-sample IMD cutoffs from the intradistance
files of the original sample and of its simulations, with and without the regional
density correction.
Inputs are synthetic intradistance files built with tmp_path (exponential background, a
block of short clustered distances in the original sample). The density correction is
replaced where a known dense region is needed. No genome reference required.
Run with: pytest tests/test_calculate_sample_imds.py -v
"""

import numpy as np
import pytest

from SigProfilerClusters import hotspot
from SigProfilerClusters.hotspot import calculateSampleIMDs

GENOME = "GRCh37"
WINDOW = 1_000_000
EDGE = 1_000_000
N_SIMS = 12  # densityCorrection samples 10 simulations
CHROM_LENGTHS = {GENOME: {"1": 0}}
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
    "regions",
    "regionsSamps",
    "avg_simCounts_corrected",
    "std_simCounts_corrected",
    "upper_CIs_refined_corrected",
    "lower_CIs_refined_corrected",
]


def _build(root, sample="S1", clustered=300, background=100, exome="", seed=5, orig=True):
    """Simulation files under sims/<sample>/ and the original sample's intradistance file."""
    rng = np.random.default_rng(seed)
    sims = root / "sims" / sample
    sims.mkdir(parents=True, exist_ok=True)
    for i in range(N_SIMS):
        distances = [max(1, int(x)) for x in rng.exponential(7.5e6, size=400)]
        positions = rng.integers(1000, 90000, size=400)
        text = "".join(
            f"{d}\t{sample}_{i}\t1\t{p}\tC\tT\t{d}\n" for d, p in zip(distances, positions)
        )
        (sims / f"{sample}_{i}_96{exome}_intradistance.txt").write_text(text)
    (root / "orig").mkdir(exist_ok=True)
    (root / "c").mkdir(exist_ok=True)
    (root / "n").mkdir(exist_ok=True)
    if orig:
        short = [int(x) for x in rng.integers(1, 100, size=clustered)]
        far = [max(1, int(x)) for x in rng.exponential(7.5e6, size=background)]
        positions = rng.integers(1000, 90000, size=clustered + background)
        text = "".join(
            f"{d}\t{sample}\t1\t{p}\tC\tT\t{d}\n" for d, p in zip(short + far, positions)
        )
        (root / "orig" / f"{sample}_96{exome}_intradistance.txt").write_text(text)


def _run(root, correction=False, folders=("S1",), exome=""):
    values = calculateSampleIMDs(
        "proj",
        list(folders),
        str(root / "sims") + "/",
        str(root / "orig") + "/",
        str(root / "c") + "/",
        str(root / "n") + "/",
        True,
        GENOME,
        WINDOW,
        False,
        "96",
        exome,
        False,
        correction,
        CHROM_LENGTHS,
    )
    return dict(zip(RESULT_NAMES, values))


def _regions(monkeypatch, regions):
    monkeypatch.setattr(hotspot, "densityCorrection", lambda *args, **kwargs: list(regions))


def _table_rows(root, name):
    path = root / name / f"proj_{'clustered' if name == 'c' else 'nonClustered'}.txt"
    return path.read_text().splitlines() if path.exists() else None


# ─── without correction ──────────────────────────────────────


class TestWithoutCorrection:
    def test_returns_all_the_result_structures(self, tmp_path):
        _build(tmp_path)
        result = _run(tmp_path)
        assert list(result) == RESULT_NAMES

    def test_sample_gets_a_cutoff_and_a_significant_excess(self, tmp_path):
        _build(tmp_path)
        result = _run(tmp_path)
        assert 0 < result["imds"]["S1"] <= 10000
        assert result["q_values"]["S1"] < 0.01

    def test_clustered_and_non_clustered_tables_are_written(self, tmp_path):
        _build(tmp_path)
        _run(tmp_path)
        assert len(_table_rows(tmp_path, "c")) == 300
        assert len(_table_rows(tmp_path, "n")) == 100

    def test_no_correction_structures_are_filled(self, tmp_path):
        _build(tmp_path)
        result = _run(tmp_path)
        assert result["imds_corrected"] == {}
        assert result["regions"] == []
        assert result["regionsSamps"] == {"S1": []}

    def test_every_per_sample_structure_has_the_sample(self, tmp_path):
        _build(tmp_path)
        result = _run(tmp_path)
        for name in ("y2s", "bincenters2s", "interval_lines", "orig_mutations_samps", "lower_CIs", "upper_CIs"):
            assert set(result[name]) == {"S1"}, name

    def test_exome_suffix_is_used_in_the_file_names(self, tmp_path):
        _build(tmp_path, exome="_exome")
        result = _run(tmp_path, exome="_exome")
        assert result["imds"]["S1"] > 0

    def test_samples_are_processed_independently(self, tmp_path):
        _build(tmp_path, sample="S1", seed=5)
        _build(tmp_path, sample="S2", seed=6)
        result = _run(tmp_path, folders=("S1", "S2"))
        assert set(result["imds"]) == {"S1", "S2"}
        samples = [row.split("\t")[1] for row in _table_rows(tmp_path, "c")]
        assert samples.count("S1") == 300 and samples.count("S2") == 300


class TestSkippedSamples:
    def test_hidden_entries_are_ignored(self, tmp_path):
        _build(tmp_path)
        result = _run(tmp_path, folders=(".DS_Store",))
        assert result["imds"] == {}

    def test_sample_without_an_original_file_is_skipped(self, tmp_path):
        _build(tmp_path, sample="S1")
        _build(tmp_path, sample="S2", orig=False)
        result = _run(tmp_path, folders=("S1", "S2"))
        assert set(result["imds"]) == {"S1"}

    def test_sample_without_distances_is_skipped(self, tmp_path):
        _build(tmp_path)
        (tmp_path / "orig" / "S1_96_intradistance.txt").write_text("0\tS1\t1\t100\tC\tT\t0\n")
        assert _run(tmp_path)["imds"] == {}

    def test_sample_without_simulation_files_is_skipped(self, tmp_path):
        _build(tmp_path, sample="S1")
        for path in (tmp_path / "sims" / "S1").iterdir():
            path.unlink()
        assert _run(tmp_path)["imds"] == {}


# ─── with the regional density correction ────────────────────


class TestWithCorrection:
    def test_without_dense_regions_only_the_sample_level_result_is_filled(self, tmp_path, monkeypatch):
        _build(tmp_path)
        _regions(monkeypatch, [])
        result = _run(tmp_path, correction=True)
        assert result["imds"]["S1"] > 0
        assert result["regionsSamps"] == {"S1": []}
        assert result["imds_corrected"] == {}

    def test_dense_region_gets_its_own_cutoff(self, tmp_path, monkeypatch):
        _build(tmp_path)
        _regions(monkeypatch, [EDGE])
        result = _run(tmp_path, correction=True)
        assert result["regionsSamps"] == {"S1": [EDGE]}
        assert set(result["imds_corrected"]["S1"]) == {EDGE}
        assert 0 < result["imds_corrected"]["S1"][EDGE] <= 10000

    def test_region_results_are_stored_per_region(self, tmp_path, monkeypatch):
        _build(tmp_path)
        _regions(monkeypatch, [EDGE])
        result = _run(tmp_path, correction=True)
        for name in ("y2s_corrected", "bincenters2s_corrected", "q_values_corrected", "interval_lines_corrected"):
            assert set(result[name]["S1"]) == {EDGE}, name

    def test_region_without_mutations_is_skipped(self, tmp_path, monkeypatch):
        _build(tmp_path)
        _regions(monkeypatch, [EDGE, 2 * EDGE])
        result = _run(tmp_path, correction=True)
        assert set(result["imds_corrected"]["S1"]) == {EDGE}
        assert result["regionsSamps"]["S1"] == [EDGE, 2 * EDGE]

    def test_regional_cutoff_does_not_duplicate_the_tables(self, tmp_path, monkeypatch):
        _build(tmp_path)
        _regions(monkeypatch, [EDGE])
        _run(tmp_path, correction=True)
        assert len(_table_rows(tmp_path, "c")) == 300
        assert len(_table_rows(tmp_path, "n")) == 100

    def test_real_density_correction_runs_with_enough_simulations(self, tmp_path):
        _build(tmp_path)
        result = _run(tmp_path, correction=True)
        assert result["imds"]["S1"] > 0
        assert "S1" in result["regionsSamps"]
