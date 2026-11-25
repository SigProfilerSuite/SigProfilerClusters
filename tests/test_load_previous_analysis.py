"""Unit tests for SigProfilerClusters.load_previous_analysis and for analysis="subClassify"
(audit finding 9: regions and imds were never defined when the hotspot analysis is skipped).
Run with: pytest tests/test_load_previous_analysis.py -v
"""

import pickle

import pytest

from SigProfilerClusters import SigProfilerClusters as spc
from SigProfilerClusters.SigProfilerClusters import load_previous_analysis, analysis

PROJECT = "proj"
REGIONS = {"S1": [1_000_000]}
IMDS = {"S1": {1_000_000: 5000}}


def _dump(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "wb") as f:
        pickle.dump(value, f)


def _previous(tmp_path, correction=True, chrom_based=False):
    (tmp_path / "input").mkdir(exist_ok=True)
    (tmp_path / "input" / "x.vcf").write_text(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "1\t255000\t.\tG\tT\t.\t.\tVAF=0.5\n"
        "1\t255001\t.\tA\tC\t.\t.\tVAF=0.5\n"
    )
    data = tmp_path / "output" / "simulations" / "data"
    suffix = "_chrom" if chrom_based else ""
    _dump(data / f"imds{suffix}.pickle", {"S1": 100})
    if correction:
        _dump(data / f"regions_corrected{suffix}.pickle", REGIONS)
        _dump(data / f"imds_corrected{suffix}.pickle", IMDS)
    folder = "vcf_files_corrected" if correction else "vcf_files"
    table = tmp_path / "output" / folder / f"{PROJECT}_clustered" / "SNV" / f"{PROJECT}_clustered.txt"
    table.parent.mkdir(parents=True, exist_ok=True)
    table.write_text("header\n")


def _load(tmp_path, contexts="96", correction=True, chrom_based=False):
    return load_previous_analysis(PROJECT, str(tmp_path) + "/", contexts, correction, chrom_based)


class TestLoadPreviousAnalysis:
    def test_loads_the_regions_and_the_corrected_cutoffs(self, tmp_path):
        _previous(tmp_path)
        assert _load(tmp_path) == (REGIONS, IMDS)

    def test_without_correction_both_are_empty(self, tmp_path):
        _previous(tmp_path, correction=False)
        assert _load(tmp_path, correction=False) == ({}, {})

    def test_missing_files_are_all_listed(self, tmp_path):
        with pytest.raises(FileNotFoundError) as error:
            _load(tmp_path)
        message = str(error.value)
        for name in ("imds.pickle", "imds_corrected.pickle", "regions_corrected.pickle", "_clustered.txt"):
            assert name in message
        assert "analysis=\"hotspot\"" in message

    def test_one_missing_file_is_named(self, tmp_path):
        _previous(tmp_path)
        (tmp_path / "output" / "simulations" / "data" / "regions_corrected.pickle").unlink()
        with pytest.raises(FileNotFoundError, match="regions_corrected.pickle"):
            _load(tmp_path)

    def test_without_correction_the_corrected_files_are_not_required(self, tmp_path):
        _previous(tmp_path, correction=False)
        assert _load(tmp_path, correction=False) == ({}, {})

    def test_indel_context_does_not_need_the_clustered_table(self, tmp_path):
        _previous(tmp_path, correction=False)
        (tmp_path / "output" / "vcf_files" / f"{PROJECT}_clustered" / "SNV" / f"{PROJECT}_clustered.txt").unlink()
        assert _load(tmp_path, contexts="ID", correction=False) == ({}, {})

    def test_chrom_based_uses_the_chrom_files(self, tmp_path):
        _previous(tmp_path, chrom_based=True)
        assert _load(tmp_path, chrom_based=True) == (REGIONS, IMDS)


class TestSubClassifyNeedsThePreviousResults:
    """The reuse of the saved results is covered by the integration test (test_variant_caller)."""

    def test_missing_previous_results_are_reported(self, tmp_path):
        (tmp_path / "input").mkdir()
        (tmp_path / "input" / "standard.vcf").write_text("")  # the input folder must not be empty
        with pytest.raises(FileNotFoundError, match="imds.pickle"):
            analysis(
                PROJECT, "GRCh37", "96", ["96"], str(tmp_path) + "/",
                analysis="subClassify", subClassify=True, includedVAFs=False,
                plotRainfall=False,
            )
