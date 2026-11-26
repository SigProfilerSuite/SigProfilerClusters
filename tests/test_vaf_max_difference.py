"""The maximum VAF difference between adjacent mutations of an event (default 0.1) can be set:
--vaf_max_difference on the command line, vaf_max_difference in analysis(), and the same
parameter of classifyFunctions.findClustersOfClusters. The default does not change anything.
The matrix generation step needs a genome and is stubbed.
Run with: pytest tests/test_vaf_max_difference.py -v
"""

import inspect
import shutil

import pytest

from clustered_table_helpers import (
    GENOME,
    HEADER_VAF,
    PROJECT,
    SUBCLASSES,
    entries,
    no_matrices,
    row,
    setup,
)
from SigProfilerClusters import SigProfilerClusters as module
from SigProfilerClusters import classifyFunctions, constants, exceptions
from SigProfilerClusters.controllers.cli_controller import CliController, parse_arguments_clusters
from test_cli_dispatch import REQUIRED, recorded  # noqa: F401  (recorded is a fixture)


@pytest.fixture(autouse=True)
def stub_matrix_generation(monkeypatch):
    monkeypatch.setattr(classifyFunctions, "generateMatrices", no_matrices)


# ─── command line ────────────────────────────────────────────


class TestCommandLine:
    def test_default_is_the_constant(self):
        assert parse_arguments_clusters(REQUIRED, "d").vaf_max_difference == constants.VAF_MAX_DIFFERENCE == 0.1

    def test_value_is_parsed_as_a_float(self):
        args = parse_arguments_clusters(REQUIRED + ["--vaf_max_difference", "0.0005"], "d")
        assert args.vaf_max_difference == 0.0005 and isinstance(args.vaf_max_difference, float)

    def test_scientific_notation_is_accepted(self):
        assert parse_arguments_clusters(REQUIRED + ["--vaf_max_difference", "5e-4"], "d").vaf_max_difference == 5e-4

    def test_a_text_is_rejected(self):
        with pytest.raises(SystemExit):
            parse_arguments_clusters(REQUIRED + ["--vaf_max_difference", "small"], "d")

    def test_the_option_is_in_the_help(self, capsys):
        with pytest.raises(SystemExit) as exit_info:
            parse_arguments_clusters(["--help"], "d")
        assert exit_info.value.code == 0
        assert "--vaf_max_difference" in capsys.readouterr().out

    def test_value_reaches_analysis(self, recorded):  # noqa: F811
        CliController().dispatch_sigProfilerClusters(REQUIRED + ["--vaf_max_difference", "0.0005"])
        assert recorded[0]["vaf_max_difference"] == 0.0005

    def test_default_reaches_analysis(self, recorded):  # noqa: F811
        CliController().dispatch_sigProfilerClusters(REQUIRED)
        assert recorded[0]["vaf_max_difference"] == 0.1


# ─── analysis() ──────────────────────────────────────────────

HERE = __import__("os").path.dirname(__import__("os").path.abspath(__file__))
STANDARD_VCF = __import__("os").path.join(HERE, "input_vcf_sim", "standard", "standard.vcf")


class TestAnalysis:
    def test_signature_default(self):
        assert inspect.signature(module.analysis).parameters["vaf_max_difference"].default == 0.1

    @pytest.mark.parametrize("value", [-0.1, "0.1", None, float("nan"), True])
    def test_invalid_values_are_rejected_at_the_start(self, tmp_path, value):
        with pytest.raises(exceptions.InvalidParameterError, match="vaf_max_difference"):
            module.analysis("proj", "GRCh37", "96", ["96"], str(tmp_path) + "/", vaf_max_difference=value)

    @pytest.mark.parametrize("value", [0, 0.0005, 0.1, 1])
    def test_valid_values_pass_the_check(self, tmp_path, value):
        with pytest.raises(Exception) as error:  # it goes on and fails later (no input folder)
            module.analysis("proj", "GRCh37", "96", ["96"], str(tmp_path) + "/", vaf_max_difference=value)
        assert not isinstance(error.value, exceptions.InvalidParameterError)

    @pytest.mark.parametrize("given, expected", [({}, 0.1), ({"vaf_max_difference": 0.0005}, 0.0005)])
    def test_the_value_reaches_the_subclassification(self, tmp_path, monkeypatch, given, expected):
        (tmp_path / "input").mkdir()
        shutil.copy(STANDARD_VCF, tmp_path / "input" / "standard.vcf")
        sims = tmp_path / "output" / "simulations" / "proj_simulations_GRCh37_96_sorted"
        sims.mkdir(parents=True)
        for number in range(100):
            (sims / f"{number}.txt").write_text("")
        seen = {}
        monkeypatch.setattr(module.hotspot, "hotSpotAnalysis", lambda *a, **k: ({}, {}))
        monkeypatch.setattr(module.classifyFunctions, "pullVaf", lambda *a, **k: None)
        monkeypatch.setattr(
            module.classifyFunctions, "findClustersOfClusters", lambda *a, **k: seen.update(k) or seen.update(positional=a)
        )
        monkeypatch.setattr(module.convertToVCF, "generateAllPaths", lambda *a, **k: None)
        monkeypatch.setattr(module.convertToVCF, "convertFiles", lambda *a, **k: None)
        module.analysis(
            "proj", "GRCh37", "96", ["96"], str(tmp_path) + "/", analysis="hotspot", sortSims=False,
            subClassify=True, plotRainfall=False, **given,
        )
        assert seen["vaf_max_difference"] == expected


# ─── findClustersOfClusters ──────────────────────────────────

CHROM_LENGTHS = {GENOME: {"1": 0, "2": 5_000_000}}


def _classify(tmp_path, vafs, **kwargs):
    """Three adjacent mutations (10 bp apart) with these VAFs; returns the subclasses that
    received chromosome 1 mutations."""
    positions = [1000, 1010, 1020]
    rows = [row(p, vaf=v, imd=1 if i == 0 else 10) for i, (p, v) in enumerate(zip(positions, vafs))]
    rows += [row(1000, chrom="2", vaf=0.4), row(1001, chrom="2", vaf=0.4)]
    clustered = setup(tmp_path, rows, HEADER_VAF, f"{PROJECT}_clustered_vaf.txt", {"S1": 100}, False, False)
    classifyFunctions.findClustersOfClusters(
        PROJECT, False, str(tmp_path) + "/", 1_000_000, CHROM_LENGTHS, {"S1": []}, None, GENOME, 1, {}, False, False,
        **kwargs,
    )
    return sorted(
        name for name in SUBCLASSES if [c for c in entries(clustered, name) if c[5] == "1"]
    )


class TestClassification:
    VAFS = [0.0001, 0.0011, 0.0021]  # adjacent differences of 0.001

    def test_with_the_default_the_vaf_difference_is_compatible(self, tmp_path):
        assert "class3" not in _classify(tmp_path, self.VAFS)

    def test_the_default_is_the_same_as_giving_0_1(self, tmp_path):
        (tmp_path / "a").mkdir()
        (tmp_path / "b").mkdir()
        assert _classify(tmp_path / "a", self.VAFS) == _classify(tmp_path / "b", self.VAFS, vaf_max_difference=0.1)

    def test_a_small_cutoff_separates_them(self, tmp_path):
        assert _classify(tmp_path, self.VAFS, vaf_max_difference=0.0005) == ["class3"]

    def test_a_cutoff_between_the_differences_keeps_them_together(self, tmp_path):
        assert "class3" not in _classify(tmp_path, self.VAFS, vaf_max_difference=0.002)
