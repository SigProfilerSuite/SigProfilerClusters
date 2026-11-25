"""Audit finding 19: analysis() reports its errors with exceptions of the package, not with
sys.exit(), so that a caller can catch them and the command line can exit with a code.
Run with: pytest tests/test_exceptions.py -v
"""

import os
import shutil

import pytest

from SigProfilerClusters import exceptions
from SigProfilerClusters.SigProfilerClusters import analysis

HERE = os.path.dirname(os.path.abspath(__file__))
STANDARD_VCF = os.path.join(HERE, "input_vcf_sim", "standard", "standard.vcf")


def _project_with_input(tmp_path):
    """A project folder with one VCF and no simulations."""
    (tmp_path / "input").mkdir()
    shutil.copy(STANDARD_VCF, tmp_path / "input" / "standard.vcf")
    return str(tmp_path) + "/"


def _analysis(path, **kwargs):
    return analysis("standard", "GRCh37", "96", ["96"], path, **kwargs)


class TestHierarchy:
    def test_all_errors_share_the_base_class(self):
        assert issubclass(exceptions.InvalidParameterError, exceptions.SigProfilerClustersError)
        assert issubclass(exceptions.MissingSimulationsError, exceptions.SigProfilerClustersError)

    def test_invalid_parameter_is_still_a_value_error(self):
        # callers that already catch ValueError (the output_type guard) keep working
        assert issubclass(exceptions.InvalidParameterError, ValueError)

    def test_base_class_is_an_ordinary_exception(self):
        # not SystemExit: "except Exception" catches it
        assert issubclass(exceptions.SigProfilerClustersError, Exception)


class TestAnalysisRaises:
    def test_chrom_based_raises_invalid_parameter(self, tmp_path):
        with pytest.raises(exceptions.InvalidParameterError, match="chrom_based"):
            _analysis(str(tmp_path) + "/", chrom_based=True)

    def test_output_type_raises_invalid_parameter(self, tmp_path):
        with pytest.raises(exceptions.InvalidParameterError, match="output_type"):
            _analysis(str(tmp_path) + "/", output_type="original")

    def test_several_simulated_contexts_raise_invalid_parameter(self, tmp_path):
        path = _project_with_input(tmp_path)
        with pytest.raises(exceptions.InvalidParameterError, match="one type of context"):
            analysis("standard", "GRCh37", "96", ["96", "ID"], path)

    def test_missing_simulations_raise_a_clear_error(self, tmp_path):
        path = _project_with_input(tmp_path)
        with pytest.raises(exceptions.MissingSimulationsError, match="no simulated data"):
            _analysis(path)

    def test_too_few_simulations_raise_a_clear_error(self, tmp_path):
        path = _project_with_input(tmp_path)
        sims = tmp_path / "output" / "simulations" / "standard_simulations_GRCh37_96"
        sims.mkdir(parents=True)
        (sims / "1.txt").write_text("")
        with pytest.raises(exceptions.MissingSimulationsError, match="minimum of 100"):
            _analysis(path)

    def test_nothing_in_analysis_calls_sys_exit(self):
        import inspect

        from SigProfilerClusters import SigProfilerClusters as module

        assert "sys.exit" not in inspect.getsource(module)


class TestCommandLine:
    def test_error_is_reported_and_the_exit_code_is_not_zero(self, tmp_path, capsys, monkeypatch):
        import sys

        from SigProfilerClusters import sigprofilerclusters_cli

        monkeypatch.setattr(
            sys, "argv",
            ["SigProfilerClusters", "analysis", "p", "GRCh37", "96", "96", str(tmp_path) + "/",
             "--output_type", "all"],
        )

        def failing_analysis(**kwargs):
            raise exceptions.MissingSimulationsError("there are no simulated data")

        monkeypatch.setattr(
            "SigProfilerClusters.SigProfilerClusters.analysis", failing_analysis
        )
        with pytest.raises(SystemExit) as exit_info:
            sigprofilerclusters_cli.main_function()
        assert exit_info.value.code == 1
        assert "there are no simulated data" in capsys.readouterr().err


class TestNoClusteredMutations:
    """Finishing without clustered mutations is not an error: analysis() returns normally"""

    def test_analysis_returns_before_the_subclassification(self, tmp_path, monkeypatch, capsys):
        path = _project_with_input(tmp_path)
        sims = tmp_path / "output" / "simulations" / "standard_simulations_GRCh37_96_sorted"
        sims.mkdir(parents=True)
        for number in range(100):
            (sims / (str(number) + ".txt")).write_text("")

        from SigProfilerClusters import SigProfilerClusters as module

        monkeypatch.setattr(module.hotspot, "hotSpotAnalysis", lambda *a, **k: (None, None))

        def must_not_run(*a, **k):
            raise AssertionError("the subclassification must not run")

        monkeypatch.setattr(module.classifyFunctions, "pullVaf", must_not_run)
        monkeypatch.setattr(module.classifyFunctions, "findClustersOfClusters", must_not_run)

        result = _analysis(path, analysis="hotspot", sortSims=False, subClassify=True)
        assert result is None
        assert "No clustered mutations were found" in capsys.readouterr().out


class TestUnsupportedGenome:
    """audit finding 23: a genome without centromere data fails at the start, with a clear
    error, not with a KeyError after the distances were calculated"""

    @pytest.mark.parametrize("genome", ["hg19", "GRCh36", "", "grch37"])
    def test_unsupported_genome_is_rejected_at_the_start(self, tmp_path, genome):
        with pytest.raises(exceptions.InvalidParameterError) as error:
            analysis("proj", genome, "96", ["96"], str(tmp_path) + "/")
        assert repr(genome) in str(error.value)

    def test_the_message_lists_the_supported_genomes(self, tmp_path):
        with pytest.raises(exceptions.InvalidParameterError) as error:
            analysis("proj", "hg19", "96", ["96"], str(tmp_path) + "/")
        for genome in ("GRCh37", "GRCh38", "mm9", "mm10", "mm39", "rn6"):
            assert genome in str(error.value)

    @pytest.mark.parametrize("genome", ["GRCh37", "GRCh38", "mm9", "mm10", "mm39", "rn6"])
    def test_supported_genomes_pass_the_check(self, tmp_path, genome):
        # they go on and fail later (no input folder), not because of the genome
        with pytest.raises(Exception) as error:
            analysis("proj", genome, "96", ["96"], str(tmp_path) + "/")
        assert not isinstance(error.value, exceptions.InvalidParameterError)
