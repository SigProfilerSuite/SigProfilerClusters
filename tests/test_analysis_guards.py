"""Unit tests for the argument guards at the start of SigProfilerClusters.analysis.
They return before any genome or simulation is read, so no reference is required.
Run with: pytest tests/test_analysis_guards.py -v
"""

import pytest

from SigProfilerClusters.SigProfilerClusters import analysis


def _analysis(**kwargs):
    return analysis("proj", "GRCh37", "96", ["96"], "/does/not/exist/", **kwargs)


class TestChromBasedGuard:
    def test_chrom_based_exits(self):
        with pytest.raises(SystemExit):
            _analysis(chrom_based=True)

    def test_chrom_based_explains_that_it_is_deprecated(self, capsys):
        with pytest.raises(SystemExit):
            _analysis(chrom_based=True)
        out = capsys.readouterr().out
        assert "chrom_based has been deprecated" in out
        assert "rerun with chrom_based=False" in out


class TestOutputTypeGuard:
    """audit finding 8: only "all" (or None) is supported; the others used to fail later"""

    @pytest.mark.parametrize("value", ["original", "simulations", "foo", ""])
    def test_unsupported_values_are_rejected_at_the_start(self, value):
        with pytest.raises(ValueError, match="output_type"):
            _analysis(output_type=value)

    def test_the_message_names_the_value_and_the_valid_ones(self):
        with pytest.raises(ValueError) as error:
            _analysis(output_type="original")
        assert "'original'" in str(error.value)
        assert "'all'" in str(error.value)

    @pytest.mark.parametrize("value", ["all", None])
    def test_supported_values_pass_the_guard(self, value):
        # they go on and fail later, reading the missing input folder, not in the guard
        with pytest.raises(Exception) as error:
            _analysis(output_type=value)
        assert "output_type" not in str(error.value)

    def test_cli_accepts_only_all(self):
        from SigProfilerClusters.controllers.cli_controller import parse_arguments_clusters

        base = ["p", "GRCh37", "96", "96", "/x/"]
        assert parse_arguments_clusters(base + ["--output_type", "all"], "d").output_type == "all"
        with pytest.raises(SystemExit):
            parse_arguments_clusters(base + ["--output_type", "original"], "d")


class TestStderrIsRestored:
    """audit finding 13: analysis() redirects sys.stderr to its .err file; it must give the
    original back, also when it stops with an error, and never leave it closed"""

    def test_stderr_restored_after_an_error(self, tmp_path):
        import sys

        original = sys.stderr
        with pytest.raises(FileNotFoundError):  # there is no input/ folder: fails after the redirection
            analysis("proj", "GRCh37", "96", ["96"], str(tmp_path) + "/")
        assert sys.stderr is original
        assert not sys.stderr.closed

    def test_stderr_restored_after_a_system_exit(self, tmp_path):
        import sys

        original = sys.stderr
        with pytest.raises(SystemExit):
            analysis("proj", "GRCh37", "96", ["96"], str(tmp_path) + "/", chrom_based=True)
        assert sys.stderr is original

    def test_the_err_file_is_still_created(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            analysis("proj", "GRCh37", "96", ["96"], str(tmp_path) + "/")
        assert len(list((tmp_path / "logs").glob("*.err"))) == 1

    def test_signature_is_unchanged(self):
        import inspect

        parameters = list(inspect.signature(analysis).parameters)
        assert parameters[:5] == ["project", "genome", "contexts", "simContext", "input_path"]
        assert "tumor_vaf_column" in parameters
