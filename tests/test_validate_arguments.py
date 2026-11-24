"""Audit finding 29: the arguments and the input folder are validated at the start of
analysis(), with a clear message, before any folder is created or deleted.
Run with: pytest tests/test_validate_arguments.py -v
"""

import os
import shutil

import pytest

from SigProfilerClusters import exceptions
from SigProfilerClusters.SigProfilerClusters import analysis

HERE = os.path.dirname(os.path.abspath(__file__))
STANDARD_VCF = os.path.join(HERE, "input_vcf_sim", "standard", "standard.vcf")


def _call(path, **kwargs):
    arguments = dict(
        project="proj", genome="GRCh37", contexts="96", simContext=["96"], input_path=str(path) + "/"
    )
    arguments.update(kwargs)
    return analysis(**arguments)


def _raises_invalid(path, word, **kwargs):
    with pytest.raises(exceptions.InvalidParameterError, match=word):
        _call(path, **kwargs)


def _passes_the_check(path, **kwargs):
    """The call goes on (and fails later, there is no input folder), not in the validation."""
    with pytest.raises(Exception) as error:
        _call(path, **kwargs)
    assert not isinstance(error.value, exceptions.InvalidParameterError), str(error.value)


class TestAnalysisValue:
    @pytest.mark.parametrize("value", ["foo", "", None, "ALL", "Hotspot"])
    def test_unsupported_values_are_rejected(self, tmp_path, value):
        _raises_invalid(tmp_path, "analysis", analysis=value)

    def test_the_message_lists_the_valid_values(self, tmp_path):
        with pytest.raises(exceptions.InvalidParameterError) as error:
            _call(tmp_path, analysis="foo")
        for value in ("all", "hotspot", "subClassify"):
            assert value in str(error.value)

    @pytest.mark.parametrize("value", ["all", "hotspot", "subClassify"])
    def test_supported_values_pass(self, tmp_path, value):
        _passes_the_check(tmp_path, analysis=value)


class TestContexts:
    @pytest.mark.parametrize("value", ["DBS", "6144", "SBS288", "", None, 96])
    def test_unsupported_contexts_are_rejected(self, tmp_path, value):
        _raises_invalid(tmp_path, "contexts", contexts=value)

    @pytest.mark.parametrize("value", ["96", "SBS96", "SBS", "ID", "ID83", "83", "INDEL"])
    def test_supported_contexts_pass(self, tmp_path, value):
        _passes_the_check(tmp_path, contexts=value)


class TestSimContext:
    @pytest.mark.parametrize("value", [None, 96, [96], [""]])
    def test_it_must_be_a_list_with_one_name(self, tmp_path, value):
        _raises_invalid(tmp_path, "simContext", simContext=value)

    # these two were already rejected before any input is read; the tests keep it that way
    @pytest.mark.parametrize("value", ["96", ()])
    def test_a_text_or_an_empty_list_is_rejected(self, tmp_path, value):
        _raises_invalid(tmp_path, "simContext", simContext=value)

    def test_several_contexts_are_rejected(self, tmp_path):
        _raises_invalid(tmp_path, "one type of context", simContext=["96", "ID"])

    @pytest.mark.parametrize("value", [["96"], ["ID"], ["6144"], ("96",)])
    def test_one_name_passes(self, tmp_path, value):
        _passes_the_check(tmp_path, simContext=value)


class TestWindowSize:
    @pytest.mark.parametrize("value", [0, -5, 1.5, "1000", None, False])
    def test_it_must_be_a_positive_integer(self, tmp_path, value):
        _raises_invalid(tmp_path, "windowSize", windowSize=value)

    @pytest.mark.parametrize("value", [500000, 1000000, 10000000])
    def test_positive_integers_pass(self, tmp_path, value):
        _passes_the_check(tmp_path, windowSize=value)


class TestVariantCaller:
    def test_unknown_caller_is_rejected_at_the_start(self, tmp_path):
        _raises_invalid(tmp_path, "variant_caller", variant_caller="foo")

    @pytest.mark.parametrize("value", ["standard", "caveman", "mutect2", "custom"])
    def test_known_callers_pass(self, tmp_path, value):
        _passes_the_check(tmp_path, variant_caller=value)


class TestInputFolder:
    def test_hierarchy(self):
        assert issubclass(exceptions.MissingInputError, exceptions.SigProfilerClustersError)
        assert issubclass(exceptions.MissingInputError, FileNotFoundError)
        assert not issubclass(exceptions.MissingInputError, exceptions.InvalidParameterError)

    def test_missing_project_folder(self, tmp_path):
        with pytest.raises(exceptions.MissingInputError, match="does not exist"):
            _call(tmp_path / "nothing")

    def test_project_folder_without_an_input_folder(self, tmp_path):
        with pytest.raises(exceptions.MissingInputError, match="input"):
            _call(tmp_path)

    def test_empty_input_folder(self, tmp_path):
        (tmp_path / "input").mkdir()
        with pytest.raises(exceptions.MissingInputError, match="no input files"):
            _call(tmp_path)

    def test_nothing_is_created_or_deleted_when_the_validation_fails(self, tmp_path):
        (tmp_path / "output" / "keep").mkdir(parents=True)
        (tmp_path / "output" / "keep" / "result.txt").write_text("x")
        with pytest.raises(exceptions.MissingInputError):
            _call(tmp_path)
        assert not (tmp_path / "logs").exists()
        assert (tmp_path / "output" / "keep" / "result.txt").read_text() == "x"

    def test_unsupported_file_format_is_an_error(self, tmp_path):
        (tmp_path / "input").mkdir()
        (tmp_path / "input" / "mutations.xyz").write_text("data")
        with pytest.raises(exceptions.InvalidParameterError, match="format"):
            _call(tmp_path)

    def test_a_valid_input_folder_passes_the_check(self, tmp_path):
        (tmp_path / "input").mkdir()
        shutil.copy(STANDARD_VCF, tmp_path / "input" / "standard.vcf")
        with pytest.raises(exceptions.MissingSimulationsError):  # the next thing that is missing
            _call(tmp_path)
