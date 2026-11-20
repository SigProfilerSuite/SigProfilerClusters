"""Unit tests for the SigProfilerClusters command line layer: argument defaults and
validation, dispatch of parsed arguments to analysis(), and the console entry point.
analysis() is replaced by a recorder, so no genome or simulation is needed.
Run with: pytest tests/test_cli_dispatch.py -v
"""

import argparse
import sys

import pytest

from SigProfilerClusters import sigprofilerclusters_cli
from SigProfilerClusters.controllers import cli_controller
from SigProfilerClusters.controllers.cli_controller import (
    CliController,
    parse_arguments_clusters,
    str2bool,
    str2list,
)

REQUIRED = ["proj", "GRCh37", "96", "96", "/data/proj/"]


@pytest.fixture
def recorded(monkeypatch):
    """Replace SigProfilerClusters.analysis with a function that records its kwargs."""
    calls = []
    monkeypatch.setattr(
        cli_controller.SigProfilerClusters, "analysis", lambda **kw: calls.append(kw)
    )
    return calls


# ─── str2bool / str2list ─────────────────────────────────────


class TestStr2Bool:
    @pytest.mark.parametrize("value", ["YES", "True", "T", "Y", "1"])
    def test_truthy_values_are_case_insensitive(self, value):
        assert str2bool(value) is True

    @pytest.mark.parametrize("value", ["NO", "False", "F", "N", "0"])
    def test_falsy_values_are_case_insensitive(self, value):
        assert str2bool(value) is False

    def test_bool_passthrough(self):
        assert str2bool(True) is True
        assert str2bool(False) is False

    @pytest.mark.parametrize("value", ["maybe", "", "2"])
    def test_invalid_raises(self, value):
        with pytest.raises(argparse.ArgumentTypeError, match="Boolean value expected"):
            str2bool(value)


class TestStr2List:
    def test_single_value(self):
        assert str2list("96") == ["96"]

    def test_comma_separated(self):
        assert str2list("96,DINUC,ID") == ["96", "DINUC", "ID"]

    def test_empty_string(self):
        assert str2list("") == [""]


# ─── parse_arguments_clusters ────────────────────────────────


class TestParseArguments:
    def test_defaults(self):
        args = parse_arguments_clusters(REQUIRED, "d")
        assert args.output_type == "all"
        assert args.analysis == "all"
        assert args.interdistance == "96"
        assert args.exome is False
        assert args.sortSims is True
        assert args.correction is True
        assert args.calculateIMD is True
        assert args.subClassify is False
        assert args.max_cpu is None
        assert args.windowSize == 1000000
        assert args.bedRanges is None
        assert args.variant_caller == "standard"
        assert args.probability is False

    def test_sim_context_is_split(self):
        args = parse_arguments_clusters(["p", "GRCh37", "96", "96,DINUC", "/d/"], "d")
        assert args.simContext == ["96", "DINUC"]

    def test_bare_boolean_flag_means_true(self):
        args = parse_arguments_clusters(REQUIRED + ["--subClassify"], "d")
        assert args.subClassify is True

    def test_boolean_flag_with_value(self):
        args = parse_arguments_clusters(REQUIRED + ["--correction", "False"], "d")
        assert args.correction is False

    def test_integer_options(self):
        args = parse_arguments_clusters(
            REQUIRED + ["--max_cpu", "4", "--windowSize", "500000"], "d"
        )
        assert args.max_cpu == 4
        assert args.windowSize == 500000

    def test_missing_required_argument_exits(self, capsys):
        with pytest.raises(SystemExit) as exc:
            parse_arguments_clusters(["proj", "GRCh37"], "d")
        assert exc.value.code == 2
        assert "required" in capsys.readouterr().err

    def test_invalid_integer_exits(self):
        with pytest.raises(SystemExit):
            parse_arguments_clusters(REQUIRED + ["--max_cpu", "many"], "d")

    def test_invalid_boolean_exits(self):
        with pytest.raises(SystemExit):
            parse_arguments_clusters(REQUIRED + ["--exome", "perhaps"], "d")

    def test_unknown_option_exits(self):
        with pytest.raises(SystemExit):
            parse_arguments_clusters(REQUIRED + ["--not_an_option"], "d")


# ─── CliController.dispatch_sigProfilerClusters ──────────────


class TestDispatch:
    def test_forwards_parsed_values_to_analysis(self, recorded):
        CliController().dispatch_sigProfilerClusters(
            REQUIRED + ["--max_cpu", "2", "--variant_caller", "caveman", "--probability"]
        )
        assert len(recorded) == 1
        kw = recorded[0]
        assert kw["project"] == "proj"
        assert kw["genome"] == "GRCh37"
        assert kw["simContext"] == ["96"]
        assert kw["input_path"] == "/data/proj/"
        assert kw["max_cpu"] == 2
        assert kw["variant_caller"] == "caveman"
        assert kw["probability"] is True

    def test_forwards_every_cli_option(self, recorded):
        CliController().dispatch_sigProfilerClusters(REQUIRED)
        parsed = vars(parse_arguments_clusters(REQUIRED, "d"))
        assert set(recorded[0]) == set(parsed)

    def test_invalid_arguments_do_not_call_analysis(self, recorded):
        with pytest.raises(SystemExit):
            CliController().dispatch_sigProfilerClusters(["only_project"])
        assert recorded == []


# ─── console entry point ─────────────────────────────────────


class TestMainFunction:
    def test_no_command_prints_usage(self, monkeypatch, capsys, recorded):
        monkeypatch.setattr(sys, "argv", ["SigProfilerClusters"])
        sigprofilerclusters_cli.main_function()
        out = capsys.readouterr().out
        assert "Usage: SigProfilerClusters <command>" in out
        assert "analysis: Analyze clustered mutations" in out
        assert recorded == []

    def test_unknown_command_prints_usage(self, monkeypatch, capsys, recorded):
        monkeypatch.setattr(sys, "argv", ["SigProfilerClusters", "bogus"])
        sigprofilerclusters_cli.main_function()
        assert "Commands:" in capsys.readouterr().out
        assert recorded == []

    def test_analysis_command_dispatches(self, monkeypatch, recorded):
        monkeypatch.setattr(sys, "argv", ["SigProfilerClusters", "analysis"] + REQUIRED)
        sigprofilerclusters_cli.main_function()
        assert len(recorded) == 1
        assert recorded[0]["project"] == "proj"

    def test_command_name_is_case_insensitive(self, monkeypatch, recorded):
        monkeypatch.setattr(sys, "argv", ["SigProfilerClusters", "ANALYSIS"] + REQUIRED)
        sigprofilerclusters_cli.main_function()
        assert len(recorded) == 1

    def test_print_usage_lists_all_commands(self, capsys):
        sigprofilerclusters_cli.print_usage({"a": "first", "b": "second"})
        out = capsys.readouterr().out
        assert "  a: first" in out
        assert "  b: second" in out
