"""Audit finding 32: the command line exposes the parameters of the custom VAF caller
(--tumor_vaf_column and --vaf_field), which analysis() already had.
Run with: pytest tests/test_cli_vaf_options.py -v
"""

import pytest

from SigProfilerClusters.controllers.cli_controller import CliController, parse_arguments_clusters

from test_cli_dispatch import REQUIRED, recorded  # noqa: F401  (recorded is a fixture)


class TestParsing:
    def test_defaults_are_those_of_analysis(self):
        args = parse_arguments_clusters(REQUIRED, "d")
        assert args.tumor_vaf_column == 10
        assert args.vaf_field == "AF"

    def test_values_are_parsed(self):
        args = parse_arguments_clusters(
            REQUIRED + ["--variant_caller", "custom", "--tumor_vaf_column", "9", "--vaf_field", "VF"],
            "d",
        )
        assert args.tumor_vaf_column == 9
        assert args.vaf_field == "VF"

    def test_the_column_must_be_an_integer(self):
        with pytest.raises(SystemExit):
            parse_arguments_clusters(REQUIRED + ["--tumor_vaf_column", "last"], "d")


class TestDispatch:
    def test_values_reach_analysis(self, recorded):  # noqa: F811
        CliController().dispatch_sigProfilerClusters(
            REQUIRED + ["--variant_caller", "custom", "--tumor_vaf_column", "9", "--vaf_field", "VF"]
        )
        assert recorded[0]["tumor_vaf_column"] == 9
        assert recorded[0]["vaf_field"] == "VF"
        assert recorded[0]["variant_caller"] == "custom"

    def test_defaults_reach_analysis(self, recorded):  # noqa: F811
        CliController().dispatch_sigProfilerClusters(REQUIRED)
        assert recorded[0]["tumor_vaf_column"] == 10
        assert recorded[0]["vaf_field"] == "AF"
