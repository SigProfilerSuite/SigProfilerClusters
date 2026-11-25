"""Audit findings 118 and 119: the names of the contexts are standardized. SBS96 and ID83 are
the standard names; 96/SBS and ID/83/INDEL are accepted as aliases (DBS is not an input:
it is an output of SPC). The context is normalized once at the start of analysis(), and the
places that ask "is this an indel analysis?" accept every alias.
Run with: pytest tests/test_context_names.py -v
"""

import os
import pickle
import re
import shutil
from pathlib import Path

import pytest

from SigProfilerClusters import constants, exceptions, hotspot
from SigProfilerClusters import SigProfilerClusters as module
from SigProfilerClusters.SigProfilerClusters import analysis, load_previous_analysis

SBS_ALIASES = ["96", "SBS96", "SBS"]
ID_ALIASES = ["ID", "ID83", "83", "INDEL"]
HERE = os.path.dirname(os.path.abspath(__file__))
STANDARD_VCF = os.path.join(HERE, "input_vcf_sim", "standard", "standard.vcf")
PACKAGE = Path(HERE).parent / "SigProfilerClusters"


def _call(path, **kwargs):
    arguments = dict(
        project="proj", genome="GRCh37", contexts="96", simContext=["96"], input_path=str(path) + "/"
    )
    arguments.update(kwargs)
    return analysis(**arguments)


class TestCanonicalContext:
    @pytest.mark.parametrize("value", SBS_ALIASES)
    def test_sbs_aliases(self, value):
        assert constants.canonical_context(value) == "SBS96"

    @pytest.mark.parametrize("value", ID_ALIASES)
    def test_indel_aliases(self, value):
        assert constants.canonical_context(value) == "ID83"

    @pytest.mark.parametrize("value", ["DBS", "6144", "SBS288", "ID415", "sbs96", "", None, 96, ["96"]])
    def test_anything_else_is_not_a_context(self, value):
        assert constants.canonical_context(value) is None

    @pytest.mark.parametrize("value", ID_ALIASES)
    def test_is_indel_context(self, value):
        assert constants.is_indel_context(value) and not constants.is_sbs_context(value)

    @pytest.mark.parametrize("value", SBS_ALIASES)
    def test_is_sbs_context(self, value):
        assert constants.is_sbs_context(value) and not constants.is_indel_context(value)

    def test_unknown_value_is_neither(self):
        assert not constants.is_indel_context("DBS") and not constants.is_sbs_context("DBS")


class TestAnalysisAcceptsTheAliases:
    @pytest.mark.parametrize("value", SBS_ALIASES + ID_ALIASES)
    def test_every_alias_passes_the_validation(self, tmp_path, value):
        with pytest.raises(exceptions.MissingInputError):  # the next thing that is missing
            _call(tmp_path, contexts=value)

    def test_dbs_is_rejected_and_explained(self, tmp_path):
        with pytest.raises(exceptions.InvalidParameterError, match="DBS") as error:
            _call(tmp_path, contexts="DBS")
        assert "output" in str(error.value)

    def test_the_message_lists_the_standard_names_and_the_aliases(self, tmp_path):
        with pytest.raises(exceptions.InvalidParameterError) as error:
            _call(tmp_path, contexts="foo")
        for name in ("SBS96", "ID83", "SBS", "INDEL"):
            assert name in str(error.value)


class TestInterdistance:
    @pytest.mark.parametrize("value", SBS_ALIASES + ID_ALIASES + ["SI"])
    def test_aliases_and_si_pass_the_validation(self, tmp_path, value):
        with pytest.raises(exceptions.MissingInputError):
            _call(tmp_path, interdistance=value)

    @pytest.mark.parametrize("value", ["foo", "", None, "DBS"])
    def test_other_values_are_rejected(self, tmp_path, value):
        with pytest.raises(exceptions.InvalidParameterError, match="interdistance"):
            _call(tmp_path, interdistance=value)


def _project_for_hotspot(tmp_path, sim_context="96"):
    """A project that goes through analysis="hotspot" up to the call of hotSpotAnalysis"""
    (tmp_path / "input").mkdir()
    shutil.copy(STANDARD_VCF, tmp_path / "input" / "standard.vcf")
    sims = tmp_path / "output" / "simulations" / f"proj_simulations_GRCh37_{sim_context}_sorted"
    sims.mkdir(parents=True)
    for number in range(100):
        (sims / f"{number}.txt").write_text("")
    return str(tmp_path) + "/"


class TestNormalizedAtTheStart:
    @pytest.mark.parametrize(
        "value, expected",
        [("96", "SBS96"), ("SBS", "SBS96"), ("SBS96", "SBS96"), ("ID", "ID83"), ("83", "ID83"), ("INDEL", "ID83")],
    )
    def test_the_analysis_receives_the_standard_name(self, tmp_path, monkeypatch, value, expected):
        sim_context = "ID" if expected == "ID83" else "96"
        path = _project_for_hotspot(tmp_path, sim_context)
        received = {}

        def stub(project, genome, contexts, *args, **kwargs):
            received["contexts"] = contexts
            return None, None

        monkeypatch.setattr(module.hotspot, "hotSpotAnalysis", stub)
        _call(path, contexts=value, simContext=[sim_context], analysis="hotspot", sortSims=False)
        assert received["contexts"] == expected


class TestLegacyResults:
    """Results of an earlier version have '96'/'ID' in the names of the distance files"""

    def test_hotspot_on_results_with_the_old_names_explains_what_to_do(self, tmp_path):
        path = _project_for_hotspot(tmp_path)
        old = tmp_path / "output" / "simulations" / "proj_intradistance_original_GRCh37_96"
        old.mkdir()
        with pytest.raises(exceptions.MissingInputError, match="earlier version"):
            _call(path, analysis="hotspot", sortSims=False)

    def test_hotspot_with_the_new_names_goes_on(self, tmp_path, monkeypatch):
        path = _project_for_hotspot(tmp_path)
        new = tmp_path / "output" / "simulations" / "proj_intradistance_original_GRCh37_SBS96"
        new.mkdir()
        monkeypatch.setattr(module.hotspot, "hotSpotAnalysis", lambda *a, **k: (None, None))
        _call(path, analysis="hotspot", sortSims=False)


class TestEverySiteAcceptsEveryAlias:
    @pytest.mark.parametrize("value", ID_ALIASES)
    def test_load_previous_analysis_does_not_ask_for_the_snv_table(self, tmp_path, value):
        data = tmp_path / "output" / "simulations" / "data"
        data.mkdir(parents=True)
        (data / "imds.pickle").write_bytes(pickle.dumps({"S1": 100}))
        # nothing else is needed for an indel analysis (there is no SNV table)
        assert load_previous_analysis("proj", str(tmp_path) + "/", value, False, False) == ({}, {})

    @pytest.mark.parametrize("value", SBS_ALIASES)
    def test_load_previous_analysis_asks_for_the_snv_table_with_sbs(self, tmp_path, value):
        data = tmp_path / "output" / "simulations" / "data"
        data.mkdir(parents=True)
        (data / "imds.pickle").write_bytes(b"x")
        with pytest.raises(FileNotFoundError, match="SNV"):
            load_previous_analysis("proj", str(tmp_path) + "/", value, False, False)

    @pytest.mark.parametrize("value", ID_ALIASES)
    def test_convert_to_vcf_paths_for_indels(self, tmp_path, value):
        from SigProfilerClusters.convertToVCF import generateAllPaths

        generateAllPaths(str(tmp_path), value)
        assert (tmp_path / "output" / "clustered" / "ID").is_dir()
        assert not (tmp_path / "output" / "clustered" / "kataegis").exists()

    @pytest.mark.parametrize("value", SBS_ALIASES)
    def test_convert_to_vcf_paths_for_sbs(self, tmp_path, value):
        from SigProfilerClusters.convertToVCF import generateAllPaths

        generateAllPaths(str(tmp_path), value)
        assert (tmp_path / "output" / "clustered" / "kataegis").is_dir()

    @pytest.mark.parametrize("value", ID_ALIASES)
    def test_clustered_mutations_check_uses_the_id_matrix(self, value):
        import pandas as pd

        matrices = {"ID": pd.DataFrame({"S": [1]}), "96": pd.DataFrame({"S": [0]})}
        assert hotspot.has_clustered_mutations(matrices, value) is True

    def test_no_module_compares_the_context_with_a_literal_name(self):
        pattern = re.compile(r"""\b(?:file_)?contexts?\s*(?:[!=]=|in)\s*\(?\s*["'](?:96|ID|INDEL)["']""")
        offenders = []
        for source in PACKAGE.rglob("*.py"):
            if source.name == "constants.py":
                continue
            for number, line in enumerate(source.read_text().splitlines(), 1):
                if pattern.search(line) and not line.lstrip().startswith("#"):
                    offenders.append(f"{source.name}:{number}: {line.strip()}")
        assert not offenders, "\n".join(offenders)
