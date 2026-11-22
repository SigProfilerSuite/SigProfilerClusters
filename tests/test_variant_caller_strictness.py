"""The integration test must fail when the pipeline fails or its output differs from the
reference (audit finding #16). SigProfilerSimulator and SigProfilerClusters.analysis are replaced
by fakes that write a chosen output, so no genome is needed.
Run with: pytest tests/test_variant_caller_strictness.py -v
"""

import random

import numpy as np
import pytest

import test_variant_caller as tvc

CALLER = "standard"
REFERENCE = "header\nrow1\nrow2\n"


@pytest.fixture
def layout(tmp_path):
    (tmp_path / "input_vcf_sim" / CALLER).mkdir(parents=True)
    (tmp_path / "output_ref_sim" / CALLER).mkdir(parents=True)
    (tmp_path / "output_ref_sim" / CALLER / f"{CALLER}_clustered_vaf.txt").write_text(REFERENCE)
    return tmp_path


def _fake_pipeline(monkeypatch, base, produced=REFERENCE, simulator_error=None, analysis_error=None):
    def simulator(**kwargs):
        if simulator_error:
            raise simulator_error

    def analysis(**kwargs):
        if analysis_error:
            raise analysis_error
        out = base / "input_vcf_sim" / CALLER / "output" / "vcf_files_corrected" / f"{CALLER}_clustered"
        out.mkdir(parents=True)
        (out / f"{CALLER}_clustered_vaf.txt").write_text(produced)

    monkeypatch.setattr(tvc.sigSim, "SigProfilerSimulator", simulator)
    monkeypatch.setattr(tvc.hp, "analysis", analysis)


def _run(base):
    tvc.run_analysis_for_variant_callers_separate_steps([CALLER], str(base), "GRCh37", 1)


class TestStrictness:
    def test_matching_output_passes(self, layout, monkeypatch):
        _fake_pipeline(monkeypatch, layout)
        _run(layout)

    def test_output_that_differs_from_the_reference_fails(self, layout, monkeypatch):
        _fake_pipeline(monkeypatch, layout, produced="header\nrow1\nCHANGED\n")
        with pytest.raises(AssertionError, match="differ"):
            _run(layout)

    def test_output_with_another_line_count_fails(self, layout, monkeypatch):
        _fake_pipeline(monkeypatch, layout, produced="header\nrow1\n")
        with pytest.raises(AssertionError):
            _run(layout)

    def test_simulator_exception_is_not_swallowed(self, layout, monkeypatch):
        _fake_pipeline(monkeypatch, layout, simulator_error=RuntimeError("simulator broke"))
        with pytest.raises(RuntimeError, match="simulator broke"):
            _run(layout)

    def test_analysis_exception_is_not_swallowed(self, layout, monkeypatch):
        _fake_pipeline(monkeypatch, layout, analysis_error=RuntimeError("analysis broke"))
        with pytest.raises(RuntimeError, match="analysis broke"):
            _run(layout)

    def test_missing_output_fails(self, layout, monkeypatch):
        monkeypatch.setattr(tvc.sigSim, "SigProfilerSimulator", lambda **kwargs: None)
        monkeypatch.setattr(tvc.hp, "analysis", lambda **kwargs: None)
        with pytest.raises(AssertionError, match="No output"):
            _run(layout)

    def test_missing_input_directory_fails(self, layout, monkeypatch):
        _fake_pipeline(monkeypatch, layout)
        with pytest.raises(AssertionError, match="Input path does not exist"):
            tvc.run_analysis_for_variant_callers_separate_steps(["other"], str(layout), "GRCh37", 1)

    def test_one_genome_test_propagates_failures(self, monkeypatch):
        def failing(*args, **kwargs):
            raise AssertionError("output differs")

        monkeypatch.setattr(tvc, "run_analysis_for_variant_callers_separate_steps", failing)
        with pytest.raises(AssertionError, match="output differs"):
            tvc.run_one_genome_test(CALLER)


class TestSimulationSeed:
    """The simulator draws the distribution of mutations over chromosomes from the global
    generators; they are seeded before every simulator call (audit register, finding 114)."""

    def _draws_seen_by_the_simulator(self, layout, monkeypatch, state):
        seen = []

        def simulator(**kwargs):
            seen.append((random.random(), float(np.random.random())))

        monkeypatch.setattr(tvc.sigSim, "SigProfilerSimulator", simulator)
        monkeypatch.setattr(tvc.hp, "analysis", lambda **kwargs: None)
        random.seed(state)
        np.random.seed(state)
        with pytest.raises(AssertionError):  # no output is generated: only the simulator call matters
            _run(layout)
        return seen[0]

    def test_the_simulator_sees_the_same_draws_whatever_the_previous_state(self, layout, monkeypatch):
        first = self._draws_seen_by_the_simulator(layout, monkeypatch, 1)
        second = self._draws_seen_by_the_simulator(layout, monkeypatch, 2)
        assert first == second

    def test_the_draws_come_from_the_documented_seed(self, layout, monkeypatch):
        seen = self._draws_seen_by_the_simulator(layout, monkeypatch, 5)
        random.seed(tvc.SIMULATION_SEED)
        np.random.seed(tvc.SIMULATION_SEED)
        assert seen == (random.random(), float(np.random.random()))
