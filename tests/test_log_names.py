"""Audit finding 31 (part 1): the log files of two runs of the same project on the same day
do not overwrite each other: their names carry the time as well as the date.
Run with: pytest tests/test_log_names.py -v
"""

import os
import re
import shutil

import pytest

from SigProfilerClusters import SigProfilerClusters as module
from SigProfilerClusters.exceptions import MissingSimulationsError

HERE = os.path.dirname(os.path.abspath(__file__))
STANDARD_VCF = os.path.join(HERE, "input_vcf_sim", "standard", "standard.vcf")


def _run(tmp_path):
    """A run that stops after the log files are created (there are no simulations)"""
    with pytest.raises(MissingSimulationsError):
        module.analysis("proj", "GRCh37", "96", ["96"], str(tmp_path) + "/")


def _project(tmp_path):
    (tmp_path / "input").mkdir()
    shutil.copy(STANDARD_VCF, tmp_path / "input" / "standard.vcf")


def test_log_stamp_has_the_date_and_the_time():
    assert re.fullmatch(r"\d{4}-\d{2}-\d{2}_\d{6}", module.log_stamp())


def test_log_names_carry_the_time(tmp_path):
    _project(tmp_path)
    _run(tmp_path)
    names = sorted(p.name for p in (tmp_path / "logs").iterdir())
    assert len(names) == 2
    for name in names:
        assert re.fullmatch(r"SigProfilerClusters_proj_GRCh37_\d{4}-\d{2}-\d{2}_\d{6}\.(out|err)", name), name


def test_two_runs_at_different_times_keep_both_logs(tmp_path, monkeypatch):
    _project(tmp_path)
    for stamp in ("2026-11-25_100000", "2026-11-25_110000"):
        monkeypatch.setattr(module, "log_stamp", lambda stamp=stamp: stamp)
        _run(tmp_path)
    names = sorted(p.name for p in (tmp_path / "logs").iterdir())
    assert len(names) == 4
    assert any("100000" in n for n in names) and any("110000" in n for n in names)


def test_the_same_stamp_still_replaces_the_log_of_that_run(tmp_path, monkeypatch):
    _project(tmp_path)
    monkeypatch.setattr(module, "log_stamp", lambda: "2026-11-25_100000")
    _run(tmp_path)
    _run(tmp_path)
    assert len(list((tmp_path / "logs").iterdir())) == 2
