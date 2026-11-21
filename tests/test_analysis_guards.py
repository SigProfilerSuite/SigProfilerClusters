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
