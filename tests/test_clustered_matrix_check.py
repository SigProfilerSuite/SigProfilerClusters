"""Unit tests for hotspot.has_clustered_mutations (audit finding 12): the matrix that decides
whether there are clustered mutations depends on the context; an ID analysis has no "96"
matrix. Matrices are small DataFrames; no genome needed.
Run with: pytest tests/test_clustered_matrix_check.py -v
"""

import pandas as pd
import pytest

from SigProfilerClusters.hotspot import has_clustered_mutations


def _matrix(total):
    return pd.DataFrame({"S1": [total, 0]}, index=["a", "b"])


class TestHasClusteredMutations:
    def test_sbs_with_mutations(self):
        assert has_clustered_mutations({"96": _matrix(3)}, "96") is True

    def test_sbs_without_mutations(self):
        assert has_clustered_mutations({"96": _matrix(0)}, "96") is False

    def test_sbs_without_matrix(self):
        assert has_clustered_mutations({}, "96") is False

    @pytest.mark.parametrize("contexts", ["ID", "INDEL"])
    def test_id_with_mutations_and_no_96_matrix(self, contexts):
        # the case of the finding: it used to be "no clustered mutations"
        assert has_clustered_mutations({"ID": _matrix(3)}, contexts) is True

    @pytest.mark.parametrize("contexts", ["ID", "INDEL"])
    def test_id_without_mutations(self, contexts):
        assert has_clustered_mutations({"ID": _matrix(0)}, contexts) is False

    def test_id_without_matrix(self):
        assert has_clustered_mutations({}, "ID") is False

    def test_id_ignores_the_96_matrix(self):
        assert has_clustered_mutations({"96": _matrix(5)}, "ID") is False

    def test_sbs_ignores_the_id_matrix(self):
        assert has_clustered_mutations({"ID": _matrix(5)}, "96") is False
