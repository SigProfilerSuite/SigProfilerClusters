"""rainfall in the SBS context when a subclass has a single mutation for a sample and a
chromosome (audit finding 113). pandas returns a Series instead of a DataFrame for a single
matching row, which used to crash the rainfall plots. The subclass files come from the real
classification (findClustersOfClusters) on a small synthetic table; the matrix generation step
needs a genome and is stubbed. No genome reference required.
Run with: pytest tests/test_rainfall_snv.py -v
"""

import pickle

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pytest

from clustered_table_helpers import (
    GENOME,
    HEADER_VAF,
    PROJECT,
    entries,
    no_matrices,
    row,
    setup,
)
from SigProfilerClusters import classifyFunctions
from SigProfilerClusters.classifyFunctions import findClustersOfClusters
from SigProfilerClusters.plottingFunctions import rainfall

CHROMOSOMES = {"1": 400000, "2": 300000, "X": 200000}
CENTROMERES = {GENOME: {"1": [150000, 160000], "2": [100000, 110000], "X": [90000, 100000]}}
CHROM_LENGTHS = {GENOME: {"1": 0, "2": 400000, "X": 700000}}
NON_CLUSTERED_HEADER = "project\tsamples\tID\tgenome\tmutType\tchr\tstart\tend\tref\talt\tmutClass\tIMDplot\tIMD\n"


@pytest.fixture(autouse=True)
def close_figures(monkeypatch):
    monkeypatch.setattr(classifyFunctions, "generateMatrices", no_matrices)
    yield
    plt.close("all")


def _mutation(pos, vaf, chrom="1"):
    return row(pos, vaf=vaf, chrom=chrom)


def _project(tmp_path, rows):
    """Classify `rows` and prepare everything rainfall reads; returns the clustered folder."""
    clustered = setup(tmp_path, rows, HEADER_VAF, f"{PROJECT}_clustered_vaf.txt", None, False, False)
    findClustersOfClusters(
        PROJECT, False, str(tmp_path) + "/", 1_000_000, {GENOME: {"1": 0, "2": 400000, "X": 700000}},
        {"S1": []}, None, GENOME, 1, {}, False, False,
    )
    # three rows: the reader drops the first data row of this table (finding 104), and two rows
    # must remain so that this table is not itself a single-row case
    non_clustered = tmp_path / "output" / "vcf_files" / f"{PROJECT}_nonClustered" / "SNV"
    non_clustered.mkdir(parents=True)
    lines = "".join(
        f"{PROJECT}\tS1\t.\t{GENOME}\tSNP\t1\t{p}\t{p}\tC\tT\tSOMATIC\t5000\t5000\n" for p in (200000, 220000, 240000)
    )
    (non_clustered / f"{PROJECT}_nonClustered.txt").write_text(NON_CLUSTERED_HEADER + lines)
    chrom_dir = tmp_path / "chroms"
    chrom_dir.mkdir()
    for chrom, length in CHROMOSOMES.items():
        (chrom_dir / f"{chrom}.txt").write_text("A" * length)
    return clustered, chrom_dir


def _rainfall(tmp_path, chrom_dir):
    rainfall(
        False, PROJECT, str(tmp_path) + "/", str(chrom_dir) + "/", CHROM_LENGTHS, CENTROMERES,
        "96", True, False, False, 10_000_000, None,
    )
    return tmp_path / "output" / "plots"


# three mutations, the middle one with another VAF: a pair (class IC) and one class III row,
# followed by another event so that the first one is not the last of the table
SINGLE_CLASS_III_ROW = [
    _mutation(1000, 0.3), _mutation(1050, 0.7), _mutation(1100, 0.3),
    _mutation(500000, 0.4), _mutation(500001, 0.4),
]


class TestSingleMutationSubclass:
    def test_the_classification_leaves_one_class_iii_row_for_the_chromosome(self, tmp_path):
        clustered, _ = _project(tmp_path, SINGLE_CLASS_III_ROW)
        rows = [r for r in entries(clustered, "class3") if r[5] == "1"]
        assert [r[6] for r in rows] == ["1050"]

    def test_rainfall_plots_a_subclass_with_a_single_row(self, tmp_path):
        _, chrom_dir = _project(tmp_path, SINGLE_CLASS_III_ROW)
        plots = _rainfall(tmp_path, chrom_dir)
        pdfs = list(plots.glob("*.pdf"))
        assert len(pdfs) == 1 and pdfs[0].stat().st_size > 0

    def test_rainfall_still_plots_subclasses_with_several_rows(self, tmp_path):
        rows = [_mutation(1000 + 50 * i, 0.3 + 0.01 * i) for i in range(4)] + [
            _mutation(500000, 0.4), _mutation(500001, 0.4),
        ]
        _, chrom_dir = _project(tmp_path, rows)
        assert len(list(_rainfall(tmp_path, chrom_dir).glob("*.pdf"))) == 1
