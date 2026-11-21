"""Unit tests for plottingFunctions.rainfall in the ID context (clustered vs non-clustered
indels along the genome, one page per sample).
Inputs are synthetic clustered and non-clustered indel tables, IMD pickles and tiny
chromosome files built with tmp_path. No genome reference required.
Run with: pytest tests/test_rainfall.py -v
"""

import pickle
import re

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pytest

from SigProfilerClusters.plottingFunctions import rainfall

PROJECT = "proj"
HEADER = "project\tsamples\tID\tgenome\tmutType\tchr\tstart\tend\tref\talt\tmutClass\tIMDplot\tIMD\n"
CHROMOSOMES = {"1": 400000, "2": 300000, "X": 200000}
CENTROMERES = {"GRCh37": {"1": [150000, 160000], "2": [100000, 110000], "X": [90000, 100000]}}
CHROM_LENGTHS = {"GRCh37": {"1": 0, "2": 400000, "X": 700000}}


@pytest.fixture(autouse=True)
def close_figures():
    yield
    plt.close("all")


def _row(sample, pos, imd, chrom="1"):
    return f"{PROJECT}\t{sample}\t.\tGRCh37\tINDEL\t{chrom}\t{pos}\t{pos}\tA\tAT\tSOMATIC\t{imd}\t{imd}\n"


def _setup(tmp_path, samples=("S1",), correction=True):
    chrom_dir = tmp_path / "chroms"
    chrom_dir.mkdir()
    for chrom, length in CHROMOSOMES.items():
        (chrom_dir / f"{chrom}.txt").write_text("A" * length)

    data = tmp_path / "output" / "simulations" / "data"
    data.mkdir(parents=True)
    with open(str(data / "imds.pickle"), "wb") as f:
        pickle.dump({sample: 1000 for sample in samples}, f)
    if correction:
        with open(str(data / "imds_corrected.pickle"), "wb") as f:
            pickle.dump({sample: {} for sample in samples}, f)

    suffix = "_corrected" if correction else ""
    base = tmp_path / "output" / f"vcf_files{suffix}"
    clustered = base / f"{PROJECT}_clustered" / "INDEL"
    non_clustered = base / f"{PROJECT}_nonClustered" / "INDEL"
    clustered.mkdir(parents=True)
    non_clustered.mkdir(parents=True)
    clustered_rows = "".join(_row(s, 1000 + i * 50, 100) for s in samples for i in range(20))
    non_clustered_rows = "".join(_row(s, 50000 + i * 3000, 3000) for s in samples for i in range(30))
    (clustered / f"{PROJECT}_clustered.txt").write_text(HEADER + clustered_rows)
    (non_clustered / f"{PROJECT}_nonClustered.txt").write_text(HEADER + non_clustered_rows)
    return chrom_dir


def _run(tmp_path, correction=True, samples=("S1",)):
    chrom_dir = _setup(tmp_path, samples, correction)
    rainfall(
        False,
        PROJECT,
        str(tmp_path) + "/",
        str(chrom_dir) + "/",
        CHROM_LENGTHS,
        CENTROMERES,
        "ID",
        False,
        False,
        correction,
        10000000,
        None,
    )
    return tmp_path / "output" / "plots"


def _pages(pdf):
    return len(re.findall(rb"/Type\s*/Page[^s]", pdf.read_bytes()))


class TestRainfallId:
    def test_a_pdf_is_written_in_the_plots_folder(self, tmp_path):
        plots = _run(tmp_path)
        pdfs = list(plots.glob("*.pdf"))
        assert [p.name for p in pdfs] == ["rainfallPlots_clustered_proj_corrected.pdf"]
        assert pdfs[0].read_bytes().startswith(b"%PDF")

    def test_uncorrected_run_has_no_corrected_suffix(self, tmp_path):
        plots = _run(tmp_path, correction=False)
        assert [p.name for p in plots.glob("*.pdf")] == ["rainfallPlots_clustered_proj.pdf"]

    def test_one_page_per_sample(self, tmp_path):
        plots = _run(tmp_path, samples=("S1", "S2", "S3"))
        assert _pages(next(plots.glob("*.pdf"))) == 3

    def test_existing_plots_folder_is_reused(self, tmp_path):
        (tmp_path / "output" / "plots").mkdir(parents=True)
        plots = _run(tmp_path)
        assert len(list(plots.glob("*.pdf"))) == 1

    def test_missing_imd_pickle_raises(self, tmp_path):
        chrom_dir = _setup(tmp_path)
        (tmp_path / "output" / "simulations" / "data" / "imds.pickle").unlink()
        with pytest.raises(FileNotFoundError):
            rainfall(
                False, PROJECT, str(tmp_path) + "/", str(chrom_dir) + "/", CHROM_LENGTHS,
                CENTROMERES, "ID", False, False, True, 10000000, None,
            )

    def test_missing_mutation_table_raises(self, tmp_path):
        chrom_dir = _setup(tmp_path)
        table = tmp_path / "output" / "vcf_files_corrected" / f"{PROJECT}_clustered" / "INDEL" / f"{PROJECT}_clustered.txt"
        table.unlink()
        with pytest.raises(FileNotFoundError):
            rainfall(
                False, PROJECT, str(tmp_path) + "/", str(chrom_dir) + "/", CHROM_LENGTHS,
                CENTROMERES, "ID", False, False, True, 10000000, None,
            )
