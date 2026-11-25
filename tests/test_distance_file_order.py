"""Audit finding 121 (and 4): the file of distances of an original sample is written by several
workers, each one appending the chromosomes it was given, so the order of its rows depended on
which worker finished first and the order of every table built from it (and the numbers of the
groups) changed from one run to the next. The rows are now sorted by chromosome and position.
The tests use synthetic mutations only.
Run with: pytest tests/test_distance_file_order.py -v
"""

import pytest

import test_distance_functions as tdf
from SigProfilerClusters import SigProfilerClusters as module
from SigProfilerClusters.SigProfilerClusters import distance_one_file

CHROMOSOMES = ["1", "2", "3", "10", "X"]
CENTROMERES = {tdf.GENOME: {chrom: [5000, 6000] for chrom in CHROMOSOMES}}


def _chromosome_files(tmp_path):
    source = tmp_path / "original"
    source.mkdir()
    for chrom in CHROMOSOMES:
        rows = [tdf._mutation(pos, chrom=chrom) for pos in (100, 110, 7000, 7100)]
        (source / f"{chrom}.txt").write_text("".join("\t".join(r) + "\n" for r in rows))
    return source


def _write(tmp_path, name, order):
    """The distance file as written by workers that finish the chromosomes in this order"""
    source = _chromosome_files(tmp_path) if not (tmp_path / "original").exists() else tmp_path / "original"
    out = tmp_path / name
    out.mkdir()
    for chrom in order:  # one call per chromosome: a worker that got only that file
        distance_one_file(
            [f"{chrom}.txt"], str(source) + "/", str(out) + "/", tdf.CONTEXT, tdf.GENOME, CENTROMERES
        )
    return out


def _read(directory):
    return {p.name: p.read_text() for p in sorted(directory.iterdir())}


class TestTheOrderDependsOnTheWorkers:
    def test_without_sorting_another_finishing_order_gives_another_file(self, tmp_path):
        a = _write(tmp_path, "a", CHROMOSOMES)
        b = _write(tmp_path, "b", list(reversed(CHROMOSOMES)))
        assert _read(a) != _read(b)  # the cause of the problem
        assert sorted(_read(a)["S1_96_intradistance.txt"].splitlines()) == sorted(
            _read(b)["S1_96_intradistance.txt"].splitlines()
        )  # the same rows


class TestSortedFiles:
    @pytest.mark.parametrize(
        "order",
        [CHROMOSOMES, list(reversed(CHROMOSOMES)), ["10", "X", "2", "1", "3"], ["3", "1", "X", "10", "2"]],
    )
    def test_the_file_does_not_depend_on_the_finishing_order(self, tmp_path, order):
        reference = _write(tmp_path, "reference", CHROMOSOMES)
        other = _write(tmp_path, "other", order)
        module.sort_distance_files(str(reference) + "/")
        module.sort_distance_files(str(other) + "/")
        assert _read(reference) == _read(other)

    def test_rows_are_in_natural_chromosome_order_then_position(self, tmp_path):
        out = _write(tmp_path, "out", ["X", "10", "2", "1", "3"])
        module.sort_distance_files(str(out) + "/")
        rows = [line.split("\t") for line in (out / "S1_96_intradistance.txt").read_text().splitlines()]
        assert [r[2] for r in rows if r[3] == "100"] == ["1", "2", "3", "10", "X"]
        for chrom in CHROMOSOMES:
            positions = [int(r[3]) for r in rows if r[2] == chrom]
            assert positions == sorted(positions)

    def test_the_content_is_not_changed(self, tmp_path):
        out = _write(tmp_path, "out", ["X", "10", "2", "1", "3"])
        before = sorted((out / "S1_96_intradistance.txt").read_text().splitlines())
        module.sort_distance_files(str(out) + "/")
        assert sorted((out / "S1_96_intradistance.txt").read_text().splitlines()) == before

    def test_other_files_of_the_folder_are_left_alone(self, tmp_path):
        out = _write(tmp_path, "out", CHROMOSOMES)
        (out / "notes.txt").write_text("b\na\n")
        module.sort_distance_files(str(out) + "/")
        assert (out / "notes.txt").read_text() == "b\na\n"

    def test_an_empty_folder_is_fine(self, tmp_path):
        (tmp_path / "empty").mkdir()
        module.sort_distance_files(str(tmp_path / "empty") + "/")


class TestChromosomeRank:
    def test_natural_order(self):
        chroms = ["X", "10", "2", "chr1", "MT", "Y", "9", "GL000", "1"]
        ordered = sorted(chroms, key=module.chromosome_rank)
        assert ordered[:4] == ["chr1", "1", "2", "9"] or ordered[:4] == ["1", "chr1", "2", "9"]
        assert ordered[4:] == ["10", "X", "Y", "MT", "GL000"]
