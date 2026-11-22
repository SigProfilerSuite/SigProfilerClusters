"""Unit tests for the inter-mutational distance (IMD) calculation of SigProfilerClusters:
distance_one_file (original samples) and distance_multiple_files_sims (simulations).
Inputs are tiny synthetic mutation files built with tmp_path; the centromere of every
chromosome is placed at 5000-6000. No genome reference required.
Tests marked xfail record behavior that looks wrong and is planned for a later fix.
Run with: pytest tests/test_distance_functions.py -v
"""

import pytest

from SigProfilerClusters.SigProfilerClusters import (
    distance_multiple_files_sims,
    distance_one_file,
)

GENOME = "GRCh37"
CENTROMERES = {GENOME: {"1": [5000, 6000], "2": [5000, 6000], "X": [5000, 6000]}}
CONTEXT = "96"


def _read(path):
    return [line.split("\t") for line in path.read_text().splitlines()]


# ─── distance_one_file (original samples) ────────────────────


def _mutation(pos, sample="S1", chrom="1", ref="C", alt="T"):
    return [sample, chrom, str(pos), ref, alt]


def _run_original(tmp_path, files):
    """files maps a file name to its list of mutations; returns {output name: rows}."""
    source, out = tmp_path / "original", tmp_path / "out"
    source.mkdir()
    out.mkdir()
    for name, rows in files.items():
        (source / name).write_text("".join("\t".join(r) + "\n" for r in rows))
    distance_one_file(list(files), str(source) + "/", str(out) + "/", CONTEXT, GENOME, CENTROMERES)
    return {path.name: _read(path) for path in sorted(out.iterdir())}


def _imds(rows):
    return [int(row[0]) for row in rows]


def _positions(rows):
    return [int(row[3]) for row in rows]


class TestDistanceOneFile:
    def test_imd_is_the_distance_to_the_closest_neighbor(self, tmp_path):
        out = _run_original(tmp_path, {"1.txt": [_mutation(100), _mutation(110), _mutation(1000)]})
        rows = out["S1_96_intradistance.txt"]
        assert _imds(rows) == [10, 10, 890]
        assert _positions(rows) == [100, 110, 1000]

    def test_adjacent_mutations_have_an_imd_of_one(self, tmp_path):
        out = _run_original(tmp_path, {"1.txt": [_mutation(100), _mutation(101)]})
        assert _imds(out["S1_96_intradistance.txt"]) == [1, 1]

    def test_q_arm_mutations_are_measured_the_same_way(self, tmp_path):
        out = _run_original(tmp_path, {"1.txt": [_mutation(7000), _mutation(7005), _mutation(7100)]})
        assert _imds(out["S1_96_intradistance.txt"]) == [5, 5, 95]

    def test_each_arm_is_measured_on_its_own(self, tmp_path):
        rows = [_mutation(100), _mutation(110), _mutation(7000), _mutation(7005)]
        out = _run_original(tmp_path, {"1.txt": rows})
        assert _imds(out["S1_96_intradistance.txt"]) == [10, 10, 5, 5]

    def test_a_pair_that_spans_the_centromere_uses_their_distance(self, tmp_path):
        out = _run_original(tmp_path, {"1.txt": [_mutation(4000), _mutation(7000)]})
        assert _imds(out["S1_96_intradistance.txt"]) == [3000, 3000]

    def test_unsorted_input_is_sorted_by_position(self, tmp_path):
        rows = [_mutation(1000), _mutation(100), _mutation(110)]
        out = _run_original(tmp_path, {"1.txt": rows})
        assert _positions(out["S1_96_intradistance.txt"]) == [100, 110, 1000]

    def test_reference_and_alternate_lengths_are_taken_into_account(self, tmp_path):
        # a deletion AC>A at 100 ends at 101: the next mutation at 110 is 9 bp away
        out = _run_original(tmp_path, {"1.txt": [_mutation(100, ref="AC", alt="A"), _mutation(110)]})
        assert _imds(out["S1_96_intradistance.txt"]) == [9, 9]

    def test_columns_of_the_output_rows(self, tmp_path):
        out = _run_original(tmp_path, {"1.txt": [_mutation(100), _mutation(110)]})
        row = out["S1_96_intradistance.txt"][0]
        assert row[1:6] == ["S1", "1", "100", "C", "T"]
        assert len(row) == 7

    def test_each_sample_gets_its_own_file(self, tmp_path):
        rows = [_mutation(100, "S1"), _mutation(110, "S1"), _mutation(100, "S2"), _mutation(150, "S2")]
        out = _run_original(tmp_path, {"1.txt": rows})
        assert sorted(out) == ["S1_96_intradistance.txt", "S2_96_intradistance.txt"]
        assert _imds(out["S2_96_intradistance.txt"]) == [50, 50]

    def test_a_single_mutation_gets_an_imd_of_one_million(self, tmp_path):
        out = _run_original(tmp_path, {"1.txt": [_mutation(100)]})
        rows = out["S1_96_intradistance.txt"]
        assert rows == [["1000000", "S1", "1", "100", "C", "T"]]

    def test_files_of_several_chromosomes_are_appended_per_sample(self, tmp_path):
        files = {
            "1.txt": [_mutation(100, chrom="1"), _mutation(110, chrom="1")],
            "2.txt": [_mutation(200, chrom="2"), _mutation(260, chrom="2")],
        }
        rows = _run_original(tmp_path, files)["S1_96_intradistance.txt"]
        assert _imds(rows) == [10, 10, 60, 60]
        assert [row[2] for row in rows] == ["1", "1", "2", "2"]

    def test_missing_file_raises(self, tmp_path):
        out = tmp_path / "out"
        out.mkdir()
        with pytest.raises(FileNotFoundError):
            distance_one_file(["nope.txt"], str(tmp_path) + "/", str(out) + "/", CONTEXT, GENOME, CENTROMERES)


class TestDistanceOneFileSameSite:
    """Several alleles at the same position (audit finding 7). The IMD of a mutation is the
    distance to the closest *other site*: a site with two alleles keeps the distance to its
    neighbors instead of getting 0, which used to hide that distance and drop the rows."""

    def _rows(self, tmp_path, mutations, name="S1_96_intradistance.txt"):
        return _run_original(tmp_path, {"1.txt": mutations})[name]

    def test_alleles_of_a_site_next_to_a_third_mutation_share_its_distance(self, tmp_path):
        rows = self._rows(tmp_path, [_mutation(7000, alt="T"), _mutation(7000, alt="G"), _mutation(7050)])
        assert _imds(rows) == [50, 50, 50]

    def test_every_allele_keeps_its_own_columns(self, tmp_path):
        rows = self._rows(tmp_path, [_mutation(7000, alt="T"), _mutation(7000, alt="G"), _mutation(7050)])
        assert [(r[3], r[4], r[5]) for r in rows] == [("7000", "C", "T"), ("7000", "C", "G"), ("7050", "C", "T")]

    def test_a_site_with_distant_neighbors_gets_the_distance_to_them_not_zero(self, tmp_path):
        rows = self._rows(tmp_path, [_mutation(1000), _mutation(4000, alt="T"), _mutation(4000, alt="G")])
        assert _imds(rows) == [3000, 3000, 3000]

    def test_no_row_has_an_imd_of_zero(self, tmp_path):
        rows = self._rows(tmp_path, [_mutation(100), _mutation(100, alt="G"), _mutation(110), _mutation(110, alt="A")])
        assert 0 not in _imds(rows)

    def test_every_input_mutation_has_an_output_row(self, tmp_path):
        mutations = [_mutation(100), _mutation(100, alt="G"), _mutation(130), _mutation(5000, alt="T"), _mutation(5000, alt="A")]
        assert len(self._rows(tmp_path, mutations)) == len(mutations)

    def test_alleles_alone_have_no_neighbor_and_get_the_singleton_distance(self, tmp_path):
        rows = self._rows(tmp_path, [_mutation(100, alt="T"), _mutation(100, alt="G")])
        assert [(r[0], r[5]) for r in rows] == [("1000000", "T"), ("1000000", "G")]

    def test_unsorted_input_with_a_repeated_site(self, tmp_path):
        rows = self._rows(tmp_path, [_mutation(7050), _mutation(7000, alt="G"), _mutation(7000, alt="T")])
        assert sorted((r[3], r[5]) for r in rows) == [("7000", "G"), ("7000", "T"), ("7050", "T")]
        assert _imds(rows) == [50, 50, 50]

    def test_the_same_position_in_two_samples_is_not_the_same_site(self, tmp_path):
        mutations = [
            _mutation(100, "S1"), _mutation(110, "S1"),
            _mutation(100, "S2"), _mutation(150, "S2"),
        ]
        out = _run_original(tmp_path, {"1.txt": mutations})
        assert _imds(out["S1_96_intradistance.txt"]) == [10, 10]
        assert _imds(out["S2_96_intradistance.txt"]) == [50, 50]

    def test_files_without_repeated_sites_are_unchanged(self, tmp_path):
        rows = self._rows(tmp_path, [_mutation(100), _mutation(110), _mutation(1000)])
        assert _imds(rows) == [10, 10, 890]


# ─── distance_multiple_files_sims (simulations) ──────────────


def _sim_row(sample, chrom, pos, ref="C", alt="T"):
    """A simulated mutation: chromosome in column 4, position 5, reference 10,
    alternate 12 and sample (with the simulation number) 15."""
    cols = ["x"] * 16
    cols[4], cols[5], cols[10], cols[12], cols[15] = chrom, str(pos), ref, alt, sample
    return cols


def _run_sims(tmp_path, rows, sort=True, name="1.txt"):
    sims, sorted_sims, out = tmp_path / "sims", tmp_path / "sorted", tmp_path / "out"
    for folder in (sims, sorted_sims, out):
        folder.mkdir()
    (sims / name).write_text("header\n" + "".join("\t".join(r) + "\n" for r in rows))
    distance_multiple_files_sims(
        str(out) + "/", [name], str(sims) + "/", str(sorted_sims) + "/", CONTEXT, GENOME, CENTROMERES, sort
    )
    outputs = {str(p.relative_to(out)): _read(p) for p in sorted(out.rglob("*.txt"))}
    return outputs, sorted_sims / name


class TestDistanceMultipleFilesSims:
    def test_imd_is_the_distance_to_the_closest_neighbor(self, tmp_path):
        rows = [_sim_row("S1_1", "1", 100), _sim_row("S1_1", "1", 110), _sim_row("S1_1", "1", 1000)]
        outputs, _ = _run_sims(tmp_path, rows)
        assert _imds(outputs["S1/S1_1_96_intradistance.txt"]) == [10, 10, 890]

    def test_output_is_written_in_a_folder_per_sample(self, tmp_path):
        rows = [_sim_row("S1_1", "1", 100), _sim_row("S1_1", "1", 110)]
        outputs, _ = _run_sims(tmp_path, rows)
        assert list(outputs) == ["S1/S1_1_96_intradistance.txt"]

    def test_each_simulation_gets_its_own_file(self, tmp_path):
        rows = [
            _sim_row("S1_1", "1", 100),
            _sim_row("S1_1", "1", 110),
            _sim_row("S1_2", "1", 100),
            _sim_row("S1_2", "1", 150),
        ]
        outputs, _ = _run_sims(tmp_path, rows)
        assert sorted(outputs) == ["S1/S1_1_96_intradistance.txt", "S1/S1_2_96_intradistance.txt"]
        assert set(_imds(outputs["S1/S1_1_96_intradistance.txt"])) == {10}
        assert set(_imds(outputs["S1/S1_2_96_intradistance.txt"])) == {50}

    def test_sorted_copy_orders_by_sample_chromosome_and_position(self, tmp_path):
        rows = [
            _sim_row("S1_1", "2", 100),
            _sim_row("S1_1", "1", 1000),
            _sim_row("S1_1", "1", 110),
            _sim_row("S1_1", "X", 50),
            _sim_row("S1_1", "1", 100),
        ]
        _, sorted_file = _run_sims(tmp_path, rows)
        order = [(cols[4], cols[5]) for cols in _read(sorted_file)]
        assert order == [("X", "50"), ("1", "100"), ("1", "110"), ("1", "1000"), ("2", "100")]

    def test_no_sorted_copy_without_sorting(self, tmp_path):
        rows = [_sim_row("S1_1", "1", 100), _sim_row("S1_1", "1", 110)]
        _, sorted_file = _run_sims(tmp_path, rows, sort=False)
        assert not sorted_file.exists()

    def test_a_chromosome_with_a_single_mutation_has_no_imd(self, tmp_path):
        rows = [
            _sim_row("S1_1", "1", 100),
            _sim_row("S1_1", "1", 110),
            _sim_row("S1_1", "2", 100),
        ]
        outputs, _ = _run_sims(tmp_path, rows)
        assert {row[2] for row in outputs["S1/S1_1_96_intradistance.txt"]} == {"1"}

    def test_a_single_mutation_writes_nothing(self, tmp_path):
        outputs, _ = _run_sims(tmp_path, [_sim_row("S1_1", "1", 100)])
        assert outputs == {}

    def test_one_row_per_mutation(self, tmp_path):
        rows = [_sim_row("S1_1", "1", p) for p in (100, 110, 1000, 1100)]
        outputs, _ = _run_sims(tmp_path, rows)
        assert len(outputs["S1/S1_1_96_intradistance.txt"]) == 4

    @pytest.mark.xfail(
        reason="the last mutation of each simulated sample is reported with the position "
        "of the previous one (distances_new repeats the last entry; distance_one_file "
        "appends the real last mutation instead)",
        strict=True,
    )
    def test_every_mutation_keeps_its_own_position(self, tmp_path):
        rows = [_sim_row("S1_1", "1", 100), _sim_row("S1_1", "1", 110), _sim_row("S1_1", "1", 1000)]
        outputs, _ = _run_sims(tmp_path, rows)
        assert _positions(outputs["S1/S1_1_96_intradistance.txt"]) == [100, 110, 1000]

    @pytest.mark.xfail(
        reason="when a simulation file holds several samples, the last sample gets one IMD "
        "row fewer than it has mutations: the boundary marker between samples consumes an entry",
        strict=True,
    )
    def test_every_simulated_sample_gets_one_row_per_mutation(self, tmp_path):
        rows = [
            _sim_row("S1_1", "1", 100),
            _sim_row("S1_1", "1", 110),
            _sim_row("S1_2", "1", 100),
            _sim_row("S1_2", "1", 150),
        ]
        outputs, _ = _run_sims(tmp_path, rows)
        assert len(outputs["S1/S1_1_96_intradistance.txt"]) == 2
        assert len(outputs["S1/S1_2_96_intradistance.txt"]) == 2
