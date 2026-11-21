"""Unit tests for SigProfilerClusters.eventProbability: the probability of every clustered
event, estimated from the density of simulated mutations in the same genomic window.
Inputs are synthetic simulation IMD files and clustered per-sample VCFs built with
tmp_path. No genome reference required.
Tests marked xfail record behavior that looks wrong and is planned for a later fix.
Run with: pytest tests/test_event_probability.py -v
"""

import matplotlib

matplotlib.use("Agg")

import numpy as np
import pytest

from SigProfilerClusters.SigProfilerClusters import eventProbability

WINDOW = 1000  # eventProbability multiplies it by 10: windows of 10 kb
K = 10  # simulation files sampled per sample
HEADER = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"


def _write_simulations(root, sample="S1", chroms=("1", "X"), n_files=K, seed=1):
    rng = np.random.default_rng(seed)
    folder = root / "simulations" / sample
    folder.mkdir(parents=True)
    for i in range(n_files):
        lines = []
        for chrom in chroms:
            positions = sorted(int(x) for x in rng.integers(1, 30000, size=300))
            lines += [f"5\t{sample}_{i}\t{chrom}\t{p}\tC\tT\t5\n" for p in positions]
        (folder / f"{sample}_{i}_96_intradistance.txt").write_text("".join(lines))
    return root / "simulations"


def _event(positions, group, chrom="1"):
    info = f"ClusteredAnalysis;VAF=0.3;groupNumber{group}"
    return [(chrom, pos, "S1", "A", "C", info) for pos in positions]


def _write_clustered(root, events, subclass="DBS", sample="S1"):
    folder = root / "output" / "clustered" / subclass
    folder.mkdir(parents=True, exist_ok=True)
    rows = [row for event in events for row in event]
    text = HEADER + "".join(f"{c}\t{p}\t{i}\t{r}\t{a}\t.\t.\t{info}\n" for c, p, i, r, a, info in rows)
    (folder / f"{sample}.vcf").write_text(text)
    return folder


def _run(root, simulations, k=K):
    eventProbability("proj", str(root) + "/", str(simulations), str(root / "orig") + "/", WINDOW, k)


def _read(path):
    return [line.split("\t") for line in path.read_text().splitlines()]


@pytest.fixture
def setup(tmp_path):
    simulations = _write_simulations(tmp_path)
    folder = _write_clustered(tmp_path, [_event([1000, 1010], 1), _event([2000, 2500], 2)])
    return tmp_path, simulations, folder


class TestEventProbability:
    def test_a_probability_file_replaces_the_clustered_vcf(self, setup):
        root, simulations, folder = setup
        _run(root, simulations)
        assert sorted(p.name for p in folder.iterdir()) == ["S1_prob.vcf"]

    def test_header_has_the_probability_column(self, setup):
        root, simulations, folder = setup
        _run(root, simulations)
        header = _read(folder / "S1_prob.vcf")[0]
        assert header[0] == "#CHROM" and header[-1] == "EVENT_PROBABILITY"

    def test_every_mutation_is_kept_with_its_original_columns(self, setup):
        root, simulations, folder = setup
        _run(root, simulations)
        rows = _read(folder / "S1_prob.vcf")[1:]
        assert sorted((r[0], int(r[1])) for r in rows) == [("1", 1000), ("1", 1010), ("1", 2000), ("1", 2500)]
        assert all(r[7].startswith("ClusteredAnalysis;VAF=0.3;groupNumber") for r in rows)

    def test_mutations_of_an_event_share_its_probability(self, setup):
        root, simulations, folder = setup
        _run(root, simulations)
        by_event = {}
        for row in _read(folder / "S1_prob.vcf")[1:]:
            by_event.setdefault(row[7], set()).add(row[-1])
        assert len(by_event) == 2
        assert all(len(values) == 1 for values in by_event.values())

    def test_probabilities_are_between_the_floor_and_one(self, setup):
        root, simulations, folder = setup
        _run(root, simulations)
        probs = [float(r[-1]) for r in _read(folder / "S1_prob.vcf")[1:]]
        assert all(1e-5 <= p <= 1 for p in probs)

    def test_very_unlikely_events_are_floored(self, tmp_path):
        simulations = _write_simulations(tmp_path)
        folder = _write_clustered(tmp_path, [_event(range(1000, 1100, 10), 1)])
        _run(tmp_path, simulations)
        probs = {float(r[-1]) for r in _read(folder / "S1_prob.vcf")[1:]}
        assert probs == {1e-5}

    def test_every_clustered_subclass_is_processed(self, tmp_path):
        simulations = _write_simulations(tmp_path)
        dbs = _write_clustered(tmp_path, [_event([1000, 1010], 1)], subclass="DBS")
        omikli = _write_clustered(tmp_path, [_event([3000, 3400], 2)], subclass="omikli")
        _run(tmp_path, simulations)
        assert (dbs / "S1_prob.vcf").exists() and (omikli / "S1_prob.vcf").exists()

    def test_subclasses_without_a_vcf_are_skipped(self, setup):
        root, simulations, _ = setup
        _run(root, simulations)
        assert not (root / "output" / "clustered" / "kataegis").exists()

    def test_hidden_entries_of_the_simulation_folder_are_ignored(self, setup):
        root, simulations, folder = setup
        (simulations / ".DS_Store").write_text("")
        _run(root, simulations)
        assert (folder / "S1_prob.vcf").exists()

    def test_fewer_simulation_files_than_requested_raises(self, tmp_path):
        simulations = _write_simulations(tmp_path, n_files=3)
        _write_clustered(tmp_path, [_event([1000, 1010], 1)])
        with pytest.raises(ValueError):
            _run(tmp_path, simulations, k=K)

    def test_only_the_requested_number_of_simulations_is_needed(self, tmp_path):
        simulations = _write_simulations(tmp_path, n_files=3)
        folder = _write_clustered(tmp_path, [_event([1000, 1010], 1)])
        _run(tmp_path, simulations, k=3)
        assert (folder / "S1_prob.vcf").exists()

    @pytest.mark.xfail(
        reason="when every chromosome name is numeric pandas reads them as integers, so "
        "the distributions are keyed by int but looked up by str(chrom): KeyError",
        raises=KeyError,
        strict=True,
    )
    def test_numeric_chromosomes_only(self, tmp_path):
        simulations = _write_simulations(tmp_path, chroms=("1",))
        folder = _write_clustered(tmp_path, [_event([1000, 1010], 1)])
        _run(tmp_path, simulations)
        assert (folder / "S1_prob.vcf").exists()
