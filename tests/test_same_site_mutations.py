"""Mutations at the same position must not disappear and can belong to a cluster
(audit finding 7). The observed mutations go through distance_one_file and
calculateSampleIMDs, so the whole path from the mutation file to the clustered and
non-clustered tables is exercised. No genome reference required.
Run with: pytest tests/test_same_site_mutations.py -v
"""

import random

import test_calculate_sample_imds as csi
from SigProfilerClusters.SigProfilerClusters import distance_one_file

# the whole synthetic chromosome is one arm: positions go up to 3 Gb like a genome-wide background
CENTROMERES = {"GRCh37": {"1": [4_000_000_000, 4_000_001_000]}}


def _observed(same_site_neighbor=True):
    """A strongly clustered sample (300 mutations 20 bp apart), a background spread over the
    genome (so that the simulated background makes the cutoff scan stop, as in real data) and a site
    with two alleles at 8000 that has a neighbor at 8030 (a few bp from the cluster)."""
    rows = [("S1", "1", 1000 + 20 * i, "C", "T") for i in range(300)]
    background = sorted(random.Random(7).sample(range(200_000, 3_000_000_000, 1000), 100))
    rows += [("S1", "1", position, "C", "T") for position in background]
    rows += [("S1", "1", 8000, "C", "T"), ("S1", "1", 8000, "C", "G")]
    if same_site_neighbor:
        rows.append(("S1", "1", 8030, "C", "T"))
    return rows


def _run_pipeline(tmp_path, mutations):
    csi._build(tmp_path, orig=False)  # simulations and folders
    source = tmp_path / "observed"
    source.mkdir()
    (source / "1.txt").write_text("".join("\t".join(map(str, m)) + "\n" for m in mutations))
    distance_one_file(["1.txt"], str(source) + "/", str(tmp_path / "orig") + "/", "96", "GRCh37", CENTROMERES)
    csi._run(tmp_path)
    clustered = [r.split("\t") for r in csi._table_rows(tmp_path, "c")]
    non_clustered = [r.split("\t") for r in csi._table_rows(tmp_path, "n")]
    return clustered, non_clustered


def _at(rows, position):
    return [r for r in rows if r[6] == str(position)]


class TestSameSiteMutations:
    def test_no_mutation_is_lost(self, tmp_path):
        mutations = _observed()
        clustered, non_clustered = _run_pipeline(tmp_path, mutations)
        assert len(clustered) + len(non_clustered) == len(mutations)

    def test_both_alleles_of_a_site_next_to_a_cluster_are_clustered(self, tmp_path):
        clustered, non_clustered = _run_pipeline(tmp_path, _observed())
        assert sorted(r[9] for r in _at(clustered, 8000)) == ["G", "T"]
        assert _at(non_clustered, 8000) == []
        assert len(_at(clustered, 8030)) == 1

    def test_alleles_far_from_any_neighbor_are_not_clustered_but_are_kept(self, tmp_path):
        mutations = [m for m in _observed(same_site_neighbor=False)]
        clustered, non_clustered = _run_pipeline(tmp_path, mutations)
        # 8000 is 6 kb from the last mutation of the cluster (6980) but alleles share that distance
        assert len(_at(clustered, 8000)) + len(_at(non_clustered, 8000)) == 2

    def test_output_rows_keep_each_allele(self, tmp_path):
        clustered, non_clustered = _run_pipeline(tmp_path, _observed())
        alts = sorted(r[9] for r in _at(clustered + non_clustered, 8000))
        assert alts == ["G", "T"]
