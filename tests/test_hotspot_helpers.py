"""Unit tests for the pure helper functions in SigProfilerClusters.hotspot.
No genome reference required.
Run with: pytest tests/test_hotspot_helpers.py -v
"""

import math

import numpy as np
import pytest

from SigProfilerClusters.hotspot import (
    catch,
    cutoffCatch,
    densityCorrectionOriginal,
    flatten,
    moving_density,
    z_test,
)


# ─── z_test ──────────────────────────────────────────────────


class TestZTest:
    def test_observed_equals_mean(self):
        z, p = z_test(10, 10, 2)
        assert z == 0
        assert p == pytest.approx(1.0)

    def test_z_score_value(self):
        z, _ = z_test(14, 10, 2)
        assert z == pytest.approx(2.0)

    def test_p_value_for_two_sigma(self):
        _, p = z_test(14, 10, 2)
        # the implementation is 2 * (1 - cdf(z)): about 0.0455 for z = 2
        assert p == pytest.approx(0.0455, abs=1e-4)

    def test_p_value_decreases_with_larger_z(self):
        _, p_small = z_test(11, 10, 2)
        _, p_large = z_test(20, 10, 2)
        assert p_large < p_small

    def test_zero_sigma_raises(self):
        with pytest.raises(ZeroDivisionError):
            z_test(5, 3, 0)

    def test_numpy_arrays(self):
        z, p = z_test(np.array([12.0, 10.0]), np.array([10.0, 10.0]), np.array([2.0, 2.0]))
        assert list(z) == [1.0, 0.0]
        assert len(p) == 2


# ─── moving_density ──────────────────────────────────────────


class TestMovingDensity:
    def test_window_is_half_open(self):
        assert moving_density([1, 5, 10, 15], 5, 15) == [5, 10]

    def test_empty_input(self):
        assert moving_density([], 0, 100) == []

    def test_no_values_in_window(self):
        assert moving_density([1, 2, 3], 10, 20) == []

    def test_preserves_order(self):
        assert moving_density([9, 3, 7], 0, 10) == [9, 3, 7]

    def test_returns_list(self):
        assert isinstance(moving_density((1, 2, 3), 0, 10), list)


# ─── flatten ─────────────────────────────────────────────────


class TestFlatten:
    def test_flat_list_unchanged(self):
        assert list(flatten([1, 2, 3])) == [1, 2, 3]

    def test_nested_lists(self):
        assert list(flatten([1, [2, [3, [4]]], 5])) == [1, 2, 3, 4, 5]

    def test_empty(self):
        assert list(flatten([])) == []

    def test_empty_nested(self):
        assert list(flatten([[], [[]]])) == []

    def test_strings_are_not_split(self):
        assert list(flatten(["ab", ["cd"]])) == ["ab", "cd"]

    def test_tuples_are_flattened(self):
        assert list(flatten([(1, 2), [3]])) == [1, 2, 3]


# ─── catch ───────────────────────────────────────────────────


class TestCatch:
    # x = [distance, sample, chrom, position]
    chrom_lengths = {"GRCh37": {"1": 0, "2": 1000}}
    regions = [0, 500, 1000, 1500]

    def test_position_on_first_chromosome(self):
        x = [10, "s1", "1", 400]
        assert catch(x, self.regions, self.chrom_lengths, "GRCh37") == 1

    def test_cumulative_offset_is_added(self):
        x = [10, "s1", "2", 100]
        # 1000 + 100 = 1100 -> first region >= 1100 is index 3 (1500)
        assert catch(x, self.regions, self.chrom_lengths, "GRCh37") == 3

    def test_exact_region_boundary(self):
        x = [10, "s1", "1", 500]
        assert catch(x, self.regions, self.chrom_lengths, "GRCh37") == 1

    def test_past_last_region_is_clamped(self):
        x = [10, "s1", "2", 5000]
        assert catch(x, self.regions, self.chrom_lengths, "GRCh37") == len(self.regions) - 1

    def test_position_is_converted_from_string(self):
        x = [10, "s1", "1", "400"]
        assert catch(x, self.regions, self.chrom_lengths, "GRCh37") == 1

    def test_unknown_chromosome_raises(self):
        with pytest.raises(KeyError):
            catch([10, "s1", "X", 1], self.regions, self.chrom_lengths, "GRCh37")


# ─── cutoffCatch ─────────────────────────────────────────────


class TestCutoffCatch:
    chroms = {"1": 0}
    regions = [0, 1000]

    def test_sample_without_corrected_imds_uses_global_cut(self):
        x = [50, "s1", "1", 10]
        assert cutoffCatch(x, {}, self.regions, 0, 100, self.chroms, 1_000_000) is True
        x = [150, "s1", "1", 10]
        assert cutoffCatch(x, {}, self.regions, 0, 100, self.chroms, 1_000_000) is False

    def test_region_without_correction_uses_global_cut(self):
        imds = {"s1": {1000: 20}}
        x = [50, "s1", "1", 10]
        # regions[0] == 0 is not corrected for s1
        assert cutoffCatch(x, imds, self.regions, 0, 100, self.chroms, 1_000_000) is True

    def test_corrected_region_uses_corrected_cut(self):
        imds = {"s1": {0: 20}}
        x = [50, "s1", "1", 10]
        assert cutoffCatch(x, imds, self.regions, 0, 100, self.chroms, 1_000_000) is False
        x = [10, "s1", "1", 10]
        assert cutoffCatch(x, imds, self.regions, 0, 100, self.chroms, 1_000_000) is True

    def test_position_beyond_window_size_is_never_clustered(self):
        imds = {"s1": {0: 20}}
        x = [1, "s1", "1", 2_000_000]
        assert cutoffCatch(x, imds, self.regions, 0, 100, self.chroms, 1_000_000) is False

    def test_distance_equal_to_cut_is_not_clustered(self):
        x = [100, "s1", "1", 10]
        assert cutoffCatch(x, {}, self.regions, 0, 100, self.chroms, 1_000_000) is False


# ─── densityCorrectionOriginal ───────────────────────────────


class TestDensityCorrectionOriginal:
    bins = np.arange(0, 101, 10)

    def _sims(self, n=3):
        # 10 mutations in the middle of every bin, in every simulation
        values = [b + 5 for b in range(0, 100, 10)] * 10
        return {i: list(values) for i in range(n)}

    def test_dense_region_is_reported(self):
        # 100 mutations in the first bin and 5 in each of the others
        original = [5] * 100 + [b + 5 for b in range(10, 100, 10)] * 5
        result = densityCorrectionOriginal(original, self._sims(), self.bins)
        assert result[0] == 0
        assert result[-1] == 6000000000
        assert len(result) == 2

    def test_no_dense_region_returns_only_sentinel(self):
        original = [b + 5 for b in range(0, 100, 10)] * 10
        result = densityCorrectionOriginal(original, self._sims(), self.bins)
        assert result == [6000000000]

    def test_dense_bin_with_few_mutations_is_ignored(self):
        # fold change is large, but the bin has fewer than 50 mutations
        original = [5] * 20 + [b + 5 for b in range(10, 100, 10)] * 1
        result = densityCorrectionOriginal(original, self._sims(), self.bins)
        assert result == [6000000000]

    def test_sentinel_is_always_last(self):
        original = [b + 5 for b in range(0, 100, 10)] * 10
        result = densityCorrectionOriginal(original, self._sims(), self.bins)
        assert result[-1] == 6000000000
        assert not math.isnan(result[-1])
