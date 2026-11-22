"""Unit tests for the IMD histogram helpers of hotspot: plot_clustered, plot_non_clustered
and plot_hist. They draw on matplotlib axes (Agg backend); the tests check the artists
that end up on each panel. No genome reference required.
Run with: pytest tests/test_hotspot_plots.py -v
"""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

from SigProfilerClusters.hotspot import plot_clustered, plot_hist, plot_non_clustered

BINS = np.array([2.0**i for i in range(1, 11)])  # 2 ... 1024
ORIG = np.array([50, 40, 30, 20, 10, 8, 6, 4, 2, 1])
SIM = np.array([5, 5, 5, 5, 5, 5, 5, 5, 5, 5])
LOWER = [3] * 10
UPPER = [8] * 10


@pytest.fixture
def panels():
    fig, axes = plt.subplots(5, 1)
    yield axes
    plt.close(fig)


def _colors(axis):
    return [line.get_color() for line in axis.get_lines()]


# ─── plot_clustered / plot_non_clustered ─────────────────────


class TestPlotClustered:
    def test_draws_simulated_and_original_lines_and_the_interval(self, panels):
        panel = panels[0]
        plot_clustered(ORIG[:4], SIM[:4], BINS, panel, LOWER, UPPER)
        assert _colors(panel) == ["red", "green"]
        assert len(panel.collections) == 1

    def test_uses_the_first_bins(self, panels):
        panel = panels[0]
        plot_clustered(ORIG[:4], SIM[:4], BINS, panel, LOWER, UPPER)
        sim, orig = panel.get_lines()
        assert list(sim.get_xdata()) == list(BINS[:4])
        assert list(orig.get_ydata()) == list(ORIG[:4])

    def test_axes_are_logarithmic(self, panels):
        panel = panels[0]
        plot_clustered(ORIG[:4], SIM[:4], BINS, panel, LOWER, UPPER)
        assert panel.get_xscale() == "log" and panel.get_yscale() == "log"


class TestPlotNonClustered:
    def test_draws_simulated_and_original_lines_and_the_interval(self, panels):
        panel = panels[0]
        plot_non_clustered(ORIG[4:], SIM[4:], BINS, panel, LOWER, UPPER)
        assert _colors(panel) == ["red", "green"]
        assert len(panel.collections) == 1

    def test_uses_the_last_bins(self, panels):
        panel = panels[0]
        plot_non_clustered(ORIG[4:], SIM[4:], BINS, panel, LOWER, UPPER)
        sim, orig = panel.get_lines()
        assert list(sim.get_xdata()) == list(BINS[4:])
        assert list(orig.get_ydata()) == list(ORIG[4:])

    def test_axes_are_logarithmic_and_labelled(self, panels):
        panel = panels[0]
        plot_non_clustered(ORIG[4:], SIM[4:], BINS, panel, LOWER, UPPER)
        assert panel.get_xscale() == "log" and panel.get_yscale() == "log"
        assert panel.get_xlabel() == "Inter-Mutational Distance (IMD)"


# ─── plot_hist ───────────────────────────────────────────────


def _plot_hist(panels, original=True, interval_line=4):
    panel2, panel3, panel4, panel5, panel6 = panels
    return plot_hist(
        ORIG,
        BINS,
        1e-5,
        interval_line,
        sum(ORIG[:interval_line]),
        12.0,
        2.0,
        BINS[interval_line],
        LOWER,
        UPPER,
        9,
        15,
        SIM,
        "S1",
        original,
        panel2,
        panel3,
        panel4,
        panel5,
        panel6,
    )


class TestPlotHist:
    def test_returns_true(self, panels):
        assert _plot_hist(panels) is True

    def test_main_panel_has_the_cutoff_and_the_simulated_and_original_lines(self, panels):
        _plot_hist(panels)
        assert _colors(panels[0]) == ["darkred", "red", "green"]
        assert panels[0].get_xscale() == "log" and panels[0].get_yscale() == "log"

    def test_cutoff_line_is_drawn_at_the_refined_imd(self, panels):
        _plot_hist(panels)
        cutoff = BINS[4]
        assert any(
            list(line.get_xdata()) == [cutoff, cutoff] for line in panels[0].get_lines()
        )

    def test_legend_names_the_series_and_the_q_value(self, panels):
        _plot_hist(panels)
        labels = [text.get_text() for text in panels[0].get_legend().get_texts()]
        assert labels == ["simulated", "real samples", "q_value = 1.00E-5"]

    def test_summary_texts_report_the_counts(self, panels):
        _plot_hist(panels)
        texts = [t.get_text() for p in panels for t in p.texts]
        assert "Total=" + str(sum(ORIG)) in texts
        assert "Clustered=" + str(sum(ORIG[:4])) + "(IMD<" + str(BINS[4]) + ")" in texts
        assert "Expected=12.0(9-15)" in texts
        assert "Non-Clustered=" + str(sum(ORIG) - sum(ORIG[:4])) in texts

    def test_clustered_and_non_clustered_panels_split_at_the_interval(self, panels):
        _plot_hist(panels, interval_line=4)
        clustered = panels[2].get_lines()[0].get_xdata()
        non_clustered = panels[4].get_lines()[0].get_xdata()
        assert list(clustered) == list(BINS[:4])
        assert list(non_clustered) == list(BINS[4:])

    def test_interval_at_the_first_bin_still_plots_one_clustered_bin(self, panels):
        _plot_hist(panels, interval_line=0)
        assert len(panels[2].get_lines()[0].get_xdata()) == 1

    def test_without_the_original_sample_only_the_simulation_is_drawn(self, panels):
        _plot_hist(panels, original=False)
        assert _colors(panels[0]) == ["red"]
        assert [t.get_text() for t in panels[0].get_legend().get_texts()] == ["simulated"]
        assert [t for p in panels for t in p.texts] == []
