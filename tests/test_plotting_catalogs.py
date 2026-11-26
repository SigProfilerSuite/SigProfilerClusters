"""Unit tests for plottingFunctions.plot96_same and plotINDEL_same: the SBS96 and ID83
catalog panels (all, clustered and non-clustered mutations).
Inputs are synthetic MatrixGenerator-style matrices built with tmp_path; the tests read the
bars drawn on each panel (Agg backend). No genome reference required.
Run with: pytest tests/test_plotting_catalogs.py -v
"""

import itertools

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pytest

from SigProfilerClusters.plottingFunctions import plot96_same, plotINDEL_same


def sbs96_labels():
    labels = []
    for ref, alts in (("C", "AGT"), ("T", "ACG")):
        for alt in alts:
            for left, right in itertools.product("ACGT", "ACGT"):
                labels.append(f"{left}[{ref}>{alt}]{right}")
    return labels


def id83_labels():
    labels = []
    for size, kind, base in [("1", "Del", "C"), ("1", "Del", "T"), ("1", "Ins", "C"), ("1", "Ins", "T")]:
        labels += [f"{size}:{kind}:{base}:{i}" for i in range(6)]
    for kind in ("Del", "Ins"):
        for size in "2345":
            labels += [f"{size}:{kind}:R:{i}" for i in range(6)]
    for size, repeats in (("2", [1]), ("3", [1, 2]), ("4", [1, 2, 3]), ("5", [1, 2, 3, 4, 5])):
        labels += [f"{size}:Del:M:{i}" for i in repeats]
    return labels


CATALOGS = {
    "SBS96": (plot96_same, sbs96_labels()),
    "ID83": (plotINDEL_same, id83_labels()),
}
PANEL_LABELS = ["All Mutations", "Clustered", "Non-Clustered"]


def _write(path, labels, columns):
    """columns maps a sample name to the list of values, one per label."""
    samples = list(columns)
    lines = ["MutationType\t" + "\t".join(samples)]
    for i, label in enumerate(labels):
        lines.append(label + "\t" + "\t".join(str(columns[s][i]) for s in samples))
    path.write_text("\n".join(lines) + "\n")
    return str(path)


def _matrices(tmp_path, labels, all_values, clustered_values, non_clustered_values, sample="S1", extra=None):
    paths = []
    for name, values in (("all", all_values), ("clustered", clustered_values), ("nonclustered", non_clustered_values)):
        columns = {sample: values}
        if extra is not None:
            columns = {"OTHER": extra, sample: values}
        paths.append(_write(tmp_path / f"{name}.txt", labels, columns))
    return paths


def _plot(kind, paths, sample="S1", percentage=False, signature=False, **kwargs):
    function, _ = CATALOGS[kind]
    fig, axes = plt.subplots(3, 1)
    function(*paths, sample, percentage, signature, axes[0], axes[1], axes[2], fig, **kwargs)
    return fig, axes


def _bars(axis):
    """Heights of the catalog bars (the only patches drawn with zorder 1000)."""
    return [float(p.get_height()) for p in axis.patches if p.get_zorder() == 1000]


@pytest.fixture(autouse=True)
def close_figures():
    yield
    plt.close("all")


@pytest.fixture(params=sorted(CATALOGS))
def catalog(request):
    return request.param, CATALOGS[request.param][1]


# ─── bars and counts ─────────────────────────────────────────


class TestCatalogBars:
    def test_one_bar_per_mutation_class_in_every_panel(self, tmp_path, catalog):
        kind, labels = catalog
        n = len(labels)
        paths = _matrices(tmp_path, labels, [1] * n, [1] * n, [1] * n)
        _, axes = _plot(kind, paths)
        assert [len(_bars(axis)) for axis in axes] == [n, n, n]

    def test_class_counts_are_the_bar_heights(self, tmp_path, catalog):
        kind, labels = catalog
        n = len(labels)
        paths = _matrices(tmp_path, labels, list(range(1, n + 1)), [2] * n, [3] * n)
        _, axes = _plot(kind, paths)
        assert _bars(axes[0]) == [float(i) for i in range(1, n + 1)]
        assert _bars(axes[1]) == [2.0] * n
        assert _bars(axes[2]) == [3.0] * n

    def test_panel_axis_is_taller_than_its_tallest_bar(self, tmp_path, catalog):
        kind, labels = catalog
        n = len(labels)
        paths = _matrices(tmp_path, labels, list(range(1, n + 1)), [2] * n, [3] * n)
        _, axes = _plot(kind, paths)
        for axis in axes:
            assert axis.get_ylim()[1] > max(_bars(axis))

    def test_panels_are_labelled(self, tmp_path, catalog):
        kind, labels = catalog
        n = len(labels)
        _, axes = _plot(kind, _matrices(tmp_path, labels, [1] * n, [1] * n, [1] * n))
        assert [axis.get_ylabel() for axis in axes] == PANEL_LABELS

    def test_figure_title_is_the_sample_name(self, tmp_path, catalog):
        kind, labels = catalog
        n = len(labels)
        fig, _ = _plot(kind, _matrices(tmp_path, labels, [1] * n, [1] * n, [1] * n))
        assert fig._suptitle.get_text() == "S1"

    def test_the_requested_sample_column_is_used(self, tmp_path, catalog):
        kind, labels = catalog
        n = len(labels)
        paths = _matrices(tmp_path, labels, [5] * n, [5] * n, [5] * n, sample="S1", extra=[9] * n)
        _, axes = _plot(kind, paths, sample="OTHER")
        assert _bars(axes[0]) == [9.0] * n

    def test_unknown_sample_raises(self, tmp_path, catalog):
        kind, labels = catalog
        n = len(labels)
        paths = _matrices(tmp_path, labels, [1] * n, [1] * n, [1] * n)
        with pytest.raises(ValueError):
            _plot(kind, paths, sample="MISSING")

    def test_missing_matrix_raises(self, tmp_path, catalog):
        kind, labels = catalog
        n = len(labels)
        paths = _matrices(tmp_path, labels, [1] * n, [1] * n, [1] * n)
        paths[1] = str(tmp_path / "nope.txt")
        with pytest.raises(FileNotFoundError):
            _plot(kind, paths)


# ─── percentage and signature options ────────────────────────


class TestPercentageAndSignature:
    def _fractions(self, n):
        return [1.0 / n] * n

    def test_percentage_matrix_keeps_the_given_fractions(self, tmp_path, catalog):
        kind, labels = catalog
        n = len(labels)
        fractions = self._fractions(n)
        paths = _matrices(tmp_path, labels, fractions, fractions, fractions)
        _, axes = _plot(kind, paths, percentage=True)
        assert _bars(axes[0]) == pytest.approx(fractions)

    def test_signature_of_counts_is_normalized_to_percentages(self, tmp_path, catalog):
        kind, labels = catalog
        n = len(labels)
        paths = _matrices(tmp_path, labels, [3] * n, [3] * n, [3] * n)
        _, axes = _plot(kind, paths, signature=True)
        for axis in axes:
            assert sum(_bars(axis)) == pytest.approx(100)

    def test_signature_of_a_percentage_matrix_is_scaled_to_percent(self, tmp_path, catalog):
        kind, labels = catalog
        n = len(labels)
        fractions = self._fractions(n)
        paths = _matrices(tmp_path, labels, fractions, fractions, fractions)
        _, axes = _plot(kind, paths, percentage=True, signature=True)
        for axis in axes:
            assert sum(_bars(axis)) == pytest.approx(100)

    def test_fractions_without_the_percentage_option_raise(self, tmp_path, catalog):
        kind, labels = catalog
        n = len(labels)
        fractions = self._fractions(n)
        paths = _matrices(tmp_path, labels, fractions, fractions, fractions)
        with pytest.raises(ValueError):
            _plot(kind, paths, percentage=False)


# ─── format specific behavior ────────────────────────────────


class TestSbs96:
    def test_chromosome_is_added_to_the_title(self, tmp_path):
        labels = sbs96_labels()
        paths = _matrices(tmp_path, labels, [1] * 96, [1] * 96, [1] * 96)
        fig, _ = _plot("SBS96", paths, chrom="7")
        assert fig._suptitle.get_text() == "S1(7)"

    def test_bars_follow_the_substitution_groups(self, tmp_path):
        labels = sbs96_labels()
        counts = [10] * 16 + [20] * 16 + [30] * 16 + [40] * 16 + [50] * 16 + [60] * 16
        paths = _matrices(tmp_path, labels, counts, [1] * 96, [1] * 96)
        _, axes = _plot("SBS96", paths)
        assert _bars(axes[0]) == [float(c) for c in counts]


class TestId83:
    def test_unknown_indel_classes_are_ignored(self, tmp_path):
        labels = id83_labels() + ["6:Del:R:0"]
        n = len(labels)
        paths = _matrices(tmp_path, labels, [1] * n, [1] * n, [1] * n)
        _, axes = _plot("ID83", paths)
        assert [len(_bars(axis)) for axis in axes] == [83, 83, 83]

    def test_microhomology_classes_start_at_one_repeat(self, tmp_path):
        labels = id83_labels()
        counts = [0] * 83
        counts[labels.index("5:Del:M:5")] = 7
        paths = _matrices(tmp_path, labels, counts, [1] * 83, [1] * 83)
        _, axes = _plot("ID83", paths)
        assert _bars(axes[0])[-1] == 7.0


# ─── matplotlib warnings (audit finding 122) ─────────────────

import warnings

from SigProfilerClusters import plottingFunctions

OPTIONS = [(False, False), (True, False), (False, True)]  # (percentage, signature)


def _draw(kind, paths, **kwargs):
    """Draw the catalog turning the two matplotlib warnings of finding 122 into errors"""
    with warnings.catch_warnings():
        warnings.filterwarnings("error", message=".*set_ticklabels.*")
        warnings.filterwarnings("error", message=".*identical low and high ylims.*")
        fig, axes = _plot(kind, paths, **kwargs)
        fig.canvas.draw()
    return fig, axes


class TestMatplotlibWarnings:
    @pytest.mark.parametrize("percentage, signature", OPTIONS)
    def test_the_labels_of_the_y_axis_are_set_after_its_ticks(self, tmp_path, catalog, percentage, signature):
        kind, labels = catalog
        n = len(labels)
        paths = _matrices(tmp_path, labels, list(range(1, n + 1)), [2] * n, [3] * n)
        _draw(kind, paths, percentage=percentage, signature=signature)

    def test_a_matrix_of_zeros_does_not_ask_for_identical_limits(self, tmp_path, catalog):
        kind, labels = catalog
        n = len(labels)
        paths = _matrices(tmp_path, labels, [0] * n, [0] * n, [0] * n)
        _draw(kind, paths)

    def test_helper_keeps_the_limits_that_matplotlib_chose_for_zero(self):
        assert plottingFunctions._panel_ylim(0) == [-0.05, 0.05]
        assert plottingFunctions._panel_ylim(48) == [0, 48]


class TestYAxisIsUnchanged:
    """What the figures show must stay as it was: these pass before and after the change"""

    @pytest.mark.parametrize("percentage, signature", OPTIONS)
    def test_every_label_is_the_text_of_its_tick(self, tmp_path, catalog, percentage, signature):
        kind, labels = catalog
        n = len(labels)
        paths = _matrices(tmp_path, labels, list(range(1, n + 1)), [2] * n, [3] * n)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fig, axes = _plot(kind, paths, percentage=percentage, signature=signature)
            fig.canvas.draw()
        for axis in axes:
            ticks = [float(t) for t in axis.get_yticks()]
            assert ticks[0] == 0
            steps = {round(b - a, 6) for a, b in zip(ticks, ticks[1:])}
            assert len(steps) == 1  # equally spaced
            texts = [t.get_text() for t in axis.get_yticklabels()]
            assert len(texts) == len(ticks)
            assert texts[0] == "0"
            assert {t.get_color() for t in axis.get_yticklabels()} == {"black"}

    def test_the_y_axis_ends_at_the_top_tick(self, tmp_path, catalog):
        kind, labels = catalog
        n = len(labels)
        paths = _matrices(tmp_path, labels, list(range(1, n + 1)), [2] * n, [3] * n)
        _, axes = _plot(kind, paths)
        for axis in axes:
            assert axis.get_ylim()[0] == 0
            assert axis.get_ylim()[1] == pytest.approx(float(axis.get_yticks()[-1]))

    def test_a_matrix_of_zeros_keeps_its_limits(self, tmp_path, catalog):
        kind, labels = catalog
        n = len(labels)
        paths = _matrices(tmp_path, labels, [0] * n, [0] * n, [0] * n)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            _, axes = _plot(kind, paths)
        for axis in axes:
            assert tuple(round(v, 6) for v in axis.get_ylim()) == (-0.05, 0.05)
