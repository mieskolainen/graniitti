# Unit tests for weighted iceplot ROC observables
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.


import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")


from core.plot import plot


# Compute a minimal generic two-score ROC observable
def make_roc_observable():
    return {
        "tag": "test_roc",
        "kind": "roc",
        "xlim": (1.0e-3, 1.0),
        "ylim": (0.0, 1.0),
        "xscale": "log",
        "yscale": "linear",
        "xlabel": "FPR",
        "ylabel": "TPR",
        "units": {"x": "unit", "y": "unit"},
        "roc": {
            "score_labels": ["score A", "score B"],
            "directions": ["higher", "higher"],
            "score_colors": ["tab:orange", "tab:green"],
            "comparisons": [
                {
                    "signal": 2,
                    "background": 0,
                    "color_sample": 0,
                    "linestyle": "-",
                    "label": "signal 2 vs background",
                },
                {
                    "signal": 1,
                    "background": 0,
                    "linestyle": "--",
                    "label": "signal 1 vs background",
                },
            ],
            "show_auc": True,
            "show_diagonal": True,
        },
    }


# Compute one sample object accepted by the generic ROC plotter
def make_roc_sample(observable, scores, weights, label, color):
    return {
        "rocdata": plot.prepare_roc_data(scores, weights, score_count=2),
        "obs": observable,
        "label": label,
        "color": color,
    }


# Verify exact weighted cumulative efficiencies and their standard AUC
def test_weighted_roc_curve_event_weights():
    fpr, tpr = plot.weighted_roc_curve(
        signal_scores=[0.9, 0.8],
        signal_weights=[1.0, 3.0],
        background_scores=[0.85, 0.1],
        background_weights=[2.0, 2.0],
        higher_is_signal=True,
    )

    np.testing.assert_allclose(fpr, [0.0, 0.0, 0.5, 0.5, 1.0])
    np.testing.assert_allclose(tpr, [0.0, 0.25, 0.25, 1.0, 1.0])
    assert plot.weighted_roc_auc(fpr, tpr) == pytest.approx(0.625)


# Verify score ties are applied simultaneously for signal and background
def test_weighted_roc_curve_groups_tied_scores():
    fpr, tpr = plot.weighted_roc_curve(
        signal_scores=[1.0, 0.0],
        signal_weights=[1.0, 1.0],
        background_scores=[1.0, 0.0],
        background_weights=[1.0, 1.0],
    )
    np.testing.assert_allclose(fpr, [0.0, 0.5, 1.0])
    np.testing.assert_allclose(tpr, [0.0, 0.5, 1.0])


# Verify lower-valued discriminants can be steered as signal-like
def test_roc_lower_signal():
    fpr, tpr = plot.weighted_roc_curve(
        signal_scores=[0.1, 0.2],
        signal_weights=[1.0, 1.0],
        background_scores=[0.8, 0.9],
        background_weights=[1.0, 1.0],
        higher_is_signal=False,
    )
    assert plot.weighted_roc_auc(fpr, tpr) == pytest.approx(1.0)


# Verify signed event weights are rejected because efficiencies are undefined
def test_weighted_roc_curve_rejects_negative_weights():
    with pytest.raises(ValueError, match="nonnegative event weights"):
        plot.weighted_roc_curve(
            signal_scores=[1.0],
            signal_weights=[-1.0],
            background_scores=[0.0],
            background_weights=[1.0],
        )


# Verify ROC score arrays survive normal histogram chunk fusion
def test_histmc_retains_and_fuses_roc_scores():
    observable = make_roc_observable()
    obs = {"test_roc": observable}
    first = plot.histmc(
        mcdata={"data": {"test_roc": [(0.1, 0.2), (0.3, 0.4)]}, "weights": [1.0, 2.0]},
        obs=obs,
        label="sample",
    )
    second = plot.histmc(
        mcdata={"data": {"test_roc": [(0.5, 0.6)]}, "weights": [3.0]}, obs=obs, label="sample"
    )

    fused = plot.fuse_worker_chunk_outputs([[first], [second]])[0]["test_roc"]["rocdata"]
    np.testing.assert_allclose(fused["scores"], [[0.1, 0.2], [0.3, 0.4], [0.5, 0.6]])
    np.testing.assert_allclose(fused["weights"], [1.0, 2.0, 3.0])


# Verify colors identify scores and line styles identify sample comparisons
def test_rocplot_config_axes_and_comparisons():
    observable = make_roc_observable()
    samples = [
        make_roc_sample(observable, [(0.1, 0.2), (0.2, 0.3)], [1.0, 1.0], "DY", "black"),
        make_roc_sample(observable, [(0.4, 0.5), (0.5, 0.6)], [1.0, 1.0], "SD on", "red"),
        make_roc_sample(observable, [(0.7, 0.8), (0.8, 0.9)], [1.0, 1.0], "SD off", "blue"),
    ]

    figure, axes = plot.rocplot(samples, observable=observable)
    axis = axes[0]
    assert axis.get_xscale() == "log"
    assert axis.get_yscale() == "linear"
    np.testing.assert_allclose(axis.get_xlim(), [1.0e-3, 1.0])
    np.testing.assert_allclose(axis.get_ylim(), [0.0, 1.0])
    assert len(axis.lines) == 5
    assert [line.get_color() for line in axis.lines[:4]] == [
        "tab:orange",
        "tab:green",
        "tab:orange",
        "tab:green",
    ]
    assert [line.get_linestyle() for line in axis.lines[:4]] == ["-", "-", "--", "--"]
    figure.clear()


# Verify score color steering has one entry per ROC discriminant
def test_roc_score_colors_require_color_per_score():
    observable = make_roc_observable()
    observable["roc"]["score_colors"] = ["tab:orange"]
    with pytest.raises(ValueError, match="score colors must match"):
        plot.roc_score_metadata(observable)
