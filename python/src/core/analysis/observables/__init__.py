# Shared observable dictionary construction
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>

# Build one histogram observable with axis, ratio, binning and projection metadata
def histogram(
    *,
    tag,
    xlim,
    xlabel,
    ylabel,
    xunit,
    label,
    bins,
    func,
    ylim=None,
    yunit="pb",
    ylim_ratio=(0.0, 2.0),
    ytick_ratio_step=0.5,
    density=False,
    **metadata,
):
    observable = {
        "tag": tag,
        "xlim": xlim,
        "ylim": ylim,
        "xlabel": xlabel,
        "ylabel": ylabel,
        "units": {"x": xunit, "y": yunit},
        "label": label,
        "ylim_ratio": ylim_ratio,
        "ytick_ratio_step": ytick_ratio_step,
        "bins": bins,
        "density": density,
        **metadata,
        "func": func,
    }
    return observable
