# Plot simulation and data histograms with fit diagnostics
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import os

import numpy as np
from matplotlib import pyplot as plt

from core.io import readers
from core.io.files import ensure_dir
from core.io.serialize import write_json_file
from core.plot import plot
from core.stats import objective as icecost
from core.tune import core as icetune


# Index canonical likelihood records by their histogram coordinates
def _likelihood_observable_map(likelihood: dict | None) -> dict[tuple[int, int, str], dict]:
    if not isinstance(likelihood, dict):
        return {}

    records = {}
    for record in likelihood.get("observables", []):
        if not isinstance(record, dict):
            continue
        try:
            key = (int(record["dataset"]), int(record["subset"]), str(record["observable"]))
        except (KeyError, TypeError, ValueError):
            continue
        records[key] = record
    return records


# Select exact fit diagnostics or compute the diagonal fallback
def _observable_loss_statistics(
    *, likelihood_record: dict | None, h_mc, h_data, fallback_valid
) -> tuple[float | None, float | None, float | None, int]:
    if likelihood_record is None:
        objective, ndf = icecost.chi2_cost(h_mc=h_mc, h_data=h_data)
        reduced = objective / float(ndf) if ndf > 0 else None
        return float(objective), float(ndf), reduced, int(bool(fallback_valid))

    if not bool(likelihood_record.get("valid", False)):
        return None, None, None, 0
    try:
        objective = float(likelihood_record["objective_value"])
        ndf = float(likelihood_record["ndf"])
    except (KeyError, TypeError, ValueError):
        return None, None, None, 0
    if not np.isfinite(objective) or not np.isfinite(ndf) or ndf < 0.0:
        return None, None, None, 0
    reduced = objective / ndf if ndf > 0.0 else None
    return objective, ndf, reduced, 1


# Format one fit diagnostic title together with the measurement cuts
def _loss_plot_title(
    *, set_title: str, tunename: str, objective: float | None, ndf: float | None, reduced: float | None
) -> str:
    quantity = "$\\chi^2$"
    if objective is None or ndf is None:
        diagnostic = f"{quantity} / ndf = undefined [{tunename}]"
    else:
        diagnostic = f"{quantity} / ndf = {objective:0.4g} / {ndf:0.4g} = " + (
            f"{reduced:0.4g} [{tunename}]" if reduced is not None else f"undefined [{tunename}]"
        )
    return "\n".join(part for part in (set_title, diagnostic) if part)


# Visualize MC and data histograms with their exact fit diagnostics
def visualize_losses(
    *, results: dict, run_name: str, cdir: str, tunename: str, summary: dict, valid_arr: list,
    likelihood: dict | None = None, output_dir: str | None = None, summary_file: str | None = None,
) -> dict:
    if output_dir is None:
        output_dir = icetune.icetune_figure_dir(cdir=cdir, run_name=run_name)

    if summary_file is None:
        summary_file = os.path.join(output_dir, "summary.json")

    datasets, mc, data = (results[key] for key in ("datasets", "mc", "data"))
    likelihood_records = _likelihood_observable_map(likelihood)

    statistics = {name: {} for name in ("chi2", "nbin", "chi2red", "valid")}

    for i in range(len(mc)):
        for j in range(len(mc[i])):
            dataset_set = datasets[i]["sets"][j]
            label = readers.clean_filename(f"{datasets[i]['type']}__{dataset_set['name']}")
            set_title = str(dataset_set.get("title") or "").strip()

            subset_statistics = {name: {} for name in statistics}

            hist_keys = icecost.get_histogram_keys(
                mc_subset=mc[i][j], data_subset=data[i][j], context=f" for dataset={i}, subset={j}, label={label}"
            )

            for obs_key in hist_keys:
                record_key = (i, j, str(obs_key))
                likelihood_record = likelihood_records.get(record_key)
                chi2_val, ndf_val, chi2red_val, valid_val = _observable_loss_statistics(
                    likelihood_record=likelihood_record,
                    h_mc=mc[i][j][obs_key]["hdata"],
                    h_data=data[i][j][obs_key]["hdata"],
                    fallback_valid=valid_arr[i][j][obs_key],
                )
                title = _loss_plot_title(
                    set_title=set_title, tunename=tunename, objective=chi2_val, ndf=ndf_val, reduced=chi2red_val
                )

                for yscale in ["linear", "log"]:
                    fullpath = os.path.join(output_dir, label, yscale)
                    ensure_dir(fullpath)

                    filename = os.path.join(fullpath, f"hplot__{obs_key}.pdf")

                    fig, ax = plot.superplot([data[i][j][obs_key], mc[i][j][obs_key]], yscale=yscale)
                    ax[0].set_title(title, fontsize=8)
                    fig.savefig(filename, bbox_inches="tight")
                    plt.close(fig)

                values = (chi2_val, ndf_val, chi2red_val, valid_val)
                for name, value in zip(statistics, values, strict=True):
                    subset_statistics[name][obs_key] = value

            for name in statistics:
                statistics[name][label] = subset_statistics[name]

    payload = icetune.figure_summary_timestamped({**summary, **statistics})

    ensure_dir(os.path.dirname(summary_file))
    write_json_file(summary_file, payload, indent=4)

    return payload
