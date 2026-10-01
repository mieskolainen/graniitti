# LHCb HEPData 2825384 exclusive charmonium reader
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from pathlib import Path

import numpy as np
from core.io.cache import cache
from core.stats.hist import doublet2binwidth, doublet2linear
from core.stats.uncertainty import (
    covariance_projection_error,
    finalize_uncertainties,
    norm,
    source_covariance,
    uncertainty_source,
)

from icepack._common.hepdata import error_pairs, read_table, table_bins, table_values

STATE_FILES = {
    "J/psi": "Table4.json",
    "psi(2S)": "Table5.json",
}

REGIONS = {"central_system"}

ERROR_LABELS = {
    "stat": "(Stat.)",
    "uncorrelated_syst": "(Uncorr.)",
    "correlated_syst": "(Corr.)",
    "luminosity": "(Lumi.)",
}


# Resolve the state explicitly or from one iceplot dataset set
def resolve_state(state, dataset):
    if state is not None:
        return state
    if dataset is None or "state" not in dataset:
        raise ValueError("ReadHEPData_2825384: dataset state is required")
    return str(dataset["state"])


# Resolve and validate one explicitly selected measurement region
def resolve_region(region, dataset):
    selected = region
    if selected is None:
        if dataset is None or "region" not in dataset:
            raise ValueError("ReadHEPData_2825384: dataset region is required")
        selected = dataset["region"]
    selected = str(selected)
    if selected not in REGIONS:
        raise ValueError(f"ReadHEPData_2825384: unsupported region {selected!r}")
    return selected


# Validate that the selected HEPData table belongs to the requested state
def validate_state_file(filename, state):
    if state not in STATE_FILES:
        raise ValueError(f"ReadHEPData_2825384: unsupported state {state!r}")
    if Path(filename).name != STATE_FILES[state]:
        raise ValueError(f"ReadHEPData_2825384: state {state!r} requires {STATE_FILES[state]}")


# Integrate one rapidity table with its published correlation structure
def integrate_rapidity(data):
    width = np.asarray(data["binwidth"], dtype=float)
    errors = {name: float(covariance_projection_error(source_covariance(source), width))
              for name, source in zip(ERROR_LABELS, data["uncertainties"], strict=True)}
    return {"value": float(np.dot(data["y"], width)), **errors,
            "total": float(norm(np.asarray(list(errors.values())))), "unit": "nb"}


# Read one resolved central-system rapidity distribution
@cache
def _read_central_system(filename, selected_state, rebin_factor):
    if rebin_factor is not None:
        raise ValueError("ReadHEPData_2825384: rebinning is not supported")

    validate_state_file(filename, selected_state)
    table = read_table(filename)
    values = table_values(table)
    x, binedges = table_bins(table)
    errors = {key: error_pairs(values, label) for key, label in ERROR_LABELS.items()}
    components = {key: np.mean(np.abs(pair), axis=1) for key, pair in errors.items()}
    bins = doublet2linear(binedges)
    ds = {
        "state": selected_state,
        "region": "central_system",
        "x": x,
        "binedges": binedges,
        "bins": bins,
        "binwidth": doublet2binwidth(binedges),
        "xlim": np.array([bins[0], bins[-1]]),
        "y": np.asarray([value["value"] for value in values], dtype=float),
        "error_components": components,
    }
    # [REFERENCE: HEPData 2825384 Tables 4 and 5 named uncertainty columns]
    finalize_uncertainties(
        ds,
        [
            uncertainty_source(
                "statistical",
                errors["stat"][:, 0],
                down=errors["stat"][:, 1],
                category="statistical",
                correlation="uncorrelated",
                provenance="HEPData stat column",
            ),
            uncertainty_source(
                "uncorrelated_systematic",
                errors["uncorrelated_syst"][:, 0],
                down=errors["uncorrelated_syst"][:, 1],
                category="systematic",
                correlation="uncorrelated",
                provenance="HEPData uncorrelated syst column",
            ),
            uncertainty_source(
                "correlated_systematic",
                errors["correlated_syst"][:, 0],
                down=errors["correlated_syst"][:, 1],
                category="systematic",
                correlation="collective",
                scope=f"LHCb_2825384:{selected_state}:correlated_systematic",
                provenance="HEPData correlated syst column",
            ),
            uncertainty_source(
                "luminosity",
                errors["luminosity"][:, 0],
                down=errors["luminosity"][:, 1],
                category="systematic",
                correlation="collective",
                scope="LHCb_2825384:luminosity",
                effect="multiplicative",
                provenance="HEPData luminosity column",
            ),
        ],
    )
    return ds


# Read the selected acceptance corrected HEPData measurement
def read(filename, rebin_factor=None, state=None, region=None, dataset=None):
    selected_state = resolve_state(state, dataset)
    resolve_region(region, dataset)
    return _read_central_system(filename, selected_state, rebin_factor)
