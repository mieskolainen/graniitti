# LHCb HEPData 1373746 combined-sample exclusive Upsilon reader
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from pathlib import Path

import numpy as np
from core.io.cache import cache
from core.stats.hist import doublet2binwidth, doublet2linear
from core.stats.uncertainty import finalize_uncertainties, uncertainty_source

from icepack._common.hepdata import error_pairs, read_table, table_bins, table_values

REGION_FILES = {
    "fiducial": "Table1.json",
    "central_system": "Table2.json",
}


# Resolve and validate the single published Upsilon state
def resolve_state(state, dataset):
    selected = state
    if selected is None:
        if dataset is None or "state" not in dataset:
            raise ValueError("ReadHEPData_1373746: dataset state is required")
        selected = dataset["state"]
    selected = str(selected)
    if selected != "Upsilon(1S)":
        raise ValueError(f"ReadHEPData_1373746: unsupported state {selected!r}")
    return selected


# Resolve and validate one explicitly selected measurement region
def resolve_region(region, dataset):
    selected = region
    if selected is None:
        if dataset is None or "region" not in dataset:
            raise ValueError("ReadHEPData_1373746: dataset region is required")
        selected = dataset["region"]
    selected = str(selected)
    if selected not in REGION_FILES:
        raise ValueError(f"ReadHEPData_1373746: unsupported region {selected!r}")
    return selected


# Require the HEPData table associated with the selected region
def validate_region_file(filename, region):
    required = REGION_FILES[region]
    if Path(filename).name != required:
        raise ValueError(f"ReadHEPData_1373746: region {region!r} requires {required}")


# Read the published Upsilon table and convert fiducial bin integrals to rapidity densities
@cache
def _read(filename, region, rebin_factor):
    if rebin_factor is not None:
        raise ValueError("ReadHEPData_1373746: rebinning is not supported")
    table = read_table(filename)
    values = table_values(table)
    x, edges = table_bins(table)
    bins = doublet2linear(edges)
    widths = doublet2binwidth(edges)
    scale = 1.0 / widths if region == "fiducial" else np.ones_like(widths)
    labels = [("statistical", "stat"), ("systematic", "sys")] if region == "fiducial" else [("combined", "error")]
    sources = []
    for category, label in labels:
        pair = error_pairs(values, label) * scale[:, None]
        sources.append(uncertainty_source(
            category, pair[:, 0], down=pair[:, 1], category=category, correlation="uncorrelated",
            provenance=f"Original HEPData {label} errors, bin covariance unavailable",
        ))
    return finalize_uncertainties({
        "state": "Upsilon(1S)", "region": region, "x": x, "binedges": edges, "bins": bins,
        "binwidth": widths, "xlim": bins[[0, -1]],
        "y": np.asarray([value["value"] for value in values], dtype=float) * scale,
    }, sources)


# Read one explicitly selected fiducial or central-system measurement
def read(filename, rebin_factor=None, state=None, region=None, dataset=None):
    resolve_state(state, dataset)
    selected_region = resolve_region(region, dataset)
    validate_region_file(filename, selected_region)
    return _read(filename, selected_region, rebin_factor)
