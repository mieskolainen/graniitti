# H1 and ZEUS exclusive photoproduction HEPData reader
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import re
from pathlib import Path

import numpy as np
from core.io.cache import cache
from core.stats.hist import bins2binwidth, doublet2binwidth, doublet2linear, edge2centerbins
from core.stats.uncertainty import finalize_uncertainties, uncertainty_source

from icepack._common.hepdata import error_pairs, normalized, number, read_table, resolve_group, table_bins, table_values


# Read positive photon fluxes from the original column or common qualifier
def photon_flux(table, rows=None):
    labels = {"PHI_T", r"$\Phi_{\gamma/e}$"}
    columns = [i for i, header in enumerate(table["headers"][:table["x_count"]]) if header["name"] in labels]
    if len(columns) == 1:
        rows = table["values"] if rows is None else rows
        flux = np.asarray([number(row["x"][columns[0]]["value"]) for row in rows])
    elif not columns:
        values = [entry["value"] for name, entries in table["qualifiers"].items() if name in labels for entry in entries]
        if len(values) != 1:
            raise ValueError("H1 photon flux requires one original column or qualifier")
        flux = np.asarray(number(values[0]))
    else:
        raise ValueError("H1 photon flux has ambiguous columns")
    if not flux.size or not np.isfinite(flux).all() or np.any(flux <= 0.0):
        raise ValueError("H1 photon flux must be finite and positive")
    return flux


# Read published intervals or infer edges using the common HEPData bin conversion
def read_bins(table, lower=None):
    coordinates = [row["x"][0] for row in table["values"]]
    centers, binedges = table_bins(table)
    if lower is not None and all("low" not in item for item in coordinates):
        binedges[0, 0] = lower
    header = table.get("headers", [{}])[0].get("name", "").upper()
    if "PHI" in header:
        # Older HERA tables label azimuth coordinates in degrees but quote densities per radian
        if "DEGREES" in header:
            centers, binedges = np.deg2rad(centers), np.deg2rad(binedges)
        elif all("value" in item and "low" not in item for item in coordinates):
            # Recover periodic bins allowing accumulated rounding of the published bin steps
            edges = np.linspace(0.0, 2.0 * np.pi, len(centers) + 1)
            precision = np.asarray([0.5 * 10.0 ** (-len(str(item["value"]).partition(".")[2]))
                                    for item in coordinates])
            if np.all(np.abs(centers - edge2centerbins(edges)) <= np.cumsum(precision)):
                binedges = np.column_stack((edges[:-1], edges[1:]))
    bins = doublet2linear(binedges)
    if len(bins) != len(binedges) + 1:
        raise ValueError("ReadHEPData_HERA: disjoint bins are not supported")
    return centers, binedges, bins


# Build uncertainty sources from the published marginal errors
def read_uncertainties(values, scope):
    count = len(values[0].get("errors", []))
    if count == 0 or any(len(value.get("errors", [])) != count for value in values):
        raise ValueError("ReadHEPData_HERA: inconsistent uncertainty fields")
    sources = []
    for index in range(count):
        first = values[0]["errors"][index]
        label = str(first.get("label", "total" if count == 1 else f"error_{index}"))
        key = normalized(label)
        category = "statistical" if "stat" in key and "syst" not in key else "combined"
        pairs = np.abs(error_pairs(values, index))
        sources.append(
            uncertainty_source(
                key or f"error_{index}",
                pairs[:, 0],
                down=pairs[:, 1],
                category=category,
                correlation="uncorrelated",
                scope=f"{scope}:{key or index}",
                provenance=f"Official HEPData {label} field, diagonal covariance assumed because bin correlations are unavailable",
            )
        )
    return sources


# Read one selected HERA dependent variable
@cache
def _read(filename, file_filter, flux, lower):
    table = read_table(filename, dimensions=1)
    group = resolve_group(table, file_filter)
    values = table_values(table, group)
    x, binedges, bins = read_bins(table, lower)
    dataset = {
        "x": x,
        "binedges": binedges,
        "bins": bins,
        "binwidth": doublet2binwidth(binedges),
        "xlim": np.asarray([bins[0], bins[-1]], dtype=float),
        "y": np.asarray([number(value["value"]) for value in values], dtype=float),
        "group": group,
        "headers": table["headers"],
        "qualifiers": table["qualifiers"],
    }
    if flux is not None:
        dataset["mc_scale"] = 1.0 / np.sum(photon_flux(read_table(Path(filename).parent / flux)))
    scope = f"{Path(filename).parent.name}:group{group}"
    sources = read_uncertainties(values, scope)
    finalize_uncertainties(dataset, sources)
    return dataset


# Read one official HERA table with an optional dependent-variable filter
def read(filename, file_filter=None, rebin_factor=None, hist=None, **_kwargs):
    if rebin_factor is not None:
        raise ValueError("ReadHEPData_HERA: rebinning is not supported")
    hist = hist or {}
    return _read(filename, file_filter, hist.get("flux"), hist.get("bin_start"))


# Transpose the H1 differential W table into absolute transfer spectra at its published W points
# [REFERENCE: H1 arXiv:hep-ex/0510016, Tables 6 and 7]
def read_wt(filename, transfer_table):
    table = read_table(filename)
    transfer = read_table(transfer_table)
    reference = read_bins(transfer)[2]
    limits = np.asarray([[number(value) for value in re.findall(r"[\d.]+", q["value"])[:2]]
                         for q in table["qualifiers"]["ABS(T)"]])
    # Table 7 has contiguous bins, the second lower edge is mistyped in HEPData
    edges = np.r_[limits[0, 0], limits[:, 1]]
    if not np.all(np.isclose(edges[:, None], reference).any(axis=1)) or np.any(np.diff(edges) <= 0):
        raise ValueError("H1 W,t boundaries disagree with the independent transfer table")
    notes = []
    if not np.allclose(limits[:, 0], edges[:-1]):
        notes.append("Contiguous t bins from paper Table 7 and HEPData Table7.json resolve the Table11.json lower-edge typo")
    spectra = []
    for row in table["values"]:
        energy = number(row["x"][0]["value"])
        values = row["y"]
        if len(values) != len(edges) - 1:
            raise ValueError("H1 W,t row does not contain every transfer bin")
        data = {"x": edge2centerbins(edges), "bins": edges,
                "binedges": np.column_stack((edges[:-1], edges[1:])), "binwidth": bins2binwidth(edges),
                "y": np.asarray([number(value["value"]) for value in values]), "W": energy,
                "source": str(filename), "bin_notes": notes}
        finalize_uncertainties(data, read_uncertainties(values, f"{filename}:W={energy}"))
        spectra.append(data)
    return spectra


# Read elastic rho W,t bins with the published statistical covariance and signed offsets
# [REFERENCE: H1 arXiv:2005.14471, Section 4.5 and Tables 12-13]
def read_rho_wt(filename, covariance):
    from icepack.PHOTOPROD._common.diss_reader import rho_stat

    table = read_table(filename)
    rows = [row for row in table["values"] if "value" in row["x"][2]]
    values = np.asarray([number(row["y"][0]["value"]) for row in rows])
    selected = [row["y"][0] for row in rows]
    offsets = np.stack([error_pairs(selected, index) for index in range(1, len(selected[0]["errors"]))], axis=1)
    return {"y": values, "binedges": np.asarray([[number(row["x"][0][key]) for key in ("low", "high")] for row in rows]),
            "W": np.asarray([np.sqrt(number(row["x"][1]["low"]) * number(row["x"][1]["high"])) for row in rows]),
            "stat_cov": rho_stat(rows, covariance, bin_index=3)["covariance"],
            "offsets": offsets, "source": str(filename)}
