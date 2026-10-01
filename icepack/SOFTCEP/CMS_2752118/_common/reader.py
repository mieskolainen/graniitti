# HEPData 2752118 reader
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import numpy as np
import pandas as pd
from core.io.cache import cache
from core.stats.hist import bins2binwidth, center2edgebins, edge2centerbins
from core.stats.uncertainty import finalize_uncertainties, uncertainty_source
from termcolor import cprint

from icepack._common.hepdata import error_pairs, read_table, rebin_common, table_values


# Report angular intervals with no detector acceptance
def acceptance_dphi_2752118(k, i, edges, cbin_2, mask, reference_cbin_2):
    """Print zero-acceptance regions for debugging"""
    if np.sum(~mask) > 0:
        delta = (reference_cbin_2[1] - reference_cbin_2[0]) / 2
        gap_indices = np.where(~mask)[0]
        segments = np.split(gap_indices, np.where(np.diff(gap_indices) != 1)[0] + 1)

        for seg in segments:
            mmin = cbin_2[seg[0]] - delta
            mmax = cbin_2[seg[-1]] + delta

            print(
                f"cut[{k}] = not (rcut(pt1, ({edges[i, 0][0]}, {edges[i, 0][1]})) and "
                f"rcut(pt2, ({edges[i, 1][0]}, {edges[i, 1][1]})) and "
                f"rcut(dphi, ({mmin}, {mmax}) ))"
            )
            k += 1

            if edges[i, 0][0] != edges[i, 1][0]:
                print(
                    f"cut[{k}] = not (rcut(pt1, ({edges[i, 1][0]}, {edges[i, 1][1]})) and "
                    f"rcut(pt2, ({edges[i, 0][0]}, {edges[i, 0][1]})) and "
                    f"rcut(dphi, ({mmin}, {mmax}) ))"
                )
                k += 1

    return k


# Compute the transverse momentum interval for one category
def category_edges_2752118(
    center: float, reference_centers: np.ndarray, reference_edges: np.ndarray
) -> np.ndarray:
    """Return category bin edges for one CMS center value"""
    index = reference_indices_2752118(
        input_array=np.asarray([center]), reference_array=reference_centers
    )[0]
    return np.asarray([reference_edges[index], reference_edges[index + 1]])


# Match rounded published coordinates to the common bin grid
def reference_indices_2752118(
    input_array: np.ndarray, reference_array: np.ndarray, tol: float = 1e-5
) -> np.ndarray:
    """Return unique reference-grid indices for one CMS category center array"""
    indices = []
    for value in np.asarray(input_array):
        matches = np.where(np.isclose(reference_array, value, atol=tol))[0]
        if len(matches) != 1:
            raise ValueError(
                f"ReadHEPData_2752118: Could not map center {value} to a unique reference bin"
            )
        indices.append(matches[0])
    return np.asarray(indices, dtype=np.int64)


# Compute the physical azimuthal grid from the published centers
def dphi_grid_2752118(published_centers: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Validate the rounded HEPData centers and return exact angular bins"""
    edges = np.linspace(0.0, np.pi, len(published_centers) + 1)
    centers = edge2centerbins(edges)
    if not len(published_centers) or not np.allclose(
        published_centers, centers, atol=5e-5, rtol=0.0
    ):
        raise ValueError("ReadHEPData_2752118: Expected uniformly spaced dPhi centers between zero and pi")
    return centers, edges


# Read the published three dimensional cross section and signed error sources
@cache
def read(filename: str, file_filter: str = None, rebin_factor: int = None) -> dict:
    """Read CMS HEPData 2752118 JSON content"""
    table = read_table(filename)
    if table["x_count"] != 3:
        raise ValueError("ReadHEPData_2752118: expected two proton transverse momenta and one central observable")
    coordinates = np.asarray([[axis["value"] for axis in row["x"]] for row in table["values"]], dtype=float)
    labels = [header["name"] for header in table["headers"][:3]]
    frame = pd.DataFrame(coordinates, columns=labels)
    if file_filter is not None:
        frame = frame.query(file_filter, engine="python")
    selected = table_values(table)
    values = [selected[index] for index in frame.index]
    X = np.column_stack((frame.to_numpy(), [float(value["value"]) for value in values],
                         error_pairs(values, "stat"), error_pairs(values, "syst (norm.)"),
                         error_pairs(values, "syst (effic.)")))
    ds = {
        "x": {},
        "bins": {},
        "binwidth": {},
        "xlim": {},
        "y": {},
        "y_err_stat": {},
        "y_err_syst": {},
        "y_err_sexp": {},
        "y_err": {},
        "valid": {},
        "uncertainties": {},
    }

    x2_label = labels[2].lower()
    use_common_x2_grid = ("\u03c6" in x2_label) or ("phi" in x2_label)

    comb = np.unique(np.round(X[:, [0, 1]], decimals=4), axis=0)
    print(__name__ + f".read: Number of (x0, x1) pair hyperbins (categories) found = {len(comb)}")

    # Infer physical bin edges from the complete table before selecting rows
    reference_cbin_0 = np.sort(np.unique(coordinates[:, 0]))
    reference_cbin_1 = np.sort(np.unique(coordinates[:, 1]))
    reference_cbin_2 = np.sort(np.unique(coordinates[:, 2]))
    reference_edges_0 = center2edgebins(reference_cbin_0)
    reference_edges_1 = center2edgebins(reference_cbin_1)
    reference_edges_2 = center2edgebins(reference_cbin_2)
    if use_common_x2_grid:
        reference_cbin_2, reference_edges_2 = dphi_grid_2752118(reference_cbin_2)

    for i, c in enumerate(comb):
        IND = (np.abs(np.round(X[:, 0], decimals=4) - c[0]) < 1e-4) & (
            np.abs(np.round(X[:, 1], decimals=4) - c[1]) < 1e-4
        )
        X_IND = X[IND]
        x0_edges = category_edges_2752118(
            center=c[0], reference_centers=reference_cbin_0, reference_edges=reference_edges_0
        )
        x1_edges = category_edges_2752118(
            center=c[1], reference_centers=reference_cbin_1, reference_edges=reference_edges_1
        )
        row_indices = reference_indices_2752118(
            input_array=X_IND[:, 2],
            reference_array=reference_cbin_2,
            tol=5e-5 if use_common_x2_grid else 1e-5,
        )
        if use_common_x2_grid:
            grid_indices = np.arange(len(reference_cbin_2), dtype=np.int64)
        else:
            grid_indices = np.arange(np.min(row_indices), np.max(row_indices) + 1, dtype=np.int64)

        cbin_2 = reference_cbin_2[grid_indices]
        edges_2 = reference_edges_2[grid_indices[0] : grid_indices[-1] + 2]
        local_indices = row_indices - grid_indices[0]
        # Missing measurement rows stay zero padded with mask false and never represent measured zeros
        # [REFERENCE: CMS and TOTEM, Phys. Rev. D 109 (2024) 112013, Sec. VI B]
        mask = np.zeros(len(cbin_2), dtype=np.bool_)
        mask[local_indices] = True

        ds["x"][i] = {0: c[0], 1: c[1], 2: cbin_2}
        ds["bins"][i] = {0: x0_edges, 1: x1_edges, 2: edges_2}
        ds["binwidth"][i] = {
            0: bins2binwidth(x0_edges),
            1: bins2binwidth(x1_edges),
            2: bins2binwidth(edges_2),
        }
        ds["xlim"][i] = [np.min(edges_2), np.max(edges_2)]
        ds["y"][i] = np.zeros_like(mask, dtype=np.float64)
        ds["y_err_stat"][i] = np.zeros_like(mask, dtype=np.float64)
        ds["y_err_syst"][i] = np.zeros_like(mask, dtype=np.float64)
        ds["y_err_sexp"][i] = np.zeros_like(mask, dtype=np.float64)
        ds["valid"][i] = mask

        ds["y"][i][local_indices] = X_IND[:, 3]
        ds["y_err_stat"][i][local_indices] = (np.abs(X_IND[:, 4]) + np.abs(X_IND[:, 5])) / 2
        ds["y_err_syst"][i][local_indices] = (np.abs(X_IND[:, 6]) + np.abs(X_IND[:, 7])) / 2
        ds["y_err_sexp"][i][local_indices] = (np.abs(X_IND[:, 8]) + np.abs(X_IND[:, 9])) / 2
        normalization_shift = np.zeros_like(mask, dtype=np.float64)
        efficiency_shift = np.zeros_like(mask, dtype=np.float64)
        normalization_shift[local_indices] = (X_IND[:, 6] - X_IND[:, 7]) / 2
        efficiency_shift[local_indices] = (X_IND[:, 8] - X_IND[:, 9]) / 2
        # [REFERENCE: HEPData 2752118 Figures 19-21 named uncertainty columns]
        ds["uncertainties"][i] = [
            uncertainty_source(
                "statistical",
                ds["y_err_stat"][i],
                category="statistical",
                correlation="uncorrelated",
                provenance="HEPData stat column",
            ),
            uncertainty_source(
                "normalization",
                ds["y_err_syst"][i],
                shift=normalization_shift,
                category="systematic", correlation="collective", effect="multiplicative",
                scope="CMS_2752118:normalization",
                provenance="HEPData syst (norm.) column",
            ),
            uncertainty_source(
                "efficiency",
                ds["y_err_sexp"][i],
                shift=efficiency_shift,
                category="systematic", correlation="collective", effect="multiplicative",
                scope="CMS_2752118:efficiency",
                provenance="HEPData syst (effic.) column",
            ),
        ]

    if rebin_factor is not None:
        cprint(__name__ + f".read: Rebin histograms with factor // {rebin_factor}", "red")
    for i in range(len(comb)):
        local = {"x": ds["x"][i][2], "bins": ds["bins"][i][2], "binwidth": ds["binwidth"][i][2],
                 "y": ds["y"][i], "valid": ds["valid"][i]}
        finalize_uncertainties(local, ds["uncertainties"][i])
        if rebin_factor is not None:
            local = rebin_common(local, rebin_factor)
        for key in ("x", "bins", "binwidth"):
            ds[key][i][2] = local[key]
        for key in ("y", "y_err", "y_err_stat", "y_err_syst", "valid", "uncertainties"):
            ds[key][i] = local[key]
        efficiency = next(source for source in local["uncertainties"] if source["name"] == "efficiency")
        ds["y_err_sexp"][i] = np.abs(np.asarray(efficiency["shift"], dtype=float))

    return ds
