# Fit the additive Pb photoabsorption continuum to HEPData
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import argparse
import json
import sys
from pathlib import Path

import json5
import numpy as np
from core.io.serialize import load_json_file
from numpy.polynomial.legendre import leggauss
from scipy.optimize import least_squares

# Resolve the icepack package for direct execution from the repository
if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from icepack._common.hepdata import error_pairs, read_table, select_rows, table_bins, table_values


# Read Pb cross sections and quoted errors, retaining finite energy-bin averages
def measurements(photo: dict, nodes: int):
    """Read original photonuclear HEPData in MeV and mb"""
    energy, weight, value, error = [], [], [], []
    rule, weights = leggauss(nodes)
    sources = (("405665", 1, 5, photo["A"] / 1000.0), ("83727", 3, 2, 1.0 / 1000.0))
    for record, table_id, group, scale in sources:
        table = read_table(f"HEPData/PHOTONUCLEAR/HEPData-ins{record}-v1-json/Table{table_id}.json", dimensions=1)
        values = table_values(table, group)
        selected = select_rows(table, [i for i, item in enumerate(values) if item["value"] != "-"])
        values = table_values(selected, group)
        centers, edges = table_bins(selected, reference=table)
        for row, center, pair in zip(selected["values"], centers, edges, strict=True):
            energy.append((np.mean(pair) + rule * np.diff(pair)[0] / 2.0) * 1000.0
                          if "low" in row["x"][0] else np.full(nodes, center * 1000.0))
            weight.append(weights / 2.0)
        value.extend(float(item["value"]) * scale for item in values)
        # Original tables supply one symmetric error, with no covariance matrix
        error.extend(np.abs(error_pairs(values, 0)[:, 0]) * scale)
    return tuple(np.asarray(item) for item in (energy, weight, value, error))


# Evaluate the analytic response in the measured energy range
def response(energy: np.ndarray, fit: np.ndarray, photo: dict):
    """Compute photoabsorption [mb] above the QD and resonance boundaries [MeV]"""
    norm, width, threshold, match = fit
    cont, qd = photo["continuum"], photo["quasi_deuteron"]
    a, z = photo["A"], photo["Z"]
    gdr = next(item for item in photo["gdr"]["isotopes"] if (item["A"], item["Z"]) == (a, z))
    res = photo["resonance"]
    peak = 2.0 * gdr["strength"] * photo["gdr"]["systematics"]["trk"] * z * (a - z) / (np.pi * a * gdr["width"])
    sigma = peak * gdr["width"]**2 * energy**2 / ((energy**2 - gdr["energy"]**2)**2 + gdr["width"]**2 * energy**2)
    sigma += (qd["levinger"] * (a - z) * z / a * qd["deuteron"]["norm"]
              * (energy - qd["deuteron"]["threshold"])**1.5 / energy**3 * np.exp(-qd["pauli"]["high_exp"] / energy))
    sigma += res["area"] / (res["width"] * np.sqrt(2.0 * np.pi)) * np.exp(-0.5 * ((energy - res["energy"]) / res["width"])**2)
    u = np.clip((energy - threshold) / (match - threshold), 0.0, 1.0)
    excess = energy - cont["mean"]
    return sigma + u**2 * (3.0 - 2.0 * u) * (cont["constant"] + cont["log2"] * np.log(energy / cont["omega0"])**2
                                            + norm * excess * np.exp(-excess / width))


# Fit four positive continuum controls while retaining the fixed low-energy response
def main():
    """Print fitted parameters without editing the source data or steering cards"""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tune", type=Path, default=Path("modeldata/TUNE0"))
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    emd = load_json_file(args.tune / "GENERAL.json", loader=json5.load)["PARAM_NUCLEAR"]["EMD"]
    isotope = next(item for item in emd["isotopes"] if (item["A"], item["Z"]) == (208, 82))
    photo = dict(isotope["photoabsorption"], A=isotope["A"], Z=isotope["Z"], gdr=emd["gdr"])
    with (args.tune / "NUMERICS.json").open() as stream:
        nodes = json5.load(stream)["NUMERICS_NUCLEAR"]["EMD"]["response"]["nodes"]
    energy, weight, value, error = measurements(photo, nodes)
    if np.min(energy) < max(photo["quasi_deuteron"]["pauli"]["high_edge"], photo["resonance"]["threshold"]):
        raise ValueError("Fit data must lie above the QD and resonance boundaries")
    cont = photo["continuum"]
    floor = max(cont["mean"], photo["resonance"]["threshold"])
    initial = [cont["norm"], cont["width"], cont["threshold"] - floor, cont["match"] - cont["threshold"]]
    if not np.isfinite(initial).all() or np.min(initial) <= 0.0:
        raise ValueError("Continuum controls require positive scales and an ordered onset")

    # Positive widths and an ordered onset avoid inadmissible trial responses
    def unpack(log_fit):
        norm, width, start, span = np.exp(log_fit)
        return np.array([norm, width, floor + start, floor + start + span])

    # Minimize residuals in the measured points and the published finite bins
    def residual(log_fit):
        prediction = np.sum(response(energy, unpack(log_fit), photo) * weight, axis=1)
        return (prediction - value) / error

    fit = least_squares(residual, np.log(initial))
    if not fit.success:
        raise RuntimeError(fit.message)
    result = dict(zip(("norm", "width", "threshold", "match"), unpack(fit.x), strict=True))
    result.update(chi2=float(fit.fun @ fit.fun), ndof=len(value) - len(fit.x),
                  errors="Quoted HEPData errors treated as independent (no covariance supplied)")
    text = json.dumps(result, indent=2)
    print(text)
    if args.output:
        args.output.write_text(text + "\n")


if __name__ == "__main__":
    main()
