# HERA differential fits and covariance diagnostics
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import math
import re
from functools import lru_cache

import numpy as np
from numba import njit
from scipy.linalg import block_diag, solve_triangular
from scipy.optimize import curve_fit, least_squares

from icepack._common import hepdata
from icepack.PHOTOPROD._common import hera_reader

from . import inputs, model


# Profile multiplicative systematic shifts without the data-dependent normalization bias
# The nuisance penalties accompany the whitened statistical residual in the least-squares fit
def profiled_residual(predicted, observed, chol, relative):
    residual = solve_triangular(chol, predicted - observed, lower=True)
    shifts = solve_triangular(chol, predicted[:, None] * relative, lower=True)
    nuisance = np.linalg.solve(np.eye(shifts.shape[1]) + shifts.T @ shifts, -shifts.T @ residual)
    return np.concatenate((residual + shifts @ nuisance, nuisance))


# Require HEPData JSON covariance before enabling this fit
def h1_spectra():
    raise ValueError('H1 covariance fit unavailable: a supported HEPData JSON covariance input is required')


# Read both ZEUS decay channels with a shared normalization uncertainty per channel
# [REFERENCE: ZEUS arXiv:hep-ex/0201043, Sections 7, 9 and 10]
def zeus_spectra():
    folder = model.REPOSITORY_ROOT / "HEPData/PHOTOPROD/HEPData-ins582237-v1-json"
    nodes, weights = np.polynomial.legendre.leggauss(32)
    spectra = []
    for channel, filename, normalization in (("muon", "Table9.json", "Table1.json"),
                                              ("electron", "Table10.json", "Table2.json")):
        path = folder / filename
        table, integrated = hepdata.read_table(path), hepdata.read_table(folder / normalization)
        for qualifier in table["qualifiers"]["W"]:
            selection = qualifier["value"]
            low, high = (float(value) for value in re.findall(r"[\d.]+", selection)[:2])
            energies = np.exp((math.log(high) + math.log(low) + nodes * math.log(high / low)) / 2.0)
            data = dict(hera_reader.read(path, file_filter=selection, hist={"bin_start": 0.0}))
            rows = [row for row in integrated["values"]
                    if hepdata.number(row["x"][0]["low"]) >= low and hepdata.number(row["x"][0]["high"]) <= high]
            if not rows or not np.allclose([hepdata.number(rows[0]["x"][0]["low"]), hepdata.number(rows[-1]["x"][0]["high"])], [low, high]):
                raise ValueError("ZEUS normalization bins do not cover the differential W range")
            # Combine the finer integrated W bins with the same logarithmic photon flux approximation
            weight = np.asarray([math.log(hepdata.number(row["x"][0]["high"]) / hepdata.number(row["x"][0]["low"])) for row in rows])
            norm = np.dot(weight, [max(hepdata.error_pair(row["y"][0]["errors"][1])) for row in rows])
            norm /= np.dot(weight, [hepdata.number(row["y"][0]["value"]) for row in rows])
            common = {"zeus_" + channel: np.full(len(data["y"]), norm)}
            # Figure 6 errors are statistical, the normalization correlation is an explicit approximation
            covariance = np.diag(data["y_err"]**2)
            data.update(stat_cov=covariance, source=str(path), experiment="ZEUS " + channel)
            data["fit_mask"] = np.ones(len(data["y"]), dtype=bool)
            if channel == "muon":
                data["fit_mask"][-1] = False
            shifts = data["y"] * norm
            spectra.append(("ZEUS " + channel + " " + selection, data, energies, weights / 2.0,
                            covariance + np.outer(shifts, shifts), common))
    return spectra


# Require HEPData JSON covariance before enabling this fit
def h1_wt_spectra(correlated=False):
    raise ValueError('H1 covariance fit unavailable: a supported HEPData JSON covariance input is required')


# Assemble independent differential measurements and retain all bins for validation
def jpsi_spectra(correlated=False):
    return h1_spectra() + zeus_spectra() + h1_wt_spectra(correlated)


# Require HEPData JSON covariance before enabling this fit
def jpsi_energy():
    raise ValueError('H1 covariance fit unavailable: a supported HEPData JSON covariance input is required')


# Assemble correlated residuals without replacing normalization information by shape fits
def jpsi_errors(observed, covariances, common):
    labels = sorted(set().union(*common))
    relative = np.concatenate([np.column_stack([sources.get(key, np.zeros(len(values))) for key in labels]) if labels else np.empty((len(values), 0))
                               for values, sources in zip(observed, common, strict=True)])
    return np.concatenate(observed), np.linalg.cholesky(block_diag(*covariances)), relative, labels


# Average the physical differential cross section over the measured transfer bins and photon flux
def jpsi_spectrum(parameters, proton: model.ProtonVertex, data, energy, flux, dlog=None):
    forward, slope, delta, alpha_prime = parameters
    if dlog is not None:
        flux = flux * model.dlog_factor(energy, 100.0, *dlog)
    nodes, weights = np.polynomial.legendre.leggauss(32)
    low, high = data["binedges"].T
    transfer = (low[:, None] + high[:, None] + (high - low)[:, None] * nodes) / 2.0
    evolution = np.sum(np.exp((delta - 4 * alpha_prime * transfer[:, :, None])
                              * np.log(energy / 100.0)) * flux, axis=-1)
    return forward * np.sum(model.proton_profile(transfer, proton)**2 * np.exp(-slope * transfer)
                            * evolution * weights, axis=1) / 2.0


# Fit absolute differential data with one physical coupling and constrained experimental nuisance shifts
# Integrated W projections and published exponential slopes are validation data, not independent fit inputs
# [REFERENCE: H1 arXiv:1304.5162, H1 arXiv:hep-ex/0510016 and ZEUS arXiv:hep-ex/0201043]
@lru_cache(maxsize=16)
def jpsi(proton: model.ProtonVertex, exclude="", correlated=False, dlog=None):
    spectra = [row for row in jpsi_spectra(correlated) if not exclude or not row[1]["experiment"].startswith(exclude)]
    observed, covariance, common = [], [], []
    for _, data, _, _, _, sources in spectra:
        use = data["fit_mask"]
        observed.append(data["y"][use])
        covariance.append(data["stat_cov"][np.ix_(use, use)])
        common.append({key: value[use] for key, value in sources.items()})
    measured, chol, relative, labels = jpsi_errors(observed, covariance, common)

    # Integrate the physical gamma-p cross section over each measured transfer bin
    def prediction(parameters):
        return np.concatenate([jpsi_spectrum(parameters, proton, data, w, flux, dlog)[data["fit_mask"]]
                               for _, data, w, flux, _, _ in spectra])

    # Profile only constrained experimental systematics, with no free spectrum normalizations
    def residual(parameters):
        return profiled_residual(prediction(parameters), measured, chol, relative)

    fit = least_squares(residual, [measured.max(), 1.0, 0.8, 0.1], bounds=(0.0, np.inf),
                        x_scale="jac", ftol=1e-9, xtol=1e-9, gtol=1e-9)
    if not fit.success or np.linalg.matrix_rank(fit.jac) != len(fit.x):
        raise ValueError("J/psi absolute W,t fit failed or is not identifiable: " + fit.message)
    parameters = fit.x
    full_cov = np.linalg.inv(fit.jac.T @ fit.jac)
    fitted_residual = residual(parameters)
    pulls = fitted_residual[len(measured):]
    report, offset = [], 0
    for name, data, w, flux, cov, _ in spectra:
        predicted = jpsi_spectrum(parameters, proton, data, w, flux, dlog)
        difference = predicted - data["y"]
        use = data["fit_mask"]
        count = int(use.sum())
        shifted = predicted[use] * (1 + relative[offset:offset + count] @ pulls)
        report.append({"spectrum": name, "experiment": data["experiment"], "source": data["source"],
                       "bin_notes": data.get("bin_notes", []), "binedges": data["binedges"].tolist(),
                       "data": data["y"].tolist(), "prediction": predicted.tolist(), "errors": np.sqrt(np.diag(cov)).tolist(),
                       "profiled_prediction": shifted.tolist(),
                       "chi2": float(difference @ np.linalg.solve(cov, difference)),
                       "fit_chi2": float(np.sum(fitted_residual[offset:offset + count]**2)),
                       "bins": len(predicted), "shape_fit_bins": use.tolist()})
        offset += count
    energy = jpsi_energy()
    values, chole, rele, elabels = jpsi_errors([row[2] for row in energy], [row[3] for row in energy],
                                             [row[4] for row in energy])
    forward, slope, delta, alpha_prime = parameters
    energy_prediction = np.concatenate([np.asarray([forward * (w / 100.0)**delta
        * model.profile_moments(proton, slope + 4 * alpha_prime * math.log(w / 100.0), tmax)[0]
        * (model.dlog_factor(w, 100.0, *dlog) if dlog is not None else 1.0) for w in wvalues]) for _, wvalues, _, _, _, tmax in energy])
    energy_residual = profiled_residual(energy_prediction, values, chole, rele)
    return {"parameters": parameters.tolist(), "dlog": dlog, "covariance": full_cov.tolist(),
            "chi2": float(2 * fit.cost), "ndf": len(measured) - len(parameters),
            "shape_chi2": float(2 * fit.cost), "energy_chi2": float(energy_residual @ energy_residual),
            "energy_used_in_fit": False, "spectra": report, "energy_data": values.tolist(),
            "energy_prediction": energy_prediction.tolist(), "excluded_experiment": exclude,
            "nuisance_pulls": dict(zip(labels, pulls.tolist(), strict=True)),
            "nuisance_chi2": float(pulls @ pulls),
            "energy_nuisance_pulls": dict(zip(elabels, energy_residual[len(values):].tolist(), strict=True)),
            "assumptions": ["One absolute gamma-p residue fits independent H1 2013, H1 2006 and ZEUS muon/electron W,t data",
                            "No free differential normalizations and no reuse of integrated projections in the fit",
                            "Fixed SOFT proton residue times exp(B_gammaPV*t/2), with the same MP/XP/GP transverse coupling",
                            "Transfer bin averages and H1 2013 published photon flux weights",
                            "ZEUS W averaging uses logarithmic quadrature, an approximation to the effective photon flux",
                            "H1 2006 uses the published W correction points",
                            "H1 2013 HE first two bins and ZEUS muon last bin excluded as in the papers, retained in reports",
                            "H1 2013 original covariance and signed multiplicative systematic sources",
                            'ZEUS t errors, with correlated normalization per decay channel',
                            "ZEUS unpublished shape covariance is not reconstructed from fitted slope errors",
                            "H1 2006 common normalization from the original H1 combination, remaining systematics are " + ("fully correlated" if correlated else "diagonal")
                            + ", compare both assumptions because the full covariance is unavailable",
                            "Parameter covariance is conditional on these experimental correlation assumptions",
                            "Integrated cross sections are correlated validation projections, their chi2 is not added to the fit",
                            "No SOFT tune changes or error rescaling"]}


# Report which experiments and missing correlation information control the extracted coupling
def jpsi_diagnostics(proton: model.ProtonVertex, dlog=None):
    nominal = jpsi(proton, dlog=dlog)
    studies = {"nominal": nominal, "H1_2006_correlated": jpsi(proton, correlated=True, dlog=dlog)}
    for experiment in sorted({row[1]["experiment"] for row in jpsi_spectra()} | {"H1", "ZEUS"}):
        studies["without_" + experiment.replace(" ", "_")] = jpsi(proton, exclude=experiment, dlog=dlog)
    return {"parameter_order": ["dsigma_dt_0_nb_per_GeV2", "B_gammaPV_GeV_minus2", "delta", "alpha_prime_GeV_minus2"],
            "amplitude": "A = g_gammaPV * g_Ppp * Fp(t) * exp(B_gammaPV*t/2) * eta(t) * W^2 * (W/W0)^(2*(alpha(t)-1))",
            "cross_section": "d sigma/d|t| = N * Fp(t)^2 * exp(B_gammaPV*t) * (W/W0)^(delta+4*alpha_prime*t)",
            "slope": "B_eff(t,W) = 2*d ln(Fp)/dt + B_gammaPV + 4*alpha_prime*ln(W/W0)",
            "normalization": "g_gammaPV = sqrt(16*pi*N/GeV2nb)/g_Ppp, Fp(0)=1, W0=100 GeV",
            "meson_residue": 'F_gammaPV(t)=exp(B_gammaPV*t/2), FF_prod depends on central mass only',
            "beam_embedding": 'Common gamma-p amplitude for ep, both pp photon directions and UPC',
            "dlog": dlog, "studies": studies}


# Fit light-vector differential spectra with the same fixed proton vertex as the generator
@lru_cache(maxsize=4)
def light(proton: model.ProtonVertex, pdg):
    folder = model.REPOSITORY_ROOT / "HEPData/PHOTOPROD"
    nodes, weights = np.polynomial.legendre.leggauss(32)
    channel = next(row for row in inputs.light(proton) if row.pdg == pdg)
    if pdg == 113:
        folder /= "HEPData-ins1798511-v1-json"
        data = hera_reader.read_rho_wt(inputs.table(1798511, "Table12-13,*"),
                                       folder / "Table12-13statisticalcorrelations.json")
        w0 = channel.w0_gev
        transfer = data["binedges"].mean(axis=1)[:, None] + np.diff(data["binedges"], axis=1) * nodes / 2
        energy = np.log(data["W"] / w0)[:, None]
        chol = np.linalg.cholesky(data["stat_cov"])
        initial = [data["y"].max(), channel.b0_per_gev2.central - 2 * proton.amplitude_slope_per_gev2,
                   4 * (channel.alpha0.central - 1), channel.alpha_prime_per_gev2.central]
        assumptions = ["H1 original W,t data, transfer bin averages and geometric W bin centres as in Section 4.5",
                       "Statistical covariance in the central fit and independent signed systematic offset refits",
                       "Linear effective trajectory and exponential vector residue with the SOFT proton form factor fixed",
                       "Integrated projections and published fit parameters are not reused as independent measurements"]
    elif pdg == 333:
        folder /= "HEPData-ins415642-v1-json"
        data = dict(hera_reader.read(folder / "Table3.json"))
        data["source"] = str(folder / "Table3.json")
        # ZEUS quotes cross sections at the effective transfer points tabulated with the bin limits
        transfer = data["x"][:, None]
        weights = np.ones(1) * 2
        energy = np.zeros_like(transfer)
        chol = np.diag(data["y_err"])
        initial = [channel.forward_dsigma_dt_ub_per_gev2.central,
                   channel.b0_per_gev2.central - 2 * proton.amplitude_slope_per_gev2]
        assumptions = ["ZEUS digitized HEPData differential points at the published effective transfer coordinates",
                       "Tabulated errors are treated as diagonal, no unreported covariance is invented",
                       'Single-W fit: normalization and slope fitted, trajectory assumed']
    else:
        raise ValueError("No light-vector differential fit for this PDG code")
    profile = model.proton_profile(transfer, proton)**2

    # Evaluate the squared factorized amplitude over each experimental transfer bin
    def prediction(parameters):
        forward, slope = parameters[:2]
        exponent = -slope * transfer
        if pdg == 113:
            exponent = exponent + (parameters[2] - 4 * parameters[3] * transfer) * energy
        return forward * np.sum(profile * np.exp(exponent) * weights, axis=1) / 2

    # Fit one measured spectrum or one signed systematic variation
    def solve(values, start):
        result = least_squares(lambda parameters: solve_triangular(chol, prediction(parameters) - values, lower=True),
                               start, bounds=(0.0, np.inf), x_scale="jac")
        if not result.success or np.linalg.matrix_rank(result.jac) != len(start):
            raise ValueError("Light-vector differential fit failed or is not identifiable")
        return result

    fit = solve(data["y"], initial)
    covariance = np.linalg.inv(fit.jac.T @ fit.jac)
    shifts = np.asarray([[solve(data["y"] + offset[:, side], fit.x).x - fit.x for side in range(2)]
                         for offset in data.get("offsets", np.empty((len(data["y"]), 0, 2))).transpose(1, 0, 2)])
    syst_up = np.sqrt(np.sum(np.maximum(shifts.max(axis=1), 0)**2, axis=0)) if len(shifts) else np.zeros(len(fit.x))
    syst_down = np.sqrt(np.sum(np.minimum(shifts.min(axis=1), 0)**2, axis=0)) if len(shifts) else np.zeros(len(fit.x))
    return {"parameters": fit.x.tolist(), "covariance": covariance.tolist(), "syst_up": syst_up.tolist(),
            "syst_down": syst_down.tolist(), "chi2": float(2 * fit.cost), "ndf": len(data["y"]) - len(fit.x),
            "source": data["source"], "data": data["y"].tolist(), "prediction": prediction(fit.x).tolist(),
            "errors": np.sqrt(np.diag(chol @ chol.T)).tolist(),
            "coordinates": {"binedges": data["binedges"].tolist(),
                            "t": data.get("x", data["binedges"].mean(axis=1)).tolist(),
                            "W": data.get("W", np.zeros(len(data["y"]))).tolist()},
            "assumptions": assumptions}


# Convert fitted values and their independent offset errors into one channel constraint
def measurement(fit, index, scale=1.0, offset=0.0):
    error = np.sqrt(fit["covariance"][index][index]) * scale
    return model.Measurement(offset + fit["parameters"][index] * scale, error, error,
                       fit.get("syst_up", [0.0] * len(fit["parameters"]))[index] * scale,
                       fit.get("syst_down", [0.0] * len(fit["parameters"]))[index] * scale, "fit")


# Compute the proton dissociation cross section integrated over the reference mass interval
@njit(cache=False)
def dissociation_spectrum(t, forward, b, n):
    return forward * np.exp(-n * np.log1p(b * t / n))


# Require HEPData JSON covariance before enabling this fit
def jpsi_dissociation():
    raise ValueError('H1 covariance fit unavailable: a supported HEPData JSON covariance input is required')


# Fit the published phi spectrum with the same transfer profile as PhotoDissKernel
# [REFERENCE: ZEUS Collaboration, arXiv:hep-ex/9910038, Table 6 and Section 10.2.2]
def phi_dissociation() -> dict[str, object]:
    path = model.REPOSITORY_ROOT / "HEPData" / "PHOTOPROD" / "HEPData-ins508770-v1-json" / "Table6.json"
    table = hepdata.read_table(path)
    data = []
    normalization = []
    for row in table["values"]:
        value = row["y"][0]
        stat, syst, mass, norm = value["errors"]
        data.append([hepdata.number(row["x"][0]["value"]), hepdata.number(value["value"]),
                     hepdata.error_pair(stat)[0], *hepdata.error_pair(syst), *hepdata.error_pair(mass)])
        normalization.append(float(norm["symerror"].rstrip("%")) / 100.0)
    if not np.allclose(normalization, normalization[0]):
        raise ValueError("Phi table requires one common relative normalization uncertainty")
    data = np.asarray(data)
    error = np.sqrt(data[:, 2] ** 2 + np.max(data[:, 3:5], axis=1) ** 2 + np.max(data[:, 5:7], axis=1) ** 2)
    fit, covariance = curve_fit(
        dissociation_spectrum,
        data[:, 0],
        data[:, 1],
        p0=[3.0, 3.0, 10.0],
        sigma=error,
        absolute_sigma=True,
        bounds=([0.0, 0.0, 1.01], [np.inf, np.inf, np.inf]),
    )
    covariance[0, 0] += (normalization[0] * fit[0]) ** 2
    return {
        "source_url": "https://arxiv.org/abs/hep-ex/9910038",
        "source": "ZEUS Table 6, W=94 GeV, M_Y^2<0.1*W^2",
        "columns": ["t_GeV2", "dsdt_ub_per_GeV2", "stat", "syst_up", "syst_down", "mass_up", "mass_down"],
        "data": data.tolist(),
        "normalization_uncertainty": normalization[0],
        "fit": {
            "parameters": ["forward_ub_per_GeV2", "b_GeV_minus2", "n"],
            "values": fit.tolist(),
            "covariance": covariance.tolist(),
            "chi2": float(np.sum(((dissociation_spectrum(data[:, 0], *fit) - data[:, 1]) / error) ** 2)),
            "ndf": len(data) - len(fit),
        },
        "assumptions": [
            "Diagonal covariance with statistical and larger asymmetric systematic and mass errors in quadrature",
            "Common normalization uncertainty propagated to the fitted forward value",
            "Transfer profile extrapolated outside the measured interval",
        ],
    }


