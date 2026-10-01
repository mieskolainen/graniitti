# HERA published constraints and explicit extrapolations
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pyjson5 as json5
from core.io.serialize import load_json_file
from scipy.optimize import least_squares

from icepack._common import hepdata
from icepack.PHOTOPROD._common import h1_psi2s, zeus_upsilon

from . import model


# Select one original JSON table, excluding saved copies and ambiguous matches
def table(record: int, pattern: str):
    folder = model.REPOSITORY_ROOT / "HEPData/PHOTOPROD" / f"HEPData-ins{record}-v1-json"
    paths = [path for path in folder.glob(pattern) if path.suffix == ".json"]
    if len(paths) != 1:
        raise ValueError(f"Expected one HEPData JSON table for {folder / pattern}, found {len(paths)}")
    return paths[0]


# Read one published parameter through the icepack HEPData decoder
# Signed offset variations are combined on their corresponding upper and lower envelopes
def measurement(record: int, pattern: str, group: int, row: int = 0, scale: float = 1.0):
    path = table(record, pattern)
    value = hepdata.table_values(hepdata.read_table(path), group)[row]
    stat = np.zeros(2)
    syst = np.zeros(2)
    for error in value["errors"]:
        target = stat if error.get("label", "").startswith("stat") else syst
        if "asymerror" in error:
            shifts = [hepdata.number(error["asymerror"][side]) for side in ("plus", "minus")]
            pair = (max(0.0, *shifts), -min(0.0, *shifts))
        else:
            pair = hepdata.error_pair(error)
        target += np.square(pair)
    return model.rescale(model.Measurement(hepdata.number(value["value"]), *np.sqrt(stat), *np.sqrt(syst)), scale)


# Compute the common heavy-vector production-channel definitions
# The TUNE0 vector cards use the SCHC approximation for direct GP helicities
# [REFERENCE: ZEUS, Eur Phys J C 12 (2000) 393, arXiv:hep-ex/9908026]
# [REFERENCE: ZEUS, Phys Lett B 377 (1996) 259, arXiv:hep-ex/9601009]
# [REFERENCE: ZEUS, Eur Phys J C 24 (2002) 345, arXiv:hep-ex/0201043]
# [REFERENCE: H1, JHEP 05 (2010) 032, arXiv:0910.5831]
def heavy_production() -> tuple[tuple[str, int, int], ...]:
    return (
        ("MP", 22, 991),
        ("XP", 22, 991),
        ("GP", 22, 990),
    )


# Build the light-vector HERA channels
def light(proton=None) -> tuple[model.HERAChannel, ...]:
    general_path = model.REPOSITORY_ROOT / "modeldata/TUNE0/GENERAL.json" if proton is None else Path(proton.source.split(":", 1)[0])
    general = load_json_file(general_path, loader=json5.load)
    phi_trajectory = next(row for row in general["PARAM_REGGE"]["photoprod"] if row[0] == 333)
    # [REFERENCE: https://doi.org/10.1140/epjc/s10052-020-08587-3]
    rho = model.HERAChannel(
        key="rho",
        name="rho(770)",
        pdg=113,
        source="H1, Eur. Phys. J. C 80 (2020) 1189, Table 14",
        source_url="https://doi.org/10.1140/epjc/s10052-020-08587-3",
        normalization_origin="H1 forward dsigma/dt",
        trajectory_origin="H1 low-|t| trajectory",
        assumptions=(),
        w0_gev=40.0,
        tmax_gev2=1.5,
        forward_dsigma_dt_ub_per_gev2=measurement(1798511, "Table14,*", 0),
        b0_per_gev2=measurement(1798511, "Table14,*", 2),
        profile_power=measurement(1798511, "Table14,*", 1).central,
        vector_cross_section_slope_per_gev2=4.0,
        normalization_scale=1.0,
        alpha0=measurement(1798511, "Table14,*", 3),
        alpha_prime_per_gev2=measurement(1798511, "Table14,*", 4),
        production_channels=(
            ("MP", 22, 991),
            ("XP", 22, 991),
            ("GP", 22, 990),
        ),
    )
    # [REFERENCE: https://arxiv.org/abs/hep-ex/9601009]
    phi = model.HERAChannel(
        key="phi",
        name="phi(1020)",
        pdg=333,
        source="ZEUS, Phys. Lett. B 377 (1996) 259-272",
        source_url="https://arxiv.org/abs/hep-ex/9601009",
        normalization_origin="ZEUS forward dsigma/dt",
        trajectory_origin="Input tune trajectory held fixed",
        assumptions=("Input tune trajectory, not independently constrained by one W spectrum",),
        w0_gev=70.0,
        tmax_gev2=0.5,
        forward_dsigma_dt_ub_per_gev2=measurement(415642, "Table2.json", 1),
        b0_per_gev2=measurement(415642, "Table2.json", 0),
        profile_power=0.0,
        vector_cross_section_slope_per_gev2=2.0,
        normalization_scale=1.0,
        alpha0=model.Measurement(phi_trajectory[3], 0.0, 0.0, 0.0, 0.0, "model"),
        alpha_prime_per_gev2=model.Measurement(phi_trajectory[4], 0.0, 0.0, 0.0, 0.0, "model"),
        production_channels=(
            ("MP", 22, 991),
            ("XP", 22, 991),
            ("GP", 22, 990),
        ),
    )
    return rho, phi


# Combine the independent published systematic and branching-ratio errors once
def psi2s_ratio():
    data = h1_psi2s.measurement()
    syst = math.hypot(data["syst"], data["branching"])
    return model.Measurement(data["ratio"], data["stat"], data["stat"], syst, syst)


# Build the ground-state and radially excited charmonium channels
def charmonium() -> tuple[model.HERAChannel, ...]:
    # [REFERENCE: https://arxiv.org/abs/hep-ex/0201043]
    jpsi = model.HERAChannel(
        key="jpsi",
        name="J/psi(1S)",
        pdg=443,
        source=("ZEUS, Eur. Phys. J. C 24 (2002) 345-360, Table 1, muon channel, 90 < W < 110 GeV"),
        source_url="https://arxiv.org/abs/hep-ex/0201043",
        normalization_origin="ZEUS forward dsigma/dt",
        trajectory_origin="ZEUS effective trajectory",
        assumptions=(),
        w0_gev=100.0,
        tmax_gev2=1.8,
        forward_dsigma_dt_ub_per_gev2=measurement(582237, "Table1.json", 1, row=4, scale=1.0e-3),
        b0_per_gev2=measurement(582237, "Table1.json", 2, row=4),
        profile_power=0.0,
        vector_cross_section_slope_per_gev2=0.0,
        normalization_scale=1.0,
        alpha0=model.Measurement(1.0, 0.0, 0.0, 0.0, 0.0, "seed"),
        alpha_prime_per_gev2=model.Measurement(0.0, 0.0, 0.0, 0.0, 0.0, "seed"),
        production_channels=heavy_production(),
    )

    # [REFERENCE: arXiv:hep-ex/0205107, Section 3.1]
    psi2s_slope = jpsi.b0_per_gev2
    ratio = psi2s_ratio()
    psi2s = model.HERAChannel(
        key="psi2s",
        name="psi(2S)",
        pdg=100443,
        source="H1 1996-2000 diffractive psi(2S)/J/psi ratio",
        source_url="https://arxiv.org/abs/hep-ex/0205107",
        normalization_origin="H1 combined diffractive ratio times the fitted J/psi amplitude",
        trajectory_origin="J/psi profile and trajectory assumed",
        assumptions=(
            "Common psi(2S) and J/psi transfer profile, energy dependence and proton transition ratio",
            "Diffractive ratio is applied to both elastic and dissociative amplitudes",
            "Equal elastic transfer profiles are compatible with the measured H1 slopes",
            "Equal energy powers are compatible with the measured H1 difference within its uncertainty",
            "independent input uncertainties",
        ),
        w0_gev=100.0,
        tmax_gev2=1.0,
        forward_dsigma_dt_ub_per_gev2=model.from_ratio(
            jpsi.forward_dsigma_dt_ub_per_gev2,
            jpsi.b0_per_gev2,
            psi2s_slope,
            ratio,
            1.0,
        ),
        b0_per_gev2=psi2s_slope,
        profile_power=0.0,
        vector_cross_section_slope_per_gev2=0.0,
        normalization_scale=1.0,
        alpha0=jpsi.alpha0,
        alpha_prime_per_gev2=jpsi.alpha_prime_per_gev2,
        production_channels=heavy_production(),
    )
    return jpsi, psi2s


# Fit independent H1 and full-statistics ZEUS bottomonium cross sections and the measured t slope
# [REFERENCE: arXiv:hep-ex/0003020, arXiv:0903.4205, arXiv:1111.2133]
def bottomonium(proton=None) -> tuple[model.HERAChannel, ...]:
    general_path = model.REPOSITORY_ROOT / "modeldata/TUNE0/GENERAL.json" if proton is None else Path(proton.source.split(":", 1)[0])
    proton = model.load_proton(general_path) if proton is None else proton
    general = load_json_file(general_path, loader=json5.load)
    _, w0, _, _, alpha_prime = next(row for row in general["PARAM_REGGE"]["photoprod"] if row[0] == 553)
    dlog = model.read_dlog(general, 553)
    zeus = zeus_upsilon.energy()
    b = model.Measurement(*zeus_upsilon.slope())
    slope, _ = model.match_slope(proton, b.central, 0.0, zeus_upsilon.slope_tmax())
    slope -= 4 * alpha_prime * math.log(zeus["W"][-1] / w0)
    integral, _ = model.profile_moments(proton, slope, None)
    h1path = table(525022, "Table5.json")
    energy = np.r_[zeus["W"][:2], hepdata.number(hepdata.read_table(h1path)["values"][0]["x"][0]["value"])]
    observations = [model.Measurement(value, stat, stat, up, down) for value, stat, up, down in
                    zip(zeus["sigma_pb"][:2] * 1e-6, zeus["stat_pb"][:2] * 1e-6,
                        zeus["syst_up_pb"][:2] * 1e-6, zeus["syst_down_pb"][:2] * 1e-6, strict=True)]
    observations.append(measurement(525022, "Table5.json", 0, scale=1e-3))
    observed = np.asarray([value.central for value in observations])
    up = np.asarray([math.hypot(value.uncertainty_up, value.syst_up) for value in observations])
    down = np.asarray([math.hypot(value.uncertainty_down, value.syst_down) for value in observations])
    transfer = np.asarray([model.profile_moments(proton, slope + 4 * alpha_prime * math.log(w / w0), None)[0]
                           / integral for w in energy])

    if dlog is not None:
        transfer *= model.dlog_factor(energy, w0, *dlog)

    # Keep the forward energy power distinct from the t-integrated cross-section power
    def prediction(parameters):
        return parameters[0] * (energy / w0)**parameters[1] * transfer

    # Keep asymmetric experimental errors on their corresponding residual sides
    def residual(parameters):
        difference = prediction(parameters) - observed
        return difference / np.where(difference >= 0, up, down)

    delta = np.log(observed[1] * transfer[0] / (observed[0] * transfer[1])) / np.log(energy[1] / energy[0])
    fit = least_squares(residual, [observed[0], max(delta, 0.0)], bounds=(0, np.inf), x_scale="jac",
                        ftol=1e-11, xtol=1e-11, gtol=1e-11)
    if not fit.success:
        raise ValueError("Upsilon W fit failed")
    covariance = np.linalg.inv(fit.jac.T @ fit.jac)
    sigma = model.Measurement(fit.x[0], math.sqrt(covariance[0, 0]), math.sqrt(covariance[0, 0]), 0.0, 0.0, "fit")
    alpha = model.Measurement(1 + fit.x[1] / 4, math.sqrt(covariance[1, 1]) / 4,
                        math.sqrt(covariance[1, 1]) / 4, 0.0, 0.0, "fit")
    constraints = {"parameters": fit.x.tolist(), "parameter_order": ["sigma_W0_ub", "delta_forward"],
                   "covariance": covariance.tolist(), "chi2": float(2 * fit.cost), "ndf": len(observed) - len(fit.x),
                   "source": [zeus["source"], str(h1path)], "W": energy.tolist(),
                   "data": observed.tolist(), "prediction": prediction(fit.x).tolist(),
                   "errors_up": up.tolist(), "errors_down": down.tolist(),
                   "assumptions": ["Independent H1 and full-statistics ZEUS resolved 1S cross sections at published effective W",
                                   "Earlier overlapping ZEUS data and the full-range projection are not independent fit constraints",
                                   'Statistical and asymmetric systematic errors in quadrature, no published covariance',
                                   'ZEUS slope matched via mean |t| over the measured interval at effective W',
                                   "Shrinkage retained from the input tune because one slope measurement cannot determine it",
                                   "Excited-state coupling ratios retained from the input tune, not independently measured at HERA"]}
    channels = []
    couplings = []
    for filename in ("Y1S.json", "Y2S.json", "Y3S.json"):
        card = load_json_file(general_path.parent / "RES" / filename, loader=json5.load)["PARAM_RES"]
        helicity = card["MODELS"]["GP"]["[22,990]"]["helicity"]
        if len(helicity) != 1 or len(helicity[0]) != 4 or helicity[0][:2] != [-1, 0]:
            raise ValueError("Upsilon extrapolation requires one transverse SCHC reference coupling")
        coupling = float(helicity[0][2])
        if not math.isfinite(coupling) or coupling < 0:
            raise ValueError("Upsilon reference couplings must be finite and nonnegative")
        couplings.append(coupling)
    if couplings[0] <= 0:
        raise ValueError("Upsilon(1S) needs a nonzero reference coupling for excited-state ratios")
    for index, pdg in enumerate((553, 100553, 200553)):
        ratio = (couplings[index] / couplings[0])**2
        channels.append(model.HERAChannel(
            key=f"upsilon{index + 1}s", name=f"Upsilon({index + 1}S)", pdg=pdg,
            source='H1 and ZEUS Upsilon(1S) data with excited-state extrapolations',
            source_url="https://arxiv.org/abs/hep-ex/0003020",
            normalization_origin="Joint published Upsilon(1S) cross-section fit",
            trajectory_origin="Independent HERA energy fit with measured transfer slope", assumptions=tuple(constraints["assumptions"]),
            w0_gev=w0, tmax_gev2=None, forward_dsigma_dt_ub_per_gev2=model.rescale(sigma, ratio / integral),
            b0_per_gev2=b, profile_power=0.0, vector_cross_section_slope_per_gev2=slope,
            normalization_scale=1.0, alpha0=alpha,
            alpha_prime_per_gev2=model.Measurement(alpha_prime, 0.0, 0.0, 0.0, 0.0, "model"),
            production_channels=heavy_production(), dlog=dlog, fit=constraints if index == 0 else {"reference_pdg": 553}))
    return tuple(channels)
