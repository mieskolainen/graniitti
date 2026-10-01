# HERA factorized amplitude, proton vertex and measured-constraint propagation
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pyjson5 as json5
from core.io.serialize import load_json_file
from numba import njit
from scipy.integrate import quad
from scipy.optimize import brentq

from .. import ROOT, soft_exchange

REPOSITORY_ROOT = ROOT


@dataclass(frozen=True)
class Measurement:
    central: float
    uncertainty_up: float
    uncertainty_down: float
    syst_up: float
    syst_down: float
    uncertainty_label: str = "stat"


@dataclass(frozen=True)
class HERAChannel:
    key: str
    name: str
    pdg: int
    source: str
    source_url: str
    normalization_origin: str
    trajectory_origin: str
    assumptions: tuple[str, ...]
    w0_gev: float
    tmax_gev2: float | None
    forward_dsigma_dt_ub_per_gev2: Measurement
    b0_per_gev2: Measurement
    profile_power: float
    vector_cross_section_slope_per_gev2: float
    normalization_scale: float
    alpha0: Measurement
    alpha_prime_per_gev2: Measurement
    production_channels: tuple[tuple[str, int, int], ...]
    fit: dict | None = None
    diss_reference: HERAChannel | None = None
    dlog: tuple[float, float] | None = None


@dataclass(frozen=True)
class DerivedChannel:
    gamma_pomeron_vector_coupling_per_gev: float
    coupling_stat_up_per_gev: float
    coupling_stat_down_per_gev: float
    coupling_syst_up_per_gev: float
    coupling_syst_down_per_gev: float
    reference_amplitude_abs: float
    reconstructed_forward_dsigma_dt_ub_per_gev2: float
    exponential_sigma_ub: float
    matched_sigma_ub: float
    forward_w_exponent: float
    shrinkage_per_gev2: float
    vector_cross_section_slope_per_gev2: float
    factorized_cross_slope_per_gev2: float


@dataclass(frozen=True)
class ProtonVertex:
    beam_residue_per_gev: float
    form_factor: str
    parameters: tuple[tuple[float, ...], ...]
    amplitude_slope_per_gev2: float
    source: str
    weights: tuple[tuple[float, ...], ...]
    transition: str


# Read an optional fixed-scale double-log term from the same physical input as the generator
def read_dlog(general, pdg):
    rows = [row for row in general["PARAM_REGGE"].get("photoprod_dlog", []) if row[0] == pdg]
    if not rows:
        return None
    w0 = next(row[1] for row in general["PARAM_REGGE"]["photoprod"] if row[0] == pdg)
    if len(rows) != 1 or len(rows[0]) != 3:
        raise ValueError("photoprod_dlog requires one [PDG, scale2, c] row per channel")
    scale2, coefficient = rows[0][1:]
    if not (0 < scale2 < w0**2 and math.isfinite(coefficient) and coefficient >= 0):
        raise ValueError("Invalid fixed-scale double-log photoproduction parameters")
    return (scale2, coefficient)


# Compute the fixed-scale double-log cross-section factor with unit value and zero local power at W0
# [REFERENCE: arXiv:1307.7099, Eq. (17)]
@njit(cache=False)
def dlog_factor(w, w0, scale2, coefficient):
    root0 = np.sqrt(np.log(w0 * w0 / scale2))
    return np.exp(coefficient * (np.sqrt(np.log(w * w / scale2)) - root0 - np.log(w / w0) / root0))


# Compute the local contribution of the double-log term to d ln(dsigma/dt) / d ln W
@njit(cache=False)
def dlog_power(w, w0, scale2, coefficient):
    return coefficient * (1 / np.sqrt(np.log(w * w / scale2)) - 1 / np.sqrt(np.log(w0 * w0 / scale2)))


# Convert one natural-unit cross section from GeV^-2 to microbarn
@njit(cache=False)
def gev2_to_microbarn() -> float:
    hbar_gev_s = 6.582119514e-25
    speed_of_light_m_per_s = 2.99792458e8
    barn_m2 = 1.0e-28
    return (hbar_gev_s * speed_of_light_m_per_s) ** 2 / barn_m2 * 1.0e6


# Derive the transition coupling against one fixed bare proton coupling
@njit(cache=False)
def coupling(
    forward_dsigma_dt_ub_per_gev2: float,
    proton_pomeron_beam_residue_per_gev: float,
) -> float:
    forward_gev4 = forward_dsigma_dt_ub_per_gev2 / gev2_to_microbarn()
    return math.sqrt(16.0 * math.pi * forward_gev4) / proton_pomeron_beam_residue_per_gev


# Reconstruct the forward differential cross section from an amplitude coupling
@njit(cache=False)
def forward(
    gamma_pomeron_vector_coupling_per_gev: float,
    proton_pomeron_beam_residue_per_gev: float,
) -> float:
    residue2 = (gamma_pomeron_vector_coupling_per_gev * proton_pomeron_beam_residue_per_gev) ** 2
    return residue2 / (16.0 * math.pi) * gev2_to_microbarn()


# Compute exponential integral factor including the zero slope limit
@njit(cache=False)
def exp_integral(
    b0_per_gev2: float,
    tmax_gev2: float,
) -> float:
    exponent = -b0_per_gev2 * tmax_gev2
    ratio = (
        1.0 + exponent / 2.0 + exponent * exponent / 6.0
        if abs(exponent) < 1.0e-8
        else math.expm1(exponent) / exponent
    )
    return tmax_gev2 * ratio


# Integrate dsigma/d|t| = dsigma/dt(0) exp(-B|t|) over the measured interval
@njit(cache=False)
def sigma_exp(
    forward_dsigma_dt_ub_per_gev2: float,
    b0_per_gev2: float,
    tmax_gev2: float,
) -> float:
    return forward_dsigma_dt_ub_per_gev2 * exp_integral(b0_per_gev2, tmax_gev2)


# Integrate the matched HERA differential profile over a finite |t| interval
@njit(cache=False)
def sigma_profile(
    forward_dsigma_dt_ub_per_gev2: float,
    b_cross_section_per_gev2: float,
    profile_power: float,
    tmax_gev2: float,
) -> float:
    if profile_power <= 0.0:
        return sigma_exp(
            forward_dsigma_dt_ub_per_gev2,
            b_cross_section_per_gev2,
            tmax_gev2,
        )
    profile = b_cross_section_per_gev2 * tmax_gev2 / profile_power
    log_profile = math.log1p(profile)
    log_ratio = (
        1.0 - profile / 2.0 + profile * profile / 3.0
        if abs(profile) < 1.0e-8
        else log_profile / profile
    )
    exponent = (1.0 - profile_power) * log_profile
    # Compute the logarithmic limit without cancellation near unit profile power
    ratio = 1.0 if abs(exponent) < 1.0e-15 else math.expm1(exponent) / exponent
    integral = tmax_gev2 * log_ratio * ratio
    return forward_dsigma_dt_ub_per_gev2 * integral


# Convert a cross section integrated over all |t| to its exponential forward value
@njit(cache=False)
def forward_from_total(
    total_cross_section_ub: float,
    b0_per_gev2: float,
) -> float:
    return total_cross_section_ub * b0_per_gev2


# Derive a daughter forward value from a parent integrated-cross-section ratio
@njit(cache=False)
def forward_ratio(
    parent_forward_dsigma_dt_ub_per_gev2: float,
    parent_b0_per_gev2: float,
    daughter_b0_per_gev2: float,
    integrated_cross_section_ratio: float,
    tmax_gev2: float,
) -> float:
    parent_sigma = parent_forward_dsigma_dt_ub_per_gev2 * exp_integral(
        parent_b0_per_gev2, tmax_gev2
    )
    return (
        integrated_cross_section_ratio
        * parent_sigma
        / exp_integral(daughter_b0_per_gev2, tmax_gev2)
    )


# Compute exact asymmetric coupling shifts induced by one data uncertainty
@njit(cache=False)
def coupling_error(
    central_dsigma: float,
    uncertainty_up: float,
    uncertainty_down: float,
    proton_pomeron_beam_residue_per_gev: float,
) -> tuple[float, float]:
    central = coupling(central_dsigma, proton_pomeron_beam_residue_per_gev)
    upper = coupling(
        central_dsigma + uncertainty_up, proton_pomeron_beam_residue_per_gev
    )
    lower = coupling(
        max(0.0, central_dsigma - uncertainty_down),
        proton_pomeron_beam_residue_per_gev,
    )
    return upper - central, central - lower


# Combine independent one-at-a-time variations in quadrature
def quadrature_shift(central: float, variations: tuple[float, ...]) -> float:
    return math.sqrt(sum((value - central) ** 2 for value in variations))


# Read the active bare proton coupling and form factor from one GRANIITTI tune
def load_proton(general_path: Path) -> ProtonVertex:
    with general_path.open(encoding="utf-8") as stream:
        general = load_json_file(stream.name, loader=json5.load)
    if general["PARAM_REGGE"]["photoprod_eta_mode"] not in {"rotating_t0", "rotating"}:
        raise ValueError("PARAM_REGGE.photoprod_eta_mode must be rotating_t0 or rotating")
    configuration = soft_exchange.load(general)
    pomeron = soft_exchange.exchange(configuration, "pomeron")
    residue = soft_exchange.proton_vertex(
        configuration.model,
        pomeron.soft_exchange,
        configuration.model_name,
    )
    exchange = configuration.model["EXCHANGE"][pomeron.soft_exchange]
    matrix = np.asarray(exchange["g"], dtype=float)
    proton = np.asarray(soft_exchange.proton_state(configuration.model["GW"]["theta"], len(matrix)))
    weights = matrix * np.outer(proton, proton)
    if exchange["transition_ff"] == "diagonal":
        weights *= np.eye(len(matrix))
    weights /= weights.sum()
    return ProtonVertex(
        beam_residue_per_gev=residue.beam_residue_per_gev,
        form_factor=residue.form_factor,
        parameters=residue.parameters,
        amplitude_slope_per_gev2=residue.amplitude_slope_per_gev2,
        source=f"{general_path}:{residue.source}",
        weights=tuple(tuple(row) for row in weights),
        transition=exchange["transition_ff"],
    )


# Evaluate the physical SOFT proton form factor at x = |t|
def proton_profile(x, proton: ProtonVertex):
    x = np.asarray(x, dtype=float)
    profiles = []
    for row in proton.parameters:
        if proton.form_factor == "EXP":
            value = np.exp(-0.5 * row[0] * x)
        elif proton.form_factor == "DPOW":
            value = np.exp(-row[2] * (np.log1p(x / row[0]) + np.log1p(x / row[1])))
        elif proton.form_factor == "EXPOW":
            value = np.exp(-(row[0] * (row[1] + x)) ** row[2] + (row[0] * row[1]) ** row[2])
        elif proton.form_factor == "MIXEXP":
            value = (1.0 - row[0]) * np.exp(-0.5 * row[1] * x) + row[0] * np.exp(-0.5 * row[2] * x)
            if row[3] > 0.0:
                value *= 1.0 - x / row[3]
        elif proton.form_factor == "ODD3G_NODE":
            value = np.exp(-0.5 * row[0] * x) * (1.0 - x / row[1]) / (1.0 + x / row[2]) ** 3
        elif proton.form_factor == "3G":
            value = (1.0 + x / row[0]) ** -2
        elif proton.form_factor == "GKERNEL":
            offset = len(row) % 4
            value = np.ones_like(x) if not offset else 1.0 + row[0] * x
            for a, p, nu, mu2 in np.asarray(row[offset:]).reshape(-1, 4):
                y = a * ((x + mu2) ** p - mu2 ** p)
                value *= np.exp(-y if nu <= 1.0e-12 else -np.log1p(nu * y) / nu)
        else:
            raise ValueError(f"unsupported SOFT proton form factor: {proton.form_factor}")
        profiles.append(value)
    result = np.zeros_like(x)
    for i, first in enumerate(profiles):
        for j, second in enumerate(profiles):
            weight = proton.weights[i][j]
            if i == j:
                transition = first
            elif proton.transition == "diagonal":
                continue
            elif proton.transition == "arithmetic":
                transition = (first + second) / 2.0
            else:
                if np.any(first * second < -1.0e-14):
                    raise ValueError("SOFT geometric transition has opposite form factor signs")
                transition = np.sign(first) * np.sqrt(np.maximum(first * second, 0.0))
            result += weight * transition
    return result


# Integrate the actual GP cross section profile and its first transfer moment
def profile_moments(proton: ProtonVertex, slope: float, tmax: float | None):
    upper = math.inf if tmax is None else tmax
    moments = [quad(lambda x, power=power: x**power * float(proton_profile(x, proton)) ** 2 * math.exp(-slope * x),
                    0.0, upper, epsabs=1.0e-10, epsrel=1.0e-9)[0] for power in (0, 1)]
    if moments[0] <= 0.0 or not np.isfinite(moments).all():
        raise ValueError("SOFT proton profile requires finite positive cross section moments")
    return moments[0], moments[1] / moments[0]


# Match the measured finite-interval mean |t| with the fixed SOFT profile
# This is the likelihood projection onto Fp(-x)^2 exp(-B_vector*x), with its normalization free
# A fitted exponential B is not the local derivative of the proton form factor at t=0
def match_slope(proton: ProtonVertex, b: float, n: float, tmax: float | None):
    upper = math.inf if tmax is None else tmax
    # Evaluate the published exponential or power profile
    def measured(x):
        return math.exp(-b * x) if n <= 0.0 else math.exp(-n * math.log1p(b * x / n))

    integral = quad(measured, 0.0, upper)[0]
    mean = quad(lambda x: x * measured(x), 0.0, upper)[0] / integral
    if profile_moments(proton, 0.0, tmax)[1] < mean:
        raise ValueError("Measured t spectrum is broader than the fixed SOFT residue permits with B_vector >= 0")
    upper_b = b
    while profile_moments(proton, upper_b, tmax)[1] > mean:
        upper_b *= 2.0
    slope = brentq(lambda value: profile_moments(proton, value, tmax)[1] - mean, 0.0, upper_b)
    return slope, integral / profile_moments(proton, slope, tmax)[0]


# Scale one measured normalization by an exact fixed factor
def rescale(value: Measurement, factor: float) -> Measurement:
    return Measurement(
        central=value.central * factor,
        uncertainty_up=value.uncertainty_up * factor,
        uncertainty_down=value.uncertainty_down * factor,
        syst_up=value.syst_up * factor,
        syst_down=value.syst_down * factor,
        uncertainty_label=value.uncertainty_label,
    )


# Shift one measurement along one asymmetric uncertainty family
def shift_value(
    value: Measurement,
    family: str,
    upward: bool,
) -> float:
    if family == "first":
        shift = value.uncertainty_up if upward else -value.uncertainty_down
    elif family == "syst":
        shift = value.syst_up if upward else -value.syst_down
    else:
        raise ValueError(f"unknown uncertainty family: {family}")
    return value.central + shift


# Build one-at-a-time variations for a total-cross-section conversion
def total_variations(
    total_cross_section_ub: Measurement,
    slope_per_gev2: Measurement,
    family: str,
    upward: bool,
) -> tuple[float, ...]:
    evaluate = forward_from_total
    return (
        evaluate(
            shift_value(total_cross_section_ub, family, upward),
            slope_per_gev2.central,
        ),
        evaluate(
            total_cross_section_ub.central,
            shift_value(slope_per_gev2, family, upward),
        ),
    )


# Propagate a total-cross-section and slope measurement to the forward value
def from_total(
    total_cross_section_ub: Measurement,
    slope_per_gev2: Measurement,
) -> Measurement:
    central = forward_from_total(
        total_cross_section_ub.central, slope_per_gev2.central
    )
    errors = [
        quadrature_shift(
            central,
            total_variations(total_cross_section_ub, slope_per_gev2, family, upward),
        )
        for family, upward in (
            ("first", True),
            ("first", False),
            ("syst", True),
            ("syst", False),
        )
    ]
    return Measurement(
        central,
        errors[0],
        errors[1],
        errors[2],
        errors[3],
        total_cross_section_ub.uncertainty_label,
    )


# Build one-at-a-time variations for a parent-to-daughter ratio conversion
def ratio_variations(
    parent_forward: Measurement,
    parent_slope: Measurement,
    daughter_slope: Measurement,
    ratio: Measurement,
    tmax_gev2: float,
    family: str,
    upward: bool,
) -> tuple[float, ...]:
    evaluate = forward_ratio
    central = (
        parent_forward.central,
        parent_slope.central,
        daughter_slope.central,
        ratio.central,
        tmax_gev2,
    )
    return (
        evaluate(shift_value(parent_forward, family, upward), *central[1:]),
        evaluate(
            central[0],
            shift_value(parent_slope, family, not upward),
            *central[2:],
        ),
        evaluate(
            *central[:2],
            shift_value(daughter_slope, family, upward),
            *central[3:],
        ),
        evaluate(
            *central[:3],
            shift_value(ratio, family, upward),
            central[4],
        ),
    )


# Propagate a measured parent-to-daughter integrated cross-section ratio
def from_ratio(
    parent_forward: Measurement,
    parent_slope: Measurement,
    daughter_slope: Measurement,
    ratio: Measurement,
    tmax_gev2: float,
) -> Measurement:
    central = forward_ratio(
        parent_forward.central,
        parent_slope.central,
        daughter_slope.central,
        ratio.central,
        tmax_gev2,
    )
    errors = [
        quadrature_shift(
            central,
            ratio_variations(
                parent_forward,
                parent_slope,
                daughter_slope,
                ratio,
                tmax_gev2,
                family,
                upward,
            ),
        )
        for family, upward in (
            ("first", True),
            ("first", False),
            ("syst", True),
            ("syst", False),
        )
    ]
    return Measurement(
        central,
        errors[0],
        errors[1],
        errors[2],
        errors[3],
        "propagated stat",
    )


# Derive all amplitude and cross-section quantities for one HERA channel
def derive(
    channel: HERAChannel,
    proton: ProtonVertex,
) -> DerivedChannel:
    measured = channel.forward_dsigma_dt_ub_per_gev2
    fitted = rescale(measured, channel.normalization_scale)
    proton_pomeron_beam_residue_per_gev = proton.beam_residue_per_gev
    fitted_forward = fitted.central
    g = coupling(fitted_forward, proton_pomeron_beam_residue_per_gev)
    stat_up, stat_down = coupling_error(
        fitted_forward,
        fitted.uncertainty_up,
        fitted.uncertainty_down,
        proton_pomeron_beam_residue_per_gev,
    )
    syst_up, syst_down = coupling_error(
        fitted_forward,
        fitted.syst_up,
        fitted.syst_down,
        proton_pomeron_beam_residue_per_gev,
    )
    reference_s = channel.w0_gev**2
    return DerivedChannel(
        gamma_pomeron_vector_coupling_per_gev=g,
        coupling_stat_up_per_gev=stat_up,
        coupling_stat_down_per_gev=stat_down,
        coupling_syst_up_per_gev=syst_up,
        coupling_syst_down_per_gev=syst_down,
        reference_amplitude_abs=(g * proton_pomeron_beam_residue_per_gev * reference_s),
        reconstructed_forward_dsigma_dt_ub_per_gev2=forward(
            g, proton_pomeron_beam_residue_per_gev
        ),
        exponential_sigma_ub=(
            measured.central / channel.b0_per_gev2.central
            if channel.tmax_gev2 is None
            else sigma_exp(
                measured.central,
                channel.b0_per_gev2.central,
                channel.tmax_gev2,
            )
        ),
        matched_sigma_ub=fitted_forward * profile_moments(
            proton, channel.vector_cross_section_slope_per_gev2, channel.tmax_gev2
        )[0],
        forward_w_exponent=4.0 * (channel.alpha0.central - 1.0),
        shrinkage_per_gev2=4.0 * channel.alpha_prime_per_gev2.central,
        vector_cross_section_slope_per_gev2=(channel.vector_cross_section_slope_per_gev2),
        factorized_cross_slope_per_gev2=(
            2.0 * proton.amplitude_slope_per_gev2
            + channel.vector_cross_section_slope_per_gev2
        ),
    )


