# Derive photoproduction couplings and steering rows from HERA measurements
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import dataclasses
import math

import numpy as np
import pyjson5 as json5
from core.io import serialize

from . import fit, inputs, model


# Derive channel couplings and finite-interval slopes with the SOFT tune held fixed
def match(channels, proton: model.ProtonVertex, config=None):
    by_pdg = {channel.pdg: channel for channel in channels}
    for pdg in (113, 333):
        result = fit.light(proton, pdg)
        channel = by_pdg[pdg]
        by_pdg[pdg] = dataclasses.replace(
            channel,
            fit=result,
            normalization_origin="Absolute differential HEPData fit with fixed SOFT residue",
            assumptions=tuple(result["assumptions"]),
            forward_dsigma_dt_ub_per_gev2=fit.measurement(result, 0),
            vector_cross_section_slope_per_gev2=result["parameters"][1],
            alpha0=fit.measurement(result, 2, scale=0.25, offset=1.0) if pdg == 113 else channel.alpha0,
            alpha_prime_per_gev2=fit.measurement(result, 3) if pdg == 113 else channel.alpha_prime_per_gev2,
        )
    general = serialize.load_json_file(proton.source.split(":", 1)[0], loader=json5.load)
    dlog = model.read_dlog(general, 443)
    result = fit.jpsi(proton, dlog=dlog)
    forward, slope, delta, alpha_prime = result["parameters"]
    error = np.sqrt(np.diag(result["covariance"]))
    by_pdg[443] = dataclasses.replace(
        by_pdg[443],
        source="H1 2013, H1 2006 and ZEUS elastic J/psi absolute differential cross sections",
        normalization_origin="Absolute HERA W,t spectra with the physical SOFT residue fixed",
        trajectory_origin="Joint W and t dependence, including shrinkage",
        assumptions=tuple(result["assumptions"]),
        fit=result,
        dlog=dlog,
        forward_dsigma_dt_ub_per_gev2=model.Measurement(
            forward / 1000.0, error[0] / 1000.0, error[0] / 1000.0, 0.0, 0.0, "fit"
        ),
        vector_cross_section_slope_per_gev2=slope,
        alpha0=model.Measurement(1.0 + delta / 4.0, error[2] / 4.0, error[2] / 4.0, 0.0, 0.0, "fit"),
        alpha_prime_per_gev2=model.Measurement(alpha_prime, error[3], error[3], 0.0, 0.0, "fit"),
    )
    if config is not None:
        from develop.tools.lib.hera import lhcb

        lhcb.fit(by_pdg, proton, config)
    # A shared heavy-vector profile preserves the measured quasi-elastic ratio for either proton final state
    parent, child = by_pdg[443], by_pdg[100443]
    ratio = inputs.psi2s_ratio()
    norm = parent.forward_dsigma_dt_ub_per_gev2
    forward = model.Measurement(
        norm.central * ratio.central,
        math.hypot(ratio.central * norm.uncertainty_up, norm.central * ratio.uncertainty_up),
        math.hypot(ratio.central * norm.uncertainty_down, norm.central * ratio.uncertainty_down),
        math.hypot(ratio.central * norm.syst_up, norm.central * ratio.syst_up),
        math.hypot(ratio.central * norm.syst_down, norm.central * ratio.syst_down),
        "fit",
    )
    by_pdg[100443] = dataclasses.replace(
        child,
        forward_dsigma_dt_ub_per_gev2=forward,
        vector_cross_section_slope_per_gev2=parent.vector_cross_section_slope_per_gev2,
        alpha0=parent.alpha0,
        alpha_prime_per_gev2=parent.alpha_prime_per_gev2,
        dlog=parent.dlog,
    )
    return tuple(
        dataclasses.replace(
            by_pdg[channel.pdg], diss_reference=parent if channel.pdg in (100443, 553, 100553, 200553) else None
        )
        for channel in channels
    )


# Build all published HERA inputs and adopted trajectory parameters
def channels(proton: model.ProtonVertex | None = None, config=None) -> tuple[model.HERAChannel, ...]:
    if proton is None:
        proton = model.load_proton(model.REPOSITORY_ROOT / "modeldata" / "TUNE0" / "GENERAL.json")
    return match(inputs.light(proton) + inputs.charmonium() + inputs.bottomonium(proton), proton, config)


# Derive dissociation rows against the elastic couplings produced in the same HERA calculation
def dissociation(channel: model.HERAChannel, derived: model.DerivedChannel) -> dict[str, object] | None:
    # Transfer the J/psi proton transition to other heavy vectors as an explicit model assumption
    # [REFERENCE: ZEUS Collaboration, arXiv:0903.4205, Section 6]
    reference = channel.diss_reference or channel
    if reference.pdg == 333:
        result = fit.phi_dissociation()
        forward, b, n = result["fit"]["values"]
        W0, delta, mass_max = 94.0, derived.forward_w_exponent, 94.0 * math.sqrt(0.1)
        result["assumptions"].append("Phi energy power from the elastic trajectory")
    elif reference.pdg in (113, 443):
        # [REFERENCE: H1 Collaboration, arXiv:2005.14471, Tables 8, 11 and Eq. (25)]
        # [REFERENCE: H1 Collaboration, arXiv:1304.5162, Tables 2, 3]
        W0, sigma, tmax, b, n, delta, source = {
            113: (
                40.0,
                inputs.measurement(1798511, "Table8,*", 2).central,
                1.5,
                inputs.measurement(1798511, "Table11,*", 5).central,
                inputs.measurement(1798511, "Table11,*", 4).central,
                inputs.measurement(1798511, "Table8,*", 3).central,
                "2005.14471",
            ),
            443: tuple(fit.jpsi_dissociation()[key] for key in ("W0", "sigma_ub"))
            + (8.0,)
            + tuple(fit.jpsi_dissociation()[key] for key in ("b", "n", "delta"))
            + ("1304.5162",),
        }[reference.pdg]
        forward = sigma / model.sigma_profile(1.0, b, n, tmax)
        mass_max = 10.0
        result = {
            "source_url": f"https://arxiv.org/abs/{source}",
            "sigma_ub": sigma,
            "tmax_GeV2": tmax,
            "assumptions": [],
        }
    else:
        return None
    if reference.pdg == 443:
        result["fit"] = fit.jpsi_dissociation()
        result["assumptions"].extend(result["fit"]["assumptions"])
    result["assumptions"].append("Smooth rho mass continuum with epsilon=0.0808 down to threshold")
    if channel.pdg != 113:
        result["assumptions"].append("Mass power borrowed from rho, not a measured channel mass spectrum")
    elastic = derived.reconstructed_forward_dsigma_dt_ub_per_gev2 * (W0 / channel.w0_gev) ** derived.forward_w_exponent
    if channel.dlog is not None:
        elastic *= model.dlog_factor(W0, channel.w0_gev, *channel.dlog)
    if reference.pdg != channel.pdg:
        power = 4.0 * (reference.alpha0.central - 1.0)
        reference_elastic = reference.normalization_scale * reference.forward_dsigma_dt_ub_per_gev2.central
        reference_elastic *= (W0 / reference.w0_gev) ** power
        if reference.dlog is not None:
            reference_elastic *= model.dlog_factor(W0, reference.w0_gev, *reference.dlog)
        scale = elastic / reference_elastic
        forward *= scale
        result["sigma_ub"] *= scale
        delta += derived.forward_w_exponent - power
        terms = {}
        for state, sign in ((channel, 1), (reference, -1)):
            if state.dlog is not None:
                scale2, coefficient = state.dlog
                delta += sign * model.dlog_power(W0, state.w0_gev, scale2, coefficient)
                terms[scale2] = terms.get(scale2, 0.0) + sign * coefficient
        result["photoprod_diss_dlog"] = [
            [channel.pdg, round(scale2, 4), round(coefficient, 4)]
            for scale2, coefficient in sorted(terms.items())
            if abs(coefficient) > 0
        ]
        result["reference_pdg"] = reference.pdg
        result["normalization_origin"] = "J/psi forward pd/el ratio at a common W and proton mass interval"
        result["assumptions"].extend(
            [
                "J/psi proton transition ratio, t profile and mass continuum, not a channel-specific measurement",
                "delta_pd(V) = delta_pd(J/psi) + delta_el(V) - delta_el(J/psi) locally at the common dissociation W0",
                "The double-log terms propagate the same proton-transition ratio to every W",
            ]
        )
        if channel.pdg == 100443:
            result["assumptions"].append(
                'The psi(2S) dissociative slope in arXiv:hep-ex/0205107 is not fitted here'
            )
    result["forward_ub_per_GeV2"] = forward
    result["elastic_forward_ub_per_GeV2"] = elastic
    result["photoprod_diss"] = [channel.pdg] + [
        round(float(value), 4) for value in (W0, forward / elastic, delta, b, n, 0.0808, mass_max)
    ]
    return result


# Convert one measurement into a machine-readable mapping
def measurement(value: model.Measurement) -> dict[str, object]:
    return {
        "central": value.central,
        "uncertainty_label": value.uncertainty_label,
        "uncertainty_up": value.uncertainty_up,
        "uncertainty_down": value.uncertainty_down,
        "syst_up": value.syst_up,
        "syst_down": value.syst_down,
    }


# Compute the direct steering fragment for one HERA coupling
def production(channel: model.HERAChannel, coupling: float) -> dict[str, object]:
    models: dict[str, object] = {}
    rounded = round(coupling, 9)
    for theory, first_pdg, second_pdg in channel.production_channels:
        key = f"[{first_pdg},{second_pdg}]"
        if theory not in {"MP", "XP", "GP"}:
            raise ValueError(f"unsupported HERA production model: {theory}")
        block: dict[str, object] = {"basis": "helicity", "helicity": [[-1, 0, rounded, 0.0]]}
        if theory in {"MP", "XP"}:
            block["Lambda"] = 1.0
        if theory == "MP":
            block["polarization"] = {"mode": "none"}
        models[theory] = {key: block}
    return models


# Build channel payload containing inputs, derived values and card values
def channel_output(channel: model.HERAChannel, derived: model.DerivedChannel) -> dict[str, object]:
    coupling = derived.gamma_pomeron_vector_coupling_per_gev
    payload = {
        "key": channel.key,
        "name": channel.name,
        "pdg": channel.pdg,
        "reference": {
            "citation": channel.source,
            "url": channel.source_url,
            "normalization": channel.normalization_origin,
            "trajectory": channel.trajectory_origin,
            "assumptions": list(channel.assumptions),
        },
        "inputs": {
            "W0_GeV": channel.w0_gev,
            "double_log": channel.dlog,
            "tmax_GeV2": channel.tmax_gev2,
            "forward_dsigma_dt_ub_per_GeV2": measurement(channel.forward_dsigma_dt_ub_per_gev2),
            "B_cross_section_GeV_minus2": measurement(channel.b0_per_gev2),
            "profile_power": channel.profile_power,
            "B_vector_cross_section_GeV_minus2": (channel.vector_cross_section_slope_per_gev2),
            "normalization_scale": channel.normalization_scale,
            "alpha0": measurement(channel.alpha0),
            "alpha_prime_GeV_minus2": measurement(channel.alpha_prime_per_gev2),
        },
        "derived": {
            "g_gammaPV_GeV_minus1": coupling,
            "card_couplings_GeV_minus1": {theory: round(coupling, 9) for theory, _, _ in channel.production_channels},
            "g_stat_up": derived.coupling_stat_up_per_gev,
            "g_stat_down": derived.coupling_stat_down_per_gev,
            "g_syst_up": derived.coupling_syst_up_per_gev,
            "g_syst_down": derived.coupling_syst_down_per_gev,
            "amplitude_abs_at_W0_t0": derived.reference_amplitude_abs,
            "reconstructed_forward_dsigma_dt_ub_per_GeV2": (derived.reconstructed_forward_dsigma_dt_ub_per_gev2),
            "exponential_sigma_ub": derived.exponential_sigma_ub,
            "matched_sigma_ub": derived.matched_sigma_ub,
            "integral_convention": "actual squared physical SOFT residue times the derived vector profile over the input t interval",
            "uncertainty_convention": "input errors scaled by the fixed normalization_scale with the SOFT proton residue held fixed",
            "forward_W_exponent": derived.forward_w_exponent,
            "dB_dlnW_GeV_minus2": derived.shrinkage_per_gev2,
            "B_vector_cross_section_GeV_minus2": (derived.vector_cross_section_slope_per_gev2),
            "B_factorized_cross_section_GeV_minus2": (derived.factorized_cross_slope_per_gev2),
        },
        "cards": {
            "photoprod": [
                channel.pdg,
                round(channel.w0_gev, 4),
                round(channel.vector_cross_section_slope_per_gev2, 4),
                round(channel.alpha0.central, 4),
                round(channel.alpha_prime_per_gev2.central, 4),
            ],
            "models": production(channel, coupling),
        },
    }
    if channel.dlog is not None:
        payload["cards"]["photoprod_dlog"] = [channel.pdg, *[round(value, 4) for value in channel.dlog]]
    if channel.fit is not None:
        payload["fit"] = channel.fit
    diss = dissociation(channel, derived)
    if diss is not None:
        payload["proton_dissociation"] = diss
        payload["cards"]["photoprod_diss"] = diss["photoprod_diss"]
        if diss.get("photoprod_diss_dlog"):
            payload["cards"]["photoprod_diss_dlog"] = diss["photoprod_diss_dlog"]
    return payload


# Build the combined JSON/card payload for selected channels
def output(
    selected: list[tuple[model.HERAChannel, model.DerivedChannel]], proton: model.ProtonVertex
) -> dict[str, object]:
    return {
        "amplitude_convention": {
            "A": ("g_gammaPV*[g_Ppp*Fp(t)]*exp(B_vector*t/2)*eta(t)*s*(s/W0^2)^(alpha(t)-1)*D(W)"),
            "dsigma_dt": "abs(A)^2/(16*pi*s^2)",
            "double_log": "D=exp[c/2*(sqrt(u)-sqrt(u0)-ln(W/W0)/sqrt(u0))], u=ln(W^2/scale2), D=1 when absent",
            "double_log_phase": "alpha_local(t)=alpha(t)+c/4*(1/sqrt(u)-1/sqrt(u0)), local derivative dispersion approximation",
            "energy_power": "4*(alpha0-1) is the local forward cross-section W power at W0",
            "signature_modulus": "abs(eta)=1 for rotating_t0 or rotating",
            "normalization": (
                'SOFT fixes g_Ppp, HERA fixes g_gammaPV'
            ),
            "proton_residue_role": "fixed bare physical SOFT proton end",
            "transition_residue_role": "HERA forward normalization",
            "Pomeron_beam_residue_t0_GeV_minus1": proton.beam_residue_per_gev,
            "proton_form_factor": proton.form_factor,
            "proton_form_factor_parameters": list(proton.parameters),
            "proton_amplitude_slope_GeV_minus2": (proton.amplitude_slope_per_gev2),
            "tune_source": proton.source,
            "GeV_minus2_to_microbarn": model.gev2_to_microbarn(),
        },
        "channels": [channel_output(channel, derived) for channel, derived in selected],
    }


# Compute only the compact steering-card rows from a full payload
def cards(payload: dict[str, object]) -> dict[str, object]:
    channels = payload["channels"]
    return {
        "PARAM_REGGE": {
            "photoprod": [channel["cards"]["photoprod"] for channel in channels],
            "photoprod_dlog": [
                channel["cards"]["photoprod_dlog"] for channel in channels if "photoprod_dlog" in channel["cards"]
            ],
            "photoprod_diss_dlog": [
                row for channel in channels for row in channel["cards"].get("photoprod_diss_dlog", [])
            ],
            "photoprod_diss": [
                channel["cards"]["photoprod_diss"] for channel in channels if "photoprod_diss" in channel["cards"]
            ],
        },
        "RES_channel_rows": {channel["name"]: channel["cards"]["models"] for channel in channels},
    }
