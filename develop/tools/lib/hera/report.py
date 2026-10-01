# Automated HERA fit tables and differential comparison PDF
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import json
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from .. import common
from . import couplings, model
from . import fit as hera_fit


# Write numerical results and physical predictions without applying fitted experimental nuisance shifts
def write(payload, cards, directory, proton, controls):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    for name, values in (("fit.json", payload), ("cards.json", cards)):
        (directory / name).write_text((common.dumps(values) if name == "cards.json" else json.dumps(values, indent=2)) + "\n")
    with PdfPages(directory / "fit_comparisons.pdf") as pdf:
        for channel in payload["channels"]:
            energy_overview(pdf, channel, proton, directory, controls)
            fit = channel.get("fit", {})
            spectra = fit.get("spectra", [])
            if "coordinates" in fit:
                spectra = [dict(fit, spectrum=channel["name"])]
            for row in spectra:
                coordinates = row.get("coordinates", {})
                bins = np.asarray(row.get("binedges", coordinates.get("binedges")))
                x = np.asarray(coordinates["t"]) if "t" in coordinates else bins.mean(axis=1)
                measured, prediction = np.asarray(row["data"]), np.asarray(row["prediction"])
                energy = np.asarray(coordinates.get("W", np.zeros(len(x))))
                for w in np.unique(energy):
                    use = np.isclose(energy, w)
                    fig, (ax, ratio) = plt.subplots(2, 1, sharex=True, gridspec_kw={"height_ratios": [3, 1]})
                    error = np.asarray(row.get("errors", np.zeros(len(x))))
                    ax.errorbar(x[use], measured[use], yerr=error[use], fmt="o", label="HERA")
                    ax.plot(x[use], prediction[use], label="Fixed proton residue fit")
                    ax.set_yscale("log")
                    unit = "nb" if channel["pdg"] == 443 else r"$\mu$b"
                    ax.set_ylabel(r"$d\sigma/d|t|$ [" + unit + r"/GeV$^2$]")
                    ax.set_title(row["spectrum"] + (f", W={w:.1f} GeV" if w > 0 else ""))
                    ax.legend()
                    ratio.errorbar(x[use], prediction[use] / measured[use],
                                   yerr=prediction[use] * error[use] / measured[use]**2, fmt="o")
                    ratio.axhline(1, color="black", linewidth=0.8)
                    ratio.set_ylabel("Fit / data")
                    ratio.set_xlabel(r"$|t|$ [GeV$^2$]")
                    fig.tight_layout()
                    pdf.savefig(fig)
                    plt.close(fig)
            energy = fit.get("energy")
            if energy is None and "W" in fit:
                energy = dict(fit, errors=(np.asarray(fit["errors_up"]) + fit["errors_down"]))
                energy["errors"] = np.asarray(energy["errors"]) / 2
            if energy is not None:
                comparison(pdf, channel["name"] + " HERA", energy["W"], energy, r"$W$ [GeV]", r"$\sigma$ [$\mu$b]")
            if "lhcb" in fit:
                row = fit["lhcb"]
                unit = "nb" if channel["pdg"] == 443 else "pb"
                comparison(pdf, channel["name"] + " LHCb", np.mean(row["binedges"], axis=1), row,
                           r"$y$", r"$d\sigma/dy$ [" + unit + "]")


# Plot independent integrated HERA measurements or screened LHCb rapidity spectra
def comparison(pdf, title, x, row, xlabel, ylabel):
    measured, prediction = np.asarray(row["data"]), np.asarray(row["prediction"])
    errors = np.asarray(row["errors"])
    fig, (ax, ratio) = plt.subplots(2, 1, sharex=True, gridspec_kw={"height_ratios": [3, 1]})
    ax.errorbar(x, measured, yerr=errors, fmt="o", label="Measurement")
    ax.plot(x, prediction, label="Common photoproduction fit")
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.legend()
    ratio.errorbar(x, prediction / measured, yerr=prediction * errors / measured**2, fmt="o")
    ratio.axhline(1, color="black", linewidth=0.8)
    ratio.set_ylabel("Fit / data")
    ratio.set_xlabel(xlabel)
    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


# Integrate the fitted elastic amplitude with the same proton vertex and transfer limits as the data
def energy_prediction(channel, proton, energy, tmax):
    inputs = channel["inputs"]
    logw = np.log(np.asarray(energy) / inputs["W0_GeV"])
    slope = inputs["B_vector_cross_section_GeV_minus2"] + 4 * inputs["alpha_prime_GeV_minus2"]["central"] * logw
    integral = np.array([model.profile_moments(proton, value, tmax)[0] for value in slope])
    prediction = inputs["forward_dsigma_dt_ub_per_GeV2"]["central"] * np.exp(4 * (inputs["alpha0"]["central"] - 1) * logw) * integral
    if inputs["double_log"] is not None:
        prediction *= model.dlog_factor(np.asarray(energy), inputs["W0_GeV"], *inputs["double_log"])
    return prediction


# Compare integrated HERA measurements with the common amplitude over the LHCb photon-energy range
# LHCb constrains two photon directions in pp, so its measured points remain on the rapidity panel
def energy_overview(pdf, channel, proton, directory, controls):
    fit = channel.get("fit", {})
    energy = fit.get("energy", fit if "W" in fit else None)
    measurements = []
    scale, unit = (1000, "nb") if channel["pdg"] == 443 else (1e6, "pb")
    if channel["pdg"] == 443:
        for name, w, data, covariance, sources, tmax in hera_fit.jpsi_energy():
            _, chol, relative, _ = hera_fit.jpsi_errors([data], [covariance], [sources])
            errors = np.sqrt(np.sum(chol**2, axis=1) + data**2 * np.sum(relative**2, axis=1))
            measurements.append((name, w, data, errors, tmax))
    elif energy is not None:
        errors = energy.get("errors", [energy.get("errors_down"), energy.get("errors_up")])
        measurements.append(("HERA", np.asarray(energy["W"]), np.asarray(energy["data"]) * scale,
                             np.asarray(errors) * scale, channel["inputs"]["tmax_GeV2"]))
    if not measurements:
        return
    lhcb = fit.get("lhcb")
    fig, axes = plt.subplots(1, 2 if lhcb else 1, figsize=(11, 5.8) if lhcb else (7, 5.8), squeeze=False)
    ax = axes[0, 0]
    low = min(np.min(row[1]) for row in measurements)
    high = max(np.max(row[1]) for row in measurements)
    if lhcb:
        bins = np.asarray(lhcb["binedges"])
        for sign, label in ((-1, r"LHCb $W_-$ range"), (1, r"LHCb $W_+$ range")):
            w = np.sqrt(lhcb["mass_GeV"] * np.asarray(lhcb["energies_GeV"])[:, None, None] * np.exp(sign * bins))
            low, high = min(low, w.min()), max(high, w.max())
            ax.axvspan(w.min(), w.max(), color="C2", alpha=0.10, label=label)
    grid = np.geomspace(low, high, controls["energy_points"])
    for name, w, data, errors, _ in measurements:
        ax.errorbar(w, data, yerr=errors, fmt="o", markersize=4, capsize=2, label=name)
    for index, tmax in enumerate(dict.fromkeys(row[-1] for row in measurements)):
        cut = r"all $|t|$" if tmax is None else rf"$|t|<{tmax:g}$ GeV$^2$"
        ax.plot(grid, scale * energy_prediction(channel, proton, grid, tmax), color="black",
                linestyle="-" if index == 0 else "--", label="Fit, " + cut)
    ax.set(xscale="log", yscale="log", xlabel=r"$W_{\gamma p}$ [GeV]", ylabel=rf"$\sigma(\gamma p\to Vp)$ [{unit}]")
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.15)
    if lhcb:
        rapidity = axes[0, 1]
        center = bins.mean(axis=1)
        rapidity.errorbar(center, lhcb["data"], xerr=(bins[:, 1] - bins[:, 0]) / 2,
                          yerr=lhcb["errors"], fmt="o", capsize=2, label="LHCb")
        order = np.argsort(lhcb["rapidity"])
        rapidity.plot(np.asarray(lhcb["rapidity"])[order], np.asarray(lhcb["curve"])[order], color="black", label="Screened fit")
        beam = "/".join(f"{value / 1000:g}" for value in lhcb["energies_GeV"])
        rapidity.set(xlabel=r"$y_V$", ylabel=rf"$d\sigma(pp\to pVp)/dy$ [{unit}]", title=rf"$\sqrt{{s}}={beam}$ TeV")
        rapidity.legend(fontsize=8)
        rapidity.grid(True, alpha=0.15)
    inputs = channel["inputs"]
    parameters = (rf"$W_0={inputs['W0_GeV']:g}$ GeV, "
                  rf"$d\sigma/dt|_0={scale * inputs['forward_dsigma_dt_ub_per_GeV2']['central']:.4f}$ {unit}/GeV$^2$, "
                  rf"$B_V={inputs['B_vector_cross_section_GeV_minus2']:.4f}$ GeV$^{{-2}}$" + "\n" +
                  rf"$\alpha_0={inputs['alpha0']['central']:.4f}$, "
                  rf"$\alpha'={inputs['alpha_prime_GeV_minus2']['central']:.4f}$ GeV$^{{-2}}$")
    if "photoprod_dlog" in channel["cards"]:
        _, scale2, coefficient = channel["cards"]["photoprod_dlog"]
        parameters += rf", $c={coefficient:.4f}$, $\mu^2={scale2:.4f}$ GeV$^2$"
    parameters += rf", $\chi^2/n_{{\rm dof}}={fit['chi2']:.2f}/{fit['ndf']}$"
    note = "HERA integrated projections are validation data" if channel["pdg"] == 443 else "Common HERA fit"
    if lhcb:
        note += r". LHCb constrains both $W_\pm^2=M_V\sqrt{s}\,e^{\pm y}$"
    fig.suptitle(channel["name"] + (": HERA and LHCb" if lhcb else ": HERA"))
    fig.text(0.5, 0.105, parameters, ha="center", fontsize=10)
    fig.text(0.5, 0.035, note, ha="center", fontsize=9)
    fig.tight_layout(rect=(0, 0.20, 1, 0.96))
    pdf.savefig(fig)
    path = directory / (channel["key"] + "_energy")
    fig.savefig(path.with_suffix(".pdf"))
    fig.savefig(path.with_suffix(".png"), dpi=160)
    plt.close(fig)


# Print compact parameter, coupling, card and reference tables
def print_tables(
    selected: list[tuple[model.HERAChannel, model.DerivedChannel]],
    proton: model.ProtonVertex,
) -> None:
    print("Vector meson photoproduction from HERA data")
    print(f"  SOFT input: {proton.source}")
    print("  A(s,t) = g_gammaPV [g_Ppp Fp(t)] exp(B_vector t / 2)")
    print("           x eta(t) s (s/W0^2)^(alpha(t)-1)")
    print("  dsigma/dt = |A|^2 / (16 pi s^2)")
    print("  |eta(t)| = 1 for the photoproduction rotating phase")
    print('  g_Ppp includes the SOFT Good Walker projection')
    print("  HERA forward data determine g_gammaPV without refitting g_Ppp")
    print(
        f"  g_Ppp(0) = {proton.beam_residue_per_gev:.4f} GeV^-1, Fp = {proton.form_factor}"
    )
    print(f"  d ln(Fp)/dt|0 = {proton.amplitude_slope_per_gev2:.4f} GeV^-2")
    print(f"  1 GeV^-2 = {model.gev2_to_microbarn():.9f} ub\n")

    _print_inputs(selected)
    _print_derived(selected)
    _print_channels(selected)


# Print measured inputs and uncertainties
def _print_inputs(selected):
    parameter_rows = []
    uncertainty_rows = []
    for channel, _ in selected:
        parameter_rows.append(
            [
                channel.name,
                channel.pdg,
                f"{channel.w0_gev:.1f}",
                "infinity" if channel.tmax_gev2 is None else f"{channel.tmax_gev2:.4f}",
                f"{channel.forward_dsigma_dt_ub_per_gev2.central:.6g}",
                f"{channel.b0_per_gev2.central:.4g}",
                f"{channel.alpha0.central:.5g}",
                f"{channel.alpha_prime_per_gev2.central:.4g}",
            ]
        )
        for label, measurement in (
            ("dsigma/dt(0) [ub/GeV2]", channel.forward_dsigma_dt_ub_per_gev2),
            ("B_cross [GeV^-2]", channel.b0_per_gev2),
            ("alpha(0)", channel.alpha0),
            ("alpha' [GeV^-2]", channel.alpha_prime_per_gev2),
        ):
            uncertainty_rows.append(
                [
                    channel.name,
                    label,
                    f"{measurement.central:.6g}",
                    measurement.uncertainty_label,
                    (f"+{measurement.uncertainty_up:.6g}/-{measurement.uncertainty_down:.6g}"),
                    f"+{measurement.syst_up:.6g}/-{measurement.syst_down:.6g}",
                ]
            )
    common.print_table(
        "HERA normalization and trajectory inputs",
        [
            "meson",
            "PDG",
            "W0 [GeV]",
            "|t| max [GeV^2]",
            "dsigma/dt(0) [ub/GeV2]",
            "B_cross [GeV^-2]",
            "alpha(0)",
            "alpha' [GeV^-2]",
        ],
        parameter_rows,
        right_align={1, 2, 3, 4, 5, 6, 7},
    )
    common.print_table(
        "\nInput uncertainties",
        ["meson", "observable", "central", "type", "first", "syst"],
        uncertainty_rows,
        right_align={2, 4, 5},
        group_by=0,
    )


# Print fitted couplings and cross sections
def _print_derived(selected):
    derived_rows = []
    coupling_error_rows = []
    for channel, derived in selected:
        derived_rows.append(
            [
                channel.name,
                f"{derived.gamma_pomeron_vector_coupling_per_gev:.9f}",
                f"{derived.reconstructed_forward_dsigma_dt_ub_per_gev2:.6g}",
                f"{derived.reference_amplitude_abs:.4f}",
                f"{derived.matched_sigma_ub:.6g}",
                f"{derived.vector_cross_section_slope_per_gev2:.4f}",
                f"{derived.factorized_cross_slope_per_gev2:.4f}",
                f"{derived.forward_w_exponent:.4f}",
                f"{derived.shrinkage_per_gev2:.4f}",
            ]
        )
        coupling_error_rows.append(
            [
                channel.name,
                f"+{derived.coupling_stat_up_per_gev:.9f}/-{derived.coupling_stat_down_per_gev:.9f}",
                f"+{derived.coupling_syst_up_per_gev:.9f}/-{derived.coupling_syst_down_per_gev:.9f}",
            ]
        )
    common.print_table(
        "\nDerived couplings and HERA reference integrals",
        [
            "meson",
            "g_gammaPV [GeV^-1]",
            "used dsigma/dt(0) [ub/GeV^2]",
            "|A(W0,0)|",
            "sigma_ref [ub]",
            "B_vector [GeV^-2]",
            "B_model(0) [GeV^-2]",
            "forward W power",
            "dB/dlnW [GeV^-2]",
        ],
        derived_rows,
        right_align={1, 2, 3, 4, 5, 6, 7, 8},
    )
    print('  sigma_ref integrates |Fp(t)|^2 exp(B_vector*t) over the input interval')
    print("  used dsigma/dt(0) includes the normalization scale listed below")
    print("  B_model(0) = 2 d ln(Fp)/dt|0 + B_vector at W0 for the SOFT proton profile")
    print('  B_vector fits the measured finite-|t| spectrum')
    common.print_table(
        "\nCoupling uncertainties [GeV^-1]", ["meson", "stat", "syst"], coupling_error_rows, right_align={1, 2}
    )
    print("  Input errors use the same fixed normalization scale as the central value")
    print('  SOFT and model uncertainties excluded')


# Print production rows and data references
def _print_channels(selected):
    card_rows = []
    for channel, derived in selected:
        for theory, first_pdg, second_pdg in channel.production_channels:
            card_rows.append(
                [
                    channel.name,
                    theory,
                    first_pdg,
                    second_pdg,
                    f"{derived.gamma_pomeron_vector_coupling_per_gev:.9f}",
                    "0.0",
                ]
            )
    common.print_table(
        "\nProduction-channel rows",
        ["meson", "model", "PDG1", "PDG2", "magnitude", "phase"],
        card_rows,
        right_align={2, 3, 4, 5},
        group_by=0,
    )

    reference_rows = [[channel.name, channel.source, channel.source_url] for channel, _ in selected]
    common.print_table("\nReferences", ["meson", "data", "URL"], reference_rows)
    trajectory_rows = [[channel.name, channel.trajectory_origin] for channel, _ in selected]
    common.print_table("\nTrajectory choices", ["meson", "origin"], trajectory_rows)
    normalization_rows = [
        [
            channel.name,
            channel.normalization_origin,
            f"{channel.normalization_scale:.4f}",
        ]
        for channel, _ in selected
    ]
    common.print_table("\nNormalization notes", ["meson", "input", "scale"], normalization_rows, right_align={2})
    for channel, _ in selected:
        if channel.assumptions:
            print(f"\n{channel.name} assumptions:")
            for assumption in channel.assumptions:
                print("  " + assumption)
    for channel, derived in selected:
        diss = couplings.dissociation(channel, derived)
        if diss is not None:
            print(f"\n{channel.name} proton dissociation: {diss['source_url']}")
            print("  [PDG, W0, forward pd/el ratio, delta, b, n, epsilon, M_max]")
            print("  " + common.dumps(diss["photoprod_diss"]))
            for assumption in diss["assumptions"]:
                print("  " + assumption)

