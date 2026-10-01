# Weighted particle flow diagnostics for the Pandora driver
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import math
from pathlib import Path

import numpy as np

from core.io.serialize import load_json_file
from core.stats.hist import bin_indices, edge2centerbins
from core.stats.uncertainty import cluster_weight_squares

PF_CLASS_NAMES = ("photon", "electron", "muon", "charged_hadron", "neutral_hadron", "other")
PF_CONFUSION_CLASS_NAMES = ("charged_hadron", "neutral_hadron", "photon", "electron", "muon", "other")
PF_PLOT_CLASS_NAMES = tuple(name for name in PF_CLASS_NAMES if name != "other")
PF_CONFUSION_PLOT_CLASS_NAMES = tuple(name for name in PF_CONFUSION_CLASS_NAMES if name != "other")
PF_RESIDUAL_CLASS_NAMES = ("charged_hadron", "photon", "neutral_hadron", "electron", "muon")
PF_CONFUSION_ROW_NAMES = (*PF_RESIDUAL_CLASS_NAMES, "lost")
PF_CONFUSION_COLUMN_NAMES = (*PF_RESIDUAL_CLASS_NAMES, "fake")
PF_OBJECT_PERFORMANCE_AXES = {
    "energy": ("energy_bins_gev", "energy", r"$E$ [GeV]", True, True),
    "pt": ("pt_bins_gev", "pt", r"$p_{\mathrm{T}}$ [GeV]", True, True),
    "costheta": ("costheta_bins", "costheta", r"$\cos\theta$", False, False),
    "eta": ("eta_bins", "eta", r"$\eta$", False, False)}
PF_ENERGY_CONFUSION_COLORBAR_LABEL = (
    r"$C_{ij} = \left\langle E_{i\to j} / E_{\mathrm{gen}}\right\rangle_{\mathrm{events}}$")
PF_CONFUSION_RESIDUAL_COLORBAR_LABEL = r"$\Delta C_{ij} = C_{ij} - C^{\mathrm{ideal}}_{ij}$"


# Compute moments in the shortest interval containing the requested event weight
def interval_stats(values, fraction=0.90, weights=None):
    x = np.asarray(values, dtype=float).reshape(-1)
    w = np.ones(x.size) if weights is None else np.asarray(weights, dtype=float).reshape(-1)
    if x.shape != w.shape or not 0.0 < fraction <= 1.0 or np.any(~np.isfinite(w) | (w < 0.0)):
        raise ValueError("Invalid central interval weights or fraction")
    valid = np.isfinite(x) & (w > 0.0)
    order = np.argsort(x[valid], kind="stable")
    x, w = x[valid][order], w[valid][order]
    if not x.size:
        return dict(valid=False, n=0, mean=None, sigma=None, low=None, high=None, fraction=fraction)
    cumulative = np.r_[0.0, np.cumsum(w)]
    target = np.nextafter(cumulative[:-1] + fraction * cumulative[-1], -np.inf)
    ends = np.searchsorted(cumulative, target, side="left")
    starts = np.flatnonzero(ends <= x.size)
    start = starts[np.argmin(x[ends[starts] - 1] - x[starts])]
    end = np.searchsorted(x, x[ends[start] - 1], side="right")
    start = np.searchsorted(x, x[start], side="left")
    mean = float(np.average(x[start:end], weights=w[start:end]))
    sigma = float(np.sqrt(np.average((x[start:end] - mean) ** 2, weights=w[start:end])))
    return dict(valid=True, n=int(end - start), mean=mean, sigma=sigma,
                low=float(x[start]), high=float(x[end - 1]), fraction=float(fraction))


# Concatenate diagnostic particle fields and expand each event weight once per collection
def _profile_arrays(profiles, weights, *, event_ids=False):
    arrays = {key: np.concatenate([p[key] for p in profiles]) for key in profiles[0]}
    for key in profiles[0]:
        if key.endswith("_energy"):
            sizes = [len(p[key]) for p in profiles]
            arrays[key.removesuffix("energy") + "weight"] = np.repeat(weights, sizes)
            if event_ids:
                arrays[key.removesuffix("energy") + "event"] = np.repeat(np.arange(len(profiles)), sizes)
    return arrays


# Extend standard energy bins to cover the observed particle range
def _energy_bins(energies, edges, factor):
    bins = list(edges)
    positive = np.asarray(energies)[np.asarray(energies) > 0.0]
    if positive.size:
        while positive.min() < bins[0]:
            bins.insert(0, bins[0] / factor)
        while positive.max() >= bins[-1]:
            bins.append(bins[-1] * factor)
    return np.asarray(bins)


# Compute particle performance and class confusion with the objective's event weights
def plot_inputs(events, weights, *, settings=None):
    settings = plot_settings() if settings is None else settings
    binning = settings["binning"]
    classes = {name: _profile_arrays([e["object_profiles"][name] for e in events], weights, event_ids=True) for name in PF_CLASS_NAMES}
    bins = {f"{field}_bins_gev": _energy_bins(np.concatenate(
                [p[f"{side}_{field}"] for p in classes.values() for side in ("truth", "reco")]),
                binning[field], binning["extend_factor"])
            for field in ("energy", "pt")}
    performance = dict(classes=classes, costheta_bins=np.asarray(binning["costheta"]),
                       eta_bins=np.asarray(binning["eta"]), **bins)
    confusion = dict(class_names=PF_CLASS_NAMES, energy_bins_gev=bins["energy_bins_gev"],
                     **_profile_arrays([e["confusion_profile"] for e in events], weights, event_ids=True))
    return dict(events=events, object_performance=performance, class_confusion=confusion, plot_settings=settings)


# Compute binned rates with approximate Wilson intervals and independent event weights
def binomial_profile(energy, passed, bins, *, weights=None, events=None, positive_only=True, geometric_centers=True, totals=None):
    energy, passed, bins = np.asarray(energy), np.asarray(passed, dtype=bool), np.asarray(bins)
    weights = np.ones(energy.size) if weights is None else np.asarray(weights)
    if energy.shape != passed.shape or weights.shape != energy.shape or not np.all(np.diff(bins) > 0.0):
        raise ValueError("Inconsistent particle efficiency inputs")
    valid = np.isfinite(energy) & (energy > 0.0 if positive_only else True)
    numerator = np.histogram(energy[valid & passed], bins, weights=weights[valid & passed])[0]
    denominator, weight2 = ([np.histogram(energy[valid], bins, weights=w[valid])[0] for w in (weights, weights**2)]
                            if totals is None else (totals["denominator"], totals["weight2"]))
    if events is not None and totals is None:
        weight2 = cluster_weight_squares(bin_indices(energy[valid], bins), weights[valid], np.asarray(events)[valid], len(bins) - 1)
    rate, lower, upper = _weighted_wilson(numerator, denominator, weight2)
    centers = np.sqrt(bins[:-1] * bins[1:]) if geometric_centers else (bins[:-1] + bins[1:]) / 2.0
    return dict(bins=bins, centers=centers, rate=rate, yerr=np.maximum(np.vstack((rate - lower, upper - rate)), 0.0),
                numerator=numerator, denominator=denominator, weight2=weight2, nonzero=denominator > 0.0)


# Approximate Wilson intervals assuming full correlation within each event and independence between events
def _weighted_wilson(numerator, denominator, weight2):
    rate = np.divide(numerator, denominator, out=np.zeros_like(numerator, dtype=float), where=denominator > 0.0)
    rate = np.clip(rate, 0.0, 1.0)
    effective = np.divide(denominator**2, weight2, out=np.zeros_like(denominator, dtype=float), where=weight2 > 0.0)
    inverse = np.divide(1.0, effective, out=np.zeros_like(effective), where=effective > 0.0)
    scale = 1.0 + inverse
    center = (rate + inverse / 2.0) / scale
    half = np.sqrt(rate * (1.0 - rate) * inverse + inverse**2 / 4.0) / scale
    return rate, np.clip(center - half, 0.0, 1.0), np.clip(center + half, 0.0, 1.0)


# Compute the finite positive or half-open energy-range selection
def energy_mask(values: np.ndarray, low: float | None, high: float | None) -> np.ndarray:
    finite = np.isfinite(values)
    if low is None or high is None:
        return finite & (values > 0.0)
    return finite & (values >= low) & (values < high)


# Compute row-normalized and column-normalized class confusion matrices for one energy range
def confusion_matrices(confusion, *, low, high):
    names = tuple(confusion["class_names"])
    count = len(names)
    result = dict(row_names=(*names, "fake"), column_names=(*names, "lost"))
    for side, axis, prefix in (("truth", 1, "row"), ("reco", 0, "column")):
        classes = np.asarray(confusion[f"{side}_class"], dtype=int)
        weight = np.asarray(confusion.get(f"{side}_weight", np.ones(classes.size)))
        mask = energy_mask(np.asarray(confusion[f"{side}_energy"]), low, high)
        denominator, weight2 = [np.bincount(classes[mask], weights=w[mask], minlength=count) for w in (weight, weight**2)]
        if f"{side}_event" in confusion:
            weight2 = cluster_weight_squares(classes[mask], weight[mask], np.asarray(confusion[f"{side}_event"])[mask], count)
        mask = energy_mask(np.asarray(confusion[f"match_{side}_energy"]), low, high)
        weight = np.asarray(confusion.get(f"match_{side}_weight", np.ones(mask.size)))
        tc, rc = (np.asarray(confusion[f"match_{s}_class"], dtype=int)[mask] for s in ("truth", "reco"))
        numerator = np.bincount(tc * count + rc, weights=weight[mask], minlength=count**2).reshape(count, count)
        unmatched = np.maximum(denominator - numerator.sum(axis=axis), 0.0)
        numerator = np.pad(numerator, ((0, 1), (0, 1)))
        if side == "truth":
            numerator[:-1, -1] = unmatched
        else:
            numerator[-1, :-1] = unmatched
        denominator, weight2 = np.r_[denominator, 0.0], np.r_[weight2, 0.0]
        rate, lower, upper = _weighted_wilson(numerator, np.expand_dims(denominator, axis), np.expand_dims(weight2, axis))
        result.update({f"{side}_denominator": denominator, f"{prefix}_numerator": numerator,
                       f"{prefix}_rate": rate, f"{prefix}_lower": lower, f"{prefix}_upper": upper})
    return result


# Compute weighted quantiles and central spread of the energy residuals
def _resolution_stats(values, weights):
    selected = weights > 0
    values, weights = values[selected], weights[selected]
    q16, median, q84 = np.quantile(values, [0.16, 0.5, 0.84], weights=weights, method="inverted_cdf")
    central = interval_stats(values, weights=weights)
    return np.array([median, (q84 - q16) / 2, central["mean"], central["sigma"]])


# Estimate statistical errors by resampling independent events with their particle weights
def _resolution_errors(values, weights, events, resampling):
    from scipy.stats import bootstrap

    groups, inverse = np.unique(events, return_inverse=True)
    if groups.size < 2:
        return np.full(4, np.nan)

    # Keep all particles from the same event together in every bootstrap sample
    def statistic(sample):
        counts = np.bincount(sample, minlength=groups.size)
        return _resolution_stats(values, weights * counts[inverse])

    return bootstrap((np.arange(groups.size),), statistic, vectorized=False, batch=1, method="percentile",
                     n_resamples=resampling["resamples"], rng=np.random.default_rng(resampling["seed"])).standard_error


# Compute weighted energy residual moments for assigned truth particles in each bin
def resolution_profile(profile, bins, *, correct_only=False, resampling=None):
    energy, residual, weights = (profile["truth_" + key] for key in ("energy", "energy_residual", "weight"))
    valid = (energy > 0) & np.isfinite(residual) & (weights > 0) & profile["truth_assigned"]
    if correct_only:
        valid &= profile["truth_matched"]
    indices = bin_indices(energy, bins)
    robust = {key: np.full(len(bins) - 1, np.nan) for key in ("median", "q68", "mean90", "sigma90")}
    errors = {key: np.full(len(bins) - 1, np.nan) for key in robust}
    for i in np.unique(indices[valid & (indices >= 0)]):
        selected = valid & (indices == i)
        values, weight = residual[selected], weights[selected]
        for key, value in zip(robust, _resolution_stats(values, weight), strict=True):
            robust[key][i] = value
        if resampling is not None:
            error = _resolution_errors(values, weight, profile["truth_event"][selected], resampling)
            for key, value in zip(errors, error, strict=True):
                errors[key][i] = value
    return dict(centers=np.sqrt(bins[:-1] * bins[1:]), errors=errors, **robust)


# Compute event-level GEN and PFO composition arrays for class summary plots
def composition_profiles(events, *, low=None, high=None):
    names = PF_CONFUSION_PLOT_CLASS_NAMES
    output = dict(class_names=names)
    profiles = [[e["object_profiles"][name] for name in names] for e in events]
    for side, prefix in (("gen", "truth"), ("pfo", "reco")):
        values = [[np.asarray(p[f"{prefix}_momentum"])[
            energy_mask(np.asarray(p[f"{prefix}_energy"]), low, high)] for p in row] for row in profiles]
        shape = (len(events), len(names))
        counts = np.array([[v.size for v in row] for row in values], dtype=float).reshape(shape)
        momentum = np.array([[np.clip(v, 0, None).sum() for v in row] for row in values]).reshape(shape)
        total = momentum.sum(axis=1, keepdims=True)
        fractions = np.divide(momentum, total, out=np.zeros_like(momentum), where=total > 0)
        output[f"{side}_total_multiplicity"] = counts.sum(axis=1)
        for field, array in (("multiplicity", counts), ("momentum_fraction", fractions)):
            output[f"{side}_{field}"] = {name: array[:, i] for i, name in enumerate(names)}
    return output


# Resolve the event losses into signed particle contributions before the logarithm
def loss_particles(truth, reco, pairs, *, pid, counts, response, momentum, spec):
    ti, ri = pairs.T
    tc, rc = truth["class_index"], reco["class_index"]
    et, er = truth["energy"], reco["energy"]
    nt, nr = len(et), len(er)
    total, reco_total = et.sum(), max(er.sum(), spec["epsilon"])
    # Use the Gen energy for assigned PFOs so correct matches cancel in the same bin
    energy = er.copy()
    energy[ri] = et[ti]
    row = np.full(nt, len(PF_RESIDUAL_CLASS_NAMES), dtype=int)
    row[ti] = rc[ri]
    fake = np.ones(nr, dtype=bool)
    fake[ri] = False
    confusion = np.r_[et / total * (pid[row, tc] - pid[tc, tc]),
                      fake * er / reco_total * pid[rc, -1]] / spec["confusion_energy_fraction_target"]**2
    delta = (counts[1] - counts[0]) / spec["multiplicity_target"]**2
    terms = dict(energy=np.r_[et, energy], classes=np.r_[tc, rc], residual=np.r_[-et, er] / total,
                 confusion=confusion, multiplicity=np.r_[-delta[tc], delta[rc]])
    for name, values in (("response", response), ("momentum", momentum)):
        terms[name] = np.zeros(nt + nr)
        terms[name][ti] = values
    return terms


# Compute an exact signed decomposition of each loss by particle type and energy
def loss_data(detail, bins):
    spec, weights = detail["loss"], detail["weights"]
    mean = detail["calibration"]["mean"]
    calibration = (mean - 1) / spec["energy_bias_target"]**2
    resolution = (math.hypot(mean, spec["epsilon"]) * spec["energy_resolution_target"])**2
    shape = (len(PF_RESIDUAL_CLASS_NAMES), len(bins) - 1)
    contributions = {name: np.zeros(shape) for name in detail["loss_components"]}
    for event, weight in zip(detail["events"], weights / weights.sum(), strict=True):
        profile = event["loss_profile"]
        indices = bin_indices(profile["energy"], bins)
        valid = indices >= 0
        indices = profile["classes"][valid] * shape[1] + indices[valid]
        terms = dict(calibration=calibration * profile["residual"],
                     resolution=(event["visible_response"] - mean) * profile["residual"] / resolution,
                     **{name: profile[name] for name in ("confusion", "response", "momentum", "multiplicity")})
        for name, values in terms.items():
            contributions[name] += weight * np.bincount(indices, weights=values[valid], minlength=np.prod(shape)).reshape(shape)
    return dict(bins=bins, classes=PF_RESIDUAL_CLASS_NAMES, energy=contributions,
                particle={name: values.sum(axis=1) for name, values in contributions.items()})


# Subdivide histogram bins uniformly in linear or logarithmic coordinates
def split_bins(edges, subdivisions, *, log_midpoints=False):
    edges = np.asarray(edges)
    space = np.geomspace if log_midpoints else np.linspace
    return np.r_[space(edges[:-1], edges[1:], subdivisions, endpoint=False, axis=1).ravel(), edges[-1]]


# Read and validate plotting controls before preparing any diagnostics
def plot_settings(path=None):
    settings = load_json_file(Path(__file__).with_name("settings.json") if path is None else path)
    bins, hist, resolution = (settings[key] for key in ("binning", "histogram", "resolution"))
    for axis, (_, _, _, positive, _) in PF_OBJECT_PERFORMANCE_AXES.items():
        edges = np.asarray(bins[axis], dtype=float)
        if (edges.ndim != 1 or edges.size < 2 or not np.isfinite(edges).all()
                or np.any(np.diff(edges) <= 0) or (positive and edges[0] <= 0)):
            raise ValueError(f"Invalid Pandora {axis} bin edges")
    integers = (bins["subdivisions"], bins["matching_bins"], hist["min_bins"], hist["max_bins"], hist["multiplier"])
    if any(type(value) is not int or value < 1 for value in integers) or hist["min_bins"] > hist["max_bins"]:
        raise ValueError("Pandora bin counts must be positive integers with min_bins <= max_bins")
    if not np.isfinite(bins["extend_factor"]) or bins["extend_factor"] <= 1:
        raise ValueError("Pandora energy bin extension factor must exceed one")
    if (type(resolution["resamples"]) is not int or resolution["resamples"] < 2
            or type(resolution["seed"]) is not int or resolution["seed"] < 0):
        raise ValueError("Resolution resampling requires at least two replicas and a nonnegative integer seed")
    return settings


# Compute efficiency, purity and matching diagnostics once for one particle variable
def _particle_axis(profile, bins, field, positive, log_x, max_angle, matching_bins):
    rates = {side: binomial_profile(profile[f"{side}_{field}"], profile[f"{side}_matched"], bins,
             weights=profile[f"{side}_weight"], events=profile[f"{side}_event"], positive_only=positive, geometric_centers=log_x)
             for side in ("truth", "reco")}
    values, weights = profile[f"truth_{field}"], profile["truth_weight"]
    assigned, correct, candidates = (profile[f"truth_{key}"] for key in ("assigned", "matched", "candidates"))
    outcomes = [rates["truth"], *[binomial_profile(values, mask, bins, weights=weights,
                 positive_only=positive, geometric_centers=log_x, totals=rates["truth"]) for mask in
                 (assigned & ~correct, candidates == 0, ~assigned & (candidates > 0))]]
    matching = {}
    for key in ("candidates", "angle", "dr"):
        data = profile[f"truth_{key}"]
        selected = np.isfinite(values) & np.isfinite(data) & (weights > 0) & ((values > 0) if positive else True)
        matching[key] = None
        if np.any(selected):
            ybins = (np.arange(-0.5, int(data[selected].max()) + 1.5, 1) if key == "candidates"
                     else np.linspace(0, max_angle if key == "angle" else max(max_angle, data[selected].max()), matching_bins + 1))
            matching[key] = np.histogram2d(values[selected], data[selected], bins=(bins, ybins), weights=weights[selected])
    return dict(rates=rates, outcomes=outcomes, matching=matching)


# Compute composition histograms and momentum fraction moments for one energy interval
def _composition_data(events, weights, low, high):
    composition = composition_profiles(events, low=low, high=high)
    names, multiplicity, fractions = composition["class_names"], {}, {}
    for name in (*names, "total"):
        values = [composition[f"{side}_total_multiplicity"] if name == "total"
                  else composition[f"{side}_multiplicity"][name] for side in ("gen", "pfo")]
        finite = np.concatenate([v[np.isfinite(v)] for v in values])
        high_count = max(math.ceil(finite.max()) if finite.size else 0, 1)
        bins = np.arange(-0.5, high_count + 1.5, 1.0)
        multiplicity[name] = dict(bins=bins, counts=[np.histogram(v, bins, weights=weights)[0] for v in values],
                                  joint=np.histogram2d(*values, bins=(bins, bins), weights=weights)[0])
    for side in ("gen", "pfo"):
        values = np.array([composition[f"{side}_momentum_fraction"][name] for name in names])
        mean = np.average(values, axis=1, weights=weights)
        fractions[side] = (mean, np.sqrt(np.average((values - mean[:, None])**2, axis=1, weights=weights)))
    return dict(class_names=names, multiplicity=multiplicity, fractions=fractions)


# Prepare all numerical plot inputs once, leaving rendering free of statistical calculations
def prepare_plots(details):
    detail = details["pflow"]
    settings = detail["plot_settings"] if "plot_settings" in detail else plot_settings()
    hist, binning = settings["histogram"], settings["binning"]
    data = dict(detail=detail, settings=settings, histograms={}, particles={}, slices=[])
    for prefix in ("", "visible_total_energy_", "truth_energy_"):
        values = detail[prefix + "values"]
        bins = min(hist["max_bins"], max(hist["min_bins"], int(math.sqrt(len(values))))) * hist["multiplier"]
        counts, edges = np.histogram(values, bins=bins, weights=detail["weights"])
        data["histograms"][prefix] = dict(counts=counts, bins=edges, n=len(values),
                                         mode=edge2centerbins(edges)[np.argmax(counts)])
    if "events" not in detail:
        return data
    performance, max_angle = detail["object_performance"], detail["loss"]["matching_max_angle_rad"]
    data["bins"] = {axis: split_bins(performance[key], binning["subdivisions"], log_midpoints=log_x)
                    for axis, (key, _, _, _, log_x) in PF_OBJECT_PERFORMANCE_AXES.items()}
    data["loss_contributions"] = loss_data(detail, data["bins"]["energy"])
    for name in PF_PLOT_CLASS_NAMES:
        profile = performance["classes"][name]
        data["particles"][name] = dict(
            resolution=[resolution_profile(profile, data["bins"]["energy"], correct_only=correct, resampling=settings["resolution"])
                        for correct in (False, True)],
            axes={axis: _particle_axis(profile, data["bins"][axis], field, positive, log_x, max_angle, binning["matching_bins"])
                  for axis, (_, field, _, positive, log_x) in PF_OBJECT_PERFORMANCE_AXES.items()})
    bins = detail["class_confusion"]["energy_bins_gev"]
    for low, high in [(None, None), *zip(bins[:-1], bins[1:], strict=True)]:
        matrices = confusion_matrices(detail["class_confusion"], low=low, high=high)
        if np.any(matrices["truth_denominator"]) or np.any(matrices["reco_denominator"]):
            data["slices"].append(dict(low=low, high=high, matrices=matrices,
                composition=_composition_data(detail["events"], detail["weights"], low, high)))
    return data
