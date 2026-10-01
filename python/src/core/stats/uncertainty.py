# Weighted statistical uncertainties and covariance propagation
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import copy

import numba
import numpy as np
from core.numerics import array


# Construct one source-aware experimental uncertainty
def uncertainty_source(
    name: str,
    values,
    *,
    category: str,
    correlation: str,
    scope: str | None = None,
    down=None,
    shift=None,
    effect: str = "additive",
    provenance: str | None = None,
) -> dict:
    if category not in {"statistical", "systematic", "combined"}:
        raise ValueError(f"Unknown uncertainty category '{category}'")
    if correlation not in {"uncorrelated", "collective", "covariance"}:
        raise ValueError(f"Unknown uncertainty correlation '{correlation}'")
    if effect not in {"additive", "multiplicative"}:
        raise ValueError(f"Unknown uncertainty effect '{effect}'")

    source = {"name": str(name), "category": category, "correlation": correlation, "effect": effect}
    if correlation == "covariance":
        covariance = np.asarray(values, dtype=float)
        if covariance.ndim != 2 or covariance.shape[0] != covariance.shape[1]:
            raise ValueError("Covariance uncertainty must be a square matrix")
        source["covariance"] = covariance
    else:
        source["up"] = np.abs(np.asarray(values, dtype=float))
        source["down"] = np.abs(np.asarray(values if down is None else down, dtype=float))
        if correlation == "collective":
            source["shift"] = 0.5 * (source["up"] + source["down"]) if shift is None else np.asarray(shift, dtype=float)
    if scope is not None:
        source["scope"] = str(scope)
    if provenance is not None:
        source["provenance"] = str(provenance)
    return source


# Sum squared cluster weights per bin, treating entries within a cluster as fully correlated
def cluster_weight_squares(indices, weights, clusters, bin_count):
    indices, weights, clusters = map(np.asarray, (indices, weights, clusters))
    if indices.ndim != 1 or indices.shape != weights.shape or indices.shape != clusters.shape:
        raise ValueError("Inconsistent cluster weight arrays")
    valid = (indices >= 0) & (indices < bin_count)
    pairs, inverse = np.unique(np.column_stack((indices[valid], clusters[valid])), axis=0, return_inverse=True)
    totals = np.bincount(inverse, weights=weights[valid])
    return np.bincount(pairs[:, 0].astype(int), weights=totals**2, minlength=bin_count)


# Compute a scaled Euclidean norm without squaring the input magnitude directly
@numba.njit
def stable_l2_norm(values: np.ndarray) -> float:
    scale = 0.0
    scaled_sum = 0.0
    for value in values:
        magnitude = abs(value)
        if magnitude == 0.0:
            continue
        if not np.isfinite(magnitude):
            return magnitude
        if scale < magnitude:
            ratio = scale / magnitude
            scaled_sum = 1.0 + scaled_sum * ratio * ratio
            scale = magnitude
        else:
            ratio = magnitude / scale
            scaled_sum += ratio * ratio
    return 0.0 if scale == 0.0 else scale * np.sqrt(scaled_sum)


# Compute the norm of event errors without losing Torch derivatives
def norm(values):
    if array.namespace(values) is np:
        return stable_l2_norm(np.asarray(values))
    scale = values.detach().abs().max().clamp_min(np.finfo(float).tiny)
    return scale * array.sqrt(((values / scale)**2).sum())


# Compute the sum of squares using the stable Euclidean norm
def stable_sum_squares(values: np.ndarray) -> float:
    norm = stable_l2_norm(np.asarray(values, dtype=np.float64).ravel())
    return float(norm * norm)


# Compute a scaled floating-point sum without overflowing intermediate additions
def stable_sum(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=np.float64).ravel()
    if len(values) == 0:
        return 0.0
    scale = float(np.max(np.abs(values)))
    if scale == 0.0:
        return 0.0
    if not np.isfinite(scale):
        return float(np.sum(values))
    return float(scale * np.sum(values / scale))


# Compute stable weighted statistical errors for fixed one-dimensional bins
@numba.njit
def binned_weight_errors(values: np.ndarray, bins: np.ndarray, weights: np.ndarray) -> np.ndarray:
    bin_count = len(bins) - 1
    scales = np.zeros(bin_count)
    scaled_sums = np.zeros(bin_count)

    for index in range(len(values)):
        value = values[index]
        if not np.isfinite(value) or value < bins[0] or value > bins[-1]:
            continue
        bin_index = np.searchsorted(bins, value, side="right") - 1
        if bin_index == bin_count and value == bins[-1]:
            bin_index = bin_count - 1
        if bin_index < 0 or bin_index >= bin_count:
            continue

        magnitude = abs(weights[index])
        if magnitude == 0.0:
            continue
        if not np.isfinite(magnitude):
            scales[bin_index] = magnitude
            scaled_sums[bin_index] = 1.0
            continue
        if not np.isfinite(scales[bin_index]):
            continue
        if scales[bin_index] < magnitude:
            ratio = scales[bin_index] / magnitude
            scaled_sums[bin_index] = 1.0 + scaled_sums[bin_index] * ratio * ratio
            scales[bin_index] = magnitude
        else:
            ratio = magnitude / scales[bin_index]
            scaled_sums[bin_index] += ratio * ratio

    return scales * np.sqrt(scaled_sums)


# Project a covariance onto one weighted sum with roundoff-safe cancellation
def covariance_projection_error(covariance, weights):
    covariance = np.asarray(covariance, dtype=float)
    weights = np.asarray(weights, dtype=float)
    if covariance.shape != (len(weights), len(weights)):
        raise ValueError("hist: covariance projection arrays have incompatible dimensions")
    terms = np.outer(weights, weights) * covariance
    variance = float(np.sum(terms))
    tolerance = 64.0 * np.finfo(float).eps * float(np.sum(np.abs(terms)))
    if variance < -tolerance:
        raise ValueError("hist: covariance projection has negative variance")
    if abs(variance) <= tolerance:
        return 0.0
    return np.sqrt(variance)


# Compute a stable square root of one positive-semidefinite covariance
def covariance_square_root(covariance: np.ndarray) -> np.ndarray:
    covariance = np.asarray(covariance, dtype=float)
    if covariance.ndim != 2 or covariance.shape[0] != covariance.shape[1]:
        raise ValueError(__name__ + ".covariance_square_root: covariance must be square")
    covariance = 0.5 * (covariance + covariance.T)
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    scale = max(float(np.max(np.abs(eigenvalues))), 1.0e-30)
    if float(np.min(eigenvalues)) < -1.0e-10 * scale:
        raise ValueError(__name__ + ".covariance_square_root: covariance is not positive semidefinite")
    return eigenvectors * np.sqrt(np.clip(eigenvalues, 0.0, None))


# Compute the largest MC to data error ratio over all resolved covariance modes
def covariance_error_ratio(mc_covariance, data_covariance):
    mc = np.asarray(mc_covariance, dtype=float)
    total = mc + np.asarray(data_covariance, dtype=float)
    if not total.size or not np.any(total):
        return None
    if not np.any(mc):
        return 0.0
    scale = float(np.max(np.abs(total)))
    root = covariance_square_root(total / scale)
    precision = np.linalg.pinv(root, rcond=1.0e-6)
    fraction = float(np.linalg.eigvalsh(precision @ (mc / scale) @ precision.T)[-1])
    if fraction >= 1.0 - 64.0 * np.finfo(float).eps:
        return None
    return float(np.sqrt(max(fraction, 0.0) / (1.0 - fraction)))


# Compute weighted Poisson integral errors using the generator exposure only as a scale
def sample_cross_section(weights, attempted, scale):
    if not np.isfinite(attempted) or attempted <= 0:
        raise ValueError("MC exposure must be finite and positive")
    factor = scale / attempted
    return factor * stable_sum(weights), abs(factor) * stable_l2_norm(np.asarray(weights))


# Compute marginal variances or the full covariance after histogram normalization
def histogram(counts, errors, scale, *, valid, denominator=None, full=False):
    counts = array.asarray(counts)
    xp = array.namespace(counts)
    valid = array.asarray(valid, like=counts, dtype=bool)
    counts = xp.where(valid, counts, 0.0)
    errors = xp.where(valid, array.asarray(errors, like=counts), 0.0)
    scale = array.asarray(scale, like=counts)
    variance = (errors * scale)**2
    result = xp.diag(variance) if full else variance
    if denominator is None:
        return result
    total, total_error = denominator
    total = array.asarray(total, like=counts)
    usable = xp.isfinite(total) & (total > 0.0)
    total = xp.where(usable, total, 1.0)
    fraction = counts / total
    error_norm = norm(errors)
    total_error = array.asarray(total_error, like=counts)
    extra = 1.0 - (error_norm / xp.where(total_error > 0.0, total_error, 1.0))**2
    shift = fraction * scale * total_error * array.sqrt(extra)
    if full:
        # Form the normalization response before squaring to preserve its null integral mode
        response = (xp.diag(xp.ones_like(counts)) - fraction[:, None]) * (xp.ones_like(counts) * scale)[:, None] * errors[None, :]
        result = response @ response.T + xp.outer(shift, shift)
    else:
        # Sum other-bin variances without subtracting a dominant bin or allocating a matrix
        squared = (errors / xp.where(error_norm > 0.0, error_norm, 1.0))**2
        zero = xp.zeros_like(squared[:1])
        before = xp.concatenate((zero, xp.cumsum(squared, 0)[:-1]))
        after = xp.concatenate((xp.flip(xp.cumsum(xp.flip(squared, (0,)), 0), (0,))[1:], zero))
        result = ((1.0 - fraction) * scale * errors)**2 + (fraction * scale * error_norm)**2 * (before + after) + shift**2
    return xp.where(usable, result, 0.0)


# Report actual MC occupancy only within the measured data acceptance
def mc_statistics(histogram, valid=None):
    measured = np.asarray(histogram.valid, dtype=bool)
    if valid is not None:
        measured = measured & np.asarray(valid, dtype=bool)
    counts, errors = (array.to_numpy(value) for value in (histogram.counts, histogram.errs))
    ratio = np.divide(counts, errors, out=np.zeros_like(counts, dtype=float), where=errors > 0)
    effective = ratio**2
    entries = histogram.entries
    empty = None if entries is None else np.flatnonzero(measured & (np.asarray(entries) == 0)).tolist()
    return {"mc_empty_bin_count": None if empty is None else len(empty), "mc_empty_bins": empty,
            "mc_valid_bin_count": int(measured.sum()),
            "mc_min_neff": float(effective[measured].min()) if np.any(measured) else None}


# Record per-observable MC precision against the fixed measured-bin masks
def comparison_statistics(results):
    records = []
    for i, dataset in enumerate(results["data"]):
        if dataset is None:
            continue
        for j, subset in enumerate(dataset):
            for name, data in subset.items():
                mc = results["mc"][i][j][name]
                stats = mc_statistics(mc["hdata"], valid=data["hdata"].valid)
                records.append({"dataset": i, "subset": j, "observable": name, **stats})
    return records


# Add occupancies while retaining the distinction between data and event histograms
def merge_entries(first, second):
    if first is None and second is None:
        return None
    if first is None or second is None:
        raise ValueError("Cannot combine known and unknown MC occupancy")
    return first + second


# Build a shared-event second moment without allocating an event by bin matrix
def cross_moment(first, second, shape):
    ids1, bins1, weights1 = first
    ids2, bins2, weights2 = second
    _, left, right = np.intersect1d(ids1, ids2, assume_unique=True, return_indices=True)
    weights1 = array.asarray(weights1, like=weights2)
    weights2 = array.asarray(weights2, like=weights1)
    index = bins1[left] * shape[1] + bins2[right]
    xp = array.namespace(weights1)
    scales = [max(float(np.max(np.abs(array.to_numpy(w)), initial=0.0)), np.finfo(float).tiny) for w in (weights1, weights2)]
    values = (weights1[left] / scales[0]) * (weights2[right] / scales[1])
    if xp is np:
        moment = np.bincount(index, weights=values, minlength=shape[0]*shape[1]).reshape(shape)
    else:
        moment = values.new_zeros(shape[0]*shape[1]).index_add(0, array.asarray(index, like=values, dtype=int), values).reshape(shape)
    return moment, scales


# Join independent chunks retaining event identities for cross-observable covariance
def merge_events(first, second):
    merged = dict(first)
    for source, value in second.items():
        if source in merged:
            if np.intersect1d(merged[source][0], value[0]).size:
                raise ValueError("Cannot merge the same MC events as independent samples")
            value = tuple(array.concatenate([a, b]) for a, b in zip(merged[source], value, strict=True))
        merged[source] = value
    return merged


# Scale event contributions when adding independent differential cross sections
def scale_events(events, scale):
    result = {}
    for source, (ids, bins, weights) in events.items():
        factor = array.asarray(scale, like=weights)
        result[source] = (ids, bins, weights * (factor if factor.ndim == 0 else factor[bins]))
    return result


# Propagate full covariance between histograms sharing sampled events
def joint_covariance(histograms, selections):
    histograms = [h if h.stat_reference is None else h.stat_reference for h in histograms]
    reference = next((h.counts for h in histograms if array.namespace(h.counts) is not np), histograms[0].counts)
    xp = array.namespace(reference)
    blocks = []
    for first, selected in zip(histograms, selections, strict=True):
        row = []
        for second, other in zip(histograms, selections, strict=True):
            if first is second:
                covariance = first.covariance_scaled
            else:
                covariance = array.asarray(np.zeros((len(first.counts), len(second.counts))), like=first.counts)
                shared = first.mc_events.keys() & second.mc_events.keys()
                for source in shared:
                    moment, scales = cross_moment(first.mc_events[source], second.mc_events[source], covariance.shape)
                    covariance = covariance + transform_cross(moment, first, second, scales)
            row.append(array.asarray(covariance, like=reference)[selected][:, other])
        blocks.append(xp.concatenate(row, axis=1))
    return xp.concatenate(blocks, axis=0)


# Apply each histogram normalization Jacobian to a cross moment
def transform_cross(covariance, first, second, scales):
    for histogram, magnitude in zip((first, second), scales, strict=True):
        xp = array.namespace(covariance)
        valid = array.asarray(histogram.valid, like=covariance, dtype=bool)
        covariance = xp.where(valid[:, None], covariance, 0.0)
        scale = histogram.density_binscale() if histogram.density else histogram.binscale
        if histogram.density and histogram.density_uncertainty == "shape":
            if histogram.density_denominator_counts is not None:
                raise ValueError("Shared-event covariance requires the histogram's own normalization")
            counts = xp.where(valid, array.asarray(histogram.counts, like=covariance), 0.0)
            total = xp.where(valid, counts, 0.0).sum()
            covariance = covariance - (counts / xp.where(total > 0, total, 1.0))[:, None] * covariance.sum(axis=0)
        scale = array.asarray(scale, like=covariance) * magnitude
        covariance = (covariance * (scale if scale.ndim == 0 else scale[:, None])).T
    return covariance


# Compute the symmetric covariance represented by one uncertainty source
def source_covariance(source: dict) -> np.ndarray:
    if source["correlation"] == "covariance":
        return np.asarray(source["covariance"], dtype=float)
    sigma = 0.5 * (np.abs(np.asarray(source["up"], dtype=float)) + np.abs(np.asarray(source["down"], dtype=float)))
    if source["correlation"] == "uncorrelated":
        return np.diag(sigma**2)
    shift = np.asarray(source.get("shift", sigma), dtype=float)
    return np.outer(shift, shift)


# Approximate a signed offset pair by a Gaussian covariance and retain its asymmetric envelope
def offset_source(name, up, down, provenance):
    positive = np.maximum(np.maximum(up, down), 0.0)
    negative = np.abs(np.minimum(np.minimum(up, down), 0.0))
    magnitude = np.sqrt(0.5 * (positive**2 + negative**2))
    direction = np.sign(up - down)
    direction[direction == 0.0] = 1.0
    shift = direction * magnitude
    return {
        "name": name,
        "category": "systematic",
        "correlation": "covariance",
        "effect": "additive",
        "covariance": np.outer(shift, shift),
        "offsets": np.column_stack((up, down)),
        "up": positive,
        "down": negative,
        "provenance": provenance,
    }


# Subtract known independent components from a published marginal error
def residual_uncertainty(total, components) -> np.ndarray:
    total = np.abs(np.asarray(total, dtype=float))
    known_variance = np.zeros_like(total)
    for component in components:
        known_variance += np.asarray(component, dtype=float) ** 2
    variance = total**2 - known_variance
    scale = max(float(np.max(total**2)), float(np.max(known_variance)), np.finfo(float).tiny)
    tolerance = 64.0 * np.finfo(float).eps * scale
    if np.any(variance < -tolerance):
        raise ValueError("Known uncertainty components exceed the published marginal")
    return np.sqrt(np.clip(variance, 0.0, None))


# Separate known common sources while preserving the quoted asymmetric marginal errors
def split_uncertainty(sources, common, source=None):
    sources = copy.deepcopy(sources)
    if source is not None:
        marginal = next(item for item in sources if item["name"] == source)
        if marginal["correlation"] != "uncorrelated":
            raise ValueError("Common sources require an independent unresolved marginal")
        for side in ("up", "down"):
            marginal[side] = residual_uncertainty(marginal[side], [item[side] for item in common])
        marginal["provenance"] += ". Remaining bin correlations are unavailable, diagonal covariance is assumed"
    return sources + copy.deepcopy(common)


# Fill weighted bins with stable errors for NumPy and differentiable Torch arrays
def binned(x, bins, weights, indices):
    if array.namespace(weights) is np:
        selected = indices >= 0
        counts = np.bincount(indices[selected], weights=np.asarray(weights, dtype=float)[selected], minlength=len(bins) - 1)
        errs = binned_weight_errors(np.asarray(x, dtype=float), bins, np.asarray(weights, dtype=float))
    else:
        selected = indices >= 0
        index = array.asarray(indices[selected], like=weights, dtype=int)
        values = weights[selected]
        counts = weights.new_zeros(len(bins) - 1).index_add(0, index, values)
        scale = counts.new_zeros(counts.shape).scatter_reduce(0, index, values.detach().abs(), reduce="amax")
        scale = scale.clamp_min(np.finfo(float).tiny)
        variance = counts.new_zeros(counts.shape).index_add(0, index, (values / scale[index])**2)
        errs = scale * array.sqrt(variance)
    return counts, errs


# Compute the Jacobian for converting bin heights to unit-integral density
def density_jacobian(values: np.ndarray, binwidth: np.ndarray, valid: np.ndarray | None = None, *, shape=True) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    binwidth = np.asarray(binwidth, dtype=float)
    mask = np.ones_like(values, dtype=bool) if valid is None else np.asarray(valid, dtype=bool)
    if values.ndim != 1 or binwidth.shape != values.shape or mask.shape != values.shape:
        raise ValueError("Density covariance arrays must be aligned vectors")
    if not np.all(np.isfinite(values[mask])) or not np.all(np.isfinite(binwidth[mask])):
        raise ValueError("Density covariance requires finite values and bin widths")
    if np.any(binwidth[mask] <= 0.0):
        raise ValueError("Density covariance requires positive bin widths")
    active_values = np.zeros_like(values)
    weighted_width = np.zeros_like(binwidth)
    active_values[mask] = values[mask]
    weighted_width[mask] = binwidth[mask]
    integral = float(np.sum(active_values * weighted_width))
    if not np.isfinite(integral) or integral <= 0.0:
        raise ValueError("Density covariance requires a positive finite integral")
    scale = np.diag(mask.astype(float)) / integral
    return scale - np.outer(active_values / integral, weighted_width / integral) if shape else scale


# Propagate source covariances and signed systematic offsets through a linear transformation
def linear_source(source, transform, full=False):
    item = copy.deepcopy(source)
    correlation = source["correlation"]
    if correlation == "covariance" or (correlation == "uncorrelated" and full):
        covariance = source_covariance(source)
        item["covariance"] = transform @ covariance @ transform.T
        item["correlation"] = "covariance"
        if correlation == "uncorrelated":
            item["mc_correlation_eligible"] = source["category"] == "statistical"
            for side in ("up", "down"):
                item[side] = np.sqrt(transform**2 @ np.asarray(source[side], dtype=float)**2)
        elif "offsets" not in source:
            # Propagate quoted margins with the same Gaussian correlation matrix
            sigma = np.sqrt(np.maximum(np.diag(covariance), 0.0))
            for side in ("up", "down"):
                if side in source:
                    scaled = transform * np.divide(source[side], sigma, out=np.zeros_like(sigma), where=sigma > 0.0)
                    item[side] = np.sqrt(np.maximum(np.sum((scaled @ covariance) * scaled, axis=1), 0.0))
    else:
        for side in ("up", "down"):
            errors = np.asarray(source[side], dtype=float)
            item[side] = np.sqrt(transform**2 @ errors**2) if correlation == "uncorrelated" else np.abs(transform @ errors)
        if correlation == "collective":
            item["shift"] = transform @ np.asarray(source["shift"], dtype=float)
    if "offsets" in source:
        offsets = item["offsets"] = transform @ source["offsets"]
        item["up"] = np.maximum(offsets.max(axis=1), 0.0)
        item["down"] = np.maximum(-offsets.min(axis=1), 0.0)
    return item


# Transform one reader uncertainty source to plotted histogram units
def transform_source(
    source: dict,
    *,
    values: np.ndarray,
    binwidth: np.ndarray,
    binscale,
    density: bool,
    density_uncertainty: str = "shape",
    valid: np.ndarray | None = None,
) -> dict:
    scale = np.asarray(binscale, dtype=float)
    if scale.ndim == 0:
        scale = np.full(len(values), float(scale), dtype=float)
    transform = np.diag(scale)
    if density:
        if density_uncertainty not in {"shape", "scaled"}:
            raise ValueError("Density uncertainty must be either 'scaled' or 'shape'")
        transform = density_jacobian(values, binwidth, valid=valid, shape=density_uncertainty == "shape")

    return linear_source(source, transform, full=density)


# Transform all sources attached to one ordinary reader histogram
def transform_sources(
    sources: list[dict],
    *,
    values: np.ndarray,
    binwidth: np.ndarray,
    binscale,
    density: bool,
    density_uncertainty: str = "shape",
    valid: np.ndarray | None = None,
) -> list[dict]:
    return [
        transform_source(
            source,
            values=np.asarray(values, dtype=float),
            binwidth=np.asarray(binwidth, dtype=float),
            binscale=binscale,
            density=density,
            density_uncertainty=density_uncertainty,
            valid=valid,
        )
        for source in sources
    ]


# Convert covariance and marginal uncertainty into one correlation matrix
def covariance_to_correlation(covariance: np.ndarray, uncertainty: np.ndarray) -> np.ndarray:
    covariance = np.asarray(covariance, dtype=float)
    uncertainty = np.asarray(uncertainty, dtype=float)
    denominator = np.outer(uncertainty, uncertainty)
    correlation = np.divide(covariance, denominator, out=np.zeros_like(covariance), where=denominator > 0.0)
    correlation = 0.5 * (correlation + correlation.T)
    return np.clip(correlation, -1.0, 1.0)


# Combine asymmetric independent error sources in quadrature
def orthogonal_error(data: dict) -> np.ndarray:
    total_error = np.zeros(len(next(iter(data.values()))))
    for errors in data.values():
        total_error = np.hypot(total_error, np.mean(np.abs(errors), axis=1))
    return total_error


# Attach covariance errors and retain the separate upper and lower marginal errors
def finalize_uncertainties(ds: dict, sources: list[dict]) -> dict:
    size = len(np.asarray(ds["y"]))
    total = np.zeros((size, size), dtype=float)
    category_covariance = {
        "statistical": np.zeros((size, size), dtype=float),
        "systematic": np.zeros((size, size), dtype=float),
        "combined": np.zeros((size, size), dtype=float),
    }
    for source in sources:
        covariance = source_covariance(source)
        if covariance.shape != (size, size):
            raise ValueError(
                f"Uncertainty source '{source['name']}' has shape {covariance.shape}, expected {(size, size)}"
            )
        total += covariance
        category_covariance[source["category"]] += covariance

    ds["uncertainties"] = sources
    ds["y_err_stat"] = np.sqrt(np.clip(np.diag(category_covariance["statistical"]), 0.0, None))
    ds["y_err_syst"] = np.sqrt(
        np.clip(np.diag(category_covariance["systematic"] + category_covariance["combined"]), 0.0, None)
    )
    ds["y_err"] = np.sqrt(np.clip(np.diag(total), 0.0, None))
    for side in ("up", "down"):
        variance = np.zeros(size)
        stat_variance = np.zeros(size)
        for source in sources:
            term = (
                np.diag(source_covariance(source))
                if side not in source
                else np.asarray(source[side], dtype=float) ** 2
            )
            variance += term
            if source["category"] == "statistical":
                stat_variance += term
        ds[f"y_err_{side}"] = np.sqrt(np.maximum(variance, 0.0))
        ds[f"y_err_stat_{side}"] = np.sqrt(np.maximum(stat_variance, 0.0))
    return ds


# Compute an uncorrelated ratio and propagated uncertainty
def cross_section_ratio(
    numerator: float,
    numerator_error: float,
    denominator: float,
    denominator_error: float,
) -> tuple[float | None, float | None]:
    if denominator == 0.0:
        return None, None
    ratio = numerator / denominator
    variance = (numerator_error / denominator) ** 2
    variance += (numerator * denominator_error / denominator**2) ** 2
    return ratio, np.sqrt(max(variance, 0.0))


# Compute the covariance represented by one plotted histogram record
def histogram_covariance(histogram, uncertainties=None):
    size = len(np.asarray(histogram.counts_scaled, dtype=float))
    if not uncertainties:
        stored = getattr(histogram, "covariance_scaled", None)
        covariance = (
            np.diag(np.asarray(histogram.errs_scaled, dtype=float) ** 2)
            if stored is None
            else np.asarray(stored, dtype=float)
        )
        if covariance.shape != (size, size):
            raise ValueError("iceplot: histogram covariance has incompatible dimensions")
        if not np.all(np.isfinite(covariance)):
            raise ValueError("iceplot: histogram covariance contains non-finite values")
        return covariance

    covariance = np.zeros((size, size), dtype=float)
    for source in uncertainties:
        source_matrix = source_covariance(source)
        if source_matrix.shape != covariance.shape:
            raise ValueError("iceplot: uncertainty source covariance has incompatible dimensions")
        covariance += source_matrix
    if not np.all(np.isfinite(covariance)):
        raise ValueError("iceplot: uncertainty source covariance contains non-finite values")
    return covariance




# Compute the combined header-mode cross section and uncertainty
def aggregate_header_xsection(chunks: list[dict], total_wsum: float) -> tuple[float, float]:
    files = {}
    for chunk in chunks:
        record = chunk["before"]
        group = files.setdefault(
            chunk["filename"],
            {
                "wsum": 0.0,
                "xsection_pb": record["xsection_pb"],
                "xsection_pb_err": record["xsection_pb_err"],
            },
        )
        group["wsum"] += record["wsum"]

    exposure = sum(
        record["wsum"] / record["xsection_pb"]
        for record in files.values()
        if record["xsection_pb"] != 0.0
    )
    if exposure == 0.0:
        return 0.0, 0.0

    xsection_pb = total_wsum / exposure
    relative_error = sum(
        abs((record["wsum"] / record["xsection_pb"]) / exposure)
        * abs(record["xsection_pb_err"] / record["xsection_pb"])
        for record in files.values()
        if record["xsection_pb"] != 0.0
    )
    return xsection_pb, abs(xsection_pb) * relative_error



# Compute the combined sample-mode cross section and uncertainty
def aggregate_sample_xsection(chunks: list[dict], wsum: float, wsum2: float) -> tuple[float, float]:
    attempted_events = sum(chunk["before"]["attempted_events"] for chunk in chunks)
    if attempted_events <= 0.0:
        return 0.0, 0.0

    rescale = next(
        (
            record["xsection_pb"] * record["attempted_events"] / record["wsum"]
            for chunk in chunks
            if (record := chunk["before"])["wsum"] != 0.0
        ),
        0.0,
    )
    return rescale * wsum / attempted_events, abs(rescale) * np.sqrt(wsum2) / attempted_events



# Compute binwise Monte Carlo precision over the selected histogram intervals
def histogram_mc_precision(counts, errors, mask):
    selected = np.asarray(mask, dtype=bool)
    values = np.asarray(counts, dtype=float)
    uncertainties = np.asarray(errors, dtype=float)
    nonzero = selected & (values != 0.0)
    zero_count = int(np.count_nonzero(selected & ~nonzero))
    maximum = None
    if np.any(nonzero):
        relative = uncertainties[nonzero] / np.abs(values[nonzero])
        maximum = float(np.max(relative))
        if not np.isfinite(maximum):
            raise ValueError("iceplot: non-finite relative MC uncertainty")
    return {
        "mc_max_rel_uncertainty": maximum,
        "mc_zero_bin_count": zero_count,
        "mc_valid_bin_count": int(np.count_nonzero(selected)),
    }
