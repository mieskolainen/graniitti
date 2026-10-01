# Histogram data, arithmetic and binning
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import copy

import numba
import numpy as np
from core.numerics import array
from core.stats import uncertainty
from core.stats.uncertainty import stable_l2_norm


# Convert ordered disjoint low and high pairs to one edge array
def doublet2linear(binedges: list) -> np.ndarray:
    edges = np.asarray(binedges, dtype=float)
    if edges.ndim != 2 or edges.shape[1] != 2 or len(edges) == 0:
        raise ValueError("HEPData bin edges must be a nonempty array of low/high pairs")
    if not np.isfinite(edges).all() or np.any(edges[:, 1] <= edges[:, 0]):
        raise ValueError("HEPData bin edges must be finite and increasing")
    gaps = edges[1:, 0] - edges[:-1, 1]
    # Allow only rounding from coordinate arithmetic, independent of physical units
    scale = np.maximum(np.abs(edges[:-1, 1]), np.abs(edges[1:, 0]))
    tolerance = 4 * np.spacing(scale)
    if np.any(gaps < -tolerance):
        raise ValueError("HEPData bins overlap or are not ordered")
    bins = [edges[0, 0]]
    for index, (low, high) in enumerate(edges):
        if index and gaps[index - 1] > tolerance[index - 1]:
            bins.append(low)
        bins.append(high)
    if np.any(np.diff(bins) <= 0.0):
        raise ValueError("HEPData bin edges are not increasing")
    return np.asarray(bins)


# Compute bin widths from HEPData low/high bin pairs
@numba.njit
def doublet2binwidth(binedges: list) -> np.ndarray:
    binwidth = np.zeros(len(binedges))
    for i in range(len(binedges)):
        binwidth[i] = binedges[i][1] - binedges[i][0]
    return binwidth


# Compute bin widths from linear bin edges
@numba.njit
def bins2binwidth(bins: np.ndarray) -> np.ndarray:
    return np.diff(bins)


# Compute bin centers from linear bin edges
@numba.njit
def edge2centerbins(bins: np.ndarray) -> np.ndarray:
    return (bins[1:] + bins[:-1]) / 2


# Infer linear bin edges from bin centers
@numba.njit
def center2edgebins(cbins: np.ndarray) -> np.ndarray:
    if cbins.ndim != 1 or not np.isfinite(cbins).all() or np.any(cbins[1:] <= cbins[:-1]):
        raise ValueError("Bin centers must be finite and strictly increasing")
    if len(cbins) == 0:
        return np.zeros(1)
    if len(cbins) == 1:
        return np.array([cbins[0] - 0.5, cbins[0] + 0.5])

    edges = np.zeros(len(cbins) + 1)
    edges[1:-1] = cbins[:-1] + (cbins[1:] - cbins[:-1]) / 2
    edges[0] = cbins[0] - (cbins[1] - cbins[0]) / 2
    edges[-1] = cbins[-1] + (cbins[-1] - cbins[-2]) / 2
    return edges


# Rebin validated histogram arrays over explicit source-bin groups
@numba.njit
def _rebin_histogram_groups(
    bin_edges: np.ndarray,
    bin_contents: np.ndarray,
    bin_errors: np.ndarray,
    groups: np.ndarray,
    differential: bool = True,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    count = len(groups) - 1
    new_edges = bin_edges[groups]
    new_contents = np.zeros(count)
    new_errors = np.zeros(count)

    if differential:
        widths = bins2binwidth(bin_edges)
        bin_contents = bin_contents * widths
        bin_errors = bin_errors * widths

    for index in range(count):
        start = groups[index]
        stop = groups[index + 1]
        new_contents[index] = np.sum(bin_contents[start:stop])
        new_errors[index] = stable_l2_norm(bin_errors[start:stop])

    if differential:
        widths = bins2binwidth(new_edges)
        new_contents /= widths
        new_errors /= widths
    return new_edges, new_contents, new_errors


# Validate and rebin a differential or integral histogram over explicit groups
def rebin_histogram_groups(
    bin_edges: np.ndarray,
    bin_contents: np.ndarray,
    bin_errors: np.ndarray,
    groups: np.ndarray,
    differential: bool = True,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    bin_edges = np.asarray(bin_edges, dtype=float)
    bin_contents = np.asarray(bin_contents, dtype=float)
    bin_errors = np.asarray(bin_errors, dtype=float)
    raw_groups = np.asarray(groups, dtype=float)
    arrays = (bin_edges, bin_contents, bin_errors, raw_groups)
    if not isinstance(differential, (bool, np.bool_)):
        raise TypeError("rebin_histogram_groups: differential must be boolean")
    if any(values.ndim != 1 for values in arrays):
        raise ValueError("rebin_histogram_groups: inputs must be one-dimensional")
    if len(bin_edges) != len(bin_contents) + 1:
        raise ValueError("rebin_histogram_groups: bin edges and contents are not aligned")
    if len(bin_contents) != len(bin_errors):
        raise ValueError("rebin_histogram_groups: bin contents and errors are not aligned")
    if len(bin_contents) == 0:
        raise ValueError("rebin_histogram_groups: histogram must contain at least one bin")
    if not np.all(np.isfinite(bin_edges)) or np.any(np.diff(bin_edges) <= 0.0):
        raise ValueError("rebin_histogram_groups: bin edges must be finite and increasing")
    if not np.all(np.isfinite(bin_contents)):
        raise ValueError("rebin_histogram_groups: bin contents must be finite")
    if not np.all(np.isfinite(bin_errors)) or np.any(bin_errors < 0.0):
        raise ValueError("rebin_histogram_groups: bin errors must be finite and non-negative")
    if not np.all(np.isfinite(raw_groups)) or np.any(raw_groups != np.floor(raw_groups)):
        raise ValueError("rebin_histogram_groups: groups must contain integer indices")
    groups = raw_groups.astype(np.int64)
    if len(groups) < 2 or groups[0] != 0 or groups[-1] != len(bin_contents):
        raise ValueError("rebin_histogram_groups: groups must span every source bin")
    if np.any(np.diff(groups) <= 0):
        raise ValueError("rebin_histogram_groups: groups must be strictly increasing")
    return _rebin_histogram_groups(bin_edges, bin_contents, bin_errors, groups, differential)


# Partition bins without crossing measured and unmeasured intervals
def rebin_partition(valid: np.ndarray, factor: int) -> np.ndarray:
    """Compute source bin groups without mixing measured and unmeasured intervals"""
    valid = np.asarray(valid, dtype=np.bool_)
    if valid.ndim != 1 or len(valid) == 0:
        raise ValueError("rebin_partition: validity mask must be a non-empty vector")
    if isinstance(factor, (bool, np.bool_)) or not isinstance(factor, (int, np.integer)):
        raise TypeError("rebin_partition: rebin factor must be an integer")
    if factor <= 0:
        raise ValueError("rebin_partition: rebin factor must be positive")
    transitions = np.flatnonzero(valid[1:] != valid[:-1]) + 1
    runs = np.concatenate((np.asarray([0]), transitions, np.asarray([len(valid)])))
    groups = [0]
    # Split at every validity change before applying the maximum group size
    for start, stop in zip(runs[:-1], runs[1:], strict=True):
        groups.extend(range(int(start) + factor, int(stop), factor))
        groups.append(int(stop))
    return np.asarray(groups, dtype=np.int64)


# Rebin a differential or integral histogram by an integer factor
def rebin_histogram(
    bin_edges: np.ndarray,
    bin_contents: np.ndarray,
    bin_errors: np.ndarray,
    rebin_factor: int,
    differential: bool = True,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    groups = rebin_partition(np.ones(len(bin_contents), dtype=bool), rebin_factor)
    return rebin_histogram_groups(bin_edges, bin_contents, bin_errors, groups, differential)


def hist_to_density(counts, errs, bins):
    """Normalize to unit integral density function over the visible histogram range"""
    return hist_to_density_fullspace(counts, errs, bins, totalweight=counts.sum())


# Normalize weighted bins without fixing their total event weight in the derivative
def hist_to_density_fullspace(counts, errs, bins, totalweight):
    scale = compute_binscale(bins2binwidth(bins), totalweight, 1.0)
    return counts * scale, errs * scale


def normalize_bin_heights(values, errors, binwidth, valid=None):
    """Normalize bin heights to a unit integral density, returning zeros for empty inputs."""
    values = np.asarray(values, dtype=float)
    errors = np.asarray(errors, dtype=float)
    binwidth = np.asarray(binwidth, dtype=float)
    mask = np.ones_like(values, dtype=bool) if valid is None else np.asarray(valid, dtype=bool)
    if errors.shape != values.shape or binwidth.shape != values.shape or mask.shape != values.shape:
        raise ValueError("hist: density arrays must be aligned vectors")
    if not np.all(np.isfinite(values[mask])) or not np.all(np.isfinite(binwidth[mask])):
        raise ValueError("hist: density requires finite values and bin widths")
    if np.any(binwidth[mask] <= 0.0):
        raise ValueError("hist: density requires positive bin widths")
    active_values = np.zeros_like(values)
    active_width = np.zeros_like(binwidth)
    active_values[mask] = values[mask]
    active_width[mask] = binwidth[mask]
    integral = np.sum(active_values * active_width)

    if not np.isfinite(integral) or integral <= 0:
        zeros = np.zeros_like(values, dtype=float)
        return zeros, zeros

    normalized_values = np.zeros_like(values)
    normalized_errors = np.zeros_like(errors)
    normalized_values[mask] = values[mask] / integral
    normalized_errors[mask] = errors[mask] / integral
    return normalized_values, normalized_errors


# Compute the bin scale with the proposal normalization included in the autograd graph
def compute_binscale(binwidth, totalweight, xsection_pb):
    total = array.asarray(totalweight)
    xp = array.namespace(total)
    width = array.asarray(binwidth, like=total)
    cross_section = array.asarray(xsection_pb, like=total)
    valid = xp.isfinite(total) & (total > 0.0)
    return xp.where(valid, cross_section / (width * xp.where(valid, total, 1.0)), 0.0)


class hobj:
    """Minimal histogram data object"""

    def __init__(
        self,
        counts=0,
        errs=0,
        bins=0,
        cbins=0,
        binscale=1.0,
        valid=None,
        density=False,
        density_uncertainty="scaled",
        density_denominator_counts=None,
        density_denominator_errs=None,
        entries=None,
        mc_events=None,
        stat_reference=None,
        differential=True,
    ):
        if density_uncertainty not in {"scaled", "shape"}:
            raise ValueError("hist: density uncertainty must be either 'scaled' or 'shape'")
        if (density_denominator_counts is None) != (density_denominator_errs is None):
            raise ValueError("hist: density denominator counts and errors must be provided together")
        self.stat_reference = stat_reference
        self.entries = entries
        self.mc_events = {} if mc_events is None else mc_events
        self.counts = counts
        self.errs = errs
        self.bins = bins
        self.cbins = cbins
        self.binscale = binscale
        self.density = bool(density)
        self.differential = bool(differential)
        self.density_uncertainty = density_uncertainty
        self.density_denominator_counts = density_denominator_counts
        self.density_denominator_errs = density_denominator_errs
        if valid is None:
            self.valid = np.ones(np.shape(counts), dtype=np.bool_)
        else:
            self.valid = np.asarray(valid, dtype=np.bool_)

    @property
    def is_empty(self):
        counts = np.asarray(self.counts_scaled, dtype=float)
        errors = np.asarray(self.errs_scaled, dtype=float)
        valid = np.asarray(self.valid, dtype=np.bool_)
        return not np.any(valid & ((counts != 0.0) | (errors != 0.0)))

    @property
    def binwidth(self):
        return bins2binwidth(self.bins)

    # Compute integration weights for densities or cross sections per bin
    @property
    def measure(self):
        return self.binwidth if self.differential or self.density else np.ones_like(self.binwidth)

    # Compute the visible range density scale after raw event counts have been fused
    def density_binscale(self):
        counts = self._density_denominator_payload()[0]
        xp = array.namespace(counts)
        valid = array.asarray(self.valid, like=counts, dtype=bool)
        total = xp.where(valid, counts, 0.0).sum()
        width = array.asarray(self.binwidth, like=counts)
        usable = xp.isfinite(total) & (total > 0.0)
        return xp.where(usable, 1.0 / (xp.where(usable, total, 1.0) * width), 0.0)

    # Compute the raw denominator histogram used by a unit density normalization
    def _density_denominator_payload(self):
        counts = self.counts if self.density_denominator_counts is None else self.density_denominator_counts
        errors = self.errs if self.density_denominator_errs is None else self.density_denominator_errs
        return array.asarray(counts, like=self.counts), array.asarray(errors, like=self.counts)

    # Compute scaled bin contents while retaining Torch parameter dependence
    @property
    def counts_scaled(self):
        counts = array.asarray(self.counts)
        xp = array.namespace(counts)
        scale = self.density_binscale() if self.density else array.asarray(self.binscale, like=counts)
        return xp.where(array.asarray(self.valid, like=counts, dtype=bool), counts * scale, 0.0)

    # Compute scaled marginal errors after the selected normalization transform
    @property
    def errs_scaled(self):
        if self.stat_reference is not None:
            return self.stat_reference.errs_scaled
        if not self.density or self.density_uncertainty == "scaled":
            return self.errs_scaled_fixed
        return array.sqrt(self._covariance())

    # Compute fixed denominator errors without the density normalization response
    @property
    def errs_scaled_fixed(self):
        errors = array.asarray(self.errs)
        xp = array.namespace(errors)
        scale = self.density_binscale() if self.density else array.asarray(self.binscale, like=errors)
        return xp.where(array.asarray(self.valid, like=errors, dtype=bool), errors * scale, 0.0)

    # Compute covariance with the same differentiable equations for NumPy and Torch histograms
    @property
    def covariance_scaled(self):
        return self._covariance(full=True)

    # Compute marginal variances by default without allocating a dense covariance
    def _covariance(self, *, full=False):
        if self.stat_reference is not None:
            return self.stat_reference._covariance(full=full)
        denominator = None
        scale = self.density_binscale() if self.density else self.binscale
        if self.density and self.density_uncertainty == "shape":
            counts, errors = self._density_denominator_payload()
            xp = array.namespace(counts)
            valid = array.asarray(self.valid, like=counts, dtype=bool)
            errors = xp.where(valid, errors, 0.0)
            denominator = (xp.where(valid, counts, 0.0).sum(), uncertainty.norm(errors))
        return uncertainty.histogram(self.counts, self.errs, scale, valid=self.valid, denominator=denominator, full=full)

    # Compute histogram integral
    # (with piece-wise differential measure and scale applied)
    def integral(self):
        dx = self.measure
        return np.sum(dx * self.counts_scaled)

    # Compute the integral error using the selected statistical reference
    def integral_error(self):
        dx = self.measure
        if not self.density or self.density_uncertainty == "scaled":
            return stable_l2_norm(dx * self.errs_scaled)
        return uncertainty.covariance_projection_error(self.covariance_scaled, dx)

    # Compute histogram count sums
    # (without differential measure or scale applied)
    def counts_total(self):
        return np.sum(self.counts)

    def errs_total(self):
        return stable_l2_norm(np.asarray(self.errs, dtype=float))

    # Check if this histogram carries normalization exposure for fusion
    def has_exposure(self):
        bins = np.asarray(self.bins)
        if bins.ndim == 0 or len(bins) < 2:
            return False

        binscale = array.to_numpy(self.binscale)
        return np.any(np.isfinite(binscale) & (binscale > 0))

    # Validate the numerical payload before combining histogram objects
    def validate_combination_payload(self, operation):
        bins = np.asarray(self.bins, dtype=float)
        counts = array.to_numpy(self.counts)
        errors = array.to_numpy(self.errs)
        valid = np.asarray(self.valid, dtype=np.bool_)
        centers = np.asarray(self.cbins, dtype=float)
        scale = array.to_numpy(self.binscale)
        denominator_counts, denominator_errors = self._density_denominator_payload()
        if bins.ndim != 1 or len(bins) < 2 or not np.all(np.diff(bins) > 0.0):
            raise ValueError(f"hist: {operation} requires strictly increasing one-dimensional bins")
        expected_shape = (len(bins) - 1,)
        if counts.shape != expected_shape or errors.shape != expected_shape:
            raise ValueError(f"hist: {operation} found incompatible bin-content shapes")
        if valid.shape != expected_shape:
            raise ValueError(f"hist: {operation} found an incompatible valid-bin mask")
        if centers.shape != expected_shape:
            raise ValueError(f"hist: {operation} found incompatible bin centers")
        if scale.ndim > 1 or (scale.ndim == 1 and scale.shape != expected_shape):
            raise ValueError(f"hist: {operation} found an incompatible bin scale")
        if denominator_counts.shape != expected_shape or denominator_errors.shape != expected_shape:
            raise ValueError(f"hist: {operation} found an incompatible density denominator")
        arrays = (bins, centers, counts, errors, scale, array.to_numpy(denominator_counts), array.to_numpy(denominator_errors))
        if any(not np.all(np.isfinite(values)) for values in arrays):
            raise ValueError(f"hist: {operation} requires finite histogram values")
        if np.any(errors < 0.0):
            raise ValueError(f"hist: {operation} requires nonnegative errors")
        if np.any(array.to_numpy(denominator_errors) < 0.0):
            raise ValueError(f"hist: {operation} requires nonnegative denominator errors")

    # Check that two histogram objects use identical bin edges
    def check_compatible_bins(self, other, operation):
        self_bins = np.asarray(self.bins)
        other_bins = np.asarray(other.bins)

        if self.differential != other.differential:
            raise ValueError(f"hist: {operation} requires matching differential conventions")
        if self_bins.shape != other_bins.shape or not np.array_equal(self_bins, other_bins):
            raise ValueError(f"hist: {operation} cannot combine different histogram bins")

    # Compute the combined inverse-exposure scale for two event chunks
    def combined_chunk_exposure_scale(self, other):
        first = array.asarray(self.binscale)
        second = array.asarray(other.binscale, like=first)
        xp = array.namespace(first)
        with np.errstate(divide="ignore", invalid="ignore"):
            inverse = 1 / first + 1 / second
            valid = xp.isfinite(inverse) & (inverse > 0.0)
            return xp.where(valid, 1 / xp.where(valid, inverse, 1.0), 0.0)

    # Fuse two statistically independent chunks from the same physics source
    def fuse_independent_chunk(self, other):
        operation = "independent chunk fusion"
        self.validate_combination_payload(operation)
        other.validate_combination_payload(operation)
        self.check_compatible_bins(other=other, operation=operation)
        if self.density != other.density:
            raise ValueError(f"hist: {operation} cannot mix density and cross-section histograms")
        if self.density:
            if self.density_uncertainty != other.density_uncertainty:
                raise ValueError(f"hist: {operation} requires matching density uncertainty policies")
            if self.density_denominator_counts is not None or other.density_denominator_counts is not None:
                raise ValueError(f"hist: {operation} does not accept an external density denominator")
            return hobj(
                self.counts + other.counts,
                array.hypot(self.errs, other.errs),
                self.bins,
                self.cbins,
                binscale=1.0,
                valid=self.valid & other.valid,
                density=True,
                differential=self.differential,
                density_uncertainty=self.density_uncertainty,
                entries=uncertainty.merge_entries(self.entries, other.entries),
                mc_events=uncertainty.merge_events(self.mc_events, other.mc_events),
            )
        if not self.has_exposure():
            return copy.deepcopy(other)
        if not other.has_exposure():
            return copy.deepcopy(self)

        counts = self.counts + other.counts
        errs = array.hypot(self.errs, other.errs)
        binscale = self.combined_chunk_exposure_scale(other=other)

        return hobj(
            counts,
            errs,
            self.bins,
            self.cbins,
            binscale,
            valid=self.valid & other.valid,
            density=self.density,
            differential=self.differential,
            density_uncertainty=self.density_uncertainty,
            entries=uncertainty.merge_entries(self.entries, other.entries),
            mc_events=uncertainty.merge_events(self.mc_events, other.mc_events),
        )

    # Sum two independent physics processes in scaled cross-section space
    def sum_independent_process(self, other):
        operation = "independent process summation"
        self.validate_combination_payload(operation)
        other.validate_combination_payload(operation)
        self.check_compatible_bins(other=other, operation=operation)
        if self.density or other.density:
            raise ValueError(f"hist: {operation} does not accept density histograms")

        counts = array.asarray(self.counts_scaled)
        counts = counts + array.asarray(other.counts_scaled, like=counts)
        errors = array.hypot(array.asarray(self.errs_scaled), array.asarray(other.errs_scaled))
        return hobj(
            counts=counts,
            errs=errors,
            bins=copy.deepcopy(self.bins),
            cbins=copy.deepcopy(self.cbins),
            binscale=1.0,
            valid=self.valid & other.valid,
            density=False,
            differential=self.differential,
            entries=uncertainty.merge_entries(self.entries, other.entries),
            mc_events=uncertainty.merge_events(
                uncertainty.scale_events(self.mc_events, self.binscale),
                uncertainty.scale_events(other.mc_events, other.binscale)),
        )

    # Fuse multiple independent chunks from one physics source
    @classmethod
    def fuse_independent_chunks(cls, histograms):
        if not histograms:
            raise ValueError("hist: independent chunk fusion requires at least one histogram")
        total = (copy.deepcopy(histograms[0]) if array.namespace(histograms[0].counts) is np
                 else copy.copy(histograms[0]))
        for histogram in histograms[1:]:
            total = total.fuse_independent_chunk(histogram)
        return total

    # Sum multiple statistically independent physics processes
    @classmethod
    def sum_independent_processes(cls, histograms):
        if not histograms:
            raise ValueError("hist: independent process summation requires at least one histogram")

        operation = "independent process summation"
        first = histograms[0]
        first.validate_combination_payload(operation)
        if first.density:
            raise ValueError(f"hist: {operation} does not accept density histograms")
        total = cls(
            counts=array.asarray(first.counts_scaled) + 0.0,
            errs=array.asarray(first.errs_scaled) + 0.0,
            bins=copy.deepcopy(first.bins),
            cbins=copy.deepcopy(first.cbins),
            binscale=1.0,
            valid=copy.deepcopy(first.valid),
            entries=copy.deepcopy(first.entries),
            mc_events=uncertainty.scale_events(first.mc_events, first.binscale),
            density=False,
            differential=first.differential,
        )
        for histogram in histograms[1:]:
            total = total.sum_independent_process(histogram)
        return total

    # Compute cumulative sums used to draw independent processes as a stack
    @classmethod
    def stack_independent_processes(cls, histograms):
        if not histograms:
            raise ValueError("hist: process stacking requires at least one histogram")
        cumulative = []
        total = None
        for histogram in histograms:
            total = (
                cls.sum_independent_processes([histogram])
                if total is None
                else total.sum_independent_process(histogram)
            )
            cumulative.append(copy.deepcopy(total))
        return cumulative


# Detach completed histogram and objective outputs for plotting and persistence
def numpy_output(value):
    if isinstance(value, hobj):
        result = copy.copy(value)
        result.__dict__ = {key: numpy_output(item) for key, item in vars(value).items()}
        return result
    if isinstance(value, dict):
        return {key: numpy_output(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return type(value)(numpy_output(item) for item in value)
    if array.namespace(value) is not np:
        return array.to_numpy(value)[()]
    return value


# Assign ordinary values to histogram bins with numpy histogram edge semantics
def bin_indices(values: np.ndarray, bins: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    bins = np.asarray(bins, dtype=float)
    indices = np.searchsorted(bins, values, side="right") - 1
    indices[values == bins[-1]] = len(bins) - 2
    valid = np.isfinite(values) & (indices >= 0) & (indices < len(bins) - 1)
    return np.where(valid, indices, -1)


# Compute a one-dimensional weighted histogram
def hist(x, bins=30, density=False, weights=None):
    x = np.asarray(x)

    weights = np.ones_like(x) if weights is None else array.asarray(weights)

    if len(weights) != len(x):
        raise Exception(f"hist.hist: len(weights) = {len(weights)} != len(x) = {len(x)}")

    bins = np.histogram_bin_edges(x, bins=bins)
    counts, errs = uncertainty.binned(x, bins, weights, bin_indices(x, bins))
    cbins = edge2centerbins(bins)

    # Density integral 1 over the histogram bins range
    if density:
        counts, errs = hist_to_density(counts=counts, errs=errs, bins=bins)

    return counts, errs, bins, cbins


def hist_obj(x, bins=30, weights=None):
    """A wrapper to return a histogram object."""
    counts, errs, bins, cbins = hist(x, bins=bins, weights=weights)
    return hobj(counts, errs, bins, cbins)


# Assign each event once to a category, including the diagonal under proton exchange
def assign_2d_hyperbins(x: np.ndarray, hyperbins: dict, symmetrized: bool = False) -> dict:
    members = {subspace: [] for subspace in hyperbins}
    for index in range(len(x)):
        x0, x1 = x[index, 0], x[index, 1]
        for subspace in hyperbins:
            x0_bins = hyperbins[subspace][0]
            x1_bins = hyperbins[subspace][1]
            accepted = (x0_bins[0] < x0 <= x0_bins[1]) and (x1_bins[0] < x1 <= x1_bins[1])
            if symmetrized and not accepted:
                accepted = (x0_bins[0] < x1 <= x0_bins[1]) and (x1_bins[0] < x0 <= x1_bins[1])
            if accepted:
                members[subspace].append(index)
                break
    return {subspace: np.asarray(members[subspace], dtype=np.int32) for subspace in hyperbins}


# Format one histogram category interval for an internal record name
def bins2txt(bins):
    low, high = (np.format_float_positional(float(value), unique=True, min_digits=2, trim="k") for value in bins)
    return f"[{low}-{high}]"
