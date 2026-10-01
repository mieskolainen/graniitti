# Statistical acceptance of published measurement comparisons
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import logging
import math

from scipy.stats import chi2


# Compute the upper tail probability for a finite chi-square with positive rank
def chi2_probability(value, ndf):
    if value is None or not math.isfinite(value) or value < 0 or ndf is None or not math.isfinite(ndf) or ndf <= 0:
        return None
    return float(chi2.sf(value, ndf))


# Yield every histogram and point prediction in one comparison report
def comparison_samples(report):
    for dataset in report["sets"]:
        if not dataset.get("mc", True):
            continue
        for observable in dataset["observables"]:
            if observable["kind"] in {"histogram", "point"}:
                for sample in observable["samples"]:
                    yield dataset["name"], observable["observable"], observable["kind"], sample


# Test every measured spectrum and rate with a Bonferroni bound per report
def measurement_assessment(report):
    criteria = report.get("validation", {}).get("measurement")
    if criteria is None:
        return {"status": "not_tested", "failures": [], "tests": 0}
    rows = [row for row in comparison_samples(report) if row[3].get("comparison_status") != "prediction"]
    absolute = report["normalization"] == "cross_section"
    tests = sum(1 + int(absolute and kind == "histogram" and observable != "cross_section")
                for _, observable, kind, _ in rows)
    threshold = criteria["alpha"] / max(tests, 1)
    failures = [] if rows else ["No measurement comparisons"]
    if report.get("reference_normalization", report["normalization"]) != report["normalization"]:
        failures.append("Plot normalization differs from the measurement criterion")
    for dataset in report["sets"]:
        if not dataset.get("mc", True):
            continue
        if not dataset["observables"] or any(not item["samples"] for item in dataset["observables"]):
            failures.append(f"{dataset['name']}: missing prediction")
    if not absolute and report.get("density_uncertainty") != "shape":
        failures.append("Shape comparisons require the normalization covariance")
    invalid = failures.copy()
    for dataset, observable, kind, sample in rows:
        context = f"{dataset}/{observable}/{sample['label']}"
        errors = []
        if sample.get("comparison_status") != "compared" or sample.get("empty", True):
            errors.append("measurement comparison is unavailable or empty")
        probabilities = {"chi2": chi2_probability(sample.get("chi2"), sample.get("ndf", 0))}
        if absolute and kind == "histogram" and observable != "cross_section":
            pull = sample.get("integral_pull")
            probabilities["integral"] = None if pull is None or not math.isfinite(pull) else math.erfc(abs(pull) / math.sqrt(2))
        for name, value in probabilities.items():
            if value is None or value < threshold:
                errors.append(f"{name} p={value} below {threshold:.4g} or unavailable")
        for field, limit in (("mc_data_error_ratio", criteria["max_mc_error_ratio"]),
                             ("mc_max_rel_uncertainty", criteria["max_mc_rel_uncertainty"])):
            value = sample.get(field)
            if value is None or not math.isfinite(value) or not 0 <= value <= limit:
                errors.append(f"{field}={value} exceeds {limit:.4g} or unavailable")
        if sample.get("data_valid_bin_count") is None or sample.get("mc_valid_bin_count") != sample["data_valid_bin_count"]:
            errors.append("MC does not cover all measured bins")
        if sample.get("mc_zero_bin_count") is None or sample["mc_zero_bin_count"] > 0:
            errors.append("MC contains unpopulated measurement bins")
        sample["measurement"] = {"status": "failed" if invalid or errors else "passed",
                                 "pvalues": probabilities, "failures": invalid + errors}
        failures.extend(f"{context}: {error}" for error in errors)
    return {"status": "failed" if failures else "passed", "tests": tests,
            "min_pvalue": threshold, "uncertainty_model": "symmetric_gaussian", "failures": failures}


# Record failed measurement criteria without interrupting the comparison
def check_measurements(report):
    report["measurement"] = measurement_assessment(report)
    for failure in report["measurement"]["failures"]:
        logging.getLogger(__name__).warning("%s", failure)


# Compute failures from optional Monte Carlo precision and population bounds
def mc_precision_failures(validation, sample, context):
    failures = []
    max_rel_uncertainty = validation.get("max_mc_rel_uncertainty")
    if max_rel_uncertainty is not None:
        zero_bins = sample.get("mc_zero_bin_count")
        relative_uncertainty = sample.get("mc_max_rel_uncertainty")
        if zero_bins is None:
            failures.append(f"{context}: MC zero-bin count is unavailable")
        elif zero_bins > 0:
            failures.append(
                f"{context}: {zero_bins} valid MC bin(s) have zero content and undefined "
                "relative uncertainty"
            )
        elif relative_uncertainty is None:
            failures.append(f"{context}: relative MC uncertainty is unavailable")
        elif (
            not math.isfinite(relative_uncertainty)
            or relative_uncertainty < 0.0
            or relative_uncertainty > max_rel_uncertainty
        ):
            failures.append(
                f"{context}: mc_max_rel_uncertainty={relative_uncertainty:.6g} "
                f"exceeds {max_rel_uncertainty:.6g}"
            )

    population_bounds = (
        ("mc_effective_events", "min_mc_effective_events"),
        ("mc_ess_fraction", "min_mc_ess_fraction"),
    )
    for field, limit_name in population_bounds:
        limit = validation.get(limit_name)
        if limit is None:
            continue
        value = sample.get(field)
        if value is None:
            failures.append(f"{context}: {field} is unavailable")
        elif not math.isfinite(value) or value < limit:
            failures.append(f"{context}: {field}={value:.6g} is below {limit:.6g}")
    return failures


# Apply optional dataset-specific physics bounds to one iceplot report
def check_physics(report):
    check_measurements(report)
    validation = report.get("validation", {})
    allow_empty = validation.get("allow_empty", False)
    comparisons = 0
    differential_comparisons = 0
    fiducial_integral_comparisons = 0
    predictions = 0
    failures = []
    comparison_required = any(validation.get(key) for key in (
        "mc_reference", "require_comparison", "require_differential_comparison",
        "require_fiducial_integral_comparison",
    )) or any(validation.get(key) is not None for key in (
        "max_chi2_ndf", "max_shape_l1", "max_abs_integral_pull", "max_integral_factor",
    ))
    for dataset, observable, kind, sample in comparison_samples(report):
        predictions += 1
        context = f"{dataset}/{observable}/{sample['label']}"
        failures.extend(mc_precision_failures(validation, sample, context))
        if sample["empty"] and not allow_empty:
            failures.append(f"{context}: empty MC histogram")
            continue
        if sample.get("comparison_status") == "prediction":
            continue
        if sample.get("comparison_status") == "reference":
            if validation.get("mc_reference") is None:
                failures.append(f"{context}: unexpected MC reference without mc_reference steering")
            continue
        if "comparison_status" not in sample:
            if comparison_required:
                failures.append(f"{context}: MC comparison status is unavailable")
            continue
        compared = sample["comparison_status"] == "compared"
        if comparison_required and not compared:
            failures.append(f"{context}: MC comparison is unavailable")
        comparisons += int(compared)
        if compared and observable != "cross_section" and sample.get("mc_valid_bin_count", 0) > 1 and all(
            value is not None and math.isfinite(value)
            for value in (sample.get("chi2_ndf"), sample.get("shape_l1"))
        ):
            differential_comparisons += 1
        if compared and kind == "histogram" and all(
            value is not None and math.isfinite(value)
            for value in (
                sample.get("integral"),
                sample.get("data_integral"),
                sample.get("integral_pull"),
            )
        ):
            fiducial_integral_comparisons += 1
        bounds = (
            ("chi2_ndf", "max_chi2_ndf", False),
            ("shape_l1", "max_shape_l1", False),
            ("integral_pull", "max_abs_integral_pull", True),
        )
        for field, limit_name, absolute in bounds:
            # The single fiducial rate uses its integral pull, not a differential chi-square bound
            if field == "chi2_ndf" and observable == "cross_section":
                continue
            limit = validation.get(limit_name)
            value = sample.get(field)
            if limit is None:
                continue
            selected = abs(value) if absolute and value is not None else value
            if selected is not None and math.isfinite(selected) and 0.0 <= selected <= limit:
                continue
            failure = (f"{context}: {field} is unavailable" if value is None else
                       f"{context}: {field}={value:.6g} exceeds {limit:.6g}")
            if field == "chi2_ndf":
                logging.getLogger(__name__).warning("%s", failure)
            else:
                failures.append(failure)
        max_factor = validation.get("max_integral_factor")
        data_integral = sample.get("data_integral")
        mc_integral = sample.get("integral")
        if max_factor is not None:
            if data_integral is None or mc_integral is None:
                failures.append(f"{context}: integral factor requires both integrals")
            elif not all(math.isfinite(value) and value > 0.0 for value in (data_integral, mc_integral)):
                failures.append(f"{context}: integral factor requires finite positive integrals")
            else:
                factor = max(mc_integral / data_integral, data_integral / mc_integral)
                if not math.isfinite(factor) or factor > max_factor:
                    failures.append(
                        f"{context}: integral factor={factor:.6g} exceeds {max_factor:.6g}"
                    )
    if predictions == 0:
        failures.append("iceplot report contains no predictions")
    if validation.get("require_comparison", False) and comparisons == 0:
        failures.append("iceplot report contains no MC-to-reference comparison")
    if validation.get("require_differential_comparison", False) and differential_comparisons == 0:
        failures.append("iceplot report contains no finite differential comparison")
    if (
        validation.get("require_fiducial_integral_comparison", False)
        and fiducial_integral_comparisons == 0
    ):
        failures.append("iceplot report contains no finite fiducial integral comparison")
    if failures:
        raise ValueError("\n".join(failures))
