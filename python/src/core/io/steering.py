# Shared resolution and loading for icepack steering references
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import argparse
import copy
import json
import math
import os
import pathlib
import re
import sys
from functools import partial
from importlib import import_module
from inspect import signature
from types import ModuleType

import pyjson5 as json5

from core import resource
from core.io.serialize import load_json_file


# Compute the normalized GRANIITTI project root
def project_root(cdir: str | os.PathLike[str]) -> pathlib.Path:
    return pathlib.Path(cdir).expanduser().resolve()


# Compute true when a reference explicitly uses dataset-relative syntax
def is_dataset_relative(reference: str) -> bool:
    return reference.startswith("./") or reference.startswith("../")


# Resolve absolute, project, dataset and package file references
def resolve_file_reference(
    reference: str,
    *,
    cdir: str | os.PathLike[str],
    dataset_path: str | os.PathLike[str] | None = None,
    suffix: str | None = None,
) -> str:
    value = os.fspath(reference).strip()
    if not value:
        raise ValueError("Empty steering file reference")

    path = resource(value.removeprefix("package:")) if value.startswith("package:") else pathlib.Path(value).expanduser()
    if path.is_absolute():
        candidate = path
    elif is_dataset_relative(value):
        if dataset_path is None:
            raise ValueError(f"Dataset-relative reference '{value}' requires a dataset card path")
        candidate = pathlib.Path(dataset_path).expanduser().resolve().parent / path
    else:
        candidate = project_root(cdir) / path

    candidate = candidate.resolve()
    if suffix is not None and candidate.suffix != suffix:
        raise ValueError(f"Steering reference '{candidate}' must use suffix '{suffix}'")
    if not candidate.is_file():
        raise FileNotFoundError(f"Steering file '{candidate}' not found")
    return str(candidate)


# Resolve a dataset card path from an explicit icepack reference
def resolve_dataset_reference(reference: str, *, cdir: str | os.PathLike[str]) -> str:
    return resolve_file_reference(reference, cdir=cdir, suffix=".json")


# Resolve a generator card relative to its dataset card or project root
def resolve_gencard_reference(
    reference: str, *, dataset_path: str | os.PathLike[str], cdir: str | os.PathLike[str]
) -> str:
    return resolve_file_reference(reference, cdir=cdir, dataset_path=dataset_path, suffix=".json")


# Resolve a data file from the shared datapath or explicitly from the dataset bundle
def resolve_data_reference(
    reference: str,
    *,
    datapath: str | os.PathLike[str],
    dataset_path: str | os.PathLike[str],
    cdir: str | os.PathLike[str],
) -> str:
    value = os.fspath(reference).strip()
    if not value:
        raise ValueError("Empty data file reference")
    if pathlib.Path(value).is_absolute() or is_dataset_relative(value):
        return resolve_file_reference(value, cdir=cdir, dataset_path=dataset_path)

    base = pathlib.Path(datapath).expanduser()
    candidate = base / value if base.is_absolute() else project_root(cdir) / base / value
    return resolve_file_reference(str(candidate), cdir=cdir)


# Resolve a built-in module name or an explicit Python definition file
def resolve_python_reference(
    reference: str, *, package: str, dataset_path: str | os.PathLike[str] | None, cdir: str | os.PathLike[str]
) -> str:
    value = os.fspath(reference).strip()
    if not value:
        raise ValueError(f"Empty {package} reference")
    if value.endswith(".py") or pathlib.Path(value).is_absolute() or is_dataset_relative(value):
        return resolve_file_reference(value, cdir=cdir, dataset_path=dataset_path, suffix=".py")
    if value.startswith(f"{package}."):
        return value
    if "/" in value or "\\" in value:
        raise ValueError(f"File reference '{value}' must end in '.py'; built-in references use '{package}.name'")
    return f"{package}.{value}"


# Load one Python module from an import name or resolved source path
def load_python_module(reference: str) -> ModuleType:
    if not reference.endswith(".py"):
        return import_module(reference)
    path = pathlib.Path(reference).resolve()
    root = next((parent for parent in path.parents if parent.name == "icepack"), path.parent)
    if str(root.parent) not in sys.path:
        sys.path.insert(0, str(root.parent))
    return import_module(".".join(path.relative_to(root.parent).with_suffix("").parts))


# Extract and validate observable dictionaries from one loaded module
def observables_from_module(module: ModuleType, *, reference: str) -> dict:
    observables = {}
    for name, value in vars(module).items():
        if not name.startswith("obs_"):
            continue
        if not isinstance(value, dict) or "tag" not in value:
            raise ValueError(f"Observable definition '{reference}:{name}' must contain 'tag'")
        kind = value.get("kind", "histogram")
        if kind not in {"histogram", "point", "roc"}:
            raise ValueError(f"Observable definition '{reference}:{name}' has unknown kind {kind!r}")
        if kind != "point" and not callable(value.get("func")):
            raise ValueError(f"Observable definition '{reference}:{name}' must contain callable 'func'")
        tag = str(value["tag"])
        if tag in observables:
            raise ValueError(f"Observable tag '{tag}' is duplicated in '{reference}'")
        observables[tag] = value
    if not observables:
        raise ValueError(f"Observable definition '{reference}' contains no obs_* dictionaries")
    return observables


# Select the observable definitions declared by one dataset set
def select_histogram_observables(dataset_set: dict, all_obs: dict) -> dict:
    selected = {}
    for histogram in dataset_set["hist"]:
        observable = histogram["obs"]
        if observable not in all_obs:
            raise KeyError(f'Observable "{observable}" not found in set observables')
        selected[observable] = copy.deepcopy(all_obs[observable])
        if "args" in histogram:
            function = selected[observable]["func"]
            signature(function).bind(None, **histogram["args"])
            selected[observable]["func"] = partial(function, **histogram["args"])
        selected[observable]["differential"] = histogram.get("differential", True)
        if not selected[observable]["differential"]:
            selected[observable]["units"]["yden"] = "1"
    return selected


# Resolve and load one set-level observable definition
def load_observables(
    reference: str, *, dataset_path: str | os.PathLike[str] | None, cdir: str | os.PathLike[str]
) -> tuple[dict, str]:
    resolved = resolve_python_reference(reference, package="core.analysis.observables", dataset_path=dataset_path, cdir=cdir)
    module = load_python_module(resolved)
    return observables_from_module(module, reference=resolved), resolved


# Validate one generator-sample name used in output filenames
def validate_sample_name(value: object, *, context: str) -> str:
    if not isinstance(value, str) or re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_+-]*", value) is None:
        raise ValueError(f"{context} name must match [A-Za-z0-9][A-Za-z0-9_+-]*")
    return value


# Validate one non-empty label without terminal control characters
def validate_sample_label(value: object, *, context: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{context} label must be a non-empty string")
    if any(ord(character) < 32 for character in value):
        raise ValueError(f"{context} label must not contain control characters")
    return value


# Validate one non-empty display string without terminal control characters
def validate_display_text(value: object, *, context: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{context} must be a non-empty string")
    if any(ord(character) < 32 for character in value):
        raise ValueError(f"{context} must not contain control characters")
    return value


# Serialize one finite JSON value for the generator override interface
def serialize_override_value(value: object, *, context: str) -> str:
    try:
        return json.dumps(value, ensure_ascii=True, allow_nan=False, separators=(",", ":"))
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{context} must be a finite JSON value") from exc


# Validate one map of generator JSON path overrides
def validate_sample_parameters(value: object, *, context: str) -> dict:
    if not isinstance(value, dict):
        raise TypeError(f"{context} parameters must be an object")
    for path, parameter in value.items():
        if (
            not isinstance(path, str)
            or not path.strip()
            or "=" in path
            or any(ord(character) < 32 for character in path)
        ):
            raise ValueError(f"{context} contains an invalid generator override path {path!r}")
        serialize_override_value(parameter, context=f"{context} parameters[{path!r}]")
    return value


# Validate one optional Pythia fragmentation stage
def validate_sample_fragmentation(
    value: object,
    *,
    context: str,
    dataset_path: str,
    cdir: str | os.PathLike[str],
) -> dict:
    if not isinstance(value, dict):
        raise TypeError(f"{context} fragmentation must be an object")
    validate_mapping_keys(
        value,
        allowed={"engine", "card", "mode", "seed", "attempts"},
        required={"engine", "card"},
        context=f"{context} fragmentation",
    )
    engine = validate_schema_choice(
        value["engine"], allowed={"pythia"}, context=f"{context} fragmentation.engine"
    )
    card = resolve_file_reference(
        value["card"], cdir=cdir, dataset_path=dataset_path, suffix=".cmnd"
    )
    mode = validate_schema_choice(
        value.get("mode", "auto"),
        allowed={"auto", "fragment", "shower"},
        context=f"{context} fragmentation.mode",
    )
    seed = value.get("seed", 0)
    attempts = value.get("attempts", 50)
    for name, number, lower in (("seed", seed, 0), ("attempts", attempts, 1)):
        if isinstance(number, bool) or not isinstance(number, int) or number < lower:
            raise ValueError(f"{context} fragmentation.{name} must be an integer at least {lower}")
    return {"engine": engine, "card": card, "mode": mode, "seed": seed, "attempts": attempts}


# Validate the allowed and required keys of one dataset mapping
def validate_mapping_keys(value: dict, *, allowed: set[str], required: set[str], context: str) -> None:
    unknown = set(value).difference(allowed)
    if unknown:
        raise KeyError(f"{context} has unknown keys {sorted(unknown)}")
    missing = required.difference(value)
    if missing:
        raise KeyError(f"{context} is missing keys {sorted(missing)}")


# Validate one string against the supported schema values
def validate_schema_choice(value: object, *, allowed: set[str], context: str) -> str:
    if not isinstance(value, str):
        raise TypeError(f"{context} must be a string")
    if value not in allowed:
        choices = ", ".join(repr(choice) for choice in sorted(allowed))
        raise ValueError(f"{context} must be one of {choices}")
    return value


# Compute the generator sample names selected by one dataset set
def set_sample_names(entry: dict) -> tuple[str, ...]:
    if entry.get("mc", True) is False:
        return ()
    if "samples" in entry:
        return tuple(entry["samples"])
    if "sample" in entry:
        return (entry["sample"],)
    return ()


# Validate the plot, fit and generator sample blocks
def validate_dataset_samples(
    dataset: dict, *, path: str, cdir: str | os.PathLike[str] = "."
) -> None:
    if "plot" not in dataset:
        raise KeyError(f"Dataset card '{path}' is missing 'plot'")
    if "fit" not in dataset:
        raise KeyError(f"Dataset card '{path}' is missing 'fit'")

    samples = dataset.get("samples")
    plot = dataset["plot"]
    fit = dataset["fit"]

    if not isinstance(plot, dict):
        raise TypeError(f"Dataset card '{path}' plot must be an object")
    validate_mapping_keys(
        plot,
        allowed={
            "normalization",
            "unit",
            "density_uncertainty",
            "ratio_plot",
            "ratio_uncertainty",
            "stack",
            "data_style",
            "data_linestyle",
            "mc_linestyle",
            "match_colors",
            "mc_first",
        },
        required={"normalization", "ratio_uncertainty", "stack", "data_style"},
        context=f"Dataset card '{path}' plot",
    )

    normalization = validate_schema_choice(
        plot["normalization"],
        allowed={"cross_section", "unit_density"},
        context=f"Dataset card '{path}' plot.normalization",
    )
    if "unit" in plot:
        validate_schema_choice(
            plot["unit"], allowed={"b", "mb", "ub", "nb", "pb", "fb"}, context=f"Dataset card '{path}' plot.unit"
        )
    if normalization == "unit_density" and "density_uncertainty" not in plot:
        raise KeyError(f"Dataset card '{path}' plot.density_uncertainty is required for unit_density")
    if normalization == "cross_section" and "density_uncertainty" in plot:
        raise KeyError(f"Dataset card '{path}' plot.density_uncertainty is only valid for unit_density")
    if "density_uncertainty" in plot:
        validate_schema_choice(
            plot["density_uncertainty"],
            allowed={"shape", "scaled"},
            context=f"Dataset card '{path}' plot.density_uncertainty",
        )

    validate_schema_choice(
        plot["ratio_uncertainty"],
        allowed={"combined", "separate", "numerator", "none"},
        context=f"Dataset card '{path}' plot.ratio_uncertainty",
    )
    if "ratio_plot" in plot and not isinstance(plot["ratio_plot"], bool):
        raise TypeError(f"Dataset card '{path}' plot.ratio_plot must be boolean")
    for key in ("match_colors", "mc_first"):
        if key in plot and not isinstance(plot[key], bool):
            raise TypeError(f"Dataset card '{path}' plot.{key} must be boolean")
    if not isinstance(plot["stack"], bool):
        raise TypeError(f"Dataset card '{path}' plot.stack must be boolean")
    validate_schema_choice(
        plot["data_style"], allowed={"hist", "errorbar"}, context=f"Dataset card '{path}' plot.data_style"
    )
    for key in ("data_linestyle", "mc_linestyle"):
        if key in plot:
            validate_schema_choice(
                plot[key], allowed={"-", "--", "-.", ":"}, context=f"Dataset card '{path}' plot.{key}"
            )

    if not isinstance(fit, dict):
        raise TypeError(f"Dataset card '{path}' fit must be an object")
    validate_mapping_keys(
        fit, allowed={"normalization"}, required={"normalization"}, context=f"Dataset card '{path}' fit"
    )
    validate_schema_choice(
        fit["normalization"],
        allowed={"cross_section", "unit_density"},
        context=f"Dataset card '{path}' fit.normalization",
    )

    if not isinstance(samples, list) or not samples:
        raise ValueError(f"Dataset card '{path}' samples must be a non-empty array")

    names = set()
    allowed_keys = {"name", "label", "gencard", "parameters", "scale", "fragmentation"}
    for index, sample in enumerate(samples):
        context = f"Dataset card '{path}' samples[{index}]"
        if not isinstance(sample, dict):
            raise TypeError(f"{context} must be an object")
        validate_mapping_keys(
            sample, allowed=allowed_keys, required={"name", "label", "gencard", "parameters"}, context=context
        )

        name = validate_sample_name(sample["name"], context=context)
        if name in names:
            raise ValueError(f"Dataset card '{path}' has duplicate sample name '{name}'")
        names.add(name)
        validate_sample_label(sample["label"], context=context)
        validate_sample_parameters(sample["parameters"], context=context)
        if "fragmentation" in sample:
            validate_sample_fragmentation(
                sample["fragmentation"], context=context, dataset_path=path, cdir=cdir
            )
        gencard = sample["gencard"]
        if not isinstance(gencard, str) or not gencard.strip():
            raise ValueError(f"{context} gencard must be a non-empty string")
        scale = sample.get("scale", 1.0)
        if (
            isinstance(scale, bool)
            or not isinstance(scale, (int, float))
            or not math.isfinite(float(scale))
            or float(scale) <= 0.0
        ):
            raise ValueError(f"{context} scale must be finite and positive")

    if "sets" in dataset:
        routed_sets = []
        set_samples = set()
        for index, entry in enumerate(dataset["sets"]):
            context = f"Dataset card '{path}' sets[{index}]"
            if "sample" in entry and "samples" in entry:
                raise ValueError(f"{context} cannot define both sample and samples")
            if "samples" in entry and (not isinstance(entry["samples"], list) or not entry["samples"]):
                raise TypeError(f"{context} samples must be a non-empty array")
            if entry.get("mc", True) is False and ("sample" in entry or "samples" in entry):
                raise ValueError(f"{context} data-only set cannot select MC samples")
            selected = set_sample_names(entry)
            if any(not isinstance(name, str) or not name.strip() for name in selected):
                raise ValueError(f"{context} sample names must be non-empty strings")
            if len(selected) != len(set(selected)):
                raise ValueError(f"{context} contains duplicate sample names")
            if entry.get("mc", True):
                routed_sets.append(bool(selected))
            set_samples.update(selected)
        if any(routed_sets) and not all(routed_sets):
            raise ValueError(f"Dataset card '{path}' routed sets must all define sample or samples")
        unknown_samples = set_samples.difference(names)
        if unknown_samples:
            raise ValueError(f"Dataset card '{path}' sets contain unknown sample values {sorted(unknown_samples)}")
        missing_samples = names.difference(set_samples) if set_samples else set()
        if missing_samples:
            raise ValueError(f"Dataset card '{path}' samples have no routed sets {sorted(missing_samples)}")

    if plot["stack"] and len(samples) < 2:
        raise ValueError(f"Dataset card '{path}' plot.stack requires at least two samples")


# Validate optional report thresholds used by the generic physics runner
def validate_dataset_validation(dataset: dict, *, path: str) -> None:
    validation = dataset.get("validation", {})
    if not isinstance(validation, dict):
        raise TypeError(f"Dataset card '{path}' validation must be an object")
    validate_mapping_keys(
        validation,
        allowed={
            "loopscreen",
            "nevents",
            "allow_empty",
            "require_comparison",
            "require_differential_comparison",
            "require_fiducial_integral_comparison",
            "mc_reference",
            "max_chi2_ndf",
            "max_shape_l1",
            "max_abs_integral_pull",
            "max_integral_factor",
            "max_mc_rel_uncertainty",
            "min_mc_effective_events",
            "min_mc_ess_fraction",
            "measurement",
        },
        required=set(),
        context=f"Dataset card '{path}' validation",
    )
    if "loopscreen" in validation and validation["loopscreen"] not in (0, 1, False, True):
        raise ValueError(f"Dataset card '{path}' validation.loopscreen must be 0 or 1")
    if "measurement" in validation:
        criteria = validation["measurement"]
        fields = {"alpha", "max_mc_error_ratio", "max_mc_rel_uncertainty"}
        if not isinstance(criteria, dict) or set(criteria) != fields:
            raise ValueError(f"Dataset card '{path}' measurement requires {sorted(fields)}")
        for key, value in criteria.items():
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not 0 < value < 1:
                raise ValueError(f"Dataset card '{path}' measurement.{key} must lie between zero and one")
        if not validation.get("require_comparison") or "mc_reference" in validation:
            raise ValueError(f"Dataset card '{path}' measurement requires a published data comparison")
        if not all(entry.get("data", True) for entry in dataset["sets"]):
            raise ValueError(f"Dataset card '{path}' measurement requires data in every set")
    if "nevents" in validation and (
        isinstance(validation["nevents"], bool)
        or not isinstance(validation["nevents"], int)
        or validation["nevents"] <= 0
    ):
        raise ValueError(f"Dataset card '{path}' validation.nevents must be a positive integer")
    comparison_requirements = (
        "require_comparison",
        "require_differential_comparison",
        "require_fiducial_integral_comparison",
    )
    for key in ("allow_empty", *comparison_requirements):
        if key in validation and not isinstance(validation[key], bool):
            raise TypeError(f"Dataset card '{path}' validation.{key} must be boolean")
    if any(validation.get(key, False) for key in comparison_requirements[1:]) and not validation.get(
        "require_comparison", False
    ):
        raise ValueError(
            f"Dataset card '{path}' differential and fiducial integral requirements need validation.require_comparison = true"
        )
    if "mc_reference" in validation:
        reference = validation["mc_reference"]
        sample_names = {sample["name"] for sample in dataset["samples"]}
        if not isinstance(reference, str) or not reference.strip():
            raise TypeError(f"Dataset card '{path}' validation.mc_reference must be a string")
        if reference not in sample_names:
            raise ValueError(f"Dataset card '{path}' validation.mc_reference names unknown sample '{reference}'")
        if len(dataset["samples"]) < 2:
            raise ValueError(f"Dataset card '{path}' validation.mc_reference requires at least two samples")
        if dataset["plot"]["stack"]:
            raise ValueError(f"Dataset card '{path}' validation.mc_reference does not support stacking")
        if any(entry.get("data", True) for entry in dataset["sets"]):
            raise ValueError(f"Dataset card '{path}' validation.mc_reference requires MC-only sets")
        if any(set_sample_names(entry) for entry in dataset["sets"]):
            raise ValueError(f"Dataset card '{path}' validation.mc_reference does not support routed sets")
    for key in (
        "max_chi2_ndf",
        "max_shape_l1",
        "max_abs_integral_pull",
        "max_integral_factor",
        "max_mc_rel_uncertainty",
        "min_mc_effective_events",
        "min_mc_ess_fraction",
    ):
        if key not in validation or validation[key] is None:
            continue
        value = validation[key]
        if (
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(float(value))
            or float(value) < 0.0
        ):
            raise ValueError(f"Dataset card '{path}' validation.{key} must be non-negative")
    max_integral_factor = validation.get("max_integral_factor")
    if max_integral_factor is not None and float(max_integral_factor) < 1.0:
        raise ValueError(f"Dataset card '{path}' validation.max_integral_factor must be at least one")
    min_ess_fraction = validation.get("min_mc_ess_fraction")
    if min_ess_fraction is not None and float(min_ess_fraction) > 1.0:
        raise ValueError(f"Dataset card '{path}' validation.min_mc_ess_fraction must not exceed one")
    source_population_bounds = ("min_mc_effective_events", "min_mc_ess_fraction")
    if dataset["plot"]["stack"] and any(validation.get(key) is not None for key in source_population_bounds):
        raise ValueError(f"Dataset card '{path}' source population bounds do not support process stacking")


# Build resolved generator commands and plot metadata for one dataset
def build_generation_plan(
    dataset: dict, *, dataset_path: str, cdir: str | os.PathLike[str], output_prefix: str
) -> dict:
    validate_sample_name(output_prefix, context="Generation output prefix")
    configured_samples = dataset["samples"]
    routed = any(set_sample_names(entry) for entry in dataset.get("sets", []))

    samples = []
    for sample in configured_samples:
        name = sample["name"]
        output = output_prefix if len(configured_samples) == 1 else f"{output_prefix}_{name}"
        gencard = resolve_gencard_reference(sample["gencard"], dataset_path=dataset_path, cdir=cdir)
        overrides = [
            f"{path}={serialize_override_value(value, context=f'Generator override {path!r}')}"
            for path, value in sample["parameters"].items()
        ]
        samples.append(
            {
                "output": output,
                "label": sample["label"],
                "gencard": gencard,
                "scale": float(sample.get("scale", 1.0)),
                "assignment": sample["name"] if routed else None,
                "fragmentation": (
                    validate_sample_fragmentation(
                        sample["fragmentation"],
                        context=f"Dataset card '{dataset_path}' sample '{name}'",
                        dataset_path=dataset_path,
                        cdir=cdir,
                    )
                    if "fragmentation" in sample
                    else None
                ),
                "overrides": overrides,
            }
        )

    plot = dataset["plot"]
    return {
        "samples": samples,
        "stack": plot["stack"],
        "density": plot["normalization"] == "unit_density",
        "loopscreen": dataset.get("validation", {}).get("loopscreen"),
        "nevents": dataset.get("validation", {}).get("nevents"),
    }


# Load one JSON5 dataset card and validate its set-level steering keys
def load_dataset(reference: str, *, cdir: str | os.PathLike[str]) -> tuple[dict, str]:
    path = resolve_dataset_reference(reference, cdir=cdir)
    dataset = load_json_file(path, loader=json5.load)

    if not isinstance(dataset, dict):
        raise TypeError(f"Dataset card '{path}' must contain an object")
    validate_mapping_keys(
        dataset,
        allowed={
            "active",
            "type",
            "reader",
            "datapath",
            "samples",
            "reference_sample",
            "sqrts",
            "sets",
            "plot",
            "fit",
            "validation",
        },
        required={"active", "type", "samples", "sets", "plot", "fit"},
        context=f"Dataset card '{path}'",
    )
    if not isinstance(dataset["sets"], list) or not dataset["sets"]:
        raise ValueError(f"Dataset card '{path}' must contain at least one sets entry")
    has_data = any(entry.get("data", True) for entry in dataset["sets"] if isinstance(entry, dict))
    if has_data and "datapath" not in dataset:
        raise KeyError(f"Dataset card '{path}' with data is missing 'datapath'")
    if has_data and dataset["type"] != "RAW_SCALAR" and "reader" not in dataset:
        raise KeyError(f"Dataset card '{path}' with data is missing 'reader'")
    plot_outputs = {}
    plot_group_observables = {}
    for index, entry in enumerate(dataset["sets"]):
        if not isinstance(entry, dict):
            raise TypeError(f"Dataset card '{path}' sets[{index}] must be an object")
        validate_mapping_keys(
            entry,
            allowed={
                "name",
                "plotname",
                "plot_group",
                "title",
                "mc_scale",
                "pid",
                "cuts",
                "obs",
                "hist",
                "state",
                "sample",
                "samples",
                "region",
                "data",
                "mc",
            },
            required={"name", "pid", "cuts", "obs", "hist"},
            context=f"Dataset card '{path}' sets[{index}]",
        )
        context = f"Dataset card '{path}' sets[{index}]"
        validate_display_text(entry["name"], context=f"{context} name")
        if "plotname" in entry:
            validate_display_text(entry["plotname"], context=f"{context} plotname")
        if "plot_group" in entry:
            validate_display_text(entry["plot_group"], context=f"{context} plot_group")
        if "title" in entry:
            validate_display_text(entry["title"], context=f"{context} title")
        if "mc_scale" in entry:
            value = entry["mc_scale"]
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{context} mc_scale must be finite and positive")
            entry["mc_scale"] = float(value)
        if "data" in entry and not isinstance(entry["data"], bool):
            raise TypeError(f"{context} data must be boolean")
        if "mc" in entry and not isinstance(entry["mc"], bool):
            raise TypeError(f"{context} mc must be boolean")
        if not entry.get("data", True) and not entry.get("mc", True):
            raise ValueError(f"{context} must enable data or MC")

        plot_output = entry.get("plot_group", entry.get("plotname", entry["name"])).replace(" ", "_")
        output_kind = "group" if "plot_group" in entry else "set"
        if plot_output in plot_outputs and (output_kind != "group" or plot_outputs[plot_output] != "group"):
            raise ValueError(f"Dataset card '{path}' has duplicate plot output '{plot_output}'")
        plot_outputs[plot_output] = output_kind
        if not isinstance(entry["hist"], list) or not entry["hist"]:
            raise ValueError(f"Dataset card '{path}' sets[{index}] must contain at least one hist entry")
        for hist_index, histogram in enumerate(entry["hist"]):
            if not isinstance(histogram, dict):
                raise TypeError(f"Dataset card '{path}' sets[{index}] hist[{hist_index}] must be an object")
            context = f"Dataset card '{path}' sets[{index}] hist[{hist_index}]"
            validate_mapping_keys(
                histogram,
                allowed={
                    "file",
                    "obs",
                    "args",
                    "scale",
                    "unit",
                    "mc_scale",
                    "differential",
                    "fitw",
                    "rebin_factor",
                    "file_filter",
                    "rows",
                    "point_axis",
                    "covariance",
                    "flux",
                    "xmin",
                    "xmax",
                    "nbins",
                },
                required={"file", "obs", "scale"} if entry.get("data", True) else {"obs"},
                context=context,
            )
            if "args" in histogram and not isinstance(histogram["args"], dict):
                raise TypeError(f"{context} args must be an object")
            if "mc_scale" in histogram:
                values = histogram["mc_scale"]
                values = values if isinstance(values, list) else [values]
                if not values or any(
                    isinstance(value, bool) or not isinstance(value, (int, float))
                    or not math.isfinite(float(value)) or float(value) <= 0.0 for value in values
                ):
                    raise ValueError(f"{context} mc_scale must contain finite positive factors")
            if "differential" in histogram and not isinstance(histogram["differential"], bool):
                raise TypeError(f"{context} differential must be boolean")
            if "rows" in histogram:
                rows = histogram["rows"]
                if (
                    not isinstance(rows, list)
                    or not rows
                    or any(isinstance(index, bool) or not isinstance(index, int) or index < 0 for index in rows)
                    or len(rows) != len(set(rows))
                ):
                    raise ValueError(f"{context} rows must contain unique nonnegative integers")
            for field in ("covariance", "flux"):
                if field in histogram and (not isinstance(histogram[field], str) or not histogram[field].strip()):
                    raise TypeError(f"{context} {field} must name a HEPData table")
            if "point_axis" in histogram and not isinstance(histogram["point_axis"], bool):
                raise TypeError(f"{context} point_axis must be boolean")
            if "unit" in histogram:
                validate_schema_choice(
                    histogram["unit"], allowed={"b", "mb", "ub", "nb", "pb", "fb"}, context=f"{context} unit"
                )
            raw_keys = {"xmin", "xmax", "nbins"}
            if entry.get("data", True) and dataset["type"] == "RAW_SCALAR":
                missing_raw_keys = raw_keys.difference(histogram)
                if missing_raw_keys:
                    raise KeyError(f"{context} is missing RAW_SCALAR keys {sorted(missing_raw_keys)}")
            elif entry.get("data", True) and raw_keys.intersection(histogram):
                raise KeyError(f"{context} RAW_SCALAR keys require type 'RAW_SCALAR'")
        if "plot_group" in entry:
            observables = tuple(item["obs"] for item in entry["hist"])
            previous = plot_group_observables.setdefault(entry["plot_group"], observables)
            if previous != observables:
                raise ValueError(
                    f"Dataset card '{path}' plot_group {entry['plot_group']!r} must use the same observables in every set"
                )
    validate_dataset_samples(dataset, path=path, cdir=cdir)
    validate_dataset_validation(dataset, path=path)
    return dataset, path


# Parse the sample and dataset commands used by the icepack runner
def parse_generation_plan_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Resolve icepack generator samples")
    parser.add_argument("command", choices=["generation-plan", "list-datasets"])
    parser.add_argument("--dataset", required=True, nargs="+")
    parser.add_argument("--cdir", required=True)
    parser.add_argument("--output-prefix")
    args = parser.parse_args()
    if args.command == "generation-plan" and (len(args.dataset) != 1 or args.output_prefix is None):
        parser.error("generation-plan requires one dataset and --output-prefix")
    return args


# Write one generation plan as unambiguous NUL-delimited shell fields
def write_generation_plan(plan: dict) -> None:
    fields = [str(len(plan["samples"])), "1" if plan["stack"] else "0", "1" if plan["density"] else "0"]
    fields.append("" if plan["loopscreen"] is None else str(int(plan["loopscreen"])))
    fields.append("" if plan["nevents"] is None else str(plan["nevents"]))
    for sample in plan["samples"]:
        fragmentation = sample.get("fragmentation")
        fields.extend(
            [
                sample["output"],
                sample["label"],
                sample["gencard"],
                repr(sample["scale"]),
                "" if sample["assignment"] is None else sample["assignment"],
                "" if fragmentation is None else fragmentation["engine"],
                "" if fragmentation is None else fragmentation["card"],
                "" if fragmentation is None else fragmentation["mode"],
                "" if fragmentation is None else str(fragmentation["seed"]),
                "" if fragmentation is None else str(fragmentation["attempts"]),
                str(len(sample["overrides"])),
                *sample["overrides"],
            ]
        )
    sys.stdout.buffer.write(b"".join(field.encode("utf-8") + b"\0" for field in fields))


# Run the icepack sample and dataset commands
def main() -> None:
    args = parse_generation_plan_args()
    if args.command == "list-datasets":
        for reference in args.dataset:
            path = resolve_dataset_reference(reference, cdir=args.cdir)
            dataset = load_json_file(path, loader=json5.load)
            if dataset["active"]:
                print(pathlib.Path(path).parent)
        return
    dataset, dataset_path = load_dataset(args.dataset[0], cdir=args.cdir)
    plan = build_generation_plan(dataset, dataset_path=dataset_path, cdir=args.cdir, output_prefix=args.output_prefix)
    write_generation_plan(plan)


if __name__ == "__main__":
    main()
