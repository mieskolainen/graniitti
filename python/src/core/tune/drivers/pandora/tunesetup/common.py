# Pandora parameter domains and dataset steering
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import copy
import math
import os
from pathlib import Path

from ray import tune

from core.tune.drivers.pandora import driver as pandora_driver


# Match a steering table to its XML algorithm scope or wrapper parameter
def _entry_matches_table(entry, name, table):
    if (entry.get("tag") or entry.get("name")) not in table or "endcap" in str(entry).lower():
        return False
    if name == "wrapper":
        return entry["kind"] == "wrapper_numeric"
    stack = entry.get("algorithm_stack", [])
    return entry["kind"] == "xml_numeric" and (
        (not stack and name == "pandora") or bool(stack) and (
            stack == name.split("/") if "/" in name else stack[-1] == name))


# Select each declared parameter once, restricting primary clustering to its first instance
def active_tune_entries(catalog, *, tuning_tables):
    entries = catalog["xml_numeric"] + catalog["wrapper_numeric"]
    selected = {}
    for name, table in tuning_tables.items():
        matches = [e for e in entries if _entry_matches_table(e, name, table)]
        if name == "ClusteringParent/ConeClustering" and matches:
            first = min(e.get("optional_parent_path", e["path"])[0] for e in matches)
            matches = [e for e in matches if e.get("optional_parent_path", e["path"])[0] == first]
        missing = set(table) - {e.get("tag") or e.get("name") for e in matches}
        if missing:
            raise ValueError(f"Pandora tuning table {name} has unmatched parameters: {sorted(missing)}")
        for entry in matches:
            if entry["key"] in selected:
                raise ValueError(f"Repeated Pandora parameter: {entry['key']}")
            selected[entry["key"]] = entry
    return list(selected.values())


# Compute a finite, nonempty sampling domain from the explicit steering bounds
def default_bounds_for_entry(entry, *, tuning_tables):
    name, spec = next((name, table[entry.get("tag") or entry.get("name")])
                      for name, table in tuning_tables.items() if _entry_matches_table(entry, name, table))
    lower, upper = map(float, spec["bounds"])
    sampling = spec["sampling"]
    if sampling not in (["between_bounds"], ["fixed"]):
        raise ValueError(f"Invalid Pandora sampling in {name}: {sampling}")
    dtype = spec["tune_dtype"]
    if (dtype not in ("int", "float") or not all(map(math.isfinite, (lower, upper)))
            or (upper != lower if sampling == ["fixed"] else upper <= lower)):
        raise ValueError(f"Invalid Pandora domain in {name}: {spec}")
    if dtype == "int" and math.floor(upper) < math.ceil(lower):
        raise ValueError(f"Empty Pandora integer domain in {name}: {spec}")
    return lower, upper


# Build the XML and wrapper catalog and annotate explicit sampling domains
def build_catalog(settings_xml, *, optional_xml_numeric_defaults, wrapper_defaults, tuning_tables, ordered_parameters=()):
    catalog = pandora_driver.build_catalog(settings_xml, optional_xml_numeric_defaults=optional_xml_numeric_defaults,
                                            wrapper_defaults=wrapper_defaults)
    for entry in pandora_driver._entry_by_key(catalog).values():
        entry.update(tune_dtype="int" if entry["kind"] == "xml_boolean" else "float", sampling_mode="inactive")
    entries = active_tune_entries(catalog, tuning_tables=tuning_tables)
    for entry in entries:
        for name, table in tuning_tables.items():
            if _entry_matches_table(entry, name, table):
                spec = table[entry.get("tag") or entry.get("name")]
                entry.update(tune_dtype=spec["tune_dtype"], sampling_mode=spec["sampling"][0],
                             tune_bounds=default_bounds_for_entry(entry, tuning_tables=tuning_tables))
                if "initial" in spec:
                    initial = float(spec["initial"])
                    low, high = entry["tune_bounds"]
                    if not low <= initial <= high or (spec["tune_dtype"] == "int" and not initial.is_integer()):
                        raise ValueError(f"Invalid Pandora initial value in {name}: {spec}")
                    entry["initial"] = initial
    catalog["ordered_parameters"] = []
    for row in ordered_parameters:
        selected = []
        if len(row) < 2:
            raise ValueError("Pandora parameter ordering requires at least two parameters")
        for reference in row:
            scope, _, tag = reference.rpartition("/")
            matches = [e for e in entries if _entry_matches_table(e, scope, {tag: None})]
            if len(matches) != 1:
                raise ValueError(f"Pandora parameter ordering has an ambiguous or missing parameter: {reference}")
            selected.append(matches[0])
        for left, right in zip(selected[:-1], selected[1:], strict=True):
            if left["tune_bounds"][1] > right["tune_bounds"][0]:
                raise ValueError(f"Pandora ordered bounds overlap: {left['tag']} <= {right['tag']}")
        catalog["ordered_parameters"].append([e["key"] for e in selected])
    return catalog


# Resolve inclusive input file ranges without repeated paths or truncated indices
def dataset_input_roots(param_paths, dataset):
    input_dir = Path(os.environ.get("PANDORA_INPUT_DIR", dataset["input_dir"])).expanduser()
    if not input_dir.is_absolute():
        input_dir = Path(param_paths["pandora_dir"]).expanduser() / input_dir
    indices = []
    for first, last in dataset["file_ranges"]:
        if type(first) is not int or type(last) is not int or first < 0 or last < first:
            raise ValueError(f"Invalid Pandora input range: {first}, {last}")
        indices.extend(range(first, last + 1))
    paths = [str(input_dir / dataset["file_template"].format(index=i)) for i in indices]
    if not paths or len(paths) != len(set(paths)):
        raise ValueError("Pandora input ranges must produce distinct, nonempty paths")
    return paths


# Build datacards from the explicit datasets and one common particle flow objective
def datacards_from_selection(*, param_paths, param_selection, param_loss, baseline_config=None):
    root = Path(param_paths["pandora_dir"]).expanduser().resolve()
    generator = param_selection["generator"]
    active = generator["active_datasets"]
    if not active or len(active) != len(set(active)):
        raise ValueError("Pandora active datasets must be distinct and nonempty")
    objective = pandora_driver.pflow_objective_spec(param_selection=param_selection, param_loss=param_loss)
    cards = []
    for name in active:
        dataset = generator["datasets"][name]
        card = dict(name=name.upper(), dataset_name=name,
            process=dataset["process"], sqrts_gev=float(dataset["sqrts_gev"]), weight=float(dataset["weight"]),
            pandora_dir=str(root), run_dir=str(root / "run"), input_roots=dataset_input_roots(param_paths, dataset),
            settings_template=str(root / param_paths["pandora_default_xml"]),
            reco_steering=str(root / param_paths["pandora_reco_py"]), nevents=dataset["nevents"],
            baseline_config=baseline_config, objectives={"pflow": copy.deepcopy(objective)})
        pandora_driver._dataset_weight(card)
        if type(card["nevents"]) is not int or card["nevents"] == 0 or card["nevents"] < -1:
            raise ValueError("Pandora nevents must be -1 or a positive integer")
        if not math.isfinite(card["sqrts_gev"]) or card["sqrts_gev"] <= 0.0:
            raise ValueError("Pandora sqrts_gev must be finite and positive")
        cards.append(card)
    if not any(card["weight"] > 0.0 for card in cards):
        raise ValueError("Pandora requires a positive weight dataset")
    return cards


# Build datacards and parameter spaces for explicit Pandora steering
def setup(*, param_paths, param_selection, param_loss, optional_xml_numeric_defaults,
                     wrapper_defaults, tuning_tables, ordered_parameters=(), baseline_config=None):
    cards = datacards_from_selection(param_paths=param_paths, param_selection=param_selection,
                                     param_loss=param_loss, baseline_config=baseline_config)
    catalog = build_catalog(cards[0]["settings_template"], optional_xml_numeric_defaults=optional_xml_numeric_defaults,
                            wrapper_defaults=wrapper_defaults, tuning_tables=tuning_tables, ordered_parameters=ordered_parameters)
    entries = active_tune_entries(catalog, tuning_tables=tuning_tables)
    param_space = {}
    for entry in entries:
        if entry["sampling_mode"] == "fixed":
            continue
        lower, upper = entry["tune_bounds"]
        param_space[entry["key"]] = (tune.randint(math.ceil(lower), math.floor(upper) + 1)
                                    if entry["tune_dtype"] == "int" else tune.uniform(lower, upper))
    aux = dict(pandora_catalog=catalog, pandora_paths=copy.deepcopy(param_paths))
    return cards, param_space, aux
