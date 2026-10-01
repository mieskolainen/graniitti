# PandoraPFA specific icetune driver and objective helpers
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import ast
import copy
import json
import math
import os
import pathlib
import re
import shlex
import shutil
import sys
import time
import xml.etree.ElementTree as ET
import xml.parsers.expat
from typing import Any

import numpy as np

from core.io import logger as logger_tools
from core.io.files import ensure_dir
from core.io.serialize import finite_or_none as _safe_float
from core.io.serialize import json_default, write_json_file
from core.tune import core as icetune
from core.tune import push as icetune_push
from core.tune import short_id
from core.tune.drivers.base import SimulatorDriver
from core.tune.drivers.pandora import diagnostics, plots, runtime
from core.tune.runtime.process import command_output_root_cause, command_output_tail, run_process

logger = logger_tools.get_logger(__name__)

PENALTY_COST = 1.0e12
CATALOG_SCHEMA_VERSION = 1

_NUMERIC_RE = re.compile(r"^[-+]?((\d+(\.\d*)?)|(\.\d+))([eE][-+]?\d+)?$")
_PARAM_SAFE_RE = re.compile(r"[^A-Za-z0-9_]+")

# Reject unstable particles in post decay visible truth
PFLOW_NONOBSERVABLE_TRUTH_ABS_PDGS = (15, 111, 311, 3112, 3212, 3222, 3312, 3322, 3334)


# Compute a path list while preserving the scalar path call convention
def _as_path_list(paths: Any) -> list[pathlib.Path]:
    return list(map(pathlib.Path, [paths] if isinstance(paths, (str, os.PathLike)) else paths))


# Compute a stable parameter key from its XML position
def _param_key(prefix, index, tag):
    return f"{prefix}_{index:04d}_{_PARAM_SAFE_RE.sub('_', str(tag)).strip('_') or 'value'}"


# Find an XML element by its catalog child indices
def _element_children_path(root, path):
    for index in path:
        root = root[int(index)]
    return root


# Identify an XML target by every ancestor tag and attribute
def _xml_identity(root, path):
    identity = [[root.tag, dict(root.attrib)]]
    for index in path:
        root = root[index]
        identity.append([root.tag, dict(root.attrib)])
    return identity


# Visit XML nodes once with their child paths and enclosing algorithm scopes
def _walk_xml_nodes(root, path=(), text=None, stack=()):
    text = root.tag if text is None else text
    if root.tag == "algorithm" and root.get("type"):
        stack = (*stack, root.get("type"))
    elif len(path) == 1 and (len(root) or not (root.text or "").strip()):
        stack = (root.tag,)
    yield root, path, text, stack
    for index, child in enumerate(root):
        yield from _walk_xml_nodes(child, (*path, index), f"{text}/{child.tag}[{index}]", stack)


# Compute numeric, boolean and optional entries with stable keys and shared XML metadata
def catalog_pandora_settings(settings_xml, *, optional_xml_numeric_defaults):
    root = ET.parse(settings_xml).getroot()
    nodes = list(_walk_xml_nodes(root))
    groups = {"xml_numeric": [], "xml_boolean": []}

    # Record the common XML location and its scalar default
    def add(node, tag, text, section, **extra):
        _, path, path_text, stack = node
        entries = groups[section]
        entries.append(dict(key=_param_key("PBOOL" if section == "xml_boolean" else "PXML", len(entries), tag),
            kind=section, tag=tag, path=list(path), path_text=path_text, algorithm_stack=list(stack),
            identity=_xml_identity(root, extra.get("optional_parent_path", path)),
            parent_algorithm=stack[-1] if stack else None, default=_xml_number(text), text_default=text, **extra))

    for node in nodes:
        elem = node[0]
        text = (elem.text or "").strip()
        if len(elem) == 0 and text:
            if _NUMERIC_RE.fullmatch(text):
                add(node, elem.tag, text, "xml_numeric")
            elif text.lower() in {"true", "false"}:
                add(node, elem.tag, text, "xml_boolean")
    existing = {(tuple(e["path"][:-1]), e["tag"]) for e in groups["xml_numeric"]}
    instances = set()
    for node in nodes:
        elem, path, path_text, stack = node
        scope = elem.get("type") if elem.tag == "algorithm" else elem.tag
        if elem.tag == "algorithm" and elem.get("instance"):
            if elem.get("instance") in instances:
                continue
            instances.add(elem.get("instance"))
        if elem is not root and elem.tag != "algorithm" and len(path) != 1:
            continue
        for tag, default in optional_xml_numeric_defaults.get(scope, {}).items():
            if (path, tag) not in existing:
                text = f"{float(default):.1f}" if float(default).is_integer() else f"{float(default):g}"
                add((elem, (), f"{path_text}/{tag}[optional]", stack), tag, text, "xml_numeric",
                    optional_parent_path=list(path), optional_insert=True)
                existing.add((path, tag))
    # Plugin settings are named root children, even when absent from the source XML
    plugins = {elem.text.strip() for elem in root if elem.tag.endswith("Plugin") and elem.text and not len(elem)}
    plugins.add("ShowerProfilePlugin")
    for scope in sorted(plugins & optional_xml_numeric_defaults.keys()):
        if root.find(scope) is not None:
            continue
        for tag, default in optional_xml_numeric_defaults[scope].items():
            add((root, (), f"{root.tag}/{scope}/{tag}[optional]", (scope,)), tag, str(default), "xml_numeric",
                optional_parent_path=[], optional_parent_tag=scope, optional_insert=True)
    return dict(schema_version=CATALOG_SCHEMA_VERSION, settings_xml=str(settings_xml),
                settings_xml_sha256=runtime.file_identity(settings_xml)["sha256"], **groups)


# Build the XML and wrapper catalog from the explicit steering
def build_catalog(settings_xml, *, optional_xml_numeric_defaults, wrapper_defaults):
    wrapper = [dict(key=f"PWRAP_{_PARAM_SAFE_RE.sub('_', name)}" + (f"_{i}" if len(values) > 1 else ""),
                    kind="wrapper_numeric", name=name, index=i, default=float(value))
               for name, values in wrapper_defaults.items() for i, value in enumerate(values)]
    return {**catalog_pandora_settings(settings_xml, optional_xml_numeric_defaults=optional_xml_numeric_defaults),
            "wrapper_numeric": wrapper}


# Index the XML and wrapper catalog by parameter key
def _entry_by_key(catalog):
    return {str(entry["key"]): entry for section in ("xml_numeric", "xml_boolean", "wrapper_numeric")
            for entry in catalog.get(section, [])}


# Decode numeric and boolean XML text
def _xml_number(text):
    boolean = {"true": 1.0, "false": 0.0}
    return boolean.get(str(text).strip().lower(), _safe_float(text))


# Format one parameter consistently for trial staging and source preserving push
def _xml_value(entry, value):
    if entry["kind"] == "xml_boolean":
        return "true" if float(value) >= 0.5 else "false"
    if entry.get("tune_dtype", "").lower() == "int":
        value = round(float(value))
        if re.fullmatch(r"[-+]?\d+", entry.get("text_default", "").strip()):
            return str(value)
    return str(float(value))


# Resolve fitted XML parameters, including absent optional children
def _xml_targets(root, catalog, config):
    for key, entry in _entry_by_key(catalog).items():
        if key not in config or not entry["kind"].startswith("xml_"):
            continue
        path = tuple(entry.get("optional_parent_path", entry["path"]))
        try:
            identity = _xml_identity(root, path)
        except IndexError as exc:
            raise ValueError(f"Pandora XML structure changed for {key}") from exc
        if identity != entry["identity"]:
            raise ValueError(f"Pandora XML structure changed for {key}")
        if entry.get("optional_insert"):
            path = tuple(entry["optional_parent_path"])
            parent = _element_children_path(root, path)
            if entry.get("optional_parent_tag"):
                parent = root.find(entry["optional_parent_tag"])
                path = (list(root).index(parent),) if parent is not None else ()
            elem = parent.find(entry["tag"]) if parent is not None else None
            if elem is not None:
                path = (*path, list(parent).index(elem))
        else:
            elem = _element_children_path(root, path)
        yield key, entry, elem, path, float(config[key])


# Compute explicit fixed settings, including values absent from the source XML
def _fixed_config(catalog):
    return {e["key"]: e["tune_bounds"][0] for e in _entry_by_key(catalog).values() if e.get("sampling_mode") == "fixed"}


# Complete fixed settings and reject invalid values before writing any reconstruction steering
def _validate_config(catalog, config):
    config = {**_fixed_config(catalog), **config}
    entries = _entry_by_key(catalog)
    unknown = config.keys() - entries.keys()
    if unknown:
        raise ValueError(f"Unknown Pandora parameters: {sorted(unknown)}")
    for key, value in config.items():
        entry, value = entries[key], float(value)
        if not math.isfinite(value):
            raise ValueError(f"Nonfinite Pandora parameter: {key}")
        if entry["kind"] == "xml_boolean" and value not in (0, 1):
            raise ValueError(f"Pandora boolean parameter {key} must be zero or one")
        if "tune_bounds" in entry:
            low, high = entry["tune_bounds"]
            if not low <= value <= high or (entry["tune_dtype"] == "int" and not value.is_integer()):
                raise ValueError(f"Pandora parameter {key} violates its {entry['tune_dtype']} bounds {low}, {high}: {value}")
    for ordered in catalog.get("ordered_parameters", []):
        for left, right in zip(ordered[:-1], ordered[1:], strict=True):
            if float(config.get(left, entries[left]["default"])) > float(config.get(right, entries[right]["default"])):
                raise ValueError(f"Pandora parameter ordering requires {left} <= {right}")
    return config


# Apply trial values to XML while preserving an unchanged source default
def _apply_config_to_xml(tree, catalog, config):
    config = _validate_config(catalog, config)
    updates = 0
    for _key, entry, elem, path, value in list(_xml_targets(tree.getroot(), catalog, config)):
        old = entry["default"] if elem is None else _xml_number(elem.text)
        if (old is not None and float(old).hex() == value.hex()
                and (elem is not None or entry.get("sampling_mode") != "fixed")):
            continue
        if elem is None:
            parent = _element_children_path(tree.getroot(), path)
            if entry.get("optional_parent_tag") and parent is tree.getroot():
                plugin = parent.find(entry["optional_parent_tag"])
                parent = ET.SubElement(parent, entry["optional_parent_tag"]) if plugin is None else plugin
            elem = ET.SubElement(parent, entry["tag"])
        text = _xml_value(entry, value)
        if elem.text != text:
            elem.text = text
            updates += 1
    return updates


# Compute contiguous wrapper overrides from fitted catalog entries
def _wrapper_overrides_from_config(catalog, config):
    config = _validate_config(catalog, config)
    overrides = {}
    for key, entry in _entry_by_key(catalog).items():
        if key in config and entry["kind"] == "wrapper_numeric":
            overrides.setdefault(entry["name"], {})[entry["index"]] = str(float(config[key]))
    for name, values in overrides.items():
        if sorted(values) != list(range(len(values))):
            raise ValueError(f'Wrapper override "{name}" has non-contiguous tuned indices {sorted(values)}')
    return {name: [values[i] for i in range(len(values))] for name, values in overrides.items()}


# Record the catalog subset needed to interpret and push this fitted configuration
def _card_config(config, catalog):
    config = {**_fixed_config(catalog), **config}
    selected = {section: [entry for entry in catalog.get(section, []) if entry["key"] in config]
                for section in ("xml_numeric", "xml_boolean", "wrapper_numeric")}
    return copy.deepcopy(dict(schema_version=1, parameters=config, tables={}, catalog=dict(
        schema_version=catalog.get("schema_version"), settings_xml=catalog.get("settings_xml"),
        ordered_parameters=[row for row in catalog.get("ordered_parameters", []) if all(key in config for key in row)],
        **selected)))


# Compute a unique existing path from a set of Pandora push candidates
def _unique_path(paths: list[pathlib.Path], label: str) -> pathlib.Path:
    unique = sorted({path.resolve() for path in paths if path.is_file()})
    if len(unique) != 1:
        raise ValueError(f"Pandora push expected one {label}, found {len(unique)}")
    return unique[0]


# Resolve the Pandora XML target using catalog metadata or matching parameter keys
def _resolve_push_xml(target: pathlib.Path, summary: dict, xml_keys: set[str]) -> pathlib.Path:
    if target.is_file():
        if target.suffix.lower() != ".xml":
            raise ValueError(f'Pandora XML parameters cannot be pushed into "{target}"')
        return target
    catalog = (summary.get("card_config") or {}).get("catalog") or {}
    recorded = pathlib.Path(str(catalog.get("settings_xml") or "")).name
    if recorded:
        candidates = list(target.rglob(recorded))
        if candidates:
            return _unique_path(candidates, f'XML file named "{recorded}"')
    direct = target / "PandoraSettings.xml"
    if direct.is_file():
        return direct.resolve()

    candidates = []
    for path in target.rglob("*.xml"):
        try:
            candidate_catalog = catalog_pandora_settings(path, optional_xml_numeric_defaults={})
        except (ET.ParseError, OSError):
            continue
        if xml_keys <= set(_entry_by_key(candidate_catalog)):
            candidates.append(path)
    return _unique_path(candidates, "XML file matching the fitted catalog")


# Resolve the Python steering file which owns Pandora wrapper parameters
def _resolve_push_python(target: pathlib.Path) -> pathlib.Path:
    if target.is_file():
        if target.suffix.lower() != ".py":
            raise ValueError(f'Pandora wrapper parameters cannot be pushed into "{target}"')
        return target
    return _unique_path(list(target.rglob("run_reco_pandora.py")), "run_reco_pandora.py")


# Compute XML element-content byte spans indexed by ElementTree child paths
def _xml_content_spans(text):
    data, parser = text.encode("utf-8"), xml.parsers.expat.ParserCreate()
    stack, spans = [], {}

    # Keep each open element's path, next child index and content start together
    def start_element(name, attributes):
        path = (*stack[-1][0], stack[-1][1]) if stack else ()
        if stack:
            stack[-1][1] += 1
        end = data.find(b">", parser.CurrentByteIndex)
        if end < 0:
            raise ValueError("Unterminated XML start tag")
        stack.append([path, 0, end + 1])

    # Finish the content span when the corresponding element closes
    def end_element(name):
        path, _, start = stack.pop()
        spans[path] = start, parser.CurrentByteIndex

    parser.StartElementHandler, parser.EndElementHandler = start_element, end_element
    parser.Parse(data, True)
    return spans


# Require recorded XML identities and infer only named scalar wrapper parameters
def _push_catalog(summary, settings_xml):
    recorded = (summary.get("card_config") or {}).get("catalog")
    if isinstance(recorded, dict):
        return copy.deepcopy(recorded)
    if settings_xml is not None:
        raise ValueError("Pandora XML push requires the recorded parameter catalog")
    catalog = dict(schema_version=CATALOG_SCHEMA_VERSION, xml_numeric=[], xml_boolean=[])
    catalog["wrapper_numeric"] = [dict(key=key, kind="wrapper_numeric", name=key.removeprefix("PWRAP_"), index=0)
                                  for key in summary["config"] if str(key).startswith("PWRAP_")]
    return catalog


# Apply non-overlapping source replacements from the final offset backwards
def _apply_source_replacements(source: Any, replacements: list[tuple]) -> Any:
    for start, end, replacement in sorted(replacements, reverse=True):
        source = source[:start] + replacement + source[end:]
    return source


# Render fitted values into an XML file while preserving all unrelated text
def _render_xml_push(path, catalog, config):
    source = path.read_text(encoding="utf-8")
    data, root = source.encode("utf-8"), ET.fromstring(source)
    spans = _xml_content_spans(source)
    replacements, insertions, plugins, old_values = [], {}, {}, {}
    for key, entry, elem, child_path, value in _xml_targets(root, catalog, config):
        old = entry["default"] if elem is None else _xml_number(elem.text)
        if old is None:
            raise ValueError(f'Pandora XML parameter "{key}" is not numeric')
        old_values[key] = old
        text = _xml_value(entry, value)
        if elem is None:
            if entry.get("sampling_mode") == "fixed" or float(old).hex() != value.hex():
                if entry.get("optional_parent_tag") and not child_path:
                    plugins.setdefault(entry["optional_parent_tag"], []).append((entry["tag"], text))
                else:
                    insertions.setdefault(child_path, []).append((entry["tag"], text))
            continue
        start, end = spans[child_path]
        content = data[start:end]
        left, right = len(content) - len(content.lstrip()), len(content.rstrip())
        if content[left:right] != text.encode():
            replacements.append((start + left, start + right, text.encode()))
    for plugin, values in plugins.items():
        text = "".join(f"\n        <{tag}>{value}</{tag}>" for tag, value in values) + "\n    "
        insertions.setdefault((), []).append((plugin, text))
    for parent_path, values in insertions.items():
        start, end = spans[parent_path]
        if start == end and data[start - 2:start] == b"/>":
            line = data[data.rfind(b"\n", 0, start) + 1:start]
            indent = re.match(rb"[ \t]*", line).group()
            newline = b"\r\n" if b"\r\n" in data else b"\n"
            addition = b"".join(newline + indent + b"    " + f"<{tag}>{value}</{tag}>".encode() for tag, value in values)
            tag = _element_children_path(root, parent_path).tag
            replacements.append((start - 2, start, b">" + addition + newline + indent + f"</{tag}>".encode()))
            continue
        content = data[start:end]
        trailing = re.search(rb"(\r?\n)([ \t]*)$", content)
        newline, indent = trailing.groups() if trailing else (b"\n", b"")
        child = re.search(rb"\r?\n([ \t]+)<", content)
        indent = child.group(1) if child else indent + b"    "
        at = end - (len(trailing.group(0)) if trailing else 0)
        addition = b"".join(newline + indent + f"<{tag}>{value}</{tag}>".encode() for tag, value in values)
        replacements.append((at, at, addition))
    rendered = _apply_source_replacements(data, replacements).decode("utf-8")
    for key, entry, elem, _, value in _xml_targets(ET.fromstring(rendered), catalog, config):
        checked = entry["default"] if elem is None else _xml_number(elem.text)
        if checked is None or float(checked).hex() != value.hex():
            raise ValueError(f'Pandora XML push validation failed for "{key}"')
    return rendered, old_values


# Compute source replacements from AST byte offsets without changing surrounding text
def _ast_replacement(source, node, value):
    lines = source.encode("utf-8").splitlines(keepends=True)
    start = sum(map(len, lines[:node.lineno - 1])) + node.col_offset
    end = sum(map(len, lines[:node.end_lineno - 1])) + node.end_col_offset
    return start, end, value.encode("utf-8")


# Render Pandora wrapper values into their existing Python dictionary literals
def _render_python_push(path, catalog, config):
    source = path.read_text(encoding="utf-8")
    by_name = {}
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, ast.Dict):
            for key, value in zip(node.keys, node.values, strict=True):
                if isinstance(key, ast.Constant) and isinstance(key.value, str) and isinstance(value, (ast.List, ast.Tuple)):
                    by_name.setdefault(key.value, []).append(value)
    replacements, old_values = [], {}
    for key, entry in _entry_by_key(catalog).items():
        if key not in config or entry["kind"] != "wrapper_numeric":
            continue
        name, index = entry["name"], entry["index"]
        matches = by_name.get(name, [])
        if len(matches) != 1:
            raise ValueError(f'Pandora wrapper parameter "{name}" occurs {len(matches)} times')
        if index >= len(matches[0].elts):
            raise ValueError(f'Pandora wrapper parameter "{name}" has no index {index}')
        node = matches[0].elts[index]
        try:
            old = ast.literal_eval(node)
        except (ValueError, TypeError, SyntaxError) as exc:
            raise ValueError(f'Pandora wrapper parameter "{name}[{index}]" is not numeric') from exc
        if isinstance(old, bool) or not isinstance(old, (str, int, float)) or _safe_float(old) is None:
            raise ValueError(f'Pandora wrapper parameter "{name}[{index}]" is not numeric')
        old_values[key] = float(old)
        old_text = ast.get_source_segment(source, node)
        value = str(float(config[key]))
        value = old_text[0] + value + old_text[-1] if old_text[:1] in {'"', "'"} else value
        replacements.append(_ast_replacement(source, node, value))
    rendered = _apply_source_replacements(source.encode("utf-8"), replacements).decode("utf-8")
    ast.parse(rendered)
    return rendered, old_values


# Compute the trial output directory
def _trial_dir(*, cdir: str, run_name: str, tunename: str) -> pathlib.Path:
    return pathlib.Path(cdir) / "runs" / "icetune" / run_name / "trials" / tunename


# Copy Pandora steering and add a trial-local override block
def _copy_and_patch_steering(*, source_path, target_path, override_path, data_dir):
    text = source_path.read_text(encoding="utf-8")
    assignments = [node for node in ast.parse(text).body if isinstance(node, ast.Assign)
                   and any(isinstance(t, ast.Name) and t.id == "dataFolder" for t in node.targets)]
    if len(assignments) != 1:
        raise RuntimeError(f"Expected one top-level Pandora dataFolder assignment in {source_path}, found {len(assignments)}")
    replacement = _ast_replacement(text, assignments[0].value, repr(str(data_dir.resolve()) + os.sep))
    text = _apply_source_replacements(text.encode("utf-8"), [replacement]).decode("utf-8")
    marker = "    TopAlg += [pandora]\n"
    if marker not in text:
        raise RuntimeError(f"Could not find Pandora TopAlg marker in {source_path}")
    patch = f"""
    # Apply icetune wrapper parameters to this trial
    import json as _icetune_json
    with open({str(override_path)!r}, "r", encoding="utf-8") as _icetune_file:
        for _icetune_key, _icetune_values in _icetune_json.load(_icetune_file).items():
            pandora.Parameters[_icetune_key] = [str(_value) for _value in _icetune_values]
"""
    target_path.write_text(text.replace(marker, patch + marker, 1), encoding="utf-8")


# Link read-only Pandora XML inputs into worker-local trial scratch
def _stage_xml_inputs(tree: ET.ElementTree, *, source_dir: pathlib.Path, target_dir: pathlib.Path) -> None:
    for elem in tree.getroot().iter("HistogramFile"):
        value = str(elem.text or "").strip()
        relative = pathlib.Path(value)
        if not value or relative.is_absolute() or ".." in relative.parts:
            raise ValueError(f"Pandora HistogramFile must be relative to the frozen run directory: {value}")
        source = source_dir / relative
        if not source.is_file():
            raise FileNotFoundError(f"Pandora XML input does not exist: {source}")
        target = target_dir / relative
        ensure_dir(target.parent)
        if not target.exists():
            target.symlink_to(source.resolve())


# Run Pandora with bounded duration and persistent command diagnostics
def _run_command(cmd, *, cwd, log_path, max_t, setup_script=None):
    ensure_dir(log_path.parent)
    started = time.time()
    run_cmd = list(map(str, cmd))
    command_text = shlex.join(run_cmd)
    environment = None
    if setup_script is not None:
        command_text = f"source {shlex.quote(str(setup_script))} && cd {shlex.quote(str(cwd))} && {command_text}"
        run_cmd = ["/bin/bash", "--noprofile", "--norc", "-c", command_text]
        environment = runtime.command_environment()
    with log_path.open("w", encoding="utf-8", errors="replace") as log:
        log.write(f"$ {command_text}\ncwd={cwd}\n\n")
        log.flush()
        result = run_process(run_cmd, cwd=str(cwd), stdout=log, timeout=max_t, env=environment)
        log.write(f"\nexit_code={result.returncode}\nwall_s={time.time() - started:.3f}\n")
    if result.returncode:
        output = log_path.read_text(encoding="utf-8", errors="replace")
        diagnostic = "\n".join(s for s in output.splitlines() if not s.startswith(("$ ", "cwd=", "exit_code=", "wall_s=")))
        cause = command_output_root_cause(diagnostic)
        tail = command_output_tail(output, max_lines=40, max_chars=8000)
        if cause and cause not in tail:
            lines = diagnostic.splitlines()
            index = next(i for i, line in enumerate(lines) if cause in line)
            context = "\n".join(lines[max(0, index - 5):index + 16])[:8000]
            tail = f"{context}\n...\n{tail}"
        raise RuntimeError(f"Command failed with exit code {result.returncode}, root_cause={cause}, log={log_path}\n{tail}")


# Validate the complete objective steering before reading or reconstructing events
def _validate_pf_spec(spec):
    for key in ("matching_max_angle_rad", "confusion_energy_fraction_target", "epsilon",
                "energy_resolution_target", "energy_bias_tolerance", "energy_bias_target", "multiplicity_target",
                "matched_momentum_response_target"):
        if not math.isfinite(float(spec[key])) or float(spec[key]) <= 0.0:
            raise ValueError(f"pflow {key} must be finite and positive")
    if spec["energy_bias_tolerance"] >= 1.0:
        raise ValueError("pflow energy_bias_tolerance must be smaller than one")
    if not 0.0 < float(spec["max_abs_parton_costheta"]) <= 1.0 or spec["matching_max_angle_rad"] > math.pi:
        raise ValueError("Invalid pflow angular selection")
    targets = spec["matched_energy_response_target"]
    if set(targets) != set(diagnostics.PF_RESIDUAL_CLASS_NAMES) or any(
        not math.isfinite(float(value)) or float(value) <= 0.0 for value in targets.values()):
        raise ValueError("Matched energy response targets must be positive for every PF class")
    for block in ("gen", "reco"):
        for field in ("energy", "theta_deg"):
            low, high = (spec[f"{block}_final_state_{edge}_{field}"] for edge in ("min", "max"))
            if not math.isfinite(low) or low < 0.0 or (high is not None and (not math.isfinite(high) or high < low)):
                raise ValueError(f"Invalid {block} {field} acceptance")
            if field == "theta_deg" and (high is None or high > 180.0 or high <= low):
                raise ValueError(f"Invalid {block} polar angle acceptance")
    if spec["gen_status"] != 1:
        raise ValueError("pflow needs stable generator truth")


# Compute the objective from the Pandora tuning card
def pflow_objective_spec(*, param_selection, param_loss):
    generator, fiducial = (param_selection[key] for key in ("generator", "fiducial"))
    spec = dict(tree="events", truth_collection="MCParticles", reco_collection="PandoraPFANewPFOs",
                gen_status=fiducial["gen"]["gen_status"], max_abs_parton_costheta=generator["max_abs_parton_costheta"],
                matching_max_angle_rad=param_selection["matching_max_angle_rad"], **copy.deepcopy(param_loss))
    for block in ("gen", "reco"):
        spec[f"{block}_excluded_abs_pdgs"] = tuple(sorted(abs(int(pdg)) for pdg in fiducial[block]["excluded_abs_pdgs"]))
        for field, key in (("energy", "final_state_energy_gev"), ("theta_deg", "final_state_theta_deg")):
            spec.update({f"{block}_final_state_{edge}_{field}": value
                         for edge, value in zip(("min", "max"), fiducial[block][key], strict=True)})
    _validate_pf_spec(spec)
    return spec


# Compute the EDM4hep fields needed for truth selection and particle flow scoring
def _pflow_branch_map(spec):
    fields = {"truth": {"pdg": "PDG", "status": "generatorStatus", "charge": "charge", "mass": "mass"},
        "reco": {"pdg": "PDG", "energy": "energy"}}
    return {f"{side}_{name}": f"{spec[side + '_collection']}/{spec[side + '_collection']}.{field}"
        for side, mapping in fields.items()
        for name, field in {**mapping, **{f"p{axis}": f"momentum.{axis}" for axis in "xyz"}}.items()}


# Compute the highest energy hard light quark direction in each reconstruction event
def _pf_parton_acceptance_mask(spec, arrays, size):
    import awkward as ak

    branches = _pflow_branch_map(spec)
    pdg, status, mass = (arrays[branches[f"truth_{key}"]] for key in ("pdg", "status", "mass"))
    px, py, pz = (arrays[branches[f"truth_p{axis}"]] for axis in "xyz")
    p2 = px * px + py * py + pz * pz
    e2 = p2 + mass * mass
    mask = (status == 23) & (abs(pdg) >= 1) & (abs(pdg) <= 3) & (p2 > 0) & np.isfinite(e2)
    selected = ak.argmax(e2[mask], axis=1, keepdims=True, mask_identity=True)
    cos_theta = ak.to_numpy(ak.fill_none(ak.firsts((pz[mask] / np.sqrt(p2[mask]))[selected]), np.nan))
    if len(cos_theta) != size:
        raise ValueError("Inconsistent MCParticles event count")
    return np.abs(cos_theta) < spec["max_abs_parton_costheta"], ~np.isfinite(cos_theta)


# Compute stable PF classes from particle labels and optional truth charge
def _pf_particle_classes(pdg, charge=None):
    labels = np.abs(np.asarray(pdg, dtype=int))
    classes = np.full(labels.size, "other", dtype=object)
    for name, ids in (("photon", (22,)), ("electron", (11,)), ("muon", (13,)), ("charged_hadron", (211, 321, 2212)),
        ("neutral_hadron", (130, 2112, 310, 3122))):
        classes[np.isin(labels, ids)] = name
    if charge is not None:
        fallback = (classes == "other") & ~np.isin(labels, PFLOW_NONOBSERVABLE_TRUTH_ABS_PDGS)
        classes[fallback & (np.abs(charge) > 1.0e-9)] = "charged_hadron"
    return classes


# Compute the fiducial selection in energy and polar angle
def _pf_final_state_acceptance(spec, block):
    low, high = (spec[f"{block}_final_state_{edge}_theta_deg"] for edge in ("min", "max"))
    return dict(min_energy_gev=spec[f"{block}_final_state_min_energy"], max_energy_gev=spec[f"{block}_final_state_max_energy"],
                min_theta_deg=low, max_theta_deg=high,
                min_costheta=math.cos(math.radians(high)), max_costheta=math.cos(math.radians(low)))


# Read and centrally validate one visible particle collection
def _pf_particles(spec, arrays, index, side):
    branches = _pflow_branch_map(spec)
    p = {key.removeprefix(side + "_"): np.asarray(arrays[branch][index]).reshape(-1) for key, branch in branches.items()
        if key.startswith(side + "_")}
    if len({value.size for value in p.values()}) != 1 or any(not np.all(np.isfinite(v)) for v in p.values()):
        raise ValueError(f"Invalid {side} particle data in event {index}")
    momentum = np.column_stack([p[f"p{axis}"] for axis in "xyz"]).astype(float)
    magnitude = np.linalg.norm(momentum, axis=1)
    energy = np.sqrt(magnitude**2 + p["mass"] ** 2) if side == "truth" else p["energy"].astype(float)
    tolerance = np.sqrt(np.finfo(momentum.dtype if side == "truth" else p["energy"].dtype).eps)
    spacelike = (energy < magnitude) & ~np.isclose(energy, magnitude, rtol=tolerance, atol=0.0)
    if np.any(energy < 0.0) or np.any(spacelike) or (side == "truth" and np.any(p["mass"] < 0.0)):
        raise ValueError(f"Unphysical {side} four momentum in event {index}")
    unit = np.divide(momentum, magnitude[:, None], out=np.zeros_like(momentum), where=magnitude[:, None] > 0.0)
    block = "gen" if side == "truth" else "reco"
    acceptance = _pf_final_state_acceptance(spec, block)
    selected = ((magnitude > 0.0) & (energy >= acceptance["min_energy_gev"]) & (unit[:, 2] > acceptance["min_costheta"])
        & (unit[:, 2] < acceptance["max_costheta"]) & ~np.isin(np.abs(p["pdg"]), spec[f"{block}_excluded_abs_pdgs"]))
    if acceptance["max_energy_gev"] is not None:
        selected &= energy <= acceptance["max_energy_gev"]
    if side == "truth":
        selected &= p["status"] == spec["gen_status"]
        if np.any(selected & np.isin(np.abs(p["pdg"]), PFLOW_NONOBSERVABLE_TRUTH_ABS_PDGS)):
            raise ValueError(f"Nonobservable truth PDGs in event {index}, use stable post decay truth")
    classes = _pf_particle_classes(p["pdg"], p.get("charge"))
    if np.any(selected & (classes == "other")):
        raise ValueError(f"Unsupported visible {side} particle PDGs in event {index}")
    energy, momentum, unit, magnitude, classes = (v[selected] for v in (energy, momentum, unit, magnitude, classes))
    pt = np.linalg.norm(momentum[:, :2], axis=1)
    return dict(energy=energy, momentum=magnitude, vector=momentum, unit=unit, costheta=unit[:, 2],
                pt=pt, eta=np.arcsinh(momentum[:, 2] / pt), classes=classes,
                class_index=np.array([diagnostics.PF_RESIDUAL_CLASS_NAMES.index(c) for c in classes], dtype=int))


# Compute class agnostic angular matches with maximum cardinality before angular cost
def _pf_angular_matched_pairs(truth, reco, max_angle):
    from scipy.optimize import linear_sum_assignment

    if not len(truth["energy"]) or not len(reco["energy"]):
        return np.empty((0, 2), dtype=int)
    # Species breaks only exact kinematic sorting ties, never changes angular costs
    ti, ri = (np.lexsort((p.get("class_index", np.zeros(len(p["energy"]))),
                         p.get("momentum", np.zeros(len(p["energy"]))), *p["unit"].T, p["energy"])) for p in (truth, reco))
    angle = np.arccos(np.clip(truth["unit"][ti] @ reco["unit"][ri].T, -1.0, 1.0))
    # One forbidden edge costs more than any sum of allowed angles, maximizing cardinality first
    forbidden = (min(angle.shape) + 1) * math.pi
    energy = [p["energy"][order] / p["energy"].sum() for p, order in ((truth, ti), (reco, ri))]
    balance = ((energy[0][:, None] - energy[1]) / (energy[0][:, None] + energy[1]))**2
    # Resolve angular ties at floating point precision using gain invariant energy fractions
    cost = angle + np.spacing(forbidden) * balance
    rows, cols = linear_sum_assignment(np.where(angle <= max_angle, cost, forbidden))
    valid = angle[rows, cols] <= max_angle
    return np.column_stack((ti[rows[valid]], ri[cols[valid]]))


# Compute truth particle matching diagnostics without changing the assigned pairs
def _pf_matching_profile(truth, reco, pairs, max_angle):
    ti, ri = pairs.T
    angles = np.arccos(np.clip(truth["unit"] @ reco["unit"].T, -1.0, 1.0))
    profile = dict(assigned=np.zeros(len(truth["energy"]), dtype=bool), candidates=np.sum(angles <= max_angle, axis=1))
    profile["assigned"][ti] = True
    for key in ("angle", "dr", "energy_residual"):
        profile[key] = np.full(len(truth["energy"]), np.nan)
    profile["angle"][ti] = angles[ti, ri]
    profile["energy_residual"][ti] = (reco["energy"][ri] - truth["energy"][ti]) / truth["energy"][ti]
    phi = [np.arctan2(p["vector"][:, 1], p["vector"][:, 0]) for p in (truth, reco)]
    dphi = phi[0][ti] - phi[1][ri]
    dphi = np.arctan2(np.sin(dphi), np.cos(dphi))
    profile["dr"][ti] = np.hypot(truth["eta"][ti] - reco["eta"][ri], dphi)
    return profile


# Compute physical energy response and PF losses independent of the event energy scale
def _evaluate_pflow_event(spec, arrays, event_index, *, collect_plot_details=True):
    truth, reco = (_pf_particles(spec, arrays, event_index, side) for side in ("truth", "reco"))
    energies = [p["energy"] for p in (truth, reco)]
    total = float(energies[0].sum())
    if total <= 0.0:
        return {"valid": False}
    four = [np.column_stack((p["energy"], p["vector"])).sum(axis=0) for p in (truth, reco)]
    pairs = _pf_angular_matched_pairs(truth, reco, spec["matching_max_angle_rad"])
    ti, ri = pairs.T
    tc, rc = truth["class_index"], reco["class_index"]
    # Penalize each class in each event so deficits cannot cancel excesses
    counts = [np.bincount(c, minlength=len(diagnostics.PF_RESIDUAL_CLASS_NAMES)) for c in (tc, rc)]
    multiplicity_loss = float(np.sum(((counts[1] - counts[0]) / spec["multiplicity_target"])**2))
    lost, fake = (np.setdiff1d(np.arange(len(c)), matched) for c, matched in ((tc, ti), (rc, ri)))
    confusion = np.zeros((len(diagnostics.PF_CONFUSION_ROW_NAMES), len(diagnostics.PF_CONFUSION_COLUMN_NAMES)))
    ideal = np.zeros_like(confusion)
    np.add.at(ideal, (tc, tc), energies[0] / total)
    np.add.at(confusion, (rc[ri], tc[ti]), energies[0][ti] / total)
    np.add.at(confusion[-1], tc[lost], energies[0][lost] / total)
    np.add.at(confusion[:, -1], rc[fake], energies[1][fake] / total)
    closure, delta = (four[1] - four[0]) / four[0][0], confusion - ideal
    reco_total = max(float(four[1][0]), spec["epsilon"])
    scale = four[0][0] / reco_total
    # Compare energy sharing at equal event energy without altering the physical response or plots
    energy = reco["energy"][ri] * scale
    pid = delta.copy()
    pid[:, -1] = np.bincount(rc[fake], weights=reco["energy"][fake], minlength=len(pid)) / reco_total
    targets = np.array([spec["matched_energy_response_target"][c] for c in truth["classes"][ti]])
    # Symmetric pair energy prevents a soft truth particle from producing an unbounded relative response
    response = (energy - truth["energy"][ti]) / targets
    pair_energy = truth["energy"][ti] + energy
    pair_weight = 2.0 / (pair_energy * total)
    response = pair_weight * response**2
    momentum = (reco["vector"][ri] * scale - truth["vector"][ti]) / spec["matched_momentum_response_target"]
    momentum = pair_weight * np.sum(momentum**2, axis=1)
    event = dict(
        valid=True, weight=1.0, confusion_loss=float(np.sum((pid / spec["confusion_energy_fraction_target"])**2)),
        response_loss=float(response.sum()), momentum_loss=float(momentum.sum()), multiplicity_loss=multiplicity_loss,
        four_momentum_residual=closure, visible_response=float(four[1][0] / four[0][0]),
        truth_energy=float(four[0][0]), reco_energy=float(four[1][0]), confusion=confusion, confusion_ideal=ideal,
        **{f"n_{name}": len(items) for name, items in (
            ("truth_visible", tc), ("reco_visible", rc), ("matched", ti), ("lost", lost), ("fake", fake))})
    if collect_plot_details:
        event["loss_profile"] = diagnostics.loss_particles(truth, reco, pairs, pid=pid, counts=counts,
                                                           response=response, momentum=momentum, spec=spec)
        event["object_profiles"] = {c: {} for c in diagnostics.PF_CLASS_NAMES}
        event["confusion_profile"] = {}
        matching = _pf_matching_profile(truth, reco, pairs, spec["matching_max_angle_rad"])
        for side, p, indices in (("truth", truth, ti), ("reco", reco, ri)):
            matched = np.zeros(len(p["energy"]), dtype=bool)
            matched[indices[tc[ti] == rc[ri]]] = True
            fields = {key: p[key] for key in ("energy", "momentum", "costheta", "pt", "eta")}
            fields["matched"] = matched
            if side == "truth":
                fields.update(matching)
            for c in diagnostics.PF_CLASS_NAMES:
                mask = p["classes"] == c
                event["object_profiles"][c].update({f"{side}_{key}": value[mask] for key, value in fields.items()})
            classes = np.array([diagnostics.PF_CLASS_NAMES.index(c) for c in p["classes"]], dtype=int)
            for prefix, selected in ((side, slice(None)), ("match_" + side, indices)):
                event["confusion_profile"].update({f"{prefix}_energy": p["energy"][selected],
                                                   f"{prefix}_class": classes[selected]})
    return event


# Separate arithmetic mean calibration from the scale invariant PF risk Q
def _pf_loss(spec, *, mean, variance, confusion, response, momentum, multiplicity):
    target, tolerance = spec["energy_bias_target"], spec["energy_bias_tolerance"]
    calibrated = 1.0 - tolerance <= mean <= 1.0 + tolerance
    bias_ratio = abs(mean - 1.0) / tolerance
    resolution = variance / (math.hypot(mean, spec["epsilon"]) * spec["energy_resolution_target"]) ** 2
    losses = dict(calibration=((mean - 1.0) / target)**2, resolution=resolution,
                  confusion=confusion, response=response, momentum=momentum, multiplicity=multiplicity)
    # L = ((mean - 1) / target)**2 + log(1 + Q), with the gain minimum at mean = 1
    cost = losses["calibration"] + math.log1p(math.fsum((resolution, confusion, response, momentum, multiplicity)))
    calibration = dict(mean=mean, bias=mean - 1.0, target=target, tolerance=tolerance, passed=bool(calibrated))
    metrics = dict(pflow=float(cost), pf_calibration_passed=int(calibrated),
                   pf_calibration_bias_ratio=float(bias_ratio), pf_calibration_tolerance=float(tolerance),
                   pf_calibration_target=float(target), **{f"pf_loss_{key}": float(value) for key, value in losses.items()})
    return metrics, calibration


# Compute weighted sample moments and plots for either one dataset or their union
def _summarize_pf(spec, events, counts, *, collect_plot_details):
    weights = np.array([e["weight"] for e in events])
    if not events or weights.sum() <= 0.0:
        raise ValueError("No positive weight visible events for pflow")
    values = {key: np.asarray([e[key] for e in events]) for key in (
        "visible_response", "truth_energy", "reco_energy", "four_momentum_residual",
        "confusion", "confusion_ideal", "confusion_loss", "response_loss", "momentum_loss", "multiplicity_loss")}
    means = {key: np.average(value, axis=0, weights=weights) for key, value in values.items()}
    response, closure = values["visible_response"], values["four_momentum_residual"]
    mean = float(means["visible_response"])
    variance = float(np.average((response - mean) ** 2, weights=weights))
    metrics, calibration = _pf_loss(spec, mean=mean, variance=variance, confusion=means["confusion_loss"],
                                    response=means["response_loss"], momentum=means["momentum_loss"],
                                    multiplicity=means["multiplicity_loss"])
    detail = dict(kind="pflow", target=1.0, weights=weights, _events=events, _spec=spec,
        parton_acceptance=dict(max_abs_costheta=spec["max_abs_parton_costheta"], **counts),
        final_state_acceptance={b: _pf_final_state_acceptance(spec, b) for b in ("gen", "reco")},
        loss=copy.deepcopy(spec), calibration=calibration, cost=metrics["pflow"],
        loss_components={name.removeprefix("pf_loss_"): value for name, value in metrics.items()
                         if name.startswith("pf_loss_")})
    for prefix, samples in (("", response), ("truth_energy_", values["truth_energy"]),
                            ("visible_total_energy_", values["reco_energy"]), ("residual_", closure[:, 0])):
        detail[prefix + "values"] = samples
        detail[prefix + "stats"] = diagnostics.interval_stats(samples, weights=weights)
    matrix, ideal = means["confusion"], means["confusion_ideal"]
    detail["confusion"] = dict(row_names=diagnostics.PF_CONFUSION_ROW_NAMES, column_names=diagnostics.PF_CONFUSION_COLUMN_NAMES,
                               matrix=matrix, ideal=ideal, residual=matrix - ideal)
    response_metrics = dict(mean=mean, bias=mean - 1, mse=variance + (mean - 1) ** 2,
                           width_rel=math.sqrt(variance) / max(abs(mean), spec["epsilon"]),
                           mean90=detail["stats"]["mean"], sigma90=detail["stats"]["sigma"])
    fractions = dict(diagonal=np.trace(matrix[:-1, :-1]), lost=matrix[-1, :-1].sum(), fake=matrix[:-1, -1].sum())
    fractions["offdiagonal"] = matrix[:-1, :-1].sum() - fractions["diagonal"]
    metrics.update(pf_total_energy_residual_mean=mean - 1, pf_total_energy_residual_variance=variance,
                   pf_n_events=float(len(events)), **{f"pf_total_energy_response_{k}": v for k, v in response_metrics.items()},
                   **{f"pf_confusion_{k}_energy_fraction": float(v) for k, v in fractions.items()})
    metrics.update({f"pf_{name}_momentum_closure_mse": float(np.average(np.sum(closure[:, start:] ** 2, axis=1), weights=weights))
                    for name, start in (("four", 0), ("three", 1))})
    metrics.update({f"pf_n_events_{name}": float(counts[key]) for name, key in (
        ("input", "n_input"), ("parton_accepted", "n_accepted"), ("rejected_parton", "n_rejected"), ("missing_parton", "n_missing"))})
    metrics.update({f"pf_n_{name}": float(sum(e[f"n_{name}"] for e in events))
                    for name in ("truth_visible", "reco_visible", "matched", "lost", "fake")})
    metrics.update({f"pf_four_momentum_closure_{name}_mean": float(value) for name, value in
                    zip(("energy", "px", "py", "pz"), means["four_momentum_residual"], strict=True)})
    if not all(math.isfinite(value) for value in metrics.values()):
        raise ValueError("Nonfinite pflow result")
    if collect_plot_details:
        detail.update(diagnostics.plot_inputs(events, weights))
    return metrics, detail


# Score each input chunk once and summarize events in their original order
def _objective_pflow(spec, chunks):
    _validate_pf_spec(spec)
    collect = spec.get("_collect_plot_details", True)
    events, counts = [], dict.fromkeys(("n_input", "n_accepted", "n_rejected", "n_missing"), 0)
    for arrays in ([chunks] if isinstance(chunks, dict) else chunks):
        size = len(arrays[_pflow_branch_map(spec)["truth_pdg"]])
        mask, missing = _pf_parton_acceptance_mask(spec, arrays, size)
        for i in np.flatnonzero(mask):
            event = _evaluate_pflow_event(spec, arrays, i, collect_plot_details=collect)
            if event["valid"]:
                events.append(event)
        for key, value in zip(counts, (size, mask.sum(), (~mask).sum(), missing.sum()), strict=True):
            counts[key] += int(value)
    return _summarize_pf(spec, events, counts, collect_plot_details=collect)


# Stream the selected EDM4hep fields using Uproot's memory bounded iteration
def _reco_chunks(reco_path, objectives, *, step_size=None):
    import uproot

    if set(objectives) != {"pflow"}:
        raise ValueError("Pandora supports the pflow objective")
    spec = objectives["pflow"]
    branches = sorted(_pflow_branch_map(spec).values())
    options = {} if step_size is None else dict(step_size=step_size)
    for path in _as_path_list(reco_path):
        with uproot.open(path, handler=uproot.source.file.MultithreadedFileSource) as source:
            yield from source[spec["tree"]].iterate(branches, library="ak", how=dict, **options)


# Write a compact JSON output atomically
def _write_json(path, payload):
    ensure_dir(path.parent)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}.{time.time_ns()}")
    write_json_file(temporary, payload, indent=4, default=json_default)
    os.replace(temporary, path)


# Retire reconstruction files after evaluating their events
def _retire_reco_roots(reco_roots: Any) -> None:
    for path in _as_path_list(reco_roots):
        try:
            if path.exists():
                path.rename(path.with_name(path.name + "._old"))
        except Exception as exc:
            logger.warning(".retire_reco_roots: could not rename %s (%s)", path, exc)


# Compute JSON-sized objective details for trial summaries
def _compact_objective_details(details):
    heavy = {"values", "residual_values", "visible_total_energy_values",
             "truth_energy_values", "object_performance", "class_confusion", "_events", "_spec", "weights", "events", "plot_settings"}
    return {name: {key: value for key, value in detail.items() if key not in heavy} for name, detail in details.items()}


# Validate a nonnegative finite dataset weight before reconstruction
def _dataset_weight(datacard):
    value = float(datacard["weight"])
    if not math.isfinite(value) or value < 0.0:
        raise ValueError("Pandora dataset weight must be finite and nonnegative")
    return value


# Resolve one reproducible configuration for all positive weight datasets
def _trial_config(catalog, config, datacards):
    configs = [_validate_config(catalog, {**card.get("baseline", {}), **_fixed_config(catalog), **config})
               for card in datacards if _dataset_weight(card) > 0]
    if not configs:
        raise ValueError("Pandora requires a positive weight dataset")
    if any(values != configs[0] for values in configs[1:]):
        raise ValueError("Pandora datasets must share the same applied parameter configuration")
    return configs[0]


# Recompute combined statistics from the same weighted events used by the objective
def _combine_datasets(records):
    first = records[0]["details"]["pflow"]
    events, counts = [], dict.fromkeys(("n_input", "n_accepted", "n_rejected", "n_missing"), 0)
    weights = [_dataset_weight(record["datacard"]) for record in records]
    for record, weight in zip(records, weights, strict=True):
        detail = record["details"]["pflow"]
        if detail["_spec"] != first["_spec"]:
            raise ValueError("Combined Pandora datasets must share objective steering")
        if weight > 0.0:
            events.extend({**event, "weight": weight * event["weight"]} for event in detail["_events"])
            counts = {key: value + detail["parton_acceptance"][key] for key, value in counts.items()}
    metrics, detail = _summarize_pf(first["_spec"], events, counts, collect_plot_details="events" in first)
    metrics.update(pf_dataset_count=float(len(records)), pf_dataset_weight_sum=sum(weights))
    detail["datasets"] = [r["dataset"] for r in records]
    return metrics, {"pflow": detail}


# Adapt reconstruction to the common icetune driver interface
class PandoraDriver(SimulatorDriver):
    BOOTSTRAP_SCHEMA_VERSION = 1
    TRIAL_PREFIX = "PANDORA_icetune"

    # Compute the stable driver identifier
    @classmethod
    def driver_name(cls):
        return "PANDORA"

    # Recognize Pandora XML, boolean and wrapper parameters
    @classmethod
    def matches_summary(cls, summary):
        config = summary.get("config")
        return isinstance(config, dict) and bool(config) and all(
            str(key).startswith(("PXML_", "PBOOL_", "PWRAP_")) for key in config)

    # Compute the default Pandora runtime root
    @classmethod
    def default_library_path(cls):
        return os.environ.get("PANDORA_IO_PATH") or os.environ["HOME"]

    # Build or apply the Pandora worker environment
    def runtime_environment(self, *, cdir, libdir, python_version, apply=False):
        variables = {key: os.environ.get(key, "") for key in ("LD_LIBRARY_PATH", "PYTHONPATH")}
        if apply:
            os.environ.update(variables)
        for path in variables["PYTHONPATH"].split(os.pathsep):
            if path and path not in sys.path:
                sys.path.append(path)
        return variables

    # Keep mutable campaign state local to this driver instance
    def __init__(self):
        self.initialized, self.datacards = False, []

    # Build the detector parameter space from the selected campaign JSON
    @classmethod
    def build_tunesetup(cls, *, config, cdir, tune_default):
        from core.tune.drivers.pandora.tunesetup import common
        paths = {**config["param_paths"], "pandora_dir": os.environ.get("PANDORA_PFA_DIR", config["param_paths"]["pandora_dir"])}
        return common.setup(param_paths=paths, param_selection=config["param_selection"], param_loss=config["param_loss"],
            optional_xml_numeric_defaults=config["param_optional_xml_numeric_defaults"],
            wrapper_defaults=config["param_wrapper_defaults"], tuning_tables=config["param_tuning_tables"],
            ordered_parameters=config["param_ordered"])

    # Initialize without integration grids or reconstruction jobs
    def prepare_init_datacards(self, datacards):
        return []

    # Validate positive weight datasets and freeze their reconstruction inputs
    def prepare_tunesetup(self, *, tunesetup, args):
        for card in tunesetup.datacards:
            _dataset_weight(card)
            _validate_pf_spec(card["objectives"]["pflow"])
        tunesetup.datacards = [card for card in tunesetup.datacards if card["weight"] > 0.0]
        if not tunesetup.datacards:
            raise ValueError("Pandora requires a positive weight dataset")
        runtime.prepare_tunesetup(tunesetup, args.cdir)

    # Identify the complete frozen campaign without producing integration outputs
    def bootstrap_payload(self, *, runtime_sha256, datacards):
        from core.tune.cache import json_fingerprint

        fingerprint = json_fingerprint(dict(datacards=datacards, runtime_sha256=runtime_sha256,
                                            schema_version=self.BOOTSTRAP_SCHEMA_VERSION, simdriver="PANDORA"))
        return dict(archive_sha256=None, archive_size=0, archive_url=None, files=[], fingerprint=fingerprint,
                    kind="pandora", reusable_vgrids=[], schema_version=self.BOOTSTRAP_SCHEMA_VERSION)

    # Supply immutable campaign identification to Ray initialization
    def initialize_backend(self, *, args, tunesetup, mc_steer):
        if args.backend != "ray":
            raise ValueError(f'Pandora does not support backend "{args.backend}"')

        # Build an independent initialization record from the selected datacards
        def bootstrap():
            return self.bootstrap_payload(runtime_sha256=args.runtime_sha256, datacards=tunesetup.datacards)

        return icetune.bootstrap_callbacks(args=args, tunesetup=tunesetup, mc_steer=mc_steer, simdriver=self,
            fingerprint_getter=lambda: bootstrap()["fingerprint"], bootstrap_builder=bootstrap)

    # Store independent campaign data without running reconstruction
    def init_data(self, run_name, datacards, obs_module=None, cdir=None, pickle_dump=True):
        self.datacards, self.initialized = copy.deepcopy(datacards), True

    # Compute the initial parameter point without trial jobs
    def initialize(self, run_name, tunesetup, mc_steer, obs_module, cdir, init_force, max_t=3600, pickle_dump=True,
                   processes=1, rngseed=0):
        self.init_data(run_name, tunesetup.datacards)
        return self.get_initial_param(tunesetup.param_space, tunesetup.aux_param_space, cdir, mc_steer.get("tune_default", "TUNE0"))

    # Load the verified shared initialization point
    def initialize_ray_bootstrap(self, *, args, tunesetup, mc_steer):
        self.init_data(args.run_name, tunesetup.datacards)
        fingerprint = self.bootstrap_payload(runtime_sha256=args.runtime_sha256, datacards=tunesetup.datacards)["fingerprint"]
        _, points = icetune.load_ray_init_state(path=args.ray_init_state, fingerprint=fingerprint)
        return points

    # Identify the physics and inputs when Ray restarts
    def physics_fingerprint(self, *, param, tunesetup):
        return self.bootstrap_payload(runtime_sha256=param.get("runtime_sha256"), datacards=param["datacards"])["fingerprint"]

    # Start inside the declared domains, retaining valid source defaults
    def get_initial_param(self, param_space, aux_param_space, cdir, tune_default="TUNE0"):
        catalog = aux_param_space.get("pandora_catalog")
        if param_space and not catalog:
            raise ValueError('Pandora aux_param_space is missing required "pandora_catalog"')
        entries, initial = _entry_by_key(catalog or {}), {}
        for key, spec in param_space.items():
            bounds = [_safe_float(getattr(spec, edge, None)) for edge in ("lower", "upper")]
            entry = entries.get(key, {})
            value = float(entry.get("initial", entry.get("default", sum(bounds) / 2 if None not in bounds else 0.0)))
            if "tune_bounds" in entry:
                low, high = entry["tune_bounds"]
                value = min(high, max(low, value))
                if entry["tune_dtype"] == "int":
                    value = min(math.floor(high), max(math.ceil(low), round(value)))
            initial[key] = value
        if catalog:
            _validate_config(catalog, initial)
        return initial

    # Write a manual steering copy using the frozen runtime
    def create_steering_card(self, param_space, tunename, cdir=None, tune_default="TUNE0"):
        aux = param_space["_aux_param_space"]
        config = _trial_config(aux["pandora_catalog"], {k: v for k, v in param_space.items() if k != "_aux_param_space"}, self.datacards)
        stage = self._stage_trial(config=config, aux_param_space=aux, datacard=self.datacards[0],
            cdir=cdir or os.getcwd(), run_name="manual", tunename=tunename)
        return {key: str(stage[key]) for key in ("trial_dir", "settings_xml")}

    # Render and validate all parameter changes before publishing the approved rows
    def push_parameters(self, *, summary, target_path, cdir, options, confirm):
        if options:
            raise ValueError(f"Pandora push does not support driver options: {sorted(options)}")
        config, target = self._resolve_push_inputs(summary, target_path, cdir)
        catalog = (summary.get("card_config") or {}).get("catalog")
        if catalog:
            config = _validate_config(catalog, config)
        if target.suffix.lower() == ".json" and not target.is_dir():
            # Export fitted parameters as a later tune baseline, retaining untuned values
            previous = icetune_push.load_summary(target, baseline=True) if target.exists() else {}
            old = previous.get("config", {})
            payload = dict(source=summary.get("_push_input_path"), config={**old, **config})
            previous_catalog = (previous.get("card_config") or {}).get("catalog", {})
            if catalog or previous_catalog:
                merged = copy.deepcopy(catalog or previous_catalog)
                for section in ("xml_numeric", "xml_boolean", "wrapper_numeric"):
                    entries = {entry["key"]: entry for entry in previous_catalog.get(section, [])}
                    for entry in merged.get(section, []):
                        if entry["key"] in entries and entry.get("identity") != entries[entry["key"]].get("identity"):
                            raise ValueError(f"Pandora baseline parameter identity changed: {entry['key']}")
                        entries[entry["key"]] = entry
                    merged[section] = list(entries.values())
                merged["ordered_parameters"] = merged.get("ordered_parameters", []) + [
                    row for row in previous_catalog.get("ordered_parameters", []) if row not in merged.get("ordered_parameters", [])]
                payload["card_config"] = _card_config(_validate_config(merged, payload["config"]), merged)
            rendered = {target: json.dumps(payload, indent=4, allow_nan=False) + "\n"}
            old_values = {key: old.get(key) for key in config}
        else:
            if not target.exists():
                raise ValueError(f'Pandora push target "{target}" does not exist')
            xml_keys = {str(k) for k in config if str(k).startswith(("PXML_", "PBOOL_"))}
            settings = _resolve_push_xml(target, summary, xml_keys) if xml_keys else None
            steering = _resolve_push_python(target) if any(str(k).startswith("PWRAP_") for k in config) else None
            catalog = _push_catalog(summary, settings)
            _validate_config(catalog, config)
            rendered, old_values = {}, {}
            for path, render in ((settings, _render_xml_push), (steering, _render_python_push)):
                if path is not None:
                    rendered[path], values = render(path, catalog, config)
                    old_values.update(values)
            missing = set(config) - old_values.keys()
            if missing:
                raise ValueError(f"Pandora push did not resolve parameters {sorted(missing)}")
        rows = [(key, old_values[key], value) for key, value in config.items()]
        if not confirm(rows):
            raise icetune_push.PushCancelled
        for path in rendered:
            ensure_dir(path.parent)
        icetune_push.atomic_write_texts(rendered)
        return rows

    # Reject reconstruction through the integration grid interface
    def compute(self, *args, **kwargs):
        if kwargs.get("datacards"):
            raise RuntimeError("Pandora reconstruction must use evaluate_trial_outputs() with local staging")
        return {"bootstrap": "noop"}

    # Keep simulation inputs on shared storage when Ray uploads the reconstruction runtime
    def shared_runtime_files(self, param):
        return [item["path"] for card in param["datacards"] for item in card["inputs"]]

    # Stage verified reconstruction and simulation inputs once in worker scratch
    def _resolve_datacard_paths(self, datacard, cdir):
        cache = pathlib.Path(cdir) / "tmp/icetune/pandora/cache"
        root = runtime.stage_runtime(datacard["runtime"], cdir, cache)
        inputs = [runtime.stage_file(item["path"], cache / "inputs" / (short_id(item["sha256"]) + ".root"),
                                    {key: item[key] for key in ("sha256", "size")}) for item in datacard["inputs"]]
        return dict(input_roots=inputs, settings_template=root / datacard["settings_template"],
                    steering_template=root / datacard["reco_steering"], run_dir=root / datacard["run_dir"],
                    setup_script=root / "setup.sh", pandora_dir=root)

    # Apply trial parameters to independent copies of the frozen steering
    def _stage_trial(self, *, config, aux_param_space, datacard, cdir, run_name, tunename, dataset_index=None):
        paths = self._resolve_datacard_paths(datacard, cdir)
        catalog = aux_param_space.get("pandora_catalog")
        if not catalog:
            raise ValueError('Pandora aux_param_space is missing required "pandora_catalog"')
        root = _trial_dir(cdir=cdir, run_name=run_name, tunename=tunename)
        if dataset_index is not None:
            root /= f"dataset_{dataset_index:03d}_{_PARAM_SAFE_RE.sub('_', datacard['dataset_name'])}"
        ensure_dir(root)
        settings, steering, wrapper = (root / name for name in (
            "PandoraSettings.xml", "run_reco_pandora.py", "pandora_wrapper_overrides.json"))
        tree = ET.parse(paths["settings_template"])
        if _apply_config_to_xml(tree, catalog, config):
            tree.write(settings, encoding="utf-8", xml_declaration=True)
        else:
            shutil.copyfile(paths["settings_template"], settings)
        _write_json(wrapper, _wrapper_overrides_from_config(catalog, config))
        _copy_and_patch_steering(source_path=paths["steering_template"], target_path=steering,
                                override_path=wrapper, data_dir=paths["run_dir"])
        _stage_xml_inputs(tree, source_dir=paths["run_dir"], target_dir=root)
        return dict(paths, trial_dir=root, settings_xml=settings, steering_py=steering, wrapper_json=wrapper,
                    logs_dir=root / "logs", reco_roots=[root / f"reco_{i:03d}.root" for i in range(len(paths["input_roots"]))])

    # Reconstruct input files sequentially within the remaining trial time
    def _run_reconstruction(self, *, stage, datacard, max_t):
        deadline = time.monotonic() + max_t
        for index, (source, target) in enumerate(zip(stage["input_roots"], stage["reco_roots"], strict=True)):
            remaining = deadline - time.monotonic()
            if remaining <= 0.0:
                raise TimeoutError("Pandora reconstruction exceeded the trial time limit")
            command = ["k4run", stage["steering_py"], "--pandoraSettings", stage["settings_xml"],
                       "--IOSvc.Input", source, "--IOSvc.Output", target, "--nevents", datacard["nevents"]]
            _run_command(list(map(str, command)), cwd=stage["trial_dir"], log_path=stage["logs_dir"] / f"reco_{index:03d}.log",
                         max_t=remaining, setup_script=stage["setup_script"])


    # Reconstruct one dataset and retain its compact numerical summary
    def _evaluate_dataset(self, config, param, tunename, index, card, deadline):
        stage = self._stage_trial(config=config, aux_param_space=param["aux_param_space"], datacard=card,
            cdir=param["cdir"], run_name=param["run_name"], tunename=tunename, dataset_index=index)
        try:
            self._run_reconstruction(stage=stage, datacard=card, max_t=deadline - time.monotonic())
            arrays = _reco_chunks(stage["reco_roots"], card["objectives"])
            spec = dict(card["objectives"]["pflow"], _collect_plot_details=param.get("plot", True))
            metrics, detail = _objective_pflow(spec, arrays)
            details = {"pflow": detail}
        finally:
            _retire_reco_roots(stage["reco_roots"])
        dataset = dict(name=card["dataset_name"], metrics=metrics, trial_dir=str(stage["trial_dir"]),
                       **{key: card[key] for key in ("process", "sqrts_gev", "weight")})
        _write_json(stage["trial_dir"] / "objective_summary.json",
                    dict(dataset=dataset, details=_compact_objective_details(details)))
        return dict(datacard=card, metrics=metrics, details=details, dataset=dataset)

    # Summarize the weighted event union and retain finite penalties for failed trials
    def evaluate_trial_outputs(self, *, config, param, trial_id, tunename):
        outputs = dict(config=copy.deepcopy(config), trial_id=trial_id, tunename=tunename, error=None, results=None)
        try:
            effective = _trial_config(param["aux_param_space"]["pandora_catalog"], config, param["datacards"])
            deadline = time.monotonic() + param.get("max_t", 7200)
            records = [self._evaluate_dataset(effective, param, tunename, i, card, deadline)
                       for i, card in enumerate(param["datacards"]) if _dataset_weight(card) > 0.0]
            metrics, details = _combine_datasets(records)
            root = _trial_dir(cdir=param["cdir"], run_name=param["run_name"], tunename=tunename)
            summary = dict(schema_version=1, trial_id=trial_id, tunename=tunename, config=effective, metrics=metrics,
                           objectives=records[0]["datacard"]["objectives"], datasets=[r["dataset"] for r in records],
                           calibration=details["pflow"]["calibration"],
                           plots={}, plots_deferred_to_best_publish=param.get("plot", True))
            _write_json(root / "objective_summary.json", summary)
            outputs.update(config=effective, card_config=_card_config(effective, param["aux_param_space"]["pandora_catalog"]), metrics=metrics,
                likelihood=dict(schema_version=1, valid=True, kind="pandora_objective",
                                objective=dict(name="loss", value=metrics[param["cost"]], source_metric=param["cost"])),
                results=dict(trial_dir=str(root), fig_dir=str(root / "figs"), objective_summary=summary,
                             details=_compact_objective_details(details), plot_details=details if param.get("plot", True) else {},
                             dataset_details=[_compact_objective_details(r["details"]) for r in records]))
        except OSError:
            raise
        except Exception as exc:
            logger.exception("Pandora trial %s failed", trial_id)
            outputs.update(error=str(exc), metrics={"pflow": PENALTY_COST},
                           likelihood=dict(schema_version=1, valid=False, kind="pandora_penalty", error=str(exc)))
        return outputs

    # Publish selected trial diagnostics from their retained numerical inputs
    def render_trial_figures_to_dir(self, *, outputs, param, summary_payload, output_dir, summary_file=None):
        out, results = pathlib.Path(output_dir), outputs.get("results") or {}
        ensure_dir(out)
        paths = plots.write_plots(diagnostics.prepare_plots(results["plot_details"]), out)
        summary = copy.deepcopy(results.get("objective_summary", {}))
        summary.update(plots={key: pathlib.Path(path).relative_to(out).as_posix() for key, path in paths.items()},
                       plots_deferred_to_best_publish=False)
        payload = dict(copy.deepcopy(summary_payload), objective_summary=summary)
        if "calibration" in summary:
            payload["calibration"] = copy.deepcopy(summary["calibration"])
            if not summary["calibration"]["passed"]:
                logger.warning("No calibrated Pandora solution in this summary: mean Reco/Gen = %.4f, tolerance = %.4f",
                               summary["calibration"]["mean"], summary["calibration"]["tolerance"])
        _write_json(pathlib.Path(summary_file) if summary_file else out / "summary.json", payload)
        return payload

    # Retain interrupted reconstruction outputs with the requested old file suffix
    def cleanup_trial_outputs(self, tunename, datacards, mc_steer, cdir=None):
        for root in (pathlib.Path(cdir or os.getcwd()) / "runs/icetune").glob("*/trials/" + str(tunename)):
            _retire_reco_roots(root.rglob("reco_*.root"))
