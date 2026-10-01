# Publish fitted parameters into steering cards without changing their formatting
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

"""Shared helpers for publishing fitted parameters into steering files."""

from __future__ import annotations

import copy
import json
import math
import os
import pathlib
import shutil
import uuid
from dataclasses import dataclass
from numbers import Real
from typing import Any

import pyjson5 as json5
from tabulate import tabulate

from core.io.serialize import load_json_file
from core.tune import summary as fit_summary


class PushCancelled(Exception):
    """Signal that the user declined a prepared parameter push."""


@dataclass(frozen=True)
class _Token:
    """Store one significant JSON5 token and its source span."""

    text: str
    start: int
    end: int


# Parse the driver-specific push options from one strict JSON object
def parse_driver_options(value: str) -> dict:
    try:
        options = json.loads(value)
    except (json.JSONDecodeError, TypeError) as exc:
        raise ValueError("Push driver options must be a valid JSON object") from exc
    if not isinstance(options, dict):
        raise ValueError("Push driver options must be a JSON object")
    return options


# Compute significant JSON5 tokens without changing their source positions
def _json5_tokens(text: str) -> list[_Token]:
    tokens: list[_Token] = []
    index = 0
    while index < len(text):
        char = text[index]
        if char.isspace():
            index += 1
            continue
        if text.startswith("//", index):
            newline = text.find("\n", index + 2)
            index = len(text) if newline < 0 else newline + 1
            continue
        if text.startswith("/*", index):
            end = text.find("*/", index + 2)
            if end < 0:
                raise ValueError("Unterminated JSON5 block comment")
            index = end + 2
            continue
        if char in "{}[]:,":
            tokens.append(_Token(char, index, index + 1))
            index += 1
            continue
        if char in "\"'":
            start = index
            quote = char
            index += 1
            while index < len(text):
                if text[index] == "\\":
                    index += 2
                    continue
                if text[index] == quote:
                    index += 1
                    break
                index += 1
            else:
                raise ValueError("Unterminated JSON5 string")
            tokens.append(_Token(text[start:index], start, index))
            continue
        start = index
        while index < len(text):
            if text[index].isspace() or text[index] in "{}[]:,":
                break
            if text.startswith("//", index) or text.startswith("/*", index):
                break
            index += 1
        if start == index:
            raise ValueError(f"Unsupported JSON5 token at byte {index}")
        tokens.append(_Token(text[start:index], start, index))
    return tokens


class _Json5SpanParser:
    """Map every JSON5 scalar path to its exact source span."""

    # Initialize one parser over significant source tokens
    def __init__(self, text: str):
        self.tokens = _json5_tokens(text)
        self.index = 0
        self.spans: dict[tuple[Any, ...], tuple[int, int]] = {}
        self.objects: dict[tuple[Any, ...], tuple[int, int]] = {}
        self.arrays: dict[tuple[Any, ...], tuple[int, int]] = {}

    # Consume and return one token, optionally checking its text
    def _take(self, expected: str | None = None) -> _Token:
        if self.index >= len(self.tokens):
            raise ValueError("Unexpected end of JSON5 input")
        token = self.tokens[self.index]
        if expected is not None and token.text != expected:
            raise ValueError(f'Expected JSON5 token "{expected}", got "{token.text}"')
        self.index += 1
        return token

    # Decode one quoted or unquoted JSON5 object key
    def _key(self) -> str:
        token = self._take()
        if token.text[:1] in {'"', "'"}:
            value = json5.loads(token.text)
            if not isinstance(value, str):
                raise ValueError("JSON5 object key is not a string")
            return value
        return token.text

    # Parse one JSON5 value and record all nested scalar spans
    def _value(self, path: tuple[Any, ...]) -> None:
        token = self.tokens[self.index]
        if token.text == "{":
            self._object(path)
            return
        if token.text == "[":
            self._array(path)
            self.arrays[path] = (token.start, self.tokens[self.index - 1].end)
            return
        scalar = self._take()
        self.spans[path] = (scalar.start, scalar.end)

    # Parse one JSON5 object
    def _object(self, path: tuple[Any, ...]) -> None:
        opening = self._take("{")
        if self.tokens[self.index].text == "}":
            closing = self._take("}")
            self.objects[path] = (opening.start, closing.start)
            return
        while True:
            key = self._key()
            self._take(":")
            self._value((*path, key))
            token = self._take()
            if token.text == "}":
                self.objects[path] = (opening.start, token.start)
                return
            if token.text != ",":
                raise ValueError(f'Expected JSON5 token "," or "}}", got "{token.text}"')
            if self.tokens[self.index].text == "}":
                closing = self._take("}")
                self.objects[path] = (opening.start, closing.start)
                return

    # Parse one JSON5 array
    def _array(self, path: tuple[Any, ...]) -> None:
        self._take("[")
        if self.tokens[self.index].text == "]":
            self._take("]")
            return
        item = 0
        while True:
            self._value((*path, item))
            item += 1
            token = self._take()
            if token.text == "]":
                return
            if token.text != ",":
                raise ValueError(f'Expected JSON5 token "," or "]", got "{token.text}"')
            if self.tokens[self.index].text == "]":
                self._take("]")
                return

    # Parse the complete document and return scalar source spans
    def parse(self) -> dict[tuple[Any, ...], tuple[int, int]]:
        self._value(())
        if self.index != len(self.tokens):
            raise ValueError("Trailing JSON5 tokens")
        return self.spans


# Compute true when two scalar values have the same steering meaning
def _same_scalar(old: Any, new: Any) -> bool:
    if isinstance(old, Real) and not isinstance(old, bool) and isinstance(new, Real) and not isinstance(new, bool):
        return float(old) == float(new)
    return type(old) is type(new) and old == new


# Collect changed scalar paths while rejecting structural card changes
def _changed_scalar_paths(old: Any, new: Any, path: tuple[Any, ...] = ()) -> list[tuple[Any, ...]]:
    if isinstance(old, dict) or isinstance(new, dict):
        if not isinstance(old, dict) or not isinstance(new, dict) or old.keys() != new.keys():
            raise ValueError(f"Push would change JSON5 object structure at {path}")
        return [changed for key in old for changed in _changed_scalar_paths(old[key], new[key], (*path, key))]
    if isinstance(old, list) or isinstance(new, list):
        if not isinstance(old, list) or not isinstance(new, list) or len(old) != len(new):
            raise ValueError(f"Push would change JSON5 array structure at {path}")
        return [
            changed
            for index, (old_item, new_item) in enumerate(zip(old, new, strict=True))
            for changed in _changed_scalar_paths(old_item, new_item, (*path, index))
        ]
    return [] if _same_scalar(old, new) else [path]


# Encode one replacement scalar as strict JSON accepted by JSON5 readers
def _scalar_text(value: Any) -> str:
    if isinstance(value, float) and not math.isfinite(value):
        raise ValueError("Cannot push a non-finite parameter value")
    return json.dumps(value, ensure_ascii=True, allow_nan=False, separators=(",", ":"))


# Collect explicitly allowed object additions while rejecting other structural changes
def _object_additions(
    old: Any, new: Any, allowed: set[tuple[Any, ...]], path: tuple[Any, ...] = ()
) -> list[tuple[tuple[Any, ...], str, dict]]:
    if isinstance(old, dict) or isinstance(new, dict):
        if not isinstance(old, dict) or not isinstance(new, dict):
            raise ValueError(f"Push would change JSON5 object structure at {path}")
        removed = old.keys() - new.keys()
        if removed:
            raise ValueError(f"Push would change JSON5 object structure at {path}")
        additions = []
        for key in new:
            if key in old:
                continue
            child = (*path, key)
            if child not in allowed or not isinstance(new[key], dict):
                raise ValueError(f"Push would change JSON5 object structure at {path}")
            additions.append((path, key, new[key]))
        for key in old:
            if key in new:
                additions.extend(_object_additions(old[key], new[key], allowed, (*path, key)))
        return additions
    if isinstance(old, list) or isinstance(new, list):
        if not isinstance(old, list) or not isinstance(new, list) or len(old) != len(new):
            raise ValueError(f"Push would change JSON5 array structure at {path}")
        return [
            addition
            for index, (old_item, new_item) in enumerate(zip(old, new, strict=True))
            for addition in _object_additions(old_item, new_item, allowed, (*path, index))
        ]
    return []


# Add one property to a multiline JSON5 object without changing surrounding text
def _insert_object(source: str, parent: tuple[Any, ...], key: str, value: Any) -> str:
    parser = _Json5SpanParser(source)
    parser.parse()
    if parent not in parser.objects:
        raise ValueError(f"Could not locate JSON5 object path {parent}")
    opening, closing = parser.objects[parent]
    line_start = source.rfind("\n", 0, closing) + 1
    indent = source[line_start:closing]
    if indent.strip():
        raise ValueError(f"Cannot add a JSON5 property to an inline object at {parent}")

    previous = next((token for token in reversed(parser.tokens) if opening <= token.start < closing), None)
    if previous is None:
        raise ValueError(f"Could not locate JSON5 object contents at {parent}")
    encoded = json.dumps(value, ensure_ascii=True, allow_nan=False, separators=(", ", ": "))
    if isinstance(value, list) and value and isinstance(value[0], list):
        encoded = "\n" + indent + "      [" + (",\n" + indent + "       ").join(
            json.dumps(row, allow_nan=False, separators=(", ", ": ")) for row in value) + "]"
    entry = f"{indent}  {json.dumps(key)}:" + ("" if encoded.startswith("\n") else " ") + encoded + "\n"
    edits = [(line_start, entry)]
    if previous.text not in {"{", ","}:
        edits.append((previous.end, ","))
    for position, text in sorted(edits, reverse=True):
        source = source[:position] + text + source[position:]
    return source


# Render desired JSON5 scalar values into the original source text
def render_json5_scalars(
    source: str, desired: Any, allowed_object_additions: set[tuple[Any, ...]] | None = None
) -> tuple[str, list[tuple[Any, ...]]]:
    old = json5.loads(source)
    allowed = set() if allowed_object_additions is None else allowed_object_additions
    additions = _object_additions(old, desired, allowed)
    rendered = source
    for parent, key, value in additions:
        rendered = _insert_object(rendered, parent, key, value)
    old = json5.loads(rendered)
    changed = _changed_scalar_paths(old, desired)
    spans = _Json5SpanParser(rendered).parse()
    missing = [path for path in changed if path not in spans]
    if missing:
        raise ValueError(f"Could not locate JSON5 scalar paths: {missing}")
    for path in sorted(changed, key=lambda item: spans[item][0], reverse=True):
        start, end = spans[path]
        value = desired
        for component in path:
            value = value[component]
        rendered = rendered[:start] + _scalar_text(value) + rendered[end:]
    if json5.loads(rendered) != desired:
        raise ValueError("JSON5 push validation failed")
    return rendered, [(*parent, key) for parent, key, _ in additions] + changed


# Render explicitly selected array parameters while preserving all surrounding JSON5 text
def render_json5_arrays(source: str, updates: dict[tuple[Any, ...], list]) -> tuple[str, list[tuple[Any, ...]]]:
    desired = copy.deepcopy(json5.loads(source))
    changed = []
    for path, value in updates.items():
        if not path or not isinstance(value, list):
            raise ValueError("Array updates require a property path and an array")
        parent = desired
        for key in path[:-1]:
            parent = parent[key]
        if not isinstance(parent, dict) or (path[-1] in parent and not isinstance(parent[path[-1]], list)):
            raise ValueError(f"Expected an array property at {path}")
        if parent.get(path[-1]) == value:
            continue
        parser = _Json5SpanParser(source)
        parser.parse()
        if path in parser.arrays:
            start, end = parser.arrays[path]
            indent = " " * (start - source.rfind("\n", 0, start) - 1)
            encoded = json.dumps(value, allow_nan=False, separators=(", ", ": "))
            if value and isinstance(value[0], list):
                encoded = "[" + (",\n" + indent + " ").join(
                    json.dumps(row, allow_nan=False, separators=(", ", ": ")) for row in value) + "]"
            source = source[:start] + encoded + source[end:]
        else:
            source = _insert_object(source, path[:-1], path[-1], value)
        parent[path[-1]] = value
        changed.append(path)
    if json5.loads(source) != desired:
        raise ValueError("JSON5 array update validation failed")
    return source, changed


# Replace all rendered text files or restore their exact prior state
def atomic_write_texts(rendered: dict[pathlib.Path, str]) -> None:
    staged: dict[pathlib.Path, pathlib.Path] = {}
    backups: dict[pathlib.Path, pathlib.Path] = {}
    try:
        for path, text in rendered.items():
            token = f"{os.getpid()}.{uuid.uuid4().hex}"
            stage = path.with_name(f".{path.name}.stage.{token}")
            stage.write_text(text, encoding="utf-8")
            if path.exists():
                os.chmod(stage, path.stat().st_mode)
                backup = path.with_name(f".{path.name}.backup.{token}")
                shutil.copy2(path, backup)
                backups[path] = backup
            staged[path] = stage
        try:
            for path, stage in staged.items():
                os.replace(stage, path)
        except BaseException as commit_error:
            rollback_errors = []
            for path in rendered:
                try:
                    if path in backups:
                        os.replace(backups[path], path)
                    else:
                        path.unlink(missing_ok=True)
                except BaseException as exc:
                    rollback_errors.append(f"{path}: {exc}")
            if rollback_errors:
                backups.clear()
                raise RuntimeError(
                    "Parameter card rollback failed and backups were retained: " + ", ".join(rollback_errors)
                ) from commit_error
            raise
    finally:
        for sidecar in (*staged.values(), *backups.values()):
            sidecar.unlink(missing_ok=True)


# Convert one flat best-fit mapping into a finite numeric parameter config
def _normalize_parameter_config(
    parameters: dict, source: str, value_records: bool, *, allow_empty: bool = False
) -> dict[str, float]:
    if not parameters and not allow_empty:
        raise ValueError(f'Push input field "{source}" is empty')
    config: dict[str, float] = {}
    for name, entry in parameters.items():
        if not str(name):
            raise ValueError(f'Push input field "{source}" contains an empty parameter name')
        value = entry
        if value_records:
            if not isinstance(entry, dict) or "value" not in entry:
                raise ValueError(f'Push input parameter "{name}" in "{source}" has no numeric value')
            value = entry["value"]
        if isinstance(value, bool) or not isinstance(value, Real):
            raise ValueError(f'Push input parameter "{name}" in "{source}" is not numeric')
        numeric = float(value)
        if not math.isfinite(numeric):
            raise ValueError(f'Push input parameter "{name}" in "{source}" is not finite')
        config[str(name)] = numeric
    return config


# Normalize native icetune, icescape and iceproxy best-fit schemas
def normalize_summary(payload: dict, summary_path: pathlib.Path, *, baseline: bool = False) -> dict:
    if not isinstance(payload, dict):
        raise ValueError(f'Push input "{summary_path}" is not a JSON object')

    if "best_fit" in payload:
        try:
            parameters = fit_summary.best_fit_values(payload["best_fit"])
        except ValueError as exc:
            raise ValueError(f'Push input "{summary_path}" has invalid best_fit: {exc}') from exc
        source = "best_fit.parameters"
        value_records = False
    else:
        for source in ("best_config", "config", "surrogate_best", "optimized_best.theta",
                       "objective_summary.config", "best_result.config"):
            keys = source.split(".")
            if keys[0] not in payload:
                continue
            parameters = payload
            for key in keys:
                parameters = parameters.get(key) if isinstance(parameters, dict) else None
            value_records = source == "surrogate_best"
            break
        else:
            if not baseline:
                raise ValueError(f'Push input "{summary_path}" has no supported fitted parameter map')
            parameters, source, value_records = payload, "parameters", False
    if not isinstance(parameters, dict):
        raise ValueError(f'Push input field "{source}" in "{summary_path}" is not a parameter map')

    normalized = dict(payload)
    normalized["config"] = _normalize_parameter_config(parameters, source, value_records, allow_empty=baseline)
    normalized["_push_parameter_source"] = source
    normalized["_push_input_path"] = str(summary_path.resolve())
    return normalized


# Load one fit summary or, for baselines, a possibly empty parameter map
def load_summary(path: str | os.PathLike[str], *, baseline: bool = False) -> dict:
    summary_path = pathlib.Path(path).expanduser()
    with summary_path.open("r", encoding="utf-8") as handle:
        payload = load_json_file(handle.name)
    return normalize_summary(payload, summary_path, baseline=baseline)


# Flatten numeric table leaves under one physical card parameter name
def _flatten_numeric(value: Any, prefix: str, output: dict[str, Any]) -> None:
    if isinstance(value, dict):
        for key, item in value.items():
            _flatten_numeric(item, f"{prefix}:{key}", output)
        return
    if isinstance(value, list):
        for index, item in enumerate(value):
            _flatten_numeric(item, f"{prefix}[{index}]", output)
        return
    if isinstance(value, Real):
        output[prefix] = value


# Flatten the physical parameters represented by one driver card config
def flatten_card_config(card_config: dict) -> dict[str, Any]:
    output = dict(card_config.get("parameters") or {})
    for base, table in (card_config.get("tables") or {}).items():
        for field in ("rows", "rho_mag", "rho_phase"):
            if field in table:
                _flatten_numeric(table[field], f"{base}:{field}", output)
    return output


# Compute true when two displayed parameter values agree numerically
def _same_display_value(old: Any, new: Any) -> bool:
    if isinstance(old, Real) and not isinstance(old, bool) and isinstance(new, Real) and not isinstance(new, bool):
        return math.isclose(float(old), float(new), rel_tol=1.0e-12, abs_tol=1.0e-12)
    return type(old) is type(new) and old == new


# Compute true when named parameters already describe one decoded card table
def _table_has_named_parameters(base: str, parameters: dict) -> bool:
    return any(key == base or key.startswith((f"{base}(", f"{base}[", f"{base}@")) for key in parameters)


# Compute display rows for changed physical fitted parameters
def card_config_rows(old: dict, new: dict) -> list[tuple[str, Any, Any]]:
    old_parameters = dict(old.get("parameters") or {})
    new_parameters = dict(new.get("parameters") or {})
    if old_parameters.keys() != new_parameters.keys():
        raise ValueError("Target and summary physical parameter sets differ")
    rows = [
        (key, old_parameters[key], new_parameters[key])
        for key in new_parameters
        if not _same_display_value(old_parameters[key], new_parameters[key])
    ]

    old_tables = old.get("tables") or {}
    new_tables = new.get("tables") or {}
    if old_tables.keys() != new_tables.keys():
        raise ValueError("Target and summary physical table sets differ")
    for base, new_table in new_tables.items():
        if _table_has_named_parameters(base, new_parameters):
            continue
        old_values: dict[str, Any] = {}
        new_values: dict[str, Any] = {}
        for field in ("rows", "rho_mag", "rho_phase"):
            if field in new_table:
                _flatten_numeric(new_table[field], f"{base}:{field}", new_values)
            if field in old_tables[base]:
                _flatten_numeric(old_tables[base][field], f"{base}:{field}", old_values)
        if old_values.keys() != new_values.keys():
            raise ValueError("Target and summary physical table values differ")
        rows.extend(
            (key, old_values[key], new_values[key])
            for key in new_values
            if not _same_display_value(old_values[key], new_values[key])
        )
    return rows


# Format one parameter value for the push table
def _display_value(value: Any, precision: int) -> str:
    if isinstance(value, Real) and not isinstance(value, bool):
        return f"{float(value):.{precision}f}"
    return str(value)


# Format the new-to-old ratio while handling zero and nonnumeric old values
def _display_ratio(old: Any, new: Any, precision: int) -> str:
    numeric = (
        isinstance(old, Real) and not isinstance(old, bool) and isinstance(new, Real) and not isinstance(new, bool)
    )
    if not numeric or float(old) == 0.0:
        return "n/a"
    return f"{float(new) / float(old):.{precision}f}"


# Compute an aligned ASCII table of old, new and ratio parameter values
def format_change_table(rows: list[tuple[str, Any, Any]], precision: int = 3) -> str:
    return tabulate(
        [
            (
                str(key),
                _display_value(old, precision),
                _display_value(new, precision),
                _display_ratio(old, new, precision),
            )
            for key, old, new in rows
        ],
        headers=("parameter", "old value", "new value", "new/old"),
        tablefmt="presto",
        colalign=("left", "right", "right", "right"),
        disable_numparse=True,
    )


# Print the prepared parameter changes and request explicit user approval
def confirm_change_table(rows: list[tuple[str, Any, Any]], precision: int = 3) -> bool:
    print(format_change_table(rows, precision=precision))
    while True:
        try:
            response = input("Apply these parameter changes? [yes,no]: ").strip().lower()
        except EOFError:
            print()
            return False
        if response in {"yes", "y"}:
            return True
        if response in {"no", "n"}:
            return False
        print('Please answer "yes" or "no".')
