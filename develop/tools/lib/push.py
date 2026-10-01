# Preserve steering card formatting while pushing derived tool parameters
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.


from __future__ import annotations

import copy
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pyjson5 as json5
from core.io import jsonref
from core.tune.push import atomic_write_texts, confirm_change_table, render_json5_arrays, render_json5_scalars
from jsonpointer import JsonPointer


@dataclass(frozen=True)
class ScalarUpdate:
    """Scalar card update"""

    path: tuple[Any, ...]
    label: str
    value: Any


# Compute the value stored at one fully resolved JSON5 path
def _path_value(payload: Any, path: tuple[Any, ...]) -> Any:
    return JsonPointer(jsonref.pointer(path)).resolve(payload)


# Assign one scalar value at a fully resolved JSON5 path
def _assign_path(payload: Any, path: tuple[Any, ...], value: Any) -> None:
    if not path:
        raise ValueError("A push update cannot replace the complete JSON5 document")
    target = JsonPointer(jsonref.pointer(path[:-1])).resolve(payload)
    key = int(path[-1]) if isinstance(target, list) else path[-1]
    old = target[key]
    if isinstance(old, (dict, list)) or isinstance(value, (dict, list)):
        raise ValueError(f"A push update must target one scalar at {path}")
    target[key] = value


# Prepare comment preserving text and preview rows for all target cards
def _prepare_updates(
    updates_by_card: dict[Path, list[ScalarUpdate]],
) -> tuple[dict[Path, str], list[tuple[str, Any, Any]], dict[Path, str]]:
    reader = jsonref.JsonReader(json5.load)
    desired: dict[Path, Any] = {}
    labels: dict[tuple[Path, str], str] = {}
    assigned: dict[tuple[Path, str], Any] = {}
    for card, updates in updates_by_card.items():
        for update in updates:
            target, parts = reader.origin(card, update.path)
            key = target, jsonref.pointer(parts)
            if key in assigned:
                if assigned[key] != update.value:
                    raise ValueError(f"Conflicting push values for shared parameter {target}#{key[1]}")
                continue
            assigned[key] = update.value
            labels[key] = update.label
            document = desired.setdefault(target, copy.deepcopy(reader.document(target)))
            _assign_path(document, tuple(parts), update.value)

    rendered: dict[Path, str] = {}
    sources: dict[Path, str] = {}
    rows: list[tuple[str, Any, Any]] = []
    for card, document in desired.items():
        source = card.read_text(encoding="utf-8")
        sources[card] = source
        output, changed = render_json5_scalars(source, document)
        if not changed:
            continue
        rendered[card] = output
        rows.extend(
            (labels[card, jsonref.pointer(path)], _path_value(reader.document(card), path), _path_value(document, path))
            for path in changed
        )
    return rendered, rows, sources


# Confirm tool changes with enough precision for small photoproduction couplings
def confirm_tool_change_table(rows: list[tuple[str, Any, Any]]) -> bool:
    return confirm_change_table(rows, precision=9)


# Preview, confirm and atomically apply scalar changes to all target cards
def push_json5_updates(
    updates_by_card: dict[Path, list[ScalarUpdate]],
    *,
    confirm: Callable[[list[tuple[str, Any, Any]]], bool] = confirm_tool_change_table,
    array_updates: dict[Path, dict[tuple[Any, ...], list]] | None = None,
) -> bool:
    rendered, rows, sources = _prepare_updates(updates_by_card)
    for card, arrays in (array_updates or {}).items():
        card = card.resolve()
        source = sources.setdefault(card, card.read_text(encoding="utf-8"))
        output, changed = render_json5_arrays(rendered.get(card, source), arrays)
        if changed:
            rendered[card] = output
        original = json5.loads(source)
        for path in changed:
            old = JsonPointer(jsonref.pointer(path)).resolve(original, default=None)
            rows.append((f"{card.name}:{'.'.join(map(str, path))}", old, arrays[path]))
    if not rows:
        print("No parameter changes to push")
        return False
    if not confirm(rows):
        print("Push cancelled")
        return False
    for card, source in sources.items():
        if card.read_text(encoding="utf-8") != source:
            raise ValueError(f"Card changed after the push preview: {card}")
    atomic_write_texts(rendered)
    print(f"Push completed: updated {len(rows)} values in {len(rendered)} card files")
    return True
