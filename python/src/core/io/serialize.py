# Shared JSON-safe numerical payload conversions
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import base64
import copy
import hashlib
import json
import math
import zlib
from os import PathLike
from pathlib import Path

import numpy as np

from core.io import jsonref
from core.io.files import ensure_dir


# Decode exact big endian IEEE 754 values from one cache payload
def decode_doubles(payload: object, field: str) -> np.ndarray:
    if not isinstance(payload, str):
        raise ValueError(f'Compressed array field "{field}" is not encoded text')
    try:
        raw = zlib.decompress(base64.b64decode(payload, validate=True))
    except (ValueError, zlib.error) as exc:
        raise ValueError(f'Compressed array field "{field}" has invalid compression') from exc
    if len(raw) % 8 != 0:
        raise ValueError(f'Compressed array field "{field}" has an invalid byte count')
    values = np.frombuffer(raw, dtype=">f8").astype(np.float64)
    if not np.all(np.isfinite(values)):
        raise ValueError(f'Compressed array field "{field}" contains a nonfinite value')
    return values


# Compute the SHA-256 digest of one file
def sha256_file(path: str | PathLike[str]) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


# Load one JSON-compatible file through the selected parser
def load_json_file(path: str | PathLike[str], *, loader=json.load):
    filename, separator, fragment = str(path).partition("#")
    if separator:
        target, parts = jsonref.reference({"$ref": "#" + fragment}, Path(filename).resolve())
        return jsonref.JsonReader(loader).node(target, parts)
    return jsonref.load(path, loader=loader)


# Write one JSON payload with caller-selected encoder options
def write_json_file(
    path: str | PathLike[str], payload, *, newline: bool = False, reference_root: str | PathLike[str] | None = None, **options
) -> None:
    from pathlib import Path

    import pyjson5

    options.setdefault("default", json_default)
    if Path(path).is_file() and not jsonref.has_references(payload):
        reader = jsonref.JsonReader(pyjson5.load)
        try:
            raw = reader.document(path)
        except (ValueError, pyjson5.Json5Exception):
            raw = None
        if jsonref.has_references(raw):
            documents = jsonref.edits(path, payload, loader=pyjson5.load)
            if reference_root is not None:
                root = Path(reference_root).resolve()
                if any(not filename.is_relative_to(root) for filename in documents):
                    raise ValueError(f"Shared fit parameters must be stored inside the tune directory {root}")
            for filename, document in documents.items():
                with filename.open("w", encoding="utf-8") as target:
                    json.dump(document, target, **options)
                    if newline:
                        target.write("\n")
            return
    ensure_dir(Path(path).parent)
    with open(path, "w", encoding="utf-8") as target:
        json.dump(payload, target, **options)
        if newline:
            target.write("\n")


# Convert one numerical or path-like value for a standard JSON encoder
def json_default(value, *, nonfinite: str = "null"):
    if isinstance(value, PathLike):
        return str(value)
    if isinstance(value, np.floating):
        value = float(value)
        return value if math.isfinite(value) else str(value) if nonfinite == "string" else None
    if isinstance(value, np.complexfloating):
        return complex(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if hasattr(value, "detach"):
        return value.detach().cpu().numpy().tolist()
    if hasattr(value, "tolist"):
        return value.tolist()
    if hasattr(value, "item"):
        return value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return str(value) if nonfinite == "string" else None
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


# Convert nested numerical values into JSON-compatible containers
def json_safe(value, *, nonfinite: str = "null"):
    if isinstance(value, dict):
        return {str(key): json_safe(item, nonfinite=nonfinite) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item, nonfinite=nonfinite) for item in value]
    try:
        return json_safe(json_default(value, nonfinite=nonfinite), nonfinite=nonfinite)
    except TypeError:
        try:
            json.dumps(value)
            return value
        except TypeError:
            return repr(value)


# Compute one finite float or the caller-selected default
def finite_or_none(value, default=None):
    try:
        output = float(value)
    except (TypeError, ValueError):
        return default
    return output if math.isfinite(output) else default


# Copy selected payload fields with optional explicit missing values
def selected_fields(payload: dict, names, *, include_missing: bool = False) -> dict:
    return {str(name): copy.deepcopy(payload.get(name)) for name in names if include_missing or name in payload}
