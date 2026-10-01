# Generic JSON reference loading and edits preserving shared parameters
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import ast
import copy
import json
import re
from pathlib import Path
from urllib.parse import unquote

from jsonpointer import JsonPointer, JsonPointerException


# Encode arbitrary object keys and array indices as an RFC 6901 pointer
def pointer(parts):
    return "".join("/" + str(part).replace("~", "~0").replace("/", "~1") for part in parts)


# Parse object keys and array indices, including quoted process names
def field_parts(field):
    tokens = []
    pattern = re.compile(r"""([^.[\]@]+)|\[("(?:\\.|[^"\\])*"|'(?:\\.|[^'\\])*')\]|\[(\d+(?:\s*,\s*\d+)*)\]""")
    pos = 0
    while pos < len(field):
        match = pattern.match(field, pos)
        if match is None:
            raise ValueError(f"Invalid card field {field}")
        key, quoted, indices = match.groups()
        if indices is not None:
            tokens.extend(int(index) for index in indices.split(","))
        else:
            tokens.append(ast.literal_eval(quoted) if quoted is not None else key)
        pos = match.end()
        if pos < len(field) and field[pos] == ".":
            pos += 1
            if pos == len(field):
                raise ValueError(f"Invalid card field {field}")
        elif pos < len(field) and field[pos] != "[":
            raise ValueError(f"Invalid card field {field}")
    if not tokens:
        raise ValueError(f"Invalid card field {field}")
    return tokens


# Decode one local URI reference while retaining the containing file as its base
def reference(node, path):
    if set(node) != {"$ref"} or not isinstance(node["$ref"], str):
        raise ValueError("A reference must contain only a string $ref")
    value = node["$ref"]
    if re.search(r"%(?![0-9a-fA-F]{2})", value):
        raise ValueError("Invalid URI percent escape")
    filename, _, fragment = value.partition("#")
    filename, fragment = (unquote(part, encoding="utf-8", errors="strict") for part in (filename, fragment))
    if any(char in filename for char in (":", "?", "\0")) or filename.startswith("//"):
        raise ValueError("$ref requires a local file path and optional JSON Pointer fragment")
    target = (path.parent / filename).resolve() if filename else path
    return target, JsonPointer(fragment).get_parts()


# Hold source documents and resolved values for one independent read or edit
class JsonReader:
    # Select the parser used for every referenced document
    def __init__(self, loader=json.load):
        self.loader = loader
        self.documents = {}
        self.values = {}
        self.active = set()

    # Resolve filenames once for the source document cache
    def path(self, path):
        path = Path(path)
        return path if path in self.documents else path.resolve()

    # Read each source file once without expanding references
    def document(self, path):
        path = self.path(path)
        if path not in self.documents:
            with path.open(encoding="utf-8") as source:
                self.documents[path] = self.loader(source)
        return self.documents[path]

    # Find the source of a value through references at any ancestor or at the value
    def origin(self, path, parts):
        path, parts = self.path(path), [str(part) for part in parts]
        seen = set()
        references = set()
        while True:
            key = path, tuple(parts)
            if key in seen:
                raise ValueError(f"Circular JSON reference at {path}#{pointer(parts)}")
            seen.add(key)
            node = self.document(path)
            for index in range(len(parts) + 1):
                if isinstance(node, dict) and "$ref" in node:
                    location = path, tuple(parts[:index])
                    if location in references:
                        raise ValueError(f"Circular JSON reference at {path}#{pointer(parts[:index])}")
                    references.add(location)
                    path, prefix = reference(node, path)
                    parts = prefix + parts[index:]
                    break
                if index == len(parts):
                    return path, parts
                if isinstance(node, list) and not re.fullmatch(r"0|[1-9][0-9]*", parts[index]):
                    raise ValueError(f"Invalid JSON array index {parts[index]}")
                node = JsonPointer(pointer([parts[index]])).resolve(node)

    # Resolve one subtree with cycle detection and independent copied values
    def node(self, path, parts):
        path, parts = self.path(path), tuple(str(part) for part in parts)
        key = path, parts
        if key in self.values:
            return copy.deepcopy(self.values[key])
        if key in self.active:
            raise ValueError(f"Circular JSON reference at {path}#{pointer(parts)}")
        self.active.add(key)
        try:
            target, tokens = self.origin(path, parts)
            if (target, tuple(tokens)) != key:
                value = self.node(target, tokens)
            else:
                value = JsonPointer(pointer(parts)).resolve(self.document(path))
                value = self.expand(value, path, parts)
            self.values[key] = copy.deepcopy(value)
            return value
        except (OSError, ValueError, TypeError, KeyError, JsonPointerException) as error:
            raise ValueError(f"{path}#{pointer(parts)}: {error}") from error
        finally:
            self.active.remove(key)

    # Traverse ordinary values directly and resolve only reference objects
    def expand(self, value, path, parts):
        if isinstance(value, dict):
            if "$ref" in value:
                return self.node(path, parts)
            return {name: self.expand(item, path, (*parts, name)) for name, item in value.items()}
        if isinstance(value, list):
            return [self.expand(item, path, (*parts, index)) for index, item in enumerate(value)]
        return value

    # Load a complete JSON or JSON5 document through the common resolver
    def read(self, path):
        return self.node(path, ())


# Load arbitrary JSON values with local and relative file references
def load(path, *, loader=json.load):
    return JsonReader(loader).read(path)


# Resolve an in-memory document using its filename for relative references
def resolve(document, path, *, loader=json.load):
    reader = JsonReader(loader)
    reader.documents[Path(path).resolve()] = copy.deepcopy(document)
    return reader.read(path)


# Detect reference objects before deciding whether a write needs source tracking
def has_references(value):
    if isinstance(value, dict):
        return "$ref" in value or any(has_references(item) for item in value.values())
    return isinstance(value, list) and any(has_references(item) for item in value)


# Compute structural edits between a resolved input and its modified value
def _changes(old, new, parts=()):
    if isinstance(old, dict) and isinstance(new, dict):
        for key in old.keys() - new.keys():
            yield (*parts, key), True, None
        for key, value in new.items():
            if key not in old:
                yield (*parts, key), False, value
            else:
                yield from _changes(old[key], value, (*parts, key))
    elif isinstance(old, list) and isinstance(new, list) and len(old) == len(new):
        for index, (first, second) in enumerate(zip(old, new, strict=True)):
            yield from _changes(first, second, (*parts, index))
    elif old != new:
        yield parts, False, new


# Compute source-card edits while preserving aliases and rejecting conflicting updates
def edits(path, payload, *, loader=json.load):
    reader = JsonReader(loader)
    original = reader.read(path)
    changes = {}
    for parts, remove, value in _changes(original, payload):
        try:
            if remove:
                target, parent = reader.origin(path, parts[:-1])
                tokens = [*parent, str(parts[-1])]
            else:
                target, tokens = reader.origin(path, parts)
        except (KeyError, IndexError, ValueError, JsonPointerException):
            if remove or not parts:
                raise
            target, parent = reader.origin(path, parts[:-1])
            tokens = [*parent, str(parts[-1])]
        key = target, tuple(tokens)
        update = remove, value
        if key in changes and changes[key] != update:
            raise ValueError(f"Conflicting edits of shared JSON value {target}#{pointer(tokens)}")
        changes[key] = update
    documents = {}
    for (target, parts), (remove, value) in changes.items():
        document = documents.setdefault(target, copy.deepcopy(reader.document(target)))
        if not parts:
            documents[target] = copy.deepcopy(value)
            continue
        parent = JsonPointer(pointer(parts[:-1])).resolve(document)
        key = int(parts[-1]) if isinstance(parent, list) else parts[-1]
        if remove:
            del parent[key]
        else:
            parent[key] = copy.deepcopy(value)
    reader.documents.update(documents)
    reader.values.clear()
    reader.read(path)
    for target in documents:
        reader.read(target)
    return documents
