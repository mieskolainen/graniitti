# Shared construction and validation helpers for command-line interfaces
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>

from __future__ import annotations

import argparse
import os


# Parse an explicit true or false command-line value
def parse_bool(value) -> bool:
    if isinstance(value, bool):
        return value
    normalized = str(value).strip().lower()
    if normalized == "true":
        return True
    if normalized == "false":
        return False
    raise argparse.ArgumentTypeError(f"expected true or false, got {value!r}")


# Add one typed value option using its default to infer the value type
def add_value(parser, name: str, default, *, value_type=None, **options):
    option = name if name.startswith("-") else f"--{name}"
    inferred_type = value_type or (str if default is None else type(default))
    return parser.add_argument(option, type=inferred_type, default=default, **options)


# Add one disabled-by-default Boolean flag
def add_flag(parser, name: str, **options):
    option = name if name.startswith("-") else f"--{name}"
    return parser.add_argument(option, action="store_true", **options)


# Add one Boolean option with matching --no-* spelling
def add_toggle(parser, name: str, default: bool | None = True, **options):
    option = name if name.startswith("-") else f"--{name}"
    return parser.add_argument(option, action=argparse.BooleanOptionalAction, default=default, **options)


# Raise the first invalid configuration rule
def validate(rules: list[tuple[bool, str]], error_type=ValueError) -> None:
    for valid, message in rules:
        if not valid:
            raise error_type(message)


# Expand one optional path and return its absolute spelling
def absolute_path(value: str | None) -> str | None:
    return None if value is None else os.path.abspath(os.path.expanduser(value))
