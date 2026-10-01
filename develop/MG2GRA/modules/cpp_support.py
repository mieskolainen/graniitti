#!/usr/bin/env python3
#
# Shared C++ source transformations for MadGraph conversion
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import re

STD_NAMES = {
    "abs",
    "array",
    "cerr",
    "complex",
    "conj",
    "cout",
    "endl",
    "ifstream",
    "ios",
    "map",
    "max",
    "min",
    "ofstream",
    "real",
    "set",
    "setiosflags",
    "setw",
    "sqrt",
    "string",
    "transform",
    "vector",
}


# Compute one complete C++ function and its source offsets
def function_block(text: str, signature: str) -> tuple[int, int, str]:
    start = cpp_code_only(text).find(signature)
    if start < 0:
        raise RuntimeError(f"Function not found: {signature}")
    opening: int | None = None
    depth = 0
    quote: str | None = None
    in_block_comment = False
    in_line_comment = False
    pos = start
    while pos < len(text):
        char = text[pos]
        if in_line_comment:
            if char == "\n":
                in_line_comment = False
            pos += 1
            continue
        if in_block_comment:
            if text.startswith("*/", pos):
                in_block_comment = False
                pos += 2
            else:
                pos += 1
            continue
        if quote is not None:
            if char == "\\" and pos + 1 < len(text):
                pos += 2
                continue
            if char == quote:
                quote = None
            pos += 1
            continue
        if text.startswith("//", pos):
            in_line_comment = True
            pos += 2
            continue
        if text.startswith("/*", pos):
            in_block_comment = True
            pos += 2
            continue
        if char in ('"', "'"):
            quote = char
            pos += 1
            continue
        if char == "{":
            if opening is None:
                opening = pos
            depth += 1
        elif char == "}" and opening is not None:
            depth -= 1
            if depth == 0:
                return start, pos + 1, text[start : pos + 1]
        pos += 1
    if opening is None:
        raise RuntimeError(f"Opening brace not found: {signature}")
    raise RuntimeError(f"Unterminated function: {signature}")


# Canonicalize generated text without retaining trailing whitespace
def normalize_whitespace(text: str) -> str:
    lines = text.replace("\r\n", "\n").replace("\r", "\n").splitlines()
    return "\n".join(line.rstrip() for line in lines).rstrip() + "\n"


# Qualify one C++ line while preserving strings and comment state
def qualify_std_line_state(
    line: str, names: set[str], in_block_comment: bool
) -> tuple[str, bool]:
    if not in_block_comment and line.lstrip().startswith(("#", "//")):
        return line.rstrip(), False
    output: list[str] = []
    index = 0
    quote: str | None = None
    while index < len(line):
        if in_block_comment:
            end = line.find("*/", index)
            if end < 0:
                output.append(line[index:])
                return "".join(output).rstrip(), True
            output.append(line[index : end + 2])
            index = end + 2
            in_block_comment = False
            continue
        char = line[index]
        if quote is not None:
            output.append(char)
            if char == "\\" and index + 1 < len(line):
                index += 1
                output.append(line[index])
            elif char == quote:
                quote = None
            index += 1
            continue
        if char in ('"', "'"):
            quote = char
            output.append(char)
            index += 1
            continue
        if line.startswith("//", index):
            output.append(line[index:])
            break
        if line.startswith("/*", index):
            output.append("/*")
            index += 2
            in_block_comment = True
            continue
        if char.isalpha() or char == "_":
            end = index + 1
            while end < len(line) and (line[end].isalnum() or line[end] == "_"):
                end += 1
            token = line[index:end]
            prefix = "".join(output)
            if token in names and not prefix.endswith((".", "->", "::")):
                output.append("std::")
            output.append(token)
            index = end
            continue
        output.append(char)
        index += 1
    return "".join(output).rstrip(), in_block_comment


# Qualify C++ text while preserving multiline comments
def qualify_std_text(text: str, names: set[str]) -> str:
    lines = []
    in_block_comment = False
    for line in text.splitlines():
        qualified, in_block_comment = qualify_std_line_state(
            line, names, in_block_comment
        )
        lines.append(qualified)
    return "\n".join(lines) + ("\n" if text.endswith("\n") else "")


# Remove C++ comments while preserving literals, newlines and token boundaries
def strip_cpp_comments(text: str) -> str:
    output: list[str] = []
    quote: str | None = None
    in_block_comment = False
    in_line_comment = False
    pos = 0
    while pos < len(text):
        char = text[pos]
        if in_line_comment:
            if char == "\n":
                output.append(char)
                in_line_comment = False
            else:
                output.append(" ")
            pos += 1
            continue
        if in_block_comment:
            if text.startswith("*/", pos):
                output.extend((" ", " "))
                in_block_comment = False
                pos += 2
            else:
                output.append("\n" if char == "\n" else " ")
                pos += 1
            continue
        if quote is not None:
            output.append(char)
            if char == "\\" and pos + 1 < len(text):
                pos += 1
                output.append(text[pos])
            elif char == quote:
                quote = None
            pos += 1
            continue
        if text.startswith("//", pos):
            output.extend((" ", " "))
            in_line_comment = True
            pos += 2
            continue
        if text.startswith("/*", pos):
            output.extend((" ", " "))
            in_block_comment = True
            pos += 2
            continue
        if char in ('"', "'"):
            quote = char
        output.append(char)
        pos += 1
    return "".join(output)


# Mask C++ comments and literals while retaining code source offsets
def cpp_code_only(text: str) -> str:
    text = strip_cpp_comments(text)
    output: list[str] = []
    quote: str | None = None
    pos = 0
    while pos < len(text):
        char = text[pos]
        if quote is not None:
            if char == "\\" and pos + 1 < len(text):
                output.extend((" ", " "))
                pos += 2
                continue
            output.append("\n" if char == "\n" else " ")
            if char == quote:
                quote = None
            pos += 1
            continue
        if char in ('"', "'"):
            quote = char
            output.append(" ")
        else:
            output.append(char)
        pos += 1
    return "".join(output)


# Split the outer body of one C++ function into top-level statements
def cpp_function_statements(text: str, signature: str) -> tuple[str, ...]:
    _, _, block = function_block(text, signature)
    block = strip_cpp_comments(block)
    statements: list[str] = []
    current: list[str] = []
    quote: str | None = None
    curly = 0
    paren = 0
    square = 0
    pos = 0
    while pos < len(block):
        char = block[pos]
        if quote is not None:
            if curly > 0:
                current.append(char)
            if char == "\\" and pos + 1 < len(block):
                pos += 1
                if curly > 0:
                    current.append(block[pos])
                pos += 1
                continue
            if char == quote:
                quote = None
            pos += 1
            continue
        if char in ('"', "'"):
            quote = char
            if curly > 0:
                current.append(char)
            pos += 1
            continue
        if char == "{":
            curly += 1
            if curly > 1:
                current.append(char)
            pos += 1
            continue
        if char == "}":
            if curly > 1:
                current.append(char)
            curly -= 1
            pos += 1
            continue
        if curly == 0:
            pos += 1
            continue
        current.append(char)
        if char == "(":
            paren += 1
        elif char == ")":
            paren -= 1
        elif char == "[":
            square += 1
        elif char == "]":
            square -= 1
        elif char == ";" and curly == 1 and paren == 0 and square == 0:
            statement = "".join(current).strip()
            if statement:
                statements.append(statement)
            current = []
        pos += 1
    if curly != 0 or paren != 0 or square != 0 or quote is not None:
        raise RuntimeError(f"Unbalanced C++ function: {signature}")
    if "".join(current).strip():
        raise RuntimeError(f"Unterminated C++ statement in function: {signature}")
    return tuple(statements)


# Replace C++ identifiers outside string and character literals
def replace_cpp_identifiers(
    text: str, replacements: dict[str, str]
) -> tuple[str, set[str]]:
    output: list[str] = []
    replaced: set[str] = set()
    quote: str | None = None
    pos = 0
    while pos < len(text):
        char = text[pos]
        if quote is not None:
            output.append(char)
            if char == "\\" and pos + 1 < len(text):
                pos += 1
                output.append(text[pos])
            elif char == quote:
                quote = None
            pos += 1
            continue
        if char in ('"', "'"):
            quote = char
            output.append(char)
            pos += 1
            continue
        if char.isalpha() or char == "_":
            end = pos + 1
            while end < len(text) and (text[end].isalnum() or text[end] == "_"):
                end += 1
            token = text[pos:end]
            if token in replacements:
                output.append(replacements[token])
                replaced.add(token)
            else:
                output.append(token)
            pos = end
            continue
        output.append(char)
        pos += 1
    return "".join(output), replaced


# Collapse C++ whitespace outside string and character literals
def compact_cpp_whitespace(text: str) -> str:
    output: list[str] = []
    quote: str | None = None
    pending_space = False
    pos = 0
    while pos < len(text):
        char = text[pos]
        if quote is not None:
            output.append(char)
            if char == "\\" and pos + 1 < len(text):
                pos += 1
                output.append(text[pos])
            elif char == quote:
                quote = None
            pos += 1
            continue
        if char in ('"', "'"):
            if pending_space and output:
                output.append(" ")
            pending_space = False
            quote = char
            output.append(char)
        elif char.isspace():
            pending_space = True
        else:
            if pending_space and output:
                output.append(" ")
            pending_space = False
            output.append(char)
        pos += 1
    return "".join(output).strip()


# Replace broad generated standard imports with explicit qualifications
def narrow_std_namespace(text: str) -> str:
    text = re.sub(r"#include <std::(\w+)>", r"#include <\1>", text)
    text = text.replace(".std::", ".").replace("std::std::", "std::")
    text = text.replace("using namespace std;", "")
    text = re.sub(r"^using std::\w+;\s*$", "", text, flags=re.M)
    return qualify_std_text(text, STD_NAMES)
