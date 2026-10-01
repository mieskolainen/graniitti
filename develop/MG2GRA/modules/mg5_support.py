#!/usr/bin/env python3
#
# Transform and validate model support from MG5 standalone C++ exports
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import re
from dataclasses import dataclass

from . import output_layout
from .cpp_support import function_block, narrow_std_namespace
from .mg5_couplings import alpha_zero_definition


@dataclass(frozen=True)
class ModelFiles:
    """Names of the generated support classes and files for one UFO model"""

    parameter_class: str
    parameter_header: str
    helas_header: str
    model_suffix: str


# Discover model-specific support names from a generated process pair
def discover_process_model(raw_header: str, raw_source: str) -> ModelFiles:
    prefix = r"(?:Graniitti/Amplitude/MG5/(?:[A-Za-z0-9_]+/)*)?"
    parameter = re.search(rf'#include "{prefix}(Parameters_([^"]+))\.h"', raw_header)
    helas = re.search(rf'#include "{prefix}(HelAmps_([^"]+))\.h"', raw_source)
    if parameter is None or helas is None:
        raise RuntimeError("Could not discover generated MadGraph model support")
    if parameter.group(2) != helas.group(2):
        raise RuntimeError(
            f"Parameter model {parameter.group(2)} differs from HELAS model {helas.group(2)}"
        )
    return ModelFiles(
        parameter_class=parameter.group(1),
        parameter_header=f"{parameter.group(1)}.h",
        helas_header=f"{helas.group(1)}.h",
        model_suffix=parameter.group(2),
    )


# Convert generated HELAS declarations to repository-safe C++
def transform_helas_header(raw_header: str) -> str:
    return narrow_std_namespace(raw_header.replace("\r\n", "\n"))


# Convert generated HELAS definitions to repository-safe C++
def transform_helas_source(
    raw_source: str, helas_header: str, include_directory: str
) -> str:
    source = raw_source.replace("\r\n", "\n")
    source = source.replace(
        f'#include "{helas_header}"',
        f'#include "{output_layout.include_path(include_directory, helas_header)}"',
    )
    # Rationalize massive and massless spinors near the negative beam axis
    for momentum in ("pp", "p[0]"):
        source = re.sub(
            rf"(?:std::)?max\({re.escape(momentum)}\s*\+\s*p\[3\],\s*0\.0+\)",
            "(p[3] < 0.0 ? (p[1] * p[1] + p[2] * p[2]) / "
            f"({momentum} - p[3]) : {momentum} + p[3])",
            source,
        )
    return narrow_std_namespace(source)


# Compute whether generated couplings use the configured electromagnetic charge
def supports_alpha_qed_zero(
    raw_source: str,
    parameter_class: str,
    charge: str | None,
    charge_square: str | None,
) -> bool:
    try:
        return (
            alpha_zero_definition(
                raw_source, parameter_class, charge, charge_square, 1.0
            )
            is not None
        )
    except RuntimeError as error:
        if "No generated coupling depends" not in str(error):
            raise
        return False


# Convert a generated parameter header to value-owned, thread-safe state
def transform_parameters_header(raw_header: str, parameter_class: str, alpha_zero: bool) -> str:
    header = raw_header.replace("\r\n", "\n")
    header = header.replace(
        '#include "read_slha.h"',
        f'#include "{output_layout.include_path(output_layout.RUNTIME, "read_slha.h")}"',
    )
    header = re.sub(
        rf"^\s*static\s+{re.escape(parameter_class)}\s*\*\s*getInstance\(\);\s*$",
        "",
        header,
        flags=re.M,
    )
    header = re.sub(
        rf"^\s*static\s+{re.escape(parameter_class)}\s*\*\s*instance;\s*$",
        "",
        header,
        flags=re.M,
    )
    header, replacements = re.subn(
        r"void\s+setDependentParameters\(\s*\)\s*;",
        "void setDependentParameters(double alpS);",
        header,
        count=1,
    )
    if replacements != 1:
        raise RuntimeError("Could not update generated dependent-parameter declaration")
    if alpha_zero:
        header, replacements = re.subn(
            r"(void\s+setDependentCouplings\(\s*\)\s*;)",
            r"\1\n    // Set electromagnetic couplings at Q2 = 0\n"
            r"    void setAlphaQEDZero();",
            header,
            count=1,
        )
        if replacements != 1:
            raise RuntimeError("Could not add alpha(0) parameter declaration")
    return narrow_std_namespace(header)


# Convert generated parameter definitions to value-owned, thread-safe state
def transform_parameters_source(
    raw_source: str,
    parameter_class: str,
    parameter_header: str,
    charge: str | None,
    charge_square: str | None,
    inverse_alpha: float,
    include_directory: str,
) -> str:
    source = raw_source.replace("\r\n", "\n")
    source = source.replace(
        f'#include "{parameter_header}"',
        f'#include "{output_layout.include_path(include_directory, parameter_header)}"',
    )
    source = re.sub(
        rf"^\s*{re.escape(parameter_class)}\s*\*\s*{re.escape(parameter_class)}::instance\s*=\s*0\s*;\s*$",
        "",
        source,
        flags=re.M,
    )
    source = re.sub(
        r"^\s*//\s*(?:Initialize|Function to get) static instance.*$",
        "",
        source,
        flags=re.M,
    )
    singleton_signature = f"{parameter_class} * {parameter_class}::getInstance()"
    if singleton_signature in source:
        start, end, _ = function_block(source, singleton_signature)
        source = source[:start] + source[end:]

    signature = f"void {parameter_class}::setDependentParameters()"
    start, end, block = function_block(source, signature)
    updated = block.replace(
        signature, f"void {parameter_class}::setDependentParameters(double alpS)", 1
    )
    opening = updated.find("{")
    assignment = "\n  aS = alpS;" if re.search(r"\baS\b", updated) else "\n  (void)alpS;"
    updated = updated[: opening + 1] + assignment + updated[opening + 1 :]
    source = source[:start] + updated + source[end:]

    alpha_definition = alpha_zero_definition(
        raw_source,
        parameter_class,
        charge,
        charge_square,
        inverse_alpha,
    )
    if alpha_definition is not None:
        dependent_signature = f"void {parameter_class}::setDependentCouplings()"
        _, insert_at, _ = function_block(source, dependent_signature)
        source = source[:insert_at] + "\n" + alpha_definition + source[insert_at:]

    if "#include <cmath>" not in source:
        source = source.replace("#include <iostream>", "#include <cmath>\n#include <iostream>", 1)
    return narrow_std_namespace(source)


# Convert the common generated SLHA reader header
def transform_slha_header(raw_header: str | None = None) -> str:
    if raw_header is not None and ("class SLHABlock" not in raw_header or "class SLHAReader" not in raw_header):
        raise RuntimeError("Generated SLHA header has an unsupported structure")
    return """// Read SLHA parameter cards for generated MG5 amplitudes
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef READ_SLHA_H
#define READ_SLHA_H

#include <cstddef>
#include <map>
#include <optional>
#include <string>
#include <vector>

class SLHABlock {
 public:
  // Construct one named SLHA block
  explicit SLHABlock(std::string name = "");

  // Store one finite value under a fixed-rank integer index
  void set_entry(const std::vector<int> &indices, double value);

  // Compute one block entry or the requested default
  double get_entry(const std::vector<int> &indices,
                   double default_value = 0.0) const;

  // Compute one block entry without emitting a missing-entry warning
  std::optional<double> find_entry(const std::vector<int> &indices) const;

  // Replace the normalized block name
  void set_name(std::string name);

  // Compute the normalized block name
  const std::string &get_name() const;

  // Compute the index rank after the first stored entry
  std::size_t get_indices() const;

 private:
  std::string _name;
  std::map<std::vector<int>, double> _entries;
  std::size_t _indices = 0;
  bool _has_indices = false;
};

class SLHAReader {
 public:
  // Construct an empty reader or load one parameter card
  explicit SLHAReader(const std::string &file_name = "");

  // Parse one complete parameter card atomically
  void read_slha_file(const std::string &file_name);

  // Compute one indexed block entry or the requested default
  double get_block_entry(const std::string &block_name,
                         const std::vector<int> &indices,
                         double default_value = 0.0) const;

  // Compute one single-index block entry or the requested default
  double get_block_entry(const std::string &block_name, int index,
                         double default_value = 0.0) const;

  // Compute one indexed block entry without emitting a missing-entry warning
  std::optional<double>
  find_block_entry(const std::string &block_name,
                   const std::vector<int> &indices) const;

  // Compute one single-index block entry without a missing-entry warning
  std::optional<double> find_block_entry(const std::string &block_name,
                                         int index) const;

  // Store one indexed block entry
  void set_block_entry(const std::string &block_name,
                       const std::vector<int> &indices, double value);

  // Store one single-index block entry
  void set_block_entry(const std::string &block_name, int index,
                       double value);

 private:
  std::map<std::string, SLHABlock> _blocks;
};

#endif
"""


# Convert the common generated SLHA reader implementation
def transform_slha_source(raw_source: str | None = None) -> str:
    required = ("SLHABlock::set_entry", "SLHAReader::read_slha_file")
    if raw_source is not None and any(signature not in raw_source for signature in required):
        raise RuntimeError("Generated SLHA source has an unsupported structure")
    return """// Read SLHA parameter cards for generated MG5 amplitudes
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "@READ_SLHA_INCLUDE@"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace {

// Compute one string without surrounding ASCII whitespace
std::string Trim(std::string text) {
  const std::size_t first = text.find_first_not_of(" \\t\\r\\n");
  if (first == std::string::npos) { return {}; }
  const std::size_t last = text.find_last_not_of(" \\t\\r\\n");
  return text.substr(first, last - first + 1);
}

// Normalize one SLHA keyword or block name to lowercase
std::string Lowercase(std::string text) {
  std::transform(text.begin(), text.end(), text.begin(), [](unsigned char value) {
    return static_cast<char>(std::tolower(value));
  });
  return text;
}

// Parse one complete finite SLHA floating-point token
double ParseValue(std::string token, std::size_t line_number) {
  std::replace(token.begin(), token.end(), 'd', 'e');
  std::replace(token.begin(), token.end(), 'D', 'E');
  std::size_t consumed = 0;
  double value = 0.0;
  try {
    value = std::stod(token, &consumed);
  } catch (const std::exception &) {
    throw std::runtime_error("Invalid SLHA value on line " +
                             std::to_string(line_number));
  }
  if (consumed != token.size() || !std::isfinite(value)) {
    throw std::runtime_error("Invalid SLHA value on line " +
                             std::to_string(line_number));
  }
  return value;
}

// Parse one complete in-range SLHA integer index
int ParseIndex(const std::string &token, std::size_t line_number) {
  std::size_t consumed = 0;
  long long value = 0;
  try {
    value = std::stoll(token, &consumed);
  } catch (const std::exception &) {
    throw std::runtime_error("Invalid SLHA index on line " +
                             std::to_string(line_number));
  }
  if (consumed != token.size() ||
      value < std::numeric_limits<int>::min() ||
      value > std::numeric_limits<int>::max()) {
    throw std::runtime_error("Invalid SLHA index on line " +
                             std::to_string(line_number));
  }
  return static_cast<int>(value);
}

} // namespace

// Construct one named SLHA block
SLHABlock::SLHABlock(std::string name) : _name(std::move(name)) {}

// Store one finite value under a fixed-rank integer index
void SLHABlock::set_entry(const std::vector<int> &indices, double value) {
  if (!std::isfinite(value)) {
    throw std::invalid_argument("SLHABlock::set_entry: invalid entry");
  }
  if (!_has_indices) {
    _indices = indices.size();
    _has_indices = true;
  } else if (indices.size() != _indices) {
    throw std::invalid_argument(
        "SLHABlock::set_entry: inconsistent index rank");
  }
  _entries[indices] = value;
}

// Compute one block entry or the requested default
double SLHABlock::get_entry(const std::vector<int> &indices,
                            double default_value) const {
  const auto value = find_entry(indices);
  if (!value.has_value()) {
    std::cerr << "Warning: No such entry in " << _name
              << ", using default value " << default_value << '\\n';
    return default_value;
  }
  return *value;
}

// Compute one block entry without emitting a missing-entry warning
std::optional<double>
SLHABlock::find_entry(const std::vector<int> &indices) const {
  const auto found = _entries.find(indices);
  return found == _entries.end() ? std::nullopt
                                 : std::optional<double>{found->second};
}

// Replace the normalized block name
void SLHABlock::set_name(std::string name) { _name = std::move(name); }

// Compute the normalized block name
const std::string &SLHABlock::get_name() const { return _name; }

// Compute the index rank after the first stored entry
std::size_t SLHABlock::get_indices() const { return _indices; }

// Construct an empty reader or load one parameter card
SLHAReader::SLHAReader(const std::string &file_name) {
  if (!file_name.empty()) { read_slha_file(file_name); }
}

// Parse one complete parameter card atomically
void SLHAReader::read_slha_file(const std::string &file_name) {
  std::ifstream input(file_name);
  if (!input) {
    throw std::runtime_error("SLHAReader: cannot open parameter card " +
                             file_name);
  }

  SLHAReader parsed;
  std::string active_block;
  std::string line;
  std::size_t line_number = 0;
  while (std::getline(input, line)) {
    ++line_number;
    const std::size_t comment = line.find('#');
    if (comment != std::string::npos) { line.erase(comment); }
    line = Trim(line);
    if (line.empty()) { continue; }

    std::istringstream stream(line);
    std::vector<std::string> fields;
    for (std::string field; stream >> field;) { fields.push_back(field); }
    if (fields.empty()) { continue; }

    const std::string keyword = Lowercase(fields.front());
    if (keyword == "block") {
      if (fields.size() < 2) {
        throw std::runtime_error("Missing SLHA block name on line " +
                                 std::to_string(line_number));
      }
      active_block = Lowercase(fields[1]);
      continue;
    }
    if (keyword == "decay") {
      if (fields.size() < 3) {
        throw std::runtime_error("Invalid SLHA decay on line " +
                                 std::to_string(line_number));
      }
      const int pdg = ParseIndex(fields[1], line_number);
      const double width = ParseValue(fields[2], line_number);
      parsed.set_block_entry("decay", pdg, width);
      active_block.clear();
      continue;
    }
    if (active_block.empty()) {
      continue;
    }
    std::vector<int> indices;
    indices.reserve(fields.size() - 1);
    for (std::size_t index = 0; index + 1 < fields.size(); ++index) {
      indices.push_back(ParseIndex(fields[index], line_number));
    }
    parsed.set_block_entry(active_block, indices,
                           ParseValue(fields.back(), line_number));
  }
  if (input.bad()) {
    throw std::runtime_error("SLHAReader: failed while reading parameter card " +
                             file_name);
  }
  if (parsed._blocks.empty()) {
    throw std::runtime_error("SLHAReader: parameter card contains no data");
  }
  _blocks.swap(parsed._blocks);
}

// Compute one indexed block entry or the requested default
double SLHAReader::get_block_entry(const std::string &block_name,
                                   const std::vector<int> &indices,
                                   double default_value) const {
  const auto found = _blocks.find(Lowercase(block_name));
  if (found == _blocks.end()) {
    std::cerr << "Warning: No such block " << block_name
              << ", using default value " << default_value << '\\n';
    return default_value;
  }
  return found->second.get_entry(indices, default_value);
}

// Compute one single-index block entry or the requested default
double SLHAReader::get_block_entry(const std::string &block_name, int index,
                                   double default_value) const {
  return get_block_entry(block_name, std::vector<int>{index}, default_value);
}

// Compute one indexed block entry without emitting a missing-entry warning
std::optional<double> SLHAReader::find_block_entry(
    const std::string &block_name, const std::vector<int> &indices) const {
  const auto found = _blocks.find(Lowercase(block_name));
  return found == _blocks.end() ? std::nullopt
                                : found->second.find_entry(indices);
}

// Compute one single-index block entry without a missing-entry warning
std::optional<double>
SLHAReader::find_block_entry(const std::string &block_name, int index) const {
  return find_block_entry(block_name, std::vector<int>{index});
}

// Store one indexed block entry
void SLHAReader::set_block_entry(const std::string &block_name,
                                 const std::vector<int> &indices,
                                 double value) {
  const std::string normalized = Lowercase(block_name);
  auto [found, inserted] = _blocks.try_emplace(normalized, normalized);
  (void)inserted;
  found->second.set_entry(indices, value);
}

// Store one single-index block entry
void SLHAReader::set_block_entry(const std::string &block_name, int index,
                                 double value) {
  set_block_entry(block_name, std::vector<int>{index}, value);
}
""".replace(
        "@READ_SLHA_INCLUDE@",
        output_layout.include_path(output_layout.RUNTIME, "read_slha.h"),
    )


# Validate that every generated process reference is provided by aggregate support
def validate_dependency_union(
    process_sources: list[str],
    helas_header: str,
    parameters_header: str,
    model_suffix: str,
) -> None:
    provided_helas = set(re.findall(r"\b(?:void|double)\s+([A-Za-z_]\w*)\s*\(", helas_header))
    provided_parameters = set(re.findall(r"\b([A-Za-z_]\w*)\b", parameters_header))
    missing_helas: set[str] = set()
    missing_parameters: set[str] = set()
    namespace = f"MG5_{model_suffix}"
    for source in process_sources:
        missing_helas.update(
            name
            for name in re.findall(rf"\b{re.escape(namespace)}::([A-Za-z_]\w*)\s*\(", source)
            if name not in provided_helas
        )
        missing_parameters.update(
            name
            for name in re.findall(r"\bpars(?:->|\.)\s*([A-Za-z_]\w*)", source)
            if name not in provided_parameters
        )
    if missing_helas or missing_parameters:
        raise RuntimeError(
            "Aggregate MG5 support is incomplete: "
            f"HELAS={sorted(missing_helas)}, parameters={sorted(missing_parameters)}"
        )


# Expose pole parameters directly from the UFO particle definitions
def add_model_particles(header: str, particles: list[dict], source: str, charge: str | None = None) -> str:
    signature = re.search(r"void \w+::setIndependentParameters\(", source)
    if signature is None:
        raise RuntimeError("Missing independent model parameter initialization")
    _, _, independent = function_block(source, signature[0][:-1])
    assigned = set(re.findall(r"\b(\w+)\s*=(?!=)", independent))
    if charge is not None and (charge not in assigned or re.search(rf"\b{re.escape(charge)}\b", header) is None):
        raise RuntimeError(f"Missing generated electromagnetic charge parameter {charge}")
    rows = []
    for particle in particles:
        values = [particle[key] for key in ("mass", "width")]
        if any(value not in assigned or re.search(rf"\b{re.escape(value)}\b", header) is None
               for value in values):
            raise RuntimeError(f"Missing generated pole parameter for PDG {particle['pdg']}")
        if particle.get("signed_mass", False):
            rows.append(f"      {{{particle['pdg']}, {{std::abs({values[0]}), std::abs({values[1]}), true}}}}")
        else:
            rows.append(f"      {{{particle['pdg']}, {{{values[0]}, {values[1]}}}}}")
    inverse = re.search(r'(\w+)\s*=\s*slha\.get_block_entry\("sminputs",\s*1\s*,', source, re.I)
    alpha = f"std::norm({charge}) / (4.0 * M_PI)" if charge else (f"1.0 / {inverse[1]}" if inverse else "0.0")
    method = (f"\n    // Compute the evaluated UFO electromagnetic coupling\n"
              f"    double AlphaQED() const {{ return {alpha}; }}\n"
              "\n    // Compute masses and widths from the evaluated UFO parameters\n"
              "    gra::mg5::ParticleMap Particles() const {\n      return {\n" +
              ",\n".join(rows) + "\n      };\n    }\n")
    header, count = re.subn(r"(public\s*:)", lambda m: m[0] + method, header, count=1)
    if count != 1:
        raise RuntimeError("Missing generated model public declarations")
    return header.replace('#include <complex>', '#include <complex>\n'
                          '#include "Graniitti/Amplitude/MG5/Runtime/AMP_MG5_Model.h"', 1)
