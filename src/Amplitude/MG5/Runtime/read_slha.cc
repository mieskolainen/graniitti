// Read SLHA parameter cards for generated MG5 amplitudes
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Amplitude/MG5/Runtime/read_slha.h"

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
  const std::size_t first = text.find_first_not_of(" \t\r\n");
  if (first == std::string::npos) { return {}; }
  const std::size_t last = text.find_last_not_of(" \t\r\n");
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
              << ", using default value " << default_value << '\n';
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
              << ", using default value " << default_value << '\n';
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
