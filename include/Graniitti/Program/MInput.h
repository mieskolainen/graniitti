// Strict numerical program input
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef PROGRAM_MINPUT_H
#define PROGRAM_MINPUT_H

#include <charconv>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include "json.hpp"

namespace gra::program {

// Parse a complete finite numerical command line value
template <typename T>
T Number(const std::string& text) {
  T          value{};
  const auto parsed = std::from_chars(text.data(), text.data() + text.size(), value);
  if (parsed.ec != std::errc{} || parsed.ptr != text.data() + text.size()) {
    throw std::invalid_argument("Invalid numerical input: " + text);
  }
  if constexpr (std::is_floating_point_v<T>) {
    if (!std::isfinite(value)) { throw std::invalid_argument("Non-finite input: " + text); }
  }
  return value;
}

// Parse a nonempty comma separated list of positive finite collision energies
inline std::vector<double> Energies(const std::string& text) {
  std::vector<double> energy;
  std::size_t begin = 0;
  while (true) {
    const auto end = text.find(',', begin);
    auto token = text.substr(begin, end == std::string::npos ? end : end - begin);
    const auto first = token.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) { throw std::invalid_argument("Empty collision energy"); }
    token = token.substr(first, token.find_last_not_of(" \t\r\n") - first + 1);
    const auto value = Number<double>(token);
    if (!(value > 0.0)) { throw std::invalid_argument("Collision energies must be positive"); }
    energy.push_back(value);
    if (end == std::string::npos) { break; }
    begin = end + 1;
  }
  return energy;
}

// Prevent an output file from truncating the input through a path or hard link alias
inline void DistinctFiles(const std::string& input, const std::string& output) {
  if (std::filesystem::exists(output) && std::filesystem::equivalent(input, output)) {
    throw std::invalid_argument("Input and output refer to the same file: " + input);
  }
}

// Read a nonnegative integer without truncation or unsigned wraparound
template <typename T>
T Count(const nlohmann::json& value, const std::string& name) {
  static_assert(std::is_integral_v<T> && !std::is_same_v<T, bool>);
  if (!value.is_number_integer() || (!value.is_number_unsigned() && value.get<std::int64_t>() < 0)) {
    throw std::invalid_argument(name + " must be a nonnegative integer");
  }
  const auto count = value.get<std::uint64_t>();
  if (count > static_cast<std::uint64_t>(std::numeric_limits<T>::max())) {
    throw std::invalid_argument(name + " exceeds integer range");
  }
  return static_cast<T>(count);
}

}  // namespace gra::program
#endif
