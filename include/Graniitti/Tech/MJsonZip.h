// Header-only compressed JSON numerical storage
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MJSONZIP_H
#define MJSONZIP_H

// C++
#include <bit>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

// Libraries
#include "json.hpp"
#include <zlib.h>

namespace gra {

// Store exact numerical data as zlib-compressed Base64 JSON strings
class MJsonZip {
public:
  using Json = nlohmann::json;

  // Compress arbitrary bytes into one JSON string
  static Json CompressBytes(const std::span<const std::uint8_t> input) {
    return Base64Encode(ZlibCompress(input));
  }

  // Decompress arbitrary bytes and validate their exact expected size
  static bool DecompressBytes(const Json &value,
                              std::vector<std::uint8_t> &output,
                              const std::size_t expected_size) {
    if (!value.is_string()) {
      return false;
    }
    std::vector<std::uint8_t> compressed;
    std::vector<std::uint8_t> decoded;
    if (!Base64Decode(value.get_ref<const std::string &>(), compressed) ||
        !ZlibDecompress(compressed, decoded, expected_size)) {
      return false;
    }
    output = std::move(decoded);
    return true;
  }

  // Compress finite doubles using their exact IEEE 754 bit patterns
  static Json CompressVector(const std::span<const double> value) {
    static_assert(std::numeric_limits<double>::is_iec559);
    if (value.size() > std::numeric_limits<std::size_t>::max() / 8) {
      throw std::invalid_argument("MJsonZip vector is too large");
    }
    std::vector<std::uint8_t> bytes;
    bytes.reserve(8 * value.size());
    for (const double entry : value) {
      if (!std::isfinite(entry)) {
        throw std::invalid_argument(
            "MJsonZip vector contains non-finite value");
      }
      const std::uint64_t bits = std::bit_cast<std::uint64_t>(entry);
      for (int shift = 56; shift >= 0; shift -= 8) {
        bytes.push_back((bits >> shift) & 0xFFU);
      }
    }
    return CompressBytes(bytes);
  }

  // Decompress finite doubles and validate their exact expected size
  static bool DecompressVector(const Json &value, std::vector<double> &output,
                               const std::size_t expected_size) {
    static_assert(std::numeric_limits<double>::is_iec559);
    if (expected_size > std::numeric_limits<std::size_t>::max() / 8) {
      return false;
    }
    std::vector<std::uint8_t> bytes;
    if (!DecompressBytes(value, bytes, 8 * expected_size)) {
      return false;
    }
    std::vector<double> decoded(expected_size);
    for (std::size_t index = 0; index < expected_size; ++index) {
      std::uint64_t bits = 0;
      for (std::size_t offset = 0; offset < 8; ++offset) {
        bits = (bits << 8U) | bytes[8 * index + offset];
      }
      decoded[index] = std::bit_cast<double>(bits);
      if (!std::isfinite(decoded[index])) {
        return false;
      }
    }
    output = std::move(decoded);
    return true;
  }

  // Compress a finite complex vector into alternating real components
  static Json
  CompressComplexVector(const std::span<const std::complex<double>> value) {
    if (value.size() > std::numeric_limits<std::size_t>::max() / 2) {
      throw std::invalid_argument("MJsonZip complex vector is too large");
    }
    std::vector<double> data;
    data.reserve(2 * value.size());
    for (const auto &entry : value) {
      data.push_back(entry.real());
      data.push_back(entry.imag());
    }
    return CompressVector(data);
  }

  // Decompress a finite complex vector with one exact expected size
  static bool DecompressComplexVector(const Json &value,
                                      std::vector<std::complex<double>> &output,
                                      const std::size_t expected_size) {
    if (expected_size > std::numeric_limits<std::size_t>::max() / 2) {
      return false;
    }
    std::vector<double> data;
    if (!DecompressVector(value, data, 2 * expected_size)) {
      return false;
    }
    std::vector<std::complex<double>> decoded(expected_size, 0.0);
    for (std::size_t index = 0; index < expected_size; ++index) {
      decoded[index] = {data[2 * index], data[2 * index + 1]};
    }
    output = std::move(decoded);
    return true;
  }

  // Compress one finite complex matrix into one JSON string
  template <typename Matrix> static Json CompressMatrix(const Matrix &matrix) {
    const std::size_t rows = matrix.size_row();
    const std::size_t cols = matrix.size_col();
    if (!ValidComplexShape(rows, cols, 1)) {
      throw std::invalid_argument("MJsonZip matrix is too large or empty");
    }
    std::vector<double> data;
    data.reserve(2 * rows * cols);
    for (std::size_t row = 0; row < rows; ++row) {
      for (std::size_t col = 0; col < cols; ++col) {
        data.push_back(std::real(matrix(row, col)));
        data.push_back(std::imag(matrix(row, col)));
      }
    }
    return CompressVector(data);
  }

  // Decompress one finite complex matrix with expected dimensions
  template <typename Matrix>
  static bool DecompressMatrix(const Json &value, Matrix &matrix,
                               const std::size_t expected_rows,
                               const std::size_t expected_cols) {
    if (!ValidComplexShape(expected_rows, expected_cols, 1)) {
      return false;
    }
    std::vector<double> data;
    if (!DecompressVector(value, data, 2 * expected_rows * expected_cols)) {
      return false;
    }
    Matrix decoded(expected_rows, expected_cols);
    std::size_t index = 0;
    for (std::size_t row = 0; row < expected_rows; ++row) {
      for (std::size_t col = 0; col < expected_cols; ++col) {
        decoded(row, col) = std::complex<double>(data[index], data[index + 1]);
        index += 2;
      }
    }
    matrix = std::move(decoded);
    return true;
  }

  // Compress a bank of equally sized finite complex matrices
  template <typename Matrix>
  static Json CompressMatrixBank(const std::vector<Matrix> &bank) {
    if (bank.empty()) {
      return CompressVector(std::span<const double>());
    }
    const std::size_t rows = bank.front().size_row();
    const std::size_t cols = bank.front().size_col();
    if (!ValidComplexShape(rows, cols, bank.size())) {
      throw std::invalid_argument("MJsonZip matrix bank is too large or empty");
    }
    const std::size_t matrix_size = rows * cols;
    std::vector<double> data;
    data.reserve(2 * bank.size() * matrix_size);
    for (const Matrix &matrix : bank) {
      if (matrix.size_row() != rows || matrix.size_col() != cols) {
        throw std::invalid_argument(
            "MJsonZip matrix bank has inconsistent dimensions");
      }
      for (std::size_t row = 0; row < rows; ++row) {
        for (std::size_t col = 0; col < cols; ++col) {
          data.push_back(std::real(matrix(row, col)));
          data.push_back(std::imag(matrix(row, col)));
        }
      }
    }
    return CompressVector(data);
  }

  // Decompress a bank of finite complex matrices with expected dimensions
  template <typename Matrix>
  static bool DecompressMatrixBank(const Json &value, std::vector<Matrix> &bank,
                                   const std::size_t expected_size,
                                   const std::size_t rows,
                                   const std::size_t cols) {
    if (!ValidComplexShape(rows, cols, expected_size)) {
      return false;
    }
    const std::size_t matrix_size = rows * cols;
    std::vector<double> data;
    if (!DecompressVector(value, data, 2 * expected_size * matrix_size)) {
      return false;
    }
    std::vector<Matrix> decoded;
    decoded.reserve(expected_size);
    std::size_t data_index = 0;
    for (std::size_t index = 0; index < expected_size; ++index) {
      decoded.emplace_back(rows, cols);
      for (std::size_t row = 0; row < rows; ++row) {
        for (std::size_t col = 0; col < cols; ++col) {
          decoded.back()(row, col) =
              std::complex<double>(data[data_index], data[data_index + 1]);
          data_index += 2;
        }
      }
    }
    bank = std::move(decoded);
    return true;
  }

private:
  // Compute whether one complex matrix bank shape has a safe nonzero size
  static bool ValidComplexShape(const std::size_t rows, const std::size_t cols,
                                const std::size_t count) {
    if (rows == 0 || cols == 0 ||
        rows > std::numeric_limits<std::size_t>::max() / cols) {
      return false;
    }
    const std::size_t matrix_size = rows * cols;
    return matrix_size <= std::numeric_limits<std::size_t>::max() / 2 &&
           count <= std::numeric_limits<std::size_t>::max() / (2 * matrix_size);
  }

  // Compute one Base64 sextet value or -1 for an invalid character
  static int Base64Value(const char character) {
    if (character >= 'A' && character <= 'Z') {
      return character - 'A';
    }
    if (character >= 'a' && character <= 'z') {
      return character - 'a' + 26;
    }
    if (character >= '0' && character <= '9') {
      return character - '0' + 52;
    }
    if (character == '+') {
      return 62;
    }
    if (character == '/') {
      return 63;
    }
    return -1;
  }

  // Encode bytes as compact Base64 text
  static std::string Base64Encode(const std::vector<std::uint8_t> &input) {
    constexpr std::string_view alphabet =
        "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    if (input.size() > std::numeric_limits<std::size_t>::max() - 2) {
      throw std::invalid_argument("MJsonZip Base64 input is too large");
    }
    const std::size_t groups = (input.size() + 2) / 3;
    if (groups > std::numeric_limits<std::size_t>::max() / 4) {
      throw std::invalid_argument("MJsonZip Base64 input is too large");
    }
    std::string output;
    output.reserve(4 * groups);
    for (std::size_t index = 0; index < input.size(); index += 3) {
      const std::uint32_t first = input[index];
      const std::uint32_t second =
          index + 1 < input.size() ? input[index + 1] : 0;
      const std::uint32_t third =
          index + 2 < input.size() ? input[index + 2] : 0;
      const std::uint32_t block = (first << 16U) | (second << 8U) | third;
      output.push_back(alphabet[(block >> 18U) & 0x3FU]);
      output.push_back(alphabet[(block >> 12U) & 0x3FU]);
      output.push_back(
          index + 1 < input.size() ? alphabet[(block >> 6U) & 0x3FU] : '=');
      output.push_back(index + 2 < input.size() ? alphabet[block & 0x3FU]
                                                : '=');
    }
    return output;
  }

  // Decode compact Base64 text into bytes
  static bool Base64Decode(const std::string &input,
                           std::vector<std::uint8_t> &output) {
    if (input.size() % 4 != 0) {
      return false;
    }
    std::vector<std::uint8_t> decoded;
    decoded.reserve(3 * (input.size() / 4));
    for (std::size_t index = 0; index < input.size(); index += 4) {
      const int first = Base64Value(input[index]);
      const int second = Base64Value(input[index + 1]);
      const bool pad2 = input[index + 2] == '=';
      const bool pad3 = input[index + 3] == '=';
      const int third = pad2 ? 0 : Base64Value(input[index + 2]);
      const int fourth = pad3 ? 0 : Base64Value(input[index + 3]);
      if (first < 0 || second < 0 || third < 0 || fourth < 0 ||
          (pad2 && !pad3) || ((pad2 || pad3) && index + 4 != input.size())) {
        return false;
      }
      const std::uint32_t block = (static_cast<std::uint32_t>(first) << 18U) |
                                  (static_cast<std::uint32_t>(second) << 12U) |
                                  (static_cast<std::uint32_t>(third) << 6U) |
                                  static_cast<std::uint32_t>(fourth);
      decoded.push_back((block >> 16U) & 0xFFU);
      if (!pad2) {
        decoded.push_back((block >> 8U) & 0xFFU);
      }
      if (!pad3) {
        decoded.push_back(block & 0xFFU);
      }
    }
    output = std::move(decoded);
    return true;
  }

  // Compress one exact byte sequence with maximum lossless compression
  static std::vector<std::uint8_t>
  ZlibCompress(const std::span<const std::uint8_t> input) {
    if (input.empty()) {
      return {};
    }
    if (input.size() > std::numeric_limits<uLong>::max()) {
      throw std::invalid_argument("MJsonZip payload is too large");
    }
    const uLong source_size = static_cast<uLong>(input.size());
    const uLongf bound = compressBound(source_size);
    if (bound > std::numeric_limits<std::size_t>::max()) {
      throw std::invalid_argument("MJsonZip payload is too large");
    }
    std::vector<std::uint8_t> compressed(bound);
    uLongf compressed_size = bound;
    if (compress2(compressed.data(), &compressed_size, input.data(),
                  source_size, Z_BEST_COMPRESSION) != Z_OK) {
      throw std::runtime_error("MJsonZip compression failed");
    }
    compressed.resize(compressed_size);
    return compressed;
  }

  // Decompress one exact byte sequence to its validated size
  static bool ZlibDecompress(const std::vector<std::uint8_t> &compressed,
                             std::vector<std::uint8_t> &output,
                             const std::size_t expected_size) {
    if (expected_size == 0) {
      output.clear();
      return compressed.empty();
    }
    if (compressed.empty() ||
        expected_size > std::numeric_limits<uLongf>::max() ||
        compressed.size() > std::numeric_limits<uLong>::max()) {
      return false;
    }
    std::vector<std::uint8_t> decoded(expected_size);
    uLongf output_size = static_cast<uLongf>(expected_size);
    uLong input_size = static_cast<uLong>(compressed.size());
    const int status =
        uncompress2(decoded.data(), &output_size, compressed.data(), &input_size);
    if (status != Z_OK || output_size != expected_size || input_size != compressed.size()) {
      return false;
    }
    output = std::move(decoded);
    return true;
  }
};

} // namespace gra

#endif // MJSONZIP_H
