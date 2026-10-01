// Unit tests for header-only compressed JSON numerical storage
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "catch.hpp"

// C++
#include <bit>
#include <complex>
#include <cstdint>
#include <limits>
#include <vector>

// Own
#include "Graniitti/Tech/MJsonZip.h"
#include "Graniitti/Math/MMatrix.h"

// Verify arbitrary byte payloads survive compression exactly
TEST_CASE("MJsonZip byte payloads round trip exactly", "[MJsonZip]") {
  std::vector<std::uint8_t> input(257);
  for (std::size_t index = 0; index < input.size(); ++index) {
    input[index] = static_cast<std::uint8_t>(index);
  }

  const auto compressed = gra::MJsonZip::CompressBytes(input);
  REQUIRE(compressed.is_string());

  std::vector<std::uint8_t> output;
  REQUIRE(gra::MJsonZip::DecompressBytes(compressed, output, input.size()));
  REQUIRE(output == input);
  REQUIRE_FALSE(
      gra::MJsonZip::DecompressBytes(compressed, output, input.size() + 1));
}

// Verify finite floating-point payloads preserve every IEEE 754 bit
TEST_CASE("MJsonZip finite double vectors round trip exactly", "[MJsonZip]") {
  const std::vector<double> input = {0.0,
                                     -0.0,
                                     1.0,
                                     -3.141592653589793,
                                     std::numeric_limits<double>::denorm_min(),
                                     std::numeric_limits<double>::min(),
                                     std::numeric_limits<double>::max()};

  const auto compressed = gra::MJsonZip::CompressVector(input);
  std::vector<double> output;
  REQUIRE(gra::MJsonZip::DecompressVector(compressed, output, input.size()));
  REQUIRE(output.size() == input.size());
  for (std::size_t index = 0; index < input.size(); ++index) {
    REQUIRE(std::bit_cast<std::uint64_t>(output[index]) ==
            std::bit_cast<std::uint64_t>(input[index]));
  }
}

// Verify complex matrix and matrix-bank storage is reusable outside eikonals
TEST_CASE("MJsonZip complex matrices round trip exactly", "[MJsonZip]") {
  using Complex = std::complex<double>;
  using Matrix = gra::MMatrix<Complex>;
  const Matrix first{{Complex(1.0, -2.0), Complex(-0.0, 0.25)},
                     {Complex(4.5, 6.0), Complex(-7.0, -8.0)}};
  const Matrix second{{Complex(-1.0, 0.5), Complex(2.0, 3.0)},
                      {Complex(0.125, -0.25), Complex(9.0, -4.0)}};

  Matrix matrix_output;
  REQUIRE(gra::MJsonZip::DecompressMatrix(gra::MJsonZip::CompressMatrix(first),
                                          matrix_output, 2, 2));
  REQUIRE(matrix_output.size_row() == first.size_row());
  REQUIRE(matrix_output.size_col() == first.size_col());
  REQUIRE(matrix_output.Flatten() == first.Flatten());

  const std::vector<Matrix> input{first, second};
  std::vector<Matrix> output;
  REQUIRE(gra::MJsonZip::DecompressMatrixBank(
      gra::MJsonZip::CompressMatrixBank(input), output, 2, 2, 2));
  REQUIRE(output.size() == input.size());
  for (std::size_t index = 0; index < input.size(); ++index) {
    REQUIRE(output[index].size_row() == input[index].size_row());
    REQUIRE(output[index].size_col() == input[index].size_col());
    REQUIRE(output[index].Flatten() == input[index].Flatten());
  }
}

// Verify malformed payloads and invalid numerical values are rejected
TEST_CASE("MJsonZip rejects invalid payloads", "[MJsonZip]") {
  std::vector<double> output{42.0};
  REQUIRE_FALSE(
      gra::MJsonZip::DecompressVector(nlohmann::json("!!!!"), output, 1));
  REQUIRE(output == std::vector<double>{42.0});
  REQUIRE_FALSE(gra::MJsonZip::DecompressVector(nlohmann::json(12), output, 1));
  REQUIRE_THROWS_AS(gra::MJsonZip::CompressVector(std::vector<double>{
                        std::numeric_limits<double>::infinity()}),
                    std::invalid_argument);

  using Matrix = gra::MMatrix<std::complex<double>>;
  const std::vector<Matrix> inconsistent{Matrix(1, 1), Matrix(2, 1)};
  REQUIRE_THROWS_AS(gra::MJsonZip::CompressMatrixBank(inconsistent),
                    std::invalid_argument);
}

// Check numerical storage rejects extra bytes after a complete compressed stream
TEST_CASE("MJsonZip consumes the complete compressed payload", "[MJsonZip]") {
  const std::vector<std::uint8_t> input{1, 2, 3, 4};
  const std::string compressed = gra::MJsonZip::CompressBytes(input).get<std::string>();
  REQUIRE(compressed.find('=') == std::string::npos);
  std::vector<std::uint8_t> output{42};
  CHECK_FALSE(gra::MJsonZip::DecompressBytes(compressed + "AAAA", output, input.size()));
  CHECK(output == std::vector<std::uint8_t>{42});
}
