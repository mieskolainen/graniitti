// Standard Model light by light amplitude tests
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include <array>
#include <cmath>
#include <complex>
#include <future>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <catch.hpp>

#include "Graniitti/Amplitude/Photon/AMP_yy_yy.h"
#include "Graniitti/Process/MSubProc.h"
#include "Graniitti/Tech/MAux.h"

namespace {
using gra::aux::indices;

// Evaluate one physical planar hard point
std::array<std::complex<double>, 16> Evaluate(const double energy,
                                              const double cosine, const gra::SMParam &sm) {
  const double s = 4.0 * energy * energy;
  const double t = -0.5 * s * (1.0 - cosine);
  const double u = -0.5 * s * (1.0 + cosine);
  std::array<std::complex<double>, 16> amplitude;
  if (!gra::lbyl::SMHelicityAmplitudes(s, t, u, 1.0 / sm.alpha_em_inv, sm, amplitude)) {
    throw std::runtime_error("analytic light-by-light evaluation failed");
  }
  return amplitude;
}

// Exchange two helicity bits in one lexicographic row
std::size_t SwapHelicities(std::size_t row, const std::size_t first,
                           const std::size_t second) {
  const std::size_t first_mask = std::size_t{1} << (3 - first);
  const std::size_t second_mask = std::size_t{1} << (3 - second);
  if (((row & first_mask) != 0) != ((row & second_mask) != 0)) {
    row ^= first_mask | second_mask;
  }
  return row;
}

} // namespace

// Check parity and identical photon exchange in the massive helicity amplitudes
TEST_CASE("Analytic light by light obeys parity and photon crossing",
          "[gra::lbyl][symmetry]") {
  const auto tune = gra::MModelTune::Load(gra::aux::ResolveProjectPath("modeldata/TUNE0/GENERAL.json"));
  const auto direct = Evaluate(50.0, 0.37, tune->SM());
  const auto crossed = Evaluate(50.0, -0.37, tune->SM());
  for (const auto &row : indices(direct)) {
    REQUIRE(direct[row].real() ==
            Approx(direct[15 - row].real()).margin(2.0e-13));
    REQUIRE(direct[row].imag() ==
            Approx(direct[15 - row].imag()).margin(2.0e-13));
    const std::size_t exchange = SwapHelicities(row, 2, 3);
    REQUIRE(direct[row].real() ==
            Approx(crossed[exchange].real()).margin(2.0e-12));
    REQUIRE(direct[row].imag() ==
            Approx(crossed[exchange].imag()).margin(2.0e-12));
  }
}

// Check finite amplitudes across the charged particle thresholds
TEST_CASE("Analytic light by light is finite across charged thresholds",
          "[gra::lbyl][threshold]") {
  const auto tune = gra::MModelTune::Load(gra::aux::ResolveProjectPath("modeldata/TUNE0/GENERAL.json"));
  const auto &sm = tune->SM();
  for (const double energy : {0.2, 1.0, sm.c * 0.9999, sm.c * 1.0001, sm.b * 0.9999, sm.b * 1.0001,
                              sm.w * 0.9999, sm.w * 1.0001, sm.t * 0.9999, sm.t * 1.0001}) {
    const auto amplitude = Evaluate(energy, -0.23, sm);
    for (const auto &value : amplitude) {
      REQUIRE(std::isfinite(value.real()));
      REQUIRE(std::isfinite(value.imag()));
    }
  }
}

// Check concurrent evaluation preserves each helicity amplitude
TEST_CASE("Analytic light by light has no shared mutable worker state",
          "[gra::lbyl][thread]") {
  const auto tune = gra::MModelTune::Load(gra::aux::ResolveProjectPath("modeldata/TUNE0/GENERAL.json"));
  const auto sm = tune->SM();
  std::vector<std::future<std::array<std::complex<double>, 16>>> workers;
  constexpr std::array<int, 8> slots = {};
  for (const auto &i : indices(slots)) {
    workers.push_back(std::async(std::launch::async,
                                 [i, sm]() { return Evaluate(30.0 + i, 0.21, sm); }));
  }
  for (const auto &i : indices(workers)) {
    const auto concurrent = workers[i].get();
    const auto serial = Evaluate(30.0 + i, 0.21, sm);
    for (const auto &row : indices(serial)) {
      REQUIRE(std::abs(concurrent[row] - serial[row]) < 1.0e-13);
    }
  }
}

