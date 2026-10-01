// Numerical and event conversion regression tests for standalone programs
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include <array>
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <limits>
#include <numeric>
#include <string>
#include <valarray>

#include "Graniitti/Program/Analysis/fitharmonic.h"
#include "Graniitti/Program/data2hepmc3.h"
#include "Graniitti/Program/hepevt2hepmc3.h"
#include "Graniitti/Program/minbias.h"
#include "Graniitti/Program/ot.h"
#include "Graniitti/Program/pathmark.h"
#include "Graniitti/Program/pdebench.h"
#include "Graniitti/Program/sommerfeld.h"
#include "Graniitti/Program/xscan.h"
#include "HepMC3/WriterHEPEVT.h"
#include "catch.hpp"
#include "support/analysis_test_support.hh"

using gra::aux::indices;

namespace {

// Keep numerical program outputs in the normal project temporary directory
std::string ProgramPath(const std::string& name) {
  const auto directory = std::filesystem::path(gra::aux::ResolveProjectPath("tmp")) / "test_programs";
  std::filesystem::create_directories(directory);
  return (directory / name).string();
}

// Write a complete input record for a converter regression
std::string InputText(const std::string& name, const std::string& text) {
  const auto    path = ProgramPath(name);
  std::ofstream output(path);
  output << text;
  output.close();
  REQUIRE_FALSE(output.fail());
  return path;
}

}  // namespace

// Verify exact event totals, positive populations and component rounding errors
TEST_CASE("Minimum bias allocation preserves the requested count", "[program][minbias]") {
  for (const int events : {0, 1, 2, 10, 10003, std::numeric_limits<int>::max()}) {
    for (const auto& xs : {std::array<double, 3>{0.25, 0.25, 0.5}, std::array<double, 3>{0.0, 0.0, 1.0},
                           std::array<double, 3>{1e-15, 0.9, 0.1}}) {
      const auto counts = gra::program::MinbiasCounts(events, xs);
      CHECK(std::accumulate(counts.begin(), counts.end(), 0LL) == events);
      const long double total = std::accumulate(xs.begin(), xs.end(), 0.0L);
      for (const auto& i : indices(counts)) {
        CHECK(counts[i] >= 0);
        CHECK(std::abs(counts[i] - events * static_cast<long double>(xs[i]) / total) <= 1.0L);
        if (!(xs[i] > 0.0)) { CHECK(counts[i] == 0); }
      }
    }
  }
  CHECK_THROWS_AS(gra::program::MinbiasCounts(1, {0, 0, 0}), std::invalid_argument);
  CHECK_THROWS_AS(gra::program::MinbiasCounts(1, {1, 1, -1}), std::invalid_argument);
  CHECK_THROWS_AS(gra::program::MinbiasCounts(-1, {1, 1, 1}), std::invalid_argument);
}

// Compare the first CNAB2 step to the analytic two Fourier mode expression
TEST_CASE("KS starts from the physical field for both FFT implementations", "[program][pde]") {
  const int                           points = 32;
  const int                           length = 16;
  const double                        k      = 2.0 * gra::math::PI / length;
  const double                        dt     = 0.025;
  std::valarray<std::complex<double>> initial(points);
  for (const auto& i : indices(initial)) { initial[i] = std::cos(2.0 * gra::math::PI * i / points); }
  const double l1 = k * k - gra::math::pow4(k);
  const double l2 = 4 * k * k - 16 * gra::math::pow4(k);
  const double a1 = (1 + dt * l1 / 2) / (1 - dt * l1 / 2);
  const double a2 = dt * k / (2 * (1 - dt * l2 / 2));
  for (const bool eigen : {false, true}) {
    gra::program::KS evolution(initial, length, dt, eigen);
    evolution.Step();
    const auto field = evolution.Field();
    for (const auto& i : indices(field)) {
      const double phase = 2.0 * gra::math::PI * i / points;
      CHECK(field[i].real() == Approx(a1 * std::cos(phase) + a2 * std::sin(2 * phase)).margin(1e-12));
      CHECK(std::abs(field[i].imag()) < 1e-12);
    }
    for (int i = 0; i < 20; ++i) { evolution.Step(); }
    CHECK(std::abs(evolution.Field().sum()) < 1e-10);
  }
}

// Reject invalid simulation inputs before allocation, indexing or division
TEST_CASE("Benchmarks validate lattice and time inputs", "[program][input]") {
  gra::program::PathParam path;
  REQUIRE_NOTHROW(path.Validate(30));
  for (const int slit : {-1, 0, 1, 3, 4}) {
    path.k_slit = slit;
    CHECK_THROWS_AS(path.Validate(30), std::invalid_argument);
  }
  path    = gra::program::PathParam{};
  path.dt = 0.0;
  CHECK_THROWS_AS(path.Validate(30), std::invalid_argument);
  CHECK_THROWS_AS(gra::program::ValidateKS(64, 128, 0.1, 1, 0, false), std::invalid_argument);
  CHECK_THROWS_AS(gra::program::ValidateKS(64, 0, 0.1, 1, 1, false), std::invalid_argument);
  CHECK_THROWS_AS(gra::program::ValidateKS(64, 6, 0.1, 1, 1, false), std::invalid_argument);
  CHECK_NOTHROW(gra::program::ValidateKS(64, 6, 0.1, 1, 1, true));
  CHECK_THROWS_AS(gra::program::Number<int>("1.5"), std::invalid_argument);
  CHECK_THROWS_AS(gra::program::Number<unsigned long>("-1"), std::invalid_argument);
  CHECK_THROWS_AS(gra::program::Number<double>("nan"), std::invalid_argument);
  const auto energies = gra::program::Energies("500, 2760.5 ,1e4");
  REQUIRE(energies.size() == 3);
  CHECK(energies[0] == Approx(500));
  CHECK(energies[1] == Approx(2760.5));
  CHECK(energies[2] == Approx(10000));
  for (const auto& text : {"", "500,", ",500", "500,,700", "500x", "nan", "inf", "0", "-1"}) {
    CHECK_THROWS_AS(gra::program::Energies(text), std::invalid_argument);
  }
}

// Compare the exact diffraction kernel with a numerical Green function derivative
TEST_CASE("Rayleigh Sommerfeld retains the near aperture Green function term", "[program][diffraction]") {
  const gra::M4Vec aperture(0.02, -0.01, 0, 0.3);
  const gra::M4Vec source(0, 0, -1000, 0);
  const gra::M4Vec normal(0, 0, 1, 0);
  const double     k = 20.0;
  for (const double z : {0.003, 0.1, 5.0}) {
    const gra::M4Vec detector(0, 0, z, 0.3);
    const double     step = 1e-7;
    const gra::M4Vec shift(0, 0, step, 0);
    // Evaluate the Green function independently on the two displaced aperture planes
    const auto green = [&](const gra::M4Vec& point) {
      const double r = (point - detector).P3mod();
      return std::exp(-gra::math::zi * k * r) / r;
    };
    const auto derivative = (green(aperture + shift) - green(aperture - shift)) / (2 * step);
    const auto expected   = gra::program::u_point(aperture, source, k) * derivative / (2 * gra::math::PI);
    const auto actual     = gra::program::RS_integrand(aperture, detector, source, normal, k);
    CHECK(std::abs(actual - expected) < 1e-7 * std::abs(expected));
  }
}

// Preserve pion four momenta, PID weights and event identities during conversion
TEST_CASE("CSV conversion preserves weighted pion events", "[program][conversion]") {
  const auto input  = InputText("pions.csv", "0.3,0,0,0,-0.3,0,0,0,0.25\n0.4,0,0,0,-0.4,0,0,0,2.5\n");
  const auto output = ProgramPath("csv/nested/pions.hepmc3");
  REQUIRE(gra::program::ConvertData(input, output, 1e-9) == 2);
  gra::MHepMCReader reader(output);
  HepMC3::GenEvent  event;
  for (int i = 0; i < 2; ++i) {
    REQUIRE(reader.Read(event));
    CHECK(event.event_number() == i);
    REQUIRE(event.weights().size() == 1);
    CHECK(event.weights()[0] == Approx(i == 0 ? 0.25 : 2.5));
    CHECK(event.cross_section()->xsec() == Approx(1000.0));
    HepMC3::FourVector pair;
    for (const auto& particle : event.particles()) {
      if (particle->status() == 1) { pair = pair + particle->momentum(); }
    }
    const double pt = i == 0 ? 0.3 : 0.4;
    CHECK(pair.m() == Approx(2 * std::sqrt(pt * pt + gra::PDG::mpi * gra::PDG::mpi)));
  }
  CHECK_FALSE(reader.Read(event));
  const auto broken = InputText("broken.csv", "1,2,3\n");
  CHECK_THROWS(gra::program::ConvertData(broken, ProgramPath("broken.hepmc3"), 1));
  CHECK_THROWS(gra::program::ConvertData(ProgramPath("absent.csv"), output, 1));
  CHECK_THROWS(gra::program::ConvertData(input, input + "/result", 1));
  if (std::filesystem::exists("/dev/full")) { CHECK_THROWS(gra::program::ConvertData(input, "/dev/full", 1)); }
  const std::array<double, 9> row{0.3, 0, 0, 0, -0.3, 0, 0, 0, 1};
  CHECK_THROWS_AS(gra::program::PionData(row, 0, std::numeric_limits<double>::max()), std::invalid_argument);
  auto large = row;
  large[0] = std::numeric_limits<double>::max();
  large[1] = std::numeric_limits<double>::max();
  CHECK_THROWS_AS(gra::program::PionData(large, 0, 1), std::invalid_argument);
}

// Preserve source bytes when converter output aliases the input through filesystem links
TEST_CASE("Converters reject output aliases before truncation", "[program][conversion]") {
  const std::string text = "0.3,0,0,0,-0.3,0,0,0,1\n";
  const auto input = InputText("alias_input.csv", text);
  const auto hard = ProgramPath("alias_hard.csv");
  const auto soft = ProgramPath("alias_soft.csv");
  if (!std::filesystem::exists(hard)) { std::filesystem::create_hard_link(input, hard); }
  if (!std::filesystem::exists(soft)) { std::filesystem::create_symlink(input, soft); }
  for (const auto& output : {input, hard, soft}) {
    CHECK_THROWS_AS(gra::program::ConvertData(input, output, 1), std::invalid_argument);
    CHECK_THROWS_AS(gra::program::ConvertHEPEVT(input, output), std::invalid_argument);
    std::ifstream source(input);
    const std::string content((std::istreambuf_iterator<char>(source)), std::istreambuf_iterator<char>());
    CHECK(content == text);
  }
}

// Convert valid HEPEVT and reject missing, malformed and truncated input
TEST_CASE("HEPEVT conversion distinguishes EOF from damaged input", "[program][conversion]") {
  const auto input = ProgramPath("valid.hepevt");
  {
    HepMC3::WriterHEPEVT writer(input);
    auto                 event = analysis_test::Event();
    writer.write_event(event);
    writer.close();
    REQUIRE_FALSE(writer.failed());
  }
  const auto output = ProgramPath("hepevt/nested/events.hepmc3");
  CHECK(gra::program::ConvertHEPEVT(input, output) == 1);
  CHECK_THROWS(gra::program::ConvertHEPEVT(ProgramPath("absent.hepevt"), output));
  const auto invalid = InputText("invalid.hepevt", "broken event\n");
  CHECK_THROWS(gra::program::ConvertHEPEVT(invalid, output));
  std::ifstream     source(input);
  const std::string text((std::istreambuf_iterator<char>(source)), std::istreambuf_iterator<char>());
  const auto no_newline = InputText("no_newline.hepevt", text.substr(0, text.find_last_not_of("\r\n") + 1));
  CHECK(gra::program::ConvertHEPEVT(no_newline, output) == 1);
  const auto        truncated = InputText("truncated.hepevt", text.substr(0, text.size() / 2));
  CHECK_THROWS(gra::program::ConvertHEPEVT(truncated, output));
  if (std::filesystem::exists("/dev/full")) { CHECK_THROWS(gra::program::ConvertHEPEVT(input, "/dev/full")); }
}

// Compare weighted transport histograms across HepMC momentum units
TEST_CASE("Transport mass probabilities preserve weights and GeV units", "[program][transport]") {
  auto first  = gra::program::PionData({0.3, 0, 0, 0, -0.3, 0, 0, 0, 1.0}, 0, 1);
  auto second = gra::program::PionData({0.6, 0, 0, 0, -0.6, 0, 0, 0, 3.0}, 1, 1);
  for (const auto unit : {HepMC3::Units::GEV, HepMC3::Units::MEV}) {
    const auto       file = analysis_test::Write("transport_" + std::to_string(unit), {first, second}, unit);
    gra::MH1<double> histogram(2, 0, 2, "mass");
    gra::program::ReadMass(file, histogram);
    const auto probability = histogram.GetProbDensity();
    CHECK(probability[0] == Approx(0.25));
    CHECK(probability[1] == Approx(0.75));
  }
  first.weights()[0]        = -1;
  const auto       negative = analysis_test::Write("negative_transport", {first});
  gra::MH1<double> histogram(2, 0, 2, "mass");
  CHECK_THROWS_AS(gra::program::ReadMass(negative, histogram), std::invalid_argument);
  const auto outside = gra::program::PionData({3, 0, 0, 0, -3, 0, 0, 0, 1}, 0, 1);
  const auto empty_range = analysis_test::Write("outside_transport", {outside});
  histogram.Fill(1.0, 1.0);
  CHECK_THROWS_AS(gra::program::ReadMass(empty_range, histogram), std::invalid_argument);
}

// Reject type and range errors before integer conversion and check output flush
TEST_CASE("Harmonic counts and measurement output reject invalid values", "[program][harmonic]") {
  using nlohmann::json;
  for (const auto& invalid : {json(-1), json(1.5), json(true), json("2"), json(nullptr)}) {
    CHECK_THROWS_AS(gra::program::Count<std::size_t>(invalid, "bins"), std::invalid_argument);
  }
  CHECK(gra::program::Count<std::size_t>(json(2), "bins") == 2);
  CHECK(gra::program::Count<std::uint64_t>(json(std::numeric_limits<std::uint64_t>::max()), "seed") ==
        std::numeric_limits<std::uint64_t>::max());
  CHECK_THROWS_AS(gra::program::Count<int>(json(std::numeric_limits<std::uint64_t>::max()), "lmax"),
                  std::invalid_argument);
  const auto output      = ProgramPath("harmonic/nested/measurement.json");
  const json measurement = {{"coefficients", {0.125, -0.75}}};
  gra::program::WriteOutput(output, measurement);
  std::ifstream input(output);
  json          read;
  input >> read;
  CHECK(read == measurement);
  CHECK_THROWS(gra::program::WriteOutput(output + "/result.json", measurement));
  if (std::filesystem::exists("/dev/full")) { CHECK_THROWS(gra::program::WriteOutput("/dev/full", measurement)); }
}

// Avoid assigning a physical zero to an unavailable eikonal cross section
TEST_CASE("Energy scan marks unavailable soft cross sections", "[program][xscan]") {
  const gra::MEikonal eikonal;
  const auto          xs = gra::program::ScanXS(eikonal);
  for (const double value : xs) { CHECK(std::isnan(value)); }
}
