// Detector aware harmonic measurement unit tests
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// Catch
#include <catch.hpp>

// C++
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>
#include <vector>

// Own
#include "Graniitti/Analysis/MHarmonic.h"
#include "Graniitti/Analysis/MHarmonicHepMC.h"
#include "Graniitti/Analysis/MHepMCReader.h"
#include "Graniitti/Analysis/MSpherical.h"
#include "support/analysis_test_support.hh"
#include "Graniitti/Math/MIntegration.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MSpecialFunctions.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Math/MStatistics.h"

using gra::aux::indices;

namespace {

using gra::harmonic::Axis;
using gra::harmonic::AngularFrame;
using gra::harmonic::Cell;
using gra::harmonic::Coordinate;
using gra::harmonic::DataEvent;
using gra::harmonic::EventKinematics;
using gra::harmonic::FiducialDecision;
using gra::harmonic::FiducialRange;
using gra::harmonic::FitConfig;
using gra::harmonic::HarmonicEstimator;
using gra::harmonic::HarmonicMeasurement;
using gra::harmonic::MeasurementMode;
using gra::harmonic::MeasurementResult;
using gra::harmonic::Observation;
using gra::harmonic::PhaseSpaceGrid;
using gra::harmonic::ResponseEvent;
using gra::harmonic::ToyResponse;
using gra::harmonic::ToyResponseConfig;

// Build one valid angular observation at selected coordinates
Observation MakeObservation(const std::vector<double> &z, double weight = 1.0) {
  return Observation{0.0, 0.0, z, weight};
}

// Build one valid observation with explicit angular coordinates
Observation MakeAngularObservation(const std::vector<double> &z,
                                   double costheta, double phi,
                                   double weight = 1.0) {
  return Observation{costheta, phi, z, weight};
}

// Build one generated and optionally reconstructed response event
ResponseEvent MakeResponse(std::uint64_t key,
                           const std::vector<double> &truth_z,
                           const std::optional<std::vector<double>> &reco_z,
                           FiducialDecision fiducial = {true, true}) {
  ResponseEvent event;
  event.event_key = key;
  event.source = "unit";
  event.truth = MakeObservation(truth_z);
  event.fiducial = fiducial;
  if (reco_z.has_value()) {
    event.reco = MakeObservation(*reco_z);
    event.reco_selected = true;
  }
  return event;
}

// Build an exclusive on-shell pp -> pp pi+ pi- state at 13 TeV
EventKinematics MakeKinematics() {
  EventKinematics event;
  event.pip.SetPxPyPzM(0.4, 0.1, 0.2, gra::PDG::mpi);
  event.pim.SetPxPyPzM(-0.4, -0.1, -0.2, gra::PDG::mpi);
  const double beam_pz = std::sqrt(6500.0 * 6500.0 - gra::PDG::mp * gra::PDG::mp);
  event.beam_plus.SetPxPyPzM(0.0, 0.0, beam_pz, gra::PDG::mp);
  event.beam_minus.SetPxPyPzM(0.0, 0.0, -beam_pz, gra::PDG::mp);
  const double energy = 6500.0 - event.pip.E();
  const double proton_pz = std::sqrt(energy * energy - gra::PDG::mp * gra::PDG::mp - 0.05);
  event.proton_plus.SetPxPyPzM(0.1, 0.2, proton_pz, gra::PDG::mp);
  event.proton_minus.SetPxPyPzM(-0.1, -0.2, -proton_pz, gra::PDG::mp);
  event.has_forward_protons = true;
  return event;
}

// Encode the exclusive state as one HepMC scattering vertex
HepMC3::GenEvent KinematicsEvent(const EventKinematics &state) {
  HepMC3::GenEvent event;
  auto vertex = std::make_shared<HepMC3::GenVertex>();
  vertex->add_particle_in(std::make_shared<HepMC3::GenParticle>(
      gra::aux::M4Vec2HepMC3(state.beam_plus), gra::PDG::PDG_p, gra::PDG::PDG_BEAM));
  vertex->add_particle_in(std::make_shared<HepMC3::GenParticle>(
      gra::aux::M4Vec2HepMC3(state.beam_minus), gra::PDG::PDG_p, gra::PDG::PDG_BEAM));
  vertex->add_particle_out(std::make_shared<HepMC3::GenParticle>(
      gra::aux::M4Vec2HepMC3(state.pip), gra::PDG::PDG_pip, gra::PDG::PDG_STABLE));
  vertex->add_particle_out(std::make_shared<HepMC3::GenParticle>(
      gra::aux::M4Vec2HepMC3(state.pim), gra::PDG::PDG_pim, gra::PDG::PDG_STABLE));
  vertex->add_particle_out(std::make_shared<HepMC3::GenParticle>(
      gra::aux::M4Vec2HepMC3(state.proton_plus), gra::PDG::PDG_p, gra::PDG::PDG_STABLE));
  vertex->add_particle_out(std::make_shared<HepMC3::GenParticle>(
      gra::aux::M4Vec2HepMC3(state.proton_minus), gra::PDG::PDG_p, gra::PDG::PDG_STABLE));
  event.add_vertex(vertex);
  event.weights() = {state.weight};
  return event;
}

// Build a zeroth-order fit configuration without response replicas
FitConfig ScalarFitConfig() {
  FitConfig config;
  config.lmax = 0;
  config.remove_odd = false;
  config.remove_negative_m = false;
  config.svd_relative_cut = 0.0;
  config.response_jackknife_bins = 0;
  config.max_parameters = 20;
  return config;
}

} // namespace

// Check inclusive fiducial boundaries and reject invalid numeric intervals
TEST_CASE("Harmonic fiducial ranges enforce finite ordered boundaries",
          "[harmonic][fiducial][config]") {
  const FiducialRange range{-2.4, 2.4};
  REQUIRE(range.Contains(-2.4));
  REQUIRE(range.Contains(2.4));
  REQUIRE_FALSE(range.Contains(2.4001));
  REQUIRE_FALSE(range.Contains(std::numeric_limits<double>::quiet_NaN()));
  REQUIRE_NOTHROW(range.Validate("pion eta interval"));
  REQUIRE_THROWS_AS((FiducialRange{1.0, 1.0}.Validate("empty interval")),
                    std::invalid_argument);
  REQUIRE_THROWS_AS((FiducialRange{2.0, 1.0}.Validate("reversed interval")),
                    std::invalid_argument);
}

// Check central and tagged HepMC configurations against their physical modes
TEST_CASE("Harmonic HepMC configuration separates central and tagged cuts",
          "[harmonic][hepmc][config]") {
  gra::harmonic::HepMCReadConfig config;
  config.mode = MeasurementMode::Central;
  config.frame = AngularFrame::CS;
  config.cuts.pion_eta = {-2.4, 2.4};
  config.cuts.pion_pt = {0.2, 100.0};
  config.sqrt_s = 13000.0;
  config.max_records = 100;
  REQUIRE_NOTHROW(config.Validate());

  config.frame = AngularFrame::GJ;
  REQUIRE_THROWS_AS(config.Validate(), std::invalid_argument);

  config.mode = MeasurementMode::Tagged;
  config.cuts.proton_xi = {0.01, 0.2};
  config.cuts.proton_abs_t = {0.03, 4.0};
  config.selection = {0.1, 0.1, 0.1};
  REQUIRE_NOTHROW(config.Validate());

  config.cuts.proton_xi.max = 1.01;
  REQUIRE_THROWS_AS(config.Validate(), std::invalid_argument);
}

// Check every explicit angular frame and reject aliases with changed case
TEST_CASE("Harmonic angular frame parsing is explicit",
          "[harmonic][hepmc][frame]") {
  REQUIRE(gra::harmonic::ParseAngularFrame("CM") == AngularFrame::CM);
  REQUIRE(gra::harmonic::ParseAngularFrame("HX") == AngularFrame::HX);
  REQUIRE(gra::harmonic::ParseAngularFrame("CS") == AngularFrame::CS);
  REQUIRE(gra::harmonic::ParseAngularFrame("AH") == AngularFrame::AH);
  REQUIRE(gra::harmonic::ParseAngularFrame("PG") == AngularFrame::PG);
  REQUIRE(gra::harmonic::ParseAngularFrame("GJ") == AngularFrame::GJ);
  REQUIRE_THROWS_AS(gra::harmonic::ParseAngularFrame("cs"),
                    std::invalid_argument);
}

// Check mass shells, beam arms and cylindrical covariance of the detector model
TEST_CASE("Harmonic toy response reconstructs physical particle momenta",
          "[harmonic][response][kinematics]") {
  ToyResponseConfig config;
  config.pion_efficiency = 0.8;
  config.proton_efficiency = 0.7;
  config.pion_logpt_sigma = 0.1;
  config.pion_eta_sigma = 0.2;
  config.pion_phi_sigma = 0.3;
  config.proton_pt_sigma = 0.1;
  config.proton_logpz_sigma = 0.0001;
  config.seed = 17;
  const ToyResponse response(config);
  EventKinematics truth = MakeKinematics();
  REQUIRE(response.Efficiency(truth) == Approx(0.8 * 0.8));
  truth.mode = MeasurementMode::Tagged;
  REQUIRE(response.Efficiency(truth) == Approx(0.8 * 0.8 * 0.7 * 0.7));
  const EventKinematics first = response.Reconstruct(truth, 42);
  const EventKinematics second = response.Reconstruct(truth, 42);
  for (std::size_t mu = 0; mu < 4; ++mu) {
    REQUIRE(first.pip[mu] == Approx(second.pip[mu]));
    REQUIRE(first.pim[mu] == Approx(second.pim[mu]));
    REQUIRE(first.proton_plus[mu] == Approx(second.proton_plus[mu]));
    REQUIRE(first.proton_minus[mu] == Approx(second.proton_minus[mu]));
  }
  REQUIRE(first.pip.M2() == Approx(truth.pip.M2()).margin(1e-12));
  REQUIRE(first.pim.M2() == Approx(truth.pim.M2()).margin(1e-12));
  REQUIRE(first.proton_plus.M2() == Approx(truth.proton_plus.M2()).margin(2e-8));
  REQUIRE(first.proton_minus.M2() == Approx(truth.proton_minus.M2()).margin(2e-8));
  REQUIRE((first.pip + first.pim).M() >= 2.0 * gra::PDG::mpi);
  REQUIRE(first.proton_plus.Pz() > 0.0);
  REQUIRE(first.proton_minus.Pz() < 0.0);
  REQUIRE(first.pip.Pt() > 0.0);

  EventKinematics rotated = truth;
  rotated.pip.RotateZ(0.7);
  rotated.pim.RotateZ(0.7);
  const EventKinematics reco_rotated = response.Reconstruct(rotated, 42);
  gra::M4Vec expected = first.pip;
  expected.RotateZ(0.7);
  for (std::size_t mu = 0; mu < 4; ++mu) {
    REQUIRE(reco_rotated.pip[mu] == Approx(expected[mu]).margin(1e-12));
  }
  const EventKinematics unchanged = ToyResponse(ToyResponseConfig{}).Reconstruct(truth, 42);
  for (std::size_t mu = 0; mu < 4; ++mu) {
    REQUIRE(unchanged.pip[mu] == Approx(truth.pip[mu]).margin(1e-12));
    REQUIRE(unchanged.proton_plus[mu] == Approx(truth.proton_plus[mu]).margin(1e-12));
  }
  config.pion_efficiency = 0.0;
  REQUIRE_FALSE(gra::harmonic::SimulateDetector(ToyResponse(config), truth, 42, 17).has_value());
}

// Check the shared inverse used by both harmonic estimators
TEST_CASE("Harmonic pseudoinverse retains rank and cutoff diagnostics",
          "[harmonic][inverse]") {
  gra::MMatrix<double> response{{1.0, 0.0}, {0.0, 1.0e-4}};
  gra::PseudoInverseDiagnostics diagnostics;
  const gra::MMatrix<double> regularized =
      response.PseudoInverse(1.0e-3, &diagnostics);
  REQUIRE(diagnostics.numerical_rank == 2);
  REQUIRE(diagnostics.retained_rank == 1);
  REQUIRE(regularized(0, 0) == Approx(1.0));
  REQUIRE(regularized(1, 1) == Approx(0.0));

  response(1, 1) = 0.0;
  const gra::MMatrix<double> inverse = response.PseudoInverse(0.0);
  REQUIRE(inverse(0, 0) == Approx(1.0));
  REQUIRE(inverse(1, 1) == Approx(0.0));
}

// Check that the two modes enforce their physically distinct coordinates
TEST_CASE("Harmonic measurement modes have canonical coordinates",
          "[harmonic][coordinates]") {
  const PhaseSpaceGrid central(MeasurementMode::Central,
                               {{Coordinate::Mass, 2, 1.0, 3.0},
                                {Coordinate::Momentum, 1, 0.0, 1.0},
                                {Coordinate::Rapidity, 1, -1.0, 1.0}});
  REQUIRE((central.Locate({2.0, 0.5, 0.0}) == Cell{1, 0, 0}));
  REQUIRE_FALSE(central.Locate({3.1, 0.5, 0.0}).has_value());

  const PhaseSpaceGrid tagged(
      MeasurementMode::Tagged,
      {{Coordinate::Mass, 1, 1.0, 3.0},
       {Coordinate::Rapidity, 1, -1.0, 1.0},
       {Coordinate::AbsT1, 1, 0.0, 1.0},
       {Coordinate::AbsT2, 1, 0.0, 1.0},
       {Coordinate::DeltaPhiPP, 1, -gra::math::PI, gra::math::PI}});
  REQUIRE((tagged.Locate({2.0, 0.0, 0.2, 0.3, 0.1}) == Cell{0, 0, 0, 0, 0}));

  REQUIRE((FiducialDecision{true, false}.Pass(MeasurementMode::Central)));
  REQUIRE_FALSE((FiducialDecision{true, false}.Pass(MeasurementMode::Tagged)));
}

// Check global inversion of migrations between conditional cells
TEST_CASE("Harmonic response unfolds conditional cell migrations",
          "[harmonic][migration]") {
  const PhaseSpaceGrid grid(MeasurementMode::Central,
                            {{Coordinate::Mass, 2, 1.0, 3.0},
                             {Coordinate::Momentum, 1, 0.0, 1.0},
                             {Coordinate::Rapidity, 1, -1.0, 1.0}});
  const std::vector<double> truth0 = {1.5, 0.5, 0.0};
  const std::vector<double> truth1 = {2.5, 0.5, 0.0};
  std::vector<ResponseEvent> response;
  for (std::uint64_t index = 0; index < 100; ++index) {
    response.push_back(
        MakeResponse(index, truth0, index < 80 ? truth0 : truth1));
  }
  for (std::uint64_t index = 0; index < 100; ++index) {
    response.push_back(
        MakeResponse(100 + index, truth1, index < 10 ? truth0 : truth1));
  }

  std::vector<DataEvent> data;
  for (std::uint64_t index = 0; index < 850; ++index) {
    data.push_back({index, "data", MakeObservation(truth0)});
  }
  for (std::uint64_t index = 0; index < 650; ++index) {
    data.push_back({850 + index, "data", MakeObservation(truth1)});
  }

  const HarmonicMeasurement measurement(grid, ScalarFitConfig());
  const auto result = measurement.Fit(response, data);
  REQUIRE(result.flat.at(Cell{0, 0, 0}).value[0] ==
          Approx(1000.0).epsilon(1e-10));
  REQUIRE(result.flat.at(Cell{1, 0, 0}).value[0] ==
          Approx(500.0).epsilon(1e-10));
  REQUIRE(result.fiducial.at(Cell{0, 0, 0}).value[0] ==
          Approx(1000.0).epsilon(1e-10));
  REQUIRE(result.fiducial.at(Cell{1, 0, 0}).value[0] ==
          Approx(500.0).epsilon(1e-10));

  FitConfig eml_config = ScalarFitConfig();
  eml_config.estimator = HarmonicEstimator::EML;
  eml_config.eml_positivity_costheta = 4;
  eml_config.eml_positivity_phi = 8;
  const MeasurementResult eml =
      HarmonicMeasurement(grid, eml_config).Fit(response, data);
  REQUIRE(eml.flat.at(Cell{0, 0, 0}).value[0] == Approx(1000.0).epsilon(1e-6));
  REQUIRE(eml.flat.at(Cell{1, 0, 0}).value[0] == Approx(500.0).epsilon(1e-6));
}

// Check that tagged forward cuts belong to F while central-only cuts do not
TEST_CASE("Tagged fiducial response requires both proton arms",
          "[harmonic][fiducial]") {
  const PhaseSpaceGrid grid(
      MeasurementMode::Tagged,
      {{Coordinate::Mass, 1, 1.0, 3.0},
       {Coordinate::Rapidity, 1, -1.0, 1.0},
       {Coordinate::AbsT1, 1, 0.0, 1.0},
       {Coordinate::AbsT2, 1, 0.0, 1.0},
       {Coordinate::DeltaPhiPP, 1, -gra::math::PI, gra::math::PI}});
  const std::vector<double> z = {2.0, 0.0, 0.2, 0.3, 0.1};
  std::vector<ResponseEvent> response;
  for (std::uint64_t index = 0; index < 100; ++index) {
    response.push_back(
        MakeResponse(index, z, z, FiducialDecision{true, index < 50}));
  }
  std::vector<DataEvent> data;
  for (std::uint64_t index = 0; index < 100; ++index) {
    data.push_back({index, "data", MakeObservation(z)});
  }

  const HarmonicMeasurement measurement(grid, ScalarFitConfig());
  MeasurementResult result = measurement.Fit(response, data);
  const Cell cell = {0, 0, 0, 0, 0};
  REQUIRE(result.flat.at(cell).value[0] == Approx(100.0).epsilon(1e-10));
  REQUIRE(result.fiducial.at(cell).value[0] == Approx(50.0).epsilon(1e-10));
  REQUIRE(result.detector.at(cell).value[0] == Approx(100.0).epsilon(1e-10));

  gra::harmonic::ApplyNormalization(result, 10.0, 0.1);
  REQUIRE(result.fiducial.at(cell).value[0] == Approx(5.0).epsilon(1e-10));
  REQUIRE(result.fiducial.at(cell).covariance[0][0] ==
          Approx(0.5).epsilon(1e-10));
  REQUIRE(result.fiducial_covariance[0][0] == Approx(0.5).epsilon(1e-10));
}

// Check the real harmonic basis normalization and response orientation
TEST_CASE("Harmonic response recovers a known angular coefficient",
          "[harmonic][basis]") {
  const PhaseSpaceGrid grid(MeasurementMode::Central,
                            {{Coordinate::Mass, 1, 1.0, 3.0},
                             {Coordinate::Momentum, 1, 0.0, 1.0},
                             {Coordinate::Rapidity, 1, -1.0, 1.0}});
  const std::vector<double> z = {2.0, 0.5, 0.0};
  constexpr std::size_t ncostheta = 40;
  constexpr std::size_t nphi = 12;
  constexpr double coefficient = 0.2;
  std::vector<ResponseEvent> response;
  std::vector<DataEvent> data;
  std::uint64_t key = 0;
  for (std::size_t i = 0; i < ncostheta; ++i) {
    const double costheta = -1.0 + (2.0 * static_cast<double>(i) + 1.0) /
                                       static_cast<double>(ncostheta);
    for (std::size_t j = 0; j < nphi; ++j) {
      const double phi = -gra::math::PI + 2.0 * gra::math::PI *
                                              (static_cast<double>(j) + 0.5) /
                                              static_cast<double>(nphi);
      const Observation truth = MakeAngularObservation(z, costheta, phi);
      ResponseEvent event;
      event.event_key = key;
      event.source = "unit";
      event.truth = truth;
      event.fiducial = {true, true};
      event.reco = truth;
      event.reco_selected = true;
      response.push_back(event);

      const double angular_weight =
          1.0 + coefficient * std::sqrt(3.0) * costheta;
      data.push_back(
          {key, "data",
           MakeAngularObservation(z, costheta, phi, angular_weight)});
      ++key;
    }
  }

  FitConfig config = ScalarFitConfig();
  config.lmax = 1;
  const HarmonicMeasurement measurement(grid, config);
  const MeasurementResult result = measurement.Fit(response, data);
  const Cell cell = {0, 0, 0};
  const double events = static_cast<double>(ncostheta * nphi);
  REQUIRE(result.flat.at(cell).value[0] == Approx(events).epsilon(1e-10));
  REQUIRE(result.flat.at(cell).value[1] == Approx(0.0).margin(1e-10));
  REQUIRE(result.flat.at(cell).value[2] ==
          Approx(events * coefficient).epsilon(1e-10));
  REQUIRE(result.flat.at(cell).value[3] == Approx(0.0).margin(1e-10));

  config.estimator = HarmonicEstimator::EML;
  config.eml_positivity_costheta = 8;
  config.eml_positivity_phi = 12;
  const MeasurementResult eml =
      HarmonicMeasurement(grid, config).Fit(response, data);
  REQUIRE(eml.flat.at(cell).value[0] == Approx(events).epsilon(2e-3));
  REQUIRE(eml.flat.at(cell).value[2] ==
          Approx(events * coefficient).epsilon(2e-3));
}

// Check that asymmetric acceptance generates moments forbidden in production
TEST_CASE("Production symmetries preserve acceptance induced moments",
          "[harmonic][basis][fiducial]") {
  const PhaseSpaceGrid grid(MeasurementMode::Central,
                           {{Coordinate::Mass, 1, 1.0, 3.0},
                            {Coordinate::Momentum, 1, 0.0, 1.0},
                            {Coordinate::Rapidity, 1, -1.0, 1.0}});
  FitConfig config = ScalarFitConfig();
  config.lmax = 1;
  std::size_t moment = 2;
  SECTION("odd l production symmetry") { config.remove_odd = true; }
  SECTION("negative m production symmetry") {
    config.remove_negative_m = true;
    moment = 1;
  }
  std::vector<ResponseEvent> response;
  std::vector<DataEvent> data;
  double selected_moment = 0.0;
  std::uint64_t key = 0;
  for (std::size_t i = 0; i < 40; ++i) {
    const double costheta = -1.0 + (2.0 * static_cast<double>(i) + 1.0) / 40.0;
    for (std::size_t j = 0; j < 12; ++j) {
      const double phi = -gra::math::PI + 2.0 * gra::math::PI *
                        (static_cast<double>(j) + 0.5) / 12.0;
      const double direction = moment == 2 ? costheta :
          std::sqrt(1.0 - costheta * costheta) * std::sin(phi);
      const bool accepted = direction > 0.0;
      ResponseEvent event = MakeResponse(key, {2.0, 0.5, 0.0}, std::nullopt,
                                        {accepted, true});
      event.truth.costheta = costheta;
      event.truth.phi = phi;
      if (accepted) {
        event.reco = event.truth;
        event.reco_selected = true;
        data.push_back({key, "data", event.truth});
        selected_moment += std::sqrt(3.0) * direction;
      }
      response.push_back(event);
      ++key;
    }
  }
  const MeasurementResult result = HarmonicMeasurement(grid, config).Fit(response, data);
  const Cell cell = {0, 0, 0};
  REQUIRE(result.flat.at(cell).value[0] == Approx(480.0).epsilon(1e-10));
  REQUIRE(result.flat.at(cell).value[moment] == Approx(0.0).margin(1e-10));
  REQUIRE(result.fiducial.at(cell).value[0] == Approx(240.0).epsilon(1e-10));
  REQUIRE(result.fiducial.at(cell).value[moment] == Approx(selected_moment).epsilon(1e-10));
  REQUIRE(result.detector.at(cell).value[moment] == Approx(selected_moment).epsilon(1e-10));
  REQUIRE(result.response_rows == 4);
  REQUIRE(result.response_columns == (moment == 2 ? 1 : 3));
  REQUIRE(result.detector.at(cell).covariance[moment][moment] > 0.0);
}

// Check rectangular EML response blocks and the covariance of fiducial dipoles
TEST_CASE("Extended likelihood retains acceptance asymmetry in each cell",
          "[harmonic][eml][fiducial][normalization]") {
  const PhaseSpaceGrid grid(MeasurementMode::Central,
                           {{Coordinate::Mass, 2, 1.0, 3.0},
                            {Coordinate::Momentum, 1, 0.0, 1.0},
                            {Coordinate::Rapidity, 1, -1.0, 1.0}});
  std::vector<ResponseEvent> response;
  std::vector<DataEvent> data;
  std::vector<double> dipole(2, 0.0);
  std::uint64_t key = 0;
  for (std::size_t cell = 0; cell < 2; ++cell) {
    const std::vector<double> z = {1.5 + static_cast<double>(cell), 0.5, 0.0};
    const double slope = cell == 0 ? 0.4 : -0.3;
    for (std::size_t i = 0; i < 24; ++i) {
      const double costheta = -1.0 + (2.0 * static_cast<double>(i) + 1.0) / 24.0;
      const double efficiency = 0.5 * (1.0 + slope * costheta);
      for (std::size_t j = 0; j < 12; ++j) {
        const double phi = -gra::math::PI + 2.0 * gra::math::PI *
                          (static_cast<double>(j) + 0.5) / 12.0;
        ResponseEvent accepted = MakeResponse(key, z, z);
        accepted.truth = MakeAngularObservation(z, costheta, phi, efficiency);
        accepted.reco = accepted.truth;
        response.push_back(accepted);
        data.push_back({key, "data", accepted.truth});
        ResponseEvent rejected = accepted;
        rejected.truth.weight = 1.0 - efficiency;
        rejected.reco.reset();
        rejected.reco_selected = false;
        rejected.fiducial = {false, false};
        response.push_back(rejected);
        dipole[cell] += efficiency * std::sqrt(3.0) * costheta;
        ++key;
      }
    }
  }
  FitConfig config = ScalarFitConfig();
  config.lmax = 1;
  config.remove_odd = true;
  config.estimator = HarmonicEstimator::EML;
  config.eml_positivity_costheta = 8;
  config.eml_positivity_phi = 12;
  MeasurementResult result = HarmonicMeasurement(grid, config).Fit(response, data);
  REQUIRE(result.response_rows == 8);
  REQUIRE(result.response_columns == 2);
  REQUIRE(result.minimum_detector_intensity > 0.0);
  REQUIRE(result.active_indices == std::vector<std::size_t>{0});
  REQUIRE(result.moment_indices == (std::vector<std::size_t>{0, 1, 2, 3}));
  for (const auto &i : indices(dipole)) {
    const Cell cell = {i, 0, 0};
    REQUIRE(result.flat.at(cell).value[0] == Approx(288.0).epsilon(2e-3));
    REQUIRE(result.fiducial.at(cell).value[0] == Approx(144.0).epsilon(2e-3));
    REQUIRE(result.fiducial.at(cell).value[2] == Approx(dipole[i]).epsilon(2e-3));
    REQUIRE(result.detector.at(cell).value[2] == Approx(dipole[i]).epsilon(2e-3));
  }
  const double cross_covariance = result.fiducial_covariance[2][6];
  const double normalization_covariance =
      0.01 * result.fiducial.at(Cell{0, 0, 0}).value[2] *
      result.fiducial.at(Cell{1, 0, 0}).value[2] / 100.0;
  gra::harmonic::ApplyNormalization(result, 10.0, 0.1);
  REQUIRE(result.fiducial_covariance[2][6] ==
          Approx(cross_covariance / 100.0 + normalization_covariance).epsilon(1e-10));
  for (const auto &i : indices(dipole)) {
    const Cell cell = {i, 0, 0};
    REQUIRE(result.fiducial_covariance[4 * i + 2][4 * i + 2] ==
            Approx(result.fiducial.at(cell).covariance[2][2]).epsilon(1e-10));
    REQUIRE(result.detector_covariance[4 * i + 2][4 * i + 2] ==
            Approx(result.detector.at(cell).covariance[2][2]).epsilon(1e-10));
  }
}

// Check event-by-event response and data weights in the algebraic estimator
TEST_CASE("Algebraic harmonic response preserves weighted sums",
          "[harmonic][weights]") {
  const PhaseSpaceGrid grid(MeasurementMode::Central,
                            {{Coordinate::Mass, 1, 1.0, 3.0},
                             {Coordinate::Momentum, 1, 0.0, 1.0},
                             {Coordinate::Rapidity, 1, -1.0, 1.0}});
  const std::vector<double> z = {2.0, 0.5, 0.0};
  std::vector<ResponseEvent> response;
  for (std::uint64_t index = 0; index < 4; ++index) {
    ResponseEvent event = MakeResponse(
        index, z,
        index < 2 ? std::optional<std::vector<double>>(z) : std::nullopt);
    event.truth.weight = static_cast<double>(index + 1);
    if (event.reco.has_value()) {
      event.reco->weight = event.truth.weight;
    }
    response.push_back(std::move(event));
  }
  const std::vector<DataEvent> data = {{0, "data", MakeObservation(z, 2.0)},
                                       {1, "data", MakeObservation(z, 3.0)}};

  const HarmonicMeasurement measurement(grid, ScalarFitConfig());
  const MeasurementResult result = measurement.Fit(response, data);
  const Cell cell = {0, 0, 0};
  REQUIRE(result.detector.at(cell).value[0] == Approx(5.0));
  REQUIRE(result.detector.at(cell).covariance[0][0] == Approx(13.0));
  REQUIRE(result.flat.at(cell).value[0] == Approx(50.0 / 3.0).epsilon(1e-12));
  REQUIRE(result.flat.at(cell).covariance[0][0] ==
          Approx(1300.0 / 9.0).epsilon(1e-12));
  REQUIRE(result.data_sum_weight == Approx(5.0));
  REQUIRE(result.data_sum_weight2 == Approx(13.0));
}

// Check signed response-MC weights when the net cell measure stays positive
TEST_CASE("Algebraic harmonic response accepts signed MC weights",
          "[harmonic][weights]") {
  const PhaseSpaceGrid grid(MeasurementMode::Central,
                            {{Coordinate::Mass, 1, 1.0, 3.0},
                             {Coordinate::Momentum, 1, 0.0, 1.0},
                             {Coordinate::Rapidity, 1, -1.0, 1.0}});
  const std::vector<double> z = {2.0, 0.5, 0.0};
  std::vector<ResponseEvent> response = {MakeResponse(0, z, z),
                                         MakeResponse(1, z, z),
                                         MakeResponse(2, z, std::nullopt)};
  response[0].truth.weight = 2.0;
  response[1].truth.weight = -1.0;
  response[2].truth.weight = 3.0;
  const std::vector<DataEvent> data = {{0, "data", MakeObservation(z, 4.0)}};

  const HarmonicMeasurement measurement(grid, ScalarFitConfig());
  const MeasurementResult result = measurement.Fit(response, data);
  const Cell cell = {0, 0, 0};
  REQUIRE(result.flat.at(cell).value[0] == Approx(16.0));
  REQUIRE(result.flat.at(cell).covariance[0][0] == Approx(256.0));
}

// Check that weighted response fluctuations enter the unfolded covariance
TEST_CASE("Algebraic harmonic response propagates jackknife uncertainty",
          "[harmonic][weights][jackknife]") {
  const PhaseSpaceGrid grid(MeasurementMode::Central,
                            {{Coordinate::Mass, 1, 1.0, 3.0},
                             {Coordinate::Momentum, 1, 0.0, 1.0},
                             {Coordinate::Rapidity, 1, -1.0, 1.0}});
  const std::vector<double> z = {2.0, 0.5, 0.0};
  const std::vector<ResponseEvent> response = {
      MakeResponse(0, z, z), MakeResponse(1, z, z),
      MakeResponse(2, z, std::nullopt), MakeResponse(3, z, z)};
  const std::vector<DataEvent> data = {{0, "data", MakeObservation(z, 4.0)}};
  FitConfig config = ScalarFitConfig();
  config.response_jackknife_bins = 2;

  const MeasurementResult result =
      HarmonicMeasurement(grid, config).Fit(response, data);
  const Cell cell = {0, 0, 0};
  REQUIRE(result.flat.at(cell).value[0] == Approx(16.0 / 3.0));
  REQUIRE(result.flat.at(cell).covariance[0][0] ==
          Approx(292.0 / 9.0));
  REQUIRE(result.response_jackknife_replicas == 2);
}

// Require estimable truth cells in every response replica while preserving the nominal fit
TEST_CASE("Harmonic response jackknife rejects cells confined to one group",
          "[harmonic][weights][jackknife]") {
  const PhaseSpaceGrid grid(MeasurementMode::Central,
                            {{Coordinate::Mass, 2, 1.0, 3.0},
                             {Coordinate::Momentum, 1, 0.0, 1.0},
                             {Coordinate::Rapidity, 1, -1.0, 1.0}});
  const std::vector<double> a = {1.5, 0.5, 0.0};
  const std::vector<double> b = {2.5, 0.5, 0.0};
  std::vector<ResponseEvent> response = {MakeResponse(0, a, a), MakeResponse(1, b, b)};
  const std::vector<DataEvent> data = {{0, "data", MakeObservation(a, 4.0)},
                                       {1, "data", MakeObservation(b, 6.0)}};
  FitConfig config = ScalarFitConfig();
  const auto nominal = HarmonicMeasurement(grid, config).Fit(response, data);
  REQUIRE(nominal.flat.at(Cell{0, 0, 0}).value[0] == Approx(4.0));
  REQUIRE(nominal.flat.at(Cell{1, 0, 0}).value[0] == Approx(6.0));

  config.response_jackknife_bins = 2;
  for (const auto key : {2U, 4U}) {
    REQUIRE_THROWS_WITH(HarmonicMeasurement(grid, config).Fit(response, data),
                        Catch::Contains("response jackknife group 0 leaves truth cell"));
    response.push_back(MakeResponse(key, a, a));
  }
  response.push_back(MakeResponse(3, a, a));
  response.push_back(MakeResponse(6, b, b));
  const auto result = HarmonicMeasurement(grid, config).Fit(response, data);
  REQUIRE(result.response_jackknife_replicas == 2);
  REQUIRE(result.flat.at(Cell{0, 0, 0}).value[0] == Approx(4.0));
  REQUIRE(result.flat.at(Cell{1, 0, 0}).value[0] == Approx(6.0));
}

// Check the weighted EML yield and sandwich covariance
TEST_CASE("Extended likelihood preserves nonnegative data weights",
          "[harmonic][weights][eml]") {
  const PhaseSpaceGrid grid(MeasurementMode::Central,
                            {{Coordinate::Mass, 1, 1.0, 3.0},
                             {Coordinate::Momentum, 1, 0.0, 1.0},
                             {Coordinate::Rapidity, 1, -1.0, 1.0}});
  const std::vector<double> z = {2.0, 0.5, 0.0};
  std::vector<ResponseEvent> response;
  for (std::uint64_t index = 0; index < 100; ++index) {
    response.push_back(MakeResponse(index, z, z));
  }
  std::vector<DataEvent> data;
  for (std::uint64_t index = 0; index < 100; ++index) {
    data.push_back({index, "data", MakeObservation(z, 2.0)});
  }

  FitConfig config = ScalarFitConfig();
  config.estimator = HarmonicEstimator::EML;
  config.eml_positivity_costheta = 4;
  config.eml_positivity_phi = 8;
  config.response_jackknife_bins = 2;
  const HarmonicMeasurement measurement(grid, config);
  const MeasurementResult result = measurement.Fit(response, data);
  const Cell cell = {0, 0, 0};
  REQUIRE(result.flat.at(cell).value[0] == Approx(200.0).epsilon(1e-6));
  REQUIRE(result.flat.at(cell).covariance[0][0] == Approx(400.0).epsilon(1e-5));
  REQUIRE(result.detector.at(cell).value[0] == Approx(200.0).epsilon(1e-6));
  REQUIRE(result.estimator == HarmonicEstimator::EML);
  REQUIRE(result.response_jackknife_replicas == 2);
}

// Check that signed data weights are not assigned an invalid likelihood
TEST_CASE("Extended likelihood rejects signed data weights",
          "[harmonic][weights][eml]") {
  const PhaseSpaceGrid grid(MeasurementMode::Central,
                            {{Coordinate::Mass, 1, 1.0, 3.0},
                             {Coordinate::Momentum, 1, 0.0, 1.0},
                             {Coordinate::Rapidity, 1, -1.0, 1.0}});
  const std::vector<double> z = {2.0, 0.5, 0.0};
  const std::vector<ResponseEvent> response = {MakeResponse(0, z, z),
                                               MakeResponse(1, z, z)};
  const std::vector<DataEvent> data = {{0, "data", MakeObservation(z, -1.0)}};
  FitConfig config = ScalarFitConfig();
  config.estimator = HarmonicEstimator::EML;
  const HarmonicMeasurement measurement(grid, config);
  REQUIRE_THROWS_AS(measurement.Fit(response, data), std::invalid_argument);
}

// Check that unresolved feed-in across the outer truth boundary cannot be
// hidden
TEST_CASE("Harmonic response rejects feed-in from outside the truth grid",
          "[harmonic][migration]") {
  const PhaseSpaceGrid grid(MeasurementMode::Central,
                            {{Coordinate::Mass, 1, 1.0, 3.0},
                             {Coordinate::Momentum, 1, 0.0, 1.0},
                             {Coordinate::Rapidity, 1, -1.0, 1.0}});
  const std::vector<double> inside = {2.0, 0.5, 0.0};
  const std::vector<double> outside = {3.5, 0.5, 0.0};
  const std::vector<ResponseEvent> response = {
      MakeResponse(0, inside, inside), MakeResponse(1, outside, inside)};
  const std::vector<DataEvent> data = {{0, "data", MakeObservation(inside)}};

  const HarmonicMeasurement measurement(grid, ScalarFitConfig());
  REQUIRE_THROWS_AS(measurement.Fit(response, data), std::invalid_argument);
}

// Check that tracking efficiency does not bias the independent momentum resolution
TEST_CASE("Toy detector efficiency is independent of momentum smearing", "[harmonic][response]") {
  ToyResponseConfig config;
  config.pion_efficiency = std::sqrt(0.5);
  config.pion_logpt_sigma = 0.1;
  config.pion_eta_sigma = 0.1;
  config.pion_phi_sigma = 0.1;
  const ToyResponse model(config);
  const EventKinematics truth = MakeKinematics();
  double sum_logpt = 0.0, sum_eta = 0.0, square_logpt = 0.0, square_eta = 0.0;
  std::size_t accepted = 0;
  for (std::uint64_t key = 0; key < 50000; ++key) {
    const auto reco = gra::harmonic::SimulateDetector(model, truth, key, config.seed);
    if (!reco) { continue; }
    ++accepted;
    const double logpt = std::log(reco->pip.Pt() / truth.pip.Pt()) + 0.005;
    const double eta = reco->pip.Eta() - truth.pip.Eta();
    sum_logpt += logpt;
    sum_eta += eta;
    square_logpt += logpt * logpt;
    square_eta += eta * eta;
  }
  REQUIRE(accepted > 24000);
  REQUIRE(accepted < 26000);
  CHECK(std::abs(sum_logpt / accepted) < 0.004);
  CHECK(std::abs(sum_eta / accepted) < 0.004);
  CHECK(square_logpt / accepted == Approx(0.01).margin(0.0005));
  CHECK(square_eta / accepted == Approx(0.01).margin(0.0005));
}

// Check physical migration into and out of pion cuts using reconstructed tracks
TEST_CASE("Toy response applies fiducial cuts after particle reconstruction",
          "[harmonic][hepmc][response][fiducial]") {
  const EventKinematics truth = MakeKinematics();
  const auto event = KinematicsEvent(truth);
  const auto path = analysis_test::Write("harmonic_cut_migrations.hepmc3",
                                        std::vector<HepMC3::GenEvent>(128, event));
  gra::harmonic::HepMCReadConfig config;
  config.frame = AngularFrame::CM;
  config.sqrt_s = 13000.0;
  config.max_records = 128;
  config.cuts.pion_eta = {-5.0, 5.0};
  bool fiducial = false;
  SECTION("feed in from below the pion threshold") {
    config.cuts.pion_pt = {1.01 * truth.pip.Pt(), 10.0};
  }
  SECTION("loss above the pion threshold") {
    fiducial = true;
    config.cuts.pion_pt = {0.99 * truth.pip.Pt(), 10.0};
  }
  ToyResponseConfig resolution;
  resolution.pion_logpt_sigma = 0.2;
  const auto events = gra::harmonic::HepMCAnalysisReader(config).ReadModelResponse(path, ToyResponse(resolution));
  REQUIRE(events.size() == 128);
  std::size_t accepted = 0;
  for (const ResponseEvent &response : events) {
    REQUIRE(response.fiducial.central == fiducial);
    REQUIRE(response.reco.has_value());
    REQUIRE(response.reco->z[0] >= 2.0 * gra::PDG::mpi);
    REQUIRE(response.reco->z[1] >= 0.0);
    if (response.DetectorAccepted()) { ++accepted; }
  }
  REQUIRE(accepted > 0);
  REQUIRE(accepted < events.size());
}

// Check tagged exclusivity on the measured proton and pion four-momenta
TEST_CASE("Toy proton resolution changes reconstructed exclusivity selection",
          "[harmonic][hepmc][response][tagged]") {
  const auto event = KinematicsEvent(MakeKinematics());
  const auto path = analysis_test::Write("harmonic_tagged_resolution.hepmc3",
                                        std::vector<HepMC3::GenEvent>(128, event));
  gra::harmonic::HepMCReadConfig config;
  config.mode = MeasurementMode::Tagged;
  config.frame = AngularFrame::CM;
  config.sqrt_s = 13000.0;
  config.max_records = 128;
  config.cuts.pion_eta = {-5.0, 5.0};
  config.cuts.pion_pt = {0.1, 10.0};
  config.cuts.proton_xi = {0.0, 1.0};
  config.cuts.proton_abs_t = {0.0, 5.0};
  config.selection = {0.1, 1.0, 1.0};
  const gra::harmonic::HepMCAnalysisReader reader(config);
  const auto identity = reader.ReadModelResponse(path, gra::harmonic::IdentityResponse());
  for (const ResponseEvent &response : identity) { REQUIRE(response.DetectorAccepted()); }
  ToyResponseConfig resolution;
  resolution.proton_pt_sigma = 0.2;
  const auto events = reader.ReadModelResponse(path, ToyResponse(resolution));
  std::size_t accepted = 0;
  for (const ResponseEvent &response : events) {
    REQUIRE(response.fiducial.Pass(config.mode));
    REQUIRE(response.reco.has_value());
    REQUIRE(response.reco->z[4] >= -gra::math::PI);
    REQUIRE(response.reco->z[4] <= gra::math::PI);
    if (response.DetectorAccepted()) { ++accepted; }
  }
  REQUIRE(accepted > 0);
  REQUIRE(accepted < events.size());
}

// Check unit conversion and malformed input through all three harmonic readers
TEST_CASE("Harmonic HepMC readers normalize units and reject broken records", "[harmonic][hepmc]") {
  gra::harmonic::HepMCReadConfig config;
  config.frame = AngularFrame::CM;
  config.sqrt_s = 13000.0;
  config.max_records = 10;
  config.cuts.pion_eta = {-5.0, 5.0};
  config.cuts.pion_pt = {0.1, 2.0};
  const gra::harmonic::HepMCAnalysisReader reader(config);
  const gra::harmonic::IdentityResponse identity;
  const auto event = analysis_test::Event();
  const auto gev = analysis_test::Write("harmonic_gev.hepmc3", {event, event});
  const auto mev = analysis_test::Write("harmonic_mev.hepmc3", {event, event}, HepMC3::Units::MEV);
  const auto data = reader.ReadData(gev);
  const auto data_mev = reader.ReadData(mev);
  const auto response = reader.ReadModelResponse(mev, identity);
  const auto paired = reader.ReadPairedResponse(gev, mev);
  REQUIRE(data.size() == 2);
  REQUIRE(data_mev.size() == 2);
  REQUIRE(response.size() == 2);
  REQUIRE(paired.size() == 2);
  for (const auto &i : indices(data)) {
    REQUIRE(paired[i].reco.has_value());
    CHECK(response[i].DetectorAccepted());
    CHECK(paired[i].DetectorAccepted());
    for (const auto &j : indices(data[i].reco.z)) {
      CHECK(data_mev[i].reco.z[j] == Approx(data[i].reco.z[j]));
      CHECK(response[i].truth.z[j] == Approx(data[i].reco.z[j]));
      CHECK(paired[i].reco->z[j] == Approx(data[i].reco.z[j]));
    }
    CHECK(data_mev[i].reco.costheta == Approx(data[i].reco.costheta));
    CHECK(data_mev[i].reco.phi == Approx(data[i].reco.phi));
  }
  for (const auto &suffix : {std::string("E 1 1 2\nP broken\n"), std::string("E 1 1 2\nU GEV MM\n")}) {
    const auto broken = analysis_test::Broken("harmonic_broken.hepmc3", suffix);
    CHECK_THROWS_AS(reader.ReadData(broken), std::invalid_argument);
    CHECK_THROWS_AS(reader.ReadModelResponse(broken, identity), std::invalid_argument);
    CHECK_THROWS_AS(reader.ReadPairedResponse(broken, gev), std::invalid_argument);
    CHECK_THROWS_AS(reader.ReadPairedResponse(gev, broken), std::invalid_argument);
  }
}

// Keep the final complete event when EOF follows its last particle without a footer
TEST_CASE("Analysis HepMC reader accepts a complete final event at EOF", "[harmonic][hepmc]") {
  const auto path = analysis_test::Broken("harmonic_no_footer.hepmc3", "");
  gra::MHepMCReader reader(path);
  HepMC3::GenEvent event;
  REQUIRE(reader.Read(event));
  CHECK(event.particles().size() == 3);
  CHECK_FALSE(reader.Read(event));
  CHECK(event.particles().empty());
  CHECK_FALSE(reader.Read(event));
  CHECK(event.particles().empty());
}

// Use the generated systematic weight when reconstructed files store only nominal weights
TEST_CASE("Paired harmonic response takes its weight from truth", "[harmonic][hepmc][weights]") {
  gra::harmonic::HepMCReadConfig config;
  config.frame = AngularFrame::CM;
  config.sqrt_s = 13000.0;
  config.max_records = 10;
  config.weight_index = 1;
  config.cuts.pion_eta = {-5.0, 5.0};
  config.cuts.pion_pt = {0.1, 2.0};
  auto truth = analysis_test::Event();
  truth.weights() = {1.0, 2.0};
  const auto truth_path = analysis_test::Write("harmonic_truth_weights.hepmc3", {truth});
  const auto reco_path = analysis_test::Write("harmonic_reco_nominal.hepmc3", {analysis_test::Event()});
  const auto paired = gra::harmonic::HepMCAnalysisReader(config).ReadPairedResponse(truth_path, reco_path);
  REQUIRE(paired.size() == 1);
  REQUIRE(paired.front().DetectorAccepted());
  CHECK(paired.front().truth.weight == Approx(2.0));
  CHECK(paired.front().reco->weight == Approx(2.0));
}

// Reject incomplete event and vertex records even when EOF is reached during parsing
TEST_CASE("Analysis HepMC reader rejects malformed final records", "[harmonic][hepmc]") {
  for (const std::string suffix : {"E broken\n", "E 1 1 0\nU GEV MM\n",
                                   "E 1 0 1\nU GEV MM\nP 1 0 211 0 0 1 1.1 0.1\n"}) {
    CAPTURE(suffix);
    const auto path = analysis_test::Broken("harmonic_final_broken.hepmc3", suffix);
    gra::MHepMCReader reader(path);
    HepMC3::GenEvent event;
    REQUIRE(reader.Read(event));
    CHECK_THROWS_AS(reader.Read(event), std::invalid_argument);
  }
}

// Check parity and azimuthal covariance through the analysis spherical basis
TEST_CASE("Analysis spherical harmonics preserve parity and rotations", "[harmonic][spherical][physics]") {
  gra::spherical::Omega point;
  point.costheta = 0.37;
  point.phi = -0.63;
  auto parity = point;
  parity.costheta = -point.costheta;
  parity.phi += gra::math::PI;
  const double angle = 0.41;
  auto rotated = point;
  rotated.phi += angle;
  const auto basis = gra::spherical::YLM({point, parity, rotated}, 4);
  for (int l = 0; l <= 4; ++l) {
    for (int m = -l; m <= l; ++m) {
      const auto index = gra::spherical::LinearInd(l, m);
      CHECK(basis[1][index] == Approx((l % 2 == 0 ? 1.0 : -1.0) * basis[0][index]).margin(1e-14));
      const double c = std::cos(std::abs(m) * angle);
      const double s = std::sin(std::abs(m) * angle);
      const double other = basis[0][gra::spherical::LinearInd(l, -m)];
      CHECK(basis[2][index] == Approx(c * basis[0][index] + (m < 0 ? s : -s) * other).margin(1e-14));
    }
  }
}

// Reject impossible angular truncations before spherical response allocation
TEST_CASE("Analysis spherical functions reject overflowing basis sizes", "[harmonic][spherical][config]") {
  for (const int lmax : {-2, std::numeric_limits<int>::max()}) {
    CAPTURE(lmax);
    const std::vector<gra::spherical::Omega> events(2);
    CHECK_THROWS_AS(gra::spherical::YLM(events, lmax), std::invalid_argument);
    CHECK_THROWS_AS(gra::spherical::SphericalMoments(events, {0, 1}, lmax, "fla"), std::invalid_argument);
    CHECK_THROWS_AS(gra::spherical::GetELM(events, {0, 1}, lmax, "fla"), std::invalid_argument);
    CHECK_THROWS_AS(gra::spherical::GetGMixing(events, {0, 1}, lmax, "fla"), std::invalid_argument);
    CHECK_THROWS_AS(gra::spherical::HarmDotProd({1.0}, {1.0}, {true}, lmax), std::invalid_argument);
    CHECK_THROWS_AS(gra::spherical::HarmDotProdError({1.0}, gra::MMatrix<double>{{1.0}}, {true}, lmax),
                    std::invalid_argument);
  }
}

// Reject impossible harmonic dimensions before allocating event response matrices
TEST_CASE("Harmonic fit validates basis size before allocating", "[harmonic][config]") {
  FitConfig config;
  SECTION("integer overflow in full basis") {
    config.lmax = std::numeric_limits<int>::max();
    CHECK_THROWS_AS(config.Validate(), std::invalid_argument);
  }
  SECTION("one cell exceeds fit size") {
    config.lmax = 2;
    config.max_parameters = 8;
    CHECK_THROWS_AS(config.Validate(), std::invalid_argument);
    config.remove_negative_m = true;
    CHECK_THROWS_AS(config.Validate(), std::invalid_argument);
    config.remove_odd = true;
    CHECK_THROWS_AS(config.Validate(), std::invalid_argument);
    config.max_parameters = 9;
    CHECK_NOTHROW(config.Validate());
  }
  SECTION("EML grid size overflow") {
    config.eml_positivity_costheta = std::numeric_limits<std::size_t>::max();
    CHECK_THROWS_AS(config.Validate(), std::invalid_argument);
  }
  SECTION("EML call count truncation") {
    config.eml_max_calls = static_cast<std::size_t>(std::numeric_limits<unsigned int>::max()) + 1;
    CHECK_THROWS_AS(config.Validate(), std::invalid_argument);
  }
}

// Check that conditional response normalization is independent of MC weight units
TEST_CASE("Harmonic response is invariant under finite weight rescaling", "[harmonic][weights]") {
  const PhaseSpaceGrid grid(MeasurementMode::Central,
                           {{Coordinate::Mass, 2, 1.0, 3.0},
                            {Coordinate::Momentum, 1, 0.0, 1.0},
                            {Coordinate::Rapidity, 1, -1.0, 1.0}});
  FitConfig config = ScalarFitConfig();
  config.response_jackknife_bins = 2;
  std::vector<ResponseEvent> response;
  std::vector<DataEvent> data;
  for (std::size_t cell = 0; cell < 2; ++cell) {
    const std::vector<double> z = {1.5 + cell, 0.5, 0.0};
    for (std::uint64_t i = 0; i < 8; ++i) {
      response.push_back(MakeResponse(i, z, i < 4 ? std::optional(z) : std::nullopt));
    }
    data.push_back({cell, "data", MakeObservation(z, 3.0)});
  }
  const HarmonicMeasurement measurement(grid, config);
  const auto expected = measurement.Fit(response, data);
  for (const double scale : {1e-310, 1e308}) {
    CAPTURE(scale);
    auto scaled = response;
    for (std::size_t i = 0; i < 8; ++i) { scaled[i].truth.weight = scale; }
    const auto result = measurement.Fit(scaled, data);
    for (const auto &cell : result.truth_cells) {
      CHECK(result.flat.at(cell).value[0] == Approx(expected.flat.at(cell).value[0]).epsilon(1e-12));
      CHECK(result.flat.at(cell).covariance[0][0] == Approx(expected.flat.at(cell).covariance[0][0]).epsilon(1e-12));
    }
  }
}

// Check zero-weight observations do not require response in otherwise unoccupied cells
TEST_CASE("Harmonic data ignore zero-weight cells", "[harmonic][weights]") {
  const PhaseSpaceGrid grid(MeasurementMode::Central,
                           {{Coordinate::Mass, 2, 1.0, 3.0},
                            {Coordinate::Momentum, 1, 0.0, 1.0},
                            {Coordinate::Rapidity, 1, -1.0, 1.0}});
  const std::vector<double> z = {1.5, 0.5, 0.0};
  const std::vector<ResponseEvent> response = {MakeResponse(0, z, z)};
  const std::vector<DataEvent> data = {{0, "data", MakeObservation(z, 2.0)},
                                      {1, "data", MakeObservation({2.5, 0.5, 0.0}, 0.0)}};
  const auto result = HarmonicMeasurement(grid, ScalarFitConfig()).Fit(response, data);
  CHECK(result.data_events == 1);
  CHECK(result.flat.at(Cell{0, 0, 0}).value[0] == Approx(2.0));
}

// Check compensated angular moments after subtraction of signed event weights
TEST_CASE("Harmonic moments preserve signed weight cancellation", "[harmonic][spherical][weights]") {
  const PhaseSpaceGrid grid(MeasurementMode::Central,
                           {{Coordinate::Mass, 1, 1.0, 3.0},
                            {Coordinate::Momentum, 1, 0.0, 1.0},
                            {Coordinate::Rapidity, 1, -1.0, 1.0}});
  const std::vector<double> z = {2.0, 0.5, 0.0};
  std::vector<DataEvent> data;
  std::vector<gra::spherical::Omega> events;
  for (const double weight : {1e20, 1.0, -1e20}) {
    data.push_back({data.size(), "data", MakeObservation(z, weight)});
    gra::spherical::Omega event;
    event.weight = weight;
    events.push_back(event);
  }
  const auto moments = gra::spherical::SphericalMoments(events, {0, 1, 2}, 0, "fla");
  CHECK(moments.sum_weight == Approx(1.0));
  CHECK(moments.value[0] == Approx(1.0));
  CHECK(moments.covariance[0][0] == Approx(2e40));
  const auto result = HarmonicMeasurement(grid, ScalarFitConfig()).Fit({MakeResponse(0, z, z)}, data);
  CHECK(result.data_sum_weight == Approx(1.0));
  CHECK(result.flat.at(Cell{0, 0, 0}).value[0] == Approx(1.0));
  CHECK(result.flat.at(Cell{0, 0, 0}).covariance[0][0] == Approx(2e40));
}

// Check synthesis normalization and independent errors without squaring their scale
TEST_CASE("Spherical synthesis and errors preserve finite scales", "[harmonic][spherical][numerics]") {
  for (const double scale : {1e-310, 1e-200, 1e200}) {
    CAPTURE(scale);
    const auto error = gra::spherical::ErrorProp(gra::MMatrix<double>{{3.0, 4.0}}, {scale, scale});
    CHECK(error[0] / scale == Approx(5.0).epsilon(1e-12));
    std::vector<double> costheta, phi;
    const auto values = gra::spherical::Y_real_synthesize({scale}, {true}, 3, costheta, phi, true);
    CHECK(values.IsFinite());
    CHECK(values[1][1] == Approx(1.0).epsilon(1e-12));
  }
  CHECK_THROWS_AS(gra::spherical::CalcError(1.0, 1e200, 2.0), std::domain_error);
  CHECK_THROWS_AS(gra::spherical::CalcError(0.0, 1e-100, 2.0), std::domain_error);
}

// Check exact internal boundaries and their neighboring representable coordinates
TEST_CASE("Harmonic grid assigns decimal boundaries to the upper cell", "[harmonic][binning]") {
  const PhaseSpaceGrid grid(MeasurementMode::Central,
                           {{Coordinate::Mass, 10, 1.0, 3.0},
                            {Coordinate::Momentum, 1, 0.0, 1.0},
                            {Coordinate::Rapidity, 1, -1.0, 1.0}});
  for (std::size_t i = 1; i < 10; ++i) {
    const double edge = 1.0 + 0.2 * i;
    CHECK(grid.Locate({edge, 0.5, 0.0}) == std::optional(Cell{i, 0, 0}));
    CHECK(grid.Locate({std::nextafter(edge, 0.0), 0.5, 0.0}) == std::optional(Cell{i - 1, 0, 0}));
    CHECK(grid.Locate({std::nextafter(edge, 3.0), 0.5, 0.0}) == std::optional(Cell{i, 0, 0}));
  }
  CHECK_THROWS_AS((Axis{Coordinate::Mass, 2, -std::numeric_limits<double>::max(),
                        std::numeric_limits<double>::max()}.Validate()), std::invalid_argument);
}

// Check covariance units when the square of the normalization is not representable
TEST_CASE("Harmonic normalization retains representable covariances", "[harmonic][normalization]") {
  const PhaseSpaceGrid grid(MeasurementMode::Central,
                           {{Coordinate::Mass, 1, 1.0, 3.0},
                            {Coordinate::Momentum, 1, 0.0, 1.0},
                            {Coordinate::Rapidity, 1, -1.0, 1.0}});
  const std::vector<double> z = {2.0, 0.5, 0.0};
  auto result = HarmonicMeasurement(grid, ScalarFitConfig()).Fit(
      {MakeResponse(0, z, z)}, {{0, "data", MakeObservation(z, 1e150)}});
  gra::harmonic::ApplyNormalization(result, 1e200, 0.0);
  const Cell cell{0, 0, 0};
  CHECK(result.flat.at(cell).value[0] / 1e-50 == Approx(1.0));
  CHECK(result.flat.at(cell).covariance[0][0] / 1e-100 == Approx(1.0));
  CHECK(result.flat_covariance[0][0] / 1e-100 == Approx(1.0));
}

// Check all real basis overlaps and the efficiency normalization by exact angular quadrature
TEST_CASE("Spherical response and synthesis share an orthonormal basis", "[harmonic][spherical][physics]") {
  const auto [node, weight] = gra::math::GaussLegendreRule(6, -1.0, 1.0);
  std::vector<gra::spherical::Omega> events;
  std::vector<std::size_t> ind;
  for (const auto &i : indices(node)) {
    for (int j = 0; j < 12; ++j) {
      gra::spherical::Omega event;
      event.costheta = node[i];
      event.phi = 2.0 * gra::math::PI * (j + 0.5) / 12.0;
      event.weight = weight[i] / 24.0;
      event.fiducial = event.selected = true;
      ind.push_back(events.size());
      events.push_back(event);
    }
  }
  const auto response = gra::spherical::GetGMixing(events, ind, 4, "det");
  const auto efficiency = gra::spherical::GetELM(events, ind, 4, "det");
  for (const auto &i : indices(efficiency.first)) {
    CHECK(efficiency.first[i] == Approx(i == 0 ? 1.0 : 0.0).margin(1e-13));
    for (const auto &j : indices(efficiency.first)) {
      CHECK(response.value[i][j] == Approx(i == j ? 1.0 : 0.0).margin(1e-13));
    }
  }
  std::vector<double> coefficients(25);
  for (const auto &i : indices(coefficients)) { coefficients[i] = std::sin(i + 0.3); }
  std::vector<double> costheta, phi;
  const auto surface = gra::spherical::Y_real_synthesize(coefficients, std::vector<bool>(25, true),
                                                        7, costheta, phi);
  for (const auto &i : indices(costheta)) {
    for (const auto &j : indices(phi)) {
      gra::spherical::Omega event;
      event.costheta = costheta[i];
      event.phi = phi[j];
      CHECK(surface[i][j] == Approx((gra::spherical::YLM({event}, 4) * coefficients)[0]).margin(1e-13));
    }
  }
  events.assign(2, gra::spherical::Omega{});
  events[1].weight = 0.0;
  CHECK_THROWS_AS(gra::spherical::GetGMixing(events, {0, 1}, 0, "fla"), std::invalid_argument);
}

// Check weighted vector covariance independently of the angular basis
TEST_CASE("Weighted vector sums retain correlations and reject overflow", "[harmonic][statistics]") {
  gra::statistics::WeightedVectorSums sums(2);
  sums.Add(std::vector<double>{1.0, 2.0}, 3.0);
  sums.Add(std::vector<double>{2.0, -1.0}, -2.0);
  const auto mean = sums.Sum();
  const auto covariance = sums.Covariance();
  CHECK(mean[0] == Approx(-1.0));
  CHECK(mean[1] == Approx(8.0));
  CHECK(covariance[0][0] == Approx(25.0));
  CHECK(covariance[0][1] == Approx(10.0));
  CHECK(covariance[1][0] == Approx(10.0));
  CHECK(covariance[1][1] == Approx(40.0));
  CHECK_THROWS_AS(sums.Add(std::vector<double>{1.0}, 1.0), std::invalid_argument);
  CHECK_THROWS_AS(sums.Add(std::vector<double>{1.0, 2.0}, std::numeric_limits<double>::infinity()),
                  std::invalid_argument);
  gra::statistics::WeightedVectorSums large(1);
  large.Add(std::vector<double>{1.0}, 1e200);
  CHECK(large.Sum()[0] / 1e200 == Approx(1.0));
  CHECK_THROWS_AS(large.Covariance(), std::overflow_error);
}

// Check a general SO(3) rotation of an anisotropic distribution through the fitted response
TEST_CASE("Harmonic fit unfolds a rotated dipole", "[harmonic][spherical][physics]") {
  const PhaseSpaceGrid grid(MeasurementMode::Central,
                           {{Coordinate::Mass, 1, 1.0, 3.0},
                            {Coordinate::Momentum, 1, 0.0, 1.0},
                            {Coordinate::Rapidity, 1, -1.0, 1.0}});
  const std::vector<double> z = {2.0, 0.5, 0.0};
  const auto [node, weight] = gra::math::GaussLegendreRule(6, -1.0, 1.0);
  std::vector<ResponseEvent> response;
  std::vector<DataEvent> data;
  const gra::M4Vec dipole(0.12, -0.08, 0.16, 0.0);
  for (const auto &i : indices(node)) {
    for (int j = 0; j < 12; ++j) {
      const double phi = -gra::math::PI + 2.0 * gra::math::PI * (j + 0.5) / 12.0;
      const double sine = std::sqrt(1.0 - node[i] * node[i]);
      gra::M4Vec direction(sine * std::cos(phi), sine * std::sin(phi), node[i], 1.0);
      const double density = 1.0 + std::sqrt(3.0) * gra::BilinearProduct(dipole.P3(), direction.P3());
      ResponseEvent event = MakeResponse(response.size(), z, z);
      event.truth = MakeAngularObservation(z, node[i], phi, weight[i] / 24.0);
      direction.RotateY(0.71);
      direction.RotateZ(-0.39);
      event.reco = MakeAngularObservation(z, direction.CosTheta(), direction.Phi(), event.truth.weight);
      response.push_back(event);
      auto reco = *event.reco;
      reco.weight *= 100.0 * density;
      data.push_back({data.size(), "data", reco});
    }
  }
  auto config = ScalarFitConfig();
  config.lmax = 1;
  const auto result = HarmonicMeasurement(grid, config).Fit(response, data);
  const auto &flat = result.flat.at(Cell{0, 0, 0}).value;
  CHECK(flat[0] == Approx(100.0).epsilon(1e-12));
  CHECK(flat[1] == Approx(100.0 * dipole.Py()).epsilon(1e-12));
  CHECK(flat[2] == Approx(100.0 * dipole.Pz()).epsilon(1e-12));
  CHECK(flat[3] == Approx(100.0 * dipole.Px()).epsilon(1e-12));
  auto rotated = dipole;
  rotated.RotateY(0.71);
  rotated.RotateZ(-0.39);
  const auto &det = result.detector.at(Cell{0, 0, 0}).value;
  CHECK(det[1] == Approx(100.0 * rotated.Py()).epsilon(1e-12));
  CHECK(det[2] == Approx(100.0 * rotated.Pz()).epsilon(1e-12));
  CHECK(det[3] == Approx(100.0 * rotated.Px()).epsilon(1e-12));
}

// Check propagated errors in covariance units and before intermediate squares overflow
TEST_CASE("Spherical covariance errors use the physical variance scale", "[harmonic][spherical][numerics]") {
  const gra::MMatrix<double> negative{{-1e-100}};
  CHECK_THROWS_AS(gra::spherical::CovarianceErrors(negative), std::domain_error);
  CHECK_THROWS_AS(gra::spherical::HarmDotProdError({1.0}, negative, {true}, 0), std::domain_error);
  const auto error = gra::spherical::HarmDotProdError({1e100}, gra::MMatrix<double>{{1e200}}, {true}, 0);
  CHECK(error / 1e200 == Approx(1.0));
  CHECK_THROWS_AS(gra::spherical::SummarizeWeights({}, {}, "unknown"), std::invalid_argument);
}
