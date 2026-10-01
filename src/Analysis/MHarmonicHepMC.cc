// HepMC input adapter for detector-aware angular measurements
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// HepMC3
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"

// Own
#include "Graniitti/Analysis/MHarmonicHepMC.h"
#include "Graniitti/Analysis/MHepMCReader.h"
#include "Graniitti/Kinematics/M4Vec.h"
#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Particle/MPDG.h"
#include "Graniitti/Tech/MAux.h"

namespace gra {
namespace harmonic {

namespace {

struct ProjectedRecord {
  Observation observation;
  FiducialDecision fiducial;
  bool detector_selected = false;
};

// Build a stable source-dependent event key
std::uint64_t EventKey(const std::string &source, std::uint64_t event_number) {
  std::uint64_t hash = 1469598103934665603ULL;
  for (const unsigned char character : source) {
    hash ^= character;
    hash *= 1099511628211ULL;
  }
  hash ^= event_number + 0x9e3779b97f4a7c15ULL + (hash << 6U) + (hash >> 2U);
  return hash;
}

// Compute the configured finite HepMC event weight
double NominalWeight(const HepMC3::GenEvent &event, std::size_t weight_index) {
  if (event.weights().empty()) {
    if (weight_index != 0) {
      throw std::invalid_argument(
          "harmonic::HepMCAnalysisReader weight index is unavailable");
    }
    return 1.0;
  }
  if (weight_index >= event.weights().size()) {
    throw std::invalid_argument(
        "harmonic::HepMCAnalysisReader weight index is unavailable");
  }
  const double weight = event.weights()[weight_index];
  if (!std::isfinite(weight)) {
    throw std::invalid_argument(
        "harmonic::HepMCAnalysisReader found a non-finite event weight");
  }
  return weight;
}

// Select exactly one stable pi+ pi- two-body central state
bool SelectCentralPair(const HepMC3::GenEvent &event, M4Vec &pip, M4Vec &pim) {
  std::size_t pip_count = 0;
  std::size_t pim_count = 0;
  std::size_t nonproton_count = 0;
  for (const auto &particle : event.particles()) {
    if (particle->status() != PDG::PDG_STABLE ||
        particle->pid() == PDG::PDG_p) {
      continue;
    }
    ++nonproton_count;
    if (particle->pid() == PDG::PDG_pip) {
      pip = gra::aux::HepMC2M4Vec(particle->momentum());
      ++pip_count;
    } else if (particle->pid() == PDG::PDG_pim) {
      pim = gra::aux::HepMC2M4Vec(particle->momentum());
      ++pim_count;
    }
  }
  return nonproton_count == 2 && pip_count == 1 && pim_count == 1;
}

// Select one final proton in each beam hemisphere
bool SelectForwardProtons(const HepMC3::GenEvent &event, M4Vec &plus,
                          M4Vec &minus) {
  std::size_t plus_count = 0;
  std::size_t minus_count = 0;
  for (const auto &particle : event.particles()) {
    if (particle->status() != PDG::PDG_STABLE ||
        particle->pid() != PDG::PDG_p) {
      continue;
    }
    const M4Vec momentum = gra::aux::HepMC2M4Vec(particle->momentum());
    if (momentum.Pz() > 0.0) {
      plus = momentum;
      ++plus_count;
    } else if (momentum.Pz() < 0.0) {
      minus = momentum;
      ++minus_count;
    }
  }
  return plus_count == 1 && minus_count == 1;
}

// Select beam records or construct nominal beams for detector-level data
void SelectBeams(const HepMC3::GenEvent &event, double sqrt_s, M4Vec &plus,
                 M4Vec &minus) {
  std::size_t plus_count = 0;
  std::size_t minus_count = 0;
  for (const auto &particle : event.particles()) {
    if (particle->status() != PDG::PDG_BEAM || particle->pid() != PDG::PDG_p) {
      continue;
    }
    const M4Vec momentum = gra::aux::HepMC2M4Vec(particle->momentum());
    if (momentum.Pz() > 0.0) {
      plus = momentum;
      ++plus_count;
    } else if (momentum.Pz() < 0.0) {
      minus = momentum;
      ++minus_count;
    }
  }
  if (plus_count == 1 && minus_count == 1) {
    return;
  }
  const double energy = sqrt_s / 2.0;
  const double momentum =
      std::sqrt(std::max(0.0, energy * energy - PDG::mp * PDG::mp));
  plus = M4Vec(0.0, 0.0, momentum, energy);
  minus = M4Vec(0.0, 0.0, -momentum, energy);
}

// Extract the topology and beam convention needed by one measurement mode
std::optional<EventKinematics> SelectEvent(const HepMC3::GenEvent &event,
                                           const HepMCReadConfig &config,
                                           std::optional<double> weight = std::nullopt) {
  EventKinematics selected;
  selected.mode = config.mode;
  if (!SelectCentralPair(event, selected.pip, selected.pim)) {
    return std::nullopt;
  }
  SelectBeams(event, config.sqrt_s, selected.beam_plus, selected.beam_minus);
  selected.has_forward_protons =
      SelectForwardProtons(event, selected.proton_plus, selected.proton_minus);
  if (config.mode == MeasurementMode::Tagged && !selected.has_forward_protons) {
    return std::nullopt;
  }
  selected.weight = weight.has_value() ? *weight : NominalWeight(event, config.weight_index);
  return selected;
}

// Test one pion against the central particle-level fiducial definition
bool AcceptPion(const M4Vec &pion, const FiducialCuts &cuts) {
  return cuts.pion_eta.Contains(pion.Eta()) && cuts.pion_pt.Contains(pion.Pt());
}

// Calculate light-cone xi and absolute momentum transfer for one beam arm
std::pair<double, double> ForwardCoordinates(const M4Vec &beam,
                                             const M4Vec &proton) {
  if (std::fpclassify(beam.Pz()) == FP_ZERO ||
      beam.Pz() * proton.Pz() <= 0.0) {
    throw std::invalid_argument(
        "harmonic::HepMCAnalysisReader found an invalid proton arm");
  }
  const double direction = beam.Pz() > 0.0 ? 1.0 : -1.0;
  const double beam_lightcone = beam.E() + direction * beam.Pz();
  const double proton_lightcone = proton.E() + direction * proton.Pz();
  if (!(beam_lightcone > 0.0)) {
    throw std::invalid_argument(
        "harmonic::HepMCAnalysisReader found an invalid beam momentum");
  }
  const double xi = 1.0 - proton_lightcone / beam_lightcone;
  const double abs_t = std::abs((beam - proton).M2());
  return {xi, abs_t};
}

// Test both forward protons against the tagged fiducial definition
bool AcceptForward(const EventKinematics &event, const FiducialCuts &cuts) {
  if (!event.has_forward_protons) {
    return false;
  }
  const auto plus = ForwardCoordinates(event.beam_plus, event.proton_plus);
  const auto minus = ForwardCoordinates(event.beam_minus, event.proton_minus);
  return cuts.proton_xi.Contains(plus.first) &&
         cuts.proton_abs_t.Contains(plus.second) &&
         cuts.proton_xi.Contains(minus.first) &&
         cuts.proton_abs_t.Contains(minus.second);
}

// Test tagged central-forward closure using reconstructed four-vectors
bool AcceptExclusiveClosure(const EventKinematics &event,
                            const ExclusiveSelection &selection) {
  if (!event.has_forward_protons) {
    return false;
  }
  const M4Vec central = event.pip + event.pim;
  const M4Vec missing = event.beam_plus + event.beam_minus - event.proton_plus -
                        event.proton_minus;
  const M4Vec residual = missing - central;
  if (!(central.M() > 0.0) || !(missing.M2() > 0.0)) {
    return false;
  }
  const double mass_match = std::abs(central.M() - missing.M()) / central.M();
  const double rapidity_match = std::abs(central.Rap() - missing.Rap());
  return residual.Pt() <= selection.pt_balance_max &&
         mass_match <= selection.mass_match_relative_max &&
         rapidity_match <= selection.rapidity_match_max;
}

// Transform the ordered pi+ pi- state into the selected angular frame
std::vector<M4Vec> TransformPair(const EventKinematics &event,
                                 AngularFrame frame) {
  std::vector<M4Vec> pair = {event.pip, event.pim};
  const M4Vec system = event.pip + event.pim;
  constexpr int direction = 1;
  if (frame == AngularFrame::CM) {
    gra::kinematics::CMframe(pair, system);
  } else if (frame == AngularFrame::HX) {
    gra::kinematics::HXframe(pair, system);
  } else if (frame == AngularFrame::CS) {
    gra::kinematics::CSframe(pair, system, event.beam_plus, event.beam_minus);
  } else if (frame == AngularFrame::AH) {
    gra::kinematics::AHframe(pair, system, event.beam_plus, event.beam_minus);
  } else if (frame == AngularFrame::PG) {
    gra::kinematics::PGframe(pair, system, direction, event.beam_plus,
                             event.beam_minus);
  } else {
    gra::kinematics::GJframe(pair, system, direction,
                             event.beam_plus - event.proton_plus,
                             event.beam_minus - event.proton_minus);
  }
  return pair;
}

// Project one accepted topology into angular and conditional coordinates
ProjectedRecord ProjectEvent(const EventKinematics &event,
                             const HepMCReadConfig &config) {
  const std::vector<M4Vec> pair = TransformPair(event, config.frame);
  const M4Vec system = event.pip + event.pim;
  ProjectedRecord output;
  output.observation.costheta = pair.front().CosTheta();
  output.observation.phi      = math::WrapAngle(pair.front().Phi());
  output.observation.weight = event.weight;
  output.fiducial.central =
      AcceptPion(event.pip, config.cuts) && AcceptPion(event.pim, config.cuts);
  output.fiducial.forward = config.mode == MeasurementMode::Central ||
                            AcceptForward(event, config.cuts);
  output.detector_selected = output.fiducial.Pass(config.mode);
  if (config.mode == MeasurementMode::Tagged) {
    output.detector_selected = output.detector_selected &&
                               AcceptExclusiveClosure(event, config.selection);
  }

  if (config.mode == MeasurementMode::Central) {
    output.observation.z = {system.M(), system.Pt(), system.Rap()};
  } else {
    const auto plus = ForwardCoordinates(event.beam_plus, event.proton_plus);
    const auto minus = ForwardCoordinates(event.beam_minus, event.proton_minus);
    output.observation.z = {system.M(), system.Rap(), plus.second, minus.second,
                            math::WrapAngle(event.proton_plus.Phi() - event.proton_minus.Phi())};
  }
  output.observation.Validate(Coordinates(config.mode).size(),
                              "harmonic::HepMCAnalysisReader projection");
  return output;
}

} // namespace

// Compute a stable printable angular-frame name
std::string ToString(AngularFrame frame) {
  if (frame == AngularFrame::CM) {
    return "CM";
  }
  if (frame == AngularFrame::HX) {
    return "HX";
  }
  if (frame == AngularFrame::CS) {
    return "CS";
  }
  if (frame == AngularFrame::AH) {
    return "AH";
  }
  if (frame == AngularFrame::PG) {
    return "PG";
  }
  return "GJ";
}

// Parse an angular frame without accepting implicit aliases
AngularFrame ParseAngularFrame(const std::string &value) {
  if (value == "CM") {
    return AngularFrame::CM;
  }
  if (value == "HX") {
    return AngularFrame::HX;
  }
  if (value == "CS") {
    return AngularFrame::CS;
  }
  if (value == "AH") {
    return AngularFrame::AH;
  }
  if (value == "PG") {
    return AngularFrame::PG;
  }
  if (value == "GJ") {
    return AngularFrame::GJ;
  }
  throw std::invalid_argument(
      "harmonic::ParseAngularFrame expects CM, HX, CS, AH, PG or GJ");
}

// Test one finite value against the inclusive interval
bool FiducialRange::Contains(double value) const {
  return std::isfinite(value) && value >= min && value <= max;
}

// Validate the finite non-empty interval
void FiducialRange::Validate(const std::string &name) const {
  if (!std::isfinite(min) || !std::isfinite(max) || !(min < max)) {
    throw std::invalid_argument("harmonic::FiducialRange has invalid " + name);
  }
}

// Validate central cuts and tagged-mode forward cuts
void FiducialCuts::Validate(MeasurementMode mode) const {
  pion_eta.Validate("pion eta interval");
  pion_pt.Validate("pion pt interval");
  if (pion_pt.min < 0.0) {
    throw std::invalid_argument(
        "harmonic::FiducialCuts requires non-negative pion pt");
  }
  if (mode == MeasurementMode::Tagged) {
    proton_xi.Validate("proton xi interval");
    proton_abs_t.Validate("proton abs(t) interval");
    if (proton_xi.min < 0.0 || proton_xi.max > 1.0 || proton_abs_t.min < 0.0) {
      throw std::invalid_argument(
          "harmonic::FiducialCuts has invalid tagged-proton bounds");
    }
  }
}

// Validate tagged-mode data-side closure cuts
void ExclusiveSelection::Validate(MeasurementMode mode) const {
  if (mode == MeasurementMode::Central) {
    return;
  }
  if (!std::isfinite(pt_balance_max) || pt_balance_max <= 0.0 ||
      !std::isfinite(mass_match_relative_max) ||
      mass_match_relative_max <= 0.0 || !std::isfinite(rapidity_match_max) ||
      rapidity_match_max <= 0.0) {
    throw std::invalid_argument(
        "harmonic::ExclusiveSelection has invalid tagged closure cuts");
  }
}

// Validate beams, cuts and frame compatibility
void HepMCReadConfig::Validate() const {
  cuts.Validate(mode);
  selection.Validate(mode);
  if (!std::isfinite(sqrt_s) || sqrt_s <= 2.0 * PDG::mp || max_records == 0) {
    throw std::invalid_argument(
        "harmonic::HepMCReadConfig has invalid beam energy or record limit");
  }
  if (mode == MeasurementMode::Central && frame == AngularFrame::GJ) {
    throw std::invalid_argument(
        "harmonic::HepMCReadConfig GJ requires measured forward protons");
  }
}

// Construct an immutable event projector
HepMCAnalysisReader::HepMCAnalysisReader(const HepMCReadConfig &config)
    : config_(config) {
  config_.Validate();
}

// Read generated events and synthesize their detector observations
std::vector<ResponseEvent> HepMCAnalysisReader::ReadModelResponse(
    const std::string &truth_path, const DetectorResponseModel &model) const {
  MHepMCReader reader(truth_path);
  std::vector<ResponseEvent> output;
  std::uint64_t records = 0;
  while (records < config_.max_records) {
    HepMC3::GenEvent event(HepMC3::Units::GEV, HepMC3::Units::MM);
    if (!reader.Read(event)) {
      break;
    }
    const std::uint64_t key = EventKey(truth_path, records);
    ++records;
    const std::optional<EventKinematics> selected = SelectEvent(event, config_);
    if (!selected.has_value()) {
      continue;
    }
    const ProjectedRecord truth = ProjectEvent(*selected, config_);
    ResponseEvent response;
    response.event_key = key;
    response.source = truth_path;
    response.truth = truth.observation;
    response.fiducial = truth.fiducial;
    const auto reconstructed = SimulateDetector(model, *selected, key, config_.seed);
    if (reconstructed.has_value()) {
      const ProjectedRecord reco = ProjectEvent(*reconstructed, config_);
      response.reco = reco.observation;
      response.reco_selected = reco.detector_selected;
    }
    output.push_back(std::move(response));
  }
  if (records == 0 || output.empty()) {
    throw std::invalid_argument(
        "harmonic::HepMCAnalysisReader found no usable truth events in " +
        truth_path);
  }
  return output;
}

// Read generated events with a sparse event-number-matched reco stream
std::vector<ResponseEvent>
HepMCAnalysisReader::ReadPairedResponse(const std::string &truth_path,
                                        const std::string &reco_path) const {
  MHepMCReader truth_reader(truth_path);
  MHepMCReader reco_reader(reco_path);
  bool reco_eof = false;
  std::optional<HepMC3::GenEvent> reco_buffer;
  std::int64_t previous_truth = std::numeric_limits<std::int64_t>::min();
  std::int64_t previous_reco = std::numeric_limits<std::int64_t>::min();
  std::vector<ResponseEvent> output;
  std::uint64_t records = 0;

  while (records < config_.max_records) {
    HepMC3::GenEvent truth_event(HepMC3::Units::GEV, HepMC3::Units::MM);
    if (!truth_reader.Read(truth_event)) {
      break;
    }
    ++records;
    const std::int64_t truth_number = truth_event.event_number();
    if (truth_number <= previous_truth) {
      throw std::invalid_argument("harmonic::HepMCAnalysisReader requires "
                                  "increasing truth event numbers");
    }
    previous_truth = truth_number;

    while (!reco_eof &&
           (!reco_buffer.has_value() ||
            reco_buffer->event_number() < truth_number)) {
      HepMC3::GenEvent next(HepMC3::Units::GEV, HepMC3::Units::MM);
      if (!reco_reader.Read(next)) {
        reco_buffer.reset();
        reco_eof = true;
        break;
      }
      if (next.event_number() <= previous_reco) {
        throw std::invalid_argument("harmonic::HepMCAnalysisReader requires "
                                    "increasing reco event numbers");
      }
      previous_reco = next.event_number();
      reco_buffer = std::move(next);
    }

    const std::optional<EventKinematics> truth_selected =
        SelectEvent(truth_event, config_);
    if (!truth_selected.has_value()) {
      continue;
    }
    const ProjectedRecord truth = ProjectEvent(*truth_selected, config_);
    ResponseEvent response;
    response.event_key = EventKey(truth_path + "|" + reco_path,
                                  static_cast<std::uint64_t>(truth_number));
    response.source = truth_path + "|" + reco_path;
    response.truth = truth.observation;
    response.fiducial = truth.fiducial;

    if (reco_buffer.has_value() &&
        reco_buffer->event_number() == truth_number) {
      // The generated event weight defines the paired response measure
      const std::optional<EventKinematics> reco_selected =
          SelectEvent(*reco_buffer, config_, response.truth.weight);
      if (reco_selected.has_value()) {
        const ProjectedRecord reco = ProjectEvent(*reco_selected, config_);
        response.reco = reco.observation;
        response.reco_selected = reco.detector_selected;
      }
    }
    output.push_back(std::move(response));
  }
  if (records == 0 || output.empty()) {
    throw std::invalid_argument("harmonic::HepMCAnalysisReader found no usable "
                                "paired truth events in " +
                                truth_path);
  }
  return output;
}

// Read selected detector-level events from data or closure MC
std::vector<DataEvent>
HepMCAnalysisReader::ReadData(const std::string &path) const {
  MHepMCReader reader(path);
  std::vector<DataEvent> output;
  std::uint64_t records = 0;
  while (records < config_.max_records) {
    HepMC3::GenEvent event(HepMC3::Units::GEV, HepMC3::Units::MM);
    if (!reader.Read(event)) {
      break;
    }
    const std::uint64_t key = EventKey(path, records);
    ++records;
    const std::optional<EventKinematics> selected = SelectEvent(event, config_);
    if (!selected.has_value()) {
      continue;
    }
    const ProjectedRecord reco = ProjectEvent(*selected, config_);
    if (!reco.detector_selected) {
      continue;
    }
    output.push_back(DataEvent{key, path, reco.observation});
  }
  if (records == 0 || output.empty()) {
    throw std::invalid_argument(
        "harmonic::HepMCAnalysisReader found no selected events in " + path);
  }
  return output;
}

} // namespace harmonic
} // namespace gra
