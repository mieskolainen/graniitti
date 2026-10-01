// GRANIITTI Les Houches Event writer helpers
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <functional>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>

// HepMC3
#include "HepMC3/Attribute.h"
#include "HepMC3/FourVector.h"
#include "HepMC3/GenCrossSection.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenPdfInfo.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/Relatives.h"

// Own
#include "Graniitti/Analysis/MHepMCReader.h"
#include "Graniitti/Particle/MPDG.h"
#include "Graniitti/Nuclear/MUPC.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MLHE.h"

using gra::aux::indices;

namespace gra {
namespace {

struct MLHEColorSlot {
  int *tag = nullptr;
  std::size_t row = 0;
};

struct MLHEDiffSide {
  bool active = false;
  bool excited = false;
  bool fragmented = false;
  bool string_skeleton = false;
  bool hard_remnant_container = false;
  bool nuclear = false;
  int nuclear_a = 0;
  int nuclear_z = 0;
  int nuclear_parent_pdg = 0;
  double neutron_threshold = 0.0;
  HepMC3::ConstGenParticlePtr beam;
  HepMC3::ConstGenParticlePtr system;
  HepMC3::ConstGenParticlePtr leading;
  HepMC3::ConstGenParticlePtr remnant;
  double xi = 0.0;
  double beta = 0.0;
  double t = 0.0;
};

struct MLHEDiffInfo {
  int side_mask = 0;
  bool elastic_cep = false;
  std::pair<int, int> parton_id = {0, 0};
  std::pair<double, double> xhard = {0.0, 0.0};
  std::pair<double, double> beam_energy = {0.0, 0.0};
  MLHEDiffSide side1;
  MLHEDiffSide side2;
};

struct MLHECrossSection {
  double value = 0.0;
  double error = 0.0;
};

// Final cross section and weight convention determined from the complete input
struct MLHEInputRun {
  MLHERunConfig weights;
  MLHECrossSection xs;
};

// Compute one particle momentum in the GeV units required by LHE
HepMC3::FourVector MomentumGeV(const HepMC3::ConstGenParticlePtr &particle) {
  auto momentum = particle->momentum();
  if (!std::isfinite(momentum.px()) || !std::isfinite(momentum.py()) ||
      !std::isfinite(momentum.pz()) || !std::isfinite(momentum.e())) {
    throw std::invalid_argument("MLHE: non-finite particle four-momentum");
  }
  if (particle->parent_event() != nullptr) {
    HepMC3::Units::convert(momentum, particle->parent_event()->momentum_unit(), HepMC3::Units::GEV);
  }
  return momentum;
}

// Compute whether two finite run-level values agree to floating-point accuracy
bool SameRunValue(const double first, const double second) {
  if (!std::isfinite(first) || !std::isfinite(second)) {
    return false;
  }
  constexpr double tolerance = 64.0 * std::numeric_limits<double>::epsilon();
  const double scale = std::max({1.0, std::abs(first), std::abs(second)});
  return std::abs(first - second) <= tolerance * scale;
}

// Read one finite positive HepMC cross section in picobarns
MLHECrossSection EventCrossSection(const HepMC3::GenEvent &ev) {
  const auto cross_section = ev.cross_section();
  if (cross_section == nullptr || cross_section->xsecs().empty()) {
    throw std::invalid_argument(
        "MLHEWriter: event has no HepMC cross section");
  }
  MLHECrossSection output{cross_section->xsec(0),
                          cross_section->xsec_err(0)};
  if (!(output.value > 0.0) || !std::isfinite(output.value) ||
      !(output.error >= 0.0) || !std::isfinite(output.error)) {
    throw std::invalid_argument(
        "MLHEWriter: event has an invalid HepMC cross section");
  }
  return output;
}

// Compute the LHA event-weight strategy for one explicit weight type
int LHEWeightStrategy(const MLHEWeightType type) {
  if (type == MLHEWeightType::Unit) {
    return 3;
  }
  if (type == MLHEWeightType::Weighted) {
    return 4;
  }
  return -4;
}

// Read one integer event attribute with a zero fallback
int ReadEventInt(const HepMC3::GenEvent &ev, const std::string &name) {
  const auto attribute = ev.attribute<HepMC3::IntAttribute>(name);
  return attribute ? attribute->value() : 0;
}

// Read one floating-point event attribute with a fallback
double ReadEventDouble(const HepMC3::GenEvent &ev, const std::string &name,
                       double fallback = 0.0) {
  const auto attribute = ev.attribute<HepMC3::DoubleAttribute>(name);
  return attribute ? attribute->value() : fallback;
}

// Read one standard HepMC color-flow attribute from a particle
int ReadColorFlow(const HepMC3::ConstGenParticlePtr &particle,
                  const std::string &name) {
  const auto attribute = particle->attribute<HepMC3::IntAttribute>(name);
  return attribute ? attribute->value() : 0;
}

// Read one integer particle attribute with a zero fallback
int ReadParticleInt(const HepMC3::ConstGenParticlePtr &particle,
                    const std::string &name) {
  const auto attribute = particle->attribute<HepMC3::IntAttribute>(name);
  return attribute ? attribute->value() : 0;
}

// Read one floating-point particle attribute with a zero fallback
double ReadParticleDouble(const HepMC3::ConstGenParticlePtr &particle,
                          const std::string &name) {
  const auto attribute = particle->attribute<HepMC3::DoubleAttribute>(name);
  return attribute ? attribute->value() : 0.0;
}

// Read one string particle attribute with an empty fallback
std::string ReadParticleString(const HepMC3::ConstGenParticlePtr &particle,
                               const std::string &name) {
  const auto attribute = particle->attribute<HepMC3::StringAttribute>(name);
  return attribute ? attribute->value() : "";
}

// Map GRANIITTI event-record statuses to Les Houches statuses
// [REFERENCE: hep-ph/0109068, ISTUP spacelike propagator and resonance statuses]
int LHEStatus(const HepMC3::ConstGenParticlePtr &particle) {
  const int status = particle->status();
  if (status == gra::PDG::PDG_BEAM) {
    return -1;
  }
  if (status == gra::PDG::PDG_STABLE) {
    return 1;
  }
  if (status == gra::PDG::PDG_INTERMEDIATE) {
    return particle->generated_mass() < 0.0 ? -2 : 2;
  }
  if (status == gra::PDG::PDG_DECAY) {
    return 2;
  }
  return status;
}

// Compute true for internal GRANIITTI container particle ids
bool IsInternalParticle(int pid) {
  return pid == gra::PDG::PDG_system || pid == gra::PDG::PDG_propagator ||
         std::abs(pid) == gra::PDG::PDG_NSTAR || pid == gra::PDG::PDG_fragment;
}

// Compute true for colored diquark remnant ids
bool IsDiquark(int pid) {
  const int apdg = std::abs(pid);
  const int spin = apdg % 10;
  return apdg >= 1000 && apdg < 6000 && (spin == 1 || spin == 3);
}

// Compute true for partons carrying QCD color
bool IsColoredLHEParton(int pid) {
  return pid == gra::PDG::PDG_gluon ||
         (std::abs(pid) >= 1 && std::abs(pid) <= 6) || IsDiquark(pid);
}

// Compute true for spectator particle ids omitted from hard-process LHE views
bool IsHardPdfSpectatorId(int pid) {
  const int apdg = std::abs(pid);
  return apdg == gra::PDG::PDG_p || apdg == gra::PDG::PDG_n || IsDiquark(pid) ||
         pid == gra::PDG::PDG_fragment;
}

// Compute true for particles descending from forward remnant containers
bool HasForwardRemnantParent(const HepMC3::ConstGenParticlePtr &particle) {
  for (const auto &parent : HepMC3::Relatives::ANCESTORS(particle)) {
    if (std::abs(parent->pid()) == gra::PDG::PDG_NSTAR) {
      return true;
    }
  }
  return false;
}

// Compute true for spectator particles omitted from hard-process LHE views
bool IsHardPdfSpectatorParticle(const HepMC3::ConstGenParticlePtr &particle) {
  if (HasForwardRemnantParent(particle)) { return true; }
  const auto parents = HepMC3::Relatives::PARENTS(particle);
  return IsHardPdfSpectatorId(particle->pid()) && parents.size() == 1 &&
         parents.front()->status() == gra::PDG::PDG_BEAM;
}

// Preserve the generated mass, including a negative sign for spacelike virtuality
// [REFERENCE: hep-ph/0109068, PUP(5,I) generated-mass convention]
double MassGeV(const HepMC3::ConstGenParticlePtr &particle) {
  double mass = particle->generated_mass();
  if (particle->parent_event() != nullptr) {
    HepMC3::Units::convert(mass, particle->parent_event()->momentum_unit(), HepMC3::Units::GEV);
  }
  if (!std::isfinite(mass)) { throw std::invalid_argument("MLHE: non-finite particle mass"); }
  return mass;
}

// Convert one HepMC four-momentum into LHE order
std::vector<double> LHEMomentum(const HepMC3::ConstGenParticlePtr &particle) {
  const HepMC3::FourVector pvec = MomentumGeV(particle);
  return {pvec.px(), pvec.py(), pvec.pz(), pvec.e(), MassGeV(particle)};
}

// Build one massless incoming parton momentum from a beam and PDF x
std::vector<double>
IncomingPartonMomentum(const HepMC3::ConstGenParticlePtr &beam, double x) {
  const HepMC3::FourVector pvec = MomentumGeV(beam);
  const double px = x * pvec.px();
  const double py = x * pvec.py();
  const double pz = x * pvec.pz();
  const double e = std::hypot(px, py, pz);
  return {px, py, pz, e, 0.0};
}

// Append one mapped HepMC particle to the LHE row list
void AppendParticleRow(std::vector<MLHERow> &rows,
                       const HepMC3::ConstGenParticlePtr &particle, long pid,
                       int status) {
  MLHERow row;
  row.particle = particle;
  row.pid = pid;
  row.status = status;
  row.colors = {ReadColorFlow(particle, "flow1"),
                ReadColorFlow(particle, "flow2")};
  row.momentum = LHEMomentum(particle);
  rows.push_back(row);
}

// Append one incoming hard parton row from HepMC PDF metadata
void AppendIncomingPartonRow(std::vector<MLHERow> &rows,
                             const HepMC3::ConstGenParticlePtr &beam, long pid,
                             double x) {
  MLHERow row;
  row.pid = pid;
  row.status = -1;
  row.momentum = IncomingPartonMomentum(beam, x);
  rows.push_back(row);
}

// Append one incoming hard parton row with an explicit four-momentum
void AppendIncomingPartonRow(std::vector<MLHERow> &rows,
                             const HepMC3::ConstGenParticlePtr &particle,
                             long pid) {
  MLHERow row;
  row.particle = particle;
  row.pid = pid;
  row.status = -1;
  row.colors = {ReadColorFlow(particle, "flow1"), ReadColorFlow(particle, "flow2")};
  row.momentum = LHEMomentum(particle);
  rows.push_back(row);
}

// Compute the beam particle on one beam side
HepMC3::ConstGenParticlePtr FindBeamParticle(const HepMC3::GenEvent &ev,
                                             bool plus_side) {
  for (const auto &particle : ev.particles()) {
    if (particle->status() != gra::PDG::PDG_BEAM) {
      continue;
    }
    if (plus_side && MomentumGeV(particle).pz() >= 0.0) {
      return particle;
    }
    if (!plus_side && MomentumGeV(particle).pz() < 0.0) {
      return particle;
    }
  }
  return nullptr;
}

// Compute the light-cone energy component for one beam side
// p+ = E+pz or p- = E-pz
double LightConeEnergy(const HepMC3::FourVector &p, bool plus_side) {
  return plus_side ? p.e() + p.pz() : p.e() - p.pz();
}

// Compute the invariant square of a HepMC four-vector
// p^2 = E^2-px^2-py^2-pz^2
double Mass2(const HepMC3::FourVector &p) {
  return p.e() * p.e() - p.px() * p.px() - p.py() * p.py() - p.pz() * p.pz();
}

// Compute the squared momentum transfer from an incoming beam to a leading proton
// t = -|(p_beam-p_leading)^2|
double MomentumTransferT(const HepMC3::ConstGenParticlePtr &beam,
                         const HepMC3::ConstGenParticlePtr &leading) {
  const HepMC3::FourVector p = MomentumGeV(beam);
  const HepMC3::FourVector q = MomentumGeV(leading);
  return -std::abs(Mass2(HepMC3::FourVector(p.px() - q.px(), p.py() - q.py(),
                                            p.pz() - q.pz(), p.e() - q.e())));
}

// Compute the fractional proton momentum loss from light-cone momenta
// xi = 1-p_leading^+/p_beam^+ or 1-p_leading^-/p_beam^-
double ReconstructXi(const HepMC3::ConstGenParticlePtr &beam,
                     const HepMC3::ConstGenParticlePtr &leading,
                     bool plus_side) {
  const double denom = LightConeEnergy(MomentumGeV(beam), plus_side);
  if (!(denom > 0.0)) {
    return 0.0;
  }
  return 1.0 - LightConeEnergy(MomentumGeV(leading), plus_side) / denom;
}

// Select the non-leading forward remnant child when available
HepMC3::ConstGenParticlePtr
FindForwardRemnant(const HepMC3::ConstGenVertexPtr &vertex,
                   const HepMC3::ConstGenParticlePtr &leading) {
  for (const auto &child : vertex->particles_out()) {
    if (child == leading) {
      continue;
    }
    return child;
  }
  return nullptr;
}

// Fill one diffractive-side metadata object from an N-star forward branch
void FillDiffSide(MLHEDiffSide &side,
                  const HepMC3::ConstGenParticlePtr &forward,
                  const HepMC3::ConstGenParticlePtr &beam, double xhard,
                  bool plus_side) {
  if (beam == nullptr || forward == nullptr) {
    return;
  }

  side.active = true;
  side.excited = true;
  side.beam = beam;
  side.system = forward;
  side.xi = ReconstructXi(beam, forward, plus_side);
  side.t = MomentumTransferT(beam, forward);
  side.beta = (side.xi > 0.0) ? xhard / side.xi : 0.0;
  if (forward->end_vertex() == nullptr) {
    return;
  }
  side.fragmented = true;

  HepMC3::ConstGenParticlePtr leading;
  for (const auto &child : forward->end_vertex()->particles_out()) {
    if (child->pid() == beam->pid()) {
      leading = child;
      break;
    }
  }
  if (leading == nullptr) {
    int colored = 0;
    std::map<int, int> balance;
    for (const auto &child : forward->end_vertex()->particles_out()) {
      if (!IsColoredLHEParton(child->pid())) {
        continue;
      }
      ++colored;
      const int flow1 = ReadColorFlow(child, "flow1");
      const int flow2 = ReadColorFlow(child, "flow2");
      if (flow1 != 0) {
        ++balance[flow1];
      }
      if (flow2 != 0) {
        --balance[flow2];
      }
    }
    side.string_skeleton = colored >= 2 && !balance.empty();
    for (const auto &[tag, value] : balance) {
      if (tag == 0 || value != 0) {
        side.string_skeleton = false;
      }
    }
    side.hard_remnant_container = colored > 0 && !side.string_skeleton;
    return;
  }

  side.leading = leading;
  side.remnant = FindForwardRemnant(forward->end_vertex(), leading);
  side.xi = ReconstructXi(beam, leading, plus_side);
  side.t = MomentumTransferT(beam, leading);
  side.beta = (side.xi > 0.0) ? xhard / side.xi : 0.0;
}

// Fill one electromagnetic nuclear-breakup side from its forward container
void FillNuclearDiffSide(MLHEDiffSide &side,
                         const HepMC3::ConstGenParticlePtr &forward,
                         const HepMC3::ConstGenParticlePtr &beam, double xhard,
                         bool plus_side) {
  if (beam == nullptr || forward == nullptr ||
      (ReadParticleString(forward, "graniitti_upc_neutrons").empty() ||
       nuclear::ParseNeutronSelection(ReadParticleString(forward, "graniitti_upc_neutrons")).Accept(0)) ||
      ReadParticleInt(forward, "graniitti_emd_realized") != 0) {
    return;
  }
  side.active = true;
  side.excited = true;
  side.nuclear = true;
  side.beam = beam;
  side.system = forward;
  side.nuclear_a = ReadParticleInt(forward, "graniitti_upc_a");
  side.nuclear_z = ReadParticleInt(forward, "graniitti_upc_z");
  side.nuclear_parent_pdg =
      ReadParticleInt(forward, "graniitti_upc_parent_pdg");
  side.neutron_threshold =
      ReadParticleDouble(forward, "graniitti_upc_neutron_threshold");
  side.xi = ReconstructXi(beam, forward, plus_side);
  side.t = MomentumTransferT(beam, forward);
  side.beta = side.xi > 0.0 ? xhard / side.xi : 0.0;
}

// Find the intact elastic beam daughter on one CEP side
HepMC3::ConstGenParticlePtr
FindElasticLeading(const HepMC3::ConstGenParticlePtr &beam) {
  if (beam == nullptr || beam->end_vertex() == nullptr) {
    return nullptr;
  }
  for (const auto &child : beam->end_vertex()->particles_out()) {
    if (child->pid() == beam->pid() &&
        child->status() == gra::PDG::PDG_STABLE) {
      return child;
    }
  }
  return nullptr;
}

// Fill one elastic CEP side without inventing a colored beam remnant
void FillElasticDiffSide(MLHEDiffSide &side,
                         const HepMC3::ConstGenParticlePtr &beam, double xhard,
                         bool plus_side) {
  const HepMC3::ConstGenParticlePtr leading = FindElasticLeading(beam);
  if (beam == nullptr || leading == nullptr) {
    return;
  }

  side.active = true;
  side.beam = beam;
  side.leading = leading;
  side.xi = ReconstructXi(beam, leading, plus_side);
  side.t = MomentumTransferT(beam, leading);
  side.beta = side.xi > 0.0 && xhard > 0.0 ? 1.0 : 0.0;
}

// Collect diffractive-side metadata from the full HepMC event record
MLHEDiffInfo BuildDiffInfo(const HepMC3::GenEvent &ev) {
  MLHEDiffInfo info;
  const int hard_sides = ReadEventInt(ev, "graniitti_hard_diffraction");
  info.elastic_cep = ReadEventInt(ev, "graniitti_elastic_colored_cep") != 0;
  const auto pdf = ev.pdf_info();
  // Exact EPA metadata must not replace the physical pp beams by effective
  // photons
  bool hard_pdf_view = false;
  if (pdf != nullptr && pdf->is_valid()) {
    info.parton_id = {pdf->parton_id[0], pdf->parton_id[1]};
    info.xhard = {pdf->x[0], pdf->x[1]};
    hard_pdf_view =
        pdf->parton_id[0] != 0 && pdf->parton_id[1] != 0 && pdf->x[0] > 0.0 &&
        pdf->x[1] > 0.0 && !info.elastic_cep &&
        ReadEventInt(ev, "graniitti_exact_forward_photon_kinematics") == 0;
  }

  const auto beam_plus = FindBeamParticle(ev, true);
  const auto beam_minus = FindBeamParticle(ev, false);
  if (beam_plus != nullptr) {
    info.beam_energy.first = MomentumGeV(beam_plus).e();
  }
  if (beam_minus != nullptr) {
    info.beam_energy.second = MomentumGeV(beam_minus).e();
  }

  if (!hard_pdf_view || info.elastic_cep) {
    FillElasticDiffSide(info.side1, beam_plus, info.xhard.first, true);
    FillElasticDiffSide(info.side2, beam_minus, info.xhard.second, false);
  }

  for (const auto &particle : ev.particles()) {
    if (std::abs(particle->pid()) != gra::PDG::PDG_NSTAR) {
      continue;
    }
    MLHEDiffSide candidate;
    if (MomentumGeV(particle).pz() >= 0.0) {
      if (hard_sides != 0 && (hard_sides & 1) == 0) { continue; }
      FillDiffSide(candidate, particle, beam_plus, info.xhard.first, true);
      if (!hard_pdf_view || !candidate.hard_remnant_container) {
        info.side1 = candidate;
      }
    } else {
      if (hard_sides != 0 && (hard_sides & 2) == 0) { continue; }
      FillDiffSide(candidate, particle, beam_minus, info.xhard.second, false);
      if (!hard_pdf_view || !candidate.hard_remnant_container) {
        info.side2 = candidate;
      }
    }
  }

  for (const auto &particle : ev.particles()) {
    if (particle->pid() != gra::PDG::PDG_fragment ||
        ReadParticleInt(particle, "graniitti_upc_parent_pdg") == 0) {
      continue;
    }
    MLHEDiffSide candidate;
    if (MomentumGeV(particle).pz() >= 0.0) {
      FillNuclearDiffSide(candidate, particle, beam_plus, info.xhard.first,
                          true);
      if (candidate.active) {
        info.side1 = candidate;
      }
    } else {
      FillNuclearDiffSide(candidate, particle, beam_minus, info.xhard.second,
                          false);
      if (candidate.active) {
        info.side2 = candidate;
      }
    }
  }

  if (info.side1.active) {
    info.side_mask |= 1;
  }
  if (info.side2.active) {
    info.side_mask |= 2;
  }
  return info;
}

// Append one particle four-vector to a GRANIITTI LHE metadata line
void AppendParticleMetadata(std::ostringstream &out, const std::string &prefix,
                            const HepMC3::ConstGenParticlePtr &particle) {
  if (particle == nullptr) {
    out << ' ' << prefix << "_pdg=0";
    return;
  }
  const HepMC3::FourVector p = MomentumGeV(particle);
  out << ' ' << prefix << "_pdg=" << particle->pid() << ' ' << prefix
      << "_px=" << p.px() << ' ' << prefix << "_py=" << p.py() << ' ' << prefix
      << "_pz=" << p.pz() << ' ' << prefix << "_e=" << p.e() << ' ' << prefix
      << "_m=" << MassGeV(particle);
}

// Append one diffractive-side block to a GRANIITTI LHE metadata line
void AppendDiffSideMetadata(std::ostringstream &out, const std::string &prefix,
                            const MLHEDiffSide &side) {
  out << ' ' << prefix << "_active=" << (side.active ? 1 : 0);
  if (!side.active) {
    return;
  }
  out << ' ' << prefix << "_excited=" << (side.excited ? 1 : 0) << ' ' << prefix
      << "_fragmented=" << (side.fragmented ? 1 : 0) << ' ' << prefix
      << "_string=" << (side.string_skeleton ? 1 : 0) << ' ' << prefix
      << "_nuclear=" << (side.nuclear ? 1 : 0) << ' ' << prefix
      << "_nuclear_a=" << side.nuclear_a << ' ' << prefix
      << "_nuclear_z=" << side.nuclear_z << ' ' << prefix
      << "_nuclear_parent_pdg=" << side.nuclear_parent_pdg << ' ' << prefix
      << "_neutron_threshold=" << side.neutron_threshold << ' ' << prefix
      << "_xi=" << side.xi << ' ' << prefix << "_beta=" << side.beta << ' '
      << prefix << "_t=" << side.t;
  AppendParticleMetadata(out, prefix + "_system", side.system);
  AppendParticleMetadata(out, prefix + "_lead", side.leading);
  AppendParticleMetadata(out, prefix + "_rem", side.remnant);
}

// Compute true when a parton-level LHE view can be built
bool HasHardPdfInfo(const HepMC3::GenEvent &ev) {
  if (ReadEventInt(ev, "graniitti_elastic_colored_cep") != 0) {
    return false;
  }
  // Exact EPA events are complete pp records even when diagnostic PDF
  // information exists
  if (ReadEventInt(ev, "graniitti_exact_forward_photon_kinematics") != 0) {
    return false;
  }
  const auto pdf = ev.pdf_info();
  if (pdf == nullptr || !pdf->is_valid() || pdf->parton_id[0] == 0 || pdf->parton_id[1] == 0) { return false; }
  for (const double x : pdf->x) {
    if (!std::isfinite(x) || !(x > 0.0) || x > 1.0) {
      throw std::invalid_argument("MLHE: incoming PDF momentum fractions must lie in (0,1]");
    }
  }
  return FindBeamParticle(ev, true) != nullptr && FindBeamParticle(ev, false) != nullptr;
}

// Compute true when one propagator enters a two-propagator fusion vertex
bool IsCentralIncomingPropagator(const HepMC3::ConstGenParticlePtr &particle) {
  if (particle->pid() != gra::PDG::PDG_propagator ||
      particle->end_vertex() == nullptr) {
    return false;
  }
  const auto &incoming = particle->end_vertex()->particles_in();
  return incoming.size() == 2 && std::all_of(incoming.begin(), incoming.end(), [](const auto &parent) {
    return parent->pid() == gra::PDG::PDG_propagator;
  });
}

// Compute the two incoming hard propagators ordered by beam side
std::pair<HepMC3::ConstGenParticlePtr, HepMC3::ConstGenParticlePtr>
FindHardPropagators(const HepMC3::GenEvent &ev) {
  std::vector<HepMC3::ConstGenParticlePtr> candidates;
  for (const auto &particle : ev.particles()) {
    if (IsCentralIncomingPropagator(particle)) {
      candidates.push_back(particle);
    }
  }
  if (candidates.size() != 2) {
    return {nullptr, nullptr};
  }

  std::sort(candidates.begin(), candidates.end(),
            [](const auto &a, const auto &b) {
              return MomentumGeV(a).pz() > MomentumGeV(b).pz();
            });
  return {candidates[0], candidates[1]};
}

// Fill row lookup from original HepMC particle ids
std::map<int, int> BuildRowIndex(const std::vector<MLHERow> &rows) {
  std::map<int, int> row_index;
  for (const auto &i : indices(rows)) {
    if (rows[i].particle != nullptr) {
      row_index[rows[i].particle->id()] = static_cast<int>(i) + 1;
    }
  }
  return row_index;
}

// Compute true if the first two rows are incoming particles
bool HasIncomingPair(const std::vector<MLHERow> &rows) {
  return rows.size() >= 2 && rows[0].status == -1 && rows[1].status == -1;
}

// Compute the largest explicit color-flow tag in the selected rows
int MaxColorTag(const std::vector<MLHERow> &rows) {
  int max_tag = 500;
  for (const auto &row : rows) {
    max_tag = std::max(max_tag, std::abs(row.colors.first));
    max_tag = std::max(max_tag, std::abs(row.colors.second));
  }
  return max_tag;
}

// Append effective LHE color slots with incoming particle flow reversed
void AppendEffectiveColorSlots(MLHERow &row, std::size_t row_index,
                               std::vector<MLHEColorSlot> &color_slots,
                               std::vector<MLHEColorSlot> &anticolor_slots) {
  if (!IsColoredLHEParton(row.pid)) {
    return;
  }

  row.colors = {0, 0};
  const bool incoming = row.status == -1;

  if (row.pid == gra::PDG::PDG_gluon) {
    if (incoming) {
      anticolor_slots.push_back({&row.colors.first, row_index});
      color_slots.push_back({&row.colors.second, row_index});
    } else {
      color_slots.push_back({&row.colors.first, row_index});
      anticolor_slots.push_back({&row.colors.second, row_index});
    }
  } else if ((row.pid > 0) != IsDiquark(row.pid)) {
    if (incoming) {
      anticolor_slots.push_back({&row.colors.first, row_index});
    } else {
      color_slots.push_back({&row.colors.first, row_index});
    }
  } else {
    if (incoming) {
      color_slots.push_back({&row.colors.second, row_index});
    } else {
      anticolor_slots.push_back({&row.colors.second, row_index});
    }
  }
}

// Preserve supplied color lines and reject incomplete colored particle records
// [REFERENCE: hep-ph/0109068, ICOLUP color flow in physical time order]
bool HasHardProcessColorFlow(const std::vector<MLHERow> &rows) {
  const bool supplied = std::any_of(rows.begin(), rows.end(), [](const auto &row) {
    return row.colors.first != 0 || row.colors.second != 0;
  });
  if (!supplied) { return false; }
  for (const auto &row : rows) {
    if (!IsColoredLHEParton(row.pid)) { continue; }
    const bool gluon = row.pid == gra::PDG::PDG_gluon;
    const bool color = gluon || ((row.pid > 0) != IsDiquark(row.pid));
    const bool anticolor = gluon || !color;
    if ((color ? row.colors.first <= 0 : row.colors.first != 0) ||
        (anticolor ? row.colors.second <= 0 : row.colors.second != 0)) {
      throw std::invalid_argument("MLHE: incomplete supplied hard-process color flow");
    }
  }
  return true;
}

// Reconstruct a unique color topology when the event has no supplied color lines
void AssignHardProcessColorFlow(std::vector<MLHERow> &rows) {
  if (HasHardProcessColorFlow(rows)) { return; }
  std::vector<MLHEColorSlot> color_slots;
  std::vector<MLHEColorSlot> anticolor_slots;

  for (const auto &i : indices(rows)) {
    AppendEffectiveColorSlots(rows[i], i, color_slots, anticolor_slots);
  }

  if (color_slots.size() != anticolor_slots.size()) {
    throw std::invalid_argument("MLHE: hard-process color representations do not balance");
  }
  // Three or more color lines admit several pairings, as do two quark lines
  const bool gluon = std::any_of(rows.begin(), rows.end(), [](const auto &row) {
    return row.pid == gra::PDG::PDG_gluon;
  });
  if (color_slots.size() > 2 || (color_slots.size() == 2 && !gluon)) {
    throw std::invalid_argument("MLHE: ambiguous hard-process color flow requires supplied color lines");
  }

  std::vector<int> anticolor_match(anticolor_slots.size(), -1);
  std::function<bool(std::size_t, std::vector<bool> &)> match_color =
      [&](std::size_t color_index, std::vector<bool> &seen) {
        for (const auto &anti_index : indices(anticolor_slots)) {
          if (seen[anti_index] ||
              color_slots[color_index].row == anticolor_slots[anti_index].row) {
            continue;
          }
          seen[anti_index] = true;
          if (anticolor_match[anti_index] < 0 ||
              match_color(static_cast<std::size_t>(anticolor_match[anti_index]),
                          seen)) {
            anticolor_match[anti_index] = static_cast<int>(color_index);
            return true;
          }
        }
        return false;
      };

  for (const auto &color_index : indices(color_slots)) {
    std::vector<bool> seen(anticolor_slots.size(), false);
    if (!match_color(color_index, seen)) {
      throw std::invalid_argument("MLHE: hard-process color lines cannot be connected");
    }
  }

  const int first_tag = MaxColorTag(rows) + 1;
  for (const auto &anti_index : indices(anticolor_slots)) {
    const std::size_t color_index =
        static_cast<std::size_t>(anticolor_match[anti_index]);
    const int tag = first_tag + static_cast<int>(color_index);
    *color_slots[color_index].tag = tag;
    *anticolor_slots[anti_index].tag = tag;
  }
}

// Assign mother indices after all rows have been selected
void AssignMotherRows(std::vector<MLHERow> &rows) {
  const std::map<int, int> row_index = BuildRowIndex(rows);

  for (auto &row : rows) {
    if (row.status == -1 || row.particle == nullptr) {
      continue;
    }

    std::vector<int> mothers;
    for (const auto &parent : HepMC3::Relatives::PARENTS(row.particle)) {
      const auto found = row_index.find(parent->id());
      if (found != row_index.end()) {
        mothers.push_back(found->second);
      }
    }
    std::sort(mothers.begin(), mothers.end());
    mothers.erase(std::unique(mothers.begin(), mothers.end()), mothers.end());

    if (mothers.empty() && HasIncomingPair(rows)) {
      row.mothers = {1, 2};
    } else if (mothers.size() == 1) {
      row.mothers = {mothers.front(), 0};
    } else if (mothers.size() >= 2) {
      row.mothers = {mothers.front(), mothers.back()};
    }
  }
}

// Compute the event weight used in one Les Houches event
double EventWeight(const HepMC3::GenEvent &ev) {
  return ev.weights().empty() ? 1.0 : ev.weights()[0];
}

// Convert one HepMC event weight into the configured LHA convention
double LHEEventWeight(const HepMC3::GenEvent &ev,
                      const MLHERunConfig &config) {
  const double input = EventWeight(ev);
  if (!std::isfinite(input)) {
    throw std::invalid_argument("MLHEWriter: event weight is not finite");
  }
  if (config.weight_type == MLHEWeightType::Unit) {
    if (!SameRunValue(input, 1.0)) {
      throw std::invalid_argument(
          "MLHEWriter: IDWTUP=3 requires positive unit event weights");
    }
    return 1.0;
  }
  if (config.weight_type == MLHEWeightType::Weighted && input < 0.0) {
    throw std::invalid_argument(
        "MLHEWriter: IDWTUP=4 does not accept negative event weights");
  }
  const double output = input * config.weight_scale;
  if (!std::isfinite(output)) {
    throw std::invalid_argument(
        "MLHEWriter: converted event weight is not finite");
  }
  return output;
}

// Infer one normalized LHA weight convention from a complete HepMC file
MLHEInputRun InspectHepMCWeights(const std::string &inputfile) {
  MHepMCReader input(inputfile);
  HepMC3::GenEvent ev(HepMC3::Units::GEV, HepMC3::Units::MM);
  std::size_t events = 0;
  long double weight_sum = 0.0L;
  bool unit = true;
  bool negative = false;
  std::shared_ptr<HepMC3::GenCrossSection> final_cross_section;
  while (input.Read(ev)) {
    const double weight = EventWeight(ev);
    if (!std::isfinite(weight)) {
      throw std::invalid_argument(
          "ConvertHepMC3ToLHE: input event weight is not finite");
    }
    if (ev.cross_section() != nullptr) {
      final_cross_section = std::make_shared<HepMC3::GenCrossSection>(*ev.cross_section());
    }
    unit = unit && SameRunValue(weight, 1.0);
    negative = negative || weight < 0.0;
    weight_sum += static_cast<long double>(weight);
    ++events;
    ev.clear();
  }

  if (events == 0) {
    throw std::invalid_argument("ConvertHepMC3ToLHE: input contains no readable events");
  }
  HepMC3::GenEvent last;
  if (final_cross_section != nullptr) { last.set_cross_section(final_cross_section); }
  MLHEInputRun run;
  run.xs = EventCrossSection(last);
  if (unit) { return run; }
  auto &config = run.weights;
  const long double mean = weight_sum / static_cast<long double>(events);
  if (!(mean > 0.0L) || !std::isfinite(mean)) {
    throw std::invalid_argument(
        "ConvertHepMC3ToLHE: weighted input has a nonpositive mean weight");
  }
  config.weight_type = negative ? MLHEWeightType::SignedWeighted
                                : MLHEWeightType::Weighted;
  config.weight_scale = run.xs.value / static_cast<double>(mean);
  if (!(config.weight_scale > 0.0) || !std::isfinite(config.weight_scale)) {
    throw std::invalid_argument(
        "ConvertHepMC3ToLHE: weighted input normalization is invalid");
  }
  return run;
}

// Compute the PDF factorization scale
double EventMuF(const HepMC3::GenEvent &ev) {
  const double mu_f = ReadEventDouble(ev, "graniitti_mu_f");
  if (std::isfinite(mu_f) && mu_f > 0.0) {
    return mu_f;
  }
  const auto pdf = ev.pdf_info();
  if (pdf == nullptr || !pdf->is_valid() || !std::isfinite(pdf->scale) || !(pdf->scale > 0.0)) {
    return -1.0;
  }
  // HepMC PDF scales are always in GeV, independently of event units
  return pdf->scale;
}

// Compute the QCD renormalization scale
double EventMuR(const HepMC3::GenEvent &ev) {
  const double mu_r = ReadEventDouble(ev, "graniitti_mu_r");
  return (std::isfinite(mu_r) && mu_r > 0.0) ? mu_r : EventMuF(ev);
}

// Compute the LHE shower starting or veto scale
double EventScalup(const HepMC3::GenEvent &ev) {
  const double scalup = ReadEventDouble(ev, "graniitti_scalup");
  return (std::isfinite(scalup) && scalup > 0.0) ? scalup : EventMuF(ev);
}

// Compute the event-local strong coupling when published by the hard process
double EventAlphaQCD(const HepMC3::GenEvent &ev) {
  const double alpha_qcd = ReadEventDouble(ev, "graniitti_alpha_qcd");
  return (std::isfinite(alpha_qcd) && alpha_qcd > 0.0) ? alpha_qcd : -1.0;
}

// Compute the incoming PDF weights from HepMC PDF information
std::pair<double, double> EventPdfWeights(const HepMC3::GenEvent &ev) {
  const auto pdf = ev.pdf_info();
  return (pdf != nullptr && pdf->is_valid())
             ? std::pair<double, double>(pdf->xf[0], pdf->xf[1])
             : std::pair<double, double>(1.0, 1.0);
}

// Compute the sampled proper decay length from the recorded flight time
// [REFERENCE: hep-ph/0109068, VTIMUP invariant lifetime c tau in mm]
double LifetimeMM(const MLHERow &row) {
  const auto &particle = row.particle;
  if (!particle || !particle->production_vertex() || !particle->end_vertex() ||
      !(row.momentum[4] > 0.0) || !(row.momentum[3] > 0.0)) { return 0.0; }
  const double time = particle->end_vertex()->position().t() - particle->production_vertex()->position().t();
  double length = time * (row.momentum[4] / row.momentum[3]);
  if (particle->parent_event() != nullptr) {
    HepMC3::Units::convert(length, particle->parent_event()->length_unit(), HepMC3::Units::MM);
  }
  if (!std::isfinite(length) || length < 0.0) {
    throw std::invalid_argument("MLHE: invalid particle proper decay length");
  }
  return length;
}

// Convert one row list into the LHE HEPEUP event block
LHEF::HEPEUP BuildHEPEUP(const HepMC3::GenEvent &ev,
                         const MLHERunConfig &config) {
  const std::vector<MLHERow> rows = BuildLHERows(ev);

  LHEF::HEPEUP hepeup;
  hepeup.resize(static_cast<int>(rows.size()));

  for (const auto &i : indices(rows)) {
    hepeup.IDUP[i] = rows[i].pid;
    hepeup.ISTUP[i] = rows[i].status;
    hepeup.MOTHUP[i] = rows[i].mothers;
    hepeup.ICOLUP[i] = rows[i].colors;
    hepeup.PUP[i] = rows[i].momentum;
    hepeup.VTIMUP[i] = LifetimeMM(rows[i]);
    hepeup.SPINUP[i] = 9.0;
  }

  hepeup.NUP = static_cast<int>(rows.size());
  hepeup.IDPRUP = 1;
  hepeup.XWGTUP = LHEEventWeight(ev, config);
  hepeup.XPDWUP = EventPdfWeights(ev);
  hepeup.SCALUP = EventScalup(ev);
  hepeup.scales = LHEF::Scales(hepeup.SCALUP, hepeup.NUP);
  hepeup.scales.muf = EventMuF(ev);
  hepeup.scales.mur = EventMuR(ev);
  hepeup.scales.mups = hepeup.SCALUP;
  const double alpha_qed = ReadEventDouble(ev, "graniitti_alpha_qed");
  hepeup.AQEDUP = std::isfinite(alpha_qed) && alpha_qed > 0.0 ? alpha_qed : -1.0;
  hepeup.AQCDUP = EventAlphaQCD(ev);
  return hepeup;
}

// Build a separate hard-process view with the selected MG5 colors for Pythia ISR
std::vector<MLHERow> HardShowerRows(const HepMC3::GenEvent &ev) {
  const auto hard = FindHardPropagators(ev);
  const auto pdf = ev.pdf_info();
  if (!pdf || !hard.first || !hard.second) {
    throw std::invalid_argument("MLHE: missing hard-diffraction propagators or PDF information");
  }
  std::vector<MLHERow> rows;
  AppendIncomingPartonRow(rows, hard.first, pdf->parton_id[0]);
  AppendIncomingPartonRow(rows, hard.second, pdf->parton_id[1]);
  rows[0].colors = {ReadEventInt(ev, "graniitti_hard1_flow1"), ReadEventInt(ev, "graniitti_hard1_flow2")};
  rows[1].colors = {ReadEventInt(ev, "graniitti_hard2_flow1"), ReadEventInt(ev, "graniitti_hard2_flow2")};
  for (const auto &particle : ev.particles()) {
    if (particle->status() == gra::PDG::PDG_BEAM || IsInternalParticle(particle->pid()) ||
        IsHardPdfSpectatorParticle(particle)) { continue; }
    AppendParticleRow(rows, particle, particle->pid(), LHEStatus(particle));
  }
  AssignMotherRows(rows);
  return rows;
}

} // namespace

// Build the generic Les Houches particle-row view of one HepMC event
std::vector<MLHERow> BuildLHERows(const HepMC3::GenEvent &ev) {
  std::vector<MLHERow> rows;

  const auto beam_plus = FindBeamParticle(ev, true);
  const auto beam_minus = FindBeamParticle(ev, false);
  // Keep generated hard-diffraction remnants and their exact MG5 color connections
  const bool hard_pdf_view = HasHardPdfInfo(ev) && ReadEventInt(ev, "graniitti_hard_diffraction") == 0;

  if (hard_pdf_view) {
    const auto pdf = ev.pdf_info();
    const auto hard = FindHardPropagators(ev);
    if (hard.first != nullptr && hard.second != nullptr) {
      AppendIncomingPartonRow(rows, hard.first, pdf->parton_id[0]);
      AppendIncomingPartonRow(rows, hard.second, pdf->parton_id[1]);
    } else {
      AppendIncomingPartonRow(rows, beam_plus, pdf->parton_id[0], pdf->x[0]);
      AppendIncomingPartonRow(rows, beam_minus, pdf->parton_id[1], pdf->x[1]);
    }
  } else {
    if (beam_plus != nullptr) {
      AppendParticleRow(rows, beam_plus, beam_plus->pid(), -1);
    }
    if (beam_minus != nullptr) {
      AppendParticleRow(rows, beam_minus, beam_minus->pid(), -1);
    }
  }

  for (const auto &particle : ev.particles()) {
    if (particle->status() == gra::PDG::PDG_BEAM) {
      continue;
    }
    if (IsInternalParticle(particle->pid())) {
      continue;
    }
    if (hard_pdf_view && IsHardPdfSpectatorParticle(particle)) {
      continue;
    }
    AppendParticleRow(rows, particle, particle->pid(),
                      LHEStatus(particle));
  }

  if (hard_pdf_view) {
    AssignHardProcessColorFlow(rows);
  }
  AssignMotherRows(rows);
  return rows;
}

// Build the optional GRANIITTI diffraction metadata LHE comment for one event
std::string BuildLHEDiffractionComment(const HepMC3::GenEvent &ev) {
  const MLHEDiffInfo info = BuildDiffInfo(ev);
  if (info.side_mask == 0) {
    return "";
  }

  long accepted_events = 0;
  long attempted_events = 0;
  const auto cross_section = ev.cross_section();
  if (cross_section != nullptr) {
    accepted_events = cross_section->get_accepted_events();
    attempted_events = cross_section->get_attempted_events();
  }

  std::ostringstream out;
  out << std::scientific << std::setprecision(17);
  out << "graniitti_diff version=1 elastic_cep=" << (info.elastic_cep ? 1 : 0)
      << " accepted_events=" << accepted_events
      << " attempted_events=" << attempted_events
      << " side_mask=" << info.side_mask << " id1=" << info.parton_id.first
      << " id2=" << info.parton_id.second << " xhard1=" << info.xhard.first
      << " xhard2=" << info.xhard.second
      << " beam1_e=" << info.beam_energy.first
      << " beam2_e=" << info.beam_energy.second;
  AppendDiffSideMetadata(out, "side1", info.side1);
  AppendDiffSideMetadata(out, "side2", info.side2);
  out << '\n';
  if (ReadEventInt(ev, "graniitti_hard_diffraction") != 0) {
    for (const auto &row : HardShowerRows(ev)) {
      out << "# graniitti_hard " << row.pid << ' ' << row.status << ' '
          << row.mothers.first << ' ' << row.mothers.second << ' '
          << row.colors.first << ' ' << row.colors.second;
      for (const auto &value : row.momentum) { out << ' ' << value; }
      out << " 0 9\n";
    }
  }
  return out.str();
}

// Open one Les Houches Event output file
MLHEWriter::MLHEWriter(const std::string &outputfile,
                       const MLHERunConfig &config)
    : config(config) {
  if (!(config.weight_scale > 0.0) ||
      !std::isfinite(config.weight_scale)) {
    throw std::invalid_argument(
        "MLHEWriter: event weight scale must be finite and positive");
  }
  if (config.weight_type == MLHEWeightType::Unit &&
      !SameRunValue(config.weight_scale, 1.0)) {
    throw std::invalid_argument(
        "MLHEWriter: unit events cannot have a weight conversion scale");
  }
  gra::aux::CreateDirectory(std::filesystem::path(outputfile).parent_path().string());
  output.open(outputfile);
  if (!output.is_open()) {
    throw std::runtime_error("MLHEWriter: cannot open " + outputfile);
  }
  writer = std::make_unique<LHEF::Writer>(output);
}

// Flush and close the Les Houches Event output file
MLHEWriter::~MLHEWriter() = default;

// Initialize the Les Houches run block from the first event
void MLHEWriter::InitializeRunBlock(const HepMC3::GenEvent &ev) {
  if (run_block_initialized) {
    return;
  }

  const auto beam_plus = FindBeamParticle(ev, true);
  const auto beam_minus = FindBeamParticle(ev, false);

  if (beam_plus == nullptr || beam_minus == nullptr) {
    throw std::invalid_argument(
        "MLHEWriter: event has no complete incoming beam pair");
  }
  const std::pair<double, double> energy = {MomentumGeV(beam_plus).e(),
                                             MomentumGeV(beam_minus).e()};
  if (!(energy.first > 0.0) || !(energy.second > 0.0) ||
      !std::isfinite(energy.first) || !std::isfinite(energy.second)) {
    throw std::invalid_argument("MLHEWriter: invalid incoming beam energy");
  }
  writer->heprup.IDBMUP = {beam_plus->pid(), beam_minus->pid()};
  writer->heprup.EBMUP = energy;
  writer->heprup.PDFGUP = {-1, -1};
  writer->heprup.PDFSUP = {-1, -1};
  writer->heprup.IDWTUP = LHEWeightStrategy(config.weight_type);
  writer->heprup.NPRUP = 1;

  const MLHECrossSection cross_section = EventCrossSection(ev);
  writer->heprup.XSECUP = {cross_section.value};
  writer->heprup.XERRUP = {cross_section.error};
  writer->heprup.XMAXUP = {
      config.weight_type == MLHEWeightType::Unit ? 1.0 : 0.0};
  writer->heprup.LPRUP = {1};

  writer->init();
  run_block_initialized = true;
}

// Require later events to retain the initialized run identity
void MLHEWriter::ValidateRunBlock(const HepMC3::GenEvent &ev) const {
  const auto beam_plus = FindBeamParticle(ev, true);
  const auto beam_minus = FindBeamParticle(ev, false);
  if (beam_plus == nullptr || beam_minus == nullptr ||
      beam_plus->pid() != writer->heprup.IDBMUP.first ||
      beam_minus->pid() != writer->heprup.IDBMUP.second ||
      !SameRunValue(MomentumGeV(beam_plus).e(), writer->heprup.EBMUP.first) ||
      !SameRunValue(MomentumGeV(beam_minus).e(), writer->heprup.EBMUP.second)) {
    throw std::invalid_argument(
        "MLHEWriter: event beam identity differs from the run block");
  }
  const MLHECrossSection cross_section = EventCrossSection(ev);
  if (!SameRunValue(cross_section.value, writer->heprup.XSECUP[0]) ||
      !SameRunValue(cross_section.error, writer->heprup.XERRUP[0])) {
    throw std::invalid_argument(
        "MLHEWriter: event cross section differs from the run block");
  }
}

// Write one HepMC event as one Les Houches event
void MLHEWriter::WriteEvent(const HepMC3::GenEvent &ev) {
  if (writer == nullptr) {
    throw std::runtime_error("MLHEWriter::WriteEvent: writer is not open");
  }

  // Validate the event before allowing it to define the fixed run block
  LHEF::HEPEUP hepeup = BuildHEPEUP(ev, config);
  hepeup.junk = BuildLHEDiffractionComment(ev);
  InitializeRunBlock(ev);
  ValidateRunBlock(ev);
  writer->hepeup = hepeup;
  writer->hepeup.heprup = &writer->heprup;
  writer->writeEvent();
  output.flush();
  if (!output.good()) { throw std::runtime_error("MLHEWriter: failed to write event output"); }
}

// Close and validate the complete LHE stream including its final XML tag
void MLHEWriter::Close() {
  if (writer == nullptr) { return; }
  writer.reset();
  output.flush();
  output.close();
  if (output.fail()) { throw std::runtime_error("MLHEWriter: failed to close event output"); }
}

// Convert one HepMC3 ASCII file into one Les Houches Event file
MLHEConversionStats ConvertHepMC3ToLHE(const std::string &inputfile,
                                       const std::string &outputfile,
                                       bool progress) {
  std::error_code error;
  if (std::filesystem::equivalent(inputfile, outputfile, error)) {
    throw std::invalid_argument("ConvertHepMC3ToLHE: input and output refer to the same file");
  }
  const MLHEInputRun run = InspectHepMCWeights(inputfile);
  MHepMCReader input(inputfile);
  MLHEWriter writer(outputfile, run.weights);
  HepMC3::GenEvent ev(HepMC3::Units::GEV, HepMC3::Units::MM);

  MLHEConversionStats stats;
  while (input.Read(ev)) {

    auto xs = ev.cross_section() != nullptr
                  ? std::make_shared<HepMC3::GenCrossSection>(*ev.cross_section())
                  : std::make_shared<HepMC3::GenCrossSection>();
    xs->set_cross_section(run.xs.value, run.xs.error, xs->get_accepted_events(), xs->get_attempted_events());
    ev.set_cross_section(xs);
    writer.WriteEvent(ev);
    ++stats.events;

    if (progress && stats.events % 10000 == 0) {
      printf("%d events processed \n", stats.events);
    }
  }

  writer.Close();
  if (stats.events > 0) {
    stats.input_size_mb = gra::aux::GetFileSize(inputfile) / 1.0e6;
    stats.output_size_mb = gra::aux::GetFileSize(outputfile) / 1.0e6;
  }
  return stats;
}

} // namespace gra
