// Shower and hadronize one Les Houches Event file with Pythia8
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Pythia8
#include "Pythia8/Pythia.h"
#include "Pythia8/MiniStringFragmentation.h"
#include "Pythia8Plugins/HepMC3.h"

// HepMC3
#include "HepMC3/Attribute.h"
#include "HepMC3/FourVector.h"
#include "HepMC3/GenCrossSection.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenPdfInfo.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/WriterAscii.h"

#ifndef PYTHIA_XML_DIR
#define PYTHIA_XML_DIR ""
#endif

using namespace Pythia8;

namespace {

constexpr int PDG_POMERON      = 990;
constexpr int PDG_STABLE       = 1;
constexpr int PDG_INTERMEDIATE = 81;

// Compute the default number of shower attempts for one hard event
int DefaultReshowerAttempts() { return 500; }

// Generate deterministic event-local Pythia seeds from one user seed
struct PythiaSeedStream {
  std::mt19937 generator;

  // Initialize the helper seed stream with Pythia's default when zero is requested
  explicit PythiaSeedStream(int seed) : generator(static_cast<std::mt19937::result_type>(seed > 0 ? seed : 19780503)) {}

  // Compute the next valid Pythia seed for an independent shower attempt
  int NextSeed() { return 1 + static_cast<int>(generator() % 900000000U); }
};

// Compute a string without leading or trailing whitespace
std::string TrimCopy(const std::string &input) {
  const auto begin = std::find_if_not(input.begin(), input.end(), [](unsigned char c) { return std::isspace(c); });
  const auto end =
      std::find_if_not(input.rbegin(), input.rend(), [](unsigned char c) { return std::isspace(c); }).base();
  if (begin >= end) { return ""; }
  return std::string(begin, end);
}

struct DiffSideInfo {
  bool   active             = false;
  bool   excited            = false;
  bool   fragmented         = false;
  bool   string_skeleton    = false;
  bool   nuclear            = false;
  int    nuclear_a          = 0;
  int    nuclear_z          = 0;
  int    nuclear_parent_pdg = 0;
  double neutron_threshold  = 0.0;
  double xi                 = 0.0;
  double beta               = 0.0;
  double t                  = 0.0;
  int    system_pdg         = 0;
  double system_px          = 0.0;
  double system_py          = 0.0;
  double system_pz          = 0.0;
  double system_e           = 0.0;
  double system_m           = 0.0;
  int    lead_pdg           = 0;
  double lead_px            = 0.0;
  double lead_py            = 0.0;
  double lead_pz            = 0.0;
  double lead_e             = 0.0;
  double lead_m             = 0.0;
  int    rem_pdg            = 0;
  double rem_px             = 0.0;
  double rem_py             = 0.0;
  double rem_pz             = 0.0;
  double rem_e              = 0.0;
  double rem_m              = 0.0;
};

struct DiffEventInfo {
  bool         active           = false;
  bool         elastic_cep      = false;
  long         accepted_events  = 0;
  long         attempted_events = 0;
  int          side_mask        = 0;
  int          id1              = 0;
  int          id2              = 0;
  double       xhard1           = 0.0;
  double       xhard2           = 0.0;
  double       beam1_e          = 0.0;
  double       beam2_e          = 0.0;
  DiffSideInfo side1;
  DiffSideInfo side2;
};

struct FourMomentumSum {
  double px = 0.0;
  double py = 0.0;
  double pz = 0.0;
  double e  = 0.0;
};

struct ClosureResult {
  FourMomentumSum initial;
  FourMomentumSum final;
  FourMomentumSum delta;
  double          scale     = 1.0;
  double          max_abs   = 0.0;
  double          tolerance = 0.0;
  bool            pass      = true;
};

struct RecoilResult {
  bool   attempted      = false;
  bool   applied        = false;
  int    particle_count = 0;
  double lambda         = 1.0;
  double current_mass   = 0.0;
  double target_mass    = 0.0;
  double mapped_mass    = 0.0;
  double target_px      = 0.0;
  double target_py      = 0.0;
  double target_pz      = 0.0;
  double target_e       = 0.0;
  double mapped_e       = 0.0;
};

struct LHEInitInfo {
  int                      beam1_id   = 2212;
  int                      beam2_id   = 2212;
  double                   beam1_e    = 6500.0;
  double                   beam2_e    = 6500.0;
  double                   beam1_mass = -1.0;
  double                   beam2_mass = -1.0;
  int                      strategy   = 3;
  std::vector<std::string> process_lines;
};

struct LHEParticle {
  int    id      = 0;
  int    status  = 0;
  int    mother1 = 0;
  int    mother2 = 0;
  int    color1  = 0;
  int    color2  = 0;
  double px      = 0.0;
  double py      = 0.0;
  double pz      = 0.0;
  double e       = 0.0;
  double m       = 0.0;
  double tau     = 0.0;
  double spin    = 9.0;
};

struct LHEEventRecord {
  bool                     valid   = false;
  int                      nup     = 0;
  int                      idprup  = 0;
  double                   weight  = 1.0;
  double                   muF     = 0.0;
  double                   muR     = 0.0;
  double                   scalup  = 0.0;
  double                   alphaEM = 0.0;
  double                   alphaS  = 0.0;
  std::vector<LHEParticle> particles;
};

struct LHEEventBlock {
  DiffEventInfo  diff;
  LHEEventRecord record;
  LHEEventRecord hard;
};

struct LHEInput {
  LHEInitInfo                init;
  std::vector<LHEEventBlock> events;
};

// Compute the invariant mass squared of a four-momentum sum
double Mass2(const FourMomentumSum &p) { return p.e * p.e - p.px * p.px - p.py * p.py - p.pz * p.pz; }

// Compute the positive invariant mass of a physical four-momentum sum
double PositiveMass(const FourMomentumSum &p) { return std::sqrt(std::max(0.0, Mass2(p))); }

// Compute the difference of two four-momentum sums
FourMomentumSum SubtractMomentum(const FourMomentumSum &a, const FourMomentumSum &b) {
  return {a.px - b.px, a.py - b.py, a.pz - b.pz, a.e - b.e};
}

// Compute the sum of two four-momentum sums
FourMomentumSum AddMomentum(const FourMomentumSum &a, const FourMomentumSum &b) {
  return {a.px + b.px, a.py + b.py, a.pz + b.pz, a.e + b.e};
}

// Compute a four-momentum sum from one HepMC four-vector
FourMomentumSum FromHepMC(const HepMC3::FourVector &p) { return {p.px(), p.py(), p.pz(), p.e()}; }

// Compute the HepMC four-vector representation of one four-momentum sum
HepMC3::FourVector ToHepMC(const FourMomentumSum &p) { return HepMC3::FourVector(p.px, p.py, p.pz, p.e); }

// Compute true when one HepMC four-vector matches an expected four-vector
bool MomentumMatches(const HepMC3::FourVector &p, const FourMomentumSum &target) {
  const double scale = std::max(1.0, std::fabs(target.e));
  const double tol   = 1.0e-8 + 1.0e-10 * scale;
  return std::fabs(p.px() - target.px) < tol && std::fabs(p.py() - target.py) < tol &&
         std::fabs(p.pz() - target.pz) < tol && std::fabs(p.e() - target.e) < tol;
}

// Compute a beta vector protected against numerical boundary excursions
std::array<double, 3> BetaVector(const FourMomentumSum &p) {
  if (std::fpclassify(p.e) == FP_ZERO) { return {0.0, 0.0, 0.0}; }
  std::array<double, 3> beta = {p.px / p.e, p.py / p.e, p.pz / p.e};
  double                b2   = beta[0] * beta[0] + beta[1] * beta[1] + beta[2] * beta[2];
  if (b2 >= 1.0) {
    constexpr double eps   = 1.0e-14;
    const double     scale = std::sqrt((1.0 - eps) / b2);
    beta[0] *= scale;
    beta[1] *= scale;
    beta[2] *= scale;
  }
  return beta;
}

// Compute the opposite beta vector
std::array<double, 3> NegativeBeta(const std::array<double, 3> &beta) { return {-beta[0], -beta[1], -beta[2]}; }

// Compute a four-vector boosted by one beta vector
FourMomentumSum BoostMomentum(const FourMomentumSum &p, const std::array<double, 3> &beta) {
  const double b2 = beta[0] * beta[0] + beta[1] * beta[1] + beta[2] * beta[2];
  if (std::fpclassify(b2) == FP_ZERO) { return p; }
  const double dot    = beta[0] * p.px + beta[1] * p.py + beta[2] * p.pz;
  const double gamma  = 1.0 / std::sqrt(1.0 - b2);
  const double gamma2 = (gamma - 1.0) / b2;
  return {p.px + gamma2 * dot * beta[0] + gamma * beta[0] * p.e, p.py + gamma2 * dot * beta[1] + gamma * beta[1] * p.e,
          p.pz + gamma2 * dot * beta[2] + gamma * beta[2] * p.e, gamma * (p.e + dot)};
}

// Compute whether one forward side requests explicit neutron fragmentation
bool NuclearBreakup(const DiffSideInfo &side) {
  return side.active && side.excited && side.nuclear && side.nuclear_a > 1 && side.nuclear_z > 0 &&
         side.nuclear_z < side.nuclear_a && side.nuclear_parent_pdg != 0 && side.system_m > 0.0 &&
         side.neutron_threshold > 0.0;
}

// Compute the minimally excited one-neutron parent mass
double NuclearExcitedMass(const DiffSideInfo &side) {
  constexpr double kinetic_release = 1.0e-3;
  return NuclearBreakup(side) ? side.system_m + side.neutron_threshold + kinetic_release : side.system_m;
}

// Compute one physical forward-side momentum from diffraction metadata
FourMomentumSum ForwardSideMomentum(const DiffSideInfo &side) {
  if (side.excited) { return {side.system_px, side.system_py, side.system_pz, side.system_e}; }
  return {side.lead_px, side.lead_py, side.lead_pz, side.lead_e};
}

// Store one remapped forward-side momentum without changing its ground mass
void SetForwardSideMomentum(DiffSideInfo &side, const FourMomentumSum &momentum) {
  if (side.excited) {
    side.system_px = momentum.px;
    side.system_py = momentum.py;
    side.system_pz = momentum.pz;
    side.system_e  = momentum.e;
    return;
  }
  side.lead_px = momentum.px;
  side.lead_py = momentum.py;
  side.lead_pz = momentum.pz;
  side.lead_e  = momentum.e;
}

// Raise tagged nuclear masses while preserving the complete forward-pair momentum
bool PrepareNuclearPair(DiffEventInfo &diff) {
  if (!NuclearBreakup(diff.side1) && !NuclearBreakup(diff.side2)) { return true; }
  if (!diff.side1.active || !diff.side2.active) { return false; }

  const FourMomentumSum original1 = ForwardSideMomentum(diff.side1);
  const FourMomentumSum original2 = ForwardSideMomentum(diff.side2);
  const FourMomentumSum total     = AddMomentum(original1, original2);
  const double          s         = Mass2(total);
  if (!(s > 0.0)) { return false; }
  const double root_s     = std::sqrt(s);
  const double mass1      = NuclearBreakup(diff.side1) ? NuclearExcitedMass(diff.side1) : PositiveMass(original1);
  const double mass2      = NuclearBreakup(diff.side2) ? NuclearExcitedMass(diff.side2) : PositiveMass(original2);
  const double sum        = mass1 + mass2;
  const double difference = mass1 - mass2;
  const double kallen     = (s - sum * sum) * (s - difference * difference);
  if (!(mass1 > 0.0) || !(mass2 >= 0.0) || !(kallen >= 0.0)) { return false; }

  const std::array<double, 3> beta      = BetaVector(total);
  const FourMomentumSum       rest1     = BoostMomentum(original1, NegativeBeta(beta));
  const double                norm      = std::sqrt(rest1.px * rest1.px + rest1.py * rest1.py + rest1.pz * rest1.pz);
  const std::array<double, 3> direction = norm > 0.0
                                              ? std::array<double, 3>{rest1.px / norm, rest1.py / norm, rest1.pz / norm}
                                              : std::array<double, 3>{0.0, 0.0, 1.0};
  const double                momentum  = std::sqrt(kallen) / (2.0 * root_s);
  const FourMomentumSum       mapped1{momentum * direction[0], momentum * direction[1], momentum * direction[2],
                                std::sqrt(mass1 * mass1 + momentum * momentum)};
  const FourMomentumSum mapped2{-mapped1.px, -mapped1.py, -mapped1.pz, std::sqrt(mass2 * mass2 + momentum * momentum)};
  const FourMomentumSum lab1 = BoostMomentum(mapped1, beta);
  const FourMomentumSum lab2 = BoostMomentum(mapped2, beta);
  if (!std::isfinite(lab1.e) || !std::isfinite(lab2.e)) { return false; }
  SetForwardSideMomentum(diff.side1, lab1);
  SetForwardSideMomentum(diff.side2, lab2);
  return true;
}

// Compute the default beam-particle mass used in closure checks
double BeamMassGeV(int pdg) {
  switch (std::abs(pdg)) {
    case 2212:
      return 0.93827208943;
    case 2112:
      return 0.9395654219;
    case 11:
      return 0.00051099895;
    case 13:
      return 0.1056583755;
    case 15:
      return 1.77693;
    default:
      return 0.0;
  }
}

// Compute one event-local beam mass from exact forward metadata when available
double BeamMassGeV(int pdg, const DiffSideInfo &side) {
  if (side.active && side.nuclear && side.system_m > 0.0 && side.nuclear_parent_pdg == pdg) { return side.system_m; }
  if (side.active && !side.excited && side.lead_m > 0.0 && side.lead_pdg == pdg) { return side.lead_m; }
  return BeamMassGeV(pdg);
}

// Compute the signed longitudinal momentum of one original beam
double BeamPz(int pdg, double energy, bool first_side, double mass = -1.0) {
  if (!(mass >= 0.0)) { mass = BeamMassGeV(pdg); }
  const double pz2 = std::max(0.0, energy * energy - mass * mass);
  const double pz  = std::sqrt(pz2);
  return first_side ? pz : -pz;
}

// Compute the original GRANIITTI pp initial-state four-momentum
FourMomentumSum OriginalInitialMomentum(const LHEInitInfo &init, const DiffEventInfo &diff) {
  const double    beam1_e = diff.beam1_e > 0.0 ? diff.beam1_e : init.beam1_e;
  const double    beam2_e = diff.beam2_e > 0.0 ? diff.beam2_e : init.beam2_e;
  FourMomentumSum sum;
  sum.pz = BeamPz(init.beam1_id, beam1_e, true,
                  init.beam1_mass >= 0.0 ? init.beam1_mass : BeamMassGeV(init.beam1_id, diff.side1)) +
           BeamPz(init.beam2_id, beam2_e, false,
                  init.beam2_mass >= 0.0 ? init.beam2_mass : BeamMassGeV(init.beam2_id, diff.side2));
  sum.e = beam1_e + beam2_e;
  return sum;
}

// Compute one original GRANIITTI pp beam four-momentum
FourMomentumSum OriginalBeamMomentum(const LHEInitInfo &init, const DiffEventInfo &diff, bool first_side) {
  const double        beam_e  = first_side ? (diff.beam1_e > 0.0 ? diff.beam1_e : init.beam1_e)
                                           : (diff.beam2_e > 0.0 ? diff.beam2_e : init.beam2_e);
  const int           beam_id = first_side ? init.beam1_id : init.beam2_id;
  const DiffSideInfo &side    = first_side ? diff.side1 : diff.side2;
  const double        mass    = first_side ? init.beam1_mass : init.beam2_mass;
  return {0.0, 0.0, BeamPz(beam_id, beam_e, first_side, mass >= 0.0 ? mass : BeamMassGeV(beam_id, side)), beam_e};
}

// Compute the exact GRANIITTI leading-proton four-momentum sum
FourMomentumSum LeadingProtonMomentum(const DiffEventInfo &diff) {
  FourMomentumSum sum;
  const auto      add_side = [&sum](const DiffSideInfo &side) {
    if (!side.active || side.lead_pdg == 0) { return; }
    sum.px += side.lead_px;
    sum.py += side.lead_py;
    sum.pz += side.lead_pz;
    sum.e += side.lead_e;
  };
  add_side(diff.side1);
  add_side(diff.side2);
  return sum;
}

// Compute the physical four-momentum left for the Pythia-produced subsystem
FourMomentumSum TargetPythiaFinalMomentum(const LHEInitInfo &init, const DiffEventInfo &diff) {
  return SubtractMomentum(OriginalInitialMomentum(init, diff), LeadingProtonMomentum(diff));
}

// Compute mutable stable HepMC3 final-state particles for kinematic projection
std::vector<HepMC3::GenParticlePtr> StableFinalParticles(HepMC3::GenEvent &event) {
  std::vector<HepMC3::GenParticlePtr> particles;
  for (const auto &particle : event.particles()) {
    if (!particle || particle->status() != PDG_STABLE) { continue; }
    particles.push_back(particle);
  }
  return particles;
}

// Compute true for colored diquark remnant ids
bool IsDiquark(int pdg) {
  const int apdg = std::abs(pdg);
  const int spin = apdg % 10;
  return apdg >= 1000 && apdg < 6000 && (spin == 1 || spin == 3);
}

// Compute true for colored parton and diquark ids
bool IsColoredPDG(int pdg) { return pdg == 21 || (std::abs(pdg) >= 1 && std::abs(pdg) <= 5) || IsDiquark(pdg); }

// Protect colorless hard decay branches while allowing connected QCD strings to take recoil
bool IsColorlessHardRow(const LHEParticle &particle) {
  return (particle.status == 1 || particle.status == 2) && !IsColoredPDG(particle.id) && particle.color1 == 0 && particle.color2 == 0;
}

// Compute true when a particle pointer has already been selected
bool ContainsParticle(const std::vector<HepMC3::GenParticlePtr> &particles, const HepMC3::GenParticlePtr &candidate) {
  return std::find(particles.begin(), particles.end(), candidate) != particles.end();
}

// Add one particle pointer if it is not already present
void AddUniqueParticle(std::vector<HepMC3::GenParticlePtr> &particles, const HepMC3::GenParticlePtr &candidate) {
  if (!candidate || ContainsParticle(particles, candidate)) { return; }
  particles.push_back(candidate);
}

// Compute true when one HepMC particle matches one colorless LHE row
bool MatchesHardFinalRoot(const HepMC3::GenParticlePtr &particle, const LHEParticle &row) {
  if (!particle || particle->pid() != row.id) { return false; }
  const FourMomentumSum target{row.px, row.py, row.pz, row.e};
  return MomentumMatches(particle->momentum(), target);
}

// Compute the HepMC hard branch root matched to one colorless LHE row
HepMC3::GenParticlePtr MatchHardFinalRoot(HepMC3::GenEvent &event, const LHEParticle &row,
                                          const std::vector<HepMC3::GenParticlePtr> &used_roots) {
  HepMC3::GenParticlePtr fallback;
  for (const auto &particle : event.particles()) {
    if (!MatchesHardFinalRoot(particle, row) || ContainsParticle(used_roots, particle)) { continue; }
    if (particle->status() != PDG_STABLE) { return particle; }
    if (!fallback) { fallback = particle; }
  }
  if (fallback) { return fallback; }
  // ISR and primordial transverse momentum can move the hard roots away from the LHE momenta
  for (const auto &particle : event.particles()) {
    if (!particle || particle->pid() != row.id || ContainsParticle(used_roots, particle)) { continue; }
    const auto vertex = particle->production_vertex();
    if (!vertex) { continue; }
    for (const auto &parent : vertex->particles_in()) {
      if (parent->status() == 21 || parent->status() == 22) { return particle; }
    }
  }
  return fallback;
}

// Recursively collect stable descendants of one hard-branch particle
void CollectStableDescendants(const HepMC3::GenParticlePtr &particle, std::vector<HepMC3::GenParticlePtr> &stable,
                              std::vector<HepMC3::GenParticlePtr> &visited) {
  if (!particle || ContainsParticle(visited, particle)) { return; }
  AddUniqueParticle(visited, particle);

  if (particle->status() == PDG_STABLE) {
    AddUniqueParticle(stable, particle);
    return;
  }

  const HepMC3::GenVertexPtr end = particle->end_vertex();
  if (!end) { return; }
  for (const auto &child : end->particles_out()) { CollectStableDescendants(child, stable, visited); }
}

// Collect stable descendants of the colorless LHE hard branches
bool ProtectedHardBranchParticles(HepMC3::GenEvent &event, const LHEEventRecord &record,
                                  std::vector<HepMC3::GenParticlePtr> &protected_particles, int &hard_rows,
                                  int &matched_roots) {
  protected_particles.clear();
  hard_rows     = 0;
  matched_roots = 0;
  std::vector<HepMC3::GenParticlePtr> roots;
  std::vector<HepMC3::GenParticlePtr> visited;

  std::vector<bool> protected_rows;
  for (const LHEParticle &row : record.particles) {
    // Protect each colorless branch once, starting at its earliest LHE ancestor
    const auto inherited = [&](int mother) {
      return mother > 0 && static_cast<std::size_t>(mother) <= protected_rows.size() && protected_rows[mother - 1];
    };
    const bool descendant = inherited(row.mother1) || inherited(row.mother2);
    const bool colorless = IsColorlessHardRow(row);
    protected_rows.push_back(descendant || colorless);
    if (!colorless || descendant) { continue; }
    ++hard_rows;
    const HepMC3::GenParticlePtr root = MatchHardFinalRoot(event, row, roots);
    if (!root) { return false; }
    AddUniqueParticle(roots, root);
    ++matched_roots;
    CollectStableDescendants(root, protected_particles, visited);
  }

  return matched_roots == hard_rows && protected_particles.size() >= static_cast<std::size_t>(hard_rows);
}

// Compute stable final-state particles split by hard-branch ancestry
std::vector<HepMC3::GenParticlePtr> StableFinalParticlesByProtection(
    HepMC3::GenEvent &event, const std::vector<HepMC3::GenParticlePtr> &protected_particles, bool select_protected) {
  std::vector<HepMC3::GenParticlePtr> particles;
  for (const auto &particle : event.particles()) {
    if (!particle || particle->status() != PDG_STABLE) { continue; }
    const bool is_protected = ContainsParticle(protected_particles, particle);
    if (is_protected == select_protected) { particles.push_back(particle); }
  }
  return particles;
}

// Compute the summed four-momentum of one particle list
FourMomentumSum ParticleListMomentum(const std::vector<HepMC3::GenParticlePtr> &particles) {
  FourMomentumSum sum;
  for (const auto &particle : particles) {
    if (!particle) { continue; }
    const HepMC3::FourVector &p = particle->momentum();
    sum.px += p.px();
    sum.py += p.py();
    sum.pz += p.pz();
    sum.e += p.e();
  }
  return sum;
}

// Compute the summed four-momentum of stable HepMC3 final-state particles
FourMomentumSum StableFinalMomentum(const HepMC3::GenEvent &event) {
  FourMomentumSum sum;
  for (const auto &particle : event.particles()) {
    if (!particle || particle->status() != PDG_STABLE) { continue; }
    const HepMC3::FourVector &p = particle->momentum();
    sum.px += p.px();
    sum.py += p.py();
    sum.pz += p.pz();
    sum.e += p.e();
  }
  return sum;
}

// Compute the closure tolerance for one event-level four-momentum check
double ClosureTolerance(const ClosureResult &result) { return 1.0e-2 + 1.0e-8 * result.scale; }

// Build the total event four-momentum closure residual
ClosureResult MomentumClosure(const HepMC3::GenEvent &event, const LHEInitInfo &init, const DiffEventInfo &diff) {
  ClosureResult result;
  result.initial = OriginalInitialMomentum(init, diff);
  result.final   = StableFinalMomentum(event);
  result.delta   = {result.final.px - result.initial.px, result.final.py - result.initial.py,
                    result.final.pz - result.initial.pz, result.final.e - result.initial.e};
  result.scale   = std::max(1.0, std::fabs(result.initial.e));
  result.max_abs = std::max(
      {std::fabs(result.delta.px), std::fabs(result.delta.py), std::fabs(result.delta.pz), std::fabs(result.delta.e)});
  result.tolerance = ClosureTolerance(result);
  result.pass      = result.max_abs <= result.tolerance;
  for (const double component : {result.delta.px, result.delta.py, result.delta.pz, result.delta.e}) {
    if (!std::isfinite(component)) { result.pass = false; }
  }
  return result;
}

// Read one positive integer command-line argument
int ReadPositiveInt(const char *value, const std::string &name) {
  std::size_t used   = 0;
  const int   parsed = std::stoi(value, &used);
  if (used != std::string(value).size()) { throw std::invalid_argument("Invalid " + name); }
  if (parsed <= 0) { throw std::invalid_argument(name + " must be positive"); }
  return parsed;
}

// Read one non-negative integer command-line argument
int ReadNonNegativeInt(const char *value, const std::string &name) {
  std::size_t used   = 0;
  const int   parsed = std::stoi(value, &used);
  if (used != std::string(value).size() || parsed > 900000000) { throw std::invalid_argument("Invalid " + name); }
  if (parsed < 0) { throw std::invalid_argument(name + " must be non-negative"); }
  return parsed;
}

// Compute one positive optional reshower-attempt count
int ReadReshowerAttempts(int argc, char *argv[]) {
  if (argc < 7) { return DefaultReshowerAttempts(); }
  return ReadPositiveInt(argv[6], "max_reshower_attempts");
}

// Compute one line stripped of leading spaces and one optional LHE comment marker
std::string StripLHECommentPrefix(const std::string &line) {
  std::size_t pos = line.find_first_not_of(" \t");
  if (pos == std::string::npos) { return ""; }
  if (line[pos] == '#') {
    ++pos;
    while (pos < line.size() && (line[pos] == ' ' || line[pos] == '\t')) { ++pos; }
  }
  return line.substr(pos);
}

// Compute true when a line is not an XML tag or comment-only line
bool IsDataLine(const std::string &line) {
  const std::size_t pos = line.find_first_not_of(" \t");
  if (pos == std::string::npos) { return false; }
  return line[pos] != '<' && line[pos] != '#';
}

// Read one floating-point XML attribute without assuming attribute order
bool ReadXMLDoubleAttribute(const std::string &line, const std::string &name, double &value) {
  const std::string key = name + "=";
  std::size_t       pos = line.find(key);
  if (pos == std::string::npos) { return false; }
  pos += key.size();
  if (pos >= line.size() || (line[pos] != '\'' && line[pos] != '"')) { return false; }
  const char        quote = line[pos++];
  const std::size_t end   = line.find(quote, pos);
  if (end == std::string::npos) { return false; }
  try {
    const double parsed = std::stod(line.substr(pos, end - pos));
    if (!std::isfinite(parsed)) { return false; }
    value = parsed;
    return true;
  } catch (...) { return false; }
}

// Read LHE3 factorization and renormalization scales with SCALUP fallbacks
void ParseLHEGeneratorScales(const std::vector<std::string> &lines, LHEEventRecord &record) {
  record.muF = record.scalup;
  record.muR = record.scalup;
  for (const auto &line : lines) {
    if (line.find("<scales") != std::string::npos) {
      ReadXMLDoubleAttribute(line, "muf", record.muF);
      ReadXMLDoubleAttribute(line, "mur", record.muR);
    }
    if (line.find("<pdfinfo") != std::string::npos) { ReadXMLDoubleAttribute(line, "scale", record.muF); }
  }
}

// Parse the fixed-width-free LHE event payload into a Pythia LHAup record
LHEEventRecord ParseLHEEventRecord(const std::vector<std::string> &lines) {
  LHEEventRecord           record;
  std::vector<std::string> data_lines;
  for (const auto &line : lines) {
    if (IsDataLine(line)) { data_lines.push_back(line); }
  }
  if (data_lines.empty()) { return record; }

  std::istringstream header(data_lines.front());
  header >> record.nup >> record.idprup >> record.weight >> record.scalup >> record.alphaEM >> record.alphaS;
  if (!header || !std::isfinite(record.weight) || !std::isfinite(record.scalup) || record.nup <= 0 ||
      static_cast<std::size_t>(record.nup) + 1 > data_lines.size()) {
    return LHEEventRecord();
  }

  record.particles.reserve(static_cast<std::size_t>(record.nup));
  for (int i = 0; i < record.nup; ++i) {
    LHEParticle        particle;
    std::istringstream row(data_lines[static_cast<std::size_t>(i + 1)]);
    row >> particle.id >> particle.status >> particle.mother1 >> particle.mother2 >> particle.color1 >>
        particle.color2 >> particle.px >> particle.py >> particle.pz >> particle.e >> particle.m >> particle.tau >>
        particle.spin;
    if (!row || !std::isfinite(particle.px) || !std::isfinite(particle.py) || !std::isfinite(particle.pz) ||
        !std::isfinite(particle.e) || !std::isfinite(particle.m) || particle.e < 0.0 || particle.mother1 < 0 ||
        particle.mother2 < 0 || particle.mother1 > record.nup || particle.mother2 > record.nup || particle.color1 < 0 ||
        particle.color2 < 0) {
      return LHEEventRecord();
    }
    record.particles.push_back(particle);
  }
  ParseLHEGeneratorScales(lines, record);
  record.valid = true;
  return record;
}

// Parse key=value tokens from one metadata line
std::map<std::string, std::string> ParseKeyValues(const std::string &line) {
  std::map<std::string, std::string> values;
  std::istringstream                 input(line);
  std::string                        token;
  input >> token;
  while (input >> token) {
    const std::size_t eq = token.find('=');
    if (eq == std::string::npos) { continue; }
    values[token.substr(0, eq)] = token.substr(eq + 1);
  }
  return values;
}

// Read one integer metadata value with a fallback
int GetIntValue(const std::map<std::string, std::string> &values, const std::string &key, int fallback) {
  const auto found = values.find(key);
  if (found == values.end()) { return fallback; }
  std::size_t used = 0;
  const auto value = std::stoi(found->second, &used);
  if (used != found->second.size()) { throw std::invalid_argument("Invalid diffraction metadata " + key); }
  return value;
}

// Read one long integer metadata value with a fallback
long GetLongValue(const std::map<std::string, std::string> &values, const std::string &key, long fallback) {
  const auto found = values.find(key);
  if (found == values.end()) { return fallback; }
  std::size_t used = 0;
  const auto value = std::stol(found->second, &used);
  if (used != found->second.size()) { throw std::invalid_argument("Invalid diffraction metadata " + key); }
  return value;
}

// Read one floating-point metadata value with a fallback
double GetDoubleValue(const std::map<std::string, std::string> &values, const std::string &key, double fallback) {
  const auto found = values.find(key);
  if (found == values.end()) { return fallback; }
  std::size_t  used  = 0;
  const double value = std::stod(found->second, &used);
  if (used != found->second.size() || !std::isfinite(value)) {
    throw std::invalid_argument("Invalid diffraction metadata " + key);
  }
  return value;
}

// Fill one diffractive-side metadata object from parsed key-value tokens
void FillDiffSide(DiffSideInfo &side, const std::map<std::string, std::string> &values, const std::string &prefix) {
  side.active = GetIntValue(values, prefix + "_active", 0) != 0;
  if (!side.active) { return; }
  side.excited            = GetIntValue(values, prefix + "_excited", 0) != 0;
  side.fragmented         = GetIntValue(values, prefix + "_fragmented", 0) != 0;
  side.string_skeleton    = GetIntValue(values, prefix + "_string", 0) != 0;
  side.nuclear            = GetIntValue(values, prefix + "_nuclear", 0) != 0;
  side.nuclear_a          = GetIntValue(values, prefix + "_nuclear_a", 0);
  side.nuclear_z          = GetIntValue(values, prefix + "_nuclear_z", 0);
  side.nuclear_parent_pdg = GetIntValue(values, prefix + "_nuclear_parent_pdg", 0);
  side.neutron_threshold  = GetDoubleValue(values, prefix + "_neutron_threshold", 0.0);
  side.xi                 = GetDoubleValue(values, prefix + "_xi", 0.0);
  side.beta               = GetDoubleValue(values, prefix + "_beta", 0.0);
  side.t                  = GetDoubleValue(values, prefix + "_t", 0.0);
  side.system_pdg         = GetIntValue(values, prefix + "_system_pdg", 0);
  side.system_px          = GetDoubleValue(values, prefix + "_system_px", 0.0);
  side.system_py          = GetDoubleValue(values, prefix + "_system_py", 0.0);
  side.system_pz          = GetDoubleValue(values, prefix + "_system_pz", 0.0);
  side.system_e           = GetDoubleValue(values, prefix + "_system_e", 0.0);
  side.system_m           = GetDoubleValue(values, prefix + "_system_m", 0.0);
  side.lead_pdg           = GetIntValue(values, prefix + "_lead_pdg", 0);
  side.lead_px            = GetDoubleValue(values, prefix + "_lead_px", 0.0);
  side.lead_py            = GetDoubleValue(values, prefix + "_lead_py", 0.0);
  side.lead_pz            = GetDoubleValue(values, prefix + "_lead_pz", 0.0);
  side.lead_e             = GetDoubleValue(values, prefix + "_lead_e", 0.0);
  side.lead_m             = GetDoubleValue(values, prefix + "_lead_m", 0.0);
  side.rem_pdg            = GetIntValue(values, prefix + "_rem_pdg", 0);
  side.rem_px             = GetDoubleValue(values, prefix + "_rem_px", 0.0);
  side.rem_py             = GetDoubleValue(values, prefix + "_rem_py", 0.0);
  side.rem_pz             = GetDoubleValue(values, prefix + "_rem_pz", 0.0);
  side.rem_e              = GetDoubleValue(values, prefix + "_rem_e", 0.0);
  side.rem_m              = GetDoubleValue(values, prefix + "_rem_m", 0.0);
}

// Parse one GRANIITTI diffractive metadata line when present
DiffEventInfo ParseDiffComment(const std::string &line) {
  DiffEventInfo     info;
  const std::string stripped = StripLHECommentPrefix(line);
  if (stripped.rfind("graniitti_diff ", 0) != 0) { return info; }

  const std::map<std::string, std::string> values = ParseKeyValues(stripped);
  info.active                                     = true;
  info.elastic_cep                                = GetIntValue(values, "elastic_cep", 0) != 0;
  info.accepted_events                            = GetLongValue(values, "accepted_events", 0);
  info.attempted_events                           = GetLongValue(values, "attempted_events", 0);
  info.side_mask                                  = GetIntValue(values, "side_mask", 0);
  info.id1                                        = GetIntValue(values, "id1", 0);
  info.id2                                        = GetIntValue(values, "id2", 0);
  info.xhard1                                     = GetDoubleValue(values, "xhard1", 0.0);
  info.xhard2                                     = GetDoubleValue(values, "xhard2", 0.0);
  info.beam1_e                                    = GetDoubleValue(values, "beam1_e", 0.0);
  info.beam2_e                                    = GetDoubleValue(values, "beam2_e", 0.0);
  FillDiffSide(info.side1, values, "side1");
  FillDiffSide(info.side2, values, "side2");
  return info;
}

struct LHEProcessInfo {
  int    id   = 0;
  double xsec = 1.0;
  double xerr = 0.0;
  double xmax = 1.0;
};

// Parse one LHE init process with finite cross section and uncertainty
LHEProcessInfo ParseLHEProcessInfo(const std::string &line) {
  LHEProcessInfo     process;
  std::istringstream input(line);
  input >> process.xsec >> process.xerr >> process.xmax >> process.id;
  if (!input || !std::isfinite(process.xsec) || !std::isfinite(process.xerr) || !std::isfinite(process.xmax) ||
      process.xerr < 0.0 || process.xmax < 0.0) {
    throw std::invalid_argument("Invalid LHE process cross section");
  }
  return process;
}

// Parse the LHE init block and event blocks needed by the Pythia bridge
LHEInput ReadLHEInput(const std::string &lhe_file, int max_events) {
  std::ifstream input(lhe_file);
  if (!input.good()) { throw std::runtime_error("failed to open LHE file " + lhe_file); }

  LHEInput                 data;
  bool                     in_init          = false;
  bool                     init_header_read = false;
  int                      process_count    = 0;
  bool                     in_event         = false;
  LHEEventBlock            current;
  std::vector<std::string> event_lines;
  std::vector<std::string> hard_lines;

  std::string line;
  while (std::getline(input, line)) {
    const std::string trimmed = TrimCopy(line);
    if (trimmed == "<init>" || trimmed.rfind("<init ", 0) == 0) {
      in_init = true;
      continue;
    }
    if (line.find("</init>") != std::string::npos) {
      in_init = false;
      continue;
    }
    if (in_init && IsDataLine(line)) {
      if (!init_header_read) {
        int                pdfg1 = 0, pdfg2 = 0, pdfs1 = 0, pdfs2 = 0, idwtup = 0;
        std::istringstream init_line(line);
        init_line >> data.init.beam1_id >> data.init.beam2_id >> data.init.beam1_e >> data.init.beam2_e >> pdfg1 >>
            pdfg2 >> pdfs1 >> pdfs2 >> idwtup >> process_count;
        if (!init_line || process_count < 1 || !std::isfinite(data.init.beam1_e) || !std::isfinite(data.init.beam2_e) ||
            !(data.init.beam1_e > 0.0) || !(data.init.beam2_e > 0.0)) {
          throw std::invalid_argument("Invalid LHE init block");
        }
        data.init.strategy = idwtup;
        init_header_read   = true;
      } else {
        const auto process = ParseLHEProcessInfo(line);
        for (const auto &previous : data.init.process_lines) {
          if (ParseLHEProcessInfo(previous).id == process.id) {
            throw std::invalid_argument("Duplicate LHE process id");
          }
        }
        data.init.process_lines.push_back(line);
      }
      continue;
    }

    if (trimmed == "<event>" || trimmed.rfind("<event ", 0) == 0) {
      if (in_event) { throw std::invalid_argument("Nested LHE event"); }
      in_event = true;
      event_lines.clear();
      hard_lines.clear();
    }
    if (in_event) {
      event_lines.push_back(line);
      const DiffEventInfo diff = ParseDiffComment(line);
      if (diff.active) { current.diff = diff; }
      const std::string comment = StripLHECommentPrefix(line);
      if (comment.rfind("graniitti_hard ", 0) == 0) { hard_lines.push_back(comment.substr(15)); }
    }
    if (in_event && line.find("</event>") != std::string::npos) {
      in_event       = false;
      current.record = ParseLHEEventRecord(event_lines);
      if (!current.record.valid) {
        throw std::invalid_argument("Invalid LHE event " + std::to_string(data.events.size() + 1));
      }
      const bool known_process = std::any_of(
          data.init.process_lines.begin(), data.init.process_lines.end(),
          [&current](const std::string &process) { return ParseLHEProcessInfo(process).id == current.record.idprup; });
      if (!known_process) { throw std::invalid_argument("LHE event has an unknown process id"); }
      if (!hard_lines.empty()) {
        std::ostringstream header;
        header << std::setprecision(17) << hard_lines.size() << ' ' << current.record.idprup << ' '
               << current.record.weight << ' ' << current.record.scalup << ' ' << current.record.alphaEM << ' '
               << current.record.alphaS;
        hard_lines.insert(hard_lines.begin(), header.str());
        current.hard     = ParseLHEEventRecord(hard_lines);
        current.hard.muF = current.record.muF;
        current.hard.muR = current.record.muR;
        if (!current.hard.valid) { throw std::invalid_argument("Invalid GRANIITTI hard-process rows"); }
      }
      if (data.events.empty() && current.record.particles.size() >= 2) {
        const auto &a = current.record.particles[0];
        const auto &b = current.record.particles[1];
        if (a.status == -1 && b.status == -1 && a.id == data.init.beam1_id && b.id == data.init.beam2_id) {
          data.init.beam1_mass = std::abs(a.m);
          data.init.beam2_mass = std::abs(b.m);
        }
      }
      data.events.push_back(std::exchange(current, LHEEventBlock{}));
      if (data.events.size() >= static_cast<std::size_t>(max_events)) { break; }
    }
  }
  if (in_event || in_init || !init_header_read ||
      data.init.process_lines.size() != static_cast<std::size_t>(process_count) || data.events.empty()) {
    throw std::invalid_argument("Incomplete LHE input " + lhe_file);
  }
  return data;
}

// Compute true when at least one LHE event carries GRANIITTI diffractive metadata
bool HasDiffractiveEvents(const LHEInput &input) {
  for (const auto &event : input.events) {
    if (event.diff.active) { return true; }
  }
  return false;
}

// Compute whether any requested event carries an Xn nuclear forward system
bool HasNuclearBreakup(const LHEInput &input, const int max_events) {
  int checked = 0;
  for (const auto &event : input.events) {
    if (!event.diff.active) { continue; }
    if (NuclearBreakup(event.diff.side1) || NuclearBreakup(event.diff.side2)) { return true; }
    if (++checked >= max_events) { break; }
  }
  return false;
}

// Require cumulative source counters for every GRANIITTI diffractive event
void ValidateDiffractiveCounters(const LHEInput &input) {
  long previous_accepted  = 0;
  long previous_attempted = 0;
  for (const auto &event : input.events) {
    if (!event.diff.active) { continue; }
    if (event.diff.accepted_events <= previous_accepted || event.diff.attempted_events < previous_attempted ||
        event.diff.attempted_events < event.diff.accepted_events) {
      throw std::runtime_error("GRANIITTI diffractive LHE event has missing or inconsistent cumulative event counters");
    }
    previous_accepted  = event.diff.accepted_events;
    previous_attempted = event.diff.attempted_events;
  }
}

// Print command-line usage
void PrintUsage(const char *argv0) {
  std::cerr << "Usage: " << argv0
            << " input.lhe output.hepmc3 nevents|all [seed] [extra.cmnd] [max_reshower_attempts] [auto|fragment|shower]\n";
}

// Compute the event-local Pythia beam id for one GRANIITTI side
int PythiaBeamId(const LHEInitInfo &init, const DiffSideInfo &side, bool first_side) {
  return side.active ? PDG_POMERON : (first_side ? init.beam1_id : init.beam2_id);
}

// Compute the original beam energy recorded for one GRANIITTI event side
double EventBeamEnergy(const LHEInitInfo &init, const DiffEventInfo &diff, bool first_side) {
  const double event_energy = first_side ? diff.beam1_e : diff.beam2_e;
  return event_energy > 0.0 ? event_energy : (first_side ? init.beam1_e : init.beam2_e);
}

// Compute the event-local Pythia beam energy for one GRANIITTI side
double PythiaBeamEnergy(const LHEInitInfo &init, const DiffEventInfo &diff, const DiffSideInfo &side, bool first_side) {
  const double beam_energy = EventBeamEnergy(init, diff, first_side);
  if (!side.active) { return beam_energy; }

  return side.xi * beam_energy;
}

// Compute the DPDF spectator energy left on one active diffractive side
double DPDFRemnantEnergy(const LHEInitInfo &init, const DiffEventInfo &diff, const DiffSideInfo &side,
                         bool first_side) {
  if (!side.active) { return 1.0e30; }
  const double beam_energy = side.xi * EventBeamEnergy(init, diff, first_side);
  return std::max(0.0, (1.0 - side.beta) * beam_energy);
}

// Compute the smallest active-side DPDF spectator energy in one event
double MinimumDPDFRemnantEnergy(const LHEInitInfo &init, const LHEEventBlock &event) {
  return std::min(DPDFRemnantEnergy(init, event.diff, event.diff.side1, true),
                  DPDFRemnantEnergy(init, event.diff, event.diff.side2, false));
}

// Compute the hard fraction relative to the event-local Pomeron or proton beam
double PythiaHardFraction(const LHEInitInfo &init, const LHEEventBlock &event, bool first_side) {
  const DiffSideInfo &side = first_side ? event.diff.side1 : event.diff.side2;
  if (event.record.particles.size() < 2) { return 0.0; }
  const double beam_energy     = PythiaBeamEnergy(init, event.diff, side, first_side);
  const double incoming_energy = event.record.particles[first_side ? 0 : 1].e;
  const double fraction        = beam_energy > 0.0 ? incoming_energy / beam_energy : 0.0;
  return (fraction > 0.0 && fraction < 1.0) ? fraction : 0.0;
}

// Compute true when all Pythia event-record momenta are finite
bool PythiaEventMomentaFinite(const Event &event) {
  for (int i = 0; i < event.size(); ++i) {
    if (!std::isfinite(event[i].px()) || !std::isfinite(event[i].py()) || !std::isfinite(event[i].pz()) ||
        !std::isfinite(event[i].e())) {
      return false;
    }
  }
  return true;
}

// Compute true for a supported incoming baryon beam
bool IsBaryonBeamPDG(int pdg) { return std::abs(pdg) == 2212 || std::abs(pdg) == 2112; }

// Compute true for a valid 10LZZZAAAI nuclear identity
bool IsNuclearPDG(int pdg) {
  const int code = std::abs(pdg);
  if (code < 1000000000) { return false; }
  const int a = (code / 10) % 1000;
  const int z = (code / 10000) % 1000;
  return a > 0 && z >= 0 && z <= a;
}

// Compute true for one colorless beam copied as a complete LHE record
bool IsFullBeamPDG(int pdg) {
  const int code = std::abs(pdg);
  return IsBaryonBeamPDG(pdg) || code == 11 || code == 13 || code == 15 || IsNuclearPDG(pdg);
}

// Compute true when the LHE event keeps the original baryon beams explicitly
bool HasFullBeamRows(const LHEEventRecord &record) {
  return record.particles.size() >= 2 && IsFullBeamPDG(record.particles[0].id) &&
         IsFullBeamPDG(record.particles[1].id) && record.particles[0].status == -1 && record.particles[1].status == -1;
}

// Recognize complete physical beam records independently of the production process
bool HasFullBeamRecords(const LHEInput &input) {
  return !input.events.empty() && std::all_of(input.events.begin(), input.events.end(),
                                              [](const LHEEventBlock &event) { return HasFullBeamRows(event.record); });
}

// Project the hard system onto collinear massless incoming partons for Pythia ISR
void PrepareHardShower(LHEEventBlock &block) {
  if (!block.diff.active) { throw std::invalid_argument("Hard shower mode requires beam-emission metadata"); }
  if (HasFullBeamRows(block.record)) {
    if (!block.hard.valid) {
      throw std::invalid_argument("Shower mode requires graniitti_hard rows in complete beam records");
    }
    block.record = std::move(block.hard);
  }
  auto &rows = block.record.particles;
  if (rows.size() < 3 || rows[0].status != -1 || rows[1].status != -1 || IsFullBeamPDG(rows[0].id) ||
      IsFullBeamPDG(rows[1].id)) {
    throw std::invalid_argument("Shower mode requires two incoming hard partons");
  }
  Vec4 total;
  for (const auto &row : rows) {
    if (row.status == 1) { total += Vec4(row.px, row.py, row.pz, row.e); }
  }
  const double mass = total.mCalc();
  if (!(mass > 0.0) || !(total.e() > std::abs(total.pz()))) {
    throw std::invalid_argument("Hard final state must have a timelike total momentum");
  }
  const double rapidity = 0.5 * std::log((total.e() + total.pz()) / (total.e() - total.pz()));
  const Vec4   collinear(0.0, 0.0, mass * std::sinh(rapidity), mass * std::cosh(rapidity));
  for (auto &row : rows) {
    if (row.status == -1) { continue; }
    Vec4 p(row.px, row.py, row.pz, row.e);
    p.bstback(total);
    p.bst(collinear);
    row.px = p.px();
    row.py = p.py();
    row.pz = p.pz();
    row.e  = p.e();
  }
  for (std::size_t i = 0; i < 2; ++i) {
    const double sign = i == 0 ? 1.0 : -1.0;
    rows[i].e         = 0.5 * (collinear.e() + sign * collinear.pz());
    rows[i].pz        = sign * rows[i].e;
    rows[i].px = rows[i].py = rows[i].m = 0.0;
  }
}

// Compute the remnant partner used for a diffractive extracted parton
int RemnantPartnerPDG(int pdg) {
  if (pdg == 21) { return 21; }
  if (std::abs(pdg) >= 1 && std::abs(pdg) <= 5) { return -pdg; }
  return 90;
}

// Compute true when one active diffractive side has complete direct metadata
bool HasDirectDiffractiveSideState(const DiffSideInfo &side) {
  if (!side.active) { return true; }
  if (side.excited) { return side.system_pdg != 0 && side.system_e > 0.0 && !side.string_skeleton; }
  return side.lead_pdg != 0 && side.lead_e > 0.0;
}

// Compute true when one forward side is complete for isolated subsystem showering
bool HasIsolatedForwardSideState(const DiffSideInfo &side) {
  if (!side.active) { return true; }
  // An unresolved excitation can be restored exactly from its complete system metadata
  if (side.excited) { return side.system_e > 0.0; }
  return side.lead_pdg != 0 && side.lead_e > 0.0;
}

// Preserve the ordinary beam identity after colorless photon emission
int ProtonRemnantPDG(int pdg, int beam_pdg) {
  if (pdg == 22) { return beam_pdg; }
  throw std::invalid_argument("Colored proton remnants require Pythia fragmentation");
}

// Compute the generated mass stored for one exact four-vector
double GeneratedMass(const FourMomentumSum &p) { return std::sqrt(std::max(0.0, Mass2(p))); }

// Compute a four-momentum from one Les Houches particle row
FourMomentumSum ParticleMomentum(const LHEParticle &particle) {
  return {particle.px, particle.py, particle.pz, particle.e};
}

// Compute true when one LHE row is an exact metadata particle
bool MatchesLHERow(const LHEParticle &row, int pdg, const FourMomentumSum &target) {
  return pdg != 0 && row.id == pdg && MomentumMatches(ToHepMC(ParticleMomentum(row)), target);
}

// Compute true when one LHE row is an exact leading forward baryon
bool IsLeadingForwardRow(const LHEParticle &row, const DiffEventInfo &diff) {
  const auto matches = [&row](const DiffSideInfo &side) {
    const FourMomentumSum lead{side.lead_px, side.lead_py, side.lead_pz, side.lead_e};
    return side.active && MatchesLHERow(row, side.lead_pdg, lead);
  };
  return matches(diff.side1) || matches(diff.side2);
}

// Compute true when one LHE row is restored separately in direct-copy mode
bool IsDirectForwardRow(const LHEParticle &row, const DiffEventInfo &diff) {
  if (IsLeadingForwardRow(row, diff)) { return true; }
  const auto matches = [&row](const DiffSideInfo &side) {
    const FourMomentumSum rem{side.rem_px, side.rem_py, side.rem_pz, side.rem_e};
    return side.active && MatchesLHERow(row, side.rem_pdg, rem);
  };
  return matches(diff.side1) || matches(diff.side2);
}

// Compute the first or second incoming hard-process row when present
const LHEParticle *IncomingParticle(const LHEEventRecord &record, bool first_side) {
  if (record.particles.size() < 2) { return nullptr; }
  return &record.particles[first_side ? 0 : 1];
}

// Compute the exact ordinary proton-remnant momentum left after parton extraction
FourMomentumSum ProtonRemnantMomentum(const LHEInitInfo &init, const DiffEventInfo &diff, const LHEParticle &incoming,
                                      bool first_side) {
  return SubtractMomentum(OriginalBeamMomentum(init, diff, first_side), ParticleMomentum(incoming));
}

// Attach one integer attribute to a HepMC particle
void AddParticleIntAttribute(const HepMC3::GenParticlePtr &particle, const std::string &name, int value) {
  particle->add_attribute(name, std::make_shared<HepMC3::IntAttribute>(value));
}

// Attach one double attribute to a HepMC particle
void AddParticleDoubleAttribute(const HepMC3::GenParticlePtr &particle, const std::string &name, double value) {
  particle->add_attribute(name, std::make_shared<HepMC3::DoubleAttribute>(value));
}

// Attach LHE-style color-flow attributes to a HepMC particle
void AddColorFlowAttributes(const HepMC3::GenParticlePtr &particle, int flow1, int flow2) {
  if (flow1 != 0) { AddParticleIntAttribute(particle, "flow1", flow1); }
  if (flow2 != 0) { AddParticleIntAttribute(particle, "flow2", flow2); }
}

// Build one HepMC particle and preserve its exact input four-vector
HepMC3::GenParticlePtr MakeParticle(const FourMomentumSum &p, int pdg, int status, double generated_mass) {
  auto particle = std::make_shared<HepMC3::GenParticle>(ToHepMC(p), pdg, status);
  particle->set_generated_mass(std::max(0.0, generated_mass));
  return particle;
}

// Build one HepMC particle from a Les Houches particle row
HepMC3::GenParticlePtr MakeParticle(const LHEParticle &row, int status) {
  return MakeParticle(ParticleMomentum(row), row.id, status, std::fabs(row.m));
}

// Build an unresolved excited forward system carried only by diffraction metadata
HepMC3::GenParticlePtr MakeUnresolvedForwardSystem(const DiffSideInfo &side) {
  if (!side.active || !side.excited || side.fragmented || side.system_pdg == 0 || !(side.system_e > 0.0)) {
    return nullptr;
  }
  const FourMomentumSum system{side.system_px, side.system_py, side.system_pz, side.system_e};
  return MakeParticle(system, side.system_pdg, PDG_STABLE, side.system_m > 0.0 ? side.system_m : GeneratedMass(system));
}

// Compute true when one LHE event has no colored final-state hard particles
bool IsDirectColorSingletEvent(const LHEEventBlock &event) {
  if (!event.record.valid || event.record.particles.size() < 2 ||
      (!event.diff.active && !HasFullBeamRows(event.record))) {
    return false;
  }
  // A reduced hard event still needs Pythia to construct and hadronize its colored remnants
  if (!HasFullBeamRows(event.record) &&
      (IsColoredPDG(event.record.particles[0].id) || IsColoredPDG(event.record.particles[1].id))) {
    return false;
  }
  if (!HasDirectDiffractiveSideState(event.diff.side1) || !HasDirectDiffractiveSideState(event.diff.side2)) {
    return false;
  }
  for (const auto &particle : event.record.particles) {
    if (particle.status != 1) { continue; }
    if (particle.color1 != 0 || particle.color2 != 0 || IsColoredPDG(particle.id)) { return false; }
  }
  return true;
}

// Compute true for an event with a closed final state suitable for isolated showering
bool IsIsolatedShowerEvent(const LHEEventBlock &event, bool shower_colorless) {
  if (!event.record.valid || !HasFullBeamRows(event.record)) { return false; }
  if (event.record.particles[0].color1 != 0 || event.record.particles[0].color2 != 0 ||
      event.record.particles[1].color1 != 0 || event.record.particles[1].color2 != 0) {
    return false;
  }
  if (!HasIsolatedForwardSideState(event.diff.side1) || !HasIsolatedForwardSideState(event.diff.side2)) {
    return false;
  }

  int                               colored         = 0;
  int                               charged_leptons = 0;
  std::map<int, std::array<int, 2>> color_counts;
  for (const auto &particle : event.record.particles) {
    if (particle.status != 1) { continue; }
    if (std::abs(particle.id) == 11 || std::abs(particle.id) == 13 || std::abs(particle.id) == 15) {
      ++charged_leptons;
    }
    if (!IsColoredPDG(particle.id)) {
      if (particle.color1 != 0 || particle.color2 != 0) { return false; }
      continue;
    }
    const bool triplet = IsDiquark(particle.id) ? particle.id < 0 : particle.id > 0;
    if (particle.id == 21) {
      if (particle.color1 == 0 || particle.color2 == 0 || particle.color1 == particle.color2) { return false; }
    } else if (triplet ? particle.color1 == 0 || particle.color2 != 0 : particle.color1 != 0 || particle.color2 == 0) {
      return false;
    }
    ++colored;
    if (particle.color1 != 0) { ++color_counts[particle.color1][0]; }
    if (particle.color2 != 0) { ++color_counts[particle.color2][1]; }
  }
  if (colored == 1 || (colored > 0 && color_counts.empty())) { return false; }
  for (const auto &[tag, counts] : color_counts) {
    if (tag == 0 || counts[0] != 1 || counts[1] != 1) { return false; }
  }
  return colored >= 2 || (shower_colorless && charged_leptons >= 2);
}

// Compute true when all requested diffractive events use isolated subsystem showering
bool CanUseIsolatedShowerMode(const LHEInput &input, int max_events, bool shower_colorless) {
  int  checked    = 0;
  bool showerable = false;
  for (const auto &event : input.events) {
    if (IsIsolatedShowerEvent(event, shower_colorless)) {
      showerable = true;
    } else if (!HasFullBeamRows(event.record) || !IsDirectColorSingletEvent(event)) {
      return false;
    }
    if (++checked >= max_events) { break; }
  }
  return checked > 0 && showerable;
}

// Compute true when the requested diffractive events can bypass Pythia exactly
bool CanUseDirectColorSingletMode(const LHEInput &input, int max_events) {
  int checked = 0;
  for (const auto &event : input.events) {
    if (!IsDirectColorSingletEvent(event)) { return false; }
    if (++checked >= max_events) { break; }
  }
  return checked > 0;
}

// Compute true when every requested direct event is already a complete pp record
bool CanUseDirectFullEventMode(const LHEInput &input, int max_events) {
  int checked = 0;
  for (const auto &event : input.events) {
    if (!IsDirectColorSingletEvent(event) || !HasFullBeamRows(event.record)) { return false; }
    if (++checked >= max_events) { break; }
  }
  return checked > 0;
}

class GraniittiEventLHAup final : public LHAup {
 public:
  // Construct a one-event LHAup source with event-specific diffractive beams
  GraniittiEventLHAup(const LHEInitInfo &init, const LHEEventBlock &event)
      : LHAup(init.strategy), init_(init), event_(event) {}

  // Initialize each diffractive side as a Pomeron beam at xi times the proton energy
  bool setInit() override {
    if (!event_.record.valid) { return false; }

    const int    beam1_id = PythiaBeamId(init_, event_.diff.side1, true);
    const int    beam2_id = PythiaBeamId(init_, event_.diff.side2, false);
    const double beam1_e  = PythiaBeamEnergy(init_, event_.diff, event_.diff.side1, true);
    const double beam2_e  = PythiaBeamEnergy(init_, event_.diff, event_.diff.side2, false);
    if (!(beam1_e > 0.0) || !(beam2_e > 0.0)) { return false; }

    setBeamA(beam1_id, beam1_e, 0, 0);
    setBeamB(beam2_id, beam2_e, 0, 0);
    setStrategy(init_.strategy);

    for (const auto &line : init_.process_lines) {
      const LHEProcessInfo process = ParseLHEProcessInfo(line);
      addProcess(process.id, process.xsec, process.xerr, process.xmax);
    }
    return true;
  }

  // Supply the same source event again after a technical shower failure
  void Rewind() { consumed_ = false; }

  // Supply the stored GRANIITTI hard event to Pythia once
  bool setEvent(int = 0) override {
    if (consumed_ || !event_.record.valid) { return false; }
    consumed_ = true;

    const LHEEventRecord &record = event_.record;
    setProcess(record.idprup, record.weight, record.scalup, record.alphaEM, record.alphaS);
    for (const auto &particle : record.particles) {
      addParticle(particle.id, particle.status, particle.mother1, particle.mother2, particle.color1, particle.color2,
                  particle.px, particle.py, particle.pz, particle.e, particle.m, particle.tau, particle.spin);
    }

    if (record.particles.size() >= 2) {
      const double x1 = PythiaHardFraction(init_, event_, true);
      const double x2 = PythiaHardFraction(init_, event_, false);
      if (!(x1 > 0.0 && x1 < 1.0) || !(x2 > 0.0 && x2 < 1.0)) { return false; }
      setIdX(record.particles[0].id, record.particles[1].id, x1, x2);
    }
    return true;
  }

 private:
  const LHEInitInfo   &init_;
  const LHEEventBlock &event_;
  bool                 consumed_ = false;
};

class GraniittiDiffractionHooks final : public UserHooks {
 public:
  // Enable ISR-emission validation when a diffractive side is present
  bool canVetoISREmission() override { return true; }

  // Veto numerically invalid ISR attempts and let Pythia reshower the event
  bool doVetoISREmission(int, const Event &event, int) override { return !PythiaEventMomentaFinite(event); }

  // Enable parton-kinematics checks before constructing beam remnants
  bool canVetoPartonLevelEarly() override { return true; }

  // Veto non-finite parton-level states while preserving the hard event
  bool doVetoPartonLevelEarly(const Event &event) override { return !PythiaEventMomentaFinite(event); }

  // Request technical reshowering instead of treating hook vetoes as physics loss
  bool retryPartonLevel() override { return true; }

  // Enable a final finite-momentum validation after hadronization
  bool canVetoAfterHadronization() override { return true; }

  // Veto non-finite hadronized events before they are converted to HepMC3
  bool doVetoAfterHadronization(const Event &event) override { return !PythiaEventMomentaFinite(event); }

};

// Apply common Pythia settings shared by file and in-memory LHA modes
void ApplyCommonPythiaSettings(Pythia &pythia, const std::string &extra_cmnd, int pythia_seed) {
  if (!extra_cmnd.empty() && !pythia.readFile(extra_cmnd)) {
    throw std::invalid_argument("Cannot read Pythia settings from " + extra_cmnd);
  }

  if (pythia.settings.flag("Merging:doKTMerging")) {
    // Retain physical merging vetoes as zero weights, never reshower them until accepted
    pythia.readString("Merging:applyVeto = off");
    pythia.readString("Merging:includeWeightInXsection = off");
  }

  pythia.readString("Beams:setProductionScalesFromLHEF = on");
  pythia.readString("PartonLevel:MPI = off");
  pythia.readString("HadronLevel:all = on");
  pythia.readString("Check:event = on");
  pythia.readString("Main:timesAllowErrors = 100");
  pythia.readString("Next:numberShowEvent = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Print:quiet = on");

  if (pythia_seed > 0) {
    pythia.readString("Random:setSeed = on");
    pythia.readString("Random:seed = " + std::to_string(pythia_seed));
  }
}

// Configure one Pythia instance for a streaming LHE input file
bool ConfigureFilePythia(Pythia &pythia, const std::string &lhe_file, const std::string &extra_cmnd, int pythia_seed) {
  ApplyCommonPythiaSettings(pythia, extra_cmnd, pythia_seed);
  pythia.readString("Beams:frameType = 4");
  pythia.readString("Beams:LHEF = " + lhe_file);
  return pythia.init();
}

// Configure one Pythia instance for an in-memory GRANIITTI LHAup event
bool ConfigureExternalPythia(Pythia &pythia, const std::shared_ptr<LHAup> &lha_up,
                             const std::string &extra_cmnd, int pythia_seed) {
  ApplyCommonPythiaSettings(pythia, extra_cmnd, pythia_seed);
  pythia.readString("Beams:frameType = 5");
  pythia.readString("PartonLevel:Remnants = on");
  pythia.readString("Check:event = off");
  pythia.setLHAupPtr(lha_up);
  return pythia.init();
}

// Let Pythia collapse low-mass strings against an existing spectator hadron
class RemnantFragmentation final : public LundFragmentation {
 public:
  // Use the native hadron recoil option only for Pythia's low-mass strings
  bool fragment(int iSub, ColConfig &colors, Event &event, bool isDiff, bool systemRecoil) override {
    if (iSub >= 0 && colors[iSub].massExcess <= parm("HadronLevel:mStringMin")) {
      return ministringFragPtr->fragment(iSub, colors, event, isDiff, false);
    }
    return LundFragmentation::fragment(iSub, colors, event, isDiff, systemRecoil);
  }
};

// Install the remnant model after Pythia creates its default fragmentation models
class RemnantHooks final : public UserHooks {
 public:
  // Bind the event generator that owns the fragmentation models
  explicit RemnantHooks(Pythia &generator) : pythia(generator) {}

  // Select native cluster collapse with spectator recoil during initialization
  bool initAfterBeams() override { return pythia.setFragmentationPtr(std::make_shared<RemnantFragmentation>()); }

 private:
  Pythia &pythia;
};

// Configure isolated closed strings for FSR and hadronization without beam remnants
bool ConfigureIsolatedShowerPythia(Pythia &pythia, const std::string &extra_cmnd, int pythia_seed) {
  ApplyCommonPythiaSettings(pythia, extra_cmnd, pythia_seed);
  pythia.readString("Merging:doKTMerging = off");
  pythia.readString("ProcessLevel:all = off");
  pythia.readString("PartonLevel:ISR = off");
  pythia.readString("PartonLevel:MPI = off");
  pythia.readString("PartonLevel:Remnants = off");
  pythia.readString("Check:event = off");
  pythia.setUserHooksPtr(std::make_shared<RemnantHooks>(pythia));
  return pythia.init();
}

// Attach one integer attribute to a HepMC event
void AddIntAttribute(HepMC3::GenEvent &event, const std::string &name, int value) {
  event.add_attribute(name, std::make_shared<HepMC3::IntAttribute>(value));
}

// Attach one double attribute to a HepMC event
void AddDoubleAttribute(HepMC3::GenEvent &event, const std::string &name, double value) {
  event.add_attribute(name, std::make_shared<HepMC3::DoubleAttribute>(value));
}

// Attach the Pythia reshower attempt metadata to a HepMC event
void AddReshowerAttributes(HepMC3::GenEvent &event, int attempt, int max_attempts, int pythia_seed) {
  AddIntAttribute(event, "graniitti_reshower_attempt", attempt);
  AddIntAttribute(event, "graniitti_reshower_max_attempts", max_attempts);
  AddIntAttribute(event, "graniitti_pythia_seed", pythia_seed);
}

// Attach the final-state recoil mapping metadata to a HepMC event
void AddRecoilAttributes(HepMC3::GenEvent &event, const RecoilResult &recoil) {
  AddIntAttribute(event, "graniitti_recoil_attempted", recoil.attempted ? 1 : 0);
  AddIntAttribute(event, "graniitti_recoil_applied", recoil.applied ? 1 : 0);
  AddIntAttribute(event, "graniitti_recoil_particle_count", recoil.particle_count);
  AddDoubleAttribute(event, "graniitti_recoil_lambda", recoil.lambda);
  AddDoubleAttribute(event, "graniitti_recoil_current_mass", recoil.current_mass);
  AddDoubleAttribute(event, "graniitti_recoil_target_mass", recoil.target_mass);
  AddDoubleAttribute(event, "graniitti_recoil_mapped_mass", recoil.mapped_mass);
  AddDoubleAttribute(event, "graniitti_recoil_target_px", recoil.target_px);
  AddDoubleAttribute(event, "graniitti_recoil_target_py", recoil.target_py);
  AddDoubleAttribute(event, "graniitti_recoil_target_pz", recoil.target_pz);
  AddDoubleAttribute(event, "graniitti_recoil_target_e", recoil.target_e);
  AddDoubleAttribute(event, "graniitti_recoil_mapped_e", recoil.mapped_e);
}

// Compute the final-state rest-frame energy after one common momentum scale
double ScaledRestEnergy(const std::vector<FourMomentumSum> &rest_momenta, double scale) {
  double energy = 0.0;
  for (const FourMomentumSum &p : rest_momenta) {
    const double mass2 = std::max(0.0, Mass2(p));
    const double p2    = p.px * p.px + p.py * p.py + p.pz * p.pz;
    energy += std::sqrt(mass2 + scale * scale * p2);
  }
  return energy;
}

// Set one selected stable-particle system to an exact invariant mass
bool SetParticleListInvariantMass(const std::vector<HepMC3::GenParticlePtr> &particles, double target_mass) {
  if (particles.empty() || !(target_mass > 0.0)) { return false; }

  const FourMomentumSum current      = ParticleListMomentum(particles);
  const double          current_mass = PositiveMass(current);
  const double          scale_mass   = std::max(1.0, target_mass);
  if (std::fabs(current_mass - target_mass) < 1.0e-10 * scale_mass) { return true; }

  const std::array<double, 3>  to_rest = NegativeBeta(BetaVector(current));
  std::vector<FourMomentumSum> rest_momenta;
  rest_momenta.reserve(particles.size());
  double min_mass = 0.0;
  for (const auto &particle : particles) {
    const FourMomentumSum rest = BoostMomentum(FromHepMC(particle->momentum()), to_rest);
    rest_momenta.push_back(rest);
    min_mass += PositiveMass(rest);
  }
  if (target_mass + 1.0e-10 * scale_mass < min_mass) { return false; }

  double lo = 0.0;
  double hi = 1.0;
  while (ScaledRestEnergy(rest_momenta, hi) < target_mass) {
    hi *= 2.0;
    if (hi > 1.0e12) { return false; }
  }

  for (int iter = 0; iter < 80; ++iter) {
    const double mid = 0.5 * (lo + hi);
    if (ScaledRestEnergy(rest_momenta, mid) > target_mass) {
      hi = mid;
    } else {
      lo = mid;
    }
  }

  const double                scale  = 0.5 * (lo + hi);
  const std::array<double, 3> to_lab = BetaVector(current);
  for (std::size_t i = 0; i < particles.size(); ++i) {
    const FourMomentumSum &rest  = rest_momenta[i];
    const double           mass2 = std::max(0.0, Mass2(rest));
    const double           p2    = rest.px * rest.px + rest.py * rest.py + rest.pz * rest.pz;
    const FourMomentumSum  scaled{scale * rest.px, scale * rest.py, scale * rest.pz,
                                 std::sqrt(mass2 + scale * scale * p2)};
    particles[i]->set_momentum(ToHepMC(BoostMomentum(scaled, to_lab)));
  }
  return true;
}

// Print one recoil warning when exact final-state remapping is not possible
void ReportRecoilWarning(const RecoilResult &recoil, std::size_t event_index) {
  if (!recoil.attempted || recoil.applied) { return; }
  std::cerr << std::scientific << std::setprecision(6) << "PYTHIA_LHE_HADRONIZE recoil-warning event=" << event_index
            << " current_mass=" << recoil.current_mass << " target_mass=" << recoil.target_mass
            << " mapped_mass=" << recoil.mapped_mass << " particle_count=" << recoil.particle_count << '\n';
}

// Transform selected stable particles while preserving protected hard leptons
RecoilResult TransformParticleListToExactTarget(const std::vector<HepMC3::GenParticlePtr> &particles,
                                                const FourMomentumSum &target, const FourMomentumSum &mapped_target,
                                                std::size_t event_index) {
  RecoilResult result;
  result.attempted = true;

  result.target_px      = target.px;
  result.target_py      = target.py;
  result.target_pz      = target.pz;
  result.target_e       = target.e;
  result.target_mass    = PositiveMass(target);
  result.mapped_mass    = PositiveMass(mapped_target);
  result.mapped_e       = mapped_target.e;
  result.particle_count = static_cast<int>(particles.size());

  if (particles.empty() || !(mapped_target.e > 0.0) || !(Mass2(mapped_target) > 0.0)) {
    ReportRecoilWarning(result, event_index);
    return result;
  }

  const FourMomentumSum current = ParticleListMomentum(particles);
  result.current_mass           = PositiveMass(current);
  if (!(current.e > 0.0) || !(Mass2(current) > 0.0)) {
    ReportRecoilWarning(result, event_index);
    return result;
  }

  const std::array<double, 3> to_rest   = NegativeBeta(BetaVector(current));
  const std::array<double, 3> to_target = BetaVector(mapped_target);
  for (const auto &particle : particles) {
    const FourMomentumSum rest = BoostMomentum(FromHepMC(particle->momentum()), to_rest);
    const FourMomentumSum lab  = BoostMomentum(rest, to_target);
    particle->set_momentum(ToHepMC(lab));
  }

  result.applied = true;
  return result;
}

// Transform the non-hard stable particles into the diffractive remainder
bool TransformPythiaSubeventPreservingHardBranches(HepMC3::GenEvent &event, const LHEInitInfo &init,
                                                   const DiffEventInfo &diff, const LHEEventRecord &record,
                                                   std::size_t event_index) {
  const FourMomentumSum available    = TargetPythiaFinalMomentum(init, diff);
  const double          available_m2 = Mass2(available);
  if (!(available_m2 > 0.0)) { return false; }

  std::vector<HepMC3::GenParticlePtr> protected_branch_particles;
  int                                 hard_rows     = 0;
  int                                 matched_roots = 0;
  if (!ProtectedHardBranchParticles(event, record, protected_branch_particles, hard_rows, matched_roots)) {
    return false;
  }

  const std::vector<HepMC3::GenParticlePtr> protected_particles =
      StableFinalParticlesByProtection(event, protected_branch_particles, true);
  const std::vector<HepMC3::GenParticlePtr> remap_particles =
      StableFinalParticlesByProtection(event, protected_branch_particles, false);

  const FourMomentumSum protected_momentum = ParticleListMomentum(protected_particles);
  const FourMomentumSum remap_available    = SubtractMomentum(available, protected_momentum);
  const double          remap_available_m2 = Mass2(remap_available);
  if (!(remap_available.e > 0.0) || !(remap_available_m2 > 0.0)) { return false; }

  AddIntAttribute(event, "graniitti_recoil_protected_particles", static_cast<int>(protected_particles.size()));
  AddIntAttribute(event, "graniitti_recoil_remapped_particles", static_cast<int>(remap_particles.size()));
  AddIntAttribute(event, "graniitti_hard_rows", hard_rows);
  AddIntAttribute(event, "graniitti_hard_roots_matched", matched_roots);
  AddDoubleAttribute(event, "graniitti_protected_mass_before_recoil", PositiveMass(protected_momentum));

  if (!SetParticleListInvariantMass(remap_particles, std::sqrt(remap_available_m2))) { return false; }

  const FourMomentumSum current = ParticleListMomentum(remap_particles);
  if (!(current.e > 0.0) || !(Mass2(current) > 0.0)) { return false; }

  const RecoilResult recoil =
      TransformParticleListToExactTarget(remap_particles, remap_available, remap_available, event_index);
  const FourMomentumSum protected_after = ParticleListMomentum(protected_particles);
  AddDoubleAttribute(event, "graniitti_protected_mass_after_recoil", PositiveMass(protected_after));
  AddRecoilAttributes(event, recoil);
  return recoil.applied;
}

// Attach the total event four-momentum closure residual to a HepMC event
void AddClosureAttributes(HepMC3::GenEvent &event, const ClosureResult &closure) {
  AddDoubleAttribute(event, "graniitti_closure_initial_e", closure.initial.e);
  AddDoubleAttribute(event, "graniitti_closure_final_e", closure.final.e);
  AddDoubleAttribute(event, "graniitti_closure_dpx", closure.delta.px);
  AddDoubleAttribute(event, "graniitti_closure_dpy", closure.delta.py);
  AddDoubleAttribute(event, "graniitti_closure_dpz", closure.delta.pz);
  AddDoubleAttribute(event, "graniitti_closure_de", closure.delta.e);
  AddDoubleAttribute(event, "graniitti_closure_max_abs", closure.max_abs);
  AddDoubleAttribute(event, "graniitti_closure_tolerance", closure.tolerance);
  AddIntAttribute(event, "graniitti_closure_pass", closure.pass ? 1 : 0);
}

// Print one closure warning when the total pp event does not close
void ReportClosureWarning(const ClosureResult &closure, std::size_t event_index) {
  if (closure.pass) { return; }
  std::cerr << std::scientific << std::setprecision(6) << "PYTHIA_LHE_HADRONIZE closure-warning event=" << event_index
            << " dpx=" << closure.delta.px << " dpy=" << closure.delta.py << " dpz=" << closure.delta.pz
            << " de=" << closure.delta.e << " max_abs=" << closure.max_abs << " tolerance=" << closure.tolerance
            << '\n';
}

// Check and annotate total pp four-momentum closure after proton reinsertion
bool CheckAndAnnotateClosure(HepMC3::GenEvent &event, const LHEInitInfo &init, const DiffEventInfo &diff,
                             std::size_t event_index) {
  const ClosureResult closure = MomentumClosure(event, init, diff);
  AddClosureAttributes(event, closure);
  ReportClosureWarning(closure, event_index);
  return closure.pass;
}

// Compute particles emitted by the converter's implicit source vertex
std::vector<HepMC3::GenParticlePtr> ImplicitSourceOutgoingParticles(HepMC3::GenEvent &event) {
  std::vector<HepMC3::GenParticlePtr> particles;
  for (const auto &particle : event.particles()) {
    if (!particle) { continue; }
    HepMC3::GenVertexPtr vertex = particle->production_vertex();
    if (vertex && vertex->particles_in().empty()) { particles.push_back(particle); }
  }
  return particles;
}

// Create one original incoming beam particle for the HepMC source vertex
HepMC3::GenParticlePtr OriginalBeamParticle(const LHEInitInfo &init, const DiffEventInfo &diff, bool first_side) {
  const int             beam_id  = first_side ? init.beam1_id : init.beam2_id;
  const DiffSideInfo   &side     = first_side ? diff.side1 : diff.side2;
  const FourMomentumSum beam     = OriginalBeamMomentum(init, diff, first_side);
  auto                  particle = std::make_shared<HepMC3::GenParticle>(ToHepMC(beam), beam_id, 4);
  const double          mass     = first_side ? init.beam1_mass : init.beam2_mass;
  particle->set_generated_mass(mass >= 0.0 ? mass : BeamMassGeV(beam_id, side));
  return particle;
}

// Add the original pp beam particles as real incoming HepMC particles
void AddOriginalBeamParticles(HepMC3::GenEvent &event, const LHEInitInfo &init, const DiffEventInfo &diff) {
  const std::vector<HepMC3::GenParticlePtr> root_particles = ImplicitSourceOutgoingParticles(event);
  auto                                      source         = std::make_shared<HepMC3::GenVertex>();
  auto                                      beam1          = OriginalBeamParticle(init, diff, true);
  auto                                      beam2          = OriginalBeamParticle(init, diff, false);
  source->add_particle_in(beam1);
  source->add_particle_in(beam2);
  for (const auto &particle : root_particles) {
    if (particle) { source->add_particle_out(particle); }
  }
  // Remove orphan vertices emptied when their particles join the beam vertex
  const auto vertices = event.vertices();
  for (const auto &vertex : vertices) {
    if (vertex->particles_in().empty() && vertex->particles_out().empty()) { event.remove_vertex(vertex); }
  }
  event.add_vertex(source);
  event.set_beam_particles(beam1, beam2);
}

// Restore or identify the tagged proton and retain its source momentum and recoil
void AddLeadingProton(HepMC3::GenEvent &event, const DiffSideInfo &side, int side_index) {
  if (!side.active || side.lead_pdg == 0) { return; }

  const FourMomentumSum original{side.lead_px, side.lead_py, side.lead_pz, side.lead_e};
  HepMC3::GenParticlePtr particle;
  for (const auto &root : event.particles()) {
    if (root->pid() != side.lead_pdg || !MomentumMatches(root->momentum(), original)) { continue; }
    std::vector<HepMC3::GenParticlePtr> stable, visited;
    CollectStableDescendants(root, stable, visited);
    for (const auto &child : stable) {
      if (child->pid() == side.lead_pdg) { particle = child; break; }
    }
    if (particle) { break; }
  }
  if (!particle) {
    particle = MakeParticle(original, side.lead_pdg, PDG_STABLE, side.lead_m);
    event.add_particle(particle);
  }
  const auto recoil = SubtractMomentum(FromHepMC(particle->momentum()), original);
  AddParticleDoubleAttribute(particle, "graniitti_recoil_dpx", recoil.px);
  AddParticleDoubleAttribute(particle, "graniitti_recoil_dpy", recoil.py);
  AddParticleDoubleAttribute(particle, "graniitti_recoil_dpz", recoil.pz);
  AddParticleDoubleAttribute(particle, "graniitti_recoil_de", recoil.e);

  particle->add_attribute("graniitti_diff_side", std::make_shared<HepMC3::IntAttribute>(side_index));
  AddParticleDoubleAttribute(particle, "graniitti_xi", side.xi);
  AddParticleDoubleAttribute(particle, "graniitti_beta", side.beta);
  AddParticleDoubleAttribute(particle, "graniitti_t", side.t);
}

// Attach GRANIITTI diffractive scalar metadata to a HepMC event
void AnnotateDiffractiveMetadata(HepMC3::GenEvent &event, const DiffEventInfo &diff) {
  if (!diff.active) { return; }

  AddIntAttribute(event, "graniitti_diff_side_mask", diff.side_mask);
  AddIntAttribute(event, "graniitti_diff_id1", diff.id1);
  AddIntAttribute(event, "graniitti_diff_id2", diff.id2);
  AddDoubleAttribute(event, "graniitti_diff_xhard1", diff.xhard1);
  AddDoubleAttribute(event, "graniitti_diff_xhard2", diff.xhard2);
  AddDoubleAttribute(event, "graniitti_diff_beam1_e", diff.beam1_e);
  AddDoubleAttribute(event, "graniitti_diff_beam2_e", diff.beam2_e);
  AddDoubleAttribute(event, "graniitti_diff_side1_xi", diff.side1.xi);
  AddDoubleAttribute(event, "graniitti_diff_side1_beta", diff.side1.beta);
  AddDoubleAttribute(event, "graniitti_diff_side1_t", diff.side1.t);
  AddDoubleAttribute(event, "graniitti_diff_side2_xi", diff.side2.xi);
  AddDoubleAttribute(event, "graniitti_diff_side2_beta", diff.side2.beta);
  AddDoubleAttribute(event, "graniitti_diff_side2_t", diff.side2.t);
}

// Attach GRANIITTI diffractive metadata to the post-Pythia HepMC event
void AnnotateDiffractiveEvent(HepMC3::GenEvent &event, const LHEInitInfo &init, const DiffEventInfo &diff) {
  AnnotateDiffractiveMetadata(event, diff);
  AddLeadingProton(event, diff.side1, 1);
  AddLeadingProton(event, diff.side2, 2);
  AddOriginalBeamParticles(event, init, diff);
}

// Attach one string attribute to a HepMC event
void AddStringAttribute(HepMC3::GenEvent &event, const std::string &name, const std::string &value) {
  event.add_attribute(name, std::make_shared<HepMC3::StringAttribute>(value));
}

// Attach LHE process cross-section metadata to a HepMC event
void AddCrossSectionInfo(HepMC3::GenEvent &event, const LHEInitInfo &init, long accepted_count = 0,
                         long attempted_count = 0) {
  double xsec = 0.0, variance = 0.0;
  for (const auto &line : init.process_lines) {
    const LHEProcessInfo process = ParseLHEProcessInfo(line);
    xsec += process.xsec;
    variance += process.xerr * process.xerr;
  }
  auto cross_section = std::make_shared<HepMC3::GenCrossSection>();
  cross_section->set_cross_section(xsec, std::sqrt(variance));
  if (accepted_count > 0) { cross_section->set_accepted_events(accepted_count); }
  if (attempted_count > 0) { cross_section->set_attempted_events(attempted_count); }
  event.set_cross_section(cross_section);
}

// Attach hard-process PDF metadata to a direct HepMC event
void AddPdfInfo(HepMC3::GenEvent &event, const LHEEventBlock &block) {
  if (!block.record.valid || block.record.particles.size() < 2) { return; }
  if (!(block.diff.xhard1 > 0.0) || !(block.diff.xhard2 > 0.0)) { return; }

  auto pdf_info = std::make_shared<HepMC3::GenPdfInfo>();
  pdf_info->set(block.diff.id1, block.diff.id2, block.diff.xhard1, block.diff.xhard2, block.record.muF, 1.0, 1.0);
  event.set_pdf_info(pdf_info);
}

// Preserve the LHE generator scales in every converted HepMC event
void AddGeneratorScaleInfo(HepMC3::GenEvent &event, const LHEEventRecord &record) {
  if (record.muF > 0.0 && std::isfinite(record.muF)) { AddDoubleAttribute(event, "graniitti_mu_f", record.muF); }
  if (record.muR > 0.0 && std::isfinite(record.muR)) { AddDoubleAttribute(event, "graniitti_mu_r", record.muR); }
  if (record.scalup > 0.0 && std::isfinite(record.scalup)) {
    AddDoubleAttribute(event, "graniitti_scalup", record.scalup);
  }
}

// Attach the event weight from one LHE record
void AddEventWeight(HepMC3::GenEvent &event, const LHEEventRecord &record) {
  event.weights().clear();
  event.weights().push_back(record.weight);
  AddDoubleAttribute(event, "graniitti_lhe_weight", record.weight);
  AddGeneratorScaleInfo(event, record);
}

// Add a diffractive remnant particle from GRANIITTI metadata
void AddDiffractiveRemnant(HepMC3::GenVertexPtr source, const DiffSideInfo &side, const LHEParticle &incoming,
                           int side_index) {
  if (!side.active || !(side.rem_e > 0.0)) { return; }

  FourMomentumSum remnant{side.rem_px, side.rem_py, side.rem_pz, side.rem_e};
  const int       pdg      = (side.rem_pdg != 0) ? side.rem_pdg : RemnantPartnerPDG(incoming.id);
  auto            particle = MakeParticle(remnant, pdg, PDG_STABLE, GeneratedMass(remnant));
  source->add_particle_out(particle);
}

// Add a normal proton-remnant particle from exact beam minus parton momentum
void AddProtonRemnant(HepMC3::GenVertexPtr source, const LHEInitInfo &init, const DiffEventInfo &diff,
                      const LHEParticle &incoming, bool first_side) {
  const int             beam_pdg = first_side ? init.beam1_id : init.beam2_id;
  const FourMomentumSum remnant  = ProtonRemnantMomentum(init, diff, incoming, first_side);
  auto particle = MakeParticle(remnant, ProtonRemnantPDG(incoming.id, beam_pdg), PDG_STABLE, GeneratedMass(remnant));
  source->add_particle_out(particle);
}

// Add one exact leading proton from GRANIITTI metadata to the direct source vertex
void AddDirectLeadingProton(HepMC3::GenVertexPtr source, const DiffSideInfo &side, int side_index) {
  if (!side.active || side.lead_pdg == 0) { return; }

  FourMomentumSum lead{side.lead_px, side.lead_py, side.lead_pz, side.lead_e};
  auto particle = MakeParticle(lead, side.lead_pdg, PDG_STABLE, side.lead_m > 0.0 ? side.lead_m : GeneratedMass(lead));
  source->add_particle_out(particle);
}

// Encode one ground or isomeric nuclear PDG identity
int NuclearPDG(const int a, const int z, const int isomer) {
  if (a < 1 || z < 0 || z > a || isomer < 0 || isomer > 9) {
    throw std::invalid_argument("invalid nuclear identity for Pythia decay");
  }
  return 1000000000 + 10000 * z + 10 * a + isomer;
}

// Configure one event-local one-neutron Pythia decay table
void ConfigureNuclearDecay(Pythia &pythia, const DiffSideInfo &side, int parent_pdg, int residual_pdg,
                           double residual_mass) {
  const int         parent_abs   = std::abs(parent_pdg);
  const int         residual_abs = std::abs(residual_pdg);
  const std::string parent_name  = "A" + std::to_string(side.nuclear_a) + "Z" + std::to_string(side.nuclear_z) + "star";
  const std::string residual_name = "A" + std::to_string(side.nuclear_a - 1) + "Z" + std::to_string(side.nuclear_z);
  pythia.particleData.addParticle(residual_abs, residual_name, residual_name + "bar", 1, 3 * side.nuclear_z, 0,
                                  residual_mass);
  pythia.particleData.isResonance(residual_abs, false);
  pythia.particleData.mayDecay(residual_abs, false);
  pythia.particleData.addParticle(parent_abs, parent_name, parent_name + "bar", 1, 3 * side.nuclear_z, 0,
                                  NuclearExcitedMass(side));
  pythia.particleData.isResonance(parent_abs, false);
  auto parent = pythia.particleData.findParticle(parent_abs);
  if (parent == nullptr) { throw std::runtime_error("Pythia nuclear parent data is unavailable"); }
  parent->clearChannels();
  const int neutron_pdg = parent_pdg > 0 ? 2112 : -2112;
  parent->addChannel(1, 1.0, 0, residual_pdg, neutron_pdg);
  pythia.particleData.mayDecay(parent_abs, true);
}

// Fragment one Xn container through an isotropic Pythia one-neutron decay
bool AddNuclearDecay(Pythia &pythia, HepMC3::GenEvent &event, const HepMC3::GenVertexPtr &source, DiffSideInfo &side) {
  if (!NuclearBreakup(side)) { return false; }
  constexpr double neutron_mass  = 0.9395654133;
  const int        sign          = side.nuclear_parent_pdg > 0 ? 1 : -1;
  const int        parent_pdg    = sign * NuclearPDG(side.nuclear_a, side.nuclear_z, 1);
  const int        residual_pdg  = sign * NuclearPDG(side.nuclear_a - 1, side.nuclear_z, 0);
  const double     residual_mass = side.system_m - neutron_mass + side.neutron_threshold;
  if (!(residual_mass > 0.0) || !(NuclearExcitedMass(side) > residual_mass + neutron_mass)) { return false; }
  ConfigureNuclearDecay(pythia, side, parent_pdg, residual_pdg, residual_mass);

  pythia.event.reset();
  const int parent_index = pythia.event.append(parent_pdg, 23, 0, 0, side.system_px, side.system_py, side.system_pz,
                                               side.system_e, NuclearExcitedMass(side));
  if (!pythia.moreDecays(parent_index, false)) { return false; }
  const int first = pythia.event[parent_index].daughter1();
  const int last  = pythia.event[parent_index].daughter2();
  if (first <= 0 || last < first) { return false; }

  const FourMomentumSum parent_momentum{side.system_px, side.system_py, side.system_pz, side.system_e};
  auto                  parent = MakeParticle(parent_momentum, parent_pdg, PDG_INTERMEDIATE, NuclearExcitedMass(side));
  source->add_particle_out(parent);
  AddParticleIntAttribute(parent, "graniitti_upc_parent_pdg", side.nuclear_parent_pdg);
  AddParticleIntAttribute(parent, "graniitti_upc_a", side.nuclear_a);
  AddParticleIntAttribute(parent, "graniitti_upc_z", side.nuclear_z);
  AddParticleDoubleAttribute(parent, "graniitti_upc_neutron_threshold", side.neutron_threshold);

  auto decay = std::make_shared<HepMC3::GenVertex>();
  decay->add_particle_in(parent);
  event.add_vertex(decay);
  int neutrons = 0;
  for (int i = first; i <= last; ++i) {
    const Particle       &daughter = pythia.event[i];
    const FourMomentumSum momentum{daughter.px(), daughter.py(), daughter.pz(), daughter.e()};
    auto                  particle = MakeParticle(momentum, daughter.id(), PDG_STABLE, daughter.m());
    decay->add_particle_out(particle);
    if (std::abs(daughter.id()) == 2112) {
      ++neutrons;
      AddParticleIntAttribute(particle, "graniitti_upc_neutron", 1);
    } else {
      AddParticleIntAttribute(particle, "graniitti_upc_residual", 1);
    }
    AddParticleIntAttribute(particle, "graniitti_upc_parent_pdg", side.nuclear_parent_pdg);
  }
  if (neutrons != 1) { return false; }
  AddParticleIntAttribute(parent, "graniitti_upc_neutron_multiplicity", neutrons);
  side.fragmented = true;
  return true;
}

// Add the complete forward-side final state for one direct diffractive event
bool AddDirectForwardSide(Pythia *pythia, HepMC3::GenEvent &event, HepMC3::GenVertexPtr source, const LHEInitInfo &init,
                          DiffEventInfo &diff, const LHEParticle &incoming, bool first_side) {
  DiffSideInfo &side = first_side ? diff.side1 : diff.side2;
  if (side.active) {
    if (NuclearBreakup(side)) { return pythia != nullptr && AddNuclearDecay(*pythia, event, source, side); }
    if (side.excited && !side.fragmented) {
      const FourMomentumSum system{side.system_px, side.system_py, side.system_pz, side.system_e};
      auto                  particle = MakeParticle(system, side.system_pdg, PDG_STABLE,
                                   side.system_m > 0.0 ? side.system_m : GeneratedMass(system));
      source->add_particle_out(particle);
      return true;
    }
    AddDirectLeadingProton(source, side, first_side ? 1 : 2);
    AddDiffractiveRemnant(source, side, incoming, first_side ? 1 : 2);
  } else {
    AddProtonRemnant(source, init, diff, incoming, first_side);
  }
  return true;
}

// Add central final-state particles from the LHE hard process
void AddDirectHardFinalState(HepMC3::GenVertexPtr hard, const LHEEventRecord &record, const DiffEventInfo &diff) {
  for (const auto &row : record.particles) {
    if (row.status != 1 || IsDirectForwardRow(row, diff)) { continue; }
    hard->add_particle_out(MakeParticle(row, PDG_STABLE));
  }
}

// Compute true when one HepMC particle matches an expected four-vector
bool MatchesParticle(const HepMC3::GenParticlePtr &particle, int pdg, const FourMomentumSum &target) {
  if (!particle || particle->pid() != pdg) { return false; }
  return MomentumMatches(particle->momentum(), target);
}

// Attach invariant-mass bookkeeping to all particles after HepMC assigns ids
void AnnotateDirectParticleMasses(HepMC3::GenEvent &event) {
  for (const auto &particle : event.particles()) {
    if (!particle) { continue; }
    AddParticleDoubleAttribute(particle, "graniitti_mass2", Mass2(FromHepMC(particle->momentum())));
  }
}

// Attach color-flow tags to the direct hard incoming partons
void AnnotateDirectHardPartonColors(HepMC3::GenEvent &event, const LHEEventRecord &record) {
  std::vector<HepMC3::GenParticlePtr> incoming;
  for (const auto &particle : event.particles()) {
    if (particle && particle->status() == 21) { incoming.push_back(particle); }
  }
  if (incoming.size() < 2 || record.particles.size() < 2) { return; }
  AddColorFlowAttributes(incoming[0], record.particles[0].color1, record.particles[0].color2);
  AddColorFlowAttributes(incoming[1], record.particles[1].color1, record.particles[1].color2);
}

// Attach metadata and color flow to a matched active diffractive remnant
void AnnotateDirectDiffractiveRemnant(const HepMC3::GenParticlePtr &particle, const DiffSideInfo &side,
                                      const LHEParticle &incoming, int side_index) {
  AddParticleIntAttribute(particle, "graniitti_remnant", 1);
  AddParticleIntAttribute(particle, "graniitti_diff_side", side_index);
  AddParticleDoubleAttribute(particle, "graniitti_xi", side.xi);
  AddParticleDoubleAttribute(particle, "graniitti_beta", side.beta);
  AddParticleDoubleAttribute(particle, "graniitti_t", side.t);
  if (incoming.color1 != 0) { AddColorFlowAttributes(particle, 0, incoming.color1); }
  if (incoming.color2 != 0) { AddColorFlowAttributes(particle, incoming.color2, 0); }
}

// Attach metadata to a matched active leading proton
void AnnotateDirectLeadingProton(const HepMC3::GenParticlePtr &particle, const DiffSideInfo &side, int side_index) {
  AddParticleIntAttribute(particle, "graniitti_diff_side", side_index);
  AddParticleDoubleAttribute(particle, "graniitti_xi", side.xi);
  AddParticleDoubleAttribute(particle, "graniitti_beta", side.beta);
  AddParticleDoubleAttribute(particle, "graniitti_t", side.t);
}

// Attach metadata and color flow to a matched ordinary proton remnant
void AnnotateDirectProtonRemnant(const HepMC3::GenParticlePtr &particle, const LHEParticle &incoming, int side_index) {
  AddParticleIntAttribute(particle, "graniitti_remnant", 1);
  AddParticleIntAttribute(particle, "graniitti_diff_side", side_index);
  if (incoming.color1 != 0) { AddColorFlowAttributes(particle, 0, incoming.color1); }
  if (incoming.color2 != 0) { AddColorFlowAttributes(particle, incoming.color2, 0); }
}

// Attach direct final-state metadata for one incoming side
void AnnotateDirectForwardSide(HepMC3::GenEvent &event, const LHEInitInfo &init, const DiffEventInfo &diff,
                               const LHEParticle &incoming, bool first_side) {
  const DiffSideInfo &side       = first_side ? diff.side1 : diff.side2;
  const int           side_index = first_side ? 1 : 2;

  if (side.active) {
    FourMomentumSum lead{side.lead_px, side.lead_py, side.lead_pz, side.lead_e};
    FourMomentumSum remnant{side.rem_px, side.rem_py, side.rem_pz, side.rem_e};
    const int       remnant_pdg = (side.rem_pdg != 0) ? side.rem_pdg : RemnantPartnerPDG(incoming.id);
    for (const auto &particle : event.particles()) {
      if (MatchesParticle(particle, side.lead_pdg, lead)) { AnnotateDirectLeadingProton(particle, side, side_index); }
      if (MatchesParticle(particle, remnant_pdg, remnant)) {
        AnnotateDirectDiffractiveRemnant(particle, side, incoming, side_index);
      }
    }
    return;
  }

  const int             beam_pdg    = first_side ? init.beam1_id : init.beam2_id;
  const FourMomentumSum remnant     = ProtonRemnantMomentum(init, diff, incoming, first_side);
  const int             remnant_pdg = ProtonRemnantPDG(incoming.id, beam_pdg);
  for (const auto &particle : event.particles()) {
    if (MatchesParticle(particle, remnant_pdg, remnant)) {
      AnnotateDirectProtonRemnant(particle, incoming, side_index);
    }
  }
}

// Attach all direct particle-level metadata after event insertion
void AnnotateDirectParticles(HepMC3::GenEvent &event, const LHEInput &input, const LHEEventBlock &block) {
  const LHEParticle *incoming1 = IncomingParticle(block.record, true);
  const LHEParticle *incoming2 = IncomingParticle(block.record, false);
  if (incoming1 == nullptr || incoming2 == nullptr) { return; }
  AnnotateDirectParticleMasses(event);
  AnnotateDirectHardPartonColors(event, block.record);
  AnnotateDirectForwardSide(event, input.init, block.diff, *incoming1, true);
  AnnotateDirectForwardSide(event, input.init, block.diff, *incoming2, false);
}

// Build one full-record particle with a remapped intact forward momentum
HepMC3::GenParticlePtr MakeRemappedFullParticle(const LHEParticle &row, const DiffEventInfo &original,
                                                const DiffEventInfo &mapped) {
  const auto remap = [&row](const DiffSideInfo &before, const DiffSideInfo &after) -> HepMC3::GenParticlePtr {
    const FourMomentumSum lead{before.lead_px, before.lead_py, before.lead_pz, before.lead_e};
    if (!before.active || before.excited || !MatchesLHERow(row, before.lead_pdg, lead)) { return nullptr; }
    return MakeParticle(ForwardSideMomentum(after), row.id, PDG_STABLE,
                        after.lead_m > 0.0 ? after.lead_m : PositiveMass(ForwardSideMomentum(after)));
  };
  if (auto particle = remap(original.side1, mapped.side1)) { return particle; }
  if (auto particle = remap(original.side2, mapped.side2)) { return particle; }
  return MakeParticle(row, PDG_STABLE);
}

// Build one direct HepMC event for a color-singlet diffractive LHE record
bool BuildDirectColorSingletEvent(const LHEInput &input, const LHEEventBlock &block, std::size_t event_index,
                                  Pythia *nuclear_pythia, HepMC3::GenEvent &event) {
  if (!IsDirectColorSingletEvent(block)) { return false; }
  const LHEParticle *incoming1 = IncomingParticle(block.record, true);
  const LHEParticle *incoming2 = IncomingParticle(block.record, false);
  if (incoming1 == nullptr || incoming2 == nullptr) { return false; }

  DiffEventInfo diff = block.diff;
  if (!PrepareNuclearPair(diff)) { return false; }

  event = HepMC3::GenEvent(HepMC3::Units::GEV, HepMC3::Units::MM);
  event.set_event_number(static_cast<int>(event_index + 1));

  auto source = std::make_shared<HepMC3::GenVertex>();
  auto hard   = std::make_shared<HepMC3::GenVertex>();
  auto beam1  = OriginalBeamParticle(input.init, diff, true);
  auto beam2  = OriginalBeamParticle(input.init, diff, false);
  source->add_particle_in(beam1);
  source->add_particle_in(beam2);
  event.set_beam_particles(beam1, beam2);

  // Full pp records are copied verbatim instead of rebuilding an effective hard vertex
  if (HasFullBeamRows(block.record)) {
    for (std::size_t i = 2; i < block.record.particles.size(); ++i) {
      const LHEParticle &row = block.record.particles[i];
      if (row.status == 1) { source->add_particle_out(MakeRemappedFullParticle(row, block.diff, diff)); }
    }
    // Add unfragmented N-stars omitted from the explicit LHE rows
    const auto add_unresolved = [&source](const DiffSideInfo &side) {
      if (side.nuclear) { return; }
      const auto particle = MakeUnresolvedForwardSystem(side);
      if (particle != nullptr) { source->add_particle_out(particle); }
    };
    add_unresolved(diff.side1);
    add_unresolved(diff.side2);

    event.add_vertex(source);
    if (NuclearBreakup(diff.side1) &&
        (nuclear_pythia == nullptr || !AddNuclearDecay(*nuclear_pythia, event, source, diff.side1))) {
      return false;
    }
    if (NuclearBreakup(diff.side2) &&
        (nuclear_pythia == nullptr || !AddNuclearDecay(*nuclear_pythia, event, source, diff.side2))) {
      return false;
    }
    AddEventWeight(event, block.record);
    AddCrossSectionInfo(event, input.init, diff.accepted_events, diff.attempted_events);
    AddPdfInfo(event, block);
    AddStringAttribute(event, "graniitti_converter_mode", "direct-full-event");
    AnnotateDiffractiveMetadata(event, diff);
    if (!CheckAndAnnotateClosure(event, input.init, diff, event_index + 1)) { return false; }
    return true;
  }

  auto hard1 = MakeParticle(*incoming1, 21);
  auto hard2 = MakeParticle(*incoming2, 21);
  source->add_particle_out(hard1);
  source->add_particle_out(hard2);
  hard->add_particle_in(hard1);
  hard->add_particle_in(hard2);

  event.add_vertex(source);
  if (!AddDirectForwardSide(nuclear_pythia, event, source, input.init, diff, *incoming1, true) ||
      !AddDirectForwardSide(nuclear_pythia, event, source, input.init, diff, *incoming2, false)) {
    return false;
  }
  AddDirectHardFinalState(hard, block.record, diff);

  event.add_vertex(hard);
  AnnotateDirectParticles(event, input, block);
  AddEventWeight(event, block.record);
  AddCrossSectionInfo(event, input.init, diff.accepted_events, diff.attempted_events);
  AddPdfInfo(event, block);
  AddStringAttribute(event, "graniitti_converter_mode", "direct-copy");
  AnnotateDiffractiveMetadata(event, diff);
  if (!CheckAndAnnotateClosure(event, input.init, diff, event_index + 1)) { return false; }
  return true;
}

struct IsolatedShowerRange {
  int    first = 0;
  int    last  = 0;
  int    count = 0;
  double scale = 0.0;
};

// Compute true for a colorless diffractive remnant hadron that Pythia should decay
bool IsPythiaRemnantHadron(const LHEParticle &particle, const DiffEventInfo &diff) {
  const int id = std::abs(particle.id);
  return diff.active && particle.status == 1 && id > 100 && id < 1000000000 && !IsColoredPDG(id) &&
         !IsLeadingForwardRow(particle, diff);
}

// Choose a physical shower scale from the hard process or forward-system masses
double IsolatedShowerScale(const LHEEventBlock &block) {
  double     scale        = block.record.scalup;
  const auto include_side = [&scale](const DiffSideInfo &side) {
    if (side.active && side.excited && side.string_skeleton) { scale = std::max(scale, side.system_m); }
  };
  include_side(block.diff.side1);
  include_side(block.diff.side2);
  if (scale > 0.0 && std::isfinite(scale)) { return scale; }

  FourMomentumSum colored;
  for (const auto &particle : block.record.particles) {
    if (particle.status == 1 && IsColoredPDG(particle.id)) {
      colored = AddMomentum(colored, ParticleMomentum(particle));
    }
  }
  return PositiveMass(colored);
}

// Fill one complete non-leading final state into a fresh Pythia event
bool FillIsolatedShowerParticles(Pythia &pythia, const LHEEventBlock &block, bool shower_colorless,
                                 IsolatedShowerRange &range) {
  const LHEEventRecord &record = block.record;
  const double          scale  = IsolatedShowerScale(block);
  if (!(scale > 0.0) || !std::isfinite(scale)) { return false; }

  pythia.event.reset();
  pythia.event.scale(scale);
  range       = IsolatedShowerRange();
  range.scale = scale;
  for (const auto &particle : record.particles) {
    if (particle.status != 1) { continue; }
    if (IsLeadingForwardRow(particle, block.diff)) {
      // Status 81 makes the on-shell spectator available to Pythia cluster collapse
      pythia.event.append(particle.id, 81, 0, 0, particle.px, particle.py, particle.pz, particle.e,
                           std::max(0.0, particle.m));
      continue;
    }
    if (!shower_colorless && !IsColoredPDG(particle.id) && !IsPythiaRemnantHadron(particle, block.diff)) { continue; }
    const int index = pythia.event.append(particle.id, 23, particle.color1, particle.color2, particle.px, particle.py,
                                          particle.pz, particle.e, std::max(0.0, particle.m), scale, particle.spin);
    if (range.count == 0) { range.first = index; }
    range.last = index;
    ++range.count;
  }
  return range.count >= 2;
}

// Shower and hadronize one event containing isolated closed color systems
bool ShowerIsolatedEvent(Pythia &pythia, const LHEEventBlock &block, bool shower_colorless) {
  IsolatedShowerRange range;
  if (!FillIsolatedShowerParticles(pythia, block, shower_colorless, range)) { return false; }
  if (pythia.flag("PartonLevel:FSR")) { pythia.forceTimeShower(range.first, range.last, range.scale); }
  return pythia.next() && PythiaEventMomentaFinite(pythia.event);
}

// Restore exact colorless LHE final states excluded from the isolated Pythia system
void AddIsolatedColorlessFinalState(HepMC3::GenEvent &event, const LHEEventBlock &block) {
  for (const auto &particle : block.record.particles) {
    if (particle.status != 1 || IsColoredPDG(particle.id) || IsLeadingForwardRow(particle, block.diff) ||
        IsPythiaRemnantHadron(particle, block.diff)) {
      continue;
    }
    event.add_particle(MakeParticle(particle, PDG_STABLE));
  }
}

// Restore unresolved excited forward systems carried only by diffraction metadata
void AddIsolatedUnresolvedForwardSystems(HepMC3::GenEvent &event, const LHEEventBlock &block) {
  const auto add_side = [&event](const DiffSideInfo &side) {
    const auto particle = MakeUnresolvedForwardSystem(side);
    if (particle != nullptr) { event.add_particle(particle); }
  };
  add_side(block.diff.side1);
  add_side(block.diff.side2);
}

// Convert one accepted isolated shower and restore exact forward baryons
bool BuildIsolatedShowerHepMC(Pythia &pythia, const LHEInput &input, const LHEEventBlock &block, int event_number,
                              int attempt, int max_attempts, int seed, bool shower_colorless,
                              HepMC3::GenEvent &output) {
  HepMC3::Pythia8ToHepMC3 converter;
  if (!converter.fill_next_event(pythia.event, output, event_number)) { return false; }
  if (!shower_colorless) { AddIsolatedColorlessFinalState(output, block); }
  AddIsolatedUnresolvedForwardSystems(output, block);

  AddEventWeight(output, block.record);
  AddCrossSectionInfo(output, input.init, block.diff.accepted_events, block.diff.attempted_events);
  AddPdfInfo(output, block);
  AddReshowerAttributes(output, attempt, max_attempts, seed);
  AddStringAttribute(output, "graniitti_converter_mode", "isolated-shower");
  AddDoubleAttribute(output, "graniitti_shower_scale", IsolatedShowerScale(block));
  AddDoubleAttribute(output, "graniitti_alpha_qcd", block.record.alphaS);
  AnnotateDiffractiveEvent(output, input.init, block.diff);
  if (!CheckAndAnnotateClosure(output, input.init, block.diff, static_cast<std::size_t>(event_number))) {
    return false;
  }
  return true;
}

// Hadronize one complete source event without dropping its weight
void ShowerIsolatedEventOrThrow(Pythia &pythia, const LHEInput &input, const LHEEventBlock &block,
                                HepMC3::GenEvent &output, int event_number, int seed, int max_attempts,
                                bool shower_colorless) {
  for (int attempt = 1; attempt <= max_attempts; ++attempt) {
    output.clear();
    if (ShowerIsolatedEvent(pythia, block, shower_colorless) &&
        BuildIsolatedShowerHepMC(pythia, input, block, event_number, attempt, max_attempts, seed,
                                 shower_colorless, output)) { return; }
  }
  pythia.stat();
  throw std::runtime_error("Pythia rejected isolated event " + std::to_string(event_number) + " after " +
                           std::to_string(max_attempts) + " attempts");
}

// Run closed color systems through isolated FSR and hadronization
int RunIsolatedShower(const LHEInput &input, const std::string &hepmc_file, int max_events, int seed,
                      const std::string &extra_cmnd, int max_reshower_attempts, bool shower_colorless) {
  Pythia pythia(PYTHIA_XML_DIR, false);
  if (!ConfigureIsolatedShowerPythia(pythia, extra_cmnd, seed)) { return 0; }

  HepMC3::WriterAscii writer(hepmc_file);
  if (writer.failed()) { throw std::runtime_error("Cannot open HepMC output " + hepmc_file); }
  int accepted = 0;
  for (const auto &block : input.events) {
    if (accepted >= max_events) { break; }
    if (!IsIsolatedShowerEvent(block, shower_colorless) && IsDirectColorSingletEvent(block) &&
        HasFullBeamRows(block.record)) {
      HepMC3::GenEvent output(HepMC3::Units::GEV, HepMC3::Units::MM);
      if (!BuildDirectColorSingletEvent(input, block, accepted, &pythia, output)) {
        throw std::runtime_error("Invalid direct event " + std::to_string(accepted + 1));
      }
      writer.write_event(output);
      if (writer.failed()) { throw std::runtime_error("HepMC write failed"); }
      ++accepted;
      continue;
    }
    if (!IsIsolatedShowerEvent(block, shower_colorless)) { throw std::runtime_error("Invalid isolated color state"); }
    HepMC3::GenEvent output(HepMC3::Units::GEV, HepMC3::Units::MM);
    ShowerIsolatedEventOrThrow(pythia, input, block, output, accepted + 1, seed, max_reshower_attempts,
                               shower_colorless);
    writer.write_event(output);
    if (writer.failed()) { throw std::runtime_error("HepMC write failed"); }
    ++accepted;
  }

  pythia.stat();
  return accepted;
}

// Run exact direct HepMC writing for color-singlet diffractive events
int RunDirectColorSinglet(const LHEInput &input, const std::string &hepmc_file, int max_events, int seed,
                          const std::string &extra_cmnd) {
  std::unique_ptr<Pythia> nuclear_pythia;
  if (HasNuclearBreakup(input, max_events)) {
    nuclear_pythia = std::make_unique<Pythia>(PYTHIA_XML_DIR, false);
    if (!ConfigureIsolatedShowerPythia(*nuclear_pythia, extra_cmnd, seed)) {
      throw std::runtime_error("Pythia nuclear fragmentation initialization failed");
    }
  }
  HepMC3::WriterAscii writer(hepmc_file);
  if (writer.failed()) { throw std::runtime_error("Cannot open HepMC output " + hepmc_file); }
  int accepted = 0;
  for (std::size_t i = 0; i < input.events.size() && accepted < max_events; ++i) {
    const LHEEventBlock &event = input.events[i];
    HepMC3::GenEvent     output(HepMC3::Units::GEV, HepMC3::Units::MM);
    if (!BuildDirectColorSingletEvent(input, event, i, nuclear_pythia.get(), output)) {
      if (NuclearBreakup(event.diff.side1) || NuclearBreakup(event.diff.side2)) {
        throw std::runtime_error("Pythia rejected nuclear fragmentation for event " + std::to_string(i + 1));
      }
      throw std::runtime_error("Invalid direct event " + std::to_string(i + 1));
    }
    writer.write_event(output);
    if (writer.failed()) { throw std::runtime_error("HepMC write failed"); }
    ++accepted;
  }
  return accepted;
}

// Run the standard streaming LHE path for non-diffractive files
int RunStandardLHE(const LHEInput &input, const std::string &lhe_file, const std::string &hepmc_file, int max_events, int seed,
                   const std::string &extra_cmnd) {
  Pythia pythia(PYTHIA_XML_DIR, false);
  if (!ConfigureFilePythia(pythia, lhe_file, extra_cmnd, seed)) { return 0; }

  Pythia8ToHepMC to_hepmc(hepmc_file);
  if (to_hepmc.output().failed()) { throw std::runtime_error("Cannot open HepMC output " + hepmc_file); }

  int accepted = 0;
  while (accepted < max_events) {
    if (!pythia.next()) {
      if (pythia.info.atEndOfFile()) { break; }
      throw std::runtime_error("Pythia rejected a source LHE event, refusing to drop its weight");
    }

    to_hepmc.setWeightNames(pythia.info.weightNameVector());
    if (!to_hepmc.fillNextEvent(pythia)) { throw std::runtime_error("Pythia HepMC conversion failed"); }
    // Preserve source normalization instead of Pythia's running estimate
    AddCrossSectionInfo(to_hepmc.event(), input.init);
    AddDoubleAttribute(to_hepmc.event(), "graniitti_lhe_weight", input.events.at(accepted).record.weight);
    AddGeneratorScaleInfo(to_hepmc.event(), input.events.at(accepted).record);
    to_hepmc.writeEvent();
    if (to_hepmc.output().failed()) { throw std::runtime_error("HepMC write failed"); }
    ++accepted;
  }

  pythia.stat();
  return accepted;
}

// Attach the event-local Pythia beam configuration used for remnant fragmentation
void AddPythiaBeamAttributes(HepMC3::GenEvent &event, const LHEInitInfo &init, const LHEEventBlock &block) {
  AddIntAttribute(event, "graniitti_pythia_beam1_id", PythiaBeamId(init, block.diff.side1, true));
  AddIntAttribute(event, "graniitti_pythia_beam2_id", PythiaBeamId(init, block.diff.side2, false));
  AddDoubleAttribute(event, "graniitti_pythia_beam1_e", PythiaBeamEnergy(init, block.diff, block.diff.side1, true));
  AddDoubleAttribute(event, "graniitti_pythia_beam2_e", PythiaBeamEnergy(init, block.diff, block.diff.side2, false));
  AddDoubleAttribute(event, "graniitti_pythia_x1", PythiaHardFraction(init, block, true));
  AddDoubleAttribute(event, "graniitti_pythia_x2", PythiaHardFraction(init, block, false));
  AddDoubleAttribute(event, "graniitti_pythia_hard1_e", block.record.particles[0].e);
  AddDoubleAttribute(event, "graniitti_pythia_hard2_e", block.record.particles[1].e);
}

// Convert one accepted event-local Pomeron shower and attach GRANIITTI metadata
bool ConvertDiffractiveEvent(Pythia &pythia, Pythia8ToHepMC &to_hepmc, const LHEInitInfo &init,
                                   const LHEEventBlock &event, std::size_t event_number, int attempt, int max_attempts,
                                   int pythia_seed, double minimum_remnant_energy) {
  to_hepmc.setWeightNames({"Weight"});
  if (!to_hepmc.fillNextEvent(pythia)) { return false; }
  AddEventWeight(to_hepmc.event(), event.record);
  AddCrossSectionInfo(to_hepmc.event(), init, event.diff.accepted_events, event.diff.attempted_events);
  AddReshowerAttributes(to_hepmc.event(), attempt, max_attempts, pythia_seed);
  AddStringAttribute(to_hepmc.event(), "graniitti_converter_mode", "full-shower");
  AddIntAttribute(to_hepmc.event(), "graniitti_low_remnant_isr_disabled", 0);
  AddDoubleAttribute(to_hepmc.event(), "graniitti_min_remnant_e", minimum_remnant_energy);
  AddPythiaBeamAttributes(to_hepmc.event(), init, event);
  if (!TransformPythiaSubeventPreservingHardBranches(to_hepmc.event(), init, event.diff, event.record, event_number)) {
    return false;
  }
  AnnotateDiffractiveEvent(to_hepmc.event(), init, event.diff);
  if (!CheckAndAnnotateClosure(to_hepmc.event(), init, event.diff, event_number)) { return false; }
  return true;
}

// Accumulate the weighted merging acceptance and its correlated numerator and denominator
struct MergingWeights {
  double sum_w = 0.0, sum_w2 = 0.0, sum_m = 0.0, sum_m2 = 0.0, sum_wm = 0.0;
  std::size_t n = 0;

  // Include every source event, including those that later receive a zero merging weight
  explicit MergingWeights(const LHEInput &input) : n(input.events.size()) {
    for (const auto &event : input.events) {
      sum_w += event.record.weight;
      sum_w2 += event.record.weight * event.record.weight;
    }
  }

  // Propagate the weighted acceptance variance and the independent source integration error
  void Apply(HepMC3::GenEvent &event, double factor) {
    if (!std::isfinite(factor) || factor < 0.0) { throw std::runtime_error("Invalid CKKW-L weight"); }
    const double source = event.weights()[0];
    event.weights()[0] *= factor;
    const double merged = event.weights()[0];
    sum_m += merged;
    sum_m2 += merged * merged;
    sum_wm += source * merged;
    const double ratio = sum_m / sum_w;
    const double residual = std::max(0.0, sum_m2 - 2.0 * ratio * sum_wm + ratio * ratio * sum_w2);
    const double variance = n > 1 ? static_cast<double>(n) / (n - 1) * residual / (sum_w * sum_w) : 0.0;
    auto xs = event.cross_section();
    AddDoubleAttribute(event, "graniitti_merging_weight", factor);
    AddDoubleAttribute(event, "graniitti_merging_source_sum", sum_w);
    AddDoubleAttribute(event, "graniitti_merging_source_xsec", xs->xsec());
    xs->set_cross_section(xs->xsec() * ratio, std::hypot(xs->xsec_err() * ratio, xs->xsec() * std::sqrt(variance)),
                         xs->get_accepted_events(), xs->get_attempted_events());
  }
};

// Shower one hard event without ever dropping it after a technical rejection
bool ShowerDiffractiveEventOrThrow(Pythia &pythia, const LHEInput &input, const LHEEventBlock &event,
                                   Pythia8ToHepMC &to_hepmc, std::size_t event_index, int max_reshower_attempts,
                                   PythiaSeedStream &seed_stream, const std::string &extra_cmnd) {
  auto         lha_up                 = std::make_shared<GraniittiEventLHAup>(input.init, event);
  const double minimum_remnant_energy = MinimumDPDFRemnantEnergy(input.init, event);
  int          pythia_seed            = seed_stream.NextSeed();
  if (!ConfigureExternalPythia(pythia, lha_up, extra_cmnd, pythia_seed)) {
    throw std::runtime_error("Pythia initialization failed for diffractive event " + std::to_string(event_index + 1));
  }
  // Query Pythia's own remnant support in the subcollision rest frame
  auto beam1 = pythia.beamA;
  auto beam2 = pythia.beamB;
  beam1.clear();
  beam2.clear();
  if (PythiaHardFraction(input.init, event, true) >= beam1.xMax(-1) ||
      PythiaHardFraction(input.init, event, false) >= beam2.xMax(-1)) {
    return false;
  }
  int showered = 0;
  for (int attempt = 1; attempt <= max_reshower_attempts; ++attempt) {
    if (attempt > 1) { pythia_seed = seed_stream.NextSeed(); }
    pythia.rndm.init(pythia_seed);
    lha_up->Rewind();
    if (!pythia.next()) {
      if (pythia.settings.flag("Merging:doKTMerging")) {
        throw std::runtime_error("CKKW-L shower failed, refusing to redraw a merging weight");
      }
      continue;
    }
    ++showered;
    if (ConvertDiffractiveEvent(pythia, to_hepmc, input.init, event, event_index + 1, attempt,
                                      max_reshower_attempts, pythia_seed, minimum_remnant_energy)) {
      return true;
    }
    if (pythia.settings.flag("Merging:doKTMerging")) {
      throw std::runtime_error("CKKW-L event conversion failed, refusing to redraw a merging weight");
    }
  }

  pythia.stat();
  throw std::runtime_error("Could not shower and convert diffractive LHAup event " + std::to_string(event_index + 1) +
                           " after " + std::to_string(max_reshower_attempts) +
                           " attempts (" + std::to_string(showered) + " showers succeeded), refusing to drop the event");
}

// Run diffractive events with exact event-local Pomeron or proton beams
int RunDiffractiveLHAup(const LHEInput &input, const std::string &hepmc_file, int max_events, int seed,
                        const std::string &extra_cmnd, int max_reshower_attempts) {
  Pythia         pythia(PYTHIA_XML_DIR, false);
  // Keep the hook alive across event-specific reinitialization
  pythia.setUserHooksPtr(std::make_shared<GraniittiDiffractionHooks>());
  Pythia8ToHepMC to_hepmc(hepmc_file);
  if (to_hepmc.output().failed()) { throw std::runtime_error("Cannot open HepMC output " + hepmc_file); }
  Pythia isolated(PYTHIA_XML_DIR, false);
  if (!ConfigureIsolatedShowerPythia(isolated, extra_cmnd, seed)) {
    throw std::runtime_error("Pythia isolated hadronization initialization failed");
  }
  PythiaSeedStream seed_stream(seed);
  int              accepted = 0;
  MergingWeights weights(input);

  for (std::size_t i = 0; i < input.events.size() && accepted < max_events; ++i) {
    const LHEEventBlock &event = input.events[i];
    if (!event.diff.active) { continue; }

    LHEEventBlock hard = event;
    PrepareHardShower(hard);
    if (!ShowerDiffractiveEventOrThrow(pythia, input, hard, to_hepmc, i, max_reshower_attempts,
                                      seed_stream, extra_cmnd)) {
      if (!HasFullBeamRows(event.record)) {
        throw std::runtime_error("Pythia remnant support requires a complete source record for fragmentation");
      }
      std::cout << "PYTHIA_LHE_HADRONIZE remnant-support fragmentation event=" << i + 1 << std::endl;
      HepMC3::GenEvent output(HepMC3::Units::GEV, HepMC3::Units::MM);
      ShowerIsolatedEventOrThrow(isolated, input, event, output, i + 1, seed,
                                 max_reshower_attempts, false);
      AddIntAttribute(output, "graniitti_pythia_remnant_fallback", 1);
      AddIntAttribute(output, "graniitti_low_remnant_isr_disabled", pythia.settings.flag("PartonLevel:ISR") ? 1 : 0);
      if (pythia.settings.flag("Merging:doKTMerging")) {
        AddIntAttribute(output, "graniitti_merging_remnant_veto", 1);
        weights.Apply(output, 0.0);
      }
      to_hepmc.output().write_event(output);
      if (to_hepmc.output().failed()) { throw std::runtime_error("HepMC write failed"); }
    } else {
      if (pythia.settings.flag("Merging:doKTMerging")) {
        weights.Apply(to_hepmc.event(), pythia.info.mergingWeight());
      }
      to_hepmc.writeEvent();
      if (to_hepmc.output().failed()) { throw std::runtime_error("HepMC write failed"); }
    }
    ++accepted;
  }
  return accepted;
}

}  // namespace

// Main program for LHE input and HepMC3 output
int main(int argc, char *argv[]) try {
  if (argc < 4 || argc > 8) {
    PrintUsage(argv[0]);
    return EXIT_FAILURE;
  }

  const std::string lhe_file              = argv[1];
  const std::string hepmc_file            = argv[2];
  const int         max_events            = std::string(argv[3]) == "all" ? std::numeric_limits<int>::max()
                                                                       : ReadPositiveInt(argv[3], "nevents");
  const int         seed                  = (argc >= 5) ? ReadNonNegativeInt(argv[4], "seed") : 0;
  const std::string extra_cmnd            = (argc >= 6) ? argv[5] : "";
  const int         max_reshower_attempts = ReadReshowerAttempts(argc, argv);
  const std::string mode                  = argc >= 8 ? argv[7] : "auto";
  if (mode != "auto" && mode != "fragment" && mode != "shower") {
    throw std::invalid_argument("Converter mode must be auto, fragment or shower");
  }
  std::error_code error;
  if (std::filesystem::equivalent(lhe_file, hepmc_file, error)) {
    throw std::invalid_argument("LHE input and HepMC output refer to the same file");
  }

  if (!extra_cmnd.empty() && std::filesystem::equivalent(extra_cmnd, hepmc_file, error)) {
    throw std::invalid_argument("Pythia steering input and HepMC output refer to the same file");
  }
  Pythia steering(PYTHIA_XML_DIR, false);
  // Isolated color singlet decays radiate photons only when explicitly requested
  steering.readString("TimeShower:QEDshowerByL = off");
  if (!extra_cmnd.empty() && !steering.readFile(extra_cmnd)) {
    throw std::invalid_argument("Cannot read Pythia settings from " + extra_cmnd);
  }
  for (const auto *setting : {"Merging:doMGMerging", "Merging:doMerging", "Merging:doPTLundMerging",
                              "Merging:doCutBasedMerging", "Merging:doDynamicMerging", "Merging:doUserMerging",
                              "Merging:doXSectionEstimate", "Merging:doUMEPSTree", "Merging:doUMEPSSubt",
                              "Merging:doNL3Tree", "Merging:doNL3Loop", "Merging:doNL3Subt", "Merging:doUNLOPSTree",
                              "Merging:doUNLOPSLoop", "Merging:doUNLOPSSubt", "Merging:doUNLOPSSubtNLO"}) {
    if (steering.settings.flag(setting)) {
      throw std::invalid_argument(std::string("Unsupported merging prescription ") + setting + ", use Merging:doKTMerging");
    }
  }
  LHEInput   input                = ReadLHEInput(lhe_file, max_events);
  const bool isr_requested        = steering.settings.flag("PartonLevel:ISR");
  const bool lepton_fsr_requested = steering.settings.flag("PartonLevel:FSR") &&
                                    steering.settings.flag("TimeShower:QEDshowerByL");
  const bool hard_shower = mode == "shower" ||
      (mode == "auto" && !input.events.empty() &&
       std::all_of(input.events.begin(), input.events.end(), [](const LHEEventBlock &event) {
         return event.diff.active && event.hard.valid;
       }));
  if (steering.settings.flag("Merging:doKTMerging") &&
      (!hard_shower || std::any_of(input.events.begin(), input.events.end(),
                                   [](const LHEEventBlock &event) { return !event.diff.active; }))) {
    throw std::invalid_argument("CKKW-L conversion requires the diffractive full shower mode");
  }
  if (steering.settings.flag("Merging:doKTMerging")) {
    const MergingWeights weights(input);
    if (!(weights.sum_w > 0.0) || !std::isfinite(weights.sum_w2) || input.init.process_lines.size() != 1) {
      throw std::invalid_argument("CKKW-L conversion requires one process and a finite positive source weight sum");
    }
  }
  int        accepted             = 0;
  if (HasDiffractiveEvents(input) || HasFullBeamRecords(input)) {
    for (const auto &event : input.events) {
      if (!event.diff.active && !HasFullBeamRows(event.record)) {
        throw std::invalid_argument("Reduced diffractive records require diffraction metadata");
      }
      if (HasFullBeamRows(event.record) && !IsIsolatedShowerEvent(event, false) && !IsDirectColorSingletEvent(event)) {
        throw std::invalid_argument("Incomplete or invalid remnant colors in the full beam record");
      }
    }
    ValidateDiffractiveCounters(input);
    const bool direct_color_singlet = CanUseDirectColorSingletMode(input, max_events);
    const bool direct_full_event    = CanUseDirectFullEventMode(input, max_events);
    // Complete pp records cannot acquire meaningful beam ISR through the reduced LHA path
    const bool direct_allowed = direct_color_singlet && (direct_full_event || !isr_requested);
    if (!hard_shower && CanUseIsolatedShowerMode(input, max_events, lepton_fsr_requested)) {
      std::cout << "PYTHIA_LHE_HADRONIZE mode=isolated-shower" << std::endl;
      accepted = RunIsolatedShower(input, hepmc_file, max_events, seed, extra_cmnd, max_reshower_attempts,
                                   lepton_fsr_requested);
    } else if (!hard_shower && direct_allowed) {
      std::cout << "PYTHIA_LHE_HADRONIZE mode=direct-copy" << std::endl;
      accepted = RunDirectColorSinglet(input, hepmc_file, max_events, seed, extra_cmnd);
    } else {
      if (mode == "fragment") {
        throw std::invalid_argument("Fragment mode requires complete color singlet beam and remnant states");
      }
      if (!hard_shower) {
        for (auto &event : input.events) {
          if (HasFullBeamRows(event.record)) {
            throw std::invalid_argument("Incomplete or invalid remnant colors in the full beam record");
          }
        }
      }
      if (std::any_of(input.events.begin(), input.events.end(),
                      [](const LHEEventBlock &event) { return !event.diff.active; })) {
        throw std::invalid_argument("Diffractive full showers require metadata for every source event");
      }
      std::cout << "PYTHIA_LHE_HADRONIZE mode=full-shower" << std::endl;
      accepted = RunDiffractiveLHAup(input, hepmc_file, max_events, seed, extra_cmnd, max_reshower_attempts);
    }
  } else {
    if (mode == "fragment") { throw std::invalid_argument("Fragment mode requires GRANIITTI diffraction metadata"); }
    if (std::abs(input.init.strategy) != 3 && std::abs(input.init.strategy) != 4) {
      throw std::invalid_argument("Streaming conversion requires LHE weight strategy +/-3 or +/-4");
    }
    std::cout << "PYTHIA_LHE_HADRONIZE mode=standard-shower" << std::endl;
    accepted = RunStandardLHE(input, lhe_file, hepmc_file, max_events, seed, extra_cmnd);
  }

  std::cout << "PYTHIA_LHE_HADRONIZE accepted " << accepted << " events from " << lhe_file << " into " << hepmc_file
            << std::endl;

  return accepted > 0 ? EXIT_SUCCESS : EXIT_FAILURE;
} catch (const std::exception &error) {
  std::cerr << "pythia_lhe_hadronize: " << error.what() << '\n';
  return EXIT_FAILURE;
}
