// YFS (Yennie-Frautschi-Suura) type QED initial and final state photon
// radiation
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

// HepMC3
#include "HepMC3/Attribute.h"
#include "HepMC3/FourVector.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"

// Libraries
#include "json.hpp"

// Own
#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Particle/MPDG.h"
#include "Graniitti/Photon/MQED.h"
#include "Graniitti/Photon/MRadiative.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MException.h"

using gra::aux::indices;

namespace gra {
namespace radiative {

namespace {

constexpr double kEuler = 0.57721566490153286061;
constexpr double kRecoilTol = 1.0e-12;

// Compute true for one charged lepton
bool IsChargedLepton(const MParticle &particle) {
  const int pdg = std::abs(particle.pdg);
  return pdg == 11 || pdg == 13 || pdg == 15;
}

// Compute true for one charged lepton
bool IsChargedLepton(int pdg) {
  pdg = std::abs(pdg);
  return pdg == 11 || pdg == 13 || pdg == 15;
}

// Compute true for one opposite-sign same-flavour charged lepton pair
bool IsLeptonPair(int pdg1, int pdg2) { return IsChargedLepton(pdg1) && pdg1 == -pdg2; }

// Compute true when a decay tree contains a directly radiating lepton pair
bool HasLeptonPair(const std::vector<MDecayBranch> &tree) {
  if (tree.size() == 2 && tree[0].legs.empty() && tree[1].legs.empty() && IsLeptonPair(tree[0].p.pdg, tree[1].p.pdg)) {
    return true;
  }
  return std::any_of(tree.begin(), tree.end(), [](const MDecayBranch &branch) { return HasLeptonPair(branch.legs); });
}

// Parse one explicit radiative mode
Mode ParseMode(const nlohmann::json &value, const std::string &name) {
  if (!value.is_string()) {
    throw std::invalid_argument("RADIATIVE." + name +
                                " must be 'none' or 'YFS'");
  }
  const std::string mode = value.get<std::string>();
  if (mode == "none") {
    return Mode::Off;
  }
  if (mode == "YFS") {
    return Mode::YFS;
  }
  throw std::invalid_argument("RADIATIVE." + name + " must be 'none' or 'YFS'");
}

// Reject misspelled or misplaced radiative steering fields
void ValidateModeKeys(const nlohmann::json &modes) {
  for (const auto &[name, value] : modes.items()) {
    (void)value;
    if (name != "ISR_QED" && name != "FSR_QED") {
      throw std::invalid_argument("RADIATIVE has unknown field " + name);
    }
  }
}

// Validate all numerical radiative controls
void Validate(const Config &config) {
  if (!std::isfinite(config.isr_param.x_min) ||
      !(config.isr_param.x_min > 0.0) || !(config.isr_param.x_min < 1.0)) {
    throw std::invalid_argument(
        "NUMERICS_RADIATIVE.ISR_QED.x_min must be between zero and one");
  }
  if (!std::isfinite(config.fsr_param.energy_min) ||
      !(config.fsr_param.energy_min > 0.0)) {
    throw std::invalid_argument(
        "NUMERICS_RADIATIVE.FSR_QED.energy_min must be positive");
  }
  if (!std::isfinite(config.fsr_param.energy_max_fraction) ||
      !(config.fsr_param.energy_max_fraction > 0.0) ||
      config.fsr_param.energy_max_fraction > 0.5) {
    throw std::invalid_argument(
        "NUMERICS_RADIATIVE.FSR_QED.energy_max_fraction must be in (0,0.5]");
  }
  if (config.fsr_param.max_photons == 0 || config.fsr_param.max_trials == 0) {
    throw std::invalid_argument(
        "NUMERICS_RADIATIVE FSR counts must be positive");
  }
}

// Read a positive size parameter without accepting negative JSON integers
std::size_t ReadSize(const nlohmann::json &value, const std::string &name) {
  if (!value.is_number_integer()) {
    throw std::invalid_argument(name + " must be a positive integer");
  }
  const auto input = value.get<long long>();
  if (input <= 0) {
    throw std::invalid_argument(name + " must be a positive integer");
  }
  return static_cast<std::size_t>(input);
}

// Construct an exact two-body recoil with a charge symmetric relative axis
bool TwoBody(const M4Vec &system, const M4Vec &source1, const M4Vec &source2, double mass1,
             double mass2, M4Vec &momentum1, M4Vec &momentum2) {
  const double mass = system.M();
  if (!std::isfinite(system.E()) || !(system.E() > 0.0) ||
      !std::isfinite(mass) || !(mass > mass1 + mass2)) {
    return false;
  }
  const double pnorm = kinematics::DecayMomentum(mass, mass1, mass2);
  if (!std::isfinite(pnorm) || !(pnorm > 0.0)) {
    return false;
  }

  std::array<M4Vec, 2> axis = {source1, source2};
  for (auto &momentum : axis) {
    kinematics::LorentzBoost(system, mass, momentum, -1);
    if (!(momentum.E() > 0.0) || !gra::AllFinite(momentum.Contravariant())) { return false; }
  }
  const M3Vec relative = (axis[0] - axis[1]).P3();
  if (!(gra::SquaredNorm(relative) > 0.0)) { return false; }
  const M3Vec direction =
      kinematics::NormalizeFrameAxis(relative, {0.0, 0.0, 1.0});
  const double energy1 =
      (mass * mass + mass1 * mass1 - mass2 * mass2) / (2.0 * mass);
  const double energy2 =
      (mass * mass + mass2 * mass2 - mass1 * mass1) / (2.0 * mass);

  momentum1 = M4Vec(pnorm * direction[0], pnorm * direction[1],
                    pnorm * direction[2], energy1);
  momentum2 = M4Vec(-pnorm * direction[0], -pnorm * direction[1],
                    -pnorm * direction[2], energy2);
  kinematics::LorentzBoost(system, mass, momentum1, 1);
  kinematics::LorentzBoost(system, mass, momentum2, 1);
  return momentum1.E() > 0.0 && momentum2.E() > 0.0 &&
         gra::AllFinite(momentum1.Contravariant()) && gra::AllFinite(momentum2.Contravariant());
}

// Sample one exponentiated electron structure-function coordinate
// 1-x = (1-x_min) u^(1/beta)
double StructurePoint(double unit, double exponent, const ISRParam &param,
                      double &weight) {
  if (!(exponent > 0.0)) {
    weight = 1.0;
    return 1.0;
  }
  const double ymax = 1.0 - param.x_min;
  const double loss = ymax * std::pow(std::max(unit, 0.0), 1.0 / exponent);
  const double x = std::clamp(1.0 - loss, param.x_min, 1.0);
  // [REFERENCE: Greco et al., arXiv:1711.00826, Eqs. (2)-(4), through second order]
  const double soft =
      std::exp(exponent * (0.75 - kEuler)) / std::tgamma(1.0 + exponent);
  weight = std::pow(ymax, exponent) * soft;
  if (!param.hard_collinear || !(loss > 0.0)) { return x; }

  // Keep the sampled loss to resolve the soft limit even when x rounds to one
  const double logx = loss < 0.5 ? std::log1p(-loss) : std::log(x);
  const double second = (1.0 + x) * (-4.0 * std::log(loss) + 3.0 * logx) -
                        4.0 * logx / loss - 5.0 - x;
  const double finite = -0.5 * (1.0 + x) + exponent * second / 8.0;
  weight += std::pow(ymax, exponent) * std::pow(loss, 1.0 - exponent) * finite;
  return x;
}

// Sample one direction from the exact massive final-state dipole
bool SampleDipoleDirection(double beta, const M3Vec &axis, MRandom &random,
                           std::size_t max_trials, M3Vec &direction) {
  const double bounded_beta = std::clamp(beta, 0.0, 1.0 - 1.0e-14);
  if (!(bounded_beta > 0.0)) {
    return false;
  }
  const double rapidity = std::atanh(bounded_beta);
  const M3Vec zaxis = kinematics::NormalizeFrameAxis(axis, {0.0, 0.0, 1.0});
  const M3Vec yaxis = kinematics::PerpendicularFrameAxis(zaxis);
  const M3Vec xaxis = kinematics::NormalizeFrameAxis(
      gra::CrossProduct(yaxis, zaxis), {1.0, 0.0, 0.0});

  for (std::size_t trial = 0; trial < max_trials; ++trial) {
    const double unit = random.U(0.0, 1.0);
    const double cosine =
        std::tanh((2.0 * unit - 1.0) * rapidity) / bounded_beta;
    const double denominator =
        1.0 - bounded_beta * bounded_beta * cosine * cosine;
    // W(c) <= 4 beta^2/(1-beta^2 c^2), with finite threshold acceptance
    const double accept = std::clamp((1.0 - cosine * cosine) / denominator, 0.0, 1.0);
    if (random.U(0.0, 1.0) > accept) {
      continue;
    }
    const double phi = 2.0 * math::PI * random.U(0.0, 1.0);
    const double sine = math::msqrt(std::max(0.0, 1.0 - cosine * cosine));
    direction = {sine * std::cos(phi) * xaxis[0] +
                     sine * std::sin(phi) * yaxis[0] + cosine * zaxis[0],
                 sine * std::cos(phi) * xaxis[1] +
                     sine * std::sin(phi) * yaxis[1] + cosine * zaxis[1],
                 sine * std::cos(phi) * xaxis[2] +
                     sine * std::sin(phi) * yaxis[2] + cosine * zaxis[2]};
    return true;
  }
  return false;
}

// Sample one exclusive soft-photon lepton-pair state
bool SamplePair(const M4Vec &lepton1, const M4Vec &lepton2,
                const FSRParam &param, MRandom &random, FSRPair &state) {
  const M4Vec parent = lepton1 + lepton2;
  const double mass = parent.M();
  const double lepton_mass = std::abs(lepton1.M());
  if (!std::isfinite(mass) || !std::isfinite(lepton_mass) ||
      !(mass > 2.0 * lepton_mass)) {
    return false;
  }

  M4Vec lepton1_rest = lepton1;
  kinematics::LorentzBoost(parent, mass, lepton1_rest, -1);
  const double beta = lepton1_rest.Beta();
  const double maximum =
      std::min(param.energy_max_fraction * mass,
               (mass * mass - 4.0 * lepton_mass * lepton_mass) / (2.0 * mass));
  if (!(maximum > param.energy_min) || !(beta > 0.0)) {
    return false;
  }

  // [REFERENCE: Schonherr and Krauss, JHEP 12 (2008) 018]
  const double mean = qed::alpha_QED() / math::PI * DipoleIntegral(beta) *
                      std::log(maximum / param.energy_min);
  if (!std::isfinite(mean) || !(mean > 0.0)) {
    return false;
  }
  for (std::size_t trial = 0; trial < param.max_trials; ++trial) {
    // Repeat the complete Poisson proposal after a recoil veto
    // [REFERENCE: Schonherr and Krauss, arXiv:0810.5071, Section 3.4, steps 1-3]
    const int photon_count = random.PoissonRandom(mean);
    if (photon_count <= 0) { return false; }
    if (static_cast<std::size_t>(photon_count) > param.max_photons) {
      throw PhaseSpaceFailure("MYFS: FSR photon multiplicity exceeds max_photons");
    }
    std::vector<M4Vec> photons;
    photons.reserve(static_cast<std::size_t>(photon_count));
    bool directions_ok = true;
    for (int photon = 0; photon < photon_count; ++photon) {
      const double energy =
          param.energy_min *
          std::pow(maximum / param.energy_min, random.U(0.0, 1.0));
      M3Vec direction{};
      if (!SampleDipoleDirection(beta, lepton1_rest.P3(), random,
                                 param.max_trials, direction)) {
        directions_ok = false;
        break;
      }
      M4Vec momentum(energy * direction[0], energy * direction[1],
                     energy * direction[2], energy);
      kinematics::LorentzBoost(parent, mass, momentum, 1);
      photons.push_back(momentum);
    }
    if (!directions_ok) {
      throw PhaseSpaceFailure("MYFS: FSR angular sampling exhausted max_trials");
    }

    if (!MapLeptonPair(lepton1, lepton2, photons, state.lepton[0],
                       state.lepton[1])) {
      continue;
    }
    state.born = {lepton1, lepton2};
    state.photons = std::move(photons);
    return true;
  }
  throw PhaseSpaceFailure("MYFS: FSR recoil sampling exhausted max_trials");
}

// Commit one sampled recoil to a HepMC lepton-pair vertex
void CommitPair(const HepMC3::GenVertexPtr &vertex, const FSRPair &state,
                bool swapped = false) {
  const auto &outgoing = vertex->particles_out();
  outgoing[0]->set_momentum(
      aux::M4Vec2HepMC3(state.lepton[swapped ? 1 : 0]));
  outgoing[1]->set_momentum(
      aux::M4Vec2HepMC3(state.lepton[swapped ? 0 : 1]));
  for (const auto &momentum : state.photons) {
    auto photon = std::make_shared<HepMC3::GenParticle>(
        aux::M4Vec2HepMC3(momentum), PDG::PDG_gamma, PDG::PDG_STABLE);
    vertex->add_particle_out(photon);
    photon->add_attribute("QED_FSR",
                          std::make_shared<HepMC3::IntAttribute>(1));
  }
}

// Compute true when one HepMC vertex is an unradiated stable lepton pair
bool LeptonVertex(const HepMC3::GenVertexPtr &vertex) {
  if (vertex == nullptr || vertex->particles_in().size() != 1 ||
      vertex->particles_out().size() != 2) {
    return false;
  }
  const auto &particle1 = vertex->particles_out()[0];
  const auto &particle2 = vertex->particles_out()[1];
  return particle1 != nullptr && particle2 != nullptr &&
         particle1->status() == PDG::PDG_STABLE &&
         particle2->status() == PDG::PDG_STABLE &&
         particle1->end_vertex() == nullptr &&
         particle2->end_vertex() == nullptr &&
         IsLeptonPair(particle1->pid(), particle2->pid());
}

// Generate and commit one exclusive soft-photon HepMC vertex
void RadiateVertex(const HepMC3::GenVertexPtr &vertex, const FSRParam &param,
                   MRandom &random) {
  if (!LeptonVertex(vertex)) {
    return;
  }
  const auto &outgoing = vertex->particles_out();
  FSRPair state;
  state.pdg = {outgoing[0]->pid(), outgoing[1]->pid()};
  if (SamplePair(aux::HepMC2M4Vec(outgoing[0]->momentum()),
                 aux::HepMC2M4Vec(outgoing[1]->momentum()), param, random,
                 state)) {
    CommitPair(vertex, state);
  }
}

// Sample one stable decay-tree lepton pair and cache its recoil
void PreparePair(std::vector<MDecayBranch> &tree, const FSRParam &param,
                 MRandom &random, std::vector<FSRPair> &pairs) {
  if (tree.size() != 2 || !tree[0].legs.empty() || !tree[1].legs.empty() ||
      !IsLeptonPair(tree[0].p.pdg, tree[1].p.pdg)) {
    return;
  }
  FSRPair state;
  state.pdg = {tree[0].p.pdg, tree[1].p.pdg};
  if (!SamplePair(tree[0].p4, tree[1].p4, param, random, state)) {
    return;
  }
  tree[0].p4 = state.lepton[0];
  tree[1].p4 = state.lepton[1];
  pairs.push_back(std::move(state));
}

// Sample all supported vertices in one recursive decay tree
void PrepareTree(std::vector<MDecayBranch> &tree, const FSRParam &param,
                 MRandom &random, std::vector<FSRPair> &pairs) {
  PreparePair(tree, param, random, pairs);
  for (auto &branch : tree) {
    PrepareTree(branch.legs, param, random, pairs);
  }
}

// Compute the largest absolute difference between two four-momenta
double MomentumDifference(const M4Vec &first, const M4Vec &second) {
  return std::max({std::abs(first.Px() - second.Px()),
                   std::abs(first.Py() - second.Py()),
                   std::abs(first.Pz() - second.Pz()),
                   std::abs(first.E() - second.E())});
}

// Compute true when two event-record momenta describe the same Born leg
bool SameMomentum(const M4Vec &first, const M4Vec &second) {
  const double scale = std::max(
      {std::abs(first.Px()), std::abs(first.Py()), std::abs(first.Pz()),
       std::abs(first.E()), std::abs(second.Px()), std::abs(second.Py()),
       std::abs(second.Pz()), std::abs(second.E()), 1.0});
  return MomentumDifference(first, second) <= 1.0e-11 * scale;
}

// Match one cached recoil to an unradiated HepMC lepton vertex
int MatchPair(const HepMC3::GenVertexPtr &vertex, const FSRPair &state) {
  if (!LeptonVertex(vertex)) {
    return -1;
  }
  const auto &outgoing = vertex->particles_out();
  const std::array<M4Vec, 2> momentum = {
      aux::HepMC2M4Vec(outgoing[0]->momentum()),
      aux::HepMC2M4Vec(outgoing[1]->momentum())};
  if (outgoing[0]->pid() == state.pdg[0] &&
      outgoing[1]->pid() == state.pdg[1] &&
      SameMomentum(momentum[0], state.born[0]) &&
      SameMomentum(momentum[1], state.born[1])) {
    return 0;
  }
  if (outgoing[0]->pid() == state.pdg[1] &&
      outgoing[1]->pid() == state.pdg[0] &&
      SameMomentum(momentum[0], state.born[1]) &&
      SameMomentum(momentum[1], state.born[0])) {
    return 1;
  }
  return -1;
}

} // namespace

// Parse and validate the gencard modes and tune numerical controls
Config ReadConfig(const nlohmann::json &modes, const nlohmann::json &numerics) {
  if (!modes.is_object() || !numerics.is_object()) {
    throw std::invalid_argument(
        "RADIATIVE and NUMERICS_RADIATIVE must be JSON objects");
  }
  ValidateModeKeys(modes);
  Config config;
  config.isr = ParseMode(modes.value("ISR_QED", "none"), "ISR_QED");
  config.fsr = ParseMode(modes.value("FSR_QED", "none"), "FSR_QED");

  const auto &isr = numerics.at("ISR_QED");
  config.isr_param.x_min = isr.at("x_min").get<double>();
  config.isr_param.hard_collinear = isr.at("hard_collinear").get<bool>();

  const auto &fsr = numerics.at("FSR_QED");
  config.fsr_param.energy_min = fsr.at("energy_min").get<double>();
  config.fsr_param.energy_max_fraction =
      fsr.at("energy_max_fraction").get<double>();
  config.fsr_param.max_photons =
      ReadSize(fsr.at("max_photons"),
               "NUMERICS_RADIATIVE.FSR_QED.max_photons");
  config.fsr_param.max_trials =
      ReadSize(fsr.at("max_trials"),
               "NUMERICS_RADIATIVE.FSR_QED.max_trials");
  Validate(config);
  return config;
}

// Compute the canonical steering name of one radiative mode
std::string ModeName(Mode mode) { return mode == Mode::YFS ? "YFS" : "none"; }

// Compute the exponent of the soft electron structure function
// beta = alpha/pi [ln(Q2/m2)-1]
double ISRExponent(double scale2, double mass) {
  if (!std::isfinite(scale2) || !std::isfinite(mass) || !(scale2 > 0.0) ||
      !(mass > 0.0)) {
    return 0.0;
  }
  const double logarithm = std::log(scale2 / (mass * mass));
  return std::max(0.0, qed::alpha_QED() / math::PI * (logarithm - 1.0));
}

// Compute the integrated massive final-state dipole factor
// I(beta) = (1+beta^2) ln[(1+beta)/(1-beta)]/beta-2
double DipoleIntegral(double beta) {
  if (!std::isfinite(beta) || !(beta > 0.0) || !(beta < 1.0)) {
    return 0.0;
  }
  if (beta < 1.0e-2) {
    const double b2 = beta * beta;
    return b2 * (8.0 / 3.0 + b2 * (16.0 / 15.0 + b2 * (24.0 / 35.0 + b2 * 32.0 / 63.0)));
  }
  return (1.0 + beta * beta) / beta * std::log1p(2.0 * beta / (1.0 - beta)) -
         2.0;
}

// Compute the massive final-state dipole angular kernel
// W(c) = 4 beta^2 (1-c^2)/(1-beta^2 c^2)^2 avoids threshold cancellation
double DipoleAngular(double beta, double cosine) {
  if (!std::isfinite(beta) || !std::isfinite(cosine) || !(beta > 0.0) ||
      !(beta < 1.0) || cosine < -1.0 || cosine > 1.0) {
    return 0.0;
  }
  const double minus = 1.0 - beta * cosine;
  const double plus = 1.0 + beta * cosine;
  return 4.0 * beta * beta * (1.0 - cosine) * (1.0 + cosine) / math::pow2(minus * plus);
}

// Map a lepton pair against explicit photons with exact four-momentum closure
bool MapLeptonPair(const M4Vec &lepton1, const M4Vec &lepton2,
                   const std::vector<M4Vec> &photons, M4Vec &mapped1,
                   M4Vec &mapped2) {
  // Reject spacelike or past directed leptons before extracting their masses
  for (const auto &lepton : {lepton1, lepton2}) {
    if (!gra::AllFinite(lepton.Contravariant()) || !(lepton.E() > 0.0) ||
        lepton.M2() < -kRecoilTol * lepton.E() * lepton.E()) {
      return false;
    }
  }
  const M4Vec parent = lepton1 + lepton2;
  M4Vec photon_sum;
  for (const auto &photon : photons) {
    if (!gra::AllFinite(photon.Contravariant()) || !(photon.E() > 0.0) ||
        std::abs(photon.M2()) > kRecoilTol * photon.E() * photon.E()) {
      return false;
    }
    photon_sum += photon;
  }
  const M4Vec recoil = parent - photon_sum;
  const double mass1 = std::sqrt(std::max(0.0, lepton1.M2()));
  const double mass2 = std::sqrt(std::max(0.0, lepton2.M2()));
  if (!std::isfinite(recoil.E()) || !(recoil.E() > 0.0) || !std::isfinite(recoil.M2()) ||
      recoil.M2() <= (mass1 + mass2) * (mass1 + mass2) * (1.0 + kRecoilTol)) {
    return false;
  }
  return TwoBody(recoil, lepton1, lepton2, mass1, mass2, mapped1, mapped2);
}

// Configure both radiative sectors
void MYFS::Configure(const Config &input) {
  Validate(input);
  config = input;
  ResetFSR();
}

// Validate that enabled radiation sectors have physical charged leptons
void MYFS::ValidateApplicability(const MParticle &beam1, const MParticle &beam2,
                                 const std::vector<MDecayBranch> &tree) const {
  if (config.isr == Mode::YFS && !IsChargedLepton(beam1) && !IsChargedLepton(beam2)) {
    throw std::invalid_argument("RADIATIVE.ISR_QED=YFS requires at least one charged-lepton beam");
  }
  if (config.fsr == Mode::YFS && !HasLeptonPair(tree)) {
    throw std::invalid_argument(
        "RADIATIVE.FSR_QED=YFS requires a direct opposite-sign same-flavour charged-lepton decay vertex");
  }
}

// Store the nominal incoming beam momenta
void MYFS::SetBeams(const M4Vec &beam1, const M4Vec &beam2) {
  isr_state.nominal = {beam1, beam2};
  isr_state.hard = isr_state.nominal;
  isr_state.photons.clear();
  isr_state.active = false;
  isr_state.weight = 1.0;
}

// Compute one sampling coordinate for each radiating lepton beam
unsigned int MYFS::ISRDim(const MParticle &beam1,
                          const MParticle &beam2) const {
  if (config.isr != Mode::YFS) {
    return 0;
  }
  return static_cast<unsigned int>(IsChargedLepton(beam1)) +
         static_cast<unsigned int>(IsChargedLepton(beam2));
}

// Restore the nominal incoming state before a new sampled point
void MYFS::ResetISR(M4Vec &beam1, M4Vec &beam2, double &s,
                    double &sqrt_s) noexcept {
  beam1 = isr_state.nominal[0];
  beam2 = isr_state.nominal[1];
  s = (beam1 + beam2).M2();
  sqrt_s = s > 0.0 ? std::sqrt(s) : 0.0;
  isr_state.hard = isr_state.nominal;
  isr_state.photons.clear();
  isr_state.active = false;
  isr_state.weight = 1.0;
}

// Generate the ISR convolution point and exact hard-beam recoil
double MYFS::GenerateISR(const std::vector<double> &random, std::size_t offset,
                         const MParticle &beam1, const MParticle &beam2,
                         M4Vec &momentum1, M4Vec &momentum2, double &s,
                         double &sqrt_s) {
  const unsigned int dimension = ISRDim(beam1, beam2);
  if (dimension == 0) {
    return 1.0;
  }
  if (random.size() < offset + dimension) {
    return 0.0;
  }

  const M4Vec parent = isr_state.nominal[0] + isr_state.nominal[1];
  const double parent_mass = parent.M();
  if (!std::isfinite(parent_mass) || !(parent_mass > beam1.mass + beam2.mass)) {
    return 0.0;
  }
  std::array<M4Vec, 2> nominal_rest = isr_state.nominal;
  for (auto &momentum : nominal_rest) {
    kinematics::LorentzBoost(parent, parent_mass, momentum, -1);
  }

  std::array<double, 2> fraction = {1.0, 1.0};
  double radiator_weight = 1.0;
  std::size_t coordinate = offset;
  const std::array<MParticle, 2> particles = {beam1, beam2};
  for (const auto &leg : indices(particles)) {
    if (!IsChargedLepton(particles[leg])) {
      continue;
    }
    double leg_weight = 1.0;
    fraction[leg] = StructurePoint(
        random[coordinate++], ISRExponent(parent.M2(), particles[leg].mass),
        config.isr_param, leg_weight);
    radiator_weight *= leg_weight;
  }
  if (!std::isfinite(radiator_weight) || !(radiator_weight > 0.0)) {
    return 0.0;
  }

  std::vector<M4Vec> photons_rest;
  for (const auto &leg : indices(particles)) {
    if (!IsChargedLepton(particles[leg])) {
      continue;
    }
    const double energy = (1.0 - fraction[leg]) * nominal_rest[leg].E();
    if (!(energy > std::numeric_limits<double>::epsilon() * parent_mass)) {
      continue;
    }
    const M3Vec direction = kinematics::NormalizeFrameAxis(
        nominal_rest[leg].P3(),
        leg == 0 ? M3Vec{0.0, 0.0, 1.0} : M3Vec{0.0, 0.0, -1.0});
    photons_rest.emplace_back(energy * direction[0], energy * direction[1],
                              energy * direction[2], energy);
  }

  M4Vec photon_sum;
  for (const auto &photon : photons_rest) {
    photon_sum += photon;
  }
  const M4Vec recoil = M4Vec(0.0, 0.0, 0.0, parent_mass) - photon_sum;
  M4Vec hard1;
  M4Vec hard2;
  if (!TwoBody(recoil, nominal_rest[0], nominal_rest[1], beam1.mass, beam2.mass, hard1, hard2)) {
    return 0.0;
  }

  isr_state.photons = photons_rest;
  for (auto &photon : isr_state.photons) {
    kinematics::LorentzBoost(parent, parent_mass, photon, 1);
  }
  kinematics::LorentzBoost(parent, parent_mass, hard1, 1);
  kinematics::LorentzBoost(parent, parent_mass, hard2, 1);

  isr_state.hard = {hard1, hard2};
  isr_state.weight = radiator_weight;
  isr_state.active = !isr_state.photons.empty();
  momentum1 = hard1;
  momentum2 = hard2;
  s = (hard1 + hard2).M2();
  sqrt_s = s > 0.0 ? std::sqrt(s) : 0.0;
  return radiator_weight;
}

// Sample FSR after the Born amplitude without changing its momenta
void MYFS::PrepareFSR(const std::vector<MDecayBranch> &tree,
                      MRandom &random) {
  ResetFSR();
  if (config.fsr != Mode::YFS) {
    return;
  }
  fsr_state.tree = tree;
  PrepareTree(fsr_state.tree, config.fsr_param, random, fsr_state.pairs);
  fsr_state.prepared = true;
}

// Compute the post-FSR tree used for central fiducial decisions
const std::vector<MDecayBranch> &MYFS::FiducialTree(
    const std::vector<MDecayBranch> &born) const noexcept {
  return fsr_state.prepared ? fsr_state.tree : born;
}

// Clear the event-local final state radiation sample
void MYFS::ResetFSR() noexcept { fsr_state = FSRState(); }

// Attach the sampled initial state to the HepMC graph
bool MYFS::AttachISR(HepMC3::GenEvent &event) {
  if (!isr_state.active) {
    return true;
  }
  const std::vector<HepMC3::GenParticlePtr> hard_beams = event.beams();
  if (hard_beams.size() != 2 || hard_beams[0] == nullptr ||
      hard_beams[1] == nullptr) {
    return false;
  }

  auto beam1 = std::make_shared<HepMC3::GenParticle>(
      aux::M4Vec2HepMC3(isr_state.nominal[0]), hard_beams[0]->pid(),
      PDG::PDG_BEAM);
  auto beam2 = std::make_shared<HepMC3::GenParticle>(
      aux::M4Vec2HepMC3(isr_state.nominal[1]), hard_beams[1]->pid(),
      PDG::PDG_BEAM);
  hard_beams[0]->set_status(PDG::PDG_INTERMEDIATE);
  hard_beams[1]->set_status(PDG::PDG_INTERMEDIATE);
  hard_beams[0]->set_momentum(aux::M4Vec2HepMC3(isr_state.hard[0]));
  hard_beams[1]->set_momentum(aux::M4Vec2HepMC3(isr_state.hard[1]));

  auto vertex = std::make_shared<HepMC3::GenVertex>();
  vertex->add_particle_in(beam1);
  vertex->add_particle_in(beam2);
  vertex->add_particle_out(hard_beams[0]);
  vertex->add_particle_out(hard_beams[1]);
  event.add_vertex(vertex);
  for (const auto &momentum : isr_state.photons) {
    auto photon = std::make_shared<HepMC3::GenParticle>(
        aux::M4Vec2HepMC3(momentum), PDG::PDG_gamma, PDG::PDG_STABLE);
    vertex->add_particle_out(photon);
    photon->add_attribute("QED_ISR", std::make_shared<HepMC3::IntAttribute>(1));
  }
  event.set_beam_particles(beam1, beam2);
  return true;
}

// Radiate all supported final state lepton-pair vertices
bool MYFS::ApplyFSR(HepMC3::GenEvent &event, MRandom &random) {
  if (config.fsr != Mode::YFS) {
    return true;
  }
  if (!fsr_state.prepared) {
    const std::vector<HepMC3::GenVertexPtr> vertices = event.vertices();
    for (const auto &vertex : vertices) {
      RadiateVertex(vertex, config.fsr_param, random);
    }
    return true;
  }

  const std::vector<HepMC3::GenVertexPtr> vertices = event.vertices();
  std::vector<bool> used(vertices.size(), false);
  for (const auto &state : fsr_state.pairs) {
    bool matched = false;
    for (std::size_t index = 0; index < vertices.size(); ++index) {
      if (used[index]) {
        continue;
      }
      const int order = MatchPair(vertices[index], state);
      if (order < 0) {
        continue;
      }
      CommitPair(vertices[index], state, order == 1);
      used[index] = true;
      matched = true;
      break;
    }
    if (!matched) {
      return false;
    }
  }
  return true;
}

// Add accepted ISR and FSR photons to one complete event record
bool MYFS::Apply(HepMC3::GenEvent &event, MRandom &random, std::exception_ptr *failure) noexcept {
  if (failure) { *failure = nullptr; }
  try {
    if (!AttachISR(event)) {
      return false;
    }
    return ApplyFSR(event, random);
  } catch (const std::exception &) {
    if (failure) { *failure = std::current_exception(); }
    return false;
  }
}

} // namespace radiative
} // namespace gra
