// Physics tests for QED initial and final state radiation
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <algorithm>
#include <array>
#include <catch.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

// HepMC3
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"

// Libraries
#include "json.hpp"

// Own
#include "Graniitti/Kinematics/M4Vec.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Particle/MPDG.h"
#include "Graniitti/Photon/MRadiative.h"
#include "Graniitti/Process/MProcessState.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MException.h"

namespace {

// Compute one complete numerical block for radiative unit tests
nlohmann::json Numerics(double energy_min = 1.0e-3) {
  return {{"ISR_QED", {{"x_min", 1.0e-6}, {"hard_collinear", true}}},
          {"FSR_QED",
           {{"energy_min", energy_min},
            {"energy_max_fraction", 0.5},
            {"max_photons", 32},
            {"max_trials", 10000}}}};
}

// Compute one on-shell collider beam along a chosen longitudinal direction
gra::M4Vec Beam(double energy, double mass, double sign) {
  return gra::M4Vec(0.0, 0.0, sign * std::sqrt(energy * energy - mass * mass),
                    energy);
}

// Compute the largest absolute component of one four-vector
double MaxComponent(const gra::M4Vec &momentum) {
  return std::max({std::abs(momentum.Px()), std::abs(momentum.Py()),
                   std::abs(momentum.Pz()), std::abs(momentum.E())});
}

// Rotate one four-vector through a fixed generic spatial rotation
gra::M4Vec Rotate(gra::M4Vec momentum) {
  momentum.RotateX(0.37);
  momentum.RotateY(-0.61);
  momentum.RotateZ(0.29);
  return momentum;
}

// Build one neutral parent to stable muon-pair event
HepMC3::GenEvent MuonEvent(double mass) {
  constexpr double muon_mass = 0.1056583745;
  const double momentum = std::sqrt(mass * mass / 4.0 - muon_mass * muon_mass);
  HepMC3::GenEvent event(HepMC3::Units::GEV, HepMC3::Units::MM);
  auto parent = std::make_shared<HepMC3::GenParticle>(
      HepMC3::FourVector(0.0, 0.0, 0.0, mass), gra::PDG::PDG_system,
      gra::PDG::PDG_INTERMEDIATE);
  auto minus = std::make_shared<HepMC3::GenParticle>(
      HepMC3::FourVector(0.0, 0.0, momentum, mass / 2.0), 13,
      gra::PDG::PDG_STABLE);
  auto plus = std::make_shared<HepMC3::GenParticle>(
      HepMC3::FourVector(0.0, 0.0, -momentum, mass / 2.0), -13,
      gra::PDG::PDG_STABLE);
  auto vertex = std::make_shared<HepMC3::GenVertex>();
  vertex->add_particle_in(parent);
  vertex->add_particle_out(minus);
  vertex->add_particle_out(plus);
  event.add_vertex(vertex);
  return event;
}

// Build the decay-tree representation of one stable muon pair
std::vector<gra::MDecayBranch> MuonTree(double mass) {
  constexpr double muon_mass = 0.1056583745;
  const double momentum = std::sqrt(mass * mass / 4.0 - muon_mass * muon_mass);
  gra::MDecayBranch minus;
  minus.p.pdg = 13;
  minus.p.mass = muon_mass;
  minus.p4 = gra::M4Vec(0.0, 0.0, momentum, mass / 2.0);
  gra::MDecayBranch plus;
  plus.p.pdg = -13;
  plus.p.mass = muon_mass;
  plus.p4 = gra::M4Vec(0.0, 0.0, -momentum, mass / 2.0);
  return {minus, plus};
}

// Compute the final photon count and energy while checking vertex closure
std::pair<std::size_t, double> CheckFSR(const HepMC3::GenEvent &event) {
  const auto &vertex = event.vertices().front();
  gra::M4Vec incoming;
  gra::M4Vec outgoing;
  std::size_t photons = 0;
  double photon_energy = 0.0;
  for (const auto &particle : vertex->particles_in()) {
    incoming += gra::aux::HepMC2M4Vec(particle->momentum());
  }
  for (const auto &particle : vertex->particles_out()) {
    const gra::M4Vec momentum = gra::aux::HepMC2M4Vec(particle->momentum());
    outgoing += momentum;
    if (particle->pid() == gra::PDG::PDG_gamma) {
      REQUIRE(particle->attribute<HepMC3::IntAttribute>("QED_FSR") != nullptr);
      ++photons;
      photon_energy += momentum.E();
    } else {
      REQUIRE(momentum.M2() ==
              Approx(0.1056583745 * 0.1056583745).margin(2.0e-12));
    }
  }
  REQUIRE(MaxComponent(incoming - outgoing) < 2.0e-11);
  return {photons, photon_energy};
}

} // namespace

// Check explicit disabled modes and reject unsupported radiative steering
TEST_CASE("Radiative steering rejects unknown modes and invalid controls",
          "[radiative][steering]") {
  const auto disabled = gra::radiative::ReadConfig({{"ISR_QED", "none"}, {"FSR_QED", "none"}}, Numerics());
  REQUIRE(disabled.isr == gra::radiative::Mode::Off);
  REQUIRE(disabled.fsr == gra::radiative::Mode::Off);
  for (const auto &key : {"ISR_QED", "FSR_QED"}) {
    REQUIRE_THROWS_AS(gra::radiative::ReadConfig({{key, "OFF"}}, Numerics()), std::invalid_argument);
  }
  REQUIRE_THROWS_AS(gra::radiative::ReadConfig(
                        {{"ISR_QED", "YFS2"}, {"FSR_QED", "none"}}, Numerics()),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(gra::radiative::ReadConfig(
                        {{"ISR_QED", "none"}, {"FSR", "YFS"}}, Numerics()),
                    std::invalid_argument);
  auto numerics = Numerics();
  numerics["FSR_QED"]["energy_min"] = -1.0;
  REQUIRE_THROWS_AS(gra::radiative::ReadConfig({{"ISR_QED", "none"}, {"FSR_QED", "YFS"}}, numerics),
                    std::invalid_argument);
}

// Check that active YFS sectors have a charged radiator in the process
TEST_CASE("YFS steering rejects processes without charged radiators", "[radiative][steering]") {
  gra::MParticle proton;
  proton.pdg = 2212;
  gra::MParticle electron;
  electron.pdg = 11;

  gra::radiative::MYFS yfs;
  yfs.Configure(gra::radiative::ReadConfig({{"ISR_QED", "YFS"}, {"FSR_QED", "none"}}, Numerics()));
  REQUIRE_THROWS_AS(yfs.ValidateApplicability(proton, proton, {}), std::invalid_argument);
  REQUIRE_NOTHROW(yfs.ValidateApplicability(electron, proton, {}));

  yfs.Configure(gra::radiative::ReadConfig({{"ISR_QED", "none"}, {"FSR_QED", "YFS"}}, Numerics()));
  REQUIRE_THROWS_AS(yfs.ValidateApplicability(proton, proton, {}), std::invalid_argument);
  REQUIRE_NOTHROW(yfs.ValidateApplicability(proton, proton, MuonTree(10.0)));
}

TEST_CASE("YFS massive dipole integrates to its analytic soft factor", "[radiative][FSR][physics]") {
  for (const double beta : {1.0e-8, 1.0e-4, 0.05, 0.5, 0.95, 0.999}) {
    constexpr int intervals = 200000;
    const double step = 2.0 / static_cast<double>(intervals);
    double integral = 0.0;
    for (int index = 0; index <= intervals; ++index) {
      const double cosine = -1.0 + step * static_cast<double>(index);
      const double coefficient = index == 0 || index == intervals ? 1.0
                                 : index % 2 == 0                 ? 2.0
                                                                  : 4.0;
      integral += coefficient * gra::radiative::DipoleAngular(beta, cosine);
    }
    integral *= step / 3.0;
    REQUIRE(integral ==
            Approx(2.0 * gra::radiative::DipoleIntegral(beta)).epsilon(2.0e-9));
  }
}

// Check dipole suppression at pair threshold without subtracting large terms
TEST_CASE("Massive dipole has the nonrelativistic radiation limit", "[radiative][FSR][threshold]") {
  for (const double beta : {1.0e-10, 1.0e-8, 1.0e-5}) {
    REQUIRE(gra::radiative::DipoleIntegral(beta) / (beta * beta) == Approx(8.0 / 3.0).epsilon(2.0e-9));
    for (const double c : {-1.0, -0.7, 0.0, 0.7, 1.0}) {
      REQUIRE(gra::radiative::DipoleAngular(beta, c) / (beta * beta) ==
              Approx(4.0 * (1.0 - c * c)).margin(2.0e-9));
    }
  }
}

TEST_CASE(
    "Final state recoil is on shell conserving and rotationally covariant",
    "[radiative][FSR][covariance]") {
  constexpr double mass = 0.1056583745;
  const double momentum = std::sqrt(25.0 - mass * mass);
  const gra::M4Vec lepton1(0.0, 0.0, momentum, 5.0);
  const gra::M4Vec lepton2(0.0, 0.0, -momentum, 5.0);
  const std::vector<gra::M4Vec> photons = {
      gra::M4Vec(0.15, 0.20, 0.10, std::sqrt(0.0725)),
      gra::M4Vec(-0.04, 0.02, -0.03, std::sqrt(0.0029))};

  gra::M4Vec mapped1;
  gra::M4Vec mapped2;
  REQUIRE(gra::radiative::MapLeptonPair(lepton1, lepton2, photons, mapped1,
                                        mapped2));
  gra::M4Vec photon_sum;
  for (const auto &photon : photons) {
    photon_sum += photon;
  }
  REQUIRE(MaxComponent(lepton1 + lepton2 - mapped1 - mapped2 - photon_sum) <
          2.0e-12);
  REQUIRE(mapped1.M2() == Approx(mass * mass).margin(2.0e-12));
  REQUIRE(mapped2.M2() == Approx(mass * mass).margin(2.0e-12));

  // Exchanging the lepton charges must exchange the mapped momenta
  gra::M4Vec swapped1, swapped2;
  REQUIRE(gra::radiative::MapLeptonPair(lepton2, lepton1, photons, swapped1, swapped2));
  REQUIRE(MaxComponent(swapped1 - mapped2) < 2.0e-12);
  REQUIRE(MaxComponent(swapped2 - mapped1) < 2.0e-12);

  // The relative recoil axis must also transform under a longitudinal boost
  const gra::M4Vec boost(0.0, 0.0, std::sinh(0.7), std::cosh(0.7));
  auto boosted1 = lepton1;
  auto boosted2 = lepton2;
  auto expected1 = mapped1;
  auto expected2 = mapped2;
  auto boosted_photons = photons;
  gra::kinematics::LorentzBoost(boost, 1.0, boosted1, 1);
  gra::kinematics::LorentzBoost(boost, 1.0, boosted2, 1);
  gra::kinematics::LorentzBoost(boost, 1.0, expected1, 1);
  gra::kinematics::LorentzBoost(boost, 1.0, expected2, 1);
  for (auto &photon : boosted_photons) { gra::kinematics::LorentzBoost(boost, 1.0, photon, 1); }
  REQUIRE(gra::radiative::MapLeptonPair(boosted1, boosted2, boosted_photons, swapped1, swapped2));
  REQUIRE(MaxComponent(swapped1 - expected1) < 2.0e-12);
  REQUIRE(MaxComponent(swapped2 - expected2) < 2.0e-12);

  std::vector<gra::M4Vec> rotated_photons;
  for (const auto &photon : photons) {
    rotated_photons.push_back(Rotate(photon));
  }
  gra::M4Vec rotated1;
  gra::M4Vec rotated2;
  REQUIRE(gra::radiative::MapLeptonPair(Rotate(lepton1), Rotate(lepton2),
                                        rotated_photons, rotated1, rotated2));
  REQUIRE(MaxComponent(rotated1 - Rotate(mapped1)) < 2.0e-12);
  REQUIRE(MaxComponent(rotated2 - Rotate(mapped2)) < 2.0e-12);
}

// Check physical input momenta independently of the available recoil mass
TEST_CASE("FSR rejects unphysical external momenta", "[radiative][FSR][kinematics]") {
  constexpr double mass = 0.1056583745;
  const double p = std::sqrt(100.0 - mass * mass);
  const gra::M4Vec first(0.0, 0.0, p, 10.0);
  const gra::M4Vec second(0.0, 0.0, -p, 10.0);
  gra::M4Vec mapped1, mapped2;
  for (const auto &photon : {gra::M4Vec(1.0, 0.0, 0.0, -1.0),
                            gra::M4Vec(0.0, 0.0, 0.0, 1.0),
                            gra::M4Vec(2.0, 0.0, 0.0, 1.0)}) {
    REQUIRE_FALSE(gra::radiative::MapLeptonPair(first, second, {photon}, mapped1, mapped2));
  }
  const gra::M4Vec spacelike(0.0, 0.0, 2.0, 1.0);
  REQUIRE_FALSE(gra::radiative::MapLeptonPair(spacelike, second, {}, mapped1, mapped2));
  REQUIRE_FALSE(gra::radiative::MapLeptonPair(first, spacelike, {}, mapped1, mapped2));
  // Empty radiation preserves a physical Born pair
  REQUIRE(gra::radiative::MapLeptonPair(first, second, {}, mapped1, mapped2));
  REQUIRE(MaxComponent(mapped1 - first) < 2.0e-12);
  REQUIRE(MaxComponent(mapped2 - second) < 2.0e-12);
}

TEST_CASE("FSR rejects past directed timelike recoil", "[radiative][FSR][kinematics]") {
  constexpr double mass = 0.1056583745;
  const double p = std::sqrt(100.0 - mass * mass);
  const gra::M4Vec first(0.0, 0.0, p, 10.0);
  const gra::M4Vec second(0.0, 0.0, -p, 10.0);
  std::vector<gra::M4Vec> photons;
  for (int i = 0; i < 3; ++i) {
    photons.emplace_back(8.0, 0.0, 0.0, 8.0);
    photons.emplace_back(-8.0, 0.0, 0.0, 8.0);
  }
  gra::M4Vec mapped1, mapped2;
  REQUIRE_FALSE(gra::radiative::MapLeptonPair(first, second, photons, mapped1, mapped2));
  // A longitudinal boost cannot make the unphysical proposal acceptable
  const gra::M4Vec boost(0.0, 0.0, std::sinh(0.7), std::cosh(0.7));
  auto boosted1 = first;
  auto boosted2 = second;
  gra::kinematics::LorentzBoost(boost, 1.0, boosted1, 1);
  gra::kinematics::LorentzBoost(boost, 1.0, boosted2, 1);
  for (auto &photon : photons) { gra::kinematics::LorentzBoost(boost, 1.0, photon, 1); }
  REQUIRE_FALSE(gra::radiative::MapLeptonPair(boosted1, boosted2, photons, mapped1, mapped2));
}

TEST_CASE("Prepared FSR controls fiducial migration and HepMC recoil",
          "[radiative][FSR][fiducial]") {
  const auto config = gra::radiative::ReadConfig(
      {{"ISR_QED", "none"}, {"FSR_QED", "YFS"}}, Numerics(1.0e-4));
  gra::radiative::MYFS yfs;
  yfs.Configure(config);
  gra::MRandom random;
  random.SetSeed(42017);

  const auto born = MuonTree(20.0);
  for (std::size_t trial = 0; trial < 1000 && yfs.GetFSR().pairs.empty();
       ++trial) {
    yfs.PrepareFSR(born, random);
  }
  REQUIRE(yfs.GetFSR().prepared);
  REQUIRE(yfs.GetFSR().pairs.size() == 1);

  const auto cached = yfs.GetFSR();
  const auto &radiated = yfs.FiducialTree(born);
  const double born_mass = (born[0].p4 + born[1].p4).M();
  const double radiated_mass = (radiated[0].p4 + radiated[1].p4).M();
  REQUIRE(radiated_mass < born_mass);

  gra::FIDPDGCUT pair_cut;
  pair_cut.pdg = {13, -13};
  pair_cut.pdg_abs = {false, false};
  pair_cut.M.active = true;
  pair_cut.M.min = 0.0;
  pair_cut.M.max = 0.5 * (born_mass + radiated_mass);
  gra::FIDCUT cuts;
  cuts.pdg_cuts = {pair_cut};
  REQUIRE_FALSE(cuts.PassSelectedParticles(born));
  REQUIRE(cuts.PassSelectedParticles(radiated));

  HepMC3::GenEvent event = MuonEvent(20.0);
  random.SetSeed(99181);
  REQUIRE(yfs.Apply(event, random));
  const auto [photons, energy] = CheckFSR(event);
  REQUIRE(photons == cached.pairs[0].photons.size());
  REQUIRE(energy > 0.0);

  HepMC3::GenParticlePtr minus;
  HepMC3::GenParticlePtr plus;
  for (const auto &particle : event.vertices().front()->particles_out()) {
    if (particle->pid() == 13) {
      minus = particle;
    } else if (particle->pid() == -13) {
      plus = particle;
    }
  }
  REQUIRE(minus != nullptr);
  REQUIRE(plus != nullptr);
  REQUIRE(MaxComponent(gra::aux::HepMC2M4Vec(minus->momentum()) -
                       cached.pairs[0].lepton[0]) < 1.0e-12);
  REQUIRE(MaxComponent(gra::aux::HepMC2M4Vec(plus->momentum()) -
                       cached.pairs[0].lepton[1]) < 1.0e-12);
}

TEST_CASE("Lepton ISR constructs an exact reduced hard collision",
          "[radiative][ISR][physics]") {
  const auto config = gra::radiative::ReadConfig(
      {{"ISR_QED", "YFS"}, {"FSR_QED", "none"}}, Numerics());
  gra::radiative::MYFS yfs;
  yfs.Configure(config);

  gra::MParticle electron;
  electron.pdg = 11;
  electron.mass = 0.0005109989461;
  gra::MParticle positron = electron;
  positron.pdg = -11;
  gra::M4Vec beam1 = Beam(45.6, electron.mass, 1.0);
  gra::M4Vec beam2 = Beam(45.6, electron.mass, -1.0);
  const gra::M4Vec nominal = beam1 + beam2;
  yfs.SetBeams(beam1, beam2);

  double s = nominal.M2();
  double sqrt_s = nominal.M();
  const double weight = yfs.GenerateISR({0.35, 0.72}, 0, electron, positron,
                                        beam1, beam2, s, sqrt_s);
  REQUIRE(weight > 0.0);
  REQUIRE(std::isfinite(weight));
  REQUIRE(beam1.M2() == Approx(electron.mass * electron.mass).margin(1.0e-12));
  REQUIRE(beam2.M2() == Approx(positron.mass * positron.mass).margin(1.0e-12));
  gra::M4Vec photons;
  for (const auto &photon : yfs.GetISR().photons) {
    REQUIRE(photon.M2() == Approx(0.0).margin(1.0e-12));
    photons += photon;
  }
  REQUIRE(MaxComponent(nominal - beam1 - beam2 - photons) < 2.0e-11);
  REQUIRE(s == Approx((beam1 + beam2).M2()).margin(1.0e-12));
  REQUIRE(sqrt_s == Approx((beam1 + beam2).M()).margin(1.0e-12));

  HepMC3::GenEvent event(HepMC3::Units::GEV, HepMC3::Units::MM);
  auto hard1 = std::make_shared<HepMC3::GenParticle>(
      gra::aux::M4Vec2HepMC3(beam1), 11, gra::PDG::PDG_BEAM);
  auto hard2 = std::make_shared<HepMC3::GenParticle>(
      gra::aux::M4Vec2HepMC3(beam2), -11, gra::PDG::PDG_BEAM);
  event.add_particle(hard1);
  event.add_particle(hard2);
  event.set_beam_particles(hard1, hard2);
  gra::MRandom random;
  REQUIRE(yfs.Apply(event, random));
  REQUIRE(event.beams().size() == 2);
  REQUIRE(event.beams()[0]->status() == gra::PDG::PDG_BEAM);
  REQUIRE(event.beams()[1]->status() == gra::PDG::PDG_BEAM);
  REQUIRE(hard1->status() == gra::PDG::PDG_INTERMEDIATE);
  REQUIRE(hard2->status() == gra::PDG::PDG_INTERMEDIATE);
  const bool has_isr_vertex = std::any_of(
      event.vertices().begin(), event.vertices().end(), [](const auto &vertex) {
        return vertex->particles_in().size() == 2 &&
               vertex->particles_out().size() >= 3;
      });
  REQUIRE(has_isr_vertex);
  for (const auto &particle : event.particles()) {
    if (particle->pid() == gra::PDG::PDG_gamma) {
      REQUIRE(particle->attribute<HepMC3::IntAttribute>("QED_ISR") != nullptr);
    }
  }
}

// Check analytic Mellin moments of the soft and second order ISR distributions
TEST_CASE("ISR sampling reproduces structure function moments",
          "[radiative][ISR][normalization]") {
  const bool hard_collinear = GENERATE(false, true);
  auto numerics = Numerics();
  numerics["ISR_QED"]["hard_collinear"] = hard_collinear;
  numerics["ISR_QED"]["x_min"] = 1.0e-10;
  const auto config = gra::radiative::ReadConfig(
      {{"ISR_QED", "YFS"}, {"FSR_QED", "none"}}, numerics);
  gra::radiative::MYFS yfs;
  yfs.Configure(config);

  gra::MParticle electron;
  electron.pdg = 11;
  electron.mass = 0.0005109989461;
  gra::MParticle photon;
  photon.pdg = 22;
  photon.mass = 0.0;
  const gra::M4Vec nominal1 = Beam(100.0, electron.mass, 1.0);
  const gra::M4Vec nominal2 = Beam(100.0, photon.mass, -1.0);
  yfs.SetBeams(nominal1, nominal2);

  constexpr std::size_t samples = 100000;
  std::array<double, 3> integral{};
  for (std::size_t sample = 0; sample < samples; ++sample) {
    gra::M4Vec beam1 = nominal1;
    gra::M4Vec beam2 = nominal2;
    double s = (beam1 + beam2).M2();
    double sqrt_s = std::sqrt(s);
    const double unit =
        (static_cast<double>(sample) + 0.5) / static_cast<double>(samples);
    const double weight = yfs.GenerateISR({unit}, 0, electron, photon, beam1, beam2, s, sqrt_s);
    const double x = s / (nominal1 + nominal2).M2();
    for (const auto &k : gra::aux::indices(integral)) {
      integral[k] += weight * std::pow(x, static_cast<int>(k)) / static_cast<double>(samples);
    }
  }

  constexpr double euler = 0.57721566490153286061;
  const double exponent =
      gra::radiative::ISRExponent((nominal1 + nominal2).M2(), electron.mass);
  // Integrate the finite splitting terms analytically over x in [0,1]
  const std::array<double, 3> second = {2.0 * gra::math::PIPI / 3.0 - 9.0 / 4.0,
                                      2.0 * gra::math::PIPI / 3.0 - 89.0 / 36.0,
                                      2.0 * gra::math::PIPI / 3.0 - 419.0 / 144.0};
  for (const auto &k : gra::aux::indices(integral)) {
    const double n = static_cast<double>(k + 1);
    double expected = std::exp(exponent * (0.75 - euler)) * std::tgamma(n) / std::tgamma(n + exponent);
    if (hard_collinear) {
      expected -= 0.5 * exponent * (1.0 / n + 1.0 / (n + 1.0));
      expected += exponent * exponent * second[k] / 8.0;
    }
    REQUIRE(integral[k] == Approx(expected).margin(5.0e-8));
  }
}

TEST_CASE("Exclusive FSR is unitary and stable under the resolved cutoff",
          "[radiative][FSR][unitarity]") {
  constexpr std::size_t events = 12000;
  std::array<double, 2> mean_energy = {0.0, 0.0};
  std::array<std::size_t, 2> radiative_events = {0, 0};
  const std::array<double, 2> cutoffs = {1.0e-2, 1.0e-4};

  for (std::size_t sample = 0; sample < cutoffs.size(); ++sample) {
    const auto config = gra::radiative::ReadConfig(
        {{"ISR_QED", "none"}, {"FSR_QED", "YFS"}}, Numerics(cutoffs[sample]));
    gra::radiative::MYFS yfs;
    yfs.Configure(config);
    gra::MRandom random;
    random.SetSeed(static_cast<std::uint32_t>(9817 + sample));

    for (std::size_t event_index = 0; event_index < events; ++event_index) {
      HepMC3::GenEvent event = MuonEvent(20.0);
      REQUIRE(yfs.Apply(event, random));
      const auto [photons, energy] = CheckFSR(event);
      radiative_events[sample] += static_cast<std::size_t>(photons > 0);
      mean_energy[sample] += energy / static_cast<double>(events);
    }
  }

  REQUIRE(radiative_events[0] > events / 10);
  REQUIRE(radiative_events[1] > radiative_events[0]);
  REQUIRE(mean_energy[1] == Approx(mean_energy[0]).epsilon(0.08));
}

// Keep numerical multiplicity limits separate from the physical no-photon probability
TEST_CASE("FSR multiplicity exhaustion remains a nonfatal phase-space failure", "[radiative][FSR][failures]") {
  auto config = gra::radiative::ReadConfig(
      {{"ISR_QED", "none"}, {"FSR_QED", "YFS"}}, Numerics(1.0e-4));
  config.fsr_param.max_photons = 1;
  gra::radiative::MYFS yfs;
  yfs.Configure(config);
  gra::MRandom random;
  random.SetSeed(73519);
  const auto born = MuonTree(20.0);
  std::size_t failures = 0;
  for (std::size_t i = 0; i < 500; ++i) {
    try {
      yfs.PrepareFSR(born, random);
      REQUIRE(yfs.GetFSR().prepared);
    } catch (const gra::PhaseSpaceFailure &) {
      ++failures;
      REQUIRE_FALSE(yfs.GetFSR().prepared);
    }
  }
  REQUIRE(failures > 0);
}
