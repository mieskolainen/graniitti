// HepMC3 nuclear reaction interface and exclusive kinematic conservation tests
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include <catch.hpp>
#include <sstream>

#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/Nuclear/MFinal.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MException.h"
#include "HepMC3/Attribute.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/WriterAscii.h"

namespace {

enum class Fault { None, Zero, Momentum, Charge, Unresolved, Hard, Orphan, Time, Weight, Throw };

// Add an isotropic two-body decay with the existing phase-space generator
HepMC3::GenParticlePtr Decay(HepMC3::GenEvent& event, const HepMC3::GenParticlePtr& parent, int id1, double m1, int id2,
                             double m2, gra::MRandom& random) {
  const auto&             q = parent->momentum();
  const gra::M4Vec        momentum(q.px(), q.py(), q.pz(), q.e());
  std::vector<gra::M4Vec> daughters;
  const auto              weight =
      gra::kinematics::TwoBodyPhaseSpace(momentum, parent->generated_mass(), {m1, m2}, daughters, random);
  REQUIRE(weight.GetW() > 0.0);
  auto vertex = std::make_shared<HepMC3::GenVertex>();
  parent->set_status(2);
  vertex->add_particle_in(parent);
  auto first  = std::make_shared<HepMC3::GenParticle>(gra::aux::M4Vec2HepMC3(daughters[0]), id1, 1);
  auto second = std::make_shared<HepMC3::GenParticle>(gra::aux::M4Vec2HepMC3(daughters[1]), id2, 1);
  first->set_generated_mass(m1);
  second->set_generated_mass(m2);
  vertex->add_particle_out(first);
  vertex->add_particle_out(second);
  event.add_vertex(vertex);
  return second;
}

// Use a known radiative deuteron breakup solely to test the public backend API
class AnalyticDecay final : public gra::nuclear::MReaction {
 public:
  // Select one controlled physical or malformed completion
  explicit AnalyticDecay(Fault fault = Fault::None) : fault_(fault) {}

  // Clone the test model independently
  std::unique_ptr<gra::nuclear::MReaction> Clone() const override { return std::make_unique<AnalyticDecay>(fault_); }

  // Identify the analytic phase-space fixture
  HepMC3::GenRunInfo::ToolInfo Tool() const override {
    return {"analytic deuteron breakup", "test", "Radiative phase-space fixture"};
  }

  // Fix the test parent mass above the proton-neutron threshold
  gra::nuclear::RecoilMass SampleMasses(const HepMC3::GenEvent&, gra::MRandom&) override { return {{2.1, 2.1}, {}}; }

  // Generate gamma, proton and neutron products at fixed parent momenta
  double Complete(HepMC3::GenEvent& event, const std::array<HepMC3::GenParticlePtr, 2>& roots,
                  gra::MRandom& random) override {
    for (const auto& parent : roots) {
      parent->set_pid(1000010020);
      auto excited = Decay(event, parent, 22, 0.0, 1000010020, 1.94, random);
      Decay(event, excited, 2212, gra::PDG::mp, 2112, 0.9395654, random);
    }
    const auto last = event.particles().back();
    switch (fault_) {
      case Fault::None:
        break;
      case Fault::Zero:
        return 0.0;
      case Fault::Momentum:
        last->set_momentum(last->momentum() + HepMC3::FourVector(0.1, 0, 0, 0));
        break;
      case Fault::Charge:
        last->set_pid(2212);
        break;
      case Fault::Unresolved:
        last->set_status(3);
        break;
      case Fault::Hard:
        event.vertices().front()->set_position(HepMC3::FourVector(1, 0, 0, 0));
        break;
      case Fault::Orphan:
        event.add_particle(std::make_shared<HepMC3::GenParticle>(HepMC3::FourVector(0, 0, 1, 1), 22, 1));
        break;
      case Fault::Time:
        last->production_vertex()->set_position(HepMC3::FourVector(10, 0, 0, 1));
        break;
      case Fault::Weight:
        return -1.0;
      case Fault::Throw:
        throw std::runtime_error("External reaction failed");
    }
    return 0.25;
  }

 private:
  Fault fault_;
};

// Construct conserving hard vertices for boosted forward mass proposals
HepMC3::GenEvent HardEvent(double pz, double phi = 0.0) {
  HepMC3::GenEvent   event(HepMC3::Units::GEV, HepMC3::Units::MM);
  auto               central = std::make_shared<HepMC3::GenVertex>();
  HepMC3::FourVector sum;
  for (const auto leg : {1, 2}) {
    const int                side = leg == 1 ? 1 : -1;
    const double             z    = side * pz;
    const HepMC3::FourVector incoming(0, 0, z, std::hypot(z, 1.8756));
    const HepMC3::FourVector outgoing(side * 0.1 * std::cos(phi), side * 0.1 * std::sin(phi), 0.99 * z,
                                      std::sqrt(0.99 * 0.99 * z * z + 0.01 + 2.1 * 2.1));
    auto                     beam = std::make_shared<HepMC3::GenParticle>(incoming, 1000010020, 4);
    beam->set_generated_mass(1.8756);
    auto forward  = std::make_shared<HepMC3::GenParticle>(outgoing, 91, 3);
    auto transfer = std::make_shared<HepMC3::GenParticle>(incoming - outgoing, 99, 3);
    auto vertex   = std::make_shared<HepMC3::GenVertex>();
    vertex->add_particle_in(beam);
    vertex->add_particle_out(forward);
    vertex->add_particle_out(transfer);
    event.add_vertex(vertex);
    forward->add_attribute("graniitti_upc_leg", std::make_shared<HepMC3::IntAttribute>(leg));
    central->add_particle_in(transfer);
    sum = sum + transfer->momentum();
  }
  central->add_particle_out(std::make_shared<HepMC3::GenParticle>(sum, 90, 3));
  event.add_vertex(central);
  event.weights().push_back(2.0);
  return event;
}

}  // namespace

TEST_CASE("External nuclear decay conserves boosted gamma proton neutron branches", "[nuclear][final]") {
  gra::MRandom random;
  random.SetSeed(724);
  gra::nuclear::MFinal model(std::make_unique<AnalyticDecay>());
  gra::MPDG            pdg;
  pdg.ReadParticleData();
  for (const auto pz : {100.0, 500000.0, -500000.0}) {
    for (const auto phi : {0.0, 0.7, 2.4}) {
      auto event = HardEvent(pz, phi);
      REQUIRE(model.SampleMasses(event, random).mass[0] == Approx(2.1));
      REQUIRE(model.Complete(event, random) == Approx(0.25));
      REQUIRE(event.weights().front() == Approx(2.0));
      std::stringstream   buffer;
      HepMC3::WriterAscii writer(buffer);
      writer.set_precision(17);
      writer.write_event(event);
      writer.close();
      HepMC3::ReaderAscii reader(buffer);
      HepMC3::GenEvent    decoded;
      reader.read_event(decoded);
      REQUIRE_FALSE(reader.failed());
      for (const auto& root : gra::nuclear::ForwardParents(decoded)) {
        REQUIRE_NOTHROW(gra::nuclear::ValidateNuclearDecay(decoded, root, pdg));
      }
      REQUIRE_THROWS_AS(model.Complete(event, random), gra::PhaseSpaceFailure);
    }
  }
}

TEST_CASE("External nuclear failures remain sample failures and clones have no proposal", "[nuclear][final]") {
  gra::MRandom random;
  random.SetSeed(725);
  for (const auto fault : {Fault::Momentum, Fault::Charge, Fault::Unresolved, Fault::Hard, Fault::Orphan, Fault::Time,
                           Fault::Weight, Fault::Throw}) {
    auto                 event = HardEvent(100.0);
    gra::nuclear::MFinal model(std::make_unique<AnalyticDecay>(fault));
    model.SampleMasses(event, random);
    auto worker(model);
    REQUIRE_THROWS_AS(worker.Complete(event, random), gra::PhaseSpaceFailure);
    REQUIRE_THROWS_AS(model.Complete(event, random), gra::PhaseSpaceFailure);
  }
}

TEST_CASE("Zero nuclear response has zero support without an invalid sample", "[nuclear][final]") {
  gra::MRandom         random;
  auto                 event = HardEvent(100.0);
  gra::nuclear::MFinal model(std::make_unique<AnalyticDecay>(Fault::Zero));
  model.SampleMasses(event, random);
  CHECK(model.Complete(event, random) == Approx(0.0));
  REQUIRE_THROWS_AS(model.Complete(event, random), gra::PhaseSpaceFailure);
}
