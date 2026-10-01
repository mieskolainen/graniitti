// ROOT analyzer normalization unit tests
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include <catch.hpp>

#include <utility>
#include <vector>

#include "Graniitti/Analysis/MAnalyzer.h"
#include "Graniitti/Program/Analysis/analyze.h"
#include "Graniitti/Math/MSpecialFunctions.h"
#include "Graniitti/Analysis/MMultiplet.h"
#include "support/analysis_test_support.hh"
#include "Graniitti/Analysis/MROOT.h"
#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/Sampling/MRandom.h"

using gra::aux::indices;

// Check that unweighted conditional histograms retain their event acceptance
TEST_CASE("ROOT analyzer preserves unweighted conditional cross sections",
          "[gra::analyzer][normalization]") {
  gra::h1Multiplet inclusive("unweighted_inclusive", "", 2, 0.0, 2.0,
                             {"sample"});
  gra::h1Multiplet conditional("unweighted_conditional", "", 2, 0.0, 2.0,
                               {"sample"});

  inclusive.h[0]->Fill(0.25);
  inclusive.h[0]->Fill(0.75);
  inclusive.h[0]->Fill(1.25);
  inclusive.h[0]->Fill(1.75);
  conditional.h[0]->Fill(0.25);
  conditional.h[0]->Fill(1.25);

  const std::vector<double> normalization = {10.0 / 4.0};
  const std::vector<double> multiplier = {1.0};
  inclusive.NormalizeAll(normalization, multiplier);
  conditional.NormalizeAll(normalization, multiplier);

  REQUIRE(inclusive.h[0]->Integral(0, 3, "width") == Approx(10.0));
  REQUIRE(conditional.h[0]->Integral(0, 3, "width") == Approx(5.0));
}

// Check weighted one- and two-dimensional conditional normalization
TEST_CASE("ROOT analyzer preserves weighted conditional cross sections",
          "[gra::analyzer][normalization]") {
  gra::h1Multiplet inclusive_1d("weighted_inclusive_1d", "", 2, 0.0, 2.0,
                                {"sample"});
  gra::h1Multiplet conditional_1d("weighted_conditional_1d", "", 2, 0.0, 2.0,
                                  {"sample"});
  gra::h2Multiplet inclusive_2d("weighted_inclusive_2d", "", 2, 0.0, 2.0, 2,
                                0.0, 2.0, {"sample"});
  gra::h2Multiplet conditional_2d("weighted_conditional_2d", "", 2, 0.0, 2.0, 2,
                                  0.0, 2.0, {"sample"});

  const std::vector<double> weights = {1.0, 2.0, 3.0};
  for (std::size_t i = 0; i < weights.size(); ++i) {
    const double x = 0.25 + 0.5 * static_cast<double>(i);
    inclusive_1d.h[0]->Fill(x, weights[i]);
    inclusive_2d.h[0]->Fill(x, x, weights[i]);
    if (i < 2) {
      conditional_1d.h[0]->Fill(x, weights[i]);
      conditional_2d.h[0]->Fill(x, x, weights[i]);
    }
  }

  const std::vector<double> normalization = {12.0 / 6.0};
  const std::vector<double> multiplier = {1.0};
  inclusive_1d.NormalizeAll(normalization, multiplier);
  conditional_1d.NormalizeAll(normalization, multiplier);
  inclusive_2d.NormalizeAll(normalization, multiplier);
  conditional_2d.NormalizeAll(normalization, multiplier);

  REQUIRE(inclusive_1d.h[0]->Integral(0, 3, "width") == Approx(12.0));
  REQUIRE(conditional_1d.h[0]->Integral(0, 3, "width") == Approx(6.0));
  REQUIRE(inclusive_2d.h[0]->Integral(0, 3, 0, 3, "width") == Approx(12.0));
  REQUIRE(conditional_2d.h[0]->Integral(0, 3, 0, 3, "width") == Approx(6.0));
}

// Check that multiplet histograms are not owned by a ROOT directory
TEST_CASE("ROOT analyzer uses explicit histogram ownership",
          "[gra::analyzer][ownership]") {
  gra::rootstyle::SetROOTStyle();
  gra::h1Multiplet histogram_1d("detached_1d", "", 1, 0.0, 1.0, {"sample"});
  gra::h2Multiplet histogram_2d("detached_2d", "", 1, 0.0, 1.0, 1, 0.0, 1.0,
                                {"sample"});
  gra::hProfMultiplet profile("detached_profile", "", 1, 0.0, 1.0, -1.0, 1.0,
                              {"sample"});

  REQUIRE_FALSE(TH1::AddDirectoryStatus());
  REQUIRE(histogram_1d.h[0]->GetDirectory() == nullptr);
  REQUIRE(histogram_2d.h[0]->GetDirectory() == nullptr);
  REQUIRE(profile.h[0]->GetDirectory() == nullptr);
}

// Check ROOT helper numerical domains and multiplet input dimensions
TEST_CASE("ROOT analysis helpers validate numerical dimensions",
          "[gra::analyzer][validation]") {
  CHECK_THROWS_AS(gra::rootstyle::CubeHelix(1, 0.5, -1.5, 1.2, 1.0),
                  std::invalid_argument);
  const auto colors = gra::rootstyle::CubeHelix(8, 0.5, -1.5, 1.2, 1.0);
  REQUIRE(colors.size() == 4);
  CHECK(colors[0].front() == Approx(0.0));
  CHECK(colors[0].back() == Approx(1.0));
  for (const auto &channel : colors) {
    REQUIRE(channel.size() == 8);
    for (const double value : channel) {
      CHECK(value >= 0.0);
      CHECK(value <= 1.0);
    }
  }

  gra::h1Multiplet histogram_1d("shape_1d", "", 2, 0.0, 1.0,
                                {"first", "second"});
  CHECK_THROWS_AS(histogram_1d.MultiFill({0.5}), std::invalid_argument);

  gra::h2Multiplet histogram_2d("shape_2d", "", 2, 0.0, 1.0, 2, 0.0, 1.0,
                                {"sample"});
  CHECK_THROWS_AS(histogram_2d.MultiFill({{0.5}}), std::invalid_argument);
}

namespace {

struct AnalyzerHistograms {
  std::map<std::string, std::shared_ptr<gra::h1Multiplet>> h1;
  std::map<std::string, std::shared_ptr<gra::h2Multiplet>> h2;
  std::map<std::string, std::shared_ptr<gra::hProfMultiplet>> hP;

  // Book the public analyzer histogram set for one sample
  AnalyzerHistograms() {
    std::vector<std::string> names1 = {"h1_1B_eta", "h1_1B_pt", "h1_PP_dphi", "h1_PP_dpt", "h1_PP_t1",
                                     "h1_S_M", "h1_S_Pt", "h1_S_Pt2", "h1_S_Y", "h1_2B_acop", "h1_2B_diffrap"};
    std::vector<std::string> names2 = {"h2_S_M_Pt", "h2_S_M_dphipp", "h2_S_M_dpt", "h2_S_M_pt", "h2_S_M_t",
                                     "h2_2B_M_dphi", "h2_2B_eta1_eta2"};
    for (const auto &frame : gra::analyzer::Frames()) {
      names1.push_back("h1_costheta_" + frame);
      names1.push_back("h1_phi_" + frame);
      names2.push_back("h2_2B_M_costheta_" + frame);
      names2.push_back("h2_2B_M_phi_" + frame);
      names2.push_back("h2_2B_costheta_phi_" + frame);
    }
    const std::vector<std::string> labels = {"sample"};
    for (const auto &name : names1) {
      h1[name] = std::make_shared<gra::h1Multiplet>(name, "", 10, -10.0, 10.0, labels);
    }
    for (const auto &name : names2) {
      h2[name] = std::make_shared<gra::h2Multiplet>(name, "", 10, -10.0, 10.0, 10, -10.0, 10.0, labels);
    }
    for (const std::string name : {"hP_S_M_Pt", "hP_2B_M_dphi", "hP_S_M_PL2_CM", "hP_S_M_PL4_CM"}) {
      hP[name] = std::make_shared<gra::hProfMultiplet>(name, "", 10, -10.0, 10.0, -10.0, 10.0, labels);
    }
  }
};

// Decay an on-shell parent using the physical phase-space API and HepMC vertices
std::array<HepMC3::GenParticlePtr, 2> Decay(HepMC3::GenEvent &event,
                                         const HepMC3::GenParticlePtr &parent,
                                         const std::array<int, 2> &pdg,
                                         gra::MPDG &table, gra::MRandom &random) {
  const auto momentum = gra::aux::HepMC2M4Vec(parent->momentum());
  std::vector<gra::M4Vec> daughters(2);
  REQUIRE(gra::kinematics::TwoBodyPhaseSpace(
              momentum, momentum.M(), {table.FindByPDG(pdg[0]).mass, table.FindByPDG(pdg[1]).mass},
              daughters, random).GetW() > 0.0);
  parent->set_status(gra::PDG::PDG_INTERMEDIATE);
  auto vertex = std::make_shared<HepMC3::GenVertex>();
  vertex->add_particle_in(parent);
  std::array<HepMC3::GenParticlePtr, 2> output;
  for (const auto &i : indices(output)) {
    output[i] = std::make_shared<HepMC3::GenParticle>(
        gra::aux::M4Vec2HepMC3(daughters[i]), pdg[i], gra::PDG::PDG_STABLE);
    vertex->add_particle_out(output[i]);
  }
  event.add_vertex(vertex);
  return output;
}

}  // namespace

// Recognize a physical central resonance produced by the two exchanged momenta
TEST_CASE("Oracle analyzer accepts physical resonance PDG records", "[gra::analyzer][hepmc]") {
  auto event = analysis_test::Event();
  const auto resonance = event.particles().front();
  resonance->set_pid(113);
  const auto momentum = gra::aux::HepMC2M4Vec(resonance->momentum());
  auto production = std::make_shared<HepMC3::GenVertex>();
  for (int i = 0; i < 2; ++i) {
    production->add_particle_in(std::make_shared<HepMC3::GenParticle>(
        gra::aux::M4Vec2HepMC3(momentum * 0.5), gra::PDG::PDG_propagator, gra::PDG::PDG_INTERMEDIATE));
  }
  production->add_particle_out(resonance);
  event.add_vertex(production);
  analysis_test::Write("oracle_resonance.hepmc3", {event});
  gra::MAnalyzer analyzer("resonance");
  AnalyzerHistograms histograms;
  CHECK(analyzer.HepMC3_OracleFill("../tmp/test_analysis/oracle_resonance", 2, 211, 10,
                                  histograms.h1, histograms.h2, histograms.hP, 0) == Approx(1.0));
  CHECK(histograms.h1.at("h1_S_M")->h[0]->GetMean() == Approx(momentum.M()));
}

// Count stable N* descendants through baryon and meson cascades
TEST_CASE("Oracle analyzer includes stable forward cascade products", "[gra::analyzer][hepmc][physics]") {
  gra::MPDG table;
  table.ReadParticleData();
  gra::MRandom random;
  random.SetSeed(314159);
  auto event = analysis_test::Event();
  const auto central = event.particles().front();
  gra::M4Vec nstar;
  nstar.SetPxPyPzM(0.0, 0.0, 3.0, 3.0);
  auto excited = std::make_shared<HepMC3::GenParticle>(
      gra::aux::M4Vec2HepMC3(nstar), gra::PDG::PDG_NSTAR, gra::PDG::PDG_INTERMEDIATE);
  gra::M4Vec proton;
  proton.SetPxPyPzM(-central->momentum().px(), 0.0, -3.0, gra::PDG::mp);
  const auto total = nstar + proton + gra::aux::HepMC2M4Vec(central->momentum());
  auto scattering = std::make_shared<HepMC3::GenVertex>();
  for (const double direction : {-1.0, 1.0}) {
    const double energy = total.E() / 2.0;
    scattering->add_particle_in(std::make_shared<HepMC3::GenParticle>(
        HepMC3::FourVector(0.0, 0.0, direction * std::sqrt(energy * energy - gra::PDG::mp * gra::PDG::mp), energy),
        gra::PDG::PDG_p, gra::PDG::PDG_BEAM));
  }
  scattering->add_particle_out(central);
  scattering->add_particle_out(excited);
  scattering->add_particle_out(std::make_shared<HepMC3::GenParticle>(
      gra::aux::M4Vec2HepMC3(proton), gra::PDG::PDG_p, gra::PDG::PDG_STABLE));
  event.add_vertex(scattering);
  const auto resonances = Decay(event, excited, {2114, 213}, table, random);
  const auto baryon = Decay(event, resonances[0], {2112, 111}, table, random);
  const auto meson = Decay(event, resonances[1], {211, 111}, table, random);
  Decay(event, baryon[1], {22, 22}, table, random);
  Decay(event, meson[1], {22, 22}, table, random);
  const double neutron_xf = 2.0 * baryon[0]->momentum().pz() / total.M();
  double pion_count = 1.0;
  SECTION("stable pion") {}
  SECTION("asymmetric beams") {
    const gra::M4Vec boost(0.0, 0.0, std::sinh(2.0), std::cosh(2.0));
    for (const auto &particle : event.particles()) {
      auto momentum = gra::aux::HepMC2M4Vec(particle->momentum());
      gra::kinematics::LorentzBoost(boost, 1.0, momentum, 1);
      particle->set_momentum(gra::aux::M4Vec2HepMC3(momentum));
    }
  }
  SECTION("antiproton excitation") {
    excited->set_pid(-gra::PDG::PDG_NSTAR);
    for (const auto &particle : HepMC3::Relatives::DESCENDANTS(excited)) {
      if (table.PDG_table.contains(-particle->pid())) { particle->set_pid(-particle->pid()); }
    }
    for (const auto &particle : event.particles()) {
      if (particle->status() == gra::PDG::PDG_BEAM && particle->momentum().pz() > 0.0) {
        particle->set_pid(-gra::PDG::PDG_p);
      }
    }
  }
  SECTION("decayed pion") {
    Decay(event, meson[0], {-13, 14}, table, random);
    pion_count = 0.0;
  }
  analysis_test::Write("oracle_cascade.hepmc3", {event});
  gra::MAnalyzer analyzer("cascade");
  AnalyzerHistograms histograms;
  REQUIRE(analyzer.HepMC3_OracleFill("../tmp/test_analysis/oracle_cascade", 2, 211, 10,
                                    histograms.h1, histograms.h2, histograms.hP, 0) == Approx(1.0));
  CHECK(analyzer.hE_Pions->Integral() == Approx(pion_count));
  CHECK(analyzer.hE_Neutron->Integral() == Approx(1.0));
  CHECK(analyzer.hXF_Neutron->GetMean() == Approx(neutron_xf));
  CHECK(analyzer.hE_Gamma->Integral() == Approx(4.0));
  CHECK(analyzer.hM_NSTAR->GetMean() == Approx(nstar.M()));
}

// Exchange the CS beam directions and check the corresponding angular reflection
TEST_CASE("Analyzer Collins-Soper angles respect beam exchange", "[gra::analyzer][frames][physics]") {
  gra::M4Vec plus, minus, pip, pim;
  plus.SetPxPyPzM(0.0, 0.0, 5.0, gra::PDG::mp);
  minus.SetPxPyPzM(0.0, 0.0, -5.0, gra::PDG::mp);
  pip.SetPxPyPzM(0.4, 0.1, 0.2, gra::PDG::mpi);
  pim.SetPxPyPzM(-0.3, -0.1, -0.2, gra::PDG::mpi);
  gra::MAnalyzer first("cs_first");
  gra::MAnalyzer exchanged("cs_exchanged");
  first.FrameObservables(1.0, plus, minus, {}, {}, {pip}, {pim});
  exchanged.FrameObservables(1.0, minus, plus, {}, {}, {pip}, {pim});
  const auto frames = gra::analyzer::Frames();
  const auto cs = std::distance(frames.begin(), std::find(frames.begin(), frames.end(), "CS"));
  CHECK(exchanged.h2CosTheta[cs][cs]->GetMean(1) == Approx(-first.h2CosTheta[cs][cs]->GetMean(1)));
  CHECK(exchanged.h2Phi[cs][cs]->GetMean(1) == Approx(-first.h2Phi[cs][cs]->GetMean(1)));
}

// Select neutral antiparticles independently of their charge and allow one neutral daughter
TEST_CASE("Oracle analyzer accepts neutral antiparticles and single photons", "[gra::analyzer][hepmc]") {
  const std::vector<std::pair<int, bool>> states = {{22, false}, {2112, false}, {2112, true}};
  for (const auto &[pdg, pair] : states) {
    const auto event = analysis_test::Event(pdg, pair);
    analysis_test::Write("oracle_neutral.hepmc3", {event});
    gra::MAnalyzer analyzer("neutral");
    AnalyzerHistograms histograms;
    const double normalization = analyzer.HepMC3_OracleFill(
        "../tmp/test_analysis/oracle_neutral", pair ? 2 : 1, pdg, 10,
        histograms.h1, histograms.h2, histograms.hP, 0);
    CHECK(normalization == Approx(1.0));
    CHECK(histograms.h1.at("h1_S_M")->h[0]->GetEntries() == Approx(1.0));
    CHECK(histograms.h1.at("h1_1B_pt")->h[0]->GetMean() == Approx(std::sqrt(0.17)));
  }
}

// Check unit conversion and fatal input parsing through the Oracle API
TEST_CASE("Oracle analyzer converts MeV input and rejects a broken second event", "[gra::analyzer][hepmc]") {
  analysis_test::Write("oracle_mev.hepmc3", {analysis_test::Event()}, HepMC3::Units::MEV);
  gra::MAnalyzer analyzer("units");
  AnalyzerHistograms histograms;
  CHECK(analyzer.HepMC3_OracleFill("../tmp/test_analysis/oracle_mev", 2, 211, 10,
                                  histograms.h1, histograms.h2, histograms.hP, 0) == Approx(1.0));
  CHECK(histograms.h1.at("h1_1B_pt")->h[0]->GetMean() == Approx(std::sqrt(0.17)));
  analysis_test::Broken("oracle_broken.hepmc3", "E 1 1 2\nP broken\n");
  CHECK_THROWS_AS(analyzer.HepMC3_OracleFill("../tmp/test_analysis/oracle_broken", 2, 211, 10,
                                           histograms.h1, histograms.h2, histograms.hP, 0), std::invalid_argument);
}

// Reject absent histograms and invalid sample slots before any histogram is filled
TEST_CASE("Oracle analyzer validates histogram inputs before filling", "[gra::analyzer][validation]") {
  analysis_test::Write("oracle_validation.hepmc3", {analysis_test::Event()});
  gra::MAnalyzer analyzer("validation");
  AnalyzerHistograms histograms;
  unsigned int sample = 0;
  SECTION("missing key") { histograms.h1.erase("h1_phi_GJ"); }
  SECTION("null multiplet") { histograms.h2.at("h2_S_M_t").reset(); }
  SECTION("invalid sample") { sample = 1; }
  SECTION("null histogram") {
    delete histograms.hP.at("hP_S_M_Pt")->h[0];
    histograms.hP.at("hP_S_M_Pt")->h[0] = nullptr;
  }
  const auto size1 = histograms.h1.size();
  const auto size2 = histograms.h2.size();
  CHECK_THROWS_AS(analyzer.HepMC3_OracleFill("../tmp/test_analysis/oracle_validation", 2, 211, 10,
                                           histograms.h1, histograms.h2, histograms.hP, sample), std::invalid_argument);
  CHECK(histograms.h1.size() == size1);
  CHECK(histograms.h2.size() == size2);
  CHECK(histograms.h1.at("h1_S_M")->h[0]->GetEntries() == Approx(0.0));
}

// Check full Legendre acceptance through the program's actual profile initializer
TEST_CASE("Analyzer profiles preserve isotropic Legendre moments", "[gra::analyzer][moments]") {
  std::map<std::string, std::shared_ptr<gra::hProfMultiplet>> profiles;
  const gra::h1Bound mass(1, 0.0, 2.0);
  gra::program::InitPrHistogram(profiles, {"isotropic"}, {2}, "", mass);
  auto *pt = profiles.at("hP_S_M_Pt")->h[0];
  pt->Fill(1.0, 0.0, 1.0);
  pt->Fill(1.0, 10.0, 3.0);
  CHECK(pt->GetBinEntries(1) == Approx(4.0));
  CHECK(pt->GetBinContent(1) == Approx(7.5));
  const int points = 10000;
  for (int i = 0; i < points; ++i) {
    const double costheta = -1.0 + (2.0 * i + 1.0) / points;
    profiles.at("hP_S_M_PL2_CM")->h[0]->Fill(1.0, gra::math::LegendrePl(2, costheta));
    profiles.at("hP_S_M_PL4_CM")->h[0]->Fill(1.0, gra::math::LegendrePl(4, costheta));
  }
  for (const auto &name : {"hP_S_M_PL2_CM", "hP_S_M_PL4_CM"}) {
    const auto *profile = profiles.at(name)->h[0];
    CHECK(profile->GetBinEntries(1) == Approx(points));
    CHECK(std::abs(profile->GetBinContent(1)) < 1e-7);
  }
  // A polar decay must retain P2 and P4 equal to one
  profiles.at("hP_S_M_PL2_CM")->h[0]->Fill(1.0, 1.0);
  profiles.at("hP_S_M_PL4_CM")->h[0]->Fill(1.0, 1.0);
  CHECK(profiles.at("hP_S_M_PL2_CM")->h[0]->GetBinEntries(1) == Approx(points + 1));
  CHECK(profiles.at("hP_S_M_PL4_CM")->h[0]->GetBinEntries(1) == Approx(points + 1));
}
