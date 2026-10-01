// Hard Pomeron PDF tests
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include <algorithm>
#include <array>
#include <atomic>
#include <catch.hpp>
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "Graniitti/Amplitude/MG5/Parton/AMP_MG5_PartonRegistry.h"
#include "Graniitti/Amplitude/MG5/Parton/AMP_MG5_pp_jj.h"
#include "Graniitti/Amplitude/MG5/Parton/AMP_MG5_pp_w.h"
#include "Graniitti/Amplitude/MG5/Parton/AMP_MG5_pp_z.h"
#include "Graniitti/Amplitude/MG5/Parton/AMP_MG5_pp_zj.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_JJ/Processes.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_W/Processes.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_Z/Processes.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_ZJ/Processes.h"
#include "Graniitti/Amplitude/MG5/Photon/AMP_MG5_yy_jj.h"
#include "Graniitti/Amplitude/MG5/Photon/AMP_MG5_yy_zjj.h"
#include "Graniitti/Amplitude/MG5/Photon/MG5_YY_JJ/Processes.h"
#include "Graniitti/Amplitude/MG5/Photon/MG5_YY_ZJJ/Processes.h"
#include "Graniitti/Kinematics/MCollinear.h"
#include "Graniitti/Kinematics/MEventKinematics.h"
#include "Graniitti/Kinematics/MHardDiffraction.h"
#include "Graniitti/Kinematics/MQuasiElastic.h"
#include "Graniitti/MGlobals.h"
#include "Graniitti/Math/MFloat.h"
#include "Graniitti/PDF/MHardPomeronPDF.h"
#include "Graniitti/Particle/MForm.h"
#include "Graniitti/Particle/MPDG.h"
#include "Graniitti/Photon/MFlux.h"
#include "Graniitti/Process/MSubProc.h"
#include "Graniitti/Spin/MHelicity.h"
#include "Graniitti/Tech/MLHE.h"
#include "HepMC3/Attribute.h"
#include "HepMC3/GenCrossSection.h"
#include "HepMC3/GenPdfInfo.h"
#include "HepMC3/GenVertex.h"
#include "json.hpp"
#include "support/models_test_support.hh"

using gra::aux::indices;

// Find one generated matrix element through its physical channel metadata
template <class ProcessBase>
ProcessBase *FindChannelProcess(std::vector<gra::mg5::Subprocess<ProcessBase>> &subprocesses,
                                const std::array<int, 2> &initial, const std::vector<int> &final) {
  for (auto &subprocess : subprocesses) {
    const auto channel = std::find_if(
        subprocess.channels.begin(), subprocess.channels.end(),
        [&](const gra::mg5::Channel &candidate) { return candidate.initial == initial && candidate.final == final; });
    if (channel != subprocess.channels.end()) { return subprocess.process.get(); }
  }
  return nullptr;
}

// Read the external pole mass from the same SLHA card as the generated family
double HardMass(int pdg, const std::string &family = "Parton/MG5_PP_Z") {
  const SLHAReader card(gra::aux::ResolveProjectPath("MG5cards/" + family + "/param_card.dat"));
  return card.get_block_entry("mass", std::abs(pdg));
}

// Close a rest-frame hard event with two massless incoming partons
void CloseHardRest(gra::LORENTZSCALAR &lts) {
  gra::M4Vec total;
  for (const auto &branch : lts.decaytree) { total += branch.p4; }
  REQUIRE(total.P3mod() < 1e-12 * total.E());
  const double energy = total.E() / 2.0;
  lts.q1 = gra::M4Vec(0.0, 0.0, energy, energy);
  lts.q2 = gra::M4Vec(0.0, 0.0, -energy, energy);
}

// Initialize the massless muon limit used by the independent analytic references
void MasslessMuon(gra::MG5Process &amplitude) {
  const auto process = amplitude.Processes().front();
  const auto path = gra::amplitude::ParameterCard(process.process_family, process.process_name);
  REQUIRE(path.has_value());
  SLHAReader card(gra::aux::ResolveProjectPath(*path));
  card.set_block_entry("mass", 13, 0.0);
  card.set_block_entry("yukawa", 13, 0.0);
  amplitude.InitParameters(std::move(card));
}

// Build one stable decay branch with fixed momentum
gra::MDecayBranch StableBranch(int pdg, const gra::M4Vec &p4, double mass = 0.0) {
  gra::MDecayBranch branch;
  branch.p.pdg = pdg;
  branch.p.mass = mass;
  branch.p4.SetPxPyPzM(p4.Px(), p4.Py(), p4.Pz(), mass);
  return branch;
}

// Rotate one generated decay branch and all stable daughters around the beam
void RotateGeneratedBranchAroundZ(gra::MDecayBranch &branch, double angle) {
  branch.p4.RotateZ(angle);
  for (auto &daughter : branch.legs) { RotateGeneratedBranchAroundZ(daughter, angle); }
}

// Rotate one generated hard event without changing its ordered beam labels
gra::LORENTZSCALAR RotateGeneratedEventAroundZ(gra::LORENTZSCALAR lts, double angle) {
  lts.q1.RotateZ(angle);
  lts.q2.RotateZ(angle);
  for (auto &branch : lts.decaytree) { RotateGeneratedBranchAroundZ(branch, angle); }
  return lts;
}

// Apply one three axis rotation used by hard frame covariance tests
void RotateCovarianceMomentum(gra::M4Vec &momentum, const std::array<double, 3> &angle = {0.41, -0.72, 1.19}) {
  momentum.RotateX(angle[0]);
  momentum.RotateY(angle[1]);
  momentum.RotateZ(angle[2]);
}

// Rotate one generated decay branch and every stable daughter
void RotateCovarianceBranch(gra::MDecayBranch &branch, const std::array<double, 3> &angle = {0.41, -0.72, 1.19}) {
  RotateCovarianceMomentum(branch.p4, angle);
  for (auto &daughter : branch.legs) { RotateCovarianceBranch(daughter, angle); }
}

// Apply the fixed noncollinear boost used by hard frame covariance tests
void BoostCovarianceMomentum(gra::M4Vec &momentum) { momentum = momentum.LorentzBoost({0.21, -0.13, 0.31}); }

// Boost one generated decay branch and every stable daughter
void BoostCovarianceBranch(gra::MDecayBranch &branch) {
  BoostCovarianceMomentum(branch.p4);
  for (auto &daughter : branch.legs) { BoostCovarianceBranch(daughter); }
}

// Transport one generated decay branch between equal-mass hard systems
void TransportGeneratedBranch(gra::MDecayBranch &branch, const gra::M4Vec &source, const gra::M4Vec &target) {
  gra::kinematics::LorentzBoost(source, source.M(), branch.p4, -1);
  gra::kinematics::LorentzBoost(target, target.M(), branch.p4, 1);
  for (auto &daughter : branch.legs) { TransportGeneratedBranch(daughter, source, target); }
}

// Rotate one complete generated hard event
gra::LORENTZSCALAR RotateCovarianceEvent(gra::LORENTZSCALAR           lts,
                                         const std::array<double, 3> &angle = {0.41, -0.72, 1.19}) {
  RotateCovarianceMomentum(lts.q1, angle);
  RotateCovarianceMomentum(lts.q2, angle);
  for (auto &branch : lts.decaytree) { RotateCovarianceBranch(branch, angle); }
  return lts;
}

// Boost one complete generated hard event
gra::LORENTZSCALAR BoostCovarianceEvent(gra::LORENTZSCALAR lts) {
  BoostCovarianceMomentum(lts.q1);
  BoostCovarianceMomentum(lts.q2);
  for (auto &branch : lts.decaytree) { BoostCovarianceBranch(branch); }
  return lts;
}

// Boost one generated decay branch and every stable daughter in rapidity
void BoostRapidityBranch(gra::MDecayBranch &branch, double rapidity) {
  branch.p4 = branch.p4.LorentzBoost({0.0, 0.0, std::tanh(rapidity)});
  for (auto &daughter : branch.legs) { BoostRapidityBranch(daughter, rapidity); }
}

// Boost one complete generated hard event in rapidity
gra::LORENTZSCALAR BoostRapidityEvent(gra::LORENTZSCALAR lts, double rapidity) {
  const gra::M3Vec boost = {0.0, 0.0, std::tanh(rapidity)};
  lts.q1                 = lts.q1.LorentzBoost(boost);
  lts.q2                 = lts.q2.LorentzBoost(boost);
  for (auto &branch : lts.decaytree) { BoostRapidityBranch(branch, rapidity); }
  return lts;
}

// Bind the immutable TUNE0 cache used by direct generated amplitude fixtures
void BindHardModelCache(gra::LORENTZSCALAR &lts) {
  if (lts.model_cache == nullptr) {
    lts.model_cache = std::make_shared<gra::MModelCache>(gra::MModelTune::Load(modelfile));
  }
}

// Evaluate one generated hard family through its result interface
template <class Amplitude>
gra::PartonMG5Evaluation EvaluateHard(Amplitude &amplitude, gra::LORENTZSCALAR &lts, double alpha_s) {
  BindHardModelCache(lts);
  return amplitude.EvaluatePrepared(lts, alpha_s);
}

// Compute one generated hard family matrix element squared
template <class Amplitude>
double HardAmp2(Amplitude &amplitude, gra::LORENTZSCALAR &lts, double alpha_s) {
  return EvaluateHard(amplitude, lts, alpha_s).amp2;
}

// Compute the fixed width LO neutral current partonic cross section
double NeutralCurrentSigmaLO(double mass, int pid, bool z_only = false) {
  constexpr double alpha        = 1.0 / 137.03599908;
  constexpr double sin2         = 0.22224648578577766;
  constexpr double z_mass       = 91.188;
  constexpr double z_width      = 2.441404;
  const double     mass2        = mass * mass;
  const double     theta        = 1.0 / (16.0 * sin2 * (1.0 - sin2));
  const double     pole         = mass2 - z_mass * z_mass;
  const double     denominator  = pole * pole + z_mass * z_mass * z_width * z_width;
  const double     gamma        = 4.0 * gra::math::PI * alpha * alpha / (3.0 * mass2);
  const double     interference = gamma * 2.0 * theta * mass2 * pole / denominator;
  const double     resonance    = gamma * gra::math::pow2(theta * mass2) / denominator;
  const double     charge       = std::abs(pid) % 2 == 0 ? 2.0 / 3.0 : -1.0 / 3.0;
  const double     axial        = std::abs(pid) % 2 == 0 ? 1.0 : -1.0;
  const double     vector       = axial - 4.0 * sin2 * charge;
  constexpr double muon_charge  = -1.0;
  constexpr double muon_axial   = -1.0;
  const double     muon_vector  = muon_axial - 4.0 * sin2 * muon_charge;
  return ((z_only ? 0.0 : charge * charge * gamma * muon_charge * muon_charge +
                          charge * vector * interference * muon_charge * muon_vector) +
          (vector * vector + axial * axial) * resonance * (muon_vector * muon_vector + muon_axial * muon_axial)) /
         3.0;
}

// Build one massless neutral current partonic point at fixed decay angle
gra::LORENTZSCALAR NeutralCurrentPoint(double mass, int pid, double cosine, bool z_only = false) {
  const double       momentum = 0.5 * mass;
  const double       sine     = std::sqrt(1.0 - cosine * cosine);
  gra::LORENTZSCALAR lts;
  lts.q1                      = gra::M4Vec(0.0, 0.0, momentum, momentum);
  lts.q2                      = gra::M4Vec(0.0, 0.0, -momentum, momentum);
  lts.id1                     = pid;
  lts.id2                     = -pid;
  lts.process.root_decay_mode = gra::RootDecayMode::Physical;
  lts.decaytree               = {StableBranch(-13, gra::M4Vec(momentum * sine, 0.0, momentum * cosine, momentum)),
                                 StableBranch(13, gra::M4Vec(-momentum * sine, 0.0, -momentum * cosine, momentum))};
  if (z_only) {
    gra::MDecayBranch z;
    z.p.pdg = 23;
    z.p4 = lts.q1 + lts.q2;
    z.legs = std::move(lts.decaytree);
    lts.decaytree = {std::move(z)};
  }
  return lts;
}

// Compute one generated photon Z plus two parton matrix element squared
double PhotonZjjAmp2(gra::AMP_MG5_yy_zjj &amplitude, gra::LORENTZSCALAR &lts, bool coherent_epa) {
  BindHardModelCache(lts);
  const auto result = amplitude.Evaluate(lts, 0.0, coherent_epa);
  return result.Valid() ? result.amp2 : 0.0;
}

// Build one four-vector from light-cone components
gra::M4Vec LightConeTestVector(double plus, double minus, double px, double py) {
  return gra::M4Vec(px, py, 0.5 * (plus - minus), 0.5 * (plus + minus));
}

// Build one plus-side leading proton for a massless beam test point
gra::M4Vec LeadingProtonTestVector(double beam_energy, double xi, double t, double phi) {
  const double pt2       = -(1.0 - xi) * t;
  const double pt        = std::sqrt(pt2);
  const double lead_plus = (1.0 - xi) * 2.0 * beam_energy;
  return LightConeTestVector(lead_plus, pt2 / lead_plus, pt * std::cos(phi), pt * std::sin(phi));
}

// Build one spacelike hard parton with a complementary lightlike Pomeron
// remnant
gra::M4Vec DiffractiveHardPartonTestVector(const gra::M4Vec &exchange, double beta, bool plus_side = true) {
  if (gra::math::IsExactEqual(beta, 1.0)) { return exchange; }

  const double remnant_px    = (1.0 - beta) * exchange.Px();
  const double remnant_py    = (1.0 - beta) * exchange.Py();
  const double remnant_lc    = (1.0 - beta) * (plus_side ? exchange.LightconePos() : exchange.LightconeNeg());
  const double conjugate_lc  = (remnant_px * remnant_px + remnant_py * remnant_py) / remnant_lc;
  const double remnant_plus  = plus_side ? remnant_lc : conjugate_lc;
  const double remnant_minus = plus_side ? conjugate_lc : remnant_lc;
  return exchange - LightConeTestVector(remnant_plus, remnant_minus, remnant_px, remnant_py);
}

// Build one minus-side collinear hard parton for a massless beam test point
gra::M4Vec MinusCollinearHardPartonTestVector(double beam_energy, double x) {
  return LightConeTestVector(0.0, x * 2.0 * beam_energy, 0.0, 0.0);
}

// Compute one HepMC four-vector invariant mass squared
double HepMCMass2(const HepMC3::FourVector &p) {
  return p.e() * p.e() - p.px() * p.px() - p.py() * p.py() - p.pz() * p.pz();
}

// Convert one HepMC momentum to the internal four-vector convention
gra::M4Vec HepMCMomentum(const HepMC3::FourVector &p) { return gra::M4Vec(p.px(), p.py(), p.pz(), p.e()); }

// Write one temporary hard-Pomeron tune with a modified Q2 floor
std::string WriteHardPomeronTune(const std::string &suffix, double q2_min) {
  const std::filesystem::path dir = "tmp/graniitti_hardpomeron_" + suffix;
  std::filesystem::create_directories(dir);

  auto j = nlohmann::json::parse(gra::aux::GetInputData(gra::ResolveModelDataFile("TUNE0", "GENERAL.json")));
  j["PARAM_HARDPOMERON"]["Q2_min"] = q2_min;

  const std::filesystem::path general_path = dir / "GENERAL.json";
  std::ofstream               out(general_path);
  if (!out.good()) { throw std::runtime_error("WriteHardPomeronTune: failed to write GENERAL.json"); }
  out << j.dump(2);
  out.close();

  std::filesystem::copy_file(gra::ResolveModelDataFile("TUNE0", "NUMERICS.json"), dir / "NUMERICS.json",
                             std::filesystem::copy_options::overwrite_existing);
  return general_path.string();
}

// Build a Z cascade branch with two muon daughters
gra::MDecayBranch ZBranch(const gra::M4Vec &mup, const gra::M4Vec &mum, double mass = HardMass(13)) {
  gra::MDecayBranch branch;
  branch.p.pdg = 23;
  branch.legs.push_back(StableBranch(-13, mup, mass));
  branch.legs.push_back(StableBranch(13, mum, mass));
  branch.p4 = branch.legs[0].p4 + branch.legs[1].p4;
  return branch;
}

// Build the fixed gamma-gamma Z+2-parton benchmark point
gra::LORENTZSCALAR YYZjjReferenceKinematics(int quark_pdg = 2) {
  gra::LORENTZSCALAR lts;
  lts.q1 = gra::M4Vec(0.0, 0.0, 220.0, 220.0);
  lts.q2 = gra::M4Vec(0.0, 0.0, -220.0, 220.0);
  lts.decaytree.push_back(ZBranch(gra::M4Vec(100.0, 0.0, 0.0, 100.0), gra::M4Vec(-100.0, 0.0, 0.0, 100.0)));
  lts.decaytree.push_back(StableBranch(quark_pdg, gra::M4Vec(0.0, 120.0, 0.0, 120.0), HardMass(quark_pdg, "Photon/MG5_YY_ZJJ")));
  lts.decaytree.push_back(StableBranch(-quark_pdg, gra::M4Vec(0.0, -120.0, 0.0, 120.0), HardMass(quark_pdg, "Photon/MG5_YY_ZJJ")));
  CloseHardRest(lts);
  return lts;
}

// Build the fixed gamma-gamma to Zjj point with photon flux kinematics
gra::LORENTZSCALAR YYZjjFluxReferenceKinematics() {
  gra::LORENTZSCALAR lts           = YYZjjReferenceKinematics();
  const double       beam_energy   = 6500.0;
  const double       beam_momentum = std::sqrt(beam_energy * beam_energy - gra::math::pow2(gra::PDG::mp));
  lts.pbeam1                       = gra::M4Vec(0.0, 0.0, beam_momentum, beam_energy);
  lts.pbeam2                       = gra::M4Vec(0.0, 0.0, -beam_momentum, beam_energy);
  lts.beam1.pdg                    = gra::PDG::PDG_p;
  lts.beam1.chargeX3               = 3;
  lts.beam1.spinX2                 = 1;
  lts.beam1.mass                   = gra::PDG::mp;
  lts.beam2                        = lts.beam1;
  const double photon_energy = lts.q1.E();
  lts.q1 = gra::M4Vec(0.05, 0.0, photon_energy, photon_energy);
  lts.q2 = gra::M4Vec(-0.05, 0.0, -photon_energy, photon_energy);
  lts.process.root_decay_mode      = gra::RootDecayMode::Physical;
  lts.pfinal.assign(6, gra::M4Vec(0.0, 0.0, 0.0, 0.0));
  lts.pfinal[0]     = lts.q1 + lts.q2;
  lts.pfinal[1]     = lts.pbeam1 - lts.q1;
  lts.pfinal[2]     = lts.pbeam2 - lts.q2;
  lts.forward_mass2 = {lts.pfinal[1].M2(), lts.pfinal[2].M2()};
  lts.s             = (lts.pbeam1 + lts.pbeam2).M2();
  lts.x1            = gra::kinematics::LongitudinalMomentumLoss(lts.pbeam1, lts.pfinal[1], true);
  lts.x2            = gra::kinematics::LongitudinalMomentumLoss(lts.pbeam2, lts.pfinal[2], false);
  lts.xi1           = lts.x1;
  lts.xi2           = lts.x2;
  lts.has_xi1       = true;
  lts.has_xi2       = true;
  lts.s_hat         = (lts.q1 + lts.q2).M2();
  lts.m2            = lts.s_hat;
  lts.t1            = lts.q1.M2();
  lts.t2            = lts.q2.M2();
  lts.qt1           = lts.q1.Pt();
  lts.qt2           = lts.q2.Pt();
  lts.LHAPDFSET     = "LUXqed17_plus_PDF4LHC15_nnlo_100";
  return lts;
}

// Load particle data needed by process-level lazy constructors
void LoadReferencePDG(gra::LORENTZSCALAR &lts) {
  gra::MODELPARAM = "TUNE0";
  lts.PDG.ReadParticleData();
}

// Build the fixed u ubar to neutral-current dimuon benchmark point
gra::LORENTZSCALAR PPZReferenceKinematics(double alpha_s, double muon_mass = HardMass(13)) {
  gra::LORENTZSCALAR lts;
  lts.q1       = gra::M4Vec(0.0, 0.0, 100.0, 100.0);
  lts.q2       = gra::M4Vec(0.0, 0.0, -100.0, 100.0);
  lts.id1      = 2;
  lts.id2      = -2;
  lts.alphaQCD = alpha_s;
  lts.decaytree.push_back(StableBranch(-13, gra::M4Vec(80.0, 0.0, 60.0, 100.0), muon_mass));
  lts.decaytree.push_back(StableBranch(13, gra::M4Vec(-80.0, 0.0, -60.0, 100.0), muon_mass));
  CloseHardRest(lts);
  return lts;
}

// Build the fixed g u to Z u benchmark point
gra::LORENTZSCALAR PPZjReferenceKinematics(double alpha_s, double muon_mass = HardMass(13)) {
  gra::LORENTZSCALAR lts;
  lts.q1       = gra::M4Vec(0.0, 0.0, 160.0, 160.0);
  lts.q2       = gra::M4Vec(0.0, 0.0, -160.0, 160.0);
  lts.id1      = 21;
  lts.id2      = 2;
  lts.alphaQCD = alpha_s;
  lts.decaytree.push_back(ZBranch(gra::M4Vec(80.0, 0.0, 60.0, 100.0), gra::M4Vec(-80.0, 0.0, 60.0, 100.0), muon_mass));
  lts.decaytree.push_back(StableBranch(2, gra::M4Vec(0.0, 0.0, -120.0, 120.0)));
  CloseHardRest(lts);
  return lts;
}

// Build a non-collinear Z plus jet point valid for every generated crossing
gra::LORENTZSCALAR PPZjAllChannelKinematics(double alpha_s) {
  gra::LORENTZSCALAR lts;
  lts.q1          = gra::M4Vec(0.0, 0.0, 140.0, 140.0);
  lts.q2          = gra::M4Vec(0.0, 0.0, -140.0, 140.0);
  lts.id1         = 21;
  lts.id2         = 2;
  lts.alphaQCD    = alpha_s;
  const double py = std::sqrt(4800.0);
  lts.decaytree.push_back(ZBranch(gra::M4Vec(80.0, 0.0, 60.0, 100.0), gra::M4Vec(-40.0, py, -60.0, 100.0)));
  lts.decaytree.push_back(StableBranch(2, gra::M4Vec(-40.0, -py, 0.0, 80.0)));
  CloseHardRest(lts);
  return lts;
}

// Build the fixed g g to u ubar dijet benchmark point
gra::LORENTZSCALAR PPJJReferenceKinematics(double alpha_s) {
  gra::LORENTZSCALAR lts;
  lts.q1       = gra::M4Vec(0.0, 0.0, 160.0, 160.0);
  lts.q2       = gra::M4Vec(0.0, 0.0, -160.0, 160.0);
  lts.id1      = 21;
  lts.id2      = 21;
  lts.alphaQCD = alpha_s;
  lts.decaytree.push_back(StableBranch(2, gra::M4Vec(160.0, 0.0, 0.0, 160.0)));
  lts.decaytree.push_back(StableBranch(-2, gra::M4Vec(-160.0, 0.0, 0.0, 160.0)));
  return lts;
}

// Build the fixed g u to g u dijet benchmark point
gra::LORENTZSCALAR PPJJReversedFlowKinematics(double alpha_s, bool mirror_initial) {
  gra::LORENTZSCALAR lts = PPJJReferenceKinematics(alpha_s);
  lts.id1                = mirror_initial ? 2 : 21;
  lts.id2                = mirror_initial ? 21 : 2;
  lts.decaytree[0].p.pdg = 21;
  lts.decaytree[1].p.pdg = 2;
  return lts;
}

// Build the fixed charged-current muon-neutrino benchmark point
gra::LORENTZSCALAR PPWReferenceKinematics(double alpha_s) {
  gra::LORENTZSCALAR lts;
  lts.q1       = gra::M4Vec(0.0, 0.0, 100.0, 100.0);
  lts.q2       = gra::M4Vec(0.0, 0.0, -100.0, 100.0);
  lts.id1      = 2;
  lts.id2      = -1;
  lts.alphaQCD = alpha_s;
  lts.decaytree.push_back(StableBranch(-13, gra::M4Vec(80.0, 0.0, 60.0, 100.0), HardMass(13)));
  lts.decaytree.push_back(StableBranch(14, gra::M4Vec(-80.0, 0.0, -60.0, 100.0)));
  CloseHardRest(lts);
  return lts;
}

// Build the fixed gamma-gamma dijet benchmark point
gra::LORENTZSCALAR YYJJReferenceKinematics(int quark_pdg = 2) {
  gra::LORENTZSCALAR lts = PPJJReferenceKinematics(0.0);
  lts.id1                = gra::PDG::PDG_gamma;
  lts.id2                = gra::PDG::PDG_gamma;
  lts.decaytree[0].p.pdg = quark_pdg;
  lts.decaytree[1].p.pdg = -quark_pdg;
  return lts;
}

// Flatten one generated external color row for exact comparison
std::vector<int> FlattenExternalColorFlow(const gra::mg5helas::ExternalColorFlow &flow) {
  std::vector<int> flattened;
  flattened.reserve(2 * flow.size());
  for (const auto &leg : flow) {
    flattened.push_back(leg.color);
    flattened.push_back(leg.anticolor);
  }
  return flattened;
}

// Check that the scalar screening amplitude stores the matrix-element norm
void RequireScalarHamp(const gra::LORENTZSCALAR &lts, double amp2) {
  REQUIRE(lts.hamp.size() == 1);
  REQUIRE(std::norm(lts.hamp.front()) == Approx(amp2).epsilon(1e-10));
}

// Check that structured MG5 helicity components retain their complex norm
void RequireComplexHelicityHamp(const gra::LORENTZSCALAR &lts, double amp2, bool require_phase = true) {
  REQUIRE(lts.hamp.size() > 1);
  double norm_sum  = 0.0;
  bool   has_phase = false;
  for (const auto amplitude : lts.hamp) {
    norm_sum += std::norm(amplitude);
    if (std::abs(std::imag(amplitude)) > 1.0e-15) { has_phase = true; }
  }
  REQUIRE(norm_sum == Approx(amp2).epsilon(1e-10));
  if (require_phase) { REQUIRE(has_phase); }
}

// Check complex components stored before the two-incoming-spin average
void RequireUnaveragedHelicityHamp(const gra::LORENTZSCALAR &lts, double amp2) {
  REQUIRE(lts.hamp.size() > 1);
  double norm_sum = 0.0;
  for (const auto amplitude : lts.hamp) { norm_sum += std::norm(amplitude); }
  REQUIRE(norm_sum == Approx(4.0 * amp2).epsilon(1e-10));
}

// Compute the normalized Hermitian overlap of two complex helicity vectors
std::complex<double> NormalizedHelicityOverlap(const std::vector<std::complex<double>> &reference,
                                               const std::vector<std::complex<double>> &shifted, double reference_norm,
                                               double shifted_norm) {
  REQUIRE(reference.size() == shifted.size());
  std::complex<double> overlap = 0.0;
  for (std::size_t i = 0; i < reference.size(); ++i) { overlap += std::conj(reference[i]) * shifted[i]; }
  return overlap / std::sqrt(reference_norm * shifted_norm);
}

// Build one fixed radius screening point while preserving hard closure
gra::LORENTZSCALAR ScreeningRingPoint(const gra::LORENTZSCALAR &reference, double transverse_momentum, double phi) {
  gra::LORENTZSCALAR shifted = reference;
  const double       kx      = transverse_momentum * std::cos(phi);
  const double       ky      = transverse_momentum * std::sin(phi);
  shifted.q1                 = gra::M4Vec(kx, ky, reference.q1.Pz(), reference.q1.E());
  shifted.q2                 = gra::M4Vec(-kx, -ky, reference.q2.Pz(), reference.q2.E());
  return shifted;
}

// Check two four-momenta componentwise at one relative numerical tolerance
void RequireMomentumClose(const gra::M4Vec &got, const gra::M4Vec &expected, double tolerance,
                          const std::string &label) {
  const gra::M4Vec delta = got - expected;
  const double     scale = std::max(
          {1.0, std::abs(expected.E()), std::abs(expected.Px()), std::abs(expected.Py()), std::abs(expected.Pz())});
  CAPTURE(label, got, expected, tolerance, scale);
  REQUIRE(std::abs(delta.E()) <= tolerance * scale);
  REQUIRE(delta.P3mod() <= tolerance * scale);
}

// Build the direct covariant incoming pair and check its defining constraints
std::array<gra::M4Vec, 2> DirectProjectedPair(const gra::LORENTZSCALAR &lts, const std::string &label) {
  const auto                    leaves      = gra::mg5::StableDecayLeaves(lts.decaytree);
  std::vector<gra::M4Vec>       final       = gra::mg5::StableLeafMomenta(leaves);
  const std::vector<gra::M4Vec> input_final = final;
  gra::M4Vec                    p1;
  gra::M4Vec                    p2;
  INFO(label);
  REQUIRE(gra::mg5helas::PrepareOnShellKinematics(lts, final, p1, p2));

  for (const auto &i : indices(final)) {
    RequireMomentumClose(final[i], input_final[i], 1.0e-15, label + " final state");
  }
  gra::M4Vec hard;
  for (const auto &particle : final) { hard += particle; }
  const double mass2 = hard.M2();
  REQUIRE(mass2 > 0.0);
  REQUIRE(std::abs(p1.M2()) <= 1.0e-10 * mass2);
  REQUIRE(std::abs(p2.M2()) <= 1.0e-10 * mass2);
  RequireMomentumClose(p1 + p2, hard, 1.0e-10, label + " closure");
  return {p1, p2};
}

// Check covariance of the direct incoming projection under one transformation
void RequireDirectProjectionCovariance(const gra::LORENTZSCALAR &reference, const gra::LORENTZSCALAR &transformed,
                                       const std::function<void(gra::M4Vec &)> &transform, double tolerance,
                                       const std::string &label) {
  const auto reference_pair   = DirectProjectedPair(reference, label);
  const auto transformed_pair = DirectProjectedPair(transformed, label);
  gra::M4Vec expected1        = reference_pair[0];
  gra::M4Vec expected2        = reference_pair[1];
  transform(expected1);
  transform(expected2);
  RequireMomentumClose(transformed_pair[0], expected1, tolerance, label + " incoming leg 1");
  RequireMomentumClose(transformed_pair[1], expected2, tolerance, label + " incoming leg 2");
}

// Embed one spacelike diffractive leg and one lightlike collinear leg
void SetHardSDPair(gra::LORENTZSCALAR &lts, bool first_diffractive) {
  constexpr double virtuality = 0.12;
  constexpr double qx         = 0.17;
  constexpr double qy         = -0.11;
  gra::M4Vec       hard;
  for (const auto &leaf : gra::mg5::StableDecayLeaves(lts.decaytree)) { hard += leaf->p4; }
  const double mass = hard.M();
  REQUIRE(mass > 0.0);
  REQUIRE(hard.P3mod() <= 1.0e-14 * mass);

  const double     momentum         = (hard.M2() + virtuality) / (2.0 * mass);
  const double     spacelike_energy = mass - momentum;
  const double     pz               = std::sqrt(momentum * momentum - qx * qx - qy * qy);
  const gra::M4Vec plus_lightlike(qx, qy, pz, momentum);
  const gra::M4Vec plus_spacelike(qx, qy, pz, spacelike_energy);
  const gra::M4Vec minus_lightlike(-qx, -qy, -pz, momentum);
  const gra::M4Vec minus_spacelike(-qx, -qy, -pz, spacelike_energy);
  lts.q1 = first_diffractive ? plus_spacelike : plus_lightlike;
  lts.q2 = first_diffractive ? minus_lightlike : minus_spacelike;

  REQUIRE(std::abs((lts.q1 + lts.q2).M2() - hard.M2()) <= 1.0e-13 * hard.M2());
  const gra::M4Vec &diffractive = first_diffractive ? lts.q1 : lts.q2;
  const gra::M4Vec &collinear   = first_diffractive ? lts.q2 : lts.q1;
  REQUIRE(std::abs(diffractive.M2() + virtuality) <= 1.0e-12 * hard.M2());
  REQUIRE(std::abs(collinear.M2()) <= 1.0e-12 * hard.M2());
}

// Check direct scalar covariance for one generated hard family
template <class Amplitude>
void RequireHardScalarCovariance(Amplitude &amplitude, gra::LORENTZSCALAR reference, double alpha_s,
                                 const std::string &label) {
  constexpr double rotation_tolerance            = 2.0e-9;
  constexpr double boost_tolerance               = 2.0e-8;
  constexpr double rapidity_projection_tolerance = 2.0e-11;
  // Direct double-precision HELAS loses a few ppm at absolute rapidity five
  constexpr double rapidity_tolerance = 5.0e-6;
  reference.process.root_decay_mode   = gra::RootDecayMode::Physical;

  const double reference_norm = HardAmp2(amplitude, reference, alpha_s);
  REQUIRE(reference_norm > 0.0);

  const std::array<std::array<double, 3>, 3> rotations = {std::array<double, 3>{0.41, -0.72, 1.19},
                                                          std::array<double, 3>{-1.07, 0.33, 2.14},
                                                          std::array<double, 3>{2.37, -1.11, -0.58}};
  for (const auto &angle : rotations) {
    gra::LORENTZSCALAR rotated = RotateCovarianceEvent(reference, angle);
    RequireDirectProjectionCovariance(
        reference, rotated, [angle](gra::M4Vec &momentum) { RotateCovarianceMomentum(momentum, angle); },
        rotation_tolerance, label + " rotated projection");
    const double rotated_norm = HardAmp2(amplitude, rotated, alpha_s);
    CAPTURE(label, angle[0], angle[1], angle[2]);
    REQUIRE(rotated_norm == Approx(reference_norm).epsilon(rotation_tolerance));
  }

  const auto         angle   = rotations[0];
  gra::LORENTZSCALAR rotated = RotateCovarianceEvent(reference, angle);
  gra::LORENTZSCALAR boosted = BoostCovarianceEvent(rotated);
  RequireDirectProjectionCovariance(
      rotated, boosted, [](gra::M4Vec &momentum) { BoostCovarianceMomentum(momentum); }, boost_tolerance,
      label + " noncollinear boost projection");
  const double boosted_norm = HardAmp2(amplitude, boosted, alpha_s);
  INFO(label << " rotated and noncollinearly boosted");
  REQUIRE(boosted_norm == Approx(reference_norm).epsilon(boost_tolerance));

  for (const double rapidity : {-5.0, 5.0}) {
    gra::LORENTZSCALAR rapidity_boosted = BoostRapidityEvent(reference, rapidity);
    RequireDirectProjectionCovariance(
        reference, rapidity_boosted,
        [rapidity](gra::M4Vec &momentum) {
          momentum = momentum.LorentzBoost({0.0, 0.0, std::tanh(rapidity)});
        },
        rapidity_projection_tolerance, label + " rapidity projection");
    const double rapidity_norm     = HardAmp2(amplitude, rapidity_boosted, alpha_s);
    const double rapidity_residual = std::abs(rapidity_norm - reference_norm) / reference_norm;
    CAPTURE(label, rapidity, rapidity_boosted.q1.E(), rapidity_boosted.q2.E());
    CAPTURE(reference_norm, rapidity_norm, rapidity_residual);
    REQUIRE(rapidity_residual <= rapidity_tolerance);
  }
}

// Check one screening ring in a common beam axis helicity section
void RequireRotatedScreeningRingSection(gra::LORENTZSCALAR                                 reference,
                                        const std::function<double(gra::LORENTZSCALAR &)> &evaluate,
                                        const std::string                                 &label) {
  constexpr int    angular_points      = 16;
  constexpr double transverse_momentum = 0.05;
  constexpr double rotation            = 0.73;
  // Resolve phases an order below the existing screening continuity bound
  constexpr double    section_tolerance = 2.0e-9;
  const double        pi                = std::acos(-1.0);
  const double        seam_offset       = 1.0e-7;
  std::vector<double> angles;
  angles.reserve(angular_points + 2);
  for (int j = 0; j < angular_points; ++j) { angles.push_back(2.0 * pi * (j + 0.5) / angular_points); }
  angles.push_back(pi - seam_offset);
  angles.push_back(-pi + seam_offset);

  reference.process.root_decay_mode = gra::RootDecayMode::Physical;
  std::vector<std::complex<double>> phases;
  std::vector<bool>                 phase_set;
  std::size_t                       phase_count = 0;

  for (const double phi : angles) {
    gra::LORENTZSCALAR point      = ScreeningRingPoint(reference, transverse_momentum, phi);
    const double       point_norm = evaluate(point);
    REQUIRE(point_norm > 0.0);

    gra::LORENTZSCALAR rotated      = RotateGeneratedEventAroundZ(point, rotation);
    const double       rotated_norm = evaluate(rotated);
    INFO(label << " phi=" << phi);
    REQUIRE(rotated_norm == Approx(point_norm).epsilon(5.0e-10));
    REQUIRE(rotated.hamp.size() == point.hamp.size());

    if (phases.empty()) {
      phases.assign(point.hamp.size(), 0.0);
      phase_set.assign(point.hamp.size(), false);
    }
    const double component_scale = std::sqrt(point_norm / static_cast<double>(point.hamp.size()));
    for (std::size_t index = 0; index < point.hamp.size(); ++index) {
      if (std::abs(point.hamp[index]) <= 1.0e-8 * component_scale) { continue; }
      if (!phase_set[index]) {
        phases[index]    = rotated.hamp[index] / point.hamp[index];
        phase_set[index] = true;
        ++phase_count;
      }
      const std::complex<double> expected  = phases[index] * point.hamp[index];
      const double               tolerance = section_tolerance * std::max(component_scale, std::abs(expected));
      CAPTURE(phi, index, phases[index], point.hamp[index], rotated.hamp[index], expected);
      REQUIRE(std::abs(rotated.hamp[index] - expected) < tolerance);
    }
  }
  REQUIRE(phase_count > 0);

  gra::LORENTZSCALAR seam_upper = ScreeningRingPoint(reference, transverse_momentum, pi - seam_offset);
  gra::LORENTZSCALAR seam_lower = ScreeningRingPoint(reference, transverse_momentum, -pi + seam_offset);
  const double       upper_norm = evaluate(seam_upper);
  const double       lower_norm = evaluate(seam_lower);
  INFO(label << " atan2 seam");
  REQUIRE(upper_norm > 0.0);
  REQUIRE(lower_norm > 0.0);
  REQUIRE(lower_norm == Approx(upper_norm).epsilon(section_tolerance));
  REQUIRE(seam_lower.hamp.size() == seam_upper.hamp.size());
  const double seam_scale = std::sqrt(upper_norm / static_cast<double>(seam_upper.hamp.size()));
  for (std::size_t index = 0; index < seam_upper.hamp.size(); ++index) {
    const double tolerance = section_tolerance * std::max(seam_scale, std::abs(seam_upper.hamp[index]));
    CAPTURE(index, seam_upper.hamp[index], seam_lower.hamp[index]);
    REQUIRE(std::abs(seam_lower.hamp[index] - seam_upper.hamp[index]) < tolerance);
  }
}

// Check azimuthal continuity of one fixed-basis MG5 screening vector
void RequireHelicitySectionContinuity(gra::LORENTZSCALAR                                 reference,
                                      const std::function<double(gra::LORENTZSCALAR &)> &evaluate,
                                      const std::string                                 &label) {
  constexpr int      angular_points      = 16;
  const double       transverse_momentum = 0.05;
  gra::LORENTZSCALAR projected_reference = reference;
  INFO(label);
  const double reference_norm       = evaluate(projected_reference);
  const auto   reference_amplitudes = projected_reference.hamp;
  REQUIRE_FALSE(reference_amplitudes.empty());
  REQUIRE(reference_norm > 0.0);
  std::vector<std::complex<double>> angular_average(reference_amplitudes.size(), 0.0);
  std::complex<double>              angular_overlap = 0.0;

  for (int j = 0; j < angular_points; ++j) {
    const double       phi          = 2.0 * std::acos(-1.0) * (j + 0.5) / angular_points;
    gra::LORENTZSCALAR shifted      = ScreeningRingPoint(projected_reference, transverse_momentum, phi);
    const double       shifted_norm = evaluate(shifted);
    REQUIRE(shifted.hamp.size() == angular_average.size());
    for (std::size_t index = 0; index < angular_average.size(); ++index) {
      angular_average[index] += shifted.hamp[index] / static_cast<double>(angular_points);
    }
    angular_overlap += NormalizedHelicityOverlap(reference_amplitudes, shifted.hamp, reference_norm, shifted_norm) /
                       static_cast<double>(angular_points);
  }

  INFO(label);
  REQUIRE(angular_overlap.real() > 0.9999);
  REQUIRE(std::abs(angular_overlap.imag()) < 1.0e-6);
  const double component_scale = std::sqrt(reference_norm / static_cast<double>(reference_amplitudes.size()));
  for (std::size_t index = 0; index < reference_amplitudes.size(); ++index) {
    CAPTURE(index, reference_amplitudes[index], angular_average[index]);
    const double tolerance = 5.0e-4 * std::max(component_scale, std::abs(reference_amplitudes[index]));
    REQUIRE(std::abs(angular_average[index] - reference_amplitudes[index]) < tolerance);
  }
}

// Build one shared MG5 momentum pointer array from stable decay leaves
void BuildGeneratedTestMomenta(const gra::LORENTZSCALAR &lts, std::vector<std::array<double, 4>> &storage,
                               std::vector<double *> &momenta) {
  const auto              leaves = gra::mg5::StableDecayLeaves(lts.decaytree);
  std::vector<gra::M4Vec> final  = gra::mg5::StableLeafMomenta(leaves);
  gra::M4Vec              p1     = lts.q1;
  gra::M4Vec              p2     = lts.q2;
  REQUIRE(gra::mg5helas::PrepareOnShellKinematics(lts, final, p1, p2));
  gra::mg5::BuildMG5Momenta(p1, p2, final, storage, momenta);
}

// Evaluate generated channels with the complete MadGraph denominator
template <class ProcessBase>
double GeneratedFullDenominatorAmp2(std::vector<gra::mg5::Subprocess<ProcessBase>> &subprocesses,
                                    const gra::LORENTZSCALAR &lts, double alpha_s) {
  std::vector<std::array<double, 4>> storage;
  std::vector<double *>              momenta;
  BuildGeneratedTestMomenta(lts, storage, momenta);
  const auto       leaves = gra::mg5::StableDecayLeaves(lts.decaytree);
  std::vector<int> final_pdgs;
  final_pdgs.reserve(leaves.size());
  for (const gra::MDecayBranch *leaf : leaves) { final_pdgs.push_back(leaf->p.pdg); }
  const gra::AmplitudeTopology topology = gra::amplitude::AmplitudeTopologyFromDecayTree(lts.decaytree);

  double amp2 = 0.0;
  for (auto &subprocess : subprocesses) {
    std::size_t available = 0;
    std::size_t selected  = 0;
    for (const auto &channel : subprocess.channels) {
      if (channel.initial != std::array<int, 2>{lts.id1, lts.id2}) { continue; }
      ++available;
      if (gra::mg5::ChannelMatches(channel, topology, final_pdgs)) { ++selected; }
    }
    if (selected == 0) { continue; }
    subprocess.process->setMomenta(momenta);
    subprocess.process->setInitial(lts.id1, lts.id2);
    subprocess.process->setAlphaS(alpha_s);
    subprocess.process->sigmaKin();
    amp2 += static_cast<double>(selected) / static_cast<double>(available) * subprocess.process->sigmaHat();
  }
  return amp2;
}

// Check every subprocess channel against its generated scalar MG5 color sum
template <class ProcessBase>
void RequireGeneratedChannelNorms(std::vector<gra::mg5::Subprocess<ProcessBase>> &subprocesses,
                                  const gra::LORENTZSCALAR &lts, double alpha_s, const std::string &family) {
  std::vector<std::array<double, 4>> storage;
  std::vector<double *>              momenta;
  const auto seed = gra::mg5::StableLeafMomenta(gra::mg5::StableDecayLeaves(lts.decaytree));
  for (auto &subprocess : subprocesses) {
    auto final = seed;
    gra::M4Vec total;
    const auto &masses = subprocess.process->getMasses();
    REQUIRE(masses.size() == final.size() + 2);
    for (const auto &i : indices(final)) {
      final[i].SetPxPyPzM(final[i].Px(), final[i].Py(), final[i].Pz(), masses[i + 2]);
      total += final[i];
    }
    REQUIRE(total.P3mod() < 1e-12 * total.E());
    REQUIRE(gra::mg5::OnShellFinal(final, masses));
    const double energy = total.E() / 2.0;
    gra::mg5::BuildMG5Momenta(gra::M4Vec(0, 0, energy, energy),
                            gra::M4Vec(0, 0, -energy, energy), final, storage, momenta);
    subprocess.process->setMomenta(momenta);
    subprocess.process->setAlphaS(alpha_s);
    subprocess.process->sigmaKin();
    for (const auto &channel : subprocess.channels) {
      subprocess.process->setInitial(channel.initial[0], channel.initial[1]);
      const double scalar              = subprocess.process->sigmaHat();
      const auto   components          = subprocess.process->helicityAmplitudes();
      const double structured          = gra::mg5helas::ComponentNorm(components);
      bool         has_flow_amplitudes = false;
      for (const auto &component : components) {
        if (!component.flow_values.empty()) {
          has_flow_amplitudes = true;
          REQUIRE(component.color == 0);
          REQUIRE(component.flow_values.size() == channel.external_color_flows.size());
        } else {
          REQUIRE(component.color != 0);
        }
      }
      INFO(family << " initial=" << channel.initial[0] << "," << channel.initial[1]);
      REQUIRE(structured == Approx(scalar).epsilon(1e-10));
      REQUIRE(has_flow_amplitudes);
      REQUIRE(channel.external_color_representations.size() == channel.final.size() + 2);
      REQUIRE_FALSE(channel.external_color_flows.empty());
    }
  }
}

// Check generated processes without external channel data
template <class ProcessBase>
void RequireGeneratedProcessNorms(std::vector<std::unique_ptr<ProcessBase>> &processes, const gra::LORENTZSCALAR &lts,
                                  const std::array<int, 2> &initial, double alpha_s, const std::string &family) {
  std::vector<std::array<double, 4>> storage;
  std::vector<double *>              momenta;
  BuildGeneratedTestMomenta(lts, storage, momenta);

  for (auto &process : processes) {
    process->setMomenta(momenta);
    process->setInitial(initial[0], initial[1]);
    process->setAlphaS(alpha_s);
    process->sigmaKin();
    const double scalar              = process->sigmaHat();
    const auto   components          = process->helicityAmplitudes();
    const double structured          = gra::mg5helas::ComponentNorm(components);
    bool         has_flow_amplitudes = false;
    for (const auto &component : components) {
      if (!component.flow_values.empty()) { has_flow_amplitudes = true; }
    }
    INFO(family);
    REQUIRE(structured == Approx(scalar).epsilon(1e-10));
    REQUIRE(has_flow_amplitudes);
  }
}

// Check one fixed matrix-element reference value
void RequireAmplitudeReference(const std::string &label, double amp2, double reference) {
  INFO(label);
  std::ostringstream values;
  values << std::scientific << std::setprecision(17) << "amp2=" << amp2 << " reference=" << reference;
  INFO(values.str());
  REQUIRE(std::isfinite(amp2));
  REQUIRE(amp2 == Approx(reference).epsilon(1e-12));
}

// Read one HepMC integer attribute with a zero fallback
int ReadIntAttribute(const HepMC3::ConstGenParticlePtr &particle, const std::string &name) {
  const auto attribute = particle->attribute<HepMC3::IntAttribute>(name);
  return attribute ? attribute->value() : 0;
}

// Compute true for shower-colored parton PDG ids
bool IsColoredParton(int pdg) {
  const int apdg = std::abs(pdg);
  const int spin = apdg % 10;
  return pdg == 21 || (apdg >= 1 && apdg <= 5) || (apdg >= 1000 && apdg < 6000 && (spin == 1 || spin == 3));
}

// Expose the protected amplitude call for process-level hard-diffraction tests
class HardDiffractionTestProbe : public gra::MHardDiffraction {
 public:
  using gra::MHardDiffraction::BuildHardPair;
  using gra::MHardDiffraction::HardLongPoint;
  using gra::MHardDiffraction::MapExp;
  using gra::MHardDiffraction::MapHardLong;
  using gra::MHardDiffraction::MapLog;
  using gra::MHardDiffraction::MHardDiffraction;
  using gra::MHardDiffraction::ScaleDoubleDiffractivePair;

  // Rebuild hard-diffraction invariants with mixed incoming virtualities
  bool ProbeHardLorentzScalars() {
    return gra::kinematics::SetLorentzScalars(state, 2, false,
                                              gra::kinematics::TransferVirtualityPolicy::HardDiffractive);
  }

  // Rebuild invariants with the strict generic spacelike prescription
  bool ProbeSpacelikeLorentzScalars() { return gra::kinematics::SetLorentzScalars(state, 2); }

  // Compute the process-level matrix element squared and its sampling
  // bookkeeping
  double ProbeAmp2(gra::MEventWeightState &aux) { return GetAmp2(true, aux); }

  // Finalize one directly prepared test amplitude before record construction
  bool FinalizeProbe() {
    amplitude_event_state_finalized = FinalizeEventState();
    return amplitude_event_state_finalized;
  }

  // Build one forward string through the common N-star implementation
  bool ProbeExciteString(const gra::M4Vec &nstar, gra::MDecayBranch &forward, const gra::MParticle &beam,
                         int color_tag) {
    return ExciteString(nstar, forward, beam, color_tag);
  }

  // Fragment the configured forward systems through the common implementation
  bool ProbeForwardFragment() { return CEPForwardFragment(); }
};

// Prepare a physical forward-remnant point for hard-color-flow process tests
void ConfigureHardColorProbe(HardDiffractionTestProbe &proc, bool double_diffraction,
                             const std::vector<gra::MDecayBranch> &decaytree, double hard_mass = 0.0) {
  proc.SetInitialState({"p+", "p+"}, {6500.0, 6500.0});
  proc.SetScreening(false);
  const double beam_pz                   = std::sqrt(gra::math::pow2(6500.0) - gra::math::pow2(gra::PDG::mp));
  proc.state.lts.pbeam1                  = gra::M4Vec(0.0, 0.0, beam_pz, 6500.0);
  proc.state.lts.pbeam2                  = gra::M4Vec(0.0, 0.0, -beam_pz, 6500.0);
  proc.state.lts.s                       = (proc.state.lts.pbeam1 + proc.state.lts.pbeam2).M2();
  proc.state.lts.sqrt_s                  = std::sqrt(proc.state.lts.s);
  proc.state.lts.LHAPDFSET               = "MMHT2014lo68cl";
  proc.state.lts.process.root_decay_mode = gra::RootDecayMode::Physical;
  proc.state.lts.decaytree               = decaytree;

  gra::M4Vec reference_hard;
  for (const auto &branch : proc.state.lts.decaytree) { reference_hard += branch.p4; }
  if (hard_mass > 0.0) {
    const double scale = hard_mass / reference_hard.M();
    // Rescale spatial momenta while retaining model masses and decay closure
    const std::function<void(gra::MDecayBranch &)> rescale = [&](gra::MDecayBranch &branch) {
      if (branch.legs.empty()) {
        branch.p4.SetPxPyPzM(scale * branch.p4.Px(), scale * branch.p4.Py(),
                            scale * branch.p4.Pz(), branch.p.mass);
      } else {
        branch.p4 = gra::M4Vec();
        for (auto &child : branch.legs) {
          rescale(child);
          branch.p4 += child.p4;
        }
      }
    };
    reference_hard = gra::M4Vec();
    for (auto &branch : proc.state.lts.decaytree) {
      rescale(branch);
      reference_hard += branch.p4;
    }
  }
  proc.state.gcuts.M_min = 0.9 * reference_hard.M();
  proc.state.gcuts.M_max = 1.1 * reference_hard.M();
  const double xhard = reference_hard.M() / proc.state.lts.sqrt_s;

  proc.state.lts.hard_diff1  = true;
  proc.state.lts.diff_xi1    = 0.04;
  proc.state.lts.diff_beta1  = xhard / proc.state.lts.diff_xi1;
  proc.state.lts.diff_t1     = -0.1;
  proc.state.lts.diff_phi1   = 0.3;
  proc.state.lts.diff_xhard1 = xhard;
  proc.state.lts.diff_xhard2 = proc.state.lts.diff_xhard1;
  if (double_diffraction) {
    proc.state.lts.hard_diff2 = true;
    proc.state.lts.diff_xi2   = proc.state.lts.diff_xi1;
    proc.state.lts.diff_beta2 = proc.state.lts.diff_beta1;
    proc.state.lts.diff_t2    = proc.state.lts.diff_t1;
    proc.state.lts.diff_phi2  = proc.state.lts.diff_phi1 + gra::math::PI;
  } else {
    proc.state.lts.hard_diff2 = false;
  }
  proc.SetModelTune(gra::MModelTune::Load(modelfile));
  proc.FinalizeProcessConfiguration();

  REQUIRE(
      proc.BuildHardPair(proc.state.lts.diff_xhard1, proc.state.lts.diff_xhard2, proc.state.lts.q1, proc.state.lts.q2));

  proc.state.lts.pfinal.assign(decaytree.size() + 3, gra::M4Vec(0.0, 0.0, 0.0, 0.0));
  proc.state.lts.pfinal[0]     = proc.state.lts.q1 + proc.state.lts.q2;
  proc.state.lts.pfinal[1]     = proc.state.lts.pbeam1 - proc.state.lts.q1;
  proc.state.lts.pfinal[2]     = proc.state.lts.pbeam2 - proc.state.lts.q2;
  proc.state.lts.forward_mass2 = {proc.state.lts.pfinal[1].M2(), proc.state.lts.pfinal[2].M2()};
  REQUIRE(proc.state.lts.pfinal[0].M() == Approx(reference_hard.M()).epsilon(1.0e-11));
  for (const auto &i : indices(proc.state.lts.decaytree)) {
    TransportGeneratedBranch(proc.state.lts.decaytree[i], reference_hard, proc.state.lts.pfinal[0]);
    proc.state.lts.pfinal[i + 3] = proc.state.lts.decaytree[i].p4;
  }
}

// Sum the stored amplitude-component norms
double HampNormSum(const gra::LORENTZSCALAR &lts) {
  double sum = 0.0;
  for (const auto &amp : lts.hamp) { sum += std::norm(amp); }
  return sum;
}

// Compute true when at least one stored component retains a physical phase
bool HasComplexHamp(const gra::LORENTZSCALAR &lts) {
  for (const auto &amp : lts.hamp) {
    if (std::abs(amp.imag()) > 1.0e-15) { return true; }
  }
  return false;
}

// Test that excited baryons become closed showerable valence strings
TEST_CASE("N-star string skeleton has baryon color closure", "[HardPomeronPDF]") {
  HardDiffractionTestProbe proc;
  proc.SetModelTune(gra::MModelTune::Load(modelfile));
  proc.state.lts.PDG.ReadParticleData();

  for (const int beam_pdg : {2212, -2212, 2112, -2112}) {
    const gra::MParticle beam = proc.state.lts.PDG.FindByPDG(beam_pdg);
    const gra::M4Vec     nstar(0.3, -0.2, beam_pdg > 0 ? 4.0 : -4.0, std::sqrt(4.0 * 4.0 + 2.0 * 2.0 + 0.13));
    gra::MDecayBranch    forward;
    REQUIRE(proc.ProbeExciteString(nstar, forward, beam, 701));
    REQUIRE(forward.legs.size() == 2);

    const auto &quark   = forward.legs[0];
    const auto &diquark = forward.legs[1];
    REQUIRE(std::abs(quark.p.pdg) >= 1);
    REQUIRE(std::abs(quark.p.pdg) <= 2);
    REQUIRE(std::abs(diquark.p.pdg) >= 1000);
    REQUIRE((quark.p.pdg > 0) == (beam_pdg > 0));
    REQUIRE((diquark.p.pdg > 0) == (beam_pdg > 0));
    if (beam_pdg > 0) {
      REQUIRE(quark.p.color_flow.flow1 == 701);
      REQUIRE(quark.p.color_flow.flow2 == 0);
      REQUIRE(diquark.p.color_flow.flow1 == 0);
      REQUIRE(diquark.p.color_flow.flow2 == 701);
    } else {
      REQUIRE(quark.p.color_flow.flow1 == 0);
      REQUIRE(quark.p.color_flow.flow2 == 701);
      REQUIRE(diquark.p.color_flow.flow1 == 701);
      REQUIRE(diquark.p.color_flow.flow2 == 0);
    }
    REQUIRE(gra::math::CheckEMC(nstar - quark.p4 - diquark.p4));
  }

  gra::MParticle unsupported_beam;
  unsupported_beam.pdg = 3122;
  gra::MDecayBranch unsupported_forward;
  REQUIRE_THROWS_AS(proc.ProbeExciteString(gra::M4Vec(0.0, 0.0, 0.0, 2.0), unsupported_forward, unsupported_beam, 701),
                    std::invalid_argument);
}

// Test that cylinder fragmentation keeps neutron charge and baryon number
TEST_CASE("N-star cylinder fragmentation conserves neutron quantum numbers", "[HardPomeronPDF]") {
  HardDiffractionTestProbe proc;
  proc.SetModelTune(gra::MModelTune::Load(modelfile));
  proc.SetBeamFrag("cylinder");
  proc.state.lts.PDG.ReadParticleData();
  proc.state.lts.beam1   = proc.state.lts.PDG.FindByPDG(gra::PDG::PDG_n);
  proc.state.lts.beam2   = proc.state.lts.PDG.FindByPDG(gra::PDG::PDG_p);
  proc.state.lts.excite1 = true;
  proc.state.lts.excite2 = false;
  proc.state.lts.pfinal.resize(3);
  const double mass        = 6.0;
  const double pz          = 10.0;
  proc.state.lts.pfinal[1] = gra::M4Vec(0.0, 0.0, pz, std::sqrt(mass * mass + pz * pz));

  REQUIRE(proc.ProbeForwardFragment());
  REQUIRE_FALSE(proc.state.lts.decayforward1.legs.empty());

  int charge = 0;
  int baryon = 0;
  for (const auto &leg : proc.state.lts.decayforward1.legs) {
    charge += leg.p.chargeX3 / 3;
    if (std::abs(leg.p.pdg) == gra::PDG::PDG_p || std::abs(leg.p.pdg) == gra::PDG::PDG_n) {
      baryon += gra::math::sign(leg.p.pdg);
    }
  }
  REQUIRE(charge == 0);
  REQUIRE(baryon == 1);
}

// Test factorized hard-Pomeron PDF access from the default model card
TEST_CASE("Hard Pomeron PDF evaluates GKG18 DPDF", "[HardPomeronPDF]") {
  const std::string          modelfile = gra::ResolveModelDataFile("TUNE0", "GENERAL.json");
  const gra::MHardPomeronPDF pdf(modelfile);

  const double xi   = 0.01;
  const double beta = 0.2;
  const double t    = -0.1;
  const double Q2   = 91.1876 * 91.1876;

  const double flux    = pdf.Flux(xi, t);
  const double gluon   = pdf.PartonDensity(21, beta, Q2);
  const double density = pdf.DiffractiveDensity(21, xi, beta, t, Q2);

  REQUIRE(std::isfinite(flux));
  REQUIRE(std::isfinite(gluon));
  REQUIRE(std::isfinite(density));
  REQUIRE(flux > 0.0);
  REQUIRE(gluon >= 0.0);
  REQUIRE(density >= 0.0);
  REQUIRE(pdf.HardX(xi, beta) == Approx(xi * beta));
  REQUIRE(pdf.Flux(xi, -0.2) / flux == Approx(std::exp(pdf.FluxTSlope(xi) * (-0.2 - t))).epsilon(1e-14));

  gra::MLHAPDFStore pdf_store;
  const auto        raw_dpdf = pdf_store.GetPDF("GKG18_DPDF_FitB_LO", 0);
  REQUIRE(pdf.AlphaS(Q2) == Approx(raw_dpdf->alphasQ2(Q2)).epsilon(1e-14));
  REQUIRE(pdf.MatchAlphaS(*pdf_store.GetPDF("MMHT2014lo68cl", 0), Q2));
  const auto proton = pdf_store.GetPDF("MMHT2014lo68cl", 0);
  REQUIRE_NOTHROW(pdf.ValidateAlphaS(*proton, 80.0 * 80.0, 100.0 * 100.0));
  REQUIRE_FALSE(pdf.MatchAlphaS(*proton, 200.0 * 200.0));
  REQUIRE_NOTHROW(pdf.ValidateAlphaS(*proton, 80.0 * 80.0, 320.0 * 320.0));
  REQUIRE_FALSE(pdf.MatchAlphaS(*pdf_store.GetPDF("NNPDF31_lo_as_0118", 0), Q2));
  REQUIRE_NOTHROW(pdf.ValidateAlphaS(*pdf_store.GetPDF("NNPDF31_lo_as_0118", 0), 80.0 * 80.0, 100.0 * 100.0));
  REQUIRE(pdf.AlphaS(Q2) == Approx(raw_dpdf->alphasQ2(Q2)).epsilon(1e-14));
  REQUIRE_THROWS_AS(pdf.ValidateAlphaS(*proton, 0.0, Q2), std::invalid_argument);
  REQUIRE_THROWS_AS(pdf.AlphaS(0.0), std::invalid_argument);

  REQUIRE(pdf.Flux(0.5, t) == Approx(0.0));
  REQUIRE(pdf.PartonDensity(21, 1.5, Q2) == Approx(0.0));
  REQUIRE(pdf.PartonDensity(21, beta, std::numeric_limits<double>::infinity()) == Approx(0.0));
  REQUIRE(pdf.DiffractiveDensity(21, xi, beta, -10.0, Q2) == Approx(0.0));
}

// Test exact one-dimensional proposal maps and their Jacobians
TEST_CASE("Hard diffraction proposal maps carry exact Jacobians", "[HardPomeronPDF][Integration]") {
  constexpr double step = 1.0e-6;
  const double     unit = 0.37;

  double value    = 0.0;
  double jacobian = 0.0;
  REQUIRE(HardDiffractionTestProbe::MapLog(unit, 1.0e-4, 0.1, value, jacobian));
  double value_plus     = 0.0;
  double jacobian_plus  = 0.0;
  double value_minus    = 0.0;
  double jacobian_minus = 0.0;
  REQUIRE(HardDiffractionTestProbe::MapLog(unit + step, 1.0e-4, 0.1, value_plus, jacobian_plus));
  REQUIRE(HardDiffractionTestProbe::MapLog(unit - step, 1.0e-4, 0.1, value_minus, jacobian_minus));
  REQUIRE((value_plus - value_minus) / (2.0 * step) == Approx(jacobian).epsilon(1.0e-9));

  for (const double slope : {0.0, 7.2, -3.4}) {
    REQUIRE(HardDiffractionTestProbe::MapExp(unit, -1.0, -0.002, slope, value, jacobian));
    REQUIRE(HardDiffractionTestProbe::MapExp(unit + step, -1.0, -0.002, slope, value_plus, jacobian_plus));
    REQUIRE(HardDiffractionTestProbe::MapExp(unit - step, -1.0, -0.002, slope, value_minus, jacobian_minus));
    INFO(slope);
    REQUIRE((value_plus - value_minus) / (2.0 * step) == Approx(jacobian).epsilon(1.0e-9));

    const double weighted_jacobian = std::exp(slope * value) * jacobian;
    const double expected = std::abs(slope) < 1.0e-12 ? 0.998 : (std::exp(-0.002 * slope) - std::exp(-slope)) / slope;
    REQUIRE(weighted_jacobian == Approx(expected).epsilon(1.0e-12));
  }

  REQUIRE_FALSE(HardDiffractionTestProbe::MapLog(unit, 0.0, 0.1, value, jacobian));
  REQUIRE_FALSE(HardDiffractionTestProbe::MapExp(unit, -0.1, -1.0, 2.0, value, jacobian));
}

// Test exact beta and x transformations through mass and rapidity
TEST_CASE("Hard diffraction beta-x map has exact Jacobian", "[HardPomeronPDF][Integration]") {
  constexpr double step          = 1.0e-6;
  const double     s             = 13000.0 * 13000.0;
  const double     mass_unit     = 0.43;
  const double     rapidity_unit = 0.37;

  const auto require_map = [&](const std::array<double, 2> &xhard1_range, const std::array<double, 2> &xhard2_range,
                               double xi_product) {
    const auto map = [&](double mass_coordinate, double rapidity_coordinate) {
      HardDiffractionTestProbe::HardLongPoint point;
      REQUIRE(HardDiffractionTestProbe::MapHardLong(mass_coordinate, rapidity_coordinate, s, 60.0 * 60.0, 120.0 * 120.0,
                                                    xhard1_range, xhard2_range, xi_product, point));
      return point;
    };

    const auto point          = map(mass_unit, rapidity_unit);
    const auto mass_plus      = map(mass_unit + step, rapidity_unit);
    const auto mass_minus     = map(mass_unit - step, rapidity_unit);
    const auto rapidity_plus  = map(mass_unit, rapidity_unit + step);
    const auto rapidity_minus = map(mass_unit, rapidity_unit - step);

    REQUIRE((mass_plus.mass2 - mass_minus.mass2) / (2.0 * step) ==
            Approx(point.mass_jacobian).epsilon(1e-9));
    REQUIRE(s * point.xhard1 * point.xhard2 == Approx(point.mass2).epsilon(1.0e-13));
    REQUIRE(0.5 * std::log(point.xhard1 / point.xhard2) == Approx(point.rapidity).epsilon(1.0e-13));
    REQUIRE(point.xhard1 >= xhard1_range[0]);
    REQUIRE(point.xhard1 <= xhard1_range[1]);
    REQUIRE(point.xhard2 >= xhard2_range[0]);
    REQUIRE(point.xhard2 <= xhard2_range[1]);

    const double dx1_mass     = (mass_plus.xhard1 - mass_minus.xhard1) / (2.0 * step);
    const double dx2_mass     = (mass_plus.xhard2 - mass_minus.xhard2) / (2.0 * step);
    const double dx1_rapidity = (rapidity_plus.xhard1 - rapidity_minus.xhard1) / (2.0 * step);
    const double dx2_rapidity = (rapidity_plus.xhard2 - rapidity_minus.xhard2) / (2.0 * step);
    const double determinant  = std::abs(dx1_mass * dx2_rapidity - dx1_rapidity * dx2_mass);
    REQUIRE(determinant / xi_product == Approx(point.jacobian).epsilon(2.0e-8));
  };

  SECTION("single diffraction beta and x map") {
    const double xi = 0.08;
    require_map({xi * 0.02, xi * 0.9}, {1.0e-5, 0.8}, xi);
  }

  SECTION("double diffraction beta and beta map") {
    const double xi1 = 0.08;
    const double xi2 = 0.06;
    require_map({xi1 * 0.02, xi1 * 0.9}, {xi2 * 0.03, xi2 * 0.95}, xi1 * xi2);
  }

  HardDiffractionTestProbe::HardLongPoint invalid;
  REQUIRE_FALSE(HardDiffractionTestProbe::MapHardLong(mass_unit, rapidity_unit, s, 120.0 * 120.0, 60.0 * 60.0,
                                                      {0.001, 0.08}, {0.001, 0.08}, 0.08, invalid));

  HardDiffractionTestProbe::HardLongPoint threshold;
  const std::array<double, 2>             low_xhard1 = {1.0e-5, 0.08};
  const std::array<double, 2>             low_xhard2 = {1.0e-6, 0.8};
  REQUIRE(HardDiffractionTestProbe::MapHardLong(1.0e-7, rapidity_unit, s, 0.0, 120.0 * 120.0, low_xhard1, low_xhard2,
                                                0.08, threshold));
  const double reachable_threshold = s * low_xhard1[0] * low_xhard2[0];
  REQUIRE(threshold.mass2 > reachable_threshold);
  REQUIRE(threshold.mass2 / reachable_threshold < 1.00001);
}

// Test mixed spacelike and lightlike hard SD transfer reconstruction
TEST_CASE("Hard SD repairs only its declared lightlike transfer", "[HardPomeronPDF][kinematics]") {
  constexpr double energy          = 6500.0;
  constexpr double x               = 0.236;
  constexpr double diff_virtuality = -0.01;
  constexpr double diff_scale      = 130.0;
  constexpr double diff_px         = 0.02;

  const auto configure = [&](HardDiffractionTestProbe &proc, bool first_diff, double ordinary_virtuality) {
    const double beam_pz  = std::sqrt(energy * energy - gra::math::pow2(gra::PDG::mp));
    proc.state.lts.pbeam1 = gra::M4Vec(0.0, 0.0, beam_pz, energy);
    proc.state.lts.pbeam2 = gra::M4Vec(0.0, 0.0, -beam_pz, energy);

    gra::M4Vec q1;
    gra::M4Vec q2;
    if (first_diff) {
      q1 = LightConeTestVector(diff_scale, (diff_virtuality + diff_px * diff_px) / diff_scale, diff_px, 0.0);
      const double minus = x * proc.state.lts.pbeam2.LightconeNeg();
      q2                 = LightConeTestVector(ordinary_virtuality / minus, minus, 0.0, 0.0);
    } else {
      const double plus = x * proc.state.lts.pbeam1.LightconePos();
      q1                = LightConeTestVector(plus, ordinary_virtuality / plus, 0.0, 0.0);
      q2 = LightConeTestVector((diff_virtuality + diff_px * diff_px) / diff_scale, diff_scale, diff_px, 0.0);
    }

    proc.state.lts.pfinal.assign(3, gra::M4Vec());
    proc.state.lts.pfinal[0]  = q1 + q2;
    proc.state.lts.pfinal[1]  = proc.state.lts.pbeam1 - q1;
    proc.state.lts.pfinal[2]  = proc.state.lts.pbeam2 - q2;
    proc.state.lts.hard_diff1 = first_diff;
    proc.state.lts.hard_diff2 = !first_diff;
    proc.state.lts.q1         = gra::M4Vec(1.0, 2.0, 3.0, 4.0);
    proc.state.lts.q2         = gra::M4Vec(-1.0, -2.0, -3.0, 5.0);
    proc.state.lts.t1         = -17.0;
    proc.state.lts.t2         = -19.0;
  };

  const auto require_repaired = [&](bool first_diff) {
    HardDiffractionTestProbe proc;
    configure(proc, first_diff, 1.0e-7);
    const gra::M4Vec reconstructed = first_diff ? proc.state.lts.pbeam2 - proc.state.lts.pfinal[2]
                                                : proc.state.lts.pbeam1 - proc.state.lts.pfinal[1];
    REQUIRE(reconstructed.M2() > 0.0);
    REQUIRE(proc.ProbeHardLorentzScalars());
    const gra::M4Vec &ordinary   = first_diff ? proc.state.lts.q2 : proc.state.lts.q1;
    const double      diff_t     = first_diff ? proc.state.lts.t1 : proc.state.lts.t2;
    const double      ordinary_t = first_diff ? proc.state.lts.t2 : proc.state.lts.t1;
    REQUIRE(diff_t == Approx(diff_virtuality).margin(1.0e-10));
    REQUIRE(ordinary_t == Approx(0.0).margin(1.0e-15));
    REQUIRE(ordinary.E() == Approx(ordinary.P3mod()).epsilon(1.0e-15));
    REQUIRE(ordinary.M2() == Approx(0.0).margin(1.0e-9));
  };

  SECTION("diffractive side one") { require_repaired(true); }
  SECTION("diffractive side two") { require_repaired(false); }

  SECTION("genuinely timelike ordinary transfer is rejected atomically") {
    HardDiffractionTestProbe proc;
    configure(proc, true, 1.0e-3);
    REQUIRE_FALSE(proc.ProbeHardLorentzScalars());
    REQUIRE(proc.state.lts.q1.E() == Approx(4.0));
    REQUIRE(proc.state.lts.q2.E() == Approx(5.0));
    REQUIRE(proc.state.lts.t1 == Approx(-17.0));
    REQUIRE(proc.state.lts.t2 == Approx(-19.0));
  }

  SECTION("generic spacelike policy stays strict") {
    HardDiffractionTestProbe proc;
    configure(proc, true, 1.0e-7);
    REQUIRE_FALSE(proc.ProbeSpacelikeLorentzScalars());
    REQUIRE(proc.state.lts.q1.E() == Approx(4.0));
    REQUIRE(proc.state.lts.q2.E() == Approx(5.0));
  }

  SECTION("double diffraction accepts two spacelike transfers") {
    HardDiffractionTestProbe proc;
    configure(proc, true, -2.0e-2);
    proc.state.lts.hard_diff2 = true;
    REQUIRE(proc.ProbeHardLorentzScalars());
    REQUIRE(proc.state.lts.t1 < 0.0);
    REQUIRE(proc.state.lts.t2 < 0.0);
  }
}

// Test the finite t recoil map against the sampled collinear hard mass
TEST_CASE("Hard SD recoil preserves the DPDF hard invariant mass", "[HardPomeronPDF][kinematics]") {
  constexpr double energy     = 6500.0;
  constexpr double xi         = 0.03;
  constexpr double beta       = 0.24;
  constexpr double t          = -0.30;
  constexpr double phi        = 0.70;
  constexpr double ordinary_x = 0.012;
  const double     diff_x     = xi * beta;

  const auto configure = [&](HardDiffractionTestProbe &proc, bool first_diff, double proton_x) {
    proc.SetModelTune(gra::MModelTune::Load(modelfile));
    proc.FinalizeProcessConfiguration();
    const double beam_pz      = std::sqrt(energy * energy - gra::math::pow2(gra::PDG::mp));
    proc.state.lts.pbeam1     = gra::M4Vec(0.0, 0.0, beam_pz, energy);
    proc.state.lts.pbeam2     = gra::M4Vec(0.0, 0.0, -beam_pz, energy);
    proc.state.lts.s          = (proc.state.lts.pbeam1 + proc.state.lts.pbeam2).M2();
    proc.state.lts.hard_diff1 = first_diff;
    proc.state.lts.hard_diff2 = !first_diff;
    if (first_diff) {
      proc.state.lts.diff_xi1    = xi;
      proc.state.lts.diff_beta1  = beta;
      proc.state.lts.diff_t1     = t;
      proc.state.lts.diff_phi1   = phi;
      proc.state.lts.diff_xhard1 = diff_x;
      proc.state.lts.diff_xhard2 = proton_x;
    } else {
      proc.state.lts.diff_xi2    = xi;
      proc.state.lts.diff_beta2  = beta;
      proc.state.lts.diff_t2     = t;
      proc.state.lts.diff_phi2   = phi;
      proc.state.lts.diff_xhard1 = proton_x;
      proc.state.lts.diff_xhard2 = diff_x;
    }
  };

  const auto require_mass_preserved = [&](bool first_diff) {
    HardDiffractionTestProbe proc;
    configure(proc, first_diff, ordinary_x);
    const double x1    = proc.state.lts.diff_xhard1;
    const double x2    = proc.state.lts.diff_xhard2;
    const double mass2 = proc.state.lts.s * x1 * x2;
    gra::M4Vec   q1;
    gra::M4Vec   q2;
    REQUIRE(proc.BuildHardPair(x1, x2, q1, q2));

    const gra::M4Vec &qdiff     = first_diff ? q1 : q2;
    const gra::M4Vec &qcol      = first_diff ? q2 : q1;
    const gra::M4Vec &beam_diff = first_diff ? proc.state.lts.pbeam1 : proc.state.lts.pbeam2;
    const gra::M4Vec &beam_col  = first_diff ? proc.state.lts.pbeam2 : proc.state.lts.pbeam1;
    REQUIRE(qdiff.M2() == Approx(beta * t).margin(1.0e-8));
    REQUIRE(qcol.M2() < 0.0);
    REQUIRE((q1 + q2).M2() == Approx(mass2).epsilon(1.0e-11));

    const double record_x =
        first_diff ? qcol.LightconeNeg() / beam_col.LightconeNeg() : qcol.LightconePos() / beam_col.LightconePos();
    REQUIRE(record_x > ordinary_x);
    REQUIRE(record_x < 1.0);

    gra::M4Vec leading;
    REQUIRE(gra::kinematics::BuildForwardParticleXiT(beam_diff, xi, t, phi, first_diff, leading));
    const gra::M4Vec pomeron_remnant = beam_diff - qdiff - leading;
    const gra::M4Vec proton_remnant  = beam_col - qcol;
    REQUIRE(pomeron_remnant.E() > 0.0);
    REQUIRE(proton_remnant.E() > 0.0);
    REQUIRE(pomeron_remnant.M2() == Approx(0.0).margin(1.0e-8));
    const auto pdf = proc.state.lts.model_cache->hard_pomeron.GetHardPomeronPDF(proc.GetSoftModel());
    REQUIRE(proton_remnant.M() == Approx(pdf->RemnantMass()).epsilon(1.0e-8));
    REQUIRE(gra::math::CheckEMC(proc.state.lts.pbeam1 + proc.state.lts.pbeam2 -
                                (leading + pomeron_remnant + proton_remnant + q1 + q2)));
  };

  SECTION("diffractive side one") { require_mass_preserved(true); }
  SECTION("diffractive side two") { require_mass_preserved(false); }

  SECTION("ordinary remnant support is enforced") {
    HardDiffractionTestProbe proc;
    const double             boundary_x = 1.0 - 1.0e-10;
    configure(proc, true, boundary_x);
    gra::M4Vec q1;
    gra::M4Vec q2;
    REQUIRE_FALSE(proc.BuildHardPair(diff_x, boundary_x, q1, q2));
  }
}

// Test the symmetric finite t DD recoil against the exact DPDF hard mass
TEST_CASE("Hard DD recoil preserves the DPDF hard invariant mass", "[HardPomeronPDF][kinematics]") {
  constexpr double energy = 6500.0;
  constexpr double xi1    = 0.030;
  constexpr double xi2    = 0.047;
  constexpr double beta1  = 0.24;
  constexpr double beta2  = 0.31;
  constexpr double t1     = -0.30;
  constexpr double t2     = -0.55;
  constexpr double phi1   = 0.70;
  constexpr double phi2   = -1.10;

  HardDiffractionTestProbe proc;
  const double             beam_pz = std::sqrt(energy * energy - gra::math::pow2(gra::PDG::mp));
  proc.state.lts.pbeam1            = gra::M4Vec(0.0, 0.0, beam_pz, energy);
  proc.state.lts.pbeam2            = gra::M4Vec(0.0, 0.0, -beam_pz, energy);
  proc.state.lts.s                 = (proc.state.lts.pbeam1 + proc.state.lts.pbeam2).M2();
  proc.state.lts.hard_diff1        = true;
  proc.state.lts.hard_diff2        = true;
  proc.state.lts.diff_xi1          = xi1;
  proc.state.lts.diff_xi2          = xi2;
  proc.state.lts.diff_beta1        = beta1;
  proc.state.lts.diff_beta2        = beta2;
  proc.state.lts.diff_t1           = t1;
  proc.state.lts.diff_t2           = t2;
  proc.state.lts.diff_phi1         = phi1;
  proc.state.lts.diff_phi2         = phi2;
  proc.state.lts.diff_xhard1       = xi1 * beta1;
  proc.state.lts.diff_xhard2       = xi2 * beta2;

  gra::M4Vec leading1;
  gra::M4Vec leading2;
  REQUIRE(gra::kinematics::BuildForwardParticleXiT(proc.state.lts.pbeam1, xi1, t1, phi1, true, leading1));
  REQUIRE(gra::kinematics::BuildForwardParticleXiT(proc.state.lts.pbeam2, xi2, t2, phi2, false, leading2));
  const gra::M4Vec exchange1 = proc.state.lts.pbeam1 - leading1;
  const gra::M4Vec exchange2 = proc.state.lts.pbeam2 - leading2;
  const gra::M4Vec base1     = DiffractiveHardPartonTestVector(exchange1, beta1, true);
  const gra::M4Vec base2     = DiffractiveHardPartonTestVector(exchange2, beta2, false);
  const double     mass2     = proc.state.lts.s * proc.state.lts.diff_xhard1 * proc.state.lts.diff_xhard2;
  const double     scale     = std::sqrt(mass2 / (base1 + base2).M2());

  gra::M4Vec q1;
  gra::M4Vec q2;
  REQUIRE(proc.BuildHardPair(proc.state.lts.diff_xhard1, proc.state.lts.diff_xhard2, q1, q2));
  REQUIRE(scale > 1.0);
  REQUIRE(scale * beta1 < 1.0);
  REQUIRE(scale * beta2 < 1.0);
  REQUIRE((q1 + q2).M2() == Approx(mass2).epsilon(1.0e-11));
  REQUIRE(q1.M2() == Approx(scale * scale * beta1 * t1).margin(1.0e-8));
  REQUIRE(q2.M2() == Approx(scale * scale * beta2 * t2).margin(1.0e-8));
  const std::array<int, 4> components = {0, 1, 2, 3};
  for (const auto &i : indices(components)) {
    REQUIRE(q1[i] == Approx(scale * base1[i]).epsilon(2.0e-12));
    REQUIRE(q2[i] == Approx(scale * base2[i]).epsilon(2.0e-12));
  }
  REQUIRE((q1 + q2).Rap() == Approx((base1 + base2).Rap()).epsilon(1.0e-12));

  const gra::M4Vec remnant1 = exchange1 - q1;
  const gra::M4Vec remnant2 = exchange2 - q2;
  REQUIRE(remnant1.E() > 0.0);
  REQUIRE(remnant2.E() > 0.0);
  REQUIRE(remnant1.M2() == Approx(-(scale - 1.0) * (1.0 - scale * beta1) * t1).margin(1.0e-8));
  REQUIRE(remnant2.M2() == Approx(-(scale - 1.0) * (1.0 - scale * beta2) * t2).margin(1.0e-8));
  REQUIRE(q1.LightconePos() / proc.state.lts.pbeam1.LightconePos() == Approx(scale * xi1 * beta1).epsilon(1.0e-12));
  REQUIRE(q2.LightconeNeg() / proc.state.lts.pbeam2.LightconeNeg() == Approx(scale * xi2 * beta2).epsilon(1.0e-12));
  REQUIRE(gra::math::CheckEMC(proc.state.lts.pbeam1 + proc.state.lts.pbeam2 -
                              (leading1 + remnant1 + leading2 + remnant2 + q1 + q2)));

  SECTION("the scalar map commutes with rotations and boosts") {
    gra::M4Vec mapped1 = base1;
    gra::M4Vec mapped2 = base2;
    REQUIRE(HardDiffractionTestProbe::ScaleDoubleDiffractivePair(mass2, beta1, beta2, mapped1, mapped2));

    gra::M4Vec transformed1 = base1;
    gra::M4Vec transformed2 = base2;
    RotateCovarianceMomentum(transformed1);
    RotateCovarianceMomentum(transformed2);
    BoostCovarianceMomentum(transformed1);
    BoostCovarianceMomentum(transformed2);
    REQUIRE(HardDiffractionTestProbe::ScaleDoubleDiffractivePair(mass2, beta1, beta2, transformed1, transformed2));

    RotateCovarianceMomentum(mapped1);
    RotateCovarianceMomentum(mapped2);
    BoostCovarianceMomentum(mapped1);
    BoostCovarianceMomentum(mapped2);
    for (const auto &i : indices(components)) {
      REQUIRE(transformed1[i] == Approx(mapped1[i]).epsilon(2.0e-12));
      REQUIRE(transformed2[i] == Approx(mapped2[i]).epsilon(2.0e-12));
    }
  }

  SECTION("beta one has no finite t local-remnant support") {
    proc.state.lts.diff_beta1  = 1.0;
    proc.state.lts.diff_xhard1 = xi1;
    gra::M4Vec rejected1(1.0, 2.0, 3.0, 4.0);
    gra::M4Vec rejected2(-1.0, -2.0, -3.0, 5.0);
    REQUIRE_FALSE(proc.BuildHardPair(proc.state.lts.diff_xhard1, proc.state.lts.diff_xhard2, rejected1, rejected2));
    REQUIRE(rejected1.E() == Approx(4.0));
    REQUIRE(rejected2.E() == Approx(5.0));
  }

  SECTION("the massless zero t limit is exactly collinear") {
    HardDiffractionTestProbe collinear;
    collinear.state.lts.pbeam1      = gra::M4Vec(0.0, 0.0, energy, energy);
    collinear.state.lts.pbeam2      = gra::M4Vec(0.0, 0.0, -energy, energy);
    collinear.state.lts.s           = (collinear.state.lts.pbeam1 + collinear.state.lts.pbeam2).M2();
    collinear.state.lts.hard_diff1  = true;
    collinear.state.lts.hard_diff2  = true;
    collinear.state.lts.diff_xi1    = xi1;
    collinear.state.lts.diff_xi2    = xi2;
    collinear.state.lts.diff_beta1  = beta1;
    collinear.state.lts.diff_beta2  = beta2;
    collinear.state.lts.diff_t1     = 0.0;
    collinear.state.lts.diff_t2     = 0.0;
    collinear.state.lts.diff_phi1   = phi1;
    collinear.state.lts.diff_phi2   = phi2;
    collinear.state.lts.diff_xhard1 = xi1 * beta1;
    collinear.state.lts.diff_xhard2 = xi2 * beta2;
    gra::M4Vec collinear1;
    gra::M4Vec collinear2;
    REQUIRE(collinear.BuildHardPair(collinear.state.lts.diff_xhard1, collinear.state.lts.diff_xhard2, collinear1,
                                    collinear2));
    REQUIRE(collinear1.M2() == Approx(0.0).margin(1.0e-10));
    REQUIRE(collinear2.M2() == Approx(0.0).margin(1.0e-10));
    REQUIRE((collinear1 + collinear2).M2() ==
            Approx(collinear.state.lts.s * xi1 * beta1 * xi2 * beta2).epsilon(1.0e-12));
  }
}

// Test exact forward-particle reconstruction used by the hard screening loop
TEST_CASE("Hard diffraction forward proton reconstruction stays on shell", "[HardPomeronPDF]") {
  const double     mass = gra::PDG::mp;
  const double     pz   = 6500.0;
  const gra::M4Vec beam(0.0, 0.0, pz, std::sqrt(pz * pz + mass * mass));
  const double     xi  = 0.03;
  const double     t   = -0.25;
  const double     phi = 0.7;

  gra::M4Vec leading;
  REQUIRE(gra::kinematics::BuildForwardParticleXiT(beam, xi, t, phi, true, leading));
  REQUIRE(leading.M2() == Approx(beam.M2()).margin(1.0e-7));
  REQUIRE((beam - leading).M2() == Approx(t).margin(1.0e-8));

  gra::M4Vec shifted;
  REQUIRE(gra::kinematics::BuildForwardParticle(beam, 1.0 - xi, leading.Px() + 0.2, leading.Py() - 0.1, true, shifted));
  REQUIRE(shifted.M2() == Approx(beam.M2()).margin(1.0e-7));
  REQUIRE(shifted.LightconePos() / beam.LightconePos() == Approx(1.0 - xi).epsilon(1.0e-12));
}

// Test shared hard-Pomeron wrapper and LHAPDF access from several threads
TEST_CASE("Hard Pomeron PDF can be evaluated from several threads", "[HardPomeronPDF]") {
  const std::string         modelfile = gra::ResolveModelDataFile("TUNE0", "GENERAL.json");
  gra::MHardPomeronPDFStore store;
  const auto                first  = store.GetHardPomeronPDF(modelfile);
  const auto                second = store.GetHardPomeronPDF(modelfile);
  REQUIRE(first == second);

  std::vector<std::shared_ptr<const gra::MHardPomeronPDF>> handles(4);
  std::atomic<int>                                         ok{0};

  std::vector<std::thread> workers;
  for (int i = 0; i < 4; ++i) {
    workers.emplace_back([i, &ok, &modelfile, &handles, &store]() {
      handles[i]         = store.GetHardPomeronPDF(modelfile);
      const double value = handles[i]->DiffractiveDensity(21, 0.01, 0.2, -0.1, 91.1876 * 91.1876);
      if (std::isfinite(value) && value >= 0.0) { ok.fetch_add(1); }
    });
  }

  for (auto &worker : workers) { worker.join(); }
  REQUIRE(ok.load() == 4);
  for (const auto &handle : handles) { REQUIRE(handle == first); }
}

// Test that the hard-Pomeron store does not return stale cards for reused paths
TEST_CASE("Hard Pomeron PDF store reloads changed model-card content", "[HardPomeronPDF]") {
  gra::MHardPomeronPDFStore store;
  const std::string         first_path = WriteHardPomeronTune("cache_reload", 2.0);
  const auto                first      = store.GetHardPomeronPDF(first_path);
  REQUIRE(first->FactorizationQ2(1.0) == Approx(2.0));

  const std::string second_path = WriteHardPomeronTune("cache_reload", 5.0);
  const auto        second      = store.GetHardPomeronPDF(second_path);
  REQUIRE(second != first);
  REQUIRE(second->FactorizationQ2(1.0) == Approx(5.0));
}

// Check PDF member and set replacement and preserve the previous state on failure
TEST_CASE("Hard Pomeron PDF reloads its LHAPDF selection atomically", "[HardPomeronPDF]") {
  const std::filesystem::path dir = "tmp/hard_pomeron_pdf_reload";
  std::filesystem::create_directories(dir);
  const std::string path    = (dir / "GENERAL.json").string();
  auto              card    = nlohmann::json::object();
  card["PARAM_HARDPOMERON"] = {{"DPDF_SET", "GKG18_DPDF_FitB_LO"},
                               {"DPDF_MEMBER", 0},
                               {"remnant_mass", gra::MHardPomeronPDFParam{}.remnant_mass}};
  // Write each selected PDF to the same card path
  const auto write_card = [&]() {
    std::ofstream out(path);
    REQUIRE(out.good());
    out << card.dump(2);
    out.close();
    REQUIRE(out.good());
  };
  write_card();
  gra::MHardPomeronPDF pdf(path);
  gra::MLHAPDFStore    store;
  for (const int member : {1, 0}) {
    card["PARAM_HARDPOMERON"]["DPDF_MEMBER"] = member;
    write_card();
    pdf.ReadParameters(path);
    const auto expected = store.GetPDF("GKG18_DPDF_FitB_LO", member);
    REQUIRE(pdf.PartonDensity(21, 0.2, 100.0) == Approx(expected->xfxQ2(21, 0.2, 100.0) / 0.2).epsilon(1e-12));
    REQUIRE(pdf.AlphaS(100.0) == Approx(expected->alphasQ2(100.0)).epsilon(1e-12));
  }
  card["PARAM_HARDPOMERON"]["DPDF_SET"] = "CT10nlo";
  write_card();
  pdf.ReadParameters(path);
  const auto   expected = store.GetPDF("CT10nlo", 0);
  const double density  = expected->xfxQ2(21, 0.2, 100.0) / 0.2;
  REQUIRE(pdf.PartonDensity(21, 0.2, 100.0) == Approx(density).epsilon(1e-12));
  REQUIRE(pdf.AlphaS(100.0) == Approx(expected->alphasQ2(100.0)).epsilon(1e-12));

  const double flux                     = pdf.Flux(0.01, -0.1);
  card["PARAM_HARDPOMERON"]["DPDF_SET"] = "../invalid";
  card["PARAM_HARDPOMERON"]["norm"]     = 2.0;
  write_card();
  REQUIRE_THROWS_AS(pdf.ReadParameters(path), std::invalid_argument);
  REQUIRE(pdf.PartonDensity(21, 0.2, 100.0) == Approx(density).epsilon(1e-12));
  REQUIRE(pdf.Flux(0.01, -0.1) == Approx(flux).epsilon(1e-12));
  REQUIRE(pdf.AlphaS(100.0) == Approx(expected->alphasQ2(100.0)).epsilon(1e-12));

  card["PARAM_HARDPOMERON"]["DPDF_SET"]    = "CT10nlo";
  card["PARAM_HARDPOMERON"]["DPDF_MEMBER"] = std::numeric_limits<int>::max();
  write_card();
  REQUIRE_THROWS_AS(pdf.ReadParameters(path), std::invalid_argument);
  REQUIRE(pdf.PartonDensity(21, 0.2, 100.0) == Approx(density).epsilon(1e-12));
  REQUIRE(pdf.Flux(0.01, -0.1) == Approx(flux).epsilon(1e-12));
  REQUIRE(pdf.AlphaS(100.0) == Approx(expected->alphasQ2(100.0)).epsilon(1e-12));
}

// Test strict support and normalization validation before loading LHAPDF
TEST_CASE("Hard Pomeron PDF parameters reject unphysical support", "[HardPomeronPDF]") {
  SECTION("xi upper limit") {
    gra::MHardPomeronPDFParam param;
    param.xi_range = {1.0e-4, 1.0};
    REQUIRE_THROWS(param.Validate());
  }

  SECTION("zero flux normalization") {
    gra::MHardPomeronPDFParam param;
    param.norm = 0.0;
    REQUIRE_THROWS(param.Validate());
  }

  SECTION("negative flux normalization") {
    gra::MHardPomeronPDFParam param;
    param.norm = -1.0;
    REQUIRE_THROWS(param.Validate());
  }
}

// Test that process-level hard diffraction samples one physical parton channel
TEST_CASE("Hard diffraction samples a channel after PDF convolution", "[HardPomeronPDF]") {
  gra::MODELPARAM = "TUNE0";
  const std::vector<gra::aux::OneCMD> syntax;
  const bool                          double_diff = GENERATE(true, false);
  const std::string                   process     = double_diff ? "IPIP[Z]<F>" : "IPp[Z]<F>";
  HardDiffractionTestProbe            proc(process, syntax);
  proc.SetHelicityConfig(gra::MModelTune::Load(modelfile));
  proc.SetInitialState({"p+", "p+"}, {6500.0, 6500.0});
  proc.SetScreening(false);
  const double beam_pz  = std::sqrt(gra::math::pow2(6500.0) - gra::math::pow2(gra::PDG::mp));
  proc.state.lts.pbeam1 = gra::M4Vec(0.0, 0.0, beam_pz, 6500.0);
  proc.state.lts.pbeam2 = gra::M4Vec(0.0, 0.0, -beam_pz, 6500.0);
  proc.state.lts.s      = (proc.state.lts.pbeam1 + proc.state.lts.pbeam2).M2();
  proc.state.lts.sqrt_s = std::sqrt(proc.state.lts.s);

  proc.state.lts.process.root_decay_mode = gra::RootDecayMode::Physical;
  proc.state.lts.LHAPDFSET               = double_diff ? "null" : "MMHT2014lo68cl";
  proc.state.lts.hard_diff1              = true;
  proc.state.lts.hard_diff2              = double_diff;
  proc.state.lts.diff_xi1                = 0.01;
  proc.state.lts.diff_beta1              = 0.2;
  proc.state.lts.diff_t1                 = -0.1;
  proc.state.lts.diff_xi2                = 0.012;
  proc.state.lts.diff_beta2              = 0.18;
  proc.state.lts.diff_t2                 = -0.12;
  proc.state.lts.diff_xhard1             = proc.state.lts.diff_xi1 * proc.state.lts.diff_beta1;
  proc.state.lts.diff_xhard2             = proc.state.lts.diff_xi2 * proc.state.lts.diff_beta2;

  proc.SetModelTune(gra::MModelTune::Load(modelfile));
  proc.FinalizeProcessConfiguration();

  proc.state.lts.diff_phi1 = 0.3;
  proc.state.lts.diff_phi2 = proc.state.lts.diff_phi1 + gra::math::PI;
  REQUIRE(
      proc.BuildHardPair(proc.state.lts.diff_xhard1, proc.state.lts.diff_xhard2, proc.state.lts.q1, proc.state.lts.q2));
  proc.state.lts.pfinal.assign(6, gra::M4Vec(0.0, 0.0, 0.0, 0.0));
  proc.state.lts.pfinal[0]     = proc.state.lts.q1 + proc.state.lts.q2;
  proc.state.lts.pfinal[1]     = proc.state.lts.pbeam1 - proc.state.lts.q1;
  proc.state.lts.pfinal[2]     = proc.state.lts.pbeam2 - proc.state.lts.q2;
  proc.state.lts.forward_mass2 = {proc.state.lts.pfinal[1].M2(), proc.state.lts.pfinal[2].M2()};
  const double hard_mass       = proc.state.lts.pfinal[0].M();
  proc.state.gcuts.M_min       = 0.9 * hard_mass;
  proc.state.gcuts.M_max       = 1.1 * hard_mass;
  const double muon_mass = HardMass(13);
  const double momentum = std::sqrt(gra::math::pow2(hard_mass / 2.0) - muon_mass * muon_mass);
  gra::M4Vec muon1(0.6 * momentum, 0.0, 0.8 * momentum, 0.5 * hard_mass);
  gra::M4Vec   muon2 = -muon1;
  muon2.SetE(muon1.E());
  gra::kinematics::LorentzBoost(proc.state.lts.pfinal[0], hard_mass, muon1, 1);
  gra::kinematics::LorentzBoost(proc.state.lts.pfinal[0], hard_mass, muon2, 1);
  proc.state.lts.decaytree.push_back(StableBranch(-13, muon1, muon_mass));
  proc.state.lts.decaytree.push_back(StableBranch(13, muon2, muon_mass));

  auto eikonal = BuildTestEikonalWithInitialState({}, "hard_diffraction_dense_spin", proc.GetInitialState(), 0.25, 1, 4,
                                                  4, proc.state.lts.s, 0.0, 0.65);
  eikonal.Numerics.LOOP.radial_integrator  = "1/3";
  eikonal.Numerics.LOOP.azimuth_integrator = "Trap";
  eikonal.Numerics.LOOP.r_min              = 0.0;
  eikonal.Numerics.LOOP.r_max              = 0.2;
  eikonal.Numerics.LOOP.radial_intervals   = 2;
  eikonal.Numerics.LOOP.azimuth_nodes      = 3;
  eikonal.InitLoopWeightMatrix();
  REQUIRE_FALSE(eikonal.GetLoopConst(proc.state.lts.s).physical_screening_is_scalar);
  proc.SetModelTune(eikonal.ModelTuneHandle());
  proc.InitializeProcessAmplitude();
  proc.FinalizeProcessConfiguration();

  gra::MEventWeightState aux;
  const double           amp2 = proc.ProbeAmp2(aux);
  REQUIRE(std::isfinite(amp2));
  REQUIRE(amp2 > 0.0);
  REQUIRE(aux.amplitude_ok);
  REQUIRE(aux.kinematics_ok);
  REQUIRE(proc.state.lts.hamp.size() > 1);
  REQUIRE(HampNormSum(proc.state.lts) == Approx(amp2).epsilon(1e-10));
  REQUIRE(HasComplexHamp(proc.state.lts));
  REQUIRE(proc.state.lts.id1 != 0);
  REQUIRE(proc.state.lts.id2 != 0);
  REQUIRE(proc.state.lts.id1 == -proc.state.lts.id2);
  REQUIRE(std::abs(proc.state.lts.id1) >= 1);
  REQUIRE(std::abs(proc.state.lts.id1) <= 5);
  REQUIRE(proc.state.lts.pdf_xf1 > 0.0);
  REQUIRE(proc.state.lts.pdf_xf2 > 0.0);
  REQUIRE(proc.state.lts.muF > 0.0);
  REQUIRE(proc.state.lts.muR == Approx(proc.state.lts.muF));
  REQUIRE(proc.state.lts.scalup == Approx(proc.state.lts.muF));
  REQUIRE_FALSE(proc.state.lts.hard_color_flows.empty());

  const std::size_t born_component_count = proc.state.lts.hamp.size();
  proc.SetEikonal(eikonal);
  proc.SetScreening(true);
  gra::MEventWeightState screened_aux;
  const double           screened_amp2 = proc.ProbeAmp2(screened_aux);
  REQUIRE(screened_aux.Valid());
  REQUIRE(screened_amp2 > 0.0);
  REQUIRE(proc.state.lts.hamp.size() == 16 * born_component_count);
  REQUIRE(proc.state.lts.hamp.metadata.spin_basis == gra::ScreeningSpinBasis::ProtonIdentity);
  REQUIRE(proc.state.lts.hamp.metadata.spin_rows == 1);
  REQUIRE(proc.state.lts.hamp.metadata.amplitude_normalization == Approx(1.0));
  // Expanding the proton identity already supplies 1/2 in each amplitude
  REQUIRE(HampNormSum(proc.state.lts) == Approx(screened_amp2).epsilon(1.0e-10));
  REQUIRE(HasComplexHamp(proc.state.lts));
  REQUIRE(proc.state.lts.id1 == -proc.state.lts.id2);
  REQUIRE(std::abs(proc.state.lts.id1) >= 1);
  REQUIRE(std::abs(proc.state.lts.id1) <= 5);
  REQUIRE(proc.state.lts.pdf_xf1 > 0.0);
  REQUIRE(proc.state.lts.pdf_xf2 > 0.0);
  const double sampled_mass2 = proc.state.lts.s * proc.state.lts.diff_xhard1 * proc.state.lts.diff_xhard2;
  REQUIRE(proc.state.lts.pfinal[0].M2() == Approx(sampled_mass2).epsilon(1.0e-11));
  for (const auto *forward : {&proc.state.lts.decayforward1, &proc.state.lts.decayforward2}) {
    REQUIRE(forward->legs.size() == 2);
    REQUIRE(forward->legs[1].p4.E() > 0.0);
    REQUIRE(forward->legs[1].p4.M2() >= -1.0e-8);
    REQUIRE(gra::math::CheckEMC(forward->p4 - forward->legs[0].p4 - forward->legs[1].p4));
  }

  HepMC3::GenEvent evt(HepMC3::Units::GEV, HepMC3::Units::MM);
  REQUIRE(proc.FinalizeProbe());
  REQUIRE(proc.EventRecord(evt));
  const auto pdf_info = evt.pdf_info();
  REQUIRE(pdf_info != nullptr);
  REQUIRE(pdf_info->is_valid());
  REQUIRE(pdf_info->x[0] == Approx(proc.state.lts.diff_xhard1));
  REQUIRE(pdf_info->x[1] == Approx(proc.state.lts.diff_xhard2));
  const std::string diff_comment = gra::BuildLHEDiffractionComment(evt);
  REQUIRE(diff_comment.find(double_diff ? "side_mask=3" : "side_mask=1") != std::string::npos);
  REQUIRE(diff_comment.find("side1_xi=") != std::string::npos);
  REQUIRE((diff_comment.find("side2_xi=") != std::string::npos) == double_diff);
  for (const auto &vertex : evt.vertices()) {
    gra::M4Vec residual;
    for (const auto &particle : vertex->particles_in()) { residual += HepMCMomentum(particle->momentum()); }
    for (const auto &particle : vertex->particles_out()) { residual -= HepMCMomentum(particle->momentum()); }
    REQUIRE(gra::math::CheckEMC(residual));
  }
}

// Test that unsupported N-star steering is rejected before integration
TEST_CASE("Hard and collinear photon processes reject N-star steering", "[HardPomeronPDF]") {
  for (const int excitation : {1, 2}) {
    gra::MHardDiffraction hard;
    hard.SetExcitation(excitation);
    REQUIRE_THROWS(hard.FinalizeProcessConfiguration());

    gra::MCollinear collinear;
    collinear.SetExcitation(excitation);
    REQUIRE_THROWS(collinear.FinalizeProcessConfiguration());

    gra::MQuasiElastic quasielastic;
    quasielastic.SetExcitation(excitation);
    REQUIRE_THROWS(quasielastic.FinalizeProcessConfiguration());
  }
}

// Check that collinear photon processes reject unsupported screening
TEST_CASE("Collinear photon processes reject unsupported screening", "[HardPomeronPDF]") {
  gra::MCollinear collinear;
  collinear.SetScreening(true);
  REQUIRE_THROWS_AS(collinear.FinalizeProcessConfiguration(), std::invalid_argument);
}

// Test process-level MG5 color-flow sampling for current and future jet
// multiplicities
TEST_CASE("Hard diffraction maps exact MG5 color candidates onto proton remnants", "[HardPomeronPDF][color]") {
  gra::MODELPARAM = "TUNE0";
  const std::vector<gra::aux::OneCMD> syntax;
  const auto                          require_closed_record = [](HardDiffractionTestProbe &probe) {
    HepMC3::GenEvent evt(HepMC3::Units::GEV, HepMC3::Units::MM);
    REQUIRE(probe.FinalizeProbe());
    REQUIRE(probe.EventRecord(evt));
    std::map<int, int> balance;
    int                colored_with_flow = 0;
    for (const auto &particle : evt.particles()) {
      if (!IsColoredParton(particle->pid())) { continue; }
      const int flow1 = ReadIntAttribute(particle, "flow1");
      const int flow2 = ReadIntAttribute(particle, "flow2");
      if (flow1 != 0 || flow2 != 0) { ++colored_with_flow; }
      if (flow1 != 0) { ++balance[flow1]; }
      if (flow2 != 0) { --balance[flow2]; }
    }
    const auto rows = gra::BuildLHERows(evt);
    REQUIRE(rows[0].pid == 2212);
    REQUIRE(rows[1].pid == 2212);
    gra::M4Vec final_sum;
    for (const auto &row : rows) {
      if (row.status != 1) { continue; }
      final_sum += gra::M4Vec(row.momentum[0], row.momentum[1], row.momentum[2], row.momentum[3]);
      REQUIRE(row.colors.first == ReadIntAttribute(row.particle, "flow1"));
      REQUIRE(row.colors.second == ReadIntAttribute(row.particle, "flow2"));
    }
    REQUIRE(gra::math::CheckEMC(probe.state.lts.pbeam1 + probe.state.lts.pbeam2 - final_sum));
    // Check that the ISR view retains hard momentum conservation and the selected color topology
    std::istringstream                hard_input(gra::BuildLHEDiffractionComment(evt));
    std::string                       line;
    gra::M4Vec                        hard_residual;
    std::map<int, std::array<int, 2>> hard_colors;
    int                               incoming = 0;
    int                               final    = 0;
    while (std::getline(hard_input, line)) {
      if (line.rfind("# graniitti_hard ", 0) != 0) { continue; }
      std::istringstream row(line.substr(16));
      int                pdg, status, mother1, mother2, flow1, flow2;
      double             px, py, pz, energy, mass, tau, spin;
      REQUIRE(static_cast<bool>(row >> pdg >> status >> mother1 >> mother2 >> flow1 >> flow2 >> px >> py >> pz >>
                                energy >> mass >> tau >> spin));
      const gra::M4Vec p(px, py, pz, energy);
      if (status == -1) {
        hard_residual += p;
        REQUIRE(pdg == evt.pdf_info()->parton_id[incoming]);
        ++incoming;
        std::swap(flow1, flow2);
      } else if (status == 1) {
        hard_residual -= p;
        ++final;
      } else {
        continue;
      }
      if (flow1 != 0) { ++hard_colors[flow1][0]; }
      if (flow2 != 0) { ++hard_colors[flow2][1]; }
    }
    REQUIRE(incoming == 2);
    REQUIRE(final >= 2);
    REQUIRE(gra::math::CheckEMC(hard_residual));
    for (const auto &[tag, count] : hard_colors) {
      REQUIRE(tag >= 501);
      REQUIRE(count[0] == 1);
      REQUIRE(count[1] == 1);
    }
    REQUIRE(colored_with_flow > 0);
    for (const auto &[tag, value] : balance) {
      REQUIRE(tag >= 501);
      REQUIRE(value == 0);
    }
  };

  {
    const std::vector<std::pair<std::string, gra::LORENTZSCALAR>> channels = {{"Z", PPZReferenceKinematics(0.118)},
                                                                              {"Zj", PPZjAllChannelKinematics(0.118)},
                                                                              {"jj", PPJJReferenceKinematics(0.118)},
                                                                              {"W", PPWReferenceKinematics(0.118)}};
    for (const auto &[channel, reference] : channels) {
      INFO(channel);
      HardDiffractionTestProbe probe("IPp[" + channel + "]<F>", syntax);
      probe.SetHelicityConfig(gra::MModelTune::Load(modelfile));
      // Use the common five-flavour PDF evolution interval for the convoluted color test
      ConfigureHardColorProbe(probe, false, reference.decaytree, 100.0);
      probe.SetModelTune(gra::MModelTune::Load(modelfile));
      probe.InitializeProcessAmplitude();
      REQUIRE_NOTHROW(probe.FinalizeProcessConfiguration());
      gra::MEventWeightState aux;
      REQUIRE(probe.ProbeAmp2(aux) > 0.0);
      REQUIRE(aux.Valid());
      REQUIRE_FALSE(probe.state.lts.hard_color_flows.empty());
      require_closed_record(probe);
      probe.state.gcuts.M_max = 320.0;
      REQUIRE_NOTHROW(probe.FinalizeProcessConfiguration());
    }
  }

  {
    std::string              process = "IPIP[jj]<F>";
    HardDiffractionTestProbe probe(process, syntax);
    probe.SetHelicityConfig(gra::MModelTune::Load(modelfile));
    ConfigureHardColorProbe(probe, true, PPJJReferenceKinematics(0.118).decaytree);
    probe.SetModelTune(gra::MModelTune::Load(modelfile));
    probe.InitializeProcessAmplitude();
    probe.FinalizeProcessConfiguration();
    gra::MEventWeightState aux;
    REQUIRE(probe.ProbeAmp2(aux) > 0.0);
    REQUIRE(aux.Valid());
    REQUIRE_FALSE(probe.state.lts.hard_color_flows.empty());
    require_closed_record(probe);
  }
}

// Test exact generated color association with raw JAMP rows and mirrors
TEST_CASE("MG5 color structure preserves raw JAMP and incoming mirror order", "[HardPomeronPDF][MG5][color]") {
  const std::array<std::vector<std::vector<int>>, 2> expected = {{
      {{1, 2, 3, 0, 3, 2, 1, 0}, {3, 2, 2, 0, 3, 1, 1, 0}},
      {{3, 0, 1, 2, 3, 2, 1, 0}, {2, 0, 3, 2, 3, 1, 1, 0}},
  }};

  gra::AMP_MG5_pp_jj amp;
  for (const bool mirror_initial : {false, true}) {
    gra::LORENTZSCALAR lts        = PPJJReversedFlowKinematics(0.118, mirror_initial);
    const auto         evaluation = EvaluateHard(amp, lts, 0.118);
    REQUIRE(evaluation.amp2 > 0.0);
    const auto &prepared = evaluation.color_flows;
    REQUIRE(prepared.size() == 2);
    for (std::size_t flow = 0; flow < prepared.size(); ++flow) {
      INFO("mirror=" << mirror_initial << " raw JAMP=" << flow);
      REQUIRE(FlattenExternalColorFlow(prepared[flow].external) == expected[mirror_initial ? 1 : 0][flow]);
      REQUIRE_FALSE(prepared[flow].amplitudes.empty());
    }
  }
}

// Test deterministic subprocess sum rejection of incomplete color structure
TEST_CASE("MG5 subprocess sum rejects unpaired raw JAMP rows", "[HardPomeronPDF][MG5][color]") {
  gra::mg5::Channel                              channel;
  std::vector<std::vector<std::complex<double>>> two_rows(2);
  std::vector<gra::mg5helas::HardColorFlow>      prepared;
  REQUIRE_FALSE(gra::mg5::PairHardColorFlows(channel, two_rows, prepared));

  channel.external_color_representations = {3, -3, 1, 1};
  channel.external_color_flows           = {{{1, 0}, {0, 1}, {0, 0}, {0, 0}}};
  REQUIRE_FALSE(gra::mg5::PairHardColorFlows(channel, two_rows, prepared));
}

// Test the generated neutral current normalization and parity even angle term
TEST_CASE("Generated neutral-current coefficient matches the LO formula", "[HardPomeronPDF][MG5][normalization]") {
  const bool z_only = GENERATE(false, true);
  gra::AMP_MG5_pp_z amplitude;
  MasslessMuon(amplitude);
  for (const double mass : {80.0, 91.188, 100.0}) {
    for (const int pid : {1, 2, 5}) {
      const double sigma = NeutralCurrentSigmaLO(mass, pid, z_only);
      REQUIRE(sigma > 0.0);
      for (const double cosine : {0.0, 0.37, 0.82}) {
        gra::LORENTZSCALAR forward  = NeutralCurrentPoint(mass, pid, cosine, z_only);
        gra::LORENTZSCALAR backward = NeutralCurrentPoint(mass, pid, -cosine, z_only);
        const double       amp2_sum = HardAmp2(amplitude, forward, 0.118) + HardAmp2(amplitude, backward, 0.118);
        const double       expected = 24.0 * gra::math::PI * mass * mass * sigma * (1.0 + cosine * cosine);
        CAPTURE(z_only, mass, pid, cosine, amp2_sum, expected);
        REQUIRE(amp2_sum == Approx(expected).epsilon(2.0e-11));
      }
    }
  }
}

// Check the explicit Z decay amplitude under boosts, rotations and screening recoil
TEST_CASE("Explicit Z decay preserves hard amplitude covariance", "[HardPomeronPDF][MG5][covariance]") {
  gra::AMP_MG5_pp_z amplitude;
  MasslessMuon(amplitude);
  const auto point = NeutralCurrentPoint(91.188, 2, 0.37, true);
  const auto evaluate = [&](gra::LORENTZSCALAR &lts) { return HardAmp2(amplitude, lts, 0.118); };
  RequireHardScalarCovariance(amplitude, point, 0.118, "explicit Z decay");
  RequireRotatedScreeningRingSection(point, evaluate, "explicit Z decay");
  RequireHelicitySectionContinuity(point, evaluate, "explicit Z decay");
}

// Test direct MG5 wrappers at fixed standalone MG5_aMC@NLO benchmark points
// Reference convention: Alwall et al., arXiv:1405.0301, standalone C++ export
TEST_CASE("Diffractive Z MG5 wrappers match standalone reference values", "[HardPomeronPDF]") {
  {
    gra::AMP_MG5_yy_zjj amp;
    std::vector<gra::mg5::Subprocess<MG5_YY_ZJJ::ProcessBase>> standalone;
    MG5_YY_ZJJ::BuildSubprocesses(standalone,
        gra::aux::ResolveProjectPath("MG5cards/Photon/MG5_YY_ZJJ/param_card.dat"));
    for (auto &subprocess : standalone) { subprocess.process->setAlphaQEDZero(true); }
    double flavour_sum = 0.0;
    double standalone_sum = 0.0;
    for (const int flavour : {1, 2, 3, 4, 5}) {
      gra::LORENTZSCALAR lts  = YYZjjReferenceKinematics(flavour);
      const double       amp2 = PhotonZjjAmp2(amp, lts, false);
      flavour_sum += amp2;
      auto reference = lts;
      reference.id1 = 22;
      reference.id2 = 22;
      const double scalar = GeneratedFullDenominatorAmp2(standalone, reference, 0.0);
      REQUIRE(scalar > 0.0);
      REQUIRE(amp2 == Approx(scalar).epsilon(1e-12));
      standalone_sum += scalar;
      RequireUnaveragedHelicityHamp(lts, amp2);
    }
    RequireAmplitudeReference("gamma gamma -> mu+ mu- plus summed q qbar", flavour_sum, standalone_sum);

    gra::LORENTZSCALAR explicit_u      = YYZjjReferenceKinematics(2);
    const double       explicit_u_amp2 = PhotonZjjAmp2(amp, explicit_u, false);
    REQUIRE(explicit_u_amp2 > 0.0);
    REQUIRE(explicit_u_amp2 < flavour_sum);
    // Equal massless up-type flavours contribute equally at fixed kinematics
    auto explicit_c = YYZjjReferenceKinematics(4);
    REQUIRE(PhotonZjjAmp2(amp, explicit_c, false) == Approx(explicit_u_amp2).epsilon(1.0e-12));
  }

  {
    gra::LORENTZSCALAR lts = PPZReferenceKinematics(0.118, 0.0);
    gra::AMP_MG5_pp_z  amp;
    MasslessMuon(amp);
    const auto         evaluation = EvaluateHard(amp, lts, 0.118);
    const double       amp2       = evaluation.amp2;
    RequireAmplitudeReference("u ubar -> mu+ mu-", amp2, 5.62803725932748008e-04);
    RequireComplexHelicityHamp(lts, amp2);
    REQUIRE(amp.SubprocessCount() == 18);
    REQUIRE(evaluation.contributing_subprocesses == 1);
    REQUIRE(evaluation.color_flows.size() == 1);

    for (const int id1 : {1, -1, 2, -2, 3, -3, 4, -4, 5, -5}) {
      gra::LORENTZSCALAR flavour_lts = PPZReferenceKinematics(0.118, 0.0);
      flavour_lts.id1                = id1;
      flavour_lts.id2                = -id1;
      const double flavour_amp2      = HardAmp2(amp, flavour_lts, 0.118);
      REQUIRE(flavour_amp2 > 0.0);
      RequireComplexHelicityHamp(flavour_lts, flavour_amp2);
    }
  }

  {
    gra::LORENTZSCALAR lts = PPZjReferenceKinematics(0.118, 0.0);
    gra::AMP_MG5_pp_zj amp;
    MasslessMuon(amp);
    const auto         evaluation = EvaluateHard(amp, lts, 0.118);
    const double       amp2       = evaluation.amp2;
    RequireAmplitudeReference("g u -> mu+ mu- u", amp2, 1.14477244095003320e-07);
    RequireComplexHelicityHamp(lts, amp2);
    REQUIRE(amp.SubprocessCount() == 27);
    REQUIRE(evaluation.contributing_subprocesses == 1);
    REQUIRE(evaluation.color_flows.size() == 1);
  }

  {
    gra::LORENTZSCALAR lts = PPJJReferenceKinematics(0.118);
    gra::AMP_MG5_pp_jj amp;
    const auto         evaluation = EvaluateHard(amp, lts, 0.118);
    const double       amp2       = evaluation.amp2;
    REQUIRE(amp2 > 0.0);
    RequireComplexHelicityHamp(lts, amp2, false);
    REQUIRE(amp.SubprocessCount() == 25);
    REQUIRE(evaluation.contributing_subprocesses == 1);
    REQUIRE(evaluation.color_flows.size() == 2);
    // Check flavour multiplicity independently of the generator's random sequence
    auto charm = lts;
    charm.decaytree[0].p.pdg = 4;
    charm.decaytree[1].p.pdg = -4;
    REQUIRE(HardAmp2(amp, charm, 0.118) == Approx(amp2).epsilon(1.0e-12));
  }

  {
    gra::LORENTZSCALAR lts = PPWReferenceKinematics(0.118);
    gra::AMP_MG5_pp_w  amp;
    const auto         evaluation = EvaluateHard(amp, lts, 0.118);
    REQUIRE(evaluation.amp2 > 0.0);
    RequireComplexHelicityHamp(lts, evaluation.amp2, false);
    REQUIRE(amp.SubprocessCount() == 6);
    REQUIRE(evaluation.contributing_subprocesses == 1);
    REQUIRE(evaluation.color_flows.size() == 1);
  }
}

// Reject massive incoming model states while retaining massive outgoing particles
TEST_CASE("MG5 families require model masses compatible with incoming kinematics", "[HardPomeronPDF][MG5][mass-shell]") {
  for (const std::string family : {"MG5_PP_Z", "MG5_PP_ZJ", "MG5_PP_JJ"}) {
    CAPTURE(family);
    auto amplitude = gra::CreatePartonMG5Process(family);
    SLHAReader card(gra::aux::ResolveProjectPath("MG5cards/Parton/" + family + "/param_card.dat"));
    REQUIRE_NOTHROW(amplitude->InitParameters(card));
    card.set_block_entry("mass", 5, 4.7);
    REQUIRE_THROWS_AS(amplitude->InitParameters(card), std::invalid_argument);
  }

  const std::string path = gra::aux::ResolveProjectPath("MG5cards/Parton/MG5_PP_Z/param_card.dat");
  std::vector<gra::mg5::Subprocess<MG5_PP_Z::ProcessBase>> subprocesses;
  MG5_PP_Z::BuildSubprocesses(subprocesses, path);
  SLHAReader massive(path);
  massive.set_block_entry("mass", 5, 4.7);
  for (auto &subprocess : subprocesses) { subprocess.process->InitParameters(massive); }
  REQUIRE_THROWS_AS(gra::mg5::SubprocessSum<MG5_PP_Z::ProcessBase>(std::move(subprocesses)), std::invalid_argument);

  gra::AMP_MG5_yy_jj photons;
  SLHAReader card(gra::aux::ResolveProjectPath("MG5cards/Photon/MG5_YY_JJ/param_card.dat"));
  card.set_block_entry("mass", 5, 4.7);
  REQUIRE_NOTHROW(photons.InitParameters(card));
  REQUIRE(photons.Particles().at(5).mass == Approx(4.7));
}

// Test the massless five-flavour hard coefficient against both incoming
// projections
TEST_CASE("Hard Pomeron bottom coefficients and projected spinors are massless", "[HardPomeronPDF][MG5][mass-shell]") {
  const std::string zj_card = gra::aux::ResolveProjectPath("MG5cards/Parton/MG5_PP_ZJ/param_card.dat");
  std::vector<gra::mg5::Subprocess<MG5_PP_ZJ::ProcessBase>> zj_subprocesses;
  MG5_PP_ZJ::BuildSubprocesses(zj_subprocesses, zj_card);
  const auto *zj_process = FindChannelProcess(zj_subprocesses, std::array<int, 2>{21, 5}, std::vector<int>{-13, 13, 5});
  REQUIRE(zj_process != nullptr);
  const auto &zj_masses = zj_process->getMasses();
  REQUIRE(zj_masses.size() == 5);
  REQUIRE(zj_masses[1] == Approx(0.0));
  REQUIRE(zj_masses[4] == Approx(0.0));

  const std::string jj_card = gra::aux::ResolveProjectPath("MG5cards/Parton/MG5_PP_JJ/param_card.dat");
  std::vector<gra::mg5::Subprocess<MG5_PP_JJ::ProcessBase>> jj_subprocesses;
  MG5_PP_JJ::BuildSubprocesses(jj_subprocesses, jj_card);
  const auto *jj_process = FindChannelProcess(jj_subprocesses, std::array<int, 2>{21, 5}, std::vector<int>{21, 5});
  REQUIRE(jj_process != nullptr);
  const auto &jj_masses = jj_process->getMasses();
  REQUIRE(jj_masses.size() == 4);
  REQUIRE(jj_masses[1] == Approx(0.0));
  REQUIRE(jj_masses[3] == Approx(0.0));

  gra::LORENTZSCALAR reference                 = PPZjAllChannelKinematics(0.118);
  reference.id1                                = 21;
  reference.id2                                = 5;
  reference.decaytree[1].p.pdg                 = 5;
  const std::vector<gra::M4Vec> physical_final = {reference.decaytree[0].legs[0].p4, reference.decaytree[0].legs[1].p4,
                                                  reference.decaytree[1].p4};

  std::vector<gra::M4Vec> projected_final = physical_final;
  gra::M4Vec              p1;
  gra::M4Vec              p2;
  REQUIRE(gra::mg5helas::PrepareOnShellKinematics(reference, projected_final, p1, p2));
  REQUIRE(p1.M2() == Approx(zj_masses[0] * zj_masses[0]).margin(1e-10));
  REQUIRE(p2.M2() == Approx(zj_masses[1] * zj_masses[1]).margin(1e-10));
}

// Test every generated subprocess against the original scalar MG5 matrix
// element
TEST_CASE("All generated MG5 hard-process components preserve scalar normalization", "[HardPomeronPDF][MG5]") {
  {
    std::vector<gra::mg5::Subprocess<MG5_PP_Z::ProcessBase>> subprocesses;
    MG5_PP_Z::BuildSubprocesses(subprocesses, gra::aux::ResolveProjectPath("MG5cards/Parton/MG5_PP_Z/param_card.dat"));
    RequireGeneratedChannelNorms(subprocesses, PPZReferenceKinematics(0.118), 0.118, "PP_Z");
  }
  {
    std::vector<gra::mg5::Subprocess<MG5_PP_ZJ::ProcessBase>> subprocesses;
    MG5_PP_ZJ::BuildSubprocesses(subprocesses,
                                 gra::aux::ResolveProjectPath("MG5cards/Parton/MG5_PP_ZJ/param_card.dat"));
    RequireGeneratedChannelNorms(subprocesses, PPZjAllChannelKinematics(0.118), 0.118, "PP_ZJ");
  }
  {
    std::vector<gra::mg5::Subprocess<MG5_PP_JJ::ProcessBase>> subprocesses;
    MG5_PP_JJ::BuildSubprocesses(subprocesses,
                                 gra::aux::ResolveProjectPath("MG5cards/Parton/MG5_PP_JJ/param_card.dat"));
    RequireGeneratedChannelNorms(subprocesses, PPJJReferenceKinematics(0.118), 0.118, "PP_JJ");
  }
  {
    std::vector<gra::mg5::Subprocess<MG5_PP_W::ProcessBase>> subprocesses;
    MG5_PP_W::BuildSubprocesses(subprocesses, gra::aux::ResolveProjectPath("MG5cards/Parton/MG5_PP_W/param_card.dat"));
    RequireGeneratedChannelNorms(subprocesses, PPWReferenceKinematics(0.118), 0.118, "PP_W");
  }
  {
    std::vector<gra::mg5::Subprocess<MG5_YY_JJ::ProcessBase>> subprocesses;
    MG5_YY_JJ::BuildSubprocesses(subprocesses,
                                 gra::aux::ResolveProjectPath("MG5cards/Photon/MG5_YY_JJ/param_card.dat"));
    RequireGeneratedChannelNorms(subprocesses, YYJJReferenceKinematics(2), 0.0, "YY_JJ");
  }
  {
    std::vector<gra::mg5::Subprocess<MG5_YY_ZJJ::ProcessBase>> subprocesses;
    MG5_YY_ZJJ::BuildSubprocesses(subprocesses,
                                  gra::aux::ResolveProjectPath("MG5cards/Photon/MG5_YY_ZJJ/param_card.dat"));
    RequireGeneratedChannelNorms(subprocesses, YYZjjReferenceKinematics(2), 0.0, "YY_ZJJ");
  }
}

// Test the split between the generated denominator and the global phase-space
// symmetry factor
TEST_CASE("MG5 hard-process channel sum removes only final-state symmetry", "[HardPomeronPDF][MG5][symmetry]") {
  std::vector<gra::mg5::Subprocess<MG5_PP_JJ::ProcessBase>> channels;
  MG5_PP_JJ::BuildSubprocesses(channels, gra::aux::ResolveProjectPath("MG5cards/Parton/MG5_PP_JJ/param_card.dat"));

  gra::LORENTZSCALAR distinct      = PPJJReferenceKinematics(0.118);
  distinct.id1                     = 2;
  distinct.id2                     = -2;
  distinct.decaytree[0].p.pdg      = 5;
  distinct.decaytree[1].p.pdg      = -5;
  const double       distinct_full = GeneratedFullDenominatorAmp2(channels, distinct, 0.118);
  gra::AMP_MG5_pp_jj wrapper;
  const double       distinct_sum = HardAmp2(wrapper, distinct, 0.118);

  gra::LORENTZSCALAR repeated = PPJJReferenceKinematics(0.118);
  repeated.id1                = 2;
  repeated.id2                = -2;
  repeated.decaytree[0].p.pdg = 21;
  repeated.decaytree[1].p.pdg = 21;
  const double repeated_full  = GeneratedFullDenominatorAmp2(channels, repeated, 0.118);
  const double repeated_sum   = HardAmp2(wrapper, repeated, 0.118);

  REQUIRE(gra::math::IsExactEqual(gra::mg5helas::FinalStateSymmetryFactor(distinct.decaytree), 1.0));
  REQUIRE(gra::math::IsExactEqual(gra::mg5helas::FinalStateSymmetryFactor(repeated.decaytree), 2.0));
  REQUIRE(distinct_full > 0.0);
  REQUIRE(repeated_full > 0.0);
  REQUIRE(distinct_sum == Approx(distinct_full).epsilon(1e-10));
  REQUIRE(repeated_sum == Approx(2.0 * repeated_full).epsilon(1e-10));
  RequireComplexHelicityHamp(distinct, distinct_sum, false);
  RequireComplexHelicityHamp(repeated, repeated_sum, false);
}

// Test that the MG5 helicity section is continuous around either incoming beam
// axis
TEST_CASE("Diffractive Z MG5 helicity phases use a fixed beam-axis section", "[HardPomeronPDF][screening]") {
  gra::AMP_MG5_pp_z amp;

  for (const int id1 : {1, -1, 2, -2, 5, -5}) {
    gra::LORENTZSCALAR reference = PPZReferenceKinematics(0.118);
    reference.id1                = id1;
    reference.id2                = -id1;
    RequireHelicitySectionContinuity(
        reference, [&](gra::LORENTZSCALAR &lts) { return HardAmp2(amp, lts, 0.118); },
        "q qbar -> Z, incoming id1=" + std::to_string(id1));
  }

  gra::AMP_MG5_pp_zj       zj;
  const auto               zj_evaluate = [&](gra::LORENTZSCALAR &lts) { return HardAmp2(zj, lts, 0.118); };
  const gra::LORENTZSCALAR gq          = PPZjAllChannelKinematics(0.118);
  RequireHelicitySectionContinuity(gq, zj_evaluate, "g u -> Z u");
  RequireRotatedScreeningRingSection(gq, zj_evaluate, "g u -> Z u");

  gra::LORENTZSCALAR qg = gq;
  qg.id1                = 2;
  qg.id2                = 21;
  RequireHelicitySectionContinuity(qg, zj_evaluate, "u g -> Z u");
  RequireRotatedScreeningRingSection(qg, zj_evaluate, "u g -> Z u");

  gra::LORENTZSCALAR qbar_g = qg;
  qbar_g.id1                = -2;
  qbar_g.decaytree[1].p.pdg = -2;
  RequireHelicitySectionContinuity(qbar_g, zj_evaluate, "ubar g -> Z ubar");
  RequireRotatedScreeningRingSection(qbar_g, zj_evaluate, "ubar g -> Z ubar");

  gra::AMP_MG5_pp_jj       jj;
  const gra::LORENTZSCALAR gg          = PPJJReferenceKinematics(0.118);
  const auto               jj_evaluate = [&](gra::LORENTZSCALAR &lts) { return HardAmp2(jj, lts, 0.118); };
  RequireHelicitySectionContinuity(gg, jj_evaluate, "g g -> u ubar");
  RequireRotatedScreeningRingSection(gg, jj_evaluate, "g g -> u ubar");

  gra::AMP_MG5_yy_zjj yy_zjj;
  RequireHelicitySectionContinuity(
      YYZjjReferenceKinematics(2), [&](gra::LORENTZSCALAR &lts) { return PhotonZjjAmp2(yy_zjj, lts, false); },
      "gamma gamma -> Z u ubar");
}

// Test every registered hard family through the direct covariant projection
TEST_CASE("Registered hard families are Lorentz stable", "[HardPomeronPDF][MG5][covariance]") {
  std::map<std::string, gra::LORENTZSCALAR> references;
  references.emplace("MG5_PP_Z", PPZReferenceKinematics(0.118));
  references.emplace("MG5_PP_ZJ", PPZjReferenceKinematics(0.118));
  references.emplace("MG5_PP_JJ", PPJJReferenceKinematics(0.118));
  references.emplace("MG5_PP_W", PPWReferenceKinematics(0.118));

  std::vector<std::string> families;
  for (const auto &info : gra::PartonMG5ProcessInfos()) { families.push_back(info.process_family); }
  REQUIRE(families.size() == references.size());
  for (const auto &family : families) {
    const auto reference = references.find(family);
    REQUIRE(reference != references.end());
    std::unique_ptr<gra::PartonMG5Process> amplitude = gra::CreatePartonMG5Process(family);
    REQUIRE(amplitude != nullptr);

    gra::LORENTZSCALAR dd = reference->second;
    dd.q1.SetPxPy(0.17, -0.11);
    dd.q2.SetPxPy(-0.17, 0.11);
    RequireHardScalarCovariance(*amplitude, dd, 0.118, family + " DD");

    for (const bool first_diffractive : {true, false}) {
      gra::LORENTZSCALAR sd = reference->second;
      SetHardSDPair(sd, first_diffractive);
      RequireHardScalarCovariance(*amplitude, sd, 0.118, family + (first_diffractive ? " SD leg 1" : " SD leg 2"));
    }
  }
}

// Test the complete generated two-body matrix in the collider helicity section
TEST_CASE("Generated neutral-current entries obey collider Rz covariance",
          "[HardPomeronPDF][MG5][helicity][covariance]") {
  gra::LORENTZSCALAR reference = PPZReferenceKinematics(0.118);

  std::vector<gra::mg5::Subprocess<MG5_PP_Z::ProcessBase>> subprocesses;
  MG5_PP_Z::BuildSubprocesses(subprocesses, gra::aux::ResolveProjectPath("MG5cards/Parton/MG5_PP_Z/param_card.dat"));
  auto *generated = FindChannelProcess(subprocesses, std::array<int, 2>{2, -2}, std::vector<int>{-13, 13});
  REQUIRE(generated != nullptr);
  std::vector<std::array<double, 4>> storage;
  std::vector<double *>              momenta;
  BuildGeneratedTestMomenta(reference, storage, momenta);
  generated->setMomenta(momenta);
  generated->setInitial(reference.id1, reference.id2);
  generated->setAlphaS(reference.alphaQCD);
  const auto labels = generated->helicityAmplitudes();

  gra::AMP_MG5_pp_z amplitude;
  const double      reference_norm       = HardAmp2(amplitude, reference, 0.118);
  const auto        reference_amplitudes = reference.hamp;
  REQUIRE(reference_norm > 0.0);
  REQUIRE(reference_amplitudes.size() == labels.size());
  const double component_scale = std::sqrt(reference_norm / static_cast<double>(labels.size()));

  for (const double angle : {-1.13, 0.47, 2.21}) {
    gra::LORENTZSCALAR rotated      = RotateGeneratedEventAroundZ(reference, angle);
    const double       rotated_norm = HardAmp2(amplitude, rotated, 0.118);
    REQUIRE(rotated_norm == Approx(reference_norm).epsilon(1.0e-11));
    REQUIRE(rotated.hamp.size() == labels.size());

    for (std::size_t index = 0; index < labels.size(); ++index) {
      REQUIRE(labels[index].outgoing.size() == 2);
      const int harmonic = gra::spin::ColliderSpinHalfHardHelicityHarmonic(
          labels[index].incoming[0], labels[index].incoming[1], labels[index].outgoing[0], labels[index].outgoing[1]);
      const std::complex<double> expected =
          std::polar(1.0, static_cast<double>(harmonic) * angle) * reference_amplitudes[index];
      const double tolerance = 2.0e-10 * std::max(component_scale, std::abs(expected));
      CAPTURE(angle, index, labels[index].incoming[0], labels[index].incoming[1], labels[index].outgoing[0],
              labels[index].outgoing[1], harmonic, rotated.hamp[index], expected);
      REQUIRE(std::abs(rotated.hamp[index] - expected) < tolerance);
    }
  }
}

TEST_CASE("Gamma-gamma Z+2-parton explicit flavours and color flow are physical", "[HardPomeronPDF][color]") {
  gra::MODELPARAM = "TUNE0";

  gra::MRandom random;
  random.SetSeed(7);
  gra::AMP_MG5_yy_zjj amp;
  for (const int flavour : {1, 2, 3, 4, 5}) {
    gra::LORENTZSCALAR lts = YYZjjReferenceKinematics(flavour);
    LoadReferencePDG(lts);
    REQUIRE(PhotonZjjAmp2(amp, lts, false) > 0.0);
    REQUIRE(amp.SubprocessCount() == 15);
    lts.decaytree[0].p.color_flow         = {610, 611};
    lts.decaytree[0].legs[0].p.color_flow = {612, 613};
    lts.decaytree[0].legs[1].p.color_flow = {614, 615};
    amp.SampleColorFlow(lts, random);
    const auto leaves = gra::mg5::StableDecayLeaves(lts.decaytree);
    REQUIRE(leaves.size() == 4);
    REQUIRE(leaves[2]->p.pdg == flavour);
    REQUIRE(leaves[3]->p.pdg == -flavour);
    REQUIRE(leaves[2]->p.color_flow.flow1 == 501);
    REQUIRE(leaves[2]->p.color_flow.flow2 == 0);
    REQUIRE(leaves[3]->p.color_flow.flow1 == 0);
    REQUIRE(leaves[3]->p.color_flow.flow2 == 501);
    REQUIRE(lts.decaytree[0].p.color_flow.empty());
    REQUIRE(lts.decaytree[0].legs[0].p.color_flow.empty());
    REQUIRE(lts.decaytree[0].legs[1].p.color_flow.empty());
  }
}

// Check unsupported colored intermediate decays fail before assigning flow
TEST_CASE("Generated color flow rejects colored intermediate decays", "[HardPomeronPDF][color]") {
  gra::MDecayBranch parent;
  parent.p.pdg   = 6;
  parent.p.color = 3;
  parent.legs    = {StableBranch(5, gra::M4Vec()), StableBranch(24, gra::M4Vec())};
  gra::LORENTZSCALAR lts;
  lts.decaytree                                = {parent};
  const std::vector<gra::MColorFlow> candidate = {{501, 0}, {0, 0}};

  CHECK_FALSE(gra::AssignHardColorFlowCandidate(lts, candidate));
  CHECK_FALSE(gra::AssignDurhamColorFlowCandidate(lts, candidate));
}

// Test that incompatible generated MG5 hard-process wrappers are skipped
TEST_CASE(
    "Diffractive Z MG5 hard subprocess channel selection rejects "
    "incompatible modes",
    "[HardPomeronPDF]") {
  {
    gra::LORENTZSCALAR lts = PPZjReferenceKinematics(0.118);
    std::swap(lts.id1, lts.id2);
    std::swap(lts.q1, lts.q2);
    gra::AMP_MG5_pp_zj amp;
    const auto         evaluation = EvaluateHard(amp, lts, 0.118);
    REQUIRE(evaluation.amp2 > 0.0);
    REQUIRE(evaluation.contributing_subprocesses == 1);
  }

  {
    gra::LORENTZSCALAR lts = PPZReferenceKinematics(0.118);
    lts.id1                = 21;
    lts.id2                = 21;
    gra::AMP_MG5_pp_z amp;
    const auto        evaluation = EvaluateHard(amp, lts, 0.118);
    REQUIRE(evaluation.amp2 == Approx(0.0));
    REQUIRE(evaluation.contributing_subprocesses == 0);
    REQUIRE(evaluation.color_flows.empty());
  }

  {
    gra::LORENTZSCALAR lts = PPZjReferenceKinematics(0.118);
    lts.id1                = 21;
    lts.id2                = 21;
    gra::AMP_MG5_pp_zj amp;
    const auto         evaluation = EvaluateHard(amp, lts, 0.118);
    REQUIRE(evaluation.amp2 == Approx(0.0));
    REQUIRE(evaluation.contributing_subprocesses == 0);
    REQUIRE(evaluation.color_flows.empty());
  }

  {
    gra::LORENTZSCALAR lts = PPJJReferenceKinematics(0.118);
    lts.decaytree[1].p.pdg = 5;
    gra::AMP_MG5_pp_jj amp;
    const auto         evaluation = EvaluateHard(amp, lts, 0.118);
    REQUIRE(evaluation.amp2 == Approx(0.0));
    REQUIRE(evaluation.contributing_subprocesses == 0);
    REQUIRE(evaluation.color_flows.empty());
  }

  {
    gra::LORENTZSCALAR lts = PPJJReferenceKinematics(0.118);
    lts.decaytree[0].p.pdg = gra::PDG::PDG_hard_jet;
    lts.decaytree[1].p.pdg = -gra::PDG::PDG_hard_jet;
    gra::AMP_MG5_pp_jj amp;
    const auto         evaluation = EvaluateHard(amp, lts, 0.118);
    REQUIRE_FALSE(evaluation.Valid());
    REQUIRE(evaluation.contributing_subprocesses == 0);
  }
}

// Test the generated stable final state and outgoing hard state
TEST_CASE("Diffractive Z MG5 channel selection validates generated final states", "[HardPomeronPDF][MG5]") {
  SECTION("ordered stable final state") {
    gra::LORENTZSCALAR lts = PPZReferenceKinematics(0.118);
    gra::AMP_MG5_pp_z  amp;
    auto               evaluation = EvaluateHard(amp, lts, 0.118);
    REQUIRE(evaluation.amp2 > 0.0);
    REQUIRE_FALSE(lts.hamp.empty());
    REQUIRE(evaluation.contributing_subprocesses == 1);

    std::swap(lts.decaytree[0].p.pdg, lts.decaytree[1].p.pdg);
    evaluation = EvaluateHard(amp, lts, 0.118);
    REQUIRE_FALSE(evaluation.Valid());
    REQUIRE(lts.hamp.empty());
    REQUIRE(evaluation.contributing_subprocesses == 0);
  }

  SECTION("explicit outgoing flavour") {
    gra::LORENTZSCALAR lts = PPZjReferenceKinematics(0.118);
    gra::AMP_MG5_pp_zj amp;
    const auto         evaluation = EvaluateHard(amp, lts, 0.118);
    REQUIRE(evaluation.amp2 > 0.0);
    REQUIRE(evaluation.contributing_subprocesses == 1);
  }

  SECTION("unresolved outgoing selector") {
    for (const int selector : {gra::PDG::PDG_hard_jet, -gra::PDG::PDG_hard_jet}) {
      gra::LORENTZSCALAR lts = PPZjReferenceKinematics(0.118);
      lts.decaytree[1].p.pdg = selector;
      gra::AMP_MG5_pp_zj amp;
      REQUIRE_FALSE(EvaluateHard(amp, lts, 0.118).Valid());
    }
  }
}

// Test that MG5 decay-chain amplitudes remain independent of phase-space
// proposals
TEST_CASE("Diffractive Z MG5 wrappers return the raw decay-chain amplitude", "[HardPomeronPDF]") {
  gra::AMP_MG5_yy_zjj amp;

  gra::LORENTZSCALAR raw_lts  = YYZjjReferenceKinematics();
  const double       raw_amp2 = PhotonZjjAmp2(amp, raw_lts, false);

  gra::LORENTZSCALAR compensated_lts               = YYZjjReferenceKinematics();
  compensated_lts.decaytree[0].p.mass              = 91.188;
  compensated_lts.decaytree[0].p.width             = 2.441404;
  compensated_lts.decaytree[0].mass_proposal = gra::MassProposal::BreitWigner;

  const double compensated_amp2 = PhotonZjjAmp2(amp, compensated_lts, false);

  REQUIRE(compensated_amp2 == Approx(raw_amp2).epsilon(1e-12));
  RequireUnaveragedHelicityHamp(compensated_lts, compensated_amp2);

  gra::LORENTZSCALAR isolated_lts      = YYZjjReferenceKinematics();
  isolated_lts.process.root_decay_mode = gra::RootDecayMode::Isolated;
  BindHardModelCache(isolated_lts);
  const auto isolated_result = amp.Evaluate(isolated_lts, 0.0, false);
  REQUIRE(isolated_result.status == gra::mg5helas::EvaluationStatus::AmplitudeFailure);
  REQUIRE(gra::math::IsZero(isolated_result.amp2));
}

// Test process-level gamma-gamma flux wrappers against standard formulae
TEST_CASE("Gamma-gamma Zjj process wrappers match flux formula references", "[HardPomeronPDF]") {
  gra::LORENTZSCALAR  raw_lts = YYZjjFluxReferenceKinematics();
  gra::AMP_MG5_yy_zjj amp;
  const double        raw_epa_amp2       = PhotonZjjAmp2(amp, raw_lts, true);
  gra::LORENTZSCALAR  raw_collinear_lts  = YYZjjFluxReferenceKinematics();
  const double        raw_collinear_amp2 = PhotonZjjAmp2(amp, raw_collinear_lts, false);
  REQUIRE(raw_epa_amp2 > 0.0);
  REQUIRE(raw_collinear_amp2 > 0.0);

  {
    gra::LORENTZSCALAR lts = YYZjjFluxReferenceKinematics();
    lts.alphaQCD           = 0.0;
    LoadReferencePDG(lts);
    BindHardModelCache(lts);
    gra::MGeneratedPhotonProc proc("yy", "Zjj", {"generated", "MG5", "yy", "", 1}, "MG5_YY_ZJJ");
    proc.BindModelTune(gra::MModelTune::Load(modelfile));
    proc.InitializeWorkerAmplitude(lts);
    const double got       = proc.Amp2(lts);
    const auto  &structure = lts.model_cache->Tune().Structure();
    const double expected  = raw_epa_amp2 * gra::flux::CohFlux(lts.x1, lts.t1, lts.qt1, structure) *
                            gra::flux::CohFlux(lts.x2, lts.t2, lts.qt2, structure) *
                            gra::flux::ExactktEPAPhaseSpaceFactor(lts);
    REQUIRE(expected > 0.0);
    REQUIRE(got == Approx(expected).epsilon(1e-12));
    RequireUnaveragedHelicityHamp(lts, got);
    REQUIRE(lts.id1 == gra::PDG::PDG_gamma);
    REQUIRE(lts.id2 == gra::PDG::PDG_gamma);
    REQUIRE(lts.exact_forward_photon_kinematics);
    REQUIRE(lts.muF == Approx(std::sqrt(lts.s_hat) / 2.0));
    REQUIRE(lts.muR == Approx(lts.muF));
    REQUIRE(lts.scalup == Approx(lts.muF));
    REQUIRE(gra::math::IsZero(lts.alphaQCD));
  }

  {
    gra::LORENTZSCALAR lts = YYZjjFluxReferenceKinematics();
    lts.alphaQCD           = 0.0;
    LoadReferencePDG(lts);
    BindHardModelCache(lts);
    gra::MGeneratedPhotonProc proc("yy_DZ", "Zjj", {"generated", "MG5", "yy", "", 1}, "MG5_YY_ZJJ");
    proc.BindModelTune(gra::MModelTune::Load(modelfile));
    proc.InitializeWorkerAmplitude(lts);
    const double got      = proc.Amp2(lts);
    const double expected = raw_collinear_amp2 * gra::flux::DZFlux(lts.x1) * gra::flux::DZFlux(lts.x2) *
                            gra::flux::CollinearPhotonPhaseSpaceFactor(lts);
    REQUIRE(expected > 0.0);
    REQUIRE(got == Approx(expected).epsilon(1e-12));
    RequireUnaveragedHelicityHamp(lts, got);
    REQUIRE(lts.id1 == gra::PDG::PDG_gamma);
    REQUIRE(lts.id2 == gra::PDG::PDG_gamma);
    REQUIRE_FALSE(lts.exact_forward_photon_kinematics);
    REQUIRE(lts.muF == Approx(std::sqrt(lts.s_hat) / 2.0));
    REQUIRE(lts.muR == Approx(lts.muF));
    REQUIRE(lts.scalup == Approx(lts.muF));
    REQUIRE(gra::math::IsZero(lts.alphaQCD));
  }

  {
    gra::LORENTZSCALAR lts = YYZjjFluxReferenceKinematics();
    lts.alphaQCD           = 0.0;
    LoadReferencePDG(lts);
    BindHardModelCache(lts);
    gra::MGeneratedPhotonProc proc("yy_LUX", "Zjj", {"generated", "MG5", "yy", "", 1}, "MG5_YY_ZJJ");
    proc.BindModelTune(gra::MModelTune::Load(modelfile));
    gra::MLHAPDFStore pdf_store;
    lts.GlobalPdfPtr = pdf_store.GetPDF(lts.LHAPDFSET, 0);
    proc.InitializeWorkerAmplitude(lts);
    const double got      = proc.Amp2(lts);
    const auto   pdf      = pdf_store.GetPDF(lts.LHAPDFSET, 0);
    const double Q2       = lts.s_hat / 4.0;
    const double f1       = pdf->xfxQ2(gra::PDG::PDG_gamma, lts.x1, Q2) / lts.x1;
    const double f2       = pdf->xfxQ2(gra::PDG::PDG_gamma, lts.x2, Q2) / lts.x2;
    const double expected = raw_collinear_amp2 * f1 * f2 * gra::flux::CollinearPhotonPhaseSpaceFactor(lts);
    REQUIRE(expected > 0.0);
    REQUIRE(got == Approx(expected).epsilon(1e-12));
    RequireUnaveragedHelicityHamp(lts, got);
    REQUIRE(lts.id1 == gra::PDG::PDG_gamma);
    REQUIRE(lts.id2 == gra::PDG::PDG_gamma);
    REQUIRE(lts.pdf_xf1 == Approx(lts.x1 * f1).epsilon(1e-13));
    REQUIRE(lts.pdf_xf2 == Approx(lts.x2 * f2).epsilon(1e-13));
    REQUIRE_FALSE(lts.exact_forward_photon_kinematics);
    REQUIRE(lts.muF == Approx(std::sqrt(Q2)));
    REQUIRE(lts.muR == Approx(lts.muF));
    REQUIRE(lts.scalup == Approx(lts.muF));
    REQUIRE(gra::math::IsZero(lts.alphaQCD));
  }
}

// Test leading-order alpha_s powers expected for generated hard channels
TEST_CASE("Diffractive MG5 wrappers obey LO alpha_s powers", "[HardPomeronPDF]") {
  {
    gra::AMP_MG5_pp_z  amp;
    gra::LORENTZSCALAR lts_half = PPZReferenceKinematics(0.0);
    gra::LORENTZSCALAR lts_nom  = PPZReferenceKinematics(0.0);
    gra::LORENTZSCALAR lts_dbl  = PPZReferenceKinematics(0.0);
    const double       half     = HardAmp2(amp, lts_half, 0.059);
    const double       nominal  = HardAmp2(amp, lts_nom, 0.118);
    const double       doubled  = HardAmp2(amp, lts_dbl, 0.236);

    REQUIRE(half == Approx(nominal).epsilon(1e-12));
    REQUIRE(doubled == Approx(nominal).epsilon(1e-12));
  }

  {
    gra::AMP_MG5_pp_zj amp;
    gra::LORENTZSCALAR lts_half = PPZjReferenceKinematics(0.0);
    gra::LORENTZSCALAR lts_nom  = PPZjReferenceKinematics(0.0);
    gra::LORENTZSCALAR lts_dbl  = PPZjReferenceKinematics(0.0);
    const double       half     = HardAmp2(amp, lts_half, 0.059);
    const double       nominal  = HardAmp2(amp, lts_nom, 0.118);
    const double       doubled  = HardAmp2(amp, lts_dbl, 0.236);

    REQUIRE(half == Approx(0.5 * nominal).epsilon(1e-12));
    REQUIRE(doubled == Approx(2.0 * nominal).epsilon(1e-12));
  }

  {
    gra::AMP_MG5_pp_jj amp;
    gra::LORENTZSCALAR lts_half = PPJJReferenceKinematics(0.0);
    gra::LORENTZSCALAR lts_nom  = PPJJReferenceKinematics(0.0);
    gra::LORENTZSCALAR lts_dbl  = PPJJReferenceKinematics(0.0);
    const double       half     = HardAmp2(amp, lts_half, 0.059);
    const double       nominal  = HardAmp2(amp, lts_nom, 0.118);
    const double       doubled  = HardAmp2(amp, lts_dbl, 0.236);

    REQUIRE(half == Approx(0.25 * nominal).epsilon(1e-12));
    REQUIRE(doubled == Approx(4.0 * nominal).epsilon(1e-12));
  }

  {
    gra::MGeneratedPartonProc proc("IPp", "Zj", {"generated", "MG5", "IPp", "", 1}, "MG5_PP_ZJ");
    gra::LORENTZSCALAR        lts_half = PPZjReferenceKinematics(0.059);
    gra::LORENTZSCALAR        lts_nom  = PPZjReferenceKinematics(0.118);
    proc.InitializeWorkerAmplitude(lts_nom);
    const double half    = proc.Amp2(lts_half);
    const double nominal = proc.Amp2(lts_nom);
    REQUIRE(half == Approx(0.5 * nominal).epsilon(1e-12));
  }
}

// Test that hard-diffraction event records expose showerable color flow
TEST_CASE("Hard diffraction event record keeps remnant color flow", "[HardPomeronPDF]") {
  {
    HardDiffractionTestProbe proc;
    proc.ProcPtr.ISTATE  = "IPp";
    proc.ProcPtr.CHANNEL = "Z";
    proc.SetModelTune(gra::MModelTune::Load(modelfile));
    proc.state.lts.PDG.ReadParticleData();

    ConfigureHardColorProbe(proc, false, PPZReferenceKinematics(0.118).decaytree);
    proc.state.lts.id1     = 2;
    proc.state.lts.id2     = -2;
    proc.state.lts.pdf_xf1 = 1.2;
    proc.state.lts.pdf_xf2 = 2.3;
    proc.state.lts.muF     = 91.1876;
    proc.state.lts.excite1 = true;
    proc.state.lts.excite2 = true;

    HepMC3::GenEvent evt(HepMC3::Units::GEV, HepMC3::Units::MM);
    REQUIRE(proc.FinalizeProbe());
    REQUIRE(proc.EventRecord(evt));
    REQUIRE_FALSE(proc.state.lts.excite1);
    REQUIRE_FALSE(proc.state.lts.excite2);
    auto cross_section = std::make_shared<HepMC3::GenCrossSection>();
    cross_section->set_cross_section(12.5, 0.25);
    cross_section->set_accepted_events(17);
    cross_section->set_attempted_events(431);
    evt.set_cross_section(cross_section);

    int    nstar_with_leading_proton = 0;
    double leading_t                 = 0.0;
    double leading_pt                = 0.0;
    for (const auto &particle : evt.particles()) {
      if (std::abs(particle->pid()) != gra::PDG::PDG_NSTAR || particle->momentum().pz() < 0.0 ||
          particle->end_vertex() == nullptr) {
        continue;
      }
      for (const auto &child : particle->end_vertex()->particles_out()) {
        if (child->pid() == 2212) {
          ++nstar_with_leading_proton;
          const HepMC3::FourVector transfer(
              proc.state.lts.pbeam1.Px() - child->momentum().px(), proc.state.lts.pbeam1.Py() - child->momentum().py(),
              proc.state.lts.pbeam1.Pz() - child->momentum().pz(), proc.state.lts.pbeam1.E() - child->momentum().e());
          leading_t  = HepMCMass2(transfer);
          leading_pt = std::sqrt(child->momentum().px() * child->momentum().px() +
                                 child->momentum().py() * child->momentum().py());
        }
      }
    }
    REQUIRE(nstar_with_leading_proton == 1);
    REQUIRE(leading_t == Approx(proc.state.lts.diff_t1).epsilon(1e-8));
    REQUIRE(leading_pt > 0.0);
    const std::string diff_comment = gra::BuildLHEDiffractionComment(evt);
    REQUIRE(diff_comment.find("graniitti_diff") != std::string::npos);
    REQUIRE(diff_comment.find("version=1") != std::string::npos);
    REQUIRE(diff_comment.find("accepted_events=17") != std::string::npos);
    REQUIRE(diff_comment.find("attempted_events=431") != std::string::npos);
    REQUIRE(diff_comment.find("side_mask=1") != std::string::npos);

    const auto pdf_info = evt.pdf_info();
    REQUIRE(pdf_info != nullptr);
    REQUIRE(pdf_info->is_valid());
    REQUIRE(pdf_info->parton_id[0] == 2);
    REQUIRE(pdf_info->parton_id[1] == -2);
    REQUIRE(pdf_info->x[0] == Approx(proc.state.lts.diff_xhard1));
    REQUIRE(pdf_info->x[1] == Approx(proc.state.lts.diff_xhard2));
    REQUIRE(pdf_info->scale == Approx(proc.state.lts.muF));

    int                colored_with_flow = 0;
    std::map<int, int> color_balance;
    for (const auto &particle : evt.particles()) {
      if (!IsColoredParton(particle->pid())) { continue; }
      const int flow1 = ReadIntAttribute(particle, "flow1");
      const int flow2 = ReadIntAttribute(particle, "flow2");
      if (flow1 != 0 || flow2 != 0) { ++colored_with_flow; }
      if (flow1 != 0) { ++color_balance[flow1]; }
      if (flow2 != 0) { --color_balance[flow2]; }
    }

    REQUIRE(colored_with_flow >= 2);
    for (const auto &[tag, balance] : color_balance) {
      REQUIRE(tag >= 501);
      REQUIRE(balance == 0);
    }

    const std::vector<gra::MLHERow> rows = gra::BuildLHERows(evt);
    REQUIRE(rows.size() >= 4);
    REQUIRE(rows[0].status == -1);
    REQUIRE(rows[1].status == -1);
    REQUIRE(rows[0].pid == 2212);
    REQUIRE(rows[1].pid == 2212);
    REQUIRE(rows[0].momentum[0] == Approx(proc.state.lts.pbeam1.Px()).epsilon(1e-10));
    REQUIRE(rows[0].momentum[1] == Approx(proc.state.lts.pbeam1.Py()).epsilon(1e-10));
    REQUIRE(rows[0].momentum[4] == Approx(gra::PDG::mp).epsilon(1e-7));
    for (const auto &row : rows) {
      REQUIRE(row.pid != gra::PDG::PDG_system);
      REQUIRE(row.pid != gra::PDG::PDG_propagator);
    }

    const std::filesystem::path lhe_path = "tmp/graniitti_hardpomeron_lhe_test.lhe";
    evt.weights()                        = {1.0};
    {
      gra::MLHEWriter writer(lhe_path.string());
      writer.WriteEvent(evt);
    }
    std::ifstream lhe_input(lhe_path);
    REQUIRE(lhe_input.good());
    const std::string lhe_text((std::istreambuf_iterator<char>(lhe_input)), std::istreambuf_iterator<char>());
    REQUIRE(lhe_text.find("<LesHouchesEvents") != std::string::npos);
    REQUIRE(lhe_text.find("<event>") != std::string::npos);
    REQUIRE(lhe_text.find(" 2212 -1 ") != std::string::npos);
    REQUIRE(lhe_text.find("graniitti_diff") != std::string::npos);
    REQUIRE(lhe_text.find("side_mask=1") != std::string::npos);
    REQUIRE(lhe_text.find("side1_xi=") != std::string::npos);
    REQUIRE(lhe_text.find("side1_lead_pdg=2212") != std::string::npos);
    REQUIRE(lhe_text.find("side1_rem_pdg=") != std::string::npos);
  }

  {
    HardDiffractionTestProbe proc;
    proc.ProcPtr.ISTATE  = "IPp";
    proc.ProcPtr.CHANNEL = "Zj";
    proc.SetModelTune(gra::MModelTune::Load(modelfile));
    proc.state.lts.PDG.ReadParticleData();

    ConfigureHardColorProbe(proc, false, PPZjReferenceKinematics(0.118).decaytree);
    proc.state.lts.id1     = 21;
    proc.state.lts.id2     = 2;
    proc.state.lts.pdf_xf1 = 1.2;
    proc.state.lts.pdf_xf2 = 2.3;
    proc.state.lts.muF     = 91.1876;

    HepMC3::GenEvent evt(HepMC3::Units::GEV, HepMC3::Units::MM);
    REQUIRE(proc.FinalizeProbe());
    REQUIRE(proc.EventRecord(evt));

    // Check the net remnant colors have the incoming gluon and quark representations
    const auto gluon_color = evt.attribute<HepMC3::IntAttribute>("graniitti_hard1_flow1");
    const auto gluon_anticolor = evt.attribute<HepMC3::IntAttribute>("graniitti_hard1_flow2");
    const auto quark_color = evt.attribute<HepMC3::IntAttribute>("graniitti_hard2_flow1");
    REQUIRE(gluon_color);
    REQUIRE(gluon_anticolor);
    REQUIRE(gluon_color->value() != gluon_anticolor->value());
    REQUIRE(quark_color);
    REQUIRE_FALSE(evt.attribute<HepMC3::IntAttribute>("graniitti_hard2_flow2"));

    const auto pdf_info = evt.pdf_info();
    REQUIRE(pdf_info != nullptr);
    REQUIRE(pdf_info->is_valid());
    REQUIRE(pdf_info->parton_id[0] == 21);
    REQUIRE(pdf_info->parton_id[1] == 2);
    REQUIRE(pdf_info->x[0] == Approx(proc.state.lts.diff_xhard1));
    REQUIRE(pdf_info->x[1] == Approx(proc.state.lts.diff_xhard2));
    REQUIRE(pdf_info->scale == Approx(proc.state.lts.muF));

    int                colored_with_flow = 0;
    std::map<int, int> color_balance;
    for (const auto &particle : evt.particles()) {
      if (!IsColoredParton(particle->pid())) { continue; }
      const int flow1 = ReadIntAttribute(particle, "flow1");
      const int flow2 = ReadIntAttribute(particle, "flow2");
      if (flow1 != 0 || flow2 != 0) { ++colored_with_flow; }
      if (flow1 != 0) { ++color_balance[flow1]; }
      if (flow2 != 0) { --color_balance[flow2]; }
    }

    REQUIRE(colored_with_flow >= 3);
    for (const auto &[tag, balance] : color_balance) {
      REQUIRE(tag >= 501);
      REQUIRE(balance == 0);
    }
  }
}

// Check constituent quantum numbers, mass shells and screening transport for every ordinary flavour
TEST_CASE("Proton remnants conserve flavour and color through screening", "[HardPomeronPDF][color][kinematics]") {
  HardDiffractionTestProbe proc;
  proc.state.lts.PDG.ReadParticleData();
  const double mass         = gra::MHardPomeronPDFParam{}.remnant_mass;
  proc.state.lts.pbeam1     = gra::M4Vec(0, 0, 100, std::hypot(100.0, gra::PDG::mp));
  proc.state.lts.pbeam2     = gra::M4Vec(0, 0, -100, std::hypot(100.0, gra::PDG::mp));
  proc.state.lts.pfinal     = {gra::M4Vec(0, 0, 0, 20), gra::M4Vec(0, 0, 70, std::hypot(70.0, mass)),
                               gra::M4Vec(0, 0, -70, std::hypot(70.0, mass))};
  proc.state.lts.beam1.pdg  = 2212;
  proc.state.lts.beam2.pdg  = -2212;
  proc.state.lts.hard_diff1 = proc.state.lts.hard_diff2 = false;

  // Compute three times baryon number directly from constituent valence content
  const auto baryon3 = [](int pdg) {
    const int id   = std::abs(pdg);
    const int sign = pdg > 0 ? 1 : -1;
    return id == 2101 || id == 2203 ? 2 * sign : (id >= 1000 && id < 6000 ? 3 * sign : (id <= 5 ? sign : 0));
  };
  // Resolve net d, u, s, c and b numbers from the PDG constituent convention
  const auto net_flavour = [](int pdg) {
    std::array<int, 5> result{};
    const int          id   = std::abs(pdg);
    const int          sign = pdg > 0 ? 1 : -1;
    if (id <= 5) {
      result[id - 1] += sign;
    } else if (id >= 1000 && id < 6000) {
      result[id / 1000 - 1] += sign;
      result[(id / 100) % 10 - 1] += sign;
      if ((id / 10) % 10 > 0) { result[(id / 10) % 10 - 1] += sign; }
    } else if (id >= 100 && id < 600) {
      const int heavy      = id / 100;
      const int light      = (id / 10) % 10;
      const int heavy_sign = heavy % 2 == 0 ? sign : -sign;
      result[heavy - 1] += heavy_sign;
      result[light - 1] -= heavy_sign;
    }
    return result;
  };
  for (const int flavour : {21, 1, 2, -1, -2, 3, -3, 4, -4, 5, -5}) {
    INFO(flavour);
    proc.state.lts.screening.active = false;
    proc.state.lts.id1              = flavour;
    proc.state.lts.id2              = flavour == 21 ? 21 : -flavour;
    REQUIRE(proc.FinalizeProbe());
    const auto                        born1 = proc.state.lts.decayforward1;
    const auto                        born2 = proc.state.lts.decayforward2;
    std::map<int, std::array<int, 2>> tags;
    for (const bool first : {true, false}) {
      const auto &branch  = first ? born1 : born2;
      const int   id      = first ? proc.state.lts.id1 : proc.state.lts.id2;
      int         charge3 = proc.state.lts.PDG.FindByPDG(id).chargeX3;
      int         baryon  = baryon3(id);
      auto        net     = net_flavour(id);
      gra::M4Vec  sum;
      REQUIRE(branch.legs.size() == 2);
      for (const auto &leg : branch.legs) {
        charge3 += leg.p.chargeX3;
        baryon += baryon3(leg.p.pdg);
        const auto constituent = net_flavour(leg.p.pdg);
        for (const auto &i : indices(net)) { net[i] += constituent[i]; }
        sum += leg.p4;
        REQUIRE(leg.p4.E() > 0.0);
        REQUIRE(leg.p4.M2() == Approx(leg.p.mass * leg.p.mass).margin(1e-9));
        if (leg.p.color_flow.flow1) { ++tags[leg.p.color_flow.flow1][0]; }
        if (leg.p.color_flow.flow2) { ++tags[leg.p.color_flow.flow2][1]; }
      }
      REQUIRE(charge3 == (first ? 3 : -3));
      REQUIRE(baryon == (first ? 3 : -3));
      REQUIRE(net == net_flavour(first ? 2212 : -2212));
      REQUIRE(gra::math::CheckEMC(branch.p4 - sum));
    }
    for (const auto &[tag, counts] : tags) {
      REQUIRE(tag >= 501);
      REQUIRE(counts[0] == 1);
      REQUIRE(counts[1] == 1);
    }

    proc.state.lts.pfinal_orig = proc.state.lts.pfinal;
    proc.state.lts.pfinal[1].RotateY(0.07);
    proc.state.lts.pfinal[2].RotateY(-0.09);
    proc.state.lts.screening.active = true;
    REQUIRE(proc.FinalizeProbe());
    for (const bool first : {true, false}) {
      const auto &born = first ? born1 : born2;
      const auto &loop = first ? proc.state.lts.decayforward1 : proc.state.lts.decayforward2;
      for (const auto &i : indices(born.legs)) {
        auto expected = born.legs[i].p4;
        gra::kinematics::LorentzBoost(born.p4, born.p4.M(), expected, -1);
        gra::kinematics::LorentzBoost(loop.p4, loop.p4.M(), expected, 1);
        REQUIRE(gra::math::CheckEMC(expected - loop.legs[i].p4));
      }
    }
    proc.state.lts.pfinal           = proc.state.lts.pfinal_orig;
    proc.state.lts.screening.active = false;
    proc.state.lts.id1              = 6;
    REQUIRE_FALSE(proc.FinalizeProbe());
    REQUIRE(proc.state.lts.decayforward1.legs[0].p.pdg == born1.legs[0].p.pdg);
  }
}
