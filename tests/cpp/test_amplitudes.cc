// Hard amplitude tests
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include <span>
#include <type_traits>

#include "Graniitti/Amplitude/MG5/Photon/AMP_MG5_PhotonRegistry.h"
#include "Graniitti/Amplitude/MG5/Photon/AMP_MG5_yy_jj.h"
#include "Graniitti/Amplitude/MG5/Photon/AMP_MG5_yy_ll.h"
#include "Graniitti/Amplitude/MG5/Photon/AMP_MG5_yy_ww.h"
#include "Graniitti/Amplitude/MG5/Photon/AMP_MG5_yy_zjj.h"
#include "Graniitti/Amplitude/MG5/Runtime/Models/sm/HelAmps_sm_lepton_masses.h"
#include "Graniitti/Amplitude/MG5/Runtime/Models/sm/Parameters_sm_lepton_masses.h"
#include "Graniitti/Amplitude/Photon/AMP_yy_yy.h"
#include "Graniitti/Kinematics/MEventKinematics.h"
#include "Graniitti/Math/MFloat.h"
#include "Graniitti/Math/MIntegration.h"
#include "Graniitti/Nuclear/MConfig.h"
#include "Graniitti/Nuclear/MNucleus.h"
#include "Graniitti/Nuclear/MUPC.h"
#include "Graniitti/Nuclear/MUPCScreen.h"
#include "Graniitti/PDF/MLHAPDF.h"
#include "Graniitti/Particle/MResonance.h"
#include "Graniitti/Photon/MFlux.h"
#include "Graniitti/Photon/MGamma.h"
#include "support/models_test_support.hh"
#include "support/nuclear_test_support.hh"

namespace {

// Complete an exact elastic two-leg state from requested transfer momenta
void CompleteElasticForwardState(gra::LORENTZSCALAR &lts, double beam_pz, double beam_mass = gra::PDG::mp) {
  lts.model_cache = std::make_shared<gra::MModelCache>(gra::MModelTune::Load(modelfile));
  lts.LHAPDFSET   = "NNPDF31_lo_as_0118";
  // Supply the physical proton charge when a compact test did not set beams
  if (lts.beam1.pdg == 0) {
    lts.beam1.pdg      = gra::PDG::PDG_p;
    lts.beam1.chargeX3 = 3;
  }
  if (lts.beam2.pdg == 0) { lts.beam2 = lts.beam1; }
  if (std::abs(lts.beam1.pdg) == gra::PDG::PDG_p) {
    lts.beam1.mass   = beam_mass;
    lts.beam1.spinX2 = 1;
  }
  if (std::abs(lts.beam2.pdg) == gra::PDG::PDG_p) {
    lts.beam2.mass   = beam_mass;
    lts.beam2.spinX2 = 1;
  }
  const double beam_energy = std::sqrt(beam_pz * beam_pz + beam_mass * beam_mass);
  lts.pbeam1               = gra::M4Vec(0.0, 0.0, beam_pz, beam_energy);
  lts.pbeam2               = gra::M4Vec(0.0, 0.0, -beam_pz, beam_energy);

  const auto elastic_transfer = [beam_mass](const gra::M4Vec &incoming, const gra::M4Vec &transfer) {
    const double     px = incoming.Px() - transfer.Px();
    const double     py = incoming.Py() - transfer.Py();
    const double     pz = incoming.Pz() - transfer.Pz();
    const gra::M4Vec outgoing(px, py, pz, std::sqrt(px * px + py * py + pz * pz + beam_mass * beam_mass));
    return incoming - outgoing;
  };
  lts.q1            = elastic_transfer(lts.pbeam1, lts.q1);
  lts.q2            = elastic_transfer(lts.pbeam2, lts.q2);
  lts.pfinal[0]     = lts.q1 + lts.q2;
  lts.pfinal[1]     = lts.pbeam1 - lts.q1;
  lts.pfinal[2]     = lts.pbeam2 - lts.q2;
  lts.forward_mass2 = {beam_mass * beam_mass, beam_mass * beam_mass};
  lts.x1            = gra::kinematics::LongitudinalMomentumLoss(lts.pbeam1, lts.pfinal[1], true);
  lts.x2            = gra::kinematics::LongitudinalMomentumLoss(lts.pbeam2, lts.pfinal[2], false);
  lts.xi1           = lts.x1;
  lts.xi2           = lts.x2;
  lts.has_xi1       = true;
  lts.has_xi2       = true;
  lts.t1            = lts.q1.M2();
  lts.t2            = lts.q2.M2();
  lts.qt1           = lts.q1.Pt();
  lts.qt2           = lts.q2.Pt();
}

// Rebuild both on-shell forward legs at one fixed-hard EPA loop node
bool ShiftEPAHardSources(gra::LORENTZSCALAR &lts, const double kx, const double ky) {
  gra::M4Vec forward1;
  gra::M4Vec forward2;
  gra::M4Vec q1;
  gra::M4Vec q2;
  double     t1 = 0.0;
  double     t2 = 0.0;
  if (!gra::kinematics::BuildEPATransfer(lts.pbeam1, lts.xi1, lts.q1.Px() + kx, lts.q1.Py() + ky, pow2(lts.beam1.mass),
                                         lts.forward_mass2[0], true, forward1, q1, t1) ||
      !gra::kinematics::BuildEPATransfer(lts.pbeam2, lts.xi2, lts.q2.Px() - kx, lts.q2.Py() - ky, pow2(lts.beam2.mass),
                                         lts.forward_mass2[1], false, forward2, q2, t2)) {
    return false;
  }
  lts.pfinal[1] = forward1;
  lts.pfinal[2] = forward2;
  lts.q1        = q1;
  lts.q2        = q2;
  lts.t1        = t1;
  lts.t2        = t2;
  lts.qt1       = q1.Pt();
  lts.qt2       = q2.Pt();
  return true;
}

// Construct a compact configuration-screened oxygen UPC model
std::shared_ptr<const gra::nuclear::MUPC> OxygenUPC(const gra::nuclear::CoherenceType emission) {
  gra::nuclear::NucleusParam nucleus_param;
  nucleus_param.pdg    = gra::nuclear::EncodeNuclearPDG(16, 8);
  nucleus_param.mass   = 14.899;
  nucleus_param.charge = {2.608, 0.513, 11.842, 64, 2048, 3.0, 4097, 1.0e-7};
  nucleus_param.matter = nucleus_param.charge;
  auto nucleus         = std::make_shared<const gra::nuclear::MNucleus>(nucleus_param);

  gra::nuclear::ConfigParam config;
  config.count             = 4;
  config.max_trials        = 10000;
  config.neutron_nodes     = 4096;
  config.density_rel_tol   = 1.0e-6;
  config.negative_norm_tol = 1.0e-7;
  config.norm_tol          = 2.0e-5;
  gra::nuclear::UPCParam param;
  param.loop.radial_integrator  = "GL";
  param.loop.azimuth_integrator = "Trap";
  param.loop.radial_map         = gra::math::RadialMap::Square;
  param.emission                = {emission, emission};
  param.structure               = gra::nuclear::StructureType::Nucleon;
  param.survival                = gra::nuclear::SurvivalType::MCGGCF;
  param.config                  = config;
  gra::test::Glauber(param.glauber, 0.0);
  param.glauber.b_max              = 20.0;
  param.glauber.q_max              = 3.0;
  param.glauber.b_nodes            = 32;
  param.glauber.q_nodes            = 32;
  param.convolution.b_max          = 20.0;
  param.convolution.smooth_b_nodes = 32;
  param.convolution.sample_b_nodes = 32;
  gra::test::Convolution(param.convolution);
  for (auto &photo : param.photo) {
    gra::test::Photo(photo);
    photo.b_nodes = 32;
    photo.z_nodes = 32;
  }
  param.loop.r_min              = 0.01;
  param.loop.r_max              = 0.12;
  param.loop.radial_intervals   = 8;
  param.loop.azimuth_nodes      = 4;
  param.convolution.b_phi_nodes = 5;
  const auto model              = std::make_shared<const gra::nuclear::MUPC>(
      std::array<gra::nuclear::BeamType, 2>{gra::nuclear::BeamType::Nucleus, gra::nuclear::BeamType::Nucleus},
      std::array<std::shared_ptr<const gra::nuclear::MNucleus>, 2>{nucleus, nucleus}, param);
  gra::MRandom random;
  random.SetSeed(9137);
  return model->Sample(random);
}

// Build one stable generated photon final-state branch
gra::MDecayBranch PhotonLeaf(int pdg, const gra::M4Vec &p4) {
  gra::MDecayBranch branch;
  branch.p.pdg = pdg;
  branch.p4    = p4;
  return branch;
}

// Build two opposite momenta with a requested common mass
std::array<gra::M4Vec, 2> PhotonPair(double mass, double particle_mass) {
  const double energy   = 0.5 * mass;
  const double momentum = std::sqrt(energy * energy - particle_mass * particle_mass);
  gra::M4Vec   first(momentum, 0.0, 0.0, energy);
  first.RotateY(0.73);
  first.RotateZ(-0.41);
  return {first, gra::M4Vec(-first.Px(), -first.Py(), -first.Pz(), energy)};
}

// Build one inclusive scalar diphoton point with proton or nuclear beams
gra::LORENTZSCALAR ScalarPhotonPoint(const int pdg, const bool nuclear) {
  gra::LORENTZSCALAR lts;
  lts.PDG                             = LoadedPDGTable();
  lts.process.root_decay_mode         = gra::RootDecayMode::None;
  lts.id1                             = gra::PDG::PDG_gamma;
  lts.id2                             = gra::PDG::PDG_gamma;

  double beam_mass = gra::PDG::mp;
  double beam_pz   = 6500.0;
  if (nuclear) {
    beam_mass          = 14.899;
    beam_pz            = 2000.0;
    lts.beam1.pdg      = gra::nuclear::EncodeNuclearPDG(16, 8);
    lts.beam1.mass     = beam_mass;
    lts.beam1.chargeX3 = 24;
    lts.beam1.spinX2   = 0;
    lts.upc_model      = OxygenUPC(gra::nuclear::CoherenceType::Coherent);
  } else {
    lts.beam1.pdg      = gra::PDG::PDG_p;
    lts.beam1.mass     = beam_mass;
    lts.beam1.chargeX3 = 3;
    lts.beam1.spinX2   = 1;
  }
  lts.beam2 = lts.beam1;

  const gra::MParticle scalar        = lts.PDG.FindByPDG(pdg);
  const double         photon_energy = 0.6 * scalar.mass;
  lts.q1                             = gra::M4Vec(0.21, -0.14, photon_energy, photon_energy);
  lts.q2                             = gra::M4Vec(-0.21, 0.14, -photon_energy, photon_energy);
  CompleteElasticForwardState(lts, beam_pz, beam_mass);
  lts.s      = (lts.pbeam1 + lts.pbeam2).M2();
  lts.s_hat  = lts.pfinal[0].M2();
  lts.m2     = lts.s_hat;
  lts.sqrt_s = std::sqrt(lts.s);

  gra::MDecayBranch root;
  root.p        = scalar;
  root.p4       = lts.pfinal[0];
  lts.decaytree = {root};
  return lts;
}

// Build the CP-even Jz=0 scalar hard tensor before EPA source contraction
std::vector<gra::mg5helas::HelicityComponent> ScalarHardTensor(const std::complex<double> same_helicity) {
  std::vector<gra::mg5helas::HelicityComponent> components;
  components.reserve(4);
  constexpr auto labels = gra::spin::BinaryHelicityLabelsX2();
  for (const int upper : labels) {
    for (const int lower : labels) {
      components.push_back({{upper, lower}, {}, 0, upper == lower ? same_helicity : 0.0});
    }
  }
  return components;
}

// Build direct fermion or nested vector decays of an analytic scalar production point
gra::LORENTZSCALAR ScalarPhotonDecayPoint(const int pdg, const bool nuclear, const bool cascade) {
  auto lts = ScalarPhotonPoint(pdg, nuclear);
  auto &res = lts.process.ROOT_RES;
  res.p = lts.PDG.FindByPDG(pdg);
  lts.process.ROOT_RES_ACTIVE = true;
  lts.process.root_resonance_pdg = pdg;
  lts.process.root_decay_mode = gra::RootDecayMode::Physical;
  lts.process.SPINDEC = true;
  lts.PS_active = true;
  lts.amplitude.DECAY_SYM = false;
  lts.decay_structure = {gra::DecayType::JacobWickIncoherent, true};
  if (cascade) {
    auto seed = ContinuumVectorCascadeLTSForTest(false);
    auto &nested = seed.decaytree[0].legs[1];
    nested.p.width = 0.03;
    const auto pair = TwoBodyRestKinematics(nested.p4.M(), 0.0, 0.0, 0.49, -0.38);
    for (const auto &i : indices(pair)) {
      gra::MDecayBranch leaf;
      leaf.p = ToyParticle("f", i == 0 ? 900401 : -900401, 1, 0.0);
      leaf.p4 = BoostFromRestFrame(pair[i], nested.p4);
      nested.legs.push_back(leaf);
    }
    gra::spin::InitTwoBodyBasis(nested.hel, 0.0, 0.5, 0.5, {0.0}, {-0.5, 0.5}, {-0.5, 0.5}, "scalar fermion decay");
    nested.hel.T = {{1.0 / std::sqrt(2.0), 0.0}, {0.0, -1.0 / std::sqrt(2.0)}};
    nested.hel.g_decay = {0.67, -0.23};
    const double scale = lts.pfinal[0].M() / seed.pfinal[0].M();
    const auto scale_branch = [&](const auto &self, gra::MDecayBranch &branch) -> void {
      branch.p4 *= scale;
      branch.p.mass *= scale;
      branch.p.width *= scale;
      for (auto &leg : branch.legs) { self(self, leg); }
      if (!branch.legs.empty()) { branch.W_event = BranchTwoBodyPhaseSpaceForTest(branch); }
    };
    lts.decaytree = seed.decaytree;
    for (auto &branch : lts.decaytree) {
      scale_branch(scale_branch, branch);
      BoostDecayBranchForTest(branch, lts.pfinal[0], lts.pfinal[0].M());
    }
    gra::spin::InitTwoBodyBasis(res.hel_decay, 0.0, 1.0, 1.0, {0.0}, {-1.0, 0.0, 1.0}, {-1.0, 0.0, 1.0}, "scalar vector decay");
    res.hel_decay.T = gra::MMatrix<std::complex<double>>::IdentityMatrix(3) / std::sqrt(3.0);
    unsigned int mask = 0;
    SetMixedMassProposalsForTest(lts.decaytree, mask);
  } else {
    const auto pair = TwoBodyRestKinematics(lts.pfinal[0].M(), 0.0, 0.0, 0.61, -0.37);
    lts.decaytree.clear();
    for (const auto &i : indices(pair)) {
      gra::MDecayBranch leaf;
      leaf.p = ToyParticle("f", i == 0 ? 900401 : -900401, 1, 0.0);
      leaf.p4 = BoostFromRestFrame(pair[i], lts.pfinal[0]);
      lts.decaytree.push_back(leaf);
    }
    gra::spin::InitTwoBodyBasis(res.hel_decay, 0.0, 0.5, 0.5, {0.0}, {-0.5, 0.5}, {-0.5, 0.5}, "scalar fermion decay");
    res.hel_decay.T = {{1.0 / std::sqrt(2.0), 0.0}, {0.0, -1.0 / std::sqrt(2.0)}};
  }
  res.hel_decay.BR = 0.23;
  res.hel_decay.g_decay = std::polar(std::sqrt(2.0 * res.p.mass * res.p.width * res.hel_decay.BR * 8.0 * gra::math::PI), 0.41);
  lts.DW = gra::kinematics::MCW(RootCascadePhaseSpaceForTest(lts, lts.decaytree));
  return lts;
}

// Build one physical monopole pair point with proton or nuclear beams
gra::LORENTZSCALAR MonopolePhotonPoint(const bool nuclear) {
  gra::LORENTZSCALAR lts;
  lts.PDG                             = LoadedPDGTable();
  lts.process.root_decay_mode         = gra::RootDecayMode::Physical;
  lts.id1                             = gra::PDG::PDG_gamma;
  lts.id2                             = gra::PDG::PDG_gamma;

  double beam_mass = gra::PDG::mp;
  double beam_pz   = 6500.0;
  if (nuclear) {
    beam_mass          = 14.899;
    beam_pz            = 5000.0;
    lts.beam1.pdg      = gra::nuclear::EncodeNuclearPDG(16, 8);
    lts.beam1.mass     = beam_mass;
    lts.beam1.chargeX3 = 24;
    lts.beam1.spinX2   = 0;
    lts.upc_model      = OxygenUPC(gra::nuclear::CoherenceType::Coherent);
  } else {
    lts.beam1.pdg      = gra::PDG::PDG_p;
    lts.beam1.mass     = beam_mass;
    lts.beam1.chargeX3 = 3;
    lts.beam1.spinX2   = 1;
  }
  lts.beam2 = lts.beam1;

  constexpr double photon_energy = 2300.0;
  lts.q1                         = gra::M4Vec(0.18, -0.12, photon_energy, photon_energy);
  lts.q2                         = gra::M4Vec(-0.18, 0.12, -photon_energy, photon_energy);
  CompleteElasticForwardState(lts, beam_pz, beam_mass);
  lts.s      = (lts.pbeam1 + lts.pbeam2).M2();
  lts.s_hat  = lts.pfinal[0].M2();
  lts.m2     = lts.s_hat;
  lts.sqrt_s = std::sqrt(lts.s);

  const gra::MParticle monopole = lts.PDG.FindByPDG(gra::PDG::PDG_monopole);
  const auto           final    = PhotonPair(lts.pfinal[0].M(), monopole.mass);
  gra::MDecayBranch    particle;
  particle.p  = monopole;
  particle.p4 = final[0];
  gra::MDecayBranch antiparticle;
  antiparticle.p  = lts.PDG.FindByPDG(-gra::PDG::PDG_monopole);
  antiparticle.p4 = final[1];
  lts.decaytree   = {particle, antiparticle};
  lts.t_hat       = (lts.q1 - particle.p4).M2();
  lts.u_hat       = (lts.q1 - antiparticle.p4).M2();
  return lts;
}

// Build one closed non-collinear massless three-body state
std::array<gra::M4Vec, 3> PhotonThreeBody(double mass) {
  std::array<gra::M4Vec, 3> out;
  const double              energy = mass / 3.0;
  for (const auto &i : indices(out)) {
    const double phi = 2.0 * gra::math::PI * static_cast<double>(i) / 3.0;
    out[i]           = gra::M4Vec(energy * std::cos(phi), energy * std::sin(phi), 0.0, energy);
    out[i].RotateX(0.37);
    out[i].RotateY(-0.29);
  }
  return out;
}

// Build one physical colored gamma-gamma point with nuclear forward legs
gra::LORENTZSCALAR NuclearColoredPoint(const gra::nuclear::CoherenceType emission, const double qx) {
  constexpr double   mass = 14.899;
  gra::LORENTZSCALAR lts;
  lts.process.root_decay_mode         = gra::RootDecayMode::Physical;
  lts.beam1.pdg                       = gra::nuclear::EncodeNuclearPDG(16, 8);
  lts.beam1.mass                      = mass;
  lts.beam1.chargeX3                  = 24;
  lts.beam1.spinX2                    = 0;
  lts.beam2                           = lts.beam1;
  lts.id1                             = gra::PDG::PDG_gamma;
  lts.id2                             = gra::PDG::PDG_gamma;
  lts.alphaQCD                        = 0.118;
  lts.q1                              = gra::M4Vec(qx, -0.04, 100.0, 100.0);
  lts.q2                              = gra::M4Vec(-qx, 0.04, -100.0, 100.0);
  CompleteElasticForwardState(lts, 2000.0, mass);
  lts.pfinal[0]    = lts.q1 + lts.q2;
  lts.s            = (lts.pbeam1 + lts.pbeam2).M2();
  lts.s_hat        = lts.pfinal[0].M2();
  lts.m2           = lts.s_hat;
  const auto final = PhotonThreeBody(lts.pfinal[0].M());
  lts.decaytree    = {PhotonLeaf(2, final[0]), PhotonLeaf(-2, final[1]), PhotonLeaf(21, final[2])};
  lts.upc_model    = OxygenUPC(emission);
  return lts;
}

// Build one closed non-collinear massless four-body state
std::array<gra::M4Vec, 4> PhotonFourBody(double mass) {
  const double              first_energy  = 0.2 * mass;
  const double              second_energy = 0.3 * mass;
  std::array<gra::M4Vec, 4> out           = {
                gra::M4Vec(first_energy, 0.0, 0.0, first_energy), gra::M4Vec(-first_energy, 0.0, 0.0, first_energy),
                gra::M4Vec(0.0, second_energy, 0.0, second_energy), gra::M4Vec(0.0, -second_energy, 0.0, second_energy)};
  for (auto &momentum : out) {
    momentum.RotateX(0.31);
    momentum.RotateY(-0.47);
    momentum.RotateZ(0.83);
  }
  return out;
}

// Build one exact elastic photon event for a registered final state
gra::LORENTZSCALAR PhotonCovariancePoint(const std::string &name) {
  gra::LORENTZSCALAR lts;
  lts.process.root_decay_mode         = gra::RootDecayMode::Physical;
  lts.beam1.pdg                       = gra::PDG::PDG_p;
  lts.beam1.chargeX3                  = 3;
  lts.beam2                           = lts.beam1;
  lts.id1                             = gra::PDG::PDG_gamma;
  lts.id2                             = gra::PDG::PDG_gamma;
  lts.alphaQCD                        = 0.118;
  const double photon_energy          = name == "yy_ww" ? 200.0 : 300.0;
  lts.q1                              = gra::M4Vec(0.35, -0.22, photon_energy, photon_energy);
  lts.q2                              = gra::M4Vec(-0.35, 0.22, -photon_energy, photon_energy);
  CompleteElasticForwardState(lts, 6500.0);
  lts.pfinal[0]     = lts.q1 + lts.q2;
  lts.s             = (lts.pbeam1 + lts.pbeam2).M2();
  lts.s_hat         = lts.pfinal[0].M2();
  lts.m2            = lts.s_hat;
  const double mass = lts.pfinal[0].M();

  if (name == "yy_yy") {
    const auto final = PhotonPair(mass, 0.0);
    lts.decaytree    = {PhotonLeaf(22, final[0]), PhotonLeaf(22, final[1])};
  } else if (name == "yy_ll") {
    const AMP_MG5_yy_ll model;
    const auto final = PhotonPair(mass, model.Particles().at(11).mass);
    lts.decaytree    = {PhotonLeaf(-11, final[0]), PhotonLeaf(11, final[1])};
  } else if (name == "yy_uubarg") {
    const auto final = PhotonThreeBody(mass);
    lts.decaytree    = {PhotonLeaf(2, final[0]), PhotonLeaf(-2, final[1]), PhotonLeaf(21, final[2])};
  } else if (name == "yy_uubargg") {
    const auto final = PhotonFourBody(mass);
    lts.decaytree    = {PhotonLeaf(2, final[0]), PhotonLeaf(-2, final[1]), PhotonLeaf(21, final[2]),
                        PhotonLeaf(21, final[3])};
  } else if (name == "yy_jj") {
    const auto final = PhotonPair(mass, 0.0);
    lts.decaytree    = {PhotonLeaf(2, final[0]), PhotonLeaf(-2, final[1])};
  } else if (name == "yy_ww") {
    const auto parameters = gra::CreatePhotonMG5Process("MG5_YY_WW")->Particles();
    const double w_mass = parameters.at(24).mass;
    const auto        w_pair      = TwoBodyRestKinematics(mass, w_mass, w_mass, 0.93, -0.41);
    const auto        plus_decay  = TwoBodyRestKinematics(w_mass, parameters.at(11).mass, 0.0, 0.73, 0.29);
    const auto        minus_decay = TwoBodyRestKinematics(w_mass, parameters.at(13).mass, 0.0, 1.18, -0.61);
    gra::MDecayBranch wplus;
    wplus.p.pdg   = 24;
    wplus.p.mass  = w_mass;
    wplus.p.width = 2.085;
    wplus.p4      = w_pair[0];
    wplus.legs    = {PhotonLeaf(-11, BoostFromRestFrame(plus_decay[0], wplus.p4)),
                     PhotonLeaf(12, BoostFromRestFrame(plus_decay[1], wplus.p4))};
    gra::MDecayBranch wminus;
    wminus.p.pdg   = -24;
    wminus.p.mass  = w_mass;
    wminus.p.width = 2.085;
    wminus.p4      = w_pair[1];
    wminus.legs    = {PhotonLeaf(13, BoostFromRestFrame(minus_decay[0], wminus.p4)),
                      PhotonLeaf(-14, BoostFromRestFrame(minus_decay[1], wminus.p4))};
    lts.decaytree  = {wplus, wminus};
  } else if (name == "yy_zjj") {
    const auto        final = PhotonFourBody(mass);
    gra::MDecayBranch z;
    z.p.pdg       = 23;
    z.p4          = final[0] + final[1];
    const auto parameters = gra::CreatePhotonMG5Process("MG5_YY_ZJJ")->Particles();
    const auto decay = TwoBodyRestKinematics(z.p4.M(), parameters.at(13).mass, parameters.at(13).mass, 0.73, 0.29);
    z.legs = {PhotonLeaf(-13, BoostFromRestFrame(decay[0], z.p4)),
              PhotonLeaf(13, BoostFromRestFrame(decay[1], z.p4))};
    lts.decaytree = {z, PhotonLeaf(2, final[2]), PhotonLeaf(-2, final[3])};
  } else {
    throw std::invalid_argument("PhotonCovariancePoint: unknown process " + name);
  }
  return lts;
}

// Transform one complete decay branch including all intermediate momenta
template <class Transform>
void TransformPhotonBranch(gra::MDecayBranch &branch, const Transform &transform) {
  transform(branch.p4);
  for (auto &leg : branch.legs) { TransformPhotonBranch(leg, transform); }
}

// Transform every physical momentum entering the coherent photon source
template <class Transform>
gra::LORENTZSCALAR TransformPhotonPoint(gra::LORENTZSCALAR lts, const Transform &transform) {
  transform(lts.pbeam1);
  transform(lts.pbeam2);
  transform(lts.q1);
  transform(lts.q2);
  for (auto &momentum : lts.pfinal) { transform(momentum); }
  for (auto &branch : lts.decaytree) { TransformPhotonBranch(branch, transform); }
  return lts;
}

// Require scalar or source-contracted covariance for one photon amplitude
template <class Evaluate>
void RequirePhotonCovariance(const gra::LORENTZSCALAR &input, const Evaluate &evaluate, bool coherent,
                             const std::string &label) {
  gra::LORENTZSCALAR reference      = input;
  INFO(label);
  const double       reference_amp2 = evaluate(reference, coherent);
  REQUIRE(reference_amp2 > 0.0);

  const auto require_near = [&](gra::LORENTZSCALAR point, const std::string &transform) {
    CAPTURE(label, coherent, transform);
    const double amp2     = evaluate(point, coherent);
    const double relative = std::abs(amp2 - reference_amp2) / reference_amp2;
    CAPTURE(label, coherent, transform, reference_amp2, amp2, relative);
    REQUIRE(amp2 == Approx(reference_amp2).epsilon(3.0e-4));
  };

  require_near(TransformPhotonPoint(input,
                                    [](gra::M4Vec &momentum) {
                                      momentum.RotateX(0.41);
                                      momentum.RotateY(-0.72);
                                      momentum.RotateZ(1.19);
                                    }),
               "rotation");
  require_near(TransformPhotonPoint(input,
                                    [](gra::M4Vec &momentum) {
                                      momentum = momentum.LorentzBoost({0.21, -0.13, 0.31});
                                    }),
               "non-collinear boost");
  for (const double rapidity : {-2.0, 2.0}) {
    require_near(TransformPhotonPoint(input,
                                      [rapidity](gra::M4Vec &momentum) {
                                        momentum = momentum.LorentzBoost({0.0, 0.0, std::tanh(rapidity)});
                                      }),
                 "longitudinal rapidity " + std::to_string(rapidity));
  }
}

// Check one public generated photon family under the requested EPA treatment
template <class Amplitude>
void RequirePhotonFamilyCovariance(const std::string &name, bool coherent) {
  const auto evaluate = [&](gra::LORENTZSCALAR &lts, bool use_coherent) {
    Amplitude  amplitude;
    const auto result = amplitude.Evaluate(lts, 0.0, use_coherent);
    REQUIRE(result.Valid());
    REQUIRE(std::isfinite(result.amp2));
    return result.amp2;
  };
  RequirePhotonCovariance(PhotonCovariancePoint(name), evaluate, coherent, name);
}

// Construct an on-shell massive Dirac pair with a prescribed scattering angle
gra::LORENTZSCALAR MassivePhotonPair(double mass, double beta, double cosine, double phi) {
  gra::LORENTZSCALAR lts;
  const double energy = mass / std::sqrt(1.0 - beta * beta);
  const double momentum = beta * energy;
  const double transverse = momentum * std::sqrt(std::max(0.0, 1.0 - cosine * cosine));
  lts.q1 = gra::M4Vec(0.0, 0.0, energy, energy);
  lts.q2 = gra::M4Vec(0.0, 0.0, -energy, energy);
  lts.pbeam1 = lts.q1;
  lts.pbeam2 = lts.q2;
  lts.pfinal[0] = lts.q1 + lts.q2;
  lts.s = lts.s_hat = lts.m2 = 4.0 * energy * energy;
  lts.id1 = lts.id2 = gra::PDG::PDG_gamma;
  const gra::M4Vec first(transverse * std::cos(phi), transverse * std::sin(phi), momentum * cosine, energy);
  lts.decaytree = {PhotonLeaf(-11, first), PhotonLeaf(11, lts.pfinal[0] - first)};
  for (auto &branch : lts.decaytree) { branch.p.mass = mass; }
  return lts;
}

// Compute the spin-averaged massive Breit-Wheeler squared amplitude
// [REFERENCE: W. Zha et al., arXiv:1804.01813, Eqs. (7) and (8)]
double BreitWheelerAmp2(double beta, double cosine) {
  const double b2 = beta * beta;
  const double c2 = cosine * cosine;
  const double denominator = 1.0 - b2 * c2;
  return 4.0 * pow2(pow2(gra::qed::e_QED())) *
         ((1.0 + b2 * c2) / denominator + 2.0 * b2 * (1.0 - b2) * (1.0 - c2) / pow2(denominator));
}

}  // namespace

// Preserve fixed model masses at high boosts and reject resolved virtuality
TEST_CASE("MG5 external momenta must satisfy the generated model masses", "[MG5][massive][kinematics]") {
  const double electron = LoadedPDGTable().FindByPDG(11).mass;
  const gra::M4Vec boosted(0.0, 0.0, 1.0e6, 1.0e6);
  REQUIRE(gra::mg5::OnShellFinal({boosted}, {0.0, 0.0, electron}));
  REQUIRE(gra::mg5::OnShellFinal({boosted}, {0.0, 0.0, 0.0}));
  REQUIRE_FALSE(gra::mg5::OnShellFinal({gra::M4Vec(0.0, 0.0, 120.0, 130.0)}, {0.0, 0.0, 40.0}));
  REQUIRE(gra::mg5::OnShellFinal({gra::M4Vec(0.0, 0.0, 120.0, 130.0)}, {0.0, 0.0, 50.0}));
  REQUIRE_FALSE(gra::mg5::OnShellFinal({gra::M4Vec(0.0, 0.0, 130.0, 120.0)}, {0.0, 0.0, 40.0}));
  REQUIRE_FALSE(gra::mg5::OnShellFinal(
      {gra::M4Vec(0.0, 0.0, 0.0, std::numeric_limits<double>::infinity())}, {0.0, 0.0, electron}));
}

// Initialize a massive Dirac model through its actual independent parameters
void ConfigureDiracMass(AMP_MG5_yy_ll &amplitude, double mass) {
  SLHAReader card(gra::aux::ResolveProjectPath("MG5cards/Photon/yy_ll/param_card.dat"));
  card.set_block_entry("mass", 11, mass);
  card.set_block_entry("yukawa", 11, mass);
  amplitude.InitParameters(std::move(card));
  REQUIRE(amplitude.getMasses()[2] == Approx(mass).epsilon(1.0e-14));
  REQUIRE(amplitude.getMasses()[3] == Approx(mass).epsilon(1.0e-14));
}

// Test the physical fermion mass in both external states and exchanged propagators
TEST_CASE("Massive photon Dirac pairs reproduce the Breit-Wheeler differential cross section",
          "[MG5][gamma-gamma][massive][Breit-Wheeler]") {
  AMP_MG5_yy_ll amplitude;
  const auto &pdg = LoadedPDGTable();
  for (const double mass : {pdg.FindByPDG(11).mass, pdg.FindByPDG(13).mass, pdg.FindByPDG(15).mass, 500.0}) {
    ConfigureDiracMass(amplitude, mass);
    for (const double beta : {0.05, 0.4, 0.9, 0.999999, 1.0 - 1.0e-10}) {
      for (const double cosine : {-1.0, -0.9, 0.0, 0.3, 0.99999, 1.0 - 1.0e-12, 1.0}) {
        auto lts = MassivePhotonPair(mass, beta, cosine, 0.37);
        const auto result = amplitude.Evaluate(lts, 0.0, false);
        CAPTURE(mass, beta, cosine);
        REQUIRE(result.Valid());
        CHECK(result.amp2 == Approx(BreitWheelerAmp2(beta, cosine)).epsilon(2.0e-9));
      }
    }
  }
}

// Integrate the real massive amplitude over the complete two-body solid angle
// [REFERENCE: W. Zha et al., arXiv:1804.01813, Eq. (7)]
TEST_CASE("Massive photon Dirac pairs integrate to the Breit-Wheeler total cross section",
          "[MG5][gamma-gamma][massive][Breit-Wheeler]") {
  AMP_MG5_yy_ll amplitude;
  const auto [nodes, weights] = gra::math::GaussLegendreRule(160, -1.0, 1.0);
  const auto &pdg = LoadedPDGTable();
  for (const double mass : {pdg.FindByPDG(11).mass, pdg.FindByPDG(13).mass, pdg.FindByPDG(15).mass, 500.0}) {
    ConfigureDiracMass(amplitude, mass);
    for (const double beta : {0.05, 0.4, 0.9, 0.99}) {
      double integral = 0.0;
      for (const auto &i : indices(nodes)) {
        auto lts = MassivePhotonPair(mass, beta, nodes[i], 0.0);
        const auto result = amplitude.Evaluate(lts, 0.0, false);
        REQUIRE(result.Valid());
        integral += weights[i] * beta * result.amp2 / (32.0 * gra::math::PI * lts.s_hat);
      }
      const double b2 = beta * beta;
      const double s = 4.0 * mass * mass / (1.0 - b2);
      const double expected = 2.0 * gra::math::PI * pow2(gra::qed::alpha_QED()) / s *
                              ((3.0 - b2 * b2) * (std::log1p(beta) - std::log1p(-beta)) -
                               2.0 * beta * (2.0 - b2));
      CAPTURE(mass, beta);
      REQUIRE(integral == Approx(expected).epsilon(2.0e-9));
    }
  }
}

// Check massive spinor amplitudes under rotations, boosts, parity and photon exchange
TEST_CASE("Massive photon Dirac pairs preserve spacetime and beam symmetries",
          "[MG5][gamma-gamma][massive][Breit-Wheeler][covariance]") {
  AMP_MG5_yy_ll amplitude;
  ConfigureDiracMass(amplitude, LoadedPDGTable().FindByPDG(15).mass);
  for (const double beta : {0.05, 0.4, 0.9}) {
    const auto born = MassivePhotonPair(LoadedPDGTable().FindByPDG(15).mass, beta, 0.3, 0.37);
    const double expected = BreitWheelerAmp2(beta, 0.3);
    std::vector<gra::LORENTZSCALAR> points = {born};
    points.push_back(TransformPhotonPoint(born, [](gra::M4Vec &p) {
      p.RotateY(0.47);
      p.RotateZ(-0.81);
    }));
    points.push_back(TransformPhotonPoint(born, [](gra::M4Vec &p) { p = p.LorentzBoost({0.21, -0.13, 0.31}); }));
    points.push_back(TransformPhotonPoint(born, [](gra::M4Vec &p) { p = gra::M4Vec(-p.Px(), -p.Py(), -p.Pz(), p.E()); }));
    auto exchanged = born;
    std::swap(exchanged.q1, exchanged.q2);
    std::swap(exchanged.pbeam1, exchanged.pbeam2);
    points.push_back(exchanged);
    for (auto &point : points) {
      const auto result = amplitude.Evaluate(point, 0.0, false);
      REQUIRE(result.Valid());
      REQUIRE(result.amp2 == Approx(expected).epsilon(2.0e-9));
    }
  }
}

// Check real generated stable-leaf amplitudes with independently varied mass proposals
TEST_CASE("MG5 WW and Zjj decay amplitudes are independent of the cascade importance proposal",
          "[MG5][MProcess][cascade][proposal][physics]") {
  ToyHelicityProcess sampler;
  for (const std::string name : {"yy_ww", "yy_zjj"}) {
    const bool ww = name == "yy_ww";
    const std::string family = ww ? "MG5_YY_WW" : "MG5_YY_ZJJ";
    gra::MGeneratedPhotonProc process("yy", ww ? "WW" : "Zjj", {"generated", "MG5", "yy", "", 1}, family);
    auto amplitude = gra::CreatePhotonMG5Process(family);
    REQUIRE(amplitude != nullptr);
    auto point = PhotonCovariancePoint(name);
    point.PS_active = true;
    gra::SynchronizeMG5DecayParameters(point.decaytree, amplitude->Particles());
    for (auto &branch : point.decaytree) {
      if (branch.legs.empty()) { continue; }
      branch.W_event = BranchTwoBodyPhaseSpaceForTest(branch);
      branch.W = gra::kinematics::MCW(branch.W_event);
      branch.hel.g_decay = 0.0;
    }
    const auto structure = process.DecayStructureFor(point);
    REQUIRE(structure.type == gra::DecayType::Full);
    const auto baseline = amplitude->Evaluate(point, 0.0, false);
    REQUIRE(baseline.Valid());
    REQUIRE(baseline.amp2 > 0.0);
    const auto helicities = point.hamp;
    for (unsigned int selection = 0; selection < (ww ? 4U : 2U); ++selection) {
      CAPTURE(name, selection);
      auto lts = point;
      unsigned int mask = selection;
      SetMixedMassProposalsForTest(lts.decaytree, mask);
      const auto result = amplitude->Evaluate(lts, 0.0, false);
      REQUIRE(result.Valid());
      REQUIRE(result.amp2 == Approx(baseline.amp2).epsilon(1e-12));
      RequireVectorNear(lts.hamp, helicities, 1e-11);
      const double density = MixedMassDensityForTest(lts.decaytree);
      const double phase = sampler.CascadePhaseSpaceForTest(lts, process.DecayStructureFor(lts));
      REQUIRE(density * phase * result.amp2 ==
              Approx(InternalCascadePhaseSpaceForTest(lts.decaytree) * baseline.amp2).epsilon(1e-11));
    }
  }
}

TEST_CASE("Gamma monopole pairs use the magnetic amplitude in either order", "[gra::MGamma][process][monopole]") {
  gra::LORENTZSCALAR forward;
  forward.PDG               = LoadedPDGTable();
  const auto   monopole     = forward.PDG.FindByPDG(gra::PDG::PDG_monopole);
  const auto   antimonopole = forward.PDG.FindByPDG(-gra::PDG::PDG_monopole);
  const double energy       = 1.2 * monopole.mass;
  const double momentum     = std::sqrt(pow2(energy) - pow2(monopole.mass));
  const double cos_theta    = 0.31;
  const double sin_theta    = std::sqrt(1.0 - pow2(cos_theta));
  forward.q1                = gra::M4Vec(0.0, 0.0, energy, energy);
  forward.q2                = gra::M4Vec(0.0, 0.0, -energy, energy);
  gra::MDecayBranch particle;
  particle.p  = monopole;
  particle.p4 = gra::M4Vec(momentum * sin_theta, 0.0, momentum * cos_theta, energy);
  gra::MDecayBranch antiparticle;
  antiparticle.p    = antimonopole;
  antiparticle.p4   = gra::M4Vec(-momentum * sin_theta, 0.0, -momentum * cos_theta, energy);
  forward.decaytree = {particle, antiparticle};
  forward.s_hat     = pow2(2.0 * energy);
  forward.t_hat     = (forward.q1 - particle.p4).M2();
  forward.u_hat     = (forward.q1 - antiparticle.p4).M2();

  gra::LORENTZSCALAR reversed = forward;
  std::swap(reversed.decaytree[0], reversed.decaytree[1]);
  std::swap(reversed.t_hat, reversed.u_hat);
  const auto  definition = gra::MGamma::ProcessDefinitionFor(gra::MGammaMode::FermionPair);
  const auto  soft_model = gra::MModelTune::Load(modelfile);
  gra::MGamma forward_amplitude(forward, soft_model, definition);
  gra::MGamma reversed_amplitude(reversed, soft_model, definition);
  const auto  forward_result  = forward_amplitude.yyffbar(forward, false);
  const auto  reversed_result = reversed_amplitude.yyffbar(reversed, false);
  REQUIRE(forward_result.Valid());
  REQUIRE(reversed_result.Valid());
  REQUIRE(forward_result.amp2 > 0.0);
  REQUIRE(reversed_result.amp2 == Approx(forward_result.amp2).epsilon(1e-12));
  REQUIRE_FALSE(forward.hamp.empty());
  REQUIRE(gra::SquaredNorm(forward.hamp) == Approx(4.0 * forward_result.amp2).epsilon(1e-12));

  const gra::MParticle monopolium = forward.PDG.FindByPDG(gra::PDG::PDG_monopolium);
  const auto   parameters   = gra::ReadMonopoleParam(*soft_model, monopole.mass, monopolium.mass, monopolium.width);
  const double beta         = std::sqrt(1.0 - 4.0 * pow2(monopole.mass) / forward.s_hat);
  const double g            = 2.0 * gra::math::PI * parameters->gn / gra::qed::e_QED();
  const auto   pair_process = gra::amplitude::FindProcess("PHOTON", "yy_ll");
  REQUIRE(pair_process.has_value());
  auto electric_kernel = gra::CreatePhotonMG5Process(*pair_process);
  REQUIRE(electric_kernel != nullptr);
  SLHAReader electric_card(gra::aux::ResolveProjectPath("MG5cards/Photon/yy_ll/param_card.dat"));
  electric_card.set_block_entry("mass", 11, monopole.mass);
  electric_card.set_block_entry("yukawa", 11, monopole.mass);
  electric_kernel->InitParameters(std::move(electric_card));
  gra::LORENTZSCALAR electric = forward;
  electric.decaytree[0].p.pdg = -11;
  electric.decaytree[1].p.pdg = 11;
  const auto electric_result  = electric_kernel->Evaluate(electric, 0.0, false);
  REQUIRE(electric_result.Valid());
  const double alpha_ref =
      gra::mg5helas::AlphaQEDAtZero(forward) ? gra::qed::alpha_QED() : electric_kernel->DefaultAlphaQED();
  const double electric2      = 4.0 * gra::math::PI * alpha_ref;
  const double magnetic_scale = pow2(g * beta) / electric2;
  REQUIRE(forward_result.amp2 == Approx(electric_result.amp2 * pow2(magnetic_scale)).epsilon(2e-11));
  REQUIRE(forward.hamp.size() == electric.hamp.size());
  for (const auto &i : indices(forward.hamp)) {
    RequireComplexNear(forward.hamp[i], magnetic_scale * electric.hamp[i], 2e-11);
  }

  gra::LORENTZSCALAR malformed = forward;
  malformed.decaytree.clear();
  malformed.hamp              = {1.0};
  const auto malformed_result = forward_amplitude.yyffbar(malformed, false);
  CHECK(malformed_result.status == gra::mg5helas::EvaluationStatus::AmplitudeFailure);
  CHECK(gra::math::IsZero(malformed_result.amp2));
  CHECK(malformed.hamp.empty());
}

TEST_CASE("Monopole parameters construct safely from one tune across threads", "[gra::MGamma][params][threading]") {
  const gra::MParticle monopole   = LoadedPDGTable().FindByPDG(gra::PDG::PDG_monopole);
  const gra::MParticle monopolium = LoadedPDGTable().FindByPDG(gra::PDG::PDG_monopolium);
  const auto            model_tune = gra::MModelTune::Load(modelfile);
  gra::MModelCache      cache(model_tune);
  const auto            first  = gra::GetMonopoleParam(cache, monopole.mass, monopolium.mass, monopolium.width);
  const auto            second = gra::GetMonopoleParam(cache, monopole.mass, monopolium.mass, monopolium.width);
  REQUIRE(first == second);
  REQUIRE(first->initialized);
  REQUIRE(second->initialized);
  REQUIRE(first->monopole_mass == Approx(monopole.mass));
  REQUIRE(first->monopolium_mass == Approx(monopolium.mass));
  REQUIRE(first->monopolium_width == Approx(monopolium.width));
  REQUIRE(first->BindingEnergy() == Approx(monopolium.mass - 2.0 * monopole.mass));

  const double binding_fraction = 2.0 - monopolium.mass / monopole.mass;
  const double expected_psi =
      std::pow(binding_fraction, 3.0 / 4.0) * std::pow(monopole.mass, 3.0 / 2.0) / std::sqrt(gra::math::PI);
  REQUIRE(first->PsiAtOrigin() == Approx(expected_psi).epsilon(1e-12));

  const double alpha_g        = 2.5;
  const double expected_gamma = 32.0 * gra::math::PI * pow2(alpha_g) / pow2(monopolium.mass) * pow2(expected_psi);
  REQUIRE(first->GammaGamma(alpha_g) == Approx(expected_gamma).epsilon(1e-12));

  const auto changed_pole = gra::GetMonopoleParam(cache, monopole.mass, monopolium.mass + 1.0, monopolium.width);
  REQUIRE(changed_pole->monopolium_mass == Approx(first->monopolium_mass + 1.0));

  constexpr std::size_t                                   nthreads = 8;
  std::vector<std::shared_ptr<const gra::MMonopoleParam>> handles(nthreads);
  std::vector<std::thread>                                workers;
  workers.reserve(nthreads);

  for (std::size_t i = 0; i < nthreads; ++i) {
    workers.emplace_back([i, &handles, &monopole, &monopolium, &cache] {
      handles[i] = gra::GetMonopoleParam(cache, monopole.mass, monopolium.mass, monopolium.width);
    });
  }
  for (auto &worker : workers) { worker.join(); }

  for (const auto &handle : handles) { REQUIRE(handle == first); }
}

TEST_CASE("MMonopoleParam validates the physical monopolium model", "[gra::MGamma][params][validation]") {
  const gra::MParticle monopole   = LoadedPDGTable().FindByPDG(gra::PDG::PDG_monopole);
  const gra::MParticle monopolium = LoadedPDGTable().FindByPDG(gra::PDG::PDG_monopolium);
  const auto            model_tune = gra::MModelTune::Load(modelfile);

  const auto above_threshold =
      gra::ReadMonopoleParam(*model_tune, monopole.mass, 2.0 * monopole.mass, monopolium.width);
  REQUIRE_THROWS(above_threshold->ValidateMonopolium());

  const auto zero_width = gra::ReadMonopoleParam(*model_tune, monopole.mass, monopolium.mass, 0.0);
  REQUIRE_THROWS(zero_width->ValidateMonopolium());

  const auto tune = WriteModifiedPhotoVMTune("monopolium_bad_wavefunction",
                                             [](auto &j) { j.at("PARAM_MONOPOLE").at("wavefunction") = "UNKNOWN"; });
  REQUIRE_THROWS(
      gra::ReadMonopoleParam(*gra::MModelTune::Load(tune.second), monopole.mass, monopolium.mass, monopolium.width));
}

TEST_CASE("MGamma monopolium amplitude uses the PDG 881 physical pole", "[gra::MGamma][monopolium][physics]") {
  gra::LORENTZSCALAR lts = ScalarPhotonPoint(gra::PDG::PDG_monopolium, false);

  gra::MParticle monopolium = lts.PDG.FindByPDG(gra::PDG::PDG_monopolium);
  monopolium.mass += 12.0;
  monopolium.width += 3.0;
  lts.PDG.PDG_table[gra::PDG::PDG_monopolium] = monopolium;
  lts.decaytree[0].p                          = monopolium;

  const auto parameters =
      gra::ReadMonopoleParam(*gra::MModelTune::Load(modelfile), lts.PDG.FindByPDG(gra::PDG::PDG_monopole).mass,
                             monopolium.mass, monopolium.width);
  const double g         = 2.0 * gra::math::PI * parameters->gn / gra::qed::e_QED();
  const double beta      = std::sqrt(1.0 - pow2(monopolium.mass) / lts.s_hat);
  const double alpha_g   = pow2(beta * g) / (4.0 * gra::math::PI);
  const double gamma_yy  = parameters->GammaGamma(alpha_g);
  const double hard_amp2 = 16.0 * gra::math::PI * pow2(monopolium.mass) * gamma_yy * monopolium.width /
                           (pow2(lts.s_hat - pow2(monopolium.mass)) + pow2(monopolium.mass * monopolium.width));
  const std::complex<double> expected_helicity =
      std::sqrt(32.0 * gra::math::PI * pow2(monopolium.mass) * gamma_yy * monopolium.width) *
      gra::resonance::FixedWidthLineShape(lts.s_hat, monopolium.mass, monopolium.width);
  const auto hard = ScalarHardTensor(expected_helicity);
  REQUIRE(hard.size() == 4);
  RequireComplexNear(hard[0].value, expected_helicity, 1e-14);
  REQUIRE(std::norm(hard[1].value) == Approx(0.0));
  REQUIRE(std::norm(hard[2].value) == Approx(0.0));
  RequireComplexNear(hard[3].value, hard[0].value, 1e-14);
  REQUIRE((std::norm(hard[0].value) + std::norm(hard[3].value)) / 4.0 == Approx(hard_amp2).epsilon(1e-12));

  gra::LORENTZSCALAR          expected_lts = lts;
  std::vector<gra::M4Vec>     final        = {expected_lts.pfinal[0]};
  gra::mg5helas::EPAHardFrame frame;
  REQUIRE(gra::mg5helas::PrepareEPAHardFrame(expected_lts, final, frame));
  const auto   expected_hamp = gra::mg5helas::ContractEPAHardPhotonSources(expected_lts, hard, frame);
  const double expected_amp2 = gra::SquaredNorm(expected_hamp) / 4.0;

  gra::MGamma  gamma(lts, gra::MModelTune::Load(modelfile),
                     gra::MGamma::ProcessDefinitionFor(gra::MGammaMode::Monopolium));
  const double amp2 = gamma.yyMP(lts);

  REQUIRE(amp2 == Approx(expected_amp2).epsilon(1e-12));
  REQUIRE(lts.hamp.size() == expected_hamp.size());
  for (const auto &i : indices(lts.hamp)) { RequireComplexNear(lts.hamp[i], expected_hamp[i], 1e-12); }
  REQUIRE(std::abs(expected_helicity.imag()) > 0.0);
}

TEST_CASE("MGamma Higgs amplitude uses the active PDG and decay tables", "[gra::MGamma][Higgs][physics]") {
  ModelParamRestoreGuard restore;
  gra::MODELPARAM = "TUNE0";

  gra::LORENTZSCALAR lts   = ScalarPhotonPoint(25, false);
  gra::MParticle     higgs = lts.PDG.FindByPDG(25);

  higgs.mass += 0.8;
  higgs.width *= 1.7;
  lts.PDG.PDG_table[25] = higgs;
  lts.decaytree[0].p    = higgs;

  const double gamma_yy  = higgs.width * gra::resonance::GammaGammaBranchingRatio(higgs.pdg);
  const double hard_amp2 = 16.0 * gra::math::PI * pow2(higgs.mass) * gamma_yy * higgs.width /
                           (pow2(lts.s_hat - pow2(higgs.mass)) + pow2(higgs.mass * higgs.width));
  const std::complex<double> expected_helicity =
      std::sqrt(32.0 * gra::math::PI * pow2(higgs.mass) * gamma_yy * higgs.width) *
      gra::resonance::FixedWidthLineShape(lts.s_hat, higgs.mass, higgs.width);
  const auto hard = ScalarHardTensor(expected_helicity);
  REQUIRE(hard.size() == 4);
  RequireComplexNear(hard[0].value, expected_helicity, 1e-14);
  REQUIRE(std::norm(hard[1].value) == Approx(0.0));
  REQUIRE(std::norm(hard[2].value) == Approx(0.0));
  RequireComplexNear(hard[3].value, hard[0].value, 1e-14);
  REQUIRE((std::norm(hard[0].value) + std::norm(hard[3].value)) / 4.0 == Approx(hard_amp2).epsilon(1e-12));

  gra::LORENTZSCALAR          expected_lts = lts;
  std::vector<gra::M4Vec>     final        = {expected_lts.pfinal[0]};
  gra::mg5helas::EPAHardFrame frame;
  REQUIRE(gra::mg5helas::PrepareEPAHardFrame(expected_lts, final, frame));
  const auto   expected_hamp = gra::mg5helas::ContractEPAHardPhotonSources(expected_lts, hard, frame);
  const double expected_amp2 = gra::SquaredNorm(expected_hamp) / 4.0;

  gra::MGamma  gamma(lts, gra::MModelTune::Load(modelfile), gra::MGamma::ProcessDefinitionFor(gra::MGammaMode::Higgs));
  const double amp2 = gamma.yyHiggs(lts);

  REQUIRE(amp2 == Approx(expected_amp2).epsilon(1e-12));
  REQUIRE(lts.hamp.size() == expected_hamp.size());
  for (const auto &i : indices(lts.hamp)) { RequireComplexNear(lts.hamp[i], expected_hamp[i], 1e-12); }
  REQUIRE(std::abs(expected_helicity.imag()) > 0.0);

  gra::LORENTZSCALAR invalid;
  invalid.PDG   = lts.PDG;
  invalid.s_hat = lts.s_hat;
  gra::MGamma invalid_gamma(invalid, gra::MModelTune::Load(modelfile),
                            gra::MGamma::ProcessDefinitionFor(gra::MGammaMode::Higgs));
  REQUIRE(invalid_gamma.yyHiggs(invalid) == Approx(0.0));
  REQUIRE(invalid.hamp.empty());
}

// Check the inclusive-to-exclusive conversion with an independently normalized partial width
TEST_CASE("analytic scalar Jacob-Wick decays integrate to production times branching ratio",
          "[gra::MGamma][scalar][cascade][physics][normalization]") {
  const auto tune = gra::MModelTune::Load(modelfile);
  for (const int pdg : {25, gra::PDG::PDG_monopolium}) {
    const auto mode = pdg == 25 ? gra::MGammaMode::Higgs : gra::MGammaMode::Monopolium;
    for (const bool nuclear : {false, true}) {
      auto inclusive = ScalarPhotonPoint(pdg, nuclear);
      gra::MGamma gamma(inclusive, tune, gra::MGamma::ProcessDefinitionFor(mode));
      const double production = pdg == 25 ? gamma.yyHiggs(inclusive) : gamma.yyMP(inclusive);
      REQUIRE(production > 0.0);
      for (const bool spin : {false, true}) {
        CAPTURE(pdg, nuclear, spin);
        auto physical = ScalarPhotonDecayPoint(pdg, nuclear, false);
        physical.process.SPINDEC = spin;
        const double value = pdg == 25 ? gamma.yyHiggs(physical) : gamma.yyMP(physical);
        REQUIRE(value / (8.0 * gra::math::PI) == Approx(production * physical.process.ROOT_RES.hel_decay.BR).epsilon(2e-11));
        REQUIRE(physical.epa_hard.amplitude.size() == 16);
        REQUIRE(physical.hamp.size() == 4 * inclusive.hamp.size());
        physical.process.ROOT_RES.hel_decay.g_decay = 0.0;
        REQUIRE((pdg == 25 ? gamma.yyHiggs(physical) : gamma.yyMP(physical)) == Approx(0.0));
      }
    }
  }
}

// Exercise the real Higgs branching reader and its identical-photon normalization
TEST_CASE("Higgs physical decay initialization supplies the Jacob-Wick coupling once",
          "[gra::MGamma][Higgs][cascade][initialization][physics]") {
  auto inclusive = ScalarPhotonPoint(25, false);
  ToyHelicityProcess sampler;
  sampler.state.lts = inclusive;
  sampler.SetModelTune(sampler.GetModelTune());
  sampler.ProcPtr.Initialize("yy", "Higgs");
  auto &lts = sampler.state.lts;
  lts.process.root_decay_mode = gra::RootDecayMode::Physical;
  lts.process.root_resonance_pdg = sampler.ProcPtr.RootResonancePDG(lts);
  lts.decaytree.clear();
  for (const auto &momentum : PhotonPair(lts.pfinal[0].M(), 0.0)) {
    gra::MDecayBranch photon;
    photon.p = lts.PDG.FindByPDG(22);
    photon.p4 = BoostFromRestFrame(momentum, lts.pfinal[0]);
    lts.decaytree.push_back(photon);
  }
  auto setup = sampler.CreateProcessSetup();
  REQUIRE_NOTHROW(sampler.ProcPtr.InitializeAmplitude(setup));
  REQUIRE(lts.process.ROOT_RES_ACTIVE);
  REQUIRE_FALSE((lts.decay_structure.type == gra::DecayType::Full));
  REQUIRE((lts.decay_structure.type == gra::DecayType::JacobWickCoherent));
  REQUIRE(sampler.state.symmetry_factor == Approx(2.0));
  gra::MGamma gamma(lts, sampler.GetModelTune(), gra::MGamma::ProcessDefinitionFor(gra::MGammaMode::Higgs));
  const double physical = gamma.yyHiggs(lts);
  const double production = gamma.yyHiggs(inclusive);
  REQUIRE(physical / (16.0 * gra::math::PI) == Approx(production * lts.process.ROOT_RES.hel_decay.BR).epsilon(2e-11));

  // A worker copy preserves the complete EPA normalization and resolves changed angular switches
  const double reference = sampler.ProcPtr.GetBareAmplitude2(lts);
  REQUIRE(reference > 0.0);
  auto worker = sampler;
  for (const bool spindec : {false, true}) {
    for (const bool coherent : {false, true}) {
      worker.SetSPINDEC(spindec);
      worker.state.lts.amplitude.DECAY_SYM = coherent;
      const auto expected = spindec && coherent ? gra::DecayType::JacobWickCoherent : gra::DecayType::JacobWickIncoherent;
      REQUIRE(worker.ProcPtr.DecayStructureFor(worker.state.lts).type == expected);
      REQUIRE(worker.ProcPtr.GetBareAmplitude2(worker.state.lts) == Approx(reference).epsilon(2e-11));
      REQUIRE(worker.state.lts.decay_structure.type == expected);
      REQUIRE(lts.decay_structure.type == gra::DecayType::JacobWickCoherent);
    }
  }
}

// Compare generic nested helicity dressing with the complete physical decay matrix
TEST_CASE("analytic scalar nested Jacob-Wick decays reverse every mass proposal once",
          "[gra::MGamma][scalar][cascade][proposal][physics]") {
  const auto tune = gra::MModelTune::Load(modelfile);
  ToyHelicityProcess sampler;
  for (const int pdg : {25, gra::PDG::PDG_monopolium}) {
    const auto mode = pdg == 25 ? gra::MGammaMode::Higgs : gra::MGammaMode::Monopolium;
    auto inclusive = ScalarPhotonPoint(pdg, false);
    gra::MGamma gamma(inclusive, tune, gra::MGamma::ProcessDefinitionFor(mode));
    const double production = pdg == 25 ? gamma.yyHiggs(inclusive) : gamma.yyMP(inclusive);
    for (const bool coherent : {false, true}) {
      auto physical = ScalarPhotonDecayPoint(pdg, false, true);
      physical.amplitude.DECAY_SYM = coherent;
      auto complete = physical;
      complete.decay_structure = {gra::DecayType::Full};
      const auto &res = complete.process.ROOT_RES;
      const auto decay = gra::spin::ResonanceDecayMatrix(complete, res, "CM") * res.hel_decay.g_decay;
      const double target = production * decay.FrobNorm2() / (2.0 * res.p.mass * res.p.width);
      REQUIRE(target > 0.0);
      for (const bool mixture : {false, true}) {
        if (mixture && !coherent) { continue; }
        for (unsigned int selection = 0; selection < 8; ++selection) {
          CAPTURE(pdg, coherent, mixture, selection);
          auto lts = physical;
          lts.decay_symmetry_proposal_active = mixture;
          unsigned int mask = selection;
          SetMixedMassProposalsForTest(lts.decaytree, mask);
          const double density = mixture ? MixedHistoryDensityForTest(lts) : MixedMassDensityForTest(lts.decaytree);
          const double amp2 = pdg == 25 ? gamma.yyHiggs(lts) : gamma.yyMP(lts);
          const double phase = sampler.CascadePhaseSpaceForTest(lts, lts.decay_structure);
          REQUIRE(amp2 * phase * density == Approx(target * InternalCascadePhaseSpaceForTest(lts.decaytree)).epsilon(2e-9));
        }
      }
    }
  }
}

// Keep all physical final helicity columns in the immutable Born tensor
TEST_CASE("scalar Jacob-Wick cascades preserve screening contraction and rotational covariance",
          "[gra::MGamma][scalar][cascade][EPA][screening][covariance]") {
  const auto tune = gra::MModelTune::Load(modelfile);
  for (const int pdg : {25, gra::PDG::PDG_monopolium}) {
    const auto mode = pdg == 25 ? gra::MGammaMode::Higgs : gra::MGammaMode::Monopolium;
    for (const bool nuclear : {false, true}) {
      CAPTURE(pdg, nuclear);
      auto born = ScalarPhotonDecayPoint(pdg, nuclear, true);
      born.amplitude.DECAY_SYM = true;
      born.decay_symmetry_proposal_active = true;
      gra::MGamma gamma(born, tune, gra::MGamma::ProcessDefinitionFor(mode));
      const double value = pdg == 25 ? gamma.yyHiggs(born) : gamma.yyMP(born);
      REQUIRE(value > 0.0);
      REQUIRE(born.epa_hard.Ready());
      auto shifted = born;
      REQUIRE(ShiftEPAHardSources(shifted, 0.017, -0.011));
      const auto fast = gra::mg5helas::ContractEPAHard(shifted, born.epa_hard);
      REQUIRE(fast.Valid());
      auto direct = shifted;
      const double exact = pdg == 25 ? gamma.yyHiggs(direct) : gamma.yyMP(direct);
      REQUIRE(fast.amp2 == Approx(exact).epsilon(2e-10));
      REQUIRE(shifted.hamp.size() == direct.hamp.size());
      for (const auto &i : indices(shifted.hamp)) { RequireComplexNear(shifted.hamp[i], direct.hamp[i], 2e-10); }
      // Test coherent rotational covariance with the isotropic mean nuclear density
      auto coherent = born;
      if (nuclear) {
        const auto& upc = *born.upc_model;
        coherent.upc_model = std::make_shared<const gra::nuclear::MUPC>(*upc.Nucleus(1), *upc.Nucleus(2), upc.Param());
      }
      const double coherent_value = pdg == 25 ? gamma.yyHiggs(coherent) : gamma.yyMP(coherent);
      REQUIRE(coherent_value > 0.0);
      auto rotated = TransformPhotonPoint(coherent, [](gra::M4Vec &p) {
        p.RotateX(0.41);
        p.RotateY(-0.72);
        p.RotateZ(1.19);
      });
      const double rotated_value = pdg == 25 ? gamma.yyHiggs(rotated) : gamma.yyMP(rotated);
      REQUIRE(rotated_value == Approx(coherent_value).epsilon(2e-9));
    }
  }
}

TEST_CASE("Prepared scalar EPA tensors reproduce shifted proton and nuclear sources",
          "[gra::MGamma][scalar][EPA][screening]") {
  const auto model_tune = gra::MModelTune::Load(modelfile);
  for (const int pdg : {25, gra::PDG::PDG_monopolium}) {
    const auto mode       = pdg == 25 ? gra::MGammaMode::Higgs : gra::MGammaMode::Monopolium;
    const auto definition = gra::MGamma::ProcessDefinitionFor(mode);
    for (const bool nuclear : {false, true}) {
      CAPTURE(pdg, nuclear);
      gra::LORENTZSCALAR born = ScalarPhotonPoint(pdg, nuclear);
      gra::MGamma        prepared(born, model_tune, definition);
      const double       born_amp2 = pdg == 25 ? prepared.yyHiggs(born) : prepared.yyMP(born);
      REQUIRE(born_amp2 > 0.0);
      REQUIRE(born.epa_hard.Ready());
      REQUIRE(gra::mg5helas::UsesProtonEPAHardSources(born) == !nuclear);
      if (!nuclear) {
        REQUIRE(born.hamp.metadata.spin_basis == gra::ScreeningSpinBasis::ProtonHelicity);
        REQUIRE(born.hamp.metadata.spin_rows == 16);
      }

      gra::LORENTZSCALAR shifted = born;
      REQUIRE(ShiftEPAHardSources(shifted, 0.017, -0.011));
      const auto fast = gra::mg5helas::ContractEPAHard(shifted, born.epa_hard);
      REQUIRE(fast.Valid());

      gra::LORENTZSCALAR direct_lts = shifted;
      gra::MGamma        direct(direct_lts, model_tune, definition);
      const double       exact = pdg == 25 ? direct.yyHiggs(direct_lts) : direct.yyMP(direct_lts);
      REQUIRE(exact > 0.0);
      REQUIRE(fast.amp2 == Approx(exact).epsilon(2.0e-11));
      REQUIRE(shifted.hamp.size() == direct_lts.hamp.size());
      for (const auto &i : indices(shifted.hamp)) { RequireComplexNear(shifted.hamp[i], direct_lts.hamp[i], 2.0e-10); }
    }
  }
}

TEST_CASE("Prepared monopole-pair EPA tensors preserve magnetic helicity amplitudes",
          "[gra::MGamma][monopole][EPA][screening]") {
  const auto definition = gra::MGamma::ProcessDefinitionFor(gra::MGammaMode::FermionPair);
  const auto model_tune = gra::MModelTune::Load(modelfile);
  for (const bool nuclear : {false, true}) {
    CAPTURE(nuclear);
    gra::LORENTZSCALAR born = MonopolePhotonPoint(nuclear);
    gra::MGamma        prepared(born, model_tune, definition);
    const auto         born_result = prepared.yyffbar(born, true);
    REQUIRE(born_result.Valid());
    REQUIRE(born_result.amp2 > 0.0);
    REQUIRE(born.epa_hard.Ready());
    REQUIRE(gra::mg5helas::UsesProtonEPAHardSources(born) == !nuclear);
    if (!nuclear) {
      REQUIRE(born.hamp.metadata.spin_basis == gra::ScreeningSpinBasis::ProtonHelicity);
      REQUIRE(born.hamp.metadata.spin_rows == 16);
    }

    gra::LORENTZSCALAR shifted = born;
    REQUIRE(ShiftEPAHardSources(shifted, -0.013, 0.019));
    const auto fast = gra::mg5helas::ContractEPAHard(shifted, born.epa_hard);
    REQUIRE(fast.Valid());

    gra::LORENTZSCALAR direct_lts = shifted;
    gra::MGamma        direct(direct_lts, model_tune, definition);
    const auto         exact = direct.yyffbar(direct_lts, true);
    REQUIRE(exact.Valid());
    REQUIRE(fast.amp2 == Approx(exact.amp2).epsilon(2.0e-11));
    REQUIRE(shifted.hamp.size() == direct_lts.hamp.size());
    for (const auto &i : indices(shifted.hamp)) { RequireComplexNear(shifted.hamp[i], direct_lts.hamp[i], 2.0e-10); }
  }
}

TEST_CASE("Scalar diphoton hard tensor is covariant under rotations", "[gra::MGamma][Higgs][covariance]") {
  const auto               definition = gra::MGamma::ProcessDefinitionFor(gra::MGammaMode::Higgs);
  const auto               model_tune = gra::MModelTune::Load(modelfile);
  const gra::LORENTZSCALAR input      = ScalarPhotonPoint(25, false);
  gra::LORENTZSCALAR       reference  = input;
  gra::MGamma              reference_gamma(reference, model_tune, definition);
  const double             reference_amp2 = reference_gamma.yyHiggs(reference);
  REQUIRE(reference_amp2 > 0.0);

  gra::LORENTZSCALAR rotated = TransformPhotonPoint(input, [](gra::M4Vec &momentum) {
    momentum.RotateX(0.41);
    momentum.RotateY(-0.72);
    momentum.RotateZ(1.19);
  });
  gra::MGamma        rotated_gamma(rotated, model_tune, definition);
  const double       rotated_amp2 = rotated_gamma.yyHiggs(rotated);
  REQUIRE(rotated_amp2 == Approx(reference_amp2).epsilon(2.0e-8));
}

TEST_CASE("Generated photon registry covers exact colored 2-to-N amplitudes",
          "[MG2GRA][gamma-gamma][registry][color]") {
  // Build one stable generated external leg
  auto stable_branch = [](int pdg, const gra::M4Vec &p4) {
    gra::MDecayBranch branch;
    branch.p.pdg = pdg;
    branch.p4    = p4;
    return branch;
  };

  gra::LORENTZSCALAR lts;
  lts.q1                              = gra::M4Vec(0.3, 0.2, 150.0, 150.0);
  lts.q2                              = gra::M4Vec(-0.3, -0.2, -150.0, 150.0);
  CompleteElasticForwardState(lts, 6500.0);
  // Match the three massless partons to the exact elastic transfer energy
  const double energy3 = lts.pfinal[0].E() / 3.0;
  const double py      = 0.5 * energy3 * std::sqrt(3.0);
  lts.s_hat            = lts.pfinal[0].M2();
  lts.decaytree = {stable_branch(2, gra::M4Vec(energy3, 0.0, 0.0, energy3)),
                   stable_branch(-2, gra::M4Vec(-0.5 * energy3, py, 0.0, energy3)),
                   stable_branch(21, gra::M4Vec(-0.5 * energy3, -py, 0.0, energy3))};

  const auto definition = gra::amplitude::FindProcess("PHOTON", "yy_uubarg");
  REQUIRE(definition.has_value());
  REQUIRE(gra::HasPhotonMG5Process(*definition));
  auto matrix_element = gra::CreatePhotonMG5Process(*definition);
  REQUIRE(matrix_element != nullptr);
  REQUIRE(matrix_element->Processes() == std::vector<gra::amplitude::Process>{*definition});
  REQUIRE(matrix_element->SubprocessCount() == 1);

  const auto incoherent = matrix_element->Evaluate(lts, 0.118, false);
  INFO("photon registry status=" << static_cast<int>(incoherent.status));
  REQUIRE(incoherent.Valid());
  auto incoherent_amplitudes = lts.hamp;
  gra::Scale(incoherent_amplitudes, 0.5);
  const auto &incoherent_flows = lts.hard_color_flows;
  REQUIRE_FALSE(incoherent_amplitudes.empty());
  REQUIRE(incoherent_flows.size() == 1);
  REQUIRE(incoherent_flows[0].amplitudes.size() == incoherent_amplitudes.size());
  double registry_amp2 = 0.0;
  for (std::size_t i = 0; i < incoherent_amplitudes.size(); ++i) {
    registry_amp2 += std::norm(incoherent_amplitudes[i]);
    RequireComplexNear(incoherent_flows[0].amplitudes[i], incoherent_amplitudes[i], 1e-12);
  }
  REQUIRE(registry_amp2 > 0.0);
  double registry_hamp2 = 0.0;
  for (const auto &amplitude : lts.hamp) { registry_hamp2 += std::norm(amplitude); }
  REQUIRE(registry_hamp2 == Approx(4.0 * registry_amp2).epsilon(1e-12));

  auto direct = gra::CreatePhotonMG5Process(*definition);
  REQUIRE(direct != nullptr);
  gra::LORENTZSCALAR direct_lts    = lts;
  const auto         direct_result = direct->Evaluate(direct_lts, 0.118, false);
  REQUIRE(direct_result.Valid());
  const double direct_amp2 = direct_result.amp2;
  REQUIRE(registry_amp2 == Approx(direct_amp2).epsilon(1e-12));
  double direct_hamp2 = 0.0;
  for (const auto &amplitude : direct_lts.hamp) { direct_hamp2 += std::norm(amplitude); }
  REQUIRE(direct_hamp2 == Approx(4.0 * direct_amp2).epsilon(1e-12));

  gra::MRandom color_random;
  color_random.SetSeed(11);
  REQUIRE(matrix_element->SampleColorFlow(lts, color_random));
  REQUIRE(lts.decaytree[0].p.color_flow.flow1 != 0);
  REQUIRE(lts.decaytree[1].p.color_flow.flow2 != 0);
  REQUIRE_FALSE(lts.decaytree[2].p.color_flow.empty());

  // Reject a changed topology before reusing the prepared exact color flow
  gra::LORENTZSCALAR wrong_color_tree = lts;
  wrong_color_tree.decaytree.push_back(stable_branch(21, gra::M4Vec(0.0, 0.0, 0.0, 0.0)));
  REQUIRE_FALSE(matrix_element->SampleColorFlow(wrong_color_tree, color_random));

  const auto coherent = matrix_element->Evaluate(lts, 0.118, true);
  REQUIRE(coherent.Valid());
  const double coherent_amp2 = coherent.amp2;
  REQUIRE(coherent_amp2 > 0.0);
  double coherent_hamp2 = 0.0;
  for (const auto &amplitude : lts.hamp) { coherent_hamp2 += std::norm(amplitude); }
  REQUIRE(coherent_hamp2 == Approx(4.0 * coherent_amp2).epsilon(1e-12));

  gra::PROC_003_QED_YY_EPA kt_epa;
  gra::LORENTZSCALAR       event_lts = lts;
  kt_epa.BindModelTune(event_lts.model_cache->TunePtr());
  kt_epa.InitializeWorkerAmplitude(event_lts);
  REQUIRE(kt_epa.GammaGammaCON(event_lts, true) == Approx(coherent_amp2).epsilon(1e-12));
  gra::PROC_030_QED_YY_DZ_EPA dz_epa;
  event_lts = lts;
  dz_epa.BindModelTune(event_lts.model_cache->TunePtr());
  dz_epa.InitializeWorkerAmplitude(event_lts);
  REQUIRE(dz_epa.GammaGammaCON(event_lts, false) == Approx(registry_amp2).epsilon(1e-12));
  gra::PROC_040_QED_YY_LUX_EPA lux_epa;
  event_lts = lts;
  lux_epa.BindModelTune(event_lts.model_cache->TunePtr());
  lux_epa.InitializeWorkerAmplitude(event_lts);
  REQUIRE(lux_epa.GammaGammaCON(event_lts, false) == Approx(registry_amp2).epsilon(1e-12));

  // Exercise a full-rank two-dimensional SU(3) basis at 2-to-4
  const double       energy4   = lts.pfinal[0].E() / 4.0;
  const double       component = energy4 / std::sqrt(3.0);
  gra::LORENTZSCALAR four_body = lts;
  four_body.decaytree          = {stable_branch(2, gra::M4Vec(component, component, component, energy4)),
                                  stable_branch(-2, gra::M4Vec(component, -component, -component, energy4)),
                                  stable_branch(21, gra::M4Vec(-component, component, -component, energy4)),
                                  stable_branch(21, gra::M4Vec(-component, -component, component, energy4))};

  const auto four_definition = gra::amplitude::FindProcess("PHOTON", "yy_uubargg");
  REQUIRE(four_definition.has_value());
  REQUIRE(gra::HasPhotonMG5Process(*four_definition));
  auto four_process = gra::CreatePhotonMG5Process(*four_definition);
  REQUIRE(four_process != nullptr);
  REQUIRE(four_process->Processes() == std::vector<gra::amplitude::Process>{*four_definition});
  REQUIRE(four_process->SubprocessCount() == 1);

  const auto four_incoherent = four_process->Evaluate(four_body, 0.118, false);
  REQUIRE(four_incoherent.Valid());
  auto four_amplitudes = four_body.hamp;
  gra::Scale(four_amplitudes, 0.5);
  const auto &four_flows = four_body.hard_color_flows;
  REQUIRE(four_flows.size() == 2);
  double four_registry_amp2 = 0.0;
  for (std::size_t i = 0; i < four_amplitudes.size(); ++i) {
    four_registry_amp2 += std::norm(four_amplitudes[i]);
    RequireComplexNear(four_flows[0].amplitudes[i] + four_flows[1].amplitudes[i], four_amplitudes[i], 1e-12);
  }
  REQUIRE(four_registry_amp2 > 0.0);
  double four_registry_hamp2 = 0.0;
  for (const auto &amplitude : four_body.hamp) { four_registry_hamp2 += std::norm(amplitude); }
  REQUIRE(four_registry_hamp2 == Approx(4.0 * four_registry_amp2).epsilon(1e-12));

  auto four_direct = gra::CreatePhotonMG5Process(*four_definition);
  REQUIRE(four_direct != nullptr);
  gra::LORENTZSCALAR four_direct_lts    = four_body;
  const auto         four_direct_result = four_direct->Evaluate(four_direct_lts, 0.118, false);
  REQUIRE(four_direct_result.Valid());
  const double four_direct_amp2 = four_direct_result.amp2;
  REQUIRE(four_registry_amp2 == Approx(four_direct_amp2).epsilon(1e-12));
  double four_direct_hamp2 = 0.0;
  for (const auto &amplitude : four_direct_lts.hamp) { four_direct_hamp2 += std::norm(amplitude); }
  REQUIRE(four_direct_hamp2 == Approx(4.0 * four_direct_amp2).epsilon(1e-12));

  gra::LORENTZSCALAR isolated_four_body      = four_body;
  isolated_four_body.process.root_decay_mode = gra::RootDecayMode::Isolated;
  const auto isolated_evaluation             = four_process->Evaluate(isolated_four_body, 0.118, false);
  REQUIRE(isolated_evaluation.Valid());
  const double isolated_registry_amp2 = isolated_evaluation.amp2;
  REQUIRE(2.0 * isolated_registry_amp2 == Approx(four_registry_amp2).epsilon(1e-12));

  gra::LORENTZSCALAR isolated_direct_lts    = isolated_four_body;
  const auto         isolated_direct_result = four_direct->Evaluate(isolated_direct_lts, 0.118, false);
  REQUIRE(isolated_direct_result.Valid());
  const double isolated_direct_amp2 = isolated_direct_result.amp2;
  REQUIRE(isolated_direct_amp2 == Approx(isolated_registry_amp2).epsilon(1e-12));

  const auto four_coherent = four_process->Evaluate(four_body, 0.118, true);
  REQUIRE(four_coherent.Valid());
  const double four_coherent_amp2  = four_coherent.amp2;
  double       four_coherent_hamp2 = 0.0;
  for (const auto &amplitude : four_body.hamp) { four_coherent_hamp2 += std::norm(amplitude); }
  REQUIRE(four_coherent_hamp2 == Approx(4.0 * four_coherent_amp2).epsilon(1e-12));
  four_direct_lts                        = four_body;
  const auto four_direct_coherent_result = four_direct->Evaluate(four_direct_lts, 0.118, true);
  REQUIRE(four_direct_coherent_result.Valid());
  const double four_direct_coherent = four_direct_coherent_result.amp2;
  REQUIRE(four_coherent_amp2 == Approx(four_direct_coherent).epsilon(1e-12));

  gra::PROC_030_QED_YY_DZ_EPA four_epa;
  event_lts = four_body;
  four_epa.BindModelTune(event_lts.model_cache->TunePtr());
  four_epa.InitializeWorkerAmplitude(event_lts);
  REQUIRE(four_epa.GammaGammaCON(event_lts, false) == Approx(four_registry_amp2).epsilon(1e-12));
  gra::MRandom photon_rng;
  photon_rng.SetSeed(17);
  four_epa.BindRandom(photon_rng);
  four_epa.SampleColorFlow(event_lts);
  for (const auto &branch : event_lts.decaytree) { REQUIRE_FALSE(branch.p.color_flow.empty()); }
}

TEST_CASE("Generated photon matrix-element norms are Lorentz scalars", "[MG2GRA][gamma-gamma][covariance][scalar]") {
  for (const auto &definition : gra::amplitude::Processes("PHOTON")) {
    CAPTURE(definition.process_name);
    REQUIRE(gra::HasPhotonMG5Process(definition));
    auto amplitude = gra::CreatePhotonMG5Process(definition);
    REQUIRE(amplitude != nullptr);
    const auto evaluate = [&](gra::LORENTZSCALAR &lts, bool coherent) {
      const auto result = amplitude->Evaluate(lts, 0.118, coherent);
      REQUIRE(result.Valid());
      REQUIRE(gra::AllFinite(lts.hamp));
      return result.amp2;
    };
    RequirePhotonCovariance(PhotonCovariancePoint(definition.process_name), evaluate, false, definition.process_name);
  }

  gra::AMP_MG5_yy_zjj amplitude;
  const auto          evaluate = [&](gra::LORENTZSCALAR &lts, bool coherent) {
    const auto result = amplitude.Evaluate(lts, 0.118, coherent);
    CAPTURE(static_cast<int>(result.status), lts.epa_hard.frame.valid);
    REQUIRE(result.Valid());
    REQUIRE(std::isfinite(result.amp2));
    return result.amp2;
  };
  RequirePhotonCovariance(PhotonCovariancePoint("yy_zjj"), evaluate, false, "yy_zjj");
  RequirePhotonFamilyCovariance<gra::AMP_MG5_yy_jj>("yy_jj", false);
  RequirePhotonFamilyCovariance<gra::AMP_MG5_yy_ww>("yy_ww", false);
}

TEST_CASE("Generated colored photon flows share ordered EPA sector weights",
          "[MG2GRA][gamma-gamma][EPA][color][nuclear]") {
  const auto definition = gra::amplitude::FindProcess("PHOTON", "yy_uubargg");
  REQUIRE(definition.has_value());

  gra::LORENTZSCALAR lts       = PhotonCovariancePoint("yy_uubargg");
  auto               amplitude = gra::CreatePhotonMG5Process(*definition);
  REQUIRE(amplitude != nullptr);
  REQUIRE(amplitude->Evaluate(lts, 0.118, false).Valid());
  const auto raw_hamp  = lts.hamp;
  const auto raw_flows = lts.hard_color_flows;
  REQUIRE(raw_flows.size() == 2);

  gra::flux::EPASectorWeights inclusive;
  inclusive.amplitude_fraction = {std::sqrt(0.1), std::sqrt(0.2), std::sqrt(0.3), std::sqrt(0.4)};
  constexpr double common      = 1.3;
  REQUIRE(gra::flux::ApplyEPAAmplitudeWeights(lts.hamp, common, inclusive));
  for (auto &flow : lts.hard_color_flows) {
    REQUIRE(gra::flux::ApplyEPAAmplitudeWeights(flow.amplitudes, common, inclusive));
  }
  REQUIRE(lts.hamp.size() == raw_hamp.size() * inclusive.amplitude_fraction.size());

  const auto &weighted_flows = lts.hard_color_flows;
  REQUIRE(weighted_flows.size() == raw_flows.size());
  for (const auto &flow : indices(weighted_flows)) {
    REQUIRE(weighted_flows[flow].amplitudes.size() == lts.hamp.size());
    for (const auto &hard : indices(raw_flows[flow].amplitudes)) {
      for (const auto &sector : indices(inclusive.amplitude_fraction)) {
        const std::size_t row = hard * inclusive.amplitude_fraction.size() + sector;
        RequireComplexNear(weighted_flows[flow].amplitudes[row],
                           common * inclusive.amplitude_fraction[sector] * raw_flows[flow].amplitudes[hard], 1.0e-12);
      }
    }
  }

  std::vector<double> probabilities(weighted_flows.size(), 0.0);
  for (const auto &flow : indices(weighted_flows)) {
    probabilities[flow] = gra::SquaredNorm(weighted_flows[flow].amplitudes);
    REQUIRE(std::isfinite(probabilities[flow]));
    REQUIRE(probabilities[flow] >= 0.0);
  }
  const double total = gra::Sum(probabilities);
  REQUIRE(total > 0.0);
  for (double &probability : probabilities) { probability /= total; }
  REQUIRE(gra::Sum(probabilities) == Approx(1.0).epsilon(1.0e-14));

  auto repeated = gra::CreatePhotonMG5Process(*definition);
  REQUIRE(repeated != nullptr);
  gra::LORENTZSCALAR repeated_lts = PhotonCovariancePoint("yy_uubargg");
  REQUIRE(repeated->Evaluate(repeated_lts, 0.118, false).Valid());
  for (auto &flow : repeated_lts.hard_color_flows) {
    REQUIRE(gra::flux::ApplyEPAAmplitudeWeights(flow.amplitudes, common, inclusive));
  }
  const auto &repeated_flows = repeated_lts.hard_color_flows;
  REQUIRE(repeated_flows.size() == weighted_flows.size());
  for (const auto &flow : indices(repeated_flows)) {
    REQUIRE(repeated_flows[flow].amplitudes.size() == weighted_flows[flow].amplitudes.size());
    for (const auto &row : indices(repeated_flows[flow].amplitudes)) {
      RequireComplexNear(repeated_flows[flow].amplitudes[row], weighted_flows[flow].amplitudes[row], 1.0e-12);
    }
  }

  auto coherent = gra::CreatePhotonMG5Process(*definition);
  REQUIRE(coherent != nullptr);
  gra::LORENTZSCALAR coherent_lts = PhotonCovariancePoint("yy_uubargg");
  REQUIRE(coherent->Evaluate(coherent_lts, 0.118, false).Valid());
  const auto coherent_rows = coherent_lts.hard_color_flows.at(0).amplitudes.size();
  for (auto &flow : coherent_lts.hard_color_flows) {
    REQUIRE(gra::flux::ApplyEPAAmplitudeWeights(flow.amplitudes, common, gra::flux::EPASectorWeights{}));
    REQUIRE(flow.amplitudes.size() == coherent_rows);
    REQUIRE(gra::AllFinite(flow.amplitudes));
  }
}

TEST_CASE("Sampled nuclear EPA keeps generated color-flow source layouts",
          "[MG2GRA][gamma-gamma][EPA][color][nuclear][config]") {
  gra::PROC_003_QED_YY_EPA process;
  gra::LORENTZSCALAR       born = NuclearColoredPoint(gra::nuclear::CoherenceType::Inclusive, 0.07);
  REQUIRE(born.upc_model != nullptr);
  REQUIRE(born.upc_model->HasSamples());
  process.BindModelTune(born.model_cache->TunePtr());
  process.InitializeWorkerAmplitude(born);
  const double born_amp2 = process.Amp2(born);
  REQUIRE(born_amp2 > 0.0);
  REQUIRE(std::isfinite(born_amp2));
  REQUIRE(born.hamp.layout.epa_sector_resolved);
  REQUIRE(born.hamp.layout.epa_sector_count == std::array<std::uint8_t, 2>{2, 2});
  const auto born_flows = born.hard_color_flows;
  REQUIRE_FALSE(born_flows.empty());
  for (const auto &flow : born_flows) {
    REQUIRE(flow.amplitudes.size() == born.hamp.size());
    REQUIRE(gra::AllFinite(flow.amplitudes));
  }

  gra::LORENTZSCALAR shifted = NuclearColoredPoint(gra::nuclear::CoherenceType::Inclusive, 0.08);
  shifted.upc_model          = born.upc_model;
  const double shifted_amp2  = process.Amp2(shifted);
  REQUIRE(shifted_amp2 > 0.0);
  REQUIRE(std::isfinite(shifted_amp2));
  REQUIRE(shifted.hamp.size() == born.hamp.size());
  REQUIRE(shifted.hamp.layout.epa_sector_count == born.hamp.layout.epa_sector_count);
  REQUIRE(shifted.hamp.layout.epa_sector_type == born.hamp.layout.epa_sector_type);
  const auto shifted_flows = shifted.hard_color_flows;
  REQUIRE(shifted_flows.size() == born_flows.size());
  for (const auto &flow : shifted_flows) {
    REQUIRE(flow.amplitudes.size() == shifted.hamp.size());
    REQUIRE(gra::AllFinite(flow.amplitudes));
  }

  gra::PROC_003_QED_YY_EPA coherent_process;
  gra::LORENTZSCALAR       coherent      = NuclearColoredPoint(gra::nuclear::CoherenceType::Coherent, 0.07);
  coherent_process.BindModelTune(coherent.model_cache->TunePtr());
  coherent_process.InitializeWorkerAmplitude(coherent);
  const double             coherent_amp2 = coherent_process.Amp2(coherent);
  REQUIRE(coherent_amp2 > 0.0);
  // Keep both source sectors until the final nuclear projection
  REQUIRE(coherent.hamp.layout.epa_sector_count == std::array<std::uint8_t, 2>{2, 2});
  REQUIRE(born.hamp.size() == coherent.hamp.size());
  CHECK(coherent_amp2 == Approx(born_amp2).epsilon(2.0e-12));
  for (const auto &flow : coherent.hard_color_flows) {
    REQUIRE(flow.amplitudes.size() == coherent.hamp.size());
    REQUIRE(gra::AllFinite(flow.amplitudes));
  }
  gra::nuclear::ScreenLayout layout;
  layout.type = gra::nuclear::ScreenType::Fusion;
  layout.fusion.rows = coherent.hamp.layout.epa_rows_per_sector;
  for (auto& sector : layout.fusion.sector) {
    sector = {gra::nuclear::CoherenceType::Coherent, gra::nuclear::CoherenceType::Incoherent};
  }
  gra::nuclear::ScreenPoint point;
  point.amplitude = coherent.hamp;
  point.transfer = {gra::nuclear::RestTransfer(coherent.pbeam1, coherent.q1, coherent.beam1.mass),
                     gra::nuclear::RestTransfer(coherent.pbeam2, coherent.q2, coherent.beam2.mass)};
  for (const auto& flow : coherent.hard_color_flows) { point.flow.push_back(flow.amplitudes); }
  const auto selected = gra::nuclear::MUPCScreen(*coherent.upc_model, layout, point).Result();
  const auto all = gra::nuclear::MUPCScreen(*born.upc_model, layout, point).Result();
  CHECK(gra::Sum(all.helicity_norm) == Approx(gra::SquaredNorm(point.amplitude)).epsilon(2.0e-12));
  CHECK(gra::Sum(selected.helicity_norm) < gra::Sum(all.helicity_norm));
  for (const auto& flow : indices(point.flow)) {
    CHECK(all.color_norm[flow] == Approx(gra::SquaredNorm(point.flow[flow])).epsilon(2.0e-12));
    CHECK(selected.color_norm[flow] == Approx(selected.color_sector[flow][0]).epsilon(2.0e-12));
  }
}

TEST_CASE("Generated photon amplitudes retain coherent EPA covariance", "[MG2GRA][gamma-gamma][EPA][covariance]") {
  for (const auto &definition : gra::amplitude::Processes("PHOTON")) {
    CAPTURE(definition.process_name);
    REQUIRE(gra::HasPhotonMG5Process(definition));
    auto amplitude = gra::CreatePhotonMG5Process(definition);
    REQUIRE(amplitude != nullptr);
    const auto evaluate = [&](gra::LORENTZSCALAR &lts, bool coherent) {
      const auto result = amplitude->Evaluate(lts, 0.118, coherent);
      REQUIRE(result.Valid());
      REQUIRE(gra::AllFinite(lts.hamp));
      return result.amp2;
    };
    RequirePhotonCovariance(PhotonCovariancePoint(definition.process_name), evaluate, true, definition.process_name);
  }

  gra::AMP_MG5_yy_zjj amplitude;
  const auto          evaluate = [&](gra::LORENTZSCALAR &lts, bool coherent) {
    const auto result = amplitude.Evaluate(lts, 0.118, coherent);
    CAPTURE(static_cast<int>(result.status), lts.epa_hard.frame.valid);
    REQUIRE(result.Valid());
    REQUIRE(std::isfinite(result.amp2));
    return result.amp2;
  };
  RequirePhotonCovariance(PhotonCovariancePoint("yy_zjj"), evaluate, true, "yy_zjj");
  RequirePhotonFamilyCovariance<gra::AMP_MG5_yy_jj>("yy_jj", true);
  RequirePhotonFamilyCovariance<gra::AMP_MG5_yy_ww>("yy_ww", true);
}

TEST_CASE("Photon QCD amplitudes use alpha_s from the configured PDF", "[MG2GRA][gamma-gamma][alpha-s]") {
  ModelParamRestoreGuard restore;
  gra::MODELPARAM = "TUNE0";

  const auto quark_gluon = gra::amplitude::FindProcess("PHOTON", "yy_uubarg");
  const auto four_parton = gra::amplitude::FindProcess("PHOTON", "yy_uubargg");
  const auto dilepton    = gra::amplitude::FindProcess("PHOTON", "yy_ll");
  REQUIRE(quark_gluon.has_value());
  REQUIRE(four_parton.has_value());
  REQUIRE(dilepton.has_value());

  auto one_parton_emission = gra::CreatePhotonMG5Process(*quark_gluon);
  auto two_parton_emission = gra::CreatePhotonMG5Process(*four_parton);
  auto pure_qed            = gra::CreatePhotonMG5Process(*dilepton);
  REQUIRE(one_parton_emission != nullptr);
  REQUIRE(two_parton_emission != nullptr);
  REQUIRE(pure_qed != nullptr);
  REQUIRE(one_parton_emission->AlphaSPower() == 1);
  REQUIRE(two_parton_emission->AlphaSPower() == 2);
  REQUIRE(pure_qed->AlphaSPower() == 0);
  REQUIRE(one_parton_emission->AlphaQEDPower() == 2);
  REQUIRE(two_parton_emission->AlphaQEDPower() == 2);
  REQUIRE(pure_qed->AlphaQEDPower() == 2);

  gra::LORENTZSCALAR qcd = PhotonCovariancePoint("yy_uubarg");
  qcd.PDG                = LoadedPDGTable();
  qcd.alphaQCD           = 0.0;
  qcd.LHAPDFSET          = "NNPDF31_lo_as_0118";
  gra::MLHAPDFStore pdf_store;
  qcd.GlobalPdfPtr = pdf_store.GetPDF(qcd.LHAPDFSET, 0);
  gra::PROC_030_QED_YY_DZ_EPA process;
  process.BindModelTune(gra::MModelTune::Load(modelfile));
  process.InitGamma(qcd);
  const double amp2 = process.GammaGammaCON(qcd, false);
  const double Q2   = qcd.s_hat / 4.0;
  REQUIRE(amp2 > 0.0);
  REQUIRE(qcd.alphaQCD == Approx(qcd.GlobalPdfPtr->alphasQ2(Q2)).epsilon(1e-14));
  REQUIRE(qcd.muF == Approx(std::sqrt(Q2)));
  REQUIRE(qcd.muR == Approx(qcd.muF));
  REQUIRE(qcd.scalup == Approx(qcd.muF));

  gra::LORENTZSCALAR qed = PhotonCovariancePoint("yy_ll");
  qed.alphaQCD           = 0.0;
  const auto qed_result  = pure_qed->Evaluate(qed, 0.0, false);
  REQUIRE(qed_result.Valid());
  REQUIRE(qed_result.amp2 > 0.0);
  REQUIRE(gra::math::IsZero(qed.alphaQCD));

  gra::AMP_MG5_yy_jj  jj;
  gra::AMP_MG5_yy_ww  ww;
  gra::AMP_MG5_yy_zjj zjj;
  REQUIRE(jj.AlphaSPower() == 0);
  REQUIRE(ww.AlphaSPower() == 0);
  REQUIRE(zjj.AlphaSPower() == 0);
  REQUIRE(jj.AlphaQEDPower() == 2);
  REQUIRE(ww.AlphaQEDPower() == 4);
  REQUIRE(zjj.AlphaQEDPower() == 4);
}

// Check the covariant source projection before any hard helicity contraction
TEST_CASE("EPA hard sources retain independent beam-transverse directions",
          "[MG2GRA][gamma-gamma][EPA][screening][source]") {
  gra::LORENTZSCALAR lts;
  lts.q1 = gra::M4Vec(0.31, 0.12, 300.0, 300.0);
  lts.q2 = gra::M4Vec(0.07, -0.23, -220.0, 220.0);
  CompleteElasticForwardState(lts, 6500.0);

  const std::array<gra::M4Vec, 2> beam = {lts.pbeam1, lts.pbeam2};
  std::array<gra::M4Vec, 2>       transverse;
  REQUIRE(gra::mg5helas::BeamTransverseSource(lts.q1, beam, transverse[0]));
  REQUIRE(gra::mg5helas::BeamTransverseSource(lts.q2, beam, transverse[1]));
  for (const auto &source : transverse) {
    const double scale = std::max({1.0, source.P3mod() * lts.pbeam1.E(), source.P3mod() * lts.pbeam2.E()});
    REQUIRE(std::abs(source * lts.pbeam1) < 2.0e-13 * scale);
    REQUIRE(std::abs(source * lts.pbeam2) < 2.0e-13 * scale);
    REQUIRE(source.M2() < 0.0);
  }
  REQUIRE(transverse[0].Px() == Approx(lts.q1.Px()).margin(1.0e-13));
  REQUIRE(transverse[0].Py() == Approx(lts.q1.Py()).margin(1.0e-13));
  REQUIRE(transverse[1].Px() == Approx(lts.q2.Px()).margin(1.0e-13));
  REQUIRE(transverse[1].Py() == Approx(lts.q2.Py()).margin(1.0e-13));

  const gra::M4Vec hard             = lts.q1 + lts.q2;
  lts.s_hat                         = hard.M2();
  const auto                  rest  = PhotonPair(hard.M(), 0.0);
  std::vector<gra::M4Vec>     final = {BoostFromRestFrame(rest[0], hard), BoostFromRestFrame(rest[1], hard)};
  gra::mg5helas::EPAHardFrame frame;
  REQUIRE(gra::mg5helas::PrepareEPAHardFrame(lts, final, frame));
  std::array<gra::M4Vec, 2> source;
  REQUIRE(gra::mg5helas::TransformEPAHardSources(frame, lts.q1, lts.q2, source));
  const auto cosine = [](const gra::M4Vec &first, const gra::M4Vec &second) {
    REQUIRE(first.M2() < 0.0);
    REQUIRE(second.M2() < 0.0);
    return -(first * second) / std::sqrt(first.M2() * second.M2());
  };
  const double lab_cosine = cosine(transverse[0], transverse[1]);
  REQUIRE(std::abs(lab_cosine) < 0.5);
  REQUIRE(cosine(source[0], source[1]) == Approx(lab_cosine).margin(2.0e-12));

  constexpr double angle   = 0.83;
  const auto       rotated = TransformPhotonPoint(lts, [angle](gra::M4Vec &momentum) { momentum.RotateZ(angle); });
  const std::array<gra::M4Vec, 2> rotated_beam = {rotated.pbeam1, rotated.pbeam2};
  for (const auto &leg : indices(transverse)) {
    gra::M4Vec        rotated_source;
    const gra::M4Vec &rotated_q = leg == 0 ? rotated.q1 : rotated.q2;
    REQUIRE(gra::mg5helas::BeamTransverseSource(rotated_q, rotated_beam, rotated_source));
    gra::M4Vec expected = transverse[leg];
    expected.RotateZ(angle);
    const auto rotated_vector  = rotated_source.Contravariant<double>();
    const auto expected_vector = expected.Contravariant<double>();
    for (const auto &component : indices(expected_vector)) {
      REQUIRE(rotated_vector[component] == Approx(expected_vector[component]).margin(1.0e-12));
    }
  }
}

TEST_CASE("Prepared dilepton EPA hard tensors reproduce shifted currents",
          "[MG2GRA][gamma-gamma][EPA][screening][yy_ll]") {
  ModelParamRestoreGuard restore;
  gra::MODELPARAM       = "TUNE0";
  const auto definition = gra::amplitude::FindProcess("PHOTON", "yy_ll");
  REQUIRE(definition.has_value());
  auto prepared = gra::CreatePhotonMG5Process(*definition);
  REQUIRE(prepared != nullptr);

  gra::LORENTZSCALAR born = PhotonCovariancePoint("yy_ll");
  REQUIRE(prepared->Evaluate(born, 0.118, true).Valid());
  REQUIRE(born.epa_hard.Ready());

  gra::LORENTZSCALAR shifted = born;
  REQUIRE(ShiftEPAHardSources(shifted, 0.027, -0.019));
  const auto fast = gra::mg5helas::ContractEPAHard(shifted, born.epa_hard);
  REQUIRE(fast.Valid());

  auto fresh = gra::CreatePhotonMG5Process(*definition);
  REQUIRE(fresh != nullptr);
  gra::LORENTZSCALAR direct_lts = shifted;
  const auto         exact      = fresh->Evaluate(direct_lts, 0.118, true);
  REQUIRE(exact.Valid());
  REQUIRE(fast.amp2 == Approx(exact.amp2).epsilon(2.0e-11));
  REQUIRE(shifted.hamp.size() == direct_lts.hamp.size());
  for (const auto &i : indices(shifted.hamp)) { RequireComplexNear(shifted.hamp[i], direct_lts.hamp[i], 2.0e-10); }

  auto incoherent = gra::CreatePhotonMG5Process(*definition);
  REQUIRE(incoherent != nullptr);
  gra::LORENTZSCALAR incoherent_lts = born;
  REQUIRE(incoherent->Evaluate(incoherent_lts, 0.118, false).Valid());
  REQUIRE_FALSE(incoherent_lts.epa_hard.Ready());

  auto muon_kernel = gra::CreatePhotonMG5Process(*definition);
  REQUIRE(muon_kernel != nullptr);
  gra::LORENTZSCALAR muon = PhotonCovariancePoint("yy_ll");
  const double muon_mass = LoadedPDGTable().FindByPDG(13).mass;
  SLHAReader muon_card(gra::aux::ResolveProjectPath("MG5cards/Photon/yy_ll/param_card.dat"));
  muon_card.set_block_entry("mass", 11, muon_mass);
  muon_card.set_block_entry("yukawa", 11, muon_mass);
  muon_kernel->InitParameters(muon_card);
  const auto muon_pair = PhotonPair(muon.pfinal[0].M(), muon_mass);
  for (const auto &i : indices(muon.decaytree)) {
    muon.decaytree[i].p.mass = muon_mass;
    muon.decaytree[i].p4 = muon_pair[i];
  }
  REQUIRE(muon_kernel->Evaluate(muon, 0.0, true).Valid());
  REQUIRE(muon.epa_hard.Ready());
  REQUIRE(ShiftEPAHardSources(muon, -0.013, 0.021));
  const auto muon_fast = gra::mg5helas::ContractEPAHard(muon, muon.epa_hard);
  REQUIRE(muon_fast.Valid());

  auto muon_direct = gra::CreatePhotonMG5Process(*definition);
  REQUIRE(muon_direct != nullptr);
  muon_direct->InitParameters(muon_card);
  gra::LORENTZSCALAR muon_check = muon;
  const auto         muon_exact = muon_direct->Evaluate(muon_check, 0.0, true);
  REQUIRE(muon_exact.Valid());
  REQUIRE(muon_fast.amp2 == Approx(muon_exact.amp2).epsilon(2.0e-11));
  REQUIRE(muon.hamp.size() == muon_check.hamp.size());
  for (const auto &i : indices(muon.hamp)) { RequireComplexNear(muon.hamp[i], muon_check.hamp[i], 2.0e-10); }

  gra::LORENTZSCALAR routed    = PhotonCovariancePoint("yy_ll");
  routed.PDG                   = LoadedPDGTable();
  routed.decaytree[0].p        = routed.PDG.FindByPDG(-13);
  routed.decaytree[1].p        = routed.PDG.FindByPDG(13);
  const auto routed_pair       = PhotonPair(routed.pfinal[0].M(), routed.decaytree[0].p.mass);
  routed.decaytree[0].p4       = routed_pair[0];
  routed.decaytree[1].p4       = routed_pair[1];
  const auto  gamma_definition = gra::MGamma::ProcessDefinitionFor(gra::MGammaMode::FermionPair);
  const auto  model_tune       = gra::MModelTune::Load(modelfile);
  gra::MGamma routed_gamma(routed, model_tune, gamma_definition);
  REQUIRE(routed_gamma.yyffbar(routed, true).Valid());
  REQUIRE(routed.epa_hard.Ready());
  REQUIRE(ShiftEPAHardSources(routed, 0.009, -0.016));
  const auto routed_fast = gra::mg5helas::ContractEPAHard(routed, routed.epa_hard);
  REQUIRE(routed_fast.Valid());

  gra::LORENTZSCALAR routed_check = routed;
  gra::MGamma        routed_direct(routed_check, model_tune, gamma_definition);
  const auto         routed_exact = routed_direct.yyffbar(routed_check, true);
  REQUIRE(routed_exact.Valid());
  REQUIRE(routed_fast.amp2 == Approx(routed_exact.amp2).epsilon(2.0e-11));
}

// Check the EPA limit against exact QED and the collider symmetries
TEST_CASE("Proton EPA dimuons approach exact QED at small photon virtuality",
          "[gamma-gamma][EPA][QED][closure]") {
  for (const double scale : {0.1, 1.0}) {
    for (const double phi : {0.3, 1.2, 2.4}) {
      auto lts = DirectCentralPairLTSForTest(-13, 13);
      lts.q1 = gra::M4Vec(0.1, 0.0, 30.0, 30.0) * scale;
      lts.q2 = gra::M4Vec(0.2 * std::cos(phi), 0.2 * std::sin(phi), -40.0, 40.0) * scale;
      CompleteElasticForwardState(lts, 6500.0);
      lts.s = (lts.pbeam1 + lts.pbeam2).M2();
      lts.s_hat = lts.pfinal[0].M2();
      lts.m2 = lts.s_hat;
      lts.process.FORWARD_NOFLIP = false;
      const auto pair = TwoBodyRestKinematics(lts.pfinal[0].M(), lts.decaytree[0].p.mass,
                                              lts.decaytree[1].p.mass, 0.91, -0.38);
      for (const auto &i : indices(lts.decaytree)) {
        lts.decaytree[i].p4 = BoostFromRestFrame(pair[i], lts.pfinal[0]);
      }
      std::vector<gra::LORENTZSCALAR> points = {lts};
      points.push_back(TransformPhotonPoint(lts, [](gra::M4Vec &p) { p.RotateZ(0.83); }));
      auto exchanged = TransformPhotonPoint(lts, [](gra::M4Vec &p) { p.Flip3(); });
      std::swap(exchanged.pbeam1, exchanged.pbeam2);
      std::swap(exchanged.pfinal[1], exchanged.pfinal[2]);
      std::swap(exchanged.q1, exchanged.q2);
      std::swap(exchanged.t1, exchanged.t2);
      std::swap(exchanged.qt1, exchanged.qt2);
      std::swap(exchanged.x1, exchanged.x2);
      std::swap(exchanged.xi1, exchanged.xi2);
      points.push_back(exchanged);
      auto conjugated = lts;
      std::swap(conjugated.decaytree[0].p4, conjugated.decaytree[1].p4);
      points.push_back(conjugated);

      std::array<double, 2> reference{};
      for (const auto &i : indices(points)) {
        auto &point = points[i];
        gra::MTensorPomeron tensor(point, gra::MModelTune::Load(modelfile),
                                   gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::QED));
        const double exact = tensor.ME4(point, gra::TensorContinuumMode::QED);
        gra::MGamma gamma(point, gra::MModelTune::Load(modelfile),
                           gra::MGamma::ProcessDefinitionFor(gra::MGammaMode::FermionPair));
        const auto result = gamma.yyffbar(point, true);
        REQUIRE(result.Valid());
        REQUIRE(point.epa_hard.Ready());
        REQUIRE(gra::mg5helas::UsesProtonEPAHardSources(point));
        const double epa = gra::flux::ApplyktEPAfluxes(result.amp2, point).amp2;
        CAPTURE(scale, phi, i, exact, epa);
        CHECK(epa == Approx(exact).epsilon(0.01));
        if (i == 0) { reference = {exact, epa}; }
        CHECK(exact == Approx(reference[0]).epsilon(3.0e-4));
        CHECK(epa == Approx(reference[1]).epsilon(3.0e-4));
      }
    }
  }
}

// Check generated families reuse only the hard tensor from the same Born state
TEST_CASE("Generated photon families cache shifted EPA hard tensors", "[MG2GRA][gamma-gamma][EPA][screening][family]") {
  gra::LORENTZSCALAR born = PhotonCovariancePoint("yy_jj");
  gra::AMP_MG5_yy_jj prepared;
  const auto         born_result = prepared.Evaluate(born, 0.118, true);
  REQUIRE(born_result.Valid());
  REQUIRE(born.epa_hard.Ready());

  gra::LORENTZSCALAR shifted = born;
  REQUIRE(ShiftEPAHardSources(shifted, 0.017, -0.011));
  shifted.exact_forward_photon_kinematics = true;
  const auto fast                         = gra::mg5helas::ContractEPAHard(shifted, born.epa_hard);
  REQUIRE(fast.Valid());

  gra::LORENTZSCALAR direct_lts = shifted;
  gra::AMP_MG5_yy_jj direct;
  const auto         exact = direct.Evaluate(direct_lts, 0.118, true);
  REQUIRE(exact.Valid());
  REQUIRE(fast.amp2 == Approx(exact.amp2).epsilon(2.0e-11));
  REQUIRE(shifted.hamp.size() == direct_lts.hamp.size());
  for (const auto &i : indices(shifted.hamp)) { RequireComplexNear(shifted.hamp[i], direct_lts.hamp[i], 2.0e-10); }
  const auto shifted_flows = shifted.hard_color_flows;
  const auto direct_flows  = direct_lts.hard_color_flows;
  REQUIRE(shifted_flows.size() == direct_flows.size());
  for (const auto &flow : indices(shifted_flows)) {
    REQUIRE(shifted_flows[flow].amplitudes.size() == direct_flows[flow].amplitudes.size());
    for (const auto &row : indices(shifted_flows[flow].amplitudes)) {
      RequireComplexNear(shifted_flows[flow].amplitudes[row], direct_flows[flow].amplitudes[row], 2.0e-10);
    }
  }

  REQUIRE(shifted.epa_hard.Ready());
  REQUIRE(gra::mg5helas::ContractEPAHard(shifted, shifted.epa_hard).Valid());

  gra::AMP_MG5_yy_jj incoherent;
  REQUIRE(incoherent.Evaluate(born, 0.118, false).Valid());
  REQUIRE_FALSE(born.epa_hard.Ready());
}

// Check that every prepared shower flow follows the shifted EPA source and normalization
TEST_CASE("Prepared colored EPA hard tensors refresh shifted shower flows",
          "[MG2GRA][gamma-gamma][EPA][screening][color]") {
  const auto definition = gra::amplitude::FindProcess("PHOTON", "yy_uubargg");
  REQUIRE(definition.has_value());

  gra::PROC_003_QED_YY_EPA prepared;
  gra::LORENTZSCALAR       born      = PhotonCovariancePoint("yy_uubargg");
  prepared.BindModelTune(born.model_cache->TunePtr());
  prepared.InitializeWorkerAmplitude(born);
  const double             born_amp2 = prepared.Amp2(born);
  REQUIRE(born_amp2 > 0.0);
  REQUIRE(born.epa_hard.Ready());
  const auto born_flows = born.hard_color_flows;
  REQUIRE(born_flows.size() == 2);

  gra::LORENTZSCALAR shifted = born;
  REQUIRE(ShiftEPAHardSources(shifted, 0.027, -0.019));
  gra::LORENTZSCALAR direct_lts = shifted;

  const double shifted_amp2 = prepared.EPAHardAmp2(shifted);
  REQUIRE(shifted_amp2 > 0.0);
  const auto shifted_flows = shifted.hard_color_flows;
  REQUIRE(shifted_flows.size() == born_flows.size());

  gra::PROC_003_QED_YY_EPA direct;
  direct.BindModelTune(direct_lts.model_cache->TunePtr());
  direct.InitializeWorkerAmplitude(direct_lts);
  double                   direct_amp2 = direct.Amp2(direct_lts);
  REQUIRE(direct_amp2 > 0.0);
  REQUIRE(shifted_amp2 == Approx(direct_amp2).epsilon(2.0e-11));
  REQUIRE(shifted.hamp.size() == direct_lts.hamp.size());
  for (const auto &i : indices(shifted.hamp)) { RequireComplexNear(shifted.hamp[i], direct_lts.hamp[i], 2.0e-10); }

  REQUIRE(gra::mg5helas::ContractEPAHard(shifted, shifted.epa_hard).Valid());

  const auto direct_flows = direct_lts.hard_color_flows;
  REQUIRE(direct_flows.size() == shifted_flows.size());
  double flow_change = 0.0;
  double flow_scale  = 0.0;
  for (const auto &flow : indices(shifted_flows)) {
    REQUIRE(shifted_flows[flow].amplitudes.size() == direct_flows[flow].amplitudes.size());
    REQUIRE(shifted_flows[flow].amplitudes.size() == born_flows[flow].amplitudes.size());
    for (const auto &row : indices(shifted_flows[flow].amplitudes)) {
      RequireComplexNear(shifted_flows[flow].amplitudes[row], direct_flows[flow].amplitudes[row], 2.0e-10);
      flow_change += std::norm(shifted_flows[flow].amplitudes[row] - born_flows[flow].amplitudes[row]);
      flow_scale += std::norm(shifted_flows[flow].amplitudes[row]) + std::norm(born_flows[flow].amplitudes[row]);
    }
  }
  REQUIRE(flow_scale > 0.0);
  REQUIRE(flow_change > 1.0e-12 * flow_scale);
}

// Check that the analytic loop keeps the hard subprocess fixed in the screening convolution
TEST_CASE("Prepared light-by-light EPA hard tensors reproduce shifted currents",
          "[gamma-gamma][light-by-light][EPA][screening]") {
  gra::LORENTZSCALAR born = PhotonCovariancePoint("yy_yy");
  born.beam1.mass         = gra::PDG::mp;
  born.beam1.spinX2       = 1;
  born.beam2              = born.beam1;
  REQUIRE(gra::mg5helas::UsesProtonEPAHardSources(born));
  gra::AMP_yy_yy prepared(born.model_cache->Tune().SM());
  REQUIRE(prepared.Evaluate(born, 0.0, true).Valid());
  REQUIRE(born.epa_hard.Ready());
  REQUIRE(born.hamp.metadata.spin_basis == gra::ScreeningSpinBasis::ProtonHelicity);
  REQUIRE(born.hamp.metadata.spin_rows == 16);

  gra::LORENTZSCALAR shifted = born;
  REQUIRE(gra::kinematics::BuildEPATransfer(born.pbeam1, born.xi1, born.q1.Px() + 0.027, born.q1.Py() - 0.019,
                                            pow2(born.beam1.mass), born.forward_mass2[0], true, shifted.pfinal[1],
                                            shifted.q1, shifted.t1));
  REQUIRE(gra::kinematics::BuildEPATransfer(born.pbeam2, born.xi2, born.q2.Px() - 0.027, born.q2.Py() + 0.019,
                                            pow2(born.beam2.mass), born.forward_mass2[1], false, shifted.pfinal[2],
                                            shifted.q2, shifted.t2));
  shifted.qt1     = shifted.q1.Pt();
  shifted.qt2     = shifted.q2.Pt();
  const auto fast = gra::mg5helas::ContractEPAHard(shifted, born.epa_hard);
  REQUIRE(fast.Valid());

  gra::LORENTZSCALAR direct_lts = shifted;
  gra::AMP_yy_yy direct(born.model_cache->Tune().SM());
  const auto         exact = direct.Evaluate(direct_lts, 0.0, true);
  REQUIRE(exact.Valid());
  REQUIRE(fast.amp2 == Approx(exact.amp2).epsilon(2.0e-11));
  REQUIRE(shifted.hamp.size() == direct_lts.hamp.size());
  for (const auto &i : indices(shifted.hamp)) { RequireComplexNear(shifted.hamp[i], direct_lts.hamp[i], 2.0e-10); }

  gra::AMP_yy_yy incoherent(born.model_cache->Tune().SM());
  gra::LORENTZSCALAR incoherent_lts = born;
  REQUIRE(incoherent.Evaluate(incoherent_lts, 0.0, false).Valid());
  REQUIRE_FALSE(incoherent_lts.epa_hard.Ready());
}

// Check the factorized nuclear source path independently of exact proton currents
TEST_CASE("Prepared nuclear light-by-light EPA tensors reproduce shifted sources",
          "[gamma-gamma][light-by-light][EPA][screening][nuclear]") {
  gra::LORENTZSCALAR born  = NuclearColoredPoint(gra::nuclear::CoherenceType::Coherent, 0.07);
  const auto         final = PhotonPair(born.pfinal[0].M(), 0.0);
  born.decaytree           = {PhotonLeaf(22, final[0]), PhotonLeaf(22, final[1])};
  REQUIRE_FALSE(gra::mg5helas::UsesProtonEPAHardSources(born));

  gra::AMP_yy_yy prepared(born.model_cache->Tune().SM());
  REQUIRE(prepared.Evaluate(born, 0.0, true).Valid());
  REQUIRE(born.epa_hard.Ready());

  gra::LORENTZSCALAR shifted = born;
  REQUIRE(gra::kinematics::BuildEPATransfer(born.pbeam1, born.xi1, born.q1.Px() + 0.017, born.q1.Py() - 0.011,
                                            pow2(born.beam1.mass), born.forward_mass2[0], true, shifted.pfinal[1],
                                            shifted.q1, shifted.t1));
  REQUIRE(gra::kinematics::BuildEPATransfer(born.pbeam2, born.xi2, born.q2.Px() - 0.017, born.q2.Py() + 0.011,
                                            pow2(born.beam2.mass), born.forward_mass2[1], false, shifted.pfinal[2],
                                            shifted.q2, shifted.t2));
  shifted.qt1     = shifted.q1.Pt();
  shifted.qt2     = shifted.q2.Pt();
  const auto fast = gra::mg5helas::ContractEPAHard(shifted, born.epa_hard);
  REQUIRE(fast.Valid());

  gra::LORENTZSCALAR direct_lts = shifted;
  gra::AMP_yy_yy direct(born.model_cache->Tune().SM());
  const auto         exact = direct.Evaluate(direct_lts, 0.0, true);
  REQUIRE(exact.Valid());
  REQUIRE(fast.amp2 == Approx(exact.amp2).epsilon(2.0e-11));
  REQUIRE(shifted.hamp.size() == direct_lts.hamp.size());
  for (const auto &i : indices(shifted.hamp)) { RequireComplexNear(shifted.hamp[i], direct_lts.hamp[i], 2.0e-10); }
}

// Check process selection rejects a mismatched Durham decay topology
TEST_CASE("Generated Durham process interfaces reject topology mismatch", "[MG2GRA][registry][decay]") {
  auto stable_branch = [](int pdg) {
    gra::MDecayBranch branch;
    branch.p.pdg = pdg;
    return branch;
  };

  gra::LORENTZSCALAR wrong_gg;
  wrong_gg.decaytree = {stable_branch(2), stable_branch(-2)};
  const auto process = gra::amplitude::FindProcess("DURHAM", "gg_gg");
  REQUIRE(process.has_value());
  auto gluons = gra::CreateDurhamMG5Process(*process);
  REQUIRE(gluons != nullptr);
  CHECK_FALSE(gluons->MatchProcess(wrong_gg.decaytree).has_value());
}

TEST_CASE("Mapped photon quark pairs use the generic color-flow assignment", "[MG2GRA][gamma-gamma][kernel][color]") {
  ModelParamRestoreGuard restore;
  gra::MODELPARAM = "TUNE0";

  // Build one physical gamma-gamma to quark-antiquark point
  auto quark_pair = []() {
    gra::LORENTZSCALAR lts;
    lts.PDG.ReadParticleData(pdgfile);
    const double energy                 = 50.0;
    const double mass                   = lts.PDG.FindByPDG(2).mass;
    const double momentum               = std::sqrt(energy * energy - mass * mass);
    lts.q1                              = gra::M4Vec(0.0, 0.0, energy, energy);
    lts.q2                              = gra::M4Vec(0.0, 0.0, -energy, energy);
    lts.s_hat                           = 4.0 * energy * energy;

    gra::MDecayBranch quark;
    quark.p  = lts.PDG.FindByPDG(2);
    quark.p4 = gra::M4Vec(momentum, 0.0, 0.0, energy);
    gra::MDecayBranch antiquark;
    antiquark.p   = lts.PDG.FindByPDG(-2);
    antiquark.p4  = gra::M4Vec(-momentum, 0.0, 0.0, energy);
    lts.decaytree = {quark, antiquark};
    return lts;
  };

  gra::LORENTZSCALAR          lts        = quark_pair();
  const auto                  model_tune = gra::MModelTune::Load(modelfile);
  gra::PROC_030_QED_YY_DZ_EPA process;
  process.BindModelTune(model_tune);
  process.InitGamma(lts);
  const double quark_amp2 = process.GammaGammaCON(lts, false);
  REQUIRE(quark_amp2 > 0.0);
  REQUIRE(gra::mg5helas::EvaluationSucceeded(process.EvaluationStatus()));
  REQUIRE(lts.decaytree[0].p.pdg == 2);
  REQUIRE(lts.decaytree[1].p.pdg == -2);
  double quark_hamp2 = 0.0;
  for (const auto &amplitude : lts.hamp) { quark_hamp2 += std::norm(amplitude); }
  REQUIRE(quark_hamp2 == Approx(4.0 * quark_amp2).epsilon(1e-12));

  // The generic Dirac kernel differs only by the quark Q^4 N_c factor
  gra::LORENTZSCALAR electron_lts = quark_pair();
  electron_lts.decaytree[0].p     = electron_lts.PDG.FindByPDG(-11);
  electron_lts.decaytree[1].p     = electron_lts.PDG.FindByPDG(11);
  for (auto &branch : electron_lts.decaytree) { branch.p.mass = lts.decaytree[0].p.mass; }
  gra::PROC_030_QED_YY_DZ_EPA electron_process;
  electron_process.BindModelTune(model_tune);
  electron_process.InitGamma(electron_lts);
  const double unit_charge_amp2 = electron_process.GammaGammaCON(electron_lts, false);
  REQUIRE(quark_amp2 == Approx((16.0 / 27.0) * unit_charge_amp2).epsilon(1e-12));
  REQUIRE(lts.decaytree[0].p.color_flow.empty());
  REQUIRE(lts.decaytree[1].p.color_flow.empty());

  process.SampleColorFlow(lts);
  REQUIRE(lts.decaytree[0].p.color_flow.flow1 == 501);
  REQUIRE(lts.decaytree[0].p.color_flow.flow2 == 0);
  REQUIRE(lts.decaytree[1].p.color_flow.flow1 == 0);
  REQUIRE(lts.decaytree[1].p.color_flow.flow2 == 501);

  // A failed mapped evaluation restores particle ids and preserves its status
  gra::LORENTZSCALAR invalid = quark_pair();
  const double       nan     = std::numeric_limits<double>::quiet_NaN();
  invalid.decaytree[0].p4    = gra::M4Vec(nan, 0.0, 0.0, nan);
  CHECK(gra::math::IsZero(process.GammaGammaCON(invalid, false)));
  CHECK(process.EvaluationStatus() == gra::mg5helas::EvaluationStatus::KinematicsFailure);
  CHECK(invalid.decaytree[0].p.pdg == 2);
  CHECK(invalid.decaytree[1].p.pdg == -2);
  CHECK(invalid.hamp.empty());

  // The shared machinery follows stable-leaf PDGs and clears stale parent tags
  gra::LORENTZSCALAR nested = quark_pair();
  gra::MDecayBranch  parent;
  parent.p.pdg        = 991;
  parent.p.color_flow = {700, 701};
  parent.legs         = {nested.decaytree[1], nested.decaytree[0]};
  nested.decaytree    = {parent};
  REQUIRE(gra::AssignPhotonSingletQuarkPairColorFlow(nested));
  REQUIRE(nested.decaytree[0].p.color_flow.empty());
  REQUIRE(nested.decaytree[0].legs[0].p.color_flow.flow1 == 0);
  REQUIRE(nested.decaytree[0].legs[0].p.color_flow.flow2 == 501);
  REQUIRE(nested.decaytree[0].legs[1].p.color_flow.flow1 == 501);
  REQUIRE(nested.decaytree[0].legs[1].p.color_flow.flow2 == 0);
}

// Check normalization and azimuthal phase against analytic half-angle spinors
TEST_CASE("Massless HELAS spinors retain normalization near the negative beam", "[MG2GRA][HELAS][helicity]") {
  constexpr double energy = 140.0;
  constexpr double tolerance = 16.0 * std::numeric_limits<double>::epsilon();
  for (const double theta : {1e-2, 1e-4, 1e-6, 1e-8, 1e-10}) {
    for (const double phi : {-2.1, 0.7, 2.4}) {
      std::array<double, 4> p = {energy, energy * std::sin(theta) * std::cos(phi),
                                energy * std::sin(theta) * std::sin(phi), -energy * std::cos(theta)};
      for (const int helicity : {-1, 1}) {
        for (const int sign : {-1, 1}) {
          std::array<std::complex<double>, 6> incoming{}, outgoing{};
          MG5_sm_lepton_masses::ixxxxx(p.data(), 0.0, helicity, sign, incoming.data());
          MG5_sm_lepton_masses::oxxxxx(p.data(), 0.0, helicity, sign, outgoing.data());
          for (const auto &spinor : {incoming, outgoing}) {
            const auto norm = gra::SquaredNorm(std::span<const std::complex<double>>(spinor).subspan(2));
            CAPTURE(theta, phi, helicity, sign, norm);
            REQUIRE(std::abs(norm / (2.0 * energy) - 1.0) < tolerance);
          }
          if (helicity == 1 && sign == 1) {
            const double small = std::sqrt(2.0 * energy) * std::sin(theta / 2.0);
            const auto large = std::polar(std::sqrt(2.0 * energy) * std::cos(theta / 2.0), phi);
            REQUIRE(std::abs(incoming[4] - small) < tolerance * small);
            REQUIRE(std::abs(incoming[5] - large) < tolerance * std::abs(large));
          }
        }
      }
    }
  }
}

// Check the collider helicity section as either beam becomes collinear
TEST_CASE("MG5 incoming HELAS states approach one collider beam section", "[MG2GRA][HELAS][helicity][covariance]") {
  constexpr double energy              = 100.0;
  constexpr double transverse_momentum = 1.0e-3;
  constexpr double collinear_tolerance = 1.0e-4;

  for (const int beam_sign : {-1, 1}) {
    for (const int helicity : {-1, 1}) {
      for (const double phi : {-2.17, -0.41, 0.73, 2.49}) {
        std::array<double, 4> exact_momentum  = {energy, 0.0, 0.0, beam_sign * energy};
        std::array<double, 4> nearby_momentum = {
            energy, transverse_momentum * std::cos(phi), transverse_momentum * std::sin(phi),
            beam_sign * std::sqrt(energy * energy - transverse_momentum * transverse_momentum)};

        std::array<std::complex<double>, 6> exact_vector  = {};
        std::array<std::complex<double>, 6> nearby_vector = {};
        MG5_sm_lepton_masses::vxxxxx(exact_momentum.data(), 0.0, helicity, -1, exact_vector.data());
        MG5_sm_lepton_masses::vxxxxx(nearby_momentum.data(), 0.0, helicity, -1, nearby_vector.data());
        const auto vector_transport =
            gra::mg5helas::IncomingTransport(nearby_momentum.data(), gra::PDG::PDG_gluon, helicity);
        for (std::size_t component = 2; component < exact_vector.size(); ++component) {
          CAPTURE(beam_sign, helicity, phi, component);
          RequireComplexNear(vector_transport * nearby_vector[component], exact_vector[component], collinear_tolerance);
        }

        std::array<std::complex<double>, 6> exact_antifermion  = {};
        std::array<std::complex<double>, 6> nearby_antifermion = {};
        MG5_sm_lepton_masses::ixxxxx(exact_momentum.data(), 0.0, helicity, 1, exact_antifermion.data());
        MG5_sm_lepton_masses::ixxxxx(nearby_momentum.data(), 0.0, helicity, 1, nearby_antifermion.data());
        const auto antifermion_transport = gra::mg5helas::IncomingTransport(nearby_momentum.data(), -2, helicity);
        for (std::size_t component = 2; component < exact_antifermion.size(); ++component) {
          CAPTURE(beam_sign, helicity, phi, component);
          RequireComplexNear(antifermion_transport * nearby_antifermion[component], exact_antifermion[component],
                             collinear_tolerance);
        }

        std::array<std::complex<double>, 6> exact_fermion  = {};
        std::array<std::complex<double>, 6> nearby_fermion = {};
        MG5_sm_lepton_masses::oxxxxx(exact_momentum.data(), 0.0, helicity, -1, exact_fermion.data());
        MG5_sm_lepton_masses::oxxxxx(nearby_momentum.data(), 0.0, helicity, -1, nearby_fermion.data());
        const auto fermion_transport = gra::mg5helas::IncomingTransport(nearby_momentum.data(), 2, helicity);
        for (std::size_t component = 2; component < exact_fermion.size(); ++component) {
          CAPTURE(beam_sign, helicity, phi, component);
          RequireComplexNear(fermion_transport * nearby_fermion[component], exact_fermion[component],
                             collinear_tolerance);
        }
      }
    }
  }

  SECTION("rounded longitudinal components retain the azimuthal helicity phase") {
    constexpr double      rounded_transverse_momentum = 1.0e-7;
    constexpr double      phi                         = 0.73;
    std::array<double, 4> exact_momentum              = {energy, 0.0, 0.0, -energy};
    std::array<double, 4> nearby_momentum             = {
                    energy, rounded_transverse_momentum * std::cos(phi), rounded_transverse_momentum * std::sin(phi),
                    -std::sqrt(energy * energy - rounded_transverse_momentum * rounded_transverse_momentum)};
    REQUIRE(nearby_momentum[0] + nearby_momentum[3] <= 0.0);

    for (const int helicity : {-1, 1}) {
      const auto transport = gra::mg5helas::IncomingFermionTransport(nearby_momentum.data(), helicity);
      RequireComplexNear(transport, -std::polar(1.0, -static_cast<double>(helicity) * phi), 1e-15);
      std::array<std::complex<double>, 6> exact  = {};
      std::array<std::complex<double>, 6> nearby = {};
      MG5_sm_lepton_masses::ixxxxx(exact_momentum.data(), 0.0, helicity, 1, exact.data());
      MG5_sm_lepton_masses::ixxxxx(nearby_momentum.data(), 0.0, helicity, 1, nearby.data());
      for (std::size_t component = 2; component < exact.size(); ++component) {
        CAPTURE(helicity, component);
        // The spinor approaches its beam value with a correction proportional to pT/sqrt(E)
        REQUIRE(std::abs(transport * nearby[component] - exact[component]) <=
                rounded_transverse_momentum / std::sqrt(energy));
      }
    }
  }
}



TEST_CASE("MG5 helicity contractions preserve complex bilinear ordering", "[MG2GRA][helicity][MMatrix]") {
  const std::array<std::complex<double>, 4> hard = {std::complex<double>(1.0, 2.0), std::complex<double>(-3.0, 4.0),
                                                    std::complex<double>(5.0, -6.0), std::complex<double>(-7.0, -8.0)};

  // One-hot sources select (--,-+,+-,++) with the lower helicity fastest
  for (std::size_t upper_index = 0; upper_index < 2; ++upper_index) {
    for (std::size_t lower_index = 0; lower_index < 2; ++lower_index) {
      std::array<std::complex<double>, 2> upper = {};
      std::array<std::complex<double>, 2> lower = {};
      upper[upper_index]                        = 1.0;
      lower[lower_index]                        = 1.0;
      RequireComplexNear(gra::mg5helas::ContractHelicitySources(hard, upper, lower),
                         hard[2 * upper_index + lower_index], 0.0);
    }
  }

  const std::array<std::complex<double>, 2> upper = {std::complex<double>(2.0, -1.0), std::complex<double>(-3.0, 4.0)};
  const std::array<std::complex<double>, 2> lower = {std::complex<double>(5.0, 2.0), std::complex<double>(-1.0, -3.0)};
  const std::complex<double> source_reference     = upper[0] * lower[0] * hard[0] + upper[0] * lower[1] * hard[1] +
                                                upper[1] * lower[0] * hard[2] + upper[1] * lower[1] * hard[3];
  RequireComplexNear(gra::mg5helas::ContractHelicitySources(hard, upper, lower), source_reference, 1.0e-14);

  const std::array<std::complex<double>, 4> kernel = {std::complex<double>(1.0, -2.0), std::complex<double>(3.0, 1.0),
                                                      std::complex<double>(-2.0, 5.0), std::complex<double>(4.0, -3.0)};
  const std::complex<double>                kernel_reference =
      kernel[0] * hard[0] + kernel[1] * hard[1] + kernel[2] * hard[2] + kernel[3] * hard[3];
  const std::complex<double> conjugating_reference = std::conj(kernel[0]) * hard[0] + std::conj(kernel[1]) * hard[1] +
                                                     std::conj(kernel[2]) * hard[2] + std::conj(kernel[3]) * hard[3];
  const auto contracted = gra::mg5helas::ContractHelicityKernel(hard, kernel);
  RequireComplexNear(contracted, kernel_reference, 1.0e-14);
  REQUIRE(std::abs(contracted - conjugating_reference) > 1.0);
}

TEST_CASE("MColorFlow defaults to colorless event tags", "[gra::MParticle][color]") {
  gra::MParticle particle;
  REQUIRE(particle.color_flow.empty());
  particle.color_flow.flow1 = 501;
  particle.color_flow.flow2 = 502;
  REQUIRE_FALSE(particle.color_flow.empty());
  particle.color_flow.clear();
  REQUIRE(particle.color_flow.empty());
}

TEST_CASE("gra::qed neutral-current helpers provide photon and Z currents", "[gra::qed][ew]") {
  gra::MDirac  dirac("DIRAC");
  const double sin2thetaW = 0.23122;
  const auto   electron   = gra::qed::FermionNeutralCurrentCouplings(11, sin2thetaW);
  const auto   up         = gra::qed::FermionNeutralCurrentCouplings(2, sin2thetaW);
  const auto   down       = gra::qed::FermionNeutralCurrentCouplings(1, sin2thetaW);

  REQUIRE(electron.chargeX1 == Approx(-1.0));
  REQUIRE(up.chargeX1 == Approx(2.0 / 3.0));
  REQUIRE(down.chargeX1 == Approx(-1.0 / 3.0));
  REQUIRE(electron.gL == Approx(-0.5 + sin2thetaW));
  REQUIRE(electron.gR == Approx(sin2thetaW));

  const auto   lts   = MakeToyPhotoZFFbar(11);
  const double e_qed = gra::qed::e_QED();
  const auto   j_gamma =
      gra::qed::PhotonFermionCurrent(dirac, lts.decaytree[0].p4, lts.decaytree[1].p4, 11, -1, 1, e_qed);
  const auto j_z =
      gra::qed::ZFermionCurrent(dirac, lts.decaytree[0].p4, lts.decaytree[1].p4, 11, -1, 1, e_qed, sin2thetaW);
  REQUIRE(j_gamma.size() == 4);
  REQUIRE(j_z.size() == 4);
  for (std::size_t i = 0; i < 4; ++i) {
    REQUIRE(std::isfinite(j_gamma[i].real()));
    REQUIRE(std::isfinite(j_gamma[i].imag()));
    REQUIRE(std::isfinite(j_z[i].real()));
    REQUIRE(std::isfinite(j_z[i].imag()));
  }

  const gra::M4Vec q = lts.decaytree[0].p4 + lts.decaytree[1].p4;
  for (const int hf : {-1, 1}) {
    for (const int ha : {-1, 1}) {
      const auto current =
          gra::qed::PhotonFermionCurrent(dirac, lts.decaytree[0].p4, lts.decaytree[1].p4, 11, hf, ha, e_qed);
      std::complex<double> ward  = 0.0;
      double               scale = 0.0;
      for (std::size_t mu = 0; mu < 4; ++mu) {
        ward += q[mu] * current[mu];
        scale += std::abs(q[mu]) * std::abs(current[mu]);
      }
      REQUIRE(std::abs(ward) <= 1e-10 * std::max(1.0, scale));
    }
  }
}

TEST_CASE("qed::alpha_QED obeys scale, threshold, and steering contracts", "[gra::qed][running-coupling]") {
  REQUIRE(qed::alpha_QED(0.0, "LL") == Approx(qed::alpha_0));
  REQUIRE(qed::alpha_QED(0.5 * math::pow2(qed::Q_ref), "LL") == Approx(qed::alpha_0));
  REQUIRE(qed::alpha_QED(1e8, "ZERO") == Approx(qed::alpha_0));
  constexpr double alpha_mg = 1.0 / 127.9;
  REQUIRE(qed::alpha_QED(1e4, "ZERO", alpha_mg) == Approx(qed::alpha_0));
  REQUIRE(qed::alpha_QED(1e4, "MG", alpha_mg) == Approx(alpha_mg));
  REQUIRE(qed::alpha_QED(1e4, "LL", alpha_mg) == Approx(qed::alpha_QED(1e4, "LL")));

  const std::vector<double> scales   = {qed::Q_ref, 0.01, 0.1, 0.2, 1.0, 2.0, 10.0, 91.2, 1000.0};
  double                    previous = 0.0;
  for (const double scale : scales) {
    CAPTURE(scale);
    const double                    q2               = scale * scale;
    constexpr std::array<double, 3> threshold        = {qed::Q_ref, 0.1056583745, 1.77686};
    double                          inverse_expected = 1.0 / qed::alpha_ref;
    for (const double mass : threshold) {
      if (scale > mass) { inverse_expected -= std::log(q2 / math::pow2(mass)) / (3.0 * math::PI); }
    }
    const double expected = 1.0 / inverse_expected;
    const double coupling = qed::alpha_QED(q2, "LL");

    REQUIRE(std::isfinite(coupling));
    REQUIRE(coupling > 0.0);
    REQUIRE(coupling == Approx(expected).margin(2e-18));
    REQUIRE(coupling >= previous);
    previous = coupling;
  }

  REQUIRE(qed::get_N_leptons(0.05) == 1);
  REQUIRE(qed::get_N_leptons(0.2) == 2);
  REQUIRE(qed::get_N_leptons(2.0) == 3);
  for (const double threshold : {0.1056583745, 1.77686}) {
    const double below = qed::alpha_QED(math::pow2(threshold * (1.0 - 1e-10)), "LL");
    const double above = qed::alpha_QED(math::pow2(threshold * (1.0 + 1e-10)), "LL");
    REQUIRE(above == Approx(below).epsilon(1e-10));
  }
  REQUIRE_THROWS_AS(qed::alpha_QED(1.0, "BAD"), std::invalid_argument);
  REQUIRE_THROWS_AS(qed::alpha_QED(1.0, "BAD", alpha_mg), std::invalid_argument);
  REQUIRE_THROWS_AS(qed::alpha_QED(-1.0, "LL"), std::invalid_argument);
  REQUIRE_THROWS_AS(qed::alpha_QED(std::numeric_limits<double>::infinity(), "LL"), std::invalid_argument);
}

// Test form factors
//

// Check reused standalone and family amplitudes restore the selected QED scheme
TEST_CASE("MG5 amplitudes can switch QED schemes between evaluations", "[MG2GRA][QED][cache]") {
  const auto directory = std::filesystem::path("tmp") / "test_mg5_qed_schemes";
  std::filesystem::create_directories(directory);
  for (const auto& filename : {"NUMERICS.json", "CON_MP.json", "CON_XP.json", "CON_GP.json", "CON_TP.json"}) {
    std::filesystem::copy_file(std::filesystem::path(modelfile).parent_path() / filename, directory / filename,
                               std::filesystem::copy_options::overwrite_existing);
  }
  // Load independent immutable tune snapshots from the same general card
  const auto tune = [&](const std::string& scheme) {
    auto card                            = nlohmann::json::parse(gra::aux::GetInputData(modelfile));
    card["PARAM_STRUCTURE"]["QED_alpha"] = scheme;
    const auto    path                   = directory / "GENERAL.json";
    std::ofstream output(path);
    output << card;
    output.close();
    return std::make_shared<gra::MModelCache>(gra::MModelTune::Load(path.string()));
  };
  const auto zero = tune("ZERO");
  const auto mg   = tune("MG");
  // Compare a reused real amplitude with a freshly constructed amplitude in each scheme
  const auto check = [&](auto &amplitude, const std::string &name, bool coherent) {
    CAPTURE(name, coherent);
    auto lts = PhotonCovariancePoint(name);
    double zero_amp2 = 0.0;
    for (const auto &cache : {zero, mg, zero}) {
      lts.model_cache = cache;
      std::decay_t<decltype(amplitude)> fresh;
      auto reference_lts = lts;
      const auto reference = fresh.Evaluate(reference_lts, lts.alphaQCD, coherent);
      const auto result = amplitude.Evaluate(lts, lts.alphaQCD, coherent);
      REQUIRE(reference.Valid());
      REQUIRE(result.Valid());
      REQUIRE(reference.amp2 > 0.0);
      REQUIRE(result.amp2 == Approx(reference.amp2).epsilon(2e-12));
      if (cache == zero) {
        zero_amp2 = result.amp2;
      } else {
        REQUIRE(std::abs(result.amp2 - zero_amp2) > 1e-3 * zero_amp2);
      }
    }
  };
  for (const bool coherent : {false, true}) {
    AMP_MG5_yy_ll standalone;
    gra::AMP_MG5_yy_jj family;
    check(standalone, "yy_ll", coherent);
    check(family, "yy_jj", coherent);
  }

  // Verify independent SLHA charge inputs through the massive QED cross section
  AMP_MG5_yy_ll massive;
  for (const double inverse_alpha : {120.0, 140.0}) {
    SLHAReader card(gra::aux::ResolveProjectPath("MG5cards/Photon/yy_ll/param_card.dat"));
    card.set_block_entry("mass", 11, 2.0);
    card.set_block_entry("yukawa", 11, 2.0);
    card.set_block_entry("sminputs", 1, inverse_alpha);
    Parameters_sm_lepton_masses parameters;
    parameters.setIndependentParameters(card);
    card.set_block_entry("mass", 24, parameters.Particles().at(24).mass);
    massive.InitParameters(card);
    auto lts = MassivePhotonPair(2.0, 0.6, 0.3, 0.4);
    lts.model_cache = mg;
    const auto result = massive.Evaluate(lts, 0.0, false);
    REQUIRE(result.Valid());
    const double expected = BreitWheelerAmp2(0.6, 0.3) / pow2(inverse_alpha * gra::qed::alpha_QED());
    REQUIRE(result.amp2 == Approx(expected).epsilon(2.0e-9));
  }
}

// Reject incompatible model restrictions before evaluating an event
TEST_CASE("MG5 model initialization preserves masses and rejects restricted parameters", "[MG5][model][massive]") {
  AMP_MG5_yy_ll amplitude;
  ConfigureDiracMass(amplitude, 2.0);
  auto valid = MassivePhotonPair(2.0, 0.6, 0.3, 0.4);
  const auto result = amplitude.Evaluate(valid, 0.0, false);
  REQUIRE(result.Valid());
  REQUIRE(result.amp2 == Approx(BreitWheelerAmp2(0.6, 0.3)).epsilon(2.0e-9));
  auto incompatible = MassivePhotonPair(3.0, 0.6, 0.3, 0.4);
  REQUIRE(amplitude.Evaluate(incompatible, 0.0, false).status == gra::mg5helas::EvaluationStatus::KinematicsFailure);
  REQUIRE(amplitude.getMasses()[2] == Approx(2.0));
  REQUIRE(amplitude.Particles().at(11).mass == Approx(2.0));

  for (const int pdg : {1, 22}) {
    SLHAReader card(gra::aux::ResolveProjectPath("MG5cards/Photon/yy_ll/param_card.dat"));
    card.set_block_entry("mass", pdg, 1.0);
    REQUIRE_THROWS_AS(amplitude.InitParameters(card), std::invalid_argument);
  }
  SLHAReader inconsistent(gra::aux::ResolveProjectPath("MG5cards/Photon/yy_ll/param_card.dat"));
  inconsistent.set_block_entry("mass", 24, 1.0);
  REQUIRE_THROWS_AS(amplitude.InitParameters(inconsistent), std::invalid_argument);
}
