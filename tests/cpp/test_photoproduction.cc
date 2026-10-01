// Photoproduction model tests
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include <limits>
#include <tuple>
#include <numeric>
#include <optional>

#include "Graniitti/Math/MIntegration.h"
#include "Graniitti/Regge/MReggeInit.h"
#include "Graniitti/Regge/MReggeGPInit.h"
#include "Graniitti/Nuclear/MConfig.h"
#include "Graniitti/Nuclear/MNucleus.h"
#include "Graniitti/Nuclear/MUPC.h"
#include "Graniitti/Nuclear/MUPCScreen.h"
#include "Graniitti/Photon/MFlux.h"
#include "Graniitti/Regge/MReggeMP.h"
#include "Graniitti/Regge/MReggeSig.h"
#include "Graniitti/Regge/MReggeSource.h"
#include "Graniitti/Regge/MReggeXP.h"
#include "HepMC3/WriterAscii.h"
#include "support/models_test_support.hh"
#include "support/nuclear_test_support.hh"

namespace {

using gra::aux::indices;

// Read the H1 energy scan through the icepack decoder without duplicating measurements
std::vector<std::pair<double, double>> H1PhotoData() {
  const auto data = gra::aux::ExecCommand(
      "python -c 'import json; from icepack.PHOTOPROD._common import hera_reader as r; "
      "t = r.read_table(\"HEPData/PHOTOPROD/HEPData-ins1228913-v1-json/Table1.json\"); "
      "print(json.dumps(list(zip(r.read_bins(t)[0].tolist(), "
      "[r.number(v[\"value\"]) for v in r.table_values(t)]))))'");
  return nlohmann::json::parse(data).get<std::vector<std::pair<double, double>>>();
}


// Integrate the skewed NLO gluon by parts independently of its analytic derivative
// [REFERENCE: arXiv:1307.7099, Eqs. (4), (7), (11) and (17)]
double JMRTSkewedKernel(const gra::MPhotoVMNumerics &num, double x, double qbar2, double kmax) {
  const auto &fit = num.jmrt_2013;
  const double q0 = fit.nlo_q02;
  const auto coupling = [&](double k) {
    const int nf = k >= pow2(num.HeavyQuarkMass(5)) ? 5 : (k >= pow2(num.HeavyQuarkMass(4)) ? 4 : 3);
    const double b0 = 11.0 - 2.0 * nf / 3.0;
    const double c = (102.0 - 38.0 * nf / 3.0) / pow2(b0);
    const double l = std::log(k / pow2(num.running_alpha_lambda_qcd.at(nf)));
    return std::pair<double, double>{4.0 * gra::math::PI / b0 * (1.0 / l - c * std::log(l) / pow2(l)),
        4.0 * gra::math::PI / b0 * (-1.0 / pow2(l) - c * (1.0 - 2.0 * std::log(l)) / gra::math::pow3(l)) / k};
  };
  const auto flux = [&](double k) {
    const double logg = std::log(std::log(k / pow2(fit.nlo_lambda_qcd)) /
                                 std::log(q0 / pow2(fit.nlo_lambda_qcd)));
    const double lambda = fit.nlo.a + 0.5 * std::sqrt(16.0 / 3.0 * logg / std::log(1.0 / x));
    const double skew = std::pow(2.0, 2.0 * lambda + 3.0) * std::tgamma(lambda + 2.5) /
                         (std::sqrt(gra::math::PI) * std::tgamma(lambda + 4.0));
    const double gluon = fit.nlo.normalization * std::pow(x, -fit.nlo.a) * std::pow(k, fit.nlo.b) *
                          std::exp(std::sqrt(16.0 / 3.0 * std::log(1.0 / x) * logg));
    return skew * gluon * std::exp(-3.0 * coupling(qbar2).first * pow2(std::max(0.0, std::log(qbar2 / k))) /
                                    (8.0 * gra::math::PI));
  };
  const auto weight = [&](double k) {
    const auto [alpha, derivative] = coupling(std::max(k, qbar2));
    return std::pair<double, double>{alpha / (qbar2 * (qbar2 + k)),
        ((k < qbar2 ? 0.0 : derivative) - alpha / (qbar2 + k)) / (qbar2 * (qbar2 + k))};
  };
  std::vector<double> bounds = {q0, kmax};
  for (double k : {qbar2, pow2(num.HeavyQuarkMass(4)), pow2(num.HeavyQuarkMass(5))}) {
    if (k > q0 && k < kmax) { bounds.push_back(k); }
  }
  std::sort(bounds.begin(), bounds.end());
  double integral = flux(kmax) * weight(kmax).first - flux(q0) * weight(q0).first;
  for (std::size_t i = 1; i < bounds.size(); ++i) {
    if (bounds[i] <= bounds[i - 1]) { continue; }
    integral -= gra::math::LogMeasureGaussIntegral(256, bounds[i - 1], bounds[i],
        [&](double k) { return flux(k) * weight(k).second * k; });
  }
  return integral + std::log1p(q0 / qbar2) * coupling(std::max(q0, qbar2)).first * flux(q0) / (qbar2 * q0);
}

// Set two on-shell decay momenta at fixed rest-frame angles
void SetPhotoDecay(gra::LORENTZSCALAR &lts, double cosine, double phi) {
  const auto parent = lts.pfinal[0];
  const double energy = parent.M() / 2.0;
  const double p = std::sqrt(pow2(energy) - pow2(lts.decaytree[0].p.mass));
  const double sine = std::sqrt(1.0 - pow2(cosine));
  gra::M4Vec first(p * sine * std::cos(phi), p * sine * std::sin(phi), p * cosine, energy);
  gra::M4Vec second(-first.Px(), -first.Py(), -first.Pz(), energy);
  gra::kinematics::LorentzBoost(parent, parent.M(), first, 1);
  gra::kinematics::LorentzBoost(parent, parent.M(), second, 1);
  lts.decaytree[0].p4 = first;
  lts.decaytree[1].p4 = second;
}

// Boost all external momenta without changing scalar invariants or light-cone fractions
gra::LORENTZSCALAR BoostPhotoEvent(gra::LORENTZSCALAR lts, double rapidity) {
  const gra::M4Vec boost(0.0, 0.0, std::sinh(rapidity), std::cosh(rapidity));
  for (auto *p : {&lts.pbeam1, &lts.pbeam2, &lts.q1, &lts.q2}) {
    gra::kinematics::LorentzBoost(boost, 1.0, *p, 1);
  }
  for (auto &p : lts.pfinal) { gra::kinematics::LorentzBoost(boost, 1.0, p, 1); }
  for (auto &branch : lts.decaytree) { gra::kinematics::LorentzBoost(boost, 1.0, branch.p4, 1); }
  return lts;
}

// Exchange the two beams by a proper rotation and relabel the forward legs
gra::LORENTZSCALAR ExchangePhotoBeams(gra::LORENTZSCALAR lts) {
  for (auto *p : {&lts.pbeam1, &lts.pbeam2, &lts.q1, &lts.q2}) { p->RotateY(gra::math::PI); }
  for (auto &p : lts.pfinal) { p.RotateY(gra::math::PI); }
  for (auto &branch : lts.decaytree) { branch.p4.RotateY(gra::math::PI); }
  std::swap(lts.pbeam1, lts.pbeam2);
  std::swap(lts.beam1, lts.beam2);
  std::swap(lts.pfinal[1], lts.pfinal[2]);
  std::swap(lts.q1, lts.q2);
  std::swap(lts.t1, lts.t2);
  std::swap(lts.qt1, lts.qt2);
  std::swap(lts.xi1, lts.xi2);
  std::swap(lts.has_xi1, lts.has_xi2);
  std::swap(lts.forward_mass2[0], lts.forward_mass2[1]);
  std::swap(lts.excite1, lts.excite2);
  return lts;
}

// Build a direct photoproduction event with a moving central system
gra::LORENTZSCALAR MovingPhotoPairEvent(const int fermion_pdg, const double nominal_mass) {
  gra::LORENTZSCALAR lts      = MakeToyPhotoZFFbar(fermion_pdg, nominal_mass);
  const double       upper_pz = lts.pfinal[1].Pz();
  const double       lower_pz = lts.pfinal[2].Pz();
  lts.pfinal[1].SetPxPyPzM(0.22, -0.09, upper_pz, gra::PDG::mp);
  lts.pfinal[2].SetPxPyPzM(-0.10, 0.04, lower_pz, gra::PDG::mp);
  RefreshToyDerivedKinematicsPreserveDecay(lts);

  const gra::M4Vec system        = lts.pfinal[0];
  const double     mass          = system.M();
  const double     daughter_mass = lts.decaytree[0].p.mass;
  const double     energy        = 0.5 * mass;
  const double     momentum      = std::sqrt(std::max(0.0, energy * energy - daughter_mass * daughter_mass));
  constexpr double cosine        = 0.23;
  const double     sine          = std::sqrt(1.0 - cosine * cosine);
  constexpr double phi           = 0.41;
  gra::M4Vec       fermion(momentum * sine * std::cos(phi), momentum * sine * std::sin(phi), momentum * cosine, energy);
  gra::M4Vec       antifermion(-fermion.Px(), -fermion.Py(), -fermion.Pz(), energy);
  gra::kinematics::LorentzBoost(system, mass, fermion, 1);
  gra::kinematics::LorentzBoost(system, mass, antifermion, 1);
  lts.decaytree[0].p4 = fermion;
  lts.decaytree[1].p4 = antifermion;
  RefreshToyDerivedKinematicsPreserveDecay(lts);
  return lts;
}

// Replace one proton beam by the corresponding antiproton
void SetPhotoAntiprotonBeam(gra::LORENTZSCALAR &lts, const int leg) {
  if (leg == 1) {
    lts.beam1 = lts.PDG.FindByPDG(-gra::PDG::PDG_p);
    return;
  }
  if (leg == 2) {
    lts.beam2 = lts.PDG.FindByPDG(-gra::PDG::PDG_p);
    return;
  }
  throw std::invalid_argument("SetPhotoAntiprotonBeam: leg should be 1 or 2");
}

// Construct a compact oxygen nucleus for photoproduction sector tests
std::shared_ptr<const gra::nuclear::MNucleus> MakePhotoOxygen() {
  gra::nuclear::NucleusParam param;
  param.pdg    = gra::nuclear::EncodeNuclearPDG(16, 8);
  param.mass   = 14.899;
  param.charge = {2.608, 0.513, 11.842, 64, 2048, 3.0, 4097, 1.0e-7};
  param.matter = param.charge;
  return std::make_shared<const gra::nuclear::MNucleus>(param);
}

// Build one symmetric ion event while retaining the direct decay kinematics
gra::LORENTZSCALAR MakeToyIonPhotoPair(const int fermion_pdg, const double central_mass,
                                       const std::shared_ptr<const gra::nuclear::MNucleus> &nucleus) {
  gra::LORENTZSCALAR lts = MakeToyPhotoZFFbar(fermion_pdg, central_mass);
  gra::MParticle     ion = lts.beam1;
  ion.pdg                = nucleus->ID().pdg;
  ion.name               = "O16";
  ion.mass               = nucleus->Mass();
  ion.chargeX3           = 3 * nucleus->Charge();
  ion.spinX2             = 0;
  ion.color              = 0;
  lts.beam1              = ion;
  lts.beam2              = ion;

  constexpr double beam_energy   = 8000.0;
  constexpr double recoil_pt     = 0.18;
  const double     beam_pz       = std::sqrt(gra::math::pow2(beam_energy) - gra::math::pow2(ion.mass));
  const double     recoil_energy = beam_energy - 0.5 * central_mass;
  const double     recoil_pz =
      std::sqrt(gra::math::pow2(recoil_energy) - gra::math::pow2(ion.mass) - gra::math::pow2(recoil_pt));
  lts.pbeam1    = gra::M4Vec(0.0, 0.0, beam_pz, beam_energy);
  lts.pbeam2    = gra::M4Vec(0.0, 0.0, -beam_pz, beam_energy);
  lts.pfinal[1] = gra::M4Vec(recoil_pt, 0.0, recoil_pz, recoil_energy);
  lts.pfinal[2] = gra::M4Vec(-recoil_pt, 0.0, -recoil_pz, recoil_energy);
  RefreshToyDerivedKinematicsPreserveDecay(lts);
  return lts;
}

// Construct one low-cost ion-ion UPC model with selected sector controls
std::shared_ptr<const gra::nuclear::MUPC> MakePhotoSectorUPC(
    const std::shared_ptr<const gra::nuclear::MNucleus> &nucleus, const gra::nuclear::CoherenceType emission,
    const gra::nuclear::CoherenceType target, const gra::nuclear::SurvivalType survival) {
  gra::nuclear::UPCParam param;
  param.loop.radial_integrator  = "GL";
  param.loop.azimuth_integrator = "Trap";
  param.loop.radial_map         = gra::math::RadialMap::Square;
  param.emission                = {emission, emission};
  param.target                  = {target, target};
  param.survival                = survival;
  gra::test::Glauber(param.glauber, 0.0);
  param.glauber.b_max              = 20.0;
  param.glauber.q_max              = 3.0;
  param.glauber.b_nodes            = 32;
  param.glauber.q_nodes            = 32;
  param.convolution.b_max          = 20.0;
  param.convolution.smooth_b_nodes = 32;
  param.convolution.sample_b_nodes = 32;
  gra::test::Convolution(param.convolution);
  param.loop.r_max                = 0.2;
  param.loop.radial_intervals     = 8;
  param.loop.azimuth_nodes        = 4;
  param.convolution.b_phi_nodes   = 5;
  param.photo[0].b_nodes          = 32;
  param.photo[1].b_nodes          = 32;
  param.photo[0].z_nodes          = 32;
  param.photo[1].z_nodes          = 32;
  for (auto &photo : param.photo) { gra::test::Photo(photo); }

  const bool sampled = survival == gra::nuclear::SurvivalType::MCGGCF ||
                       (emission != gra::nuclear::CoherenceType::Coherent &&
                        target != gra::nuclear::CoherenceType::Coherent);
  if (sampled) {
    param.structure = gra::nuclear::StructureType::Nucleon;
    gra::nuclear::ConfigParam config;
    config.count             = 4;
    config.max_trials        = 10000;
    config.neutron_nodes     = 4096;
    config.density_rel_tol   = 1.0e-6;
    config.negative_norm_tol = 1.0e-7;
    config.norm_tol          = 2.0e-5;
    param.config             = config;
  }
  const auto model = std::make_shared<const gra::nuclear::MUPC>(
      std::array{gra::nuclear::BeamType::Nucleus, gra::nuclear::BeamType::Nucleus}, std::array{nucleus, nucleus},
      param);
  if (!sampled) { return model; }
  gra::MRandom random;
  random.SetSeed(99173);
  return model->Sample(random);
}

enum class ToyPhotoType { VM, Z };

class ToyNuclearPhotoScreeningProcess : public gra::MFactorized {
 public:
  // Construct one prepared ion-ion ygg process with a UPC convolution
  ToyNuclearPhotoScreeningProcess(gra::MModelTunePtr                                   model_tune,
                                  const std::shared_ptr<const gra::nuclear::MNucleus> &nucleus,
                                  std::shared_ptr<const gra::nuclear::MUPC>            upc,
                                  const ToyPhotoType                                   type = ToyPhotoType::VM) {
    state.screening = true;
    SetModelTune(std::move(model_tune));
    ProcPtr.CHANNEL           = type == ToyPhotoType::VM ? "jpsi" : "Z";
    ProcPtr.ISTATE            = "ygg";
    const double central_mass = type == ToyPhotoType::VM ? 3.096900 : 91.1876;
    state.lts                 = MakeToyIonPhotoPair(13, central_mass, nucleus);
    state.lts.upc_model       = std::move(upc);
    state.lts.upc_event       = state.lts.upc_model;
    PrepareScreeningPoint(state);
    state.lts.pfinal_orig     = state.lts.pfinal;
    if (type == ToyPhotoType::VM) {
      photovm = std::make_unique<gra::MPhotoVM>(state.lts, state.model_tune, "jpsi",
                                                gra::MPhotoVM::ProcessDefinitionFor("jpsi"));
    } else {
      photoz = std::make_unique<gra::MPhotoZ>(state.lts, state.model_tune, gra::MPhotoZ::ProcessDefinitionFor());
    }
  }

  // Evaluate the ygg amplitude with optional UPC screening
  double ScreenedAmp2(const bool include_screening = true) { return ScreenedAmplitudeSquared(include_screening); }

 protected:
  // Evaluate the selected ygg hard amplitude at the current recoil point
  double EvaluateBareAmplitude() override {
    return photovm != nullptr ? photovm->Amp2(state.lts) : photoz->Amp2(state.lts);
  }

 private:
  std::unique_ptr<gra::MPhotoVM> photovm;
  std::unique_ptr<gra::MPhotoZ>  photoz;
};

}  // namespace

TEST_CASE("Analytic photon amplitudes expose exact processes", "[gra::MGamma][gra::MPhotoZ][gra::MPhotoVM][process]") {
  // Build a direct two-particle topology from bare PDG codes
  const auto direct_pair = [](const int first_pdg, const int second_pdg) {
    gra::MDecayBranch first;
    first.p.pdg    = first_pdg;
    first.p.spinX2 = 1;
    gra::MDecayBranch second;
    second.p.pdg    = second_pdg;
    second.p.spinX2 = 1;
    return std::vector<gra::MDecayBranch>{first, second};
  };

  const auto gamma_pair  = gra::MGamma::ProcessDefinitionFor(gra::MGammaMode::FermionPair);
  const auto muon_pair   = direct_pair(-13, 13);
  const auto gamma_match = gamma_pair->MatchProcess(muon_pair);
  REQUIRE(gamma_match.has_value());
  REQUIRE(gamma_match->decay_structure.type == gra::DecayType::Full);
  REQUIRE(gamma_pair->MatchProcess(direct_pair(6, -6)).has_value());
  REQUIRE(gamma_pair->MatchProcess(direct_pair(gra::PDG::PDG_monopole, -gra::PDG::PDG_monopole)).has_value());
  REQUIRE_FALSE(gamma_pair->MatchProcess(direct_pair(14, -14)).has_value());
  REQUIRE_FALSE(gamma_pair->MatchProcess(direct_pair(13, -11)).has_value());
  auto wrong_spin        = muon_pair;
  wrong_spin[0].p.spinX2 = 0;
  REQUIRE_FALSE(gamma_pair->MatchProcess(wrong_spin).has_value());
  auto nested_pair = muon_pair;
  nested_pair[0].legs.push_back(muon_pair[0]);
  REQUIRE_FALSE(gamma_pair->MatchProcess(nested_pair).has_value());

  gra::MDecayBranch leaf;
  leaf.p.pdg                                         = 22;
  const std::vector<gra::MDecayBranch> nonempty_tree = {leaf};
  for (const auto mode :
       {gra::MGammaMode::Generic, gra::MGammaMode::Higgs, gra::MGammaMode::Monopolium, gra::MGammaMode::Flux}) {
    const auto definition = gra::MGamma::ProcessDefinitionFor(mode);
    const auto match      = definition->MatchProcess(nonempty_tree);
    REQUIRE(match.has_value());
    REQUIRE(match->decay_structure.type == gra::DecayType::None);
    REQUIRE_FALSE(definition->MatchProcess({}).has_value());
  }

  const auto photoz = gra::MPhotoZ::ProcessDefinitionFor();
  REQUIRE(photoz->MatchProcess(muon_pair).has_value());
  REQUIRE(photoz->MatchProcess(direct_pair(5, -5)).has_value());
  REQUIRE_FALSE(photoz->MatchProcess(direct_pair(6, -6)).has_value());
  REQUIRE_FALSE(photoz->MatchProcess(nested_pair).has_value());
  REQUIRE(photoz->MatchProcess(muon_pair)->decay_structure == gra::MPhotoZ::DirectDecayStructure());

  const auto photovm = gra::MPhotoVM::ProcessDefinitionFor("jpsi");
  REQUIRE(photovm->MatchProcess(direct_pair(-15, 15)).has_value());
  REQUIRE_FALSE(photovm->MatchProcess(direct_pair(5, -5)).has_value());
  const auto photovm_patterns = photovm->Processes();
  REQUIRE(photovm_patterns.size() == 1);
  REQUIRE(photovm_patterns[0].process_name.find("jpsi") != std::string::npos);
  REQUIRE_THROWS_AS(gra::MPhotoVM::ProcessDefinitionFor("unsupported"), std::invalid_argument);
}

// Check continuum process caches before the amplitudes index their rows
TEST_CASE("MRegge coherently sums the two automatic mixed-exchange orderings",
          "[gra::MRegge][photoproduction][symmetry][physics]") {
  ToyHelicityProcess proc;
  ConfigureToyProductionProcess(proc, "MP", "RES", "pi+ pi-");
  MRandom        rng;
  gra::PARAM_RES rho = gra::resonance::Read("RES/rho_770.json", rng, gra::ReggeProductionModel::MP);
  proc.SetResonances({{"rho_770", rho}});
  REQUIRE_NOTHROW(proc.InitializeProcessAmplitude());
  const gra::PARAM_RES configured = proc.GetResonances().at("rho_770");
  REQUIRE(configured.production.size() == 2);

  // Keep exactly one expanded beam ordering while preserving its matching spin
  // tensor
  auto select_channel = [](const gra::PARAM_RES &input, std::size_t index) {
    gra::PARAM_RES output  = input;
    output.production = {input.production.at(index)};
    return output;
  };

  const auto         model_tune = gra::MModelTune::Load(modelfile);
  gra::LORENTZSCALAR full_lts   = MakeToyCoherentPhotonLTS();
  gra::RequireModelCache(full_lts.model_cache, model_tune, "test coherent Regge sum");
  gra::LORENTZSCALAR upper_lts = full_lts;
  gra::LORENTZSCALAR lower_lts = full_lts;
  gra::PARAM_RES     full      = configured;
  gra::PARAM_RES     upper     = select_channel(configured, 0);
  gra::PARAM_RES     lower     = select_channel(configured, 1);
  gra::MRegge        regge(full_lts, model_tune, gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));

  TestReggeRes(regge, full_lts, full, gra::ReggeProductionModel::MP);
  TestReggeRes(regge, upper_lts, upper, gra::ReggeProductionModel::MP);
  TestReggeRes(regge, lower_lts, lower, gra::ReggeProductionModel::MP);
  REQUIRE(full_lts.hamp.size() == upper_lts.hamp.size());
  REQUIRE(full_lts.hamp.size() == lower_lts.hamp.size());

  double upper_norm2 = 0.0;
  double lower_norm2 = 0.0;
  for (std::size_t i = 0; i < full_lts.hamp.size(); ++i) {
    const std::complex<double> expected = upper_lts.hamp[i] + lower_lts.hamp[i];
    const double               scale    = std::max({1.0, std::abs(full_lts.hamp[i]), std::abs(expected)});
    CHECK(std::abs(full_lts.hamp[i] - expected) < 1e-11 * scale);
    upper_norm2 += std::norm(upper_lts.hamp[i]);
    lower_norm2 += std::norm(lower_lts.hamp[i]);
  }
  CHECK(upper_norm2 > 0.0);
  CHECK(lower_norm2 > 0.0);
}

// Check each production basis carries the same physical transverse HERA limit
TEST_CASE("Identical rho input has one MP XP GP HERA residue",
          "[gra::MRegge][photoproduction][normalization][physics]") {
  const std::array models            = {gra::ReggeProductionModel::MP, gra::ReggeProductionModel::XP,
                                        gra::ReggeProductionModel::GP};
  MRandom rng;
  auto input = gra::resonance::Read("RES/rho_770.json", rng, gra::ReggeProductionModel::XP);
  input.XP.phi = 0.31;
  for (auto &channel : input.XP.channels) {
    for (auto &coupling : channel.g_helicity) { coupling = {0.37, 0.11}; }
  }
  static_cast<gra::RES_PRODUCTION_MODEL &>(input.MP) = input.XP;
  input.GP = input.XP;
  for (auto &channel : input.GP.channels) {
    for (auto &exchange : channel.exchange) {
      if (exchange == 991) { exchange = 990; }
    }
  }
  double physical_residue2 = -1.0;
  for (const auto model : models) {
    const std::string name = gra::ReggeProductionModelName(model);
    CAPTURE(name);
    ToyHelicityProcess process;
    ConfigureToyProductionProcess(process, name, "RES", "pi+ pi-");
    auto resonance = input;
    resonance.production_model = model;
    process.SetResonances({{"rho_770", resonance}});
    REQUIRE_NOTHROW(process.InitializeProcessAmplitude());

    const auto &rho = process.GetResonances().at("rho_770");
    REQUIRE(rho.production.size() == 2);
    for (const auto &channel : indices(rho.production)) {
      double reference_density = 0.0;
      if (model == gra::ReggeProductionModel::GP) {
        const auto &hel = rho.production.at(channel).hel;
        REQUIRE(hel.UsesHelicityCouplings());
        const std::size_t zero         = gra::gpom::AnalyticMIndex(0, hel.analytic_MMAX, "GP HERA transverse limit");
        const bool        upper_photon = rho.production[channel].tree[0].p.pdg == gra::PDG::PDG_gamma;
        for (const int photon_m : {-1, 1}) {
          const std::size_t transverse =
              gra::gpom::AnalyticMIndex(photon_m, hel.analytic_MMAX, "GP HERA transverse limit");
          reference_density += upper_photon ? std::norm(hel.T[transverse][zero]) : std::norm(hel.T[zero][transverse]);
        }
        reference_density *= 0.5;
      } else {
        reference_density = ConfiguredPoleReferenceDensity(rho, model, channel);
      }
      const double residue2 = reference_density;
      CAPTURE(channel, reference_density, residue2);
      if (physical_residue2 < 0.0) { physical_residue2 = residue2; }
      CHECK(residue2 == Approx(physical_residue2).epsilon(4.0e-6));
    }
  }
}

// Compare complete vector amplitudes, including configured helicities, spin projections and decay
TEST_CASE("Vector MP and XP amplitudes have the same physical residue",
          "[gra::MRegge][photoproduction][normalization][regression]") {
  const auto [card, nstars] = GENERATE(
      std::make_pair("rho_770", 0), std::make_pair("phi_1020", 0), std::make_pair("jpsi", 0),
      std::make_pair("psi2S", 0), std::make_pair("Y1S", 0), std::make_pair("Y2S", 0), std::make_pair("Y3S", 0),
      std::make_pair("rho_770_odd", 0), std::make_pair("phi_1020_odd", 0),
      std::make_pair("rho_770", 1), std::make_pair("phi_1020", 1), std::make_pair("jpsi", 1));
  const std::string photon_vertex = GENERATE("EPA", "QED");
  MRandom rng;
  auto input = gra::resonance::Read("RES/" + std::string(card) + ".json", rng, gra::ReggeProductionModel::XP);
  input.XP.phi = 0.31;
  for (auto &channel : input.XP.channels) {
    for (auto &coupling : channel.g_helicity) { coupling = {0.37, 0.11}; }
  }
  // Share the mass, width, couplings, phase and form factors between spin constructions
  static_cast<gra::RES_PRODUCTION_MODEL &>(input.MP) = input.XP;
  const auto amplitude = [&](gra::ReggeProductionModel model, const std::string& frame, double angle) {
    ToyHelicityProcess process;
    ConfigureToyProductionProcess(process, gra::ReggeProductionModelName(model), "RES", "mu+ mu-");
    auto resonance = input;
    resonance.production_model = model;
    auto lts       = MovingPhotoPairEvent(13, resonance.p.mass);
    if (std::string(card).find("_odd") == std::string::npos) {
      lts.beam1 = lts.PDG.FindByPDG(-11);
      lts.pbeam1.SetPxPyPzM(0.0, 0.0, lts.pbeam1.Pz(), lts.beam1.mass);
      lts.pfinal[1].SetPxPyPzM(lts.pfinal[1].Px(), lts.pfinal[1].Py(), lts.pfinal[1].Pz(), lts.beam1.mass);
    }
    if (nstars) { SetToyPhotoForwardExcitation(lts, 2, 2.0); }
    SetPhotoDecay(lts, 0.23, 0.41);
    RefreshToyDerivedKinematicsPreserveDecay(lts);
    process.state.lts = lts;
    process.SetModelTune(gra::MModelTune::Load(modelfile));
    process.SetExcitation(nstars);
    process.SetResonances({{card, resonance}});
    process.PrepareRun();
    auto event                   = RotateToyEventAroundZ(process.state.lts, angle);
    event.process.PHOTON_VERTEX = photon_vertex;
    event.process.MP_FRAME      = frame;
    REQUIRE(process.ProcPtr.GetBareAmplitude2(event) > 0.0);
    return event.hamp;
  };
  const auto reference = amplitude(gra::ReggeProductionModel::XP, "CM", 0.0);
  for (const std::string frame : {"CM", "CS", "HX"}) {
    for (const double angle : {0.0, 0.57}) {
      CAPTURE(card, nstars, photon_vertex, frame, angle);
      RequireVectorNear(amplitude(gra::ReggeProductionModel::MP, frame, angle), reference, 1.0e-7);
    }
  }
}

// Integrate the local double-log trajectory independently of the native energy primitive
double PhotoDLogIntegral(const gra::regge::PhotoDLog &term, double s, double s0) {
  if (std::abs(std::log(s / s0)) < std::numeric_limits<double>::epsilon()) { return 0.0; }
  const auto derivative = [&](double w2) {
    return term.c / 4.0 * (1.0 / std::sqrt(std::log(w2 / term.scale2)) - 1.0 / term.root0);
  };
  return (s > s0 ? 1.0 : -1.0) *
      gra::math::LogMeasureGaussIntegral(64, std::min(s, s0), std::max(s, s0), derivative);
}

// Check the channel kernels against the forward HERA normalization
TEST_CASE("MRegge photoproduction amplitudes use factorized HERA normalizations",
          "[gra::MRegge][photoproduction][physics]") {
  gra::LORENTZSCALAR        lts = MakeToyCoherentPhotonLTS();
  gra::MRegge               regge(lts, gra::MModelTune::Load(modelfile),
                                  gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));
  const auto               &soft    = *regge.SoftModelHandle();
  const gra::SoftExchangeId pomeron = soft.ExchangeId("P");
  const auto                param   = gra::regge::ReadParam({211, -211}, LoadedPDGTable(), *regge.ModelTuneHandle());
  for (const auto &mapping : param.exchanges) {
    REQUIRE(regge.SoftModelHandle()->Exchange(mapping.soft_exchange).Eta() == gra::EtaMode::Rotating);
  }
  const auto tune_dir = regge.ModelTuneHandle()->Directory();
  MRandom              rng;
  const gra::PARAM_RES rho   = gra::resonance::Read("RES/rho_770.json", rng, gra::ReggeProductionModel::GP, tune_dir);
  const gra::PARAM_RES phi   = gra::resonance::Read("RES/phi_1020.json", rng, gra::ReggeProductionModel::GP, tune_dir);
  const gra::PARAM_RES jpsi  = gra::resonance::Read("RES/jpsi.json", rng, gra::ReggeProductionModel::GP, tune_dir);
  const gra::PARAM_RES psi2s = gra::resonance::Read("RES/psi2S.json", rng, gra::ReggeProductionModel::GP, tune_dir);
  const gra::PARAM_RES y1s   = gra::resonance::Read("RES/Y1S.json", rng, gra::ReggeProductionModel::GP, tune_dir);
  const gra::PARAM_RES y2s   = gra::resonance::Read("RES/Y2S.json", rng, gra::ReggeProductionModel::GP, tune_dir);
  const gra::PARAM_RES y3s   = gra::resonance::Read("RES/Y3S.json", rng, gra::ReggeProductionModel::GP, tune_dir);

  struct HERAChannel {
    int                   pdg;
    const gra::PARAM_RES *resonance;
  };
  const std::array<HERAChannel, 7> channels = {{{113, &rho}, {333, &phi}, {443, &jpsi},
                                                {100443, &psi2s}, {553, &y1s}, {100553, &y2s}, {200553, &y3s}}};
  const double                     gev2_to_ub      = gra::PDG::GeV2barn * 1.0e6;
  const double                     proton_coupling = soft.PhysicalResidue(pomeron, 0.0);
  for (const auto &data : channels) {
    CAPTURE(data.pdg);
    const auto  &channel = gra::regge::Photo(param, data.pdg);
    REQUIRE(data.resonance->GP.channels.size() == 1);
    const auto& input = data.resonance->GP.channels.front();
    REQUIRE(input.g_helicity.size() == 1);
    const double expected_coupling = std::abs(input.g_helicity.front());
    const double forward = pow2(expected_coupling * proton_coupling) * gev2_to_ub / (16.0 * gra::math::PI);

    const double               s0 = pow2(channel.W0);
    const std::complex<double> forward_amplitude =
        expected_coupling * soft.PhysicalResidue(pomeron, 0.0) * regge.PhotoKernel(s0, 0.0, data.pdg);
    const double reconstructed_forward = std::norm(forward_amplitude) / (16.0 * gra::math::PI * pow2(s0)) * gev2_to_ub;
    CHECK(std::abs(soft.PhysicalResidue(pomeron, 0.0) * regge.PhotoKernel(s0, 0.0, data.pdg)) ==
          Approx(s0 * proton_coupling).epsilon(1.0e-12));
    const std::complex<double> forward_phase =
        regge.PhotoKernel(s0, 0.0, data.pdg) / std::abs(regge.PhotoKernel(s0, 0.0, data.pdg));
    const std::complex<double> expected_forward_phase = gra::regge::EtaPhase(
        channel.a0, soft.Exchange(param.exchanges.at(param.pomeron_trajectory).soft_exchange).Signature());
    CHECK(std::real(forward_phase) == Approx(std::real(expected_forward_phase)).epsilon(1.0e-12));
    CHECK(std::imag(forward_phase) == Approx(std::imag(expected_forward_phase)).epsilon(1.0e-12));
    CHECK(reconstructed_forward == Approx(forward).epsilon(1.0e-12));

    // Extract the HERA coupling from the full GP production spin sum
    ToyHelicityProcess process;
    ConfigureToyProductionProcess(process, "GP", "RES", "mu+ mu-");
    process.SetResonances({{"vector", *data.resonance}});
    process.InitializeProcessAmplitude();
    lts.process.MMAX = process.state.lts.process.MMAX;
    gra::gpom::AmpCache cache;
    const auto          matrices = gra::gpom::Resonance(lts, param, process.GetResonances().at("vector"), &cache);
    const double        flux =
        gra::spin::SourceSpinAveragedDensity(cache.sources.front().upper_residue, 2, "HERA forward source");
    const double gp_coupling2 = 0.25 * matrices.front().FrobNorm2() / flux;
    const double gp_forward   = gp_coupling2 * pow2(proton_coupling) * gev2_to_ub / (16.0 * gra::math::PI);
    CHECK(gp_forward == Approx(forward).epsilon(0.002));

    const double               s = 1.44 * s0;
    const double               t = -0.15;
    const std::complex<double> shifted_amplitude =
        expected_coupling * soft.PhysicalResidue(pomeron, t) * regge.PhotoKernel(s, t, data.pdg);
    const double dsigma_ratio = (std::norm(shifted_amplitude) / pow2(s)) / (std::norm(forward_amplitude) / pow2(s0));
    const double alpha_t      = channel.a0 + channel.ap * t;
    const double proton_form_factor = soft.PhysicalResidue(pomeron, t) / proton_coupling;
    const double profile2 = gra::math::pow2(proton_form_factor * gra::form::ExpSlopeAmplitude(channel.B_gammaPV, t));
    const double evolution = channel.dlog ? std::exp(2 * PhotoDLogIntegral(*channel.dlog, s, s0)) : 1.0;
    const double expected_ratio = std::pow(s / s0, 2.0 * (alpha_t - 1.0)) * profile2 * evolution;
    CHECK(dsigma_ratio == Approx(expected_ratio).epsilon(1.0e-12));

  }
}

// Compute the spin averaged photon and proton vertex normalization in microbarn
double PhotoResidueNorm(gra::MRegge& regge, gra::LORENTZSCALAR& lts, const gra::PARAM_RES& res) {
  const auto& production = res.production.front();
  double      coupling   = 0.0;
  if (res.production_model == gra::ReggeProductionModel::GP) {
    const auto          param = ReggeParametersForTest(regge, lts);
    gra::gpom::AmpCache cache;
    const auto          matrices = gra::gpom::Resonance(lts, *param, res, &cache);
    const auto&         source   = cache.sources.front();
    const auto&         photon   = production.tree[0].p.pdg == 22 ? source.upper_residue : source.lower_residue;
    const double        flux     = gra::spin::SourceSpinAveragedDensity(photon, 2, "HERA photon source");
    coupling                     = std::sqrt(0.25 * matrices.front().FrobNorm2() / flux);
  } else {
    coupling = gra::rspin::PhotoCoupling(production);
  }
  const auto& soft = *regge.SoftModelHandle();
  return pow2(coupling * soft.PhysicalResidue(soft.ExchangeId("P"), 0.0)) * gra::PDG::GeV2barn * 1.0e6 /
         (16.0 * gra::math::PI);
}

// Check double-log evolution through independent energy integration and the complex signature phase
TEST_CASE("Photoproduction double-log evolution preserves HERA normalization and dispersion phase",
          "[gra::MRegge][photoproduction][normalization][phase][physics]") {
  const auto tune = WriteModifiedPhotoVMTune("photo_dlog", [](auto &card) {
    for (const auto &row : card["PARAM_REGGE"]["photoprod"]) {
      if (row[0].template get<int>() == 443) {
        card["PARAM_REGGE"]["photoprod_dlog"] = {{443, pow2(row[1].template get<double>()) / std::exp(5.0), 2.3}};
      }
    }
  });
  auto lts = MakeToyCoherentPhotonLTS();
  const auto model = gra::MModelTune::Load(tune.second);
  gra::MRegge regge(lts, model, gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "double log"));
  const auto param = gra::regge::ReadParam({211, -211}, LoadedPDGTable(), *model);
  const auto &photo = gra::regge::Photo(param, 443);
  REQUIRE(photo.dlog.has_value());
  const auto &evolution = *photo.dlog;
  const double s0 = pow2(photo.W0);
  const auto &exchange = regge.SoftModelHandle()->Exchange(param.exchanges[param.pomeron_trajectory].soft_exchange);
  const auto anchor = static_cast<double>(exchange.ResidueSign()) *
      gra::regge::EtaFactor(photo.a0, photo.a0, exchange.Signature(), param.photoprod_eta_mode);
  RequireComplexNear(regge.PhotoKernel(s0, 0.0, 443) / s0, anchor, 1.0e-12);
  for (const double ratio : {0.2, 1.0 - 1.0e-6, 1.0 + 1.0e-6, 2.0, 20.0}) {
    const double s = s0 * ratio;
    const auto derivative = [&](double w2) {
      return photo.a0 - 1.0 + evolution.c / 4.0 *
          (1.0 / std::sqrt(std::log(w2 / evolution.scale2)) - 1.0 / evolution.root0);
    };
    const double integral = (s > s0 ? 1.0 : -1.0) *
        gra::math::LogMeasureGaussIntegral(64, std::min(s, s0), std::max(s, s0), derivative);
    const auto amplitude = regge.PhotoKernel(s, 0.0, 443) / s;
    CHECK(std::log(std::abs(amplitude)) == Approx(integral).margin(1.0e-11));
    const double step = 1.0e-5;
    const double power = (std::log(std::abs(regge.PhotoKernel(s * std::exp(step), 0.0, 443)) / (s * std::exp(step))) -
                          std::log(std::abs(regge.PhotoKernel(s * std::exp(-step), 0.0, 443)) / (s * std::exp(-step)))) / (2 * step);
    const auto phase = static_cast<double>(exchange.ResidueSign()) *
        gra::regge::EtaFactor(1 + power, 1 + power, exchange.Signature(), param.photoprod_eta_mode);
    RequireComplexNear(amplitude / std::abs(amplitude), phase, 1.0e-9);
    CHECK(std::abs(regge.PhotoKernel(s, -0.2, 443) / regge.PhotoKernel(s, 0.0, 443)) ==
          Approx(std::exp(-0.1 * photo.B_gammaPV - 0.2 * photo.ap * std::log(s / s0))).epsilon(1.0e-12));
  }
  CHECK_THROWS_AS(regge.PhotoKernel(evolution.scale2 / 2, 0, 443), gra::AmplitudeFailure);
}

// Reject malformed double-log input before sampling and omit pure-zero corrections
TEST_CASE("Photoproduction double-log input is validated once", "[photoproduction][input]") {
  const std::vector<nlohmann::json> invalid = {
      {{443, -1.0, 1.0}}, {{443, 1.0, -1.0}}, {{443, 1.0}},
      {{443, 1.0, 1.0}, {443, 2.0, 1.0}}, {{999999, 1.0, 1.0}}};
  for (const auto &i : indices(invalid)) {
    const auto path = WriteModifiedPhotoVMTune("dlog_invalid_" + std::to_string(i), [&](auto &card) {
      card["PARAM_REGGE"]["photoprod_dlog"] = invalid[i];
    });
    CHECK_THROWS_AS(gra::regge::ReadParam({211, -211}, LoadedPDGTable(), *gra::MModelTune::Load(path.second)), std::invalid_argument);
  }
  const auto path = WriteModifiedPhotoVMTune("dlog_zero", [](auto &card) {
    card["PARAM_REGGE"]["photoprod_dlog"] = {{443, 1.0, 0.0}};
    card["PARAM_REGGE"]["photoprod_diss_dlog"] = {{443, 1.0, 0.0}};
  });
  const auto param = gra::regge::ReadParam({211, -211}, LoadedPDGTable(), *gra::MModelTune::Load(path.second));
  CHECK_FALSE(gra::regge::Photo(param, 443).dlog.has_value());
  CHECK(gra::regge::Photo(param, 443).diss_dlog.empty());
}

// Integrate every configured vector transition and compare to the analytic normalized HERA density
TEST_CASE("HERA proton dissociation has normalized MP XP GP cross sections",
          "[gra::MRegge][photoproduction][dissociation][physics]") {
  const std::array<std::pair<std::string, std::string>, 7> channels = {{{"rho_770", "pi+ pi-"}, {"phi_1020", "K+ K-"},
      {"jpsi", "mu+ mu-"}, {"psi2S", "mu+ mu-"}, {"Y1S", "mu+ mu-"}, {"Y2S", "mu+ mu-"}, {"Y3S", "mu+ mu-"}}};
  const auto tune = gra::MModelTune::Load(modelfile);
  for (const auto model : {gra::ReggeProductionModel::MP, gra::ReggeProductionModel::XP, gra::ReggeProductionModel::GP}) {
    for (const auto& [card, decay] : channels) {
      CAPTURE(gra::ReggeProductionModelName(model), card);
      ToyHelicityProcess process;
      ConfigureToyProductionProcess(process, gra::ReggeProductionModelName(model), "RES", decay);
      MRandom rng;
      process.SetResonances({{card, gra::resonance::Read("RES/" + card + ".json", rng, model)}});
      process.InitializeProcessAmplitude();
      const auto& res = process.GetResonances().at(card);
      auto lts = MakeToyCoherentPhotonLTS();
      lts.process.MMAX = process.state.lts.process.MMAX;
      gra::MRegge regge(lts, tune, gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "HERA dissociation"));
      const auto param = gra::regge::ReadParam({211, -211}, LoadedPDGTable(), *tune);
      const auto& photo = gra::regge::Photo(param, res.p.pdg);
      REQUIRE(photo.diss.has_value());
      const auto& diss = *photo.diss;
      const double W = diss.W0;
      const double residue = PhotoResidueNorm(regge, lts, res);
      const double norm = residue / pow2(W * W);
      const double evolution = photo.dlog ? std::exp(2 * PhotoDLogIntegral(*photo.dlog, W * W, pow2(photo.W0))) : 1.0;
      const double forward = residue * diss.ratio * std::pow(W / photo.W0, 4.0 * (photo.a0 - 1.0)) * evolution;
      const double mass_min2 = pow2(gra::PDG::mp + gra::PDG::mpi) * (1.0 + 1.0e-12);
      const double log_span = std::log(pow2(diss.mass_max) / mass_min2);
      const auto dsdt = [&](const double t) {
        std::vector<double> density(257);
        for (const auto& i : indices(density)) {
          const double mass2 = mass_min2 * std::exp(log_span * i / 256.0);
          density[i] = norm * std::norm(regge.PhotoDissKernel(W * W, t, mass2, res.p.pdg, gra::DissociationType::Hera)) * mass2;
        }
        return gra::math::CS13Integral(density, log_span / 256.0);
      };
      constexpr double tmax = 1.3;
      std::vector<double> spectrum(513);
      for (const auto& i : indices(spectrum)) { spectrum[i] = dsdt(-tmax * i / 512.0); }
      const double integral = diss.n / (diss.b * (diss.n - 1.0)) *
          (1.0 - std::pow(1.0 + diss.b * tmax / diss.n, 1.0 - diss.n));
      CHECK(gra::math::CS13Integral(spectrum, tmax / 512.0) == Approx(forward * integral).epsilon(1.0e-6));
      CHECK(dsdt(-0.2) == Approx(forward * std::pow(1.0 + diss.b * 0.2 / diss.n, -diss.n)).epsilon(1.0e-6));
      const auto amp = regge.PhotoDissKernel(W * W, -0.2, 4.0, res.p.pdg, gra::DissociationType::Hera);
      const auto higher = regge.PhotoDissKernel(4.0 * W * W, -0.2, 4.0, res.p.pdg, gra::DissociationType::Hera);
      double log_evolution = 0.0, phase = 0.0;
      for (const auto &term : photo.diss_dlog) {
        log_evolution += PhotoDLogIntegral(term, 4.0 * W * W, W * W);
        phase -= gra::math::PI / 2.0 * term.c / 4.0 *
            (1.0 / std::sqrt(std::log(4.0 * W * W / term.scale2)) - 1.0 / term.root0);
      }
      const auto expected = std::polar(4.0 * std::pow(2.0, diss.delta / 2.0) * std::exp(log_evolution), phase);
      RequireComplexNear(higher / amp, expected, 1.0e-12);
      CHECK(std::abs(regge.PhotoDissKernel(W * W, 0.0, 4.0, res.p.pdg, gra::DissociationType::Hera)) > 0.0);
      CHECK(std::abs(regge.PhotoDissKernel(W * W, -0.2, pow2(gra::PDG::mp), res.p.pdg, gra::DissociationType::Hera)) == Approx(0.0));
      CHECK_THROWS_AS(regge.PhotoDissKernel(W * W, 0.2, 4.0, res.p.pdg, gra::DissociationType::Hera), gra::AmplitudeFailure);
    }
  }
}

// Check the HERA target density is independent of the triple Pomeron profile and refreshed for shifted transfers
TEST_CASE("HERA dissociation sources preserve beam symmetry and screening shifts",
          "[gra::MRegge][photoproduction][dissociation][screening]") {
  const int         pdg        = GENERATE(113, 333);
  const std::string transition = GENERATE("geometric", "arithmetic", "diagonal");
  auto              lts        = MakeToyCoherentPhotonLTS();
  lts.process.PHOTO_DISSOCIATION = gra::DissociationType::Hera;
  const auto        path       = WriteModifiedPhotoVMTune("hera_projection_" + transition, [&transition](auto& card) {
    card["PARAM_SOFT"]["active_model"] = "double";
    auto& model = card["PARAM_SOFT"]["MODEL"]["double"];
    model["GW"]["theta"]                  = {0.41};
    model["EXCHANGE"]["P"]["g"]            = {{4.0, 1.25}, {1.25, 2.0}};
    model["EXCHANGE"]["P"]["transition_ff"] = transition;
  });
  const auto tune = gra::MModelTune::Load(path.second);
  gra::MRegge regge(lts, tune, gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "HERA source"));
  const auto   param     = gra::regge::ReadParam({211, -211}, LoadedPDGTable(), *tune);
  const auto&  soft      = *regge.SoftModelHandle();
  const auto   exchange  = param.exchanges[param.pomeron_trajectory].soft_exchange;
  const double reference = pow2(soft.PhysicalResidue(exchange, 0.0));
  for (const auto leg : {gra::ForwardBeamLeg::Upper, gra::ForwardBeamLeg::Lower}) {
    auto  upper        = gra::ResolveForwardLegState(lts, gra::ForwardBeamLeg::Upper);
    auto  lower        = gra::ResolveForwardLegState(lts, gra::ForwardBeamLeg::Lower);
    auto& target       = leg == gra::ForwardBeamLeg::Upper ? upper : lower;
    target.final_state = gra::ForwardFinalState::InclusiveExcitation;
    target.mass2       = 4.0;
    for (const double t : {-0.1, -0.4}) {
      target.t = t;
      gra::ReggeSourceCache cache(lts, regge, param, upper, lower);
      const auto&           source = cache.Sources(leg, 990, 1600.0, true, pdg);
      const auto            direct = gra::LegSources(lts, target, regge, param, 990, 1600.0, true, pdg);
      REQUIRE(source.size() == 1);
      CHECK(source.front().sector == gra::ProtonGoodWalkerSector::PhotoDiss);
      CHECK(gra::SquaredNorm(source.front().flip) == Approx(0.0));
      const double norm = gra::SquaredNorm(source.front().nonflip);
      CHECK(norm == Approx(gra::SquaredNorm(direct.front().nonflip)).epsilon(1.0e-12));
      const double reduced = norm / std::norm(regge.PhotoDissKernel(1600.0, t, 4.0, pdg, gra::DissociationType::Hera));
      CHECK(reduced == Approx(reference).epsilon(1.0e-12));
    }
  }
}

// Exercise initialized pp SD and DD and ep SD amplitudes with either physical photon source
TEST_CASE("MP XP GP pp and ep proton dissociation preserves rotations and beam exchange",
          "[gra::MRegge][photoproduction][dissociation][symmetry]") {
  const auto [lepton, nstars] = GENERATE(std::make_pair(false, 1), std::make_pair(false, 2), std::make_pair(true, 1));
  const std::string card      = GENERATE("jpsi", "phi_1020");
  const std::string prescription = GENERATE("soft", "hera");
  MRandom           rng;
  const bool        phi       = card == "phi_1020";
  for (const std::string model : {"MP", "XP", "GP"}) {
    const auto selected = WriteModifiedPhotoVMTune("diss_symmetry_" + model + prescription, [&](auto &input) {
      input["PARAM_NSTAR"]["MODEL"][model] = prescription;
    });
    const auto resonance = gra::resonance::Read("RES/" + card + ".json", rng, gra::ParseReggeProductionModel(model));
    for (const std::string photon_vertex : {"EPA", "QED"}) {
      CAPTURE(model, prescription, photon_vertex, lepton, nstars, card);
      ToyHelicityProcess process;
      ConfigureToyProductionProcess(process, model, "RES", phi ? "e+ e-" : "mu+ mu-");
      auto lts = MovingPhotoPairEvent(phi ? 11 : 13, resonance.p.mass);
      if (lepton) {
        lts.beam1 = lts.PDG.FindByPDG(11);
        lts.pbeam1.SetPxPyPzM(0.0, 0.0, lts.pbeam1.Pz(), lts.beam1.mass);
        lts.pfinal[1].SetPxPyPzM(lts.pfinal[1].Px(), lts.pfinal[1].Py(), lts.pfinal[1].Pz(), lts.beam1.mass);
      }
      SetToyPhotoForwardExcitation(lts, 2, 2.0);
      if (nstars == 2) { SetToyPhotoForwardExcitation(lts, 1, 2.0); }
      SetPhotoDecay(lts, 0.23, 0.41);
      RefreshToyDerivedKinematicsPreserveDecay(lts);
      process.state.lts = lts;
      process.SetModelTune(gra::MModelTune::Load(selected.second));
      process.SetExcitation(nstars);
      process.SetResonances({{card, resonance}});
      REQUIRE_NOTHROW(process.PrepareRun());
      process.state.lts.process.PHOTON_VERTEX = photon_vertex;
      auto         original                   = process.state.lts;
      const double norm                       = process.ProcPtr.GetBareAmplitude2(original);
      REQUIRE(std::isfinite(norm));
      REQUIRE(norm > 0.0);
      CHECK(original.excite1 == (nstars == 2));
      CHECK(original.excite2);
      auto rotated = RotateToyEventAroundZ(process.state.lts, 0.57);
      CHECK(process.ProcPtr.GetBareAmplitude2(rotated) == Approx(norm).epsilon(1.0e-7));
      auto reflected = ReflectToyEventInXZ(process.state.lts);
      CHECK(process.ProcPtr.GetBareAmplitude2(reflected) == Approx(norm).epsilon(1.0e-7));
      auto exchanged = ExchangePhotoBeams(process.state.lts);
      RefreshToyDerivedKinematicsPreserveDecay(exchanged);
      CHECK(process.ProcPtr.GetBareAmplitude2(exchanged) == Approx(norm).epsilon(1.0e-7));
      CHECK(exchanged.excite1);
      CHECK(exchanged.excite2 == (nstars == 2));
    }
  }
}

// Check the configured dissociation density acts on the proton in ep and pe amplitudes
TEST_CASE("Vector dissociation follows the proton leg in ep and pe",
          "[gra::MRegge][gra::MPhotoVM][photoproduction][dissociation][symmetry][regression]") {
  const std::string model = GENERATE("MP", "XP", "GP", "ygg");
  const int lepton = GENERATE(11, -11);
  const std::string prescription = GENERATE("soft", "hera");
  const auto selected = WriteModifiedPhotoVMTune("photo_diss_" + model + prescription, [&](auto& card) {
    card["PARAM_NSTAR"]["MODEL"][model] = {{"[*]", "soft"}, {"[22]", "structure"}, {"[22,P]", prescription}};
  });
  const auto tune = gra::MModelTune::Load(selected.second);
  const auto changed = WriteModifiedPhotoVMTune("photo_diss_varied_" + model + prescription, [&](auto& card) {
    card["PARAM_NSTAR"]["MODEL"][model] = {{"[*]", "soft"}, {"[22]", "structure"}, {"[22,P]", prescription}};
    auto& soft = card["PARAM_SOFT"]["FORWARD_EXCITATION"];
    soft["s0"] = soft["s0"].template get<double>() * 1.3;
    soft["a"] = soft["a"].template get<double>() * 1.7;
    for (auto& row : card["PARAM_REGGE"]["photoprod_diss"]) {
      row[2] = row[2].template get<double>() * 1.3;
      row[3] = row[3].template get<double>() + 0.2;
      row[4] = row[4].template get<double>() * 1.7;
    }
  });
  const auto varied = gra::MModelTune::Load(changed.second);
  const auto density = gra::flux::ReadPhotoDiss(tune->General()["PARAM_REGGE"]["photoprod_diss"]).at(443);
  const auto shifted = gra::flux::ReadPhotoDiss(varied->General()["PARAM_REGGE"]["photoprod_diss"]).at(443);
  auto event = MovingPhotoPairEvent(13, LoadedPDGTable().FindByPDG(443).mass);
  event.beam1 = event.PDG.FindByPDG(lepton);
  event.pbeam1.SetPxPyPzM(0.0, 0.0, event.pbeam1.Pz(), event.beam1.mass);
  event.pfinal[1].SetPxPyPzM(event.pfinal[1].Px(), event.pfinal[1].Py(), event.pfinal[1].Pz(), event.beam1.mass);
  SetToyPhotoForwardExcitation(event, 2, 2.0);
  SetPhotoDecay(event, 0.23, 0.41);
  RefreshToyDerivedKinematicsPreserveDecay(event);
  const auto amplitude = [&](gra::LORENTZSCALAR lts, gra::MModelTunePtr parameters) {
    ToyHelicityProcess process;
    ConfigureToyProductionProcess(process, model, model == "ygg" ? "jpsi" : "RES", "mu+ mu-");
    process.state.lts = lts;
    process.state.lts.model_cache = std::make_shared<gra::MModelCache>(parameters);
    process.SetModelTune(parameters);
    process.SetExcitation(1);
    if (model != "ygg") {
      process.SetResonances({{"jpsi", gra::resonance::Read("RES/jpsi.json", process.state.random,
                                                         gra::ParseReggeProductionModel(model))}});
    }
    process.PrepareRun();
    REQUIRE(process.ProcPtr.GetBareAmplitude2(process.state.lts) > 0.0);
    return process.state.lts.hamp;
  };
  double reference = 0.0;
  for (const bool reverse : {false, true}) {
    auto lts = reverse ? ExchangePhotoBeams(event) : event;
    RefreshToyDerivedKinematicsPreserveDecay(lts);
    lts.model_cache = std::make_shared<gra::MModelCache>(tune);
    const auto target = gra::ResolveForwardLegState(lts, reverse ? gra::ForwardBeamLeg::Upper : gra::ForwardBeamLeg::Lower);
    const double w2 = ((reverse ? lts.q2 : lts.q1) + target.incoming).M2();
    const auto exchange = tune->Soft()->ForwardExcitationExchange();
    const double factor = prescription == "soft"
        ? varied->Soft()->ForwardExcitationFactor(exchange, target.t, target.mass2) /
          tune->Soft()->ForwardExcitationFactor(exchange, target.t, target.mass2)
        : gra::flux::PhotoDissFactor(shifted, w2, target.t, target.mass2) /
          gra::flux::PhotoDissFactor(density, w2, target.t, target.mass2);
    CAPTURE(model, prescription, lepton, reverse, target.t, target.mass2);
    auto expected = amplitude(lts, tune);
    const double norm = gra::SquaredNorm(expected);
    if (!reverse) { reference = norm; }
    CHECK(norm == Approx(reference).epsilon(1.0e-7));
    std::complex<double> weight = factor;
    if (prescription == "hera" && model != "ygg") {
      const auto param = gra::regge::ReadParam({211, -211}, LoadedPDGTable(), *tune);
      const auto signature = tune->Soft()->Exchange(exchange).Signature();
      const double alpha = 1.0 + density.delta / 4.0;
      const double shifted_alpha = 1.0 + shifted.delta / 4.0;
      weight *= gra::regge::EtaFactor(shifted_alpha, shifted_alpha, signature, param.photoprod_eta_mode) /
          gra::regge::EtaFactor(alpha, alpha, signature, param.photoprod_eta_mode);
    }
    gra::Scale(expected, weight);
    RequireVectorNear(amplitude(lts, varied), expected, 1.0e-7);
  }
}

// Compare the full ygg amplitudes with the independently evaluated elastic forward normalization
TEST_CASE("HERA ygg dissociation uses one forward normalization at W0", "[gra::MPhotoVM][dissociation][normalization]") {
  const auto soft_path = WriteModifiedPhotoVMTune("diss_anchor_soft", [](auto &card) {
    card["PARAM_NSTAR"]["MODEL"]["ygg"] = "soft";
  });
  const auto hera_path = WriteModifiedPhotoVMTune("diss_anchor_hera", [](auto &card) {
    card["PARAM_NSTAR"]["MODEL"]["ygg"] = {{"[*]", "soft"}, {"[22,P]", "hera"}};
  });
  const auto soft = gra::MModelTune::Load(soft_path.second);
  const auto hera = gra::MModelTune::Load(hera_path.second);
  const auto density = gra::flux::ReadPhotoDiss(hera->General()["PARAM_REGGE"]["photoprod_diss"]).at(443);
  for (const double scale : {0.8, 1.0, 1.2}) {
    auto first = MovingPhotoPairEvent(13, scale * LoadedPDGTable().FindByPDG(443).mass);
    first.beam1 = first.PDG.FindByPDG(11);
    first.pbeam1.SetPxPyPzM(0.0, 0.0, first.pbeam1.Pz(), first.beam1.mass);
    first.pfinal[1].SetPxPyPzM(first.pfinal[1].Px(), first.pfinal[1].Py(), first.pfinal[1].Pz(), first.beam1.mass);
    RefreshToyDerivedKinematicsPreserveDecay(first);
    auto second = first;
    first.model_cache = std::make_shared<gra::MModelCache>(soft);
    second.model_cache = std::make_shared<gra::MModelCache>(hera);
    gra::MPhotoVM a(first, soft, "jpsi", gra::MPhotoVM::ProcessDefinitionFor("jpsi"));
    gra::MPhotoVM b(second, hera, "jpsi", gra::MPhotoVM::ProcessDefinitionFor("jpsi"));
    REQUIRE(a.Amp2(first) > 0.0);
    REQUIRE(b.Amp2(second) > 0.0);
    RequireVectorNear(first.hamp, second.hamp, 1.0e-12);
    for (auto *event : {&first, &second}) {
      SetToyPhotoForwardExcitation(*event, 2, 2.0);
      SetPhotoDecay(*event, 0.23, 0.41);
      RefreshToyDerivedKinematicsPreserveDecay(*event);
    }
    const auto target = gra::ResolveForwardLegState(first, gra::ForwardBeamLeg::Lower);
    const double w2 = (first.q1 + target.incoming).M2();
    const double factor = gra::flux::PhotoDissFactor(density, w2, target.t, target.mass2) /
        soft->Soft()->ForwardExcitationFactor(soft->Soft()->ForwardExcitationExchange(), target.t, target.mass2);
    const double ratio = pow2(factor) * a.GammaPDSigmaDt(first, density.W0, 0.0) /
        a.GammaPDSigmaDt(first, std::sqrt(w2), 0.0);
    CHECK(b.Amp2(second) / a.Amp2(first) == Approx(ratio).epsilon(1.0e-10));
    auto rotated = RotateToyEventAroundZ(second, 0.57);
    CHECK(b.Amp2(rotated) == Approx(b.Amp2(second)).epsilon(1.0e-10));
  }
}

// Propagate one proton transition across different elastic energy laws at complex amplitude level
TEST_CASE("Double-log dissociation preserves the common proton transition", "[photoproduction][dissociation][phase]") {
  constexpr double w0 = 120.0, wd = 85.0, scale2 = 9.0, coefficient = 2.0;
  const double delta = 0.4 + 4.0 * (1.23 - 1.19) - coefficient *
      (1.0 / std::sqrt(std::log(wd * wd / scale2)) - 1.0 / std::sqrt(std::log(w0 * w0 / scale2)));
  const auto selected = WriteModifiedPhotoVMTune("diss_double_log", [&](auto &card) {
    auto &regge = card["PARAM_REGGE"];
    regge["photoprod"] = {{443, w0, 1.4, 1.19, 0.09}, {553, w0, 1.7, 1.23, 0.08}};
    regge["photoprod_dlog"] = {{443, scale2, coefficient}};
    regge["photoprod_diss"] = {{443, wd, 0.3, 0.4, 2.0, 4.0, 0.07, 10.0}, {553, wd, 0.3, delta, 2.0, 4.0, 0.07, 10.0}};
    regge["photoprod_diss_dlog"] = {{553, scale2, -coefficient}};
    regge["photoprod_eta_mode"] = "rotating_t0";
  });
  auto lts = MakeToyCoherentPhotonLTS();
  const auto tune = gra::MModelTune::Load(selected.second);
  gra::MRegge regge(lts, tune, gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "double-log dissociation"));
  for (const double w : {0.4 * wd, wd, 3.0 * wd, 12.0 * wd}) {
    const auto transition = [&](int pdg) {
      return regge.PhotoDissKernel(w * w, -0.2, 4.0, pdg, gra::DissociationType::Hera) / regge.PhotoKernel(w * w, 0.0, pdg);
    };
    const auto ratio = transition(553) / transition(443);
    CHECK(ratio.real() == Approx(1.0).epsilon(1e-12));
    CHECK(ratio.imag() == Approx(0.0).margin(1e-12));
  }
}

// Reject malformed HERA input and require a measured or explicitly supplied channel
TEST_CASE("HERA dissociation steering fails during input parsing",
          "[gra::MRegge][photoproduction][dissociation][input]") {
  for (const std::size_t column : {1U, 2U, 4U, 5U, 6U, 7U}) {
    const auto path = WriteModifiedPhotoVMTune("hera_diss_invalid_" + std::to_string(column), [column](auto& card) {
      card["PARAM_REGGE"]["photoprod_diss"][0][column] = -1.0;
    });
    const auto tune = gra::MModelTune::Load(path.second);
    CHECK_THROWS_AS(gra::regge::ReadParam({211, -211}, LoadedPDGTable(), *tune), std::invalid_argument);
  }
  const auto path = WriteModifiedPhotoVMTune(
      "hera_diss_missing", [](auto& card) {
        card["PARAM_REGGE"]["photoprod_diss"] = nlohmann::json::array();
        card["PARAM_REGGE"]["photoprod_diss_dlog"] = nlohmann::json::array();
        for (const auto model : {"MP", "XP", "GP", "ygg"}) { card["PARAM_NSTAR"]["MODEL"][model] = "hera"; }
      });
  for (const std::string model : {"MP", "XP", "GP", "ygg"}) {
    ToyHelicityProcess process;
    ConfigureToyProductionProcess(process, model, model == "ygg" ? "jpsi" : "RES", "mu+ mu-");
    process.SetModelTune(gra::MModelTune::Load(path.second));
    process.SetExcitation(1);
    if (model != "ygg") {
      process.SetResonances({{"jpsi", gra::resonance::Read("RES/jpsi.json", process.state.random,
                                                        gra::ParseReggeProductionModel(model))}});
    }
    CHECK_THROWS_WITH(process.PrepareRun(), Catch::Contains("photoprod_diss"));
  }
}

// Check shared HERA factors act once on the TP target, including double dissociation
TEST_CASE("Tensor HERA targets preserve complex normalization and beam symmetry",
          "[MTensorPomeron][photoproduction][dissociation][symmetry][regression]") {
  const std::string name = GENERATE("rho_770", "phi_1020");
  const bool leptonic = GENERATE(true, false);
  const bool noflip = GENERATE(true, false);
  const bool phi = name == "phi_1020";
  const int pdg = phi ? 333 : 113;
  const auto make_tune = [&](bool varied) {
    return WriteModifiedPhotoVMTune("tensor_hera_" + name + std::to_string(varied) + std::to_string(noflip), [&](auto &card) {
      card["PARAM_NSTAR"]["MODEL"]["TP"] = {{"[*]", "soft"}, {"[22]", "structure"}, {"[22,P]", "hera"}};
      card["PARAM_TENSORPOM"]["FORWARD_NOFLIP"] = noflip;
      if (varied) {
        for (auto &row : card["PARAM_REGGE"]["photoprod_diss"]) {
          if (row[0].template get<int>() == pdg) { row[2] = row[2].template get<double>() * 1.44; }
        }
        card["PARAM_SOFT"]["FORWARD_EXCITATION"]["s0"] = 2.3;
      }
    });
  };
  const auto tune = gra::MModelTune::Load(make_tune(false).second);
  const auto changed = gra::MModelTune::Load(make_tune(true).second);
  auto event = MovingPhotoPairEvent(phi ? 321 : 211, LoadedPDGTable().FindByPDG(pdg).mass);
  if (leptonic) {
    event.beam1 = event.PDG.FindByPDG(11);
    event.pbeam1.SetPxPyPzM(0.0, 0.0, event.pbeam1.Pz(), event.beam1.mass);
    event.pfinal[1].SetPxPyPzM(event.pfinal[1].Px(), event.pfinal[1].Py(), event.pfinal[1].Pz(), event.beam1.mass);
  } else {
    SetToyPhotoForwardExcitation(event, 1, 2.0);
  }
  SetToyPhotoForwardExcitation(event, 2, 2.0);
  SetPhotoDecay(event, 0.23, 0.41);
  RefreshToyDerivedKinematicsPreserveDecay(event);
  const auto amplitude = [&](gra::LORENTZSCALAR lts, gra::MModelTunePtr model, bool transfer) {
    ToyHelicityProcess process;
    ConfigureToyProductionProcess(process, "TP", "RES", phi ? "K+ K-" : "pi+ pi-");
    process.state.lts = lts;
    process.state.lts.model_cache = std::make_shared<gra::MModelCache>(model);
    process.SetModelTune(model);
    process.SetExcitation(leptonic ? 1 : 2);
    auto res = gra::resonance::Read("RES/" + name + ".json", process.state.random, gra::ReggeProductionModel::TP);
    if (transfer) {
      for (auto &channel : res.TP.channels) {
        channel.ff_transfer = gra::regge::FFParam{};
        channel.ff_transfer.type = gra::regge::FFType::Exponential;
        channel.ff_transfer.param = {13.0};
      }
    }
    process.SetResonances({{name, res}});
    process.PrepareRun();
    REQUIRE(process.state.lts.process.DISSOCIATION == gra::DissociationType::Soft);
    REQUIRE(process.state.lts.process.PHOTO_DISSOCIATION == gra::DissociationType::Hera);
    REQUIRE(process.ProcPtr.GetBareAmplitude2(process.state.lts) > 0.0);
    return process.state.lts.hamp;
  };
  const auto expected = amplitude(event, tune, false);
  const auto varied = amplitude(event, changed, true);
  REQUIRE(expected.size() == varied.size());
  for (const auto &i : indices(expected)) { RequireComplexNear(varied[i], 1.2 * expected[i], 3.0e-9); }
  const double norm = gra::SquaredNorm(expected);
  const std::array transformations = {RotateToyEventAroundZ(event, 0.57), ReflectToyEventInXZ(event),
                                      ExchangePhotoBeams(event), BoostPhotoEvent(event, 0.37)};
  for (const auto i : indices(transformations)) {
    auto transformed = transformations[i];
    RefreshToyDerivedKinematicsPreserveDecay(transformed);
    // Preserve light-cone fractions under a longitudinal boost
    transformed.xi1 = transformations[i].xi1;
    transformed.xi2 = transformations[i].xi2;
    CAPTURE(name, leptonic, noflip, i);
    CHECK(gra::SquaredNorm(amplitude(transformed, tune, false)) == Approx(norm).epsilon(3.0e-7));
  }
}

// Compare the TP HERA tensor current to the elastic forward current with the same complex phase
TEST_CASE("Tensor HERA current uses the common W0 anchor and density",
          "[MTensorPomeron][photoproduction][dissociation][normalization]") {
  auto lts = MovingPhotoPairEvent(211, LoadedPDGTable().FindByPDG(113).mass);
  SetToyPhotoForwardExcitation(lts, 2, 2.0);
  const auto tune = gra::MModelTune::Load(modelfile);
  gra::MTensorPomeron tensor(lts, tune, gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Resonance));
  const auto parameters = gra::GetTensorParam(*lts.model_cache, lts.PDG, lts.process.RESONANCES);
  const auto &diss = parameters->photo_diss.at(113);
  auto target = gra::ResolveForwardLegState(lts, gra::ForwardBeamLeg::Lower);
  auto forward = target;
  forward.final_state = gra::ForwardFinalState::Elastic;
  forward.t = 0.0;
  const auto vertex = tensor.iG_TForwardHE(forward, 995);
  const auto anchor = tensor.PomeronPropagatorCurrent(vertex, tensor.TensorPropagatorFactor(995, pow2(diss.W0), 0.0));
  for (const double scale : {0.7, 1.0, 1.4}) {
    const double w2 = pow2(scale * diss.W0);
    for (const double t : {-0.02, -0.4, -2.0}) {
      target.t = t;
      const auto value = tensor.PhotoDissCurrent(target, 995, 113, w2, w2);
      const double factor = pow2(diss.W0) / w2 * gra::flux::PhotoDissFactor(diss, w2, t, target.mass2);
      for (const auto &mu : tensor.LI) {
        for (const auto &nu : tensor.LI) { RequireComplexNear(value(mu, nu), factor * anchor(mu, nu), 1.0e-12); }
      }
    }
  }
}

// Require a measured vector row and reject a continuum interpretation of the HERA fit
TEST_CASE("Tensor HERA selection validates its physical channel before sampling",
          "[MTensorPomeron][photoproduction][dissociation][input]") {
  const auto missing = WriteModifiedPhotoVMTune("tensor_hera_missing", [](auto &card) {
    card["PARAM_NSTAR"]["MODEL"]["TP"] = {{"[*]", "soft"}, {"[22,P]", "hera"}};
    card["PARAM_REGGE"]["photoprod_diss"] = nlohmann::json::array();
    card["PARAM_REGGE"]["photoprod_diss_dlog"] = nlohmann::json::array();
  });
  ToyHelicityProcess vector;
  ConfigureToyProductionProcess(vector, "TP", "RES", "pi+ pi-");
  vector.state.lts = MovingPhotoPairEvent(211, LoadedPDGTable().FindByPDG(113).mass);
  vector.SetModelTune(gra::MModelTune::Load(missing.second));
  vector.SetExcitation(1);
  vector.SetResonances({{"rho_770", gra::resonance::Read("RES/rho_770.json", vector.state.random, gra::ReggeProductionModel::TP)}});
  CHECK_THROWS_WITH(vector.PrepareRun(), Catch::Contains("photoprod_diss"));
  ToyHelicityProcess continuum;
  ConfigureToyProductionProcess(continuum, "TP", "PHOTO", "pi+ pi-");
  continuum.SetExcitation(1);
  CHECK_THROWS_WITH(continuum.PrepareRun(), Catch::Contains("continuum"));
}

// A photon-target override must not modify double-Pomeron scalar production
TEST_CASE("Photoproduction HERA overrides leave hadronic resonance amplitudes unchanged",
          "[photoproduction][dissociation][normalization]") {
  for (const std::string model : {"MP", "XP", "GP", "TP"}) {
    const auto point = MovingPhotoPairEvent(211, LoadedPDGTable().FindByPDG(9010221).mass);
    const auto amplitude = [&](const std::string &photo) {
      const auto path = WriteModifiedPhotoVMTune("diss_hadronic_rule_" + model + photo, [&](auto &card) {
        card["PARAM_NSTAR"]["MODEL"][model] = {{"[*]", "soft"}, {"[22,P]", photo}};
      });
      const auto tune = gra::MModelTune::Load(path.second);
      ToyHelicityProcess process;
      ConfigureToyProductionProcess(process, model, "RES", "pi+ pi-");
      process.state.lts = point;
      process.state.lts.model_cache = std::make_shared<gra::MModelCache>(tune);
      process.SetModelTune(tune);
      process.SetExcitation(2);
      SetToyPhotoForwardExcitation(process.state.lts, 1, 2.0);
      SetToyPhotoForwardExcitation(process.state.lts, 2, 2.0);
      SetPhotoDecay(process.state.lts, 0.23, 0.41);
      RefreshToyDerivedKinematicsPreserveDecay(process.state.lts);
      process.SetResonances({{"f0_980", gra::resonance::Read("RES/f0_980.json", process.state.random,
                                                         gra::ParseReggeProductionModel(model))}});
      process.PrepareRun();
      REQUIRE(process.ProcPtr.GetBareAmplitude2(process.state.lts) > 0.0);
      return process.state.lts.hamp;
    };
    CAPTURE(model);
    RequireVectorNear(amplitude("hera"), amplitude("soft"), 1.0e-12);
  }
}

// HERA fits apply only to vector photoproduction, while soft excitation also covers hadronic resonances
TEST_CASE("HERA target steering rejects hadronic production before sampling", "[photoproduction][dissociation][input]") {
  for (const std::string prescription : {"soft", "hera"}) {
    const auto selected = WriteModifiedPhotoVMTune("diss_hadronic_" + prescription, [&](auto &card) {
      card["PARAM_NSTAR"]["MODEL"]["MP"] = prescription;
    });
    ToyHelicityProcess process;
    ConfigureToyProductionProcess(process, "MP", "RES", "pi+ pi-");
    process.state.lts = MovingPhotoPairEvent(211, LoadedPDGTable().FindByPDG(9010221).mass);
    process.SetModelTune(gra::MModelTune::Load(selected.second));
    process.SetExcitation(1);
    process.SetResonances({{"f0_980", gra::resonance::Read("RES/f0_980.json", process.state.random, gra::ReggeProductionModel::MP)}});
    if (prescription == "soft") { CHECK_NOTHROW(process.PrepareRun()); }
    else { CHECK_THROWS_WITH(process.PrepareRun(), Catch::Contains("hera requires vector")); }
  }
}

TEST_CASE("gra::spin:: photoproduction branches use photon source not proton residue", "[gra::spin]") {
  gra::LORENTZSCALAR lts                  = MakeToyProductionLTSAsymmetric(0.3, 4.2, -4.6);
  gra::PARAM_RES     res                  = MakeToyPhotoMPResonance(false);
  res.production[0].tree[0].hel           = RealisticProtonLegHelicityMatrix(22, 2, -1, 1, 2);
  res.production[0].tree[0].hel.Jz_values = {-1.0, 1.0};
  res.production.front().pole = gra::spin::PreparePoleLS(
      res.p, res.production[0].tree[0].p, res.production[0].tree[1].p, {{0, 2, 1.0}}, 1.0, true, false, true);
  res.production.front().hel = res.production.front().pole->helicity;

  const auto up_rows        = HelicityConservingRows(res.production[0].tree[0].hel);
  const auto dn_rows        = HelicityConservingRows(res.production[0].tree[1].hel);

  const auto f_up_photon = SelectRows(PhotonLegMatrixForTest(lts, res.production[0].tree[0].hel, 1, false), up_rows);
  const auto f_dn_regge  = NormalizeForwardRowsForTest(
       SelectRows(gra::spin::Forward(lts, res.production[0].tree[1], lts.pbeam2, lts.pfinal[2], true,
                                     gra::spin::Rows(res.production[0].tree[1], false), lts.process.PHOTON_VERTEX,
                                     gra::spin::ForwardSpec{}),
                  dn_rows));
  const auto f_x      = gra::mpom::Fusion(lts, res.production.front().pole.value());
  const auto expected = f_up_photon.Kronecker(f_dn_regge) * f_x;
  const auto actual   = gra::rspin::Resonance(lts, res, (res.production_model == gra::ReggeProductionModel::MP ? gra::mpom::Fusion : gra::xpom::Fusion), 1.0);

  REQUIRE(actual.size() == 1);
  REQUIRE(gra::spin::SourceSpinAveragedDensity(f_up_photon, 2, "test upper photon") == Approx(1.0).epsilon(1e-13));
  REQUIRE(gra::spin::ForwardSourceColumnSpinAveragedDensity(f_dn_regge, f_dn_regge.size_col() / 2, 2,
                                                            "test lower Regge m0") == Approx(1.0).epsilon(1e-13));
  RequireMatrixNear(actual[0], expected);

  std::swap(res.production[0].tree[0], res.production[0].tree[1]);
  res.production.front().pole = gra::spin::PreparePoleLS(
      res.p, res.production[0].tree[0].p, res.production[0].tree[1].p, {{0, 2, 1.0}}, 1.0, true, false, true);
  res.production.front().hel = res.production.front().pole->helicity;
  const auto up_rows_lower_gamma = HelicityConservingRows(res.production[0].tree[0].hel);
  const auto dn_rows_lower_gamma = HelicityConservingRows(res.production[0].tree[1].hel);
  const auto f_up_regge          = NormalizeForwardRowsForTest(
               SelectRows(gra::spin::Forward(lts, res.production[0].tree[0], lts.pbeam1, lts.pfinal[1], false,
                                             gra::spin::Rows(res.production[0].tree[0], false), lts.process.PHOTON_VERTEX,
                                             gra::spin::ForwardSpec{}),
                          up_rows_lower_gamma));
  const auto f_dn_photon =
      SelectRows(PhotonLegMatrixForTest(lts, res.production[0].tree[1].hel, 2, true), dn_rows_lower_gamma);
  const auto f_x_lower = gra::mpom::Fusion(lts, res.production.front().pole.value());
  const auto expected_lower_gamma = f_up_regge.Kronecker(f_dn_photon) * f_x_lower;
  const auto actual_lower_gamma   = gra::rspin::Resonance(lts, res, (res.production_model == gra::ReggeProductionModel::MP ? gra::mpom::Fusion : gra::xpom::Fusion), 1.0);

  REQUIRE(actual_lower_gamma.size() == 1);
  REQUIRE(gra::spin::ForwardSourceColumnSpinAveragedDensity(f_up_regge, f_up_regge.size_col() / 2, 2,
                                                            "test upper Regge m0") == Approx(1.0).epsilon(1e-13));
  REQUIRE(gra::spin::SourceSpinAveragedDensity(f_dn_photon, 2, "test lower photon") == Approx(1.0).epsilon(1e-13));
  RequireMatrixNear(actual_lower_gamma[0], expected_lower_gamma);
}

TEST_CASE(
    "gra::spin:: FORWARD_VERTEX scalar-Pomeron full-helicity "
    "photoproduction works in both beam orderings",
    "[gra::spin]") {
  gra::LORENTZSCALAR lts     = MakeToyProductionLTSAsymmetric(0.3, 4.2, -4.6);
  gra::PARAM_RES     res     = MakeToyScalarPomeronPhotoXP();
  lts.process.FORWARD_VERTEX = gra::ForwardVertexMode::HelicityResidue;
  lts.process.FORWARD_NOFLIP = false;

  const gra::MDecayBranch gamma_branch   = res.production[0].tree[0];
  const gra::MDecayBranch pomeron_branch = res.production[0].tree[1];
  res.production.push_back({{pomeron_branch, gamma_branch}});
  PrepareToyXPOperators(res, {{{0, 2, 1.0}}, {{0, 2, 1.0}}}, true, true);

  const auto full = gra::rspin::Resonance(lts, res, (res.production_model == gra::ReggeProductionModel::MP ? gra::mpom::Fusion : gra::xpom::Fusion), 1.0);
  REQUIRE(full.size() == 2);
  for (const auto &channel_matrix : full) {
    CAPTURE(MatrixNorm2(channel_matrix));
    REQUIRE(channel_matrix.size_row() == 16);
    REQUIRE(channel_matrix.size_col() == 3);
    REQUIRE(MatrixNorm2(channel_matrix) > 1e-12);
  }
}

// Reject malformed quadrature counts before constructing integration rules
TEST_CASE("PhotoZ rejects fractional and overflowing integration counts", "[gra::MPhotoZ][params]") {
  const auto source = gra::ResolveModelDataFile("TUNE0", "NUMERICS.json");
  const auto document = nlohmann::json::parse(gra::aux::GetInputData(source));
  const nlohmann::json invalid = {-1, 0, 1.5, "4", 4294967297ULL, -4294967295LL};
  for (const auto &name : {"N_kappa", "N_k"}) {
    for (const auto &value : invalid) {
      auto input = document;
      input["NUMERICS_PHOTOZ"][name] = value;
      gra::MPhotoQCDNumerics numerics;
      INFO("count=" << name << " value=" << value);
      REQUIRE_THROWS_AS(numerics.Configure(input, source, "NUMERICS_PHOTOZ"), std::invalid_argument);
    }
  }
}

TEST_CASE("MPhotoQCD reads separate PhotoZ physics and numerical blocks", "[gra::MPhotoZ][params]") {
  ModelParamRestoreGuard restore;
  gra::MODELPARAM                 = "TUNE0";
  const std::string numerics_file = gra::ResolveModelDataFile("TUNE0", "NUMERICS.json");

  gra::MPhotoQCDParam default_param;
  REQUIRE_NOTHROW(default_param.ConfigureFromJson(modelfile, gra::aux::GetInputData(modelfile), "PARAM_PHOTOZ"));
  REQUIRE(default_param.initialized);

  gra::MPhotoQCDNumerics default_numerics;
  REQUIRE_NOTHROW(
      default_numerics.ConfigureFromJson(numerics_file, gra::aux::GetInputData(numerics_file), "NUMERICS_PHOTOZ"));
  REQUIRE(default_numerics.initialized);
  REQUIRE(default_numerics.kappa2_nodes.size() == default_numerics.N_kappa);
  REQUIRE(default_numerics.kappa2_weights.size() == default_numerics.N_kappa);
  REQUIRE(default_numerics.K2Max(91.1876 * 91.1876) == Approx(91.1876 * 91.1876));

  const double w2    = gra::math::pow2(190.0);
  const double t     = -0.35;
  const double slope = default_param.t_slope_B0 +
                       4.0 * default_param.t_slope_alpha_prime * std::log(std::sqrt(w2) / default_param.t_slope_W0);
  REQUIRE(gra::MPhotoQCD::XGluon(default_param, 100.0, w2) == Approx(default_param.xg_coefficient * 100.0 / w2));
  REQUIRE(gra::MPhotoQCD::HardScale(default_param, 1.0) == Approx(default_param.mu_MIN));
  REQUIRE(gra::MPhotoQCD::HardScale(default_param, 100.0) == Approx(default_param.mu_over_m * 10.0));
  REQUIRE(gra::MPhotoQCD::RealPartFactor(default_param).real() ==
          Approx(std::tan(0.5 * gra::math::PI * default_param.real_part_delta)));
  REQUIRE(gra::MPhotoQCD::RealPartFactor(default_param).imag() == Approx(1.0));
  REQUIRE(gra::MPhotoQCD::TSlope(default_param, w2) == Approx(slope));
  REQUIRE(gra::MPhotoQCD::TSlopeFactor(default_param, w2, t) == Approx(std::exp(0.5 * slope * t)).epsilon(1e-14));

  const auto steered = WriteModifiedPhotoZTune(
      "valid_steering", [](auto &j) { j["PARAM_PHOTOZ"]["real_part_delta"] = 0.12; },
      [](auto &j) {
        j["NUMERICS_PHOTOZ"]["N_kappa"]       = 3;
        j["NUMERICS_PHOTOZ"]["N_k"]           = 4;
        j["NUMERICS_PHOTOZ"]["k2_MAX_use_q2"] = false;
        j["NUMERICS_PHOTOZ"]["k2_MAX"]        = 10.0;
      });
  gra::MODELPARAM = steered.first;
  gra::MPhotoQCDParam steered_param;
  REQUIRE_NOTHROW(
      steered_param.ConfigureFromJson(steered.second, gra::aux::GetInputData(steered.second), "PARAM_PHOTOZ"));
  gra::MPhotoQCDNumerics steered_numerics;
  const std::string      steered_numerics_file =
      (std::filesystem::path(steered.second).parent_path() / "NUMERICS.json").string();
  REQUIRE_NOTHROW(steered_numerics.ConfigureFromJson(steered_numerics_file,
                                                     gra::aux::GetInputData(steered_numerics_file), "NUMERICS_PHOTOZ"));
  REQUIRE(steered_numerics.N_kappa == 3);
  REQUIRE(steered_numerics.N_k == 4);
  REQUIRE(steered_numerics.kappa2_nodes.size() == steered_numerics.N_kappa);
  REQUIRE(steered_numerics.kappa2_weights.size() == steered_numerics.N_kappa);
  REQUIRE_FALSE(steered_numerics.k2_MAX_use_q2);
  REQUIRE(steered_numerics.K2Max(91.1876 * 91.1876) == Approx(10.0));
  REQUIRE(steered_param.real_part_delta == Approx(0.12));

  auto require_invalid_numerics = [](const std::string &suffix, const std::function<void(nlohmann::json &)> &mutate,
                                     const std::string &message) {
    const auto tune = WriteModifiedPhotoZTune(
        "invalid_numerics_" + suffix, [](auto &) {}, mutate);
    gra::MODELPARAM             = tune.first;
    const std::string      path = (std::filesystem::path(tune.second).parent_path() / "NUMERICS.json").string();
    gra::MPhotoQCDNumerics numerics;
    REQUIRE_THROWS(numerics.ConfigureFromJson(path, gra::aux::GetInputData(path), "NUMERICS_PHOTOZ"));
  };

  require_invalid_numerics(
      "nkappa", [](auto &j) { j["NUMERICS_PHOTOZ"]["N_kappa"] = 0; }, "N_kappa");
  require_invalid_numerics(
      "kappa_range", [](auto &j) { j["NUMERICS_PHOTOZ"]["kappa2_MAX"] = 1000.0; }, "must not exceed NUMERICS_SUDAKOV");

  auto require_invalid_param = [](const std::string &suffix, const std::function<void(nlohmann::json &)> &mutate,
                                  const std::string &message) {
    const auto tune = WriteModifiedPhotoZTune("invalid_param_" + suffix, mutate, [](auto &) {});
    gra::MODELPARAM = tune.first;
    gra::MPhotoQCDParam param;
    REQUIRE_THROWS(param.ConfigureFromJson(tune.second, gra::aux::GetInputData(tune.second), "PARAM_PHOTOZ"));
  };

  require_invalid_param(
      "xg", [](auto &j) { j["PARAM_PHOTOZ"]["xg_coefficient"] = 0.0; }, "xg_coefficient");
  require_invalid_param(
      "mu", [](auto &j) { j["PARAM_PHOTOZ"]["mu_MIN"] = 0.0; }, "mu_MIN and mu_over_m");
  require_invalid_param(
      "mass", [](auto &j) { j["PARAM_PHOTOZ"]["light_quark_mass_MIN"] = -0.1; }, "light_quark_mass_MIN");
  require_invalid_param(
      "real_delta", [](auto &j) { j["PARAM_PHOTOZ"]["real_part_delta"] = 1.0; }, "real_part_delta");
  require_invalid_param(
      "slope", [](auto &j) { j["PARAM_PHOTOZ"]["t_slope_alpha_prime"] = -0.1; }, "t-slope parameters");
}

TEST_CASE("Photoproduction excitation keeps EPA and target factors separate", "[gra::MPhotoQCD][EPA][dissociation]") {
  ModelParamRestoreGuard restore;
  gra::MODELPARAM       = "TUNE0";
  const auto model_tune = gra::MModelTune::Load(modelfile);

  gra::LORENTZSCALAR emitter = MakeToyPhotoZFFbar(13);
  gra::RequireModelCache(emitter.model_cache, model_tune, "test emitter");
  SetToyPhotoForwardExcitation(emitter, 1, 2.0);
  const auto                 emitter_state  = gra::ResolveForwardLegState(emitter, gra::ForwardBeamLeg::Upper);
  const std::complex<double> emitter_factor = gra::qed::PhotonSourceAmplitude(emitter, emitter_state, "EPA");
  const double expected_flux = gra::flux::IncohFlux(emitter.xi1, emitter.t1, emitter.qt1, emitter.pfinal[1].M2());
  REQUIRE(std::norm(emitter_factor) * emitter.xi1 == Approx(expected_flux).epsilon(1.0e-12));
  REQUIRE_FALSE(emitter.proton_good_walker.has_value());

  gra::LORENTZSCALAR target = MakeToyPhotoZFFbar(13);
  gra::RequireModelCache(target.model_cache, model_tune, "test target");
  SetToyPhotoForwardExcitation(target, 2, 2.0);
  const auto   photo_param  = gra::ReadPhotoQCDParam(*model_tune, "PARAM_PHOTOZ");
  const double w2           = (target.q1 + target.pbeam2).M2();
  const auto   target_state = gra::ResolveForwardLegState(target, gra::ForwardBeamLeg::Lower);
  const double target_factor =
      gra::MPhotoQCD::TargetTransitionFactor(target_state, *photo_param, w2, *model_tune->Soft());
  const double expected_target = model_tune->Soft()->ForwardExcitationFactor(
      model_tune->Soft()->ForwardExcitationExchange(), target.q2.M2(), target.pfinal[2].M2());
  REQUIRE(target_factor == Approx(expected_target).epsilon(1.0e-12));
  REQUIRE_FALSE(target.proton_good_walker.has_value());
}

TEST_CASE("MPhotoQCD crosses the lower photon into the global boson section",
          "[gra::MPhotoQCD][EPA][helicity][phase]") {
  gra::LORENTZSCALAR lts = MakeToyPhotoZFFbar(13);
  gra::RequireModelCache(lts.model_cache, gra::MModelTune::Load(modelfile), "test photon source");
  constexpr double         rotation      = -0.64;
  const gra::LORENTZSCALAR rotated       = RotateToyEventAroundZ(lts, rotation);
  const auto               upper         = gra::ResolveForwardLegState(lts, gra::ForwardBeamLeg::Upper);
  const auto               lower         = gra::ResolveForwardLegState(lts, gra::ForwardBeamLeg::Lower);
  const auto               upper_rotated = gra::ResolveForwardLegState(rotated, gra::ForwardBeamLeg::Upper);
  const auto               lower_rotated = gra::ResolveForwardLegState(rotated, gra::ForwardBeamLeg::Lower);
  const std::vector<std::pair<double, double>> noflip                   = {{-0.5, -0.5}};
  const std::vector<int>                       local_photon_projections = {-1, 1};
  const auto                                   upper_exchange =
      gra::qed::PhotonSourceMatrixTransitions(lts, 1, noflip, local_photon_projections, "EPA", false);
  const auto lower_exchange =
      gra::qed::PhotonSourceMatrixTransitions(lts, 2, noflip, local_photon_projections, "EPA", true);

  for (const int projection : {-1, 1}) {
    for (const auto &[state, rotated_state] : std::array<std::pair<gra::ForwardLegState, gra::ForwardLegState>, 2>{
             std::pair{upper, upper_rotated}, std::pair{lower, lower_rotated}}) {
      const std::complex<double> expected =
          gra::qed::PhotonSourceAmplitude(lts, state, "EPA") *
          std::exp(gra::math::zi * static_cast<double>(projection) * state.transfer.Phi()) / std::sqrt(2.0);
      REQUIRE(std::abs(gra::qed::TransversePhotonSourceAmplitude(lts, state, projection) - expected) < 1.0e-14);
      const std::complex<double> reference = gra::qed::TransversePhotonSourceAmplitude(lts, state, projection);
      const std::complex<double> actual = gra::qed::TransversePhotonSourceAmplitude(rotated, rotated_state, projection);
      const std::complex<double> rotation_phase = std::exp(gra::math::zi * static_cast<double>(projection) * rotation);
      REQUIRE(std::abs(actual - rotation_phase * reference) < 2.0e-13);
    }

    const std::size_t          upper_local  = projection > 0 ? 1 : 0;
    const std::size_t          lower_local  = projection > 0 ? 0 : 1;
    const std::complex<double> upper_direct = gra::qed::TransversePhotonSourceAmplitude(lts, upper, projection);
    const std::complex<double> lower_direct = gra::qed::TransversePhotonSourceAmplitude(lts, lower, projection);
    // The direct dual circular basis and the JW exchange basis have different helicity phases
    REQUIRE(std::abs(upper_direct / std::abs(upper_direct) + static_cast<double>(projection) *
                                                                 upper_exchange[0][upper_local] /
                                                                 std::abs(upper_exchange[0][upper_local])) < 1.0e-13);
    REQUIRE(std::abs(lower_direct / std::abs(lower_direct) + static_cast<double>(projection) *
                                                                 lower_exchange[0][lower_local] /
                                                                 std::abs(lower_exchange[0][lower_local])) < 1.0e-13);
  }
}

TEST_CASE("Direct photoproduction preserves coherent beam-direction phases",
          "[gra::MPhotoZ][gra::MPhotoVM][helicity][phase][interference]") {
  ModelParamRestoreGuard restore;
  gra::MODELPARAM           = "TUNE0";
  constexpr double rotation = 0.57;

  const auto validate = [rotation](const std::string &label, const double mass, const auto &evaluate) {
    const gra::LORENTZSCALAR pp_event         = MovingPhotoPairEvent(13, mass);
    gra::LORENTZSCALAR       upper_antiproton = pp_event;
    gra::LORENTZSCALAR       lower_antiproton = pp_event;
    SetPhotoAntiprotonBeam(upper_antiproton, 1);
    SetPhotoAntiprotonBeam(lower_antiproton, 2);

    const auto pp             = evaluate(pp_event);
    const auto pbar_p         = evaluate(upper_antiproton);
    const auto p_pbar         = evaluate(lower_antiproton);
    const auto pp_rotated     = evaluate(RotateToyEventAroundZ(pp_event, rotation));
    const auto pbar_p_rotated = evaluate(RotateToyEventAroundZ(upper_antiproton, rotation));

    CAPTURE(label, pp.first.size(), pbar_p.first.size(), p_pbar.first.size());
    REQUIRE(pp.first.size() == 16);
    REQUIRE(pbar_p.first.size() == pp.first.size());
    REQUIRE(p_pbar.first.size() == pp.first.size());
    REQUIRE(pp_rotated.first.size() == pp.first.size());
    REQUIRE(pbar_p_rotated.first.size() == pp.first.size());
    REQUIRE(pp.second == Approx(gra::SquaredNorm(pp.first) / 4.0).epsilon(1.0e-12));
    REQUIRE(pp_rotated.second == Approx(pp.second).epsilon(2.0e-10));
    REQUIRE(evaluate(ExchangePhotoBeams(pp_event)).second == Approx(pp.second).epsilon(2.0e-9));

    std::vector<int> harmonics(pp.first.size(), 0);
    RequireAzimuthalCovariance(pp_rotated.first, pp.first, harmonics, rotation, 1.0e-9);

    double upper_norm2     = 0.0;
    double lower_norm2     = 0.0;
    double interference_l1 = 0.0;
    for (std::size_t index = 0; index < pp.first.size(); ++index) {
      const std::complex<double> upper                 = 0.5 * (pp.first[index] - pbar_p.first[index]);
      const std::complex<double> lower                 = 0.5 * (pp.first[index] + pbar_p.first[index]);
      const std::complex<double> upper_rotated         = 0.5 * (pp_rotated.first[index] - pbar_p_rotated.first[index]);
      const std::complex<double> lower_rotated         = 0.5 * (pp_rotated.first[index] + pbar_p_rotated.first[index]);
      // Only external fermion states remain after the internal spin sum
      const std::complex<double> phase = 1.0;
      const double scale =
          std::max({1.0, std::abs(upper), std::abs(lower), std::abs(upper_rotated), std::abs(lower_rotated)});

      CAPTURE(label, index, upper, lower, upper_rotated, lower_rotated);
      CHECK(std::abs(pp.first[index] - upper - lower) < 2.0e-12 * scale);
      CHECK(std::abs(p_pbar.first[index] + pbar_p.first[index]) < 2.0e-12 * scale);
      CHECK(std::abs(upper_rotated - phase * upper) < 3.0e-10 * scale);
      CHECK(std::abs(lower_rotated - phase * lower) < 3.0e-10 * scale);
      CHECK(std::abs(pp_rotated.first[index] - phase * pp.first[index]) < 3.0e-10 * scale);

      upper_norm2 += std::norm(upper);
      lower_norm2 += std::norm(lower);
      const double row_interference          = std::norm(pp.first[index]) - std::norm(upper) - std::norm(lower);
      const double opposite_row_interference = std::norm(pbar_p.first[index]) - std::norm(upper) - std::norm(lower);
      CHECK(std::abs(opposite_row_interference + row_interference) < 2.0e-11 * scale * scale);
      interference_l1 += std::abs(row_interference);
    }

    CAPTURE(label, upper_norm2, lower_norm2, interference_l1);
    REQUIRE(upper_norm2 > 0.0);
    REQUIRE(lower_norm2 > 0.0);
    REQUIRE(interference_l1 > 1.0e-10 * (upper_norm2 + lower_norm2));
  };

  const auto evaluate_photoz = [](gra::LORENTZSCALAR lts) {
    gra::MPhotoZ amplitude(lts, gra::MModelTune::Load(modelfile), gra::MPhotoZ::ProcessDefinitionFor());
    const double amp2 = amplitude.Amp2(lts);
    return std::make_pair(lts.hamp, amp2);
  };
  validate("PhotoZ", 91.1876, evaluate_photoz);
  const std::array<std::pair<const char *, int>, 5> vector_mesons = {{
      {"jpsi", 443},
      {"psi(2S)", 100443},
      {"Upsilon(1S)", 553},
      {"Upsilon(2S)", 100553},
      {"Upsilon(3S)", 200553},
  }};
  for (const auto &[channel, pdg] : vector_mesons) {
    const double mass             = LoadedPDGTable().FindByPDG(pdg).mass;
    const auto   evaluate_channel = [channel](gra::LORENTZSCALAR lts) {
      gra::MPhotoVM amplitude(lts, gra::MModelTune::Load(modelfile), channel,
                                gra::MPhotoVM::ProcessDefinitionFor(channel));
      const double  amp2 = amplitude.Amp2(lts);
      return std::make_pair(lts.hamp, amp2);
    };
    validate(std::string("PhotoVM ") + channel, mass, evaluate_channel);
  }
}

TEST_CASE("Direct photoproduction retains polarized decays and boost invariance",
          "[gra::MPhotoVM][gra::MPhotoZ][physics][covariance][polarization]") {
  ModelParamRestoreGuard restore;
  const auto tune = gra::MModelTune::Load(modelfile);
  const auto check = [&](const std::string &channel, int pdg) {
    // An electron emitter and proton target isolate one production direction
    auto lts = MakeToyPhotoZFFbar(11, LoadedPDGTable().FindByPDG(pdg).mass);
    lts.beam1 = lts.PDG.FindByPDG(11);
    lts.pbeam1.SetPxPyPzM(0.0, 0.0, lts.pbeam1.Pz(), lts.beam1.mass);
    lts.pfinal[1].SetPxPyPzM(0.18, 0.0, lts.pfinal[1].Pz(), lts.beam1.mass);
    RefreshToyDerivedKinematicsPreserveDecay(lts);
    lts.forward_mass2 = {pow2(lts.beam1.mass), pow2(lts.beam2.mass)};
    std::unique_ptr<gra::MPhotoVM> vm;
    std::unique_ptr<gra::MPhotoZ> z;
    if (pdg == 23) { z = std::make_unique<gra::MPhotoZ>(lts, tune, gra::MPhotoZ::ProcessDefinitionFor()); }
    else { vm = std::make_unique<gra::MPhotoVM>(lts, tune, channel, gra::MPhotoVM::ProcessDefinitionFor(channel)); }
    const auto evaluate = [&](gra::LORENTZSCALAR event) { return z ? z->Amp2(event) : vm->Amp2(event); };

    SetPhotoDecay(lts, 0.4, gra::math::PI / 2.0);
    const double norm = evaluate(lts);
    REQUIRE(norm > 0.0);
    // [REFERENCE: Zha et al., Phys. Rev. D 103 (2021) 033007, Eq. (2)]
    // The electron mass correction is below the requested angular accuracy
    for (double cosine : {-0.7, 0.4, 0.8}) {
      for (double phi : {0.0, 0.37, gra::math::PI / 2.0}) {
        SetPhotoDecay(lts, cosine, phi);
        const double expected = 1.0 - (1.0 - pow2(cosine)) * pow2(std::cos(phi));
        CAPTURE(channel, cosine, phi);
        REQUIRE(evaluate(lts) / norm == Approx(expected).epsilon(2.0e-6));
        auto reversed = lts;
        std::swap(reversed.decaytree[0], reversed.decaytree[1]);
        REQUIRE(evaluate(reversed) == Approx(evaluate(lts)).epsilon(1.0e-12));
      }
    }
    // Integrating the polarized decay retains the transverse partial-width normalization
    const auto rule = gra::math::GaussLegendreRule(4, -1.0, 1.0);
    double integral = 0.0;
    for (const auto &i : indices(rule.first)) {
      for (int j = 0; j < 8; ++j) {
        SetPhotoDecay(lts, rule.first[i], 2.0 * gra::math::PI * (j + 0.5) / 8.0);
        integral += rule.second[i] * evaluate(lts) / norm / 16.0;
      }
    }
    REQUIRE(integral == Approx(2.0 / 3.0).epsilon(2.0e-6));

    // Nonzero central transverse momentum exposes use of a lab helicity axis
    lts.pfinal[2].SetPxPyPzM(-0.10, 0.04, lts.pfinal[2].Pz(), lts.beam2.mass);
    RefreshToyDerivedKinematicsPreserveDecay(lts);
    SetPhotoDecay(lts, 0.23, 0.41);
    const double original = evaluate(lts);
    REQUIRE(original > 0.0);
    for (double rapidity : {-1.0, -0.2, 0.2, 1.0}) {
      const auto boosted = BoostPhotoEvent(lts, rapidity);
      CAPTURE(channel, rapidity);
      REQUIRE(evaluate(boosted) == Approx(original).epsilon(2.0e-6));
    }
  };
  check("Z", 23);
  for (const auto &[channel, pdg] : std::array<std::pair<const char *, int>, 5>{
           std::pair{"jpsi", 443}, {"psi(2S)", 100443}, {"Upsilon(1S)", 553},
           {"Upsilon(2S)", 100553}, {"Upsilon(3S)", 200553}}) { check(channel, pdg); }
}

TEST_CASE("PhotoVM rejects invalid integration counts before quadrature", "[gra::MPhotoVM][validation]") {
  const auto tune = gra::MModelTune::Load(modelfile);
  const std::array<nlohmann::json, 6> invalid = {-1, 0, 3, 4.5, "48", 4294967296ULL};
  for (const auto &count : invalid) {
    auto input = tune->Numerics();
    input["NUMERICS_PHOTOVM"]["N_k"] = count;
    gra::MPhotoVMNumerics parameters;
    CAPTURE(count);
    REQUIRE_THROWS_AS(parameters.Configure(tune->General(), "GENERAL", input, "NUMERICS"), std::invalid_argument);
  }
}

TEST_CASE("Photo Prod3 and decay keep both beam orderings in one collider section",
          "[gra::rspin][photoproduction][helicity][phase][interference]") {
  gra::LORENTZSCALAR lts = MakeToyCoherentPhotonLTS();
  gra::RequireModelCache(lts.model_cache, gra::MModelTune::Load(modelfile), "test resonance photoproduction");
  lts.process.FORWARD_NOFLIP = false;
  lts.process.FORWARD_VERTEX = gra::ForwardVertexMode::HelicityResidue;
  lts.process.PHOTON_VERTEX  = "QED";
  ConfigureToyPionPair(lts);

  gra::PARAM_RES upper = MakeToyCovariantScalarPomeronPhotoXP();
  gra::PARAM_RES lower = upper;
  std::swap(lower.production[0].tree[0], lower.production[0].tree[1]);
  PrepareToyXPOperators(lower, {{{0, 2, 1.0}}}, true, true);

  const auto full_amplitude = [](gra::LORENTZSCALAR event, const gra::PARAM_RES &resonance) {
    const auto production = gra::rspin::Resonance(event, resonance, (resonance.production_model == gra::ReggeProductionModel::MP ? gra::mpom::Fusion : gra::xpom::Fusion), 1.0);
    const auto decay      = gra::spin::ResonanceDecayMatrix(event, resonance, "CM");
    REQUIRE(production.size() == 1);
    return production[0] * decay;
  };

  constexpr double         rotation        = -0.73;
  const gra::LORENTZSCALAR rotated_lts     = RotateToyEventAroundZ(lts, rotation);
  const auto               upper_reference = full_amplitude(lts, upper);
  const auto               lower_reference = full_amplitude(lts, lower);
  const auto               upper_rotated   = full_amplitude(rotated_lts, upper);
  const auto               lower_rotated   = full_amplitude(rotated_lts, lower);

  REQUIRE(upper_reference.size_row() == 16);
  REQUIRE(lower_reference.size_row() == upper_reference.size_row());
  REQUIRE(upper_reference.size_col() == lower_reference.size_col());
  REQUIRE(upper_rotated.size_row() == upper_reference.size_row());
  REQUIRE(lower_rotated.size_row() == lower_reference.size_row());

  constexpr auto helicity_x2      = gra::spin::BinaryHelicityLabelsX2();
  double         upper_norm2      = 0.0;
  double         lower_norm2      = 0.0;
  double         incoherent_norm2 = 0.0;
  double         interference_l1  = 0.0;
  for (std::size_t in1 = 0; in1 < 2; ++in1) {
    for (std::size_t in2 = 0; in2 < 2; ++in2) {
      for (std::size_t out1 = 0; out1 < 2; ++out1) {
        for (std::size_t out2 = 0; out2 < 2; ++out2) {
          const std::size_t row      = gra::spin::CanonicalProtonPairSpinLayout::HardRow(in1, in2, out1, out2);
          const int         harmonic = gra::spin::ColliderSpinHalfHelicityHarmonic(helicity_x2[in1], helicity_x2[in2],
                                                                                   helicity_x2[out1], helicity_x2[out2]);
          const std::complex<double> phase = std::exp(gra::math::zi * static_cast<double>(harmonic) * rotation);
          for (std::size_t col = 0; col < upper_reference.size_col(); ++col) {
            const std::complex<double> up       = upper_reference[row][col];
            const std::complex<double> down     = lower_reference[row][col];
            const std::complex<double> coherent = up + down;
            const double               scale    = std::max({1.0, std::abs(up), std::abs(down), std::abs(coherent)});
            CAPTURE(in1, in2, out1, out2, row, col, harmonic, up, down, upper_rotated[row][col],
                    lower_rotated[row][col]);
            CHECK(std::abs(upper_rotated[row][col] - phase * up) < 5.0e-10 * scale);
            CHECK(std::abs(lower_rotated[row][col] - phase * down) < 5.0e-10 * scale);
            CHECK(std::abs(upper_rotated[row][col] + lower_rotated[row][col] - phase * coherent) < 7.0e-10 * scale);
            upper_norm2 += std::norm(up);
            lower_norm2 += std::norm(down);
            incoherent_norm2 += std::norm(up) + std::norm(down);
            interference_l1 += std::abs(std::norm(coherent) - std::norm(up) - std::norm(down));
          }
        }
      }
    }
  }

  CAPTURE(upper_norm2, lower_norm2, interference_l1, incoherent_norm2);
  REQUIRE(upper_norm2 > 0.0);
  REQUIRE(lower_norm2 > 0.0);
  REQUIRE(interference_l1 > 1.0e-10 * incoherent_norm2);
}

TEST_CASE("Nuclear ygg amplitudes resolve emitter and target sectors",
          "[gra::MPhotoZ][gra::MPhotoVM][nuclear][good-walker]") {
  ModelParamRestoreGuard restore;
  gra::MODELPARAM    = "TUNE0";
  const auto nucleus = MakePhotoOxygen();
  const auto inclusive =
      MakePhotoSectorUPC(nucleus, gra::nuclear::CoherenceType::Inclusive, gra::nuclear::CoherenceType::Inclusive,
                         gra::nuclear::SurvivalType::Optical);
  const auto coherent_config =
      MakePhotoSectorUPC(nucleus, gra::nuclear::CoherenceType::Coherent, gra::nuclear::CoherenceType::Coherent,
                         gra::nuclear::SurvivalType::MCGGCF);

  const auto require_channels = [](const gra::LORENTZSCALAR &lts, const std::size_t channel_count) {
    REQUIRE(lts.hamp.size() % channel_count == 0);
    const auto &metadata = lts.hamp.layout;
    REQUIRE(metadata.photo_sector_resolved);
    REQUIRE(metadata.photo_channel_count == channel_count);
    const auto require_channel =
        [&](const std::size_t index, const gra::nuclear::PhotoDirection direction,
            const gra::nuclear::CoherenceType emission, const gra::nuclear::CoherenceType target,
            const gra::nuclear::CoherenceType upper, const gra::nuclear::CoherenceType lower) {
          const auto &channel = metadata.photo_channel[index];
          CHECK(channel.direction == direction);
          CHECK(channel.emission == emission);
          CHECK(channel.target == target);
          CHECK(channel.Pair() == gra::nuclear::FinalPairIndex(upper, lower));
        };
    using gra::nuclear::CoherenceType;
    using gra::nuclear::PhotoDirection;
    require_channel(0, PhotoDirection::Upper, CoherenceType::Coherent, CoherenceType::Coherent,
                    CoherenceType::Coherent, CoherenceType::Coherent);
    if (channel_count == 2) {
      require_channel(1, PhotoDirection::Lower, CoherenceType::Coherent, CoherenceType::Coherent,
                      CoherenceType::Coherent, CoherenceType::Coherent);
    } else {
      REQUIRE(channel_count == 8);
      require_channel(1, PhotoDirection::Upper, CoherenceType::Coherent, CoherenceType::Incoherent,
                      CoherenceType::Coherent, CoherenceType::Incoherent);
      require_channel(2, PhotoDirection::Upper, CoherenceType::Incoherent, CoherenceType::Coherent,
                      CoherenceType::Incoherent, CoherenceType::Coherent);
      require_channel(3, PhotoDirection::Upper, CoherenceType::Incoherent, CoherenceType::Incoherent,
                      CoherenceType::Incoherent, CoherenceType::Incoherent);
      require_channel(4, PhotoDirection::Lower, CoherenceType::Coherent, CoherenceType::Coherent,
                      CoherenceType::Coherent, CoherenceType::Coherent);
      require_channel(5, PhotoDirection::Lower, CoherenceType::Coherent, CoherenceType::Incoherent,
                      CoherenceType::Incoherent, CoherenceType::Coherent);
      require_channel(6, PhotoDirection::Lower, CoherenceType::Incoherent, CoherenceType::Coherent,
                      CoherenceType::Coherent, CoherenceType::Incoherent);
      require_channel(7, PhotoDirection::Lower, CoherenceType::Incoherent, CoherenceType::Incoherent,
                      CoherenceType::Incoherent, CoherenceType::Incoherent);
    }
    std::vector<double> norm(channel_count, 0.0);
    for (const auto &i : gra::aux::indices(lts.hamp)) { norm[i % channel_count] += std::norm(lts.hamp[i]); }
    for (const double value : norm) { REQUIRE(value > 0.0); }
    if (channel_count != 8) { return; }

    const auto upper        = gra::ResolveForwardLegState(lts, gra::ForwardBeamLeg::Upper);
    const auto lower        = gra::ResolveForwardLegState(lts, gra::ForwardBeamLeg::Lower);
    const auto source_ratio = [](const gra::ForwardLegState &state) {
      const double coherent =
          gra::flux::NuclearPhotonFluxTransverse(state, gra::nuclear::CoherenceType::Coherent).Trace();
      const double incoherent =
          gra::flux::NuclearPhotonFluxTransverse(state, gra::nuclear::CoherenceType::Incoherent).Trace();
      REQUIRE(coherent > 0.0);
      REQUIRE(incoherent > 0.0);
      return std::sqrt(incoherent / coherent);
    };
    const double upper_ratio = source_ratio(upper);
    const double lower_ratio = source_ratio(lower);
    for (std::size_t base = 0; base < lts.hamp.size(); base += channel_count) {
      for (const auto &[coherent, incoherent, ratio] : std::array<std::tuple<std::size_t, std::size_t, double>, 4>{
               std::tuple{base, base + 2, upper_ratio}, std::tuple{base + 1, base + 3, upper_ratio},
               std::tuple{base + 4, base + 6, lower_ratio}, std::tuple{base + 5, base + 7, lower_ratio}}) {
        const std::complex<double> expected = ratio * lts.hamp[coherent];
        const double               scale    = std::max({1.0, std::abs(expected), std::abs(lts.hamp[incoherent])});
        CHECK(std::abs(lts.hamp[incoherent] - expected) < 2.0e-12 * scale);
      }
    }
  };

  const auto vm_size = [&](const std::shared_ptr<const gra::nuclear::MUPC> &upc, const std::size_t channel_count) {
    gra::LORENTZSCALAR lts = MakeToyIonPhotoPair(13, 3.096900, nucleus);
    lts.upc_model          = upc;
    gra::MPhotoVM amplitude(lts, gra::MModelTune::Load(modelfile), "jpsi", gra::MPhotoVM::ProcessDefinitionFor("jpsi"));
    const double  amp2 = amplitude.Amp2(lts);
    REQUIRE(std::isfinite(amp2));
    REQUIRE(amp2 > 0.0);
    REQUIRE(gra::AllFinite(lts.hamp));
    require_channels(lts, channel_count);
    if (upc->Bank(1) != nullptr) {
      for (const auto &current : lts.screening.photo.current) {
        REQUIRE(current.has_value());
        REQUIRE(current->sample.size() == upc->Bank(1)->Size());
      }
    }
    return lts.hamp.size();
  };
  const auto z_size = [&](const std::shared_ptr<const gra::nuclear::MUPC> &upc, const std::size_t channel_count) {
    gra::LORENTZSCALAR lts = MakeToyIonPhotoPair(13, 91.1876, nucleus);
    lts.upc_model          = upc;
    gra::MPhotoZ amplitude(lts, gra::MModelTune::Load(modelfile), gra::MPhotoZ::ProcessDefinitionFor());
    const double amp2 = amplitude.Amp2(lts);
    REQUIRE(std::isfinite(amp2));
    REQUIRE(amp2 > 0.0);
    REQUIRE(gra::AllFinite(lts.hamp));
    require_channels(lts, channel_count);
    if (upc->Bank(1) != nullptr) {
      for (const auto &current : lts.screening.photo.current) {
        REQUIRE(current.has_value());
        REQUIRE(current->sample.size() == upc->Bank(1)->Size());
      }
    }
    return lts.hamp.size();
  };

  // Four final spin rows carry four sectors per direction for inclusive AA
  REQUIRE(vm_size(inclusive, 8) == 32);
  REQUIRE(z_size(inclusive, 8) == 32);
  // Configuration survival retains all sources before selecting the coherent final state
  REQUIRE(vm_size(coherent_config, 8) == 32);
  REQUIRE(z_size(coherent_config, 8) == 32);
}

TEST_CASE("Photonuclear embedding preserves the nuclear target invariant flux",
          "[gra::flux][photoproduction][nuclear][normalization][physics]") {
  const auto         nucleus = MakePhotoOxygen();
  const auto         aa      = MakePhotoSectorUPC(nucleus, gra::nuclear::CoherenceType::Coherent,
                                                  gra::nuclear::CoherenceType::Coherent, gra::nuclear::SurvivalType::Optical);
  gra::LORENTZSCALAR lts     = MakeToyIonPhotoPair(13, 3.096900, nucleus);
  gra::RequireModelCache(lts.model_cache, gra::MModelTune::Load(modelfile), "test photonuclear embedding");
  lts.upc_model                              = aa;
  const auto                          target = gra::ResolveForwardLegState(lts, gra::ForwardBeamLeg::Lower);
  const gra::flux::PhotoTargetProfile profile{0.0, 4.0, 0.0};
  const auto                          embedded = gra::flux::ResolvePhotoTarget(target, profile, 1.0, 1.0);
  REQUIRE(embedded.factor.size() == 1);

  const gra::M3Vec transfer = gra::nuclear::RestTransfer(target.incoming, target.transfer, target.emitter.mass);
  const gra::nuclear::PhotoProfile nuclear_profile{0.0, 4.0, 0.0, 0.0};
  const auto                 transition = aa->Photo(2)->Factors(nuclear_profile, transfer[0], transfer[1], transfer[2],
                                                                gra::nuclear::TargetPhotonDirection(2));
  const std::complex<double> expected   = static_cast<double>(nucleus->A()) * transition.coherent;
  CHECK(embedded.factor[0].real() == Approx(expected.real()).epsilon(2.0e-12));
  CHECK(embedded.factor[0].imag() == Approx(expected.imag()).epsilon(2.0e-12));

  gra::LORENTZSCALAR pp = MakeToyPhotoZFFbar(13, 3.096900);
  gra::RequireModelCache(pp.model_cache, gra::MModelTune::Load(modelfile), "test elementary embedding");
  const auto proton     = gra::ResolveForwardLegState(pp, gra::ForwardBeamLeg::Lower);
  const auto elementary = gra::flux::ResolvePhotoTarget(proton, profile, 2.0, 3.0);
  REQUIRE(elementary.factor.size() == 1);
  CHECK(elementary.factor[0].real() == Approx(3.0));
  CHECK(elementary.factor[0].imag() == Approx(0.0).margin(1.0e-15));
}

TEST_CASE("Configuration ygg screening combines matching photon directions",
          "[gra::MProcess][gra::MPhotoVM][nuclear][good-walker]") {
  ModelParamRestoreGuard restore;
  gra::MODELPARAM                         = "TUNE0";
  const auto                      nucleus = MakePhotoOxygen();
  const auto                      upc     = MakePhotoSectorUPC(nucleus, gra::nuclear::CoherenceType::Inclusive,
                                                               gra::nuclear::CoherenceType::Inclusive, gra::nuclear::SurvivalType::MCGGCF);
  ToyNuclearPhotoScreeningProcess process(gra::MModelTune::Load(modelfile), nucleus, upc);

  const double screened = process.ScreenedAmp2();
  REQUIRE(std::isfinite(screened));
  REQUIRE(screened > 0.0);
  REQUIRE(process.state.lts.hamp.layout.photo_sector_resolved);
  REQUIRE(process.state.lts.hamp.layout.photo_channel_count == 8);
  REQUIRE(process.state.lts.hamp.size() == 32);

  double projected_norm = 0.0;
  for (std::size_t hard = 0; hard < process.state.lts.hamp.size() / 8; ++hard) {
    const std::size_t base = hard * 8;
    for (std::size_t channel = 0; channel < 4; ++channel) {
      projected_norm += std::norm(process.state.lts.hamp[base + channel]);
    }
    // Lower-direction entries have been coherently folded into the matching
    // final-leg sector represented by the first four entries
    for (std::size_t channel = 4; channel < 8; ++channel) {
      CHECK(std::fpclassify(std::norm(process.state.lts.hamp[base + channel])) == FP_ZERO);
    }
  }
  CHECK(projected_norm > 0.0);
}

TEST_CASE("Bare and screened ygg amplitudes project sampled final nuclear sectors",
          "[gra::MProcess][gra::MPhotoZ][gra::MPhotoVM][nuclear]"
          "[interference]") {
  ModelParamRestoreGuard restore;
  gra::MODELPARAM       = "TUNE0";
  const auto model_tune = gra::MModelTune::Load(modelfile);
  const auto nucleus    = MakePhotoOxygen();
  const auto inclusive =
      MakePhotoSectorUPC(nucleus, gra::nuclear::CoherenceType::Inclusive, gra::nuclear::CoherenceType::Inclusive,
                         gra::nuclear::SurvivalType::Optical);

  for (const ToyPhotoType type : {ToyPhotoType::VM, ToyPhotoType::Z}) {
    const std::string label = type == ToyPhotoType::VM ? "vector meson" : "Z";
    DYNAMIC_SECTION(label) {
      ToyNuclearPhotoScreeningProcess bare(model_tune, nucleus, inclusive, type);
      gra::LORENTZSCALAR raw_lts = bare.state.lts;
      if (type == ToyPhotoType::VM) {
        gra::MPhotoVM amplitude(raw_lts, model_tune, "jpsi", gra::MPhotoVM::ProcessDefinitionFor("jpsi"));
        REQUIRE(amplitude.Amp2(raw_lts) > 0.0);
      } else {
        gra::MPhotoZ amplitude(raw_lts, model_tune, gra::MPhotoZ::ProcessDefinitionFor());
        REQUIRE(amplitude.Amp2(raw_lts) > 0.0);
      }
      const std::vector<std::complex<double>> raw = raw_lts.hamp;
      REQUIRE(raw.size() == 32);

      const auto channels = gra::nuclear::PhotoChannels(*inclusive);
      const std::array<gra::M3Vec, 2> transfer = {
          gra::nuclear::RestTransfer(raw_lts.pbeam1, raw_lts.q1, raw_lts.beam1.mass),
          gra::nuclear::RestTransfer(raw_lts.pbeam2, raw_lts.q2, raw_lts.beam2.mass)};
      using Ensemble = std::vector<std::complex<double>>;
      std::array<std::vector<Ensemble>, 2> direction;
      for (auto &bank : direction) { bank.assign(4, Ensemble(inclusive->SampleCount(), 0.0)); }
      // Reconstruct each complete photon direction before averaging nuclear configurations
      const auto accumulate = [&](const std::vector<std::complex<double>> &amplitude,
                                  const std::array<std::optional<gra::nuclear::PhotoCurrent>, 2> &current,
                                  std::optional<gra::nuclear::PhotoDirection> selected) {
        for (const auto &c : indices(channels)) {
          const auto &channel = channels[c];
          if (selected.has_value() && *selected != channel.direction) { continue; }
          const int emitter = channel.direction == gra::nuclear::PhotoDirection::Upper ? 1 : 2;
          const int target = 3 - emitter;
          const auto source = inclusive->EmissionRatios(emitter, channel.emission, transfer[emitter - 1]);
          const auto response = current[target - 1].has_value()
              ? inclusive->TargetCurrentRatios(target, channel.target, *current[target - 1])
              : inclusive->TargetRatios(target, channel.target, transfer[target - 1]);
          for (std::size_t hard = 0; hard < 4; ++hard) {
            auto &values = direction[emitter - 1][hard];
            for (const auto &sample : indices(values)) {
              values[sample] += amplitude[8 * hard + c] * source[sample] * response[sample];
            }
          }
        }
      };
      if (raw_lts.screening.photo.term.empty()) {
        accumulate(raw, raw_lts.screening.photo.current, std::nullopt);
      } else {
        for (const auto &term : raw_lts.screening.photo.term) {
          accumulate(term.amplitude, term.photo_current, term.direction);
        }
      }
      std::vector<std::complex<double>> expected(raw.size(), 0.0);
      std::array<double, 4> coherent_scale{};
      double expected_norm = 0.0, diagonal_norm = 0.0;
      for (std::size_t hard = 0; hard < 4; ++hard) {
        for (std::size_t pair = 0; pair < 4; ++pair) {
          const auto &upper = direction[0][hard];
          const auto &lower = direction[1][hard];
          Ensemble total(upper.size());
          for (const auto &sample : indices(total)) {
            total[sample] = upper[sample] + lower[sample];
            if (pair == 0) {
              coherent_scale[hard] += (std::abs(upper[sample]) + std::abs(lower[sample])) / total.size();
            }
          }
          expected_norm += gra::nuclear::NuclearGoodWalkerProject(total, inclusive->SampleShape())[pair / 2][pair % 2];
          diagonal_norm += gra::nuclear::NuclearGoodWalkerProject(upper, inclusive->SampleShape())[pair / 2][pair % 2] +
                           gra::nuclear::NuclearGoodWalkerProject(lower, inclusive->SampleShape())[pair / 2][pair % 2];
          if (pair == 0) {
            expected[8 * hard] = std::accumulate(total.begin(), total.end(), std::complex<double>{}) /
                                 static_cast<double>(total.size());
          }
        }
      }
      const double bare_amp2 = bare.ScreenedAmp2(false);
      REQUIRE(bare_amp2 == Approx(expected_norm).epsilon(1.0e-12));
      REQUIRE(bare.state.lts.hamp.size() == expected.size());
      for (const auto &h : indices(expected)) {
        // Bound roundoff relative to the interfering currents before their coherent cancellation
        const double scale = h % 8 == 0 ? std::max(1.0, coherent_scale[h / 8]) : 1.0;
        CAPTURE(h, expected[h], bare.state.lts.hamp[h], raw[8 * (h / 8)], raw[8 * (h / 8) + 4]);
        CHECK(std::abs(bare.state.lts.hamp[h] - expected[h]) < 2.0e-12 * scale);
      }
      REQUIRE(std::abs(expected_norm - diagonal_norm) > 1.0e-14 * diagonal_norm);

      for (const auto survival : {gra::nuclear::SurvivalType::Optical, gra::nuclear::SurvivalType::OpticalGGCF}) {
        const auto                      upc = MakePhotoSectorUPC(nucleus, gra::nuclear::CoherenceType::Inclusive,
                                                                 gra::nuclear::CoherenceType::Inclusive, survival);
        ToyNuclearPhotoScreeningProcess screened(model_tune, nucleus, upc, type);
        const double                    amp2 = screened.ScreenedAmp2();
        REQUIRE(std::isfinite(amp2));
        REQUIRE(amp2 > 0.0);
        for (std::size_t hard = 0; hard < screened.state.lts.hamp.size() / 8; ++hard) {
          for (std::size_t channel = 4; channel < 8; ++channel) {
            CHECK(std::fpclassify(std::norm(screened.state.lts.hamp[8 * hard + channel])) == FP_ZERO);
          }
        }
      }
    }
  }
}

TEST_CASE("Forward leg states resolve upper and lower generated systems", "[gra::ForwardLegState][kinematics]") {
  gra::LORENTZSCALAR lts = MakeToyPhotoZFFbar(13);
  gra::RequireModelCache(lts.model_cache, gra::MModelTune::Load(modelfile), "test forward leg");
  SetToyPhotoForwardExcitation(lts, 1, 2.0);
  SetToyPhotoForwardExcitation(lts, 2, 3.0);

  const auto upper = gra::ResolveForwardLegState(lts, gra::ForwardBeamLeg::Upper);
  const auto lower = gra::ResolveForwardLegState(lts, gra::ForwardBeamLeg::Lower);
  REQUIRE(upper.Index() == 1);
  REQUIRE(upper.FinalIndex() == 1);
  REQUIRE(upper.IsExcited());
  CHECK(upper.xi == Approx(lts.xi1));
  CHECK(upper.t == Approx(lts.t1));
  CHECK(upper.qt == Approx(lts.qt1));
  CHECK(upper.mass2 == Approx(lts.pfinal[1].M2()));
  CHECK(upper.incoming == lts.pbeam1);
  CHECK(upper.outgoing == lts.pfinal[1]);
  CHECK(upper.transfer == lts.q1);

  REQUIRE(lower.Index() == 2);
  REQUIRE(lower.FinalIndex() == 2);
  REQUIRE(lower.IsExcited());
  CHECK(lower.xi == Approx(lts.xi2));
  CHECK(lower.t == Approx(lts.t2));
  CHECK(lower.qt == Approx(lts.qt2));
  CHECK(lower.mass2 == Approx(lts.pfinal[2].M2()));
  CHECK(lower.incoming == lts.pbeam2);
  CHECK(lower.outgoing == lts.pfinal[2]);
  CHECK(lower.transfer == lts.q2);

  gra::LORENTZSCALAR missing = lts;
  missing.pfinal.resize(2);
  CHECK_THROWS(gra::ResolveForwardLegState(missing, gra::ForwardBeamLeg::Lower));

  gra::LORENTZSCALAR invalid = lts;
  invalid.qt1                = -1.0;
  CHECK_THROWS(gra::ResolveForwardLegState(invalid, gra::ForwardBeamLeg::Upper));
}

TEST_CASE("PhotoZ parameters construct safely from one tune across threads", "[gra::MPhotoZ][threading]") {
  ModelParamRestoreGuard restore;
  gra::MODELPARAM = "TUNE0";

  const auto       model_tune = gra::MModelTune::Load(modelfile);
  gra::MModelCache cache(model_tune);
  const auto       first        = gra::GetPhotoQCDNumerics(cache, "NUMERICS_PHOTOZ");
  const auto       second       = gra::GetPhotoQCDNumerics(cache, "NUMERICS_PHOTOZ");
  const auto       first_param  = gra::GetPhotoQCDParam(cache, "PARAM_PHOTOZ");
  const auto       second_param = gra::GetPhotoQCDParam(cache, "PARAM_PHOTOZ");
  REQUIRE(first == second);
  REQUIRE(first_param == second_param);
  REQUIRE(first->initialized);
  REQUIRE(second->initialized);
  REQUIRE(first_param->initialized);
  REQUIRE(second_param->initialized);

  constexpr std::size_t                                      nthreads = 8;
  std::vector<std::shared_ptr<const gra::MPhotoQCDNumerics>> handles(nthreads);
  std::vector<std::shared_ptr<const gra::MPhotoQCDParam>>    param_handles(nthreads);
  std::vector<std::thread>                                   workers;
  workers.reserve(nthreads);

  for (std::size_t i = 0; i < nthreads; ++i) {
    workers.emplace_back([i, &handles, &param_handles, &cache] {
      handles[i]       = gra::GetPhotoQCDNumerics(cache, "NUMERICS_PHOTOZ");
      param_handles[i] = gra::GetPhotoQCDParam(cache, "PARAM_PHOTOZ");
    });
  }
  for (auto &worker : workers) { worker.join(); }

  for (const auto &handle : handles) { REQUIRE(handle == first); }
  for (const auto &handle : param_handles) { REQUIRE(handle == first_param); }
}

// Require photoproduction parameters to come from complete tune blocks
TEST_CASE("Photoproduction readers reject incomplete tune blocks", "[gra::MPhotoZ][gra::MPhotoVM][snapshot]") {
  const auto photoz_tune = WriteModifiedPhotoZTune(
      "photoz_missing_numerics", [](auto &) {}, [](auto &j) { j.erase("NUMERICS_PHOTOZ"); });
  const auto photoz_model = gra::MModelTune::Load(photoz_tune.second);
  REQUIRE_THROWS(gra::ReadPhotoQCDNumerics(*photoz_model, "NUMERICS_PHOTOZ"));

  const auto photovm_tune =
      WriteModifiedPhotoVMTune("photovm_missing_parameters", [](auto &j) { j.erase("PARAM_PHOTOVM"); });
  const auto photovm_model = gra::MModelTune::Load(photovm_tune.second);
  REQUIRE_THROWS(gra::ReadPhotoVMNumerics(*photovm_model));
}

TEST_CASE("PhotoZ readers use changed card content", "[gra::MPhotoZ][threading]") {
  ModelParamRestoreGuard restore;

  const auto first_tune = WriteModifiedPhotoZTune(
      "photoqcd_cache_reload", [](auto &j) { j["PARAM_PHOTOZ"]["real_part_delta"] = 0.11; },
      [](auto &j) { j["NUMERICS_PHOTOZ"]["N_k"] = 3; });
  gra::MODELPARAM           = first_tune.first;
  const auto first_model    = gra::MModelTune::Load(first_tune.second);
  const auto first_param    = gra::ReadPhotoQCDParam(*first_model, "PARAM_PHOTOZ");
  const auto first_numerics = gra::ReadPhotoQCDNumerics(*first_model, "NUMERICS_PHOTOZ");
  REQUIRE(first_param->real_part_delta == Approx(0.11));
  REQUIRE(first_numerics->N_k == 3);

  const auto second_tune = WriteModifiedPhotoZTune(
      "photoqcd_cache_reload", [](auto &j) { j["PARAM_PHOTOZ"]["real_part_delta"] = 0.17; },
      [](auto &j) { j["NUMERICS_PHOTOZ"]["N_k"] = 5; });
  gra::MODELPARAM         = second_tune.first;
  const auto second_model = gra::MModelTune::Load(second_tune.second);
  REQUIRE(second_model->GeneralFile() == first_model->GeneralFile());
  REQUIRE(second_model != first_model);
  const auto second_param    = gra::ReadPhotoQCDParam(*second_model, "PARAM_PHOTOZ");
  const auto second_numerics = gra::ReadPhotoQCDNumerics(*second_model, "NUMERICS_PHOTOZ");
  REQUIRE(second_param->real_part_delta == Approx(0.17));
  REQUIRE(second_numerics->N_k == 5);
}

TEST_CASE("MPhotoZ initializes through the shared Sudakov UGD store", "[gra::MPhotoZ][sudakov]") {
  ModelParamRestoreGuard restore;
  gra::MODELPARAM = "nonexistent_global_tune";

  gra::LORENTZSCALAR first_lts = MakeToyPhotoZFFbar(13);
  REQUIRE(first_lts.GlobalSudakovPtr == nullptr);
  REQUIRE(first_lts.GlobalPdfPtr == nullptr);
  gra::MPhotoZ first(first_lts, gra::MModelTune::Load(modelfile), gra::MPhotoZ::ProcessDefinitionFor());
  const gra::amplitude::ProcessDefinition &process_definition = first;
  REQUIRE(process_definition.DecayStructureFor(first_lts) == gra::MPhotoZ::DirectDecayStructure());
  REQUIRE(first_lts.GlobalSudakovPtr != nullptr);
  REQUIRE(first_lts.GlobalPdfPtr == nullptr);

  gra::LORENTZSCALAR second_lts = MakeToyPhotoZFFbar(13);
  second_lts.model_cache        = first_lts.model_cache;
  gra::MPhotoZ second(second_lts, gra::MModelTune::Load(modelfile), gra::MPhotoZ::ProcessDefinitionFor());
  REQUIRE(second_lts.GlobalSudakovPtr == first_lts.GlobalSudakovPtr);
}

TEST_CASE("MPhotoZ rejects a Sudakov from a different SOFT snapshot", "[gra::MPhotoZ][sudakov][snapshot]") {
  ModelParamRestoreGuard restore;
  gra::MODELPARAM = "TUNE0";

  const auto first_tune =
      WriteModifiedPhotoVMTune("photoz_first_soft", [](auto &j) { j["PARAM_REGGE"]["omega"]["MP"] = 0.81; });
  const auto second_tune =
      WriteModifiedPhotoVMTune("photoz_second_soft", [](auto &j) { j["PARAM_REGGE"]["omega"]["MP"] = 0.97; });
  const auto first_model  = gra::MModelTune::Load(first_tune.second);
  const auto second_model = gra::MModelTune::Load(second_tune.second);
  REQUIRE(first_model != second_model);

  gra::LORENTZSCALAR lts = MakeToyPhotoZFFbar(13);
  lts.GlobalSudakovPtr   = lts.model_cache->sudakov.GetSudakov(lts.sqrt_s, lts.LHAPDFSET, first_model->Soft());
  REQUIRE(lts.GlobalSudakovPtr->SoftModelHandle() == first_model->Soft());
  REQUIRE_THROWS(gra::MPhotoZ(lts, second_model, gra::MPhotoZ::ProcessDefinitionFor()));
}

// Reject the whole impact amplitude when the infrared gluon fails at an integration point
TEST_CASE("MPhotoZ preserves gluon amplitude failures", "[gra::MPhotoZ][sudakov][regression]") {
  ModelParamRestoreGuard restore;
  const auto tune = WriteModifiedPhotoZTune(
      "failed_gluon", [](auto &j) {
        j["PARAM_SKEWED_UGD"]["mode"] = "IR_RKHS";
        j["PARAM_SKEWED_UGD"]["rkhs_radius"] = 100.0;
      }, [](auto &j) {
        j["NUMERICS_SUDAKOV"]["SUDA"]["N"] = {12, 12};
        j["NUMERICS_SUDAKOV"]["SHUV"]["N"] = {12, 12};
        // Sample the infrared region where the configured RKHS ball violates positivity
        j["NUMERICS_PHOTOZ"]["kappa2_MIN"] = 0.1;
      });
  gra::MODELPARAM = tune.first;
  const auto model = gra::MModelTune::Load(tune.second);
  gra::LORENTZSCALAR lts = MakeToyPhotoZFFbar(13);
  gra::MPhotoZ photoz(lts, model, gra::MPhotoZ::ProcessDefinitionFor());
  REQUIRE_THROWS_AS(photoz.GammaPTotalCrossSectionNb(lts, 200.0), gra::AmplitudeFailure);
  REQUIRE_THROWS_AS(photoz.Amp2(lts), gra::AmplitudeFailure);
}

// Keep optional nuclear PDF failures separate from the proton impact amplitude
TEST_CASE("PhotoZ low-scale nuclear profiles reject sampling points safely", "[gra::MPhotoZ][regression]") {
  const auto files = WriteModifiedPhotoZTune(
      "photoz_low_scale",
      [](auto &j) {
        j["PARAM_PHOTOZ"]["mu_MIN"]    = 0.2;
        j["PARAM_PHOTOZ"]["mu_over_m"] = 0.1;
      },
      [](auto &j) {
        j["NUMERICS_SUDAKOV"]["SHUV"]["N"] = {32, 32};
        j["NUMERICS_SUDAKOV"]["SUDA"]["N"] = {32, 32};
      });
  const auto   tune   = gra::MModelTune::Load(files.second);
  auto         proton = MovingPhotoPairEvent(13, 0.8);
  gra::MPhotoZ amplitude(proton, tune, gra::MPhotoZ::ProcessDefinitionFor());
  REQUIRE(std::isfinite(amplitude.Amp2(proton)));
  const auto nucleus = MakePhotoOxygen();
  auto       ion     = MakeToyIonPhotoPair(13, 0.8, nucleus);
  ion.upc_model      = MakePhotoSectorUPC(nucleus, gra::nuclear::CoherenceType::Coherent,
                                          gra::nuclear::CoherenceType::Coherent, gra::nuclear::SurvivalType::Optical);
  gra::MPhotoZ nuclear(ion, tune, gra::MPhotoZ::ProcessDefinitionFor());
  REQUIRE_THROWS_AS(nuclear.Amp2(ion), gra::AmplitudeFailure);
}

// Check the analytic z integration against the original CSS integral and its Plemelj discontinuity
TEST_CASE("CSS flavour integral retains the timelike pole phase", "[gra::MPhotoZ][physics][pole]") {
  const auto   tune = gra::MModelTune::Load(modelfile);
  auto         lts  = MakeToyPhotoZFFbar(13);
  gra::MPhotoZ photoz(lts, tune, gra::MPhotoZ::ProcessDefinitionFor());
  const auto   param = gra::ReadPhotoQCDParam(*tune, "PARAM_PHOTOZ");
  auto         num   = *gra::ReadPhotoQCDNumerics(*tune, "NUMERICS_PHOTOZ");
  num.N_k            = 96;
  const double m2    = 1.0;
  const double x     = 0.001;
  const double mu    = 5.0;
  // Integrate the published impact kernel with the actual Shuvaev gluon
  const auto impact = [&](double z, double k2) {
    double out = 0.0;
    for (const auto &j : indices(num.kappa2_nodes)) {
      const double kappa2 = num.kappa2_nodes[j];
      const double root   = std::sqrt(pow2(k2 - m2 - kappa2) + 4.0 * m2 * k2);
      const double w0     = 1.0 / (k2 + m2) - 1.0 / root;
      const double w1     = 1.0 - (k2 + m2) / (2.0 * k2) * (1.0 + (k2 - m2 - kappa2) / root);
      out += num.kappa2_weights[j] / pow2(kappa2) *
             lts.GlobalSudakovPtr->AlphaSFlux_xQ2Mu(x, kappa2, mu, std::max(k2 + m2, kappa2)) *
             (m2 * w0 + (pow2(z) + pow2(1.0 - z)) * k2 / (k2 + m2) * w1);
    }
    return gra::math::pow3(gra::math::PI) * out;
  };
  for (const double q2 : {1.0e-12, 0.5, 25.0, 4.0 * (m2 + num.kappa2_nodes.front())}) {
    const auto   value        = gra::MPhotoQCD::CSSFlavorIntegral(lts, *param, num, x, std::sqrt(m2), q2, mu);
    const bool   timelike_cut = q2 > 4.0 * (m2 + num.k2_MIN);
    const double zmin         = timelike_cut ? 0.5 * (1.0 - std::sqrt(1.0 - 4.0 * (m2 + num.k2_MIN) / q2)) : 0.0;
    const auto [z, weight]    = gra::math::GaussLegendreRule(256, zmin, 1.0 - zmin);
    double reference          = 0.0;
    for (const auto &i : indices(z)) {
      if (timelike_cut) {
        reference += weight[i] * gra::math::PI * impact(z[i], z[i] * (1.0 - z[i]) * q2 - m2);
      } else {
        reference += weight[i] * gra::math::LogGaussIntegral(256, num.k2_MIN, num.K2Max(q2), [&](double k2) {
                       return impact(z[i], k2) / (k2 + m2 - z[i] * (1.0 - z[i]) * q2);
                     });
      }
    }
    CAPTURE(q2, value, reference);
    CHECK((timelike_cut ? value.imag() : value.real()) == Approx(reference).epsilon(2.0e-4));
    if (!timelike_cut) { CHECK(std::abs(value.imag()) < 1.0e-14); }
    auto coarse      = num;
    coarse.N_k       = 48;
    const auto other = gra::MPhotoQCD::CSSFlavorIntegral(lts, *param, coarse, x, std::sqrt(m2), q2, mu);
    CHECK(std::abs(other - value) < 2.0e-4 * std::abs(value));
  }
}

// Require stable complex amplitudes as both transverse quadratures are refined
TEST_CASE("CSS light and heavy quark loops converge", "[gra::MPhotoZ][physics][convergence]") {
  const auto   tune = gra::MModelTune::Load(modelfile);
  auto         lts  = MakeToyPhotoZFFbar(13);
  gra::MPhotoZ photoz(lts, tune, gra::MPhotoZ::ProcessDefinitionFor());
  const auto   param = gra::ReadPhotoQCDParam(*tune, "PARAM_PHOTOZ");
  auto         num   = *gra::ReadPhotoQCDNumerics(*tune, "NUMERICS_PHOTOZ");
  const double q2    = pow2(lts.PDG.FindByPDG(23).mass);
  for (const double mass : {0.001, 1.5, 4.5}) {
    std::complex<double> previous;
    for (const unsigned int n : {24U, 48U, 96U}) {
      num.N_k     = n;
      num.N_kappa = n;
      num.PrepareIntegrationRules();
      const auto value = gra::MPhotoQCD::CSSFlavorIntegral(lts, *param, num, 0.001, mass, q2, std::sqrt(q2) / 2.0);
      CAPTURE(mass, n, value, previous);
      if (n > 24) { CHECK(std::abs(value - previous) < 0.003 * std::abs(value)); }
      previous = value;
    }
  }
}

TEST_CASE("MPhotoZ gamma-p benchmark has the CSS physical scale and energy growth", "[gra::MPhotoZ][benchmark]") {
  ModelParamRestoreGuard restore;
  gra::MODELPARAM = "TUNE0";

  gra::LORENTZSCALAR lts = MakeToyPhotoZFFbar(13);
  gra::MPhotoZ       photoz(lts, gra::MModelTune::Load(modelfile), gra::MPhotoZ::ProcessDefinitionFor());
  const double       z_threshold = lts.PDG.FindByPDG(23).mass + gra::PDG::mp;
  REQUIRE(photoz.GammaPTotalCrossSectionNb(lts, z_threshold - 1.0e-6) == Approx(0.0));
  // The CSS impact factor already includes the colour trace
  // [REFERENCE: Cisek, Schafer, Szczurek, arXiv:0906.1739, Eqs. (2.1)-(2.3)]
  const auto           tune   = gra::MModelTune::Load(modelfile);
  const auto           param  = gra::ReadPhotoQCDParam(*tune, "PARAM_PHOTOZ");
  const auto           num    = gra::ReadPhotoQCDNumerics(*tune, "NUMERICS_PHOTOZ");
  const double         mass2  = pow2(lts.PDG.FindByPDG(23).mass);
  const double         sw2    = 1.0 - pow2(lts.PDG.FindByPDG(24).mass) / mass2;
  const double         alpha  = gra::qed::alpha_QED(mass2, tune->Structure().QED_alpha, 1.0 / 132.507);
  std::complex<double> impact = 0.0;
  for (const int pdg : {1, 2, 3, 4, 5}) {
    const auto  &quark  = lts.PDG.FindByPDG(pdg);
    const double charge = quark.chargeX3 / 3.0;
    const double t3     = charge > 0.0 ? 0.5 : -0.5;
    impact += 2.0 * alpha / gra::math::PI * charge * (t3 - 2.0 * charge * sw2) / (2.0 * std::sqrt(sw2 * (1.0 - sw2))) *
              gra::MPhotoQCD::CSSFlavorIntegral(lts, *param, *num, gra::MPhotoQCD::XGluon(*param, mass2, pow2(200.0)),
                                                quark.mass, mass2, gra::MPhotoQCD::HardScale(*param, mass2));
  }
  const double expected = std::norm(gra::MPhotoQCD::RealPartFactor(*param) * impact) /
                          (16.0 * gra::math::PI * gra::MPhotoQCD::TSlope(*param, pow2(200.0))) * gra::PDG::GeV2barn *
                          1.0e9;
  CHECK(photoz.GammaPTotalCrossSectionNb(lts, 200.0) == Approx(expected).epsilon(1.0e-12));
  const double sigma_hera_nb  = photoz.GammaPTotalCrossSectionNb(lts, 200.0);
  const double sigma_10tev_nb = photoz.GammaPTotalCrossSectionNb(lts, 10000.0);

  INFO("sigma(gamma p -> Z p, W=200 GeV) = " << sigma_hera_nb << " nb");
  INFO("sigma(gamma p -> Z p, W=10 TeV) = " << sigma_10tev_nb << " nb");
  // The CSS anchors are order 1e-5 nb and 1e-3 nb for their different IN UGD
  // MMHT2014lo68cl gives stronger small-x growth at the upper benchmark energy
  REQUIRE(sigma_hera_nb > 5.0e-7);
  REQUIRE(sigma_hera_nb < 2.0e-4);
  REQUIRE(sigma_10tev_nb > 5.0e-5);
  REQUIRE(sigma_10tev_nb < 1.0e-1);
  REQUIRE(sigma_10tev_nb > 10.0 * sigma_hera_nb);
}

TEST_CASE(
    "MPhotoZ transverse SCHC model supports direct dilepton and "
    "quark-pair final states",
    "[gra::MPhotoZ][amplitude][SCHC]") {
  ModelParamRestoreGuard restore;
  gra::MODELPARAM = "TUNE0";

  auto require_finite = [](int pdg, double mass, std::size_t expected_hamps) {
    gra::LORENTZSCALAR lts = MakeToyPhotoZFFbar(pdg, mass);
    gra::MPhotoZ       photoz(lts, gra::MModelTune::Load(modelfile), gra::MPhotoZ::ProcessDefinitionFor());
    const double       amp2 = photoz.Amp2(lts);
    REQUIRE_FALSE(lts.proton_good_walker.has_value());
    REQUIRE(lts.exact_forward_photon_kinematics);
    REQUIRE(std::isfinite(amp2));
    REQUIRE(amp2 >= 0.0);
    REQUIRE(amp2 > 0.0);
    REQUIRE(lts.hamp.size() == expected_hamps);
    REQUIRE(std::isfinite(gra::SquaredNorm(lts.hamp)));
  };

  require_finite(11, 60.0, 16);
  require_finite(13, 91.1876, 16);
  require_finite(2, 91.1876, 48);
  require_finite(5, 91.1876, 48);

  gra::LORENTZSCALAR lts = MakeToyPhotoZFFbar(13, 91.1876);
  gra::RequireModelCache(lts.model_cache, gra::MModelTune::Load(modelfile), "test PhotoZ source");
  const auto emitter_state = gra::ResolveForwardLegState(lts, gra::ForwardBeamLeg::Upper);
  REQUIRE_THROWS(gra::qed::TransversePhotonSourceAmplitude(lts, emitter_state, 0));
}

TEST_CASE(
    "MPhotoZ keeps explicit proton-spin normalization and supports "
    "either dissociative role",
    "[gra::MPhotoZ][amplitude][dissociation]") {
  ModelParamRestoreGuard restore;
  gra::MODELPARAM = "TUNE0";

  auto evaluate = [](bool excite1, bool excite2) {
    gra::LORENTZSCALAR lts = MakeToyPhotoZFFbar(13);
    if (excite1) { SetToyPhotoForwardExcitation(lts, 1, 2.0); }
    if (excite2) { SetToyPhotoForwardExcitation(lts, 2, 2.0); }
    gra::MPhotoZ photoz(lts, gra::MModelTune::Load(modelfile), gra::MPhotoZ::ProcessDefinitionFor());
    const double amp2 = photoz.Amp2(lts);
    REQUIRE_FALSE(lts.proton_good_walker.has_value());
    const std::size_t direction_count = (excite1 || excite2) ? 2 : 1;
    const std::size_t spin_block      = 4 * direction_count;
    REQUIRE(lts.hamp.size() == 4 * spin_block);
    for (std::size_t spin = 1; spin < 4; ++spin) {
      for (std::size_t i = 0; i < spin_block; ++i) { REQUIRE(lts.hamp[spin * spin_block + i] == lts.hamp[i]); }
    }
    double direct_norm = 0.0;
    for (const auto &amplitude : lts.hamp) { direct_norm += std::norm(amplitude); }
    REQUIRE(amp2 == Approx(0.25 * direct_norm).epsilon(1e-12));
    REQUIRE(std::isfinite(amp2));
    REQUIRE(amp2 > 0.0);
    return amp2;
  };

  const double elastic        = evaluate(false, false);
  const double upper_excited  = evaluate(true, false);
  const double lower_excited  = evaluate(false, true);
  const double double_excited = evaluate(true, true);
  REQUIRE(elastic > 0.0);
  REQUIRE(upper_excited == Approx(lower_excited).epsilon(1e-10));
  REQUIRE(double_excited > 0.0);

  // Require the generated target mass in both directional thresholds
  auto require_closed_direction = [](int target_leg) {
    gra::LORENTZSCALAR lts = MakeToyPhotoZFFbar(13);
    SetToyPhotoForwardExcitation(lts, target_leg, 120.0);
    const int         photon_leg = target_leg == 1 ? 2 : 1;
    const gra::M4Vec &q_photon   = photon_leg == 1 ? lts.q1 : lts.q2;
    const gra::M4Vec &p_target   = target_leg == 1 ? lts.pbeam1 : lts.pbeam2;
    const double      W          = std::sqrt((q_photon + p_target).M2());
    REQUIRE(W < std::sqrt(gra::MPhotoQCD::CentralMass2(lts)) + lts.pfinal[static_cast<std::size_t>(target_leg)].M());
    if (photon_leg == 1) {
      lts.x2  = 0.0;
      lts.xi2 = 0.0;
    } else {
      lts.x1  = 0.0;
      lts.xi1 = 0.0;
    }
    gra::MPhotoZ photoz(lts, gra::MModelTune::Load(modelfile), gra::MPhotoZ::ProcessDefinitionFor());
    REQUIRE(photoz.Amp2(lts) == Approx(0.0));
  };
  require_closed_direction(1);
  require_closed_direction(2);
}

TEST_CASE("MPhotoZ remains finite under steered PhotoZ numerical kernels", "[gra::MPhotoZ][amplitude]") {
  ModelParamRestoreGuard restore;

  auto evaluate_with_tune = [](const std::pair<std::string, std::string> &tune) {
    gra::MODELPARAM        = tune.first;
    gra::LORENTZSCALAR lts = MakeToyPhotoZFFbar(13, 91.1876);
    gra::MPhotoZ       photoz(lts, gra::MModelTune::Load(tune.second), gra::MPhotoZ::ProcessDefinitionFor());
    const double       amp2 = photoz.Amp2(lts);
    REQUIRE_FALSE(lts.hamp.empty());
    REQUIRE(std::isfinite(amp2));
    REQUIRE(amp2 > 0.0);
    return amp2;
  };

  const auto coarse = WriteModifiedPhotoZTune(
      "coarse_kernel", [](auto &) {},
      [](auto &j) {
        j["NUMERICS_PHOTOZ"]["N_kappa"]       = 3;
        j["NUMERICS_PHOTOZ"]["N_k"]           = 3;
        j["NUMERICS_PHOTOZ"]["k2_MAX"]        = 30.0;
        j["NUMERICS_PHOTOZ"]["k2_MAX_use_q2"] = false;
      });
  const auto refined = WriteModifiedPhotoZTune(
      "refined_kernel", [](auto &) {},
      [](auto &j) {
        j["NUMERICS_PHOTOZ"]["N_kappa"]       = 5;
        j["NUMERICS_PHOTOZ"]["N_k"]           = 5;
        j["NUMERICS_PHOTOZ"]["k2_MAX"]        = 60.0;
        j["NUMERICS_PHOTOZ"]["k2_MAX_use_q2"] = false;
      });

  const double coarse_amp2  = evaluate_with_tune(coarse);
  const double refined_amp2 = evaluate_with_tune(refined);
  const double ratio        = refined_amp2 / coarse_amp2;
  REQUIRE(std::isfinite(ratio));
  REQUIRE(ratio > 1e-3);
  REQUIRE(ratio < 1e3);
}

TEST_CASE("MPhotoZ rejects root syntax and unsupported direct final states", "[gra::MPhotoZ][validation]") {
  ModelParamRestoreGuard restore;
  gra::MODELPARAM = "TUNE0";

  gra::LORENTZSCALAR root_lts = MakeToyPhotoZFFbar(13);
  gra::MDecayBranch  z_branch;
  z_branch.p         = root_lts.PDG.FindByPDG(23);
  z_branch.p4        = root_lts.decaytree[0].p4 + root_lts.decaytree[1].p4;
  z_branch.legs      = root_lts.decaytree;
  root_lts.decaytree = {z_branch};
  root_lts.process.root_decay_mode = gra::RootDecayMode::Physical;
  REQUIRE_THROWS_AS(gra::MPhotoZ(root_lts, gra::MModelTune::Load(modelfile), gra::MPhotoZ::ProcessDefinitionFor()),
                    std::invalid_argument);

  gra::LORENTZSCALAR mixed_lts = MakeToyPhotoZFFbar(11);
  mixed_lts.decaytree[1].p = mixed_lts.PDG.FindByPDG(-13);
  REQUIRE_THROWS_AS(gra::MPhotoZ(mixed_lts, gra::MModelTune::Load(modelfile), gra::MPhotoZ::ProcessDefinitionFor()),
                    std::invalid_argument);

  gra::LORENTZSCALAR top_lts = MakeToyPhotoZFFbar(6, 360.0);
  REQUIRE_THROWS_AS(gra::MPhotoZ(top_lts, gra::MModelTune::Load(modelfile), gra::MPhotoZ::ProcessDefinitionFor()),
                    std::invalid_argument);

  gra::MSubProc subproc({"ygg"}, "F");
  subproc.Initialize("ygg", "Z");
  REQUIRE(subproc.RootResonancePDG(root_lts) == 0);
}

TEST_CASE("MPhotoZ assigns HepMC-compatible color flow through MSubProc", "[gra::MPhotoZ][color]") {
  ModelParamRestoreGuard restore;
  gra::MODELPARAM = "TUNE0";

  gra::LORENTZSCALAR quark_lts = MakeToyPhotoZFFbar(2);
  ToyHelicityProcess master;
  master.state.lts = quark_lts;
  master.ProcPtr   = gra::MSubProc({"ygg"}, "F");
  master.ProcPtr.Initialize("ygg", "Z");
  master.InitializeProcessAmplitude();
  gra::MSubProc subproc = master.ProcPtr;
  quark_lts             = master.state.lts;
  const double amp2     = subproc.GetBareAmplitude2(quark_lts);
  REQUIRE(std::isfinite(amp2));
  REQUIRE(quark_lts.decaytree[0].p.color_flow.empty());
  subproc.SampleColorFlow(quark_lts);
  REQUIRE(quark_lts.decaytree[0].p.color_flow.flow1 == 501);
  REQUIRE(quark_lts.decaytree[0].p.color_flow.flow2 == 0);
  REQUIRE(quark_lts.decaytree[1].p.color_flow.flow1 == 0);
  REQUIRE(quark_lts.decaytree[1].p.color_flow.flow2 == 501);

  std::map<int, int> balance;
  for (const auto &branch : quark_lts.decaytree) {
    if (branch.p.color_flow.flow1 != 0) { ++balance[branch.p.color_flow.flow1]; }
    if (branch.p.color_flow.flow2 != 0) { --balance[branch.p.color_flow.flow2]; }
  }
  REQUIRE(balance[501] == 0);

  ToyPhotoZRecordProcess recorder;
  recorder.SetModelTune(gra::MModelTune::Load(modelfile));
  recorder.state.lts = quark_lts;
  HepMC3::GenEvent evt(HepMC3::Units::GEV, HepMC3::Units::MM);
  REQUIRE(recorder.Record(evt));
  bool found_quark     = false;
  bool found_antiquark = false;
  for (const auto &particle : evt.particles()) {
    if (particle->pid() == 2) {
      found_quark = true;
      REQUIRE(ReadHepMCIntAttributeForTest(particle, "flow1") == 501);
      REQUIRE(ReadHepMCIntAttributeForTest(particle, "flow2") == 0);
    }
    if (particle->pid() == -2) {
      found_antiquark = true;
      REQUIRE(ReadHepMCIntAttributeForTest(particle, "flow1") == 0);
      REQUIRE(ReadHepMCIntAttributeForTest(particle, "flow2") == 501);
    }
  }
  REQUIRE(found_quark);
  REQUIRE(found_antiquark);

  gra::LORENTZSCALAR lepton_lts = MakeToyPhotoZFFbar(13);
  ToyHelicityProcess lepton_master;
  lepton_master.state.lts = lepton_lts;
  lepton_master.ProcPtr = gra::MSubProc({"ygg"}, "F");
  lepton_master.ProcPtr.Initialize("ygg", "Z");
  lepton_master.InitializeProcessAmplitude();
  subproc = lepton_master.ProcPtr;
  lepton_lts = lepton_master.state.lts;
  REQUIRE(std::isfinite(subproc.GetBareAmplitude2(lepton_lts)));
  subproc.SampleColorFlow(lepton_lts);
  REQUIRE(lepton_lts.decaytree[0].p.color_flow.empty());
  REQUIRE(lepton_lts.decaytree[1].p.color_flow.empty());
}

TEST_CASE("Exact EPA photon records keep proton beams in the LHE bridge view", "[gra::MPhotoZ][gra::MLHE]") {
  gra::LORENTZSCALAR lts              = MakeToyPhotoZFFbar(13);
  lts.id1                             = gra::PDG::PDG_gamma;
  lts.id2                             = gra::PDG::PDG_gamma;
  lts.muF                             = 91.1876;
  lts.muR                             = 80.0;
  lts.scalup                          = 70.0;
  lts.pdf_xf1                         = 1.0;
  lts.pdf_xf2                         = 1.0;
  lts.exact_forward_photon_kinematics = true;

  ToyPhotoZRecordProcess recorder;
  recorder.SetModelTune(gra::MModelTune::Load(modelfile));
  recorder.state.lts = lts;
  HepMC3::GenEvent evt(HepMC3::Units::GEV, HepMC3::Units::MM);
  REQUIRE(recorder.Record(evt));
  REQUIRE(evt.attribute<HepMC3::IntAttribute>("graniitti_exact_forward_photon_kinematics")->value() == 1);
  REQUIRE(evt.attribute<HepMC3::DoubleAttribute>("graniitti_mu_f")->value() == Approx(lts.muF));
  REQUIRE(evt.attribute<HepMC3::DoubleAttribute>("graniitti_mu_r")->value() == Approx(lts.muR));
  REQUIRE(evt.attribute<HepMC3::DoubleAttribute>("graniitti_scalup")->value() == Approx(lts.scalup));
  REQUIRE(evt.pdf_info() != nullptr);
  REQUIRE(evt.pdf_info()->parton_id[0] == gra::PDG::PDG_gamma);
  REQUIRE(evt.pdf_info()->parton_id[1] == gra::PDG::PDG_gamma);

  const std::vector<gra::MLHERow> rows = gra::BuildLHERows(evt);
  REQUIRE(rows.size() >= 6);
  REQUIRE(rows[0].pid == 2212);
  REQUIRE(rows[1].pid == 2212);
  REQUIRE(rows[0].status == -1);
  REQUIRE(rows[1].status == -1);
  int stable_protons = 0;
  for (const auto &row : rows) {
    if (row.pid == 2212 && row.status == 1) { ++stable_protons; }
  }
  REQUIRE(stable_protons == 2);
  const std::string diff_comment = gra::BuildLHEDiffractionComment(evt);
  REQUIRE(diff_comment.find("side_mask=3") != std::string::npos);

  const std::string lhe_path      = "tmp/graniitti_exact_photo_lhe_test.lhe";
  auto              cross_section = std::make_shared<HepMC3::GenCrossSection>();
  cross_section->set_cross_section(1.0, 0.0);
  cross_section->set_accepted_events(1);
  cross_section->set_attempted_events(1);
  evt.set_cross_section(cross_section);
  evt.weights() = {1.0};
  {
    gra::MLHEWriter writer(lhe_path);
    writer.WriteEvent(evt);
  }
  std::ifstream input(lhe_path);
  REQUIRE(input.good());
  const std::string lhe_text((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
  REQUIRE(lhe_text.find(" 2212 -1 ") != std::string::npos);
  REQUIRE(lhe_text.find("side_mask=3") != std::string::npos);
  REQUIRE(lhe_text.find("<scales") != std::string::npos);
  REQUIRE(lhe_text.find("muf=") != std::string::npos);
  REQUIRE(lhe_text.find("mur=") != std::string::npos);
  LHEF::Reader unit_reader(lhe_path);
  REQUIRE(unit_reader.heprup.IDWTUP == 3);
  REQUIRE(unit_reader.heprup.XSECUP.size() == 1);
  REQUIRE(unit_reader.heprup.XSECUP[0] == Approx(1.0));
  REQUIRE(unit_reader.heprup.XMAXUP[0] == Approx(1.0));
  REQUIRE(unit_reader.readEvent());
  REQUIRE(unit_reader.hepeup.IDPRUP == 1);
  REQUIRE(unit_reader.hepeup.XWGTUP == Approx(1.0));

  const std::string weighted_lhe_path = "tmp/graniitti_exact_photo_weighted_lhe_test.lhe";
  evt.weights()                       = {2.0};
  {
    gra::MLHERunConfig config;
    config.weight_type  = gra::MLHEWeightType::Weighted;
    config.weight_scale = 2.5;
    gra::MLHEWriter writer(weighted_lhe_path, config);
    writer.WriteEvent(evt);
  }
  LHEF::Reader weighted_reader(weighted_lhe_path);
  REQUIRE(weighted_reader.heprup.IDWTUP == 4);
  REQUIRE(weighted_reader.heprup.XSECUP[0] == Approx(1.0));
  REQUIRE(weighted_reader.heprup.XMAXUP[0] == Approx(0.0));
  REQUIRE(weighted_reader.readEvent());
  REQUIRE(weighted_reader.hepeup.XWGTUP == Approx(5.0));

  const std::string signed_lhe_path = "tmp/graniitti_exact_photo_signed_lhe_test.lhe";
  evt.weights()                     = {-2.0};
  {
    gra::MLHERunConfig config;
    config.weight_type  = gra::MLHEWeightType::SignedWeighted;
    config.weight_scale = 2.5;
    gra::MLHEWriter writer(signed_lhe_path, config);
    writer.WriteEvent(evt);
  }
  LHEF::Reader signed_reader(signed_lhe_path);
  REQUIRE(signed_reader.heprup.IDWTUP == -4);
  REQUIRE(signed_reader.readEvent());
  REQUIRE(signed_reader.hepeup.XWGTUP == Approx(-5.0));

  const std::string invalid_signed_lhe_path = "tmp/graniitti_exact_photo_invalid_signed_lhe_test.lhe";
  {
    gra::MLHERunConfig config;
    config.weight_type = gra::MLHEWeightType::Weighted;
    gra::MLHEWriter writer(invalid_signed_lhe_path, config);
    REQUIRE_THROWS_AS(writer.WriteEvent(evt), std::invalid_argument);
  }

  const std::string invalid_lhe_path = "tmp/graniitti_exact_photo_invalid_lhe_test.lhe";
  evt.weights()                      = {2.0};
  {
    gra::MLHEWriter writer(invalid_lhe_path);
    REQUIRE_THROWS_AS(writer.WriteEvent(evt), std::invalid_argument);
  }
  evt.weights() = {1.0};

  const std::string hepmc_path         = "tmp/graniitti_exact_photo_weighted_lhe_test.hepmc3";
  const std::string converted_lhe_path = "tmp/graniitti_exact_photo_converted_lhe_test.lhe";
  {
    HepMC3::WriterAscii writer(hepmc_path);
    writer.write_event(evt);
    evt.weights() = {3.0};
    writer.write_event(evt);
    writer.close();
  }
  const gra::MLHEConversionStats conversion = gra::ConvertHepMC3ToLHE(hepmc_path, converted_lhe_path, false);
  REQUIRE(conversion.events == 2);
  LHEF::Reader converted_reader(converted_lhe_path);
  REQUIRE(converted_reader.heprup.IDWTUP == 4);
  REQUIRE(converted_reader.readEvent());
  REQUIRE(converted_reader.hepeup.XWGTUP == Approx(0.5));
  REQUIRE(converted_reader.readEvent());
  REQUIRE(converted_reader.hepeup.XWGTUP == Approx(1.5));
  evt.weights() = {1.0};

  const std::string changed_run_path = "tmp/graniitti_exact_photo_changed_run_lhe_test.lhe";
  {
    gra::MLHEWriter writer(changed_run_path);
    writer.WriteEvent(evt);
    cross_section->set_cross_section(2.0, 0.0);
    REQUIRE_THROWS_AS(writer.WriteEvent(evt), std::invalid_argument);
  }
  cross_section->set_cross_section(1.0, 0.0);

  gra::LORENTZSCALAR dissociative_lts = lts;
  SetToyPhotoForwardExcitation(dissociative_lts, 1, 2.0);
  SetToyPhotoForwardExcitation(dissociative_lts, 2, 2.0);
  ToyPhotoZRecordProcess dissociative_recorder;
  dissociative_recorder.SetModelTune(gra::MModelTune::Load(modelfile));
  dissociative_recorder.state.lts = dissociative_lts;
  HepMC3::GenEvent dissociative_evt(HepMC3::Units::GEV, HepMC3::Units::MM);
  REQUIRE(dissociative_recorder.Record(dissociative_evt));
  const std::string dissociative_comment = gra::BuildLHEDiffractionComment(dissociative_evt);
  CAPTURE(dissociative_comment);
  REQUIRE(dissociative_comment.find("side1_excited=1") != std::string::npos);
  REQUIRE(dissociative_comment.find("side2_excited=1") != std::string::npos);

  auto dissociative_cross_section = std::make_shared<HepMC3::GenCrossSection>();
  dissociative_cross_section->set_cross_section(1.0, 0.0);
  dissociative_cross_section->set_accepted_events(1);
  dissociative_cross_section->set_attempted_events(1);
  dissociative_evt.set_cross_section(dissociative_cross_section);
  dissociative_evt.weights()              = {1.0};
  const std::string dissociative_lhe_path = "tmp/graniitti_exact_photo_nstar_lhe_test.lhe";
  {
    gra::MLHEWriter writer(dissociative_lhe_path);
    writer.WriteEvent(dissociative_evt);
  }
  std::ifstream dissociative_input(dissociative_lhe_path);
  REQUIRE(dissociative_input.good());

  gra::LORENTZSCALAR colored_lts              = MakeToyPhotoZFFbar(2);
  colored_lts.id1                             = gra::PDG::PDG_gamma;
  colored_lts.id2                             = gra::PDG::PDG_gamma;
  colored_lts.muF                             = 91.1876;
  colored_lts.muR                             = 80.0;
  colored_lts.scalup                          = 70.0;
  colored_lts.pdf_xf1                         = 1.0;
  colored_lts.pdf_xf2                         = 1.0;
  colored_lts.exact_forward_photon_kinematics = true;
  SetToyPhotoForwardExcitation(colored_lts, 1, 2.0);
  SetToyPhotoForwardExcitation(colored_lts, 2, 2.0);
  colored_lts.decaytree[0].p.color_flow.flow1 = 501;
  colored_lts.decaytree[1].p.color_flow.flow2 = 501;
  ToyPhotoZRecordProcess colored_recorder;
  colored_recorder.SetModelTune(gra::MModelTune::Load(modelfile));
  colored_recorder.state.lts = colored_lts;
  HepMC3::GenEvent colored_evt(HepMC3::Units::GEV, HepMC3::Units::MM);
  REQUIRE(colored_recorder.Record(colored_evt));
  auto colored_cross_section = std::make_shared<HepMC3::GenCrossSection>();
  colored_cross_section->set_cross_section(1.0, 0.0);
  colored_cross_section->set_accepted_events(1);
  colored_cross_section->set_attempted_events(1);
  colored_evt.set_cross_section(colored_cross_section);
  colored_evt.weights()              = {1.0};
  const std::string colored_lhe_path = "tmp/graniitti_exact_photo_colored_nstar_lhe_test.lhe";
  {
    gra::MLHEWriter writer(colored_lhe_path);
    writer.WriteEvent(colored_evt);
  }
  std::ifstream colored_input(colored_lhe_path);
  REQUIRE(colored_input.good());
}

TEST_CASE("MPhotoZ remains compatible with the existing screening loop", "[gra::MPhotoZ][screening]") {
  ModelParamRestoreGuard model_restore;
  gra::MODELPARAM = "TUNE0";

  const gra::LORENTZSCALAR event   = MakeToyPhotoZFFbar(13);
  gra::MEikonal            eikonal = BuildTestEikonalWithInitialState({}, "screen_photoz", {event.beam1, event.beam2},
                                                                                0.25, 1, 4, 4, event.s);
  eikonal.Numerics.LOOP.radial_integrator  = "1/3";
  eikonal.Numerics.LOOP.azimuth_integrator = "Trap";
  eikonal.Numerics.LOOP.r_min              = 0.0;
  eikonal.Numerics.LOOP.r_max              = 0.1;
  eikonal.Numerics.LOOP.radial_intervals   = 2;
  eikonal.Numerics.LOOP.azimuth_nodes      = 2;
  eikonal.InitLoopWeightMatrix();
  ToyPhotoZScreeningProcess process(eikonal.ModelTuneHandle());
  process.SetEikonal(eikonal);

  const double bare = process.BareAmp2();
  REQUIRE(std::isfinite(bare));
  REQUIRE_FALSE(process.state.lts.hamp.empty());
  REQUIRE_FALSE(process.state.lts.proton_good_walker.has_value());
  const double screened = process.ScreenedAmp2();
  REQUIRE(std::isfinite(screened));
  REQUIRE(screened >= 0.0);
  REQUIRE_FALSE(process.state.lts.hamp.empty());
  REQUIRE_FALSE(process.state.lts.proton_good_walker.has_value());
}

TEST_CASE("MPhotoVMNumerics rejects unsupported physics lookups", "[gra::MPhotoVM][validation]") {
  ModelParamRestoreGuard restore;
  gra::MODELPARAM = "TUNE0";

  gra::MPhotoVMNumerics numerics;
  const std::string     numerics_file = gra::ResolveModelDataFile("TUNE0", "NUMERICS.json");
  REQUIRE_NOTHROW(numerics.ConfigureFromJson(modelfile, gra::aux::GetInputData(modelfile), numerics_file,
                                             gra::aux::GetInputData(numerics_file)));
  REQUIRE(numerics.initialized);
  REQUIRE_THROWS_AS(numerics.HeavyQuarkMass(6), std::invalid_argument);
  REQUIRE_THROWS_AS(numerics.Channel(999999), std::invalid_argument);
}

TEST_CASE("PhotoVM parameters construct safely from one tune across threads", "[gra::MPhotoVM][threading]") {
  ModelParamRestoreGuard restore;
  gra::MODELPARAM = "TUNE0";

  const auto       model_tune = gra::MModelTune::Load(modelfile);
  gra::MModelCache cache(model_tune);
  const auto       first  = gra::GetPhotoVMNumerics(cache);
  const auto       second = gra::GetPhotoVMNumerics(cache);
  REQUIRE(first == second);
  REQUIRE(first->initialized);
  REQUIRE(second->initialized);

  constexpr std::size_t                                     nthreads = 8;
  std::vector<std::shared_ptr<const gra::MPhotoVMNumerics>> handles(nthreads);
  std::vector<std::thread>                                  workers;
  workers.reserve(nthreads);

  for (std::size_t i = 0; i < nthreads; ++i) {
    workers.emplace_back([i, &handles, &cache] { handles[i] = gra::GetPhotoVMNumerics(cache); });
  }
  for (auto &worker : workers) { worker.join(); }

  for (const auto &handle : handles) { REQUIRE(handle == first); }
}

TEST_CASE("PhotoVM reader uses changed NUMERICS content", "[gra::MPhotoVM][threading]") {
  ModelParamRestoreGuard restore;

  const auto first_tune = WriteModifiedPhotoZTune(
      "photovm_cache_reload", [](auto &) {}, [](auto &j) { j["NUMERICS_PHOTOVM"]["jmrt_lambda_step"] = 0.06; });
  gra::MODELPARAM        = first_tune.first;
  const auto first_model = gra::MModelTune::Load(first_tune.second);
  const auto first       = gra::ReadPhotoVMNumerics(*first_model);
  REQUIRE(first->jmrt_lambda_step == Approx(0.06));

  const auto second_tune = WriteModifiedPhotoZTune(
      "photovm_cache_reload", [](auto &) {}, [](auto &j) { j["NUMERICS_PHOTOVM"]["jmrt_lambda_step"] = 0.09; });
  gra::MODELPARAM         = second_tune.first;
  const auto second_model = gra::MModelTune::Load(second_tune.second);
  REQUIRE(second_model != first_model);
  REQUIRE(second_model->Numerics() != first_model->Numerics());
  const auto second = gra::ReadPhotoVMNumerics(*second_model);
  REQUIRE(second->jmrt_lambda_step == Approx(0.09));
}

TEST_CASE("PhotoVM reader uses changed GENERAL physics content", "[gra::MPhotoVM][threading]") {
  ModelParamRestoreGuard restore;

  const auto first_tune =
      WriteModifiedPhotoVMTune("physics_cache_reload", [](auto &j) { j["PARAM_PHOTOVM"]["channels"][0][3] = 0.95; });
  gra::MODELPARAM        = first_tune.first;
  const auto first_model = gra::MModelTune::Load(first_tune.second);
  const auto first       = gra::ReadPhotoVMNumerics(*first_model);
  REQUIRE(first->Channel(443).wavefunction_correction == Approx(0.95));

  const auto second_tune =
      WriteModifiedPhotoVMTune("physics_cache_reload", [](auto &j) { j["PARAM_PHOTOVM"]["channels"][0][3] = 0.90; });
  gra::MODELPARAM         = second_tune.first;
  const auto second_model = gra::MModelTune::Load(second_tune.second);
  REQUIRE(second_model->GeneralFile() == first_model->GeneralFile());
  REQUIRE(second_model != first_model);
  const auto second = gra::ReadPhotoVMNumerics(*second_model);
  REQUIRE(second->Channel(443).wavefunction_correction == Approx(0.90));
}

TEST_CASE("MPhotoVM gamma-p cross sections have consistent normalization",
          "[gra::MPhotoVM][physics][literature]") {
  ModelParamRestoreGuard restore;
  gra::MODELPARAM = "TUNE0";

  gra::LORENTZSCALAR lts = MakeToyPhotoZFFbar(13, 3.096900);
  gra::MPhotoVM      jpsi(lts, gra::MModelTune::Load(modelfile), "jpsi", gra::MPhotoVM::ProcessDefinitionFor("jpsi"));
  gra::MPhotoVM psi2s(lts, gra::MModelTune::Load(modelfile), "psi(2S)", gra::MPhotoVM::ProcessDefinitionFor("psi(2S)"));
  const double  jpsi_threshold = lts.PDG.FindByPDG(443).mass + gra::PDG::mp;
  REQUIRE(jpsi.GammaPCrossSection(lts, jpsi_threshold - 1.0e-6, 1.2) == Approx(0.0));
  REQUIRE(jpsi.GammaPDSigmaDt(lts, jpsi_threshold - 1.0e-6, 0.0) == Approx(0.0));

  // [REFERENCE: https://arxiv.org/abs/1307.7099]
  REQUIRE(jpsi.GammaPTSlope(90.0) == Approx(4.9).epsilon(1e-14));
  REQUIRE(jpsi.GammaPTSlope(180.0) == Approx(4.9 + 4.0 * 0.06 * std::log(2.0)).epsilon(1e-14));

  // The exponential differential spectrum must integrate with the same measured
  // slope
  const double forward   = jpsi.GammaPDSigmaDt(lts, 90.0, 0.0);
  const double half_gev2 = jpsi.GammaPDSigmaDt(lts, 90.0, -0.5);
  REQUIRE(forward > 0.0);
  REQUIRE(half_gev2 / forward == Approx(std::exp(-0.5 * 4.9)).epsilon(1e-12));

  // [REFERENCE: https://arxiv.org/abs/1304.5162]
  // Elastic gamma-p data integrated over |t|<1.2 GeV2
  const auto h1_integrated_nb = H1PhotoData();
  for (const auto &[W, measured_nb] : h1_integrated_nb) {
    CAPTURE(W, measured_nb);
    // This broad envelope for rounded published fit coefficients does not assert fiducial data agreement
    REQUIRE(jpsi.GammaPCrossSection(lts, W, 1.2) == Approx(measured_nb).epsilon(0.35));
  }
  // Integrate the differential prediction at the same W as its gamma-p cross section
  // The H1 t spectrum is averaged over W and cannot be compared at W=90 GeV
  const double slope = jpsi.GammaPTSlope(90.0);
  REQUIRE(jpsi.GammaPCrossSection(lts, 90.0, 1.2) ==
          Approx(forward * (-std::expm1(-slope * 1.2)) / slope).epsilon(1e-12));

  // The electronic-width normalization contains photon-pole alpha(0), not
  // running alpha(MV^2)
  const auto running_tune = WriteModifiedPhotoVMTune(
      "running_qed_photon_pole", [](auto &general) { general["PARAM_STRUCTURE"]["QED_alpha"] = "LL"; });
  gra::LORENTZSCALAR running_lts = MakeToyPhotoZFFbar(13, 3.096900);
  gra::MPhotoVM      running_jpsi(running_lts, gra::MModelTune::Load(running_tune.second), "jpsi",
                                  gra::MPhotoVM::ProcessDefinitionFor("jpsi"));
  const double       running_steered_forward = running_jpsi.GammaPDSigmaDt(running_lts, 90.0, 0.0);
  REQUIRE(running_steered_forward == Approx(forward).epsilon(1e-14));

  // [REFERENCE: https://arxiv.org/abs/2206.13343]
  // Fiducial range 30<W<180 GeV, Q2<1 GeV2 and |t|<1 GeV2
  // R=0.146 +/-0.010(stat) +0.016/-0.020(syst), with no significant W
  // dependence
  for (const double W : {50.0, 90.0, 150.0}) {
    const double sigma_jpsi  = jpsi.GammaPCrossSection(lts, W, 1.0);
    const double sigma_psi2s = psi2s.GammaPCrossSection(lts, W, 1.0);
    CAPTURE(W, sigma_jpsi, sigma_psi2s);
    REQUIRE(sigma_jpsi > 0.0);
    REQUIRE(sigma_psi2s / sigma_jpsi == Approx(0.146).margin(0.026));
  }

  // ZEUS used b(J/psi)=4.6+/-0.3 GeV^-2 and b(psi(2S))=4.3+/-0.7 GeV^-2
  REQUIRE(jpsi.GammaPTSlope(90.0) == Approx(4.6).margin(0.6));
  REQUIRE(psi2s.GammaPTSlope(90.0) == Approx(4.3).margin(0.7));
}

TEST_CASE("MPhotoVM implements the published JMRT 2013 fitted gluons", "[gra::MPhotoVM][physics][literature]") {
  ModelParamRestoreGuard                         restore;
  const auto h1_integrated_nb = H1PhotoData();

  auto fitted_prediction = [&](const std::string &source, const std::string &suffix, const std::string &pdf) {
    const auto tune = WriteModifiedPhotoVMTune(suffix, [&](auto &j) { j["PARAM_PHOTOVM"]["gluon_source"] = source; });
    gra::MODELPARAM = tune.first;
    gra::LORENTZSCALAR lts = MakeToyPhotoZFFbar(13, 3.096900);
    lts.LHAPDFSET          = pdf;
    gra::MPhotoVM jpsi(lts, gra::MModelTune::Load(tune.second), "jpsi", gra::MPhotoVM::ProcessDefinitionFor("jpsi"));

    std::vector<double> out;
    for (const auto &[W, measured_nb] : h1_integrated_nb) {
      const double sigma = jpsi.GammaPCrossSection(lts, W, 1.2);
      CAPTURE(source, pdf, W, sigma, measured_nb);
      REQUIRE(std::isfinite(sigma));
      REQUIRE(sigma > 0.0);
      REQUIRE(sigma == Approx(measured_nb).epsilon(0.35));
      out.push_back(sigma);
    }
    return out;
  };

  // [REFERENCE: S.P. Jones et al., arXiv:1307.7099, Eqs. (1), (4), (5) and (16)]
  const auto lo = fitted_prediction("JMRT_2013_LO", "jmrt_2013_lo", "CT10nlo");

  // [REFERENCE: S.P. Jones et al., arXiv:1307.7099, Eqs. (4), (5), (7), (9)--(11) and (17)]
  const auto nlo_ct10  = fitted_prediction("JMRT_2013_NLO", "jmrt_2013_nlo_ct10", "CT10nlo");
  const auto nlo_nnpdf = fitted_prediction("JMRT_2013_NLO", "jmrt_2013_nlo_nnpdf", "NNPDF31_lo_as_0118");

  REQUIRE(lo != nlo_ct10);
  REQUIRE(nlo_ct10.size() == nlo_nnpdf.size());
  for (std::size_t i = 0; i < nlo_ct10.size(); ++i) { REQUIRE(nlo_ct10[i] == Approx(nlo_nnpdf[i]).epsilon(1e-13)); }
}

// Resolve the integrable derivative singularity of the published NLO gluon at Q0
TEST_CASE("JMRT fitted NLO quadrature converges", "[gra::MPhotoVM][physics][convergence]") {
  ModelParamRestoreGuard restore;
  std::vector<double>    previous;
  for (const unsigned int n : {24U, 48U, 96U}) {
    const auto tune = WriteModifiedPhotoZTune(
        "jmrt_convergence_" + std::to_string(n), [](auto &j) { j["PARAM_PHOTOVM"]["gluon_source"] = "JMRT_2013_NLO"; },
        [n](auto &j) { j["NUMERICS_PHOTOVM"]["N_k"] = n; });
    auto          lts = MakeToyPhotoZFFbar(13, 3.0969);
    gra::MPhotoVM vm(lts, gra::MModelTune::Load(tune.second), "jpsi", gra::MPhotoVM::ProcessDefinitionFor("jpsi"));
    std::vector<double> current;
    for (const double w : {30.0, 90.0, 1000.0}) { current.push_back(vm.GammaPCrossSection(lts, w, 1.2)); }
    if (!previous.empty()) {
      for (const auto &i : indices(current)) {
        CAPTURE(n, i, current[i], previous[i]);
        CHECK(current[i] == Approx(previous[i]).epsilon(0.002));
      }
    }
    previous = current;
  }
}

// Compare complex amplitudes with independent fitted and signed LHAPDF gluon integrals
TEST_CASE("JMRT complex amplitudes reproduce independent gluon integrals", "[gra::MPhotoVM][physics][phase][signed]") {
  ModelParamRestoreGuard restore;
  const bool fitted = GENERATE(false, true);
  const int pdg = fitted ? 443 : 553;
  const auto             files = WriteModifiedPhotoZTune(
                  "jmrt_signed",
                  [fitted](auto &j) {
        j["PARAM_PHOTOVM"]["gluon_source"] = fitted ? "JMRT_2013_NLO" : "LHAPDF_SHUVAEV";
        j["PARAM_SKEWED_UGD"]["mode"]      = "PERTURBATIVE_ONLY";
      },
                  [](auto &j) {
        j["NUMERICS_SUDAKOV"]["SHUV"]["N"] = {40, 40};
        j["NUMERICS_SUDAKOV"]["SUDA"]["N"] = {40, 40};
        j["NUMERICS_PHOTOVM"]["N_k"]       = 256;
      });
  const auto tune     = gra::MModelTune::Load(files.second);
  const auto num      = gra::ReadPhotoVMNumerics(*tune);
  bool       negative = false;
  for (const double beam_energy : fitted ? std::vector<double>{200.0, 650.0, 6500.0} : std::vector<double>{6.5, 8.0, 12.0}) {
    auto         lts      = MakeToyPhotoZFFbar(13, LoadedPDGTable().FindByPDG(pdg).mass);
    const auto  &particle = lts.PDG.FindByPDG(pdg);
    const double pz       = std::sqrt(pow2(beam_energy) - pow2(gra::PDG::mp));
    const double energy   = beam_energy - particle.mass / 2.0;
    const double recoil   = std::sqrt(pow2(energy) - pow2(gra::PDG::mp) - lts.pfinal[1].Pt2());
    lts.pbeam1.SetPzE(pz, beam_energy);
    lts.pbeam2.SetPzE(-pz, beam_energy);
    lts.pfinal[1].SetPzE(recoil, energy);
    lts.pfinal[2].SetPzE(-recoil, energy);
    SetPhotoAntiprotonBeam(lts, 2);
    RefreshToyDerivedKinematicsPreserveDecay(lts);
    const std::string channel = fitted ? "jpsi" : "Upsilon(1S)";
    gra::MPhotoVM vm(lts, tune, channel, gra::MPhotoVM::ProcessDefinitionFor(channel));
    REQUIRE(vm.Amp2(lts) > 0.0);
    const double w2     = (lts.q1 + lts.pbeam2).M2();
    const double qbar2  = pow2(num->HeavyQuarkMass(num->Channel(pdg).quark_pdg));
    const double x      = 4.0 * qbar2 / w2;
    const double kmax = fitted ? (w2 - 4.0 * qbar2) / 4.0 : std::min(num->jmrt_k2_max, (w2 - 4.0 * qbar2) / 4.0);
    const auto   kernel = [&](double fraction) {
      if (fitted) { return JMRTSkewedKernel(*num, fraction, qbar2, kmax); }
      return gra::math::LogMeasureGaussIntegral(256, lts.GlobalSudakovPtr->GetQ2Min(), kmax, [&](double k2) {
        return lts.GlobalSudakovPtr->AlphaSFlux_xQ2Mu(
                     fraction, k2, std::max(std::sqrt(std::max(k2, qbar2)), lts.GlobalSudakovPtr->GetMuMin()),
                     std::max(k2, qbar2)) /
               (qbar2 * (qbar2 + k2));
      });
    };
    const double k = kernel(x);
    negative       = negative || k < 0.0;
    const double h = num->jmrt_lambda_step;
    const double rho =
        gra::math::PI / (4.0 * h) * std::log(std::abs(kernel(x * std::exp(-h)) / kernel(x * std::exp(h))));
    const double width   = num->LeptonicWidthEE(pdg);
    const auto   forward = w2 *
                         std::sqrt(width * gra::math::pow3(particle.mass) * gra::math::pow4(gra::math::PI) /
                                   (3.0 * gra::qed::alpha_QED())) *
                         k * std::complex<double>(rho, 1.0) * num->Channel(pdg).wavefunction_correction;
    const auto        prop   = 1.0 / std::complex<double>(gra::MPhotoQCD::CentralMass2(lts) - pow2(particle.mass),
                                                 particle.mass * particle.width);
    const auto        states = gra::MPhotoQCD::TransverseStates(lts, lts.pfinal[0]);
    const gra::MDirac dirac("DIRAC");
    std::size_t       row = 0;
    for (const int hf : {-1, 1}) {
      for (const int ha : {-1, 1}) {
        const double g       = std::sqrt(12.0 * gra::math::PI * width / particle.mass);
        const auto   current = gra::qed::FFVCurrent(dirac, lts.decaytree[0].p4, lts.decaytree[1].p4, hf, ha, g, g);
        std::complex<double> expected = 0.0;
        for (const int lambda : {-1, 1}) {
          for (const auto leg : {gra::ForwardBeamLeg::Upper, gra::ForwardBeamLeg::Lower}) {
            const auto source = gra::flux::PhotoSourceAmplitudes(lts, gra::ResolveForwardLegState(lts, leg), lambda);
            expected +=
                source[0] * forward * std::exp(vm.GammaPTSlope(std::sqrt(w2)) * lts.t1 / 2.0) * prop *
                gra::MPhotoQCD::ContractCurrentPolarization(current, states[gra::spin::BinaryHelicityIndexX2(lambda)]);
          }
        }
        CAPTURE(beam_energy, k, row, expected, lts.hamp[row]);
        CHECK(std::abs(lts.hamp[row] - expected) < (fitted ? 2.0e-4 : 1.0e-9) * std::sqrt(gra::SquaredNorm(lts.hamp)));
        ++row;
      }
    }
  }
  if (!fitted) { REQUIRE(negative); }
}

// Reject overflowing JMRT contributions during HERA initialization or event sampling
TEST_CASE("MPhotoVM rejects overflowing JMRT contributions", "[gra::MPhotoVM][regression]") {
  ModelParamRestoreGuard restore;
  const auto source = GENERATE("JMRT_2013_LO", "JMRT_2013_NLO");
  const std::string profile = GENERATE("soft", "hera");
  const auto tune = WriteModifiedPhotoVMTune(std::string("overflow_") + source + profile, [&](auto &j) {
    j["PARAM_NSTAR"]["MODEL"]["ygg"]["[22,P]"] = profile;
    j["PARAM_PHOTOVM"]["gluon_source"] = source;
    j["PARAM_PHOTOVM"]["jmrt_2013"]["lo_fit"][0] = std::numeric_limits<double>::max();
    j["PARAM_PHOTOVM"]["jmrt_2013"]["nlo_fit"][0] = std::numeric_limits<double>::max();
  });
  gra::MODELPARAM = tune.first;
  gra::LORENTZSCALAR lts = MakeToyPhotoZFFbar(13, 3.096900);
  if (profile == "hera") {
    REQUIRE_THROWS_AS(gra::MPhotoVM(lts, gra::MModelTune::Load(tune.second), "jpsi", gra::MPhotoVM::ProcessDefinitionFor("jpsi")),
                      gra::AmplitudeFailure);
    return;
  }
  gra::MPhotoVM jpsi(lts, gra::MModelTune::Load(tune.second), "jpsi", gra::MPhotoVM::ProcessDefinitionFor("jpsi"));
  REQUIRE_THROWS_AS(jpsi.GammaPDSigmaDt(lts, 90.0, 0.0), gra::AmplitudeFailure);
  REQUIRE_THROWS_AS(jpsi.Amp2(lts), gra::AmplitudeFailure);
}

TEST_CASE("MPhotoVM perturbative-only real part remains finite", "[gra::MPhotoVM][physics][regression]") {
  ModelParamRestoreGuard restore;
  const auto             tune = WriteModifiedPhotoVMTune("perturbative_real_part",
                                                         [](auto &j) { j["PARAM_SKEWED_UGD"]["mode"] = "PERTURBATIVE_ONLY"; });
  gra::MODELPARAM             = tune.first;
  gra::LORENTZSCALAR lts      = MakeToyPhotoZFFbar(13, 3.096900);
  gra::MPhotoVM      jpsi(lts, gra::MModelTune::Load(tune.second), "jpsi", gra::MPhotoVM::ProcessDefinitionFor("jpsi"));
  const double       sigma = jpsi.GammaPCrossSection(lts, 700.0, 1.2);
  REQUIRE(std::isfinite(sigma));
  REQUIRE(sigma > 0.0);
  REQUIRE(sigma < 1.0e5);
}

TEST_CASE("MPhotoVM supports heavy vector meson direct dilepton final states", "[gra::MPhotoVM][amplitude]") {
  ModelParamRestoreGuard restore;
  gra::MODELPARAM = "TUNE0";

  auto require_finite = [](const std::string &channel, double mass) {
    gra::LORENTZSCALAR lts = MakeToyPhotoZFFbar(13, mass);
    gra::MPhotoVM photovm(lts, gra::MModelTune::Load(modelfile), channel, gra::MPhotoVM::ProcessDefinitionFor(channel));
    const gra::amplitude::ProcessDefinition &process_definition = photovm;
    REQUIRE(process_definition.DecayStructureFor(lts) == gra::MPhotoVM::DirectDecayStructure());
    const double amp2 = photovm.Amp2(lts);
    REQUIRE_FALSE(lts.proton_good_walker.has_value());
    REQUIRE(lts.exact_forward_photon_kinematics);
    REQUIRE(std::isfinite(amp2));
    REQUIRE(amp2 >= 0.0);
    REQUIRE(amp2 > 0.0);
    REQUIRE(lts.hamp.size() == 16);
    for (std::size_t spin = 1; spin < 4; ++spin) {
      for (std::size_t i = 0; i < 4; ++i) { REQUIRE(lts.hamp[spin * 4 + i] == lts.hamp[i]); }
    }
    REQUIRE(amp2 == Approx(gra::SquaredNorm(std::vector<std::complex<double>>(lts.hamp.begin(), lts.hamp.begin() + 4)))
                        .epsilon(1e-12));
    for (const auto &amp : lts.hamp) {
      REQUIRE(std::isfinite(amp.real()));
      REQUIRE(std::isfinite(amp.imag()));
    }
    photovm.SampleColorFlow(lts);
    REQUIRE(lts.decaytree[0].p.color_flow.empty());
    REQUIRE(lts.decaytree[1].p.color_flow.empty());
  };

  require_finite("jpsi", 3.096900);
  require_finite("psi(2S)", 3.68610);
  require_finite("Upsilon(1S)", 9.46030);
  require_finite("Upsilon(2S)", 10.02326);
  require_finite("Upsilon(3S)", 10.3552);
}

TEST_CASE("MPhotoVM supports photon-emitter and UGD-target dissociation", "[gra::MPhotoVM][dissociation]") {
  ModelParamRestoreGuard restore;
  gra::MODELPARAM = "TUNE0";

  auto evaluate = [](bool excite1, bool excite2) {
    gra::LORENTZSCALAR lts = MakeToyPhotoZFFbar(13, 3.096900);
    if (excite1) { SetToyPhotoForwardExcitation(lts, 1, 2.0); }
    if (excite2) { SetToyPhotoForwardExcitation(lts, 2, 2.0); }
    gra::MPhotoVM photovm(lts, gra::MModelTune::Load(modelfile), "jpsi", gra::MPhotoVM::ProcessDefinitionFor("jpsi"));
    const double  amp2 = photovm.Amp2(lts);
    REQUIRE_FALSE(lts.proton_good_walker.has_value());
    const std::size_t direction_count = (excite1 || excite2) ? 2 : 1;
    const std::size_t spin_block      = 4 * direction_count;
    REQUIRE(lts.hamp.size() == 4 * spin_block);
    for (std::size_t spin = 1; spin < 4; ++spin) {
      for (std::size_t i = 0; i < spin_block; ++i) { REQUIRE(lts.hamp[spin * spin_block + i] == lts.hamp[i]); }
    }
    double direct_norm = 0.0;
    for (const auto &amplitude : lts.hamp) { direct_norm += std::norm(amplitude); }
    REQUIRE(amp2 == Approx(0.25 * direct_norm).epsilon(1e-12));
    REQUIRE(std::isfinite(amp2));
    REQUIRE(amp2 > 0.0);
    return amp2;
  };

  const double upper_excited = evaluate(true, false);
  const double lower_excited = evaluate(false, true);
  REQUIRE(upper_excited == Approx(lower_excited).epsilon(1e-10));
  REQUIRE(evaluate(true, true) > 0.0);

  // Require the generated target mass in both directional thresholds
  auto require_closed_direction = [](int target_leg) {
    gra::LORENTZSCALAR lts = MakeToyPhotoZFFbar(13, 3.096900);
    SetToyPhotoForwardExcitation(lts, target_leg, 120.0);
    const int         photon_leg = target_leg == 1 ? 2 : 1;
    const gra::M4Vec &q_photon   = photon_leg == 1 ? lts.q1 : lts.q2;
    const gra::M4Vec &p_target   = target_leg == 1 ? lts.pbeam1 : lts.pbeam2;
    const double      W          = std::sqrt((q_photon + p_target).M2());
    REQUIRE(W < lts.PDG.FindByPDG(443).mass + lts.pfinal[static_cast<std::size_t>(target_leg)].M());
    if (photon_leg == 1) {
      lts.x2  = 0.0;
      lts.xi2 = 0.0;
    } else {
      lts.x1  = 0.0;
      lts.xi1 = 0.0;
    }
    gra::MPhotoVM photovm(lts, gra::MModelTune::Load(modelfile), "jpsi", gra::MPhotoVM::ProcessDefinitionFor("jpsi"));
    REQUIRE(photovm.Amp2(lts) == Approx(0.0));
  };
  require_closed_direction(1);
  require_closed_direction(2);
}

TEST_CASE("MSubProc registers ygg photoproduction processes", "[gra::MPhotoVM][gra::MSubProc]") {
  ModelParamRestoreGuard restore;
  gra::MODELPARAM = "TUNE0";

  const gra::MSubProc proc_f({"ygg"}, "F");
  const gra::MSubProc proc_c({"ygg"}, "C");

  REQUIRE(proc_f.ProcessExist("ygg[Z]<F>"));
  REQUIRE(proc_f.ProcessExist("ygg[jpsi]<F>"));
  REQUIRE(proc_f.ProcessExist("ygg[psi(2S)]<F>"));
  REQUIRE(proc_f.ProcessExist("ygg[Upsilon(1S)]<F>"));
  REQUIRE(proc_f.ProcessExist("ygg[Upsilon(2S)]<F>"));
  REQUIRE(proc_f.ProcessExist("ygg[Upsilon(3S)]<F>"));
  REQUIRE(proc_c.ProcessExist("ygg[Z]<C>"));
  REQUIRE(proc_c.ProcessExist("ygg[jpsi]<C>"));
  REQUIRE(proc_c.ProcessExist("ygg[psi(2S)]<C>"));
  REQUIRE(proc_c.ProcessExist("ygg[Upsilon(1S)]<C>"));
  REQUIRE(proc_c.ProcessExist("ygg[Upsilon(2S)]<C>"));
  REQUIRE(proc_c.ProcessExist("ygg[Upsilon(3S)]<C>"));

  ToyHelicityProcess master;
  master.state.lts = MakeToyPhotoZFFbar(13, 3.096900);
  master.ProcPtr   = gra::MSubProc({"ygg"}, "F");
  master.ProcPtr.Initialize("ygg", "jpsi");
  master.InitializeProcessAmplitude();
  gra::MSubProc      active = master.ProcPtr;
  gra::LORENTZSCALAR lts    = master.state.lts;
  const double       amp2   = active.GetBareAmplitude2(lts);
  REQUIRE(std::isfinite(amp2));
  REQUIRE(amp2 > 0.0);
  active.SampleColorFlow(lts);
  REQUIRE(lts.decaytree[0].p.color_flow.empty());
  REQUIRE(lts.decaytree[1].p.color_flow.empty());
}

// Match the continued forward vertex to physical spin sums and to the nuclear VMD profile
TEST_CASE("GP sliding photo reference matches production and target profile",
          "[gra::MRegge][GP][photoproduction][normalization]") {
  const auto tune = gra::MModelTune::Load(modelfile);
  for (const bool reverse : {false, true}) {
    for (const bool derivative : {false, true}) {
      CAPTURE(reverse, derivative);
      auto lts = MakeToyCoherentPhotonLTS();
      lts.process.MMAX = 2;
      lts.process.FORWARD_NOFLIP = true;
      lts.process.FORWARD_VERTEX = gra::ForwardVertexMode::HelicityResidue;
      lts.process.SPINGEN = true;
      lts.process.SPINDEC = false;
      lts.process.DERIVATIVE_FACTOR = derivative;
      gra::MRegge regge(lts, tune, gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "GP photo profile"));
      const auto param = ReggeParametersForTest(regge, lts);
      auto res = MakeToyPhotoGPResonance(lts.process.MMAX);
      res.production_model = gra::ReggeProductionModel::GP;
      auto &production = res.production.front();
      if (reverse) {
        std::swap(production.tree[0], production.tree[1]);
      }
      production.hel.alpha_ls.Set(2, 2, {0.5, 0.25});
      gra::gpom::InitResonanceLS(production.hel, 1, reverse ? 4 : 2, reverse ? 2 : 4);
      const auto exchange = param->exchanges.at(gra::regge::TrajectoryIndex(*param, 990)).soft_exchange;

      // With vanishing target transverse transfer, only its m=0 column survives
      auto forward = lts;
      const int target = reverse ? 1 : 2;
      forward.pfinal[target].SetPxPyPzM(0.0, 0.0, forward.pfinal[target].Pz(), gra::PDG::mp);
      UpdateToyDerivedKinematics(forward);
      gra::gpom::AmpCache cache;
      const auto matrices = gra::gpom::Resonance(forward, *param, res, &cache);
      const auto &source = cache.sources.front();
      const auto &photon = reverse ? source.lower_residue : source.upper_residue;
      const double flux = gra::spin::SourceSpinAveragedDensity(photon, 2, "GP forward spin sum");
      const double t = reverse ? forward.t1 : forward.t2;
      const double momentum = gra::kinematics::PairRest(forward.q1_in_X, forward.q2_in_X).momentum;
      const auto coupling =
          gra::gpom::PhotoCoupling(production, param->soft_model->Alpha(exchange, t), momentum, derivative);
      CHECK(0.25 * matrices.front().FrobNorm2() / flux == Approx(std::norm(coupling)).epsilon(1.0e-10));
      const int mmax = lts.process.MMAX;
      const auto residue = source.upper_residue[0][mmax - (reverse ? 0 : 1)] *
                           source.lower_residue[0][mmax - (reverse ? 1 : 0)];
      REQUIRE(std::abs(residue) > 0.0);
      CHECK(std::abs(matrices.front()[0][reverse ? 0 : 2] / residue - coupling) < 1.0e-12);

      lts.process.RESONANCES = {{"rho", res}};
      const std::array<gra::ForwardLegState, 2> state = {gra::ResolveForwardLegState(lts, gra::ForwardBeamLeg::Upper),
                                                         gra::ResolveForwardLegState(lts, gra::ForwardBeamLeg::Lower)};
      const std::array<double, 2> subenergy = {(lts.q2 + lts.pbeam1).M2(), (lts.q1 + lts.pbeam2).M2()};
      const auto terms = regge.PhotoAmplitudes(lts, gra::ReggeProductionModel::GP, state, subenergy, {true, true});
      REQUIRE(terms.size() == 1);
      const auto &term = terms.front();
      const auto amplitude = [&](double transfer) {
        const double mass2 = lts.pfinal[0].M2();
        const double photon_t = reverse ? lts.t2 : lts.t1;
        const double q = gra::kinematics::SqrtKallenLambda(mass2, photon_t, transfer) / (2.0 * std::sqrt(mass2));
        return gra::gpom::PhotoCoupling(production, param->soft_model->Alpha(exchange, transfer), q, derivative) *
               param->soft_model->PhysicalResidue(exchange, transfer) * regge.PhotoKernel(term.w2, transfer, 113);
      };
      const auto expected = amplitude(0.0);
      CHECK(std::abs(term.target_forward - expected) < 1.0e-12 * std::abs(expected));
      CHECK(term.target_eta == Approx(expected.real() / expected.imag()).epsilon(1.0e-10));
      // An independent smaller difference checks the local slope including the continued vertex
      const double dt = 1.0e-6;
      CHECK(term.target_slope == Approx(std::log(std::norm(expected) / std::norm(amplitude(-dt))) / dt).epsilon(0.001));
      REQUIRE(gra::SquaredNorm(term.amplitude) > 0.0);

      // A pure target helicity flip keeps finite-transfer events despite zero VMD absorption
      gra::RES_PRODUCTION_CHANNEL flip;
      flip.exchange = {22, 990};
      flip.basis = gra::ReggeVertexBasis::Helicity;
      flip.helicity = {{-1, -1}};
      flip.g_helicity = {{1.0, 0.2}};
      const auto pdg = LoadedPDGTable();
      auto &zero = lts.process.RESONANCES.at("rho").production.front();
      zero.hel = gra::gpom::PrepareResonance(res.p, {zero.tree[0].p, zero.tree[1].p}, flip, *param, pdg, lts.process.MMAX, 0.0);
      const auto flipped = regge.PhotoAmplitudes(lts, gra::ReggeProductionModel::GP, state, subenergy, {true, true});
      REQUIRE(flipped.size() == 1);
      CHECK(std::abs(flipped.front().target_forward) < 1.0e-15);
      REQUIRE(gra::SquaredNorm(flipped.front().amplitude) > 0.0);
    }
  }
}

// Integrate the actual single-photon source in the exact lepton phase-space measure
// [REFERENCE: Frixione et al., arXiv:hep-ph/9310350, Eqs. (17)-(20)]
TEST_CASE("Single-photon EPA recoil integrates to the Weizsacker-Williams spectrum", "[photoproduction][EPA][normalization]") {
  auto lts = MakeToyPhotoZFFbar(13);
  gra::RequireModelCache(lts.model_cache, gra::MModelTune::Load(modelfile), "EPA recoil test");
  lts.process.PHOTON_VERTEX = "EPA";
  for (const auto leg : {gra::ForwardBeamLeg::Upper, gra::ForwardBeamLeg::Lower}) {
    auto state = gra::ResolveForwardLegState(lts, leg);
    state.emitter = lts.PDG.FindByPDG(11);
    if (leg == gra::ForwardBeamLeg::Upper) { lts.beam1 = state.emitter; }
    else { lts.beam2 = state.emitter; }
    for (const double x : {0.03, 0.2, 0.7}) {
      state.xi = x;
      state.has_xi = true;
      const double qmin = pow2(state.emitter.mass * x) / (1-x), qmax = 0.8;
      const auto [nodes, weights] = gra::math::GaussLegendreRule(64, std::log(qmin), std::log(qmax));
      double integral = 0.0;
      for (const auto &i : indices(nodes)) {
        const double q2 = std::exp(nodes[i]);
        state.t = -q2;
        state.qt = std::sqrt((1-x) * (q2-qmin));
        const auto scalar = gra::flux::PhotoSourceScalarAmplitudes(lts,state).front();
        double transverse = 0.0;
        for (const int h : {-1,1}) {
          const auto amplitude = gra::flux::PhotoSourceAmplitudes(lts,state,h).front();
          const auto expected = std::sqrt(1-x) * gra::qed::TransversePhotonSourceAmplitude(lts,state,h);
          RequireComplexNear(amplitude,expected,1e-12);
          transverse += std::norm(amplitude);
        }
        CHECK(transverse == Approx(std::norm(scalar)).epsilon(1e-12));
        integral += weights[i] * q2 * std::norm(scalar) * x / (16*gra::math::PIPI);
      }
      const double analytic = gra::qed::alpha_QED() / (2*gra::math::PI) *
          ((1+pow2(1-x))/x * std::log(qmax/qmin) - 2*(1-x)/x * (1-qmin/qmax));
      CHECK(integral == Approx(analytic).epsilon(2e-12));
    }
  }
}
