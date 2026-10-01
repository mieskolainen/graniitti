// Gamma-Gamma Amplitudes
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <bit>
#include <cmath>
#include <complex>
#include <cstdint>
#include <random>
#include <limits>
#include <utility>
#include <vector>

// Own
#include "Graniitti/Amplitude/MG5/Runtime/AMP_MG5_Helicity.h"
#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/Math/MCombinatorics.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Process/MDecayTree.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Particle/MForm.h"
#include "Graniitti/Particle/MPDG.h"
#include "Graniitti/Particle/MResonance.h"
#include "Graniitti/Photon/MGamma.h"
#include "Graniitti/Photon/MQED.h"
#include "Graniitti/Regge/MRegge.h"
#include "Graniitti/Spin/MHelicity.h"
#include "Graniitti/Spin/MHelicityDecay.h"
#include "Graniitti/Tech/MException.h"

using gra::aux::indices;
using gra::math::abs2;
using gra::math::msqrt;
using gra::math::PI;
using gra::math::pow2;
using gra::math::zi;

using namespace gra::form;

namespace gra {
namespace {

// Compute the complex CP-even scalar diphoton amplitude from its pole coefficient
// H_-- = H_++ = sqrt(32 pi M^2 Gamma_gg Gamma)/(s-M^2+i M Gamma)
bool ScalarYYResonanceAmplitude(const double shat, const double mass, const double gamma_yy, const double gamma_tot,
                                std::complex<double> &amplitude) {
  amplitude = 0.0;
  if (!std::isfinite(shat) || !std::isfinite(mass) || !(mass > 0.0) || !std::isfinite(gamma_yy) || gamma_yy < 0.0 ||
      !std::isfinite(gamma_tot) || !(gamma_tot > 0.0)) {
    return false;
  }
  const double residue2 = 32.0 * math::PI * pow2(mass) * gamma_yy * gamma_tot;
  if (!std::isfinite(residue2) || residue2 < 0.0) { return false; }
  amplitude = std::sqrt(residue2) * resonance::FixedWidthLineShape(shat, mass, gamma_tot);
  return std::isfinite(amplitude.real()) && std::isfinite(amplitude.imag());
}

// Collect doubled external helicities in the same leaf order as the Jacob-Wick matrix
void ScalarLeafHelicities(const std::vector<MDecayBranch> &tree, std::vector<std::vector<int>> &helicities) {
  for (const auto &branch : tree) {
    if (!branch.legs.empty()) {
      ScalarLeafHelicities(branch.legs, helicities);
      continue;
    }
    std::vector<int> labels;
    for (const double helicity : spin::FinalStateHelicities(branch.p, "MGamma scalar decay")) {
      labels.push_back(static_cast<int>(std::lround(2.0 * helicity)));
    }
    helicities.push_back(std::move(labels));
  }
}

// Dress inclusive scalar production with a normalized physical Jacob-Wick decay
// Integral dPhi |g_decay D|^2 / S = 2 M Gamma_partial at the pole
// [REFERENCE: PDG 2025, Kinematics review, Eq. (49.11)]
std::vector<mg5helas::HelicityComponent> ScalarDecayTensor(LORENTZSCALAR &lts, const std::complex<double> production,
                                                         const double mass, const double width) {
  MMatrix<std::complex<double>> decay = {{1.0}};
  std::vector<std::vector<int>> final = {{}};
  if (lts.process.root_decay_mode == RootDecayMode::Physical) {

    const auto &res = lts.process.ROOT_RES;
    decay = spin::ResonanceDecayMatrix(lts, res, "CM") * (res.hel_decay.g_decay / std::sqrt(2.0 * mass * width));
    std::vector<std::vector<int>> helicities;
    ScalarLeafHelicities(lts.decaytree, helicities);
    final.clear();
    math::IndexComb(helicities, 0, std::vector<int>{}, final);

  }

  std::vector<mg5helas::HelicityComponent> components;
  components.reserve(4 * final.size());
  constexpr auto labels = spin::BinaryHelicityLabelsX2();
  for (const auto &i : indices(final)) {
    for (const int upper : labels) {
      for (const int lower : labels) {
        components.push_back({{upper, lower}, final[i], 0, upper == lower ? production * decay[0][i] : 0.0});
      }
    }
  }
  return components;
}

// Compute true for a fermion species implemented by the gamma-gamma pair
// amplitude
bool IsGammaFermionPairParticle(const int apdg) {
  return (apdg >= 1 && apdg <= 6) || apdg == 11 || apdg == 13 || apdg == 15 || apdg == std::abs(PDG::PDG_monopole);
}

// Compute true for one direct same-flavour particle-antiparticle pair
bool GammaDirectFermionPair(const std::vector<MDecayBranch> &tree) {
  if (tree.size() != 2 || !tree[0].legs.empty() || !tree[1].legs.empty()) { return false; }
  const int first_pdg = tree[0].p.pdg;
  return first_pdg == -tree[1].p.pdg && tree[0].p.spinX2 == 1 && tree[1].p.spinX2 == 1 &&
         IsGammaFermionPairParticle(std::abs(first_pdg));
}

// Map one physical direct pair to the generated electron process for one call
class PairPDGMap {
 public:
  explicit PairPDGMap(std::vector<MDecayBranch> &tree) : first(tree[0].p.pdg), second(tree[1].p.pdg), branches(tree) {
    branches[0].p.pdg = -11;
    branches[1].p.pdg = 11;
  }
  ~PairPDGMap() {
    branches[0].p.pdg = first;
    branches[1].p.pdg = second;
  }
  PairPDGMap(const PairPDGMap &)            = delete;
  PairPDGMap &operator=(const PairPDGMap &) = delete;

 private:
  int                        first;
  int                        second;
  std::vector<MDecayBranch> &branches;
};

// Compute true for one physical quark-pair alias mode in kernel order
bool GammaDirectPartonMode(const std::vector<MDecayBranch> &tree) {
  return GammaDirectFermionPair(tree) && tree[0].p.pdg > 0 && tree[0].p.pdg <= 6;
}

// Compute true when one decay tree contains the requested particle
bool GammaDecayTreeContains(const std::vector<MDecayBranch> &tree, const int pdg) {
  for (const auto &branch : tree) {
    if (std::abs(branch.p.pdg) == std::abs(pdg) || GammaDecayTreeContains(branch.legs, pdg)) { return true; }
  }
  return false;
}

// Compute whether one gamma amplitude mode needs monopole parameters
bool GammaModeNeedsMonopoleParameters(const MGammaMode mode, const LORENTZSCALAR &lts) {
  if (mode == MGammaMode::Monopolium) { return true; }
  return mode == MGammaMode::FermionPair && GammaDecayTreeContains(lts.decaytree, PDG::PDG_monopole);
}

// Load immutable monopole parameters for one process state
MMonopoleParamPtr GammaMonopoleParameters(const LORENTZSCALAR &lts, const MModelTunePtr &model_tune) {
  if (model_tune == nullptr) { throw std::invalid_argument("MGamma::InitializeParameters: missing model tune"); }
  if (lts.model_cache == nullptr) {
    throw std::logic_error("MGamma::InitializeParameters: missing run owned model cache");
  }
  lts.model_cache->Bind(model_tune);
  return GetMonopoleParam(*lts.model_cache, lts.PDG.FindByPDG(PDG::PDG_monopole).mass,
                          lts.PDG.FindByPDG(PDG::PDG_monopolium).mass, lts.PDG.FindByPDG(PDG::PDG_monopolium).width);
}

// Compute the fixed process identity for one gamma-gamma process mode
std::string GammaProcessName(const MGammaMode mode) {
  switch (mode) {
    case MGammaMode::Generic:
      return "yy_generic";
    case MGammaMode::Higgs:
      return "yy_Higgs";
    case MGammaMode::Monopolium:
      return "yy_monopolium0";
    case MGammaMode::FermionPair:
      return "yy_ffbar";
    case MGammaMode::Flux:
      return "yy_flux";
  }
  throw std::invalid_argument("MGamma::ProcessDefinitionFor: unknown process mode");
}

// Compute the final-state pattern for one gamma-gamma process mode
std::string GammaProcessPattern(const MGammaMode mode) {
  switch (mode) {
    case MGammaMode::Generic:
      return "event gamma-gamma final state";
    case MGammaMode::Higgs:
      return "Higgs decay tree";
    case MGammaMode::Monopolium:
      return "monopolium decay tree";
    case MGammaMode::FermionPair:
      return "l+l-, qqbar, monopole pairs";
    case MGammaMode::Flux:
      return "arbitrary central state";
  }
  throw std::invalid_argument("MGamma::ProcessDefinitionFor: unknown process mode");
}

}  // namespace

// Build one immutable gamma-gamma process definition for an analytic mode
std::shared_ptr<const amplitude::ProcessDefinition> MGamma::ProcessDefinitionFor(const MGammaMode mode) {
  return std::make_shared<amplitude::AnalyticProcess>(
      "MGAMMA", GammaProcessName(mode), GammaProcessPattern(mode),
      DecayStructure{mode == MGammaMode::FermionPair ? DecayType::Full : DecayType::None},
      [mode](const std::vector<MDecayBranch> &tree) {
        if (mode == MGammaMode::FermionPair) { return GammaDirectFermionPair(tree); }
        if (mode == MGammaMode::Flux) { return true; }
        return !tree.empty();
      },
      [mode](const LORENTZSCALAR &lts) {
        const bool scalar_decay = (mode == MGammaMode::Higgs || mode == MGammaMode::Monopolium) &&
                                  lts.process.root_decay_mode == RootDecayMode::Physical;
        if (scalar_decay && lts.decaytree.size() != 2) {
          throw std::invalid_argument("MGamma: physical scalar decays require a binary Jacob-Wick root");
        }
        if (mode == MGammaMode::FermionPair) { return DecayStructure{DecayType::Full}; }
        return scalar_decay ? decay::JacobWickStructure(lts) : DecayStructure{};
      },
      mode == MGammaMode::FermionPair ? amplitude::TopologyCondition(GammaDirectPartonMode)
                                      : amplitude::TopologyCondition{});
}

// Initialize mode-dependent monopole parameters before worker process copies
void MGamma::InitializeParameters(const MProcessSetup &setup, const MGammaMode mode) {
  if (!GammaModeNeedsMonopoleParameters(mode, setup.lts)) { return; }
  const auto parameters = GammaMonopoleParameters(setup.lts, setup.model_tune);
  if (mode == MGammaMode::Monopolium) {
    parameters->ValidateMonopolium();
    if (parameters->monopolium_mass > setup.lts.sqrt_s) {
      throw std::invalid_argument("MGamma::InitializeParameters: monopolium mass exceeds the collision energy");
    }
  }
}

// Construct the gamma-gamma amplitude with one immutable process definition
MGamma::MGamma(gra::LORENTZSCALAR &lts, MModelTunePtr model_tune_snapshot,
               std::shared_ptr<const amplitude::ProcessDefinition> definition)
    : amplitude::ProcessFamily(std::move(definition)),
      model_tune(RequireModelCache(lts.model_cache, model_tune_snapshot, "MGamma").TunePtr()),
      higgs(lts.PDG.FindByPDG(25)),
      higgs_gamma_gamma_width(gra::resonance::GammaGammaPartialWidth(higgs, model_tune->Directory())) {
  std::map<int, double> species;
  if (GammaDirectFermionPair(lts.decaytree)) {
    const auto &first = lts.decaytree[0].p;
    const auto &second = lts.decaytree[1].p;
    if (std::abs(first.mass - second.mass) > 64.0 * std::numeric_limits<double>::epsilon() *
                                               std::max({1.0, first.mass, second.mass})) {
      throw std::invalid_argument("MGamma: a Dirac pair requires the same particle and antiparticle mass");
    }
    species.emplace(std::abs(first.pdg), first.mass);
  }
  for (const auto &mode : lts.final_state_parton_modes) {
    if (mode.size() != 2 || mode[0] != -mode[1]) { continue; }
    const int pdg = std::abs(mode[0]);
    const auto mass = lts.final_state_parton_masses.find(pdg);
    species.emplace(pdg, mass != lts.final_state_parton_masses.end() ? mass->second : lts.PDG.FindByPDG(pdg).mass);
  }
  if (species.empty()) { return; }
  const auto process = amplitude::FindProcess("PHOTON", "yy_ll");
  if (!process.has_value()) { throw std::invalid_argument("MGamma: missing massive Dirac kernel"); }
  for (const auto &[pdg, mass] : species) {
    SLHAReader card(aux::ResolveProjectPath(*amplitude::ParameterCard("PHOTON", "yy_ll")));
    // The analytic Dirac kernel uses one mass parameter at both photon vertices
    card.set_block_entry("mass", 11, mass);
    card.set_block_entry("yukawa", 11, mass);
    auto kernel = CreatePhotonMG5Process(*process);
    kernel->InitParameters(std::move(card));
    pair_amplitudes.emplace(pdg, std::move(kernel));
  }
}

// Compute the physical coupling multiplier of one direct fermion pair
bool MGamma::PairScale(const gra::LORENTZSCALAR &lts, double &scale) const {
  scale = 1.0;
  if (!GammaDirectFermionPair(lts.decaytree) || pair_amplitude == nullptr) { return false; }

  const int apdg = std::abs(lts.decaytree[0].p.pdg);
  if (apdg <= 6) {
    // Two photon vertices give Q^2 and three orthogonal colors give sqrt(Nc)
    const double charge = lts.decaytree[0].p.chargeX3 / 3.0;
    scale               = std::sqrt(3.0) * pow2(charge);
    return std::isfinite(scale);
  }
  if (apdg != std::abs(PDG::PDG_monopole)) { return true; }

  const MMonopoleParam &param     = MonopoleParameters(lts);
  const double          threshold = 4.0 * pow2(param.monopole_mass);
  if (!std::isfinite(lts.s_hat) || !(lts.s_hat >= threshold)) { return false; }
  double beta = 1.0;
  if (param.coupling == "beta-dirac") {
    const double beta2 = 1.0 - threshold / lts.s_hat;
    if (!std::isfinite(beta2) || beta2 < 0.0) { return false; }
    beta = std::sqrt(beta2);
  } else if (param.coupling != "dirac") {
    throw std::invalid_argument("MGamma::PairScale: Unknown PARAM_MONOPOLE.coupling " + param.coupling);
  }

  // The generated kernel contains two electric vertices, so replace e^2 by
  // (g beta)^2 at complex amplitude level
  const double alpha_ref = mg5helas::AlphaQEDAtZero(lts) ? qed::alpha_QED() : pair_amplitude->DefaultAlphaQED();
  const double electric2 = 4.0 * math::PI * alpha_ref;
  const double g         = 2.0 * math::PI * param.gn / qed::e_QED();
  if (!std::isfinite(electric2) || !(electric2 > 0.0) || !std::isfinite(g)) { return false; }
  scale = pow2(g * beta) / electric2;
  return std::isfinite(scale);
}

// Build and contract one complex CP-even scalar diphoton hard tensor
double MGamma::ScalarYY(gra::LORENTZSCALAR &lts, const double mass, const double gamma_yy, const double gamma_tot) {
  lts.epa_hard.Clear();
  lts.hard_color_flows.clear();
  lts.hamp.clear();
  std::complex<double> same_helicity;
  if (!ScalarYYResonanceAmplitude(lts.s_hat, mass, gamma_yy, gamma_tot, same_helicity)) {
    lts.hamp.clear();
    return 0.0;
  }

  if (lts.pfinal.empty()) {
    lts.hamp.clear();
    return 0.0;
  }
  auto components = ScalarDecayTensor(lts, same_helicity, mass, gamma_tot);

  // The central state stays fixed while only the forward transfers move in
  // the screening loop
  std::vector<M4Vec> final = {lts.pfinal[0]};
  if (mg5helas::HasEPAHardBeamState(lts)) {
    mg5helas::EPAHardFrame frame;
    if (!mg5helas::PrepareEPAHardFrame(lts, final, frame)) {
      lts.hamp.clear();
      return 0.0;
    }
    auto amplitudes = mg5helas::ContractEPAHardPhotonSources(lts, components, frame);
    if (amplitudes.empty() || !gra::AllFinite(amplitudes)) {
      lts.hamp.clear();
      return 0.0;
    }
    for (auto &component : components) { component.value *= 0.5; }
    lts.epa_hard.frame         = frame;
    lts.epa_hard.amplitude     = std::move(components);
    lts.epa_hard.normalization = 2.0;
    lts.hamp.assign(amplitudes.begin(), amplitudes.end());
    return gra::SquaredNorm(amplitudes) / 4.0;
  }

  M4Vec k1;
  M4Vec k2;
  if (!mg5helas::PrepareOnShellKinematics(lts, final, k1, k2)) {
    lts.hamp.clear();
    return 0.0;
  }
  lts.hamp = mg5helas::ContractEPAPhotonSources(lts, components, k1, k2);
  if (lts.hamp.empty() || !gra::AllFinite(lts.hamp)) {
    lts.hamp.clear();
    return 0.0;
  }
  return gra::SquaredNorm(lts.hamp) / 4.0;
}

// Compute immutable parameters required by a monopole amplitude
const MMonopoleParam &MGamma::MonopoleParameters(const gra::LORENTZSCALAR &lts) const {
  std::call_once(monopole_parameter_load_once,
                 [this, &lts] { monopole_param_handle = GammaMonopoleParameters(lts, model_tune); });
  return *monopole_param_handle;
}

// ============================================================================
// (yy -> fermion-antifermion pair)
// yy -> e+e-, mu+mu-, tau+tau-, qqbar or Monopole-Antimonopole (spin-1/2
// monopole)
// production
//
//
// This is the same amplitude as yy -> e+e- (two diagrams)
// obtained by crossing e+e- -> yy annihilation, see:
//
// http://theory.sinp.msu.ru/comphep_old/tutorial/QED/node4.html
//
// beta = v/c of the monopole (or antimonopole) in their CM system
//
//
// *************************************************************
// Dirac quantization condition:
//
// g = 2\pi \hbar n / (\mu_0 e),    where n = 1,2,3,...
// where alpha = e^2/(4*PI) is the running QED coupling
//
//
// Numerical values:
// alpha_g  = g^2 / (4*PI) ~ 34  (when n = 1)
// alpha_em = e^2 / (4*PI) ~ 1/137
// *************************************************************
//
// From lepton pair to monopole pair:
// replace e -> g*beta
//
// 16 * pow2(PI * alpha_EM) -> pow(g*beta, 4)
//
// *************************************************************
// Coupling schemes:
//
// Dirac:     alpha_g = g^2/(4pi)
// Beta-dirac alpha_g = (g*beta)^2 / (4pi)
//
// [REFERENCE: Rajantie, physicstoday.scitation.org/doi/pdf/10.1063/PT.3.3328]
// [REFERENCE: Dougall, Wick, arxiv.org/abs/0706.1042]
// [REFERENCE: Rels, Sauter, arxiv.org/abs/1707.04170v1]
//
mg5helas::MatrixElementEvaluation MGamma::yyffbar(gra::LORENTZSCALAR &lts, bool coherent_epa) {
  lts.epa_hard.Clear();
  lts.hard_color_flows.clear();
  if (!GammaDirectFermionPair(lts.decaytree)) {
    lts.hamp.clear();
    return {mg5helas::EvaluationStatus::AmplitudeFailure, 0.0};
  }

  // Select a model initialized before sampling with this species' fixed mass
  const auto kernel = pair_amplitudes.find(std::abs(lts.decaytree[0].p.pdg));
  pair_amplitude = kernel != pair_amplitudes.end() ? kernel->second.get() : nullptr;
  if (pair_amplitude == nullptr) {
    lts.hamp.clear();
    return {mg5helas::EvaluationStatus::AmplitudeFailure, 0.0};
  }
  mg5helas::MatrixElementEvaluation result;
  {
    PairPDGMap mapped(lts.decaytree);
    result = pair_amplitude->Evaluate(lts, 0.0, coherent_epa);
  }
  if (!result.Valid()) { return {result.status, 0.0}; }

  double scale = 1.0;
  if (!PairScale(lts, scale)) {
    lts.hamp.clear();
    return {mg5helas::EvaluationStatus::KinematicsFailure, 0.0};
  }
  const double factor = pow2(scale);
  if (!std::isfinite(factor) || !std::isfinite(result.amp2 * factor)) {
    lts.hamp.clear();
    return {mg5helas::EvaluationStatus::AmplitudeFailure, 0.0};
  }
  gra::Scale(lts.hamp, scale);
  if (lts.epa_hard.Ready()) {
    for (auto &component : lts.epa_hard.amplitude) { component.value *= scale; }
    for (auto &flow : lts.epa_hard.color) {
      for (auto &component : flow.components) { component.value *= scale; }
    }
  }
  result.amp2 *= factor;
  return result;
}

// Assign the unique color-singlet flow of a direct quark pair
bool MGamma::SampleColorFlow(gra::LORENTZSCALAR &lts) {
  if (!GammaDirectFermionPair(lts.decaytree)) { return true; }
  return std::abs(lts.decaytree[0].p.pdg) > 6 || AssignPhotonSingletQuarkPairColorFlow(lts);
}

// --------------------------------------------------------------------------------------------
// A Monopolium (monopole-antimonopole) bound state process
//
// See e.g
//
// [REFERENCE: Preskill, http://www.theory.caltech.edu/~preskill/pubs/preskill-1984-monopoles.pdf]
// [REFERENCE: Epele, Franchiotti, Garcia, Canal, Vento, https://arxiv.org/abs/hep-ph/0701133v2]
// [REFERENCE: Barrie, Sugamoto, Yamashita, https://arxiv.org/abs/1607.03987v3]
// [REFERENCE: Fanchiotti, Canal, Vento, https://arxiv.org/pdf/1703.06649.pdf]
// [REFERENCE: Reis, Sauter, https://arxiv.org/abs/1707.04170v1]
//
// *************************************************************
// Dirac quantization condition:
//
// g = 2\pi \hbar / (\mu_0 e) n,    where n = 1,2,3,...
// where alpha = e^2/(4*PI) is the running QED coupling
//
//
// Numerical values:
// alpha_g  = g^2 / (4*PI) ~ 34  (when n = 1)
// alpha_em = e^2 / (4*PI) ~ 1/137
// *************************************************************
//
double MGamma::yyMP(gra::LORENTZSCALAR &lts) {
  lts.epa_hard.Clear();
  lts.hard_color_flows.clear();
  const MMonopoleParam &monopole_param = MonopoleParameters(lts);

  // Lazy printing
  PrintMonopoleParameters(monopole_param, lts.sqrt_s);

  // Physical monopolium pole parameters from PDG 881
  const double M       = monopole_param.monopolium_mass;
  const double Gamma_M = monopole_param.monopolium_width;

  // Two coupling scenarios:
  const double g    = 2.0 * math::PI * monopole_param.gn / qed::e_QED();
  double       beta = 0.0;

  if (monopole_param.coupling == "beta-dirac") {
    if (!std::isfinite(lts.s_hat) || !(lts.s_hat >= pow2(M))) {
      lts.hamp.clear();
      return 0.0;
    }
    const double beta2 = 1.0 - pow2(M) / lts.s_hat;
    if (!std::isfinite(beta2) || beta2 < 0.0) {
      lts.hamp.clear();
      return 0.0;
    }
    beta = std::sqrt(beta2);
  } else if (monopole_param.coupling == "dirac") {
    beta = 1.0;
  } else {
    throw std::invalid_argument("MGamma::yyMP: Unknown PARAM_MONOPOLE.coupling " + monopole_param.coupling);
  }

  // Magnetic coupling
  const double alpha_g = pow2(beta * g) / (4.0 * math::PI);

  // Energy-dependent diphoton partial width
  const double Gamma_E = monopole_param.GammaGamma(alpha_g);

  // [REFERENCE: Epele et al., arXiv:0809.0272, Eqs. (4)-(7)]
  return ScalarYY(lts, M, Gamma_E, Gamma_M);
}

// --------------------------------------------------------------------------------------------
// Gamma-Gamma to SM-Higgs 0++ helicity amplitudes
//
// Generic narrow width yy -> X cross section in terms of partial decay widths:
//
// \sigma(yy -> X) = 8\pi^2/M_X (2J+1) \Gamma(X -> yy) \delta(shat - M_X^2) (1 +
// h1h2)
//                 = (8 * \pi)  (2J+1) \Gamma(X -> yy) \Gamma_X (1 + h1h2) /
//                   ((shat - M_X^2)^2 + M_X^2\Gamma_X^2),
//
// The second line uses the narrow-width replacement
// \delta(shat-M_X^2) -> M_X Gamma_X / [\pi((shat-M_X^2)^2+M_X^2 Gamma_X^2)]
// With the exact gamma-gamma flux 2*shat, the fixed on-shell pole coefficient below
// retains an additional M_X^2/shat away from the pole; both forms agree in the
// narrow-width limit
//
// where h1,h2 = +- gamma helicities
//
// [REFERENCE: Khoze, Martin, Ryskin, https://arxiv.org/abs/hep-ph/0111078]
// [REFERENCE: Bernal, Lopez-Val, Sola, https://arxiv.org/pdf/0903.4978.pdf]
// [REFERENCE: Enterria, Lansberg, https://www.slac.stanford.edu/pubs/slacpubs/13750/slac-pub-13786.pdf]
//
double MGamma::yyHiggs(gra::LORENTZSCALAR &lts) {
  // The active PDG and decay tables are the single source for the pole and
  // partial widths
  return ScalarYY(lts, higgs.mass, higgs_gamma_gamma_width, higgs.width);
}

// Validate parameters shared by monopole and monopolium production
void MMonopoleParam::Validate() const {
  if (!std::isfinite(monopole_mass) || monopole_mass <= 0.0) {
    throw std::invalid_argument("MMonopoleParam: PDG 882 monopole mass should be finite and positive");
  }
  if (!std::isfinite(monopolium_mass) || monopolium_mass <= 0.0) {
    throw std::invalid_argument(
        "MMonopoleParam: PDG 881 monopolium mass "
        "should be finite and positive");
  }
  if (!std::isfinite(monopolium_width) || monopolium_width < 0.0) {
    throw std::invalid_argument(
        "MMonopoleParam: PDG 881 monopolium width "
        "should be finite and nonnegative");
  }
  if (gn < 1) { throw std::invalid_argument("MMonopoleParam: PARAM_MONOPOLE.gn should be at least one"); }
  if (coupling != "beta-dirac" && coupling != "dirac") {
    throw std::invalid_argument("MMonopoleParam: unknown PARAM_MONOPOLE.coupling " + coupling);
  }
  if (wavefunction != "COULOMB") {
    throw std::invalid_argument("MMonopoleParam: unknown PARAM_MONOPOLE.wavefunction " + wavefunction);
  }
}

// Validate the physical monopolium bound-state pole
void MMonopoleParam::ValidateMonopolium() const {
  Validate();
  if (monopolium_mass >= 2.0 * monopole_mass) {
    throw std::invalid_argument(
        "MMonopoleParam: PDG 881 mass should be below "
        "the two-monopole threshold");
  }
  if (monopolium_width <= 0.0) {
    throw std::invalid_argument(
        "MMonopoleParam: PDG 881 total width should be "
        "positive for yy -> monopolium");
  }
}

// Compute the monopolium binding energy relative to two free monopoles
// E_bind = M_bound-2m with the signed threshold convention
double MMonopoleParam::BindingEnergy() const { return monopolium_mass - 2.0 * monopole_mass; }

// Compute the Coulombic monopolium wavefunction at the origin
// psi(0) = (2-M/m)^(3/4) m^(3/2)/sqrt(pi)
double MMonopoleParam::PsiAtOrigin() const {
  const double binding_fraction = 2.0 - monopolium_mass / monopole_mass;

  // [REFERENCE: Epele et al., arXiv:0809.0272, Eq. (6)]
  return std::pow(binding_fraction, 3.0 / 4.0) * std::pow(monopole_mass, 3.0 / 2.0) / msqrt(PI);
}

// Compute the monopolium diphoton width for a magnetic coupling
// Gamma_gg = 32 pi alpha_g^2 |psi(0)|^2/M^2
double MMonopoleParam::GammaGamma(double alpha_g) const {


  // [REFERENCE: Epele et al., arXiv:0809.0272, Eq. (5)]
  return 32.0 * PI * pow2(alpha_g) / pow2(monopolium_mass) * math::abs2(PsiAtOrigin());
}

// Print monopole parameters once for this process helper instance
//
void MGamma::PrintMonopoleParameters(const MMonopoleParam &parameters, double) const {
  std::call_once(monopole_parameters_print_once, [&parameters] {
    aux::PrintBar("*");
    std::cout << rang::style::bold << "Monopolium process parameters:" << rang::style::reset << std::endl << std::endl;

    printf("- Monopole mass (PDG 882)        = %0.3f GeV \n", parameters.monopole_mass);
    printf("- Monopolium pole mass (PDG 881) = %0.3f GeV \n", parameters.monopolium_mass);
    printf("- Monopolium total width         = %0.3f GeV \n", parameters.monopolium_width);
    printf("- Binding energy                 = %0.3f GeV \n\n", parameters.BindingEnergy());

    std::cout << "Dirac charge = " << parameters.gn << std::endl;
    std::cout << "Coupling scheme = " << parameters.coupling << std::endl;
    std::cout << "Wavefunction = " << parameters.wavefunction << std::endl;
    aux::PrintBar("*");
    std::cout << std::endl;
  });
}

// Construct one immutable monopole parameter block
MMonopoleParamPtr ReadMonopoleParam(const MModelTune &tune, double monopole_mass, double monopolium_mass,
                                    double monopolium_width) {
  if (!std::isfinite(monopole_mass) || !std::isfinite(monopolium_mass) || !std::isfinite(monopolium_width)) {
    throw std::invalid_argument("ReadMonopoleParam: non-finite PDG mass or width");
  }

  auto param              = std::make_shared<MMonopoleParam>();
  param->monopole_mass    = monopole_mass;
  param->monopolium_mass  = monopolium_mass;
  param->monopolium_width = monopolium_width;

  try {
    const auto &block   = tune.General("PARAM_MONOPOLE");
    param->coupling     = block.at("coupling");
    param->wavefunction = block.at("wavefunction");
    param->gn           = block.at("gn");
  } catch (const nlohmann::json::exception &e) {
    throw std::invalid_argument("ReadMonopoleParam: Error parsing " + tune.GeneralFile() + ": " + e.what());
  }

  param->Validate();
  param->initialized = true;
  return param;
}

// Compute one run owned immutable monopole parameter block
MMonopoleParamPtr GetMonopoleParam(MModelCache &cache, double monopole_mass, double monopolium_mass,
                                   double monopolium_width) {
  const std::string key = "monopole:" + std::to_string(std::bit_cast<std::uint64_t>(monopole_mass)) + ":" +
                          std::to_string(std::bit_cast<std::uint64_t>(monopolium_mass)) + ":" +
                          std::to_string(std::bit_cast<std::uint64_t>(monopolium_width));
  return cache.Get<MMonopoleParam>(key, [&cache, monopole_mass, monopolium_mass, monopolium_width] {
    return ReadMonopoleParam(cache.Tune(), monopole_mass, monopolium_mass, monopolium_width);
  });
}

}  // namespace gra
