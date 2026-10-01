// GP production vertex preparation
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Regge/MReggeGPInit.h"
#include "Graniitti/Regge/MReggeGP.h"
#include "Graniitti/Regge/MReggeInit.h"
#include "Graniitti/Process/MHelicityConfig.h"
#include "Graniitti/Tech/MAux.h"
#include "rang.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <set>
#include <stdexcept>

using gra::aux::indices;

namespace gra::gpom {

constexpr double kComparisonTolerance = 1e-9;

// Initialize a GP analytic trajectory H(lambda1,lambda2;m) matrix
// Rows are physical leg helicities and columns are exchange m labels
void InitHelicity(HELMatrix& hc, double s1, double s2, int MMAX, const std::string& context) {
  const std::size_t n1 = spin::SpinStateCount(s1, context + " leg 1");
  const std::size_t n2 = spin::SpinStateCount(s2, context + " leg 2");
  (void)AnalyticMIndex(0, MMAX, context);
  const std::size_t nm = static_cast<std::size_t>(2 * MMAX + 1);

  hc.coupling_basis = CouplingBasis::Unspecified;
  hc.domain         = HelicityDomain::ReggeTrajectory;
  hc.exchange_basis = ExchangeBasisType::HelicityTransport;
  hc.analytic_MMAX  = MMAX;
  hc.T_active.clear();
  hc.gp_orbital = {};
  hc.m_ls.clear();
  hc.m_orbital.clear();
  hc.alpha_ls.Clear();
  spin::InitTwoBodyBasis(hc, -1.0, s1, s2, spin::SpinProjections(static_cast<double>(MMAX)), spin::SpinProjections(s1),
                         spin::SpinProjections(s2), context);
  hc.T     = HelAmp(n1 * n2, nm, 0.0);
  hc.T_set = MMatrix<bool>(n1 * n2, nm, false);
}

// Initialize one angle-free crossed Regge-helicity vertex
void InitCrossed(HELMatrix& hc, const double s1, const double s2, const int MMAX, const std::string& context) {
  InitHelicity(hc, s1, s2, MMAX, context);
  hc.exchange_basis = ExchangeBasisType::ReggeHelicity;
}

// Initialize one crossed photon vertex in its physical spin one basis
void InitPhotonCrossed(HELMatrix& hc, const double s1, const double s2, const std::string& context) {
  InitHelicity(hc, s1, s2, 1, context);
  // Prepare the physical transverse-photon rotation used after pole projection
  auto pole      = hc;
  pole.J         = 1.0;
  pole.Jz_values = {-1.0, 1.0};
  spin::InitJWRotation(pole);
  hc.jw_rotation = std::move(pole.jw_rotation);
}

// Validate a GP analytic trajectory H(lambda1,lambda2;m) matrix
// Direct Regge couplings fix the tensor scale without normalization
void CheckHelicity(HELMatrix& hc, const std::string& context, bool verbose) {
  (void)verbose;
  if (!hc.UsesReggeDomain()) {
    throw std::invalid_argument(context + ": HELMatrix is not an analytic trajectory basis");
  }
  if (hc.analytic_MMAX < 0) { throw std::invalid_argument(context + ": analytic trajectory MMAX is not initialized"); }
  (void)AnalyticMIndex(0, hc.analytic_MMAX, context);
  const std::size_t nm = static_cast<std::size_t>(2 * hc.analytic_MMAX + 1);
  if (hc.T.size_row() != hc.lambda_values.size_row() || hc.T.size_col() != nm) {
    throw std::invalid_argument(context + ": analytic trajectory T matrix has invalid dimensions");
  }
  if (hc.T_set.isEmpty()) {
    hc.T_set = hc.T.Transform([](const auto& value) { return math::abs2(value) > 1e-24; });
  }
  if (hc.T_set.size_row() != hc.T.size_row() || hc.T_set.size_col() != hc.T.size_col()) {
    throw std::invalid_argument(context + ": analytic trajectory T_set matrix has invalid dimensions");
  }

  if (!hc.T.IsFinite() || hc.T.MaskedSquaredNorm(hc.T_set) <= 1e-24) {
    throw std::invalid_argument(context +
                                ": analytic trajectory tensor is empty or "
                                "non-finite");
  }
  PruneHelicityCouplings(hc, 0.0, context);
}

// Format an integer or half-integer helicity projection for diagnostics
std::string FormatProjection(const double value) {
  const int spin2 = static_cast<int>(std::llround(2.0 * std::abs(value)));
  return std::string(value < 0.0 ? "-" : "") + gra::aux::Spin2XtoString(spin2);
}

// Format fixed-spin particles and analytic trajectory metadata for diagnostics
std::string FormatParticleSpin(const MParticle &p) { return p.name.find("a(t)") != std::string::npos ? "a(t)" : gra::aux::NullableSpin2XtoString(p.spinX2); }

// Print the GP proton-side analytic vertex table in the m1,m2 basis
void PrintContinuum(const LORENTZSCALAR &state, const std::vector<int> &exchange_channel, const std::vector<MDecayBranch> &tree) {
  if (exchange_channel.size() != 2 || tree.size() != 2) { throw std::invalid_argument("Amplitude initialization: GP print expects two exchange legs"); }
  const int MMAX = state.process.MMAX;
  if (MMAX < 0) { throw std::invalid_argument("Amplitude initialization: GP MMAX must be non-negative"); }
  std::cout << std::endl << "Amplitude initialization: GP forward helicity residues" << std::endl << std::endl << "[state 0 + state 1 <- beam legs]" << std::endl;
  gra::aux::PrintTable(
      {"ID", "Role", "PDG basis", "Type", "J", "P", "C", "PDG", "Name"},
      {{"0", "upper exchange", "trajectory metadata", "boson", FormatParticleSpin(tree[0].p), gra::aux::ParityToString(tree[0].p.P), gra::aux::ParityToString(tree[0].p.C), std::to_string(exchange_channel[0]), tree[0].p.name},
       {"1", "lower exchange", "trajectory metadata", "boson", FormatParticleSpin(tree[1].p), gra::aux::ParityToString(tree[1].p.P), gra::aux::ParityToString(tree[1].p.C), std::to_string(exchange_channel[1]), tree[1].p.name}});
  const bool barrier = state.process.FORWARD_VERTEX == ForwardVertexMode::HelicityResidue;
  std::cout << std::endl
            << "Analytic exchange labels: m1,m2 = [" << -MMAX << ", " << MMAX << "]" << std::endl
            << "Hadronic exchange residue: R_i(m_i) = N_i(m_i) " << (barrier ? "(qT_i/sqrt(s0))^|m_i| " : "") << "exp(+- i m_i phi_i), N_i = product_{k=0}^{|m_i|-1} [alpha_i(t_i)-k]" << std::endl
            << "Intact hadron SOFT residue: V_nf,ij = G_ij F_ij(t_i), V_f,ij = V_nf,ij kappa_ij exp(B_kappa,ij t_i/2) sqrt(-t_i)/(2 m_hel)" << std::endl
            << "Forward proton helicities: lambda_in,lambda_out = +-1/2" << (state.process.FORWARD_NOFLIP ? " with flip rows dropped" : " with flip rows kept") << std::endl
            << std::endl
            << "Photon legs use PARAM_REGGE.PHOTON_VERTEX = " << state.process.PHOTON_VERTEX << " (inclusive dissociation remains EPA)" << std::endl
            << std::endl;
  std::vector<std::vector<std::string>> rows;
  const auto                            upper_rows = spin::SpinHalfTransitions(state.process.FORWARD_NOFLIP);
  const auto                            lower_rows = spin::SpinHalfTransitions(state.process.FORWARD_NOFLIP);
  for (const auto &[h1, h3] : upper_rows) {
    for (const auto &[h2, h4] : lower_rows) {
      for (int m1 = -MMAX; m1 <= MMAX; ++m1) {
        for (int m2 = -MMAX; m2 <= MMAX; ++m2) {
          rows.push_back({FormatProjection(h1), FormatProjection(h3), FormatProjection(h2), FormatProjection(h4), std::to_string(m1), std::to_string(m2), FormatProjection(h1 - h3), FormatProjection(h2 - h4),
                          std::to_string(barrier ? std::abs(m1) : 0), std::to_string(barrier ? std::abs(m2) : 0)});
        }
      }
    }
  }
  gra::aux::PrintTable({"lambda1", "lambda3", "lambda2", "lambda4", "m1", "m2", "delta1", "delta2", "qT1 power", "qT2 power"}, rows);
  std::cout << std::endl;
}

// Compute the Jacob-Wick parity phase with analytic spin phases set by tau
int GPParityPhase(const MParticle &mother, const std::array<const MParticle *, 2> &legs, const regge::Param &param) {
  if (mother.P != -1 && mother.P != 1) { throw std::invalid_argument("Amplitude initialization: GP resonance has undefined parity"); }
  int phase = mother.P * IntegerPhaseSign(0.5 * static_cast<double>(mother.spinX2), "Amplitude initialization: GP parity");
  for (const MParticle *leg : legs) {
    if (leg->P != -1 && leg->P != 1) { throw std::invalid_argument("Amplitude initialization: GP exchange has undefined parity"); }
    phase *= leg->P;
    if (leg->pdg == PDG::PDG_gamma) {
      phase *= IntegerPhaseSign(0.5 * static_cast<double>(leg->spinX2), "Amplitude initialization: photon parity");
      continue;
    }
    const std::size_t      trajectory = regge::TrajectoryIndex(param, leg->pdg);
    const regge::Signature signature  = param.soft_model->Exchange(param.exchanges.at(trajectory).soft_exchange).Signature();
    phase *= regge::Tau(signature);
  }
  return phase;
}

// Distinguish helicity reversal under parity from exchange of identical legs
enum class SpinSymmetry { Parity, Exchange };

// Check H(partner) = phase H with the same relative tolerance for both symmetries
void ValidateRawSymmetry(const HELMatrix &hc, int phase, SpinSymmetry symmetry) {
  const int MMAX = hc.analytic_MMAX;
  const bool parity = symmetry == SpinSymmetry::Parity;
  const std::string context = parity ? "GP helicity Jacob-Wick parity" : "GP helicity identical-leg exchange";
  for (int m1 = -MMAX; m1 <= MMAX; ++m1) {
    for (int m2 = -MMAX; m2 <= MMAX; ++m2) {
      const auto i1 = AnalyticMIndex(m1, MMAX, context);
      const auto i2 = AnalyticMIndex(m2, MMAX, context);
      if (!hc.T_set[i1][i2]) { continue; }
      const auto p1 = AnalyticMIndex(parity ? -m1 : m2, MMAX, context);
      const auto p2 = AnalyticMIndex(parity ? -m2 : m1, MMAX, context);
      if (!hc.T_set[p1][p2]) { throw std::invalid_argument(context + ": symmetry partner is missing"); }
      const double scale = std::max({1.0, std::abs(hc.T[i1][i2]), std::abs(hc.T[p1][p2])});
      if (std::abs(hc.T[p1][p2] - static_cast<double>(phase) * hc.T[i1][i2]) > kComparisonTolerance * scale) {
        throw std::invalid_argument(context + ": inconsistent symmetry partner");
      }
    }
  }
}

// Check one analytic direct tensor against the physical pole LS subspace
void ValidateGPPoleSubspace(const HELMatrix &hc, const MParticle &mother, const MParticle &physical1, const MParticle &physical2, bool reverse, int exchange_phase) {
  const MParticle &pole1 = reverse ? physical2 : physical1;
  const MParticle &pole2 = reverse ? physical1 : physical2;
  if (pole1.spinX2 < 0 || pole2.spinX2 < 0 || pole1.spinX2 % 2 != 0 || pole2.spinX2 % 2 != 0) {
    throw std::invalid_argument(
        "Amplitude initialization: GP helicity requires integer pole "
        "spins");
  }
  HELMatrix pole;
  pole.C_symmetry      = hc.C_symmetry;
  pole.P_symmetry      = hc.P_symmetry;
  pole.coupling_basis  = gra::CouplingBasis::Helicity;
  const std::size_t n1 = static_cast<std::size_t>(pole1.spinX2 + 1);
  const std::size_t n2 = static_cast<std::size_t>(pole2.spinX2 + 1);
  pole.T               = MMatrix<std::complex<double>>(n1, n2, 0.0);
  pole.T_set           = MMatrix<bool>(n1, n2, false);
  const int j1         = pole1.spinX2 / 2;
  const int j2         = pole2.spinX2 / 2;
  for (int m1 = -j1; m1 <= j1; ++m1) {
    for (int m2 = -j2; m2 <= j2; ++m2) {
      if (std::abs(m1) > hc.analytic_MMAX || std::abs(m2) > hc.analytic_MMAX) { continue; }
      const std::size_t source1    = gpom::AnalyticMIndex(reverse ? m2 : m1, hc.analytic_MMAX, "GP pole subspace");
      const std::size_t source2    = gpom::AnalyticMIndex(reverse ? m1 : m2, hc.analytic_MMAX, "GP pole subspace");
      const std::size_t target1    = static_cast<std::size_t>(m1 + j1);
      const std::size_t target2    = static_cast<std::size_t>(m2 + j2);
      pole.T[target1][target2]     = static_cast<double>(reverse ? exchange_phase : 1) * hc.T[source1][source2];
      pole.T_set[target1][target2] = hc.T_set[source1][source2];
    }
  }
  const auto   expansion = gra::spin::DirectHelicityToLSCoefficients(pole, mother, pole1, pole2, true, gra::spin::VertexContext::Auto);
  const double norm2     = pole.T.FrobNorm2();
  if (expansion.residual_norm2 > 1.0e-10 * std::max(1.0, norm2)) {
    throw std::invalid_argument(
        "Amplitude initialization: GP helicity input is outside the "
        "canonical pole LS subspace");
  }
}

// Build one typed GP resonance steering matrix
HELMatrix PrepareResonance(const MParticle &mother, const std::vector<MParticle> &legs, const RES_PRODUCTION_CHANNEL &channel, const regge::Param &regge_param, const MPDG &pdg_table, int MMAX, double coupling_min) {
  if (legs.size() != 2 || MMAX < 0 || mother.spinX2 < 0 || mother.spinX2 % 2 != 0) {
    throw std::invalid_argument(
        "Amplitude initialization: GP resonance needs two analytic exchanges, "
        "non-negative MMAX and integer central spin");
  }
  (void)gpom::AnalyticMIndex(0, MMAX, "PrepareResonance");
  for (const auto &leg : legs) {
    if (leg.pdg != PDG::PDG_gamma && leg.spinX2 != gra::aux::kNullSpinX2) {
      throw std::invalid_argument(
          "Amplitude initialization: GP accepts only null-spin analytic PDG "
          "aliases or photons");
    }
  }
  if (channel.basis != ReggeVertexBasis::LS && channel.basis != ReggeVertexBasis::Helicity) {
    throw std::invalid_argument(
        "Amplitude initialization: GP basis should be g_ls or "
        "helicity");
  }

  HELMatrix hc;
  hc.BR            = 1.0;
  hc.P_symmetry    = channel.P_symmetry;
  hc.C_symmetry    = channel.C_symmetry;
  hc.J             = mother.spinX2 / 2.0;
  hc.domain        = gra::HelicityDomain::ReggeTrajectory;
  hc.analytic_MMAX = MMAX;
  hc.alpha_ls.Clear();
  const bool direct  = channel.exchange[0] == legs[0].pdg && channel.exchange[1] == legs[1].pdg;
  const bool reverse = !direct && channel.exchange[0] == legs[1].pdg && channel.exchange[1] == legs[0].pdg;
  if (!direct && !reverse) {
    throw std::invalid_argument(
        "Amplitude initialization: GP exchange ordering does not match "
        "the card channel");
  }
  const int j1x2 = gpom::PoleSpinX2(regge_param, pdg_table, legs[0].pdg);
  const int j2x2 = gpom::PoleSpinX2(regge_param, pdg_table, legs[1].pdg);
  if (channel.basis == ReggeVertexBasis::LS) {
    hc.coupling_basis  = gra::CouplingBasis::LS;
    hc.analytic_Lambda = channel.Lambda;
    hc.alpha_ls        = channel.g_ls;
    if (reverse) {
      // Transport the card coupling through the Jacob Wick daughter exchange
      for (spin::LSTerm &term : hc.alpha_ls) {
        const double exponent = static_cast<double>(term.l) + 0.5 * static_cast<double>(j1x2) + 0.5 * static_cast<double>(j2x2) - 0.5 * static_cast<double>(term.two_s);
        term.coefficient *= IntegerPhaseSign(exponent, "Amplitude initialization: GP leg exchange");
      }
    }
    if (!hc.alpha_ls.IsFinite()) { throw std::invalid_argument("Amplitude initialization: GP g_ls has non-finite coupling"); }
    const auto &pole1 = legs[0].pdg == PDG::PDG_gamma ? pdg_table.FindByPDG(PDG::PDG_gamma) : regge::PoleRepresentative(regge_param, pdg_table, legs[0].pdg);
    const auto &pole2 = legs[1].pdg == PDG::PDG_gamma ? pdg_table.FindByPDG(PDG::PDG_gamma) : regge::PoleRepresentative(regge_param, pdg_table, legs[1].pdg);
    // The complete canonical pole basis enforces angular momentum, P, C and identical-leg symmetry
    spin::ValidateCanonicalPoleTerms(mother, pole1, pole2,
        std::vector<spin::LSTerm>(channel.g_ls.cbegin(), channel.g_ls.cend()), true, channel.C_symmetry, channel.P_symmetry);
    hc.alpha_ls.RemoveBelow(coupling_min);
    // Retain validation of zero rows before omitting this production channel
    if (hc.alpha_ls.Empty()) { return hc; }
    gpom::InitResonanceLS(hc, mother.spinX2 / 2, pole1.spinX2, pole2.spinX2);
    return hc;
  }

  hc.coupling_basis = gra::CouplingBasis::Helicity;
  if (channel.helicity.empty() || channel.helicity.size() != channel.g_helicity.size()) {
    throw std::invalid_argument(
        "Amplitude initialization: GP helicity labels and couplings "
        "have inconsistent dimensions");
  }
  const std::size_t nm                 = static_cast<std::size_t>(2 * MMAX + 1);
  hc.T                                 = MMatrix<std::complex<double>>(nm, nm, 0.0);
  hc.T_set                             = MMatrix<bool>(nm, nm, false);
  const auto &leg1                     = pdg_table.FindByPDG(legs[0].pdg);
  const auto &leg2                     = pdg_table.FindByPDG(legs[1].pdg);
  const bool  identical                = leg1.pdg == leg2.pdg;
  const int   parity_phase             = hc.P_symmetry ? GPParityPhase(mother, {&leg1, &leg2}, regge_param) : 1;
  const int   identical_exchange_phase = identical ? IntegerPhaseSign(0.5 * static_cast<double>(j1x2 + j2x2 - mother.spinX2), "Amplitude initialization: GP identical-leg exchange") : 1;
  // Transport direct helicities through the same daughter exchange
  const int                     exchange_phase = reverse ? IntegerPhaseSign(0.5 * static_cast<double>(j1x2 + j2x2 - mother.spinX2), "Amplitude initialization: GP helicity leg exchange") : 1;
  std::set<std::pair<int, int>> input_orbits;
  const auto                    add = [&](int m1, int m2, std::complex<double> value, const std::string &origin) {
    if (2 * std::abs(m1) > j1x2 || 2 * std::abs(m2) > j2x2) {
      throw std::invalid_argument(
                             "Amplitude initialization: GP helicity row exceeds the "
                                                "exchange pole spin");
    }
    const std::size_t i1 = gpom::AnalyticMIndex(m1, MMAX, "PrepareResonance");
    const std::size_t i2 = gpom::AnalyticMIndex(m2, MMAX, "PrepareResonance");
    const long long   mu = static_cast<long long>(m1) - m2;
    if (std::abs(mu) > mother.spinX2 / 2) {
      throw std::invalid_argument(
                             "Amplitude initialization: GP helicity row violates "
                                                "|m1-m2| <= JX");
    }
    if ((legs[0].pdg == PDG::PDG_gamma && std::abs(m1) != 1) || (legs[1].pdg == PDG::PDG_gamma && std::abs(m2) != 1)) { throw std::invalid_argument("Amplitude initialization: GP photon helicity must be transverse"); }
    if (hc.T_set[i1][i2]) {
      const double scale = std::max({1.0, std::abs(hc.T[i1][i2]), std::abs(value)});
      if (std::abs(hc.T[i1][i2] - value) > kComparisonTolerance * scale) {
        throw std::invalid_argument(
                               "Amplitude initialization: GP helicity symmetry orbit "
                                                  "is inconsistent at " +
                               origin);
      }
      return;
    }
    hc.T[i1][i2]     = value;
    hc.T_set[i1][i2] = true;
  };
  for (const auto &i : indices(channel.helicity)) {
    int m1 = static_cast<int>(std::llround(channel.helicity[i][0]));
    int m2 = static_cast<int>(std::llround(channel.helicity[i][1]));
    if (reverse) { std::swap(m1, m2); }
    if (hc.P_symmetry || identical) {
      std::pair<int, int> orbit = std::make_pair(m1, m2);
      if (hc.P_symmetry) { orbit = std::min(orbit, std::make_pair(-m1, -m2)); }
      if (identical) {
        orbit = std::min(orbit, std::make_pair(m2, m1));
        if (hc.P_symmetry) { orbit = std::min(orbit, std::make_pair(-m2, -m1)); }
      }
      if (!input_orbits.insert(orbit).second) {
        throw std::invalid_argument(
            "Amplitude initialization: GP helicity rows repeat one "
            "symmetry orbit");
      }
    }
    const std::complex<double> coupling = static_cast<double>(exchange_phase) * channel.g_helicity[i];
    add(m1, m2, coupling, "input row");
    if (hc.P_symmetry) { add(-m1, -m2, static_cast<double>(parity_phase) * coupling, "parity partner"); }
    if (identical) {
      const std::complex<double> exchanged = static_cast<double>(identical_exchange_phase) * coupling;
      add(m2, m1, exchanged, "identical-leg exchange partner");
      if (hc.P_symmetry) { add(-m2, -m1, static_cast<double>(parity_phase) * exchanged, "identical-leg parity-exchange partner"); }
    }
  }
  if (hc.P_symmetry) { ValidateRawSymmetry(hc, parity_phase, SpinSymmetry::Parity); }
  if (identical) { ValidateRawSymmetry(hc, identical_exchange_phase, SpinSymmetry::Exchange); }
  const auto &pole1 = legs[0].pdg == PDG::PDG_gamma ? leg1 : regge::PoleRepresentative(regge_param, pdg_table, legs[0].pdg);
  const auto &pole2 = legs[1].pdg == PDG::PDG_gamma ? leg2 : regge::PoleRepresentative(regge_param, pdg_table, legs[1].pdg);
  ValidateGPPoleSubspace(hc, mother, pole1, pole2, reverse, exchange_phase);
  PruneHelicityCouplings(hc, coupling_min, "Amplitude initialization: GP helicity");
  return hc;
}

// Normalize analytic diphoton couplings at the physical transverse pole
void NormalizeGammaGamma(HELMatrix &hel, const PARAM_RES &res, const std::vector<MParticle> &legs,
                         bool derivative_factor, bool isolated_decay) {
  double norm2 = 0.0;
  if (hel.UsesHelicityCouplings()) {
    norm2 = hel.T.MaskedSquaredNorm(hel.T_set);
  } else {
    HELMatrix physical = hel;
    spin::ApplyCanonicalPoleLSCoefficients(physical, res.p, legs[0], legs[1]);
    const double q_ratio = derivative_factor ? 0.5 * res.p.mass / hel.analytic_Lambda : 1.0;
    for (auto &term : physical.alpha_ls) { term.coefficient *= math::IntegerPower(q_ratio, static_cast<unsigned int>(term.l)); }
    spin::InitTMatrix(physical, res.p, legs[0], legs[1], true, "GP diphoton transverse normalization", false, false);
    norm2 = physical.T.FrobNorm2();
  }
  const double scale = GammaGammaScale(res, norm2, isolated_decay);
  if (!hel.UsesHelicityCouplings()) { hel.alpha_ls.Scale(scale); }
  hel.T *= scale;
}

// Build one ordered GP continuum pair in its analytic trajectory basis
ReggeContinuumPole PreparePair(MProcessSetup &setup, const regge::Param &param, const MParticle &upper_exchange, const MParticle &lower_exchange, const MParticle &first, const MParticle &second, const std::string &context) {
  for (const auto *exchange : {&upper_exchange, &lower_exchange}) {
    if (exchange->pdg != PDG::PDG_gamma && exchange->spinX2 != gra::aux::kNullSpinX2) {
      throw std::invalid_argument(
          "Amplitude initialization: GP continuum requires null-spin "
          "analytic exchange aliases or photons");
    }
  }
  ReggeContinuumPole out;
  out.gp_vertex[0] = setup.ProcessHelicityStructure(upper_exchange, {first, second}, true, true, context + " upper", true, spin::VertexContext::SubTUChannelExchange);
  out.gp_vertex[1] = setup.ProcessHelicityStructure(lower_exchange, {second, first}, true, true, context + " lower", true, spin::VertexContext::SubTUChannelExchange);
  if (out.gp_vertex[0].UsesLSCouplings()) { gpom::InitCrossedLS(out.gp_vertex[0], gpom::PoleSpinX2(param, setup.lts.PDG, upper_exchange.pdg) / 2); }
  if (out.gp_vertex[1].UsesLSCouplings()) { gpom::InitCrossedLS(out.gp_vertex[1], gpom::PoleSpinX2(param, setup.lts.PDG, lower_exchange.pdg) / 2); }
  const bool upper_photon = upper_exchange.pdg == PDG::PDG_gamma;
  const bool lower_photon = lower_exchange.pdg == PDG::PDG_gamma;
  if (upper_photon != lower_photon && setup.lts.process.MMAX < 1) {
    throw std::invalid_argument("Amplitude initialization: mixed photon GP continuum requires MMAX >= 1");
  }
  return out;
}

// Build the four ordered analytic GP continuum subvertices
std::vector<HELMatrix> PrepareContinuum(MProcessSetup &setup, const std::vector<MDecayBranch> &production_tree, const regge::Param &param) {
  if (production_tree.size() != 2 || setup.lts.decaytree.size() != 2) { throw std::invalid_argument("Amplitude initialization: GP continuum expects two branches"); }
  const auto &central  = setup.lts.decaytree;
  auto        t        = PreparePair(setup, param, production_tree[0].p, production_tree[1].p, central[0].p, central[1].p, "GP continuum t");
  auto        u        = PreparePair(setup, param, production_tree[0].p, production_tree[1].p, central[1].p, central[0].p, "GP continuum u");
  return {std::move(t.gp_vertex[0]), std::move(t.gp_vertex[1]), std::move(u.gp_vertex[0]), std::move(u.gp_vertex[1])};
}

// Retain m=0 in the same GP representation used by two-particle continuum vertices
void PrepareLadder(HELMatrix &hel, const std::string &context) {
  const std::size_t zero = gpom::AnalyticMIndex(0, hel.analytic_MMAX, context);
  for (const auto &column : indices(hel.Jz_values)) {
    if (column == zero) { continue; }
    const bool active = hel.UsesLSCouplings()
                            ? !hel.m_ls[column].Empty()
                            : std::any_of(hel.T_active.cbegin(), hel.T_active.cend(),
                                          [column](const auto &entry) { return entry.second == column; });
    if (active) { throw std::invalid_argument(context + ": finite nonzero-m GP ladder transport is not defined"); }
  }
  hel.T             = hel.T.SelectColumns({zero});
  hel.T_set         = hel.T_set.SelectColumns({zero});
  hel.analytic_MMAX = 0;
  hel.Jz_values     = {0.0};
  if (hel.UsesLSCouplings()) {
    auto terms = std::move(hel.m_ls[zero]);
    hel.m_ls   = {std::move(terms)};
    gpom::InitCrossedLS(hel, static_cast<int>(std::llround(hel.J)));
  } else {
    gpom::CheckHelicity(hel, context, false);
  }
}

}  // namespace gra::gpom
