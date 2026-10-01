// Shared finite spin amplitudes for the MP and XP models
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// For Pomeron spin discussions, see:
//
// [REFERENCE: Close, Schuler, https://arxiv.org/abs/hep-ph/9902243v1]
// [REFERENCE: Close, Schuler, https://arxiv.org/abs/hep-ph/9905305]
// [REFERENCE: Kaidalov, KMR, https://arxiv.org/abs/hep-ph/0307064]

// C++
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

// Own
#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/Regge/MReggeMPXP.h"
#include "Graniitti/Regge/MReggeInit.h"
#include "Graniitti/Process/MHelicityConfig.h"
#include "Graniitti/Spin/MHelicity.h"
#include "Graniitti/Spin/MHelicityScatter.h"
#include "Graniitti/Spin/MPoleLS.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MException.h"

using gra::aux::indices;

namespace gra::rspin {

// Compute the pole-normalized MP/XP photon transition coupling
double PhotoCoupling(const RES_PRODUCTION& production) {
  if (!production.pole) { throw AmplitudeFailure("MReggeMPXP::PhotoCoupling: missing physical-pole LS vertex"); }
  const double density = spin::LeadingPoleDensity(*production.pole, production.pole->Lambda);
  if (!std::isfinite(density) || !(density > 0.0)) {
    throw AmplitudeFailure("MReggeMPXP::PhotoCoupling: invalid transition");
  }
  return std::sqrt(density);
}

// Compute the local MP or XP continuum contraction basis
std::vector<double> LocalBasis(const ReggeContinuumPole& pole, const std::size_t leg) {
  const auto& vertex = pole.pole_operator[leg].Pole();
  return vertex.helicity.exchange_basis == ExchangeBasisType::HelicityTransport ? vertex.helicity.Jz_values
                                                                                : std::vector<double>{0.0};
}

namespace {

// Compute the initial forward helicity dimension of one two-leg channel
std::size_t InitialDim(const std::vector<MDecayBranch>& tree, bool drop_flip) {
  return spin::Rows(tree[0], drop_flip).size() * spin::Rows(tree[1], drop_flip).size();
}

// Boost the two continuum particles once into their common X rest frame
std::array<M4Vec, 2> ContinuumFinalPairInX(const LORENTZSCALAR& lts) {
  return {kinematics::BoostToRestFrame(lts.decaytree[0].p4, lts.pfinal[0], "MReggeMPXP::ContinuumFinalPairInX"),
          kinematics::BoostToRestFrame(lts.decaytree[1].p4, lts.pfinal[0], "MReggeMPXP::ContinuumFinalPairInX")};
}

// Evaluate a continuum pole in its reduced Regge or physical photon basis
const spin::EvaluatedPoleSubvertex& ContinuumPoleVertex(const spin::PoleResidue& vertex, const M4Vec& final_in_X,
                                                        const M4Vec& parent_dir_in_X, bool second_exchange_daughter,
                                                        spin::EvaluatedPoleSubvertex& local) {
  if (vertex.Pole().helicity.exchange_basis == ExchangeBasisType::ReducedRegge) { return vertex.Reduced(); }
  local = spin::EvaluatePoleSubvertex(vertex.Reduced().helicity, final_in_X, parent_dir_in_X, second_exchange_daughter);
  return local;
}

}  // namespace

// Compute the pole-normalized t and u scales of four ordered subvertices
PoleScales ContinuumPoleScales(const std::vector<spin::PoleResidue>& pole) {
  return {std::sqrt(pole[0].Density() * pole[1].Density()), std::sqrt(pole[2].Density() * pole[3].Density())};
}

// Build per-channel finite-spin resonance production matrices
std::vector<HelAmp> Resonance(const LORENTZSCALAR& lts, const PARAM_RES& res, FusionTensor fusion, double s0,
                              const HelAmp* spin_filter, const std::array<ForwardLegState, 2>* forward_state) {
  std::vector<HelAmp> out;
  out.reserve(res.production.size());
  const spin::ForwardSpec forward = {lts.process.FORWARD_VERTEX, s0};
  if (!lts.process.SPINGEN) {
    const double      momentum = lts.q1_in_X.P3mod();
    const std::size_t spin_states =
        spin::SpinStateCount(res.p.spinX2 / 2.0, "MReggeMPXP::Resonance resonance spin");
    for (const auto& production : res.production) {
      const double      pole_scale     = std::sqrt(spin::LeadingPoleDensity(*production.pole, momentum));
      const std::size_t initial_states = InitialDim(production.tree, lts.process.FORWARD_NOFLIP);
      out.push_back(spin::Blind(initial_states, spin_states, pole_scale));
    }
    return out;
  }

  if (!kinematics::FiniteFourVector(lts.q1_in_X) || !kinematics::FiniteFourVector(lts.q2_in_X)) {
    throw AmplitudeFailure("MReggeMPXP::Resonance: non-finite generated exchange momentum");
  }
  if (forward_state != nullptr &&
      (forward_state->at(0).leg != ForwardBeamLeg::Upper || forward_state->at(1).leg != ForwardBeamLeg::Lower)) {
    throw AmplitudeFailure("MReggeMPXP::Resonance: invalid forward state ordering");
  }
  const M4Vec& upper_in  = forward_state != nullptr ? forward_state->at(0).incoming : lts.pbeam1;
  const M4Vec& upper_out = forward_state != nullptr ? forward_state->at(0).outgoing : lts.pfinal[1];
  const M4Vec& lower_in  = forward_state != nullptr ? forward_state->at(1).incoming : lts.pbeam2;
  const M4Vec& lower_out = forward_state != nullptr ? forward_state->at(1).outgoing : lts.pfinal[2];
  for (const auto& production : res.production) {
    const auto& tree          = production.tree;
    const auto  up_rows       = spin::Rows(tree[0], lts.process.FORWARD_NOFLIP);
    const auto  dn_rows       = spin::Rows(tree[1], lts.process.FORWARD_NOFLIP);
    const auto  upper_factors = spin::ForwardFactors(tree[0], upper_in, upper_out, false, up_rows, forward);
    const auto  lower_factors = spin::ForwardFactors(tree[1], lower_in, lower_out, true, dn_rows, forward);
    const auto& vertex        = *production.pole;
    HelAmp      central       = fusion(lts, vertex);
    if (spin_filter != nullptr) { central = central * *spin_filter; }
    if (upper_factors.has_value() && lower_factors.has_value()) {
      out.push_back(spin::Contract(*upper_factors, *lower_factors, central, 1.0));
      continue;
    }

    const auto upper =
        spin::Forward(lts, tree[0], upper_in, upper_out, false, up_rows, lts.process.PHOTON_VERTEX, forward);
    const auto lower =
        spin::Forward(lts, tree[1], lower_in, lower_out, true, dn_rows, lts.process.PHOTON_VERTEX, forward);
    out.push_back(spin::Contract(upper, lower, central, 1.0));
  }
  return out;
}

// Build each continuum channel as sources, pole projection, then proton rows
std::vector<HelPair> Continuum(const LORENTZSCALAR& lts, double s0, PoleMetric numerator) {
  const auto& left  = lts.decaytree[0].p;
  const auto& right = lts.decaytree[1].p;
  const auto  n_final =
      spin::FinalStateHelicityCount(left, "Continuum left") * spin::FinalStateHelicityCount(right, "Continuum right");
  const spin::ForwardSpec forward = {lts.process.FORWARD_VERTEX, s0, ExchangeBasisType::ReducedRegge};
  std::vector<HelPair>    out;
  out.reserve(lts.process.CONT_PRODUCTIONTREE.size());
  const auto final   = lts.process.SPINGEN ? ContinuumFinalPairInX(lts) : std::array<M4Vec, 2>{};
  M4Vec      q2_axis = lts.q2_in_X;
  q2_axis.Flip3();
  for (const auto channel : indices(lts.process.CONT_PRODUCTIONTREE)) {
    const auto& tree = lts.process.CONT_PRODUCTIONTREE[channel];
    const auto& pole = lts.process.CONTINUUM_POLE[channel];
    if (!lts.process.SPINGEN) {
      const auto scale     = ContinuumPoleScales(pole);
      const auto n_initial = InitialDim(tree, lts.process.FORWARD_NOFLIP);
      out.emplace_back(spin::Blind(n_initial, n_final, scale.t), spin::Blind(n_initial, n_final, scale.u));
      continue;
    }
    const auto up_rows = spin::Rows(tree[0], lts.process.FORWARD_NOFLIP);
    const auto dn_rows = spin::Rows(tree[1], lts.process.FORWARD_NOFLIP);
    const auto up      = spin::ForwardFactors(tree[0], lts.pbeam1, lts.pfinal[1], false, up_rows, forward);
    const auto dn      = spin::ForwardFactors(tree[1], lts.pbeam2, lts.pfinal[2], true, dn_rows, forward);

    // Evaluate both final-particle orders with the same projection and sewing
    const auto project = [&](const auto& upper_source, const auto& lower_source) {
      using Matrix = decltype(spin::ProjectedSubchannel(pole[0].Reduced(), pole[1].Reduced(), upper_source,
                                                        lower_source, left, right, false));
      std::array<Matrix, 2> sub;
      for (const auto order : indices(sub)) {
        std::array<spin::EvaluatedPoleSubvertex, 2> local;
        const auto& upper = ContinuumPoleVertex(pole[2 * order], final[order], lts.q1_in_X, false, local[0]);
        const auto& lower = ContinuumPoleVertex(pole[2 * order + 1], final[1 - order], q2_axis, true, local[1]);
        const auto  metric =
            numerator ? numerator(lts.decaytree[order].p, lts.q1_in_X - final[order]) : spin::InternalHelicityMetric{};
        sub[order] = spin::ProjectedSubchannel(upper, lower, upper_source, lower_source, left, right, order == 1, false,
                                               numerator ? &metric : nullptr);
      }
      return std::pair{std::move(sub[0]), std::move(sub[1])};
    };
    if (up && dn) {
      const auto sub = project(up->exchange_helicity, dn->exchange_helicity);
      out.push_back(spin::Contract(*up, *dn, sub.first, sub.second));
    } else {
      out.push_back(project(
          spin::Forward(lts, tree[0], lts.pbeam1, lts.pfinal[1], false, up_rows, lts.process.PHOTON_VERTEX, forward),
          spin::Forward(lts, tree[1], lts.pbeam2, lts.pfinal[2], true, dn_rows, lts.process.PHOTON_VERTEX, forward)));
    }
  }
  return out;
}

// Evaluate one local final-pair kernel from the two prepared subvertices
HelAmp PairKernel(const ReggeContinuumPole& pole, const MDecayBranch& first, const MDecayBranch& second,
                  const M4Vec& upper_exchange, const M4Vec& lower_exchange) {
  std::array<spin::EvaluatedPoleSubvertex, 2> local;
  M4Vec                                       first_rest, second_rest, upper_rest, lower_rest;
  const bool                                  physical_pole =
      pole.pole_operator[0].Pole().helicity.exchange_basis == ExchangeBasisType::HelicityTransport ||
      pole.pole_operator[1].Pole().helicity.exchange_basis == ExchangeBasisType::HelicityTransport;
  if (physical_pole) {
    const M4Vec total = first.p4 + second.p4;
    first_rest        = kinematics::BoostToRestFrame(first.p4, total, "MReggeMPXP::PairKernel first");
    second_rest       = kinematics::BoostToRestFrame(second.p4, total, "MReggeMPXP::PairKernel second");
    upper_rest        = kinematics::BoostToRestFrame(upper_exchange, total, "MReggeMPXP::PairKernel upper");
    lower_rest        = kinematics::BoostToRestFrame(lower_exchange, total, "MReggeMPXP::PairKernel lower");
  }
  const auto&       upper      = ContinuumPoleVertex(pole.pole_operator[0], first_rest, upper_rest, false, local[0]);
  const auto&       lower      = ContinuumPoleVertex(pole.pole_operator[1], second_rest, lower_rest, true, local[1]);
  const std::size_t upper_dim  = upper.helicity.Jz_values.size();
  const std::size_t lower_dim  = lower.helicity.Jz_values.size();
  auto              subchannel = spin::Subchannel(upper, lower, first.p, second.p, false, true);
  subchannel.Reshape(upper_dim, lower_dim);
  return subchannel;
}

// Prepare the shared finite-spin pole coefficients before event generation

// Normalize one physical diphoton pole using its absolute LS coefficients
void NormalizeGammaGamma(spin::PoleLS &vertex, const PARAM_RES &res, bool isolated_decay) {
  const double scale = GammaGammaScale(res, spin::PoleLSReduced(vertex, 0.5 * res.p.mass).FrobNorm2(), isolated_decay);
  for (auto &term : vertex.terms) { term.coefficient *= scale; }
  vertex.helicity.alpha_ls.Scale(scale);
  vertex.helicity.T *= scale;
}

// Compute whether one ordered topology reverses its typed production channel
bool ReversesChannel(const RES_PRODUCTION_CHANNEL &channel, const std::vector<MParticle> &legs, const std::string &context) {
  if (legs.size() != 2) { throw std::invalid_argument(context + " expects two production legs"); }
  const std::array<int, 2> pair = {legs[0].pdg, legs[1].pdg};
  if (pair == channel.exchange) { return false; }
  if (pair[0] == channel.exchange[1] && pair[1] == channel.exchange[0]) { return true; }
  throw std::invalid_argument(context + " topology does not match its production channel");
}

// Order canonical raw STF LS terms with the phase (-1)^(L+j1+j2-S)
std::vector<spin::LSTerm> OrderedLSTerms(const RES_PRODUCTION_CHANNEL &channel, const std::vector<MParticle> &legs,
                                       std::vector<spin::LSTerm> terms) {
  if (ReversesChannel(channel, legs, "Amplitude initialization: fixed-pole leg exchange")) {
    for (auto &row : terms) {
      const double exponent = static_cast<double>(row.l) + 0.5 * legs[0].spinX2 + 0.5 * legs[1].spinX2 - 0.5 * row.two_s;
      row.coefficient *= IntegerPhaseSign(exponent, "Amplitude initialization: fixed-pole leg exchange");
    }
  }
  return terms;
}

// Map one generated basis into the shared spin builder mode
gra::spin::AutoCentralCouplingMode ReggeAutoMode(ReggeVertexBasis basis) {
  if (basis == ReggeVertexBasis::AutoMinL) { return gra::spin::AutoCentralCouplingMode::MinL; }
  if (basis == ReggeVertexBasis::AutoMinS) { return gra::spin::AutoCentralCouplingMode::MinS; }
  if (basis == ReggeVertexBasis::AutoEqualLS) { return gra::spin::AutoCentralCouplingMode::FlatLS; }
  if (basis == ReggeVertexBasis::AutoEqualHelicity) { return gra::spin::AutoCentralCouplingMode::FlatHelicity; }
  throw std::invalid_argument("Amplitude initialization: invalid automatic resonance basis");
}

// Convert one generated helicity matrix into canonical raw pole coefficients
std::vector<spin::LSTerm> AutomaticPoleTerms(const HELMatrix &helicity, const MParticle &mother, const std::vector<MParticle> &legs, const std::complex<double> coupling, gra::spin::VertexContext context, bool check = false) {
  std::vector<spin::LSTerm> terms;
  if (!helicity.UsesHelicityCouplings()) {
    terms.reserve(helicity.alpha_ls.Size());
    for (const spin::LSTerm &term : helicity.alpha_ls) { terms.push_back({term.l, term.two_s, term.coefficient * coupling}); }
    return terms;
  }

  const auto expansion = gra::spin::DirectHelicityToLSCoefficients(helicity, mother, legs[0], legs[1], true, context);
  if (check && expansion.residual_norm2 > 1e-18) { throw std::invalid_argument("Amplitude initialization: resonance vertex is outside the canonical LS subspace"); }
  terms.reserve(expansion.rows.size());
  for (const auto &i : indices(expansion.rows)) {
    const auto                &row         = expansion.rows[i];
    const std::complex<double> coefficient = expansion.coefficients[i] / gra::spin::RawPoleLSNormalization(legs[0].spinX2, legs[1].spinX2, mother.spinX2, row.l, static_cast<int>(row.two_s));
    terms.push_back({row.l, row.two_s, coefficient * coupling});
  }
  return terms;
}

// Convert one typed fixed-pole channel into ordered canonical LS terms
std::vector<spin::LSTerm> PoleTerms(MProcessSetup &setup, const PARAM_RES &RES, const RES_PRODUCTION_CHANNEL &channel, const std::vector<MParticle> &legs) {
  if (legs.size() != 2) { throw std::invalid_argument("Amplitude initialization: pole production expects two legs"); }
  if (channel.basis == ReggeVertexBasis::LS) { return OrderedLSTerms(channel, legs, {channel.g_ls.cbegin(), channel.g_ls.cend()}); }
  if (channel.basis != ReggeVertexBasis::Helicity) { throw std::invalid_argument("Amplitude initialization: unsupported fixed-pole resonance basis"); }
  const auto        vertex_legs = JacobWickVertexLegs(legs, true, gra::spin::VertexContext::Auto, setup.lts.PDG);
  const std::string context     = std::string("Regge resonance vertex");
  const auto        helicity    = BuildDirectHelicityCoupling(RES.p, vertex_legs, channel.helicity, channel.g_helicity, channel.C_symmetry, channel.P_symmetry, ReversesChannel(channel, vertex_legs, context), context);
  return AutomaticPoleTerms(helicity, RES.p, legs, 1.0, spin::VertexContext::Auto, true);
}

// Prepare a finite-spin pole with the common absolute LS normalization
spin::PoleLS PreparePole(MProcessSetup &setup, const PARAM_RES &res, const RES_PRODUCTION_CHANNEL &channel,
                         const std::vector<MParticle> &legs, const HELMatrix *automatic) {
  if (!setup.model_tune) { throw std::logic_error("Amplitude initialization: model tune is not configured"); }
  auto canonical = legs;
  if (ReversesChannel(channel, legs, "Amplitude initialization: fixed-pole production")) { std::swap(canonical[0], canonical[1]); }
  const auto terms = automatic ? OrderedLSTerms(channel, legs, AutomaticPoleTerms(*automatic, res.p, canonical, channel.g, spin::VertexContext::Auto))
                               : PoleTerms(setup, res, channel, legs);
  if (!automatic) { spin::ValidateCanonicalPoleTerms(res.p, legs[0], legs[1], terms, true, channel.C_symmetry, channel.P_symmetry); }
  return spin::PreparePoleLS(res.p, legs[0], legs[1], terms, channel.Lambda, true, channel.C_symmetry,
                             channel.P_symmetry, spin::VertexContext::Auto, setup.model_tune->Global().coupling_min, setup.lts.process.DERIVATIVE_FACTOR);
}

// Build one ordered MP or XP pair from its representative pole coefficient
ReggeContinuumPole PreparePolePair(MProcessSetup &setup, const ReggeProductionModel model, const MParticle &upper_exchange, const MParticle &lower_exchange, const MParticle &first, const MParticle &second) {
  ReggeContinuumPole out;
  out.pole_operator[0] = spin::PoleResidue(setup.ProcessPoleOperatorStructure(
      model, upper_exchange, {first, second}, spin::VertexContext::SubTUChannelExchange));
  out.pole_operator[1] = spin::PoleResidue(setup.ProcessPoleOperatorStructure(
      model, lower_exchange, {second, first}, spin::VertexContext::SubTUChannelExchange));
  return out;
}

// Build the four ordered MP or XP canonical pole subvertices
std::vector<spin::PoleResidue> PreparePoleContinuum(MProcessSetup                   &setup,
                                                  const std::vector<MDecayBranch> &production_tree,
                                                  ReggeProductionModel model) {
  if (production_tree.size() != 2 || setup.lts.decaytree.size() != 2) { throw std::invalid_argument("Amplitude initialization: continuum pole expects two branches"); }
  if (model != ReggeProductionModel::MP && model != ReggeProductionModel::XP) { throw std::invalid_argument("Amplitude initialization: continuum pole requires MP or XP"); }
  const auto &central = setup.lts.decaytree;
  auto        t       = PreparePolePair(setup, model, production_tree[0].p, production_tree[1].p, central[0].p, central[1].p);
  auto        u       = PreparePolePair(setup, model, production_tree[0].p, production_tree[1].p, central[1].p, central[0].p);
  return {std::move(t.pole_operator[0]), std::move(t.pole_operator[1]), std::move(u.pole_operator[0]), std::move(u.pole_operator[1])};
}

}  // namespace gra::rspin
