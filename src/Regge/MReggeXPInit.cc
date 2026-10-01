// XP production vertex preparation
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Regge/MReggeXPInit.h"
#include "Graniitti/Regge/MReggeXP.h"
#include "Graniitti/Process/MHelicityConfig.h"
#include "Graniitti/Tech/MAux.h"
#include "rang.hpp"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <utility>

using gra::aux::indices;

namespace gra::xpom {

// Format one topology coupling for diagnostics
std::string FormatCoupling(const std::complex<double> coupling) {
  std::ostringstream ss;
  ss << std::fixed << std::setprecision(4) << "[mag=" << std::abs(coupling) << ", phase=" << std::arg(coupling) << "]";
  return ss.str();
}

// Print one complete canonical XP pole operator and its helicity projection
void PrintPole(const std::string &label, const MParticle &mother, const std::vector<MParticle> &legs,
                         const spin::PoleLS &vertex) {
  if (legs.size() != 2 || vertex.terms.size() != vertex.raw_normalization.size()) { throw std::invalid_argument("Amplitude initialization: invalid XP pole diagnostic metadata"); }
  std::cout << rang::fg::yellow << label << ": " << mother.pdg << " < [" << legs[0].pdg << ", " << legs[1].pdg << "]" << rang::fg::reset << std::endl;
  std::vector<std::vector<std::string>> rows;
  rows.reserve(vertex.terms.size());
  for (const auto &i : indices(vertex.terms)) {
    const auto &term = vertex.terms[i];
    rows.push_back({std::to_string(term.l), gra::aux::Spin2XtoString(static_cast<int>(term.two_s)), FormatCoupling(term.coefficient), std::to_string(vertex.raw_normalization[i]),
                    vertex.context == spin::VertexContext::SubTUChannelExchange ? "DL reduced residue" : "(p/Lambda)^" + std::to_string(term.l)});
  }
  gra::aux::PrintTable({"L", "S", "g_LS", "raw STF", "power"}, rows);
  const auto reduced = spin::PoleLSReduced(vertex, vertex.Lambda);
  rows.clear();
  rows.reserve(vertex.helicity.lambda_values.size_row());
  for (std::size_t row = 0; row < vertex.helicity.lambda_values.size_row(); ++row) {
    const auto i1 = vertex.helicity.lambda_idx[row][0];
    const auto i2 = vertex.helicity.lambda_idx[row][1];
    rows.push_back({std::to_string(vertex.helicity.lambda_values[row][0]), std::to_string(vertex.helicity.lambda_values[row][1]), FormatCoupling(reduced[i1][i2])});
  }
  gra::aux::PrintTable({"lambda1", "lambda2", "H(lambda1,lambda2)"}, rows);
}

// Prepare one XP fusion vertex from absolute LS or helicity couplings
spin::PoleLS PrepareResonance(MProcessSetup &setup, const PARAM_RES &res, const RES_PRODUCTION_CHANNEL &channel,
                              const std::vector<MParticle> &legs) {
  if (!IsAutomaticReggeVertexBasis(channel.basis)) { return rspin::PreparePole(setup, res, channel, legs); }
  auto canonical = legs;
  if (canonical.size() == 2 && canonical[0].pdg != channel.exchange[0]) { std::swap(canonical[0], canonical[1]); }
  const auto automatic = spin::BuildAutomaticCentralCoupling(res.p, canonical, rspin::ReggeAutoMode(channel.basis), true,
      "XP " + std::string(ReggeVertexBasisName(channel.basis, ReggeVertexRole::Resonance)) + " auto-selection", "", true,
      channel.P_symmetry, channel.C_symmetry, spin::VertexContext::Auto);
  return rspin::PreparePole(setup, res, channel, legs, &automatic);
}

// Prepare the ordered XP continuum vertices in the common pole normalization
std::vector<spin::PoleResidue> PrepareContinuum(MProcessSetup &setup, const std::vector<MDecayBranch> &tree) {
  for (const auto &branch : setup.lts.decaytree) { ClassifyPole(branch.p); }
  return rspin::PreparePoleContinuum(setup, tree, ReggeProductionModel::XP);
}

// Prepare one local XP pair for multiparticle production
ReggeContinuumPole PreparePair(MProcessSetup &setup, const MParticle &upper, const MParticle &lower,
                               const MParticle &first, const MParticle &second) {
  ClassifyPole(first);
  ClassifyPole(second);
  return rspin::PreparePolePair(setup, ReggeProductionModel::XP, upper, lower, first, second);
}

}  // namespace gra::xpom
