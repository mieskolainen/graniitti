// MP production vertex preparation
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Regge/MReggeMPInit.h"

#include <cmath>
#include <iostream>
#include <utility>

#include "Graniitti/Process/MHelicityConfig.h"
#include "Graniitti/Regge/MReggeMPXP.h"
#include "Graniitti/Spin/MSpinDensity.h"
#include "Graniitti/Tech/MAux.h"
#include "rang.hpp"

using gra::aux::indices;

namespace gra::mpom {

// Fill a missing resonance name from the PDG table when available
MParticle NamedMPMother(const MProcessSetup& setup, const MParticle& mother) {
  MParticle named_mother = mother;
  if (!named_mother.name.empty()) { return named_mother; }
  try {
    named_mother.name = setup.lts.PDG.FindByPDG(named_mother.pdg).name;
  } catch (...) {
    // Keep the resonance-card particle unchanged if no PDG-table name exists
  }
  return named_mother;
}

// Compute the spin-projection label for one vector or matrix index
std::string SpinLabel(const std::size_t index, const std::size_t dimension) {
  return gra::spin::FormatSpinLabel(-0.5 * static_cast<double>(dimension - 1) + static_cast<double>(index));
}

// Print the active MP spin steering
void PrintPolarization(const PARAM_RES& res, const bool metrics_enabled) {
  constexpr double kZeroTolerance = 1e-12;
  std::cout << rang::fg::yellow << "MP spin basis: " << res.spin_basis << rang::fg::reset << std::endl;
  if (metrics_enabled && !res.UsesUnrestrictedSpinBasis()) {
    const auto metrics = spin::DensityMetrics(res.rho);
    std::cout << "MP steering density (input): purity = " << gra::aux::ToString(metrics.purity, 4)
              << ", entropy = " << gra::aux::ToString(metrics.entropy, 4) << " bits" << std::endl;
  }
  std::vector<std::vector<std::string>> rows;
  if (res.UsesCoherentSpinBasis()) {
    for (const auto& i : indices(res.a_Jz)) {
      rows.push_back({SpinLabel(i, res.a_Jz.size()), gra::aux::ToString(std::abs(res.a_Jz[i]), 3),
                      gra::aux::ToString(std::arg(res.a_Jz[i]), 3), gra::aux::ToString(std::norm(res.a_Jz[i]), 3)});
    }
    std::cout << "MP coherent spin steering a_Jz:" << std::endl;
    gra::aux::PrintTable({"Jz", "|a_Jz|", "arg(a_Jz)", "|a_Jz|^2"}, rows);
    std::cout << std::endl;
    return;
  }
  if (!res.UsesDensitySpinBasis()) {
    std::cout << "MP resonance spin state is determined by the production amplitude" << std::endl << std::endl;
    return;
  }
  for (std::size_t i = 0; i < res.rho.size_row(); ++i) {
    for (std::size_t j = 0; j < res.rho.size_col(); ++j) {
      if (std::abs(res.rho[i][j]) <= kZeroTolerance) { continue; }
      rows.push_back({SpinLabel(i, res.rho.size_row()), SpinLabel(j, res.rho.size_col()),
                      gra::aux::ToString(std::abs(res.rho[i][j]), 3), gra::aux::ToString(std::arg(res.rho[i][j]), 3)});
    }
  }
  std::cout << "MP spin density matrix rho:" << std::endl;
  gra::aux::PrintTable({"Jz(row)", "Jz(col)", "|rho|", "arg(rho)"}, rows);
  std::cout << "Tr(rho) = " << gra::aux::ToString(std::real(res.rho.Trace()), 6) << " + i "
            << gra::aux::ToString(std::imag(res.rho.Trace()), 6) << std::endl
            << std::endl;

}

// Prepare one MP fusion vertex from absolute LS or helicity couplings
spin::PoleLS PrepareResonance(MProcessSetup& setup, const PARAM_RES& res, const RES_PRODUCTION_CHANNEL& channel,
                              const std::vector<MParticle>& legs) {
  const auto basis = channel.basis;
  if (!IsAutomaticReggeVertexBasis(basis)) { return rspin::PreparePole(setup, res, channel, legs); }
  auto canonical = legs;
  if (canonical.size() == 2 && canonical[0].pdg != channel.exchange[0]) { std::swap(canonical[0], canonical[1]); }
  const auto automatic = spin::BuildAutomaticCentralCoupling(
      NamedMPMother(setup, res.p), JacobWickVertexLegs(canonical, true, spin::VertexContext::Auto, setup.lts.PDG),
      rspin::ReggeAutoMode(basis), true,
      "MP " + std::string(ReggeVertexBasisName(basis, ReggeVertexRole::Resonance)) + " auto-selection", "",
      true, channel.P_symmetry, channel.C_symmetry, spin::VertexContext::Auto);
  return rspin::PreparePole(setup, res, channel, legs, &automatic);
}

// Project a random density into the physical frame symmetries
void PrepareDensity(PARAM_RES& res, const LORENTZSCALAR& lts) {
  if (!res.MP.random_rho || !res.UsesDensitySpinBasis()) { return; }
  if (lts.process.MP_FRAME == "CM") { res.rho = HelAmp::DiagonalMatrix(res.rho.GetDiag()); }
  if (lts.beam1.pdg == lts.beam2.pdg) {
    HelAmp image;
    if (lts.process.MP_FRAME == "HX") {
      HelVec phase(res.rho.size_row());
      for (const auto i : indices(phase)) { phase[i] = i % 2 == 0 ? 1.0 : -1.0; }
      image = res.rho.LeftDiagonalProduct(phase).RightDiagonalProduct(phase);
    } else {
      std::vector<std::size_t> reversed(res.rho.size_row());
      for (const auto i : indices(reversed)) { reversed[i] = reversed.size() - 1 - i; }
      image = res.rho.SelectRows(reversed).SelectColumns(reversed);
    }
    // Averaging unitary images preserves positivity and unit trace
    res.rho = (res.rho + image) * 0.5;
  }
}

// Prepare S^T for production columns, with S = sqrt((2J+1) rho)
void PreparePolarization(PARAM_RES& res, const LORENTZSCALAR& lts) {
  res.MP.filter = {};
  if (res.UsesUnrestrictedSpinBasis()) { return; }
  if (res.UsesCoherentSpinBasis()) { res.rho = RankOneProjector(res.a_Jz); }
  PrepareDensity(res, lts);
  // Beam exchange rotates CS axes by pi about x and HX axes by pi about z
  // CM azimuthal symmetry requires diagonal rho in the fixed collider axes
  const auto& frame           = lts.process.MP_FRAME;
  const bool  identical_beams = lts.beam1.pdg == lts.beam2.pdg;
  const auto& rho             = res.rho;
  // Beam exchange constrains rho, since an overall phase of a_Jz cancels in its projector
  for (std::size_t i = 0; i < rho.size_row(); ++i) {
    for (std::size_t j = 0; j < rho.size_col(); ++j) {
      if (frame == "CM" && i != j && std::abs(rho[i][j]) > 1e-10) {
        throw std::invalid_argument("MP polarization in fixed CM axes requires diagonal rho for azimuthal symmetry");
      }
      const auto image = frame == "HX" ? ((i + j) % 2 == 0 ? 1.0 : -1.0) * rho[i][j]
                                       : rho[rho.size_row() - 1 - i][rho.size_col() - 1 - j];
      if (identical_beams && std::abs(rho[i][j] - image) > 1e-10) {
        throw std::invalid_argument("MP polarization requires a beam-exchange invariant rho for identical beams");
      }
    }
  }
  // A pure projector has an exact square root, avoiding numerical null eigenvalues
  const auto root = res.UsesCoherentSpinBasis() ? rho : rho.PrincipalPositiveSemidefiniteSquareRoot();
  res.MP.filter = root.Transpose() * std::sqrt(static_cast<double>(res.p.spinX2 + 1));
}

// Prepare the ordered MP continuum vertices in the common pole normalization
std::vector<spin::PoleResidue> PrepareContinuum(MProcessSetup& setup, const std::vector<MDecayBranch>& tree) {
  return rspin::PreparePoleContinuum(setup, tree, ReggeProductionModel::MP);
}

// Prepare one local MP pair for multiparticle production
ReggeContinuumPole PreparePair(MProcessSetup& setup, const MParticle& upper, const MParticle& lower,
                               const MParticle& first, const MParticle& second) {
  return rspin::PreparePolePair(setup, ReggeProductionModel::MP, upper, lower, first, second);
}

}  // namespace gra::mpom
