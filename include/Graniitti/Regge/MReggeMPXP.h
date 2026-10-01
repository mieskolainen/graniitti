// Shared finite spin amplitudes for the MP and XP models
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MREGGEMPXP_H
#define MREGGEMPXP_H

#include <array>
#include <complex>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include "Graniitti/Kinematics/MEventKinematics.h"
#include "Graniitti/Process/MProcessSetup.h"
#include "Graniitti/Regge/MReggeParam.h"
#include "Graniitti/Regge/MReggeProductionModel.h"
#include "Graniitti/Spin/MHelicityScatter.h"
#include "Graniitti/Spin/MPoleLS.h"

namespace gra::rspin {

// Compute the shared MP/XP photon pole strength
double PhotoCoupling(const RES_PRODUCTION& production);

// Store the pole-normalized continuum scales of the t and u subchannels
struct PoleScales {
  double t = 1.0;
  double u = 1.0;
};

// Compute the local MP or XP continuum contraction basis
std::vector<double> LocalBasis(const ReggeContinuumPole& pole, std::size_t leg);

// Select the model production tensor without importing its implementation
using FusionTensor = HelAmp (*)(const LORENTZSCALAR&, const spin::PoleLS&);

// Select an internal numerator, or null for the physical helicity sum
using PoleMetric = spin::InternalHelicityMetric (*)(const MParticle&, const M4Vec&);

// Build per-channel finite-spin resonance production matrices
std::vector<HelAmp> Resonance(const LORENTZSCALAR& lts, const PARAM_RES& res, FusionTensor fusion, double s0, const HelAmp* spin_filter = nullptr, const std::array<ForwardLegState, 2>* forward_state = nullptr);

// Build per-channel crossed Regge continuum production matrices
std::vector<HelPair> Continuum(const LORENTZSCALAR& lts, double s0, PoleMetric numerator = nullptr);

// Compute the pole-normalized t and u scales of four ordered subvertices
PoleScales ContinuumPoleScales(const std::vector<spin::PoleResidue>& pole);

// Evaluate one local final-pair kernel from the two prepared subvertices
HelAmp PairKernel(const ReggeContinuumPole& pole, const MDecayBranch& first, const MDecayBranch& second, const M4Vec& upper_exchange, const M4Vec& lower_exchange);

// Normalize a physical diphoton pole
void NormalizeGammaGamma(spin::PoleLS &vertex, const PARAM_RES &res, bool isolated_decay);

// Select the automatic LS basis
spin::AutoCentralCouplingMode ReggeAutoMode(ReggeVertexBasis basis);

// Prepare a pole with the common absolute LS normalization
spin::PoleLS PreparePole(MProcessSetup &setup, const PARAM_RES &res, const RES_PRODUCTION_CHANNEL &channel, const std::vector<MParticle> &legs, const HELMatrix *automatic = nullptr);

// Prepare one ordered pair of finite-spin continuum vertices
ReggeContinuumPole PreparePolePair(MProcessSetup &setup, ReggeProductionModel model, const MParticle &upper, const MParticle &lower, const MParticle &first, const MParticle &second);

// Prepare the four ordered finite-spin continuum vertices
std::vector<spin::PoleResidue> PreparePoleContinuum(MProcessSetup &setup, const std::vector<MDecayBranch> &tree, ReggeProductionModel model);

}  // namespace gra::rspin

#endif
