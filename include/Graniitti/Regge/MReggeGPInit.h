// GP production vertex preparation
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MREGGEGPINIT_H
#define MREGGEGPINIT_H

#include "Graniitti/Process/MProcessSetup.h"
#include "Graniitti/Regge/MReggeGP.h"

namespace gra::gpom {

// Initialize a GP analytic trajectory H(lambda1,lambda2;m) matrix
void InitHelicity(HELMatrix& hc, double s1, double s2, int MMAX, const std::string& context);

// Initialize one angle-free crossed Regge-helicity vertex
void InitCrossed(HELMatrix& hc, double s1, double s2, int MMAX, const std::string& context);

// Initialize one crossed photon vertex in its physical spin one basis
void InitPhotonCrossed(HELMatrix& hc, double s1, double s2, const std::string& context);

// Validate one GP analytic trajectory helicity matrix without rescaling it
void CheckHelicity(HELMatrix& hc, const std::string& context, bool verbose);

// Prepare one analytic resonance fusion vertex
HELMatrix PrepareResonance(const MParticle &mother, const std::vector<MParticle> &legs, const RES_PRODUCTION_CHANNEL &channel, const regge::Param &param, const MPDG &pdg_table, int MMAX, double coupling_min);
// Normalize the analytic diphoton operator at its physical pole
void NormalizeGammaGamma(HELMatrix &hel, const PARAM_RES &res, const std::vector<MParticle> &legs, bool derivative_factor, bool isolated_decay);
// Prepare four ordered analytic continuum vertices
std::vector<HELMatrix> PrepareContinuum(MProcessSetup &setup, const std::vector<MDecayBranch> &tree, const regge::Param &param);
// Prepare one local analytic continuum pair
ReggeContinuumPole PreparePair(MProcessSetup &setup, const regge::Param &param, const MParticle &upper, const MParticle &lower, const MParticle &first, const MParticle &second, const std::string &context);
// Restrict a ladder vertex to its supported m=0 transport
void PrepareLadder(HELMatrix &hel, const std::string &context);
// Print the analytic continuum input basis
void PrintContinuum(const LORENTZSCALAR &lts, const std::vector<int> &exchange, const std::vector<MDecayBranch> &tree);

}  // namespace gra::gpom

#endif
