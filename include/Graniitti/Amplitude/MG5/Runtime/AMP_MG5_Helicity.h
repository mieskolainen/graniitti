// Shared MG5 helicity and color amplitudes
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef AMP_MG5_HELICITY_H
#define AMP_MG5_HELICITY_H

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Process/MProcessState.h"
#include "Graniitti/Spin/MHelicityBasis.h"

namespace gra {
namespace mg5helas {

// Classify one event-local generated matrix-element evaluation
enum class EvaluationStatus { Success, KinematicsFailure, AmplitudeFailure };

// Compute true only for a successful generated matrix-element evaluation
inline bool EvaluationSucceeded(EvaluationStatus status) { return status == EvaluationStatus::Success; }

// Compute true when the event uses the on-shell QED coupling
inline bool AlphaQEDAtZero(const LORENTZSCALAR &lts) {
  return lts.model_cache == nullptr || lts.model_cache->Tune().Structure().QED_alpha == "ZERO";
}

// Own the status and value from one scalar matrix-element evaluation
struct MatrixElementEvaluation {
  EvaluationStatus status = EvaluationStatus::Success;
  double           amp2   = 0.0;

  // Compute true when event kinematics and amplitudes were evaluated
  bool Valid() const { return EvaluationSucceeded(status); }
};


// Collect stable final-state PDG ids from one decay branch
inline void CollectStableFinalPDGs(const MDecayBranch &branch, std::vector<int> &pdgs) {
  if (branch.legs.empty()) {
    pdgs.push_back(branch.p.pdg);
    return;
  }
  for (const auto &leg : branch.legs) { CollectStableFinalPDGs(leg, pdgs); }
}

// Compute the factorial symmetry factor for stable final-state PDG ids
inline double FinalStateSymmetryFactor(const std::vector<int> &final_pdgs) {
  std::vector<int> sorted_pdgs = final_pdgs;
  std::sort(sorted_pdgs.begin(), sorted_pdgs.end());

  double symmetry = 1.0;
  for (std::size_t begin = 0; begin < sorted_pdgs.size();) {
    std::size_t end = begin + 1;
    while (end < sorted_pdgs.size() && sorted_pdgs[end] == sorted_pdgs[begin]) { ++end; }
    for (std::size_t factor = 2; factor <= end - begin; ++factor) { symmetry *= static_cast<double>(factor); }
    begin = end;
  }
  return symmetry;
}

// Compute the factorial symmetry factor for one stable decay tree
inline double FinalStateSymmetryFactor(const std::vector<MDecayBranch> &decaytree) {
  std::vector<int> final_pdgs;
  for (const auto &branch : decaytree) { CollectStableFinalPDGs(branch, final_pdgs); }
  return FinalStateSymmetryFactor(final_pdgs);
}

// Compute the symmetry factor applied by the outer GRANIITTI phase space
inline double AppliedFinalStateSymmetryFactor(const LORENTZSCALAR &lts) {
  if (lts.process.root_decay_mode == RootDecayMode::Isolated) { return 1.0; }
  return FinalStateSymmetryFactor(lts.decaytree);
}

// Compute true for a massless vector parton supported by the MG5 amplitudes
inline bool IsMasslessVector(int pdg) { return pdg == 21 || pdg == 22; }

// Transport one incoming massless HELAS fermion into the fixed beam-axis
// section
inline std::complex<double> IncomingFermionTransport(const double *momentum, int helicity) {
  const double pt2 = momentum[1] * momentum[1] + momentum[2] * momentum[2];
  if (momentum[3] >= 0.0 || !(pt2 > 0.0)) { return 1.0; }

  /*
   * HELAS uses a north-pole massless-spinor chart. Near the negative beam axis
   * its chi[1] component contains -exp(i h phi) relative to the exactly
   * collinear real spinor. Applying the inverse keeps screening-loop momenta in
   * one fixed little-group section.
   */
  const double phi = std::atan2(momentum[2], momentum[1]);
  return -std::polar(1.0, -static_cast<double>(helicity) * phi);
}

// Transport one incoming massless HELAS vector into the fixed transverse
// section
inline std::complex<double> IncomingVectorTransport(const double *momentum, int helicity) {
  const double pt2 = momentum[1] * momentum[1] + momentum[2] * momentum[2];
  if (pt2 <= 1.0e-30) { return 1.0; }

  /*
   * For the upper and lower incoming legs HELAS contributes exp(-i h phi) and
   * exp(+i h phi), respectively. The lower chart also differs by a common
   * minus sign from the exactly collinear south-pole branch.
   */
  const double phi = std::atan2(momentum[2], momentum[1]);
  if (momentum[3] >= 0.0) { return std::polar(1.0, static_cast<double>(helicity) * phi); }
  return -std::polar(1.0, -static_cast<double>(helicity) * phi);
}

// Transport one incoming MG5 external state according to its physical PDG id
inline std::complex<double> IncomingTransport(const double *momentum, int pdg, int helicity) {
  return IsMasslessVector(pdg) ? IncomingVectorTransport(momentum, helicity)
                               : IncomingFermionTransport(momentum, helicity);
}

// Transport both incoming MG5 legs into one fixed screening basis
inline std::complex<double> IncomingPairTransport(const double *momentum1, int pdg1, int helicity1,
                                                  const double *momentum2, int pdg2, int helicity2) {
  return IncomingTransport(momentum1, pdg1, helicity1) * IncomingTransport(momentum2, pdg2, helicity2);
}

// Compute orthogonal complex color components whose norm is the MG5 color sum
// Compute an empty vector when the generated color metric is invalid
// [REFERENCE: Lifson and Mattelaer, Eur. Phys. J. C 82, 1144 (2022), arXiv:2210.07267]
inline std::vector<std::complex<double>> ColorMetricAmplitudes(const std::vector<std::complex<double>> &jamp,
                                                               const std::vector<double>               &denominators,
                                                               const std::vector<std::vector<double>>  &color_factors) {
  const std::size_t n = jamp.size();
  if (n == 0 || denominators.size() != n || color_factors.size() != n) { return {}; }
  if (!gra::AllFinite(jamp)) { return {}; }

  MMatrix<double> metric(n, n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    if (!std::isfinite(denominators[i]) || !(denominators[i] > 0.0) || color_factors[i].size() != n) { return {}; }
    for (std::size_t j = 0; j < n; ++j) {
      if (!std::isfinite(color_factors[i][j])) { return {}; }
      metric[i][j] = color_factors[i][j] / denominators[i];
    }
  }

  // Cholesky factorize the symmetric positive semidefinite MG5 color metric
  // without concealing inconsistent row denominators by symmetrization
  MMatrix<double> lower;
  try {
    lower = metric.CholeskyLower(1.0e-12, CholeskyMode::AllowSemidefinite);
  } catch (const std::invalid_argument &) { return {}; } catch (const std::domain_error &) {
    return {};
  } catch (const std::runtime_error &) { return {}; }

  const std::vector<std::complex<double>> out = lower.LeftMultiply(jamp);
  if (!gra::AllFinite(out)) { return {}; }
  return out;
}

// Compute the incoherent norm of structured complex helicity components
inline double ComponentNorm(const std::vector<HelicityComponent> &components) {
  double norm = 0.0;
  for (const auto &component : components) { norm += std::norm(component.value); }
  return norm;
}

// Compute true when a physical beam pair defines the canonical EPA hard frame
inline bool HasEPAHardBeamState(const LORENTZSCALAR &lts) {
  const auto finite = [](const M4Vec &momentum) {
    return std::isfinite(momentum.E()) && std::isfinite(momentum.Px()) && std::isfinite(momentum.Py()) &&
           std::isfinite(momentum.Pz());
  };
  return finite(lts.pbeam1) && finite(lts.pbeam2) && lts.pbeam1.E() > 0.0 && lts.pbeam2.E() > 0.0;
}

// Prepare covariant massless incoming momenta while preserving the final state
bool PrepareOnShellKinematics(const LORENTZSCALAR &lts, const std::vector<M4Vec> &final, M4Vec &p1, M4Vec &p2);

// Prepare one reusable EPA hard rest frame and transform all final momenta
bool PrepareEPAHardFrame(const LORENTZSCALAR &lts, std::vector<M4Vec> &final, EPAHardFrame &frame);

// Remove the components of one EPA source in the physical two-beam plane
bool BeamTransverseSource(const M4Vec &q, const std::array<M4Vec, 2> &beam, M4Vec &source);

// Transform one node-local beam-transverse source pair into a prepared hard frame
bool TransformEPAHardSources(const EPAHardFrame &frame, const M4Vec &q1, const M4Vec &q2, std::array<M4Vec, 2> &source);

// Compute whether one fixed hard point uses exact elastic proton currents
bool UsesProtonEPAHardSources(const LORENTZSCALAR &lts);

// Remove longitudinal components along a physical on-shell incoming pair
M4Vec TransverseSourceVector(const M4Vec &source, const M4Vec &k1, const M4Vec &k2);

// Compute the HELAS massless incoming-vector polarization in (t,x,y,z) order
std::array<std::complex<double>, 4> IncomingVectorPolarization(const M4Vec &momentum, int helicity);

// Project one real source vector onto the physical HELAS helicities of an on-shell pair
std::array<std::complex<double>, 2> TransverseHelicityCoefficients(const M4Vec &source, const M4Vec &k1,
                                                                   const M4Vec &k2, bool normalize);

// Store the linear helicity projection of a transverse Cartesian source
struct TransverseHelicityProjector {
  bool                                               valid       = false;
  std::array<std::array<std::complex<double>, 2>, 2> coefficient = {};
  double                                             norm_xx     = 0.0;
  double                                             norm_xy     = 0.0;
  double                                             norm_yy     = 0.0;

  // Project qx e_x + qy e_y without rebuilding the covariant basis
  std::array<std::complex<double>, 2> Project(double qx, double qy) const;

  // Test the same finite transverse-source threshold used by Project
  bool Accepts(double qx, double qy) const;
};

// Prepare the covariant transverse-source map for one on-shell incoming leg
TransverseHelicityProjector PrepareTransverseHelicityProjector(const M4Vec &k1, const M4Vec &k2);

// Contract two prepared source-helicity coefficient vectors with one hard
// amplitude
std::complex<double> ContractHelicitySources(const std::array<std::complex<double>, 4> &hard,
                                             const std::array<std::complex<double>, 2> &source1,
                                             const std::array<std::complex<double>, 2> &source2);

// Contract one integrated incoming-helicity kernel with one hard amplitude
std::complex<double> ContractHelicityKernel(const std::array<std::complex<double>, 4> &hard,
                                            const std::array<std::complex<double>, 4> &kernel);

// Contract a Durham transverse-source pair with hard amplitudes ordered as
// (--,-+,+-,++)
std::complex<double> ContractTransverseSources(const std::array<std::complex<double>, 4> &hard, const M4Vec &source1,
                                               const M4Vec &source2, const M4Vec &k1, const M4Vec &k2);

// Contract fixed-basis gamma-gamma hard amplitudes with resolved EPA sources
std::vector<std::complex<double>> ContractEPAPhotonSources(LORENTZSCALAR                        &lts,
                                                           const std::vector<HelicityComponent> &components,
                                                           const M4Vec &k1, const M4Vec &k2);

// Contract a transfer-independent hard tensor with node-local transverse sources
std::vector<std::complex<double>> ContractEPAHardPhotonSources(LORENTZSCALAR                        &lts,
                                                               const std::vector<HelicityComponent> &components,
                                                               const M4Vec &source1, const M4Vec &source2,
                                                               const M4Vec &k1, const M4Vec &k2);

// Contract one fixed hard tensor with frame-covariant node-local EPA currents
std::vector<std::complex<double>> ContractEPAHardPhotonSources(LORENTZSCALAR                        &lts,
                                                               const std::vector<HelicityComponent> &components,
                                                               const EPAHardFrame                   &frame);

// Contract the Born event hard tensor with node-local photon currents
MatrixElementEvaluation ContractEPAHard(LORENTZSCALAR &lts, const EPAHardTensor &hard);

}  // namespace mg5helas
}  // namespace gra

#endif
