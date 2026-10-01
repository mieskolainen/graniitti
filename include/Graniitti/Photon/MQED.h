// QED couplings, currents and photon sources
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MQED_H
#define MQED_H

#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Own
#include "Graniitti/Kinematics/MEventKinematics.h"
#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Process/MProcessState.h"
#include "Graniitti/Spin/MDirac.h"

namespace gra {
namespace qed {

// Store one resolved elastic spin-half QED emitter
struct ElasticEmitter {
  gra::MDirac::FermionKind fermion_kind = gra::MDirac::FermionKind::Particle;
  double mass = 0.0;
  double charge = 0.0;
  double F1 = 0.0;
  double F2 = 0.0;
  double electric_form_factor = 0.0;
  double magnetic_form_factor = 0.0;
  std::size_t incoming_spin_states = 0;
};

// Store neutral-current couplings for one fermion species
struct NeutralCurrentCouplings {
  double chargeX1 = 0.0;
  double isospin3 = 0.0;
  double gL = 0.0;
  double gR = 0.0;
  double gV = 0.0;
  double gA = 0.0;
};

// Store one row-major photon matrix in the common (--,-+,+-,++) pair order
using ElasticSpinHalfPhotonMatrix = std::array<std::complex<double>, 16>;

constexpr double alpha_ref =
    0.0072973525693;                      // reference value of alpha at scale Q
constexpr double Q_ref = 0.0005109989461; // reference scale (GeV)
constexpr double alpha_0 =
    0.0072973525693; // value of alpha at Q = 0 (fine structure constant)^{-1}

// -----------------------------------------------------------------------

// Get number of charged leptons at given scale
inline int get_N_leptons(double Q) {
  if (Q > 1.77686E+00) {
    return 3;
  }
  if (Q > 1.056583745E-01) {
    return 2;
  }
  return 1;
}

// Compute the coefficient of the one-loop QED beta function
// beta_0^QED = -N_l/(3 pi) for d(1/alpha)/dln(Q2)
inline double beta0_QED(int N_lf) { return -1.0 / (3.0 * math::PI) * N_lf; }

// Compute the QED running coupling
// alpha(Q2)^-1 = alpha_ref^-1-sum_(m_l<Q) ln(Q2/m_l^2)/(3 pi)
inline double alpha_QED(double Q2, const std::string &order = "LL") {
  if (!std::isfinite(Q2) || Q2 < 0.0) {
    throw std::invalid_argument(
        "alpha_QED: Q2 must be finite and non-negative");
  }
  if (order != "ZERO" && order != "LL") {
    throw std::invalid_argument("alpha_QED: Unknown order parameter: " + order);
  }

  // Freeze the coupling at Q < Q_ref
  if (order == "ZERO" || Q2 < (Q_ref * Q_ref)) {
    return alpha_0;
  }

  // Match each charged-lepton vacuum-polarization logarithm at its threshold
  const double Q = math::msqrt(Q2);
  constexpr std::array<double, 3> mass = {Q_ref, 0.1056583745, 1.77686};
  double inverse_alpha = 1.0 / alpha_ref;
  for (const double threshold : mass) {
    if (Q > threshold) {
      inverse_alpha -=
          std::log(Q2 / math::pow2(threshold)) / (3.0 * math::PI);
    }
  }
  if (!(inverse_alpha > 0.0) || !std::isfinite(inverse_alpha)) {
    throw std::invalid_argument("alpha_QED: scale exceeds the LL domain");
  }
  return 1.0 / inverse_alpha;
}

// Compute alpha_QED according to the shared generator steering convention
inline double alpha_QED(const double Q2, const std::string &scheme,
                        const double alpha_mg) {
  if (!std::isfinite(alpha_mg) || !(alpha_mg > 0.0)) {
    throw std::invalid_argument(
        "alpha_QED: fixed model coupling must be finite and positive");
  }
  if (scheme == "ZERO") {
    return alpha_QED(0.0, "ZERO");
  }
  if (scheme == "MG") {
    return alpha_mg;
  }
  if (scheme == "LL") {
    return alpha_QED(Q2, "LL");
  }
  throw std::invalid_argument("alpha_QED: Unknown steering scheme: " + scheme);
}

// QED coupling at Q = 0
inline double alpha_QED() { return alpha_QED(0, "ZERO"); }

// Compute the electric charge in natural units
// e(Q2) = sqrt[4 pi alpha(Q2)]
inline double e_QED(double Q2) {
  return math::msqrt(alpha_QED(Q2) * 4.0 * math::PI);
}
inline double e_QED() { return math::msqrt(alpha_QED() * 4.0 * math::PI); }

// Compute a lower-index fermion-antifermion vector current
gra::MDirac::Current FFVCurrent(const gra::MDirac &dirac,
                                const gra::M4Vec &fermion,
                                const gra::M4Vec &antifermion,
                                int helicity_fermion, int helicity_antifermion,
                                std::complex<double> cL,
                                std::complex<double> cR);

// Compute a lower-index photon to fermion-antifermion current
gra::MDirac::Current PhotonFermionCurrent(const gra::MDirac &dirac,
                                          const gra::M4Vec &fermion,
                                          const gra::M4Vec &antifermion,
                                          int fermion_pdg, int helicity_fermion,
                                          int helicity_antifermion,
                                          double electric_charge);

// Compute neutral-current couplings for charged leptons and quarks
NeutralCurrentCouplings FermionNeutralCurrentCouplings(int pdg,
                                                       double sin2thetaW);

// Compute a lower-index Z to fermion-antifermion current
gra::MDirac::Current ZFermionCurrent(const gra::MDirac &dirac,
                                     const gra::M4Vec &fermion,
                                     const gra::M4Vec &antifermion,
                                     int fermion_pdg, int helicity_fermion,
                                     int helicity_antifermion,
                                     double electric_charge, double sin2thetaW);

// Validate the photon source steering mode
void ValidatePhotonMode(const std::string &mode, const std::string &context);

// Resolve one supported elastic spin-half QED emitter
ElasticEmitter Emitter(const gra::MParticle &particle, double t,
                       const std::string &context,
                       const gra::form::ParamStore &structure = {});

// Resolve one supported elastic spin-half QED emitter from a beam leg
ElasticEmitter Emitter(const gra::LORENTZSCALAR &lts, int leg, double t,
                       const std::string &context);

// Compute the exact elastic one-photon pair-helicity matrix
ElasticSpinHalfPhotonMatrix ElasticSpinHalfPhotonExchange(
    const gra::MDirac &dirac, const gra::MParticle &beam1,
    const gra::MParticle &beam2, const gra::M4Vec &p1_in,
    const gra::M4Vec &p2_in, const gra::M4Vec &p1_out,
    const gra::M4Vec &p2_out,
    const gra::form::ParamStore &structure = {});

// Validate one exclusive beam leg for the elastic QED current
void ValidateEmitter(const gra::LORENTZSCALAR &lts, int leg,
                     const std::string &context);

// Compute the electric charge of one photon-emitting beam leg
double Charge(const gra::LORENTZSCALAR &lts, int leg);

// Compute the elastic EPA density used to match the QED current
double ElasticPhotonDensity(const gra::LORENTZSCALAR &lts, int leg);

// Build unnormalized elastic photon currents with propagator and collider phases
std::vector<MDirac::Current> PhotonCurrentTransitions(
    const gra::LORENTZSCALAR &lts, int leg,
    const std::vector<std::pair<double, double>> &transitions);

// Compute the EPA or matched QED photon-source amplitude normalization
std::complex<double> PhotonSourceAmplitude(const gra::LORENTZSCALAR &lts,
                                           const gra::ForwardLegState &state,
                                           const std::string &mode);

// Compute the transverse EPA source amplitude for one global spin projection
std::complex<double>
TransversePhotonSourceAmplitude(const LORENTZSCALAR &lts,
                                const ForwardLegState &state,
                                int boson_spin_projection);

// Build a matched EPA or QED photon source for explicit helicity transitions
MMatrix<std::complex<double>> PhotonSourceMatrixTransitions(
    const gra::LORENTZSCALAR &lts, int leg,
    const std::vector<std::pair<double, double>> &transitions,
    const std::vector<int> &m_values, const std::string &mode,
    bool second_exchange_daughter = false);

} // namespace qed
} // namespace gra

#endif
