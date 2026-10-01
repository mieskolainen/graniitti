// Durham QCD amplitudes for central production
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <mutex>
#include <random>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

// Own
#include "Graniitti/Amplitude/MG5/Runtime/AMP_MG5_Helicity.h"
#include "Graniitti/Amplitude/Photon/AMP_gg_yy.h"
#include "Graniitti/Amplitude/MG5/Runtime/AMP_MG5_Kinematics.h"
#include "Graniitti/Kinematics/MEventKinematics.h"
#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/Math/MIntegration.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Process/MDecayTree.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/PDF/MSudakov.h"
#include "Graniitti/Particle/MForm.h"
#include "Graniitti/Particle/MPDG.h"
#include "Graniitti/QCD/MDurham.h"
#include "Graniitti/Spin/MDirac.h"
#include "Graniitti/Spin/MHelicityBasis.h"
#include "Graniitti/Spin/MHelicityDecay.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MException.h"

// Libraries
#include "json.hpp"
#include "rang.hpp"

using gra::aux::indices;

namespace gra {
using math::msqrt;
using math::PI;
using math::pow2;
using math::zi;
using PDG::mp;

// Central system J^P = 0^+, 0^-, +2^+, -2^+
enum SPINPARITY { P0, M0, P2, M2 };  // Implicit conversion to int

namespace {

// Compute one Durham or generated MG5 pair index in (--,-+,+-,++) order
constexpr std::size_t DurhamHelicityIndex(const int first, const int second) {
  return spin::BinaryPairHelicityIndexX2(first, second);
}

struct DurhamCharmoniumState {
  double mass  = 0.0;
  double width = 0.0;
};

// Compute the PDG code associated with one Durham charmonium subprocess
//
int DurhamCharmoniumPDG(const std::string &process) {
  if (process == "chic(0)") { return 10441; }
  if (process == "chic(1)") { return 20443; }
  if (process == "chic(2)") { return 445; }
  throw std::invalid_argument("DurhamCharmonium: unknown process " + process);
}

// Compute the active pole parameters for one Durham charmonium subprocess
//
DurhamCharmoniumState DurhamCharmonium(const gra::LORENTZSCALAR &lts, const std::string &process) {
  const int             pdg      = DurhamCharmoniumPDG(process);
  const gra::MParticle *particle = nullptr;

  if (lts.process.ROOT_RES_ACTIVE) {
    if (lts.process.ROOT_RES.p.pdg != pdg) {
      throw std::invalid_argument("DurhamCharmonium: active root resonance does not match subprocess " + process);
    }
    particle = &lts.process.ROOT_RES.p;
  } else {
    particle = &lts.PDG.FindByPDG(pdg);
  }

  if (!std::isfinite(particle->mass) || !std::isfinite(particle->width) || particle->mass <= 0.0 ||
      particle->width <= 0.0) {
    throw std::invalid_argument("DurhamCharmonium: invalid pole parameters for subprocess " + process);
  }
  return {particle->mass, particle->width};
}

// Compute the chi_c0 pole parameters used to normalize every chi_cJ hard vertex
//
DurhamCharmoniumState DurhamCharmoniumNormalizationState(const gra::LORENTZSCALAR &lts) {
  const gra::MParticle *particle = lts.process.ROOT_RES_ACTIVE && lts.process.ROOT_RES.p.pdg == 10441
                                       ? &lts.process.ROOT_RES.p
                                       : &lts.PDG.FindByPDG(10441);

  if (!std::isfinite(particle->mass) || !std::isfinite(particle->width) || particle->mass <= 0.0 ||
      particle->width <= 0.0) {
    throw std::invalid_argument("DurhamCharmoniumNormalizationState: invalid chi_c0 pole parameters");
  }
  return {particle->mass, particle->width};
}

// Extract the common alpha_s |phi'_P(0)| product from the chi_c0 two-gluon
// width Gamma_gg = 96 alpha_s^2 |phi'_P(0)|^2 / M_chic0^4 SuperChic uses
// Gamma_gg = Gamma_tot / K_NLO with K_NLO = 1.5
// [REFERENCE: Harland-Lang et al., arxiv.org/abs/1508.02718v2, Eq. (35)]
//
double DurhamCharmoniumAlphaPhiPrime(const DurhamCharmoniumState &normalization_state, double gg_width_fraction) {
  const double gamma_gg = gg_width_fraction * normalization_state.width;
  return pow2(normalization_state.mass) * msqrt(gamma_gg / 96.0);
}

// Compute the root decay coupling matched to Durham's normalized delta-BW kernel
// coupling = g_decay sqrt[pi/(M Gamma)]
std::complex<double> DurhamRootDecayCoupling(const gra::PARAM_RES &res) {
  if (!std::isfinite(res.p.mass) || !std::isfinite(res.p.width) || res.p.mass <= 0.0 || res.p.width <= 0.0) {
    throw std::invalid_argument(
        "DurhamRootDecayCoupling: Durham resonance "
        "decay needs a finite mass and width");
  }
  return res.hel_decay.g_decay * std::sqrt(PI / (res.p.mass * res.p.width));
}

struct DurhamEtaMixing {
  double octet   = 0.0;
  double singlet = 0.0;
};

// Compute the number of orthogonal amplitude channels for one Durham subprocess
std::size_t DurhamChannelCount(const std::string &process) {
  if (process == "gg") { return 4; }
  if (process == "chic(1)") { return 3; }
  if (process == "chic(2)") { return 5; }
  if (process == "MMbar" || process == "chic(0)" || process == "FLUX") { return 1; }
  throw std::invalid_argument("DurhamChannelCount: unknown subprocess " + process);
}

// Compute true for Durham subprocesses whose central state carries QCD color
bool IsDurhamColoredProcess(const std::string &process) { return process == "gg"; }

// Compute true when one four-vector has finite physical components
//
bool DurhamFourVectorFinite(const gra::M4Vec &p) {
  return std::isfinite(p.E()) && std::isfinite(p.Px()) && std::isfinite(p.Py()) && std::isfinite(p.Pz());
}

// Zero color amplitudes while preserving the screening channel layout
void ZeroDurhamColorAmplitudes(gra::LORENTZSCALAR &lts, std::size_t channels) {
  for (auto &flow : lts.hard_color_flows) {
    flow.amplitudes.assign(channels, 0.0);
    flow.screened_weight.reset();
  }
}

// Compute true when one event lies inside the perturbative and tabulated Durham
// domain
bool DurhamEventDomainValid(const gra::LORENTZSCALAR &lts, const MDurhamParam &param) {
  if (lts.GlobalSudakovPtr == nullptr || !std::isfinite(lts.s_hat) || lts.s_hat <= 0.0) { return false; }

  if (lts.pfinal.size() < 3 || !DurhamFourVectorFinite(lts.pbeam1) || !DurhamFourVectorFinite(lts.pbeam2)) {
    return false;
  }
  for (std::size_t i = 0; i < 3; ++i) {
    if (!DurhamFourVectorFinite(lts.pfinal[i])) { return false; }
  }
  if (!(lts.pfinal[0].E() > 0.0) || !(lts.pfinal[0].M2() > 0.0)) { return false; }

  const MDurhamScales scales = DurhamCentralScales(std::sqrt(lts.s_hat), param);
  return std::isfinite(scales.muR) && scales.muR > 0.0 && lts.GlobalSudakovPtr->SupportsMu(scales.muF);
}

// Compute the hard coupling and propagate failures to the sampling boundary
double DurhamEventAlphaS(const gra::LORENTZSCALAR& lts, const MDurhamParam& param) {
  try {
    const auto   scales  = DurhamCentralScales(std::sqrt(lts.s_hat), param);
    const double alpha_s = lts.GlobalSudakovPtr->AlphaS_Q2(scales.muR * scales.muR);
    if (!std::isfinite(alpha_s) || !(alpha_s > 0.0)) { throw AmplitudeFailure("Durham: invalid alpha_s"); }
    return alpha_s;
  } catch (const std::logic_error& error) { throw AmplitudeFailure(error.what()); }
}

// Reuse the central hard-scale coupling already evaluated by DurhamQCD
double DurhamHardAlphaS(const gra::LORENTZSCALAR& lts, const MDurhamParam& param) {
  const MDurhamScales scales          = DurhamCentralScales(std::sqrt(lts.s_hat), param);
  const double        scale_tolerance = 1.0e-12 * std::max(1.0, std::abs(scales.muR));
  if (std::isfinite(lts.alphaQCD) && lts.alphaQCD > 0.0 && std::isfinite(lts.muR) &&
      std::abs(lts.muR - scales.muR) <= scale_tolerance) {
    return lts.alphaQCD;
  }
  return DurhamEventAlphaS(lts, param);
}

// Compute the eta/eta-prime two-angle octet-singlet decay-constant coefficients
//
// [REFERENCE: Feldmann, Kroll and Stech, arXiv:hep-ph/9802409]
//
DurhamEtaMixing DurhamEtaMixingCoefficients(int pdg, double eta_theta8_deg, double eta_theta1_deg) {
  const double theta8 = math::Deg2Rad(eta_theta8_deg);
  const double theta1 = math::Deg2Rad(eta_theta1_deg);

  const int apdg = std::abs(pdg);
  if (apdg == 221) { return {std::cos(theta8), -std::sin(theta1)}; }
  if (apdg == 331) { return {std::sin(theta8), std::cos(theta1)}; }
  return {};
}

// Compute the sign of the four-dimensional Levi-Civita tensor epsilon^{0123}
//
int LeviCivita4(std::size_t a, std::size_t b, std::size_t c, std::size_t d) {
  const std::array<std::size_t, 4> index = {a, b, c, d};
  for (const auto &i : indices(index)) {
    for (std::size_t j = i + 1; j < index.size(); ++j) {
      if (index[i] == index[j]) { return 0; }
    }
  }

  int sign = 1;
  for (const auto &i : indices(index)) {
    for (std::size_t j = i + 1; j < index.size(); ++j) {
      if (index[i] > index[j]) { sign *= -1; }
    }
  }
  return sign;
}

// Compute the squared norm of one transverse Durham momentum
//
double DurhamTransverseMomentum2(const MDurham::DurhamTransverseMomentum &q) { return pow2(q[0]) + pow2(q[1]); }

// Compute true when both transverse loop momenta are finite
//
bool DurhamLoopMomentaValid(const MDurham::DurhamTransverseMomentum &q1, const MDurham::DurhamTransverseMomentum &q2) {
  return std::isfinite(q1[0]) && std::isfinite(q1[1]) && std::isfinite(q2[0]) && std::isfinite(q2[1]);
}

// Minkowski product of the Durham transverse embeddings
//
//   q_i^mu = (0, q_ix, q_iy, 0)
//
// in the pp working frame.  This is a Lorentz scalar and therefore may be used
// directly in the covariant chi_cJ vertices
//
double DurhamTransverseMinkowskiDot(const MDurham::DurhamTransverseMomentum &q1,
                                    const MDurham::DurhamTransverseMomentum &q2) {
  return -(q1[0] * q2[0] + q1[1] * q2[1]);
}

// Rotate one vector so the chosen axis defines the local hard-scattering z-axis
//
gra::M4Vec DurhamRotateToAxisFrame(const gra::M4Vec &p, const gra::M4Vec &axis) {
  gra::M4Vec out = p;
  out.RotateZ(-axis.Phi());
  out.RotateY(-axis.Theta());
  return out;
}

// Embed parity-even meson helicities in the physical incoming HELAS basis
// V = (A0+A2) g_perp + 2 A2 d_perp d_perp / (-d_perp^2)
// In the planar hard frame this gives M_++ = A0 and M_+- = A2
MDurham::DurhamInitialHelicity DurhamMesonTensor(const M4Vec &k1, const M4Vec &k2, const M4Vec &meson, double A0,
                                                 double A2) {
  MDurham::DurhamInitialHelicity out{};
  const M4Vec                    transverse = mg5helas::TransverseSourceVector(meson, k1, k2);
  const double                   pt2        = -transverse.M2();
  const double                   kdot       = k1 * k2;
  if (!std::isfinite(pt2) || !(pt2 > 0.0) || !std::isfinite(kdot) || !(kdot > 0.0)) { return out; }
  const auto first     = k1.Contravariant<double>();
  const auto second    = k2.Contravariant<double>();
  const auto direction = transverse.Contravariant<double>();
  for (const int h1 : {-1, 1}) {
    const auto eps1 = mg5helas::IncomingVectorPolarization(k1, h1);
    for (const int h2 : {-1, 1}) {
      const auto eps2 = mg5helas::IncomingVectorPolarization(k2, h2);
      const auto metric =
          MinkowskiProduct(eps1, eps2) - (MinkowskiProduct(eps1, first) * MinkowskiProduct(eps2, second) +
                                          MinkowskiProduct(eps1, second) * MinkowskiProduct(eps2, first)) /
                                             kdot;
      out[DurhamHelicityIndex(h1, h2)] =
          (A0 + A2) * metric + 2.0 * A2 * MinkowskiProduct(eps1, direction) * MinkowskiProduct(eps2, direction) / pt2;
    }
  }
  return out;
}

// Compute the lowered real Lorentz components of one four-vector
//
std::array<double, 4> DurhamLowerComponents(const gra::M4Vec &p) { return {p % 0, p % 1, p % 2, p % 3}; }

using DurhamSpin1Basis = std::array<FTensor::Tensor1<std::complex<double>, 4>, 3>;
using DurhamSpin2Basis = std::array<FTensor::Tensor2<std::complex<double>, 4, 4>, 5>;

// Only the x and y components of the axial projector are needed in the Durham
// loop, because q_i^mu=(0,q_ix,q_iy,0) in the pp working frame
using DurhamAxialProjector = std::array<std::array<std::complex<double>, 2>, 3>;

struct DurhamTensorProjector {
  // Only the transverse 2x2 block enters q1_mu q2_nu eps^{*mu nu}
  std::array<std::array<std::array<std::complex<double>, 2>, 2>, 5> transverse = {};

  // p1_mu p2_nu eps^{*mu nu}, event-local and loop independent
  std::array<std::complex<double>, 5> beam = {};
};

// Canonically (rotationlessly) boost one complex Lorentz vector from the
// central-system rest frame into the pp working frame
//
// The CM spin labels m are therefore preserved: this is the inverse of the
// old event-by-event boost of all momenta into the central rest frame, but is
// applied only once to the polarization basis
//
FTensor::Tensor1<std::complex<double>, 4> DurhamCanonicalBoostSpin1(
    const gra::M4Vec &system, const FTensor::Tensor1<std::complex<double>, 4> &rest) {
  const double mass = system.M();
  const double e    = system.E();
  if (!(mass > 0.0) || !(e > 0.0) || !DurhamFourVectorFinite(system)) {
    throw AmplitudeFailure(
        "DurhamCanonicalBoostSpin1: central system must be "
        "finite, timelike and positive-energy");
  }

  const std::array<double, 3> beta = {system.Px() / e, system.Py() / e, system.Pz() / e};

  const double beta2 = pow2(beta[0]) + pow2(beta[1]) + pow2(beta[2]);
  if (!std::isfinite(beta2) || beta2 >= 1.0) {
    throw AmplitudeFailure("DurhamCanonicalBoostSpin1: invalid central-system boost");
  }

  const double gamma = e / mass;

  // beta . epsilon_rest (ordinary Euclidean spatial dot product)
  const std::complex<double> beta_dot = beta[0] * rest(1) + beta[1] * rest(2) + beta[2] * rest(3);

  // (gamma - 1)/beta^2 = gamma^2/(gamma + 1), with a stable beta -> 0 limit
  const double boost_coeff = gamma * gamma / (gamma + 1.0);

  FTensor::Tensor1<std::complex<double>, 4> out;
  out(0) = gamma * (rest(0) + beta_dot);
  out(1) = rest(1) + boost_coeff * beta_dot * beta[0] + gamma * rest(0) * beta[0];
  out(2) = rest(2) + boost_coeff * beta_dot * beta[1] + gamma * rest(0) * beta[1];
  out(3) = rest(3) + boost_coeff * beta_dot * beta[2] + gamma * rest(0) * beta[2];
  return out;
}

// Compute the canonical spin-1 basis in the pp working frame
//
// Do not call EpsMassiveSpin1(system,m) here: that constructor produces the
// helicity basis quantized along the moving system momentum.  Durham resonance
// decays use the CM spin-projection basis, so we construct that basis at rest
// and canonically boost it once
//
DurhamSpin1Basis DurhamSpin1Polarizations(const gra::MDirac &dirac, const gra::M4Vec &system) {
  const double mass = system.M();
  if (!(mass > 0.0) || !std::isfinite(mass)) {
    throw AmplitudeFailure("DurhamSpin1Polarizations: central system must be timelike");
  }

  const gra::M4Vec rest_system(0.0, 0.0, 0.0, mass);

  DurhamSpin1Basis eps;
  for (int m = -1; m <= 1; ++m) {
    const auto rest_eps                  = dirac.EpsMassiveSpin1(rest_system, m);
    eps[static_cast<std::size_t>(m + 1)] = DurhamCanonicalBoostSpin1(system, rest_eps);
  }
  return eps;
}

// Build the spin-2 basis directly from the already boosted spin-1 basis
//
// Lorentz transformation is linear, so coupling the boosted spin-1 vectors is
// exactly equivalent to constructing epsilon^{mu nu} in the CM and boosting
// both Lorentz indices. The spin labels remain the original CM labels
//
DurhamSpin2Basis DurhamSpin2Polarizations(const DurhamSpin1Basis &eps1) {
  const auto &em = eps1[0];
  const auto &e0 = eps1[1];
  const auto &ep = eps1[2];

  constexpr double invsqrt2 = 0.707106781186547524400844362104849039;
  constexpr double invsqrt6 = 0.408248290463863016366214012450981899;

  DurhamSpin2Basis out;
  for (std::size_t mu = 0; mu < 4; ++mu) {
    for (std::size_t nu = 0; nu < 4; ++nu) {
      out[0](mu, nu) = em(mu) * em(nu);
      out[1](mu, nu) = (em(mu) * e0(nu) + e0(mu) * em(nu)) * invsqrt2;
      out[2](mu, nu) = (ep(mu) * em(nu) + em(mu) * ep(nu) + 2.0 * e0(mu) * e0(nu)) * invsqrt6;
      out[3](mu, nu) = (ep(mu) * e0(nu) + e0(mu) * ep(nu)) * invsqrt2;
      out[4](mu, nu) = ep(mu) * ep(nu);
    }
  }
  return out;
}

// Lower one complex Lorentz-vector component with the (+,-,-,-) metric
//
std::complex<double> LowerComplexVectorComponent(const FTensor::Tensor1<std::complex<double>, 4> &v, std::size_t mu) {
  return (mu == 0) ? v(mu) : -v(mu);
}

// Precontract the covariant chi_c1 beam/polarization pseudovector
//
// A_m^mu = eps^{mu nu alpha beta} p1_nu p2_alpha eps^*_{m,beta}
//
// Only mu=x,y is retained because the Durham loop vectors have only transverse
// components in the pp working frame
//
DurhamAxialProjector DurhamAxialBeamProjector(const gra::M4Vec &p1, const gra::M4Vec &p2, const DurhamSpin1Basis &eps) {
  DurhamAxialProjector projector = {};
  const auto           p1_lo     = DurhamLowerComponents(p1);
  const auto           p2_lo     = DurhamLowerComponents(p2);

  for (const auto &m : indices(projector)) {
    for (std::size_t transverse = 0; transverse < 2; ++transverse) {
      const std::size_t mu = transverse + 1;  // x,y Lorentz indices

      for (std::size_t nu = 0; nu < 4; ++nu) {
        for (std::size_t alpha = 0; alpha < 4; ++alpha) {
          for (std::size_t beta = 0; beta < 4; ++beta) {
            const int eps4 = LeviCivita4(mu, nu, alpha, beta);
            if (eps4 == 0) { continue; }

            projector[m][transverse] += static_cast<double>(eps4) * p1_lo[nu] * p2_lo[alpha] *
                                        std::conj(LowerComplexVectorComponent(eps[m], beta));
          }
        }
      }
    }
  }
  return projector;
}

// Precontract the covariant chi_c2 polarization tensors in the pp working
// frame
//
// The beam contraction needs all 4x4 components once per event.  The Q_t loop
// needs only the transverse 2x2 block
//
DurhamTensorProjector DurhamSpin2Projector(const gra::M4Vec &p1, const gra::M4Vec &p2, const DurhamSpin2Basis &eps) {
  DurhamTensorProjector projector;
  const auto            p1_lo = DurhamLowerComponents(p1);
  const auto            p2_lo = DurhamLowerComponents(p2);

  for (const auto &m : indices(projector.beam)) {
    for (std::size_t mu = 0; mu < 4; ++mu) {
      for (std::size_t nu = 0; nu < 4; ++nu) {
        const std::complex<double> spin = std::conj(eps[m](mu, nu));
        projector.beam[m] += p1_lo[mu] * p2_lo[nu] * spin;

        if (mu >= 1 && mu <= 2 && nu >= 1 && nu <= 2) { projector.transverse[m][mu - 1][nu - 1] = spin; }
      }
    }
  }
  return projector;
}

// Compute the full off-shell gluon scalar product
// q1.q2 = (M_X^2+q1^2+q2^2)/2
double DurhamGluonDot(const DurhamCharmoniumState &state, double q1_2, double q2_2) {
  return 0.5 * (pow2(state.mass) + q1_2 + q2_2);
}

// Compute the common chi_c off-shell vertex coupling c_chi
//
// The Breit-Wigner factor is the normalized complex delta-BW amplitude,
// with propagator 1/(s-M^2+iM Gamma) and negative imaginary pole phase
//
std::complex<double> DurhamCharmoniumCoupling(const std::complex<double>  &breit_wigner,
                                              const DurhamCharmoniumState &state, double q1_2, double q2_2,
                                              double alpha_phi_prime) {
  const double NC   = 3.0;
  const double q1q2 = DurhamGluonDot(state, q1_2, q2_2);

  std::complex<double> coupling =
      16.0 * PI / (2.0 * msqrt(NC) * pow2(q1q2)) * msqrt(6.0 / (4.0 * PI * state.mass)) * alpha_phi_prime;
  coupling *= breit_wigner;
  return coupling;
}

// Convert a reduced Durham Mbar kernel to the contracted q_i V_ij q_j kernel
//
std::complex<double> DurhamMbarToQV(double central_mass2, const std::complex<double> &mbar) {
  return 0.5 * central_mass2 * mbar;
}

// Compute true for Durham charmonium subprocesses with resonance decay support
//
bool IsDurhamCharmonium(const std::string &process) {
  return process == "chic(0)" || process == "chic(1)" || process == "chic(2)";
}

// Compute true when the physical arrow activates Durham root spin decays
//
bool DurhamDecaySpinActive(const gra::LORENTZSCALAR &lts, const std::string &process) {
  return lts.process.ROOT_RES_ACTIVE && lts.process.root_decay_mode == RootDecayMode::Physical &&
         IsDurhamCharmonium(process);
}

// Compute the Durham root decay matrix bound to the current Born event
//
const MMatrix<std::complex<double>> &DurhamDecayMatrixForEvent(gra::LORENTZSCALAR            &lts,
                                                               MMatrix<std::complex<double>> &local) {
  constexpr const char *kDurhamNativeFrame = "CM";
  constexpr const char *kDurhamDecayKey    = "root";

  if (!lts.process.ROOT_RES_ACTIVE) {
    throw std::invalid_argument(
        "MDurham::DurhamDecayMatrixForEvent: Durham "
        "resonance decay is not initialized");
  }
  if (lts.screening.active) {
    return lts.amplitude.durham_decay.Get(lts.amplitude.central, kDurhamDecayKey, "MDurham::DurhamDecayMatrixForEvent");
  }

  local = gra::spin::ResonanceDecayMatrix(lts, lts.process.ROOT_RES, kDurhamNativeFrame);
  if (lts.amplitude.central.Active()) {
    return lts.amplitude.durham_decay.Store(lts.amplitude.central, kDurhamDecayKey, std::move(local),
                                            "MDurham::DurhamDecayMatrixForEvent");
  }
  return local;
}

// Validate Durham production channels are in resonance spin-projection order
//
std::vector<std::complex<double>> DurhamSpinProjectionAmplitudes(const std::string                       &process,
                                                                 const std::vector<std::complex<double>> &production) {
  if (process == "chic(0)") {
    if (production.size() != 1) {
      throw AmplitudeFailure("DurhamSpinProjectionAmplitudes: chic(0) expects one channel");
    }
    return production;
  }

  if (process == "chic(1)") {
    if (production.size() != 3) {
      throw AmplitudeFailure("DurhamSpinProjectionAmplitudes: chic(1) expects three channels");
    }
    return production;
  }

  if (process == "chic(2)") {
    if (production.size() != 5) {
      throw AmplitudeFailure("DurhamSpinProjectionAmplitudes: chic(2) expects five channels");
    }
    return production;
  }

  throw std::invalid_argument("DurhamSpinProjectionAmplitudes: unknown process " + process);
}

// Contract Durham production amplitudes with MHelicity root decay machinery
//
double DurhamDecayAmp2(gra::LORENTZSCALAR &lts, const std::string &process,
                       const std::vector<std::complex<double>> &production) {
  MMatrix<std::complex<double>> local_decay;
  const auto                   &decay_f = DurhamDecayMatrixForEvent(lts, local_decay);

  const std::vector<std::complex<double>> spin_amp = DurhamSpinProjectionAmplitudes(process, production);
  if (spin_amp.size() != decay_f.size_row()) {
    throw AmplitudeFailure("MDurham::DurhamDecayAmp2: production/decay basis mismatch");
  }

  const std::complex<double> decay_coupling = DurhamRootDecayCoupling(lts.process.ROOT_RES);

  lts.hamp.assign(decay_f.size_col(), 0.0);
  for (std::size_t col = 0; col < decay_f.size_col(); ++col) {
    for (std::size_t row = 0; row < decay_f.size_row(); ++row) {
      lts.hamp[col] += spin_amp[row] * decay_f[row][col] * decay_coupling;
    }
  }

  return gra::SquaredNorm(lts.hamp);
}

// Compute the conversion from a normalized gluon singlet to the Durham colour average
// factor = 1/sqrt(Nc^2-1)
double DurhamIncomingColorAverageFactor() {
  constexpr double NC = 3.0;
  return 1.0 / std::sqrt(pow2(NC) - 1.0);
}

// Match generated incoming phases to the selected Durham helicity projector
MDurham::DurhamInitialHelicity DurhamHelicityPhases(const gra::LORENTZSCALAR &lts, bool jzp) {
  MDurham::DurhamInitialHelicity phase;
  phase.fill(1.0);
  if (!jzp) { return phase; }
  const auto                 &k1 = lts.screening.durham.hard_k1;
  const auto                 &k2 = lts.screening.durham.hard_k2;
  const std::array<double, 4> p1 = {k1.E(), k1.Px(), k1.Py(), k1.Pz()};
  const std::array<double, 4> p2 = {k2.E(), k2.Px(), k2.Py(), k2.Pz()};
  for (const int h1 : {-1, 1}) {
    for (const int h2 : {-1, 1}) {
      const auto transport =
          mg5helas::IncomingPairTransport(p1.data(), PDG::PDG_gluon, h1, p2.data(), PDG::PDG_gluon, h2);
      phase[DurhamHelicityIndex(h1, h2)] = transport;
    }
  }
  return phase;
}

// Convert one flattened MG5 projector row to Durham final-by-initial layout
MDurham::DurhamProjectedAmp BuildDurhamProjectedChannels(const std::complex<double> *projected,
                                                         std::size_t helicity_count, std::size_t final_helicity_count,
                                                         double                                normalization,
                                                         const MDurham::DurhamInitialHelicity &phase) {
  if (projected == nullptr || helicity_count != 4 * final_helicity_count) {
    throw AmplitudeFailure("BuildDurhamProjectedChannels: incompatible helicity dimensions");
  }

  MDurham::DurhamProjectedAmp out(final_helicity_count);
  for (auto &channel : out) { channel.fill(0.0); }
  for (std::size_t initial = 0; initial < 4; ++initial) {
    for (std::size_t final = 0; final < final_helicity_count; ++final) {
      out[final][initial] = normalization * phase[initial] * projected[final_helicity_count * initial + final];
    }
  }
  return out;
}

// Integrate one square matrix sampled on a tensor-product quadrature rule
double TensorProductIntegral(const MMatrix<double> &values, const std::vector<double> &weights) {
  if (values.size_row() != weights.size() || values.size_col() != weights.size()) {
    throw std::invalid_argument("TensorProductIntegral: matrix and quadrature dimensions differ");
  }
  return values.MatrixElement(weights, weights);
}

}  // namespace

// Construct one immutable Durham model and loop configuration
std::shared_ptr<const MDurham::DurhamConfig> ReadDurhamConfig(const MModelTune &tune) {
  auto config = std::make_shared<MDurham::DurhamConfig>();
  config->param.ConfigureFromJson(tune.GeneralFile(), tune.General().dump());
  config->numerics.ConfigureFromJson(tune.NumericsFile(), tune.Numerics().dump());
  const double cut = config->param.loop_q2_cut;
  if (!(cut < config->numerics.qt2_MAX)) {
    throw std::invalid_argument("ReadDurhamConfig: loop_q2_cut must be below qt2_MAX");
  }
  if (!(cut > 0.0) && config->param.PDF_scale != "MIN") {
    throw std::invalid_argument("ReadDurhamConfig: only MIN permits zero loop_q2_cut");
  }
  config->numerics.loop.r_min = std::sqrt(cut);
  config->loop_const.emplace(config->numerics.loop);
  return config;
}

// Compute the run owned immutable Durham configuration
std::shared_ptr<const MDurham::DurhamConfig> GetDurhamConfig(MModelCache &cache) {
  return cache.Get<MDurham::DurhamConfig>("durham", [&cache] { return ReadDurhamConfig(cache.Tune()); });
}

// Initialize shared Durham configuration and Sudakov state for one process
// instance
//
std::shared_ptr<const MDurham::DurhamConfig> InitializeDurhamConfig(gra::LORENTZSCALAR  &lts,
                                                                    const MModelTunePtr &tune) {
  MModelCache        &cache      = RequireModelCache(lts.model_cache, tune, "InitializeDurhamConfig");
  const SoftModelPtr &soft_model = tune->Soft();
  const auto          config     = GetDurhamConfig(cache);
  if (lts.GlobalSudakovPtr == nullptr) {
    if (lts.model_cache == nullptr) {
      throw std::invalid_argument("InitializeDurhamConfig: missing run owned model caches");
    }
    lts.GlobalSudakovPtr = lts.model_cache->sudakov.GetSudakov(lts.sqrt_s, lts.LHAPDFSET, soft_model);
  }
  if (lts.GlobalSudakovPtr->SoftModelHandle() != soft_model) {
    throw std::invalid_argument("InitializeDurhamConfig: Sudakov and Durham SOFT models differ");
  }
  if (lts.GlobalSudakovPtr->GetIRMode() == gra::MSudakovIRMode::PERTURBATIVE_ONLY &&
      config->param.loop_q2_cut < lts.GlobalSudakovPtr->GetQ2Min()) {
    throw std::invalid_argument(
        "InitializeDurhamConfig: PERTURBATIVE_ONLY "
        "requires loop_q2_cut >= Q0^2");
  }
  return config;
}

// Compute true for one direct meson pair supported by the Durham hard kernel
bool DurhamMesonPairSupported(const std::vector<MDecayBranch> &tree) {
  if (tree.size() != 2 || !tree[0].legs.empty() || !tree[1].legs.empty()) { return false; }

  const int  raw_pdg0   = tree[0].p.pdg;
  const int  raw_pdg1   = tree[1].p.pdg;
  const int  pdg0       = std::abs(raw_pdg0);
  const int  pdg1       = std::abs(raw_pdg1);
  const auto is_octet   = [](const int pdg) { return pdg == 111 || pdg == 211 || pdg == 321 || pdg == 311; };
  const auto is_singlet = [](const int pdg) { return pdg == 221 || pdg == 331; };

  if (is_singlet(pdg0) && is_singlet(pdg1)) { return true; }
  if (!(is_octet(pdg0) && is_octet(pdg1))) { return false; }
  const bool neutral_pion_pair          = raw_pdg0 == raw_pdg1 && pdg0 == 111 && pdg1 == 111;
  const bool particle_antiparticle_pair = raw_pdg0 == -raw_pdg1 && pdg0 == pdg1;
  return neutral_pion_pair || particle_antiparticle_pair;
}

// Compute the final-state pattern exposed by one analytic Durham mode
std::string DurhamProcessPattern(MDurhamMode mode) {
  switch (mode) {
    case MDurhamMode::Generic:
      return "event Durham final state";
    case MDurhamMode::Resonance:
      return "chic decay tree";
    case MDurhamMode::Parton:
      return "generated partonic states";
    case MDurhamMode::PhotonPair:
      return "photon pair";
    case MDurhamMode::MesonPair:
      return "pseudoscalar pairs";
    case MDurhamMode::Flux:
      return "arbitrary central state";
  }
  throw std::invalid_argument("MDurham::ProcessDefinitionFor: unknown process mode");
}

// Build one immutable Durham process definition for an analytic channel
std::shared_ptr<const amplitude::ProcessDefinition> MDurham::ProcessDefinitionFor(MDurhamMode        mode,
                                                                                  const std::string &process_name) {
  if (mode == MDurhamMode::Parton) {
    return std::make_shared<amplitude::ProcessRegistry>(amplitude::Processes("DURHAM"));
  }
  if (mode == MDurhamMode::PhotonPair) {
    return std::make_shared<amplitude::ProcessRegistry>(std::vector<amplitude::Process>{AMP_gg_yy::Definition()});
  }
  return std::make_shared<amplitude::AnalyticProcess>(
      "MDURHAM", process_name, DurhamProcessPattern(mode),
      DecayStructure{mode == MDurhamMode::Generic || mode == MDurhamMode::MesonPair ? DecayType::Full : DecayType::None},
      [mode](const std::vector<MDecayBranch> &tree) {
        if (mode == MDurhamMode::Generic || mode == MDurhamMode::Flux) { return true; }
        if (mode == MDurhamMode::Resonance) { return !tree.empty(); }
        return DurhamMesonPairSupported(tree);
      },
      [mode](const LORENTZSCALAR &lts) {
        if (mode == MDurhamMode::Generic || mode == MDurhamMode::MesonPair) { return DecayStructure{DecayType::Full}; }
        return mode == MDurhamMode::Resonance && lts.process.ROOT_RES_ACTIVE
                   ? decay::JacobWickStructure(lts) : DecayStructure{};
      });
}

// Build the tune-independent Durham continuum process union
std::shared_ptr<const amplitude::ProcessDefinition> MDurham::ContinuumProcesses() {
  std::vector<amplitude::Process> generated = amplitude::Processes("DURHAM");
  generated.push_back(AMP_gg_yy::Definition());

  std::vector<std::shared_ptr<const amplitude::ProcessDefinition>> definitions;
  definitions.push_back(std::make_shared<amplitude::ProcessRegistry>(std::move(generated)));
  definitions.push_back(ProcessDefinitionFor(MDurhamMode::MesonPair, "gg_MM"));
  return std::make_shared<amplitude::ProcessAlternatives>(std::move(definitions));
}

namespace {

// Compute the model tune after validating an injected immutable snapshot
MModelTunePtr RequireDurhamTune(MModelTunePtr tune) {
  if (tune == nullptr) { throw std::invalid_argument("MDurham: null model tune"); }
  return tune;
}

}  // namespace

// Initialize Durham and Sudakov parameters before worker process copies
void MDurham::InitializeParameters(MProcessSetup &setup) {
  InitializeDurhamConfig(setup.lts, RequireDurhamTune(setup.model_tune));
}

// Constructor
//
MDurham::MDurham(gra::LORENTZSCALAR& lts, MModelTunePtr model_tune_snapshot, gra::MRandom& rng_in,
                 std::shared_ptr<const amplitude::ProcessDefinition> definition)
    : amplitude::ProcessFamily(std::move(definition)),
      model_tune(
          RequireModelCache(lts.model_cache, RequireDurhamTune(std::move(model_tune_snapshot)), "MDurham").TunePtr()),
      soft_model(model_tune->Soft()),
      rng(rng_in),
      config(InitializeDurhamConfig(lts, model_tune)),
      param(config->param),
      numerics(config->numerics),
      loop_const(*config->loop_const),
      q2_flavour([&lts] {
        const auto pdf = lts.model_cache->pdf.GetPDF(lts.LHAPDFSET, 0);
        return std::array<double, 2>{pow2(pdf->info().get_entry_as<double>("MCharm")),
                                      pow2(pdf->info().get_entry_as<double>("MBottom"))};
      }()),
      generated_durham_registry([] {
        auto processes = amplitude::Processes("DURHAM");
        processes.push_back(AMP_gg_yy::Definition());
        return processes;
      }()) {}

// Clear event color-flow tags from the central decay tree
//
void MDurham::ClearCentralColorFlow(gra::LORENTZSCALAR &lts) const {
  std::function<void(gra::MDecayBranch &)> clear_branch = [&](gra::MDecayBranch &branch) {
    branch.p.color_flow.clear();
    for (auto &leg : branch.legs) { clear_branch(leg); }
  };

  for (auto &branch : lts.decaytree) { clear_branch(branch); }
}

// Compute generalized-kt jet power for the selected MadGraph jet algorithm
int MDurham::MadGraphJetPower() const {
  if (param.JET_ALGO == "anti-kt") { return -1; }
  if (param.JET_ALGO == "CA") { return 0; }
  if (param.JET_ALGO == "kt") { return 1; }
  throw std::invalid_argument("MDurham::MadGraphJetPower: unsupported JET_ALGO = " + param.JET_ALGO);
}

// Cluster partons with the selected generalized-kt MadGraph jet algorithm
// d_iB = pTi^(2p), d_ij = min(d_iB,d_jB) DeltaR_ij^2/R^2
// [REFERENCE: Cacciari, Salam and Soyez, arXiv:0802.1189, Eqs. (1)-(2)]
std::vector<gra::M4Vec> MDurham::ClusterMadGraphJets(const std::vector<gra::M4Vec> &partons) const {
  if (param.JET_ALGO == "none") { return partons; }

  const int    power = MadGraphJetPower();
  const double R2    = pow2(param.JET_R);

  auto beam_distance = [power](const gra::M4Vec &p4) {
    const double pt2 = p4.Pt2();
    if (power == 0) { return 1.0; }
    if (power > 0) { return pt2; }
    return (pt2 > 0.0) ? 1.0 / pt2 : std::numeric_limits<double>::infinity();
  };

  std::vector<gra::M4Vec> active = partons;
  std::vector<gra::M4Vec> jets;
  jets.reserve(partons.size());

  while (!active.empty()) {
    double      best_distance = std::numeric_limits<double>::infinity();
    bool        best_is_pair  = false;
    std::size_t best_i        = 0;
    std::size_t best_j        = 0;

    for (const auto &i : indices(active)) {
      const double distance = beam_distance(active[i]);
      if (distance < best_distance) {
        best_distance = distance;
        best_is_pair  = false;
        best_i        = i;
      }
    }

    for (const auto &i : indices(active)) {
      for (std::size_t j = i + 1; j < active.size(); ++j) {
        const double dy       = active[i].Rap() - active[j].Rap();
        const double dphi     = active[i].DeltaPhi(active[j]);
        const double dR2      = pow2(dy) + pow2(dphi);
        const double distance = std::min(beam_distance(active[i]), beam_distance(active[j])) * dR2 / R2;
        if (std::isfinite(distance) && distance < best_distance) {
          best_distance = distance;
          best_is_pair  = true;
          best_i        = i;
          best_j        = j;
        }
      }
    }

    if (best_is_pair) {
      active[best_i] += active[best_j];
      active.erase(active.begin() + static_cast<std::ptrdiff_t>(best_j));
    } else {
      jets.push_back(active[best_i]);
      active.erase(active.begin() + static_cast<std::ptrdiff_t>(best_i));
    }
  }

  return jets;
}

// Apply resolved generalized-kt parton-level jet cuts to MadGraph final states
bool MDurham::PassMadGraphJetCuts(const gra::LORENTZSCALAR &lts,
                                  const std::vector<int>   *final_color_representations) const {
  if (param.JET_ALGO == "none") { return true; }

  const auto leaves = mg5::StableDecayLeaves(lts.decaytree);
  if (final_color_representations != nullptr && final_color_representations->size() != leaves.size()) {
    throw std::invalid_argument("MDurham::PassMadGraphJetCuts: final-state color structure mismatch");
  }

  std::vector<gra::M4Vec> partons;
  partons.reserve(leaves.size());

  for (const auto &i : indices(leaves)) {
    if (final_color_representations != nullptr && std::abs((*final_color_representations)[i]) == 1) { continue; }
    const double pt = leaves[i]->p4.Pt();
    const double y  = leaves[i]->p4.Rap();
    if (!std::isfinite(pt) || pt < param.JET_pt_min) { return false; }
    if (!std::isfinite(y) || std::abs(y) > param.JET_rap_max) { return false; }
    partons.push_back(leaves[i]->p4);
  }

  const std::vector<gra::M4Vec> jets = ClusterMadGraphJets(partons);
  return jets.size() == partons.size();
}

// Reuse one event-local jet-cut decision across all screening-loop momenta
bool MDurham::PassCachedMadGraphJetCuts(gra::LORENTZSCALAR     &lts,
                                        const std::vector<int> *final_color_representations) const {
  if (!lts.screening.active || !lts.screening.durham.jet_cuts_cached) {
    lts.screening.durham.jet_cuts_pass   = PassMadGraphJetCuts(lts, final_color_representations);
    lts.screening.durham.jet_cuts_cached = true;
  }
  return lts.screening.durham.jet_cuts_pass;
}

// Reset projected amplitudes after a hard partonic event fails jet cuts
void MDurham::ZeroProjectedAmplitudes(DurhamProjectedAmp &Amp, std::vector<DurhamProjectedAmp> *ColorAmp) const {
  for (auto &channel : Amp) { channel.fill(0.0); }
  if (ColorAmp != nullptr) { ColorAmp->clear(); }
}

// Resolve and cache the generated matrix element matching the exact process
DurhamMG5Process &MDurham::ResolveGeneratedDurhamProcess(const gra::LORENTZSCALAR &lts) {
  const auto process = generated_durham_registry.MatchProcess(lts.decaytree);
  if (!process.has_value()) {
    throw std::invalid_argument(
        "MDurham::ResolveGeneratedDurhamProcess: no "
        "generated matrix element for "
        "the exact decay topology");
  }

  if (generated_durham_process == nullptr || !generated_durham_match.has_value() ||
      *generated_durham_match != *process) {
    if (process->process_name == "gg_yy") {
      generated_durham_process = std::make_unique<AMP_gg_yy>(model_tune->SM());
    } else {
      generated_durham_process = CreateDurhamMG5Process(*process);
    }
    if (generated_durham_process == nullptr) {
      throw std::invalid_argument(
          "MDurham::ResolveGeneratedDurhamProcess: no "
          "generated matrix element for "
          "the matched process");
    }
    generated_durham_match = process;
  }

  const auto matrix_element_process = generated_durham_process->MatchProcess(lts.decaytree);
  if (!matrix_element_process.has_value() || *matrix_element_process != *process) {
    throw std::invalid_argument(
        "MDurham::ResolveGeneratedDurhamProcess: "
        "matrix element process mismatch");
  }
  return *generated_durham_process;
}

// Sample one Durham color-flow candidate after the final event amplitude is
// known
//
bool MDurham::SampleColorFlow(gra::LORENTZSCALAR &lts) {
  ClearCentralColorFlow(lts);
  if (lts.hard_color_flows.empty()) { return true; }
  const auto candidate = [](const mg5helas::ExternalColorFlow &external) {
    std::vector<MColorFlow> out;
    out.reserve(external.size());
    for (const auto &leg : external) { out.push_back({leg.color, leg.anticolor}); }
    return out;
  };

  const auto selected = mg5helas::SelectHardColorFlow(lts.hard_color_flows, rng);
  return selected.has_value() &&
         AssignDurhamColorFlowCandidate(lts, candidate(lts.hard_color_flows[*selected].external));
}

// Durham QCD / KMR model
//
// [REFERENCE: Pumplin, Phys.Rev.D 52 (1995)]
// [REFERENCE: Khoze, Kaidalov, Martin, Ryskin, Stirling, arxiv.org/abs/hep-ph/0507040]
// [REFERENCE: Khoze, Martin, Ryskin, arxiv.org/abs/hep-ph/0605113]
// [REFERENCE: Harland-Lang, Khoze, Ryskin, Stirling, arxiv.org/abs/1005.0695]
// [REFERENCE: Harland-Lang, Khoze, Ryskin, arxiv.org/abs/1409.4785]
//
double MDurham::DurhamQCD(gra::LORENTZSCALAR& lts, const std::string& process) try {
  evaluation_status                          = mg5helas::EvaluationStatus::Success;
  const bool        top_level_call           = !lts.screening.active;
  DurhamMG5Process *generated_matrix_element = nullptr;
  if (process == "MG5" || process == "gg") { generated_matrix_element = &ResolveGeneratedDurhamProcess(lts); }
  const bool        colored_process = generated_matrix_element != nullptr || IsDurhamColoredProcess(process);
  const std::size_t channel_count =
      generated_matrix_element != nullptr
          ? generated_matrix_element->ColorRank() * (generated_matrix_element->HelicityCount() / 4)
          : DurhamChannelCount(process);
  if (top_level_call) {
    ClearCentralColorFlow(lts);
    lts.hard_color_flows.clear();
    lts.screening.durham.jet_cuts_cached = false;
    if (colored_process) {
      lts.id1      = 0;
      lts.id2      = 0;
      lts.muF      = 0.0;
      lts.muR      = 0.0;
      lts.scalup   = 0.0;
      lts.alphaQCD = 0.0;
    }
  }
  lts.screening.durham.hard_cached = false;

  if (!DurhamEventDomainValid(lts, param)) {
    evaluation_status = mg5helas::EvaluationStatus::KinematicsFailure;
    lts.hamp.assign(channel_count, 0.0);
    ZeroDurhamColorAmplitudes(lts, channel_count);
    return 0.0;
  }

  const double event_alpha_s = colored_process ? DurhamEventAlphaS(lts, param) : 0.0;

  if (top_level_call && colored_process) {
    lts.id1                          = PDG::PDG_gluon;
    lts.id2                          = PDG::PDG_gluon;
    const double        central_mass = std::sqrt(lts.s_hat);
    const MDurhamScales scales       = DurhamCentralScales(central_mass, param);
    lts.muF                          = scales.muF;
    lts.muR                          = scales.muR;
    lts.scalup                       = scales.muF;
    lts.alphaQCD                     = event_alpha_s;
  }

  // Any MG5 generated finite-Nc hard process
  if (generated_matrix_element != nullptr) {
    DurhamProjectedAmp Amp(channel_count);
    for (auto &channel : Amp) { channel.fill(0.0); }

    std::vector<DurhamProjectedAmp> color_amp;
    Dgg2Generated(lts, *generated_matrix_element, Amp, &color_amp);
    if (!lts.screening.durham.hard_cached || !mg5helas::EvaluationSucceeded(evaluation_status)) {
      lts.hamp.assign(channel_count, 0.0);
      ZeroDurhamColorAmplitudes(lts, channel_count);
      return 0.0;
    }
    if (top_level_call) {
      for (const auto &candidate : generated_matrix_element->FlowCandidates()) {
        mg5helas::ExternalColorFlow external;
        external.reserve(candidate.size());
        for (const auto &leg : candidate) { external.push_back({leg.flow1, leg.flow2}); }
        lts.hard_color_flows.push_back({{}, std::move(external)});
      }
    }
    return DQtloop(lts, Amp, &color_amp);

    // Meson pair continuum
  } else if (process == "MMbar") {
    // [Final state helicities/polarizations x 4 initial state helicities]
    DurhamProjectedAmp Amp(1);
    Amp[0].fill(0.0);

    // Amplitude evaluated outside the Qt-loop (approximation)
    Dgg2MMbar(lts, Amp);

    // Run loop
    return DQtloop(lts, Amp);

    // chic(0), chic(1), chic(2) resonances
  } else if (process == "chic(0)") {
    // [Spin projection m = 0 x dynamic contracted q_i V_ij q_j vertex]
    DurhamProjectedAmp Amp(1);
    Amp[0].fill(0.0);

    const DurhamLoopAmp loop_amp = Dgg2chic0(lts);

    if (DurhamDecaySpinActive(lts, process)) {
      return DurhamDecayAmp2(lts, process, DQtloopAmplitudes(lts, Amp, nullptr, &loop_amp));
    }

    // Run loop
    return DQtloop(lts, Amp, nullptr, &loop_amp);

  } else if (process == "chic(1)") {
    // [Final state polarizations x dynamic loop kernel]
    DurhamProjectedAmp Amp(3);
    for (auto &channel : Amp) { channel.fill(0.0); }

    // The axial-vector contracted kernel depends on the off-shell transverse
    // momenta
    const DurhamLoopAmp loop_amp = Dgg2chic1(lts);

    if (DurhamDecaySpinActive(lts, process)) {
      return DurhamDecayAmp2(lts, process, DQtloopAmplitudes(lts, Amp, nullptr, &loop_amp));
    }

    // Run loop
    return DQtloop(lts, Amp, nullptr, &loop_amp);

  } else if (process == "chic(2)") {
    // [Spin projections m = -2,-1,0,+1,+2 x dynamic contracted q_i V_ij q_j
    // vertex]
    DurhamProjectedAmp Amp(5);
    for (auto &channel : Amp) { channel.fill(0.0); }

    const DurhamLoopAmp loop_amp = Dgg2chic2(lts);

    if (DurhamDecaySpinActive(lts, process)) {
      return DurhamDecayAmp2(lts, process, DQtloopAmplitudes(lts, Amp, nullptr, &loop_amp));
    }

    // Run loop
    return DQtloop(lts, Amp, nullptr, &loop_amp);

    // Flux with |A| = 1
  } else if (process == "FLUX") {
    // [Final state helicities/polarizations x 4 initial state helicities]
    DurhamProjectedAmp Amp(1);
    Amp[0].fill(0.0);

    // Initial state gluon helicity combinations,
    // ++ and -- give contribution for 0+ state (see DHelProj function)
    const double A = 1.0;

    DurhamInitialHelicity is{};
    is.fill(0.0);
    is[DurhamHelicityIndex(-1, -1)] = A;
    is[DurhamHelicityIndex(1, 1)]   = A;

    // No polarization for the final state to loop over, only one index
    Amp[0] = is;

    // Run loop
    return DQtloop(lts, Amp);

    // ------------------------------------------------------------
    //
    // Implement more processes here ...
    //
    // ------------------------------------------------------------

  } else {
    std::string str = "MDurham::Durham: Unknown subprocess: " + process;
    throw std::invalid_argument(str);
  }
} catch (const AmplitudeFailure&) {
  evaluation_status = mg5helas::EvaluationStatus::AmplitudeFailure;
  throw;
}

// Compute the KMR transverse projector tensor components
//
// The four entries are
//   P0 = -1/2 (q1x q2x + q1y q2y)
//   M0 = -i/2 (q1x q2y - q1y q2x)
//   P2 =  1/2 [(q1x q2x - q1y q2y) + i(q1x q2y + q1y q2x)]
//   M2 =  1/2 [(q1x q2x - q1y q2y) - i(q1x q2y + q1y q2x)]
//
// In the forward limit q1_t = - q2_t = Q_t
// => gives gluon polarization vectors eps_1 = -eps_2 => central system J_z = 0
//
void MDurham::DHelicity(const DurhamTransverseMomentum &q1, const DurhamTransverseMomentum &q2,
                        std::vector<std::complex<double>> &JzP) const {
  JzP.resize(4);
  const unsigned int X = 0;  // component for readability
  const unsigned int Y = 1;

  // 1/2  q1_t dot q2_t
  JzP[P0] = -0.5 * (q1[X] * q2[X] + q1[Y] * q2[Y]);

  // Durham papers often print |q1_t x q2_t| for the 0- spin-factor magnitude
  // Here the signed transverse orientation is kept because DHelProj coherently
  // sums helicity phases and no compensating azimuthal phase is applied
  // elsewhere
  JzP[M0] = -0.5 * math::zi * (q1[X] * q2[Y] - q1[Y] * q2[X]);

  const double Re = 0.5 * (q1[X] * q2[X] - q1[Y] * q2[Y]);
  const double Im = 0.5 * (q1[X] * q2[Y] + q1[Y] * q2[X]);

  // +1/2[ (xx - yy) + i(xy + yx) ]
  JzP[P2] = Re + math::zi * Im;

  // +1/2[ (xx - yy) - i(xy + yx) ]
  JzP[M2] = Re - math::zi * Im;
}

// Contract hard g(lambda1) g(lambda2) -> X helicity amplitudes with J_z
// projectors
//
// Compute sum_J J(q1,q2) M_{\lambda1 lambda2}; the initial helicity order is
// (--,-+,+-,++), and the coherent sum follows SuperChic Eq. (20)
//
// [REFERENCE: arxiv.org/abs/1405.0018v2, formula (20)]
//
std::complex<double> MDurham::DHelProj(const std::vector<std::complex<double>> &A,
                                       const std::vector<std::complex<double>> &JzP) const {
  if (A.size() != 4 || JzP.size() != 4) {
    throw AmplitudeFailure("MDurham::DHelProj: expected four initial helicities and projectors");
  }
  // M_{++} + M_{--}
  // (this term gives J_z^PC = 0^++ selection rule in the forward pt->0 limit)
  const std::complex<double> aP0 = JzP[P0] * (A[DurhamHelicityIndex(1, 1)] + A[DurhamHelicityIndex(-1, -1)]);
  // M_{++} - M_{--}
  const std::complex<double> aM0 = JzP[M0] * (A[DurhamHelicityIndex(1, 1)] - A[DurhamHelicityIndex(-1, -1)]);
  // M_{-+}
  const std::complex<double> aP2 = JzP[P2] * A[DurhamHelicityIndex(-1, 1)];
  // M_{+-}
  const std::complex<double> aM2 = JzP[M2] * A[DurhamHelicityIndex(1, -1)];

  /*
  // DEBUG
  std::cout << "0+ : " << aP0 << std::endl;
  std::cout << "0- : " << aM0 << std::endl;
  std::cout << "2+ : " << aP2 << std::endl;
  std::cout << "2- : " << aM2 << std::endl;
  std::cout << std::endl << std::endl;
  */

  // Coherent sum!
  return aP0 + aM0 + aP2 + aM2;
}

// Convenience overload for fixed-size initial helicity arrays
std::complex<double> MDurham::DHelProj(const DurhamInitialHelicity             &A,
                                       const std::vector<std::complex<double>> &JzP) const {
  if (JzP.size() != 4) { throw AmplitudeFailure("MDurham::DHelProj: expected four projectors"); }
  const std::complex<double> aP0 = JzP[P0] * (A[DurhamHelicityIndex(1, 1)] + A[DurhamHelicityIndex(-1, -1)]);
  const std::complex<double> aM0 = JzP[M0] * (A[DurhamHelicityIndex(1, 1)] - A[DurhamHelicityIndex(-1, -1)]);
  const std::complex<double> aP2 = JzP[P2] * A[DurhamHelicityIndex(-1, 1)];
  const std::complex<double> aM2 = JzP[M2] * A[DurhamHelicityIndex(1, -1)];
  return aP0 + aM0 + aP2 + aM2;
}

// Alternative (semi-ad-hoc) scenarios for the scale choise
void MDurham::DScaleChoise(double qt2, double q1_2, double q2_2, double &Q1_2_scale, double &Q2_2_scale) const {
  if (param.PDF_scale == "MIN") {
    Q1_2_scale = std::min(qt2, q1_2);
    Q2_2_scale = std::min(qt2, q2_2);
  } else if (param.PDF_scale == "MAX") {
    Q1_2_scale = std::max(qt2, q1_2);
    Q2_2_scale = std::max(qt2, q2_2);
  } else if (param.PDF_scale == "IN") {
    Q1_2_scale = q1_2;
    Q2_2_scale = q2_2;
  } else if (param.PDF_scale == "EX") {
    Q1_2_scale = qt2;
    Q2_2_scale = qt2;
  } else if (param.PDF_scale == "AVG") {
    Q1_2_scale = (qt2 + q1_2) / 2.0;
    Q2_2_scale = (qt2 + q2_2) / 2.0;
  } else {
    throw std::invalid_argument("MDurham::DScaleChoise: Unknown 'Durham::PDF_scale' option!");
  }
}

// Evaluate the reduced Durham loop and match it to the invariant pp amplitude
//
// For each final helicity channel h this computes
//
//   T_Mbar^h = pi^2 int d^2Q_t f_g1 f_g2 [2 q_i V_h,ij q_j / M_X^2]
//              / (Q_t^2 q1_t^2 q2_t^2)
//
// Generic helicity amplitudes use the steered transverse-basis or JzP
// contraction. Dynamic off-shell kernels, such as the chi_cJ vertices, return
// q_i V_ij q_j directly. After applying the common Mbar normalization and doing
// the numerical integral this function returns
//
// A_pp[h] = F_pp(t_1) F_pp(t_2) s T_Mbar^h,
//
// where F_pp is the proton form factor. DQtloop stores this directly in
// lts.hamp, while Durham resonance decays first contract it with decay_f
// There is no initial 1/4 helicity average,
// because protons are treated as unpolarized (also screening should not apply
// that)
//
// The integration measure is d^2Q_t = q_t dq_t dphi over phi in [0,2pi]; no
// additional symmetry factor is present in the KMR normalization
// No extra global i is inserted here: the Durham Born amplitude is kept in the
// hard-subprocess phase convention, while the screening loop supplies the
// i/(8 pi^2 s) prefactor through MEikonal loop weights
//
// [REFERENCE: Khoze, Martin, Ryskin, journals.aps.org/prd/pdf/10.1103/PhysRevD.56.5867]
// [REFERENCE: Harland-Lang, Khoze, Martin, Ryskin, Stirling, arxiv.org/abs/1405.0018v2]
//
// See also:
// [REFERENCE: Lonnblad, Zlebcik, arxiv.org/abs/1608.03765]
//
std::vector<std::complex<double>> MDurham::DQtloopAmplitudes(gra::LORENTZSCALAR& lts, const DurhamProjectedAmp& Amp,
                                                             const std::vector<DurhamProjectedAmp>* ColorAmp,
                                                             const DurhamLoopAmp*                   LoopAmp) try {
  // Amp.size == 4 for example for the gg -> gg process [initial states summed
  // coherently]
  std::vector<std::complex<double>> sum(Amp.size(), 0.0);
  if (!DurhamEventDomainValid(lts, param)) {
    evaluation_status = mg5helas::EvaluationStatus::KinematicsFailure;
    ZeroDurhamColorAmplitudes(lts, Amp.size());
    return sum;
  }

  // Forward proton transverse momenta after validating the final state
  const DurhamTransverseMomentum pt1 = {lts.pfinal[1].Px(), lts.pfinal[1].Py()};
  const DurhamTransverseMomentum pt2 = {lts.pfinal[2].Px(), lts.pfinal[2].Py()};

  M4Vec      hard_k1;
  M4Vec      hard_k2;
  const bool use_transverse_projector = numerics.helicity_projector == "transverse";
  if (use_transverse_projector && LoopAmp == nullptr) {
    if (lts.screening.durham.hard_cached) {
      hard_k1 = lts.screening.durham.hard_k1;
      hard_k2 = lts.screening.durham.hard_k2;
    } else {
      std::vector<M4Vec> hard_final;
      hard_final.reserve(lts.decaytree.size());
      for (const auto &branch : lts.decaytree) { hard_final.push_back(branch.p4); }
      if (!mg5helas::PrepareOnShellKinematics(lts, hard_final, hard_k1, hard_k2)) {
        evaluation_status = mg5helas::EvaluationStatus::KinematicsFailure;
        ZeroDurhamColorAmplitudes(lts, Amp.size());
        return sum;
      }
    }
  }

  // *************************************************************************
  // ** Sudakov suppression kt^2 integral upper bound in GeV **
  const MDurhamScales scales = DurhamCentralScales(msqrt(lts.s_hat), param);
  // *************************************************************************
  const double           qv_to_mbar = 2.0 / lts.s_hat;
  auto                   flux_1     = lts.GlobalSudakovPtr->PrepareFlux(lts.x1, scales.muF);
  auto                   flux_2     = lts.GlobalSudakovPtr->PrepareFlux(lts.x2, scales.muF);

  const bool color_matches_total = ColorAmp != nullptr && ColorAmp->size() == 1 && ColorAmp->front() == Amp;
  std::vector<std::vector<std::complex<double>>> color_sum;
  if (ColorAmp != nullptr && !color_matches_total) {
    color_sum.resize(ColorAmp->size());
    for (std::size_t c = 0; c < ColorAmp->size(); ++c) {
      if ((*ColorAmp)[c].size() != Amp.size()) {
        throw AmplitudeFailure("MDurham::DQtloop: Color amplitude helicity dimension mismatch");
      }
      color_sum[c].assign((*ColorAmp)[c].size(), 0.0);
    }
  }

  // Spin-Parity
  std::vector<std::complex<double>>     JzP(4, 0.0);
  std::array<std::complex<double>, 4>   integrated_projector         = {};
  std::array<double, 4>                 integrated_transverse_tensor = {};
  mg5helas::TransverseHelicityProjector source_projector_1;
  mg5helas::TransverseHelicityProjector source_projector_2;
  if (use_transverse_projector && LoopAmp == nullptr) {
    source_projector_1 = mg5helas::PrepareTransverseHelicityProjector(hard_k1, hard_k2);
    source_projector_2 = mg5helas::PrepareTransverseHelicityProjector(hard_k2, hard_k1);
  }

  // 2D-loop integral
  //
  // \int d^2 \vec{qt} [...] = \int dphi \int dqt qt [...]
  //

  const bool   pdf_scale_min = param.PDF_scale == "MIN";
  const double q2_max        = lts.GlobalSudakovPtr->GetQ2Max();
  // Tie the azimuth origin to the recoil plane to preserve rotations and beam exchange
  const double dx            = pt1[0] - pt2[0];
  const double dy            = pt1[1] - pt2[1];
  const double phi0          = (dx * dx + dy * dy > 0.0) ? std::atan2(dy, dx) : std::atan2(pt1[1], pt1[0]);
  // Split at infrared matching, flavour thresholds, radiation and PDF table boundaries
  std::vector<std::array<double, 3>> splits;
  const double q2_min = lts.GlobalSudakovPtr->GetQ2Min();
  for (const double scale2 : {q2_min, q2_flavour[0], q2_flavour[1], scales.muF * scales.muF, q2_max}) {
    if (scale2 < q2_min || (pdf_scale_min && scale2 >= numerics.qt2_MAX)) { continue; }
    if (param.PDF_scale != "IN") { splits.push_back({0.0, 0.0, std::sqrt(scale2)}); }
    if (param.PDF_scale == "EX") { continue; }
    for (const auto& centre : {pt1, DurhamTransverseMomentum{-pt2[0], -pt2[1]}}) {
      if (param.PDF_scale == "AVG") {
        const double r2 = scale2 - DurhamTransverseMomentum2(centre) / 4.0;
        if (r2 > 0.0) { splits.push_back({centre[0] / 2.0, centre[1] / 2.0, std::sqrt(r2)}); }
      } else {
        splits.push_back({centre[0], centre[1], std::sqrt(scale2)});
      }
    }
  }
  const auto nodes = loop_const.Nodes({pt1, {-pt2[0], -pt2[1]}}, std::sqrt(param.loop_q2_cut), phi0, splits);
  unsigned int ring          = 0;
  double       qt2           = 0.0;
  double       radial_flux_1 = 0.0;
  double       radial_flux_2 = 0.0;
  for (const auto& node : nodes) {
    if (node.radial != ring) {
      ring = node.radial;
      qt2  = std::min(node.x * node.x + node.y * node.y, numerics.qt2_MAX);
      if (pdf_scale_min && qt2 <= q2_max) {
        radial_flux_1 = flux_1.OverQ2(qt2);
        radial_flux_2 = flux_2.OverQ2(qt2);
      }
    }
    const double qt_x = node.x;
    const double qt_y = node.y;
    // Fusing gluon pt-vectors
    const DurhamTransverseMomentum q1 = {qt_x - pt1[0], qt_y - pt1[1]};    //  qt - pt1
    const DurhamTransverseMomentum q2 = {-qt_x - pt2[0], -qt_y - pt2[1]};  // -qt - pt2

    const double q1_2 = DurhamTransverseMomentum2(q1);
    const double q2_2 = DurhamTransverseMomentum2(q2);

    // Omit only measure-zero singular nodes, never a finite infrared interval
    if (!std::isfinite(qt2) || !std::isfinite(q1_2) || !std::isfinite(q2_2) || qt2 <= 0.0 || q1_2 <= 0.0 ||
        q2_2 <= 0.0) {
      continue;
    }

    if (!use_transverse_projector || LoopAmp != nullptr) {
      // Get fusing gluon spin-parity components [q1,q2] -> [0+,0-,+2,-2]
      DHelicity(q1, q2, JzP);
    }

    // ** Durham scale choise **
    const bool q1_uses_qt = pdf_scale_min && qt2 <= q1_2;
    const bool q2_uses_qt = pdf_scale_min && qt2 <= q2_2;
    double     Q1_2_scale = q1_uses_qt ? qt2 : q1_2;
    double     Q2_2_scale = q2_uses_qt ? qt2 : q2_2;
    if (!pdf_scale_min) { DScaleChoise(qt2, q1_2, q2_2, Q1_2_scale, Q2_2_scale); }

    // Enforce only the tabulated upper scale boundary
    if (!std::isfinite(Q1_2_scale) || !std::isfinite(Q2_2_scale) || Q1_2_scale < 0.0 || Q2_2_scale < 0.0 ||
        Q1_2_scale > q2_max || Q2_2_scale > q2_max) {
      continue;
    }

    // --------------------------------------------------------------------------
    // Factor the color-neutral Q2 zero from each amplitude-level gluon
    const double fg_over_q2_1 = q1_uses_qt ? radial_flux_1 : flux_1.OverQ2(Q1_2_scale);
    const double fg_over_q2_2 = q2_uses_qt ? radial_flux_2 : flux_2.OverQ2(Q2_2_scale);

    double propagator_factor = 0.0;
    if (pdf_scale_min) {
      if (q1_uses_qt && q2_uses_qt) {
        propagator_factor = qt2 / (q1_2 * q2_2);
      } else if (q1_uses_qt) {
        propagator_factor = 1.0 / q1_2;
      } else if (q2_uses_qt) {
        propagator_factor = 1.0 / q2_2;
      } else {
        propagator_factor = 1.0 / qt2;
      }
    } else {
      propagator_factor = Q1_2_scale * Q2_2_scale / (qt2 * q1_2 * q2_2);
    }

    // Amplitude weight:
    // * \pi^2 : KMR Durham loop normalization
    // * weight: quadrature and jacobian of d^2qt -> dphi dqt qt
    const double weight = fg_over_q2_1 * fg_over_q2_2 * propagator_factor * math::PIPI * node.weight;

    if (LoopAmp != nullptr) {
      const std::vector<std::complex<double>> dynamic_amp = (*LoopAmp)(q1, q2, JzP);
      if (dynamic_amp.size() != sum.size()) {
        throw AmplitudeFailure("MDurham::DQtloop: dynamic loop amplitude size mismatch");
      }
      for (const auto& h : indices(sum)) { sum[h] += weight * qv_to_mbar * dynamic_amp[h]; }
    } else {
      if (use_transverse_projector) {
        if (!source_projector_1.Accepts(q1[0], q1[1]) || !source_projector_2.Accepts(q2[0], q2[1])) { continue; }
        const double tensor_weight = weight * qv_to_mbar;
        integrated_transverse_tensor[0] += tensor_weight * q1[0] * q2[0];
        integrated_transverse_tensor[1] += tensor_weight * q1[0] * q2[1];
        integrated_transverse_tensor[2] += tensor_weight * q1[1] * q2[0];
        integrated_transverse_tensor[3] += tensor_weight * q1[1] * q2[1];
      } else {
        for (const auto& component : indices(integrated_projector)) {
          integrated_projector[component] += weight * qv_to_mbar * JzP[component];
        }
      }
    }
  }

  if (LoopAmp == nullptr) {
    if (use_transverse_projector) {
      for (std::size_t h1 = 0; h1 < 2; ++h1) {
        for (std::size_t h2 = 0; h2 < 2; ++h2) {
          integrated_projector[spin::BinaryPairHelicityIndex(h1, h2)] =
              integrated_transverse_tensor[0] * source_projector_1.coefficient[0][h1] *
                  source_projector_2.coefficient[0][h2] +
              integrated_transverse_tensor[1] * source_projector_1.coefficient[0][h1] *
                  source_projector_2.coefficient[1][h2] +
              integrated_transverse_tensor[2] * source_projector_1.coefficient[1][h1] *
                  source_projector_2.coefficient[0][h2] +
              integrated_transverse_tensor[3] * source_projector_1.coefficient[1][h1] *
                  source_projector_2.coefficient[1][h2];
        }
      }
    }
    if (use_transverse_projector) {
      for (const auto &h : indices(sum)) { sum[h] = mg5helas::ContractHelicityKernel(Amp[h], integrated_projector); }
      for (const auto &c : indices(color_sum)) {
        for (const auto &h : indices(color_sum[c])) {
          color_sum[c][h] = mg5helas::ContractHelicityKernel((*ColorAmp)[c][h], integrated_projector);
        }
      }
    } else {
      const std::vector<std::complex<double>> integrated_jzp(integrated_projector.begin(), integrated_projector.end());
      for (const auto &h : indices(sum)) { sum[h] = DHelProj(Amp[h], integrated_jzp); }
      for (const auto &c : indices(color_sum)) {
        for (const auto &h : indices(color_sum[c])) { color_sum[c][h] = DHelProj((*ColorAmp)[c][h], integrated_jzp); }
      }
    }
  }

  // Match the reduced Durham loop to the invariant pp amplitude
  std::complex<double> external       = 1.0;
  const SoftExchangeId pomeron        = soft_model->ForwardExcitationExchange();
  const auto           forward_factor = [&](const ForwardLegState &state) {
    return state.IsExcited() ? soft_model->ForwardExcitationFactor(pomeron, state.t, state.mass2)
                                       : soft_model->NormalizedPhysicalResidue(pomeron, state.t);
  };
  external *= forward_factor(ResolveForwardLegState(lts, ForwardBeamLeg::Upper));
  external *= forward_factor(ResolveForwardLegState(lts, ForwardBeamLeg::Lower));
  external *= lts.s;

  // Outgoing helicity combinations
  for (auto &amplitude : sum) { amplitude *= external; }

  if (color_matches_total) {
    if (lts.hard_color_flows.size() != 1) { throw AmplitudeFailure("MDurham::DQtloop: color-flow layout is missing"); }
    lts.hard_color_flows.front().amplitudes = sum;
    lts.hard_color_flows.front().screened_weight.reset();
  } else if (ColorAmp != nullptr) {
    if (lts.hard_color_flows.size() != color_sum.size()) {
      throw AmplitudeFailure("MDurham::DQtloop: color-flow count changed");
    }
    for (const auto &c : indices(color_sum)) {
      auto &amplitudes = lts.hard_color_flows[c].amplitudes;
      amplitudes.resize(color_sum[c].size(), 0.0);
      lts.hard_color_flows[c].screened_weight.reset();
      for (const auto &h : indices(color_sum[c])) { amplitudes[h] = external * color_sum[c][h]; }
    }
  }
  // Initial state helicity average 1/4 not needed here

  return sum;
} catch (const std::logic_error& error) {
  evaluation_status = mg5helas::EvaluationStatus::AmplitudeFailure;
  throw AmplitudeFailure(error.what());
}

// Evaluate the reduced Durham loop and store channel amplitudes in lts.hamp
//
double MDurham::DQtloop(gra::LORENTZSCALAR &lts, const DurhamProjectedAmp &Amp,
                        const std::vector<DurhamProjectedAmp> *ColorAmp, const DurhamLoopAmp *LoopAmp) {
  lts.hamp = DQtloopAmplitudes(lts, Amp, ColorAmp, LoopAmp);

  return gra::SquaredNorm(lts.hamp);
}

// ======================================================================
// Covariant chi_cJ off-shell vertices
//
// The spin basis is the same CM spin-projection basis used by the decay
// machinery, but it is canonically boosted once into the pp working frame
// All Lorentz contractions are then evaluated directly in that frame
//
// The Durham loop momenta are embedded as
//
//   q_i^mu = (0, q_ix, q_iy, 0),
//
// so no Q_t-dependent Lorentz boost is needed

// Compute the q_i V_ij q_j representation of the reduced chi_c0 amplitude
//
// [REFERENCE: Harland-Lang et al., arxiv.org/abs/1405.0018v2, Eq. (26)]
// [REFERENCE: Harland-Lang et al., arxiv.org/abs/1508.02718v2, Eq. (30)]
//
MDurham::DurhamLoopAmp MDurham::Dgg2chic0(const gra::LORENTZSCALAR &lts) const {
  const DurhamCharmoniumState state               = DurhamCharmonium(lts, "chic(0)");
  const DurhamCharmoniumState normalization_state = DurhamCharmoniumNormalizationState(lts);
  const double alpha_phi_prime = DurhamCharmoniumAlphaPhiPrime(normalization_state, param.chic0_gg_width_fraction);
  const double central_mass2   = lts.s_hat;
  const std::complex<double> breit_wigner =
      gra::resonance::ComplexBreitWignerAmplitude(lts.s_hat, state.mass, state.width);

  return [state, alpha_phi_prime, central_mass2, breit_wigner](const DurhamTransverseMomentum &q1,
                                                               const DurhamTransverseMomentum &q2,
                                                               const std::vector<std::complex<double>> &) {
    if (!DurhamLoopMomentaValid(q1, q2)) { return std::vector<std::complex<double>>(1, 0.0); }

    const double q1_2 = DurhamTransverseMomentum2(q1);
    const double q2_2 = DurhamTransverseMomentum2(q2);

    // Exact Lorentz invariants of the transverse embeddings
    const double q1_sq = -q1_2;
    const double q2_sq = -q2_2;
    const double q1q2  = DurhamTransverseMinkowskiDot(q1, q2);

    const std::complex<double> c_chi = DurhamCharmoniumCoupling(breit_wigner, state, q1_2, q2_2, alpha_phi_prime);

    const std::complex<double> mbar = std::sqrt(1.0 / 6.0) * c_chi / state.mass *
                                      (3.0 * pow2(state.mass) * q1q2 - q1q2 * (q1_sq + q2_sq) - 2.0 * q1_sq * q2_sq);

    return std::vector<std::complex<double>>{DurhamMbarToQV(central_mass2, mbar)};
  };
}

// Compute the q_i V_ij q_j representation of the reduced chi_c1 amplitude
//
// [REFERENCE: Harland-Lang et al., arxiv.org/abs/1405.0018v2, Eq. (27)]
// [REFERENCE: Harland-Lang et al., arxiv.org/abs/1508.02718v2, Eq. (31)]
//
MDurham::DurhamLoopAmp MDurham::Dgg2chic1(const gra::LORENTZSCALAR &lts) const {
  const DurhamCharmoniumState state               = DurhamCharmonium(lts, "chic(1)");
  const DurhamCharmoniumState normalization_state = DurhamCharmoniumNormalizationState(lts);
  const double                s                   = lts.s;
  const double alpha_phi_prime = DurhamCharmoniumAlphaPhiPrime(normalization_state, param.chic0_gg_width_fraction);
  const double central_mass2   = lts.s_hat;
  const std::complex<double> breit_wigner =
      gra::resonance::ComplexBreitWignerAmplitude(lts.s_hat, state.mass, state.width);

  const gra::MDirac dirac;

  // Preserve the CM spin-projection basis, but express it covariantly in the
  // current pp frame once per event
  const DurhamSpin1Basis spin_basis = DurhamSpin1Polarizations(dirac, lts.pfinal[0]);

  const DurhamAxialProjector axial_projector = DurhamAxialBeamProjector(lts.pbeam1, lts.pbeam2, spin_basis);

  return [state, s, alpha_phi_prime, central_mass2, breit_wigner, axial_projector](
             const DurhamTransverseMomentum &q1, const DurhamTransverseMomentum &q2,
             const std::vector<std::complex<double>> &) {
    if (!std::isfinite(s) || s <= 0.0 || !DurhamLoopMomentaValid(q1, q2)) {
      return std::vector<std::complex<double>>(3, 0.0);
    }

    const double q1_2 = DurhamTransverseMomentum2(q1);
    const double q2_2 = DurhamTransverseMomentum2(q2);

    // q_i^2 = -|q_iT|^2 and q_{i,x/y} = -q_i[x/y]
    //
    // [q2_mu q1^2 - q1_mu q2^2] for mu=x,y therefore becomes
    //
    //   q2_i |q1T|^2 - q1_i |q2T|^2
    const double qterm_x = q2[0] * q1_2 - q1[0] * q2_2;
    const double qterm_y = q2[1] * q1_2 - q1[1] * q2_2;

    const std::complex<double> c_chi = DurhamCharmoniumCoupling(breit_wigner, state, q1_2, q2_2, alpha_phi_prime);

    std::vector<std::complex<double>> out(3, 0.0);
    for (const auto &index : indices(out)) {
      const std::complex<double> contraction =
          qterm_x * axial_projector[index][0] + qterm_y * axial_projector[index][1];

      const std::complex<double> mbar = -2.0 * zi * c_chi * contraction / s;

      out[index] = DurhamMbarToQV(central_mass2, mbar);
    }
    return out;
  };
}

// Compute the q_i V_ij q_j representation of the reduced chi_c2 amplitude
//
// [REFERENCE: Harland-Lang et al., arxiv.org/abs/1405.0018v2, Eq. (28)]
// [REFERENCE: Harland-Lang et al., arxiv.org/abs/1508.02718v2, Eq. (32)]
//
MDurham::DurhamLoopAmp MDurham::Dgg2chic2(const gra::LORENTZSCALAR &lts) const {
  const DurhamCharmoniumState state               = DurhamCharmonium(lts, "chic(2)");
  const DurhamCharmoniumState normalization_state = DurhamCharmoniumNormalizationState(lts);
  const double                s                   = lts.s;
  const double alpha_phi_prime = DurhamCharmoniumAlphaPhiPrime(normalization_state, param.chic0_gg_width_fraction);
  const double central_mass2   = lts.s_hat;
  const std::complex<double> breit_wigner =
      gra::resonance::ComplexBreitWignerAmplitude(lts.s_hat, state.mass, state.width);

  const gra::MDirac dirac;

  // Canonically boost the CM spin basis once.  Building spin-2 from these
  // boosted vectors is exactly the Lorentz transform of the CM spin-2 tensor
  const DurhamSpin1Basis spin1_basis = DurhamSpin1Polarizations(dirac, lts.pfinal[0]);
  const DurhamSpin2Basis spin2_basis = DurhamSpin2Polarizations(spin1_basis);

  const DurhamTensorProjector tensor_projector = DurhamSpin2Projector(lts.pbeam1, lts.pbeam2, spin2_basis);

  // Beam contraction is Q_t independent
  std::array<std::complex<double>, 5> beam_over_s = {};
  if (std::isfinite(s) && s > 0.0) {
    for (const auto &m : indices(beam_over_s)) { beam_over_s[m] = tensor_projector.beam[m] / s; }
  }

  return [state, s, alpha_phi_prime, central_mass2, breit_wigner, tensor_projector, beam_over_s](
             const DurhamTransverseMomentum &q1, const DurhamTransverseMomentum &q2,
             const std::vector<std::complex<double>> &) {
    if (!std::isfinite(s) || s <= 0.0 || !DurhamLoopMomentaValid(q1, q2)) {
      return std::vector<std::complex<double>>(5, 0.0);
    }

    const double q1_2 = DurhamTransverseMomentum2(q1);
    const double q2_2 = DurhamTransverseMomentum2(q2);
    const double q1q2 = DurhamTransverseMinkowskiDot(q1, q2);

    const std::complex<double> c_chi = DurhamCharmoniumCoupling(breit_wigner, state, q1_2, q2_2, alpha_phi_prime);

    std::vector<std::complex<double>> out(5, 0.0);
    for (const auto &index : indices(out)) {
      const auto &T = tensor_projector.transverse[index];

      // q1_mu q2_nu eps^{*mu nu}.  Since both lowered transverse components
      // carry a minus sign, the two signs cancel
      const std::complex<double> q_contraction =
          q1[0] * (q2[0] * T[0][0] + q2[1] * T[0][1]) + q1[1] * (q2[0] * T[1][0] + q2[1] * T[1][1]);

      const std::complex<double> contraction = q_contraction + 2.0 * q1q2 * beam_over_s[index];

      const std::complex<double> mbar = std::sqrt(2.0) * c_chi * state.mass * contraction;

      out[index] = DurhamMbarToQV(central_mass2, mbar);
    }
    return out;
  };
}

// ======================================================================
// gg -> gg tree-level helicity amplitudes
//
//
// Basic result: d\hat{\sigma}/dt = 9/4 \pi \alpha_s^2 / E_T^4
//
// In pure gluon amplitudes: when gluon helicities are the same,
// or at most one is different from the rest, vanish for any n >= 4,
// where n is the total number of gluons (in+out)
//
// [REFERENCE: Dixon, arxiv.org/abs/1310.5353v1]
//
// Fill any registry-generated finite-Nc Durham hard process
//
// The generated projector P obeys P^dagger P = G_singlet, where G_singlet
// is MadGraph's color Gram matrix after contraction with delta_ab/sqrt(8)
// [REFERENCE: Alwall et al., JHEP 07 (2014) 079, arXiv:1405.0301]
void MDurham::Dgg2Generated(gra::LORENTZSCALAR &lts, DurhamMG5Process &matrix_element, DurhamProjectedAmp &Amp,
                            std::vector<DurhamProjectedAmp> *ColorAmp) {
  evaluation_status                = mg5helas::EvaluationStatus::Success;
  lts.screening.durham.hard_cached = false;
  const auto leaves = mg5::StableDecayLeaves(lts.decaytree);
  if (leaves.size() != matrix_element.FinalPDGs().size()) {
    evaluation_status = mg5helas::EvaluationStatus::AmplitudeFailure;
    ZeroProjectedAmplitudes(Amp, ColorAmp);
    return;
  }
  if (!PassCachedMadGraphJetCuts(lts, &matrix_element.FinalColorRepresentations())) {
    ZeroProjectedAmplitudes(Amp, ColorAmp);
    return;
  }

  gra::LORENTZSCALAR work = lts;
  work.GlobalSudakovPtr   = nullptr;
  work.GlobalPdfPtr       = nullptr;

  const double alpha_s              = DurhamHardAlphaS(lts, param);
  generated_durham_evaluation       = matrix_element.Evaluate(work, alpha_s, &lts.screening.durham.hard_k1, &lts.screening.durham.hard_k2);
  evaluation_status                 = generated_durham_evaluation.status;
  if (!generated_durham_evaluation.Valid()) {
    ZeroProjectedAmplitudes(Amp, ColorAmp);
    return;
  }

  const std::size_t helicity_count       = matrix_element.HelicityCount();
  const std::size_t final_helicity_count = helicity_count / 4;
  const std::size_t rank                 = matrix_element.ColorRank();
  const auto       &evaluation           = generated_durham_evaluation;
  if (helicity_count == 0 || helicity_count % 4 != 0 || evaluation.projected.size() != rank * helicity_count) {
    evaluation_status = mg5helas::EvaluationStatus::AmplitudeFailure;
    ZeroProjectedAmplitudes(Amp, ColorAmp);
    return;
  }

  const double normalization = DurhamIncomingColorAverageFactor();
  const auto   phase         = DurhamHelicityPhases(lts, numerics.helicity_projector == "jzp");
  Amp.clear();
  Amp.reserve(rank * final_helicity_count);
  for (std::size_t color = 0; color < rank; ++color) {
    const DurhamProjectedAmp channels =
        BuildDurhamProjectedChannels(evaluation.projected.data() + color * helicity_count, helicity_count,
                                     final_helicity_count, normalization, phase);
    Amp.insert(Amp.end(), channels.begin(), channels.end());
  }

  if (ColorAmp == nullptr) {
    lts.screening.durham.hard_cached = true;
    return;
  }
  if (evaluation.flow_projected.size() != matrix_element.FlowCandidates().size()) {
    evaluation_status = mg5helas::EvaluationStatus::AmplitudeFailure;
    ZeroProjectedAmplitudes(Amp, ColorAmp);
    return;
  }
  ColorAmp->assign(evaluation.flow_projected.size(), {});
  for (const auto &flow : indices(evaluation.flow_projected)) {
    if (evaluation.flow_projected[flow].size() != rank * helicity_count) {
      evaluation_status = mg5helas::EvaluationStatus::AmplitudeFailure;
      ZeroProjectedAmplitudes(Amp, ColorAmp);
      return;
    }
    auto &candidate = (*ColorAmp)[flow];
    candidate.reserve(Amp.size());
    for (std::size_t color = 0; color < rank; ++color) {
      const DurhamProjectedAmp channels =
          BuildDurhamProjectedChannels(evaluation.flow_projected[flow].data() + color * helicity_count, helicity_count,
                                       final_helicity_count, normalization, phase);
      candidate.insert(candidate.end(), channels.begin(), channels.end());
    }
  }
  lts.screening.durham.hard_cached = true;
}

// ======================================================================
// gg -> meson pair (gg -> q\bar{q} q\bar{q})
//
// [REFERENCE: Harland-Lang, Khoze, Ryskin, Stirling, arxiv.org/pdf/1105.1626.pdf]

// Compute the fixed low-scale Chernyak-Zhitnitsky meson wavefunction
// phi_CZ(x) = 5 sqrt(3) f_M x(1-x)(1-2x)^2
// x the longitudinal momentum fraction of a parton within the meson
// fM the meson decay constant
// The fixed low-scale CZ ansatz lacks ERBL evolution, limiting high-mass predictions
// Kaons share this symmetric shape without SU(3)-breaking odd moments or narrowing
// Eta and eta-prime amplitudes omit valence-gluon Fock components
//
// Normalization: \int_0^1 dx phi(x) = fM/(2 sqrt(3))
//
double MDurham::phi_CZ(double x, double fM) const {
  return 5.0 * std::sqrt(3.0) * fM * x * (1.0 - x) * pow2(1.0 - 2.0 * x);
}

// Tabulate the quark-Fock meson wave function
//
// ----------------------------------------------------------------------
// Three charge neutral state are in the SU(3)_F quark model nonet
// (octet+singlet):
// pi0, eta8, eta0
//
// The quark-Fock decay constants use the two-angle physical basis:
//
// f_eta^8  =  f_8 * cos(theta_8)
// f_eta'^8 =  f_8 * sin(theta_8)
// f_eta^1  = -f_1 * sin(theta_1)
// f_eta'^1 =  f_1 * cos(theta_1)
//
// Valence-gluon eta/eta-prime Fock components are outside this quark kernel
//
std::vector<double> MDurham::EvalPhi(const std::vector<double> &xval, int pdg) const {
  if (xval.empty()) { throw std::invalid_argument("MDurham::EvalPhi: empty integration grid"); }

  const int apdg = std::abs(pdg);

  // Meson decay constants
  double fM = 0.0;
  if (apdg == 111 || apdg == 211) {
    fM = param.f_pi;
  } else if (apdg == 321 || apdg == 311) {
    fM = PDG::fM_meson.at(apdg);
  } else if (apdg != 221 && apdg != 331) {
    throw std::invalid_argument("MDurham::EvalPhi: unsupported meson PDG " + std::to_string(pdg));
  } else {
    fM = param.f_pi;
  }

  const double fM_eta8 = param.f_pi * param.f_eta8_over_fpi;
  const double fM_eta0 = param.f_pi * param.f_eta0_over_fpi;

  // Use f_eta^8 = f8 cos(theta8), f_eta^1 = -f1 sin(theta1),
  // f_eta'^8 = f8 sin(theta8), f_eta'^1 = f1 cos(theta1)
  std::vector<double>   f(xval.size());
  const DurhamEtaMixing mixing = DurhamEtaMixingCoefficients(apdg, param.eta_theta8_deg, param.eta_theta1_deg);
  for (const auto &i : aux::indices(f)) {
    const double x = xval[i];
    if (apdg == 221 || apdg == 331) {  // eta or eta-prime via two-angle mixing
      f[i] = mixing.octet * phi_CZ(x, fM_eta8) + mixing.singlet * phi_CZ(x, fM_eta0);
    } else {
      f[i] = phi_CZ(x, fM);
    }
  }
  return f;
}

// Compute cached meson light-cone wave functions for the active integration grid
//
const std::vector<double> &MDurham::CachedMesonWaveFunction(const std::vector<double> &xval, int pdg) const {

  const MesonWaveCacheKey     key{pdg, xval};
  std::lock_guard<std::mutex> lock(meson_wave_cache_mutex);
  auto                        it = meson_wave_cache.find(key);
  if (it == meson_wave_cache.end()) { it = meson_wave_cache.emplace(key, EvalPhi(xval, pdg)).first; }
  return it->second;
}

// Fill gg -> M Mbar hard amplitudes from the Durham meson-pair kernels
//
// This is a leading-twist massless hard kernel and should only be used above
// the hard pT cut Its unevolved meson distributions limit quantitative
// extrapolation beyond moderate masses
//
// Compute T_{\lambda1\lambda2} = (64 pi^2 alpha_s^2)/(N_C shat)
// int dx dy phi_M(x) phi_Mbar(y) T_{\lambda1\lambda2}(x,y,theta)
void MDurham::Dgg2MMbar(const gra::LORENTZSCALAR &lts, DurhamProjectedAmp &Amp) {
  Amp.assign(1, DurhamInitialHelicity{});
  // Open Gauss-Legendre nodes avoid removable boundary singularities in the
  // hard kernels
  const unsigned int Nx      = numerics.N_x;
  const auto [xval, xweight] = math::GaussLegendreRule(Nx, 0.0, 1.0);

  const double NC = 3;                              // Three colors
  const double CF = (pow2(NC) - 1.0) / (NC * 2.0);  // SU(3) algebra

  const double alpha_s = DurhamHardAlphaS(lts, param);

  const bool         transverse = numerics.helicity_projector == "transverse";
  std::vector<M4Vec> final      = {lts.decaytree[0].p4, lts.decaytree[1].p4};
  M4Vec              k1, k2;
  M4Vec              axis = lts.q1_in_X;
  M4Vec              d0   = lts.d0_in_X;
  M4Vec              d1   = lts.d1_in_X;
  if (transverse) {
    if (!mg5helas::PrepareOnShellKinematics(lts, final, k1, k2)) { return; }
    const M4Vec  hard = final[0] + final[1];
    const double mass = hard.M();
    axis              = k1;
    d0                = final[0];
    d1                = final[1];
    kinematics::LorentzBoost(hard, mass, axis, -1);
    kinematics::LorentzBoost(hard, mass, d0, -1);
    kinematics::LorentzBoost(hard, mass, d1, -1);
  }

  const double shat        = transverse ? (k1 + k2).M2() : lts.s_hat;
  const double q1_mod      = axis.P3mod();
  const double d0_mod      = d0.P3mod();
  const double angle_denom = q1_mod * d0_mod;
  if (!std::isfinite(angle_denom) || angle_denom <= 0.0) {
    Amp[0].fill(0.0);
    return;
  }
  const gra::M4Vec d0_in_hard_frame = DurhamRotateToAxisFrame(d0, axis);
  const double     costheta_raw     = d0_in_hard_frame.Pz() / d0_mod;
  if (!std::isfinite(costheta_raw)) {
    Amp[0].fill(0.0);
    return;
  }
  const double costheta  = std::clamp(costheta_raw, -1.0, 1.0);
  const double costheta2 = pow2(costheta);

  // ------------------------------------------------------------------
  // ** Hard angular cut-off **
  // some sub-amplitudes are singular when |costheta| -> 1

  if (std::abs(costheta) > param.MAXCOS) {
    Amp[0].fill(0.0);  // Compute zero
    return;
  }

  // Keep the leading-twist meson kernel inside its perturbative hard-pT domain
  const gra::M4Vec d1_in_hard_frame = DurhamRotateToAxisFrame(d1, axis);
  const double     hard_pt0         = d0_in_hard_frame.Pt();
  const double     hard_pt1         = d1_in_hard_frame.Pt();
  if (!std::isfinite(hard_pt0) || !std::isfinite(hard_pt1)) {
    Amp[0].fill(0.0);
    return;
  }
  if (hard_pt0 < param.MESON_pt_min || hard_pt1 < param.MESON_pt_min) {
    Amp[0].fill(0.0);
    return;
  }
  // ------------------------------------------------------------------

  // Restore the beam-axis azimuth for JzP after rotating into the hard frame
  // Spinless two-body transport: M_+- -> e^{+2i phi} M_+-
  const double               phi      = d0_in_hard_frame.Phi() + axis.Phi();
  const std::complex<double> posphase = std::exp(2.0 * zi * phi);
  const std::complex<double> negphase = std::exp(-2.0 * zi * phi);

  // ------------------------------------------------------------------
  // Mesons scalar flavor octet (non-singlet) amplitude:
  // |\pi0\pi0>, |\pi+\pi->, |K+K->, |K0\bar{K0}>
  //
  // T_+- = T_-+
  auto T_SFO_PM = [&](double x, double y) {
    const double a = (1.0 - x) * (1.0 - y) + x * y;  // +
    const double b = (1.0 - x) * (1.0 - y) - x * y;  // -

    return 1.0 / (x * y * (1.0 - x) * (1.0 - y)) * (x * (1.0 - x) + y * (1.0 - y)) / (pow2(a) - pow2(b) * costheta2) *
           (NC / 2.0) * (costheta2 - 2.0 * CF / NC * a);
  };
  // ------------------------------------------------------------------

  // ------------------------------------------------------------------
  // SU(3)_F scalar flavor-singlet hard kernel for eta0 eta0 components:
  //
  // Three normalized quark flavours contribute coherently to each ladder amplitude
  // [REFERENCE: Harland-Lang et al., arXiv:1105.1626v2, Eqs. (4.1)-(4.6)]
  constexpr double singlet_flavour_multiplicity = 3.0;

  // T_++ = T_--
  auto T_SFS_PP = [&](double x, double y) {
    return singlet_flavour_multiplicity / (x * y * (1.0 - x) * (1.0 - y)) * (1.0 + costheta2) / pow2(1.0 - costheta2);
  };

  // T_+- = T_-+
  auto T_SFS_PM = [&](double x, double y) {
    return singlet_flavour_multiplicity / (x * y * (1.0 - x) * (1.0 - y)) * (1.0 + 3.0 * costheta2) /
           (2.0 * pow2(1.0 - costheta2));
  };
  // ------------------------------------------------------------------

  const int raw_pdg0 = lts.decaytree[0].p.pdg;
  const int raw_pdg1 = lts.decaytree[1].p.pdg;
  const int pdg0     = std::abs(raw_pdg0);
  const int pdg1     = std::abs(raw_pdg1);

  const auto is_sfo_pdg = [](int code) { return code == 111 || code == 211 || code == 321 || code == 311; };
  const auto is_sfs_pdg = [](int code) { return code == 221 || code == 331; };

  const bool is_sfo  = is_sfo_pdg(pdg0) && is_sfo_pdg(pdg1);
  const bool is_sfs  = is_sfs_pdg(pdg0) && is_sfs_pdg(pdg1);
  // ------------------------------------------------------------------
  // Integral over meson wave functions:
  // M\int_0^1 dx dy \phi_M(x) \phi_\bar{M}(y) T_{\lambda\lambda'}
  // (x,y,\hat{s},\theta)

  // Normalization factor after normalized incoming singlet contraction:
  // (delta^AB / (N_C^2-1)) * (delta^AB / N_C) = 1/N_C
  const double norm = (1.0 / NC) * (64.0 * math::PIPI * pow2(alpha_s)) / shat;
  double       A0   = 0.0;
  double       A2   = 0.0;

  // Scalar flavor octet
  if (is_sfo) {
    // Evaluate or reuse scalar-octet meson wave functions on the active grid
    const std::vector<double> &wfphi0 = CachedMesonWaveFunction(xval, pdg0);
    const std::vector<double> &wfphi1 = CachedMesonWaveFunction(xval, pdg1);

    MMatrix<double> f_PM(Nx, Nx, 0.0);

    for (std::size_t i = 0; i < Nx; ++i) {
      for (std::size_t j = 0; j < Nx; ++j) { f_PM[i][j] = wfphi0[i] * wfphi1[j] * T_SFO_PM(xval[i], xval[j]); }
    }

    A2 = norm * TensorProductIntegral(f_PM, xweight);
  }

  // Physical eta/eta-prime quark-Fock pairs: ordinary OO and SS plus singlet
  // ladder SS
  else if (is_sfs) {
    MMatrix<double> f_PP(Nx, Nx, 0.0);
    MMatrix<double> f_PM(Nx, Nx, 0.0);

    const DurhamEtaMixing mix0    = DurhamEtaMixingCoefficients(pdg0, param.eta_theta8_deg, param.eta_theta1_deg);
    const DurhamEtaMixing mix1    = DurhamEtaMixingCoefficients(pdg1, param.eta_theta8_deg, param.eta_theta1_deg);
    const double          fM_eta8 = param.f_pi * param.f_eta8_over_fpi;
    const double          fM_eta0 = param.f_pi * param.f_eta0_over_fpi;
    std::vector<double>   wf8_0(Nx, 0.0);
    std::vector<double>   wf8_1(Nx, 0.0);
    std::vector<double>   wf0_0(Nx, 0.0);
    std::vector<double>   wf0_1(Nx, 0.0);

    for (std::size_t i = 0; i < Nx; ++i) {
      wf8_0[i] = mix0.octet * phi_CZ(xval[i], fM_eta8);
      wf8_1[i] = mix1.octet * phi_CZ(xval[i], fM_eta8);
      wf0_0[i] = mix0.singlet * phi_CZ(xval[i], fM_eta0);
      wf0_1[i] = mix1.singlet * phi_CZ(xval[i], fM_eta0);
    }

    for (std::size_t i = 0; i < Nx; ++i) {
      for (std::size_t j = 0; j < Nx; ++j) {
        const double octet_octet     = wf8_0[i] * wf8_1[j];
        const double singlet_singlet = wf0_0[i] * wf0_1[j];
        f_PP[i][j]                   = singlet_singlet * T_SFS_PP(xval[i], xval[j]);
        f_PM[i][j] =
            (octet_octet + singlet_singlet) * T_SFO_PM(xval[i], xval[j]) + singlet_singlet * T_SFS_PM(xval[i], xval[j]);
      }
    }

    A0 = norm * TensorProductIntegral(f_PP, xweight);
    A2 = norm * TensorProductIntegral(f_PM, xweight);
  }

  if (transverse) {
    Amp[0] = DurhamMesonTensor(k1, k2, final[0], A0, A2);
  } else {
    Amp[0][DurhamHelicityIndex(1, 1)]   = A0;
    Amp[0][DurhamHelicityIndex(-1, -1)] = A0;
    Amp[0][DurhamHelicityIndex(1, -1)]  = A2 * posphase;
    Amp[0][DurhamHelicityIndex(-1, 1)]  = A2 * negphase;
  }
}

}  // namespace gra
