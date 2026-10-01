// Coherent and incoherent nuclear transverse photon densities
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Nuclear/MPhoton.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

#include "Graniitti/Math/MMath.h"
#include "Graniitti/Particle/MForm.h"
#include "Graniitti/Particle/MPDG.h"
#include "Graniitti/Photon/MQED.h"

namespace gra::nuclear {
namespace {

// Compute whether one photon momentum fraction is physical
bool ValidXi(const double xi) { return std::isfinite(xi) && xi > 0.0 && xi < 1.0; }

}  // namespace

// Compute the unpolarized photon density
// n = n_parallel + n_perpendicular
double PhotonDensity::Trace() const { return parallel + perpendicular; }

// Compute one transverse nuclear photon-source helicity amplitude
// A_lambda = sign(Q_beam) a_parallel exp(i lambda phi) / sqrt(2 xi)
std::complex<double> PhotonAmp(const double source, const double xi, const double charge, const double phi,
                               const int helicity) {
  if (!ValidXi(xi)) { return 0.0; }
  const double               charge_phase = charge / std::abs(charge);
  const std::complex<double> phase        = std::exp(math::zi * static_cast<double>(helicity) * phi);
  return charge_phase * source * phase / std::sqrt(2.0 * xi);
}

// Compute the point charge EPA kernel with its small-x limit and floating-point tail
double PhotonKernel(const double x, const double gamma) {
  if (!(x > std::numeric_limits<double>::epsilon())) { return 1.0; }
  if (x > -0.5 * std::log(std::numeric_limits<double>::min())) { return 0.0; }
  const double xk0 = x * std::cyl_bessel_k(0.0, x), xk1 = x * std::cyl_bessel_k(1.0, x);
  return xk1 * xk1 + xk0 * xk0 / (gamma * gamma);
}

// Compute n = Z^2 alpha C(b)^2 K(omega b/(gamma hbarc))/(pi^2 omega b^2)
// [REFERENCE: Baltz et al., Phys. Rept. 458 (2008) 1, impact dependent equivalent photon density]
double ImpactPhotonDensity(double omega, double b, double gamma, double charge, double fraction) {
  if (std::fpclassify(b) == FP_ZERO) { return fraction > 0.0 ? std::numeric_limits<double>::infinity() : 0.0; }
  return qed::alpha_QED() * math::pow2(charge * fraction / b) * PhotonKernel(omega * b / (gamma * PDG::GeV2fm), gamma) /
         (math::PIPI * omega);
}

// Construct one immutable nuclear photon source
MPhoton::MPhoton(const MNucleus& nucleus, form::ParamStore structure)
    : MPhoton(std::make_shared<const MNucleus>(nucleus), std::move(structure)) {}

// Construct one photon source sharing an immutable nuclear model
MPhoton::MPhoton(std::shared_ptr<const MNucleus> nucleus, form::ParamStore structure)
    : nucleus_(std::move(nucleus)), structure_(std::move(structure)) {
  if (nucleus_ == nullptr) { throw std::invalid_argument("MPhoton: null nuclear model"); }
}

// Compute the signed coherent transverse photon source
// a_parallel = 4 sqrt[pi alpha (1-xi)/xi] p_T Z F_A(Q)/(p_T^2 + xi^2 M_A^2), Q^2 = -t
// The real diffraction sign belongs to F_A(Q) before any density is formed
// [REFERENCE: Harland-Lang et al., Eur. Phys. J. C 79 (2019) 39, Eqs. (11)-(14)]
double MPhoton::CoherentAmp(const double xi, const double t, const double pt) const {
  if (!ValidXi(xi) || !std::isfinite(t) || t > 0.0 || !std::isfinite(pt) || !(pt > 0.0)) { return 0.0; }
  const double pt2    = pt * pt;
  const double mass2  = nucleus_->Mass() * nucleus_->Mass();
  const double denom  = pt2 + xi * xi * mass2;
  const double form   = nucleus_->ChargeDensity().Form(std::sqrt(-t));
  const double norm   = 4.0 * std::sqrt(math::PI * qed::alpha_QED() * (1.0 - xi) / xi);
  const double source = norm * (pt / denom) * static_cast<double>(nucleus_->Z()) * form;
  return std::isfinite(source) ? source : 0.0;
}

// Compute the elastic bound-proton transverse density
// n_E = C (1 - xi) p_T^2 F_E / D^2, n_M = C xi^2 G_M^2 / (4 D)
// C = 16 pi alpha / xi, D = p_T^2 + xi^2 m_p^2, F_E = (4 m_p^2 G_E^2 + Q^2 G_M^2) / (4 m_p^2 + Q^2)
// n_parallel = n_E + n_M, n_perpendicular = n_M
// [REFERENCE: Budnev et al., Phys. Rept. 15 (1975) 181]
PhotonDensity MPhoton::Proton(const double xi, const double t, const double pt) const {
  if (!ValidXi(xi) || !std::isfinite(t) || t > 0.0 || !std::isfinite(pt) || pt < 0.0) { return {}; }
  const double pt2           = pt * pt;
  const double q2            = -t;
  const double ge            = form::G_E(q2, structure_);
  const double gm            = form::G_M(q2, structure_);
  const double mass2         = PDG::mp * PDG::mp;
  const double electric_part = (4.0 * mass2 * ge * ge + q2 * gm * gm) / (4.0 * mass2 + q2);
  const double denom         = pt2 + xi * xi * mass2;
  const double common        = 16.0 * math::PI * qed::alpha_QED() / xi;
  const double electric      = common * (1.0 - xi) * math::pow2(pt / denom) * electric_part;
  const double magnetic      = 4.0 * math::PI * qed::alpha_QED() * xi * gm * gm / denom;
  if (!std::isfinite(electric) || !std::isfinite(magnetic) || electric < 0.0 || magnetic < 0.0) { return {}; }
  return {electric + magnetic, magnetic};
}

// Compute the charge-current variance in proton-count units
// Var[J] = <|J|^2> - |<J>|^2 = Z [1 - F_A(q)^2]
// [REFERENCE: Good and Walker, Phys. Rev. 120 (1960) 1857]
double MPhoton::Variance(const M3Vec& q, const MConfigBank* bank) const {
  if (!std::isfinite(q[0]) || !std::isfinite(q[1]) || !std::isfinite(q[2])) {
    throw std::invalid_argument("MPhoton::Variance: non-finite rest-frame transfer");
  }
  if (bank != nullptr) {
    if (bank->Nucleus().ID().pdg != nucleus_->ID().pdg) {
      throw std::invalid_argument("MPhoton::Variance: configuration bank nucleus does not match");
    }
    return bank->ChargeStat(q[0], q[1], q[2]).variance;
  }
  const double form = nucleus_->ChargeDensity().Form(std::hypot(q[0], q[1], q[2]));
  return static_cast<double>(nucleus_->Z()) * std::max(0.0, 1.0 - form * form);
}

// Compute the incoherent bound-proton fluctuation density
// n_inc(xi,t) = A Var[J(q)] n_p(A xi,t_p)
// [REFERENCE: Budnev et al., Phys. Rept. 15 (1975) 181]
PhotonDensity MPhoton::Incoherent(const double xi, const double t, const double pt, const M3Vec& q,
                                  const MConfigBank* bank) const {
  if (!ValidXi(xi) || !std::isfinite(t) || t > 0.0 || !std::isfinite(pt) || pt < 0.0) { return {}; }
  const double proton_xi = static_cast<double>(nucleus_->A()) * xi;
  if (!ValidXi(proton_xi)) { return {}; }
  const double  proton_q2 = (pt * pt + proton_xi * proton_xi * PDG::mp * PDG::mp) / (1.0 - proton_xi);
  PhotonDensity density   = Proton(proton_xi, -proton_q2, pt);
  const double  scale     = static_cast<double>(nucleus_->A()) * Variance(q, bank);
  density.parallel *= scale;
  density.perpendicular *= scale;
  if (!std::isfinite(density.Trace()) || density.parallel < 0.0 || density.perpendicular < 0.0) { return {}; }
  return density;
}

// Compute a smooth nuclear density from the invariant momentum transfer
PhotonDensity MPhoton::Density(const CoherenceType type, const double xi, const double t, const double pt) const {
  return Density(type, xi, t, pt, M3Vec{std::sqrt(std::max(0.0, -t)), 0.0, 0.0}, nullptr);
}

// Compute one coherent, incoherent or inclusive nuclear photon density
// n_coh = (|a_parallel|^2,0), n_incl = n_coh + n_inc
// [REFERENCE: Budnev et al., Phys. Rept. 15 (1975) 181]
PhotonDensity MPhoton::Density(const CoherenceType type, const double xi, const double t, const double pt,
                               const M3Vec& q, const MConfigBank* bank) const {
  if (type == CoherenceType::Incoherent) { return Incoherent(xi, t, pt, q, bank); }
  const double  source = CoherentAmp(xi, t, pt);
  PhotonDensity density{source * source, 0.0};
  if (!std::isfinite(density.parallel)) { return {}; }
  if (type == CoherenceType::Coherent) { return density; }
  if (type != CoherenceType::Inclusive) { throw std::invalid_argument("MPhoton::Density: invalid sector"); }
  const auto incoherent = Incoherent(xi, t, pt, q, bank);
  density.parallel += incoherent.parallel;
  density.perpendicular = incoherent.perpendicular;
  return density;
}

}  // namespace gra::nuclear
