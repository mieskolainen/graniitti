// Analytic quark box amplitudes for gg -> gamma gamma
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Amplitude/Photon/AMP_gg_yy.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <numbers>
#include <stdexcept>

#include "Graniitti/Amplitude/MG5/Runtime/AMP_MG5_Kinematics.h"
#include "Graniitti/Amplitude/MG5/Runtime/Models/sm/HelAmps_sm_lepton_masses.h"
#include "Graniitti/MModelCache.h"
#include "Graniitti/Photon/MQED.h"
#include "Graniitti/Tech/MAux.h"

namespace gra {
using gra::aux::indices;
namespace {
using Complex = std::complex<long double>;
using Helicity5 = std::array<Complex, 5>;
constexpr long double pi = std::numbers::pi_v<long double>;

// Evaluate the rapidly convergent Spence series after analytic mapping
// Li2(x) = z-z^2/4+sum_(n>=1) B_(2n) z^(2n+1)/(2n+1)!, z=-log(1-x)
Complex SpenceSeries(const Complex &x) {
  const Complex z   = -std::log(Complex{1.0L, 0.0L} - x);
  const Complex z2  = z * z;
  Complex       sum = 43867.0L / 798.0L;
  sum               = sum * z2 / 342.0L - 3617.0L / 510.0L;
  sum               = sum * z2 / 272.0L + 7.0L / 6.0L;
  sum               = sum * z2 / 210.0L - 691.0L / 2730.0L;
  sum               = sum * z2 / 156.0L + 5.0L / 66.0L;
  sum               = sum * z2 / 110.0L - 1.0L / 30.0L;
  sum               = sum * z2 / 72.0L + 1.0L / 42.0L;
  sum               = sum * z2 / 42.0L - 1.0L / 30.0L;
  sum               = sum * z2 / 20.0L + 1.0L / 6.0L;
  return z2 * z * sum / 6.0L - 0.25L * z2 + z;
}

// Evaluate the principal branch of the complex dilogarithm
// Use Li2 reflection and inversion identities outside the series domain
Complex Li2(const Complex &x) {
  const Complex one{1.0L, 0.0L};
  if (std::abs(one - x) < 1.0e-13L) { return {(pi * pi / 6.0L), 0.0L}; }
  if (std::abs(one - x) < 0.5L) { return (pi * pi / 6.0L) - std::log(one - x) * std::log(x) - SpenceSeries(one - x); }
  if (std::abs(x) > 1.0L) { return -(pi * pi / 6.0L) - 0.5L * std::log(-x) * std::log(-x) - SpenceSeries(one / x); }
  return SpenceSeries(x);
}

// Compute the finite equal-mass two-point scalar function
// B0(q,m2) = -beta log[(beta+1)^2/(-4m2/q)]
Complex B0(const long double q, const Complex &m2) {
  if (std::abs(q) < 0.5L * m2.real()) {
    const long double x = q / m2.real();
    long double term = x / 6.0L;
    long double sum = -2.0L;
    for (unsigned int n = 1; n <= 24; ++n) {
      sum += term;
      term *= x * n / (2.0L * (2 * n + 3));
    }
    return sum;
  }
  const Complex beta = std::sqrt(1.0L - 4.0L * m2 / q);
  return -beta * std::log((beta + 1.0L) * (beta + 1.0L) / (-4.0L * m2 / q));
}

// Compute the equal-mass triangle with two real photon legs
Complex C0(const long double q, const Complex &m2) {
  if (std::abs(q) < 0.5L * m2.real()) {
    const long double x = q / m2.real();
    long double term = -0.5L;
    long double sum = 0.0L;
    for (unsigned int n = 0; n <= 24; ++n) {
      sum += term;
      term *= x * (n + 1) * (n + 1) / ((2.0L * n + 3) * (2.0L * n + 4));
    }
    return sum / m2.real();
  }
  const long double  q2   = -q;
  const Complex root = std::sqrt(1.0L + 4.0L * m2 / q2);
  const Complex x    = (1.0L + root) / 2.0L;
  const Complex c01  = -(Li2(1.0L / x) + Li2((-q2 / m2) * x)) / q2;
  return -c01;
}

// Compute the equal-mass box with four real photon legs
Complex D0(const long double q, const long double p, const Complex &m2) {
  if (std::max(std::abs(q), std::abs(p)) < 0.5L * m2.real()) {
    // Expand the Feynman parameter denominator before cancellations of heavy quark terms
    // D0 = sum_(k,l>=0) (k+l+1)! k! l! x^k y^l / [(2k+2l+3)! m^4]
    const long double x = q / m2.real();
    const long double y = p / m2.real();
    std::array<long double, 25> xp;
    std::array<long double, 25> yp;
    xp[0] = yp[0] = 1.0L;
    for (std::size_t n = 1; n < xp.size(); ++n) {
      xp[n] = xp[n - 1] * x;
      yp[n] = yp[n - 1] * y;
    }
    long double coefficient = 1.0L / 6.0L;
    long double sum = 0.0L;
    for (const auto &n : indices(xp)) {
      long double weight = coefficient;
      long double bound = 0.0L;
      for (std::size_t k = 0; k <= n; ++k) {
        const long double term = weight * xp[k] * yp[n - k];
        sum += term;
        bound += std::abs(term);
        if (k < n) { weight *= static_cast<long double>(k + 1) / (n - k); }
      }
      if (n > 1 && bound < std::numeric_limits<long double>::epsilon() * std::abs(sum)) { break; }
      coefficient *= (n + 1) / (2.0L * (2 * n + 5));
    }
    return sum / (m2.real() * m2.real());
  }
  const Complex a1    = std::sqrt(1.0L - 4.0L * m2 / q);
  const Complex a2    = std::sqrt(1.0L - 4.0L * m2 / p);
  const Complex a3    = std::sqrt(1.0L - 4.0L * m2.real() * (q + p) / (q * p));
  const Complex la3   = 4.0L * m2 * (q + p) / (q * p) / (1.0L + a3);
  const Complex a2a3  = 4.0L * m2 / q / (a2 + a3);
  const Complex a1a3  = 4.0L * m2 / p / (a1 + a3);
  const Complex one3  = 1.0L + a3;
  Complex       value = Li2(one3 / (a1 + a3)) - Li2(la3 / a1a3) - Li2(-a1a3 / one3) -
                  0.5L * std::log(a1a3 / one3) * std::log(a1a3 / one3) - (pi * pi / 6.0L) - Li2(-la3 / (a1 + a3)) +
                  Li2(one3 / (a2 + a3)) - Li2(la3 / a2a3) - Li2(-a2a3 / one3) -
                  0.5L * std::log(a2a3 / one3) * std::log(a2a3 / one3) - (pi * pi / 6.0L) - Li2(-la3 / (a2 + a3));
  if ((a1 + a3).imag() < 0.0L) { value += Complex{0.0L, 2.0L * pi} * std::log(1.0L + 1.0L / a3); }
  value *= -2.0L / (q * p) / a3;
  if (q * p > 0.0L || 4.0L * m2.real() > q) { value = {value.real(), 0.0L}; }
  return value;
}

// Add one massive quark box to the five independent amplitudes
void AddFermion(const long double s, const long double t, const long double u, const long double mass, const long double charge2,
                Helicity5 &amplitude) {
  const long double  mr2 = mass * mass;
  const Complex m2{mr2, -1.0e-30L};
  const Complex bs   = B0(s, m2);
  const Complex bt   = B0(t, m2);
  const Complex bu   = B0(u, m2);
  const Complex cs   = C0(s, m2);
  const Complex ct   = C0(t, m2);
  const Complex cu   = C0(u, m2);
  const Complex dst  = D0(s, t, m2);
  const Complex dsu  = D0(s, u, m2);
  const Complex dut  = D0(u, t, m2);
  const Complex dsum = dst + dsu + dut;
  const long double  r2   = s * s + t * t + u * u;

  amplitude[0] +=
      charge2 * (-1.0L + (u - t) / s * (bu - bt) + (4.0L * mr2 / s + 2.0L * (t * u / (s * s) - 0.5L)) * (u * cu + t * ct) -
                 2.0L * mr2 * s * (mr2 / s - 0.5L) * dsum - t * u * (4.0L * mr2 / s + t * u / (s * s) - 0.5L) * dut);
  amplitude[1] +=
      charge2 * (1.0L - mr2 * r2 * (cs / (u * t) + ct / (s * u) + cu / (s * t)) -
                 mr2 * ((2.0L * mr2 + s * t / u) * dst + (2.0L * mr2 + s * u / t) * dsu + (2.0L * mr2 + t * u / s) * dut));
  amplitude[2] += charge2 * (1.0L - 2.0L * mr2 * mr2 * dsum);
  amplitude[3] +=
      charge2 * (-1.0L + (u - s) / t * (bu - bs) + (4.0L * mr2 / t + 2.0L * (s * u / (t * t) - 0.5L)) * (u * cu + s * cs) -
                 2.0L * mr2 * t * (mr2 / t - 0.5L) * dsum - s * u * (4.0L * mr2 / t + s * u / (t * t) - 0.5L) * dsu);
  amplitude[4] +=
      charge2 * (-1.0L + (s - t) / u * (bs - bt) + (4.0L * mr2 / u + 2.0L * (t * s / (u * u) - 0.5L)) * (s * cs + t * ct) -
                 2.0L * mr2 * u * (mr2 / u - 0.5L) * dsum - t * s * (4.0L * mr2 / u + t * s / (u * u) - 0.5L) * dst);
}

// Compute the HELAS phase relative to a covariantly transported planar polarization
std::complex<double> PolarizationPhase(const M4Vec &p, int helicity, int direction,
                                      const M4Vec &real, const M4Vec &imaginary) {
  auto momentum = p.Contravariant();
  std::array<std::complex<double>, 6> wave;
  MG5_sm_lepton_masses::vxxxxx(momentum.data(), 0.0, helicity, direction, wave.data());
  const std::array<std::complex<double>, 4> epsilon = {wave[2], wave[3], wave[4], wave[5]};
  return -MinkowskiProduct(real.Contravariant(), epsilon) +
         math::zi * MinkowskiProduct(imaginary.Contravariant(), epsilon);
}

// Construct the scattering-plane tetrad and transport every external helicity phase
std::array<std::array<std::complex<double>, 2>, 4> BoxPhases(const std::array<M4Vec, 4> &p,
                                                           double s, double t, double u) {
  const double mass = std::sqrt(s);
  const double cosine = (t - u) / s;
  const double sine = 2.0 * std::sqrt(t * u) / s;
  const M4Vec hard = p[0] + p[1];
  const M4Vec time = hard / mass;
  const M4Vec z = (p[0] - p[1]) / mass;
  const M4Vec x = (p[2] * (2.0 / mass) - time - z * cosine) / sine;
  M4Vec rest_z = z;
  M4Vec rest_x = x;
  kinematics::LorentzBoost(hard, mass, rest_z, -1);
  kinematics::LorentzBoost(hard, mass, rest_x, -1);
  const auto normal = CrossProduct(std::array<double, 3>{rest_z.Px(), rest_z.Py(), rest_z.Pz()},
                                   std::array<double, 3>{rest_x.Px(), rest_x.Py(), rest_x.Pz()});
  M4Vec y(normal[0], normal[1], normal[2], 0.0);
  kinematics::LorentzBoost(hard, mass, y, 1);
  const M4Vec transverse = x * cosine - z * sine;
  const std::array<M4Vec, 4> real_axis = {x, x, transverse, transverse};
  const std::array<double, 4> imaginary_sign = {-1.0, 1.0, 1.0, -1.0};
  std::array<std::array<std::complex<double>, 2>, 4> phase;
  for (const auto &leg : indices(phase)) {
    for (const auto &hel : indices(phase[leg])) {
      const int h = hel == 0 ? -1 : 1;
      phase[leg][hel] = PolarizationPhase(p[leg], h, leg < 2 ? -1 : 1,
                                          real_axis[leg] * (-h / std::sqrt(2.0)),
                                          y * (imaginary_sign[leg] / std::sqrt(2.0)));
    }
  }
  return phase;
}

}  // namespace

// Define the native quark-box process independently of MG2GRA
amplitude::Process AMP_gg_yy::Definition() {
  const AmplitudeTopology topology = {AmplitudeTopologyNode{{22}, {}}, AmplitudeTopologyNode{{22}, {}}};
  return amplitude::Process{"DURHAM", "gg_yy", "gamma gamma", "gamma gamma", "a,a", "g g > a a",
                            topology, {topology}, {22, 22},
                            DecayStructure{DecayType::Full, true},
                            amplitude::MatrixElementForm::Analytic, amplitude::TopologyMode::Exact};
}

// Prepare all six massive quark loops from the validated Standard Model inputs
AMP_gg_yy::AMP_gg_yy(const SMParam &sm)
    : DurhamMG5Process(std::vector<amplitude::Process>{Definition()}),
      quarks_{{{sm.d, 1.0 / 9.0}, {sm.u, 4.0 / 9.0}, {sm.s, 1.0 / 9.0},
               {sm.c, 4.0 / 9.0}, {sm.b, 1.0 / 9.0}, {sm.t, 4.0 / 9.0}}},
      alpha_(1.0 / sm.alpha_em_inv) {}

// Compute A_ab = delta_ab 4 alpha_s alpha sum_q Q_q^2 H_q
// Tr(T_a T_b) = delta_ab/2 replaces the QED closed-color trace
// [REFERENCE: Bardin et al., arXiv:0911.5634]
std::array<std::complex<double>, 16> AMP_gg_yy::Helicity(double s, double t, double u,
                                                       double alpha_s, double alpha) const {
  Helicity5 independent = {};
  for (const auto &quark : quarks_) { AddFermion(s, t, u, quark.mass, quark.charge2, independent); }
  constexpr std::array<std::size_t, 16> map = {0, 1, 1, 2, 1, 4, 3, 1, 1, 3, 4, 1, 2, 1, 1, 0};
  const long double coupling = 4.0L * alpha_s * alpha;
  std::array<std::complex<double>, 16> amplitude;
  for (const auto &row : indices(amplitude)) {
    amplitude[row] = static_cast<std::complex<double>>(coupling * independent[map[row]]);
  }
  return amplitude;
}

// Evaluate the quark box and its exact normalized incoming gluon singlet projection
DurhamMG5Evaluation AMP_gg_yy::Evaluate(LORENTZSCALAR &lts, double alpha_s, M4Vec *hard_k1, M4Vec *hard_k2) {
  DurhamMG5Evaluation result;
  lts.hamp.clear();
  if (hard_k1 != nullptr) { *hard_k1 = M4Vec(); }
  if (hard_k2 != nullptr) { *hard_k2 = M4Vec(); }
  const auto final = mg5::StableLeafMomenta(mg5::StableDecayLeaves(lts.decaytree));
  std::array<M4Vec, 4> p;
  if (final.size() != 2 || !mg5helas::PrepareOnShellKinematics(lts, final, p[0], p[1])) {
    result.status = mg5helas::EvaluationStatus::KinematicsFailure;
    return result;
  }
  p[2] = final[0];
  p[3] = final[1];
  const double s = (p[0] + p[1]).M2();
  const double t = -2.0 * MinkowskiProduct(p[0].Contravariant<long double>(), p[2].Contravariant<long double>());
  const double u = -s - t;
  if (!(s > 0.0) || !(t < 0.0) || !(u < 0.0)) {
    result.status = mg5helas::EvaluationStatus::KinematicsFailure;
    return result;
  }
  const double alpha = qed::alpha_QED(s, lts.model_cache->Tune().Structure().QED_alpha, alpha_);
  const auto hard = Helicity(s, t, u, alpha_s, alpha);
  const auto phase = BoxPhases(p, s, t, u);
  result.projected.resize(16);
  for (const auto &row : indices(result.projected)) {
    std::complex<double> value = std::sqrt(8.0) * hard[row];
    for (const auto &leg : indices(phase)) { value *= phase[leg][(row >> (3 - leg)) & 1U]; }
    result.projected[row] = value;
  }
  if (hard_k1 != nullptr) { *hard_k1 = p[0]; }
  if (hard_k2 != nullptr) { *hard_k2 = p[1]; }
  return result;
}

}  // namespace gra
