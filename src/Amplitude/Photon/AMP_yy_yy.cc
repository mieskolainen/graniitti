// Standard Model yy -> yy one loop amplitude
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.
//
// [REFERENCE: Bardin et al., arXiv:0911.5634]
// [REFERENCE: S. Navas et al., Phys. Rev. D 110 (2024) 030001, 2025 update]

#include "Graniitti/Amplitude/Photon/AMP_yy_yy.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <utility>

#include "Graniitti/Amplitude/MG5/Runtime/AMP_MG5_Helicity.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Particle/MForm.h"
#include "Graniitti/Particle/MPDG.h"
#include "Graniitti/Photon/MFlux.h"
#include "Graniitti/Photon/MQED.h"
#include "Graniitti/Tech/MAux.h"

namespace gra::lbyl {
using gra::aux::indices;
namespace {

using Complex   = std::complex<double>;
using Helicity5 = std::array<Complex, 5>;

// Evaluate the rapidly convergent Spence series after analytic mapping
// Li2(x) = z-z^2/4+sum_(n>=1) B_(2n) z^(2n+1)/(2n+1)!, z=-log(1-x)
Complex SpenceSeries(const Complex &x) {
  const Complex z   = -std::log(Complex{1.0, 0.0} - x);
  const Complex z2  = z * z;
  Complex       sum = 43867.0 / 798.0;
  sum               = sum * z2 / 342.0 - 3617.0 / 510.0;
  sum               = sum * z2 / 272.0 + 7.0 / 6.0;
  sum               = sum * z2 / 210.0 - 691.0 / 2730.0;
  sum               = sum * z2 / 156.0 + 5.0 / 66.0;
  sum               = sum * z2 / 110.0 - 1.0 / 30.0;
  sum               = sum * z2 / 72.0 + 1.0 / 42.0;
  sum               = sum * z2 / 42.0 - 1.0 / 30.0;
  sum               = sum * z2 / 20.0 + 1.0 / 6.0;
  return z2 * z * sum / 6.0 - 0.25 * z2 + z;
}

// Evaluate the principal branch of the complex dilogarithm
// Use Li2 reflection and inversion identities outside the series domain
Complex Li2(const Complex &x) {
  const Complex one{1.0, 0.0};
  if (std::abs(one - x) < 1.0e-13) { return {gra::math::ZETA2, 0.0}; }
  if (std::abs(one - x) < 0.5) { return gra::math::ZETA2 - std::log(one - x) * std::log(x) - SpenceSeries(one - x); }
  if (std::abs(x) > 1.0) { return -gra::math::ZETA2 - 0.5 * std::log(-x) * std::log(-x) - SpenceSeries(one / x); }
  return SpenceSeries(x);
}

// Compute the finite equal-mass two-point scalar function
// B0(q,m2) = -beta log[(beta+1)^2/(-4m2/q)]
Complex B0(const double q, const Complex &m2) {
  const Complex beta = std::sqrt(1.0 - 4.0 * m2 / q);
  return -beta * std::log((beta + 1.0) * (beta + 1.0) / (-4.0 * m2 / q));
}

// Compute the equal-mass triangle with two real photon legs
Complex C0(const double q, const Complex &m2) {
  const double  q2   = -q;
  const Complex root = std::sqrt(1.0 + 4.0 * m2 / q2);
  const Complex x    = (1.0 + root) / 2.0;
  const Complex c01  = -(Li2(1.0 / x) + Li2((-q2 / m2) * x)) / q2;
  return -c01;
}

// Compute the equal-mass box with four real photon legs
Complex D0(const double q, const double p, const Complex &m2) {
  const Complex a1    = std::sqrt(1.0 - 4.0 * m2 / q);
  const Complex a2    = std::sqrt(1.0 - 4.0 * m2 / p);
  const Complex a3    = std::sqrt(1.0 - 4.0 * m2.real() * (q + p) / (q * p));
  const Complex la3   = 4.0 * m2 * (q + p) / (q * p) / (1.0 + a3);
  const Complex a2a3  = 4.0 * m2 / q / (a2 + a3);
  const Complex a1a3  = 4.0 * m2 / p / (a1 + a3);
  const Complex one3  = 1.0 + a3;
  Complex       value = Li2(one3 / (a1 + a3)) - Li2(la3 / a1a3) - Li2(-a1a3 / one3) -
                  0.5 * std::log(a1a3 / one3) * std::log(a1a3 / one3) - gra::math::ZETA2 - Li2(-la3 / (a1 + a3)) +
                  Li2(one3 / (a2 + a3)) - Li2(la3 / a2a3) - Li2(-a2a3 / one3) -
                  0.5 * std::log(a2a3 / one3) * std::log(a2a3 / one3) - gra::math::ZETA2 - Li2(-la3 / (a2 + a3));
  if ((a1 + a3).imag() < 0.0) { value += Complex{0.0, 2.0 * gra::math::PI} * std::log(1.0 + 1.0 / a3); }
  value *= -2.0 / (q * p) / a3;
  if (q * p > 0.0 || 4.0 * m2.real() > q) { value = {value.real(), 0.0}; }
  return value;
}

// Add one charged fermion loop to the five independent amplitudes
void AddFermion(const double s, const double t, const double u, const double mass, const double charge4,
                Helicity5 &amplitude) {
  const double  mr2 = mass * mass;
  const Complex m2{mr2, -1.0e-30};
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
  const double  r2   = s * s + t * t + u * u;

  amplitude[0] +=
      charge4 * (-1.0 + (u - t) / s * (bu - bt) + (4.0 * mr2 / s + 2.0 * (t * u / (s * s) - 0.5)) * (u * cu + t * ct) -
                 2.0 * mr2 * s * (mr2 / s - 0.5) * dsum - t * u * (4.0 * mr2 / s + t * u / (s * s) - 0.5) * dut);
  amplitude[1] +=
      charge4 * (1.0 - mr2 * r2 * (cs / (u * t) + ct / (s * u) + cu / (s * t)) -
                 mr2 * ((2.0 * mr2 + s * t / u) * dst + (2.0 * mr2 + s * u / t) * dsu + (2.0 * mr2 + t * u / s) * dut));
  amplitude[2] += charge4 * (1.0 - 2.0 * mr2 * mr2 * dsum);
  amplitude[3] +=
      charge4 * (-1.0 + (u - s) / t * (bu - bs) + (4.0 * mr2 / t + 2.0 * (s * u / (t * t) - 0.5)) * (u * cu + s * cs) -
                 2.0 * mr2 * t * (mr2 / t - 0.5) * dsum - s * u * (4.0 * mr2 / t + s * u / (t * t) - 0.5) * dsu);
  amplitude[4] +=
      charge4 * (-1.0 + (s - t) / u * (bs - bt) + (4.0 * mr2 / u + 2.0 * (t * s / (u * u) - 0.5)) * (s * cs + t * ct) -
                 2.0 * mr2 * u * (mr2 / u - 0.5) * dsum - t * s * (4.0 * mr2 / u + t * s / (u * u) - 0.5) * dst);
}

// Compute one crossed W-loop helicity-conserving amplitude
Complex WConserving(const double x, const double y, const double z, const double mr2, const Complex &bx,
                    const Complex &by, const Complex &bz, const Complex &cx, const Complex &cy, const Complex &cz,
                    const Complex &dxy, const Complex &dxz, const Complex &dyz) {
  const Complex dsum = dxy + dxz + dyz;
  return -1.0 + (z - y) / x * (bz - by) + (4.0 * mr2 / x + 2.0 * (y * z / (x * x) - 4.0 / 3.0)) * (z * cz + y * cy) -
         (2.0 * mr2 * x * (mr2 / x - 4.0 / 3.0) + 2.0 * x * x / 3.0) * dsum -
         y * z * (4.0 * mr2 / x + y * z / (x * x) - 4.0 / 3.0) * dyz;
}

// Subtract the charged W loop from the five independent amplitudes
void AddW(const double s, const double t, const double u, const double mass, Helicity5 &amplitude) {
  const double     mr2 = mass * mass;
  const Complex    m2{mr2, -1.0e-30};
  const Complex    bs     = B0(s, m2);
  const Complex    bt     = B0(t, m2);
  const Complex    bu     = B0(u, m2);
  const Complex    cs     = C0(s, m2);
  const Complex    ct     = C0(t, m2);
  const Complex    cu     = C0(u, m2);
  const Complex    dst    = D0(s, t, m2);
  const Complex    dsu    = D0(s, u, m2);
  const Complex    dut    = D0(u, t, m2);
  const Complex    dsum   = dst + dsu + dut;
  const double     r2     = s * s + t * t + u * u;
  constexpr double weight = -1.5;

  amplitude[0] += weight * WConserving(s, t, u, mr2, bs, bt, bu, cs, ct, cu, dst, dsu, dut);
  amplitude[1] +=
      weight * (1.0 - mr2 * r2 * (cs / (u * t) + ct / (s * u) + cu / (s * t)) -
                mr2 * ((2.0 * mr2 + s * t / u) * dst + (2.0 * mr2 + s * u / t) * dsu + (2.0 * mr2 + t * u / s) * dut));
  amplitude[2] += weight * (1.0 - 2.0 * mr2 * mr2 * dsum);
  amplitude[3] += weight * WConserving(t, s, u, mr2, bt, bs, bu, ct, cs, cu, dst, dut, dsu);
  amplitude[4] += weight * WConserving(u, t, s, mr2, bu, bt, bs, cu, ct, cs, dut, dsu, dst);
}

// Compute true when every complex amplitude component is finite
bool Finite(const std::array<Complex, 16> &amplitude) {
  return std::all_of(amplitude.cbegin(), amplitude.cend(),
                     [](const Complex &value) { return std::isfinite(value.real()) && std::isfinite(value.imag()); });
}

}  // namespace

// Evaluate all helicity amplitudes for real photon scattering
bool SMHelicityAmplitudes(const double s, const double t, const double u, const double alpha, const SMParam &mass,
                          std::array<std::complex<double>, 16> &amplitude) {
  amplitude.fill({0.0, 0.0});
  const double scale = std::max({std::abs(s), std::abs(t), std::abs(u)});
  if (!std::isfinite(scale) || scale <= 0.0 || !std::isfinite(alpha) || alpha <= 0.0 ||
      std::abs(s + t + u) > 1.0e-8 * scale || s <= 0.0 || t >= 0.0 || u >= 0.0) {
    return false;
  }

  Helicity5        independent = {};
  // PDG light-quark masses are MSbar at 2 GeV and c,b at their own mass
  AddFermion(s, t, u, mass.e, 1.0, independent);
  AddFermion(s, t, u, mass.mu, 1.0, independent);
  AddFermion(s, t, u, mass.tau, 1.0, independent);
  AddFermion(s, t, u, mass.d, 1.0 / 27.0, independent);
  AddFermion(s, t, u, mass.u, 16.0 / 27.0, independent);
  AddFermion(s, t, u, mass.s, 1.0 / 27.0, independent);
  AddFermion(s, t, u, mass.c, 16.0 / 27.0, independent);
  AddFermion(s, t, u, mass.b, 1.0 / 27.0, independent);
  // Use the pole extraction for the top propagator mass
  AddFermion(s, t, u, mass.t, 16.0 / 27.0, independent);
  AddW(s, t, u, mass.w, independent);

  constexpr std::array<std::size_t, 16> map  = {0, 1, 1, 2, 1, 4, 3, 1, 1, 3, 4, 1, 2, 1, 1, 0};
  const double                          norm = 8.0 * alpha * alpha;
  for (const auto &row : indices(amplitude)) { amplitude[row] = norm * independent[map[row]]; }
  return Finite(amplitude);
}

}  // namespace gra::lbyl

namespace gra {
using gra::aux::indices;
namespace {

// Compute the helicity label for one lexicographic row and leg
int LightByLightHelicity(const std::size_t row, const std::size_t leg) {
  const std::size_t divisor = std::size_t{1} << (3 - leg);
  return ((row / divisor) % 2) == 0 ? -1 : 1;
}

// Compute the invariant hard scattering cosine
bool HardScatteringCosine(const M4Vec &incoming, const M4Vec &outgoing, const double s, double &cosine) {
  if (!std::isfinite(s) || s <= 0.0) { return false; }
  const auto            p1        = incoming.Contravariant<long double>();
  const auto            p3        = outgoing.Contravariant<long double>();
  const long double     product   = gra::MinkowskiProduct(p1, p3);
  const long double     value     = 1.0L - 4.0L * product / s;
  constexpr long double tolerance = 1.0e-12L;
  if (!std::isfinite(value) || value < -1.0L - tolerance || value > 1.0L + tolerance) { return false; }
  cosine = std::clamp(static_cast<double>(value), -1.0, 1.0);
  return true;
}

}  // namespace

// Compute the exact analytic light-by-light process definition
std::vector<amplitude::Process> LightByLightProcesses() {
  const AmplitudeTopology topology = {AmplitudeTopologyNode{{PDG::PDG_gamma}, {}},
                                      AmplitudeTopologyNode{{PDG::PDG_gamma}, {}}};
  return {amplitude::Process{"SM_YY_YY", "yy_yy", "gamma gamma", "gamma gamma", "a,a",
                             "Standard Model charged-particle one-loop amplitude", topology,
                             std::vector<AmplitudeTopology>{topology}, std::vector<int>{PDG::PDG_gamma, PDG::PDG_gamma},
                             DecayStructure{DecayType::Full, true},
                             amplitude::MatrixElementForm::Generated, amplitude::TopologyMode::Exact}};
}

// Construct the analytic light-by-light amplitude
AMP_yy_yy::AMP_yy_yy(const SMParam &sm)
    : PhotonMG5Process(LightByLightProcesses(), 1.0 / sm.alpha_em_inv), sm_(sm) {}

// Compute the number of subprocess matrix elements
std::size_t AMP_yy_yy::SubprocessCount() const { return 1; }

// The charged-particle loop has zero power of the strong coupling
int AMP_yy_yy::AlphaSPower() const noexcept { return 0; }

// Compute the Born matrix-element power of alpha_QED
int AMP_yy_yy::AlphaQEDPower() const noexcept { return 4; }

// Evaluate the complete charged Standard Model loop amplitude
mg5helas::MatrixElementEvaluation AMP_yy_yy::Evaluate(LORENTZSCALAR &lts, const double alpha_s,
                                                         const bool coherent_epa) {
  (void)alpha_s;
  amplitudes_.clear();
  lts.epa_hard.Clear();
  lts.hard_color_flows.clear();
  lts.hamp.clear();
  if (!MatchProcess(lts.decaytree).has_value() || lts.decaytree.size() != 2) {
    return {mg5helas::EvaluationStatus::AmplitudeFailure, 0.0};
  }

  final_buffer_.clear();
  final_buffer_.reserve(2);
  for (const auto &branch : lts.decaytree) {
    const double energy2 = branch.p4.E() * branch.p4.E();
    const double mass2   = branch.p4.M2();
    if (!std::isfinite(energy2) || !std::isfinite(mass2) || std::abs(mass2) > 1.0e-8 * std::max(1.0, energy2)) {
      return {mg5helas::EvaluationStatus::KinematicsFailure, 0.0};
    }
    final_buffer_.push_back(branch.p4);
  }

  mg5helas::EPAHardFrame event_frame;
  std::array<M4Vec, 2>   incoming;
  const bool             epa_hard =
      coherent_epa && mg5helas::HasEPAHardBeamState(lts);
  if (epa_hard) {
    if (!mg5helas::PrepareEPAHardFrame(lts, final_buffer_, event_frame)) {
      return {mg5helas::EvaluationStatus::KinematicsFailure, 0.0};
    }
    incoming = event_frame.incoming;
  } else if (!mg5helas::PrepareOnShellKinematics(lts, final_buffer_, incoming[0], incoming[1])) {
    return {mg5helas::EvaluationStatus::KinematicsFailure, 0.0};
  }
  double cosine = 0.0;
  if (!HardScatteringCosine(incoming[0], final_buffer_[0], lts.s_hat, cosine)) {
    return {mg5helas::EvaluationStatus::KinematicsFailure, 0.0};
  }
  const double                         t = -0.5 * lts.s_hat * (1.0 - cosine);
  const double                         u = -0.5 * lts.s_hat * (1.0 + cosine);
  std::array<std::complex<double>, 16> hard;
  const double alpha = qed::alpha_QED(lts.s_hat, lts.model_cache->Tune().Structure().QED_alpha, DefaultAlphaQED());
  if (!lbyl::SMHelicityAmplitudes(lts.s_hat, t, u, alpha, sm_, hard)) {
    return {mg5helas::EvaluationStatus::AmplitudeFailure, 0.0};
  }

  const double                             azimuth = final_buffer_[0].Pt() > 1.0e-14 * std::sqrt(lts.s_hat)
                                                         ? std::atan2(final_buffer_[0].Py(), final_buffer_[0].Px())
                                                         : 0.0;
  std::vector<mg5helas::HelicityComponent> components;
  components.reserve(hard.size());
  for (const auto &row : indices(hard)) {
    const int h1 = LightByLightHelicity(row, 0);
    const int h2 = LightByLightHelicity(row, 1);
    const int h3 = LightByLightHelicity(row, 2);
    const int h4 = LightByLightHelicity(row, 3);
    components.push_back(
        {{h1, h2}, {h3, h4}, 0, std::polar(1.0, static_cast<double>(h1 - h2) * azimuth) * hard[row], {}});
  }

  const double normalization = std::sqrt(mg5helas::AppliedFinalStateSymmetryFactor(lts) / 8.0);
  amplitudes_                = epa_hard ? mg5helas::ContractEPAHardPhotonSources(lts, components, event_frame)
                                        : mg5::PhotonAmplitudes(lts, components, incoming[0], incoming[1], coherent_epa);
  Scale(amplitudes_, normalization);
  const double amp2 = SquaredNorm(amplitudes_);
  if (amplitudes_.empty() || !std::isfinite(amp2) || !AllFinite(amplitudes_)) {
    amplitudes_.clear();
    return {mg5helas::EvaluationStatus::AmplitudeFailure, 0.0};
  }
  if (epa_hard) {
    for (auto &component : components) { component.value *= normalization; }
    lts.epa_hard.frame         = event_frame;
    lts.epa_hard.amplitude     = std::move(components);
    lts.epa_hard.normalization = 2.0;
  }
  lts.hamp = amplitudes_;
  Scale(lts.hamp, 2.0);
  return {mg5helas::EvaluationStatus::Success, amp2};
}

// Clear color tags for the colorless final state
bool AMP_yy_yy::SampleColorFlow(LORENTZSCALAR &lts, MRandom &random) {
  (void)random;
  return MatchProcess(lts.decaytree).has_value() && ClearHardColorFlow(lts);
}

// Both final-state photons are color singlets
bool AMP_yy_yy::HasFinalStateColor() const { return false; }

}  // namespace gra
