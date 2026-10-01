// Electromagnetic absorption and deposited excitation for nuclear decay
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Nuclear/MBreakup.h"

#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <limits>
#include <numeric>
#include <set>
#include <stdexcept>
#include <utility>

#include "Graniitti/Math/MAlgorithms.h"
#include "Graniitti/Math/MIntegration.h"
#include "Graniitti/Math/MInterpolation.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MStatistics.h"
#include "Graniitti/Nuclear/MPhoton.h"
#include "Graniitti/Particle/MForm.h"
#include "Graniitti/Particle/MPDG.h"
#include "Graniitti/Tech/MAux.h"

using gra::aux::indices;
using gra::math::ParallelFor;

namespace gra::nuclear {
namespace {

// Find the measured response for one isotope
template <class T>
const typename T::value_type* FindIsotope(const T& isotopes, const unsigned int a, const unsigned int z) {
  const auto value =
      std::find_if(isotopes.cbegin(), isotopes.cend(), [a, z](const auto& item) { return item.a == a && item.z == z; });
  return value == isotopes.cend() ? nullptr : &*value;
}

// Compute the Levinger quasi-deuteron photoabsorption component [mb]
// sigma_QD = L NZ/A sigma_d(E) P_Pauli(E)
// [REFERENCE: Klusek-Gawenda et al., Phys. Rev. C89 (2014) 054907]
double QuasiDeuteron(const unsigned int a, const unsigned int z, const double energy, const QDParam& param) {
  if (!(energy > param.deuteron.threshold)) { return 0.0; }
  const double n       = static_cast<double>(a - z);
  const double mass    = static_cast<double>(a);
  const double charge  = static_cast<double>(z);
  const double excess  = 1.0 - param.deuteron.threshold / energy;
  const double sigma_d = param.deuteron.norm * excess * std::sqrt(excess) * std::pow(energy, -1.5);
  double       pauli   = 0.0;
  if (energy < param.pauli.low_edge) {
    pauli = std::exp(-param.pauli.low_exp / energy);
  } else if (energy < param.pauli.high_edge) {
    pauli = param.pauli.poly.back();
    for (std::size_t i = param.pauli.poly.size() - 1; i > 0; --i) { pauli = param.pauli.poly[i - 1] + energy * pauli; }
  } else {
    pauli = std::exp(-param.pauli.high_exp / energy);
  }
  return param.levinger * n * charge / mass * sigma_d * std::max(0.0, pauli);
}

// Compute one Gaussian nucleon-resonance photoabsorption component [mb]
// sigma_R(E) = A_R exp[-(E - E_R)^2 / (2 Gamma_R^2)] / (sqrt(2 pi) Gamma_R)
// [REFERENCE: Klusek-Gawenda et al., Phys. Rev. C89 (2014) 054907]
double Resonance(const double energy, const double threshold, const ResonanceIsotope& param) {
  if (energy < threshold) { return 0.0; }
  const double pull = (energy - param.energy) / param.width;
  return param.area * std::exp(-0.5 * pull * pull) / (param.width * std::sqrt(2.0 * math::PI));
}

// Compute a smoothly switched additive photonuclear continuum [mb]
// The fitted exponential term decays into the squared-log asymptote without a hard splice
// [REFERENCE: arXiv:1311.1938, Eq. (3.7), high-energy term only]
double Continuum(const double energy, const ContinuumParam& param) {
  if (!(energy > param.threshold)) { return 0.0; }
  const double u         = std::min(1.0, (energy - param.threshold) / (param.match - param.threshold));
  const double excess    = energy - param.mean;
  const double logarithm = std::log(energy / param.omega0);
  return u * u * (3.0 - 2.0 * u) *
         (param.constant + param.log2 * logarithm * logarithm + param.norm * excess * std::exp(-excess / param.width));
}

// Resolve the dimensionless photon kernel down to double precision relative to K(0) = 1
double KernelLimit(const double gamma) {
  constexpr double tolerance = std::numeric_limits<double>::epsilon();
  double upper = 1.0;
  while (PhotonKernel(upper, gamma) > tolerance && upper < 1024.0) { upper *= 2.0; }
  if (!(upper < 1024.0) || !std::isfinite(upper)) {
    throw std::invalid_argument("MBreakup: EPA kernel support is unresolved");
  }
  double lower = upper > 1.0 ? 0.5 * upper : 0.0;
  for (unsigned int iteration = 0; iteration < 64; ++iteration) {
    const double middle = 0.5 * (lower + upper);
    if (PhotonKernel(middle, gamma) > tolerance) {
      lower = middle;
    } else {
      upper = middle;
    }
  }
  return upper;
}

// Validate one isotope list and return its unique identities
template <typename T>
bool ValidIsotopes(const std::vector<T>& isotope) {
  std::set<std::pair<unsigned int, unsigned int>> identity;
  for (const auto& item : isotope) {
    if (item.a <= 1 || item.z == 0 || item.z >= item.a || !identity.insert({item.a, item.z}).second) { return false; }
  }
  return true;
}

// Compute normalized absorption weights after cancelling the energy independent charge factor
std::discrete_distribution<std::size_t> PhotonSpectrum(const std::vector<double>& omega,
                                                       const std::vector<double>& weight,
                                                       const std::vector<double>& absorption, double b, double gamma,
                                                       double limit) {
  return {omega.size(), 0.0, static_cast<double>(omega.size()), [&](double x) {
            const auto   i = static_cast<std::size_t>(x);
            const double k = omega[i] * b / (gamma * PDG::GeV2fm);
            return k > limit ? 0.0 : weight[i] * absorption[i] * PhotonKernel(k, gamma) / omega[i];
          }};
}

}  // namespace

// Validate the same GDR coefficients and isotope overrides for absorption and decay
void ValidateGDR(const GDRParam& param) {
  const auto& g = param.systematics;
  for (const double value : {g.energy_a, g.energy_b, g.width_norm, g.width_power, g.trk, g.strength}) {
    if (!std::isfinite(value) || !(value > 0.0)) { throw std::invalid_argument("GDR: invalid systematics"); }
  }
  if (!ValidIsotopes(param.isotope)) { throw std::invalid_argument("GDR: invalid isotope identities"); }
  for (const auto& item : param.isotope) {
    for (const double value : {item.energy, item.width, item.strength}) {
      if (!std::isfinite(value) || !(value > 0.0)) { throw std::invalid_argument("GDR: invalid isotope response"); }
    }
  }
}

// Resolve the common E1 systematics and measured isotope overrides
// [REFERENCE: Capote et al., Nucl. Data Sheets 110 (2009) 3107]
Dipole::Dipole(unsigned int a, unsigned int z, const GDRParam& param) {
  const auto& g   = param.systematics;
  energy          = g.energy_a * std::pow(a, -1.0 / 3.0) + g.energy_b * std::pow(a, -1.0 / 6.0);
  width           = g.width_norm * std::pow(energy, g.width_power);
  double strength = g.strength;
  if (const auto* isotope = FindIsotope(param.isotope, a, z)) {
    energy   = isotope->energy;
    width    = isotope->width;
    strength = isotope->strength;
  }
  peak = 2.0 * strength * g.trk * z * (a - z) / (math::PI * a * width);
  energy *= 1.0e-3;
  width *= 1.0e-3;
}

// Compute the standard Lorentzian giant dipole component [mb]
// sigma_GDR(E) = sigma_0 Gamma^2 E^2 / [(E^2 - E_0^2)^2 + Gamma^2 E^2]
// [REFERENCE: Berman and Fultz, Rev. Mod. Phys. 47 (1975) 713]
double Dipole::Sigma(const double omega) const {
  const double ratio        = omega / energy;
  const double reduced      = ratio > 1.0 ? 1.0 / ratio : ratio;
  const double scaled_width = (width / energy) * reduced;
  const double offset       = reduced * reduced - 1.0;
  return peak * scaled_width * scaled_width / (offset * offset + scaled_width * scaled_width);
}

// Construct one immutable electromagnetic breakup model
// E_GDR = a A^(-1/3) + b A^(-1/6), Gamma = c E_GDR^d
// sigma_GDR^peak = 2 s S_TRK / (pi Gamma), S_TRK = k NZ/A
// [REFERENCE: Capote et al., Nucl. Data Sheets 110 (2009) 3107]
MBreakup::MBreakup(BreakupParam param, std::shared_ptr<const MNucleus> emitter, form::ParamStore structure)
    : param_(std::move(param)), emitter_(std::move(emitter)), structure_(std::move(structure)) {
  if (!std::isfinite(param_.z_emit) || !std::isfinite(param_.gamma)) {
    throw std::invalid_argument("MBreakup: runtime controls must be finite");
  }
  if (param_.a == 0) {
    if (param_.z != 0 || std::fpclassify(param_.z_emit) != FP_ZERO || std::fpclassify(param_.gamma) != FP_ZERO ||
        param_.emitter != EmitterType::Point || emitter_ != nullptr) {
      throw std::invalid_argument("MBreakup: a disabled model must use default runtime controls");
    }
    return;
  }

  ValidateGDR(param_.photo.gdr);
  const auto& qd        = param_.photo.qd;
  const auto& pauli     = qd.pauli;
  const auto& continuum = param_.photo.continuum;
  if (continuum.a != param_.a || continuum.z != param_.z) {
    throw std::invalid_argument("MBreakup: photonuclear continuum fit does not match the target isotope");
  }
  const bool emitter_valid = param_.emitter == EmitterType::Point || param_.emitter == EmitterType::Proton ||
                             param_.emitter == EmitterType::Nuclear;
  const auto& response = param_.response;
  const auto& profile  = param_.profile;
  // Check positive coefficients together within each physical response
  const auto positive = [](std::initializer_list<double> values) {
    return std::all_of(values.begin(), values.end(), [](double value) { return std::isfinite(value) && value > 0.0; });
  };
  if (!positive({qd.deuteron.threshold, qd.deuteron.norm, qd.levinger, pauli.low_edge, pauli.high_edge, pauli.low_exp,
                 pauli.high_exp}) ||
      !positive({param_.photo.energy_min, param_.photo.resonance.threshold, continuum.threshold, continuum.match,
                 continuum.omega0, continuum.constant, continuum.norm, continuum.width}) ||
      !positive({param_.transfer.E0, param_.transfer.match}) ||
      !positive({response.rel_tol, response.tail_rel_tol, profile.b_min, param_.proton.q_max}) ||
      !std::isfinite(continuum.log2) || continuum.log2 < 0.0 || !std::isfinite(continuum.mean) ||
      continuum.mean < 0.0 || !gra::AllFinite(pauli.poly)) {
    throw std::invalid_argument("MBreakup: invalid photonuclear coefficients");
  }
  if (!emitter_valid || param_.a <= 1 || param_.z == 0 || param_.z >= param_.a || !(param_.z_emit > 0.0) ||
      !(param_.gamma > 1.0) || !ValidIsotopes(param_.photo.resonance.isotope)) {
    throw std::invalid_argument("MBreakup: invalid isotope response");
  }
  if (!(pauli.low_edge > qd.deuteron.threshold) || !(pauli.high_edge > pauli.low_edge) ||
      !(continuum.threshold > param_.photo.resonance.threshold) || !(continuum.match > continuum.threshold) ||
      !(continuum.omega0 > continuum.match) || !(continuum.mean < continuum.threshold)) {
    throw std::invalid_argument("MBreakup: unordered photoabsorption energy ranges");
  }
  if (response.nodes < 32 || response.nodes > 4096 || !(response.rel_tol < 1.0) || !(response.tail_rel_tol < 1.0) ||
      profile.nodes < 32 || profile.nodes > 65536 || param_.proton.q_nodes < 32 ||
      param_.proton.q_nodes > 4096) {
    throw std::invalid_argument("MBreakup: invalid quadrature controls");
  }
  for (const auto& item : param_.photo.resonance.isotope) {
    if (!std::isfinite(item.energy) || !std::isfinite(item.width) || !std::isfinite(item.area) ||
        !(item.energy > 0.0) || !(item.width > 0.0) || !(item.area > 0.0)) {
      throw std::invalid_argument("MBreakup: invalid resonance isotope entry");
    }
  }
  if ((param_.emitter == EmitterType::Nuclear) != (emitter_ != nullptr)) {
    throw std::invalid_argument("MBreakup: emitter form factor does not match");
  }

  gdr_          = Dipole(param_.a, param_.z, param_.photo.gdr);
  omega_min_    = 1.0e-3 * param_.photo.energy_min;
  kernel_limit_ = KernelLimit(param_.gamma);

  omega_edge_ = {1.0e-3 * pauli.low_edge,      1.0e-3 * param_.photo.resonance.threshold,
                 1.0e-3 * pauli.high_edge,     1.0e-3 * param_.transfer.match,
                 1.0e-3 * continuum.threshold, 1.0e-3 * continuum.match};
  std::sort(omega_edge_.begin(), omega_edge_.end());

  // Four widths is a validation margin and does not enter accepted weights
  if (!(omega_min_ < gdr_.energy) || !(1.0e-3 * continuum.match > gdr_.energy + 4.0 * gdr_.width) ||
      !(gdr_.peak > 0.0) || !std::isfinite(gdr_.energy) || !std::isfinite(gdr_.width) || !std::isfinite(gdr_.peak) ||
      !std::isfinite(kernel_limit_)) {
    throw std::invalid_argument("MBreakup: photoabsorption range does not resolve the GDR");
  }

  PrepareEnergy();
  PrepareMean();
}

// Compute whether electromagnetic breakup is active
bool MBreakup::Enabled() const { return param_.a > 0; }

// Compute the total photoabsorption cross section [mb]
// Add the fitted continuum to the GDR, quasi-deuteron and resonance response
// [REFERENCE: Klusek-Gawenda et al., Phys. Rev. C89 (2014) 054907]
double MBreakup::PhotoAbsorption(const double omega) const {
  if (!std::isfinite(omega) || omega < 0.0) {
    throw std::invalid_argument("MBreakup::PhotoAbsorption: energy must be finite and nonnegative");
  }
  if (!Enabled() || omega < omega_min_ || std::fpclassify(omega) == FP_ZERO) { return 0.0; }
  const double energy = 1000.0 * omega;
  double       sigma  = Continuum(energy, param_.photo.continuum);
  sigma += gdr_.Sigma(omega) + QuasiDeuteron(param_.a, param_.z, energy, param_.photo.qd);
  if (const auto* resonance = FindIsotope(param_.photo.resonance.isotope, param_.a, param_.z); resonance != nullptr) {
    sigma += Resonance(energy, param_.photo.resonance.threshold, *resonance);
  }
  return std::isfinite(sigma) && sigma > 0.0 ? sigma : 0.0;
}

// Compute the stable dimensionless Bessel-function EPA kernel
// K(x) = x^2 [K_1(x)^2 + K_0(x)^2 / gamma^2]
double MBreakup::BesselKernel(const double x) const { return x > kernel_limit_ ? 0.0 : PhotonKernel(x, param_.gamma); }

// Compute the emitter charge fraction inside one transverse cylinder
// C(b) = int_{r_T < b} rho_ch(r) d^3r
double MBreakup::ChargeFraction(const double b) const {
  if (!std::isfinite(b) || b < 0.0) { throw std::invalid_argument("MBreakup::ChargeFraction: invalid radius"); }
  if (param_.emitter == EmitterType::Point) { return 1.0; }
  if (param_.emitter == EmitterType::Proton) {
    if (std::fpclassify(b) == FP_ZERO) { return 0.0; }
    if (charge_b_node_.empty() || charge_fraction_.empty()) {
      throw std::logic_error("MBreakup: proton charge table is not prepared");
    }
    // Continue the enclosed charge as C(b) proportional to b^2 at the origin
    if (b <= charge_b_node_.front()) { return charge_fraction_.front() * math::pow2(b / charge_b_node_.front()); }
    if (b >= charge_b_node_.back()) { return 1.0; }
    return std::clamp(math::LinearInterpolateValidatedGrid(charge_b_node_, charge_fraction_, b), 0.0, 1.0);
  }
  return emitter_->ChargeDensity().Cylinder(b);
}

// Compute the finite-size impact-space EPA density [GeV^-1 fm^-2]
// n(omega,b) = Z^2 alpha C(b)^2 K(x) / (pi^2 omega b^2), x = omega b/(gamma hbarc)
// [REFERENCE: Baltz et al., Phys. Rept. 458 (2008) 1]
double MBreakup::PhotonDensity(const double omega, const double b) const {
  if (!std::isfinite(omega) || !(omega > 0.0) || !std::isfinite(b) || b < 0.0) {
    throw std::invalid_argument("MBreakup::PhotonDensity: expected positive energy and nonnegative finite radius");
  }
  if (!Enabled()) { return 0.0; }
  const double x = omega * b / (param_.gamma * PDG::GeV2fm);
  const double density =
      x > kernel_limit_ ? 0.0 : ImpactPhotonDensity(omega, b, param_.gamma, param_.z_emit, ChargeFraction(b));
  if (std::isnan(density) || density < 0.0) {
    throw std::runtime_error("MBreakup::PhotonDensity: invalid EPA density");
  }
  return density;
}

// Fold the photon density and photoabsorption cross section directly
// mu(b) = int n(omega,b) sigma_gammaA(omega) d omega
double MBreakup::FoldMean(const double b) const {
  if (std::fpclassify(b) == FP_ZERO && param_.emitter == EmitterType::Point) {
    return std::numeric_limits<double>::infinity();
  }
  statistics::CompensatedSum mean;
  for (const auto& i : indices(omega_node_)) {
    constexpr double mb_to_fm2 = 0.1;
    const double     density   = PhotonDensity(omega_node_[i], b);
    const double     term      = omega_weight_[i] * mb_to_fm2 * absorption_[i] * density;
    mean.Add(term);
  }
  const double value = static_cast<double>(mean.Value());
  if (std::isnan(value) || value < 0.0) { throw std::runtime_error("MBreakup::Mean: invalid folded excitation"); }
  return value;
}

// Prepare and validate the piecewise photon energy quadrature
void MBreakup::PrepareEnergy() {
  constexpr double       hbarc      = PDG::GeV2fm;
  constexpr unsigned int tail_nodes = 32;
  const double           b_ref      = 1.2 * std::cbrt(static_cast<double>(param_.a));
  const double           tail_start = 1.0e-3 * param_.photo.continuum.match;

  // Integrate one logarithmic tail interval at the target nuclear radius
  const auto tail_integral = [&](const double lower, const double upper) {
    const auto                 rule = math::GaussLegendreRule(tail_nodes, std::log(lower), std::log(upper));
    statistics::CompensatedSum sum;
    for (const auto& i : indices(rule.first)) {
      const double omega = std::exp(rule.first[i]);
      const double x     = omega * b_ref / (param_.gamma * hbarc);
      sum.Add(rule.second[i] * PhotoAbsorption(omega) * PhotonKernel(x, param_.gamma));
    }
    return static_cast<double>(sum.Value());
  };

  double       total     = tail_integral(omega_min_, tail_start);
  double       lower     = tail_start;
  unsigned int converged = 0;
  for (unsigned int iteration = 0; iteration < 128; ++iteration) {
    const double upper = 2.0 * lower;
    const double shell = tail_integral(lower, upper);
    total += shell;
    omega_edge_.push_back(upper);
    const double scale = std::max(total, std::numeric_limits<double>::min());
    converged          = shell <= param_.response.tail_rel_tol * scale ? converged + 1U : 0U;
    lower              = upper;
    if (converged >= 2U) { break; }
  }
  if (converged < 2U || !std::isfinite(total) || !(total > 0.0)) {
    throw std::runtime_error("MBreakup: photon-energy tail does not meet tail_rel_tol");
  }
  omega_max_ = lower;

  std::vector<double> interval = {omega_min_};
  std::sort(omega_edge_.begin(), omega_edge_.end());
  for (const double value : omega_edge_) {
    if (value > interval.back() && value < omega_max_) { interval.push_back(value); }
  }
  interval.push_back(omega_max_);

  const auto build = [this, &interval](const unsigned int nodes, std::vector<double>& energy,
                                       std::vector<double>& weight, std::vector<double>& absorption) {
    energy.clear();
    weight.clear();
    for (std::size_t i = 1; i < interval.size(); ++i) {
      const auto rule = math::GaussLegendreRule(nodes, std::log(interval[i - 1]), std::log(interval[i]));
      for (const auto& j : indices(rule.first)) {
        const double omega = std::exp(rule.first[j]);
        energy.push_back(omega);
        weight.push_back(rule.second[j] * omega);
      }
    }
    absorption.resize(energy.size(), 0.0);
    ParallelFor(energy.size(), [&](const std::size_t i) { absorption[i] = PhotoAbsorption(energy[i]); });
  };
  std::vector<double> coarse_energy, coarse_weight, coarse_absorption;
  build(param_.response.nodes / 2, coarse_energy, coarse_weight, coarse_absorption);
  // Check absorption at three resolved EPA energy scales
  const auto moment = [this](const std::vector<double>& energy, const std::vector<double>& weight,
                             const std::vector<double>& absorption) {
    std::array<double, 3> value{};
    for (const auto& i : indices(energy)) {
      const double                measure = weight[i] * absorption[i] / energy[i];
      const std::array<double, 3> kernel{1.0, BesselKernel(energy[i] / gdr_.energy),
                                         BesselKernel(energy[i] / omega_min_)};
      for (const auto& k : indices(kernel)) { value[k] += measure * kernel[k]; }
    }
    return value;
  };
  auto coarse = moment(coarse_energy, coarse_weight, coarse_absorption);
  for (unsigned int nodes = param_.response.nodes;; nodes = std::min(4096U, 2U * nodes)) {
    build(nodes, omega_node_, omega_weight_, absorption_);
    auto fine      = moment(omega_node_, omega_weight_, absorption_);
    bool converged = gra::AllFinite(fine) && gra::AllFinite(coarse);
    for (const auto& i : indices(fine)) {
      const double bound = param_.response.rel_tol * std::max(fine[i], coarse[i]) +
                           64.0 * std::numeric_limits<double>::epsilon() * fine[0];
      converged = converged && std::abs(fine[i] - coarse[i]) <= bound;
    }
    if (converged) { return; }
    if (nodes == 4096U) { break; }
    coarse = fine;
  }
  throw std::runtime_error("MBreakup: photoabsorption quadrature does not meet rel_tol");
}

// Tabulate the projected charge of the shared proton form factor
// C_p(b) = b / hbarc int_0^infinity J_1(qb / hbarc) G_E(q^2) dq
void MBreakup::PrepareCharge() {
  if (param_.emitter != EmitterType::Proton) { return; }
  constexpr double hbarc      = PDG::GeV2fm;
  const auto       q_rule     = math::GaussLegendreRule(param_.proton.q_nodes, 0.0, param_.proton.q_max);
  const double     resolved_b = math::PI * hbarc * static_cast<double>(param_.proton.q_nodes - 1) / param_.proton.q_max;
  const auto       stop       = std::upper_bound(log_b_node_.cbegin(), log_b_node_.cend(), std::log(resolved_b));
  const std::size_t count     = static_cast<std::size_t>(stop - log_b_node_.cbegin());
  charge_b_node_.resize(count, 0.0);
  charge_fraction_.resize(count, 0.0);
  ParallelFor(count, [&](const std::size_t i) {
    const double b        = std::exp(log_b_node_[i]);
    double       integral = 0.0;
    for (const auto& iq : indices(q_rule.first)) {
      const double q = q_rule.first[iq];
      integral += q_rule.second[iq] * std::cyl_bessel_j(1.0, q * b / hbarc) * form::G_E(q * q, structure_);
    }
    charge_b_node_[i]   = b;
    charge_fraction_[i] = std::clamp(b * integral / hbarc, 0.0, 1.0);
  });
  if (charge_b_node_.size() < 2) { throw std::runtime_error("MBreakup: proton charge grid is unresolved"); }
  for (std::size_t i = 1; i < charge_fraction_.size(); ++i) {
    charge_fraction_[i] = std::max(charge_fraction_[i], charge_fraction_[i - 1]);
  }
  math::ValidateInterpolationGrid(charge_b_node_);
}

// Prepare the validated logarithmic impact-parameter interpolation table
void MBreakup::PrepareMean() {
  constexpr double hbarc     = PDG::GeV2fm;
  const double     log_b_min = std::log(param_.profile.b_min);
  const double log_b_max = std::log(kernel_limit_) + std::log(param_.gamma) + std::log(hbarc) - std::log(omega_min_);
  const double log_limit = std::log(std::numeric_limits<double>::max());
  if (!std::isfinite(log_b_min) || !std::isfinite(log_b_max) || !(log_b_max > log_b_min) || !(log_b_max < log_limit)) {
    throw std::invalid_argument("MBreakup: invalid impact-parameter interpolation range");
  }

  log_b_node_.resize(param_.profile.nodes, 0.0);
  mean_.resize(param_.profile.nodes, 0.0);
  const double step = (log_b_max - log_b_min) / static_cast<double>(param_.profile.nodes - 1);
  for (const auto& i : indices(log_b_node_)) { log_b_node_[i] = log_b_min + static_cast<double>(i) * step; }
  PrepareCharge();
  ParallelFor(log_b_node_.size(), [&](const std::size_t i) {
    const double mean = FoldMean(std::exp(log_b_node_[i]));
    if (!std::isfinite(mean) || mean < 0.0) {
      throw std::runtime_error("MBreakup: invalid impact-parameter interpolation value");
    }
    mean_[i] = mean;
  });
  math::ValidateInterpolationGrid(log_b_node_);
}

// Compute the mean number of absorbed photons
// mu(b) = int n(omega,b) sigma_gammaA(omega) d omega
// [REFERENCE: Baltz et al., Phys. Rept. 458 (2008) 1]
double MBreakup::Mean(const double b) const {
  if (!std::isfinite(b) || b < 0.0) {
    throw std::invalid_argument("MBreakup::Mean: impact parameter must be finite and nonnegative");
  }
  if (!Enabled()) { return 0.0; }
  if (std::fpclassify(b) == FP_ZERO) { return FoldMean(b); }
  const double log_b = std::log(b);
  if (log_b < log_b_node_.front()) { return FoldMean(b); }
  if (log_b > log_b_node_.back()) { return 0.0; }
  return std::max(0.0, math::LinearInterpolateValidatedGrid(log_b_node_, mean_, log_b));
}

// Compute the derived upper impact-parameter support in fm
double MBreakup::ImpactMax() const { return Enabled() ? std::exp(log_b_node_.back()) : 0.0; }

// Compute the same normalized photon spectrum used by the compound Poisson sampler
std::vector<double> MBreakup::Spectrum(const double b) const {
  if (!Enabled() || !(Mean(b) > 0.0)) { return std::vector<double>(omega_node_.size(), 0.0); }
  return PhotonSpectrum(omega_node_, omega_weight_, absorption_, b, param_.gamma, kernel_limit_).probabilities();
}

// Sample one absorbed photon from the impact-dependent EMD spectrum
// p(omega|b) = n(omega,b) sigma_gammaA(omega) / mu(b)
double MBreakup::SampleEnergy(const double b, MRandom& random) const {
  if (!Enabled() || !std::isfinite(b) || b < 0.0) { return 0.0; }
  if (!(Mean(b) > 0.0)) { return 0.0; }
  auto spectrum = PhotonSpectrum(omega_node_, omega_weight_, absorption_, b, param_.gamma, kernel_limit_);
  return omega_node_[spectrum(random.rng)];
}

// Sample one physical deposited energy from the continuous TCM distribution
// w = min(omega,omega_match)
// p(E|w) = exp(-w/E0) delta(E - w) + [1 - exp(-w/E0)] / w for 0 <= E <= w
// [REFERENCE: Jucha et al., Phys. Rev. C111 (2025) 034901]
double MBreakup::SampleTransfer(const double omega, MRandom& random) const {
  if (!Enabled() || !std::isfinite(omega) || !(omega > 0.0)) { return 0.0; }
  const double transfer = std::min(omega, 1.0e-3 * param_.transfer.match);
  const double full     = std::exp(-transfer / (1.0e-3 * param_.transfer.E0));
  return random.U(0.0, 1.0) < full ? transfer : random.U(0.0, transfer);
}

// Sample independent photon absorptions and sum their deposited energies before one remnant decay
double MBreakup::SampleExcitation(const double b, MRandom& random, bool absorbed) const {
  const double mean = Mean(b);
  if (!(mean > 0.0)) { return 0.0; }
  const int count = random.PoissonRandom(mean, absorbed);
  if (count == 0) { return 0.0; }
  auto   spectrum   = PhotonSpectrum(omega_node_, omega_weight_, absorption_, b, param_.gamma, kernel_limit_);
  double excitation = 0.0;
  for (int i = 0; i < count; ++i) { excitation += SampleTransfer(omega_node_[spectrum(random.rng)], random); }
  return excitation;
}

}  // namespace gra::nuclear
