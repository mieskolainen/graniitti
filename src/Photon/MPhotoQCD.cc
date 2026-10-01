// QCD photoproduction utilities
//
// [REFERENCE: Cisek, Schafer, Szczurek, PRD 80 (2009), arXiv:0906.1739]
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <mutex>
#include <stdexcept>

// Own
#include "Graniitti/MGlobals.h"
#include "Graniitti/Math/MIntegration.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Particle/MForm.h"
#include "Graniitti/Photon/MPhotoQCD.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MException.h"
#include "json.hpp"

using gra::aux::indices;

using gra::math::msqrt;
using gra::math::pow2;

namespace gra {

namespace {

// Compute the CSS W0 kernel after azimuthal integration
// W0 = 1/(k2+m2)-1/sqrt[(k2-m2-kappa2)^2+4m2 k2]
double CSSW0(double k2, double kappa2, double m2) {
  const double radicand = pow2(k2 - m2 - kappa2) + 4.0 * m2 * k2;
  return 1.0 / (k2 + m2) - 1.0 / std::sqrt(radicand);
}

// Compute the CSS W1 kernel after azimuthal integration
// W1 = 1-(k2+m2)[1+(k2-m2-kappa2)/D]/(2k2)
double CSSW1(double k2, double kappa2, double m2) {
  const double radicand = pow2(k2 - m2 - kappa2) + 4.0 * m2 * k2;
  const double root = std::sqrt(radicand);
  return 1.0 - (k2 + m2) / (2.0 * k2) * (1.0 + (k2 - m2 - kappa2) / root);
}

// Integrate 1 and z^2+(1-z)^2 over the timelike denominator a-q2*z*(1-z)-i0
// [REFERENCE: Cisek, Schafer, Szczurek, arXiv:0906.1739, Eqs. (2.1)-(2.5)]
std::pair<std::complex<double>, std::complex<double>> CSSZMoments(double a, double q2) {
  if (q2 < std::sqrt(std::numeric_limits<double>::epsilon()) * a) {
    return {(1.0 + q2 / (6.0 * a)) / a, (2.0 / 3.0 + q2 / (10.0 * a)) / a};
  }
  const double               r = 4.0 * a / q2;
  const double               v = std::sqrt(std::abs(1.0 - r));
  const std::complex<double> i0 =
      r < 1.0 ? 2.0 / (q2 * v) * std::complex<double>(std::log(r) - 2.0 * std::log1p(v), math::PI)
              : 4.0 / (q2 * v) * std::atan(1.0 / v);
  return {i0, (1.0 - 2.0 * a / q2) * i0 + 2.0 / q2};
}

} // namespace

// Validate photoproduction physics steering
void MPhotoQCDParam::Validate(const std::string &block_name) const {
  const std::string context = "MPhotoQCDParam::Validate(" + block_name + ")";
  if (!std::isfinite(xg_coefficient) || xg_coefficient <= 0.0) {
    throw std::invalid_argument(context + ": xg_coefficient must be positive");
  }
  if (!std::isfinite(mu_MIN) || mu_MIN <= 0.0 || !std::isfinite(mu_over_m) ||
      mu_over_m <= 0.0) {
    throw std::invalid_argument(context +
                                ": mu_MIN and mu_over_m must be positive");
  }
  if (!std::isfinite(light_quark_mass_MIN) || light_quark_mass_MIN < 0.0) {
    throw std::invalid_argument(context +
                                ": light_quark_mass_MIN must be non-negative");
  }
  if (!std::isfinite(real_part_delta) || real_part_delta < 0.0 ||
      real_part_delta >= 1.0) {
    throw std::invalid_argument(context + ": real_part_delta must be in [0,1)");
  }
  if (!std::isfinite(t_slope_B0) || !std::isfinite(t_slope_alpha_prime) ||
      !std::isfinite(t_slope_W0) || t_slope_B0 < 0.0 ||
      t_slope_alpha_prime < 0.0 || t_slope_W0 <= 0.0) {
    throw std::invalid_argument(context + ": invalid t-slope parameters");
  }
}

// Read photoproduction physics steering from one immutable GENERAL JSON text
void MPhotoQCDParam::ConfigureFromJson(const std::string &source_file,
                                       const std::string &json_text,
                                       const std::string &block_name) {
  using json = nlohmann::json;

  try {
    Configure(json::parse(json_text), source_file, block_name);
  } catch (const json::exception &e) {
    const std::string message =
        "MPhotoQCDParam::ConfigureFromJson: Error parsing " + source_file +
        " (Check " + block_name + "): " + std::string(e.what());
    throw std::invalid_argument(message);
  }
}

// Read photoproduction physics steering from one parsed GENERAL document
void MPhotoQCDParam::Configure(const nlohmann::json &document,
                               const std::string &source_file,
                               const std::string &block_name) {
  using json = nlohmann::json;

  try {
    const auto &block = document.at(block_name);
    xg_coefficient = block.at("xg_coefficient");
    mu_MIN = block.at("mu_MIN");
    mu_over_m = block.at("mu_over_m");
    light_quark_mass_MIN = block.at("light_quark_mass_MIN");
    real_part_delta = block.at("real_part_delta");
    t_slope_B0 = block.at("t_slope_B0");
    t_slope_alpha_prime = block.at("t_slope_alpha_prime");
    t_slope_W0 = block.at("t_slope_W0");
    Validate(block_name);

    {
      std::lock_guard<std::mutex> lock(gra::g_mutex);
      std::cout << "MPhotoQCDParam::ReadParameters: [" << block_name << "]"
                << std::endl;
      std::cout << block << std::endl;
      std::cout << std::endl;
    }
  } catch (const json::exception &e) {
    const std::string message =
        "MPhotoQCDParam::Configure: Error parsing " + source_file +
        " (Check " + block_name + "): " + std::string(e.what());
    throw std::invalid_argument(message);
  }

  initialized = true;
}

// Compute the active k2 upper integration bound
double MPhotoQCDNumerics::K2Max(double q2) const {
  return k2_MAX_use_q2 ? std::max(k2_MAX, q2) : k2_MAX;
}

// Prepare fixed integration rules for repeated impact-factor calls
void MPhotoQCDNumerics::PrepareIntegrationRules() {
  kappa2_nodes.clear();
  kappa2_weights.clear();
  if (!(kappa2_MIN > 0.0) || !(kappa2_MAX > kappa2_MIN) || N_kappa == 0) {
    return;
  }

  const auto rule = gra::math::GaussLegendreRule(N_kappa, std::log(kappa2_MIN),
                                                 std::log(kappa2_MAX));
  const std::vector<double> &nodes = rule.first;
  const std::vector<double> &weights = rule.second;

  kappa2_nodes.reserve(nodes.size());
  kappa2_weights.reserve(weights.size());
  for (const auto &i : indices(nodes)) {
    const double kappa2 = std::exp(nodes[i]);
    kappa2_nodes.push_back(kappa2);
    kappa2_weights.push_back(weights[i] * kappa2);
  }
}

// Validate photoproduction numerical steering against the Sudakov grid
void MPhotoQCDNumerics::Validate(double sudakov_q2_max,
                                 const std::string &block_name) const {
  const std::string context = "MPhotoQCDNumerics::Validate(" + block_name + ")";
  if (!std::isfinite(kappa2_MIN) || !std::isfinite(kappa2_MAX) ||
      kappa2_MIN <= 0.0 || kappa2_MIN >= kappa2_MAX) {
    throw std::invalid_argument(context +
                                ": require 0 < kappa2_MIN < kappa2_MAX");
  }
  if (!std::isfinite(k2_MIN) || !std::isfinite(k2_MAX) || k2_MIN <= 0.0 ||
      k2_MIN >= k2_MAX) {
    throw std::invalid_argument(context + ": require 0 < k2_MIN < k2_MAX");
  }
  if (N_kappa == 0 || N_k == 0) { throw std::invalid_argument(context + ": N_kappa and N_k must be positive"); }
  if (!std::isfinite(sudakov_q2_max) || sudakov_q2_max <= 0.0) {
    throw std::invalid_argument(context + ": invalid NUMERICS_SUDAKOV q2_MAX");
  }
  if (kappa2_MAX > sudakov_q2_max * (1.0 + 1e-12)) {
    throw std::invalid_argument(
        context + ": kappa2_MAX must not exceed NUMERICS_SUDAKOV q2_MAX");
  }
}

// Read photoproduction numerical steering from one immutable NUMERICS JSON text
void MPhotoQCDNumerics::ConfigureFromJson(const std::string &source_file,
                                          const std::string &json_text,
                                          const std::string &block_name) {
  using json = nlohmann::json;

  try {
    Configure(json::parse(json_text), source_file, block_name);
  } catch (const json::exception &e) {
    std::string str = "MPhotoQCDNumerics::ConfigureFromJson: Error parsing " +
                      source_file + " (Check " + block_name +
                      "): " + std::string(e.what());
    throw std::invalid_argument(str);
  }
}

// Read numerical steering from one parsed NUMERICS document
void MPhotoQCDNumerics::Configure(const nlohmann::json &document,
                                  const std::string &source_file,
                                  const std::string &block_name) {
  using json = nlohmann::json;

  try {
    const auto &block = document.at(block_name);
    kappa2_MIN = block.at("kappa2_MIN");
    kappa2_MAX = block.at("kappa2_MAX");
    k2_MIN = block.at("k2_MIN");
    k2_MAX = block.at("k2_MAX");
    k2_MAX_use_q2 = block.at("k2_MAX_use_q2");

    // Validate integer counts before narrowing the JSON representation
    const auto count = [&](const std::string &name) {
      const auto &value = block.at(name);
      if (!value.is_number_integer() || value.get<long double>() < 1.0L ||
          value.get<long double>() > std::numeric_limits<unsigned int>::max()) {
        throw std::invalid_argument("MPhotoQCDNumerics::Configure(" + block_name +
                                    "): " + name + " must be a positive integer in the unsigned range");
      }
      return value.get<unsigned int>();
    };
    N_kappa = count("N_kappa");
    N_k = count("N_k");

    const double sudakov_q2_max = document.at("NUMERICS_SUDAKOV").at("q2_MAX");
    Validate(sudakov_q2_max, block_name);
    PrepareIntegrationRules();

    {
      std::lock_guard<std::mutex> lock(gra::g_mutex);
      std::cout << "MPhotoQCDNumerics::ReadParameters: [" << block_name << "]"
                << std::endl;
      std::cout << block << std::endl;
      std::cout << std::endl;
    }
  } catch (const json::exception &e) {
    std::string str = "MPhotoQCDNumerics::Configure: Error parsing " +
                      source_file + " (Check " + block_name +
                      "): " + std::string(e.what());
    throw std::invalid_argument(str);
  }

  initialized = true;
}

// Construct one immutable photoproduction physics block
MPhotoQCDParamPtr ReadPhotoQCDParam(const MModelTune &tune,
                                    const std::string &block_name) {
  auto param = std::make_shared<MPhotoQCDParam>();
  param->Configure(tune.General(), tune.GeneralFile(), block_name);
  return param;
}

// Construct one immutable photoproduction numerical block
MPhotoQCDNumericsPtr ReadPhotoQCDNumerics(const MModelTune &tune,
                                          const std::string &block_name) {
  auto numerics = std::make_shared<MPhotoQCDNumerics>();
  numerics->Configure(tune.Numerics(), tune.NumericsFile(), block_name);
  return numerics;
}

// Compute one run owned photoproduction physics block
MPhotoQCDParamPtr GetPhotoQCDParam(MModelCache &cache,
                                   const std::string &block_name) {
  return cache.Get<MPhotoQCDParam>(
      "photoqcd:param:" + block_name, [&cache, &block_name] {
        return ReadPhotoQCDParam(cache.Tune(), block_name);
      });
}

// Compute one run owned photoproduction numerical block
MPhotoQCDNumericsPtr GetPhotoQCDNumerics(MModelCache &cache,
                                         const std::string &block_name) {
  return cache.Get<MPhotoQCDNumerics>(
      "photoqcd:numerics:" + block_name, [&cache, &block_name] {
        return ReadPhotoQCDNumerics(cache.Tune(), block_name);
      });
}

// Acquire shared read-only Sudakov and UGD tables
void MPhotoQCD::EnsureSudakov(gra::LORENTZSCALAR &lts,
                              const SoftModelPtr &soft_model,
                              const std::string &context) {
  if (!soft_model) {
    throw std::invalid_argument(context + ": missing SOFT model");
  }
  if (lts.s <= 0.0) {
    lts.s = (lts.pbeam1 + lts.pbeam2).M2();
  }
  if (lts.sqrt_s <= 0.0 && lts.s > 0.0) {
    lts.sqrt_s = std::sqrt(lts.s);
  }
  if (lts.LHAPDFSET.empty() || lts.LHAPDFSET == "null") {
    throw std::invalid_argument(
        context + ": LHAPDFSET must be configured for UGD access");
  }
  if (lts.GlobalSudakovPtr == nullptr) {
    if (lts.model_cache == nullptr) {
      throw std::invalid_argument(context + ": missing run owned model caches");
    }
    lts.GlobalSudakovPtr = lts.model_cache->sudakov.GetSudakov(
        lts.sqrt_s, lts.LHAPDFSET, soft_model);
  }
  if (lts.GlobalSudakovPtr->SoftModelHandle() != soft_model) {
    throw std::invalid_argument(
        context + ": Sudakov and photoproduction SOFT models differ");
  }
}

// Compute a positive invariant mass squared for the generated central system
double MPhotoQCD::CentralMass2(const gra::LORENTZSCALAR &lts) {
  return (lts.decaytree[0].p4 + lts.decaytree[1].p4).M2();
}

// Compute the steerable real-part correction for high-energy vector production
// factor = rho+i with rho = Re(A)/Im(A)
std::complex<double> MPhotoQCD::RealPartFactor(const MPhotoQCDParam &param) {
  const double rho = std::tan(0.5 * gra::math::PI * param.real_part_delta);
  return std::complex<double>(rho, 1.0);
}

// Compute the elastic momentum-transfer slope factor around the UGD impact factor
// factor = exp[B(W)t/2]
double MPhotoQCD::TSlopeFactor(const MPhotoQCDParam &param, double w2,
                               double t) {
  const double slope = TSlope(param, w2);
  return gra::form::ExpSlopeAmplitude(slope, t);
}

// Compute the elastic slope or inclusive dissociative target transition factor
double MPhotoQCD::TargetTransitionFactor(const ForwardLegState &state,
                                         const MPhotoQCDParam &param,
                                         const double w2,
                                         const SoftModel &soft_model) {
  if (!state.IsExcited()) {
    return TSlopeFactor(param, w2, state.t);
  }

  const double factor = soft_model.ForwardExcitationFactor(
      soft_model.ForwardExcitationExchange(), state.t, state.mass2);
  return std::isfinite(factor) && factor > 0.0 ? factor : 0.0;
}

// Compute the cross-section t slope associated with the amplitude factor
double MPhotoQCD::TSlope(const MPhotoQCDParam &param, double w2) {
  const double log_arg = std::max(w2 / pow2(param.t_slope_W0), 1e-12);
  return std::max(0.0, param.t_slope_B0 +
                           2.0 * param.t_slope_alpha_prime * std::log(log_arg));
}

// Expand one spin-averaged source over explicit incoming proton helicity rows
std::vector<std::complex<double>> MPhotoQCD::InitialProtonSpinCopies(
    const std::vector<std::complex<double>> &amplitudes) {
  std::vector<std::complex<double>> out;
  out.reserve(4 * amplitudes.size());
  for (int spin_row = 0; spin_row < 4; ++spin_row) {
    out.insert(out.end(), amplitudes.begin(), amplitudes.end());
  }
  return out;
}

// Compute the small-x value used by the impact-factor model
double MPhotoQCD::XGluon(const MPhotoQCDParam &param, double q2, double w2) {
  return (w2 > 0.0) ? param.xg_coefficient * q2 / w2 : 0.0;
}

// Compute the UGD hard scale used by the impact-factor model
double MPhotoQCD::HardScale(const MPhotoQCDParam &param, double q2) {
  return std::max(param.mu_MIN, param.mu_over_m * msqrt(q2));
}

// Contract a lower-index current with an upper-index polarization vector
std::complex<double> MPhotoQCD::ContractCurrentPolarization(
    const MDirac::Current &current,
    const FTensor::Tensor1<std::complex<double>, 4> &eps) {
  std::complex<double> out = 0.0;
  for (std::size_t mu = 0; mu < 4; ++mu) {
    out += current[mu] * eps(mu);
  }
  return out;
}

// Transport the beam transverse rest axes through the collider CM to the lab
// The dual circular states obey sum_m exp(i m phi) eps_m/sqrt(2) = ex cos(phi)+ey sin(phi)
// [REFERENCE: Zha et al., Phys. Rev. D 103 (2021) 033007, Eq. (2), beam axis approximation]
std::array<FTensor::Tensor1<std::complex<double>, 4>, 2>
MPhotoQCD::TransverseStates(const LORENTZSCALAR &lts, const M4Vec &boson) {
  const M4Vec collider = lts.pbeam1 + lts.pbeam2;
  const double mass = boson.M();
  const double energy = collider.M();
  if (!std::isfinite(mass) || !std::isfinite(energy) || !(mass > 0.0) || !(energy > 0.0) ||
      !(boson.E() > 0.0) || !(collider.E() > 0.0)) {
    throw AmplitudeFailure("MPhotoQCD::TransverseStates: invalid timelike production state");
  }
  M4Vec rest = boson;
  kinematics::LorentzBoost(collider, energy, rest, -1);
  if (!(rest.E() > 0.0)) { throw AmplitudeFailure("MPhotoQCD::TransverseStates: invalid collider boost"); }
  M4Vec ex(1.0, 0.0, 0.0, 0.0);
  M4Vec ey(0.0, 1.0, 0.0, 0.0);
  kinematics::LorentzBoost(rest, mass, ex, 1);
  kinematics::LorentzBoost(rest, mass, ey, 1);
  kinematics::LorentzBoost(collider, energy, ex, 1);
  kinematics::LorentzBoost(collider, energy, ey, 1);
  const auto x = ex.Contravariant();
  const auto y = ey.Contravariant();
  if (!gra::AllFinite(x) || !gra::AllFinite(y) ||
      std::abs(ex.M2() + 1.0) > 1.0e-10 * std::max(1.0, gra::SquaredNorm(x)) ||
      std::abs(ey.M2() + 1.0) > 1.0e-10 * std::max(1.0, gra::SquaredNorm(y))) {
    throw AmplitudeFailure("MPhotoQCD::TransverseStates: invalid polarization boost");
  }
  std::array<FTensor::Tensor1<std::complex<double>, 4>, 2> states;
  for (const auto &i : indices(states)) {
    const int m = i == 0 ? -1 : 1;
    for (const auto &mu : indices(x)) { states[i](mu) = (x[mu] - math::zi * static_cast<double>(m) * y[mu]) / std::sqrt(2.0); }
  }
  return states;
}

// Compute the CSS flavour integral with analytic z moments and the timelike pole phase
std::complex<double> MPhotoQCD::CSSFlavorIntegral(const gra::LORENTZSCALAR &lts, const MPhotoQCDParam &param,
                                                  const MPhotoQCDNumerics &num, double xg, double quark_mass, double q2,
                                                  double mu) {
  const double m2               = pow2(std::max(quark_mass, param.light_quark_mass_MIN));
  const double pole             = q2 / 4.0 - m2;
  const double maximum          = num.K2Max(q2);
  const auto [nodes, weights]   = math::GaussLegendreRule(num.N_k, 0.0, 1.0);
  std::complex<double> integral = 0.0;
  try {
    for (const auto &j : indices(num.kappa2_nodes)) {
      const double kappa2 = num.kappa2_nodes[j];
      const double fg     = lts.GlobalSudakovPtr->fg_xQ2Mu(xg, kappa2, std::max(mu, lts.GlobalSudakovPtr->GetMuMin()));
      std::vector<double> bounds = {num.k2_MIN, maximum};
      for (const double split : {pole, kappa2}) {
        if (split > num.k2_MIN && split < maximum) { bounds.push_back(split); }
      }
      std::sort(bounds.begin(), bounds.end());
      for (std::size_t b = 1; b < bounds.size(); ++b) {
        if (!(bounds[b] > bounds[b - 1])) { continue; }
        const double low    = std::log(bounds[b - 1]);
        const double length = std::log(bounds[b] / bounds[b - 1]);

        // A sine map removes the square-root singularity at the pair threshold
        for (const auto &i : indices(nodes)) {
          const double angle  = math::PI * nodes[i];
          const double k2     = std::exp(low + length * pow2(std::sin(0.5 * angle)));
          const double jac    = k2 * length * math::PI / 2.0 * std::sin(angle);
          const auto [i0, i1] = CSSZMoments(k2 + m2, q2);
          const auto   kernel = m2 * CSSW0(k2, kappa2, m2) * i0 + k2 / (k2 + m2) * CSSW1(k2, kappa2, m2) * i1;
          const double alpha_fg =
              fg * lts.GlobalSudakovPtr->AlphaS_Q2(std::max({k2 + m2, kappa2, lts.GlobalSudakovPtr->GetQ2Min()}));
          integral += num.kappa2_weights[j] * weights[i] * jac * alpha_fg * kernel / pow2(kappa2);
        }
      }
    }
  } catch (const std::exception &error) {
    throw AmplitudeFailure(std::string("MPhotoQCD::CSSFlavorIntegral: ") + error.what());
  }
  integral *= math::pow3(math::PI);
  if (!std::isfinite(integral.real()) || !std::isfinite(integral.imag())) {
    throw AmplitudeFailure("MPhotoQCD::CSSFlavorIntegral: nonfinite quark impact factor");
  }
  return integral;
}

} // namespace gra
