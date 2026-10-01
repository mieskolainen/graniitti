// Regge form factors and exchanged-meson kernels
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Regge/MReggeForm.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <set>
#include <stdexcept>

#include "Graniitti/MGlobals.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MException.h"

using gra::aux::indices;
using gra::math::pow2;

namespace gra::regge {
namespace {

// Compute 2 log cosh(d/2) for a finite absolute rapidity gap without cancellation or overflow
double SmoothReggeInterval(const double d) {
  if (d < 1.0) { return 2.0 * std::log1p(2.0 * pow2(std::sinh(0.25 * d))); }
  return d + 2.0 * (std::log1p(std::exp(-d)) - std::log(2.0));
}

// Require an exact field inventory for one form factor object
void CheckFields(const nlohmann::json &block, const std::set<std::string> &fields, const std::string &context) {
  if (!block.is_object()) { throw std::invalid_argument(context + " must be an object"); }
  for (const auto &[field, value] : block.items()) {
    (void)value;
    if (!fields.contains(field)) { throw std::invalid_argument(context + " has unknown field " + field); }
  }
  for (const auto &field : fields) {
    if (!block.contains(field)) { throw std::invalid_argument(context + " is missing field " + field); }
  }
}

// Read one finite numeric form factor field
double ReadNumber(const nlohmann::json &block, const std::string &field, const std::string &context) {
  if (!block.at(field).is_number()) { throw std::invalid_argument(context + "." + field + " must be numeric"); }
  const double value = block.at(field).get<double>();
  if (!std::isfinite(value)) { throw std::invalid_argument(context + "." + field + " must be finite"); }
  return value;
}

}  // namespace

// Parse one steering card form factor type
FFType ParseFFType(const std::string &type) {
  if (type == "none") { return FFType::None; }
  if (type == "exp") { return FFType::Exponential; }
  if (type == "power") { return FFType::Power; }
  if (type == "orear") { return FFType::Orear; }
  if (type == "gaussian") { return FFType::Gaussian; }
  if (type == "vector") { return FFType::Vector; }
  if (type == "gkernel") { return FFType::GKernel; }
  if (type == "dirac") { return FFType::Dirac; }
  if (type == "logexp") { return FFType::LogExp; }
  throw std::invalid_argument("unknown form-factor type " + type);
}

// Parse one steering card form factor normalization point
FFNorm ParseFFNorm(const std::string &norm) {
  if (norm == "zero") { return FFNorm::Zero; }
  if (norm == "pole") { return FFNorm::Pole; }
  throw std::invalid_argument("unknown form-factor normalization " + norm);
}

// Read one minimal named form factor object
FFParam ReadFF(const nlohmann::json &block, const std::string &context) {
  if (!block.is_object() || !block.contains("type") || !block.at("type").is_string()) { throw std::invalid_argument(context + " must contain string type"); }

  FFParam ff;
  ff.type = ParseFFType(block.at("type").get<std::string>());
  if (ff.type == FFType::None) {
    CheckFields(block, {"type"}, context);
    return ff;
  }

  if (!block.contains("norm") || !block.at("norm").is_string()) { throw std::invalid_argument(context + " must contain string norm"); }
  ff.norm = ParseFFNorm(block.at("norm").get<std::string>());

  if (ff.type == FFType::Dirac) {
    CheckFields(block, {"type", "norm", "EM"}, context);
    ff.structure.EM = block.at("EM").get<std::string>();
  } else if (ff.type == FFType::Exponential) {
    CheckFields(block, {"type", "norm", "b"}, context);
    ff.param = {ReadNumber(block, "b", context)};
  } else if (ff.type == FFType::LogExp) {
    CheckFields(block, {"type", "norm", "b", "Lambda2"}, context);
    ff.param = {ReadNumber(block, "b", context), ReadNumber(block, "Lambda2", context)};
  } else if (ff.type == FFType::Power || ff.type == FFType::Vector) {
    CheckFields(block, {"type", "norm", "Lambda2", "n"}, context);
    ff.param = {ReadNumber(block, "Lambda2", context), ReadNumber(block, "n", context)};
  } else if (ff.type == FFType::Orear) {
    CheckFields(block, {"type", "norm", "b", "a"}, context);
    ff.param = {ReadNumber(block, "b", context), ReadNumber(block, "a", context)};
  } else if (ff.type == FFType::Gaussian) {
    CheckFields(block, {"type", "norm", "Lambda2"}, context);
    ff.param = {ReadNumber(block, "Lambda2", context)};
  } else {
    CheckFields(block, {"type", "norm", "terms"}, context);
    const auto &terms = block.at("terms");
    if (!terms.is_array() || terms.empty()) { throw std::invalid_argument(context + ".terms must be a nonempty array"); }
    for (const auto &i : indices(terms)) {
      const auto       &term         = terms.at(i);
      const std::string term_context = context + ".terms[" + std::to_string(i) + "]";
      CheckFields(term, {"a", "p", "nu", "mu2"}, term_context);
      ff.param.push_back(ReadNumber(term, "a", term_context));
      ff.param.push_back(ReadNumber(term, "p", term_context));
      ff.param.push_back(ReadNumber(term, "nu", term_context));
      ff.param.push_back(ReadNumber(term, "mu2", term_context));
    }
  }
  CheckFF(ff, context);
  return ff;
}

// Validate one positive generated invariant mass squared
void CheckMass2(const double mass2, const std::string &context) {
  if (!std::isfinite(mass2) || !(mass2 > 0.0)) { throw AmplitudeFailure(context + ": invalid generated mass squared"); }
}

// Validate one parsed form factor
void CheckFF(const FFParam &ff, const std::string &context) {
  if (!gra::AllFinite(ff.param)) { throw std::invalid_argument(context + " parameters must be finite"); }
  const auto require_size = [&](const std::size_t expected) {
    if (ff.param.size() != expected) { throw std::invalid_argument(context + " has an invalid parameter count"); }
  };
  if (ff.type == FFType::None) {
    require_size(0);
    return;
  }
  if (ff.type == FFType::Dirac) {
    require_size(0);
    if (ff.norm != FFNorm::Zero || (ff.structure.EM != "DIPOLE" && ff.structure.EM != "KELLY")) {
      throw std::invalid_argument(context + " dirac requires norm zero and EM DIPOLE or KELLY");
    }
    return;
  }
  if (ff.type == FFType::Exponential) {
    require_size(1);
    if (ff.param[0] < 0.0) { throw std::invalid_argument(context + " b must be nonnegative"); }
    return;
  }
  if (ff.type == FFType::LogExp) {
    require_size(2);
    if (ff.param[0] < 0.0 || !(ff.param[1] > 0.0)) {
      throw std::invalid_argument(context + " logexp b must be nonnegative and Lambda2 positive");
    }
    return;
  }
  if (ff.type == FFType::Power) {
    require_size(2);
    if (!(ff.param[0] > 0.0) || !(ff.param[1] > 0.0)) { throw std::invalid_argument(context + " Lambda2 and n must be positive"); }
    return;
  }
  if (ff.type == FFType::Orear) {
    require_size(2);
    if (ff.param[0] < 0.0 || ff.param[1] < 0.0) { throw std::invalid_argument(context + " b and a must be nonnegative"); }
    return;
  }
  if (ff.type == FFType::Gaussian) {
    require_size(1);
    if (!(ff.param[0] > 0.0)) { throw std::invalid_argument(context + " Lambda2 must be positive"); }
    return;
  }
  if (ff.type == FFType::Vector) {
    require_size(2);
    if (ff.norm != FFNorm::Pole || !(ff.param[0] > 0.0) || !(ff.param[1] > 0.0)) { throw std::invalid_argument(context + " vector Lambda2 and n must be positive with norm pole"); }
    return;
  }
  if (ff.type == FFType::GKernel) {
    if (ff.param.empty() || ff.param.size() % 4 != 0) { throw std::invalid_argument(context + " gkernel requires terms"); }
    for (std::size_t i = 0; i < ff.param.size(); i += 4) {
      if (ff.param[i] < 0.0 || !(ff.param[i + 1] > 0.0) || ff.param[i + 2] < 0.0 || ff.param[i + 3] < 0.0) { throw std::invalid_argument(context + " has invalid gkernel terms"); }
    }
    return;
  }
}

namespace {

// Evaluate one validated off-shell form-factor parameter block
// F(q2) is evaluated at x = q_ref^2-q2 with the selected normalized profile
double EvalFF(double q2, double pole2, const FFParam &ff) {
  if (!std::isfinite(q2) || (ff.norm == FFNorm::Pole && (!std::isfinite(pole2) || pole2 < 0.0))) { throw AmplitudeFailure("FormFactor: invalid virtuality or pole squared"); }
  if (ff.type == FFType::None) { return 1.0; }
  if (ff.type == FFType::Dirac) { return form::F1(q2, ff.structure); }
  const double x     = (ff.norm == FFNorm::Zero ? 0.0 : pole2) - q2;
  double       value = 0.0;

  if (ff.type == FFType::Exponential) {
    value = std::exp(-ff.param[0] * x);
  } else if (ff.type == FFType::LogExp) {
    value = LogExpFF(x, ff.param[0], ff.param[1]);
  } else if (ff.type == FFType::Power) {
    const double base = 1.0 + x / ff.param[0];
    if (base <= 0.0) { throw AmplitudeFailure("FormFactor: power event domain is nonpositive"); }
    value = std::pow(base, -ff.param[1]);
  } else if (ff.type == FFType::Orear) {
    const double b = ff.param[0];
    const double a = ff.param[1];
    if (x + pow2(a) < 0.0) { throw AmplitudeFailure("FormFactor: orear event domain is negative"); }
    value = std::exp(-b * (std::sqrt(x + pow2(a)) - a));
  } else if (ff.type == FFType::Gaussian) {
    value = std::exp(-pow2(x) / pow2(ff.param[0]));
  } else if (ff.type == FFType::Vector) {
    // [REFERENCE: arXiv:2508.06334v1, Eq. (2.29)]
    const double base = 1.0 + q2 * (q2 - pole2) / pow2(ff.param[0]);
    if (base <= 0.0) { throw AmplitudeFailure("FormFactor: vector event domain is nonpositive"); }
    value = std::pow(base, -ff.param[1]);
  } else {
    value = GKernel(x, ff.param, "FormFactor");
  }
  return value;
}

}  // namespace

// For data, see:
//
// [REFERENCE: Kirk for WA102, https://arxiv.org/abs/hep-ph/9908253v1]

// Evaluate one form factor at q^2 with its pole mass squared
double FormFactor(const double q2, const double pole2, const FFParam &ff) {
  return EvalFF(q2, pole2, ff);
}

// Evaluate one transfer form factor normalized at zero
double TransferFF(const double t, const FFParam &ff) { return EvalFF(t, 0.0, ff); }

// Evaluate one invariant mass form factor normalized at its pole
double MassFF(const double mass2, const double pole2, const FFParam &ff) { return EvalFF(mass2, pole2, ff); }

// Compute the analytic divided difference without subtracting nearly equal form factors
// The caller validates the none, exp or gaussian profile during initialization
std::array<double, 3> MassFFPair(const double q1, const double q2, const double pole2, const FFParam &ff) {
  if (ff.type == FFType::None) { return {1.0, 1.0, 0.0}; }
  const double ref = ff.norm == FFNorm::Pole ? pole2 : 0.0;
  const double x1 = q1 - ref, x2 = q2 - ref;
  const bool exponential = ff.type == FFType::Exponential;
  const double a = exponential ? ff.param[0] * x1 : -pow2(x1 / ff.param[0]);
  const double b = exponential ? ff.param[0] * x2 : -pow2(x2 / ff.param[0]);
  const double slope = exponential ? ff.param[0] : -(x1 + x2) / pow2(ff.param[0]);
  const double gap = std::abs(slope * (q1 - q2));
  const double ratio = std::fpclassify(gap) == FP_ZERO ? 1.0 : -std::expm1(-gap) / gap;
  return {std::exp(a), std::exp(b), slope * std::exp(std::max(a, b)) * ratio};
}

// Compute the off-shell hadron propagator with optional subchannel reggeization
// propagator = [t_hat-M2]^-1 exp[(alpha(t_hat)-J)L], L = 2 log cosh(Delta y/2)
// The smooth gap is a phenomenological variation of the absolute (hard) gap
// 
// [REFERENCE: https://arxiv.org/abs/1103.1203]
// [REFERENCE: Harland-Lang, Khoze, Ryskin, arxiv.org/abs/1312.4553]
//
std::complex<double> OffshellProp(double t_hat, double M2, const ReggeizeParam &reggeize, const MesonTraj &trajectory, const gra::M4Vec &left, const gra::M4Vec &right) {
  if (!std::isfinite(t_hat) || !std::isfinite(M2) || M2 < 0.0) { throw AmplitudeFailure("OffshellProp: invalid virtuality or nominal mass squared"); }

  const double denominator = t_hat - M2;
  if (std::fpclassify(denominator) == FP_ZERO) { throw AmplitudeFailure("OffshellProp: sampled the exchanged-meson pole"); }
  const std::complex<double> F = 1.0 / denominator;

  if (!reggeize.active) { return F; }

  // Anchor alpha(M^2)=J and freeze trajectory at the configured spacelike virtuality
  const double alpha = trajectory.spin + trajectory.ap * (std::max(-reggeize.freeze_scale2, t_hat) - M2);
  const double dY    = std::abs(left.Rap() - right.Rap());

  if (!std::isfinite(dY)) { throw AmplitudeFailure("OffshellProp: reggeized rapidity gap is not finite"); }
  
  return F * std::exp((alpha - trajectory.spin) * SmoothReggeInterval(dY));
}

// Evaluate one pair-specific exchanged-meson propagator
// [REFERENCE: Kycia et al., Phys. Rev. D 95 (2017) 094020, arXiv:1702.07572]
//
std::complex<double> MesonPropagator(const Param &param, double t_hat, double M2, const PairParam &entry, const gra::MDecayBranch &left, const gra::MDecayBranch &right) {
  const MesonTraj trajectory = !entry.reggeize.active ? MesonTraj{} : MesonTrajectory(param, left.p.pdg);
  return OffshellProp(t_hat, M2, entry.reggeize, trajectory, left.p4, right.p4);
}

// Evaluate the two vertex form factors of one exchanged-meson line
double MesonVertexFactor(const double t_hat, const double M2, const VertexParam &vertex) { return EvalFF(t_hat, M2, vertex.offshell[0]) * EvalFF(t_hat, M2, vertex.offshell[1]); }

// Evaluate one off-shell meson exchange without mixing in a Regge virtuality
std::complex<double> MesonExchange(const Param &param, double t_hat, double M2, const PairParam &entry, const VertexParam &vertex, const gra::MDecayBranch &left, const gra::MDecayBranch &right) {
  return MesonVertexFactor(t_hat, M2, vertex) * MesonPropagator(param, t_hat, M2, entry, left, right);
}

// Evaluate the charge and transfer factor of one ordered continuum pair vertex
double PairVertexFactor(const Param &param, const VertexParam &vertex, const gra::MDecayBranch &left, const gra::MDecayBranch &right, const double upper_t, const double lower_t, const bool upper_ff, const bool lower_ff) {
  if (!std::isfinite(upper_t) || !std::isfinite(lower_t)) { throw AmplitudeFailure("PairVertexFactor: non-finite exchange transfer"); }
  double factor = VertexSign(vertex, left, right, param);
  if (upper_ff && vertex.first != PDG::PDG_gamma) {
    factor *= TransferFF(upper_t, vertex.transfer[0]);
  }
  if (lower_ff && vertex.second != PDG::PDG_gamma) {
    factor *= TransferFF(lower_t, vertex.transfer[1]);
  }
  return factor;
}

// Evaluate one complete ordered continuum pair scalar factor
std::complex<double> PairExchange(const Param &param, const PairParam &entry, const VertexParam &vertex, const gra::MDecayBranch &left, const gra::MDecayBranch &right, const double t_hat,
                                  const double mass2, const double upper_t, const double lower_t) {
  const std::complex<double> value = MesonExchange(param, t_hat, mass2, entry, vertex, left, right) * PairVertexFactor(param, vertex, left, right, upper_t, lower_t, true, true);
  return value;
}

// Compute the amplitude factor for Poisson-distributed zero-secondary production
// [REFERENCE: Harland-Lang, Khoze, Ryskin, arXiv:1312.4553]
//
double Veto(const VetoParam &veto, double central_mass) {
  if (!veto.active) { return 1.0; }
  if (!std::isfinite(central_mass) || central_mass <= 0.0) { throw AmplitudeFailure("Veto: central mass must be finite and positive"); }
  if (central_mass <= veto.M0) { return 1.0; }

  return std::pow(veto.M0 / central_mass, veto.c);
}

}  // namespace gra::regge
