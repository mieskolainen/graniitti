// Regge eta (signature) factors and steering modes
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <cmath>
#include <complex>
#include <stdexcept>
#include <string>

// Own
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Regge/MReggeSig.h"
#include "Graniitti/Tech/MException.h"

namespace gra {
namespace regge {

namespace {

struct HalfAngle {
  double sine;
  double cosine;
};

// Validate one finite complex angular momentum
void CheckJ(const std::complex<double> J) {
  if (!std::isfinite(J.real()) || !std::isfinite(J.imag())) { throw std::invalid_argument("gra::regge: angular momentum J must be finite"); }
}

// Evaluate sin and cos of pi alpha over two with exact integer values
HalfAngle EvaluateHalfAngle(const double alpha_t) {
  const double reduced_alpha = std::remainder(alpha_t, 4.0);
  const double half_phase    = math::PI * reduced_alpha / 2.0;
  HalfAngle    value{std::sin(half_phase), std::cos(half_phase)};

  const double nearest_integer = std::nearbyint(reduced_alpha);
  if (std::fpclassify(reduced_alpha - nearest_integer) == FP_ZERO) {
    const int residue = (static_cast<int>(nearest_integer) % 4 + 4) % 4;
    if (residue == 0) {
      value = {0.0, 1.0};
    } else if (residue == 1) {
      value = {1.0, 0.0};
    } else if (residue == 2) {
      value = {0.0, -1.0};
    } else {
      value = {-1.0, 0.0};
    }
  }
  return value;
}

// Compute the rim dependent constant term of the signature factor
std::complex<double> RimTerm(const Signature signature, const Rim rim) {
  const double tau = signature == Signature::Positive ? 1.0 : -1.0;
  return rim == Rim::Lower ? math::zi * tau : -math::zi * tau;
}

}  // namespace

// Parse one exact integer Regge signature
Signature ParseSignature(const int tau, const std::string &path) {
  if (tau == 1) { return Signature::Positive; }
  if (tau == -1) { return Signature::Negative; }
  throw std::invalid_argument(path + " must be -1 or 1");
}

// Compute tau = +1 or -1 for one Regge signature
int Tau(const Signature signature) {
  if (signature == Signature::Positive) { return 1; }
  if (signature == Signature::Negative) { return -1; }
  throw std::invalid_argument("gra::regge::Tau: invalid Regge signature");
}

// Classify the signature at one integer angular momentum
IntegerType Classify(const int J, const Signature signature) {
  Tau(signature);
  const Signature integer_signature = J % 2 == 0 ? Signature::Positive : Signature::Negative;
  return signature == integer_signature ? IntegerType::Allowed : IntegerType::Opposite;
}

// Evaluate the full unregulated Regge signature factor
// eta_+ = -cot(pi J/2)+i, eta_- = -tan(pi J/2)-i on the lower rim
std::complex<double> EtaRaw(const double J, const Signature signature, const Rim rim) {

  const HalfAngle angle       = EvaluateHalfAngle(J);
  const double    meromorphic = signature == Signature::Positive ? -angle.cosine / angle.sine : -angle.sine / angle.cosine;
  return meromorphic + RimTerm(signature, rim);
}

// Evaluate the full unregulated signature factor at complex angular momentum
std::complex<double> EtaRaw(const std::complex<double> J, const Signature signature, const Rim rim) {
  if (std::fpclassify(J.imag()) == FP_ZERO) { return EtaRaw(J.real(), signature, rim); }

  const bool                 positive = signature == Signature::Positive;
  const std::complex<double> reduced(std::remainder(J.real(), 2.0), J.imag());
  if (J.imag() > 0.0) {
    const std::complex<double> q = std::exp(math::zi * math::PI * reduced);
    if (rim == Rim::Lower) { return positive ? 2.0 * math::zi / (1.0 - q) : -2.0 * math::zi / (1.0 + q); }
    return positive ? 2.0 * math::zi * q / (1.0 - q) : 2.0 * math::zi * q / (1.0 + q);
  }

  const std::complex<double> r = std::exp(-math::zi * math::PI * reduced);
  if (rim == Rim::Lower) { return positive ? -2.0 * math::zi * r / (1.0 - r) : -2.0 * math::zi * r / (1.0 + r); }
  return positive ? -2.0 * math::zi / (1.0 - r) : 2.0 * math::zi / (1.0 + r);
}

// Evaluate the full Regge signature factor with regulated physical poles
std::complex<double> Eta(const double alpha_t, const Signature signature, const double pole_epsilon) {
  CheckJ({alpha_t, 0.0});
  Tau(signature);
  if (!std::isfinite(pole_epsilon) || pole_epsilon <= 0.0) { throw std::invalid_argument("gra::regge::Eta: pole epsilon must be positive and finite"); }

  const HalfAngle            angle = EvaluateHalfAngle(alpha_t);
  const std::complex<double> phase(angle.cosine, -angle.sine);
  const double               regulator = std::tanh(pole_epsilon / 2.0);

  if (signature == Signature::Positive) {
    const std::complex<double> denominator(angle.sine, angle.cosine * regulator);
    return -phase / denominator;
  }
  const std::complex<double> denominator(angle.cosine, -angle.sine * regulator);
  return -math::zi * phase / denominator;
}

// Evaluate the pole-stripped unit-modulus rotating signature phase
// eta_+ = -exp(-i pi alpha/2), eta_- = -i exp(-i pi alpha/2)
std::complex<double> EtaPhase(const double alpha_t, const Signature signature) {

  const HalfAngle            angle = EvaluateHalfAngle(alpha_t);
  const std::complex<double> phase(angle.cosine, -angle.sine);
  const std::complex<double> factor = signature == Signature::Positive ? -1.0 : -math::zi;
  return factor * phase;
}

// Evaluate the positive energy power using its real logarithm
// power = (s/s0)^J
std::complex<double> Power(const double s, const double s0, const std::complex<double> J) {
  if (!std::isfinite(s) || !(s > 0.0)) { throw AmplitudeFailure("gra::regge::Power: invalid generated subenergy"); }
  return std::exp(J * (std::log(s) - std::log(s0)));
}

// Evaluate one signatured four-point Regge-pole kernel
// P_tau(s,J) = eta_tau(J)(s/s0)^J
// The trajectory vertex and particle pole normalization remain external
std::complex<double> Pole(const double s, const double s0, const std::complex<double> J, const Signature signature, const Rim rim) { return EtaRaw(J, signature, rim) * Power(s, s0, J); }

// Parse one exact Regge eta factor mode
EtaMode ParseEta(const std::string &mode, const std::string &path) {
  if (mode == "raw") { return EtaMode::Raw; }
  if (mode == "rotating_t0") { return EtaMode::RotatingT0; }
  if (mode == "rotating") { return EtaMode::Rotating; }
  throw std::invalid_argument(path + " has unknown mode " + mode);
}

// Compute the steering card name of one Regge eta factor mode
std::string EtaName(const EtaMode mode) {
  if (mode == EtaMode::Raw) { return "raw"; }
  if (mode == EtaMode::RotatingT0) { return "rotating_t0"; }
  if (mode == EtaMode::Rotating) { return "rotating"; }
  throw std::invalid_argument("gra::regge::EtaName: invalid mode");
}

// Validate one trajectory and eta factor prescription
void CheckEta(const double alpha0, const double alpha_prime, const Signature signature, const EtaMode mode, const std::string &path) {
  CheckJ({alpha0, 0.0});
  Tau(signature);
  if (!std::isfinite(alpha_prime)) { throw std::invalid_argument(path + " alpha_prime must be finite"); }
  if (mode != EtaMode::Raw && mode != EtaMode::RotatingT0 && mode != EtaMode::Rotating) { throw std::invalid_argument(path + " has an invalid eta mode"); }
}

// Evaluate one complete Regge eta factor
std::complex<double> EtaFactor(const double alpha_t, const double alpha0, const Signature signature, const EtaMode mode) {
  if (!std::isfinite(alpha_t)) { throw AmplitudeFailure("MReggeSig::EtaFactor: non-finite generated trajectory"); }

  std::complex<double> value = 0.0;
  if (mode == EtaMode::Raw) {
    value = EtaRaw(alpha_t, signature, Rim::Lower);
  } else if (mode == EtaMode::RotatingT0) {
    value = EtaPhase(alpha0, signature);
  } else if (mode == EtaMode::Rotating) {
    value = EtaPhase(alpha_t, signature);
  } else {
    throw std::invalid_argument("MReggeSig::EtaFactor: invalid eta mode");
  }

  return value;
}

}  // namespace regge
}  // namespace gra
