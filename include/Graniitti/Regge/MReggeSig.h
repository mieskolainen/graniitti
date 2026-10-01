// Regge eta factors and steering modes
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MREGGESIG_H
#define MREGGESIG_H

// C++
#include <complex>
#include <string>

namespace gra {

// Select one complete Regge eta factor prescription
enum class EtaMode { RotatingT0, Rotating, Raw };

namespace regge {

// Select the positive or negative Regge signature
enum class Signature { Positive, Negative };

// Select the crossed energy boundary value of the s channel cut
enum class Rim { Lower, Upper };

// Classify an integer angular momentum under one Regge signature
enum class IntegerType { Allowed, Opposite };

// Parse one exact integer Regge signature
Signature ParseSignature(int tau, const std::string &path);

// Compute tau = +1 or -1 for one Regge signature
int Tau(Signature signature);

// Classify the signature at one integer angular momentum
IntegerType Classify(int J, Signature signature);

// Evaluate the full unregulated Regge signature factor
std::complex<double> EtaRaw(double J, Signature signature, Rim rim);

// Evaluate the full unregulated signature factor at complex angular momentum
std::complex<double> EtaRaw(std::complex<double> J, Signature signature, Rim rim);

// Evaluate the full Regge signature factor with regulated physical poles
std::complex<double> Eta(double alpha_t, Signature signature, double pole_epsilon = 1e-8);

// Evaluate the pole stripped unit modulus rotating signature phase
std::complex<double> EtaPhase(double alpha_t, Signature signature);

// Evaluate the positive s power using its real logarithm
std::complex<double> Power(double s, double s0, std::complex<double> J);

// Evaluate one signatured four point Regge pole kernel
// The trajectory vertex and particle pole normalization remain external
std::complex<double> Pole(double s, double s0, std::complex<double> J, Signature signature, Rim rim);

// Parse one exact Regge eta factor mode
EtaMode ParseEta(const std::string &mode, const std::string &path);

// Compute the steering card name of one Regge eta factor mode
std::string EtaName(EtaMode mode);

// Validate one trajectory and eta factor prescription
void CheckEta(double alpha0, double alpha_prime, Signature signature, EtaMode mode, const std::string &path);

// Evaluate one complete Regge eta factor
std::complex<double> EtaFactor(double alpha_t, double alpha0, Signature signature, EtaMode mode);

}  // namespace regge
}  // namespace gra

#endif  // MREGGESIG_H
