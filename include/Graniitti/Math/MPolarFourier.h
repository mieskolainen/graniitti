// Non-radial polar Fourier-Bessel transforms and convolutions
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MPOLARFOURIER_H
#define MPOLARFOURIER_H

#include <complex>
#include <cstddef>
#include <vector>

#include "Graniitti/Math/MMatrix.h"

namespace gra::math {

// Store radial nodes and weights which discretize the measure r dr
struct PolarMeasureRule {
  std::vector<double> node;
  std::vector<double> measure_weight;
};

// Store angular Fourier coefficients at every impact parameter node
struct PolarHarmonicField {
  std::size_t max_harmonic = 0;
  MMatrix<std::complex<double>> coefficient;
};

// Store one prepared inverse Fourier-Bessel quadrature kernel
struct PolarInverseKernel {
  std::size_t max_harmonic = 0;
  MMatrix<double> radial_weight;
  std::vector<std::complex<double>> angular_phase;
};

// Compute one signed harmonic from a periodic discrete Fourier grid
int SignedHarmonic(unsigned int mode, unsigned int count);

// Compute a finite radial order-zero Fourier-Bessel transform
double RadialBesselTransform0(const std::vector<double> &node,
                              const std::vector<double> &weight,
                              const std::vector<double> &field, double momentum,
                              double coordinate_scale);

// Prepare the exact Hankel transform of a field linear in log radius
// The field is constant below the first node and zero at and beyond the last node
MMatrix<double> LogBesselKernel(const std::vector<double>& radius, const std::vector<double>& momentum,
                                double coordinate_scale);

// Construct a linear Gauss-Legendre rule including the radial measure r dr
PolarMeasureRule GaussLegendrePolarMeasure(unsigned int node_count,
                                           double minimum, double maximum);

// Transform non-radial complex fields between polar momentum and impact space
class MPolarFourier {
public:
  // Initialize fixed momentum, impact parameter, and angular quadrature grids
  MPolarFourier(PolarMeasureRule momentum_rule, PolarMeasureRule impact_rule,
                std::vector<double> azimuth_node, std::size_t max_harmonic);

  // Compute the momentum radial quadrature rule
  const PolarMeasureRule &MomentumRule() const;

  // Compute the impact parameter radial quadrature rule
  const PolarMeasureRule &ImpactRule() const;

  // Compute the uniform angular quadrature nodes
  const std::vector<double> &AzimuthNodes() const;

  // Compute the largest retained harmonic of each input field
  std::size_t MaxInputHarmonic() const;

  // Apply F(b)=1/(2 pi) int d2k exp(+i k.b) f(k)
  PolarHarmonicField
  Forward(const MMatrix<std::complex<double>> &momentum_field) const;

  // Multiply two impact parameter fields by discrete angular convolution
  PolarHarmonicField Multiply(const PolarHarmonicField &first,
                              const PolarHarmonicField &second) const;

  // Multiply three impact parameter fields by discrete angular convolution
  PolarHarmonicField Multiply(const PolarHarmonicField &first,
                              const PolarHarmonicField &second,
                              const PolarHarmonicField &third) const;

  // Prepare the inverse kernel once for a fixed output momentum
  PolarInverseKernel PrepareInverseKernel(double q, double azimuth,
                                          std::size_t max_harmonic) const;

  // Evaluate one impact parameter field with a prepared inverse kernel
  std::complex<double> InverseEvaluate(const PolarHarmonicField &impact_field,
                                       const PolarInverseKernel &kernel) const;

  // Evaluate several impact parameter fields with one prepared inverse kernel
  std::vector<std::complex<double>>
  InverseEvaluateBatch(const std::vector<PolarHarmonicField> &impact_field,
                       const PolarInverseKernel &kernel) const;

  // Evaluate f(q)=1/(2 pi) int d2b exp(-i q.b) F(b) at arbitrary q
  std::complex<double> InverseEvaluate(const PolarHarmonicField &impact_field,
                                       double q, double azimuth) const;

  // Evaluate the direct two-field transverse convolution at arbitrary q
  std::complex<double> ConvolutionEvaluate(const PolarHarmonicField &first,
                                           const PolarHarmonicField &second,
                                           double q, double azimuth) const;

  // Evaluate the direct three-field transverse convolution at arbitrary q
  std::complex<double> ConvolutionEvaluate(const PolarHarmonicField &first,
                                           const PolarHarmonicField &second,
                                           const PolarHarmonicField &third,
                                           double q, double azimuth) const;

private:
  PolarMeasureRule momentum_rule_;
  PolarMeasureRule impact_rule_;
  std::vector<double> azimuth_node_;
  std::size_t max_harmonic_ = 0;
  MMatrix<std::complex<double>> angular_kernel_;
  std::vector<MMatrix<double>> forward_kernel_;
};

// Build one transform from an existing polar momentum rule
MPolarFourier BuildPolarFourierTransform(PolarMeasureRule momentum_rule,
                                         unsigned int azimuth_count,
                                         unsigned int impact_count,
                                         double impact_max,
                                         std::size_t max_harmonic);

} // namespace gra::math

#endif
