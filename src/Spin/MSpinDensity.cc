// Spin density and bipartite entanglement diagnostics
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Spin/MSpinDensity.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

#include "Graniitti/Spin/MHelicity.h"
#include "Graniitti/Tech/MAux.h"

using gra::aux::indices;

namespace gra::spin {

namespace {

// Estimate accumulated floating point roundoff for a normalized density matrix
double DensityTolerance(const std::size_t dimension) {
  const double n = static_cast<double>(dimension);
  return n * n * std::numeric_limits<double>::epsilon();
}

// Normalize a finite spin density without changing its relative matrix elements
MMatrix<std::complex<double>> NormalizeDensity(const MMatrix<std::complex<double>>& rho) {
  if (rho.isEmpty() || !rho.IsFinite()) { throw std::invalid_argument("Spin density is empty or non-finite"); }
  const auto trace = rho.Trace();
  if (!(trace.real() > 0.0) || !std::isfinite(trace.real()) ||
      std::abs(trace.imag()) > DensityTolerance(rho.size_row()) * trace.real()) {
    throw std::invalid_argument("Spin density must have a real positive trace");
  }
  return rho / trace.real();
}

}  // namespace

// Compute purity and von Neumann entropy in bits for a normalized spin density
SpinMetrics DensityMetrics(const MMatrix<std::complex<double>>& rho) {
  const auto   normalized = NormalizeDensity(rho);
  const double entropy    = normalized.SelfAdjointEntropy(DensityTolerance(rho.size_row())) / std::log(2.0);
  return {normalized.FrobNorm2(), std::max(0.0, entropy)};
}

// Compute N = (||rho^T2||_1 - 1)/2 and E_N = log2(1 + 2N)
// [REFERENCE: Vidal and Werner, Phys. Rev. A 65 (2002) 032314, arXiv:quant-ph/0102117]
PairMetrics EntanglementMetrics(const MMatrix<std::complex<double>>& rho, std::size_t n1, std::size_t n2) {
  const auto  normalized = NormalizeDensity(rho);
  PairMetrics out;
  out.pair               = DensityMetrics(normalized);
  out.entropy1           = DensityMetrics(normalized.PartialTrace(n1, n2, TensorFactor::Second)).entropy;
  out.entropy2           = DensityMetrics(normalized.PartialTrace(n1, n2, TensorFactor::First)).entropy;
  const double tolerance = DensityTolerance(rho.size_row());
  for (const double eigenvalue :
       normalized.PartialTranspose(n1, n2, TensorFactor::Second).SelfAdjointEigenvalues(tolerance)) {
    if (eigenvalue < -tolerance) { out.negativity -= eigenvalue; }
  }
  out.log_negativity = std::log2(1.0 + 2.0 * out.negativity);
  return out;
}

// Embed physical helicity rows in the full product basis before tracing the mother spin
MMatrix<std::complex<double>> DecayDensity(const HELMatrix& hel, const MMatrix<std::complex<double>>& parent,
                                           double theta, double phi) {
  (void)DensityMetrics(parent);
  const auto               decay = fDecayMatrix(hel, theta, phi);
  std::vector<std::size_t> rows(hel.lambda_idx.size_row());
  for (const auto i : indices(rows)) { rows[i] = hel.lambda_idx[i][0] * hel.T.size_col() + hel.lambda_idx[i][1]; }
  MMatrix<std::complex<double>> amplitude(hel.T.Elements().size(), decay.size_col(), 0.0);
  amplitude.SetRows(rows, decay);
  auto         rho  = amplitude * parent * amplitude.Dagger();
  const double norm = std::real(rho.Trace());
  if (!(norm > 0.0) || !std::isfinite(norm)) {
    throw std::invalid_argument("DecayDensity: zero or invalid decay probability at the selected angles");
  }
  return rho / norm;
}

// Print conditional pair metrics with the reference state and entropy units explicit
void PrintDecayEntanglement(const HELMatrix& hel) {
  const std::size_t n1 = hel.T.size_row(), n2 = hel.T.size_col();
  if (n1 <= 1 || n2 <= 1) { return; }
  const std::size_t n       = hel.Jz_values.size();
  const auto        parent  = MMatrix<std::complex<double>>::IdentityMatrix(n) / static_cast<double>(n);
  const auto        metrics = EntanglementMetrics(DecayDensity(hel, parent, 0.0, 0.0), n1, n2);
  std::cout << "Decay spin entanglement: reduced vertex, unpolarized mother, fixed daughter direction" << std::endl;
  gra::aux::PrintTable({"Tr(rho12^2)", "S12 [bits]", "S1 [bits]", "S2 [bits]", "Negativity", "Log negativity [bits]"},
                       {{gra::aux::ToString(metrics.pair.purity, 4), gra::aux::ToString(metrics.pair.entropy, 4),
                         gra::aux::ToString(metrics.entropy1, 4), gra::aux::ToString(metrics.entropy2, 4),
                         gra::aux::ToString(metrics.negativity, 4), gra::aux::ToString(metrics.log_negativity, 4)}});
  std::cout << "N > 0 certifies entanglement. S1 and S2 are entanglement entropies only for a pure pair" << std::endl;
}

}  // namespace gra::spin
