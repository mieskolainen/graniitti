// Entropy regularized optimal transport algorithms
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MTRANSPORT_H
#define MTRANSPORT_H

#include <cstddef>
#include <vector>

#include "Graniitti/Math/MMatrix.h"

namespace gra {
namespace opt {

// Construct a Gaussian convolution kernel on two regular one-dimensional grids
void ConvKernel(std::size_t n, std::size_t m, double lambda, MMatrix<double> &K);

// Construct an entropy-regularized Gibbs kernel from a nonnegative cost matrix
void GibbsKernel(double lambda, const MMatrix<double> &C, MMatrix<double> &K);

// Compute the final source-marginal L1 residual after Sinkhorn iteration
double SinkHorn(MMatrix<double> &Pi, const MMatrix<double> &K, const std::vector<double> &p,
                const std::vector<double> &q, std::size_t iter);

}  // namespace opt
}  // namespace gra

#endif
