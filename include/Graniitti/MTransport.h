// Optimal Transport Class
//
//
// (c) 2017-2021 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MTRANSPORT_H
#define MTRANSPORT_H

// C++
#include <algorithm>
#include <complex>
#include <random>
#include <stdexcept>
#include <vector>

// Eigen
#include <Eigen/Dense>

// Own
#include "Graniitti/MMatrix.h"
#include "Graniitti/MTensor.h"

namespace gra {
namespace opt {

void ConvKernel(std::size_t n, std::size_t m, double lambda, MMatrix<double>& K);
void GibbsKernel(double lambda, const MMatrix<double>& C, MMatrix<double>& K);

double SinkHorn(MMatrix<double>& P, const MMatrix<double>& K, std::vector<double>& p,
                std::vector<double>& q, std::size_t iter);

}  // namespace opt
}  // namespace gra

#endif