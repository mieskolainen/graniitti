// Optimal Transport Class
//
//
// (c) 2017-2019 Mikael Mieskolainen
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

  double SinkHorn(MMatrix<double>& P, const MMatrix<double>& C,
                  std::vector<double>& p, std::vector<double>& q, double lambda, unsigned int iter);

}  // math namespace ends
}  // gra namespace ends

#endif