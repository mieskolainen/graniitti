// Optimal Transport Algorithms
//
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++ standard
#include <algorithm>
#include <complex>
#include <iostream>
#include <random>
#include <vector>

// Own
#include "Graniitti/MTransport.h"
#include "Graniitti/MMatrix.h"
#include "Graniitti/MMatOper.h"

// Eigen
#include <Eigen/Dense>

namespace gra {
namespace opt {

  // Sinkhorn-Knopp non-linear but convex optimization algorithm
  // for Optimal Transport problems
  // 
  // Input:
  // 
  // Cost matrix C between elements of p and q              (n x m)   e.g. ||x_i - x_j||^2
  // Probability density p (histogram / dirac point mass)   (n x 1)   with \sum = 1
  // Probability density q (histogram / dirac point mass)   (m x 1)   with \sum = 1
  // Entropic regularization lambda                         (scalar)
  // Iterations                                             (scalar)
  // 
  // Output:
  // 
  // Optimal Transport Matrix gamma                         (n x m)
  // Distance                                               (scalar)
  // 
	double SinkHorn(MMatrix<double>& gamma, const MMatrix<double>& C,
									std::vector<double>& p, std::vector<double>& q, double lambda, unsigned int iter) {
		
		const std::size_t n = C.size_row();
		const std::size_t m = C.size_col();
		
		// Calculate entropy regularized distance matrix
		MMatrix<double> K(n,m,0.0);
		for (std::size_t i = 0; i < n; ++i) {
			for (std::size_t j = 0; j < m; ++j) {
				K[i][j] = std::exp(- C[i][j] / lambda);
			}
		}
		
		// Transposed version
		MMatrix<double> KT = K; KT.Transpose();
		
		// Sinkhorn iterations
		std::vector<double> u(n, 1.0/n);
		std::vector<double> v(m, 1.0/m);
		
		for (unsigned int k = 0; k < iter; ++k) {
			
			const std::vector<double> Kv  = K  * v;
			for (std::size_t i = 0; i < n; ++i) { u[i] = p[i] /  Kv[i]; } // Possible division by zero!
			
			const std::vector<double> KTu = KT * u;
			for (std::size_t i = 0; i < m; ++i) { v[i] = q[i] / KTu[i]; } // Possible division by zero!
		}
		
		// Transport coupling matrix
		gamma = matoper::OuterProd(v,u);
		
		// Transport cost
		const double cost = gamma.Sum();
		
		return cost;
	}

} // opt namespace
} // gra namespace
