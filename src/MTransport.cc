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
  // Optimal Transport Matrix Pi                            (n x m)
  // Distance measure                                       (scalar)
  // 
	double SinkHorn(MMatrix<double>& Pi, const MMatrix<double>& C,
									std::vector<double>& p, std::vector<double>& q, double lambda, unsigned int iter) {
		
		const double SAFE_EPS = 1e-30;

		std::cout << "opt::SinkHorn optimization:" << std::endl;

		// ============================================================
		// Test the normalization
		auto fsum = [] (const std::vector<double>& x) {
			double value = 0.0;
			for (std::size_t i = 0; i < x.size(); ++i) { value += x[i]; }
			return value;
		};

		const double EPS = 1e-5;
		if (std::abs(fsum(p) - 1) > EPS) { 
			const std::string str = "SinkHorn:: Input p elements sum != 1";
			throw std::invalid_argument(str);
		}
		if (std::abs(fsum(q) - 1) > EPS) {
			const std::string str = "SinkHorn:: Input q elements sum != 1";
			throw std::invalid_argument(str);
		}
		// ============================================================

		const std::size_t n = C.size_row();
		const std::size_t m = C.size_col();
		
		// Gibbs Kernel: Entropy regularized distance matrix
		MMatrix<double> K(n,m,0.0);
		for (std::size_t i = 0; i < n; ++i) {
			for (std::size_t j = 0; j < m; ++j) {
				K[i][j] = std::exp(- C[i][j] / lambda);
			}
		}
		
		// Transposed version
		MMatrix<double> KT = K.Transpose();

		// Initialization
		std::vector<double> u(n, 1.0/n);
		std::vector<double> v(m, 1.0/m);

		// || Pi * 1 - p || or || Pi^T 1 - q ||
		std::vector<double> ones_u(n, 1.0);
		std::vector<double> ones_v(m, 1.0);

		auto resnorm = [&] (const MMatrix<double>& Pi) {
			double value = 0.0;
			const std::vector<double> lhs = Pi * ones_v;
			for (std::size_t i = 0; i < lhs.size(); ++i) {
				value += math::pow2(lhs[i] - p[i]);
			}
			return math::msqrt(value);
		};

		// Sinkhorn iterations
		double cost = 0.0;

		for (unsigned int k = 0; k < iter; ++k) {

			const std::vector<double> Kv  = K  * v; // >= 0
			for (std::size_t i = 0; i < n; ++i) { u[i] = p[i] /  std::max(SAFE_EPS, Kv[i]); }
			
			const std::vector<double> KTu = KT * u; // >= 0
			for (std::size_t i = 0; i < m; ++i) { v[i] = q[i] / std::max(SAFE_EPS, KTu[i]); }
			
			// Calculate metrics
			if ((k+1) % (iter / 10) == 0 || k == 0 || k == iter - 1) {
			
				// Transport coupling matrix
				Pi = matoper::diagAdiag(u, K, v);

				// Transport cost
				cost = resnorm(Pi);

				printf("iter = %3d / %3d : cost = %0.5E, log(cost) = %0.1f \n", k+1, iter, cost, std::log(cost));
			}
		}
		
		return cost;
	}

} // opt namespace
} // gra namespace
