// GRANIITTI - Monte Carlo event generator for high energy diffraction
// https://github.com/mieskolainen/graniitti
//
// <NeuroJacobian neural net prototype>
//
//
// (c) 2017-2020 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MNEUROJACOBIAN_H
#define MNEUROJACOBIAN_H

// C++ includes
#include <future>
#include <iostream>

// Own
#include <Graniitti/MAux.h>
#include <Graniitti/MH2.h>
#include <Graniitti/MMath.h>
#include <Graniitti/MProcess.h>
#include <Graniitti/MRandom.h>

// Eigen includes
#include <Eigen/Core>

// L-BFGS algorithms
#include <LBFGS/include/LBFGS.h>

// autodiff library
#include <autodiff/forward.hpp>
#include <autodiff/forward/eigen.hpp>
#include <autodiff/reverse.hpp>
#include <autodiff/reverse/eigen.hpp>

using namespace gra;
using gra::aux::indices;
using gra::math::PI;
using gra::math::pow2;

using Eigen::VectorXd;
using namespace Eigen;
using namespace autodiff;
using namespace LBFGSpp;

namespace gra {
namespace neurojac {
MRandom      randx;
unsigned int BATCHSIZE   = 0;
bool         FLAT_TARGET = false;

// declare a pointer to member function
MProcess *procptr;
double (MProcess::*ptfptr)(const std::vector<double> &randvec,
                           AuxIntData &               aux) = &MProcess::EventWeight;

// Function to be integrated
double func(std::vector<double> &x) {
  if (FLAT_TARGET) { return 1.0; }

  gra::AuxIntData aux;
  aux.vegasweight  = 1.0;
  aux.burn_in_mode = false;
  double y         = (gra::neurojac::procptr->*gra::neurojac::ptfptr)(x, aux);

  // Numerical protection
  if (std::isnan(y) || std::isinf(y)) { y = 0.0; }
  return y;
}

// Exponential linear unit
inline dual elu(dual z) {
  if (z >= 0.0) {
    return z;
  } else {
    return exp(z) - 1.0;
  }
}

// Hyperbolic functions
inline dual hypersinh(dual z) { return (exp(z) - exp(-z)) / 2; }
inline dual hypertanh(dual z) { return (exp(2.0 * z) - 1.0) / (exp(2.0 * z) + 1.0); }

// Sigmoid
inline dual sigmoid(dual z) { return 1.0 / (1.0 + exp(-z)); }
// Modified Sigmoid
inline dual logexp(dual z) {
  const double alpha = 25.0;
  return 1.0 / alpha * log((1.0 + exp(alpha * z)) / (1.0 + exp(alpha * (z - 1.0))));
}

// Gaussian
dual gaussprob(const VectorXdual &z, double mu, double sigma) {
  dual prod = 1.0;
  for (int i = 0; i < z.size(); ++i) {
    prod *= (1.0 / sqrt(2 * PI * sigma * sigma)) * exp(-(z[i] - mu) * (z[i] - mu) / 2.0);
  }
  return prod;
}

// Feedforward Network layer
class Layer {
 public:
  Layer(int M, int N) {
    w = MatrixXdual(M, N);
    b = VectorXdual(M);
  }
  ~Layer() {}

  MatrixXdual w;
  VectorXdual b;
};

// Parameter struct
struct NetParams {
  std::vector<Layer> L;

  // Return number of parameters
  int size_param() {
    int k = 0;
    for (std::size_t l = 0; l < L.size(); ++l) { k += L[l].w.rows() * L[l].w.cols(); }
    for (std::size_t l = 0; l < L.size(); ++l) { k += L[l].b.size(); }
    return k;
  }

  int  D = 0;  // Integral dimension
  dual alpha;
};

// Global parameters
NetParams par;

// Neural network function G: D^N -> D^N, D = [0 ... 1]
//
// This depends on global parameter struct p
//
//
VectorXdual G_net(const VectorXdual &x) {
  const double beta = 0.5;

  // Network input
  VectorXdual in(x.size());
  for (std::size_t i = 0; i < (unsigned int)x.size(); ++i) { in[i] = x[i]; }

  // Network layers
  for (std::size_t l = 0; l < par.L.size() - 1; ++l) {
    VectorXdual out(par.L[l].w.rows());
    for (std::size_t i = 0; i < (unsigned int)par.L[l].w.rows(); ++i) {
      for (std::size_t j = 0; j < (unsigned int)par.L[l].w.cols(); ++j) {
        out[i] += par.L[l].w(i, j) * in[j];
      }
      out[i] = hypertanh(out[i] + par.L[l].b[i]) * beta + (1.0 - beta) * out[i];
    }
    in = out;  // Next layer input is this layer output
  }

  // Output layer compression to [0,1] range
  const unsigned int last = par.L.size() - 1;
  VectorXdual        out(par.L[last].w.rows());
  for (std::size_t i = 0; i < (unsigned int)par.L[last].w.rows(); ++i) {
    for (std::size_t j = 0; j < (unsigned int)par.L[last].w.cols(); ++j) {
      out[i] += par.L[last].w(i, j) * in[j];
    }
    out[i] = logexp(out[i] + par.L[last].b[i]);
  }
  return out;
}

class MNeuroJacobian {
 public:
  MNeuroJacobian() {}
  ~MNeuroJacobian() {}

  // Vector format to parameters
  static void Vec2Par(const std::vector<double> &w, NetParams par) {
    std::size_t k = 0;
    for (std::size_t l = 0; l < par.L.size(); ++l) {
      for (std::size_t i = 0; i < (unsigned int)par.L[l].w.rows(); ++i) {
        for (std::size_t j = 0; j < (unsigned int)par.L[l].w.cols(); ++j) {
          par.L[l].w(i, j) = w[k];
          ++k;
        }
      }
    }

    for (std::size_t l = 0; l < par.L.size(); ++l) {
      for (std::size_t j = 0; j < (unsigned int)par.L[l].b.size(); ++j) {
        par.L[l].b(j) = w[k];
        ++k;
      }
    }
  }

  // Random init for parameter vector
  std::vector<double> RandomInit(NetParams &par, double sigma) {
    std::vector<double> w(par.size_param());
    unsigned int        k = 0;

    // Connection matrices with random initialization
    for (std::size_t l = 0; l < par.L.size(); ++l) {
      for (std::size_t i = 0; i < (unsigned int)par.L[l].w.rows(); ++i) {
        for (std::size_t j = 0; j < (unsigned int)par.L[l].w.cols(); ++j) {
          w[k] = randx.G(0, sigma) * sqrt(2.0 / par.L[l].w.cols());  // He-et all style
          ++k;
        }
      }
    }
    // Bias initialized to zero
    for (std::size_t l = 0; l < par.L.size(); ++l) {
      for (std::size_t i = 0; i < (unsigned int)par.L[l].b.size(); ++i) {
        w[k] = 0.0;
        ++k;
      }
    }
    return w;
  }

  // Kullback-Leibler divergence loss function
  //
  //
  // For KL-divergence and algorithms with Jacobians, see:
  // [REFERENCE: Leow et al., Statistical Properties of Jacobian Maps ..., 2007]
  // http://www-sop.inria.fr/asclepios/cours/MVA/Module2/Papers/Leow_TMI07_LogUnbiased.pdf
  //
  // For neural networks and Jacobians, see:
  // [REFERENCE: Dillon et al, Tensorflow Bijector API,
  // https://arxiv.org/abs/1711.10604v1]
  //
  // Matrix derivative results, see:
  // [REFERENCE: M. Giles, http://eprints.maths.ox.ac.uk/1079/1/NA-08-01.pdf]
  //
  // MC Integration:
  // [REFERENCE: J. Bendavid, https://arxiv.org/pdf/1707.00028.pdf]
  //
  //
  static double singleloss(VectorXdual &z, const NetParams &par) {
    // Input likelihood value
    double p = val(gaussprob(z, 0, 1));

    // Evaluate network
    VectorXdual u = G_net(z);

    // Jacobian matrix du/dz [AUTODIFF]
    MatrixXdual J = jacobian(G_net, u, z);

    // Absolute Jacobian determinant
    double absDetJ = abs(val(J.determinant()));

    // ==============================================================
    // Evaluate integrand function
    std::vector<double> u_(u.size());
    for (std::size_t i = 0; i < (unsigned int)u.size(); ++i) { u_[i] = val(u[i]); }
    double fG = func(u_);
    // ==============================================================

    // Prior term
    double KL = 0.0;
    if (p > 0) { KL = log(p); }

    // Jacobian determinant term
    if (absDetJ > 0) { KL = KL - log(absDetJ); }

    // Target function fidelity term
    if (fG > 0) { KL = KL - log(fG); }

    return KL;
  }

  static dual loss(VectorXdual &Z, const NetParams &par) {
    VectorXdual z(par.D);
    dual        sumloss = 0.0;

    // Loop over all samples
    int k = 0;
    while (true) {
      // Pick prior p(z) distribution samples
      for (std::size_t i = 0; i < (unsigned int)par.D; ++i) {
        z[i] = Z[k];
        ++k;
      }
      sumloss = sumloss + singleloss(z, par);
      if (k == Z.size()) break;
    }
    return sumloss;
  }

  // Evaluate loss and gradient vector
  static double CostGrad(const std::vector<double> &w, std::vector<double> &gradv) {
    // Generate new virtual batch from noise distribution (prior)
    VectorXdual Z(BATCHSIZE * par.D);
    for (std::size_t i = 0; i < (unsigned int)Z.size(); ++i) { Z[i] = randx.G(0, 1); }

    // Update vector to global parameters
    Vec2Par(w, par);

    // ------------------------------------------------------------------
    // Evaluate gradient vector components numerically (not very optimal..)
    std::size_t k = 0;
    gradv.resize(w.size());

    // Numerical derivative epsilon
    const double h = 1e-4;

    for (std::size_t l = 0; l < par.L.size(); ++l) {
      for (std::size_t i = 0; i < (unsigned int)par.L[l].w.rows(); ++i) {
        for (std::size_t j = 0; j < (unsigned int)par.L[l].w.cols(); ++j) {
          par.L[l].w(i, j)   = w[k] + h;
          const double f_pos = val(loss(Z, par));
          par.L[l].w(i, j)   = w[k] - h;
          const double f_neg = val(loss(Z, par));

          gradv[k] = (f_pos - f_neg) / (2 * h);

          par.L[l].w(i, j) = w[k];  // back to previous value
          ++k;
        }
      }
    }

    for (std::size_t l = 0; l < par.L.size(); ++l) {
      for (std::size_t j = 0; j < (unsigned int)par.L[l].b.size(); ++j) {
        par.L[l].b(j)      = w[k] + h;
        const double f_pos = val(loss(Z, par));
        par.L[l].b(j)      = w[k] - h;
        const double f_neg = val(loss(Z, par));

        gradv[k]      = (f_pos - f_neg) / (2 * h);
        par.L[l].b(j) = w[k];  // back to previous value
        ++k;
      }
    }

    const double cost = val(loss(Z, par));
    std::cout << "loss = " << cost << std::endl;
    return cost;
  }

  // Naive gradient descent
  void NaiveGradDesc(std::vector<double> &w, unsigned int MAXITER, double rate) {
    // -----------------------------------------------------------------------
    // Gradient Descent
    // w_n+1 = w_n - gamma*Grad(F(w_n))
    unsigned int iter = 0;

    // gradvient vector
    std::vector<double> gradv(w.size(), 0.0);

    while (true) {
      CostGrad(w, gradv);
      // Naive gradient descent update
      for (std::size_t k = 0; k < gradv.size(); ++k) {
        w[k] = w[k] - rate * gradv[k];
        printf("w[%3lu] = %6.3f    <gradv[%3lu] = %10.3E> \n", k, w[k], k, gradv[k]);
      }
      ++iter;

      if (iter > MAXITER) break;
    }
  }

  // Test network dimensions
  void TestNetworkDimensions() {
    for (std::size_t l = 0; l < par.L.size() - 1; ++l) {
      if (par.L[l + 1].w.cols() != par.L[l].w.rows()) {
        throw std::invalid_argument("Network layer dimension mismatch!");
      }
    }
  }

  // Optimize network
  void Optimize() {
    // ------------------------------------------------------------------
    // Random init (too high sigma -> may get trapped to a bad local minima)
    double              sigma = 1e-4;
    std::vector<double> w     = RandomInit(par, sigma);

    // ------------------------------------------------------------------
    // By Jacobi formula: d/dt ln det A(t) = tr[A(t)^{-1} d/dt A]
    // ------------------------------------------------------------------

    const unsigned int n = w.size();

    LBFGSParam<double> param;

    // Initial guess
    VectorXd x = VectorXd::Zero(n);
    for (std::size_t i = 0; i < n; ++i) { x[i] = w[i]; }

    // Function object
    Neurocost fun(n);

    // x will be overwritten to be the best point found
    double fx;
    int    niter;

    FLAT_TARGET = true;
    std::cout << "Solving flat target: " << std::endl;

    param.epsilon        = 1e-4;
    param.max_iterations = 10;

    LBFGSSolver<double> solver0(param);
    niter = solver0.minimize(fun, x, fx);

    FLAT_TARGET = false;
    std::cout << "Solving real target: " << std::endl;

    param.epsilon        = 1e-4;
    param.max_iterations = 20;

    LBFGSSolver<double> solver1(param);
    niter = solver1.minimize(fun, x, fx);

    std::cout << niter << " iterations" << std::endl;
    std::cout << "x = \n" << x.transpose() << std::endl;
    std::cout << "f(x) = " << fx << std::endl;
  }

 private:
  class Neurocost {

  public:
    Neurocost(int _n) {
      n = _n;
    }
    
    double fx = 0.0;
    double operator()(const VectorXd &x, VectorXd &grad) {
      // Map to vectors
      std::vector<double> w(x.size(), 0.0);
      std::vector<double> gradv(x.size(), 0.0);

      for (std::size_t i = 0; i < (unsigned int)x.size(); ++i) { w[i] = x[i]; }

      fx = CostGrad(w, gradv);

      for (std::size_t i = 0; i < (unsigned int)x.size(); ++i) { grad[i] = gradv[i]; }

      return fx;
    }

    private:
      int n;
  };

};  // MNeuroJacobian

}  // Namespace neurojac
}  // Namespace gra

#endif