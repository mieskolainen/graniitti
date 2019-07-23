// GRANIITTI - Monte Carlo event generator for high energy diffraction
// https://github.com/mieskolainen/graniitti
//
// <Sommerfeld-Rayleigh diffraction benchmark>
//
// We work in the point source rest frame
// 
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <cassert>
#include <complex>
#include <future>
#include <iostream>
#include <random>
#include <stdexcept>
#include <valarray>
#include <vector>
#include <fstream>
#include <string>

// Own
#include "Graniitti/MAux.h"
#include "Graniitti/MH1.h"
#include "Graniitti/MMath.h"
#include "Graniitti/MTimer.h"

// Libraries
#include "cxxopts.hpp"
#include "json.hpp"
#include "rang.hpp"

using gra::aux::indices;
using namespace gra;

const double c = 1.0; // Wavespeed
double    kmod = 0.0;

// Source position (x,y,z,t)
M4Vec xs(0, 0, -1000, 0);

// Unit normal inwards (x,y,z,t)
const M4Vec n(0, 0, 1.0, 0);

const double EPSILON = 1e-9;


// Spherical point source field at 4-position x
//
// Source 4-position: xs
//
inline std::complex<double> u_point(const M4Vec& x) {

  double s = (x - xs).P3mod();
  if (s < EPSILON) { return 1.0; }

  const double omega = kmod * c;
  const double     t = x.T();

  return std::exp(math::zi * (omega * t - kmod * s)) / s;
}


// Sommerfeld-Rayleigh aperture boundary integrand
inline std::complex<double> RS_integrand(const M4Vec& x, const M4Vec& x0) {

  // Detector -> Aperture vector
  const double x_x0 = (x - x0).P3mod();

  if (x_x0 < EPSILON) { return u_point(x); }

  // SR-integral
  // Incoming field u(x)
  // Green function G(x-x0)
  // Inclination factor: cos delta(x,x0)
  const std::complex<double> out =
        -math::zi * 
         kmod / (2.0*math::PI) *
         u_point(x) *
         std::exp(-math::zi * kmod * x_x0) / x_x0 *
         n.Dot3(x - x0) / x_x0;

  return out;
}

struct LIM {
  std::vector<double> xlim = {0.0, 0.0};
  std::vector<double> ylim = {0.0, 0.0};
};


// Main
int main(int argc, char *argv[]) {

  MTimer global_timer;
  MTimer timer;

  std::vector<std::future<std::complex<double>>> futures;  // std::async return values

  // Simpson weight matrix (too small N will give visible artifacts at the boundary)
  const int M = 33;
  const MMatrix<double> SW = math::Simpson38Weight2D(M,M);

  // ------------------------------------------------------------
  // Open aperture (A) definition

  std::vector<LIM> A;
  {
    LIM a;
    a.xlim = {-0.2,  0.2};
    a.ylim = {-0.2,  0.2};
    A.push_back(a);
  }
  {
    LIM a;
    a.xlim = {-0.2,  0.2};
    a.ylim = { 0.6,  0.8};
    A.push_back(a);
  }
  {
    LIM a;
    a.xlim = {-0.2,  0.2};
    a.ylim = {-0.8, -0.6};
    A.push_back(a);
  }
  
  // ------------------------------------------------------------

  MMatrix<std::complex<double>> f(M+1, M+1);

  const double kmag = 20 / std::sqrt(3.0);
  const M4Vec k4(kmag, kmag, kmag, 0.0);
  kmod = k4.P3mod();

  // 4D-Grid discretization
  const std::vector<double> xval = math::linspace( 0.0, 0.0, 1);
  const std::vector<double> yval = math::linspace(-1.0, 1.0, 50);
  const std::vector<double> zval = math::linspace(-1.0, 5.0, 150);
  const std::vector<double> tval = math::linspace( 0.0, 2.0*math::PI/kmod, 20); // One period
  
  // Aperture vector
  M4Vec  x(0, 0, 0, 0);
  M4Vec x0(0, 0, 0, 0);

  MTensor<std::complex<double>> tensor = MTensor({xval.size(), yval.size(), zval.size(), tval.size()}, std::complex<double>(0.0));

  // 3D-grid
  for (std::size_t i = 0; i < xval.size(); ++i) {
  for (std::size_t j = 0; j < yval.size(); ++j) {
  for (std::size_t k = 0; k < zval.size(); ++k) {
  for (std::size_t l = 0; l < tval.size(); ++l) {

    // Detector 4-position
    x0.Set(xval[i], yval[j], zval[k], tval[l]);

    // Source side
    if (zval[k] <= 0.0) {

      // Spherical wave only
      tensor({i,j,k,l}) = u_point(x0);

    // Scattering side
    } else {

      std::complex<double> I = 0.0;

      // Loop over apertures
      for (std::size_t a = 0; a < A.size(); ++a) {

        const double hstepX = (A[a].xlim[1] - A[a].xlim[0]) / M;
        const double hstepY = (A[a].ylim[1] - A[a].ylim[0]) / M;

        for (std::size_t u = 0; u < f.size_row(); ++u) {
          for (std::size_t v = 0; v < f.size_col(); ++v) {

            // Aperture vector
            x.Set(A[a].xlim[0] + u * hstepX, A[a].ylim[0] + v * hstepY, 0, tval[l]);
            f[u][v] = RS_integrand(x, x0);
          }
        }
        // Integrate
        I += math::Simpson38Integral2D(f, SW, hstepX, hstepY);
      }
      // Save real part
      tensor({i,j,k,l}) = I;
    }

    if (timer.ElapsedSec() > 0.1) {
      timer.Reset();
      gra::aux::PrintProgress(j / static_cast<double>(yval.size()));
    }

  }}}}

  gra::aux::ClearProgress();  // Clear progressbar

  const std::vector<std::vector<size_t>> allind = {math::linspace(size_t(0), xval.size()-1, xval.size()), 
                                                   math::linspace(size_t(0), yval.size()-1, yval.size()),
                                                   math::linspace(size_t(0), zval.size()-1, zval.size()),
                                                   math::linspace(size_t(0), tval.size()-1, tval.size())};

  std::vector<size_t> current;
  std::vector<std::vector<std::size_t>> output;
  math::IndexComb(allind, 0, current, output);

  std::ofstream textout("3D.ascii");

  for (std::size_t i = 0; i < output.size(); ++i) {
    for (std::size_t j = 0; j < output[i].size(); ++j) {
      textout << output[i][j] << ",";
    }
    textout << std::real(tensor(output[i])) << "," << std::imag(tensor(output[i]));
    textout << std::endl;
  }

  printf("Elapsed time %0.1f sec \n", global_timer.ElapsedSec());

  return EXIT_SUCCESS;
}
