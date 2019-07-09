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

// Source position (x,y,z,t)
M4Vec xs(0, 0, -5, 0);

const double EPSILON = 1e-10;


// Spherical point source field at 4-position x
//
// Source 4-position: xs
//
inline std::complex<double> u_point(const M4Vec& x, const M4Vec& k4) {

  const M4Vec d = x - xs;

  double s = d.P3mod();
  if (s < EPSILON) { s = EPSILON; }

  const double  kmod = k4.P3mod();
  const double omega = kmod * c;
  const double     t = x.T();

  return std::exp(math::zi * (omega * t - kmod * s)) / s;
}


// Sommerfeld-Rayleigh aperture boundary integral
inline std::complex<double> RS_integral(const M4Vec& x, const M4Vec& x0, const M4Vec& k4) {

  // Unit normal inwards (x,y,z,t)
  const M4Vec n(0, 0, 1.0, 0);

  // Detector -> Aperture vector
  const M4Vec x_x0vec = x - x0;
  double x_x0 = x_x0vec.P3mod();

  if (x_x0 < EPSILON) { x_x0 = EPSILON; }

  const double kmod = k4.P3mod();

  // SR-integral
  // Incoming field u(x)
  // Green function G(x-x0)
  // Inclination factor: cos delta(x,x0)
  const std::complex<double> out =
        -math::zi * kmod / (2.0*math::PI) *
         u_point(x, k4) *
         std::exp(-math::zi * kmod * x_x0) / x_x0 *
         n.Dot3(x_x0vec) / x_x0;

  return out;
}


// Main
int main(int argc, char *argv[]) {

  MTimer timer(true);

  std::vector<std::future<std::complex<double>>> futures;  // std::async return values

  // Simpson weight matrix
  const int M = 20;
  const int N = 20;
  const MMatrix<double> SW = math::Simpson13Weight2D(M,N);

  // Aperture limits
  const double a = -0.3;
  const double b =  0.3;
  const double c = -0.3;
  const double d =  0.3;

  const double hstepM = (b-a) / M;
  const double hstepN = (d-c) / N;
  MMatrix<std::complex<double>> f(M+1, N+1);

  // 4D-Grid discretization
  const std::vector<double> xval = math::linspace(-1.0, 1.0, 150);
  const std::vector<double> yval = math::linspace(-1.0, 1.0, 150);
  const std::vector<double> zval = math::linspace(-1.0, 1.0, 150);
  const std::vector<double> tval = math::linspace( 0.0, 1.0, 20);

  const double kmag = 40 / std::sqrt(3.0);
  const M4Vec k4(kmag, kmag, kmag, 0.0);

  // Aperture vector
  M4Vec  x(0, 0, 0, 0);
  M4Vec x0(0, 0, 0, 0);

  MTensor<double> tensor = MTensor({xval.size(), yval.size(), zval.size(), tval.size()}, 0.0);

  // 3D-grid
  for (std::size_t i = 0; i < xval.size(); ++i) {
  for (std::size_t j = 0; j < yval.size(); ++j) {
  for (std::size_t k = 0; k < zval.size(); ++k) {
  for (std::size_t l = 0; l < tval.size(); ++l) {

    // Detector 4-position
    x0.Set(xval[i], yval[j], zval[k], tval[l]);

    if (zval[k] < 0.0) {
      // Spherical wave only
      tensor({i,j,k,l}) = std::real( u_point(x0, k4) );
    } else {

      // Sample function values for the Sommerfeld-Rayleigh integral
      for (int u = 0; u < M+1; ++u) {
        for (int v = 0; v < N+1; ++v) {

          // Aperture vector
          x.Set(a + u * hstepM, c + v * hstepN, 0, tval[l]);

          f[u][v] = RS_integral(x, x0, k4);
        }
      }
      // Calculate aperture surface 2D-integral    
      tensor({i,j,k,l}) = std::real( math::Simpson13Integral2D(f, SW, hstepM, hstepN) );
    }

    if (timer.ElapsedSec() > 0.1) {
      timer.Reset();
      gra::aux::PrintProgress(i / static_cast<double>(xval.size()));
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
    textout << tensor(output[i]);
    textout << std::endl;
  }

  return EXIT_SUCCESS;
}
