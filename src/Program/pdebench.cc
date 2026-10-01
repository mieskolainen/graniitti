// Evolve the periodic Kuramoto Sivashinsky equation
//
// (c) 2026 Mikael Mieskolainen
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

// Own
#include "Graniitti/Math/MAlgorithms.h"
#include "Graniitti/Math/MFFT.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Program/MInput.h"
#include "Graniitti/Program/pdebench.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MH1.h"
#include "Graniitti/Tech/MTimer.h"

// Libraries
#include "rang.hpp"

// Eigen
#include <Eigen/unsupported/Eigen/FFT>

using gra::aux::indices;
using namespace gra;

// Evolve and display the coordinate space field at the requested intervals
void CNAB2(std::valarray<std::complex<double>>& u, int Lx, double dt, int Nt, int nplot, bool use_eigen) {
  gra::program::ValidateKS(Lx, static_cast<int>(u.size()), dt, Nt, nplot, use_eigen);
  gra::program::KS evolution(u, Lx, dt, use_eigen);
  for (int step = 0; step < Nt; ++step) {
    evolution.Step();
    if (step % nplot != 0) { continue; }
    const auto field = evolution.Field();
    double     peak  = 0.0;
    for (const auto& i : indices(field)) { peak = std::max(peak, std::norm(field[i])); }
    const std::string shades = " .:-=+*#%@";
    std::cout << step << " ";
    for (const auto& i : indices(field)) {
      const double value = std::norm(field[i]) / (peak + 1e-9);
      const auto   index = static_cast<std::size_t>(std::clamp(std::round(value * 9), 0.0, 9.0));
      std::cout << shades[index];
    }
    std::cout << std::endl;
  }
  u = evolution.Field();
}

// Display a forward and inverse FFT example
void testfft() {
  // FFT test
  std::valarray<float> xxx = gra::math::linspace<std::valarray>((float)0.0, (float)10.0, (unsigned int)16);

  std::valarray<std::complex<double>> X(xxx.size());
  for (const auto& i : indices(xxx)) { X[i] = xxx[i]; }
  MFFT::fft(X);
  std::cout << "FFT:" << std::endl;
  for (const auto& i : indices(X)) { std::cout << X[i] << std::endl; }
  MFFT::ifft(X);
  std::cout << "IFFT:" << std::endl;
  for (const auto& i : indices(X)) { std::cout << X[i] << std::endl; }
}

// Main
int main(int argc, char* argv[]) {
  aux::PrintArgv(argc, argv);

  if (argc != 7) {
    std::cout << "Kuramoto-Sivashinsky 4th order PDE FFT-testbench" << std::endl;
    std::cout << "Usage: ./pdebench 64 128 16 4 2400 0" << std::endl;
    return EXIT_FAILURE;
  }

  try {
    int       Lx       = gra::program::Number<int>(argv[1]);
    int       Nx       = gra::program::Number<int>(argv[2]);
    double    dt       = 1.0 / static_cast<double>(gra::program::Number<int>(argv[3]));
    int       nplot    = gra::program::Number<int>(argv[4]);
    int       Nt       = gra::program::Number<int>(argv[5]);
    const int fft_type = gra::program::Number<int>(argv[6]);
    if (fft_type != 0 && fft_type != 1) { throw std::invalid_argument("pdebench: FFT type must be 0 or 1"); }
    const bool USE_EIGEN = fft_type == 1;
    gra::program::ValidateKS(Lx, Nx, dt, Nt, nplot, USE_EIGEN);
    testfft();

    std::valarray<std::complex<double>> x(Nx);
    for (const auto& i : indices(x)) { x[i] = static_cast<double>(i) * Lx / static_cast<double>(Nx); }
    std::valarray<std::complex<double>> u(Nx);
    for (const auto& i : indices(u)) {
      // Initialize with single wave
      // u[i] = std::cos(x[i]/(double)Lx)*(1.0 + std::sin(x[i]/(double)Lx));

      // Initialize with a canonic one
      u[i] = std::cos(x[i]) + 0.1 * std::cos(x[i] / 16.0) * (1.0 + 2.0 * std::sin(x[i] / 16.0));
    }
    gra::PrintArray(x, "x");
    gra::PrintArray(u, "u");

    // Solve the equation
    CNAB2(u, Lx, dt, Nt, nplot, USE_EIGEN);

    std::cout << "[pdebench: done]" << std::endl;

    return EXIT_SUCCESS;
  } catch (const std::exception& error) {
    std::cerr << "pdebench: " << error.what() << std::endl;
    return EXIT_FAILURE;
  }
}
