// Compute diffraction from a spherical source
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <cassert>
#include <complex>
#include <fstream>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <valarray>
#include <vector>

// Own
#include "Graniitti/Math/MAlgorithms.h"
#include "Graniitti/Math/MCombinatorics.h"
#include "Graniitti/Math/MIntegration.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Program/sommerfeld.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MH1.h"
#include "Graniitti/Tech/MTimer.h"

// Libraries
#include "cxxopts.hpp"
#include "json.hpp"
#include "rang.hpp"

using gra::aux::indices;
using namespace gra;

struct LIM {
  std::vector<double> xlim = {0.0, 0.0};
  std::vector<double> ylim = {0.0, 0.0};
};


// Main
int main(int argc, char* argv[]) {
  aux::PrintArgv(argc, argv);

  try {
    std::ofstream textout;
    textout.exceptions(std::ios::failbit | std::ios::badbit);
    textout.open("3D.ascii");

    MTimer global_timer;
    MTimer timer;

    // Simpson weight matrix for the aperture quadrature
    const int             M  = 33;
    const MMatrix<double> SW = math::Simpson38Weight2D(M, M);

    // ------------------------------------------------------------
    // Open aperture (A) definition

    std::vector<LIM> A;
    {
      LIM a;
      a.xlim = {-0.2, 0.2};
      a.ylim = {-0.2, 0.2};
      A.push_back(a);
    }
    {
      LIM a;
      a.xlim = {-0.2, 0.2};
      a.ylim = {0.6, 0.8};
      A.push_back(a);
    }
    {
      LIM a;
      a.xlim = {-0.2, 0.2};
      a.ylim = {-0.8, -0.6};
      A.push_back(a);
    }

    // ------------------------------------------------------------

    MMatrix<std::complex<double>> f(M + 1, M + 1);

    const double kmag = 20 / std::sqrt(3.0);
    const M4Vec  k4(kmag, kmag, kmag, 0.0);
    const double kmod = k4.P3mod();
    const M4Vec  source(0, 0, -1000, 0);
    const M4Vec  normal(0, 0, 1, 0);

    // 4D-Grid discretization
    const std::vector<double> xval = math::linspace(0.0, 0.0, 1);
    const std::vector<double> yval = math::linspace(-1.0, 1.0, 100);
    const std::vector<double> zval = math::linspace(-1.0, 5.0, 300);
    const std::vector<double> tval = math::linspace(0.0, 2.0 * math::PI / kmod, 20);  // One period

    // Aperture vector
    M4Vec x(0, 0, 0, 0);
    M4Vec x0(0, 0, 0, 0);

    MTensor<std::complex<double>> tensor =
        MTensor({xval.size(), yval.size(), zval.size(), tval.size()}, std::complex<double>(0.0));

    // 3D-grid
    for (const auto& i : indices(xval)) {
      for (const auto& j : indices(yval)) {
        for (const auto& k : indices(zval)) {
          for (const auto& l : indices(tval)) {
            // Detector 4-position
            x0.Set(xval[i], yval[j], zval[k], tval[l]);

            // Source side
            if (zval[k] <= 0.0) {
              // Spherical wave only
              tensor({i, j, k, l}) = gra::program::u_point(x0, source, kmod);

              // Scattering side
            } else {
              std::complex<double> I = 0.0;

              // Loop over apertures
              for (const auto& a : indices(A)) {
                const double hstepX = (A[a].xlim[1] - A[a].xlim[0]) / M;
                const double hstepY = (A[a].ylim[1] - A[a].ylim[0]) / M;

                for (std::size_t u = 0; u < f.size_row(); ++u) {
                  for (std::size_t v = 0; v < f.size_col(); ++v) {
                    // Aperture vector
                    x.Set(A[a].xlim[0] + u * hstepX, A[a].ylim[0] + v * hstepY, 0, tval[l]);
                    f[u][v] = gra::program::RS_integrand(x, x0, source, normal, kmod);
                  }
                }
                // Integrate
                I += math::Simpson38Integral2D(f, SW, hstepX, hstepY);
              }
              // Save real part
              tensor({i, j, k, l}) = I;
            }

            if (timer.ElapsedSec() > 0.1) {
              timer.Reset();
              gra::aux::PrintProgress(j / static_cast<double>(yval.size()));
            }
          }
        }
      }
    }

    gra::aux::ClearProgress();  // Clear progressbar

    const std::vector<std::vector<size_t>> allind = {
        math::linspace(size_t(0), xval.size() - 1, xval.size()), math::linspace(size_t(0), yval.size() - 1, yval.size()),
        math::linspace(size_t(0), zval.size() - 1, zval.size()), math::linspace(size_t(0), tval.size() - 1, tval.size())};

    std::vector<size_t>                   current;
    std::vector<std::vector<std::size_t>> output;
    math::IndexComb(allind, 0, current, output);

    for (const auto& i : indices(output)) {
      for (const auto& j : indices(output[i])) { textout << output[i][j] << ","; }
      textout << std::real(tensor(output[i])) << "," << std::imag(tensor(output[i]));
      textout << std::endl;
    }

    textout.close();

    printf("Elapsed time %0.1f sec \n", global_timer.ElapsedSec());
  } catch (const std::exception& error) {
    std::cerr << "sommerfeld: " << error.what() << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
