// Sample free particle lattice paths through a slit
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <cassert>
#include <complex>
#include <iostream>
#include <random>
#include <stdexcept>
#include <valarray>
#include <vector>

// Own
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Program/MInput.h"
#include "Graniitti/Program/pathmark.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MH1.h"
#include "Graniitti/Tech/MTimer.h"

// Libraries
#include "cxxopts.hpp"
#include "json.hpp"
#include "rang.hpp"

using gra::aux::indices;
using namespace gra;

// Main
int main(int argc, char* argv[]) {
  aux::PrintArgv(argc, argv);

  try {
    gra::program::PathParam param;
    auto& [N, dt, m, k_slit, SLIT, NS] = param;

    // Detector observable histogram
    const double              XWIDTH = 30;
    MH1<std::complex<double>> hx(100, -XWIDTH, XWIDTH, "Path Integral MC");

    // Seed here
    std::random_device rd;
    // Standard mersenne_twister_engine seeded with rd()
    std::mt19937_64 gen(rd());

    // double offset = 2.0;

    if (argc != 1 && argc != 7) { throw std::invalid_argument("Usage: pathmark N dt m k_slit slit samples"); }
    if (argc == 1) {  // Default

      N      = 4;
      dt     = 1.0;
      m      = 1.5;
      k_slit = 2;
      SLIT   = 0.02;
      NS     = 10000000;

      printf("Example input ./pathmark %d %0.3f %0.3f %d %0.3f %lu \n", N, dt, m, k_slit, SLIT, NS);
      std::cout << std::endl;
    } else {
      // Collect input variables
      N      = gra::program::Number<int>(argv[1]);
      dt     = gra::program::Number<double>(argv[2]);
      m      = gra::program::Number<double>(argv[3]);
      k_slit = gra::program::Number<int>(argv[4]);  // Slit position
      SLIT   = gra::program::Number<double>(argv[5]);
      NS     = gra::program::Number<unsigned long>(argv[6]);
    }

    param.Validate(XWIDTH);

    // Random variables
    std::uniform_real_distribution<> dis(-XWIDTH, XWIDTH);
    std::uniform_real_distribution<> slit(-SLIT, SLIT);

    // Imaginary unit
    const std::complex<double> ii(0, 1);

    // Path configuration vector
    std::vector<double> x(N, 0.0);

    MTimer globaltimer(true);
    MTimer localtimer(true);

    // MC samples
    for (unsigned long i = 0; i < NS; ++i) {
      // Draw a single path configuration
      for (int k = 1; k < N; ++k) { x[k] = dis(gen); }

      // Apply double slit boundary condition to two time steps
      // (one is not enough)
      x[k_slit]     = slit(gen);
      x[k_slit - 1] = slit(gen);

      // ** Get complex action weight **
      std::complex<double> W = std::exp(ii * gra::program::S_lat(x, param));

      // Fill histogram
      hx.Fill(x[N - 1], W);

      if (localtimer.ElapsedSec() > 0.1) {
        localtimer.Reset();
        gra::aux::PrintProgress(i / static_cast<double>(NS));
      }
    }
    gra::aux::ClearProgress();  // Clear progressbar

    // Print histogram
    hx.Print();

    double runtime = globaltimer.ElapsedSec();
    std::cout << rang::style::bold << "<pathmark - complex path integral CPU benchmark>" << rang::style::reset
              << std::endl
              << std::endl;
    printf("- Runtime         = %0.2E sec\n", runtime);
    printf("- Action integral = %0.3E eval/sec\n", NS / runtime);
    printf("\n");


    std::cout << std::endl;

    hx.RawOutput();

    std::cout << "[pathmark: done]" << std::endl;

    return EXIT_SUCCESS;
  } catch (const std::exception& error) {
    std::cerr << "pathmark: " << error.what() << std::endl;
    return EXIT_FAILURE;
  }
}
