// Optimal transport demonstration program
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <fstream>
#include <iomanip>

// HepMC3
#include "HepMC3/FourVector.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/LHEFAttributes.h"
#include "HepMC3/Print.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/Relatives.h"
#include "HepMC3/Selector.h"
#include "HepMC3/WriterAscii.h"

// Own
#include "Graniitti/Math/MAlgorithms.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Math/MTransport.h"
#include "Graniitti/Particle/MPDG.h"
#include "Graniitti/Program/MInput.h"
#include "Graniitti/Program/ot.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MH1.h"


// Libraries
#include "cxxopts.hpp"

using gra::aux::indices;

using namespace gra;


// Main function
int main(int argc, char* argv[]) {
  aux::PrintArgv(argc, argv);

  std::cout << std::endl;
  std::cout << "OT" << std::endl;

  // Parameters
  std::vector<std::string> input;
  double                   lambda = 1.0;
  unsigned int             iter   = 10;

  // Save the number of input arguments
  const int NARGC = argc - 1;
  try {
    cxxopts::Options options(argv[0], "");
    options.add_options("")("i, input", "Input files            <input1.hepmc3,input2.hepmc3>",
                            cxxopts::value<std::string>())("a, lambda", "Regularization lambda  <double>",
                                                           cxxopts::value<std::string>()->default_value("1.0"))(
        "r, iter", "Number of iterations   <integer>", cxxopts::value<std::string>()->default_value("10"))("H, help", "Help");

    auto r = options.parse(argc, argv);

    if (r.count("help") || NARGC == 0) {
      std::cout << options.help({""}) << std::endl;
      std::cout << rang::style::bold << "Example:" << rang::style::reset << std::endl;
      std::cout << "  " << argv[0] << " -i input1.hepmc3,input2.hepmc3 -a 1.5 -r 15" << std::endl << std::endl;
      return EXIT_FAILURE;
    }

    // Read input
    input  = gra::aux::SplitStr2Str(r["i"].as<std::string>(), ',');
    lambda = gra::program::Number<double>(r["a"].as<std::string>());
    iter   = gra::program::Number<unsigned int>(r["r"].as<std::string>());
    if (!(lambda > 0.0) || iter == 0) { throw std::invalid_argument("lambda and iterations must be positive"); }

  } catch (const std::exception& error) {
    std::cerr << "ot: Error reading command-line input: " << error.what() << std::endl;
    return EXIT_FAILURE;
  } catch (...) {
    std::cerr << "ot: Non-standard exception while reading command-line input" << std::endl;
    return EXIT_FAILURE;
  }

  if (input.size() != 2) {
    std::cerr << "ot: exactly two input files are required" << std::endl;
    return EXIT_FAILURE;
  }

  std::string inputfile1(input[0]);
  std::string inputfile2(input[1]);

  try {
    // ------------------------------------------------------------
    // Create histograms
    const int BINS1 = 60;
    const int BINS2 = 60;


    int MODE = 1;

    // Transport matrix to be calculated by optimization
    gra::MMatrix<double> Pi;


    if (MODE == 1) {
      MH1<double> h1(BINS1, 0.28, 3.0, "data_1");
      MH1<double> h2(BINS2, 0.28, 3.0, "data_2");

      gra::program::ReadMass(inputfile1, h1);
      gra::program::ReadMass(inputfile2, h2);

      std::vector<double> p = h1.GetProbDensity();
      std::vector<double> q = h2.GetProbDensity();

      // ------------------------------------------------------------
      /*
      // l2-metric squared || x - y ||^2
      auto l2metric = [] (const std::vector<double>& x, const std::vector<double>& y) {

        double value = 0.0;
        for (std::size_t i = 0; i < x.size(); ++i) {
          value += gra::math::pow2(x[i] - y[i]);
        }
        return value;
      };

      // Distance matrix
      gra::MMatrix<double> C(S1.size(), S2.size());

      for (std::size_t i = 0; i < C.size_row(); ++i) {
        for (std::size_t j = 0; j < C.size_col(); ++j) {
          //C[i][j] = l2metric(S1[i], S2[j]);

          C[i][j] = math::pow2(p[i] - q[j]);
        }
      }
      */

      // Kernel matrix for comparing histograms
      MMatrix<double> K;
      gra::opt::ConvKernel(BINS1, BINS2, lambda, K);

      /*
      gra::opt::GibbsKernel(lambda, C, K);

      // Target histograms (for unweighted set == 1/size)
      std::vector<double> p(S1.size(), 1.0 / S1.size());
      std::vector<double> q(S2.size(), 1.0 / S2.size());
      */

      gra::opt::SinkHorn(Pi, K, p, q, iter);
      Pi.Print();
    }

    // ---------------------------------------------------------------
    // Gaussian histograms
    if (MODE == 2) {
      MMatrix<double> K;
      gra::opt::ConvKernel(BINS1, BINS2, lambda, K);

      // Normal function
      auto gaussfunc = [](double mu, double sigma, std::size_t N) {
        std::vector<double> y(N);

        // Create x-value
        std::vector<double> x = math::linspace(0.0, N - 1.0, N);
        gra::Scale(x, 1.0 / (N - 1.0));

        for (std::size_t i = 0; i < N; ++i) { y[i] = std::exp(-math::pow2(x[i] - mu) / (2 * math::pow2(sigma))); }
        return y;
      };

      auto normfunc = [](std::vector<double>& x) {
        const double alpha  = 0.02;
        const double maxval = *std::max_element(x.begin(), x.end());
        for (const auto& i : indices(x)) { x[i] += maxval * alpha; }
        x = gra::NormalizedSum(x);
      };

      /*
      auto printfunc = [] (const std::vector<double>& x) {
        for (std::size_t i = 0; i < x.size(); ++i) {
          std::cout << x[i] << std::endl;
        }
      };
      */

      const double        sigma = 0.06;
      std::vector<double> p     = gaussfunc(0.25, sigma, BINS1);
      std::vector<double> q     = gaussfunc(0.8, sigma, BINS2);
      normfunc(p);
      normfunc(q);

      // ---------------------------------------------------------------

      gra::opt::SinkHorn(Pi, K, p, q, iter);
      Pi.Print();
    }

    return EXIT_SUCCESS;
  } catch (const std::exception& error) {
    std::cerr << "ot: " << error.what() << std::endl;
    return EXIT_FAILURE;
  }
}
