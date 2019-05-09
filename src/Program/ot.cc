// GRANIITTI - Monte Carlo event generator for high energy diffraction
// https://github.com/mieskolainen/graniitti
//
// <OT>
//
// (c) 2017-2019 Mikael Mieskolainen
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
#include "Graniitti/MPDG.h"
#include "Graniitti/MTransport.h"
#include "Graniitti/MMatrix.h"
#include "Graniitti/MMath.h"

// Libraries
#include "cxxopts.hpp"


void ReadEvents(const std::string& inputfile, std::vector<std::vector<double>>& S) {

  // Input and output
  HepMC3::ReaderAscii input(inputfile);

  int events = 0;

  // HepMC3 event object
  HepMC3::GenEvent ev(HepMC3::Units::GEV, HepMC3::Units::MM);

  // Event loop over all events in HepMC3 file
  while (!input.failed()) {

    input.read_event(ev);
    if (input.failed()) break;

    // Event vector
    std::vector<double> x;

    // Loop over all particles
    for (HepMC3::ConstGenParticlePtr p1 : ev.particles()) {

      // Check is final state and pion
      if (p1->status() == gra::PDG::PDG_STABLE && std::abs(p1->id()) < 1000) {

        // Take 4-momentum
        HepMC3::FourVector pvec = p1->momentum();

        //x.push_back(pvec.e());
        x.push_back(pvec.px());
        x.push_back(pvec.py());
        x.push_back(pvec.pz());
      }
    }

    S.push_back(x);
    ++events;

    if (events % 10000 == 0) {
      std::cout << events << " events processed" << std::endl;
    }
  }
  
  printf("From input: %s found %lu events \n", inputfile.c_str(), S.size());

}

// Main function
int main(int argc, char *argv[]) {

  std::cout << std::endl;
  std::cout << "OT" << std::endl;

  // Parameters
  std::vector<std::string> input;
  double lambda = 1.0;
  unsigned int iter = 10;

  // Save the number of input arguments
  const int NARGC = argc - 1;
  try {
    cxxopts::Options options(argv[0], "");
    options.add_options("")(
        "i, input",   "Input files            <input1.hepmc3,input2.hepmc3,...>", cxxopts::value<std::string>())(
        "a, lambda",  "Regularization lambda  <double>",  cxxopts::value<double>())(
        "r, iter",    "Number of iterations   <integer>", cxxopts::value<unsigned int>())(
        "H, help",    "Help");

    auto r = options.parse(argc, argv);

    if (r.count("help") || NARGC == 0) {
      std::cout << options.help({""}) << std::endl;
      std::cout << rang::style::bold << "Example:" << rang::style::reset << std::endl;
      std::cout << "  " << argv[0] << " -i input1.hepmc3,input2.hepmc3 -a 1.5 -r 15" << std::endl << std::endl;
      return EXIT_FAILURE;
    }

    // Read input
    input = gra::aux::SplitStr2Str(r["i"].as<std::string>(), ',');
    lambda = r["a"].as<double>();
    iter   = r["r"].as<unsigned int>();

  } catch (...) {
    std::cout << "unknown error" << std::endl;
    return EXIT_FAILURE;
  }

  if (input.size() < 2) {
    std::cout << "input should be >= 2" << std::endl;
    return EXIT_FAILURE;
  }

  std::string inputfile1(input[0]);
  std::string inputfile2(input[1]);

  // ------------------------------------------------------------

  std::vector<std::vector<double>> S1;
  ReadEvents(inputfile1, S1);

  std::vector<std::vector<double>> S2;
  ReadEvents(inputfile2, S2);

  // ------------------------------------------------------------

  // l2-metric || x - y ||^2
  auto l2metric = [] (const std::vector<double>& x, const std::vector<double>& y) {

    double value = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
      value += gra::math::pow2(x[i] - y[i]);
    }
    return value;
  };

  // Minkowski
  auto minkowski = [] (const std::vector<double>& x, const std::vector<double>& y) {

    double value = 0.0;
    for (std::size_t i = 0; i < x.size(); i += 4) { // jump in 4
      // e^2 - pvec^2 (+,---)
      value += gra::math::pow2(x[i] - y[i]) - (gra::math::pow2(x[i+1] - y[i+1]) + gra::math::pow2(x[i+2] - y[i+2]) + gra::math::pow2(x[i+3] - y[i+3])); 
    }
    return value;
  };

  // Distance matrix
  gra::MMatrix<double> C(S1.size(), S2.size());

  for (std::size_t i = 0; i < C.size_row(); ++i) {
    for (std::size_t j = 0; j < C.size_col(); ++j) {
      //C[i][j] = minkowski(S1[i], S2[j]);
      C[i][j] = l2metric(S1[i], S2[j]); 
    }
  }

  // Transport matrix
  gra::MMatrix<double> gamma;

  //C.Print();

  // Event vector weights (for unweighted set == 1/size)
  std::vector<double> p1(S1.size(), 1.0 / S1.size());
  std::vector<double> p2(S2.size(), 1.0 / S2.size());

  double cost = gra::opt::SinkHorn(gamma, C, p1, p2, lambda, iter);

  printf("Sinkhorn cost = %0.5E \n", cost);

  return EXIT_SUCCESS;
}

