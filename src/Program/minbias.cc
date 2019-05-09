// GRANIITTI - Monte Carlo event generator for high energy diffraction
// https://github.com/mieskolainen/graniitti
//
// <minimum bias processes combined>
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <math.h>
#include <algorithm>
#include <chrono>
#include <complex>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <memory>
#include <memory>
#include <mutex>
#include <random>
#include <stdexcept>
#include <thread>
#include <vector>

// HepMC3
#include "HepMC3/WriterAscii.h"
#include "HepMC3/WriterAsciiHepMC2.h"

// Own
#include "Graniitti/MGraniitti.h"
#include "Graniitti/MTimer.h"
#include "Graniitti/MAux.h"

// Libraries
#include "json.hpp"

using gra::aux::indices;
using namespace gra;

// Main
int main(int argc, char *argv[]) {

  MTimer timer(true);

  try {
    if (argc < 3) {  // We expect > 2 arguments: the program name, energies,
      // total number
      // of events
      std::stringstream ss;
      ss << "Usage: ./minbias <ENERGY_0,ENERGY_1,...,ENERGY_K> <EVENTS>";

      aux::CheckUpdate();

      throw std::invalid_argument(ss.str());
    }

    // Input energy list
    std::string         input(argv[1]);
    std::vector<double> sqrtsvec = gra::aux::SplitStr(input, double(0.0));

    // Number of events
    int EVENTS = std::atoi(argv[2]);

    std::cout << "EVENTS: " << EVENTS << std::endl;

    std::vector<std::string> json_in = {"./tests/run_minbias/sd.json", "./tests/run_minbias/dd.json",
                                        "./tests/run_minbias/nd.json"};

    const std::vector<std::string> beam = {"p+", "p+"};

    // Loop over energies
    for (const auto &e : indices(sqrtsvec)) {

      std::vector<double> xs0     = {0, 0};
      std::vector<double> xs0_err = {0, 0};

      double xs_tot = 0.0;
      double xs_el  = 0.0;
      double xs_in  = 0.0;
      
      // Beam and energy
      const std::vector<double> energy = {sqrtsvec[e] / 2, sqrtsvec[e] / 2};

      // Then calculate screened SD and DD integrated cross section
      for (const auto& p : indices(json_in)) {

        // Create generator object first
        std::unique_ptr<MGraniitti> gen = std::make_unique<MGraniitti>();

        // Read process input from file
        gen->ReadInput(json_in[p]);
        //gen->SetNumberOfEvents(0);

        gen->proc->SetInitialState(beam, energy);
        gen->proc->SetScreening(true);

        // ** ALWAYS LAST **
        gen->Initialize();

        // Get process cross sections
        gen->GetXS(xs0[p], xs0_err[p]);

        // Total inelastic
        if (p == 0) {
          gen->proc->Eikonal.GetTotXS(xs_tot, xs_el, xs_in);
        }
      }

      // Non-diffractive = Total_inelastic - (screened_SD + screened_DD);
      const double xs_nd = xs_in - (xs0[0] + xs0[1]);

      // Events for each process
      const int SD_EVT = std::ceil(EVENTS * xs0[0] / xs_in);
      const int DD_EVT = std::ceil(EVENTS * xs0[1] / xs_in);
      int       ND_EVT = std::ceil(EVENTS * xs_nd / xs_in);

      // Make it sure we have exact amount of events
      const int D = EVENTS - (SD_EVT + DD_EVT + ND_EVT);
      ND_EVT -= D;
      std::vector<int> NEVT = {SD_EVT, DD_EVT, ND_EVT};

      // HepMC33
      // outputHepMC33 =
      //    std::make_shared<HepMC3::WriterAscii>("./output/" + OUTPUT +
      //    ".hepmc3");
      //} else if (FORMAT.compare("hepmc2") == 0) {
      // HepMC3::outputHepMC32 =
      // std::shared_ptr<HepMC3::WriterAscii> outputHepMC33;

      // HepMC32
      const std::string OUTPUTNAME =
          "minbias_" + std::to_string(static_cast<int>(sqrtsvec[e]));  // Note x 2
      const std::string                          outputstr = "./output/" + OUTPUTNAME + ".hepmc2";
      std::shared_ptr<HepMC3::WriterAsciiHepMC2> outputHepMC2 =
          std::make_shared<HepMC3::WriterAsciiHepMC2>(outputstr);

      // Loop over processes
      for (std::size_t i = 0; i < NEVT.size(); ++i) {
        
        // Create generator object first
        std::unique_ptr<MGraniitti> gen = std::make_unique<MGraniitti>();

        // Read the process input from a file
        gen->ReadInput(json_in[i]);
        gen->SetNumberOfEvents(NEVT[i]);

        // Set beam and energy (same as above)
        gen->proc->SetInitialState(beam, energy);

        // External HepMC2 output
        gen->SetHepMC2Output(outputHepMC2, OUTPUTNAME);

        // g->proc->SetScreening(false);

        // ** Always Last! **
        gen->Initialize();

        // Force the cross section (total inelastic)
        gen->ForceXS(xs_in);

        // Generate events
        gen->Generate();
      }
      
      // Finalize
      outputHepMC2->close();

      std::cout << "CMS-energy: " << sqrtsvec[e] << " GeV : Generated in total " << EVENTS << " minimum bias events" << std::endl << std::endl;

    }
  } catch (const std::invalid_argument &e) {
    gra::aux::PrintGameOver();
    std::cerr << rang::fg::red << "Exception catched: " << rang::fg::reset << e.what() << std::endl;
    return EXIT_FAILURE;
  } catch (const std::ios_base::failure &e) {
    gra::aux::PrintGameOver();
    std::cerr << rang::fg::red << "Exception catched: std::ios_base::failure: " << rang::fg::reset
              << e.what() << std::endl;
    return EXIT_FAILURE;
  } catch (const nlohmann::json::exception &e) {
    gra::aux::PrintGameOver();
    std::cerr << rang::fg::red << "Exception catched: JSON input: " << rang::fg::reset
              << e.what() << std::endl;
    return EXIT_FAILURE;
  } catch (...) {
    gra::aux::PrintGameOver();
    std::cerr << rang::fg::red << "Exception catched: Unspecified (...) (Probably JSON input)"
              << rang::fg::reset << std::endl;
    return EXIT_FAILURE;
  }
  
  std::cout << "[minbias: done]" << std::endl;
  aux::CheckUpdate();

  return EXIT_SUCCESS;
}
