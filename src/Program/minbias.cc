// Generate a screened minimum bias mixture
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <math.h>

#include <algorithm>
#include <chrono>
#include <complex>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <iomanip>
#include <iostream>
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
#include "Graniitti/Program/MInput.h"
#include "Graniitti/Program/minbias.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MTimer.h"

// Libraries
#include "json.hpp"

using gra::aux::indices;
using namespace gra;


// Main
int main(int argc, char* argv[]) {
  aux::PrintArgv(argc, argv);

  MTimer timer(true);

  try {
    if (argc != 3) {
      std::stringstream ss;
      ss << "Usage: ./minbias <ENERGY_0,ENERGY_1,...,ENERGY_K> <EVENTS>";

      aux::CheckUpdate();

      throw std::invalid_argument(ss.str());
    }

    // Input energy list
    const auto sqrtsvec = gra::program::Energies(argv[1]);

    // Number of events
    const int EVENTS = gra::program::Number<int>(argv[2]);
    if (EVENTS < 0) { throw std::invalid_argument("minbias: event count must be nonnegative"); }

    std::cout << "EVENTS: " << EVENTS << std::endl;

    std::vector<std::string> json_in = {"./tests/physics/studies/minbias/gencard_sd.json",
                                        "./tests/physics/studies/minbias/gencard_dd.json",
                                        "./tests/physics/studies/minbias/gencard_nd.json"};

    const std::vector<std::string> beam = {"p+", "p+"};

    // Loop over energies
    for (const auto& e : indices(sqrtsvec)) {
      std::vector<double> xs0     = {0, 0};
      std::vector<double> xs0_err = {0, 0};

      double xs_tot = 0.0;
      double xs_el  = 0.0;
      double xs_in  = 0.0;

      // Beam and energy
      const std::vector<double> energy = {sqrtsvec[e] / 2, sqrtsvec[e] / 2};

      // Use the same screened cards in integration and event generation
      std::vector<json> cards;
      for (const auto& file : json_in) {
        auto card                     = json::parse(gra::aux::GetInputData(file));
        card["SCATTERING"]["BEAM"]    = beam;
        card["SCATTERING"]["ENERGY"]  = energy;
        card["SCATTERING"]["LOOPSCREEN"] = true;
        cards.push_back(std::move(card));
      }

      // Then calculate screened SD and DD integrated cross section
      for (const auto& p : indices(xs0)) {
        // Create generator object first
        std::unique_ptr<MGraniitti> gen = std::make_unique<MGraniitti>();

        gen->ReadInput(cards[p]);
        gen->SetNumberOfEvents(0);

        // ** ALWAYS LAST **
        gen->Initialize();

        // Get process cross sections
        gen->GetXS(xs0[p], xs0_err[p]);

        // Total inelastic
        if (p == 0) { gen->proc->eikonal.GetTotXS(xs_tot, xs_el, xs_in); }
      }

      // Non-diffractive = Total_inelastic - (screened_SD + screened_DD);
      const double xs_nd = xs_in - (xs0[0] + xs0[1]);

      // Preserve the total event count and nonnegative component populations
      const auto NEVT = gra::program::MinbiasCounts(EVENTS, {xs0[0], xs0[1], xs_nd});

      // Preserve fractional energies in output names without an integer conversion
      std::ostringstream name;
      name << "minbias_" << std::setprecision(std::numeric_limits<double>::max_digits10) << sqrtsvec[e];
      const std::string OUTPUTNAME = name.str();
      gra::aux::CreateDirectory("output");
      const std::string outputstr = "./output/" + OUTPUTNAME + ".hepmc2";
      std::ofstream target(outputstr);
      if (!target.is_open()) { throw std::runtime_error("minbias: cannot open " + outputstr); }
      auto outputHepMC2 = std::make_shared<HepMC3::WriterAsciiHepMC2>(target);

      // Loop over processes
      for (const auto& p : indices(NEVT)) {
        if (NEVT[p] == 0) { continue; }
        // Create generator object first
        std::unique_ptr<MGraniitti> gen = std::make_unique<MGraniitti>();

        gen->ReadInput(cards[p]);
        gen->SetNumberOfEvents(NEVT[p]);

        // External HepMC2 output
        gen->SetHepMC2Output(outputHepMC2, OUTPUTNAME);

        // ** Always Last! **
        gen->Initialize();

        // Force the cross section (total inelastic)
        gen->ForceXS(xs_in);

        // Generate events
        gen->Generate();
      }

      // Finalize
      outputHepMC2->close();
      target.close();
      if (outputHepMC2->failed() || target.fail()) { throw std::runtime_error("minbias: failed event output"); }

      aux::PrintBar("=");
      std::cout << "CMS-energy: " << sqrtsvec[e] << " GeV" << std::endl;
      std::cout << std::endl;
      std::cout << "Screened Cross Sections within space phase applied:" << std::endl << std::endl;

      printf(" Single Diffractive (forward+backward): %0.3f mb\n", xs0[0] * 1e3);
      printf(" Double Diffractive:                    %0.3f mb\n", xs0[1] * 1e3);
      printf(" Non-Diffractive:                       %0.3f mb\n", xs_nd * 1e3);
      printf(" Total inelastic:                       %0.3f mb\n", xs_in * 1e3);

      std::cout << std::endl;
      std::cout << "Generated in total " << EVENTS << " minimum bias events according to xs above" << std::endl
                << std::endl;
      aux::PrintBar("=");
    }
  } catch (const std::invalid_argument& e) {
    gra::aux::PrintGameOver();
    std::cerr << rang::fg::red << "Exception catched: " << rang::fg::reset << e.what() << std::endl;
    return EXIT_FAILURE;
  } catch (const std::ios_base::failure& e) {
    gra::aux::PrintGameOver();
    std::cerr << rang::fg::red << "Exception catched: std::ios_base::failure: " << rang::fg::reset << e.what()
              << std::endl;
    return EXIT_FAILURE;
  } catch (const nlohmann::json::exception& e) {
    gra::aux::PrintGameOver();
    std::cerr << rang::fg::red << "Exception catched: JSON input: " << rang::fg::reset << e.what() << std::endl;
    return EXIT_FAILURE;
  } catch (const std::exception& e) {
    gra::aux::PrintGameOver();
    std::cerr << rang::fg::red << "Exception catched: std::exception: " << rang::fg::reset << e.what() << std::endl;
    return EXIT_FAILURE;
  } catch (...) {
    gra::aux::PrintGameOver();
    std::cerr << rang::fg::red << "Exception catched: non-standard exception" << rang::fg::reset << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "[minbias: done]" << std::endl;
  aux::CheckUpdate();

  return EXIT_SUCCESS;
}
