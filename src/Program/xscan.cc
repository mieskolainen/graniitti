// GRANIITTI - Monte Carlo event generator for high energy diffraction
// https://github.com/mieskolainen/graniitti
//
// <Integrated cross section energy evolution scanner>
//
// (c) 2017-2021 Mikael Mieskolainen
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
#include <mutex>
#include <random>
#include <stdexcept>
#include <thread>
#include <vector>

// OWN
#include "Graniitti/MAux.h"
#include "Graniitti/MGraniitti.h"
#include "Graniitti/MTimer.h"

// Libraries
#include "cxxopts.hpp"
#include "json.hpp"

using gra::aux::indices;
using namespace gra;

// Main
int main(int argc, char *argv[]) {
  aux::PrintArgv(argc, argv);

  MTimer timer;

  // Save the number of input arguments
  const int NARGC = argc - 1;

  try {
    cxxopts::Options options(argv[0], "");

    options.add_options("")("i,input", "Input cards            <card1.json,card2.json,...>",
                            cxxopts::value<std::string>())(
        "e,ENERGY", "CMS energies           <energy0,energy1,...>", cxxopts::value<std::string>())(
        "l,POMLOOP", "Pomeron loop screening <true|false>", cxxopts::value<std::string>())("H,help",
                                                                                           "Help");

    auto r = options.parse(argc, argv);

    if (r.count("help") || NARGC == 0) {
      std::unique_ptr<MGraniitti> gen = std::make_unique<MGraniitti>();
      gen->GetProcessNumbers();

      std::cout << options.help({""}) << std::endl;
      std::cout << rang::style::bold << "Example:" << rang::style::reset << std::endl;
      std::cout << "  " << argv[0] << " -i ./input/test.json -e 500,2760,7000,13000,100000 -l false"
                << std::endl
                << std::endl;
      aux::CheckUpdate();

      return EXIT_FAILURE;
    }

    // Input energy list
    std::vector<double> energy = gra::aux::SplitStr(r["e"].as<std::string>(), double(0.0), ',');

    // Input file list
    std::vector<std::string> jsinput = gra::aux::SplitStr2Str(r["i"].as<std::string>(), ',');

    // Output text file
    FILE *fout;
    fout = fopen("scan.csv", "w");
    if (fout == NULL) {
      printf("Scan:: Error opening scan.csv file! \n");
      return EXIT_FAILURE;
    }
    FILE *fout_latex;
    fout_latex = fopen("scan.tex", "w");
    if (fout_latex == NULL) {
      printf("Scan:: Error opening scan.tex file! \n");
      return EXIT_FAILURE;
    }

    fprintf(fout, "sqrts\txstot\txsin\txsel");
    for (unsigned int k = 0; k < jsinput.size(); ++k) { fprintf(fout, "\txs%d", k); }
    fprintf(fout, "\n");
    fflush(fout);

    // LOOP over energy
    for (const auto &i : indices(energy)) {
      double xs_tot = 0;
      double xs_el  = 0;
      double xs_in  = 0;

      // LOOP over processes
      std::vector<double> xs0(jsinput.size(), 0.0);
      std::vector<double> xs0_err(jsinput.size(), 0.0);
      for (const auto &k : indices(jsinput)) {
        // Create generator object first
        std::unique_ptr<MGraniitti> gen = std::make_unique<MGraniitti>();

        // Read process input from file
        nlohmann::json js = json::parse(gra::aux::GetInputData(jsinput[k]));

        // Re-set parameters
        js["PROCESSPARAM"]["ENERGY"]  = std::vector<double>(2, energy[i] / 2);
        js["PROCESSPARAM"]["POMLOOP"] = (r["l"].as<std::string>() == "true");

        gen->ReadInput(js);

        // Always last!
        gen->Initialize();

        if (k == 0) {  // One process is enough
          gen->proc->Eikonal.GetTotXS(xs_tot, xs_el, xs_in);
        }

        // > Get process cross section and error
        gen->GetXS(xs0[k], xs0_err[k]);
      }

      // Write out
      fprintf(fout, "%0.3E\t%0.3E\t%0.3E\t%0.3E", energy.at(i), xs_tot, xs_in, xs_el);
      fprintf(fout_latex, "%0.3E & %0.3E & %0.3E & %0.3E", energy.at(i), xs_tot, xs_in, xs_el);

      for (unsigned int k = 0; k < jsinput.size(); ++k) {
        fprintf(fout, "\t%0.3E", xs0.at(k));
        fprintf(fout_latex, " & %0.3E", xs0.at(k));
      }
      fprintf(fout, "\n");
      fprintf(fout_latex, " \\\\ \n");

      fflush(fout);        // flush it out
      fflush(fout_latex);  // flush it out
    }
    fclose(fout);
    fclose(fout_latex);
  } catch (const std::invalid_argument &e) {
    std::unique_ptr<MGraniitti> gen = std::make_unique<MGraniitti>();
    gen->GetProcessNumbers();

    gra::aux::PrintGameOver();
    std::cerr << rang::fg::red << "Exception catched: " << rang::fg::reset << e.what() << std::endl;
    return EXIT_FAILURE;
  } catch (const std::ios_base::failure &e) {
    gra::aux::PrintGameOver();
    std::cerr << rang::fg::red << "Exception catched: std::ios_base::failure: " << rang::fg::reset
              << e.what() << std::endl;
    return EXIT_FAILURE;
  } catch (const cxxopts::OptionException &e) {
    gra::aux::PrintGameOver();
    std::cerr << rang::fg::red << "Exception catched: Commandline options: " << rang::fg::reset
              << e.what() << std::endl;
    return EXIT_FAILURE;
  } catch (const nlohmann::json::exception &e) {
    gra::aux::PrintGameOver();
    std::cerr << rang::fg::red << "Exception catched: JSON input: " << rang::fg::reset << e.what()
              << std::endl;
    return EXIT_FAILURE;
  } catch (...) {
    std::unique_ptr<MGraniitti> gen = std::make_unique<MGraniitti>();
    gen->GetProcessNumbers();

    gra::aux::PrintGameOver();
    std::cerr << rang::fg::red << "Exception catched: Unspecified (...) (Probably JSON input)"
              << rang::fg::reset << std::endl;
    return EXIT_FAILURE;
  }

  printf("\n");
  printf("scan:: Finished in %0.1f sec \n", timer.ElapsedSec());
  printf("scan:: Output created to scan{.csv,.tex} \n\n");

  std::cout << "[xscan: done]" << std::endl;

  return EXIT_SUCCESS;
}
