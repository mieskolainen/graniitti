// Scan integrated cross sections across collision energies
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
#include <iomanip>
#include <iostream>
#include <memory>
#include <mutex>
#include <random>
#include <stdexcept>
#include <thread>
#include <vector>

// OWN
#include "Graniitti/MGraniitti.h"
#include "Graniitti/Program/MInput.h"
#include "Graniitti/Program/xscan.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MJsonOverride.h"
#include "Graniitti/Tech/MTimer.h"

// Libraries
#include "cxxopts.hpp"
#include "json.hpp"

using gra::aux::indices;
using namespace gra;

namespace {

// Collect all repeated JSON override option values in original command line order
std::vector<std::string> CollectJsonOverrideSpecs(const cxxopts::ParseResult& r) {
  std::vector<std::string> specs;
  for (const auto& arg : r.arguments()) {
    if (arg.key() == "set") { specs.push_back(arg.value()); }
  }
  return specs;
}

}  // namespace

// Main
int main(int argc, char* argv[]) {
  aux::PrintArgv(argc, argv);

  MTimer timer;

  // Save the number of input arguments
  const int NARGC = argc - 1;

  try {
    cxxopts::Options options(argv[0], "");

    options.add_options("")("i,input", "Input cards            <card1.json,card2.json,...>",
                            cxxopts::value<std::string>())("e,ENERGY", "CMS energies           <energy0,energy1,...>",
                                                           cxxopts::value<std::string>())(
        "l,LOOPSCREEN", "Soft survival screening <true|false|1|0>", cxxopts::value<std::string>())(
        "set", "Override JSON card entry <[card:]path=json>", cxxopts::value<std::string>())("H,help", "Help");

    auto                                                r = options.parse(argc, argv);
    const std::vector<gra::json_override::OverrideSpec> json_overrides =
        gra::json_override::ParseSpecs(CollectJsonOverrideSpecs(r));

    if (r.count("help") || NARGC == 0) {
      std::unique_ptr<MGraniitti> gen = std::make_unique<MGraniitti>();
      gen->GetProcessNumbers();

      std::cout << options.help({""}) << std::endl;
      std::cout << rang::style::bold << "Example:" << rang::style::reset << std::endl;
      std::cout << "  " << argv[0] << " -i gencard/test.json -e 500,2760,7000,13000,100000 -l false" << std::endl;
      std::cout << "  " << argv[0]
                << " -i gencard/test.json -e 13000 -l true --set 'GENERAL.json:PARAM_SOFT.active_model=\"double\"'"
                << std::endl
                << std::endl;
      aux::CheckUpdate();

      return EXIT_FAILURE;
    }

    // Input energy list
    std::vector<double> energy = gra::program::Energies(r["e"].as<std::string>());

    // Input file list
    std::vector<std::string> jsinput = gra::aux::SplitStr2Str(r["i"].as<std::string>(), ',');

    if (jsinput.empty()) { throw std::invalid_argument("xscan: at least one input card is required"); }
    const bool loopscreen = aux::ParseBool(r["l"].as<std::string>(), "--LOOPSCREEN");

    // Check writes, flushes and close operations on both scan outputs
    for (const auto& file : jsinput) {
      gra::program::DistinctFiles(file, "scan.csv");
      gra::program::DistinctFiles(file, "scan.tex");
    }
    std::ofstream fout, fout_latex;
    fout.exceptions(std::ios::failbit | std::ios::badbit);
    fout_latex.exceptions(std::ios::failbit | std::ios::badbit);
    fout.open("scan.csv");
    fout_latex.open("scan.tex");
    fout << std::scientific << std::uppercase << std::setprecision(4);
    fout_latex << std::scientific << std::uppercase << std::setprecision(4);
    fout << "sqrts\txstot\txsin\txsel";
    for (const auto& k : indices(jsinput)) { fout << "\txs" << k; }
    for (const auto& k : indices(jsinput)) { fout << "\txs" << k << "_err"; }
    fout << std::endl;

    // LOOP over energy
    for (const auto& i : indices(energy)) {
      double xs_tot = 0;
      double xs_el  = 0;
      double xs_in  = 0;

      // LOOP over processes
      std::vector<double> xs0(jsinput.size(), 0.0);
      std::vector<double> xs0_err(jsinput.size(), 0.0);
      for (const auto& k : indices(jsinput)) {
        gra::json_override::ClearCardOverrides();

        // Create generator object first
        std::unique_ptr<MGraniitti> gen = std::make_unique<MGraniitti>();

        // Read process input from file
        nlohmann::json js = json::parse(gra::aux::GetInputData(jsinput[k]));

        // Process generic JSON command line overrides for the input card and model cards
        gra::json_override::ApplyInputOverrides(js, json_overrides);
        const std::vector<std::string> override_cards = gra::json_override::ResolveCardOverrideTargets(
            json_overrides, gra::ResolveModelTuneDir(js.at("GENERIC").at("MODELPARAM")));
        gra::json_override::RegisterCardOverrides(json_overrides);
        for (const auto& card : override_cards) { (void)gra::aux::GetInputData(card); }
        gra::json_override::RequireAllCardOverridesApplied();

        // Re-set parameters
        js["SCATTERING"]["ENERGY"]  = std::vector<double>(2, energy[i] / 2);
        js["SCATTERING"]["LOOPSCREEN"] = loopscreen;

        gen->ReadInput(js);

        // Always last!
        gen->Initialize();
        gra::json_override::ClearCardOverrides();

        if (k == 0) {  // One process is enough
          const auto xs = gra::program::ScanXS(gen->proc->GetEikonal());
          xs_tot        = xs[0];
          xs_el         = xs[1];
          xs_in         = xs[2];
          if (!std::isfinite(xs_tot)) {
            std::cerr << "xscan: soft total cross sections unavailable for the first process, writing nan" << std::endl;
          }
        }

        // > Get process cross section and error
        gen->GetXS(xs0[k], xs0_err[k]);
      }

      // Write each completed energy point before starting the next integration
      fout << energy[i] << '\t' << xs_tot << '\t' << xs_in << '\t' << xs_el;
      fout_latex << energy[i] << " & " << xs_tot << " & " << xs_in << " & " << xs_el;
      for (const auto& k : indices(jsinput)) {
        fout << '\t' << xs0[k];
        fout_latex << " & " << xs0[k];
      }
      for (const auto& k : indices(jsinput)) { fout << '\t' << xs0_err[k]; }
      fout << std::endl;
      fout_latex << " \\\\ " << std::endl;
    }
    fout.close();
    fout_latex.close();
  } catch (const std::invalid_argument& e) {
    std::unique_ptr<MGraniitti> gen = std::make_unique<MGraniitti>();
    gen->GetProcessNumbers();

    gra::aux::PrintGameOver();
    std::cerr << rang::fg::red << "Exception catched: " << rang::fg::reset << e.what() << std::endl;
    return EXIT_FAILURE;
  } catch (const std::ios_base::failure& e) {
    gra::aux::PrintGameOver();
    std::cerr << rang::fg::red << "Exception catched: std::ios_base::failure: " << rang::fg::reset << e.what()
              << std::endl;
    return EXIT_FAILURE;
  } catch (const cxxopts::OptionException& e) {
    gra::aux::PrintGameOver();
    std::cerr << rang::fg::red << "Exception catched: Commandline options: " << rang::fg::reset << e.what()
              << std::endl;
    return EXIT_FAILURE;
  } catch (const nlohmann::json::exception& e) {
    gra::aux::PrintGameOver();
    std::cerr << rang::fg::red << "Exception catched: JSON input: " << rang::fg::reset << e.what() << std::endl;
    return EXIT_FAILURE;
  } catch (const std::exception& e) {
    std::unique_ptr<MGraniitti> gen = std::make_unique<MGraniitti>();
    gen->GetProcessNumbers();

    gra::aux::PrintGameOver();
    std::cerr << rang::fg::red << "Exception catched: std::exception: " << rang::fg::reset << e.what() << std::endl;
    return EXIT_FAILURE;
  } catch (...) {
    std::unique_ptr<MGraniitti> gen = std::make_unique<MGraniitti>();
    gen->GetProcessNumbers();

    gra::aux::PrintGameOver();
    std::cerr << rang::fg::red << "Exception catched: non-standard exception" << rang::fg::reset << std::endl;
    return EXIT_FAILURE;
  }

  printf("\n");
  printf("scan:: Finished in %0.1f sec \n", timer.ElapsedSec());
  printf("scan:: Output created to scan{.csv,.tex} \n\n");

  std::cout << "[xscan: done]" << std::endl;

  return EXIT_SUCCESS;
}
