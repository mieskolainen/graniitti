// Main event generator program
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <chrono>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

// Own
#include "Graniitti/MGraniitti.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MJsonOverride.h"

// Libraries
#include "cxxopts.hpp"
#include "json.hpp"
#include "rang.hpp"

using gra::aux::indices;

using namespace gra;

namespace {

// Collect all repeated JSON override option values in original command line
// order
std::vector<std::string>
CollectJsonOverrideSpecs(const cxxopts::ParseResult &r) {
  std::vector<std::string> specs;
  for (const auto &arg : r.arguments()) {
    if (arg.key() == "set") {
      specs.push_back(arg.value());
    }
  }
  return specs;
}

} // namespace

// Main
int main(int argc, char *argv[]) {
  aux::PrintArgv(argc, argv);

  // Save the number of input arguments
  const int NARGC = argc - 1;

  // Create generator object first
  std::unique_ptr<MGraniitti> gen = std::make_unique<MGraniitti>();

  try {
    cxxopts::Options options(argv[0], "");

    options.add_options("")("i,INPUT", "Input card                <string>",
                            cxxopts::value<std::string>())(
        "d,VGRID", "Use pre-computed integration proposal <string>",
        cxxopts::value<std::string>())("H,help", "Help");

    // Add global options
    gen->ConstructTerminal(options);

    auto r = options.parse(argc, argv);
    const std::vector<gra::json_override::OverrideSpec> json_overrides =
        gra::json_override::ParseSpecs(CollectJsonOverrideSpecs(r));

    if (r.count("help") || NARGC == 0) {
      gen->GetProcessNumbers(r.count("m") ? r["m"].as<std::string>() : gra::MODELPARAM);
      std::cout << options.help({"", "GENERIC", "SCATTERING"}) << std::endl;
      std::cout << rang::style::bold
                << " Arrow operators:" << rang::style::reset << std::endl;
      std::cout << "  use -> between initial and final states" << std::endl;
      std::cout << "  use &> (instead of ->) for a decoupled central system "
                   "phase space (use with "
                   "<F>, <P> class, e.g. for s-channel resonances)"
                << std::endl;
      std::cout << "  use  > for recursive decaytree branchings with curly "
                   "brackets { grand daughters }"
                << std::endl;
      std::cout << std::endl;
      std::cout << rang::style::bold
                << " Inline 'on-the-flight' parameters to concatenate with "
                   "PROCESS string:"
                << rang::style::reset << std::endl
                << std::endl;

      std::cout << rang::style::bold << "  [Generic]" << rang::style::reset
                << std::endl;

      std::cout << "  @FLATAMP:N                              flat matrix "
                   "element for 'pure' phase "
                   "space generation, set N to -1 for more info"
                << std::endl;
      std::cout << "  @FLATMASS2:true                         flat sampling in "
                   "M^2 instead of "
                   "relativistic Breit-Wigner f(M^2) in decay trees "
                   "(true/false or 1/0)"
                << std::endl;
      std::cout << "  @OFFSHELL:X                             how many +- full "
                   "widths particles "
                   "off-shell in decay trees (X = 0 on-shell, default from "
                   "NUMERICS_MC.offshell_max)"
                << std::endl;
      std::cout << "  @PDG[X]{M:350.0, W:5.0}                 new mass and "
                   "width for pdg particle id X"
                << std::endl;
      std::cout << "  @j={u,d,s,c,b,g}                        outgoing generic "
                   "parton species"
                << std::endl;
      std::cout << "  @R[f0_980]{M:0.990, W:0.065}            set new central "
                   "resonance mass and width"
                << std::endl;
      std::cout << "  @RES{rho_770:1, ..., f2_1270:1}         set active "
                   "resonances in the "
                   "amplitude (true/false or 1/0)"
                << std::endl;

      std::cout << std::endl;
      std::cout << rang::style::bold << "  [Pomeron amplitudes]"
                << rang::style::reset << std::endl;

      std::cout << "  @SPINGEN:true                           set generation "
                   "2->1 spin "
                   "correlations active (true/false or 1/0)"
                << std::endl;
      std::cout << "  @SPINDEC:true                           set decay 1->2 "
                   "spin correlations "
                   "active (true/false or 1/0)"
                << std::endl;

      std::cout << "  @QMETRICS:true                          print integrated quantum metrics (true/false)"
                << std::endl;
      std::cout << "  @MP_FRAME:X                             set MP central "
                   "spin basis frame "
                   "(X = HX, CS, CM)"
                << std::endl;
      std::cout << "  @R[f2_1270]{JZ0:0.5, JZ1:0.0, JZ2:0.5}  set Jz populations "
                   "for MP[RES] / MP[RES+CON] "
                   "(polarization.mode selects a_Jz or diagonal rho)"
                << std::endl;
      std::cout << "  @MMAX:X                                 set maximum "
                   "analytic Regge helicity for GP"
                << std::endl;

      std::cout << std::endl;
      std::cout << rang::style::bold << "  [Tensor pomeron amplitudes]"
                << rang::style::reset << std::endl;

      std::cout
          << "  @R[f0_980]{g0:1.0, g1:0.2, ...}         set new production "
             "couplings {g0,g1} "
             "[scalar/pseudoscalar/vector/axial-vector] {g0,...,g6} [tensor]"
          << std::endl;

      std::cout << std::endl;
      std::cout << rang::style::bold
                << " PROCESS string examples:" << rang::style::reset
                << std::endl;

      std::cout << "  yy[EPA]<C> -> mu+ mu-" << std::endl;
      std::cout << "  yy[EPA]<F> -> mu+ mu-" << std::endl;
      std::cout << "  gg[QCD]<F> -> j ~j @j={u,d,s,c,b,g}" << std::endl;
      std::cout << "  IPp[Zj]<F> -> Z > {mu+ mu-} j @j={u,d,s,c,b,g}"
                << std::endl;
      std::cout << "  yy[Higgs]<F> &> 22 22" << std::endl;
      std::cout << "  yy[EPA]<C> -> 882 -882 @PDG[882]{M:1500,W:0}"
                << std::endl;
      std::cout << "  TP[RES]<F> -> pi+ pi- @RES{rho_770:1}" << std::endl;
      std::cout << "  MP[CON]<F> -> rho(770)0 > {pi+ pi-} rho(770)0 > {pi+ pi-}"
                << std::endl;
      std::cout << "  MP[CON]<F> -> pi+ pi- pi+ pi-" << std::endl;
      std::cout << "  MP[CON]<F> -> p+ p- @FLATAMP:2" << std::endl;
      std::cout << "  MP[RES+CON]<F> -> pi+ pi- "
                   "@RES{f0_500:0,rho_770:1,f0_980:1,f2_1270:1} "
                   "@R[f0_980]{M:0.98,W:0.065}"
                << std::endl;
      std::cout << std::endl;

      std::cout << rang::style::bold << rang::fg::green
                << " A steering card example with no additional input:"
                << rang::style::reset << std::endl;
      std::cout << "  " << argv[0] << " -i gencard/test.json" << std::endl
                << std::endl;

      std::cout << rang::style::bold << rang::fg::red
                << " A steering card example with commandline input override:"
                << rang::style::reset << std::endl;
      std::cout << "  " << argv[0]
                << " -i gencard/test.json -p \"yy[QED]<F> -> e+ e-\""
                << std::endl
                << std::endl;

      std::cout << rang::style::bold << rang::fg::red
                << " Generic JSON card override examples:" << rang::style::reset
                << std::endl;
      std::cout << "  " << argv[0]
                << " -i gencard/test.json --set 'GENCUTS.<F>.M=[0.5,2.0]'"
                << std::endl;
      std::cout
          << "  " << argv[0]
          << " -i gencard/test.json --set "
             "'GENERAL.json:PARAM_SOFT.MODEL.single.EXCHANGE.P.g[0,0]=8.4'"
          << std::endl
          << std::endl;

      std::cout << rang::style::bold << rang::fg::green
                << " A steering card example with pre-computed MC integration "
                   "array and custom "
                   "random seed (for e.g. GRID computing):"
                << rang::style::reset << std::endl;
      std::cout << "  " << argv[0]
                << " -i gencard/test.json -d vgrid/test.vgrid -r 12345"
                << std::endl
                << std::endl;

      aux::CheckUpdate();

      return EXIT_SUCCESS;
    }

    // ===================================================================

    // Read and parse json input
    const std::string inputfile = r["i"].as<std::string>();
    std::cout << "gr: Reading input: " << inputfile << std::endl;
    json j = json::parse(gra::aux::GetInputData(inputfile));

    // Process terminal input
    gen->ProcessTerminal(j, r);

    // Process generic JSON command line overrides for the input card
    gra::json_override::ApplyInputOverrides(j, json_overrides);
    const std::vector<std::string> override_cards =
        gra::json_override::ResolveCardOverrideTargets(
            json_overrides,
            gra::ResolveModelTuneDir(j.at("GENERIC").at("MODELPARAM")));
    gra::json_override::RegisterCardOverrides(json_overrides);
    for (const auto &card : override_cards) {
      (void)gra::aux::GetInputData(card);
    }
    gra::json_override::RequireAllCardOverridesApplied();

    // Read input parameters
    gen->ReadInput(j);
    gen->proc->SetDebug(r.count("debug") != 0);

    // -------------------------------------------------------------------
    // Use a pre-computed integration proposal
    if (r.count("d")) {
      gen->SetVgridFile(r["d"].as<std::string>());
      gen->proc->SetHistograms(0); // Always off when using pre-computed grid!
    }
    // -------------------------------------------------------------------

    // Initialize generator (always last!)
    gen->Initialize();

    // Generate events
    gen->Generate();

    // Print histograms
    gen->PrintHistograms();

    // Report the original failure without querying partially configured processes
  } catch (const std::invalid_argument &e) {
    gra::aux::PrintGameOver();
    std::cerr << rang::fg::red << "Exception catched: " << rang::fg::reset
              << e.what() << std::endl;
    return EXIT_FAILURE;
  } catch (const std::ios_base::failure &e) {
    gra::aux::PrintGameOver();
    std::cerr << rang::fg::red << "Exception catched: std::ios_base::failure: "
              << rang::fg::reset << e.what() << std::endl;
    return EXIT_FAILURE;
  } catch (const cxxopts::OptionException &e) {
    gra::aux::PrintGameOver();
    std::cerr << rang::fg::red
              << "Exception catched: Commandline options: " << rang::fg::reset
              << e.what() << std::endl;
    return EXIT_FAILURE;
  } catch (const nlohmann::json::exception &e) {
    gra::aux::PrintGameOver();
    std::cerr << rang::fg::red
              << "Exception catched: JSON input: " << rang::fg::reset
              << e.what() << std::endl;
    return EXIT_FAILURE;
  } catch (const std::exception &e) {
    gra::aux::PrintGameOver();
    std::cerr << rang::fg::red
              << "Exception catched: std::exception: " << rang::fg::reset
              << e.what() << std::endl;
    return EXIT_FAILURE;
  } catch (...) {
    gra::aux::PrintGameOver();
    std::cerr << rang::fg::red << "Exception catched: non-standard exception"
              << rang::fg::reset << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "[gr: done]" << std::endl;

  return EXIT_SUCCESS;
}
