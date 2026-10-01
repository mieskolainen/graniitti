// ROOT-based fiducial observable analyzer
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <math.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <vector>

// ROOT
#include "TBranch.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TMinuit.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TRandom3.h"
#include "TString.h"
#include "TStyle.h"
#include "TTree.h"

// Own
#include "Graniitti/Analysis/MAnalyzer.h"
#include "Graniitti/Program/Analysis/analyze.h"
#include "Graniitti/Analysis/MMultiplet.h"
#include "Graniitti/Analysis/MROOT.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MMatrix.h"

// Libraries
#include "cxxopts.hpp"
#include "json.hpp"
#include "rang.hpp"

using gra::aux::indices;
using namespace gra;

// Initialize 1D-histograms
//
void Init1DHistogram(std::map<std::string, std::shared_ptr<h1Multiplet>> &h,
                     const std::vector<std::string> &legendtext, std::vector<int> multiplicity,
                     const std::string &title, const std::string &units, const h1Bound &bM,
                     const h1Bound &bP, const h1Bound &bY) {
  std::string       name = "null";
  const std::string U    = (units != "events") ? "#sigma" : "N";

  // Central system observables
  name    = "h1_S_M";
  h[name] = std::make_shared<h1Multiplet>(
      name, title + ";System M  (GeV);d" + U + "/dM  (" + units + "/GeV)", bM.N, bM.min, bM.max,
      legendtext);

  name    = "h1_S_Pt";
  h[name] = std::make_shared<h1Multiplet>(
      name, title + ";System P_{T} (GeV);d" + U + "/dP_{T}  (" + units + "/GeV)", bP.N, bP.min,
      bP.max, legendtext);

  name    = "h1_S_Pt2";
  h[name] = std::make_shared<h1Multiplet>(
      name, title + ";System P_{T}^{2} (GeV^{2});d" + U + "/dP_{T}^{2}  (" + units + "/GeV^{2})",
      bP.N, math::pow2(bP.min), math::pow2(bP.max), legendtext);

  name    = "h1_S_Y";
  h[name] = std::make_shared<h1Multiplet>(name, title + ";System Y;d" + U + "/dY  (" + units + ")",
                                          bY.N, bY.min, bY.max, legendtext);

  // Central track observables
  name    = "h1_1B_pt";
  h[name] = std::make_shared<h1Multiplet>(
      name, title + ";Central final state p_{T} (GeV);d" + U + "/dp_{T}  (" + units + "/GeV)", bP.N,
      bP.min, bP.max, legendtext);

  name    = "h1_1B_eta";
  h[name] = std::make_shared<h1Multiplet>(
      name, title + ";Central final state #eta;d" + U + "/d#eta  (" + units + ")", bY.N, bY.min,
      bY.max, legendtext);

  // 2-Body observables
  if (std::find(multiplicity.begin(), multiplicity.end(), 2) != multiplicity.end()) {
    for (const auto &i : analyzer::Frames()) {
      name = "h1_costheta_" + i;
      h[name] =
          std::make_shared<h1Multiplet>(name,
                                        title + ";Central final state cos(#theta) [" + i +
                                            " frame];d" + U + "/dcos(#theta)  (" + units + ")",
                                        bP.N, -1.0, 1.0, legendtext);

      name    = "h1_phi_" + i;
      h[name] = std::make_shared<h1Multiplet>(
          name,
          title + ";Central final state #phi [" + i + " frame];d" + U + "/d#phi  (" + units + ")",
          bP.N, -math::PI, math::PI, legendtext);
    }

    name    = "h1_2B_acop";
    h[name] = std::make_shared<h1Multiplet>(name,
                                            title +
                                                ";Central final state acoplanarity #rho = 1 - "
                                                "|#delta#phi|/#pi;d" +
                                                U + "/d#rho  (" + units + "/rad)",
                                            100, 0.0, 1.0, legendtext);

    name    = "h1_2B_diffrap";
    h[name] = std::make_shared<h1Multiplet>(
        name, title + ";#deltay #equiv y_{1} - y_{2};d" + U + "/d#deltay  (" + units + ")", bY.N,
        bY.min - bY.max, bY.max - bY.min, legendtext);
  }

  // 4-Body observables
  if (std::find(multiplicity.begin(), multiplicity.end(), 4) != multiplicity.end()) {
    // ...
  }

  // Forward proton observables
  name    = "h1_PP_dphi";
  h[name] = std::make_shared<h1Multiplet>(
      name, title + ";Proton pair #delta#phi (rad);d" + U + "/d#delta#phi  (" + units + "/rad)",
      100, 0.0, 3.14159, legendtext);

  name    = "h1_PP_t1";
  h[name] = std::make_shared<h1Multiplet>(
      name, title + ";Mandelstam -t_{1} (GeV^{2});d" + U + "/dt  (" + units + "/GeV^{2})", bP.N,
      bP.min, bP.max, legendtext);

  name    = "h1_PP_dpt";
  h[name] = std::make_shared<h1Multiplet>(name,
                                          title +
                                              ";Proton pair |#delta#bar{p}_{T}| "
                                              "(GeV);d" +
                                              U + "/d|#delta#bar{p}_{T}|  (" + units + "/GeV)",
                                          bP.N, bP.min, bP.max, legendtext);
}

// Initialize 2D-histograms
//
void Init2DHistogram(std::map<std::string, std::shared_ptr<h2Multiplet>> &h,
                     const std::vector<std::string> &legendtext, std::vector<int> multiplicity,
                     const std::string &title, const std::string &units, const h1Bound &bM,
                     const h1Bound &bP, const h1Bound &bY) {
  std::string       name = "null";
  const std::string U    = (units != "events") ? "#sigma" : "N";

  // Central system observables
  name    = "h2_S_M_Pt";
  h[name] = std::make_shared<h2Multiplet>(name,
                                          "d" + U + "^2/dMdP_{T}  (" + units + "/GeV/GeV) | " +
                                              title + ";System M (GeV); System P_{T} (GeV)",
                                          bM.N, bM.min, bM.max, bP.N, bP.min, bP.max, legendtext);

  name    = "h2_S_M_t";
  h[name] = std::make_shared<h2Multiplet>(
      name,
      "d" + U + "^2/dMd#t  (" + units + "/GeV/GeV) | " + title + ";System M (GeV); |t| (GeV^{2})",
      bM.N, bM.min, bM.max, bP.N, math::pow2(bP.min), math::pow2(bP.max), legendtext);

  name = "h2_S_M_pt";
  h[name] =
      std::make_shared<h2Multiplet>(name,
                                    "d" + U + "^2/dMdp_{T}  (" + units + "/GeV/GeV) | " + title +
                                        ";System M (GeV); Central final state p_{T} (GeV)",
                                    bM.N, bM.min, bM.max, bP.N, bP.min, bP.max, legendtext);

  name = "h2_S_M_dphipp";
  h[name] =
      std::make_shared<h2Multiplet>(name,
                                    "d" + U + "^2/dMd#delta#phi_{pp}  (" + units + "/GeV/rad) | " +
                                        title + ";System M (GeV); Forward proton #delta#phi_{pp}",
                                    bM.N, bM.min, bM.max, 100, 0.0, gra::math::PI, legendtext);

  name    = "h2_S_M_dpt";
  h[name] = std::make_shared<h2Multiplet>(
      name,
      "d" + U + "^2/dMd|#delta#bar{p}_{T}|  (" + units + "/GeV/GeV) | " + title +
          ";System M (GeV); Proton pair |#delta#bar{p}_{T}| (GeV)",
      bM.N, bM.min, bM.max, 100, 0.0, 2.0, legendtext);

  // 2-Body
  if (std::find(multiplicity.begin(), multiplicity.end(), 2) != multiplicity.end()) {
    name    = "h2_2B_M_dphi";
    h[name] = std::make_shared<h2Multiplet>(
        name,
        "d" + U + "^2/dMd#delta#phi  (" + units + "/GeV/rad) | " + title +
            ";System M (GeV); Central final state #delta#phi (rad)",
        bM.N, bM.min, bM.max, 100, 0.0, gra::math::PI, legendtext);

    name    = "h2_2B_eta1_eta2";
    h[name] = std::make_shared<h2Multiplet>(
        name, "d" + U + "^2/d#eta_{1}d#eta_{2}  (" + units + ") | " + title + ";#eta_{1}; #eta_{2}",
        bY.N, bY.min, bY.max, bY.N, bY.min, bY.max, legendtext);

    // 2D (costheta, phi) in different rest frames
    for (const auto &i : analyzer::Frames()) {
      name    = "h2_2B_costheta_phi_" + i;
      h[name] = std::make_shared<h2Multiplet>(
          name,
          "d" + U + "^2/dcos(#theta)d#phi  (" + units + ") | " + title +
              ";daughter cos(#theta); daughter #phi (rad) [" + i + " FRAME]",
          100, -1, 1, 100, -gra::math::PI, gra::math::PI, legendtext);
    }

    // 2D (M, costheta) in different rest frames
    for (const auto &i : analyzer::Frames()) {
      name    = "h2_2B_M_costheta_" + i;
      h[name] = std::make_shared<h2Multiplet>(
          name,
          "d" + U + "^2/dMdcos(#theta)  (" + units + "/GeV) | " + title +
              ";M (GeV); daughter cos(#theta) [" + i + " FRAME]",
          bM.N, bM.min, bM.max, 100, -1, 1, legendtext);
    }

    // 2D (M, phi) in different rest frames
    for (const auto &i : analyzer::Frames()) {
      name    = "h2_2B_M_phi_" + i;
      h[name] = std::make_shared<h2Multiplet>(
          name,
          "d" + U + "^2/dMd#phi  (" + units + "/GeV/rad) | " + title +
              ";M (GeV); daughter #phi (rad) [" + i + " FRAME]",
          bM.N, bM.min, bM.max, 100, -gra::math::PI, gra::math::PI, legendtext);
    }
  }

  // 4-Body observables
  if (std::find(multiplicity.begin(), multiplicity.end(), 4) != multiplicity.end()) {
    // ...
  }
}

// Main program
int main(int argc, char *argv[]) {
  gra::aux::PrintArgv(argc, argv);
  gra::rootstyle::SetROOTStyle();

  gra::aux::PrintFlashScreen(rang::fg::magenta);
  std::cout << rang::style::bold << "GRANIITTI - Fast Analyzer" << rang::style::reset << std::endl
            << std::endl;
  gra::aux::PrintVersion();

  // Save the number of input arguments
  const int NARGC = argc - 1;

  try {
    cxxopts::Options options(argv[0], "");
    options.add_options()("i,input",
                          "input HepMC3 file                <input1,input2,...> (without .hepmc3)",
                          cxxopts::value<std::string>())(
        "g,pdg", "central final state PDG          <input1,input2,...>",
        cxxopts::value<std::string>())("n,number",
                                       "central final state multiplicity <input1,input2,...>",
                                       cxxopts::value<std::string>())(
        "l,labels", "plot legend string               <input1,input2,...>",
        cxxopts::value<std::string>())("t,title",
                                       "plot title string                <input>            ",
                                       cxxopts::value<std::string>())(
        "u,units", "plot unit                        <barn|mb|ub|nb|pb|fb>",
        cxxopts::value<std::string>())("M,mass", "plot mass binning                <bins,min,max>",
                                       cxxopts::value<std::string>())(
        "Y,rapidity", "plot rapidity binning            <bins,min,max>",
        cxxopts::value<std::string>())("P,momentum",
                                       "plot momentum binning            <bins,min,max>",
                                       cxxopts::value<std::string>())(
        "L,luminosity", "integrated luminosity (opt.)     <inverse barn>",
        cxxopts::value<double>())("X,maximum", "max nr. events to process (opt.) <value>",
                                  cxxopts::value<unsigned int>())(
        "S,scale", "scale plots                      <scale1,scale2,...>",
        cxxopts::value<std::string>())("R,ratio", "ratio plotting on                <true|false|1|0>",
                                       cxxopts::value<std::string>())("H,help", "Help");

    auto r = options.parse(argc, argv);

    if (r.count("help") || NARGC == 0) {
      std::cout << options.help({""}) << std::endl;
      std::cout << rang::style::bold << "Example:" << rang::style::reset << std::endl;
      std::cout << "  " << argv[0]
                << " -i ALICE_2pi,ALICE_2K -g 211,321 -n 2,2 -l "
                   "'#pi+#pi-','K+K-' -M 100,0.0,3.0 -Y 100,-1.5,1.5 -P 100,0.0,2.0 -u ub"
                << std::endl
                << std::endl;

      aux::CheckUpdate();
      return EXIT_FAILURE;
    }

    // Create Analysis Objects for each data input
    // NOTE HERE THAT THESE MUST BE POINTER TYPE; OTHERWISE WE RUN OUT OF
    // MEMORY
    std::vector<std::shared_ptr<MAnalyzer>> analysis;
    std::map<std::string, std::shared_ptr<h1Multiplet>>    h1;
    std::map<std::string, std::shared_ptr<h2Multiplet>>    h2;
    std::map<std::string, std::shared_ptr<hProfMultiplet>> hP;

    // Input list
    std::vector<std::string> inputfile     = gra::aux::SplitStr2Str(r["input"].as<std::string>());
    std::vector<std::string> labels        = gra::aux::SplitStr2Str(r["labels"].as<std::string>());
    std::vector<int>         finalstatePDG = gra::aux::SplitStr2Int(r["pdg"].as<std::string>());
    std::vector<int>         multiplicity  = gra::aux::SplitStr2Int(r["number"].as<std::string>());

    // Scaling
    std::vector<double> scale(inputfile.size(), 1.0);  // Default 1.0 for all
    if (r.count("scale")) {
      const std::vector<std::string> str_vals =
          gra::aux::SplitStr2Str(r["scale"].as<std::string>());
      if (str_vals.size() == inputfile.size()) {
        for (auto const &i : indices(str_vals)) {
          std::size_t consumed = 0;
          scale[i] = std::stod(str_vals[i], &consumed);
          if (consumed != str_vals[i].size() || !std::isfinite(scale[i])) {
            throw std::invalid_argument("analyze:: scale values must be finite numbers");
          }
        }
      } else {
        throw std::invalid_argument("analyzer::scale input list needs to be of length 0 or N");
      }
    }

    // Title string
    std::string title = "";
    if (r.count("title")) { title = r["title"].as<std::string>(); }

    unsigned int MAXEVENTS = std::numeric_limits<unsigned int>::max();
    if (r.count("maximum")) { MAXEVENTS = r["maximum"].as<unsigned int>(); }
    if (MAXEVENTS == 0) {
      throw std::invalid_argument("analyze:: maximum must be positive");
    }

    std::string units      = r["units"].as<std::string>();
    double      multiplier = 0.0;

    // Integrated luminosity given -> change units from dsigma/dx to dN/dx
    double luminosity = 1.0;
    if (r.count("luminosity")) {
      luminosity = r["luminosity"].as<double>();
      if (!std::isfinite(luminosity) || luminosity <= 0.0) {
        throw std::invalid_argument("analyze:: luminosity must be finite and positive");
      }
      gra::Scale(scale, luminosity);
      units = "events";
    } else if (units == "events") {
      throw std::invalid_argument(
          "analyze:: event yields require an integrated luminosity");
    }

    if (units == "events") {
      multiplier = 1.0;
    } else if (units == "barn") {
      multiplier = 1.0;
    } else if (units == "mb") {
      multiplier = 1E3;
    } else if (units == "ub") {
      units      = "#mub";
      multiplier = 1E6;
    } else if (units == "nb") {
      multiplier = 1E9;
    } else if (units == "pb") {
      multiplier = 1E12;
    } else if (units == "fb") {
      multiplier = 1E15;
    } else {
      throw std::invalid_argument("Unknown 'units' parameter: " + units);
    }

    // ---------------------------------------------------------------------
    // Check input lengths
    auto CheckInputLength = [](const std::vector<std::size_t> &x) {
      for (const auto &i : indices(x)) {
        for (const auto &j : indices(x)) {
          if (x[i] != x[j]) { return false; }
        }
      }
      return true;
    };
    const std::vector<std::size_t> lengths = {inputfile.size(), finalstatePDG.size(), labels.size(),
                                              multiplicity.size()};

    if (!CheckInputLength(lengths)) {
      throw std::invalid_argument(
          "Commandline input lengths (|inputfile| == |PDG| == |labels| == |multiplicity|) do not "
          "match!");
    }
    if (inputfile.empty()) {
      throw std::invalid_argument("analyze:: at least one input sample is required");
    }
    for (const auto &i : indices(multiplicity)) {
      if (multiplicity[i] <= 0) {
        throw std::invalid_argument("analyze:: multiplicities must be positive");
      }
      if (finalstatePDG[i] == 0) {
        throw std::invalid_argument("analyze:: central final-state PDG codes must be non-zero");
      }
    }

    // Scale each data source
    for (const auto &i : indices(scale)) {
      scale[i] *= multiplier;
      if (!std::isfinite(scale[i])) {
        throw std::invalid_argument("analyze:: combined scale is non-finite");
      }
    }

    // ---------------------------------------------------------------------
    // Create histogram and add pointer to the map

    auto tripletfunc = [&](const std::string &str) {
      const std::vector<std::string> tokens =
          gra::aux::SplitStr2Str(r[str].as<std::string>());
      if (tokens.size() != 3 || tokens[0].empty() || tokens[0][0] == '-') {
        throw std::invalid_argument(
            "analyze:: " + str + " discretization must be <bins,min,max>");
      }
      std::size_t bins_consumed = 0;
      std::size_t min_consumed  = 0;
      std::size_t max_consumed  = 0;
      const unsigned long long bins = std::stoull(tokens[0], &bins_consumed);
      const double min = std::stod(tokens[1], &min_consumed);
      const double max = std::stod(tokens[2], &max_consumed);
      if (bins_consumed != tokens[0].size() ||
          min_consumed != tokens[1].size() ||
          max_consumed != tokens[2].size() || bins == 0 ||
          bins > std::numeric_limits<unsigned int>::max() ||
          !std::isfinite(min) || !std::isfinite(max) || min >= max) {
        throw std::invalid_argument(
            "analyze:: invalid " + str + " discretization <bins,min,max>");
      }
      return h1Bound(static_cast<unsigned int>(bins), min, max);
    };

    h1Bound bM = tripletfunc("M");
    h1Bound bY = tripletfunc("Y");
    h1Bound bP = tripletfunc("P");

    Init1DHistogram(h1, labels, multiplicity, title, units, bM, bP, bY);
    Init2DHistogram(h2, labels, multiplicity, title, units, bM, bP, bY);
    gra::program::InitPrHistogram(hP, labels, multiplicity, title, bM);

    // Analyze them
    printf("Analyze:: \n");
    std::vector<double> normalization(inputfile.size(), 0.0);
    for (const auto &i : indices(inputfile)) {
      std::cout << i << " :: input:" << inputfile[i] << std::endl;

      analysis.push_back(std::make_shared<MAnalyzer>("ID" + std::to_string(i)));
      normalization[i] = analysis[i]->HepMC3_OracleFill(
          inputfile[i], static_cast<unsigned int>(multiplicity[i]), finalstatePDG[i],
          MAXEVENTS, h1, h2, hP, i);
    }
    // Name
    std::string fullpath = gra::aux::GetBasePath(2) + "/figs/";
    for (const auto &i : indices(inputfile)) {
      fullpath += inputfile[i];
      if (i < inputfile.size() - 1) { fullpath += "+"; }
    }
    fullpath += "/";  // important

    bool RATIOPLOT = true;
    if (r.count("ratio")) {
      RATIOPLOT = aux::ParseBool(r["ratio"].as<std::string>(), "--ratio");
    }

    // Iterate over all 1D-histograms
    for (const auto &x : h1) {
      x.second->NormalizeAll(normalization, scale);
      x.second->SaveFig(fullpath, RATIOPLOT);  // Select second member of map
    }

    // Iterate over all 2D-histograms
    for (const auto &x : h2) {
      x.second->NormalizeAll(normalization, scale);
      x.second->SaveFig(fullpath, RATIOPLOT);
    }

    // Iterate over all Profile-histograms
    for (const auto &x : hP) { x.second->SaveFig(fullpath, RATIOPLOT); }

    gra::rootstyle::MergePdfDirectory(fullpath, fullpath + "/merged.pdf", "h");

    // Plot all separate histograms
    for (const auto &i : indices(analysis)) { analysis[i]->PlotAll(title); }

  } catch (const std::invalid_argument &e) {
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
  } catch (const std::exception &e) {
    gra::aux::PrintGameOver();
    std::cerr << rang::fg::red << "Exception catched: std::exception: " << rang::fg::reset
              << e.what() << std::endl;
    return EXIT_FAILURE;
  } catch (...) {
    gra::aux::PrintGameOver();
    std::cerr << rang::fg::red << "Exception catched: non-standard exception"
              << rang::fg::reset << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "[analyze: done]" << std::endl;

  return EXIT_SUCCESS;
}
