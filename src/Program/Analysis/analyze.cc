// GRANIITTI - Monte Carlo event generator for high energy diffraction
// https://github.com/mieskolainen/graniitti
//
// <Fiducial Observables Analyzer>
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <math.h>
#include <algorithm>
#include <complex>
#include <cstdlib>
#include <iomanip>
#include <iostream>
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
#include "Graniitti/Analysis/MMultiplet.h"
#include "Graniitti/Analysis/MROOT.h"
#include "Graniitti/MAux.h"
#include "Graniitti/MMath.h"

// Libraries
#include "json.hpp"
#include "cxxopts.hpp"
#include "rang.hpp"

using gra::aux::indices;
using namespace gra;

// Initialize 1D-histograms
void Init1DHistogram(std::map<std::string, std::unique_ptr<h1Multiplet>> &h,
                     const std::vector<std::string> &legendtext, std::vector<int> multiplicity,
                     const std::string &title, const std::string &units, const h1Bound &bM,
                     const h1Bound &bP, const h1Bound &bY) {
  std::string       name = "null";
  const std::string U    = (units != "events") ? "#sigma" : "N";

  // Central system observables
  name    = "h1_S_M";
  h[name] = std::make_unique<h1Multiplet>(
      name, title + ";System M  (GeV);d" + U + "/dM  (" + units + "/GeV)", bM.N, bM.min, bM.max,
      legendtext);

  name    = "h1_S_Pt";
  h[name] = std::make_unique<h1Multiplet>(
      name, title + ";System P_{T} (GeV);d" + U + "/dP_{T}  (" + units + "/GeV)", bP.N, bP.min,
      bP.max, legendtext);

  name    = "h1_S_Y";
  h[name] = std::make_unique<h1Multiplet>(name, title + ";System Y;d" + U + "/dY  (" + units + ")",
                                          bY.N, bY.min, bY.max, legendtext);

  // Central track observables
  name    = "h1_1B_pt";
  h[name] = std::make_unique<h1Multiplet>(
      name, title + ";Central final state p_{T} (GeV);d" + U + "/dp_{T}  (" + units + "/GeV)", bP.N,
      bP.min, bP.max, legendtext);

  name    = "h1_1B_eta";
  h[name] = std::make_unique<h1Multiplet>(
      name, title + ";Central final state #eta;d" + U + "/d#eta  (" + units + ")", bY.N, bY.min,
      bY.max, legendtext);

  // 2-Body observables
  if (std::find(multiplicity.begin(), multiplicity.end(), 2) != multiplicity.end()) {

    for (const auto& i : analyzer::FRAMES) {
      name    = "h1_costheta_" + i;
      h[name] = std::make_unique<h1Multiplet>(
        name, title + ";Central final state cos(#theta) [" + i + " frame];d" + U + "/dcos(#theta)  (" + units + ")", bP.N, -1.0, 1.0, legendtext);

      name    = "h1_phi_" + i;
      h[name] = std::make_unique<h1Multiplet>(
        name, title + ";Central final state #phi [" + i + " frame];d" + U + "/dcos #phi  (" + units + ")", bP.N, -math::PI, math::PI, legendtext);
    }
    
    name = "h1_2B_acop";
    h[name] =
        std::make_unique<h1Multiplet>(name, title +
                                                ";Central final state acoplanarity #rho = 1 - "
                                                "|#Delta#phi|/#pi;d" +
                                                U + "/d#rho  (" + units + "/rad)",
                                      100, 0.0, 1.0, legendtext);

    name    = "h1_2B_diffrap";
    h[name] = std::make_unique<h1Multiplet>(
        name, title + ";#deltay #equiv y_{1} - y_{2};d" + U + "/d#Deltay  (" + units + ")", bY.N,
        bY.min, bY.max, legendtext);
  }

  // 4-Body observables
  if (std::find(multiplicity.begin(), multiplicity.end(), 4) != multiplicity.end()) {
    // ...
  }

  // Forward proton observables
  name    = "h1_PP_dphi";
  h[name] = std::make_unique<h1Multiplet>(
      name, title + ";Proton pair #Delta#phi (rad);d" + U + "/d#Delta#phi  (" + units + "/rad)", 100,
      0.0, 3.14159, legendtext);

  name    = "h1_PP_t1";
  h[name] = std::make_unique<h1Multiplet>(
      name, title + ";Mandelstam -t_{1} (GeV^{2});d" + U + "/dt  (" + units + "/GeV^{2})", bP.N,
      bP.min, bP.max, legendtext);

  name    = "h1_PP_dpt";
  h[name] = std::make_unique<h1Multiplet>(name, title +
                                                    ";Proton pair |#Delta#bar{p}_{T}| "
                                                    "(GeV);d" +
                                                    U + "/d|#Delta#bar{p}_{T}|  (" + units + "/GeV)",
                                          bP.N, bP.min, bP.max, legendtext);
}

// Initialize 2D-histograms
void Init2DHistogram(std::map<std::string, std::unique_ptr<h2Multiplet>> &h,
                     const std::vector<std::string> &legendtext, std::vector<int> multiplicity,
                     const std::string &title, const std::string &units, const h1Bound &bM,
                     const h1Bound &bP, const h1Bound &bY) {
  std::string       name = "null";
  const std::string U    = (units != "events") ? "#sigma" : "N";

  // Central system observables
  name = "h2_S_M_Pt";
  h[name] =
      std::make_unique<h2Multiplet>(name, "d" + U + "^2/dMdP_{T}  (" + units + "/GeV/GeV) | " +
                                              title + ";System M (GeV); System P_{T} (GeV)",
                                    bM.N, bM.min, bM.max, bP.N, bP.min, bP.max, legendtext);

  name    = "h2_S_M_t";
  h[name] = std::make_unique<h2Multiplet>(name, "d" + U + "^2/dMd#t  (" + units +
                                                    "/GeV/GeV) | " + title +
                                                    ";System M (GeV); |t| (GeV^{2})",
                                          bM.N, bM.min, bM.max, bP.N, math::pow2(bP.min), math::pow2(bP.max), legendtext);

  name    = "h2_S_M_pt";
  h[name] = std::make_unique<h2Multiplet>(
      name, "d" + U + "^2/dMdp_{T}  (" + units + "/GeV/GeV) | " + title +
                ";System M (GeV); Central final state p_{T} (GeV)",
      bM.N, bM.min, bM.max, bP.N, bP.min, bP.max, legendtext);

  name    = "h2_S_M_dphipp";
  h[name] = std::make_unique<h2Multiplet>(
      name, "d" + U + "^2/dMd#Delta#phi_{pp}  (" + units + "/GeV/rad) | " + title +
                ";System M (GeV); Forward proton #Delta#phi_{pp}",
      bM.N, bM.min, bM.max, 100, 0.0, gra::math::PI, legendtext);

  name    = "h2_S_M_dpt";
  h[name] = std::make_unique<h2Multiplet>(
      name, "d" + U + "^2/dMd|#Delta#bar{p}_{T}|  (" + units + "/GeV/GeV) | " + title +
                ";System M (GeV); Proton pair |#Delta#bar{p}_{T}| (GeV)",
      bM.N, bM.min, bM.max, 100, 0.0, 2.0, legendtext);

  // 2-Body
  if (std::find(multiplicity.begin(), multiplicity.end(), 2) != multiplicity.end()) {
    name    = "h2_2B_M_dphi";
    h[name] = std::make_unique<h2Multiplet>(
        name, "d" + U + "^2/dMd#Delta#phi  (" + units + "/GeV/rad) | " + title +
                  ";System M (GeV); Central final state #Delta#phi (rad)",
        bM.N, bM.min, bM.max, 100, 0.0, gra::math::PI, legendtext);

    name    = "h2_2B_eta1_eta2";
    h[name] = std::make_unique<h2Multiplet>(
        name, "d" + U + "^2/d#eta_{1}d#eta_{2}  (" + units + ") | " + title + ";#eta_{1}; #eta_{2}",
        bY.N, bY.min, bY.max, bY.N, bY.min, bY.max, legendtext);
    
    // 2D (costheta, phi) in different rest frames
    for (const auto& i : analyzer::FRAMES) {

      name = "h2_2B_costheta_phi_" + i;
      h[name] = std::make_unique<h2Multiplet>(
          name, "d" + U + "^2/cos(#theta)#phi  (" + units + ") | " + title + ";daughter cos(#theta); daughter #phi (rad) [" + i + " FRAME]",
          100, -1, 1, 100, -gra::math::PI, gra::math::PI, legendtext);
    }
    
    // 2D (M, costheta) in different rest frames
    for (const auto& i : analyzer::FRAMES) {
      
      name = "h2_2B_M_costheta_" + i;
      h[name] = std::make_unique<h2Multiplet>(
          name, "d" + U + "^2/dMdcos(#theta)  (" + units + "/GeV) | " + title + ";M (GeV); daughter cos(#theta) [" + i + " FRAME]",
          bM.N, bM.min, bM.max, 100, -1, 1, legendtext);
    }
    
    // 2D (M, phi) in different rest frames
    for (const auto& i : analyzer::FRAMES) {
      
      name = "h2_2B_M_phi_" + i;
      h[name] = std::make_unique<h2Multiplet>(
          name, "d" + U + "^2/dMd#phi  (" + units + "/GeV/rad) | " + title + ";M (GeV); daughter #phi (rad) [" + i + " FRAME]",
          bM.N, bM.min, bM.max, 100, -gra::math::PI, gra::math::PI, legendtext);
    }
  }

  // 4-Body observables
  if (std::find(multiplicity.begin(), multiplicity.end(), 4) != multiplicity.end()) {
    // ...
  }
}

// Initialize Profile histograms
void InitPrHistogram(std::map<std::string, std::unique_ptr<hProfMultiplet>> &h,
                     const std::vector<std::string> &legendtext, std::vector<int> multiplicity,
                     const std::string &title, const h1Bound &bM, const h1Bound &bP,
                     const h1Bound &bY) {
  std::string name = "null";

  // Central system observables
  name = "hP_S_M_Pt";
  h[name] =
      std::make_unique<hProfMultiplet>(name, title + ";System M  (GeV); System #LTP_{T}#GT (GeV)",
                                       bM.N, bM.min, bM.max, bP.min, bP.max, legendtext);

  name    = "hP_S_M_PL2_SR";
  h[name] = std::make_unique<hProfMultiplet>(
      name, title + ";System M  (GeV); Legendre #LTP_{l=2}(cos #theta)#GT [SR frame]", bM.N, bM.min,
      bM.max, -0.4, 0.4, legendtext);

  name    = "hP_S_M_PL4_SR";
  h[name] = std::make_unique<hProfMultiplet>(
      name, title + ";System M  (GeV); Legendre #LTP_{l=4}(cos #theta)#GT [SR frame]", bM.N, bM.min,
      bM.max, -0.4, 0.4, legendtext);

  // 2-body
  if (std::find(multiplicity.begin(), multiplicity.end(), 2) != multiplicity.end()) {
    name    = "hP_2B_M_dphi";
    h[name] = std::make_unique<hProfMultiplet>(
        name, title + ";System M (GeV); Central final state pair #LT#delta#phi#GT  (rad)", bM.N,
        bM.min, bM.max, 0.0, math::PI, legendtext);
  }

  // 4-Body observables
  if (std::find(multiplicity.begin(), multiplicity.end(), 4) != multiplicity.end()) {
    // ...
  }
}

// Histogram collected here
std::map<std::string, std::unique_ptr<h1Multiplet>>    h1;
std::map<std::string, std::unique_ptr<h2Multiplet>>    h2;
std::map<std::string, std::unique_ptr<hProfMultiplet>> hP;

// -----------------------------------------------------------------------

// Main program
int main(int argc, char *argv[]) {
  gra::rootstyle::SetROOTStyle();

  gra::aux::PrintFlashScreen(rang::fg::green);
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
        cxxopts::value<double>())("X,maximum", "max number of events to process (opt.)  <value>",
                                  cxxopts::value<int>())(
        "S,scale", "scale plots                      <scale1,scale2,...>",
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
    std::vector<std::unique_ptr<MAnalyzer>> analysis;

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
        for (auto const &i : indices(str_vals)) { scale[i] = std::stod(str_vals[i]); }
      } else {
        throw std::invalid_argument("analyzer::scale input list needs to be of length 0 or N");
      }
    }

    // Title string
    std::string title = "";
    if (r.count("title")) { title = r["title"].as<std::string>(); }

    int MAXEVENTS = 1e9;
    if (r.count("X")) { MAXEVENTS = r["X"].as<int>(); }

    std::string units      = r["units"].as<std::string>();
    double      multiplier = 0.0;

    // Integrated luminosity given -> change units from dsigma/dx to dN/dx
    double luminosity = 1.0;
    if (r.count("luminosity")) {
      luminosity = r["luminosity"].as<double>();
      for (auto const &i : indices(scale)) { scale[i] *= luminosity; }
      units = "events";
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

    if (inputfile.size() != finalstatePDG.size() || inputfile.size() != multiplicity.size()) {
      throw std::invalid_argument("Commandline input length do not match!");
    }

    // Scale each data source
    for (const auto &i : indices(scale)) { scale[i] *= multiplier; }

    // ---------------------------------------------------------------------
    // Create histogram and add pointer to the map

    auto tripletfunc = [&](const std::string &str) {

      const std::string         vecstr = r[str].as<std::string>();
      const std::vector<double> vec    = gra::aux::SplitStr(vecstr, double(0), ',');
      if (vec.size() != 3) {
        throw std::invalid_argument("analyze:: " + str +
                                    " discretization not size 3 <bins,min,max>");
      }
      return h1Bound(vec[0], vec[1], vec[2]);
    };

    h1Bound bM = tripletfunc("M");
    h1Bound bY = tripletfunc("Y");
    h1Bound bP = tripletfunc("P");

    Init1DHistogram(h1, labels, multiplicity, title, units, bM, bP, bY);
    Init2DHistogram(h2, labels, multiplicity, title, units, bM, bP, bY);
    InitPrHistogram(hP, labels, multiplicity, title, bM, bP, bY);

    // Analyze them
    printf("Analyze:: \n");
    std::vector<double> cross_section(inputfile.size(), 0.0);
    for (const auto &i : indices(inputfile)) {

      std::cout << i << " :: input:" << inputfile[i] << std::endl;

      analysis.push_back(std::make_unique<MAnalyzer>("ID" + std::to_string(i)));
      cross_section[i] = analysis[i]->HepMC3_OracleFill(
          inputfile[i], (uint)multiplicity[i], finalstatePDG[i], (uint)MAXEVENTS, h1, h2, hP, i);
    }
    // Name
    std::string fullpath = gra::aux::GetBasePath(2) + "/figs/";
    for (const auto &i : indices(inputfile)) {
      fullpath += inputfile[i];
      if (i < inputfile.size() - 1) { fullpath += "+"; }
    }
    fullpath += "/";  // important

    // Iterate over all 1D-histograms
    for (const auto &x : h1) {
      x.second->NormalizeAll(cross_section, scale);
      x.second->SaveFig(fullpath);  // Select second member of map
    }

    // Iterate over all 2D-histograms
    for (const auto &x : h2) {
      x.second->NormalizeAll(cross_section, scale);
      x.second->SaveFig(fullpath);
    }

    // Iterate over all Profile-histograms
    for (const auto &x : hP) { x.second->SaveFig(fullpath); }

    // Merge pdfs using Ghostscript (gs)
    const std::string cmd = "gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=" + fullpath +
                            "/merged.pdf " + fullpath + "/h*.pdf";
    if (system(cmd.c_str()) == -1) {
      throw std::invalid_argument("Error: Problem executing Ghostscript merge on pdfs!");
    }

    // Plot all separate histograms
    for (const auto &i : indices(analysis)) {
      analysis[i]->PlotAll(title);
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
  } catch (const cxxopts::OptionException &e) {
    gra::aux::PrintGameOver();
    std::cerr << rang::fg::red << "Exception catched: Commandline options: " << rang::fg::reset
              << e.what() << std::endl;
    return EXIT_FAILURE;
  } catch (const nlohmann::json::exception &e) {
    gra::aux::PrintGameOver();
    std::cerr << rang::fg::red << "Exception catched: JSON input: " << rang::fg::reset
              << e.what() << std::endl;
    return EXIT_FAILURE;
  } catch (...) {
    gra::aux::PrintGameOver();
    std::cerr << rang::fg::red << "Exception catched: Unspecified (...) (Probably input)"
              << rang::fg::reset << std::endl;
    return EXIT_FAILURE;
  }
  
  std::cout << "[analyze: done]" << std::endl;
  aux::CheckUpdate();

  return EXIT_SUCCESS;
}
