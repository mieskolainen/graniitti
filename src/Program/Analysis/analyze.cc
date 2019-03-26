// GRANIITTI - Monte Carlo event generator for high energy diffraction
// https://github.com/mieskolainen/graniitti
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
#include "Graniitti/MAux.h"
#include "Graniitti/Analysis/MAnalyzer.h"
#include "Graniitti/Analysis/MMultiplet.h"

// Libraries
#include "cxxopts.hpp"
#include "rang.hpp"

using gra::aux::indices;

using namespace gra;

// Set "nice" 2D-plot style
// Read here more about problems with the Rainbow
void set_plot_style() {
    
    // Set Smooth color gradients
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;

    Double_t stops[NRGBs] = {0.00, 0.34, 0.61, 0.84, 1.00};
    Double_t red[NRGBs]   = {0.00, 0.00, 0.87, 1.00, 0.51};
    Double_t green[NRGBs] = {0.00, 0.81, 1.00, 0.20, 0.00};
    Double_t blue[NRGBs]  = {0.51, 1.00, 0.12, 0.00, 0.00};
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);

    // Black-Red palette
    gStyle->SetPalette(53); // 53/56 for inverted
    gStyle->SetTitleOffset(1.6, "x"); // X-axis title offset from axis
    gStyle->SetTitleOffset(1.0, "y"); // X-axis title offset from axis
    gStyle->SetTitleSize(0.03,  "x");  // X-axis title size
    gStyle->SetTitleSize(0.035, "y");
    gStyle->SetTitleSize(0.03,  "z");
    gStyle->SetLabelOffset(0.025);
}


// Global Style Setup
void setROOTstyle() {
    gStyle->SetOptStat(0); // Statistics BOX OFF [0,1]

    gStyle->SetOptFit(); // Fit parameters

    gStyle->SetTitleSize(0.0475, "t"); // Title with "t" (or anything else than xyz)
    gStyle->SetStatY(1.0);
    gStyle->SetStatX(1.0);
    gStyle->SetStatW(0.15);
    gStyle->SetStatH(0.09);
    
    // See below
    set_plot_style();
}


// Initialize 1D-histograms
void Init1DHistogram(std::map<std::string, std::unique_ptr<h1Multiplet>>& h,
                     const std::vector<std::string>& legendtext, std::vector<int> multiplicity, const std::string& units,
                     const h1Bound& bM, const h1Bound& bP, const h1Bound& bY) {

    std::string name = "null";

    // Central system observables
    name = "h1_S_M";
    h[name] = std::make_unique<h1Multiplet>(name, ";System M  (GeV);d#sigma/dM  (" + units + "/GeV)", bM.N, bM.min, bM.max, legendtext);

    name = "h1_S_Pt";
    h[name] = std::make_unique<h1Multiplet>(name, ";System P_{t} (GeV);d#sigma/dP_{t}  (" + units + "/GeV)", bP.N, bP.min, bP.max, legendtext);

    name = "h1_S_Y";
    h[name] = std::make_unique<h1Multiplet>(name, ";System Y;d#sigma/dY  (" + units + ")", bY.N, bY.min, bY.max, legendtext);

    // Central track observables
    name = "h1_1B_pt";
    h[name] = std::make_unique<h1Multiplet>(name, ";Final state p_{t} (GeV);d#sigma/dp_{t}  (" + units + "/GeV)", bP.N, bP.min, bP.max, legendtext);

    name = "h1_1B_eta";
    h[name] = std::make_unique<h1Multiplet>(name, ";Final state #eta;d#sigma/d#eta  (" + units + ")", bY.N, bY.min, bY.max, legendtext);


    // 2-Body observables
    if (std::find(multiplicity.begin(), multiplicity.end(), 2) != multiplicity.end()) {
    name = "h1_2B_acop";
    h[name] = std::make_unique<h1Multiplet>(
        name, ";Final state acoplanarity #rho = 1 - |#delta#phi|/#pi;d#sigma/d#rho  (" + units + "/rad)", 100, 0.0, 1.0, legendtext);
    
    name = "h1_2B_diffrap";
    h[name] = std::make_unique<h1Multiplet>(
        name, ";#deltay;d#sigma/d#deltay  (" + units + ")", bY.N, bY.min, bY.max, legendtext);
    
    }

    // 4-Body observables
    if (std::find(multiplicity.begin(), multiplicity.end(), 4) != multiplicity.end()) {
        // ...
    }

    // Forward proton observables
    name = "h1_PP_dphi";
    h[name] = std::make_unique<h1Multiplet>(name, ";Proton pair #delta#phi (rad);d#sigma/#delta#phi  ("  + units + "/rad)", 100, 0.0, 3.14159, legendtext);

    name = "h1_PP_t1";
    h[name] = std::make_unique<h1Multiplet>(name, ";Mandelstam -t_{1} (GeV^{2});d#sigma/dt  ("   + units + "/GeV^{2})", bP.N, bP.min, bP.max, legendtext);

    name = "h1_PP_dpt";
    h[name] = std::make_unique<h1Multiplet>(name, ";Proton pair |#delta#bar{p}_{t}| (GeV);d#sigma/|#delta#bar{p}_{t}|  (" + units + "/GeV)", bP.N, bP.min, bP.max, legendtext);
}


// Initialize 2D-histograms
void Init2DHistogram(std::map<std::string, std::unique_ptr<h2Multiplet>>& h,
                     const std::vector<std::string>& legendtext, std::vector<int> multiplicity, const std::string& units,
                     const h1Bound& bM, const h1Bound& bP, const h1Bound& bY) {

    std::string name  = "null";

    // Central system observables
    name = "h2_S_M_Pt";
    h[name] = std::make_unique<h2Multiplet>(name, "d#sigma^2/dMdP_{t}  (" + units + "/GeV/GeV)" + ";System M (GeV); System P_{t} (GeV)",
        bM.N, bM.min, bM.max, bP.N, bP.min, bP.max, legendtext);
    
    // 2-Body
    if (std::find(multiplicity.begin(), multiplicity.end(), 2) != multiplicity.end()) {
        name = "h2_2B_M_dphi";
        h[name] = std::make_unique<h2Multiplet>(name, "d#sigma^2/dMd#delta#phi  (" + units + "/GeV/rad)" + ";System M (GeV); Final state #delta#phi (rad)",
            bM.N, bM.min, bM.max, 100, 0.0, 3.14159, legendtext);

        name = "h2_2B_eta1_eta2";
        h[name] = std::make_unique<h2Multiplet>(name, "d#sigma^2/d#eta_{1}d#eta_{2}  (" + units + ")" + ";#eta_{1}; #eta_{2}",
            bY.N, bY.min, bY.max, bY.N, bY.min, bY.max, legendtext);
    }
    
    // 4-Body observables
    if (std::find(multiplicity.begin(), multiplicity.end(), 4) != multiplicity.end()) {
        // ...
    }
}


// Initialize Profile histograms
void InitPrHistogram(std::map<std::string, std::unique_ptr<hProfMultiplet>>& h,
    const std::vector<std::string>& legendtext, std::vector<int> multiplicity,
    const h1Bound& bM, const h1Bound& bP, const h1Bound& bY) {

    std::string name   = "null";

    // Central system observables
    name = "hP_S_M_Pt";
    h[name] = std::make_unique<hProfMultiplet>(
        name, ";System M  (GeV); System #LTP_{t}#GT  (GeV)",  bM.N, bM.min, bM.max, bP.min, bP.max, legendtext);

    // 2-body
    if (std::find(multiplicity.begin(), multiplicity.end(), 2) != multiplicity.end()) {
    name = "hP_2B_M_dphi";
    h[name] = std::make_unique<hProfMultiplet>(
        name, ";System M (GeV); Final state pair #LT#delta#phi#GT  (rad)",  bM.N, bM.min, bM.max, 0.0, 3.14159, legendtext);
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
int main(int argc, char* argv[]) {
    setROOTstyle();

    gra::aux::PrintFlashScreen(rang::fg::green);
    std::cout << rang::style::bold
              << "GRANIITTI - Fast Analyzer"
              << rang::style::reset << std::endl
              << std::endl;
    gra::aux::PrintVersion();

    // Save the number of input arguments
    const int NARGC = argc - 1;

try {

    cxxopts::Options options(argv[0], "");
    options.add_options()
        ("i,input",     "input HepMC3 file        <input1,input2,...> (without .hepmc3)", cxxopts::value<std::string>() )
        ("g,pdg",       "final state PDG          <input1,input2,...>",            cxxopts::value<std::string>() )
        ("n,number",    "final state multiplicity <input1,input2,...>",   cxxopts::value<std::string>() )
        ("l,labels",    "plot legend string       <input1,input2,...>",      cxxopts::value<std::string>() )
        ("u,units",     "cross section unit       <barn|mb|ub|nb|pb|fb>",    cxxopts::value<std::string>() )
        ("M,mass",      "plot mass limit",     cxxopts::value<double>() )
        ("Y,rapidity",  "plot rapidity limit", cxxopts::value<double>() )
        ("P,momentum",  "plot momentum limit", cxxopts::value<double>() )
        ("X,maximum",   "maximum number of events", cxxopts::value<int>() )
        ("S,scale",     "scale histograms", cxxopts::value<double>() )
        ("H,help",      "Help")
        ;
        
    auto r = options.parse(argc, argv);

    if (r.count("help") || NARGC == 0) {
        std::cout << options.help({""}) << std::endl;
        std::cout << "Example:" << std::endl;
        std::cout << "  " << argv[0] << " -i ALICE_2pi,ALICE_2K -g 211,321 -n 2,2 -l '#pi+#pi-','K+K-' -M 3.0 -Y 1.5 -P 2.0 -u ub"
                  << std::endl << std::endl;

        return EXIT_FAILURE;
    }
    

    // Create Analysis Objects for each data input
    // NOTE HERE THAT THESE MUST BE POINTER TYPE; OTHERWISE WE RUN OUT OF
    // MEMORY
    std::vector<MAnalyzer*> analysis;

    // Input list
    std::vector<std::string> inputfile = gra::aux::SplitStr2Str(r["input"].as<std::string>());
    std::vector<std::string> labels    = gra::aux::SplitStr2Str(r["labels"].as<std::string>());
    std::vector<int> finalstatePDG     = gra::aux::SplitStr2Int(r["pdg"].as<std::string>());
    std::vector<int> multiplicity      = gra::aux::SplitStr2Int(r["number"].as<std::string>());
    
    // Scaling
    double scale = 1.0;
    if (r.count("scale")) { scale = r["scale"].as<double>(); }

    const double MMAX  = r["mass"].as<double>();
    const double YMAX  = r["rapidity"].as<double>();
    const double PTMAX = r["momentum"].as<double>();

    int MAXEVENTS = 1e9;
    if (r.count("X")) { MAXEVENTS = r["X"].as<int>(); }

    std::string units = r["units"].as<std::string>();
    double multiplier = 0.0;

    if      (units == "barn") {
        multiplier = 1.0;
    }
    else if (units == "mb") {
        multiplier = 1E3;
    }
    else if (units == "ub") {
        units = "#mub";
        multiplier = 1E6;
    }
    else if (units == "nb") {
        multiplier = 1E9;       
    }
    else if (units == "pb") {
        multiplier = 1E12;
    }
    else if (units == "fb") {
        multiplier = 1E15;       
    } else {
        throw std::invalid_argument("Unknown 'units' parameter: " + units);
    }
    
    if (inputfile.size() != finalstatePDG.size() || inputfile.size() != multiplicity.size()) {
	   throw std::invalid_argument("Commandline input length do not match!");
    }
    
    // ---------------------------------------------------------------------
    // Create histogram and add pointer to the map

    h1Bound bM(95, 0, MMAX);
    h1Bound bP(95, 0, PTMAX);
    h1Bound bY(95, -YMAX, YMAX);

    Init1DHistogram(h1, labels, multiplicity, units, bM, bP, bY);
    Init2DHistogram(h2, labels, multiplicity, units, bM, bP, bY);
    InitPrHistogram(hP, labels, multiplicity, bM, bP, bY);

    // Analyze them
    printf("Analyze:: \n");
    std::vector<double> cross_section(inputfile.size(), 0.0);
    for (const auto& SID : indices(inputfile)) {
    	analysis.push_back(new MAnalyzer());
    	cross_section[SID] = analysis[SID]->HepMC3_OracleFill(inputfile[SID], (uint)multiplicity[SID],
    	                          finalstatePDG[SID], (uint)MAXEVENTS, h1, h2, hP, SID);
    }
    // Name
    std::string fullpath = gra::aux::GetBasePath(2) + "/figs/";
    for (const auto& i : indices(inputfile)) {
        fullpath += inputfile[i];
        if (i < inputfile.size()-1) { fullpath += "+"; }
    }
    fullpath += "/"; // important

    // Iterate over all 1D-histograms
    for (const auto& x : h1) {
    	x.second->NormalizeAll(cross_section, multiplier * scale);
    	x.second->SaveFig(fullpath); // Select second member of map
    }

    // Iterate over all 2D-histograms
    for (const auto& x : h2) {
    	x.second->NormalizeAll(cross_section, multiplier * scale);
    	x.second->SaveFig(fullpath);
    }

    // Iterate over all Profile-histograms
    for (const auto& x : hP) {
	   x.second->SaveFig(fullpath);
    }

    // Merge pdfs using Ghostscript (gs)
    const std::string cmd = "gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile="+fullpath+"/merged.pdf "+fullpath+"/h*.pdf";
    if (system(cmd.c_str()) == -1) {
        throw std::invalid_argument("Error: Problem executing Ghostscript merge on pdfs!");
    }

    // Print all separate histograms
    for (const auto& i : indices(analysis)) {
        analysis[i]->PlotAll();
    }

} catch (const std::invalid_argument& e) {
    gra::aux::PrintGameOver();
    std::cerr << rang::fg::red
              << "Exception catched: " << rang::fg::reset << e.what()
              << std::endl;
    return EXIT_FAILURE;
} catch (const std::ios_base::failure& e) {
    gra::aux::PrintGameOver();
    std::cerr << rang::fg::red
              << "Exception catched: std::ios_base::failure: "
              << rang::fg::reset << e.what() << std::endl;
    return EXIT_FAILURE;
} catch (const cxxopts::OptionException& e) {
    gra::aux::PrintGameOver();
    std::cerr << rang::fg::red
              << "Exception catched: Commandline options: " << rang::fg::reset << e.what() << std::endl;
    return EXIT_FAILURE;
} catch (...) {
    gra::aux::PrintGameOver();
    std::cerr << rang::fg::red << "Exception catched: Unspecified (...) (Probably input)" << rang::fg::reset
              << std::endl;
    return EXIT_FAILURE;
}

    return EXIT_SUCCESS;

}
