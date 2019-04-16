// ROOT style namespace
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MROOT_H
#define MROOT_H

#include <complex>
#include <memory>
#include <vector>
#include <tuple>

// ROOT
#include "TCanvas.h"
#include "TColor.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TPad.h"


namespace gra {
namespace rootstyle {


// Create grid canvas with N subpads
inline TCanvas* AutoGridCanvas(unsigned int N) {
    
    TCanvas* c1;
    
    unsigned int ADD = 0;
    while (true) { // Adjust grid size
        const unsigned int val = std::sqrt(N + ADD); 
        if ( val*val == (N + ADD) ) { break; }
        ++ADD;
    }
    // Calculate if we have a full empty row -> remove that
    unsigned int DEL = 0;
    if (ADD*ADD - ADD == N) {DEL = 1;}
    
    const unsigned int COLS = std::sqrt(N+ADD);
    const unsigned int ROWS = std::sqrt(N+ADD)-DEL;
    
    // Adjust aspect ratio
    if (ROWS == COLS) {
    c1 = new TCanvas("c1", "c1", 700, 600); // horizontal, vertical
    }
    if (ROWS != COLS) {
    c1 = new TCanvas("c1", "c1", 700, 450); // horizontal, vertical
    }
    
    c1->Divide(COLS, ROWS, 0.002, 0.001);
    
    // This is needed
    gStyle->SetPadLeftMargin(0.15);

    return c1;
}

// Set "nice" 2D-plot style
inline void SetPlotStyle() {
    
    // Set Smooth color gradients
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;

    double stops[NRGBs] = {0.00, 0.34, 0.61, 0.84, 1.00};
    double red[NRGBs]   = {0.00, 0.00, 0.87, 1.00, 0.51};
    double green[NRGBs] = {0.00, 0.81, 1.00, 0.20, 0.00};
    double blue[NRGBs]  = {0.51, 1.00, 0.12, 0.00, 0.00};
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
    
    // Black-Red palette
    gStyle->SetPalette(53); // 53/56 for inverted
    gStyle->SetTitleOffset(1.6, "x");  // title offset from axis
    gStyle->SetTitleOffset(1.0, "y");  // 
    gStyle->SetTitleSize(0.03,  "x");  // title size
    gStyle->SetTitleSize(0.035, "y");
    gStyle->SetTitleSize(0.03,  "z");
    gStyle->SetLabelOffset(0.025);
    
    // Necessary with multiple plots per canvas
    gStyle->SetPadTopMargin(0.1);
    gStyle->SetPadBottomMargin(0.15);

    gStyle->SetPadLeftMargin(0.1);
    gStyle->SetPadRightMargin(0.09);
}

// Global Style Setup
inline void SetROOTStyle() {
    gStyle->SetOptStat(0); // Statistics BOX OFF [0,1]

    gStyle->SetOptFit(); // Fit parameters

    gStyle->SetTitleSize(0.04, "t"); // Title with "t" (or anything else than xyz)
    gStyle->SetStatY(1.0);
    gStyle->SetStatX(1.0);
    gStyle->SetStatW(0.15);
    gStyle->SetStatH(0.09);
    
    SetPlotStyle();
}

// Before calling this, call mother TCanvas cd->()
inline TPad* TransparentPad() {

    TPad* newpad = new TPad("newpad", "a transparent pad", 0,0,1,1);
    newpad->SetFillStyle(4000);
    newpad->Draw();
    newpad->cd();
    return newpad;
}

// Create GRANIITTI Text
inline std::tuple<TLatex*,TLatex*> MadeInFinland(double xpos = 0.93) {
    
    TLatex* l1 = new TLatex(xpos, 0.03, gra::aux::GetVersionTLatex().c_str());
    l1->SetNDC(); // Normalized coordinates
    l1->SetTextAngle(90);
    l1->Draw();

    TLatex* l2 = new TLatex(xpos, 0.68, gra::aux::GetWebTLatex().c_str());
    l2->SetNDC(); // Normalized coordinates
    l2->SetTextAngle(90);
    l2->Draw();

    return std::make_tuple(l1, l2);
}

} // rootstyle namespace ends
} // gra namespace ends

#endif