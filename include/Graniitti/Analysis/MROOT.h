// ROOT style namespace
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MROOT_H
#define MROOT_H

#include <complex>
#include <memory>
#include <vector>

// ROOT
#include "TColor.h"
#include "TROOT.h"
#include "TStyle.h"


namespace gra {

namespace rootstyle {

// Set "nice" 2D-plot style
void SetPlotStyle() {
    
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
    gStyle->SetTitleOffset(1.6, "x");  // X-axis title offset from axis
    gStyle->SetTitleOffset(1.0, "y");  // X-axis title offset from axis
    gStyle->SetTitleSize(0.03,  "x");  // X-axis title size
    gStyle->SetTitleSize(0.035, "y");
    gStyle->SetTitleSize(0.03,  "z");
    gStyle->SetLabelOffset(0.025);

    // Necessary with multiple plots per canvas
    gStyle->SetPadTopMargin(0.1);
    gStyle->SetPadRightMargin(0.1);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadLeftMargin(0.1);
}


// Global Style Setup
void SetROOTStyle() {
    gStyle->SetOptStat(0); // Statistics BOX OFF [0,1]

    gStyle->SetOptFit(); // Fit parameters

    gStyle->SetTitleSize(0.04, "t"); // Title with "t" (or anything else than xyz)
    gStyle->SetStatY(1.0);
    gStyle->SetStatX(1.0);
    gStyle->SetStatW(0.15);
    gStyle->SetStatH(0.09);
    
    SetPlotStyle();
}

} // rootstyle namespace ends
} // gra namespace ends

#endif