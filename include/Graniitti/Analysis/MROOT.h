// ROOT style namespace
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MROOT_H
#define MROOT_H

#include <complex>
#include <memory>
#include <tuple>
#include <vector>

// Own
#include "Graniitti/MMath.h"

// ROOT
#include "TCanvas.h"
#include "TColor.h"
#include "TLatex.h"
#include "TPad.h"
#include "TROOT.h"
#include "TStyle.h"

// Own
#include "Graniitti/MAux.h"

namespace gra {
namespace rootstyle {

// Returns a new colormap
inline void CreateColorMap(std::vector<int>& color, std::vector<std::shared_ptr<TColor>>& rootcolor, int COLORSCHEME = 1) {
  std::vector<std::vector<double>> colormap(150);

  if (COLORSCHEME == 1) {
    // "Modern colormap"
    std::vector<std::vector<double>> cm = {{0, 0.4470, 0.7410},      {0.8500, 0.3250, 0.0980},
                                           {0.9290, 0.6940, 0.1250}, {0.4940, 0.1840, 0.5560},
                                           {0.4660, 0.6740, 0.1880}, {0.3010, 0.7450, 0.9330},
                                           {0.6350, 0.0780, 0.1840}};
    colormap = cm;
  }

  if (COLORSCHEME == 2) {
    // "Classic colormap"
    std::vector<std::vector<double>> cm = {{0, 0, 0.9},       {0, 0.5, 0},     {0.9, 0, 0},
                                           {0, 0.75, 0.75},   {0.75, 0, 0.75}, {0.75, 0.75, 0},
                                           {0.25, 0.25, 0.25}};
    colormap = cm;
  }

  color     = std::vector<int>(colormap.size(), 0);
  rootcolor = std::vector<std::shared_ptr<TColor>> (color.size(), nullptr);

  for (const auto &i : aux::indices(color)) {
    color[i] = TColor::GetFreeColorIndex();
    //color[i] = 9000 + i;  // some big number not used

    // ROOT style, we need to create some hidden memory part
    rootcolor[i] = std::make_shared<TColor>(color[i], colormap[i][0], colormap[i][1], colormap[i][2]);
  }
}


// Create grid canvas with N subpads
inline void AutoGridCanvas(std::shared_ptr<TCanvas>& c1, unsigned int N) {
  
  unsigned int ADD = 0;
  while (true) {  // Adjust grid size
    const unsigned int val = std::sqrt(N + ADD);
    if (val * val == (N + ADD)) { break; }
    ++ADD;
  }
  // Calculate if we have a full empty row -> remove that
  unsigned int DEL = 0;
  if (ADD * ADD - ADD == N) { DEL = 1; }

  const unsigned int COLS = std::sqrt(N + ADD);
  const unsigned int ROWS = std::sqrt(N + ADD) - DEL;

  // Adjust aspect ratio
  if (ROWS == COLS) {
    c1 = std::make_shared<TCanvas>("c1", "c1", 700, 600);  // horizontal, vertical
  }
  if (ROWS != COLS) {
    c1 = std::make_shared<TCanvas>("c1", "c1", 700, 450);  // horizontal, vertical
  }

  c1->Divide(COLS, ROWS, 0.002, 0.001);

  // This is needed
  gStyle->SetPadLeftMargin(0.15);
}

// CubeHelix colormap generator
//
// N     = number of discrete steps
// start = start color (1 = red ... 2 = green ... 3 = red)
// R     = number of helix rotations
// hue   = hue, with 0 gives black&white
// gamma = intensity correction
// 
// Default CubeHelix(256, 0.5, -1.5, 1.2, 1.0);
// 
// [REFERENCE: D.A. Green, https://arxiv.org/abs/1108.5083]
// https://www.mrao.cam.ac.uk/~dag/CUBEHELIX
//
inline std::vector<std::vector<double>> CubeHelix(int N, double start, double R, double hue, double gamma) {

  auto limitfunc = [] (std::vector<double>& x) {
    for (std::size_t i = 0; i < x.size(); ++i) {
      if (x[i] < 0.0) { x[i] = 0.0; }
      if (x[i] > 1.0) { x[i] = 1.0; }
    }
  };

  // Color matrix
  const std::vector<std::vector<double>> A = {{-0.14861,  1.78277},
                                              {-0.29227, -0.90649},
                                              { 1.97294,  0}};
  const double PI = 3.14159265359;

  // Steps, Red, Green, Blue
  std::vector<std::vector<double>> M(4, std::vector<double>(N, 0.0));

  // Lightning
  const double maxlight = 1.0;
  const double minlight = 0.0;
  const double lightstep = (maxlight - minlight) / N;

  for (int i = 1; i <= N; ++i) {

    // Rotation angle and shift
    double alpha = (i - 1.0)/(N - 1.0);
    const double phi = 2.0*PI*(start/3.0 + 1.0 + R*alpha);

    // Apply gamma-correction
    alpha = std::pow(alpha, gamma);
    const double a = hue*alpha*(1.0 - alpha)/2.0;

    // Affine Map
    const std::vector<double> x = {a*std::cos(phi), a*std::sin(phi)};
    std::vector<double> y = {A[0][0]*x[0] + A[0][1]*x[1] + alpha,
                             A[1][0]*x[0] + A[1][1]*x[1] + alpha,
                             A[2][0]*x[0] + A[2][1]*x[1] + alpha};
    limitfunc(y); // Limit values to [0,1]

    // Save values
    M[0][i-1] = minlight + lightstep * (i-1);
    for (std::size_t j = 0; j < 3; ++j) { M[j+1][i-1] = y[j]; }
  }
  
  return M;
}


// Set "nice" 2D-plot style
inline void SetPlotStyle() {

  // Set smooth color gradients
  const int NCont = 256;

  const std::string style = "default";

  if      (style == "default") {
  const int NRGBs = 5;

  double stops[NRGBs] = {0.00, 0.34, 0.61, 0.84, 1.00};
  double red[NRGBs]   = {0.00, 0.00, 0.87, 1.00, 0.51};
  double green[NRGBs] = {0.00, 0.81, 1.00, 0.20, 0.00};
  double blue[NRGBs]  = {0.51, 1.00, 0.12, 0.00, 0.00};
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);

  // See https://root.cern.ch/doc/master/classTColor.html
  //gStyle->SetPalette(53);     // kDarkBodyRadiator
  //gStyle->SetPalette(51);       // kDeepSea
  //gStyle->SetPalette(87);     // kLightTemperature
  //gStyle->SetPalette(105);    // kThermometer
  //gStyle->SetPalette(71);     // kBlueGreenYellow
  //gStyle->SetPalette(57);     // kBird
  //gStyle->SetPalette(75);       // kCherry
  gStyle->SetPalette(112);      // kViridis

  //TColor::InvertPalette();      // Palette inversion
  }
  else if (style == "gray") {
  const int NRGBs = 5;

  double stops[NRGBs] = {0.00, 0.34, 0.61, 0.84, 1.00};
  double red[NRGBs]   = {1.00, 0.84, 0.61, 0.34, 0.00};
  double green[NRGBs] = {1.00, 0.84, 0.61, 0.34, 0.00};
  double blue[NRGBs]  = {1.00, 0.84, 0.61, 0.34, 0.00};
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  }
  else if (style == "cubehelix") {

  const int NRGBs = 256;
  const std::vector<std::vector<double>> M = CubeHelix(NRGBs, 0.5, -1.5, 1.2, 1.0);

  double stops[NRGBs];
  double red[NRGBs];
  double green[NRGBs];
  double blue[NRGBs];

  for (std::size_t i = 0; i < NRGBs; ++i) {
    stops[i] = M[0][i];
    red[i]   = M[1][i];
    green[i] = M[2][i];
    blue[i]  = M[3][i];
  }
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  }

  gStyle->SetNumberContours(NCont);

  gStyle->SetTitleOffset(1.6, "x");  // title offset from axis
  gStyle->SetTitleOffset(1.0, "y");  //
  gStyle->SetTitleSize(0.03, "x");   // title size
  gStyle->SetTitleSize(0.035, "y");
  gStyle->SetTitleSize(0.03, "z");
  gStyle->SetLabelOffset(0.025);

  // Necessary with multiple plots per canvas
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.15);

  gStyle->SetPadLeftMargin(0.1);
  gStyle->SetPadRightMargin(0.09);
}

// Global Style Setup
inline void SetROOTStyle() {
  gStyle->SetOptStat(0);  // Statistics BOX OFF [0,1]

  gStyle->SetOptFit();  // Fit parameters

  gStyle->SetTitleSize(0.04, "t");  // Title with "t" (or anything else than xyz)
  gStyle->SetStatY(1.0);
  gStyle->SetStatX(1.0);
  gStyle->SetStatW(0.15);
  gStyle->SetStatH(0.09);

  SetPlotStyle();
}

// Before calling this, call mother TCanvas cd->()
inline void TransparentPad(std::shared_ptr<TPad>& pad) {
  pad = std::make_shared<TPad>("newpad", "a transparent pad", 0, 0, 1, 1);
  pad->SetFillStyle(4000);
  pad->Draw();
  pad->cd();
}

// Create GRANIITTI Text
inline void MadeInFinland(std::shared_ptr<TLatex>& l1, std::shared_ptr<TLatex>& l2, double xpos = 0.935, std::vector<double> ypos = {0.03, 0.58}) {
    
  if (ypos.size() != 2) {
    throw std::invalid_argument("MROOT::MadeInFinland: argument ypos should be of size 2");
  }

  l1 = std::make_shared<TLatex>(xpos, ypos[0], gra::aux::GetVersionTLatex().c_str());
  l1->SetNDC();  // Normalized coordinates
  l1->SetTextAngle(90);
  l1->Draw();

  l2 = std::make_shared<TLatex>(xpos, ypos[1], gra::aux::GetWebTLatex().c_str());
  l2->SetNDC();  // Normalized coordinates
  l2->SetTextAngle(90);
  l2->Draw();
}

}  // rootstyle namespace ends
}  // gra namespace ends

#endif
