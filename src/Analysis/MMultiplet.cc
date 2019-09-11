// (Histo0, Histo1, Histo2, ..., Histo N-1) Multiplet ROOT histograms
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <string>
#include <tuple>

// ROOT
#include "TColor.h"
#include "TGaxis.h"
#include "TH1.h"
#include "TH2.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TProfile.h"

// Own
#include "Graniitti/Analysis/MMultiplet.h"
#include "Graniitti/Analysis/MROOT.h"
#include "Graniitti/MAux.h"

using gra::aux::indices;

namespace gra {

void GetLegendPosition(unsigned int N, double &x1, double &x2, double &y1, double &y2,
                       const std::string &legendposition) {
  // North-East
  x1 = 0.70;
  x2 = x1 + 0.18;
  y1 = 0.75 - 0.02 * N;

  // South-East
  if (legendposition.compare("southeast") == 0) { y1 = 0.10 - 0.01 * N; }
  y2 = y1 + 0.04 * N;  // Scale by the number of histograms
}

h1Multiplet::h1Multiplet(const std::string &name, const std::string &labeltext, int N,
                         double minval, double maxval, const std::vector<std::string> &legendtext) {
  N_              = N;
  name_           = name;
  minval_         = minval;
  maxval_         = maxval;
  legendtext_     = legendtext;
  legendposition_ = "northeast";

  // Initialize histogram vector container size
  h.resize(legendtext.size());

  for (const auto &i : indices(h)) {
    h[i] = new TH1D(Form("%s_%lu", name.c_str(), i), labeltext.c_str(), N, minval, maxval);
    h[i]->Sumw2();  // Error saving on
  }
}

// Histogram normalization
void h1Multiplet::NormalizeAll(const std::vector<double> &cross_section,
                               const std::vector<double> &multiplier) {
  for (const auto &i : indices(h)) {
    double scale = 1.0;

    // Takes into account weighted and unweight (weight=1) filling
    // We take into account also under/overflow here with macro: (0, Nbins+1)
    const double integral = h[i]->Integral(0, h[i]->GetNbinsX() + 1);
    if (integral > 0) { scale /= integral; }

    // Binwidth
    const int    bin      = 1;
    const double binwidth = h[i]->GetXaxis()->GetBinWidth(bin);
    if (binwidth > 0) { scale /= binwidth; }

    // Cross Section
    if (cross_section[i] > 0) { scale *= cross_section[i]; }

    // Additional scale factor
    scale *= multiplier[i];

    h[i]->Scale(scale);
  }
}

std::vector<double> h1Multiplet::SaveFig(const std::string &fullpath) const {
  
  // New colorscheme
  std::vector<int> color;
  std::vector<std::shared_ptr<TColor>> rootcolor;
  rootstyle::CreateColorMap(color, rootcolor);

  // ----------------------------------------------------
  // Apply the chi2 test and retrieve the residuals
  gra::aux::PrintBar("*");
  std::vector<double> chi2ndf(h.size(), 0.0);

  for (const auto &i : indices(h)) {
    double res[N_];
    printf("%s [%lu] :: \n", legendtext_[i].c_str(), i);
    double c2ndf = h[0]->Chi2Test(h[i], "WW P CHI2/NDF", res);
    chi2ndf[i]   = c2ndf;

    printf("chi2/ndf = %0.3f \n\n", c2ndf);
  }
  gra::aux::PrintBar("*");

  TCanvas c0("c", "c", 750, 800);

  // Upper plot will be in pad1
  std::shared_ptr<TPad> pad1 = std::make_shared<TPad>("pad1", "pad1", 0, 0.3, 1, 1.0);

  pad1->SetBottomMargin(0.015);  // Upper and lower plot are joined
  // pad1->SetGridx();          // Vertical grid
  pad1->Draw();       // Draw the upper pad: pad1
  pad1->cd();         // pad1 becomes the current pad
  h[0]->SetStats(0);  // No statistics on upper plot

  // Find maximum value for y-range limits
  double MAXVAL = 0.0;
  for (const auto &i : indices(h)) {
    const double max = h[i]->GetMaximum();
    if (max > MAXVAL) { MAXVAL = max; }
  }

  // Find minimum (non-zero) value for y-range limits
  double MINVAL = 1e32;
  for (const auto &i : indices(h)) {
    for (int k = 0; k < h[i]->GetNbinsX(); ++k) {
      double value = h[i]->GetBinContent(k);
      if (value < MINVAL && value > 0) { MINVAL = value; }
    }
  }

  // Loop over histograms
  for (const auto &i : indices(h)) {
    h[i]->SetLineColor(color[i]);
    h[i]->SetLineWidth(2);
    h[i]->SetMarkerColor(color[i]);
    h[i]->SetMarkerStyle(20 + i);
    h[i]->SetMarkerSize(0.73);

    std::cout << MINVAL << std::endl;
    h[i]->GetYaxis()->SetRangeUser(0, MAXVAL * 1.5);

    // h[i]->Draw("e2 same");   // Needs this combo to draw box-lines
    h[i]->Draw("hist same");

    // h[i]->SetFillStyle(4050);
    // h[i]->SetFillColor(color[i]);
    // h[i]->SetFillColorAlpha(color[i], 0.1);
    // h[i]->Draw("f hist same"); // Filled
  }

  pad1->RedrawAxis();  // Fix overlapping histogram lines on borders

  // -------------------------------------------------------------------
  // Create dynamically sized legend (depending on number of legend entries)

  double x1, x2, y1, y2 = 0.0;
  GetLegendPosition(h.size(), x1, x2, y1, y2, legendposition_);
  std::shared_ptr<TLegend> legend = std::make_shared<TLegend>(x1, y1, x2, y2);
  legend->SetFillColor(0);  // White background
  // legend->SetBorderSize(0); // No box

  // Add legend entries
  for (const auto &i : indices(h)) { legend->AddEntry(h[i], legendtext_[i].c_str()); }
  legend->Draw();

  // -------------------------------------------------------------------
  // Create labels
  std::shared_ptr<TLatex> l1;
  std::shared_ptr<TLatex> l2;
  gra::rootstyle::MadeInFinland(l1, l2);

  // -------------------------------------------------------------------
  // Ratio plots
  c0.cd();
  std::shared_ptr<TPad> pad2 = std::make_shared<TPad>("pad2", "pad2", 0, 0.05, 1, 0.3);
  pad2->SetTopMargin(0.025);
  pad2->SetBottomMargin(0.25);
  pad2->SetGridx();  // vertical grid
  pad2->Draw();
  pad2->cd();  // pad2 becomes the current pad

  // *** Ratio histograms ***
  std::vector<TH1D *> hR(h.size(), nullptr);

  for (const auto &i : indices(h)) {
    hR[i] = static_cast<TH1D*>( h[i]->Clone(Form("ratio_%lu", i)) );
    hR[i]->Divide(h[0]);

    hR[i]->SetMinimum(0.0);  // y-axis range
    hR[i]->SetMaximum(2.0);  //
    hR[i]->SetStats(0);      // statistics box off
    hR[i]->Draw("same");     // ratio plot

    // Ratio plot (h3) settings
    hR[i]->SetTitle("");  // Remove the ratio title

    // Y axis ratio plot settings
    hR[i]->GetYaxis()->SetTitle("Ratio");
    hR[i]->GetYaxis()->CenterTitle();
    hR[i]->GetYaxis()->SetNdivisions(505);
    hR[i]->GetYaxis()->SetTitleSize(22);
    hR[i]->GetYaxis()->SetTitleFont(43);
    hR[i]->GetYaxis()->SetTitleOffset(1.55);
    hR[i]->GetYaxis()->SetLabelFont(45);  // Absolute font size in pixel (precision 3)
    hR[i]->GetYaxis()->SetLabelSize(15);

    // X axis ratio plot settings
    hR[i]->GetXaxis()->SetTitleSize(22);
    hR[i]->GetXaxis()->SetTitleFont(43);
    hR[i]->GetXaxis()->SetTitleOffset(4.);
    hR[i]->GetXaxis()->SetLabelFont(45);  // Absolute font size in pixel (precision 3)
    hR[i]->GetXaxis()->SetLabelSize(18);
  }

  // Draw horizontal line
  const double           ymax = 1.0;
  std::shared_ptr<TLine> line = std::make_shared<TLine>(minval_, ymax, maxval_, ymax);
  line->SetLineColor(15);
  line->SetLineWidth(2.0);
  line->Draw();

  // -------------------------------------------------------------------
  // Remove x-axis of UPPER PLOT
  h[0]->GetXaxis()->SetLabelOffset(999);
  h[0]->GetXaxis()->SetLabelSize(0);

  // Give y-axis title some offset to avoid overlapping with numbers
  h[0]->GetYaxis()->SetTitleOffset(1.45);
  // -------------------------------------------------------------------

  // Create output directory if it does not exist
  aux::CreateDirectory(fullpath);

  // Save pdf
  std::string fullfile = fullpath + name_ + ".pdf";
  c0.SaveAs(fullfile.c_str());

  // Save logscale pdf
  if (MINVAL > 0) {
    pad1->cd()->SetLogy();  // pad2 becomes the current pad
    for (const auto &i : indices(h)) { h[i]->GetYaxis()->SetRangeUser(std::max(MINVAL, MAXVAL * 1e-8), MAXVAL * 5); }
    fullfile = fullpath + name_ + "_logy" + ".pdf";
    c0.SaveAs(fullfile.c_str());
  }

  // Remove histograms
  for (const auto &i : indices(hR)) { delete hR[i]; }

  return chi2ndf;
}

h2Multiplet::h2Multiplet(const std::string &name, const std::string &labeltext, int N1,
                         double minval1, double maxval1, int N2, double minval2, double maxval2,
                         const std::vector<std::string> &legendtext) {
  name_ = name;

  N1_ = N1;
  N2_ = N2;

  minval1_ = minval1;
  maxval1_ = maxval1;
  minval2_ = minval2;
  maxval2_ = maxval2;

  legendtext_ = legendtext;

  // Initialize histogram vector container size
  h.resize(legendtext.size());

  for (const auto &i : indices(h)) {
    h[i] = new TH2D(Form("%s_%lu", name.c_str(), i), labeltext.c_str(), N1, minval1, maxval1, N2,
                    minval2, maxval2);
    h[i]->Sumw2();  // Error saving on
  }
}

void h2Multiplet::NormalizeAll(const std::vector<double> &cross_section,
                               const std::vector<double> &multiplier) {
  for (const auto &i : indices(h)) {
    double scale = 1.0;

    // Takes into account weighted and unweight (weight=1) filling
    // We take into account also under/overflow here with macro: (0, Nbins+1)
    const double integral = h[i]->Integral(0, h[i]->GetNbinsX() + 1, 0, h[i]->GetNbinsY() + 1);
    if (integral > 0) { scale /= integral; }

    // Binwidth
    const int    bin      = 1;
    const double binwidth = h[i]->GetXaxis()->GetBinWidth(bin) * h[i]->GetYaxis()->GetBinWidth(bin);
    if (binwidth > 0) { scale /= binwidth; }

    // Cross Section
    if (cross_section[i] > 0) { scale *= cross_section[i]; }

    // Additional scale factor
    scale *= multiplier[i];

    h[i]->Scale(scale);
  }
}

double h2Multiplet::SaveFig(const std::string &fullpath) const {

  TCanvas c0("c", "c", h.size() == 1 ? 525 : 750 / 3.0 * h.size(), 525); // scale canvas according to number of sources
  c0.Divide(h.size(), h.size() == 1 ? 1 : 2, 0.01, 0.02);                // [columns] x [rows]
  
  // Normalize by the first histogram
  //const double ZMAX = h[0]->GetMaximum();
  
  // Histograms on TOP ROW
  for (const auto &i : indices(h)) {
    c0.cd(i + 1);  // choose position
    c0.cd(i + 1)->SetRightMargin(0.13);

    h[i]->SetStats(0);
    h[i]->Draw("COLZ");
    h[i]->GetYaxis()->SetTitleOffset(1.3);
    //h[i]->GetZaxis()->SetRangeUser(0.0, ZMAX);
    h[i]->SetTitle(legendtext_[i].c_str());
  }

  // Ratio histograms on BOTTOM ROW
  std::vector<TH2D *> hR(h.size(), nullptr);

  if (h.size() > 1) { // Only if we have more than 1

  for (const auto &i : indices(h)) {
    c0.cd(i + 1 + h.size()); // Choose position
    c0.cd(i + 1 + h.size())->SetRightMargin(0.13);

    hR[i] = static_cast<TH2D*>(h[i]->Clone(Form("h2R_%lu", i)));

    hR[i]->Divide(h[0]);        // Divide by 0-th histogram
    hR[i]->GetYaxis()->SetTitleOffset(1.3);
    hR[i]->SetStats(0);         // No statistics on upper plot
    hR[i]->Draw("COLZ");
    hR[i]->GetZaxis()->SetRangeUser(0.0, 2.0);
    hR[i]->SetTitle(Form("Ratio: %s / %s", legendtext_[i].c_str(), legendtext_[0].c_str()));
  }

  }

  // -------------------------------------------------------------------
  // Create labels
  c0.cd();  // Important!
  std::shared_ptr<TPad> tpad;
  gra::rootstyle::TransparentPad(tpad);

  std::shared_ptr<TLatex> l1;
  std::shared_ptr<TLatex> l2;
  gra::rootstyle::MadeInFinland(l1, l2, 0.99);
  // -------------------------------------------------------------------

  // Create output directory if it does not exist
  aux::CreateDirectory(fullpath);

  // Save pdf
  std::string fullfile = fullpath + name_ + ".pdf";
  c0.SaveAs(fullfile.c_str());

  // Delete ratio histograms from memory
  if (h.size() > 1) {
    for (const auto &i : indices(hR)) { delete hR[i]; }
  }  
  return 0.0;
}


hProfMultiplet::hProfMultiplet(const std::string &name, const std::string &labeltext, int N,
                               double minval1, double maxval1, double minval2, double maxval2,
                               const std::vector<std::string> &legendtext) {
  name_ = name;

  N_       = N;
  minval1_ = minval1;
  maxval1_ = maxval1;

  minval2_ = minval2;
  maxval2_ = maxval2;

  legendtext_     = legendtext;
  legendposition_ = "northeast";

  // Initialize histogram vector container size
  h.resize(legendtext.size());

  for (const auto &i : indices(h)) {
    h[i] = new TProfile(Form("%s_%lu", name.c_str(), i), labeltext.c_str(), N, minval1, maxval1,
                        minval2, maxval2);
    h[i]->Sumw2();  // Error saving on
  }
}

double hProfMultiplet::SaveFig(const std::string &fullpath) const {
  
  // New colorscheme
  std::vector<int> color;
  std::vector<std::shared_ptr<TColor>> rootcolor;
  rootstyle::CreateColorMap(color, rootcolor);

  TCanvas c0("c", "c", 750, 800);

  // Upper plot will be in pad1
  std::shared_ptr<TPad> pad1 = std::make_shared<TPad>("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetBottomMargin(0.015);  // Upper and lower plot are joined
  // pad1->SetGridx();           // Vertical grid
  pad1->Draw();       // Draw the upper pad: pad1
  pad1->cd();         // pad1 becomes the current pad
  h[0]->SetStats(0);  // No statistics on upper plot

  // Find maximum value for y-range limits
  double MAXVAL = 0.0;
  double MINVAL = 1e32;
  for (const auto &i : indices(h)) {
    if (h[i]->GetMaximum() > MAXVAL) { MAXVAL = h[i]->GetMaximum(); }
    if (h[i]->GetMinimum() < MINVAL) { MINVAL = h[i]->GetMinimum(); }
  }
  
  // Loop over histograms
  for (const auto &i : indices(h)) {
    h[i]->SetLineColor(color[i]);
    h[i]->SetMarkerColor(color[i]);
    h[i]->SetMarkerStyle(20);
    h[i]->SetMarkerSize(0.5);
    h[i]->GetYaxis()->SetRangeUser(MINVAL, MAXVAL * 1.25);

    h[i]->Draw("L SAME");
  }

  // -------------------------------------------------------------------
  // LEGEND
  // North-East
  double x1, x2, y1, y2 = 0.0;
  GetLegendPosition(h.size(), x1, x2, y1, y2, legendposition_);
  std::shared_ptr<TLegend> legend = std::make_shared<TLegend>(x1, y1, x2, y2);

  legend->SetFillColor(0);  // White background
  // legend->SetBorderSize(0); // No box

  // Add legend entries
  for (const auto &i : indices(h)) { legend->AddEntry(h[i], legendtext_[i].c_str()); }
  legend->Draw();

  // -------------------------------------------------------------------
  // GRANIITTI text
  std::shared_ptr<TLatex> l1;
  std::shared_ptr<TLatex> l2;
  gra::rootstyle::MadeInFinland(l1, l2);
  // -------------------------------------------------------------------

  // -------------------------------------------------------------------
  // Ratio plots
  c0.cd();
  std::shared_ptr<TPad> pad2 = std::make_shared<TPad>("pad2", "pad2", 0, 0.05, 1, 0.3);
  pad2->SetTopMargin(0.025);
  pad2->SetBottomMargin(0.25);
  pad2->SetGridx();  // vertical grid
  pad2->Draw();
  pad2->cd();  // pad2 becomes the current pad

  // *** Ratio histograms ***
  std::vector<TProfile *> hxP(h.size(), nullptr);
  std::vector<TH1D *>     hR(h.size(),  nullptr);
  TH1D* hx0 = nullptr;

  for (const auto &i : indices(h)) {
    hxP[i] = static_cast<TProfile*>(h[i]->Clone(Form("ratio_%lu", i)));

    // --------------------------------------------------------------
    // To get proper errors, we need to use TH1 instead of TProfile
    hR[i] = hxP[i]->ProjectionX();
    
    if (i == 0) {  // this is needed as the denominator
    hx0 = static_cast<TH1D*>(hxP[i]->Clone("denominator"));
    }

    // Set same colors as above
    hR[i]->SetLineColor(color[i]);
    hR[i]->SetMarkerColor(color[i]);
    hR[i]->SetMarkerStyle(20);
    hR[i]->SetMarkerSize(0.5);

    // --------------------------------------------------------------

    // Get ratio to 0-th histogram
    hR[i]->Divide(hx0);

    hR[i]->SetMinimum(0.0);  // y-range
    hR[i]->SetMaximum(2.0);  //
    hR[i]->SetStats(0);      // no statistics box
    hR[i]->Draw("same");     // ratio plot

    // Ratio plot (h3) settings
    hR[i]->SetTitle("");  // remove title

    // Y axis ratio plot settings
    hR[i]->GetYaxis()->SetTitle("Ratio");
    hR[i]->GetYaxis()->CenterTitle();
    hR[i]->GetYaxis()->SetNdivisions(505);
    hR[i]->GetYaxis()->SetTitleSize(20);
    hR[i]->GetYaxis()->SetTitleFont(43);
    hR[i]->GetYaxis()->SetTitleOffset(1.55);
    hR[i]->GetYaxis()->SetLabelFont(43);  // in pixel (precision 3)
    hR[i]->GetYaxis()->SetLabelSize(15);

    // X axis ratio plot settings
    hR[i]->GetXaxis()->SetTitleSize(20);
    hR[i]->GetXaxis()->SetTitleFont(43);
    hR[i]->GetXaxis()->SetTitleOffset(4.);
    hR[i]->GetXaxis()->SetLabelFont(43);  // in pixel (precision 3)
    hR[i]->GetXaxis()->SetLabelSize(15);
  }

  // Draw horizontal line
  const double           ymax = 1.0;
  std::shared_ptr<TLine> line = std::make_shared<TLine>(minval1_, ymax, maxval1_, ymax);
  line->SetLineColor(15);
  line->SetLineWidth(2.0);
  line->Draw();

  // -------------------------------------------------------------------
  // Remove x-axis of UPPER PLOT
  h[0]->GetXaxis()->SetLabelOffset(999);
  h[0]->GetXaxis()->SetLabelSize(0);

  // Give y-axis title some offset to avoid overlapping with numbers
  h[0]->GetYaxis()->SetTitleOffset(1.45);
  // -------------------------------------------------------------------

  // Create output directory if it does not exist
  aux::CreateDirectory(fullpath);

  // Save pdf
  std::string fullfile = fullpath + name_ + ".pdf";
  c0.SaveAs(fullfile.c_str());

  // Save logscale pdf
  if (MINVAL > 0) {
    pad1->cd()->SetLogy();  // pad2 becomes the current pad
    fullfile = fullpath + name_ + "_logy" + ".pdf";
    c0.SaveAs(fullfile.c_str());
  }

  // Remove histograms
  for (const auto &i : indices(hR)) {
    delete hxP[i];
    delete hR[i];
  }
  delete hx0;

  return 0.0;
}

}  // gra namespace ends
