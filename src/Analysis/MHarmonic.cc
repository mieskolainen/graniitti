// "Meta class" for Spherical Harmonics Expansion - contains
// dataprocessing/Minuit/ROOT plotting functions.
//
// All spherical harmonic processing is separated to MSpherical namespace.
//
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <complex>
#include <iostream>
#include <random>
#include <regex>
#include <string>
#include <vector>

// ROOT
#include "TCanvas.h"
#include "TColor.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLorentzVector.h"
#include "TMinuit.h"
#include "TMultiGraph.h"
#include "TProfile.h"
#include "TStyle.h"

// Own
#include "Graniitti/Analysis/MHarmonic.h"
#include "Graniitti/Analysis/MROOT.h"
#include "Graniitti/M4Vec.h"
#include "Graniitti/MAux.h"
#include "Graniitti/MKinematics.h"
#include "Graniitti/MMath.h"
#include "Graniitti/MPDG.h"
#include "Graniitti/MSpherical.h"

// Eigen
#include <Eigen/Dense>

using gra::aux::indices;
using gra::math::msqrt;
using gra::math::zi;
using gra::math::PI;
using gra::math::pow2;

namespace gra {
// Constructor
MHarmonic::MHarmonic() {}

// Initialize
void MHarmonic::Init(const HPARAM &hp) {
  param = hp;
  NCOEF = (param.LMAX + 1) * (param.LMAX + 1);

  // ------------------------------------------------------------------
  // SETUP parameters of the fit
  ACTIVE.resize(NCOEF);  // Active moments
  ACTIVENDF = 0;

  for (int l = 0; l <= param.LMAX; ++l) {
    for (int m = -l; m <= l; ++m) {
      const int index = gra::spherical::LinearInd(l, m);
      ACTIVE[index]   = true;

      // FIX ODD MOMENTS TO ZERO
      if (param.REMOVEODD && ((l % 2) != 0)) { ACTIVE[index] = false; }
      // FIX NEGATIVE M TO ZERO
      if (param.REMOVENEGATIVEM && (m < 0)) { ACTIVE[index] = false; }
      if (ACTIVE[index]) { ++ACTIVENDF; }
    }
  }

  // ------------------------------------------------------------------
  std::cout << std::endl;
  param.Print();
  std::cout << std::endl;
  // ------------------------------------------------------------------

  // Init arrays used by Minuit function
  t_lm       = std::vector<double>(NCOEF, 0.0);
  t_lm_error = std::vector<double>(NCOEF, 0.0);
  errmat     = MMatrix<double>(NCOEF, NCOEF, 0.0);
  covmat     = MMatrix<double>(NCOEF, NCOEF, 0.0);

  // Test functions
  gra::spherical::TestSphericalIntegrals(param.LMAX);

  std::cout << rang::fg::yellow << "<Spherical Harmonic Based (costheta,phi)_r.f. "
                                   "Decomposition and Efficiency inversion>"
            << std::endl
            << std::endl;
  std::cout << "TERMINOLOGY:" << rang::fg::reset << std::endl;
  std::cout << "  {G} Generated == Events in the angular flat phase space (no "
               "cuts on final "
               "states, only on the system)"
            << std::endl;
  std::cout << "  {F} Fiducial  == Events within the strict fiducial "
               "(geometric-kinematic) "
               "final state phase space (cuts on final states)"
            << std::endl;
  std::cout << "  {S} Selected  == Events after the detector efficiency losses" << std::endl;
  std::cout << std::endl;
  std::cout << "{S} subset of {F} subset of {G} (this strict hierarchy might "
               "be violated in "
               "some special cases)"
            << std::endl;
  std::cout << std::endl;
  std::cout << "  The basic idea is to define {G} such that minimal extrapolation " << std::endl;
  std::cout << "  is required from the strict fiducial (geometric) phase space {F}. " << std::endl;
  std::cout << "  Flatness requirement of {G} is strictly required to represent "
               "moments in "
               "an unmixed basis (non-flat phase space <=> geometric moment mixing)."
            << std::endl;
  std::cout << std::endl;
  std::cout << rang::fg::yellow << "EXAMPLE OF A FORMALLY VALID DEFINITION:" << rang::fg::reset
            << std::endl;
  std::cout << "  G = {|Y(system)| < 0.9}" << std::endl;
  std::cout << "  F = {|eta(pi)|   < 0.9 && pt(pi) > 0.1 GeV}" << std::endl;
  std::cout << std::endl << std::endl;
  std::cout << "Note also the rotation between inactive coefficients and "
               "active one, due to "
               "moment mixing."
            << std::endl;
  std::cout << std::endl << std::endl;

  // pause(5);
}

// Destructor
MHarmonic::~MHarmonic() {}

bool MHarmonic::PrintLoop(const std::string &output) const {
  /*
  // ------------------------------------------------------------------
  // Print out results for EXTERNAL analysis
  FILE* fp;

  char buff[100];
  snprintf(buff, sizeof(buff), "./fits/%s_%s.m", output.c_str(), FRAME.c_str());
  std::string outputfile = buff;

  fp = fopen(outputfile.c_str(), "w");

  if (fp == NULL) {
                  fprintf(stderr, "Cannot open outputfile %s !\n", outputfile.c_str());
                  return false;
  }

  fprintf(fp, "FRAME = %s; \n",     FRAME.c_str());
  fprintf(fp, "param.LMAX = %d; \n",      param.LMAX);
  fprintf(fp, "LAMBDA = %0.3f; \n", LAMBDA);
  fprintf(fp, "ACTIVE = [");

  for (std::size_t i = 0; i < ACTIVE.size(); ++i) {
                  fprintf(fp, "%d ", (int)ACTIVE[i]);
  }
  fprintf(fp, "]; \n\n");
  fprintf(fp, "M_CELL = [");

  for (std::size_t i = 0; i < masscenter.size(); ++i) {
                  fprintf(fp, "%0.3f ", masscenter[i]);
  }
  fprintf(fp, "]; \n\n");

  fprintf(fp, "A{1} = [");
  PrintMatrix(fp, fla.t_lm_EML);
  fprintf(fp, "];\n");
  fprintf(fp, "A{2} = [");
  PrintMatrix(fp, fid.t_lm_EML);
  fprintf(fp, "];\n");
  fprintf(fp, "A{3} = [");
  PrintMatrix(fp, det.t_lm_EML);
  fprintf(fp, "];\n");

  fprintf(fp, "A{4} = [");
  PrintMatrix(fp, fla.t_lm_MPP);
  fprintf(fp, "];\n");
  fprintf(fp, "A{5} = [");
  PrintMatrix(fp, fid.t_lm_MPP);
  fprintf(fp, "];\n");
  fprintf(fp, "A{6} = [");
  PrintMatrix(fp, det.t_lm_MPP);
  fprintf(fp, "];\n");

  fclose(fp);
  */
  return true;
}

void MHarmonic::PlotAll(const std::string &outputpath) const {
  for (const auto &OBSERVABLE : {0, 1, 2}) {
    // **** EFFICIENCY DECOMPOSITION ****
    Plot1DEfficiency(OBSERVABLE, outputpath);

    // **** ALGEBRAIC INVERSE MOMENTS ****
    PlotFigures(det, OBSERVABLE, "h_{Moments}[MPP]<det>", 33, outputpath);
    PlotFigures(fid, OBSERVABLE, "h_{Moments}[MPP]<fid>", 33, outputpath);
    PlotFigures(fla, OBSERVABLE, "h_{Moments}[MPP]<fla>", 33, outputpath);

    // **** EXTENDED MAXIMUM LIKELIHOOD INVERSE MOMENTS ****
    if (param.EML) {
      PlotFigures(det, OBSERVABLE, "h_{Moments}[EML]<det>", 33, outputpath);
      PlotFigures(fid, OBSERVABLE, "h_{Moments}[EML]<fid>", 33, outputpath);
      PlotFigures(fla, OBSERVABLE, "h_{Moments}[EML]<fla>", 33, outputpath);
    }
  }

  // 2D
  std::vector<std::vector<int>> OBSERVABLES = {{0, 1}, {0, 2}, {1, 2}};  // Different pairs

  for (const auto &OBSERVABLE2 : OBSERVABLES) {
    // **** EFFICIENCY DECOMPOSITION ****
    PlotFigures2D(fla, OBSERVABLE2, "{Response}[FLAT_REFERENCE]<fla>", 17, outputpath);
    PlotFigures2D(fid, OBSERVABLE2, "{Response}[FIDUCIAL_ACCEPTANCE]<fid>", 33, outputpath);
    PlotFigures2D(det, OBSERVABLE2, "{Response}[ACCEPTANCE_x_EFFICIENCY]<det>", 29, outputpath);

    // **** ALGEBRAIC INVERSE MOMENTS ****
    PlotFigures2D(fla, OBSERVABLE2, "{Moments}[MPP]<fla>", 17, outputpath);
    PlotFigures2D(fid, OBSERVABLE2, "{Moments}[MPP]<fid>", 33, outputpath);
    PlotFigures2D(det, OBSERVABLE2, "{Moments}[MPP]<det>", 29, outputpath);

    // **** EXTENDED MAXIMUM LIKELIHOOD INVERSE MOMENTS ****
    if (param.EML) {
      PlotFigures2D(fla, OBSERVABLE2, "{Moments}[EML]<fla>", 17, outputpath);
      PlotFigures2D(fid, OBSERVABLE2, "{Moments}[EML]<fid>", 33, outputpath);
      PlotFigures2D(det, OBSERVABLE2, "{Moments}[EML]<det>", 29, outputpath);
    }
  }

  const int OBSERVABLE = 0;
  Plot2DExpansion(fid, OBSERVABLE, "h_{Moments}[MPP]<fid>", 33, outputpath);
  Plot2DExpansion(fla, OBSERVABLE, "h_{Moments}[MPP]<fla>", 33, outputpath);
  if (param.EML) {
    Plot2DExpansion(fid, OBSERVABLE, "h_{Moments}[EML]<fid>", 33, outputpath);
    Plot2DExpansion(fla, OBSERVABLE, "h_{Moments}[EML]<fla>", 33, outputpath);
  }
}

// Synthesized (costheta,phi) plots
//
void MHarmonic::Plot2DExpansion(
    const std::map<gra::spherical::Meta, MTensor<gra::spherical::SH>> &tensor,
    unsigned int OBSERVABLE, const std::string &TYPESTRING, int barcolor,
    const std::string &outputpath) const {
  // ------------------------------------------------------------------
  // Extract name strings

  // find {string}
  std::smatch sma;
  std::regex_search(TYPESTRING, sma, std::regex(R"(\{.*?\})"));  // R"()" for Raw string literals
  std::string DATAMODE = sma[0];
  DATAMODE             = DATAMODE.substr(1, DATAMODE.size() - 2);

  // find <string>
  std::smatch smb;
  std::regex_search(TYPESTRING, smb, std::regex(R"(\<.*?\>)"));  // R"()" for Raw string literals
  std::string SPACE = smb[0];
  SPACE             = SPACE.substr(1, SPACE.size() - 2);

  // find [string]
  std::smatch smc;
  std::regex_search(TYPESTRING, smc, std::regex(R"(\[.*?\])"));  // R"()" for Raw string literals
  std::string ALGO = smc[0];
  ALGO             = ALGO.substr(1, ALGO.size() - 2);

  // ------------------------------------------------------------------
  // f(cos(theta),phi; M) synthesis

  std::size_t N = 50;

  // Cos(theta) and phi
  std::vector<double> costheta;
  std::vector<double> phi;

  std::size_t BINS = 0;

  // TH2 for each
  std::vector<std::vector<TH2D *>> h2;

  // Loop over fla
  std::size_t source_ind = 0;
  for (const auto &source : tensor) {
    // legendstrs[ind] = source.first.LEGEND;
    BINS = source.second.size(OBSERVABLE);

    // Add histograms
    h2.push_back(std::vector<TH2D *>(BINS, NULL));

    // Loop over observable

    for (std::size_t bin = 0; bin < BINS; ++bin) {

      // Get x-axis point
      //const double value = grid[OBSERVABLE][bin].center();

      // Set indices {0,0,0, ..., 0}
      std::vector<std::size_t> cell(grid.size(), 0);
      cell[OBSERVABLE] = bin;

      // Synthesize distribution
      MMatrix<double> Z;

      if (ALGO == "MPP") {
        Z = spherical::Y_real_synthesize(source.second(cell).t_lm_MPP, ACTIVE, N, costheta, phi,
                                         true);
      } else if (ALGO == "EML") {
        Z = spherical::Y_real_synthesize(source.second(cell).t_lm_EML, ACTIVE, N, costheta, phi,
                                         true);
      }

      // Create histogram
      const double      EPS  = 1e-3;
      const std::string name = "h2_" + std::to_string(source_ind) + "_" + std::to_string(bin);
      h2[source_ind][bin]    = new TH2D(name.c_str(), "; cos #theta; #phi (rad)", N, -(1 + EPS),
                                     1 + EPS, N, -(math::PI + EPS), math::PI + EPS);

      for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = 0; j < N; ++j) {
          h2[source_ind][bin]->Fill(costheta[i], phi[j], Z[i][j]);
        }
      }
    }

    ++source_ind;
  }  // loop over sources

  // Loop over bins
  for (std::size_t bin = 0; bin < BINS; ++bin) {
    TCanvas *c1 = new TCanvas("c1", "c1", 200 + tensor.size() * 400, 500);  // horizontal, vertical
    c1->Divide(tensor.size(), 1, 0.002, 0.001);

    std::size_t source_ind = 0;
    std::string FRAME;
    for (const auto &source : tensor) {
      c1->cd(source_ind + 1);
      c1->cd(source_ind + 1)->SetRightMargin(0.13);
      //
      FRAME = source.first.FRAME;
      h2[source_ind][bin]->SetTitle(source.first.LEGEND.c_str());
      //
      h2[source_ind][bin]->GetXaxis()->CenterTitle();
      h2[source_ind][bin]->GetZaxis()->SetRangeUser(0.1, 1.0);
      h2[source_ind][bin]->Draw("COLZ");
      //

      // --------------------------------------------------------------
      // Draw Lorentz FRAME string on left
      TText *t2 = new TText(0.4, 0.85,
                            Form("%s : [%0.2f, %0.2f] %s", FRAME.c_str(), grid[OBSERVABLE][bin].min,
                                 grid[OBSERVABLE][bin].max, xlabels[OBSERVABLE].c_str()));
      t2->SetNDC();
      t2->SetTextAlign(22);
      t2->SetTextColor(kBlack);
      t2->SetTextFont(43);
      t2->SetTextSize(18);
      // t2->SetTextAngle(45);
      t2->Draw("same");
      // --------------------------------------------------------------

      ++source_ind;
    }
    //
    aux::CreateDirectory("./figs");
    aux::CreateDirectory("./figs/harmonicfit");
    aux::CreateDirectory("./figs/harmonicfit/" + outputpath + "/synthesis");
    //
    const std::string subpath =
        SPACE + "_OBS_" + std::to_string(OBSERVABLE) + "_" + FRAME + "_" + ALGO;
    const std::string fullpath = "./figs/harmonicfit/" + outputpath + "/synthesis/" + subpath;
    aux::CreateDirectory(fullpath);
    c1->Print(Form("%s/%04lu.pdf", fullpath.c_str(), bin));
    //
    delete c1;
  }

  // Delete all histograms
  for (std::size_t i = 0; i < h2.size(); ++i) {
    for (std::size_t j = 0; j < h2[i].size(); ++j) { delete h2[i][j]; }
  }
}

void GetLegendPosition2(unsigned int N, double &x1, double &x2, double &y1, double &y2,
                        const std::string &legendposition) {
  // North-East
  x1 = 0.63;
  x2 = x1 + 0.18;
  y1 = 0.75 - 0.01 * N;

  // South-East
  if (legendposition.compare("southeast") == 0) { y1 = 0.10 - 0.01 * N; }
  y2 = y1 + 0.05 * N;  // Scale by the number of histograms
}

// LM-Expansion plots as a function of observables
//
void MHarmonic::PlotFigures(
    const std::map<gra::spherical::Meta, MTensor<gra::spherical::SH>> &tensor,
    unsigned int OBSERVABLE, const std::string &TYPESTRING, int barcolor,
    const std::string &outputpath) const {
  // ------------------------------------------------------------------
  std::shared_ptr<TCanvas> c1;
  gra::rootstyle::AutoGridCanvas(c1, ACTIVENDF);
  // ------------------------------------------------------------------

  if (OBSERVABLE > xlabels.size() - 1) {
    throw std::invalid_argument("MHarmonic::PlotFigures: Unknown observable " +
                                std::to_string(OBSERVABLE));
  }
  const std::string xlabel = xlabels[OBSERVABLE];

  // ------------------------------------------------------------------
  // Extract name strings

  // find {string}
  std::smatch sma;
  std::regex_search(TYPESTRING, sma, std::regex(R"(\{.*?\})"));  // R"()" for Raw string literals
  std::string DATAMODE = sma[0];
  DATAMODE             = DATAMODE.substr(1, DATAMODE.size() - 2);

  // find <string>
  std::smatch smb;
  std::regex_search(TYPESTRING, smb, std::regex(R"(\<.*?\>)"));  // R"()" for Raw string literals
  std::string SPACE = smb[0];
  SPACE             = SPACE.substr(1, SPACE.size() - 2);

  // find [string]
  std::smatch smc;
  std::regex_search(TYPESTRING, smc, std::regex(R"(\[.*?\])"));  // R"()" for Raw string literals
  std::string ALGO = smc[0];
  ALGO             = ALGO.substr(1, ALGO.size() - 2);

  // ------------------------------------------------------------------

  // Loop over data sources
  std::size_t BINS = 0;
  int         ind  = 0;

  // Graphs for each data [source] x [lm-moment]
  std::vector<std::vector<TGraphErrors *>> gr(tensor.size(),
                                              std::vector<TGraphErrors *>(ACTIVENDF, NULL));

  // Legend titles
  std::vector<std::string> legendstrs(tensor.size());

  // Minimum and maximum y-values for each plot
  std::vector<double> MINVAL(ACTIVENDF, 1e32);
  std::vector<double> MAXVAL(ACTIVENDF, -1e32);

  // Turn of horizontal errors
  gStyle->SetErrorX(0);

  std::vector<TMultiGraph *> mg(ACTIVENDF, NULL);
  for (std::size_t k = 0; k < ACTIVENDF; ++k) { mg[k] = new TMultiGraph(); }

  std::string              FRAME;
  std::vector<std::string> TITLES;

  // Y-axis title
  std::string yaxis_label = "Events / bin";

  for (const auto &source : tensor) {
    legendstrs[ind] = source.first.LEGEND;
    BINS            = source.second.size(OBSERVABLE);

    // Loop over moments
    int k = 0;

    double SCALE = source.first.SCALE;

    if (SCALE < 0 && source.first.YAXIS == "") {
      yaxis_label = "Normalized to 1";
    }
    if (source.first.YAXIS != "") {
      yaxis_label  = source.first.YAXIS;
    }

    for (int l = 0; l <= param.LMAX; ++l) {
      for (int m = -l; m <= l; ++m) {
        const int index = gra::spherical::LinearInd(l, m);

        if (!ACTIVE[index]) { continue; }  // Not active

        // Set canvas position
        c1->cd(k + 1);

        // Loop over observable
        double x[BINS]     = {0.0};
        double y[BINS]     = {0.0};
        double x_err[BINS] = {0.0};
        double y_err[BINS] = {0.0};

        for (std::size_t bin = 0; bin < BINS; ++bin) {
        
          // Get x-axis point
          x[bin] = grid[OBSERVABLE][bin].center();

          // Set indices {0,0,0, ..., 0}
          std::vector<std::size_t> cell(grid.size(), 0);
          cell[OBSERVABLE] = bin;

          // CHOOSE DATAMODE
          if        (ALGO == "MPP") {
            y[bin]     = source.second(cell).t_lm_MPP[index];
            y_err[bin] = source.second(cell).t_lm_MPP_error[index];
          } else if (ALGO == "EML") {
            y[bin]     = source.second(cell).t_lm_EML[index];
            y_err[bin] = source.second(cell).t_lm_EML_error[index];
          } else {
            throw std::invalid_argument(
                "MHarmonic::PlotFigures: Unknown input: "
                "DATAMODE = " +
                DATAMODE + " ALGO = " + ALGO);
          }
        }

        // ------------------------------------------------------
        // Scaling (e.g. luminosity)
        // Normalize lm = 00 to sum to 1 -> then apply to all
        if ((l == 0 && m == 0) && SCALE < 0.0) {
          double sum = 0.0;
          for (std::size_t bin = 0; bin < BINS; ++bin) { sum += y[bin]; }
          if (sum > 0) { SCALE = 1.0 / sum; }
        }

        // Apply scale by user
        for (std::size_t bin = 0; bin < BINS; ++bin) {
          y[bin] *= SCALE;
          y_err[bin] *= SCALE;
        }
        // ------------------------------------------------------

        // Save maximum for visualization
        for (std::size_t bin = 0; bin < BINS; ++bin) {
          MAXVAL[k] = y[bin] > MAXVAL[k] ? y[bin] : MAXVAL[k];
          MINVAL[k] = y[bin] < MINVAL[k] ? y[bin] : MINVAL[k];
        }
        // Data displayed using TGraphErrors
        gr[ind][k] = new TGraphErrors(BINS, x, y, x_err, y_err);

        // Colors
        gr[ind][k]->SetMarkerColor(colors[ind]);
        gr[ind][k]->SetLineColor(colors[ind]);
        gr[ind][k]->SetFillColor(colors[ind]);
        gr[ind][k]->SetLineWidth(2.0);
        gr[ind][k]->SetMarkerStyle(21);  // square
        gr[ind][k]->SetMarkerSize(0.3);

        FRAME  = source.first.FRAME;
        TITLES = source.first.TITLES;

        // Add to the multigraph
        mg[k]->Add(gr[ind][k]);

        ++k;

      }  // over m
    }    // over l

    ++ind;
  }  // Loop over sources

  // Draw multigraph
  unsigned int k = 0;

  // Aux variables
  std::shared_ptr<TPad> tpad;
  std::shared_ptr<TLatex> l1;
  std::shared_ptr<TLatex> l2;

  for (int l = 0; l <= param.LMAX; ++l) {
    for (int m = -l; m <= l; ++m) {
      const int index = gra::spherical::LinearInd(l, m);

      if (!ACTIVE[index]) { continue; }  // Not active

      // Set canvas position
      c1->cd(k + 1);

      // First draw, then setup (otherwise crash)
      mg[k]->Draw("AC*");

      // Title and y-axis
      if (k == 0) {

        // Title
        if (SPACE == "det") {
          mg[k]->SetTitle(Form("%s | #it{lm} = <%d,%d>", TITLES[0].c_str(), l, m));
        }
        if (SPACE == "fid") {
          mg[k]->SetTitle(Form("%s | #it{lm} = <%d,%d>", TITLES[1].c_str(), l, m));
        }
        if (SPACE == "fla") {
          mg[k]->SetTitle(Form("%s | #it{lm} = <%d,%d>", TITLES[2].c_str(), l, m));
        }

        // y-axis (for some ROOT reason, this needs to be always after SetTitle)
        mg[k]->GetYaxis()->SetTitle(yaxis_label.c_str());
        gPad->SetLeftMargin(0.15); // 15 per cent of pad for left margin, default is 10%
        mg[k]->GetYaxis()->SetTitleOffset(1.25);

      } else {
        mg[k]->SetTitle(Form("Moment #it{lm} = <%d,%d>", l, m));
      }

      // Set x-axis
      mg[k]->GetXaxis()->SetTitle(xlabel.c_str());
      // gStyle->SetBarWidth(0.5);
      // mg[k]->SetFillStyle(0);

      mg[k]->GetXaxis()->SetTitleSize(0.05);
      mg[k]->GetXaxis()->SetLabelSize(0.05);

      mg[k]->GetYaxis()->SetTitleSize(0.05);
      mg[k]->GetYaxis()->SetLabelSize(0.05);

      // Set y-axis
      // mg[k]->GetYaxis()->SetTitle("Intensity");

      gStyle->SetTitleFontSize(0.08);

      if (k == 0) {
        gStyle->SetTitleW(0.95);  // width percentage
      }

      // Y-axis range
      if (k == 0) {
        mg[k]->GetHistogram()->SetMaximum(MAXVAL[k] * 1.1);
        mg[k]->GetHistogram()->SetMinimum(0.0);
      } else {
        // Skip the first element (lm=00)
        const double maxval = *std::max_element(std::begin(MAXVAL) + 1, std::end(MAXVAL));
        const double minval = *std::min_element(std::begin(MINVAL) + 1, std::end(MINVAL));
        const double bound  = std::max(std::abs(minval), std::abs(maxval));

        mg[k]->GetHistogram()->SetMaximum(bound * 1.1);
        mg[k]->GetHistogram()->SetMinimum(-bound * 1.1);
      }

      // --------------------------------------------------------------
      // Draw Lorentz FRAME string on left
      TText *t2 = new TText(0.225, 0.825, FRAME.c_str());
      t2->SetNDC();
      t2->SetTextAlign(22);
      t2->SetTextColor(kRed + 2);
      t2->SetTextFont(43);
      t2->SetTextSize(std::ceil(1.0 / msqrt(ACTIVENDF)) * 16);
      // t2->SetTextAngle(45);
      t2->Draw("same");
      // --------------------------------------------------------------

      // --------------------------------------------------------------
      // Draw horizontal line
      if (k != 0) {
        TLine *line = new TLine(grid[OBSERVABLE][0].min, 0.0,
                                grid[OBSERVABLE][grid[OBSERVABLE].size() - 1].max, 0.0);
        line->SetLineColor(kBlack);
        line->SetLineWidth(1.0);
        line->Draw("same");
      }
      // --------------------------------------------------------------

      // --------------------------------------------------------------
      // Who made it
      if (k == ACTIVENDF - 1) {
        // New pad on top of all
        c1->cd();  // Important!
        gra::rootstyle::TransparentPad(tpad);

        const double xpos = 0.99;
        gra::rootstyle::MadeInFinland(l1, l2, xpos);
      }
      // --------------------------------------------------------------

      ++k;
    }
  }

  // ------------------------------------------------------------------
  // Draw legend to the upper left most
  c1->cd(1);

  // Create legend
  double            x1, x2, y1, y2 = 0.0;
  const std::string legendposition = "northeast";
  GetLegendPosition2(tensor.size(), x1, x2, y1, y2, legendposition);
  TLegend *legend = new TLegend(x1, y1, x2, y2);
  legend->SetFillColor(0);   // White background
  legend->SetBorderSize(0);  // No box
  legend->SetTextSize(0.035);

  // Add legend entries
  for (const auto &i : indices(gr)) { legend->AddEntry(gr[i][0], legendstrs[i].c_str()); }

  // Draw legend
  legend->Draw("same");
  // ------------------------------------------------------------------

  // Require that we have data
  if (BINS > 1) {
    aux::CreateDirectory("./figs");
    aux::CreateDirectory("./figs/harmonicfit");
    aux::CreateDirectory("./figs/harmonicfit/" + outputpath);
    const std::string subpath  = "OBS_" + std::to_string(OBSERVABLE) + "_" + FRAME;
    const std::string fullpath = "./figs/harmonicfit/" + outputpath + "/" + subpath;
    aux::CreateDirectory(fullpath);
    c1->Print(Form("%s/%s.pdf", fullpath.c_str(), TYPESTRING.c_str()));

    // Merge pdfs using Ghostscript (gs)
    const std::string cmd = "gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=" + fullpath +
                            "/" + FRAME + "_merged.pdf " + fullpath + "/h*.pdf";
    if (system(cmd.c_str()) == -1) {
      throw std::invalid_argument("Error: Problem executing Ghostscript merge on pdfs!");
    }
  }

  /*
  for (std::size_t i = 0; i < gr.size(); ++i) {
                  delete gr[i];
  }
  */
}

// LM-Efficiency figures
//
void MHarmonic::Plot1DEfficiency(unsigned int OBSERVABLE, const std::string &outputpath) const {
  // ------------------------------------------------------------------
  std::shared_ptr<TCanvas> c1;
  gra::rootstyle::AutoGridCanvas(c1, ACTIVENDF);
  // ------------------------------------------------------------------

  if (OBSERVABLE > xlabels.size() - 1) {
    throw std::invalid_argument("MHarmonic::PlotFigures: Unknown observable " +
                                std::to_string(OBSERVABLE));
  }
  const std::string xlabel = xlabels[OBSERVABLE];

  // Graphs for each efficiency [level] x [lm-moment]
  std::vector<std::vector<TGraphErrors *>> gr(3, std::vector<TGraphErrors *>(ACTIVENDF, NULL));

  // Legend titles
  std::vector<std::string> legendstrs(3);

  // Minimum and maximum y-values for each plot
  std::vector<double> MINVAL(ACTIVENDF, 1e32);
  std::vector<double> MAXVAL(ACTIVENDF, -1e32);

  // Turn of horizontal errors
  gStyle->SetErrorX(0);

  std::vector<TMultiGraph *> mg(ACTIVENDF, NULL);
  for (std::size_t k = 0; k < ACTIVENDF; ++k) { mg[k] = new TMultiGraph(); }

  std::size_t BINS = det_DET.size(OBSERVABLE);

  // Read data
  std::string              FRAME;
  std::vector<std::string> TITLES;

  for (const auto &source : det) {
    FRAME  = source.first.FRAME;
    TITLES = source.first.TITLES;
    break;
  }

  for (std::size_t level = 0; level < 3; ++level) {
    // Loop over moments
    int k = 0;
    for (int l = 0; l <= param.LMAX; ++l) {
      for (int m = -l; m <= l; ++m) {
        const int index = gra::spherical::LinearInd(l, m);

        if (!ACTIVE[index]) { continue; }  // Not active

        // Set canvas position
        c1->cd(k + 1);

        // Loop over mass
        double x[BINS]     = {0.0};
        double y[BINS]     = {0.0};
        double x_err[BINS] = {0.0};
        double y_err[BINS] = {0.0};

        for (std::size_t bin = 0; bin < BINS; ++bin) {
          // Get x-axis point
          x[bin] = grid[OBSERVABLE][bin].center();

          // Set indices {0,0,0, ..., 0}
          std::vector<std::size_t> cell(grid.size(), 0);
          cell[OBSERVABLE] = bin;

          // CHOOSE DATAMODE
          if (level == 0) {
            y[bin]     = det_DET(cell).E_lm[index];
            y_err[bin] = det_DET(cell).E_lm_error[index];
          } else if (level == 1) {
            y[bin]     = fid_DET(cell).E_lm[index];
            y_err[bin] = fid_DET(cell).E_lm_error[index];
          } else if (level == 2) {
            y[bin]     = fla_DET(cell).E_lm[index];
            y_err[bin] = fla_DET(cell).E_lm_error[index];
          }

          // Save maximum for visualization
          MAXVAL[k] = y[bin] > MAXVAL[k] ? y[bin] : MAXVAL[k];
          MINVAL[k] = y[bin] < MINVAL[k] ? y[bin] : MINVAL[k];
        }

        // Data displayed using TGraphErrors
        gr[level][k] = new TGraphErrors(BINS, x, y, x_err, y_err);

        // Colors
        gr[level][k]->SetMarkerColor(colors[level]);
        gr[level][k]->SetLineColor(colors[level]);
        gr[level][k]->SetFillColor(colors[level]);
        gr[level][k]->SetLineWidth(2.0);
        gr[level][k]->SetMarkerStyle(21);  // square
        gr[level][k]->SetMarkerSize(0.3);

        // Add to the multigraph
        mg[k]->Add(gr[level][k]);

        ++k;

      }  // over m
    }    // over l

  }  // Loop over level

  // Draw multigraph
  unsigned int k = 0;

  std::shared_ptr<TPad> tpad;
  std::shared_ptr<TLatex>  l1;
  std::shared_ptr<TLatex>  l2;

  for (int l = 0; l <= param.LMAX; ++l) {
    for (int m = -l; m <= l; ++m) {
      const int index = gra::spherical::LinearInd(l, m);

      if (!ACTIVE[index]) { continue; }  // Not active

      // Set canvas position
      c1->cd(k + 1);

      // First draw, then setup (otherwise crash)
      mg[k]->Draw("AC*");

      // Title
      if (k == 0) {
        mg[k]->SetTitle(Form("Acceptance decomposition: #it{lm} = <%d,%d>", l, m));
      } else {
        mg[k]->SetTitle(Form("Moment #it{lm} = <%d,%d>", l, m));
      }

      // Set x-axis
      mg[k]->GetXaxis()->SetTitle(xlabel.c_str());
      // gStyle->SetBarWidth(0.5);
      // mg[k]->SetFillStyle(0);

      mg[k]->GetXaxis()->SetTitleSize(0.05);
      mg[k]->GetXaxis()->SetLabelSize(0.05);

      mg[k]->GetYaxis()->SetTitleSize(0.05);
      mg[k]->GetYaxis()->SetLabelSize(0.05);

      // Set y-axis
      // mg[k]->GetYaxis()->SetTitle("");

      gStyle->SetTitleFontSize(0.08);

      if (k == 0) {
        gStyle->SetTitleW(0.95);  // width percentage
      }
      // Y-axis range
      if (k == 0) {
        mg[k]->GetHistogram()->SetMaximum(1.2);
        mg[k]->GetHistogram()->SetMinimum(0.0);
      } else {
        // Skip the first element (lm=00)
        const double maxval = *std::max_element(std::begin(MAXVAL) + 1, std::end(MAXVAL));
        const double minval = *std::min_element(std::begin(MINVAL) + 1, std::end(MINVAL));
        const double bound  = std::max(std::abs(minval), std::abs(maxval));

        mg[k]->GetHistogram()->SetMaximum(bound * 1.1);
        mg[k]->GetHistogram()->SetMinimum(-bound * 1.1);
      }

      // Draw Lorentz FRAME string on left
      TText *t2 = new TText(0.225, 0.825, FRAME.c_str());
      t2->SetNDC();
      t2->SetTextAlign(22);
      t2->SetTextColor(kRed + 2);
      t2->SetTextFont(43);
      t2->SetTextSize(std::ceil(1.0 / msqrt(ACTIVENDF)) * 16);
      // t2->SetTextAngle(45);
      t2->Draw("same");

      // Draw horizontal line
      if (k != 0) {
        TLine *line = new TLine(grid[OBSERVABLE][0].min, 0.0,
                                grid[OBSERVABLE][grid[OBSERVABLE].size() - 1].max, 0.0);
        line->SetLineColor(kBlack);
        line->SetLineWidth(1.0);
        line->Draw("same");
      }

      // Who made it
      if (k == ACTIVENDF - 1) {
        // New pad on top of all
        c1->cd();  // Important!
        gra::rootstyle::TransparentPad(tpad);

        const double xpos = 0.99;
        gra::rootstyle::MadeInFinland(l1, l2, xpos);
      }
      ++k;
    }
  }

  // ------------------------------------------------------------------
  // Draw legend to the upper left most

  c1->cd(1);

  // Create legend
  double            x1, x2, y1, y2 = 0.0;
  const std::string legendposition = "southeast";
  GetLegendPosition2(3, x1, x2, y1, y2, legendposition);
  TLegend *legend = new TLegend(x1 - 0.3, y1 + 0.1, x2 - 0.3, y2 + 0.1);
  legend->SetFillColor(0);   // White background
  legend->SetBorderSize(0);  // No box
  legend->SetTextSize(0.04);

  // Add legend entries
  legend->AddEntry(gr[0][0], TITLES[0].c_str());
  legend->AddEntry(gr[1][0], TITLES[1].c_str());
  legend->AddEntry(gr[2][0], TITLES[2].c_str());

  // Draw legend
  legend->Draw("same");

  // ------------------------------------------------------------------

  // Require that we have data
  if (BINS > 1) {
    aux::CreateDirectory("./figs");
    aux::CreateDirectory("./figs/harmonicfit");
    aux::CreateDirectory("./figs/harmonicfit/" + outputpath);
    const std::string subpath  = "OBS_" + std::to_string(OBSERVABLE) + "_" + FRAME;
    const std::string fullpath = "./figs/harmonicfit/" + outputpath + "/" + subpath;
    aux::CreateDirectory(fullpath);
    c1->Print(Form("%s/h_Response.pdf", fullpath.c_str()));
  }

  /*
  for (std::size_t i = 0; i < gr.size(); ++i) {
                  delete gr[i];
  }
  */
}

// Print 2D-figures
void MHarmonic::PlotFigures2D(
    const std::map<gra::spherical::Meta, MTensor<gra::spherical::SH>> &tensor,
    const std::vector<int> &OBSERVABLE2, const std::string &TYPESTRING, int barcolor,
    const std::string &outputpath) const {
  /*
                  TCanvas* c1 = new TCanvas("c1", "c1", 700, 500); // horizontal,
    vertical
                  c1->Divide(sqrt(NCOEF), std::sqrt(NCOEF), 0.002, 0.001);

    // Loop over moments
    std::vector<TH2D*> gr(NCOEF, NULL);

    std::vector<std::string> label(2);
    for (std::size_t i = 0; i < 2; ++i) {

          if      (OBSERVABLE[i] == 0) {
            label[i] = "M (GeV)";
          }
          else if (OBSERVABLE[i] == 1) {
            label[i] = "P_{T} (GeV)";
          }
          else if (OBSERVABLE[i] == 2) {
            label[i] = "Y";
          } else {
            throw std::invalid_argument("MHarmonic::PlotFigures2D: Unknown
    observable " +
    std::to_string(OBSERVABLE[i]));
          }
    }
    int k = 0;

    std::vector<std::size_t> BINS = {tensor.size(OBSERVABLE[0]),
    tensor.size(OBSERVABLE[1])};

    for (int l = 0; l <= param.LMAX; ++l) {
          for (int m = -l; m <= l; ++m) {
            const int index = gra::spherical::LinearInd(l, m);

            // Set canvas position
            c1->cd(k + 1);

            // Loop over mass
            double x[BINS[0]] = {0};
            double y[BINS[1]] = {0};
            double z[BINS[0]][BINS[1]] = {{0}};
            //double x_err[BINS] = {1e-2};
            //double y_err[BINS] = {1e-2};

            for (std::size_t k1 = 0; k1 < BINS[0]; ++k1) {
                  x[k1] = grid[OBSERVABLE[0]][k1].center();

                  for (std::size_t k2 = 0; k2 < BINS[1]; ++k2) {
                    y[k2] = grid[OBSERVABLE[1]][k2].center();

                    // Set indices {0,0,0, ..., 0}
                    std::vector<std::size_t> cell(xcenter.size(),0);
                    cell[OBSERVABLE[0]] = k1;
                    cell[OBSERVABLE[1]] = k2;

                    // CHOOSE DATAMODE
                    if      (DATAMODE == "t_lm_MPP") {
                          z[k1][k2] = tensor(cell).t_lm_MPP[index];
                    }
                    else if (DATAMODE == "t_lm_EML") {
                          z[k1][k2] = tensor(cell).t_lm_EML[index];
                    }
                    else if (DATAMODE == "E_lm") {
                          z[k1][k2] = tensor(cell).E_lm[index];
                    } else {
                          throw std::invalid_argument("MHarmonic::PlotFigures: Unknown MODE
    = " +
    DATAMODE);
                    }
                  }
            }

            // Boundaries
            double xw = (x[1]-x[0])/2;
            double yw = (y[1]-y[0])/2;

            gr[k] = new TH2D(("h2_" + std::to_string(k)).c_str(),"title",
                                   BINS[0], x[0]-xw, x[BINS[0]-1]+xw,
                                   BINS[1], y[0]-yw, y[BINS[1]-1]+yw);


            const std::string titletxt = (DATAMODE == "E_lm") ? "MC" : param.TYPE;
    // Efficiency
    always MC
            if (!(l == 0 && m == 0)) {
                  gr[k]->SetTitle(Form("%s: #it{lm} = <%d,%d>", titletxt.c_str(), l,
    m));
            } else {
                  gr[k]->SetTitle(Form("%s: %s: #it{lm} = <%d,%d>", titletxt.c_str(),
    outputfile.c_str(), l, m));
            }

            // Set x-axis label
            gr[k]->GetXaxis()->SetTitle(label[0].c_str());
            gr[k]->GetYaxis()->SetTitle(label[1].c_str());

            gr[k]->GetXaxis()->SetTitleSize(0.05);
            gr[k]->GetXaxis()->SetLabelSize(0.05);

            gr[k]->GetYaxis()->SetTitleSize(0.05);
            gr[k]->GetYaxis()->SetLabelSize(0.05);

            gStyle->SetTitleFontSize(0.08);

            // Plot
            for (std::size_t i = 0; i < BINS[0]; ++i) {
                   for (std::size_t j = 0; j < BINS[1]; ++j) {
                          gr[k]->Fill(x[i], y[j], z[i][j]);
                   }
            }
            gr[k]->Draw("COLZ");

            //gr[k]->Draw("A B");

            // gr[k]->Draw("ALP");
            // c1->Update();
            ++k;
          }
    }

    // Require that we have truly 2D data
    if (BINS[0] > 1 && BINS[1] > 1) {
    aux::CreateDirectory("./figs");
    aux::CreateDirectory("./figs/harmonicfit");
    aux::CreateDirectory("./figs/harmonicfit/" + outputpath);
    const std::string subpath = "2D_OBS_" + std::to_string(OBSERVABLE[0]) +
    "_OBS_" +
    std::to_string(OBSERVABLE[1]) + "_" + legendstr;
    aux::CreateDirectory("./figs/harmonicfit/" + outputpath + "/" + subpath);
    c1->Print(Form("./figs/harmonicfit/%s/%s/%s.pdf", outputpath.c_str(),
    subpath.c_str(),
    outputfile.c_str()));
    }

    for (std::size_t i = 0; i < gr.size(); ++i) {
          delete gr[i];
    }
    delete c1;
    */
}

// Loop over system mass, pt, rapidity
//
void MHarmonic::HyperLoop(void (*fitfunc)(int &, double *, double &, double *, int),
                          const std::vector<gra::spherical::Omega> &MC,
                          const std::vector<gra::spherical::Data> &DATA, const HPARAM &hp) {
  // Initialize detector expansion tensors
  fla_DET = MTensor<gra::spherical::SH_DET>(
      {(unsigned int)hp.M[0], (unsigned int)hp.PT[0], (unsigned int)hp.Y[0]});
  fid_DET = fla_DET;
  det_DET = fla_DET;

  // ------------------------------------------------------------------
  // Create grid discretization

  grid.resize(3);
  grid[0].resize(hp.M[0]);
  grid[1].resize(hp.PT[0]);
  grid[2].resize(hp.Y[0]);

  // Mass steps
  const double M_STEP  = (hp.M[2] - hp.M[1]) / hp.M[0];
  const double PT_STEP = (hp.PT[2] - hp.PT[1]) / hp.PT[0];
  const double Y_STEP  = (hp.Y[2] - hp.Y[1]) / hp.Y[0];

  for (std::size_t i = 0; i < hp.M[0]; ++i) {
    grid[0][i].min = i * M_STEP + hp.M[1];
    grid[0][i].max = grid[0][i].min + M_STEP;
  }
  for (std::size_t j = 0; j < hp.PT[0]; ++j) {
    grid[1][j].min = j * PT_STEP + hp.PT[1];
    grid[1][j].max = grid[1][j].min + PT_STEP;
  }
  for (std::size_t k = 0; k < hp.Y[0]; ++k) {
    grid[2][k].min = k * Y_STEP + hp.Y[1];
    grid[2][k].max = grid[2][k].min + Y_STEP;
  }

  // ------------------------------------------------------------------
  // Expand the detector transfer function

  for (const auto &i : indices(grid[0])) {
    for (const auto &j : indices(grid[1])) {
      for (const auto &k : indices(grid[2])) {
        const std::vector<std::size_t> MC_ind = gra::spherical::GetIndices(
            MC, {grid[0][i].min, grid[0][i].max}, {grid[1][j].min, grid[1][j].max},
            {grid[2][k].min, grid[2][k].max});

        // Acceptance mixing matrices
        fla_DET({i, j, k}).MIXlm = gra::spherical::GetGMixing(MC, MC_ind, param.LMAX, "fla");
        fid_DET({i, j, k}).MIXlm = gra::spherical::GetGMixing(MC, MC_ind, param.LMAX, "fid");
        det_DET({i, j, k}).MIXlm = gra::spherical::GetGMixing(MC, MC_ind, param.LMAX, "det");

        // Get efficiency decomposition for this interval
        std::pair<std::vector<double>, std::vector<double>> E0 =
            gra::spherical::GetELM(MC, MC_ind, param.LMAX, "fla");
        std::pair<std::vector<double>, std::vector<double>> E1 =
            gra::spherical::GetELM(MC, MC_ind, param.LMAX, "fid");
        std::pair<std::vector<double>, std::vector<double>> E2 =
            gra::spherical::GetELM(MC, MC_ind, param.LMAX, "det");

        fla_DET({i, j, k}).E_lm       = E0.first;
        fla_DET({i, j, k}).E_lm_error = E0.second;

        fid_DET({i, j, k}).E_lm       = E1.first;
        fid_DET({i, j, k}).E_lm_error = E1.second;

        det_DET({i, j, k}).E_lm       = E2.first;
        det_DET({i, j, k}).E_lm_error = E2.second;
      }
    }
  }

  // ------------------------------------------------------------------

  // Loop over data sources
  for (const auto &ind : indices(DATA)) {
    // --------------------------------------------------------
    // Pre-Calculate once Spherical Harmonics for the MINUIT fit
    DATA_events = DATA[ind].EVENTS;
    Y_lm        = gra::spherical::YLM(DATA_events, param.LMAX);
    // --------------------------------------------------------

    double chi2 = 0.0;

    // Data source identifier string
    const gra::spherical::Meta META = DATA[ind].META;

    // Initialize tensor arrays
    MTensor<gra::spherical::SH> temp(
        {(unsigned int)hp.M[0], (unsigned int)hp.PT[0], (unsigned int)hp.Y[0]});
    fla[META] = temp;
    fid[META] = temp;
    det[META] = temp;

    // Expand Data
    for (const auto &i : indices(grid[0])) {
      for (const auto &j : indices(grid[1])) {
        for (const auto &k : indices(grid[2])) {
          // Data indices
          DATA_ind = gra::spherical::GetIndices(DATA_events, {grid[0][i].min, grid[0][i].max},
                                                {grid[1][j].min, grid[1][j].max},
                                                {grid[2][k].min, grid[2][k].max});

          const unsigned int MINEVENTS = 75;
          if (DATA_ind.size() < MINEVENTS) {
            std::cout << rang::fg::red << "WARNING: Less than " << MINEVENTS << " in the cell!"
                      << rang::fg::reset << std::endl;
          }

          // ==============================================================
          // ALGORITHM 1: DIRECT / OBSERVED / ALGEBRAIC decomposition

          fid[META]({i, j, k}).t_lm_MPP =
              gra::spherical::SphericalMoments(DATA_events, DATA_ind, param.LMAX, "fid");
          det[META]({i, j, k}).t_lm_MPP =
              gra::spherical::SphericalMoments(DATA_events, DATA_ind, param.LMAX, "det");

          // Forward matrix
          Eigen::MatrixXd M = gra::aux::Matrix2Eigen(det_DET({i, j, k}).MIXlm);

          // Matrix pseudoinverse
          Eigen::VectorXd b = gra::aux::Vector2Eigen(det[META]({i, j, k}).t_lm_MPP);
          Eigen::VectorXd x = gra::math::PseudoInverse(M, param.SVDREG) * b;

          // Collect the inversion result
          fla[META]({i, j, k}).t_lm_MPP = gra::aux::Eigen2Vector(x);

          // ==============================================================

          // Update these to proper errors (TBD. PUT NAIVE POISSON
          // COUNTING FOR NOW)
          fla[META]({i, j, k}).t_lm_MPP_error =
              std::vector<double>(fla[META]({i, j, k}).t_lm_MPP.size(), 0.0);
          for (const auto &z : indices(fla[META]({i, j, k}).t_lm_MPP)) {
            fla[META]({i, j, k}).t_lm_MPP_error[z] =
                msqrt(std::abs(fla[META]({i, j, k}).t_lm_MPP[z]));
          }

          // Propagate errors assuming no error on mixing matrix
          // (infinite reference MC statistics limit)
          fid[META]({i, j, k}).t_lm_MPP_error =
              spherical::ErrorProp(fid_DET({i, j, k}).MIXlm, fla[META]({i, j, k}).t_lm_MPP_error);
          det[META]({i, j, k}).t_lm_MPP_error =
              spherical::ErrorProp(det_DET({i, j, k}).MIXlm, fla[META]({i, j, k}).t_lm_MPP_error);

          std::cout << "Algebraic Moore-Penrose/SVD inverted "
                       "(unmixed) moments in the (angular flat) "
                       "reference phase space:"
                    << std::endl;
          gra::spherical::PrintOutMoments(fla[META]({i, j, k}).t_lm_MPP,
                                          fla[META]({i, j, k}).t_lm_MPP_error, ACTIVE, param.LMAX);

          std::cout << "Algebraic (mixed) moments in the fiducial "
                       "phase space:"
                    << std::endl;
          gra::spherical::PrintOutMoments(fid[META]({i, j, k}).t_lm_MPP,
                                          fid[META]({i, j, k}).t_lm_MPP_error, ACTIVE, param.LMAX);

          std::cout << "Algebraic (mixed) moments in the detector space:" << std::endl;
          gra::spherical::PrintOutMoments(det[META]({i, j, k}).t_lm_MPP,
                                          det[META]({i, j, k}).t_lm_MPP_error, ACTIVE, param.LMAX);

          // ==============================================================
          // ALGORITHM 2: Extended Maximum Likelihood Fit

          if (param.EML) {
            // *** Get TMinuit based fit decomposition ***
            MomentFit(META, {i, j, k}, fitfunc);

            // Save result
            fla[META]({i, j, k}).t_lm_EML = t_lm;
            fid[META]({i, j, k}).t_lm_EML =
                fid_DET({i, j, k}).MIXlm * t_lm;  // Moment forward rotation: Matrix *
            // Vector
            det[META]({i, j, k}).t_lm_EML =
                det_DET({i, j, k}).MIXlm * t_lm;  // Moment forward rotation: Matrix *
            // Vector

            // Uncertanties
            fla[META]({i, j, k}).t_lm_EML_error = t_lm_error;

            // Propagate errors assuming no error on mixing
            // matrix (infinite reference MC statistics limit)
            fid[META]({i, j, k}).t_lm_EML_error =
                spherical::ErrorProp(fid_DET({i, j, k}).MIXlm, t_lm_error);
            det[META]({i, j, k}).t_lm_EML_error =
                spherical::ErrorProp(det_DET({i, j, k}).MIXlm, t_lm_error);
            // ==============================================================

            std::cout << "Extended Maximum-Likelihood inverse "
                         "fitted (unmixed) moments in the "
                         "(angular flat)  phase space:"
                      << std::endl;
            gra::spherical::PrintOutMoments(fla[META]({i, j, k}).t_lm_EML,
                                            fla[META]({i, j, k}).t_lm_EML_error, ACTIVE,
                                            param.LMAX);

            std::cout << "Extended Maximum-Likelihood "
                         "Re-Back-Projected (mixed) moments in "
                         "the fiducial phase space:"
                      << std::endl;
            gra::spherical::PrintOutMoments(fid[META]({i, j, k}).t_lm_EML,
                                            fid[META]({i, j, k}).t_lm_EML_error, ACTIVE,
                                            param.LMAX);

            std::cout << "Extended Maximum-Likelihood "
                         "Re-Back-Projected (mixed) moments in "
                         "the detector space:"
                      << std::endl;
            gra::spherical::PrintOutMoments(det[META]({i, j, k}).t_lm_EML,
                                            det[META]({i, j, k}).t_lm_EML_error, ACTIVE,
                                            param.LMAX);
          }

          // --------------------------------------------------------------
          // Print results
          chi2 += PrintOutHyperCell(META, {i, j, k});

          // --------------------------------------------------------------
          // Make comparison of synthetic MC data vs estimate

          if (DATA[ind].META.MODE == "MC") {
            int fiducial = 0;
            int selected = 0;
            for (const auto &l : DATA_ind) {
              if (DATA_events[l].fiducial) { ++fiducial; }
              if (DATA_events[l].fiducial && DATA_events[l].selected) { ++selected; }
            }
            std::cout << std::endl;

            std::cout << rang::fg::yellow << "MC GROUND TRUTH: " << rang::fg::reset << std::endl;
            printf(
                "   MC 'synthetic data' events generated       "
                "      = %u \n",
                (unsigned int)DATA_ind.size());
            printf(
                "   MC 'synthetic data' events fiducial        "
                "      = %d (acceptance %0.1f percent) \n",
                fiducial, fiducial / (double)DATA_ind.size() * 100);
            printf(
                "   MC 'synthetic data' events fiducial and "
                "selected = %d (efficiency %0.1f percent) \n",
                selected, selected / (double)fiducial * 100);
            std::cout << std::endl;
          }
        }
      }
    }

    if (param.EML) {
      gra::aux::PrintBar("=");
      double reducedchi2 = chi2 / (double)(ACTIVENDF * hp.M[0] * hp.PT[0] * hp.Y[0]);
      if (reducedchi2 < 3) {
        std::cout << rang::fg::green;
      } else {
        std::cout << rang::fg::red;
      }
      printf(
          "Total chi2(MPP - EML) / (ACTIVENDF x BINS) = %0.2f / (%d x %d) = "
          "%0.2f \n",
          chi2, ACTIVENDF, (int)(hp.M[0] * hp.PT[0] * hp.Y[0]), reducedchi2);

      std::cout << rang::fg::reset;
      gra::aux::PrintBar("=");
      std::cout << std::endl << std::endl;
    }

  }  // Loop over data sources
}

// Print out results for the hypercell
//
double MHarmonic::PrintOutHyperCell(const gra::spherical::Meta &    META,
                                    const std::vector<std::size_t> &cell) {
  double chi2 = 0;

  // Print information
  META.Print();

  // Extended Maximum Likelihood
  if (param.EML == true) {
    // Loop over moments
    for (int l = 0; l <= param.LMAX; ++l) {
      for (int m = -l; m <= l; ++m) {
        const int index = gra::spherical::LinearInd(l, m);

        const double obs = det[META](cell).t_lm_MPP[index];
        const double fit = det[META](cell).t_lm_EML[index];

        // Moment active
        if (ACTIVE[index]) { chi2 += pow2(obs - fit) / pow2(obs); }
      }
    }
    std::cout << std::endl;

    const double reducedchi2 = chi2 / (double)ACTIVENDF;
    if (reducedchi2 < 3) {
      std::cout << rang::fg::green;
    } else {
      std::cout << rang::fg::red;
    }
    printf("chi2(MPP - EML) / ndf = %0.3f / %d = %0.3f \n", chi2, ACTIVENDF, reducedchi2);
    std::cout << rang::fg::reset << std::endl;

    const double sum_fla = fla[META](cell).t_lm_EML[gra::spherical::LinearInd(0, 0)];
    printf(
        "EML: Estimate of events in this hyperbin in the flat phase space = "
        "%0.1f +- %0.1f "
        "\n",
        sum_fla, msqrt(sum_fla)); // Poisson error

    const double sum_FID = gra::spherical::HarmDotProd(fid_DET(cell).E_lm, fla[META](cell).t_lm_EML,
                                                       ACTIVE, param.LMAX);
    printf(
        "EML: Estimate of events in this hyperbin in the fiducial  phase "
        "space = %0.1f "
        "+- %0.1f \n",
        sum_FID, msqrt(sum_FID)); // Poisson error

    const double sum_DET = gra::spherical::HarmDotProd(det_DET(cell).E_lm, fla[META](cell).t_lm_EML,
                                                       ACTIVE, param.LMAX);
    printf(
        "EML: Estimate of events in this hyperbin in the detector        "
        "space = %0.1f "
        "+- %0.1f \n",
        sum_DET, msqrt(sum_DET)); // Poisson error
  }
  std::cout << std::endl;

  // Algebraic inverse
  const double sum_fla = fla[META](cell).t_lm_MPP[gra::spherical::LinearInd(0, 0)];
  printf(
      "MPP: Estimate of events in this hyperbin in the flat phase space = %0.1f "
      "+- %0.1f \n",
      sum_fla, msqrt(sum_fla)); // Poisson error

  const double sum_FID =
      gra::spherical::HarmDotProd(fid_DET(cell).E_lm, fla[META](cell).t_lm_MPP, ACTIVE, param.LMAX);
  printf(
      "MPP: Estimate of events in this hyperbin in the fiducial  phase "
      "space = %0.1f +- "
      "%0.1f \n",
      sum_FID, msqrt(sum_FID)); // Poisson error

  const double sum_DET =
      gra::spherical::HarmDotProd(det_DET(cell).E_lm, fla[META](cell).t_lm_MPP, ACTIVE, param.LMAX);
  printf(
      "MPP: Estimate of events in this hyperbin in the detector        "
      "space = %0.1f +- "
      "%0.1f \n",
      sum_DET, msqrt(sum_DET)); // Poisson error

  return chi2;
}

// MINUIT based fit routine for the Extended Maximum Likelihood formalism
//
void MHarmonic::MomentFit(const gra::spherical::Meta &META, const std::vector<std::size_t> &cell,
                          void (*fitfunc)(int &, double *, double &, double *, int)) {
  std::cout << "MomentFit: Starting ..." << std::endl << std::endl;

  // **** This must be set for the loss function ****
  activecell = cell;

  // Init TMinuit
  TMinuit *gMinuit = new TMinuit(NCOEF);  // initialize TMinuit with a maximum of N params
  gMinuit->SetFCN(fitfunc);

  // Set Print Level
  // -1 no output
  // 1 standard output
  gMinuit->SetPrintLevel(-1);

  // Set error Definition
  // 1 for Chi square
  // 0.5 for negative log likelihood
  //    gMinuit->SetErrorDef(0.5);
  double arglist[2];
  int    ierflg = 0;
  arglist[0]    = 0.5;  // 0.5 <=> We use negative log likelihood cost function
  gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);

  // double fminbest = 1e32;
  // Try different initial values to find out true minimum
  const std::size_t TRIALMAX = 1;

  for (std::size_t trials = 0; trials < TRIALMAX; ++trials) {
    for (int l = 0; l <= param.LMAX; ++l) {
      for (int m = -l; m <= l; ++m) {
        const int index = gra::spherical::LinearInd(l, m);

        const std::string str = "t_" + std::to_string(l) + std::to_string(m);

        const double max = 1e9;  // Physically bounded ultimately by the number of events
        const double min = -1e5;

        // ======================================================
        // *** Use the algebraic inverse solution as a good starting value
        // ***
        const double start_value = fla[META](activecell).t_lm_MPP[index];
        // ======================================================

        const double step_value = 0.1;  // in units of Events

        gMinuit->mnparm(index, str, start_value, step_value, min, max, ierflg);

        // After first trial, fix t_00 (error are not estimated with
        // constant parameters)
        // if (trials > 0) {
        //	gMinuit->mnparm(0, "t_00", t_lm[0], 0, 0, 0, ierflg);
        //	gMinuit->FixParameter(0);
        //}
        // FIX ODD MOMENTS TO ZERO
        if (param.REMOVEODD && ((l % 2) != 0)) {
          gMinuit->mnparm(index, str, 0, 0, 0, 0, ierflg);
          gMinuit->FixParameter(index);
        }
        // FIX NEGATIVE M TO ZERO
        if (param.REMOVENEGATIVEM && (m < 0)) {
          gMinuit->mnparm(index, str, 0, 0, 0, 0, ierflg);
          gMinuit->FixParameter(index);
        }
      }
    }
    // Scan main parameter
    // gMinuit->mnscan();
    arglist[0] = 100000;   // Minimum number of function calls
    arglist[1] = 0.00001;  // Minimum tolerance

    // First simplex to find approximate answers
    gMinuit->mnexcm("SIMPLEX", arglist, 2, ierflg);

    // Numerical Hessian (2nd derivatives matrix), inverse of this -> covariance
    // matrix
    gMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);

    // Confidence intervals based on the profile likelihood ratio
    gMinuit->mnexcm("MINOS", arglist, 2, ierflg);

    // Calculate error & covariance matrix
    gMinuit->mnemat(&errmat[0][0], NCOEF);

    for (int i = 0; i < NCOEF; ++i) {
      for (int j = 0; j < NCOEF; ++j) {
        covmat[i][j] = sqrt(errmat[i][i] * errmat[j][j]);
        if (covmat[i][j] > 1E-80) {
          covmat[i][j] = errmat[i][j] / covmat[i][j];
        } else
          covmat[i][j] = 0.0;
      }
    }
    errmat.Print("Error matrix");
    covmat.Print("Covariance matrix");

    // Print results
    double fmin, fedm, errdef  = 0.0;
    int    nvpar, nparx, istat = 0;
    gMinuit->mnstat(fmin, fedm, errdef, nvpar, nparx, istat);

    // Collect fit result
    for (int l = 0; l <= param.LMAX; ++l) {
      for (int m = -l; m <= l; ++m) {
        const int index = gra::spherical::LinearInd(l, m);

        // Set results into t_lm[], t_lm_error[]
        gMinuit->GetParameter(index, t_lm[index], t_lm_error[index]);
      }
    }

  }  // TRIALS LOOP END

  // Sum of N random variables Y = X_1 + X_2 + ... X_N
  // sigma_Y^2 = \sum_i^N \sigma_i^2 + 2 \sum_i \sum_i < j
  // covariance(X_i, X_j)

  // Print out results (see MINUIT manual for these parameters)
  double fmin, fedm, errdef  = 0.0;
  int    nvpar, nparx, istat = 0;
  gMinuit->mnstat(fmin, fedm, errdef, nvpar, nparx, istat);
  gMinuit->mnprin(4, fmin);

  std::cout << "<MINUIT (MIGRAD+MINOS)> done." << std::endl << std::endl;

  delete gMinuit;
}

// Unbinned Extended Maximum Likelihood function
// Extended means that the number of events (itself) is a Poisson distributed
// random variable and that is incorporated to the fit.
void MHarmonic::logLfunc(int &npar, double *gin, double &f, double *par, int iflag) const {
  // Collect fit t_LM coefficients from MINUIT
  std::vector<double> T(NCOEF, 0.0);

  for (const auto &i : indices(T)) {
    T[i] = par[i];
    // printf("T[%d] = %0.5f \n", i);
  }
  // This is the number of events in the current phase space point at detector
  // level
  const int METHOD = 2;
  double    nhat   = 0.0;

  // Equivalent estimator 1
  if (METHOD == 1) {
    const std::vector<double> t_lm_det = det_DET(activecell).MIXlm * T;  // Matrix * Vector
    nhat                               = t_lm_det[0];
  }
  // Equivalent estimator 2
  if (METHOD == 2) {
    nhat = gra::spherical::HarmDotProd(det_DET(activecell).E_lm, T, ACTIVE, param.LMAX);
  }

  // For each event, calculate \sum_{LM} t_{lm}
  // Re[Y_{lm}(costheta,phi_k)], k is the event index
  std::vector<double> I0;
  const double        V = msqrt(4.0 * PI);  // Normalization volume

  for (const auto &k : DATA_ind) {
    // Event is accepted
    if (DATA_events[k].fiducial && DATA_events[k].selected) {
      // fine
    } else {
      continue;
    }

    // Loop over (l,m) terms
    double sum = 0.0;
    for (int l = 0; l <= param.LMAX; ++l) {
      for (int m = -l; m <= l; ++m) {
        const int index = gra::spherical::LinearInd(l, m);

        if (ACTIVE[index]) {
          // Calculate here
          // const std::complex<double> Y =
          //	gra::math::Y_complex_basis(DATA_events[k].costheta,DATA_events[k].phi,
          // l,m);
          // const double ReY = gra::math::NReY(Y,l,m);

          // Pre-calculated for speed
          const double ReY = Y_lm[k][index];

          // Add
          sum += T[index] * ReY;
        }
      }
    }
    I0.push_back(V * sum);
  }
  // Fidelity term
  double logL = 0;

  for (const auto &k : indices(I0)) {
    if (I0[k] > 0) {
      logL += std::log(I0[k]);
      // Positivity is not satisfied, make strong penalty
    } else {
      logL = -1e32;
      // The Re[Y_lm] (purely real) formulation here does
      // not, explicitly, enforce non-negativity.
      // Where as |A|^2 formulation would enforce it by
      // construction,
      // but that would have -l <= m <= l parameters, i.e.,
      // over redundant representation problem
      // on the other hand.
    }
  }

  // Note the minus sign on logL term
  f = -logL + nhat;

  // High penalty, total event count estimate has gone out of physical
  if (nhat < 0) { f = 1e32; }

  // L1-norm regularization (Laplace prior), -> + ln(P_laplace)
  double             l1term = 0.0;
  const unsigned int START  = 1;
  for (std::size_t i = START; i < T.size(); ++i) { l1term += std::abs(T[i]); }
  f += l1term * param.L1REG * nhat;  // nhat for scale normalization

  // printf("MHarmonic:: cost-functional = %0.4E, nhat = %0.1f \n", f, nhat);
}

}  // gra namespace ends
