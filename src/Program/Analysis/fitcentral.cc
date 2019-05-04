// GRANIITTI - Monte Carlo event generator for high energy diffraction
// https://github.com/mieskolainen/graniitti
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <math.h>
#include <algorithm>
#include <chrono>
#include <complex>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <memory>
#include <mutex>
#include <random>
#include <stdexcept>
#include <thread>
#include <vector>

// Own
#include "Graniitti/MAux.h"
#include "Graniitti/MGraniitti.h"
#include "Graniitti/MH1.h"
#include "Graniitti/MMath.h"
#include "Graniitti/MRegge.h"

// ROOT
#include "TMinuit.h"
#include "TString.h"

// Libraries
#include "json.hpp"
#include "cxxopts.hpp"
#include "rang.hpp"

using gra::aux::indices;

using namespace gra;

namespace fitcentral {
// Single dataset / observable
struct DATASET {
  std::vector<double> data;  // Events

  std::string INPUTFILE;
  std::string DATAFILE;
  int         DATATYPE = 0;

  int    OBSERVABLE = 0;
  uint   NBINS      = 0;
  double MMIN       = 0.0;
  double MMAX       = 0.0;

  MH1<double> *h1DATA;
};

// All datasets here
std::vector<DATASET> datasets;

// Function prototypes
void RunMC(double *par);
bool ReadData(DATASET &dataset);
void Chi2Func(int &npar, double *gin, double &f, double *par, int iflag);
void SetSoftParam(std::map<std::string, gra::PARAM_RES> &RESONANCES, double *par);

// Choose continuum form factor
int  offshellFF = 1;
bool POMLOOP    = false;

// Resonance fixlist and zerolist
std::vector<std::string> afixlist;
std::vector<std::string> pfixlist;
std::vector<std::string> zerolist;

// ===============================================================
// Note that TMinuit uses (internally) Fortran indexing =>
// The parameter 0 in C++ is the parameter 1 in TMinuit output.
//
// Main function

// Set soft parameters back to the model
void SetSoftParam(std::map<std::string, gra::PARAM_RES> &RESONANCES, double *par) {
  printf("fitcentral::SetSoftParam: \n\n");

  PARAM_REGGE::offshellFF = offshellFF;
  int k                   = -1;

  // Loop over resonances
  for (const auto &x : RESONANCES) {
    double A              = par[++k];
    double phi            = par[++k];
    RESONANCES[x.first].g = A * std::exp(gra::math::zi * phi);

    printf("[%s] \t A = %0.8f \t phi = %0.6f (%0.0f deg) ", x.first.c_str(), A, phi,
           phi / (2 * gra::math::PI) * 360);

    bool afix = std::find(afixlist.begin(), afixlist.end(), x.first) != afixlist.end();
    bool pfix = std::find(pfixlist.begin(), pfixlist.end(), x.first) != pfixlist.end();

    if (afix && pfix) {
      std::cout << rang::fg::red << "\t[AMPLITUDE & PHASE FIXED]" << rang::fg::reset << std::endl;
    } else if (afix) {
      std::cout << rang::fg::green << "\t[AMPLITUDE FIXED]" << rang::fg::reset << std::endl;
    } else if (pfix) {
      std::cout << rang::fg::blue << "\t[PHASE FIXED]" << rang::fg::reset << std::endl;
    } else if (std::find(zerolist.begin(), zerolist.end(), x.first) != zerolist.end()) {
      std::cout << rang::fg::yellow << "\t[SET TO ZERO]" << rang::fg::reset << std::endl;
    } else {
      std::cout << std::endl;
    }
  }
  std::cout << std::endl;

  // Continuum form factors
  if (offshellFF == 0) {
    PARAM_REGGE::b_EXP = par[++k];
    printf("[Offshell-FF] \t b_EXP = %0.3f", PARAM_REGGE::b_EXP);
  }
  if (offshellFF == 1) {
    PARAM_REGGE::a_OREAR = par[++k];
    PARAM_REGGE::b_OREAR = par[++k];
    printf("[Offshell-FF] \t a/b_OREAR = [%0.3f,%0.3f]", PARAM_REGGE::a_OREAR,
           PARAM_REGGE::b_OREAR);
  }
  if (offshellFF == 2) {
    PARAM_REGGE::b_POW = par[++k];
    printf("[Offshell-FF] \t b_POW = %0.3f", PARAM_REGGE::b_POW);
  }

  bool fixcontinuum =
      (std::find(afixlist.begin(), afixlist.end(), "continuum") != afixlist.end()) ||
      (std::find(pfixlist.begin(), pfixlist.end(), "continuum") != pfixlist.end());

  if (fixcontinuum) {
    std::cout << rang::fg::red << "\t[FIXED]" << rang::fg::reset << std::endl;
  } else {
    std::cout << std::endl;
  }
}

// Read input data event by event
bool ReadData(DATASET &dataset) {
  printf("fitcentral::ReadData: with inputfile %s \n", dataset.DATAFILE.c_str());

  FILE * fp;
  double M  = 0.0;
  double Px = 0.0;
  double Py = 0.0;
  double Pz = 0.0;

  if ((fp = fopen(dataset.DATAFILE.c_str(), "r+")) == NULL) {
    printf("No inputfile %s found \n", dataset.DATAFILE.c_str());
    return false;
  }

  while (true) {
    // System: mass,Px,Py,Pz
    int ret = -1;

    if (dataset.DATATYPE == 4) ret = fscanf(fp, "%lf,%lf,%lf,%lf", &M, &Px, &Py, &Pz);

    if (dataset.DATATYPE == 1) ret = fscanf(fp, "%lf", &M);

    if (ret == dataset.DATATYPE) {  // number of values per line
                                    // Fine
    } else if (ret == EOF) {
      break;
    } else {
      printf("Error in the file structure of %s!\n", dataset.DATAFILE.c_str());
      return false;
    }

    // Fill out vector
    dataset.data.push_back(M);

    // DEBUG print
    // printf("x = %0.8f \n", xval);

    // Fill out mass histogram
    dataset.h1DATA->Fill(M);
  }
  fclose(fp);

  dataset.h1DATA->Print();

  return true;
}

MEikonal eik;  // For speed, calculate only once
int      iter = 0;

// Cost function
void Chi2Func(int &npar, double *gin, double &f, double *par, int iflag) {
  std::cout << rang::fg::red << "fitcentral::Chi2Func: iter = " << iter << rang::fg::reset
            << std::endl;

  double cost   = 0;
  double chi2   = 0;
  int    points = 0;

  // Loop over different datasets
  for (const auto &k : indices(datasets)) {
    double local_cost   = 0.0;
    double local_chi2   = 0.0;
    int    local_points = 0;

    // Generator
    std::unique_ptr<MGraniitti> gen = std::make_unique<MGraniitti>();

    try {
      if (iter > 1) {
        gen->HILJAA = true;  // SILENT output
      }

      // Read process input from file
      gen->ReadInput(datasets[k].INPUTFILE);
      std::map<std::string, gra::PARAM_RES> RESONANCES = gen->proc->GetResonances();

      // Set new parameters
      fitcentral::SetSoftParam(RESONANCES, par);
      gen->proc->SetResonances(RESONANCES);

      // Set screening
      gen->proc->SetScreening(fitcentral::POMLOOP);

      // Force histogramming level 1
      gen->proc->SetHistograms(1);

      // Reset histogram bounds according to what we use here
      // (MUST BE DONE BEFORE INITIALIZATION!)
      gen->proc->h1["M"].ResetBounds(datasets[k].NBINS, datasets[k].MMIN, datasets[k].MMAX);

      // Initialize (always last)
      if (iter == 0) {
        gen->Initialize();
        eik = gen->proc->GetEikonal();
      } else {
        gen->Initialize(eik);
      }
      ++iter;

      // Generate events
      // gen->Generate();

      // Histogram fusion due to multithreading!
      gen->HistogramFusion();

      printf("MC: \n");
      gen->proc->h1["M"].Print();
      // ---------------------------------------------------------

      // Get normalization from data
      double S = datasets[k].h1DATA->SumWeights() / gen->proc->h1["M"].SumWeights();

      // Chi^2 and log-likelihood
      for (std::size_t i = 0; i < datasets[k].h1DATA->XBins(); ++i) {
        if (datasets[k].h1DATA->GetBinWeight(i) > 0) {
          // Chi^2
          local_chi2 += gra::math::pow2(S * gen->proc->h1["M"].GetBinWeight(i) -
                                        datasets[k].h1DATA->GetBinWeight(i)) /
                        datasets[k].h1DATA->GetBinWeight(i);

          // Log-likelihood
          double nb = datasets[k].h1DATA->GetBinWeight(i);
          double fb = S * gen->proc->h1["M"].GetBinWeight(i);
          if (fb > 0) {
            // make it negative (we use minimizer)
            local_cost -= nb * std::log(fb);
          }
          ++local_points;
        }
      }

      // Make it compatible with chi^2 uncertainty calculation
      local_cost *= 2.0;

      // UPDATE GLOBAL
      chi2 += local_chi2;
      cost += local_cost;
      points += local_points;

      printf("DATA: \n");
      datasets[k].h1DATA->Print();

      if (chi2 / (double)points < 2.0) {
        std::cout << rang::fg::green;
      } else {
        std::cout << rang::fg::red;
      }
      printf(
          "============================================================="
          " \n");
      printf("DATASET %lu/%lu \n", k + 1, datasets.size());
      printf("LOCAL \n");
      printf(" chi^2 / points = %0.2f / %d = %0.2f \n", local_chi2, local_points,
             local_chi2 / (double)local_points);
      printf(" -logL = %0.2f \n", local_cost);
      printf("GLOBAL \n");
      printf(" chi^2 / points = %0.2f / %d = %0.2f \n", chi2, points, chi2 / (double)points);
      printf(" -logL = %0.2f \n", cost);
      printf(
          "============================================================="
          " \n\n");
      std::cout << rang::fg::reset;
    } catch (const std::invalid_argument &e) {
      gra::aux::PrintGameOver();
      std::cerr << rang::fg::red << "Exception catched: " << rang::fg::reset << e.what()
                << std::endl;
      exit(0);
    } catch (const std::ios_base::failure &e) {
      gra::aux::PrintGameOver();
      std::cerr << rang::fg::red << "Exception catched: std::ios_base::failure: " << rang::fg::reset
                << e.what() << std::endl;
      exit(0);
    } catch (...) {
      gra::aux::PrintGameOver();
      std::cerr << rang::fg::red << "Exception catched: Unspecified (...) (Probably JSON input)"
                << rang::fg::reset << std::endl;
      exit(0);
    }

  }  // Over different datasets

  f = cost;
}
}  // Namespace fitcentral ends

int main(int argc, char *argv[]) {
  gra::aux::PrintFlashScreen(rang::fg::red);
  std::cout << rang::style::bold << "GRANIITTI - Soft Central Production Fitter"
            << rang::style::reset << std::endl
            << std::endl;
  gra::aux::PrintVersion();

  // Save the number of input arguments
  const int NARGC = argc - 1;

  try {
    cxxopts::Options options(argv[0], "");
    options.add_options()("i,input", "Input dataset collection numbers               <0,1,2,3,...>",
                          cxxopts::value<std::string>())(
        "c,continuum", "Continuum form factor                          <0|1|2>",
        cxxopts::value<int>())(
        "a,afix", "Fix resonance amplitude to their input values  <f0_980,f2_1270,...>",
        cxxopts::value<std::string>())(
        "p,pfix", "Fix resonance phase to their input values      <f0_980,f2_1270,...>",
        cxxopts::value<std::string>())(
        "x,fix", "Fix resonance amplitude and phase              <f0_980,f2_1270,...>",
        cxxopts::value<std::string>())(
        "z,zero", "Fix resonances to zero                         <f2_1525,...>",
        cxxopts::value<std::string>())("A,allfixed",
                                       "Fix all resonance parameters               <true|false>",
                                       cxxopts::value<std::string>())(
        "l,POMLOOP", "Screening Pomeron loop", cxxopts::value<std::string>())("H,help", "Help");
    auto r = options.parse(argc, argv);

    if (r.count("help") || NARGC == 0) {
      std::cout << options.help({""}) << std::endl;
      std::cout << rang::style::bold << "Example:" << rang::style::reset << std::endl;
      std::cout << "  " << argv[0] << " -i 0,1 -c 2 -x 'rho_770,continuum' " << std::endl
                << std::endl;
      return EXIT_FAILURE;
    }

    // Input dataset numbers
    std::vector<int> input = gra::aux::SplitStr2Int(r["input"].as<std::string>());

    // Off-shell form factor
    fitcentral::offshellFF = r["c"].as<int>();

    // Amplitude
    if (r.count("a")) { fitcentral::afixlist = gra::aux::SplitStr2Str(r["a"].as<std::string>()); }
    // Phase
    if (r.count("p")) { fitcentral::pfixlist = gra::aux::SplitStr2Str(r["p"].as<std::string>()); }
    // Both
    if (r.count("x")) {
      std::vector<std::string> out1 = gra::aux::SplitStr2Str(r["x"].as<std::string>());
      fitcentral::afixlist.insert(fitcentral::afixlist.end(), out1.begin(), out1.end());

      std::vector<std::string> out2 = gra::aux::SplitStr2Str(r["x"].as<std::string>());
      fitcentral::pfixlist.insert(fitcentral::pfixlist.end(), out2.begin(), out2.end());
    }
    // Zero
    if (r.count("z")) { fitcentral::zerolist = gra::aux::SplitStr2Str(r["z"].as<std::string>()); }

    // Fix-all resonances
    bool FIXALLRES = false;
    if (r.count("A")) { FIXALLRES = (r["A"].as<std::string>() == "true" ? true : false); }

    // Pomeron loop
    if (r.count("l")) {
      const std::string val = r["l"].as<std::string>();
      fitcentral::POMLOOP   = (val == "true");
    }

    // ==================================================================
    // Push datasets
    // ------------------------------------------------------------------

    const std::string BASEPATH = gra::aux::GetBasePath(2);

    fitcentral::DATASET temp;

    // ------------------------------------------------------------------
    // 0

    if (std::find(input.begin(), input.end(), 0) != input.end()) {
      temp.INPUTFILE = BASEPATH + "/fitcard/" + "ALICE_7_pipi.json";
      temp.DATAFILE  = BASEPATH + "/HEPdata/CD/" + "pipi_7_TeV_eta_09_pt_150.csv";
      temp.DATATYPE  = 4;

      temp.OBSERVABLE = 1;
      temp.MMIN       = 0.3;
      temp.MMAX       = 3.0;
      temp.NBINS      = 100;
      temp.h1DATA     = new MH1<double>(temp.NBINS, temp.MMIN, temp.MMAX, "ALICE 7 TeV pi+pi-");

      fitcentral::datasets.push_back(temp);
    }

    // ------------------------------------------------------------------
    // 1

    if (std::find(input.begin(), input.end(), 1) != input.end()) {
      temp.INPUTFILE = BASEPATH + "/fitcard/" + "ALICE_7_KK.json";
      temp.DATAFILE  = BASEPATH + "/HEPdata/CD/" + "KK_7_TeV.csv";
      temp.DATATYPE  = 1;

      temp.OBSERVABLE = 1;
      temp.MMIN       = 0.98;

      temp.MMAX   = 2.5;
      temp.NBINS  = 50;
      temp.h1DATA = new MH1<double>(temp.NBINS, temp.MMIN, temp.MMAX, "ALICE 7 TeV K+K-");

      fitcentral::datasets.push_back(temp);
    }

    // ------------------------------------------------------------------
    // 2

    if (std::find(input.begin(), input.end(), 2) != input.end()) {
      temp.INPUTFILE = BASEPATH + "/fitcard/" + "ATLAS_13_pipi.json";
      temp.DATAFILE  = BASEPATH + "/HEPdata/CD/" + "ATLAS13pipi.csv";
      temp.DATATYPE  = 1;

      temp.OBSERVABLE = 1;
      temp.MMIN       = 0.3;
      temp.MMAX       = 2.4;
      temp.NBINS      = 100;
      temp.h1DATA     = new MH1<double>(temp.NBINS, temp.MMIN, temp.MMAX, "ATLAS 13 TeV pi+pi-");

      fitcentral::datasets.push_back(temp);
    }

    if (fitcentral::datasets.size() == 0) {
      throw std::invalid_argument("fitcentral::datasets is empty (did not find any)");
    }

    // Read input data
    for (const auto &k : indices(fitcentral::datasets)) {
      if (!fitcentral::ReadData(fitcentral::datasets[k])) return EXIT_FAILURE;
    }

    // Initialize TMinuit with a maximum of 50 params
    std::unique_ptr<TMinuit> gMinuit = std::make_unique<TMinuit>(50);
    gMinuit->SetFCN(fitcentral::Chi2Func);

    double arglist[10];
    int    ierflg = 0;

    // Chi^2 type cost function, set 1
    // -log(likelihood), set 0.5
    arglist[0] = 1;
    gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);

    // Create generator object first to be able to read model parameters
    std::unique_ptr<MGraniitti> gen = std::make_unique<MGraniitti>();

    // Read process steering input (just the first one is enough)
    gen->ReadInput(fitcentral::datasets[0].INPUTFILE);
    std::map<std::string, gra::PARAM_RES> RESONANCES = gen->proc->GetResonances();

    // SET SOFT PARAMETERS, phase angles between [-pi,pi] but for optimization
    // reasons we allow larger domain [-2pi, 2pi]
    // due to phase wrapping.

    const double STEPSIZE_A   = 0.05;  // starting scale
    const double STEPSIZE_phi = 0.15;

    // Loop over resonances
    int k = -1;

    if (FIXALLRES) {
      for (const auto &x : RESONANCES) {
        fitcentral::afixlist.push_back(x.first);
        fitcentral::pfixlist.push_back(x.first);
      }
    }

    for (const auto &x : RESONANCES) {
      bool setzero = false;
      if (std::find(fitcentral::zerolist.begin(), fitcentral::zerolist.end(), x.first) !=
          fitcentral::zerolist.end()) {
        setzero = true;
      }
      gMinuit->mnparm(++k, Form("RESONANCES[%s].g_A", x.first.c_str()),
                      !setzero ? std::abs(x.second.g) : 0.0, STEPSIZE_A, 1e-9, 20.0,
                      ierflg);  // couplings [1e-9,...,2.0]
      gMinuit->mnparm(++k, Form("RESONANCES[%s].g_phi", x.first.c_str()),
                      !setzero ? std::arg(x.second.g) : 0.0, STEPSIZE_phi, -2 * gra::math::PI,
                      2 * gra::math::PI, ierflg);

      // Amplitude and Phase
      bool afix = std::find(fitcentral::afixlist.begin(), fitcentral::afixlist.end(), x.first) !=
                  fitcentral::afixlist.end();
      bool pfix = std::find(fitcentral::pfixlist.begin(), fitcentral::pfixlist.end(), x.first) !=
                  fitcentral::pfixlist.end();

      // Is set to zero (then fix always)
      if (setzero) {
        afix = true;
        pfix = true;
      }
      if (afix) { gMinuit->FixParameter(k - 1); }
      if (pfix) { gMinuit->FixParameter(k); }
    }

    // Continuum off-shell form factor
    if (fitcentral::offshellFF == 0) {
      gMinuit->mnparm(++k, "PARAM_REGGE::a_EXP", PARAM_REGGE::b_EXP, STEPSIZE_A, 0.0, 2.0, ierflg);
    } else if (fitcentral::offshellFF == 1) {
      gMinuit->mnparm(++k, "PARAM_REGGE::a_OREAR", PARAM_REGGE::a_OREAR, STEPSIZE_A, 0.0, 2.0,
                      ierflg);
      gMinuit->mnparm(++k, "PARAM_REGGE::b_OREAR", PARAM_REGGE::b_OREAR, STEPSIZE_A, 0.0, 2.0,
                      ierflg);
    } else if (fitcentral::offshellFF == 2) {
      gMinuit->mnparm(++k, "PARAM_REGGE::b_POW", PARAM_REGGE::b_POW, STEPSIZE_A, 0.0, 2.0, ierflg);
    } else {
      throw std::invalid_argument("fitcentral::offshellFF should be 0,1,2 (input is " +
                                  std::to_string(fitcentral::offshellFF) + ")");
    }

    // Do we fix it
    bool afix = std::find(fitcentral::afixlist.begin(), fitcentral::afixlist.end(), "continuum") !=
                fitcentral::afixlist.end();
    bool pfix = std::find(fitcentral::pfixlist.begin(), fitcentral::pfixlist.end(), "continuum") !=
                fitcentral::pfixlist.end();

    if (afix || pfix) {
      gMinuit->FixParameter(k);  // continuum
    }

    arglist[0] = 10000;  // Maximum number of iterations
    arglist[1] = 0.001;  // Cost function stop tolerance

    // Release all parameters
    // gMinuit->mnfree(0);

    // First MC search (necessary for global search)
    // gMinuit->mnexcm("SEEK", arglist, 2, ierflg);

    // First simplex (necessary for global search)
    gMinuit->mnexcm("SIMPLEX", arglist, 2, ierflg);

    // Then migrad near the optimum
    // gMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);

    // Get proper errors with MINOS
    // gMinuit->mnexcm("MINOS", arglist, 2, ierflg);

    // Print results
    double amin, edm, errdef;
    int    nvpar, nparx, icstat;
    gMinuit->mnstat(amin, edm, errdef, nvpar, nparx, icstat);
    // gMinuit->mnprin(3,amin);

    for (const auto &i : indices(fitcentral::datasets)) { delete fitcentral::datasets[i].h1DATA; }

    // Done
    std::cout << "[fitcentral:: done]" << std::endl;

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

  return EXIT_SUCCESS;
}
