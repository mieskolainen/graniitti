// GRANIITTI - Monte Carlo event generator for high energy diffraction
// https://github.com/mieskolainen/graniitti

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
#include "Graniitti/MMath.h"
#include "Graniitti/MPDG.h"

// ROOT
#include "TMinuit.h"
#include "TString.h"

// Libraries
#include "rang.hpp"


using gra::aux::indices;
using namespace gra::form;
using namespace gra;

// Function prototypes
namespace fitsoft {

const std::string INPUTFILE = gra::aux::GetBasePath(2) + "/fitcard/" + "EL.json";

int iter = 0;

// Global for the differential fit
namespace DIFFMEAS {

    std::vector<double> sqrts;
    std::vector<std::vector<std::string>> initialstate;

    std::vector<std::vector<double>> x;
    std::vector<std::vector<double>> y;
    std::vector<std::vector<double>> err;    
}

// Measurements and errors for each energy
namespace INTEGMEAS {

    // Which datasets are active in the fit
    std::vector<bool> active        = {false, false, true, true, false};

    std::vector<double> energy      = {31, 62, 1800, 7000, 8000};
    std::vector<std::vector<std::string>> beam = {{"p+","p-"}, {"p+","p-"}, {"p+","p-"}, {"p+","p+"}, {"p+","p+"}}; 
    
    std::vector<double> sigma_tot   = {40.2, 43.9, 80.03, 98.3, 96.07};
    std::vector<double> sigma_tot_e = {0.2, 0.3, 2.24, 2.8, 0.9225};

    std::vector<double> sigma_el    = {7.12, 7.61, 19.7, 24.3, 24.33};
    std::vector<double> sigma_el_e  = {0.34, 0.22, 0.85, 1.2, 0.392};
}


std::vector<double> dsigma_el_dt(double* par, double sqrts, const std::vector<string>& beam,
                                 std::vector<double>& x);

void SigmaXS(double* par, double& xs_tot, double& xs_el, double& xs_inel,
             const double sqrts, const std::vector<std::string>& beam);

bool ReadData(std::string inputfile, std::vector<double>& x,
              std::vector<double>& y, std::vector<double>& err);

void PushNull();
bool InitData();
void Chi2Func(int& npar, double* gin, double& f, double* par, int iflag);
void SetSoftParam(double* par);


// Differential elastic
std::vector<double> dsigma_el_dt(double* par, const double sqrts,
                                 const std::vector<std::string>& beam,
                                 std::vector<double>& x) {
    if (iter > 0) {
	   HILJAA = true; // SILENT output
    }
    
    // Mandelstam s
    const double s = sqrts * sqrts;

    // Create generator object first
    std::unique_ptr<MGraniitti> gen = std::make_unique<MGraniitti>();
    std::vector<double> values(x.size(), 0.0);

    try {
	// Read process input from file
	gen->ReadInput(fitsoft::INPUTFILE);

    // Set beam and energy
    const std::vector<double> energy = {sqrts/2, sqrts/2};
    gen->proc->SetInitialState(beam, energy);

	// ** Set soft parameters here **
	fitsoft::SetSoftParam(par);

	// Initialize generator, always last
	gen->InitEikonal();

	// Construct dsigma_el/dt
	for (const auto& i : indices(values)) {
	    // Get amplitude
	    const double kt2 = x[i]; // |t| (Mandelstam)
	    const std::complex<double> A =
	        gen->proc->Eikonal.MSA.Interpolate1D(kt2);

	    // dsigma/dt in units of millibarns
	    values[i] = gra::math::abs2(A) / (16.0 * gra::math::PI * s *
	                                  (s - 4.0 * gra::math::pow2(PDG::mp))) * PDG::GeV2mb;
	}
    } catch (const std::invalid_argument& e) {
    gra::aux::PrintGameOver();
    std::cerr << rang::fg::red
              << "Exception catched: " << rang::fg::reset << e.what()
              << std::endl;
    exit(0);
    } catch (const std::ios_base::failure& e) {
    gra::aux::PrintGameOver();
    std::cerr << rang::fg::red
              << "Exception catched: std::ios_base::failure: "
              << rang::fg::reset << e.what() << std::endl;
    exit(0);
    } catch (...) {
    gra::aux::PrintGameOver();
    std::cerr << rang::fg::red << "Exception catched: Unspecified (...) (Probably JSON input)" << rang::fg::reset
              << std::endl;
    exit(0);
    }

    ++iter;

    return values;
}


// Integrated cross sections
void SigmaXS(double* par, double& xs_tot, double& xs_el, double& xs_inel,
             const double sqrts, const std::vector<std::string>& beam) {
    // Create generator object first
    std::unique_ptr<MGraniitti> gen = std::make_unique<MGraniitti>();

    try {
	// Read process input from file and set initial state
	gen->ReadInput(fitsoft::INPUTFILE);
    const std::vector<double> energy = {sqrts/2, sqrts/2};
    gen->proc->SetInitialState(beam, energy);

	// ** Set soft parameters **
	fitsoft::SetSoftParam(par);

	// Initialize generator, always last
	gen->InitEikonalOnly();
	gen->proc->Eikonal.GetTotXS(xs_tot, xs_el, xs_inel);

    } catch (const std::invalid_argument& e) {
	gra::aux::PrintGameOver();
	std::cerr << "Exception catched: " << e.what() << std::endl;
	exit(0);
    } catch (const std::ios_base::failure& e) {
	gra::aux::PrintGameOver();
	std::cerr << "Exception catched: Caught std::ios_base::failure: "
	          << e.what() << std::endl;
	exit(0);
    }
}


// Read differential input
bool ReadData(std::string inputfile, std::vector<double>& x,
              std::vector<double>& y, std::vector<double>& err) {
    printf("fitsoft::ReadData: using HepDATA inputfile %s \n",
           inputfile.c_str());

    double xval        = 0.0;
    double xval_low    = 0.0;
    double xval_high   = 0.0;
    double yval        = 0.0;

    double yerr_stat_p = 0.0;
    double yerr_stat_m = 0.0;

    double yerr_syst_p = 0.0;
    double yerr_syst_m = 0.0;

    FILE* fp;

    if ((fp = fopen(inputfile.c_str(), "r+")) == NULL) {
	printf("No inputfile %s found \n", inputfile.c_str());
	return false;
    }

    while (true) {
	int ret = fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf", &xval,
	                 &xval_low, &xval_high, &yval, &yerr_stat_p,
	                 &yerr_stat_m, &yerr_syst_p, &yerr_syst_m);
	if (ret == 8) {
	    // Fine
	} else if (ret == EOF) {
	    break;
	} else {
	    printf("Error in the file structure of %s!\n", inputfile.c_str());
	    return false;
	}

	x.push_back(xval);
	y.push_back(yval);

	// Construct error
	double error = std::sqrt(gra::math::pow2((yerr_stat_p - yerr_stat_m) / 2.0) +
	                         gra::math::pow2((yerr_syst_p - yerr_syst_m) / 2.0));
	err.push_back(error);

	printf("x = %0.8f, y = %0.8f, err = %0.8f, err/y = %0.8f \n", xval,
	       yval, error, error / yval);
    }

    printf("fitsoft::ReadData: [done] \n");

    return true;
}

// Call this for every dataset
void PushNull() {
    // Make empty vectors
    std::vector<double> nullvec;
    nullvec.clear();

    DIFFMEAS::sqrts.push_back(0.0);
    DIFFMEAS::initialstate.push_back({"",""});

    DIFFMEAS::x.push_back(nullvec);
    DIFFMEAS::y.push_back(nullvec);
    DIFFMEAS::err.push_back(nullvec);
}


// Differential dsigma/dt dataset initialization
//
// Add here new datasets as you wish
bool InitData() {
    printf("fitsoft::InitData: \n");

    uint N = 0;
    const std::string BASEPATH = gra::aux::GetBasePath(2) + "/HEPdata/EL/";
    
    // 7 TeV DATA
    fitsoft::PushNull();
    DIFFMEAS::sqrts[N] = 7000;
    DIFFMEAS::initialstate[N] = {"p+","p+"};
    if (!fitsoft::ReadData(BASEPATH + "TOTEM_7_low_t.csv", DIFFMEAS::x[N], DIFFMEAS::y[N], DIFFMEAS::err[N])) {
	   return false;
    }
    if (!fitsoft::ReadData(BASEPATH + "TOTEM_7.csv",       DIFFMEAS::x[N], DIFFMEAS::y[N], DIFFMEAS::err[N])) {
	   return false;
    }
    ++N;
    
    // 1.96 TeV DATA
    /*
    fitsoft::PushNull();
    DIFFMEAS::sqrts[N] = 1960;
    DIFFMEAS::initialstate[N] = {"p+","p-"};
    if (!fitsoft::ReadData(BASEPATH + "D0_1960.csv", DIFFMEAS::x[N], DIFFMEAS::y[N], DIFFMEAS::err[N]) ) {
        return false;
    }
    ++N;
    */

    // 1.8 TeV DATA
    fitsoft::PushNull();
    DIFFMEAS::sqrts[N] = 1800;
    DIFFMEAS::initialstate[N] = {"p+","p-"};
    if (!fitsoft::ReadData(BASEPATH + "ABE_1994_1800.csv", DIFFMEAS::x[N], DIFFMEAS::y[N], DIFFMEAS::err[N])) {
	   return false;
    }
    ++N;

    // 546 GeV DATA
    fitsoft::PushNull();
    DIFFMEAS::sqrts[N] = 546;
    DIFFMEAS::initialstate[N] = {"p+","p-"};
    if (!fitsoft::ReadData(BASEPATH + "FNAL_546.csv", DIFFMEAS::x[N], DIFFMEAS::y[N], DIFFMEAS::err[N])){
	   return false;
    }
    if (!fitsoft::ReadData(BASEPATH + "SPS_546.csv",  DIFFMEAS::x[N], DIFFMEAS::y[N], DIFFMEAS::err[N])){
	   return false;
    }
    ++N;
    /*
      // 62.5 GeV DATA (too low energy, needs Reggeons)
      PushNull();
      DIFFMEAS::sqrts[N] = 62.5;
      DIFFMEAS::initialstate[N] = {"p+","p-"};
      if (!ReadData(BASEPATH + "ISR_62.csv", DIFFMEAS::x[N], DIFFMEAS::y[N], DIFFMEAS::err[N]) )
      return false;
      ++N;
    */
    return true;
}

// Chi^2 function
void Chi2Func(int& npar, double* gin, double& f, double* par, int iflag) {

    printf("fitsoft::Chi2Func: [iter = %d] \n", iter);
    double chi2 = 0;
    int points = 0;

    // Loop over at different energies
    for (const auto& i : indices(DIFFMEAS::sqrts)) {

    	// Get dsigma_el/dt values
    	std::vector<double> values = fitsoft::dsigma_el_dt(
    	    par, DIFFMEAS::sqrts[i], DIFFMEAS::initialstate[i], DIFFMEAS::x[i]);

    	for (const auto& j : indices(values)) {
    	    chi2 += gra::math::pow2((DIFFMEAS::y[i][j] - values[j]) / DIFFMEAS::err[i][j]);
    	    ++points;
    	}
    }

    // Integrated cross sections
    double xs_tot  = 0.0;
    double xs_el   = 0.0;
    double xs_inel = 0.0;

    // Loop over energies
    for (const auto& i : indices(INTEGMEAS::active)) {

      if (INTEGMEAS::active[i] == true) {
        SigmaXS(par, xs_tot, xs_el, xs_inel, INTEGMEAS::energy[i], INTEGMEAS::beam[i]);
        
        chi2 += gra::math::pow2((INTEGMEAS::sigma_tot[i] - xs_tot)/INTEGMEAS::sigma_tot_e[i]);
        chi2 += gra::math::pow2((INTEGMEAS::sigma_el[i] - xs_el)/INTEGMEAS::sigma_el_e[i]);
        ++points;
      }
    }
    
    // Print out parameter values
    PARAM_SOFT::PrintParam();

    if (chi2 / (double)points < 2.0) {
	   std::cout << rang::fg::green;
    } else {
	   std::cout << rang::fg::red;
    }
    printf("============================================================= \n");
    printf(" chi^2 / points = %0.2f / %d = %0.2f \n", chi2, points, chi2 / (double)points);
    printf("============================================================= \n\n");
    std::cout << rang::fg::reset;

    f = chi2;
}


// Set soft parameters
// These need to follow the same order as below in main()
void SetSoftParam(double* par) {

    printf("fitsoft::SetSoftParam: \n\n");

    PARAM_SOFT::DELTA_P = par[0];
    PARAM_SOFT::ALPHA_P = par[1];
    PARAM_SOFT::gN_P    = par[2];
    PARAM_SOFT::g3P     = par[3] * par[2]; // Convention
    PARAM_SOFT::gamma   = par[4];

    PARAM_SOFT::fc1     = par[5];
    PARAM_SOFT::fc2     = par[6];
    PARAM_SOFT::fc3     = par[7];

    // Initialization loops -> put just enough for the event by event accuracy
    MEikonalNumerics::NumberBT  = 500;
    MEikonalNumerics::NumberKT2 = 500;
}

} // Namespace fitsoft ends


// ===============================================================
// Note that MINUIT uses Fortran indexing, i.e. parameter 0 in C++
// is parameter 1 in MINUIT.
//
// Main function
int main(int argc, char* argv[]) {

    gra::aux::PrintFlashScreen(rang::fg::red);
    std::cout << rang::style::bold
              << "GRANIITTI - Eikonal Pomeron Fitter"
              << rang::style::reset << std::endl
              << std::endl;
    gra::aux::PrintVersion();

    // Initialize input data
    if (!fitsoft::InitData()) {
	   return EXIT_FAILURE;
    }

    // Maximum number of parameters of the soft model
    const int NPARAM = 20;

    // Init TMinuit
    std::unique_ptr<TMinuit> gMinuit = std::make_unique<TMinuit>(
        NPARAM); // initialize TMinuit with a maximum of NPARAM params
    gMinuit->SetFCN(fitsoft::Chi2Func);

    double arglist[10];
    int ierflg = 0;

    // Chi^2 type cost function (for -log(likelihood) type put 0.5)
    arglist[0] = 1; 
    gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);

    // Create generator object first to be able to read model parameters
    MGraniitti gen;

    // Read process input from file
    gen.ReadInput(fitsoft::INPUTFILE);

    // SET SOFT PARAMETERS
    // Pomeron trajectors
    gMinuit->mnparm(0, "DELTA_P", PARAM_SOFT::DELTA_P, 0.01, -0.5, 0.5, ierflg);
    gMinuit->mnparm(1, "ALPHA_P", PARAM_SOFT::ALPHA_P, 0.01, 1e-9, 0.5, ierflg);

    // Couplings
    gMinuit->mnparm(2, "gN_P",  PARAM_SOFT::gN_P,   0.1, 1e-9, 12.0, ierflg);
    gMinuit->mnparm(3, "g3P",   PARAM_SOFT::g3P / PARAM_SOFT::gN_P, 0.1, 1e-9, 1.0, ierflg); // Convention
    gMinuit->mnparm(4, "gamma", PARAM_SOFT::gamma, 0.01, 1e-9, 1.0 - 1e-9, ierflg);

    // Proton (+ pion loop) form factor parameters (put maximum at least 15 GeV)
    gMinuit->mnparm(5, "fc1", PARAM_SOFT::fc1, 0.1, 1e-9, 100.0, ierflg);
    gMinuit->mnparm(6, "fc2", PARAM_SOFT::fc2, 0.1, 1e-9, 100.0, ierflg);
    gMinuit->mnparm(7, "fc3", PARAM_SOFT::fc3, 0.1, 1e-9, 100.0, ierflg);

    // Fix triple Pomeron coupling and eikonal gamma (we do not fit inelastic
    // diffraction here)
    {
    	const std::vector<int> fixed = {3, 4};
    	for (const auto& x : fixed) {
    	    gMinuit->FixParameter(x);
    	}
    }
    arglist[0] = 10000; // Maximum number of iterations
    arglist[1] = 0.01;  // Cost function tolerance
    
    // Now ready for minimization step

    // First simplex
    // gMinuit->mnexcm("SIMPLEX", arglist, 2, ierflg);

    // Then migrad
    gMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);

    // Get proper errors with MINOS
    gMinuit->mnexcm("MINOS", arglist, 2, ierflg);

    // Print results
    double amin, edm, errdef;
    int nvpar, nparx, icstat;
    gMinuit->mnstat(amin, edm, errdef, nvpar, nparx, icstat);
    // gMinuit->mnprin(3,amin);

    return EXIT_SUCCESS;
}

