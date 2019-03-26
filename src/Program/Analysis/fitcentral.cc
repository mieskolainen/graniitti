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
#include "Graniitti/MRegge.h"
#include "Graniitti/MAux.h"
#include "Graniitti/MGraniitti.h"
#include "Graniitti/MH1.h"
#include "Graniitti/MMath.h"

// ROOT
#include "TMinuit.h"
#include "TString.h"

// Libraries
#include "cxxopts.hpp"
#include "rang.hpp"


using gra::aux::indices;

using namespace gra;

namespace fitcentral {


// Single dataset / observable
struct DATASET {
    std::vector<double> data; // Events

    std::string INPUTFILE;
    std::string DATAFILE;
    int DATATYPE = 0;

    int OBSERVABLE = 0;
    uint NBINS  = 0;
    double MMIN = 0.0;
    double MMAX = 0.0;

    MH1<double>* h1DATA;
};


// All datasets here
std::vector<DATASET> datasets;


// Function prototypes
void RunMC(double* par);
bool ReadData(DATASET& dataset);
void Chi2Func(int& npar, double* gin, double& f, double* par, int iflag);
void SetSoftParam(std::map<std::string, gra::PARAM_RES>& RESONANCES, double* par);


// Choose continuum form factor parametrization here (1,2,3)
int offshellFF = 1;

bool POMLOOP = false;

// ===============================================================
// Note that TMinuit uses Fortran indexing => The parameter 0 in C++
// is the parameter 1 in TMinuit.
//
// Main function


// Set soft parameters back to the model
void SetSoftParam(std::map<std::string, gra::PARAM_RES>& RESONANCES, double* par) {

    printf("fitcentral::SetSoftParam: \n\n");
    
    PARAM_REGGE::offshellFF = offshellFF;
    int k = -1;

    // Loop over resonances
    for (const auto& x : RESONANCES) {
    	double A = par[++k];
    	double phi = par[++k];
    	RESONANCES[x.first].g = A * std::exp(gra::math::zi * phi);

    	printf("[%s] \t A = %0.8f \t phi = %0.6f (%0.0f deg)\n",
    	       x.first.c_str(), A, phi, phi / (2 * gra::math::PI) * 360);
    }
    printf("\n");
    
    // Continuum form factors
    if (offshellFF == 1) {
	PARAM_REGGE::b_EXP = par[++k];
	printf("[Offshell-FF] \t b_EXP = %0.3f \n", PARAM_REGGE::b_EXP);
    }
    if (offshellFF == 2) {
	PARAM_REGGE::a_OREAR = par[++k];
	PARAM_REGGE::b_OREAR = par[++k];
	printf("[Offshell-FF] \t a/b_OREAR = [%0.3f,%0.3f] \n",
	       PARAM_REGGE::a_OREAR, PARAM_REGGE::b_OREAR);
    }
    if (offshellFF == 3) {
	PARAM_REGGE::b_POW = par[++k];
	printf("[Offshell-FF] \t b_POW = %0.3f \n", PARAM_REGGE::b_POW);
    }
}


// Read input data event by event
bool ReadData(DATASET& dataset) {
    printf("fitcentral::ReadData: with inputfile %s \n",
           dataset.DATAFILE.c_str());

    FILE* fp;
    double M = 0.0;
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

	if (dataset.DATATYPE == 4)
	    ret = fscanf(fp, "%lf,%lf,%lf,%lf", &M, &Px, &Py, &Pz);

	if (dataset.DATATYPE == 1)
	    ret = fscanf(fp, "%lf", &M);

	if (ret == dataset.DATATYPE) { // number of values per line
	                               // Fine
	} else if (ret == EOF) {
	    break;
	} else {
	    printf("Error in the file structure of %s!\n",
	           dataset.DATAFILE.c_str());
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

MEikonal eik; // For speed, calculate only once
int iter = 0;


// Cost function
void Chi2Func(int& npar, double* gin, double& f, double* par, int iflag) {

    std::cout << rang::fg::red << "fitcentral::Chi2Func: iter = " << iter << rang::fg::reset << std::endl;

    double cost = 0;
    double chi2 = 0;
    int points = 0;

    // Loop over different datasets
    for (const auto& k : indices(datasets)) {

    	double local_cost = 0.0;
    	double local_chi2 = 0.0;
    	int local_points = 0;

    	// Generator
        std::unique_ptr<MGraniitti> gen = std::make_unique<MGraniitti>();

    	try {

    	    if (iter > 1) {
    		HILJAA = true; // SILENT output
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
    	    //gen->Generate();

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
        		                    datasets[k].h1DATA->GetBinWeight(i)) / datasets[k].h1DATA->GetBinWeight(i);

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
    	    chi2   += local_chi2;
    	    cost   += local_cost;
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
    	    printf(" chi^2 / points = %0.2f / %d = %0.2f \n", local_chi2,
    	           local_points, local_chi2 / (double)local_points);
    	    printf(" -logL = %0.2f \n", local_cost);
    	    printf("GLOBAL \n");
    	    printf(" chi^2 / points = %0.2f / %d = %0.2f \n", chi2, points,
    	           chi2 / (double)points);
    	    printf(" -logL = %0.2f \n", cost);
    	    printf(
    	        "============================================================="
    	        " \n\n");
    	    std::cout << rang::fg::reset;

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
        std::cerr << rang::fg::red << "Exception catched: Unspecified (...) (Probably JSON input)"
        		  << rang::fg::reset << std::endl;
        exit(0);
        }

    } // Over different datasets

    f = cost;
}
} // Namespace fitcentral ends


int main(int argc, char* argv[]) {

    gra::aux::PrintFlashScreen(rang::fg::red);
    std::cout << rang::style::bold
              << "GRANIITTI - Soft Central Production Fitter"
              << rang::style::reset << std::endl
              << std::endl;
    gra::aux::PrintVersion();

    // Save the number of input arguments
    const int NARGC = argc - 1;

try {

    cxxopts::Options options(argv[0], "");
    options.add_options()
        ("i,input",     "Input dataset collection numbers         <0,1,2,3,...>", cxxopts::value<std::string>() )
        ("c,continuum", "Continuum form factor                    <0|1|2>",       cxxopts::value<int>() )
        ("x,fix",       "Fix all resonances to their input values <true|false>",  cxxopts::value<std::string>() )
        ("l,POMLOOP",   "Screening Pomeron loop", cxxopts::value<std::string>() )
        ("H,help",      "Help")
        ;
    auto r = options.parse(argc, argv);

    if (r.count("help") || NARGC == 0) {
        std::cout << options.help({""}) << std::endl;
        std::cout << "Example:" << std::endl;
        std::cout << "  " << argv[0] << " -i 0,1 -c 2"
                  << std::endl << std::endl;
        return EXIT_FAILURE;
    }

    // Input dataset numbers
    std::vector<int> input = gra::aux::SplitStr2Int(r["input"].as<std::string>());

    // Off-shell form factor
    fitcentral::offshellFF = r["c"].as<int>(); 

    // Fix-all resonances
    bool FIXALLRES = false;
    if (r.count("x")) {
        FIXALLRES = (r["x"].as<std::string>() == "true" ? true : false);
    }

    // Pomeron loop
    if (r.count("l"))  { 
        const std::string val = r["l"].as<std::string>();
        fitcentral::POMLOOP = (val == "true");
    }

    // ==================================================================
    // Push datasets
    // ------------------------------------------------------------------

    const std::string BASEPATH = gra::aux::GetBasePath(2);

    fitcentral::DATASET temp;


    // ------------------------------------------------------------------
    // 0

    if (std::find(input.begin(), input.end(), 0) != input.end() ) {

        temp.INPUTFILE = BASEPATH + "/fitcard/"    + "ALICE_7_pipi.json";
        temp.DATAFILE  = BASEPATH + "/HEPdata/CD/" + "pipi_7_TeV_eta_09_pt_150.csv";
        temp.DATATYPE  = 4;

        temp.OBSERVABLE = 1;
        temp.MMIN   = 0.3;
        temp.MMAX   = 3.0;
        temp.NBINS  = 100;
        temp.h1DATA = new MH1<double>(temp.NBINS, temp.MMIN, temp.MMAX, "ALICE 7 TeV pi+pi-");

        fitcentral::datasets.push_back(temp);
    }

    // ------------------------------------------------------------------
    // 1

    if (std::find(input.begin(), input.end(), 1) != input.end() ) {

        temp.INPUTFILE = BASEPATH + "/fitcard/"    + "ALICE_7_KK.json";
        temp.DATAFILE  = BASEPATH + "/HEPdata/CD/" + "KK_7_TeV.csv";
        temp.DATATYPE  = 1;

        temp.OBSERVABLE = 1;
        temp.MMIN = 0.98;

        temp.MMAX = 2.5;
        temp.NBINS = 50;
        temp.h1DATA = new MH1<double>(temp.NBINS, temp.MMIN, temp.MMAX, "ALICE 7 TeV K+K-");

        fitcentral::datasets.push_back(temp);
    }

    // ------------------------------------------------------------------
    // 2

    if (std::find(input.begin(), input.end(), 2) != input.end() ) {

        temp.INPUTFILE = BASEPATH + "/fitcard/"    + "ATLAS_13_pipi.json";
        temp.DATAFILE  = BASEPATH + "/HEPdata/CD/" + "ATLAS13pipi.csv";
        temp.DATATYPE  = 1;

        temp.OBSERVABLE = 1;
        temp.MMIN   = 0.3;
        temp.MMAX   = 2.4;
        temp.NBINS  = 100;
        temp.h1DATA = new MH1<double>(temp.NBINS, temp.MMIN, temp.MMAX, "ATLAS 13 TeV pi+pi-");

        fitcentral::datasets.push_back(temp);
    }

    if (fitcentral::datasets.size() == 0) {
        throw std::invalid_argument("fitcentral::datasets is empty (did not find any)");
    }

    // Read input data
    for (const auto& k : indices(fitcentral::datasets)) {
	if (!fitcentral::ReadData(fitcentral::datasets[k]))
	    return EXIT_FAILURE;
    }

    // Initialize TMinuit with a maximum of 50 params
    std::unique_ptr<TMinuit> gMinuit = std::make_unique<TMinuit>(50);
    gMinuit->SetFCN(fitcentral::Chi2Func);

    double arglist[10];
    int ierflg = 0;
    
    
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

    const double STEPSIZE_A   = 0.1;
    const double STEPSIZE_phi = 0.15;

    
    // Loop over resonances
    int k = -1;

    for (const auto& x : RESONANCES) {
	gMinuit->mnparm(++k, Form("RESONANCES[%s].g_A", x.first.c_str()),
	                std::abs(x.second.g), STEPSIZE_A, 1e-9, 20.0, ierflg); // couplings [1e-9,...,2.0]
	gMinuit->mnparm(++k, Form("RESONANCES[%s].g_phi", x.first.c_str()),
	                std::arg(x.second.g), STEPSIZE_phi, -2 * gra::math::PI,
	                2 * gra::math::PI, ierflg);
    }
    
    // Continuum off-shell form factor
    if (fitcentral::offshellFF == 1) {
	gMinuit->mnparm(++k, "PARAM_REGGE::a_EXP", PARAM_REGGE::b_EXP,
	                STEPSIZE_A, 0.0, 2.0, ierflg);
    }
    if (fitcentral::offshellFF == 2) {
	gMinuit->mnparm(++k, "PARAM_REGGE::a_OREAR", PARAM_REGGE::a_OREAR,
	                STEPSIZE_A, 0.0, 2.0, ierflg);
	gMinuit->mnparm(++k, "PARAM_REGGE::b_OREAR", PARAM_REGGE::b_OREAR,
	                STEPSIZE_A, 0.0, 2.0, ierflg);
    }
    if (fitcentral::offshellFF == 3) {
	gMinuit->mnparm(++k, "PARAM_REGGE::b_POW", PARAM_REGGE::b_POW,
	                STEPSIZE_A, 0.0, 2.0, ierflg);
    }

    arglist[0] = 10000; // Maximum number of iterations
    arglist[1] = 0.001; // Cost function stop tolerance

    // FIX ALL RESONANCES
    if (FIXALLRES) {
        for (std::size_t i = 0; i < RESONANCES.size() * 2; ++i) { // 2 x for amplitude and phase
            gMinuit->FixParameter(i);
        }
    }

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
    int nvpar, nparx, icstat;
    gMinuit->mnstat(amin, edm, errdef, nvpar, nparx, icstat);
    // gMinuit->mnprin(3,amin);

    for (const auto& i : indices(fitcentral::datasets)) {
	   delete fitcentral::datasets[i].h1DATA;
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
