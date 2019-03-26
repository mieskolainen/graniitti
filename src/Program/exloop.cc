// GRANIITTI - Monte Carlo event generator for high energy diffraction
// https://github.com/mieskolainen/graniitti
//
// <Fiducial Measurements vs Monte Carlo>
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
#include <memory>
#include <mutex>
#include <random>
#include <stdexcept>
#include <thread>
#include <vector>

// Own
#include "Graniitti/M4Vec.h"
#include "Graniitti/MAux.h"
#include "Graniitti/MGraniitti.h"
#include "Graniitti/MTimer.h"

// Libraries
#include "rang.hpp"
#include "cxxopts.hpp"

using gra::aux::indices;
using namespace gra;


struct MEASUREMENT {
    MEASUREMENT(std::string card_, double value_, double stat_, double syst_) {
	card = card_;
	value = value_;
	stat = stat_;
	syst = syst_;
    }
    std::string card;
    double value = 0.0;
    double stat  = 0.0;
    double syst  = 0.0;
};


void experiment(bool screening);


// Main
int main(int argc, char* argv[]) {

    gra::aux::PrintFlashScreen(rang::fg::green);
    std::cout << rang::style::bold
              << "GRANIITTI - Fiducial Measurements vs Monte Carlo Processes"
              << rang::style::reset << std::endl
              << std::endl;
    gra::aux::PrintVersion();

    // Save the number of input arguments
    const int NARGC = argc - 1;
    try {

    cxxopts::Options options(argv[0], "");
    options.add_options()
        ("l,POMLOOP", "Screening Pomeron loop (slow)", cxxopts::value<std::string>() )
        ("H,help",    "Help")
        ;
    auto r = options.parse(argc, argv);
    
    if (r.count("help") || NARGC == 0) {
        std::cout << options.help({""}) << std::endl;
        std::cout << "Example:" << std::endl;
        std::cout << "  " << argv[0] << " -l false"
                  << std::endl << std::endl;
        return EXIT_FAILURE;
    }

    bool screening = false;
	if (r.count("l"))  { 
		const std::string val = r["l"].as<std::string>();
		screening = (val == "true");
	}

	// Project
    experiment(screening);

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
		std::cerr << rang::fg::red << "Exception catched: Unspecified (...) (Probably JSON input)" << rang::fg::reset
		          << std::endl;
		return EXIT_FAILURE;
    }

    //autoloop();
    return EXIT_SUCCESS;
}


void experiment(bool screening) {

    MTimer timer;
    std::vector<MEASUREMENT> input;

    /*
      input.push_back(MEASUREMENT("./tests/processes/el.json",            1.6e-12, 0.5e-12, 0.3e-12));

      // <https://cds.cern.ch/record/1495761/files/CERN-PH-EP-2012-353.pdf>
      input.push_back(MEASUREMENT("./tests/processes/TOTEM12_.json",      1.6e-12, 0.5e-12, 0.3e-12));

      // <https://arxiv.org/abs/1712.06153>
      input.push_back(MEASUREMENT("./tests/processes/TOTEM17_total.json", 1.6e-12, 0.5e-12, 0.3e-12));
    */

    const std::string BASEPATH = "./tests/processes/";


    // ...
    input.push_back(MEASUREMENT("gg2MMbar.json",     0, 0, 0));
    

    // https://arxiv.org/pdf/hep-ex/0611040.pdf
    input.push_back(MEASUREMENT("CDF07_ee.json",     1.6e-12,   0.5e-12,   0.3e-12));

    // <https://arxiv.org/pdf/1112.0858.pdf>
    input.push_back(MEASUREMENT("CDF11_ee.json",     2.88e-12,  0.57e-12,  0.63e-12));

    // <https://arxiv.org/abs/1111.5536>
    input.push_back(MEASUREMENT("CMS11_mumu.json",   3.38e-12,  0.58e-12,  0.21e-12));

    // <https://arxiv.org/abs/1506.07098>
    input.push_back(MEASUREMENT("ATLAS15_ee.json",   0.428e-12, 0.035e-12, 0.018e-12));

    // <https://arxiv.org/abs/1506.07098>
    input.push_back(MEASUREMENT("ATLAS15_mumu.json", 0.628e-12, 0.032e-12, 0.021e-12));

    // <https://arxiv.org/abs/1708.04053>
    input.push_back(MEASUREMENT("ATLAS17_mumu.json", 3.12e-12,  0.07e-12,  0.14e-12));

    // <https://discoverycenter.nbi.ku.dk/teaching/thesis_page/MasterEmilBolsFinal.pdf>
    input.push_back(MEASUREMENT("ATLAS17_2pi.json",  18.75e-6,  0.048e-6,  0.770e-6));

    // Comparison of 2->4 vs 2->2 cf. MadGraph EPA (1.31e-12)
    input.push_back(MEASUREMENT("yy2ee.json",       0, 0, 0));
    input.push_back(MEASUREMENT("yy2ee_DZ.json",    0, 0, 0));
    

    // ALICE fiducial
    input.push_back(MEASUREMENT("ALICE_2pi.json",    0, 0, 0));

    input.push_back(MEASUREMENT("ALICE_2K.json",     0, 0, 0));

    input.push_back(MEASUREMENT("ALICE_4pi.json",    0, 0, 0));

    input.push_back(MEASUREMENT("ALICE_4K.json",     0, 0, 0));
    
    
    // <ALICE 7 TeV PWA, |y(pi+pi-)| < 0.9
    //input.push_back(MEASUREMENT("./tests/experiment/ALICE_2pi_PWA.json", 31e-6, 0.5e-6, 2e-6));
            
    // <https://arxiv.org/pdf/1806.04079.pdf>
    //input.push_back(MEASUREMENT("./tests/experiment/LHCb18_jpsi.json", 399e-12, 16e-12, 19e-12));


    // <ATLAS W+W->
    input.push_back(MEASUREMENT("ATLAS_WW.json",     0, 0, 0));

	// <https://arxiv.org/pdf/1806.04079.pdf>
    input.push_back(MEASUREMENT("gg2chic0.json",     0, 0, 0));

    // ...
    input.push_back(MEASUREMENT("gg2gg.json",        0, 0, 0));

    // ...
    input.push_back(MEASUREMENT("yy2Higgs.json",     0, 0, 0));

    // ...
    input.push_back(MEASUREMENT("monopolepair.json",     0, 0, 0));

    // ...
    input.push_back(MEASUREMENT("monopolepair_LUX.json", 0, 0, 0));

    // ...
    input.push_back(MEASUREMENT("monopolium.json",       0, 0, 0));


    // ------------------------------------------------------------------
    /*
    // Output text file
    FILE* fout;
    fout = fopen("burn.csv", "w");
    if (fout == NULL) {
		printf("Scan:: Error opening burn.csv file! \n");
		return;
    }
    FILE* fout_latex;
    fout_latex = fopen("burn.tex", "w");
    if (fout_latex == NULL) {
		printf("Scan:: Error opening burn.tex file! \n");
		return;
    }

    fprintf(fout, "sqrts\t\txstot\t\txsin\t\txsel");
    for (unsigned int k = 0; k < input.size(); ++k) {
		fprintf(fout, "\txs%d", k);
    }
    fprintf(fout, "\n");
    fflush(fout);
    */

    std::vector<std::vector<double>> xs0(input.size(), std::vector<double>(2, 0.0));
    std::vector<std::vector<double>> xs0_err = xs0;

    int MODEMAX = 1;
    if (screening) {
		MODEMAX = 2;
    }
	
    // LOOP over ALL input
    for (auto const& p : indices(input)) {

		// Screening off/on
		for (int mode = 0; mode < MODEMAX; ++mode) {
	    	
		    try {
                
		    	MGraniitti* gen = new MGraniitti;

				HILJAA = false;

				// Read process input from file
				gen->ReadInput(BASEPATH + input[p].card);

				// Set spesific parameters
				gen->proc->SetScreening(mode);
				gen->proc->SetExcitation(false);
                gen->proc->SetHistograms(0);

				// Initialize (always last!)
				gen->Initialize();
                
				// Get total, elastic and inelastic cross section
				double xs_tot = 0.0;
				double xs_el  = 0.0;
				double xs_in  = 0.0;
				gen->proc->Eikonal.GetTotXS(xs_tot, xs_el, xs_in);

				// > Get process cross section and error
				gen->GetXS(xs0[p][mode], xs0_err[p][mode]);

                delete gen;

		    } catch (const std::invalid_argument& err) {
				std::cerr << "PROCESS:: " << input[p].card
				          << " : Exception catched: " << err.what()
				          << std::endl;
		        exit(0);
		    } catch (const std::ios_base::failure& err) {
				std::cerr << "PROCESS:: " << input[p].card
				          << " : Exception catched: " << err.what()
				          << std::endl;
	            exit(0);
		    } catch (...) {
				std::cerr << "Exception (...) catched!" << std::endl;
				exit(0);
		    }
		} // Screening on/off
    } // Process loop

    std::cout << std::endl;
    printf("Finished in %0.1f sec \n\n", timer.ElapsedSec());

    // Loop over processes
    for (const auto& i : indices(input)) {
	printf(
	    "[%22s] DATA = %0.2E +- %0.2E (stat) +- %0.2E (syst) | MC BARE = "
	    "%0.3E (ratio = %0.2f) | MC SCRN = %0.3E (ratio = %0.2f) | <S^2> = %0.3f \n",
	    input[i].card.c_str(), input[i].value, input[i].stat,
	    input[i].syst, xs0[i][0], xs0[i][0] / input[i].value, xs0[i][1],
	    xs0[i][1] / input[i].value, xs0[i][1] / xs0[i][0]);
    }

    printf("\n");
}
