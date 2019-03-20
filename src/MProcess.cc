// Abstract process class
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <algorithm>
#include <complex>
#include <future>
#include <iostream>
#include <iterator>
#include <map>
#include <random>
#include <regex>
#include <stdexcept>
#include <vector>
#include <cmath>


// Own
#include "Graniitti/MAux.h"
#include "Graniitti/MForm.h"
#include "Graniitti/MH1.h"
#include "Graniitti/MH2.h"
#include "Graniitti/MKinematics.h"
#include "Graniitti/MMath.h"
#include "Graniitti/MMatrix.h"
#include "Graniitti/MProcess.h"
#include "Graniitti/MFactorized.h"
#include "Graniitti/MSpin.h"
#include "Graniitti/MTimer.h"
#include "Graniitti/MUserCuts.h"
#include "Graniitti/MFragment.h"


// HepMC33
#include "HepMC3/FourVector.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/Print.h"


// Libraries
#include "json.hpp"
#include "rang.hpp"


using gra::aux::indices;
using gra::math::msqrt;
using gra::math::pow2;
using gra::math::zi;
using gra::math::abs2;


namespace gra {


void MProcess::PrintSetup() const {

	std::cout << std::endl;
	std::cout << rang::style::bold
	          << "Process setup:" << rang::style::reset << std::endl
	          << std::endl;
	std::cout << "- Random seed:      " << random.GetSeed() << std::endl;
	std::cout << "- Initial state:    " << beam1.name << " " << beam2.name << std::endl;
	printf("- Beam energies:   [%0.1f %0.1f] GeV \n", lts.pbeam1.E(), lts.pbeam2.E());
	printf("- CMS energy:       %0.1f GeV\n", lts.sqrt_s);
	std::cout << "- Process:          " << PROCESS << rang::fg::green
	          << "  <" << GetProcessDescriptor(PROCESS) << ">"
	          << rang::fg::reset << std::endl
	          << std::endl;

	// Subprocess
	std::cout << rang::style::bold
	          << "Subprocess parameters:" << rang::style::reset
	          << std::endl << std::endl;
  	std::cout << "- Pomeron loop screening:  " << std::boolalpha
	          << SCREENING << std::endl;
	
	// All other than inclusive processes
    if (ProcPtr.ISTATE != "X") {
		std::cout << "- Final state:             " << DECAYMODE << std::endl;
		std::cout << "- Proton N* excitation:    " << std::boolalpha << EXCITATION;
		
		if (EXCITATION == 0){
			std::cout << rang::fg::green << "  <elastic>"
			          << rang::fg::reset << std::endl;
		}
		if (EXCITATION == 1){
			std::cout << rang::fg::green << "  <single>"
			          << rang::fg::reset << std::endl;
		}
		if (EXCITATION == 2){
			std::cout << rang::fg::green << "  <double>"
			          << rang::fg::reset << std::endl;
		}
    }
	if (FLATAMP != 0){
		std::cout << "- Flat amplitude mode:     " << FLATAMP << std::endl;
	}

	std::cout << std::endl << std::endl;
}


// Print out process lists
std::vector<std::string> MProcess::PrintProcesses() const {

	// Iterate through processes
	std::vector<std::string> processes;
	std::map<std::string, std::string>::const_iterator it = Processes.begin();

	while (it != Processes.end()) {
		processes.push_back(it->first);
		printf("%25s  =  ", it->first.c_str());
		std::cout << it->second << std::endl;
		++it;
	}
	return processes;
}

// Check if process exists
bool MProcess::ProcessExist(std::string str) const {
	if (Processes.find(str) != Processes.end()) {
		return true;
	}
	return false;
}

// Return process description string
std::string MProcess::GetProcessDescriptor(std::string str) const {
	if (!ProcessExist(str)) {
		throw std::invalid_argument("MProcess::GetProcessDescriptor: Process by name " + str + " does not exist");
	}
	return Processes.find(str)->second;
}


// Pomeron Loop Screened Amplitude
//
// ======xxxxxxxxx======>
//       *        *
//       *     -p1t + kt = qt1
//       *          *
//       * kt        x==>
//       *          *
//       *     -p2t - kt = qt2
//       *        *
// ======xxxxxxxxx======>


complex<double> MProcess::S3ScreenedAmp() {

	// Eikonal Loop Screening not on
	if (SCREENING == false) {
		return ProcPtr.GetBareAmplitude(lts);
	}

	if (ProcPtr.CHANNEL == "EL") {
		// Elastic scattering, return directly eikonalized amplitude itself
		return Eikonal.MSA.Interpolate1D(-lts.t);
	}
	
	// Proton pt vectors
	const std::vector<double> p1T = {lts.pfinal[1].Px(), lts.pfinal[1].Py()};
	const std::vector<double> p2T = {lts.pfinal[2].Px(), lts.pfinal[2].Py()};

	// Save born amplitudes
	const std::complex<double> A_0 = ProcPtr.GetBareAmplitude(lts);
	const std::vector<std::complex<double>> hamp_0 = lts.hamp;

	// Save old kinematics
	lts.pfinal_orig = lts.pfinal;


	// @@ C++11 handles multithreaded static initialization @@
	// Discretization steps
	const static double StepKT = (MEikonalNumerics::MaxLoopKT - MEikonalNumerics::MinLoopKT) /
	    						  MEikonalNumerics::NumberLoopKT;    
	const double MinPhi  = 0.0;
	const double MaxPhi  = 2.0 * gra::math::PI;
	const static double StepPhi = (MaxPhi - MinPhi) / MEikonalNumerics::NumberLoopPHI;


	// Init 2D-Simpson weight matrix (will be calculated only once, being
	// static)
	const static MMatrix<double> WSimpson =
	    gra::math::Simpson38Weight2D(MEikonalNumerics::NumberLoopPHI,
	                             MEikonalNumerics::NumberLoopKT);
    
	// NOTE N + 1, init with zero!
	MMatrix<std::complex<double>> f(MEikonalNumerics::NumberLoopPHI + 1,
									MEikonalNumerics::NumberLoopKT  + 1, 0.0);

	std::vector<MMatrix<std::complex<double>>> f_hamp(lts.hamp.size(),
	    MMatrix<std::complex<double>>(MEikonalNumerics::NumberLoopPHI + 1,
	    							  MEikonalNumerics::NumberLoopKT  + 1, 0.0));
	
	
	// 2D-integral
	//
	// \int d^2kt A_eik(kt^2) A(kt1,kt2)
	// = \int d\phi \int dkt * kt * A_el(kt^2) * A(kt1,kt2)
	//
	//        *2
	//    *3  |   *1
	//        |/phi
	//  *4----------*0/8  ---> kt --->
	//        |
	//   *5   |   *7
	//        *6
	//
	
	for (std::size_t i = 0; i < MEikonalNumerics::NumberLoopPHI+1; ++i) {
		const double phi = MinPhi + i * StepPhi;

		for (std::size_t j = 0; j < MEikonalNumerics::NumberLoopKT+1; ++j) {
			const double kt = MEikonalNumerics::MinLoopKT + j * StepKT;
			const double kt2 = gra::math::pow2(kt);

			// -------------------------------------------------------------------

			// Get screening amplitude
			// (compiler will optimize this outside the loop)
			const std::complex<double> A_eik = Eikonal.MSA.Interpolate1D(kt2);
			
			// 0. Construct loop 2D kT-vector
			const std::vector<double> kt_{kt * std::cos(phi),
			                              kt * std::sin(phi)};

			// 1. New proton pt vectors
			const std::vector<double> p1p = {p1T[0] - kt_[0],
			                                 p1T[1] - kt_[1]};
			const std::vector<double> p2p = {p2T[0] + kt_[0],
			                                 p2T[1] + kt_[1]};

			// 2. Update kinematics
			if (!LoopKinematics(p1p, p2p)){
				continue; // not valid kinematically
			}

			// 3. Get new amplitude
			const std::complex<double> A = ProcPtr.GetBareAmplitude(lts);
			// -------------------------------------------------------------------

			if (lts.hamp.size() > 1) { // Helicity amps
				for (const auto& h : indices(lts.hamp)) {
					f_hamp[h][i][j] = A_eik * lts.hamp[h] * kt; // note x kt (jacobian)
				}
			} else {
				f[i][j] = A_eik * A * kt;
			} // note x kt (jacobian)

		} // inner-loop
	} // outer-loop

	std::complex<double> A_loop(0, 0);
	std::vector<std::complex<double>> hamp_loop(lts.hamp.size(), 0.0);

	// Normalization
	const std::complex<double> norm = zi / (8.0 * gra::math::PIPI * lts.s);

	// If the process has helicity amplitudes
	if (lts.hamp.size() > 1) {
		for (const auto& h : indices(lts.hamp)) {
			hamp_loop[h] =
			    norm * gra::math::Simpson38Integral2D(f_hamp[h], WSimpson,
			                               		  StepPhi, StepKT);
		}
	} else {
		A_loop =
		    norm * gra::math::Simpson38Integral2D(f, WSimpson, StepPhi, StepKT);
	}

	// Update back to tree level kinematics
	LoopKinematics(p1T, p2T);

	// ------------------------------------------------------------
	// Final amplitude (squared)

	// Helicity amplitude process
	if (lts.hamp.size() > 1) {
		if (ProcPtr.ISTATE != "gg") { // Not Durham-QCD

			// Separate (incoherent) sum
			double A2 = 0.0;
			for (const auto& h : indices(lts.hamp)) {
				A2 += abs2(hamp_0[h] + hamp_loop[h]);
			}
			A2 /= 4; // Initial state helicity average
			return msqrt(A2);

		} else { // Durham-QCD

			// Coherent sum
			std::complex<double> A = 0.0;
			for (const auto& h : indices(lts.hamp)) {
				A += hamp_0[h] + hamp_loop[h];
			}
			// Helicity average already taken care of
			return A;
		}
	}

	// All the other amplitude processes
	return A_0 + A_loop;
}


/*
std::complex<double> MProcess::S3screeningSD(double kt2) {

  // Local discretization
  const int N = 2 * MEikonalNumerics::IntegralN; // even number
  const double step = (MEikonalNumerics::MaxBT-1)/(double) N;

  std::vector<std::complex<double>> f(N, 0.0);

  // Integral independent part
  std::vector<std::complex<double>> A =
    Beta_j(0)*pow2(Beta_i(0))*g_zij(0)*
    std::pow(s/M2,  2*alpha_i(tmin)-2)*
    std::pow(M2/s0, alpha_j(0)-1);


  // OUTERMOST INTEGRAL b_2

  // Numerical integral loop over impact parameter (b_t) space
  for (int n1 = 1; n1 <= N; ++n1) {

      const double bt1 = n*step;

      f[n-1] = Fi(bt1) * gra::aux::BESSJ0(bt1*msqrt(qt2) ) * bt1;

    // MIDDLE INTEGRAL b_3

    // Numerical integral loop over impact parameter (b_t) space
    for (int n2 = 1; n2 <= N; ++n2) {

        const double bt2 = n2*step;

        f[n-1] = Fi(bt2) * gra::aux::BESSJ0(bt2*msqrt(qt2) ) * bt2;

      // INNERMOST INTEGRAL b_1

      // Numerical integral loop over impact parameter (b_t) space
      for (int n3 = 1; n3 <= N; ++n3) {

          const double bt3 = n3*step;

          f[n-1] = Fj(bt3) * gra::aux::BESSJ0(bt3*0 ) * bt3;
      }
    }
  }


  // Composite Simpson's rule
  return A * CSIntegral(f, step) / (2*gra::math::PI * 2*gra::math::PI * 2*gra::math::PI);
}

// Calculate here the amplitudes inside
std::complex<double> MProcess::Fi(double bt) {

// Discretization of qt
const int    nqt     = 2*7000; // Even number
const double qtstep  = 4.0 / (double)(nqt);

std::vector<std::complex<double>> f(nqt, 0.0);

// Loop over
for (int n = 0; n < nqt; ++n) {

    const double qt = n*qtstep;
    double qt2 = qt*qt;
    double alpha_i = PARAM_SOFT::ALPHAP; // Pomeron slope

    // Amplitude_, and its Fourier-Bessel transform
    double A = Beta_i(qt) * std::pow(s/M2, -alphap_i*qt2) * std::exp(bprime_zij
* q2);
    f[n] = A * gra::aux::BESSJ0(bt*qt)*qt;
}
return CSIntegral(f, ktstep) / (2.0*gra::math::PI) / Beta_i(0);
}


// Calculate here the amplitudes inside
std::complex<double> MProcess::Fj(double bt) {

// Discretization of kt
const int    N_kt     = 2*7000; // Even number
const double ktstep  = 4.0 / (double)(N_kt);

std::vector<std::complex<double>> f(N_kt, 0.0);

// Loop over
for (int n = 0; n < N_kt; ++n) {

    const double kt = n*ktstep;
    double kt2 = kt*kt;
    double alpha_j = PARAM_SOFT::ALPHAP; // Pomeron slope

    // Amplitude_, and its Fourier-Bessel transform
    double A = Beta_j(kt) * std::pow(M2/s0, -alphap_j*kt2) *
std::exp(-bprime_zij * kt2);
    f[n] = A * gra::aux::BESSJ0(bt*kt)*kt;
}
return CSIntegral(f, ktstep) / (2.0*gra::math::PI) / Beta_j(0);
}
*/


// Set CMS energy and beam particle 4-vectors
void MProcess::SetInitialState(const std::vector<std::string>& beam,
							   const std::vector<double>& energy) {
	
	printf("MProcess::SetInitialState: beam: [%s, %s], energy = [%0.1f, %0.1f] \n",
		beam[0].c_str(), beam[1].c_str(), energy[0], energy[1]);

	if (beam.size() != 2){
		throw std::invalid_argument("MProcess::SetInitialState: Input BEAM vector not dim 2!");
	}
	if (energy.size() != 2){
		throw std::invalid_argument("MProcess::SetInitialState: Input ENERGY vector not dim 2!");
	}

	beam1 = PDG.FindByPDGName(beam[0]);
	beam2 = PDG.FindByPDGName(beam[1]);
	
	// Beam particle 4-momenta re-setup with safety threshold if fixed target setup
	const double E1 = std::max(energy[0], 1.001*beam1.mass);
	const double E2 = std::max(energy[1], 1.001*beam2.mass);
	
	lts.pbeam1 = M4Vec(0, 0,  msqrt(pow2(E1) - pow2(beam1.mass)), E1); // positive z-axis
	lts.pbeam2 = M4Vec(0, 0, -msqrt(pow2(E2) - pow2(beam2.mass)), E2); // negative z-axis

	// Mandelstam s
	lts.s      = (lts.pbeam1 + lts.pbeam2).M2();
	lts.sqrt_s = msqrt(lts.s);

	if (lts.sqrt_s < (beam1.mass + beam2.mass)) {
		std::string str = "MProcess::SetInitialState: Error with input CMS energy: " +
			std::to_string(lts.sqrt_s) + " GeV < initial state masses!";
		throw std::invalid_argument(str);
	}
}


// Return "flat" matrix element squared |A|^2 -> for evaluating the phase space
double MProcess::GetFlatAmp2(const gra::LORENTZSCALAR& lts) const {
	double W = 1.0;

	// Peripheral phase space: |A|^2 ~ exp(bt1) exp(bt2)
	if (FLATAMP == 1) {
		W = std::exp(PARAM_FLAT::b * lts.t1) *
		    std::exp(PARAM_FLAT::b * lts.t2);
	}
	// Peripheral phase space: |A|^2 ~ exp(bt1) exp(bt2) / sqrt(shat)
	else if (FLATAMP == 2) {
		W = std::exp(PARAM_FLAT::b * lts.t1) *
		    std::exp(PARAM_FLAT::b * lts.t2) / msqrt(lts.s_hat);
	}
	// Peripheral phase space: |A|^2 ~ exp(bt1) exp(bt2) / shat
	else if (FLATAMP == 3) {
		W = std::exp(PARAM_FLAT::b * lts.t1) *
		    std::exp(PARAM_FLAT::b * lts.t2) / lts.s_hat;
	}
	// Constant
	else if (FLATAMP == 4) {
		W = 1.0;
	} else {
		// Throw an error, unknown mode
		std::string str =
		    "MProcess::GetFlatAmp2: Unknown FLATAMP mode (|A|^2 :: 0 = off, 1 = exp{b(t1+t2)}, 2 = exp{b(t1+t2)}/M, 3 = exp{b(t1+t2)}/M^2, 4 = 1.0) : input was " + std::to_string(FLATAMP);
		throw std::invalid_argument(str);
	}
	return W;
}


// Set decaymode
void MProcess::SetDecayMode(std::string str) {

	// Clear decaytree
	lts.decaytree.clear();

	// ------------------------------------------------------------------
	// 0. Check if contains only spaces (or is empty)
	std::string check = str;
	std::string::iterator end_pos = std::remove(check.begin(), check.end(), ' ');
	check.erase(end_pos, check.end());

	// Can be empty only for <Q>-class (2->2/inclusive) processes or <F> factorized phase space
	if  (check.size() == 0 &&
		!(PROCESS.find("<Q>") != std::string::npos) && 
		!(PROCESS.find("<F>") != std::string::npos)) {
		throw std::invalid_argument("MProcess::SetDecayMode: DECAY string input is empty (can be only for <Q> or <F> class)!");
	}

	// ------------------------------------------------------------------
	// 1. Syntax check, we must have same amount of left { and right } brackets
	std::vector<std::string> words = gra::aux::Extract(str);

	unsigned int arrows = 0;
	unsigned int L = 0;
	unsigned int R = 0;
	for (std::string::size_type i = 0; i < str.size(); ++i) {
		if (str[i] == '{')
			++L;
		if (str[i] == '}')
			++R;
		if (str[i] == '>')
			++arrows;
	}
	if (L != R) {
		std::string strerr =
		    "MProces::SetDecayMode: ERROR: Decay tree has " +
		    std::to_string(L) + " left brackets { and " +
		    std::to_string(R) + " right brackets } !";
		throw std::invalid_argument(strerr);
	}
	if (arrows != L) {
		std::string strerr =
		    "MProces::SetDecayMode: ERROR: Decay tree syntax with arrows > not equal to { } brackets!";
		throw std::invalid_argument(strerr);
	}

	// ===================================================================
	// ** Read decay by recursion **
	PDG.TokenizeProcess(str, 0, lts.decaytree);
	// ===================================================================

	// Save decaymode string
	DECAYMODE = str;

	// Calculate symmetry factor
	CalculateSymmetryFactor();

	// Remind the user
	if  (lts.decaytree.size() > 2 ) {
		gra::aux::PrintWarning();
		std::cout << "Reminder: Resonance decay |matrix element 1->K|^2 is non-factorizable from the phase space for K = "
					+ std::to_string(lts.decaytree.size()) +" > 2!" << std::endl;
		std::cout << "Use &> arrow for fully separating 2->3 [+] 1->K" << std::endl;
		std::cout << std::endl;
	}
}


// Count QFT statistical (combinatorial) factor for identical final states
// 
// BASED here counting the first level final states
// Example: yyyy               -> 4!
//          pi+pi-pi+pi-       -> 2!2! (CHECK THIS!)
//          pi+pi-pi+pi-pi+pi- -> 3!3!3!
//          yy                 -> 2!
//          pi+pi-K+K-         -> 1
//
void MProcess::CalculateSymmetryFactor() {

	S_factor = 1.0; // Init with 1

	std::vector<bool > marked(lts.decaytree.size(), false);

	// Collect here multiplicities
	std::vector<int> multiplicities;

	for (const auto& i : indices(lts.decaytree)) {
	  int n = 1;
	  marked[i] = true;
	  for (const auto& j : indices(lts.decaytree)) {
	    if (i != j && lts.decaytree[i].p.pdg == lts.decaytree[j].p.pdg && marked[j] == false) {
	      marked[j] = true;
	      ++n;
	    }
	  }
	  S_factor *= gra::math::factorial(n);
	}
}


// Read resonance branching fractions
void MProcess::SetupBranching() {

	aux::PrintBar(".");
	
	if (lts.RESONANCES.size() == 0) { // Empty one is not treated
		std::cout << rang::fg::red << "MProcess::SetupBranching: lts.RESONANCES.size() == 0 !" << rang::fg::reset << std::endl;
		return;
	}
	if (lts.decaytree.size() == 0) {
		std::cout << rang::fg::red << "MProcess::SetupBranching: decaytree.size() == 0 !" << rang::fg::reset << std::endl;
		return;
	}
	
	// Get particle listing
	std::vector<int> refdecay;
	for (const auto& i : indices(lts.decaytree)) {
		refdecay.push_back(lts.decaytree[i].p.pdg);
	}

	// Setup branching fractions for resonances
	using json = nlohmann::json;

	const std::string inputfile = gra::aux::GetBasePath(2) + "/modeldata/" + "BRANCHING.json";
	const std::string data      = gra::aux::GetInputData(inputfile);
	json j;
	try {
		j = json::parse(data);
	} catch (...) {
		std::string str =
		    "MProcess::SetupBranching: Error parsing " + inputfile +
		    " (Check for extra/missing commas)";
		throw std::invalid_argument(str);
	}

	// Loop over resonances
	for (auto const& xpoint : lts.RESONANCES) {
		// Take resonance
		gra::PARAM_RES res = xpoint.second;

		// Find out if we have this resonance in tables
		std::string PDG_STR = std::to_string(res.p.pdg);

		// Did we find that resonance written down it in the Branching Tables
		bool found = j.count(PDG_STR);

		if (found == true) {

			std::cout << "MProcess::SetupBranching: Resonance PDG = " + std::to_string(res.p.pdg) << std::endl;

			// Try to find the decaymode
			bool found_decay_mode = false;
			const unsigned int MAXDECAYS = 500;
			for (std::size_t i = 0; i < MAXDECAYS; ++i) {

				const std::string SID = std::to_string(i);
				std::vector<int> decay;
				try {
					std::vector<int> vec =
					    j[PDG_STR][SID]["PDG"];
					decay = vec;
				} catch (...) {
					continue; // Did not found anything with this
					          // index ["i"]
				}

				if (std::is_permutation(refdecay.begin(),
				                        refdecay.end(),
				                        decay.begin())) {

					// Found matching decay
					res.BR = j[PDG_STR][SID]["BR"];
					res.hc.P_conservation = j[PDG_STR][SID]["P_conservation"];
					found_decay_mode = true;

					// ----------------------------------------------------------------------------
					// Construct helicity ls-coupling matrix (put default 1.0)
					// indexed directly by ls values
					MMatrix<std::complex<double>> alpha(20, 20, 1.0);
					MMatrix<bool> alpha_set(20, 20, false);

					// Allow 50 coupling parameter (l,s) pairs (more than enough)
					for (std::size_t a = 0; a < 50; ++a) {
						try {
							const unsigned int    l = j[PDG_STR][SID]["alpha_ls"][a][0];
							const unsigned int    s = j[PDG_STR][SID]["alpha_ls"][a][1];
							const double Re = j[PDG_STR][SID]["alpha_ls"][a][2];
							const double Im = j[PDG_STR][SID]["alpha_ls"][a][3];

							// Set value
							alpha[l][s] = Re + zi*Im;

							// Mark it set
							alpha_set[l][s] = true;

							std::cout << "Found ls-coupling: "
									  << "[l=" << l << ",s=" << s << "] = " << alpha[l][s] << " (Re,Im)" << std::endl;

						} catch (...) {
							continue;
						}
					}
					res.hc.alpha     = alpha;
					res.hc.alpha_set = alpha_set;
					break;
				}
			}

			if (!found_decay_mode) {
				gra::aux::PrintWarning();
				std::cout << rang::fg::red
						  << "WARNING: BRANCHING.json contains no information on this decaymode '" + DECAYMODE + "' for resonance PDG = "
							 + std::to_string(res.p.pdg) + " (setting up BR = 1.0, alpha_ls == 1.0, P_conservation = true) "
						  << rang::fg::reset << std::endl;
			  	res.BR = 1.0;
			  	res.hc.P_conservation = true;
				res.hc.alpha     = MMatrix<std::complex<double>>(20, 20, 1.0);
				res.hc.alpha_set = MMatrix<bool>(20, 20, true);
			}
			
			// ----------------------------------------------------------------------------
			// Init T-matrix
			if (lts.decaytree.size() == 2) {
				try {
					gra::spin::InitTMatrix(res.hc, res.p, lts.decaytree[0].p, lts.decaytree[1].p);
				} catch (std::invalid_argument& e) {
					throw std::invalid_argument("Problem with resonance PDG = "
						+ std::to_string(res.p.pdg) + " : " + e.what());
				}
			}
			// ----------------------------------------------------------------------------

			// ----------------------------------------------------------------------------
			// Calculate equivalent decay coupling for 2-body cases
			//
			// In 2-body case, the decay phase space and decay amplitude squared factorize

			if (lts.decaytree.size() == 2) {

				// Calculate phase space part:
				// Gamma = PS * |A_decay|^2 then with |A_decay|^2 = 1 <=> Gamma = PS

				const double amp2 = 1.0;
				const double sym  = 1.0;

				double PS = gra::kinematics::PDW2body(
				    pow2(res.p.mass),
				    pow2(lts.decaytree[0].p.mass),
				    pow2(lts.decaytree[1].p.mass), amp2, sym);

				if (std::abs(PS) < 1e-9) {

					// Try again with higher mother mass, we might
					// be trying purely off-shell decay (such as f0(980) -> K+K- )
					const unsigned int N_width = 3;
					for (std::size_t i = 1; i <= N_width; ++i) {
						PS = gra::kinematics::PDW2body(
							pow2(res.p.mass + i * res.p.width),
						    pow2(lts.decaytree[0].p.mass),
						    pow2(lts.decaytree[1].p.mass),
						    amp2, sym);

						if (PS > 1e-9) { break; }
					}
				}

				// BR \equiv Gamma/Gamma_tot = (PS * |A_decay|^2) / Gamma_tot
				// <=>
				// |A|^2 = BR * Gamma_tot / Gamma_PS and sqrt to get 'amplitude level'
				//
				// ** GeV unit can be obtained by inspecting PDW2body function **
				res.g_decay = msqrt(res.BR * res.p.width / PS);

				if (std::isnan(res.g_decay) || std::isinf(res.g_decay)) {
					std::string str =
					    "MProcess::SetupBranching:: Kinematic coupling problem to '"
					    + DECAYMODE + "' with resonance PDG "
					    + res.p.name + "(" + std::to_string(res.p.pdg) + ")" +
					    " (daughters too heavy?) (check RESONANCE, DECAYMODE and BRANCHING tables)";
					throw std::invalid_argument(str);
				}

			// 3,4,... body cases [NOT TREATED YET, phase space and matrix element do not factorize there]
			} else {
				res.g_decay = 1.0; // For the rest, put 1.0
			}

			printf("(Mass, Full width):                (%0.3E, %0.3E GeV) \n",  res.p.mass, res.p.width);
			printf("Branching ratio || Partial width:   %0.3E || %0.3E GeV \n", res.BR, res.BR * res.p.width);
			printf("=> Effective decay vertex:          %0.3E GeV \n", res.g_decay);

			// Set resonance
			lts.RESONANCES[xpoint.first] = res;

			aux::PrintBar(".");

		} else {
			std::string str =
			    "MProcess::SetupBranching:: Did not find any "
			    "branching ratio data for resonance with PDG: " +
			    std::to_string(res.p.pdg);
			throw std::invalid_argument(str);
		}
	}
}


// Get intermediate off-shell mass
void MProcess::GetOffShellMass(const gra::MDecayBranch& branch, double& mass) {

	// If zero-width particle (electron, gamma ...)
	if (branch.p.width < 1e-9) {
		mass = branch.p.mass;
		return;
	}

	const unsigned int OUTERMAXTRIAL = 1e4;
	const unsigned int INNERMAXTRIAL = 1e4;

	unsigned int outertrials = 0;
	
	while (true) {

		// We have decay daughters
		double daughter_masses = 0;
		if (branch.legs.size() > 0) {

			// Find daughter offshell masses
			for (std::size_t i = 0; i < branch.legs.size(); ++i) {

				const double M = branch.legs[i].p.mass;
				const double W = branch.legs[i].p.width;

				if (!FLATMASS2) {
					daughter_masses += std::max(0.0, random.BreitWignerRandom(M, W, OFFSHELL));
				} else {
					daughter_masses += msqrt( random.U(std::max(0.0,pow2(M - OFFSHELL*W)), std::min(lts.s, pow2(M + OFFSHELL*W)) ));
				}
			}
			const double safe_margin = 1e-4; // GeV
			daughter_masses += safe_margin;
		}

		// Pick mother offshell mass
		unsigned int innertrials = 0;
		while (true) {

			const double M = branch.p.mass;
			const double W = branch.p.width;
			
			if (!FLATMASS2) {
				mass = std::max(0.0, random.BreitWignerRandom(M, W, OFFSHELL));
			} else {
				mass = msqrt( random.U(std::max(0.0,pow2(M - OFFSHELL*W)), std::min(lts.s, pow2(M + OFFSHELL*W)) ));
			}

			++innertrials;
			if (mass > daughter_masses) {
				return; // all done
			}
			if (innertrials > INNERMAXTRIAL) {
				break; // try again with different daughter masses
			}
		}
		++outertrials;

		// Impossible
		if (outertrials > OUTERMAXTRIAL) {
			std::string str;
			for (const auto& i : indices(branch.legs)) {
				str += branch.legs[i].p.name + " ";
			}
			str = "MProcess::GetOffShellMass: Kinematically impossible decay: " + branch.p.name + " > " + str;
			throw std::invalid_argument(str);
		}
	}
}


// Proton continuum excitation (this is not yet written in the same format as central system - TBD!)
//
int MProcess::ExciteContinuum(const M4Vec& nstar, 
	const HepMC3::GenParticlePtr gen_nstar, HepMC3::GenEvent& evt, double Q2_scale, int B, int Q) {

	// Sanity check (2 x pion mass)
	if (msqrt(Q2_scale) < 0.3) {
		//return 1; // not valid
	}
	const double QProb   = 2.0/3.0; // Probability for a charged particle (isospin)
	const int OUTERTRIAL = 100;

	// Average multiplicity
    const double a0 = 0.0;
    double b0 = 2.0;
    double avgN = a0 + b0 * (1/QProb) * std::log(Q2_scale);	
    int outertrials = 0;
    
	const double M = nstar.M();

    while (true) { // OUTER
    	
	    // Draw fluctuating number of particles
	    int N = 0;
	    while (true) {

			N = random.PoissonRandom(avgN);
			//if (M < 10) { // low-mass smearing treatment
				for (std::size_t r = 0; r < 3; ++r) {
			//		N = PoissonRND(N);
				}
			//}
		 	if (nstar.M() > N*0.14) break; // Minimal mass threshold
		}

		// Boundary conditions
	    N = (N < 1)  ? 1 : N;   // At least 1
	    N = (N >= B) ? N : B;   // Baryon number
	    N = (N >= Q) ? N : Q;   // Charge conservation

		// --------------------------------------------------------------------
		// Pick particles

		std::vector<double> mass;
		std::vector<int> pdgcode;
		if (MFragment::PickParticles(M, N, B, 0, Q, mass, pdgcode, PDG, random) == 1) {
			++outertrials;
			continue;
		}

		// --------------------------------------------------------------------
		// Now decay >>
	    std::vector<bool> isstable(N, false);
	    MFragment::GetDecayStatus(pdgcode, isstable);
	    
	    // Decay p
	    std::vector<M4Vec> products;

		double q = std::pow(N, 0.11/3);    // Powerlaw high-pt slope
		double T = 0.065;                   // Temperature
		if (etree.size() > 0) {            // "Temperature" rises ~ 1/impact parameter squared
			T *= std::pow(1/(bt*bt), 0.10);
		}
		const double maxpt = 15.0;         // Maximum Pt per particle
		const double W = MFragment::TubeFragment(nstar, M, mass, products, q, T, maxpt, random);
		
	    if (W  <= 0 ) {
	    	++outertrials;
	    	if (outertrials > OUTERTRIAL) {
	    		//std::cout << "MProcess::ExciteContinuum: Outer loop failure!" << std::endl;
	    		return 1; // too many failures
	    	} else {
	    		continue; // try again!
	    	}
	    } else {

	    	// Re-check kinematics
	    	bool fail = false;
	    	for (const auto& i : indices(products)) {
	    		if ( std::isnan(products[i].Rap()) ) {
	    			//std::cout << rang::fg::yellow << "MProcess::ExciteContinuum: Kinematics failure!" << rang::style::reset << std::endl;
			    	++outertrials;
			    	if (outertrials > OUTERTRIAL) {
			    		//std::cout << "MProcess::ExciteContinuum: Outer loop failure!" << std::endl;
			    		return 1; // too many failures
			    	}
			    	fail = true;
			    	break;
	    		}
	    	}
	    	if (fail) {
	    		continue;
	    	}
		    // Save event output
		    BranchForwardSystem(products, pdgcode, isstable, gen_nstar, evt);
		    break;
		}
	}
	return 0;
}


// This is used with continuum excitation
//
void MProcess::BranchForwardSystem(const std::vector<M4Vec>& p,
								   const std::vector<int>& pdgcode, const std::vector<bool>& isstable,
								   const HepMC3::GenParticlePtr gen_nstar, HepMC3::GenEvent& evt) {

	// Create N* decay vertex
	HepMC3::GenVertexPtr v = std::make_shared<HepMC3::GenVertex>();
	evt.add_vertex(v);

	// Add N* in to the decay vertex
	v->add_particle_in(gen_nstar);

	// Add decay products to the decay vertex
	std::vector<HepMC3::GenParticlePtr> gen_part;
	for (const auto& i : indices(p)) {
		
		int status = (isstable[i] == true) ? PDG::PDG_STABLE : PDG::PDG_DECAY;
		HepMC3::GenParticlePtr gen_p =
		    std::make_shared<HepMC3::GenParticle>(gra::aux::M4Vec2HepMC3(p[i]),
		                                         pdgcode[i], status);
		gen_part.push_back(gen_p);

		// Add to the vertex
		v->add_particle_out(gen_p);
	}

	// --------------------------------------------------------------
	// Do the (pi0 -> yy) decays

	for (const auto& i : indices(pdgcode)) {

		if (pdgcode[i] == PDG::PDG_pi0) { // pi0

			// Decay
			std::vector<double> m = {0.0, 0.0};
			std::vector<M4Vec> pgamma;
			gra::kinematics::TwoBodyPhaseSpace(p[i], p[i].M(), m, pgamma, random);

			// Gammas
			HepMC3::GenParticlePtr gen_p1 =
			    std::make_shared<HepMC3::GenParticle>(gra::aux::M4Vec2HepMC3(pgamma[0]), PDG::PDG_gamma, PDG::PDG_STABLE);
			HepMC3::GenParticlePtr gen_p2 =
			    std::make_shared<HepMC3::GenParticle>(gra::aux::M4Vec2HepMC3(pgamma[1]), PDG::PDG_gamma, PDG::PDG_STABLE);

			// Create a decay vertex and put particles in and out
			HepMC3::GenVertexPtr vert = std::make_shared<HepMC3::GenVertex>();
			evt.add_vertex(vert);

			// In
			vert->add_particle_in(gen_part[i]);

			// Out
			vert->add_particle_out(gen_p1);
			vert->add_particle_out(gen_p2);
		}
	}
}


// Proton low mass excitation (this function is now written in same format as central system)
//
void MProcess::ExciteNstar(const M4Vec& nstar, gra::MDecayBranch& forward) {

	// Find random decaymode
	std::vector<int> pdgcode;
	MFragment::NstarDecayTable(nstar.M(), pdgcode, random);

	// Get corresponding masses
	std::vector<gra::MParticle> p(pdgcode.size());
	std::vector<double> mass(pdgcode.size(), 0.0);

	for (const auto& i : indices(pdgcode)) {
		p[i]    = PDG.FindByPDG(pdgcode[i]);
		mass[i] = p[i].mass;
	}

	// 4-momenta
	std::vector<M4Vec> p4;

	// Do the 2 or 3-body isotropic decay
	if (mass.size() == 2) {
		gra::kinematics::TwoBodyPhaseSpace(nstar, nstar.M(), mass, p4, random);
	}
	if (mass.size() == 3) {
		const bool UNWEIGHT = true;
		gra::kinematics::ThreeBodyPhaseSpace(nstar, nstar.M(), mass, p4, UNWEIGHT, random);
	}

	// ------------------------------------------------------------------

	// Construct decaytree
	forward       = gra::MDecayBranch(); // Initialize!!
	forward.p.pdg = PDG::PDG_NSTAR;
	forward.p4    = nstar;

	// Daughters
	forward.legs.resize(p.size());
	forward.depth = 0;

	for (const auto& i : indices(p)) {

		// Decay particle
		gra::MDecayBranch branch;
		branch.p  = p[i];
		branch.p4 = p4[i];

		// Treat pi0 -> yy
		if (branch.p.pdg == PDG::PDG_pi0) {

			// Decay
			std::vector<double> mgamma = {0.0, 0.0};
			std::vector<M4Vec> pgamma;
			gra::kinematics::TwoBodyPhaseSpace(branch.p4, branch.p4.M(), mgamma, pgamma, random);

			// Add gamma legs
			branch.legs.resize(2);
			for (std::size_t k = 0; k < 2; ++k) {
				gra::MDecayBranch gammabranch;
				gammabranch.p  = PDG.FindByPDG(22);
				gammabranch.p4 = pgamma[k];

				branch.legs[k] = gammabranch;
				branch.legs[k].depth = 2;
			}
		}

		// Add leg
		forward.legs[i] = branch;
		forward.legs[i].depth = 1;
	}
}

bool MProcess::VetoCuts() const {

	bool ok = true;

	// Veto cuts
	if (vetocuts.active == true) {

		// Forward system
		FindVetoCuts(lts.decayforward1, ok);
		FindVetoCuts(lts.decayforward2, ok);

		// Central system
		for (const auto& i : indices(lts.decaytree)) {
			FindVetoCuts(lts.decaytree[i], ok);
		}
	}
	return ok;
}


// Common cuts for MResonance and MContinuum classes
bool MProcess::CommonCuts() const {

	bool ok = true;

	// Fiducial cuts
	if (fcuts.active == true) {

		// Check custom user cuts (do not substitute to kinematics_ok = UserCut...)
		if (!UserCut(USERCUTS, lts)) {
			return false; // not fine
		}

		// Check forward system variables
		if (CID != "P") { // Collinear class does not support these
					
			if (std::abs(lts.t1) >= fcuts.forward_t_min && std::abs(lts.t1) <= fcuts.forward_t_max &&
			    std::abs(lts.t2) >= fcuts.forward_t_min && std::abs(lts.t2) <= fcuts.forward_t_max) {
					// fine
			} else { return false; }

			if (lts.excite1) {
			    if (lts.pfinal[1].M() >= fcuts.forward_M_min && lts.pfinal[1].M() <= fcuts.forward_M_max) {
					// fine
			    } else { return false; }
			}

			if (lts.excite2) {
				if (lts.pfinal[2].M() >= fcuts.forward_M_min && lts.pfinal[2].M() <= fcuts.forward_M_max) {
					// fine
				} else { return false; }
			}
		}

		// Check system variables
		if (msqrt(lts.m2) >= fcuts.M_min && msqrt(lts.m2) <= fcuts.M_max &&
		    lts.Y  >= fcuts.Y_min  && lts.Y  <= fcuts.Y_max &&
		    lts.Pt >= fcuts.Pt_min && lts.Pt <= fcuts.Pt_max) {	
			// fine, do not touch
		} else {
			return false; // not fine
		}
		
		// Check fiducial cuts of the central final state particles
		for (const auto& i : indices(lts.decaytree)) {
			FindDecayCuts(lts.decaytree[i], ok);
		}
	}

	return ok;
}


// Recursive find out final daughters and check if they pass fiducial cuts
void MProcess::FindDecayCuts(const gra::MDecayBranch& branch, bool& ok) const {

	// We are at the end -> must be a final state
	if (branch.legs.size() == 0) {

		// Check cuts
		if (branch.p4.Pt()  >= fcuts.pt_min &&
		    branch.p4.Pt()  <= fcuts.pt_max &&
		    branch.p4.Et()  >= fcuts.Et_min &&
		    branch.p4.Et()  <= fcuts.Et_max &&
		    branch.p4.Eta() >= fcuts.eta_min &&
		    branch.p4.Eta() <= fcuts.eta_max &&
		    branch.p4.Rap() >= fcuts.rap_min &&
		    branch.p4.Rap() <= fcuts.rap_max) {
			// Event passes cuts
		} else {
			ok = false; // does not pass
		}
	}

	// ** RECURSION here **
	for (const auto& i : indices(branch.legs)) {
		FindDecayCuts(branch.legs[i], ok);
	}
}


// Recursive find out final daughters and check if they trigget veto cuts
void MProcess::FindVetoCuts(const gra::MDecayBranch& branch, bool& ok) const {

	// We are at the end -> must be a final state
	if (branch.legs.size() == 0) {
		
		// Loop over veto domains
		for (const auto& i : indices(vetocuts.cuts)) {

			// Check cuts
			if (branch.p4.Pt()  >= vetocuts.cuts[i].pt_min &&
			    branch.p4.Pt()  <= vetocuts.cuts[i].pt_max &&
			    branch.p4.Eta() >= vetocuts.cuts[i].eta_min &&
			    branch.p4.Eta() <= vetocuts.cuts[i].eta_max) {

				ok = false; // VETO, does not pass
			} else {
				// does not trigger veto, do nothing
			}
		}
	}

	// ** RECURSION here **
	for (const auto& i : indices(branch.legs)) {
		FindVetoCuts(branch.legs[i], ok);
	}
}


// Recursive function to plot out the decay tree
void MProcess::PrintDecayTree(const gra::MDecayBranch& branch) const {
	std::string spaces(branch.depth * 2, ' '); // Give empty space
	std::string lines = spaces + "|--->  ";

	if ((branch.depth + 1) % 3 == 0) {
		std::cout << rang::fg::blue << lines << rang::fg::reset;
	} else if ((branch.depth + 1) % 2 == 0) {
		std::cout << rang::fg::yellow << lines << rang::fg::reset;
	} else {
		std::cout << rang::fg::green << lines << rang::fg::reset;
	}

	printf(
	    "M = %0.2E, W = %0.2E (GeV) [<tau> = %0.1E s] PDG = [%s%d | %s] "
	    "[Q=%s, J=%s] \n",
	    branch.p.mass, branch.p.width, branch.p.tau,
	    branch.p.pdg > 0 ? " " : "", branch.p.pdg, branch.p.name.c_str(),
	    gra::aux::Charge3XtoString(branch.p.chargeX3).c_str(),
	    gra::aux::Spin2XtoString(branch.p.spinX2).c_str());

	// ** RECURSION here **
	for (const auto& i : indices(branch.legs)) {
		PrintDecayTree(branch.legs[i]);
	}
}


// Recursive function to plot out the decay tree
void MProcess::PrintPhaseSpace(const gra::MDecayBranch& branch) const {
	std::string spaces(branch.depth * 2, ' '); // Give empty space
	std::string lines = spaces + "|---";
	
	if (branch.W.GetW() > 0) { // There is decay information
		printf("%s> \t {1->%lu LIPS}:  %0.3E +- %0.3E \n",
			lines.c_str(), branch.legs.size(), branch.W.Integral(), branch.W.IntegralError());
	}
	for (const auto& i : indices(branch.legs)) {
		PrintPhaseSpace(branch.legs[i]);
	}
}


// Recursive decay tree kinematics (called event by event from inhereting
// classes)
bool MProcess::ConstructDecayKinematics(gra::MDecayBranch& branch) {
	
	// This leg has any daughters
	if (branch.legs.size() > 0) {
		std::vector<double> masses;
		
		// Generate decay product masses
		for (const auto& i : indices(branch.legs)) {

			GetOffShellMass(branch.legs[i], branch.legs[i].m_offshell);
			masses.push_back(branch.legs[i].m_offshell);
		}
		std::vector<M4Vec> p;
		
		// For now, keep always unweighted
		const bool UNWEIGHT = true;
		
		// 2-body
		gra::kinematics::MCW w;
		if (branch.legs.size() == 2) {
			w = gra::kinematics::TwoBodyPhaseSpace(
				branch.p4, branch.p4.M(), masses, p, random);
		// 3-body
		} else if (branch.legs.size() == 3) {
			w = gra::kinematics::ThreeBodyPhaseSpace(
			    branch.p4, branch.p4.M(), masses, p, UNWEIGHT, random);
		// N-body
		} else {
			w = gra::kinematics::NBodyPhaseSpace(
			    branch.p4, branch.p4.M(), masses, p, UNWEIGHT, random);
		}
		if (w.GetW() < 0) {
			std::string str = "MProcess::ConstructDecayKinematics: Fatal error: Weight < 0 (Check your decay tree)";
			std::cout << str << std::endl;
			return false;
			//throw std::invalid_argument();
		}
		
		// Collect weight
		branch.W += w;
		
		// Collect decay product 4-momenta
		for (const auto& i : indices(branch.legs)) {
			branch.legs[i].p4 = p[i];
		}
		
		//  ** Now the INNER recursion **
		for (const auto& i : indices(branch.legs)) {
			if (!ConstructDecayKinematics(branch.legs[i])) { return false; }
		}
	}
	return true;
}


// Recursively add final states to the event structure
void MProcess::WriteDecayKinematics(gra::MDecayBranch& branch,
                                    HepMC3::GenParticlePtr& mother,
                                    HepMC3::GenEvent& evt) {
	// This particle has daughters
	if (branch.legs.size() > 0) {

		// Create new vertex with decay 4-position
		HepMC3::GenVertexPtr vertex =
		    std::make_shared<HepMC3::GenVertex>(gra::aux::M4Vec2HepMC3(branch.decay_position));
		evt.add_vertex(vertex);

		// The decaying particle
		vertex->add_particle_in(mother);

		// Add daughters
		for (const auto& i : indices(branch.legs)) {
			const int STATE = (branch.legs[i].legs.size() > 0)
			                ? PDG::PDG_DECAY : PDG::PDG_STABLE;
			                
			// ADD HERE THE ctau > 1.0 cm definition for the status
			// code [TBD]

			HepMC3::GenParticlePtr particle =
			    std::make_shared<HepMC3::GenParticle>(gra::aux::M4Vec2HepMC3(branch.legs[i].p4), branch.legs[i].p.pdg, STATE);
			vertex->add_particle_out(particle);

			// ** RECURSION here **
			WriteDecayKinematics(branch.legs[i], particle, evt);
		}
	}
}


// Forward excitation mass sampling
void MProcess::SampleForwardMasses(std::vector<double>& mvec, const std::vector<double>& randvec) {

	mvec = {beam1.mass, beam2.mass};
	
	lts.excite1 = false;
	lts.excite2 = false;

	M2_f_min = pow2(1.07); // proton + pion threshold
	M2_f_max = gcuts.XI_max * lts.s;
	
	if      (EXCITATION == 1) {
		const double mforward = msqrt( M2_f_min + (M2_f_max - M2_f_min) * randvec[0] );
		if (random.U(0,1) < 0.5) {
			mvec[0]     = mforward;
			lts.excite1 = true;
		} else {
			mvec[1]     = mforward;
			lts.excite2 = true;
		}
	}
	else if (EXCITATION == 2) {
		mvec[0]     = msqrt( M2_f_min + (M2_f_max - M2_f_min) * randvec[0] );
		mvec[1]     = msqrt( M2_f_min + (M2_f_max - M2_f_min) * randvec[1] );
		lts.excite1 = true;
		lts.excite2 = true;
	}
}


// Print fiducial cuts set by user
void MProcess::PrintFiducialCuts() const {

	gra::aux::PrintBar("-");
	std::cout << std::endl;
	std::cout << rang::style::bold;
	std::cout << "Fiducial cuts:" << std::endl << std::endl;
	std::cout << rang::style::reset;
	
	if (fcuts.active == true) {

		std::cout << "Central final states" << std::endl;
		printf("- Eta   [min, max] = [%0.2f, %0.2f] \n",     fcuts.eta_min, fcuts.eta_max);
		printf("- Pt    [min, max] = [%0.2f, %0.2f] GeV \n", fcuts.pt_min,  fcuts.pt_max);
		printf("- Et    [min, max] = [%0.2f, %0.2f] GeV \n", fcuts.Et_min,  fcuts.Et_max);
		printf("- Rap   [min, max] = [%0.2f, %0.2f] \n",     fcuts.rap_min, fcuts.rap_max);
		std::cout << std::endl;

		std::cout << "Central system" << std::endl;
		printf("- M     [min, max] = [%0.2f, %0.2f] GeV \n", fcuts.M_min,  fcuts.M_max);
		printf("- Pt    [min, max] = [%0.2f, %0.2f] GeV \n", fcuts.Pt_min, fcuts.Pt_max);
		printf("- Rap   [min, max] = [%0.2f, %0.2f] \n",     fcuts.Y_min,  fcuts.Y_max);
		std::cout << std::endl;

		std::cout << "Forward kinematics" << std::endl;
		printf("- M     [min, max] = [%0.2f, %0.2f] GeV \n",   fcuts.forward_M_min, fcuts.forward_M_max);
		printf("- |t_i| [min, max] = [%0.2f, %0.2f] GeV^2 \n", fcuts.forward_t_min, fcuts.forward_t_max);
		
	} else {
		std::cout << "- Not active" << std::endl;
	}
	gra::aux::PrintBar("-");
    std::cout << std::endl;
	std::cout << rang::style::bold;
	std::cout << "Extra custom fiducial cuts:" << std::endl << std::endl;
	std::cout << rang::style::reset;
	std::cout << rang::fg::red;
	if (USERCUTS != 0) {
		printf("- Active ID = %d (see MUserCuts.cc) \n", USERCUTS);
	} else {
		std::cout << "- Not active" << std::endl;
	}
	std::cout << rang::fg::reset;

	gra::aux::PrintBar("-");
	std::cout << std::endl;
	std::cout << rang::style::bold;
	std::cout << "VETO cuts:" << std::endl << std::endl;
	std::cout << rang::style::reset;

	if (vetocuts.active == true) {
		for (std::size_t i = 0; i < vetocuts.cuts.size(); ++i) {
			std::cout << "DOMAIN: " << i << std::endl;
			printf("- Eta   [min, max] = [%0.2f, %0.2f] \n", vetocuts.cuts[i].eta_min, vetocuts.cuts[i].eta_max);
			printf("- Pt    [min, max] = [%0.2f, %0.2f] GeV \n", vetocuts.cuts[i].pt_min, vetocuts.cuts[i].pt_max);
			std::cout << std::endl;
		}
	} else {
		std::cout << "- Not active" << std::endl;
	}

	gra::aux::PrintBar("-");
	std::cout << std::endl;
	std::cout << rang::style::bold;
	std::cout << "Central system decay tree:" << std::endl << std::endl;
	std::cout << rang::style::reset;
	// Print out decaytree recursively
	for (std::size_t i = 0; i < lts.decaytree.size(); ++i) {
		PrintDecayTree(lts.decaytree[i]);
	}

	std::cout << std::endl;
	printf("Final state symmetry factor (1/S = 1/%0.0f) applied at cross section level \n", S_factor);

	gra::aux::PrintBar("-");
	std::cout << std::endl;
}


/*
bool MProcess::B5Getgra::LORENTZSCALARs() {

  // 5 Lorentz scalars:
  //
  //      t1
  // ----------->              lts.pf[1]
  //      |  s10 (s1)
  //  s   |-------> shat (s0)  lts.pf[0]
  //      |  s20 (s2)
  // ----------->              lts.pf[2]
  //      t2

  // s-type scalars
  lts.ss[1][0] = (lts.pfinal[1] + lts.pfinal[0]).M2();
  lts.ss[2][0] = (lts.pfinal[2] + lts.pfinal[0]).M2();
  lts.s1 = lts.ss[1][0];
  lts.s2 = lts.ss[2][0];
  
  lts.m2    = lts.pfinal[0].M2();
  lts.s_hat = lts.m2;
  lts.Y     = lts.pfinal[0].Rap();
  lts.Pt    = lts.pfinal[0].Pt();
  
  // t-type scalars
  lts.t1 = (lts.pbeam1 - lts.pfinal[1]).M2();
  lts.t2 = (lts.pbeam2 - lts.pfinal[2]).M2();

  // For 2-body central process
  lts.t_hat = (lts.pbeam1 - lts.pfinal[1] - lts.decaytree[0].p4).M2(  ); // that
  lts.u_hat = (lts.pbeam1 - lts.pfinal[1] - lts.decaytree[1].p4).M2(  ); // uhat
  
  // Bjorken-x
  lts.x1 = (1.0 - lts.pfinal[1].Pz() / lts.pbeam1.Pz() );
  lts.x2 = (1.0 - lts.pfinal[2].Pz() / lts.pbeam2.Pz() );

  // Propagator vectors
  lts.q1 = lts.pbeam1 - lts.pfinal[1];
  lts.q2 = lts.pbeam2 - lts.pfinal[2];

  // Propagator pt
  lts.qt1 = lts.q1.Pt();
  lts.qt2 = lts.q2.Pt();

  // Test scalars
  if ( lts.s1 > lts.s || lts.s1 < 0) return false;
  if ( lts.s2 > lts.s || lts.s2 < 0) return false;
  if ( lts.m2  > lts.s) return false;
  if ( lts.t1  > 0) return false;
  if ( lts.t2  > 0) return false;

  return true;
}
*/

// Build and check Lorentz scalars
// Input as the number of final states
bool MProcess::GetLorentzScalars(unsigned int Nf) {
	
	// Example Nf = 4 gives 8 Lorentz scalars:
	//
	//        t1
	// ---------------->                 lts.pfinal[1]
	//         |   s13
	//         |------------->           lts.pfinal[3]
	//         |
	// s   that,uhat       shat (s34)    lts.pfinal[0]
	//         |
	//         |------------->           lts.pfinal[4]
	//         |   s24
	// ---------------->                 lts.pfinal[2]
	//        t2

	// s-type Lorentz scalars -->
	const int offset = 3; // central system indexing starts from offset

	// Upper right triangle
	for (std::size_t i = 0; i <= Nf; ++i) { // start at 0
		for (std::size_t j = 0; j <= Nf; ++j) {
			if (i < j){
				lts.ss[i][j] = (lts.pfinal[i] + lts.pfinal[j]).M2();
			}
			if (lts.ss[i][j] < 0) { return false; }
		}
	}
	// Copy the i<->j permutated to the left bottom triangle
	// (faster than calculating twice)
	for (std::size_t i = 0; i <= Nf; ++i) { // start at 0
		for (std::size_t j = i; j <= Nf; ++j) {
			if (i != j) {
				lts.ss[j][i] = lts.ss[i][j];
			}
		}
	}
	
	lts.s1 = lts.ss[1][0]; // Sub invariants w.r.t central system
	lts.s2 = lts.ss[2][0];
	if (lts.s1 < 0 || lts.s2 < 0) { return false; }

	// ------------------------------------------------------------------

	// t-type Lorentz scalars -->

	// Propagator vectors
	lts.q1 = lts.pbeam1 - lts.pfinal[1];
	lts.q2 = lts.pbeam2 - lts.pfinal[2];

	lts.t1 = lts.q1.M2();
	lts.t2 = lts.q2.M2();

	if (lts.t1 > 0 || lts.t2 > 0) { return false; }

	// For 2-body central processes
	if (lts.decaytree.size() == 2) {
		lts.t_hat = (lts.q1 - lts.decaytree[0].p4).M2(); // note q1 on both
		lts.u_hat = (lts.q1 - lts.decaytree[1].p4).M2(); // in t and u!
	}
	
	for (std::size_t i = offset; i <= Nf; ++i) {
		lts.tt_1[i] = (lts.q1 - lts.pfinal[i]).M2();
	}
	for (std::size_t i = offset; i <= Nf; ++i) {
		for (std::size_t j = offset; j <= Nf; ++j) {
			lts.tt_xy[i][j] = (lts.q1 - lts.pfinal[i] - lts.pfinal[j]).M2();
		}
	}
	for (std::size_t i = offset; i <= Nf; ++i) {
		lts.tt_2[i] = (lts.q2 - lts.pfinal[i]).M2();
	}
	
	// Fractional longitudinal momentum loss [0,1]
	lts.x1 = (1 - lts.pfinal[1].Pz() / lts.pbeam1.Pz());
	lts.x2 = (1 - lts.pfinal[2].Pz() / lts.pbeam2.Pz());

	// Bjorken-x [0,1] (this Lorentz invariant expression
	// gives 1 always for elastic central production forward leg)
	lts.xbj1 = lts.t1 / (2 * (lts.pbeam1 * lts.q1));
	lts.xbj2 = lts.t2 / (2 * (lts.pbeam2 * lts.q2));

	// Propagator pt
	lts.qt1 = lts.q1.Pt();
	lts.qt2 = lts.q2.Pt();
	
	// Often used system variables
	lts.m2    = lts.pfinal[0].M2();
	lts.s_hat = lts.m2;
	lts.Y     = lts.pfinal[0].Rap();
	lts.Pt    = lts.pfinal[0].Pt();

	// DEBUG
	/*
	std::cout << std::endl;
	printf("sqrt[ lts.m2 ]  = %0.9f \n", msqrt(lts.m2));
	printf("lts.Y           = %0.9f \n", lts.Y);
	printf("lts.Pt          = %0.9f \n", lts.Pt);
	printf("lts.t1          = %0.9f \n", lts.t1);
	printf("lts.t2          = %0.9f \n", lts.t2);
	printf("lts.t_hat - decaytree[0].M2() = %0.9f \n", lts.t_hat - lts.decaytree[0].p4.M2());
	printf("lts.u_hat - decaytree[0].M2() = %0.9f \n", lts.u_hat - lts.decaytree[0].p4.M2());
	
	for (std::size_t i = 0; i < randvec.size(); ++i) {
		printf("randvec[%d] = %0.9f \n", i, randvec[i]);
	}
	for (std::size_t i = 0; i <= lts.decaytree.size() + 2; ++i) {
		for (std::size_t j = 0; j <= lts.decaytree.size() + 2; ++j) {
			printf("sqrt[ ss[%lu][%lu] ] = %0.9f \n", i, j, msqrt(lts.ss[i][j]));
		}
	}
	for (std::size_t i = 0; i <= lts.decaytree.size() + 2; ++i) {
		for (std::size_t j = 0; j <= lts.decaytree.size() + 2; ++j) {
			printf("tt_xy[%lu][%lu] = %0.9f \n", i, j, lts.tt_xy[i][j]);
		}
	}
	std::cout << std::endl;
	*/

	return true;
}

bool MProcess::CommonRecord(HepMC3::GenEvent& evt) {

	// Initial state protons (4-momentum, pdg-id, status code)
	HepMC3::GenParticlePtr gen_p1 =
	    std::make_shared<HepMC3::GenParticle>(gra::aux::M4Vec2HepMC3(lts.pbeam1), beam1.pdg, PDG::PDG_BEAM);
	HepMC3::GenParticlePtr gen_p2 =
	    std::make_shared<HepMC3::GenParticle>(gra::aux::M4Vec2HepMC3(lts.pbeam2), beam2.pdg, PDG::PDG_BEAM);

	// Final state protons/N*
	int PDG_ID1 = beam1.pdg;
	int PDG_ID2 = beam2.pdg;

	int PDG_status1 = PDG::PDG_STABLE;
	int PDG_status2 = PDG::PDG_STABLE;

	if (lts.excite1 == true) {
		PDG_ID1 = PDG::PDG_NSTAR;
		PDG_status1 = PDG::PDG_INTERMEDIATE;
	}
	if (lts.excite2 == true) {
		PDG_ID2 = PDG::PDG_NSTAR;
		PDG_status2 = PDG::PDG_INTERMEDIATE;
	}

	HepMC3::GenParticlePtr gen_p1f = std::make_shared<HepMC3::GenParticle>(
	    gra::aux::M4Vec2HepMC3(lts.pfinal[1]), PDG_ID1, PDG_status1);
	HepMC3::GenParticlePtr gen_p2f = std::make_shared<HepMC3::GenParticle>(
	    gra::aux::M4Vec2HepMC3(lts.pfinal[2]), PDG_ID2, PDG_status2);

	// -------------------------------------------------------------------

	// Propagator 1 and 2
	HepMC3::GenParticlePtr gen_q1 = std::make_shared<HepMC3::GenParticle>(
	    gra::aux::M4Vec2HepMC3(lts.q1), PDG::PDG_pomeron, PDG::PDG_INTERMEDIATE);
	HepMC3::GenParticlePtr gen_q2 = std::make_shared<HepMC3::GenParticle>(
	    gra::aux::M4Vec2HepMC3(lts.q2), PDG::PDG_pomeron, PDG::PDG_INTERMEDIATE);

	// -------------------------------------------------------------------

	// Central system / resonance
	HepMC3::GenParticlePtr gen_q = std::make_shared<HepMC3::GenParticle>(
	    gra::aux::M4Vec2HepMC3(lts.pfinal[0]), PDG::PDG_system, PDG::PDG_INTERMEDIATE);

	// ====================================================================
	// Construct vertices

	// Upper proton-pomeron-proton
	HepMC3::GenVertexPtr v1 = std::make_shared<HepMC3::GenVertex>();
	v1->add_particle_in(gen_p1);
	v1->add_particle_out(gen_p1f);
	v1->add_particle_out(gen_q1);

	// Lower proton-pomeron-proton
	HepMC3::GenVertexPtr v2 = std::make_shared<HepMC3::GenVertex>();
	v2->add_particle_in(gen_p2);
	v2->add_particle_out(gen_p2f);
	v2->add_particle_out(gen_q2);

	// Pomeron-Pomeron-System vertex
	HepMC3::GenVertexPtr v3 = std::make_shared<HepMC3::GenVertex>();
	v3->add_particle_in(gen_q1);
	v3->add_particle_in(gen_q2);
	v3->add_particle_out(gen_q);

	evt.add_vertex(v1);
	evt.add_vertex(v2);
	evt.add_vertex(v3);

	// ====================================================================
	// System->Decay products vertex

	HepMC3::GenVertexPtr v4 = std::make_shared<HepMC3::GenVertex>();
	evt.add_vertex(v4);

	// Add resonance in
	v4->add_particle_in(gen_q);

	// Add direct daughters
	for (const auto& i : indices(lts.decaytree)) {
		const int STATE = (lts.decaytree[i].legs.size() > 0) ? PDG::PDG_DECAY : PDG::PDG_STABLE;
		
			// TBD: ADD HERE THE ctau > 1.0 cm definition for the status code

			if (STATE == PDG::PDG_DECAY) {
			// ----------------------------------------------------------
			// Pick exponential lifetime in the rest frame
			const double tau = random.ExpRandom(1.0 / lts.decaytree[i].p.tau);

			// Get decay vertex coordinates in the lab
			const double MM = 1E3; // meters to millimeters
			lts.decaytree[i].decay_position = lts.decaytree[i].p4.PropagatePosition(tau, MM);
			// ----------------------------------------------------------
			}

		HepMC3::GenParticlePtr particle = std::make_shared<HepMC3::GenParticle>(gra::aux::M4Vec2HepMC3(lts.decaytree[i].p4), lts.decaytree[i].p.pdg, STATE);
		v4->add_particle_out(particle);

		WriteDecayKinematics(lts.decaytree[i], particle, evt);
	}

	// ----------------------------------------------------------------
	// Upper proton excitation
	if (lts.excite1 == true) {

		// Create vertex
		HepMC3::GenVertexPtr v1X = std::make_shared<HepMC3::GenVertex>();
		evt.add_vertex(v1X);

		// Add N* in
		v1X->add_particle_in(gen_p1f);	

		// Add direct daughters
		for (const auto& i : indices(lts.decayforward1.legs)) {

			const int STATE = (lts.decayforward1.legs[i].legs.size() > 0) ? PDG::PDG_DECAY : PDG::PDG_STABLE;
			// TBD: ADD HERE THE ctau > 1.0 cm definition for the status code

			HepMC3::GenParticlePtr particle = std::make_shared<HepMC3::GenParticle>(
				gra::aux::M4Vec2HepMC3(lts.decayforward1.legs[i].p4), lts.decayforward1.legs[i].p.pdg, STATE);
			v1X->add_particle_out(particle);

			WriteDecayKinematics(lts.decayforward1.legs[i], particle, evt);
		}
	}

	// Lower proton excitation
	if (lts.excite2 == true) {

		// Create vertex
		HepMC3::GenVertexPtr v2X = std::make_shared<HepMC3::GenVertex>();
		evt.add_vertex(v2X);

		// Add N* in
		v2X->add_particle_in(gen_p2f);

		// Add direct daughters
		for (const auto& i : indices(lts.decayforward2.legs)) {
			
			const int STATE = (lts.decayforward2.legs[i].legs.size() > 0) ? PDG::PDG_DECAY : PDG::PDG_STABLE;
			// TBD: ADD HERE THE ctau > 1.0 cm definition for the status code
			
			HepMC3::GenParticlePtr particle = std::make_shared<HepMC3::GenParticle>(
				gra::aux::M4Vec2HepMC3(lts.decayforward2.legs[i].p4), lts.decayforward2.legs[i].p.pdg, STATE);
			v2X->add_particle_out(particle);
			
			WriteDecayKinematics(lts.decayforward2.legs[i], particle, evt);
		}
	}

	return true;
}

// Set process
void MProcess::SetProcess(std::string& str, const std::vector<aux::OneCMD>& syntax) {

	// SET IT HERE!
	PROCESS = str;
		
	// Call this always first
	PDG.ReadParticleData(gra::aux::GetBasePath(2) + "/modeldata/mass_width_2018.mcd");
	
	// @SYNTAX Read and set new PDG input
	for (const auto& i : indices(syntax)) {
		if (syntax[i].id.find("PDG") != std::string::npos) {
			std::size_t left  = syntax[i].id.find("[");
			std::size_t right = syntax[i].id.find("]");
			if (left == std::string::npos || right == std::string::npos) {
				throw std::invalid_argument("MProcess::SetProcess: invalid @PDG[]{} syntax");
			}
			// Read pdg id / number
			const int pdg = std::stoi( syntax[i].id.substr(left+1, right-left-1) );
			
			MParticle p = PDG.FindByPDG(pdg);

			// Try to find the particle from PDG table
			try {
				MParticle p = PDG.FindByPDG(pdg);
				MParticle p_anti;
				bool found_anti = false;
				
				// Check if it has an anti-particle
				try {
					p_anti = PDG.FindByPDG(-pdg);
					found_anti = true;
				} catch (...) {} // do nothing
				
				// Set new properties
				for (auto& [key,val] : syntax[i].arg) {
					if (key == "M") {
						p.mass = std::stod(std::any_cast<std::string>(val));
						if (found_anti) {
							p_anti.mass = p.mass;
						}
					}
					if (key == "W") {
						p.width = std::stod(std::any_cast<std::string>(val));
						if (found_anti) {
							p_anti.width = p.width;
						}
					}
				}
				// Set new modified to the PDG table
				PDG.PDG_table[pdg]  = p;
				if (found_anti) {
				PDG.PDG_table[-pdg] = p_anti; 
				}

				std::cout << rang::fg::red <<
					"MProcess::SetProcess: New particle properties set with @PDG[ID]{key:val} syntax:" << rang::fg::reset << std::endl;
				p.print();

			} catch (...) {
				throw std::invalid_argument("MProcess::SetProcess: syntax problem in @PDG[ID] with ID = " + std::to_string(pdg));
			}
		}
	}

	// Remove whitespace
	std::string::iterator end_pos = std::remove(str.begin(), str.end(), ' ');
	str.erase(end_pos, str.end());

	//ID = "PP[RES]<C>"; // EXAMPLE of valid
	std::string first  = "";
	std::string second = "";
	bool found_first = false;

	for (const auto& i : indices(str)) {
	
	    if (str[i] != '[' && !found_first) {
	    	first += str[i]; continue;
	    }
	    else if (str[i] == '[' && !found_first) {
	    	found_first = true; continue;
	    }
	    if (found_first && str[i] != ']') {
	    	second += str[i];
	    }
	    if (found_first && str[i] == ']'){
	    	break;
	    }
	}
	// Setup subprocess
	ProcPtr = MSubProc(first, second, PDG);
	ConstructProcesses();
	
	// Check do we find the process
	if (Processes.count(str)) {
		// fine
	} else {
		throw std::invalid_argument(
		    "MProcess::SetProcess: Unknown PROCESS: " + str);
	}
}

} // gra namespace
