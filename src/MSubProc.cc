// (Sub)-Processes and Amplitudes
//
// Next step: update to use functor binding.
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <complex>
#include <random>
#include <vector>

// Own
#include "Graniitti/MPDG.h"
#include "Graniitti/MSubProc.h"
#include "Graniitti/MAux.h"
#include "Graniitti/MKinematics.h"
#include "Graniitti/MSpin.h"
#include "Graniitti/MForm.h"
#include "Graniitti/MGlobals.h"

// LHAPDF
#include "LHAPDF/LHAPDF.h"

// Libraries
#include "rang.hpp"


using gra::math::msqrt;
using gra::math::pow2;
using gra::math::pow4;
using gra::math::zi;


namespace gra {

// Initialize
MSubProc::MSubProc(const std::string& _ISTATE, const std::string& _CHANNEL, const MPDG& _PDG) {

	ISTATE  = _ISTATE;
	CHANNEL = _CHANNEL;	
	PDG     = _PDG;

	ConstructDescriptions(ISTATE);

	// Construct 4 and 6 body Regge Ladder leg permutations
	const int mode = 1; // charged permutations
	permutations4_ = gra::math::GetAmpPerm(4, mode);
	permutations6_ = gra::math::GetAmpPerm(6, mode);

	// ** Read in monopole mass (needed by monopolium process) **
	PARAM_MONOPOLE::M0 = _PDG.FindByPDG(PDG::PDG_monopole).mass;
}

MSubProc::MSubProc(const std::vector<std::string>& first) {
	for (std::size_t i = 0; i < first.size(); ++i) {
		ConstructDescriptions(first[i]);
	}
}

void MSubProc::ConstructDescriptions(const std::string& first) {

	if      (first == "X") {

		std::map<std::string, std::string> channels;
		channels.insert(std::pair<std::string, std::string>("EL",           "Elastic                                  [Eikonal Pomeron]"));
		channels.insert(std::pair<std::string, std::string>("SD",           "Single Diffractive                       [Triple Pomeron]      @@ UNDER CONSTRUCTION @@"));
		channels.insert(std::pair<std::string, std::string>("DD",           "Double Diffractive                       [Triple Pomeron]      @@ UNDER CONSTRUCTION @@"));
		channels.insert(std::pair<std::string, std::string>("ND",           "Non-Diffractive                          [N-cut soft Pomerons] @@ UNDER CONSTRUCTION @@"));
			descriptions.insert(std::pair<std::string, std::map<std::string,std::string>>("X", channels));
	}
	else if (first == "PP") {

		std::map<std::string, std::string> channels;
		channels.insert(std::pair<std::string, std::string>("CONTENSOR",    "Regge continuum 2-body                   [Tensor Pomeron] x [Tensor Pomeron] @@ UNDER TESTING @@"));
		channels.insert(std::pair<std::string, std::string>("CON",          "Regge continuum 2/4/6-body               [Pomeron] x [Pomeron]"));
		channels.insert(std::pair<std::string, std::string>("CON-",         "Regge continuum 2-body negative subamp   [Pomeron] x [Pomeron]"));
		channels.insert(std::pair<std::string, std::string>("RES+CON",      "Regge cont. 2-body + complex resonances  [Pomeron] x [Gamma/Pomeron]"));
		channels.insert(std::pair<std::string, std::string>("RES",          "Regge parametric resonance               [Pomeron] x [Pomeron]"));
		channels.insert(std::pair<std::string, std::string>("RESHEL",       "Regge sliding helicity amplitudes        [Pomeron] x [Pomeron] @@ UNDER CONSTRUCTION @@"));
			descriptions.insert(std::pair<std::string, std::map<std::string,std::string>>("PP", channels));
	}
	else if (first == "yP") {
		
		std::map<std::string, std::string> channels;
		channels.insert(std::pair<std::string, std::string>("RES",          "Photoproduced parametric resonance       [EPA] x [Pomeron]"));
			descriptions.insert(std::pair<std::string, std::map<std::string,std::string>>("yP", channels));
	}
	else if (first == "yy") {

		std::map<std::string, std::string> channels;
		channels.insert(std::pair<std::string, std::string>("RES",          "Gamma-Gamma to parametric resonance      [EPA] x [EPA]"));
		channels.insert(std::pair<std::string, std::string>("Higgs",        "Gamma-Gamma to SM Higgs                  [EPA] x [EPA]"));
		channels.insert(std::pair<std::string, std::string>("monopolium(0)","Gamma-Gamma to Monopolium (J=0)          [EPA] x [EPA]         @@ UNDER TESTING @@"));
		channels.insert(std::pair<std::string, std::string>("CON",          "Gamma-Gamma to l+l-, w+w-, monopoles     [EPA] x [EPA]"));
		channels.insert(std::pair<std::string, std::string>("QED",          "Gamma-Gamma qq -> q l+ l- q              [FULL QED]            @@ TEST PROCESS @@"));
			descriptions.insert(std::pair<std::string, std::map<std::string,std::string>>("yy", channels));
	}
	else if (first == "gg") {

		std::map<std::string, std::string> channels;
		channels.insert(std::pair<std::string, std::string>("chic(0)",      "QCD resonance chic(0)                    [Durham QCD]"));
		channels.insert(std::pair<std::string, std::string>("CON",          "QCD continuum to gg, 2 x pseudoscalar    [Durham QCD]          @@ UNDER TESTING @@"));
			descriptions.insert(std::pair<std::string, std::map<std::string,std::string>>("gg", channels));
	}	
	else if (first == "yy_DZ") {

		std::map<std::string, std::string> channels;
		channels.insert(std::pair<std::string, std::string>("CON",          "Collinear yy to l+l, w+w- or monopoles   [DZ] x [DZ]"));
			descriptions.insert(std::pair<std::string, std::map<std::string,std::string>>("yy_DZ", channels));	
	}
	else if (first == "yy_LUX") {

		std::map<std::string, std::string> channels;
		channels.insert(std::pair<std::string, std::string>("CON",          "Collinear yy to l+l, w+w- or monopoles   [LUX] x [LUX]         @@ UNDER TESTING @@"));
			descriptions.insert(std::pair<std::string, std::map<std::string,std::string>>("yy_LUX", channels));	
	}
}


// This is called last by the initialization routines
// as the last step before event generation.
void MSubProc::SetTechnicalBoundaries(gra::GENCUT& gcuts, unsigned int EXCITATION) {
	
	if (gcuts.forward_pt_min < 0.0) { // Not set yet
		gcuts.forward_pt_min = 0.0;
	}
	
	if (gcuts.forward_pt_max < 0.0) { // Not set yet
		
		if      (EXCITATION == 0) {   // Fully elastic
			if      (ISTATE == "PP") {
			gcuts.forward_pt_max = 1.75;
			}
			else if (ISTATE == "yP") {
			gcuts.forward_pt_max = 1.75;
			}
			else if (ISTATE == "yy") {
			gcuts.forward_pt_max = 1.5;
			}
			else if (ISTATE == "gg") {
			gcuts.forward_pt_max = 2.0;
			}
		}

		// Here, we cannot put too high if we use on-shell matrix (e.g. in EPA)
		// because t ~ -pt^2
		else if (EXCITATION == 1) {   // Single
			gcuts.forward_pt_max = 3.0;
		}
		else if (EXCITATION == 2) {   // Double
			gcuts.forward_pt_max = 3.0;
		}
	}
}


std::complex<double> MSubProc::GetBareAmplitude(gra::LORENTZSCALAR& lts) {

	if      (ISTATE == "X") {
		return GetBareAmplitude_X(lts);
	}
	else if (ISTATE == "PP") {
		return GetBareAmplitude_PP(lts);
	}
	else if (ISTATE == "yP") {
		return GetBareAmplitude_yP(lts);
	}
	else if (ISTATE == "yy") {
		return GetBareAmplitude_yy(lts);
	}
	else if (ISTATE == "gg") {
		return GetBareAmplitude_gg(lts);
	}
	else if (ISTATE == "yy_DZ") {
		return GetBareAmplitude_yy_DZ(lts);
	}
	else if (ISTATE == "yy_LUX") {
		return GetBareAmplitude_yy_LUX(lts);
	} else {
		throw std::invalid_argument("MSubProc::GetBareAmplitude: Unknown ISTATE '" + ISTATE + '"');
	}
}


// Inclusive processes
inline std::complex<double> MSubProc::GetBareAmplitude_X(gra::LORENTZSCALAR& lts) {

	std::complex<double> A(0, 0);
	if (CHANNEL == "EL" || CHANNEL == "SD" || CHANNEL == "DD") {
		A = ME2(lts, LIPSDIM);
	} else if (CHANNEL == "ND") {
		A = 1.0;
	} else {
		throw std::invalid_argument("MSubProc::GetBareAmplitude: Unknown CHANNEL = " + CHANNEL);
	}
	return A;
}


// Pomeron-Pomeron
inline std::complex<double> MSubProc::GetBareAmplitude_PP(gra::LORENTZSCALAR& lts) {
	
	// ------------------------------------------------------------------
	// First run, init parameters
	gra::g_mutex.lock();
	if (!PARAM_REGGE::initialized) { 
		const int PDG = std::abs(lts.decaytree[0].p.pdg);
		InitReggeAmplitude(PDG, gra::MODELPARAM);
		PARAM_REGGE::initialized = true;
	}
	gra::g_mutex.unlock();
	// ------------------------------------------------------------------
	
	std::complex<double> A(0, 0);


	if      (CHANNEL == "RES") {
		A = ME3(lts, lts.RESONANCES.begin()->second);
	}
	else if (CHANNEL == "RESHEL") {
		A = ME3HEL(lts, lts.RESONANCES.begin()->second);
	}

	else if (CHANNEL == "CONTENSOR") {
		if      (lts.decaytree.size() == 2) {
			static MTensorPomeron TensorPomeron;
			A = TensorPomeron.ME4(lts);
		} else {
			throw std::invalid_argument("MSubProc: Only 2-body + sequential final states for [CONTENSOR] process");
		}
	}
	
	else if (CHANNEL == "CON") {
		if      (lts.decaytree.size() == 2) {
			A = ME4(lts,1);
		}
		else if (lts.decaytree.size() == 4) {
			A = ME6(lts);
		}
		else if (lts.decaytree.size() == 6) {
			A = ME8(lts);
		} else {
			throw std::invalid_argument("MSubProc: Only 2/4/6-body + sequential final states for [CON] process");
		}
	}
	
	else if (CHANNEL == "CON-") {
		if      (lts.decaytree.size() == 2) {
			A = ME4(lts, -1);
		} else {
			throw std::invalid_argument("MSubProc: Only 2-body final states for [CON-] process");
		}
	}
	
	else if (CHANNEL == "RES+CON") {
		A = ME4RES(lts, lts.RESONANCES, 1);
	} else {
		throw std::invalid_argument("MSubProc::GetBareAmplitude: Unknown CHANNEL = " + CHANNEL);
	}
	return A;
}

// Gamma-Pomeron
inline std::complex<double> MSubProc::GetBareAmplitude_yP(gra::LORENTZSCALAR& lts) {
	std::complex<double> A(0, 0);

	if (CHANNEL == "RES") {
		A = PhotoME3(lts, lts.RESONANCES.begin()->second);
	} else {
		throw std::invalid_argument("MSubProc::GetBareAmplitude: Unknown CHANNEL = " + CHANNEL);
	}
	return A;
}

// Gamma-Gamma
inline std::complex<double> MSubProc::GetBareAmplitude_yy(gra::LORENTZSCALAR& lts) {
	std::complex<double> A(0, 0);
	
	if      (CHANNEL == "RES") {
		A = yyX(lts, lts.RESONANCES.begin()->second);
	}
	else if (CHANNEL == "Higgs") {
		A = yyHiggs(lts);
	}
	else if (CHANNEL == "monopolium(0)") {
		A = yyMP(lts);
	}
	else if (CHANNEL == "CON") {
		if (std::abs(lts.decaytree[0].p.pdg) == 24 &&
			std::abs(lts.decaytree[1].p.pdg) == 24) {   // W+W-
			A = AmpMG5_yy_ww.CalcAmp(lts);
		} else { // ffbar
			A = yyffbar(lts);
		}
	}
	else if (CHANNEL == "QED") {
		A = AmpMG5_yy_ll_2to4.CalcAmp(lts);
	} else {
		throw std::invalid_argument("MSubProc::GetBareAmplitude: Unknown CHANNEL = " + CHANNEL);
	}
	
 	// Apply non-collinear EPA fluxes
	const double gammaflux1 = lts.excite1 ? gra::form::IncohFlux(lts.x1, lts.t1, lts.qt1, lts.pfinal[1].M2()) : form::CohFlux(lts.x1, lts.t1, lts.qt1);
	const double gammaflux2 = lts.excite2 ? gra::form::IncohFlux(lts.x2, lts.t2, lts.qt2, lts.pfinal[2].M2()) : form::CohFlux(lts.x2, lts.t2, lts.qt2);
	const double phasespace = lts.s / lts.s_hat; // Moller flux replacement
	
	// To "amplitude level"
	const double flux = msqrt( gammaflux1 * gammaflux2 * phasespace);

	// Helicity amplitudes
	for (const auto& i : aux::indices(lts.hamp)) {
		lts.hamp[i] *= flux;
	}

	return A * flux;
}


// Gamma-Gamma collinear Drees-Zeppenfeld (coherent flux)
inline std::complex<double> MSubProc::GetBareAmplitude_yy_DZ(gra::LORENTZSCALAR& lts) { 

	//throw std::invalid_argument("yy_DZ::Not active in this version");

	// Amplitude
	std::complex<double> A = yyffbar(lts);

	// Evaluate gamma pdfs
	const double   f1 = form::DZFlux(lts.x1);
	const double   f2 = form::DZFlux(lts.x2);
	const double amp2 = f1 * f2 * math::abs2(A) * (lts.s / lts.s_hat);

	return msqrt(amp2); // We take square later
}


// Gamma-Gamma LUX-pdf (use at \mu > 10 GeV)
// *** UNDER TESTING ***
inline std::complex<double> MSubProc::GetBareAmplitude_yy_LUX(gra::LORENTZSCALAR& lts) {

	// @@ MULTITHREADING LOCK @@
	gra::g_mutex.lock();
	if (gra::GlobalPdfPtr == nullptr) {
	std::string pdfname = gra::LHAPDF;

retry:
	try {
		gra::GlobalPdfPtr = LHAPDF::mkPDF(pdfname, 0);
		gra::pdf_trials = 0; // fine
	} catch (...) {
		++gra::pdf_trials;
		std::string str = "MSubProc::InitLHAPDF: Problem with reading '" + pdfname + "'";
		aux::AutoDownloadLHAPDF(pdfname); // Try autodownload
		gra::g_mutex.unlock();       // Remember before throw, otherwise deadlock
		if (gra::pdf_trials >= 2) { // too many failures
			throw std::invalid_argument(str);
		} else {
			goto retry;
		}
	}
	}
	gra::g_mutex.unlock();
	// @@ MULTITHREADING UNLOCK @@

	// Amplitude
	std::complex<double> A(0, 0);
	if (CHANNEL == "CON") {
		if (std::abs(lts.decaytree[0].p.pdg) == 24 &&
			std::abs(lts.decaytree[1].p.pdg) == 24) {   // W+W-
			A = AmpMG5_yy_ww.CalcAmp(lts);
		} else { // ffbar
			A = yyffbar(lts);
		}
	} else {
		throw std::invalid_argument("MSubProc::GetBareAmplitude: Unknown CHANNEL = " + CHANNEL);
	}

	// @@ MULTITHREADING LOCK @@
	gra::g_mutex.lock();

	// pdf factorization scale
	const double Q2 = lts.s_hat / 4.0; 
	
	// Evaluate gamma pdfs
	double f1 = 0.0;
	double f2 = 0.0;
	try {
	// Divide x out
	f1 = gra::GlobalPdfPtr->xfxQ2(PDG::PDG_gamma, lts.x1, Q2) / lts.x1;
	f2 = gra::GlobalPdfPtr->xfxQ2(PDG::PDG_gamma, lts.x2, Q2) / lts.x2;
	} catch (...) {
		gra::g_mutex.unlock(); // remember before throw, otherwise deadlock		
		throw std::invalid_argument("MSubProc:yy_LUX: Failed evaluating LHAPDF");
	}
	gra::g_mutex.unlock();

	const double tot = f1 * f2 * math::abs2(A) * (lts.s / lts.s_hat);
	return msqrt(tot);
}


// Durham gg
inline std::complex<double> MSubProc::GetBareAmplitude_gg(gra::LORENTZSCALAR& lts) {
	std::complex<double> A(0,0);

	if (CHANNEL == "chic(0)") {
		A = DurhamQCD(lts, CHANNEL);
	}
	else if (CHANNEL == "CON") {
		if (std::abs(lts.decaytree[0].p.pdg) == 21 && 
			std::abs(lts.decaytree[1].p.pdg) == 21) {
			A = DurhamQCD(lts, "gg");
		} else {
			A = DurhamQCD(lts, "MMbar");
		}
	} else {
		throw std::invalid_argument("MSubProc::GetBareAmplitude: Unknown CHANNEL = " + CHANNEL);
	}
	return A;
}

} // gra namespace ends
