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
#include "Graniitti/MSubProc.h"
#include "Graniitti/MAux.h"
#include "Graniitti/MKinematics.h"
#include "Graniitti/MSpin.h"
#include "Graniitti/MForm.h"

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

	// 4 and 6 body Regge Ladder leg permutations
	const int mode = 1;
	permutations4_ = gra::math::GetAmpPerm(4, mode);
	permutations6_ = gra::math::GetAmpPerm(6, mode);

	// ** Read in monopole mass (needed by monopolium process) **
	PARAM_MONOPOLE::M0 = _PDG.FindByPDG(992).mass;
}

MSubProc::MSubProc(const std::vector<std::string>& first) {
	for (std::size_t i = 0; i < first.size(); ++i) {
		ConstructDescriptions(first[i]);
	}
}

void MSubProc::ConstructDescriptions(const std::string& first) {

	if      (first == "X") {

		std::map<std::string, std::string> channels;
		channels.insert(std::pair<std::string, std::string>("EL",  "Elastic [Eikonal Pomeron]"));
		channels.insert(std::pair<std::string, std::string>("SD",  "Single Diffractive [Triple Pomeron] (TEST / UNDER CONSTRUCTION)"));
		channels.insert(std::pair<std::string, std::string>("DD",  "Double Diffractive [Triple Pomeron] (TEST / UNDER CONSTRUCTION)"));
		channels.insert(std::pair<std::string, std::string>("ND",  "Non-Diffractive [Multiple Soft Pomeron cuts] (TEST / UNDER CONSTRUCTION)"));
			descriptions.insert(std::pair<std::string, std::map<std::string,std::string>>("X", channels));
	}
	else if (first == "PP") {

		std::map<std::string, std::string> channels;
		channels.insert(std::pair<std::string, std::string>("CONTENSOR", "Continuum 2-body [Tensor Pomeron] x [Tensor Pomeron]"));
		channels.insert(std::pair<std::string, std::string>("CON",       "Continuum 2/4/6-body [Pomeron] x [Pomeron] (WARNING: 4/6-body are singular with certain form factors)"));
		channels.insert(std::pair<std::string, std::string>("CON-",      "Continuum 2-body with negative subamplitude (2-body) [Pomeron] x [Pomeron]"));
		channels.insert(std::pair<std::string, std::string>("RES+CON",   "Resonances + Continuum 2-body (interfering) [Pomeron] x [Gamma/Pomeron]"));
		channels.insert(std::pair<std::string, std::string>("RES",       "Resonance [Pomeron] x [Pomeron]"));
		channels.insert(std::pair<std::string, std::string>("RESHEL",    "Resonance with 'sliding helicity amplitudes' (TEST / UNDER CONSTRUCTION)"));
			descriptions.insert(std::pair<std::string, std::map<std::string,std::string>>("PP", channels));
	}
	else if (first == "yP") {

		std::map<std::string, std::string> channels;
		channels.insert(std::pair<std::string, std::string>("RES",       "Photoproduction of resonance [EPA] x [Pomeron]"));
			descriptions.insert(std::pair<std::string, std::map<std::string,std::string>>("yP", channels));
	}
	else if (first == "yy") {

		std::map<std::string, std::string> channels;
		channels.insert(std::pair<std::string, std::string>("RES",           "Gamma-Gamma to resonance [EPA] x [EPA]"));
		channels.insert(std::pair<std::string, std::string>("monopolium(0)", "Gamma-Gamma to Monopolium J=0 bound state [EPA] x [EPA]"));
		channels.insert(std::pair<std::string, std::string>("CON",           "Gamma-Gamma to l+l-, w+w-, monopole pair continuum [EPA] x [EPA]"));
		channels.insert(std::pair<std::string, std::string>("QED",           "Gamma-Gamma qq -> q l+ l- q [FULL QED]"));
			descriptions.insert(std::pair<std::string, std::map<std::string,std::string>>("yy", channels));
	}
	else if (first == "gg") {

		std::map<std::string, std::string> channels;
		channels.insert(std::pair<std::string, std::string>("chi_c(0)", "QCD resonance gg -> chi_c(0) [Durham QCD]"));
		channels.insert(std::pair<std::string, std::string>("CON",      "QCD continuum [Durham QCD] (-> g g currently)"));
			descriptions.insert(std::pair<std::string, std::map<std::string,std::string>>("gg", channels));
	}
}

// This is called last by the initialization routines
// as the last step before event generation.
void MSubProc::SetTechnicalBoundaries(gra::GENCUT& gcuts, unsigned int EXCITATION) {
	
	if (gcuts.forward_pt_min < 0.0) { // Not set yet
		gcuts.forward_pt_min = 0.0;
	}

	if (gcuts.forward_pt_max < 0.0) { // Not set yet

		if      (EXCITATION == 0) {  // Fully elastic
			gcuts.forward_pt_max = 2.0;
		}
		else if (EXCITATION == 1) {  // Single
			gcuts.forward_pt_max = 4.0;
		}
		else if (EXCITATION == 2) {  // Double
			gcuts.forward_pt_max = 6.0;
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
	gra::aux::g_mutex.lock();
	if (!PARAM_REGGE::initialized) { 
		const int PDG = std::abs(lts.decaytree[0].p.pdg);
		InitReggeAmplitude(PDG, gra::aux::MODELPARAM);
		PARAM_REGGE::initialized = true;
	}
	gra::aux::g_mutex.unlock();
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
			throw std::invalid_argument("MSubProc: Only 2-body+sequential final states for [CONTENSOR] process");
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
			throw std::invalid_argument("MSubProc: Only 2/4/6-body+sequential final states for [CON] process");
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

	if (CHANNEL == "RES") {
		A = yyX(lts, lts.RESONANCES.begin()->second);
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
	return A;
}

// Durham gg
inline std::complex<double> MSubProc::GetBareAmplitude_gg(gra::LORENTZSCALAR& lts) {
	std::complex<double> A(0, 0);

	if (CHANNEL == "chi_c(0)") {
		A = DurhamQCD(lts, "chi0");
	}
	else if (CHANNEL == "CON") {
		if (std::abs(lts.decaytree[0].p.pdg) == 21 && 
			std::abs(lts.decaytree[1].p.pdg) == 21) {
			A = DurhamQCD(lts, "gg");
		} else {
			A = DurhamQCD(lts, "qqbar");			
		}
	} else {
		throw std::invalid_argument("MSubProc::GetBareAmplitude: Unknown CHANNEL = " + CHANNEL);
	}
	return A;
}

} // gra namespace ends
