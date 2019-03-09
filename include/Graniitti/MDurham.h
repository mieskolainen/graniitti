// "Durham QCD" Processes and Amplitudes
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MDURHAM_H
#define MDURHAM_H

// C++
#include <complex>
#include <random>
#include <vector>

// MADGRAPH
#include "Graniitti/Amplitude/MAmpMG5_gg_qqbar.h"
#include "Graniitti/Amplitude/MAmpMG5_gg_ggg.h"
#include "Graniitti/Amplitude/MAmpMG5_gg_gg.h"


// Own
#include "Graniitti/MAux.h"
#include "Graniitti/MKinematics.h"
#include "Graniitti/MSudakov.h"
#include "Graniitti/M4Vec.h"


namespace gra {


// Durham loop integral discretization technicals
namespace Durham {

	extern unsigned int N_qt;  // (> 30)
	extern unsigned int N_phi; // (> 10)

	extern double qt2_MIN; // Loop momentum qt^2 minimum (GeV^2)
	extern double qt2_MAX; // Loop momentum qt^2 maximum (GeV^2)

	extern std::string PDF_scale; // PDF scale scheme
	extern double alphas_scale;   // alphas scale scheme

	// THESE MUST BE SET LAST
	extern double qt_MIN;
	extern double qt_MAX;
	extern double qt_STEP;
	extern double phi_STEP;

	extern bool initialized;

	void ReadParameters();
}


// Matrix element dimension: " GeV^" << -(2*external_legs - 8)
class MDurham {

public:

	MDurham() {}
	~MDurham() {}

	std::complex<double> DurhamQCD(gra::LORENTZSCALAR& lts, const std::string& process);
	std::complex<double> DQtloop(gra::LORENTZSCALAR& lts, std::vector<std::vector<std::complex<double>>> Amp);
	inline void DScaleChoise(double qt2, double q1_2, double q2_2, double& Q1_2_scale, double& Q2_2_scale) const;
	
	inline void Dgg2chi0(const gra::LORENTZSCALAR& lts,
	              std::vector<std::vector<std::complex<double>>>& Amp,
	              const std::vector<double>& qt1,
	              const std::vector<double>& qt2) const;

	inline void DHelicity(const std::vector<double>& q1,
	               const std::vector<double>& q2,
	               std::vector<std::complex<double>>& JzP) const;

	inline std::complex<double> DHelProj(
	    const std::vector<std::complex<double>>& A,
	    const std::vector<std::complex<double>>& JzP) const;

	void Dgg2gg(const gra::LORENTZSCALAR& lts,
	            std::vector<std::vector<std::complex<double>>>& Amp);
	
	void Dgg2qqbar(const gra::LORENTZSCALAR& lts,
	               std::vector<std::vector<std::complex<double>>>& Amp);

	double Asum = 0.0;
	double Nsum = 0.0;

	// ------------------------------------------------------------------
	// MadGraph amplitudes here
	MAmpMG5_gg_gg       AmpMG5_gg_gg;
	MAmpMG5_gg_ggg      AmpMG5_gg_ggg;
	MAmpMG5_gg_qqbar    AmpMG5_gg_qqbar;

private:

};


} // gra namespace ends

#endif