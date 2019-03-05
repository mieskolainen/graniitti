// Gamma amplitudes
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MGAMMA_H
#define MGAMMA_H

// C++
#include <complex>
#include <random>
#include <vector>

// MadGraph
#include "Graniitti/Amplitude/MAmpMG5_yy_ll.h"
#include "Graniitti/Amplitude/MAmpMG5_yy_ll_2to4.h"
#include "Graniitti/Amplitude/MAmpMG5_yy_ww.h"

// Own
#include "Graniitti/MKinematics.h"
#include "Graniitti/M4Vec.h"


namespace gra {


// "Functionoid class"
// Matrix element dimension: " GeV^" << -(2*external_legs - 8)
class MGamma {

public:

	MGamma(){}
	~MGamma(){}

	// yy->resonance X
	std::complex<double> yyX(const gra::LORENTZSCALAR& lts, gra::PARAM_RES& resonance) const;
	
	// yy->monopolium
	std::complex<double> yyMP(const gra::LORENTZSCALAR& lts) const;

	// yy->lepton pair, or monopole antimonopole amplitude
	std::complex<double> yyffbar(gra::LORENTZSCALAR& lts);
	
protected:

	// MADGRAPH amplitudes added here
	MAmpMG5_yy_ll_2to4  AmpMG5_yy_ll_2to4;
	MAmpMG5_yy_ll       AmpMG5_yy_ll;
	MAmpMG5_yy_ww       AmpMG5_yy_ww;

};

} // gra namespace ends

#endif
