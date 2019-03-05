// Functional methods for Spherical Harmonic Expansions
// 
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MSPHERICAL_H
#define MSPHERICAL_H

// C++
#include <iostream>
#include <vector>
#include <complex>

// Own
#include "Graniitti/MMatrix.h"
#include "Graniitti/MMath.h"

namespace gra {
namespace spherical {

// Microevent structure
struct Omega {
	
	// Decay daughter
	// rest frame variables
	double costheta = 0.0;
	double phi      = 0.0;
	
	// System lab frame
	double M        = 0.0;
	double Pt       = 0.0;
	double Y        = 0.0;
	
	bool fiducial   = false;
	bool selected   = false;
};


// Data for one hypercell (M,Pt,Y)
struct SH {

	// Moment mixing matrix
	MMatrix<double>     MIXlm;

	// Efficiency decomposition coefficients
	std::vector<double> E_lm;
	std::vector<double> E_lm_error;

	// Diractly (algebraic) observed moments
	std::vector<double> t_lm_OBS;
	std::vector<double> t_lm_OBS_error;		

	// Maximum Likelihood fitted moments
	std::vector<double> t_lm_FIT;
	std::vector<double> t_lm_FIT_error;
};


MMatrix<double>     GetGMixing(const std::vector<Omega>& events,
							   const std::vector<std::size_t>& ind, int LMAX, int mode);

std::pair<std::vector<double>,std::vector<double>>
                    GetELM(const std::vector<Omega>& MC,
                    	   const std::vector<std::size_t>& ind, int LMAX, int mode);

std::vector<double> SphericalMoments(const std::vector<Omega>& input,
	                                 const std::vector<std::size_t>& ind, int LMAX, int mode);

MMatrix<double>     YLM(const std::vector<Omega>& events, int LMAX);

double HarmDotProd(const std::vector<double>& G, const std::vector<double>& x,
	               const std::vector<bool>& ACTIVE, int LMAX);

void   PrintOutMoments(const std::vector<double>& x, const std::vector<double>& x_error,
	                   const std::vector<bool>& ACTIVE, int LMAX);

std::vector<std::size_t>   GetIndices(const std::vector<Omega>& events,
				  					  const std::vector<double>& M,
				  					  const std::vector<double>& Pt,
				  					  const std::vector<double>& Y);

void   TestSphericalIntegrals(int LMAX);
int    LinearInd(int l, int m);
double CalcError(double f2, double f, double N);
void   PrintMatrix(FILE* fp, const std::vector<std::vector<double>>& A);


} // spherical namespace
} // gra namespace


#endif