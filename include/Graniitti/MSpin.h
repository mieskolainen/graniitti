// Spin polarization and correlation functions
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.


#ifndef MSPIN_H
#define MSPIN_H

// C++
#include <complex>
#include <random>
#include <vector>

// Own
#include "Graniitti/MKinematics.h"

namespace gra {
namespace spin {

	// Amplitude squared
	std::complex<double> SpinAmp(const gra::LORENTZSCALAR& lts, gra::PARAM_RES& resonance);
	
	MMatrix<std::complex<double>> DMatrix(unsigned int J, double theta_mother, double phi_mother);
	MMatrix<std::complex<double>> fMatrix(const MMatrix<std::complex<double>>& T,
										  unsigned int J, double s1, double s2, double theta, double phi);
	std::vector<double> SpinProjections(double J);
	void InitTMatrix(gra::HELMatrix& hc, const gra::MParticle& p, const gra::MParticle& p1, const gra::MParticle& p2);
	
	// Spin-Statistics
	bool  BoseSymmetry(int l, int s);
	bool FermiSymmetry(int l, int s);
	
	// Spin algebra functions
	double ClebschGordan(double j1, double j2, double m1, double m2, double j, double m);
	double CGrules(double j1, double j2, double m1, double m2, double j, double m);
	bool chalfint(double x);
	bool cint(double x);
	bool isequal(double x, double y);
	double W3j(double j1, double j2, double j3, double m1, double m2, double m3);
	bool W3jrules(double j1, double j2, double j3, double m1, double m2, double m3);
	double TriangleCoeff(double j1, double j2, double j3);
	std::complex<double> WignerD(double theta, double phi, double m, double mp, double J);
	double WignerSmalld(double theta, double m, double mp, double J);

	// Density matrix functions
	bool Positivity(const MMatrix<std::complex<double>>& rho, unsigned int J);
	MMatrix<std::complex<double>> RandomRho(unsigned int J, bool parity, std::mt19937_64& rng);
	double VonNeumannEntropy(const MMatrix<std::complex<double>>& rho);

} // spin namespace
} // gra namespace


#endif