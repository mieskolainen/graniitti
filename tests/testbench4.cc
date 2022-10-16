// GRANIITTI Testbench 4 (unit tests)
//
// (c) 2017-2022 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include <catch.hpp>
#include <random>

#include "Graniitti/MRandom.h"
#include "Graniitti/MMath.h"
#include "Graniitti/MMatrix.h"
#include "Graniitti/MKinematics.h"
#include "Graniitti/M4Vec.h"
#include "Graniitti/MDirac.h"
#include "Graniitti/MTensorPomeron.h"
#include "Graniitti/MForm.h"
#include "Graniitti/MRegge.h"
#include "Graniitti/MDurham.h"
#include "Graniitti/MGamma.h"
#include "Graniitti/MQED.h"


using namespace gra;

using gra::aux::indices;
using gra::math::pow2;

const std::string modelfile = gra::aux::GetBasePath(2) + "/modeldata/TUNE0/GENERAL.json";


// Test initializing Tensor Pomeron amplitudes
// 
TEST_CASE("MTensorPomeron", "[gra::MTensorPomeron]") {

	gra::LORENTZSCALAR lts;
	MTensorPomeron a(lts, modelfile);
}

// Test initializing Regge amplitudes
//
TEST_CASE("MRegge", "[gra::MRegge]") {

	gra::LORENTZSCALAR lts;
	MRegge a(lts, modelfile);
}

// Test initializing Durham model amplitudes
//
TEST_CASE("MDurham", "[gra::MDurham]") {

	gra::LORENTZSCALAR lts;
	lts.sqrt_s = 13000;
	lts.LHAPDFSET = "CT10nlo";
	MDurham a(lts, modelfile);	
}

// Test initializing gamma-gamma amplitudes
//
/*
TEST_CASE("MGamma", "[gra::MGamma]") {

	gra::LORENTZSCALAR lts;
	MGamma a(lts, modelfile);
}
*/

// Test running QED coupling
//
TEST_CASE("qed::alpha_QED", "[gra::qed]") {

	std::vector<double> Qval = {0.0, 1e-12, 1e-9, 1e-6, 0.0005109989461, 1.0, 10.0, 91.2, 100, 1000, 10000};

	double Q2 = 0.0;
	for (std::size_t i = 0; i < Qval.size(); ++i) {
		Q2 = Qval[i]*Qval[i];
		printf("Q = %9.6E GeV \t 1/alpha_QED(LL) = %0.6f \n", Qval[i], 1.0/qed::alpha_QED(Q2, "LL"));
	}

	// Zero momentum scale
	Q2 = 0;
	REQUIRE( 1.0/qed::alpha_QED(Q2, "LL") == Approx(137.0359).epsilon(0.01) );
	
	// Z-pole scale
	Q2 = math::pow2(91.2);
	REQUIRE( 1.0/qed::alpha_QED(Q2, "LL") == Approx(128.5).epsilon(0.01) );	
}

// Test form factors
//
TEST_CASE("form:: Form factors", "[gra::form]") {

	const double EPS_DIPOLE = 1e-5;
	const double EPS_KELLY  = 0.4;

	gra::LORENTZSCALAR lts;
	MTensorPomeron tensor(lts, modelfile);

	std::cout << "F1(0): " << form::F1(0)  << std::endl;
	std::cout << "F2(0): " << form::F2(0)  << std::endl;
	std::cout << "F1_(0): " << tensor.F1_(0) << std::endl;
	std::cout << "F2_(0): " << tensor.F2_(0) << std::endl;

	gra::PARAM_STRUCTURE::EM = "DIPOLE";
	REQUIRE( form::F1(0) == Approx( tensor.F1_(0) ).epsilon(1e-5) );
	REQUIRE( form::F2(0) == Approx( tensor.F2_(0) ).epsilon(1e-5) );
	
	gra::PARAM_STRUCTURE::EM = "KELLY";
	REQUIRE( form::F1(0) == Approx( tensor.F1_(0) ).epsilon(1e-5) );
	REQUIRE( form::F2(0) == Approx( tensor.F2_(0) ).epsilon(1e-5) );	

	MRandom rng;
	
	for (std::size_t i = 0; i < 1e2; ++i) {
		const double t = rng.U(-5,0);

		gra::PARAM_STRUCTURE::EM = "DIPOLE";
		REQUIRE( form::F1(t) == Approx( tensor.F1_(t) ).epsilon(EPS_DIPOLE) );
		REQUIRE( form::F2(t) == Approx( tensor.F2_(t) ).epsilon(EPS_DIPOLE) );

		gra::PARAM_STRUCTURE::EM = "KELLY";
		REQUIRE( form::F1(t) == Approx( tensor.F1_(t) ).epsilon(EPS_KELLY) );
		REQUIRE( form::F2(t) == Approx( tensor.F2_(t) ).epsilon(EPS_KELLY) );
	}
}

// Regge trajectory phase
//
TEST_CASE("gra::form:: ReggeEta and ReggeEtaLinear", "[Form]") {

	const double EPS = 0.05;
	MRandom rng;
	
	for (std::size_t i = 0; i < 1e2; ++i) {
		
		// Mandelstam t
		const double t = rng.U(-1.0, 0);

		// Pomeron trajectory
		const double alpha_t0 = 1.08;
		const double ap       = 0.25;
		const double sigma    = 1; // even signature
		const double alpha    = alpha_t0 + ap * t;

		const std::complex<double> eta1 = gra::form::ReggeEtaLinear(t, alpha_t0, ap, sigma);
		const std::complex<double> eta2 = gra::form::ReggeEta(alpha, sigma);

		if (i < 2) {
		std::cout << "Mandelstam t:   " << t    << std::endl;
		std::cout << "ReggeEtaLinear: " << eta1 << std::endl;
		std::cout << "ReggeEta:       " << eta2 << std::endl;
		}

		REQUIRE( std::real(eta1) == Approx(std::real(eta2)).epsilon(EPS) );
		REQUIRE( std::imag(eta1) == Approx(std::imag(eta2)).epsilon(EPS) );
	}
}
