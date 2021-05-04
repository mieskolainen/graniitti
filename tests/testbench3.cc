// GRANIITTI Testbench 3 (unit tests)
//
// (c) 2017-2021 Mikael Mieskolainen
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


// Matrix trace and dagger
//
TEST_CASE("MMatrix:: Dagger and Trace", "[MMatrix]") {

	const double EPS = 1e-5;

	const MMatrix<std::complex<double>> A =
		 { {std::complex<double>(0.4),  std::complex<double>(0.0)},
		   {std::complex<double>(0.0), -std::complex<double>(0.4)}};
	A.Print();

	const MMatrix<std::complex<double>> A_dagger = A.Transpose();
	A_dagger.Print();

	const MMatrix<std::complex<double>> AAT = A*A_dagger;
	AAT.Print();

	const double trace = std::real( AAT.Trace() );

	REQUIRE( trace == Approx(0.32).epsilon(EPS) );
}


// Test Lorentz frame transformation functions
// 
TEST_CASE("gra::kinematics::LorentzFrame versus PGframe, HEframe, CSframe", "[Lorentz Frames]") {
	
	const bool DEBUG = false;

	const int N = 1e5; // Number of trials
	const double EPS = 1e-5;

    MRandom rng;
    rng.SetSeed(123456);

	const double mproton = 0.938;
	const double mpion = 0.14; // pion
	std::vector<M4Vec> pf(2);

	int repeat = 0;
	do {
		const int direction = 1; // GJ and PG frames

		if (repeat % 10000 == 0)
		std::cout << repeat << " / " << N << std::endl;

		// Random mother mass from the threshold up to 100 GeV
		const double M0 = rng.U(mpion * 2 + 0.01, 100);

		// Random mother 3-momentum
		const double px = rng.U(-100,100);
		const double py = rng.U(-100,100);
		const double pz = rng.U(-100,100);
		
		M4Vec mother;
		mother.SetPxPyPzM(px, py, pz, M0);

		// Random 1 -> 2 process
		gra::kinematics::MCW x;
 		x = gra::kinematics::TwoBodyPhaseSpace(mother, M0, {mpion, mpion}, pf, rng);

		// Proton beams
	    const double PZ = rng.U(50, 7000); // GeV
		const M4Vec p1(0,0, PZ, sqrt( pow2(mproton) + pow2(PZ)));
		const M4Vec p2(0,0,-PZ, sqrt( pow2(mproton) + pow2(PZ)));

		// --------------------------------------------------------------------
		// ** TRANSFORM TO DIFFERENT LORENTZ FRAMES **
		const M4Vec X = pf[0] + pf[1];

		M4Vec              pb1boost;  // beam1 particle boosted
		M4Vec              pb2boost;  // beam2 particle boosted
		std::vector<M4Vec> pfboost;   // particles boosted
		gra::kinematics::LorentFramePrepare(pf, X, p1, p2, pb1boost, pb2boost, pfboost);

		std::vector<std::string> frametype = {"CS", "HX", "AH", "PG", "CM"};
		
		for (const auto &k : indices(frametype)) {

			// Transform using the generic routine
			std::vector<M4Vec> pfout;
			gra::kinematics::LorentzFrame(pfout, pb1boost, pb2boost, pfboost, frametype[k], direction);

			// Pseudo-Gottfried-Jackson frame test
			if (frametype[k] == "PG") {
				//printf("PG[0]:      "); pfout[0].Print();
				//printf("PG[1]:      "); pfout[1].Print();				

		  		std::vector<M4Vec> pfPG = pf;
		  		gra::kinematics::PGframe(pfPG, X, direction, p1, p2, DEBUG);
				//printf("PGframe[0]: "); pfPG[0].Print();
				//printf("PGframe[1]: "); pfPG[1].Print();

				// Difference
				REQUIRE( gra::math::CheckEMC(pfout[0] - pfPG[0], EPS) );
				REQUIRE( gra::math::CheckEMC(pfout[1] - pfPG[1], EPS) );		
			}

			// Helicity frame test
			if (frametype[k] == "HX") {
				//printf("HX[0]:      "); pfout[0].Print();
				//printf("HX[1]:      "); pfout[1].Print();				

		  		std::vector<M4Vec> pfHX = pf;
		  		gra::kinematics::HXframe(pfHX, X, DEBUG);
				//printf("HXframe[0]: "); pfHX[0].Print();
				//printf("HXframe[1]: "); pfHX[1].Print();				

				// Difference
				REQUIRE( gra::math::CheckEMC(pfout[0] - pfHX[0], EPS) );
				REQUIRE( gra::math::CheckEMC(pfout[1] - pfHX[1], EPS) );	
			}

			// Collins-Soper frame test
			if (frametype[k] == "CS") {
				//printf("CS[0]:      "); pfout[0].Print();
				//printf("CS[1]:      "); pfout[1].Print();				

		  		std::vector<M4Vec> pfCS = pf;
		  		gra::kinematics::CSframe(pfCS, X, p1, p2, DEBUG);
				//printf("CSframe[0]: "); pfCS[0].Print();
				//printf("CSframe[1]: "); pfCS[1].Print();
				
				// Difference
				REQUIRE( gra::math::CheckEMC(pfout[0] - pfCS[0], EPS) );
				REQUIRE( gra::math::CheckEMC(pfout[1] - pfCS[1], EPS) );	
			}
		}

	} while ((++repeat) < N);

}


// Test Lorentz frame transformation functions
// 
TEST_CASE("gra::kinematics::GJframe", "[Lorentz Frames]") {
	
	const bool DEBUG = false;

	const int N = 1e5; // Number of trials
	const double EPS = 1e-5;

    MRandom rng;
    rng.SetSeed(123456);

	const double mproton = 0.938;
	const double mpion = 0.14; // pion
	std::vector<M4Vec> pf4(4);	

	int repeat = 0;
	do {
		const int direction = 1; // GJ and PG frames
		
		if (repeat % 10000 == 0)
		std::cout << repeat << " / " << N << std::endl;

		// Proton beams
	    const double PZ = rng.U(50, 7000); // GeV
		const M4Vec p1(0,0, PZ, sqrt( pow2(mproton) + pow2(PZ)));
		const M4Vec p2(0,0,-PZ, sqrt( pow2(mproton) + pow2(PZ)));

		// Random 2 -> 4 process
		gra::kinematics::NBodyPhaseSpace(p1 + p2, (p1+p2).M(), {mproton, mpion, mpion, mproton}, pf4, false, rng);

 		const M4Vec X  = pf4[1] + pf4[2];
  		const M4Vec q1 = p1 - pf4[0];
  		const M4Vec q2 = p2 - pf4[3];

  		// Transform q1 and q2
  		std::vector<M4Vec> out = {q1, q2};
  		gra::kinematics::GJframe(out, X, direction, q1, q2, DEBUG);

  		// Propagators should be aligned back-to-back on z-axis
		REQUIRE( std::abs(out[0].Px()) < EPS ); // ~ 0
		REQUIRE( std::abs(out[0].Py()) < EPS ); // ~ 0
		REQUIRE( std::abs(out[1].Px()) < EPS ); // ~ 0
		REQUIRE( std::abs(out[1].Py()) < EPS ); // ~ 0
		REQUIRE( std::abs(out[0].Pz() + out[1].Pz()) < EPS ); // pz_1 + pz_2 ~ 0

	} while ((++repeat) < N);
}

