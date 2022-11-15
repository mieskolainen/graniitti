// GRANIITTI Testbench 3 (unit tests)
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

using namespace gra;

using gra::aux::indices;
using gra::math::pow2;


// Matrix trace and dagger
//
TEST_CASE("MMatrix:: Dagger and Trace", "[MMatrix]") {

	const double EPS = 1e-5;

	const MMatrix<std::complex<double>> A =
		 { {std::complex<double>(0.4),  std::complex<double>(0.0)},
		   {std::complex<double>(0.0), -std::complex<double>(0.4)}};
	A.Print();

	const MMatrix<std::complex<double>> A_dagger = A.Dagger();
	A_dagger.Print();
	
	const MMatrix<std::complex<double>> AAT = A*A_dagger;
	AAT.Print();

	const double trace = std::real( AAT.Trace() );

	REQUIRE( trace == Approx(0.32).epsilon(EPS) );
}


// Test Lorentz frame transformation functions
//
TEST_CASE("gra::kinematics::LorentzFrame versus PGframe, HXframe, CSframe", "[Lorentz Frames]") {
	
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
		const int direction = -1; // GJ and PG frames

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
	    const double PZ1 = rng.U(50, 7000); // z-momentum
		M4Vec p1; p1.SetPxPyPzM(0,0,PZ1, mproton);

		const double PZ2 = -PZ1;
		M4Vec p2; p2.SetPxPyPzM(0,0,PZ2, mproton);
		
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

			// pb1boost.Print("pb1boost");
			// pb2boost.Print("pb2boost");

			// Pseudo-Gottfried-Jackson frame standalone function test
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

			// Helicity frame standalone function test
			if (frametype[k] == "HX") {
				// printf("HX[0]:      "); pfout[0].Print();
				// printf("HX[1]:      "); pfout[1].Print();				

		  		std::vector<M4Vec> pfHX = pf;
		  		gra::kinematics::HXframe(pfHX, X, DEBUG);
				// printf("HXframe[0]: "); pfHX[0].Print();
				// printf("HXframe[1]: "); pfHX[1].Print();				

				// Difference
				REQUIRE( gra::math::CheckEMC(pfout[0] - pfHX[0], EPS) );
				REQUIRE( gra::math::CheckEMC(pfout[1] - pfHX[1], EPS) );	
			}

			// Collins-Soper frame function standalone test
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

