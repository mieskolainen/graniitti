// Testbench 3
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.


#include <catch.hpp>
#include <random>

#include "Graniitti/MRandom.h"
#include "Graniitti/MMath.h"
#include "Graniitti/MMatrix.h"
#include "Graniitti/MKinematics.h"
#include "Graniitti/M4Vec.h"
#include "Graniitti/MDirac.h"

using namespace gra;

using gra::aux::indices;
using gra::math::pow2;


// Test Lorentz frame transformation functions
// 
TEST_CASE("gra::kinematics::LorentzFrame, PGframe, HEframe, CSframe", "[Lorentz Frames]") {
	
	const int N = 1e5; // Number of trials
	const double EPS = 1e-6;
	
    MRandom rng;
    rng.SetSeed(123456);

    const double sqrts = 7000; // GeV

	// Massless beams
	M4Vec p1(0,0, sqrts/2, sqrts/2);
	M4Vec p2(0,0,-sqrts/2, sqrts/2);

	const double mpion = 0.14; // pion
	std::vector<double> m(2, mpion);
	std::vector<M4Vec> pf(2);

	int repeat = 0;
	do {

		std::cout << repeat << " / " << N << std::endl;

		// Random mother mass from the threshold up to 100 GeV
		const double M0 = rng.U(mpion * 2 + 0.01, 100);

		// Random mother 3-momentum
		const double px = rng.U(-100,100);
		const double py = rng.U(-100,100);
		const double pz = rng.U(-100,100);
		
		M4Vec mother;
		mother.SetPxPyPzM(px, py, pz, M0);

		// Random decay
		gra::kinematics::MCW x;
 		x = gra::kinematics::TwoBodyPhaseSpace(mother, M0, m, pf, rng);

		printf("pf[0]: "); pf[0].Print();
		printf("pf[1]: "); pf[1].Print();

		std::cout << " ***************************************************** " << std::endl;

 		// Do Lorentz transforms

		// --------------------------------------------------------------------
		// ** TRANSFORM TO DIFFERENT LORENTZ FRAMES **
		
		M4Vec              pb1boost;  // beam1 particle boosted
		M4Vec              pb2boost;  // beam2 particle boosted
		std::vector<M4Vec> pfboost;   // central particles boosted
		gra::kinematics::LorentFramePrepare(pf, p1, p2, pb1boost, pb2boost, pfboost);

		std::vector<std::string> frametype = {"CS", "HE", "AH", "PG", "SR"};
		const int direction = 1;
		for (const auto &k : indices(frametype)) {

			// Transform using the generic routine
			std::vector<M4Vec> pfout;
			gra::kinematics::LorentzFrame(pfout, pb1boost, pb2boost, pfboost, frametype[k], direction);

			// Pseudo-Gottfried-Jackson frame test
			if (frametype[k] == "PG") {
				printf("PG[0]:      "); pfout[0].Print();
				printf("PG[1]:      "); pfout[1].Print();				

		  		std::vector<M4Vec> pfPG = pf;
		  		gra::kinematics::PGframe(pfPG, direction, p1, p2);
				printf("PGframe[0]: "); pfPG[0].Print();
				printf("PGframe[1]: "); pfPG[1].Print();

				// Difference
				REQUIRE( gra::math::CheckEMC(pfout[0] - pfPG[0], EPS) );
				REQUIRE( gra::math::CheckEMC(pfout[1] - pfPG[1], EPS) );		
			}

			// Helicity frame test
			if (frametype[k] == "HE") {
				printf("HE[0]:      "); pfout[0].Print();
				printf("HE[1]:      "); pfout[1].Print();				

		  		std::vector<M4Vec> pfHE = pf;
		  		gra::kinematics::HEframe(pfHE);
				printf("HEframe[0]: "); pfHE[0].Print();
				printf("HEframe[1]: "); pfHE[1].Print();				

				// Difference
				REQUIRE( gra::math::CheckEMC(pfout[0] - pfHE[0], EPS) );
				REQUIRE( gra::math::CheckEMC(pfout[1] - pfHE[1], EPS) );	
			}

			// Collins-Soper frame test
			if (frametype[k] == "CS") {
				printf("CS[0]:      "); pfout[0].Print();
				printf("CS[1]:      "); pfout[1].Print();				

		  		std::vector<M4Vec> pfCS = pf;
		  		gra::kinematics::CSframe(pfCS);
				printf("CSframe[0]: "); pfCS[0].Print();
				printf("CSframe[1]: "); pfCS[1].Print();
				
				// Difference
				REQUIRE( gra::math::CheckEMC(pfout[0] - pfCS[0], EPS) );
				REQUIRE( gra::math::CheckEMC(pfout[1] - pfCS[1], EPS) );	
			}
		}
		// --------------------------------------------------------------------
		std::cout << " -----------------------------------------------------" << std::endl;


	} while ((++repeat) < N);

}

