// Testbench 2 (unit tests)
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

// Dirac algebra tests
//
//
TEST_CASE("MDirac::TestGammaAntiCommutation", "[MDirac]") {
	const double EPS = 1e-5;
	MDirac dirac;
	REQUIRE(dirac.TestGammaAntiCommutation() < EPS);
}


// Dirac algebra tests
// 
// 
TEST_CASE("MDirac multiple tests", "[MDirac]") {

	const double EPS = 1E-5;

	// First create some random particles
	MRandom rng;
	rng.SetSeed(123456);

	const double m0 = 100;
	const double mdaughter = 0.139;
	const int N = 2;

	M4Vec mother(0,0,0, m0);
	std::vector<double> m(N, mdaughter);
	
	// Create random 2-body decays
	for (std::size_t trial = 0; trial < 10; ++trial) {
		std::vector<M4Vec> p(N);
		gra::kinematics::TwoBodyPhaseSpace(mother, m0, m, p, rng);

		// Different basis
		std::vector<std::string> BASIS = {"DIRAC", "CHIRAL"};
		for (std::size_t k = 0; k < BASIS.size(); ++k) {

			MDirac dirac(BASIS[k]);
			
			// For each particle, do tests
			for (std::size_t i = 0; i < N; ++i) {

	  			//REQUIRE(dirac.TestSpinorHELimit(const M4Vec &p1, const M4Vec &p2));
				
				std::vector<std::string> MODE = {"helicity", "gauge", "spin"};
				for (const auto& mode : MODE) {
	  				REQUIRE(dirac.TestSpinorComplete(p[i], "u", mode) < EPS);
	  				REQUIRE(dirac.TestSpinorComplete(p[i], "v", mode) < EPS);
	  			}

		    	REQUIRE(dirac.TestMassiveSpin1Complete(p[i]) < EPS);
				REQUIRE(dirac.TestFSlashFSlash(p[i]) < EPS);
			}
		}
	}
}

