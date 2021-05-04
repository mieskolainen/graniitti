// GRANIITTI Testbench 0 (unit tests)
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

using namespace gra;

using gra::aux::indices;
using gra::math::pow2;


gra::M4Vec psum(const std::vector<gra::M4Vec>& p) {
	gra::M4Vec ptot;
	for (const auto& i : indices(p)) { ptot += p[i]; }
	return ptot;
}


// Generic phase space algorithms
// 
// 
bool PhaseSpaceTest(uint N, const M4Vec& mother, double mpart, double& a, double& b, double& c, double& d) {

    MRandom rng;
    rng.SetSeed(123456);

    const double EPS = 1e-4;

	// Test routine
	printf("Phase space volume test \n");

	const double m0 = mother.M();
	
	// Daughter masses
	printf("<m = %0.3f GeV> : ", mpart);


	std::vector<double> m(N, mpart); // n massless particles
	std::vector<M4Vec> p(N);
	
	bool unweight = false;
	uint repeat   = 0;
	gra::kinematics::MCW Xbod;
	gra::kinematics::MCW Nbod;
	gra::kinematics::MCW Rbod;	

	do {
		gra::kinematics::MCW x;
		if (N == 2) {
	 		x = gra::kinematics::TwoBodyPhaseSpace(mother, m0, m, p, rng);
	 		REQUIRE(math::CheckEMC(mother - psum(p), EPS));
			Xbod += x;
		}
		if (N == 3) {
	 		x = gra::kinematics::ThreeBodyPhaseSpace(mother, m0, m, p, unweight, rng);
	 		REQUIRE(math::CheckEMC(mother - psum(p), EPS));
	 		Xbod += x;
		}

		x = gra::kinematics::NBodyPhaseSpace(mother, m0, m, p, unweight, rng);
		REQUIRE(math::CheckEMC(mother - psum(p), EPS));
		Nbod += x;

		x = gra::kinematics::RamboMassive(mother, m0, m, p, rng);
		REQUIRE(math::CheckEMC(mother - psum(p), EPS));
		Rbod += x;

	} while ((++repeat) < 1E4);
	
	if (N == 2) {

		a = Xbod.Integral();
		b = Nbod.Integral();
		c = Rbod.Integral();
		d = gra::kinematics::PS2Massive(m0*m0, m[0]*m[0], m[1]*m[1]);
		printf("N = %d :: M0 = %6.1f :: \t 2BodyPhaseSpace = %0.5E  NBodyPhaseSpace = %0.5E  RamboMassive = %0.5E  \t[PS2Massive = %0.5E] \n", 
				N, m0, a, b, c, d);
	}
	else if (N == 3) {

		a = Xbod.Integral();
		b = Nbod.Integral();
		c = Rbod.Integral();
		d = a;
		printf("N = %d :: M0 = %6.1f :: \t 3BodyPhaseSpace = %0.5E  NBodyPhaseSpace = %0.5E  RamboMassive = %0.5E  \t[PSnMassless = %0.5E] \n", 
				N, m0, a, b, c, d);
	} else {

		a = 0.0;
		b = Nbod.Integral();
		c = Rbod.Integral();
		d = gra::kinematics::PSnMassless(m0*m0, N);
		printf("N = %d :: M0 = %6.1f :: \t\t\t\t\t  NBodyPhaseSpace = %0.5E  RamboMassive = %0.5E  \t[PSnMassless = %0.5E] \n",
			    N, m0, b, c, d);
	}
	std::cout << std::endl;

	return true;
}

// Generic phase space algorithms
//
//
TEST_CASE("gra::kinematics::TwoBodyPhaseSpace, ThreeBodyPhaseSpace, NBodyPhaseSpace, RamboMassive", "[PhaseSpace]") {

    MRandom rng;
    rng.SetSeed(123456);

	const double EPS = 0.05;

	double  a = 0.0;
	double  b = 0.0;
	double  c = 0.0;
	double  d = 0.0;	

	// Sample different mother momentum
	int trials = 0;
	while (trials < 30) {
		++trials;

		// --------------------------------------------------
		// Random momentum for the mother

		const double Px = rng.U(-5, 5);
		const double Py = rng.U(-5, 5);
		const double Pz = rng.U(-5, 5);
		const double M  = rng.U(1.5, 30);

		M4Vec mother;
		mother.SetPxPyPzM(Px,Py,Pz,M);

		// --------------------------------------------------

		double mpart = 0.0;

		//SECTION("N = 2 | m = 0") {

			PhaseSpaceTest(2, mother, mpart, a,b,c,d);
			REQUIRE( a == Approx(d).epsilon(EPS) );
			REQUIRE( b == Approx(d).epsilon(EPS) );
			REQUIRE( c == Approx(d).epsilon(EPS) );	
		//}
		//SECTION("N = 3 | m = 0") {

			PhaseSpaceTest(3, mother, mpart, a,b,c,d);
			REQUIRE( a == Approx(d).epsilon(EPS) );
			REQUIRE( b == Approx(d).epsilon(EPS) );
			REQUIRE( c == Approx(d).epsilon(EPS) );	
		//}
		//SECTION("N = 4 | m = 0") {

			PhaseSpaceTest(4, mother, mpart, a,b,c,d);
			REQUIRE( b == Approx(d).epsilon(EPS) );
			REQUIRE( c == Approx(d).epsilon(EPS) );
		//}
		
		//SECTION("N = 5 | m = 0") {

			PhaseSpaceTest(5, mother, mpart, a,b,c,d);
			REQUIRE( b == Approx(d).epsilon(EPS) );
			REQUIRE( c == Approx(d).epsilon(EPS) );
		//}
		/*
		SECTION("N = 6 | m = 0") {

			PhaseSpaceTest(6, mother, mpart, a,b,c,d);
			REQUIRE( b == Approx(d).epsilon(EPS) );
			REQUIRE( c == Approx(d).epsilon(EPS) );
		}
		*/

		// ===================================================================
		mpart = 0.139;


		//SECTION("N = 2 | m != 0") {

			PhaseSpaceTest(2, mother, mpart, a,b,c,d);
			REQUIRE( b == Approx(a).epsilon(EPS) );
			REQUIRE( c == Approx(a).epsilon(EPS) );
		//}
		//SECTION("N = 3 | m != 0") {

			PhaseSpaceTest(3, mother, mpart, a,b,c,d);
			REQUIRE( b == Approx(a).epsilon(EPS) );
			REQUIRE( c == Approx(a).epsilon(EPS) );	
		//}
		//SECTION("N = 4 | m != 0") {

			PhaseSpaceTest(4, mother, mpart, a,b,c,d);
			REQUIRE( b == Approx(c).epsilon(EPS) );
		//}
		
		//SECTION("N = 5 | m != 0") {

			PhaseSpaceTest(5, mother, mpart, a,b,c,d);
			REQUIRE( b == Approx(c).epsilon(EPS) );
		//}
		/*
		SECTION("N = 6 | m != 0") {

			PhaseSpaceTest(6, mother, mpart, a,b,c,d);
			REQUIRE( b == Approx(c).epsilon(EPS) );
		}
		*/
	}
}

// Basic operations of 4-vectors
//
//
TEST_CASE("M4Vec: Basic kinematic operations", "[M4Vec]") {

	const double EPS = 1e-5;
	const double mpi = 0.139; // GeV

	M4Vec p(0,0,0, mpi);
	SECTION("Invariant operators") {
		REQUIRE( p.M()  == Approx(mpi).epsilon(EPS));
		REQUIRE( p.M2() == Approx(mpi*mpi).epsilon(EPS));	
	}

	M4Vec k(1.0,2.0,3.0, mpi);
	SECTION("3-Momentum components") {
		REQUIRE( k.Px() == Approx(1.0).epsilon(EPS));
		REQUIRE( k.Py() == Approx(2.0).epsilon(EPS));
		REQUIRE( k.Pz() == Approx(3.0).epsilon(EPS));
	}

	SECTION("Bracket indexing: k[1] == k.Px()") {
		REQUIRE( k[1] == Approx(k.Px()).epsilon(EPS));
		REQUIRE( k[2] == Approx(k.Py()).epsilon(EPS));
		REQUIRE( k[3] == Approx(k.Pz()).epsilon(EPS));
		REQUIRE( k[0] == Approx(k.E()).epsilon(EPS));
	}

	SECTION("Covariant indexing with %% operator: ") {
		for (std::size_t i = 1; i <= 3; ++i) {
		REQUIRE( (k%i) == Approx(-k[i]).epsilon(EPS));	
		}
		REQUIRE( (k%0) == Approx(k[0]).epsilon(EPS));
	}
}

// Factorials
//
//
TEST_CASE("gra::math::factorial: with values of 0,1,2,3,10", "[gra::math::factorial]") {
	
    REQUIRE( gra::math::factorial(0) == 1 );
    REQUIRE( gra::math::factorial(1) == 1 );
    REQUIRE( gra::math::factorial(2) == 2 );
    REQUIRE( gra::math::factorial(3) == 6 );
    REQUIRE( gra::math::factorial(10) == 3628800 );
}

// Legendre polynomials
//
//
TEST_CASE("gra::math::sf_legendre: numerically l = 0 ... 8, m = 0", "[gra::math::sf_legendre]") {

	const double EPS = 1e-5;

	std::vector<double> costheta = {-0.7, -0.3, 0.0, 0.3, 0.7};

	for (std::size_t i = 0; i < costheta.size(); ++i) {
		for (int l = 0; l < 8; ++l) {
			REQUIRE( gra::math::sf_legendre(l, 0, costheta[i]) == Approx(gra::math::LegendrePl(l, costheta[i])).epsilon(EPS) );
		}
	}
}

// Complex spherical harmonics
//
//
TEST_CASE("gra::math::Y_complex_basis: numerically", "[gra::math::Y_complex_basis]") {

	const double EPS = 1e-5;

	std::vector<double> costheta = {-0.7, -0.3, 0.0, 0.3, 0.7};
	std::vector<double> phi      = {-2.5, 2.0, 1.0, 1.5, 2.5};	

	for (std::size_t i = 0; i < costheta.size(); ++i) {
	for (std::size_t j = 0; j < phi.size(); ++j) {
		for (int l = 0;  l <= 4; ++l) {
		for (int m = -l; m <= l; ++m) {
			
			SECTION("l = " + std::to_string(l) + " m = " + std::to_string(m)) {
				REQUIRE( std::real(gra::math::Y_complex_basis(costheta[i], phi[j], l, m)) == Approx(std::real(gra::math::Y_complex_basis_ref(costheta[i], phi[j], l, m))).epsilon(EPS) );
				REQUIRE( std::imag(gra::math::Y_complex_basis(costheta[i], phi[j], l, m)) == Approx(std::imag(gra::math::Y_complex_basis_ref(costheta[i], phi[j], l, m))).epsilon(EPS) );
			}
		}
		}
	}}
}


// Matrix initialization
//
//
TEMPLATE_TEST_CASE("MMatrix:: Initialization with initialization list", "[MMatrix][template]", int) {

	const MMatrix<TestType> A = {{1, 2, 3},
		        			     {4, 5, 6},
						         {7, 8, 9}};

	SECTION("Two-fold initialization list") {

		REQUIRE( A[0][1] == 2);
		REQUIRE( A[2][0] == 7);
		REQUIRE( A[1][2] == 6);
	}

	const std::vector<TestType> row0 = {1, 2, 3};
	const std::vector<TestType> row1 = {4, 5, 6};
	const std::vector<TestType> row2 = {7, 8, 9};
	const MMatrix<TestType> B = {row0, row1, row2};

	SECTION("One-fold initialization list with vector") {

		REQUIRE( B[0][1] == 2);
		REQUIRE( B[2][0] == 7);
		REQUIRE( B[1][2] == 6);
	}

	const MMatrix<TestType> C = A * B;

	SECTION("Matrix x Matrix multiplication") {

		REQUIRE( C[0][0] == 30);
		REQUIRE( C[0][1] == 36);
		REQUIRE( C[2][1] == 126);
	}

	MMatrix<TestType> D = {{1,2,3,4},
						   {5,6,7,8}};
	D = D.Transpose();

	SECTION("Matrix Transpose") {

		REQUIRE( D[0][1] == 5);
		REQUIRE( D[1][1] == 6);
		REQUIRE( D[2][1] == 7);
		REQUIRE( D[3][1] == 8);	
	}
}


