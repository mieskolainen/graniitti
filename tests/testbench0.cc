// Testbench 0
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

using namespace gra;

using gra::aux::indices;
using gra::math::pow2;


// Phase space functions
bool PhaseSpaceTest(uint N, double m0, double& a, double& b, double& c) {
	
    MRandom rng;
    rng.SetSeed(123456);

	// Test routine
	printf("Phase space volume test \n");
	M4Vec mother(0,0,0, m0);

	// Daughter masses
	double mpart = 0.0;
	if (N == 2 || N == 3) {
		mpart = 0.0;
	}
	printf("<m = %0.3f GeV> : ", mpart);

	std::vector<double> m(N, mpart); // n massless particles
	std::vector<M4Vec> p(N);

	bool unweight = true;
	uint repeat   = 0;
	gra::kinematics::MCW Xbod;
	gra::kinematics::MCW Nbod;

	do {
		gra::kinematics::MCW x;
		if (N == 2) {
	 		x = gra::kinematics::TwoBodyPhaseSpace(mother, m0, m, p, rng);
			Xbod += x;
		}
		if (N == 3) {
	 		x = gra::kinematics::ThreeBodyPhaseSpace(mother, m0, m, p, unweight, rng);
	 		Xbod += x;
		}
		x = gra::kinematics::NBodyPhaseSpace(mother, m0, m, p, unweight, rng);
		Nbod += x;

	} while ((++repeat) < 1E4);

	double exact = gra::kinematics::PSnMassless(m0*m0, N);

	if (N == 2) {
		exact = gra::kinematics::PS2Massive(m0*m0, m[0]*m[0], m[1]*m[1]);
	}
	if (N == 2 || N == 3) {

		a = Xbod.Integral();
		b = Nbod.Integral();
		c = exact;
		printf("N = %d :: M0 = %6.1f :: \t XBodyPhaseSpace = %0.5E (%0.5f)  NBodyPhaseSpace = %0.5E (%0.5f) \t[PSnMassless = %0.5E] \n", 
				N, m0, a, a/c, b, b/c, c);
	} else {

		a = 0.0;
		b = Nbod.Integral();
		c = exact;
		printf("N = %d :: M0 = %6.1f :: \t\t\t\t\t\t  NBodyPhaseSpace = %0.5E (%0.5f) \t[PSnMassless = %0.5E] \n",
			    N, m0, b/c, b/c, c);
	}
	std::cout << std::endl;

	return true;
}


TEST_CASE("gra::kinematics::TwoBodyPhaseSpace, ThreeBodyPhaseSpace, NBodyPhaseSpace", "[PhaseSpace]") {

	double M0 = 1500;
	double  a = 0.0;
	double  b = 0.0;
	double  c = 0.0;

	const double EPS = 0.01;

	SECTION("N = 2") {

		PhaseSpaceTest(2, M0, a,b,c);
		REQUIRE( a == Approx(c).epsilon(EPS) );
		REQUIRE( b == Approx(c).epsilon(EPS) );
	}
	SECTION("N = 3") {

		PhaseSpaceTest(3, M0, a,b,c);
		REQUIRE( a == Approx(c).epsilon(EPS) );
		REQUIRE( b == Approx(c).epsilon(EPS) );
	}
	SECTION("N = 4") {

		PhaseSpaceTest(4, M0, a,b,c);
		REQUIRE( b == Approx(c).epsilon(EPS) );
	}
	SECTION("N = 5") {

		PhaseSpaceTest(5, M0, a,b,c);
		REQUIRE( b == Approx(c).epsilon(EPS) );
	}
	SECTION("N = 6") {

		PhaseSpaceTest(6, M0, a,b,c);
		REQUIRE( b == Approx(c).epsilon(EPS) );
	}
}


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


TEST_CASE("gra::math::factorial: with values of 0,1,2,3,10", "[gra::math::factorial]") {
	
    REQUIRE( gra::math::factorial(0) == 1 );
    REQUIRE( gra::math::factorial(1) == 1 );
    REQUIRE( gra::math::factorial(2) == 2 );
    REQUIRE( gra::math::factorial(3) == 6 );
    REQUIRE( gra::math::factorial(10) == 3628800 );
}


TEST_CASE("gra::math::sf_legendre: numerically l = 0 ... 8, m = 0", "[gra::math::sf_legendre]") {

	const double EPS = 1e-5;

	std::vector<double> costheta = {-0.7, -0.3, 0.0, 0.3, 0.7};

	for (std::size_t i = 0; i < costheta.size(); ++i) {
		for (int l = 0; l < 8; ++l) {
			REQUIRE( gra::math::sf_legendre(l, 0, costheta[i]) == Approx(gra::math::LegendrePl(l, costheta[i])).epsilon(EPS) );
		}
	}
}


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


