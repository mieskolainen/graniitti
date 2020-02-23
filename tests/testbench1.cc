// Testbench 1 (unit tests)
//
// (c) 2017-2020 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include <catch.hpp>
#include <random>

#include "Graniitti/MMath.h"
#include "Graniitti/MMatrix.h"
#include "Graniitti/MKinematics.h"
#include "Graniitti/M4Vec.h"
#include "Graniitti/MRandom.h"
#include "Graniitti/MPDG.h"

// LHAPDF
#include "LHAPDF/LHAPDF.h"

// LUXqed / LHAPDF read
//
//
TEST_CASE("LUXqed parton fractions", "[LUXqed]") {

	const double EPS = 0.03; // Assert relative error criteria

	std::unique_ptr<LHAPDF::PDF> GlobalPdfPtr = nullptr;

	std::string pdfname = "LUXqed17_plus_PDF4LHC15_nnlo_100";
	try {
		GlobalPdfPtr = std::unique_ptr<LHAPDF::PDF>(LHAPDF::mkPDF(pdfname, 0));
	} catch (...) {
		std::string str = "TEST_CASE::LUXqed: Problem with reading a pdfset '" + pdfname + "'";
		throw std::invalid_argument(str);
	}
	
	// Evaluate pdf values for all supported parton PIDs (return to map)
	// pids: -5,-4,-3,-2,-1, 1, 2, 3, 4, 5, 21, 22
	std::vector<double> Q2val = {10*10, 100*100, 1000*1000, 10000*10000}; // pdf factorization scale

	for (const auto& Q2 :  Q2val) {

	std::cout << "Q   = " << sqrt(Q2) << " GeV" << std::endl;
	std::cout << "----------------------------" << std::endl;
	std::cout << "(1 ~ d, 2 ~ u, 3 ~ s, 4 ~ c, 5 ~ b, 21 ~ g, 22 ~ gamma)" << std::endl;
	std::cout << std::endl;

	gra::MRandom rng;
	std::map<int,double> sums;
	const int N = 1e6;

	std::map<int,double> rtn;

	for (int i = 0; i < N; ++i) {

		double x = rng.U(0,1);
		GlobalPdfPtr->xfxQ2(x, Q2, rtn);
		
		// Collect pdf values
		for (const auto& [key, val] : rtn) {
			sums[key] += val;
		}
	}
	double total = 0;
	GlobalPdfPtr->xfxQ2(0.5, Q2, rtn);

	for (const auto& [key,val] : sums) {
		printf("%2d  = %0.5f \n", key, val / N);
		total += val / N;
	}
	printf("sum = %0.5f \n\n", total);

	REQUIRE(total == Approx(1.0).epsilon(EPS));
	}
}

// Drell-Yan test
//
//
TEST_CASE("Drell-Yan u u~ -> mu+ mu- & LHAPDF test", "[Drell-Yan]") {

	std::unique_ptr<LHAPDF::PDF> GlobalPdfPtr = nullptr;

	const double EPS = 0.5; // Assert relative error criteria

	if (GlobalPdfPtr == nullptr) {
		std::string pdfname = "LUXqed17_plus_PDF4LHC15_nnlo_100";
		//std::string pdfname = "CT10nlo";
		try {
			GlobalPdfPtr = std::unique_ptr<LHAPDF::PDF>(LHAPDF::mkPDF(pdfname, 0));
		} catch (...) {
			std::string str = "TEST_CASE::Drell-Yan: Problem with reading a pdfset '" + pdfname + "'";
			throw std::invalid_argument(str);
		}
	}
	using namespace gra;
	using gra::math::pow2;
	using gra::math::pow4;

	// Constants
	const double alpha = 1.0/137.0;                   // Take it at q^2 = 0
	//const double   m_e = 0.0005109989461;           // Electron mass
	const double  m_mu = 0.1056583745;                // Muon mass
	const double     e = sqrt(alpha * 4 * math::PI);
	const double  pico = 1e12;

	// CMS energy^2
	const double s = 13000*13000;

	MRandom rng;
	double sumw = 0.0;	
	const int N = 1e5;

	// generate u u~ > mu+ mu- / z
	std::vector<int> flavor = {2}; // 1=down, 2=up, 3=strange, 4=charm

	// Final states
	std::vector<M4Vec> pf;
	const std::vector<double> mf = {m_mu, m_mu};

	for (std::size_t f = 0; f < flavor.size(); ++f) {
	for (std::size_t i = 0; i < N; ++i) {

		const double x1 = rng.U(0,1);
		const double x2 = rng.U(0,1);
		
		// We-work in CMS frame
		// Initial state (px,py,pz,E)
		M4Vec p(0, 0,  x1 * sqrt(s)/2, x1 * sqrt(s)/2);
		M4Vec k(0, 0, -x2 * sqrt(s)/2, x2 * sqrt(s)/2);

		// Propagator
		M4Vec q = p + k;

		// Factorization scale
		const double Q2 = q.M2();

		// Divide x out
		const double f1 = GlobalPdfPtr->xfxQ2( flavor[f], x1, Q2) / x1;
		const double f2 = GlobalPdfPtr->xfxQ2(-flavor[f], x2, Q2) / x2;	

		// Generate final state kinematics
		gra::kinematics::MCW PS = gra::kinematics::TwoBodyPhaseSpace(q, q.M(), mf, pf, rng);

		// s-channel amplitude squared
		// 1/4 \sum_{spins} |M_s|_i^2
		//
		// 1/3 for averaging over initial state colors
		// (2/3)^2 for up type charge
		const double Amp2 = (1.0/3.0) * pow2(2.0/3.0) * 32.0 * pow4(e) / (4.0 * pow2(q.M2())) *
							( (k * pf[0])*(p * pf[1]) + (p * pf[0])*(k * pf[1]) );

		// Total event weight
		if (q.M() > 60) // invariant mass cut
		sumw += PS.GetW() * Amp2 * f1 * f2 / (2.0*Q2) * PDG::GeV2barn;
	}
	}
	std::cout << std::endl;

	printf("Check the MadGraph factorization scale Q^2 definition - big impact here with LO-matrix element \n");

	const double integral_mc       = sumw / N * pico;
	const double integral_madgraph = 19.17; // pb, u u~ > mu+ mu- MadGraph, ct10nlo, M > 60, 13 TeV
	printf("sqrts = %0.1f GeV \t Integral = %0.10f pb \t MadGraph = %0.10f pb \n",
			sqrt(s), integral_mc, integral_madgraph);

	REQUIRE(integral_mc == Approx(integral_madgraph).epsilon(EPS));
}


// generate a a > e+ e- / z
// Compare with MadGraph (remember alpha_qed scale setup in MadGraph)
// generate e+ e- > mu+ mu- / z
//
TEST_CASE("gra::M4Vec & gra::kinematics::TwoBodyPhaseSpace & e+e- -> mu+mu- tree-level", "[M4Vec & TwoBodyPhaseSpace]") {

	using namespace gra;
	using gra::math::pow2;
	using gra::math::pow4;

	const double EPS = 1e-3; // Assert relative error criteria

	// Analytic integrated cross section
	auto xs_analytic = [&] (double s, double alpha, double Nc, double Qf) {
		return Nc * pow2(Qf) * 4.0 * math::PI * pow2(alpha) / (3*s) * PDG::GeV2barn;
	};

	// Constants
	const double alpha = 1.0/137.0;                   // Take it at q^2 = 0
	const double   m_e = 0.0005109989461;             // Electron mass
	const double  m_mu = 0.1056583745;                // Muon mass
	const double     e = sqrt(alpha * 4 * math::PI);
	const double pico = 1e12;

	// CMS-energies
	const std::vector<double> sqrt_s = {1000,13000};

	for (const auto & i : gra::aux::indices(sqrt_s)) {

		// CMS energy^2
		const double s = pow2(sqrt_s[i]);

		// Moller flux (high energy limit)
		const double Flux = 2*s;

		// We-work in CMS frame
		// Initial state (px,py,pz,E)
		M4Vec p(0, 0,  sqrt(s/4 - pow2(m_e)), sqrt(s)/2);
		M4Vec k(0, 0, -sqrt(s/4 - pow2(m_e)), sqrt(s)/2);

		// Propagator
		M4Vec q = p + k;

		// Final states
		std::vector<M4Vec> pf;
		const std::vector<double> mf = {m_mu, m_mu};

		// Number of MC samples
		const int N = 1e6;
		double sumw = 0.0;
		MRandom rng;
		for (int i = 0; i < N; ++i) {

			// Generate final state kinematics
			gra::kinematics::MCW PS = gra::kinematics::TwoBodyPhaseSpace(q, q.M(), mf, pf, rng);

			// s-channel amplitude squared
			// 1/4 \sum_{spins} |M_s|_i^2
			const double Amp2 = 32.0 * pow4(e) / (4.0 * pow2((p+k).M2())) *
								( (k * pf[0])*(p * pf[1]) + (p * pf[0])*(k * pf[1]) );

			// Total event weight
			sumw += PS.GetW() * Amp2 / Flux * PDG::GeV2barn;
		}
		
		double integral_mc       = sumw / N;
		double integral_analytic = xs_analytic(s, alpha, 1, 1);
		printf("sqrts = %0.1f GeV \t Integral = %0.10f pb \t Analytic = %0.10f pb \n",
			sqrt(s), integral_mc * pico, integral_analytic * pico);

		REQUIRE(integral_mc == Approx(integral_analytic).epsilon(EPS));
	}
}
