// GRANIITTI - Monte Carlo event generator for high energy diffraction
// https://github.com/mieskolainen/graniitti
//
// <Path integral CPU performance benchmark>
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <cassert>
#include <complex>
#include <future>
#include <iostream>
#include <random>
#include <stdexcept>
#include <valarray>
#include <vector>

// Own
#include "Graniitti/MAux.h"
#include "Graniitti/MH1.h"
#include "Graniitti/MMath.h"
#include "Graniitti/MTimer.h"

// Libraries
#include "cxxopts.hpp"
#include "json.hpp"
#include "rang.hpp"

using gra::aux::indices;
using namespace gra;

// Simulation variables
int N = 0;
double dt = 0;
double m = 0;
int k_slit = 0;
double SLIT = 0;
unsigned long NS = 0;

// Potential
double V(double x) {
	// return gra::aux::pow2(x)/2.0; // Harmonic (quadratic)
	return 0.0; // no external potential
}

// 1D-lattice action
double S_lat(std::vector<double>& x) {
	double S = 0.0;
	for(int j = 0; j < N - 1; ++j) {
		S += (m / (2.0 * dt)) * gra::math::pow2(x[j + 1] - x[j]) - dt * V(x[j]);
	}
	return S;
}

// Print variables
void printvar() {
	printf("  N = %d, dt = %0.3f, m = %0.3f, k_slit = %d, SLIT = %0.3f, samples "
		   "= %lu \n",
		   N, dt, m, k_slit, SLIT, NS);
}

// Main
int main(int argc, char* argv[]) {
	// Detector observable histogram
	const double XWIDTH = 18;
	MH1<std::complex<double>> hx(50, -XWIDTH, XWIDTH, "Path Integral MC");

	// Seed here
	std::random_device rd;
	// Standard mersenne_twister_engine seeded with rd()
	std::mt19937_64 gen(rd());

	// double offset = 2.0;

	if(argc != 7) { // Default

		N = 4;
		dt = 0.8;
		m = 1;
		k_slit = 2;
		SLIT = 0.03;
		NS = 10000000;

		printf("Example input ./pathmark %d %0.3f %0.3f %d %0.3f %lu \n", N, dt, m, k_slit, SLIT,
			   NS);
		std::cout << std::endl;
	} else {
		// Collect input variables
		N = atoi(argv[1]);
		dt = atof(argv[2]);
		m = atof(argv[3]);
		k_slit = atoi(argv[4]); // Slit position
		SLIT = atof(argv[5]);
		NS = atoi(argv[6]);

		printvar();
	}

	// Random variables
	std::uniform_real_distribution<> dis(-XWIDTH, XWIDTH);
	std::uniform_real_distribution<> slit(-SLIT, SLIT);
	std::uniform_real_distribution<> outL(-XWIDTH, 0 - SLIT);
	std::uniform_real_distribution<> outR(0 + SLIT, XWIDTH);

	// Imaginary unit
	const std::complex<double> ii(0, 1);

	// Path configuration vector
	std::vector<double> x(N, 0.0);

	MTimer globaltimer(true);
	MTimer localtimer(true);

	std::vector<std::future<std::complex<double>>> futures; // std::async return values

	// MC samples
	for(unsigned long i = 0; i < NS; ++i) {
		// Draw a single path configuration
		for(int k = 1; k < N; ++k) {
			x[k] = dis(gen);
		}

		// Apply double slit boundary condition to two time steps
		// (one is not enough)
		x[k_slit] = slit(gen);
		x[k_slit - 1] = slit(gen);

		// ** Get complex action weight **
		std::complex<double> W = std::exp(ii * S_lat(x));

		// Fill histogram
		hx.Fill(x[N - 1], W);

		if(localtimer.ElapsedSec() > 0.1) {
			localtimer.Reset();
			gra::aux::PrintProgress(i / static_cast<double>(NS));
		}
	}
	gra::aux::ClearProgress(); // Clear progressbar

	// Print histogram
	hx.Print();

	double runtime = globaltimer.ElapsedSec();
	std::cout << rang::style::bold << "<pathmark - complex path integral CPU benchmark>"
			  << rang::style::reset << std::endl
			  << std::endl;
	printf("- Runtime         = %0.2E sec\n", runtime);
	printf("- Action integral = %0.3E eval/sec\n", NS / runtime);
	printf("\n");

	printvar();
	std::cout << std::endl;

	return EXIT_SUCCESS;
}
