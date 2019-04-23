// GRANIITTI - Monte Carlo event generator for high energy diffraction
// https://github.com/mieskolainen/graniitti
//
// <Kuramoto-Sivashinsky (KS) PDE "snakes" / FFT testbench>
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
#include "Graniitti/MFFT.h"
#include "Graniitti/MH1.h"
#include "Graniitti/MMath.h"
#include "Graniitti/MTimer.h"

// Libraries
#include "rang.hpp"

// Eigen
#include <Eigen/unsupported/Eigen/FFT>

using namespace gra;

// ---------------------------------------------------------------
// Solve Kuramoto-Sivashinsky (KS)
// 4-th order partial differential equation
// u_t + uu_x + u_xx + u_xxxx = 0 on [-pi L,pi L]
//
// via Crank-Nicolson/Adams-Bashforth (CNAB2) iteration
//
// We follow MATLAB/Julia implementation from:
// http://online.kitp.ucsb.edu/online/transturb17/gibson/html/5-kuramoto-sivashinksy.html
//
void CNAB2(std::valarray<std::complex<double>>& u, int Lx, double dt, int Nt, int nplot,
		   bool USE_EIGEN) {
	// number of gridpoints
	const int Nx = u.size();

	// Wavenumbers: exp(2*pi*i*kx*x/L)
	std::valarray<std::complex<double>> kx(Nx);
	unsigned int j = 0;
	for(int i = 0; i <= Nx / 2 - 1; ++i) {
		kx[j] = i;
		++j;
	}
	kx[j] = 0.0;
	++j;
	for(int i = -Nx / 2 + 1; i <= -1; ++i) {
		kx[j] = i;
		++j;
	}

	// Wavenumbers: exp(i*alpha*x)
	std::valarray<std::complex<double>> alpha =
		std::complex<double>(2.0 * gra::math::PI / Lx, 0.0) * kx;
	gra::math::PrintArray(alpha, "alpha");

	// Fourier space: D = d/dx operator
	std::valarray<std::complex<double>> D = std::complex<double>(0, 1) * alpha;
	gra::math::PrintArray(D, "D");

	// Fourier space: -1/2 D operator
	std::valarray<std::complex<double>> G = std::complex<double>(-0.5, 0.0) * D;
	gra::math::PrintArray(G, "G");

	// Fourier space: linear operator -D^2 - D^4
	std::valarray<std::complex<double>> L(alpha.size());
	for(std::size_t i = 0; i < L.size(); ++i) {
		L[i] = alpha[i] * alpha[i] - alpha[i] * alpha[i] * alpha[i] * alpha[i];
	}
	gra::math::PrintArray(L, "L");

	// Create x-vector
	std::valarray<std::complex<double>> x(Nx);
	for(int i = 0; i < (int)x.size(); ++i) {
		x[i] = i * Lx / static_cast<double>(Nx);
	}
	gra::math::PrintArray(x, "x");

	// Time step variables
	const std::complex<double> dt2 = dt / 2;
	const std::complex<double> dt32 = 3 * dt / 2;

	std::valarray<std::complex<double>> A = 1.0 + dt2 * L;
	std::valarray<std::complex<double>> B = 1.0 / (1.0 - dt2 * L);
	gra::math::PrintArray(A, "A");
	gra::math::PrintArray(B, "B");

	// ------------------------------------------------------------------
	Eigen::FFT<double> fft;
	std::vector<std::complex<double>> timevec;
	std::vector<std::complex<double>> freqvec;

	std::valarray<std::complex<double>> uDOTu = u * u;
	if(USE_EIGEN) {
		timevec = gra::math::valarray2vector(uDOTu);
		fft.fwd(freqvec, timevec);
		uDOTu = gra::math::vector2valarray(freqvec);
	} else {
		MFFT::fft(uDOTu);
	}
	std::valarray<std::complex<double>> Nn = G * uDOTu;
	std::valarray<std::complex<double>> Nn1 = Nn;
	gra::math::PrintArray(Nn, "Nn");

	// To spectral domain
	// FFT
	if(USE_EIGEN) {
		timevec = gra::math::valarray2vector(u);
		fft.fwd(freqvec, timevec);
		u = gra::math::vector2valarray(freqvec);
	} else {
		MFFT::fft(u);
	}

	gra::math::PrintArray(u, "fft(u)");
	std::valarray<std::complex<double>> temp;

	// timestepping loop, index for saving
	for(int n = 0; n < Nt; ++n) {
		// IFFT
		if(USE_EIGEN) {
			freqvec = gra::math::valarray2vector(Nn);
			fft.inv(timevec, freqvec);
			Nn = gra::math::vector2valarray(timevec);
		} else {
			MFFT::ifft(Nn);
		}

		// Square
		for(std::size_t k = 0; k < Nn.size(); ++k) {
			Nn[k] = gra::math::pow2(std::real(Nn[k]));
		}

		// To spectral domain
		// FFT
		if(USE_EIGEN) {
			timevec = gra::math::valarray2vector(Nn);
			fft.fwd(freqvec, timevec);
			Nn = gra::math::vector2valarray(freqvec);
		} else {
			MFFT::fft(Nn);
		}

		// Fourier space: Nn = -u u_x
		Nn = G * Nn;
		u = B * (A * u + dt32 * Nn - dt2 * Nn1);

		// Non-Linear term shifted in time: N^{n-1} <- N^n
		Nn1 = Nn;
		Nn = u;

		// Plot
		if((n % nplot) == 0) {
			std::valarray<std::complex<double>> temp2 = u;

			// IFFT
			if(USE_EIGEN) {
				freqvec = gra::math::valarray2vector(temp2);
				fft.inv(timevec, freqvec);
				temp2 = gra::math::vector2valarray(timevec);
			} else {
				MFFT::ifft(temp2);
			}

			const std::vector<std::string> asciiart = {" ", ".", ":", "-", "=",
													   "+", "*", "#", "%", "@"};
			double maxval = 0;
			for(unsigned int i = 0; i < temp2.size(); ++i) {
				maxval = gra::math::abs2(temp2[i]) > maxval ? gra::math::abs2(temp2[i]) : maxval;
			}
			std::cout << n << " ";
			for(unsigned int i = 0; i < temp2.size(); ++i) {
				const double w = gra::math::abs2(temp2[i]) / (maxval + 1e-9);
				unsigned int ind = std::round(w * 9);
				if(ind < 0)
					ind = 0;
				if(ind > 9)
					ind = 9;
				std::cout << rang::fg::red << asciiart[ind];
			}
			std::cout << rang::fg::reset << std::endl;
			//++np;
		}
	}
}

void testfft() {
	// FFT test
	std::valarray<float> xxx =
		gra::math::linspace<std::valarray>((float)0.0, (float)10.0, (unsigned int)16);

	std::valarray<std::complex<double>> X(xxx.size());
	for(unsigned int i = 0; i < xxx.size(); ++i) {
		X[i] = xxx[i];
	}
	MFFT::fft(X);
	std::cout << "FFT:" << std::endl;
	for(unsigned int i = 0; i < X.size(); ++i) {
		std::cout << X[i] << std::endl;
	}
	MFFT::ifft(X);
	std::cout << "IFFT:" << std::endl;
	for(unsigned int i = 0; i < X.size(); ++i) {
		std::cout << X[i] << std::endl;
	}
}

// Main
int main(int argc, char* argv[]) {
	if(argc != 7) {
		std::cout << "Kuramoto-Sivashinsky 4th order PDE FFT-testbench" << std::endl;
		std::cout << "Usage: ./pdebench 64 128 16 4 2400 0" << std::endl;
		return EXIT_FAILURE;
	}

	// Test FFT
	testfft();

	int Lx = std::atoi(argv[1]);
	int Nx = std::atoi(argv[2]);
	double dt = 1.0 / static_cast<double>(std::atoi(argv[3]));
	int nplot = std::atoi(argv[4]);
	int Nt = std::atoi(argv[5]);
	bool USE_EIGEN = std::atoi(argv[6]);

	std::valarray<std::complex<double>> x(Nx);
	for(int i = 0; i < (int)x.size(); ++i) {
		x[i] = i * Lx / static_cast<double>(Nx);
	}
	std::valarray<std::complex<double>> u(Nx);
	for(int i = 0; i < (int)u.size(); ++i) {
		// Initialize with single wave
		// u[i] = std::cos(x[i]/(double)Lx)*(1.0 + std::sin(x[i]/(double)Lx));

		// Initialize with a canonic one
		u[i] = std::cos(x[i]) + 0.1 * std::cos(x[i] / 16.0) * (1.0 + 2.0 * std::sin(x[i] / 16.0));
	}
	gra::math::PrintArray(x, "x");
	gra::math::PrintArray(u, "u");

	// Solve the equation
	CNAB2(u, Lx, dt, Nt, nplot, USE_EIGEN);

	return EXIT_SUCCESS;
}
