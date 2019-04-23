// "Brute force" Dirac Gamma Algebra for amplitudes
//
// Next step, optimize this class using more efficient spinor-products.
// For those, see e.g.:
//
// http://www.phys.nthu.edu.tw/~class/Helicity/YPYAO_Note_01.pdf
// http://www.physics.wayne.edu/~ablechman/main/Research_files/shm.pdf
// R. Kleiss, W.J. Stirling, Spinor Techniques, Nucl.Phys. B262 (1985) 235-262
// https://arxiv.org/pdf/hep-ph/0110108.pdf
//
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <complex>
#include <iostream>
#include <vector>

// Tensor algebra
#include "FTensor.hpp"

// Own
#include "Graniitti/M4Vec.h"
#include "Graniitti/MAux.h"
#include "Graniitti/MDirac.h"
#include "Graniitti/MKinematics.h"
#include "Graniitti/MMath.h"
#include "Graniitti/MSpin.h"

using FTensor::Tensor1;
using FTensor::Tensor2;

using gra::math::msqrt;
using gra::math::pow2;
using gra::math::zi;

namespace gra {
MDirac::MDirac() {
	// Default
	InitGammaMatrices("DIRAC");
}

void MDirac::InitGammaMatrices(const std::string& basis) {
	// Set gamma basis
	if(basis == "DIRAC") {
		BASIS = "D";
	} else if(basis == "CHIRAL") {
		BASIS = "C";
	} else {
		throw std::invalid_argument("MDirac::InitGammaMatrices: Unknown gamma basis chosen: " +
									basis);
	}

	// PARITY OPERATOR = y^0
	// Intrinsic parity of fermions +1, anti-fermion -1

	const MMatrix<std::complex<double>> y0_chiral{
		std::vector<std::complex<double>>{0.0, 0.0, 1.0, 0.0},
		std::vector<std::complex<double>>{0.0, 0.0, 0.0, 1.0},
		std::vector<std::complex<double>>{1.0, 0.0, 0.0, 0.0},
		std::vector<std::complex<double>>{0.0, 1.0, 0.0, 0.0}};

	const MMatrix<std::complex<double>> y0_dirac{
		std::vector<std::complex<double>>{1.0, 0.0, 0.0, 0.0},
		std::vector<std::complex<double>>{0.0, 1.0, 0.0, 0.0},
		std::vector<std::complex<double>>{0.0, 0.0, -1.0, 0.0},
		std::vector<std::complex<double>>{0.0, 0.0, 0.0, -1.0}};

	// \equiv i \gamma^0\gamma^1\gamma^2\gamma^3
	const MMatrix<std::complex<double>> y5_chiral{
		std::vector<std::complex<double>>{-1.0, 0.0, 0.0, 0.0},
		std::vector<std::complex<double>>{0.0, -1.0, 0.0, 0.0},
		std::vector<std::complex<double>>{0.0, 0.0, 1.0, 0.0},
		std::vector<std::complex<double>>{0.0, 0.0, 0.0, 1.0}};

	// \equiv i \gamma^0\gamma^1\gamma^2\gamma^3
	const MMatrix<std::complex<double>> y5_dirac{
		std::vector<std::complex<double>>{0.0, 0.0, 1.0, 0.0},
		std::vector<std::complex<double>>{0.0, 0.0, 0.0, 1.0},
		std::vector<std::complex<double>>{1.0, 0.0, 0.0, 0.0},
		std::vector<std::complex<double>>{0.0, 1.0, 0.0, 0.0}};

	// ------------------------------------------------------------------
	// Both chiral and dirac basis
	// Contravariant matrices (upper index)

	MMatrix<std::complex<double>> y1_up{std::vector<std::complex<double>>{0.0, 0.0, 0.0, 1.0},
										std::vector<std::complex<double>>{0.0, 0.0, 1.0, 0.0},
										std::vector<std::complex<double>>{0.0, -1.0, 0.0, 0.0},
										std::vector<std::complex<double>>{-1.0, 0.0, 0.0, 0.0}};

	MMatrix<std::complex<double>> y2_up{std::vector<std::complex<double>>{0.0, 0.0, 0.0, -zi},
										std::vector<std::complex<double>>{0.0, 0.0, zi, 0.0},
										std::vector<std::complex<double>>{0.0, zi, 0.0, 0.0},
										std::vector<std::complex<double>>{-zi, 0.0, 0.0, 0.0}};

	MMatrix<std::complex<double>> y3_up{std::vector<std::complex<double>>{0.0, 0.0, 1.0, 0.0},
										std::vector<std::complex<double>>{0.0, 0.0, 0.0, -1.0},
										std::vector<std::complex<double>>{-1.0, 0.0, 0.0, 0.0},
										std::vector<std::complex<double>>{0.0, 1.0, 0.0, 0.0}};

	// ------------------------------------------------------------------
	// Contravariant and covariant
	if(BASIS == "D") {
		gamma_up = {y0_dirac, y1_up, y2_up, y3_up, y5_dirac};
		gamma_lo = {y0_dirac, -y1_up, -y2_up, -y3_up, y5_dirac};
	}
	if(BASIS == "C") {
		gamma_up = {y0_chiral, y1_up, y2_up, y3_up, y5_chiral};
		gamma_lo = {y0_chiral, -y1_up, -y2_up, -y3_up, y5_chiral};
	}

	// sigma = i/2 [\gamma^\mu, \gamma^\nu]
	// \sigma_{\mu\nu} each [mu][nu] contains one MMatrix
	sigma_up = std::vector<std::vector<MMatrix<std::complex<double>>>>(
		4, std::vector<MMatrix<std::complex<double>>>(4, MMatrix<std::complex<double>>(4, 4)));
	sigma_lo = std::vector<std::vector<MMatrix<std::complex<double>>>>(
		4, std::vector<MMatrix<std::complex<double>>>(4, MMatrix<std::complex<double>>(4, 4)));

	for(const auto& mu : LI) {
		for(const auto& nu : LI) {
			sigma_up[mu][nu] =
				(gamma_up[mu] * gamma_up[nu] - gamma_up[nu] * gamma_up[mu]) * (zi / 2.0);
			sigma_lo[mu][nu] =
				(gamma_lo[mu] * gamma_lo[nu] - gamma_lo[nu] * gamma_lo[mu]) * (zi / 2.0);
		}
	}
}

// Two component Weyl (chiral) spinor
//
std::vector<std::complex<double>> MDirac::XiSpinor(const M4Vec& p, int helicity) const {
	const double theta2 = p.Theta() / 2.0;
	const double phi = p.Phi();
	const double c = std::cos(theta2);
	const double s = std::sin(theta2);

	if(helicity == 1) {
		return {c, s * std::exp(zi * phi)};
	}
	if(helicity == -1) {
		return {-s * std::exp(-zi * phi), c};
	}
	throw std::invalid_argument("MDirac::XiSpinor: helicity is not -1 or 1");
}

// Helicity eigenstate spinors for particles
//
// <@@ DEFINED IN CHIRAL GAMMA MATRIX REPRESENTATION @@>
//
std::vector<std::complex<double>> MDirac::uHelChiral(const M4Vec& p, int helicity) const {
	if(BASIS != "C") {
		throw std::invalid_argument("MDirac::uHelChiral: Wrong gamma basis in use");
	}

	const double E = p.E();
	const double p3 = p.P3mod();
	const double m = p.M();
	const double theta2 = p.Theta() / 2.0;
	const double phi = p.Phi();

	const double c = std::cos(theta2);
	const double s = std::sin(theta2);
	const std::complex<double> phase = std::exp(zi * phi);

	const double neg = E + m - p3;
	const double pos = E + m + p3;
	const double N = 1.0 / (std::sqrt(2.0) * msqrt(E + m)); // Volume normalization to 2E

	switch(helicity) {
		case 1:
			return {c * neg * N, s * phase * neg * N, c * pos * N, s * phase * pos * N};
		case -1:
			return {-s * pos * N, c * phase * pos * N, -s * neg * N, c * phase * neg * N};
	}

	throw std::invalid_argument("MDirac::uHelChiral: helicity is not -1 or 1");
}

// Helicity eigenstate spinors for Anti-Particles
//
// <@@ DEFINED IN CHIRAL GAMMA MATRIX REPRESENTATION @@>
//
std::vector<std::complex<double>> MDirac::vHelChiral(const M4Vec& p, int helicity) const {
	if(BASIS != "C") {
		throw std::invalid_argument("MDirac::vHelChiral: Wrong gamma basis in use");
	}

	if(helicity != 1 && helicity != -1) {
		throw std::invalid_argument("MDirac::vHelChiral: helicity is not -1 or 1");
	}
	// Flip the helicity, so we can use particle solution permutated
	helicity = -helicity;
	const std::vector<std::complex<double>> v = uHelChiral(p, helicity);
	return {-v[0], -v[1], v[2], v[3]};

	throw std::invalid_argument("MDirac::vHelChiral: helicity is not -1 or 1");
}

// Chirality projectors, valid for all gamma representations
// Right handed
MMatrix<std::complex<double>> MDirac::PR() const {
	return (I4 + gamma_up[4]) * 0.5;
}
// Left handed
MMatrix<std::complex<double>> MDirac::PL() const {
	return (I4 - gamma_up[4]) * 0.5;
}

// Helicity eigenstate spinors
//
// <@@ DEFINED IN DIRAC GAMMA MATRIX REPRESENTATION @@>
//
// [REFERENCE: Thomson, Modern Particle Physics, Cambridge University Press]
// [https://www.hep.phy.cam.ac.uk/~thomson/partIIIparticles/handouts/Handout_2_2011.pdf]
//
// In high energy limit E >> m, helicity eigenstates are eigenstates of gamma^5
// gamma^5 u_+ = +u_+
// gamma^5 u_- = -u_-
// gamma^5 v_+ = -v_+
// gamma^5 v_- = +v_-
//
// Remember: In the limit E >> m (only then)
// -> left and right handed chiral states == helicity states.
//
std::vector<std::complex<double>> MDirac::uHelDirac(const M4Vec& p, int helicity) const {
	if(BASIS != "D") {
		throw std::invalid_argument("MDirac::uHelDirac: Wrong gamma basis in use");
	}

	const double E = p.E();
	const double p3 = p.P3mod();
	const double m = p.M();
	const double N = msqrt(E + m); // Volume normalization to 2E

	const double theta2 = p.Theta() / 2.0;
	const double phi = p.Phi();

	const double c = std::cos(theta2);
	const double s = std::sin(theta2);
	const std::complex<double> phase = std::exp(zi * phi);

	switch(helicity) {
		case 1:
			return {N * c, N * phase * s, N * p3 / (E + m) * c, N * p3 / (E + m) * phase * s};
		case -1:
			return {-N * s, N * phase * c, N * p3 / (E + m) * s, -N * p3 / (E + m) * phase * c};
	}

	throw std::invalid_argument("MDirac::uHelDirac: helicity is not -1 or 1");
}

// Helicity eigenstate spinors for Anti-Fermions
//
// <@@ DEFINED IN DIRAC GAMMA-MATRIX REPRESENTATION @@>
//
std::vector<std::complex<double>> MDirac::vHelDirac(const M4Vec& p, int helicity) const {
	if(BASIS != "D") {
		throw std::invalid_argument("MDirac::vHelDirac: Wrong gamma basis in use");
	}
	if(helicity != 1 && helicity != -1) {
		throw std::invalid_argument("MDirac::vHelDirac: helicity is not -1 or 1");
	}
	// Flip the helicity, so we can use particle solution permutated
	helicity = -helicity;
	const std::vector<std::complex<double>> v = uHelDirac(p, helicity);
	return {v[2], v[3], v[0], v[1]};
}

// Dirac-Spinor "spin eigenstate" fermions: (\gamma^\mu p_\mu - m)u = 0
//
// <@@ DEFINED IN DIRAC GAMMA MATRIX REPRESENTATION @@>
//
std::vector<std::complex<double>> MDirac::uDirac(const M4Vec& p, int spin) const {
	if(BASIS != "D") {
		throw std::invalid_argument("MDirac::uDirac: Wrong gamma basis in use");
	}
	if(spin != -1 && spin != 1) {
		throw std::invalid_argument("MDirac::uDirac: spin state argument != -1 or 1");
	}
	const double E = p.E();
	const double m = p.M();
	const double N = msqrt(E + m); // Normalization to 2E per unit volume

	if(spin == 1) {
		return {N * 1.0, 0.0, N * p[3] / (E + m), N * (p[1] + zi * p[2]) / (E + m)};
	} else { // spin == -1
		return {0.0, N * 1.0, N * (p[1] - zi * p[2]) / (E + m), N * (-p[3] / (E + m))};
	}
}

// Dirac-Spinor "spin eigenstate" for anti-fermions: (\gamma^\mu p_\mu + m)v = 0
//
// <@@ DEFINED IN DIRAC GAMMA MATRIX REPRESENTATION @@>
//
std::vector<std::complex<double>> MDirac::vDirac(const M4Vec& p, int spin) const {
	if(BASIS != "D") {
		throw std::invalid_argument("MDirac::vDirac: Wrong gamma basis in use");
	}
	if(!(spin == -1 || spin == 1)) {
		throw std::invalid_argument("MDirac::vDirac: spin state argument != -1 or 1");
	}

	// Flip the spin, then use the u-particle solution permutated
	spin = -spin;
	const std::vector<std::complex<double>> v = uDirac(p, spin);
	return {v[2], v[3], v[0], v[1]};
}

// Return Hermitician angular momentum operators J_i = 1/2 \sigma_i, i = 1,2,3
MMatrix<std::complex<double>> MDirac::J_operator(unsigned int i) const {
	if(i == 1) {
		return sigma_x * 0.5;
	}
	if(i == 2) {
		return sigma_y * 0.5;
	}
	if(i == 3) {
		return sigma_z * 0.5;
	}

	throw std::invalid_argument("MDirac::J_operator: i is invalid (not 1,2,3)");
}

// Polarization vector eps^{(m)\mu}(k) for massless spin-1 with helicity m =
// -1,1
//
// Spatial dependence only on the direction.
//
// [REFERENCE: http://scipp.ucsc.edu/~haber/ph218/polsum.pdf]
//
// [t,x,y,z] order convention!
//
Tensor1<std::complex<double>, 4> MDirac::EpsSpin1(const M4Vec& k, int m) const {
	const double theta = k.Theta();
	const double phi = k.Phi();

	// const std::vector<double> e1 = {1.0, 0.0, 0.0};
	// const std::vector<double> e2 = {0.0, 1.0, 0.0};

	// First define in \vec{z}-direction, then rotate spatial part such that
	// \vec{k} = R \vec{z}

	// Normalization
	const double N = 1.0 / std::sqrt(2.0);

	// \eps^{\mu,-1}
	if(m == -1) {
		// e = {0.0,  (e1[0] + zi*e2[0]) / msqrt(2.0),  (e1[1] + zi*e2[1]) /
		// msqrt(2.0),
		// 0.0};
		Tensor1<std::complex<double>, 4> e;
		e(0) = 0.0;
		e(1) = (std::cos(theta) * std::cos(phi) + zi * std::sin(phi)) * N;
		e(2) = (std::cos(theta) * std::sin(phi) - zi * cos(phi)) * N;
		e(3) = -std::sin(theta) * N;
		return e;
	}
	// \eps^{\mu,+1}
	if(m == 1) {
		// e = {0.0, -(e1[0] - zi*e2[0]) / msqrt(2.0), -(e1[1] - zi*e2[1]) /
		// msqrt(2.0),
		// 0.0};
		Tensor1<std::complex<double>, 4> e;
		e(0) = 0.0;
		e(1) = -(std::cos(theta) * std::cos(phi) + zi * std::sin(phi)) * N;
		e(2) = -(std::cos(theta) * std::sin(phi) - zi * cos(phi)) * N;
		e(3) = std::sin(theta) * N;
		return e;
	}
	throw std::invalid_argument("MDirac::EpsSpin1: helicity m should be -1 or 1");
}

// Transformation bilinears:
//
// Vector current: \bar{u}_1 \gamma^\mu u_2
// Axial current:  \bar{u}_1 \gamma^\mu \gamma_5 \u_2
// Tensor current: \bar{u}_1 \sigma^{\mu\nu} \gamma_5 \u_2
// Pseudoscalar current:  \bar{u}_1 \gamma_5 \u_2
//
//
// Adjoing gamma matrix: gamma^\mu\dagger = gamma^0 gamma^\mu gamma^0
//
//
// Charge conjugation operator matrix is basis dependent.
// u = C \bar{v}^T <=> \bar{v}^T = C^{-1} u
//
//
// Basis independent definition:
// C^{-1} \gamma_\mu C = -\gamma_\mu^T   <=>  C \gamma^\mu C^{-1} = -\gamma^\mu
// C^\dagger = C^{-1}
// C^T = -C
//
MMatrix<std::complex<double>> MDirac::C_up() const {
	if(BASIS == "D" || BASIS == "C") {
		return -gamma_up[2] * gamma_up[0] * zi;
	} else {
		throw std::invalid_argument("MDirac::COP: Unknown gamma basis");
	}
}

// --------------------------------------------------------------------
// Basic notes:
//
// 1. Helicity is conserved in fermion lines
// 2. By gauge invariance, polarization vectors can be chosen with explicit
//    representation. Invidual sub-amplitudes are gauge dependent, their sum is
//    not.
// 3. Helicity amplitudes do not interfere, sum incoherently \sum_{config}
// |M_i|^2
//
//
// --------------------------------------------------------------------
// Conventions:
//
// Spin-1/2:
// incoming particle:      u(p)
// outgoing particle:      \bar{u}(p)
// incoming anti-particle: \bar{v}(p)
// outgoing anti-particle: v(p)
//
// Propagators:
// spin 1:                 -ig_\mu\nu / q^2
// spin 1/2:                i(\gamma^\mu q_\mu + m)/(q^2 - m^2)
//
// Vertex factor:
//                         ie\gamma^\mu
//
// Matrix element: -iM     product of rules
// --------------------------------------------------------------------

// Photon propagator: -ig_{\mu\nu} / q^2
//
// Remember:
// \sum_{4 virtual polarizations} \eps_\mu^\lambda (\eps_{\nu}^\lambda)^* =
// -g_{\mu\nu}
//
Tensor2<std::complex<double>, 4, 4> MDirac::iD_y(const M4Vec& q) const {
	const double q2 = q.M2();
	Tensor2<std::complex<double>, 4, 4> T;

	for(const auto& u : LI) {
		for(const auto& v : LI) {
			T(u, v) = -zi * g[u][v] / q2;
		}
	}
	return T;
}

// Internal Fermion propagator (matrix with spinor indices, no Lorentz indices)
//
// Input as contravariant (upper) index 4-vector
//
MMatrix<std::complex<double>> MDirac::iD_F(const M4Vec& q, double m) const {
	const double q2 = q.M2();

	// gamma_\mu q^\mu contraction
	MMatrix<std::complex<double>> M(4, 4, 0.0); // Init with zero!

	for(const auto& mu : LI) {
		M += (gamma_lo[mu] * q[mu] + I4 * m);
	}
	M = M * (zi / (q2 - m * m));
	return M;
}

// Polarization vector eps^{(m)\mu}(k) for massive spin-1 with helicity m =
// -1,0,1
//
// m = +-1 state depends only on the direction of the momentum
// m = 0   state depends also on magnitude
//
// [t,x,y,z] order convention!
//
Tensor1<std::complex<double>, 4> MDirac::EpsMassiveSpin1(const M4Vec& k, int m) const {
	// \eps^{(0),-+1} (massless case applies here too)
	if(m == -1 || m == 1) {
		return EpsSpin1(k, m);
	}
	// Should be 0 at this point
	if(m != 0) {
		throw std::invalid_argument("MDirac::EpsMassiveSpin1: helicity m should be -1, 0 or 1");
	}

	// \eps^{(0),\mu} (longitudinal case)
	const double theta = k.Theta();
	const double phi = k.Phi();
	const double E = k.E();
	const double M = k.M();

	Tensor1<std::complex<double>, 4> e;
	e(0) = k.P3mod() / M;
	e(1) = E * std::sin(theta) * std::cos(phi) / M;
	e(2) = E * std::sin(theta) * std::sin(phi) / M;
	e(3) = E * std::cos(theta) / M;

	return e;
}

// Polarization tensor eps^{(m)\mu\nu}(k) for massive spin-2 with helicity m =
// -2,-1,0,1,2
//
// [t,x,y,z] convention!
//
// Should obey: (eps_{\mu\nu}^{(m)}(k))^* eps^{(n)\mu\nu}(k) ) = \delta_{mn}
//
Tensor2<std::complex<double>, 4, 4> MDirac::EpsMassiveSpin2(const M4Vec& k, int m) const {
	if(!(m == -2 || m == -1 || m == 0 || m == 1 || m == 2)) {
		throw std::invalid_argument("MDirac::EpsMassiveSpin2: m is not -2,-1,0,1,2");
	}
	Tensor2<std::complex<double>, 4, 4> epsmat;

	const int offset = 1; // array indexing (-1 means 0)

	// Get Massive Spin-1 basis vectors
	std::array<Tensor1<std::complex<double>, 4>, 3> epsvec;
	for(const int& m : {-1, 0, 1}) {
		epsvec[m + offset] = EpsMassiveSpin1(k, m);
	}

	// Loop over two Lorentz indices
	for(const auto& mu : LI) {
		for(const auto& nu : LI) {
			// Clebsch-Gordan decomposition
			for(const int& m1 : {-1, 0, 1}) {
				for(const int& m2 : {-1, 0, 1}) {
					// ClebschGordan(double j1, double j2, double m1, double m2,
					// double j, double m)
					// Get coefficient <1,1, m1, m2| 2,m>
					epsmat(mu, nu) += gra::spin::ClebschGordan(1.0, 1.0, m1, m2, 2.0, m) *
									  epsvec[m1 + offset](mu) * epsvec[m2 + offset](nu);
				}
			}
		}
	}
	return epsmat;
}

// Adjoint Dirac spinor: \bar{u} = u^dagger * gamma^0
std::vector<std::complex<double>>
	MDirac::Bar(const std::vector<std::complex<double>>& spinor) const {
	// First conjugate elements, then a matrix product with gamma^0 matrix
	return gra::matoper::VecMatMultiply(gra::matoper::VecDagger(spinor), gamma_up[0]);
}

// Feynman slash matrix operator: \slash{a} = \gamma_\mu a^\mu = \gamma^\mu
// a_\mu
//
// Input assumed contravariant (upper) index 4-vector
MMatrix<std::complex<double>> MDirac::FSlash(const M4Vec& a) const {
	MMatrix<std::complex<double>> aslash(4, 4, 0.0); // Init with zero!
	for(const auto& mu : LI) {
		aslash += gamma_up[mu] * (a % mu);
	}
	return aslash;
}

// ----------------------------------------------------------------------
// High Energy (E >> m) relation for the vertex current when p = p' (e.g.
// forward scattering)
//
// \bar{u}_s'(p) \gamma^\mu u_s(p) = 2p^\mu \delta_{s',s}
// \bar{v}_s'(p) \gamma^\mu v_s(p) = 2p^\mu \delta_{s',s}
//
// Other relation:
// \bar{u}_s'(p) \gamma^0 v_s(-p) = 0
// \bar{v}_s'(p) \gamma^0 u_s(-p) = 0
//
bool MDirac::SpinorHELimit(const M4Vec& pi, const M4Vec& pf) const {
	std::cout << "MDirac::SpinorHELimit:" << std::endl;

	const M4Vec psum = pi + pf;

	// Helicities
	for(const auto& hi : {1, 2}) {
		const std::vector<std::complex<double>> u = uDirac(pi, hi);

		// Helicities
		for(const auto& hf : {1, 2}) {
			const std::vector<std::complex<double>> ubar = Bar(uDirac(pf, hf));

			if(hi != hf) {
				std::cout << std::endl;
				continue;
			} // Helicity conservation

			for(const auto& mu : LI) {
				for(const auto& nu : LI) {
					const std::vector<std::complex<double>> prod = (gamma_up[mu] * (psum % nu)) * u;

					const std::complex<double> lhs = gra::matoper::VecVecMultiply(ubar, prod);
					const double rhs = (psum % mu) * (psum % nu) * Delta(hi, hf);

					const double absratio = std::abs(std::real(lhs)) / std::abs(rhs);

					if(absratio < 0.9 || absratio > 1.1) {
						std::cout << rang::fg::red;
					} else {
						std::cout << rang::fg::green;
					}

					printf(" hi = %2d, hf = %2d | (mu,nu) = (%lu,%lu) | Re[lhs] = "
						   "%0.3E \t rhs = %0.3E \t |lhs| / |rhs| = %0.5f \n",
						   hi, hf, mu, nu, std::real(lhs), rhs, absratio);

					std::cout << rang::fg::reset;
				}
			}
		}
	}
	std::cout << std::endl << std::endl;
	return true;
}

// Test gamma matrix anticommutation relation {\gamma^\mu, \gamma^nu} = 2
// g^{\mu\nu} I_4
bool MDirac::GammaAntiCommutation() const {
	for(const auto& mu : LI) {
		for(const auto& nu : LI) {
			const MMatrix<std::complex<double>> AC_lo =
				gamma_lo[mu] * gamma_lo[nu] + gamma_lo[nu] * gamma_lo[mu];
			std::cout << "gamma_lo:: mu:" << mu << " nu: " << nu << std::endl;
			AC_lo.Print();

			const MMatrix<std::complex<double>> AC_up =
				gamma_up[mu] * gamma_up[nu] + gamma_up[nu] * gamma_up[mu];
			std::cout << "gamma_up:: mu:" << mu << " nu: " << nu << std::endl;
			AC_up.Print();
		}
	}
	return true;
}

// ----------------------------------------------------------------------
// [REFERENCE: https://arxiv.org/pdf/hep-ph/0110108.pdf]

// Spinor product: s_\lambda(p1,p2)
std::complex<double> MDirac::sProd(const M4Vec& p1, const M4Vec& p2, int helicity) const {
	return gra::matoper::VecVecMultiply(Bar(uGauge(p1, helicity)), uGauge(p2, -helicity));
}

// Helicity u-spinor via massless gauge vector
//
std::vector<std::complex<double>> MDirac::uGauge(const M4Vec& p, int helicity) const {
	const M4Vec l(100, 0, 0, 100);
	return ((FSlash(p) + I4 * p.M()) / msqrt(2.0 * (p * l))) * uHelDirac(l, -helicity);
}

// Helicity u-spinor via massless gauge vector
//
std::vector<std::complex<double>> MDirac::vGauge(const M4Vec& p, int helicity) const {
	const M4Vec l(100, 0, 0, 100);
	return (-(FSlash(p) - I4 * p.M()) / msqrt(2.0 * (p * l))) * vHelDirac(l, -helicity);
}

// ----------------------------------------------------------------------
// Test functions

// Test \slash{p}\slash{p} = p^2 I_4
double MDirac::FSlashFSlash(const M4Vec& p) const {
	const MMatrix<std::complex<double>> A = FSlash(p) * FSlash(p);
	const MMatrix<std::complex<double>> B = I4 * p.M2();

	const double norm = (A - B).FrobNorm();
	std::cout << "FSlashFSlash:: Frobenius norm |A-B|_F: " << norm << std::endl;
	return norm;
}

// Check Dirac u- and v-spinor normalization and completeness relations
//
// 1. Normalization of the Lorentz invariant inner product:
//
//   \bar{u}_s'(p) u_s(p) = +2m \delta_{s's}
//   \bar{v}_s'(p) v_s(p) = -2m \delta_{s's}
//   \bar{u}_s'(p) v_s(p) = 0
//   \bar{v}_s'(p) u_s(p) = 0
//
// Note that \dagger{u}_s'(p) u_s(p) = 2E \delta_{s's}
// (standard relativistic normalization of a wavefunction)
//
//
// 2. Completeness relation (Projection operators):
//
//   \sum_{s=+-} u_s(p) \bar{u}_s(p) = \slash{p} + m
//   \sum_{s=+-} v_s(p) \bar{v}_s(p) = \slash{p} - m
//
//
// Set type = "u" or "v"
//
// Normalization test (add to unit tests!)
//
// SPINOR      CHIRAL  |  DIRAC
// ------------------------------------------
// uDirac    | -       |  OK
// vDirac    | -       |  OK
// uHelDirac | -       |  OK
// vHelDirac | -       |  OK
// uHel      | OK      |  -
// vHel      | OK      |  -
//
//
bool MDirac::DiracSpinorComplete(const M4Vec& p, const std::string& type,
								 const std::string& basis) const {
	std::cout << "DiracSpinorComplete:: Type: " << type << std::endl;
	// InitGammaMatrices(basis);
	MMatrix<std::complex<double>> lhs(4, 4, 0.0); // Init with zero!
	const double SIGN = ((type == "u") ? 1.0 : -1.0);

	for(const auto& lambda : SPINORSTATE) {
		std::vector<std::complex<double>> spinor;
		if(type == "u") {
			// spinor = uDirac(p,lambda);
			// spinor = uHelDirac(p,lambda);

			spinor = uHelChiral(p, lambda);
			// spinor = uGauge(p,lambda);
		} else if(type == "v") {
			// spinor = vDirac(p,lambda);
			// spinor = vHelDirac(p,lambda);

			spinor = vHelChiral(p, lambda);
			// spinor = vGauge(p,lambda);
		} else {
			throw std::invalid_argument("DiracSpinorComplete: Unknown type (set u or v)");
		}
		// Adjoint
		const std::vector<std::complex<double>> spinorbar = Bar(spinor);

		// Check normalization
		const double nlhs = std::real(
			gra::matoper::VecVecMultiply(spinorbar, spinor)); // Take real to cast to double
		const double nrhs = SIGN * 2. * p.M();

		printf("s = %2d : Normalization = ", lambda);
		if(gra::aux::AssertRatio(nlhs, nrhs, 5E-3)) {
			std::cout << rang::fg::green << nlhs << " | 2 x p.M = " << nrhs << " OK"
					  << rang::fg::reset << std::endl;
		} else {
			std::cout << rang::fg::red << nlhs << " | 2 x p.M = " << nrhs << " NOT OK"
					  << rang::fg::reset << std::endl;
		}

		// Take outerproduct, sum
		lhs += gra::matoper::OuterProd(spinor, spinorbar);
	}
	std::cout << std::endl;

	// Completeness relation
	const MMatrix<std::complex<double>> rhs = FSlash(p) + I4 * p.M() * SIGN;

	// Compare
	lhs.Print("sum_{s}{spinor_s(p) barspinor_s(p)}");
	std::cout << std::endl;
	rhs.Print("slash{p} + I4*m*SIGN");
	std::cout << std::endl;

	const double frobnorm = (lhs - rhs).FrobNorm();
	if(frobnorm < 0.1) {
		std::cout << rang::fg::green << "Completeness relation OK, Frobenius norm = " << frobnorm
				  << rang::fg::reset << std::endl;
	} else {
		std::cout << rang::fg::red << "Completeness relation NOT OK, Frobenius norm = " << frobnorm
				  << rang::fg::reset << std::endl;
	}
	std::cout << std::endl;

	return true;
}

// Check massive spin1-polarization sum (completeness relation)
//
// \sum_{\lambda = -1,0.,1} eps^\mu_\lambda eps^{*\mu}_\lambda = -g^{\mu\nu} +
// k^\mu k^\nu / M^2
//
bool MDirac::MassiveSpin1Complete(const M4Vec& k) const {
	for(const auto& mu : LI) {
		for(const auto& nu : LI) {
			// Massive spin-1 helicities
			std::complex<double> sum = 0.;
			for(const int& lambda : {-1, 0, 1}) {
				const Tensor1<std::complex<double>, 4> eps = EpsMassiveSpin1(k, lambda);
				sum += eps(mu) * std::conj(eps(nu));
			}
			// Right hand side
			const double rhs = -g[mu][nu] + k[mu] * k[nu] / k.M2();
			printf("CheckEps1Complete:: lhs = %0.2E, rhs = %0.2E \n", std::real(sum), rhs);
		}
	}
	return true;
}

} // gra namespace ends
