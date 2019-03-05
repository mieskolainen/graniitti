// Tensor Pomeron amplitudes
// 
// See:
// [REFERENCE: Ewerz, Maniatis, Nachtmann, https://arxiv.org/pdf/1309.3478.pdf]
// [REFERENCE: Lebiedowicz, Nachtmann, Szczurek, https://arxiv.org/pdf/1801.03902.pdf]
// [REFERENCE: Lebiedowicz, Nachtmann, Szczurek, https://arxiv.org/pdf/1901.11490.pdf] 
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <iostream>
#include <vector>
#include <complex>

// Own
#include "Graniitti/MTensorPomeron.h"
#include "Graniitti/MDirac.h"
#include "Graniitti/MKinematics.h"
#include "Graniitti/M4Vec.h"
#include "Graniitti/MMath.h"
#include "Graniitti/MForm.h"
#include "Graniitti/MRegge.h"

// FTensor algebra
#include "FTensor.hpp"


using gra::math::msqrt;
using gra::math::pow2;
using gra::math::zi;
using gra::math::PI;


using FTensor::Tensor0;
using FTensor::Tensor1;
using FTensor::Tensor2;
using FTensor::Tensor3;
using FTensor::Tensor4;


namespace gra {


// 2 -> 4 amplitudes
std::complex<double> MTensorPomeron::ME4(gra::LORENTZSCALAR& lts) const {

	// Kinematics
	const M4Vec pa = lts.pbeam1;
	const M4Vec pb = lts.pbeam2;

	const M4Vec p1 = lts.pfinal[1];
	const M4Vec p2 = lts.pfinal[2];
	const M4Vec p3 = lts.decaytree[0].p4;
	const M4Vec p4 = lts.decaytree[1].p4;
	
	// Intermediate boson/fermion mass
	const double M_mes     = lts.decaytree[0].p.mass;
	const double Gamma_mes = lts.decaytree[0].p.width;

	// Momentum convention of sub-diagrams
	// 
	// ------<------ anti-particle p3
	// |
	// | \hat{t} (arrow down)
	// |
	// ------>------ particle p4
	//
	// ------<------ particle p4
	// |
	// | \hat{u} (arrow up)
	// |
	// ------>------ anti-particle p3
	//
	const M4Vec pt = pa - p1 - p3; // => q1 = pt + p3, q2 = pt - p4
	const M4Vec pu = p4 - pa + p1; // => q2 = pu + p3, q1 = pu - p3

	// ------------------------------------------------------------------

	// Spinors (2 helicities)
	const std::array<std::vector<std::complex<double>>, 2> u_a      = SpinorStates(pa, "u");
	const std::array<std::vector<std::complex<double>>, 2> u_b      = SpinorStates(pb, "u");

	const std::array<std::vector<std::complex<double>>, 2> ubar_1   = SpinorStates(p1, "ubar");
	const std::array<std::vector<std::complex<double>>, 2> ubar_2   = SpinorStates(p2, "ubar");

	// Massive polarization vectors (3 helicities)
	const std::array<Tensor1<std::complex<double>,4>, 3> eps_3_conj = Spin1States(p3, "conj");
	const std::array<Tensor1<std::complex<double>,4>, 3> eps_4_conj = Spin1States(p4, "conj");

	// ------------------------------------------------------------------
	
	// t-channel Pomeron propagators
	const Tensor4<std::complex<double>,4,4,4,4>  iDP_13  = iD_P(lts.ss[1][3], lts.t1);
	const Tensor4<std::complex<double>,4,4,4,4>  iDP_24  = iD_P(lts.ss[2][4], lts.t2);

	// u-channel Pomeron propagators
	const Tensor4<std::complex<double>,4,4,4,4>  iDP_14  = iD_P(lts.ss[1][4], lts.t1);
	const Tensor4<std::complex<double>,4,4,4,4>  iDP_23  = iD_P(lts.ss[2][3], lts.t2);

	// ------------------------------------------------------------------

	// Fermion propagator
	const double m = lts.decaytree[0].p.mass;
	const MMatrix<std::complex<double>> iSF_t = iD_F(pt, m);
	const MMatrix<std::complex<double>> iSF_u = iD_F(pu, m);

	// Central spinors (2 helicities)
	const std::array<std::vector<std::complex<double>>, 2> v_3    = SpinorStates(p3, "v");
	const std::array<std::vector<std::complex<double>>, 2> ubar_4 = SpinorStates(p4, "ubar");

	// Reset
	lts.hamp.clear();


	// Free Lorentz indices [second parameter denotes the range of index]
	FTensor::Index<'a', 4> mu1;
	FTensor::Index<'b', 4> nu1;
	FTensor::Index<'c', 4> rho1;
	FTensor::Index<'d', 4> rho2;
	FTensor::Index<'g', 4> alpha1;
	FTensor::Index<'h', 4> beta1;
	FTensor::Index<'i', 4> alpha2;
	FTensor::Index<'j', 4> beta2;
	FTensor::Index<'k', 4> mu2;
	FTensor::Index<'l', 4> nu2;
	
	
	// DEDUCE mode
	uint SPINMODE = 0;
	if      (lts.decaytree[0].p.spinX2 == 0 && lts.decaytree[1].p.spinX2 == 0) {
		SPINMODE = 0;
	}
	else if (lts.decaytree[0].p.spinX2 == 1 && lts.decaytree[1].p.spinX2 == 1) {
		SPINMODE = 1;
	}
	else if (lts.decaytree[0].p.spinX2 == 2 && lts.decaytree[1].p.spinX2 == 2) {
		SPINMODE = 2;
	} else {
		throw std::invalid_argument("MTensorPomeron::ME4: Invalid daughter spin (J=0,1/2,1 pairs supported)");
	}


	// Two helicity states for incoming and outgoing protons
	for (const auto& ha : {0,1}) {
	for (const auto& hb : {0,1}) {
	for (const auto& h1 : {0,1}) {
	for (const auto& h2 : {0,1}) {

		// Apply proton leg helicity conservation / No helicity flip (high energy limit)
		// This gives at least 4 x speed improvement
		if (ha != h1 || hb != h2) { continue; }

		// Full proton-Pomeron-proton spinor structure (upper and lower vertex)
		//const Tensor2<std::complex<double>,4,4> iG_1a = iG_Ppp(p1, pa, ubar_1[h1], u_a[ha]);
		//const Tensor2<std::complex<double>,4,4> iG_2b = iG_Ppp(p2, pb, ubar_2[h2], u_b[hb]);

		// High Energy Limit proton-Pomeron-proton spinor structure
		const Tensor2<std::complex<double>,4,4> iG_1a = iG_PppHE(p1,pa);
		const Tensor2<std::complex<double>,4,4> iG_2b = iG_PppHE(p2,pb);

		// ==============================================================
		// 2 x pseudoscalar
		if (SPINMODE == 0) {

			// t-channel blocks
			const Tensor2<std::complex<double>,4,4>  iG_ta   = iG_Ppsps(pt, -p3);
			const std::complex<double>               iDMES_t = iD_MES0(pt, M_mes);
			const Tensor2<std::complex<double>,4,4>  iG_tb   = iG_Ppsps(p4,  pt);

			// u-channel blocks
			const Tensor2<std::complex<double>,4,4>  iG_ua   = iG_Ppsps(p4,  pu);
			const std::complex<double>               iDMES_u = iD_MES0(pu, M_mes);
			const Tensor2<std::complex<double>,4,4>  iG_ub   = iG_Ppsps(pu, -p3);


			std::complex<double> M_t;
			std::complex<double> M_u;
			// t-channel
			{
				// Upper block
				const std::complex<double> A =
				                iG_1a(mu1,nu1) *
				       		    iDP_13(mu1,nu1, alpha1, beta1) *
				       	 	    iG_ta(alpha1, beta1) *
								iG_tb(alpha2, beta2) *
								iDP_24(alpha2,beta2, mu2, nu2) *
								iG_2b(mu2, nu2);
				// Apply off-shell form factors
				M_t = A * iDMES_t * pow2(PARAM_REGGE::Meson_FF(pt.M2(), pow2(M_mes) ));
			}

			// u-channel
			{
				// Upper block
				const std::complex<double> A =
								iG_1a(mu1,nu1) *
				       		    iDP_14(mu1,nu1, alpha1, beta1) *
				       	 	    iG_ua(alpha1, beta1) *
								iG_ub(alpha2, beta2) *
								iDP_23(alpha2,beta2, mu2, nu2) *
								iG_2b(mu2, nu2);
				// Apply off-shell form factors
				M_u = A * iDMES_u * pow2(PARAM_REGGE::Meson_FF(pu.M2(), pow2(M_mes) ));
			}

			// Full amplitude: iM = [ ... ]  <-> M = (-i)*[ ... ]
			std::complex<double> M = (-zi)*(M_t + M_u);
			lts.hamp.push_back(M);
		}
		
		
		// ==============================================================
		// 2 x fermion
		if (SPINMODE == 1) {

			for (const auto& h3 : {0,1}) {
			for (const auto& h4 : {0,1}) {

				// Fermion propagator and connected parts
				const Tensor4<std::complex<double>,4,4,4,4> iG_t = iG_PppbarP(p4, ubar_4[h4], pt, iSF_t, v_3[h3], -p3);
				const Tensor4<std::complex<double>,4,4,4,4> iG_u = iG_PppbarP(p4, ubar_4[h4], pu, iSF_u, v_3[h3], -p3);

				std::complex<double> M_t;
				std::complex<double> M_u;
				// t-channel
				{
					// Upper block
					const std::complex<double> A =
					                iG_1a(mu1,nu1) *
					       		    iDP_13(mu1,nu1, alpha1, beta1) *
					       	 	    iG_t(alpha2, beta2, alpha1, beta1) *
									iDP_24(alpha2,beta2, mu2, nu2) * 
									iG_2b(mu2, nu2);
					// Apply off-shell form factors
					M_t = A * pow2(PARAM_REGGE::Baryon_FF(pt.M2(), pow2(M_mes) ));
				}

				// u-channel
				{
					// Upper block
					const std::complex<double> A =
									iG_1a(mu1,nu1) *
					       		    iDP_14(mu1,nu1, alpha1, beta1) *
					       	 	    iG_u(alpha1, beta1, alpha2, beta2) *
									iDP_23(alpha2,beta2, mu2, nu2) *
									iG_2b(mu2, nu2);
					// Apply off-shell form factors
					M_u = A * pow2(PARAM_REGGE::Baryon_FF(pu.M2(), pow2(M_mes) ));
				}

				// Full amplitude: iM = [ ... ]  <-> M = (-i)*[ ... ]
				std::complex<double> M = (-zi)*(M_t + M_u);
				lts.hamp.push_back(M);
			}}
		}


		// ==============================================================
		// 2 x vector meson
		if (SPINMODE == 2) {

			// t-channel blocks
			const Tensor4<std::complex<double>,4,4,4,4>  iG_tA   = iG_Pvv(pt, -p3);
			const Tensor2<std::complex<double>,4,4>      iDMES_t = iD_VMES(pt, M_mes, Gamma_mes);
			const Tensor4<std::complex<double>,4,4,4,4>  iG_tB   = iG_Pvv(p4, pt);

			// u-channel blocks
			const Tensor4<std::complex<double>,4,4,4,4>  iG_uA   = iG_Pvv(p4, pu);
			const Tensor2<std::complex<double>,4,4>      iDMES_u = iD_VMES(pu, M_mes, Gamma_mes);
			const Tensor4<std::complex<double>,4,4,4,4>  iG_uB   = iG_Pvv(pu, -p3);

			FTensor::Index<'e', 4> rho3;
			FTensor::Index<'f', 4> rho4;

			// -------------------------------------------------------------------------
			// TENSOR LORENTZ INDEX CONTRACTION BLOCK
			// This is done in pieces, due to template <> argument deduction constraints

			Tensor2<std::complex<double>, 4,4> M_t;
			Tensor2<std::complex<double>, 4,4> M_u;

			// t-channel
			{
				// Upper block
				Tensor2<std::complex<double>, 4,4> A;
				A(rho3, rho1) = iG_1a(mu1,nu1) *
				       		    iDP_13(mu1,nu1, alpha1, beta1) *
				       	 	    iG_tA(rho3,rho1, alpha1, beta1);
				// Lower block
				Tensor2<std::complex<double>, 4,4> B;
				B(rho2, rho4) = iG_2b(mu2, nu2) *
								iDP_24(alpha2,beta2, mu2, nu2) *
								iG_tB(rho2,rho4, alpha2, beta2);
				// Apply off-shell form factors
				M_t(rho3,rho4) = A(rho3,rho1) * iDMES_t(rho1,rho2) * B(rho2,rho4) * pow2(PARAM_REGGE::Meson_FF(pt.M2(), pow2(M_mes)) );
			}

			// u-channel
			{
				// Upper block
				Tensor2<std::complex<double>, 4,4> A;
				A(rho4, rho1) = iG_1a(mu1,nu1) *
				       		    iDP_14(mu1,nu1, alpha1, beta1) *
				       	 	    iG_uA(rho4,rho1, alpha1, beta1);
				// Lower block
				Tensor2<std::complex<double>, 4,4> B;
				B(rho2, rho3) = iG_2b(mu2, nu2) *
								iDP_23(alpha2,beta2, mu2, nu2) *
								iG_uB(rho2,rho3, alpha2, beta2);
				// Apply off-shell form factors
				M_u(rho3,rho4) = A(rho4,rho1) * iDMES_u(rho1,rho2) * B(rho2,rho3) * pow2(PARAM_REGGE::Meson_FF(pu.M2(), pow2(M_mes)) );
			}

			// Total amplitude: iM = [ ... ]  <-> M = (-i)*[ ... ]
			Tensor2<std::complex<double>, 4,4> M;
			for (const auto& mu : LI) {
			for (const auto& nu : LI) {
				M(mu,nu) = (-zi)*(M_t(mu,nu) + M_u(mu,nu));
			}}

			const int OPTION = 1;
			
			// Polarization sum
			if (OPTION == 1) {

				double Amp2 = 0.0;
				for (const auto& sigma3 : LI) {
				for (const auto& sigma4 : LI) {
				for (const auto& rho3   : LI) {
				for (const auto& rho4   : LI) {

					const std::complex<double> contract = 
						std::conj(M(sigma3,sigma4)) * M(rho3,rho4) * g[sigma3][rho3] * g[sigma4][rho4];

						Amp2 += std::real(contract); // real for casting to double
				}}}}

				lts.hamp.push_back(msqrt(Amp2));
			}
			// --------------------------------------------------

			// Explicit polarization vectors
			if (OPTION == 2) {

				// Loop over massive Spin-1 helicity states	(-1,0,1)
				for (const auto& h3 : {0,1,2}) {
				for (const auto& h4 : {0,1,2}) {

					// Contract Lorentz indices to get the helicity amplitude (ha,hb,h1,h2,h3,h4)
					const std::complex<double> amp =
						eps_3_conj[h3](rho3) * eps_4_conj[h4](rho4) * M(rho3,rho4);

					lts.hamp.push_back(amp);
				}}
			}
		}
	}}}}

	// Get total amplitude squared 1/4 \sum_h |A_h|^2
	double SumAmp2 = 0.0;
	for (std::size_t i = 0; i < lts.hamp.size(); ++i) {
		SumAmp2 += gra::math::abs2(lts.hamp[i]);
	}
	SumAmp2 /= 4; // Initial state helicity average

	return msqrt(SumAmp2); // We expect amplitude
}



// Construct spin-1/2 helicity spinors (-1,1) [indexing with 0,1]
std::array<std::vector<std::complex<double>>, 2> MTensorPomeron::SpinorStates(const M4Vec& p, std::string type) const {

	std::array<std::vector<std::complex<double>>, 2> spinor;

	for (const auto& m : {0,1}) {
		if (type == "u")    { spinor[m] =     MDirac::uHelDirac(p, SPINORSTATE[m]); }
		if (type == "ubar") { spinor[m] = Bar(MDirac::uHelDirac(p, SPINORSTATE[m])); }

		if (type == "v")    { spinor[m] =     MDirac::vHelDirac(p, SPINORSTATE[m]); }
		if (type == "vbar") { spinor[m] = Bar(MDirac::vHelDirac(p, SPINORSTATE[m])); }	
	}
	return spinor;
}


// Construct Massive Spin-1 polarization vectors (m = -1,0,1) [indexing with 0,1,2]
std::array<Tensor1<std::complex<double>,4>, 3> MTensorPomeron::Spin1States(const M4Vec& p, std::string type) const {

	std::array<Tensor1<std::complex<double>,4>, 3> eps;
	const int lambda[3] = {-1,0,1};

	for (const auto& m : {0,1,2}) { // loop over helicities
		eps[m] = MDirac::EpsMassiveSpin1(p, lambda[m]);
		
		if (type == "conj") {     // Take complex conjugate per element
			for (const auto& mu : LI) {
				eps[m](mu) = std::conj(eps[m](mu));
			}
		}
	}
	return eps;
}



// -------------------------------------------------------------------------
// Vertex functions


// Gamma-Lepton-Lepton vertex
// contracted in \bar{spinor} G_\mu \bar{spinor}
//
// iGamma_\mu (p',p)
//
// High Energy Limit:
// ubar(p',\lambda') \gamma_\mu u(p,\lambda ~= (p' + p)_\mu \delta_{\lambda',\lambda}
// 
// Input as contravariant (upper index) 4-vectors
//
//
Tensor1<std::complex<double>,4> MTensorPomeron::iG_yee(const M4Vec& prime, const M4Vec& p,
	const std::vector<std::complex<double>>& ubar, const std::vector<std::complex<double>>& u) const {

	//const double q2 = (prime-p).M2();
	const double e = msqrt(alpha_QED * 4.0 * PI); // ~ 0.3

	Tensor1<std::complex<double>, 4> T;
	for (const auto& mu : LI) {
		// \bar{spinor} [Gamma Matrix] \spinor product
		T(mu) = zi * e * gra::matoper::VecVecMultiply(ubar, gamma_lo[mu] * u);
	}
	return T;
}


// Gamma-Proton-Proton vertex: iGamma_\mu (p', p)
//
//
// contracted in \bar{spinor} iG_\mu \bar{spinor}
//
// Input as contravariant (upper index) 4-vectors
//
Tensor1<std::complex<double>,4> MTensorPomeron::iG_ypp(const M4Vec& prime, const M4Vec& p,
	const std::vector<std::complex<double>>& ubar, const std::vector<std::complex<double>>& u) const {

	const double t = (prime-p).M2();
	const double e = msqrt(alpha_QED * 4.0 * PI); // ~ 0.3
	const M4Vec psum = prime - p;

	Tensor1<std::complex<double>, 4> T;
	for (const auto& mu : LI) {

		MMatrix<std::complex<double>> SUM(4,4, 0.0);
		for (const auto& nu : LI) {
			SUM += sigma_lo[mu][nu] * psum[nu];
		}
		const MMatrix<std::complex<double>> A = gamma_lo[mu] * F1(t);
		const MMatrix<std::complex<double>> B = SUM * zi/(2*mp) * F2(t);
		const MMatrix<std::complex<double>> G = (A + B) * ((-1.0) * zi * e);

		// \bar{spinor} [Gamma Matrix] \spinor product
		T(mu) = gra::matoper::VecVecMultiply(ubar, G * u);
	}
	return T;
}


// High-Energy limit proton-Pomeron-proton vertex times \delta^{lambda_prime, \lambda} (helicity conservation)
// (~1.5 faster evaluation than the exact spinor structure)
Tensor2<std::complex<double>,4,4> MTensorPomeron::iG_PppHE(const M4Vec& prime, const M4Vec p) const {

	Tensor2<std::complex<double>,4,4> T;
	const M4Vec psum = prime + p;

	const std::complex<double> FACTOR = - zi * 3.0 * betaPNN * F1((prime-p).M2());

	for (auto& mu : LI) {
		for (auto& nu : LI) {
			T(mu,nu) = FACTOR * (psum % mu) * (psum % nu);
		}
	}
	return T;
}


// Pomeron-Proton-Proton [-Neutron-Neutron, -Antiproton-Antiproton] vertex function:
// contracted in \bar{spinor} G_{\mu\nu} \bar{spinor}
//
// i\Gamma_{\mu\nu} (p', p)
//
// Input as contravariant (upper index) 4-vectors
//
Tensor2<std::complex<double>,4,4> MTensorPomeron::iG_Ppp(const M4Vec& prime, const M4Vec& p, 
	const std::vector<std::complex<double>>& ubar, const std::vector<std::complex<double>>& u) const {

	const   double t = (prime-p).M2();
	const M4Vec psum = prime + p;

	Tensor2<std::complex<double>,4,4> T;
	const std::complex<double> FACTOR = - zi * 3.0 * betaPNN * F1(t);

	// Lorentz indices
	MMatrix<std::complex<double>> A;
	MMatrix<std::complex<double>> MAT;

	const MMatrix<std::complex<double>> slash = FSlash(psum);

	for (const auto& mu : LI) {
	for (const auto& nu : LI) {

		A = (gamma_lo[mu] * (psum%nu) + gamma_lo[nu] * (psum%mu)) * 0.5;
		if (mu == nu) { // Speed it up
			MAT = (A - slash * 0.25 * g[mu][nu]) * FACTOR;
		} else {
			MAT = A * FACTOR;
		}

		// \bar{spinor} [Gamma Matrix] \spinor product
		T(mu,nu) = gra::matoper::VecVecMultiply(ubar, MAT * u);
	}}
	return T;
}


// Double vertex function in High Energy Limit
// Computationally very heavy, should optimize this [TBD]
//
// Pomeron - outgoing proton spinor ubar -
//           <antiproton/proton fermion propagator> 
//         - outgoing antiproton spinor v - Pomeron
//
Tensor4<std::complex<double>,4,4,4,4> MTensorPomeron::iG_PppbarP(const M4Vec& prime, 
	const std::vector<std::complex<double>>& ubar, const M4Vec& pt, const MMatrix<std::complex<double>>& iSF,
	const std::vector<std::complex<double>>& v, const M4Vec& p) const {

	const std::complex<double> lhs_FACTOR = - zi * 3.0 * betaPNN * F1((prime - pt).M2());
	const std::complex<double> rhs_FACTOR = - zi * 3.0 * betaPNN * F1((pt - p).M2());

	const M4Vec lhs_psum = prime + pt;
	const M4Vec rhs_psum = pt + p;

	// Feynman slashes
	const MMatrix<std::complex<double>> lhs_slash = FSlash(lhs_psum);
	const MMatrix<std::complex<double>> rhs_slash = FSlash(rhs_psum);

	// Aux matrices
	MMatrix<std::complex<double>> lhs_A;
	MMatrix<std::complex<double>> lhs_MAT;

	MMatrix<std::complex<double>> rhs_A;
	MMatrix<std::complex<double>> rhs_MAT;

	Tensor4<std::complex<double>,4,4,4,4> T;

	for (const auto& mu2 : LI) {
	for (const auto& nu2 : LI) {

		lhs_A = (gamma_lo[mu2] * (lhs_psum%nu2) + gamma_lo[nu2] * (lhs_psum%mu2)) * 0.5;
		if (mu2 == nu2) { // Speed it up
			lhs_MAT = (lhs_A - lhs_slash * 0.25 * g[mu2][nu2]) * lhs_FACTOR;
		} else {
			lhs_MAT = lhs_A * lhs_FACTOR;
		}

		// Fermion propagator applied in the middle
		const std::vector<std::complex<double>> ubarM =
			gra::matoper::VecMatMultiply(ubar, lhs_MAT * iSF);

		for (const auto& mu1 : LI) {
		for (const auto& nu1 : LI) {

			rhs_A = (gamma_lo[mu1] * (rhs_psum%nu1) + gamma_lo[nu1] * (rhs_psum%mu1)) * 0.5;
			if (mu1 == nu1) { // Speed it up
				rhs_MAT = (rhs_A - rhs_slash * 0.25 * g[mu1][nu1]) * rhs_FACTOR;
			} else {
				rhs_MAT = rhs_A * rhs_FACTOR;
			}

			// \bar{spinor} [Gamma Matrix] \spinor product
			T(mu2,nu2,mu1,nu1) = gra::matoper::VecVecMultiply(ubarM, rhs_MAT * v);
		}}
	}}

	return T;
}

/*

[NOT IMPLEMENTED]

// Pomeron (\mu\nu) - Pomeron (\kappa\lambda) - f2 (\rho\sigma)
// vertex function: i\Gamma_{\mu\nu,\kappa\lambda,\rho\sigma}
//
// M0 (GeV) is the nominal mass
//
// Input as contravariant (upper index) 4-vector
//
MTensor<std::complex<double>> MTensorPomeron::iG_PPf2(const M4Vec& q1, const M4Vec& q2, double M_f2) const {

	MTensor<std::complex<double>> T = MTensor({4,4,4,4,4,4}, 0.0);
	auto FM_ = [] (double t) -> double {
		const double LAMBDA = 1.0; // GeV
		return 1.0 / (1.0 - t/pow2(LAMBDA));
	};
	auto F_ = [] (double q2, double M_f2) -> double {
		const double LAMBDA = 1.0; // GeV^2
		return std::exp(-pow2(q2 - pow2(M_f2)) / pow4(LAMBDA));
	};

	const double F_tilde = FM_(q1.M2()) * FM_(q2.M2()) * F_((q1+q2).M2(), M_f2);

	...
}

MTensor<std::complex<double>> MTensorPomeron::iG_PPf2_1(const M4Vec& q1, const M4Vec& q2) {

	MTensor<std::complex<double>> T = MTensor({4,4,4,4,4,4}, 0.0);

	const double M0 = 1.0;
	const double  g = 1.0; // Coupling

	for (auto const& u : LI) {
	for (auto const& v : LI) {
	for (auto const& k : LI) {
	for (auto const& l : LI) {
	for (auto const& r : LI) {
	for (auto const& s : LI) {

		for (auto const& v1 : LI) {
		for (auto const& a1 : LI) {
		for (auto const& l1 : LI) {
		for (auto const& r1 : LI) {
		for (auto const& s1 : LI) {
		for (auto const& u1 : LI) {

		if (n1 != a1 || l1 != r1 || s1 != u1) { continue; } // speed it up

		T({u,n,k,l,r,s}) = (2.0*zi) * g * M0 * R(u,v,u1,v1) * R(k,l,a1,l1) * R(r,s,r1,s1) * g[n1][a1] * g[l1][r1] * g[s1][u1];

		}} }} }}

	}}} }}}
}
*/


// Rho-Pi-Pi vertex function: i\Gamma_\mu(k1,k2)
//
// Input as contravariant (upper index) 4-vectors
//
Tensor1<std::complex<double>,4> MTensorPomeron::iG_rhopipi(const M4Vec& k1, const M4Vec& k2) const {

	const double g_rhopipi = 11.51;
	const M4Vec p = k1 - k2;

	Tensor1<std::complex<double>,4> T;
	for (const auto& mu : LI) {
		T(mu) = -0.5 * zi * g_rhopipi * (p % mu);
	}
	return T;
}


// f2-Pi-Pi vertex function: i\Gamma_{\kappa\lambda}(k1,k2)
//
// Input as contravariant (upper index) 4-vectors
//
Tensor2<std::complex<double>,4,4> MTensorPomeron::iG_f2pipi(const M4Vec& k1, const M4Vec& k2) const {

	// Form factor
	auto F_f2pipi = [&] (double q2) {
		return 1.0;
	};
	const double g_f2pipi = 9.26;    // Coupling
	const double M0 = 1.0;           // Mass scale (GeV)
	const double ksq = (k1+k2).M2();
	const M4Vec p = k1 - k2;

	Tensor2<std::complex<double>,4,4> T;
	for (const auto& k : LI) {
	for (const auto& l : LI) {
		T(k,l) = -zi * g_f2pipi / (2*M0) * ((p%k) * (p%l) - 0.25 * g[k][l] * p.M2()) * F_f2pipi(ksq);
	}}
	return T;
}


// f2-gamma-gamma vertex function: i\Gamma_{\mu\nu\kappa\lambda}(k1,k2)
//
// Input as contravariant (upper index) 4-vectors
//
Tensor4<std::complex<double>,4,4,4,4> MTensorPomeron::iG_f2yy(const M4Vec& k1, const M4Vec& k2) const {

	const double e = msqrt(alpha_QED * 4.0 * PI); // ~ 0.3

	// Form factor
	auto F_f2yy = [&] (double q2) {
		return 1.0;
	};

	// Couplings
	const double a_f2yy = pow2(e) / (4*PI) * 1.45; // GeV^{-3}
	const double b_f2yy = pow2(e) / (4*PI) * 2.49; // GeV^{-1}

	const Tensor4<std::complex<double>,4,4,4,4> G0 = Gamma0(k1,k2);
	const Tensor4<std::complex<double>,4,4,4,4> G2 = Gamma2(k1,k2);	

	const std::complex<double> FACTOR = zi * FM(k1.M2()) * FM(k2.M2()) * F_f2yy((k1+k2).M2());

	Tensor4<std::complex<double>,4,4,4,4> T;
	for (const auto& u : LI) {
	for (const auto& v : LI) {
	for (const auto& k : LI) {
	for (const auto& l : LI) {
		T(u,v,k,l) = FACTOR * (2.0 * a_f2yy * G0(u,v,k,l) - b_f2yy * G2(u,v,k,l)); 
	}}}}
	return T;
}


// Pomeron-Pseudoscalar-Pseudoscalar vertex function: i\Gamma_{\mu\nu}(p',p)
//
// Input as contravariant (upper index) 4-vectors
//
Tensor2<std::complex<double>,4,4> MTensorPomeron::iG_Ppsps(const M4Vec& prime, const M4Vec& p) const {

	double beta_P = 0.0;
	if (prime.M() < 0.2) { // Choose based on particle mass
		// Coupling (pion)		
		beta_P = 1.76; // GeV^{-1}
	} else {
		// Coupling (kaon)
		beta_P = 1.54; // GeV^{-1}
	}

	const M4Vec psum = prime + p;
	const double M2  = psum.M2();

	Tensor2<std::complex<double>,4,4> T;
	const std::complex<double> FACTOR = -zi * 2.0 * beta_P * FM((prime - p).M2());

	for (const auto& mu : LI) {
	for (const auto& nu : LI) {
		T(mu,nu) = FACTOR * ((psum%mu) * (psum%nu) - 0.25 * g[mu][nu] * M2);
	}}
	return T;
}


// Pomeron-Vector(massive)-Vector(massive) vertex function: i\Gamma_{\alpha\beta\gamma\delta} (p', p)
//
// Input M0 is the vector meson on-shell mass
//
// Input as contravariant (upper index) 4-vectors
//
Tensor4<std::complex<double>,4,4,4,4> MTensorPomeron::iG_Pvv(const M4Vec& prime, const M4Vec& p) const {

	// Couplings (phi meson)
	const double a = 0.49; // GeV^{-3}
	const double b = 4.27; // GeV^{-1}

	const Tensor4<std::complex<double>,4,4,4,4> G0 = Gamma0(prime, -p);
	const Tensor4<std::complex<double>,4,4,4,4> G2 = Gamma2(prime, -p);

	const std::complex<double> FACTOR = zi * FM((prime - p).M2());

	Tensor4<std::complex<double>,4,4,4,4> T;
	for (const auto& alpha : LI) {
	for (const auto& beta  : LI) {
	for (const auto& gamma : LI) {
	for (const auto& delta : LI) {

		T(alpha,beta,gamma,delta) = FACTOR * (2.0*a*G0(alpha,beta,gamma,delta) - b*G2(alpha,beta,gamma,delta));
	}}}}

	return T;
}


// Gamma-Vector meson vertex
//
// type = "rho", "omega", "phi"
//
Tensor2<std::complex<double>,4,4> MTensorPomeron::yV(std::string type) const {

	const double e = msqrt(alpha_QED * 4.0 * PI); // ~ 0.3
	double gammaV = 0.0;
	double mV = 0.0;

	if (type == "rho") {
		gammaV = msqrt( 4*PI / 0.496);
		mV = 0.770;
	}
	if (type == "omega") {
		gammaV = msqrt( 4*PI / 0.042);
		mV = 0.785;
	}
	if (type == "phi") {
		gammaV = (-1.0) * msqrt( 4*PI / 0.0716);
		mV = 1.020;
	}

	Tensor2<std::complex<double>,4,4> T;

	for (const auto& mu : LI) {
	for (const auto& nu : LI) {
		T(mu,nu) = -zi * e * pow2(mV) / gammaV * g[mu][nu];
	}}
	return T;
}


// ----------------------------------------------------------------------
// Propagators


// Tensor Pomeron propagator: (covariant == contravariant)
//
// i\Delta_{\mu\nu,\kappa\lambda}(s,t) = i\Delta^{\mu\nu,\kappa\lambda}(s,t)
//
// Input as (sub)-Mandelstam invariants: s,t
//
// 
// Symmetry relations:
// \Delta_{\mu\nu,\kappa\lambda} = \Delta_{\nu\mu,\kappa\lambda} = 
// \Delta_{\mu\nu,\lambda\kappa} = \Delta_{\kappa\lambda,\mu\nu}
//
// Contraction with g^{\mu\nu} or g^{\kappa\lambda} gives 0
//
Tensor4<std::complex<double>,4,4,4,4> MTensorPomeron::iD_P(double s, double t) const {
	
	Tensor4<std::complex<double>,4,4,4,4> T;

	for (const auto& u : LI) {
	for (const auto& v : LI) {
	for (const auto& k : LI) {
	for (const auto& l : LI) {
		T(u,v,k,l) =
			1.0 / (4.0*s) *
		    (g[u][k] * g[v][l] + g[u][l] * g[v][k] - 0.5 * g[u][v] * g[k][l]) *
			std::pow(-zi * s * ap_P, alpha_P(t) - 1.0);
	}}}}
	return T;
}


// Tensor Reggeon propagator iD (f2-a2-reggeon)
//
// Input as (sub)-Mandelstam invariants: s,t
//
Tensor4<std::complex<double>,4,4,4,4> MTensorPomeron::iD_2R(double s, double t) const {
	
	Tensor4<std::complex<double>,4,4,4,4> T;

	for (const auto& u : LI) {
	for (const auto& v : LI) {
	for (const auto& k : LI) {
	for (const auto& l : LI) {
		T(u,v,k,l) =
			1.0 / (4.0*s) *
		    (g[u][k] * g[v][l] + g[u][l] * g[v][k] - 0.5 * g[u][v] * g[k][l]) *
			std::pow(-zi * s * ap_2R, alpha_2R(t) - 1.0);
	}}}}
	return T;
}


// Vector Reggeon propagator iD (rho-reggeon, omega-reggeon)
//
// Input as (sub)-Mandelstam invariants s,t
//
Tensor2<std::complex<double>,4,4> MTensorPomeron::iD_1R(double s, double t) const {
	
	Tensor2<std::complex<double>,4,4> T;

	for (const auto& u : LI) {
	for (const auto& v : LI) {
		T(u,v) = zi * g[u][v] / pow2(M_1R) * std::pow(-zi * s * ap_1R, alpha_1R(t) - 1.0);
	}}
	return T;
}


// Odder propagator iD
//
// Input as (sub)-Mandelstam invariants s,t
//
Tensor2<std::complex<double>,4,4> MTensorPomeron::iD_O(double s, double t) const {
	
	Tensor2<std::complex<double>,4,4> T;

	for (const auto& u : LI) {
	for (const auto& v : LI) {
		T(u,v) = - zi * g[u][v] * eta_O / pow2(M_O) * std::pow(-zi * s * ap_O, alpha_O(t) - 1.0);
	}}
	return T;
}


// Scalar propagator iD with zero width
//
// Input as 4-vector
//
std::complex<double> MTensorPomeron::iD_MES0(const M4Vec& p, double M0) const {

	return zi / (p.M2() - pow2(M0));
}


// Scalar propagator iD with finite width
//
// Input as 4-vector
//
std::complex<double> MTensorPomeron::iD_MES(const M4Vec& p, double M0, double Gamma) const {

	const double m2 = p.M2();

	// Breit-Wigner
	const std::complex<double> Delta = 1.0 / (m2 - pow2(M0) + zi*M0*Gamma);
	return zi * Delta;
}


// Massive vector propagator iD_{\mu \nu}
//
// Input as contravariant 4-vector
//
Tensor2<std::complex<double>,4,4> MTensorPomeron::iD_VMES(const M4Vec& p, double M0, double Gamma) const {
	
	Tensor2<std::complex<double>,4,4> T;

	// Transverse part
	const double m2 = p.M2();
	const std::complex<double> Delta_T = 1.0 / (m2 - pow2(M0) + zi*M0*Gamma);

	// Longitudinal part (does not enter here)
	//const double Delta_L = 0.0;

	for (const auto& mu : LI) {
	for (const auto& nu : LI) {
		T(mu,nu) = zi *( -g[mu][nu] + (p%mu)*(p%nu)/m2 ) * Delta_T;
		// - zi*(p%u)*(p%v)/m2 * Delta_L; // this part does not enter
	}}
	return T;
}


// Tensor meson propagator iD_{\mu \nu \kappa \lambda}
//
// Input as contravariant (upper index) 4-vector
//
Tensor4<std::complex<double>,4,4,4,4> MTensorPomeron::iD_TMES(const M4Vec& p, double M0, double Gamma) const {

	const double m2 = p.M2();
	Tensor4<std::complex<double>, 4,4,4,4> T;

	// Inline lambda function
	auto ghat = [&] (int mu, int nu) {
		return -g[mu][nu] + (p%mu)*(p%nu)/m2;
	};

	// Breit-Wigner
	const std::complex<double> Delta = 1.0 / (m2 - pow2(M0) + zi*M0*Gamma);

	for (const auto& u : LI) {
	for (const auto& v : LI) {
	for (const auto& k : LI) {
	for (const auto& l : LI) {
		T(u,v,k,l) = zi * Delta *
			(0.5 * (ghat(u,k)*ghat(v,l) + ghat(u,l)*ghat(v,k)) - 1.0/3.0 *  ghat(u,v)*ghat(k,l));
	}}}}
	return T;
}


// ----------------------------------------------------------------------
// Trajectories (simple linear/affine here)

// Pomeron trajectory
double MTensorPomeron::alpha_P(double t) const {
	return 1.0 + delta_P + ap_P*t;
}
// Odderon trajectory
double MTensorPomeron::alpha_O(double t) const {
	return 1.0 + delta_O + ap_O*t;
}
// Reggeon (rho,omega) trajectory
double MTensorPomeron::alpha_1R(double t) const {
	return 1.0 + delta_1R + ap_1R*t;
}
// Reggeon (f2,a2) trajectory
double MTensorPomeron::alpha_2R(double t) const {
	return 1.0 + delta_2R + ap_2R*t;
}


// ----------------------------------------------------------------------
// Form factors

// Proton electromagnetic Dirac form factor (electric)
double MTensorPomeron::F1(double t) const {
	return (1 - (t/(4*mp*mp)) * mu_ratio) / (1 - t/(4*mp*mp)) * GD(t);
}

// Proton electromagnetic Pauli form factor (magnetic)
double MTensorPomeron::F2(double t) const {
	return (mu_ratio - 1)  / (1 - t/(4*mp*mp)) * GD(t); 
}

// Dipole form factor
double MTensorPomeron::GD(double t) const {
	const double m2D = 0.71; // GeV^2
	return 1.0/pow2(1 - t/m2D);
}

// Meson form factor (pion electromagnetic type)
double MTensorPomeron::FM(double q2) const {
	const double A2 = 0.5; // GeV^2
	return 1.0 / (1.0 - q2/A2);
}


// ----------------------------------------------------------------------
// Tensor functions
 
// Output: \Gamma^{(0)}_{\mu\nu\kappa\lambda}
//
// Input must be contravariant (upper index) 4-vectors
//
Tensor4<std::complex<double>,4,4,4,4> MTensorPomeron::Gamma0(const M4Vec& k1, const M4Vec& k2) const {

	const double k1k2 = k1 * k2; // 4-dot product
	Tensor4<std::complex<double>,4,4,4,4> T;

	for (const auto& u : LI) {
	for (const auto& v : LI) {
	for (const auto& k : LI) {
	for (const auto& l : LI) {
		T(u,v,k,l) = (k1k2 * g[u][v] -  (k2%u) * (k1%v)) *
				     ((k1%k) * (k2%l) + (k2%k) * (k1%l) - 0.5 * k1k2 * g[k][l]);
	}}}}
	return T;
}


// Output: \Gamma^{(2)}_{\mu\nu\kappa\lambda}
//
// Input must be contravariant (upper index) 4-vectors
//
Tensor4<std::complex<double>,4,4,4,4> MTensorPomeron::Gamma2(const M4Vec& k1, const M4Vec& k2) const {

	const double k1k2 = k1 * k2; // 4-dot product
	Tensor4<std::complex<double>,4,4,4,4> T;

	for (const auto& u : LI) {
	for (const auto& v : LI) {
	for (const auto& k : LI) {
	for (const auto& l : LI) {
		T(u,v,k,l) = 
			k1k2 * (g[u][k]*g[v][l] + g[u][l]*g[v][k]) + g[u][v]*((k1%k) * (k2%l) + (k2%k) * (k1%l))
		 - (k1%v) * (k2%l)*g[u][k] - (k1%v) * (k2%k) *g[u][l] - (k2%u) * (k1%l) *g[v][k] - (k2%u) * (k1%k)*g[v][l]
		 - (k1k2 * g[u][v] - (k2%u)*(k1%v)) * g[k][l]; 
	}}}}
	return T;
}

// Output: Auxialary (helper) tensor
// 1/2 g_{\mu\kappa} g_{\nu\lambda} + 1/2 g_{\mu\lambda}g_{\nu\kappa} - 1/4 g_{\mu\nu}g_{kappa\lambda}
//
void MTensorPomeron::CalcRTensor() {

	Tensor4<std::complex<double>,4,4,4,4> T;

	for (const auto& mu     : LI) {
	for (const auto& nu     : LI) {
	for (const auto& kappa  : LI) {
	for (const auto& lambda : LI) {

		T(mu,nu,kappa,lambda) =  0.5 * g[mu][kappa]*g[nu][lambda] +
								 0.5 * g[mu][lambda]*g[nu][kappa] +
								0.25 * g[mu][nu]*g[kappa][lambda];
	}}}}

	TR = T;
}



} // gra namespace ends
