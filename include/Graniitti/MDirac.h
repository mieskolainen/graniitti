// Dirac gamma algebra routines
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MDIRAC_H
#define MDIRAC_H

// C++
#include <complex>
#include <random>
#include <vector>

// Own
#include "Graniitti/MKinematics.h"
#include "Graniitti/M4Vec.h"
#include "Graniitti/MMatrix.h"
#include "Graniitti/MMath.h"

// FTensor
#include "FTensor.hpp"


namespace gra {


class MDirac {

public:

	MDirac();
	~MDirac() {}

	void InitGammaMatrices(const std::string& basis);
	bool GammaAntiCommutation() const;

	// Polarization vectors/tensors
	FTensor::Tensor1<std::complex<double>,4>          EpsSpin1(const M4Vec& k, int m) const;
	FTensor::Tensor1<std::complex<double>,4>   EpsMassiveSpin1(const M4Vec& k, int m) const;
	FTensor::Tensor2<std::complex<double>,4,4> EpsMassiveSpin2(const M4Vec& k, int m) const;
	
	// Helicity spinors
	std::vector<std::complex<double>> XiSpinor(const M4Vec& p, int helicity) const;
	std::vector<std::complex<double>> uHelChiral(const M4Vec& p, int helicity) const;
	std::vector<std::complex<double>> vHelChiral(const M4Vec& p, int helicity) const;
	
	std::vector<std::complex<double>> uHelDirac(const M4Vec& p, int helicity) const;
	std::vector<std::complex<double>> vHelDirac(const M4Vec& p, int helicity) const;

	// Dirac adjoint
	std::vector<std::complex<double>> Bar(const std::vector<std::complex<double>>& spinor) const;
	
	// Feynman slash matrix operator
	MMatrix<std::complex<double>>     FSlash(const M4Vec& a) const;
	
	// Propagators
	FTensor::Tensor2<std::complex<double>,4,4> iD_y(const M4Vec& q) const;
	MMatrix<std::complex<double>>     iD_F(const M4Vec& q, double m) const;

	// Dirac spinors
	std::vector<std::complex<double>> uDirac(const M4Vec& p, int spin) const;
	std::vector<std::complex<double>> vDirac(const M4Vec& p, int spin) const;
	
	// Spinor-Helicity style methods
	std::complex<double>              sProd(const M4Vec& p1, const M4Vec& p2, int helicity) const;
	std::vector<std::complex<double>> uGauge(const M4Vec& p, int helicity) const;
	std::vector<std::complex<double>> vGauge(const M4Vec& p, int helicity) const;
	
	// Charge conjugate operator
	MMatrix<std::complex<double>> C_up() const;

	// Right and left chiral projectors
	MMatrix<std::complex<double>> PR() const;
	MMatrix<std::complex<double>> PL() const;

	// Angular Momentum operators
	MMatrix<std::complex<double>> J_operator(unsigned int i) const;

    bool   SpinorHELimit(const M4Vec& p1, const M4Vec& p2) const;
	bool   DiracSpinorComplete(const M4Vec& p, const std::string& type, const std::string& basis) const;
	bool   MassiveSpin1Complete(const M4Vec& k) const;
	double FSlashFSlash(const M4Vec& p) const;
	
	// Kronecker delta	
	double Delta(int i, int j) const {
		return (i == j) ? 1.0 : 0.0; 
	}
	
	// Transform matrix from the Weyl to Dirac basis by an operator sandwich:
	//
	// \gamma^\mu_Dirac = S \gamma_{Weyl}^\mu S^{-1}
	//
	// S = S^{-1} <=> So this matrix works also from Dirac to Weyl basis
	// 
	// Spinor Transformation:
	// 
	// u_Dirac = S u_Weyl
	// v_Dirac = S v_Weyl
	//
	// S = 1/\sqrt{2}(1 + y^5y^0)
	//
	MMatrix<std::complex<double>> S_basis =
			{std::vector<std::complex<double>>{1.0/std::sqrt(2.0), 0.0,  1.0/std::sqrt(2.0),  0.0},
		 	 std::vector<std::complex<double>>{0.0, 1.0/std::sqrt(2.0),  0.0,  1.0/std::sqrt(2.0)},
		 	 std::vector<std::complex<double>>{1.0/std::sqrt(2.0), 0.0, -1.0/std::sqrt(2.0),  0.0},
		 	 std::vector<std::complex<double>>{0.0, 1.0/std::sqrt(2.0),  0.0, -1.0/std::sqrt(2.0)}};

	// Pauli matrices
	MMatrix<std::complex<double>> sigma_x =
		{std::vector<std::complex<double>> {0.0, 1.0},
		 std::vector<std::complex<double>> {1.0, 0.0}};

	MMatrix<std::complex<double>> sigma_y =
		{std::vector<std::complex<double>>{0.0, -gra::math::zi},
		 std::vector<std::complex<double>>{gra::math::zi,  0.0}};

	MMatrix<std::complex<double>> sigma_z =
		{std::vector<std::complex<double>>{1.0, 0.0},
		 std::vector<std::complex<double>>{0.0, 1.0}};


	// Contravariant (up) and covariant (lo) gamma matrix set
	std::vector<MMatrix<std::complex<double>> > gamma_up;
	std::vector<MMatrix<std::complex<double>> > gamma_lo;

	// Contravariant (up) and covariant (lo) \sigma_{\mu\nu} matrices
	std::vector<std::vector<MMatrix<std::complex<double>>>> sigma_up;
	std::vector<std::vector<MMatrix<std::complex<double>>>> sigma_lo;

	// Indices etc.
	std::vector<std::size_t> LI = {0,1,2,3};
	std::vector<int> SPINORSTATE = {-1,1};

	// Minkowski metric tensor (+,-,-,-) in convention (t,px,py,pz)
	MMatrix<double> g = MMatrix<double>(4,4,"minkowski");
	
	// Identity matrix
	MMatrix<std::complex<double>> I4 = MMatrix<std::complex<double>>(4,4,"eye");

protected:

	std::string BASIS = ""; // D for Dirac, C for Chiral
	
};


} // gra namespace ends


#endif
