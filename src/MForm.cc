// Form factors, structure functions, Regge trajectories etc. parametrizations
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <complex>
#include <cstdlib>
#include <iostream>
#include <map>
#include <vector>

// Own
#include "Graniitti/MForm.h"
#include "Graniitti/MKinematics.h"
#include "Graniitti/MMath.h"
#include "Graniitti/MProcess.h"
#include "Graniitti/MSpin.h"

// Libraries
#include "json.hpp"
#include "rang.hpp"

using gra::math::msqrt;
using gra::math::pow2;
using gra::math::pow3;
using gra::math::zi;
using gra::math::PI;
using gra::math::abs2;

using gra::PDG::mp;
using gra::PDG::mpi;


namespace gra {


// Model parameters
namespace PARAM_SOFT {

	// Pomeron trajectory
	double DELTA_P = 0.0;
	double ALPHA_P = 0.0;

	// Couplings
	double gN_P = 0.0;
	double g3P = 0.0;
	double gamma = 0.0;

	// Proton form factor
	double fc1 = 0.0;
	double fc2 = 0.0;

	// Pion loop
	double fc3 = 0.0;

	std::string GetHashString() {
		std::string str = std::to_string(PARAM_SOFT::DELTA_P) +
		                  std::to_string(PARAM_SOFT::ALPHA_P) +
		                  std::to_string(PARAM_SOFT::gN_P) +
		                  std::to_string(PARAM_SOFT::fc1) +
		                  std::to_string(PARAM_SOFT::fc2) +
		                  std::to_string(PARAM_SOFT::fc3);
		return str;
	}

	void PrintParam() {
		printf("PARAM_SOFT:: Soft model parameters: \n\n");
		printf("- DELTA_P = %0.4f \n", DELTA_P);
		printf("- ALPHA_P = %0.4f [GeV^{-2}] \n", ALPHA_P);
		printf("- gN_P    = %0.4f [GeV^{-1}] \n", gN_P);
		printf("- g3P     = %0.4f \n", g3P / gN_P); // Convention
		printf("- gamma   = %0.4f \n", gamma);
		printf("- fc1     = %0.4f [GeV^2] \n", fc1);
		printf("- fc2     = %0.4f [GeV^2] \n", fc2);
		printf("- fc3     = %0.4f [GeV^2] \n", fc3);
		std::cout << std::endl << std::endl;
	}
}

// Flat amplitude parameters
namespace PARAM_FLAT {

double b = 0.0;
}

// Monopole parameters
namespace PARAM_MONOPOLE {

int En = 0;                 // Bound state energy level
double M0 = 0.0;            // Monopole mass
double Gamma0 = 0.0;        // Monopolium width
std::string coupling = "";  // Coupling scenarios
int gn = 0;                 // Dirac charge 1,2,3,...

bool initialized = false;

// Functions as a solution to radial Schroedinger equation with Coulomb type
// potential
// V(r) ~= -g^2/(4pi) * 1/r
//
// Monopolium binding energy eigenvalues
// -> Monopolium system sits at lower energy than monopole + antimonopole

// Monopolium running width
double GammaMP(double n, double alpha_g) {
	return (8.0 * PI * pow2(alpha_g)) / pow2(PARAM_MONOPOLE::M0) * math::abs2(PsiMP(n));
}

// Binding energy
double EnergyMP(double n) {
	return -pow2(1 / (8.0 * form::alpha_EM(0))) * PARAM_MONOPOLE::M0 / (n*n);
	//return -2.0*M0/15.0; // DEBUG
}

// Monopolium wavefunction in the origin of the bound system
double PsiMP(double n) {
	return (1.0 / msqrt(PI)) * std::pow(PARAM_MONOPOLE::M0 / (8.0 * form::alpha_EM(0) * n), 3.0/2.0);
}

void PrintParam(double sqrts) {
	std::cout << rang::style::bold
	          << "Monopolium process parameters:" << rang::style::reset
	          << std::endl
	          << std::endl;

	printf("- M0     = %0.3f GeV \n", PARAM_MONOPOLE::M0);
	printf("- Gamma0 = %0.3f GeV \n", PARAM_MONOPOLE::Gamma0);
	printf("- En     = %d \n",        PARAM_MONOPOLE::En);
	
	const double M = 2*PARAM_MONOPOLE::M0 + EnergyMP(PARAM_MONOPOLE::En);
	printf("\nGives monopolium mass: M = %0.3f GeV (Binding Energy = %0.3f GeV) \n\n",
		M, M - 2*PARAM_MONOPOLE::M0);
	
	// Check we have enough energy
	if (M > sqrts) {
		std::string str =
		    "gra::form::PARAM_MONOPOLE::PrintParam: FATAL error "
		    "Monopolium mass > CMS energy!";
		throw std::invalid_argument(str);
	}

	std::cout << "Dirac charge = "    << PARAM_MONOPOLE::gn << std::endl;
	std::cout << "Coupling scheme = " << PARAM_MONOPOLE::coupling << std::endl;
	
}
}


// N*, N**, N** excitation couplings
namespace PARAM_NSTAR {
std::vector<double> rc = {0.0, 0.0, 0.0};
}

namespace form {

// LHAPDFset name
std::string LHAPDF;

// Read resonance parameters
gra::PARAM_RES ReadResonance(const std::string& resparam_str, MRandom& rng) {

	std::cout << "gra::form::ReadResonance: Reading " + resparam_str + " ";

	// Create a JSON object from file
	std::string inputfile = gra::aux::GetBasePath(2) + "/modeldata/" + resparam_str;

	// Read and parse
	std::string data;
	nlohmann::json j;

	try {
		data = gra::aux::GetInputData(inputfile);
		j    = nlohmann::json::parse(data);
	} catch (...) {
		throw std::invalid_argument("form::ReadResonance: Error parsing '" + resparam_str + "'");
	}

	// Resonance parameters
	gra::PARAM_RES res;

	try {

	std::complex<double> g(0, 0);
	double g_A   = j.at("PARAM_RES").at("g_A");
	double g_phi = j.at("PARAM_RES").at("g_phi");
	res.g = g_A * std::exp(std::complex<double>(0,1) * g_phi); // Complex coupling
	
	// PDG code
	res.p.pdg = j.at("PARAM_RES").at("PDG");

	// "Glueball state"
	res.p.glue = j.at("PARAM_RES").at("glue");
	
	// mass
	res.p.mass = j.at("PARAM_RES").at("M");
	if (res.p.mass < 0) {
		std::string str = "MAux:ReadResonance:: <" + resparam_str + "> Invalid M < 0 !";
		throw std::invalid_argument(str);
	}

	// width
	res.p.width = j.at("PARAM_RES").at("W");
	if (res.p.width < 0) {
		std::string str = "MAux:ReadResonance:: <" + resparam_str + "> Invalid W < 0 !";
		throw std::invalid_argument(str);
	}

	const double J = j.at("PARAM_RES").at("J");
	res.p.spinX2 = J * 2;
	if (res.p.spinX2 < 0) {
		std::string str = "MAux:ReadResonance:: <" + resparam_str + "> Invalid J < 0 !";
		throw std::invalid_argument(str);
	}

	//
	res.p.P = j.at("PARAM_RES").at("P");
	if (!(res.p.P == -1 || res.p.P == 1)) {
		std::string str = "MAux:ReadResonance:: <" + resparam_str + "> Invalid P (not -1 or 1) !";
		throw std::invalid_argument(str);
	}

	// Validity of these is taken care of in the functions
	res.BW             = j.at("PARAM_RES").at("BW");

	// If we have spin
	if (res.p.spinX2 != 0) {

		// Spin dependent
		res.hc.FRAME               = j.at("PARAM_RES").at("FRAME");
		const bool P_conservation  = true;
		
		const int n = res.p.spinX2 + 1; // n = 2J + 1
		MMatrix<std::complex<double>> rho(n,n);

		// Draw random density matrices (until the number set by user)
		if (j.at("PARAM_RES").at("random_rho") > 0) {
			for (std::size_t k = 0; k < j.at("PARAM_RES").at("random_rho"); ++k) {
				rho = gra::spin::RandomRho(res.p.spinX2/2.0, P_conservation, rng);
			}

		// Construct spin density matrix from the input
		} else {
			for (std::size_t a = 0; a < rho.size_row(); ++a) {
				for (std::size_t b = 0; b < rho.size_col(); ++b) {
					const double Re = j.at("PARAM_RES").at("rho_real").at(a).at(b);
					const double Im = j.at("PARAM_RES").at("rho_imag").at(a).at(b);
					rho[a][b] = Re + zi * Im;
				}
			}
			// Check positivity conditions
			if (gra::spin::Positivity(rho, res.p.spinX2/2.0) == false) {
				std::string str = "gra::form::ReadResonance: <" + resparam_str +
				                  "> Input density matrix not positive definite!";
				throw std::invalid_argument(str);
			}
		}
		res.hc.rho = rho;
	}
	std::cout << " [DONE]" << std::endl;

	} catch (nlohmann::json::exception& e) {
		throw std::invalid_argument("form::ReadResonance: Missing parameter in '"
			+ resparam_str + "' : " + e.what());
	}

	return res;
}


// Regge signature factor, alpha_t is alpha(t), and signature sigma = +-1
//
// Remember, there are poles (also) at t < 0 (scattering domain) with linear trajectory
//
// sin(pi x) = pi / (-x Gamma(x) Gamma(-x))
//
std::complex<double> ReggeEta(double alpha_t, double sigma) {
	const double denom = std::sin(PI * alpha_t);
	return -(1.0 + sigma * std::exp(-zi * PI * alpha_t)) / denom;
}


// ----------------------------------------------------------------------
// Non-linear Pomeron trajectory functional form:
//
// [REFERENCE: Khoze, Martin, Ryskin, https://arxiv.org/abs/hep-ph/0007359]
double S3PomAlpha(double t) {
	// Additive quark model Beta_pi / Beta_p = 2/3 ->
	const double BETApi = (2.0 / 3.0) * PARAM_SOFT::gN_P;
	const double PIC = (pow2(BETApi) * pow2(mpi)) / (32.0 * pow3(PI));

	// Pomeron trajectory value
	const double h =
	    PIC * S3HPL(4.0 * pow2(mpi) / std::abs(t), t); // pion loop insert (ADD with minus sign)
	const double alpha_P =
	    1.0 + PARAM_SOFT::DELTA_P + PARAM_SOFT::ALPHA_P * t - h;

	return alpha_P;
}

// Pion loop insert to the non-linear Pomeron trajectory
double S3HPL(double tau, double t) {
	const double m = 1.0; // fixed scale (GeV)
	const double sqrtau = msqrt(1 + tau);

	// Pion-Pomeron form factor parametrization
	const double F_pi = 1.0 / (1.0 - t / PARAM_SOFT::fc3);

	return (4.0 / tau) * pow2(F_pi) *
	       (2.0 * tau -
	        std::pow(1.0 + tau, 3.0 / 2.0) *
	            std::log((sqrtau + 1.0) / (sqrtau - 1.0)) +
	        std::log((m * m) / (mpi * mpi)));
}
// ----------------------------------------------------------------------


// Elastic proton form factor parametrization
//
// <apply at amplitude level>
// [REFERENCE: Khoze, Martin, Ryskin, https://arxiv.org/abs/hep-ph/0007359]
double S3F(double t) {
	return (1.0 / (1 - t / PARAM_SOFT::fc1)) *
	       (1.0 / (1 - t / PARAM_SOFT::fc2));
}


// Proton inelastic form factor / structure function
// parametrization for Pomeron processes [THIS FUNCTION IS ANSATZ - IMPROVE!]
// 
// <apply at amplitude level>
double S3FINEL(double t, double M2) {

	constexpr double DELTA_P = 0.0808;
	constexpr       double a = 0.5616; // GeV^{2}

	double f = std::pow(std::abs(t) / (M2*(std::abs(t) + a)), 0.5 * (1 + DELTA_P));

	// Coupling ansatz
	f *= msqrt(PARAM_SOFT::g3P / PARAM_SOFT::gN_P);

	return f;
}


// Proton inelastic structure function F2(x,Q^2) parametrization 
// 
// The basic idea is that at low-Q^2, a fully non-perturbative description (parametrization) is needed.
// At high Q^2, DGLAP evolution could be done in log(Q^2) starting from the input description.
// 
// Now, some (very) classic ones have been implemented. 
// TBD: Can we "invert" F_2 and F_L e.g from LUXqed pdfs, or interface to other library.
// 
// [REFERENCE: Donnachie, Landshoff, https://arxiv.org/abs/hep-ph/9305319]
// [REFERENCE: Capella, Kaidalov, Merino, Tran Tranh Van, https://arxiv.org/abs/hep-ph/9405338v1]
//
double F2xQ2(double xbj, double Q2) {

	const std::string F2TYPE = "CKMT";

	if (F2TYPE == "DL") {

		constexpr double A = 0.324;
		constexpr double B = 0.098;

		constexpr double DELTA_P = 0.0808;
		constexpr double DELTA_R = 0.5475;

		constexpr double a = 0.561991692786383;
		constexpr double b = 0.011133;
		
		const double F2 = 
		  A * std::pow(xbj, - DELTA_P) * std::pow(Q2 / (Q2 + a), 1 + DELTA_P)
		+ B * std::pow(xbj, 1-DELTA_R) * std::pow(Q2 / (Q2 + b), DELTA_R);

		return F2;
	}

	else if (F2TYPE == "CKMT") {
		constexpr double A        = 0.1502;
		constexpr double B_u      = 1.2064;
		constexpr double B_d      = 0.1798;
		constexpr double alpha_R  = 0.4150;
		constexpr double DELTA_0  = 0.0800;

		constexpr double a        = 0.2631;
		constexpr double b        = 0.6452;
		constexpr double c        = 3.5489;
		constexpr double d        = 1.1170;

		const double n_Q2     = (3.0/2.0) * (1 + Q2/(Q2 + c));
		const double DELTA_Q2 = DELTA_0 * (1 + (2*Q2) / (Q2 + d));
		
		const double C1 = std::pow(Q2/(Q2 + a), 1.0 + DELTA_Q2);
		const double C2 = std::pow(Q2/(Q2 + b), alpha_R);

		const double F2 = A * std::pow(xbj, -DELTA_Q2)    * std::pow(1-xbj, n_Q2 + 4.0)  * C1 + 
			   				  std::pow(xbj, 1.0-alpha_R)  *
			   		   (B_u * std::pow(1-xbj, n_Q2) + B_d * std::pow(1-xbj, n_Q2 + 1.0)) * C2;

	   	return F2;
    } else {
    	throw std::invalid_argument("gra::form::F2xQ2: Unknown F2TYPE = " + F2TYPE);
    }
}


// "Purely magnetic structure function"
//
// Callan-Gross relation for spin-1/2: F_2(x) = 2xF_1(x) under Bjorken scaling
// For spin-0, F_1(x) = 0
//
double F1xQ2(double xbj, double Q2) {

	// F_L(xbj,Q2) = (1 + 4*pow2(xbj*mp)/Q2) * F2(xbj, Q2) - 2xbj * F1(xbj,Q2)

	return F2xQ2(xbj, Q2) / (2.0*xbj);
}


// ============================================================================
// Photon flux densities and form factors, input Q^2 as positive

// Alpha EM(Q^2 = 0)
double alpha_EM(double Q2) { // no running here
	return 1.0 / 137.035999139;
}
// Electric charge in natural units ~ 0.3
double e_EM(double Q2) {  // no running here
	return msqrt(alpha_EM(Q2) * 4.0 * PI);
}
double e_EM() {
	return msqrt(alpha_EM(0.0) * 4.0 * PI);
}

// kT unintegrated coherent EPA photon flux as in:
// 
// [REFERENCE: Harland-Lang, Khoze, Ryskin, https://arxiv.org/abs/1601.03772]
// [REFERENCE: Luszczak, Schaefer, Szczurek, https://arxiv.org/abs/1802.03244]
// 
// Form factors:
// [REFERENCE, Punjabi et al., https://arxiv.org/abs/1503.01452v4]
// 
// Proton electromagnetic form factors: Basic notions, present
// achievements and future perspectives, Physics Reports, 2015
// <https://www.sciencedirect.com/science/article/pzi/S0370157314003184>
// 
// [REFERENCE: Budnev, Ginzburg, Meledin, Serbo, EPA paper, Physics Reports, 1976]
// <https://www.sciencedirect.com/science/article/pzi/0370157375900095>
// 
// 
// Proton electric form factor (F_electric == F_2 == Pauli)
double F_E(double Q2) {
	return (4.0 * pow2(mp) * pow2(G_E(Q2)) + Q2 * pow2(G_M(Q2))) /
	       (4.0 * pow2(mp) + Q2);
}
// Proton magnetic form factor (F_magnetic == F_1  == Dirac)
double F_M(double Q2) {
	return pow2(G_M(Q2));
}
// Rosenbluth separation:
// low-Q^2 dominated by G_E, high-Q^2 dominated by G_M

// "Sachs" form factors:
// <http://www.scholarpedia.org/article/Nucleon_Form_factors>
constexpr double mup = 2.792847337; // Proton magnetic moment

/*
// The simplest possible: Dipole parametrization

double G_E(double Q2) {
  return G_M(Q2) / mup; // Scaling assumption
}
const double lambda2 = 0.71;
double G_M(double Q2) {
  return mup / pow2(1.0 + Q2/lambda2);
}
*/

// "Sachs Form Factor" goes as follows:
// G_E(0) = 1 for proton, 0 for neutron
// G_M(0) = mu_p for proton, mu_n for neutrons

// Simple parametrization of nucleon EM-form factors
// 	
// [REFERENCE: JJ Kelly, https://journals.aps.org/prc/pdf/10.1103/PhysRevC.70.068202]
double G_E(double Q2) {
	static const std::vector<double> a = {1, -0.24};
	static const std::vector<double> b = {10.98, 12.82, 21.97};

	const double tau = Q2 / (4.0 * mp * mp);

	// Numerator
	double num = 0.0;  // 0
	num += a[0];       // a_0 tau^0
	num += a[1] * tau; // a_1 tau^1

	// Denominator
	double den = 1.0;        // 1.0
	den += b[0] * tau;       // b_1 tau^1
	den += b[1] * pow2(tau); // b_2 tau^2
	den += b[2] * pow3(tau); // b_3 tau^3

	return num / den;
}

double G_M(double Q2) {
	static const std::vector<double> a = {1, 0.12};
	static const std::vector<double> b = {10.97, 18.86, 6.55};

	const double tau = Q2 / (4.0 * mp * mp);

	// Numerator
	double num = 0.0;  // 0
	num += a[0];       // a_0 tau^0
	num += a[1] * tau; // a_1 tau^1

	// Denominator
	double den = 1.0;        // 1.0
	den += b[0] * tau;       // b_1 tau^1
	den += b[1] * pow2(tau); // b_2 tau^2
	den += b[2] * pow3(tau); // b_3 tau^3

	return mup * num / den;
}

// Coherent photon flux from proton
// xi ~ longitudinal momentum loss [0,1]
// t  ~ Mandelstam t
// pt ~ proton transverse momentum
// 
// 
//  p ---------F_E---------> p' with xi = 1 - p^*_z'/p_z
//               $
//                $
//                 $
//
// Factors applied here:
//  
//  1/xi    [~ sub Moller flux]
//  1/pt2   [~ kt-factorization] (cancels with pt2 from numerator)
//  16pi^2  [~ kinematics volume factor]
//
double CohFlux(double xi, double t, double pt) {

	const double pt2 = pow2(pt);
	const double xi2 = pow2(xi);
	const double mp2 = pow2(mp);
	const double Q2  = std::abs(t);

	double f = alpha_EM(0) / PI  * 
	           (pt2 / (pt2 + xi2 * mp2)) *
	           ((1.0 - xi) * (pt2 / (pt2 + xi2 * mp2)) * F_E(Q2) +
	           	(xi2 / 4.0) * F_M(Q2));
    
	// Factors
	f /= xi;
	f /= pt2;
	f *= 16.0 * gra::math::PIPI;

	return f; // Use at cross section level
}


// Incoherent photon flux from a dissociated proton with mass M.
// When M -> mp, this reproduces CohFlux() if
// F2(x,Q^2) does reproduce the elastic limit (not all parametrizations do)
// 
// 
//  p ---------F2(x,Q^2)----->  p* with xi = 1 - p^*_z / p_z
//             -------------->
//             ---x---------->
//                 $
//                  $
//                   $
//
// Factors applied as with CohFlux() above.
//
double IncohFlux(double xi, double t, double pt, double M2) {
	
	constexpr double mp2 = pow2(mp);

	const double pt2 = pow2(pt);
	const double xi2 = pow2(xi);
	const double Q2  = std::abs(t);
	const double xbj = Q2 / (Q2 + M2 - mp2); // Bjorken-x

	double f = alpha_EM(0) / PI *
			   (pt2 / (pt2 + xi*(M2 - mp2) + xi2*mp2)) *
	           ((1.0 - xi) * (pt2 / (pt2 + xi*(M2 - mp2) + xi2*mp2)) * F2xQ2(xbj,Q2) / (Q2 + M2 - mp2) +
	           	(xi2 / (4.0*pow2(xbj))) * 2.0 * xbj * F1xQ2(xbj,Q2) / (Q2 + M2 - mp2) );
    
	// Factors
	f /= xi;
	f /= pt2;
	f *= 16.0 * gra::math::PIPI;

	return f; // Use at cross section level
}


// Drees-Zeppenfeld proton coherent gamma flux (collinear)
double DZFlux(double x) {
	const double Q2min = (pow2(mp) * pow2(x)) / (1.0 - x);
	const double A = 1.0 + 0.71 / Q2min;

	double f = alpha_EM(0) / (2.0 * PI * x) * (1.0 + pow2(1.0 - x)) *
	           (std::log(A) - 11.0 / 6.0 + 3.0 / A -
	            3.0 / (2.0 * pow2(A)) + 1.0 / (3 * pow3(A)));

	return f; // Use at cross section level
}


// Breit-Wigner propagators / form factors
//
//
//
// Useful identity for normalization:
//
// \int dm^2 \frac{1}{(m^2 - M0^2)^2 + M0^2Gamma^2} \equiv \frac{\pi}{M0 Gamma}
//
// based on squaring the complex propagator
// D(m^2) = 1 / (m^2 - M0^2 + iM0*Gamma), and integrating.
//
std::complex<double> CBW(const gra::LORENTZSCALAR& lts,
                         const gra::PARAM_RES& resonance) {
	switch (resonance.BW) {
		case 1:
			return CBW_FW(lts.m2, resonance.p.mass,
			              resonance.p.width);
		case 2:
			return CBW_RW(lts.m2, resonance.p.mass,
			              resonance.p.width);
		case 3:
			return CBW_BF(lts.m2, resonance.p.mass,
			              resonance.p.width, resonance.p.spinX2/2.0,
			              lts.decaytree[0].p4.M(),
			              lts.decaytree[1].p4.M());
		case 4:
			return CBW_JR(lts.m2, resonance.p.mass,
						  resonance.p.width, resonance.p.spinX2/2.0);

		default:
			throw std::invalid_argument(
			    "CBW: Unknown BW (Breit-Wigner) parameter: " +
			    std::to_string(resonance.BW));
	}
}

// See e.g.
// [REFERENCE: TASI Lectures on propagators, http://users.ictp.it/~smr2244/tait-supplemental.pdf]
// [REFERENCE: Cacciapaglia, Deandrea, Curtis, https://arxiv.org/abs/0906.3417v2]
// [REFERENCE: http://www.t2.ucsd.edu/twiki2/pub/UCSDTier2/Physics214Spring2015/ajw-breit-wigner-cbx99-55.pdf]

// ----------------------------------------------------------------------
// Delta function \delta(\hat{s} - M_0^2) replacement function:
//    \int d\hat{s} \delta(\hat{s} - M_0^2) -> \int d\hat{s} deltaBW(\hat{s},M0,Gamma)
// 
// To be applied at cross section level
double deltaBWxsec(double shat, double M0, double Gamma) {
	return M0 * Gamma / PI / (pow2(shat - pow2(M0)) + pow2(M0 * Gamma));
}

// To be applied at amplitude level
double deltaBWamp(double shat, double M0, double Gamma) {
	return msqrt(deltaBWxsec(shat, M0, Gamma));
}
// ----------------------------------------------------------------------


// ----------------------------------------------------------------------
// Breit-Wigner propagator parametrizations

// 1. Complex Fixed Width Relativistic Breit-Wigner
std::complex<double> CBW_FW(double m2, double M0, double Gamma) {
	return -1.0 / (m2 - M0 * M0 + zi * M0 * Gamma);
}

// 2. Complex Running Width Relativistic Breit-Wigner
std::complex<double> CBW_RW(double m2, double M0, double Gamma) {
	return -1.0 / (m2 - M0 * M0 + zi * std::sqrt(m2) * Gamma);
}

// 3. J = 0,1,2 Complex Relativistic Breit-Wigner
// with angular barrier effects type of Bleit-Weiskopf
// mA and mB are the masses of daughters (GeV)
std::complex<double> CBW_BF(double m2, double M0, double Gamma, int J,
                            double mA, double mB) {
	
	const double u = msqrt((m2 - pow2(mA + mB)) *
						   (m2 - pow2(mA - mB))) / (2 * msqrt(m2));
	const double d = msqrt((M0 * M0 - pow2(mA + mB)) *
						   (M0 * M0 - pow2(mA - mB))) / (2 * M0);
	const double Bfactor = pow(u / d, 2 * J + 1);
	return -1.0 / (m2 - M0 * M0 + zi * Gamma * M0 * M0 / msqrt(m2) * Bfactor);
}

// 4. Spin dependent Relativistic Breit-Wigner
//
// [REFERENCE: https://arxiv.org/pdf/1402.1178.pdf]
std::complex<double> CBW_JR(double m2, double M0, double Gamma, double J) {

	const std::complex<double> denom = (m2 - M0*M0 + zi*M0*Gamma);

	if      (static_cast<int>(J) == 0)   { // J = 0
		return -1.0 / denom;
	}
	else if (static_cast<int>(J*2) == 1) { // J = 1/2
		return -2.0*msqrt(m2) / denom;
	}
	else if (static_cast<int>(J) == 1)   { // J = 1
		return -(1.0 - m2/(M0*M0)) / denom;
	}
	else if (static_cast<int>(J*2) == 3) { // J = 3/2
		return -(2.0/3.0)*msqrt(m2)*(1.0-m2/(M0*M0)) / denom;
	}
	else if (static_cast<int>(J) == 2)   { // J = 2
		return -(7.0/6.0 - (4.0/3.0) * (m2/(M0*M0)) +
			(2.0/3.0) * (m2*m2)/(gra::math::pow4(M0))) / denom;
	}
	else {
		return -1.0 / denom; // Too high spin
	}
}


} // form namespace ends
} // gra namespace ends
