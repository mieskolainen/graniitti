// 'Durham QCD' Processes and Amplitudes
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.


// C++
#include <complex>
#include <random>
#include <vector>

// Own
#include "Graniitti/MAux.h"
#include "Graniitti/MForm.h"
#include "Graniitti/MKinematics.h"
#include "Graniitti/MSudakov.h"
#include "Graniitti/MDurham.h"
#include "Graniitti/MPDG.h"

// Libraries
#include "json.hpp"
#include "rang.hpp"

using gra::math::msqrt;
using gra::math::pow2;
using gra::math::zi;
using gra::math::PI;
using gra::PDG::mp;


namespace gra {


// 2-body incoming or outgoing helicity combinations
// (keep it algebraic order to match with MadGraph)
const int MM = 0; // --
const int MP = 1; // -+
const int PM = 2; // +-
const int PP = 3; // ++

// 2-body final state helicity configurations
const std::vector<int> fs2 = {MM,MP,PM,PP};

// Central system J^P = 0^+, 0^-, +2^+, -2^+
enum SPINPARITY { P0, M0, P2, M2 }; // Implicit conversion to int


// Durham loop integral discretization technical parameters
namespace Durham {

	unsigned int N_qt = 0;      // Number of qt discretization intervals
	unsigned int N_phi = 0;     // Number of phi discretization intervals
	
	double qt2_MIN = 0; // Loop momentum qt^2 minimum (GeV^2)
	double qt2_MAX = 0; // Loop momentum qt^2 maximum (GeV^2)
	
	std::string PDF_scale = "MIN"; // Scheme
	double alphas_scale = 4.0;

	void ReadParameters() {

		using json = nlohmann::json;

		const std::string inputfile = gra::aux::GetBasePath(2) + "/modeldata/" + gra::aux::MODELPARAM + ".json";
		const std::string data = gra::aux::GetInputData(inputfile);
		json j;
		
		try {
			j = json::parse(data);
			
			// JSON block identifier
			const std::string XID = "PARAM_DURHAM_QCD";
			N_qt        = j[XID]["N_qt"];
			N_phi       = j[XID]["N_phi"];
			qt2_MIN     = j[XID]["qt2_MIN"];
			qt2_MAX     = j[XID]["qt2_MAX"];

			PDF_scale   = j[XID]["PDF_scale"];
			alphas_scale = j[XID]["alphas_scale"];

			// Now calculate rest
			qt_MIN = msqrt(Durham::qt2_MIN);
			qt_MAX = msqrt(Durham::qt2_MAX);

			qt_STEP = (Durham::qt_MIN - Durham::qt_MAX) / Durham::N_qt;
			phi_STEP = (2.0 * gra::math::PI) / Durham::N_phi;

			initialized = true;

		} catch (...) {
			std::string str =
			    "Durham::ReadParameters: Error parsing " +
			    inputfile + " (Check for extra/missing commas)";
			throw std::invalid_argument(str);
		}
	}
	
	// THESE ARE CALCULATED FROM ABOVE
	double qt_MIN = 0;
	double qt_MAX = 0;
	double qt_STEP = 0;
	double phi_STEP = 0;
	bool initialized = false;
}

// Durham QCD / KMR model
//
// [REFERENCE: Pumplin, Phys.Rev.D 52 (1995)]
// [REFERENCE: Khoze, Kaidalov, Martin, Ryskin, Stirling, https://arxiv.org/abs/hep-ph/0507040]
// [REFERENCE: Khoze, Martin, Ryskin, https://arxiv.org/abs/hep-ph/0605113]
// [REFERENCE: Harland-Lang, Khoze, Ryskin, Stirling, https://arxiv.org/abs/1005.0695]
// [REFERENCE: Harland-Lang, Khoze, Ryskin, https://arxiv.org/abs/1409.4785]
//
std::complex<double> MDurham::DurhamQCD(gra::LORENTZSCALAR& lts, const std::string& process) {
	
	// First run, init parameters
	gra::aux::g_mutex.lock();
	if (gra::aux::sudakov.initialized == false) {
		try {
			gra::aux::sudakov.Init(lts.sqrt_s, gra::form::LHAPDF, true);
		} catch (std::invalid_argument e) {
			gra::aux::g_mutex.unlock(); // need to release here,
			                        // otherwise get infinite lock
			throw(e);
		}
	}
	if (!Durham::initialized) {
		Durham::ReadParameters();
	}
	gra::aux::g_mutex.unlock();

	if (process == "gg") {
		std::vector<std::vector<std::complex<double>>> Amp(
		    4, std::vector<std::complex<double>>(4, 0.0));

		// Madgraph
		gra::aux::g_mutex.lock();
		const double alpha_s = gra::aux::sudakov.AlphaS_Q2(lts.s_hat);
		gra::aux::g_mutex.unlock();
		AmpMG5_gg_gg.CalcAmp(lts, alpha_s);

		// Amplitude evaluated outside the Qt-loop (approximation)
		Dgg2gg(lts, Amp);

		// Run loop
		return DQtloop(lts, Amp);
	} else if (process == "qqbar") {

		std::vector<std::vector<std::complex<double>>> Amp(
		    4, std::vector<std::complex<double>>(4, 0.0));
		
		// Madgraph
		gra::aux::g_mutex.lock();
		const double alpha_s = gra::aux::sudakov.AlphaS_Q2(lts.s_hat);
		gra::aux::g_mutex.unlock();
		AmpMG5_gg_qqbar.CalcAmp(lts, alpha_s);

		// Amplitude evaluated outside the Qt-loop (approximation)
		Dgg2qqbar(lts, Amp);

		// Run loop
		return DQtloop(lts, Amp);
	} else if (process == "chi0") {
		std::vector<std::vector<std::complex<double>>> Amp(
		    1, std::vector<std::complex<double>>(4, 0.0));

		// Run loop
		std::complex<double> A = DQtloop(lts, Amp);
		return A;

	// ------------------------------------------------------------
	//
	// Implement more processes here ...
	//
	// ------------------------------------------------------------

	} else {
		std::string str =
		    "MDurham::Durham: Unknown subprocess: " + process;
		throw std::invalid_argument(str);
	}
}


// Helicity basis decomposition
//
// In the forward limit q1_t = - q2_t = Q_t
// => gives gluon polarization vectors eps_1 = -eps_2 => central system J_z = 0
//
inline void MDurham::DHelicity(const std::vector<double>& q1,
                          	   const std::vector<double>& q2,
                          	   std::vector<std::complex<double>>& JzP) const {
	
	const unsigned int X = 0; // component for readability
	const unsigned int Y = 1;

	// 1/2  q1_t dot q2_t
	JzP[P0] = - 0.5 * (q1[X] * q2[X] + q1[Y] * q2[Y]);

	// 1/2 i |q1_t x q2_t|
	JzP[M0] = - 0.5 * gra::math::zi * std::abs(q1[X] * q2[Y] - q1[Y] * q2[X]);

	const double Re = 0.5 * (q1[X] * q2[X] - q1[Y] * q2[Y]);
	const double Im = 0.5 * (q1[X] * q2[Y] + q1[Y] * q2[X]);

	// +1/2[ (xx - yy) + i(xy + yx) ]
	JzP[P2] = Re + gra::math::zi * Im;

	// +1/2[ (xx - yy) - i(xy + yx) ]
	JzP[M2] = Re - gra::math::zi * Im;
}


// M_{++,--,+-,-+} denotes M_{\lambda_1,\lambda_2}
// i.e. g(\lambda_1) g(\lambda_2) -> X (system) helicity amplitudes
//
// where lambda1, lambda2 are the initial state gluon helicities
//
// [REFERENCE: https://arxiv.org/abs/1405.0018v2, formula (20)]
//
inline std::complex<double> MDurham::DHelProj(
    const std::vector<std::complex<double>>& A,
    const std::vector<std::complex<double>>& JzP) const {
	
	// M_{++} + M_{--}
	// (this term gives J_z^PC = 0^++ selection rule in the forward pt->0 limit)
	const std::complex<double> aP0 = JzP[P0] * (A[PP] + A[MM]);
	// M_{++} - M_{--}
	const std::complex<double> aM0 = JzP[M0] * (A[PP] - A[MM]);
	// M_{-+}
	const std::complex<double> aP2 = JzP[P2] * A[MP];
	// M_{+-}
	const std::complex<double> aM2 = JzP[M2] * A[PM];

	/*
	// DEBUG
	std::cout << "0+ : " << aP0 << std::endl;
	std::cout << "0- : " << aM0 << std::endl;
	std::cout << "2+ : " << aP2 << std::endl;
	std::cout << "2- : " << aM2 << std::endl;
	std::cout << std::endl << std::endl;
	*/

	return aP0 + aM0 + aP2 + aM2;
}

/*
// MADGRAPH HELICITY FORMAT [direct algebraic order]
const static int helicities[ncomb][nexternal] =
  {{-1, -1, -1, -1}, 0
   {-1, -1, -1,  1}, 1
   {-1, -1,  1, -1}, 2
   {-1, -1,  1,  1}, 3

   {-1,  1, -1, -1}, 4
   {-1,  1, -1,  1}, 5
   {-1,  1,  1, -1}, 6 
   {-1,  1,  1,  1}, 7

   { 1, -1, -1, -1}, 8
   { 1, -1, -1,  1}, 9
   { 1, -1,  1, -1}, 10
   { 1, -1,  1,  1}, 11

   { 1,  1, -1, -1}, 12
   { 1,  1, -1,  1}, 13
   { 1,  1,  1, -1}, 14
   { 1,  1,  1,  1}}; 15
*/


// Alternative (semi-ad-hoc) scenarios for the scale choise
inline void MDurham::DScaleChoise(double qt2, double q1_2, double q2_2, double& Q1_2_scale, double& Q2_2_scale) const {

	if      (Durham::PDF_scale == "MIN") {
		Q1_2_scale = std::min(qt2, q1_2);
		Q2_2_scale = std::min(qt2, q2_2);
	}
	else if (Durham::PDF_scale == "MAX") {
		Q1_2_scale = std::max(qt2, q1_2);
		Q2_2_scale = std::max(qt2, q2_2);
	}
	else if (Durham::PDF_scale == "IN") {
		Q1_2_scale = q1_2;
		Q2_2_scale = q2_2;
	}
	else if (Durham::PDF_scale == "EX") {
		Q1_2_scale = qt2;
		Q2_2_scale = qt2;
	}
	else if (Durham::PDF_scale == "AVG") {
		Q1_2_scale = (qt2 + q1_2) / 2.0;
		Q2_2_scale = (qt2 + q2_2) / 2.0;
	} else {
		throw std::invalid_argument("MDurham::DScaleChoise: Unknown 'Durham::PDF_scale' option!");
	}
}



// [REFERENCE: Khoze, Martin, Ryskin, https://journals.aps.org/prd/pdf/10.1103/PhysRevD.56.5867]
// [REFERENCE: Harland-Lang, Khoze, Martin, Ryskin, Stirling, https://arxiv.org/abs/1405.0018v2]
//
// See also: <https://arxiv.org/pdf/1608.03765.pdf> (LUND)
//
//
// Durham loop integral amplitude:
// A = pi^2 \int \frac{d^2 Q_t M(gg->X)}{Q_t^2(Q_t-{p_1t})^2(Q_t+p_{2t})^2}
//               x f_g(x_1,x_1',Q_1^2,\mu_2;t_1) x
//               f_g(x_2,x_2',Q_2^2,\mu_2;t_2)
//
std::complex<double> MDurham::DQtloop(
    gra::LORENTZSCALAR& lts, std::vector<std::vector<std::complex<double>>> Amp) {
	
	// Forward (proton) system pt-vectors
	const std::vector<double> pt1 = {lts.pfinal[1].Px(), lts.pfinal[1].Py()};
	const std::vector<double> pt2 = {lts.pfinal[2].Px(), lts.pfinal[2].Py()};

	// *************************************************************************
	// ** Process scale (GeV) / Sudakov suppression kt^2 integral upper bound **
	const double MU = msqrt(lts.s_hat / Durham::alphas_scale);
	// *************************************************************************

	// Init 2D-Simpson weight matrix (will be calculated only once, being
	// static), C++11 handles multithreaded static initialization
	const static MMatrix<double> WSimpson =
	    gra::math::Simpson38Weight2D(Durham::N_qt, Durham::N_phi);

	// NOTE N + 1, init with zero!
	std::vector<MMatrix<std::complex<double>>> f(Amp.size(),
	    MMatrix<std::complex<double>>(Durham::N_qt + 1, Durham::N_phi + 1, 0.0));

	// Spin-Parity
	std::vector<std::complex<double>> JzP(4, 0.0);

	// 2D-loop integral
	//
	// \int d^2 \vec{qt} [...] = \int dphi \int dqt qt [...]
	//

	// Linearly discretized qt-loop, N+1!
	for (std::size_t i = 0; i < Durham::N_qt+1; ++i) {
		const double qt = Durham::qt_MIN + i * Durham::qt_STEP;
		const double qt2 = pow2(qt);

		// Linearly discretized phi in [0,2pi), N+1!
		for (std::size_t j = 0; j < Durham::N_phi+1; ++j) {
			const double qphi = j * Durham::phi_STEP;

			// --------------------------------------------------------------------------
			// Loop vector
			const std::vector<double> qt_ = {qt * std::cos(qphi),
			                                 qt * std::sin(qphi)};

			// Fusing gluon pt-vectors
			const std::vector<double> q1 = {qt_[0] - pt1[0],
			                                qt_[1] - pt1[1]};
			const std::vector<double> q2 = {qt_[0] + pt2[0],
			                                qt_[1] + pt2[1]};

			const double q1_2 = gra::math::vpow2(q1);
			const double q2_2 = gra::math::vpow2(q2);

			// Get fusing gluon spin-parity (J_z^P) components
			// [q1,q2] -> [0^+,0^-,+2^+,-2^-]
			DHelicity(q1, q2, JzP);

			// ** Durham scale choise **
			double Q1_2_scale = 0.0;
			double Q2_2_scale = 0.0;
			DScaleChoise(qt2, q1_2, q2_2, Q1_2_scale, Q2_2_scale);

			// Minimum scale cutoff
			if (Q1_2_scale < Durham::qt2_MIN || Q2_2_scale < Durham::qt2_MIN) { continue; }

			// --------------------------------------------------------------------------
			// Get amplitude level pdfs
			const double fg_1 = gra::aux::sudakov.fg_xQ2M(lts.x1, Q1_2_scale, MU);
			const double fg_2 = gra::aux::sudakov.fg_xQ2M(lts.x2, Q2_2_scale, MU);
			
			// Amplitude weight:
			// * \pi^2 : see original KMR papers
			// *    2  : factor of from initial state
			// boson-statistics, check it:!
			// *    qt : jacobian of d^2qt -> dphi dqt qt
			std::complex<double> weight = fg_1 * fg_2 / (qt2 * q1_2 * q2_2);
			weight *= gra::math::PIPI * 2.0 * qt;

			// Update sub-amplitude
			if (Amp.size() == 1) { // Calculate here
				Dgg2chi0(lts, Amp, q1, q2);
			}
			
			// Loop over final state helicity combinations
			for (std::size_t h = 0; h < f.size(); ++h) {
			  f[h][i][j] = weight * DHelProj(Amp[h], JzP);
			}

		} // phi-loop
	} // |q_t|-loop
	
	// Evaluate the total numerical integral for each helicity amplitude
	std::vector<std::complex<double>> sum(f.size(), 0.0);
	for (std::size_t h = 0; h < f.size(); ++h) {
		sum[h] = gra::math::Simpson38Integral2D(f[h], WSimpson, Durham::qt_STEP, Durham::phi_STEP);
	}

	// *****Make sure it is of right size!*****
	lts.hamp.resize(f.size());

	// Final outcome as a coherent sum over incoming gluons
	// That is, sum at amplitude level -
	// in Contrast to the usual \sum_h |A_h|^2
	std::complex<double> A(0, 0);
	for (std::size_t h = 0; h < f.size(); ++h) {

		lts.hamp[h] = sum[h];

		// Apply proton form factors
		lts.hamp[h] *= lts.excite1 ? gra::form::S3FINEL(lts.t1) : gra::form::S3F(lts.t1);
		lts.hamp[h] *= lts.excite2 ? gra::form::S3FINEL(lts.t2) : gra::form::S3F(lts.t2);

		// Apply phase space factors
		lts.hamp[h] *= msqrt(16.0 * gra::math::PIPI);
		lts.hamp[h] *= msqrt(16.0 * gra::math::PIPI);
		lts.hamp[h] *= msqrt(lts.s / lts.s_hat);

		// For the total complex amplitude (for no eikonal screening applied)
		A += lts.hamp[h];
	}
	return A;
}


// ======================================================================
// Exclusive \chi_c(0+) sub-amplitude
//
// [REFERENCE: Khoze, Martin, Ryskin, Stirling, https://arxiv.org/abs/hep-ph/0403218]
// [REFERENCE: Pasechnik, Szczurek, Teryaev, https://arxiv.org/abs/0709.0857v1]
//
void MDurham::Dgg2chi0(const gra::LORENTZSCALAR& lts, std::vector<std::vector<std::complex<double>>>& Amp,
                       const std::vector<double>& qt1, const std::vector<double>& qt2) const {
	
	// @@ MUTEX LOCK @@
	gra::aux::g_mutex.lock();
	const double alpha_s = gra::aux::sudakov.AlphaS_Q2(lts.s_hat / Durham::alphas_scale);
	gra::aux::g_mutex.unlock();
	// @@ MUTEX LOCK @@

	const double gs2 = 4.0 * PI * alpha_s; // coupling
	const double K_NLO = 1.68;             // NLO correction
	const double M0 = 3.41475;             // chi_c(0+) mass (GeV)
	const double W0 = 0.0108;              // chi_c(0+) width (GeV)
	const double NC = 3.0;                 // #colors

	// Gluonic width \Gamma(\chi_c(0+) -> gg), see references
	std::complex<double> A = K_NLO * 8.0 * gra::math::zi * gs2 / M0 * msqrt(0.075) / msqrt(PI * M0 * NC);

	// Apply Breit-Wigner propagator shape delta-function
	A *= gra::form::deltaBWamp(lts.s_hat, M0, W0);
	
	// Initial state gluon combinations
	std::vector<std::complex<double>> is(4, 0.0);
	is[MM] = A;
	is[MP] = A;
	is[PM] = A;
	is[PP] = A;

	// 0+ final state
	Amp[P0] = is;
}


// ======================================================================
// gg -> gg tree-level helicity amplitudes
// 
//
// Basic result: d\hat{\sigma}/dt = 9/4 \pi \alpha_s^2 / E_T^4
//
// In pure gluon amplitudes: when gluon helicities are the same,
// or at most one is different from the rest, vanish for any n >= 4,
// where n is the total number of gluons (in+out)
//
// [REFERENCE: Dixon, https://arxiv.org/abs/1310.5353v1]
// 
// Note that incoming gluons (a,b) are in a color singlet state,
// i.e., a = b or \delta^{ab} Tr [ .color algebra. ] has been applied
// 
void MDurham::Dgg2gg(const gra::LORENTZSCALAR& lts,
                       std::vector<std::vector<std::complex<double>>>& Amp) {
	
	// @@ MUTEX LOCK @@
	gra::aux::g_mutex.lock();
	const double alpha_s = gra::aux::sudakov.AlphaS_Q2(lts.s_hat / Durham::alphas_scale);
	gra::aux::g_mutex.unlock();
	// @@ MUTEX LOCK @@

	// Vertex factor coupling, gs^2 = 4 pi alpha_s
	double norm = 4.0 * PI * alpha_s;
	
	// Color part
	const double NC = 3.0;
	norm *= NC / msqrt(NC*NC - 1);
	norm *= 2.0;

	// Helicity phase
	const double phi = lts.decaytree[0].p4.Phi();
	const std::complex<double> negphase = std::exp(-zi * phi);
	const std::complex<double> posphase = std::exp( zi * phi);

	const double aux01 = pow2(lts.s_hat) / (lts.u_hat * lts.t_hat);
	const double aux23 = lts.u_hat / lts.t_hat;

	// We use natural binary order / Madgraph order here
	std::vector<std::complex<double>> dualAmp(16, 0.0);

	// Loop over final state gluon pair helicity combinations
	for (std::size_t h = 0; h < fs2.size(); ++h) {
		
		std::vector<std::complex<double>> is(4, 0.0);
		if (fs2[h] == PP || fs2[h] == MM) {

			is[MM] = (fs2[h] == PP) ? 0 : norm*aux01;                         // -- => ++ or --
			is[MP] = 0;                                                       // -+ => ++ or --
			is[PM] = 0;                                                       // +- => ++ or --
			is[PP] = (fs2[h] == PP) ? norm*aux01 : 0;                         // ++ => ++ or --

		} else if (fs2[h] == PM || fs2[h] == MP) {

			is[MM] = 0;                                                       // -- => +- or -+
			is[MP] = negphase * norm * ((fs2[h] == PM) ? 1/aux23 :   aux23);  // -+ => +- or -+
			is[PM] = posphase * norm * ((fs2[h] == PM) ?   aux23 : 1/aux23);  // +- => +- or -+
			is[PP] = 0;                                                       // ++ => +- or -+			
		}
		Amp[h] = is;

		for (std::size_t k = 0; k < 4; ++k) {
			dualAmp[4*k+h] = is[k];
		}
	}

	// ------------------------------------------------------------------------
	const bool DEBUG = false;
	if (DEBUG) {
	// TEST amplitude squared WITH MADGRAPH (UNDER IMPLEMENTATION)
	double  madsum = 0.0;
	double thissum = 0.0;
	for (unsigned int i = 0; i < 16; ++i) {
		printf("i=%2d : |mad|^2 = %0.3f , |amp|^2 = %0.3f \n", i, gra::math::abs2(lts.hamp[i]), gra::math::abs2(dualAmp[i]));
		madsum  += gra::math::abs2(lts.hamp[i]);
		thissum += gra::math::abs2(dualAmp[i]); 
	}
	Asum += thissum / madsum;
	Nsum += 1.0;
	printf("THIS/MADGRAPH = %0.10f \n\n", Asum / Nsum);
	}
	// ------------------------------------------------------------------------
}


// ======================================================================
// gg -> qqbar tree-level helicity amplitudes
//
//
void MDurham::Dgg2qqbar(const gra::LORENTZSCALAR& lts,
						std::vector<std::vector<std::complex<double>>>& Amp) {

	throw std::invalid_argument("MDurham::DurhamQCD: qqbar amplitude in the next version");
}

} // gra namespace ends
