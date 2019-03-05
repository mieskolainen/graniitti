// Factorized type phase space class
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.


// C++
#include <algorithm>
#include <complex>
#include <iostream>
#include <random>
#include <vector>

// Own
#include "Graniitti/MAux.h"
#include "Graniitti/MForm.h"
#include "Graniitti/MKinematics.h"
#include "Graniitti/MMath.h"
#include "Graniitti/MMatrix.h"
#include "Graniitti/MProcess.h"
#include "Graniitti/MFactorized.h"
#include "Graniitti/MSpin.h"
#include "Graniitti/MUserCuts.h"
#include "Graniitti/MFragment.h"

// Libraries
#include "rang.hpp"

using gra::aux::indices;

using gra::math::msqrt;
using gra::math::pow2;
using gra::math::pow3;
using gra::math::pow4;
using gra::math::pow5;
using gra::math::CheckEMC;
using gra::math::zi;
using gra::math::PI;
using gra::math::abs2;

using gra::PDG::GeV2barn;


namespace gra {


double auxw = 0.0;

// This is needed by construction
MFactorized::MFactorized() {
	std::vector<std::string> supported = {"PP","yP","yy","gg"};
	ProcPtr = MSubProc(supported);
	ConstructProcesses();
}

// Constructor
MFactorized::MFactorized(std::string process) {

	InitHistograms();
	SetProcess(process);

	// Init final states
	M4Vec zerovec(0, 0, 0, 0);
	for (std::size_t i = 0; i < 10; ++i) {
		lts.pfinal.push_back(zerovec);
	}
	std::cout << "MFactorized:: [Constructor done]" << std::endl;
}

// Destructor
MFactorized::~MFactorized() {
}

// This is manually constructed and updated here
void MFactorized::ConstructProcesses() {

	Processes.clear();
	const std::string CID = "F";
	for (auto const& x : ProcPtr.descriptions) {
	    std::map<std::string, std::string> value = x.second;
	    for (auto const& y : value) {
	    	Processes.insert(std::make_pair(x.first + "[" + y.first + "]<"+CID+">", y.second));
		}
	}
}


// Initialize cut and process spesific postsetup
void MFactorized::post_Constructor() {

	if (ProcPtr.CHANNEL == "RES") {
		// Here we support only single resonances
		if (lts.RESONANCES.size() != 1) {
			std::string str = "MFactorized::post_Constructor: Only single resonance supported for this process (RESPARAM.size() != 1)";
			throw std::invalid_argument(str);
		}
	}

	// Set sampling boundaries
	ProcPtr.SetTechnicalBoundaries(gcuts, EXCITATION);

	// Initialize phase space dimension
	ProcPtr.LIPSDIM = 5 + 1; // All processes, +1 from central system mass
}


// Update kinematics (screening kT loop calls this)
// Exact 4-momentum conservation at loop vertices, that is, by using this one
// does not assume vanishing external momenta in the screening loop calculation.

bool MFactorized::LoopKinematics(const std::vector<double>& p1p,
                                   const std::vector<double>& p2p) {

	const uint Nf = lts.decaytree.size() + 2; // Number of final states

	// SET new final states pT degrees of freedom
	lts.pfinal[1].SetPxPy(p1p[0], p1p[1]);
	lts.pfinal[2].SetPxPy(p2p[0], p2p[1]);

	const double pt1 = lts.pfinal[1].Pt();
	const double pt2 = lts.pfinal[2].Pt();

	// Read original variables
	const double m1 = lts.pfinal_orig[1].M();
	const double m2 = lts.pfinal_orig[2].M();

	const double mX = lts.pfinal_orig[0].M();
	const double yX = lts.pfinal_orig[0].Rap();
	
	// Central system momentum
	lts.pfinal[0].SetPxPy(-(lts.pfinal[1].Px() + lts.pfinal[2].Px()),
	                      -(lts.pfinal[1].Py() + lts.pfinal[2].Py()));
	
	// Central system pz and E
	const double mtX = msqrt(pow2(mX) + lts.pfinal[0].Pt2());
	lts.pfinal[0].SetPzE(mtX * std::sinh(yX), mtX * std::cosh(yX));
	
	// Energy overflow
	if (lts.pfinal[0].E() > (lts.sqrt_s - (m1 + m2))) { return false; }

	double p1z = gra::kinematics::SolvepPZ1(m1, m2, pt1, pt2, lts.pfinal[0].Pz(), lts.pfinal[0].E(), lts.sqrt_s);
	double p2z = -(lts.pfinal[0].Pz() + p1z); // by momentum conservation
	if (std::isnan(p1z)) { return false; }
	
	// *** Enforce scattering direction +p -> +p, -p -> -p ***
	p1z = std::abs(p1z);
	p2z = -std::abs(p2z);

	// Pz and E of protons/N*
	lts.pfinal[1].SetPzE(p1z, msqrt(pow2(m1) + pow2(pt1) + pow2(p1z)) );
	lts.pfinal[2].SetPzE(p2z, msqrt(pow2(m2) + pow2(pt2) + pow2(p2z)) );

	if (!CheckEMC((lts.pbeam1 + lts.pbeam2) - (lts.pfinal[1] + lts.pfinal[2] + lts.pfinal[0]))) { return false; }

	return GetLorentzScalars(Nf);
}


// Return Monte Carlo integrand weight
double MFactorized::EventWeight(const std::vector<double>& randvec, AuxIntData& aux) {

	double W = 0.0;

	// Kinematics and cuts
	aux.kinematics_ok = B51RandomKin(randvec);
	aux.fidcuts_ok    = FiducialCuts();
	aux.vetocuts_ok   = VetoCuts();
	
	if (aux.Valid()) {
		
		// Matrix element squared
		const double MatESQ = (FLATAMPLITUDE == 0) ? abs2(S3ScreenedAmp()) : GetFlatAmp2(lts);
		
		// Calculate central system Phase Space volume
		double exact = 0.0;
		DecayWidthPS(exact);
		lts.DW_sum_exact.Add(gra::kinematics::MCW(exact, pow2(exact), 1), aux.vegasweight);
		
		// Add to the integral sum and take into account the VEGAS weight
		lts.DW_sum.Add(lts.DW, aux.vegasweight);

		double C_space = 1.0;
		if (lts.decaytree.size() != 0 && !ISOLATE) { // We have some legs in the central system
			C_space = (lts.DW.Integral()/(2*PI));
		}
		
		// ** EVENT WEIGHT **
		// /(2*PI) from phase space factorization
		W = C_space * (1.0 / S_factor) * B51PhaseSpaceWeight() * B51IntegralVolume() * MatESQ *
		    GeV2barn / MollerFlux(); // Total weight: phase-space x |M|^2 x barn units
	}
	
	aux.amplitude_ok = CheckInfNan(W);
	
	// As the last STEP: Histograms
	if (!aux.burn_in_mode) {
		const double totalweight = W * aux.vegasweight;
		FillHistograms(totalweight, lts);
	}

	return W;
}

// Fiducial cuts
bool MFactorized::FiducialCuts() const {

	return CommonCuts();
}

// Record event
bool MFactorized::EventRecord(HepMC::GenEvent& evt) {
	
	return CommonRecord(evt);
}


void MFactorized::PrintInit(bool silent) const {
	if (!silent) {
		PrintSetup();

		// Construct prettyprint diagram
		std::string proton1 = "-----------EL--------->";
		std::string proton2 = "-----------EL--------->";

		if (EXCITATION == 1) {
			proton1 = "-----------F2-xxxxxxxx>";
		}
		if (EXCITATION == 2) {
			proton1 = "-----------F2-xxxxxxxx>";
			proton2 = "-----------F2-xxxxxxxx>";
		}

		std::vector<std::string> feynmangraph;
		feynmangraph = {"||          ",
		                "||          ",
		                "xx--------->",
		                "||          ",
		                "||          "};

		// Print diagram
		std::cout << proton1 << std::endl;
		for (const auto& i : indices(feynmangraph)) {
			if (SCREENING) { // Put red
				std::cout << rang::fg::red << "     **    " << rang::style::reset;
			} else {
				std::cout << rang::fg::red << "           " << rang::style::reset;
			}
			std::cout << feynmangraph[i] << std::endl;
		}
		std::cout << proton2 << std::endl;
		std::cout << std::endl;
		
		// --------------------------------------------------------------
		// Amplitude setup printout
		// Monopolium
		if (ProcPtr.CHANNEL == "monopolium(0)") {
			PARAM_MONOPOLE::PrintParam(lts.sqrt_s);
		}

		// Custom resonances processes
		if (ProcPtr.CHANNEL == "RES") {
			for (auto& x : lts.RESONANCES) {
				x.second.PrintParam(lts.sqrt_s);
			}
		}
		// --------------------------------------------------------------

		std::cout << std::endl;
		std::cout << rang::style::bold
		          << "Generation cuts:" << rang::style::reset
		          << std::endl
		          << std::endl;
		printf("- System rapidity (Rap) [min, max] = [%0.2f, %0.2f]     (user) \n"
			   "- System mass (M)       [min, max] = [%0.2f, %0.2f] GeV (user) \n"
		       "- Forward leg (Pt)      [min, max] = [%0.2f, %0.2f] GeV (fixed/user) \n",
		    gcuts.Y_min, gcuts.Y_max,
		    gcuts.M_min, gcuts.M_max,
		    gcuts.forward_pt_min, gcuts.forward_pt_max);
		
		PrintFiducialCuts();
	}
}


// 5+1-dimensional phase space vector initialization
bool MFactorized::B51RandomKin(const std::vector<double>& randvec) {

	const double pt1  = gcuts.forward_pt_min + (gcuts.forward_pt_max - gcuts.forward_pt_min) * randvec[0];
	const double pt2  = gcuts.forward_pt_min + (gcuts.forward_pt_max - gcuts.forward_pt_min) * randvec[1];
	const double phi1 = 2.0 * gra::math::PI * randvec[2];
	const double phi2 = 2.0 * gra::math::PI * randvec[3];
	const double yX   = gcuts.Y_min + (gcuts.Y_max - gcuts.Y_min) * randvec[4];

	// Pick daughter masses, can fail due to off-shelliness, then
	// retry
	uint trials = 0;
	const uint MAXTRIAL = 1e5;
	while (true) {
		double M_sum = 0.0;
		
		// ==============================================================
		for (const auto& i : indices(lts.decaytree)) {
			if (FIXONSHELL) { // Use only on-shell value
				lts.decaytree[i].m_offshell = lts.decaytree[i].p.mass;
			} else {
				GetOffShellMass(lts.decaytree[i], lts.decaytree[i].m_offshell);
			}
			M_sum += lts.decaytree[i].m_offshell;
		}
		// ==============================================================

		// Apply absolute boundary conditions first
		const double MARGIN = 1e-4; // GeV
		M_MIN = M_sum + MARGIN;
		M_MAX = lts.sqrt_s - (lts.pbeam1.M() + lts.pbeam2.M());
		
		/*
		// If a single resonance, Now tighten >>
		if (ProcPtr.CHANNEL == "RES") {
			PARAM_RES res;
			for (auto const& x : lts.RESONANCES) {
				res = x.second;
			}
			if (res.p.mass > 0) {
				M_MIN = std::max(M_MIN, res.p.mass - res.p.width * 5);
				M_MAX = std::min(M_MAX, res.p.mass + res.p.width * 5);
			}
		}
		*/
		// Apply generator cuts (either automatic or user provided)
		M_MIN = std::max(M_MIN, gcuts.M_min);
		M_MAX = std::min(M_MAX, gcuts.M_max);

		if (M_MIN < M_MAX) { break; }
		++trials;
		if (trials > MAXTRIAL) {
			throw std::invalid_argument("MFactorized::B51RandomKin: Infinite loop in kinematics. Check the decaymode and cuts!");
		}
	};

	// Mass squared
	const double M2 = pow2(M_MIN) + (pow2(M_MAX) - pow2(M_MIN)) * randvec[5];

	return B51BuildKin(pt1, pt2, phi1, phi2, yX, M2);
}


// Build kinematics for 2->3 skeleton
bool MFactorized::B51BuildKin(double pt1, double pt2, double phi1, double phi2, double yX, double m2X) {

	// Forward N* system masses
	double m1 = 0;
	double m2 = 0;
	MFragment::GetForwardMass(m1, m2, lts.excite1, lts.excite2, EXCITATION, random);

	// Final state 4-momenta, set px,py first
	M4Vec p1(pt1 * std::cos(phi1), pt1 * std::sin(phi1), 0, 0);
	M4Vec p2(pt2 * std::cos(phi2), pt2 * std::sin(phi2), 0, 0);
	M4Vec pX(-(p1.Px() + p2.Px()), -(p1.Py() + p2.Py()), 0, 0);

	// Central system pz and E
	const double mtX = msqrt(m2X + pX.Pt2());
	pX.SetPzE(mtX * std::sinh(yX), mtX * std::cosh(yX));

	// Energy overflow
	if (pX.E() > (lts.sqrt_s - (m1 + m2))) { return false; }

	double p1z = gra::kinematics::SolvepPZ1(m1, m2, pt1, pt2, pX.Pz(), pX.E(), lts.sqrt_s);
	double p2z = -(pX.Pz() + p1z); // by momentum conservation

	// *** Enforce scattering direction +p -> +p, -p -> -p *** (RARE POLYNOMIAL BRANCH FLIP)
	p1z = std::abs(p1z);
	p2z = -std::abs(p2z);

	// Pz and E of forward protons/N*
	p1.SetPzE(p1z, msqrt(pow2(m1) + pow2(pt1) + pow2(p1z)));
	p2.SetPzE(p2z, msqrt(pow2(m2) + pow2(pt2) + pow2(p2z)));

	// Save
	lts.pfinal[1] = p1;
	lts.pfinal[2] = p2;
	lts.pfinal[0] = pX; // Central system

	// -------------------------------------------------------------------
	// Kinematic checks

	// Total 4-momentum conservation
	if (!CheckEMC((lts.pbeam1 + lts.pbeam2) - (lts.pfinal[1] + lts.pfinal[2] + lts.pfinal[0]))) { return false; }

	// ==============================================================================
	// Central system decay tree first branch kinematics set up here, the
	// rest is done recursively

	// Mother mass
	std::vector<double> masses;

	// Collect decay product masses
	for (const auto& i : indices(lts.decaytree)) {
		// @@ Note, we need to take offshell masses here @@
		masses.push_back(lts.decaytree[i].m_offshell);
	}
	std::vector<M4Vec> products;
	
	// false if amplitude has dependence on the final state legs (generic),
	// true if amplitude is a function of central system kinematics only (limited)
	const bool UNWEIGHT = ISOLATE; // Use ISOLATE tag
	
	gra::kinematics::MCW w;
	// 2-body
	if        (lts.decaytree.size() == 2) {
		w = gra::kinematics::TwoBodyPhaseSpace(
		    lts.pfinal[0], msqrt(m2X), masses, products, random.rng);
	// 3-body
	} else if (lts.decaytree.size() == 3) {
		w = gra::kinematics::ThreeBodyPhaseSpace(
		    lts.pfinal[0], msqrt(m2X), masses, products, UNWEIGHT, random.rng);
	// N-body
	} else if (lts.decaytree.size() > 3) {
		w = gra::kinematics::NBodyPhaseSpace(
		    lts.pfinal[0], msqrt(m2X), masses, products, UNWEIGHT, random.rng);
	}
	
	if (w.GetW() < 0){
		return false; // Kinematically impossible
	}
	lts.DW = w;

	// Collect decay products
	const uint offset = 3;
	for (const auto& i : indices(lts.decaytree)) {
		lts.decaytree[i].p4  = products[i];
		lts.pfinal[i+offset] = products[i]; 
	}

	// Treat decaytree recursively
	for (const auto& i : indices(lts.decaytree)) {
		if (!ConstructDecayKinematics(lts.decaytree[i])) { return false; }
	}

	// Forward excitation
	if (lts.excite1) { ExciteNstar(lts.pfinal[1], lts.decayforward1); }
	if (lts.excite2) { ExciteNstar(lts.pfinal[2], lts.decayforward2); }

	// ==============================================================================
	// Check that we are above mass threshold -> not necessary, this is
	// done in mass sampling function
	const uint Nf = lts.decaytree.size() + 2;
	return GetLorentzScalars(Nf);
}


// Calculate Pure Phase Space Decay Width (Volume)
void MFactorized::DecayWidthPS(double& exact) const {
	
	exact = 0;
	// Massive exact closed form for 2-body case
	if (lts.decaytree.size() == 2) {
		exact = gra::kinematics::PS2Massive(lts.m2,
										pow2(lts.decaytree[0].p4.M()),
		    							pow2(lts.decaytree[1].p4.M()) );
	}
	// Massless case
	if (lts.decaytree.size() > 2) {
		exact = gra::kinematics::PSnMassless(lts.m2, lts.decaytree.size());
	}
}


// 5-Dim Integral Volume \int_{M^2_MIN}^{M^2_MAX} { dM^2 } [phi1] x [phi2] x
// [pt1] x [pt2] x [y]
// Integral over central mass^2 is separate [phase-space factorization], but
// encapsulated here (+1 dimension)
double MFactorized::B51IntegralVolume() const {
	return (pow2(M_MAX) - pow2(M_MIN)) *
		   (2.0 * gra::math::PI) *
	       (2.0 * gra::math::PI) *
	       (gcuts.forward_pt_max - gcuts.forward_pt_min) *
	       (gcuts.forward_pt_max - gcuts.forward_pt_min) *
	       (gcuts.Y_max - gcuts.Y_min);
}


// 5-Dim phase space weight
double MFactorized::B51PhaseSpaceWeight() const {
	const double J = 1.0 /
	    std::abs(lts.pfinal[1].Pz() / lts.pfinal[1].E() -
	             lts.pfinal[2].Pz() / lts.pfinal[2].E()); // Jacobian, close to 0.5
	const double factor = (1.0 / 2.0) *
						  (1.0 / pow5(2.0 * gra::math::PI)) *
	                      lts.pfinal[1].Pt() *
	                      lts.pfinal[2].Pt() *
	                      (1.0 / (2.0 * lts.pfinal[1].E())) *
	                      (1.0 / (2.0 * lts.pfinal[2].E()));
	return J * factor;
}

// For high mass limit kinematics, see e.g. [arxiv.org/pdf/hep-ph/9903279.pdf]

} // gra namespace ends
