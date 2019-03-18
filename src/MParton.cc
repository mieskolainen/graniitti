// Collinear parton model 2->2 or (2->N) type phase space class
// for calculations without proton (remnant) considerations
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
#include "Graniitti/MParton.h"
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

// This is needed by construction
MParton::MParton() {
	std::vector<std::string> supported = {"yy_LUX", "yy_DZ"};
	ProcPtr = MSubProc(supported);
	ConstructProcesses();
}

// Constructor
MParton::MParton(std::string process) {

	InitHistograms();
	SetProcess(process);

	// Init final states
	M4Vec zerovec(0, 0, 0, 0);
	for (std::size_t i = 0; i < 10; ++i) {
		lts.pfinal.push_back(zerovec);
	}
	std::cout << "MParton:: [Constructor done]" << std::endl;
}

// Destructor
MParton::~MParton() {
}

// This is manually constructed and updated here
void MParton::ConstructProcesses() {

	Processes.clear();
	CID = "P";
	for (auto const& x : ProcPtr.descriptions) {
	    std::map<std::string, std::string> value = x.second;
	    for (auto const& y : value) {
	    	Processes.insert(std::make_pair(x.first + "[" + y.first + "]<"+CID+">", y.second));
		}
	}
}


// Initialize cut and process spesific postsetup
void MParton::post_Constructor() {
	
	if (ProcPtr.CHANNEL == "RES") {
		// Here we support only single resonances
		if (lts.RESONANCES.size() != 1) {
			std::string str = "MParton::post_Constructor: Only single resonance supported for this process (RESPARAM.size() != 1)";
			throw std::invalid_argument(str);
		}
	}
	
	// Set sampling boundaries
	ProcPtr.SetTechnicalBoundaries(gcuts, EXCITATION);
	
	// Initialize phase space dimension
	ProcPtr.LIPSDIM = 2; // All processes
	
	//if (EXCITATION == 1) { ProcPtr.LIPSDIM += 1; }
	//if (EXCITATION == 2) { ProcPtr.LIPSDIM += 2; }
}


// Update kinematics (screening kT loop calls this)
// Exact 4-momentum conservation at loop vertices, that is, by using this one
// does not assume vanishing external momenta in the screening loop calculation.

bool MParton::LoopKinematics(const std::vector<double>& p1p,
                             const std::vector<double>& p2p) {
	return false;
}


// Return Monte Carlo integrand weight
double MParton::EventWeight(const std::vector<double>& randvec, AuxIntData& aux) {

	double W = 0.0;

	// Kinematics and cuts
	aux.kinematics_ok = B2RandomKin(randvec);
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
		 // We have some legs in the central system
		if (lts.decaytree.size() != 0 && !ISOLATE) {
			C_space = lts.DW.Integral();
		}
		
		// ** EVENT WEIGHT **
		W = C_space * (1.0 / S_factor) * MatESQ * GeV2barn *
			B2IntegralVolume() * B2PhaseSpaceWeight() / MollerFlux();
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
bool MParton::FiducialCuts() const {

	return CommonCuts();
}

// Record event
bool MParton::EventRecord(HepMC3::GenEvent& evt) {
	
	return CommonRecord(evt);
}

void MParton::PrintInit(bool silent) const {
	if (!silent) {
		PrintSetup();

		// Construct prettyprint diagram
		std::string proton1 = "-----------pdf-------->";
		std::string proton2 = "-----------pdf-------->";

		/*
		if (EXCITATION == 1) {
			proton1 = "-----------F2-xxxxxxxx>";
		}
		if (EXCITATION == 2) {
			proton1 = "-----------F2-xxxxxxxx>";
			proton2 = "-----------F2-xxxxxxxx>";
		}
		*/

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
		
		std::cout << "- None" << std::endl;
		/*
		printf("- System rapidity (Rap) [min, max] = [%0.2f, %0.2f]     \t(user) \n"
			   "- System mass (M)       [min, max] = [%0.2f, %0.2f] GeV \t(user) \n"
		       "- Forward leg (Pt)      [min, max] = [%0.2f, %0.2f] GeV \t(fixed/user) \n",
		    gcuts.Y_min, gcuts.Y_max,
		    gcuts.M_min, gcuts.M_max,
		    gcuts.forward_pt_min, gcuts.forward_pt_max);
		
		if (EXCITATION != 0) {
		printf("- Forward leg (Xi)      [min, max] = [%0.2E, %0.2E]     \t(fixed/user) \n", gcuts.XI_min, gcuts.XI_max);
		}
		*/

		PrintFiducialCuts();
	}
}

// 2-dimensional phase space vector initialization
bool MParton::B2RandomKin(const std::vector<double>& randvec) {

    const double x1 = randvec[0];
    const double x2 = randvec[1];

	return B2BuildKin(x1, x2);
}

// Build kinematics for 2->2 skeleton
bool MParton::B2BuildKin(double x1, double x2) {

	static const M4Vec beamsum = lts.pbeam1 + lts.pbeam2;

	// We-work in CMS-frame

	// Initial state collinear (pt=0) and massless (m=0) parton 4-momentum
	M4Vec q1(0, 0,  x1 * lts.sqrt_s/2, x1 * lts.sqrt_s/2);
	M4Vec q2(0, 0, -x2 * lts.sqrt_s/2, x2 * lts.sqrt_s/2);

	// ------------------------------------------------------------------
	// Now boost if asymmetric beams
	if (std::abs(beamsum.Pz()) > 1e-9) {
		constexpr int sign = 1; // positive -> boost to the lab
		kinematics::LorentzBoost(beamsum, lts.sqrt_s, q1, sign);
		kinematics::LorentzBoost(beamsum, lts.sqrt_s, q2, sign);
	}
	// ------------------------------------------------------------------
	
	M4Vec p1 = lts.pbeam1 - q1; // Remnant
	M4Vec p2 = lts.pbeam2 - q2; // Remnant
	M4Vec pX = q1 + q2;         // System

	// Save
	lts.pfinal[1] = p1;
	lts.pfinal[2] = p2;
	lts.pfinal[0] = pX; // Central system

	// -------------------------------------------------------------------
	// Kinematic checks

	// Total 4-momentum conservation
	if (!CheckEMC(beamsum - (lts.pfinal[1] + lts.pfinal[2] + lts.pfinal[0]))) { return false; }

	// ==============================================================================
	// Central system decay tree first branch kinematics set up here, the
	// rest is done recursively

	// Mother mass
	std::vector<double> masses;

	// Collect decay product masses
	for (const auto& i : indices(lts.decaytree)) {
		// @@ Note, we need to take offshell masses here @@
		//masses.push_back(lts.decaytree[i].m_offshell);
		masses.push_back(lts.decaytree[i].p.mass);	
	}
	std::vector<M4Vec> products;
	
	// false if amplitude has dependence on the final state legs (generic),
	// true if amplitude is a function of central system kinematics only (limited)
	const bool UNWEIGHT = ISOLATE; // Use ISOLATE tag
	
	gra::kinematics::MCW w;
	// 2-body
	if        (lts.decaytree.size() == 2) {
		w = gra::kinematics::TwoBodyPhaseSpace(
		    lts.pfinal[0], lts.pfinal[0].M(), masses, products, random);
	// 3-body
	} else if (lts.decaytree.size() == 3) {
		w = gra::kinematics::ThreeBodyPhaseSpace(
		    lts.pfinal[0], lts.pfinal[0].M(), masses, products, UNWEIGHT, random);
	// N-body
	} else if (lts.decaytree.size() > 3) {
		w = gra::kinematics::NBodyPhaseSpace(
		    lts.pfinal[0], lts.pfinal[0].M(), masses, products, UNWEIGHT, random);
	}
	
	if (w.GetW() < 0){
		return false; // Kinematically impossible
	}
	lts.DW = w;

	// Collect decay products
	const unsigned int offset = 3;
	for (const auto& i : indices(lts.decaytree)) {
		lts.decaytree[i].p4  = products[i];
		lts.pfinal[i+offset] = products[i]; 
	}

	// Treat decaytree recursively
	for (const auto& i : indices(lts.decaytree)) {
		if (!ConstructDecayKinematics(lts.decaytree[i])) { return false; }
	}

	// ==============================================================================
	// Check that we are above mass threshold -> not necessary, this is
	// done in mass sampling function
	const unsigned int Nf = lts.decaytree.size() + 2;
	return GetLorentzScalars(Nf);
}

// Calculate Pure Phase Space Decay Width (Volume)
void MParton::DecayWidthPS(double& exact) const {
	
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

// 2-Dim Integral Volume
double MParton::B2IntegralVolume() const {

	return 1.0;
}

// 2-Dim phase space weight
double MParton::B2PhaseSpaceWeight() const {
	
	return 1.0;
}

} // gra namespace ends