// Continuum type phase space class
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.


// C++
#include <algorithm>
#include <complex>
#include <iostream>
#include <random>
#include <vector>

// OWN
#include "Graniitti/MAux.h"
#include "Graniitti/MContinuum.h"
#include "Graniitti/MForm.h"
#include "Graniitti/MKinematics.h"
#include "Graniitti/MMath.h"
#include "Graniitti/MMatrix.h"
#include "Graniitti/MProcess.h"
#include "Graniitti/MSpin.h"
#include "Graniitti/MUserCuts.h"
#include "Graniitti/MFragment.h"

// Libraries
#include "rang.hpp"

using gra::aux::indices;
using gra::PDG::GeV2barn;
using gra::math::abs2;
using gra::math::msqrt;
using gra::math::pow2;
using gra::math::pow3;
using gra::math::pow4;
using gra::math::zi;
using gra::math::PI;


namespace gra {


// This is needed by construction
MContinuum::MContinuum() {
	std::vector<std::string> supported = {"PP","yP","yy","gg"};
	ProcPtr = MSubProc(supported);
	ConstructProcesses();
}


// Constructor
MContinuum::MContinuum(std::string process) {
	
	InitHistograms();
	SetProcess(process);

	// Init final states
	M4Vec zerovec(0, 0, 0, 0);
	for (std::size_t i = 0; i < 10; ++i) {
		lts.pfinal.push_back(zerovec);
	}
	std::cout << "MContinuum:: [Constructor done]" << std::endl;
}

// Destructor
MContinuum::~MContinuum() {
}

// This is manually constructed and updated here
void MContinuum::ConstructProcesses() {

	Processes.clear();

	const std::string CID = "C";
	for (auto const& x : ProcPtr.descriptions) {
	    std::map<std::string, std::string> value = x.second;
	    for (const auto& y : value) {
	    	Processes.insert(std::make_pair(x.first + "[" + y.first + "]<"+CID+">", y.second));
		}
	}
}


// Initialize cut and process spesific postsetup
void MContinuum::post_Constructor() {

	std::cout << "MContinuum::post_Constructor: " << EXCITATION << std::endl;
	
	// Set sampling boundaries
	ProcPtr.SetTechnicalBoundaries(gcuts, EXCITATION);
	
	// Technical cuts here
	if (gcuts.kt_max < 0.0) { // not set by user

		gcuts.kt_min = 0.0;

		// Pomeron-Pomeron continuum (low-pt guaranteed by form factors)
		gcuts.kt_max = 10.0 / lts.decaytree.size();

		// QED/EW/Durham-QCD processes (high pt)
		if (ProcPtr.ISTATE != "PP" && ProcPtr.ISTATE != "yP") {
			gcuts.kt_max = lts.sqrt_s / 2.0;
		}
	}

	// Initialize phase space dimension (3*Nf - 4)
	ProcPtr.LIPSDIM = 3*(lts.decaytree.size()+2)-4;

	if (EXCITATION == 1) { ProcPtr.LIPSDIM += 1; }
	if (EXCITATION == 2) { ProcPtr.LIPSDIM += 2; }
}


// Update kinematics (screening kT loop calls this)
// Exact 4-momentum conservation at loop vertices, that is, by using this one
// does not assume vanishing external momenta in the screening loop calculation
bool MContinuum::LoopKinematics(const std::vector<double>& p1p,
                                const std::vector<double>& p2p) {

	const unsigned int Kf = lts.decaytree.size(); // Number of central particles
	const unsigned int Nf = Kf + 2;               // Number of final states

	// SET new proton pT degrees of freedom
	lts.pfinal[1].SetPxPy(p1p[0],p1p[1]);
	lts.pfinal[2].SetPxPy(p2p[0],p2p[1]);

	// Get central final states pT degrees of freedom
	std::vector<M4Vec> p(Kf, M4Vec(0,0,0,0));
	BLinearSystem(p, pkt_, lts.pfinal[1], lts.pfinal[2]);

	// Set central particles px,py,pz,e
	const unsigned int offset = 3; // indexing
	M4Vec sumP(0,0,0,0);
	for (std::size_t i = 0; i < Kf; ++i) {

		const double y = lts.pfinal_orig[i+offset].Rap();
		const double m = lts.pfinal_orig[i+offset].M();

		// px,py
		lts.pfinal[i+offset].SetPxPy(p[i].Px(), p[i].Py());

		// pz,E
		const double mt = msqrt(pow2(m) + lts.pfinal[i+offset].Pt2());
		lts.pfinal[i+offset].SetPzE(mt * std::sinh(y), mt * std::cosh(y));

		sumP  += lts.pfinal[i+offset];
	}

	// pz and E of forward protons/N*
	const double pt1 = lts.pfinal[1].Pt();
	const double pt2 = lts.pfinal[2].Pt();
	const double m1  = lts.pfinal_orig[1].M();
	const double m2  = lts.pfinal_orig[2].M();
	
	double p1z = gra::kinematics::SolvepPZ1(m1, m2, pt1, pt2, sumP.Pz(), sumP.E(), lts.sqrt_s);
	double p2z = -(sumP.Pz() + p1z); // by momentum conservation
	
	// Enforce scattering direction +p -> +p, -p -> -p (VERY RARE POLYNOMIAL BRANCH FLIP)
	if (p1z < 0 || p2z > 0) { return false; }
	
	lts.pfinal[1].SetPzE(p1z, msqrt(pow2(m1) + pow2(pt1) + pow2(p1z)));
	lts.pfinal[2].SetPzE(p2z, msqrt(pow2(m2) + pow2(pt2) + pow2(p2z)));
	lts.pfinal[0] = sumP;

	// Check Energy-Momentum
	if (!gra::math::CheckEMC((lts.pbeam1 + lts.pbeam2) - (lts.pfinal[1] + lts.pfinal[2] + lts.pfinal[0]))) { return false; }

	return GetLorentzScalars(Nf);
}


// Get weight
double MContinuum::EventWeight(const std::vector<double>& randvec, AuxIntData& aux) {
	
	double W = 0.0;

	// Kinematics
	aux.kinematics_ok = BNRandomKin(lts.decaytree.size() + 2, randvec);
	aux.fidcuts_ok    = FiducialCuts();
	aux.vetocuts_ok   = VetoCuts();

	if (aux.Valid()) {
		
		// ** EVENT WEIGHT **
		double MatESQ = (FLATAMPLITUDE == 0) ? abs2(S3ScreenedAmp()) : GetFlatAmp2(lts); // Matrix element squared

		W = (1.0 / S_factor) * BNPhaseSpaceWeight() * BNIntegralVolume() * MatESQ * GeV2barn /
		    MollerFlux(); // Total weight: phase-space x |M|^2 x barn units
	}
	
	aux.amplitude_ok = CheckInfNan(W);
	
	// As the last step: Histograms
	if (!aux.burn_in_mode) {
		const double totalweight = W * aux.vegasweight;
		FillHistograms(totalweight, lts);
	}

	return W;
}

// Fiducial cuts
bool MContinuum::FiducialCuts() const {

	return CommonCuts();
	
}

// Record event
bool MContinuum::EventRecord(HepMC3::GenEvent& evt) {

	return CommonRecord(evt);
}

void MContinuum::PrintInit(bool silent) const {
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

		std::string pomeronloop =           
		    (SCREENING == true) ? "     **    " : "           ";

		std::vector<std::string> legs;
		for (const auto& i : indices(lts.decaytree)) {
			char buff[250];
			snprintf(
			    buff, sizeof(buff),
			    "x---------> %d (%s) [Q=%s, J=%s]",
			    lts.decaytree[i].p.pdg,
			    lts.decaytree[i].p.name.c_str(),
			    gra::aux::Charge3XtoString(lts.decaytree[i].p.chargeX3)
			        .c_str(),
			    gra::aux::Spin2XtoString(lts.decaytree[i].p.spinX2)
			        .c_str());
			std::string leg = buff;
			legs.push_back(leg);
		}

		std::vector<std::string> feynmangraph;
		feynmangraph.push_back("||         ");
		for (const auto& i : indices(lts.decaytree)) {
			feynmangraph.push_back(legs[i]);
			if (i < lts.decaytree.size() - 1)
				feynmangraph.push_back("|          ");
		}
		feynmangraph.push_back("||         ");

		// Print diagram
		std::cout << proton1 << std::endl;
		for (const auto& i : indices(feynmangraph)) {
			if (SCREENING) { // Put red
				std::cout << rang::fg::red << pomeronloop << rang::style::reset;
			} else {
				std::cout << rang::fg::red << pomeronloop << rang::style::reset;
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
		
		// Generation cuts
		std::cout << std::endl;
		std::cout << rang::style::bold
		          << "Generation cuts:" << rang::style::reset
		          << std::endl
		          << std::endl;
		printf("- Final state rapidity (Rap) [min, max] = [%0.2f, %0.2f]     \t(user)       \n"
		       "- Intermediate (Kt)          [min, max] = [%0.2f, %0.2f] GeV \t(fixed/user) \n"
		       "- Forward leg (Pt)           [min, max] = [%0.2f, %0.2f] GeV \t(fixed/user) \n",
		       gcuts.rap_min, gcuts.rap_max, gcuts.kt_min, gcuts.kt_max, gcuts.forward_pt_min, gcuts.forward_pt_max);

		if (EXCITATION != 0) {
		printf("- Forward leg (Xi)           [min, max] = [%0.2f, %0.2f]     \t(fixed/user) \n", gcuts.XI_min, gcuts.XI_max);
		}

		PrintFiducialCuts();
	}
}


// (3*Nf-4)-dimensional phase space vector initialization
bool MContinuum::BNRandomKin(unsigned int Nf, const std::vector<double>& randvec) {

	const unsigned int Kf = Nf - 2; // Central system multiplicity

	const double pt1  = gcuts.forward_pt_min + (gcuts.forward_pt_max - gcuts.forward_pt_min) * randvec[0];
	const double pt2  = gcuts.forward_pt_min + (gcuts.forward_pt_max - gcuts.forward_pt_min) * randvec[1];
	const double phi1 = 2.0 * gra::math::PI * randvec[2];
	const double phi2 = 2.0 * gra::math::PI * randvec[3];
	const unsigned int offset = 4; // 4 variables above

	// Decay product masses
	// ==============================================================
	for (const auto& i : indices(lts.decaytree)) {
		if (FIXONSHELL) { // Use only on-shell value
			lts.decaytree[i].m_offshell = lts.decaytree[i].p.mass;
		} else {
			GetOffShellMass(lts.decaytree[i], lts.decaytree[i].m_offshell);
		}
	}
	// ==============================================================

	// Intermediate kt
	std::vector<double> kt(Kf - 1, 0.0);  // Kf-1
	size_t ind = offset;
	for (const auto& i : indices(kt)) {
		kt[i] = gcuts.kt_min + (gcuts.kt_max - gcuts.kt_min) * randvec[ind];
		++ind;
	}

	// Intermediate phi
	std::vector<double> phi(Kf - 1, 0.0); // Kf-1
	for (const auto& i : indices(phi)) {
		phi[i] = 2.0 * PI * randvec[ind];
		++ind;
	}

	// Final state rapidity
	std::vector<double> y(Kf, 0.0);       // Kf
	for (const auto& i : indices(y)) {
		y[i] = gcuts.rap_min + (gcuts.rap_max - gcuts.rap_min) * randvec[ind];
		++ind;
	}

	// Forward N* system masses
	double m1 = beam1.mass;
	double m2 = beam2.mass;
	lts.excite1 = false;
	lts.excite2 = false;

	M2_f_min = pow2(1.07);
	M2_f_max = gcuts.XI_max * lts.s;
	
	if (EXCITATION == 1) {
		const double mforward = msqrt(M2_f_min + (M2_f_max - M2_f_min) * randvec[ind]);
		if (random.U(0,1) < 0.5) {
			m1 = mforward;
			lts.excite1 = true;
		} else {
			m2 = mforward;
			lts.excite2 = true;
		}
	}
	if (EXCITATION == 2) {
		m1 = msqrt(M2_f_min + (M2_f_max - M2_f_min) * randvec[ind]  );
		m2 = msqrt(M2_f_min + (M2_f_max - M2_f_min) * randvec[ind+1]);
		lts.excite1 = true;
		lts.excite2 = true;
	}

	return BNBuildKin(Nf, pt1, pt2, phi1, phi2, kt, phi, y, m1, m2);
}


// Build kinematics of 2->N
bool MContinuum::BNBuildKin(unsigned int Nf, double pt1, double pt2, double phi1, double phi2,
	const std::vector<double>& kt, const std::vector<double>& phi, const std::vector<double>& y, double m1, double m2) {

	const unsigned int Kf = Nf - 2; // Central system multiplicity

	if (lts.decaytree.size() != Kf) {
		std::string str =
		    "MContinuum::BNBuildKin: Decaytree first level topology "
		    "= " + std::to_string(lts.decaytree.size()) + 
		    " (should be " + std::to_string(Kf) +" for this process)!";
		throw std::invalid_argument(str);
	}

	//double m1 = 0;
	//double m2 = 0;
	//MFragment::GetForwardMass(m1, m2, lts.excite1, lts.excite2, EXCITATION, random);
	
	// Forward protons px,py
	M4Vec p1(pt1 * std::cos(phi1), pt1 * std::sin(phi1), 0, 0);
	M4Vec p2(pt2 * std::cos(phi2), pt2 * std::sin(phi2), 0, 0);

	// Auxialary "difference momentum" q0 = p0 - p1 ...
	pkt_.resize(Kf-1);
	for (const auto& i : indices(pkt_)) {
		pkt_[i] = M4Vec(kt[i] * std::cos(phi[i]), kt[i] * std::sin(phi[i]), 0, 0);
	}

	// Apply linear system to get p
	std::vector<M4Vec> p(Kf, M4Vec(0,0,0,0));
	BLinearSystem(p, pkt_, p1, p2);

	// Set pz and E for central final states
	M4Vec sumP(0,0,0,0);
	for (const auto& i : indices(p)) {

		const double m = lts.decaytree[i].m_offshell; // Note offshell!
		const double mt = msqrt(pow2(m) + p[i].Pt2());
		p[i].SetPzE(mt*std::sinh(y[i]), mt*std::cosh(y[i]));

		sumP += p[i];
	}

	// Check crude energy overflow
	if (sumP.E() > lts.sqrt_s) { return false; }

	double p1z = gra::kinematics::SolvepPZ1(m1, m2, pt1, pt2, sumP.Pz(), sumP.E(), lts.sqrt_s);
	double p2z = -(sumP.Pz() + p1z); // by momentum conservation

	// Enforce scattering direction +p -> +p, -p -> -p (VERY RARE POLYNOMIAL BRANCH FLIP)
	if (p1z < 0 || p2z > 0) { return false; }

	// pz and E of protons/N*
	p1.SetPzE(p1z, msqrt(pow2(m1) + pow2(pt1) + pow2(p1z)));
	p2.SetPzE(p2z, msqrt(pow2(m2) + pow2(pt2) + pow2(p2z)));

	// First branch kinematics
	lts.pfinal[1] = p1;    // Forward systems
	lts.pfinal[2] = p2;
	lts.pfinal[0] = sumP;  // Central system

	double sumM = 0;
	const unsigned int offset = 3;
	for (const auto& i : indices(p)) {
		lts.decaytree[i].p4  = p[i];
		lts.pfinal[i+offset] = p[i];
		sumM += p[i].M();
	}

	// -------------------------------------------------------------------
	// Kinematic checks

	// Check we are above mass threshold
	if (sumP.M() < sumM) { return false; }

	// Total 4-momentum conservation
	if (!gra::math::CheckEMC((lts.pbeam1 + lts.pbeam2) - (lts.pfinal[1] + lts.pfinal[2] + lts.pfinal[0]))) { return false; }

	// -------------------------------------------------------------------

	// Treat decaytree recursively
	for (const auto& i : indices(lts.decaytree)) {
		if (!ConstructDecayKinematics(lts.decaytree[i])) { return false; }
	}
	
	// Forward excitation
	if (lts.excite1) { ExciteNstar(lts.pfinal[1], lts.decayforward1); }
	if (lts.excite2) { ExciteNstar(lts.pfinal[2], lts.decayforward2); }

	return GetLorentzScalars(Nf);
}


// Construct the linear system
void MContinuum::BLinearSystem(std::vector<M4Vec>& p, const std::vector<M4Vec>& q,
	const M4Vec& p1f, const M4Vec& p2f) const {

	/*

	kk = p1f+p2f
	
	-------------------------------
	KF = 4

	2p0 = ( q0 - kk - p2 - p3)
	2p1 = (-q0 - kk - p2 - p3)
	2p2 = (-q1 - kk - p0 - p3)
	2p3 = (-q2 - kk - p0 - p1)
	
	[2 0 1 1
	 0 2 1 1
	 1 0 2 1
	 1 1 0 2]

	-------------------------------
	KF = 3

	2p0 = ( q0 - kk - p2)
	2p1 = (-q0 - kk - p2)
	2p2 = (-q1 - kk - p0)

	[2 0 1
	 0 2 1
	 1 0 2]

	-------------------------------
	KF = 2

	2p0 = ( q0 - kk)
	2p1 = (-q0 - kk)

	[2 0
	 0 2]

	-------------------------------
	*/

	static const std::vector<std::vector<std::vector<double>>> A = 
	{ 
	{{1.0/2.0, 0.0},
	{0.0, 1.0/2.0}},

	{{2.0/3.0, 0.0, -1.0/3.0},
	{1.0/6.0, 1.0/2.0, -1.0/3.0},
	{-1.0/3.0, 0.0, 2.0/3.0}},

	{{7.0/8.0, 1.0/8.0, -1.0/2.0, -1.0/4.0},
	{3.0/8.0, 5.0/8.0, -1.0/2.0, -1.0/4.0},
	{-1.0/8.0, 1.0/8.0, 1.0/2.0, -1.0/4.0},
	{-5.0/8.0, -3.0/8.0, 1.0/2.0, 3.0/4.0}},

	{{11.0/10.0, 3.0/10.0, -3.0/5.0, -2.0/5.0, -1.0/5.0},
	{3.0/5.0, 4.0/5.0, -3.0/5.0, -2.0/5.0, -1.0/5.0},
	{1.0/10.0, 3.0/10.0, 2.0/5.0, -2.0/5.0, -1.0/5.0},
	{-2.0/5.0, -1.0/5.0, 2.0/5.0, 3.0/5.0, -1.0/5.0},
	{-9.0/10.0, -7.0/10.0, 2.0/5.0, 3.0/5.0, 4.0/5.0}},

	{{4.0/3.0, 1.0/2.0, -2.0/3.0, -1.0/2.0, -1.0/3.0, -1.0/6.0},
	{5.0/6.0, 1.0/1.0, -2.0/3.0, -1.0/2.0, -1.0/3.0, -1.0/6.0},
	{1.0/3.0, 1.0/2.0, 1.0/3.0, -1.0/2.0, -1.0/3.0, -1.0/6.0},
	{-1.0/6.0, 0.0, 1.0/3.0, 1.0/2.0, -1.0/3.0, -1.0/6.0},
	{-2.0/3.0, -1.0/2.0, 1.0/3.0, 1.0/2.0, 2.0/3.0, -1.0/6.0},
	{-7.0/6.0, -1.0/1.0, 1.0/3.0, 1.0/2.0, 2.0/3.0, 5.0/6.0}},

	{{11.0/7.0, 5.0/7.0, -5.0/7.0, -4.0/7.0, -3.0/7.0, -2.0/7.0, -1.0/7.0},
	{15.0/14.0, 17.0/14.0, -5.0/7.0, -4.0/7.0, -3.0/7.0, -2.0/7.0, -1.0/7.0},
	{4.0/7.0, 5.0/7.0, 2.0/7.0, -4.0/7.0, -3.0/7.0, -2.0/7.0, -1.0/7.0},
	{1.0/14.0, 3.0/14.0, 2.0/7.0, 3.0/7.0, -3.0/7.0, -2.0/7.0, -1.0/7.0},
	{-3.0/7.0, -2.0/7.0, 2.0/7.0, 3.0/7.0, 4.0/7.0, -2.0/7.0, -1.0/7.0},
	{-13.0/14.0, -11.0/14.0, 2.0/7.0, 3.0/7.0, 4.0/7.0, 5.0/7.0, -1.0/7.0},
	{-10.0/7.0, -9.0/7.0, 2.0/7.0, 3.0/7.0, 4.0/7.0, 5.0/7.0, 6.0/7.0}},

	{{29.0/16.0, 15.0/16.0, -3.0/4.0, -5.0/8.0, -1.0/2.0, -3.0/8.0, -1.0/4.0, -1.0/8.0},
	{21.0/16.0, 23.0/16.0, -3.0/4.0, -5.0/8.0, -1.0/2.0, -3.0/8.0, -1.0/4.0, -1.0/8.0},
	{13.0/16.0, 15.0/16.0, 1.0/4.0, -5.0/8.0, -1.0/2.0, -3.0/8.0, -1.0/4.0, -1.0/8.0},
	{5.0/16.0, 7.0/16.0, 1.0/4.0, 3.0/8.0, -1.0/2.0, -3.0/8.0, -1.0/4.0, -1.0/8.0},
	{-3.0/16.0, -1.0/16.0, 1.0/4.0, 3.0/8.0, 1.0/2.0, -3.0/8.0, -1.0/4.0, -1.0/8.0},
	{-11.0/16.0, -9.0/16.0, 1.0/4.0, 3.0/8.0, 1.0/2.0, 5.0/8.0, -1.0/4.0, -1.0/8.0},
	{-19.0/16.0, -17.0/16.0, 1.0/4.0, 3.0/8.0, 1.0/2.0, 5.0/8.0, 3.0/4.0, -1.0/8.0},
	{-27.0/16.0, -25.0/16.0, 1.0/4.0, 3.0/8.0, 1.0/2.0, 5.0/8.0, 3.0/4.0, 7.0/8.0}}
	};

	// Construct vector b
	const unsigned int Kf = p.size(); // Number of central system particles
	std::vector<M4Vec> b(Kf);
	const M4Vec p1p2sum = p1f + p2f;

	for (const auto& i : indices(b)) {
		if (i == 0){
			b[i] = q[0] - p1p2sum;
		} else {
			b[i] = q[i-1]*(-1.0) - p1p2sum;
		}
	}

	// Apply linear system p = Ab to get px,py components for each p[i]
	const unsigned int index = Kf-2; // -2 because of C++
	for (const auto& i : indices(p)) {
		for (const auto& j : indices(b)) {
			p[i] += b[j] * A[index][i][j]; // notice plus
		}
	}

}

// ----------------------------------------------------------------------
// MATLAB symbolic code to print out the system
// ----------------------------------------------------------------------
/*
% Number of particles in the central system
KFMAX = 8;

fprintf('static const std::vector<std::vector<std::vector<double>>> A = \n{ ');
for KF = 2:KFMAX

% Create system matrix A
A = ones(KF);
A(1:2,1:2) = [2 0; 0 2]; % Starting corner

% Lower triangle
for i = 1:KF-2
   m = ones(1,KF);
   m(i+1) = 0;
   m(i+2) = 2;
   A(2+i,:) = m;
end

% Invert
I = inv(sym(A));

% Plot C++ style
fprintf('\n');
fprintf('{');
for i = 1:size(A,1)
    fprintf('{');
    for j = 1:size(A,2)
        if (j < size(A,2))
            mark = ', ';
        else
            mark = '';
        end
        [N,D] = numden(I(i,j));
        if (N == 0)
        fprintf('0.0%s', mark);
        else
        fprintf('%0.1f/%0.1f%s', N, D, mark);
        end
    end
    if (i < size(A,2))
        mark = '},';
    else
        if (KF < KFMAX)
            emark = ',';
        else
            emark = '';
        end
        mark = sprintf('}}%s', emark);
    end
    fprintf('%s\n', mark);
end

end
fprintf('};\n');
*/

// (3*Nf-4)-Dim Integral Volume [phi1] x [phi2] x [phi_k1] x [phi_k2] x ... x [phi_{Kf-1}]
//                        [pt1] x [pt2] x [kt1] x [kt2] x ... x [kt{Kf-1}]
//                        x [y3] x [y4] x ... x [yN]
//
// For reference, Nf = 4 special case in:
// [REFERENCE: Lebiedowicz, Szczurek, arxiv.org/abs/0912.0190]
double MContinuum::BNIntegralVolume() const {
	
	// Number of central states
	const unsigned int Kf = pkt_.size() + 1;
	
	// Forward leg integration
	double M2_forward_volume = 1.0;
	
	if      (EXCITATION == 1) {
		M2_forward_volume = M2_f_max - M2_f_min;
	}
	else if (EXCITATION == 2) {
		M2_forward_volume = pow2(M2_f_max - M2_f_min);
	}
	
	return pow2(2.0 * PI) * 
		   pow2(gcuts.forward_pt_max - gcuts.forward_pt_min) *
		   std::pow(2.0 * PI, Kf-1) *
	       std::pow(gcuts.kt_max  - gcuts.kt_min, Kf-1) *
	       std::pow(gcuts.rap_max - gcuts.rap_min, Kf) *
	       M2_forward_volume;
}


double MContinuum::BNPhaseSpaceWeight() const {
	
	// Jacobian, close to constant 0.5
	const double J = 1.0 /
	    std::abs(lts.pfinal[1].Pz() / lts.pfinal[1].E() -
	             lts.pfinal[2].Pz() / lts.pfinal[2].E());

    // Number of final states
    const unsigned int Nf = lts.decaytree.size()+2;

	// Intermediate "difference"
	// kt factors from \prod_i d^2 k_i = \prod_i dphi_i kt_i dkt_i
    double PROD = 1.0;
	for (const auto& i : indices(pkt_)) {
		PROD *= pkt_[i].Pt();
	}

	const double factor =
	    pow4(2.0 * PI) *
	    (1.0 / std::pow(2.0 * pow3(2.0 * PI), Nf) ) *
	    lts.pfinal[1].Pt() *
	    lts.pfinal[2].Pt() *
	    (1.0 / (2.0 * lts.pfinal[1].E())) *
	    (1.0 / (2.0 * lts.pfinal[2].E())) * 
	    PROD * (1.0 / std::pow(2,Nf-4));
	return J * factor;
}

} // gra namespace ends
