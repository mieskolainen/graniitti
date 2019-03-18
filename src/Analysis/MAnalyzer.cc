// Fast MC analysis class
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.


// C++
#include <algorithm>
#include <iostream>
#include <memory>
#include <vector>

// C-file processing
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

// ROOT
#include "TBranch.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TRandom3.h"
#include "TString.h"
#include "TStyle.h"
#include "TTree.h"

// HepMC3 3
#include "HepMC3/GenCrossSection.h"
#include "HepMC3/FourVector.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/Print.h"
#include "HepMC3/Selector.h"
#include "HepMC3/Relatives.h"


// Own
#include "Graniitti/Analysis/MAnalyzer.h"
#include "Graniitti/Analysis/MMultiplet.h"
#include "Graniitti/MMath.h"
#include "Graniitti/M4Vec.h"
#include "Graniitti/MKinematics.h"
#include "Graniitti/MPDG.h"

const bool DEBUG = false;

using gra::math::msqrt;
using gra::aux::indices;


namespace gra {


// Constructor
MAnalyzer::MAnalyzer() {

	// Initialize histograms
	const int NBINS = 150;

	// Energy
	hE_Pions        = std::make_unique<TH1D>("Energy #pi (GeV)",    ";Energy (GeV);Events", NBINS, 0, sqrts / 2.0);
	hE_Gamma        = std::make_unique<TH1D>("Energy #gamma (GeV)", ";Energy (GeV);Events", NBINS, 0, sqrts / 2.0);
	hE_Neutron      = std::make_unique<TH1D>("Energy n (GeV)",      ";Energy (GeV);Events", NBINS, 0, sqrts / 2.0);
	hE_GammaNeutron = std::make_unique<TH1D>("Energy y+n (GeV)",    ";Energy (GeV);Events", NBINS, 0, sqrts / 2.0);
	
	// Feynman-x
	hXF_Pions    = std::make_unique<TH1D>("xF #pi",    ";Feynman-x;Events", NBINS, -1.0, 1.0);
	hXF_Gamma    = std::make_unique<TH1D>("xF #gamma", ";Feynman-x;Events", NBINS, -1.0, 1.0);
	hXF_Neutron  = std::make_unique<TH1D>("xF n",      ";Feynman-x;Events", NBINS, -1.0, 1.0);

	// Forward systems
	hEta_Pions   = std::make_unique<TH1D>("#eta pi", ";#eta;Events", NBINS, -12, 12);
	hEta_Gamma   = std::make_unique<TH1D>("#eta y",  ";#eta;Events", NBINS, -12, 12);
	hEta_Neutron = std::make_unique<TH1D>("#eta n",  ";#eta;Events", NBINS, -12, 12);
	hM_NSTAR     = std::make_unique<TH1D>("M (GeV)", ";M (GeV);Events", NBINS, 0, 10);


	// Legendre polynomials, DO NOT CHANGE THE Y-RANGE [-1,1]
	for (std::size_t i = 0; i < 8; ++i) {
		hPl[i] = std::make_unique<TProfile>(Form("hPl%lu", i + 1), "", 100, 0.0, 4.0, -1, 1);
		hPl[i]->Sumw2(); // Error saving on
		hPl[i]->SetXTitle(Form("System M (GeV)"));
		hPl[i]->SetYTitle(Form("#LTP_{l}(cos(#theta)#GT |_{ r.f.}"));
	}

	// Costheta correlations between different frames
	for (std::size_t i = 0; i < NFR; ++i) {
		for (std::size_t j = 0; j < NFR; ++j) {
			h2CosTheta[i][j] = std::make_unique<TH2D>(
			    Form("%s^{+} cos(theta) %s vs %s", pstr.c_str(),
			         frame_labels[i].Data(),
			         frame_labels[j].Data()),
			    Form(
			        ";%s^{+} cos(#theta) %s;%s^{+} cos(#theta) %s",
			        pstr.c_str(), frame_labels[i].Data(),
			        pstr.c_str(), frame_labels[j].Data()),
			    NBINS, -1, 1, NBINS, -1, 1);
		}
	}

	// Phi correlations between different frames
	for (std::size_t i = 0; i < NFR; ++i) {
		for (std::size_t j = 0; j < NFR; ++j) {
			h2Phi[i][j] = std::make_unique<TH2D>(
			    Form("%s^{+} #phi %s vs %s", pstr.c_str(),
			         frame_labels[i].Data(),
			         frame_labels[j].Data()),
			    Form(";%s^{+} #phi %s (rad);%s^{+} #phi %s (rad)",
			         pstr.c_str(), frame_labels[i].Data(),
			         pstr.c_str(), frame_labels[j].Data()),
			    NBINS, -gra::math::PI, gra::math::PI, NBINS, -gra::math::PI,
			    gra::math::PI);
		}
	}

	// 1D costheta
	for (std::size_t i = 0; i < NFR; ++i) {
		hCosTheta_Meson_p[i] = std::make_unique<TH1D>(
		    Form("%s^{+} cos(#theta) %s", pstr.c_str(),
		         frame_labels[i].Data()),
		    Form(";%s^{+} cos(#theta) %s; Events", pstr.c_str(),
		         frame_labels[i].Data()),
		    NBINS, -1, 1);
		hCosTheta_Meson_m[i] = std::make_unique<TH1D>(
		    Form("%s^{-} cos(#theta) %s", pstr.c_str(),
		         frame_labels[i].Data()),
		    Form(";%s^{+} cos(#theta) %s; Events", pstr.c_str(),
		         frame_labels[i].Data()),
		    NBINS, -1, 1);
	}
	
	// 1D phi
	for (std::size_t i = 0; i < NFR; ++i) {
		hPhi_Meson_p[i] = std::make_unique<TH1D>(
		    Form("%s^{+} #phi %s", pstr.c_str(),
		         frame_labels[i].Data()),
		    Form(";%s^{+} #phi %s; Events", pstr.c_str(),
		         frame_labels[i].Data()),
		    NBINS, -gra::math::PI, gra::math::PI);
		hPhi_Meson_m[i] = std::make_unique<TH1D>(
		    Form("%s^{-} #phi %s", pstr.c_str(),
		         frame_labels[i].Data()),
		    Form(";%s^{+} #phi %s; Events", pstr.c_str(),
		         frame_labels[i].Data()),
		    NBINS, -gra::math::PI, gra::math::PI);
	}

	// 2D (costheta, phi)
	for (std::size_t i = 0; i < NFR; ++i) {
		h2CosTheta_Phi[i] = std::make_unique<TH2D>(
		    Form("cos(#theta) vs #phi %s", frame_labels[i].Data()),
		    Form(" ;%s^{+} cos(#theta) %s;%s^{+} #phi %s (rad)",
		         pstr.c_str(), frame_labels[i].Data(), pstr.c_str(),
		         frame_labels[i].Data()),
		    NBINS / 2, -1, 1, NBINS / 2, -gra::math::PI, gra::math::PI);
	}
}

// Destructor
MAnalyzer::~MAnalyzer() {
}


// Fiducial cuts
// return true if event passes, else false
/*
bool MAnalyzer::FiducialCuts(HepMC3::GenEvent& ) {

	return true;
}
*/


// "Oracle" histogram filler:
//
// Oracle here means that in this function we use event tree information (unphysical),
// not just pure fiducial final state information (physical).
//
double MAnalyzer::HepMC3_OracleFill(const std::string input, unsigned int multiplicity, int finalPDG, unsigned int MAXEVENTS,
    std::map<std::string, std::unique_ptr<h1Multiplet>>& h1,
    std::map<std::string, std::unique_ptr<h2Multiplet>>& h2,
    std::map<std::string, std::unique_ptr<hProfMultiplet>>& hP, unsigned int SID) {
	
	inputfile = input;
	const std::string totalpath = gra::aux::GetBasePath(2) + "/output/" + input + ".hepmc3";
	HepMC3::ReaderAscii input_file(totalpath);
	
	if (input_file.failed()) {
		throw std::invalid_argument("MAnalyzer::HepMC3Read: Cannot open file " + totalpath);
	}

	// Event loop
	unsigned int events_read  = 0;

	// Variables for calculating selection efficiency
	double totalW = 0;
	double selecW = 0;

	while (true) {

		// Read event from input file
		HepMC3::GenEvent evt(HepMC3::Units::GEV, HepMC3::Units::MM);
		input_file.read_event(evt);
		
		// Reading failed
		if (input_file.failed()) {
			if (events_read == 0) {
				throw std::invalid_argument("MAnalyzer::HepMC3Read: File " + totalpath + " is empty!");
			} else { break; }
		}
		if (events_read == 0) {
			HepMC3::Print::listing(evt);
			HepMC3::Print::content(evt);
		}
		++events_read;

		// *** Get generator cross section (in picobarns by HepMC3 convention) ***
		std::shared_ptr<HepMC3::GenCrossSection> cs =
		    evt.attribute<HepMC3::GenCrossSection>("GenCrossSection");
		if (cs) {
			cross_section = 1E-12 * cs->xsec(0); // turn into barns
		} else {
			std::cout << "Problem accessing 'GenCrossSection' attribute!" << std::endl;
		}
		// --------------------------------------------------------------
		
		// *** Get event weight (always in barn units) ***
		double W = 1.0;
		if (evt.weights().size() != 0) { // check do we have weights saved
			W = evt.weights()[0];        // take the first one
		}
		totalW += W;
		// --------------------------------------------------------------

		int sign = -1;
		if (finalPDG == PDG::PDG_gamma || finalPDG == PDG::PDG_gluon || finalPDG == 113) {
			sign = 0; // Do not double count neutral
		}
		int NEGfinalPDG = sign * finalPDG;

		// ==============================================================
		sqrts = CheckEnergyMomentum(evt);
		// ==============================================================

		// Central particles
		std::vector<M4Vec> pip;
		std::vector<M4Vec> pim;

		for (HepMC3::ConstGenParticlePtr p1: HepMC3::applyFilter(*abs(HepMC3::Selector::PDG_ID) == finalPDG, evt.particles())) {
			M4Vec pvec = gra::aux::HepMC2M4Vec(p1->momentum());

			// Check that ancestor is a central system
			std::vector<HepMC3::ConstGenParticlePtr> results =
				HepMC3::applyFilter(*abs(HepMC3::Selector::PDG_ID) == PDG::PDG_system, HepMC3::Relatives::ANCESTORS(p1));
			if (results.size() != 0) { pip.push_back(pvec); }
		}
		for (HepMC3::ConstGenParticlePtr p1: HepMC3::applyFilter(*abs(HepMC3::Selector::PDG_ID) == NEGfinalPDG, evt.particles())) {
			M4Vec pvec = gra::aux::HepMC2M4Vec(p1->momentum());

			// Check that ancestor is a central system
			std::vector<HepMC3::ConstGenParticlePtr> results =
				HepMC3::applyFilter(*abs(HepMC3::Selector::PDG_ID) == PDG::PDG_system, HepMC3::Relatives::ANCESTORS(p1));
			if (results.size() != 0) { pim.push_back(pvec); }
		}

		// CHECK CONDITION
		if (pip.size() + pim.size() != (unsigned int)multiplicity) {
			printf("MAnalyzer::ReadHepMC3:: Multiplicity condition not filled %lu %lu %d! \n", pip.size(), pim.size(), multiplicity);
			continue; // skip event
		}

		// ---------------------------------------------------------------
		// CENTRAL SYSTEM plots
		M4Vec system;
		for (const auto& x : pip) { system += x; }
		for (const auto& x : pim) { system += x; }

		std::vector<HepMC3::GenParticlePtr> beam_protons =
			HepMC3::applyFilter(HepMC3::Selector::STATUS == PDG::PDG_BEAM   && *abs(HepMC3::Selector::PDG_ID) == PDG::PDG_p, evt.particles());
		
		std::vector<HepMC3::GenParticlePtr> final_protons =
			HepMC3::applyFilter(HepMC3::Selector::STATUS == PDG::PDG_STABLE && *abs(HepMC3::Selector::PDG_ID) == PDG::PDG_p, evt.particles());

		M4Vec p_beam_plus;
		M4Vec p_beam_minus;
		M4Vec p_final_plus;
		M4Vec p_final_minus;

		// Beam (initial state ) protons
		for (const HepMC3::GenParticlePtr& p1 : beam_protons) {
			M4Vec pvec = gra::aux::HepMC2M4Vec(p1->momentum());
			if (pvec.Rap() > 0) {
				p_beam_plus = pvec;
			} else {
				p_beam_minus = pvec;
			}
		}

		// Final state protons
		for (const HepMC3::GenParticlePtr& p1 : final_protons) {

			M4Vec pvec = gra::aux::HepMC2M4Vec(p1->momentum());

			// Check that ancestor is NOT excited forward system or central system
			std::vector<HepMC3::GenParticlePtr> results =
				HepMC3::applyFilter(*abs(HepMC3::Selector::PDG_ID) == PDG::PDG_NSTAR ||
					                *abs(HepMC3::Selector::PDG_ID) == PDG::PDG_system, HepMC3::Relatives::ANCESTORS(p1));

			if (results.size() == 0) {
				if (pvec.Rap() > 0) {
					p_final_plus = pvec;
				} else {
					p_final_minus = pvec;
				}
			}
		}

		// Observables for 2-body case only
		if (multiplicity == 2) {
			FrameObservables(W, evt, p_beam_plus, p_beam_minus, p_final_plus, p_final_minus, pip, pim);
		}

		// Observables for N stars
		NStarObservables(W, evt);

		// **************************************************************
		// SUPERPLOTTER >>
		try {

			M4Vec a = pip[0];
			M4Vec b;
			if (pim.size() != 0) { // Charged pair
				b = pim[0];
			} else { // Neutral pair
				b = pip[1];
			}

			// Mandelstam -t_1,2
			const double t1 = -(p_beam_plus - p_final_plus).M2();
			//const double t2 = -(p_beam_minus - p_final_minus).M2();

			// Deltaphi
			const double deltaphi_pp = p_final_plus.DeltaPhiAbs(p_final_minus);
			M4Vec pp_diff = p_final_plus - p_final_minus;
			const double pp_dpt = pp_diff.Pt();
			
			// 1D: System
			h1["h1_S_M"]->h[SID]->Fill(system.M(), W);
			h1["h1_S_Pt"]->h[SID]->Fill(system.Pt(), W);
			h1["h1_S_Y"]->h[SID]->Fill(system.Rap(), W);
			hP["hP_S_M_Pt"]->h[SID]->Fill(system.M(), system.Pt(), W);

			// 1D: 1-Body
			h1["h1_1B_pt"]->h[SID]->Fill(a.Pt(), W);
			h1["h1_1B_eta"]->h[SID]->Fill(a.Eta(), W);

			// 1D: Forward proton pair
			h1["h1_PP_dphi"]->h[SID]->Fill(deltaphi_pp, W);
			h1["h1_PP_t1"]->h[SID]->Fill(t1, W);
			h1["h1_PP_dpt"]->h[SID]->Fill(pp_dpt, W);

			// 2D: System
			h2["h2_S_M_Pt"]->h[SID]->Fill(
			    system.M(), system.Pt(), W);

			// 2-Body only
			if (multiplicity == 2) {

				hP["hP_2B_M_dphi"]->h[SID]->Fill(system.M(), a.DeltaPhi(b), W);
				h1["h1_2B_acop"]->h[SID]->Fill(1.0 - a.DeltaPhi(b) / gra::math::PI, W);
				h1["h1_2B_diffrap"]->h[SID]->Fill(b.Rap() - a.Rap(), W);
				
				h2["h2_2B_M_dphi"]->h[SID]->Fill(system.M(), a.DeltaPhi(b), W);
				h2["h2_2B_eta1_eta2"]->h[SID]->Fill(a.Eta(), b.Eta(), W);	
			}

			// 4-Body only 
			if (multiplicity == 4) {
				// ...
			}

		} catch (...) {
			throw std::invalid_argument("MAnalyzer::HepMC3Read: Problem filling histogram!");
		}

		// << SUPERPLOTTER
		// **************************************************************

		if (events_read >= MAXEVENTS){
			std::cout << "MAnalyzer::HepMC3Read: Maximum event count " << MAXEVENTS << " reached!";
			break; // Enough events
		}

		if (events_read % 10000 == 0) {
			std::cout << std::endl
			          << "Events processed: " << events_read
			          << std::endl;
		}
		// [THIS AS LAST!] Sum selected event weights
		selecW += W;
	}
	std::cout << std::endl;
	std::cout << "MAnalyzer::HepMC3Read: Events processed in total: "
			  << events_read << std::endl;

	// Close HepMC3 file
	input_file.close();

	if (selecW == 0.0) {
		throw std::invalid_argument("MAnalyzer::HepMC3Read:: Valid events in <"
			+ totalpath + ">"+ " == 0 out of " + std::to_string(events_read));
	}
	// Take into account extra fiducial cut efficiency here
	double efficiency = selecW / totalW;
	printf("MAnalyzer::HepMC3Read: Fiducial cut efficiency: %0.3f \n", efficiency);
	std::cout << std::endl;

	return cross_section * efficiency;
}


// Sanity check
double MAnalyzer::CheckEnergyMomentum(HepMC3::GenEvent& evt) const {

	std::vector<HepMC3::GenParticlePtr> all_init  =
		HepMC3::applyFilter(HepMC3::Selector::STATUS == PDG::PDG_BEAM, evt.particles());    // Beam
	
	std::vector<HepMC3::GenParticlePtr> all_final =
		HepMC3::applyFilter(HepMC3::Selector::STATUS == PDG::PDG_STABLE, evt.particles());  // Final state
	
	M4Vec beam(0,0,0,0);
	for (const HepMC3::GenParticlePtr& p1 : all_init)  {
		beam += gra::aux::HepMC2M4Vec(p1->momentum());
	}
	M4Vec final(0,0,0,0);
	for (const HepMC3::GenParticlePtr& p1 : all_final) {
		final += gra::aux::HepMC2M4Vec(p1->momentum());
	}
	if (!gra::math::CheckEMC(beam-final)) {
		std::cout << rang::fg::red << 
		"Energy-Momentum not conserved!" << rang::fg::reset << std::endl;
		(beam - final).Print();
		HepMC3::Print::listing(evt);
		HepMC3::Print::content(evt);
	}
	return beam.M();
}


// 2-body angular observables
void MAnalyzer::FrameObservables(double W, HepMC3::GenEvent& evt, const M4Vec& p_beam_plus,
	const M4Vec& p_beam_minus, const M4Vec& p_final_plus, const M4Vec& p_final_minus,
	const std::vector<M4Vec>& pip, const std::vector<M4Vec>& pim) {

	// GJ-frame direction
	const unsigned int direction = 1;
	M4Vec propagator;

	// Calculate propagator (=Pomeron) 4-vector
	if (direction == 1) {
		propagator = p_beam_plus  - p_final_plus;
	} else {
		propagator = p_beam_minus - p_final_minus;
	}
	std::vector<M4Vec> twopions;
	
	if (pip.size() != 0 && pim.size() != 0) { // Charged pair
		twopions = {pip[0], pim[0]};
	}
	if (pip.size() == 2 && pim.size() == 0) { // Neutral pair
		twopions = {pip[0], pip[1]};
	}

	// Make copies
	std::vector<std::vector<M4Vec>> pions;
	for (std::size_t i = 0; i < NFR; ++i) {
		pions.push_back(twopions);
	}

	// Frame transformations
	gra::kinematics::CSframe(pions[0]);
	gra::kinematics::HEframe(pions[1]);
	// LABframe(pions[2]), already there, do nothing
	gra::kinematics::GJframe(pions[3], propagator);
	gra::kinematics::PGframe(pions[4], direction, p_beam_plus, p_beam_minus);
	gra::kinematics::SRframe(pions[5]);
	

	// Legendre polynomials P_l cos(theta)
	const M4Vec system = twopions[0] + twopions[1];
	const unsigned int FRAMENUMBER = 5; // Non-rotated frame
	for (std::size_t l = 0; l < 8; ++l) { // note l+1
		double value = gra::math::LegendrePl((l + 1), pions[FRAMENUMBER][0].CosTheta()); // cos(theta)
		hPl[l]->Fill(system.M(), value, W);
	}

	// FILL HISTOGRAMS -->
	// 2D
	for (std::size_t i = 0; i < NFR; ++i) {
		for (std::size_t j = 0; j < NFR; ++j) {
			h2CosTheta[i][j]->Fill(pions[i][0].CosTheta(), pions[j][0].CosTheta(), W);
			h2Phi[i][j]->Fill(pions[i][0].Phi(), pions[j][0].Phi(), W);
		}
	}
	// 1D, 2D
	for (std::size_t i = 0; i < NFR; ++i) {
		hCosTheta_Meson_p[i]->Fill(pions[i][0].CosTheta(), W);
		hCosTheta_Meson_m[i]->Fill(pions[i][1].CosTheta(), W);

		hPhi_Meson_p[i]->Fill(pions[i][0].Phi(), W);
		hPhi_Meson_m[i]->Fill(pions[i][1].Phi(), W);

		h2CosTheta_Phi[i]->Fill(pions[i][0].CosTheta(), pions[i][0].Phi(), W);
	}	
}

// Forward system observables
void MAnalyzer::NStarObservables(double W, HepMC3::GenEvent& evt) {

	// Excite forward system particles
	std::vector<HepMC3::GenParticlePtr> search_gammas   =
		HepMC3::applyFilter(HepMC3::Selector::STATUS == PDG::PDG_STABLE && *abs(HepMC3::Selector::PDG_ID) == PDG::PDG_gamma, evt.particles());

	std::vector<HepMC3::GenParticlePtr> search_neutrons =
		HepMC3::applyFilter(HepMC3::Selector::STATUS == PDG::PDG_STABLE && *abs(HepMC3::Selector::PDG_ID) == PDG::PDG_n,     evt.particles());

	std::vector<HepMC3::GenParticlePtr> search_pip   =
		HepMC3::applyFilter(HepMC3::Selector::STATUS == PDG::PDG_STABLE && *abs(HepMC3::Selector::PDG_ID) == PDG::PDG_pip,   evt.particles());

	std::vector<HepMC3::GenParticlePtr> search_pim   =
		HepMC3::applyFilter(HepMC3::Selector::STATUS == PDG::PDG_STABLE && *abs(HepMC3::Selector::PDG_ID) == PDG::PDG_pim,   evt.particles());

	std::vector<HepMC3::GenParticlePtr> search_nstar =
		HepMC3::applyFilter(HepMC3::Selector::STATUS == PDG::PDG_STABLE && *abs(HepMC3::Selector::PDG_ID) == PDG::PDG_NSTAR, evt.particles());

	// Find out if we excited one or two protons
	bool excited_plus  = false;
	bool excited_minus = false;
	for (const HepMC3::GenParticlePtr& p1 : search_nstar) {

		M4Vec pvec = gra::aux::HepMC2M4Vec(p1->momentum());
		hM_NSTAR->Fill(pvec.M());

		if (pvec.Rap() > 0) {
			excited_plus = true;
		}
		if (pvec.Rap() < 0) {
			excited_minus = true;
		}
	}
	// Excited system found
	if (excited_plus || excited_minus) {
		N_STAR_ON = true;
	}

	// N* system decay products

	// Gammas
	double gamma_e_plus = 0;
	double gamma_e_minus = 0;

	// Gammas
	for (const HepMC3::GenParticlePtr& p1 : search_gammas) {
		
		M4Vec pvec = gra::aux::HepMC2M4Vec(p1->momentum());

		// Check that ancestor is excited forward system
		std::vector<HepMC3::GenParticlePtr> ancestor =
			HepMC3::applyFilter(*abs(HepMC3::Selector::PDG_ID) == PDG::PDG_NSTAR, HepMC3::Relatives::ANCESTORS(p1));

		if (ancestor.size() != 0) {
			hEta_Gamma->Fill(pvec.Eta(), W);
			hE_Gamma->Fill(pvec.E(), W);
			hXF_Gamma->Fill(pvec.Pz() / (sqrts/2), W);

			if (excited_plus && pvec.Rap() > 0)  {
				gamma_e_plus += pvec.E();
			}
			if (excited_minus && pvec.Rap() < 0) {
				gamma_e_minus += pvec.E();
			}
		}
	}

	// Pi+
	for (const HepMC3::GenParticlePtr& p1 : search_pip) {
		// HepMC3::Print::line(p1);

		// Check that parent is the excited system
		std::vector<HepMC3::GenParticlePtr> parents = p1->parents();
		bool found = false;
		for (const auto& k : indices(parents)) {
			if (parents[k]->pid() == PDG::PDG_NSTAR) { found = true; }
		}
		if (found) {
			M4Vec pvec = gra::aux::HepMC2M4Vec(p1->momentum());
			hEta_Pions->Fill(pvec.Eta(), W);
			hE_Pions->Fill(pvec.E(), W);
			hXF_Pions->Fill(pvec.Pz() / (sqrts/2), W);
		}
	}

	// Pi-
	for (const HepMC3::GenParticlePtr& p1 : search_pim) {
		// HepMC3::Print::line(p1);

		// Check that parent is the excited system
		std::vector<HepMC3::GenParticlePtr> parents = p1->parents();
		bool found = false;
		for (const auto& k : indices(parents)) {
			if (parents[k]->pid() == PDG::PDG_NSTAR) { found = true; }
		}
		if (found) {
			M4Vec pvec = gra::aux::HepMC2M4Vec(p1->momentum());
			hEta_Pions->Fill(pvec.Eta(), W);
			hE_Pions->Fill(pvec.E(), W);
			hXF_Pions->Fill(pvec.Pz() / (sqrts/2), W);
		}
	}
	
	// Neutrons
	double neutron_e_plus  = 0;
	double neutron_e_minus = 0;

	for (const HepMC3::GenParticlePtr& p1 : search_neutrons) {
		M4Vec pvec = gra::aux::HepMC2M4Vec(p1->momentum());

		// Check that parent is the excited system
		std::vector<HepMC3::GenParticlePtr> parents = p1->parents();
		bool found = false;
		for (const auto& k : indices(parents)) {
			if (parents[k]->pid() == PDG::PDG_NSTAR) { found = true; }
		}
		if (found) {
			hEta_Neutron->Fill(pvec.Eta(), W);
			hE_Neutron->Fill(pvec.E(), W);
			hXF_Neutron->Fill(pvec.Pz() / (sqrts/2), W);

			if (excited_plus && pvec.Rap() > 0) {
				neutron_e_plus += pvec.E();
			}
			if (excited_minus && pvec.Rap() < 0) {
				neutron_e_minus += pvec.E();
			}
		}
	}
	
	// Gamma+Neutron energy histogram
	if (excited_plus){
		hE_GammaNeutron->Fill(gamma_e_plus + neutron_e_plus, W);
	}
	if (excited_minus){
		hE_GammaNeutron->Fill(gamma_e_minus + neutron_e_minus, W);
	}
}


double powerlaw(double* x, double* par) {
	return par[0] / std::pow(1.0 + std::pow(x[0], 2) / (std::pow(par[1], 2) * par[2]), par[2]);
}

double exponential(double* x, double* par) {
	return par[0] * exp(par[1] * x[0]);
}


// Plotter
void MAnalyzer::PlotAll() {

    // Create output directory if it does not exist
	const std::string FOLDER = gra::aux::GetBasePath(2) + "/figs/" + inputfile;
   	aux::CreateDirectory(FOLDER);
   	
	// FIT FUNCTIONS
	TF1* fb = new TF1("exp_fit", exponential, 0.05, 0.5, 2);
	fb->SetParameter(0, 10.0); // A
	fb->SetParameter(1, -8.0); // b

	TF1* fa = new TF1("pow_fit", powerlaw, 0.5, 3.0, 3);
	fa->SetParameter(0, 10.0); // A
	fa->SetParameter(1, 0.15); // T
	fa->SetParameter(2, 0.1);  // n

	//        hpT2->Fit("exp_fit","R");
	//        hpT_Meson_p->Fit("pow_fit","R"); // "R" for range

	// *************** FORWARD EXCITATION ***************
	if (N_STAR_ON) {

		// Draw histograms
		TCanvas c1("c", "c", 800, 600);
		c1.Divide(2, 2, 0.0002, 0.0002);

		c1.cd(1);
		//gPad->SetLogy();
		hEta_Gamma->SetLineColor(2);
		hEta_Pions->SetLineColor(4);
		hEta_Neutron->SetLineColor(8);

		hEta_Gamma->Draw("same");
		hEta_Neutron->Draw("same");
		hEta_Pions->Draw("same");
		
		// Set x-axis range
		const double eta_min = 3;
		const double eta_max = 15;
		hEta_Pions->SetAxisRange(eta_min, eta_max, "X");
		hEta_Gamma->SetAxisRange(eta_min, eta_max, "X");
		hEta_Neutron->SetAxisRange(eta_min, eta_max, "X");

		c1.cd(2);
		//gPad->SetLogy();
		hXF_Gamma->SetLineColor(2);
		hXF_Pions->SetLineColor(4);
		hXF_Neutron->SetLineColor(8);

		hXF_Gamma->Draw("same");
		hXF_Pions->Draw("same");
		hXF_Neutron->Draw("same");
		/*
		hE_Gamma->SetLineColor(2);
		hE_Pions->SetLineColor(4);
		hE_Neutron->SetLineColor(8);
		hE_GammaNeutron->SetLineColor(kBlack);

		hE_Gamma->Draw("same");
		hE_Pions->Draw("same");
		hE_Neutron->Draw("same");
		hE_GammaNeutron->Draw("same");
		*/

		c1.cd(3);
		gPad->SetLogy();
		//gPad->SetLogx();
		hM_NSTAR->Draw();

		c1.cd(4);
		gPad->SetLogy();
		hXF_Gamma->SetLineColor(2);
		hXF_Pions->SetLineColor(4);
		hXF_Neutron->SetLineColor(8);

		hXF_Gamma->Draw("same");
		hXF_Pions->Draw("same");
		hXF_Neutron->Draw("same");

		c1.SaveAs(Form("%s/figs/%s/forward.pdf", gra::aux::GetBasePath(2).c_str(), inputfile.c_str()));
	}
	// *************** *************** ***************

	TCanvas c2("c", "c", 800, 800);
	c2.Divide(NFR, NFR, 0.0001, 0.0002);

	int k = 1;
	for (std::size_t i = 0; i < NFR; ++i) {
		for (std::size_t j = 0; j < NFR; ++j) {
			c2.cd(k);
			++k;
			if (j >= i)
				h2CosTheta[i][j]->Draw("COLZ");
		}
	}
	c2.SaveAs(Form("%s/figs/%s/QA_matrix2.pdf", gra::aux::GetBasePath(2).c_str(), inputfile.c_str()));


	TCanvas c3("c", "c", 800, 800);
	c3.Divide(NFR, NFR, 0.0001, 0.0002);

	k = 1;
	for (std::size_t i = 0; i < NFR; ++i) {
		for (std::size_t j = 0; j < NFR; ++j) {
			c3.cd(k);
			++k;
			if (j >= i)
				h2Phi[i][j]->Draw("COLZ");
		}
	}
	c3.SaveAs(Form("%s/figs/%s/QA_matrix3.pdf", gra::aux::GetBasePath(2).c_str(), inputfile.c_str()));


	TCanvas c4("c", "c", 800, 800);
	c4.Divide(NFR, NFR, 0.0001, 0.0002);

	k = 1;
	for (std::size_t i = 0; i < NFR; ++i) {
		c4.cd(k);
		++k;
		// gPad->SetLogy();
		hCosTheta_Meson_p[i]->SetLineColor(9);
		hCosTheta_Meson_m[i]->SetLineColor(46);
		hCosTheta_Meson_p[i]->Draw("same");
		hCosTheta_Meson_m[i]->Draw("same");
	}
	for (std::size_t i = 0; i < NFR; ++i) {
		c4.cd(k);
		++k;
		gPad->SetLogy();
		hPhi_Meson_p[i]->SetLineColor(9);
		hPhi_Meson_m[i]->SetLineColor(46);
		hPhi_Meson_p[i]->Draw("same");
		hPhi_Meson_m[i]->Draw("same");
	}
	for (std::size_t i = 0; i < NFR; ++i) {
		c4.cd(k);
		++k;
		h2CosTheta_Phi[i]->Draw("COLZ");
	}
	c4.SaveAs(Form("%s/figs/%s/QA_matrix4.pdf", gra::aux::GetBasePath(2).c_str(), inputfile.c_str()));

	// -------------------------------------------------------------------------------------
	// -------------------------------------------------------------------------------------
	// Legendre polynomials in the Rest Frame

	const int colors[4] = {48, 53, 98, 32};

	// 1...4
	{
		TCanvas* c115 = new TCanvas("c115", "Legendre polynomials", 600, 400);
		TLegend* leg[4];
		c115->Divide(2, 2, 0.001, 0.001);

		for (std::size_t l = 0; l < 4; ++l) {
			c115->cd(l + 1); // note (l+1)

			leg[l] = new TLegend(0.15, 0.75, 0.4, 0.85); // x1,y1,x2,y2
			hPl[l]->SetLineColor(colors[l]);
			hPl[l]->Draw();
			hPl[l]->SetMinimum(-0.4); // Y-axis minimum
			hPl[l]->SetMaximum(0.4);  // Y-axis maximum

			leg[l]->SetFillColor(0);  // White background
			leg[l]->SetBorderSize(0); // No box
			leg[l]->AddEntry(hPl[l].get(), Form("l = %lu", l + 1),
			                 "l");
			leg[l]->Draw();
		}
		c115->SaveAs(
		    Form("%s/figs/%s/hPl_1to4.pdf", gra::aux::GetBasePath(2).c_str(), inputfile.c_str()));
	}

	// 5...8
	{
		TCanvas* c115 = new TCanvas("c115", "Legendre polynomials", 600, 400);
		TLegend* leg[4];
		c115->Divide(2, 2, 0.001, 0.001);

		const int colors[4] = {48, 53, 98, 32};
		for (std::size_t l = 4; l < 8; ++l) {
			c115->cd(l - 3); // note (l-3)

			leg[l] = new TLegend(0.15, 0.75, 0.4, 0.85); // x1,y1,x2,y2
			hPl[l]->SetLineColor(colors[l - 4]);
			hPl[l]->Draw();
			hPl[l]->SetMinimum(-0.4); // Y-axis minimum
			hPl[l]->SetMaximum(0.4);  // Y-axis maximum

			leg[l]->SetFillColor(0);  // White background
			leg[l]->SetBorderSize(0); // No box
			leg[l]->AddEntry(hPl[l].get(), Form("l = %lu", l + 1),
			                 "l");
			leg[l]->Draw();
		}
		c115->SaveAs(
		    Form("%s/figs/%s/hPl_5to8.pdf", gra::aux::GetBasePath(2).c_str(), inputfile.c_str()));
	}
}

} // gra namespace ends
