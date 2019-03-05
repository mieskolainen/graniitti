// Fast analysis class
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MANALYZER_H
#define MANALYZER_H

#include <complex>
#include <memory>
#include <vector>


// ROOT
#include "TBranch.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TRandom3.h"
#include "TString.h"
#include "TStyle.h"
#include "TTree.h"


// HepMC3
#include "HepMC/FourVector.h"
#include "HepMC/GenEvent.h"
#include "HepMC/GenParticle.h"
#include "HepMC/GenVertex.h"
#include "HepMC/Print.h"
#include "HepMC/ReaderAscii.h"
#include "HepMC/Search/FindParticles.h"
#include "HepMC/WriterAscii.h"


// Own
#include "Graniitti/MAux.h"
#include "Graniitti/Analysis/MMultiplet.h"


namespace gra {


class MAnalyzer {

public:
	// Constructor, destructor
	MAnalyzer();
	~MAnalyzer();

	// Default particle name string
	std::string pstr = "daughter";

	// Different Lorentz frame labels
	const std::vector<TString> frame_labels = {"CS", "HE",  "LAB",
	                                           "GJ", "PG", "SR"};

	// ----------------------------------------------------------
	// Forward system quantities
	std::unique_ptr<TH1D> hE_Pions;
	std::unique_ptr<TH1D> hE_Gamma;
	std::unique_ptr<TH1D> hE_Neutron;
	std::unique_ptr<TH1D> hE_GammaNeutron;

	std::unique_ptr<TH1D> hXF_Pions;
	std::unique_ptr<TH1D> hXF_Gamma;
	std::unique_ptr<TH1D> hXF_Neutron;

	std::unique_ptr<TH1D> hEta_Pions;
	std::unique_ptr<TH1D> hEta_Gamma;
	std::unique_ptr<TH1D> hEta_Neutron;
	std::unique_ptr<TH1D> hM_NSTAR;


	// ----------------------------------------------------------
	// Angular observables
	std::unique_ptr<TProfile> hPl[8];

	const uint NFR = 6; // number of frames

	std::unique_ptr<TH2D> h2CosTheta[6][6];
	std::unique_ptr<TH2D> h2Phi[6][6];
	std::unique_ptr<TH1D> hCosTheta_Meson_p[6];
	std::unique_ptr<TH1D> hCosTheta_Meson_m[6];
	std::unique_ptr<TH1D> hPhi_Meson_p[6];
	std::unique_ptr<TH1D> hPhi_Meson_m[6];
	std::unique_ptr<TH2D> h2CosTheta_Phi[6];
	
	
	// HepMC3 reader
	double HepMC3_OracleFill(
	    const std::string inputfile, uint multiplicity, int finalPDG, uint MAXEVENTS,
	    std::map<std::string, std::unique_ptr<h1Multiplet> >& h1,
	    std::map<std::string, std::unique_ptr<h2Multiplet> >& h2,
	    std::map<std::string, std::unique_ptr<hProfMultiplet> >& hP,
	    uint SID);
	
	
	// Plot out all local histograms
	void PlotAll();

	double cross_section = 0;

	double CheckEnergyMomentum(HepMC::GenEvent& evt) const;
	void FrameObservables(double W, HepMC::GenEvent& evt, const M4Vec& p_beam_plus,
							const M4Vec& p_beam_minus, const M4Vec& p_final_plus, const M4Vec& p_final_minus,
							const std::vector<M4Vec>& pip, const std::vector<M4Vec>& pim);
	void NStarObservables(double W, HepMC::GenEvent& evt);
	
private:
	
	double sqrts = 0.0;
	
	// Name of the HepMC3 input
	std::string inputfile;
	
	// Proton excitation is turned on
	bool N_STAR_ON = false;
};


} // gra namespace ends

#endif