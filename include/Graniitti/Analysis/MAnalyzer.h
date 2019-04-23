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

// HepMC33
#include "HepMC3/FourVector.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/Print.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/Relatives.h"
#include "HepMC3/Selector.h"
#include "HepMC3/WriterAscii.h"

// Own
#include "Graniitti/Analysis/MMultiplet.h"
#include "Graniitti/MAux.h"

namespace gra {

class MAnalyzer {
       public:
	// Constructor, destructor
	MAnalyzer();
	~MAnalyzer();

	// Default particle name string
	std::string pstr = "daughter";

	// Different Lorentz frame labels
	const std::vector<TString> frame_labels = {"CS", "HE", "LAB", "GJ", "PG", "SR"};

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

	static constexpr unsigned int NFR = 6; // number of frames

	// Correlations between frames
	std::unique_ptr<TH2D> h2CosTheta[NFR][NFR];
	std::unique_ptr<TH2D> h2Phi[NFR][NFR];

	std::unique_ptr<TH1D> hCosTheta_Meson_p[NFR];
	std::unique_ptr<TH1D> hCosTheta_Meson_m[NFR];
	std::unique_ptr<TH1D> hPhi_Meson_p[NFR];
	std::unique_ptr<TH1D> hPhi_Meson_m[NFR];
	std::unique_ptr<TH2D> h2CosTheta_Phi[NFR];

	std::unique_ptr<TH2D> h2M_CosTheta[NFR];
	std::unique_ptr<TH2D> h2M_Phi[NFR];

	// HepMC3 reader
	double HepMC3_OracleFill(const std::string inputfile, unsigned int multiplicity,
	                         int finalPDG, unsigned int MAXEVENTS,
	                         std::map<std::string, std::unique_ptr<h1Multiplet>> &h1,
	                         std::map<std::string, std::unique_ptr<h2Multiplet>> &h2,
	                         std::map<std::string, std::unique_ptr<hProfMultiplet>> &hP,
	                         unsigned int SID);

	// Plot out all local histograms
	void PlotAll(const std::string &titlestr);

	double cross_section = 0;

	double CheckEnergyMomentum(HepMC3::GenEvent &evt) const;
	void FrameObservables(double W, HepMC3::GenEvent &evt, const M4Vec &p_beam_plus,
	                      const M4Vec &p_beam_minus, const M4Vec &p_final_plus,
	                      const M4Vec &p_final_minus, const std::vector<M4Vec> &pip,
	                      const std::vector<M4Vec> &pim);
	void NStarObservables(double W, HepMC3::GenEvent &evt);

       private:
	double sqrts = 0.0;

	// Name of the HepMC33 input
	std::string inputfile;

	// Proton excitation is turned on
	bool N_STAR_ON = false;
};

} // gra namespace ends

#endif