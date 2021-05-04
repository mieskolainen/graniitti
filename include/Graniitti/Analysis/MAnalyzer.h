// Fast analysis class
//
// (c) 2017-2021 Mikael Mieskolainen
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

namespace analyzer {

// Different Lorentz frame labels
const std::vector<std::string> FRAMES = {"CM", "HX", "CS", "PG", "GJ", "LAB"};

}  // namespace analyzer

class MAnalyzer {
 public:
  // Constructor, destructor
  MAnalyzer(const std::string &ID);
  ~MAnalyzer();

  // Default daughter particle name string
  std::string pstr = "daughter";

  // ----------------------------------------------------------
  // Forward system quantities
  std::shared_ptr<TH1D> hE_Pions;
  std::shared_ptr<TH1D> hE_Gamma;
  std::shared_ptr<TH1D> hE_Neutron;
  std::shared_ptr<TH1D> hE_GammaNeutron;

  std::shared_ptr<TH1D> hXF_Pions;
  std::shared_ptr<TH1D> hXF_Gamma;
  std::shared_ptr<TH1D> hXF_Neutron;

  std::shared_ptr<TH1D> hEta_Pions;
  std::shared_ptr<TH1D> hEta_Gamma;
  std::shared_ptr<TH1D> hEta_Neutron;
  std::shared_ptr<TH1D> hM_NSTAR;

  // ----------------------------------------------------------
  // Angular observables
  std::shared_ptr<TProfile> hPl[8];

  static constexpr unsigned int NFR = 6;  // number of frames

  // Correlations between frames
  std::shared_ptr<TH2D> h2CosTheta[NFR][NFR];
  std::shared_ptr<TH2D> h2Phi[NFR][NFR];
  // ----------------------------------------------------------

  // HepMC3 reader
  double HepMC3_OracleFill(const std::string inputfile, unsigned int multiplicity, int finalPDG,
                           unsigned int                                            MAXEVENTS,
                           std::map<std::string, std::shared_ptr<h1Multiplet>> &   h1,
                           std::map<std::string, std::shared_ptr<h2Multiplet>> &   h2,
                           std::map<std::string, std::shared_ptr<hProfMultiplet>> &hP,
                           unsigned int                                            SID);

  // Plot out all local histograms
  void PlotAll(const std::string &titlestr);

  double cross_section = 0;

  double CheckEnergyMomentum(HepMC3::GenEvent &evt) const;
  void   FrameObservables(double W, HepMC3::GenEvent &evt, const M4Vec &p_beam_plus,
                          const M4Vec &p_beam_minus, const M4Vec &p_final_plus,
                          const M4Vec &p_final_minus, const std::vector<M4Vec> &pip,
                          const std::vector<M4Vec> &pim);
  void   NStarObservables(double W, HepMC3::GenEvent &evt);

 private:
  double sqrts = 0.0;

  // Name of the HepMC33 input
  std::string inputfile;

  // Proton excitation is turned on
  bool N_STAR_ON = false;
};

}  // namespace gra

#endif