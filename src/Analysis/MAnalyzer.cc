// Fast MC analysis class
//
// (c) 2017-2020 Mikael Mieskolainen
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
#include "TLine.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TRandom3.h"
#include "TString.h"
#include "TStyle.h"
#include "TTree.h"

// HepMC3 3
#include "HepMC3/FourVector.h"
#include "HepMC3/GenCrossSection.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/Print.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/Relatives.h"
#include "HepMC3/Selector.h"

// Own
#include "Graniitti/Analysis/MAnalyzer.h"
#include "Graniitti/Analysis/MMultiplet.h"
#include "Graniitti/M4Vec.h"
#include "Graniitti/MKinematics.h"
#include "Graniitti/MMath.h"
#include "Graniitti/MPDG.h"

const bool DEBUG = false;

using gra::aux::indices;
using gra::math::msqrt;

namespace gra {

// Constructor with unique ID string for ROOT bookkeeping reasons
MAnalyzer::MAnalyzer(const std::string &ID) {
  // Initialize histograms
  const int NBINS = 150;

  // Energy
  hE_Pions        = std::make_shared<TH1D>(Form("%s_%s", "Energy #pi (GeV)", ID.c_str()),
                                    ";Energy (GeV);Events", NBINS, 0, sqrts / 2.0);
  hE_Gamma        = std::make_shared<TH1D>(Form("%s_%s", "Energy #gamma (GeV)", ID.c_str()),
                                    ";Energy (GeV);Events", NBINS, 0, sqrts / 2.0);
  hE_Neutron      = std::make_shared<TH1D>(Form("%s_%s", "Energy n (GeV)", ID.c_str()),
                                      ";Energy (GeV);Events", NBINS, 0, sqrts / 2.0);
  hE_GammaNeutron = std::make_shared<TH1D>(Form("%s_%s", "Energy y+n (GeV)", ID.c_str()),
                                           ";Energy (GeV);Events", NBINS, 0, sqrts / 2.0);

  // Feynman-x
  hXF_Pions   = std::make_shared<TH1D>(Form("%s_%s", "xF #pi", ID.c_str()), ";Feynman-x;Events",
                                     NBINS, -1.0, 1.0);
  hXF_Gamma   = std::make_shared<TH1D>(Form("%s_%s", "xF #gamma", ID.c_str()), ";Feynman-x;Events",
                                     NBINS, -1.0, 1.0);
  hXF_Neutron = std::make_shared<TH1D>(Form("%s_%s", "xF n", ID.c_str()), ";Feynman-x;Events",
                                       NBINS, -1.0, 1.0);

  // Forward systems
  hEta_Pions =
      std::make_shared<TH1D>(Form("%s_%s", "#eta pi", ID.c_str()), ";#eta;Events", NBINS, -12, 12);
  hEta_Gamma =
      std::make_shared<TH1D>(Form("%s_%s", "#eta y", ID.c_str()), ";#eta;Events", NBINS, -12, 12);
  hEta_Neutron =
      std::make_shared<TH1D>(Form("%s_%s", "#eta n", ID.c_str()), ";#eta;Events", NBINS, -12, 12);
  hM_NSTAR =
      std::make_shared<TH1D>(Form("%s_%s", "M (GeV)", ID.c_str()), ";M (GeV);Events", NBINS, 0, 10);

  // Legendre polynomials, DO NOT CHANGE THE Y-RANGE [-1,1]
  for (std::size_t i = 0; i < 8; ++i) {
    hPl[i] =
        std::make_shared<TProfile>(Form("hPl%lu_%s", i + 1, ID.c_str()), "", 100, 0.0, 4.0, -1, 1);
    hPl[i]->Sumw2();  // Error saving on
    hPl[i]->SetXTitle(Form("System M (GeV)"));
    hPl[i]->SetYTitle(Form("Legendre #LTP_{l}(cos #theta)#GT [CM frame]"));
  }

  // Costheta correlations between different frames
  for (std::size_t i = 0; i < analyzer::FRAMES.size(); ++i) {
    for (std::size_t j = 0; j < analyzer::FRAMES.size(); ++j) {
      h2CosTheta[i][j] = std::make_shared<TH2D>(
          Form("%s^{+} cos(theta) %s vs %s_%s", pstr.c_str(), analyzer::FRAMES[i].c_str(),
               analyzer::FRAMES[j].c_str(), ID.c_str()),
          Form(";%s^{+} cos(#theta) %s;%s^{+} cos(#theta) %s", pstr.c_str(),
               analyzer::FRAMES[i].c_str(), pstr.c_str(), analyzer::FRAMES[j].c_str()),
          NBINS, -1, 1, NBINS, -1, 1);
    }
  }

  // Phi correlations between different frames
  for (std::size_t i = 0; i < analyzer::FRAMES.size(); ++i) {
    for (std::size_t j = 0; j < analyzer::FRAMES.size(); ++j) {
      h2Phi[i][j] = std::make_shared<TH2D>(
          Form("%s^{+} #phi %s vs %s_%s", pstr.c_str(), analyzer::FRAMES[i].c_str(),
               analyzer::FRAMES[j].c_str(), ID.c_str()),
          Form(";%s^{+} #phi %s (rad);%s^{+} #phi %s (rad)", pstr.c_str(),
               analyzer::FRAMES[i].c_str(), pstr.c_str(), analyzer::FRAMES[j].c_str()),
          NBINS, -gra::math::PI, gra::math::PI, NBINS, -gra::math::PI, gra::math::PI);
    }
  }
}

// Destructor
MAnalyzer::~MAnalyzer() {}

// "Oracle" histogram filler:
//
// Oracle here means that in this function we (may) use event tree information,
// not just pure fiducial final state information based on purely physical observables.
//
double MAnalyzer::HepMC3_OracleFill(const std::string input, unsigned int multiplicity,
                                    int finalPDG, unsigned int MAXEVENTS,
                                    std::map<std::string, std::shared_ptr<h1Multiplet>> &   h1,
                                    std::map<std::string, std::shared_ptr<h2Multiplet>> &   h2,
                                    std::map<std::string, std::shared_ptr<hProfMultiplet>> &hP,
                                    unsigned int                                            SID) {
  inputfile                     = input;
  const std::string   totalpath = gra::aux::GetBasePath(2) + "/output/" + input + ".hepmc3";
  HepMC3::ReaderAscii input_file(totalpath);

  if (input_file.failed()) {
    throw std::invalid_argument("MAnalyzer::HepMC3Read: Cannot open file " + totalpath);
  }

  // Event loop
  unsigned int events_read = 0;

  // Variables for calculating selection efficiency
  double totalW = 0;
  double selecW = 0;

  // ---------------------------------------------------------------------
  // Set final state [charged pair or neutral pair]
  MPDG PDG;
  PDG.ReadParticleData(gra::aux::GetBasePath(2) + "/modeldata/mass_width_2018.mcd");

  // Try to find the particle from PDG table, will throw exception if fails
  MParticle p           = PDG.FindByPDG(finalPDG);
  const int NEGfinalPDG = (p.chargeX3 != 0) ? -finalPDG : 0;  // Do not double count neutral
  // ---------------------------------------------------------------------

  while (true) {
    // Read event from input file
    HepMC3::GenEvent evt(HepMC3::Units::GEV, HepMC3::Units::MM);
    input_file.read_event(evt);

    // Reading failed
    if (input_file.failed()) {
      if (events_read == 0) {
        throw std::invalid_argument("MAnalyzer::HepMC3Read: File " + totalpath + " is empty!");
      } else {
        break;
      }
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
      cross_section = 1E-12 * cs->xsec(0);  // turn into barns
    } else {
      std::cout << "Problem accessing 'GenCrossSection' attribute!" << std::endl;
    }
    // --------------------------------------------------------------

    // *** Get event weight (always in barn units) ***
    double W = 1.0;
    if (evt.weights().size() != 0) {  // check do we have weights saved
      W = evt.weights()[0];           // take the first one
    }
    totalW += W;
    // --------------------------------------------------------------

    // Central particles
    std::vector<M4Vec> pip;
    std::vector<M4Vec> pim;

    for (HepMC3::ConstGenParticlePtr p1 :
         HepMC3::applyFilter(HepMC3::StandardSelector::PDG_ID == finalPDG, evt.particles())) {
      M4Vec pvec = gra::aux::HepMC2M4Vec(p1->momentum());

      // Check that ancestor is a central system
      std::vector<HepMC3::ConstGenParticlePtr> results =
          HepMC3::applyFilter(*abs(HepMC3::StandardSelector::PDG_ID) == PDG::PDG_system,
                              HepMC3::Relatives::ANCESTORS(p1));
      if (results.size() != 0) { pip.push_back(pvec); }
    }
    for (HepMC3::ConstGenParticlePtr p1 :
         HepMC3::applyFilter(HepMC3::StandardSelector::PDG_ID == NEGfinalPDG, evt.particles())) {
      M4Vec pvec = gra::aux::HepMC2M4Vec(p1->momentum());

      // Check that ancestor is a central system
      std::vector<HepMC3::ConstGenParticlePtr> results = HepMC3::applyFilter(
          HepMC3::StandardSelector::PDG_ID == PDG::PDG_system, HepMC3::Relatives::ANCESTORS(p1));
      if (results.size() != 0) { pim.push_back(pvec); }
    }

    // CHECK CONDITION
    if (pip.size() + pim.size() != (unsigned int)multiplicity) {
      printf(
          "MAnalyzer::ReadHepMC3:: Multiplicity condition not filled +[%lu] "
          "-[%lu] %d! \n",
          pip.size(), pim.size(), multiplicity);
      continue;  // skip event
    }

    // ---------------------------------------------------------------
    // CENTRAL SYSTEM plots
    M4Vec system;
    for (const auto &x : pip) { system += x; }
    for (const auto &x : pim) { system += x; }

    std::vector<HepMC3::GenParticlePtr> beam_protons =
        HepMC3::applyFilter(HepMC3::StandardSelector::STATUS == PDG::PDG_BEAM &&
                                HepMC3::StandardSelector::PDG_ID == PDG::PDG_p,
                            evt.particles());

    std::vector<HepMC3::GenParticlePtr> final_protons =
        HepMC3::applyFilter(HepMC3::StandardSelector::STATUS == PDG::PDG_STABLE &&
                                HepMC3::StandardSelector::PDG_ID == PDG::PDG_p,
                            evt.particles());

    M4Vec p_beam_plus;
    M4Vec p_beam_minus;
    M4Vec p_final_plus;
    M4Vec p_final_minus;

    // If we have full event, check energy-momentum
    if (beam_protons.size()) {
      // ==============================================================
      sqrts = CheckEnergyMomentum(evt);
      // ==============================================================
    }

    // Beam (initial state ) protons
    for (const HepMC3::GenParticlePtr &p1 : beam_protons) {
      M4Vec pvec = gra::aux::HepMC2M4Vec(p1->momentum());
      if (pvec.Rap() > 0) {
        p_beam_plus = pvec;
      } else {
        p_beam_minus = pvec;
      }
    }

    // Final state protons
    for (const HepMC3::GenParticlePtr &p1 : final_protons) {
      M4Vec pvec = gra::aux::HepMC2M4Vec(p1->momentum());

      // Check that ancestor is NOT excited forward system or central system
      std::vector<HepMC3::GenParticlePtr> results =
          HepMC3::applyFilter(HepMC3::StandardSelector::PDG_ID == PDG::PDG_NSTAR ||
                                  HepMC3::StandardSelector::PDG_ID == -PDG::PDG_NSTAR ||
                                  HepMC3::StandardSelector::PDG_ID == PDG::PDG_system,
                              HepMC3::Relatives::ANCESTORS(p1));

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
      if (pim.size() != 0) {  // Charged pair
        b = pim[0];
      } else {  // Neutral pair
        b = pip[1];
      }

      const double M  = system.M();
      const double Pt = system.Pt();
      const double Y  = system.Rap();

      // 1D: System
      h1["h1_S_M"]->h[SID]->Fill(M, W);
      h1["h1_S_Pt"]->h[SID]->Fill(Pt, W);
      h1["h1_S_Pt2"]->h[SID]->Fill(math::pow2(Pt), W);
      h1["h1_S_Y"]->h[SID]->Fill(Y, W);
      hP["hP_S_M_Pt"]->h[SID]->Fill(M, Pt, W);

      // 1D: 1-Body
      h1["h1_1B_pt"]->h[SID]->Fill(a.Pt(), W);
      h1["h1_1B_eta"]->h[SID]->Fill(a.Eta(), W);

      // 1D: Forward proton pair
      double deltaphi_pp = -1.0;
      if (p_final_plus.M() > 0) {
        // Mandelstam -t_1,2
        const double t1 = -(p_beam_plus - p_final_plus).M2();
        // const double t2 = -(p_beam_minus - p_final_minus).M2();

        // Deltaphi
        deltaphi_pp          = p_final_plus.DeltaPhiAbs(p_final_minus);
        M4Vec        pp_diff = p_final_plus - p_final_minus;
        const double pp_dpt  = pp_diff.Pt();

        h1["h1_PP_dphi"]->h[SID]->Fill(deltaphi_pp, W);
        h1["h1_PP_t1"]->h[SID]->Fill(t1, W);
        h1["h1_PP_dpt"]->h[SID]->Fill(pp_dpt, W);

        h2["h2_S_M_dphipp"]->h[SID]->Fill(M, deltaphi_pp, W);
        h2["h2_S_M_dpt"]->h[SID]->Fill(M, pp_dpt, W);
        h2["h2_S_M_t"]->h[SID]->Fill(M, std::abs(t1), W);
      }

      // 2D
      h2["h2_S_M_Pt"]->h[SID]->Fill(M, Pt, W);
      h2["h2_S_M_pt"]->h[SID]->Fill(M, a.Pt(), W);

      // 2-Body only
      if (multiplicity == 2) {
        hP["hP_2B_M_dphi"]->h[SID]->Fill(M, a.DeltaPhi(b), W);
        h1["h1_2B_acop"]->h[SID]->Fill(1.0 - a.DeltaPhi(b) / gra::math::PI, W);
        h1["h1_2B_diffrap"]->h[SID]->Fill(b.Rap() - a.Rap(), W);
        h2["h2_2B_M_dphi"]->h[SID]->Fill(M, a.DeltaPhi(b), W);
        h2["h2_2B_eta1_eta2"]->h[SID]->Fill(a.Eta(), b.Eta(), W);


        // Frame transform
        const int   direction = 1;
        const M4Vec X         = a + b;

        std::vector<M4Vec> CM = {a, b};
        gra::kinematics::CMframe(CM, X);

        std::vector<M4Vec> HX = {a, b};
        gra::kinematics::HXframe(HX, X);

        std::vector<M4Vec> CS = {a, b};
        gra::kinematics::CSframe(CS, X, p_beam_plus, p_beam_minus);

        std::vector<M4Vec> GJ = {a, b};
        gra::kinematics::GJframe(GJ, X, direction, p_beam_plus - p_final_plus,
                                 p_beam_minus - p_final_minus);

        std::vector<M4Vec> PG = {a, b};
        gra::kinematics::PGframe(PG, X, direction, p_beam_plus, p_beam_minus);


        h1["h1_costheta_CM"]->h[SID]->Fill(CM[0].CosTheta(), W);
        h1["h1_costheta_HX"]->h[SID]->Fill(HX[0].CosTheta(), W);
        h1["h1_costheta_CS"]->h[SID]->Fill(CS[0].CosTheta(), W);
        h1["h1_costheta_GJ"]->h[SID]->Fill(GJ[0].CosTheta(), W);
        h1["h1_costheta_PG"]->h[SID]->Fill(PG[0].CosTheta(), W);
        h1["h1_costheta_LAB"]->h[SID]->Fill(a.CosTheta(), W);


        h1["h1_phi_CM"]->h[SID]->Fill(CM[0].Phi(), W);
        h1["h1_phi_HX"]->h[SID]->Fill(HX[0].Phi(), W);
        h1["h1_phi_CS"]->h[SID]->Fill(CS[0].Phi(), W);
        h1["h1_phi_GJ"]->h[SID]->Fill(GJ[0].Phi(), W);
        h1["h1_phi_PG"]->h[SID]->Fill(PG[0].Phi(), W);
        h1["h1_phi_LAB"]->h[SID]->Fill(a.Phi(), W);


        h2["h2_2B_costheta_phi_CM"]->h[SID]->Fill(CM[0].CosTheta(), CM[0].Phi(), W);
        h2["h2_2B_costheta_phi_HX"]->h[SID]->Fill(HX[0].CosTheta(), HX[0].Phi(), W);
        h2["h2_2B_costheta_phi_CS"]->h[SID]->Fill(CS[0].CosTheta(), CS[0].Phi(), W);
        h2["h2_2B_costheta_phi_GJ"]->h[SID]->Fill(GJ[0].CosTheta(), GJ[0].Phi(), W);
        h2["h2_2B_costheta_phi_PG"]->h[SID]->Fill(PG[0].CosTheta(), PG[0].Phi(), W);
        h2["h2_2B_costheta_phi_LAB"]->h[SID]->Fill(a.CosTheta(), a.Phi(), W);


        h2["h2_2B_M_costheta_CM"]->h[SID]->Fill(M, CM[0].CosTheta(), W);
        h2["h2_2B_M_costheta_HX"]->h[SID]->Fill(M, HX[0].CosTheta(), W);
        h2["h2_2B_M_costheta_CS"]->h[SID]->Fill(M, CS[0].CosTheta(), W);
        h2["h2_2B_M_costheta_GJ"]->h[SID]->Fill(M, GJ[0].CosTheta(), W);
        h2["h2_2B_M_costheta_PG"]->h[SID]->Fill(M, PG[0].CosTheta(), W);
        h2["h2_2B_M_costheta_LAB"]->h[SID]->Fill(M, a.CosTheta(), W);


        h2["h2_2B_M_phi_CM"]->h[SID]->Fill(M, CM[0].Phi(), W);
        h2["h2_2B_M_phi_HX"]->h[SID]->Fill(M, HX[0].Phi(), W);
        h2["h2_2B_M_phi_CS"]->h[SID]->Fill(M, CS[0].Phi(), W);
        h2["h2_2B_M_phi_GJ"]->h[SID]->Fill(M, GJ[0].Phi(), W);
        h2["h2_2B_M_phi_PG"]->h[SID]->Fill(M, PG[0].Phi(), W);
        h2["h2_2B_M_phi_LAB"]->h[SID]->Fill(M, a.Phi(), W);


        // ---------------------------------------------------------------------------
        hP["hP_S_M_PL2_CM"]->h[SID]->Fill(M, math::LegendrePl(2, CM[0].CosTheta()), W);
        hP["hP_S_M_PL4_CM"]->h[SID]->Fill(M, math::LegendrePl(4, CM[0].CosTheta()), W);
        h2["h2_2B_eta1_eta2"]->h[SID]->Fill(a.Eta(), b.Eta(), W);
        // ---------------------------------------------------------------------------
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

    if (events_read >= MAXEVENTS) {
      std::cout << "MAnalyzer::HepMC3Read: Maximum event count " << MAXEVENTS << " reached!";
      break;  // Enough events
    }

    if (events_read % 10000 == 0) {
      std::cout << std::endl << "Events processed: " << events_read << std::endl;
    }
    // [THIS AS LAST!] Sum selected event weights
    selecW += W;
  }
  std::cout << std::endl;
  std::cout << "MAnalyzer::HepMC3Read: Events processed in total: " << events_read << std::endl;

  // Close HepMC3 file
  input_file.close();

  if (selecW == 0.0) {
    throw std::invalid_argument("MAnalyzer::HepMC3Read:: Valid events in <" + totalpath + ">" +
                                " == 0 out of " + std::to_string(events_read));
  }
  // Take into account extra fiducial cut efficiency here
  double efficiency = selecW / totalW;
  printf("MAnalyzer::HepMC3Read: Fiducial cut efficiency: %0.3f \n", efficiency);
  std::cout << std::endl;

  return cross_section * efficiency;
}

// Sanity check
double MAnalyzer::CheckEnergyMomentum(HepMC3::GenEvent &evt) const {
  std::vector<HepMC3::GenParticlePtr> all_init = HepMC3::applyFilter(
      HepMC3::StandardSelector::STATUS == PDG::PDG_BEAM, evt.particles());  // Beam

  std::vector<HepMC3::GenParticlePtr> all_final = HepMC3::applyFilter(
      HepMC3::StandardSelector::STATUS == PDG::PDG_STABLE, evt.particles());  // Final state

  M4Vec beam(0, 0, 0, 0);
  for (const HepMC3::GenParticlePtr &p1 : all_init) {
    beam += gra::aux::HepMC2M4Vec(p1->momentum());
  }
  M4Vec final(0, 0, 0, 0);
  for (const HepMC3::GenParticlePtr &p1 : all_final) {
    final += gra::aux::HepMC2M4Vec(p1->momentum());
  }
  if (!gra::math::CheckEMC(beam - final)) {
    gra::aux::PrintWarning();
    std::cout << rang::fg::red << "Energy-Momentum not conserved!" << rang::fg::reset << std::endl;
    (beam - final).Print();
    HepMC3::Print::listing(evt);
    HepMC3::Print::content(evt);
  }
  return beam.M();
}

// 2-body angular observables
void MAnalyzer::FrameObservables(double W, HepMC3::GenEvent &evt, const M4Vec &p_beam_plus,
                                 const M4Vec &p_beam_minus, const M4Vec &p_final_plus,
                                 const M4Vec &p_final_minus, const std::vector<M4Vec> &pip,
                                 const std::vector<M4Vec> &pim) {
  // Find index
  const auto ind = [&](const std::string str) {
    for (const auto &i : indices(analyzer::FRAMES)) {
      if (analyzer::FRAMES[i] == str) { return i; }
    }
    throw std::invalid_argument("MAnalyzer::FrameObservables: Unknown Lorentz frame: " + str);
  };

  std::vector<M4Vec> pf;

  if (pip.size() != 0 && pim.size() != 0) {  // Charged pair
    pf = {pip[0], pim[0]};
  }
  if (pip.size() == 2 && pim.size() == 0) {  // Neutral pair
    pf = {pip[0], pip[1]};
  }

  // ---------------------------------------------------------------------
  // Lorentz frame transformations

  // Make copies
  std::vector<std::vector<M4Vec>> pions;
  for (std::size_t i = 0; i < analyzer::FRAMES.size(); ++i) { pions.push_back(pf); }

  // System
  const M4Vec X         = pf[0] + pf[1];
  const int   direction = 1;  // PG and GJ

  gra::kinematics::CMframe(pions[ind("CM")], X);
  gra::kinematics::HXframe(pions[ind("HX")], X);
  gra::kinematics::CSframe(pions[ind("CS")], X, p_beam_plus, p_beam_minus);
  gra::kinematics::GJframe(pions[ind("GJ")], X, direction, p_beam_plus - p_final_plus,
                           p_beam_minus - p_final_minus);
  gra::kinematics::PGframe(pions[ind("PG")], X, direction, p_beam_plus, p_beam_minus);
  pions[ind("LAB")] = pions[ind("LAB")];  // already there, do nothing

  // No forward protons -> set zero
  if (p_final_plus.M() < 0.5) { pions[ind("GJ")] = {M4Vec(0, 0, 0, 0), M4Vec(0, 0, 0, 0)}; }

  // ---------------------------------------------------------------------

  // FILL HISTOGRAMS -->

  // Legendre polynomials P_l cos(theta)
  for (std::size_t l = 0; l < 8; ++l) {  // note l+1
    // Take first daughter [0]
    double value = gra::math::LegendrePl((l + 1), pions[ind("CM")][0].CosTheta());
    hPl[l]->Fill(X.M(), value, W);
  }

  // FRAME correlations
  for (std::size_t i = 0; i < analyzer::FRAMES.size(); ++i) {
    for (std::size_t j = 0; j < analyzer::FRAMES.size(); ++j) {
      h2CosTheta[i][j]->Fill(pions[i][0].CosTheta(), pions[j][0].CosTheta(), W);
      h2Phi[i][j]->Fill(pions[i][0].Phi(), pions[j][0].Phi(), W);
    }
  }
}

// Forward system observables
void MAnalyzer::NStarObservables(double W, HepMC3::GenEvent &evt) {
  // Excite forward system particles
  std::vector<HepMC3::GenParticlePtr> search_gammas =
      HepMC3::applyFilter(HepMC3::StandardSelector::PDG_ID == PDG::PDG_gamma, evt.particles());

  std::vector<HepMC3::GenParticlePtr> search_neutrons =
      HepMC3::applyFilter(HepMC3::StandardSelector::PDG_ID == PDG::PDG_n, evt.particles());

  std::vector<HepMC3::GenParticlePtr> search_pip =
      HepMC3::applyFilter(HepMC3::StandardSelector::PDG_ID == PDG::PDG_pip, evt.particles());

  std::vector<HepMC3::GenParticlePtr> search_pim =
      HepMC3::applyFilter(HepMC3::StandardSelector::PDG_ID == PDG::PDG_pim, evt.particles());

  std::vector<HepMC3::GenParticlePtr> search_nstar =
      HepMC3::applyFilter(HepMC3::StandardSelector::PDG_ID == PDG::PDG_NSTAR ||
                              HepMC3::StandardSelector::PDG_ID == -PDG::PDG_NSTAR,
                          evt.particles());

  // Find out if we excited one or two protons
  bool excited_plus  = false;
  bool excited_minus = false;
  for (const HepMC3::GenParticlePtr &p1 : search_nstar) {
    M4Vec pvec = gra::aux::HepMC2M4Vec(p1->momentum());
    hM_NSTAR->Fill(pvec.M(), W);

    if (pvec.Rap() > 0) { excited_plus = true; }
    if (pvec.Rap() < 0) { excited_minus = true; }
  }
  // Excited system found
  if (excited_plus || excited_minus) { N_STAR_ON = true; }

  // N* system decay products

  // Gammas
  double gamma_e_plus  = 0;
  double gamma_e_minus = 0;

  // Gammas
  for (const HepMC3::GenParticlePtr &p1 : search_gammas) {
    M4Vec pvec = gra::aux::HepMC2M4Vec(p1->momentum());

    // Check that ancestor is excited forward system
    std::vector<HepMC3::GenParticlePtr> ancestor =
        HepMC3::applyFilter(HepMC3::StandardSelector::PDG_ID == PDG::PDG_NSTAR ||
                                HepMC3::StandardSelector::PDG_ID == -PDG::PDG_NSTAR,
                            HepMC3::Relatives::ANCESTORS(p1));

    if (ancestor.size() != 0) {
      hEta_Gamma->Fill(pvec.Eta(), W);
      hE_Gamma->Fill(pvec.E(), W);
      hXF_Gamma->Fill(pvec.Pz() / (sqrts / 2), W);

      if (excited_plus && pvec.Rap() > 0) { gamma_e_plus += pvec.E(); }
      if (excited_minus && pvec.Rap() < 0) { gamma_e_minus += pvec.E(); }
    }
  }

  // Pi+
  for (const HepMC3::GenParticlePtr &p1 : search_pip) {
    // HepMC3::Print::line(p1);

    // Check that parent is the excited system
    std::vector<HepMC3::GenParticlePtr> parents = p1->parents();
    bool                                found   = false;
    for (const auto &k : indices(parents)) {
      if (std::abs(parents[k]->pid()) == PDG::PDG_NSTAR) { found = true; }
    }
    if (found) {
      M4Vec pvec = gra::aux::HepMC2M4Vec(p1->momentum());
      hEta_Pions->Fill(pvec.Eta(), W);
      hE_Pions->Fill(pvec.E(), W);
      hXF_Pions->Fill(pvec.Pz() / (sqrts / 2), W);
    }
  }

  // Pi-
  for (const HepMC3::GenParticlePtr &p1 : search_pim) {
    // HepMC3::Print::line(p1);

    // Check that parent is the excited system
    std::vector<HepMC3::GenParticlePtr> parents = p1->parents();
    bool                                found   = false;
    for (const auto &k : indices(parents)) {
      if (std::abs(parents[k]->pid()) == PDG::PDG_NSTAR) { found = true; }
    }
    if (found) {
      M4Vec pvec = gra::aux::HepMC2M4Vec(p1->momentum());
      hEta_Pions->Fill(pvec.Eta(), W);
      hE_Pions->Fill(pvec.E(), W);
      hXF_Pions->Fill(pvec.Pz() / (sqrts / 2), W);
    }
  }

  // Neutrons
  double neutron_e_plus  = 0;
  double neutron_e_minus = 0;

  for (const HepMC3::GenParticlePtr &p1 : search_neutrons) {
    M4Vec pvec = gra::aux::HepMC2M4Vec(p1->momentum());

    // Check that parent is the excited system
    std::vector<HepMC3::GenParticlePtr> parents = p1->parents();
    bool                                found   = false;
    for (const auto &k : indices(parents)) {
      if (std::abs(parents[k]->pid()) == PDG::PDG_NSTAR) { found = true; }
    }
    if (found) {
      hEta_Neutron->Fill(pvec.Eta(), W);
      hE_Neutron->Fill(pvec.E(), W);
      hXF_Neutron->Fill(pvec.Pz() / (sqrts / 2), W);

      if (excited_plus && pvec.Rap() > 0) { neutron_e_plus += pvec.E(); }
      if (excited_minus && pvec.Rap() < 0) { neutron_e_minus += pvec.E(); }
    }
  }

  // Gamma+Neutron energy histogram
  if (excited_plus) { hE_GammaNeutron->Fill(gamma_e_plus + neutron_e_plus, W); }
  if (excited_minus) { hE_GammaNeutron->Fill(gamma_e_minus + neutron_e_minus, W); }
}

double powerlaw(double *x, double *par) {
  return par[0] / std::pow(1.0 + std::pow(x[0], 2) / (std::pow(par[1], 2) * par[2]), par[2]);
}

double exponential(double *x, double *par) { return par[0] * exp(par[1] * x[0]); }

// Custom plotter
void MAnalyzer::PlotAll(const std::string &titlestr) {
  // Create output directory if it does not exist
  const std::string FOLDER = gra::aux::GetBasePath(2) + "/figs/" + inputfile;
  aux::CreateDirectory(FOLDER);

  /*
  // FIT FUNCTIONS
  std::shared_ptr<TF1> fb = std::make_shared<TF1>("exp_fit", exponential, 0.05,
  0.5, 2);
  fb->SetParameter(0, 10.0); // A
  fb->SetParameter(1, -8.0); // b

  std::shared_ptr<TF1> fa = std::make_shared<TF1>("pow_fit", powerlaw, 0.5, 3.0,
  3);
  fa->SetParameter(0, 10.0); // A
  fa->SetParameter(1, 0.15); // T
  fa->SetParameter(2, 0.1);  // n
  */
  //        hpT2->Fit("exp_fit","R");
  //        hpT_Meson_p->Fit("pow_fit","R"); // "R" for range

  // *************** FORWARD EXCITATION ***************
  if (N_STAR_ON) {
    // Draw histograms
    TCanvas c1("c", "c", 800, 600);
    c1.Divide(2, 2, 0.0002, 0.0002);

    c1.cd(1);
    // gPad->SetLogy();
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
    // gPad->SetLogy();
    hXF_Gamma->SetLineColor(2);
    hXF_Pions->SetLineColor(4);
    hXF_Neutron->SetLineColor(8);

    hXF_Gamma->Draw("same");
    hXF_Pions->Draw("same");
    hXF_Neutron->Draw("same");

    c1.cd(3);
    gPad->SetLogy();
    // gPad->SetLogx();
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

  // -------------------------------------------------------------------------------------
  // FRAME correlations

  TCanvas c2("c", "c", 800, 800);
  c2.Divide(analyzer::FRAMES.size(), analyzer::FRAMES.size(), 0.0001, 0.0002);

  int k = 1;
  for (std::size_t i = 0; i < analyzer::FRAMES.size(); ++i) {
    for (std::size_t j = 0; j < analyzer::FRAMES.size(); ++j) {
      c2.cd(k);
      ++k;
      if (j >= i) { h2CosTheta[i][j]->Draw("COLZ"); }

      // Titlestr
      if ((i == 0) & (j == 0)) { h2CosTheta[i][j]->SetTitle(titlestr.c_str()); }
    }
  }
  c2.SaveAs(Form("%s/figs/%s/h2_frame_correlations_costheta.pdf", gra::aux::GetBasePath(2).c_str(),
                 inputfile.c_str()));

  // -------------------------------------------------------------------------------------

  TCanvas c3("c", "c", 800, 800);
  c3.Divide(analyzer::FRAMES.size(), analyzer::FRAMES.size(), 0.0001, 0.0002);

  k = 1;
  for (std::size_t i = 0; i < analyzer::FRAMES.size(); ++i) {
    for (std::size_t j = 0; j < analyzer::FRAMES.size(); ++j) {
      c3.cd(k);
      ++k;
      if (j >= i) { h2Phi[i][j]->Draw("COLZ"); }

      // Titlestr
      if ((i == 0) & (j == 0)) { h2Phi[i][j]->SetTitle(titlestr.c_str()); }
    }
  }
  c3.SaveAs(Form("%s/figs/%s/h2_frame_correlations_phi.pdf", gra::aux::GetBasePath(2).c_str(),
                 inputfile.c_str()));

  // -------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------
  // Legendre polynomials in the Rest Frame (non-rotated one)

  const int    colors[8] = {48, 53, 98, 32, 48, 53, 98, 32};
  const double XBOUND[2] = {0.0, 4.0};

  // 1...8
  TCanvas *c115 = new TCanvas("c115", "Legendre polynomials 1-4", 600, 400);
  TCanvas *c116 = new TCanvas("c116", "Legendre polynomials 5-8", 600, 400);

  TLegend *leg[8];
  TLine *  line[8];
  c115->Divide(2, 2, 0.001, 0.001);
  c116->Divide(2, 2, 0.001, 0.001);

  for (std::size_t l = 0; l < 8; ++l) {
    if (l < 4) {
      c115->cd(l + 1);
    } else {
      c116->cd(l - 3);
    }

    // Legend
    leg[l] = new TLegend(0.15, 0.75, 0.4, 0.85);  // x1,y1,x2,y2

    // Data
    hPl[l]->SetLineColor(colors[l]);
    hPl[l]->Draw("SAME");
    hPl[l]->SetMinimum(-0.4);  // Y-axis minimum
    hPl[l]->SetMaximum(0.4);   // Y-axis maximum

    // Adjust legend
    leg[l]->SetFillColor(0);   // White background
    leg[l]->SetBorderSize(0);  // No box
    leg[l]->AddEntry(hPl[l].get(), Form("l = %lu", l + 1), "l");
    leg[l]->Draw("SAME");

    // Horizontal line
    line[l] = new TLine(XBOUND[0], 0, XBOUND[1], 0);
    line[l]->Draw("SAME");

    // Titlestr
    if (l == 0 || l == 4) { hPl[l]->SetTitle(titlestr.c_str()); }
  }
  c115->SaveAs(
      Form("%s/figs/%s/hPl_1to4.pdf", gra::aux::GetBasePath(2).c_str(), inputfile.c_str()));
  c116->SaveAs(
      Form("%s/figs/%s/hPl_6to8.pdf", gra::aux::GetBasePath(2).c_str(), inputfile.c_str()));

  for (std::size_t i = 0; i < 8; ++i) {
    delete leg[i];
    delete line[i];
  }
  delete c115;
  delete c116;
}

}  // namespace gra
