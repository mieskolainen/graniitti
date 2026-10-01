// Fast MC analysis class
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <algorithm>
#include <cmath>
#include <exception>
#include <filesystem>
#include <iostream>
#include <limits>
#include <memory>
#include <utility>
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
#include "HepMC3/Relatives.h"
#include "HepMC3/Selector.h"

// Own
#include "Graniitti/Analysis/MAnalyzer.h"
#include "Graniitti/Analysis/MHepMCReader.h"
#include "Graniitti/Analysis/MMultiplet.h"
#include "Graniitti/Kinematics/M4Vec.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MSpecialFunctions.h"
#include "Graniitti/Particle/MPDG.h"
#include "Graniitti/Math/MStatistics.h"

using gra::aux::indices;
using gra::math::msqrt;

namespace {

// Recognize central decays from a system record or the two exchanged momenta
bool IsCentralDecay(const HepMC3::ConstGenParticlePtr &particle) {
  for (const auto &ancestor : HepMC3::Relatives::ANCESTORS(particle)) {
    if (ancestor->pid() == gra::PDG::PDG_system || ancestor->pid() == -gra::PDG::PDG_system) {
      return true;
    }
    const auto parents = ancestor->parents();
    if (parents.size() == 2 &&
        std::all_of(parents.begin(), parents.end(), [](const auto &parent) {
          return parent->pid() == gra::PDG::PDG_propagator;
        })) {
      return true;
    }
  }
  return false;
}

// Identify forward excitation descendants through the full decay cascade
bool IsNStarDecay(const HepMC3::ConstGenParticlePtr &particle) {
  const auto ancestors = HepMC3::Relatives::ANCESTORS(particle);
  return std::any_of(ancestors.begin(), ancestors.end(), [](const auto &ancestor) {
    return ancestor->pid() == gra::PDG::PDG_NSTAR || ancestor->pid() == -gra::PDG::PDG_NSTAR;
  });
}

// Create a histogram with explicit C++ ownership outside ROOT directories
template <typename Histogram, typename... Args>
std::shared_ptr<Histogram> MakeDetachedHistogram(Args &&...args) {
  auto histogram = std::make_shared<Histogram>(std::forward<Args>(args)...);
  histogram->SetDirectory(nullptr);
  return histogram;
}

// Validate every histogram required by the selected final-state multiplicity
template <typename Multiplet>
void RequireHistograms(const std::map<std::string, std::shared_ptr<Multiplet>> &histograms,
                       const std::vector<std::string> &names, unsigned int sample) {
  for (const auto &name : names) {
    const auto found = histograms.find(name);
    if (found == histograms.end() || !found->second ||
        sample >= found->second->h.size() || !found->second->h[sample]) {
      throw std::invalid_argument("MAnalyzer::HepMC3Read: missing histogram or invalid sample for " + name);
    }
  }
}

// Check the histogram set before reading or filling any events
void ValidateHistograms(const std::map<std::string, std::shared_ptr<gra::h1Multiplet>> &h1,
                        const std::map<std::string, std::shared_ptr<gra::h2Multiplet>> &h2,
                        const std::map<std::string, std::shared_ptr<gra::hProfMultiplet>> &hP,
                        unsigned int sample, unsigned int multiplicity) {
  RequireHistograms(h1, {"h1_1B_eta", "h1_1B_pt", "h1_PP_dphi", "h1_PP_dpt", "h1_PP_t1",
                         "h1_S_M", "h1_S_Pt", "h1_S_Pt2", "h1_S_Y"}, sample);
  RequireHistograms(h2, {"h2_S_M_Pt", "h2_S_M_dphipp", "h2_S_M_dpt", "h2_S_M_pt", "h2_S_M_t"}, sample);
  RequireHistograms(hP, {"hP_S_M_Pt"}, sample);
  if (multiplicity != 2) { return; }
  RequireHistograms(h1, {"h1_2B_acop", "h1_2B_diffrap"}, sample);
  RequireHistograms(h2, {"h2_2B_M_dphi", "h2_2B_eta1_eta2"}, sample);
  RequireHistograms(hP, {"hP_2B_M_dphi", "hP_S_M_PL2_CM", "hP_S_M_PL4_CM"}, sample);
  for (const auto &frame : gra::analyzer::Frames()) {
    RequireHistograms(h1, {"h1_costheta_" + frame, "h1_phi_" + frame}, sample);
    RequireHistograms(h2, {"h2_2B_M_costheta_" + frame, "h2_2B_M_phi_" + frame,
                           "h2_2B_costheta_phi_" + frame}, sample);
  }
}

}  // namespace

namespace gra {

// Constructor with unique ID string for ROOT bookkeeping reasons
MAnalyzer::MAnalyzer(const std::string &ID) {
  // Initialize histograms
  const int NBINS = 150;

  // Energy
  hE_Pions = MakeDetachedHistogram<TH1D>(Form("%s_%s", "Energy #pi (GeV)", ID.c_str()),
                                         ";Energy (GeV);Events", NBINS, 0, 1.0);
  hE_Gamma = MakeDetachedHistogram<TH1D>(Form("%s_%s", "Energy #gamma (GeV)", ID.c_str()),
                                         ";Energy (GeV);Events", NBINS, 0, 1.0);
  hE_Neutron = MakeDetachedHistogram<TH1D>(Form("%s_%s", "Energy n (GeV)", ID.c_str()),
                                           ";Energy (GeV);Events", NBINS, 0, 1.0);
  hE_GammaNeutron =
      MakeDetachedHistogram<TH1D>(Form("%s_%s", "Energy y+n (GeV)", ID.c_str()),
                                  ";Energy (GeV);Events", NBINS, 0, 1.0);

  // Feynman-x
  hXF_Pions = MakeDetachedHistogram<TH1D>(Form("%s_%s", "xF #pi", ID.c_str()),
                                          ";Feynman-x;Events", NBINS, -1.0, 1.0);
  hXF_Gamma = MakeDetachedHistogram<TH1D>(Form("%s_%s", "xF #gamma", ID.c_str()),
                                          ";Feynman-x;Events", NBINS, -1.0, 1.0);
  hXF_Neutron = MakeDetachedHistogram<TH1D>(Form("%s_%s", "xF n", ID.c_str()),
                                            ";Feynman-x;Events", NBINS, -1.0, 1.0);

  // Forward systems
  hEta_Pions = MakeDetachedHistogram<TH1D>(Form("%s_%s", "#eta pi", ID.c_str()), ";#eta;Events",
                                           NBINS, -12, 12);
  hEta_Gamma = MakeDetachedHistogram<TH1D>(Form("%s_%s", "#eta y", ID.c_str()), ";#eta;Events",
                                           NBINS, -12, 12);
  hEta_Neutron = MakeDetachedHistogram<TH1D>(Form("%s_%s", "#eta n", ID.c_str()), ";#eta;Events",
                                             NBINS, -12, 12);
  hM_NSTAR = MakeDetachedHistogram<TH1D>(Form("%s_%s", "M (GeV)", ID.c_str()),
                                         ";M (GeV);Events", NBINS, 0, 10);

  // Legendre polynomials, DO NOT CHANGE THE Y-RANGE [-1,1]
  for (std::size_t i = 0; i < 8; ++i) {
    hPl[i] = MakeDetachedHistogram<TProfile>(Form("hPl%lu_%s", i + 1, ID.c_str()), "", 100, 0.0,
                                             4.0, -1, 1);
    hPl[i]->Sumw2();  // Error saving on
    hPl[i]->SetXTitle(Form("System M (GeV)"));
    hPl[i]->SetYTitle(Form("Legendre #LTP_{l}(cos #theta)#GT [CM frame]"));
  }

  // Costheta correlations between different frames
  const auto frames = analyzer::Frames();
  for (const auto &i : indices(frames)) {
    for (const auto &j : indices(frames)) {
      h2CosTheta[i][j] = MakeDetachedHistogram<TH2D>(
          Form("%s^{+} cos(theta) %s vs %s_%s", pstr.c_str(), frames[i].data(),
               frames[j].data(), ID.c_str()),
          Form(";%s^{+} cos(#theta) %s;%s^{+} cos(#theta) %s", pstr.c_str(),
               frames[i].data(), pstr.c_str(), frames[j].data()),
          NBINS, -1, 1, NBINS, -1, 1);
    }
  }

  // Phi correlations between different frames
  for (const auto &i : indices(frames)) {
    for (const auto &j : indices(frames)) {
      h2Phi[i][j] = MakeDetachedHistogram<TH2D>(
          Form("%s^{+} #phi %s vs %s_%s", pstr.c_str(), frames[i].data(),
               frames[j].data(), ID.c_str()),
          Form(";%s^{+} #phi %s (rad);%s^{+} #phi %s (rad)", pstr.c_str(),
               frames[i].data(), pstr.c_str(), frames[j].data()),
          NBINS, -gra::math::PI, gra::math::PI, NBINS, -gra::math::PI, gra::math::PI);
    }
  }
}

// Destructor
MAnalyzer::~MAnalyzer() {}

// Configure energy axes once the collider energy is known
void MAnalyzer::ConfigureColliderEnergy(const M4Vec &collision) {
  const double collider_energy = collision.M();
  if (!std::isfinite(collider_energy) || collider_energy <= 0.0 ||
      !std::isfinite(collision.E()) || collision.E() <= 0.0) {
    throw std::invalid_argument("MAnalyzer::ConfigureColliderEnergy: invalid collider energy");
  }
  if (energy_range_initialized) {
    const double tolerance =
        1e-9 * std::max({1.0, std::abs(sqrts), std::abs(collider_energy)});
    if (std::abs(collider_energy - sqrts) > tolerance) {
      throw std::invalid_argument(
          "MAnalyzer::ConfigureColliderEnergy: collider energy changes within one sample");
    }
    return;
  }
  sqrts = collider_energy;
  constexpr int energy_bins = 150;
  for (const auto &histogram : {hE_Pions, hE_Gamma, hE_Neutron, hE_GammaNeutron}) {
    histogram->SetBins(energy_bins, 0.0, collision.E());
  }
  energy_range_initialized = true;
}

// Fill histograms using generator ancestry and final-state momenta
double MAnalyzer::HepMC3_OracleFill(const std::string input, unsigned int multiplicity,
                                    int finalPDG, unsigned int MAXEVENTS,
                                    std::map<std::string, std::shared_ptr<h1Multiplet>> &   h1,
                                    std::map<std::string, std::shared_ptr<h2Multiplet>> &   h2,
                                    std::map<std::string, std::shared_ptr<hProfMultiplet>> &hP,
                                    unsigned int                                            SID) {
  if (multiplicity == 0 || finalPDG == 0 || finalPDG == std::numeric_limits<int>::min()) {
    throw std::invalid_argument("MAnalyzer::HepMC3Read: invalid multiplicity or PDG code");
  }
  ValidateHistograms(h1, h2, hP, SID, multiplicity);

  inputfile                     = input;
  const std::string totalpath = (std::filesystem::path(gra::aux::GetBasePath(2)) / "output" /
                                 (input + ".hepmc3")).lexically_normal().string();
  MHepMCReader input_file(totalpath);

  // Event loop
  unsigned int events_read = 0;

  // Variables for calculating selection efficiency
  statistics::ScaledWeightSums total_weights;
  statistics::ScaledWeightSums selected_weights;
  std::size_t selected_events = 0;
  bool cross_section_seen = false;

  // ---------------------------------------------------------------------
  // Set final state [charged pair or neutral pair]
  MPDG PDG;
  PDG.ReadParticleData();

  // Try to find the particle from PDG table, will throw exception if fails
  PDG.FindByPDG(finalPDG);
  const bool has_antiparticle = PDG.PDG_table.contains(-finalPDG);
  const int NEGfinalPDG = has_antiparticle ? -finalPDG : 0;
  // ---------------------------------------------------------------------

  while (events_read < MAXEVENTS) {
    // Read event from input file
    HepMC3::GenEvent evt(HepMC3::Units::GEV, HepMC3::Units::MM);
    if (!input_file.Read(evt)) {
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
      const double event_cross_section = 1E-12 * cs->xsec(0);
      if (!std::isfinite(event_cross_section) || event_cross_section < 0.0) {
        throw std::invalid_argument(
            "MAnalyzer::HepMC3Read: invalid GenCrossSection in " + totalpath);
      }
      if (cross_section_seen) {
        const double tolerance = 1e-2 * std::max(event_cross_section, cross_section);
        if (std::abs(event_cross_section - cross_section) > tolerance) {
          throw std::invalid_argument(
              "MAnalyzer::HepMC3Read: inconsistent GenCrossSection in " + totalpath);
        }
      } else {
        cross_section = event_cross_section;
      }
      cross_section_seen = true;
    } else {
      throw std::invalid_argument(
          "MAnalyzer::HepMC3Read: missing GenCrossSection in " + totalpath);
    }
    // --------------------------------------------------------------

    // Get the nominal event weight before cross-section normalization
    double W = 1.0;
    if (evt.weights().size() != 0) {  // check do we have weights saved
      W = evt.weights()[0];           // take the first one
    }
    if (!std::isfinite(W)) {
      throw std::invalid_argument(
          "MAnalyzer::HepMC3Read: non-finite event weight in " + totalpath);
    }
    total_weights.Add(W);
    // --------------------------------------------------------------

    // Central particles
    std::vector<M4Vec> pip;
    std::vector<M4Vec> pim;

    for (HepMC3::ConstGenParticlePtr p1 :
         HepMC3::applyFilter(HepMC3::StandardSelector::STATUS == PDG::PDG_STABLE &&
                                 HepMC3::StandardSelector::PDG_ID == finalPDG,
                             evt.particles())) {
      M4Vec pvec = gra::aux::HepMC2M4Vec(p1->momentum());

      if (IsCentralDecay(p1)) { pip.push_back(pvec); }
    }
    for (HepMC3::ConstGenParticlePtr p1 :
         HepMC3::applyFilter(HepMC3::StandardSelector::STATUS == PDG::PDG_STABLE &&
                                 HepMC3::StandardSelector::PDG_ID == NEGfinalPDG,
                             evt.particles())) {
      M4Vec pvec = gra::aux::HepMC2M4Vec(p1->momentum());

      if (IsCentralDecay(p1)) { pim.push_back(pvec); }
    }

    // CHECK CONDITION
    const bool valid_topology =
        (multiplicity == 1 || !has_antiparticle)
            ? pip.size() == multiplicity && pim.empty()
            : (multiplicity == 2
                   ? pip.size() == 1 && pim.size() == 1
                   : pip.size() + pim.size() == multiplicity && !pip.empty() && !pim.empty());
    if (!valid_topology) {
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
                                *abs(HepMC3::StandardSelector::PDG_ID) == PDG::PDG_p,
                            evt.particles());

    std::vector<HepMC3::GenParticlePtr> final_protons =
        HepMC3::applyFilter(HepMC3::StandardSelector::STATUS == PDG::PDG_STABLE &&
                                *abs(HepMC3::StandardSelector::PDG_ID) == PDG::PDG_p,
                            evt.particles());

    M4Vec p_beam_plus;
    M4Vec p_beam_minus;
    M4Vec p_final_plus;
    M4Vec p_final_minus;
    std::size_t beam_plus_count  = 0;
    std::size_t beam_minus_count = 0;
    std::size_t final_plus_count  = 0;
    std::size_t final_minus_count = 0;

    // Beam (initial state ) protons
    for (const HepMC3::GenParticlePtr &p1 : beam_protons) {
      M4Vec pvec = gra::aux::HepMC2M4Vec(p1->momentum());
      if (pvec.Rap() > 0) {
        p_beam_plus = pvec;
        ++beam_plus_count;
      } else {
        p_beam_minus = pvec;
        ++beam_minus_count;
      }
    }

    // Final state protons
    for (const HepMC3::GenParticlePtr &p1 : final_protons) {
      M4Vec pvec = gra::aux::HepMC2M4Vec(p1->momentum());

      // Select intact forward protons outside central and N* decays
      if (!IsNStarDecay(p1) && !IsCentralDecay(p1)) {
        if (pvec.Rap() > 0) {
          p_final_plus = pvec;
          ++final_plus_count;
        } else {
          p_final_minus = pvec;
          ++final_minus_count;
        }
      }
    }
    const bool has_beams = beam_plus_count == 1 && beam_minus_count == 1;
    const bool has_forward_protons =
        final_plus_count == 1 && final_minus_count == 1;
    if (has_beams) {
      CheckEnergyMomentum(evt);
      ConfigureColliderEnergy(p_beam_plus + p_beam_minus);
    }
    if (!has_beams) {
      p_beam_plus = M4Vec();
      p_beam_minus = M4Vec();
    }
    if (!has_forward_protons) {
      p_final_plus = M4Vec();
      p_final_minus = M4Vec();
    }

    // Observables for 2-body case only
    if (multiplicity == 2) {
      FrameObservables(W, p_beam_plus, p_beam_minus, p_final_plus, p_final_minus, pip, pim);
    }

    // Observables for N stars
    if (sqrts > 0.0) { NStarObservables(W, evt); }

    // **************************************************************
    // SUPERPLOTTER >>
    try {
      const M4Vec &a = pip.front();

      const double M  = system.M();
      const double Pt = system.Pt();
      const double Y  = system.Rap();

      // 1D: System
      h1.at("h1_S_M")->h.at(SID)->Fill(M, W);
      h1.at("h1_S_Pt")->h.at(SID)->Fill(Pt, W);
      h1.at("h1_S_Pt2")->h.at(SID)->Fill(math::pow2(Pt), W);
      h1.at("h1_S_Y")->h.at(SID)->Fill(Y, W);
      hP.at("hP_S_M_Pt")->h.at(SID)->Fill(M, Pt, W);

      // 1D: 1-Body
      h1.at("h1_1B_pt")->h.at(SID)->Fill(a.Pt(), W);
      h1.at("h1_1B_eta")->h.at(SID)->Fill(a.Eta(), W);

      // 1D: Forward proton pair
      double deltaphi_pp = -1.0;
      if (has_beams && has_forward_protons) {
        // Mandelstam -t_1,2
        const double t1 = -(p_beam_plus - p_final_plus).M2();
        // const double t2 = -(p_beam_minus - p_final_minus).M2();

        // Deltaphi
        deltaphi_pp          = p_final_plus.DeltaPhiAbs(p_final_minus);
        M4Vec        pp_diff = p_final_plus - p_final_minus;
        const double pp_dpt  = pp_diff.Pt();

        h1.at("h1_PP_dphi")->h.at(SID)->Fill(deltaphi_pp, W);
        h1.at("h1_PP_t1")->h.at(SID)->Fill(t1, W);
        h1.at("h1_PP_dpt")->h.at(SID)->Fill(pp_dpt, W);

        h2.at("h2_S_M_dphipp")->h.at(SID)->Fill(M, deltaphi_pp, W);
        h2.at("h2_S_M_dpt")->h.at(SID)->Fill(M, pp_dpt, W);
        h2.at("h2_S_M_t")->h.at(SID)->Fill(M, std::abs(t1), W);
      }

      // 2D
      h2.at("h2_S_M_Pt")->h.at(SID)->Fill(M, Pt, W);
      h2.at("h2_S_M_pt")->h.at(SID)->Fill(M, a.Pt(), W);

      // 2-Body only
      if (multiplicity == 2) {
        const M4Vec &b = pim.empty() ? pip[1] : pim.front();
        const double delta_phi = a.DeltaPhiAbs(b);
        hP.at("hP_2B_M_dphi")->h.at(SID)->Fill(M, delta_phi, W);
        h1.at("h1_2B_acop")->h.at(SID)->Fill(1.0 - delta_phi / gra::math::PI, W);
        h1.at("h1_2B_diffrap")->h.at(SID)->Fill(a.Rap() - b.Rap(), W);
        h2.at("h2_2B_M_dphi")->h.at(SID)->Fill(M, delta_phi, W);
        h2.at("h2_2B_eta1_eta2")->h.at(SID)->Fill(a.Eta(), b.Eta(), W);


        // Frame transform
        const int   direction = 1;
        const M4Vec X         = a + b;

        std::vector<M4Vec> CM = {a, b};
        gra::kinematics::CMframe(CM, X);

        std::vector<M4Vec> HX = {a, b};
        gra::kinematics::HXframe(HX, X);

        h1.at("h1_costheta_CM")->h.at(SID)->Fill(CM[0].CosTheta(), W);
        h1.at("h1_costheta_HX")->h.at(SID)->Fill(HX[0].CosTheta(), W);
        h1.at("h1_costheta_LAB")->h.at(SID)->Fill(a.CosTheta(), W);


        h1.at("h1_phi_CM")->h.at(SID)->Fill(CM[0].Phi(), W);
        h1.at("h1_phi_HX")->h.at(SID)->Fill(HX[0].Phi(), W);
        h1.at("h1_phi_LAB")->h.at(SID)->Fill(a.Phi(), W);


        h2.at("h2_2B_costheta_phi_CM")->h.at(SID)->Fill(CM[0].CosTheta(), CM[0].Phi(), W);
        h2.at("h2_2B_costheta_phi_HX")->h.at(SID)->Fill(HX[0].CosTheta(), HX[0].Phi(), W);
        h2.at("h2_2B_costheta_phi_LAB")->h.at(SID)->Fill(a.CosTheta(), a.Phi(), W);


        h2.at("h2_2B_M_costheta_CM")->h.at(SID)->Fill(M, CM[0].CosTheta(), W);
        h2.at("h2_2B_M_costheta_HX")->h.at(SID)->Fill(M, HX[0].CosTheta(), W);
        h2.at("h2_2B_M_costheta_LAB")->h.at(SID)->Fill(M, a.CosTheta(), W);


        h2.at("h2_2B_M_phi_CM")->h.at(SID)->Fill(M, CM[0].Phi(), W);
        h2.at("h2_2B_M_phi_HX")->h.at(SID)->Fill(M, HX[0].Phi(), W);
        h2.at("h2_2B_M_phi_LAB")->h.at(SID)->Fill(M, a.Phi(), W);

        if (has_beams) {
          std::vector<M4Vec> CS = {a, b};
          gra::kinematics::CSframe(CS, X, p_beam_plus, p_beam_minus);
          std::vector<M4Vec> PG = {a, b};
          gra::kinematics::PGframe(PG, X, direction, p_beam_plus, p_beam_minus);
          h1.at("h1_costheta_CS")->h.at(SID)->Fill(CS[0].CosTheta(), W);
          h1.at("h1_costheta_PG")->h.at(SID)->Fill(PG[0].CosTheta(), W);
          h1.at("h1_phi_CS")->h.at(SID)->Fill(CS[0].Phi(), W);
          h1.at("h1_phi_PG")->h.at(SID)->Fill(PG[0].Phi(), W);
          h2.at("h2_2B_costheta_phi_CS")->h.at(SID)->Fill(CS[0].CosTheta(), CS[0].Phi(), W);
          h2.at("h2_2B_costheta_phi_PG")->h.at(SID)->Fill(PG[0].CosTheta(), PG[0].Phi(), W);
          h2.at("h2_2B_M_costheta_CS")->h.at(SID)->Fill(M, CS[0].CosTheta(), W);
          h2.at("h2_2B_M_costheta_PG")->h.at(SID)->Fill(M, PG[0].CosTheta(), W);
          h2.at("h2_2B_M_phi_CS")->h.at(SID)->Fill(M, CS[0].Phi(), W);
          h2.at("h2_2B_M_phi_PG")->h.at(SID)->Fill(M, PG[0].Phi(), W);
        }
        if (has_beams && has_forward_protons) {
          std::vector<M4Vec> GJ = {a, b};
          gra::kinematics::GJframe(GJ, X, direction, p_beam_plus - p_final_plus,
                                   p_beam_minus - p_final_minus);
          h1.at("h1_costheta_GJ")->h.at(SID)->Fill(GJ[0].CosTheta(), W);
          h1.at("h1_phi_GJ")->h.at(SID)->Fill(GJ[0].Phi(), W);
          h2.at("h2_2B_costheta_phi_GJ")->h.at(SID)->Fill(GJ[0].CosTheta(), GJ[0].Phi(), W);
          h2.at("h2_2B_M_costheta_GJ")->h.at(SID)->Fill(M, GJ[0].CosTheta(), W);
          h2.at("h2_2B_M_phi_GJ")->h.at(SID)->Fill(M, GJ[0].Phi(), W);
        }


        // ---------------------------------------------------------------------------
        hP.at("hP_S_M_PL2_CM")->h.at(SID)->Fill(M, math::LegendrePl(2, CM[0].CosTheta()), W);
        hP.at("hP_S_M_PL4_CM")->h.at(SID)->Fill(M, math::LegendrePl(4, CM[0].CosTheta()), W);
        // ---------------------------------------------------------------------------
      }

      // 4-Body only
      if (multiplicity == 4) {
        // ...
      }
    } catch (const std::exception &e) {
      throw std::invalid_argument("MAnalyzer::HepMC3Read: Problem filling histogram: " +
                                  std::string(e.what()));
    } catch (...) {
      std::throw_with_nested(
          std::runtime_error("MAnalyzer::HepMC3Read: Unknown problem filling histogram"));
    }

    // << SUPERPLOTTER
    // **************************************************************

    if (events_read % 10000 == 0) {
      std::cout << std::endl << "Events processed: " << events_read << std::endl;
    }

    // Sum weights only after the event has passed the analysis selection
    selected_weights.Add(W);
    ++selected_events;
  }
  if (events_read >= MAXEVENTS) {
    std::cout << "MAnalyzer::HepMC3Read: Maximum event count " << MAXEVENTS << " reached!";
  }
  std::cout << std::endl;
  std::cout << "MAnalyzer::HepMC3Read: Events processed in total: " << events_read << std::endl;

  if (selected_events == 0) {
    throw std::invalid_argument("MAnalyzer::HepMC3Read:: Valid events in <" + totalpath + ">" +
                                " == 0 out of " + std::to_string(events_read));
  }
  if (!cross_section_seen ||
      !total_weights.HasSignificantSignedSum(
          std::numeric_limits<double>::epsilon())) {
    throw std::invalid_argument(
        "MAnalyzer::HepMC3Read: invalid cross-section normalization weight sum in <" +
        totalpath + ">");
  }
  // Take into account extra fiducial cut efficiency here
  const long double total_weight = total_weights.Sum();
  const double efficiency =
      static_cast<double>(selected_weights.Sum() / total_weight);
  printf("MAnalyzer::HepMC3Read: Fiducial cut efficiency: %0.3f \n", efficiency);
  std::cout << std::endl;

  return static_cast<double>(
      static_cast<long double>(cross_section) / total_weight);
}

// Check event four-momentum conservation and compute the collision mass
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
void MAnalyzer::FrameObservables(double W, const M4Vec &p_beam_plus,
                                 const M4Vec &p_beam_minus, const M4Vec &p_final_plus,
                                 const M4Vec &p_final_minus, const std::vector<M4Vec> &pip,
                                 const std::vector<M4Vec> &pim) {
  const auto frames = analyzer::Frames();
  // Find index
  const auto ind = [&](const std::string &str) {
    for (const auto &i : indices(frames)) {
      if (frames[i] == str) { return i; }
    }
    throw std::invalid_argument(
        "MAnalyzer::FrameObservables: unknown Lorentz frame " + std::string(str));
  };

  std::vector<M4Vec> pf;

  if (pip.size() != 0 && pim.size() != 0) {  // Charged pair
    pf = {pip[0], pim[0]};
  }
  if (pip.size() == 2 && pim.size() == 0) {  // Neutral pair
    pf = {pip[0], pip[1]};
  }
  if (pf.size() != 2) {
    throw std::invalid_argument("MAnalyzer::FrameObservables: invalid two-body topology");
  }

  // ---------------------------------------------------------------------
  // Lorentz frame transformations

  // Make copies
  std::vector<std::vector<M4Vec>> pions(frames.size(), pf);
  std::vector<bool> valid(frames.size(), false);

  // System
  const M4Vec X         = pf[0] + pf[1];
  const int   direction = 1;  // PG and GJ

  gra::kinematics::CMframe(pions[ind("CM")], X);
  valid[ind("CM")] = true;
  gra::kinematics::HXframe(pions[ind("HX")], X);
  valid[ind("HX")] = true;
  valid[ind("LAB")] = true;
  const bool has_beams = p_beam_plus.E() > 0.0 && p_beam_minus.E() > 0.0;
  const bool has_forward =
      p_final_plus.E() > 0.0 && p_final_minus.E() > 0.0;
  if (has_beams) {
    gra::kinematics::CSframe(pions[ind("CS")], X, p_beam_plus, p_beam_minus);
    gra::kinematics::PGframe(pions[ind("PG")], X, direction, p_beam_plus, p_beam_minus);
    valid[ind("CS")] = true;
    valid[ind("PG")] = true;
  }
  if (has_beams && has_forward) {
    gra::kinematics::GJframe(pions[ind("GJ")], X, direction,
                             p_beam_plus - p_final_plus,
                             p_beam_minus - p_final_minus);
    valid[ind("GJ")] = true;
  }

  // ---------------------------------------------------------------------

  // FILL HISTOGRAMS -->

  // Legendre polynomials P_l cos(theta)
  for (std::size_t l = 0; l < 8; ++l) {  // note l+1
    // Take first daughter [0]
    double value = gra::math::LegendrePl((l + 1), pions[ind("CM")][0].CosTheta());
    hPl[l]->Fill(X.M(), value, W);
  }

  // FRAME correlations
  for (const auto &i : indices(frames)) {
    for (const auto &j : indices(frames)) {
      if (!valid[i] || !valid[j]) { continue; }
      h2CosTheta[i][j]->Fill(pions[i][0].CosTheta(), pions[j][0].CosTheta(), W);
      h2Phi[i][j]->Fill(pions[i][0].Phi(), pions[j][0].Phi(), W);
    }
  }
}

// Forward system observables
void MAnalyzer::NStarObservables(double W, HepMC3::GenEvent &evt) {
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
  if (!excited_plus && !excited_minus) { return; }

  // Define Feynman x in the collision CM even for asymmetric beam energies
  M4Vec collision;
  for (const auto &particle : evt.particles()) {
    if (particle->status() == PDG::PDG_BEAM) {
      collision += gra::aux::HepMC2M4Vec(particle->momentum());
    }
  }
  ConfigureColliderEnergy(collision);

  // Count stable particles from every stage of the forward decay cascade
  double gamma_e_plus = 0.0;
  double gamma_e_minus = 0.0;
  double neutron_e_plus = 0.0;
  double neutron_e_minus = 0.0;
  for (const auto &particle : evt.particles()) {
    if (particle->status() != PDG::PDG_STABLE) { continue; }
    const int pdg = particle->pid();
    const bool neutron = pdg == PDG::PDG_n || pdg == -PDG::PDG_n;
    if (pdg != PDG::PDG_gamma && !neutron &&
        pdg != PDG::PDG_pip && pdg != PDG::PDG_pim) {
      continue;
    }
    if (!IsNStarDecay(particle)) { continue; }
    const M4Vec p = gra::aux::HepMC2M4Vec(particle->momentum());
    M4Vec cm = p;
    gra::kinematics::LorentzBoost(collision, sqrts, cm, -1);
    const double xf = 2.0 * cm.Pz() / sqrts;
    if (pdg == PDG::PDG_gamma) {
      hEta_Gamma->Fill(p.Eta(), W);
      hE_Gamma->Fill(p.E(), W);
      hXF_Gamma->Fill(xf, W);
      if (excited_plus && p.Pz() > 0.0) { gamma_e_plus += p.E(); }
      if (excited_minus && p.Pz() < 0.0) { gamma_e_minus += p.E(); }
    } else if (neutron) {
      hEta_Neutron->Fill(p.Eta(), W);
      hE_Neutron->Fill(p.E(), W);
      hXF_Neutron->Fill(xf, W);
      if (excited_plus && p.Pz() > 0.0) { neutron_e_plus += p.E(); }
      if (excited_minus && p.Pz() < 0.0) { neutron_e_minus += p.E(); }
    } else {
      hEta_Pions->Fill(p.Eta(), W);
      hE_Pions->Fill(p.E(), W);
      hXF_Pions->Fill(xf, W);
    }
  }

  // Gamma+Neutron energy histogram
  if (excited_plus) { hE_GammaNeutron->Fill(gamma_e_plus + neutron_e_plus, W); }
  if (excited_minus) { hE_GammaNeutron->Fill(gamma_e_minus + neutron_e_minus, W); }
}

// Evaluate the transverse-momentum power-law fit function
double powerlaw(double *x, double *par) {
  return par[0] / std::pow(1.0 + std::pow(x[0], 2) / (std::pow(par[1], 2) * par[2]), par[2]);
}

// Evaluate an exponential momentum-transfer fit function
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

  const auto frames = analyzer::Frames();
  TCanvas c2("c", "c", 800, 800);
  c2.Divide(frames.size(), frames.size(), 0.0001, 0.0002);

  int k = 1;
  for (const auto &i : indices(frames)) {
    for (const auto &j : indices(frames)) {
      c2.cd(k);
      ++k;
      if (j >= i) { h2CosTheta[i][j]->Draw("COLZ"); }

      // Titlestr
      if (i == 0 && j == 0) { h2CosTheta[i][j]->SetTitle(titlestr.c_str()); }
    }
  }
  c2.SaveAs(Form("%s/figs/%s/h2_frame_correlations_costheta.pdf", gra::aux::GetBasePath(2).c_str(),
                 inputfile.c_str()));

  // -------------------------------------------------------------------------------------

  TCanvas c3("c", "c", 800, 800);
  c3.Divide(frames.size(), frames.size(), 0.0001, 0.0002);

  k = 1;
  for (const auto &i : indices(frames)) {
    for (const auto &j : indices(frames)) {
      c3.cd(k);
      ++k;
      if (j >= i) { h2Phi[i][j]->Draw("COLZ"); }

      // Titlestr
      if (i == 0 && j == 0) { h2Phi[i][j]->SetTitle(titlestr.c_str()); }
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
