// Container class for different type of histograms
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <complex>
#include <iostream>
#include <vector>

// Own
#include "Graniitti/MAux.h"
#include "Graniitti/MForm.h"
#include "Graniitti/MKinematics.h"
#include "Graniitti/MMath.h"
#include "Graniitti/MSpin.h"
#include "Graniitti/MUserHistograms.h"

using gra::aux::indices;
using gra::math::PI;

namespace gra {
void MUserHistograms::InitHistograms() {
  // Level 1
  unsigned int Nbins = 40;

  h1["M"]   = MH1<double>(Nbins, "Central System M (GeV)");
  h1["Rap"] = MH1<double>(Nbins, "Central System Rap");
  h1["Rap"].SetAutoSymmetry(true);
  h1["Pt"]       = MH1<double>(Nbins, 0.0, 2.5, "Central System Pt (GeV)");
  h1["dPhi_pp"]  = MH1<double>(Nbins, 0.0, gra::math::PI, "Forward deltaphi (rad)");
  h1["pPt"]      = MH1<double>(Nbins, 0.0, 2.0, "Forward Pt (GeV)");
  h2["rap1rap2"] = MH2(Nbins, Nbins, "Rapidity1 vs Rapidity2");
  h2["rap1rap2"].SetAutoSymmetry({true, true});

  h1["FM"] = MH1<double>(Nbins, "Forward M (GeV)");
  h1["m0"] = MH1<double>(Nbins, "Intermediate daughter M (GeV)");

  // Level 2
  Nbins = 40;

  h2["costhetaphi_CS"] =
      MH2(Nbins, -1.0, 1.0, Nbins, -PI, PI, "(cos theta, phi) [CS] [Collins-Soper frame]");
  h2["costhetaphi_HX"] =
      MH2(Nbins, -1.0, 1.0, Nbins, -PI, PI, "(cos theta, phi) [HX] [Helicity frame]");
  h2["costhetaphi_AH"] =
      MH2(Nbins, -1.0, 1.0, Nbins, -PI, PI, "(cos theta, phi) [AH] [Anti-Helicity frame]");
  h2["costhetaphi_PG"] =
      MH2(Nbins, -1.0, 1.0, Nbins, -PI, PI, "(cos theta, phi) [PG] [Pseudo-GJ frame]");
  h2["costhetaphi_GJ"] =
      MH2(Nbins, -1.0, 1.0, Nbins, -PI, PI, "(cos theta, phi) [GJ] [Gottfried-Jackson frame]");
  h2["costhetaphi_CM"] =
      MH2(Nbins, -1.0, 1.0, Nbins, -PI, PI, "(cos theta, phi) [CM] [Direct system rest frame]");

  h2["costhetaphi_LA"] =
      MH2(Nbins, -1.0, 1.0, Nbins, -PI, PI, "(cos theta, phi) [LAB] [Laboratory frame]");

  /*
    h2["costhetaphi_G1"] =
        MH2(Nbins, -1.0, 1.0, Nbins, -PI, PI, "(cos(theta), phi) [G1 frame]");
    h2["costhetaphi_G2"] =
        MH2(Nbins, -1.0, 1.0, Nbins, -PI, PI, "(cos(theta), phi) [G2 frame]");
    h2["costhetaphi_G3"] =
        MH2(Nbins, -1.0, 1.0, Nbins, -PI, PI, "(cos(theta), phi) [G3 frame]");
    h2["costhetaphi_G4"] =
        MH2(Nbins, -1.0, 1.0, Nbins, -PI, PI, "(cos(theta), phi) [G4 frame]");
  */
  // Level 3
  /*
  h2["Mphi_CS"]        = MH2(Nbins,  1.0, 1.5, Nbins, -PI, PI, "(M, phi) [Collins-Soper frame]");
*/
}

// Input as the total event weight
void MUserHistograms::FillHistograms(double totalweight, const gra::LORENTZSCALAR &lts) {
  // Level 1
  if (HIST >= 1) {
    h1["M"].Fill(gra::math::msqrt(lts.m2), totalweight);
    h1["Rap"].Fill(lts.Y, totalweight);
    h1["Pt"].Fill(lts.Pt, totalweight);
    h1["dPhi_pp"].Fill(std::abs(lts.pfinal[1].DeltaPhi(lts.pfinal[2])), totalweight);
    h1["pPt"].Fill(lts.pfinal[1].Pt(), totalweight);

    // Cascade decay
    if (lts.decaytree[0].legs.size() > 0) { h1["m0"].Fill(lts.decaytree[0].p4.M(), totalweight); }

    // Dissociated proton
    if (lts.pfinal[1].M() > 1.0) { h1["FM"].Fill(lts.pfinal[1].M(), totalweight); }
  }

  // Level 2
  if (HIST >= 2) {
    h2["rap1rap2"].Fill(lts.decaytree[0].p4.Rap(), lts.decaytree[1].p4.Rap(), totalweight);
    FillCosThetaPhi(totalweight, lts);
  }
}

// (costheta,phi) of daughter in different Lorentz frames, input as the total
// event weight
void MUserHistograms::FillCosThetaPhi(double totalweight, const gra::LORENTZSCALAR &lts) {
  // Central final states
  std::vector<M4Vec> pf;
  pf.push_back(lts.decaytree[0].p4);
  pf.push_back(lts.decaytree[1].p4);

  // System 4-momentum
  M4Vec X;
  for (const auto &i : indices(lts.decaytree)) { X += lts.decaytree[i].p4; }

  // ------------------------------------------------------------------
  // Choose Pseudo-Gottfried-Jackson beam direction (-1,1)
  const int direction = 1;

  {
    // ** PREPARE LORENTZ TRANSFORMATION COMMON VARIABLES **
    M4Vec              pb1boost;  // beam1 particle boosted
    M4Vec              pb2boost;  // beam2 particle boosted
    std::vector<M4Vec> pfboost;   // central particles boosted
    gra::kinematics::LorentFramePrepare(pf, X, lts.pbeam1, lts.pbeam2, pb1boost, pb2boost, pfboost);

    // ** TRANSFORM TO DIFFERENT LORENTZ FRAMES **
    std::vector<std::string> frametype = {"CS", "HX", "AH", "PG", "CM"};
    for (const auto &k : indices(frametype)) {
      // Transform and histogram
      std::vector<M4Vec> pfout;
      gra::kinematics::LorentzFrame(pfout, pb1boost, pb2boost, pfboost, frametype[k], direction);
      h2["costhetaphi_" + frametype[k]].Fill(std::cos(pfout[0].Theta()), pfout[0].Phi(),
                                             totalweight);
      // h2["Mphi_"+frametype[k]].Fill(gra::math::msqrt(lts.m2), pfout[0].Phi(),
      // totalweight);
    }
  }

  /*
  {
    M4Vec              q1boost;  // beam1 particle boosted
    M4Vec              q2boost;  // beam2 particle boosted
    std::vector<M4Vec> pfboost;   // central particles boosted
    gra::kinematics::LorentFramePrepare(pf, lts.q1, lts.q2, q1boost, q2boost, pfboost);

    // ** TRANSFORM TO DIFFERENT LORENTZ FRAMES **
    std::vector<std::string> frametype = {"G1", "G2", "G3", "G4"};
    for (const auto &k : indices(frametype)) {
      // Transform and histogram
      std::vector<M4Vec> pfout;
      gra::kinematics::LorentzFrame(pfout, q1boost, q2boost, pfboost, frametype[k], direction);
      h2["costhetaphi_" + frametype[k]].Fill(std::cos(pfout[0].Theta()), pfout[0].Phi(),
  totalweight);
      // h2["Mphi_"+frametype[k]].Fill(gra::math::msqrt(lts.m2), pfout[0].Phi(),
      // totalweight);
    }
  }
  */

  // ------------------------------------------------------------------
  // Gottfried-Jackson frame
  std::vector<M4Vec> pfGJ = pf;
  gra::kinematics::GJframe(pfGJ, lts.pfinal[0], direction, lts.q1, lts.q2, false);
  h2["costhetaphi_GJ"].Fill(std::cos(pfGJ[0].Theta()), pfGJ[0].Phi(), totalweight);

  // Laboratory frame
  h2["costhetaphi_LA"].Fill(std::cos(pf[0].Theta()), pf[0].Phi(), totalweight);
}

// Print all histograms out
void MUserHistograms::PrintHistograms() {
  if (HIST >= 1) {
    for (auto const &x : h1) { h1[x.first].Print(); }
  }
  if (HIST >= 2) {
    for (auto const &x : h2) { h2[x.first].Print(); }
  }
}

}  // namespace gra
