// Container class for different type of histograms
//
// (c) 2017-2021 Mikael Mieskolainen
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
#include "json.hpp"


using gra::aux::indices;
using gra::math::PI;

namespace gra {


// Lorentz frames
const std::vector<std::string> frames = {"CS", "HX", "AH", "PG", "GJ", "CM", "LA"};

const std::vector<std::string> description = {
    "Collins-Soper rest",     "Helicity rest", "Anti-Helicity rest", "Pseudo-GJ rest",
    "Gottfried-Jackson rest", "Direct rest",   "Laboratory"};


void MUserHistograms::InitHistograms() {
  // *** Level 1 ***
  unsigned int Nbins = 40;

  h1["M"]   = MH1<double>(Nbins, "Central System M (GeV)");
  h1["Rap"] = MH1<double>(Nbins, "Central System Rap");
  h1["Rap"].SetAutoSymmetry(true);
  h1["Pt"]      = MH1<double>(Nbins, 0.0, 2.5, "Central System Pt (GeV)");
  h1["dPhi_pp"] = MH1<double>(Nbins, 0.0, 180, "Forward deltaphi (deg)");
  h1["pPt"]     = MH1<double>(Nbins, 0.0, 2.0, "Forward Pt (GeV)");
  h1["FM"]      = MH1<double>(Nbins, "Forward M (GeV)");
  h1["m0"]      = MH1<double>(Nbins, "Intermediate daughter M (GeV)");

  // *** Level 2 ***
  h1["|t1+t2|"]  = MH1<double>(Nbins, "|t1 + t2| (GeV^2)");
  h2["rap1rap2"] = MH2(Nbins, Nbins, "Rapidity1 vs Rapidity2");
  h2["rap1rap2"].SetAutoSymmetry({true, true});

  Nbins = 40;
  for (std::size_t i = 0; i < frames.size(); ++i) {
    const std::string desc = frames[i] + "] [" + description[i] + " frame]";

    h1["costheta_" + frames[i]] = MH1<double>(Nbins, -1, 1, "(cos theta)[+] [" + desc);
    h1["phi_" + frames[i]]      = MH1<double>(Nbins, -180, 180, "(phi)[+] (deg) [" + desc);
    h2["costhetaphi_" + frames[i]] =
        MH2(Nbins, -1.0, 1.0, Nbins, -180, 180, "(cos theta, phi)[+] [" + desc);
  }
}


// Input as the total event weight
void MUserHistograms::FillHistograms(double totalweight, const gra::LORENTZSCALAR &lts) {
  // Level 1
  if (HIST >= 1) {
    h1["M"].Fill(gra::math::msqrt(lts.m2), totalweight);
    h1["Rap"].Fill(lts.Y, totalweight);
    h1["Pt"].Fill(lts.Pt, totalweight);
    h1["dPhi_pp"].Fill(gra::math::Rad2Deg(std::abs(lts.pfinal[1].DeltaPhi(lts.pfinal[2]))),
                       totalweight);
    h1["pPt"].Fill(lts.pfinal[1].Pt(), totalweight);

    // Cascade decay
    if (lts.decaytree[0].legs.size() > 0) { h1["m0"].Fill(lts.decaytree[0].p4.M(), totalweight); }

    // Dissociated proton
    if (lts.pfinal[1].M() > 1.0) { h1["FM"].Fill(lts.pfinal[1].M(), totalweight); }
  }

  // Level 2
  if (HIST >= 2) {
    h1["|t1+t2|"].Fill(std::abs(lts.t1 + lts.t2), totalweight);
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

    for (const auto &k : indices(frames)) {
      if (frames[k] == "GJ" || frames[k] == "LA") {
        continue;  // Treated outside this loop
      }

      // Transform and histogram
      std::vector<M4Vec> pfout;
      gra::kinematics::LorentzFrame(pfout, pb1boost, pb2boost, pfboost, frames[k], direction);

      const double costheta = pfout[0].CosTheta();
      const double phi      = gra::math::Rad2Deg(pfout[1].Phi());

      // 1D
      h1["costheta_" + frames[k]].Fill(costheta, totalweight);
      h1["phi_" + frames[k]].Fill(phi, totalweight);

      // 2D
      h2["costhetaphi_" + frames[k]].Fill(costheta, phi, totalweight);
    }
  }

  // ------------------------------------------------------------------
  // Laboratory frame

  {
    const double costheta = pf[0].CosTheta();
    const double phi      = gra::math::Rad2Deg(pf[0].Phi());

    h1["costheta_LA"].Fill(costheta, totalweight);
    h1["phi_LA"].Fill(phi, totalweight);
    h2["costhetaphi_LA"].Fill(costheta, phi, totalweight);
  }

  // ------------------------------------------------------------------
  // Gottfried-Jackson frame

  {
    std::vector<M4Vec> pfGJ = pf;
    gra::kinematics::GJframe(pfGJ, lts.pfinal[0], direction, lts.q1, lts.q2, false);

    const double costheta = pfGJ[0].CosTheta();
    const double phi      = gra::math::Rad2Deg(pfGJ[0].Phi());

    h1["costheta_GJ"].Fill(costheta, totalweight);
    h1["phi_GJ"].Fill(phi, totalweight);
    h2["costhetaphi_GJ"].Fill(costheta, phi, totalweight);
  }
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

// Save all histograms out
void MUserHistograms::SaveHistograms(const std::string filename) {
  std::ofstream file(filename);

  if (HIST >= 1) {
    nlohmann::json j;

    for (auto const &x : h1) {
      nlohmann::json jthis;
      h1[x.first].struct2json(jthis);
      j["h1"][x.first] = jthis;
    }

    file << j << std::endl;
    std::cout << "MUserHistograms::SaveHistograms: JSON to " << filename << std::endl;
  }
  /*
  // To be implemented
  if (HIST >= 2) {
    for (auto const &x : h2) {
      h2[x.first].RawOutput();
    }
  }
  */
  std::cout << std::endl;
}

}  // namespace gra
