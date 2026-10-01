// Container class for different type of histograms
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <array>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string_view>
#include <vector>

// Own
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Particle/MForm.h"
#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Spin/MHelicity.h"
#include "Graniitti/MUserHistograms.h"
#include "json.hpp"


using gra::aux::indices;
using gra::math::PI;

namespace gra {

// Set the supported histogram detail level
void MUserHistograms::SetHistograms(unsigned int in) {
  if (in > 2) { throw std::invalid_argument("MUserHistograms::SetHistograms: HIST must be 0, 1 or 2"); }
  HIST = in;
}

namespace {

// Compute the Lorentz-frame names used by the fast angular histograms
std::array<std::string_view, 7> HistogramFrames() {
  return {"CS", "HX", "AH", "PG", "GJ", "CM", "LA"};
}

// Compute the descriptions corresponding to the fast angular frames
std::array<std::string_view, 7> HistogramFrameDescriptions() {
  return {"Collins-Soper rest",     "Helicity rest",
          "Anti-Helicity rest",     "Pseudo-GJ rest",
          "Gottfried-Jackson rest", "Direct rest",
          "Laboratory"};
}

}  // namespace


// Initialize the fast event histograms
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
  const auto frames = HistogramFrames();
  const auto descriptions = HistogramFrameDescriptions();
  for (const auto &i : indices(frames)) {
    const std::string frame(frames[i]);
    const std::string desc =
        frame + "] [" + std::string(descriptions[i]) + " frame]";

    h1["costheta_" + frame] =
        MH1<double>(Nbins, -1, 1, "(cos theta)[+] [" + desc);
    h1["phi_" + frame] =
        MH1<double>(Nbins, -180, 180, "(phi)[+] (deg) [" + desc);
    h2["costhetaphi_" + frame] =
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
    if (lts.pfinal.size() >= 3) {
      h1["dPhi_pp"].Fill(
          gra::math::Rad2Deg(
              std::abs(lts.pfinal[1].DeltaPhi(lts.pfinal[2]))),
          totalweight);
      h1["pPt"].Fill(lts.pfinal[1].Pt(), totalweight);

      // Dissociated proton
      if (lts.pfinal[1].M() > 1.0) {
        h1["FM"].Fill(lts.pfinal[1].M(), totalweight);
      }
    }

    // Cascade decay
    if (!lts.decaytree.empty() && !lts.decaytree[0].legs.empty()) {
      h1["m0"].Fill(lts.decaytree[0].p4.M(), totalweight);
    }
  }

  // Level 2
  if (HIST >= 2) {
    h1["|t1+t2|"].Fill(std::abs(lts.t1 + lts.t2), totalweight);
    if (lts.decaytree.size() >= 2) {
      h2["rap1rap2"].Fill(lts.decaytree[0].p4.Rap(),
                           lts.decaytree[1].p4.Rap(), totalweight);
    }
    if (!lts.decaytree.empty() && lts.pfinal.size() >= 1) {
      FillCosThetaPhi(totalweight, lts);
    }
  }
}


// (costheta,phi) of daughter in different Lorentz frames, input as the total
// event weight
void MUserHistograms::FillCosThetaPhi(double totalweight, const gra::LORENTZSCALAR &lts) {
  if (lts.decaytree.empty() || lts.pfinal.empty()) {
    return;
  }

  // Count zero-weight trials without interpreting rejected momenta
  if (std::fpclassify(totalweight) == FP_ZERO) {
    const double missing = std::numeric_limits<double>::quiet_NaN();
    for (const auto frame : HistogramFrames()) {
      const std::string name(frame);
      h1["costheta_" + name].Fill(missing, 0.0);
      h1["phi_" + name].Fill(missing, 0.0);
      h2["costhetaphi_" + name].Fill(missing, missing, 0.0);
    }
    return;
  }

  // The angular pair must describe the same selected daughter
  const std::vector<M4Vec> pf = {lts.decaytree[0].p4};

  // System 4-momentum
  M4Vec X;
  for (const auto &i : indices(lts.decaytree)) { X += lts.decaytree[i].p4; }

  // Rejected phase space points have no rest-frame angular observables
  if (!std::isfinite(X.E()) || !(X.E() > 0.0) ||
      !std::isfinite(X.M2()) || !(X.M2() > 0.0)) {
    return;
  }

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

    const auto frames = HistogramFrames();
    for (const auto &k : indices(frames)) {
      const std::string frame(frames[k]);
      if (frame == "GJ" || frame == "LA") {
        continue;  // Treated outside this loop
      }
      
      // Transform and histogram
      std::vector<M4Vec> pfout;
      gra::kinematics::LorentzFrame(pfout, pb1boost, pb2boost, pfboost,
                                    frame, direction);
      
      const double costheta = pfout[0].CosTheta();
      const double phi      = gra::math::Rad2Deg(pfout[0].Phi());

      // 1D
      h1["costheta_" + frame].Fill(costheta, totalweight);
      h1["phi_" + frame].Fill(phi, totalweight);

      // 2D
      h2["costhetaphi_" + frame].Fill(costheta, phi, totalweight);
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
    const M4Vec &system = lts.pfinal[0];
    if (!std::isfinite(system.E()) || !(system.E() > 0.0) ||
        !std::isfinite(system.M2()) || !(system.M2() > 0.0)) {
      return;
    }
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
    for (auto &entry : h1) {
      entry.second.FlushBuffer();
      entry.second.Print();
    }
  }
  if (HIST >= 2) {
    for (auto &entry : h2) {
      entry.second.FlushBuffer();
      entry.second.Print();
    }
  }
}

// Save all histograms out
void MUserHistograms::SaveHistograms(const std::string filename) {
  if (HIST < 1) {
    return;
  }
  gra::aux::CreateDirectory(std::filesystem::path(filename).parent_path().string());
  std::ofstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error(
        "MUserHistograms::SaveHistograms: cannot open " + filename);
  }

  nlohmann::json j;

  for (auto &entry : h1) {
    entry.second.FlushBuffer();
    nlohmann::json jthis;
    entry.second.struct2json(jthis);
    j["h1"][entry.first] = jthis;
  }

  file << j << std::endl;
  if (!file.good()) {
    throw std::runtime_error(
        "MUserHistograms::SaveHistograms: failed to write " + filename);
  }
  std::cout << "MUserHistograms::SaveHistograms: JSON to " << filename
            << std::endl;
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
