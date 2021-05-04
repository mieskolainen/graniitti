// Custom user cuts which cannot be implemented directly via json steering files
//
// (c) 2017-2021 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <complex>
#include <iostream>
#include <vector>

// Own
#include "Graniitti/MForm.h"
#include "Graniitti/MKinematics.h"
#include "Graniitti/MMath.h"
#include "Graniitti/MUserCuts.h"

namespace gra {

// USERCUTS (implement custom cuts here; cuts which cannot be implemented
// in .json steering file). Label these by unique integer.

bool UserCut(int id, const gra::LORENTZSCALAR &lts) {
  // ** NO CUTS CASE, THIS SHOULD BE FIRST **
  if (id == 0) {
    return true;
  }

  // -------------------------------------------------------------------
  // "Spin-filter" cut ('Glueball' filter)

  // Forward proton |dpt| < 0.3 GeV
  else if (id == -3) {
    const double dpt = (lts.pfinal[1] - lts.pfinal[2]).Pt();
    if (dpt < 0.3) {
      // fine
    } else {
      return false;  // did not pass
    }
  }

  // Forward proton |dpt| > 0.3 GeV
  else if (id == 3) {
    const double dpt = (lts.pfinal[1] - lts.pfinal[2]).Pt();
    if (dpt > 0.3) {
      // fine
    } else {
      return false;  // did not pass
    }
  }

  // -------------------------------------------------------------------
  // "Spin-filter" cut

  // Forward proton pt1 dot pt2 < 0
  else if (id == -111) {
    if (lts.pfinal[1].DotPt(lts.pfinal[2]) < 0) {
      // fine
    } else {
      return false;  // did not pass
    }
  }

  // Forward proton pt1 dot pt2 > 0
  else if (id == 111) {
    if (lts.pfinal[1].DotPt(lts.pfinal[2]) > 0) {
      // fine
    } else {
      return false;  // did not pass
    }
  }


  // -------------------------------------------------------------------
  // "Spin-filter" cut

  // Forward proton |deltaphi| in (90, 180]
  else if (id == 90180) {
    const double deltaphiabs = lts.pfinal[1].DeltaPhiAbs(lts.pfinal[2]);
    if (math::Deg2Rad(90) < deltaphiabs && deltaphiabs <= math::Deg2Rad(180)) {
      // fine
    } else {
      return false;  // did not pass
    }
  }

  // Forward proton |deltaphi| in (0, 90]
  else if (id == 90) {
    const double deltaphiabs = lts.pfinal[1].DeltaPhiAbs(lts.pfinal[2]);
    if (0 < deltaphiabs && deltaphiabs <= math::Deg2Rad(90)) {
      // fine
    } else {
      return false;  // did not pass
    }
  }

  // --------------------------------------------------------------------
  // STAR/RHIC \sqrt{s} = 200 GeV pi+pi- / K+K-, ppbar (+ other cuts needed in .json file)
  // https://rivet.hepforge.org/analyses/STAR_2020_I1792394
  //
  else if (id == 1792394000 || id == 1792394001 || id == 1792394002 || id == 1792394010 ||
           id == 1792394011 || id == 1792394012 || id == 1792394020 || id == 1792394021 ||
           id == 1792394022) {
    const std::vector<int> indices = {1, 2};

    // Loop over forward protons
    for (const auto &i : indices) {
      if (lts.pfinal[i].Px() > -0.2 && std::abs(lts.pfinal[i].Py()) > 0.2 &&
          std::abs(lts.pfinal[i].Py()) < 0.4 &&
          (math::pow2(lts.pfinal[i].Px() + 0.3) + math::pow2(lts.pfinal[i].Py())) < 0.25) {
        // fine
      } else {
        return false;  // not passed
      }

      // ============================================================
      // *** Extra forward cuts for all pi+pi- / K+K- / ppbar ***
      if (id == 1792394001 || id == 1792394011 || id == 1792394021) {
        const double deltaphiabs = lts.pfinal[1].DeltaPhiAbs(lts.pfinal[2]);
        if (0 < deltaphiabs && deltaphiabs <= math::Deg2Rad(90)) {
          // fine
        } else {
          return false;  // did not pass
        }
      } else if (id == 1792394002 || id == 1792394012 || id == 1792394022) {
        const double deltaphiabs = lts.pfinal[1].DeltaPhiAbs(lts.pfinal[2]);
        if (math::Deg2Rad(90) < deltaphiabs && deltaphiabs <= math::Deg2Rad(180)) {
          // fine
        } else {
          return false;  // did not pass
        }
      }

      // ============================================================
      // ** Extra central K+, K- cuts ***
      if (id == 1792394010 || id == 1792394011 || id == 1792394012) {
        if (std::min(lts.decaytree[0].p4.Pt(), lts.decaytree[1].p4.Pt()) < 0.7) {
          // fine
        } else {
          return false;  // did not pass
        }
      }

      // ============================================================
      // ** Extra central p, pbar cuts ***
      if (id == 1792394020 || id == 1792394021 || id == 1792394022) {
        if (std::min(lts.decaytree[0].p4.Pt(), lts.decaytree[1].p4.Pt()) < 1.1) {
          // fine
        } else {
          return false;  // did not pass
        }
      }
    }

    // --------------------------------------------------------------------
    // [arxiv.org/abs/1608.03765]

  } else if (id == 160803765) {
    const double xi1 = (lts.pbeam1.Pz() - lts.pfinal[1].Pz()) / lts.pbeam1.Pz();
    const double xi2 = (lts.pbeam2.Pz() - lts.pfinal[2].Pz()) / lts.pbeam2.Pz();

    const double XI_MAX = 0.03;

    if (xi1 < XI_MAX && xi2 < XI_MAX) {
      // fine
    } else {
      return false;  // not passed
    }
  }

  // --------------------------------------------------------------------
  // CDF exclusive dijets
  // [arxiv.org/abs/0712.0604]

  else if (id == 7120604) {
    // antiproton longitudinal momentum loss fraction
    const double xi_pbar = (lts.pbeam2.Pz() - lts.pfinal[2].Pz()) / lts.pbeam2.Pz();

    if (0.03 < xi_pbar && xi_pbar < 0.08) {
      // fine
    } else {
      return false;  // not passed
    }
  }

  // --------------------------------------------------------------------
  // ATLAS yy->mu+mu- 13 TeV fiducial cuts (+other cuts needed in .json file)
  // [arxiv.org/abs/hep-ex/170804053]

  else if (id == 170804053) {
    const double M = math::msqrt(lts.m2);

    if (12 <= M && M < 30) {  // GeV
      if (lts.decaytree[0].p4.Pt() > 6 && lts.decaytree[1].p4.Pt() > 6) {
        // fine
      } else {
        return false;  // not passed
      }
    } else if (30 <= M && M <= 70) {  // GeV
      if (lts.decaytree[0].p4.Pt() > 10 && lts.decaytree[1].p4.Pt() > 10) {
        // fine
      } else {
        return false;  // not passed
      }
    }
  }

  // --------------------------------------------------------------------
  // ATLAS pi+pi- 13 TeV roman pot fiducial cuts
  // (+other cuts needed in .json file)
  // from R. Sikora, ATLAS poster, Bad Honnef QCD School 2017
  //
  // (N.B. check the implementation
  // c.f. cut |t| > 0.03 GeV^2 seems to give more physical results)

  else if (id == 1230123) {
    // Forward protons |py| and |phi|
    const std::vector<double> pyabs  = {std::abs(lts.pfinal[1].Py()), std::abs(lts.pfinal[2].Py())};
    const std::vector<double> phiabs = {std::abs(lts.pfinal[1].Phi()),
                                        std::abs(lts.pfinal[2].Phi())};

    for (const auto &i : aux::indices(pyabs)) {
      if ((0.17 < pyabs[i]) && (pyabs[i] < 0.5)) {  // GeV
                                                    // fine
      } else {
        return false;  // not passed
      }
      if ((math::PI / 4 < phiabs[i]) && (phiabs[i] < 3.0 * math::PI / 4)) {
        // fine
      } else {
        return false;  // not passed
      }
    }

    // Roman pot geometry
    const double deltaphiabs = lts.pfinal[1].DeltaPhiAbs(lts.pfinal[2]);
    if ((deltaphiabs < math::Deg2Rad(40.0)) || (deltaphiabs > math::Deg2Rad(140.0))) {
      // fine
    } else {
      return false;  // not passed
    }
  } else {  // Poor input
    std::string str = "MUserCuts:: Unknown cut ID: " + std::to_string(id);
    throw std::invalid_argument(str);
  }

  return true;
}

}  // namespace gra
