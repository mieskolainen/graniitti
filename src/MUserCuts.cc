// Custom user cuts which cannot be implemented directly via json steering files
//
// (c) 2017-2019 Mikael Mieskolainen
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

  // Forward proton |deltaphi| in (90, 180]
  else if (id == 90180) {
    const double deltaphiabs = lts.pfinal[1].DeltaPhiAbs(lts.pfinal[2]);
    if (gra::math::Deg2Rad(90) < deltaphiabs && deltaphiabs <= gra::math::Deg2Rad(180)) {
      // fine
    } else {
      return false;  // did not pass
    }
  }

  // Forward proton |deltaphi| in (0, 90]
  else if (id == 90) {
    const double deltaphiabs = lts.pfinal[1].DeltaPhiAbs(lts.pfinal[2]);
    if (0 < deltaphiabs && deltaphiabs <= gra::math::Deg2Rad(90)) {
      // fine
    } else {
      return false;  // did not pass
    }
  }

  // --------------------------------------------------------------------
  // STAR/RHIC \sqrt{s} = 200 GeV pi+pi- (+ other cuts needed in .json file)
  // https://indico.cern.ch/event/713101/contributions/3102315/attachments/1705771/2748440/Diffraction2018_RafalSikora.pdf
  else if (id == 280818) {
    
    // Loop over forward protons
    for (std::size_t i = 0; i < 2; ++i) {
      if (gra::math::pow2(lts.pfinal[i].Px() + 0.3) + gra::math::pow2(lts.pfinal[i].Py()) < 0.25) {  // GeV^2
                                                                             // fine
      } else {
        return false;  // not passed
      }

      if (0.2 < std::abs(lts.pfinal[i].Py()) && std::abs(lts.pfinal[i].Py()) < 0.4) {  // GeV
                                                             // fine
      } else {
        return false;  // not passed
      }

      if (lts.pfinal[i].Px() > -0.2) {  // GeV
                           // fine
      } else {
        return false;  // not passed
      }
    }

    // https://arxiv.org/pdf/1608.03765.pdf
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
  // https://arxiv.org/pdf/0712.0604.pdf

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
  // https://arxiv.org/abs/hep-ex/170804053

  else if (id == 170804053) {
    const double M = gra::math::msqrt(lts.m2);

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
  // ATLAS pi+pi- 13 TeV roman pot fiducial cuts (+other cuts needed in .json
  // file)
  // from R. Sikora, ATLAS poster, Bad Honnef QCD School 2017
  //
  // (not sure if correct here/implemented right way
  // c.f. cut |t| > 0.03 GeV^2 seems to give more physical results)

  else if (id == 1230123) {
    // Forward protons |py| and |phi|
    std::vector<double> pyabs  = {std::abs(lts.pfinal[1].Py()), std::abs(lts.pfinal[2].Py())};
    std::vector<double> phiabs = {std::abs(lts.pfinal[1].Phi()), std::abs(lts.pfinal[2].Phi())};
    for (std::size_t i = 0; i < 2; ++i) {
      if ((0.17 < pyabs[i]) && (pyabs[i] < 0.5)) {  // GeV
                                                    // fine
      } else {
        return false;  // not passed
      }
      if ((gra::math::PI / 4 < phiabs[i]) && (phiabs[i] < 3.0 * gra::math::PI / 4)) {
        // fine
      } else {
        return false;  // not passed
      }
    }

    // Roman pot geometry
    const double deltaphiabs = lts.pfinal[1].DeltaPhiAbs(lts.pfinal[2]);
    if ((deltaphiabs < gra::math::Deg2Rad(40.0)) || (deltaphiabs > gra::math::Deg2Rad(140.0))) {
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

}  // gra namespace ends
