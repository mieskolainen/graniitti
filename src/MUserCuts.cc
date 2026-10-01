// Custom user cuts which cannot be implemented directly via json steering files
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <array>
#include <cmath>
#include <complex>
#include <iostream>

// Own
#include "Graniitti/Particle/MForm.h"
#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/MUserCuts.h"

namespace gra {

namespace {

// Compute whether a momentum slot contains finite positive-energy data
bool HasMomentum(const M4Vec &momentum) noexcept {
  return std::isfinite(momentum.Px()) && std::isfinite(momentum.Py()) &&
         std::isfinite(momentum.Pz()) && std::isfinite(momentum.E()) &&
         momentum.E() > 0.0;
}

// Test one outgoing proton against the STAR Roman-pot fiducial region
bool PassStarForwardProton(const M4Vec &proton) {
  const double px = proton.Px();
  const double py = proton.Py();

  return px > -0.2 && 0.2 < std::abs(py) && std::abs(py) < 0.4 &&
         math::pow2(px + 0.3) + math::pow2(py) < 0.25;
}

// Test both outgoing protons against the shared STAR Roman-pot acceptance
bool PassStarForwardAcceptance(const LORENTZSCALAR &lts) {
  return PassStarForwardProton(lts.pfinal[1]) && PassStarForwardProton(lts.pfinal[2]);
}

// Apply the STAR 510 GeV branch geometry in the STAR laboratory axes
// [REFERENCE: STAR Collaboration, arXiv:2510.27482, Eq. (4.1), Table 1 and Fig. 3]
bool PassStar510Proton(const M4Vec &proton) {
  // EU, ED, WU, WD rows: px_min, py_min, py_max, px_center, R2, px_bar, R2_bar
  constexpr std::array<std::array<double, 7>, 4> rp = {{
      {-0.2300, 0.4200, 0.8600, 0.6400, 1.3600,  0.0000, 0.0000},
      {-0.2500, 0.4800, 0.8400, 0.7000, 1.5000, -0.2500, 0.9590},
      {-0.2100, 0.4600, 0.8400, 0.6000, 1.3000, -0.2800, 0.9460},
      {-0.1900, 0.4600, 0.8800, 0.7000, 1.5000,  0.0000, 0.0000}}};
  const bool west = proton.Pz() > 0.0;
  const bool up = proton.Py() > 0.0;
  const auto &[px_min, py_min, py_max, px_center, r2, px_bar, r2_bar] =
      rp[(west ? 2 : 0) + (up ? 0 : 1)];
  const double px = proton.Px(), py = proton.Py();
  return px > px_min && std::abs(py) > py_min &&
         math::pow2(px + px_center) + math::pow2(py) < r2 &&
         (west == up ? math::pow2(px + px_bar) + math::pow2(py) < r2_bar
                     : std::abs(py) < py_max);
}

}  // namespace

// Compute whether an integer selects an implemented user cut
bool IsKnownUserCut(std::int64_t id) noexcept {
  switch (id) {
    case 0:
    case -3:
    case 3:
    case -111:
    case 111:
    case 1792394000:
    case 1792394010:
    case 1792394020:
    case 3075716000LL:
    case 3075716010LL:
    case 3075716020LL:
    case 160803765:
    case 7120604:
    case 170804053:
    case 1230123:
      return true;
    default:
      return false;
  }
}

namespace {

// Apply custom cuts with the selected central tree and mass definition
bool ApplyUserCut(std::int64_t id, const gra::LORENTZSCALAR &lts,
                  const std::vector<MDecayBranch> &central,
                  bool radiated) noexcept {
  // ** NO CUTS CASE, THIS SHOULD BE FIRST **
  if (id == 0) {
    return true;
  }

  const bool has_forward_pair =
      lts.pfinal.size() >= 3 && HasMomentum(lts.pfinal[1]) &&
      HasMomentum(lts.pfinal[2]);
  const bool has_central_pair = central.size() >= 2 &&
                                HasMomentum(central[0].p4) && HasMomentum(central[1].p4);

  // -------------------------------------------------------------------
  // "Spin-filter" cut ('Glueball' filter)

  // Forward proton |dpt| < 0.3 GeV
  if (id == -3) {
    if (!has_forward_pair) {
      return false;
    }
    const double dpt = (lts.pfinal[1] - lts.pfinal[2]).Pt();
    if (dpt < 0.3) {
      // fine
    } else {
      return false;  // did not pass
    }
  }

  // Forward proton |dpt| > 0.3 GeV
  else if (id == 3) {
    if (!has_forward_pair) {
      return false;
    }
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
    if (!has_forward_pair) {
      return false;
    }
    if (lts.pfinal[1].DotPt(lts.pfinal[2]) < 0) {
      // fine
    } else {
      return false;  // did not pass
    }
  }

  // Forward proton pt1 dot pt2 > 0
  else if (id == 111) {
    if (!has_forward_pair) {
      return false;
    }
    if (lts.pfinal[1].DotPt(lts.pfinal[2]) > 0) {
      // fine
    } else {
      return false;  // did not pass
    }
  }


  // --------------------------------------------------------------------
  // STAR/RHIC sqrt(s) = 200 GeV pi+pi- / K+K- / ppbar non-factorizable cuts
  // [REFERENCE: STAR Collaboration, arXiv:2004.11078]
  //
  else if (id == 1792394000 || id == 1792394010 || id == 1792394020) {
    if (!has_forward_pair ||
        ((id == 1792394010 || id == 1792394020) && !has_central_pair)) {
      return false;
    }
    if (!PassStarForwardAcceptance(lts)) {
      return false;
    }

    // Keep the non-factorizable K+K- transverse-momentum condition in C++
    if (id == 1792394010 &&
        std::min(central[0].p4.Pt(), central[1].p4.Pt()) >= 0.7) {
      return false;
    }

    // Keep the non-factorizable ppbar transverse-momentum condition in C++
    if (id == 1792394020 &&
        std::min(central[0].p4.Pt(), central[1].p4.Pt()) >= 1.1) {
      return false;
    }

  // STAR/RHIC sqrt(s) = 510 GeV pi+pi-, K+K- and ppbar fiducial cuts
  // [REFERENCE: STAR Collaboration, arXiv:2510.27482, Eqs. (4.1), (4.10) and (4.11)]
  } else if (id == 3075716000LL || id == 3075716010LL || id == 3075716020LL) {
    if (!has_forward_pair || !has_central_pair) { return false; }
    if (!PassStar510Proton(lts.pfinal[1]) || !PassStar510Proton(lts.pfinal[2])) { return false; }
    if (id != 3075716000LL && std::min(central[0].p4.Pt(), central[1].p4.Pt()) >=
        (id == 3075716010LL ? 0.7000 : 1.1000)) { return false; }
  
  // --------------------------------------------------------------------
  // [arxiv.org/abs/1608.03765]
  
  } else if (id == 160803765) {
    const double xi1 = lts.xi1;
    const double xi2 = lts.xi2;

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
    const double xi_pbar = lts.xi2;

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
    if (!has_central_pair) {
      return false;
    }
    const double M = radiated ? (central[0].p4 + central[1].p4).M()
                              : math::msqrt(lts.m2);

    if (12 <= M && M < 30) {  // GeV
      if (central[0].p4.Pt() > 6 && central[1].p4.Pt() > 6) {
        // fine
      } else {
        return false;  // not passed
      }
    } else if (30 <= M && M <= 70) {  // GeV
      if (central[0].p4.Pt() > 10 && central[1].p4.Pt() > 10) {
        // fine
      } else {
        return false;  // not passed
      }
    } else {
      return false;
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
    if (!has_forward_pair) {
      return false;
    }
    // Forward protons |py| and |phi|
    const std::array<double, 2> pyabs = {
        std::abs(lts.pfinal[1].Py()), std::abs(lts.pfinal[2].Py())};
    const std::array<double, 2> phiabs = {
        std::abs(lts.pfinal[1].Phi()), std::abs(lts.pfinal[2].Phi())};

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
  } else {
    return false;
  }

  return true;
}

}  // namespace

// Apply custom cuts to the event kinematics stored in the Lorentz state
bool UserCut(std::int64_t id, const gra::LORENTZSCALAR &lts) noexcept {
  return ApplyUserCut(id, lts, lts.decaytree, false);
}

// Apply custom cuts with an explicit post-radiation central momentum tree
bool UserCut(std::int64_t id, const gra::LORENTZSCALAR &lts,
             const std::vector<MDecayBranch> &central) noexcept {
  return ApplyUserCut(id, lts, central, true);
}

}  // namespace gra
