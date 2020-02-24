// Photon and other fluxes
//
// (c) 2017-2020 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <complex>
#include <random>
#include <vector>

// Own
#include "Graniitti/MFlux.h"

namespace gra {
namespace flux {

using math::abs2;
using math::msqrt;

// Apply non-collinear EPA fluxes at cross section level
// Use with full 2 -> N kinematics
double ApplyktEPAfluxes(double amp2, gra::LORENTZSCALAR& lts) {
  const double gammaflux1 = lts.excite1
                                ? gra::form::IncohFlux(lts.x1, lts.t1, lts.qt1, lts.pfinal[1].M2())
                                : form::CohFlux(lts.x1, lts.t1, lts.qt1);
  const double gammaflux2 = lts.excite2
                                ? gra::form::IncohFlux(lts.x2, lts.t2, lts.qt2, lts.pfinal[2].M2())
                                : form::CohFlux(lts.x2, lts.t2, lts.qt2);

  // const double phasespace = lts.s / lts.s_hat;    // This gives problems at low masses
  const double phasespace = 1.0 / (lts.x1 * lts.x2);  // Consistent with kt-factorization

  // Combine all
  const double fluxes = gammaflux1 * gammaflux2 * phasespace;
  const double tot    = fluxes * amp2;

  // --------------------------------------------------------------------
  // Apply fluxes to helicity amplitudes
  const double sqrt_fluxes = msqrt(fluxes);  // to "amplitude level"
  for (const auto& i : aux::indices(lts.hamp)) { lts.hamp[i] *= sqrt_fluxes; }
  // --------------------------------------------------------------------

  return tot;
}

// Apply Gamma-Gamma collinear Drees-Zeppenfeld (coherent flux) at cross section level
// Use with collinear kinematics
double ApplyDZfluxes(double amp2, gra::LORENTZSCALAR& lts) {
  // Evaluate gamma pdfs
  const double f1         = form::DZFlux(lts.x1);
  const double f2         = form::DZFlux(lts.x2);
  const double phasespace = 1.0 / (lts.x1 * lts.x2);
  const double tot        = f1 * f2 * amp2 * phasespace;

  return tot;
}

// Apply Gamma-Gamma LUX-pdf (use at \mu > 10 GeV) at cross section level
// Use with collinear kinematics
double ApplyLUXfluxes(double amp2, gra::LORENTZSCALAR& lts) {
  // @@ MULTITHREADING LOCK NEEDED FOR INITIALIZATION @@
  gra::g_mutex.lock();

  // Not created yet
  if (lts.GlobalPdfPtr == nullptr) {
    std::string pdfname = lts.LHAPDFSET;

  retry:
    try {
      lts.GlobalPdfPtr = LHAPDF::mkPDF(pdfname, 0);
      lts.pdf_trials   = 0;  // fine
    } catch (...) {
      ++lts.pdf_trials;
      std::string str = "ApplyLUXfluxes: Problem with reading LHAPDF '" + pdfname + "'";
      aux::AutoDownloadLHAPDF(pdfname);  // Try autodownload
      gra::g_mutex.unlock();             // Remember before throw, otherwise deadlock
      if (lts.pdf_trials >= 2) {         // too many failures
        throw std::invalid_argument(str);
      } else {
        goto retry;
      }
    }
  }
  gra::g_mutex.unlock();
  // @@ MULTITHREADING UNLOCK @@

  // pdf factorization scale
  const double Q2 = lts.s_hat / 4.0;

  // Evaluate gamma pdfs
  double f1 = 0.0;
  double f2 = 0.0;
  try {
    // Divide x out
    f1 = lts.GlobalPdfPtr->xfxQ2(PDG::PDG_gamma, lts.x1, Q2) / lts.x1;
    f2 = lts.GlobalPdfPtr->xfxQ2(PDG::PDG_gamma, lts.x2, Q2) / lts.x2;
  } catch (...) { throw std::invalid_argument("ApplyLUXfluxes: Failed evaluating LHAPDF"); }

  const double phasespace = 1.0 / (lts.x1 * lts.x2);
  const double tot        = f1 * f2 * amp2 * phasespace;

  return tot;
}

}  // namespace flux
}  // namespace gra