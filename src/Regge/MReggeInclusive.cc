// Inclusive elastic and diffractive Regge amplitudes
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Regge/MReggeInclusive.h"

#include <cmath>
#include <complex>
#include <stdexcept>

#include "Graniitti/Math/MMath.h"
#include "Graniitti/Regge/MRegge.h"
#include "Graniitti/Regge/MReggeForm.h"
#include "Graniitti/Regge/MReggeSig.h"
#include "Graniitti/Regge/MReggeSource.h"
#include "Graniitti/Tech/MException.h"

namespace gra {

// Triple-Regge (Pomeron) limits
//
// discontinuity line (cut line)
//        |
//        |
//        |
//        .
//        .
// SD:
//
// a =gN= a a =gN= a
//     *      *
//  i   *    *  j
//       *3P*
//      k *  M_X^2
//        *
// a =====gN===== a
//
//
// DD:
//
// a =====gN===== a
//    k1  *
//        *  M_Y^2
//       *3P*
//   i  *    *  j
//      *    *
//       *3P*
//        *  M_X^2
//    k2  *
// a =====gN===== a
//
//
// CD:
//
// a =gN= a  a =gN= a
//     *      *
//  i1  *    *  j1
//       *3P*
//     k  *  M_X^2
//        *
//       *3P*
//  i2  *    *  j2
//     *      *
// a =gN= a  a =gN= a
//
//
//
// [REFERENCE: Gribov, A Reggeon Diagram Technique, Soviet JETP, 1968, jetp.ac.ru/cgi-bin/dn/e_026_02_0414.pdf]
//
// Strong coupling in the Pomeranchuk pole problem
// [REFERENCE: Gribov, Migdal, Sov. Phys. JETP 28(4), 784-795 (1968)]
// [REFERENCE: Muller, 1972]
//
// For different forms of triple Pomeron coupling (weak/strong, scalar/vector):
// [REFERENCE: Luna, Khoze, Martin, Ryskin, arxiv.org/abs/1005.4864v1]

// Evaluate one inclusive EL, SD or DD Regge matrix element
std::complex<double> MRegge::ME2(LORENTZSCALAR &lts, const MReggeInclusive mode) const {
  
  if (mode != MReggeInclusive::EL && mode != MReggeInclusive::SD && mode != MReggeInclusive::DD) {
    throw std::invalid_argument("MRegge::ME2: unknown inclusive mode");
  }
  if (!std::isfinite(lts.s) || !(lts.s > 0.0) || !std::isfinite(lts.t)) {
    throw AmplitudeFailure("MRegge::ME2: invalid generated s or t");
  }
  if (lts.t > 1.0e-12) { throw AmplitudeFailure("MRegge::ME2: generated transfer is timelike"); }
  
  const SoftExchangeId pomeron_id = param.exchanges.at(param.pomeron_trajectory).soft_exchange;
  const auto          &pomeron    = soft_model->Exchange(pomeron_id);
  
  // EL
  if (mode == MReggeInclusive::EL) {
    StoreElasticGoodWalker(lts, pomeron_id, PropagatorForExchange(lts.s, lts.t, pomeron_id));
    const std::complex<double> amplitude = math::msqrt(BornNorm(lts));
    lts.proton_good_walker.reset();
    return amplitude;
  }
  const double         alpha0 = pomeron.Alpha0();
  const double         alpha  = soft_model->Alpha(pomeron_id, lts.t);
  const auto           eta    = regge::EtaFactor(alpha, alpha0, pomeron.Signature(), soft_model->TriplePomeronEtaMode());
  std::complex<double> kernel;

  // SD
  if (mode == MReggeInclusive::SD) {
    if (lts.excite1 == lts.excite2) { throw AmplitudeFailure("MRegge::ME2(SD): exactly one proton must be excited"); }
    const double mass2 = lts.excite1 ? lts.ss[1][1] : lts.ss[2][2];
    regge::CheckMass2(mass2, "MRegge::ME2(SD)");
    kernel = eta * math::msqrt(std::pow(lts.s / mass2, 2.0 * alpha) *
                               std::pow(mass2 / param.s0, alpha0));
  
  // DD
  } else {
    const double upper_mass2 = lts.ss[1][1];
    const double lower_mass2 = lts.ss[2][2];
    regge::CheckMass2(upper_mass2, "MRegge::ME2(DD) upper");
    regge::CheckMass2(lower_mass2, "MRegge::ME2(DD) lower");
    kernel = eta * math::msqrt(std::pow((lts.s * param.s0) / (upper_mass2 * lower_mass2), 2.0 * alpha) *
                               std::pow(upper_mass2 / param.s0, alpha0) * std::pow(lower_mass2 / param.s0, alpha0));
  
  }

  kernel = static_cast<double>(pomeron.ResidueSign()) * kernel;
  StoreTripleGoodWalker(lts, mode, kernel);
  if (GoodWalkerOnly(lts, "MRegge::ME2")) { return 0.0; }
  ProjectGoodWalkerBorn(lts);
  
  return math::msqrt(BornNorm(lts));
}

} // namespace gra
