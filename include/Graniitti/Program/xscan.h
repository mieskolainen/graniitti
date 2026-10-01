// Available integrated soft cross sections for energy scans
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef PROGRAM_XSCAN_H
#define PROGRAM_XSCAN_H

#include <array>
#include <limits>

#include "Graniitti/Eikonal/MEikonal.h"

namespace gra::program {

// Mark unavailable total, elastic and inelastic cross sections explicitly
inline std::array<double, 3> ScanXS(const MEikonal& eikonal) {
  const double          missing = std::numeric_limits<double>::quiet_NaN();
  std::array<double, 3> xs{missing, missing, missing};
  if (eikonal.IsInitialized()) { eikonal.GetTotXS(xs[0], xs[1], xs[2]); }
  return xs;
}

}  // namespace gra::program
#endif
