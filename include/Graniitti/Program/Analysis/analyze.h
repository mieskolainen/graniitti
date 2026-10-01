// Fiducial observable profile initialization
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef PROGRAM_ANALYZE_H
#define PROGRAM_ANALYZE_H
#include <algorithm>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "Graniitti/Analysis/MMultiplet.h"
#include "Graniitti/Math/MMath.h"
namespace gra::program {
// Initialize Profile histograms
//
inline void InitPrHistogram(std::map<std::string, std::shared_ptr<hProfMultiplet>>& h,
                            const std::vector<std::string>& legendtext, std::vector<int> multiplicity,
                            const std::string& title, const h1Bound& bM) {
  std::string name = "null";

  // Central system observables
  name    = "hP_S_M_Pt";
  h[name] = std::make_shared<hProfMultiplet>(name, title + ";System M  (GeV); System #LTP_{T}#GT (GeV)", bM.N, bM.min,
                                             bM.max, 0.0, std::numeric_limits<double>::max(), legendtext);

  name = "hP_S_M_PL2_CM";
  h[name] =
      std::make_shared<hProfMultiplet>(name, title + ";System M  (GeV); Legendre #LTP_{l=2}(cos #theta)#GT [SR frame]",
                                       bM.N, bM.min, bM.max, -1.0, 1.0, legendtext);

  name = "hP_S_M_PL4_CM";
  h[name] =
      std::make_shared<hProfMultiplet>(name, title + ";System M  (GeV); Legendre #LTP_{l=4}(cos #theta)#GT [SR frame]",
                                       bM.N, bM.min, bM.max, -1.0, 1.0, legendtext);

  // 2-body
  if (std::find(multiplicity.begin(), multiplicity.end(), 2) != multiplicity.end()) {
    name    = "hP_2B_M_dphi";
    h[name] = std::make_shared<hProfMultiplet>(
        name, title + ";System M (GeV); Central final state pair #LT#delta#phi#GT  (rad)", bM.N, bM.min, bM.max, 0.0,
        math::PI, legendtext);
  }

  // 4-Body observables
  if (std::find(multiplicity.begin(), multiplicity.end(), 4) != multiplicity.end()) {
    // ...
  }
}

}  // namespace gra::program
#endif
