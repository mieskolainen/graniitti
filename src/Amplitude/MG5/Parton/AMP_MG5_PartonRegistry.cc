// Generated parton MadGraph process registry
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Amplitude/MG5/Parton/AMP_MG5_PartonRegistry.h"

#include <memory>
#include <string>
#include <vector>

#include "Graniitti/Amplitude/MG5/Parton/AMP_MG5_pp_z.h"
#include "Graniitti/Amplitude/MG5/Parton/AMP_MG5_pp_zj.h"
#include "Graniitti/Amplitude/MG5/Parton/AMP_MG5_pp_jj.h"
#include "Graniitti/Amplitude/MG5/Parton/AMP_MG5_pp_w.h"

namespace gra {

// Compute generated parton families and their public process channels
std::vector<MG5ProcessInfo> PartonMG5ProcessInfos() {
  return {
      MG5ProcessInfo{"MG5_PP_Z", "Z", "MG5cards/Parton/MG5_PP_Z/param_card.dat"},
      MG5ProcessInfo{"MG5_PP_ZJ", "Zj", "MG5cards/Parton/MG5_PP_ZJ/param_card.dat"},
      MG5ProcessInfo{"MG5_PP_JJ", "jj", "MG5cards/Parton/MG5_PP_JJ/param_card.dat"},
      MG5ProcessInfo{"MG5_PP_W", "W", "MG5cards/Parton/MG5_PP_W/param_card.dat"}
  };
}

// Construct one generated parton process family
std::unique_ptr<PartonMG5Process> CreatePartonMG5Process(
    const std::string &process_family) {
  if (process_family == "MG5_PP_Z") {
    return std::make_unique<AMP_MG5_pp_z>();
  }
  if (process_family == "MG5_PP_ZJ") {
    return std::make_unique<AMP_MG5_pp_zj>();
  }
  if (process_family == "MG5_PP_JJ") {
    return std::make_unique<AMP_MG5_pp_jj>();
  }
  if (process_family == "MG5_PP_W") {
    return std::make_unique<AMP_MG5_pp_w>();
  }
  return nullptr;
}

}  // namespace gra
