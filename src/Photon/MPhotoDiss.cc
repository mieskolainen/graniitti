// HERA proton dissociation profiles for vector meson photoproduction
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Photon/MPhotoDiss.h"

#include <cmath>
#include <stdexcept>
#include <string>

#include "Graniitti/Particle/MPDG.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MException.h"

using gra::aux::indices;

namespace gra::flux {

// Read HERA rows and normalize the continuum in the reference mass interval
std::map<int, PhotoDissParam> ReadPhotoDiss(const nlohmann::json& rows) {
  if (!rows.is_array()) { throw std::invalid_argument("PARAM_REGGE.photoprod_diss must be a row array"); }
  std::map<int, PhotoDissParam> output;
  for (const auto i : indices(rows)) {
    const auto& row = rows[i];
    const std::string context = "PARAM_REGGE.photoprod_diss[" + std::to_string(i) + "]";
    if (!row.is_array() || row.size() != 8 || !row[0].is_number_integer()) {
      throw std::invalid_argument(context + " must be [PDG, W0, ratio, delta, b, n, epsilon, M_max]");
    }
    for (std::size_t k = 1; k < row.size(); ++k) {
      if (!row[k].is_number() || !std::isfinite(row[k].get<double>())) {
        throw std::invalid_argument(context + " contains an invalid number");
      }
    }
    const int pdg = row[0].get<int>();
    PhotoDissParam param{row[1], row[2], row[3], row[4], row[5], row[6], row[7]};
    const double threshold = PDG::mp + PDG::mpi;
    if (pdg <= 0 || !(param.W0 > 0.0) || param.ratio < 0.0 || !(param.delta > -4.0) ||
        !(param.b > 0.0) || !(param.n > 1.0) || param.epsilon < 0.0 || !(param.mass_max > threshold)) {
      throw std::invalid_argument(context + " contains an invalid physical parameter");
    }
    const double span = 2.0 * std::log(param.mass_max / threshold);
    param.mass_norm = threshold * threshold *
        (param.epsilon > 1.0e-12 ? -std::expm1(-param.epsilon * span) / param.epsilon : span);
    if (!std::isfinite(param.mass_norm) || !(param.mass_norm > 0.0) || !output.emplace(pdg, param).second) {
      throw std::invalid_argument(context + " has invalid mass normalization or duplicate PDG");
    }
  }
  return output;
}

// Compute the smooth mass continuum and the measured nonzero forward t profile
// [REFERENCE: H1 Collaboration, arXiv:2005.14471, Eqs. (25), (41) and Tables 8, 11]
// [REFERENCE: H1 Collaboration, arXiv:1304.5162, Sections 2.2, 3.1, 3.2 and Tables 2, 3]
double PhotoDissFactor(const PhotoDissParam& param, double w2, double t, double mass2) {
  if (!std::isfinite(w2) || !(w2 > 0.0) || !std::isfinite(t) || t > 0.0 ||
      !std::isfinite(mass2) || !(mass2 > 0.0)) {
    throw AmplitudeFailure("PhotoDissFactor: invalid generated invariants");
  }
  const double threshold2 = (PDG::mp + PDG::mpi) * (PDG::mp + PDG::mpi);
  if (mass2 <= threshold2) { return 0.0; }
  const double density = std::exp(-(1.0 + param.epsilon) * std::log(mass2 / threshold2)) / param.mass_norm;
  const double shape = std::exp(-0.5 * param.n * std::log1p(-param.b * t / param.n));
  return std::pow(w2 / (param.W0 * param.W0), param.delta / 4.0) * std::sqrt(param.ratio * density) * shape;
}

}  // namespace gra::flux
