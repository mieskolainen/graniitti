// Numerical parameters for parallel Regge production
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Regge/MReggeNumerics.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include "Graniitti/Particle/MPDG.h"
#include "json.hpp"

namespace gra {

// Validate all parallel production integration parameters
void MReggeNumerics::Validate() const {
  (void)math::PolarRadialRule(loop);
  (void)math::PolarAzimuthRule(loop);
  if (!std::isfinite(impact_max) || impact_max <= 0.0 || impact_count == 0) { throw std::invalid_argument("MReggeNumerics: invalid parallel impact parameter rule"); }
  if (!std::isfinite(photo_dt) || photo_dt <= 0.0) { throw std::invalid_argument("MReggeNumerics: photo_dt must be positive"); }
}

// Build the polar transverse production rule
void MReggeNumerics::BuildParallelRule() {
  const math::PolarRule1D radial  = math::PolarRadialRule(loop);
  const math::PolarRule1D angular = math::PolarAzimuthRule(loop);

  parallel_rule.kt             = radial.node;
  parallel_rule.radial_weight  = radial.measure;
  parallel_rule.measure_weight = OuterProduct(radial.measure, angular.measure);

  std::vector<double> cos_phi(angular.node.size());
  std::vector<double> sin_phi(angular.node.size());
  std::transform(angular.node.cbegin(), angular.node.cend(), cos_phi.begin(), [](const double phi) { return std::cos(phi); });
  std::transform(angular.node.cbegin(), angular.node.cend(), sin_phi.begin(), [](const double phi) { return std::sin(phi); });
  parallel_rule.kt_x = OuterProduct(radial.node, cos_phi);
  parallel_rule.kt_y = OuterProduct(radial.node, sin_phi);
}

// Compute the run owned immutable Regge numerical block
MReggeNumericsPtr GetReggeNumerics(MModelCache &cache) {
  return cache.Get<MReggeNumerics>("regge:numerics", [&cache] {
    auto numerics = std::make_shared<MReggeNumerics>();
    numerics->ConfigureFromJson(cache.Tune().NumericsFile(), cache.Tune().Numerics().dump());
    return numerics;
  });
}

// Configure the immutable parallel production rule from NUMERICS JSON
void MReggeNumerics::ConfigureFromJson(const std::string &source, const std::string &json_text) {
  if (initialized) { throw std::logic_error("MReggeNumerics: parameters are already initialized"); }
  try {
    const nlohmann::json document    = nlohmann::json::parse(json_text);
    photo_dt                        = document.at("NUMERICS_REGGE").at("photo_dt").get<double>();
    const auto          &integration = document.at("NUMERICS_REGGE").at("LOOP_INTEGRAL");
    // Validate discrete counts before narrowing or allocating quadrature nodes
    for (const char *name : {"NumberKT", "NumberPHI", "NumberBT"}) {
      const auto &value = integration.at(name);
      if (!value.is_number_integer() || value.get<long double>() < 1.0L || value.get<long double>() > std::numeric_limits<unsigned int>::max()) {
        throw std::invalid_argument(std::string(name) + " must be a positive integer in the unsigned range");
      }
    }
    loop.radial_integrator           = integration.at("kT_integrator").get<std::string>();
    loop.azimuth_integrator          = integration.at("phi_integrator").get<std::string>();
    loop.radial_map                  = integration.at("log_kT").get<bool>() ? math::RadialMap::Log : math::RadialMap::Linear;
    loop.r_min                       = integration.at("MinKT").get<double>();
    loop.r_max                       = integration.at("MaxKT").get<double>();
    loop.radial_intervals            = integration.at("NumberKT").get<unsigned int>();
    loop.azimuth_nodes               = integration.at("NumberPHI").get<unsigned int>();
    impact_max                       = integration.at("MaxBT").get<double>() / PDG::GeV2fm;
    impact_count                     = integration.at("NumberBT").get<unsigned int>();
    Validate();
    BuildParallelRule();
    initialized = true;
  } catch (const std::exception &error) { throw std::invalid_argument("MReggeNumerics::ConfigureFromJson: error reading " + source + " NUMERICS_REGGE: " + error.what()); }
}

// Compute the polar transverse production rule
const ReggeParallelRule &MReggeNumerics::ParallelRule() const {
  if (!initialized) { throw std::logic_error("MReggeNumerics: parameters are not initialized"); }
  return parallel_rule;
}

// Compute the impact parameter range used by the triple convolution
double MReggeNumerics::ParallelImpactMax() const {
  if (!initialized) { throw std::logic_error("MReggeNumerics: parameters are not initialized"); }
  return impact_max;
}

// Compute the impact parameter node count used by the triple convolution
unsigned int MReggeNumerics::ParallelImpactCount() const {
  if (!initialized) { throw std::logic_error("MReggeNumerics: parameters are not initialized"); }
  return impact_count;
}

}  // namespace gra
