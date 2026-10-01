// Numerical parameters for parallel Regge production
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MREGGENUMERICS_H
#define MREGGENUMERICS_H

#include <memory>
#include <string>
#include <vector>

#include "Graniitti/MModelCache.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Math/MPolarQuadrature.h"

namespace gra {

// Store the polar transverse rule used by simultaneous Regge ladders
struct ReggeParallelRule {
  std::vector<double> kt;
  std::vector<double> radial_weight;
  MMatrix<double>     measure_weight;
  MMatrix<double>     kt_x;
  MMatrix<double>     kt_y;
};

// Own numerical steering for simultaneous Regge production
class MReggeNumerics {
 public:
  // Configure the immutable parallel production rule from NUMERICS JSON
  void ConfigureFromJson(const std::string &source, const std::string &json_text);

  // Compute the polar transverse production rule
  const ReggeParallelRule &ParallelRule() const;

  // Compute the impact parameter range used by the triple convolution
  double ParallelImpactMax() const;

  // Compute the impact parameter node count used by the triple convolution
  unsigned int ParallelImpactCount() const;

  // Compute the transfer step for the elementary photoproduction slope
  double PhotoStep() const { return photo_dt; }

 private:
  math::PolarParam  loop;
  double            impact_max   = 0.0;
  unsigned int      impact_count = 0;
  double            photo_dt     = 0.0;
  bool              initialized  = false;
  ReggeParallelRule parallel_rule;

  // Validate all parallel production integration parameters
  void Validate() const;

  // Build the polar transverse production rule
  void BuildParallelRule();
};

using MReggeNumericsPtr = std::shared_ptr<const MReggeNumerics>;

// Compute the run owned immutable Regge numerical block
MReggeNumericsPtr GetReggeNumerics(MModelCache &cache);

}  // namespace gra

#endif
