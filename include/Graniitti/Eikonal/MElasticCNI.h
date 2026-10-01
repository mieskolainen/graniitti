// Elastic Coulomb and nuclear interference runtime
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MELASTIC_CNI_H
#define MELASTIC_CNI_H

// C++
#include <complex>
#include <memory>
#include <vector>

// Own
#include "Graniitti/Kinematics/M4Vec.h"
#include "Graniitti/Eikonal/MEikonalHelicity.h"
#include "Graniitti/Eikonal/MEikonalMatrix.h"
#include "Graniitti/Particle/MForm.h"
#include "Graniitti/Particle/MParticle.h"
#include "Graniitti/Regge/MSoftModel.h"

namespace gra {

// Numerical controls for the finite Coulomb and nuclear interference table
struct ElasticCNINumerics {
  unsigned int NumberKT2 = 2048;
  bool logKT2 = true;
  double FBIntegralMaxKT = 15.0;
  unsigned int FBIntegralN = 8000;
  double tail_max_qb = 2048.0;
  unsigned int TailIntegralN = 8192;
  int tail_order = 0;
  double tail_split_qb = 0.0;
  double interp_rel_tol = 1.0e-5;
};

// Immutable physical elastic Coulomb and nuclear interference runtime
class MElasticCNI {
public:
  using Complex = std::complex<double>;
  using Matrix = MMatrix<Complex>;

  // Construct or load the pole-regular total physical helicity table
  static std::shared_ptr<const MElasticCNI>
  Build(const MEikonalMatrix &strong_runtime, const SoftModelPtr &model,
        double s, const std::vector<MParticle> &initialstate,
        const ElasticCNINumerics &numerics, bool strong_log_b, double min_abs_t,
        double max_abs_t, const form::ParamStore &structure);

  // Compute the complete physical helicity amplitude for one elastic event
  ProtonHelicityMatrix PhysicalHelicityMatrix(const M4Vec &p1_in,
                                              const M4Vec &p2_in,
                                              const M4Vec &p1_out,
                                              const M4Vec &p2_out) const;

  // Compute the analytic pure point-Coulomb correction beyond one photon
  ProtonHelicityMatrix PointHigherOrderHelicityMatrix(double abs_t) const;

  // Compute the symmetric residual electromagnetic sandwich of the strong S
  static Matrix SymmetricShortRangeSMatrix(const Matrix &strong_s,
                                           const Matrix &chi_residual);

  // Compute the finite impact-parameter remainder after both Born subtractions
  static Matrix FiniteRemainderProfile(const Matrix &strong_s,
                                       const Matrix &chi_residual,
                                       double chi_point);

  // Compute the IR-renormalized analytic point-Coulomb phase ratio
  static Complex RenormalizedPointCoulombPhase(double eta,
                                               double point_ir_scale,
                                               double abs_t);

  // Compute the lower physical momentum-transfer bound
  double MinimumAbsT() const noexcept { return min_abs_t_; }

  // Compute the upper physical momentum-transfer bound
  double MaximumAbsT() const noexcept { return max_abs_t_; }

private:
  // Construct an empty mutable runtime before immutable publication
  MElasticCNI(double s, std::vector<MParticle> initialstate, double min_abs_t,
              double max_abs_t, double point_ir_scale, double sqrt_lambda,
              double eta);

  double s_ = 0.0;
  std::vector<MParticle> initialstate_;
  double min_abs_t_ = 0.0;
  double max_abs_t_ = 0.0;
  double point_ir_scale_ = 1.0;
  double eta_ = 0.0;
  double point_born_coefficient_ = 0.0;
  std::vector<double> q2_node_;
  std::vector<ProtonHelicityMatrix> physical_total_residue_;
};

} // namespace gra

#endif
