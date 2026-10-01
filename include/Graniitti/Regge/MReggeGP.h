// Analytic GP trajectory helicity functions
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MREGGEGP_H
#define MREGGEGP_H

// C++
#include <complex>
#include <cstddef>
#include <optional>
#include <string>
#include <utility>
#include <vector>

// Own
#include "Graniitti/Regge/MReggeParam.h"
#include "Graniitti/Spin/MHelicityScatter.h"

namespace gra {
namespace gpom {

// Compute the photon transition at a specified transfer
std::complex<double> PhotoCoupling(const LORENTZSCALAR& lts, const regge::Param& param, const RES_PRODUCTION& production, double t);

// Compute the complex transverse helicity-conserving GP vertex at the requested trajectory spin
std::complex<double> PhotoCoupling(const RES_PRODUCTION& production, double alpha, double momentum, bool derivative_factor);

// Key one call-local analytic source by its exchange pair
struct SourceKey {
  int MMAX   = 0;
  int up_pdg = 0;
  int dn_pdg = 0;

  bool operator==(const SourceKey& other) const = default;
};

// Store the analytic exchange basis for one amplitude and exchange pair
struct ReggeBasis {
  int                                    MMAX                          = 0;
  std::size_t                            nm                            = 0;
  bool                                   upper_photon                  = false;
  bool                                   lower_photon                  = false;
  bool                                   use_exchange_helicity_barrier = true;
  double                                 alpha1                        = 1.0;
  double                                 alpha2                        = 1.0;
  std::complex<double>                   exchange_section              = 1.0;
  std::vector<int>                       m_values;
  std::vector<double>                    upper_nonsense;
  std::vector<double>                    lower_nonsense;
  std::vector<std::pair<double, double>> upper_rows;
  std::vector<std::pair<double, double>> lower_rows;
};

// Store one dynamically sized pole reflection projected Gamma spin block
struct ReggeCGBlock {
  std::size_t                                      two_s            = 0;
  int                                              reflection_phase = 1;
  std::vector<std::optional<std::complex<double>>> coefficient;
};

// Store one amplitude-local analytic forward source
struct ForwardSource {
  SourceKey                 key;
  ReggeBasis                basis;
  HelAmp                    upper_residue;
  HelAmp                    lower_residue;
  HelVec                    upper_row_factor;
  HelVec                    upper_column_factor;
  HelVec                    lower_row_factor;
  HelVec                    lower_column_factor;
  std::vector<ReggeCGBlock> regge_cg;
};

// Share GP[CON] and GP[RES] work at one kinematic point, recreate at every screening node
struct AmpCache {
  std::vector<ForwardSource> sources;
};

// Compute the central tensor with exchange-pair rows and resonance-spin columns
HelAmp Fusion(const HELMatrix& hel, ForwardSource& source, double momentum, bool derivative_factor);

// Build per-channel analytic resonance production matrices
std::vector<HelAmp> Resonance(const LORENTZSCALAR& lts, const regge::Param& param, const PARAM_RES& res, AmpCache* amp_cache);

// Build per-channel analytic continuum production matrices
std::vector<HelPair> Continuum(const LORENTZSCALAR& lts, const regge::Param& param, AmpCache* amp_cache = nullptr);

// Compute the local GP reduced Regge contraction basis
std::vector<double> LocalBasis(const ReggeContinuumPole& pole, std::size_t leg);

// Evaluate one analytic GP local final-pair kernel
HelAmp PairKernel(const ReggeContinuumPole& pole, const regge::Param& param, const MDecayBranch& first, const MDecayBranch& second, double alpha_upper, double alpha_lower);

// Compute the analytic trajectory m basis matrix column index
std::size_t AnalyticMIndex(int m, int MMAX, const std::string& context);

// Compute the canonical pole spin for a photon or mapped analytic trajectory
int PoleSpinX2(const regge::Param& param, const MPDG& pdg_table, int pdg);

// Precompute the pole-normalized central-resonance orbital SU(2) factors
void InitResonanceLS(HELMatrix& hc, int JX, int j1x2, int j2x2);

// Precompute the pole normalization of continuous-J crossed LS vertices
void InitCrossedLS(HELMatrix& hc, int JX);

// Compute the physical photon pole helicity tensor
HELMatrix PhotonPoleHelicity(const HELMatrix& input);

// Evaluate one angle-free crossed Regge-helicity vertex
spin::EvaluatedPoleSubvertex Crossed(const HELMatrix& input, double alpha, const M4Vec& axis, bool second_exchange_daughter);

// Compute one forward Regge helicity vertex factor for an integer m label
std::complex<double> Residue(int m, double delta_lambda, double qt, double phi, double s0, bool second_exchange_daughter, bool use_exchange_helicity_barrier);

// Compute the anchor-free analytic nonsense-zero factor
double NonsenseZero(double alpha_t, int m);

// Build R_(i,m) = f_i r_m in the retained analytic helicity basis
HelAmp ResidueMatrix(const std::vector<std::pair<double, double>>& rows, const std::vector<int>& m_values, const std::vector<double>& nonsense, const M4Vec& q, double s0, bool second_exchange_daughter, bool use_exchange_helicity_barrier,
                     const std::string& context);

// Internal contractions used by the production equations in MReggeGP.cc
namespace detail {
// Share forward vertices within one kinematic point
ForwardSource& SourceFor(const LORENTZSCALAR& lts, const regge::Param& param, SourceKey key, AmpCache& cache);
// Construct the exchange basis at the physical pole
ForwardSource PoleSource(const HELMatrix& hel, bool upper_photon, bool lower_photon);
// Compute the physical pole density for spin averaging
double ResonanceScale(const HELMatrix& hel, const std::vector<MDecayBranch>& tree, double momentum, bool derivative_factor);
// Contract active fusion terms without allocating the full tensor
HelAmp ContractFusion(const HELMatrix& hel, ForwardSource& source, double momentum, bool derivative_factor);
// Contract the ordered continuum vertices
HelPair ConProd(const LORENTZSCALAR& lts, const regge::Param& param, const std::vector<HELMatrix>& vertex, const std::vector<int>& exchange_channel, const std::array<M4Vec, 2>& final_in_X, AmpCache& cache);
}  // namespace detail

}  // namespace gpom
}  // namespace gra

#endif
