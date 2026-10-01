// Shared physical-pole LS algebra and immutable continuum vertices
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MPOLELS_H
#define MPOLELS_H

// C++
#include <complex>
#include <cstddef>
#include <memory>
#include <vector>

// Own
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Particle/MParticle.h"
#include "Graniitti/Spin/MHELMatrix.h"
#include "Graniitti/Spin/MHelicity.h"

namespace gra::spin {

// Describe one complete canonical physical-pole LS operator
struct CanonicalPoleOperator {
  LSCoupling           coupling;
  double               raw_normalization  = 0.0;
  std::complex<double> leg_exchange_phase = 1.0;
};

// Store physical-pole LS coefficients and their fixed helicity recoupling data
// The production model separately selects the angular and momentum continuation
// Photon flags denote a transverse pole projection, not a massive spin-one
// field or a complete gauge-invariant photon operator
struct PoleLS {
  HELMatrix                               helicity;
  std::vector<LSTerm>                     terms;
  std::vector<double>                     raw_normalization;
  std::vector<HelAmp>                     reduced_basis;
  // Pair M with -M through the spherical spin metric
  std::vector<double>                     dual_spin;
  std::shared_ptr<const wigner::Rotation> rotation;
  std::vector<double>                     spin_metric;
  std::size_t                             rows                 = 0;
  std::size_t                             cols                 = 0;
  double                                  Lambda               = 1.0;
  VertexContext                           context              = VertexContext::Auto;
  bool                                    ready                = false;
  bool                                    derivative_factor    = true;
  bool                                    leg1_transverse_pole = false;
  bool                                    leg2_transverse_pole = false;
};

// Compute the raw Cartesian STF coupling relative to normalized SU(2) CG
// coefficients
double RawSTFCouplingNormalization(std::size_t rank1, std::size_t rank2, std::size_t output_rank);

// Compute the complete raw STF normalization relative to the JW LS matrix
double RawLSOperatorNormalization(std::size_t j1, std::size_t j2, std::size_t J, std::size_t l, std::size_t s);

// Compute the canonical pole normalization for boson or spin-half leg pairs
double RawPoleLSNormalization(int j1x2, int j2x2, int Jx2, std::size_t l, int two_s);

// Enumerate every independent canonical physical-pole operator
std::vector<CanonicalPoleOperator> CanonicalPoleOperators(const MParticle& mother, const MParticle& leg1,
                                                          const MParticle& leg2, bool production_mode = true,
                                                          bool C_symmetry = true, bool P_symmetry = true,
                                                          VertexContext context = VertexContext::Auto);

// Validate a complete canonical pole coefficient table
void ValidateCanonicalPoleTerms(const MParticle& mother, const MParticle& leg1, const MParticle& leg2,
                                const std::vector<LSTerm>& terms, bool production_mode = true, bool C_symmetry = true,
                                bool P_symmetry = true, VertexContext context = VertexContext::Auto);

// Apply the canonical raw pole normalization to LS coefficients
void ApplyCanonicalPoleLSCoefficients(HELMatrix& helicity, const MParticle& mother, const MParticle& leg1,
                                      const MParticle& leg2);

// Prepare immutable pole-normalized STF LS matrices and helicity metadata
PoleLS PreparePoleLS(const MParticle& mother, const MParticle& leg1, const MParticle& leg2,
                     const std::vector<LSTerm>& terms, double Lambda, bool production_mode = true,
                     bool C_symmetry = true, bool P_symmetry = true, VertexContext context = VertexContext::Auto,
                     double coupling_min = 0.0, bool derivative_factor = true);

// Average the squared coherent pole amplitudes over the lowest active helicity sector
double LeadingPoleDensity(const PoleLS& vertex, double momentum);

// Sum canonical operators using timelike or crossed-vertex kinematics
HelAmp PoleLSReduced(const PoleLS& vertex, double momentum);

// Evaluate the raw LS operators into one reduced-helicity matrix cache
HELMatrix PoleLSHelicity(const PoleLS& vertex, double momentum);

// Evaluate one fixed pole as a reduced crossed Regge vertex
EvaluatedPoleSubvertex EvaluateCrossedPole(const PoleLS& vertex);

// Evaluate one canonical fixed-spin pole subvertex in its local frame
EvaluatedPoleSubvertex EvaluatePoleSubvertex(const PoleLS& vertex, const M4Vec& final_in_X,
                                             const M4Vec& parent_dir_in_X, bool second_exchange_daughter);

// Store a fixed continuum vertex with its reduced helicity amplitudes
// Coupling changes construct a new value, while worker copies share only const data
class PoleResidue {
 public:
  // Construct an unused slot for a process without a fixed-spin continuum
  PoleResidue() = default;

  // Prepare the pole density and reduced vertex amplitude before event generation
  explicit PoleResidue(const PoleLS& vertex);

  // Compute whether the continuum vertex has been prepared
  bool Ready() const { return data != nullptr; }

  // Access the immutable LS vertex used to prepare the helicity amplitudes
  const PoleLS& Pole() const { return data->pole; }

  // Access the reduced Regge frame or the physical photon helicity tensor
  const EvaluatedPoleSubvertex& Reduced() const { return data->reduced; }

  // Access the mean squared leading helicity amplitude at the fixed momentum Lambda
  double Density() const { return data->density; }

 private:
  // Keep every prepared value in the same immutable allocation as its couplings
  struct Data {
    PoleLS                 pole;
    EvaluatedPoleSubvertex reduced;
    double                 density = 0.0;
  };

  std::shared_ptr<const Data> data;
};

}  // namespace gra::spin

#endif
