// Fundamental spin bases, couplings and density matrix algebra
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MHELICITY_H
#define MHELICITY_H

// C++
#include <complex>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Own
#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/Particle/MParticle.h"
#include "Graniitti/Sampling/MRandom.h"
#include "Graniitti/Spin/MHELMatrix.h"
#include "Graniitti/Spin/MHelicityBasis.h"
#include "Graniitti/Spin/MSpin.h"

namespace gra {
namespace spin {

enum class VertexContext { Auto, CrossedBeamLeg, SubTUChannelExchange };

enum class AutoCentralCouplingMode { MinL, MinS, FlatLS, FlatHelicity };

struct DirectHelicityBasis {
  std::vector<std::vector<std::complex<double>>> basis;
  MMatrix<bool> support;
};

struct DirectHelicityPhase {
  bool found = false;
  std::complex<double> phase = 1.0;
};

// Store one direct-helicity vector expanded in the canonical JW LS basis
struct DirectHelicityLSExpansion {
  std::vector<LSCoupling> rows;
  std::vector<std::complex<double>> coefficients;
  double residual_norm2 = 0.0;
};

// Store one canonical pole helicity tensor and its evaluated local frame
struct EvaluatedPoleSubvertex {
  HELMatrix                     helicity;
  MMatrix<std::complex<double>> frame;
};

// Store one internal pole numerator in its physical helicity basis
struct InternalHelicityMetric {
  std::vector<double>           projections;
  MMatrix<std::complex<double>> matrix;
};

// Infer a daughter spin from its explicit helicity column
double SpinFromHelicityColumn(const HELMatrix &hel, std::size_t column);

// Find one row with the requested pair of helicities
std::size_t FindHelicityPairRow(const HELMatrix &hel, double lambda1, double lambda2);

// Convert crossed second-leg rows into physical helicities
MMatrix<std::complex<double>> CrossedSecondRowsToPhysical(const HELMatrix                     &hel,
                                                          const MMatrix<std::complex<double>> &aux);

// Rotate helicity columns between azimuthal spin sections
void ApplyColumnHelicityPhase(MMatrix<std::complex<double>> &f, const HELMatrix &hel, double phi, double coefficient);

// Build one reduced crossed Regge frame without an angular rotation
MMatrix<std::complex<double>> ReducedCrossedFrame(const gra::HELMatrix &hel);

// Build one angle-free crossed Regge-helicity frame with its spin section
MMatrix<std::complex<double>> ReggeCrossedFrame(const gra::HELMatrix &hel, const gra::M4Vec &axis,
                                                bool second_exchange_daughter);

// Evaluate one prepared physical-pole helicity tensor in its local frame
EvaluatedPoleSubvertex EvaluatePoleSubvertex(HELMatrix helicity, const gra::M4Vec &final_in_X,
                                             const gra::M4Vec &parent_dir_in_X, bool second_exchange_daughter);

// Evaluate one virtual subchannel vertex in the common X rest frame
MMatrix<std::complex<double>> VirtualSubchannelFrame(const gra::HELMatrix &hel, gra::M4Vec final_in_X,
                                                     const gra::M4Vec &parent_dir_in_X, bool second_exchange_daughter);

// Compute the physical spin-half helicity transitions
std::vector<std::pair<double, double>> SpinHalfTransitions(bool drop_flip);

// Compute one spherical forward helicity barrier and azimuthal factor
std::complex<double> HelicityFactor(int m, double delta_lambda, double qt,
                                    double phi, double exchange_mass,
                                    double flip_mass, bool reverse_phase,
                                    bool use_exchange_barrier);

// Compute the external helicity phase of one ordered collider beam source
std::complex<double>
ForwardHelicitySectionPhase(double lambda_in, double lambda_out,
                            double collider_transfer_azimuth,
                            bool second_exchange_daughter);

// Compute the complete external section factor of a factorized spin-half source
std::complex<double>
SpinHalfForwardHelicitySectionFactor(double lambda_in, double lambda_out,
                                     double collider_transfer_azimuth,
                                     bool second_exchange_daughter);

// Compute the Jacob-Wick phase from reversing the ordered second leg
double JacobWickSecondLegReversalPhase(double spin, double physical_helicity);

// Compute the azimuth harmonic of one ordered spin-half pair transition
constexpr int ColliderSpinHalfHelicityHarmonic(int h1_in, int h2_in, int h1_out,
                                               int h2_out) {
  if ((h1_in != -1 && h1_in != 1) || (h2_in != -1 && h2_in != 1) ||
      (h1_out != -1 && h1_out != 1) || (h2_out != -1 && h2_out != 1)) {
    throw std::invalid_argument(
        "ColliderSpinHalfHelicityHarmonic: doubled helicities must be +-1");
  }
  return (h1_in - h2_in - h1_out + h2_out) / 2;
}

// Compute the azimuth harmonic for a central spin-half two-to-two process
constexpr int ColliderSpinHalfHardHelicityHarmonic(int h1_in, int h2_in,
                                                   int h1_out, int h2_out) {
  if ((h1_in != -1 && h1_in != 1) || (h2_in != -1 && h2_in != 1) ||
      (h1_out != -1 && h1_out != 1) || (h2_out != -1 && h2_out != 1)) {
    throw std::invalid_argument(
        "ColliderSpinHalfHardHelicityHarmonic: doubled helicities must be "
        "+-1");
  }
  return (h1_in - h2_in - h1_out - h2_out) / 2;
}

// Compute the signed reciprocity factor of one spin-half pair transition
constexpr double ColliderSpinHalfReciprocitySign(int h1_in, int h2_in,
                                                 int h1_out, int h2_out) {
  const int harmonic =
      ColliderSpinHalfHelicityHarmonic(h1_in, h2_in, h1_out, h2_out);
  return ((harmonic < 0 ? -harmonic : harmonic) % 2 == 0) ? 1.0 : -1.0;
}

// Compute a validated direct helicity basis matrix
MMatrix<std::complex<double>> DirectFrame(const gra::HELMatrix &hel,
                                          const std::string &context);

// Compute physical final-state helicities, excluding longitudinal real photons
// and gluons
std::vector<double> FinalStateHelicities(const gra::MParticle &particle,
                                         const std::string &context);

// Compute the physical final-state helicity dimension of one particle
std::size_t FinalStateHelicityCount(const gra::MParticle &particle,
                                    const std::string &context);

// Compute the compact physical final-state index of one helicity
std::size_t FinalStateHelicityIndex(const gra::MParticle &particle,
                                    double helicity,
                                    const std::string &context);

// Initialize one two-body helicity basis and its Jacob-Wick row cache
void InitTwoBodyBasis(gra::HELMatrix &hel, double J, double s1, double s2,
                      const std::vector<double> &parent_helicities,
                      const std::vector<double> &helicities1,
                      const std::vector<double> &helicities2,
                      const std::string &context);

// Compute a compact integer or half-integer spin label for diagnostics
std::string FormatSpinLabel(double value);

// Precompute the finite Jacob-Wick helicity row map
void InitJWRotation(gra::HELMatrix &hel);

MMatrix<std::complex<double>> fDecayMatrix(const gra::HELMatrix &hel,
                                           double theta, double phi);

// Compute the pole-normalized helicity matrix with LS threshold factors
MMatrix<std::complex<double>>
DecayLSHelicityMatrix(const gra::HELMatrix &hel, double q_ratio, bool barrier);

// Compute the angle-integrated LS threshold intensity
double DecayLSIntensity(const gra::HELMatrix &hel, double q_ratio,
                        bool barrier);

void InitTMatrix(gra::HELMatrix &hc, const gra::MParticle &p,
                 const gra::MParticle &p1, const gra::MParticle &p2,
                 bool production_mode = false,
                 const std::string &source_name = "",
                 bool require_complete_l2s = true, bool verbose = true,
                 const std::string &verbose_label = "",
                 VertexContext context = VertexContext::Auto);

// Validate direct reduced helicity amplitudes against the allowed Jacob-Wick
// subspace
void ValidateDirectTMatrix(gra::HELMatrix &hc, const gra::MParticle &p,
                           const gra::MParticle &p1, const gra::MParticle &p2,
                           bool production_mode = false,
                           const std::string &source_name = "",
                           bool require_complete_helicity = true,
                           bool verbose = true,
                           const std::string &verbose_label = "",
                           VertexContext context = VertexContext::Auto);

// Compute true when one two-body LS row passes the requested C-parity constraint
bool TwoBodyCParityAllowed(const gra::MParticle &mother,
                           const gra::MParticle &p1, const gra::MParticle &p2,
                           std::size_t l, double s, bool c_conservation);

// Check if an integer stores a defined C-parity eigenvalue
bool IsValidCParity(int C);

// Compute whether the two physical legs are identical in the current vertex
// context
bool PhysicalIdenticalPair(const gra::MParticle &p1, const gra::MParticle &p2,
                           bool production_mode, VertexContext context);

// Compute the flattened row-major direct-helicity coordinate index
std::size_t DirectHelicityCoordinateIndex(double lambda1, double lambda2,
                                          double s1, double s2,
                                          const std::string &context);

// Compute the direct-helicity leg-exchange phase for swapped table lookup
double DirectHelicityLegExchangePhase(const gra::MParticle &mother,
                                      const gra::MParticle &leg1,
                                      const gra::MParticle &leg2,
                                      const std::string &context);

// Compute the active LS rows in minimal orbital-angular-momentum order
std::vector<LSCoupling>
AllowedLSCouplings(const gra::HELMatrix &hc, const gra::MParticle &p,
                   const gra::MParticle &p1, const gra::MParticle &p2,
                   bool production_mode, const std::string &source_name,
                   VertexContext context);

// Normalize one non-zero direct helicity vector in place
void NormalizeDirectHelicityVector(std::vector<std::complex<double>> &v);

// Fix an arbitrary direct-vector phase by making the first non-zero element
// real positive
void CanonicalizeDirectHelicityVectorPhase(
    std::vector<std::complex<double>> &v);

// Add one vector to an orthonormal direct-helicity basis
void AddDirectHelicityBasisVector(
    std::vector<std::vector<std::complex<double>>> &basis,
    const std::vector<std::complex<double>> &candidate, double tol = 1e-12);

// Compute the orthonormal JW direct-helicity basis generated by all allowed LS
// rows
DirectHelicityBasis
BuildJWDirectHelicityBasis(const gra::HELMatrix &hc, const gra::MParticle &p,
                           const gra::MParticle &p1, const gra::MParticle &p2,
                           bool production_mode, const std::string &source_name,
                           const std::string &verbose_label,
                           VertexContext context);

// Compute the raw LS-generated direct-helicity vectors before orthonormalization
std::vector<std::vector<std::complex<double>>> BuildJWDirectHelicityRawBasis(
    const gra::HELMatrix &hc, const gra::MParticle &p, const gra::MParticle &p1,
    const gra::MParticle &p2, bool production_mode,
    const std::string &source_name, const std::string &verbose_label,
    VertexContext context);

// Convert one direct-helicity matrix into its complete JW LS coefficients
DirectHelicityLSExpansion DirectHelicityToLSCoefficients(
    const gra::HELMatrix &helicity, const gra::MParticle &mother,
    const gra::MParticle &leg1, const gra::MParticle &leg2,
    bool production_mode, VertexContext context);

// Compute the fixed phase relating two direct-helicity coordinates in a basis
DirectHelicityPhase DirectHelicityCoordinatePhase(
    const std::vector<std::vector<std::complex<double>>> &basis_vectors,
    std::size_t from_index, std::size_t to_index, double tolerance = 1e-9);

// Project a direct helicity vector into one orthonormal basis
std::vector<std::complex<double>> ProjectDirectHelicityVector(
    const std::vector<std::vector<std::complex<double>>> &basis,
    const std::vector<std::complex<double>> &seed);

// Build a flat direct-helicity target projected into an allowed JW basis
std::vector<std::complex<double>>
BuildFlatDirectHelicityVector(const DirectHelicityBasis &basis, double J,
                              double s1, double s2, const std::string &context);

// Store a flattened direct-helicity vector into a prepared HELMatrix
void StoreDirectHelicityVector(gra::HELMatrix &hc,
                               const std::vector<std::complex<double>> &vector);

// Build one automatic two-body central coupling in LS or direct-helicity basis
gra::HELMatrix BuildAutomaticCentralCoupling(
    const gra::MParticle &mother,
    const std::vector<gra::MParticle> &vertex_legs,
    AutoCentralCouplingMode mode, bool production_mode,
    const std::string &source_name, const std::string &verbose_label,
    bool verbose_output, bool p_symmetry, bool c_symmetry,
    VertexContext context);

// Spin-Statistics
bool BoseSymmetry(int l, int s);
bool FermiSymmetry(int l, int s);

// Density matrix functions
// Compute a source density averaged over the incoming spin states
double SourceSpinAveragedDensity(const MMatrix<std::complex<double>> &matrix,
                                 std::size_t incoming_spin_states,
                                 const std::string &context);
// Compute the spin-averaged density of one source helicity column
double ForwardSourceColumnSpinAveragedDensity(
    const MMatrix<std::complex<double>> &matrix, std::size_t column,
    std::size_t incoming_spin_states, const std::string &context);
bool Positivity(const MMatrix<std::complex<double>> &rho, double J);
// Decompose a density matrix into weighted pure-state spin vectors
std::vector<std::vector<std::complex<double>>>
SpectralPureStates(const MMatrix<std::complex<double>> &rho);
// Generate a random density matrix for the supplied doubled spin label
MMatrix<std::complex<double>> RandomRho(int spinX2, bool parity, MRandom &rng);
double VonNeumannEntropy(const MMatrix<std::complex<double>> &rho);

} // namespace spin
} // namespace gra

#endif
