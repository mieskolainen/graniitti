// Resonance parameters and card parsing
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MRESONANCE_H
#define MRESONANCE_H

// C++
#include <algorithm>
#include <array>
#include <complex>
#include <cstddef>
#include <optional>
#include <string>
#include <vector>

// Own
#include "Graniitti/MGlobals.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Particle/MParticle.h"
#include "Graniitti/Regge/MFormFactor.h"
#include "Graniitti/Regge/MReggeProductionModel.h"
#include "Graniitti/Sampling/MRandom.h"
#include "Graniitti/Spin/MPoleLS.h"
#include "Graniitti/Spin/MHELMatrix.h"

// Libraries

namespace gra {

// Explicit resonance model form factors
struct RES_MODEL_FORM {
  regge::FFParam ff_transfer;
  regge::FFParam ff_prod;
};

// One fully parsed MP, XP or GP resonance production channel
struct RES_PRODUCTION_CHANNEL {
  std::array<int, 2>                 exchange      = {0, 0};
  bool                               width_derived = false;
  ReggeVertexBasis                   basis         = ReggeVertexBasis::LS;
  double                             Lambda        = 1.0;
  bool                               C_symmetry    = true;
  bool                               P_symmetry    = true;
  std::complex<double>               g             = 0.0;
  spin::LSCoefficients               g_ls;
  std::vector<std::array<double, 2>> helicity;
  std::vector<std::complex<double>>  g_helicity;

  // Compute whether the selected basis has a nonzero production coupling
  bool Active(double minimum) const {
    if (width_derived) { return true; }
    if (IsAutomaticReggeVertexBasis(basis)) { return std::abs(g) > minimum; }
    if (basis == ReggeVertexBasis::LS) {
      return std::any_of(g_ls.cbegin(), g_ls.cend(), [minimum](const auto &term) {
        return std::abs(term.coefficient) > minimum;
      });
    }
    return std::any_of(g_helicity.cbegin(), g_helicity.cend(), [minimum](const auto value) {
      return std::abs(value) > minimum;
    });
  }
};

// One typed MP, XP or GP resonance production model
struct RES_PRODUCTION_MODEL : public RES_MODEL_FORM {
  double                              phi = 0.0;
  std::vector<RES_PRODUCTION_CHANNEL> channels;
};

struct RES_MP_PRODUCTION_MODEL : public RES_PRODUCTION_MODEL {
  bool random_rho = false;
  MMatrix<std::complex<double>> filter;  // S^T acts on production spin columns
};

struct RES_TENSOR_VMD_MIXING {
  int                      pdg      = 0;
  std::array<double, 2>    coupling = {0.0, 0.0};
  std::vector<std::size_t> active_coupling;
};

// One unordered covariant rank-two exchange pair and central coupling
struct RES_TENSOR_CHANNEL {
  std::array<int, 2>                 exchange = {0, 0};
  std::vector<double>                g_tensor;
  std::vector<std::size_t>           active_g_tensor;
  std::vector<RES_TENSOR_VMD_MIXING> VMD_MIXING;
  std::vector<std::size_t>           active_vmd;
  bool                               active_ready         = false;
  regge::FFParam                     ff_transfer;
  regge::FFParam                     ff_prod;
};

// One ordered runtime resonance production channel
struct RES_PRODUCTION {
  std::vector<MDecayBranch>               tree;
  HELMatrix                              hel;
  std::optional<spin::PoleLS> pole;
  std::complex<double>                    g = 0.0;  // Direct polarization coupling
};

// Covariant resonance production channels keyed by exchange pair
struct RES_TENSOR_MODEL : public RES_MODEL_FORM {
  double                          phi = 0.0;
  std::vector<RES_TENSOR_CHANNEL> channels;
};

// Select the reduced relativistic resonance denominator
enum class BreitWigner { FixedWidth, KinematicWidth, RunningWidth };

// Select the active resonance production dynamics without overloaded flags
// Resonance parameters filled from steering cards by MResonance.cc
class PARAM_RES {
 public:
  PARAM_RES() = default;

  // Particle class
  MParticle p;
  bool      C_from_card = false;

  // Per-model resonance-card configuration
  RES_PRODUCTION_MODEL    XP;
  RES_PRODUCTION_MODEL    GP;
  RES_MP_PRODUCTION_MODEL MP;
  RES_TENSOR_MODEL        TP;

  // Model directory used for parity and measured decay widths
  std::string modelparam;

  // Prepared production-density index, absent when metrics are disabled
  std::optional<std::size_t> spin_index;

  // MP resonance polarization steering
  std::string                                    spin_basis = "none";
  std::vector<std::complex<double>>              a_Jz;

  // Compute true when the MP model block selects coherent a_Jz steering
  bool UsesCoherentSpinBasis() const { return spin_basis == "a_Jz"; }

  // Compute true when the MP model block selects spin-density rho steering
  bool UsesDensitySpinBasis() const { return spin_basis == "rho"; }

  // Compute true when production dynamics determines the resonance spin state
  bool UsesUnrestrictedSpinBasis() const { return spin_basis == "none"; }

  // Compute true when the MP model block selects a supported polarization mode
  bool UsesMPPolarizationMode() const {
    return UsesUnrestrictedSpinBasis() || UsesCoherentSpinBasis() || UsesDensitySpinBasis();
  }

  // Reduced relativistic resonance denominator
  BreitWigner BW = BreitWigner::FixedWidth;

  // Active resonance production dynamics selected by the process string
  ReggeProductionModel production_model = ReggeProductionModel::None;

  // Spin-density matrix (2J+1) (constant)
  MMatrix<std::complex<double>> rho;

  // --------------------------------------------------------------------
  // Production and decay amplitude MHelicity information (dynamic)
  MMatrix<std::complex<double>> prod_f;
  MMatrix<std::complex<double>> decay_f;

  // --------------------------------------------------------------------
  // Resonance decay helicity amplitude information
  HELMatrix hel_decay;
  // --------------------------------------------------------------------

  // Ordered runtime production channels
  std::vector<RES_PRODUCTION> production;
};

namespace resonance {

// Read a resonance card with the selected model pole and GP for non Regge processes
PARAM_RES Read(const std::string &resonance_file, MRandom &rng, ReggeProductionModel model,
               const std::string &modelparam = gra::MODELPARAM);

// Build coherent spin amplitudes from diagonal Jz probabilities and phases
std::vector<std::complex<double>> CoherentAJzFromDiagonalWeights(const PARAM_RES           &res,
                                                                 const std::vector<double> &probabilities);

// Compute the diphoton branching ratio from the active decay table
double GammaGammaBranchingRatio(int resonance_pdg, const std::string &modelparam = gra::MODELPARAM);

// Compute the electronic partial width from particle and decay tables
double ElectronicPartialWidth(const MParticle &resonance, const std::string &modelparam = gra::MODELPARAM);

// Compute the on shell diphoton partial width from particle and decay tables
double GammaGammaPartialWidth(const MParticle &resonance, const std::string &modelparam = gra::MODELPARAM);

// Compute the width normalized reduced gamma gamma resonance coupling
double GammaGammaResonanceCoupling(const MParticle &resonance, const std::string &modelparam = gra::MODELPARAM);

// Evaluate the reduced line shape selected by one resonance configuration
std::complex<double> LineShape(double mass2, const PARAM_RES &resonance, double running_profile = 1.0);

// Evaluate a reduced fixed width relativistic Breit Wigner line shape
std::complex<double> FixedWidthLineShape(double mass2, double pole_mass, double width);

// Evaluate a reduced kinematic width relativistic Breit Wigner line shape
std::complex<double> KinematicWidthLineShape(double mass2, double pole_mass, double width);

// Evaluate a reduced LS running width relativistic Breit Wigner line shape
std::complex<double> RunningWidthLineShape(double mass2, double pole_mass, double width, double running_profile);

// Evaluate the MadWidth propagator estimate for approximate decay widths
// This reference helper does not replace a physical spin propagator
// A_J(s) = N_J(s)/(s-M0^2+i M0 Gamma), 0 <= J <= 2
std::complex<double> JacksonLineShape(double mass2, double pole_mass, double width, double spin);

// Compute the normalized relativistic Breit Wigner density
double BreitWignerDensity(double mass2, double pole_mass, double width);

// Compute the positive square root of the Breit Wigner density
double BreitWignerAmplitude(double mass2, double pole_mass, double width);

// Compute the normalized complex reduced Breit Wigner amplitude
std::complex<double> ComplexBreitWignerAmplitude(double mass2, double pole_mass, double width);

}  // namespace resonance

}  // namespace gra

#endif
