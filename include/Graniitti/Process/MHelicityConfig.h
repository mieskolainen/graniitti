// Helicity steering-card loading and vertex construction
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef PROCESS_MHELICITYCONFIG_H
#define PROCESS_MHELICITYCONFIG_H

// C++
#include <array>
#include <complex>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

// Own
#include "Graniitti/MModelTune.h"
#include "Graniitti/Particle/MPDG.h"
#include "Graniitti/Process/MProcessState.h"
#include "Graniitti/Process/MSubProc.h"
#include "Graniitti/Regge/MReggeSig.h"
#include "Graniitti/Spin/MPoleLS.h"
#include "Graniitti/Spin/MHELMatrix.h"
#include "Graniitti/Spin/MHelicity.h"
#include "Graniitti/Spin/MHelicityScatter.h"

// Libraries
#include "json.hpp"

namespace gra {

using json = nlohmann::json;

// Report missing helicity card data without classifying diagnostic text
class MissingHelicityData final : public std::invalid_argument {
 public:
  using std::invalid_argument::invalid_argument;
};

// Store the result of one helicity channel lookup
struct HelicityChannelMatch {
  bool             found = false;
  std::string      channel_id;
  std::vector<int> matched_decay;
  bool             matched_by_leg_exchange       = false;
  bool             matched_by_charge_conjugation = false;
  const json      *channel_block                 = nullptr;
};

// Read decay-only parameters in a two-body channel block
void ParseDecayParameters(HELMatrix &hc, const json &channel_block, const std::string &context, bool production_mode,
                          const std::string &model);

// Read one two-body coupling block in a supported helicity basis
void ParseTwoBodyCouplings(HELMatrix &hc, const json &channel_block, const std::string &context,
                           const std::string &table_name, const MParticle &mother,
                           const std::vector<MParticle> &vertex_legs, const HelicityChannelMatch &match,
                           bool production_mode, gra::spin::VertexContext vertex_context,
                           regge::Signature trajectory_signature, bool analytic_production, int analytic_mmax);

// Build one direct two-body helicity matrix from typed couplings
HELMatrix BuildDirectHelicityCoupling(const MParticle &mother, const std::vector<MParticle> &legs,
                                      const std::vector<std::array<double, 2>> &helicity,
                                      const std::vector<std::complex<double>> &couplings, bool C_symmetry,
                                      bool P_symmetry, bool leg_exchange, const std::string &source_name,
                                      bool verbose_output = true,
                                      const std::string &verbose_label = "",
                                      gra::spin::VertexContext context = gra::spin::VertexContext::Auto);

// Remove validated helicity couplings at or below one magnitude threshold
void PruneHelicityCouplings(HELMatrix &hc, double coupling_min, const std::string &context);

// Compute the PDG ids carried by one particle vector
std::vector<int> ParticlePDGList(const std::vector<MParticle> &legs);

// Format a compact PDG list for diagnostics
std::string FormatPDGList(const std::vector<int> &values);

// Compute a diagnostic label for one JSON helicity channel
std::string ChannelContext(const std::string &table_name, const std::string &pdg_str, const std::string &channel_id);

// Test if a two-body vertex has exchangeable physical legs
bool AllowTwoBodyLegExchange(bool production_mode, gra::spin::VertexContext context, const std::vector<int> &refdecay);

// Build auxiliary JW metadata with the crossed second production leg
std::vector<MParticle> JacobWickVertexLegs(const std::vector<MParticle> &legs, bool production_mode,
                                           gra::spin::VertexContext context, const MPDG &pdg_table);

// Build beam-side photon source metadata without a steering-card row
HELMatrix ForwardPhotonSourceHelicityStructure(const MParticle &particle, const std::vector<MParticle> &legs);

// Build beam-side fixed-spin exchange metadata without a continuum card row
HELMatrix ForwardHadronSourceHelicityStructure(const MParticle &particle, const std::vector<MParticle> &legs);

// Compute the best matching helicity channel from a resonance block
HelicityChannelMatch FindHelicityChannelMatch(const json &resonance_block, const std::vector<int> &refdecay,
                                              const std::string &table_name, const std::string &pdg_str,
                                              bool allow_charge_conjugate_match, bool allow_leg_exchange_match,
                                              bool                          production_pair_schema     = false,
                                              bool                          analytic_production_schema = false,
                                              const std::vector<MParticle> *requested_legs             = nullptr);

// Apply the two-body LS-basis leg-exchange phase
void ApplyTwoBodyLegExchangePhase(HELMatrix &hc, const MParticle &leg1, const MParticle &leg2,
                                  const std::vector<int> &requested_decay, const HelicityChannelMatch &match,
                                  const std::string &table_name);

// Immutable helicity steering cards shared by copied process workers
class MHelicityConfig {
 public:
  // Construct an unloaded helicity configuration
  MHelicityConfig() = default;

  // Load immutable helicity steering cards for one model tune
  explicit MHelicityConfig(MModelTunePtr tune);

  // Load immutable helicity steering cards for one model tune
  void Load(MModelTunePtr tune);

  // Compute whether the helicity steering cards have been loaded
  bool IsLoaded() const noexcept;

  // Construct one helicity vertex from the immutable steering cards
  HELMatrix ProcessHelicityStructure(const MParticle &particle, const std::vector<MParticle> &legs,
                                     const LORENTZSCALAR &lts, const MSubProc &subprocess,
                                     const MModelTunePtr &model_tune, const std::string &decay_mode,
                                     bool production_mode = false, bool strict_mode = false,
                                     const std::string &verbose_label = "", bool verbose_output = true,
                                     gra::spin::VertexContext context = gra::spin::VertexContext::Auto) const;

  // Construct one canonical MP or XP pole operator from its continuum card
  gra::spin::PoleLS ProcessPoleOperatorStructure(ReggeProductionModel model, const MParticle &particle,
                                                            const std::vector<MParticle> &legs,
                                                            const LORENTZSCALAR &lts, const MSubProc &subprocess,
                                                            const MModelTunePtr     &model_tune,
                                                            gra::spin::VertexContext context) const;

  // Construct helicity vertices recursively through one decay branch
  void ProcessHelicityTree(MDecayBranch &branch, const LORENTZSCALAR &lts, const MSubProc &subprocess,
                           const MModelTunePtr &model_tune, const std::string &decay_mode, bool production_mode = false,
                           bool strict_mode = false) const;

 private:
  // Validate every continuum card before process-specific lookup
  void ValidateContinuum() const;

  // Immutable parsed steering JSON objects
  MModelTunePtr                                 tune;
  std::map<std::string, const nlohmann::json *> continuum;
  std::shared_ptr<const nlohmann::json>         decays;

  std::string decays_path;

  // Trajectory signatures keyed by every mapped analytic and fixed PDG alias
  std::map<int, regge::Signature> trajectory_signature;

  // Canonical pole spins keyed by every mapped analytic and fixed PDG alias
  std::map<int, int> trajectory_pole_spinX2;
};

}  // namespace gra

#endif
