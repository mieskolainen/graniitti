// Regge Amplitudes
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MREGGE_H
#define MREGGE_H

// C++
#include <array>
#include <complex>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

// Own
#include "Graniitti/Kinematics/MEventKinematics.h"
#include "Graniitti/Math/MPolarFourier.h"
#include "Graniitti/Process/MAmpMatch.h"
#include "Graniitti/Process/MProcessSetup.h"
#include "Graniitti/Regge/MReggeInclusive.h"
#include "Graniitti/Regge/MReggeMulti.h"
#include "Graniitti/Regge/MReggeNumerics.h"
#include "Graniitti/Regge/MReggeParam.h"
#include "Graniitti/Regge/MReggeProductionModel.h"

namespace gra {

// Select the Regge amplitude family exposed through the process definition
enum class MReggeMode { Generic, Resonance, ContinuumTwoBody, ContinuumTwoFourSixBody, ResonanceContinuumTwoBody, Soft };

// One model-specific photoproduction channel before beam contraction
struct ReggePhotoAmplitude {
  ForwardBeamLeg                    photon_leg   = ForwardBeamLeg::Upper;
  int                               central_pdg  = 0;
  double                            central_mass = 0.0;
  double                            w2           = 0.0;
  std::vector<std::complex<double>> amplitude;
  std::complex<double>              target_forward = 0.0;
  std::complex<double>              target_elastic = 0.0;
  double                            target_slope   = 0.0;
  double                            target_eta     = 0.0;
};

// Compute the beam particle attached to one forward Regge leg
const gra::MParticle &ReggeBeamParticleForLeg(const gra::LORENTZSCALAR &lts, int leg);

// Matrix element dimension: " GeV^" << -(2*external_legs - 8)
class MRegge : public amplitude::ProcessFamily {
 public:
  MRegge(gra::LORENTZSCALAR &lts, MModelTunePtr model_tune, std::shared_ptr<const amplitude::ProcessDefinition> definition);
  ~MRegge() = default;

  // Build one immutable Regge process definition for a selected process mode
  static std::shared_ptr<const amplitude::ProcessDefinition> ProcessDefinitionFor(MReggeMode mode, const std::string &process_name);

  // Initialize immutable Regge parameters before worker process copies
  static void InitializeParameters(const MProcessSetup &setup);

  // Initialize Regge resonance and continuum branching structures
  static void InitializeBranching(MProcessSetup &setup, ReggeProductionModel model, MReggeMode mode);

  // Validate all Regge caches and steering fixed during initialization
  static void ValidateInitializedProcess(MReggeMode mode, const LORENTZSCALAR &lts, const std::string &process_name);

  // Cache the continuum t/u interference signs fixed by model steering
  static void InitializeContinuumInterference(LORENTZSCALAR &lts);

  // Compute sqrt(BornNorm) and store coherent six-particle amplitudes in lts
  std::complex<double> EvalContinuum6(gra::LORENTZSCALAR &lts, ReggeProductionModel model);

  // Compute sqrt(BornNorm) and store coherent four-particle amplitudes in lts
  std::complex<double> EvalContinuum4(gra::LORENTZSCALAR &lts, ReggeProductionModel model);

  // Compute sqrt(BornNorm) and store inclusive elastic or diffractive amplitudes in lts
  std::complex<double> ME2(gra::LORENTZSCALAR &lts, MReggeInclusive mode) const;

  // Evaluate one complete MP, XP or GP central process
  double Amp2(gra::LORENTZSCALAR &lts, ReggeProductionModel model, MReggeMode mode);

  // Resolve selected model-specific photo channels with caller supplied beam states
  std::vector<ReggePhotoAmplitude> PhotoAmplitudes(gra::LORENTZSCALAR &lts, ReggeProductionModel model, const std::array<ForwardLegState, 2> &state, const std::array<double, 2> &subenergy, const std::array<bool, 2> &direction) const;

  // Compute the propagator of the explicitly mapped Pomeron trajectory
  std::complex<double> PomeronKernel(double s, double t) const;

  // Compute the propagator of one central trajectory PDG alias
  std::complex<double> ExchangeKernel(double s, double t, int exchange_pdg) const;

  // Compute the channel-specific photoproduction kernel without vertex couplings
  std::complex<double> PhotoKernel(double s, double t, int central_pdg) const;

  // Compute the selected proton target transition kernel
  std::complex<double> PhotoDissKernel(double s, double t, double mass2, int central_pdg, DissociationType type) const;

  // Access the immutable triple Pomeron coupling matrix square root
  const MMatrix<double> &TriplePomeronCouplingRoot() const { return triple_pomeron_root; }

  // Access the immutable SOFT model used by every Regge beam source
  const SoftModelPtr &SoftModelHandle() const noexcept { return soft_model; }

  // Access the immutable complete tune used by the Regge amplitude
  const MModelTunePtr &ModelTuneHandle() const noexcept { return model_tune; }

 private:
  // Project the completed proton Good Walker source into the Born amplitude
  void ProjectGoodWalkerBorn(gra::LORENTZSCALAR &lts) const;

  // Store one elastic source for an arbitrary proton Good Walker channel count
  void StoreElasticGoodWalker(gra::LORENTZSCALAR &lts, SoftExchangeId exchange, std::complex<double> kernel) const;

  // Store one triple-Regge source without an additional inelastic profile
  void StoreTripleGoodWalker(gra::LORENTZSCALAR &lts, MReggeInclusive mode, std::complex<double> kernel) const;

  // Add one direct multi-Regge ladder term to the proton Good Walker source
  void AddLadderGoodWalker(gra::LORENTZSCALAR &lts, const ForwardLegState &upper_state, const ForwardLegState &lower_state, int upper_exchange, double upper_s, int lower_exchange, double lower_s, std::complex<double> central) const;

  // Accumulate one parallel multi-Regge production topology
  void AccumulateParallelTopology(gra::LORENTZSCALAR &lts, ReggeProductionModel model, const regge::Topology &topology, const std::vector<std::vector<int>> &permutations) const;

  // Evaluate one prepared four- or six-particle continuum ladder
  std::complex<double> EvalContinuumLadder(gra::LORENTZSCALAR &lts, ReggeProductionModel model, std::size_t central_count);

  // Compute one mapped SOFT propagator without beam vertices
  std::complex<double> PropagatorForExchange(double s, double t, SoftExchangeId exchange) const;

  // Compute one prepared mapped SOFT propagator without repeated lookup
  std::complex<double> PreparedPropagator(double s, double t, SoftExchangeId exchange) const;

  MModelTunePtr                        model_tune;
  SoftModelPtr                         soft_model;
  regge::ParamPtr                      param_handle;
  const regge::Param                  &param;
  const regge::ContinuumPlan           continuum_plan;
  const bool                           parallel_required;
  MReggeNumericsPtr                    regge_numerics;
  const bool                           parallel_fourier_required;
  std::unique_ptr<math::MPolarFourier> parallel_fourier;
  MMatrix<double>                      triple_pomeron_root;
};

}  // namespace gra

#endif
