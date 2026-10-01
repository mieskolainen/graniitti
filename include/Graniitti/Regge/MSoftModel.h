// Immutable soft Regge model parameters
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MSOFTMODEL_H
#define MSOFTMODEL_H

// C++
#include <complex>
#include <cstddef>
#include <map>
#include <memory>
#include <span>
#include <string>
#include <vector>

// Own
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Regge/MReggeSig.h"

namespace gra {

// Select one cached physical final-state basis
enum class GoodWalkerFinalBasis { Proton, Excited, Complete };

// Immutable Good Walker mixing geometry shared by hard and soft amplitudes
class GoodWalkerSpace {
public:
  // Construct and validate one arbitrary-channel Good Walker space
  GoodWalkerSpace(std::size_t channel_count, std::vector<double> mixing_angles,
                  std::vector<double> resolved_coefficients);

  // Compute the number of Good Walker eigenstates
  std::size_t ChannelCount() const noexcept { return channel_count_; }

  // Compute the dimension of the two-proton eigenstate pair space
  std::size_t PairDimension() const noexcept { return pair_dimension_; }

  // Compute the row-major pair-space index for eigenstates i and k
  std::size_t PairIndex(std::size_t i, std::size_t k) const;

  // Compute the lexicographically ordered Good Walker mixing angles
  const std::vector<double> &MixingAngles() const noexcept {
    return mixing_angles_;
  }

  // Compute the coefficients defining the normalized resolved direction
  const std::vector<double> &ResolvedCoefficients() const noexcept {
    return resolved_coefficients_;
  }

  // Compute the physical-state to eigenstate orthogonal mixing matrix
  const MMatrix<double> &MixingMatrix() const noexcept { return mixing_; }

  // Compute the physical proton vector in the eigenstate basis
  const std::vector<double> &ProtonVector() const noexcept { return proton_; }

  // Compute the normalized resolved excitation vector in the eigenstate basis
  const std::vector<double> &ResolvedVector() const noexcept {
    return resolved_;
  }

  // Compute the physical proton projector
  const MMatrix<double> &ProtonProjector() const noexcept {
    return proton_projector_;
  }

  // Compute the complete physical excited-state projector
  const MMatrix<double> &ExcitedProjector() const noexcept {
    return excited_projector_;
  }

  // Compute the resolved-direction projector
  const MMatrix<double> &ResolvedProjector() const noexcept {
    return resolved_projector_;
  }

  // Compute the inclusive complement of the resolved-direction projector
  const MMatrix<double> &InclusiveProjector() const noexcept {
    return inclusive_projector_;
  }

  // Compute one cached orthonormal physical final-state basis
  const MMatrix<double> &FinalBasis(GoodWalkerFinalBasis basis) const;

  // Project one pair-space source onto two cached physical final-state bases
  std::vector<std::complex<double>>
  ProjectPair(std::span<const std::complex<double>> source,
              GoodWalkerFinalBasis upper, GoodWalkerFinalBasis lower) const;

  // Project one pair-space source onto two explicit orthonormal bases
  std::vector<std::complex<double>>
  ProjectPair(std::span<const std::complex<double>> source,
              const MMatrix<double> &upper_basis,
              const MMatrix<double> &lower_basis) const;

private:
  std::size_t channel_count_ = 0;
  std::size_t pair_dimension_ = 0;
  std::vector<double> mixing_angles_;
  std::vector<double> resolved_coefficients_;
  MMatrix<double> mixing_;
  std::vector<double> proton_;
  std::vector<double> resolved_;
  MMatrix<double> proton_projector_;
  MMatrix<double> excited_projector_;
  MMatrix<double> resolved_projector_;
  MMatrix<double> inclusive_projector_;
  MMatrix<double> proton_basis_;
  MMatrix<double> excited_basis_;
  MMatrix<double> complete_basis_;
};

// Strong index identifying one exchange inside a SoftModel snapshot
class SoftExchangeId {
public:
  // Construct one exchange identifier from its immutable model index
  explicit constexpr SoftExchangeId(const std::size_t value) : value_(value) {}

  // Access the immutable model index
  constexpr std::size_t Value() const noexcept { return value_; }

  // Compare two exchange identifiers
  friend constexpr bool operator==(const SoftExchangeId first,
                                   const SoftExchangeId second) noexcept {
    return first.value_ == second.value_;
  }

  // Order two exchange identifiers by their immutable model index
  friend constexpr bool operator<(const SoftExchangeId first,
                                  const SoftExchangeId second) noexcept {
    return first.value_ < second.value_;
  }

private:
  std::size_t value_ = 0;
};

// Physical exchange family used for mapping central Regge trajectories
enum class SoftExchangeRole { Pomeron, Reggeon, Odderon };

// Analytic trajectory form selected by one soft exchange
enum class SoftTrajectoryMode { Linear, PionLoop };

// Unitarization map selected for the soft eikonal
enum class SoftUnitarization { Exponential, QExponential };

// Transition form factor prescription for off diagonal eigenstates
enum class SoftTransitionFormFactor { Diagonal, Arithmetic, Geometric };

// Supported eigenstate form factor families
enum class SoftFormFactor {
  Exponential,
  ExponentialPower,
  DualPower,
  MixedExponential,
  OddThreeGluonNode,
  ThreeGluon,
  GeneralizedKernel
};

// Compute the steering name of one soft exchange role
std::string SoftExchangeRoleName(SoftExchangeRole role);

// Compute the steering name of one trajectory mode
std::string SoftTrajectoryModeName(SoftTrajectoryMode mode);

// Immutable complete parameters of one soft Regge exchange
class SoftExchange {
public:
  // Construct one validated soft exchange value
  SoftExchange(std::string name, SoftExchangeRole role,
               SoftTrajectoryMode trajectory_mode, bool enabled,
               int crossing_parity, regge::Signature signature, double alpha0,
               double alpha_prime, MMatrix<double> coupling, int residue_sign,
               EtaMode eta_mode, std::string form_factor_name,
               SoftFormFactor form_factor,
               std::vector<std::vector<double>> form_factor_parameters,
               SoftTransitionFormFactor transition_form_factor,
               MMatrix<double> helicity_flip_coupling,
               MMatrix<double> helicity_flip_slope);

  // Compute the unique exchange name inside the model snapshot
  const std::string &Name() const noexcept { return name_; }

  // Compute the physical exchange family
  SoftExchangeRole Role() const noexcept { return role_; }

  // Compute the analytic trajectory mode
  SoftTrajectoryMode TrajectoryMode() const noexcept {
    return trajectory_mode_;
  }

  // Check whether the exchange participates in amplitudes
  bool Enabled() const noexcept { return enabled_; }

  // Compute the beam crossing parity
  int CrossingParity() const noexcept { return crossing_parity_; }

  // Compute the Regge signature
  regge::Signature Signature() const noexcept { return signature_; }

  // Compute the trajectory intercept
  double Alpha0() const noexcept { return alpha0_; }

  // Compute the linear trajectory slope
  double AlphaPrime() const noexcept { return alpha_prime_; }

  // Compute the complete Good Walker coupling matrix
  const MMatrix<double> &CouplingMatrix() const noexcept { return coupling_; }

  // Compute one Good Walker transition coupling
  double Coupling(std::size_t i, std::size_t j) const;

  // Compute the real exchange coupling sign
  int ResidueSign() const noexcept { return residue_sign_; }

  // Compute the Regge eta factor prescription
  EtaMode Eta() const noexcept { return eta_mode_; }

  // Compute the named form factor bank
  const std::string &FormFactorName() const noexcept {
    return form_factor_name_;
  }

  // Compute the form factor family
  SoftFormFactor FormFactorType() const noexcept { return form_factor_; }

  // Compute the eigenstate form factor parameters
  const std::vector<std::vector<double>> &
  FormFactorParameters() const noexcept {
    return form_factor_parameters_;
  }

  // Compute the off diagonal transition form factor prescription
  SoftTransitionFormFactor TransitionFormFactor() const noexcept {
    return transition_form_factor_;
  }

  // Compute one Pauli flip coupling
  double HelicityFlipCoupling(std::size_t i, std::size_t j) const;

  // Compute one Pauli flip slope
  double HelicityFlipSlope(std::size_t i, std::size_t j) const;

private:
  std::string name_;
  SoftExchangeRole role_ = SoftExchangeRole::Reggeon;
  SoftTrajectoryMode trajectory_mode_ = SoftTrajectoryMode::Linear;
  bool enabled_ = true;
  int crossing_parity_ = 1;
  regge::Signature signature_ = regge::Signature::Positive;
  double alpha0_ = 0.0;
  double alpha_prime_ = 0.0;
  MMatrix<double> coupling_;
  int residue_sign_ = 1;
  EtaMode eta_mode_ = EtaMode::Rotating;
  std::string form_factor_name_;
  SoftFormFactor form_factor_ = SoftFormFactor::DualPower;
  std::vector<std::vector<double>> form_factor_parameters_;
  SoftTransitionFormFactor transition_form_factor_ =
      SoftTransitionFormFactor::Geometric;
  MMatrix<double> helicity_flip_coupling_;
  MMatrix<double> helicity_flip_slope_;
};

// Immutable matrix eikonal controls from one soft model snapshot
class SoftEikonalSettings {
public:
  // Construct one validated set of matrix eikonal controls
  SoftEikonalSettings(SoftUnitarization unitarization, double q,
                      bool helicity_enabled, double helicity_mass_scale,
                      std::vector<SoftExchangeId> screening_exchanges);

  // Compute the selected eikonal unitarization map
  SoftUnitarization Unitarization() const noexcept { return unitarization_; }

  // Compute the q exponential deformation parameter
  double Q() const noexcept { return q_; }

  // Check whether proton helicity transitions are enabled
  bool HelicityEnabled() const noexcept { return helicity_enabled_; }

  // Compute the proton helicity mass scale
  double HelicityMassScale() const noexcept { return helicity_mass_scale_; }

  // Compute the exchanges included in event screening
  const std::vector<SoftExchangeId> &ScreeningExchanges() const noexcept {
    return screening_exchanges_;
  }

private:
  SoftUnitarization unitarization_ = SoftUnitarization::Exponential;
  double q_ = 1.0;
  bool helicity_enabled_ = false;
  double helicity_mass_scale_ = 0.0;
  std::vector<SoftExchangeId> screening_exchanges_;
};

// Immutable forward excitation profile from PARAM_SOFT.FORWARD_EXCITATION
class SoftForwardExcitationProfile {
public:
  // Construct one validated inelastic proton transition profile
  SoftForwardExcitationProfile(double s0, double a);

  // Compute the inelastic profile reference scale squared
  double S0() const noexcept { return s0_; }

  // Compute the inelastic momentum transfer scale squared
  double A() const noexcept { return a_; }

private:
  double s0_ = 0.0;
  double a_ = 0.0;
};

// Immutable complete soft exchange and Good Walker physics snapshot
class SoftModel {
public:
  // Parse one complete model snapshot from GENERAL JSON text
  static std::shared_ptr<const SoftModel>
  LoadFromJson(const std::string &source_file, const std::string &json_text);

  // Parse one model snapshot with GENERAL and NUMERICS JSON text
  static std::shared_ptr<const SoftModel>
  LoadFromJson(const std::string &source_file, const std::string &json_text,
               const std::string &numerics_source_file,
               const std::string &numerics_json);

  // Compute the canonical source file path
  const std::string &SourceFile() const noexcept { return source_file_; }

  // Compute the exact immutable GENERAL JSON text used by this snapshot
  const std::string &SourceJson() const noexcept { return source_json_; }

  // Compute the NUMERICS source path captured with this model snapshot
  const std::string &NumericsSourceFile() const noexcept {
    return numerics_source_file_;
  }

  // Compute the exact immutable NUMERICS JSON text captured on the master
  const std::string &NumericsSourceJson() const noexcept {
    return numerics_source_json_;
  }

  // Compute true when the master captured a complete NUMERICS snapshot
  bool HasNumericsSnapshot() const noexcept {
    return !numerics_source_file_.empty() && !numerics_source_json_.empty();
  }

  // Compute the selected PARAM_SOFT model name
  const std::string &ActiveModel() const noexcept { return active_model_; }

  // Compute the deterministic physics snapshot fingerprint
  const std::string &Fingerprint() const noexcept { return fingerprint_; }

  // Compute the number of configured soft exchanges
  std::size_t ExchangeCount() const noexcept { return exchanges_.size(); }

  // Resolve one named exchange inside this snapshot
  SoftExchangeId ExchangeId(const std::string &name) const;

  // Compute one exchange inside this snapshot
  const SoftExchange &Exchange(SoftExchangeId exchange) const;

  // Compute every configured exchange in immutable model order
  const std::vector<SoftExchange> &Exchanges() const noexcept {
    return exchanges_;
  }

  // Compute the Pomeron selected for precontracted forward excitation
  SoftExchangeId ForwardExcitationExchange() const noexcept {
    return forward_excitation_exchange_;
  }

  // Access the immutable arbitrary-channel Good Walker geometry
  const GoodWalkerSpace &GoodWalker() const noexcept { return good_walker_; }

  // Access the immutable matrix eikonal controls
  const SoftEikonalSettings &Eikonal() const noexcept { return eikonal_; }

  // Compute the Pomeron pion loop form factor scale squared
  double PionLoopScale2() const noexcept { return pion_loop_scale2_; }

  // Compute the raw triple Pomeron coupling ratio
  double TriplePomeronRatio() const noexcept { return triple_pomeron_ratio_; }

  // Compute the triple Pomeron eta factor prescription
  EtaMode TriplePomeronEtaMode() const noexcept {
    return triple_pomeron_eta_mode_;
  }

  // Access the immutable full-range forward excitation profile
  const SoftForwardExcitationProfile &
  ForwardExcitationProfile() const noexcept {
    return forward_excitation_;
  }

  // Evaluate one mapped linear or pion loop trajectory
  double Alpha(SoftExchangeId exchange, double t) const;

  // Evaluate one diagonal eigenstate form factor
  double FormFactor(SoftExchangeId exchange, double t, std::size_t i) const;

  // Evaluate one symmetric eigenstate transition form factor
  double TransitionFormFactor(SoftExchangeId exchange, double t, std::size_t i,
                              std::size_t j) const;

  // Compute one complete exchange vertex matrix at fixed momentum transfer
  MMatrix<double> ResidueMatrix(SoftExchangeId exchange, double t) const;

  // Compute one complete Pauli flip vertex matrix at fixed momentum transfer
  MMatrix<double> HelicityFlipResidueMatrix(SoftExchangeId exchange,
                                            double t) const;

  // Evaluate only the diagonal Good Walker vertices into caller storage
  void DiagonalResidues(SoftExchangeId exchange, double t,
                        std::span<double> output) const;

  // Compute one exchange coupling projected onto the physical proton
  double PhysicalCoupling(SoftExchangeId exchange) const;

  // Compute one exchange vertex projected onto the physical proton
  double PhysicalResidue(SoftExchangeId exchange, double t) const;

  // Compute the physical exchange form factor normalized at zero transfer
  double NormalizedPhysicalResidue(SoftExchangeId exchange, double t) const;

  // Compute the effective mapped triple Pomeron coupling
  double EffectiveTriplePomeronCoupling(SoftExchangeId exchange) const;

  // Compute the real symmetric principal triple Pomeron matrix square root
  MMatrix<double> TriplePomeronCouplingRoot(SoftExchangeId exchange) const;

  // Evaluate the full-range mapped inelastic proton transition profile
  double ForwardExcitationFactor(SoftExchangeId exchange, double t,
                                 double mass2) const;

private:
  // Construct one already parsed and validated immutable model snapshot
  SoftModel(std::string source_file, std::string source_json,
            std::string numerics_source_file, std::string numerics_source_json,
            std::string active_model, std::string fingerprint,
            std::vector<SoftExchange> exchanges,
            std::map<std::string, SoftExchangeId> exchange_ids,
            GoodWalkerSpace good_walker, SoftEikonalSettings eikonal,
            SoftExchangeId forward_excitation_exchange, double pion_loop_scale2,
            double triple_pomeron_ratio, EtaMode triple_pomeron_eta_mode,
            SoftForwardExcitationProfile forward_excitation);

  std::string source_file_;
  std::string source_json_;
  std::string numerics_source_file_;
  std::string numerics_source_json_;
  std::string active_model_;
  std::string fingerprint_;
  std::vector<SoftExchange> exchanges_;
  std::map<std::string, SoftExchangeId> exchange_ids_;
  GoodWalkerSpace good_walker_;
  SoftEikonalSettings eikonal_;
  SoftExchangeId forward_excitation_exchange_;
  double pion_loop_scale2_ = 0.0;
  double triple_pomeron_ratio_ = 0.0;
  EtaMode triple_pomeron_eta_mode_ = EtaMode::Rotating;
  SoftForwardExcitationProfile forward_excitation_;
};

using SoftModelPtr = std::shared_ptr<const SoftModel>;

} // namespace gra

#endif
