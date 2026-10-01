// YFS type QED initial and final state photon radiation
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MRADIATIVE_H
#define MRADIATIVE_H

// C++
#include <array>
#include <cstddef>
#include <exception>
#include <string>
#include <vector>

// HepMC3
#include "HepMC3/GenEvent.h"

// Libraries
#include "json.hpp"

// Own
#include "Graniitti/Kinematics/M4Vec.h"
#include "Graniitti/Particle/MParticle.h"
#include "Graniitti/Sampling/MRandom.h"

namespace gra {
namespace radiative {

// Available QED radiation algorithms
enum class Mode { Off, YFS };

// Numerical controls for lepton beam radiation
struct ISRParam {
  double x_min = 1.0e-6;
  // Add finite collinear terms through O(alpha^2) to the soft structure function
  bool hard_collinear = true;
};

// Numerical controls for radiation from final lepton pairs
struct FSRParam {
  double energy_min = 1.0e-3;
  double energy_max_fraction = 0.5;
  std::size_t max_photons = 32;
  std::size_t max_trials = 10000;
};

// Complete immutable radiative configuration copied to each worker
struct Config {
  Mode isr = Mode::Off;
  Mode fsr = Mode::Off;
  ISRParam isr_param;
  FSRParam fsr_param;
};

// Event-local initial state radiation and hard-beam recoil
struct ISRState {
  bool active = false;
  std::array<M4Vec, 2> nominal;
  std::array<M4Vec, 2> hard;
  std::vector<M4Vec> photons;
  double weight = 1.0;
};

// One sampled final state lepton-pair recoil
struct FSRPair {
  std::array<int, 2> pdg{};
  std::array<M4Vec, 2> born;
  std::array<M4Vec, 2> lepton;
  std::vector<M4Vec> photons;
};

// Event-local final state radiation and fiducial decay tree
struct FSRState {
  bool prepared = false;
  std::vector<MDecayBranch> tree;
  std::vector<FSRPair> pairs;
};

// Parse and validate the gencard modes and tune numerical controls
Config ReadConfig(const nlohmann::json &modes, const nlohmann::json &numerics);

// Compute the canonical steering name of one radiative mode
std::string ModeName(Mode mode);

// Compute the exponent of the soft electron structure function
double ISRExponent(double scale2, double mass);

// Compute the integrated massive final state dipole factor
double DipoleIntegral(double beta);

// Compute the massive final state dipole angular kernel
double DipoleAngular(double beta, double cosine);

// Map a lepton pair against explicit photons with exact four-momentum closure
bool MapLeptonPair(const M4Vec &lepton1, const M4Vec &lepton2,
                   const std::vector<M4Vec> &photons, M4Vec &mapped1,
                   M4Vec &mapped2);

// Worker-local YFS type radiation service
class MYFS {
 public:
  // Configure both radiative sectors
  void Configure(const Config &config);

  // Validate that enabled radiation sectors have physical charged leptons
  void ValidateApplicability(const MParticle &beam1, const MParticle &beam2,
                             const std::vector<MDecayBranch> &tree) const;

  // Store the nominal incoming beam momenta
  void SetBeams(const M4Vec &beam1, const M4Vec &beam2);

  // Compute one sampling coordinate for each radiating lepton beam
  unsigned int ISRDim(const MParticle &beam1, const MParticle &beam2) const;

  // Restore the nominal incoming state before a new sampled point
  void ResetISR(M4Vec &beam1, M4Vec &beam2, double &s, double &sqrt_s) noexcept;

  // Generate the ISR convolution point and exact hard-beam recoil
  double GenerateISR(const std::vector<double> &random, std::size_t offset,
                     const MParticle &beam1, const MParticle &beam2,
                     M4Vec &momentum1, M4Vec &momentum2, double &s,
                     double &sqrt_s);

  // Sample FSR after the Born amplitude without changing its momenta
  void PrepareFSR(const std::vector<MDecayBranch> &tree, MRandom &random);

  // Compute the post-FSR tree used for central fiducial decisions
  const std::vector<MDecayBranch> &
  FiducialTree(const std::vector<MDecayBranch> &born) const noexcept;

  // Clear the event-local final state radiation sample
  void ResetFSR() noexcept;

  // Add accepted ISR and FSR photons to one complete event record
  bool Apply(HepMC3::GenEvent &event, MRandom &random, std::exception_ptr *failure = nullptr) noexcept;

  // Compute the event-local ISR state
  const ISRState &GetISR() const noexcept { return isr_state; }

  // Compute the event-local FSR state
  const FSRState &GetFSR() const noexcept { return fsr_state; }

  // Access the configured modes
  const Config &GetConfig() const noexcept { return config; }

private:
  // Attach the sampled initial state to the HepMC graph
  bool AttachISR(HepMC3::GenEvent &event);

  // Radiate all supported final state lepton-pair vertices
  bool ApplyFSR(HepMC3::GenEvent &event, MRandom &random);

  Config config;
  ISRState isr_state;
  FSRState fsr_state;
};

} // namespace radiative
} // namespace gra

#endif
