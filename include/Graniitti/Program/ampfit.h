// Physical helicity amplitudes on fixed exclusive two-body HepMC3 events
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef PROGRAM_AMPFIT_H
#define PROGRAM_AMPFIT_H

#include <array>
#include <cmath>
#include <complex>
#include <vector>

#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "Graniitti/Kinematics/MEventKinematics.h"
#include "Graniitti/Kinematics/MFactorized.h"
#include "Graniitti/Tech/MException.h"
#include "Graniitti/Tech/MAux.h"

namespace gra::program {

// Evaluate the normal generator amplitude and screening on a common event sample
class AmplitudeProcess final : public MFactorized {
 public:
  // Preserve the fully initialized process and its physical decay operators
  explicit AmplitudeProcess(const MFactorized &process) : MFactorized(process) {
    const auto &lts = state.lts;
    if (lts.process.REGGE_MODEL == ReggeProductionModel::MP) {
      for (const auto &[name, res] : lts.process.RESONANCES) {
        if (res.UsesDensitySpinBasis()) {
          throw std::invalid_argument("ampfit requires coherent MP a_Jz or fusion production, found rho for " + name);
        }
      }
    }
    if (lts.decaytree.size() != 2 || !lts.decaytree[0].legs.empty() || !lts.decaytree[1].legs.empty() ||
        lts.decaytree[0].p.pdg == lts.decaytree[1].p.pdg || lts.excite1 || lts.excite2) {
      throw std::invalid_argument("ampfit requires exclusive production of two distinct stable central particles");
    }
  }
  
  // Import measured four-momenta and use the generators's common amplitude failure boundary
  std::vector<std::complex<double>> Evaluate(const HepMC3::GenEvent &event, MEventWeightState &weight, double &intensity) {
    auto &lts = state.lts;
    std::array<bool, 5> found{};
    for (const auto &particle : event.particles()) {
      if (particle->status() != PDG::PDG_STABLE) { continue; }
      const auto &p = particle->momentum();
      std::size_t index = 0;
      if (particle->pid() == lts.beam1.pdg && p.pz() > 0.0) { index = 1; }
      if (particle->pid() == lts.beam2.pdg && p.pz() < 0.0) { index = 2; }
      // In native generator records the recoil is a direct daughter of a beam
      // Its PDG may also occur in the central pair (for example p pbar)
      const auto parents = particle->parents();
      const bool beam_recoil = parents.size() == 1 && parents.front()->status() == PDG::PDG_BEAM;
      if (!beam_recoil) {
        if (particle->pid() == lts.decaytree[0].p.pdg) { index = 3; }
        if (particle->pid() == lts.decaytree[1].p.pdg) { index = 4; }
      }
      if (index == 0) { continue; }
      if (found[index]) { throw PhaseSpaceFailure("ampfit found repeated final-state particles"); }
      lts.pfinal[index] = aux::HepMC2M4Vec(p);
      found[index] = true;
    }
    for (std::size_t index = 1; index < found.size(); ++index) {
      if (!found[index]) { throw PhaseSpaceFailure("ampfit event is missing a final-state particle"); }
    }
    lts.pfinal[0] = lts.pfinal[3] + lts.pfinal[4];
    lts.decaytree[0].p4 = lts.pfinal[3];
    lts.decaytree[1].p4 = lts.pfinal[4];
    lts.proton_good_walker.reset();
    lts.amplitude.BeginCentral();
    lts.screening.Clear();
    if (!kinematics::SetLorentzScalars(state, 4)) { throw PhaseSpaceFailure("ampfit has invalid event kinematics"); }
    lts.pfinal_orig = lts.pfinal;
    intensity = GetAmp2(GetScreening(), weight);
    if (!weight.Valid()) { return {}; }
    std::vector<std::complex<double>> amplitude(lts.hamp.begin(), lts.hamp.end());
    const double scale = std::sqrt(lts.hamp.metadata.amplitude_normalization);
    for (auto &value : amplitude) { value *= scale; }
    return amplitude;
  }
};

}  // namespace gra::program
#endif
