// Weighted central mass distributions for optimal transport
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef PROGRAM_OT_H
#define PROGRAM_OT_H

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

#include "Graniitti/Analysis/MHepMCReader.h"
#include "Graniitti/Particle/MPDG.h"
#include "Graniitti/Tech/MH1.h"
#include "HepMC3/GenParticle.h"

namespace gra::program {

// Fill the central mass in GeV with a finite nonnegative nominal event weight
inline void ReadMass(const std::string& path, MH1<double>& histogram) {
  MHepMCReader     reader(path);
  HepMC3::GenEvent event;
  double           accepted = 0.0;
  while (reader.Read(event)) {
    const double weight = event.weights().empty() ? 1.0 : event.weights()[0];
    if (!std::isfinite(weight) || weight < 0.0) {
      throw std::invalid_argument("ot: nominal weights must be finite and nonnegative");
    }
    HepMC3::FourVector system;
    for (const auto& particle : event.particles()) {
      if (particle->status() == PDG::PDG_STABLE && particle->pid() > -1000 && particle->pid() < 1000) {
        system = system + particle->momentum();
      }
    }
    const double mass = system.m();
    if (!std::isfinite(mass)) { throw std::invalid_argument("ot: non-finite central mass"); }
    int bin = 0;
    histogram.GetBinIdx(mass, bin);
    histogram.Fill(mass, weight);
    if (bin >= 0) { accepted += weight; }
  }
  const auto density = histogram.GetProbDensity();
  if (!(accepted > 0.0) || !std::isfinite(accepted) ||
      !std::all_of(density.begin(), density.end(), [](double p) { return std::isfinite(p) && p >= 0.0; }) ||
      std::none_of(density.begin(), density.end(), [](double p) { return p > 0.0; })) {
    throw std::invalid_argument("ot: no positive weight in the selected mass range");
  }
}

}  // namespace gra::program
#endif
