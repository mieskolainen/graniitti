// HepMC event samples for analysis regression tests
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef ANALYSIS_TEST_SUPPORT_HH
#define ANALYSIS_TEST_SUPPORT_HH

#include <cmath>
#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

#include "Graniitti/Particle/MPDG.h"
#include "Graniitti/Tech/MAux.h"
#include "HepMC3/GenCrossSection.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/WriterAscii.h"

namespace analysis_test {

using gra::aux::indices;

// Keep generated event samples under the project temporary directory
inline std::string Path(const std::string &name) {
  const auto directory = std::filesystem::path(gra::aux::GetBasePath(2)) / "tmp" / "test_analysis";
  std::filesystem::create_directories(directory);
  return (directory / name).string();
}

// Build a central decay with a distinct particle-antiparticle pair or one photon
inline HepMC3::GenEvent Event(int pdg = 211, bool pair = true) {
  HepMC3::GenEvent event;
  gra::MPDG table;
  table.ReadParticleData();
  const double mass = table.FindByPDG(pdg).mass;
  const HepMC3::FourVector first(0.4, 0.1, 0.2, std::sqrt(0.21 + mass * mass));
  const HepMC3::FourVector second(-0.3, -0.1, -0.2, std::sqrt(0.14 + mass * mass));
  auto vertex = std::make_shared<HepMC3::GenVertex>();
  vertex->add_particle_in(std::make_shared<HepMC3::GenParticle>(
      pair ? first + second : first, gra::PDG::PDG_system, 2));
  vertex->add_particle_out(std::make_shared<HepMC3::GenParticle>(first, pdg, 1));
  if (pair) { vertex->add_particle_out(std::make_shared<HepMC3::GenParticle>(second, -pdg, 1)); }
  event.add_vertex(vertex);
  event.weights() = {1.0};
  auto cross_section = std::make_shared<HepMC3::GenCrossSection>();
  cross_section->set_cross_section(1.0e12, 0.0);
  event.set_cross_section(cross_section);
  return event;
}

// Write numbered events with the requested HepMC momentum units
inline std::string Write(const std::string &name, std::vector<HepMC3::GenEvent> events,
                         HepMC3::Units::MomentumUnit unit = HepMC3::Units::GEV) {
  const auto path = Path(name);
  HepMC3::WriterAscii writer(path);
  for (const auto &i : indices(events)) {
    events[i].set_event_number(static_cast<int>(i));
    events[i].set_units(unit, HepMC3::Units::MM);
    writer.write_event(events[i]);
  }
  writer.close();
  return path;
}

// Append malformed or truncated input after a valid event
inline std::string Broken(const std::string &name, const std::string &suffix) {
  const auto good = Write(name + "_good", {Event()});
  const auto path = Path(name);
  std::ifstream input(good);
  std::ofstream output(path);
  std::string line;
  while (std::getline(input, line)) {
    if (line.find("END_EVENT_LISTING") != std::string::npos) { break; }
    output << line << '\n';
  }
  output << suffix;
  return path;
}

}  // namespace analysis_test

#endif
