// Convert weighted pion pair CSV data to HepMC3
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef PROGRAM_DATA2HEPMC3_H
#define PROGRAM_DATA2HEPMC3_H

#include <array>
#include <filesystem>
#include <fstream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>

#include "Graniitti/Kinematics/M4Vec.h"
#include "Graniitti/Particle/MPDG.h"
#include "Graniitti/Program/MInput.h"
#include "Graniitti/Tech/MAux.h"
#include "HepMC3/GenCrossSection.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/WriterAscii.h"

namespace gra::program {

// Build the pion pair and retain the supplied PID correction weight
inline HepMC3::GenEvent PionData(const std::array<double, 9>& row, int number, double xs) {
  using gra::aux::indices;
  if (!std::isfinite(xs) || xs < 0.0 || !std::isfinite(xs * 1e12)) {
    throw std::invalid_argument("data2hepmc3: invalid cross section in pb");
  }
  for (const auto& i : indices(row)) {
    if (!std::isfinite(row[i])) { throw std::invalid_argument("data2hepmc3: non-finite CSV value"); }
  }
  M4Vec first, second;
  first.SetPxPyPzM(row[0], row[1], row[2], PDG::mpi);
  second.SetPxPyPzM(row[4], row[5], row[6], PDG::mpi);
  if (!std::isfinite(first.E() + second.E())) {
    throw std::invalid_argument("data2hepmc3: pion energy overflow");
  }
  HepMC3::GenEvent event(HepMC3::Units::GEV, HepMC3::Units::MM);
  event.set_event_number(number);
  // Retain the container's dummy beams for the existing analysis convention
  auto collision = std::make_shared<HepMC3::GenVertex>();
  collision->add_particle_in(
      std::make_shared<HepMC3::GenParticle>(HepMC3::FourVector(0, 0, 1000, 1000), PDG::PDG_p, PDG::PDG_BEAM));
  collision->add_particle_in(
      std::make_shared<HepMC3::GenParticle>(HepMC3::FourVector(0, 0, -1000, 1000), PDG::PDG_p, PDG::PDG_BEAM));
  auto system =
      std::make_shared<HepMC3::GenParticle>(aux::M4Vec2HepMC3(first + second), PDG::PDG_system, PDG::PDG_INTERMEDIATE);
  collision->add_particle_out(system);
  auto decay = std::make_shared<HepMC3::GenVertex>();
  decay->add_particle_in(system);
  decay->add_particle_out(
      std::make_shared<HepMC3::GenParticle>(aux::M4Vec2HepMC3(first), PDG::PDG_pip, PDG::PDG_STABLE));
  decay->add_particle_out(
      std::make_shared<HepMC3::GenParticle>(aux::M4Vec2HepMC3(second), PDG::PDG_pim, PDG::PDG_STABLE));
  event.add_vertex(collision);
  event.add_vertex(decay);
  auto cross_section = std::make_shared<HepMC3::GenCrossSection>();
  event.set_cross_section(cross_section);
  cross_section->set_cross_section(xs * 1e12, 0.0);
  event.weights() = {row[8]};
  return event;
}

// Convert complete CSV rows and check the input and the closed output stream
inline int ConvertData(const std::string& input, const std::string& output, double xs) {
  if (!std::isfinite(xs) || xs < 0.0 || !std::isfinite(xs * 1e12)) {
    throw std::invalid_argument("data2hepmc3: invalid cross section in pb");
  }
  std::ifstream source(input);
  if (!source.is_open()) { throw std::invalid_argument("data2hepmc3: cannot open " + input); }
  DistinctFiles(input, output);
  gra::aux::CreateDirectory(std::filesystem::path(output).parent_path().string());
  std::ofstream target(output);
  if (!target.is_open()) { throw std::invalid_argument("data2hepmc3: cannot open " + output); }
  int events = 0;
  {
    auto run = std::make_shared<HepMC3::GenRunInfo>();
    run->tools().push_back({"data2hepmc3", std::to_string(aux::GetVersion()), "data"});
    run->set_weight_names({"PID"});
    HepMC3::WriterAscii writer(target, run);
    std::string         line;
    while (std::getline(source, line)) {
      if (line.find_first_not_of(" \t\r") == std::string::npos) { continue; }
      std::istringstream    record(line);
      std::array<double, 9> row{};
      for (const auto& i : gra::aux::indices(row)) {
        if (!(record >> row[i]) || !std::isfinite(row[i])) {
          throw std::invalid_argument("data2hepmc3: invalid CSV number");
        }
        if (i + 1 < row.size()) {
          char comma = 0;
          if (!(record >> comma) || comma != ',') { throw std::invalid_argument("data2hepmc3: expected comma"); }
        }
      }
      record >> std::ws;
      if (!record.eof()) { throw std::invalid_argument("data2hepmc3: extra CSV fields"); }
      if (events == std::numeric_limits<int>::max()) { throw std::invalid_argument("data2hepmc3: too many events"); }
      auto event = PionData(row, events, xs);
      writer.write_event(event);
      if (writer.failed() || target.fail()) { throw std::runtime_error("data2hepmc3: failed output " + output); }
      ++events;
    }
    if (source.bad() || !source.eof()) { throw std::runtime_error("data2hepmc3: failed input " + input); }
    writer.close();
  }
  // WriterAscii closes the supplied file stream when writing its footer
  if (target.is_open()) { target.close(); }
  if (target.fail()) { throw std::runtime_error("data2hepmc3: failed output " + output); }
  return events;
}

}  // namespace gra::program
#endif
