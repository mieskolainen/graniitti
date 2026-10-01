// Convert complete HEPEVT records with explicit stream error checks
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef PROGRAM_HEPEVT2HEPMC3_H
#define PROGRAM_HEPEVT2HEPMC3_H

#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

#include "Graniitti/Program/MInput.h"
#include "Graniitti/Tech/MAux.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/ReaderHEPEVT.h"
#include "HepMC3/WriterAscii.h"

namespace gra::program {

// Distinguish EOF before a record from errors inside a record
inline std::size_t ConvertHEPEVT(const std::string& input, const std::string& output) {
  std::ifstream source(input);
  if (!source.is_open()) { throw std::invalid_argument("hepevt2hepmc3: cannot open " + input); }
  DistinctFiles(input, output);
  gra::aux::CreateDirectory(std::filesystem::path(output).parent_path().string());
  std::ofstream target(output);
  if (!target.is_open()) { throw std::invalid_argument("hepevt2hepmc3: cannot open " + output); }
  std::size_t count = 0;
  {
    HepMC3::ReaderHEPEVT reader(source);
    HepMC3::WriterAscii  writer(target);
    while (true) {
      source >> std::ws;
      if (source.bad()) { throw std::runtime_error("hepevt2hepmc3: failed input " + input); }
      if (source.eof()) { break; }
      if (source.fail()) { throw std::invalid_argument("hepevt2hepmc3: invalid input " + input); }
      HepMC3::GenEvent event(HepMC3::Units::GEV, HepMC3::Units::MM);
      if (!reader.read_event(event) || source.fail() || source.bad()) {
        throw std::invalid_argument("hepevt2hepmc3: malformed or truncated event in " + input);
      }
      event.set_run_info(reader.run_info());
      writer.write_event(event);
      if (writer.failed() || target.fail()) { throw std::runtime_error("hepevt2hepmc3: failed output " + output); }
      ++count;
    }
    reader.close();
    writer.close();
  }
  // WriterAscii closes the supplied file stream when writing its footer
  if (target.is_open()) { target.close(); }
  if (target.fail()) { throw std::runtime_error("hepevt2hepmc3: failed output " + output); }
  return count;
}

}  // namespace gra::program
#endif
