// Read analysis events with explicit HepMC units and input errors
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef ANALYSIS_MHEPMCREADER_H
#define ANALYSIS_MHEPMCREADER_H

#include <fstream>
#include <stdexcept>
#include <string>

#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/ReaderAscii.h"

namespace gra {

class MHepMCReader {
 public:
  // Open one HepMC stream and retain access to its EOF and error flags
  explicit MHepMCReader(const std::string &path)
      : path_(path), stream_(path), reader_(stream_) {
    if (!stream_.is_open() || stream_.fail()) {
      throw std::invalid_argument("MHepMCReader: cannot open " + path_);
    }
  }

  // Read one event in GeV and mm or reject malformed and truncated records
  bool Read(HepMC3::GenEvent &event) {
    event.clear();
    if (stream_.bad() || (stream_.fail() && !stream_.eof())) {
      throw std::invalid_argument("MHepMCReader: input error in " + path_);
    }
    if (stream_.eof()) { return false; }
    const bool parsed = reader_.read_event(event);
    if (stream_.bad() || !parsed ||
        (stream_.fail() && !stream_.eof())) {
      event.clear();
      throw std::invalid_argument("MHepMCReader: malformed event or input error in " + path_);
    }
    // HepMC can retain default vertices when their records are missing at EOF
    for (const auto &vertex : event.vertices()) {
      if (vertex->particles_in().empty() && vertex->particles_out().empty()) {
        event.clear();
        throw std::invalid_argument("MHepMCReader: uninitialized vertex in " + path_);
      }
    }
    if (stream_.eof() && event.particles().empty()) { return false; }
    // HepMC can retain default particles when the final record is truncated
    for (const auto &particle : event.particles()) {
      if (particle->pid() == 0 && particle->status() == 0) {
        event.clear();
        throw std::invalid_argument("MHepMCReader: uninitialized particle in " + path_);
      }
    }
    event.set_units(HepMC3::Units::GEV, HepMC3::Units::MM);
    return true;
  }

 private:
  std::string path_;
  std::ifstream stream_;
  HepMC3::ReaderAscii reader_;
};

}  // namespace gra

#endif
