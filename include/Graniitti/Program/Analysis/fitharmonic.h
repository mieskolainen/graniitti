// Checked harmonic measurement output
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef PROGRAM_FITHARMONIC_H
#define PROGRAM_FITHARMONIC_H

#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

#include "Graniitti/Program/MInput.h"
#include "Graniitti/Tech/MAux.h"
#include "json.hpp"

namespace gra::program {

// Write and close the full measurement before reporting success
inline void WriteOutput(const std::string& path, const nlohmann::json& output) {
  gra::aux::CreateDirectory(std::filesystem::path(path).parent_path().string());
  std::ofstream file(path);
  if (!file.is_open()) { throw std::invalid_argument("fitharmonic: cannot open output " + path); }
  file << output.dump(2) << '\n';
  file.close();
  if (file.fail()) { throw std::runtime_error("fitharmonic: failed output " + path); }
}

}  // namespace gra::program
#endif
