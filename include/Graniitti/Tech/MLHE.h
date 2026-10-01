// GRANIITTI Les Houches Event writer helpers
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MLHE_H
#define MLHE_H

// C++
#include <fstream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

// HepMC3
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/LHEF.h"

namespace gra {

// Les Houches event-weight interpretation
enum class MLHEWeightType { Unit, Weighted, SignedWeighted };

// Fixed run-level weight conversion into the selected LHA convention
struct MLHERunConfig {
  MLHEWeightType weight_type = MLHEWeightType::Unit;
  double         weight_scale = 1.0;
};

// One particle row after HepMC to LHE mapping
struct MLHERow {
  HepMC3::ConstGenParticlePtr particle = nullptr;
  long                        pid      = 0;
  int                         status   = 0;
  std::pair<int, int>         mothers  = {0, 0};
  std::pair<int, int>         colors   = {0, 0};
  std::vector<double>         momentum;
};

// Summary of one HepMC3 to LHE conversion run
struct MLHEConversionStats {
  int    events         = 0;
  double input_size_mb  = 0.0;
  double output_size_mb = 0.0;
};

// Build the generic Les Houches particle-row view of one HepMC event
std::vector<MLHERow> BuildLHERows(const HepMC3::GenEvent &ev);

// Build the optional GRANIITTI diffraction metadata LHE comment for one event
std::string BuildLHEDiffractionComment(const HepMC3::GenEvent &ev);

// Convert one HepMC3 ASCII file into one Les Houches Event file
MLHEConversionStats ConvertHepMC3ToLHE(const std::string &inputfile,
                                       const std::string &outputfile,
                                       bool progress);

class MLHEWriter {
 public:
  // Open one Les Houches Event output file
  explicit MLHEWriter(const std::string &outputfile,
                      const MLHERunConfig &config = {});

  // Flush and close the Les Houches Event output file
  ~MLHEWriter();

  // Copy, assignment and move disabled
  MLHEWriter(const MLHEWriter &)            = delete;
  MLHEWriter &operator=(const MLHEWriter &) = delete;
  MLHEWriter(MLHEWriter &&)                 = delete;
  MLHEWriter &operator=(MLHEWriter &&)      = delete;

  // Write one HepMC event as one Les Houches event
  void WriteEvent(const HepMC3::GenEvent &ev);

  // Flush the complete LHE output and report stream failures
  void Close();

 private:
  // Initialize the Les Houches run block from the first event
  void InitializeRunBlock(const HepMC3::GenEvent &ev);

  // Require later events to retain the initialized run identity
  void ValidateRunBlock(const HepMC3::GenEvent &ev) const;

  std::ofstream output;
  std::unique_ptr<LHEF::Writer> writer;
  MLHERunConfig                 config;
  bool run_block_initialized = false;
};

}  // namespace gra

#endif
