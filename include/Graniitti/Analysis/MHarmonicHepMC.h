// HepMC input adapter for detector-aware angular measurements
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef ANALYSIS_MHARMONICHEPMC_H
#define ANALYSIS_MHARMONICHEPMC_H

// C++
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

// Own
#include "Graniitti/Analysis/MHarmonic.h"

namespace gra {
namespace harmonic {

// Identify the angular rest-frame convention independently of measurement mode
enum class AngularFrame { CM, HX, CS, AH, PG, GJ };

// Compute a stable printable angular-frame name
std::string ToString(AngularFrame frame);

// Parse an angular frame without accepting implicit aliases
AngularFrame ParseAngularFrame(const std::string &value);

// Hold one inclusive finite fiducial interval
struct FiducialRange {
  double min = 0.0;
  double max = 0.0;

  // Test one finite value against the inclusive interval
  bool Contains(double value) const;

  // Validate the finite non-empty interval
  void Validate(const std::string &name) const;
};

// Define particle-level fiducial cuts applied equally at truth and reco
struct FiducialCuts {
  FiducialRange pion_eta;
  FiducialRange pion_pt;
  FiducialRange proton_xi;
  FiducialRange proton_abs_t;

  // Validate central cuts and tagged-mode forward cuts
  void Validate(MeasurementMode mode) const;
};

// Define reconstructed exclusivity cuts without changing truth fiducial space
struct ExclusiveSelection {
  double pt_balance_max = 0.0;
  double mass_match_relative_max = 0.0;
  double rapidity_match_max = 0.0;

  // Validate tagged-mode data-side closure cuts
  void Validate(MeasurementMode mode) const;
};

// Configure the standalone HepMC measurement adapter
struct HepMCReadConfig {
  MeasurementMode mode = MeasurementMode::Central;
  AngularFrame frame = AngularFrame::CS;
  FiducialCuts cuts;
  ExclusiveSelection selection;
  double sqrt_s = 0.0;
  std::uint64_t max_records = 0;
  std::uint64_t seed = 1;
  std::size_t weight_index = 0;

  // Validate beams, cuts and frame compatibility
  void Validate() const;
};

// Read truth, paired detector simulation and measured HepMC streams
class HepMCAnalysisReader {
public:
  // Construct an immutable event projector
  explicit HepMCAnalysisReader(const HepMCReadConfig &config);

  // Read generated events and synthesize their detector observations
  std::vector<ResponseEvent>
  ReadModelResponse(const std::string &truth_path,
                    const DetectorResponseModel &model) const;

  // Read generated events with a sparse event-number-matched reco stream
  std::vector<ResponseEvent>
  ReadPairedResponse(const std::string &truth_path,
                     const std::string &reco_path) const;

  // Read selected detector-level events from data or closure MC
  std::vector<DataEvent> ReadData(const std::string &path) const;

private:
  HepMCReadConfig config_;
};

} // namespace harmonic
} // namespace gra

#endif
