// Generate inclusive gamma*/Z -> mu+ mu- events with Pythia8 and write HepMC3
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

// HepMC3
#include "HepMC3/Attribute.h"
#include "HepMC3/GenCrossSection.h"

// Pythia8
#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC3.h"

#ifndef PYTHIA_XML_DIR
#define PYTHIA_XML_DIR ""
#endif

using namespace Pythia8;

namespace {

struct HardMuonCuts {
  bool enabled = false;
  double pt_min = 20.0;
  double abs_eta_max = 2.5;
  double mass_min = 80.0;
  double mass_max = 100.0;
};

struct HardMuonKinematics {
  bool valid = false;
  double mu_minus_pt = 0.0;
  double mu_minus_eta = 0.0;
  double mu_plus_pt = 0.0;
  double mu_plus_eta = 0.0;
  double mumu_mass = 0.0;
};

struct CommandLineOptions {
  std::string steering;
  std::string hepmc_file;
  int nevents = 0;
  int seed = 0;
  HardMuonCuts hard_muon_cuts;
};

struct RunSummary {
  int written = 0;
  int aborted = 0;
  long tried = 0;
  long selected = 0;
  long accepted = 0;
  double sigma_pb = 0.0;
  double sigma_err_pb = 0.0;
};

// Compute a positive integer parsed from a command-line argument
int ReadPositiveInt(const char *value, const std::string &name) {
  std::size_t parsed = 0;
  const int result = std::stoi(value, &parsed);
  if (parsed != std::string(value).size() || result <= 0) {
    throw std::runtime_error(name + " must be a positive integer");
  }
  return result;
}

// Compute a non-negative integer parsed from a command-line argument
int ReadNonNegativeInt(const char *value, const std::string &name) {
  std::size_t parsed = 0;
  const int result = std::stoi(value, &parsed);
  if (parsed != std::string(value).size() || result < 0 || result > 900000000) {
    throw std::runtime_error(name + " must be an integer between 0 and 900000000");
  }
  return result;
}

// Compute an on/off command-line value as a boolean
bool ReadOnOff(const std::string &value, const std::string &name) {
  if (value == "on") { return true; }
  if (value == "off") { return false; }
  throw std::runtime_error(name + " must be either 'on' or 'off'");
}

// Print command-line usage
void PrintUsage(std::ostream &output, const char *argv0) {
  output << "Usage: " << argv0
         << " steering.cmnd output.hepmc3 nevents [seed]"
         << " [--hard-muon-cuts on|off]\n"
         << "\n"
         << "Hard-process filter when enabled:\n"
         << "  pT(mu+/-) > 20 GeV, |eta(mu+/-)| < 2.5, "
         << "80 <= m(mu+mu-) <= 100 GeV\n";
}

// Parse the native Pythia driver command line
CommandLineOptions ParseCommandLine(int argc, char *argv[]) {
  if (argc < 4) {
    throw std::runtime_error("missing required command-line arguments");
  }

  CommandLineOptions options;
  options.steering = argv[1];
  options.hepmc_file = argv[2];
  options.nevents = ReadPositiveInt(argv[3], "nevents");

  bool seed_seen = false;
  bool hard_muon_cuts_seen = false;
  for (int index = 4; index < argc; ++index) {
    const std::string argument = argv[index];
    if (argument == "--hard-muon-cuts") {
      if (hard_muon_cuts_seen) {
        throw std::runtime_error("--hard-muon-cuts was specified more than once");
      }
      if (++index >= argc) {
        throw std::runtime_error("--hard-muon-cuts requires on or off");
      }
      options.hard_muon_cuts.enabled =
          ReadOnOff(argv[index], "--hard-muon-cuts");
      hard_muon_cuts_seen = true;
      continue;
    }

    const std::string prefix = "--hard-muon-cuts=";
    if (argument.rfind(prefix, 0) == 0) {
      if (hard_muon_cuts_seen) {
        throw std::runtime_error("--hard-muon-cuts was specified more than once");
      }
      options.hard_muon_cuts.enabled =
          ReadOnOff(argument.substr(prefix.size()), "--hard-muon-cuts");
      hard_muon_cuts_seen = true;
      continue;
    }

    if (!argument.empty() && argument.front() == '-') {
      throw std::runtime_error("unknown command-line option: " + argument);
    }
    if (seed_seen) {
      throw std::runtime_error("unexpected positional argument: " + argument);
    }
    options.seed = ReadNonNegativeInt(argv[index], "seed");
    seed_seen = true;
  }
  std::error_code error;
  if (std::filesystem::equivalent(options.steering, options.hepmc_file, error)) {
    throw std::invalid_argument("Steering input and HepMC output refer to the same file");
  }
  return options;
}

// Find one direct hard gamma*/Z decay into an opposite-sign muon pair
bool FindHardProcessMuonPair(const Event &process, int &mu_minus_index,
                             int &mu_plus_index) {
  int candidate_count = 0;
  int candidate_minus = -1;
  int candidate_plus = -1;

  for (int index = 0; index < process.size(); ++index) {
    const Particle &particle = process[index];
    if (particle.id() != 23 || particle.status() != -22) { continue; }

    int local_minus = -1;
    int local_plus = -1;
    bool duplicate_charge = false;
    for (const int daughter_index : particle.daughterList()) {
      if (daughter_index <= 0 || daughter_index >= process.size()) { continue; }
      const int daughter_id = process[daughter_index].id();
      if (daughter_id == 13) {
        duplicate_charge = duplicate_charge || local_minus >= 0;
        local_minus = daughter_index;
      } else if (daughter_id == -13) {
        duplicate_charge = duplicate_charge || local_plus >= 0;
        local_plus = daughter_index;
      }
    }
    if (duplicate_charge || local_minus < 0 || local_plus < 0) { continue; }

    ++candidate_count;
    candidate_minus = local_minus;
    candidate_plus = local_plus;
  }

  if (candidate_count != 1) { return false; }
  mu_minus_index = candidate_minus;
  mu_plus_index = candidate_plus;
  return true;
}

// Check one extracted hard-process muon pair against the fiducial thresholds
bool PassHardMuonKinematics(const HardMuonKinematics &kinematics,
                            const HardMuonCuts &cuts) {
  return kinematics.valid &&
         kinematics.mu_minus_pt > cuts.pt_min &&
         kinematics.mu_plus_pt > cuts.pt_min &&
         std::abs(kinematics.mu_minus_eta) < cuts.abs_eta_max &&
         std::abs(kinematics.mu_plus_eta) < cuts.abs_eta_max &&
         kinematics.mumu_mass >= cuts.mass_min &&
         kinematics.mumu_mass <= cuts.mass_max;
}

// Apply the fiducial cuts to the direct hard-process muon legs
bool PassHardProcessMuonCuts(const Event &process, const HardMuonCuts &cuts,
                             HardMuonKinematics &kinematics) {
  int mu_minus_index = -1;
  int mu_plus_index = -1;
  if (!FindHardProcessMuonPair(process, mu_minus_index, mu_plus_index)) {
    return false;
  }

  const Particle &mu_minus = process[mu_minus_index];
  const Particle &mu_plus = process[mu_plus_index];
  kinematics.mu_minus_pt = mu_minus.pT();
  kinematics.mu_minus_eta = mu_minus.eta();
  kinematics.mu_plus_pt = mu_plus.pT();
  kinematics.mu_plus_eta = mu_plus.eta();
  kinematics.mumu_mass = (mu_minus.p() + mu_plus.p()).mCalc();
  kinematics.valid = std::isfinite(kinematics.mu_minus_pt) &&
                     std::isfinite(kinematics.mu_minus_eta) &&
                     std::isfinite(kinematics.mu_plus_pt) &&
                     std::isfinite(kinematics.mu_plus_eta) &&
                     std::isfinite(kinematics.mumu_mass);
  if (!kinematics.valid) { return false; }

  return PassHardMuonKinematics(kinematics, cuts);
}

class HardProcessMuonCutHook : public UserHooks {
 public:
  // Construct a process-level hook with fixed fiducial thresholds
  explicit HardProcessMuonCutHook(const HardMuonCuts &cuts) : cuts_(cuts) {}

  // Advertise the process-level veto to Pythia
  bool canVetoProcessLevel() override { return true; }

  // Reject hard events whose direct gamma*/Z decay legs fail the cuts
  bool doVetoProcessLevel(Event &process) override {
    ++seen_;
    last_passed_kinematics_ = HardMuonKinematics{};
    HardMuonKinematics kinematics;
    if (PassHardProcessMuonCuts(process, cuts_, kinematics)) {
      last_passed_kinematics_ = kinematics;
      ++passed_;
      return false;
    }
    if (!kinematics.valid) {
      ++malformed_;
      throw std::runtime_error(
          "hard-muon process filter requires exactly one finite direct "
          "gamma*/Z -> mu+ mu- pair");
    }
    ++kinematic_vetoed_;
    return true;
  }

  // Compute the number of hard candidates inspected by the hook
  long seen() const { return seen_; }

  // Compute the number of hard candidates accepted by the hook
  long passed() const { return passed_; }

  // Compute the number rejected by the requested kinematic thresholds
  long kinematicVetoed() const { return kinematic_vetoed_; }

  // Compute the number of malformed hard-process records encountered
  long malformed() const { return malformed_; }

  // Compute the kinematics accepted by the most recent process-level decision
  const HardMuonKinematics &lastPassedKinematics() const {
    return last_passed_kinematics_;
  }

 private:
  const HardMuonCuts cuts_;
  HardMuonKinematics last_passed_kinematics_;
  long seen_ = 0;
  long passed_ = 0;
  long kinematic_vetoed_ = 0;
  long malformed_ = 0;
};

// Configure Pythia from one steering file and an optional command-line seed
bool ConfigurePythia(Pythia &pythia, const std::string &steering, int seed,
                     bool process_veto_enabled) {
  if (!pythia.readFile(steering)) {
    throw std::invalid_argument("Cannot read Pythia settings from " + steering);
  }
  pythia.readString("Main:timesAllowErrors = 100");
  pythia.readString("Next:numberShowEvent = 0");
  pythia.readString("Next:numberShowProcess = 0");
  if (process_veto_enabled) {
    pythia.readString("Check:abortIfVeto = off");
  }

  if (seed > 0) {
    pythia.readString("Random:setSeed = on");
    pythia.readString("Random:seed = " + std::to_string(seed));
  }
  return pythia.init();
}

// Attach generator counters and hard-muon-filter provenance to one event
void AddGeneratorMetadata(HepMC3::GenEvent &event, const Info &info,
                          const HardMuonCuts &cuts,
                          const HardMuonKinematics &kinematics) {
  const auto cross_section = event.cross_section();
  if (!cross_section) {
    throw std::runtime_error("Pythia HepMC3 event has no GenCrossSection");
  }
  cross_section->set_accepted_events(info.nAccepted());
  cross_section->set_attempted_events(info.nTried());

  event.add_attribute("pythia_n_selected",
                      std::make_shared<HepMC3::LongAttribute>(info.nSelected()));
  event.add_attribute(
      "hard_process_muon_cuts_enabled",
      std::make_shared<HepMC3::IntAttribute>(cuts.enabled ? 1 : 0));
  event.add_attribute("hard_process_muon_pt_min_GeV",
                      std::make_shared<HepMC3::DoubleAttribute>(cuts.pt_min));
  event.add_attribute(
      "hard_process_muon_abs_eta_max",
      std::make_shared<HepMC3::DoubleAttribute>(cuts.abs_eta_max));
  event.add_attribute("hard_process_mumu_mass_min_GeV",
                      std::make_shared<HepMC3::DoubleAttribute>(cuts.mass_min));
  event.add_attribute("hard_process_mumu_mass_max_GeV",
                      std::make_shared<HepMC3::DoubleAttribute>(cuts.mass_max));

  if (!kinematics.valid) { return; }
  event.add_attribute(
      "hard_process_mu_minus_pt_GeV",
      std::make_shared<HepMC3::DoubleAttribute>(kinematics.mu_minus_pt));
  event.add_attribute(
      "hard_process_mu_minus_eta",
      std::make_shared<HepMC3::DoubleAttribute>(kinematics.mu_minus_eta));
  event.add_attribute(
      "hard_process_mu_plus_pt_GeV",
      std::make_shared<HepMC3::DoubleAttribute>(kinematics.mu_plus_pt));
  event.add_attribute(
      "hard_process_mu_plus_eta",
      std::make_shared<HepMC3::DoubleAttribute>(kinematics.mu_plus_eta));
  event.add_attribute(
      "hard_process_mumu_mass_GeV",
      std::make_shared<HepMC3::DoubleAttribute>(kinematics.mumu_mass));
}

// Print the process-filter accounting after generation
void PrintFilterSummary(const HardMuonCuts &cuts,
                        const std::shared_ptr<HardProcessMuonCutHook> &hook) {
  std::cout << "PYTHIA_ZMUMU_HARD_MUON_CUTS "
            << (cuts.enabled ? "on" : "off")
            << " pt_min_GeV " << cuts.pt_min
            << " abs_eta_max " << cuts.abs_eta_max
            << " mass_min_GeV " << cuts.mass_min
            << " mass_max_GeV " << cuts.mass_max << '\n';
  if (!hook) { return; }
  std::cout << "PYTHIA_ZMUMU_HARD_MUON_FILTER seen " << hook->seen()
            << " passed " << hook->passed()
            << " kinematic_vetoed " << hook->kinematicVetoed()
            << " malformed " << hook->malformed() << '\n';
}

// Run Pythia and stream accepted events to a HepMC3 file
RunSummary RunPythia(const CommandLineOptions &options) {
  Pythia pythia(PYTHIA_XML_DIR, false);
  std::shared_ptr<HardProcessMuonCutHook> hook;
  if (options.hard_muon_cuts.enabled) {
    hook = std::make_shared<HardProcessMuonCutHook>(options.hard_muon_cuts);
    if (!pythia.setUserHooksPtr(hook)) {
      throw std::runtime_error("failed to install the hard-muon process hook");
    }
  }
  if (!ConfigurePythia(pythia, options.steering, options.seed,
                       options.hard_muon_cuts.enabled)) {
    throw std::runtime_error("Pythia initialization failed");
  }

  Pythia8ToHepMC to_hepmc(options.hepmc_file);
  if (to_hepmc.output().failed()) { throw std::runtime_error("Cannot open HepMC output " + options.hepmc_file); }
  RunSummary summary;
  while (summary.written < options.nevents) {
    if (!pythia.next()) {
      if (pythia.info.atEndOfFile()) { break; }
      if (++summary.aborted > pythia.mode("Main:timesAllowErrors")) { break; }
      continue;
    }

    to_hepmc.setWeightNames(pythia.info.weightNameVector());
    if (!to_hepmc.fillNextEvent(pythia)) {
      throw std::runtime_error("Pythia-to-HepMC3 conversion failed");
    }

    HardMuonKinematics kinematics;
    if (hook) {
      kinematics = hook->lastPassedKinematics();
      if (!PassHardMuonKinematics(kinematics, options.hard_muon_cuts)) {
        throw std::runtime_error(
            "accepted event has no matching hard-muon hook decision");
      }
    } else {
      PassHardProcessMuonCuts(
          pythia.process, options.hard_muon_cuts, kinematics);
    }
    AddGeneratorMetadata(to_hepmc.event(), pythia.info,
                         options.hard_muon_cuts, kinematics);
    to_hepmc.writeEvent();
    if (to_hepmc.output().failed()) {
      throw std::runtime_error("failed while writing the HepMC3 output");
    }
    ++summary.written;
  }

  to_hepmc.output().close();
  pythia.stat();
  summary.tried = pythia.info.nTried();
  summary.selected = pythia.info.nSelected();
  summary.accepted = pythia.info.nAccepted();
  summary.sigma_pb = pythia.info.sigmaGen() * 1e9;
  summary.sigma_err_pb = pythia.info.sigmaErr() * 1e9;
  PrintFilterSummary(options.hard_muon_cuts, hook);
  return summary;
}

}  // namespace

// Main program for inclusive Pythia gamma*/Z -> mu+ mu- HepMC3 generation
int main(int argc, char *argv[]) {
  if (argc == 2 && std::string(argv[1]) == "--help") {
    PrintUsage(std::cout, argv[0]);
    return EXIT_SUCCESS;
  }

  try {
    const CommandLineOptions options = ParseCommandLine(argc, argv);
    const RunSummary summary = RunPythia(options);
    std::cout << "PYTHIA_ZMUMU_HEPMC3 written " << summary.written
              << " events into " << options.hepmc_file
              << " tried " << summary.tried
              << " selected " << summary.selected
              << " accepted " << summary.accepted
              << " aborted " << summary.aborted
              << " sigma_pb " << summary.sigma_pb
              << " sigma_err_pb " << summary.sigma_err_pb << std::endl;
    return summary.written == options.nevents ? EXIT_SUCCESS : EXIT_FAILURE;
  } catch (const std::exception &error) {
    std::cerr << "PYTHIA_ZMUMU_HEPMC3 error: " << error.what() << '\n';
    PrintUsage(std::cerr, argv[0]);
    return EXIT_FAILURE;
  }
}
