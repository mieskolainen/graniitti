// Generate Pythia8 events from a steering card and write HepMC3
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

// Pythia8
#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC3.h"

using namespace Pythia8;

namespace {

// Compute a positive integer parsed from a command-line argument
int ReadPositiveInt(const char *value, const std::string &name) {
  std::size_t parsed = 0;
  const int   result = std::stoi(value, &parsed);
  if (parsed != std::string(value).size() || result <= 0) {
    throw std::runtime_error(name + " must be a positive integer");
  }
  return result;
}

// Compute a non-negative integer parsed from a command-line argument
int ReadNonNegativeInt(const char *value, const std::string &name) {
  std::size_t parsed = 0;
  const int   result = std::stoi(value, &parsed);
  if (parsed != std::string(value).size() || result < 0) {
    throw std::runtime_error(name + " must be a non-negative integer");
  }
  return result;
}

// Compute the Pythia XML data path from the environment or the local default
std::string PythiaDataPath() {
  const char *env_path = std::getenv("PYTHIA8DATA");
  if (env_path != nullptr && std::string(env_path).size() > 0) {
    return std::string(env_path);
  }
  return "../pythia8317/share/Pythia8/xmldoc";
}

// Print command-line usage
void PrintUsage(const char *argv0) {
  std::cerr << "Usage: " << argv0
            << " steering.cmnd output.hepmc3 nevents [seed]\n";
}

// Configure one Pythia instance from a steering file and optional seed
bool ConfigurePythia(Pythia &pythia, const std::string &steering, int seed) {
  if (!pythia.readFile(steering)) { return false; }
  pythia.readString("Main:timesAllowErrors = 100");
  pythia.readString("Next:numberShowEvent = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowInfo = 0");

  if (seed > 0) {
    pythia.readString("Random:setSeed = on");
    pythia.readString("Random:seed = " + std::to_string(seed));
  }
  return pythia.init();
}

// Run Pythia and stream accepted events to a HepMC3 file
int RunPythia(const std::string &steering, const std::string &hepmc_file,
              int nevents, int seed) {
  Pythia pythia(PythiaDataPath(), false);
  if (!ConfigurePythia(pythia, steering, seed)) { return 0; }

  Pythia8ToHepMC to_hepmc(hepmc_file);

  int accepted = 0;
  int aborted  = 0;
  while (accepted < nevents) {
    if (!pythia.next()) {
      if (pythia.info.atEndOfFile()) { break; }
      if (++aborted > pythia.mode("Main:timesAllowErrors")) { break; }
      continue;
    }

    to_hepmc.setWeightNames(pythia.info.weightNameVector());
    if (!to_hepmc.writeNextEvent(pythia)) { return 0; }
    ++accepted;
  }

  to_hepmc.output().close();
  pythia.stat();
  return to_hepmc.output().failed() ? 0 : accepted;
}

}  // namespace

// Main program for steering-card driven Pythia HepMC3 generation
int main(int argc, char *argv[]) {
  if (argc < 4 || argc > 5) {
    PrintUsage(argv[0]);
    return EXIT_FAILURE;
  }

  const std::string steering   = argv[1];
  const std::string hepmc_file = argv[2];
  const int         nevents    = ReadPositiveInt(argv[3], "nevents");
  const int         seed       = (argc == 5) ? ReadNonNegativeInt(argv[4], "seed") : 0;

  const int accepted = RunPythia(steering, hepmc_file, nevents, seed);
  std::cout << "PYTHIA_HEPMC3 accepted " << accepted << " events into "
            << hepmc_file << std::endl;

  return accepted == nevents ? EXIT_SUCCESS : EXIT_FAILURE;
}
