// Convert HepMC3 ASCII events to HepMC2 ASCII for Delphes pile-up tools
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

// HepMC3
#include "HepMC3/GenEvent.h"
#include "Graniitti/Analysis/MHepMCReader.h"
#include "HepMC3/WriterAsciiHepMC2.h"

namespace {

struct Arguments {
  std::string input_file;
  std::string output_file;
};

// Print command-line usage
void PrintUsage(const char *argv0) {
  std::cerr << "Usage: " << argv0 << " input.hepmc3 output.hepmc2\n";
}

// Parse and validate command-line arguments
Arguments ParseArguments(int argc, char *argv[]) {
  if (argc != 3) {
    PrintUsage(argv[0]);
    throw std::runtime_error("expected exactly two file arguments");
  }

  Arguments arguments{argv[1], argv[2]};
  if (std::filesystem::exists(arguments.output_file) &&
      std::filesystem::equivalent(arguments.input_file, arguments.output_file)) {
    throw std::runtime_error("input and output must be distinct files");
  }
  return arguments;
}

// Convert all HepMC3 events in the input file to HepMC2 ASCII
int ConvertFile(const Arguments &arguments) {
  gra::MHepMCReader        input(arguments.input_file);
  HepMC3::WriterAsciiHepMC2 output(arguments.output_file);
  HepMC3::GenEvent          event;

  int converted = 0;
  while (input.Read(event)) {
    output.write_event(event);
    if (output.failed()) { throw std::runtime_error("HepMC2 output write failed"); }
    ++converted;
  }

  output.close();
  if (output.failed()) { throw std::runtime_error("HepMC2 output close failed"); }
  return converted;
}

}  // namespace

// Main program for HepMC3 to HepMC2 conversion
int main(int argc, char *argv[]) {
  try {
    const Arguments arguments = ParseArguments(argc, argv);
    const int       converted = ConvertFile(arguments);
    std::cout << "HEPMC3_TO_HEPMC2 converted " << converted << " events into "
              << arguments.output_file << std::endl;
    return converted > 0 ? EXIT_SUCCESS : EXIT_FAILURE;
  } catch (const std::exception &error) {
    std::cerr << "hepmc3_to_hepmc2: " << error.what() << std::endl;
    return EXIT_FAILURE;
  }
}
