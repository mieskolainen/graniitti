// Convert weighted pion pair CSV data to HepMC3
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Program/data2hepmc3.h"

#include <filesystem>
#include <iostream>

// Convert one data file and report success only after closing the output
int main(int argc, char* argv[]) {
  if (argc != 2) {
    std::cerr << "Usage: data2hepmc3 data.csv" << std::endl;
    return EXIT_FAILURE;
  }
  try {
    const std::filesystem::path input(argv[1]);
    const auto output = std::filesystem::path("output") / (input.filename().string() + ".hepmc3");
    const int  events = gra::program::ConvertData(input.string(), output.string(), 5.5e-6);
    std::cout << "Converted " << events << " events to " << output << std::endl;
    return EXIT_SUCCESS;
  } catch (const std::exception& error) {
    std::cerr << "data2hepmc3: " << error.what() << std::endl;
    return EXIT_FAILURE;
  }
}
