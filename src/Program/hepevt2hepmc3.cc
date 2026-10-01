// Convert HEPEVT to HepMC3
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Program/hepevt2hepmc3.h"

#include <iostream>

// Convert every complete input event and propagate read and write errors
int main(int argc, char* argv[]) {
  if (argc != 2) {
    std::cerr << "Usage: hepevt2hepmc3 filename.hepevt" << std::endl;
    return EXIT_FAILURE;
  }
  try {
    const std::string input(argv[1]);
    const auto        count = gra::program::ConvertHEPEVT(input, input + ".hepmc3");
    std::cout << "Converted " << count << " events" << std::endl;
    return EXIT_SUCCESS;
  } catch (const std::exception& error) {
    std::cerr << "hepevt2hepmc3: " << error.what() << std::endl;
    return EXIT_FAILURE;
  }
}
