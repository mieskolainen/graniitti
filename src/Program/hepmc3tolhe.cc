// Convert HepMC3 events to Les Houches Event format
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <cstdlib>
#include <iostream>
#include <string>

// Own
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Program/hepmc3tolhe.h"

// For LHE event format see:
// [REFERENCE: http://home.thep.lu.se/~torbjorn/talks/fnal04lha.pdf]
int main(int argc, char *argv[]) {
  gra::aux::PrintArgv(argc, argv);

  if (argc != 2) {
    std::cout << std::endl;
    std::cout << "[HepMC3 to LHEF converter]" << std::endl << std::endl;
    std::cout << "Example: ./hepmc3tolhe filename.hepmc3" << std::endl;

    gra::aux::CheckUpdate();
    return EXIT_FAILURE;
  }

  try {
    std::string inputfile(argv[1]);
    std::string outputfile = inputfile + ".lhe";

    const gra::MLHEConversionStats stats = gra::ConvertHepMC3ToLHE(inputfile, outputfile, true);

    if (stats.events > 0) {
      printf("HepMC3: input  (%0.1f MB, %0.5f MB/event) %s \n", stats.input_size_mb,
             stats.input_size_mb / stats.events, inputfile.c_str());
      printf("LHEF:   output (%0.1f MB, %0.5f MB/event) %s \n", stats.output_size_mb,
             stats.output_size_mb / stats.events, outputfile.c_str());
      printf("Total %d events converted from HepMC3 to LHE \n", stats.events);
    }

    std::cout << "[hepmc3tolhe: done]" << std::endl;
    gra::aux::CheckUpdate();
  } catch (const std::exception& error) {
    std::cerr << "hepmc3tolhe: " << error.what() << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
