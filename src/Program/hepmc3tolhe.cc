// GRANIITTI - Monte Carlo event generator for high energy diffraction
// https://github.com/mieskolainen/graniitti
//
// <HepMC3 to LHE (.xml) format converter>
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.


// C++
#include <fstream>
#include <iomanip>

// HepMC3
#include "HepMC/FourVector.h"
#include "HepMC/GenEvent.h"
#include "HepMC/GenParticle.h"
#include "HepMC/GenVertex.h"
#include "HepMC/LHEFAttributes.h"
#include "HepMC/Print.h"
#include "HepMC/ReaderAscii.h"
#include "HepMC/Search/FindParticles.h"
#include "HepMC/WriterAscii.h"


// Return filesize for statistics
long GetFileSize(const std::string &filename) {
    std::ifstream file(filename.c_str(),
                       std::ifstream::in | std::ifstream::binary);
    if (!file.is_open())
	return -1;

    file.seekg(0, std::ios::end);
    long fs = file.tellg();
    file.close();

    return fs;
}


// For LHE event format see:
// [REFERENCE: http://home.thep.lu.se/~torbjorn/talks/fnal04lha.pdf]
//
int main(int argc, char *argv[]) {
    if (argc != 2) {
	std::cerr << "<< HepMC3 to LHEF converter>>" << std::endl;
	std::cerr << "Example: ./hepmc3tolhe filename.hepmc3" << std::endl;
	return EXIT_FAILURE;
    }

    std::string inputfile(argv[1]);
    std::string outputfile = inputfile + ".lhe";

    // Input and output
    HepMC::ReaderAscii input(inputfile);
    LHEF::Writer writer(outputfile);

    // TODO: Check what this does actually?
    writer.init();

    int events = 0;

    // HepMC3 event object
    HepMC::GenEvent ev(HepMC::Units::GEV, HepMC::Units::MM);

    // Event loop over all events in HepMC3 file
    while (!input.failed()) {
	// Read event from input file
	input.read_event(ev);
	if (input.failed())
	    break;

	// Create LHE HEPEUP event
	LHEF::HEPEUP hepeup;

	// ** Particle properties **
	std::vector<long> IDUP; // PDG code
	std::vector<int> ISTUP; // status code (-1 in state, 1 = final state, 2
	                        // = intermediate)
	std::vector<std::pair<int, int> > MOTHUP; // position of 1 or 2 mothers
	std::vector<std::pair<int, int> >
	    ICOLUP; // color charges (1 for quarks, 2 for gluons)
	std::vector<std::vector<double> > PUP; // (px,py,pz,E,m)
	std::vector<double> VTIMUP;            // invariant lifetime ctau
	std::vector<double> SPINUP;            // helicity (spin) information

	// All particles
	HepMC::FindParticles search(ev, HepMC::FIND_ALL);

	// Loop over particles
	int NUP = 0;
	FOREACH(const HepMC::GenParticlePtr &p1, search.results()) {
	    // IDUP
	    IDUP.push_back(p1->pid());

	    // ISTUP
	    ISTUP.push_back(p1->status());

	    // Find mother ids
	    HepMC::FindParticles search2(p1, HepMC::FIND_ALL_ANCESTORS);

	    std::vector<int> mother_ids;
	    FOREACH(const HepMC::GenParticlePtr &k, search2.results()) {
		// HepMC::Print::line(k);
		mother_ids.push_back(k->id());
	    }
	    std::pair<int, int> MOTHUP_this(0, 0);
	    const int offset = 0; // convention
	    if (mother_ids.size() == 1) {
		MOTHUP_this.first = mother_ids.at(0) + offset;
	    }
	    if (mother_ids.size() >= 2) {
		MOTHUP_this.first = mother_ids.at(0) + offset;
		MOTHUP_this.second =
		    mother_ids.at(mother_ids.size() - 1) + offset;
	    }

	    // MOTHUP
	    MOTHUP.push_back(MOTHUP_this);

	    // ICOLUP
	    std::pair<int, int> ICOLUP_this(0, 0);
	    ICOLUP.push_back(ICOLUP_this);

	    // PUP
	    HepMC::FourVector pvec = p1->momentum();
	    std::vector<double> PUP_this = {pvec.px(), pvec.py(), pvec.pz(),
	                                    pvec.e(), pvec.m()};
	    PUP.push_back(PUP_this);

	    // VTIMUP
	    VTIMUP.push_back(0);

	    // SPINUP
	    SPINUP.push_back(0);

	    ++NUP;
	}

	hepeup.resize(NUP); // <- Important, otherwise segfault

	// ** Global Event Properties **

	// Number of particle entries
	hepeup.NUP = NUP;
	// Subprocess code (as given in LPRUP)
	hepeup.IDPRUP = 0;
	// Event weight
	hepeup.XWGTUP = 1.0;
	// PDF weights for incoming partons
	hepeup.XPDWUP = std::pair<double, double>(1.0, 1.0);
	// PDF evaluation scale of the event (GeV)
	hepeup.SCALUP = 0;
	// QED coupling of the event
	hepeup.AQEDUP = 0;
	// QCD coupling of the event
	hepeup.AQCDUP = 0;


	// ** Particles **

	// PDG Ids
	hepeup.IDUP = IDUP;
	// Status codes
	hepeup.ISTUP = ISTUP;
	// First and last mother indices
	hepeup.MOTHUP = MOTHUP;
	// Color-flow indices first(second) is (anti)colour
	hepeup.ICOLUP = ICOLUP;
	// Lab frame (px,py,pz,e,m) (GeV)
	hepeup.PUP = PUP;
	// Lifetime
	hepeup.VTIMUP = VTIMUP;
	// Spin (polarization) information
	hepeup.SPINUP = SPINUP;
	
	
	// Write the event out
	writer.hepeup = hepeup;
	writer.hepeup.heprup = &writer.heprup;
	writer.writeEvent();
	++events;

	if (events % 10000 == 0)
	    printf("%d events processed \n", events);
    }

    if (events > 0) {
	double input_size = GetFileSize(inputfile) / 1.0e6;
	double output_size = GetFileSize(outputfile) / 1.0e6;

	printf("HepMC3:: input  (%0.1f MB, %0.5f MB/event) %s \n", input_size,
	       input_size / events, inputfile.c_str());
	printf("LHEF::   output (%0.1f MB, %0.5f MB/event) %s \n", output_size,
	       output_size / events, outputfile.c_str());
	printf("Total %d events converted from HepMC3 to LHE \n", events);
    }

    return EXIT_SUCCESS;
}
