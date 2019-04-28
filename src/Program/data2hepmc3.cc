// GRANIITTI - Monte Carlo event generator for high energy diffraction
// https://github.com/mieskolainen/graniitti
//
// <HepMC33 to LHE (.xml) format converter>
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <fstream>
#include <iomanip>
#include <iostream>

// Own
#include "Graniitti/M4Vec.h"
#include "Graniitti/MAux.h"
#include "Graniitti/MH1.h"
#include "Graniitti/MPDG.h"

// HepMC3
#include "HepMC3/FourVector.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/LHEFAttributes.h"
#include "HepMC3/Print.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/Relatives.h"
#include "HepMC3/Selector.h"
#include "HepMC3/WriterAscii.h"

std::shared_ptr<HepMC3::GenRunInfo>  runinfo      = nullptr;
std::shared_ptr<HepMC3::WriterAscii> outputHepMC3 = nullptr;

int main(int argc, char *argv[]) {
  using namespace gra;

  if (argc != 2) {
    std::cerr << "<Data to HepMC3 container>" << std::endl;
    std::cerr << "Example: ./data2hepmc3 datafile.csv" << std::endl;
    return EXIT_FAILURE;
  }

  std::string DATAFILE(argv[1]);
  std::string outputfile = "./output/" + DATAFILE + ".hepmc3";

  // Force cross-section (barn)
  const double xsforced = 5.5e-6;

  // Read events
  FILE *fp;

  if ((fp = fopen(DATAFILE.c_str(), "r+")) == NULL) {
    printf("No inputfile %s found \n", DATAFILE.c_str());
    return false;
  }

  // --------------------------------------------------------------
  // Generator info
  runinfo = std::make_shared<HepMC3::GenRunInfo>();

  struct HepMC3::GenRunInfo::ToolInfo generator = {
      std::string("data2hepmc3"), std::to_string(gra::aux::GetVersion()).substr(0, 5),
      std::string("data")};
  runinfo->tools().push_back(generator);

  // struct HepMC3::GenRunInfo::ToolInfo config = {
  //	FULL_INPUT_STR, "1.0", std::string("Steering card")};
  // runinfo->tools().push_back(config);

  outputHepMC3 = std::make_shared<HepMC3::WriterAscii>(outputfile, runinfo);
  // --------------------------------------------------------------
  double dummy     = 0;
  double P1_gen[3] = {0.0};
  double P2_gen[3] = {0.0};
  //  double P1_rec[3]  = {0.0};
  //  double P2_rec[3]  = {0.0};
  //  int    PDG1       = 0;
  //  int    PDG2       = 0;
  //  int    REC        = 0;
  double PID_weight = 0;

  const int DATATYPE = 9;  // # number of variables

  // Diagnostics
  MH1<double> hM(40, 0, 2.5, "System Mass (GeV)");

  // Loop over events
  unsigned int events = 0;
  while (true) {
    int ret = -1;

    ret = fscanf(fp, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf", &P1_gen[0], &P1_gen[1], &P1_gen[2],
                 &dummy, &P2_gen[0], &P2_gen[1], &P2_gen[2], &dummy, &PID_weight);

    /*
    ret = fscanf(fp, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%d,%d,%d",
                                             &P1_gen[0], &P1_gen[1], &P1_gen[2],
                                             &P2_gen[0], &P2_gen[1], &P2_gen[2],
                                             &P1_rec[0], &P1_rec[1], &P1_rec[2],
                                             &P2_rec[0], &P2_rec[1], &P2_rec[2],
                                             &PDG1, &PDG2, &REC);
    */

    if (ret == DATATYPE) {  // number of values per line
                            // Fine
    } else if (ret == EOF) {
      break;
    } else {
      printf("Error in the file structure of %s!\n", DATAFILE.c_str());
      return false;
    }

    const double m1 = PDG::mpi;  // Pion mass
    const double m2 = PDG::mpi;

    // Create 4-vectors
    gra::M4Vec pf1;
    pf1.SetPxPyPzM(P1_gen[0], P1_gen[1], P1_gen[2], m1);
    gra::M4Vec pf2;
    pf2.SetPxPyPzM(P2_gen[0], P2_gen[1], P2_gen[2], m2);

    gra::M4Vec system = pf1 + pf2;

    hM.Fill(system.M());

    // Event weight (NOT USED CURRENTLY)
    //    double weight = 1.0;

    // ******************************************************************
    // Create HepMC3 event
    HepMC3::GenEvent evt(HepMC3::Units::GEV, HepMC3::Units::MM);

    // Dummy beams
    gra::M4Vec beam1(0, 0, 1000, 1000);
    gra::M4Vec beam2(0, 0, -1000, 1000);

    HepMC3::GenParticlePtr gen_beam1 = std::make_shared<HepMC3::GenParticle>(
        gra::aux::M4Vec2HepMC3(beam1), PDG::PDG_p, PDG::PDG_BEAM);
    HepMC3::GenParticlePtr gen_beam2 = std::make_shared<HepMC3::GenParticle>(
        gra::aux::M4Vec2HepMC3(beam2), PDG::PDG_p, PDG::PDG_BEAM);

    HepMC3::GenParticlePtr gen_system = std::make_shared<HepMC3::GenParticle>(
        gra::aux::M4Vec2HepMC3(system), PDG::PDG_system, PDG::PDG_INTERMEDIATE);
    HepMC3::GenParticlePtr gen_p1f = std::make_shared<HepMC3::GenParticle>(
        gra::aux::M4Vec2HepMC3(pf1), PDG::PDG_pip, PDG::PDG_STABLE);
    HepMC3::GenParticlePtr gen_p2f = std::make_shared<HepMC3::GenParticle>(
        gra::aux::M4Vec2HepMC3(pf2), PDG::PDG_pim, PDG::PDG_STABLE);
    HepMC3::GenVertexPtr v1 = std::make_shared<HepMC3::GenVertex>();
    HepMC3::GenVertexPtr v2 = std::make_shared<HepMC3::GenVertex>();

    v1->add_particle_in(gen_beam1);
    v1->add_particle_in(gen_beam2);
    v1->add_particle_out(gen_system);

    v2->add_particle_in(gen_system);
    v2->add_particle_out(gen_p1f);
    v2->add_particle_out(gen_p2f);

    // Finally add all vertices to the event
    evt.add_vertex(v1);
    evt.add_vertex(v2);

    // Save cross section information (HepMC3 format wants it event by event)
    std::shared_ptr<HepMC3::GenCrossSection> xsobj = std::make_shared<HepMC3::GenCrossSection>();
    evt.add_attribute("GenCrossSection", xsobj);

    // Now add the value in picobarns [HepMC3 convention]
    xsobj->set_cross_section(xsforced * 1E12, 0);  // external fixed one

    // Save event weight (unweighted events with weight 1)
    // const double HepMC3_weight = PID_weight;
    const double HepMC3_weight = 1.0;
    evt.weights().push_back(HepMC3_weight);  // add more weights with .push_back()

    // Write event out
    outputHepMC3->write_event(evt);
    // ******************************************************************

    ++events;

    if (events % 100000 == 0) { std::cout << "Event " << events << " processed" << std::endl; }
  }
  fclose(fp);

  hM.Print();

  std::cout << "Conversion done for " << events << " events" << std::endl;


  // Done
  std::cout << "[data2hepmc3:: done]" << std::endl;

  return EXIT_SUCCESS;
}
