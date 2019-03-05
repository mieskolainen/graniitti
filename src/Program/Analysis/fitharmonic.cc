// GRANIITTI - Monte Carlo event generator for high energy diffraction
// https://github.com/mieskolainen/graniitti
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.
//
//
// Spherical Harmonic Moment t_LM Based (M,costheta,phi) Decomposition
// based Efficiency inversion and expansion for 2-body central (exclusive)
// production at the LHC.


// C++
#include <math.h>
#include <algorithm>
#include <chrono>
#include <complex>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <random>
#include <stdexcept>
#include <thread>
#include <vector>

// Own
#include "Graniitti/Analysis/MHarmonic.h"
#include "Graniitti/MAux.h"
#include "Graniitti/MPDG.h"

// Libraries
#include "cxxopts.hpp"
#include "rang.hpp"

// HepMC 3
#include "HepMC/FourVector.h"
#include "HepMC/GenEvent.h"
#include "HepMC/GenParticle.h"
#include "HepMC/GenVertex.h"
#include "HepMC/Print.h"
#include "HepMC/ReaderAscii.h"
#include "HepMC/Search/FindParticles.h"
#include "HepMC/WriterAscii.h"

using gra::aux::indices;

using namespace gra;

// Create Harmonic moment expansion object (needs to be global for TMinuit reasons)
MHarmonic ha;

// Random numbers for fast detector simulation
MRandom rrand;

const double ETACUT = 0.9;
const double PTCUT  = 0.1;

// TMinuit fit wrapper calling the global object
void fitwrapper(int& npar, double* gin, double& f, double* par, int iflag) {
    ha.logLfunc(npar, gin, f, par, iflag);
}

void ReadIn(const std::string inputfile, std::vector<gra::spherical::Omega>& events, const std::string& FRAME, int MAXEVENTS);

// Main function
int main(int argc, char* argv[]) {

    gra::aux::PrintFlashScreen(rang::fg::blue);
    std::cout << rang::style::bold
              << "GRANIITTI - Spherical Harmonics Inverse Expansion"
              << rang::style::reset << std::endl
              << std::endl;
    gra::aux::PrintVersion();
    
    // Save the number of input arguments
    const int NARGC = argc - 1;
    try {

    cxxopts::Options options(argv[0], "");
    options.add_options()
        ("r,ref",     "Reference MC (angular flat MC sample)", cxxopts::value<std::string>() )
        ("i,input",   "Input sample MC or DATA",               cxxopts::value<std::string>() )
        ("l,lmax",    "Maximum angular order (1,2,3,...)",     cxxopts::value<uint>() )
        ("f,frame",   "Lorentz rest frame (HE,CS,GJ,PG,AH)",   cxxopts::value<std::string>() )
        ("m,mass",    "Mass binning [bins,min,max]",           cxxopts::value<std::string>() )
        ("X,maximum", "Maximum number of events",              cxxopts::value<uint>() )
        ("H,help",    "Help")
        ;
    auto r = options.parse(argc, argv);
    
    if (r.count("help") || NARGC == 0) {
        std::cout << options.help({""}) << std::endl;
        std::cout << "Example:" << std::endl;
        std::cout << "  " << argv[0] << " -r SPHERICAL_2pi_FLAT -i SPHERICAL_2pi -l 4 -f HE -m 40,0.4,1.5"
                  << std::endl << std::endl;
        return EXIT_FAILURE;
    }
    
    const uint LMAX          = r["lmax"].as<uint>();
    const std::string FRAME  = r["frame"].as<std::string>();

    const std::string Mstr   = r["mass"].as<std::string>();
    std::vector<double> M    = gra::aux::SplitStr(Mstr,double(0),',');
  
    /*
    const std::string Pstr   = r["pt"].as<std::string>();
    std::vector<double> PT   = gra::aux::SplitStr(Pstr,double(0),',');

    const std::string Ystr   = r["Y"].as<std::string>();
    std::vector<double> Y    = gra::aux::SplitStr(Ystr,double(0),',');
    */

    // ------------------------------------------------------------------
    // INITIALIZE HARMONIC EXPANSION

    MHarmonic::HPARAM hparam;

    // Discretization
    //hparam.M  = {M_BIN,   M_MIN,  M_MAX};
    //hparam.PT = {PT_BIN, PT_MIN, PT_MAX};
    //hparam.Y  = {Y_BIN,   Y_MIN,  Y_MAX};

    // Moments
    hparam.LMAX            = LMAX;
    hparam.REMOVEODD       = false;
    hparam.REMOVENEGATIVEM = true;
    hparam.LAMBDA          = 1e-6;

    ha.Init(hparam);


    // ------------------------------------------------------------------
    // DATA

    const std::string refinput = r["ref"].as<std::string>();
    const std::string datinput = r["input"].as<std::string>();

    // Read in MC
    int MAXEVENTS = 1E9;
    if (r.count("X")) {
      MAXEVENTS = r["X"].as<uint>();
    }
    
    // Read in reference MC
    std::vector<gra::spherical::Omega> MC_events;
    ReadIn(refinput, MC_events, FRAME, MAXEVENTS);
    
    // Read in DATA
    std::vector<gra::spherical::Omega> DATA_events;
    ReadIn(datinput, DATA_events, FRAME, MAXEVENTS);


    // ------------------------------------------------------------------
    // LOOP OVER HYPERBINS

    ha.HyperLoop(fitwrapper, MC_events, DATA_events);

    // Print out results for external analysis
    std::string outputname(argv[2]);
    ha.PrintLoop(outputname);

    // Plot all
    ha.PlotAll();


  } catch (const std::invalid_argument& e) {
      gra::aux::PrintGameOver();
      std::cerr << rang::fg::red
                << "Exception catched: " << rang::fg::reset << e.what()
                << std::endl;
      return EXIT_FAILURE;
  } catch (const std::ios_base::failure& e) {
      gra::aux::PrintGameOver();
      std::cerr << rang::fg::red
                << "Exception catched: std::ios_base::failure: "
                << rang::fg::reset << e.what() << std::endl;
      return EXIT_FAILURE;
  } catch (const cxxopts::OptionException& e) {
      gra::aux::PrintGameOver();
      std::cerr << rang::fg::red
                << "Exception catched: Commandline options: " << rang::fg::reset << e.what() << std::endl;
      return EXIT_FAILURE;
  } catch (...) {
      gra::aux::PrintGameOver();
      std::cerr << rang::fg::red << "Exception catched: Unspecified (...) (Probably input)" << rang::fg::reset
                << std::endl;
      return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}


// Fast toy simulation pt-efficiency parameter (put infinite for perfect pt-efficiency)
const double pt_scale = 4.5;


// Read events in
void ReadIn(const std::string inputfile, std::vector<gra::spherical::Omega>& events, const std::string& FRAME, int MAXEVENTS) {

  const std::string total_input = gra::aux::GetBasePath(2) + "/output/" + inputfile + ".hepmc3";
  printf("ReadIn:: Maximum event count %d \n", MAXEVENTS);
  printf("Reading %s \n", total_input.c_str());
  
  HepMC::ReaderAscii input_file(total_input);
  // Reading failed
  if (input_file.failed()) {
    throw std::invalid_argument("MHarmonic::ReadIn: Failed to open <" + total_input + ">");
  }
  // Event loop
  int events_read = 0;

  // Dummy vector
  M4Vec pvec;
  HepMC::GenEvent evt(HepMC::Units::GEV, HepMC::Units::MM);

  // Allocate memory here for speed
  events.resize((int) 1e7);
  printf("Memory allocated \n");

  // Event loop
  while (!input_file.failed()) {

    if (events_read >= MAXEVENTS) { break; }

    // Read event from input file
    input_file.read_event(evt);

    // FIND SYSTEM PARTICLES
    HepMC::FindParticles search_pip(evt, HepMC::FIND_ALL, HepMC::PDG_ID == PDG::PDG_pip && HepMC::STATUS == PDG::PDG_STABLE);
    HepMC::FindParticles search_pim(evt, HepMC::FIND_ALL, HepMC::PDG_ID == PDG::PDG_pim && HepMC::STATUS == PDG::PDG_STABLE);

    // FIND BEAM PROTONS
    HepMC::FindParticles  search_beam_p(evt, HepMC::FIND_ALL, HepMC::PDG_ID == PDG::PDG_p && HepMC::STATUS == PDG::PDG_BEAM);
    HepMC::FindParticles search_final_p(evt, HepMC::FIND_ALL, HepMC::PDG_ID == PDG::PDG_p && HepMC::STATUS == PDG::PDG_STABLE);

    if (events_read == 0) {
      //    Print::listing(evt);
      //    Print::content(evt);
    }
    
    // Mesons
    M4Vec pip;
    M4Vec pim;

    // Find out central particles
    int found = 0;
    FOREACH(const HepMC::GenParticlePtr& p1, search_pip.results()) {
      pip = gra::aux::HepMC2M4Vec(p1->momentum());
      ++found;
    }
    FOREACH(const HepMC::GenParticlePtr& p1, search_pim.results()) {
      pim = gra::aux::HepMC2M4Vec(p1->momentum());
      ++found;
    }

    // CHECK CONDITION
    if (found != 2)
      continue; // skip event, something is wrong

    // ------------------------------------------------------------------
    // Valid events

    // System
    M4Vec sys_ = pip + pim;

    // Find out beam protons
    M4Vec p_beam_plus;
    M4Vec p_beam_minus;
    M4Vec p_final_plus;
    M4Vec p_final_minus;

    if (FRAME == "GJ" || FRAME == "PG") {
      FOREACH(const HepMC::GenParticlePtr& p1, search_beam_p.results()) {
        // Print::line(p1);
        pvec = gra::aux::HepMC2M4Vec(p1->momentum());
        if (pvec.Pz() > 0) {
          p_beam_plus = pvec;
        } else {
          p_beam_minus = pvec;
        }
      }
    }

    if (FRAME == "GJ") {
      FOREACH(const HepMC::GenParticlePtr& p1, search_final_p.results()) {
        // Print::line(p1);
        pvec = gra::aux::HepMC2M4Vec(p1->momentum());
        if (pvec.Pz() > 0) {
          p_final_plus = pvec;
        } else {
          p_final_minus = pvec;
        }
      }
    }

    // ------------------------------------------------------------------
    // Do the frame transformation
    std::vector<M4Vec> boosted;
    boosted.push_back(pip);
    boosted.push_back(pim);

    // Helicity frame
    if (FRAME == "HE") {
      gra::kinematics::HEframe(boosted);
    }
    // Pseudo GJ-frame
    else if (FRAME == "GJ") {
      const uint direction = 1; // (1,-1) which proton direction
      gra::kinematics::PGframe(boosted, direction, p_beam_plus, p_beam_minus);
    }
    // Collins-Soper frame
    else if (FRAME == "CS") {
      gra::kinematics::CSframe(boosted);
    }
    // Gottfried-Jackson frame
    else if (FRAME == "GJ") {
      M4Vec propagator = p_beam_plus - p_final_plus;
      gra::kinematics::GJframe(boosted, propagator);
    }
    else {
      throw std::invalid_argument("MHarmonic::ReadIn: Unknown Lorentz frame: " + FRAME);
    }

    // ------------------------------------------------------------------
    // Construct microevent structure
    gra::spherical::Omega evt;

    // Event data
    evt.costheta = std::cos(boosted[0].Theta());
    evt.phi      = boosted[0].Phi();
    evt.M        = sys_.M();
    evt.Pt       = sys_.Pt();
    evt.Y        = sys_.Rap();
    evt.fiducial = false;
    evt.selected = false;

    // ------------------------------------------------------------------
    // FIDUCIAL CUTS
    
    if ((std::abs(pip.Eta()) < ETACUT) &&
        (std::abs(pim.Eta()) < ETACUT) &&
        pip.Pt() > PTCUT &&
        pim.Pt() > PTCUT) {
      evt.fiducial = true;
    }

    // ------------------------------------------------------------------
    // EFFICIENCY (FAST) SIMULATION

    // This block can be replaced with FULL GEANT SIMULATION / DELPHES style fast simulation
    // This is a simple parametrization to mimick pt-efficiency effects.

    // Fast simulate smooth detector pt-efficiency curve here by hyperbolic tangent per particle
    if ((std::tanh(pip.Pt() * pt_scale) > rrand.U(0,1)) && (std::tanh(pim.Pt() * pt_scale) > rrand.U(0,1)) ) {
      evt.selected = true;
    }
    // ------------------------------------------------------------------

    events[events_read] = evt;

    ++events_read;
    if (events_read % 500000 == 0) { printf("Event %d \n", events_read); }
  }

  // Remove empty memory
  events.resize(events_read);
  printf("Events read: %d \n\n", events_read);

  // Close HepMC file
  input_file.close();
}

