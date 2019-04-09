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
#include "HepMC3/FourVector.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/Print.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC3/Selector.h"
#include "HepMC3/Relatives.h"


using gra::aux::indices;

using namespace gra;

// Create Harmonic moment expansion object (needs to be global for TMinuit reasons)
MHarmonic ha;

// Random numbers for fast detector simulation
MRandom rrand;

// Final state fiducial cuts
struct FIDCUTS {

  std::vector<double> ETA = {0.0, 0.0};
  std::vector<double> PT  = {0.0, 0.0};
};

FIDCUTS fidcuts;

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
        ("r,ref",      "Reference MC (angular flat MC sample) <filename without .hepmc3>",     cxxopts::value<std::string>()  )
        ("i,input",    "Input sample MC or DATA               <filename without .hepmc3>",     cxxopts::value<std::string>()  )
        ("l,lmax",     "Maximum angular order                 <1|2|3|...>",                    cxxopts::value<unsigned int>() )
        ("f,frame",    "Lorentz rest frame                    <HE|CS|GJ|PG|SR>",               cxxopts::value<std::string>()  )
        ("c,cuts",     "Fiducial cuts                         <ETAMIN, ETAMAX, PTMIN, PTMAX>", cxxopts::value<std::string>()  )
        ("M,mass",     "System mass binning                   <bins,min,max>",                 cxxopts::value<std::string>()  )
        ("P,pt",       "System pt binning                     <bins,min,max>",                 cxxopts::value<std::string>()  )
        ("Y,rapidity", "System rapidity binning               <bins,min,max>",                 cxxopts::value<std::string>()  )
        ("X,maximum",  "Maximum number of events              <N>",                            cxxopts::value<unsigned int>() )
        ("H,help",     "Help")
        ;
    auto r = options.parse(argc, argv);
          
    if (r.count("help") || NARGC == 0) {
        std::cout << options.help({""}) << std::endl;
        std::cout << "Example:" << std::endl;
        std::cout << "  " << argv[0] << " -r SH_2pi_REF -i SH_2pi -l 4 -f HE -M 40,0.4,1.5"
                  << std::endl << std::endl;
        return EXIT_FAILURE;
    }
    
    const uint LMAX         = r["lmax"].as<unsigned int>();
    const std::string FRAME = r["frame"].as<std::string>();

    // Fiducial cuts
    if (r.count("cuts")) {
      const std::string str = r["cuts"].as<std::string>();
      std::vector<double> cuts = gra::aux::SplitStr(str,double(0),',');
      if (cuts.size() != 4)  { throw std::invalid_argument("fitharmonic:: cuts not size 4 <ETAMIN,ETAMAX,PTMIN,PTMAX>"); }
      fidcuts.ETA = {cuts[0],cuts[1]};
      fidcuts.PT  = {cuts[2],cuts[3]};     
    } else {
      throw std::invalid_argument("fitharmonic:: cuts parameter not given");
    }

    MHarmonic::HPARAM hparam;

    // Default discretization
    hparam.M  = {20, 0.0,   2.5}; // system mass
    hparam.PT = {1,  0.0, 100.0}; // system pt
    hparam.Y  = {1, -10.0, 10.0}; // system rapidity

    if (r.count("mass")) {
      const std::string Mstr = r["M"].as<std::string>();
      hparam.M = gra::aux::SplitStr(Mstr,double(0),',');
      if (hparam.M.size() != 3)  { throw std::invalid_argument("fitharmonic:: mass discretization not size 3 <bins,min,max>"); }
    }
    if (r.count("pt")) {
      const std::string Pstr = r["P"].as<std::string>();
      hparam.PT = gra::aux::SplitStr(Pstr,double(0),',');
      if (hparam.PT.size() != 3) { throw std::invalid_argument("fitharmonic:: pt discretization not size 3 <bins,min,max>"); }
    }
    if (r.count("rapidity")) {
      const std::string Ystr = r["Y"].as<std::string>();
      hparam.Y  = gra::aux::SplitStr(Ystr,double(0),',');
      if (hparam.Y.size() != 3)  { throw std::invalid_argument("fitharmonic:: rapidity discretization not size 3 <bins,min,max>"); }
    }

    // ------------------------------------------------------------------
    // INITIALIZE HARMONIC EXPANSION PARAMETERS

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
      MAXEVENTS = r["X"].as<unsigned int>();
    }
    
    // Read in reference MC
    std::vector<gra::spherical::Omega> MC_events;
    ReadIn(refinput, MC_events, FRAME, MAXEVENTS);
    
    // Read in DATA
    std::vector<gra::spherical::Omega> DATA_events;
    ReadIn(datinput, DATA_events, FRAME, MAXEVENTS);

    // ------------------------------------------------------------------
    // LOOP OVER HYPERBINS

    ha.HyperLoop(fitwrapper, MC_events, DATA_events, hparam);

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
  
  HepMC3::ReaderAscii input_file(total_input);
  // Reading failed
  if (input_file.failed()) {
    throw std::invalid_argument("MHarmonic::ReadIn: Failed to open <" + total_input + ">");
  }
  // Event loop
  int events_read = 0;

  // Dummy vector
  M4Vec pvec;
  HepMC3::GenEvent evt(HepMC3::Units::GEV, HepMC3::Units::MM);

  // Allocate memory here for speed
  events.resize((int) 1e7);
  printf("Memory allocated \n");

  // Event loop
  while (!input_file.failed()) {

    if (events_read >= MAXEVENTS) { break; }

    // Read event from input file
    input_file.read_event(evt);

    if (events_read == 0) {
      HepMC3::Print::listing(evt);
      HepMC3::Print::content(evt);
    }
    
    // Mesons
    std::vector<M4Vec> pip;
    std::vector<M4Vec> pim;

    // Find out central particles
    std::vector<HepMC3::GenParticlePtr> find_pip = 
      HepMC3::applyFilter(HepMC3::Selector::PDG_ID == PDG::PDG_pip, evt.particles());

    std::vector<HepMC3::GenParticlePtr> find_pim = 
      HepMC3::applyFilter(HepMC3::Selector::PDG_ID == PDG::PDG_pim, evt.particles());

    for (const HepMC3::GenParticlePtr& p1 : find_pip) {
      pip.push_back(gra::aux::HepMC2M4Vec(p1->momentum()));
    }
    for (const HepMC3::GenParticlePtr& p1 : find_pim) {
      pim.push_back(gra::aux::HepMC2M4Vec(p1->momentum()));
    }

    // CHECK CONDITION
    if (!(pip.size() == 1 && pim.size() == 1)) {

      //HepMC3::Print::listing(evt);
      //HepMC3::Print::content(evt);
      continue; // skip event, something is wrong
    } else {
      //printf("Found valid \n");
    }

    // ------------------------------------------------------------------
    // Valid events

    // System
    M4Vec sys_ = pip[0] + pim[0];

    // Find out beam protons
    M4Vec p_beam_plus;
    M4Vec p_beam_minus;
    M4Vec p_final_plus;
    M4Vec p_final_minus;

    // Beam protons
    if (FRAME == "GJ" || FRAME == "PG") {
    for (HepMC3::ConstGenParticlePtr p1 :
          HepMC3::applyFilter(HepMC3::Selector::STATUS == PDG::PDG_BEAM &&
                            HepMC3::Selector::PDG_ID == PDG::PDG_p, evt.particles())) {
        // Print::line(p1);
        pvec = gra::aux::HepMC2M4Vec(p1->momentum());
        if (pvec.Pz() > 0) { p_beam_plus = pvec; } else { p_beam_minus = pvec; }
      }
    }

    // Final state protons
    if (FRAME == "GJ") {
    for (HepMC3::ConstGenParticlePtr p1 :
          HepMC3::applyFilter(HepMC3::Selector::STATUS == PDG::PDG_STABLE &&
                              HepMC3::Selector::PDG_ID == PDG::PDG_p, evt.particles())) {
        // Print::line(p1);
        pvec = gra::aux::HepMC2M4Vec(p1->momentum());
        if (pvec.Pz() > 0) { p_final_plus = pvec; } else { p_final_minus = pvec; }
      }
    }

    // ------------------------------------------------------------------
    // Do the frame transformations
    std::vector<M4Vec> rf = {pip[0], pim[0]};

    // Non-rotated rest frame
    if        (FRAME == "SR") {
      gra::kinematics::SRframe(rf);
    // Helicity frame
    } else if (FRAME == "HE") {
      gra::kinematics::HEframe(rf);
    // Pseudo GJ-frame
    } else if (FRAME == "PG") {
      const unsigned int direction = 1; // (1,-1) which proton direction
      gra::kinematics::PGframe(rf, direction, p_beam_plus, p_beam_minus);
    // Collins-Soper frame
    } else if (FRAME == "CS") {
      gra::kinematics::CSframe(rf);
    // Gottfried-Jackson frame
    } else if (FRAME == "GJ") {
      M4Vec propagator = p_beam_plus - p_final_plus;
      gra::kinematics::GJframe(rf, propagator);
    }
    else {
      throw std::invalid_argument("MHarmonic::ReadIn: Unknown Lorentz frame: " + FRAME);
    }

    // ------------------------------------------------------------------
    // Construct microevent structure
    gra::spherical::Omega evt;

    // Event data
    evt.costheta = std::cos(rf[0].Theta());
    evt.phi      = rf[0].Phi();
    evt.M        = sys_.M();
    evt.Pt       = sys_.Pt();
    evt.Y        = sys_.Rap();
    evt.fiducial = false;
    evt.selected = false;

    // ------------------------------------------------------------------
    // FIDUCIAL CUTS
    
    if (pip[0].Eta() > fidcuts.ETA[0] && pip[0].Eta() < fidcuts.ETA[1] &&
        pip[0].Pt()  > fidcuts.PT[0]  && pip[0].Pt()  < fidcuts.PT[1]  &&

        pim[0].Eta() > fidcuts.ETA[0] && pim[0].Eta() < fidcuts.ETA[1] &&
        pim[0].Pt()  > fidcuts.PT[0]  && pim[0].Pt()  < fidcuts.PT[1]) {
      evt.fiducial = true;
    }

    // ------------------------------------------------------------------
    // EFFICIENCY (FAST) SIMULATION

    // This block can be replaced with FULL GEANT SIMULATION / DELPHES style fast simulation
    // This is a simple parametrization to mimick pt-efficiency effects.

    // Fast simulate smooth detector pt-efficiency curve here by hyperbolic tangent per particle
    if ((std::tanh(pip[0].Pt() * pt_scale) > rrand.U(0,1)) &&
        (std::tanh(pim[0].Pt() * pt_scale) > rrand.U(0,1)) ) {
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

