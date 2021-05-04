// GRANIITTI Monte Carlo main class
//
// (c) 2017-2021 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <algorithm>
#include <complex>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <mutex>
#include <random>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

// Own
#include "Graniitti/MAux.h"
#include "Graniitti/MContinuum.h"
#include "Graniitti/MFactorized.h"
#include "Graniitti/MGraniitti.h"
#include "Graniitti/MKinematics.h"
#include "Graniitti/MMath.h"
#include "Graniitti/MNeuroJacobian.h"
#include "Graniitti/MParton.h"
#include "Graniitti/MProcess.h"
#include "Graniitti/MQuasiElastic.h"
#include "Graniitti/MTimer.h"

// HepMC3
#include "HepMC3/GenEvent.h"
#include "HepMC3/WriterAsciiHepMC2.h"
#include "HepMC3/WriterHEPEVT.h"

// Libraries
#include "json.hpp"
#include "rang.hpp"


namespace gra {

// ******************************************************************
// GLOBALS from MGlobals.h: initialized by MGraniitti()

// Model tune
std::string MODELPARAM;

// Multithreading
std::mutex         g_mutex;
std::exception_ptr globalExceptionPtr;

// ******************************************************************

using json = nlohmann::json;

using gra::aux::indices;
using gra::math::msqrt;
using gra::math::PI;
using gra::math::pow2;
using gra::math::pow3;
using gra::math::zi;


// Constructor
MGraniitti::MGraniitti() {
  // ******************************************************************
  // Init program globals

  gra::MODELPARAM         = "null";
  gra::globalExceptionPtr = nullptr;
  // ******************************************************************

  // Print general layout
  PrintInit();
}


// Destructor
MGraniitti::~MGraniitti() {
  // Destroy processes
  for (const auto &i : indices(pvec)) { delete pvec[i]; }

  std::cout << "~MGraniitti [DONE]" << std::endl;
}


void MGraniitti::Initialize() {
  // ** ALWAYS HERE ONLY AS LAST! **
  proc->post_Constructor();

  // If Eikonals are not yet initialized and pomeron screening loop is on
  if (proc->Eikonal.IsInitialized() == false && proc->GetScreening() == true) {
    proc->Eikonal.S3Constructor(proc->GetMandelstam_s(), proc->GetInitialState(), false);
  }

  // Integrate
  CallIntegrator(0);
}


// Initialize with external Eikonal
void MGraniitti::Initialize(const MEikonal &eikonal_in) {
  // ** ALWAYS HERE ONLY AS LAST! **
  proc->post_Constructor();

  // Set input eikonal
  proc->SetEikonal(eikonal_in);

  // Integrate
  CallIntegrator(0);
}


// Save VEGAS grids to a file
void MGraniitti::SaveVGRID() const {
  // VEGAS data
  json jvegas;
  VD.struct2json(jvegas);

  // MC stats
  json jstat;
  stat.struct2json(jstat);

  // Combined
  json j;
  j["PROCESS"] = FULL_PROCESS_STR;
  j["VEGAS"]   = jvegas;
  j["STAT"]    = jstat;
  j["POMLOOP"] = proc->GetScreening();
  j["SQRTS"]   = proc->lts.sqrt_s;

  const std::string filename = gra::aux::GetBasePath(2) + "/vgrid/" + OUTPUT + ".vgrid";
  std::ofstream     file(filename);
  file << j << std::endl;

  std::cout << "MGraniitti::SaveVGRID: Saved MC integration grids to " + filename << std::endl;
}


// Read VEGAS grids from a file
void MGraniitti::ReadVGRID(const std::string &inputfile) {
  std::cout << "MGraniitti::ReadVGRID: Reading pre-computed MC integration grids from " + inputfile
            << std::endl;
  std::cout << "** Re-cycling of these grids requires that the process energy, parameters and cuts "
               "match! **"
            << std::endl;

  // Read and parse
  const std::string data = gra::aux::GetInputData(inputfile);
  json              j;
  try {
    j = json::parse(data);
  } catch (...) {
    std::string str = "MGraniitti::ReadVGRID: Error parsing " + inputfile;
    throw std::invalid_argument(str);
  }

  // Read in VEGAS data
  VD.json2struct(j.at("VEGAS"));

  // Check that the phase space dimension matches
  if (VD.FDIM != proc->GetdLIPSDim()) {
    throw std::invalid_argument(
        "MGraniitti::ReadVGRID: Input .vgrid with VEGAS phase space dim = " +
        std::to_string(VD.FDIM) + " but process with dim = " + std::to_string(proc->GetdLIPSDim()));
  }

  // Check that the screening matches
  const bool pomloop = j.at("POMLOOP");
  if (pomloop != proc->GetScreening()) {
    throw std::invalid_argument(
        "MGraniitti::ReadVGRID: Input .vgrid with POMLOOP = " + aux::bool_cast(pomloop) +
        " but process with POMLOOP = " + aux::bool_cast(proc->GetScreening()));
  }

  // Check that the CMS-energy matches
  const double sqrts = j.at("SQRTS");
  if (!aux::AssertRatio(sqrts, proc->lts.sqrt_s, 0.01)) {
    throw std::invalid_argument(
        "MGraniitti::ReadVGRID: Input .vgrid with SQRTS = " + std::to_string(sqrts) +
        " but process with SQRTS = " + std::to_string(proc->lts.sqrt_s));
  }

  // Read in MC stats
  stat.json2struct(j["STAT"]);
}


// Print fast histograms
void MGraniitti::PrintHistograms() {
  if (proc->HIST != 0) {
    HistogramFusion();

    // First print to screen (flushes buffers)
    proc->PrintHistograms();

    // Then save
    const std::string filename = gra::aux::GetBasePath(2) + "/output/" + OUTPUT + ".hfast";
    proc->SaveHistograms(filename);
  }
}


// Unify histogram bounds across threads, so the histogram fusion is possible
void MGraniitti::UnifyHistogramBounds() {
  // Loop over all 1D histograms
  for (auto const &xpoint : proc->h1) {
    int    xbins = 0;
    double xmin  = 0;
    double xmax  = 0;

    // Fuse buffer values [start index 1]
    for (std::size_t p = 1; p < pvec.size(); ++p) {
      proc->h1[xpoint.first].FuseBuffer(pvec[p]->h1[xpoint.first]);
    }

    proc->h1[xpoint.first].FlushBuffer();
    proc->h1[xpoint.first].GetBounds(xbins, xmin, xmax);
    proc->h1[xpoint.first].ResetBounds(xbins, xmin, xmax);  // New start

    // Loop over processes and set histogram bounds [start index = 1]
    for (std::size_t p = 1; p < pvec.size(); ++p) {
      pvec[p]->h1[xpoint.first].ResetBounds(xbins, xmin, xmax);
    }
  }
  // Loop over all 2D histograms
  for (auto const &xpoint : proc->h2) {
    int    xbins = 0;
    double xmin  = 0;
    double xmax  = 0;

    int    ybins = 0;
    double ymin  = 0;
    double ymax  = 0;

    // Fuse buffer values [start index 1]
    for (std::size_t p = 1; p < pvec.size(); ++p) {
      proc->h2[xpoint.first].FuseBuffer(pvec[p]->h2[xpoint.first]);
    }

    proc->h2[xpoint.first].FlushBuffer();
    proc->h2[xpoint.first].GetBounds(xbins, xmin, xmax, ybins, ymin, ymax);
    proc->h2[xpoint.first].ResetBounds(xbins, xmin, xmax, ybins, ymin, ymax);  // New start

    // Loop over processes and set histogram bounds [start index = 1]
    for (std::size_t p = 1; p < pvec.size(); ++p) {
      pvec[p]->h2[xpoint.first].ResetBounds(xbins, xmin, xmax, ybins, ymin, ymax);
    }
  }
}

// Fuse histograms for N-fold statistics
void MGraniitti::HistogramFusion() {
  if (hist_fusion_done == false) {
    // START with process index 1, because 0 is the base
    for (std::size_t i = 1; i < pvec.size(); ++i) {
      // Loop over all 1D histograms
      for (auto const &xpoint : proc->h1) {
        proc->h1[xpoint.first] = proc->h1[xpoint.first] + pvec[i]->h1[xpoint.first];
      }
      // Loop over all 2D histograms
      for (auto const &xpoint : proc->h2) {
        proc->h2[xpoint.first] = proc->h2[xpoint.first] + pvec[i]->h2[xpoint.first];
      }
    }
    hist_fusion_done = true;
  } else {
    std::cout << "MGraniitti::HistogramFusion: Multithreaded histograms have been "
                 "fused already"
              << std::endl;
  }
}

// Generate events
void MGraniitti::Generate() {
  if (NEVENTS > 0) {  // Generate events

    // Just before event generation
    // (in order not to generate unnecessary empty files
    // if file name / output type is changed)
    InitFileOutput();
    CallIntegrator(NEVENTS);
  }
}


// Return process numbers available
std::vector<std::string> MGraniitti::GetProcessNumbers() const {
  std::cout << rang::style::bold;
  std::cout << "Available processes: ";
  std::cout << rang::style::reset << std::endl << std::endl;
  std::cout << std::endl;

  std::cout << rang::style::bold << "2->3 x 1->(N-2) LIPS:" << rang::style::reset << std::endl;
  std::vector<std::string> list1 = proc_F.ProcPtr.PrintProcesses();
  std::cout << std::endl;

  std::cout << rang::style::bold << "2->N LIPS:" << rang::style::reset << std::endl;
  std::vector<std::string> list2 = proc_C.ProcPtr.PrintProcesses();
  std::cout << std::endl;

  std::cout << rang::style::bold << "Quasielastic:" << rang::style::reset << std::endl;
  std::vector<std::string> list3 = proc_Q.ProcPtr.PrintProcesses();
  std::cout << std::endl;

  std::cout << rang::style::bold << "2->1 x (1->M) collinear:" << rang::style::reset << std::endl;
  std::vector<std::string> list4 = proc_P.ProcPtr.PrintProcesses();
  std::cout << std::endl;

  // Concatenate all
  list1.insert(list1.end(), list2.begin(), list2.end());
  list1.insert(list1.end(), list3.begin(), list3.end());
  list1.insert(list1.end(), list4.begin(), list4.end());

  return list1;
}

// (Re-)assign the pointers to the local memory space
void MGraniitti::InitProcessMemory(std::string process, unsigned int seed) {
  // These must be here!
  PROCESS = process;

  // <Q> process
  if (proc_Q.ProcPtr.ProcessExist(process)) {
    proc_Q = MQuasiElastic(process, syntax);
    proc   = &proc_Q;

    // <F> processes
  } else if (proc_F.ProcPtr.ProcessExist(process)) {
    proc_F = MFactorized(process, syntax);
    proc   = &proc_F;

    // <C> processes
  } else if (proc_C.ProcPtr.ProcessExist(process)) {
    proc_C = MContinuum(process, syntax);
    proc   = &proc_C;

    // <P> processes
  } else if (proc_P.ProcPtr.ProcessExist(process)) {
    proc_P = MParton(process, syntax);
    proc   = &proc_P;
  } else {
    std::string str = "MGraniitti::InitProcessMemory: Unknown PROCESS: " + process;
    // GetProcessNumbers(); // This is done by main program
    throw std::invalid_argument(str);
  }

  // Set random seed last!
  proc->random.SetSeed(seed);
}


void MGraniitti::InitMultiMemory() {
  // ** Init multithreading memory by making copies of the process **
  for (const auto &i : indices(pvec)) { delete pvec[i]; }

  pvec.resize(CORES, nullptr);

  // Create new process objects for each thread
  for (const auto &i : indices(pvec)) {
    if (proc_Q.ProcPtr.ProcessExist(PROCESS)) {
      pvec[i] = new MQuasiElastic(proc_Q);
    } else if (proc_F.ProcPtr.ProcessExist(PROCESS)) {
      pvec[i] = new MFactorized(proc_F);
    } else if (proc_C.ProcPtr.ProcessExist(PROCESS)) {
      pvec[i] = new MContinuum(proc_C);
    } else if (proc_P.ProcPtr.ProcessExist(PROCESS)) {
      pvec[i] = new MParton(proc_P);
    }
  }

  // RANDOM SEED PER THREAD (IMPORTANT!)
  for (const auto &i : indices(pvec)) {
    // Deterministic seed
    const uint32_t thread_seed = static_cast<uint32_t>(pvec[i]->random.GetSeed() + i * 10);
    pvec[i]->random.SetSeed(thread_seed);

    printf("Thread [%3lu] : RNG seed = %9d \n", i, pvec[i]->random.GetSeed());
  }
  std::cout << std::endl;

  // SET main control pointer to the first one
  // (needed for printing etc.)
  proc = pvec[0];
  proc->PrintInit(HILJAA);
}


// Set simple MC parameters
void MGraniitti::SetMCParam(MCPARAM &in) { mcparam = in; }

// Set VEGAS parameters
void MGraniitti::SetVegasParam(const VEGASPARAM &in) { vparam = in; }

// Read parameters from a single JSON file
void MGraniitti::ReadInput(const json &j) {
  ReadGeneralParam(j);
  ReadProcessParam(j);
  ReadModelParam(gra::MODELPARAM);
}

// General parameter initialization
void MGraniitti::ReadGeneralParam(const json &j) {
  // JSON block identifier
  const std::string XID = "GENERALPARAM";

  // Setup parameters (order is important)
  SetNumberOfEvents(j.at(XID).at("NEVENTS"));
  SetOutput(j.at(XID).at("OUTPUT"));
  SetFormat(j.at(XID).at("FORMAT"));
  SetWeighted(j.at(XID).at("WEIGHTED"));
  SetIntegrator(j.at(XID).at("INTEGRATOR"));
  SetCores(j.at(XID).at("CORES"));

  // Save for later use
  gra::MODELPARAM = j.at(XID).at("MODELPARAM");
}

// General model parameters initialized from .json file
//
void MGraniitti::ReadModelParam(const std::string &tune) const {
  const std::string fullpath = gra::aux::GetBasePath(2) + "/modeldata/" + tune + "/GENERAL.json";

  // Read generic blocks
  PARAM_SOFT::ReadParameters(fullpath);
  PARAM_STRUCTURE::ReadParameters(fullpath);
  PARAM_FLAT::ReadParameters(fullpath);
  PARAM_NSTAR::ReadParameters(fullpath);

  // The rest are handled by spesific amplitude classes
}

// Process parameter initialization, Call proc->post_Constructor() after this
void MGraniitti::ReadProcessParam(const json &j) {
  const std::string XID = "PROCESSPARAM";

  // ----------------------------------------------------------------
  // Initialize process

  std::string fullstring = j.at(XID).at("PROCESS");
  FULL_PROCESS_STR       = fullstring;  // Save it for later

  // ----------------------------------------------------------------
  // First separate possible extra arguments by @... ...
  std::vector<std::size_t> markerpos = aux::FindOccurance(fullstring, "@");

  if (markerpos.size() != 0) {
    syntax     = aux::SplitCommands(fullstring.substr(markerpos[0]));
    fullstring = fullstring.substr(0, markerpos[0] - 1);
  }
  // ----------------------------------------------------------------

  // Now separate process and decay parts by "->"
  std::string PROCESS_PART = "";
  std::string DECAY_PART   = "";
  std::size_t pos          = 0;
  std::size_t pos1         = fullstring.find("->");
  std::size_t pos2         = fullstring.find("&>");  // ISOLATEd Phase-Space

  pos = (pos1 != std::string::npos) ? pos1 : std::string::npos;  // Try to find ->
  if (pos == std::string::npos) {
    pos = (pos2 != std::string::npos) ? pos2 : std::string::npos;  // Try to find &>
  }

  if (pos != std::string::npos) {
    PROCESS_PART = fullstring.substr(0, pos - 1);  // beginning to the pos-1
    DECAY_PART   = fullstring.substr(pos + 2);     // from pos+2 to the end
  } else {
    PROCESS_PART = fullstring;  // No decay defined
  }

  // Trim extra spaces away
  gra::aux::TrimExtraSpace(PROCESS_PART);
  gra::aux::TrimExtraSpace(DECAY_PART);

  InitProcessMemory(PROCESS_PART, j.at(XID).at("RNDSEED"));

  // ----------------------------------------------------------------
  // SETUP process: process memory needs to be initialized before!

  proc->SetLHAPDF(j.at(XID).at("LHAPDF"));
  proc->SetScreening(j.at(XID).at("POMLOOP"));
  proc->SetExcitation(j.at(XID).at("NSTARS"));
  proc->SetHistograms(j.at(XID).at("HIST"));

  if (pos2 != std::string::npos) {
    if (PROCESS_PART.find("<F>") == std::string::npos &&
        PROCESS_PART.find("<P>") == std::string::npos) {
      throw std::invalid_argument(
          "MGraniitti::ReadProcessParam: Phase space "
          "isolation arrow '&>' to be "
          "used only with <F> or <P> class!");
    }
    proc->SetISOLATE(true);
  }

  // ----------------------------------------------------------------
  // Decaymode setup
  proc->SetDecayMode(DECAY_PART);

  // Read free resonance parameters (only if needed)
  std::map<std::string, gra::PARAM_RES> RESONANCES;
  if (PROCESS_PART.find("RES") != std::string::npos) {
    // From .json input
    std::vector<std::string> RES = j.at(XID).at("RES");

    // @ Command syntax override
    for (const auto &i : indices(syntax)) {
      if (syntax[i].id == "RES") {
        std::vector<std::string> temp;
        for (const auto &x : syntax[i].arg) {
          if (x.second == "true" || x.second == "1") {
            temp.push_back(x.first);  // f0_980, ...
          }
        }
        RES = temp;
      }
    }

    // Read resonance data
    for (const auto &i : indices(RES)) {
      const std::string str = "RES_" + RES[i] + ".json";
      RESONANCES[RES[i]]    = gra::form::ReadResonance(str, proc->random);
    }

    // @ Command syntax override "on-the-flight parameters"
    for (const auto &i : indices(syntax)) {
      if (syntax[i].id == "R") {
        // Take target string
        std::string RESNAME;
        if (syntax[i].target.size() == 1) {
          RESNAME = syntax[i].target[0];
        } else if (syntax[i].target.size() == 0) {
          throw std::invalid_argument("@Syntax error: invalid R[] without any target []");
        } else if (syntax[i].target.size() > 1) {
          throw std::invalid_argument("@Syntax error: invalid R[] with multiple targets inside []");
        }

        // Do we find target
        if (RESONANCES.find(RESNAME) == RESONANCES.end()) {
          throw std::invalid_argument("@Syntax error: invalid R[" + RESNAME + "] not found");
        }
        bool couplings_touched = false;
        bool density_touched   = false;

        MMatrix<std::complex<double>> newrho(1, 1, 0.0);
        if (RESONANCES[RESNAME].p.spinX2 != 0) {
          const int N = RESONANCES[RESNAME].rho.size_row();
          newrho      = MMatrix<std::complex<double>>(N, N, 0.0);
        }

        // Loop over key:val arguments
        for (const auto &x : syntax[i].arg) {
          if (x.first == "M") {
            RESONANCES[RESNAME].p.mass = std::stod(x.second);
            std::cout << rang::fg::green << "@R[" << RESNAME
                      << "] new mass: " << RESONANCES[RESNAME].p.mass << rang::fg::reset
                      << std::endl;
          }
          if (x.first == "W") {
            RESONANCES[RESNAME].p.width = std::stod(x.second);
            std::cout << rang::fg::green << "@R[" << RESNAME
                      << "] new width: " << RESONANCES[RESNAME].p.width << rang::fg::reset
                      << std::endl;
          }

          // Set couplings
          for (std::size_t n = 0; n < RESONANCES[RESNAME].g_Tensor.size(); ++n) {
            if (x.first == ("g" + std::to_string(n))) {
              RESONANCES[RESNAME].g_Tensor[n] = std::stod(x.second);
              couplings_touched               = true;
            }
          }

          // Set diagonal spin density elements
          for (std::size_t n = 0; n <= (newrho.size_row() - 1) / 2; ++n) {
            const int J = (newrho.size_row() - 1) / 2;
            if (x.first == ("JZ" + std::to_string(n))) {
              const double value = std::stod(x.second);

              if (value < 0) {
                throw std::invalid_argument("MGraniitti::ReadProcessParam: @R[] JZ value < 0");
              }

              newrho[newrho.size_row() - 1 - J - n][newrho.size_row() - 1 - J - n] = value;

              newrho[newrho.size_row() - 1 - J + n][newrho.size_row() - 1 - J + n] = value;
              density_touched                                                      = true;
            }
          }
        }

        // Normalize the trace and save it
        if (density_touched) {
          const std::complex<double> trace = newrho.Trace();
          if (std::abs(trace) > 0) {
            newrho = newrho * (1.0 / trace);
          } else {
            throw std::invalid_argument(
                "MGraniitti::ReadProcessParam: @R[] Spin density error with <= 0 trace");
          }
          RESONANCES[RESNAME].rho = newrho;

          std::cout << rang::fg::green << "@R[" << RESNAME
                    << "] new spin density set:" << std::endl;
          RESONANCES[RESNAME].rho.Print();
          std::cout << rang::fg::reset << std::endl;
        }

        // Print new coupling array
        if (couplings_touched) {
          std::cout << rang::fg::green << "@R[" << RESNAME
                    << "] new production couplings set: g_Tensor = [";
          for (std::size_t k = 0; k < RESONANCES[RESNAME].g_Tensor.size(); ++k) {
            printf("%0.3E", RESONANCES[RESNAME].g_Tensor[k]);
            if (k < RESONANCES[RESNAME].g_Tensor.size() - 1) { std::cout << ", "; }
          }
          std::cout << "]" << rang::fg::reset << std::endl;
        }
      }
    }
    proc->SetResonances(RESONANCES);

    // Setup resonance branching (final step!)
    proc->SetupBranching();
  }
  // ----------------------------------------------------------------

  // Command syntax parameters
  for (const auto &i : indices(syntax)) {
    if (syntax[i].id == "FLATAMP") { proc->SetFLATAMP(std::stoi(syntax[i].arg["_SINGLET_"])); }
    if (syntax[i].id == "FLATMASS2") { proc->SetFLATMASS2(syntax[i].arg["_SINGLET_"] == "true"); }
    if (syntax[i].id == "OFFSHELL") { proc->SetOFFSHELL(std::stod(syntax[i].arg["_SINGLET_"])); }

    if (syntax[i].id == "SPINGEN") { proc->SetSPINGEN(syntax[i].arg["_SINGLET_"] == "true"); }
    if (syntax[i].id == "SPINDEC") { proc->SetSPINDEC(syntax[i].arg["_SINGLET_"] == "true"); }

    if (syntax[i].id == "FRAME") { proc->SetFRAME(syntax[i].arg["_SINGLET_"]); }
    if (syntax[i].id == "JMAX") { proc->SetJMAX(std::stoi(syntax[i].arg["_SINGLET_"])); }
  }

  // Always last (we need PDG tables)
  std::vector<std::string> beam   = j.at(XID).at("BEAM");
  std::vector<double>      energy = j.at(XID).at("ENERGY");
  proc->SetInitialState(beam, energy);

  // Now rest of the parameters
  ReadIntegralParam(j);
  ReadGenCuts(j);
  ReadFidCuts(j);
  ReadVetoCuts(j);
}

// MC integrator parameter initialization
void MGraniitti::ReadIntegralParam(const json &j) {
  const std::string XID = "INTEGRALPARAM";

  // Numerical loop integration
  const int ND = j.at(XID).at("POMLOOP").at("ND");
  gra::aux::AssertRange(ND, {-10, 10}, "POMLOOP::ND", true);
  proc->Eikonal.Numerics.SetLoopDiscretization(ND);

  // FLAT (naive) MC parameters
  MCPARAM mpam;
  mpam.PRECISION = j.at(XID).at("FLAT").at("PRECISION");
  gra::aux::AssertRange(mpam.PRECISION, {0.0, 1.0}, "FLAT::PRECISION", true);
  mpam.MIN_EVENTS = j.at(XID).at("FLAT").at("MIN_EVENTS");
  gra::aux::AssertRange(mpam.MIN_EVENTS, {10, (unsigned int)1e9}, "FLAT::MIN_EVENTS", true);
  SetMCParam(mpam);

  try {
    j.at(XID).at("VEGAS").at("BINS");
  } catch (...) {
    std::cout << "MGraniitti::ReadIntegralParam: Did not found VEGAS parameter block "
                 "from json input, using default"
              << std::endl;
    return;  // Did not find VEGAS parameter block at all, use default
  }

  // VEGAS parameters
  VEGASPARAM vpam;
  vpam.BINS = j.at(XID).at("VEGAS").at("BINS");
  gra::aux::AssertRange(vpam.BINS, {0, (unsigned int)1e9}, "VEGAS::BINS", true);
  if ((vpam.BINS % 2) != 0) {
    throw std::invalid_argument("VEGAS::BINS = " + std::to_string(vpam.BINS) +
                                " should be even number");
  }

  vpam.LAMBDA = j.at(XID).at("VEGAS").at("LAMBDA");
  gra::aux::AssertRange(vpam.LAMBDA, {1.0, 10.0}, "VEGAS::LAMBDA", true);
  vpam.NCALL = j.at(XID).at("VEGAS").at("NCALL");
  gra::aux::AssertRange(vpam.NCALL, {50, (unsigned int)1e9}, "VEGAS::NCALL", true);
  vpam.ITER = j.at(XID).at("VEGAS").at("ITER");
  gra::aux::AssertRange(vpam.ITER, {1, (unsigned int)1e9}, "VEGAS::ITER", true);
  vpam.CHI2MAX = j.at(XID).at("VEGAS").at("CHI2MAX");
  gra::aux::AssertRange(vpam.CHI2MAX, {0.0, 1e3}, "VEGAS::CHI2MAX", true);
  vpam.PRECISION = j.at(XID).at("VEGAS").at("PRECISION");
  gra::aux::AssertRange(vpam.PRECISION, {0.0, 1.0}, "VEGAS::PRECISION", true);
  vpam.DEBUG = j.at(XID).at("VEGAS").at("DEBUG");
  gra::aux::AssertSet(vpam.DEBUG, {-1, 0, 1}, "VEGAS::DEBUG", true);

  SetVegasParam(vpam);
}

// Generator cuts
void MGraniitti::ReadGenCuts(const json &j) {
  gra::GENCUT       gcuts;
  const std::string XID = "GENCUTS";

  // Continuum phase space class
  if (PROCESS.find("<C>") != std::string::npos) {
    // Daughter rapidity
    std::vector<double> rap;
    try {
      std::vector<double> temp = j.at(XID).at("<C>").at("Rap");
      rap                      = temp;
    } catch (...) {
      throw std::invalid_argument(
          "MGraniitti::ReadGenCuts: <C> phase space class requires from user: "
          "\"GENCUTS\" : { \"<C>\" : { \"Rap\" : [min, max] }}");
    }
    gcuts.rap_min = rap[0];
    gcuts.rap_max = rap[1];
    gra::aux::AssertCut(rap, "GENCUTS::<C>::Rap", true);

    // This is optional, intermediate kt
    std::vector<double> kt;
    try {
      std::vector<double> temp = j.at(XID).at("<C>").at("Kt");
      kt                       = temp;
    } catch (...) {
      // Do nothing
    }
    if (kt.size() != 0) {
      gcuts.kt_min = kt[0];
      gcuts.kt_max = kt[1];
      gra::aux::AssertCutRange(kt, {0.0, 1e32}, "GENCUTS::<C>::Kt", true);
    }

    // This is optional, forward leg pt
    std::vector<double> pt;
    try {
      std::vector<double> temp = j.at(XID).at("<C>").at("Pt");
      pt                       = temp;
    } catch (...) {
      // Do nothing
    }
    if (pt.size() != 0) {
      gcuts.forward_pt_min = pt[0];
      gcuts.forward_pt_max = pt[1];
      gra::aux::AssertCutRange(pt, {0.0, 1e32}, "GENCUTS::<C>::Pt", true);
    }

    // This is optional, forward leg excitation
    std::vector<double> Xi;
    try {
      std::vector<double> temp = j.at(XID).at("<C>").at("Xi");
      Xi                       = temp;
    } catch (...) {
      // Do nothing
    }
    if (Xi.size() != 0) {
      gcuts.XI_min = Xi[0];
      gcuts.XI_max = Xi[1];
      gra::aux::AssertCutRange(Xi, {0.0, 1.0}, "GENCUTS::<C>::Xi", true);
    }
  }

  // Factorized phase space class
  if (PROCESS.find("<F>") != std::string::npos) {
    // System rapidity
    std::vector<double> Y;
    try {
      std::vector<double> temp = j.at(XID).at("<F>").at("Rap");
      Y                        = temp;
    } catch (...) {
      throw std::invalid_argument(
          "MGraniitti::ReadGenCuts: <F> phase space class requires from user: "
          "\"GENCUTS\" : { \"<F>\" : { \"Rap\" : [min, max] }}");
    }
    gcuts.Y_min = Y[0];
    gcuts.Y_max = Y[1];
    gra::aux::AssertCut(Y, "GENCUTS::<F>::Rap", true);

    // System mass
    std::vector<double> M;
    try {
      std::vector<double> temp = j.at(XID).at("<F>").at("M");
      M                        = temp;
    } catch (...) {
      throw std::invalid_argument(
          "MGraniitti::ReadGenCuts: <F> phase space class requires from user: "
          "\"GENCUTS\" : { \"<F>\" : { \"M\" : [min, max] }}");
    }
    gcuts.M_min = M[0];
    gcuts.M_max = M[1];
    gra::aux::AssertCutRange(M, {0.0, 1e32}, "GENCUTS::<F>::M", true);

    // This is optional, forward leg pt
    std::vector<double> pt;
    try {
      std::vector<double> temp = j.at(XID).at("<F>").at("Pt");
      pt                       = temp;
    } catch (...) {
      // Do nothing
    }
    if (pt.size() != 0) {
      gcuts.forward_pt_min = pt[0];
      gcuts.forward_pt_max = pt[1];
      gra::aux::AssertCutRange(pt, {0.0, 1e32}, "GENCUTS::<F>::Pt", true);
    }

    // This is optional, forward leg excitation
    std::vector<double> Xi;
    try {
      std::vector<double> temp = j.at(XID).at("<F>").at("Xi");
      Xi                       = temp;
    } catch (...) {
      // Do nothing
    }
    if (Xi.size() != 0) {
      gcuts.XI_min = Xi[0];
      gcuts.XI_max = Xi[1];
      gra::aux::AssertCutRange(Xi, {0.0, 1.0}, "GENCUTS::<F>::Xi", true);
    }
  }

  // Quasielastic phase space class
  if (PROCESS.find("<Q>") != std::string::npos) {
    std::vector<double> Xi;
    try {
      std::vector<double> temp = j.at(XID).at("<Q>").at("Xi");
      Xi                       = temp;
    } catch (...) {
      throw std::invalid_argument(
          "MGraniitti::ReadGenCuts: <Q> phase space class requires from user: "
          "\"GENCUTS\" : { \"<Q>\" : { \"Xi\" : [min, max] }}");
    }
    gcuts.XI_min = Xi[0];
    gcuts.XI_max = Xi[1];
    gra::aux::AssertCutRange(Xi, {0.0, 1.0}, "GENCUTS::<Q>::Xi", true);
  }

  proc->SetGenCuts(gcuts);
}

// Fiducial cuts
void MGraniitti::ReadFidCuts(const json &j) {
  // Fiducial cuts
  gra::FIDCUT       fcuts;
  const std::string XID = "FIDCUTS";

  try {
    fcuts.active = j.at(XID).at("active");
  } catch (...) {
    return;  // Did not find cut block at all
  }

  // Particles
  {
    std::vector<double> eta = j.at(XID).at("PARTICLE").at("Eta");
    fcuts.eta_min           = eta[0];
    fcuts.eta_max           = eta[1];
    gra::aux::AssertCut(eta, "FIDCUTS::PARTICLE::Eta", true);

    std::vector<double> rap = j.at(XID).at("PARTICLE").at("Rap");
    fcuts.rap_min           = rap[0];
    fcuts.rap_max           = rap[1];
    gra::aux::AssertCut(rap, "FIDCUTS::PARTICLE::Rap", true);

    std::vector<double> pt = j.at(XID).at("PARTICLE").at("Pt");
    fcuts.pt_min           = pt[0];
    fcuts.pt_max           = pt[1];
    gra::aux::AssertCut(pt, "FIDCUTS::PARTICLE::Pt", true);

    std::vector<double> Et = j.at(XID).at("PARTICLE").at("Et");
    fcuts.Et_min           = Et[0];
    fcuts.Et_max           = Et[1];
    gra::aux::AssertCut(Et, "FIDCUTS::PARTICLE::Et", true);
  }

  // System
  {
    std::vector<double> M = j.at(XID).at("SYSTEM").at("M");
    fcuts.M_min           = M[0];
    fcuts.M_max           = M[1];
    gra::aux::AssertCut(M, "FIDCUTS::SYSTEM::M", true);

    std::vector<double> Y = j.at(XID).at("SYSTEM").at("Rap");
    fcuts.Y_min           = Y[0];
    fcuts.Y_max           = Y[1];
    gra::aux::AssertCut(Y, "FIDCUTS::SYSTEM::Rap", true);

    std::vector<double> Pt = j.at(XID).at("SYSTEM").at("Pt");
    fcuts.Pt_min           = Pt[0];
    fcuts.Pt_max           = Pt[1];
    gra::aux::AssertCut(Pt, "FIDCUTS::SYSTEM::Pt", true);
  }

  // Forward
  {
    std::vector<double> t = j.at(XID).at("FORWARD").at("t");
    fcuts.forward_t_min   = t[0];
    fcuts.forward_t_max   = t[1];
    gra::aux::AssertCut(t, "FIDCUTS::FORWARD::t", true);

    std::vector<double> M = j.at(XID).at("FORWARD").at("M");
    fcuts.forward_M_min   = M[0];
    fcuts.forward_M_max   = M[1];
    gra::aux::AssertCut(M, "FIDCUTS::FORWARD::M", true);
  }

  // Set fiducial cuts
  proc->SetFidCuts(fcuts);

  // Set user cuts
  proc->SetUserCuts(j.at(XID).at("USERCUTS"));
}

// Fiducial cuts
void MGraniitti::ReadVetoCuts(const json &j) {
  // Fiducial cuts
  gra::VETOCUT      veto;
  const std::string XID = "VETOCUTS";

  try {
    veto.active = j.at(XID).at("active");
  } catch (...) {
    return;  // Did not find cut block at all
  }

  // Find domains
  for (std::size_t i = 0; i < 100; ++i) {
    const std::string NUMBER = std::to_string(i);

    gra::VETODOMAIN     domain;
    std::vector<double> eta;
    std::vector<double> pt;

    // try to find new domain
    try {
      std::vector<double> temp1 = j.at(XID).at(NUMBER).at("Eta");
      eta                       = temp1;
      std::vector<double> temp2 = j.at(XID).at(NUMBER).at("Pt");
      pt                        = temp2;
    } catch (...) {
      break;  // Not found
    }

    domain.eta_min = eta[0];
    domain.eta_max = eta[1];
    gra::aux::AssertCut(eta, "VETOCUTS::" + NUMBER + "::Eta", true);

    domain.pt_min = pt[0];
    domain.pt_max = pt[1];
    gra::aux::AssertCut(pt, "VETOCUTS::" + NUMBER + "::Pt", true);

    veto.cuts.push_back(domain);
  }

  // Set fiducial cuts
  proc->SetVetoCuts(veto);
}

// Get maximum integration weight
double MGraniitti::GetMaxweight() const {
  if (INTEGRATOR == "VEGAS") {
    return stat.maxf;
  } else {
    return stat.maxW;
  }
}

// Set maximum weight for the integration process
void MGraniitti::SetMaxweight(double weight) {
  if (weight > 0) {
    if (INTEGRATOR == "VEGAS") {
      stat.maxf = weight;
    } else {
      stat.maxW = weight;
    }
  } else {
    std::string str =
        "MGraniitti::SetMaxweight: Maximum weight: " + std::to_string(weight) + " not valid!";
    throw std::invalid_argument(str);
  }
}


void MGraniitti::PrintInit() const {
  if (!HILJAA) {
    gra::aux::PrintFlashScreen(rang::fg::green);
    std::cout << rang::style::bold
              << "GRANIITTI - Monte Carlo event generator for "
                 "high energy diffraction"
              << rang::style::reset << std::endl
              << std::endl;
    gra::aux::PrintVersion();
    gra::aux::PrintBar("-");

    const double GB = pow3(1024.0);
    printf("Running on %s (%d CORE / %0.2f GB RAM) at %s \n", gra::aux::HostName().c_str(),
           std::thread::hardware_concurrency(), gra::aux::TotalSystemMemory() / GB,
           gra::aux::DateTime().c_str());
    int64_t size = 0;
    int64_t free = 0;
    int64_t used = 0;
    gra::aux::GetDiskUsage("/", size, free, used);
    printf("Path '/': size | used | free = %0.1f | %0.1f | %0.1f GB \n", size / GB, used / GB,
           free / GB);
    std::cout << "Program path: " << gra::aux::GetBasePath(2) << std::endl;
    std::cout << gra::aux::SystemName() << std::endl;
    gra::aux::PrintBar("-");
  }
}


void MGraniitti::CallIntegrator(unsigned int N) {
  // Print out basic information
  std::cout << std::endl;
  std::cout << rang::style::bold << "General setup:" << rang::style::reset << std::endl
            << std::endl;

  std::cout << "Output file:            " << OUTPUT << std::endl;
  std::cout << "Output format:          " << FORMAT << std::endl;
  std::cout << "Multithreading:         " << CORES << std::endl;
  std::cout << "Integrator:             " << INTEGRATOR << std::endl;
  std::cout << "Number of events:       " << NEVENTS << std::endl;
  std::cout << "Parameter setup:        " << gra::MODELPARAM << std::endl;

  std::string str = (WEIGHTED == true) ? "weighted" : "unweighted";
  std::cout << rang::fg::green << "Generation mode:        " << str << rang::fg::reset << std::endl;
  std::cout << std::endl;


  // =====================================================================

  // Initialize global clock
  if (N == 0) { global_tictoc = MTimer(true); }

  // Sample the phase space
  if (INTEGRATOR == "VEGAS") {
    SampleVegas(N);
  } else if (INTEGRATOR == "FLAT") {
    SampleFlat(N);
  } else if (INTEGRATOR == "NEURO") {
    SampleNeuro(N);
  } else {
    throw std::invalid_argument("MGraniitti::CallIntegrator: Unknown INTEGRATOR = " + INTEGRATOR +
                                " (use VEGAS, FLAT, NEURO)");
  }
}

// Initialize and generate events using VEGAS MC
void MGraniitti::SampleVegas(unsigned int N) {
  if (N == 0) {
    InitMultiMemory();
    GMODE = 0;  // Pure integration

    // ******************************************************
    // Create integral sampling region array with given phase space dimension
    VD.InitRegion(proc->GetdLIPSDim());
    // ******************************************************
  }
  if (N > 0) {
    GMODE = 1;  // Event generation
  }

  // A. Initalize from a pre-computed file
  if (MC_VGRID_INPUT != "null" && GMODE == 0) {
    ReadVGRID(MC_VGRID_INPUT);

    // Print pre-computed statistics
    PrintStatistics(N);
  }

  // B. On-the-flight initialization
  else if (GMODE == 0) {
    const double MINTIME     = 0.1;  // Seconds
    unsigned int BURNIN_ITER = 3;    // BURN-IN iterations (default)!

    // Initialize GRID on-the-flight
    do {
      unsigned int init   = 0;
      int          factor = 0;

      do {  // Loop until stable
        factor       = VEGAS(init, vparam.NCALL, BURNIN_ITER, N);
        vparam.NCALL = vparam.NCALL * factor;
        if (factor == 1) { break; }  // We are done
      } while (true);

      // Increase CALLS and re-run, if too fast
      if (itertime < MINTIME) {
        const double time_per_iter = itertime / vparam.NCALL;

        // Max, because this is only minimum condition
        vparam.NCALL =
            std::max((unsigned int)vparam.NCALL, (unsigned int)(MINTIME / time_per_iter));
        VEGAS(init, vparam.NCALL, BURNIN_ITER, N);
      }

      // Now re-calculate by skipping burn-in (init = 0) iterations
      // because they detoriate the integral value by bad grid
      init   = 1;
      factor = VEGAS(init, vparam.NCALL, vparam.ITER, N);

      if (factor == 1) {
        break;
      } else {
        BURNIN_ITER = 2 * BURNIN_ITER;
      }
    } while (true);

    // Save VEGAS grid
    SaveVGRID();
  }

  // Event generation mode
  else if (GMODE == 1) {
    const unsigned int init    = 2;
    const unsigned int itermin = 1E9;
    VEGAS(init, vparam.NCALL * 10, itermin, N);
  }

  else {
    throw std::invalid_argument("MGraniitti::SampleVegas: Unknown operation mode!");
  }
}


// Multithreaded VEGAS integrator
// [close to optimal importance sampling iff
//  integrand factorizes dimension by dimension]
//
// Original algorithm from:
// [REFERENCE: Lepage, G.P. Journal of Computational Physics, 1978]
// en.wikipedia.org/wiki/VEGAS_algorithm

int MGraniitti::VEGAS(unsigned int init, unsigned int calls, unsigned int itermin, unsigned int N) {
  // First the initialization
  VD.Init(init, vparam);

  if (init == 0 && !HILJAA) {
    gra::aux::ClearProgress();
    std::cout << rang::style::bold;
    printf("VEGAS burn-in iterations: \n\n");
    std::cout << rang::style::reset;
  }

  if (init == 1 && !HILJAA) {
    gra::aux::ClearProgress();  // Progressbar clear
    gra::aux::PrintBar("-");
    std::cout << rang::style::bold;
    std::cout << "VEGAS importance sampling:" << std::endl << std::endl;
    std::cout << rang::style::reset;
    printf("<Multithreading CORES = %d> \n\n", CORES);

    printf("- dimension = %u \n", VD.FDIM);
    printf(" \n");
    printf("- itermin   = %d \n", itermin);
    printf("- calls     = %d \n", calls);
    printf("- BINS      = %u \n", vparam.BINS);
    printf("- PRECISION = %0.4f \n", vparam.PRECISION);
    printf("- CHI2MAX   = %0.4f \n", vparam.CHI2MAX);
    printf("- LAMBDA    = %0.4f \n", vparam.LAMBDA);
    printf("- DEBUG     = %d \n", vparam.DEBUG);
    gra::aux::PrintBar("-");
  }

  // Reset local timers
  local_tictoc = MTimer(true);

  MTimer stime = MTimer(true);  // For statusprint
  atime        = MTimer(true);  // For progressbar

  // -------------------------------------------------------------------
  // Create number of calls per thread, their sum == calls
  const std::vector<unsigned int> LOCALcalls = VD.GetLocalCalls(calls, CORES);


  MTimer gridtic;

  // VEGAS grid iterations
  for (std::size_t iter = 0; iter < itermin; ++iter) {
    if (init == 0 && iter == 2) {           // Save time for one iteration
      itertime = gridtic.ElapsedSec() / 3;  // Average over 3 iter
    }

    for (std::size_t j = 0; j < VD.FDIM; ++j) {
      for (std::size_t i = 0; i < vparam.BINS; ++i) {
        VD.fmat[i][j]  = 0.0;
        VD.f2mat[i][j] = 0.0;
      }
    }
    // Init here outside parallel processing
    VD.fsum  = 0.0;
    VD.f2sum = 0.0;

    // --------------------------------------------------------------
    // SPAWN PARALLEL PROCESSING HERE

    std::vector<std::thread> threads;
    for (int tid = 0; tid < CORES; ++tid) {
      threads.push_back(
          std::thread([=, this] { VEGASMultiThread(N, tid, init, LOCALcalls[tid]); }));
    }

    // Loop again to join the threads once they complete
    for (auto &t : threads) { t.join(); }
    if (gra::globalExceptionPtr) {  // Exception handling of threads
      std::rethrow_exception(gra::globalExceptionPtr);
    }

    // --------------------------------------------------------------
    // Estimates based on this iteration

    VD.fsum *= 1.0 / static_cast<double>(calls);
    VD.f2sum *= pow2(1.0 / static_cast<double>(calls));

    // Local integral estimate
    const double I_this = VD.fsum;

    VD.f2sum = msqrt(VD.f2sum * calls);
    VD.f2sum = (VD.f2sum - VD.fsum) * (VD.f2sum + VD.fsum);  // - +
    if (VD.f2sum <= 0.0) { VD.f2sum = vparam.EPS; }

    // Local error squared estimate
    const double I2_this = VD.f2sum * 1.0 / static_cast<double>(calls);

    // --------------------------------------------------------------
    // Update global estimates
    const double alpha = 1.0 / I2_this;

    VD.sumdata += alpha * I_this;
    VD.sumchi2 += alpha * pow2(I_this);
    VD.sweight += alpha;

    // Integral and its error estimate
    stat.sigma      = VD.sumdata / VD.sweight;
    stat.sigma_err  = msqrt(1.0 / VD.sweight);
    stat.sigma_err2 = pow2(stat.sigma_err);

    // Chi2
    double chi2this = (VD.sumchi2 - pow2(VD.sumdata) / VD.sweight) / (iter + 1E-5);
    chi2this        = (chi2this < 0) ? 0.0 : chi2this;
    stat.chi2       = chi2this;
    // --------------------------------------------------------------

    // Got enough events generated
    if (GMODE == 1 && stat.generated >= N) { goto stop; }

    // Fatal error and convergence restart treatment
    if (GMODE == 0 && iter > 0) {
      if (std::isnan(stat.sigma) || std::isinf(stat.sigma)) {
        gra::aux::ClearProgress();
        throw std::invalid_argument(
            "VEGAS:: Integral inf/nan: FATAL ERROR, cuts or parameters "
            "need fixing. Otherwise, increase 'NCALL' or try <F> phase space class if not in use.");
      }
      if (stat.sigma < 1e-60) {
        gra::aux::ClearProgress();
        printf("VEGAS:: integral = %0.1f < 1e-60: Convergence problem, increasing 10 x calls. \n",
               stat.sigma);
        printf("Try <F> phase space class if not in use. \n");
        return 10;  // 10 times more calls
      }
      if (stat.chi2 > 50) {
        gra::aux::ClearProgress();
        printf("VEGAS:: chi2 = %0.1f > 50: Convergence problem, increasing 10 x calls. \n",
               stat.chi2);
        printf("Try <F> phase space class if not in use. \n");
        return 10;  // 10 times more calls
      }
    }

    if (vparam.DEBUG >= 0) {
      gra::aux::ClearProgress();
      printf(
          "VEGAS:: local iter = %4lu integral = %0.5E +- std = "
          "%0.5E \t [global integral = %0.5E +- std = %0.5E] \t "
          "chi2this = %0.2f \n",
          iter + 1, I_this, gra::math::msqrt(I2_this), stat.sigma, stat.sigma_err, chi2this);

      if (vparam.DEBUG > 0) {
        for (std::size_t j = 0; j < VD.FDIM; ++j) {
          printf("VEGAS:: data for dimension j = %lu (VD.FDIM = %u) \n", j, VD.FDIM);

          for (std::size_t i = 0; i < vparam.BINS; ++i) {
            printf(
                "xmat[%3lu][j] = %0.5E, fmat[%3lu][j] = %0.5E "
                "\n",
                i, VD.xmat[i][j], i, VD.fmat[i][j]);
          }
        }
      }
    }

    // We are not below the chi2 or precision condition -> increase iterations
    if (init > 0) {  // do not consider in burn-in (init = 0) mode
      if ((iter == (itermin - 1) && stat.chi2 > vparam.CHI2MAX) ||
          (iter == (itermin - 1) && stat.sigma_err / stat.sigma > vparam.PRECISION)) {
        ++itermin;
      }
    }

    // Status and grid optimization
    if (GMODE == 0) {
      if (stime.ElapsedSec() > 2.0 || init == 0) {
        PrintStatus(stat.evaluations, N, local_tictoc, -1.0);
        stime.Reset();
      }
      if (atime.ElapsedSec() > 0.01) {
        gra::aux::PrintProgress((iter + 1) / static_cast<double>(itermin));
        atime.Reset();
      }

      // ==========================================================
      // **** Update VEGAS grid (not during event generation) ****
      VD.OptimizeGrid(vparam);
      // ==========================================================
    }

    // Unify histogram boundaries after burn-in across different threads
    // (due to adaptive histogramming)
    if (init == 0 && iter == itermin - 1) { UnifyHistogramBounds(); }
  }  // Main grid iteration loop

stop:  // We jump here once finished (GOTO point)

  // Final statistics
  if (init > 0) { PrintStatistics(N); }

  return 1;  // Return 1 for good
}

// This is called once for every VEGAS grid iteration
void MGraniitti::VEGASMultiThread(unsigned int N, unsigned int THREAD_ID, unsigned int init,
                                  unsigned int LOCALcalls) {
  double zn = 0.0;
  double zo = 0.0;
  double ac = 0.0;

  try {
    for (std::size_t k = 0; k < LOCALcalls; ++k) {  // ** LOCALcalls ~= calls / CORES **

      double vegasweight = 1.0;

      // Phase space point vector
      std::vector<double> xpoint(pvec[THREAD_ID]->GetdLIPSDim(), 0.0);

      // Loop over dimensions and construct random vector xpoint
      std::vector<double> indvec(VD.FDIM, 0);
      for (std::size_t j = 0; j < VD.FDIM; ++j) {
        // Draw random number
        zn = pvec[THREAD_ID]->random.U(0, 1) * vparam.BINS + 1.0;
        indvec[j] =
            std::max((unsigned int)1, std::min((unsigned int)zn, (unsigned int)vparam.BINS));

        if (indvec[j] > 1) {
          zo = VD.xmat[indvec[j] - 1][j] - VD.xmat[indvec[j] - 2][j];
          ac = VD.xmat[indvec[j] - 2][j] + (zn - indvec[j]) * zo;
        } else {
          zo = VD.xmat[indvec[j] - 1][j];
          ac = (zn - indvec[j]) * zo;
        }

        // Multidim space vector component
        xpoint[j] = VD.region[j] + ac * VD.dxvec[j];
        vegasweight *= zo * vparam.BINS;

      }  // VD.FDIM loop

      // *******************************************************************
      // ****** Call the process under integration to get the weight *******

      gra::AuxIntData aux;
      aux.vegasweight  = vegasweight;
      aux.burn_in_mode = (init == 0) ? true : false;

      const double W = pvec[THREAD_ID]->EventWeight(xpoint, aux);

      // *******************************************************************

      // @@ Multithreading lock (FAST SECTION) @@
      gra::g_mutex.lock();

      // Increase statistics
      stat.Accumulate(aux);

      // *** Importance weighting ***
      const double f  = W * vegasweight;
      const double f2 = pow2(f);

      VD.fsum += f;
      VD.f2sum += f2;

      // Loop over dimensions and add importance weighted results
      for (std::size_t j = 0; j < VD.FDIM; ++j) {
        VD.fmat[indvec[j] - 1][j] += f;
        VD.f2mat[indvec[j] - 1][j] += f2;
      }

      // ----------------------------------------------------------
      // Initialization (integration) mode
      if (GMODE == 0) {
        // Do not consider burn-in phase weights (unstable)
        if (init != 0) {
          // Maximum raw weight (for general information, not used
          // here)
          if (W > stat.maxW) { stat.maxW = W; }
          // Maximum total VEGAS importance weighted
          if (f > stat.maxf) { stat.maxf = f; }
        }
      }
      gra::g_mutex.unlock();
      // @@ Multithreading unlock (FAST SECTION) @@

      // ----------------------------------------------------------
      // Event generation mode
      if (GMODE == 1) {
        // Enough events
        if (stat.generated == (unsigned int)GetNumberOfEvents()) { break; }

        // Event trial
        SaveEvent(pvec[THREAD_ID], f, stat.maxf, aux);

        if (THREAD_ID == 0 && atime.ElapsedSec() > 0.5) {
          PrintStatus(stat.generated, N, local_tictoc, 10.0);
          gra::aux::PrintProgress(stat.generated / static_cast<double>(N));
          atime.Reset();
        }
      }
    }  // calls loop
  } catch (...) {
    // Set the global exception pointer if exception arises
    // This is because of multithreading
    gra::globalExceptionPtr = std::current_exception();
  }
}

// Generate events using plain simple MC (for reference/DEBUG purposes)
void MGraniitti::SampleFlat(unsigned int N) {
  // Integration mode
  if (N == 0) {
    GMODE = 0;
    proc->PrintInit(HILJAA);
  }
  // Event generation mode
  if (N > 0) { GMODE = 1; }

  // Get dimension of the phase space
  const unsigned int  dim = proc->GetdLIPSDim();
  std::vector<double> randvec(dim, 0.0);


  // Reset local timer
  local_tictoc = MTimer(true);

  // Progressbar
  atime = MTimer(true);

  // Reservation (test)
  // using namespace std::placeholders; // _1, _2, ... come from here
  // std::function<double(const std::vector<double>&, int, double,
  // std::vector<double>&,
  // bool&)> fifth =
  //  std::bind(&MQuasiElastic::EventWeight, &proc_Q, _1, _2, _3, _4, _5);

  // Event loop
  while (true) {
    // ** Generate new random numbers **
    for (const auto &i : indices(randvec)) { randvec[i] = proc->random.U(0, 1); }

    // Generate event
    // Used for in-out control of the process
    gra::AuxIntData aux;
    aux.vegasweight  = 1.0;
    aux.burn_in_mode = false;
    const double W   = proc->EventWeight(randvec, aux);

    // Increase statistics
    stat.Accumulate(aux);
    stat.Wsum += W;
    stat.W2sum += pow2(W);

    // Update cross section estimate
    stat.CalculateCrossSection();

    // Initialization
    if (GMODE == 0) {
      stat.maxW = (W > stat.maxW) ? W : stat.maxW;  // Update maximum weight
      PrintStatus(stat.evaluations, N, local_tictoc, 10.0);

      if ((stat.sigma_err / stat.sigma) < mcparam.PRECISION &&
          stat.evaluations > mcparam.MIN_EVENTS) {
        break;
      }

      // Progressbar
      if (atime.ElapsedSec() > 0.1) {
        gra::aux::PrintProgress(stat.evaluations / static_cast<double>(mcparam.MIN_EVENTS));
        atime.Reset();
      }
    }

    // Event generation mode
    if (GMODE == 1) {
      SaveEvent(proc, W, stat.maxW, aux);
      PrintStatus(stat.generated, N, local_tictoc, 10.0);
      if (stat.generated >= N) { break; }

      // Progressbar
      if (atime.ElapsedSec() > 0.1) {
        gra::aux::PrintProgress(stat.generated / static_cast<double>(N));
        atime.Reset();
      }
    }
  }
  PrintStatus(stat.generated, N, local_tictoc, -1.0);
  PrintStatistics(N);
}

// Neural net Monte Carlo (small scale PROTOTYPE)
void MGraniitti::SampleNeuro(unsigned int N) {
  // Integration mode
  if (N == 0) {
    GMODE = 0;
    proc->PrintInit(HILJAA);
  }
  // Event generation mode
  if (N > 0) { GMODE = 1; }

  // Get dimension of the phase space
  const unsigned int D = proc->GetdLIPSDim();

  // Reset timers
  local_tictoc = MTimer(true);
  atime        = MTimer(true);

  // Reservation (test)
  // using namespace std::placeholders; // _1, _2, ... come from here
  // std::function<double(const std::vector<double>&, int, double,
  // std::vector<double>&,
  // bool&)> fifth =
  // std::bind(&MQuasiElastic::EventWeight, &proc_Q, _1, _2, _3, _4, _5);

  if (N == 0) {
    // -------------------------------------------------------
    // Bind the object and call it
    gra::neurojac::procptr = proc;  // First set address
    // -------------------------------------------------------

    gra::neurojac::MNeuroJacobian neurojac;
    gra::neurojac::BATCHSIZE = 100;

    // Set network layer dimensions [first, ..., output]
    gra::neurojac::par.D = D;  // Integrand dimension

    gra::neurojac::par.L.push_back(gra::neurojac::Layer(2, D));  // Input
    gra::neurojac::par.L.push_back(gra::neurojac::Layer(2, 2));
    gra::neurojac::par.L.push_back(gra::neurojac::Layer(2, 2));
    gra::neurojac::par.L.push_back(gra::neurojac::Layer(D, 2));  // Output

    neurojac.Optimize();
  }

  // -------------------------------------------------------------------
  // Now to the event generation

  // Lambda capture
  std::vector<double> u(D);

  auto NeuroSample = [&](gra::AuxIntData &aux) {
    // Prior p(z) distribution sampling
    VectorXdual z(u.size());
    for (std::size_t i = 0; i < D; ++i) { z[i] = proc->random.G(0, 1); }
    const double p = val(gra::neurojac::gaussprob(z, 0, 1));

    // Evaluate network map
    VectorXdual u_ = gra::neurojac::G_net(z);

    // Evaluate the Jacobian matrix du/dz
    MatrixXd dudz = jacobian(gra::neurojac::G_net, u_, z);

    // Abs Jacobian determinant and inverse prior
    const double jacweight = abs(dudz.determinant()) / p;

    for (std::size_t i = 0; i < D; ++i) { u[i] = val(u_[i]); }

    // Evaluate event weight
    aux.vegasweight     = jacweight;
    aux.burn_in_mode    = false;
    const double weight = proc->EventWeight(u, aux) * jacweight;

    return weight;
  };

  // Event loop
  while (true) {
    // aux used for in-out control of the process
    gra::AuxIntData aux;
    const double    W = NeuroSample(aux);

    // Increase statistics
    stat.Accumulate(aux);
    stat.Wsum += W;
    stat.W2sum += pow2(W);

    // Update cross section estimate
    stat.CalculateCrossSection();

    // Initialization
    if (GMODE == 0) {
      stat.maxW = (W > stat.maxW) ? W : stat.maxW;  // Update maximum weight
      PrintStatus(stat.evaluations, N, local_tictoc, 10.0);

      if ((stat.sigma_err / stat.sigma) < mcparam.PRECISION &&
          stat.evaluations > mcparam.MIN_EVENTS) {
        break;
      }

      // Progressbar
      if (atime.ElapsedSec() > 0.1) {
        gra::aux::PrintProgress(stat.evaluations / static_cast<double>(mcparam.MIN_EVENTS));
        atime.Reset();
      }
    }

    // Event generation mode
    if (GMODE == 1) {
      SaveEvent(proc, W, stat.maxW, aux);
      PrintStatus(stat.generated, N, local_tictoc, 10.0);
      if (stat.generated >= N) { break; }

      // Progressbar
      if (atime.ElapsedSec() > 0.1) {
        gra::aux::PrintProgress(stat.generated / static_cast<double>(N));
        atime.Reset();
      }
    }
  }
  PrintStatus(stat.generated, N, local_tictoc, -1.0);
  PrintStatistics(N);
}

// Save unweighted or weighted event
int MGraniitti::SaveEvent(MProcess *pr, double weight, double MAXWEIGHT,
                          const gra::AuxIntData &aux) {
  gra::g_mutex.lock();
  stat.trials += 1;  // This is one trial more

  // 0. Hit-Miss
  bool hit_in = proc->random.U(0, 1) < (weight / MAXWEIGHT);

  // 2. We see weight larger than maxweight
  if (!WEIGHTED && weight > MAXWEIGHT) {
    ++stat.N_overflow;
    hit_in = false;
  }
  gra::g_mutex.unlock();

  // Three ways to accept the event. N.B. weight > 0 needed if the amplitude
  // fails numerically
  if ((hit_in && aux.Valid()) || (WEIGHTED && aux.Valid()) || aux.forced_accept) {
    // Create HepMC3 event (do not lock yet for speed)
    HepMC3::GenEvent evt(HepMC3::Units::GEV, HepMC3::Units::MM);

    // ** Construct event record **
    if (!pr->EventRecord(evt)) {  // Event not ok!
      // std::cout << "MGraniitti::SaveEvent: Last moment rare veto!" <<
      // std::endl;
      return 2;
    }

    // Event ok, continue >>

    // @@@ THIS IS THREAD-NON-SAFE -> LOCK IT
    gra::g_mutex.lock();

    // ** This is a multithreading race-condition treatment **
    if (stat.generated == (unsigned int)GetNumberOfEvents()) {
      stat.trials -= 1;  // Correct statistics, unnecessary trial
      gra::g_mutex.unlock();
      return 1;
    }

    // Set event number
    evt.set_event_number(stat.generated);


    // ** Save event weight (unweighted events with weight 1)
    const double HepMC3_weight = WEIGHTED ? weight : 1.0;
    evt.weights().push_back(HepMC3_weight);  // add more weights with .push_back()


    // ** Save cross section information (HepMC3 wants event by event)
    std::shared_ptr<HepMC3::GenCrossSection> xsobj = std::make_shared<HepMC3::GenCrossSection>();
    evt.add_attribute("GenCrossSection", xsobj);

    // Now add the value in picobarns [HepMC3 convention]
    if (xsforced > 0) {
      xsobj->set_cross_section(xsforced * OUTPUT_XS_SCALE, 0);  // external fixed one
    } else {
      xsobj->set_cross_section(stat.sigma * OUTPUT_XS_SCALE, stat.sigma_err * OUTPUT_XS_SCALE);
    }

    // ** Save PDF info
    std::shared_ptr<HepMC3::GenPdfInfo> pdf = std::make_shared<HepMC3::GenPdfInfo>();
    evt.add_attribute("GenPdfInfo", pdf);
    pdf->set(pr->lts.id1, pr->lts.id2, pr->lts.x1, pr->lts.x2, pr->lts.scale, pr->lts.xf1,
             pr->lts.xf2, 0, 0);


    if (FORMAT == "hepmc3") {
      outputHepMC3->write_event(evt);
    } else if (FORMAT == "hepmc2") {
      outputHepMC2->write_event(evt);
    } else if (FORMAT == "hepevt") {
      outputHEPEVT->write_event(evt);
    } else {
      throw std::invalid_argument("MGraniitti::SaveEvent: Unknown output FORMAT " + FORMAT);
    }

    // LAST STEP
    stat.generated += 1;  // +1 event generated

    // Cleanup up the event from memory
    evt.clear();

    gra::g_mutex.unlock();
    // @@ THIS IS THREAD-NON-SAFE <- LOCK IT @@
    return 0;
  } else {
    return 1;
  }
}

// This is called just before event generation
// Check against nullptr is done below, because if the output is already set
// externally,
// this function is not used.
void MGraniitti::InitFileOutput() {
  if (NEVENTS > 0) {
    if (OUTPUT == "") {  // OUTPUT must be set
      throw std::invalid_argument("MGraniitti::InitFileOutput: OUTPUT filename not set!");
    }

    FULL_OUTPUT_STR = gra::aux::GetBasePath(2) + "/output/" + OUTPUT + "." + FORMAT;

    // --------------------------------------------------------------
    // Generator info
    runinfo = std::make_shared<HepMC3::GenRunInfo>();

    struct HepMC3::GenRunInfo::ToolInfo generator = {
        std::string("GRANIITTI (" + gra::MODELPARAM + ")"),
        std::to_string(aux::GetVersion()).substr(0, 5), std::string("Generator")};
    runinfo->tools().push_back(generator);

    struct HepMC3::GenRunInfo::ToolInfo config = {FULL_PROCESS_STR, "1.0",
                                                  std::string("Steering card")};
    runinfo->tools().push_back(config);


    // ** Cross section is added also event by event later **
    // Now add the value in picobarns [HepMC3 convention]
    if (xsforced > 0) {
      // external fixed one
      runinfo->add_attribute("xs",
                             std::make_shared<HepMC3::FloatAttribute>(xsforced * OUTPUT_XS_SCALE));
      runinfo->add_attribute("xs_err", std::make_shared<HepMC3::FloatAttribute>(0));
    } else {
      runinfo->add_attribute(
          "xs", std::make_shared<HepMC3::FloatAttribute>(stat.sigma * OUTPUT_XS_SCALE));
      runinfo->add_attribute(
          "xs_err", std::make_shared<HepMC3::FloatAttribute>(stat.sigma_err * OUTPUT_XS_SCALE));
    }

    // --------------------------------------------------------------

    if (FORMAT == "hepmc3" && outputHepMC3 == nullptr) {
      outputHepMC3 = std::make_shared<HepMC3::WriterAscii>(FULL_OUTPUT_STR, runinfo);
    } else if (FORMAT == "hepmc2" && outputHepMC2 == nullptr) {
      outputHepMC2 = std::make_shared<HepMC3::WriterAsciiHepMC2>(FULL_OUTPUT_STR, runinfo);
    } else if (FORMAT == "hepevt" && outputHEPEVT == nullptr) {
      outputHEPEVT = std::make_shared<HepMC3::WriterHEPEVT>(FULL_OUTPUT_STR);
    }
  }
}


// Intermediate statistics
void MGraniitti::PrintStatus(unsigned int events, unsigned int N, MTimer &tictoc, double timercut) {
  if (tictoc.ElapsedSec() > timercut) {
    tictoc.Reset();
    gra::aux::ClearProgress();

    double       peak_use     = 0.0;
    double       resident_use = 0.0;
    const double MB           = 1024 * 1024;
    const double GB           = MB * 1024;
    gra::aux::GetProcessMemory(peak_use, resident_use);
    peak_use /= MB;
    resident_use /= MB;

    if (GMODE == 0) {
      const double global_lap = global_tictoc.ElapsedSec();
      printf(
          "[%0.1f MB] xs: %9.3E, er: %7.5f [chi2 = %2.1f], %4.1f "
          "min ~ %0.1E Hz \n",
          resident_use, stat.sigma, stat.sigma_err / stat.sigma, stat.chi2, global_lap / 60.0,
          events / global_lap);
    }
    if (GMODE == 1) {
      const double global_lap     = global_tictoc.ElapsedSec() - time_t0;
      double       outputfilesize = gra::aux::GetFileSize(FULL_OUTPUT_STR) / GB;

      printf(
          "[%0.1f MB/%0.2f GB] E: %9d, xs: %9.3E, er: %7.5f, %0.1f/%0.1f min ~ "
          "%0.1E Hz \n",
          resident_use, outputfilesize, events, stat.sigma, stat.sigma_err / stat.sigma,
          global_lap / 60.0, (N - events) * global_lap / (double)events / 60.0,
          events / global_lap);
    }
  }
}


// Final statistics
void MGraniitti::PrintStatistics(unsigned int N) {
  gra::aux::ClearProgress();  // Clear progressbar

  if (GMODE == 0) {
    time_t0 = global_tictoc.ElapsedSec();
    gra::aux::PrintBar("=");
    std::cout << rang::style::bold;
    printf("Monte Carlo integration: \n\n\n");
    std::cout << rang::style::reset;

    if (proc->GetISOLATE()) {
      std::cout << rang::fg::red << "NOTE: Central leg phase space isolation tag &> in use!";
      std::cout << rang::fg::reset;
      std::cout << std::endl << std::endl;
    }

    // Check if we have cascaded phase space turned on in the x-section calculation
    double prod    = 1.0;
    double prod2pi = 1.0;
    double volume  = 1.0;
    int    Nf      = 0;
    for (const auto &i : indices(proc->lts.decaytree)) {
      proc->CalculatePhaseSpace(proc->lts.decaytree[i], prod, prod2pi, volume, Nf);
    }
    unsigned int N_leg = std::max((int)proc->lts.decaytree.size(), Nf) + 2;  // +2 forward legs

    // Special cases
    if (proc->GetISOLATE()) { N_leg = 3; }        // ISOLATEd 2->3 process with <F> phase space
    if (proc->GetdLIPSDim() == 2) { N_leg = 2; }  // <P> and <Q> class

    printf("{2->%d cross section}:             [%0.3E +- %0.3E] barn \n", N_leg, stat.sigma,
           stat.sigma_err);
    printf("Sampling uncertainty:             %0.3f %% \n", 100 * stat.sigma_err / stat.sigma);
    std::cout << std::endl;

    if (proc->lts.DW_sum.Integral() > 0) {  // Recursive phase space on

      std::cout << std::endl;
      const double PSvolume = proc->lts.DW_sum.Integral();
      // const double PSvolume_error       = proc->lts.DW_sum.IntegralError();

      const double PSvolume_exact       = proc->lts.DW_sum_exact.Integral();
      const double PSvolume_exact_error = proc->lts.DW_sum_exact.IntegralError();

      // Get sum of decay daughter masses
      double MSUM = 0.0;
      for (std::size_t i = 0; i < proc->lts.decaytree.size(); ++i) {
        MSUM += proc->lts.decaytree[i].m_offshell;
      }
      // only 2-body exact or massless case
      if (proc->lts.decaytree.size() == 2 || MSUM < 1e-6) {
        printf("Analytic phase space volume:      %0.3E +- %0.3E \n", PSvolume_exact,
               PSvolume_exact_error);
        printf("RATIO: MC/analytic:               %0.6f \n", PSvolume / PSvolume_exact);
      }
      std::cout << std::endl;

      // Print out phase space weight
      printf("{1->%lu LIPS}:                      %0.3E +- %0.3E  ", proc->lts.decaytree.size(),
             proc->lts.DW_sum.Integral(), proc->lts.DW_sum.IntegralError());

      if (proc->GetISOLATE()) {
        std::cout << "[" << rang::fg::red << "INACTIVE " << rang::fg::reset
                  << "part of integral  <=> apply manual BR]" << std::endl;
      } else {
        std::cout << "[" << rang::fg::green << "ACTIVE " << rang::fg::reset
                  << "part of integral / (2PI)]" << std::endl;
      }
      // Recursion relation based (phase space factorization):
      // d^N PS(s; p_1, p2, ...p_N)
      // = 1/(2*PI) * d^3 PS(s; p1,p2,pX) d^{N-2} PS(M^2; pX,p3,p4,..,pN) dM^2
      //
      // Decaywidth = 1/(2M S) \int dPS |M_decay|^2, where M
      // = mother mass, S = final state symmetry factor

      if (!proc->GetISOLATE() && proc->GetdLIPSDim() != 2) {  // Special cases
        printf("\n\n{2->3 cross section ~=~ [2->%u / (1->%lu LIPS) x 2PI]}:       %0.3E \n", N_leg,
               proc->lts.decaytree.size(), stat.sigma / proc->lts.DW_sum.Integral() * (2 * PI));
        std::cout << std::endl;
        printf("** Remember to use &> operator instead of -> for phase space isolation ** \n\n");
      }
    }

    // Print recursively
    double product    = 1.0;
    double product2pi = 1.0;
    int    N_final    = 2;  // forward legs == 2
    for (const auto &i : indices(proc->lts.decaytree)) {
      proc->PrintPhaseSpace(proc->lts.decaytree[i], product, product2pi, N_final);
      if (proc->lts.decaytree[i].legs.size() != 0) { std::cout << std::endl; }
    }
    gra::aux::PrintBar("-");

    std::cout << std::endl;
    printf("Integration runtime:              %0.2f sec \n", time_t0);
    printf("Integrand sampling frequency:     %0.2E Hz \n",
           (stat.evaluations - stat.trials) / time_t0);
    std::cout << std::endl;
    printf("<Statistics> \n");
    printf("MAX raw integrand weight:         %0.3E \n", stat.maxW);
    if (INTEGRATOR == "VEGAS") { printf("MAX weight (vegas x integrand):   %0.3E \n", stat.maxf); }

    printf("\n");
    printf("Amplitude passing rate:           %0.3E \n", stat.amplitude_ok / stat.evaluations);
    printf("Kinematics passing rate:          %0.3E \n", stat.kinematics_ok / stat.evaluations);
    printf("Fiducial cuts passing rate:       %0.3E \n", stat.fidcuts_ok / stat.evaluations);
    printf("Veto cuts passing rate:           %0.3E \n", stat.vetocuts_ok / stat.evaluations);
    printf("\n");

    std::cout << std::endl;
    printf(
        "** All values include phase space generation and "
        "fiducial (+ veto) cuts ** \n");
    gra::aux::PrintBar("=");
  }
  if (GMODE == 1) {
    const double lap = global_tictoc.ElapsedSec() - time_t0;

    gra::aux::PrintBar("=");
    if (WEIGHTED) {
      gra::aux::PrintNotice();
      std::cout << rang::fg::red << "You did WEIGHTED event generation:" << std::endl
                << std::endl
                << std::endl;
    } else {
      std::cout << rang::fg::green
                << "You did UNWEIGHTED (acceptance-rejection) event generation:" << std::endl
                << std::endl
                << std::endl;
    }
    std::cout << rang::style::reset;

    printf("Generation efficiency:    %0.3E (%d / %0.0f) \n", N / (double)stat.trials, N,
           stat.trials);
    printf("Weight overflow:          %0.3E (%d / %0.0f) \n", stat.N_overflow / (double)stat.trials,
           stat.N_overflow, stat.trials);
    printf("Generation runtime:       %0.2f sec \n", lap);
    printf("Generation frequency:     %0.2E Hz \n", N / lap);

    double outputfilesize = gra::aux::GetFileSize(FULL_OUTPUT_STR) / (1024.0 * 1024.0 * 1024.0);

    printf("Outputfile size:          %0.3f GB [%s] \n", outputfilesize, FULL_OUTPUT_STR.c_str());
    gra::aux::PrintBar("=");
    std::cout << std::endl;
  }
}


// Construct terminal input parameters
void MGraniitti::ConstructTerminal(cxxopts::Options &options) const {
  options.add_options("GENERALPARAM")("o,OUTPUT", "Output name             <string>",
                                      cxxopts::value<std::string>())(
      "f,FORMAT", "Output format           <hepmc3|hepmc2|hepevt>", cxxopts::value<std::string>())(
      "c,CORES", "Number of CPU threads   <integer>", cxxopts::value<unsigned int>())(
      "n,NEVENTS", "Number of events        <integer32>", cxxopts::value<unsigned int>())(
      "g,INTEGRATOR", "Integrator              <VEGAS|FLAT>", cxxopts::value<std::string>())(
      "w,WEIGHTED", "Weighted events         <true|false>", cxxopts::value<std::string>())(
      "m,MODELPARAM", "Model tune              <string>", cxxopts::value<std::string>());

  options.add_options("PROCESSPARAM")("p,PROCESS", "Process                 <string>",
                                      cxxopts::value<std::string>())(
      "e,ENERGY", "CMS energy              <double>", cxxopts::value<double>())(
      "l,POMLOOP", "Screening Pomeron loop  <true|false>", cxxopts::value<std::string>())(
      "s,NSTARS", "Excite protons          <0|1|2>", cxxopts::value<unsigned int>())(
      "q,LHAPDF", "Set LHAPDF              <string>", cxxopts::value<std::string>())(
      "h,HIST", "Histogramming           <0|1|2>", cxxopts::value<unsigned int>())(
      "r,RNDSEED", "Random seed             <integer32>", cxxopts::value<unsigned int>());
}


// Override json object parameters from the command line object r
void MGraniitti::ProcessTerminal(json &j, cxxopts::ParseResult const &r) const {
  // GENERALPARAM

  if (r.count("n")) { j.at("GENERALPARAM").at("NEVENTS") = r["n"].as<unsigned int>(); }
  if (r.count("o")) { j.at("GENERALPARAM").at("OUTPUT") = r["o"].as<std::string>(); }
  if (r.count("f")) { j.at("GENERALPARAM").at("FORMAT") = r["f"].as<std::string>(); }
  if (r.count("m")) { j.at("GENERALPARAM").at("MODELPARAM") = r["m"].as<std::string>(); }
  if (r.count("g")) { j.at("GENERALPARAM").at("INTEGRATOR") = r["g"].as<std::string>(); }
  if (r.count("w")) { j.at("GENERALPARAM").at("WEIGHTED") = (r["w"].as<std::string>() == "true"); }
  if (r.count("c")) { j.at("GENERALPARAM").at("CORES") = r["c"].as<unsigned int>(); }


  // PROCESSPARAM

  if (r.count("p")) { j.at("PROCESSPARAM").at("PROCESS") = r["p"].as<std::string>(); }
  if (r.count("l")) { j.at("PROCESSPARAM").at("POMLOOP") = (r["l"].as<std::string>() == "true"); }
  if (r.count("e")) {
    j.at("PROCESSPARAM").at("ENERGY") = std::vector<double>(2, r["e"].as<double>() / 2.0);
  }
  if (r.count("s")) { j.at("PROCESSPARAM").at("NSTARS") = r["s"].as<unsigned int>(); }
  if (r.count("q")) { j.at("PROCESSPARAM").at("LHAPDF") = r["q"].as<std::string>(); }
  if (r.count("h")) { j.at("PROCESSPARAM").at("HIST") = r["h"].as<unsigned int>(); }
  if (r.count("r")) { j.at("PROCESSPARAM").at("RNDSEED") = r["r"].as<unsigned int>(); }
}

}  // namespace gra
