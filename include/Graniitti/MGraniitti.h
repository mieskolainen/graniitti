// GRANIITTI Monte Carlo main class
//
// (c) 2017-2021 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MGRANIITTI_H
#define MGRANIITTI_H

// C++
#include <complex>
#include <mutex>
#include <random>
#include <thread>
#include <vector>

// HepMC3
#include "HepMC3/GenEvent.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC3/WriterAsciiHepMC2.h"
#include "HepMC3/WriterHEPEVT.h"

// Own
#include "Graniitti/M4Vec.h"
#include "Graniitti/MAux.h"
#include "Graniitti/MContinuum.h"
#include "Graniitti/MEikonal.h"
#include "Graniitti/MFactorized.h"
#include "Graniitti/MKinematics.h"
#include "Graniitti/MMatrix.h"
#include "Graniitti/MParton.h"
#include "Graniitti/MProcess.h"
#include "Graniitti/MQuasiElastic.h"
#include "Graniitti/MSpin.h"
#include "Graniitti/MTimer.h"
#include "Graniitti/MVEGAS.h"

// Other
#include "cxxopts.hpp"
#include "json.hpp"

namespace gra {

// Simple MC parameters
struct MCPARAM {
  double       PRECISION  = 0.05;     // Integral relative precision
  unsigned int MIN_EVENTS = 1000000;  // Minimum number of events to be sampled
};

// Integration statistics
class Stats {
 public:
  void Accumulate(const gra::AuxIntData &aux) {
    evaluations += 1.0;

    amplitude_ok += (aux.amplitude_ok ? 1.0 : 0.0);
    kinematics_ok += (aux.kinematics_ok ? 1.0 : 0.0);
    fidcuts_ok += (aux.fidcuts_ok ? 1.0 : 0.0);
    vetocuts_ok += (aux.vetocuts_ok ? 1.0 : 0.0);
  }

  // Calculate integrated cross section sigma and its
  // error for direct (simple) sampling. NOT TO BE USED WITH VEGAS!
  void CalculateCrossSection() {
    // <f>
    sigma = Wsum / evaluations;
    // err^2 = <f^2> - <f>^2
    sigma_err2 = W2sum / evaluations - gra::math::pow2(sigma);
    // err_trials = err^2 / sqrt(trials) (standard error of the mean)
    sigma_err = gra::math::msqrt(sigma_err2 / evaluations);
  }

  void struct2json(json &j) const {
    j["amplitude_ok"]  = amplitude_ok;
    j["kinematics_ok"] = kinematics_ok;
    j["fidcuts_ok"]    = fidcuts_ok;
    j["vetocuts_ok"]   = vetocuts_ok;

    j["evaluations"] = evaluations;
    j["trials"]      = trials;

    j["generated"]  = generated;
    j["N_overflow"] = N_overflow;

    j["sigma"]      = sigma;
    j["sigma_err"]  = sigma_err;
    j["sigma_err2"] = sigma_err2;

    j["Wsum"]  = Wsum;
    j["W2sum"] = W2sum;
    j["maxW"]  = maxW;

    j["maxf"] = maxf;
    j["chi2"] = chi2;
  }

  void json2struct(json &j) {
    j.at("amplitude_ok").get_to(amplitude_ok);
    j.at("kinematics_ok").get_to(kinematics_ok);
    j.at("fidcuts_ok").get_to(fidcuts_ok);
    j.at("vetocuts_ok").get_to(vetocuts_ok);

    j.at("evaluations").get_to(evaluations);
    j.at("trials").get_to(trials);

    j.at("generated").get_to(generated);
    j.at("N_overflow").get_to(N_overflow);

    j.at("sigma").get_to(sigma);
    j.at("sigma_err").get_to(sigma_err);
    j.at("sigma_err2").get_to(sigma_err2);

    j.at("Wsum").get_to(Wsum);
    j.at("W2sum").get_to(W2sum);
    j.at("maxW").get_to(maxW);

    j.at("maxf").get_to(maxf);
    j.at("chi2").get_to(chi2);
  }


  double amplitude_ok  = 0.0;
  double kinematics_ok = 0.0;
  double fidcuts_ok    = 0.0;
  double vetocuts_ok   = 0.0;

  // Keep as double to avoid overflow of range
  double evaluations = 0.0;  // Integrand evaluations
  double trials      = 0.0;  // Event generation trials

  unsigned int generated  = 0.0;  // Event generation
  unsigned int N_overflow = 0.0;  // Weight overflows

  // Cross section and its error
  double sigma      = 0.0;
  double sigma_err  = 0.0;
  double sigma_err2 = 0.0;

  // Weight statistics FLAT MC
  double Wsum  = 0.0;
  double W2sum = 0.0;
  double maxW  = 0.0;

  // Weight statistics VEGAS MC
  double maxf = 0.0;
  double chi2 = 0.0;
};

class MGraniitti {
 public:
  // Destructor & Constructor
  MGraniitti();
  ~MGraniitti();

  // -------------------------------------------------------------------
  // Copy, assignment and move disabled
  MGraniitti(const MGraniitti &) = delete;
  MGraniitti &operator=(const MGraniitti &) = delete;
  MGraniitti(MGraniitti &&)                 = delete;
  MGraniitti &operator=(MGraniitti &&) = delete;
  // -------------------------------------------------------------------

  // For cross-section for HepMC3 output
  void ForceXS(double xs) { xsforced = xs; }

  // Read parameters
  void ReadInput(const json &j);

  // Initialize memory
  void InitProcessMemory(std::string process, unsigned int RNDSEED);
  void InitMultiMemory();


  // Set simple MC parameters
  void SetMCParam(MCPARAM &in);

  // Set VEGAS parameters
  void SetVegasParam(const VEGASPARAM &in);

  // Set external file handle
  void SetHepMC2Output(std::shared_ptr<HepMC3::WriterAsciiHepMC2> &hepmc,
                       const std::string &                         OUTPUTNAME) {
    OUTPUT       = OUTPUTNAME;
    FORMAT       = "hepmc2";
    outputHepMC2 = hepmc;
  }
  void SetHepMC3Output(std::shared_ptr<HepMC3::WriterAscii> &hepmc, const std::string &OUTPUTNAME) {
    OUTPUT       = OUTPUTNAME;
    FORMAT       = "hepmc3";
    outputHepMC3 = hepmc;
  }

  // Maximum weight set/get
  void   SetMaxweight(double w);
  double GetMaxweight() const;

  // Get number of CPU cores
  void SetCores(int N) {
    CORES = N;
    if (CORES == 0) {  // SETUP number of threads automatically
      CORES = std::round(std::thread::hardware_concurrency() * 1.5);

      // If autodetection fails, set 1
      if (CORES < 1) { CORES = 1; }
    }
    if (CORES < 0) {
      std::string str = "MGraniitti::SetCORES: CORES < 0";
      throw std::invalid_argument(str);
    }
  }
  int  GetCores() const { return CORES; }
  void SetIntegrator(const std::string &integrator) { INTEGRATOR = integrator; }
  void SetWeighted(bool weighted) { WEIGHTED = weighted; }
  // Output file name
  void SetOutput(const std::string &output) { OUTPUT = output; }
  // Output file format
  void SetFormat(const std::string &format) {
    if (format == "hepmc3" || format == "hepmc2" || format == "hepevt") {
      FORMAT = format;
    } else {
      throw std::invalid_argument("MGraniitti::SetFormat: Unknown output format: " + format +
                                  " (valid: hepmc3, hepmc2, hepevt)");
    }
  }

  // Special method for combining histograms from each thread
  void HistogramFusion();

  // Get cross section and error
  void GetXS(double &xs, double &xs_err) const {
    xs     = stat.sigma;
    xs_err = stat.sigma_err;
  }

  // Set number of (un)weighted events to be generated
  void SetNumberOfEvents(int n) {
    if (n < 0) {
      throw std::invalid_argument(
          "MGraniitti::SetNumberOfEvents():: Error: Number "
          "of events < 0!");
    }
    NEVENTS = n;
  }
  int GetNumberOfEvents() const { return NEVENTS; }

  // Return processes
  std::vector<std::string> GetProcessNumbers() const;

  // Initialize generator
  void Initialize();
  void Initialize(const MEikonal &eikonal_in);

  // Process object pointer, public so methods can be accessed
  MProcess *proc = nullptr;

  void SaveVGRID() const;
  void ReadVGRID(const std::string &inputfile);

  void PrintHistograms();
  void ReadGeneralParam(const json &j);
  void ReadProcessParam(const json &j);
  void ReadIntegralParam(const json &j);
  void ReadGenCuts(const json &j);
  void ReadFidCuts(const json &j);
  void ReadVetoCuts(const json &j);
  void ReadModelParam(const std::string &tune) const;

  void Generate();

  std::string FULL_PROCESS_STR = "null";  // Full process string
  std::string PROCESS          = "null";  // Physics process identifier

  bool        WEIGHTED   = false;   // Unweighted or weighted event generation
  int         NEVENTS    = 0;       // Number of events to be generated
  int         CORES      = 0;       // Number of CPU cores (threads) in use
  std::string INTEGRATOR = "null";  // Integrator (VEGAS, FLAT, ...)

  // SILENT OUTPUT
  bool HILJAA = false;

  void SetVgridFile(const std::string &inputfile) { MC_VGRID_INPUT = inputfile; }

  // Terminal input
  void ConstructTerminal(cxxopts::Options &options) const;
  void ProcessTerminal(json &j, cxxopts::ParseResult const &r) const;

 private:
  const double OUTPUT_XS_SCALE = 1E12;  // Output in picobarns


  // MC initialized from a file
  std::string MC_VGRID_INPUT = "null";

  // Integration statistics
  Stats stat;

  // Histogram fusion
  bool hist_fusion_done = false;


  // HepMC outputfile
  std::string FULL_OUTPUT_STR = "null";
  std::string OUTPUT          = "null";
  std::string FORMAT          = "null";  // hepmc3 or hepmc2 or hepevt

  std::shared_ptr<HepMC3::GenRunInfo>        runinfo      = nullptr;
  std::shared_ptr<HepMC3::WriterAscii>       outputHepMC3 = nullptr;
  std::shared_ptr<HepMC3::WriterAsciiHepMC2> outputHepMC2 = nullptr;
  std::shared_ptr<HepMC3::WriterHEPEVT>      outputHEPEVT = nullptr;

  // VEGAS creates copies here
  std::vector<MProcess *> pvec;
  MContinuum              proc_C;
  MFactorized             proc_F;
  MQuasiElastic           proc_Q;
  MParton                 proc_P;

  // Forced cross-section for HepMC3 output
  double xsforced = -1;

  // Global timing
  MTimer global_tictoc;
  MTimer local_tictoc;
  MTimer atime;
  double time_t0  = 0.0;
  double itertime = 0.0;

  // -----------------------------------------------
  // Process MC integration generic variables

  // Generation mode (integration = 0, event generation = 1)
  unsigned int GMODE = 0;

  // -----------------------------------------------
  // FLAT MC
  MCPARAM mcparam;

  // -----------------------------------------------
  // VEGAS MC

  // Parameters
  VEGASPARAM vparam;

  // DATA
  VEGASData VD;

  void VEGASInit(unsigned int init, unsigned int calls);
  int  VEGAS(unsigned int init, unsigned int calls, unsigned int iter, unsigned int N);
  void VEGASMultiThread(unsigned int N, unsigned int tid, unsigned int init,
                        unsigned int LOCALcalls);

  // -----------------------------------------------

  void UnifyHistogramBounds();

  // Calculate cross section
  void CalculateCrossSection();

  // Vegas wrapper
  double VegasWrapper(std::vector<double> &randvec, double wgt);

  // Event sampling/generation
  void CallIntegrator(unsigned int N);
  void SampleVegas(unsigned int N);
  void SampleFlat(unsigned int N);
  void SampleNeuro(unsigned int N);

  // Helper functions
  void           PrintInit() const;
  int            SaveEvent(MProcess *pr, double W, double MAXW, const gra::AuxIntData &aux);
  void           PrintStatus(unsigned int events, unsigned int N, MTimer &tictoc, double timercut);
  void           PrintStatistics(unsigned int N);
  gra::PARAM_RES ReadFactorized(const std::string &resparam_str);
  void           InitFileOutput();

  // Interpreter commands
  std::vector<aux::OneCMD> syntax;
};

// Launcher functions
void MLaunch(MGraniitti gen, int randomseed, int tid, int events);
void MThreader(MGraniitti &gen);

}  // namespace gra

#endif