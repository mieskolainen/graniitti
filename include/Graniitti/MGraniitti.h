// GRANIITTI Monte Carlo main class
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.


#ifndef MGRANIITTI_H
#define MGRANIITTI_H

// C++
#include <complex>
#include <mutex>
#include <random>
#include <thread>
#include <vector>

// HepMC
#include "HepMC/WriterHEPEVT.h"
#include "HepMC/WriterAscii.h"
#include "HepMC/WriterAsciiHepMC2.h"

// Own
#include "Graniitti/MAux.h"
#include "Graniitti/MContinuum.h"
#include "Graniitti/MEikonal.h"
#include "Graniitti/MKinematics.h"
#include "Graniitti/MMatrix.h"
#include "Graniitti/MProcess.h"
#include "Graniitti/MQuasiElastic.h"
#include "Graniitti/MFactorized.h"
#include "Graniitti/MSpin.h"
#include "Graniitti/MTimer.h"
#include "Graniitti/M4Vec.h"


namespace gra {


// Global to silent output
extern bool HILJAA;


// Simple MC parameters
struct MCPARAM {
	double PRECISION = 0.05;    // Integral relative precision
	unsigned int MIN_EVENTS  = 1000000; // Minimum number of events to be sampled
};


// Vegas MC default parameters
struct VEGASPARAM {

	unsigned int BINS = 128;            // Maximum number of bins per dimension (EVEN NUMBER!)
	double LAMBDA = 1.5;        // Regularization parameter

	// Initialization
	unsigned int NCALL = 20000;         // Number of calls per iteration
	unsigned int ITER = 15;             // Number of iterations
	double CHI2MAX = 10.0;      // Maximum chi2 in initialization
	double PRECISION = 0.01;    // Maximum relative error of cross section integral
	int DEBUG = -1;             // Debug mode

	// User cannot set these
	unsigned int MAXFDIM = 100;         // Maximum integral dimension
	double EPS = 1.0e-30;       // Epsilon parameter
};


// Vegas MC adaptation data
struct VEGASData {

	// VEGAS initialization function
	void Init(unsigned int init, const VEGASPARAM& param) {

		// First initialization: Create the grid and initial data
		if (init == 0) {
			ClearAll(param);
			InitGridDependent(param);
		}
		// Use the previous grid but NOT integral data
		else if (init == 1) {
			ClearIntegral();
			InitGridDependent(param);
		}
		// Use the previous grid and its integral data
		else if (init == 2) {
			InitGridDependent(param);	
		} else {
			throw std::invalid_argument("VEGASData::Init: Unknown init parameter = " + std::to_string(init));
		}
	}

	// Initialize 
	void InitGridDependent(const VEGASPARAM& param) {

		// Create grid spacing
		for (std::size_t j = 0; j < FDIM; ++j) {
			dxvec[j]  = region[j + FDIM] - region[j];
		}

		// If binning parameter changed from previous call
		if (param.BINS != BINS_prev) {
			for (std::size_t i = 0; i < std::max(param.BINS, BINS_prev); ++i) {
				rvec[i] = 1.0;
			}
			for (std::size_t j = 0; j < FDIM; ++j) {
				Rebin(BINS_prev / static_cast<double>(param.BINS), j, param);
			}
			BINS_prev = param.BINS;
		}
	}

	// VEGAS rebinning function (algorithm adapted from Numerical Recipes)
	void Rebin(double ac, unsigned int j, const VEGASPARAM& param) {

		unsigned int    k = 0;
		double dr = 0.0;
		double zn = 0.0;
		double zo = 0.0;

		for (std::size_t i = 0; i < param.BINS - 1; ++i) {
			while (ac > dr) {
				dr += rvec[(++k) - 1];
			}
			if (k > 1) {
				zo = xmat[k - 2][j];
			}
			zn = xmat[k - 1][j];
			dr -= ac;
			xcache[i] = zn - (zn-zo) * dr / rvec[k - 1];
		}

		for (std::size_t i = 0; i < param.BINS - 1; ++i) {
			xmat[i][j] = xcache[i];
		}

		xmat[param.BINS - 1][j] = 1.0;
	}

	// VEGAS grid optimizing function (algorithm adapted from Numerical Recipes)
	void OptimizeGrid(const VEGASPARAM& param) {

		double zo = 0;
		double zn = 0;
		double ac = 0;

		for (std::size_t j = 0; j < FDIM; ++j) {
			zo = f2mat[0][j];
			zn = f2mat[1][j];
			f2mat[0][j] = (zo + zn) / 2.0;
			dcache[j] = f2mat[0][j];

			for (std::size_t i = 2; i < param.BINS; ++i) {
				ac = zo + zn;
				zo = zn;
				zn = f2mat[i][j];
				f2mat[i - 1][j] = (ac + zn) / 3.0;
				dcache[j] += f2mat[i - 1][j];
			}
			f2mat[param.BINS - 1][j] = (zo + zn) / 2.0;
			dcache[j] += f2mat[param.BINS - 1][j];
		}

		for (std::size_t j = 0; j < FDIM; ++j) {
			ac = 0.0;
			for (std::size_t i = 0; i < param.BINS; ++i) {
				f2mat[i][j] = (f2mat[i][j] < param.EPS) ? param.EPS : f2mat[i][j];
				rvec[i] = std::pow((1.0 - f2mat[i][j] / dcache[j]) /
				         (std::log(dcache[j]) - std::log(f2mat[i][j])), param.LAMBDA);
				ac += rvec[i];
			}
			Rebin(ac / param.BINS, j, param);
		}
	}

	// Initialize sampling region [0,1] x [0,1] x ... x [0,1]
	void InitRegion(unsigned int fdim) {

		FDIM = fdim;
		region.resize(2 * FDIM, 0.0);

		// Lower bound [zero]
		for (std::size_t i = 0; i < FDIM; ++i) {
			region[i] = 0.0;
		}
		// Upper bound [one]
		for (std::size_t i = FDIM; i < 2 * FDIM; ++i) {
			region[i] = 1.0;
		}
	}

	// Clear integrated data
	void ClearIntegral() {
		sumdata = 0.0;
		sumchi2 = 0.0;
		sweight = 0.0;
	}

	// Full init
	void ClearAll(const VEGASPARAM& param) {

		// Vectors [MAXFDIM]
		dcache = std::vector<double>(param.MAXFDIM, 0.0);
		dxvec  = std::vector<double>(param.MAXFDIM, 0.0);

		// Vectors [BINS]
		rvec   = std::vector<double>(param.BINS, 0.0);
		xcache = std::vector<double>(param.BINS, 0.0);

		// Matrices [BINS x MAXFDIM]
		fmat = std::vector<std::vector<double>>(
		    param.BINS, std::vector<double>(param.MAXFDIM, 0.0));
		f2mat = std::vector<std::vector<double>>(
		    param.BINS, std::vector<double>(param.MAXFDIM, 0.0));
		xmat = std::vector<std::vector<double>>(
		    param.BINS, std::vector<double>(param.MAXFDIM, 0.0));

		// Init with 1!
		for (std::size_t j = 0; j < FDIM; ++j) {
			xmat[0][j] = 1.0;
		}

		// VEGAS integraldata
		sumdata     = 0.0;
		sumchi2     = 0.0;
		sweight     = 0.0;

		// Previous binning
		BINS_prev   = 1; 
	}

	// VEGAS scalars
	unsigned int BINS_prev   = 0;
	unsigned int FDIM        = 0;

	double fsum      = 0.0;
	double f2sum     = 0.0;

	double sumdata   = 0.0;
	double sumchi2   = 0.0;
	double sweight   = 0.0;

	// Vectors
	std::vector<double> region;

	// Vectors
	std::vector<double> dcache;
	std::vector<double> dxvec;

	// Vectors
	std::vector<double> rvec;
	std::vector<double> xcache;

	// Matrices
	std::vector<std::vector<double>> fmat;
	std::vector<std::vector<double>> f2mat;
	std::vector<std::vector<double>> xmat;
};


// Integration statistics
class Stats {

public:

	void Accumulate(const gra::AuxIntData& aux) {

		evaluations   += 1.0;

		amplitude_ok  += (aux.amplitude_ok  ? 1.0 : 0.0);
		kinematics_ok += (aux.kinematics_ok ? 1.0 : 0.0);
		fidcuts_ok    += (aux.fidcuts_ok    ? 1.0 : 0.0);
		vetocuts_ok   += (aux.vetocuts_ok   ? 1.0 : 0.0);
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

	double amplitude_ok   = 0.0;
	double kinematics_ok  = 0.0;
	double fidcuts_ok     = 0.0;
	double vetocuts_ok    = 0.0;

	// Keep as double to avoid overflow of range
	double evaluations    = 0.0; // Integrand evaluations
	double trials         = 0.0; // Event generation trials
	
	unsigned int   generated      = 0.0; // Event generation
	unsigned int   N_overflow     = 0.0; // Weight overflows

	// Cross section and its error
	double sigma       = 0.0;
	double sigma_err   = 0.0;
	double sigma_err2  = 0.0;

	// Weight statistics FLAT MC
	double Wsum        = 0.0;
	double W2sum       = 0.0;
	double maxW        = 0.0;

	// Weight statistics VEGAS MC
	double maxf        = 0.0;
	double chi2        = 0.0;
};


class MGraniitti {

public:
	// Destructor & Constructor
	MGraniitti();
	~MGraniitti();

	// For cross-section for HepMC output
	void ForceXS(double xs) {
		xsforced = xs;
	}
	
	// Read parameters
	void ReadInput(const std::string& inputfile, const std::string cmd_PROCESS = "null");

	// Initialize memory
	void InitProcessMemory(std::string process, unsigned int RNDSEED);
	void InitMultiMemory();
	
	// Set simple MC parameters
	void SetMCParam(MCPARAM& in);

	// Set VEGAS parameters
	void SetVegasParam(const VEGASPARAM& in);

	// Set external file handle
	void SetHepMC2Output(std::shared_ptr<HepMC::WriterAsciiHepMC2>& hepmc, const std::string& OUTPUTNAME) {
		OUTPUT = OUTPUTNAME;
		FORMAT = "hepmc2";
		outputHepMC2 = hepmc;
	}	
	void SetHepMC3Output(std::shared_ptr<HepMC::WriterAscii>& hepmc, const std::string& OUTPUTNAME) {
		OUTPUT = OUTPUTNAME;
		FORMAT = "hepmc3";
		outputHepMC3 = hepmc;
	}

	// Maximum weight set/get
	void   SetMaxweight(double w);
	double GetMaxweight() const;
	
	// Get number of CPU cores
	void SetCores(int N) {
		CORES = N;
		if (CORES == 0) { // SETUP number of threads automatically
			CORES = std::round(std::thread::hardware_concurrency() * 1.5);

			// If autodetection fails, set 1
			if (CORES < 1) {
				CORES = 1;
			}
		}
		if (CORES < 0) {
			std::string str = "MGraniitti::SetCORES: CORES < 0";
			throw std::invalid_argument(str);
		}
	}
	int GetCores() const {
		return CORES;
	}
	void SetIntegrator(const std::string& integrator) {
		INTEGRATOR = integrator;
	}
	void SetWeighted(bool weighted) {
		WEIGHTED = weighted;
	}
	// Output file name
	void SetOutput(const std::string& output) {
		OUTPUT = output;
	}
	// Output file format
	void SetFormat(const std::string& format) {
		if (format == "hepmc3" ||
			format == "hepmc2" ||
			format == "hepevt") {
			FORMAT = format;
		} else {
			throw std::invalid_argument("MGraniitti::SetFormat: Unknown output format: "
										+ format + " (valid: hepmc3, hepmc2, hepevt)");
		}
	}

	// Special method for combining histograms from each thread
	void HistogramFusion();

	// Get cross section and error
	void GetXS(double& xs, double& xs_err) const {
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
	int GetNumberOfEvents() const {
		return NEVENTS;
	}

	// Return processes
	std::vector<std::string> GetProcessNumbers() const;

	// Initialize generator
	void InitEikonal();
	void InitEikonalOnly();
	void Initialize();
	void Initialize(const MEikonal& eikonal_in);

	// Process object pointer, public so methods can be accessed
	MProcess* proc = nullptr;


	void PrintHistograms();
	void ReadGeneralParam(const std::string& inputfile);
	void ReadModelParam(const std::string& inputfile);
	void ReadProcessParam(const std::string& inputfile, const std::string cmd_PROCESS = "null");
	void ReadIntegralParam(const std::string& inputfile);
	void ReadGenCuts(const std::string& inputfile);
	void ReadFidCuts(const std::string& inputfile);
	void ReadVetoCuts(const std::string& inputfile);
	
	void Generate();

	std::string PROCESS = "";    // Physics process identifier
	
	bool WEIGHTED = false;       // Unweighted or weighted event generation
	int NEVENTS = 0;             // Number of events to be generated
	int CORES   = 0;             // Number of CPU cores (threads) in use
	std::string INTEGRATOR = ""; // Integrator (VEGAS, FLAT, ...)

private:

	// -------------------------------------------------------------------
	// Copy and assignment disabled by making the private
	MGraniitti(const MGraniitti&);
	MGraniitti& operator=(const MGraniitti&);
	// -------------------------------------------------------------------

	// Integration statistics
	Stats stat;

	// Histogram fusion
	bool hist_fusion_done = false;

	// HepMC outputfile
	std::string FULL_OUTPUT_STR = "";
	std::string OUTPUT = "";
	std::string FORMAT = "";   // hepmc3 or hepmc2 or hepevt
	
	// HepMC3 output
	std::shared_ptr<HepMC::WriterAscii>       outputHepMC3 = nullptr;
	std::shared_ptr<HepMC::WriterAsciiHepMC2> outputHepMC2 = nullptr;
	std::shared_ptr<HepMC::WriterHEPEVT>      outputHEPEVT = nullptr;

	// VEGAS creates copies here
	std::vector<MProcess*> pvec;
	MContinuum    proc_C;
	MFactorized   proc_F;
	MQuasiElastic proc_Q;

	// Forced cross-section for HepMC output
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

	int  Vegas(unsigned int init, unsigned int calls, unsigned int iter, unsigned int N);
	void VEGASInit(unsigned int init, unsigned int calls);
	void VEGASMultiThread(unsigned int N, unsigned int tid, unsigned int init, unsigned int LOCALcalls);
	std::vector<unsigned int> VEGASGetLocalCalls(unsigned int calls);

	// -----------------------------------------------

	// Calculate cross section
	void CalculateCrossSection();

	// Vegas wrapper
	double VegasWrapper(std::vector<double>& randvec, double wgt);

	// Event sampling/generation
	void CallIntegrator(unsigned int N);
	void SampleVegas(unsigned int N);
	void SampleFlat(unsigned int N);
	void SampleNeuro(unsigned int N);

	// Helper functions
	void PrintInit() const;
	int  SaveEvent(MProcess* pr, double W, double MAXW, const gra::AuxIntData& aux);
	void PrintStatus(unsigned int events, unsigned int N, MTimer& tictoc, double timercut);
	void PrintStatistics(unsigned int N);
	gra::PARAM_RES ReadFactorized(const std::string& resparam_str);
	void InitFileOutput();

};

// Launcher functions
void MLaunch(MGraniitti gen, int randomseed, int tid, int events);
void MThreader(MGraniitti& gen);


} // gra namespace ends

#endif