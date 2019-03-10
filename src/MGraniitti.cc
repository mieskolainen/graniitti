// GRANIITTI Monte Carlo main class
//
// (c) 2017-2019 Mikael Mieskolainen
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
#include "Graniitti/MNeuroJacobian.h"
#include "Graniitti/MAux.h"
#include "Graniitti/MContinuum.h"
#include "Graniitti/MGraniitti.h"
#include "Graniitti/MKinematics.h"
#include "Graniitti/MMath.h"
#include "Graniitti/MProcess.h"
#include "Graniitti/MQuasiElastic.h"
#include "Graniitti/MFactorized.h"
#include "Graniitti/MTimer.h"

// HepMC3
#include "HepMC/WriterHEPEVT.h"
#include "HepMC/WriterAsciiHepMC2.h"
#include "HepMC/GenEvent.h"

// Libraries
#include "json.hpp"
#include "rang.hpp"


namespace gra {

using json = nlohmann::json;

bool HILJAA = false;

using gra::aux::indices;
using gra::math::msqrt;
using gra::math::pow2;
using gra::math::pow3;
using gra::math::zi;
using gra::math::PI;


// Constructor
MGraniitti::MGraniitti() {

	// Print general layout
	PrintInit();
}

// Destructor
MGraniitti::~MGraniitti() {

	// Destroy processes
	//for (unsigned int i = 0; i < pvec.size(); ++i) {
	//	delete pvec[i];
	//}
}

void MGraniitti::PrintHistograms() {
	HistogramFusion();
	proc->PrintHistograms();
}

// Fuse VEGAS histograms for NCORES-fold statistics
void MGraniitti::HistogramFusion() {

	if (hist_fusion_done == false) {

		// pvec[0] points to base instance -> fuse histograms to that

		// START with process index 1
		for (std::size_t i = 1; i < pvec.size(); ++i) {
			// Loop over all 1D histograms
			for (auto const& xpoint : pvec[0]->h1) {
				pvec[0]->h1[xpoint.first] = pvec[0]->h1[xpoint.first] +
				                           	pvec[i]->h1[xpoint.first];
			}
			// Loop over all 2D histograms
			for (auto const& xpoint : pvec[0]->h2) {
				pvec[0]->h2[xpoint.first] = pvec[0]->h2[xpoint.first] +
				                           	pvec[i]->h2[xpoint.first];
			}
		}
		hist_fusion_done = true;

	} else {
		std::cout << "MGraniitti::HistogramFusion: Multithreaded histograms have been fused already" << std::endl;
	}
}

// Generate events
void MGraniitti::Generate() {
	if (NEVENTS > 0) { // Generate events

		// Just before event generation
		// (in order not to generate unnecessary empty files
		// if file name / output type is changed)
		InitFileOutput();

		CallIntegrator(NEVENTS);
	}
}


// This is called just before event generation
// Check against nullptr is done below, because if the output is already set externally,
// this function is not used.
void MGraniitti::InitFileOutput() {

	if (NEVENTS > 0) { 

		if (OUTPUT == "") { // OUTPUT must be set
			throw std::invalid_argument("MGraniitti::InitFileOutput: OUTPUT filename not set!");
		}

		FULL_OUTPUT_STR = gra::aux::GetBasePath(2) + "/output/" + OUTPUT + "." + FORMAT;

		if        (FORMAT == "hepmc3" && outputHepMC3 == nullptr) {
			outputHepMC3  = std::make_shared<HepMC::WriterAscii>(FULL_OUTPUT_STR);
	    } else if (FORMAT == "hepmc2" && outputHepMC2 == nullptr) {    	
			outputHepMC2  = std::make_shared<HepMC::WriterAsciiHepMC2>(FULL_OUTPUT_STR);
		} else if (FORMAT == "hepevt" && outputHEPEVT == nullptr) {
			outputHEPEVT  = std::make_shared<HepMC::WriterHEPEVT>(FULL_OUTPUT_STR);
		}
	}
}


// Return process numbers available
std::vector<std::string> MGraniitti::GetProcessNumbers() const {

	std::cout << rang::style::bold;
	std::cout << "Available processes: ";
	std::cout << rang::style::reset << std::endl << std::endl;
	std::cout << std::endl;

	std::cout << rang::style::bold << "2->3 x 1->(N-2) PS:" << rang::style::reset << std::endl;
	std::vector<std::string> list2 = proc_F.PrintProcesses();
	std::cout << std::endl;

	std::cout << rang::style::bold << "2->N PS:" << rang::style::reset << std::endl;
	std::vector<std::string> list3 = proc_C.PrintProcesses();
	std::cout << std::endl;

	std::cout << rang::style::bold << "Quasielastic:" << rang::style::reset << std::endl;
	std::vector<std::string> list1 = proc_Q.PrintProcesses();
	std::cout << std::endl;

	// Concatenate all
	list1.insert(list1.end(), list2.begin(), list2.end());
	list1.insert(list1.end(), list3.begin(), list3.end());

	return list1;
}

// (Re-)assign the pointers to the local memory space
void MGraniitti::InitProcessMemory(std::string process, unsigned int seed) {
	
	// These must be here!
	PROCESS = process;

	// <Q> process
	if (proc_Q.ProcessExist(process)) {
		proc_Q = MQuasiElastic(process);
		proc = &proc_Q;

	// <F> processes
	} else if (proc_F.ProcessExist(process)) {
		proc_F = MFactorized(process);
		proc = &proc_F;

	// <C> processes
	} else if (proc_C.ProcessExist(process)) {
		proc_C = MContinuum(process);
		proc = &proc_C;

	} else {
		std::string str =
		    "MGraniitti::InitProcessMemory: Unknown PROCESS: " + process;
		//GetProcessNumbers(); // This is done by main program
		throw std::invalid_argument(str);
	}

	// Set random seed last!
	proc->random.SetSeed(seed);
}


void MGraniitti::InitMultiMemory() {

	// ** Init multithreading memory by making copies of the process **
	if (pvec.size() != 0) {
		for (std::size_t i = 0; i < pvec.size(); ++i) {
			delete pvec[i];
		}
	}
	pvec.resize(CORES, nullptr);

	// Copy process objects for each thread
	for (int i = 0; i < CORES; ++i) {
		if        (proc_Q.ProcessExist(PROCESS)) {
			pvec[i] = new MQuasiElastic(proc_Q);
		} else if (proc_F.ProcessExist(PROCESS)) {
			pvec[i] = new MFactorized(proc_F);
		} else if (proc_C.ProcessExist(PROCESS)) {
			pvec[i] = new MContinuum(proc_C);
		}
	}

	// RANDOM SEED PER THREAD (IMPORTANT!)
	for (int i = 0; i < CORES; ++i) {
		const unsigned int tid = i + 1;
		// Deterministic seed
		const unsigned int thread_seed = static_cast<unsigned int>(
			std::max(0, (int)pvec[i]->random.GetSeed()) * tid);
		pvec[i]->random.SetSeed(thread_seed);
	}

	// SET main control pointer to the first one
	// (needed for printing etc.)
	proc = pvec[0];
	proc->PrintInit(HILJAA);
}


// Set simple MC parameters
void MGraniitti::SetMCParam(MCPARAM& in) {
	mcparam = in;
}


// Set VEGAS parameters
void MGraniitti::SetVegasParam(const VEGASPARAM& in) {
	vparam = in;
}


// Read parameters from a single JSON file
void MGraniitti::ReadInput(const std::string& inputfile, const std::string cmd_PROCESS) {

	std::cout << rang::fg::green << "MGraniitti::ReadInput: "
		+ inputfile << rang::fg::reset << std::endl << std::endl;

	ReadGeneralParam(inputfile);
	ReadProcessParam(inputfile, cmd_PROCESS);
}

// General parameter initialization
void MGraniitti::ReadGeneralParam(const std::string& inputfile) {
	// Read and parse
	const std::string data = gra::aux::GetInputData(inputfile);
	json j;
	try {
		j = json::parse(data);
	} catch (...) {
		std::string str =
		    "MGraniitti::ReadGeneralParam: Error parsing " +
		    inputfile + " (Check for extra/missing commas)";
		throw std::invalid_argument(str);
	}

	// JSON block identifier
	const std::string XID = "GENERALPARAM";

	// Setup parameters (order is important)
	SetNumberOfEvents(j[XID]["NEVENTS"]);
	SetOutput(j[XID]["OUTPUT"]);
	SetFormat(j[XID]["FORMAT"]);
	SetWeighted(j[XID]["WEIGHTED"]);
	SetIntegrator(j[XID]["INTEGRATOR"]);
	SetCores(j[XID]["CORES"]);

	// Model parameter file
	ReadModelParam(j[XID]["MODELPARAM"]);
}


// Soft model parameter initialization
void MGraniitti::ReadModelParam(const std::string& inputfile) {

	// Setup for later use
	gra::aux::MODELPARAM = inputfile;

	// Read and parse
	const std::string fullpath = gra::aux::GetBasePath(2) + "/modeldata/" + inputfile + ".json";
	const std::string data     = gra::aux::GetInputData(fullpath);

	json j;
	try {
		j = json::parse(data);
	} catch (...) {
		std::string str =
		    "MGraniitti::ReadModelParam: Error parsing " + fullpath +
		    " (Check for extra/missing commas)";
		throw std::invalid_argument(str);
	}

	// Soft model parameters
	PARAM_SOFT::DELTA_P = j["PARAM_SOFT"]["DELTA_P"];
	PARAM_SOFT::ALPHA_P = j["PARAM_SOFT"]["ALPHA_P"];
	PARAM_SOFT::gN_P    = j["PARAM_SOFT"]["gN_P"];

	double triple3P     = j["PARAM_SOFT"]["g3P"];
	PARAM_SOFT::g3P     = triple3P * PARAM_SOFT::gN_P; // Convention
	PARAM_SOFT::gamma   = j["PARAM_SOFT"]["gamma"];

	PARAM_SOFT::fc1     = j["PARAM_SOFT"]["fc1"];
	PARAM_SOFT::fc2     = j["PARAM_SOFT"]["fc2"];
	PARAM_SOFT::fc3     = j["PARAM_SOFT"]["fc3"];

	// Flat amplitude parameters
	PARAM_FLAT::b = j["PARAM_FLAT"]["b"];

	// Monopole production
	PARAM_MONOPOLE::En       = j["PARAM_MONOPOLE"]["En"];
	PARAM_MONOPOLE::Gamma0   = j["PARAM_MONOPOLE"]["Gamma0"];
	PARAM_MONOPOLE::coupling = j["PARAM_MONOPOLE"]["coupling"];

	// Regge amplitude parameters
	std::vector<double> a0   = j["PARAM_REGGE"]["a0"];
	std::vector<double> ap   = j["PARAM_REGGE"]["ap"];
	std::vector<double> sgn  = j["PARAM_REGGE"]["sgn"];
	PARAM_REGGE::a0          = a0;
	PARAM_REGGE::ap          = ap;
	PARAM_REGGE::sgn         = sgn;

	PARAM_REGGE::s0          = j["PARAM_REGGE"]["s0"];

	PARAM_REGGE::offshellFF  = j["PARAM_REGGE"]["offshellFF"];
	PARAM_REGGE::b_EXP       = j["PARAM_REGGE"]["b_EXP"];
	PARAM_REGGE::a_OREAR     = j["PARAM_REGGE"]["a_OREAR"];
	PARAM_REGGE::b_OREAR     = j["PARAM_REGGE"]["b_OREAR"];
	PARAM_REGGE::b_POW       = j["PARAM_REGGE"]["b_POW"];
	PARAM_REGGE::reggeize    = j["PARAM_REGGE"]["reggeize"];

	// Proton (Good-Walker) resonances
	std::vector<double> rc   = j["PARAM_NSTAR"]["rc"];
	PARAM_NSTAR::rc = rc;

	// Make sure they sum to one
	const double sum_rc = std::accumulate(rc.begin(), rc.end(), 0);
	for (std::size_t i = 0; i < PARAM_NSTAR::rc.size(); ++i) {
		PARAM_NSTAR::rc[i] /= sum_rc;
	}

	PARAM_SOFT::PrintParam();
}


// Process parameter initialization, Call proc->post_Constructor() after this
void MGraniitti::ReadProcessParam(const std::string& inputfile, const std::string cmd_PROCESS) {

	// Read and parse	
	const std::string data = gra::aux::GetInputData(inputfile);
	json j;
	try {
		j = json::parse(data);
	} catch (...) {
		std::string str =
		    "MGraniitti::ReadProcessParam: Error parsing " +
		    inputfile + " (Check for extra/missing commas)";
		throw std::invalid_argument(str);
	}
	const std::string XID = "PROCESSPARAM";

	// ----------------------------------------------------------------
	// Initialize process

	std::string fullstring;
	if (cmd_PROCESS == "null") {
		fullstring = j[XID]["PROCESS"];
	} else { // commandline override
		fullstring = cmd_PROCESS;
	}

	// Now separate process and decay parts by "->"
	std::string PROCESS_PART = "";
	std::string DECAY_PART   = "";
	std::size_t pos  = 0;
	std::size_t pos1 = fullstring.find("->");
	std::size_t pos2 = fullstring.find("&>"); // Isolated Phase-Space

	pos = (pos1 != std::string::npos) ? pos1 : std::string::npos; // Try to find ->
	if (pos == std::string::npos) {
		pos = (pos2 != std::string::npos) ? pos2 : std::string::npos; // Try to find &>
	}

	if (pos != std::string::npos) {
		PROCESS_PART = fullstring.substr(0, pos-1); // beginning to the pos-1
		DECAY_PART   = fullstring.substr(pos + 2);  // from pos+2 to the end
	} else {
		PROCESS_PART = fullstring; // No decay defined
	}

	// Trim extra spaces away
	gra::aux::TrimExtraSpace(PROCESS_PART);
	gra::aux::TrimExtraSpace(DECAY_PART);

	InitProcessMemory(PROCESS_PART, j[XID]["RNDSEED"]);

	// ----------------------------------------------------------------
	// SETUP process: process memory needs to be initialized before!
	
	proc->SetLHAPDF(j[XID]["LHAPDF"]);
	proc->SetScreening(j[XID]["POMLOOP"]);
	proc->SetExcitation(j[XID]["NSTARS"]);
	proc->SetHistograms(j[XID]["HIST"]);
	proc->SetFlatAmp(j[XID]["FLATAMP"]);

	if (pos2 != std::string::npos) {

		if (PROCESS_PART.find("<F>") == std::string::npos) {
			throw std::invalid_argument("MGraniitti::ReadProcessParam: Phase space isolation arrow '&>' to be used only with <F> class!");
		}
		proc->SetIsolate(true);
	}
	
	// ----------------------------------------------------------------
	// Decaymode setup
	proc->SetDecayMode(DECAY_PART);

	// Read free resonance parameters (only if needed)
	std::map<std::string, gra::PARAM_RES> RESONANCES;
	if (PROCESS_PART.find("RES") != std::string::npos) {

		const std::vector<std::string> RES = j[XID]["RES"];
		for (std::size_t i = 0; i < RES.size(); ++i) {
			const std::string str = "RES_" + RES[i] + ".json";
			RESONANCES[RES[i]] = gra::form::ReadResonance(str, proc->random.rng);
		}
		proc->SetResonances(RESONANCES);
		
		// Setup resonance branching (final step!)
		proc->SetupBranching();
	}
	// ----------------------------------------------------------------

	// Always last (we need PDG tables)
	std::vector<std::string> beam = j[XID]["BEAM"];
	std::vector<double> energy    = j[XID]["ENERGY"];
	proc->SetInitialState(beam,energy);

	// Now rest of the parameters
	ReadIntegralParam(inputfile);
	ReadGenCuts(inputfile);
	ReadFidCuts(inputfile);
	ReadVetoCuts(inputfile);
	
	// ** ALWAYS LAST AFTER READING ALL CUTS! **
	proc->post_Constructor();
}


// MC integrator parameter initialization
void MGraniitti::ReadIntegralParam(const std::string& inputfile) {

	using namespace gra::aux;

	// Read and parse
	const std::string data = gra::aux::GetInputData(inputfile);
	json j;
	try {
		j = json::parse(data);
	} catch (...) {
		std::string str =
		    "MGraniitti::ReadIntegralParam: Error parsing " +
		    inputfile + " (Check for extra/missing commas)";
		throw std::invalid_argument(str);
	}
	const std::string XID = "INTEGRALPARAM";

	// Numerical loop integration
	const int ND = j[XID]["POMLOOP"]["ND"];         AssertRange(ND,  {1, 1000}, "POMLOOP::ND",  true);
	MEikonalNumerics::SetLoopDiscretization(ND);


	// FLAT (naive) MC parameters
	MCPARAM mpam;
	mpam.PRECISION  = j[XID]["FLAT"]["PRECISION"];  AssertRange(mpam.PRECISION,  {0.0, 1.0},      "FLAT::PRECISION",  true);
	mpam.MIN_EVENTS = j[XID]["FLAT"]["MIN_EVENTS"]; AssertRange(mpam.MIN_EVENTS, {10, (unsigned int)1e9}, "FLAT::MIN_EVENTS", true);
	SetMCParam(mpam);

	try {
		j[XID]["VEGAS"]["BINS"];
	} catch (...) {
		std::cout << "MGraniitti::ReadIntegralParam: Did not found VEGAS parameter block from json input, using default" << std::endl;
		return; // Did not find VEGAS parameter block at all, use default
	}

	// VEGAS parameters
	VEGASPARAM vpam;
	vpam.BINS      = j[XID]["VEGAS"]["BINS"];      AssertRange(vpam.BINS,      {0, (unsigned int)1e9}, "VEGAS::BINS",      true);
	if ((vpam.BINS % 2) != 0) { throw std::invalid_argument("VEGAS::BINS = " + std::to_string(vpam.BINS) + " should be even number"); }

	vpam.LAMBDA    = j[XID]["VEGAS"]["LAMBDA"];    AssertRange(vpam.LAMBDA,    {1.0, 10.0},    "VEGAS::LAMBDA",    true);
	vpam.NCALL     = j[XID]["VEGAS"]["NCALL"];     AssertRange(vpam.NCALL,     {50,(unsigned int)1e9}, "VEGAS::NCALL",     true);
	vpam.ITER      = j[XID]["VEGAS"]["ITER"];      AssertRange(vpam.ITER,      {1, (unsigned int)1e9}, "VEGAS::ITER",      true);
	vpam.CHI2MAX   = j[XID]["VEGAS"]["CHI2MAX"];   AssertRange(vpam.CHI2MAX,   {0.0,  1e3},    "VEGAS::CHI2MAX",   true);
	vpam.PRECISION = j[XID]["VEGAS"]["PRECISION"]; AssertRange(vpam.PRECISION, {0.0,  1.0},    "VEGAS::PRECISION", true);
	vpam.DEBUG     = j[XID]["VEGAS"]["DEBUG"];     AssertSet(vpam.DEBUG,       {-1, 0, 1},     "VEGAS::DEBUG",     true);

	SetVegasParam(vpam);
}


// Generator cuts
void MGraniitti::ReadGenCuts(const std::string& inputfile) {
	// Read and parse
	const std::string data = gra::aux::GetInputData(inputfile);
	json j;
	try {
		j = json::parse(data);
	} catch (...) {
		std::string str = "MGraniitti::ReadGenCuts: Error parsing " + inputfile + " (Check for extra/missing commas)";
		throw std::invalid_argument(str);
	}
	const std::string XID = "GENCUTS";

	gra::GENCUT gcuts;

	// Continuum phase space class
	if (PROCESS.find("<C>") != std::string::npos) {

		// Daughter rapidity
		std::vector<double> rap;
		try {
		std::vector<double> temp = j[XID]["<C>"]["Rap"]; rap = temp;
		} catch (...) {
		throw std::invalid_argument("MGraniitti::ReadGenCuts: <C> phase space class requires from user: \"GENCUTS\" : { \"<C>\" : { \"Rap\" : [min, max] }}");
		}
		gcuts.rap_min = rap[0];
		gcuts.rap_max = rap[1];
		gra::aux::AssertCut(rap, "GENCUTS::<C>::Rap", true);

		// This is optional, intermediate kt
		std::vector<double> kt;
		try {
		std::vector<double> temp = j[XID]["<C>"]["Kt"]; kt = temp;
		} catch (...) {
			// Do nothing
		}
		if (kt.size() != 0) {
			gcuts.kt_min   = kt[0];
			gcuts.kt_max   = kt[1];
			gra::aux::AssertCut(kt, "GENCUTS::<C>::Kt", true);
		}

		// This is optional, forward leg pt
		std::vector<double> pt;
		try {
		std::vector<double> temp = j[XID]["<C>"]["Pt"]; pt = temp;
		} catch (...) {
			// Do nothing
		}
		if (pt.size() != 0) {
			gcuts.forward_pt_min   = pt[0];
			gcuts.forward_pt_max   = pt[1];
			gra::aux::AssertCut(pt, "GENCUTS::<C>::Pt", true);
		}
	}

	// Factorized phase space class
	if (PROCESS.find("<F>") != std::string::npos) {

		// System rapidity
		std::vector<double> Y;
		try {
		std::vector<double> temp = j[XID]["<F>"]["Rap"]; Y = temp;
		} catch (...) {
		throw std::invalid_argument("MGraniitti::ReadGenCuts: <F> phase space class requires from user: \"GENCUTS\" : { \"<F>\" : { \"Rap\" : [min, max] }}");
		}
		gcuts.Y_min = Y[0];
		gcuts.Y_max = Y[1];
		gra::aux::AssertCut(Y, "GENCUTS::<F>::Rap", true);

		// System mass
		std::vector<double> M;
		try {
		std::vector<double> temp = j[XID]["<F>"]["M"]; M = temp;
		} catch (...) {
		throw std::invalid_argument("MGraniitti::ReadGenCuts: <F> phase space class requires from user: \"GENCUTS\" : { \"<F>\" : { \"M\" : [min, max] }}");
		}
		gcuts.M_min = M[0];
		gcuts.M_max = M[1];
		gra::aux::AssertCut(M, "GENCUTS::<F>::M", true);
		
		// This is optional, forward leg pt
		std::vector<double> pt;
		try {
		std::vector<double> temp = j[XID]["<F>"]["Pt"]; pt = temp;
		} catch (...) {
			// Do nothing
		}
		if (pt.size() != 0) {
			gcuts.forward_pt_min = pt[0];
			gcuts.forward_pt_max = pt[1];
			gra::aux::AssertCut(pt, "GENCUTS::<F>::Pt", true);
		}
	}
	
	// Quasielastic phase space class
	if (PROCESS.find("<Q>") != std::string::npos) {	

		std::vector<double> Xi;
		try {
		std::vector<double> temp = j[XID]["<Q>"]["Xi"]; Xi = temp;
		} catch (...) {
		throw std::invalid_argument("MGraniitti::ReadGenCuts: <Q> phase space class requires from user: \"GENCUTS\" : { \"<Q>\" : { \"Xi\" : [min, max] }}");	
		}
		gcuts.XI_min = Xi[0];
		gcuts.XI_max = Xi[1];
		gra::aux::AssertCut(Xi, "GENCUTS::<Q>::Xi", true);
	}

	proc->SetGenCuts(gcuts);
}

// Fiducial cuts
void MGraniitti::ReadFidCuts(const std::string& inputfile) {
	// Read and parse
	const std::string data = gra::aux::GetInputData(inputfile);
	json j;
	try {
		j = json::parse(data);
	} catch (...) {
		std::string str = "MGraniitti::ReadFidCuts: Error parsing " +
		                  inputfile +
		                  " (Check for extra/missing commas)";
		throw std::invalid_argument(str);
	}

	// Fiducial cuts
	gra::FIDCUT fcuts;
	const std::string XID = "FIDCUTS";

	try {
		fcuts.active = j[XID]["active"];
	} catch (...) {
		return; // Did not find cut block at all
	}

	// Particles
	{
		std::vector<double> eta = j[XID]["PARTICLE"]["Eta"];
		fcuts.eta_min = eta[0];
		fcuts.eta_max = eta[1];
		gra::aux::AssertCut(eta, "FIDCUTS::PARTICLE::Eta", true);

		std::vector<double> rap = j[XID]["PARTICLE"]["Rap"];
		fcuts.rap_min = rap[0];
		fcuts.rap_max = rap[1];
		gra::aux::AssertCut(rap, "FIDCUTS::PARTICLE::Rap", true);

		std::vector<double> pt = j[XID]["PARTICLE"]["Pt"];
		fcuts.pt_min = pt[0];
		fcuts.pt_max = pt[1];
		gra::aux::AssertCut(pt,  "FIDCUTS::PARTICLE::Pt", true);

		std::vector<double> Et = j[XID]["PARTICLE"]["Et"];
		fcuts.Et_min = Et[0];
		fcuts.Et_max = Et[1];
		gra::aux::AssertCut(Et,  "FIDCUTS::PARTICLE::Et", true);
	}

	// System
	{
		std::vector<double> M = j[XID]["SYSTEM"]["M"];
		fcuts.M_min = M[0];
		fcuts.M_max = M[1];
		gra::aux::AssertCut(M, "FIDCUTS::SYSTEM::M", true);

		std::vector<double> Y = j[XID]["SYSTEM"]["Rap"];
		fcuts.Y_min = Y[0];
		fcuts.Y_max = Y[1];
		gra::aux::AssertCut(Y, "FIDCUTS::SYSTEM::Rap", true);

		std::vector<double> Pt = j[XID]["SYSTEM"]["Pt"];
		fcuts.Pt_min = Pt[0];
		fcuts.Pt_max = Pt[1];
		gra::aux::AssertCut(Pt, "FIDCUTS::SYSTEM::Pt", true);
	}

	// Forward
	{
		std::vector<double> t = j[XID]["FORWARD"]["t"];
		fcuts.forward_t_min = t[0];
		fcuts.forward_t_max = t[1];
		gra::aux::AssertCut(t, "FIDCUTS::FORWARD::t", true);

		std::vector<double> M = j[XID]["FORWARD"]["M"];
		fcuts.forward_M_min = M[0];
		fcuts.forward_M_max = M[1];
		gra::aux::AssertCut(M, "FIDCUTS::FORWARD::M", true);
	}

	// Set fiducial cuts
	proc->SetFidCuts(fcuts);

	// Set user cuts
	proc->SetUserCuts(j[XID]["USERCUTS"]);
}


// Fiducial cuts
void MGraniitti::ReadVetoCuts(const std::string& inputfile) {

	// Read and parse
	const std::string data = gra::aux::GetInputData(inputfile);
	json j;
	try {
		j = json::parse(data);
	} catch (...) {
		std::string str = "MGraniitti::ReadVetoCuts: Error parsing " +
		                  inputfile + " (Check for extra/missing commas)";
		throw std::invalid_argument(str);
	}

	// Fiducial cuts
	gra::VETOCUT veto;
	const std::string XID = "VETOCUTS";

	try {
		veto.active = j[XID]["active"];
	} catch (...) {
		return; // Did not find cut block at all
	}

	// Find domains
	for (std::size_t i = 0; i < 100; ++i) {
		const std::string NUMBER = std::to_string(i);

		gra::VETODOMAIN domain;
		std::vector<double> eta;
		std::vector<double> pt;
		
		// try to find new domain
		try {
			std::vector<double> temp1 = j[XID][NUMBER]["Eta"]; eta = temp1;
			std::vector<double> temp2 = j[XID][NUMBER]["Pt"];  pt  = temp2;
		} catch (...) {
			break; // Not found
		}

		domain.eta_min = eta[0];
		domain.eta_max = eta[1];
		gra::aux::AssertCut(eta, "VETOCUTS::" + NUMBER + "::Eta", true);

		domain.pt_min  = pt[0];
		domain.pt_max  = pt[1];
		gra::aux::AssertCut(pt,  "VETOCUTS::" + NUMBER + "::Pt", true);

		veto.cuts.push_back(domain);
	}

	// Set fiducial cuts
	proc->SetVetoCuts(veto);
}

// Get maximum vegasweight
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
		    "MGraniitti::SetMaxweight: Maximum weight: " +
		    std::to_string(weight) + " not valid!";
		throw std::invalid_argument(str);
	}
}

void MGraniitti::PrintInit() const {

	if (!HILJAA) {

		gra::aux::PrintFlashScreen(rang::fg::yellow);
		std::cout << rang::style::bold
		          << "GRANIITTI - Monte Carlo event generator for "
		             "high energy diffraction"
		          << rang::style::reset << std::endl
		          << std::endl;
		gra::aux::PrintVersion();
		gra::aux::PrintBar("-");
		
 		const double GB  = pow3(1024.0);
		printf("Running on %s (%d CORE / %0.2f GB RAM) at %s \n",
			gra::aux::HostName().c_str(), std::thread::hardware_concurrency(), 
			gra::aux::TotalSystemMemory() / GB, gra::aux::DateTime().c_str());
		int64_t size = 0; int64_t free = 0; int64_t used = 0;
		gra::aux::GetDiskUsage("/", size, free, used);
		printf("Path '/': size | used | free = %0.1f | %0.1f | %0.1f GB \n", size/GB, used/GB, free/GB);
		std::cout << "Program path: " << gra::aux::GetBasePath(2) << std::endl;
		std::cout << gra::aux::SystemName() << std::endl;
		gra::aux::PrintBar("-");
		
	}
}

void MGraniitti::InitEikonal() {
	// Print general and process init
	// PrintInit();
	const bool onlyeikonal = false;
	proc->Eikonal.S3Constructor(proc->GetMandelstam_s(),
		proc->GetInitialState(), onlyeikonal);
}

void MGraniitti::InitEikonalOnly() {
	// Print general and process init
	// PrintInit();
	const bool onlyeikonal = true;
	proc->Eikonal.S3Constructor(proc->GetMandelstam_s(),
		proc->GetInitialState(), onlyeikonal);
}

void MGraniitti::Initialize() {
	// Print out basic information
	std::cout << std::endl;
	std::cout << rang::style::bold
	          << "General setup:" << rang::style::reset << std::endl
	          << std::endl;
	std::cout << "Output file:            " << OUTPUT     << std::endl;
	std::cout << "Output format:          " << FORMAT     << std::endl;
	std::cout << "Multithreading:         " << CORES      << std::endl;
	std::cout << "Integrator:             " << INTEGRATOR << std::endl;
	std::cout << "Number of events:       " << NEVENTS    << std::endl;
	
	std::string str = (WEIGHTED == true) ? "weighted" : "unweighted";
	std::cout << rang::fg::green << "Generation mode:        " << str
		<< rang::fg::reset << std::endl;
	std::cout << std::endl;
	
	// If Eikonals are not yet initialized and pomeron screening loop is on
	if (proc->Eikonal.IsInitialized() == false &&
	    proc->GetScreening() == true) {
		proc->Eikonal.S3Constructor(proc->GetMandelstam_s(),
									proc->GetInitialState(), false);
	}
	CallIntegrator(0);
}


// Initialize with external Eikonal
void MGraniitti::Initialize(const MEikonal& eikonal_in) {

	proc->SetEikonal(eikonal_in);
	CallIntegrator(0);
}


void MGraniitti::CallIntegrator(unsigned int N) {

	// Initialize global clock
	if (N == 0) { global_tictoc = MTimer(true); }

	// Sample the phase space
	if      (INTEGRATOR == "VEGAS") {
		SampleVegas(N);
	}
	else if (INTEGRATOR == "FLAT") {
	 	SampleFlat(N);
	}
	else if (INTEGRATOR == "NEURO") {
	 	SampleNeuro(N);
	}
	else {
		throw std::invalid_argument("MGraniitti::CallIntegrator: Unknown INTEGRATOR = " 
			+ INTEGRATOR + " (use VEGAS, FLAT, NEURO)");
	}
}


// Initialize and generate events using VEGAS MC
void MGraniitti::SampleVegas(unsigned int N) {

	if (N == 0) {
		InitMultiMemory();
		GMODE = 0; // Pure integration
	}
	if (N > 0) {
		GMODE = 1; // Event generation
	}

	// ******************************************************
	// Create integral sampling region array

	VD.InitRegion(proc->GetdLIPSDim());

	// ******************************************************

	// Pure integration mode
	if (GMODE == 0) {

		unsigned int BURNIN_ITER = 3; // BURN-IN iterations (default)!

		// Initialize GRID
		do {
			unsigned int init = 0;
			int factor = 0;

			do { // Loop until stable
				factor = Vegas(init, vparam.NCALL, BURNIN_ITER, N);
				vparam.NCALL = vparam.NCALL * factor;
				if (factor == 1) { break; } // We are done
			} while (true);

			// Increase CALLS and re-run, if too fast
			const double MINTIME = 0.3; // seconds
			if (itertime < MINTIME) {
				const double time_per_iter = itertime / vparam.NCALL;

				// Max, because this is only minimum condition
				vparam.NCALL = std::max((unsigned int) vparam.NCALL, (unsigned int)(MINTIME / time_per_iter));
				Vegas(init, vparam.NCALL, BURNIN_ITER, N);				
			}

			// Now re-calculate by skipping burn-in (init = 0) iterations
			// because they detoriate the integral value by bad grid
			init = 1;
			factor = Vegas(init, vparam.NCALL, vparam.ITER, N);

			if (factor == 1) {
				break;
			} else {
				BURNIN_ITER = 2*BURNIN_ITER;
			}
		} while (true);
	}

	// Event generation mode
	if (GMODE == 1) {
		const unsigned int init = 2;
		const unsigned int itermin = 1E9;
		Vegas(init, vparam.NCALL * 10, itermin, N);
	}
}


// Create number of calls per thread, they need to sum to calls
std::vector<unsigned int> MGraniitti::VEGASGetLocalCalls(unsigned int calls) {

	std::vector<unsigned int> LOCALcalls(CORES, 0.0);
	unsigned int sum = 0;
	for (int k = 0; k < CORES; ++k) {
		LOCALcalls[k] = std::floor(calls / CORES);
		sum += LOCALcalls[k];
	}
	// Add remainder to the thread number[0]
	LOCALcalls[0] += calls - sum;
	return LOCALcalls;
}


// Multithreaded VEGAS integrator code
// [close to optimal importance sampling iff
//  integrand factorizes dimension by dimension]
//
// Original algorithm from:
// [REFERENCE: Lepage, G.P. Journal of Computational Physics, 1978]
// https://en.wikipedia.org/wiki/VEGAS_algorithm

int MGraniitti::Vegas(unsigned int init, unsigned int calls, unsigned int itermin, unsigned int N) {

	// First the initialization
	VD.Init(init, vparam);

	if (init == 0 && !HILJAA) {
		gra::aux::ClearProgress();
		std::cout << rang::style::bold;
		printf("VEGAS burn-in iterations: \n\n");
		std::cout << rang::style::reset;
	}

	if (init == 1 && !HILJAA) {
		gra::aux::ClearProgress(); // Progressbar clear
		gra::aux::PrintBar("-");
		std::cout << rang::style::bold;
		std::cout << "VEGAS importance sampling:" << std::endl << std::endl;
		std::cout << rang::style::reset;
		printf("<Multithreading CORES = %d> \n\n", CORES);

		printf("- dimension = %u \n",    VD.FDIM);
		printf(" \n");
		printf("- itermin   = %d \n",    itermin);
		printf("- calls     = %d \n",    calls);
		printf("- BINS      = %u \n",    vparam.BINS);
		printf("- PRECISION = %0.4f \n", vparam.PRECISION);
		printf("- CHI2MAX   = %0.4f \n", vparam.CHI2MAX);
		printf("- LAMBDA    = %0.4f \n", vparam.LAMBDA);
		printf("- DEBUG     = %d \n",    vparam.DEBUG);
		gra::aux::PrintBar("-");
	}

	// Reset local timers
	local_tictoc = MTimer(true);
	atime        = MTimer(true); // For progressbar

	// -------------------------------------------------------------------
	// Create number of calls per thread, their sum == calls
	std::vector<unsigned int> LOCALcalls = VEGASGetLocalCalls(calls);

	MTimer gridtic;

	// VEGAS grid iterations
	for (std::size_t iter = 0; iter < itermin; ++iter) {

		if (init == 0 && iter == 2) { // Save time for one iteration
			itertime = gridtic.ElapsedSec() / 3; // Average over 3 iter
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
			threads.push_back(std::thread([=] {
				VEGASMultiThread(N, tid, init, LOCALcalls[tid]);
			}));
		}

		// Loop again to join the threads
		for (auto& t : threads) {
			t.join();
		}
		if (gra::aux::globalExceptionPtr) { // Exception handling of threads
			std::rethrow_exception(gra::aux::globalExceptionPtr);
		}

		// --------------------------------------------------------------
		// Estimates based on this iteration

		VD.fsum  *= 1.0 / static_cast<double>(calls);
		VD.f2sum *= pow2(1.0 / static_cast<double>(calls));

		// Local integral estimate
		const double I_this  = VD.fsum;

		VD.f2sum = msqrt(VD.f2sum * calls);
		VD.f2sum = (VD.f2sum - VD.fsum) * (VD.f2sum + VD.fsum); // - +
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
		stat.sigma       = VD.sumdata / VD.sweight;
		stat.sigma_err   = msqrt(1.0 / VD.sweight);
		stat.sigma_err2  = pow2(stat.sigma_err);

		// Chi2
		double chi2this = (VD.sumchi2 - pow2(VD.sumdata)/VD.sweight) / (iter + 1E-5);
		chi2this        = (chi2this < 0) ? 0.0 : chi2this;
		stat.chi2       = chi2this;
		// --------------------------------------------------------------

		// Status
		if (GMODE == 0) {
			PrintStatus(stat.evaluations, N, local_tictoc, -1.0);

			if (atime.ElapsedSec() > 0.1) {
				gra::aux::PrintProgress((iter + 1) / static_cast<double>(itermin));
				atime.Reset();
			}

			// ==========================================================
			// **** Update VEGAS grid (not during event generation) ****
			VD.OptimizeGrid(vparam);
			// ==========================================================
		}

		// Got enough events generated
		if (GMODE == 1 && stat.generated >= N) {
			goto stop;
		}

		// Fatal error and convergence restart treatment
		if (GMODE == 0 && iter > 0) {

			if (std::isnan(stat.sigma) || std::isinf(stat.sigma)) {
				gra::aux::ClearProgress();
				throw std::invalid_argument("VEGAS:: Integral inf/nan: FATAL ERROR, code or parameters needs fixing.");
			}
			if (stat.sigma < 1e-60) {
				gra::aux::ClearProgress();
				throw std::invalid_argument("VEGAS:: Integral < 1E-60: Check VEGAS parameters, decaymode sanity, generation and fiducial cuts!");
			}
			if (stat.chi2 > 50) {
				gra::aux::ClearProgress();
				printf("VEGAS:: chi2 = %0.1f > 50: Convergence problem, increasing 10 x calls \n", stat.chi2);
				return 10; // 10 times more calls
			}
		}

		if (vparam.DEBUG >= 0) {
			gra::aux::ClearProgress();
			printf("VEGAS:: local iter = %4lu integral = %0.5E +- std = "
			       "%0.5E \t [global integral = %0.5E +- std = %0.5E] \t "
			       "chi2this = %0.2f \n", iter + 1, I_this, gra::math::msqrt(I2_this), stat.sigma, stat.sigma_err, chi2this);

			if (vparam.DEBUG > 0) {
				for (std::size_t j = 0; j < VD.FDIM; ++j) {
					printf("VEGAS:: data for dimension j = %lu (VD.FDIM = %u) \n", j, VD.FDIM);
					
					for (std::size_t i = 0; i < vparam.BINS; ++i) {
						printf("xmat[%3lu][j] = %0.5E, fmat[%3lu][j] = %0.5E \n", i, VD.xmat[i][j], i, VD.fmat[i][j]);
					}
				}
			}
		}

		// We are not below the chi2 or precision condition -> increase iterations
		if (init > 0) { // do not consider in burn-in (init = 0) mode
			if ((iter == (itermin - 1) && stat.chi2 > vparam.CHI2MAX) ||
			    (iter == (itermin - 1) && stat.sigma_err / stat.sigma > vparam.PRECISION)) {
				++itermin;
			}
		}

	} // Main grid iteration loop

stop: // We jump here once finished (GOTO point)

	// Final statistics
	if (init > 0) { PrintStatistics(N); }

	return 1; // Return 1 for good
}

// This is called once for every VEGAS grid iteration
void MGraniitti::VEGASMultiThread(unsigned int N, unsigned int THREAD_ID, unsigned int init, unsigned int LOCALcalls) {

	double zn = 0.0;
	double zo = 0.0;
	double ac = 0.0;

	try {

		for (std::size_t k = 0; k < LOCALcalls; ++k) { // ** LOCALcalls ~= calls / CORES **

			double vegasweight = 1.0;

			// Phase space point vector
			std::vector<double> xpoint(pvec[THREAD_ID]->GetdLIPSDim(), 0.0);

			// Loop over dimensions and construct random vector xpoint
			std::vector<double> indvec(VD.FDIM, 0);
			for (std::size_t j = 0; j < VD.FDIM; ++j) {

				// Draw random number
				zn        = pvec[THREAD_ID]->random.U(0,1) * vparam.BINS + 1.0;
				indvec[j] = std::max((unsigned int)1, std::min((unsigned int)zn, (unsigned int)vparam.BINS));

				if (indvec[j] > 1) {
					zo = VD.xmat[indvec[j] - 1][j] - VD.xmat[indvec[j] - 2][j];
					ac = VD.xmat[indvec[j] - 2][j] + (zn - indvec[j]) * zo;
				} else {
					zo = VD.xmat[indvec[j] - 1][j];
					ac = (zn - indvec[j]) * zo;
				}
				
				// Multidim space vector component
				xpoint[j]    = VD.region[j] + ac * VD.dxvec[j];
				vegasweight *= zo * vparam.BINS;

			} // VD.FDIM loop

			// *******************************************************************
			// ****** Call the process under integration to get the weight *******

			gra::AuxIntData aux;
			aux.vegasweight = vegasweight;
			aux.burn_in_mode = (init == 0) ? true : false;

			const double W = pvec[THREAD_ID]->EventWeight(xpoint, aux);

			// *******************************************************************

			// @@ Multithreading lock (FAST SECTION) @@
			gra::aux::g_mutex.lock();

			// Increase statistics
			stat.Accumulate(aux);

			// *** Importance weighting ***
			const double f  = W * vegasweight;
			const double f2 = pow2(f);
			
			VD.fsum  += f;
			VD.f2sum += f2;

			// Loop over dimensions and add importance weighted results
			for (std::size_t j = 0; j < VD.FDIM; ++j) {
				VD.fmat[indvec[j] - 1][j]  += f;
				VD.f2mat[indvec[j] - 1][j] += f2;
			}

			// ----------------------------------------------------------
			// Initialization (integration) mode
			if (GMODE == 0) {
				// Do not consider burn-in phase weights (unstable)
				if (init != 0) {
					// Maximum raw weight (for general information, not used here)
					if (W > stat.maxW) { stat.maxW = W; }
					// Maximum total VEGAS importance weighted
					if (f > stat.maxf) { stat.maxf = f; }
				}
			}
			gra::aux::g_mutex.unlock();
			// @@ Multithreading unlock (FAST SECTION) @@

			// ----------------------------------------------------------
			// Event generation mode
			if (GMODE == 1) {

				// Enough events
				if (stat.generated == (unsigned int)GetNumberOfEvents()) { break; }
				
				// Event trial
				SaveEvent(pvec[THREAD_ID], f, stat.maxf, aux);

				if (THREAD_ID == 0 &&
				    atime.ElapsedSec() > 0.5) {
					PrintStatus(stat.generated, N, local_tictoc, 10.0);
					gra::aux::PrintProgress(stat.generated / static_cast<double>(N));
					atime.Reset();
				}
			}
		} // calls loop

	} catch (...) {
		// Set the global exception pointer if exception arises
		// This is because of multithreading
		gra::aux::globalExceptionPtr = std::current_exception();
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
	if (N > 0) {
		GMODE = 1;
	}

	// Get dimension of the phase space
	const unsigned int dim = proc->GetdLIPSDim();
	std::vector<double> randvec(dim, 0.0);

	// Reset local timer
	local_tictoc = MTimer(true);

	// Progressbar
	atime = MTimer(true);

	// Reservation (test)
	//using namespace std::placeholders; // _1, _2, ... come from here
 	//std::function<double(const std::vector<double>&, int, double, std::vector<double>&, bool&)> fifth = 
 	//	std::bind(&MQuasiElastic::EventWeight, &proc_Q, _1, _2, _3, _4, _5);
 	
	// Event loop
	while (true) {

		// ** Generate new random numbers **
		for (const auto& i : indices(randvec)) {
			randvec[i] = proc->random.U(0,1);
		}

		// Generate event
		// Used for in-out control of the process
		gra::AuxIntData aux;
		aux.vegasweight  = 1.0;
		aux.burn_in_mode = false;
		const double W = proc->EventWeight(randvec, aux);

		// Increase statistics
		stat.Accumulate(aux);
		stat.Wsum  += W;
		stat.W2sum += pow2(W);
		
		// Update cross section estimate
		stat.CalculateCrossSection();
		
		// Initialization
		if (GMODE == 0) {
			stat.maxW = (W > stat.maxW) ? W : stat.maxW; // Update maximum weight
			PrintStatus(stat.evaluations, N, local_tictoc, 10.0);

			if ((stat.sigma_err / stat.sigma) < mcparam.PRECISION &&
			    stat.evaluations > mcparam.MIN_EVENTS) { break; }

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
	if (N > 0) {
		GMODE = 1;
	}

	// Get dimension of the phase space
	const unsigned int D = proc->GetdLIPSDim();

	// Reset timers
	local_tictoc = MTimer(true);
	atime        = MTimer(true);
	
	
	// Reservation (test)
	// using namespace std::placeholders; // _1, _2, ... come from here
 	// std::function<double(const std::vector<double>&, int, double, std::vector<double>&, bool&)> fifth = 
 	// std::bind(&MQuasiElastic::EventWeight, &proc_Q, _1, _2, _3, _4, _5);

	if ( N == 0) {

		// -------------------------------------------------------
		// Bind the object and call it
		gra::neurojac::procptr = proc; // First set address
		// -------------------------------------------------------

		gra::neurojac::MNeuroJacobian neurojac;
		gra::neurojac::BATCHSIZE = 100;

	    // Set network layer dimensions [first, ..., output]
		gra::neurojac::par.D = D; // Integrand dimension

		gra::neurojac::par.L.push_back(gra::neurojac::Layer(2,D)); // Input
		gra::neurojac::par.L.push_back(gra::neurojac::Layer(2,2));
		gra::neurojac::par.L.push_back(gra::neurojac::Layer(2,2));
		gra::neurojac::par.L.push_back(gra::neurojac::Layer(D,2)); // Output

		neurojac.Optimize();
	}

	// -------------------------------------------------------------------
	// Now to the event generation


	// Lambda capture
    std::vector<double> u(D);

    auto NeuroSample = [&] (gra::AuxIntData& aux)  {

		// Prior p(z) distribution sampling
		VectorXdual z(u.size());
		for (std::size_t i = 0; i < D; ++i) {
			z[i] = proc->random.G(0,1);
		}
	    const double p = val(gra::neurojac::gaussprob(z,0,1));
	    
		// Evaluate network map
	    VectorXdual u_ = gra::neurojac::G_net(z);

	    // Evaluate the Jacobian matrix du/dz
	    MatrixXd dudz = jacobian(gra::neurojac::G_net, u_, z);
	    
	    // Abs Jacobian determinant and inverse prior
	    const double jacweight = abs(dudz.determinant()) / p;
	    
	    for (std::size_t i = 0; i < D; ++i) {
	    	u[i] = val(u_[i]);
	    }

	    // Evaluate event weight
		aux.vegasweight  = jacweight;
		aux.burn_in_mode = false;
	    const double weight = proc->EventWeight(u, aux) * jacweight;

	    return weight;
    };


	// Event loop
	while (true) {

		// aux used for in-out control of the process
		gra::AuxIntData aux;
		const double W = NeuroSample(aux);

		// Increase statistics
		stat.Accumulate(aux);
		stat.Wsum  += W;
		stat.W2sum += pow2(W);
		
		// Update cross section estimate
		stat.CalculateCrossSection();
		
		// Initialization
		if (GMODE == 0) {
			stat.maxW = (W > stat.maxW) ? W : stat.maxW; // Update maximum weight
			PrintStatus(stat.evaluations, N, local_tictoc, 10.0);

			if ((stat.sigma_err / stat.sigma) < mcparam.PRECISION &&
			    stat.evaluations > mcparam.MIN_EVENTS) { break; }

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
int MGraniitti::SaveEvent(MProcess* pr, double weight, double MAXWEIGHT, const gra::AuxIntData& aux) {

	gra::aux::g_mutex.lock();
	stat.trials += 1; // This is one trial more
	
	// 0. Hit-Miss
	bool hit_in = proc->random.U(0,1) < (weight / MAXWEIGHT);

	// 2. We see weight larger than maxweight
	if (!WEIGHTED && weight > MAXWEIGHT) {
		++stat.N_overflow;
		hit_in = false;
	}
	gra::aux::g_mutex.unlock();
	
	// Three ways to accept the event. N.B. weight > 0 needed if the amplitude fails numerically
	if ((hit_in && aux.Valid()) || (WEIGHTED && aux.Valid()) || aux.forced_accept) {
		
		// Create HepMC event (do not lock yet for speed)
		HepMC::GenEvent evt(HepMC::Units::GEV, HepMC::Units::MM);
		if (!pr->EventRecord(evt)) { // Event not ok!
			//std::cout << "MGraniitti::SaveEvent: Last moment rare veto!" << std::endl;
			return 2;
		}
		// Event ok, continue >>
		
		// @@@ THIS IS THREAD-NON-SAFE -> LOCK IT
		gra::aux::g_mutex.lock();   
		
		// ** This is a multithreading race-condition treatment **
		if (stat.generated == (unsigned int) GetNumberOfEvents()) {
			stat.trials -= 1; // Correct statistics, unnecessary trial
			gra::aux::g_mutex.unlock();
			return 1;
		}

		// Save cross section information (HepMC format wants it event
		// by event)
		std::shared_ptr<HepMC::GenCrossSection> xsobj = std::make_shared<HepMC::GenCrossSection>();
		evt.add_attribute("GenCrossSection", xsobj);
		
		// Save event weight (unweighted events with weight 1)
		const double HepMC_weight = WEIGHTED ? weight : 1.0;
		evt.weights()[0] = HepMC_weight; // add more weights with .push_back()

		// Now add the value in picobarns [HepMC convention]
		if (xsforced > 0) {
			xsobj->set_cross_section(xsforced*1E12, 0); // external fixed one
		} else {
			xsobj->set_cross_section(stat.sigma*1E12, stat.sigma_err*1E12);
		}
		if      (FORMAT == "hepmc3") {
			outputHepMC3->write_event(evt);
		}
		else if (FORMAT == "hepmc2") {
			outputHepMC2->write_event(evt);
		}
		else if (FORMAT == "hepevt") {
			outputHEPEVT->write_event(evt);
		}
		else {
			throw std::invalid_argument("MGraniitti::SaveEvent: Unknown output FORMAT " + FORMAT);
		}

		// LAST STEP
		stat.generated += 1; // +1 event generated

		gra::aux::g_mutex.unlock();
		// @@ THIS IS THREAD-NON-SAFE <- LOCK IT @@
		return 0;
	} else {
		return 1;
	}
}


// Intermediate statistics
void MGraniitti::PrintStatus(unsigned int events, unsigned int N, MTimer& tictoc, double timercut) {

	if (tictoc.ElapsedSec() > timercut) {
		tictoc.Reset();
		gra::aux::ClearProgress();

		double peak_use     = 0.0;
		double resident_use = 0.0;
		const double MB     = 1024*1024;
		const double GB     = MB * 1024;
		gra::aux::GetProcessMemory(peak_use, resident_use);
		peak_use     /= MB;
		resident_use /= MB;

		if (GMODE == 0) {
			const double global_lap = global_tictoc.ElapsedSec();
			printf(
			    "[%0.1f MB] xs: %9.3E, er: %7.5f [chi2 = %2.1f], %4.1f "
			    "min ~ %0.1E Hz \n", resident_use, stat.sigma, stat.sigma_err / stat.sigma, stat.chi2,
			    					 global_lap / 60.0, events / global_lap);
		}
		if (GMODE == 1) {
			const double global_lap = global_tictoc.ElapsedSec() - time_t0;
			double outputfilesize = gra::aux::GetFileSize(FULL_OUTPUT_STR) / GB;

			printf(
			    "[%0.1f MB/%0.2f GB] E: %9d, xs: %9.3E, er: %7.5f, %0.1f/%0.1f min ~ %0.1E Hz \n",
			    resident_use, outputfilesize, events, stat.sigma, stat.sigma_err / stat.sigma,
			    global_lap / 60.0,
			    (N - events) * global_lap / (double)events / 60.0,
			    events / global_lap);
		}
	}
}


// Final statistics
void MGraniitti::PrintStatistics(unsigned int N) {

	gra::aux::ClearProgress(); // Clear progressbar

	if (GMODE == 0) {
		time_t0 = global_tictoc.ElapsedSec();
		gra::aux::PrintBar("=");
		std::cout << rang::style::bold;
		printf("Monte Carlo integration: \n\n\n");
		std::cout << rang::style::reset;

		if (proc->GetIsolate()) {
			std::cout << rang::fg::red << "NOTE: Central leg phase space isolation tag &> in use!";
			std::cout << rang::fg::reset;
			std::cout << std::endl << std::endl;
		}

		unsigned int N_leg = proc->lts.decaytree.size() + 2;
		if (proc->GetIsolate()) { N_leg = 3; } // Isolated 2->3 process with <F> phase space

		printf("{2->%d cross section}:             %0.3E +- %0.3E barn \n", N_leg, stat.sigma, stat.sigma_err);
		printf("Sampling uncertainty:             %0.3f %% \n", 100 * stat.sigma_err / stat.sigma);
		std::cout << std::endl;

		if (proc->lts.DW_sum.Integral() > 0) { // Recursive phase space on

			std::cout << std::endl;
			const double PSvolume             = proc->lts.DW_sum.Integral();
			//const double PSvolume_error       = proc->lts.DW_sum.IntegralError();
			
			const double PSvolume_exact       = proc->lts.DW_sum_exact.Integral();
			const double PSvolume_exact_error = proc->lts.DW_sum_exact.IntegralError();

			// Get sum of decay daughter masses
			double MSUM = 0.0;
			for (std::size_t i = 0; i < proc->lts.decaytree.size(); ++i) {
				MSUM += proc->lts.decaytree[i].m_offshell;
			}
			// only 2-body exact or massless case
			if (proc->lts.decaytree.size() == 2 || MSUM < 1e-6) {
				printf(
				    "Analytic phase space volume:      %0.3E +- "
				    "%0.3E \n",
				    PSvolume_exact, PSvolume_exact_error);
				printf(
				    "RATIO: MC/analytic:               %0.6f \n",
				    PSvolume / PSvolume_exact);
			}
			std::cout << std::endl;

			// Collect out phase space weight recursively
			printf("{1->%lu LIPS}:                      %0.3E +- %0.3E \n",
			       	proc->lts.decaytree.size(),
			       	proc->lts.DW_sum.Integral(),
			       	proc->lts.DW_sum.IntegralError());

			for (std::size_t i = 0; i < proc->lts.decaytree.size(); ++i) {
					proc->PrintPhaseSpace(proc->lts.decaytree[i]);
					std::cout << std::endl;
			}

			// Recursion relation based (phase space factorization):
			// d^N PS(s; p_1, p2, ...p_N) = (1/(2*PI)) * d^3 PS(s;
			// p1,p2,p3) d^{N-2} PS(M^2; p4,p5,..,pN) dM^2
			
			// Decaywidth = 1/(2M S) \int dPS |M_decay|^2, where M
			// = mother mass, S = final state symmetry factor

			if (!proc->GetIsolate()) { // Not isolated
			printf("{2->3 cross section ~=~ [2->%u / (1->%lu LIPS) x 2PI]}:       %0.3E barn \n",
				N_leg, proc->lts.decaytree.size(), stat.sigma / proc->lts.DW_sum.Integral() * (2*PI) );
			std::cout << std::endl;
			printf("** Remember to use &> operator instead of -> for phase space isolation ** \n");
			}
		}
		gra::aux::PrintBar("-");

		std::cout << std::endl;
		printf("Integration runtime:              %0.2f sec \n",
		       time_t0);
		printf("Integrand sampling frequency:     %0.2E Hz \n",
		       (stat.evaluations - stat.trials) / time_t0);
		std::cout << std::endl;
		printf("<Statistics> \n");
		printf("MAX raw integrand weight:         %0.3E \n", stat.maxW);
		if (INTEGRATOR == "VEGAS") {
		printf("MAX weight (vegas x integrand):   %0.3E \n", stat.maxf);
		}

		printf("\n");
		printf("Amplitude passing rate:           %0.3E \n", stat.amplitude_ok / stat.evaluations);
		printf("Kinematics passing rate:          %0.3E \n", stat.kinematics_ok / stat.evaluations);
		printf("Fiducial cuts passing rate:       %0.3E \n", stat.fidcuts_ok / stat.evaluations);	
		printf("Veto cuts passing rate:           %0.3E \n", stat.vetocuts_ok / stat.evaluations);	
		printf("\n");

		std::cout << std::endl;
		printf(
		    "** All values include phase space generation and "
		    "fiducial cuts ** \n");
		gra::aux::PrintBar("=");
	}
	if (GMODE == 1) {
		const double lap = global_tictoc.ElapsedSec() - time_t0;

		gra::aux::PrintBar("=");
		if (WEIGHTED) {
			std::cout << rang::fg::green << "Weighted event generation:"
				<< std::endl << std::endl << std::endl;
		} else {
			std::cout << rang::fg::red   << "Unweighted event generation:"
				<< std::endl << std::endl << std::endl;
		}
		std::cout << rang::style::reset;

		printf("Generation efficiency:    %0.3E (%d / %0.0f) \n", N / (double) stat.trials, N, stat.trials);
		printf("Weight overflow:          %0.3E (%d / %0.0f) \n", stat.N_overflow / (double) stat.trials, stat.N_overflow, stat.trials);
		printf("Generation runtime:       %0.2f sec \n", lap);
		printf("Generation frequency:     %0.2E Hz \n", N / lap);

		double outputfilesize = gra::aux::GetFileSize(FULL_OUTPUT_STR) / (1024.0*1024.0*1024.0);

		printf("Outputfile size:          %0.3f GB [%s] \n", outputfilesize, FULL_OUTPUT_STR.c_str());
		gra::aux::PrintBar("=");
		std::cout << std::endl;
	}
}

} // gra namespace
