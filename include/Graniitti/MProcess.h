// Abstract process class
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.


#ifndef MPROCESS_H
#define MPROCESS_H

// C++
#include <complex>
#include <random>
#include <vector>

// HepMC33
#include "HepMC3/FourVector.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"

// Own
#include "Graniitti/MSubProc.h"
#include "Graniitti/MAux.h"
#include "Graniitti/MEikonal.h"
#include "Graniitti/MH1.h"
#include "Graniitti/MH2.h"
#include "Graniitti/MKinematics.h"
#include "Graniitti/MSudakov.h"
#include "Graniitti/MRandom.h"
#include "Graniitti/MUserHistograms.h"
#include "Graniitti/M4Vec.h"
#include "Graniitti/MPDG.h"


namespace gra {


// Event-by-event auxialary data for integration
struct AuxIntData {

	// Aux weight
	double vegasweight = 1.0;
	bool burn_in_mode  = false;

	// Event-by-event assertations (init all with true!)
	bool amplitude_ok  = true;
	bool kinematics_ok = true;
	bool fidcuts_ok    = true;
	bool vetocuts_ok   = true;
	
	// Forced acceptance of the event
	bool forced_accept = false;

	bool Valid() const {
		return kinematics_ok && fidcuts_ok && vetocuts_ok;
	}
};


// Multipomeron kinematics
struct MPI {

	M4Vec p1i;
	M4Vec p2i;

	M4Vec q1;
	M4Vec q2;

	M4Vec k;

	M4Vec p1f;
	M4Vec p2f;

	M4Vec p3;
	M4Vec p4;
};


// Abstract process class
class MProcess : public MUserHistograms {

public:

	// Of polymorphic type
	MProcess(){};
	virtual ~MProcess(){}; // MUST HAVE IT HERE as virtual

	// Pure virtual functions, without definitions here
	// WITHOUT =0, undefined reference to vtable for MProcess will occur in
	// compilation
	virtual   void post_Constructor() = 0;
	virtual   void PrintInit(bool silent) const = 0;

	virtual double operator ()(const std::vector<double>& randvec, AuxIntData& aux) = 0;
	virtual double EventWeight(const std::vector<double>& randvec, AuxIntData& aux) = 0;
	virtual   bool EventRecord(HepMC3::GenEvent& evt) = 0;

	std::vector<std::string> PrintProcesses() const;
	std::string GetProcessDescriptor(std::string process) const;
	
	// Set central system decay structure
	void SetDecayMode(std::string str);
	void SetupBranching();

	// Set initial state
	void SetInitialState(const std::vector<std::string>& beam, const std::vector<double>& energy);

	// ISOLATE phase space in <F> class processes
	void SetISOLATE(bool in) {
		ISOLATE = in;
	}
	bool GetISOLATE() {
		return ISOLATE;
	}
	void SetFLATMASS2(bool in) {
		if (in == true) {
		aux::PrintWarning();
		std::cout << rang::fg::red << "MProcess::SetFLATMASS2: Set flat in mass^2 sampling in decay trees : true" << rang::fg::reset << std::endl;
		}
		FLATMASS2 = in;
	}
	bool GetFLATMASS2() {
		return FLATMASS2;
	}
	void SetOFFSHELL(double in) {
		if (in == true) {
		aux::PrintWarning();
		std::cout << rang::fg::red << "MProcess::SetOFFSHELL: Set number of decay widths in decay trees : " << in << rang::fg::reset << std::endl;
		}
		OFFSHELL = in;
	}
	double GetOFFSHELL() {
		return OFFSHELL;
	}

	// Get initial state
	std::vector<gra::MParticle> GetInitialState() {
		std::vector<gra::MParticle> beams = {beam1, beam2};
		return beams;
	}

	// Phase space dimension
	unsigned int GetdLIPSDim() {
		return ProcPtr.LIPSDIM;
	}

	// Pomeron loop screening
	void SetScreening(bool value) {
		SCREENING = value;
	}
	bool GetScreening() {
		return SCREENING;
	}
	// Set/Get input eikonal
	void SetEikonal(const MEikonal& in) {
		Eikonal = in;
	}
	MEikonal GetEikonal() const {
		return Eikonal;
	}

	// Set LHAPDF
	void SetLHAPDF(const std::string& in) {
		gra::form::LHAPDF = in;
	}
	
	// Set cuts
	void SetGenCuts(const gra::GENCUT& in) {
		gcuts = in;
	}
	void SetFidCuts(const gra::FIDCUT& in) {
		fcuts = in;
	}
	void SetUserCuts(int in) {
		USERCUTS = in;
	}
	void SetVetoCuts(const gra::VETOCUT& in) {
		vetocuts = in;
	}

	double GetMandelstam_s() const {
		return lts.s;
	}

	// Set proton excitation to low-mass N*
	void SetExcitation(int in) {
		if (in > 2 || in < 0) {
			std::string str =
			    "MProcess::SetExcitation: Not valid input "
			    "(0,1,2): " +
			    std::to_string(in) + " !";
			throw std::invalid_argument(str);
		}
		EXCITATION = in;
	}

	// Set flat matrix element mode
	void SetFLATAMP(int in) {

		if (in > 0) {
		aux::PrintWarning();
		std::cout << rang::fg::red << "MProcess::SetFLATAMP: Flat matrix element FLATAMP : " << in << rang::fg::reset << std::endl;
		}
		FLATAMP = in;
	}

	// pp invariant Moller flux (high energy limit)
	double MollerFlux() {
		return 2.0 * lts.s;
	}

	// Flat amplitudes (for DEBUG)
	double GetFlatAmp2(const gra::LORENTZSCALAR& lts) const;

	// Set input resonances
	void SetResonances(const std::map<std::string, gra::PARAM_RES>& in) {
		lts.RESONANCES = in;
	}
	
	// Get input resonances
	std::map<std::string, gra::PARAM_RES> GetResonances() const {
		return lts.RESONANCES;
	}

	// Eikonal (screening) functions
	MEikonal Eikonal;

	// Check processes
	bool ProcessExist(std::string process) const;

	// Lets keep these public for easy access
	gra::LORENTZSCALAR lts; // Lorentz scalars and others for kinematics

	// Cut structures
	gra::GENCUT gcuts;      // Generator sampling cuts (phase space boundaries)
	gra::FIDCUT fcuts;      // Fiducial cuts (phase space boundaries)
	gra::VETOCUT vetocuts;  // Veto cuts

	void PrintDecayTree(const gra::MDecayBranch& branch) const;
	void PrintPhaseSpace(const gra::MDecayBranch& branch) const;

	// Random numbers (keep it public for seeding)
	MRandom random;


protected:

	// Copy and assignment made private
	//MProcess(const MProcess& other);
	//MProcess& operator=(const MProcess& rhs);

	// Internal virtual functions, without definitions here
	virtual void ConstructProcesses() = 0;
	virtual bool LoopKinematics(const std::vector<double>& p1p,
	                            const std::vector<double>& p2p) = 0;
	virtual bool FiducialCuts() const = 0;
	
	// Eikonal screening loop
	std::complex<double> S3ScreenedAmp();
	
	// First print
	void PrintSetup() const;

	// Setup process
	void SetProcess(std::string& process, const std::vector<aux::OneCMD>& syntax);

	// QFT symmetry factor
	void CalculateSymmetryFactor();


	// -------------------------------------------------------
	// Recursive function to treat decay trees

	bool CommonRecord(HepMC3::GenEvent& evt);
	bool VetoCuts() const;
	bool CommonCuts() const;
	void FindDecayCuts(const gra::MDecayBranch& branch, bool& ok) const;
	void FindVetoCuts(const gra::MDecayBranch& branch, bool& ok) const;
	bool ConstructDecayKinematics(gra::MDecayBranch& branch);
	void WriteDecayKinematics(gra::MDecayBranch& branch,
	                          HepMC3::GenParticlePtr& mother,
	                          HepMC3::GenEvent& evt);
	void PrintFiducialCuts() const;
	
	
	// Offshell mass pick
	void GetOffShellMass(const gra::MDecayBranch& branch, double& mass);
	
	// Lorentz scalars
	bool GetLorentzScalars(unsigned int Nf);
	
	void SampleForwardMasses(std::vector<double>& mvec, const std::vector<double>& randvec);

	// --------------------------------------------------------
	// System fragmentation

	void ExciteNstar(const M4Vec& nstar, gra::MDecayBranch& forward);

	int  ExciteContinuum(const M4Vec& nstar, const HepMC3::GenParticlePtr gen_nstar,
						HepMC3::GenEvent& evt, double Q2_scale, int B_sum, int Q_sum);

	void BranchForwardSystem(const std::vector<M4Vec>& products,
							 const std::vector<int>& pdgcode, const std::vector<bool>& isstable,
							 const HepMC3::GenParticlePtr gen_nstar, HepMC3::GenEvent& evt);

	// ---------------------------------------------------------
	
	// Check std::nan/std::inf
	bool CheckInfNan(double& W) {
		if        (std::isnan(W)) {
			++N_nan;
			W = 0;
			return false;
		} else if (std::isinf(W)) {
			++N_inf;
			W = 0;
			return false;
		}
		if (N_nan > 50) {
			throw std::invalid_argument("MProcess::CheckInfNan: Too many NaN weights - Check model parameters and cuts!");
		}
		if (N_inf > 50) {
			throw std::invalid_argument("MProcess::CheckInfNan: Too many Inf weights - Check model parameters and cuts!");
		}
		return true;
	}
	unsigned int N_inf = 0;
	unsigned int N_nan = 0;
	
	// Initial states
	gra::MParticle beam1;
	gra::MParticle beam2;
	double S_factor = 0.0;        // Cross section statistical 1/S symmetry factor
	                              // (for identical final states)	
	// ----------------------------------------------------------------------
	// Subprocesses
	
	// Process ID lists
	std::map<std::string, std::string> Processes;
	MSubProc ProcPtr;
	
	// ----------------------------------------------------------------------
	// Steering parameters

	std::string PROCESS;          // Process identifier string
	std::string CID;              // Phase space sampler identifier <F>,<C>,...
	std::string DECAYMODE;        // Decaymode identifier string
	bool SCREENING = false;       // Pomeron loop on/off
	unsigned int EXCITATION = 0;  // Forward proton excitation (0 = off, 1 = single, 2 = double)
	unsigned int USERCUTS = 0;    // User custom cuts identifier
	unsigned int FLATAMP = 0;     // Flat matrix element mode
	
	// ----------------------------------------------------------------------
	// Phase-space control

	bool ISOLATE    = false;      // ISOLATE the decay phase space (<F> class)
	bool FLATMASS2  = false;      // Flat in M^2 instead of Breit-Wigner sampling
	double OFFSHELL = 15;         // How many full widths to sample particles in trees
	                              // (note that small values can create artifacts!)
	
	// Forward excitation minimum/maximum M^2 boundaries
	double M2_f_min = 0.0;
	double M2_f_max = 0.0;
	// ----------------------------------------------------------------------
	
	// Non-Diffractive
	std::vector<MPI> etree;
	double bt = 0.0;
	
	// Particle database
	MPDG PDG;
};


} // gra namespace ends


#endif