// Abstract process class
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MPROCESS_H
#define MPROCESS_H

// C++
#include <array>
#include <cmath>
#include <complex>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

// HepMC3
#include "HepMC3/FourVector.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"

// Own
#include "Graniitti/Nuclear/MUPCScreen.h"
#include "Graniitti/Eikonal/MEikonal.h"
#include "Graniitti/Kinematics/M4Vec.h"
#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/MUserHistograms.h"
#include "Graniitti/Math/MFloat.h"
#include "Graniitti/Math/MStatistics.h"
#include "Graniitti/Process/MHelicityConfig.h"
#include "Graniitti/Process/MProcessState.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MException.h"

namespace gra {

// Abstract process interface for concrete event generators
class MProcess : public MUserHistograms {
 public:
  // Of polymorphic type
  virtual ~MProcess() {}

  // Enable failure diagnostics before making worker copies
  void SetDebug(bool enabled) noexcept { debug = enabled; }

  // Complete configuration and mutable data copied independently per worker
  MProcessState state;

  // Selected physical subprocess and worker-local amplitude implementation
  MSubProc ProcPtr;

  // Worker-local soft-screening service
  MEikonal eikonal;

  // Load process-local immutable helicity steering cards
  void SetHelicityConfig(MModelTunePtr tune);

  // Set central system structure
  void SetDecayMode(std::string str);

  // Initialize the complete process runtime before worker copies are made
  void PrepareRun();

  // Initialize the complete process runtime with an external eikonal
  void PrepareRun(const MEikonal &eikonal);

  // Initialize only the selected physical amplitude and its immutable data
  void InitializeProcessAmplitude();

  // Construct one helicity vertex from the process-local steering cards
  HELMatrix ProcessHelicityStructure(const MParticle &particle, const std::vector<MParticle> &legs,
                                     bool production_mode = false, bool strict_mode = false,
                                     const std::string &verbose_label = "", bool verbose_output = true,
                                     gra::spin::VertexContext context = gra::spin::VertexContext::Auto) const;

  // Construct one shared physical-pole LS vertex from its model card
  gra::spin::PoleLS ProcessPoleOperatorStructure(ReggeProductionModel model, const MParticle &particle,
                                                            const std::vector<MParticle> &legs,
                                                            gra::spin::VertexContext      context) const;

  // Construct helicity vertices recursively through one decay branch
  void ProcessHelicityTree(MDecayBranch &branch, bool production_mode = false, bool strict_mode = false);

  // Set isolated phase-space sampling
  void SetISOLATE(bool value);
  bool GetISOLATE() const;

  // Set root decay semantics selected by the process arrow
  void SetRootDecayMode(RootDecayMode mode);

  // Set cascade phase-space controls
  void   SetFLATMASS2(bool value);
  bool   GetFLATMASS2() const;
  void   SetOFFSHELL(double value);
  double GetOFFSHELL() const;
  void   SetDefaultOFFSHELL(double value);
  void   SetWIDTHMIN(double value);
  double GetWIDTHMIN() const;

  // Set helicity-amplitude controls
  void SetSPINGEN(bool value);
  void SetSPINDEC(bool value);
  // Enable spin metrics independently of physical spin correlations
  void SetQMetrics(bool value);
  void SetMPFrame(const std::string &frame);
  void SetMMAX(int value);

  // Compute the configured initial state
  std::vector<MParticle> GetInitialState() const;

  // Compute the phase-space dimension
  unsigned int GetdLIPSDim() const;

  // Set screening and coupled eikonal data
  void            SetScreening(bool value);
  bool            GetScreening() const;
  void            SetEikonal(const MEikonal &eikonal);
  const MEikonal &GetEikonal() const;

  // Set the immutable SOFT model shared by amplitudes and screening
  void                SetSoftModel(SoftModelPtr model);
  const SoftModelPtr &GetSoftModel() const noexcept;

  // Set the immutable complete model tune shared by all workers
  void SetModelTune(MModelTunePtr tune);

  // Access the immutable complete model tune
  const MModelTunePtr &GetModelTune() const noexcept;

  // Set the process PDF set
  void SetLHAPDF(const std::string &name);

  // Set process cuts
  void SetGenCuts(const GENCUT &cuts);
  void SetFidCuts(const FIDCUT &cuts);
  void SetUserCuts(std::int64_t id);
  void SetVetoCuts(const VETOCUT &cuts);

  // Select beam remnant fragmentation
  void SetBeamFrag(const std::string &mode);

  // Set forward excitation and debug amplitude controls
  void SetExcitation(int value);
  void SetFLATAMP(int value);

  // Set and return input resonances
  void                                    SetResonances(const std::map<std::string, PARAM_RES> &resonances);
  const std::map<std::string, PARAM_RES> &GetResonances() const;

  // Pure virtual process hooks
  virtual void PrintInit(bool silent) const = 0;

  // Evaluate and classify one process weight through the common sampling
  // interface
  double EventWeight(const std::vector<double> &randvec, MEventWeightState &aux);

  // Construct one accepted event record through the common cleanup boundary
  bool EventRecord(HepMC3::GenEvent &evt);

  // Set the initial state using per-nucleon energies for nuclear beams
  void SetInitialState(const std::vector<std::string> &beam, const std::vector<double> &energy);

  // Set the heavy-ion UPC model after both beam particles are known
  void SetNuclear(const nuclear::UPCParam &param, std::unique_ptr<nuclear::MReaction> model = nullptr);

  // Attach a worker-cloned nuclear reaction model to factorized or central phase space
  void SetNuclearFinalState(std::unique_ptr<nuclear::MReaction> model);

  // Compute whether a heavy-ion UPC model is active
  bool HasNuclear() const noexcept;

  // Set full beam-particle energies
  void SetBeamEnergies(double E1, double E2);

  // Set the QED initial and final state radiation configuration
  void SetRadiative(const radiative::Config &config);

  // Add configured QED photons to one accepted event record
  bool ApplyRadiation(HepMC3::GenEvent &evt) noexcept;

  // Validate fiducial cuts against the configured phase-space class
  void ValidateFiducialCuts(const gra::FIDCUT &cuts, std::int64_t usercuts) const;

  // Flat amplitudes for debug
  double GetFlatAmp2(const gra::LORENTZSCALAR &lts) const;

 protected:
  static constexpr double ZERO_EPS = 1e-12;

  // Construct the process-independent amplitude initialization context
  MProcessSetup CreateProcessSetup();

  // Calculate the QFT symmetry factor for the configured stable final state
  void CalculateSymmetryFactor();

  // Finalize phase-space-class state after amplitude initialization
  virtual void FinalizeProcessConfiguration() = 0;

  // Compute the minimum eikonal kT2 table maximum required by this process
  virtual double EikonalMaxKT2() const;

  // Compute one process-specific weight before common numerical classification
  virtual double ComputeEventWeight(const std::vector<double> &randvec, MEventWeightState &aux) = 0;

  // Construct one process-specific event record
  virtual bool BuildEventRecord(HepMC3::GenEvent &evt);

  // Copy and assignment made private
  // MProcess(const MProcess& other);
  // MProcess& operator=(const MProcess& rhs);

  // Internal virtual functions
  virtual bool LoopKinematics(const std::array<double, 2> &p1p, const std::array<double, 2> &p2p) = 0;
  virtual bool FiducialCuts() const;

  // Restore the accepted Born four-momenta and refresh process-specific
  // kinematics
  bool RestoreBornKinematics();

  // Refresh process-specific quantities from the restored Born four-momenta
  virtual bool RefreshBornKinematics();

  // Reconstruct and evaluate one common pp or UPC screening node
  bool EvaluateScreeningNode(const std::array<double, 2> &p1p, const std::array<double, 2> &p2p);

  // Compute one screening node through the common failure boundary
  bool ComputeScreeningNode(const std::array<double, 2> &p1p, const std::array<double, 2> &p2p);

  // Cascade phase-space factor
  double CascadePS();

  // Cache stable-leaf assignment plan for symmetrized or diagnostic decay
  // handling
  void CacheStableLeafSymmetryAssignments();

  // Compute incoherent stable-leaf assignment compensation for diagnostic
  // unsymmetrized decays
  double DecaySymmetryCompensationFactor();

  // Select one stable-leaf proposal channel for symmetrized cascade decays
  void PrepareDecaySymmetryProposal();

  // Apply the selected stable-leaf proposal channel to generated kinematics
  bool ApplyDecaySymmetryProposal();

  // Amplitude squared with optional eikonal screening
  double GetAmp2(bool include_screening, MEventWeightState &aux);

  // Sample nuclear fluctuations once for the active event evaluation
  void PrepareUPCEvent(bool include_screening);

  // Apply generic failure bookkeeping around one configured amplitude
  // evaluation
  double EvaluateAmplitudeBoundary(bool include_screening, MEventWeightState &aux);

  // Eikonal screening loop
  double ScreenedAmplitudeSquared(bool include_screening);

  // Screen one Good Walker pair-space amplitude with the pp eikonal
  double ScreenPPGoodWalker(ProtonGoodWalkerAmplitude born, const ScreeningMetadata &metadata,
                            const std::array<double, 2> &upper_transverse,
                            const std::array<double, 2> &lower_transverse, const MEikonal::LoopConst &loop_const);

  // Screen physical helicity amplitudes with the pp eikonal
  double ScreenPPAmplitude(const std::array<double, 2> &upper_transverse, const std::array<double, 2> &lower_transverse,
                           const MEikonal::LoopConst &loop_const);

  // Screen one heavy-ion UPC amplitude with a scalar Glauber kernel
  double ScreenUPCAmplitude(const std::array<double, 2> &upper_transverse,
                            const std::array<double, 2> &lower_transverse, const nuclear::MUPC &upc);

  // Normalize a resolved nuclear projection and finalize event weights
  double FinalizeUPCAmplitude(const nuclear::ScreenResult &result, const nuclear::ScreenLayout &layout,
                              const ScreeningMetadata &metadata);

  // Evaluate one amplitude and convert status-based failures to the generic
  // exception
  double EvaluateAmplitude();

  // Validate one Born or screening amplitude at the shared sampling boundary
  void ValidateAmplitude(double value) const;

  // Restore Born kinematics after screening without allowing cleanup to throw
  bool RestoreScreeningBornState() noexcept;

  // Roll back one failed screening evaluation through the generic cleanup
  // interface
  void AbortScreeningLoopState() noexcept;

  // Restore Born kinematics automatically unless one screening loop succeeds
  class ScreeningLoopGuard {
   public:
    // Mark the owner as evaluating shifted screening kinematics
    explicit ScreeningLoopGuard(MProcess &owner) : owner_(owner) { owner_.state.lts.screening.active = true; }

    // Restore and clear an incomplete screening evaluation
    ~ScreeningLoopGuard() {
      if (!accepted_) { owner_.AbortScreeningLoopState(); }
    }

    // Preserve the completed Born state after successful projection
    void Accept() noexcept { accepted_ = true; }

    ScreeningLoopGuard(const ScreeningLoopGuard &)            = delete;
    ScreeningLoopGuard &operator=(const ScreeningLoopGuard &) = delete;

   private:
    MProcess &owner_;
    bool      accepted_ = false;
  };

  // Clear amplitude outputs after their event-local consumers have finished
  void ResetTransientAmplitudeState() noexcept;

  // Clear incomplete amplitude outputs and every transient screening cache
  void ResetFailedAmplitudeState() noexcept;

  // Validate and normalize orthogonal nuclear final-sector weights
  void SetUPCFinalWeights(const nuclear::FinalWeights &weight);

  // Sample the event-local nuclear final sector at most once
  bool SampleUPCFinalState();

  // Finalize nuclear, color-flow and process-specific event state atomically
  bool FinalizeEventState();

  // Clear event-local nuclear final-sector probabilities and selection
  void ResetUPCFinalState() noexcept;

  // Receive final amplitude-component weights after optional screening
  virtual void PostScreeningAmplitude(const std::vector<double> &component_amp_squared) { (void)component_amp_squared; }

  // Finalize event-local amplitude state through the common failure boundary
  virtual bool FinalizeAmplitudeEventState() { return ProcPtr.SampleColorFlow(state.lts); }

  // Evaluate bare subprocess amplitude
  virtual double EvaluateBareAmplitude() {
    ProcPtr.BindRandom(state.random);
    const double value =
        state.lts.screening.active ? ProcPtr.ScreeningAmp2(state.lts) : ProcPtr.GetBareAmplitude2(state.lts);
    evaluation_status = ProcPtr.EvaluationStatus();
    return value;
  }

  // First print
  void PrintSetup() const;

  // Setup process
  void SetProcess(std::string &process, const std::vector<aux::OneCMD> &syntax, MModelTunePtr tune = {});

  // Prepare common state after one physical phase-space construction
  void PreparePhaseSpacePoint(bool kinematics_ok, MEventWeightState &aux);

  // Apply FSR and the final fiducial and veto classification
  void ClassifyFinalState(MEventWeightState &aux);

  // Apply common cuts
  bool CommonCuts() const;

  // Construct recursive decay kinematics
  bool ConstructDecayKinematics(gra::MDecayBranch &branch);

  // Print fiducial cuts
  void PrintFiducialCuts() const;

  // Sample an offshell mass for a branch
  void GetOffShellMass(gra::MDecayBranch &branch, double &mass);

  // Prepare top-level branch masses without independently sampling a sole root
  bool PrepareCentralBranchMasses(double &mass_sum, double *mass_max = nullptr);

  // Bind one sole production root to the complete central-system momentum
  bool SetSingleCentralRootKinematics();

  // Compute the direct central phase-space volume for the current root count
  double CentralDecayWidthPS() const;

  // Set technical phase-space boundaries
  void SetTechnicalBoundaries(gra::GENCUT &gcuts, unsigned int EXCITATION);

  // Compute forward phase-space volume
  double ForwardVolume() const;

  // Compute the inverse probability of the sampled single dissociation side
  double DissociationCrossSectionFactor() const;

  // Sample forward masses
  void SampleForwardMasses(std::vector<double> &mvec, const std::vector<double> &randvec);

  // Excite low-mass N* system
  bool ExciteNstar(const M4Vec &nstar, gra::MDecayBranch &forward, const MParticle &pbeam);

  // Build one showerable color-singlet quark-diquark forward system
  bool ExciteString(const M4Vec &nstar, gra::MDecayBranch &forward, const MParticle &pbeam, int color_tag);

  // Excite continuum forward system
  bool ExciteContinuum(const M4Vec &nstar, gra::MDecayBranch &forward, double Q2_scale, int B_sum, int Q_sum);

  // Branch one forward excitation system transactionally
  bool BranchForwardSystem(const std::vector<M4Vec> &p4, const std::vector<MParticle> &p, const M4Vec &nstar,
                           gra::MDecayBranch &forward);

  // Clear event-local color and parton metadata
  void ResetEventRecordState() noexcept;

  // Initialize all generic event-local sampling state
  void BeginEventEvaluation(MEventWeightState &aux) noexcept;

  // Clear generic state belonging to a rejected or interrupted event
  void ResetRejectedEventState() noexcept;

  // Construct CEP forward fragmentation
  bool CEPForwardFragment();

  // Parse process command syntax
  void ParseCMD(const std::string &str, std::string &first, std::string &second, std::string &third) const;

  // Check std::nan/std::inf
  bool CheckInfNan(double &W) {
    if (std::isnan(W)) {
      ++N_nan;
      W = 0;
      return false;
    } else if (std::isinf(W)) {
      ++N_inf;
      W = 0;
      return false;
    }
    return true;
  }

  // Separate invalid numerical weights from zero physical amplitude support
  void BookkeepAmplitudeWeight(double &W, MEventWeightState &aux);

  unsigned int               N_inf                           = 0;
  unsigned int               N_nan                           = 0;
  mg5helas::EvaluationStatus evaluation_status               = mg5helas::EvaluationStatus::Success;
  bool                       amplitude_event_state_finalized = false;
  // Forward excitation minimum/maximum M^2 boundaries
  double M2_f_min     = 0.0;
  double M2_f_max     = 0.0;
  double log_M2_f_min = 0.0;
  double log_M2_f_max = 0.0;

 private:
  // Print one nonfatal failure without changing sampling or throwing
  void PrintFailure(const char *type, const char *message) const noexcept;

  bool debug = false;

  // Immutable parsed helicity cards shared safely across worker copies
  MHelicityConfig helicity_config;
};

}  // namespace gra

#endif
