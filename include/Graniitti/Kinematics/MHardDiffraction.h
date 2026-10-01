// Hard diffraction with factorized or adaptive central phase space
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MHARDDIFFRACTION_H
#define MHARDDIFFRACTION_H

// C++
#include <array>
#include <map>
#include <memory>
#include <random>
#include <string>
#include <vector>

// Own
#include "Graniitti/Amplitude/MG5/Runtime/AMP_MG5_Helicity.h"
#include "Graniitti/Kinematics/M4Vec.h"
#include "Graniitti/PDF/MHardPomeronPDF.h"
#include "Graniitti/Process/MProcess.h"
#include "Graniitti/Tech/MAux.h"

// HepMC3
#include "HepMC3/GenEvent.h"

namespace gra {

// Driver for hard Pomeron PDF based diffraction
class MHardDiffraction : public MProcess {
 public:
  // Construct the process list
  MHardDiffraction();

  // Construct one selected hard-diffraction process
  MHardDiffraction(std::string process, const std::vector<aux::OneCMD> &syntax, MModelTunePtr tune = {});

  // Destroy one selected hard-diffraction process
  virtual ~MHardDiffraction();

  // Initialize process-specific phase-space dimensionality
  void FinalizeProcessConfiguration() override;

  // Print the hard-diffraction setup
  void PrintInit(bool silent) const;

 private:
  // Initialize supported hard-diffraction process names
  void Initialize(const std::string &mode = "F");

 protected:
  // Compute one hard diffraction process weight
  double ComputeEventWeight(const std::vector<double> &randvec, MEventWeightState &aux) override;

  // Construct one hard diffraction event record
  bool BuildEventRecord(HepMC3::GenEvent &evt) override;

  // Update loop kinematics for the eikonal screening integration
  bool LoopKinematics(const std::array<double, 2> &p1p, const std::array<double, 2> &p2p);

  // Refresh hard-diffraction kinematics from the restored Born four-momenta
  bool RefreshBornKinematics() override;

  // Build t-dependent hard-scattering kinematics
  bool HardBuildKin(double xhard1, double xhard2);

  // Build one hard parton pair with the exact sampled DPDF mass
  bool BuildHardPair(double xhard1, double xhard2, M4Vec &q1, M4Vec &q2) const;

  // Scale two diffractive hard legs symmetrically to the exact DPDF mass
  static bool ScaleDoubleDiffractivePair(double mass2, double beta1, double beta2, M4Vec &q1, M4Vec &q2);

  // One collinear hard-system point and its longitudinal Jacobian
  struct HardLongPoint {
    double mass2    = 0.0;
    double rapidity = 0.0;
    double xhard1   = 0.0;
    double xhard2   = 0.0;
    double jacobian = 0.0;
    double mass_jacobian = 0.0;
  };

  // Map one unit coordinate logarithmically onto a positive interval
  static bool MapLog(double unit, double min, double max, double &value, double &jacobian);

  // Map one unit coordinate with density proportional to exp(slope value)
  static bool MapExp(double unit, double min, double max, double slope, double &value, double &jacobian);

  // Map collinear hard fractions through mass squared and rapidity
  static bool MapHardLong(double mass_unit, double rapidity_unit, double s, double mass2_min, double mass2_max,
                          const std::array<double, 2> &xhard1_range, const std::array<double, 2> &xhard2_range,
                          double xi_product, HardLongPoint &point);

 private:
  // Apply common fiducial cuts
  bool FiducialCuts() const;

  // Apply hard-diffraction generation cuts to the sampled hard system
  bool HardGenerationCuts() const;

  // Evaluate the PDF-weighted hard subprocess amplitude
  double EvaluateBareAmplitude() override;

  // Skip default helicity averaging because hard-channel components are already
  // MG5 averaged
  // Sample one incoming hard-parton channel from final component weights
  void PostScreeningAmplitude(const std::vector<double> &component_amp_squared) override;

  // Finalize the selected hard channel and both remnant branches atomically
  bool FinalizeAmplitudeEventState() override;

  // Build one random hard-diffraction phase-space point
  bool HardRandomKin(const std::vector<double> &randvec);

  // Sample one diffractive external leg and return its xi-t Jacobian
  bool SampleDiffractiveSide(bool first_side, double xi_unit, double t_unit, double phi_unit, double &jacobian);

  // Compute the hard-fraction support for one current beam side
  std::array<double, 2> HardXRange(bool first_side) const;

  // Store one sampled hard fraction as beta or ordinary proton x
  bool SetHardX(bool first_side, double xhard);

  // Sample both hard fractions through mass squared and rapidity
  bool SampleHardLong(double mass_unit, double rapidity_unit, double mass_sum, double mass_max, double &jacobian);

  // Build one incoming hard parton for a selected beam side
  bool BuildHardPartonForSide(bool first_side, M4Vec &parton) const;

  // Build the tagged particle and physical exchange remnant on one side
  bool BuildTaggedRemnant(bool first_side, const M4Vec &parton, M4Vec &leading, M4Vec &remnant) const;

  // Check one hard parton against its physical beam-remnant support
  bool AcceptHardParton(bool first_side, const M4Vec &parton) const;

  // Rebuild deterministic forward remnant decay branches
  bool BuildHardRemnantBranches(bool assign_colors);

  // Build one four-vector from light-cone components and transverse momentum
  M4Vec LightConeVector(double plus, double minus, double px, double py) const;

  // Build one spacelike hard parton and a physical lightlike Pomeron remnant
  bool BuildDiffractiveHardParton(const M4Vec &beam, const M4Vec &exchange, double xi, double beta, bool plus_side,
                                  M4Vec &parton) const;

  // Build the ordinary hard leg with a massive remnant and the sampled hard mass
  bool BuildProtonHardParton(bool first_side, const M4Vec &qdiff, double mass2, M4Vec &parton) const;

  // Resolve the ordinary proton remnant mass for this process
  double ProtonRemnantMass() const;

  // Split one ordinary remnant using the fixed event decay angles
  bool BuildProtonRemnant(bool first_side, MDecayBranch &branch) const;

  // Rebuild nominal t-dependent forward remnant branches
  bool BuildNominalHardRemnantBranches(std::array<MDecayBranch, 2> &branches);

  // Rebuild one nominal forward remnant side
  bool BuildNominalHardRemnantSide(bool first_side, MDecayBranch &branch);

  // Rebuild loop-deformed forward remnant branches by momentum fractions
  bool BuildLoopHardRemnantBranches(std::array<MDecayBranch, 2> &branches);

  // Rebuild one loop-deformed forward remnant side
  bool BuildLoopHardRemnantSide(bool first_side, MDecayBranch &branch);

  // Assign shower-compatible color tags to hard partons and remnants
  bool AssignHardColorFlow();

  // Enumerate leading-color assignments for the current hard event
  std::vector<std::vector<MColorFlow>> BuildHardColorFlowCandidates();

  // Map one exact generated external flow onto central and remnant leaves
  bool BuildGeneratedHardColorFlowCandidate(const mg5helas::ExternalColorFlow &external_flow,
                                            std::vector<MColorFlow>           &candidate);

  // Convert one generated symbolic color pair to event color tags
  MColorFlow ConvertGeneratedColorFlowLeg(const mg5helas::ColorFlowLeg &leg, std::map<int, int> &tag_map) const;

  // Distribute the crossed incoming color across one forward remnant
  bool MapRemnantColor(const std::vector<MDecayBranch *> &leaves, const MColorFlow &crossed,
                       std::size_t offset, int &next_tag, std::vector<MColorFlow> &candidate) const;

  // Compute true when one event color pair matches the particle representation
  bool HardColorFlowMatchesParticle(const MColorFlow &flow, const MParticle &particle) const;

  // Validate one complete candidate without changing the event color state
  bool ValidateHardColorFlowCandidate(const std::vector<MColorFlow> &candidate);

  // Compute true when two hard color banks have identical channels and tags
  bool HardColorFlowsMatch(const std::vector<mg5helas::HardColorFlow> &first,
                           const std::vector<mg5helas::HardColorFlow> &second) const;

  // Apply one generated color-flow assignment in stable-leaf order
  bool ApplyHardColorFlowCandidate(const std::vector<MColorFlow> &candidate);

  // Collect stable leaves from one mutable decay branch
  void CollectStableLeaves(MDecayBranch &branch, std::vector<MDecayBranch *> &leaves) const;

  // Compute true for partons carrying QCD color
  bool IsColoredParton(int pdg) const;

  // Compute true for colored diquark remnant ids
  bool IsDiquark(int pdg) const;

  // Compute a standard colored remnant partner for an extracted parton
  int RemnantPartnerPDG(int pdg) const;

  // Compute the ordinary remnant flavours with baryon number and charge conservation
  std::array<int, 2> ProtonRemnantIDs(int pdg, int beam_pdg) const;

  // Build one stable remnant branch
  MDecayBranch MakeRemnantBranch(int pdg, const M4Vec &p4) const;

  // Compute true when one remnant four-vector is usable in an event record
  bool AcceptRemnantMomentum(const M4Vec &p4) const;

  // Compute a particle definition or a local fallback for remnant bookkeeping
  MParticle HardParticle(int pdg) const;

  // Compute the sampled hard-diffraction integral volume
  double HardIntegralVolume() const;

  // Compute the hard-scattering phase-space Jacobian
  double HardPhaseSpaceWeight() const;

  // Calculate pure phase-space decay width
  void DecayWidthPS(double &exact) const;

  // Compute the supported incoming hard-parton flavours
  std::vector<int> HardPartonFlavours() const;

  // One PDF-positive incoming hard-parton channel
  struct HardPartonChannel {
    int    id1     = 0;
    int    id2     = 0;
    double f1      = 0.0;
    double f2      = 0.0;
    double pdf_xf1 = 0.0;
    double pdf_xf2 = 0.0;
  };

  // Publish the tagged leading proton losses independently of hard parton x
  void SetTaggedForwardXi();

  // Build the incoming hard-parton basis for the current phase-space point
  std::vector<HardPartonChannel> BuildHardPartonChannels(double Q2);

  // Recalculate channel densities without changing the cached channel ordering
  std::vector<HardPartonChannel> RefreshHardPartonChannelDensities(const std::vector<HardPartonChannel> &channels,
                                                                   double                                Q2);

  // Compute the current diffractive proton momentum transfer for one beam side
  double CurrentDiffractiveT(bool first_side) const;

  // Compute the PDF or DPDF density for one incoming beam side
  double HardPartonDensityForSide(int pid, bool first_side, double Q2);

  // Select and publish one hard-parton channel to lts metadata
  bool SelectHardPartonChannel(const std::vector<double> &component_amp2);

  // Compute a proton parton density f_i/p(x,Q2)
  double ProtonPartonDensity(int pid, double x, double Q2);

  // Compute alpha_s for the current hard scale
  double AlphaQCD(double Q2);

  std::shared_ptr<const MHardPomeronPDF> hard_pomeron_pdf = nullptr;
  std::vector<HardPartonChannel>         hard_parton_channels;
  std::vector<std::size_t>               hard_channel_component_counts;
  std::vector<std::size_t>               hard_component_channel_indices;
  // Born leading-color bank used to validate every screening node
  std::vector<mg5helas::HardColorFlow> hard_color_flows;
  // Incoming hard channel selected for the accepted event
  std::size_t hard_selected_channel_index = 0;
  double      hard_integral_volume        = 0.0;
  // Isotropic rest-frame decay coordinates shared by channels and screening nodes
  std::array<double, 2> remnant_angles = {0.5, 0.0};
  // Central production coordinates for adaptive N-body integration
  std::vector<double> central_coordinates;
};

}  // namespace gra

#endif
