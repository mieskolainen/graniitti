// Process kinematics unit tests
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include <algorithm>
#include <array>
#include <catch.hpp>
#include <cmath>
#include <complex>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/Math/MFloat.h"
#include "Graniitti/Math/MIntegration.h"
#include "Graniitti/Particle/MForm.h"
#include "Graniitti/Nuclear/MBeam.h"
#include "Graniitti/Kinematics/MCentral.h"
#include "Graniitti/Kinematics/MCollinear.h"
#include "Graniitti/Kinematics/MFactorized.h"
#include "Graniitti/Kinematics/MHardDiffraction.h"
#include "Graniitti/Process/MProcess.h"
#include "Graniitti/Process/MDecayTree.h"
#include "Graniitti/MGraniitti.h"
#include "Graniitti/Process/MForwardExcitation.h"
#include "Graniitti/Kinematics/MEventKinematics.h"
#include "Graniitti/Kinematics/MQuasiElastic.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"

namespace {

using gra::aux::indices;

// Expose protected process kinematics without selecting a physics amplitude
class ProcessKinematicsProbe : public gra::MFactorized {
public:
  // Initialize deterministic storage and random numbers for unit tests
  ProcessKinematicsProbe() {
    state.lts.pfinal.assign(11, gra::M4Vec());
    state.lts.beam1.pdg = state.lts.beam2.pdg = gra::PDG::PDG_p;
    state.random.SetSeed(13579);
  }

  // Configure the common forward-mass integration domain
  void ConfigureForward(int excitation, double s, double xi_min, double xi_max,
                        double pt_min, double pt_max) {
    state.excitation = excitation;
    state.lts.s = s;
    state.lts.beam1.mass = gra::PDG::mp;
    state.lts.beam2.mass = gra::PDG::mp;
    state.gcuts.XI_min = xi_min;
    state.gcuts.XI_max = xi_max;
    state.gcuts.forward_pt_min = pt_min;
    state.gcuts.forward_pt_max = pt_max;
  }

  // Sample the forward masses through the production implementation
  std::vector<double> SampleMasses(const std::vector<double> &randoms) {
    std::vector<double> masses;
    SampleForwardMasses(masses, randoms);
    return masses;
  }

  // Install sampled forward masses and fixed transverse momenta
  void SetForwardPoint(const std::vector<double> &masses, double pt1,
                       double pt2) {
    state.lts.pfinal[1].SetPxPyPzM(pt1, 0.0, 0.0, masses[0]);
    state.lts.pfinal[2].SetPxPyPzM(0.0, pt2, 0.0, masses[1]);
  }

  // Compute the production forward phase-space volume
  double ProbeForwardVolume() const { return ForwardVolume(); }

  // Configure one process topology for dissociation normalization tests
  void ConfigureDissociationTopology(int excitation,
                                     const std::string &initial_state,
                                     const std::string &channel) {
    state.excitation = excitation;
    ProcPtr.ISTATE = initial_state;
    ProcPtr.CHANNEL = channel;
  }

  // Compute the common dissociation cross-section multiplicity
  double ProbeDissociationCrossSectionFactor() const {
    return DissociationCrossSectionFactor();
  }

  // Compute the configured invariant-mass-squared integration bounds
  std::pair<double, double> ForwardMassBounds() const {
    return {M2_f_min, M2_f_max};
  }

  // Select whether generated cascade phase space belongs to the integrand
  void SetActiveDecayPhaseSpace(bool active) { state.lts.PS_active = active; }

  // Select the sampled resonance-width window
  void SetOffShellRange(double offshell) { state.offshell_widths = offshell; }

  // Replace the top-level decay tree
  void SetDecayTree(const gra::MDecayBranch &root) {
    state.lts.decaytree = {root};
  }

  // Generate all kinematics below the only top-level branch
  bool GenerateDecayTree() {
    return ConstructDecayKinematics(state.lts.decaytree.at(0));
  }

  // Compute the current event cascade factor
  double ProbeCascadePS() { return CascadePS(); }

  // Compute the direct central phase-space volume for the current roots
  double ProbeCentralDecayWidthPS() const { return CentralDecayWidthPS(); }

  // Prepare the production-root mass threshold through the common helper
  bool ProbePrepareCentralBranchMasses(double &mass_sum,
                                       double *mass_max = nullptr) {
    return PrepareCentralBranchMasses(mass_sum, mass_max);
  }

  // Bind a sole production root through the common kinematic helper
  bool ProbeSetSingleCentralRootKinematics() {
    return SetSingleCentralRootKinematics();
  }

  // Exercise scalar-storage range validation
  bool ProbeLorentzScalars(unsigned int final_count) {
    return gra::kinematics::SetLorentzScalars(state, final_count);
  }

  // Build one forward branch through the transactional process helper
  bool ProbeForwardBranch(const std::vector<gra::M4Vec> &p4,
                          const std::vector<gra::MParticle> &particles,
                          const gra::M4Vec &parent, gra::MDecayBranch &branch) {
    return BranchForwardSystem(p4, particles, parent, branch);
  }

  // Install veto domains for direct source and charge-selection tests
  void ConfigureVetoCuts(const gra::VETOCUT &cuts) { SetVetoCuts(cuts); }

  // Apply configured veto cuts to explicit forward and central branches
  bool ProbeVetoCuts(const gra::MDecayBranch &forward,
                     const gra::MDecayBranch &central) {
    state.lts.decayforward1 = forward;
    state.lts.decayforward2 = gra::MDecayBranch();
    state.lts.decaytree = {central};
    return state.vetocuts.Pass(state.lts);
  }

  using gra::MProcess::BookkeepAmplitudeWeight;
  using gra::MProcess::ExciteString;
};

// Expose the central phase space transverse map helpers for unit tests
class CentralKinematicsProbe : public gra::MCentral {
public:
  // Select a physical process owner for amplitude normalization
  CentralKinematicsProbe() { ProcPtr.Initialize("MP", "RES"); }

  // Compute the central mass ceiling for one forward point
  static double MassCeiling(double sqrt_s, double mass1, double mass2,
                            const gra::M4Vec &forward1,
                            const gra::M4Vec &forward2) {
    return CentralMassKinematicMaximum(sqrt_s, mass1, mass2, forward1,
                                       forward2);
  }

  // Compute the squared transverse Helmert ball radius
  static double TransverseRadius2(double mass_max, double recoil_pt2,
                                  unsigned int central_multiplicity) {
    return MassConditionedTransverseRadius2(mass_max, recoil_pt2,
                                            central_multiplicity);
  }

  // Map one logarithmic Helmert hyperradius and return its Jacobian
  static std::pair<std::vector<gra::M4Vec>, double>
  MapTransverseLogarithmic(const std::vector<double> &radial_units,
                           const std::vector<double> &angle_units,
                           const gra::M4Vec &central_transverse_momentum,
                           double radius2, double scale2) {
    double jacobian = 0.0;
    const auto differences =
        MapTransverseLog(radial_units, angle_units, central_transverse_momentum,
                         radius2, scale2, jacobian);
    return {differences, jacobian};
  }

  // Build one deterministic central phase space point
  bool BuildBornForTest(unsigned int final_count, double pt1, double pt2,
                        double phi1, double phi2, const std::vector<double> &kt,
                        const std::vector<double> &phi,
                        const std::vector<double> &rapidity, double mass1,
                        double mass2) {
    state.gcuts.M_min = 0.0;
    state.gcuts.M_max = state.lts.sqrt_s;
    return BNBuildKin(final_count, pt1, pt2, phi1, phi2, kt, phi, rapidity,
                      mass1, mass2);
  }

  // Rebuild one central phase space screening point
  bool BuildLoopForTest(const std::array<double, 2> &p1t,
                        const std::array<double, 2> &p2t) {
    return LoopKinematics(p1t, p2t);
  }

  // Restore the saved central Born point through the production path
  bool RestoreBornForTest() { return RestoreBornKinematics(); }

};

// Expose the factorized central mass map without selecting an amplitude
class FactorizedKinematicsProbe : public gra::MFactorized {
public:
  // Select a physical process owner for amplitude normalization
  FactorizedKinematicsProbe() { ProcPtr.Initialize("MP", "RES"); }

  // Sample one central mass squared point through the production map
  static double SampleCentralMassSquaredForTest(double unit, double mass_min,
                                                double mass_max) {
    return SampleCentralMassSquared(unit, mass_min, mass_max);
  }

  // Compute the production Jacobian for one sampled central mass squared
  static double CentralMassSquaredJacobianForTest(double mass_squared,
                                                  double mass_min,
                                                  double mass_max) {
    return CentralMassSquaredJacobian(mass_squared, mass_min, mass_max);
  }

  // Build one deterministic factorized point
  bool BuildBornForTest(double pt1, double pt2, double phi1, double phi2,
                        double rapidity, double mass_squared, double mass1,
                        double mass2) {
    return B51BuildKin(pt1, pt2, phi1, phi2, rapidity, mass_squared, mass1,
                       mass2);
  }

  // Rebuild one factorized screening point
  bool BuildLoopForTest(const std::array<double, 2> &p1t,
                        const std::array<double, 2> &p2t) {
    return LoopKinematics(p1t, p2t);
  }

  // Restore the saved factorized Born point through the production path
  bool RestoreBornForTest() { return RestoreBornKinematics(); }

};

// Expose the collinear central mass map without selecting an amplitude
class CollinearKinematicsProbe : public gra::MCollinear {
public:
  // Exercise the complete collinear random map through its production API
  using gra::MCollinear::ComputeEventWeight;

  // Initialize the stable final-state normalization for the sampled tree
  using gra::MProcess::CalculateSymmetryFactor;

  // Select a physical process owner for amplitude normalization
  CollinearKinematicsProbe() { ProcPtr.Initialize("MP", "CON"); }

  // Sample one central mass squared point through the production map
  static double SampleCentralMassSquaredForTest(double unit, double mass_min,
                                                double mass_max) {
    return SampleCentralMassSquared(unit, mass_min, mass_max);
  }

  // Compute the production Jacobian for one sampled central mass squared
  static double CentralMassSquaredJacobianForTest(double mass_squared,
                                                  double mass_min,
                                                  double mass_max) {
    return CentralMassSquaredJacobian(mass_squared, mass_min, mass_max);
  }

  // Build one deterministic collinear parton point
  bool BuildBornForTest(double x1, double x2) { return B2BuildKin(x1, x2); }

  // Compute the unsupported parton screening-loop decision
  bool BuildLoopForTest(const std::array<double, 2> &p1t,
                        const std::array<double, 2> &p2t) {
    return LoopKinematics(p1t, p2t);
  }

};

// Expose quasielastic Born and screening kinematics for unit tests
class QuasiElasticKinematicsProbe : public gra::MQuasiElastic {
public:
  using gra::MQuasiElastic::BuildEventRecord;
  // Select a physical process owner for amplitude normalization
  QuasiElasticKinematicsProbe() { ProcPtr.Initialize("MP", "RES"); }

  // Seed the quasielastic azimuth generator
  void SeedForTest(std::uint64_t seed) { state.random.SetSeed(seed); }

  // Build one deterministic quasielastic point
  bool BuildBornForTest(double mass_squared1, double mass_squared2, double t) {
    return B3BuildKin(mass_squared1, mass_squared2, t);
  }

  // Rebuild one quasielastic screening point
  bool BuildLoopForTest(const std::array<double, 2> &p1t,
                        const std::array<double, 2> &p2t) {
    return LoopKinematics(p1t, p2t);
  }

  // Restore the saved quasielastic Born point through the production path
  bool RestoreBornForTest() { return RestoreBornKinematics(); }

  // Configure the quasielastic absolute momentum-transfer maximum
  void SetMaximumSampledAbsTForTest(double configured_abs_t_max,
                                    double loop_max_kt) {
    state.gcuts.q_t_abs_max = configured_abs_t_max;
    eikonal.Numerics.LOOP.r_max = loop_max_kt;
  }

  // Configure one X<Q> channel and its generated momentum-transfer range
  void SetQuasiElasticRangeForTest(const std::string &channel, double min_abs_t,
                                   double max_abs_t) {
    ProcPtr.Initialize("X", channel);
    state.gcuts.q_t_abs_min = min_abs_t;
    state.gcuts.q_t_abs_max = max_abs_t;
  }

  // Compute the active quasielastic absolute momentum-transfer maximum
  double MaximumSampledAbsTForTest() const { return MaximumSampledAbsT(); }

  // Compute the process requirement for the eikonal momentum table
  double EikonalMaxKT2ForTest() const { return EikonalMaxKT2(); }

  // Sample elastic absolute momentum transfer through the production map
  static double SampleElasticAbsTForTest(double unit, double min_abs_t,
                                         double max_abs_t) {
    return SampleElasticAbsT(unit, min_abs_t, max_abs_t);
  }

  // Compute the production elastic mixed-proposal inverse density
  static double ElasticAbsTJacobianForTest(double abs_t, double min_abs_t,
                                           double max_abs_t) {
    return ElasticAbsTJacobian(abs_t, min_abs_t, max_abs_t);
  }

};

TEST_CASE("quasielastic t sampling is independent of the loop momentum range",
          "[Process][PhaseSpace][QuasiElastic]") {
  QuasiElasticKinematicsProbe probe;

  probe.SetQuasiElasticRangeForTest("EL", 1.0e-4, 0.8);
  probe.SetMaximumSampledAbsTForTest(0.8, 2.25);
  REQUIRE(probe.MaximumSampledAbsTForTest() == Approx(0.8));

  probe.SetMaximumSampledAbsTForTest(12.0, 0.25);
  REQUIRE(probe.MaximumSampledAbsTForTest() == Approx(12.0));
}

TEST_CASE("only elastic X<Q> extends the eikonal momentum table",
          "[Process][PhaseSpace][QuasiElastic][Eikonal]") {
  QuasiElasticKinematicsProbe probe;

  probe.SetQuasiElasticRangeForTest("EL", 1.0e-4, 12.0);
  REQUIRE_NOTHROW(probe.FinalizeProcessConfiguration());
  CHECK(probe.EikonalMaxKT2ForTest() == Approx(12.0));

  probe.SetQuasiElasticRangeForTest("SD", 0.0, 50.0);
  REQUIRE_NOTHROW(probe.FinalizeProcessConfiguration());
  CHECK(probe.EikonalMaxKT2ForTest() == Approx(0.0));

  probe.SetQuasiElasticRangeForTest("DD", 0.0, 50.0);
  REQUIRE_NOTHROW(probe.FinalizeProcessConfiguration());
  CHECK(probe.EikonalMaxKT2ForTest() == Approx(0.0));
}

TEST_CASE("non-diffractive soft sampling rejects an unavailable eikonal",
          "[Process][QuasiElastic][Failure]") {
  gra::MQuasiElastic process;
  process.ProcPtr.ISTATE = "X";
  process.ProcPtr.CHANNEL = "ND";
  process.ProcPtr.LIPSDIM = 1;

  gra::MEventWeightState status;
  double weight = -1.0;
  CHECK_NOTHROW(weight = process.EventWeight({0.5}, status));
  CHECK(gra::math::IsZero(weight));
  CHECK_FALSE(status.kinematics_ok);
  CHECK(status.technical_failure);
  CHECK_FALSE(status.Valid());
  CHECK(process.state.multipomeron_chain.empty());
}

TEST_CASE("elastic t sampling mixes inverse-variable and logarithmic maps",
          "[Process][PhaseSpace][QuasiElastic][Elastic]") {
  constexpr double min_abs_t = 0.0006;
  constexpr double max_abs_t = 5.0;

  REQUIRE(QuasiElasticKinematicsProbe::SampleElasticAbsTForTest(
              0.0, min_abs_t, max_abs_t) == Approx(min_abs_t));
  REQUIRE(QuasiElasticKinematicsProbe::SampleElasticAbsTForTest(
              0.5, min_abs_t, max_abs_t) == Approx(min_abs_t));
  REQUIRE(QuasiElasticKinematicsProbe::SampleElasticAbsTForTest(
              1.0, min_abs_t, max_abs_t) == Approx(max_abs_t));

  const double inverse_mid =
      QuasiElasticKinematicsProbe::SampleElasticAbsTForTest(0.25, min_abs_t,
                                                            max_abs_t);
  const double logarithmic_mid =
      QuasiElasticKinematicsProbe::SampleElasticAbsTForTest(0.75, min_abs_t,
                                                            max_abs_t);
  REQUIRE(inverse_mid ==
          Approx(2.0 * min_abs_t * max_abs_t / (min_abs_t + max_abs_t)));
  REQUIRE(logarithmic_mid == Approx(std::sqrt(min_abs_t * max_abs_t)));

  constexpr std::size_t samples = 20000;
  double integral = 0.0;
  for (std::size_t i = 0; i < samples; ++i) {
    const double unit =
        (static_cast<double>(i) + 0.5) / static_cast<double>(samples);
    const double abs_t = QuasiElasticKinematicsProbe::SampleElasticAbsTForTest(
        unit, min_abs_t, max_abs_t);
    REQUIRE(abs_t >= min_abs_t);
    REQUIRE(abs_t <= max_abs_t);
    integral += QuasiElasticKinematicsProbe::ElasticAbsTJacobianForTest(
        abs_t, min_abs_t, max_abs_t);
  }
  integral /= static_cast<double>(samples);
  REQUIRE(integral == Approx(max_abs_t - min_abs_t).epsilon(2.0e-4));
}

// Expose hard diffraction Born and screening kinematics for unit tests
class HardDiffractionKinematicsProbe : public gra::MHardDiffraction {
public:
  // Select a physical process owner for amplitude normalization
  HardDiffractionKinematicsProbe() {
    ProcPtr.Initialize("MP", "RES");
    state.lts.PDG.ReadParticleData();
    SetModelTune(gra::MModelTune::Load(gra::ResolveModelDataFile("TUNE0", "GENERAL.json")));
    FinalizeProcessConfiguration();
  }

  // Build one deterministic hard diffraction point
  bool BuildBornForTest(double x1, double x2) {
    state.gcuts.M_min = 0.0;
    state.gcuts.M_max = state.lts.sqrt_s;
    state.gcuts.Y_min = -30.0;
    state.gcuts.Y_max = 30.0;
    state.lts.diff_xhard1 = x1;
    state.lts.diff_xhard2 = x2;
    return HardBuildKin(x1, x2);
  }

  // Rebuild one hard diffraction screening point
  bool BuildLoopForTest(const std::array<double, 2> &p1t,
                        const std::array<double, 2> &p2t) {
    return LoopKinematics(p1t, p2t);
  }

  // Restore the saved hard-diffraction Born point through the production path
  bool RestoreBornForTest() { return RestoreBornKinematics(); }

  // Apply the shared forward xi cut to the prepared hard diffraction point
  bool PassForwardXiForTest(const gra::FIDCUT &cuts) {
    state.fcuts = cuts;
    return state.fcuts.PassForwardXi(state.lts);
  }

};

// Initialize the actual EPA dimuon process and its normal phase-space configuration
std::unique_ptr<gra::MGraniitti> DimuonGenerator(bool fsr = false) {
  auto card = nlohmann::json::parse(gra::aux::GetInputData(gra::aux::ResolveProjectPath("gencard/test.json")));
  card["GENERIC"]["HIST"] = 0;
  card["GENERIC"]["CORES"] = 1;
  card["GENERIC"]["RNDSEED"] = 13579;
  card["GENERIC"]["INTEGRATOR"] = "VEGAS";
  card["SCATTERING"]["PROCESS"] = "yy[EPA]<F> -> mu+ mu-";
  card["SCATTERING"]["ENERGY"] = {80.0, 20.0};
  card["SCATTERING"]["LOOPSCREEN"] = false;
  card["GENCUTS"]["<F>"]["M"] = {19.9, 20.1};
  card["GENCUTS"]["<F>"]["Rap"] = {-1.0, 1.0};
  card["GENCUTS"]["<F>"]["Pt"] = {0.1, 0.6};
  card["FIDCUTS"] = {{"active", false}};
  auto generator = std::make_unique<gra::MGraniitti>();
  generator->ReadInput(card);
  if (fsr) {
    gra::radiative::Config config;
    config.fsr = gra::radiative::Mode::YFS;
    config.fsr_param.energy_min = 1e-4;
    generator->proc->SetRadiative(config);
  }
  generator->proc->PrepareRun();
  REQUIRE(generator->proc->GetdLIPSDim() > 0);
  return generator;
}

// Find a positive event using the real process and independent unit coordinates
std::vector<double> DimuonPoint(gra::MProcess &process, gra::MRandom &random) {
  std::vector<double> point(process.GetdLIPSDim());
  for (std::size_t trial = 0; trial < 1000; ++trial) {
    for (auto &unit : point) { unit = random.U(0.0, 1.0); }
    gra::MEventWeightState aux;
    const double weight = process.EventWeight(point, aux);
    REQUIRE_FALSE(aux.technical_failure);
    if (weight > 0.0) {
      REQUIRE(aux.Valid());
      return point;
    }
  }
  FAIL("EPA dimuon sampling found no physical support");
  return {};
}

// Construct one stable massless decay leaf
gra::MDecayBranch MasslessLeaf() {
  gra::MDecayBranch leaf;
  leaf.p.mass = 0.0;
  leaf.p.width = 0.0;
  return leaf;
}

// Construct one stable particle with fixed charge and cylindrical momentum
gra::MDecayBranch StableParticle(int charge_x3, double pt, double eta) {
  gra::MDecayBranch branch;
  branch.p.chargeX3 = charge_x3;
  branch.p4.SetPxPyPzM(pt, 0.0, pt * std::sinh(eta), 0.0);
  return branch;
}

// Construct one direct on-shell muon pair at rest as a system
std::vector<gra::MDecayBranch> MuonPair(double mass) {
  constexpr double muon_mass = 0.1056583745;
  const double momentum =
      std::sqrt(mass * mass / 4.0 - muon_mass * muon_mass);
  gra::MDecayBranch minus;
  minus.p.pdg = 13;
  minus.p.mass = muon_mass;
  minus.p4 = gra::M4Vec(0.0, 0.0, momentum, mass / 2.0);
  gra::MDecayBranch plus;
  plus.p.pdg = -13;
  plus.p.mass = muon_mass;
  plus.p4 = gra::M4Vec(0.0, 0.0, -momentum, mass / 2.0);
  return {minus, plus};
}

// Construct one fixed-mass parent with the requested daughter branches
gra::MDecayBranch FixedDecay(double mass,
                             const std::vector<gra::MDecayBranch> &daughters) {
  gra::MDecayBranch branch;
  branch.p.mass = mass;
  branch.p.width = 0.0;
  branch.p4 = gra::M4Vec(0.0, 0.0, 0.0, mass);
  branch.legs = daughters;
  return branch;
}

// Configure physical proton beams with unequal lab energies
template <typename Probe>
void ConfigureAsymmetricProtonBeams(Probe &probe, double energy1,
                                    double energy2) {
  probe.state.lts.PDG.ReadParticleData();
  if (probe.state.model_tune == nullptr) {
    probe.SetModelTune(gra::MModelTune::Load(gra::ResolveModelDataFile("TUNE0", "GENERAL.json")));
  }
  probe.SetInitialState({"p+", "p+"}, {energy1, energy2});
  probe.state.lts.pfinal.assign(11, gra::M4Vec());
}

// Require two four-vectors to agree component by component
void RequireFourVectorNear(const gra::M4Vec &actual, const gra::M4Vec &expected,
                           double tolerance = 1.0e-10) {
  REQUIRE(actual.Px() == Approx(expected.Px()).margin(tolerance));
  REQUIRE(actual.Py() == Approx(expected.Py()).margin(tolerance));
  REQUIRE(actual.Pz() == Approx(expected.Pz()).margin(tolerance));
  REQUIRE(actual.E() == Approx(expected.E()).margin(tolerance));
}

// Require one restored four-vector to retain its exact saved representation
void RequireFourVectorExact(const gra::M4Vec &actual,
                            const gra::M4Vec &expected) {
  REQUIRE(gra::math::IsExactEqual(actual.Px(), expected.Px()));
  REQUIRE(gra::math::IsExactEqual(actual.Py(), expected.Py()));
  REQUIRE(gra::math::IsExactEqual(actual.Pz(), expected.Pz()));
  REQUIRE(gra::math::IsExactEqual(actual.E(), expected.E()));
}

// Construct one unbound finite-width root with a two-body stable decay
gra::MDecayBranch SingleCascadeRoot(double pole_mass = 2.0,
                                    double width = 0.2) {
  gra::MDecayBranch root;
  root.p.pdg = 23;
  root.p.mass = pole_mass;
  root.p.width = width;
  root.m_offshell = 1.7;
  root.mass_proposal_norm = 0.4;
  root.mass_proposal_min2 = 1.0;
  root.mass_proposal_max2 = 4.0;
  root.mass_proposal = gra::MassProposal::BreitWigner;
  root.legs = {MasslessLeaf(), MasslessLeaf()};
  return root;
}

// Require one production root to own the full central-system momentum
void RequireSingleRootKinematics(const gra::LORENTZSCALAR &lts) {
  REQUIRE(lts.decaytree.size() == 1);
  const auto &root = lts.decaytree.front();
  RequireFourVectorNear(root.p4, lts.pfinal[0], 1.0e-10);
  RequireFourVectorNear(lts.pfinal[3], lts.pfinal[0], 1.0e-10);
  REQUIRE(root.m_offshell == Approx(lts.pfinal[0].M()).epsilon(1.0e-12));
  REQUIRE(root.mass_proposal_norm == Approx(1.0));
  REQUIRE(root.mass_proposal_min2 == Approx(0.0));
  REQUIRE(root.mass_proposal_max2 == Approx(0.0));
  REQUIRE_FALSE((root.mass_proposal == gra::MassProposal::BreitWigner));
  REQUIRE_FALSE((root.mass_proposal == gra::MassProposal::Uniform));
  REQUIRE_FALSE((root.mass_proposal == gra::MassProposal::Fixed));
  REQUIRE(lts.DW.GetN() == Approx(1.0));
  REQUIRE(lts.DW.Integral() == Approx(1.0));
  REQUIRE(root.W_event > 0.0);
  REQUIRE(root.legs.size() == 2);
  RequireFourVectorNear(root.legs[0].p4 + root.legs[1].p4, root.p4, 1.0e-8);
}

// Compute one lab four-vector in the collision CM
gra::M4Vec CollisionCMVector(const gra::LORENTZSCALAR &lts,
                             const gra::M4Vec &vector_lab) {
  const gra::M4Vec beamsum = lts.pbeam1 + lts.pbeam2;
  gra::M4Vec vector_cm = vector_lab;
  gra::kinematics::LorentzBoost(beamsum, beamsum.M(), vector_cm, -1);
  return vector_cm;
}

// Require lab-frame energy-momentum conservation for one phase-space skeleton
void RequirePhaseSpaceClosure(const gra::LORENTZSCALAR &lts,
                              bool central_system_active) {
  const gra::M4Vec central =
      central_system_active ? lts.pfinal[0] : gra::M4Vec();
  REQUIRE(gra::math::CheckEMC(lts.pbeam1 + lts.pbeam2 - lts.pfinal[1] -
                              lts.pfinal[2] - central));
}

// Require one forward branch to reproduce its parent momentum
void RequireForwardBranchClosure(const gra::MDecayBranch &branch,
                                 const gra::M4Vec &parent,
                                 std::size_t expected_legs) {
  REQUIRE(branch.legs.size() == expected_legs);
  gra::M4Vec sum;
  for (const auto &leg : branch.legs) {
    sum += leg.p4;
  }
  RequireFourVectorNear(sum, parent, 1.0e-8);
}

// Require every central decay momentum to remain fixed through screening
void RequireDecayTreeMomentaNear(const std::vector<gra::MDecayBranch> &actual,
                                 const std::vector<gra::MDecayBranch> &born) {
  REQUIRE(actual.size() == born.size());
  for (const auto &i : gra::aux::indices(born)) {
    RequireFourVectorNear(actual[i].p4, born[i].p4);
    RequireDecayTreeMomentaNear(actual[i].legs, born[i].legs);
  }
}

// Exercise bare channel selection and screening-loop kinematics at one
// asymmetric point
template <typename Probe>
void RequireAsymmetricScreeningKinematics(Probe &probe,
                                          bool central_system_active,
                                          bool require_hard_remnants = false) {
  const std::vector<gra::M4Vec> born = probe.state.lts.pfinal;
  const std::vector<gra::MDecayBranch> born_decay = probe.state.lts.decaytree;
  const std::array<gra::MDecayBranch, 2> born_forward = {
      probe.state.lts.decayforward1, probe.state.lts.decayforward2};
  REQUIRE(probe.state.lts.pbeam1.E() !=
          Approx(probe.state.lts.pbeam2.E()).epsilon(1.0e-12));

  probe.SetScreening(false);
  REQUIRE_FALSE(probe.GetScreening());
  RequirePhaseSpaceClosure(probe.state.lts, central_system_active);
  for (std::size_t i = 0; i < 3; ++i) {
    RequireFourVectorNear(probe.state.lts.pfinal[i], born[i]);
  }

  probe.SetScreening(true);
  REQUIRE(probe.GetScreening());
  for (std::size_t i = 0; i < 3; ++i) {
    RequireFourVectorNear(probe.state.lts.pfinal[i], born[i]);
  }
  probe.state.lts.pfinal_orig = born;
  probe.state.lts.screening.active = true;

  constexpr double shift_x = 0.04;
  constexpr double shift_y = -0.03;
  const std::array<double, 2> shifted1 = {born[1].Px() - shift_x,
                                          born[1].Py() - shift_y};
  const std::array<double, 2> shifted2 = {born[2].Px() + shift_x,
                                          born[2].Py() + shift_y};
  REQUIRE(probe.BuildLoopForTest(shifted1, shifted2));
  RequirePhaseSpaceClosure(probe.state.lts, central_system_active);
  REQUIRE(probe.state.lts.pfinal[1].M2() ==
          Approx(born[1].M2()).margin(1.0e-8));
  REQUIRE(probe.state.lts.pfinal[2].M2() ==
          Approx(born[2].M2()).margin(1.0e-8));
  if (central_system_active) {
    RequireFourVectorNear(probe.state.lts.pfinal[0], born[0]);
    RequireDecayTreeMomentaNear(probe.state.lts.decaytree, born_decay);
  }
  if (require_hard_remnants) {
    RequireForwardBranchClosure(probe.state.lts.decayforward1,
                                probe.state.lts.pfinal[1], 2);
    RequireForwardBranchClosure(probe.state.lts.decayforward2,
                                probe.state.lts.pfinal[2], 2);
    REQUIRE(probe.state.lts.decayforward1.legs[1].p4.E() > 0.0);
    REQUIRE(probe.state.lts.decayforward2.legs[1].p4.E() > 0.0);
    REQUIRE(probe.state.lts.decayforward1.legs[1].p4.M2() ==
            Approx(born_forward[0].legs[1].p4.M2()).margin(1.0e-8));
    REQUIRE(probe.state.lts.decayforward2.legs[1].p4.M2() ==
            Approx(born_forward[1].legs[1].p4.M2()).margin(1.0e-8));
  }

  const gra::M4Vec beamsum = probe.state.lts.pbeam1 + probe.state.lts.pbeam2;
  gra::M4Vec forward1_cm = probe.state.lts.pfinal[1];
  gra::M4Vec forward2_cm = probe.state.lts.pfinal[2];
  gra::kinematics::LorentzBoost(beamsum, beamsum.M(), forward1_cm, -1);
  gra::kinematics::LorentzBoost(beamsum, beamsum.M(), forward2_cm, -1);
  REQUIRE(forward1_cm.Pz() > 0.0);
  REQUIRE(forward2_cm.Pz() < 0.0);

  probe.state.lts.screening.active = false;
  REQUIRE(probe.RestoreBornForTest());
  for (std::size_t i = 0; i < 3; ++i) {
    RequireFourVectorExact(probe.state.lts.pfinal[i], born[i]);
  }
  if (require_hard_remnants) {
    RequireForwardBranchClosure(probe.state.lts.decayforward1,
                                probe.state.lts.pfinal[1], 2);
    RequireForwardBranchClosure(probe.state.lts.decayforward2,
                                probe.state.lts.pfinal[2], 2);
    REQUIRE(probe.state.lts.decayforward1.legs[1].p4.M2() ==
            Approx(born_forward[0].legs[1].p4.M2()).margin(1.0e-8));
    REQUIRE(probe.state.lts.decayforward2.legs[1].p4.M2() ==
            Approx(born_forward[1].legs[1].p4.M2()).margin(1.0e-8));
  }
}

} // namespace

TEST_CASE("phase-space Born and supported screening-loop kinematics "
          "support asymmetric beam energies",
          "[Process][PhaseSpace][AsymmetricBeams][screening]") {
  const auto beam_energies =
      GENERATE(std::make_pair(80.0, 20.0), std::make_pair(20.0, 80.0));
  CAPTURE(beam_energies.first, beam_energies.second);

  SECTION("central C") {
    CentralKinematicsProbe probe;
    ConfigureAsymmetricProtonBeams(probe, beam_energies.first,
                                   beam_energies.second);
    gra::MDecayBranch left = MasslessLeaf();
    gra::MDecayBranch right = MasslessLeaf();
    left.p.mass = 0.3;
    right.p.mass = 0.4;
    left.m_offshell = left.p.mass;
    right.m_offshell = right.p.mass;
    probe.state.lts.decaytree = {left, right};
    probe.state.lts.central_phase_space_mode =
        gra::CentralPhaseSpaceMode::Central;

    REQUIRE(probe.BuildBornForTest(4, 0.25, 0.35, 0.2, 2.4, {0.45}, {1.1},
                                   {0.15, -0.25}, gra::PDG::mp, gra::PDG::mp));
    REQUIRE(CollisionCMVector(probe.state.lts, probe.state.lts.decaytree[0].p4)
                .Rap() == Approx(0.15).margin(1.0e-11));
    REQUIRE(CollisionCMVector(probe.state.lts, probe.state.lts.decaytree[1].p4)
                .Rap() == Approx(-0.25).margin(1.0e-11));
    RequireAsymmetricScreeningKinematics(probe, true);
  }

  SECTION("factorized F") {
    FactorizedKinematicsProbe probe;
    ConfigureAsymmetricProtonBeams(probe, beam_energies.first,
                                   beam_energies.second);
    probe.state.lts.decaytree.clear();

    REQUIRE(probe.BuildBornForTest(0.25, 0.35, 0.2, 2.4, 0.1, 9.0, gra::PDG::mp,
                                   gra::PDG::mp));
    REQUIRE(
        CollisionCMVector(probe.state.lts, probe.state.lts.pfinal[0]).Rap() ==
        Approx(0.1).margin(1.0e-11));
    RequireAsymmetricScreeningKinematics(probe, true);
  }

  SECTION("quasielastic Q") {
    QuasiElasticKinematicsProbe probe;
    ConfigureAsymmetricProtonBeams(probe, beam_energies.first,
                                   beam_energies.second);
    probe.SeedForTest(13579);

    REQUIRE(probe.BuildBornForTest(gra::math::pow2(gra::PDG::mp),
                                   gra::math::pow2(1.3), -0.25));
    RequireAsymmetricScreeningKinematics(probe, false);
    REQUIRE_FALSE(probe.BuildBornForTest(
        std::numeric_limits<double>::quiet_NaN(), gra::math::pow2(1.3), -0.25));
  }

  SECTION("hard diffraction D") {
    HardDiffractionKinematicsProbe probe;
    ConfigureAsymmetricProtonBeams(probe, beam_energies.first,
                                   beam_energies.second);
    probe.state.lts.decaytree.clear();
    probe.state.lts.hard_diff1 = true;
    probe.state.lts.hard_diff2 = true;
    probe.state.lts.diff_xi1 = 0.12;
    probe.state.lts.diff_xi2 = 0.10;
    probe.state.lts.diff_beta1 = 0.65;
    probe.state.lts.diff_beta2 = 0.55;
    probe.state.lts.diff_t1 = -0.16;
    probe.state.lts.diff_t2 = -0.09;
    probe.state.lts.diff_phi1 = 0.3;
    probe.state.lts.diff_phi2 = 0.3 + gra::math::PI;
    probe.state.lts.id1 = gra::PDG::PDG_gluon;
    probe.state.lts.id2 = gra::PDG::PDG_gluon;

    REQUIRE(probe.BuildBornForTest(
        probe.state.lts.diff_xi1 * probe.state.lts.diff_beta1,
        probe.state.lts.diff_xi2 * probe.state.lts.diff_beta2));
    REQUIRE(probe.state.lts.pfinal[0].M2() ==
            Approx(probe.state.lts.s * probe.state.lts.diff_xhard1 *
                   probe.state.lts.diff_xhard2)
                .epsilon(1.0e-11));
    REQUIRE(probe.state.lts.x1 > probe.state.lts.diff_xhard1);
    REQUIRE(probe.state.lts.x2 > probe.state.lts.diff_xhard2);
    REQUIRE(probe.state.lts.x1 < probe.state.lts.diff_xi1);
    REQUIRE(probe.state.lts.x2 < probe.state.lts.diff_xi2);
    REQUIRE(probe.state.lts.xi1 ==
            Approx(probe.state.lts.diff_xi1).margin(1.0e-12));
    REQUIRE(probe.state.lts.xi2 ==
            Approx(probe.state.lts.diff_xi2).margin(1.0e-12));
    REQUIRE(probe.state.lts.has_xi1);
    REQUIRE(probe.state.lts.has_xi2);
    REQUIRE(probe.state.lts.x1 != Approx(probe.state.lts.xi1));
    REQUIRE(probe.state.lts.x2 != Approx(probe.state.lts.xi2));
    gra::FIDCUT xi_cut;
    xi_cut.forward_xi_active = true;
    xi_cut.forward_xi_min = 0.09;
    xi_cut.forward_xi_max = 0.13;
    REQUIRE(probe.PassForwardXiForTest(xi_cut));
    REQUIRE(
        std::abs(CollisionCMVector(probe.state.lts, probe.state.lts.pfinal[0])
                     .Rap()) < 0.25);
    RequireAsymmetricScreeningKinematics(probe, true, true);
    REQUIRE(probe.state.lts.pfinal[0].M2() ==
            Approx(probe.state.lts.s * probe.state.lts.diff_xhard1 *
                   probe.state.lts.diff_xhard2)
                .epsilon(1.0e-11));
  }

  SECTION("collinear parton P") {
    CollinearKinematicsProbe probe;
    ConfigureAsymmetricProtonBeams(probe, beam_energies.first,
                                   beam_energies.second);
    gra::MDecayBranch left = MasslessLeaf();
    gra::MDecayBranch right = MasslessLeaf();
    left.m_offshell = 0.0;
    right.m_offshell = 0.0;
    probe.state.lts.decaytree = {left, right};

    const bool parton_born = probe.BuildBornForTest(0.20, 0.15);
    CAPTURE(probe.state.lts.pfinal[0].M2(), probe.state.lts.pfinal[1].M2(),
            probe.state.lts.pfinal[2].M2(), probe.state.lts.t1,
            probe.state.lts.t2);
    REQUIRE(parton_born);
    REQUIRE(probe.state.lts.x1 == Approx(0.20).epsilon(1.0e-14));
    REQUIRE(probe.state.lts.x2 == Approx(0.15).epsilon(1.0e-14));
    REQUIRE(probe.state.lts.pfinal[0].M2() ==
            Approx(probe.state.lts.s * 0.20 * 0.15).epsilon(1.0e-11));
    REQUIRE_FALSE(probe.state.lts.has_xi1);
    REQUIRE_FALSE(probe.state.lts.has_xi2);
    REQUIRE(probe.state.lts.t1 <= 0.0);
    REQUIRE(probe.state.lts.t2 <= 0.0);
    REQUIRE(std::abs(probe.state.lts.t1) < 1.0e-9);
    REQUIRE(std::abs(probe.state.lts.t2) < 1.0e-9);
    REQUIRE(probe.state.lts.q1.M2() == Approx(0.0).margin(1.0e-12));
    REQUIRE(probe.state.lts.q2.M2() == Approx(0.0).margin(1.0e-12));
    RequireFourVectorNear(probe.state.lts.q1,
                          probe.state.lts.pbeam1 - probe.state.lts.pfinal[1],
                          1.0e-11);
    RequireFourVectorNear(probe.state.lts.q2,
                          probe.state.lts.pbeam2 - probe.state.lts.pfinal[2],
                          1.0e-11);
    RequireFourVectorNear(probe.state.lts.q1 + probe.state.lts.q2,
                          probe.state.lts.pfinal[0], 1.0e-11);
    REQUIRE(
        CollisionCMVector(probe.state.lts, probe.state.lts.pfinal[0]).Rap() ==
        Approx(0.5 * std::log(0.20 / 0.15)).margin(1.0e-11));
    const std::vector<gra::M4Vec> born = probe.state.lts.pfinal;
    probe.SetScreening(false);
    REQUIRE_FALSE(probe.GetScreening());
    RequirePhaseSpaceClosure(probe.state.lts, true);

    probe.SetScreening(true);
    REQUIRE(probe.GetScreening());
    REQUIRE_THROWS_AS(probe.FinalizeProcessConfiguration(), std::invalid_argument);
    probe.SetScreening(false);
    REQUIRE_NOTHROW(probe.FinalizeProcessConfiguration());
    REQUIRE_FALSE(probe.GetScreening());
    probe.state.lts.pfinal_orig = born;
    REQUIRE_FALSE(
        probe.BuildLoopForTest({born[1].Px() - 0.04, born[1].Py() + 0.03},
                               {born[2].Px() + 0.04, born[2].Py() - 0.03}));
    for (std::size_t i = 0; i < 3; ++i) {
      RequireFourVectorNear(probe.state.lts.pfinal[i], born[i]);
    }
  }
}

TEST_CASE("one production root inherits the central-system phase space",
          "[Process][PhaseSpace][Topology]") {
  SECTION("factorized") {
    FactorizedKinematicsProbe probe;
    ConfigureAsymmetricProtonBeams(probe, 80.0, 20.0);
    probe.state.lts.PS_active = true;
    probe.state.lts.decaytree = {SingleCascadeRoot(3.0)};

    REQUIRE(probe.BuildBornForTest(0.25, 0.35, 0.2, 2.4, 0.1, 9.0, gra::PDG::mp,
                                   gra::PDG::mp));
    RequireSingleRootKinematics(probe.state.lts);
  }

  SECTION("collinear parton") {
    CollinearKinematicsProbe probe;
    ConfigureAsymmetricProtonBeams(probe, 80.0, 20.0);
    probe.state.lts.PS_active = true;
    probe.state.lts.decaytree = {
        SingleCascadeRoot(std::sqrt(0.20 * 0.15 * probe.state.lts.s))};

    REQUIRE(probe.BuildBornForTest(0.20, 0.15));
    RequireSingleRootKinematics(probe.state.lts);
  }

  SECTION("hard diffraction") {
    HardDiffractionKinematicsProbe probe;
    ConfigureAsymmetricProtonBeams(probe, 80.0, 20.0);
    probe.state.lts.PS_active = true;
    probe.state.lts.decaytree = {SingleCascadeRoot(10.0, 10.0)};
    probe.state.lts.hard_diff1 = true;
    probe.state.lts.hard_diff2 = true;
    probe.state.lts.diff_xi1 = 0.12;
    probe.state.lts.diff_xi2 = 0.10;
    probe.state.lts.diff_beta1 = 0.65;
    probe.state.lts.diff_beta2 = 0.55;
    probe.state.lts.diff_t1 = -0.16;
    probe.state.lts.diff_t2 = -0.09;
    probe.state.lts.diff_phi1 = 0.3;
    probe.state.lts.diff_phi2 = 0.3 + gra::math::PI;
    probe.state.lts.id1 = gra::PDG::PDG_gluon;
    probe.state.lts.id2 = gra::PDG::PDG_gluon;

    REQUIRE(probe.BuildBornForTest(
        probe.state.lts.diff_xi1 * probe.state.lts.diff_beta1,
        probe.state.lts.diff_xi2 * probe.state.lts.diff_beta2));
    RequireSingleRootKinematics(probe.state.lts);
  }

  SECTION("direct phase-space volume") {
    ProcessKinematicsProbe probe;
    probe.state.lts.decaytree = {SingleCascadeRoot(3.0)};
    probe.state.lts.pfinal[0] = gra::M4Vec(0.0, 0.0, 0.0, 3.0);
    probe.state.lts.decaytree[0].p4 = probe.state.lts.pfinal[0];
    probe.state.lts.m2 = probe.state.lts.pfinal[0].M2();
    double mass_sum = 0.0;
    double mass_max = 100.0;
    REQUIRE(probe.ProbePrepareCentralBranchMasses(mass_sum, &mass_max));
    REQUIRE(mass_sum == Approx(2.0));
    REQUIRE(mass_max == Approx(4.0));
    REQUIRE(probe.ProbeCentralDecayWidthPS() == Approx(1.0));
  }

  SECTION("stable one-root topology is rejected") {
    ProcessKinematicsProbe probe;
    probe.state.lts.decaytree = {MasslessLeaf()};
    probe.state.lts.pfinal[0] = gra::M4Vec(0.0, 0.0, 0.0, 3.0);
    double mass_sum = 0.0;
    REQUIRE_FALSE(probe.ProbePrepareCentralBranchMasses(mass_sum));
    REQUIRE_FALSE(probe.ProbeSetSingleCentralRootKinematics());
    REQUIRE(probe.ProbeCentralDecayWidthPS() == Approx(0.0));
  }

  SECTION("out-of-window root virtuality is rejected") {
    ProcessKinematicsProbe probe;
    probe.state.lts.decaytree = {SingleCascadeRoot(2.0, 0.2)};
    probe.state.lts.pfinal[0] = gra::M4Vec(0.0, 0.0, 0.0, 4.0);
    REQUIRE_FALSE(probe.ProbeSetSingleCentralRootKinematics());
    probe.state.lts.decaytree[0].p4 = probe.state.lts.pfinal[0];
    REQUIRE(probe.ProbeCentralDecayWidthPS() == Approx(0.0));
  }
}

TEST_CASE("central transverse map follows the event mass ceiling",
          "[Process][CentralPhaseSpace]") {
  const gra::M4Vec forward1(3.0, 4.0, 0.0, 0.0);
  const gra::M4Vec forward2(-1.0, 2.0, 0.0, 0.0);
  constexpr double sqrt_s = 20.0;
  constexpr double mass1 = 3.0;
  constexpr double mass2 = 4.0;

  const double available_energy = sqrt_s -
                                  std::sqrt(mass1 * mass1 + forward1.Pt2()) -
                                  std::sqrt(mass2 * mass2 + forward2.Pt2());
  const double expected_mass = std::sqrt(available_energy * available_energy -
                                         (forward1 + forward2).Pt2());
  const double mass_ceiling = CentralKinematicsProbe::MassCeiling(
      sqrt_s, mass1, mass2, forward1, forward2);
  REQUIRE(mass_ceiling == Approx(expected_mass).epsilon(1.0e-14));

  constexpr double radius2 = 36.0;
  REQUIRE(CentralKinematicsProbe::TransverseRadius2(6.0, 24.0, 4) ==
          Approx(radius2));
  REQUIRE(std::sqrt(2.0 * CentralKinematicsProbe::TransverseRadius2(
                              6.0, 24.0, 4)) == Approx(std::sqrt(72.0)));
  const gra::M4Vec center(0.5, -0.25, 0.0, 0.0);
  const std::vector<double> radial_units = {0.37, 0.50, 0.75};
  const std::vector<double> angle_units = {0.20, 0.40, 0.60};
  constexpr double scale2 = 4.0;
  const gra::M4Vec central_momentum(2.0, -1.0, 0.0, 0.0);
  const auto [log_differences, log_jacobian] =
      CentralKinematicsProbe::MapTransverseLogarithmic(
          radial_units, angle_units, central_momentum, radius2, scale2);
  REQUIRE(log_differences.size() == 3);
  std::vector<gra::M4Vec> log_momenta;
  REQUIRE(gra::kinematics::BuildCentralTransverseMomenta(
      -central_momentum, gra::M4Vec(), log_differences, log_momenta));
  const double log_range = std::log1p(radius2 / scale2);
  const double rho2 = scale2 * std::expm1(radial_units[0] * log_range);
  double log_centered_norm2 = 0.0;
  gra::M4Vec momentum_sum;
  for (const auto &momentum : log_momenta) {
    momentum_sum += momentum;
    log_centered_norm2 += (momentum - center).Pt2();
  }
  REQUIRE(momentum_sum.Px() == Approx(central_momentum.Px()));
  REQUIRE(momentum_sum.Py() == Approx(central_momentum.Py()));
  REQUIRE(log_centered_norm2 == Approx(rho2).epsilon(1.0e-13));
  const double expected_log_jacobian = 4.0 * std::pow(gra::math::PI, 3) * rho2 *
                                       rho2 * log_range * (scale2 + rho2) / 2.0;
  REQUIRE(log_jacobian == Approx(expected_log_jacobian).epsilon(1.0e-13));

  // Rotate the full transverse construction without changing its Jacobian
  constexpr double angle_shift = 0.125;
  const double rotation = 2.0 * gra::math::PI * angle_shift;
  const double cosine = std::cos(rotation);
  const double sine = std::sin(rotation);
  auto rotated_angles = angle_units;
  for (auto &angle : rotated_angles) {
    angle += angle_shift;
  }
  const gra::M4Vec rotated_central(
      cosine * central_momentum.Px() - sine * central_momentum.Py(),
      sine * central_momentum.Px() + cosine * central_momentum.Py(), 0.0, 0.0);
  const auto [rotated_differences, rotated_jacobian] =
      CentralKinematicsProbe::MapTransverseLogarithmic(
          radial_units, rotated_angles, rotated_central, radius2, scale2);
  REQUIRE(rotated_jacobian == Approx(log_jacobian).epsilon(1.0e-14));
  for (const auto &i : indices(log_differences)) {
    const double rotated_x =
        cosine * log_differences[i].Px() - sine * log_differences[i].Py();
    const double rotated_y =
        sine * log_differences[i].Px() + cosine * log_differences[i].Py();
    REQUIRE(rotated_differences[i].Px() == Approx(rotated_x).margin(1.0e-13));
    REQUIRE(rotated_differences[i].Py() == Approx(rotated_y).margin(1.0e-13));
  }

  // Integrate the logarithmic radial Jacobian back to the full Helmert ball
  constexpr std::size_t radial_nodes = 4096;
  double mapped_volume = 0.0;
  for (std::size_t i = 0; i < radial_nodes; ++i) {
    auto units = radial_units;
    units[0] =
        (static_cast<double>(i) + 0.5) / static_cast<double>(radial_nodes);
    mapped_volume += CentralKinematicsProbe::MapTransverseLogarithmic(
                         units, angle_units, central_momentum, radius2, scale2)
                         .second;
  }
  mapped_volume /= static_cast<double>(radial_nodes);
  const double ball_volume = 4.0 * std::pow(gra::math::PI * radius2, 3) / 6.0;
  REQUIRE(mapped_volume == Approx(ball_volume).epsilon(1.0e-6));

  auto invalid_radial_units = radial_units;
  invalid_radial_units[0] = -0.1;
  REQUIRE_THROWS_AS(
      CentralKinematicsProbe::MapTransverseLogarithmic(
          invalid_radial_units, angle_units, central_momentum, radius2, scale2),
      gra::PhaseSpaceFailure);
  REQUIRE_THROWS_AS(
      CentralKinematicsProbe::MapTransverseLogarithmic(
          radial_units, angle_units, central_momentum, radius2, 0.0),
      gra::PhaseSpaceFailure);
  REQUIRE_THROWS_AS(CentralKinematicsProbe::TransverseRadius2(6.0, 24.0, 1),
                    gra::PhaseSpaceFailure);
  REQUIRE(gra::math::IsZero(CentralKinematicsProbe::MassCeiling(
      2.0, mass1, mass2, forward1, forward2)));
}

TEST_CASE("factorized central mass map has the exact logarithmic Jacobian",
          "[Process][FactorizedPhaseSpace]") {
  constexpr double mass_min = 20.0;
  constexpr double mass_max = 7000.0;
  constexpr double unit = 0.37;
  constexpr double step = 1.0e-6;

  const double lower =
      FactorizedKinematicsProbe::SampleCentralMassSquaredForTest(0.0, mass_min,
                                                                 mass_max);
  const double upper =
      FactorizedKinematicsProbe::SampleCentralMassSquaredForTest(1.0, mass_min,
                                                                 mass_max);
  const double midpoint =
      FactorizedKinematicsProbe::SampleCentralMassSquaredForTest(0.5, mass_min,
                                                                 mass_max);
  REQUIRE(lower == Approx(mass_min * mass_min).epsilon(1.0e-14));
  REQUIRE(upper == Approx(mass_max * mass_max).epsilon(1.0e-14));
  REQUIRE(midpoint == Approx(mass_min * mass_max).epsilon(1.0e-14));

  const double mass_squared =
      FactorizedKinematicsProbe::SampleCentralMassSquaredForTest(unit, mass_min,
                                                                 mass_max);
  const double jacobian =
      FactorizedKinematicsProbe::CentralMassSquaredJacobianForTest(
          mass_squared, mass_min, mass_max);
  const double plus =
      FactorizedKinematicsProbe::SampleCentralMassSquaredForTest(
          unit + step, mass_min, mass_max);
  const double minus =
      FactorizedKinematicsProbe::SampleCentralMassSquaredForTest(
          unit - step, mass_min, mass_max);
  REQUIRE((plus - minus) / (2.0 * step) == Approx(jacobian).epsilon(1.0e-9));
}

// Check light lepton mass shells through asymmetric collinear photon kinematics
TEST_CASE("collinear photon pairs preserve mass shells with asymmetric beams",
          "[Process][CollinearPhaseSpace][mass-shell]") {
  auto card = nlohmann::json::parse(gra::aux::GetInputData(gra::aux::ResolveProjectPath("gencard/test.json")));
  card["GENERIC"]["CORES"] = 1;
  card["GENERIC"]["HIST"] = 0;
  card["GENERIC"]["INTEGRATOR"] = "VEGAS";
  card["SCATTERING"]["LOOPSCREEN"] = false;
  card["SCATTERING"]["RES"] = nlohmann::json::array();
  card["GENCUTS"]["<P>"] = {{"M", {1.0, 2.5}}, {"Rap", {-4.0, 4.0}}};
  card["FIDCUTS"] = {{"active", false}};
  for (const std::string leptons : {"e+ e-", "mu+ mu-"}) {
    for (const double ratio : {0.01, 1.0, 100.0}) {
      card["SCATTERING"]["PROCESS"] = "yy_DZ[EPA]<P> -> " + leptons;
      card["SCATTERING"]["ENERGY"] = {6500.0, 6500.0 * ratio};
      gra::MGraniitti generator;
      generator.ReadInput(card);
      generator.proc->PrepareRun();
      for (const double mass : {0.15, 0.5, 0.85}) {
        for (const double rapidity : {0.05, 0.275, 0.5, 0.725, 0.95}) {
          for (unsigned int trial = 0; trial < 32; ++trial) {
            generator.proc->state.random.SetSeed(7103 + trial);
            gra::MEventWeightState state;
            state.include_screening = false;
            const double weight = generator.proc->EventWeight({mass, rapidity}, state);
            CAPTURE(leptons, ratio, mass, rapidity, trial, state.technical_failure);
            REQUIRE(state.Valid());
            REQUIRE(weight > 0.0);
            const auto &lts = generator.proc->state.lts;
            CHECK(lts.pfinal[0].M() == Approx(std::pow(2.5, mass)).epsilon(1e-10));
            CHECK(lts.pfinal[0].Rap() == Approx(-4.0 + 8.0 * rapidity).margin(1e-11));
            REQUIRE(gra::math::CheckEMC(lts.decaytree[0].p4 + lts.decaytree[1].p4 - lts.pfinal[0]));
            for (const auto &branch : lts.decaytree) {
              const double scale = std::max(1.0, branch.p4.E() * branch.p4.E());
              CHECK(std::abs(branch.p4.M2() - branch.p.mass * branch.p.mass) <
                    64.0 * std::numeric_limits<double>::epsilon() * scale);
            }
          }
        }
      }
    }
  }
}

// Check the stored inverse density belongs to the generated collinear mass map
TEST_CASE("collinear cascades record their conditional central mass proposal",
          "[Process][CollinearPhaseSpace][Cascade][proposal]") {
  auto card = nlohmann::json::parse(gra::aux::GetInputData(gra::aux::ResolveProjectPath("gencard/test.json")));
  card["GENERIC"]["CORES"] = 1;
  card["GENERIC"]["HIST"] = 0;
  card["GENERIC"]["INTEGRATOR"] = "VEGAS";
  card["SCATTERING"]["PROCESS"] = "yy_DZ[FLUX]<P> &> rho(770)0 > {pi+ pi-} rho(770)0 > {pi+ pi-}";
  card["SCATTERING"]["ENERGY"] = {80.0, 20.0};
  card["SCATTERING"]["LOOPSCREEN"] = false;
  card["GENCUTS"]["<P>"] = {{"M", {1.5, 5.0}}, {"Rap", {-0.2, 0.2}}};
  card["FIDCUTS"] = {{"active", false}};
  gra::MGraniitti generator;
  generator.ReadInput(card);
  auto &process = *generator.proc;
  process.SetOFFSHELL(0.0);
  process.PrepareRun();
  gra::MEventWeightState aux;
  aux.include_screening = false;
  const double weight = process.EventWeight({0.37, 0.62}, aux);
  REQUIRE(std::isfinite(weight));
  REQUIRE(weight > 0.0);
  REQUIRE(aux.Valid());
  const auto &lts = process.state.lts;
  REQUIRE(lts.central_phase_space_mode == gra::CentralPhaseSpaceMode::Collinear);
  const double lower = lts.decaytree[0].p.mass + lts.decaytree[1].p.mass + lts.central_phase_space_mass_margin;
  const double log_range = std::log(5.0 / lower);
  REQUIRE(std::log(lts.pfinal[0].M() / lower) / log_range == Approx(0.37).epsilon(1e-11));
  REQUIRE(lts.central_phase_space_generated_jacobian ==
          Approx(2.0 * lts.pfinal[0].M2() * log_range).epsilon(1e-11));
}

TEST_CASE("parton central mass map has the exact logarithmic Jacobian",
          "[Process][CollinearPhaseSpace]") {
  constexpr double mass_min = 1.0e-4;
  constexpr double mass_max = 7000.0;
  constexpr double unit = 0.63;
  constexpr double step = 1.0e-6;

  const double lower =
      CollinearKinematicsProbe::SampleCentralMassSquaredForTest(0.0, mass_min,
                                                                mass_max);
  const double upper =
      CollinearKinematicsProbe::SampleCentralMassSquaredForTest(1.0, mass_min,
                                                                mass_max);
  const double midpoint =
      CollinearKinematicsProbe::SampleCentralMassSquaredForTest(0.5, mass_min,
                                                                mass_max);
  REQUIRE(lower == Approx(mass_min * mass_min).epsilon(1.0e-14));
  REQUIRE(upper == Approx(mass_max * mass_max).epsilon(1.0e-14));
  REQUIRE(midpoint == Approx(mass_min * mass_max).epsilon(1.0e-14));

  const double mass_squared =
      CollinearKinematicsProbe::SampleCentralMassSquaredForTest(unit, mass_min,
                                                                mass_max);
  const double jacobian =
      CollinearKinematicsProbe::CentralMassSquaredJacobianForTest(
          mass_squared, mass_min, mass_max);
  const double plus = CollinearKinematicsProbe::SampleCentralMassSquaredForTest(
      unit + step, mass_min, mass_max);
  const double minus =
      CollinearKinematicsProbe::SampleCentralMassSquaredForTest(
          unit - step, mass_min, mass_max);
  REQUIRE((plus - minus) / (2.0 * step) == Approx(jacobian).epsilon(1.0e-9));
}

TEST_CASE("veto domains select scattering source and particle charge category",
          "[Process][VetoCuts]") {
  ProcessKinematicsProbe probe;
  gra::VETODOMAIN domain;
  domain.eta_min = -2.5;
  domain.eta_max = 2.5;
  domain.pt_min = 0.4;
  domain.pt_max = 100000.0;
  domain.source_forward = true;
  domain.source_central = false;

  gra::VETOCUT cuts;
  cuts.active = true;

  const auto charged_track = StableParticle(3, 1.0, 0.5);
  const auto neutral_track = StableParticle(0, 1.0, 0.5);

  REQUIRE(domain.charge == gra::VetoCharge::Any);
  domain.charge = gra::VetoCharge::Charged;
  cuts.cuts = {domain};
  probe.ConfigureVetoCuts(cuts);
  REQUIRE_FALSE(probe.ProbeVetoCuts(charged_track, neutral_track));
  REQUIRE(probe.ProbeVetoCuts(neutral_track, charged_track));

  domain.charge = gra::VetoCharge::Neutral;
  cuts.cuts = {domain};
  probe.ConfigureVetoCuts(cuts);
  REQUIRE_FALSE(probe.ProbeVetoCuts(neutral_track, charged_track));
  REQUIRE(probe.ProbeVetoCuts(charged_track, neutral_track));

  domain.charge = gra::VetoCharge::Any;
  cuts.cuts = {domain};
  probe.ConfigureVetoCuts(cuts);
  REQUIRE_FALSE(probe.ProbeVetoCuts(charged_track, neutral_track));
  REQUIRE_FALSE(probe.ProbeVetoCuts(neutral_track, charged_track));
}

TEST_CASE(
    "single-diffraction forward volume describes one sampled beam assignment",
    "[Process][ForwardPhaseSpace]") {
  ProcessKinematicsProbe probe;
  probe.SetModelTune(gra::MModelTune::Load(
      gra::aux::ResolveProjectPath("modeldata/TUNE0/GENERAL.json")));
  const double s = 10000.0;
  const double xi_min = 0.02;
  const double xi_max = 0.20;
  const double pt_min = 0.10;
  const double pt_max = 1.40;
  const double pt1 = 0.37;
  const double pt2 = 0.62;
  probe.ConfigureForward(1, s, xi_min, xi_max, pt_min, pt_max);
  const std::vector<double> masses = probe.SampleMasses({0.35});
  probe.SetForwardPoint(masses, pt1, pt2);

  const auto bounds = probe.ForwardMassBounds();
  REQUIRE(bounds.first == Approx(xi_min * s).margin(1e-13));
  REQUIRE(bounds.second == Approx(xi_max * s).margin(1e-13));
  const double excited_mass2 = probe.state.lts.excite1
                                   ? probe.state.lts.pfinal[1].M2()
                                   : probe.state.lts.pfinal[2].M2();
  REQUIRE(probe.state.lts.excite1 != probe.state.lts.excite2);
  const double phi_volume = gra::math::pow2(2.0 * gra::math::PI);
  const double pt_volume =
      pt1 * pt2 * gra::math::pow2(std::log(pt_max) - std::log(pt_min + 1e-12));
  const double expected = phi_volume * pt_volume * excited_mass2 *
                          std::log(bounds.second / bounds.first);
  REQUIRE(probe.ProbeForwardVolume() == Approx(expected).epsilon(2e-13));
  const double p = probe.state.nstar_param->single_side_prob;
  const double side_probability = probe.state.lts.excite1 ? p : 1.0 - p;
  REQUIRE(probe.ProbeDissociationCrossSectionFactor() *
              probe.ProbeForwardVolume() ==
          Approx(expected / side_probability).epsilon(2e-13));
}

// Check the lepton or nucleus remains unchanged and the sole proton assignment has unit probability
TEST_CASE("ep and pA NSTAR samples only the proton in either beam ordering", "[Process][Dissociation][Normalization]") {
  const auto partner = GENERATE(std::make_pair(11, 0.000511), std::make_pair(1000822080, 193.7));
  ProcessKinematicsProbe probe;
  probe.SetModelTune(gra::MModelTune::Load(gra::aux::ResolveProjectPath("modeldata/TUNE0/GENERAL.json")));
  for (const bool upper_proton : {false, true}) {
    probe.ConfigureForward(1, 1.0e8, 2.0e-8, 1.0e-6, 0.001, 2.0);
    probe.state.lts.beam1.pdg  = upper_proton ? gra::PDG::PDG_p : partner.first;
    probe.state.lts.beam2.pdg  = upper_proton ? partner.first : gra::PDG::PDG_p;
    probe.state.lts.beam1.mass = upper_proton ? gra::PDG::mp : partner.second;
    probe.state.lts.beam2.mass = upper_proton ? partner.second : gra::PDG::mp;
    for (const double u : {0.0, 0.2, 0.8, 1.0}) {
      const auto masses = probe.SampleMasses({u});
      CHECK(probe.state.lts.excite1 == upper_proton);
      CHECK(probe.state.lts.excite2 == !upper_proton);
      CHECK(masses[upper_proton ? 1 : 0] == Approx(partner.second));
      CHECK(masses[upper_proton ? 0 : 1] == Approx(std::sqrt(2.0 * std::pow(50.0, u))));
      CHECK(probe.ProbeDissociationCrossSectionFactor() == Approx(1.0));
    }
    probe.state.lts.beam1.chargeX3 = upper_proton ? 3 : gra::nuclear::BeamChargeX3(partner.first);
    probe.state.lts.beam2.chargeX3 = upper_proton ? gra::nuclear::BeamChargeX3(partner.first) : 3;
    probe.ConfigureDissociationTopology(2, "yy", "EPA");
    CHECK_THROWS_AS(probe.PrepareRun(), std::invalid_argument);
  }
}

TEST_CASE("single dissociation sums both physical beam assignments once",
          "[Process][Dissociation][Normalization]") {
  struct Topology {
    int excitation;
    std::string initial_state;
    std::string channel;
    double expected;
  };

  const std::vector<Topology> topologies = {
      {0, "MP", "CON", 1.0},  {1, "MP", "CON", 2.0}, {2, "MP", "CON", 1.0},
      {0, "X", "EL", 1.0},    {0, "X", "SD", 2.0},   {0, "X", "DD", 1.0},
      {0, "X", "ND", 1.0},    {0, "IPp", "Z", 2.0},  {0, "IPIP", "Z", 1.0},
      {0, "yy_DZ", "Z", 1.0},
  };

  ProcessKinematicsProbe probe;
  auto nstar = std::make_shared<gra::MNstarParam>();
  nstar->single_side_prob = 0.5;
  probe.state.nstar_param = nstar;
  probe.state.lts.excite1 = true;
  probe.state.lts.excite2 = false;
  for (const auto &topology : topologies) {
    probe.ConfigureDissociationTopology(
        topology.excitation, topology.initial_state, topology.channel);
    INFO("topology = " << topology.initial_state << "[" << topology.channel
                       << "], excitation = " << topology.excitation);
    REQUIRE(probe.ProbeDissociationCrossSectionFactor() ==
            Approx(topology.expected));
  }
}

// Preserve the nominal on-shell beams when a direct energy update is invalid
TEST_CASE("beam energy updates reject unphysical input before changing the collision",
          "[Process][Beams][Validation]") {
  ProcessKinematicsProbe probe;
  ConfigureAsymmetricProtonBeams(probe, 7.0, 3.0);
  const auto beam1 = probe.state.lts.pbeam1;
  const auto beam2 = probe.state.lts.pbeam2;
  const double s = probe.state.lts.s;
  for (const double energy : {-1.0, 0.0, 0.5 * gra::PDG::mp,
                              std::numeric_limits<double>::infinity(),
                              std::numeric_limits<double>::quiet_NaN(), 1e300}) {
    CAPTURE(energy);
    CHECK_THROWS_AS(probe.SetBeamEnergies(energy, 3.0), std::invalid_argument);
    CHECK_THROWS_AS(probe.SetBeamEnergies(7.0, energy), std::invalid_argument);
    RequireFourVectorNear(probe.state.lts.pbeam1, beam1);
    RequireFourVectorNear(probe.state.lts.pbeam2, beam2);
    CHECK(probe.state.lts.s == Approx(s));
  }
  for (const double energy : {-1.0, std::numeric_limits<double>::infinity(),
                              std::numeric_limits<double>::quiet_NaN()}) {
    CHECK_THROWS_AS(probe.SetInitialState({"p+", "p+"}, {energy, 3.0}), std::invalid_argument);
    CHECK_THROWS_AS(probe.SetInitialState({"p+", "p+"}, {7.0, energy}), std::invalid_argument);
  }
  probe.SetBeamEnergies(3.0, 7.0);
  CHECK(probe.state.lts.s == Approx(s));
  const double mass = probe.state.lts.beam2.mass;
  for (const double energy : {mass, std::nextafter(mass, 2.0)}) {
    probe.SetBeamEnergies(7.0, energy);
    CHECK(probe.state.lts.pbeam2.M2() == Approx(mass * mass));
    CHECK(probe.state.lts.pbeam2.E() == Approx(energy));
    CHECK(probe.state.lts.pbeam2.Pz() <= 0.0);
  }
}

// Check stationary targets against the invariant collision energy in either beam order
TEST_CASE("beam initialization preserves stationary targets and threshold energies",
          "[Process][Beams][physics][regression]") {
  gra::MFactorized process;
  process.state.lts.PDG.ReadParticleData();
  const double mass = process.state.lts.PDG.FindByPDG(gra::PDG::PDG_p).mass;
  const double energy = 7.0;
  // [REFERENCE: PDG 2025, Kinematics, Eqs. (49.2)-(49.3)]
  const double s = 2.0 * mass * (mass + energy);
  for (const bool upper_target : {false, true}) {
    const std::vector<double> beams = upper_target ? std::vector<double>{mass, energy}
                                                   : std::vector<double>{energy, mass};
    process.SetInitialState({"p+", "p+"}, beams);
    const auto &target = upper_target ? process.state.lts.pbeam1 : process.state.lts.pbeam2;
    CHECK(target.E() == Approx(mass).epsilon(1.0e-14));
    CHECK(target.Pz() == Approx(0.0).margin(1.0e-14));
    CHECK(process.state.lts.s == Approx(s).epsilon(1.0e-14));
  }
  const double near_mass = std::nextafter(mass, energy);
  process.SetInitialState({"p+", "p+"}, {energy, near_mass});
  CHECK(process.state.lts.pbeam2.E() == Approx(near_mass).epsilon(1.0e-14));
  CHECK(process.state.lts.pbeam2.Pz() ==
        Approx(-std::sqrt((near_mass - mass) * (near_mass + mass))).epsilon(1.0e-14));
  for (const double invalid : {0.0, 0.5 * mass, std::nextafter(mass, 0.0)}) {
    CHECK_THROWS_AS(process.SetInitialState({"p+", "p+"}, {energy, invalid}), std::invalid_argument);
    CHECK_THROWS_AS(process.SetInitialState({"p+", "p+"}, {invalid, energy}), std::invalid_argument);
  }
}

// Conserve baryon charge in either valence channel and apply vetoes to the resulting diquark
TEST_CASE("forward strings conserve charge and retain diquark veto selection",
          "[Process][ForwardPhaseSpace][VetoCuts][physics][regression]") {
  ProcessKinematicsProbe process;
  process.state.lts.PDG.ReadParticleData();
  process.SetModelTune(gra::MModelTune::Load(gra::aux::ResolveProjectPath("modeldata/TUNE0/GENERAL.json")));
  auto param = std::make_shared<gra::MNstarParam>(*process.state.nstar_param);
  process.state.nstar_param = param;
  const int pdg = GENERATE(2212, -2212, 2112, -2112);
  const auto beam = process.state.lts.PDG.FindByPDG(pdg);
  const gra::M4Vec parent(0.0, 0.0, 0.0, 2.0);
  for (const double probability : {0.0, 1.0}) {
    CAPTURE(pdg, probability);
    param->string.first_valence_prob = probability;
    gra::MDecayBranch branch;
    REQUIRE(process.ExciteString(parent, branch, beam, 701));
    REQUIRE(branch.legs.size() == 2);
    CHECK(branch.legs[0].p.chargeX3 + branch.legs[1].p.chargeX3 == beam.chargeX3);
    CHECK(gra::math::CheckEMC(parent - branch.legs[0].p4 - branch.legs[1].p4));

    gra::VETODOMAIN domain;
    domain.source_forward = true;
    domain.source_central = false;
    domain.eta_min = -100.0;
    domain.eta_max = 100.0;
    domain.pt_min = 0.0;
    domain.pt_max = parent.M();
    gra::VETOCUT cuts;
    cuts.active = true;
    for (const auto charge : {gra::VetoCharge::Charged, gra::VetoCharge::Neutral}) {
      domain.charge = charge;
      cuts.cuts = {domain};
      process.ConfigureVetoCuts(cuts);
      CHECK(process.ProbeVetoCuts(branch.legs[1], {}) == (charge == gra::VetoCharge::Neutral));
    }
  }
}

// Reject fixed forward input at initialization and sampled coordinates during generation
TEST_CASE("forward mass initialization validates input and sampling validates coordinates",
          "[Process][ForwardPhaseSpace][Validation]") {
  ProcessKinematicsProbe probe;
  auto nstar = std::make_shared<gra::MNstarParam>();
  nstar->single_side_prob = 0.5;
  probe.state.nstar_param = nstar;
  probe.ConfigureForward(2, 10000.0, 0.02, 0.2, 0.1, 1.4);
  for (const auto &coordinates : std::vector<std::vector<double>>{
           {}, {0.5}, {-0.1, 0.5}, {0.5, 1.1},
           {std::numeric_limits<double>::quiet_NaN(), 0.5},
           {0.5, std::numeric_limits<double>::infinity()}}) {
    CHECK_THROWS_AS(probe.SampleMasses(coordinates), gra::PhaseSpaceFailure);
    CHECK_FALSE(probe.state.lts.excite1);
    CHECK_FALSE(probe.state.lts.excite2);
  }
  for (const double coordinate : {0.0, 1.0}) {
    const auto masses = probe.SampleMasses({coordinate, coordinate});
    const double expected2 = coordinate < 0.5 ? 200.0 : 2000.0;
    for (const double mass : masses) { CHECK(mass * mass == Approx(expected2)); }
  }
  for (const double energy2 : {-1.0, 0.0, 1.0,
                               std::numeric_limits<double>::infinity(),
                               std::numeric_limits<double>::quiet_NaN()}) {
    probe.ConfigureForward(2, energy2, 0.0, 0.2, 0.1, 1.4);
    CHECK_THROWS_AS(gra::forward::MassBounds(probe.state, probe.state.gcuts), std::invalid_argument);
  }
  for (const double xi : {-0.1, 0.0, 0.01,
                          std::numeric_limits<double>::infinity(),
                          std::numeric_limits<double>::quiet_NaN()}) {
    probe.ConfigureForward(2, 10000.0, 0.02, xi, 0.1, 1.4);
    CHECK_THROWS_AS(gra::forward::MassBounds(probe.state, probe.state.gcuts), std::invalid_argument);
  }
  probe.ConfigureForward(1, 10000.0, 0.02, 0.2, 0.1, 1.4);
  for (const double probability : {0.0, 1.0, std::numeric_limits<double>::quiet_NaN()}) {
    nstar->single_side_prob = probability;
    CHECK_THROWS_AS(gra::forward::MassBounds(probe.state, probe.state.gcuts), std::invalid_argument);
  }
}

// Distinguish invalid input from a forward interval lost at lower collision energy
TEST_CASE("factorized dissociation rejects invalid input and records sampling failures",
          "[Process][ForwardPhaseSpace][Failure]") {
  gra::MFactorized process;
  process.ProcPtr.Initialize("MP", "RES");
  ConfigureAsymmetricProtonBeams(process, 50.0, 50.0);
  process.state.lts.decaytree = MuonPair(2.0);
  auto nstar = std::make_shared<gra::MNstarParam>();
  nstar->single_side_prob = 0.5;
  process.state.nstar_param = nstar;
  process.SetExcitation(1);
  process.state.gcuts.forward_pt_min = 0.1;
  process.state.gcuts.forward_pt_max = 1.0;
  process.state.gcuts.M_min = 2.0;
  process.state.gcuts.M_max = 3.0;
  process.state.gcuts.Y_min = -1.0;
  process.state.gcuts.Y_max = 1.0;
  process.state.gcuts.XI_min = 0.0;
  for (const double xi : {-0.1, 0.0, 1e-8, std::numeric_limits<double>::infinity(),
                          std::numeric_limits<double>::quiet_NaN()}) {
    process.state.gcuts.XI_max = xi;
    REQUIRE_THROWS_AS(process.FinalizeProcessConfiguration(), std::invalid_argument);
  }
  process.state.gcuts.XI_max = 0.02;
  REQUIRE_NOTHROW(process.FinalizeProcessConfiguration());
  process.SetBeamEnergies(2.0, 2.0);
  gra::MEventWeightState status;
  status.adaptation_mode = true;
  double weight = -1.0;
  REQUIRE_NOTHROW(weight = process.EventWeight(std::vector<double>(process.GetdLIPSDim(), 0.5), status));
  CHECK(gra::math::IsZero(weight));
  CHECK_FALSE(status.kinematics_ok);
  CHECK(status.technical_failure);
  CHECK_FALSE(process.state.lts.excite1);
  CHECK_FALSE(process.state.lts.excite2);
}

// Integrate both discrete beam assignments with unequal sampling probabilities
TEST_CASE("single dissociation weights remove the beam sampling probability",
          "[Process][Dissociation][Normalization]") {
  for (const double p : {0.2, 0.5, 0.8}) {
    CAPTURE(p);
    ProcessKinematicsProbe probe;
    auto nstar = std::make_shared<gra::MNstarParam>();
    nstar->single_side_prob = p;
    probe.state.nstar_param = nstar;
    probe.ConfigureForward(1, 10000.0, 0.02, 0.2, 0.1, 1.4);
    std::array<double, 2> integral = {0.0, 0.0};
    constexpr std::size_t samples = 20000;
    for (std::size_t i = 0; i < samples; ++i) {
      probe.SampleMasses({0.35});
      REQUIRE(probe.state.lts.excite1 != probe.state.lts.excite2);
      const std::size_t side = probe.state.lts.excite1 ? 0 : 1;
      const double probability = side == 0 ? p : 1.0 - p;
      const double weight = probe.ProbeDissociationCrossSectionFactor();
      REQUIRE(probability * weight == Approx(1.0).epsilon(1e-14));
      integral[side] += weight / samples;
    }
    for (const auto &side : indices(integral)) {
      const double probability = side == 0 ? p : 1.0 - p;
      const double error = std::sqrt((1.0 - probability) / (probability * samples));
      CHECK(integral[side] == Approx(1.0).margin(6.0 * error));
    }
    // Unequal physical weights also retain the sum of the two cross sections
    const double error = std::abs(3.0 / p - 0.7 / (1.0 - p)) * std::sqrt(p * (1.0 - p) / samples);
    CHECK(3.0 * integral[0] + 0.7 * integral[1] == Approx(3.7).margin(6.0 * error));
    probe.state.lts.excite1 = probe.state.lts.excite2 = false;
    REQUIRE_THROWS_AS(probe.ProbeDissociationCrossSectionFactor(), gra::PhaseSpaceFailure);
  }
}

TEST_CASE("forward excitation Xi minimum defaults to the physical threshold",
          "[Process][ForwardPhaseSpace]") {
  ProcessKinematicsProbe probe;
  probe.SetModelTune(gra::MModelTune::Load(
      gra::aux::ResolveProjectPath("modeldata/TUNE0/GENERAL.json")));
  const double s = 10000.0;
  probe.ConfigureForward(1, s, 0.0, 0.20, 0.1, 1.0);
  probe.SampleMasses({0.5});
  REQUIRE(probe.state.nstar_param != nullptr);
  const double threshold =
      gra::PDG::mp + gra::PDG::mpi + probe.state.nstar_param->mass_margin;
  REQUIRE(probe.ForwardMassBounds().first ==
          Approx(gra::math::pow2(threshold)).margin(1e-14));

  probe.ConfigureForward(1, s, 0.03, 0.20, 0.1, 1.0);
  probe.SampleMasses({0.5});
  REQUIRE(probe.ForwardMassBounds().first == Approx(0.03 * s).margin(1e-13));
}

TEST_CASE(
    "active three-body cascade weights reproduce the James phase-space volume",
    "[Process][Cascade][Literature]") {
  // F James, CERN 68-15, https://cds.cern.ch/record/275743/files/CERN-68-15.pdf
  // The tabulated normalization gives Phi_3(M^2)=M^2/(256 pi^3)
  ProcessKinematicsProbe probe;
  probe.SetActiveDecayPhaseSpace(true);
  const double mass = 6.0;
  probe.SetDecayTree(
      FixedDecay(mass, {MasslessLeaf(), MasslessLeaf(), MasslessLeaf()}));

  gra::kinematics::MCW event_weights;
  constexpr std::size_t events = 50000;
  for (std::size_t event = 0; event < events; ++event) {
    REQUIRE(probe.GenerateDecayTree());
    const auto &root = probe.state.lts.decaytree.front();
    REQUIRE(root.W_event > 0.0);
    // The externally fixed root adds no independent invariant-mass integral
    REQUIRE(probe.ProbeCascadePS() == Approx(root.W_event).epsilon(1e-14));
    event_weights.Push(root.W_event);
  }

  const double exact = mass * mass / (256.0 * gra::math::pow3(gra::math::PI));
  const auto &root = probe.state.lts.decaytree.front();
  REQUIRE(root.W.GetN() == Approx(static_cast<double>(events)));
  REQUIRE(std::abs(event_weights.Integral() - exact) <
          std::max(6.0 * event_weights.IntegralError(), 0.005 * exact));
  REQUIRE(root.W.Integral() == Approx(event_weights.Integral()).epsilon(1e-14));
}

TEST_CASE(
    "nested fixed-mass decay trees multiply every local invariant phase space",
    "[Process][Cascade]") {
  ProcessKinematicsProbe probe;
  probe.SetActiveDecayPhaseSpace(true);
  const double parent_mass = 8.0;
  const double intermediate_mass = 3.0;
  const gra::MDecayBranch intermediate =
      FixedDecay(intermediate_mass, {MasslessLeaf(), MasslessLeaf()});
  probe.SetDecayTree(FixedDecay(parent_mass, {intermediate, MasslessLeaf()}));
  REQUIRE(probe.GenerateDecayTree());

  const auto &root = probe.state.lts.decaytree.front();
  REQUIRE(root.legs.front().W_event > 0.0);
  const double expected =
      gra::kinematics::PS2Massive(parent_mass * parent_mass,
                                  intermediate_mass * intermediate_mass, 0.0) *
      gra::kinematics::PS2Massive(intermediate_mass * intermediate_mass, 0.0,
                                  0.0) /
      (2.0 * gra::math::PI);
  REQUIRE(probe.ProbeCascadePS() == Approx(expected).epsilon(1e-13));
  REQUIRE(
      gra::math::CheckEMC(root.p4 - root.legs[0].p4 - root.legs[1].p4, 1e-11));
  REQUIRE(gra::math::CheckEMC(root.legs[0].p4 - root.legs[0].legs[0].p4 -
                                  root.legs[0].legs[1].p4,
                              1e-11));

  probe.SetActiveDecayPhaseSpace(false);
  REQUIRE(probe.GenerateDecayTree());
  REQUIRE(probe.ProbeCascadePS() == Approx(1.0));
}

TEST_CASE("internal cascade branches reject invalid local phase-space weights",
          "[Process][Cascade][Failure]") {
  ProcessKinematicsProbe probe;
  probe.SetActiveDecayPhaseSpace(true);

  gra::MDecayBranch root = FixedDecay(3.0, {MasslessLeaf(), MasslessLeaf()});
  for (const double invalid_weight :
       {0.0, -1.0, std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::infinity()}) {
    root.W_event = invalid_weight;
    probe.SetDecayTree(root);
    CAPTURE(invalid_weight);
    CHECK_THROWS_AS(probe.ProbeCascadePS(), gra::PhaseSpaceFailure);
  }
}

TEST_CASE(
    "cascade phase space cancels Z and W BW proposals for full MG5 amplitudes",
    "[Process][Cascade][MG5]") {
  const std::array<std::tuple<int, double, double>, 2> resonances = {
      std::tuple{23, 91.1876, 2.4952}, std::tuple{24, 80.379, 2.085}};
  for (const auto &[pdg, mass, width] : resonances) {
    CAPTURE(pdg, mass, width);
    ProcessKinematicsProbe probe;
    probe.SetActiveDecayPhaseSpace(true);

    gra::MDecayBranch resonance =
        FixedDecay(mass, {MasslessLeaf(), MasslessLeaf()});
    resonance.p.pdg = pdg;
    resonance.p.width = width;
    probe.SetDecayTree(FixedDecay(2.5 * mass, {resonance, MasslessLeaf()}));
    REQUIRE(probe.GenerateDecayTree());

    const auto &sampled = probe.state.lts.decaytree.front().legs.front();
    REQUIRE((sampled.mass_proposal == gra::MassProposal::BreitWigner));
    REQUIRE(sampled.mass_proposal_norm > 0.0);
    const double bw2 = std::norm(gra::resonance::FixedWidthLineShape(
        sampled.p4.M2(), sampled.p.mass, sampled.p.width));
    REQUIRE(bw2 > 0.0);

    probe.state.lts.decay_structure = {gra::DecayType::JacobWickCoherent};
    const double proposal_weight = probe.ProbeCascadePS();
    probe.state.lts.decay_structure = {gra::DecayType::Full};
    const double full_amplitude_weight = probe.ProbeCascadePS();

    REQUIRE(proposal_weight > 0.0);
    REQUIRE(full_amplitude_weight ==
            Approx(proposal_weight).epsilon(1e-13));
    const auto &root = probe.state.lts.decaytree.front();
    REQUIRE(full_amplitude_weight * bw2 ==
            Approx(root.W_event * sampled.W_event * sampled.mass_proposal_norm /
                   (2.0 * gra::math::PI)).epsilon(1e-13));
  }
}

TEST_CASE("unresolved cascade poles reproduce the physical branching ratio",
          "[Process][Cascade][NWA][physics]") {
  ProcessKinematicsProbe probe;
  probe.SetActiveDecayPhaseSpace(true);
  probe.SetOffShellRange(0.0);

  gra::MDecayBranch resonance =
      FixedDecay(3.0, {MasslessLeaf(), MasslessLeaf()});
  resonance.p.width = GENERATE(0.2, 1e-9, 1e-20);
  probe.SetDecayTree(FixedDecay(8.0, {resonance, MasslessLeaf()}));
  REQUIRE(probe.GenerateDecayTree());

  const auto &sampled = probe.state.lts.decaytree.front().legs.front();
  REQUIRE((sampled.mass_proposal == gra::MassProposal::Fixed));
  const double bw2 = std::norm(gra::resonance::FixedWidthLineShape(
      sampled.p4.M2(), sampled.p.mass, sampled.p.width));
  REQUIRE(std::isfinite(bw2));
  REQUIRE(bw2 > 0.0);

  probe.state.lts.decay_structure = {gra::DecayType::JacobWickCoherent};
  const double factorized_weight = probe.ProbeCascadePS();
  probe.state.lts.decay_structure = {gra::DecayType::Full};
  const double full_weight = probe.ProbeCascadePS();
  REQUIRE(full_weight == Approx(factorized_weight).epsilon(1e-13));
  // Gamma_partial = |g|^2 Phi2 / (2M), so the integrated pole must give BR
  // [REFERENCE: PDG, Review of Particle Physics, Kinematics, Eqs. (49.10)-(49.13)]
  constexpr double branching_ratio = 0.37;
  const double coupling2 = 2.0 * sampled.p.mass * sampled.p.width * branching_ratio / sampled.W_event;
  REQUIRE(full_weight * bw2 * coupling2 ==
          Approx(probe.state.lts.decaytree.front().W_event * branching_ratio).epsilon(1e-12));
}

// A single decay history still has a nontrivial invariant-mass proposal density
TEST_CASE("one cascade history retains its mass density", "[Process][Cascade][proposal][physics]") {
  ProcessKinematicsProbe probe;
  probe.SetActiveDecayPhaseSpace(true);
  auto resonance = FixedDecay(3.0, {MasslessLeaf(), MasslessLeaf()});
  resonance.p.width = 0.2;
  resonance.legs[0].p.pdg = 11;
  resonance.legs[1].p.pdg = -11;
  auto spectator = MasslessLeaf();
  spectator.p.pdg = 22;
  probe.SetDecayTree(FixedDecay(8.0, {resonance, spectator}));
  REQUIRE(probe.GenerateDecayTree());
  auto &lts = probe.state.lts;
  lts.pfinal[0] = lts.decaytree.front().p4;
  lts.central_phase_space_mode = gra::CentralPhaseSpaceMode::Unknown;
  const auto &branch = lts.decaytree.front().legs.front();
  const double density = 1.0 / (branch.mass_proposal_norm *
      (gra::math::pow2(branch.p4.M2() - gra::math::pow2(branch.p.mass)) +
       gra::math::pow2(branch.p.mass * branch.p.width)));
  REQUIRE(gra::decay::MixtureDensity(lts, lts.decaytree) == Approx(density).epsilon(1e-11));
  const double weight = probe.ProbeCascadePS();
  lts.decay_symmetry_proposal_active = true;
  REQUIRE(probe.ProbeCascadePS() == Approx(weight).epsilon(1e-11));
}

// Check generated nested masses against an independent invariant-mass integral
TEST_CASE("nested cascade construction reproduces invariant-mass integrals and moments",
          "[Process][Cascade][proposal][physics][closure]") {
  // X(2) -> R(1.4) a, R -> S(0.7) b, S -> c d with four massless stable scalars
  // dPhi4 = dsR dsS Phi2(X,R,a) Phi2(R,S,b) Phi2(S,c,d)/(2pi)^2
  constexpr double parent_mass = 2.0;
  const auto bw2 = [](double s, double mass, double width) {
    return 1.0 / (gra::math::pow2(s - mass * mass) + gra::math::pow2(mass * width));
  };
  std::array<double, 3> exact = {};
  for (const auto &range : std::array<std::array<double, 2>, 2>{{{0.64, 1.21}, {1.21, 4.0}}}) {
    const auto [outer, outer_weight] = gra::math::GaussLegendreRule(96, range[0], range[1]);
    for (const auto &i : indices(outer)) {
      const double sr = outer[i];
      const auto [inner, inner_weight] = gra::math::GaussLegendreRule(96, 0.09, std::min(1.21, sr));
      for (const auto &j : indices(inner)) {
        const double ss = inner[j];
        const double phase = (1.0 - sr / 4.0) * (1.0 - ss / sr) /
                             (std::pow(8.0 * gra::math::PI, 3) * gra::math::pow2(2.0 * gra::math::PI));
        const double weight = outer_weight[i] * inner_weight[j] * phase * bw2(sr, 1.4, 0.3) * bw2(ss, 0.7, 0.2);
        exact[0] += weight;
        exact[1] += weight * sr / 4.0;
        exact[2] += weight * ss / 4.0;
      }
    }
  }
  for (const bool flat : {false, true}) {
    CAPTURE(flat);
    ProcessKinematicsProbe probe;
    probe.SetActiveDecayPhaseSpace(true);
    probe.state.lts.decay_structure = {gra::DecayType::Full};
    probe.SetOffShellRange(2.0);
    probe.state.lts.s = parent_mass * parent_mass;
    probe.state.flat_mass2 = flat;
    auto inner = FixedDecay(0.7, {MasslessLeaf(), MasslessLeaf()});
    inner.p.width = 0.2;
    auto outer = FixedDecay(1.4, {inner, MasslessLeaf()});
    outer.p.width = 0.3;
    probe.SetDecayTree(FixedDecay(parent_mass, {outer, MasslessLeaf()}));
    std::array<gra::kinematics::MCW, 3> moments;
    unsigned int rejected = 0;
    constexpr unsigned int events = 40000;
    for (unsigned int event = 0; event < events; ++event) {
      if (!probe.GenerateDecayTree()) {
        ++rejected;
        for (auto &moment : moments) { moment.Push(0.0); }
        continue;
      }
      const auto &r = probe.state.lts.decaytree.front().legs.front();
      const auto &s = r.legs.front();
      const double sr = r.p4.M2();
      const double ss = s.p4.M2();
      const double weight = probe.ProbeCascadePS() * bw2(sr, 1.4, 0.3) * bw2(ss, 0.7, 0.2);
      REQUIRE(std::isfinite(weight));
      moments[0].Push(weight);
      moments[1].Push(weight * sr / 4.0);
      moments[2].Push(weight * ss / 4.0);
      REQUIRE(gra::math::CheckEMC(r.p4 - s.p4 - r.legs.back().p4, 1e-10));
      REQUIRE(gra::math::CheckEMC(s.p4 - s.legs.front().p4 - s.legs.back().p4, 1e-10));
    }
    REQUIRE(rejected > 0);
    for (const auto &i : indices(moments)) {
      CAPTURE(i, rejected, exact[i], moments[i].Integral(), moments[i].IntegralError());
      REQUIRE(moments[i].GetN() == Approx(events));
      REQUIRE(std::abs(moments[i].Integral() - exact[i]) < std::max(6.0 * moments[i].IntegralError(), 0.005 * exact[i]));
    }
  }
}

TEST_CASE("all integrators share storage for eight direct central particles",
          "[Process][Topology]") {
  gra::LORENTZSCALAR scalars;
  REQUIRE(scalars.pfinal.size() == 11);
  REQUIRE(sizeof(scalars.ss) / sizeof(scalars.ss[0]) == 11);
  REQUIRE(sizeof(scalars.tt_1) / sizeof(scalars.tt_1[0]) == 11);

  ProcessKinematicsProbe probe;
  REQUIRE_THROWS_AS(probe.ProbeLorentzScalars(11), std::invalid_argument);

  gra::MCentral central;
  central.state.lts.decaytree.resize(9);
  REQUIRE_THROWS(central.FinalizeProcessConfiguration());
  gra::MFactorized factorized;
  factorized.state.lts.decaytree.resize(9);
  REQUIRE_THROWS(factorized.FinalizeProcessConfiguration());
  gra::MCollinear collinear;
  collinear.state.lts.decaytree.resize(9);
  REQUIRE_THROWS(collinear.FinalizeProcessConfiguration());
  gra::MHardDiffraction hard;
  hard.state.lts.decaytree.resize(9);
  REQUIRE_THROWS(hard.FinalizeProcessConfiguration());
  gra::MQuasiElastic quasi_elastic;
  REQUIRE(quasi_elastic.state.lts.pfinal.size() == 11);
}

TEST_CASE("nonpartonic transfers reject every positive virtuality",
          "[Process][PhaseSpace][virtuality]") {
  ProcessKinematicsProbe probe;
  probe.state.lts.pbeam1 = gra::M4Vec(0.0, 0.0, 10.0, 10.0);
  probe.state.lts.pbeam2 = gra::M4Vec(0.0, 0.0, -10.0, 10.0);
  probe.state.lts.pfinal[0] = gra::M4Vec(0.0, 0.0, 0.0, 5.0);
  probe.state.lts.pfinal[1] = gra::M4Vec(0.0, 0.0, 10.0, 10.0 - 1.0e-6);
  probe.state.lts.pfinal[2] = gra::M4Vec(0.0, 0.0, -9.0, 9.0);

  REQUIRE((probe.state.lts.pbeam1 - probe.state.lts.pfinal[1]).M2() > 0.0);
  REQUIRE_FALSE(probe.ProbeLorentzScalars(2));
}

// Reject non-finite four-momenta in the common invariant construction
TEST_CASE("common Lorentz scalars reject non-finite generated momenta", "[Process][PhaseSpace][failure]") {
  ProcessKinematicsProbe probe;
  ConfigureAsymmetricProtonBeams(probe, 10.0, 10.0);
  const double energy = std::hypot(9.0, gra::PDG::mp);
  probe.state.lts.pfinal[0] = gra::M4Vec(0.0, 0.0, 0.0, 20.0 - 2.0 * energy);
  probe.state.lts.pfinal[1] = gra::M4Vec(0.0, 0.0, 9.0, energy);
  probe.state.lts.pfinal[2] = gra::M4Vec(0.0, 0.0, -9.0, energy);
  REQUIRE(probe.ProbeLorentzScalars(2));
  const auto finite = probe.state.lts.pfinal;
  for (const double invalid : {std::numeric_limits<double>::quiet_NaN(),
                               std::numeric_limits<double>::infinity()}) {
    for (std::size_t leg = 0; leg < 3; ++leg) {
      probe.state.lts.pfinal = finite;
      probe.state.lts.pfinal[leg].SetE(invalid);
      REQUIRE_FALSE(probe.ProbeLorentzScalars(2));
    }
  }
}

// Check numerical weight bookkeeping independently of phase-space sampling
TEST_CASE("weight bookkeeping separates zero support from numerical failures",
          "[Process][Weight]") {
  ProcessKinematicsProbe probe;

  gra::MEventWeightState positive_aux;
  gra::MEventWeightState zero_aux;
  gra::MEventWeightState negative_aux;
  gra::MEventWeightState infinity_aux;
  gra::MEventWeightState nan_aux;
  gra::MEventWeightState reported_failure_aux;

  double positive = 1.0;
  probe.BookkeepAmplitudeWeight(positive, positive_aux);
  double zero = 0.0;
  probe.BookkeepAmplitudeWeight(zero, zero_aux);
  double negative = -1.0;
  probe.BookkeepAmplitudeWeight(negative, negative_aux);
  double infinity = std::numeric_limits<double>::infinity();
  probe.BookkeepAmplitudeWeight(infinity, infinity_aux);
  double nan = std::numeric_limits<double>::quiet_NaN();
  probe.BookkeepAmplitudeWeight(nan, nan_aux);
  double reported_failure = 2.0;
  reported_failure_aux.technical_failure = true;
  probe.BookkeepAmplitudeWeight(reported_failure, reported_failure_aux);
  REQUIRE(positive_aux.amplitude_ok);
  REQUIRE_FALSE(zero_aux.amplitude_ok);
  REQUIRE(negative_aux.amplitude_ok);
  REQUIRE(infinity_aux.amplitude_ok);
  REQUIRE(nan_aux.amplitude_ok);
  REQUIRE_FALSE(positive_aux.technical_failure);
  REQUIRE_FALSE(zero_aux.technical_failure);
  REQUIRE(negative_aux.technical_failure);
  REQUIRE(infinity_aux.technical_failure);
  REQUIRE(nan_aux.technical_failure);
  REQUIRE(reported_failure_aux.technical_failure);
  REQUIRE(positive_aux.Valid());
  REQUIRE_FALSE(zero_aux.Valid());
  REQUIRE_FALSE(negative_aux.Valid());
  REQUIRE_FALSE(infinity_aux.Valid());
  REQUIRE_FALSE(nan_aux.Valid());
  REQUIRE_FALSE(reported_failure_aux.Valid());
  REQUIRE(gra::math::IsExactEqual(positive, 1.0));
  REQUIRE(gra::math::IsZero(zero));
  REQUIRE(gra::math::IsZero(negative));
  REQUIRE(gra::math::IsZero(infinity));
  REQUIRE(gra::math::IsZero(nan));
  REQUIRE(gra::math::IsZero(reported_failure));
}

// Final-state radiation moves real dimuons across the measured pair-mass cut
TEST_CASE("Physical EPA dimuons apply fiducial cuts after YFS radiation", "[Process][radiative][FSR]") {
  auto generator = DimuonGenerator(true);
  auto inclusive = DimuonGenerator(true);
  auto &process = *generator->proc;
  auto &reference = *inclusive->proc;
  process.state.fcuts.active = true;
  process.state.fcuts.pdg_cuts.clear();
  auto &pair = process.state.fcuts.pdg_cuts.emplace_back();
  pair.pdg = {13, -13};
  pair.pdg_abs = {false, false};
  pair.M = {true, 19.0, 20.2};
  gra::MRandom random;
  random.SetSeed(67531);
  std::vector<double> point(process.GetdLIPSDim());
  std::size_t accepted = 0;
  std::size_t migrated = 0;
  for (std::size_t event = 0; event < 1000; ++event) {
    for (auto &unit : point) { unit = random.U(0.0, 1.0); }
    // Reproduce the same physical radiation before applying the pair cut
    reference.state.random = process.state.random;
    gra::MEventWeightState uncut;
    const double inclusive_weight = reference.EventWeight(point, uncut);
    gra::MEventWeightState aux;
    const double weight = process.EventWeight(point, aux);
    REQUIRE_FALSE(uncut.technical_failure);
    REQUIRE_FALSE(aux.technical_failure);
    if (!(inclusive_weight > 0.0)) { continue; }
    REQUIRE(uncut.Valid());
    const auto &lts = reference.state.lts;
    const auto &radiated = reference.state.radiative.FiducialTree(lts.decaytree);
    const double mass = (radiated[0].p4 + radiated[1].p4).M();
    CAPTURE(event, mass, weight, inclusive_weight);
    REQUIRE(mass <= lts.pfinal[0].M() + 1e-10);
    REQUIRE(aux.fidcuts_ok == (mass >= pair.M.min && mass <= pair.M.max));
    if (aux.fidcuts_ok) {
      REQUIRE(aux.Valid());
      REQUIRE(weight == Approx(inclusive_weight).epsilon(1e-12));
      REQUIRE(process.state.radiative.GetFSR().prepared);
      ++accepted;
    } else {
      REQUIRE(gra::math::IsZero(weight));
      ++migrated;
    }
  }
  REQUIRE(accepted > 0);
  REQUIRE(migrated > 0);
}

// Malformed coordinates must fail before sampling and leave the real process reusable
TEST_CASE("event sampling rejects malformed unit hypercube coordinates", "[Process][Weight][PhaseSpace]") {
  auto generator = DimuonGenerator();
  auto &process = *generator->proc;
  process.SetDebug(GENERATE(0, 1) != 0);
  for (const double invalid : {-0.1, 1.1, std::numeric_limits<double>::quiet_NaN(),
                               std::numeric_limits<double>::infinity()}) {
    std::vector<double> coordinates(process.GetdLIPSDim(), 0.5);
    coordinates[0] = invalid;
    gra::MEventWeightState rejected;
    REQUIRE(process.EventWeight(coordinates, rejected) == Approx(0.0));
    CHECK_FALSE(rejected.kinematics_ok);
    CHECK(rejected.technical_failure);
  }
  gra::MEventWeightState missing;
  CHECK(process.EventWeight({}, missing) == Approx(0.0));
  CHECK(missing.technical_failure);
  gra::MRandom random;
  random.SetSeed(314159);
  DimuonPoint(process, random);
  const auto &lts = process.state.lts;
  CHECK(gra::math::CheckEMC(lts.pbeam1 + lts.pbeam2 - lts.pfinal[1] - lts.pfinal[2] -
                           lts.decaytree[0].p4 - lts.decaytree[1].p4));
}

// Accepted and rejected physical events must clear the preceding forward state
TEST_CASE("event boundaries own both forward excitation sides", "[Process][ForwardPhaseSpace][state]") {
  auto generator = DimuonGenerator();
  auto &process = *generator->proc;
  const auto seed = [&]() {
    process.state.lts.excite1 = process.state.lts.excite2 = true;
    process.state.lts.decayforward1.legs.resize(1);
    process.state.lts.decayforward2.legs.resize(1);
  };
  const auto empty = [&]() {
    CHECK_FALSE(process.state.lts.excite1);
    CHECK_FALSE(process.state.lts.excite2);
    CHECK(process.state.lts.decayforward1.legs.empty());
    CHECK(process.state.lts.decayforward2.legs.empty());
  };
  seed();
  gra::MEventWeightState invalid;
  CHECK(process.EventWeight({}, invalid) == Approx(0.0));
  empty();
  seed();
  gra::MRandom random;
  random.SetSeed(314159);
  DimuonPoint(process, random);
  empty();
  HepMC3::GenEvent event(HepMC3::Units::GEV, HepMC3::Units::MM);
  REQUIRE(process.EventRecord(event));
  REQUIRE_FALSE(event.particles().empty());
  REQUIRE_FALSE(process.EventRecord(event));
  empty();
}

TEST_CASE("forward branch construction commits only a complete decay",
          "[Process][ForwardPhaseSpace][state]") {
  ProcessKinematicsProbe probe;
  gra::MDecayBranch branch;
  branch.p.pdg = 990001;
  branch.legs.resize(1);
  branch.legs.front().p.pdg = 990002;
  const gra::MDecayBranch original = branch;

  const gra::M4Vec parent(0.0, 0.0, 0.0, 1.0);
  gra::MParticle daughter;
  daughter.pdg = 211;
  REQUIRE_THROWS_AS(probe.ProbeForwardBranch({}, {daughter}, parent, branch),
                    std::invalid_argument);
  REQUIRE(branch.p.pdg == original.p.pdg);
  REQUIRE(branch.legs.size() == original.legs.size());
  REQUIRE(branch.legs.front().p.pdg == original.legs.front().p.pdg);

  const gra::M4Vec invalid_daughter(0.0, 0.0, 0.0, -1.0);
  REQUIRE_FALSE(
      probe.ProbeForwardBranch({invalid_daughter}, {daughter}, parent, branch));
  REQUIRE(branch.p.pdg == original.p.pdg);
  REQUIRE(branch.legs.size() == original.legs.size());
  REQUIRE(branch.legs.front().p.pdg == original.legs.front().p.pdg);
}

TEST_CASE("process parser accepts MadGraph names and generic antiparton alias",
          "[Process][Parser][MG5]") {
  gra::MPDG pdg;
  for (const int code : {12, -12, 14, -14, 16, -16, 89, -89}) {
    gra::MParticle particle;
    particle.pdg = code;
    pdg.PDG_table.emplace(code, particle);
  }

  REQUIRE(pdg.FindByPDGName("ve").pdg == 12);
  REQUIRE(pdg.FindByPDGName("ve~").pdg == -12);
  REQUIRE(pdg.FindByPDGName("vm").pdg == 14);
  REQUIRE(pdg.FindByPDGName("vm~").pdg == -14);
  REQUIRE(pdg.FindByPDGName("vt").pdg == 16);
  REQUIRE(pdg.FindByPDGName("vt~").pdg == -16);
  REQUIRE(pdg.FindByPDGName("~j").pdg == -89);
}

// Preserve tiny elastic momentum transfers in the physical production builder
TEST_CASE("elastic recoil resolves tiny t and rejects positive t", "[Process][PhaseSpace][Elastic]") {
  QuasiElasticKinematicsProbe probe;
  ConfigureAsymmetricProtonBeams(probe, 6500.0, 6500.0);
  const double mass2 = probe.state.lts.pbeam1.M2();
  for (const double t : {-1e-3, -1e-6, -1e-9, -1e-12}) {
    REQUIRE(probe.BuildBornForTest(mass2, mass2, t));
    REQUIRE(probe.state.lts.t == Approx(t).epsilon(1e-9).margin(1e-18));
    REQUIRE(gra::math::CheckEMC(probe.state.lts.pbeam1 + probe.state.lts.pbeam2 -
                               probe.state.lts.pfinal[1] - probe.state.lts.pfinal[2]));
  }
  REQUIRE_FALSE(probe.BuildBornForTest(mass2, mass2, 1.0));
  REQUIRE_FALSE(probe.BuildBornForTest(mass2, mass2, -probe.state.lts.s));
}

// Serialize actual soft scatterings and each terminal remnant with local momentum closure
TEST_CASE("soft event vertices conserve momentum through both remnant decays", "[Process][SoftChain][HepMC]") {
  QuasiElasticKinematicsProbe probe;
  ConfigureAsymmetricProtonBeams(probe, 100.0, 100.0);
  const int lower_pdg = GENERATE(2212, -2212);
  probe.state.lts.beam2.pdg = lower_pdg;
  probe.ProcPtr.CHANNEL = "ND";
  probe.SetModelTune(gra::MModelTune::Load(gra::aux::ResolveProjectPath("modeldata/TUNE0/GENERAL.json")));
  probe.state.lts.PDG.ReadParticleData();
  probe.SeedForTest(74531);
  probe.state.multipomeron_impact_parameter = 1.0;
  auto &chain = probe.state.multipomeron_chain;
  chain.resize(3);
  auto &point = chain[0];
  point.p1i = probe.state.lts.pbeam1;
  point.p2i = probe.state.lts.pbeam2;
  point.p1f = gra::M4Vec(0, 0, 60, std::sqrt(4500.0));
  point.p2f = gra::M4Vec(0, 0, -60, std::sqrt(4500.0));
  point.q1 = point.p1i - point.p1f;
  point.q2 = point.p2i - point.p2f;
  point.k = point.q1 + point.q2;
  chain[1].k = point.p1f;
  chain[2].k = point.p2f;
  // Terminal entries carry only k and must never be interpreted as scatterings
  chain[1].q1 = gra::M4Vec(1, 2, 3, 4);
  chain[2].p1f = gra::M4Vec(5, 6, 7, 8);
  HepMC3::GenEvent event(HepMC3::Units::GEV, HepMC3::Units::MM);
  REQUIRE(probe.BuildEventRecord(event));
  for (const auto &vertex : event.vertices()) {
    gra::M4Vec balance;
    for (const auto &particle : vertex->particles_in()) {
      const auto &p = particle->momentum();
      balance += gra::M4Vec(p.px(), p.py(), p.pz(), p.e());
    }
    for (const auto &particle : vertex->particles_out()) {
      const auto &p = particle->momentum();
      REQUIRE(p.e() > 0.0);
      balance -= gra::M4Vec(p.px(), p.py(), p.pz(), p.e());
    }
    REQUIRE(gra::math::CheckEMC(balance));
  }
  int charge = 0;
  for (const auto &particle : event.particles()) {
    if (particle->status() == gra::PDG::PDG_STABLE) {
      charge += probe.state.lts.PDG.FindByPDG(particle->pid()).chargeX3;
    }
  }
  REQUIRE(charge == (lower_pdg > 0 ? 6 : 0));
}
