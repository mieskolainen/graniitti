// QuasiElastic (EL,SD,DD) and soft ND (simplified) class
// with phase space class <Q>
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <exception>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

// HepMC3
#include "HepMC3/FourVector.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"

// Own
#include "Graniitti/Kinematics/MEventKinematics.h"
#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/Kinematics/MQuasiElastic.h"
#include "Graniitti/MUserCuts.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Particle/MForm.h"
#include "Graniitti/Process/MEventRecord.h"
#include "Graniitti/Process/MProcess.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MException.h"

using gra::aux::indices;
using gra::math::abs2;
using gra::math::msqrt;
using gra::math::PI;
using gra::math::pow2;
using gra::math::zi;
using gra::PDG::GeV2barn;
using gra::PDG::mp;

namespace gra {
// This is needed by construction
MQuasiElastic::MQuasiElastic() { Initialize(); }

// Constructor
MQuasiElastic::MQuasiElastic(std::string process, const std::vector<aux::OneCMD> &syntax, MModelTunePtr tune) {
  Initialize();
  InitHistograms();
  SetProcess(process, syntax, std::move(tune));

  // Init final states
  M4Vec zerovec(0, 0, 0, 0);
  state.lts.pfinal.assign(11, zerovec);

  // Pomeron weights
  MAXPOMW = std::vector<double>(100, 0.0);
  std::cout << "MQuasiElastic:: [Constructor done]" << std::endl;
}

void MQuasiElastic::Initialize() {
  const std::vector<std::string> supported = {"X"};
  state.phase_space_class                  = "Q";
  ProcPtr                                  = MSubProc(supported, state.phase_space_class);
}

// Compute true when one value is inside an inclusive fiducial range
bool PassQuasiElasticFiducialRange(double value, double min, double max) { return min <= value && value <= max; }

// Destructor
MQuasiElastic::~MQuasiElastic() {}

// Initialize cut and process spesific postsetup
void MQuasiElastic::FinalizeProcessConfiguration() {
  if (state.excitation != 0) {
    throw std::invalid_argument(
        "MQuasiElastic::FinalizeProcessConfiguration: X[EL/SD/DD/ND]<Q> encodes the forward topology and requires "
        "NSTARS=0");
  }

  // Set phase space dimension
  if (ProcPtr.CHANNEL == "EL") { ProcPtr.LIPSDIM = 1; };
  if (ProcPtr.CHANNEL == "SD") { ProcPtr.LIPSDIM = 2; };
  if (ProcPtr.CHANNEL == "DD") { ProcPtr.LIPSDIM = 3; };
  if (ProcPtr.CHANNEL == "ND") { ProcPtr.LIPSDIM = 1; };  // Keep it 1

  const bool elastic       = ProcPtr.ISTATE == "X" && ProcPtr.CHANNEL == "EL";
  const bool diffraction   = elastic || (ProcPtr.ISTATE == "X" && (ProcPtr.CHANNEL == "SD" || ProcPtr.CHANNEL == "DD"));
  const bool valid_minimum = elastic ? state.gcuts.q_t_abs_min > 0.0 : state.gcuts.q_t_abs_min >= 0.0;
  if (diffraction && (!std::isfinite(state.gcuts.q_t_abs_min) || !std::isfinite(state.gcuts.q_t_abs_max) ||
                      !valid_minimum || state.gcuts.q_t_abs_max <= state.gcuts.q_t_abs_min)) {
    throw std::invalid_argument(
        "MQuasiElastic::FinalizeProcessConfiguration: invalid generated "
        "absolute momentum-transfer range");
  }
  if (ProcPtr.ISTATE == "X" && (ProcPtr.CHANNEL == "SD" || ProcPtr.CHANNEL == "DD")) {
    aux::PrintWarning(false);
    std::cout << rang::style::bold << rang::fg::red << "Minbias X[" << ProcPtr.CHANNEL
              << "]<Q> is work-in-progress" << rang::style::reset << std::endl;
  }
}

// Fiducial user cuts
bool MQuasiElastic::FiducialCuts() const {
  if (state.fcuts.active == true) {
    if (!state.fcuts.forward_t1.Contains(std::abs(state.lts.t)) ||
        !state.fcuts.forward_t2.Contains(std::abs(state.lts.t))) { return false; }
    if (!state.fcuts.PassForwardDeltaPhi(state.lts)) { return false; }
    if (!state.fcuts.PassForwardXi(state.lts)) { return false; }

    if (state.fcuts.forward_t_active &&
        !PassQuasiElasticFiducialRange(std::abs(state.lts.t), state.fcuts.forward_t_min, state.fcuts.forward_t_max)) {
      return false;
    }

    // EL cuts
    if (ProcPtr.CHANNEL == "EL") {
      // no mass cuts for elastic final states
    }

    // SD cuts
    if (ProcPtr.CHANNEL == "SD") {
      if (state.fcuts.forward_M_active) {
        const M4Vec &excited = (state.lts.pfinal[1].M() > 1.0) ? state.lts.pfinal[1] : state.lts.pfinal[2];
        if (!PassQuasiElasticFiducialRange(excited.M(), state.fcuts.forward_M_min, state.fcuts.forward_M_max)) {
          return false;
        }
      }
    }

    // DD cuts
    if (ProcPtr.CHANNEL == "DD") {
      if (state.fcuts.forward_M_active) {
        if (!PassQuasiElasticFiducialRange(state.lts.pfinal[1].M(), state.fcuts.forward_M_min,
                                           state.fcuts.forward_M_max) ||
            !PassQuasiElasticFiducialRange(state.lts.pfinal[2].M(), state.fcuts.forward_M_min,
                                           state.fcuts.forward_M_max)) {
          return false;
        }
      }
    }

    // Check user cuts (do not substitute to kinematics =
    // UserCut...)
    if (!UserCut(state.usercuts, state.lts)) {
      return false;  // not fine
    }
  }
  return true;  // fine
}

bool MQuasiElastic::LoopKinematics(const std::array<double, 2> &p1p, const std::array<double, 2> &p2p) {
  if (!kinematics::RebuildScreeningKinematics(state.lts, p1p, p2p, false)) { return false; }
  return B3GetLorentzScalars();
}

// Recompute quasielastic invariants from the restored Born four-momenta
bool MQuasiElastic::RefreshBornKinematics() { return B3GetLorentzScalars(); }

// Compute the elastic momentum-table support required by the generation range
double MQuasiElastic::EikonalMaxKT2() const {
  if (ProcPtr.ISTATE == "X" && ProcPtr.CHANNEL == "EL") { return MaximumSampledAbsT(); }
  return 0.0;
}

// Access the configured minimum generated |t|
double MQuasiElastic::MinimumSampledAbsT() const { return state.gcuts.q_t_abs_min; }

// Access the configured maximum generated |t|
double MQuasiElastic::MaximumSampledAbsT() const { return state.gcuts.q_t_abs_max; }

// Sample |t| from equal inverse-variable and logarithmic proposals
// p(|t|) = [p_inverse(|t|)+p_log(|t|)]/2
double MQuasiElastic::SampleElasticAbsT(double unit, double min_abs_t, double max_abs_t) {
  if (unit < 0.5) {
    const double branch_unit   = 2.0 * unit;
    const double inverse_abs_t = (1.0 - branch_unit) / min_abs_t + branch_unit / max_abs_t;
    return 1.0 / inverse_abs_t;
  }
  const double branch_unit = 2.0 * unit - 1.0;
  return min_abs_t * std::exp(branch_unit * std::log(max_abs_t / min_abs_t));
}

// Compute the inverse density of the equal elastic proposal mixture
// 1/p(|t|), p = (p_inverse+p_log)/2
double MQuasiElastic::ElasticAbsTJacobian(double abs_t, double min_abs_t, double max_abs_t) {
  const double inverse_variable_density = 1.0 / (pow2(abs_t) * (1.0 / min_abs_t - 1.0 / max_abs_t));
  const double logarithmic_density      = 1.0 / (abs_t * std::log(max_abs_t / min_abs_t));
  return 1.0 / (0.5 * inverse_variable_density + 0.5 * logarithmic_density);
}

// Get weight
double MQuasiElastic::ComputeEventWeight(const std::vector<double> &randvec, MEventWeightState &aux) {
  double W = 0.0;

  if (ProcPtr.CHANNEL != "ND") {  // Diffractive
    PreparePhaseSpacePoint(B3RandomKin(randvec), aux);

    if (aux.Valid()) {
      // ** EVENT WEIGHT **
      const double LIPS   = B3PhaseSpaceWeight();  // Phase-space weight
      const double MatESQ = GetAmp2(aux.include_screening, aux);

      // Total weight: phase-space x |M|^2 x barn units
      W = DissociationCrossSectionFactor() * LIPS * B3IntegralVolume() * MatESQ * GeV2barn;
    }

  } else {  // Non-Diffractive

    if (!eikonal.IsInitialized()) {
      // This should be rejected during setup, but sampling must remain alive
      aux.kinematics_ok     = false;
      aux.technical_failure = true;
      return 0.0;
    }

    aux.fidcuts_ok    = true;
    aux.vetocuts_ok   = true;

    const bool kinematics_ok = BuildSoftChain();
    aux.kinematics_ok        = kinematics_ok;
    if (!kinematics_ok) { return 0.0; }

    const double LIPS   = B3PhaseSpaceWeight();
    const double MatESQ = GetAmp2(false, aux);

    // W = LIPS * MatESQ * GeV2barn;  // Total weight: phase-space x |M|^2 x
    // barn units
    W = DissociationCrossSectionFactor() * LIPS * MatESQ;  // Total weight:

    // Bypass the outer rejection step for the internally sampled cut-Pomeron state
    aux.forced_accept = W > 0.0;
  }

  return W;
}

// Record HepMC3 event
bool MQuasiElastic::BuildEventRecord(HepMC3::GenEvent &evt) {
  if (ProcPtr.CHANNEL == "ND") { return BuildSoftEventRecord(evt); }

  // Diffractive processes

  // Initial states (4-momentum, pdg-id, status code)
  HepMC3::GenParticlePtr gen_p1 = std::make_shared<HepMC3::GenParticle>(gra::aux::M4Vec2HepMC3(state.lts.pbeam1),
                                                                        state.lts.beam1.pdg, PDG::PDG_BEAM);
  HepMC3::GenParticlePtr gen_p2 = std::make_shared<HepMC3::GenParticle>(gra::aux::M4Vec2HepMC3(state.lts.pbeam2),
                                                                        state.lts.beam2.pdg, PDG::PDG_BEAM);

  // Propagator 4-vector and generator particle
  M4Vec                  q1(state.lts.pbeam1 - state.lts.pfinal[1]);
  HepMC3::GenParticlePtr gen_q1 =
      std::make_shared<HepMC3::GenParticle>(gra::aux::M4Vec2HepMC3(q1), PDG::PDG_propagator, PDG::PDG_INTERMEDIATE);

  // Final state protons/excited systems
  int PDG_ID1     = state.lts.beam1.pdg;
  int PDG_ID2     = state.lts.beam2.pdg;
  int PDG_status1 = PDG::PDG_STABLE;
  int PDG_status2 = PDG::PDG_STABLE;

  // EL
  // already fine

  // SD
  if (ProcPtr.CHANNEL == "SD") {
    if (state.lts.excite1) {  // proton 1 excited
      PDG_ID1     = std::abs(PDG::PDG_NSTAR) * math::sign(state.lts.beam1.pdg);
      PDG_status1 = PDG::PDG_INTERMEDIATE;
    } else {
      PDG_ID2     = std::abs(PDG::PDG_NSTAR) * math::sign(state.lts.beam2.pdg);
      PDG_status2 = PDG::PDG_INTERMEDIATE;
    }
  }
  // DD
  if (ProcPtr.CHANNEL == "DD") {
    PDG_ID1     = std::abs(PDG::PDG_NSTAR) * math::sign(state.lts.beam1.pdg);
    PDG_status1 = PDG::PDG_INTERMEDIATE;
    PDG_ID2     = std::abs(PDG::PDG_NSTAR) * math::sign(state.lts.beam2.pdg);
    PDG_status2 = PDG::PDG_INTERMEDIATE;
  }

  HepMC3::GenParticlePtr gen_p1f =
      std::make_shared<HepMC3::GenParticle>(gra::aux::M4Vec2HepMC3(state.lts.pfinal[1]), PDG_ID1, PDG_status1);
  HepMC3::GenParticlePtr gen_p2f =
      std::make_shared<HepMC3::GenParticle>(gra::aux::M4Vec2HepMC3(state.lts.pfinal[2]), PDG_ID2, PDG_status2);

  // Construct vertices

  // Upper proton-pomeron-proton
  HepMC3::GenVertexPtr v1 = std::make_shared<HepMC3::GenVertex>();
  v1->add_particle_in(gen_p1);
  v1->add_particle_out(gen_p1f);
  v1->add_particle_out(gen_q1);

  // Lower proton-pomeron-proton
  HepMC3::GenVertexPtr v2 = std::make_shared<HepMC3::GenVertex>();
  v2->add_particle_in(gen_p2);
  v2->add_particle_out(gen_p2f);
  v2->add_particle_in(gen_q1);

  // Finally add all vertices
  evt.add_vertex(v1);
  evt.add_vertex(v2);

  // Fragment every excited diffractive side through the common N-star
  // dispatcher
  if (!CEPForwardFragment()) { throw PhaseSpaceFailure("MQuasiElastic::BuildEventRecord: forward excitation failed"); }
  if (state.lts.excite1 && !state.lts.decayforward1.legs.empty()) {
    record::WriteBranch(state.lts.decayforward1, gen_p1f, evt, state.random, false);
  }
  if (state.lts.excite2 && !state.lts.decayforward2.legs.empty()) {
    record::WriteBranch(state.lts.decayforward2, gen_p2f, evt, state.random, false);
  }
  return true;
}

// Fragment one positive-energy soft system at its existing production particle
bool MQuasiElastic::WriteSoftDecay(const HepMC3::GenParticlePtr &particle, const M4Vec &momentum,
                                    int baryon, int charge, HepMC3::GenEvent &evt) {
  if (!(momentum.E() > 0.0) || !std::isfinite(momentum.M2()) || !(momentum.M2() > 0.0)) { return false; }
  MDecayBranch branch;
  if (!ExciteContinuum(momentum, branch, momentum.M2(), baryon, charge)) { return false; }
  record::WriteBranch(branch, particle, evt, state.random, false);
  return true;
}

// Write the cut Pomeron chain and decay its two final remnants independently
bool MQuasiElastic::BuildSoftEventRecord(HepMC3::GenEvent &evt) {
  const auto &chain = state.multipomeron_chain;
  if (chain.size() < 3) { return false; }
  const std::size_t count = chain.size() - 2;
  const auto particle = [](const M4Vec &p, int pdg, int status) {
    return std::make_shared<HepMC3::GenParticle>(aux::M4Vec2HepMC3(p), pdg, status);
  };
  auto upper = particle(chain[0].p1i, state.lts.beam1.pdg, PDG::PDG_BEAM);
  auto lower = particle(chain[0].p2i, state.lts.beam2.pdg, PDG::PDG_BEAM);
  for (std::size_t i = 0; i < count; ++i) {
    const auto &point = chain[i];
    if (!(point.k.E() > 0.0) ||
        !math::CheckEMC(point.p1i - point.p1f - point.q1) ||
        !math::CheckEMC(point.p2i - point.p2f - point.q2) ||
        !math::CheckEMC(point.q1 + point.q2 - point.k)) { return false; }
    if (i > 0 && (!math::CheckEMC(chain[i - 1].p1f - point.p1i) ||
                  !math::CheckEMC(chain[i - 1].p2f - point.p2i))) { return false; }
    auto upper_out = particle(point.p1f, PDG::PDG_fragment, PDG::PDG_INTERMEDIATE);
    auto lower_out = particle(point.p2f, PDG::PDG_fragment, PDG::PDG_INTERMEDIATE);
    auto q1 = particle(point.q1, PDG::PDG_propagator, PDG::PDG_INTERMEDIATE);
    auto q2 = particle(point.q2, PDG::PDG_propagator, PDG::PDG_INTERMEDIATE);
    auto system = particle(point.k, PDG::PDG_system, PDG::PDG_INTERMEDIATE);
    auto v1 = std::make_shared<HepMC3::GenVertex>();
    v1->add_particle_in(upper);
    v1->add_particle_out(upper_out);
    v1->add_particle_out(q1);
    auto v2 = std::make_shared<HepMC3::GenVertex>();
    v2->add_particle_in(lower);
    v2->add_particle_out(lower_out);
    v2->add_particle_out(q2);
    auto vx = std::make_shared<HepMC3::GenVertex>();
    vx->add_particle_in(q1);
    vx->add_particle_in(q2);
    vx->add_particle_out(system);
    evt.add_vertex(v1);
    evt.add_vertex(v2);
    evt.add_vertex(vx);
    if (!WriteSoftDecay(system, point.k, 0, 0, evt)) { return false; }
    upper = upper_out;
    lower = lower_out;
  }
  if (!math::CheckEMC(chain[count].k - chain[count - 1].p1f) ||
      !math::CheckEMC(chain[count + 1].k - chain[count - 1].p2f)) { return false; }
  const int b1 = math::sign(state.lts.beam1.pdg);
  const int b2 = math::sign(state.lts.beam2.pdg);
  return WriteSoftDecay(upper, chain[count].k, b1, b1, evt) &&
         WriteSoftDecay(lower, chain[count + 1].k, b2, b2, evt);
}

// Print out setup
void MQuasiElastic::PrintInit(bool silent) const {
  if (!silent) {
    PrintSetup();

    // Diffractive processes
    if (ProcPtr.CHANNEL != "ND") {
      std::string proton1 = "-----------EL--------->";
      std::string proton2 = "-----------EL--------->";

      if (ProcPtr.CHANNEL == "SD") { proton1 = "-----------F2-xxxxxxxx>"; }
      if (ProcPtr.CHANNEL == "DD") {
        proton1 = "-----------F2-xxxxxxxx>";
        proton2 = "-----------F2-xxxxxxxx>";
      }

      std::vector<std::string> feynmangraph;
      feynmangraph = {"||          ", "||          ", "||          ", "||          ", "||          "};

      // Print diagram
      std::cout << proton1 << std::endl;
      for (const auto &i : indices(feynmangraph)) {
        if (state.screening) {  // Put red
          std::cout << rang::fg::red << "     **    " << rang::style::reset;
        } else {
          std::cout << rang::fg::red << "           " << rang::style::reset;
        }
        std::cout << feynmangraph[i] << std::endl;
      }
      std::cout << proton2 << std::endl << std::endl;

      // Generation cuts
      std::cout << rang::style::bold << "Generation cuts:" << rang::style::reset << std::endl << std::endl;

      if (ProcPtr.CHANNEL != "EL") {
        printf("- Xi  [min, max] = [%0.3E, %0.3E] (Xi == M^2/s) \n", state.gcuts.XI_min, state.gcuts.XI_max);
      }
      printf("- |t| [min, max] = [%0.3f, %0.3f] GeV^2 \n", MinimumSampledAbsT(), MaximumSampledAbsT());
      std::cout << std::endl;

      // Fiducial cuts
      std::cout << rang::style::bold << "Fiducial cuts:" << rang::style::reset << std::endl << std::endl;
      if (state.fcuts.active) {
        if (state.fcuts.HasForwardCuts()) {
          if (state.fcuts.forward_M_active) {
            printf("- M   [min, max]   = [%0.2f, %0.2f] GeV^2 \n", state.fcuts.forward_M_min,
                   state.fcuts.forward_M_max);
          }
          if (state.fcuts.forward_t_active) {
            printf("- |t| [min, max]   = [%0.2f, %0.2f] GeV^2 \n", state.fcuts.forward_t_min,
                   state.fcuts.forward_t_max);
          }
          if (state.fcuts.forward_xi_active) {
            printf("- Xi_i [min, max]  = [%0.3E, %0.3E] \n", state.fcuts.forward_xi_min, state.fcuts.forward_xi_max);
          }
          if (state.fcuts.forward_dPhi_active) {
            printf("- dPhi [min, max]  = [%0.2f, %0.2f] deg \n", state.fcuts.forward_dPhi_min,
                   state.fcuts.forward_dPhi_max);
          }
        } else {
          std::cout << "- no cuts" << std::endl;
        }
        std::cout << std::endl;
      } else {
        std::cout << "- Not active" << std::endl;
      }

      // Non-Diffractive
    } else {
      for (std::size_t i = 0; i < 7; ++i) {
        if (i > 0 && i < 6) {
          if (state.screening) {
            std::cout << rang::fg::red << "     **    " << rang::style::reset;
          } else {
            std::cout << rang::fg::red << "           " << rang::style::reset;
          }
        } else {
          std::cout << "-----------";
        }
        std::cout << "|x|--------->" << std::endl;
      }
    }
    std::cout << std::endl;
  }
}

// N semi-independent (soft) Pomeron exchanges, 4-momentum conserved
//
// ------------------------------------> remnant 1
//      |      |       |          |
//      -- s0  |       -- s2      -- s(N-1)
// s    |      |       |          |
//      |      -- s1   |     ...  |
//      |      |       |          |
// ------------------------------------> remnant 2
//
//
// Beam fragmentation uses Regge-like momentum fractions without explicit color degrees of freedom
// Multiplicity shapes depend on the fragmentation prescription and Poisson eikonal ansatz
// Their energy evolution requires validation against several observables
bool MQuasiElastic::BuildSoftChain() {
  // Boundary conditions
  const double XMIN = 1e-6;
  const double XMAX = 1.0;

  const double MMIN       = 3.0;
  const double REM_M2_MIN = pow2(10.0);
  const double DELTA      = 0.10;

  int       outertrials = 0;
  const int MAXTRIAL    = 1e4;

  // Sample the number of inelastic cut Pomerons
  unsigned int N = 0;

  try {
    eikonal.S3GetRandomCutsBt(N, state.multipomeron_impact_parameter, state.random.rng);
  } catch (const std::exception &error) {
    throw PhaseSpaceFailure("MQuasiElastic::BuildSoftChain: cut-Pomeron sampling failed: " +
                            std::string(error.what()));
  }
  if (N == 0) { throw PhaseSpaceFailure("MQuasiElastic::BuildSoftChain: no cut Pomeron was sampled"); }

  state.multipomeron_chain.assign(N + 2, MultipomeronKinematics{});

  // Beams in
  state.multipomeron_chain[0].p1i = state.lts.pbeam1;  // @@@@@@
  state.multipomeron_chain[0].p2i = state.lts.pbeam2;  // @@@@@@

  while (true) {
    bool faulty = false;
    for (std::size_t i = 0; i < N; ++i) {
      const double s = (state.multipomeron_chain[i].p1i + state.multipomeron_chain[i].p2i).M2();

      // Draw Bjorken-x
      double x1 = 0;
      double x2 = 0;

      const double p1_m2 = state.multipomeron_chain[i].p1i.M2();
      const double p2_m2 = state.multipomeron_chain[i].p2i.M2();
      if (p1_m2 < 0.0 || p2_m2 < 0.0 || s <= 0.0) {
        faulty = true;
        break;
      }

      const double sqrt_s   = msqrt(s);
      const double m1       = msqrt(p1_m2);
      const double m2       = msqrt(p2_m2);
      const double pnorm_cm = gra::kinematics::DecayMomentum(sqrt_s, m1, m2);
      if (!(pnorm_cm > 0.0)) {
        faulty = true;
        break;
      }


      // No Q^2 dependence taken into account here
      // We assume xG(x) ~ 1/x^{DELTA} <=> G(x) = 1/x^{1+DELTA}
      x1 = state.random.PowerBoundedRandom(XMIN, XMAX, DELTA);
      x2 = state.random.PowerBoundedRandom(XMIN, XMAX, DELTA);

      // Pick gaussian Fermi pt of Pomerons
      const double sigma = 0.4;  // GeV
      const double pt1   = msqrt(pow2(state.random.G(0, sigma)) + pow2(state.random.G(0, sigma)));
      const double pt2   = msqrt(pow2(state.random.G(0, sigma)) + pow2(state.random.G(0, sigma)));
      const double phi1  = state.random.U(0, 2.0 * PI);
      const double phi2  = state.random.U(0, 2.0 * PI);

      // Pick remnant mass^2
      double p3_m2 = p1_m2;
      double p4_m2 = p2_m2;

      if (i == N - 1) {  // Excite remnants (one could excite also intermediate)
        // Max operator for low energies
        const double p1_E_cm = 0.5 * (s + p1_m2 - p2_m2) / sqrt_s;
        const double p2_E_cm = 0.5 * (s + p2_m2 - p1_m2) / sqrt_s;
        if (pow2(p1_E_cm) > REM_M2_MIN) {
          p3_m2 = state.random.PowerBoundedRandom(REM_M2_MIN, pow2(p1_E_cm), (1 + DELTA));
        }
        if (pow2(p2_E_cm) > REM_M2_MIN) {
          p4_m2 = state.random.PowerBoundedRandom(REM_M2_MIN, pow2(p2_E_cm), (1 + DELTA));
        }
      }

      // Preserve the incoming remnant directions and positive production energy
      auto &chain = state.multipomeron_chain[i];
      const std::array<M4Vec, 2> transverse = {
          M4Vec(pt1 * std::cos(phi1), pt1 * std::sin(phi1), 0.0, 0.0),
          M4Vec(pt2 * std::cos(phi2), pt2 * std::sin(phi2), 0.0, 0.0)};
      std::array<M4Vec, 2> remnant;
      if (!kinematics::BuildSoftRemnants({chain.p1i, chain.p2i}, {x1, x2},
                                          {p3_m2, p4_m2}, transverse, remnant, chain.k) ||
          chain.k.M() < MMIN) {
        faulty = true;
        break;
      }
      chain.p1f = remnant[0];
      chain.p2f = remnant[1];
      chain.q1 = chain.p1i - chain.p1f;
      chain.q2 = chain.p2i - chain.p2f;

      // Remnants out
      if (i < N - 1) {
        state.multipomeron_chain[i + 1].p1i = state.multipomeron_chain[i].p1f;  // @@@@@@
        state.multipomeron_chain[i + 1].p2i = state.multipomeron_chain[i].p2f;  // @@@@@@
      }
    }

    // Take forward and backward remnants energy-momentum
    state.multipomeron_chain[N].k     = state.multipomeron_chain[N - 1].p1f;
    state.multipomeron_chain[N + 1].k = state.multipomeron_chain[N - 1].p2f;

    // Check Energy-Momentum Conservation here explicitly !!!
    if (!faulty) {
      M4Vec p_sum;
      for (const auto &i : indices(state.multipomeron_chain)) {
        p_sum += state.multipomeron_chain[i].k;
      }
      if (!gra::math::CheckEMC(p_sum - (state.lts.pbeam1 + state.lts.pbeam2))) {
        faulty = true;
      }
    }

    if (faulty) {
      ++outertrials;
      if (outertrials > MAXTRIAL) {
        throw PhaseSpaceFailure("MQuasiElastic::BuildSoftChain: failed to construct a physical soft chain");
      }
      continue;
    } else {
      break;
    }
  }

  return true;
}

// 3-dimensional phase space vector initialization
// The adaptive phase space boundaries here are taken into account by
// B3IntegralVolume() function
bool MQuasiElastic::B3RandomKin(const std::vector<double> &randvec) {
  t_sampling_jacobian = 0.0;
  state.lts.excite1   = false;
  state.lts.excite2   = false;

  // Elastic case: s1 = s3, s2 = s4
  double s3 = state.lts.pbeam1.M2();
  double s4 = state.lts.pbeam2.M2();

  // Set Diffractive mass boundaries
  // neutron + piplus + safe margin
  M2_f_min = std::max(state.gcuts.XI_min * state.lts.s, gra::math::pow2(1.3));
  M2_f_max = state.gcuts.XI_max * state.lts.s;

  if (ProcPtr.CHANNEL != "EL" && M2_f_max <= M2_f_min) { return false; }

  log_M2_f_min = std::log(M2_f_min);
  log_M2_f_max = (M2_f_max > 0.0) ? std::log(M2_f_max) : 0.0;

  // Sample diffractive system masses
  if (ProcPtr.CHANNEL == "SD") {
    // Log-change of variable
    const double u = log_M2_f_min + (log_M2_f_max - log_M2_f_min) * randvec[1];
    const double r = std::exp(u);

    // Choose state.random permutation
    if (state.random.U(0, 1) < 0.5) {
      s3                = r;
      state.lts.excite1 = true;
      state.lts.excite2 = false;
    } else {
      s4                = r;
      state.lts.excite1 = false;
      state.lts.excite2 = true;
    }
  } else if (ProcPtr.CHANNEL == "DD") {
    // Product bound M_1^2 M_2^2 <= xi_max m_p^2 s defines a triangle in log
    // masses
    const double DD_M2_product_max = state.gcuts.XI_max * pow2(mp) * state.lts.s;

    if (DD_M2_product_max <= pow2(M2_f_min)) { return false; }

    DD_M2_1_max     = std::min(M2_f_max, DD_M2_product_max / M2_f_min);
    log_DD_M2_1_max = std::log(DD_M2_1_max);

    // Log-change of variable
    const double u1 = log_M2_f_min + (log_DD_M2_1_max - log_M2_f_min) * randvec[1];
    const double r1 = std::exp(u1);

    DD_M2_max     = std::min(M2_f_max, DD_M2_product_max / r1);
    log_DD_M2_max = std::log(DD_M2_max);
    if (DD_M2_max <= M2_f_min) { return false; }

    // Log-change of variable
    const double u2 = log_M2_f_min + (log_DD_M2_max - log_M2_f_min) * randvec[2];
    const double r2 = std::exp(u2);

    // Choose state.random permutations
    if (state.random.U(0, 1) < 0.5) {
      s3 = r1;
      s4 = r2;
    } else {
      s4 = r1;
      s3 = r2;
    }
    state.lts.excite1 = true;
    state.lts.excite2 = true;
  }

  // Calculate kinematically valid t-range
  const double s1 = state.lts.pbeam1.M2();
  const double s2 = state.lts.pbeam2.M2();

  // Mandelstam t-range calculation
  gra::kinematics::Two2TwoLimit(state.lts.s, s1, s2, s3, s4, t_min, t_max);

  // Then limit the sampled absolute momentum-transfer range
  t_min = std::max(-MaximumSampledAbsT(), t_min);
  t_max = std::min(-MinimumSampledAbsT(), t_max);
  if (t_min > t_max) { return false; }

  // Logarithmic change of variable sampling
  const double A = std::abs(t_max);
  const double B = std::abs(t_min);

  double abs_t = 0.0;
  if (ProcPtr.CHANNEL == "EL") {
    if (!(A > 0.0) || !(B > A) || !std::isfinite(A) || !std::isfinite(B)) { return false; }
    abs_t               = SampleElasticAbsT(randvec[0], A, B);
    t_sampling_jacobian = ElasticAbsTJacobian(abs_t, A, B);
  } else {
    const double r = std::log(A + ZERO_EPS) + (std::log(B + ZERO_EPS) - std::log(A + ZERO_EPS)) * randvec[0];
    const double shifted_t = std::exp(r);
    abs_t               = shifted_t - ZERO_EPS;
    t_sampling_jacobian = (std::log(B + ZERO_EPS) - std::log(A + ZERO_EPS)) * shifted_t;
  }
  if (!(t_sampling_jacobian > 0.0) || !std::isfinite(t_sampling_jacobian)) { return false; }

  return B3BuildKin(s3, s4, -abs_t);
}

// Build kinematics for elastic, single and double diffractive 2->2 quasielastic
bool MQuasiElastic::B3BuildKin(double s3, double s4, double t) {
  if (!std::isfinite(s3) || !std::isfinite(s4) || !std::isfinite(t) || s3 < 0.0 || s4 < 0.0 ||
      state.lts.pfinal.size() <= 2 || !std::isfinite(state.lts.s) || !(state.lts.s > 0.0) ||
      !std::isfinite(state.lts.sqrt_s) || !(state.lts.sqrt_s > 0.0)) {
    return false;
  }
  const double s1      = state.lts.pbeam1.M2();
  const double s2      = state.lts.pbeam2.M2();
  const M4Vec  beamsum = state.lts.pbeam1 + state.lts.pbeam2;

  // Resolve the transverse recoil directly from t, including the forward limit
  double pz = 0.0;
  double pt = 0.0;
  if (!kinematics::ForwardScattering(state.lts.s, s1, s2, s3, s4, t, pz, pt)) { return false; }
  const double phi = state.random.U(0.0, 2.0 * PI);
  M4Vec p3(pt * std::cos(phi), pt * std::sin(phi), pz,
            0.5 * (state.lts.s + s3 - s4) / state.lts.sqrt_s);
  M4Vec p4(-p3.Px(), -p3.Py(), -pz,
            0.5 * (state.lts.s + s4 - s3) / state.lts.sqrt_s);

  // ------------------------------------------------------------------
  // Now boost if asymmetric beams
  if (std::abs(beamsum.Pz()) > 1e-6) {
    constexpr int sign = 1;  // positive -> boost to the lab
    kinematics::LorentzBoost(beamsum, state.lts.sqrt_s, p3, sign);
    kinematics::LorentzBoost(beamsum, state.lts.sqrt_s, p4, sign);
  }
  // ------------------------------------------------------------------

  state.lts.pfinal[1]     = p3;
  state.lts.pfinal[2]     = p4;
  state.lts.forward_mass2 = {s3, s4};

  // Check Energy-Momentum conservation
  if (!gra::math::CheckEMC(beamsum - (state.lts.pfinal[1] + state.lts.pfinal[2]))) { return false; }

  return B3GetLorentzScalars();
}

// Build and check scalars
bool MQuasiElastic::B3GetLorentzScalars() {
  if (state.lts.pfinal.size() <= 2) { return false; }

  // Calculate Lorentz scalars
  state.lts.ss[1][1] = state.lts.pfinal[1].M2();
  state.lts.ss[2][2] = state.lts.pfinal[2].M2();

  state.lts.t = (state.lts.pbeam1 - state.lts.pfinal[1]).M2();
  state.lts.u = (state.lts.pbeam1 - state.lts.pfinal[2]).M2();

  // Exact exchange and forward-system losses coincide for quasielastic
  // kinematics
  state.lts.x1      = kinematics::LongitudinalMomentumLoss(state.lts.pbeam1, state.lts.pfinal[1], true);
  state.lts.x2      = kinematics::LongitudinalMomentumLoss(state.lts.pbeam2, state.lts.pfinal[2], false);
  state.lts.xi1     = state.lts.x1;
  state.lts.xi2     = state.lts.x2;
  state.lts.has_xi1 = true;
  state.lts.has_xi2 = true;

  // Test scalars
  const std::array<double, 8> scalars = {state.lts.ss[1][1], state.lts.ss[2][2], state.lts.s,  state.lts.t,
                                         state.lts.u,        state.lts.x1,       state.lts.x2, state.lts.sqrt_s};
  if (!std::all_of(scalars.begin(), scalars.end(), [](const double value) { return std::isfinite(value); }) ||
      !(state.lts.s > 0.0) || !(state.lts.sqrt_s > 0.0) || state.lts.ss[1][1] < 0.0 || state.lts.ss[2][2] < 0.0 ||
      state.lts.ss[1][1] > state.lts.s || state.lts.ss[2][2] > state.lts.s || state.lts.t > 0.0 || state.lts.u > 0.0) {
    return false;
  }
  /*
  printf("s = %E, s^1/2 = %E \n", state.lts.s, msqrt(state.lts.s));
  printf("t = %E, u = %E \n", state.lts.t, state.lts.u);
  printf("s3 = %E, s4 = %E \n", state.lts.ss[1][1], state.lts.ss[2][2]);
  printf("\n");
  */
  return true;
}

// 1/2/3-Dim Integral Volume [t] x [M^2] x [M^2]
//
double MQuasiElastic::B3IntegralVolume() const {
  if (ProcPtr.CHANNEL == "EL") {
    return t_sampling_jacobian;

  } else if (ProcPtr.CHANNEL == "SD") {
    // Jacobian from log-change of variable:
    // \int_a^b f(M2) dM2 = \int_{ln(a)}^{ln(b)} f(exp(u)) * exp(u) du, where u
    // = ln(M2)
    const double J      = (state.lts.excite1) ? state.lts.ss[1][1] : state.lts.ss[2][2];
    const double M2_VOL = (log_M2_f_max - log_M2_f_min) * J;
    return t_sampling_jacobian * M2_VOL;

  } else if (ProcPtr.CHANNEL == "DD") {
    // Jacobian from log-change of variables
    const double J      = state.lts.ss[1][1] * state.lts.ss[2][2];
    const double M2_VOL = (log_DD_M2_1_max - log_M2_f_min) * (log_DD_M2_max - log_M2_f_min) * J;
    return t_sampling_jacobian * M2_VOL;

  } else {
    return 0;
  }
}

// Compute the standard two-body scattering phase-space and flux factor
// dPhi_2/flux = 1/[16 pi s^2 beta12(s)^2] per dt
//
double MQuasiElastic::B3PhaseSpaceWeight() const {
  // expression -> 16 * pi * s^2 (if s >> m1,m2)
  const double norm =
      16.0 * gra::math::PI *
      pow2(state.lts.s * gra::kinematics::beta12(state.lts.s, state.lts.beam1.mass, state.lts.beam2.mass));

  if (ProcPtr.CHANNEL == "EL") {
    return 1.0 / norm;
  } else if (ProcPtr.CHANNEL == "SD") {
    return 1.0 / norm;
  } else if (ProcPtr.CHANNEL == "DD") {
    return 1.0 / norm;
  } else if (ProcPtr.CHANNEL == "ND") {
    return 1.0;
  } else {
    return 0;
  }
}

}  // namespace gra
