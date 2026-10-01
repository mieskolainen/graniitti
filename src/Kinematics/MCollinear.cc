// Collinear parton model 2->2 or (2->N) type phase space class <P>
// for calculations without proton (remnant) considerations
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <algorithm>
#include <complex>
#include <iostream>
#include <random>
#include <vector>

// Own
#include "Graniitti/Kinematics/MCollinear.h"
#include "Graniitti/Kinematics/MEventKinematics.h"
#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/MUserCuts.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Particle/MForm.h"
#include "Graniitti/Process/MProcess.h"
#include "Graniitti/QCD/MPartonProposal.h"
#include "Graniitti/Regge/MFragment.h"
#include "Graniitti/Spin/MHelicity.h"
#include "Graniitti/Tech/MAux.h"

// Libraries
#include "rang.hpp"

using gra::aux::indices;

using gra::math::abs2;
using gra::math::CheckEMC;
using gra::math::msqrt;
using gra::math::PI;
using gra::math::pow2;
using gra::math::pow3;
using gra::math::pow4;
using gra::math::pow5;
using gra::math::zi;

using gra::PDG::GeV2barn;

namespace gra {
// This is needed by construction
MCollinear::MCollinear() { Initialize(); }

// Constructor
MCollinear::MCollinear(std::string process, const std::vector<aux::OneCMD> &syntax, MModelTunePtr tune) {
  Initialize();
  InitHistograms();
  SetProcess(process, syntax, std::move(tune));

  // Init final states
  M4Vec zerovec(0, 0, 0, 0);
  state.lts.pfinal.assign(11, zerovec);

  std::cout << "MCollinear:: [Constructor done]" << std::endl;
}

// Initialize the collinear central mass and rapidity proposal
void MCollinear::Initialize() {
  const std::vector<std::string> supported = {"yy_LUX", "yy_DZ"};
  state.phase_space_class                  = "P";
  ProcPtr                                  = MSubProc(supported, state.phase_space_class);
  state.lts.central_phase_space_mode       = CentralPhaseSpaceMode::Collinear;
}

// Destructor
MCollinear::~MCollinear() {}

// Map one unit coordinate logarithmically onto the central mass squared
// M2 = Mmin^2 exp[u log(Mmax^2/Mmin^2)]
double MCollinear::SampleCentralMassSquared(double unit, double mass_min, double mass_max) {
  const double log_mass2_min   = 2.0 * std::log(mass_min);
  const double log_mass2_range = 2.0 * (std::log(mass_max) - std::log(mass_min));
  return std::exp(log_mass2_min + unit * log_mass2_range);
}

// Compute the local Jacobian of the logarithmic central-mass map
// dM2/du = M2 log(Mmax^2/Mmin^2)
double MCollinear::CentralMassSquaredJacobian(double mass_squared, double mass_min, double mass_max) {
  const double log_mass2_range = 2.0 * (std::log(mass_max) - std::log(mass_min));
  return mass_squared * log_mass2_range;
}

// Initialize cut and process spesific postsetup
void MCollinear::FinalizeProcessConfiguration() {
  if (state.lts.decaytree.size() > 8) {
    throw std::invalid_argument(
        "MCollinear::FinalizeProcessConfiguration: direct central "
        "multiplicity cannot exceed 8");
  }

  if (state.screening) {
    throw std::invalid_argument(
        "MCollinear::FinalizeProcessConfiguration: LOOPSCREEN is not supported for collinear photon fluxes");
  }

  // Reject dissociation when the selected flux has no physical remnant state
  if (state.excitation != 0) {
    if (ProcPtr.ISTATE == "yy_DZ") {
      throw std::invalid_argument(
          "MCollinear::FinalizeProcessConfiguration: yy_DZ is an elastic proton flux and requires NSTARS=0; use "
          "yy<F> for physical proton dissociation");
    }
    throw std::invalid_argument(
        "MCollinear::FinalizeProcessConfiguration: yy_LUX is an inclusive photon PDF without event-level remnant "
        "kinematics and requires NSTARS=0; use yy<F> for physical proton dissociation");
  }

  if (ProcPtr.ISTATE == "MP" && ProcPtr.CHANNEL == "RES") {
    // Here we support only single resonances
    if (state.lts.process.RESONANCES.size() != 1) {
      std::string str =
          "MCollinear::FinalizeProcessConfiguration: Only single "
          "resonance supported for this "
          "process (RESPARAM.size() != 1)";
      throw std::invalid_argument(str);
    }
  }

  // Set sampling boundaries
  SetTechnicalBoundaries(state.gcuts, state.excitation);

  // Initialize phase space dimension
  ProcPtr.LIPSDIM = 2;  // All processes

  // Not applicable here
}

// No screening loop kinematics considered here
bool MCollinear::LoopKinematics(const std::array<double, 2> &p1p, const std::array<double, 2> &p2p) { return false; }

// Compute Monte Carlo integrand weight
double MCollinear::ComputeEventWeight(const std::vector<double> &randvec, MEventWeightState &aux) {
  double W = 0.0;

  if (parton::PrepareFinalStateProposal(state.lts, state.random)) { CalculateSymmetryFactor(); }

  // Kinematics and cuts
  PreparePhaseSpacePoint(B2RandomKin(randvec), aux);

  if (aux.Valid()) {
    // Matrix element squared
    const double MatESQ = GetAmp2(aux.include_screening, aux);

    // Calculate central system Phase Space volume
    double exact = 0.0;
    DecayWidthPS(exact);
    if (!aux.adaptation_mode) {
      state.lts.DW_sum_exact.AddLogWeight(gra::kinematics::MCW(exact), aux.log_inverse_density);

      // Add to the integral sum and take into account the proposal weight
      state.lts.DW_sum.AddLogWeight(state.lts.DW, aux.log_inverse_density);
    }

    double C_space = 1.0;
    // We have some legs in the central system
    if (state.lts.decaytree.size() != 0 && state.lts.PS_active) {
      C_space = state.lts.DW.Integral();

      // --------------------------------------------------------------------
      // Cascade resonances phase-space
      C_space *= CascadePS();
      // --------------------------------------------------------------------
    }

    // ** EVENT WEIGHT **
    W = DissociationCrossSectionFactor() * DecaySymmetryCompensationFactor() * C_space * (1.0 / state.symmetry_factor) *
        MatESQ * B2IntegralVolume() * B2PhaseSpaceWeight() * parton::ProposalWeight(state.lts) * GeV2barn /
        kinematics::MollerFlux(state.lts.pbeam1, state.lts.pbeam2);
  }

  return W;
}

void MCollinear::PrintInit(bool silent) const {
  if (!silent) {
    PrintSetup();

    // Construct prettyprint diagram
    std::string proton1 = "-----------pdf-------->";
    std::string proton2 = "-----------pdf-------->";

    /*
    if (state.excitation == 1) {
                    proton1 = "-----------F2-xxxxxxxx>";
    }
    if (state.excitation == 2) {
                    proton1 = "-----------F2-xxxxxxxx>";
                    proton2 = "-----------F2-xxxxxxxx>";
    }
    */

    std::vector<std::string> feynmangraph;
    feynmangraph = {"||          ", "||          ", "xx--------->", "||          ", "||          "};

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
    std::cout << proton2 << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << rang::style::bold << "Generation cuts:" << rang::style::reset << std::endl << std::endl;

    std::cout << "- System rapidity: [" << state.gcuts.Y_min << ", " << state.gcuts.Y_max << "]" << std::endl;
    std::cout << "- System mass:     [" << state.gcuts.M_min << ", " << state.gcuts.M_max << "] GeV" << std::endl;

    PrintFiducialCuts();
  }
}

// 2-dimensional phase space vector initialization
bool MCollinear::B2RandomKin(const std::vector<double> &randvec) {
  // Pick daughter masses once so rejected configurations retain zero MC weight
  PrepareDecaySymmetryProposal();
  double M_sum    = 0.0;
  double mass_max = std::min(state.lts.sqrt_s, state.gcuts.M_max);
  if (!PrepareCentralBranchMasses(M_sum, &mass_max)) { return false; }

  // Apply the central threshold, beam energy and generator boundaries
  constexpr double MASS_MARGIN = 1.0e-4;
  const double     mass_min    = std::max(M_sum + MASS_MARGIN, state.gcuts.M_min);
  if (!(mass_min < mass_max)) { return false; }

  // Sample mass squared logarithmically and rapidity uniformly
  const double mass2         = SampleCentralMassSquared(randvec[0], mass_min, mass_max);
  state.lts.central_phase_space_mass_cut_min = state.gcuts.M_min;
  state.lts.central_phase_space_mass_max = mass_max;
  state.lts.central_phase_space_mass_margin = MASS_MARGIN;
  state.lts.central_phase_space_generated_jacobian = CentralMassSquaredJacobian(mass2, mass_min, mass_max);
  const double rapidity_lab  = state.gcuts.Y_min + randvec[1] * (state.gcuts.Y_max - state.gcuts.Y_min);
  const double rapidity      = rapidity_lab - kinematics::CollisionRapidity(state.lts.pbeam1, state.lts.pbeam2);
  const double mass_fraction = msqrt(mass2 / state.lts.s);
  const double x1            = mass_fraction * std::exp(rapidity);
  const double x2            = mass_fraction * std::exp(-rapidity);

  collinear_integral_volume =
      state.lts.central_phase_space_generated_jacobian * (state.gcuts.Y_max - state.gcuts.Y_min) / state.lts.s;
  if (!(x1 > 0.0 && x1 < 1.0 && x2 > 0.0 && x2 < 1.0)) { return false; }

  return B2BuildKin(x1, x2);
}

// Build kinematics for 2 -> 1 x 1 -> N collinear
bool MCollinear::B2BuildKin(double x1, double x2) {
  const M4Vec beamsum = state.lts.pbeam1 + state.lts.pbeam2;

  // We-work in CMS-frame

  // Initial state collinear (pt=0) and massless (m=0) parton 4-momentum
  M4Vec q1(0, 0, x1 * state.lts.sqrt_s / 2.0, x1 * state.lts.sqrt_s / 2.0);
  M4Vec q2(0, 0, -x2 * state.lts.sqrt_s / 2.0, x2 * state.lts.sqrt_s / 2.0);

  // ------------------------------------------------------------------
  // Now boost if asymmetric beams
  if (std::abs(beamsum.Pz()) > 1e-6) {
    constexpr int sign = 1;  // positive -> boost to the lab
    kinematics::LorentzBoost(beamsum, state.lts.sqrt_s, q1, sign);
    kinematics::LorentzBoost(beamsum, state.lts.sqrt_s, q2, sign);
  }
  q1.SetE(q1.P3mod());
  q2.SetE(q2.P3mod());
  // ------------------------------------------------------------------

  M4Vec p1 = state.lts.pbeam1 - q1;  // Remnant
  M4Vec p2 = state.lts.pbeam2 - q2;  // Remnant
  M4Vec pX = q1 + q2;                // System

  // Save
  state.lts.pfinal[1] = p1;
  state.lts.pfinal[2] = p2;
  state.lts.pfinal[0] = pX;  // Central system

  // -------------------------------------------------------------------
  // Kinematic checks

  // Total 4-momentum conservation
  if (!CheckEMC(beamsum - (state.lts.pfinal[1] + state.lts.pfinal[2] + state.lts.pfinal[0]))) { return false; }

  // ==============================================================================
  // Central system decay tree first branch kinematics set up here, the
  // rest is done recursively

  // Mother mass
  std::vector<double> masses;

  // Collect decay product masses
  for (const auto &i : indices(state.lts.decaytree)) {
    // Use the sampled off-shell daughter masses
    masses.push_back(state.lts.decaytree[i].m_offshell);
  }
  std::vector<M4Vec> products;

  // false if amplitude has dependence on the final state legs (generic),
  // true if amplitude is a function of central system kinematics only (limited)
  const bool UNWEIGHT = !state.lts.PS_active;

  gra::kinematics::MCW w;
  // One production root inherits the complete central-system momentum
  if (state.lts.decaytree.size() == 1) {
    if (!SetSingleCentralRootKinematics()) { return false; }
  } else if (state.lts.decaytree.size() == 2) {
    w = gra::kinematics::TwoBodyPhaseSpace(state.lts.pfinal[0], state.lts.pfinal[0].M(), masses, products,
                                           state.random);
    // 3-body
  } else if (state.lts.decaytree.size() == 3) {
    w = gra::kinematics::ThreeBodyPhaseSpace(state.lts.pfinal[0], state.lts.pfinal[0].M(), masses, products, UNWEIGHT,
                                             state.random);
    // N-body
  } else if (state.lts.decaytree.size() > 3) {
    w = gra::kinematics::NBodyPhaseSpace(state.lts.pfinal[0], state.lts.pfinal[0].M(), masses, products, UNWEIGHT,
                                         state.random);
  }

  if (state.lts.decaytree.size() != 1) {
    if (w.GetW() < 0) {
      return false;  // Kinematically impossible
    }
    state.lts.DW = w;

    // Collect decay products
    const unsigned int offset = 3;
    for (const auto &i : indices(state.lts.decaytree)) {
      state.lts.decaytree[i].p4    = products[i];
      state.lts.pfinal[i + offset] = products[i];
    }
  }

  // Treat decaytree recursively
  for (const auto &i : indices(state.lts.decaytree)) {
    if (!ConstructDecayKinematics(state.lts.decaytree[i])) { return false; }
  }
  if (!ApplyDecaySymmetryProposal()) { return false; }

  // ==============================================================================
  // Check that we are above mass threshold -> not necessary, this is
  // done in mass sampling function
  const unsigned int Nf = state.lts.decaytree.size() + 2;
  if (!kinematics::SetLorentzScalars(state, Nf, false, kinematics::TransferVirtualityPolicy::CollinearLightlike)) {
    return false;
  }

  // Keep the sampled collinear PDF fractions instead of massive-beam losses
  state.lts.x1 = x1;
  state.lts.x2 = x2;
  return true;
}

// Calculate Pure Phase Space Decay Width (Volume)
void MCollinear::DecayWidthPS(double &exact) const { exact = CentralDecayWidthPS(); }

// 2-Dim Integral Volume
double MCollinear::B2IntegralVolume() const { return collinear_integral_volume; }

// 2-Dim phase space weight
double MCollinear::B2PhaseSpaceWeight() const { return 1.0; }

}  // namespace gra
