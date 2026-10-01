// Factorized type phase space class <F> for 2->3 x 1->N-2 (central system)
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
#include "Graniitti/Kinematics/MEventKinematics.h"
#include "Graniitti/Kinematics/MFactorized.h"
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
MFactorized::MFactorized() { Initialize(); }

// Constructor
MFactorized::MFactorized(std::string process, const std::vector<aux::OneCMD> &syntax, MModelTunePtr tune) {
  Initialize();
  InitHistograms();
  SetProcess(process, syntax, std::move(tune));

  // Init final states
  M4Vec zerovec(0, 0, 0, 0);
  state.lts.pfinal.assign(11, zerovec);
  std::cout << "MFactorized:: [Constructor done]" << std::endl;
}

// Construct a factorized phase space behind one compatible F or C selector
MFactorized::MFactorized(std::string process, const std::vector<aux::OneCMD> &syntax, const std::string &mode,
                         MModelTunePtr tune) {
  Initialize(mode);
  InitHistograms();
  SetProcess(process, syntax, std::move(tune));

  M4Vec zerovec(0, 0, 0, 0);
  state.lts.pfinal.assign(11, zerovec);
  std::cout << "MFactorized:: [Constructor done]" << std::endl;
}

// Initialize the factorized phase space with one compatible selector tag
void MFactorized::Initialize(const std::string &mode) {
  if (mode != "F" && mode != "C") { throw std::invalid_argument("MFactorized::Initialize: mode must be F or C"); }
  const std::vector<std::string> supported = {"MP", "XP", "GP", "TP", "ygg", "yy", "gg"};
  state.phase_space_class                  = mode;
  ProcPtr                                  = MSubProc(supported, state.phase_space_class);
  state.lts.central_phase_space_mode       = CentralPhaseSpaceMode::Factorized;
}

// Destructor
MFactorized::~MFactorized() {}

// Initialize cut and process spesific postsetup
void MFactorized::FinalizeProcessConfiguration() {
  if (state.lts.decaytree.size() > 8) {
    throw std::invalid_argument(
        "MFactorized::FinalizeProcessConfiguration: direct central "
        "multiplicity cannot exceed 8");
  }

  // Set sampling boundaries
  SetTechnicalBoundaries(state.gcuts, state.excitation);

  // Initialize phase space dimension
  ProcPtr.LIPSDIM = 5 + 1;  // All processes, +1 from central system mass

  if (state.excitation == 1) { ProcPtr.LIPSDIM += 1; }
  if (state.excitation == 2) { ProcPtr.LIPSDIM += 2; }
}

// Map one unit coordinate logarithmically onto the central mass squared
// M2 = Mmin^2 exp[u log(Mmax^2/Mmin^2)]
double MFactorized::SampleCentralMassSquared(double unit, double mass_min, double mass_max) {
  const double log_mass2_min   = 2.0 * std::log(mass_min);
  const double log_mass2_range = 2.0 * (std::log(mass_max) - std::log(mass_min));
  return std::exp(log_mass2_min + unit * log_mass2_range);
}

// Compute the local Jacobian of the logarithmic central-mass map
// dM2/du = M2 log(Mmax^2/Mmin^2)
double MFactorized::CentralMassSquaredJacobian(double mass_squared, double mass_min, double mass_max) {
  const double log_mass2_range = 2.0 * (std::log(mass_max) - std::log(mass_min));
  return mass_squared * log_mass2_range;
}

// Update kinematics (screening kT loop calls this)
// Exact 4-momentum conservation at loop vertices, that is, by using this one
// does not assume vanishing external momenta in the screening loop calculation

bool MFactorized::LoopKinematics(const std::array<double, 2> &p1p, const std::array<double, 2> &p2p) {
  if (!kinematics::RebuildScreeningKinematics(state.lts, p1p, p2p, true)) { return false; }
  return kinematics::SetLorentzScalars(state, state.lts.decaytree.size() + 2, true);
}

// Refresh factorized invariants from the restored Born four-momenta
bool MFactorized::RefreshBornKinematics() {
  return kinematics::SetLorentzScalars(state, state.lts.decaytree.size() + 2, true);
}

// Compute Monte Carlo integrand weight
double MFactorized::ComputeEventWeight(const std::vector<double> &randvec, MEventWeightState &aux) {
  double W = 0.0;

  if (parton::PrepareFinalStateProposal(state.lts, state.random)) { CalculateSymmetryFactor(); }

  // Kinematics and cuts
  PreparePhaseSpacePoint(B51RandomKin(randvec), aux);

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
    if (state.lts.decaytree.size() != 0 && state.lts.PS_active) {  // We have some legs in the central system
      C_space = (state.lts.DW.Integral() / (2 * PI));              // /(2*PI) from phase space factorization

      // --------------------------------------------------------------------
      // Cascade resonances phase-space
      C_space *= CascadePS();
      // --------------------------------------------------------------------
    }

    // ** EVENT WEIGHT **
    const double phase_weight = DissociationCrossSectionFactor() * DecaySymmetryCompensationFactor() * C_space * (1.0 / state.symmetry_factor) *
        B51PhaseSpaceWeight() * B51IntegralVolume() * parton::ProposalWeight(state.lts) * GeV2barn /
        kinematics::MollerFlux(state.lts.pbeam1, state.lts.pbeam2);
    W = phase_weight * MatESQ;
    if (aux.qmetrics) { aux.phase_weight = phase_weight; }
  }

  return W;
}

// Print factorized phase-space initialization information
void MFactorized::PrintInit(bool silent) const {
  if (!silent) {
    PrintSetup();

    // Construct prettyprint diagram
    std::string proton1 = "-----------EL--------->";
    std::string proton2 = "-----------EL--------->";

    if (state.excitation == 1) { proton1 = "-----------F2-xxxxxxxx>"; }
    if (state.excitation == 2) {
      proton1 = "-----------F2-xxxxxxxx>";
      proton2 = "-----------F2-xxxxxxxx>";
    }

    std::vector<std::string> feynmangraph;
    feynmangraph = {"||          ", "||          ", "xx--------->", "||          ", "||          "};

    // Print diagram
    const bool nuclear_survival = state.lts.upc_model != nullptr && state.lts.upc_model->HadronicConvolution();
    std::cout << proton1 << std::endl;
    for (const auto &i : indices(feynmangraph)) {
      if (nuclear_survival) {
        std::cout << rang::fg::yellow << "     **    " << rang::style::reset;
      } else if (state.lts.upc_model == nullptr && state.screening) {
        std::cout << rang::fg::red << "     **    " << rang::style::reset;
      } else {
        std::cout << "           ";
      }
      std::cout << feynmangraph[i] << std::endl;
    }
    std::cout << proton2 << std::endl;
    std::cout << std::endl;

    std::cout << std::endl;
    std::cout << rang::style::bold << "Generation cuts:" << rang::style::reset << std::endl << std::endl;
    printf(
        "- Rap : System rapidity     [min, max] = [%0.2f, %0.2f]     "
        "\t(user) \n"
        "- M   : System mass         [min, max] = [%0.2f, %0.2f] GeV "
        "\t(user) \n"
        "- Pt  : Forward leg         [min, max] = [%0.2f, %0.2f] GeV "
        "\t(fixed/user) \n",
        state.gcuts.Y_min, state.gcuts.Y_max, state.gcuts.M_min, state.gcuts.M_max, state.gcuts.forward_pt_min,
        state.gcuts.forward_pt_max);

    if (state.excitation != 0) {
      printf(
          "- Xi  : Forward leg (M^2/s) [min, max] = [%0.2E, %0.2E]     "
          "\t(fixed/user) \n",
          state.gcuts.XI_min, state.gcuts.XI_max);
    }

    PrintFiducialCuts();
  }
}

// 5+1-dimensional phase space vector initialization
bool MFactorized::B51RandomKin(const std::vector<double> &randvec) {
  // log-change of variables for pt
  const double u1 =
      std::log(state.gcuts.forward_pt_min + ZERO_EPS) +
      (std::log(state.gcuts.forward_pt_max) - std::log(state.gcuts.forward_pt_min + ZERO_EPS)) * randvec[0];
  const double u2 =
      std::log(state.gcuts.forward_pt_min + ZERO_EPS) +
      (std::log(state.gcuts.forward_pt_max) - std::log(state.gcuts.forward_pt_min + ZERO_EPS)) * randvec[1];

  const double pt1 = std::exp(u1);
  const double pt2 = std::exp(u2);

  const double phi1  = 2.0 * gra::math::PI * randvec[2];
  const double phi2  = 2.0 * gra::math::PI * randvec[3];
  const double y_lab = state.gcuts.Y_min + (state.gcuts.Y_max - state.gcuts.Y_min) * randvec[4];
  const double yX    = y_lab - kinematics::CollisionRapidity(state.lts.pbeam1, state.lts.pbeam2);

  // Pick daughter masses once so rejected configurations retain zero MC weight
  PrepareDecaySymmetryProposal();
  double M_sum = 0.0;
  M_MAX        = std::min(state.lts.sqrt_s - (state.lts.pbeam1.M() + state.lts.pbeam2.M()), state.gcuts.M_max);
  if (!PrepareCentralBranchMasses(M_sum, &M_MAX)) { return false; }

  // Apply absolute boundary conditions and generator cuts
  const double MARGIN = 1e-4;  // GeV
  M_MIN               = std::max(M_sum + MARGIN, state.gcuts.M_min);
  if (!(M_MIN < M_MAX)) { return false; }

  // Sample mass squared logarithmically to resolve broad low-mass continua
  const double m2X                                      = SampleCentralMassSquared(randvec[5], M_MIN, M_MAX);
  state.lts.central_phase_space_mass_cut_min            = state.gcuts.M_min;
  state.lts.central_phase_space_mass_max                = M_MAX;
  state.lts.central_phase_space_mass_margin             = MARGIN;
  state.lts.central_phase_space_generated_jacobian = CentralMassSquaredJacobian(m2X, M_MIN, M_MAX);

  // Forward N* system masses
  std::vector<double> mvec;
  std::vector<double> rvec;
  if (state.excitation == 1) { rvec = {randvec[6]}; }
  if (state.excitation == 2) { rvec = {randvec[6], randvec[7]}; }
  SampleForwardMasses(mvec, rvec);

  return B51BuildKin(pt1, pt2, phi1, phi2, yX, m2X, mvec[0], mvec[1]);
}

// Build kinematics for 2->3 skeleton
bool MFactorized::B51BuildKin(double pt1, double pt2, double phi1, double phi2, double yX, double m2X, double m1,
                              double m2) {
  const M4Vec beamsum     = state.lts.pbeam1 + state.lts.pbeam2;
  state.lts.forward_mass2 = {pow2(m1), pow2(m2)};

  // Final state 4-momenta, set px,py first
  M4Vec p1(pt1 * std::cos(phi1), pt1 * std::sin(phi1), 0, 0);
  M4Vec p2(pt2 * std::cos(phi2), pt2 * std::sin(phi2), 0, 0);
  M4Vec pX(-(p1.Px() + p2.Px()), -(p1.Py() + p2.Py()), 0, 0);

  // Central system pz and E
  const double mtX = msqrt(m2X + pX.Pt2());
  pX.SetPzE(mtX * std::sinh(yX), mtX * std::cosh(yX));

  // Energy overflow
  if (pX.E() > (state.lts.sqrt_s - (m1 + m2))) { return false; }

  double p1z = gra::kinematics::SolvePz(m1, m2, pt1, pt2, pX.Pz(), pX.E(), state.lts.s);
  double p2z = -(pX.Pz() + p1z);  // by momentum conservation

  // Enforce scattering direction +p -> +p, -p -> -p (VERY RARE POLYNOMIAL
  // BRANCH FLIP)
  if (p1z < 0 || p2z > 0) { return false; }

  // Pz and E of forward protons/N*
  p1.SetPzE(p1z, msqrt(pow2(m1) + pow2(pt1) + pow2(p1z)));
  p2.SetPzE(p2z, msqrt(pow2(m2) + pow2(pt2) + pow2(p2z)));

  // ------------------------------------------------------------------
  // Now boost if asymmetric beams
  if (std::abs(beamsum.Pz()) > 1e-6) {
    constexpr int sign = 1;  // positive -> boost to the lab
    kinematics::LorentzBoost(beamsum, state.lts.sqrt_s, p1, sign);
    kinematics::LorentzBoost(beamsum, state.lts.sqrt_s, p2, sign);
    // Preserve the generated invariant mass when changing the longitudinal frame
    const double y_lab = yX + kinematics::CollisionRapidity(state.lts.pbeam1, state.lts.pbeam2);
    pX.SetPzE(mtX * std::sinh(y_lab), mtX * std::cosh(y_lab));
  }
  // ------------------------------------------------------------------

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

  // Collect decay product masses
  std::vector<double> masses;
  for (const auto &i : indices(state.lts.decaytree)) {
    // Use the off-shell daughter masses sampled by B51RandomKin
    masses.push_back(state.lts.decaytree[i].m_offshell);
  }
  std::vector<M4Vec> products(state.lts.decaytree.size());

  // false if amplitude has dependence on the final state legs (generic),
  // true if amplitude is a function of central system kinematics only (limited)
  const bool UNWEIGHT = !state.lts.PS_active;

  gra::kinematics::MCW w;

  // One production root inherits the complete central-system momentum
  if (state.lts.decaytree.size() == 1) {
    if (!SetSingleCentralRootKinematics()) { return false; }
  } else if (state.lts.decaytree.size() == 2) {
    w = gra::kinematics::TwoBodyPhaseSpace(state.lts.pfinal[0], msqrt(m2X), masses, products, state.random);

    // 3-body
  } else if (state.lts.decaytree.size() == 3) {
    w = gra::kinematics::ThreeBodyPhaseSpace(state.lts.pfinal[0], msqrt(m2X), masses, products, UNWEIGHT, state.random);

    // N-body
  } else if (state.lts.decaytree.size() > 3) {
    if (UNWEIGHT) {
      w = gra::kinematics::NBodyPhaseSpace(state.lts.pfinal[0], msqrt(m2X), masses, products, UNWEIGHT, state.random);
    } else {
      w = gra::kinematics::RamboMassive(state.lts.pfinal[0], msqrt(m2X), masses, products, state.random);
    }
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
  return kinematics::SetLorentzScalars(state, Nf);
}

// Calculate Pure Phase Space Decay Width (Volume)
void MFactorized::DecayWidthPS(double &exact) const { exact = CentralDecayWidthPS(); }

// 5-Dim Integral Volume \int_{M^2_MIN}^{M^2_MAX} { dM^2 } [phi1] x [phi2] x
// [pt1] x [pt2] x [y]
// Integral over central mass^2 is separate [phase-space factorization], but
// encapsulated here (+1 dimension)
double MFactorized::B51IntegralVolume() const {
  // Forward leg integration
  const double forward_volume = ForwardVolume();
  const double mass2_volume   = CentralMassSquaredJacobian(state.lts.m2, M_MIN, M_MAX);

  return mass2_volume * (state.gcuts.Y_max - state.gcuts.Y_min) * forward_volume;
}

// Compute the five-dimensional phase-space weight
// weight = dPhi_3/[dphi1 dphi2 dpT1 dpT2 dy] after the pz delta function
double MFactorized::B51PhaseSpaceWeight() const {
  const double J = 1.0 / std::abs(state.lts.pfinal[1].Pz() / state.lts.pfinal[1].E() -
                                  state.lts.pfinal[2].Pz() / state.lts.pfinal[2].E());  // Jacobian, close to 0.5

  const double factor = (1.0 / 2.0) * (1.0 / pow5(2.0 * gra::math::PI)) *
                        (state.lts.pfinal[1].Pt() / (2.0 * state.lts.pfinal[1].E())) *
                        (state.lts.pfinal[2].Pt() / (2.0 * state.lts.pfinal[2].E())) * J;

  return factor;
}

// For high mass limit kinematics, see e.g. [arxiv.org/pdf/hep-ph/9903279.pdf]

}  // namespace gra
