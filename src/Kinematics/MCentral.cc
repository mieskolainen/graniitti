// Central multiparticle type phase space class <C>
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <random>
#include <vector>

// OWN
#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/MUserCuts.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Particle/MForm.h"
#include "Graniitti/Kinematics/MCentral.h"
#include "Graniitti/Kinematics/MEventKinematics.h"
#include "Graniitti/Process/MProcess.h"
#include "Graniitti/QCD/MPartonProposal.h"
#include "Graniitti/Regge/MFragment.h"
#include "Graniitti/Spin/MHelicity.h"
#include "Graniitti/Tech/MAux.h"

// Libraries
#include "rang.hpp"

using gra::aux::indices;
using gra::math::abs2;
using gra::math::msqrt;
using gra::math::PI;
using gra::math::pow2;
using gra::math::pow3;
using gra::math::pow4;
using gra::math::zi;
using gra::PDG::GeV2barn;

namespace gra {
// This is needed by construction
MCentral::MCentral() { Initialize(); }

// Constructor
MCentral::MCentral(std::string process,
                   const std::vector<aux::OneCMD> &syntax, MModelTunePtr tune) {
  Initialize();
  InitHistograms();
  SetProcess(process, syntax, std::move(tune));

  // Init final states
  M4Vec zerovec(0, 0, 0, 0);
  state.lts.pfinal.assign(11, zerovec);
  std::cout << "MCentral:: [Constructor done]" << std::endl;
}

// Initialize central multiparticle phase space support
void MCentral::Initialize() {
  const std::vector<std::string> supported = {"MP",  "XP", "GP", "TP",
                                              "ygg", "yy", "gg"};
  state.phase_space_class = "C";
  ProcPtr = MSubProc(supported, state.phase_space_class);
  state.lts.central_phase_space_mode = CentralPhaseSpaceMode::Central;
}

// Destructor
MCentral::~MCentral() {}

// Initialize cut and process spesific postsetup
void MCentral::FinalizeProcessConfiguration() {
  std::cout << "MCentral::FinalizeProcessConfiguration: " << state.excitation
            << std::endl;

  if (state.lts.decaytree.size() < 2 || state.lts.decaytree.size() > 8) {
    throw std::invalid_argument(
        "MCentral::FinalizeProcessConfiguration: direct central "
        "multiplicity must be between 2 and 8");
  }

  if (!(state.gcuts.M_min >= 0.0) || !(state.gcuts.M_min < state.gcuts.M_max)) {
    throw std::invalid_argument(
        "MCentral::FinalizeProcessConfiguration: central mass "
        "range must satisfy 0 <= min < max");
  }

  if (!std::isfinite(state.gcuts.rap_min) ||
      !std::isfinite(state.gcuts.rap_max) ||
      !(state.gcuts.rap_min < state.gcuts.rap_max)) {
    throw std::invalid_argument(
        "MCentral::FinalizeProcessConfiguration: central rapidity "
        "range must satisfy finite min < max");
  }
  const double y_collision =
      kinematics::CollisionRapidity(state.lts.pbeam1, state.lts.pbeam2);
  state.lts.central_phase_space_rap_min_cm = state.gcuts.rap_min - y_collision;
  state.lts.central_phase_space_rap_max_cm = state.gcuts.rap_max - y_collision;

  // Set sampling boundaries
  SetTechnicalBoundaries(state.gcuts, state.excitation);

  // Initialize phase space dimension (3*Nf - 4)
  ProcPtr.LIPSDIM = 3 * (state.lts.decaytree.size() + 2) - 4;

  if (state.excitation == 1) {
    ProcPtr.LIPSDIM += 1;
  }
  if (state.excitation == 2) {
    ProcPtr.LIPSDIM += 2;
  }
}

// Compute the central mass ceiling allowed by one sampled forward point
double MCentral::CentralMassKinematicMaximum(double sqrt_s,
                                             double forward_mass1,
                                             double forward_mass2,
                                             const M4Vec &forward1,
                                             const M4Vec &forward2) {
  if (!std::isfinite(sqrt_s) || !std::isfinite(forward_mass1) ||
      !std::isfinite(forward_mass2) || !(sqrt_s > 0.0) ||
      !(forward_mass1 >= 0.0) || !(forward_mass2 >= 0.0)) {
    throw PhaseSpaceFailure("MCentral::CentralMassKinematicMaximum: "
                                "invalid mass or collision energy");
  }

  const double mt1 = msqrt(pow2(forward_mass1) + forward1.Pt2());
  const double mt2 = msqrt(pow2(forward_mass2) + forward2.Pt2());
  const double available_energy = sqrt_s - mt1 - mt2;
  if (!(available_energy > 0.0)) {
    return 0.0;
  }
  const double recoil_pt2 = (forward1 + forward2).Pt2();
  const double mass2_max = pow2(available_energy) - recoil_pt2;
  return mass2_max > 0.0 ? msqrt(mass2_max) : 0.0;
}

// Compute the squared Helmert ball radius for one mass and recoil ceiling
double
MCentral::MassConditionedTransverseRadius2(double mass_max, double recoil_pt2,
                                           unsigned int central_multiplicity) {
  if (!std::isfinite(mass_max) || !std::isfinite(recoil_pt2) ||
      !(mass_max > 0.0) || !(recoil_pt2 >= 0.0) || central_multiplicity < 2) {
    throw PhaseSpaceFailure(
        "MCentral::MassConditionedTransverseRadius2: invalid mapping domain");
  }

  // Two antiparallel legs maximize the internal norm at fixed Mmax and P_T
  // For rho^2 = sum_i |p_iT - P_T/K|^2 this gives the exact ball envelope
  return 0.5 * pow2(mass_max) +
         (1.0 - 1.0 / static_cast<double>(central_multiplicity)) * recoil_pt2;
}

// Map a logarithmic Helmert hyperradius onto transverse differences
std::vector<M4Vec>
MCentral::MapTransverseLog(const std::vector<double> &radial_units,
                           const std::vector<double> &angle_units,
                           const M4Vec &central_transverse_momentum,
                           const double radius2, const double scale2,
                           double &jacobian) {
  if (radial_units.empty() || radial_units.size() != angle_units.size() ||
      !std::isfinite(radius2) || !std::isfinite(scale2) || !(radius2 > 0.0) ||
      !(scale2 > 0.0)) {
    throw PhaseSpaceFailure(
        "MCentral::MapTransverseLog: invalid mapping domain");
  }
  for (const auto &i : indices(radial_units)) {
    if (!std::isfinite(radial_units[i]) || !std::isfinite(angle_units[i]) ||
        !(radial_units[i] >= 0.0 && radial_units[i] <= 1.0) ||
        !(angle_units[i] >= 0.0 && angle_units[i] <= 1.0)) {
      throw PhaseSpaceFailure(
          "MCentral::MapTransverseLog: unit coordinate outside [0, 1]");
    }
  }

  const unsigned int mode_count = radial_units.size();
  const unsigned int multiplicity = mode_count + 1;

  // Sample t = rho^2 with density proportional to 1 / (scale2 + t)
  // This is logarithmic above scale2 and finite at the physical origin
  const double log_range = std::log1p(radius2 / scale2);
  const double rho2 = scale2 * std::expm1(radial_units[0] * log_range);

  // Split rho^2 symmetrically as |h_a|^2 = rho^2 z_a with sum_a z_a = 1
  // The inverse stick map generates Dirichlet(1,...,1) mode fractions
  std::vector<double> fraction(mode_count, 0.0);
  double remaining = 1.0;
  for (unsigned int i = 0; i + 1 < mode_count; ++i) {
    const double remaining_modes = static_cast<double>(mode_count - i - 1);
    const double stick =
        1.0 - std::pow(1.0 - radial_units[i + 1], 1.0 / remaining_modes);
    fraction[i] = remaining * stick;
    remaining *= 1.0 - stick;
  }
  fraction.back() = remaining;

  // Give each orthonormal two-dimensional Helmert mode a uniform azimuth
  std::vector<M4Vec> modes(mode_count);
  for (const auto &i : indices(modes)) {
    const double magnitude = msqrt(rho2 * fraction[i]);
    const double angle = 2.0 * PI * angle_units[i];
    modes[i] = M4Vec(magnitude * std::cos(angle), magnitude * std::sin(angle),
                     0.0, 0.0);
  }

  // Invert the Helmert basis around P_T / K, preserving exact momentum closure
  std::vector<M4Vec> momenta(multiplicity,
                             central_transverse_momentum /
                                 static_cast<double>(multiplicity));
  for (const auto &j : indices(modes)) {
    const double norm = msqrt(static_cast<double>((j + 1) * (j + 2)));
    for (unsigned int i = 0; i <= j; ++i) {
      momenta[i] += modes[j] / norm;
    }
    momenta[j + 1] -= modes[j] * (static_cast<double>(j + 1) / norm);
  }

  std::vector<M4Vec> differences(mode_count);
  for (const auto &i : indices(differences)) {
    differences[i] = momenta[i] - momenta[i + 1];
  }

  // The unit-cube Jacobian is
  // K pi^n t^(n-1) [dt/du] / (n-1)!, where n = K-1
  // The factor K is the two-dimensional Helmert-to-difference determinant
  // Its integral is K (pi R^2)^n / n!, the exact difference-ball volume
  // For K = 2, q^2 = 2t and scale2 = q0^2 / 2, giving the two-body log map
  jacobian = kinematics::TransverseLogJacobian(rho2, radius2, scale2, multiplicity);
  if (!std::isfinite(jacobian) || jacobian < 0.0) {
    throw PhaseSpaceFailure(
        "MCentral::MapTransverseLog: non-finite Jacobian");
  }
  return differences;
}

// Update kinematics (screening kT loop calls this)
// Exact 4-momentum conservation at loop vertices, that is, by using this one
// does not assume vanishing external momenta in the screening loop calculation
bool MCentral::LoopKinematics(const std::array<double, 2> &p1p,
                              const std::array<double, 2> &p2p) {
  if (!kinematics::RebuildScreeningKinematics(state.lts, p1p, p2p, true)) {
    return false;
  }
  return kinematics::SetLorentzScalars(state, state.lts.decaytree.size() + 2, true);
}

// Refresh central phase space invariants from the restored Born four-momenta
bool MCentral::RefreshBornKinematics() {
  return kinematics::SetLorentzScalars(state, state.lts.decaytree.size() + 2, true);
}

// Get weight
double MCentral::ComputeEventWeight(const std::vector<double> &randvec,
                                    MEventWeightState &aux) {
  double W = 0.0;

  if (parton::PrepareFinalStateProposal(state.lts, state.random)) {
    CalculateSymmetryFactor();
  }

  // Kinematics and sequential cuts
  PreparePhaseSpacePoint(BNRandomKin(state.lts.decaytree.size() + 2, randvec),
                         aux);

  if (aux.Valid()) {
    // Matrix element squared
    const double MatESQ = GetAmp2(aux.include_screening, aux);

    // --------------------------------------------------------------------
    // Cascade resonances phase-space
    const double C_space = CascadePS();
    // --------------------------------------------------------------------

    // ** EVENT WEIGHT **
    const double phase_weight = DissociationCrossSectionFactor() * DecaySymmetryCompensationFactor() *
        C_space * (1.0 / state.symmetry_factor) * BNPhaseSpaceWeight() *
        BNIntegralVolume() * parton::ProposalWeight(state.lts) *
        GeV2barn / kinematics::MollerFlux(state.lts.pbeam1, state.lts.pbeam2);
    W = phase_weight * MatESQ;
    if (aux.qmetrics) { aux.phase_weight = phase_weight; }
  }

  return W;
}

// Print central phase space initialization information
void MCentral::PrintInit(bool silent) const {
  if (!silent) {
    PrintSetup();

    // Construct prettyprint diagram
    std::string proton1 = "-----------EL--------->";
    std::string proton2 = "-----------EL--------->";

    if (state.excitation == 1) {
      proton1 = "-----------F2-xxxxxxxx>";
    }
    if (state.excitation == 2) {
      proton1 = "-----------F2-xxxxxxxx>";
      proton2 = "-----------F2-xxxxxxxx>";
    }

    std::vector<std::string> legs;
    for (const auto &i : indices(state.lts.decaytree)) {
      char buff[250];
      snprintf(
          buff, sizeof(buff), "x---------> %d (%s) [Q=%s, J=%s]",
          state.lts.decaytree[i].p.pdg, state.lts.decaytree[i].p.name.c_str(),
          gra::aux::Charge3XtoString(state.lts.decaytree[i].p.chargeX3).c_str(),
          gra::aux::NullableSpin2XtoString(state.lts.decaytree[i].p.spinX2)
              .c_str());
      std::string leg = buff;
      legs.push_back(leg);
    }

    std::vector<std::string> feynmangraph;
    feynmangraph.push_back("||         ");
    for (const auto &i : indices(state.lts.decaytree)) {
      feynmangraph.push_back(legs[i]);
      if (i < state.lts.decaytree.size() - 1)
        feynmangraph.push_back("|          ");
    }
    feynmangraph.push_back("||         ");

    // Print diagram
    const bool nuclear_survival =
        state.lts.upc_model != nullptr && state.lts.upc_model->HadronicConvolution();
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

    // Generation cuts
    std::cout << std::endl;
    std::cout << rang::style::bold << "Generation cuts:" << rang::style::reset
              << std::endl
              << std::endl;
    const double recoil_pt_max = 2.0 * state.gcuts.forward_pt_max;
    const double radius2_max = MassConditionedTransverseRadius2(
        state.gcuts.M_max, pow2(recoil_pt_max),
        static_cast<unsigned int>(state.lts.decaytree.size()));
    const double kt_max = msqrt(2.0 * radius2_max);

    printf("- Rap : Final state rapidity [min, max] = [%0.2f, %0.2f]     "
           "\t(user)       \n"
           "- M   : Central system       [min, max] = [%0.2f, %0.2f] GeV "
           "\t(user)       \n"
           "- Kt  : Adjacent difference  [min, max] = [0.00, %0.2f] GeV "
           "\t(derived bounds)   \n"
           "- Pt  : Forward leg          [min, max] = [%0.2f, %0.2f] GeV "
           "\t(fixed/user) "
           "\n",
           state.gcuts.rap_min, state.gcuts.rap_max, state.gcuts.M_min,
           state.gcuts.M_max, kt_max, state.gcuts.forward_pt_min,
           state.gcuts.forward_pt_max);

    if (state.excitation != 0) {
      printf("- Xi  : Forward leg (M^2/s)  [min, max] = [%0.2E, %0.2E]     "
             "\t(fixed/user) \n",
             state.gcuts.XI_min, state.gcuts.XI_max);
    }

    PrintFiducialCuts();
  }
}

// (3*Nf-4)-dimensional phase space vector initialization
bool MCentral::BNRandomKin(unsigned int Nf,
                           const std::vector<double> &randvec) {
  const unsigned int Kf = Nf - 2; // Central system multiplicity
  state.lts.central_phase_space_generated_jacobian = 0.0;

  // log-change of variables for pt
  const double u1 = std::log(state.gcuts.forward_pt_min + ZERO_EPS) +
                    (std::log(state.gcuts.forward_pt_max) -
                     std::log(state.gcuts.forward_pt_min + ZERO_EPS)) *
                        randvec[0];
  const double u2 = std::log(state.gcuts.forward_pt_min + ZERO_EPS) +
                    (std::log(state.gcuts.forward_pt_max) -
                     std::log(state.gcuts.forward_pt_min + ZERO_EPS)) *
                        randvec[1];

  const double pt1 = std::exp(u1);
  const double pt2 = std::exp(u2);

  const double phi1 = 2.0 * gra::math::PI * randvec[2];
  const double phi2 = 2.0 * gra::math::PI * randvec[3];
  const unsigned int offset = 4; // 4 variables above

  // Decay product masses
  // ==============================================================
  PrepareDecaySymmetryProposal();
  for (const auto &i : indices(state.lts.decaytree)) {
    GetOffShellMass(state.lts.decaytree[i], state.lts.decaytree[i].m_offshell);
    if (state.lts.decaytree[i].m_offshell < 0.0) {
      return false;
    }
  }
  // ==============================================================

  // Sample forward system masses before constructing the transverse envelope
  const std::size_t forward_mass_index = offset + 2 * (Kf - 1) + Kf;
  std::vector<double> mvec;
  std::vector<double> rvec;
  if (state.excitation == 1) {
    rvec = {randvec[forward_mass_index]};
  }
  if (state.excitation == 2) {
    rvec = {randvec[forward_mass_index], randvec[forward_mass_index + 1]};
  }
  SampleForwardMasses(mvec, rvec);

  // Construct the mass-conditioned transverse difference envelope
  const M4Vec forward1(pt1 * std::cos(phi1), pt1 * std::sin(phi1), 0.0, 0.0);
  const M4Vec forward2(pt2 * std::cos(phi2), pt2 * std::sin(phi2), 0.0, 0.0);
  const double mass_kinematic_max = CentralMassKinematicMaximum(
      state.lts.sqrt_s, mvec[0], mvec[1], forward1, forward2);
  const double mass_max = std::min(state.gcuts.M_max, mass_kinematic_max);
  double offshell_mass_sum = 0.0;
  for (const auto &branch : state.lts.decaytree) {
    offshell_mass_sum += branch.m_offshell;
  }
  const double mass_threshold = std::max(state.gcuts.M_min, offshell_mass_sum);
  if (!(mass_threshold < mass_max)) {
    return false;
  }

  const M4Vec central_transverse_momentum = -(forward1 + forward2);
  const double recoil_pt2 = central_transverse_momentum.Pt2();
  const double radius2 =
      MassConditionedTransverseRadius2(mass_max, recoil_pt2, Kf);

  // Map the transverse unit coordinates through one correlated Helmert ball
  std::vector<double> radial_units(Kf - 1);
  std::vector<double> angle_units(Kf - 1);
  size_t ind = offset;
  for (auto &unit : radial_units) {
    unit = randvec[ind++];
  }
  for (auto &unit : angle_units) {
    unit = randvec[ind++];
  }
  // Regularize the logarithm at the physical mass scale or one percent of Mmax
  state.lts.central_phase_space_mass_max = mass_max;
  state.lts.central_phase_space_transverse_radius2 = radius2;
  const double scale = std::max(offshell_mass_sum, 0.01 * mass_max);
  const std::vector<M4Vec> differences = MapTransverseLog(
      radial_units, angle_units, central_transverse_momentum, radius2,
      pow2(scale) / static_cast<double>(Kf), state.lts.central_phase_space_generated_jacobian);

  // Convert the difference vectors for the existing central momentum builder
  std::vector<double> kt(Kf - 1, 0.0);  // Kf-1
  std::vector<double> phi(Kf - 1, 0.0); // Kf-1
  for (const auto &i : indices(differences)) {
    kt[i] = differences[i].Pt();
    phi[i] = std::atan2(differences[i].Py(), differences[i].Px());
  }

  // Map laboratory rapidity steering into the collision CM
  const double y_collision =
      kinematics::CollisionRapidity(state.lts.pbeam1, state.lts.pbeam2);
  std::vector<double> y(Kf, 0.0); // Kf
  for (const auto &i : indices(y)) {
    y[i] = state.gcuts.rap_min +
           (state.gcuts.rap_max - state.gcuts.rap_min) * randvec[ind] -
           y_collision;
    ++ind;
  }

  return BNBuildKin(Nf, pt1, pt2, phi1, phi2, kt, phi, y, mvec[0], mvec[1]);
}

// Build kinematics of 2->N
bool MCentral::BNBuildKin(unsigned int Nf, double pt1, double pt2, double phi1,
                          double phi2, const std::vector<double> &kt,
                          const std::vector<double> &phi,
                          const std::vector<double> &y, double m1, double m2) {
  const unsigned int Kf = Nf - 2; // Central system multiplicity
  const M4Vec beamsum = state.lts.pbeam1 + state.lts.pbeam2;
  state.lts.forward_mass2 = {pow2(m1), pow2(m2)};

  // Forward protons px,py
  M4Vec p1(pt1 * std::cos(phi1), pt1 * std::sin(phi1), 0, 0);
  M4Vec p2(pt2 * std::cos(phi2), pt2 * std::sin(phi2), 0, 0);

  // Auxialary "difference momentum" q0 = p0 - p1 ...
  pkt_.resize(Kf - 1);
  for (const auto &i : indices(pkt_)) {
    pkt_[i] = M4Vec(kt[i] * std::cos(phi[i]), kt[i] * std::sin(phi[i]), 0, 0);
  }

  // Apply linear system to get p
  std::vector<M4Vec> p;
  if (!gra::kinematics::BuildCentralTransverseMomenta(p1, p2, pkt_, p)) { return false; }

  // Set pz and E for central final states
  M4Vec sumP(0, 0, 0, 0);
  for (const auto &i : indices(p)) {
    const double m = state.lts.decaytree[i].m_offshell; // Note offshell!
    const double mt = msqrt(pow2(m) + p[i].Pt2());
    p[i].SetPzE(mt * std::sinh(y[i]), mt * std::cosh(y[i]));
    sumP += p[i];
  }

  // Check crude energy overflow
  if (sumP.E() > state.lts.sqrt_s) {
    return false;
  }

  // Apply the central system generation mass window
  const double central_mass = sumP.M();
  if (!std::isfinite(central_mass) || central_mass < state.gcuts.M_min ||
      central_mass > state.gcuts.M_max) {
    return false;
  }

  double p1z = gra::kinematics::SolvePz(m1, m2, pt1, pt2, sumP.Pz(), sumP.E(),
                                        state.lts.s);
  double p2z = -(sumP.Pz() + p1z); // by momentum conservation

  // Enforce scattering direction +p -> +p, -p -> -p (VERY RARE POLYNOMIAL
  // BRANCH FLIP)
  if (p1z < 0 || p2z > 0) {
    return false;
  }

  // pz and E of protons/N*
  p1.SetPzE(p1z, msqrt(pow2(m1) + pow2(pt1) + pow2(p1z)));
  p2.SetPzE(p2z, msqrt(pow2(m2) + pow2(pt2) + pow2(p2z)));

  // ------------------------------------------------------------------
  // Now boost if asymmetric beams
  if (std::abs(beamsum.Pz()) > 1e-6) {
    constexpr int sign = 1; // positive -> boost to the lab
    kinematics::LorentzBoost(beamsum, state.lts.sqrt_s, p1, sign);
    kinematics::LorentzBoost(beamsum, state.lts.sqrt_s, p2, sign);
    // Preserve each generated mass and sum the daughters in the laboratory frame
    const double y_collision = kinematics::CollisionRapidity(state.lts.pbeam1, state.lts.pbeam2);
    sumP = M4Vec();
    for (const auto &i : indices(p)) {
      const double mt = std::hypot(state.lts.decaytree[i].m_offshell, p[i].Pt());
      const double y_lab = y[i] + y_collision;
      p[i].SetPzE(mt * std::sinh(y_lab), mt * std::cosh(y_lab));
      sumP += p[i];
    }
  }
  // ------------------------------------------------------------------

  // First branch kinematics
  state.lts.pfinal[1] = p1; // Forward systems
  state.lts.pfinal[2] = p2;
  state.lts.pfinal[0] = sumP; // Central system

  double sumM = 0;
  const unsigned int offset = 3;
  for (const auto &i : indices(p)) {
    state.lts.decaytree[i].p4 = p[i];
    state.lts.pfinal[i + offset] = p[i];
    sumM += p[i].M();
  }

  // -------------------------------------------------------------------
  // Kinematic checks

  // Check we are above mass threshold
  if (sumP.M() < sumM) {
    return false;
  }

  // Total 4-momentum conservation
  if (!gra::math::CheckEMC(
          beamsum -
          (state.lts.pfinal[1] + state.lts.pfinal[2] + state.lts.pfinal[0]))) {
    return false;
  }

  // -------------------------------------------------------------------

  // Treat decaytree recursively
  for (const auto &i : indices(state.lts.decaytree)) {
    if (!ConstructDecayKinematics(state.lts.decaytree[i])) {
      return false;
    }
  }
  if (!ApplyDecaySymmetryProposal()) {
    return false;
  }

  return kinematics::SetLorentzScalars(state, Nf);
}

/*
The inverse-matrix representation of the transverse linear system

For K central particles, w = p1f + p2f and the auxiliary vector was

  b[0] =  q[0] - w
  b[i] = -q[i - 1] - w,  i = 1,...,K-1

The solution p = A[K - 2] b below is algebraically equivalent to the
adjacent-difference recurrence used by BuildCentralTransverseMomenta.

static const std::vector<std::vector<std::vector<double>>> A = {
    {{1.0 / 2.0, 0.0}, {0.0, 1.0 / 2.0}},

    {{2.0 / 3.0, 0.0, -1.0 / 3.0},
     {1.0 / 6.0, 1.0 / 2.0, -1.0 / 3.0},
     {-1.0 / 3.0, 0.0, 2.0 / 3.0}},

    {{7.0 / 8.0, 1.0 / 8.0, -1.0 / 2.0, -1.0 / 4.0},
     {3.0 / 8.0, 5.0 / 8.0, -1.0 / 2.0, -1.0 / 4.0},
     {-1.0 / 8.0, 1.0 / 8.0, 1.0 / 2.0, -1.0 / 4.0},
     {-5.0 / 8.0, -3.0 / 8.0, 1.0 / 2.0, 3.0 / 4.0}},

    {{11.0 / 10.0, 3.0 / 10.0, -3.0 / 5.0, -2.0 / 5.0, -1.0 / 5.0},
     {3.0 / 5.0, 4.0 / 5.0, -3.0 / 5.0, -2.0 / 5.0, -1.0 / 5.0},
     {1.0 / 10.0, 3.0 / 10.0, 2.0 / 5.0, -2.0 / 5.0, -1.0 / 5.0},
     {-2.0 / 5.0, -1.0 / 5.0, 2.0 / 5.0, 3.0 / 5.0, -1.0 / 5.0},
     {-9.0 / 10.0, -7.0 / 10.0, 2.0 / 5.0, 3.0 / 5.0, 4.0 / 5.0}},

    {{4.0 / 3.0, 1.0 / 2.0, -2.0 / 3.0, -1.0 / 2.0, -1.0 / 3.0, -1.0 / 6.0},
     {5.0 / 6.0, 1.0 / 1.0, -2.0 / 3.0, -1.0 / 2.0, -1.0 / 3.0, -1.0 / 6.0},
     {1.0 / 3.0, 1.0 / 2.0, 1.0 / 3.0, -1.0 / 2.0, -1.0 / 3.0, -1.0 / 6.0},
     {-1.0 / 6.0, 0.0, 1.0 / 3.0, 1.0 / 2.0, -1.0 / 3.0, -1.0 / 6.0},
     {-2.0 / 3.0, -1.0 / 2.0, 1.0 / 3.0, 1.0 / 2.0, 2.0 / 3.0, -1.0 / 6.0},
     {-7.0 / 6.0, -1.0 / 1.0, 1.0 / 3.0, 1.0 / 2.0, 2.0 / 3.0, 5.0 / 6.0}},

    {{11.0 / 7.0, 5.0 / 7.0, -5.0 / 7.0, -4.0 / 7.0, -3.0 / 7.0, -2.0 / 7.0,
      -1.0 / 7.0},
     {15.0 / 14.0, 17.0 / 14.0, -5.0 / 7.0, -4.0 / 7.0, -3.0 / 7.0, -2.0 / 7.0,
      -1.0 / 7.0},
     {4.0 / 7.0, 5.0 / 7.0, 2.0 / 7.0, -4.0 / 7.0, -3.0 / 7.0, -2.0 / 7.0,
      -1.0 / 7.0},
     {1.0 / 14.0, 3.0 / 14.0, 2.0 / 7.0, 3.0 / 7.0, -3.0 / 7.0, -2.0 / 7.0,
      -1.0 / 7.0},
     {-3.0 / 7.0, -2.0 / 7.0, 2.0 / 7.0, 3.0 / 7.0, 4.0 / 7.0, -2.0 / 7.0,
      -1.0 / 7.0},
     {-13.0 / 14.0, -11.0 / 14.0, 2.0 / 7.0, 3.0 / 7.0, 4.0 / 7.0, 5.0 / 7.0,
      -1.0 / 7.0},
     {-10.0 / 7.0, -9.0 / 7.0, 2.0 / 7.0, 3.0 / 7.0, 4.0 / 7.0, 5.0 / 7.0,
      6.0 / 7.0}},

    {{29.0 / 16.0, 15.0 / 16.0, -3.0 / 4.0, -5.0 / 8.0, -1.0 / 2.0,
      -3.0 / 8.0, -1.0 / 4.0, -1.0 / 8.0},
     {21.0 / 16.0, 23.0 / 16.0, -3.0 / 4.0, -5.0 / 8.0, -1.0 / 2.0,
      -3.0 / 8.0, -1.0 / 4.0, -1.0 / 8.0},
     {13.0 / 16.0, 15.0 / 16.0, 1.0 / 4.0, -5.0 / 8.0, -1.0 / 2.0,
      -3.0 / 8.0, -1.0 / 4.0, -1.0 / 8.0},
     {5.0 / 16.0, 7.0 / 16.0, 1.0 / 4.0, 3.0 / 8.0, -1.0 / 2.0,
      -3.0 / 8.0, -1.0 / 4.0, -1.0 / 8.0},
     {-3.0 / 16.0, -1.0 / 16.0, 1.0 / 4.0, 3.0 / 8.0, 1.0 / 2.0,
      -3.0 / 8.0, -1.0 / 4.0, -1.0 / 8.0},
     {-11.0 / 16.0, -9.0 / 16.0, 1.0 / 4.0, 3.0 / 8.0, 1.0 / 2.0,
      5.0 / 8.0, -1.0 / 4.0, -1.0 / 8.0},
     {-19.0 / 16.0, -17.0 / 16.0, 1.0 / 4.0, 3.0 / 8.0, 1.0 / 2.0,
      5.0 / 8.0, 3.0 / 4.0, -1.0 / 8.0},
     {-27.0 / 16.0, -25.0 / 16.0, 1.0 / 4.0, 3.0 / 8.0, 1.0 / 2.0,
      5.0 / 8.0, 3.0 / 4.0, 7.0 / 8.0}}};
*/

/*
SymPy code used to generate the inverse matrices

import sympy as sp

KFMAX = 8

def cpp_rational(value):
    value = sp.Rational(value)
    if value == 0:
        return "0.0"
    numerator, denominator = value.as_numer_denom()
    return f"{int(numerator)}.0 / {int(denominator)}.0"

blocks = []
for multiplicity in range(2, KFMAX + 1):
    system = sp.ones(multiplicity)
    system[:2, :2] = sp.Matrix([[2, 0], [0, 2]])

    for row_index in range(multiplicity - 2):
        row = [1] * multiplicity
        row[row_index + 1] = 0
        row[row_index + 2] = 2
        system[row_index + 2, :] = sp.Matrix(1, multiplicity, row)

    inverse = system.inv()
    rows = ["{" + ", ".join(cpp_rational(value) for value in inverse.row(row)) +
"}" for row in range(multiplicity)] blocks.append("    {" + ",\n ".join(rows) +
"}")

print("static const std::vector<std::vector<std::vector<double>>> A = {")
print(",\n\n".join(blocks))
print("};")
*/

// (3*Nf-4)-Dim Integral Volume [phi1] x [phi2] x [phi_k1] x [phi_k2] x ... x
// [phi_{Kf-1}]
//                        [pt1] x [pt2] x [kt1] x [kt2] x ... x [kt{Kf-1}]
//                        x [y3] x [y4] x ... x [yN]
//
// For reference, Nf = 4 special case in:
// [REFERENCE: Lebiedowicz, Szczurek, arxiv.org/abs/0912.0190]
double MCentral::BNIntegralVolume() const {
  // Number of central states
  const unsigned int Kf = state.lts.decaytree.size();

  // Forward leg integration
  const double forward_volume = ForwardVolume();

  // Mass-conditioned transverse and rapidity map volume
  double KT_vol = state.lts.central_phase_space_generated_jacobian;
  KT_vol *= std::pow(state.gcuts.rap_max - state.gcuts.rap_min, Kf);

  return KT_vol * forward_volume;
}

// Compute direct 2 -> N Lorentz invariant phase-space density
double MCentral::BNPhaseSpaceWeight() const {
  // Longitudinal energy-delta Jacobian after solving the two forward pz values
  const double J =
      1.0 / std::abs(state.lts.pfinal[1].Pz() / state.lts.pfinal[1].E() -
                     state.lts.pfinal[2].Pz() / state.lts.pfinal[2].E());

  // Number of final states
  const unsigned int Nf = state.lts.decaytree.size() + 2;
  const unsigned int Kf = Nf - 2;

  // With q_i = p_i - p_{i+1} and fixed sum_i p_i, the map closes through K
  // momenta This gives d^{K-1}p = d^{K-1}q / K per x- and y-axis, thus 1/K^2
  const double central_transverse_map_det = 1.0 / pow2(static_cast<double>(Kf));

  // Each central rapidity measure contributes d^3p/(2E) = d^2pT dy/2
  const double central_rapidity_measure = 1.0 / std::pow(2.0, Kf);

  const double factor =
      central_transverse_map_det * central_rapidity_measure *
      (1.0 / std::pow(2.0 * PI, 3 * Nf - 4)) *
      (state.lts.pfinal[1].Pt() / (2.0 * state.lts.pfinal[1].E())) *
      (state.lts.pfinal[2].Pt() / (2.0 * state.lts.pfinal[2].E())) * J;
  return factor;
}

} // namespace gra
