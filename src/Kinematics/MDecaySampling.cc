// Cascade phase-space proposal densities and identical-particle mixtures
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Kinematics/MDecaySampling.h"

#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/Spin/MHelicityDecay.h"
#include "Graniitti/Tech/MException.h"

namespace gra::decay {
// Compute the invariant-mass proposal density independently of the amplitude
// Fixed massive poles carry the integrated narrow-width normalization
// [REFERENCE: PDG, Review of Particle Physics, Kinematics, Eqs. (49.10)-(49.13)]
double MassDensity(const MDecayBranch &branch) {
  if (branch.mass_proposal == MassProposal::None) { return 1.0; }
  double density = 1.0 / branch.mass_proposal_norm;
  if (branch.mass_proposal != MassProposal::Uniform && branch.p.width > 0.0) {
    density /= math::pow2(branch.p4.M2() - math::pow2(branch.p.mass)) + math::pow2(branch.p.mass * branch.p.width);
  }
  if (!std::isfinite(density) || density <= 0.0) {
    throw PhaseSpaceFailure("decay::MassDensity: invalid invariant-mass proposal density");
  }
  return density;
}

namespace {

// Compute the generated invariant-mass density for one cascade branch
double InternalMassProposalDensity(const MDecayBranch &branch) {
  if (branch.legs.empty()) { return 1.0; }

  double out = MassDensity(branch);
  for (const auto &leg : branch.legs) { out *= InternalMassProposalDensity(leg); }
  return out;
}

// Compute the mass proposal density without propagators from unsampled roots
double InternalMassProposalDensity(const std::vector<MDecayBranch> &tree) {
  double out = 1.0;
  for (const auto &branch : tree) { out *= InternalMassProposalDensity(branch); }
  return out;
}

// Compute a numerical tolerance for one invariant mass squared support
double InvariantMassSupportTolerance(const MDecayBranch &branch, double lower2, double upper2) {
  const double support_scale  = std::max({1.0, std::abs(lower2), std::abs(upper2)});
  const double roundoff_scale = std::max({support_scale, branch.p4.E() * branch.p4.E(), branch.p4.P3mod2()});
  return std::max(1e-12 * support_scale, 128.0 * std::numeric_limits<double>::epsilon() * roundoff_scale);
}

// Compute whether one crossed branch remains inside every generated mass
// proposal support
bool InternalMassProposalSupports(const MDecayBranch &branch) {
  const double mass2 = branch.p4.M2();
  if (!std::isfinite(mass2)) { return false; }

  const bool continuous =
      (branch.mass_proposal == gra::MassProposal::BreitWigner) || (branch.mass_proposal == gra::MassProposal::Uniform);
  if (continuous) {
    if (!std::isfinite(branch.mass_proposal_min2) || !std::isfinite(branch.mass_proposal_max2) ||
        !(branch.mass_proposal_max2 > branch.mass_proposal_min2)) {
      throw PhaseSpaceFailure(
          "Stable-leaf symmetrization: invalid "
          "intermediate mass-proposal support");
    }
    const double tolerance =
        InvariantMassSupportTolerance(branch, branch.mass_proposal_min2, branch.mass_proposal_max2);
    if (mass2 < branch.mass_proposal_min2 - tolerance || mass2 > branch.mass_proposal_max2 + tolerance) {
      return false;
    }
  } else if ((branch.mass_proposal == gra::MassProposal::Fixed)) {
    if (!std::isfinite(branch.mass_proposal_min2) || !std::isfinite(branch.mass_proposal_max2)) {
      throw PhaseSpaceFailure(
          "Stable-leaf symmetrization: invalid fixed "
          "intermediate mass support");
    }
    const double tolerance =
        InvariantMassSupportTolerance(branch, branch.mass_proposal_min2, branch.mass_proposal_max2);
    if (std::abs(branch.mass_proposal_max2 - branch.mass_proposal_min2) > tolerance) {
      throw PhaseSpaceFailure(
          "Stable-leaf symmetrization: fixed "
          "intermediate mass support is not a point");
    }
    if (std::abs(mass2 - branch.mass_proposal_min2) > tolerance) { return false; }
  }
  for (const auto &leg : branch.legs) {
    if (!InternalMassProposalSupports(leg)) { return false; }
  }
  return true;
}

// Compute whether all crossed intermediate masses lie inside the sampled
// proposal support
bool InternalMassProposalSupports(const std::vector<MDecayBranch> &tree) {
  for (const auto &branch : tree) {
    if (!InternalMassProposalSupports(branch)) { return false; }
  }
  return true;
}

// Compute whether crossed top-level branches remain inside the central phase
// space rapidity range
bool CentralRapidityProposalSupports(const gra::LORENTZSCALAR &lts, const std::vector<MDecayBranch> &tree) {
  if (lts.central_phase_space_mode != gra::CentralPhaseSpaceMode::Central) { return true; }
  if (!std::isfinite(lts.central_phase_space_rap_min_cm) || !std::isfinite(lts.central_phase_space_rap_max_cm) ||
      !(lts.central_phase_space_rap_min_cm < lts.central_phase_space_rap_max_cm)) {
    throw PhaseSpaceFailure(
        "Stable-leaf symmetrization: invalid "
        "central phase space rapidity support");
  }

  const gra::M4Vec beamsum        = lts.pbeam1 + lts.pbeam2;
  const double     collision_mass = beamsum.M();
  if (!std::isfinite(collision_mass) || !(collision_mass > 0.0)) {
    throw PhaseSpaceFailure(
        "Stable-leaf symmetrization: invalid "
        "collision CM boost");
  }

  const double scale =
      std::max({1.0, std::abs(lts.central_phase_space_rap_min_cm), std::abs(lts.central_phase_space_rap_max_cm)});
  const double tolerance = 1e-12 * scale;
  for (const auto &branch : tree) {
    gra::M4Vec branch_cm = branch.p4;
    gra::kinematics::LorentzBoost(beamsum, collision_mass, branch_cm, -1);
    const double rapidity_cm = branch_cm.Rap();
    if (!std::isfinite(rapidity_cm) || rapidity_cm < lts.central_phase_space_rap_min_cm - tolerance ||
        rapidity_cm > lts.central_phase_space_rap_max_cm + tolerance) {
      return false;
    }
  }
  return true;
}

// Compute whether the central map has a history-dependent Jacobian
bool HasCentralProposal(const gra::LORENTZSCALAR &lts) {
  return lts.central_phase_space_mode == gra::CentralPhaseSpaceMode::Central ||
         lts.central_phase_space_mode == gra::CentralPhaseSpaceMode::Factorized ||
         lts.central_phase_space_mode == gra::CentralPhaseSpaceMode::Collinear ||
         lts.central_phase_space_mode == gra::CentralPhaseSpaceMode::HardDiffraction;
}

// Compute the inverse-density Jacobian of the central proposal in each history
double CentralProposalJacobian(const gra::LORENTZSCALAR &lts, const std::vector<MDecayBranch> &tree) {
  if (lts.central_phase_space_mode == CentralPhaseSpaceMode::Central) {
    const double multiplicity = tree.size();
    const M4Vec  mean         = lts.pfinal[0] / multiplicity;
    double       rho2         = 0.0;
    double       mass_sum     = 0.0;
    for (const auto &branch : tree) {
      rho2 += (branch.p4 - mean).Pt2();
      mass_sum += branch.p4.M();
    }
    const double scale = std::max(mass_sum, 0.01 * lts.central_phase_space_mass_max);
    return kinematics::TransverseLogJacobian(rho2, lts.central_phase_space_transverse_radius2,
                                             math::pow2(scale) / multiplicity, tree.size());
  }
  if (!HasCentralProposal(lts)) { return 1.0; }
  if (!std::isfinite(lts.central_phase_space_mass_cut_min) || !std::isfinite(lts.central_phase_space_mass_max) ||
      !std::isfinite(lts.central_phase_space_mass_margin) || lts.central_phase_space_mass_cut_min < 0.0 ||
      !(lts.central_phase_space_mass_max > 0.0) || lts.central_phase_space_mass_margin < 0.0) {
    throw PhaseSpaceFailure(
        "Stable-leaf symmetrization: invalid "
        "central-mass parameters");
  }

  // A single root keeps the generated central mass and its original mass window
  if (tree.size() == 1) {
    const double jacobian = lts.central_phase_space_generated_jacobian;
    if (!std::isfinite(jacobian) || jacobian <= 0.0) {
      throw PhaseSpaceFailure(
          "Stable-leaf symmetrization: invalid generated "
          "central proposal Jacobian");
    }
    return jacobian;
  }

  double branch_mass_sum = 0.0;
  for (const auto &branch : tree) {
    double mass = 0.0;
    if (!branch.p4.PhysicalMass(mass)) {
      throw PhaseSpaceFailure("Stable-leaf symmetrization: invalid crossed top-level mass");
    }
    branch_mass_sum += mass;
  }

  const double mass_min =
      std::max(branch_mass_sum + lts.central_phase_space_mass_margin, lts.central_phase_space_mass_cut_min);
  const double mass_max = lts.central_phase_space_mass_max;
  if (!(mass_min < mass_max)) { return 0.0; }

  const double mass2 = lts.pfinal[0].M2();
  if (!std::isfinite(mass2) || !(mass2 > 0.0)) {
    throw PhaseSpaceFailure("Stable-leaf symmetrization: invalid central invariant mass");
  }
  const double lower2    = gra::math::pow2(mass_min);
  const double upper2    = gra::math::pow2(mass_max);
  const double tolerance = 64.0 * std::numeric_limits<double>::epsilon() * std::max({1.0, std::abs(mass2), upper2});
  if (mass2 < lower2 - tolerance || mass2 > upper2 + tolerance) { return 0.0; }
  return 2.0 * mass2 * (std::log(mass_max) - std::log(mass_min));
}

// Compute the deterministic phase-space weight for one generated branch split
// Recursive N-body decays are expressed as chained two-body weights
//
double BranchDecayPhaseSpace(const MDecayBranch &branch) {
  if (branch.legs.empty()) { return 1.0; }

  double M0 = 0.0;
  if (!branch.p4.PhysicalMass(M0)) { return 0.0; }
  std::vector<double> masses;
  masses.reserve(branch.legs.size());
  for (const auto &leg : branch.legs) {
    double mass = 0.0;
    if (!leg.p4.PhysicalMass(mass)) { return 0.0; }
    masses.push_back(mass);
  }

  if (branch.legs.size() == 2) {
    const double pnorm = gra::kinematics::DecayMomentum(M0, masses[0], masses[1]);
    return gra::kinematics::dPhi2(M0, pnorm);
  }

  if (branch.legs.size() == 3) {
    double m12 = 0.0;
    if (!(branch.legs[1].p4 + branch.legs[2].p4).PhysicalMass(m12)) { return 0.0; }
    const double volume = M0 - masses[0] - masses[1] - masses[2];
    if (volume <= 0.0) { return 0.0; }
    const double pnorm0 = gra::kinematics::DecayMomentum(M0, masses[0], m12);
    const double pnorm1 = gra::kinematics::DecayMomentum(m12, masses[1], masses[2]);
    return volume * (m12 / gra::math::PI) * gra::kinematics::dPhi2(M0, pnorm0) * gra::kinematics::dPhi2(m12, pnorm1);
  }

  const std::size_t N = branch.legs.size();
  if (N < 2) { return 0.0; }

  std::vector<double> M_eff(N, 0.0);
  M_eff[0]     = masses[0];
  M_eff[N - 1] = M0;
  for (std::size_t i = 1; i + 1 < N; ++i) {
    M4Vec sum(0.0, 0.0, 0.0, 0.0);
    for (std::size_t j = 0; j <= i; ++j) { sum += branch.legs[j].p4; }
    if (!sum.PhysicalMass(M_eff[i])) { return 0.0; }
  }

  double mass_sum = 0.0;
  for (const double m : masses) { mass_sum += m; }
  const double delta = M0 - mass_sum;
  if (delta <= 0.0) { return 0.0; }

  double pnorm_product = 1.0;
  for (std::size_t i = 1; i < N; ++i) {
    const double pnorm = gra::kinematics::DecayMomentum(M_eff[i], M_eff[i - 1], masses[i]);
    if (pnorm <= 0.0) { return 0.0; }
    pnorm_product *= pnorm;
  }

  return pnorm_product / (2.0 * std::pow(2.0 * gra::math::PI, 2 * N - 3)) * std::pow(delta, N - 2) /
         gra::math::factorial(N - 2) / M0;
}

// Compute the product of deterministic recursive phase-space weights for one branch
// Every internal split contributes its generated phase-space factor
//
double InternalDecayPhaseSpaceProduct(const MDecayBranch &branch) {
  if (branch.legs.empty()) { return 1.0; }

  double out = BranchDecayPhaseSpace(branch);
  for (const auto &leg : branch.legs) { out *= InternalDecayPhaseSpaceProduct(leg); }
  return out;
}

// Compute the product of deterministic recursive phase-space weights for a tree
// The full decay proposal is the product over branch splits
//
double InternalDecayPhaseSpaceProduct(const std::vector<MDecayBranch> &tree) {
  double out = 1.0;
  for (const auto &branch : tree) { out *= InternalDecayPhaseSpaceProduct(branch); }
  return out;
}

// Compute the deterministic phase-space weight for the top-level central split
// The X -> top-level branches split is part of the proposal density
//
double RootDecayPhaseSpace(const gra::LORENTZSCALAR &lts, const std::vector<MDecayBranch> &tree) {
  if (lts.central_phase_space_mode == gra::CentralPhaseSpaceMode::Central) { return 1.0; }
  if (!lts.PS_active || lts.DW.GetN() <= 0.0 || tree.size() <= 1) { return 1.0; }

  // Active F kinematics uses massive RAMBO for the direct split above three bodies
  if (lts.central_phase_space_mode == gra::CentralPhaseSpaceMode::Factorized && tree.size() > 3) {
    std::vector<M4Vec> momenta;
    momenta.reserve(tree.size());
    for (const auto &branch : tree) {
      momenta.push_back(kinematics::BoostToRestFrame(branch.p4, lts.pfinal[0], "RootDecayPhaseSpace RAMBO"));
    }
    return kinematics::RamboWeight(lts.pfinal[0].M(), momenta);
  }

  MDecayBranch root;
  root.p4   = lts.pfinal[0];
  root.legs = tree;
  return BranchDecayPhaseSpace(root);
}

// Compute all deterministic phase-space factors included in the active proposal
// This is the denominator for coherent leaf-permutation mixtures
//
double StableLeafProposalPhaseSpace(const gra::LORENTZSCALAR &lts, const std::vector<MDecayBranch> &tree) {
  return RootDecayPhaseSpace(lts, tree) * InternalDecayPhaseSpaceProduct(tree);
}


}  // namespace

// Compute the stable-leaf mixture density used by coherent cascade symmetrization
// The phase-space weight divides by the mean of q_mass J_generated W_generated / (J_history W_history)
//
double MixtureDensity(gra::LORENTZSCALAR &lts, const std::vector<MDecayBranch> &reference_tree) {
  const auto   terms = spin::StableLeafAmplitudeTrees(lts, reference_tree);
  const double reference_phase_space =
      (lts.decay_symmetry_proposal_active && lts.decay_symmetry_proposal_phase_space > 0.0)
          ? lts.decay_symmetry_proposal_phase_space
          : StableLeafProposalPhaseSpace(lts, reference_tree);
  if (!std::isfinite(reference_phase_space) || reference_phase_space <= 0.0) {
    throw PhaseSpaceFailure("Stable-leaf symmetrization: invalid reference proposal phase-space");
  }
  double generated_jacobian = 1.0;
  if (HasCentralProposal(lts)) {
    generated_jacobian = lts.central_phase_space_generated_jacobian;
    if (!std::isfinite(generated_jacobian) || generated_jacobian <= 0.0) {
      throw PhaseSpaceFailure(
          "Stable-leaf symmetrization: invalid generated "
          "central proposal Jacobian");
    }
  }

  double mixture_density = 0.0;
  for (const auto &term_tree : terms) {
    // Crossed histories outside the exact generated window have zero proposal
    // density
    if (!CentralRapidityProposalSupports(lts, term_tree.tree) || !InternalMassProposalSupports(term_tree.tree)) {
      continue;
    }

    const double term_jacobian = CentralProposalJacobian(lts, term_tree.tree);
    if (term_jacobian <= 0.0) { continue; }

    const double term_phase_space = StableLeafProposalPhaseSpace(lts, term_tree.tree);
    if (!std::isfinite(term_phase_space) || term_phase_space <= 0.0) {
      throw PhaseSpaceFailure("Stable-leaf symmetrization: invalid mixture proposal phase-space");
    }

    const double spectral_density = InternalMassProposalDensity(term_tree.tree);
    mixture_density += spectral_density * reference_phase_space / term_phase_space * generated_jacobian / term_jacobian;
  }
  mixture_density /= static_cast<double>(terms.size());

  if (!std::isfinite(mixture_density) || mixture_density <= 0.0) {
    throw PhaseSpaceFailure("Stable-leaf symmetrization: invalid mixture proposal density");
  }

  return mixture_density;
}


}  // namespace gra::decay
