// Forward proton excitation and fragmentation state
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Process/MForwardExcitation.h"

#include <array>
#include <cmath>
#include <vector>

#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/MGlobals.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Nuclear/MBeam.h"
#include "Graniitti/Process/MProcess.h"
#include "Graniitti/Regge/MFragment.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MException.h"

using gra::aux::indices;
using gra::math::msqrt;
using gra::math::pow2;

namespace gra {

namespace forward {

// Compute and validate the invariant mass interval of an excited proton
std::array<double, 2> MassBounds(const MProcessState &state, const GENCUT &cuts) {
  if (!state.nstar_param) { throw std::invalid_argument("forward::MassBounds: missing PARAM_NSTAR"); }
  if (!std::isfinite(state.lts.s) || !(state.lts.s > 0.0) ||
      !std::isfinite(cuts.XI_min) || cuts.XI_min < 0.0 || !std::isfinite(cuts.XI_max) ||
      !std::isfinite(state.nstar_param->mass_margin) || state.nstar_param->mass_margin < 0.0) {
    throw std::invalid_argument("forward::MassBounds: invalid forward mass bounds");
  }
  const double p = state.nstar_param->single_side_prob;
  if (state.excitation == 1 && (!std::isfinite(p) || !(p > 0.0 && p < 1.0))) {
    throw std::invalid_argument("forward::MassBounds: invalid beam sampling probability");
  }
  const double threshold = PDG::mp + PDG::mpi + state.nstar_param->mass_margin;
  const double lower = std::max(pow2(threshold), cuts.XI_min * state.lts.s);
  const double upper = cuts.XI_max * state.lts.s;
  if (!std::isfinite(lower) || !std::isfinite(upper) || upper <= lower) {
    throw std::invalid_argument("forward::MassBounds: empty or non-finite forward mass interval");
  }
  return {lower, upper};
}

// Clear generated forward branches after event serialization
void ClearBranches(LORENTZSCALAR &lts) noexcept {
  lts.excite1 = lts.excite2 = false;
  lts.decayforward1 = MDecayBranch();
  lts.decayforward2 = MDecayBranch();
}

// Clear all event-local forward excitation and diffractive variables
void Reset(LORENTZSCALAR &lts) noexcept {
  ClearBranches(lts);
  lts.hard_diff1 = lts.hard_diff2 = false;
  lts.has_xi1 = lts.has_xi2 = false;
  lts.xi1 = lts.xi2 = 0.0;
  lts.diff_xi1 = lts.diff_xi2 = 0.0;
  lts.diff_beta1 = lts.diff_beta2 = 0.0;
  lts.diff_t1 = lts.diff_t2 = 0.0;
  lts.diff_phi1 = lts.diff_phi2 = 0.0;
  lts.diff_xhard1 = lts.diff_xhard2 = 0.0;
}

}  // namespace forward

// -----------------------------------------------------------------------

// Proton continuum excitation
//
// Compute true if OK
// Compute false if FAILS
//
bool MProcess::ExciteContinuum(const M4Vec &nstar, gra::MDecayBranch &forward, double Q2_scale, int B, int Q) {
  if (!state.nstar_param || !std::isfinite(Q2_scale) || Q2_scale <= 0.0) { return false; }
  const auto &param    = *state.nstar_param;
  const auto &cylinder = param.cylinder;

  // Average multiplicity
  const double avgN =
      cylinder.mult_a + cylinder.mult_b / cylinder.charged_prob * std::log1p(Q2_scale / cylinder.q2_ref);
  unsigned int outertrials = 0;

  const double M = nstar.M();

  while (true) {  // OUTER

    // Draw fluctuating number of particles
    int          N                   = 0;
    unsigned int multiplicity_trials = 0;
    while (multiplicity_trials < cylinder.mult_trials) {
      N = state.random.PoissonRandom(avgN);
      N = std::max({N, 2, std::abs(B), std::abs(Q)});
      if (nstar.M() > static_cast<double>(N) * cylinder.mass_unit) { break; }
      ++multiplicity_trials;
    }
    if (multiplicity_trials == cylinder.mult_trials) { return false; }

    // Boundary conditions
    N = (N < 2) ? 2 : N;                       // Tube kinematics requires at least 2
    N = (N >= std::abs(B)) ? N : std::abs(B);  // Baryon number
    N = (N >= std::abs(Q)) ? N : std::abs(Q);  // Charge conservation

    // --------------------------------------------------------------------
    // Pick particles

    std::vector<double> mass;
    std::vector<int>    pdgcode;
    if (!MFragment::PickParticles(M, N, B, 0, Q, mass, pdgcode, state.lts.PDG, state.random, param)) {
      ++outertrials;
      if (outertrials > cylinder.outer_trials) {
        return false;  // too many failures
      } else {
        continue;  // try again!
      }
    }

    // --------------------------------------------------------------------
    // Now decay >>

    std::vector<M4Vec> p4;

    const double q = std::pow(N, cylinder.q_power);
    double       T = cylinder.T;
    if (state.multipomeron_chain.size() > 0) {  // ND temperature scales as (1/b^2)^impact_power
      const double impact2 = state.multipomeron_impact_parameter * state.multipomeron_impact_parameter;
      if (!std::isfinite(impact2) || !(impact2 > 0.0)) { return false; }
      T *= std::pow(1.0 / impact2, cylinder.impact_power);
    }
    const double W =
        MFragment::TubeFragment(nstar, M, mass, p4, q, T, cylinder.exp_lambda, cylinder.max_pt, state.random, param);

    if (W <= 0) {
      ++outertrials;
      if (outertrials > cylinder.outer_trials) {
        return false;  // too many failures
      } else {
        continue;  // try again!
      }
    } else {
      // Re-check kinematics
      bool fail = false;
      for (const auto &i : indices(p4)) {
        if (std::isnan(p4[i].Rap())) {
          ++outertrials;
          if (outertrials > cylinder.outer_trials) {
            return false;  // too many failures
          }
          fail = true;
          break;
        }
      }
      if (fail) { continue; }

      // Get corresponding PDG particles
      std::vector<gra::MParticle> p(p4.size());
      for (const auto &i : indices(p)) { p[i] = state.lts.PDG.FindByPDG(pdgcode[i]); }

      // Create decaytree
      if (!BranchForwardSystem(p4, p, nstar, forward)) {
        ++outertrials;
        if (outertrials > cylinder.outer_trials) { return false; }
        continue;
      }
      break;
    }
  }
  return true;
}

// Proton low mass (N* like) 2-body and 3-body excitation
//
//
bool MProcess::ExciteNstar(const M4Vec &nstar, gra::MDecayBranch &forward, const MParticle &pbeam) {
  if (!state.nstar_param) { return false; }
  // Find state.random decaymode
  std::vector<int> pdgcodes;
  MFragment::NstarDecayTable(pbeam.chargeX3 / 3, nstar.M(), pdgcodes, state.random, *state.nstar_param);

  // Get corresponding PDG particles
  std::vector<gra::MParticle> p(pdgcodes.size());
  std::vector<double>         mass(pdgcodes.size(), 0.0);
  for (const auto &i : indices(pdgcodes)) {
    p[i]    = state.lts.PDG.FindByPDG(pdgcodes[i]);
    mass[i] = p[i].mass;
  }

  // Products 4-momenta
  std::vector<M4Vec> p4;

  // Do the 2 or 3-body isotropic decay
  gra::kinematics::MCW weight = gra::kinematics::InvalidPhaseSpacePoint();
  if (mass.size() == 2) { weight = gra::kinematics::TwoBodyPhaseSpace(nstar, nstar.M(), mass, p4, state.random); }
  if (mass.size() == 3) {
    weight = gra::kinematics::ThreeBodyPhaseSpace(nstar, nstar.M(), mass, p4, state.nstar_param->fewbody.unweight,
                                                  state.random);
  }
  if (weight.GetW() < 0.0 || p4.size() != p.size()) { return false; }

  // Create decaytree
  return BranchForwardSystem(p4, p, nstar, forward);
}

// Build one showerable color-singlet quark-diquark forward system
bool MProcess::ExciteString(const M4Vec &nstar, gra::MDecayBranch &forward, const MParticle &pbeam, int color_tag) {
  if (!state.nstar_param) { return false; }
  const auto &param        = state.nstar_param->string;
  const int   beam_sign    = math::sign(pbeam.pdg);
  const int   beam_abs_pdg = std::abs(pbeam.pdg);
  if (beam_abs_pdg != PDG::PDG_p && beam_abs_pdg != PDG::PDG_n) {
    throw std::invalid_argument("MProcess::ExciteString: unsupported beam PDG " + std::to_string(pbeam.pdg));
  }
  const bool proton        = beam_abs_pdg == PDG::PDG_p;
  const bool first_valence = state.random.U(0.0, 1.0) < param.first_valence_prob;

  int    quark        = 0;
  int    diquark      = 0;
  double diquark_mass = 0.0;
  if (proton && first_valence) {
    quark        = 2;
    diquark      = 2101;
    diquark_mass = param.scalar_diquark_mass;
  } else if (proton) {
    quark        = 1;
    diquark      = 2203;
    diquark_mass = param.vector_diquark_mass;
  } else if (first_valence) {
    quark        = 1;
    diquark      = 2101;
    diquark_mass = param.scalar_diquark_mass;
  } else {
    quark        = 2;
    diquark      = 1103;
    diquark_mass = param.vector_diquark_mass;
  }
  quark *= beam_sign;
  diquark *= beam_sign;

  MParticle quark_particle = state.lts.PDG.FindByPDG(quark);
  MParticle diquark_particle;
  diquark_particle.pdg      = diquark;
  diquark_particle.name     = "forward-diquark";
  diquark_particle.mass     = diquark_mass;
  diquark_particle.chargeX3 = pbeam.chargeX3 - quark_particle.chargeX3;
  diquark_particle.spinX2   = (std::abs(diquark) % 10 == 1) ? 0 : 2;
  diquark_particle.color    = (diquark > 0) ? -3 : 3;

  if (quark > 0) {
    quark_particle.color_flow.flow1   = color_tag;
    diquark_particle.color_flow.flow2 = color_tag;
  } else {
    quark_particle.color_flow.flow2   = color_tag;
    diquark_particle.color_flow.flow1 = color_tag;
  }

  const std::vector<double> masses = {quark_particle.mass, diquark_particle.mass};
  std::vector<M4Vec>        momenta;
  const auto weight = gra::kinematics::TwoBodyPhaseSpace(nstar, nstar.M(), masses, momenta, state.random);
  if (weight.GetW() < 0.0 || momenta.size() != 2) { return false; }

  return BranchForwardSystem(momenta, {quark_particle, diquark_particle}, nstar, forward);
}

// Build one complete forward branch and publish it only after validation
bool MProcess::BranchForwardSystem(const std::vector<M4Vec> &p4, const std::vector<MParticle> &p, const M4Vec &nstar,
                                   gra::MDecayBranch &forward) {
  if (p4.size() != p.size()) {
    throw std::invalid_argument(
        "MProcess::BranchForwardSystem: momentum and "
        "particle dimensions differ");
  }
  /*
    std::vector<bool> isstable(N, false);
    MFragment::GetDecayStatus(pdgcode, isstable);
  */

  // Construct one local decay tree
  MDecayBranch candidate;
  candidate.p4 = nstar;

  // Daughters
  candidate.legs.resize(p.size());
  candidate.depth = 0;

  for (const auto &i : indices(p)) {
    if (!std::isfinite(p4[i].Px()) || !std::isfinite(p4[i].Py()) || !std::isfinite(p4[i].Pz()) ||
        !std::isfinite(p4[i].E()) || !(p4[i].E() > 0.0)) {
      return false;
    }

    // Decay particle
    gra::MDecayBranch branch;
    branch.p  = p[i];
    branch.p4 = p4[i];

    // Treat pi0 -> yy (BR ~ 100%)
    if (branch.p.pdg == PDG::PDG_pi0) {
      // Decay
      const std::vector<double> md = {0.0, 0.0};
      std::vector<M4Vec>        pd;
      const auto weight = gra::kinematics::TwoBodyPhaseSpace(branch.p4, branch.p4.M(), md, pd, state.random);
      if (weight.GetW() < 0.0 || pd.size() != 2) { return false; }

      // Add gamma legs
      branch.legs.resize(2);
      for (std::size_t k = 0; k < 2; ++k) {
        gra::MDecayBranch decaybranch;
        decaybranch.p  = state.lts.PDG.FindByPDG(PDG::PDG_gamma);
        decaybranch.p4 = pd[k];

        branch.legs[k]       = decaybranch;
        branch.legs[k].depth = 2;
      }
    }

    // Treat rho0 -> pi+ pi- (BR ~ 100%)
    else if (branch.p.pdg == PDG::PDG_rho0) {
      // Decay
      const std::vector<double> md  = {PDG::mpi, PDG::mpi};
      const std::vector<int>    pdg = {PDG::PDG_pip, PDG::PDG_pim};
      std::vector<M4Vec>        pd;
      const auto weight = gra::kinematics::TwoBodyPhaseSpace(branch.p4, branch.p4.M(), md, pd, state.random);
      if (weight.GetW() < 0.0 || pd.size() != 2) { return false; }

      // Add pion legs
      branch.legs.resize(2);
      for (std::size_t k = 0; k < 2; ++k) {
        gra::MDecayBranch decaybranch;
        decaybranch.p  = state.lts.PDG.FindByPDG(pdg[k]);
        decaybranch.p4 = pd[k];

        branch.legs[k]       = decaybranch;
        branch.legs[k].depth = 2;
      }
    }

    // Here, we could add other MAJOR resonant decays
    // ...

    // Add leg
    candidate.legs[i]       = branch;
    candidate.legs[i].depth = 1;
  }

  if (!candidate.legs.empty()) {
    M4Vec daughter_sum;
    for (const auto &leg : candidate.legs) { daughter_sum += leg.p4; }
    if (!math::CheckEMC(nstar - daughter_sum)) { return false; }
  }

  forward = std::move(candidate);
  return true;
}


// Fragment forward systems in central production
//
bool MProcess::CEPForwardFragment() {
  if (!state.nstar_param) { throw std::logic_error("MProcess::CEPForwardFragment: missing PARAM_NSTAR snapshot"); }
  const auto fragment = state.beamfrag;

  std::array<MDecayBranch, 2>    branches   = {state.lts.decayforward1, state.lts.decayforward2};
  const std::array<bool, 2>      excited    = {state.lts.excite1, state.lts.excite2};
  const std::array<M4Vec, 2>     systems    = {state.lts.pfinal[1], state.lts.pfinal[2]};
  const std::array<MParticle, 2> beams      = {state.lts.beam1, state.lts.beam2};
  constexpr std::array<int, 2>   color_tags = {701, 702};

  for (const auto &side : indices(branches)) {
    if (!excited[side] || !branches[side].legs.empty()) { continue; }

    const int baryon  = math::sign(beams[side].pdg);
    const int charge  = beams[side].chargeX3 / 3;
    bool      success = false;
    if (fragment == BeamFragType::None) {
      success = BranchForwardSystem({}, {}, systems[side], branches[side]);
    } else if (fragment == BeamFragType::FewBody) {
      success = ExciteNstar(systems[side], branches[side], beams[side]);
    } else if (fragment == BeamFragType::Cylinder) {
      success = ExciteContinuum(systems[side], branches[side], systems[side].M2(), baryon, charge);
    } else {
      success = ExciteString(systems[side], branches[side], beams[side], color_tags[side]);
    }
    if (!success) { return false; }
  }

  state.lts.decayforward1 = std::move(branches[0]);
  state.lts.decayforward2 = std::move(branches[1]);
  return true;
}


// Forward excitation mass sampling
void MProcess::SampleForwardMasses(std::vector<double> &mvec, const std::vector<double> &randvec) {
  mvec = {state.lts.beam1.mass, state.lts.beam2.mass};
  state.lts.forward_emd = {};
  state.lts.upc_excitation.reset();

  if (state.nuclear_final.has_value()) {
    HepMC3::GenEvent beams(HepMC3::Units::GEV, HepMC3::Units::MM);
    auto first = std::make_shared<HepMC3::GenParticle>(aux::M4Vec2HepMC3(state.lts.pbeam1), state.lts.beam1.pdg, 4);
    auto second = std::make_shared<HepMC3::GenParticle>(aux::M4Vec2HepMC3(state.lts.pbeam2), state.lts.beam2.pdg, 4);
    first->set_generated_mass(state.lts.beam1.mass);
    second->set_generated_mass(state.lts.beam2.mass);
    beams.set_beam_particles(first, second);
    const auto mass = state.nuclear_final->SampleMasses(beams, state.random);
    mvec.assign(mass.mass.cbegin(), mass.mass.cend());
    state.lts.forward_emd = mass.emd;
    state.lts.upc_excitation = mass.channel;
  }

  state.lts.excite1 = false;
  state.lts.excite2 = false;

  if (state.excitation != 0) {
    if (randvec.size() < static_cast<std::size_t>(state.excitation)) {
      throw PhaseSpaceFailure("MProcess::SampleForwardMasses: invalid excitation coordinates");
    }
    for (int i = 0; i < state.excitation; ++i) {
      if (!std::isfinite(randvec[i]) || randvec[i] < 0.0 || randvec[i] > 1.0) {
        throw PhaseSpaceFailure("MProcess::SampleForwardMasses: mass coordinate outside [0,1]");
      }
    }

    // Compute the mass interval at the sampled collision energy
    const double threshold = PDG::mp + PDG::mpi + state.nstar_param->mass_margin;
    M2_f_min = std::max(pow2(threshold), state.gcuts.XI_min * state.lts.s);
    M2_f_max = state.gcuts.XI_max * state.lts.s;
    if (!(M2_f_max > M2_f_min)) {
      throw PhaseSpaceFailure("MProcess::SampleForwardMasses: empty generated forward mass interval");
    }
    log_M2_f_min = std::log(M2_f_min);
    log_M2_f_max = std::log(M2_f_max);

    if (state.excitation == 1) {  // Single

      // log-change of variable
      const double u = log_M2_f_min + (log_M2_f_max - log_M2_f_min) * randvec[0];
      const double r = std::exp(u);

      const bool proton1 = nuclear::IsProton(state.lts.beam1.pdg);
      const bool proton2 = nuclear::IsProton(state.lts.beam2.pdg);
      if (proton1 && (!proton2 || state.random.U(0, 1) < state.nstar_param->single_side_prob)) {
        mvec[0]           = msqrt(r);
        state.lts.excite1 = true;
      } else {
        mvec[1]           = msqrt(r);
        state.lts.excite2 = true;
      }
    } else if (state.excitation == 2) {  // Double

      // log-change of variable
      const double u1 = log_M2_f_min + (log_M2_f_max - log_M2_f_min) * randvec[0];
      const double r1 = std::exp(u1);

      const double u2 = log_M2_f_min + (log_M2_f_max - log_M2_f_min) * randvec[1];
      const double r2 = std::exp(u2);

      mvec[0]           = msqrt(r1);
      mvec[1]           = msqrt(r2);
      state.lts.excite1 = true;
      state.lts.excite2 = true;
    }
  }
}

// Compute the forward-leg integration volume
// V = V_phi V_log(pT) V_log(M2)
double MProcess::ForwardVolume() const {
  // Forward leg phi1,phi2 volumes gives (2pi)^2
  const double PHI_vol = pow2(2.0 * math::PI);

  // Forward leg pt volume with log-change of variables
  const double J_pt = state.lts.pfinal[1].Pt() * state.lts.pfinal[2].Pt();
  const double PT_vol =
      J_pt * pow2(std::log(state.gcuts.forward_pt_max) - std::log(state.gcuts.forward_pt_min + ZERO_EPS));

  if (state.excitation == 0) {
    return PHI_vol * PT_vol;

  } else if (state.excitation == 1) {
    // log-change of variable jacobian
    // \int_a^b f(M2) dM2 = \int_{ln(a)}^{ln(b)} f(exp(u)) * exp(u) du, where u
    // = ln(M2)
    const double J_M2 = (state.lts.excite1) ? state.lts.pfinal[1].M2() : state.lts.pfinal[2].M2();
    return PHI_vol * PT_vol * J_M2 * (log_M2_f_max - log_M2_f_min);

  } else if (state.excitation == 2) {
    // log-change of variable jacobian
    const double J_M2 = state.lts.pfinal[1].M2() * state.lts.pfinal[2].M2();
    return PHI_vol * PT_vol * J_M2 * pow2(log_M2_f_max - log_M2_f_min);
  } else {
    throw std::invalid_argument("MProcess::ForwardVolume: state.excitation variable in invalid state");
  }
}

// Compute the inverse probability of the sampled single dissociation side
double MProcess::DissociationCrossSectionFactor() const {
  const bool soft_single_dissociation     = ProcPtr.ISTATE == "X" && ProcPtr.CHANNEL == "SD";
  const bool hard_single_dissociation     = ProcPtr.ISTATE == "IPp";

  if (soft_single_dissociation || hard_single_dissociation) { return 2.0; }
  if (state.excitation != 1) { return 1.0; }
  if (!state.nstar_param || state.lts.excite1 == state.lts.excite2) {
    throw PhaseSpaceFailure("MProcess::DissociationCrossSectionFactor: missing single dissociation side");
  }
  if (!nuclear::IsProton(state.lts.beam1.pdg) || !nuclear::IsProton(state.lts.beam2.pdg)) { return 1.0; }
  const double p = state.nstar_param->single_side_prob;
  return 1.0 / (state.lts.excite1 ? p : 1.0 - p);
}


}  // namespace gra
