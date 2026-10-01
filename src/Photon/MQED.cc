// QED and EW functions
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <stdexcept>
#include <string>

// Own
#include "Graniitti/Particle/MForm.h"
#include "Graniitti/Particle/MPDG.h"
#include "Graniitti/Photon/MFlux.h"
#include "Graniitti/Photon/MQED.h"
#include "Graniitti/Spin/MHelicity.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MException.h"

using gra::aux::indices;

namespace gra {
namespace qed {

namespace {

constexpr double kHelicityTolerance = 1e-9;

// Compute the particle metadata for one photon-emitting beam leg
const gra::MParticle &Beam(const gra::LORENTZSCALAR &lts, int leg) {
  if (leg == 1) { return lts.beam1; }
  if (leg == 2) { return lts.beam2; }
  throw std::invalid_argument("gra::qed::Beam: leg should be 1 or 2");
}

// Compute the electric charge Q/e for one Standard Model fermion species
double FermionChargeX1(int pdg) {
  const int apdg = std::abs(pdg);
  if (apdg == 11 || apdg == 13 || apdg == 15) { return -1.0; }
  if (apdg == 12 || apdg == 14 || apdg == 16) { return 0.0; }
  if (apdg == 2 || apdg == 4 || apdg == 6) { return 2.0 / 3.0; }
  if (apdg == 1 || apdg == 3 || apdg == 5) { return -1.0 / 3.0; }
  throw std::invalid_argument("gra::qed: unsupported fermion PDG " + std::to_string(pdg));
}

// Compute true for a photon source with an exclusive outgoing spinor
bool ExclusivePhotonEmitter(const gra::LORENTZSCALAR &lts, int leg) {
  if (leg == 1) { return !lts.excite1; }
  if (leg == 2) { return !lts.excite2; }
  throw std::invalid_argument("gra::qed: photon leg should be 1 or 2");
}

// Rotate photon columns into the common exchange-helicity section
void ApplyPhotonSourceSectionPhase(MMatrix<std::complex<double>> &matrix, const std::vector<int> &m_values, double phi,
                                   bool second_exchange_daughter) {
  if (matrix.size_col() != m_values.size()) {
    throw std::invalid_argument(
        "gra::qed::ApplyPhotonSourceSectionPhase: photon m basis has "
        "inconsistent dimensions");
  }
  const double phase_sign = second_exchange_daughter ? -1.0 : 1.0;
  for (const auto &col : indices(m_values)) {
    const double angle = phase_sign * static_cast<double>(m_values[col]) * phi;
    matrix.ScaleColumn(col, std::complex<double>(std::cos(angle), std::sin(angle)));
  }
}

// Build the electric EPA source in the Condon-Shortley photon basis
// The opposite helicity signs follow from epsilon_m dot q_T
// Reversing the beam axis also reverses the transverse photon basis
MMatrix<std::complex<double>> PhotonEPAMatrix(const gra::LORENTZSCALAR &lts, int leg,
                                              const std::vector<std::pair<double, double>> &transitions,
                                              const std::vector<int> &m_values, bool second_exchange_daughter) {
  MMatrix<std::complex<double>>     out(transitions.size(), m_values.size(), 0.0);
  const std::complex<double>        transverse_weight = (leg == 1 ? -1.0 : 1.0) / gra::math::msqrt(2.0);
  std::vector<std::complex<double>> column_factor(m_values.size(), 0.0);
  const double                      phi = ((leg == 1) ? lts.q1 : lts.q2).Phi();
  for (const auto &col : indices(m_values)) {
    if (std::abs(m_values[col]) == 1) { column_factor[col] = static_cast<double>(m_values[col]) * transverse_weight; }
  }
  for (const auto &row : indices(transitions)) {
    if (std::abs(transitions[row].first - transitions[row].second) > kHelicityTolerance) { continue; }
    std::copy(column_factor.cbegin(), column_factor.cend(), out.Row(row).begin());
  }
  ApplyPhotonSourceSectionPhase(out, m_values, phi, second_exchange_daughter);
  return out;
}

// Align each exact photon transition with the signed EPA photon field section
void AlignPhotonSourceRows(MMatrix<std::complex<double>> &source, const MMatrix<std::complex<double>> &reference) {
  if (source.size_row() != reference.size_row() || source.size_col() != reference.size_col()) {
    throw std::invalid_argument("gra::qed::AlignPhotonSourceRows: source dimensions differ");
  }
  for (std::size_t row = 0; row < source.size_row(); ++row) {
    std::complex<double> overlap = 0.0;
    for (const auto &col : indices(source.Row(row))) { overlap += std::conj(reference[row][col]) * source[row][col]; }
    const double magnitude = std::abs(overlap);
    if (magnitude > 0.0 && std::isfinite(magnitude)) { source.ScaleRow(row, std::conj(overlap) / magnitude); }
  }
}

// Build the unnormalized QED photon source from an elastic spin-half current
MMatrix<std::complex<double>> RawQEDPhotonMatrix(const gra::LORENTZSCALAR &lts, int leg,
                                                 const std::vector<std::pair<double, double>> &transitions,
                                                 const std::vector<int>                       &m_values) {
  const std::size_t outgoing_index = (leg == 1) ? 1 : 2;
  if (lts.pfinal.size() <= outgoing_index) {
    throw AmplitudeFailure("gra::qed::RawQEDPhotonMatrix: missing forward final-state leg");
  }
  const gra::M4Vec &q = (leg == 1) ? lts.q1 : lts.q2;

  const gra::MDirac             dirac("DIRAC");
  MMatrix<std::complex<double>> out(transitions.size(), m_values.size(), 0.0);
  const auto                    eps_conj = dirac.MasslessSpin1States(q, "conj", true);
  const auto currents = PhotonCurrentTransitions(lts, leg, transitions);

  for (const auto &row : indices(transitions)) {
    for (const auto &col : indices(m_values)) {
      const int m = m_values[col];
      if (std::abs(m) != 1) { continue; }
      const auto          &eps       = eps_conj[(m + 1) / 2];
      std::complex<double> projected = 0.0;
      for (const auto &mu : dirac.LI) { projected += currents[row][mu] * eps(mu); }
      out[row][col] = projected;
    }
  }
  return out;
}

// Select requested helicity rows from a complete spin-half photon source
MMatrix<std::complex<double>> SelectPhotonTransitionRows(
    const MMatrix<std::complex<double>> &full, const std::vector<std::pair<double, double>> &full_transitions,
    const std::vector<std::pair<double, double>> &selected_transitions) {
  if (full.size_row() != full_transitions.size()) {
    throw std::invalid_argument("gra::qed::SelectPhotonTransitionRows: full matrix row mismatch");
  }

  MMatrix<std::complex<double>> out(selected_transitions.size(), full.size_col(), 0.0);
  for (const auto &row : indices(selected_transitions)) {
    bool found = false;
    for (const auto &full_row : indices(full_transitions)) {
      const bool same_in =
          std::abs(selected_transitions[row].first - full_transitions[full_row].first) < kHelicityTolerance;
      const bool same_out =
          std::abs(selected_transitions[row].second - full_transitions[full_row].second) < kHelicityTolerance;
      if (!same_in || !same_out) { continue; }
      std::copy(full.Row(full_row).begin(), full.Row(full_row).end(), out.Row(row).begin());
      found = true;
      break;
    }
    if (!found) {
      throw std::invalid_argument(
          "gra::qed::SelectPhotonTransitionRows: requested transition is not "
          "spin-half");
    }
  }
  return out;
}

}  // namespace

// Build unnormalized elastic photon currents with propagator and collider phases
std::vector<MDirac::Current> PhotonCurrentTransitions(
    const LORENTZSCALAR &lts, const int leg,
    const std::vector<std::pair<double, double>> &transitions) {
  const std::size_t outgoing_index = leg == 1 ? 1 : leg == 2 ? 2 : 0;
  if (outgoing_index == 0 || lts.pfinal.size() <= outgoing_index) {
    throw AmplitudeFailure(
        "gra::qed::PhotonCurrentTransitions: missing forward final-state leg");
  }
  const M4Vec &incoming = leg == 1 ? lts.pbeam1 : lts.pbeam2;
  const M4Vec &outgoing = lts.pfinal[outgoing_index];
  const M4Vec     &transfer = leg == 1 ? lts.q1 : lts.q2;
  const double     t        = leg == 1 ? lts.t1 : lts.t2;
  if (!std::isfinite(t) || !(t < 0.0)) { throw AmplitudeFailure("PhotonCurrentTransitions: invalid transfer"); }

  const ElasticEmitter emitter   = Emitter(lts, leg, t, "gra::qed::PhotonCurrentTransitions");
  const double inverse_t = 1.0 / t;
  if (!std::isfinite(inverse_t)) {
    throw AmplitudeFailure(
        "gra::qed::PhotonCurrentTransitions: photon propagator overflow");
  }

  const double               azimuth = leg == 1 ? (-transfer).Phi() : transfer.Phi();
  const std::complex<double> propagator =
      -math::zi * e_QED() * emitter.charge * inverse_t;
  const MDirac dirac("DIRAC");
  const auto basis = dirac.DiracPauliCurrentBasis(incoming, outgoing, -transfer,
                                                  emitter.fermion_kind, emitter.mass,
                                                  emitter.F1, emitter.F2);
  std::vector<MDirac::Current> out(transitions.size());
  for (const auto &row : indices(transitions)) {
    const std::size_t in = spin::BinaryHelicityIndexX2(static_cast<int>(2.0 * transitions[row].first));
    const std::size_t final = spin::BinaryHelicityIndexX2(static_cast<int>(2.0 * transitions[row].second));
    out[row] = basis[2 * final + in];
    const std::complex<double> phase =
        MDirac::ElasticSpinHalfColliderCurrentPhase(
            emitter.fermion_kind, leg, transitions[row].first,
            transitions[row].second, azimuth);
    for (auto &component : out[row]) {
      component *= phase * propagator;
    }
  }
  return out;
}

// Compute a lower-index chiral vector current for a fermion-antifermion pair
// J_mu = ubar gamma_mu(cL P_L+cR P_R)v
gra::MDirac::Current FFVCurrent(const gra::MDirac &dirac, const gra::M4Vec &fermion, const gra::M4Vec &antifermion,
                                int helicity_fermion, int helicity_antifermion, std::complex<double> cL,
                                std::complex<double> cR) {

  const auto          ubar = dirac.Bar(dirac.uHelDirac(fermion, helicity_fermion));
  const auto          v    = dirac.vHelDirac(antifermion, helicity_antifermion);
  const auto          vL   = dirac.PL() * v;
  const auto          vR   = dirac.PR() * v;
  gra::MDirac::Spinor chiral_v{};
  for (const auto &i : indices(chiral_v)) { chiral_v[i] = cL * vL[i] + cR * vR[i]; }

  gra::MDirac::Current current{};
  for (const auto &mu : dirac.LI) { current[mu] = gra::BilinearProduct(ubar, dirac.gamma_lo[mu] * chiral_v); }
  return current;
}

// Compute a lower-index photon current for f fbar
gra::MDirac::Current PhotonFermionCurrent(const gra::MDirac &dirac, const gra::M4Vec &fermion,
                                          const gra::M4Vec &antifermion, int fermion_pdg, int helicity_fermion,
                                          int helicity_antifermion, double electric_charge) {
  const std::complex<double> coupling = electric_charge * FermionChargeX1(fermion_pdg);
  return FFVCurrent(dirac, fermion, antifermion, helicity_fermion, helicity_antifermion, coupling, coupling);
}

// Compute neutral-current charge and weak couplings for one fermion species
// gL = T3-Q sin2(thetaW), gR = -Q sin2(thetaW), gV,A = (gL +/- gR)/2
NeutralCurrentCouplings FermionNeutralCurrentCouplings(int pdg, double sin2thetaW) {
  const int apdg = std::abs(pdg);
  if (!std::isfinite(sin2thetaW) || sin2thetaW <= 0.0 || sin2thetaW >= 1.0) {
    throw std::invalid_argument("gra::qed::FermionNeutralCurrentCouplings: invalid sin2thetaW");
  }

  NeutralCurrentCouplings out;
  out.chargeX1 = FermionChargeX1(pdg);
  if (apdg == 11 || apdg == 13 || apdg == 15) {
    out.isospin3 = -0.5;
  } else if (apdg == 12 || apdg == 14 || apdg == 16 || apdg == 2 || apdg == 4 || apdg == 6) {
    out.isospin3 = 0.5;
  } else if (apdg == 1 || apdg == 3 || apdg == 5) {
    out.isospin3 = -0.5;
  } else {
    throw std::invalid_argument("gra::qed::FermionNeutralCurrentCouplings: unsupported PDG " + std::to_string(pdg));
  }

  out.gL = out.isospin3 - out.chargeX1 * sin2thetaW;
  out.gR = -out.chargeX1 * sin2thetaW;
  out.gV = 0.5 * (out.gL + out.gR);
  out.gA = 0.5 * (out.gL - out.gR);
  return out;
}

// Compute a lower-index Z current for f fbar
gra::MDirac::Current ZFermionCurrent(const gra::MDirac &dirac, const gra::M4Vec &fermion, const gra::M4Vec &antifermion,
                                     int fermion_pdg, int helicity_fermion, int helicity_antifermion,
                                     double electric_charge, double sin2thetaW) {
  const NeutralCurrentCouplings coupling = FermionNeutralCurrentCouplings(fermion_pdg, sin2thetaW);
  const double                  sw       = std::sqrt(sin2thetaW);
  const double                  cw       = std::sqrt(1.0 - sin2thetaW);
  const double                  norm     = electric_charge / (sw * cw);
  return FFVCurrent(dirac, fermion, antifermion, helicity_fermion, helicity_antifermion, norm * coupling.gL,
                    norm * coupling.gR);
}

// Validate the photon source steering mode
void ValidatePhotonMode(const std::string &mode, const std::string &context) {
  if (mode == "EPA" || mode == "QED") { return; }
  throw std::invalid_argument(context +
                              ": PARAM_REGGE.PHOTON_VERTEX should be 'EPA' or "
                              "'QED'");
}

// Resolve one supported elastic spin-half QED emitter from particle metadata
ElasticEmitter Emitter(const gra::MParticle &particle, double t, const std::string &context,
                       const gra::form::ParamStore &structure) {
  if (!std::isfinite(t) || t > 0.0) { throw AmplitudeFailure(context + ": photon transfer is not spacelike"); }

  ElasticEmitter emitter;
  emitter.fermion_kind = particle.pdg < 0 ? gra::MDirac::FermionKind::Antiparticle : gra::MDirac::FermionKind::Particle;
  emitter.mass         = particle.mass;
  emitter.charge       = particle.chargeX3 / 3.0;
  emitter.incoming_spin_states = 2;

  const int absolute_pdg = std::abs(particle.pdg);
  if (absolute_pdg == gra::PDG::PDG_p) {
    const double Q2              = -t;
    emitter.F1                   = gra::form::F1(t, structure);
    emitter.F2                   = gra::form::F2(t, structure);
    emitter.electric_form_factor = gra::form::G_E(Q2, structure);
    emitter.magnetic_form_factor = gra::form::G_M(Q2, structure);
    return emitter;
  }
  emitter.F1                   = 1.0;
  emitter.F2                   = 0.0;
  emitter.electric_form_factor = 1.0;
  emitter.magnetic_form_factor = 1.0;
  return emitter;
}

// Resolve one supported elastic spin-half QED emitter from a beam leg
ElasticEmitter Emitter(const gra::LORENTZSCALAR &lts, int leg, double t, const std::string &context) {
  if (lts.model_cache == nullptr) { throw std::logic_error(context + ": model cache is not set"); }
  return Emitter(Beam(lts, leg), t, context, lts.model_cache->Tune().Structure());
}

// Compute the exact elastic one-photon pair-helicity matrix
// M = (4 pi alpha Q1 Q2/t) J1_mu J2^mu
ElasticSpinHalfPhotonMatrix ElasticSpinHalfPhotonExchange(const gra::MDirac &dirac, const gra::MParticle &beam1,
                                                          const gra::MParticle &beam2, const gra::M4Vec &p1_in,
                                                          const gra::M4Vec &p2_in, const gra::M4Vec &p1_out,
                                                          const gra::M4Vec            &p2_out,
                                                          const gra::form::ParamStore &structure) {
  const double beam_scale = std::max({1.0, std::abs(p1_in.E()), std::abs(p2_in.E()), p1_in.P3mod(), p2_in.P3mod()});
  const double transverse_tolerance = 1.0e-10 * beam_scale;
  if (std::abs(p1_in.Px()) > transverse_tolerance || std::abs(p1_in.Py()) > transverse_tolerance ||
      std::abs(p2_in.Px()) > transverse_tolerance || std::abs(p2_in.Py()) > transverse_tolerance ||
      !(p1_in.Pz() > 0.0) || !(p2_in.Pz() < 0.0)) {
    throw std::invalid_argument("ElasticSpinHalfPhotonExchange: ordered collinear beams are required");
  }
  const M4Vec  q1    = p1_out - p1_in;
  const M4Vec  q2    = p2_out - p2_in;
  const double t1    = q1.M2();
  const double t2    = q2.M2();
  const double scale = std::max({1.0, std::abs(t1), std::abs(t2)});
  if (!std::isfinite(t1) || !std::isfinite(t2) || !(t1 < 0.0) || !(t2 < 0.0) ||
      std::abs(t1 - t2) > 1.0e-9 * scale) {
    throw std::invalid_argument("ElasticSpinHalfPhotonExchange: inconsistent spacelike transfers");
  }

  const ElasticEmitter                emitter1 = Emitter(beam1, t1, "ElasticSpinHalfPhotonExchange beam 1", structure);
  const ElasticEmitter                emitter2 = Emitter(beam2, t2, "ElasticSpinHalfPhotonExchange beam 2", structure);
  const double                        transfer_azimuth = q1.Phi();
  const auto current1 = dirac.DiracPauliCurrentBasis(p1_in, p1_out, q1, emitter1.fermion_kind,
                                                     emitter1.mass, emitter1.F1, emitter1.F2);
  const auto current2 = dirac.DiracPauliCurrentBasis(p2_in, p2_out, q2, emitter2.fermion_kind,
                                                     emitter2.mass, emitter2.F1, emitter2.F2);

  const double                coefficient = 4.0 * math::PI * alpha_0 * emitter1.charge * emitter2.charge / t1;
  ElasticSpinHalfPhotonMatrix amplitude{};
  for (std::size_t row = 0; row < 4; ++row) {
    const auto outgoing = gra::spin::BinaryPairHelicityLabelsX2(row);
    for (std::size_t col = 0; col < 4; ++col) {
      const auto                 incoming = gra::spin::BinaryPairHelicityLabelsX2(col);
      const std::size_t          out1     = gra::spin::BinaryHelicityIndexX2(outgoing[0]);
      const std::size_t          out2     = gra::spin::BinaryHelicityIndexX2(outgoing[1]);
      const std::size_t          in1      = gra::spin::BinaryHelicityIndexX2(incoming[0]);
      const std::size_t          in2      = gra::spin::BinaryHelicityIndexX2(incoming[1]);
      const std::complex<double> section_phase =
          gra::MDirac::ElasticSpinHalfColliderCurrentPhase(emitter1.fermion_kind, 1, incoming[0] / 2.0,
                                                           outgoing[0] / 2.0, transfer_azimuth) *
          gra::MDirac::ElasticSpinHalfColliderCurrentPhase(emitter2.fermion_kind, 2, incoming[1] / 2.0,
                                                           outgoing[1] / 2.0, transfer_azimuth);
      amplitude[gra::spin::PairHelicityMatrixIndex(row, col)] =
          section_phase * coefficient *
          dirac.g.BilinearForm(current1[gra::spin::BinaryPairHelicityIndex(out1, in1)],
                               current2[gra::spin::BinaryPairHelicityIndex(out2, in2)]);
    }
  }
  return amplitude;
}

// Validate one exclusive beam leg for the elastic QED current
void ValidateEmitter(const gra::LORENTZSCALAR &lts, int leg, const std::string &context) {
  const auto &particle = Beam(lts, leg);
  if (particle.spinX2 != 1) {
    throw std::invalid_argument(
        context + ": elastic QED current requires a spin-half beam, PDG = " + std::to_string(particle.pdg));
  }
  if (!(particle.mass > 0.0) || !std::isfinite(particle.mass)) {
    throw std::invalid_argument(context + ": invalid beam mass for PDG = " + std::to_string(particle.pdg));
  }

  if (particle.chargeX3 == 0) {
    throw std::invalid_argument(
        context + ": neutral beam cannot emit an elastic QED photon, PDG = " + std::to_string(particle.pdg));
  }
  const int absolute_pdg = std::abs(particle.pdg);
  if (absolute_pdg != PDG::PDG_p && absolute_pdg != 11 && absolute_pdg != 13 && absolute_pdg != 15) {
    throw std::invalid_argument(context + ": no elastic QED current model for beam PDG = " + std::to_string(particle.pdg));
  }
}

// Compute the electric charge of one photon-emitting beam leg
double Charge(const gra::LORENTZSCALAR &lts, int leg) { return Beam(lts, leg).chargeX3 / 3.0; }

// Compute the elastic EPA density used to match the QED current
double ElasticPhotonDensity(const gra::LORENTZSCALAR &lts, int leg) {
  const double         xi      = (leg == 1) ? lts.xi1 : (leg == 2) ? lts.xi2 : 0.0;
  const double         t       = (leg == 1) ? lts.t1 : (leg == 2) ? lts.t2 : 0.0;
  const double         qt      = (leg == 1) ? lts.qt1 : (leg == 2) ? lts.qt2 : 0.0;
  const ElasticEmitter emitter = Emitter(lts, leg, t, "gra::qed::ElasticPhotonDensity");
  if (!(xi > 0.0)) { return 0.0; }
  return gra::flux::ElasticSpinHalfFluxTransverse(xi, t, qt, emitter.mass, emitter.electric_form_factor,
                                                  emitter.magnetic_form_factor, emitter.charge)
             .Trace() /
         xi;
}

// Compute the EPA or matched QED photon-source amplitude normalization
// A_source = sqrt(n_gamma/xi), or the density-matched QED current
std::complex<double> PhotonSourceAmplitude(const gra::LORENTZSCALAR &lts, const gra::ForwardLegState &state,
                                           const std::string &mode) {
  const double charge = Charge(lts, state.Index());
  if (mode == "EPA") {
    const double flux = gra::flux::ForwardPhotonFlux(state);
    if (!(flux > 0.0) || !(state.xi > 0.0)) { return 0.0; }
    const double sign = charge < 0.0 ? -1.0 : 1.0;
    return sign * gra::math::msqrt(flux / state.xi);
  }
  if (mode == "QED" && lts.process.SPINGEN && !state.IsExcited()) { return 1.0; }

  if (!state.IsExcited()) {
    const double density = ElasticPhotonDensity(lts, state.Index());
    const double sign    = (charge < 0.0) ? -1.0 : 1.0;
    return sign * gra::math::msqrt(density);
  }

  const double inclusive_flux = gra::flux::ForwardPhotonFlux(state);
  if (!(inclusive_flux > 0.0) || !(state.xi > 0.0)) { return 0.0; }
  return charge * gra::math::msqrt(inclusive_flux / state.xi);
}

// Compute the transverse EPA source amplitude for one global spin projection
// A_m = sqrt(n_gamma/xi) exp(i m phi)/sqrt(2), including charge sign
std::complex<double> TransversePhotonSourceAmplitude(const LORENTZSCALAR &lts, const ForwardLegState &state,
                                                     const int boson_spin_projection) {
  if (boson_spin_projection != -1 && boson_spin_projection != 1) {
    throw std::invalid_argument("qed::TransversePhotonSourceAmplitude: spin projection must be +/-1");
  }
  const std::complex<double> scalar = PhotonSourceAmplitude(lts, state, "EPA");
  if (std::fpclassify(std::abs(scalar)) == FP_ZERO) { return 0.0; }

  // The lower local photon label is minus the shared produced-boson projection
  const std::complex<double> phase =
      std::exp(math::zi * static_cast<double>(boson_spin_projection) * state.transfer.Phi());
  return scalar * phase / std::sqrt(2.0);
}

// Build a matched EPA or QED photon source for explicit helicity transitions
MMatrix<std::complex<double>> PhotonSourceMatrixTransitions(const gra::LORENTZSCALAR &lts, int leg,
                                                            const std::vector<std::pair<double, double>> &transitions,
                                                            const std::vector<int> &m_values, const std::string &mode,
                                                            bool second_exchange_daughter) {
  if (mode == "EPA" || !ExclusivePhotonEmitter(lts, leg)) {
    return PhotonEPAMatrix(lts, leg, transitions, m_values, second_exchange_daughter);
  }

  const auto                          full_transitions = gra::spin::SpinHalfTransitions(false);
  const MMatrix<std::complex<double>> raw_full         = RawQEDPhotonMatrix(lts, leg, full_transitions, m_values);
  MMatrix<std::complex<double>>       phased_full      = raw_full;
  ApplyPhotonSourceSectionPhase(phased_full, m_values, ((leg == 1) ? lts.q1 : lts.q2).Phi(), second_exchange_daughter);
  const double         t       = (leg == 1) ? lts.t1 : lts.t2;
  const ElasticEmitter emitter = Emitter(lts, leg, t, "gra::qed::PhotonSourceMatrixTransitions");
  const double         target  = ElasticPhotonDensity(lts, leg);
  const double         rho_qed = gra::spin::SourceSpinAveragedDensity(
              raw_full, emitter.incoming_spin_states, "gra::qed::PhotonSourceMatrixTransitions raw QED current");
  if (target <= 0.0 || rho_qed <= 0.0) {
    return MMatrix<std::complex<double>>(transitions.size(), m_values.size(), 0.0);
  }
  phased_full *= gra::math::msqrt(target / rho_qed);

  MMatrix<std::complex<double>> reference =
      PhotonEPAMatrix(lts, leg, full_transitions, m_values, second_exchange_daughter);
  const double charge_sign = (emitter.charge < 0.0) ? -1.0 : 1.0;
  reference *= charge_sign * gra::math::msqrt(target);
  AlignPhotonSourceRows(phased_full, reference);
  return SelectPhotonTransitionRows(phased_full, full_transitions, transitions);
}

}  // namespace qed
}  // namespace gra
