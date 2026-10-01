// XP pole amplitudes and covariant numerator projections
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Regge/MReggeXP.h"

#include "Graniitti/Spin/MDirac.h"
#include "Graniitti/Spin/MWigner.h"
#include "Graniitti/Tech/MException.h"

namespace gra::xpom {

// Compute exchange-helicity rows by resonance-spin columns in the cached X rest frame
HelAmp Fusion(const LORENTZSCALAR& lts, const spin::PoleLS& vertex) {
  const auto& q1_in_X = lts.q1_in_X;
  const std::size_t rows = vertex.helicity.lambda_values.size_row();
  const std::size_t cols = vertex.helicity.Jz_values.size();

  const auto reduced = spin::PoleLSReduced(vertex, q1_in_X.P3mod());
  const double phi = q1_in_X.Phi();

  HelAmp     central(rows, cols, 0.0);
  const auto rotation = vertex.rotation->Evaluate(q1_in_X.Theta());
  for (std::size_t row = 0; row < rows; ++row) {
    const std::size_t i1 = vertex.helicity.lambda_idx[row][0];
    const std::size_t i2 = vertex.helicity.lambda_idx[row][1];

    const std::size_t rotation_row = vertex.helicity.jw_rotation.row[row];
    const double lambda = vertex.helicity.jw_rotation.difference[rotation_row];
    for (std::size_t col = 0; col < cols; ++col) {
      // C_{lambda,M} = G_M d^J_{lambda,-M}(theta) exp[-i(M+lambda)phi] H_lambda
      const auto phase = std::exp(std::complex<double>(0.0, (-vertex.helicity.Jz_values[col] - lambda) * phi));
      central[row][col] = vertex.spin_metric[col] * rotation(rotation_row, col) * phase * reduced[i1][i2];
    }
  }
  return central;
}

// Build the XP resonance production matrices
std::vector<HelAmp> Resonance(const LORENTZSCALAR& lts, const PARAM_RES& res, const double s0,
                              const std::array<ForwardLegState, 2>* forward_state) {
  return rspin::Resonance(lts, res, Fusion, s0, nullptr, forward_state);
}

// Build the XP crossed continuum production matrices
std::vector<HelPair> Continuum(const LORENTZSCALAR& lts, const double s0) {
  // Reduced LS vertices sew the physical pole coefficient with no internal spin average
  return rspin::Continuum(lts, s0);
}

namespace {

// Project (slash k + m) using ubar_h u_h' = 2m delta_hh' and ubar_h gamma^0 u_h' = 2E_pole delta_hh'
// The pole value is 4m^2, fixing the same unit pole coefficient used by the DL couplings
HelAmp DiracProjection(const M4Vec& momentum, double mass) {
  const double p           = std::hypot(momentum.Px(), momentum.Py(), momentum.Pz());
  const double pole_energy = std::hypot(p, mass);
  const double value       = std::fma(0.5 * (pole_energy / mass), (momentum.E() - pole_energy) / mass, 1.0);
  if (!std::isfinite(value)) { throw AmplitudeFailure("MReggeXP Dirac: non-finite projected numerator"); }
  return HelAmp::IdentityMatrix(2) * value;
}

// Evaluate the Proca projection in the physical pole basis without boost cancellations
// epsilon_0.k/m = |p|(E-E_pole)/m^2 and -epsilon_h*.epsilon_h' = delta_hh'
HelAmp ProcaProjection(const M4Vec& momentum, double mass) {
  const double p           = std::hypot(momentum.Px(), momentum.Py(), momentum.Pz());
  const double pole_energy = std::hypot(p, mass);
  const double overlap     = (p / mass) * ((momentum.E() - pole_energy) / mass);
  auto         output      = HelAmp::IdentityMatrix(3);
  output[1][1] += overlap * overlap;
  return output;
}

}  // namespace

// Classify one supported scalar, fermion, vector or photon pole field
PoleType ClassifyPole(const MParticle& particle) {
  if (particle.pdg == PDG::PDG_gamma) { return PoleType::Photon; }
  if (particle.spinX2 == 0) { return PoleType::Scalar; }
  if (particle.spinX2 == 1 || particle.spinX2 == 2) {
    if (!std::isfinite(particle.mass) || !(particle.mass > 0.0)) {
      throw std::invalid_argument("MReggeXP::ClassifyPole: massive pole requires a positive finite mass");
    }
    return particle.spinX2 == 1 ? PoleType::Dirac : PoleType::Proca;
  }
  throw std::invalid_argument("MReggeXP::ClassifyPole: unsupported XP field");
}

// Compute the unit pole coefficient supported by reduced continuum vertices
PoleNumerator ReducedNumerator(const MParticle& particle, const M4Vec& momentum) {
  if (!kinematics::FiniteFourVector(momentum)) {
    throw AmplitudeFailure("MReggeXP::ReducedNumerator: non-finite momentum");
  }
  // Off-shell Dirac and longitudinal Proca terms require covariant vertex currents
  // The reduced LS vertices specify only the unit pole coefficient, as in MP and GP
  PoleNumerator output;
  output.type       = ClassifyPole(particle);
  output.helicities = spin::FinalStateHelicities(particle, "MReggeXP::ReducedNumerator");
  output.matrix     = HelAmp::IdentityMatrix(output.helicities.size());
  return output;
}

// Project a covariant numerator into the positive-energy pole basis for diagnostics
PoleNumerator Numerator(const MParticle& particle, const M4Vec& momentum) {
  auto output = ReducedNumerator(particle, momentum);
  if (output.type == PoleType::Dirac) {
    output.matrix = DiracProjection(momentum, particle.mass);
  } else if (output.type == PoleType::Proca) {
    output.matrix = ProcaProjection(momentum, particle.mass);
  }
  return output;
}

// Compute the antisymmetric photon field-strength tensor k^mu eps^nu-k^nu eps^mu
HelAmp FieldStrength(const M4Vec& momentum, int helicity) {
  if (!kinematics::FiniteFourVector(momentum)) {
    throw AmplitudeFailure("MReggeXP::FieldStrength: non-finite momentum");
  }
  const MDirac dirac("DIRAC");
  const auto   polarization = dirac.EpsSpin1(momentum, helicity);
  HelAmp       field(4, 4, 0.0);
  for (std::size_t mu = 0; mu < 4; ++mu) {
    for (std::size_t nu = 0; nu < 4; ++nu) {
      field[mu][nu] = momentum[mu] * polarization(nu) - momentum[nu] * polarization(mu);
    }
  }
  return field;
}

// Express the reduced XP pole coefficient in the shared continuum helicity metric
spin::InternalHelicityMetric PoleMetric(const MParticle& particle, const M4Vec& momentum) {
  auto numerator = ReducedNumerator(particle, momentum);
  return {std::move(numerator.helicities), std::move(numerator.matrix)};
}

}  // namespace gra::xpom
