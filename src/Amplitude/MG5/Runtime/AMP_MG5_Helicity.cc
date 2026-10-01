// Shared MG5 photon source helicity contraction
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Amplitude/MG5/Runtime/AMP_MG5_Helicity.h"

#include <cstdint>
#include <limits>

#include "Graniitti/Photon/MFlux.h"
#include "Graniitti/Photon/MQED.h"
#include "Graniitti/Spin/MHelicityScatter.h"
#include "Graniitti/Tech/MAux.h"

using gra::aux::indices;

namespace gra {
namespace mg5helas {

namespace {

// Select the incoming source basis for one EPA hard contraction
enum class EPASourceMode { Projected, Factorized };

// Compute one resolved forward photon leg
ForwardLegState PhotonLegState(const LORENTZSCALAR &lts, const int leg) {
  if (leg == 1) { return ResolveForwardLegState(lts, ForwardBeamLeg::Upper); }
  if (leg == 2) { return ResolveForwardLegState(lts, ForwardBeamLeg::Lower); }
  throw std::invalid_argument("mg5helas::PhotonLegState: photon leg should be 1 or 2");
}

// Rotate one normalized real linear-polarization source by 90 degrees
// e_perp = (-e_+^*,e_-^*) from e_parallel = (e_-,e_+)
std::array<std::complex<double>, 2> PerpendicularHelicitySource(const std::array<std::complex<double>, 2> &parallel) {
  return {-std::conj(parallel[1]), std::conj(parallel[0])};
}

// Compute the normalized parallel EPA source in the active MG5 helicity chart
std::array<std::complex<double>, 2> ParallelEPASource(const LORENTZSCALAR &lts, int leg, const M4Vec &source1,
                                                      const M4Vec &source2, const M4Vec &k1, const M4Vec &k2,
                                                      EPASourceMode mode) {
  if (mode == EPASourceMode::Factorized) {
    return leg == 1 ? TransverseHelicityCoefficients(source1, k1, k2, true)
                    : TransverseHelicityCoefficients(source2, k2, k1, true);
  }
  return leg == 1 ? TransverseHelicityCoefficients(lts.pbeam1, k1, k2, true)
                  : TransverseHelicityCoefficients(lts.pbeam2, k2, k1, true);
}

// Build transverse EPA eigen-sources including the MG5 photon spin average
// source_a = sqrt(2) a_a e_a / sqrt[Tr(rho)], rho_a = a_a^2
MMatrix<std::complex<double>> PhotonDensityEigenSources(const LORENTZSCALAR &lts, const int leg,
                                                        const ForwardLegState                     &state,
                                                        const std::vector<flux::PhotonFluxSector> &sector,
                                                        const M4Vec &k1, const M4Vec &k2, const M4Vec &source1,
                                                        const M4Vec &source2, const EPASourceMode mode) {
  flux::TransversePhotonFlux total;
  for (const auto &item : sector) {
    total.parallel += item.density.parallel;
    total.perpendicular += item.density.perpendicular;
  }
  const double                  trace = total.Trace();
  MMatrix<std::complex<double>> sources(2 * sector.size(), 2, 0.0);
  if (!(trace > 0.0) || !std::isfinite(trace)) { return sources; }

  const auto   parallel      = ParallelEPASource(lts, leg, source1, source2, k1, k2, mode);
  const auto   perpendicular = PerpendicularHelicitySource(parallel);
  const double charge_phase  = qed::Charge(lts, leg) < 0.0 ? -1.0 : 1.0;
  const double inverse_norm  = std::sqrt(2.0 / trace);
  for (const auto &s : indices(sector)) {
    const double parallel_scale      = inverse_norm * sector[s].source.parallel;
    const double perpendicular_scale = inverse_norm * sector[s].source.perpendicular;
    for (std::size_t h = 0; h < 2; ++h) {
      sources[2 * s][h]     = charge_phase * parallel_scale * parallel[h];
      sources[2 * s + 1][h] = charge_phase * perpendicular_scale * perpendicular[h];
    }
  }
  return sources;
}

// Contract one hard tensor with resolved EPA sources in the selected basis
// M_ab = sum_(h1,h2) e1_h1 e2_h2 H_h1h2
std::vector<std::complex<double>> ContractPhotonSources(LORENTZSCALAR                                  &lts,
                                                        const std::vector<mg5helas::HelicityComponent> &components,
                                                        const M4Vec &source1, const M4Vec &source2, const M4Vec &k1,
                                                        const M4Vec &k2, EPASourceMode mode) {
  struct HardGroup {
    std::vector<int>                    outgoing;
    std::size_t                         color = 0;
    std::array<std::complex<double>, 4> hard  = {};
  };

  std::vector<HardGroup> groups;
  for (const auto &component : components) {
    if (std::abs(component.incoming[0]) != 1 || std::abs(component.incoming[1]) != 1) {
      throw std::invalid_argument("mg5helas::ContractPhotonSources: photon helicity is not transverse");
    }
    auto found = std::find_if(groups.begin(), groups.end(), [&](const HardGroup &group) {
      return group.color == component.color && group.outgoing == component.outgoing;
    });
    if (found == groups.end()) {
      groups.push_back({component.outgoing, component.color, {}});
      found = std::prev(groups.end());
    }
    const std::size_t pair = spin::BinaryPairHelicityIndexX2(component.incoming[0], component.incoming[1]);
    found->hard[pair] += component.value;
  }

  const std::array<ForwardLegState, 2>                     state  = {PhotonLegState(lts, 1), PhotonLegState(lts, 2)};
  const std::array<std::vector<flux::PhotonFluxSector>, 2> sector = {flux::ForwardPhotonFluxSectors(state[0]),
                                                                     flux::ForwardPhotonFluxSectors(state[1])};
  const MMatrix<std::complex<double>>                      upper =
      PhotonDensityEigenSources(lts, 1, state[0], sector[0], k1, k2, source1, source2, mode);
  const MMatrix<std::complex<double>> lower =
      PhotonDensityEigenSources(lts, 2, state[1], sector[1], k1, k2, source1, source2, mode);
  lts.hamp.layout.nuclear_type        = nuclear::ScreenType::Fusion;
  lts.hamp.layout.epa_sector_resolved = true;
  lts.hamp.layout.epa_rows_per_sector = 2;
  for (const auto &leg : indices(state)) {
    lts.hamp.layout.epa_sector_count[leg] = static_cast<std::uint8_t>(sector[leg].size());
    lts.hamp.layout.epa_sector_type[leg]  = {255, 255};
    for (const auto &i : indices(sector[leg])) { lts.hamp.layout.epa_sector_type[leg][i] = sector[leg][i].code; }
  }

  std::vector<std::complex<double>> out;
  out.reserve(groups.size() * upper.size_row() * lower.size_row());
  for (const auto &group : groups) {
    for (std::size_t row1 = 0; row1 < upper.size_row(); ++row1) {
      for (std::size_t row2 = 0; row2 < lower.size_row(); ++row2) {
        out.push_back(gra::KroneckerBilinearProduct(upper.Row(row1), lower.Row(row2), group.hard));
      }
    }
  }
  return out;
}

// Validate the stored hard invariant against the summed final state
bool EPAHardInvariant(const LORENTZSCALAR &lts, const M4Vec &hard, long double &mass2) {
  const auto        h          = hard.Contravariant<long double>();
  const long double calculated = h[0] * h[0] - h[1] * h[1] - h[2] * h[2] - h[3] * h[3];
  const long double stored     = lts.s_hat;
  const long double scale      = 1.0L + h[0] * h[0] + h[1] * h[1] + h[2] * h[2] + h[3] * h[3] + std::abs(stored);
  const long double tolerance  = 64.0L * std::numeric_limits<double>::epsilon() * scale;
  if (!std::isfinite(calculated) || !std::isfinite(stored) || !(stored > 0.0L) ||
      std::abs(calculated - stored) > tolerance || !std::isfinite(hard.E()) || !(hard.E() > 0.0)) {
    mass2 = 0.0L;
    return false;
  }
  mass2 = stored;
  return true;
}

// Construct the minimal rotation of the hard upper-photon axis onto +z
bool EPAHardAxisRotation(const M4Vec &direction, MMatrix<double> &rotation) {
  const auto        p     = direction.Contravariant<long double>();
  const long double norm2 = p[1] * p[1] + p[2] * p[2] + p[3] * p[3];
  if (!std::isfinite(norm2) || !(norm2 > 0.0L)) { return false; }
  const long double inverse_norm = 1.0L / std::sqrt(norm2);
  const long double nx           = p[1] * inverse_norm;
  const long double ny           = p[2] * inverse_norm;
  const long double nz           = std::clamp(p[3] * inverse_norm, -1.0L, 1.0L);

  // The antiparallel chart needs one deterministic transverse axis
  const long double denominator = 1.0L + nz;
  if (denominator <= 64.0L * std::numeric_limits<long double>::epsilon()) {
    rotation = direction.RotationTo(M3Vec{0.0, 0.0, 1.0});
    return true;
  }

  // Rodrigues rotation with v = n x z remains stable as n approaches +z
  const long double vx  = ny;
  const long double vy  = -nx;
  const long double inv = 1.0L / denominator;
  rotation              = {
                   {static_cast<double>(1.0L - vy * vy * inv), static_cast<double>(vx * vy * inv), static_cast<double>(vy)},
                   {static_cast<double>(vx * vy * inv), static_cast<double>(1.0L - vx * vx * inv), static_cast<double>(-vx)},
                   {static_cast<double>(-vy), static_cast<double>(vx), static_cast<double>(1.0L - (vx * vx + vy * vy) * inv)}};
  for (std::size_t row = 0; row < rotation.size_row(); ++row) {
    for (std::size_t column = 0; column < rotation.size_col(); ++column) {
      if (!std::isfinite(rotation[row][column])) { return false; }
    }
  }
  return true;
}

}  // namespace

// Prepare covariant massless incoming momenta while preserving the final state
bool PrepareOnShellKinematics(const LORENTZSCALAR &lts, const std::vector<M4Vec> &final, M4Vec &p1, M4Vec &p2) {
  return kinematics::CovariantOnShellInitialState(lts.q1, lts.q2, final, p1, p2);
}

namespace {

// Prepare the transfer-independent leading-power kT-EPA hard incoming pair
bool PrepareEPAHardKinematics(const LORENTZSCALAR &lts, const std::vector<M4Vec> &final, M4Vec &p1, M4Vec &p2) {
  if (final.empty()) { return false; }
  M4Vec hard;
  for (const auto &momentum : final) { hard += momentum; }
  long double hard_mass2 = 0.0L;
  if (!EPAHardInvariant(lts, hard, hard_mass2)) { return false; }
  const double upper_product = lts.pbeam1 * hard;
  const double lower_product = lts.pbeam2 * hard;
  if (!std::isfinite(upper_product) || !std::isfinite(lower_product) || !(upper_product > 0.0) ||
      !(lower_product > 0.0)) {
    return false;
  }
  // Use the symmetric beam direction in the hard rest frame
  // [REFERENCE: Collins and Soper, Phys. Rev. D16 (1977) 2219]
  const M4Vec axis       = (lts.pbeam1 / upper_product - lts.pbeam2 / lower_product) * static_cast<double>(hard_mass2);
  const M4Vec reference1 = hard * 0.5 + axis;
  const M4Vec reference2 = hard - reference1;
  return kinematics::CovariantOnShellInitialState(reference1, reference2, final, p1, p2);
}

}  // namespace

// Prepare one reusable EPA hard rest frame and transform all final momenta
bool PrepareEPAHardFrame(const LORENTZSCALAR &lts, std::vector<M4Vec> &final, EPAHardFrame &frame) {
  frame = {};
  if (!PrepareEPAHardKinematics(lts, final, frame.incoming[0], frame.incoming[1])) { return false; }
  frame.beam = {lts.pbeam1, lts.pbeam2};

  for (const auto &momentum : final) { frame.hard += momentum; }
  long double mass2 = 0.0L;
  if (!EPAHardInvariant(lts, frame.hard, mass2)) { return false; }

  frame.mass = static_cast<double>(std::sqrt(mass2));
  for (auto &momentum : frame.incoming) { kinematics::LorentzBoost(frame.hard, frame.mass, momentum, -1); }
  for (auto &momentum : final) { kinematics::LorentzBoost(frame.hard, frame.mass, momentum, -1); }

  const std::array<M4Vec, 2> direction = frame.incoming;
  if (!kinematics::CovariantOnShellInitialState(direction[0], direction[1], final, frame.incoming[0],
                                                frame.incoming[1])) {
    return false;
  }

  // Align the physical incoming axis with the stable planar beam chart
  if (!EPAHardAxisRotation(frame.incoming[0], frame.rotation)) { return false; }
  for (auto &momentum : final) { momentum.Rotate(frame.rotation); }

  // HELAS selects the collinear polarization chart from the exact pt = 0 case
  const double photon_energy = 0.5 * frame.mass;
  frame.incoming = {M4Vec(0.0, 0.0, photon_energy, photon_energy), M4Vec(0.0, 0.0, -photon_energy, photon_energy)};

  const auto finite = [](const M4Vec &momentum) {
    return std::isfinite(momentum.E()) && std::isfinite(momentum.Px()) && std::isfinite(momentum.Py()) &&
           std::isfinite(momentum.Pz());
  };
  frame.valid = frame.incoming[0].E() > 0.0 && frame.incoming[1].E() > 0.0 &&
                std::all_of(frame.incoming.cbegin(), frame.incoming.cend(), finite) &&
                std::all_of(final.cbegin(), final.cend(),
                            [&](const M4Vec &momentum) { return momentum.E() > 0.0 && finite(momentum); });
  return frame.valid;
}

// Remove the components of one EPA source in the physical two-beam plane
bool BeamTransverseSource(const M4Vec &q, const std::array<M4Vec, 2> &beam, M4Vec &source) {
  // q_T = q - a p_1 - b p_2 with (a,b) fixed by q_T.p_1 = q_T.p_2 = 0
  const auto dot = [](const M4Vec &a, const M4Vec &b) {
    const auto av = a.Contravariant<long double>();
    const auto bv = b.Contravariant<long double>();
    return av[0] * bv[0] - av[1] * bv[1] - av[2] * bv[2] - av[3] * bv[3];
  };
  const long double beam11 = dot(beam[0], beam[0]);
  const long double beam12 = dot(beam[0], beam[1]);
  const long double beam22 = dot(beam[1], beam[1]);
  const long double qbeam1 = dot(q, beam[0]);
  const long double qbeam2 = dot(q, beam[1]);
  const long double gram   = beam12 * beam12 - beam11 * beam22;
  const long double scale  = std::max({1.0L, std::abs(beam12 * beam12), std::abs(beam11 * beam22)});
  if (!std::isfinite(gram) || !(gram > 64.0L * std::numeric_limits<long double>::epsilon() * scale)) {
    source = {};
    return false;
  }

  // This is the exact inverse of the massive two-beam Gram matrix
  const long double upper  = (beam12 * qbeam2 - beam22 * qbeam1) / gram;
  const long double lower  = (beam12 * qbeam1 - beam11 * qbeam2) / gram;
  const auto        qv     = q.Contravariant<long double>();
  const auto        p1     = beam[0].Contravariant<long double>();
  const auto        p2     = beam[1].Contravariant<long double>();
  const long double energy = qv[0] - upper * p1[0] - lower * p2[0];
  const long double px     = qv[1] - upper * p1[1] - lower * p2[1];
  const long double py     = qv[2] - upper * p1[2] - lower * p2[2];
  const long double pz     = qv[3] - upper * p1[3] - lower * p2[3];
  source =
      M4Vec(static_cast<double>(px), static_cast<double>(py), static_cast<double>(pz), static_cast<double>(energy));
  return std::isfinite(source.E()) && std::isfinite(source.Px()) && std::isfinite(source.Py()) &&
         std::isfinite(source.Pz());
}

// Transform one node-local beam-transverse source pair into a prepared hard frame
bool TransformEPAHardSources(const EPAHardFrame &frame, const M4Vec &q1, const M4Vec &q2,
                             std::array<M4Vec, 2> &source) {
  if (!frame.valid) { return false; }
  if (!BeamTransverseSource(q1, frame.beam, source[0]) || !BeamTransverseSource(q2, frame.beam, source[1])) {
    return false;
  }
  for (auto &momentum : source) {
    kinematics::LorentzBoost(frame.hard, frame.mass, momentum, -1);
    momentum.Rotate(frame.rotation);
  }
  return std::all_of(source.cbegin(), source.cend(), [](const M4Vec &momentum) {
    return std::isfinite(momentum.E()) && std::isfinite(momentum.Px()) && std::isfinite(momentum.Py()) &&
           std::isfinite(momentum.Pz());
  });
}

// Remove longitudinal components along a physical on-shell incoming pair
// q_perp = q-k1(q.k2)/(k1.k2)-k2(q.k1)/(k1.k2)
M4Vec TransverseSourceVector(const M4Vec &source, const M4Vec &k1, const M4Vec &k2) {
  const double pair_product = k1 * k2;
  if (!std::isfinite(pair_product) || pair_product <= 0.0) { return M4Vec(); }
  return source - k1 * ((source * k2) / pair_product) - k2 * ((source * k1) / pair_product);
}

// Compute the HELAS massless incoming-vector polarization in (t,x,y,z) order
std::array<std::complex<double>, 4> IncomingVectorPolarization(const M4Vec &momentum, int helicity) {
  const double spatial_norm = momentum.P3mod();
  if (!std::isfinite(spatial_norm) || spatial_norm <= 0.0 || std::abs(helicity) != 1) { return {}; }

  const double inverse_sqrt_two = 1.0 / std::sqrt(2.0);
  const double transverse_norm  = momentum.Pt();
  const double cos_theta        = momentum.Pz() / spatial_norm;
  const double sin_theta        = transverse_norm / spatial_norm;
  const double cos_phi          = transverse_norm > 0.0 ? momentum.Px() / transverse_norm : 1.0;
  const double sin_phi          = transverse_norm > 0.0 ? momentum.Py() / transverse_norm : 0.0;
  const double h                = static_cast<double>(helicity);

  return {0.0, inverse_sqrt_two * (-h * cos_theta * cos_phi + math::zi * sin_phi),
          inverse_sqrt_two * (-h * cos_theta * sin_phi - math::zi * cos_phi), inverse_sqrt_two * h * sin_theta};
}

// Project one real source vector onto the physical HELAS helicities
// c_lambda = -q_perp.epsilon_lambda^*/sqrt(-q_perp^2)
std::array<std::complex<double>, 2> TransverseHelicityCoefficients(const M4Vec &source, const M4Vec &k1,
                                                                   const M4Vec &k2, bool normalize) {
  const M4Vec  transverse = TransverseSourceVector(source, k1, k2);
  const double norm2      = -transverse.M2();
  if (!std::isfinite(norm2) || norm2 <= 1.0e-24) { return {}; }

  const double                        scale        = normalize ? 1.0 / std::sqrt(norm2) : 1.0;
  std::array<std::complex<double>, 2> coefficients = {};
  for (const auto &index : indices(coefficients)) {
    const int                  helicity = spin::BinaryHelicityLabelX2(index);
    const auto                 epsilon  = IncomingVectorPolarization(k1, helicity);
    const std::complex<double> contraction =
        transverse.E() * std::conj(epsilon[0]) - transverse.Px() * std::conj(epsilon[1]) -
        transverse.Py() * std::conj(epsilon[2]) - transverse.Pz() * std::conj(epsilon[3]);
    coefficients[index] = -scale * contraction;
  }
  return coefficients;
}

// Prepare the linear transverse Cartesian map into the physical helicity basis
TransverseHelicityProjector PrepareTransverseHelicityProjector(const M4Vec &k1, const M4Vec &k2) {
  TransverseHelicityProjector projector;
  const double                pair_product = k1 * k2;
  if (!std::isfinite(pair_product) || pair_product <= 0.0) { return projector; }

  const M4Vec basis_x = TransverseSourceVector(M4Vec(1.0, 0.0, 0.0, 0.0), k1, k2);
  const M4Vec basis_y = TransverseSourceVector(M4Vec(0.0, 1.0, 0.0, 0.0), k1, k2);
  projector.norm_xx   = -basis_x.M2();
  projector.norm_xy   = -(basis_x * basis_y);
  projector.norm_yy   = -basis_y.M2();

  for (std::size_t index = 0; index < 2; ++index) {
    const int  helicity = spin::BinaryHelicityLabelX2(index);
    const auto epsilon  = IncomingVectorPolarization(k1, helicity);
    const auto project  = [&](const M4Vec &basis) {
      const std::complex<double> contraction = basis.E() * std::conj(epsilon[0]) - basis.Px() * std::conj(epsilon[1]) -
                                               basis.Py() * std::conj(epsilon[2]) - basis.Pz() * std::conj(epsilon[3]);
      return -contraction;
    };
    projector.coefficient[0][index] = project(basis_x);
    projector.coefficient[1][index] = project(basis_y);
  }
  projector.valid =
      std::isfinite(projector.norm_xx) && std::isfinite(projector.norm_xy) && std::isfinite(projector.norm_yy);
  return projector;
}

// Project one transverse Cartesian vector through a prepared covariant basis
// c_lambda = qx c_x,lambda+qy c_y,lambda
std::array<std::complex<double>, 2> TransverseHelicityProjector::Project(double qx, double qy) const {
  if (!Accepts(qx, qy)) { return {}; }
  return {qx * coefficient[0][0] + qy * coefficient[1][0], qx * coefficient[0][1] + qy * coefficient[1][1]};
}

// Test the prepared transverse norm against the source-projection threshold
bool TransverseHelicityProjector::Accepts(double qx, double qy) const {
  if (!valid || !std::isfinite(qx) || !std::isfinite(qy)) { return false; }
  const double norm2 = qx * qx * norm_xx + 2.0 * qx * qy * norm_xy + qy * qy * norm_yy;
  return std::isfinite(norm2) && norm2 > 1.0e-24;
}

// Contract two prepared source-helicity vectors with one hard amplitude
// output = sum_(h1,h2) c1_h1 c2_h2 H_h1h2
std::complex<double> ContractHelicitySources(const std::array<std::complex<double>, 4> &hard,
                                             const std::array<std::complex<double>, 2> &source1,
                                             const std::array<std::complex<double>, 2> &source2) {
  return gra::KroneckerBilinearProduct(source1, source2, hard);
}

// Contract one integrated incoming-helicity kernel with one hard amplitude
// output = sum_(h1,h2) K_h1h2 H_h1h2
std::complex<double> ContractHelicityKernel(const std::array<std::complex<double>, 4> &hard,
                                            const std::array<std::complex<double>, 4> &kernel) {
  return gra::BilinearProduct(kernel, hard);
}

// Contract a Durham transverse-source pair with hard amplitudes ordered as
// (--,-+,+-,++)
std::complex<double> ContractTransverseSources(const std::array<std::complex<double>, 4> &hard, const M4Vec &source1,
                                               const M4Vec &source2, const M4Vec &k1, const M4Vec &k2) {
  const auto upper = TransverseHelicityCoefficients(source1, k1, k2, false);
  const auto lower = TransverseHelicityCoefficients(source2, k2, k1, false);
  return ContractHelicitySources(hard, upper, lower);
}

// Contract fixed-basis gamma-gamma hard amplitudes with resolved EPA sources
std::vector<std::complex<double>> ContractEPAPhotonSources(LORENTZSCALAR                        &lts,
                                                           const std::vector<HelicityComponent> &components,
                                                           const M4Vec &k1, const M4Vec &k2) {
  return ContractPhotonSources(lts, components, lts.q1, lts.q2, k1, k2, EPASourceMode::Projected);
}

// Contract a transfer-independent hard tensor with node-local transverse sources
std::vector<std::complex<double>> ContractEPAHardPhotonSources(LORENTZSCALAR                        &lts,
                                                               const std::vector<HelicityComponent> &components,
                                                               const M4Vec &source1, const M4Vec &source2,
                                                               const M4Vec &k1, const M4Vec &k2) {
  return ContractPhotonSources(lts, components, source1, source2, k1, k2, EPASourceMode::Factorized);
}

// Compute whether a fixed-hard EPA point has two exclusive proton spin sources
bool UsesProtonEPAHardSources(const LORENTZSCALAR &lts) {
  return lts.upc_model == nullptr && lts.upc_event == nullptr && !lts.excite1 && !lts.excite2 &&
         std::abs(lts.beam1.pdg) == PDG::PDG_p && std::abs(lts.beam2.pdg) == PDG::PDG_p && lts.beam1.spinX2 == 1 &&
         lts.beam2.spinX2 == 1;
}

namespace {

// Project all elastic Dirac-Pauli transitions into one fixed hard frame
MMatrix<std::complex<double>> ProtonEPASources(const LORENTZSCALAR &lts, const int leg, const EPAHardFrame &frame) {
  const auto                    transitions = spin::SpinHalfTransitions(false);
  auto                          currents    = qed::PhotonCurrentTransitions(lts, leg, transitions);
  MMatrix<std::complex<double>> source(transitions.size(), 2, 0.0);
  const M4Vec                  &photon = frame.incoming[static_cast<std::size_t>(leg - 1)];
  const M4Vec                  &q = leg == 1 ? lts.q1 : lts.q2;
  const M4Vec                  &n = leg == 1 ? lts.pbeam2 : lts.pbeam1;
  const std::array<double, 4>   q_lo = {q % 0, q % 1, q % 2, q % 3};
  for (const auto &row : indices(currents)) {
    // Use the Ward identity before projecting onto the on-shell hard photons
    // J -> J - q (J.n)/(q.n) retains each emitting leg's transverse direction
    AddScaled(currents[row], q_lo, -BilinearProduct(currents[row], n.Contravariant<double>()) / (q * n));
    // Keep collider helicity labels fixed while transforming the four-current
    currents[row] = MDirac::BoostCurrent(currents[row], frame.hard, frame.mass, -1);
    currents[row] = MDirac::RotateCurrent(currents[row], frame.rotation);
    for (std::size_t column = 0; column < 2; ++column) {
      const int  helicity = spin::BinaryHelicityLabelX2(column);
      const auto epsilon  = IncomingVectorPolarization(photon, helicity);
      for (std::size_t mu = 0; mu < 4; ++mu) { source[row][column] += currents[row][mu] * std::conj(epsilon[mu]); }
    }
  }

  // Normalize only the positive spin trace so the outer EPA flux enters once
  // (1/2) sum_(lambda_i,lambda_f,lambda_gamma) |J_hat|^2 = 1
  const double density = spin::SourceSpinAveragedDensity(source, 2, "mg5helas::ProtonEPASources");
  if (!(density > 0.0) || !std::isfinite(density)) { return MMatrix<std::complex<double>>(transitions.size(), 2, 0.0); }
  source *= 1.0 / std::sqrt(density);
  return source;
}

// Contract two proton transition banks in canonical initial-final pair order
std::vector<std::complex<double>> ContractProtonEPAHardSources(LORENTZSCALAR                        &lts,
                                                               const std::vector<HelicityComponent> &components,
                                                               const EPAHardFrame                   &frame) {
  struct HardGroup {
    std::vector<int>                    outgoing;
    std::size_t                         color = 0;
    std::array<std::complex<double>, 4> hard  = {};
  };

  std::vector<HardGroup> groups;
  for (const auto &component : components) {
    if (std::abs(component.incoming[0]) != 1 || std::abs(component.incoming[1]) != 1) {
      throw std::invalid_argument(
          "mg5helas::ContractProtonEPAHardSources: photon helicity is not "
          "transverse");
    }
    auto found = std::find_if(groups.begin(), groups.end(), [&](const HardGroup &group) {
      return group.color == component.color && group.outgoing == component.outgoing;
    });
    if (found == groups.end()) {
      groups.push_back({component.outgoing, component.color, {}});
      found = std::prev(groups.end());
    }
    const std::size_t pair = spin::BinaryPairHelicityIndexX2(component.incoming[0], component.incoming[1]);
    found->hard[pair] += component.value;
  }

  const auto        upper        = ProtonEPASources(lts, 1, frame);
  const auto        lower        = ProtonEPASources(lts, 2, frame);
  ScreeningMetadata metadata     = lts.hamp.metadata;
  metadata.spin_basis            = ScreeningSpinBasis::ProtonHelicity;
  metadata.amplitude_normalization = 0.25;
  metadata.forward_noflip        = false;
  metadata.spin_rows             = 16;
  lts.hamp.Configure(metadata);
  lts.hamp.layout.Clear();

  std::vector<std::complex<double>> out;
  out.reserve(16 * groups.size());
  for (std::size_t initial = 0; initial < 4; ++initial) {
    const auto initial_helicity = spin::BinaryPairHelicityLabelsX2(initial);
    for (std::size_t final = 0; final < 4; ++final) {
      const auto        final_helicity = spin::BinaryPairHelicityLabelsX2(final);
      const std::size_t row1           = spin::BinaryPairHelicityIndex(spin::BinaryHelicityIndexX2(initial_helicity[0]),
                                                                       spin::BinaryHelicityIndexX2(final_helicity[0]));
      const std::size_t row2           = spin::BinaryPairHelicityIndex(spin::BinaryHelicityIndexX2(initial_helicity[1]),
                                                                       spin::BinaryHelicityIndexX2(final_helicity[1]));
      for (const auto &group : groups) {
        out.push_back(KroneckerBilinearProduct(upper.Row(row1), lower.Row(row2), group.hard));
      }
    }
  }
  return out;
}

}  // namespace

// Contract one fixed hard tensor with frame-covariant node-local EPA currents
std::vector<std::complex<double>> ContractEPAHardPhotonSources(LORENTZSCALAR                        &lts,
                                                               const std::vector<HelicityComponent> &components,
                                                               const EPAHardFrame                   &frame) {
  if (!frame.valid) { return {}; }
  if (UsesProtonEPAHardSources(lts)) { return ContractProtonEPAHardSources(lts, components, frame); }
  std::array<M4Vec, 2> source;
  if (!TransformEPAHardSources(frame, lts.q1, lts.q2, source)) { return {}; }
  return ContractEPAHardPhotonSources(lts, components, source[0], source[1], frame.incoming[0], frame.incoming[1]);
}

// Contract the event-local hard tensor and preserve its physical color rows
MatrixElementEvaluation ContractEPAHard(LORENTZSCALAR &lts, const EPAHardTensor &hard) {
  lts.hamp.clear();
  lts.hard_color_flows.clear();
  if (!hard.Ready()) { return {EvaluationStatus::AmplitudeFailure, 0.0}; }
  auto amplitudes = ContractEPAHardPhotonSources(lts, hard.amplitude, hard.frame);
  if (amplitudes.empty() || !gra::AllFinite(amplitudes)) { return {EvaluationStatus::AmplitudeFailure, 0.0}; }
  for (const auto &flow : hard.color) {
    auto values = ContractEPAHardPhotonSources(lts, flow.components, hard.frame);
    if (values.empty() || !gra::AllFinite(values)) {
      lts.hard_color_flows.clear();
      return {EvaluationStatus::AmplitudeFailure, 0.0};
    }
    lts.hard_color_flows.push_back({std::move(values), flow.external});
  }
  const double amp2 = gra::SquaredNorm(amplitudes);
  if (!std::isfinite(amp2)) { return {EvaluationStatus::AmplitudeFailure, 0.0}; }
  lts.hamp = std::move(amplitudes);
  gra::Scale(lts.hamp, hard.normalization);
  return {EvaluationStatus::Success, amp2};
}

}  // namespace mg5helas
}  // namespace gra
