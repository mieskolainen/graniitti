// Conditional nuclear impulse fragmentation through the common HepMC3 reaction API
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Nuclear/MIncoherent.h"

#include <cmath>

#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/Math/MIntegration.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MException.h"
#include "HepMC3/Attribute.h"
#include "HepMC3/GenHeavyIon.h"
#include "HepMC3/GenVertex.h"

namespace gra::nuclear {

using gra::aux::indices;

namespace {

// Compute the continuum energy-transfer interval from the on-shell struck nucleon
std::array<double, 2> Interval(double ma, double m, double mr, double p, double t, double fermi) {
  const double eps = ma - std::hypot(mr, p), a = eps * eps - p * p;
  if (!(a > 0.0) || !(t < 0.0)) { return {}; }
  const double delta = (m * m - a - t) / 2.0, d = p * std::sqrt(delta * delta - a * t);
  const double minimum = std::max(ma, m + mr);
  // [REFERENCE: Bodek and Cai, Eur. Phys. J. C 79 (2019) 293, Appendix B, occupied Fermi sphere and Pauli blocking]
  return {
      std::max({(eps * delta - d) / a, ((minimum - ma) * (minimum + ma) - t) / (2.0 * ma), std::hypot(m, fermi) - eps}),
      (eps * delta + d) / a};
}

}  // namespace

// Match a uniform Fermi sphere to the production matter rms radius
MIncoherent::MIncoherent(std::shared_ptr<const MUPC> upc) {
  if (upc == nullptr) { throw std::invalid_argument("MIncoherent: missing UPC model"); }
  const auto& param = upc->Param();
  if (param.impulse_nodes < 16 || param.impulse_nodes > 4096) {
    throw std::invalid_argument("MIncoherent: invalid impulse quadrature");
  }
  rule_ = math::GaussLegendreRule(param.impulse_nodes, 0.0, 1.0);
  for (const auto& i : indices(decay_)) {
    const auto* nucleus = upc->Nucleus(i + 1);
    if (nucleus == nullptr) { continue; }
    const auto& id = nucleus->ID();
    if (id.a <= 4 || id.z == 0 || id.z == id.a || id.anti || id.lambda || id.isomer) {
      throw std::invalid_argument("MIncoherent: requires an ordinary ground state nucleus with A > 4");
    }
    charge_[i]        = upc->Photon(i + 1) && param.emission[i] != CoherenceType::Coherent;
    const bool matter = upc->Photo(i + 1) && param.target[i] != CoherenceType::Coherent;
    if (charge_[i] && matter) {
      throw std::invalid_argument(
          "MIncoherent: simultaneous charge and matter excitation needs a joint external response");
    }
    active_[i]          = charge_[i] || matter;
    const double radius = std::sqrt(5.0 / 3.0) * nucleus->MatterDensity().Rms() / std::cbrt(id.a);
    decay_[i] = std::make_shared<MEvaporation>(param.mass, radius, param.emd.gdr, param.decay_nodes,
                                               param.decay_steps);
    // Check source masses at initialization rather than failing every sampled point
    try {
      for (const auto z : {id.z, id.z - 1}) { (void)decay_[i]->Mass(id.a - 1, z); }
      if (std::abs(nucleus->Mass() - decay_[i]->Mass(id.a, id.z)) > 1.0e-8) {
        throw std::invalid_argument(
            "MIncoherent: beam masses must use the same bare AME masses as the decay calculation");
      }
    } catch (const PhaseSpaceFailure& error) { throw std::invalid_argument(error.what()); }
  }
  if (upc->Breakup(1)->Enabled() || upc->Breakup(2)->Enabled()) {
    if (charge_[0] || charge_[1]) {
      throw std::invalid_argument(
          "MIncoherent: EMD with an incoherent photon source requires an external joint response");
    }
    emd_ = std::make_shared<const MEMD>(std::move(upc));
  }
}

// Clone only the immutable nuclear calculations and current selectors
std::unique_ptr<MReaction> MIncoherent::Clone() const {
  auto model         = std::make_unique<MIncoherent>(*this);
  model->state_      = {};
  model->excitation_ = {};
  return model;
}

// Describe the physical scope without claiming a transport calculation
HepMC3::GenRunInfo::ToolInfo MIncoherent::Tool() const {
  return {"GRANIITTI nuclear impulse", "1",
          emd_ ? "Fermi gas knockout, resolved amplitude EMD, prompt evaporation, no FSI"
               : "Fermi gas knockout and prompt evaporation, no FSI or additional EMD"};
}

// Propose a normalized continuum mass density with an unbounded rational tail
RecoilMass MIncoherent::SampleMasses(const HepMC3::GenEvent& beams, MRandom& random) {
  state_      = {};
  excitation_ = emd_ ? emd_->Sample(random) : MEMD::State{};
  std::array<double, 2> mass{};
  for (const auto& i : indices(mass)) {
    auto&      state = state_[i];
    const auto beam  = beams.beams().at(i);
    mass[i]          = beam->generated_mass();
    if (excitation_.excitation[i] > 0.0) {
      const auto ion = DecodeNuclearPDG(beam->pid());
      mass[i]        = decay_[i]->Mass(ion.a, ion.z) + excitation_.excitation[i];
    }
    // The half probability is a compensated sector proposal, not a neutron probability
    if (!active_[i] || random.U(0.0, 1.0) < 0.5) { continue; }
    const auto ion     = DecodeNuclearPDG(beam->pid());
    state.excited      = true;
    state.proton       = charge_[i] || random.U(0.0, ion.a) < ion.z;
    state.mass         = decay_[i]->Mass(1, state.proton);
    const double fermi = decay_[i]->Fermi(ion.a, ion.z, state.proton);
    state.p            = fermi * std::cbrt(random.U(0.0, 1.0));
    const double hole =
        (fermi * fermi - state.p * state.p) / (std::hypot(state.mass, fermi) + std::hypot(state.mass, state.p));
    state.residue        = decay_[i]->Mass(ion.a - 1, ion.z - state.proton) + hole + excitation_.excitation[i];
    const double scale   = fermi * fermi / (std::hypot(state.mass, fermi) + state.mass);
    const double u       = random.U(0.0, 1.0);
    const double minimum = std::max(mass[i], state.mass + state.residue);
    mass[i]              = minimum + scale * u / (1.0 - u);
    state.pdf            = (1.0 - u) * (1.0 - u) / scale;
  }
  return {mass, excitation_.excitation, excitation_.channel};
}

// Integrate the angular delta function before normalizing the joint species and hole response
double MIncoherent::Norm(std::size_t leg, const NucleusID& ion, double ma, double t) const {
  double norm = 0.0;
  for (const bool proton : {false, true}) {
    if (charge_[leg] && !proton) { continue; }
    const double abundance = charge_[leg] ? 1.0 : static_cast<double>(proton ? ion.z : ion.a - ion.z) / ion.a;
    if (!(abundance > 0.0)) { continue; }
    const auto&  decay = *decay_[leg];
    const double m = decay.Mass(1, proton), fermi = decay.Fermi(ion.a, ion.z, proton);
    for (const auto& k : indices(rule_.first)) {
      const double p = fermi * rule_.first[k];
      const double mr =
          decay.Mass(ion.a - 1, ion.z - proton) + std::hypot(m, fermi) - std::hypot(m, p) + excitation_.excitation[leg];
      const auto range = Interval(ma, m, mr, p, t, fermi);
      if (!(range[1] > range[0])) { continue; }
      const double eps = ma - std::hypot(mr, p), qt = std::sqrt(-t);
      const double integral = eps * (std::asinh(range[1] / qt) - std::asinh(range[0] / qt)) + std::hypot(range[1], qt) -
                              std::hypot(range[0], qt);
      norm += abundance * 3.0 * rule_.first[k] * rule_.first[k] * rule_.second[k] * integral / (2.0 * p);
    }
  }
  return norm;
}

// Fold the angular energy delta function at fixed t with the sampled mass proposal
// eps = M_A - sqrt(M_R*^2+p^2), omega = (M_X^2-M_A^2-t)/(2 M_A)
// S(omega|p,t) = (eps+omega)/(2 p sqrt(omega^2-t)), normalized over its physical interval
double MIncoherent::Complete(HepMC3::GenEvent& event, const std::array<HepMC3::GenParticlePtr, 2>& forward,
                             MRandom& random) {
  double                     weight = 1.0;
  for (const auto& i : indices(forward)) {
    const auto sector = event.attribute<HepMC3::StringAttribute>("graniitti_upc_final_sector_" + std::to_string(i + 1));
    if (sector == nullptr || (sector->value() != "coherent" && sector->value() != "incoherent")) {
      throw PhaseSpaceFailure("MIncoherent: missing selected hard nuclear sector");
    }
    const bool incoherent = sector->value() == "incoherent";
    if (incoherent != state_[i].excited) { return 0.0; }
    if (active_[i]) { weight *= 2.0; }
  }
  for (const auto& i : indices(forward)) {
    const auto& root = forward[i];
    const auto  beam = root->production_vertex()->particles_in().front();
    if (decay_[i] == nullptr) { continue; }
    root->set_pid(beam->pid());
    root->set_status(1);
    const auto& state = state_[i];
    if (!state.excited) {
      if (excitation_.excitation[i] > 0.0) { decay_[i]->Decay(event, root, random); }
      continue;
    }
    const double ma = beam->generated_mass(), mx = root->generated_mass();
    const auto&  b = beam->momentum();
    const auto&  f = root->momentum();
    const M4Vec  incoming(b.px(), b.py(), b.pz(), b.e());
    const M4Vec  outgoing(f.px(), f.py(), f.pz(), f.e());
    M4Vec        recoil = outgoing;
    kinematics::LorentzBoost(incoming, ma, recoil, -1);
    const double t     = (outgoing - incoming).M2();
    const double omega = ((mx - ma) * (mx + ma) - t) / (2.0 * ma);
    const double q = std::sqrt(omega * omega - t), p = state.p;
    const double er = std::hypot(state.residue, p), eps = ma - er;
    const auto   ion   = DecodeNuclearPDG(beam->pid());
    const auto   range = Interval(ma, state.mass, state.residue, p, t, decay_[i]->Fermi(ion.a, ion.z, state.proton));
    if (!(omega > range[0] && omega < range[1]) || !(eps + omega > state.mass) || !(p > 0.0)) { return 0.0; }
    const double norm = Norm(i, ion, ma, t);
    if (!(norm > 0.0)) { return 0.0; }
    weight *= (eps + omega) / (2.0 * p * q * norm) * mx / ma / state.pdf;
    const double cosine = (std::pow(eps + omega, 2) - state.mass * state.mass - p * p - q * q) / (2.0 * p * q);
    if (std::abs(cosine) > 1.0) { return 0.0; }
    const double phi = random.U(0.0, 2.0 * math::PI), transverse = p * std::sqrt(1.0 - cosine * cosine);
    M4Vec        spectator(-transverse * std::cos(phi), -transverse * std::sin(phi), -p * cosine, er);
    spectator.Rotate(M4Vec(0.0, 0.0, 1.0, 0.0).RotationTo(recoil.P3()));
    kinematics::LorentzBoost(incoming, ma, spectator, 1);
    auto vertex = std::make_shared<HepMC3::GenVertex>(root->production_vertex()->position());
    root->set_status(2);
    vertex->add_particle_in(root);
    auto nucleon =
        std::make_shared<HepMC3::GenParticle>(aux::M4Vec2HepMC3(outgoing - spectator), state.proton ? 2212 : 2112, 1);
    auto residue = std::make_shared<HepMC3::GenParticle>(aux::M4Vec2HepMC3(spectator),
                                                         EncodeNuclearPDG(ion.a - 1, ion.z - state.proton), 1);
    nucleon->set_generated_mass(state.mass);
    residue->set_generated_mass(state.residue);
    vertex->add_particle_out(nucleon);
    vertex->add_particle_out(residue);
    event.add_vertex(vertex);
    decay_[i]->Decay(event, residue, random);
  }
  event.add_attribute("graniitti_nuclear_additional_emd",
                      std::make_shared<HepMC3::StringAttribute>(emd_ ? "resolved_amplitude_tcm" : "absent"));
  return weight;
}

}  // namespace gra::nuclear
