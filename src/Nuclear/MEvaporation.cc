// Statistical nuclear particle and photon emission with exact recoil kinematics
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Nuclear/MEvaporation.h"

#include <algorithm>
#include <cmath>
#include <vector>

#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/Math/MSpecialFunctions.h"
#include "Graniitti/Photon/MQED.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MException.h"
#include "HepMC3/GenVertex.h"

namespace gra::nuclear {

// Prepare immutable nuclear masses and validate independent evaporation inputs
MEvaporation::MEvaporation(const MassParam& mass, double radius, const GDRParam& gdr, unsigned int nodes,
                           unsigned int steps)
    : mass_(mass), radius_(radius), gdr_(gdr), nodes_(nodes), steps_(steps) {
  ValidateGDR(gdr_);
  if (!std::isfinite(radius) || !(radius > 0.0) || nodes < 16 || nodes > 4096 || steps == 0) {
    throw std::invalid_argument("MEvaporation: invalid radius or quadrature");
  }
}

// Compute a mass using the same evaluated and theoretical prescription for every decay channel
double MEvaporation::Mass(unsigned int a, unsigned int z) const { return mass_.Mass(a, z); }

// Compute the occupied species momentum sphere from the same nuclear volume
double MEvaporation::Fermi(unsigned int a, unsigned int z, bool proton) const {
  const double volume = 4.0 * math::PI * std::pow(radius_, 3) * a / 3.0;
  return PDG::GeV2fm * std::cbrt(3.0 * math::PI * math::PI * (proton ? z : a - z) / volume);
}

// Compute the logarithm of the inverse Laplace transform of exp(a/beta), excluding the discrete ground state
double MEvaporation::LogDensity(unsigned int a, unsigned int z, double energy) const {
  if (a <= 4 || energy < 0.0) { return -std::numeric_limits<double>::infinity(); }
  double density = 0.0;
  for (const bool proton : {false, true}) {
    const auto count = proton ? z : a - z;
    if (count == 0) { continue; }
    const double p = Fermi(a, z, proton), m = Mass(1, proton);
    const double ef = p * p / (std::hypot(m, p) + m);
    density += math::PI * math::PI * count / (4.0 * ef);
  }
  const double root = std::sqrt(density * energy);
  return std::log(density) +
         (root * root > std::numeric_limits<double>::epsilon() ? math::LogBesselI1(2.0 * root) - std::log(root) : 0.0);
}

// Append one isotropic two-body emission with exact recoil and explicit daughter masses
HepMC3::GenParticlePtr MEvaporation::Emit(HepMC3::GenEvent& event, const HepMC3::GenParticlePtr& parent, int id1,
                                          double mass1, int id2, double mass2, MRandom& random) const {
  const double       mass = parent->generated_mass();
  const auto&        q    = parent->momentum();
  std::vector<M4Vec> daughter;
  const auto         weight =
      kinematics::TwoBodyPhaseSpace(M4Vec(q.px(), q.py(), q.pz(), q.e()), mass, {mass1, mass2}, daughter, random);
  if (!(weight.GetW() > 0.0)) { throw PhaseSpaceFailure("MEvaporation: invalid decay phase space"); }
  auto vertex = std::make_shared<HepMC3::GenVertex>(
      parent->production_vertex() ? parent->production_vertex()->position() : HepMC3::FourVector());
  vertex->add_particle_in(parent);
  parent->set_status(2);
  auto emitted = std::make_shared<HepMC3::GenParticle>(aux::M4Vec2HepMC3(daughter[0]), id1, 1);
  auto residue = std::make_shared<HepMC3::GenParticle>(aux::M4Vec2HepMC3(daughter[1]), id2, 1);
  emitted->set_generated_mass(mass1);
  residue->set_generated_mass(mass2);
  vertex->add_particle_out(emitted);
  vertex->add_particle_out(residue);
  event.add_vertex(vertex);
  return residue;
}

// Sample normalized detailed balance spectra, including the daughter ground state
// [REFERENCE: Weisskopf, Phys. Rev. 52 (1937) 295, statistical evaporation and inverse capture]
void MEvaporation::Decay(HepMC3::GenEvent& event, HepMC3::GenParticlePtr parent, MRandom& random) const {
  if (parent == nullptr || parent->parent_event() != &event || !IsNuclearPDG(parent->pid()) || parent->end_vertex()) {
    throw PhaseSpaceFailure("MEvaporation: require an undecayed nucleus in the supplied event");
  }
  const auto id = DecodeNuclearPDG(parent->pid());
  if (id.anti || id.lambda || id.isomer) { throw PhaseSpaceFailure("MEvaporation: unsupported nuclear species"); }
  struct Channel {
    unsigned int        a, z;
    int                 id, spin;
    double              mass, ground, maximum;
    std::vector<double> rate;
  };
  // Channels share one density convention and one inverse capture normalization
  const std::array<std::array<int, 4>, 7> ejectile = {{{0, 0, 22, 2},
                                                       {1, 0, 2112, 2},
                                                       {1, 1, 2212, 2},
                                                       {2, 1, 1000010020, 3},
                                                       {3, 1, 1000010030, 2},
                                                       {3, 2, 1000020030, 2},
                                                       {4, 2, 1000020040, 1}}};
  for (unsigned int step = 0; step < steps_; ++step) {
    const auto   ion  = DecodeNuclearPDG(parent->pid());
    const double mass = parent->generated_mass(), ground = Mass(ion.a, ion.z);
    parent->set_status(1);
    const double tolerance = 64.0 * std::numeric_limits<double>::epsilon() * mass;
    if (mass < ground - tolerance) { throw PhaseSpaceFailure("MEvaporation: parent below its ground state"); }
    const bool ground_state = mass - ground <= tolerance;
    // Resolve prompt particle unbound He-5, Li-5 and Be-8 ground states
    // [REFERENCE: Tilley et al., Nucl. Phys. A 708 (2002) 3, A = 5 ground states]
    // [REFERENCE: Tilley et al., Nucl. Phys. A 745 (2004) 155, Be-8 ground state]
    if (ground_state && ((ion.a == 5 && (ion.z == 2 || ion.z == 3)) || (ion.a == 8 && ion.z == 4))) {
      const auto& emitted = ejectile[ion.a == 8 ? 6 : ion.z - 1];
      Emit(event, parent, emitted[2], Mass(emitted[0], emitted[1]), 1000020040, Mass(4, 2), random);
      return;
    }
    if (ground_state) { return; }
    std::vector<Channel> channels;
    std::vector<double>  weights;
    // Integrate in residual excitation, with an explicit unit-weight ground state
    for (const auto& item : ejectile) {
      const auto da = static_cast<unsigned int>(item[0]), dz = static_cast<unsigned int>(item[1]);
      if (ion.a <= da || ion.z < dz || ion.a - ion.z < da - dz) { continue; }
      const unsigned int a = ion.a - da, z = ion.z - dz;
      // Free neutron or proton systems do not form bound evaporation remnants
      if (a > 1 && (z == 0 || z == a)) { continue; }
      Channel      c{a, z, item[2], item[3], da ? Mass(da, dz) : 0.0, Mass(a, z), 0.0, {}};
      const double radius  = radius_ * (std::cbrt(a) + std::cbrt(da));
      const double barrier = da ? qed::alpha_QED() * PDG::GeV2fm * dz * z / radius : 0.0;
      c.maximum            = mass - c.mass - c.ground - barrier;
      if (!(c.maximum > 0.0)) { continue; }
      const Dipole dipole(a, z, gdr_);
      // Sharp absorption and classical Coulomb capture, with no branching parameters
      const auto rate = [&](double excitation) {
        const double daughter = c.ground + excitation;
        if (da == 0) {
          const double photon = (mass - daughter) * (mass + daughter) / (2.0 * mass);
          return photon * photon * 0.1 * dipole.Sigma(photon) * daughter / mass;
        }
        const double kinetic = mass - c.mass - daughter;
        const double mu      = c.mass * daughter / (c.mass + daughter);
        return c.spin * mu * math::PI * radius * radius * std::max(0.0, kinetic - barrier);
      };
      weights.push_back(std::log(rate(0.0)));
      const double width = c.maximum / nodes_;
      for (unsigned int k = 0; k <= nodes_; ++k) {
        const double energy = k * width;
        c.rate.push_back(std::log(rate(energy)) + LogDensity(a, z, energy));
        if (k > 0) {
          const double high = std::max(c.rate[k - 1], c.rate[k]);
          weights.push_back(std::isfinite(high) ? std::log(0.5 * width) + high +
                                                      std::log1p(std::exp(std::min(c.rate[k - 1], c.rate[k]) - high))
                                                : high);
        }
      }
      channels.push_back(std::move(c));
    }
    const double maximum =
        weights.empty() ? -std::numeric_limits<double>::infinity() : *std::max_element(weights.begin(), weights.end());
    for (double& weight : weights) { weight = std::exp(weight - maximum); }
    const double total = gra::Sum(weights);
    if (!(total > 0.0) || !std::isfinite(total)) { throw PhaseSpaceFailure("MEvaporation: no finite decay width"); }
    double      select = random.U(0.0, total);
    std::size_t bin    = 0;
    while (bin + 1 < weights.size() && select >= weights[bin]) { select -= weights[bin++]; }
    const auto&        c          = channels[bin / (nodes_ + 1)];
    const unsigned int k          = bin % (nodes_ + 1);
    double             excitation = 0.0;
    if (k > 0) {
      const double high = std::max(c.rate[k - 1], c.rate[k]);
      const double left = std::exp(c.rate[k - 1] - high), right = std::exp(c.rate[k] - high), u = random.U(0.0, 1.0);
      // Invert the linear density within the selected integration interval
      const double denominator = left + std::sqrt((1.0 - u) * left * left + u * right * right);
      const double fraction    = denominator > 0.0 ? u * (left + right) / denominator : u;
      excitation               = (k - 1 + fraction) * c.maximum / nodes_;
    }
    const int residue = c.a == 1 ? (c.z ? 2212 : 2112) : EncodeNuclearPDG(c.a, c.z);
    parent            = Emit(event, parent, c.id, c.mass, residue, c.ground + excitation, random);
    if (c.a == 1) { return; }
  }
  throw PhaseSpaceFailure("MEvaporation: decay did not terminate");
}

}  // namespace gra::nuclear
