// HepMC event-record construction shared by physical processes
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Process/MEventRecord.h"

#include <array>
#include <cmath>

#include "Graniitti/MGlobals.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Nuclear/MRecord.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MException.h"
#include "HepMC3/Attribute.h"
#include "HepMC3/GenPdfInfo.h"
#include "HepMC3/GenVertex.h"

namespace gra::record {

using aux::indices;

namespace {

// Compute the partonic momentum fraction to publish in HepMC PDF metadata
double OutputPartonX(const gra::LORENTZSCALAR &lts, int leg) {
  if (lts.diff_xhard1 > 0.0 && lts.diff_xhard2 > 0.0) { return (leg == 1) ? lts.diff_xhard1 : lts.diff_xhard2; }
  return (leg == 1) ? lts.x1 : lts.x2;
}

// Attach hard-process PDF metadata when a partonic event view is available
void AttachPdfInfo(const gra::LORENTZSCALAR &lts, HepMC3::GenEvent &evt) {
  // Mark kT-EPA records so the LHE writer retains the physical pp event
  if (lts.exact_forward_photon_kinematics) {
    evt.add_attribute("graniitti_exact_forward_photon_kinematics", std::make_shared<HepMC3::IntAttribute>(1));
  }
  const double x1 = OutputPartonX(lts, 1);
  const double x2 = OutputPartonX(lts, 2);
  if (lts.id1 == 0 || lts.id2 == 0) { return; }
  if (!(x1 > 0.0) || !(x2 > 0.0) || !(lts.muF > 0.0)) { return; }
  if (!std::isfinite(lts.pdf_xf1) || !std::isfinite(lts.pdf_xf2)) { return; }

  auto pdf_info = std::make_shared<HepMC3::GenPdfInfo>();
  pdf_info->set(lts.id1, lts.id2, x1, x2, lts.muF, lts.pdf_xf1, lts.pdf_xf2);
  evt.set_pdf_info(pdf_info);
}

// Attach standard event-generator scales needed by the LHE writer
void AttachGeneratorScales(const gra::LORENTZSCALAR &lts, HepMC3::GenEvent &evt) {
  if (std::isfinite(lts.muF) && lts.muF > 0.0) {
    evt.add_attribute("graniitti_mu_f", std::make_shared<HepMC3::DoubleAttribute>(lts.muF));
  }
  if (std::isfinite(lts.muR) && lts.muR > 0.0) {
    evt.add_attribute("graniitti_mu_r", std::make_shared<HepMC3::DoubleAttribute>(lts.muR));
  }
  if (std::isfinite(lts.scalup) && lts.scalup > 0.0) {
    evt.add_attribute("graniitti_scalup", std::make_shared<HepMC3::DoubleAttribute>(lts.scalup));
  }
}

// Compute true for a Durham colored central-exclusive event
bool IsDurhamColoredCEP(const gra::LORENTZSCALAR &lts) {
  return !lts.hard_diff1 && !lts.hard_diff2 && lts.id1 == PDG::PDG_gluon && lts.id2 == PDG::PDG_gluon &&
         !lts.hard_color_flows.empty();
}

// Attach Durham CEP metadata without enabling the generic hard-PDF LHE view
void AttachDurhamCEPInfo(const gra::LORENTZSCALAR &lts, HepMC3::GenEvent &evt) {
  if (!IsDurhamColoredCEP(lts)) { return; }

  evt.add_attribute("graniitti_elastic_colored_cep", std::make_shared<HepMC3::IntAttribute>(1));
  if (std::isfinite(lts.alphaQCD) && lts.alphaQCD > 0.0) {
    evt.add_attribute("graniitti_alpha_qcd", std::make_shared<HepMC3::DoubleAttribute>(lts.alphaQCD));
  }
}

// Propagate one unstable particle from its production vertex in event units
M4Vec DecayPosition(const MDecayBranch &branch, const HepMC3::GenParticlePtr &mother,
                   const HepMC3::GenEvent &evt, MRandom &random, bool displaced) {
  const auto production = mother->production_vertex();
  M4Vec position = aux::HepMC2M4Vec(production ? production->position() : evt.event_pos());
  if (!displaced) { return position; }
  if (!std::isfinite(branch.p.tau) || branch.p.tau < 0.0) {
    throw std::invalid_argument("record::DecayPosition: invalid particle lifetime");
  }
  if (branch.p.tau > 0.0) {
    if (!std::isfinite(branch.p4.M2()) || !(branch.p4.M2() > 0.0) || !(branch.p4.E() > 0.0)) {
      throw PhaseSpaceFailure("record::DecayPosition: invalid unstable particle momentum");
    }
    const double tau = branch.p.tau * random.ExpRandom(1.0);
    const double scale = evt.length_unit() == HepMC3::Units::MM ? 1e3 : 1e2;
    position += branch.p4.PropagatePosition(tau, scale);
  }
  if (!std::isfinite(position.X()) || !std::isfinite(position.Y()) ||
      !std::isfinite(position.Z()) || !std::isfinite(position.E())) {
    throw PhaseSpaceFailure("record::DecayPosition: non-finite decay position");
  }
  return position;
}

}  // namespace

// Attach standard HepMC color-flow attributes to a generated particle
void AttachColorFlow(const MParticle &particle, const HepMC3::GenParticlePtr &generated) {
  if (particle.color_flow.empty()) { return; }
  generated->add_attribute("flow1", std::make_shared<HepMC3::IntAttribute>(particle.color_flow.flow1));
  generated->add_attribute("flow2", std::make_shared<HepMC3::IntAttribute>(particle.color_flow.flow2));
}

// Sample lifetimes and write a decay branch with an optional prompt root
void WriteBranch(MDecayBranch &branch, const HepMC3::GenParticlePtr &mother, HepMC3::GenEvent &evt,
                 MRandom &random, bool displaced) {
  if (branch.legs.empty()) { return; }
  branch.decay_position = DecayPosition(branch, mother, evt, random, displaced);
  auto vertex = std::make_shared<HepMC3::GenVertex>(aux::M4Vec2HepMC3(branch.decay_position));
  vertex->add_particle_in(mother);
  evt.add_vertex(vertex);
  for (auto &leg : branch.legs) {
    const int status   = leg.legs.empty() ? PDG::PDG_STABLE : PDG::PDG_DECAY;
    auto      particle = std::make_shared<HepMC3::GenParticle>(aux::M4Vec2HepMC3(leg.p4), leg.p.pdg, status);
    vertex->add_particle_out(particle);
    AttachColorFlow(leg.p, particle);
    WriteBranch(leg, particle, evt, random);
  }
}

// Event output recording (to HepMC containers)
//
bool WriteCentral(MProcessState &state, HepMC3::GenEvent &evt) {
  AttachGeneratorScales(state.lts, evt);
  AttachDurhamCEPInfo(state.lts, evt);
  AttachPdfInfo(state.lts, evt);
  const auto &upc = state.lts.upc_event != nullptr ? state.lts.upc_event : state.lts.upc_model;

  std::array<M4Vec, 2> forward    = {state.lts.pfinal[1], state.lts.pfinal[2]};
  std::array<M4Vec, 2> propagator = {state.lts.q1, state.lts.q2};
  if (upc != nullptr) {
    // Keep inclusive recoils unresolved until the selected reaction completes their decay
    nuclear::AttachRecord(*upc, state.upc_final, evt);
    if (state.upc_param.additional_emd && state.nuclear_final) {
      for (const auto& i : indices(state.upc_photo)) {
        evt.add_attribute("graniitti_upc_photo_" + std::to_string(i),
                          std::make_shared<HepMC3::DoubleAttribute>(state.upc_photo[i]));
      }
    }
  }

  // Initial states (4-momentum, pdg-id, status code)
  HepMC3::GenParticlePtr gen_p1 = std::make_shared<HepMC3::GenParticle>(gra::aux::M4Vec2HepMC3(state.lts.pbeam1),
                                                                        state.lts.beam1.pdg, PDG::PDG_BEAM);
  HepMC3::GenParticlePtr gen_p2 = std::make_shared<HepMC3::GenParticle>(gra::aux::M4Vec2HepMC3(state.lts.pbeam2),
                                                                        state.lts.beam2.pdg, PDG::PDG_BEAM);
  gen_p1->set_generated_mass(state.lts.beam1.mass);
  gen_p2->set_generated_mass(state.lts.beam2.mass);

  // Final states (protons/excited system)
  int PDG_ID1 = state.lts.beam1.pdg;
  int PDG_ID2 = state.lts.beam2.pdg;

  int        PDG_status1 = PDG::PDG_STABLE;
  int        PDG_status2 = PDG::PDG_STABLE;
  const bool forward_branch1 = state.lts.excite1 || !state.lts.decayforward1.legs.empty();
  const bool forward_branch2 = state.lts.excite2 || !state.lts.decayforward2.legs.empty();

  if (forward_branch1) {
    PDG_ID1     = std::abs(PDG::PDG_NSTAR) * math::sign(state.lts.beam1.pdg);
    PDG_status1 = PDG::PDG_INTERMEDIATE;
  }
  if (forward_branch2) {
    PDG_ID2     = std::abs(PDG::PDG_NSTAR) * math::sign(state.lts.beam2.pdg);
    PDG_status2 = PDG::PDG_INTERMEDIATE;
  }
  // Inclusive nuclear systems are documentation entries until a backend resolves them
  if (upc != nullptr && upc->Nucleus(1) != nullptr) { PDG_status1 = 3; }
  if (upc != nullptr && upc->Nucleus(2) != nullptr) { PDG_status2 = 3; }

  HepMC3::GenParticlePtr gen_p1f =
      std::make_shared<HepMC3::GenParticle>(gra::aux::M4Vec2HepMC3(forward[0]), PDG_ID1, PDG_status1);
  HepMC3::GenParticlePtr gen_p2f =
      std::make_shared<HepMC3::GenParticle>(gra::aux::M4Vec2HepMC3(forward[1]), PDG_ID2, PDG_status2);
  if (state.lts.forward_mass2[0] >= 0.0) { gen_p1f->set_generated_mass(std::sqrt(state.lts.forward_mass2[0])); }
  if (state.lts.forward_mass2[1] >= 0.0) { gen_p2f->set_generated_mass(std::sqrt(state.lts.forward_mass2[1])); }

  // -------------------------------------------------------------------

  // Propagator 1 and 2
  // it is ill-posed to try classify pomeron/gamma/gluon etc. here, thus,
  // we tag only a generic propagator ID in the record
  HepMC3::GenParticlePtr gen_q1 = std::make_shared<HepMC3::GenParticle>(gra::aux::M4Vec2HepMC3(propagator[0]),
                                                                        PDG::PDG_propagator, PDG::PDG_INTERMEDIATE);
  HepMC3::GenParticlePtr gen_q2 = std::make_shared<HepMC3::GenParticle>(gra::aux::M4Vec2HepMC3(propagator[1]),
                                                                        PDG::PDG_propagator, PDG::PDG_INTERMEDIATE);

  // -------------------------------------------------------------------

  // Central system / resonance
  int central_pdg = state.lts.process.root_resonance_pdg != 0 ? state.lts.process.root_resonance_pdg : PDG::PDG_system;
  if (state.lts.process.root_decay_mode == RootDecayMode::Isolated && state.lts.process.RESONANCES.size() == 1) {
    central_pdg = state.lts.process.RESONANCES.begin()->second.p.pdg;
  }
  HepMC3::GenParticlePtr gen_q = std::make_shared<HepMC3::GenParticle>(gra::aux::M4Vec2HepMC3(state.lts.pfinal[0]),
                                                                       central_pdg, PDG::PDG_INTERMEDIATE);

  // ====================================================================
  // Construct vertices

  // Upper proton-propagator-proton
  HepMC3::GenVertexPtr v1 = std::make_shared<HepMC3::GenVertex>();
  v1->add_particle_in(gen_p1);
  v1->add_particle_out(gen_p1f);
  v1->add_particle_out(gen_q1);

  // Lower proton-propagator-proton
  HepMC3::GenVertexPtr v2 = std::make_shared<HepMC3::GenVertex>();
  v2->add_particle_in(gen_p2);
  v2->add_particle_out(gen_p2f);
  v2->add_particle_out(gen_q2);

  // Propagator-Propagator-System vertex
  HepMC3::GenVertexPtr v3 = std::make_shared<HepMC3::GenVertex>();
  v3->add_particle_in(gen_q1);
  v3->add_particle_in(gen_q2);
  v3->add_particle_out(gen_q);

  evt.add_vertex(v1);
  evt.add_vertex(v2);
  evt.add_vertex(v3);
  if (upc != nullptr) {
    const std::array roots = {gen_p1f, gen_p2f};
    for (const auto& i : indices(roots)) {
      roots[i]->add_attribute("graniitti_upc_leg", std::make_shared<HepMC3::IntAttribute>(i + 1));
      roots[i]->add_attribute("graniitti_upc_final_sector", std::make_shared<HepMC3::StringAttribute>(
          state.upc_final.valid ? nuclear::CoherenceName(state.upc_final.leg[i]) : "unspecified"));
    }
  }

  // ====================================================================
  // System->Decay products vertex

  HepMC3::GenVertexPtr v4 = std::make_shared<HepMC3::GenVertex>();
  evt.add_vertex(v4);

  // Add resonance in
  v4->add_particle_in(gen_q);

  // Add direct daughters
  for (const auto &i : indices(state.lts.decaytree)) {
    const int STATE = (state.lts.decaytree[i].legs.size() > 0) ? PDG::PDG_DECAY : PDG::PDG_STABLE;

    HepMC3::GenParticlePtr particle = std::make_shared<HepMC3::GenParticle>(
        gra::aux::M4Vec2HepMC3(state.lts.decaytree[i].p4), state.lts.decaytree[i].p.pdg, STATE);
    v4->add_particle_out(particle);
    AttachColorFlow(state.lts.decaytree[i].p, particle);

    WriteBranch(state.lts.decaytree[i], particle, evt, state.random);
  }

  // Upper proton excitation
  if (forward_branch1 && state.lts.decayforward1.legs.size() != 0) {
    WriteBranch(state.lts.decayforward1, gen_p1f, evt, state.random, false);
  }

  // Lower proton excitation
  if (forward_branch2 && state.lts.decayforward2.legs.size() != 0) {
    WriteBranch(state.lts.decayforward2, gen_p2f, evt, state.random, false);
  }

  return true;
}


}  // namespace gra::record
