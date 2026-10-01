// Tensor Pomeron process classification and decay validation
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Tensor/MTensorProcess.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "Graniitti/Math/MMath.h"
#include "Graniitti/Process/MDecayTree.h"
#include "Graniitti/Tensor/MTensorPomeron.h"

using gra::math::pow2;

namespace gra::tensor {

// Classify resonance quantum numbers supported by the implemented amplitudes
TensorResonanceType ClassifyTensorResonance(const MParticle &particle) {
  if (particle.chargeX3 == 0 && particle.C == 1) {
    if (particle.spinX2 == 0 && particle.P == 1) {
      return TensorResonanceType::Scalar;
    }
    if (particle.spinX2 == 0 && particle.P == -1) {
      return TensorResonanceType::Pseudoscalar;
    }
    if (particle.spinX2 == 2 && particle.P == 1) {
      return TensorResonanceType::AxialVector;
    }
    if (particle.spinX2 == 4 && particle.P == 1) {
      return TensorResonanceType::Tensor;
    }
  }
  if (particle.chargeX3 == 0 && particle.spinX2 == 2 && particle.P == -1 &&
      particle.C == -1) {
    return TensorResonanceType::Vector;
  }
  if (particle.chargeX3 == 0 && particle.spinX2 == 6 && particle.P == -1 && particle.C == -1) {
    return TensorResonanceType::Spin3;
  }
  throw std::invalid_argument(
      "MTensorPomeron::ME3: unsupported resonance quantum numbers for PDG = " +
      std::to_string(particle.pdg) +
      "; supported neutral J^PC states are 0++, 0-+, 1--, 1++, 2++ and 3--");
}

// Compute true when the central state is one direct two-body topology
bool TensorDirectTwoBody(const std::vector<MDecayBranch> &tree) {
  return tree.size() == 2 && tree[0].legs.empty() && tree[1].legs.empty();
}

// Compute true when both top-level branches contain exact two-body cascades
bool TensorPairedTwoBodyCascades(const std::vector<MDecayBranch> &tree) {
  if (tree.size() != 2 || tree[0].legs.size() != 2 ||
      tree[1].legs.size() != 2) {
    return false;
  }
  for (const auto &branch : tree) {
    if (!branch.legs[0].legs.empty() || !branch.legs[1].legs.empty()) {
      return false;
    }
  }
  return true;
}

// Require finite nonnegative equal stable daughter masses before event
// kinematics
void ValidateStableVectorDaughterMasses(const MParticle &first,
                                        const MParticle &second,
                                        const std::string &context) {
  if (!std::isfinite(first.mass) || !std::isfinite(second.mass) ||
      first.mass < 0.0 || second.mass < 0.0) {
    throw std::invalid_argument(
        context + ": vector daughters require finite nonnegative masses");
  }

  const double first_mass2 = pow2(first.mass);
  const double second_mass2 = pow2(second.mass);
  const double scale = std::max({1.0, first_mass2, second_mass2});
  if (std::abs(first_mass2 - second_mass2) > 1e-10 * scale) {
    throw std::invalid_argument(
        context + ": vector daughters require equal particle masses");
  }
}

// Compute true for a direct charged-fermion pair supported by QED ME4
bool TensorQEDDirectPair(const std::vector<MDecayBranch> &tree) {
  if (!TensorDirectTwoBody(tree)) {
    return false;
  }
  const int first_pdg = tree[0].p.pdg;
  const int abs_pdg = std::abs(first_pdg);
  const bool is_quark = abs_pdg >= 1 && abs_pdg <= 6;
  const bool is_charged_lepton =
      abs_pdg == 11 || abs_pdg == 13 || abs_pdg == 15;
  return first_pdg == -tree[1].p.pdg && (is_quark || is_charged_lepton);
}

// Compute true for direct charged-pion and charged-kaon photoproduction
bool TensorPhotoPair(const std::vector<MDecayBranch> &tree) {
  if (!TensorDirectTwoBody(tree)) {
    return false;
  }
  const int pdg = std::abs(tree[0].p.pdg);
  return tree[0].p.pdg == -tree[1].p.pdg && (pdg == PDG::PDG_pip || pdg == PDG::PDG_Kp);
}

// Compute true for a direct hadron pair supported by Tensor Pomeron ME4
bool TensorHadronicDirectPair(const std::vector<MDecayBranch> &tree) {
  if (!TensorDirectTwoBody(tree)) {
    return false;
  }
  const int first_pdg = tree[0].p.pdg;
  const int second_pdg = tree[1].p.pdg;
  const int abs_pdg = std::abs(first_pdg);
  const bool conjugate_pair = first_pdg == -second_pdg;
  const bool spin_zero_pair = tree[0].p.spinX2 == 0 && tree[1].p.spinX2 == 0;
  const bool pion_pair =
      spin_zero_pair && ((conjugate_pair && abs_pdg == 211) ||
                         (first_pdg == 111 && second_pdg == 111));
  const bool kaon_pair =
      spin_zero_pair && conjugate_pair && (abs_pdg == 311 || abs_pdg == 321);
  const bool baryon_pair = tree[0].p.spinX2 == 1 && tree[1].p.spinX2 == 1 &&
                           conjugate_pair && abs_pdg >= 1000;
  const bool vector_pair = tree[0].p.spinX2 == 2 && tree[1].p.spinX2 == 2 &&
                           first_pdg == second_pdg &&
                           (first_pdg == 113 || first_pdg == 333);
  return pion_pair || kaon_pair || baryon_pair || vector_pair;
}

// Compute true for paired vectors with direct pseudoscalar decays
bool TensorHadronicVectorCascades(const std::vector<MDecayBranch> &tree) {
  if (!TensorPairedTwoBodyCascades(tree) || tree[0].p.pdg != tree[1].p.pdg) {
    return false;
  }
  for (const auto &branch : tree) {
    if (branch.p.spinX2 != 2 || branch.legs[0].p.spinX2 != 0 ||
        branch.legs[1].p.spinX2 != 0 ||
        branch.legs[0].p.pdg != -branch.legs[1].p.pdg) {
      return false;
    }
  }
  return true;
}

// Compute true when both direct top-level states have the requested spin
bool TensorDirectSpinPair(const std::vector<MDecayBranch> &tree, int spin_x2) {
  return TensorDirectTwoBody(tree) && tree[0].p.spinX2 == spin_x2 &&
         tree[1].p.spinX2 == spin_x2;
}

// Compute true when one absolute PDG id is present in resolved model data
bool TensorModelContains(const std::vector<int> &pdgs, int pdg) {
  return std::find(pdgs.begin(), pdgs.end(), std::abs(pdg)) != pdgs.end();
}

// Compute the configured daughter id for one resolved vector-meson channel
int TensorVectorDaughterPDG(const LORENTZSCALAR &lts, int vector_pdg) {
  if (!lts.process.TENSOR_MODEL_READY) {
    throw std::invalid_argument(
        "MTensorPomeron process requires resolved model data");
  }
  const auto channel =
      lts.process.TENSOR_VECTOR_DECAY_PDGS.find(std::abs(vector_pdg));
  if (channel == lts.process.TENSOR_VECTOR_DECAY_PDGS.end()) {
    throw std::invalid_argument(
        "MTensorPomeron process has no configured vector channel for PDG " +
        std::to_string(vector_pdg));
  }
  return channel->second;
}

// Validate tune-specific direct or cascaded Tensor Pomeron continuum support
void ValidateTensorContinuumModel(const LORENTZSCALAR &lts) {
  if (!lts.process.TENSOR_MODEL_READY) {
    throw std::invalid_argument(
        "MTensorPomeron continuum process requires resolved model data");
  }
  const auto &tree = lts.decaytree;
  if (TensorDirectTwoBody(tree)) {
    const int pdg = tree[0].p.pdg;
    const int spin_x2 = tree[0].p.spinX2;
    bool supported = false;
    if (spin_x2 == 0) {
      supported =
          TensorModelContains(lts.process.TENSOR_PSEUDOSCALAR_PDGS, pdg);
    } else if (spin_x2 == 1) {
      supported = TensorModelContains(lts.process.TENSOR_BARYON_PDGS, pdg);
    } else if (spin_x2 == 2) {
      static_cast<void>(TensorVectorDaughterPDG(lts, pdg));
      supported = true;
    }
    if (!supported) {
      throw std::invalid_argument(
          "MTensorPomeron continuum process has no configured channel for "
          "PDG " +
          std::to_string(pdg));
    }
    return;
  }

  const int daughter_pdg = TensorVectorDaughterPDG(lts, tree[0].p.pdg);
  for (const auto &branch : tree) {
    const bool daughters_match =
        branch.legs[0].p.pdg == -branch.legs[1].p.pdg &&
        std::abs(branch.legs[0].p.pdg) == daughter_pdg;
    if (!daughters_match || branch.hel.g_decay_TP.empty() ||
        !std::isfinite(branch.hel.g_decay_TP[0]) ||
        !std::isfinite(branch.p.mass) || !std::isfinite(branch.p.width) ||
        !(branch.p.mass > 0.0) || !(branch.p.width > 0.0)) {
      throw std::invalid_argument(
          "MTensorPomeron cascade process does not match the configured "
          "vector decay channel");
    }
    ValidateStableVectorDaughterMasses(
        branch.legs[0].p, branch.legs[1].p,
        "MTensorPomeron continuum cascade process");
  }
}

// Require the vector cascade structure implemented by resonance ME3
void ValidateTensorResonanceVectorCascades(const LORENTZSCALAR &lts) {
  const auto &tree = lts.decaytree;
  if (!TensorPairedTwoBodyCascades(tree)) {
    throw std::invalid_argument(
        "MTensorPomeron resonance process requires two direct vector "
        "cascades");
  }
  for (const auto &branch : tree) {
    (void)TensorVectorDaughterPDG(lts, branch.p.pdg);
    if (branch.p.spinX2 != 2 || branch.legs[0].p.spinX2 != 0 ||
        branch.legs[1].p.spinX2 != 0) {
      throw std::invalid_argument(
          "MTensorPomeron resonance process requires spin-1 vectors "
          "decaying to stable spin-0 pairs");
    }
    if (branch.hel.g_decay_TP.empty() ||
        !std::isfinite(branch.hel.g_decay_TP[0]) ||
        !std::isfinite(branch.p.mass) || !std::isfinite(branch.p.width) ||
        !(branch.p.mass > 0.0) || !(branch.p.width > 0.0)) {
      throw std::invalid_argument(
          "MTensorPomeron resonance cascade process requires a finite "
          "decay coupling and positive vector mass and width");
    }
    ValidateStableVectorDaughterMasses(
        branch.legs[0].p, branch.legs[1].p,
        "MTensorPomeron resonance cascade process");
  }
}

// Validate the particle spins required by the QED continuum amplitude
void ValidateTensorQEDProcess(const LORENTZSCALAR &lts) {
  if (lts.decaytree.size() != 2 || lts.decaytree[0].p.spinX2 != 1 ||
      lts.decaytree[1].p.spinX2 != 1) {
    throw std::invalid_argument(
        "MTensorPomeron QED process requires a spin-1/2 particle pair");
  }
}

// Validate one non-axial Tensor Pomeron resonance decay topology
void ValidateTensorResonanceDecay(const PARAM_RES &resonance,
                                  const LORENTZSCALAR &lts) {
  const auto &tree = lts.decaytree;
  const TensorResonanceType type = ClassifyTensorResonance(resonance.p);
  if (type == TensorResonanceType::AxialVector) {
    return;
  }
  if (tree.size() != 2) {
    throw std::invalid_argument(
        "MTensorPomeron non-axial resonance process requires exactly two "
        "top-level central branches");
  }
  if (resonance.hel_decay.g_decay_TP.empty()) {
    throw std::invalid_argument(
        "MTensorPomeron resonance process requires tensor decay "
        "couplings for PDG " +
        std::to_string(resonance.p.pdg));
  }

  const bool direct_spin_zero = TensorDirectSpinPair(tree, 0);
  const bool direct_spin_one = TensorDirectSpinPair(tree, 2);
  const bool photons = direct_spin_one && tree[0].p.pdg == PDG::PDG_gamma &&
                       tree[1].p.pdg == PDG::PDG_gamma;
  const bool cascaded = !tree[0].legs.empty() || !tree[1].legs.empty();
  if (cascaded) {
    ValidateTensorResonanceVectorCascades(lts);
  }
  const bool vector_pair = (direct_spin_one && !photons) || cascaded;

  if (type == TensorResonanceType::Scalar) {
    if (!direct_spin_zero && !vector_pair) {
      throw std::invalid_argument(
          "MTensorPomeron scalar resonance process supports a direct "
          "spin-0 pair or a massive vector pair");
    }
    if (vector_pair && resonance.hel_decay.g_decay_TP.size() != 2) {
      throw std::invalid_argument(
          "MTensorPomeron scalar to vector-pair process requires two "
          "tensor decay couplings");
    }
    return;
  }
  if (type == TensorResonanceType::Vector || type == TensorResonanceType::Spin3) {
    if (!direct_spin_zero) {
      throw std::invalid_argument(
          "MTensorPomeron vector resonance process requires a direct "
          "spin-0 pair");
    }
    ValidateStableVectorDaughterMasses(
        tree[0].p, tree[1].p, "MTensorPomeron vector resonance process");
    if (type == TensorResonanceType::Vector) { static_cast<void>(TensorVectorDaughterPDG(lts, resonance.p.pdg)); }
    return;
  }
  if (type == TensorResonanceType::Pseudoscalar) {
    if (!photons) {
      throw std::invalid_argument(
          "MTensorPomeron pseudoscalar resonance process requires a "
          "direct photon pair");
    }
    return;
  }
  if (!direct_spin_zero && !vector_pair && !photons) {
    throw std::invalid_argument(
        "MTensorPomeron tensor resonance process has an unsupported decay "
        "topology");
  }
  if ((vector_pair || photons) &&
      resonance.hel_decay.g_decay_TP.size() != 2) {
    throw std::invalid_argument(
        "MTensorPomeron tensor to vector-pair process requires two tensor "
        "decay couplings");
  }
}

// Validate event-dependent resonance data used by the Tensor Pomeron ME3
void ValidateTensorResonanceProcess(const LORENTZSCALAR &lts) {
  if (lts.decaytree.empty()) {
    return;
  }
  if (lts.decaytree.size() < 2) {
    throw std::invalid_argument(
        "MTensorPomeron resonance process requires at least two central "
        "branches");
  }
  if (!lts.process.SPINGEN || !lts.process.SPINDEC) {
    throw std::invalid_argument(
        "MTensorPomeron resonance process requires SPINGEN and SPINDEC");
  }
  if (lts.process.RESONANCES.empty()) {
    throw std::invalid_argument(
        "MTensorPomeron resonance process requires an active resonance");
  }

  bool all_axial = true;
  for (const auto &entry : lts.process.RESONANCES) {
    const TensorResonanceType type = ClassifyTensorResonance(entry.second.p);
    all_axial = all_axial && type == TensorResonanceType::AxialVector;
    ValidateTensorResonanceDecay(entry.second, lts);
  }
  if (lts.process.root_decay_mode == RootDecayMode::Isolated && !all_axial) {
    throw std::invalid_argument(
        "MTensorPomeron isolated resonance sampling requires the generic "
        "axial-vector decay amplitude");
  }
}

// Match the decay-tree structure accepted by one Tensor Pomeron mode
bool ProcessAccepts(const std::vector<MDecayBranch> &tree,
                    MTensorPomeronMode mode) {
  switch (mode) {
  case MTensorPomeronMode::Generic:
    return !tree.empty();
  case MTensorPomeronMode::Resonance:
    return tree.size() >= 2;
  case MTensorPomeronMode::Continuum:
    return TensorHadronicDirectPair(tree) || TensorHadronicVectorCascades(tree);
  case MTensorPomeronMode::ResonanceContinuum:
    return TensorHadronicDirectPair(tree);
  case MTensorPomeronMode::Photo:
    return TensorPhotoPair(tree);
  case MTensorPomeronMode::QED:
    return TensorQEDDirectPair(tree);
  }
  return false;
}

// Compute the stable process identity of one Tensor Pomeron process mode
std::string ProcessName(MTensorPomeronMode mode) {
  switch (mode) {
  case MTensorPomeronMode::Generic:
    return "generic";
  case MTensorPomeronMode::Resonance:
    return "resonance";
  case MTensorPomeronMode::Continuum:
    return "continuum";
  case MTensorPomeronMode::ResonanceContinuum:
    return "resonance_continuum";
  case MTensorPomeronMode::Photo:
    return "photo";
  case MTensorPomeronMode::QED:
    return "qed";
  }
  throw std::invalid_argument(
      "MTensorPomeron::ProcessDefinitionFor: unknown process mode");
}

// Compute the final-state pattern exposed by one Tensor Pomeron mode
std::string ProcessPattern(MTensorPomeronMode mode) {
  switch (mode) {
  case MTensorPomeronMode::Generic:
    return "event Tensor Pomeron final state";
  case MTensorPomeronMode::Resonance:
    return "resonance decays";
  case MTensorPomeronMode::Continuum:
    return "2 hadrons and vector pair decays";
  case MTensorPomeronMode::ResonanceContinuum:
    return "2 hadrons";
  case MTensorPomeronMode::Photo:
    return "pi+ pi-, K+ K-";
  case MTensorPomeronMode::QED:
    return "l+l-, qqbar";
  }
  throw std::invalid_argument(
      "MTensorPomeron::ProcessDefinitionFor: unknown process mode");
}

// Resolve the decay structure for one accepted Tensor Pomeron topology
DecayStructure ProcessDecayStructure(MTensorPomeronMode mode,
                                     const LORENTZSCALAR &lts) {
  if (mode == MTensorPomeronMode::Resonance ||
      (mode == MTensorPomeronMode::ResonanceContinuum &&
       !lts.process.RESONANCES.empty()) ||
      (mode == MTensorPomeronMode::Photo && !lts.process.RESONANCES.empty())) {
    ValidateTensorResonanceProcess(lts);
  }
  if ((mode == MTensorPomeronMode::Continuum ||
       mode == MTensorPomeronMode::ResonanceContinuum ||
       mode == MTensorPomeronMode::Photo || mode == MTensorPomeronMode::QED) &&
      !lts.process.SPINGEN) {
    throw std::invalid_argument(
        "MTensorPomeron continuum process requires SPINGEN");
  }
  if (mode == MTensorPomeronMode::QED) {
    ValidateTensorQEDProcess(lts);
  }
  if (mode == MTensorPomeronMode::Continuum &&
      TensorHadronicVectorCascades(lts.decaytree) && !lts.process.SPINDEC) {
    throw std::invalid_argument(
        "MTensorPomeron cascade process requires SPINDEC");
  }
  if (mode == MTensorPomeronMode::Continuum ||
      mode == MTensorPomeronMode::ResonanceContinuum) {
    ValidateTensorContinuumModel(lts);
  }
  if (lts.process.root_decay_mode == RootDecayMode::Isolated) { return {DecayType::None}; }
  if ((mode == MTensorPomeronMode::Resonance || mode == MTensorPomeronMode::Generic) &&
      !lts.process.RESONANCES.empty() &&
      std::all_of(lts.process.RESONANCES.begin(), lts.process.RESONANCES.end(), [](const auto &entry) {
        return ClassifyTensorResonance(entry.second.p) == TensorResonanceType::AxialVector;
      })) {
    return decay::JacobWickStructure(lts);
  }
  switch (mode) {
  case MTensorPomeronMode::Generic:
  case MTensorPomeronMode::Resonance:
  case MTensorPomeronMode::Continuum:
  case MTensorPomeronMode::ResonanceContinuum:
  case MTensorPomeronMode::Photo:
  case MTensorPomeronMode::QED:
    return MTensorPomeron::DirectDecayStructure();
  }
  throw std::invalid_argument(
      "MTensorPomeron::ProcessDefinitionFor: unknown process mode");
}

} // namespace gra::tensor
