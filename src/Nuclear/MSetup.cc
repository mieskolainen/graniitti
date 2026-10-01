// Nuclear UPC model construction
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Nuclear/MSetup.h"

// C++
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <iterator>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Own
#include "Graniitti/Math/MIntegration.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Nuclear/MConfig.h"
#include "Graniitti/Nuclear/MNucleus.h"
#include "Graniitti/Particle/MPDG.h"
#include "Graniitti/Tech/MAux.h"

namespace gra::nuclear {

using gra::aux::indices;

namespace {

// Compute the full unitarity inelastic pp probability at one impact parameter
// Abar = (A_00 + A_11 + A_22 + A_33) / 4
// S_el = I + i A_el, P_inel = 1 - Tr[S_el^dagger S_el]/4
double InelasticProb(const MEikonalMatrix &runtime, const double b, const double unitarity_margin) {
  const auto                 amplitude   = runtime.ImpactHelicityMatrix(b, 0, 0);
  const std::complex<double> forward     = 0.25 * (amplitude[0] + amplitude[5] + amplitude[10] + amplitude[15]);
  const double               total       = 2.0 * std::imag(forward);
  const double               elastic     = 0.25 * gra::SquaredNorm(amplitude);
  const double               probability = total - elastic;
  if (!std::isfinite(probability) || probability < -unitarity_margin || probability > 1.0 + unitarity_margin) {
    throw std::runtime_error("InelasticProb: probability " + std::to_string(probability) +
                             " is not physical at b = " + std::to_string(b));
  }
  return std::clamp(probability, 0.0, 1.0);
}

}  // namespace

// Derive the elementary inelastic profile from one pp eikonal table
// sigma_inel = 2 pi int b db P_inel(b)
NNProfile BuildNNProfile(const MEikonal &eikonal, const std::optional<double> sigma_nn) {
  if (!eikonal.IsInitialized() || eikonal.InitialState().size() != 2 ||
      std::abs(eikonal.InitialState()[0].pdg) != PDG::PDG_p || std::abs(eikonal.InitialState()[1].pdg) != PDG::PDG_p) {
    throw std::invalid_argument("BuildNNProfile: initialized pp eikonal required");
  }
  const auto &runtime = eikonal.GetMatrixRuntime();
  const auto &impact  = runtime.ImpactParameterNodes();
  if (impact.size() < 3) { throw std::invalid_argument("BuildNNProfile: pp impact grid is incomplete"); }

  double sigma_tot  = 0.0;
  double sigma_el   = 0.0;
  double sigma_inel = 0.0;
  eikonal.GetTotXS(sigma_tot, sigma_el, sigma_inel);
  if (!std::isfinite(sigma_tot) || !std::isfinite(sigma_el) || !std::isfinite(sigma_inel) || !(sigma_tot > 0.0) ||
      sigma_el < 0.0 || !(sigma_inel > 0.0) || std::abs(sigma_tot - sigma_el - sigma_inel) > 1.0e-10 * sigma_tot) {
    throw std::runtime_error("BuildNNProfile: invalid pp eikonal cross sections");
  }
  const double unitarity_tolerance = runtime.UnitarityTolerance();
  if (!std::isfinite(unitarity_tolerance) || !(unitarity_tolerance > 0.0)) {
    throw std::runtime_error("BuildNNProfile: invalid pp unitarity tolerance");
  }
  // sigma_max(S_el) <= 1 + eps gives P_inel >= -(2 eps + eps^2)
  const double unitarity_margin = unitarity_tolerance * (2.0 + unitarity_tolerance);

  NNProfile profile;
  profile.sigma = 1.0e3 * sigma_inel;
  profile.b_node.resize(impact.size(), 0.0);
  profile.inelastic.resize(impact.size(), 0.0);
  profile.fingerprint = runtime.RuntimeFingerprint() + ";NN_PROFILE_V6_ONE_PROJECTILE";
  for (const auto &i : indices(impact)) {
    profile.b_node[i]    = impact[i] * PDG::GeV2fm;
    profile.inelastic[i] = InelasticProb(runtime, impact[i], unitarity_margin);
  }

  const auto last = std::find_if(profile.inelastic.rbegin(), profile.inelastic.rend(),
                                 [](const double probability) { return probability > 0.0; });
  if (last == profile.inelastic.rend()) { throw std::runtime_error("BuildNNProfile: pp inelastic profile is empty"); }
  const std::size_t last_nonzero =
      profile.inelastic.size() - 1 - static_cast<std::size_t>(std::distance(profile.inelastic.rbegin(), last));
  const std::size_t support = std::min(profile.inelastic.size(), last_nonzero + 2);
  profile.b_node.resize(support);
  profile.inelastic.resize(support);

  const double eikonal_sigma = profile.sigma;
  const double table_sigma   = 10.0 * 2.0 * math::PI * math::LinearRadialIntegral(profile.b_node, profile.inelastic);
  if (!std::isfinite(table_sigma) || !(table_sigma > 0.0) ||
      std::abs(table_sigma - eikonal_sigma) > 0.01 * eikonal_sigma) {
    throw std::runtime_error("BuildNNProfile: pp profile does not reproduce sigma_inel");
  }
  profile.sigma = sigma_nn.value_or(eikonal_sigma);
  if (!std::isfinite(profile.sigma) || !(profile.sigma > 0.0)) {
    throw std::invalid_argument("BuildNNProfile: sigma_NN override must be finite and positive");
  }
  // Preserve P_inel while normalizing its tabulated integral to sigma_NN
  const double radius_scale = std::sqrt(profile.sigma / table_sigma);
  for (auto &b : profile.b_node) { b *= radius_scale; }

  profile.omega = ForwardGW(runtime).omega;
  return profile;
}

// Build one validated UPC runtime from physical beam states
// gamma_rel = p_1.p_2 / (m_1 m_2)
std::shared_ptr<const MUPC> BuildUPC(const std::array<MParticle, 2> &beam, const std::array<M4Vec, 2> &momentum,
                                     const int excitation, const UPCParam &param, const UPCMode mode) {
  std::array<BeamType, 2>                        type;
  std::array<std::shared_ptr<const MNucleus>, 2> nucleus;
  bool                                           has_nucleus = false;
  for (const auto &leg : indices(beam)) {
    if (IsNuclearPDG(beam[leg].pdg)) {
      if (beam[leg].spinX2 < 0) { throw std::invalid_argument("BuildUPC: nuclear beams require a physical spin"); }
      // Spin-independent nuclear currents have Tr(I)/(2J+1) = 1 for unpolarized ions
      // Magnetic and spin-dependent nuclear multipoles are omitted
      type[leg] = BeamType::Nucleus;
      if (leg > 0 && nucleus[0] != nullptr && beam[leg].pdg == beam[0].pdg) {
        nucleus[leg] = nucleus[0];
      } else {
        nucleus[leg] =
            std::make_shared<const MNucleus>(DefaultNucleusParam(beam[leg].pdg, beam[leg].mass, param.geometry));
      }
      has_nucleus = true;
    } else if (std::abs(beam[leg].pdg) == PDG::PDG_p) {
      type[leg] = BeamType::Proton;
    } else {
      const int pdg = std::abs(beam[leg].pdg);
      if (pdg != 11 && pdg != 13 && pdg != 15) {
        throw std::invalid_argument(
            "BuildUPC: UPC supports proton, charged-lepton and nuclear "
            "beams");
      }
      type[leg] = BeamType::Lepton;
    }
  }
  if (!has_nucleus) { throw std::invalid_argument("BuildUPC: UPC requires a nuclear beam"); }
  if (excitation != 0 && !(excitation == 1 && (type[0] == BeamType::Proton || type[1] == BeamType::Proton))) {
    throw std::invalid_argument("BuildUPC: NSTARS requires single proton excitation in pA collisions");
  }

  UPCParam     configured   = param;
  const double mass_product = beam[0].mass * beam[1].mass;
  if (!(mass_product > 0.0)) { throw std::invalid_argument("BuildUPC: beam masses must be positive"); }
  const double gamma = (momentum[0] * momentum[1]) / mass_product;
  if (!std::isfinite(gamma) || !(gamma > 1.0)) {
    throw std::invalid_argument("BuildUPC: invalid relative beam Lorentz factor");
  }
  // Additional absorption supplies excitation to the explicit remnant decay
  const bool realize_emd = configured.additional_emd;
  for (const auto &leg : indices(nucleus)) {
    if (nucleus[leg] == nullptr || !realize_emd) {
      configured.breakup[leg] = {};
      continue;
    }
    const std::size_t emitter = 1 - leg;
    const auto response = std::find_if(configured.emd.isotope.begin(), configured.emd.isotope.end(),
        [&](const auto &item) { return item.a == nucleus[leg]->A() && item.z == nucleus[leg]->Z(); });
    if (response == configured.emd.isotope.end()) {
      throw std::invalid_argument("BuildUPC: no EMD response for the target isotope");
    }
    configured.breakup[leg]   = *response;
    auto &breakup             = configured.breakup[leg];
    breakup.photo.gdr         = configured.emd.gdr;
    breakup.z_emit            = std::abs(static_cast<double>(beam[emitter].chargeX3)) / 3.0;
    breakup.gamma             = gamma;
    if (nucleus[emitter] != nullptr) {
      breakup.emitter = EmitterType::Nuclear;
    } else if (type[emitter] == BeamType::Proton) {
      breakup.emitter = EmitterType::Proton;
    } else {
      breakup.emitter = EmitterType::Point;
    }
  }

  return std::make_shared<const MUPC>(type, std::move(nucleus), configured,
                                      std::array<std::shared_ptr<const MConfigBank>, 2>{}, mode);
}

// Print the active nuclear UPC setup fields
void PrintUPCSetup(const MUPC &upc, std::ostream &stream) {
  const UPCParam &param    = upc.Param();
  const bool      emission = upc.Photon(1) != nullptr || upc.Photon(2) != nullptr;
  const bool      photo    = upc.Photo(1) != nullptr || upc.Photo(2) != nullptr;

  if (emission) {
    stream << "- Nuclear emission:       [" << CoherenceName(param.emission[0]) << ", "
           << CoherenceName(param.emission[1]) << ']' << std::endl;
  }
  if (photo) {
    stream << "- Nuclear photo target:   [" << CoherenceName(param.target[0]) << ", " << CoherenceName(param.target[1])
           << ']' << std::endl;
    stream << "- Nuclear photo model:     '" << PhotoModelName(param.photo_model) << '\'' << std::endl;
  }

  stream << "- Nuclear screening:      " << std::boolalpha << upc.HadronicConvolution() << " ['"
         << SurvivalName(param.survival) << "']" << std::endl;
  if (param.glauber.profile.sigma > 0.0) {
    stream << "- Nuclear sigma_NN:       " << param.glauber.profile.sigma << " mb ["
           << (param.sigma_nn.has_value() ? "gencard override" : "eikonal") << ']' << std::endl;
  } else {
    stream << "- Nuclear sigma_NN:       inactive" << std::endl;
  }
  if (!param.survival_eikonal.empty()) {
    stream << "- Nuclear NN eikonal:     '" << param.survival_eikonal << "', sqrt(s_NN) = " << std::sqrt(param.s_nn)
           << " GeV, omega_N = " << param.glauber.profile.omega << std::endl;
  }
  stream << "- Nuclear additional EMD: " << std::boolalpha << param.additional_emd << std::endl;
  stream << "- Nuclear neutron class:  [" << NeutronName(param.neutron[0]) << ", " << NeutronName(param.neutron[1])
         << ']' << std::endl;
  stream << "- Nuclear structure:       '" << StructureName(param.structure) << '\'' << std::endl;
}

}  // namespace gra::nuclear
