// Nuclear UPC steering and screening loops
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Nuclear/MUPC.h"

#include <algorithm>
#include <cmath>
#include <charconv>
#include <regex>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <utility>

#include "Graniitti/Math/MAlgorithms.h"
#include "Graniitti/Math/MIntegration.h"
#include "Graniitti/Math/MInterpolation.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MPolarFourier.h"
#include "Graniitti/Math/MStatistics.h"
#include "Graniitti/Nuclear/MCache.h"
#include "Graniitti/Particle/MForm.h"
#include "Graniitti/Particle/MPDG.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MException.h"
#include "Graniitti/Tech/MJsonZip.h"
#include "json.hpp"

using gra::aux::indices;
using gra::math::ParallelFor;

namespace gra::nuclear {
namespace {

// Transform the finite-range configuration profiles with the shared harmonic quadrature
std::vector<std::complex<double>> ConfigTransform(const ConfigGrid& grid,
                                                  const std::vector<std::vector<std::complex<double>>>& profile) {
  const auto modes = grid.b_unit.size(), samples = profile.size(), radii = grid.b_node.size();
  std::vector<std::complex<double>> result(grid.k_node.size() * samples * modes, 0.0);
  for (const auto& ik : indices(grid.k_node)) {
    for (const auto& sample : indices(profile)) {
      for (std::size_t mode = 0; mode < modes; ++mode) {
        for (const auto& ib : indices(grid.b_node)) {
          result[(ik * samples + sample) * modes + mode] +=
              grid.bessel[(ik * modes + mode) * radii + ib] * profile[sample][ib * modes + mode];
        }
      }
    }
  }
  return result;
}

// Compute one checked zero-based two-leg index
std::size_t LegIndex(const int leg, const std::string &source) {
  if (leg < 1 || leg > 2) { throw std::out_of_range(source + ": leg must be one or two"); }
  return static_cast<std::size_t>(leg - 1);
}

// Construct an open logarithmic rule for the polar measure r dr
std::vector<std::pair<double, double>> LogRadialRule(const unsigned int count, const double r_min, const double r_max,
                                                     const double scale) {
  if (count == 0 || !(r_max > r_min) || !(scale > 0.0)) { return {}; }
  const auto                             rule   = math::GaussLegendreRule(count, 0.0, 1.0);
  const double                           scale2 = scale * scale;
  const double                           base   = scale2 + r_min * r_min;
  const double                           range  = std::log((scale2 + r_max * r_max) / base);
  std::vector<std::pair<double, double>> output;
  output.reserve(count);
  for (const auto &i : indices(rule.first)) {
    const double radius2 = base * std::exp(rule.first[i] * range) - scale2;
    const double radius  = std::sqrt(std::max(0.0, radius2));
    const double measure = 0.5 * rule.second[i] * range * (scale2 + radius2);
    output.emplace_back(radius, measure);
  }
  return output;
}

// Compute the allowed radial intervals of one ray inside a centered annulus
std::vector<std::pair<double, double>> DiskRayIntervals(const double cx, const double cy, const double angle,
                                                        const double r_min, const double r_max) {
  const double cosine     = std::cos(angle);
  const double sine       = std::sin(angle);
  const double center2    = cx * cx + cy * cy;
  const double projection = cx * cosine + cy * sine;
  const double outer_disc = projection * projection + r_max * r_max - center2;
  if (!(outer_disc > 0.0)) { return {}; }
  const double outer = -projection + std::sqrt(outer_disc);
  if (!(outer > 0.0) || !(r_min > 0.0)) {
    return outer > 0.0 ? std::vector<std::pair<double, double>>{{0.0, outer}}
                       : std::vector<std::pair<double, double>>{};
  }

  const double inner_disc = projection * projection + r_min * r_min - center2;
  if (!(inner_disc > 0.0)) { return {{0.0, outer}}; }
  const double root  = std::sqrt(inner_disc);
  const double lower = std::max(0.0, -projection - root);
  const double upper = std::min(outer, -projection + root);
  if (!(upper > lower)) { return {{0.0, outer}}; }
  std::vector<std::pair<double, double>> interval;
  if (lower > 0.0) { interval.emplace_back(0.0, lower); }
  if (outer > upper) { interval.emplace_back(upper, outer); }
  return interval;
}

// Allocate a fixed radial node count between at most two ray intervals
std::vector<unsigned int> RadialCounts(const std::vector<std::pair<double, double>> &interval,
                                       const unsigned int                            count) {
  if (interval.empty() || count < interval.size()) { return {}; }
  std::vector<unsigned int> output(interval.size(), 1);
  unsigned int              remaining = count - static_cast<unsigned int>(interval.size());
  if (interval.size() == 1) {
    output[0] += remaining;
    return output;
  }
  const double       area0    = math::pow2(interval[0].second) - math::pow2(interval[0].first);
  const double       area1    = math::pow2(interval[1].second) - math::pow2(interval[1].first);
  const double       fraction = area0 / (area0 + area1);
  const unsigned int first    = static_cast<unsigned int>(std::llround(fraction * static_cast<double>(remaining)));
  output[0] += std::min(first, remaining);
  output[1] += remaining - std::min(first, remaining);
  return output;
}

// Compute the current persistent UPC table representation version
constexpr std::size_t CacheVersion() { return 2; }

// Append one floating-point value in an exact cache-key representation
void AppendDouble(std::ostringstream &stream, const double value) { stream << std::hexfloat << value << ';'; }

// Append one density definition to a cache key
void AppendDensity(std::ostringstream &stream, const DensityParam &param) {
  AppendDouble(stream, param.radius);
  AppendDouble(stream, param.skin);
  AppendDouble(stream, param.r_max);
  stream << param.nodes << ';' << param.cdf_nodes << ';';
  AppendDouble(stream, param.form_q_max);
  stream << param.form_nodes << ';';
  AppendDouble(stream, param.form_abs_tol);
}

// Serialize the effective convolution and Glauber parameters independently of input parsing
nlohmann::json ConvolutionKey(const UPCParam &p) {
  const auto &l = p.loop;
  const auto &c = p.convolution;
  const auto &g = p.glauber;
  const auto &f = g.fluctuation;
  return {"one-projectile-omega-v1", {l.radial_integrator, l.azimuth_integrator, static_cast<unsigned int>(l.radial_map), l.r_min, l.r_max,
           l.radial_intervals, l.azimuth_nodes},
          {c.b_max, c.smooth_b_nodes, c.sample_b_nodes, c.b_phi_nodes, c.smooth_kt_nodes, c.sample_kt_nodes},
          {g.profile.sigma, g.profile.omega, g.profile.b_node, g.profile.inelastic, g.b_max, g.q_max, g.b_nodes,
           g.q_nodes, g.profile_nodes},
          {f.omega, f.nodes, f.cdf_tol, f.cdf_max_iter, f.normal_shape_min},
          {p.structure_param.EM, p.structure_param.F2, p.structure_param.QED_alpha}};
}

// Compute the persistent filename for one complete UPC cache key
std::string CacheFilename(const std::string &key, const std::array<BeamType, 2> &type,
                          const std::array<std::shared_ptr<const MNucleus>, 2> &nucleus) {
  std::ostringstream name;
  name << "UPC";
  for (const auto &leg : indices(type)) {
    name << '_' << BeamName(type[leg]);
    if (nucleus[leg] != nullptr) { name << nucleus[leg]->ID().pdg; }
  }
  return NuclearCacheFilename(name.str(), key);
}

// Resolve the mean and normalized fluctuation components of one nuclear current
// R_coh,c = 1, R_inc,c = (J_c - <J>)/sqrt[Var(J)]
// R_incl,c = J_c/sqrt[<|J|^2>]
std::vector<std::complex<double>> CurrentRatios(const std::vector<std::complex<double>> &current,
                                                const CoherenceType type) {
  if (current.empty()) { throw std::logic_error("MUPC: empty current ensemble"); }
  if (type == CoherenceType::Coherent) { return std::vector<std::complex<double>>(current.size(), 1.0); }
  const auto stat = statistics::ComplexMoments(current);

  const bool   incoherent = type == CoherenceType::Incoherent;
  const double second     = incoherent ? stat.variance : stat.second;
  const auto   center     = incoherent ? stat.mean : std::complex<double>{};
  if (!std::isfinite(second)) { throw AmplitudeFailure("MUPC: invalid current normalization"); }
  std::vector<std::complex<double>> ratio(current.size(), 0.0);
  // Normalize every representable positive moment, independent of the absolute current scale
  if (!(second > 0.0)) { return ratio; }
  const double norm = std::sqrt(second);
  for (const auto &sample : indices(ratio)) { ratio[sample] = (current[sample] - center) / norm; }
  return ratio;
}

// Parse one exact steering selector
template <typename T, std::size_t N>
T Parse(const std::string &name, const std::array<std::string_view, N> &value, const std::string_view source) {
  for (const auto &i : indices(value)) {
    if (name == value[i]) { return static_cast<T>(i); }
  }
  throw std::invalid_argument(std::string(source) + ": unsupported value " + name);
}

// Compute the exact steering name of one selector
template <typename T, std::size_t N>
std::string Name(const T type, const std::array<std::string_view, N> &value, const std::string_view source) {
  const auto index = static_cast<std::size_t>(type);
  if (index >= value.size()) { throw std::logic_error(std::string(source) + ": invalid selector"); }
  return std::string(value[index]);
}

constexpr std::array<std::string_view, 3> beam_name      = {"lepton", "proton", "nucleus"};
constexpr std::array<std::string_view, 3> coherence_name = {"coherent", "incoherent", "inclusive"};
constexpr std::array<std::string_view, 7> neutron_op = {"*", "==", "!=", ">", ">=", "<", "<="};
constexpr std::array<std::string_view, 3> structure_name = {"smooth", "nucleon", "hotspot"};
constexpr std::array<std::string_view, 3> survival_name  = {"optical", "optical_ggcf", "mc_ggcf"};
constexpr std::array<std::string_view, 3> photo_name     = {"impulse", "glauber", "LTA"};

}  // namespace

// Parse one lowercase beam name
BeamType ParseBeamType(const std::string &name) { return Parse<BeamType>(name, beam_name, "ParseBeamType"); }

// Parse one lowercase nuclear-sector name
CoherenceType ParseCoherenceType(const std::string &name) {
  return Parse<CoherenceType>(name, coherence_name, "ParseCoherenceType");
}

// Check an actual emitted neutron count without parsing in the event loop
bool NeutronSelection::Accept(const std::size_t n) const {
  switch (type) {
    case Type::Any:          return true;
    case Type::Equal:        return n == count;
    case Type::NotEqual:     return n != count;
    case Type::Greater:      return n > count;
    case Type::GreaterEqual: return n >= count;
    case Type::Less:         return n < count;
    case Type::LessEqual:    return n <= count;
  }
  return false;
}

// Parse a wildcard or one nonnegative neutron-count comparison at initialization
NeutronSelection ParseNeutronSelection(const std::string &name) {
  const std::regex grammar(R"(^\s*(\*|n\s*(==|=|!=|>=|<=|>|<)\s*([0-9]+))\s*$)");
  std::smatch match;
  if (!std::regex_match(name, match, grammar)) {
    throw std::invalid_argument("ParseNeutronSelection: expected * or n <operator> <nonnegative integer>: " + name);
  }
  if (match[1] == "*") { return {}; }
  NeutronSelection result;
  result.type = Parse<NeutronSelection::Type>(match[2] == "=" ? "==" : match[2].str(), neutron_op, "ParseNeutronSelection");
  const auto digits = match[3].str();
  const auto parsed = std::from_chars(digits.data(), digits.data() + digits.size(), result.count);
  if (parsed.ec != std::errc{}) { throw std::invalid_argument("ParseNeutronSelection: neutron count overflow"); }
  return result;
}

// Parse one lowercase nuclear-structure name
StructureType ParseStructureType(const std::string &name) {
  return Parse<StructureType>(name, structure_name, "ParseStructureType");
}

// Parse one lowercase survival-model name
SurvivalType ParseSurvivalType(const std::string &name) {
  return Parse<SurvivalType>(name, survival_name, "ParseSurvivalType");
}

// Parse one lowercase photonuclear model name
PhotoModel ParsePhotoModel(const std::string &name) { return Parse<PhotoModel>(name, photo_name, "ParsePhotoModel"); }

// Compute one lowercase beam name
std::string BeamName(const BeamType type) { return Name(type, beam_name, "BeamName"); }

// Compute one lowercase nuclear-sector name
std::string CoherenceName(const CoherenceType type) { return Name(type, coherence_name, "CoherenceName"); }

// Compute the canonical neutron-count expression
std::string NeutronName(const NeutronSelection selection) {
  if (selection.type == NeutronSelection::Type::Any) { return "*"; }
  return "n " + Name(selection.type, neutron_op, "NeutronName") + " " + std::to_string(selection.count);
}

// Compute one lowercase nuclear-structure name
std::string StructureName(const StructureType type) { return Name(type, structure_name, "StructureName"); }

// Compute one lowercase survival-model name
std::string SurvivalName(const SurvivalType type) { return Name(type, survival_name, "SurvivalName"); }

// Compute one lowercase photonuclear model name
std::string PhotoModelName(const PhotoModel model) { return Name(model, photo_name, "PhotoModelName"); }

// Construct and precompute one heterogeneous UPC convolution
MUPC::MUPC(std::array<BeamType, 2> type, std::array<std::shared_ptr<const MNucleus>, 2> nucleus, UPCParam param,
           std::array<std::shared_ptr<const MConfigBank>, 2> bank, const UPCMode mode)
    : param_(param),
      type_(type),
      nucleus_(std::move(nucleus)),
      breakup_{std::make_shared<const MBreakup>(param_.breakup[0], param_.breakup[0].a > 0 ? nucleus_[1] : nullptr,
                                                param_.structure_param),
               std::make_shared<const MBreakup>(param_.breakup[1], param_.breakup[1].a > 0 ? nucleus_[0] : nullptr,
                                                param_.structure_param)},
      bank_(std::move(bank)),
      screening_(mode.screening) {
  param_.glauber.fluctuation.omega = param_.glauber.profile.omega;
  for (const auto &i : indices(type_)) {
    if ((type_[i] == BeamType::Nucleus) != (nucleus_[i] != nullptr)) {
      throw std::invalid_argument("MUPC: beam type and nuclear identity disagree");
    }
    if (type_[i] == BeamType::Nucleus) {
      if (mode.emission[i]) { photon_[i] = std::make_shared<const MPhoton>(nucleus_[i], param_.structure_param); }
      if (mode.target[i]) {
        photo_[i] = std::make_shared<const MPhoto>(nucleus_[i], param_.photo[i], param_.photo_model);
      }
      const bool config_survival = HadronicConvolution() && param_.survival == SurvivalType::MCGGCF;
      if (param_.config.count > 0 && (mode.current[i] || config_survival)) {
        sampler_[i] = std::make_shared<const MConfigSampler>(*nucleus_[i], param_.config);
      }
    }
  }
  if (HadronicConvolution()) {
    std::array<HadronType, 2> hadron;
    for (const auto &i : indices(type_)) {
      hadron[i] = type_[i] == BeamType::Proton ? HadronType::Proton : HadronType::Nucleus;
    }
    glauber_ = std::make_shared<const MGlauber>(hadron, nucleus_, param_.glauber);
  }
  Validate();
  Prepare();
}

// Construct and precompute one nuclear screening convolution
MUPC::MUPC(const MNucleus &beam1, const MNucleus &beam2, UPCParam param, std::shared_ptr<const MConfigBank> bank1,
           std::shared_ptr<const MConfigBank> bank2, const UPCMode mode)
    : MUPC({BeamType::Nucleus, BeamType::Nucleus},
           {std::make_shared<const MNucleus>(beam1), std::make_shared<const MNucleus>(beam2)}, param,
           {std::move(bank1), std::move(bank2)}, mode) {}

// Construct one event-local runtime sharing all deterministic model tables
MUPC::MUPC(const MUPC &model, std::array<std::shared_ptr<const MConfigBank>, 2> bank,
           std::array<std::optional<MHotSpot>, 2> hotspot, const bool prepare_convolution)
    : param_(model.param_),
      type_(model.type_),
      nucleus_(model.nucleus_),
      glauber_(model.glauber_),
      breakup_(model.breakup_),
      bank_(std::move(bank)),
      sampler_(model.sampler_),
      hotspot_(std::move(hotspot)),
      screening_(model.screening_),
      config_grid_(model.config_grid_),
      photon_(model.photon_),
      photo_(model.photo_) {
  ValidateBanks();
  const bool config_profile = glauber_ != nullptr && param_.survival == SurvivalType::MCGGCF;
  if (config_profile && prepare_convolution) {
    Prepare();
  } else {
    PrepareSamples();
    b_node_        = model.b_node_;
    b_weight_      = model.b_weight_;
    profile_       = model.profile_;
    loop_k_        = model.loop_k_;
    loop_profile_  = model.loop_profile_;
    loop_harmonic_ = model.loop_harmonic_;
  }
}

// Validate all UPC controls and optional configuration banks
void MUPC::Validate() const {
  if (!std::isfinite(param_.emd_focus) || param_.emd_focus < 0.0 || !(param_.emd_focus < 1.0)) {
    throw std::invalid_argument("MUPC: EMD focus fraction must be in [0,1)");
  }
  if (!param_.reaction && (param_.additional_emd || std::any_of(param_.neutron.begin(), param_.neutron.end(),
                                                                [](auto n) { return n.type != NeutronSelection::Type::Any; }))) {
    throw std::invalid_argument("MUPC: EMD and neutron classes require internal or external fragmentation");
  }
  if (param_.sigma_nn.has_value() &&
      (!HadronicConvolution() || !std::isfinite(*param_.sigma_nn) || *param_.sigma_nn <= 0.0)) {
    throw std::invalid_argument("MUPC: sigma_NN override requires active hadronic survival");
  }
  if (param_.structure == StructureType::Smooth) {
    std::array<bool, 2> coherent{};
    std::array<bool, 2> incoherent{};
    bool                both_directions = true;
    for (const auto &leg : indices(type_)) {
      const bool nuclear = type_[leg] == BeamType::Nucleus;
      both_directions &= (type_[leg] == BeamType::Proton || (nuclear && photo_[leg] != nullptr)) &&
                         (!nuclear || photon_[leg] != nullptr);
      coherent[leg] = !nuclear || (param_.emission[leg] != CoherenceType::Incoherent &&
                                   param_.target[leg] != CoherenceType::Incoherent);
      incoherent[leg] =
          nuclear && param_.emission[leg] != CoherenceType::Coherent && param_.target[leg] != CoherenceType::Coherent;
    }
    if (both_directions &&
        ((incoherent[0] && (coherent[1] || incoherent[1])) || (incoherent[1] && (coherent[0] || incoherent[0])))) {
      throw std::invalid_argument("MUPC: interfering incoherent photon directions require sampled nuclear structure");
    }
  }
  const bool sampled_structure = param_.structure != StructureType::Smooth;
  if (sampled_structure != (param_.config.count > 0) || (!sampled_structure && param_.current_count > 0)) {
    throw std::invalid_argument("MUPC: structure and sampling disagree");
  }
  if (param_.config.count == 1) {
    throw std::invalid_argument("MUPC: Good-Walker sampling requires at least two nuclear configurations");
  }
  if (!std::isfinite(param_.loop.r_min) || !std::isfinite(param_.loop.r_max) || param_.loop.r_min < 0.0 ||
      !(param_.loop.r_max > param_.loop.r_min) ||
      (param_.loop.radial_map == math::RadialMap::Log && !(param_.loop.r_min > 0.0))) {
    throw std::invalid_argument("MUPC: invalid screening momentum range");
  }
  const bool convolution = HadronicConvolution() || param_.additional_emd;
  if (convolution && (param_.loop.radial_integrator != "GL" || param_.loop.azimuth_integrator != "Trap")) {
    throw std::invalid_argument(
        "MUPC: focused screening requires GL radial and Trap azimuthal "
        "rules");
  }
  if (convolution && (param_.loop.radial_intervals < 2 || param_.loop.azimuth_nodes == 0)) {
    throw std::invalid_argument("MUPC: focused screening requires at least two radial nodes and a nonempty azimuth");
  }
  if (convolution && (!std::isfinite(param_.convolution.b_max) || !(param_.convolution.b_max > 0.0) ||
                                param_.convolution.smooth_b_nodes == 0 || param_.convolution.sample_b_nodes == 0 ||
                                param_.convolution.b_phi_nodes == 0 || param_.convolution.b_phi_nodes % 2 == 0 ||
                                param_.convolution.smooth_kt_nodes == 0 || param_.convolution.sample_kt_nodes == 0)) {
    throw std::invalid_argument("MUPC: invalid impact screening rule");
  }
  if ((param_.current_count == 1) ||
      (param_.current_count > 0 && (param_.config.count == 0 || param_.current_count > param_.config.count))) {
    throw std::invalid_argument("MUPC: invalid current configuration count");
  }
  for (const auto& i : indices(type_)) {
    if (sampler_[i] && (param_.structure == StructureType::Hotspot || param_.photo[i].hotspot.count > 0)) {
      param_.photo[i].hotspot.Validate();
    }
    if (type_[i] != BeamType::Nucleus &&
        (param_.emission[i] != CoherenceType::Coherent || param_.target[i] != CoherenceType::Coherent)) {
      throw std::invalid_argument("MUPC: incoherent nuclear modes require a nuclear beam");
    }
    if (type_[i] != BeamType::Nucleus && param_.neutron[i].type != NeutronSelection::Type::Any) {
      throw std::invalid_argument("MUPC: neutron classes apply only to nuclear beams");
    }
    const bool hard = (photon_[i] && param_.emission[i] != CoherenceType::Coherent) ||
                      (photo_[i] && param_.target[i] != CoherenceType::Coherent);
    if (!param_.neutron[i].Accept(0) && !breakup_[i]->Enabled() && !hard) {
      throw std::invalid_argument("MUPC: xn requires EMD or an incoherent hard nuclear transition");
    }
  }
  if (HadronicConvolution() && param_.survival == SurvivalType::MCGGCF && param_.config.count == 0 &&
      std::any_of(type_.begin(), type_.end(), [](BeamType type) { return type == BeamType::Nucleus; })) {
    throw std::invalid_argument("MUPC: configuration survival requires sampled configurations");
  }
  ValidateBanks();
}

// Validate only the new event configuration banks against the immutable beam model
void MUPC::ValidateBanks() const {
  const bool paired  = HadronicConvolution() && param_.survival == SurvivalType::MCGGCF;
  const bool sampled = bank_[0] != nullptr || bank_[1] != nullptr;
  for (const auto &i : indices(type_)) {
    if (bank_[i] && (!nucleus_[i] || bank_[i]->Nucleus().ID().pdg != nucleus_[i]->ID().pdg)) {
      throw std::invalid_argument("MUPC: configuration bank nucleus does not match its beam");
    }
    if (hotspot_[i] && (!bank_[i] || hotspot_[i]->ConfigCount() != bank_[i]->Size())) {
      throw std::invalid_argument("MUPC: hotspot and nucleon banks disagree");
    }
    if (paired && sampled && nucleus_[i] && !bank_[i]) {
      throw std::invalid_argument("MUPC: event-local configuration ensemble is incomplete");
    }
  }
}

// Prepare the complete Cartesian configuration ensemble
void MUPC::PrepareSamples() {
  bool sampled = false;
  for (const auto &i : indices(type_)) {
    if (type_[i] == BeamType::Nucleus && bank_[i] != nullptr) {
      sample_shape_[i] = bank_[i]->Size();
      sampled          = true;
    }
  }
  if (!sampled) { return; }
  if (sample_shape_[0] > std::numeric_limits<std::size_t>::max() / sample_shape_[1]) {
    throw std::invalid_argument("MUPC: configuration ensemble is too large");
  }
  sample_count_ = sample_shape_[0] * sample_shape_[1];
}

// Construct the full physics and numerical disk-cache key
std::string MUPC::CacheKey() const {
  std::ostringstream stream;
  stream << CacheVersion() << ';' << param_.table_fingerprint << ';' << param_.structure_param.EM << ';';
  stream << ConvolutionKey(param_).dump() << ';';
  for (const auto &leg : indices(type_)) {
    stream << static_cast<unsigned int>(type_[leg]) << ';' << static_cast<unsigned int>(param_.emission[leg]) << ';'
           << static_cast<unsigned int>(param_.target[leg]) << ';';
    if (nucleus_[leg] != nullptr) {
      const auto &nucleus = nucleus_[leg]->Param();
      stream << nucleus.pdg << ';';
      AppendDouble(stream, nucleus.mass);
      AppendDensity(stream, nucleus.charge);
      AppendDensity(stream, nucleus.matter);
    } else {
      stream << 0 << ';';
    }
  }
  stream << static_cast<unsigned int>(screening_) << ';' << static_cast<unsigned int>(param_.survival) << ';';
  if (HadronicConvolution()) {
    stream << param_.glauber.profile.fingerprint << ';';
    AppendDouble(stream, param_.glauber.profile.sigma);
    AppendDouble(stream, param_.glauber.profile.omega);
  }
  return stream.str();
}

// Load and validate one compressed JSON screening table
bool MUPC::ReadCache(const std::string &filename, const std::string &key) {
  try {
    nlohmann::json cache;
    if (!ReadNuclearCache(filename, CacheVersion(), "scalar", key, cache)) { return false; }
    const std::size_t b_count         = param_.convolution.smooth_b_nodes;
    auto              transform_loop  = param_.loop;
    transform_loop.radial_intervals   = std::max(param_.convolution.smooth_kt_nodes, param_.loop.radial_intervals);
    const std::size_t transform_count = math::PolarNodeCount(transform_loop) + 2;
    if (!MJsonZip::DecompressVector(cache.at("b_node"), b_node_, b_count) ||
        !MJsonZip::DecompressVector(cache.at("b_weight"), b_weight_, b_count) ||
        !MJsonZip::DecompressVector(cache.at("profile"), profile_, b_count)) {
      return false;
    }
    if (!MJsonZip::DecompressVector(cache.at("loop_k"), loop_k_, transform_count) ||
        !MJsonZip::DecompressVector(cache.at("loop_profile"), loop_profile_, transform_count)) {
      return false;
    }
    math::ValidateInterpolationGrid(b_node_);
    math::ValidateInterpolationGrid(loop_k_);
    sample_profile_.clear();
    return true;
  } catch (const std::exception &) { return false; }
}

// Save one complete compressed JSON screening table
void MUPC::WriteCache(const std::string &filename, const std::string &key) const {
  auto cache = NuclearCache(CacheVersion(), "scalar", key);
  cache.update({{"b_node", MJsonZip::CompressVector(b_node_)},
                {"b_weight", MJsonZip::CompressVector(b_weight_)},
                {"profile", MJsonZip::CompressVector(profile_)},
                {"loop_k", MJsonZip::CompressVector(loop_k_)},
                {"loop_profile", MJsonZip::CompressVector(loop_profile_)}});
  PublishNuclearCache(filename, cache);
}

// Load one deterministic config-screening grid
bool MUPC::ReadConfigCache(const std::string &filename, const std::string &key) {
  try {
    nlohmann::json cache;
    if (!ReadNuclearCache(filename, CacheVersion(), "config_grid", key, cache)) { return false; }
    const std::size_t b_count        = param_.convolution.sample_b_nodes;
    const std::size_t b_phi_count    = param_.convolution.b_phi_nodes;
    auto              transform_loop = param_.loop;
    transform_loop.radial_intervals  = std::max(param_.convolution.sample_kt_nodes, param_.loop.radial_intervals);
    const std::size_t k_count        = math::PolarNodeCount(transform_loop) + 2;
    auto              grid           = std::make_shared<ConfigGrid>();
    if (!MJsonZip::DecompressVector(cache.at("b_node"), grid->b_node, b_count) ||
        !MJsonZip::DecompressVector(cache.at("b_weight"), grid->b_weight, b_count) ||
        !MJsonZip::DecompressComplexVector(cache.at("b_unit"), grid->b_unit, b_phi_count) ||
        !MJsonZip::DecompressComplexVector(cache.at("b_phase"), grid->b_phase, b_phi_count * b_phi_count) ||
        !MJsonZip::DecompressVector(cache.at("k_node"), grid->k_node, k_count) ||
        !MJsonZip::DecompressVector(cache.at("bessel"), grid->bessel, k_count * b_phi_count * b_count)) {
      return false;
    }
    math::ValidateInterpolationGrid(grid->b_node);
    math::ValidateInterpolationGrid(grid->k_node);
    config_grid_ = std::move(grid);
    return true;
  } catch (const std::exception &) { return false; }
}

// Save one deterministic config-screening grid
void MUPC::WriteConfigCache(const std::string &filename, const std::string &key) const {
  if (config_grid_ == nullptr) { throw std::logic_error("MUPC config cache has no deterministic grid"); }
  auto cache = NuclearCache(CacheVersion(), "config_grid", key);
  cache.update({{"b_node", MJsonZip::CompressVector(config_grid_->b_node)},
                {"b_weight", MJsonZip::CompressVector(config_grid_->b_weight)},
                {"b_unit", MJsonZip::CompressComplexVector(config_grid_->b_unit)},
                {"b_phase", MJsonZip::CompressComplexVector(config_grid_->b_phase)},
                {"k_node", MJsonZip::CompressVector(config_grid_->k_node)},
                {"bessel", MJsonZip::CompressVector(config_grid_->bessel)}});
  PublishNuclearCache(filename, cache);
}

// Prepare deterministic config-screening quadrature and transform tables
// F_m(k) = 2 pi / hbarc^2 int b db J_m(kb/hbarc) P_m(b)
void MUPC::PrepareConfigGrid() {
  auto grid      = std::make_shared<ConfigGrid>();
  auto b_rule    = math::GaussLegendreRule(param_.convolution.sample_b_nodes, 0.0, param_.convolution.b_max);
  grid->b_node   = std::move(b_rule.first);
  grid->b_weight = std::move(b_rule.second);

  const double inverse_b_phi = 1.0 / static_cast<double>(param_.convolution.b_phi_nodes);
  grid->b_unit.resize(param_.convolution.b_phi_nodes);
  grid->b_phase.resize(static_cast<std::size_t>(param_.convolution.b_phi_nodes) * param_.convolution.b_phi_nodes, 0.0);
  for (const auto &angle : indices(grid->b_unit)) {
    const double phi    = 2.0 * math::PI * (static_cast<double>(angle) + 0.5) * inverse_b_phi;
    grid->b_unit[angle] = std::polar(1.0, phi);
    for (unsigned int mode = 0; mode < param_.convolution.b_phi_nodes; ++mode) {
      grid->b_phase[angle * param_.convolution.b_phi_nodes + mode] =
          inverse_b_phi * std::polar(1.0, -math::SignedHarmonic(mode, param_.convolution.b_phi_nodes) * phi);
    }
  }

  auto transform_loop             = param_.loop;
  transform_loop.radial_intervals = std::max(param_.convolution.sample_kt_nodes, param_.loop.radial_intervals);
  auto kt_rule                    = math::PolarRadialRule(transform_loop);
  grid->k_node                    = std::move(kt_rule.node);
  // Add interpolation bounds without changing the fixed quadrature
  if (grid->k_node.front() > param_.loop.r_min) { grid->k_node.insert(grid->k_node.begin(), param_.loop.r_min); }
  if (grid->k_node.back() < param_.loop.r_max) { grid->k_node.push_back(param_.loop.r_max); }

  const std::size_t mode_count = param_.convolution.b_phi_nodes;
  const std::size_t b_count    = grid->b_node.size();
  grid->bessel.resize(grid->k_node.size() * mode_count * b_count, 0.0);
  constexpr double hbarc = PDG::GeV2fm;
  const double     scale = 2.0 * math::PI / (hbarc * hbarc);
  ParallelFor(
      grid->k_node.size() * mode_count,
      [&](const std::size_t task) {
        const std::size_t  ik       = task / mode_count;
        const unsigned int mode     = static_cast<unsigned int>(task % mode_count);
        const int          harmonic = math::SignedHarmonic(mode, param_.convolution.b_phi_nodes);
        const std::size_t  first    = task * b_count;
        for (const auto &ib : indices(grid->b_node)) {
          double bessel = std::cyl_bessel_j(std::abs(harmonic), grid->k_node[ik] * grid->b_node[ib] / hbarc);
          if (harmonic < 0 && std::abs(harmonic) % 2 != 0) { bessel = -bessel; }
          grid->bessel[first + ib] = scale * grid->b_weight[ib] * grid->b_node[ib] * bessel;
        }
      },
      "Nuclear UPC deterministic Bessel grid");
  config_grid_ = std::move(grid);
}

// Prepare the dense radial transform used by smooth focused quadrature
// F(k) = 2 pi / hbarc^2 int b db J_0(kb/hbarc) P(b)
void MUPC::PrepareScalarTransform() {
  auto transform_loop             = param_.loop;
  transform_loop.radial_intervals = std::max(param_.convolution.smooth_kt_nodes, param_.loop.radial_intervals);
  auto rule                       = math::PolarRadialRule(transform_loop);
  loop_k_                         = std::move(rule.node);
  // Add both bounds once so event-local interpolation never transforms again
  if (loop_k_.front() > param_.loop.r_min) { loop_k_.insert(loop_k_.begin(), param_.loop.r_min); }
  if (loop_k_.back() < param_.loop.r_max) { loop_k_.push_back(param_.loop.r_max); }
  loop_profile_.resize(loop_k_.size(), 0.0);
  ParallelFor(
      loop_k_.size(), [&](const std::size_t i) { loop_profile_[i] = ProfileTransform(loop_k_[i]); },
      "Nuclear UPC smooth radial transform");
  math::ValidateInterpolationGrid(loop_k_);
  if (!gra::AllFinite(loop_profile_)) {
    throw std::runtime_error("MUPC::PrepareScalarTransform: non-finite radial transform");
  }
}

// Interpolate the smooth radial profile transform
double MUPC::ScalarTransform(const double kt) const {
  const bool separate = channel_ && !loop_harmonic_.empty();
  const auto& momentum = separate ? channel_->momentum : loop_k_;
  const auto& profile = separate ? channel_->bare : loop_profile_;
  if (!std::isfinite(kt) || momentum.size() < 2 || profile.size() != momentum.size()) {
    throw std::logic_error("MUPC::ScalarTransform: transform grid is invalid");
  }
  if (kt < momentum.front() || kt > momentum.back()) { return ProfileTransform(kt); }
  const auto        upper    = std::upper_bound(momentum.cbegin(), momentum.cend(), kt);
  const std::size_t hi       = std::clamp<std::size_t>(std::distance(momentum.cbegin(), upper), 1, momentum.size() - 1);
  const std::size_t lo       = hi - 1;
  const double      fraction = (separate || param_.loop.radial_map == math::RadialMap::Log)
                                   ? std::clamp(std::log(kt / momentum[lo]) / std::log(momentum[hi] / momentum[lo]), 0.0, 1.0)
                                   : std::clamp((kt - momentum[lo]) / (momentum[hi] - momentum[lo]), 0.0, 1.0);
  return profile[lo] * (1.0 - fraction) + profile[hi] * fraction;
}

// Prepare the impact-parameter and momentum-space convolution
// P_c(b) = 1 - S_c(b)
// w_c(k) = -F_c(k) d^2k / (4 pi^2)
void MUPC::Prepare() {
  PrepareSamples();
  if (!HadronicConvolution()) { return; }
  const bool                  config_profile = glauber_ != nullptr && param_.survival == SurvivalType::MCGGCF;
  std::string                 cache_key;
  std::string                 cache_filename;
  std::unique_ptr<MCacheLock> cache_lock;
  if ((!config_profile || config_grid_ == nullptr) && !param_.table_fingerprint.empty()) {
    cache_key      = CacheKey();
    cache_filename = CacheFilename(cache_key, type_, nucleus_);
    cache_lock     = std::make_unique<MCacheLock>(cache_filename + ".lock");
    const bool loaded =
        config_profile ? ReadConfigCache(cache_filename, cache_key) : ReadCache(cache_filename, cache_key);
    if (loaded) {
      std::cout << "Loaded nuclear UPC cache: " << cache_filename << std::endl;
      if (!config_profile) { return; }
    }
    if (config_profile) {
      if (config_grid_ == nullptr) {
        PrepareConfigGrid();
        if (cache_lock->Acquired()) {
          try {
            WriteConfigCache(cache_filename, cache_key);
            std::cout << "Saved nuclear UPC cache: " << cache_filename << std::endl;
          } catch (const std::exception &error) {
            std::cerr << "WARNING: MUPC cache not saved: " << error.what() << std::endl;
          }
        }
      }
      // The cache contains only immutable quadrature tables. Event-dependent
      // configuration amplitudes are constructed below and are never cached
      cache_lock.reset();
    }
  }
  if (config_profile && config_grid_ == nullptr) { PrepareConfigGrid(); }
  if (config_profile) {
    b_node_   = config_grid_->b_node;
    b_weight_ = config_grid_->b_weight;
    if (!HasSamples()) { return; }
  }
  if (!config_profile) {
    auto b_rule = math::GaussLegendreRule(param_.convolution.smooth_b_nodes, 0.0, param_.convolution.b_max);
    b_node_     = std::move(b_rule.first);
    b_weight_   = std::move(b_rule.second);
  }
  profile_.resize(b_node_.size(), 0.0);
  if (glauber_ != nullptr && param_.survival == SurvivalType::MCGGCF) {
    sample_profile_.assign(sample_count_,
                           std::vector<std::complex<double>>(b_node_.size() * param_.convolution.b_phi_nodes, 0.0));
    if (sample_count_ > std::numeric_limits<std::size_t>::max() / b_node_.size()) {
      throw std::invalid_argument("MUPC: screening ensemble is too large");
    }
    std::vector<std::array<double, 2>> impact;
    impact.reserve(b_node_.size() * param_.convolution.b_phi_nodes);
    for (const double b : b_node_) {
      for (unsigned int angle = 0; angle < param_.convolution.b_phi_nodes; ++angle) {
        impact.push_back({b * config_grid_->b_unit[angle].real(), b * config_grid_->b_unit[angle].imag()});
      }
    }
    // Build configuration profiles serially inside the existing VEGAS worker
    for (const auto &sample : indices(sample_profile_)) {
      const std::size_t upper     = SampleIndex(sample, 1);
      const std::size_t lower     = SampleIndex(sample, 2);
      const auto        amplitude = glauber_->ConfigPairAmp(impact, bank_[0].get(), bank_[1].get(), upper, lower);
      for (const auto &i : indices(b_node_)) {
        for (unsigned int angle = 0; angle < param_.convolution.b_phi_nodes; ++angle) {
          const std::size_t node  = i * static_cast<std::size_t>(param_.convolution.b_phi_nodes) + angle;
          const double      value = 1.0 - amplitude[node];
          for (unsigned int mode = 0; mode < param_.convolution.b_phi_nodes; ++mode) {
            sample_profile_[sample][i * param_.convolution.b_phi_nodes + mode] +=
                value * config_grid_->b_phase[angle * param_.convolution.b_phi_nodes + mode];
          }
        }
      }
    }
    for (const auto &i : indices(profile_)) {
      for (const auto &sample : indices(sample_profile_)) {
        profile_[i] += sample_profile_[sample][i * param_.convolution.b_phi_nodes].real();
      }
      profile_[i] /= static_cast<double>(sample_count_);
    }
  } else {
    ParallelFor(
        b_node_.size(), [&](const std::size_t i) { profile_[i] = 1.0 - SurvivalAmp(b_node_[i]); },
        "Nuclear UPC impact profile");
  }

  loop_profile_.clear();
  loop_harmonic_.clear();
  if (config_profile) {
    loop_k_ = config_grid_->k_node;
    loop_harmonic_ = ConfigTransform(*config_grid_, sample_profile_);
  } else {
    PrepareScalarTransform();
  }
  if (cache_lock && cache_lock->Acquired()) {
    try {
      WriteCache(cache_filename, cache_key);
      std::cout << "Saved nuclear UPC cache: " << cache_filename << std::endl;
    } catch (const std::exception &error) {
      std::cerr << "WARNING: MUPC cache not saved: " << error.what() << std::endl;
    }
  }
}

// Compute one checked beam profile type
BeamType MUPC::Type(const int leg) const { return type_[LegIndex(leg, "MUPC::Type")]; }

// Compute one checked beam nucleus
const MNucleus *MUPC::Nucleus(const int leg) const { return nucleus_[LegIndex(leg, "MUPC::Nucleus")].get(); }

// Compute one optional immutable configuration bank
const MConfigBank *MUPC::Bank(const int leg) const { return bank_[LegIndex(leg, "MUPC::Bank")].get(); }

// Compute one optional event-local hotspot bank
const MHotSpot *MUPC::HotSpot(const int leg) const {
  const auto &hotspot = hotspot_[LegIndex(leg, "MUPC::HotSpot")];
  return hotspot.has_value() ? &*hotspot : nullptr;
}

// Compute whether both beams are hadrons
bool MUPC::HasHadronicPair() const { return type_[0] != BeamType::Lepton && type_[1] != BeamType::Lepton; }

// Compute one optional preconstructed nuclear photon source
const MPhoton *MUPC::Photon(const int leg) const { return photon_[LegIndex(leg, "MUPC::Photon")].get(); }

// Compute one optional preconstructed photonuclear target model
const MPhoto *MUPC::Photo(const int leg) const { return photo_[LegIndex(leg, "MUPC::Photo")].get(); }

// Compute one electromagnetic breakup model
const MBreakup *MUPC::Breakup(const int leg) const { return breakup_[LegIndex(leg, "MUPC::Breakup")].get(); }

// Compute event-local soft and shifted-current multichannel screening nodes
//
// Every channel covers the complete convolution annulus in its own polar
// coordinates.  The soft channel is centered at k = 0, while the two EPA
// channels are centered at k = -q1T and k = +q2T.  Logarithmic radial maps
// resolve the nuclear kernel and the local coherent-current scale |qz|.
// Inverse-distance partitions sum to one pointwise, so the channel sum remains
// the original two-dimensional integral without double counting
// w_a(k) = -F(k) p_a(k) d^2k / (4 pi^2), sum_a p_a(k) = 1
std::vector<LoopNode> MUPC::Nodes(const std::array<M3Vec, 2> &transfer) const {
  const bool scalar  = !loop_profile_.empty();
  const bool sampled = !loop_harmonic_.empty();
  if (!scalar && !sampled) { return {}; }
  if (scalar == sampled || loop_k_.size() < 2 ||
      (sampled && (sample_count_ == 0 ||
                   loop_harmonic_.size() != loop_k_.size() * sample_count_ * param_.convolution.b_phi_nodes))) {
    throw std::logic_error("MUPC::Nodes: screening harmonics are incomplete");
  }

  struct Channel {
    double x     = 0.0;
    double y     = 0.0;
    double scale = 0.0;
  };
  const std::array<Channel, 2> cancellation = {Channel{-transfer[0][0], -transfer[0][1], std::abs(transfer[0][2])},
                                               Channel{transfer[1][0], transfer[1][1], std::abs(transfer[1][2])}};
  double                       axis         = 0.0;
  for (const auto &item : cancellation) {
    if (std::hypot(item.x, item.y) > 0.0) {
      axis = std::atan2(item.y, item.x);
      break;
    }
  }

  // Each importance width follows a physical transverse scale
  std::vector<Channel> channel = {{0.0, 0.0, PDG::GeV2fm / param_.convolution.b_max}};
  if (channel_) { channel.push_back({0.0, 0.0, PDG::GeV2fm / channel_->radius}); }
  for (const auto &item : cancellation) {
    if (!(item.scale > 0.0) || !(std::hypot(item.x, item.y) < param_.loop.r_max)) { continue; }
    channel.push_back(item);
  }
  const unsigned int    radial_count = math::PolarNodeCount(param_.loop);
  const unsigned int    phi_count    = param_.loop.azimuth_nodes;
  std::vector<unsigned int> counts(channel.size(), radial_count);
  if (channel_) {
    // Preserve the configured node density in log(r^2 + width^2) across the longer EMD support
    const double range = std::log1p(math::pow2(param_.loop.r_max / channel.front().scale));
    for (const auto& i : indices(counts)) {
      const double extended = std::log1p(math::pow2(param_.loop.r_max / channel[i].scale));
      counts[i] = std::max(radial_count, static_cast<unsigned int>(std::ceil(radial_count * extended / range)));
    }
  }
  std::vector<LoopNode> output;
  output.reserve(static_cast<std::size_t>(channel.size()) * radial_count * phi_count);

  // Compute one channel's pointwise partition of unity
  const auto partition = [&](const std::size_t selected, const double kx, const double ky) {
    double sum    = 0.0;
    double weight = 0.0;
    for (const auto &i : indices(channel)) {
      const double dx    = kx - channel[i].x;
      const double dy    = ky - channel[i].y;
      const double score = 1.0 / (dx * dx + dy * dy + math::pow2(channel[i].scale));
      sum += score;
      if (i == selected) { weight = score; }
    }
    return weight / sum;
  };

  // Append one Cartesian node using interpolated profile harmonics
  const auto append = [&](const std::size_t selected, const double kx, const double ky, const double measure) {
    LoopNode node;
    node.kt            = std::hypot(kx, ky);
    node.phi           = std::atan2(ky, kx);
    node.kx            = kx;
    node.ky            = ky;
    const double scale = -measure * partition(selected, kx, ky) / (4.0 * math::PIPI);
    if (scalar) {
      node.weight = scale * ScalarTransform(node.kt);
      output.push_back(std::move(node));
      return;
    }
    node.sample_weight.assign(sample_count_, 0.0);
    const auto        upper = std::upper_bound(loop_k_.cbegin(), loop_k_.cend(), node.kt);
    const std::size_t hi    = std::clamp<std::size_t>(std::distance(loop_k_.cbegin(), upper), 1, loop_k_.size() - 1);
    const std::size_t lo    = hi - 1;
    const double      fraction =
        param_.loop.radial_map == math::RadialMap::Log
                 ? std::clamp(std::log(node.kt / loop_k_[lo]) / std::log(loop_k_[hi] / loop_k_[lo]), 0.0, 1.0)
                 : std::clamp((node.kt - loop_k_[lo]) / (loop_k_[hi] - loop_k_[lo]), 0.0, 1.0);
    std::vector<std::complex<double>> phase(param_.convolution.b_phi_nodes, 0.0);
    for (unsigned int mode = 0; mode < param_.convolution.b_phi_nodes; ++mode) {
      phase[mode] =
          std::polar(1.0, math::SignedHarmonic(mode, param_.convolution.b_phi_nodes) * (node.phi + 0.5 * math::PI));
    }
    const double excitation = channel_ ? ScalarTransform(node.kt) : 0.0;
    for (const auto &sample : indices(node.sample_weight)) {
      std::complex<double> transform = 0.0;
      for (unsigned int mode = 0; mode < param_.convolution.b_phi_nodes; ++mode) {
        const std::size_t index0   = (lo * sample_count_ + sample) * param_.convolution.b_phi_nodes + mode;
        const std::size_t index1   = (hi * sample_count_ + sample) * param_.convolution.b_phi_nodes + mode;
        const auto        harmonic = loop_harmonic_[index0] * (1.0 - fraction) + loop_harmonic_[index1] * fraction;
        transform += phase[mode] * harmonic;
      }
      node.sample_weight[sample] = scale * (transform + excitation);
      node.weight += node.sample_weight[sample];
    }
    node.weight /= static_cast<double>(sample_count_);
    output.push_back(std::move(node));
  };

  const double phi_weight = 2.0 * math::PI / static_cast<double>(phi_count);
  for (const auto &selected : indices(channel)) {
    for (unsigned int iphi = 0; iphi < phi_count; ++iphi) {
      const double angle = axis + phi_weight * (static_cast<double>(iphi) + 0.5);
      const auto   interval =
          DiskRayIntervals(channel[selected].x, channel[selected].y, angle, param_.loop.r_min, param_.loop.r_max);
      const auto interval_count = RadialCounts(interval, counts[selected]);
      if (interval_count.size() != interval.size()) {
        throw std::logic_error("MUPC::Nodes: insufficient local radial channel nodes");
      }
      for (const auto &part : indices(interval)) {
        const auto radial =
            LogRadialRule(interval_count[part], interval[part].first, interval[part].second, channel[selected].scale);
        for (const auto &[radius, radial_weight] : radial) {
          const double kx = channel[selected].x + radius * std::cos(angle);
          const double ky = channel[selected].y + radius * std::sin(angle);
          append(selected, kx, ky, radial_weight * phi_weight);
        }
      }
    }
  }
  return output;
}

// Compute whether a configuration-current ensemble is available
bool MUPC::HasSamples() const { return sample_count_ > 0; }

// Compute one leg's configuration index in the paired ensemble
std::size_t MUPC::SampleIndex(const std::size_t sample, const int leg) const {
  const std::size_t index = LegIndex(leg, "MUPC::SampleIndex");
  if (sample >= sample_count_) { throw std::out_of_range("MUPC::SampleIndex: sample is out of range"); }
  return index == 0 ? sample / sample_shape_[1] : sample % sample_shape_[1];
}

// Expand one leg-local current bank over the paired ensemble
std::vector<std::complex<double>> MUPC::ExpandRatios(const std::vector<std::complex<double>> &ratio,
                                                     const int                                leg) const {
  const std::size_t index = LegIndex(leg, "MUPC::ExpandRatios");
  if (ratio.size() != sample_shape_[index] || !gra::AllFinite(ratio)) {
    throw std::invalid_argument("MUPC::ExpandRatios: current bank changed");
  }
  std::vector<std::complex<double>> expanded(sample_count_, 0.0);
  for (const auto &sample : indices(expanded)) { expanded[sample] = ratio[SampleIndex(sample, leg)]; }
  return expanded;
}

// Compute one charge or matter current bank
std::vector<std::complex<double>> MUPC::Currents(const int leg, const M3Vec &q, const bool charge) const {
  const std::size_t index = LegIndex(leg, "MUPC::Currents");
  if (!gra::AllFinite(q) || !HasSamples()) { throw std::logic_error("MUPC::Currents: invalid ensemble or momentum"); }
  if (type_[index] != BeamType::Nucleus) { return {}; }
  if (bank_[index] == nullptr || bank_[index]->Size() == 0) {
    throw std::logic_error("MUPC::Currents: missing configuration bank");
  }
  std::vector<std::complex<double>> current(bank_[index]->Size(), 0.0);
  for (const auto &sample : indices(current)) {
    current[sample] = charge ? bank_[index]->At(sample).ChargeCurrent(q[0], q[1], q[2])
                             : bank_[index]->At(sample).MatterCurrent(q[0], q[1], q[2]);
  }
  return current;
}

// Normalize and expand one leg-local current bank
std::vector<std::complex<double>> MUPC::Ratios(const int leg, const CoherenceType type,
                                               const std::vector<std::complex<double>> &current) const {
  const std::size_t index = LegIndex(leg, "MUPC::Ratios");
  if (!HasSamples()) { throw std::logic_error("MUPC::Ratios: ensemble is unavailable"); }
  if (type_[index] != BeamType::Nucleus) { return std::vector<std::complex<double>>(sample_count_, 1.0); }
  if (bank_[index] == nullptr || current.size() != bank_[index]->Size() || !gra::AllFinite(current)) {
    throw std::invalid_argument("MUPC::Ratios: current bank changed");
  }
  return ExpandRatios(CurrentRatios(current, type), leg);
}

// Include mean and fluctuation sources before projecting a sampled nuclear final state
CoherenceType MUPC::SourceSector(const int leg, const CoherenceType selected) const {
  return Bank(leg) != nullptr ? CoherenceType::Inclusive : selected;
}

// Compute normalized charge currents for one explicit emission sector
// R_coh,c = 1, R_inc,c = (J_c - <J>)/sqrt[Var(J)], R_incl,c = J_c/sqrt[<|J|^2>]
// [REFERENCE: Good and Walker, Phys. Rev. 120 (1960) 1857]
std::vector<std::complex<double>> MUPC::EmissionRatios(const int leg, const CoherenceType emission,
                                                       const M3Vec &q) const {
  return Ratios(leg, emission, Currents(leg, q, true));
}

// Compute coherent and incoherent charge-current ratios in one bank traversal
// (R_coh,c, R_inc,c) = (1, (J_c - <J>)/sqrt[Var(J)])
MUPC::SectorRatios MUPC::EmissionComponents(const int leg, const M3Vec &q) const {
  SectorRatios ratio;
  const auto   current = Currents(leg, q, true);
  ratio[0]             = Ratios(leg, CoherenceType::Coherent, current);
  ratio[1]             = Ratios(leg, CoherenceType::Incoherent, current);
  return ratio;
}

// Compute normalized matter currents for one explicit target sector
// R_coh,c = 1, R_inc,c = (J_c - <J>)/sqrt[Var(J)], R_incl,c = J_c/sqrt[<|J|^2>]
// [REFERENCE: Good and Walker, Phys. Rev. 120 (1960) 1857]
std::vector<std::complex<double>> MUPC::TargetRatios(const int leg, const CoherenceType target, const M3Vec &q) const {
  return Ratios(leg, target, Currents(leg, q, false));
}

// Compute target ratios from one event-local configuration-current bank
// R_coh,c = 1, R_inc,c = (J_c - <J>)/sqrt[Var(J)], R_incl,c = J_c/sqrt[<|J|^2>]
// [REFERENCE: Good and Walker, Phys. Rev. 120 (1960) 1857]
std::vector<std::complex<double>> MUPC::TargetCurrentRatios(const int leg, const CoherenceType target,
                                                            const PhotoCurrent &current) const {
  return Ratios(leg, target, current.sample);
}

// Attach an immutable excitation transform to this event's nuclear configuration
std::shared_ptr<const MUPC> MUPC::WithExcitation(std::shared_ptr<const ExcitationChannel> channel,
                                               const bool survival) const {
  auto event = std::make_shared<MUPC>(*this);
  event->channel_ = std::move(channel);
  event->screening_ = screening_ && survival;
  event->param_.loop.r_min = event->channel_->momentum.front();
  const bool sampled = event->HadronicConvolution() && param_.survival == SurvivalType::MCGGCF;
  if (sampled) {
    // T_h(b) is configuration independent, while every S_c(b,phi) retains its correlations
    // B_h - S_c T_h = (B_h - T_h) + T_h (1 - S_c)
    auto log_radius = event->channel_->impact;
    for (double& b : log_radius) { b = std::log(b); }
    for (const auto& ib : indices(b_node_)) {
      const double excitation = math::LinearInterpolateValidatedGrid(log_radius, event->channel_->amplitude,
                                                                     std::log(b_node_[ib]));
      for (auto& profile : event->sample_profile_) {
        for (unsigned int mode = 0; mode < param_.convolution.b_phi_nodes; ++mode) {
          profile[ib * param_.convolution.b_phi_nodes + mode] *= excitation;
        }
      }
    }
    event->loop_harmonic_ = ConfigTransform(*config_grid_, event->sample_profile_);
    event->loop_profile_.clear();
  } else {
    event->loop_k_ = event->channel_->momentum;
    event->loop_profile_ = survival ? event->channel_->screened : event->channel_->bare;
    event->loop_harmonic_.clear();
    event->param_.loop.radial_map = math::RadialMap::Log;
  }
  return event;
}

// Compute whether hadronic survival requires an impact-space convolution
bool MUPC::HadronicConvolution() const { return screening_ && HasHadronicPair(); }

// Compute the selected no-additional-interaction amplitude
// S(b) = 1, S_opt(b), S_GGCF(b), or <S_c(b)> according to the selected model
double MUPC::SurvivalAmp(const double b) const {
  if (!HadronicConvolution()) { return 1.0; }
  switch (param_.survival) {
    case SurvivalType::Optical:
      return glauber_->OpticalAmp(b);
    case SurvivalType::OpticalGGCF:
      return glauber_->GGCFAmp(b);
    case SurvivalType::MCGGCF:
      return glauber_->ConfigAmp(b, bank_[0].get(), bank_[1].get());
  }
  throw std::logic_error("MUPC::SurvivalAmp: invalid survival type");
}

// Sample one complete event-local nuclear fluctuation ensemble
std::shared_ptr<const MUPC> MUPC::Sample(MRandom &random, const bool prepare_convolution,
                                         const std::size_t count) const {
  std::array<std::shared_ptr<const MConfigBank>, 2> bank;
  std::array<std::optional<MHotSpot>, 2>            hotspot;
  for (const auto &leg : indices(sampler_)) {
    if (sampler_[leg] == nullptr) { continue; }
    const std::size_t samples = count == 0 ? sampler_[leg]->Param().count : count;
    bank[leg]                 = std::make_shared<const MConfigBank>(*sampler_[leg], random, samples);
    if (param_.photo[leg].hotspot.count > 0) { hotspot[leg].emplace(*bank[leg], param_.photo[leg].hotspot, random); }
  }
  return std::shared_ptr<const MUPC>(new MUPC(*this, std::move(bank), std::move(hotspot), prepare_convolution));
}

// Compute the momentum-space hadronic profile in GeV^-2
// F(k) = 2 pi / hbarc^2 int b db J_0(kb/hbarc) [1 - S(b)]
double MUPC::ProfileTransform(const double kt) const {
  if (!std::isfinite(kt) || kt < 0.0) {
    throw std::invalid_argument("MUPC::ProfileTransform: momentum must be finite and nonnegative");
  }
  if (channel_) {
    std::vector<double> field(channel_->impact.size());
    for (const auto& i : indices(field)) {
      field[i] = channel_->born - (loop_harmonic_.empty() ? SurvivalAmp(channel_->impact[i]) : 1.0) *
                                      channel_->amplitude[i];
    }
    field.back() = 0.0;
    return (math::LogBesselKernel(channel_->impact, {kt}, PDG::GeV2fm) * field).front();
  }
  if (!HadronicConvolution()) { return 0.0; }
  const double transform = math::RadialBesselTransform0(b_node_, b_weight_, profile_, kt, PDG::GeV2fm);
  if (!std::isfinite(transform)) { throw std::runtime_error("MUPC::ProfileTransform: non-finite transform"); }
  return transform;
}

}  // namespace gra::nuclear
