// Regge forward sources and Good Walker pair amplitudes
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Regge/MReggeSource.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <compare>
#include <stdexcept>
#include <utility>

#include "Graniitti/Photon/MQED.h"
#include "Graniitti/Nuclear/MBeam.h"
#include "Graniitti/Regge/MReggeForm.h"
#include "Graniitti/Spin/MHelicityScatter.h"
#include "Graniitti/Tech/MException.h"

namespace gra {
namespace {

// Build the inclusive proton target source with its physical forward coupling
std::vector<LegSource> PhotoDissProfile(const ForwardLegState &state, const MRegge &regge, const regge::Param &param, const int exchange_pdg) {
  if (!nuclear::IsProton(state.emitter.pdg)) { throw AmplitudeFailure("MRegge: photoproduction dissociation requires a proton target"); }
  const auto trajectory = regge::TrajectoryIndex(param, exchange_pdg);
  if (trajectory != param.pomeron_trajectory) { throw AmplitudeFailure("MRegge: photoproduction dissociation requires the mapped Pomeron"); }
  const auto &model = *regge.SoftModelHandle();
  const auto &space = model.GoodWalker();
  std::vector<std::complex<double>> proton(space.ProtonVector().begin(), space.ProtonVector().end());
  gra::Scale(proton, model.PhysicalResidue(param.exchanges[trajectory].soft_exchange, 0.0));
  return {{ProtonGoodWalkerSector::PhotoDiss, std::move(proton), std::vector<std::complex<double>>(space.ChannelCount(), 0.0)}};
}

// Fill one pair-source row directly in canonical upper-major order
void FillPairRow(MMatrix<std::complex<double>> &matrix, std::size_t row, const std::vector<std::complex<double>> &upper, const std::vector<std::complex<double>> &lower, std::complex<double> scale) {
  if (row >= matrix.size_row() || matrix.size_col() != upper.size() * lower.size()) { throw AmplitudeFailure("MRegge: invalid pair-source row dimensions"); }
  auto output = matrix.Row(row);
  for (const auto &i : aux::indices(upper)) {
    for (const auto &j : aux::indices(lower)) { output[i * lower.size() + j] = scale * upper[i] * lower[j]; }
  }
}

}  // namespace

// Bind one call-local cache to immutable forward amplitude inputs
ReggeSourceCache::ReggeSourceCache(const LORENTZSCALAR &lts, const MRegge &regge, const regge::Param &param, const ForwardLegState &upper, const ForwardLegState &lower)
    : ReggeSourceCache(lts, regge, param, upper, lower, {lts.s1, lts.s2}, {}) {}

// Bind explicit subenergies and optional photon-source amplitudes
ReggeSourceCache::ReggeSourceCache(const LORENTZSCALAR &lts, const MRegge &regge, const regge::Param &param, const ForwardLegState &upper, const ForwardLegState &lower, std::array<double, 2> subenergy,
                                   std::array<std::optional<std::complex<double>>, 2> photon_source)
    : lts(lts), regge(regge), param(param), upper(upper), lower(lower), subenergy(subenergy), photon_source(photon_source) {
  if (upper.leg != ForwardBeamLeg::Upper || lower.leg != ForwardBeamLeg::Lower) { throw AmplitudeFailure("ReggeSourceCache: forward leg ordering is invalid"); }
  std::size_t expected_kernels = 2 * lts.process.CONT_PRODUCTION.size() + 2;
  for (const auto &[name, resonance] : lts.process.RESONANCES) {
    static_cast<void>(name);
    expected_kernels += 2 * resonance.production.size();
  }
  max_profiles = 2 * (param.PDG_TO_INDEX.size() + 1);
  profiles.reserve(max_profiles);
  kernels.reserve(expected_kernels);
}

// Compute the immutable forward state selected at construction
const ForwardLegState &ReggeSourceCache::State(const ForwardBeamLeg leg) const {
  if (leg == ForwardBeamLeg::Upper) { return upper; }
  if (leg == ForwardBeamLeg::Lower) { return lower; }
  throw AmplitudeFailure("ReggeSourceCache: invalid forward leg");
}

// Compute one caller supplied forward subenergy
double ReggeSourceCache::Subenergy(const ForwardBeamLeg leg) const {
  if (leg != ForwardBeamLeg::Upper && leg != ForwardBeamLeg::Lower) { throw AmplitudeFailure("ReggeSourceCache: invalid forward leg"); }
  return subenergy[leg == ForwardBeamLeg::Upper ? 0 : 1];
}

// Canonicalize arguments which are absent from the selected kernel equation
ReggeSourceCache::KernelKey ReggeSourceCache::Key(const ForwardBeamLeg leg, const int exchange_pdg, const double s_forward, const bool photoproduction, const int central_pdg) const {
  if (exchange_pdg == PDG::PDG_gamma) { return {leg, exchange_pdg, 0, false, 0}; }
  if (!std::isfinite(s_forward)) { throw AmplitudeFailure("ReggeSourceCache: non-finite subenergy"); }
  return {leg, exchange_pdg, std::bit_cast<std::uint64_t>(s_forward), photoproduction, photoproduction ? central_pdg : 0};
}

// Compute one Good Walker profile computed at most once in this call
const std::vector<LegSource> &ReggeSourceCache::Profile(const ForwardBeamLeg leg, const int exchange_pdg) {
  const ProfileKey key{leg, exchange_pdg};
  const auto       found = std::find_if(profiles.cbegin(), profiles.cend(), [&key](const auto &entry) { return entry.first == key; });
  if (found != profiles.cend()) { return found->second; }
  if (profiles.size() >= max_profiles) { throw std::logic_error("ReggeSourceCache: profile capacity invariant failed"); }
  const std::size_t index = leg == ForwardBeamLeg::Upper ? 0 : 1;
  profiles.emplace_back(key, LegProfile(lts, State(leg), regge, param, exchange_pdg, photon_source[index]));
  return profiles.back().second;
}

// Compute one scalar Regge kernel computed at most once in this call
std::complex<double> ReggeSourceCache::Kernel(const ForwardBeamLeg leg, const int exchange_pdg, const double s_forward, const bool photoproduction, const int central_pdg) {
  const KernelKey key   = Key(leg, exchange_pdg, s_forward, photoproduction, central_pdg);
  const auto      found = std::find_if(kernels.cbegin(), kernels.cend(), [&key](const auto &entry) { return entry.first == key; });
  if (found != kernels.cend()) { return found->second; }
  kernels.emplace_back(key, LegKernel(lts, State(leg), regge, param, exchange_pdg, s_forward, photoproduction, central_pdg));
  return kernels.back().second;
}

// Compute one complete forward source computed at most once in this call
const std::vector<LegSource> &ReggeSourceCache::Sources(const ForwardBeamLeg leg, const int exchange_pdg, const double s_forward, const bool photoproduction, const int central_pdg) {
  const KernelKey key   = Key(leg, exchange_pdg, s_forward, photoproduction, central_pdg);
  const auto      found = std::find_if(sources.cbegin(), sources.cend(), [&key](const auto &entry) { return entry.first == key; });
  if (found != sources.cend()) { return found->second; }
  const bool diss = photoproduction && exchange_pdg != PDG::PDG_gamma && State(leg).IsExcited();
  auto source = diss ? PhotoDissProfile(State(leg), regge, param, exchange_pdg) : Profile(leg, exchange_pdg);
  const std::complex<double> kernel = Kernel(leg, exchange_pdg, s_forward, photoproduction, central_pdg);
  for (auto &component : source) {
    gra::Scale(component.nonflip, kernel);
    gra::Scale(component.flip, kernel);
  }
  sources.emplace_back(key, std::move(source));
  return sources.back().second;
}

// Compute the real QFT sign left by joining reduced amplitudes
// Each internal line gives (iM_a)(i/D)(iM_b)/i = -M_a M_b/D
// This is i-factor bookkeeping, not a statistics sign or fitted phase
double SewingSign(const std::size_t internal_lines) { return (internal_lines % 2 == 0) ? 1.0 : -1.0; }

// Apply one common factor to every proton Good Walker component
void ScaleGoodWalker(LORENTZSCALAR &lts, const double scale, const std::string &context) {
  if (!lts.proton_good_walker) { throw AmplitudeFailure(context + ": missing Good Walker source"); }
  if (std::is_eq(scale <=> 1.0)) { return; }
  for (auto &component : lts.proton_good_walker->components) { component.source *= scale; }
}

// Compute the Born norm with the configured initial spin average
double BornNorm(const LORENTZSCALAR &lts) {
  const double average = lts.hamp.metadata.amplitude_normalization;
  if (!std::isfinite(average) || average <= 0.0) { throw AmplitudeFailure("MRegge: invalid initial-state spin average"); }
  return average * SquaredNorm(lts.hamp);
}

// Retain only the proton Good Walker source during shifted screening
bool GoodWalkerOnly(LORENTZSCALAR &lts, const std::string &context) {
  if (!lts.screening.active) { return false; }
  if (!lts.proton_good_walker) { throw AmplitudeFailure(context + ": missing shifted Good Walker source"); }
  lts.hamp.clear();
  return true;
}

// Compute the configured hard spin layout or the standalone Regge default
ScreeningMetadata ProtonSpinLayout(const gra::LORENTZSCALAR &lts) {
  if (lts.hamp.metadata.spin_basis == ScreeningSpinBasis::ProtonIdentity || lts.hamp.metadata.spin_basis == ScreeningSpinBasis::ProtonHelicity) { return lts.hamp.metadata; }
  ScreeningMetadata layout;
  layout.spin_basis     = ScreeningSpinBasis::ProtonHelicity;
  layout.forward_noflip = lts.process.FORWARD_NOFLIP;
  layout.spin_rows      = layout.forward_noflip ? 4 : 16;
  layout.PrepareSpinTransitions();
  return layout;
}

// Initialize or validate the event-local proton Good Walker amplitude
void PrepareGoodWalker(gra::LORENTZSCALAR &lts, const SoftModelPtr &model) {
  if (model == nullptr) { throw AmplitudeFailure("MRegge: missing immutable SOFT model"); }
  const std::size_t n       = model->GoodWalker().ChannelCount();
  auto             &payload = lts.proton_good_walker;
  if (!payload) {
    payload.emplace();
    auto &state         = *payload;
    state.model         = model;
    state.channel_count = n;
    return;
  }
  const auto &state = *payload;
  if (state.model != model || state.channel_count != n) { throw AmplitudeFailure("MRegge: event Good Walker layout changed"); }
}

// Build the Good-Walker profile carried by one forward beam leg
std::vector<LegSource> LegProfile(const gra::LORENTZSCALAR &lts, const ForwardLegState &state, const gra::MRegge &regge, const regge::Param &param, const int exchange_pdg, const std::optional<std::complex<double>> photon_source) {
  if (!std::isfinite(state.t) || state.t > 0.0) { throw AmplitudeFailure("MRegge::LegProfile: generated transfer is invalid"); }
  const auto                       &model       = *regge.SoftModelHandle();
  const auto                       &good_walker = model.GoodWalker();
  const std::size_t                 n           = good_walker.ChannelCount();
  std::vector<std::complex<double>> proton(good_walker.ProtonVector().begin(), good_walker.ProtonVector().end());
  if (exchange_pdg == PDG::PDG_gamma) {
    const std::complex<double> flux = photon_source.has_value() ? *photon_source : gra::qed::PhotonSourceAmplitude(lts, state, lts.process.PHOTON_VERTEX);
    gra::Scale(proton, flux);
    const auto common = proton;
    return {{state.IsExcited() ? ProtonGoodWalkerSector::InelasticEPA : ProtonGoodWalkerSector::Elastic, common, std::move(proton)}};
  }

  const std::size_t    trajectory = regge::TrajectoryIndex(param, exchange_pdg);
  const SoftExchangeId exchange   = param.exchanges.at(trajectory).soft_exchange;

  if (!state.IsExcited()) {
    const double helicity_mass = model.Eikonal().HelicityMassScale();
    if (!std::isfinite(helicity_mass) || helicity_mass <= 0.0) {
      throw AmplitudeFailure(
          "MRegge::LegProfile: invalid fitted helicity mass "
          "scale");
    }
    auto         flip         = model.HelicityFlipResidueMatrix(exchange, state.t) * proton;
    const double flip_barrier = std::sqrt(-state.t) / (2.0 * helicity_mass);
    if (!std::isfinite(flip_barrier)) {
      throw AmplitudeFailure(
          "MRegge::LegProfile: non-finite generated flip "
          "barrier");
    }
    gra::Scale(flip, flip_barrier);
    return {{ProtonGoodWalkerSector::Elastic, model.ResidueMatrix(exchange, state.t) * proton, std::move(flip)}};
  }
  if (model.Exchange(exchange).Role() != SoftExchangeRole::Pomeron) { return {}; }

  auto         complete           = regge.TriplePomeronCouplingRoot() * proton;
  const double excitation_profile = model.ForwardExcitationFactor(exchange, state.t, state.mass2);
  gra::Scale(complete, excitation_profile);
  auto resolved  = good_walker.ResolvedProjector() * complete;
  auto inclusive = good_walker.InclusiveProjector() * complete;
  if (n == 1) { return {{ProtonGoodWalkerSector::TripleInclusive, std::move(inclusive), std::vector<std::complex<double>>(n, 0.0)}}; }
  return {{ProtonGoodWalkerSector::TripleResolved, std::move(resolved), std::vector<std::complex<double>>(n, 0.0)}, {ProtonGoodWalkerSector::TripleInclusive, std::move(inclusive), std::vector<std::complex<double>>(n, 0.0)}};
}

// Compute the scalar propagator and beam-sign factor of one forward leg
std::complex<double> LegKernel(const gra::LORENTZSCALAR &lts, const ForwardLegState &state, const gra::MRegge &regge, const regge::Param &param, const int exchange_pdg, const double s_forward, const bool photoproduction,
                               const int central_pdg) {
  if (exchange_pdg == PDG::PDG_gamma) { return 1.0; }
  const std::complex<double> propagator = photoproduction
      ? (state.IsExcited() ? regge.PhotoDissKernel(s_forward, state.t, state.mass2, central_pdg, lts.process.PHOTO_DISSOCIATION) : regge.PhotoKernel(s_forward, state.t, central_pdg))
      : regge.ExchangeKernel(s_forward, state.t, exchange_pdg);
  return regge::AntiparticleSign(param, exchange_pdg, state.emitter.pdg) * propagator;
}

// Build all orthogonal Regge sources carried by one forward beam leg
std::vector<LegSource> LegSources(const gra::LORENTZSCALAR &lts, const ForwardLegState &state, const gra::MRegge &regge, const regge::Param &param, const int exchange_pdg, const double s_forward, const bool photoproduction,
                                  const int central_pdg) {
  const bool diss = photoproduction && exchange_pdg != PDG::PDG_gamma && state.IsExcited();
  auto source = diss ? PhotoDissProfile(state, regge, param, exchange_pdg) : LegProfile(lts, state, regge, param, exchange_pdg);
  const std::complex<double> kernel = LegKernel(lts, state, regge, param, exchange_pdg, s_forward, photoproduction, central_pdg);
  for (auto &component : source) {
    gra::Scale(component.nonflip, kernel);
    gra::Scale(component.flip, kernel);
  }
  return source;
}

// Form every upper-lower Cartesian Regge-source combination
std::vector<PairSource> PairSources(const std::vector<LegSource> &upper, const std::vector<LegSource> &lower, const GoodWalkerSpace &good_walker, const std::complex<double> scale) {
  ScreeningMetadata scalar;
  scalar.spin_basis     = ScreeningSpinBasis::ProtonIdentity;
  scalar.spin_rows      = 1;
  scalar.forward_noflip = true;
  return PairSources(upper, lower, good_walker, scalar, scale);
}

// Form every upper-lower source combination in the hard proton spin basis
std::vector<PairSource> PairSources(const std::vector<LegSource> &upper, const std::vector<LegSource> &lower, const GoodWalkerSpace &good_walker, const ScreeningMetadata &spin_layout, const std::complex<double> scale) {
  std::vector<PairSource> out;
  out.reserve(upper.size() * lower.size());
  const std::size_t n = good_walker.ChannelCount();
  const std::size_t d = good_walker.PairDimension();
  for (const auto &up : upper) {
    for (const auto &dn : lower) {
      if (up.nonflip.size() != n || up.flip.size() != n || dn.nonflip.size() != n || dn.flip.size() != n) { throw AmplitudeFailure("MRegge: unequal Good Walker leg dimensions"); }
      PairSource pair;
      pair.upper_sector  = up.sector;
      pair.lower_sector  = dn.sector;
      pair.amplitude     = MMatrix<std::complex<double>>(spin_layout.spin_rows, d, std::complex<double>(0.0, 0.0));
      const auto set_row = [&](const std::size_t row, const auto &upper_profile, const auto &lower_profile) { FillPairRow(pair.amplitude, row, upper_profile, lower_profile, scale); };

      if (spin_layout.spin_basis == ScreeningSpinBasis::ProtonIdentity) {
        if (spin_layout.spin_rows != 1 && spin_layout.spin_rows != 4) { throw AmplitudeFailure("MRegge: invalid proton identity pair-source layout"); }
        for (std::size_t row = 0; row < spin_layout.spin_rows; ++row) { set_row(row, up.nonflip, dn.nonflip); }
      } else if (spin_layout.spin_basis == ScreeningSpinBasis::ProtonHelicity) {
        if (spin_layout.spin_rows != (spin_layout.forward_noflip ? 4 : 16)) { throw AmplitudeFailure("MRegge: invalid proton helicity pair-source layout"); }
        for (const std::size_t i1 : spin::BinaryHelicityIndices()) {
          for (const std::size_t i2 : spin::BinaryHelicityIndices()) {
            if (spin_layout.forward_noflip) {
              set_row(spin::CanonicalProtonPairSpinLayout::CompactNoFlipHardRow(i1, i2), up.nonflip, dn.nonflip);
              continue;
            }
            for (const std::size_t f1 : spin::BinaryHelicityIndices()) {
              for (const std::size_t f2 : spin::BinaryHelicityIndices()) { set_row(spin::CanonicalProtonPairSpinLayout::HardRow(i1, i2, f1, f2), i1 == f1 ? up.nonflip : up.flip, i2 == f2 ? dn.nonflip : dn.flip); }
            }
          }
        }
      } else {
        throw AmplitudeFailure("MRegge: unsupported spin basis for a Good Walker source");
      }
      out.push_back(std::move(pair));
    }
  }
  return out;
}

namespace {

// Build pair-space beam sources for one ordered central-production row
template <bool collect>
void AddPairSource(gra::LORENTZSCALAR &lts, const std::vector<PairSource> &pairs,
                   const MMatrix<std::complex<double>> &central, const std::size_t group, const SoftModelPtr &model,
                   const std::optional<std::size_t> spin_index) {
  if (pairs.empty()) { return; }
  if (central.isEmpty()) { throw AmplitudeFailure("MRegge: empty central helicity matrix"); }
  PrepareGoodWalker(lts, model);
  auto       &state       = *lts.proton_good_walker;
  const auto &good_walker = state.model->GoodWalker();
  if (state.channel_count != good_walker.ChannelCount()) {
    throw AmplitudeFailure("MRegge: event Good Walker channel count changed");
  }
  const std::size_t d = good_walker.PairDimension();
  for (const auto &pair : pairs) {
    if (pair.amplitude.size_col() != d || pair.amplitude.size_row() == 0 ||
        central.size_row() % pair.amplitude.size_row() != 0) {
      throw AmplitudeFailure("MRegge: incompatible central and Good Walker spin rows");
    }
    const std::size_t spectators  = central.size_row() / pair.amplitude.size_row();
    const std::size_t source_rows = central.size_row() * central.size_col();
    auto match = std::find_if(state.components.begin(), state.components.end(), [&](const auto &component) {
      return (!collect || component.spin_index == spin_index) && component.coherence_group == group &&
             component.upper_sector == pair.upper_sector && component.lower_sector == pair.lower_sector;
    });
    if (match == state.components.end()) {
      ProtonGoodWalkerComponent component;
      if constexpr (collect) { component.spin_index = spin_index; }
      component.coherence_group = group;
      component.upper_sector    = pair.upper_sector;
      component.lower_sector    = pair.lower_sector;
      component.source          = MMatrix<std::complex<double>>(source_rows, d, 0.0);
      state.components.push_back(std::move(component));
      match = std::prev(state.components.end());
    }
    if (match->source.size_row() != source_rows || match->source.size_col() != d) {
      throw AmplitudeFailure("MRegge: incompatible coherent Good Walker source layout");
    }
    for (std::size_t row = 0; row < central.size_row(); ++row) {
      const std::size_t hard_row = row / spectators;
      for (std::size_t col = 0; col < central.size_col(); ++col) {
        auto target = match->source.Row(row * central.size_col() + col);
        gra::AddScaled(target, pair.amplitude.Row(hard_row), central[row][col]);
      }
    }
  }
}

// Project every Born pair source onto its physical final-state sector
template <bool collect>
std::vector<std::complex<double>> ProjectBorn(const ProtonGoodWalkerAmplitude &state, spin::MQMetrics *metrics,
                                              const double normalization, const std::size_t proton_rows) {
  if (state.model == nullptr || state.channel_count == 0 || state.components.empty()) {
    throw AmplitudeFailure("MRegge: inactive proton Good Walker amplitude");
  }
  struct CoherentSource {
    std::size_t                   group = 0;
    ProtonGoodWalkerSector        upper = ProtonGoodWalkerSector::Elastic;
    ProtonGoodWalkerSector        lower = ProtonGoodWalkerSector::Elastic;
    MMatrix<std::complex<double>> source;
    std::optional<std::size_t>    spin_index;
  };
  std::vector<CoherentSource> coherent;
  const auto                 &good_walker = state.model->GoodWalker();
  if (state.channel_count != good_walker.ChannelCount()) {
    throw AmplitudeFailure("MRegge: event Good Walker channel count changed");
  }
  const std::size_t d = good_walker.PairDimension();
  for (const auto &component : state.components) {
    if (component.source.size_col() != d) { throw AmplitudeFailure("MRegge: invalid pair-space Born source"); }
    const auto match = std::find_if(coherent.begin(), coherent.end(), [&](const auto &item) {
      return (!collect || item.spin_index == component.spin_index) && item.group == component.coherence_group &&
             item.upper == component.upper_sector && item.lower == component.lower_sector;
    });
    if (match == coherent.end()) {
      coherent.push_back({component.coherence_group, component.upper_sector, component.lower_sector, component.source,
                          collect ? component.spin_index : std::nullopt});
    } else {
      if (match->source.size_row() != component.source.size_row() ||
          match->source.size_col() != component.source.size_col()) {
        throw AmplitudeFailure("MRegge: incompatible coherent pair-space Born sources");
      }
      match->source += component.source;
    }
  }
  std::vector<std::complex<double>> out;
  for (const auto &item : coherent) {
    const auto index = collect ? (item.spin_index ? item.spin_index : metrics->FinalPair()) : std::nullopt;
    MMatrix<std::complex<double>> amplitude;
    for (std::size_t row = 0; row < item.source.size_row(); ++row) {
      const auto projected =
          good_walker.ProjectPair(item.source.Row(row), SectorFinalBasis(item.upper), SectorFinalBasis(item.lower));
      if (!collect || !item.spin_index) { out.insert(out.end(), projected.begin(), projected.end()); }
      if constexpr (collect) {
        if (index) {
          if (row == 0) { amplitude.Resize(item.source.size_row(), projected.size()); }
          std::copy(projected.begin(), projected.end(), amplitude.Row(row).begin());
        }
      }
    }
    if constexpr (collect) {
      if (index) { metrics->Add(*index, amplitude, normalization, proton_rows); }
    }
  }
  return out;
}

}  // namespace

// Add the physical source without constructing optional diagnostic arguments
void AddGoodWalker(LORENTZSCALAR &lts, const std::vector<PairSource> &pairs,
                   const MMatrix<std::complex<double>> &central, const std::size_t group, const SoftModelPtr &model) {
  if (lts.qmetrics.Active()) {
    AddPairSource<true>(lts, pairs, central, group, model, std::nullopt);
  } else {
    AddPairSource<false>(lts, pairs, central, group, model, std::nullopt);
  }
}

// Add one production source inside the enabled density-collection branch
void AddGoodWalker(LORENTZSCALAR &lts, const std::vector<PairSource> &pairs,
                   const MMatrix<std::complex<double>> &central, const std::size_t group, const SoftModelPtr &model,
                   const std::size_t spin_index) {
  if (lts.qmetrics.Active()) { AddPairSource<true>(lts, pairs, central, group, model, spin_index); }
}

// Select density collection once before entering the amplitude loops
std::vector<std::complex<double>> ProjectGoodWalker(const ProtonGoodWalkerAmplitude &state, spin::MQMetrics *metrics,
                                                    const double normalization, const std::size_t proton_rows) {
  return metrics != nullptr ? ProjectBorn<true>(state, metrics, normalization, proton_rows)
                            : ProjectBorn<false>(state, metrics, normalization, proton_rows);
}

// Project the completed proton Good Walker source into the Born amplitude
void MRegge::ProjectGoodWalkerBorn(gra::LORENTZSCALAR &lts) const {
  if (!lts.proton_good_walker) {
    throw AmplitudeFailure("MRegge::ProjectGoodWalkerBorn: missing Good Walker amplitude");
  }
  const auto &state = *lts.proton_good_walker;
  if (state.model != soft_model) {
    throw AmplitudeFailure("MRegge::ProjectGoodWalkerBorn: SOFT model identity mismatch");
  }
  auto *metrics = lts.qmetrics.Active() ? &lts.qmetrics : nullptr;
  if (metrics != nullptr) { metrics->ClearEvent(); }
  lts.hamp = ProjectGoodWalker(state, metrics, lts.hamp.metadata.amplitude_normalization,
                               metrics != nullptr ? ProtonSpinLayout(lts).spin_rows : 0);
}

// Store one elastic pair source for an arbitrary Good Walker channel count
void MRegge::StoreElasticGoodWalker(gra::LORENTZSCALAR &lts, const SoftExchangeId exchange,
                                    const std::complex<double> kernel) const {
  lts.proton_good_walker.reset();
  const auto                             &good_walker = soft_model->GoodWalker();
  const std::vector<std::complex<double>> proton(good_walker.ProtonVector().begin(), good_walker.ProtonVector().end());
  const auto                              residue = soft_model->ResidueMatrix(exchange, lts.t);
  const LegSource                         elastic{ProtonGoodWalkerSector::Elastic, residue * proton, std::vector<std::complex<double>>(good_walker.ChannelCount(), 0.0)};
  const auto                              pairs = PairSources({elastic}, {elastic}, good_walker, kernel);
  MMatrix<std::complex<double>>           central(1, 1, 1.0);
  AddGoodWalker(lts, pairs, central, 0, soft_model);
  ProjectGoodWalkerBorn(lts);
}

// Store one triple-Good Walker source without an additional inelastic profile
void MRegge::StoreTripleGoodWalker(gra::LORENTZSCALAR &lts, const MReggeInclusive mode, const std::complex<double> kernel) const {
  if (mode == MReggeInclusive::EL) { throw AmplitudeFailure("MRegge::StoreTripleGoodWalker: mode must be SD or DD"); }
  lts.proton_good_walker.reset();
  const auto                       &good_walker = soft_model->GoodWalker();
  const std::size_t                 n           = good_walker.ChannelCount();
  std::vector<std::complex<double>> proton(good_walker.ProtonVector().begin(), good_walker.ProtonVector().end());
  const SoftExchangeId              pomeron           = param.exchanges.at(param.pomeron_trajectory).soft_exchange;
  const auto                        elastic_matrix    = soft_model->ResidueMatrix(pomeron, lts.t);
  auto                              elastic_amplitude = elastic_matrix * proton;
  const auto                        triple_complete   = triple_pomeron_root * proton;
  auto                              triple_resolved   = good_walker.ResolvedProjector() * triple_complete;
  auto                              triple_inclusive  = good_walker.InclusiveProjector() * triple_complete;

  const LegSource        elastic{ProtonGoodWalkerSector::Elastic, std::move(elastic_amplitude), std::vector<std::complex<double>>(n, 0.0)};
  std::vector<LegSource> triple;
  if (n > 1) { triple.push_back({ProtonGoodWalkerSector::TripleResolved, std::move(triple_resolved), std::vector<std::complex<double>>(n, 0.0)}); }
  triple.push_back({ProtonGoodWalkerSector::TripleInclusive, std::move(triple_inclusive), std::vector<std::complex<double>>(n, 0.0)});

  std::vector<PairSource> pairs;
  if (mode == MReggeInclusive::SD) {
    pairs = lts.excite1 ? PairSources(triple, {elastic}, good_walker, 1.0) : PairSources({elastic}, triple, good_walker, 1.0);
  } else {
    pairs = PairSources(triple, triple, good_walker, 1.0);
  }
  MMatrix<std::complex<double>> central(1, 1, 0.0);
  central(0, 0) = kernel;
  AddGoodWalker(lts, pairs, central, 0, soft_model);
}

// Add one direct multi-Regge ladder term to the proton Good Walker source
void MRegge::AddLadderGoodWalker(gra::LORENTZSCALAR &lts, const ForwardLegState &upper_state, const ForwardLegState &lower_state, const int upper_exchange, const double upper_s, const int lower_exchange, const double lower_s,
                                 const std::complex<double> central) const {
  const auto              upper  = LegSources(lts, upper_state, *this, param, upper_exchange, upper_s, false, 0);
  const auto              lower  = LegSources(lts, lower_state, *this, param, lower_exchange, lower_s, false, 0);
  const ScreeningMetadata layout = ProtonSpinLayout(lts);
  const auto              pairs  = PairSources(upper, lower, soft_model->GoodWalker(), layout, 1.0);
  const std::size_t       rows   = layout.spin_rows;
  if (rows != 1 && rows != 4 && rows != 16) { throw AmplitudeFailure("MRegge::AddLadderGoodWalker: invalid proton spin row count"); }
  MMatrix<std::complex<double>> central_matrix(rows, 1, 0.0);
  if (layout.spin_basis == ScreeningSpinBasis::ProtonHelicity && rows == 16) {
    const double upper_phi = (upper_state.outgoing - upper_state.incoming).Phi();
    const double lower_phi = (lower_state.incoming - lower_state.outgoing).Phi();
    for (const std::size_t i1 : spin::BinaryHelicityIndices()) {
      for (const std::size_t i2 : spin::BinaryHelicityIndices()) {
        for (const std::size_t f1 : spin::BinaryHelicityIndices()) {
          for (const std::size_t f2 : spin::BinaryHelicityIndices()) {
            const double h1i = spin::BinaryHelicityLabelX2(i1) / 2.0;
            const double h2i = spin::BinaryHelicityLabelX2(i2) / 2.0;
            const double h1f = spin::BinaryHelicityLabelX2(f1) / 2.0;
            const double h2f = spin::BinaryHelicityLabelX2(f2) / 2.0;
            central_matrix(spin::CanonicalProtonPairSpinLayout::HardRow(i1, i2, f1, f2), 0) =
                central * spin::SpinHalfForwardHelicitySectionFactor(h1i, h1f, upper_phi, false) * spin::SpinHalfForwardHelicitySectionFactor(h2i, h2f, lower_phi, true);
          }
        }
      }
    }
  } else {
    for (std::size_t row = 0; row < rows; ++row) { central_matrix(row, 0) = central; }
  }
  AddGoodWalker(lts, pairs, central_matrix, 0, soft_model);
}

}  // namespace gra
