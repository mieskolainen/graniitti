// Universal pp screening of helicity and Good Walker amplitudes
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Eikonal/MProtonScreen.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <span>

#include "Graniitti/Math/MMath.h"
#include "Graniitti/Spin/MHelicityScatter.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MException.h"

using gra::aux::indices;

namespace gra::eikonal {
namespace {

// Compute the amplitude scale that preserves the unpolarized norm when one identity row is expanded
double DenseScale(const ScreeningMetadata &metadata) {
  return metadata.spin_basis == ScreeningSpinBasis::ProtonIdentity && metadata.spin_rows == 1 ? 0.5 : 1.0;
}

// Compute the number of non-proton-spin amplitudes in one source bank
std::size_t SpectatorCount(const ScreeningMetadata &metadata, const std::size_t size) {
  if (metadata.spin_rows == 0 || size % metadata.spin_rows != 0) {
    throw AmplitudeFailure("MProtonScreen: amplitude size does not match its proton-spin basis");
  }
  return size / metadata.spin_rows;
}

// Expand compact process amplitudes into the dense proton transition basis
std::vector<std::complex<double>> DenseAmplitude(const std::vector<std::complex<double>> &source,
                                                 const ScreeningMetadata &metadata, const std::size_t spectators) {
  if (metadata.spin_transition_count == 0) {
    throw AmplitudeFailure("MProtonScreen: process has no coherent proton-spin basis");
  }
  std::vector<std::complex<double>> amplitude(16 * spectators, 0.0);
  const double                      scale = DenseScale(metadata);
  for (std::size_t i = 0; i < metadata.spin_transition_count; ++i) {
    const auto &spin = metadata.spin_transition[i];
    for (std::size_t h = 0; h < spectators; ++h) {
      amplitude[spin::PairHelicityTransitionIndex(spin.initial, spin.intermediate) * spectators + h] =
          scale * source[spin.source_row * spectators + h];
    }
  }
  return amplitude;
}

// Add one proton-helicity screening matrix times one compact hard amplitude
void AddHelicity(std::span<const std::complex<double>> source, const ScreeningMetadata &metadata,
                 const ProtonHelicityMatrix &soft, const std::size_t spectators,
                 std::vector<std::complex<double>> &amplitude, const std::size_t offset) {
  const double scale = DenseScale(metadata);
  for (std::size_t i = 0; i < metadata.spin_transition_count; ++i) {
    const auto &spin = metadata.spin_transition[i];
    for (std::size_t final = 0; final < 4; ++final) {
      const auto weight = soft[spin::PairHelicityMatrixIndex(final, spin.intermediate)];
      if (math::IsZero(weight)) { continue; }
      for (std::size_t h = 0; h < spectators; ++h) {
        amplitude[offset + spin::PairHelicityTransitionIndex(spin.initial, final) * spectators + h] +=
            weight * scale * source[spin.source_row * spectators + h];
      }
    }
  }
}

// Multiply one pair-space matrix into a source block and add the result
void AddPair(const MMatrix<std::complex<double>> &matrix, const std::complex<double> *source,
             const std::size_t dimension, std::vector<std::complex<double>> &amplitude, const std::size_t offset,
             const std::complex<double> weight) {
  matrix.MultiplyAdd(std::span<const std::complex<double>>(source, dimension),
                     std::span<std::complex<double>>(amplitude).subspan(offset, dimension), weight);
}

// Compare one shifted Good Walker component with its Born layout
bool SameComponent(const ProtonGoodWalkerComponent &born, const ProtonGoodWalkerComponent &shifted) {
  return born.coherence_group == shifted.coherence_group &&
         born.upper_sector == shifted.upper_sector &&
         born.lower_sector == shifted.lower_sector && born.source.size_row() == shifted.source.size_row() &&
         born.source.size_col() == shifted.source.size_col();
}

}  // namespace

// Expand Born amplitudes into the physical proton and Good Walker basis
MProtonScreen::MProtonScreen(const std::vector<std::vector<std::complex<double>>> &born,
                             const ScreeningMetadata &metadata, const MEikonal::LoopConst &loop,
                             const GWFinalClass final_class)
    : metadata_(metadata), loop_(loop) {
  channel_     = &loop_.good_walker_channels[static_cast<std::size_t>(final_class)];
  good_walker_ = metadata.proton_mode == ProtonScreeningMode::ForwardExcitation &&
                 final_class != GWFinalClass::Elastic && !channel_->empty();
  const std::size_t channel_count = good_walker_ ? channel_->size() : 1;
  if (channel_count == 0) { throw AmplitudeFailure("MProtonScreen: empty Good Walker final-state class"); }
  scalar_good_walker_ =
      good_walker_ && std::all_of(channel_->begin(), channel_->end(), [](const auto &c) { return c.spin_scalar; });
  helicity_ = !loop_.physical_screening_is_scalar || (good_walker_ && !scalar_good_walker_);
  if (helicity_ && metadata.spin_basis != ScreeningSpinBasis::ProtonIdentity &&
      metadata.spin_basis != ScreeningSpinBasis::ProtonHelicity) {
    throw AmplitudeFailure("MProtonScreen: process has no coherent proton-spin basis");
  }

  bank_.reserve(born.size());
  for (const auto &source : born) {
    Bank bank;
    bank.source_size = source.size();
    if (helicity_) {
      bank.spectators  = SpectatorCount(metadata_, source.size());
      const auto dense = DenseAmplitude(source, metadata_, bank.spectators);
      bank.born.assign(channel_count * dense.size(), 0.0);
      for (std::size_t c = 0; c < channel_count; ++c) {
        const double coefficient = good_walker_ ? (*channel_)[c].born_coefficient : 1.0;
        for (const auto &h : indices(dense)) { bank.born[c * dense.size() + h] = coefficient * dense[h]; }
      }
    } else if (scalar_good_walker_) {
      bank.born.assign(channel_count * source.size(), 0.0);
      for (std::size_t c = 0; c < channel_count; ++c) {
        for (const auto &h : indices(source)) {
          bank.born[c * source.size() + h] = (*channel_)[c].born_coefficient * source[h];
        }
      }
    } else {
      bank.born = source;
    }
    bank.loop.assign(bank.born.size(), 0.0);
    bank_.push_back(std::move(bank));
  }
}

// Accumulate A_loop(fi,h) += sum_m S(fi,m;kT) A_hard(i,m;kT) d2kT for every bank
void MProtonScreen::Add(const std::size_t radial, const std::size_t azimuth,
                        const std::vector<std::span<const std::complex<double>>> &amplitude) {
  if (amplitude.size() != bank_.size()) { throw AmplitudeFailure("MProtonScreen: amplitude bank count changed"); }
  const std::size_t n_phi = loop_.physical_screening_weight.size_col();
  const std::size_t node  = radial * n_phi + azimuth;
  for (const auto &b : indices(bank_)) {
    auto &bank = bank_[b];
    if (amplitude[b].size() != bank.source_size) {
      throw AmplitudeFailure("MProtonScreen: external-helicity layout changed");
    }
    if (helicity_) {
      const std::size_t channel_count = good_walker_ ? channel_->size() : 1;
      const std::size_t block         = 16 * bank.spectators;
      for (std::size_t c = 0; c < channel_count; ++c) {
        if (good_walker_ && (*channel_)[c].spin_scalar) {
          ProtonHelicityMatrix soft{};
          for (const std::size_t diagonal : {0U, 5U, 10U, 15U}) {
            soft[diagonal] = (*channel_)[c].scalar_screening_weight[node];
          }
          AddHelicity(amplitude[b], metadata_, soft, bank.spectators, bank.loop, c * block);
          continue;
        }
        const auto &soft =
            good_walker_ ? (*channel_)[c].screening_weight[node] : loop_.physical_screening_helicity_weight[node];
        AddHelicity(amplitude[b], metadata_, soft, bank.spectators, bank.loop, c * block);
      }
    } else if (scalar_good_walker_) {
      for (const auto &c : indices(*channel_)) {
        const std::size_t offset = c * bank.source_size;
        for (const auto &h : indices(amplitude[b])) {
          bank.loop[offset + h] += (*channel_)[c].scalar_screening_weight[node] * amplitude[b][h];
        }
      }
    } else {
      gra::AddScaled(bank.loop, amplitude[b], loop_.physical_screening_weight[radial][azimuth]);
    }
  }
}

// Compute A_screened = A_Born + A_loop for every physical or color-projection bank
std::vector<std::vector<std::complex<double>>> MProtonScreen::Result() const {
  std::vector<std::vector<std::complex<double>>> result;
  result.reserve(bank_.size());
  for (const auto &bank : bank_) {
    result.push_back(bank.born);
    gra::AddScaled(result.back(), bank.loop, std::complex<double>(1.0, 0.0));
  }
  return result;
}

// Initialize dense pair-space blocks from the Born amplitude
MProtonGoodWalkerScreen::MProtonGoodWalkerScreen(const ProtonGoodWalkerAmplitude &born,
                                                 const ScreeningMetadata &metadata, const MEikonal::LoopConst &loop,
                                                 const SoftModelPtr &soft_model, const std::size_t eikonal_channels)
    : born_(born), metadata_(metadata), loop_(loop) {
  if (born.model == nullptr || born.model != soft_model || born.channel_count != eikonal_channels ||
      born.channel_count != born.model->GoodWalker().ChannelCount() || born.components.empty()) {
    throw AmplitudeFailure("MProtonGoodWalkerScreen: hard and eikonal Good Walker spaces disagree");
  }
  const std::size_t dimension = born.model->GoodWalker().PairDimension();
  block_.reserve(born.components.size());
  spectators_.reserve(born.components.size());
  for (const auto &component : born.components) {
    const std::size_t spectators = SpectatorCount(metadata_, component.source.size_row());
    if (component.source.size_col() != dimension) {
      throw AmplitudeFailure("MProtonGoodWalkerScreen: invalid Born pair-space source");
    }
    spectators_.push_back(spectators);
    block_.emplace_back(16 * spectators * dimension, 0.0);
    auto &output = block_.back();
    for (std::size_t i = 0; i < metadata_.spin_transition_count; ++i) {
      const auto &spin = metadata_.spin_transition[i];
      for (std::size_t h = 0; h < spectators; ++h) {
        const std::size_t row = spin.source_row * spectators + h;
        const std::size_t offset =
            (spin::PairHelicityTransitionIndex(spin.initial, spin.intermediate) * spectators + h) * dimension;
        for (std::size_t a = 0; a < dimension; ++a) {
          output[offset + a] = DenseScale(metadata_) * component.source(row, a);
        }
      }
    }
  }
}

// Accumulate the pp eikonal convolution directly in proton spin and Good Walker pair space
void MProtonGoodWalkerScreen::Add(const std::size_t radial, const std::size_t azimuth,
                                  const ProtonGoodWalkerAmplitude &amplitude) {
  if (amplitude.model != born_.model || amplitude.channel_count != born_.channel_count ||
      amplitude.components.size() != born_.components.size()) {
    throw AmplitudeFailure("MProtonGoodWalkerScreen: shifted pair-space layout changed");
  }
  const std::size_t dimension = born_.model->GoodWalker().PairDimension();
  const std::size_t n_phi     = loop_.node_weight.size_col();
  // Jacob-Wick covariance gives exp(i m phi), m=(h1i-h2i-h1f+h2f)/2=-2,...,+2
  const auto    &harmonic   = loop_.pair_screening_harmonic_weight[radial * n_phi + azimuth];
  constexpr auto transition = CanonicalProtonHelicityTransitions();

  for (const auto &c : indices(born_.components)) {
    const auto &born_component = born_.components[c];
    const auto &component      = amplitude.components[c];
    if (!SameComponent(born_component, component)) {
      throw AmplitudeFailure("MProtonGoodWalkerScreen: shifted component layout changed");
    }
    const auto       &soft        = loop_.pair_screening_spin.at(radial);
    const auto       *diagonal    = loop_.pair_screening_pair_diagonal ? &loop_.pair_screening_diagonal.at(radial) : nullptr;
    const bool        spin_scalar = loop_.pair_screening_spin_scalar;
    const std::size_t spectators  = spectators_[c];
    const double      hard_scale  = DenseScale(metadata_);
    for (std::size_t i = 0; i < metadata_.spin_transition_count; ++i) {
      const auto &spin = metadata_.spin_transition[i];
      for (std::size_t final = 0; final < 4; ++final) {
        if (spin_scalar && final != spin.intermediate) { continue; }
        const std::size_t op = spin::PairHelicityMatrixIndex(final, spin.intermediate);
        const auto        weight =
            spin_scalar ? loop_.node_weight[radial][azimuth] : harmonic[transition[op].azimuth_harmonic + 2];
        for (std::size_t h = 0; h < spectators; ++h) {
          const std::size_t row = spin.source_row * spectators + h;
          const std::size_t offset =
              (spin::PairHelicityTransitionIndex(spin.initial, final) * spectators + h) * dimension;
          if (spin_scalar && diagonal != nullptr) {
            gra::DiagonalMultiplyAdd(std::span<const std::complex<double>>(*diagonal),
                                     std::span<const std::complex<double>>(component.source[row], dimension),
                                     std::span<std::complex<double>>(block_[c]).subspan(offset, dimension),
                                     hard_scale * weight);
          } else {
            AddPair(soft[op], component.source[row], dimension, block_[c], offset, hard_scale * weight);
          }
        }
      }
    }
  }
}

namespace {

// Project screened sources with optional density collection removed at compile time
template <bool collect>
std::vector<std::complex<double>> ProjectScreened(const ProtonGoodWalkerAmplitude                      &born,
                                                  const std::vector<std::vector<std::complex<double>>> &block,
                                                  const ScreeningMetadata &metadata, spin::MQMetrics *metrics) {
  const std::size_t dimension = born.model->GoodWalker().PairDimension();
  const auto       &space     = born.model->GoodWalker();
  struct CoherentSource {
    std::size_t                       group;
    ProtonGoodWalkerSector            upper;
    ProtonGoodWalkerSector            lower;
    std::vector<std::complex<double>> amplitude;
    std::optional<std::size_t>        spin_index;
  };
  std::vector<CoherentSource> coherent;
  for (const auto &c : indices(born.components)) {
    const auto &component = born.components[c];
    const auto  found     = std::find_if(coherent.begin(), coherent.end(), [&](const auto &source) {
      return (!collect || source.spin_index == component.spin_index) && source.group == component.coherence_group &&
             source.upper == component.upper_sector && source.lower == component.lower_sector;
    });
    if (found == coherent.end()) {
      coherent.push_back({component.coherence_group, component.upper_sector, component.lower_sector, block[c],
                          collect ? component.spin_index : std::nullopt});
    } else {
      gra::AddScaled(found->amplitude, block[c], std::complex<double>(1.0, 0.0));
    }
  }
  std::vector<std::complex<double>> result;
  for (const auto &source : coherent) {
    const auto index = collect ? (source.spin_index ? source.spin_index : metrics->FinalPair()) : std::nullopt;
    MMatrix<std::complex<double>> amplitude;
    for (std::size_t offset = 0; offset < source.amplitude.size(); offset += dimension) {
      const auto physical =
          space.ProjectPair(std::span<const std::complex<double>>(source.amplitude.data() + offset, dimension),
                            SectorFinalBasis(source.upper), SectorFinalBasis(source.lower));
      if (!collect || !source.spin_index) { result.insert(result.end(), physical.begin(), physical.end()); }
      if constexpr (collect) {
        if (index) {
          if (offset == 0) { amplitude.Resize(source.amplitude.size() / dimension, physical.size()); }
          std::copy(physical.begin(), physical.end(), amplitude.Row(offset / dimension).begin());
        }
      }
    }
    if constexpr (collect) {
      if (index) { metrics->Add(*index, amplitude, metadata.amplitude_normalization, 16); }
    }
  }
  return result;
}

}  // namespace

// Project the completed amplitude and select density collection outside the spin loops
std::vector<std::complex<double>> MProtonGoodWalkerScreen::Result(spin::MQMetrics *metrics) const {
  return metrics != nullptr ? ProjectScreened<true>(born_, block_, metadata_, metrics)
                            : ProjectScreened<false>(born_, block_, metadata_, metrics);
}

}  // namespace gra::eikonal
