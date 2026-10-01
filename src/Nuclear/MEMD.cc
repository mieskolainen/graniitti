// Resolved electromagnetic excitation amplitudes and normalized history sampling
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Nuclear/MEMD.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>

#include "Graniitti/Math/MPolarFourier.h"
#include "Graniitti/Tech/MAux.h"

namespace gra::nuclear {

using gra::aux::indices;

// Prepare a common logarithmic grid spanning the complete photon absorption support
MEMD::MEMD(std::shared_ptr<const MUPC> upc) : upc_(std::move(upc)) {
  if (!upc_) { throw std::invalid_argument("MEMD: missing UPC model"); }
  const auto& param = upc_->Param();
  double minimum = param.convolution.b_max, maximum = minimum;
  std::size_t count = 0;
  for (const auto& leg : indices(mean_)) {
    const auto& breakup = *upc_->Breakup(leg + 1);
    if (upc_->Type(leg + 1) == BeamType::Lepton) {
      throw std::invalid_argument("MEMD: internal EMD requires AA, pA or Ap beams");
    }
    if (upc_->Nucleus(leg + 1) && !breakup.Enabled()) {
      throw std::invalid_argument("MEMD: missing nuclear EMD response");
    }
    if (!breakup.Enabled()) { continue; }
    minimum = std::min(minimum, breakup.Param().profile.b_min);
    maximum = std::max(maximum, breakup.ImpactMax());
    count = std::max(count, breakup.Param().profile.nodes);
    const bool hard = (upc_->Photon(leg + 1) && param.emission[leg] != CoherenceType::Coherent) ||
                      (upc_->Photo(leg + 1) && param.target[leg] != CoherenceType::Coherent);
    absorbed_[leg] = param.emd_condition && !param.neutron[leg].Accept(0) && !hard;
  }
  if (count < 2 || !(minimum > 0.0) || !(maximum > minimum)) {
    throw std::invalid_argument("MEMD: no finite absorption support");
  }
  radius_.resize(count);
  survival_.resize(count);
  proposal_.resize(count);
  log_condition_.assign(count, 0.0);
  const double step = std::log(maximum / minimum) / (count - 1);
  for (const auto& i : indices(radius_)) {
    radius_[i] = minimum * std::exp(i * step);
    // For sampled survival this smooth amplitude only focuses the history proposal
    survival_[i] = upc_->HadronicConvolution() && param.survival == SurvivalType::MCGGCF
                       ? upc_->Glauber()->GGCFAmp(radius_[i]) : upc_->SurvivalAmp(radius_[i]);
    proposal_[i] = survival_[i] * survival_[i];
  }
  for (const auto& leg : indices(mean_)) {
    const auto& breakup = *upc_->Breakup(leg + 1);
    mean_[leg].resize(count);
    log_spectrum_[leg] = MMatrix<double>(count, breakup.Energies().size(), 0.0);
    for (const auto& i : indices(radius_)) {
      mean_[leg][i] = breakup.Mean(radius_[i]);
      const auto spectrum = breakup.Spectrum(radius_[i]);
      for (const auto& j : indices(spectrum)) { log_spectrum_[leg][i][j] = std::log(spectrum[j]); }
      if (absorbed_[leg]) {
        const double probability = -std::expm1(-mean_[leg][i]);
        proposal_[i] *= probability;
        log_condition_[i] += std::log(probability);
      }
    }
  }
  // Focus inclusive histories on absorption while the uniform mixture retains the vacuum support
  if (!absorbed_[0] && !absorbed_[1]) {
    for (const auto& i : indices(proposal_)) {
      proposal_[i] *= -std::expm1(-mean_[0][i] - mean_[1][i]);
    }
  }
  const double norm = gra::Sum(proposal_);
  if (!(norm > 0.0)) { throw std::invalid_argument("MEMD: selected absorption has no support"); }
  for (const auto& i : indices(proposal_)) {
    proposal_[i] = std::isfinite(log_condition_[i])
                      ? (1.0 - param.emd_focus) / count + param.emd_focus * proposal_[i] / norm : 0.0;
  }
  const double total = gra::Sum(proposal_);
  for (double& weight : proposal_) { weight /= total; }
  auto loop = param.loop;
  loop.r_min = loop.r_min > 0.0 ? std::min(loop.r_min, PDG::GeV2fm / maximum) : PDG::GeV2fm / maximum;
  loop.radial_map = math::RadialMap::Log;
  loop.radial_intervals = std::max(loop.radial_intervals, param.convolution.smooth_kt_nodes);
  momentum_ = math::PolarRadialRule(loop).node;
  momentum_.insert(momentum_.begin(), loop.r_min);
  momentum_.push_back(loop.r_max);
  kernel_ = math::LogBesselKernel(radius_, momentum_, PDG::GeV2fm);
}

// Use real independent-exchange amplitudes and trace absorption histories incoherently
// [REFERENCE: arXiv:1311.1938, Eqs. (2.4)-(2.9), independent photon absorption probabilities]
// The proposal is a normalized mixture over impact nodes, not the density at the sampled node
MEMD::State MEMD::Sample(MRandom& random) const {
  State state;
  std::discrete_distribution<std::size_t> impact(proposal_.begin(), proposal_.end());
  const auto selected = impact(random.rng);
  std::vector<double> log_probability(radius_.size(), 0.0);
  bool vacuum = true;
  for (const auto& leg : indices(mean_)) {
    const auto& breakup = *upc_->Breakup(leg + 1);
    const int count = random.PoissonRandom(mean_[leg][selected], absorbed_[leg]);
    vacuum &= count == 0;
    for (const auto& i : indices(log_probability)) {
      log_probability[i] -= mean_[leg][i] + std::lgamma(count + 1.0);
      if (count > 0) { log_probability[i] += count * std::log(mean_[leg][i]); }
    }
    if (count == 0) { continue; }
    std::vector<double> probability(breakup.Energies().size());
    for (const auto& j : indices(probability)) { probability[j] = std::exp(log_spectrum_[leg][selected][j]); }
    std::discrete_distribution<std::size_t> photon(probability.begin(), probability.end());
    for (int n = 0; n < count; ++n) {
      const auto j = photon(random.rng);
      state.excitation[leg] += breakup.SampleTransfer(breakup.Energies()[j], random);
      for (const auto& i : indices(log_probability)) { log_probability[i] += log_spectrum_[leg][i][j]; }
    }
  }
  std::vector<double> terms(radius_.size(), -std::numeric_limits<double>::infinity());
  for (const auto& i : indices(terms)) {
    if (proposal_[i] > 0.0) { terms[i] = std::log(proposal_[i]) + log_probability[i] - log_condition_[i]; }
  }
  const double high = *std::max_element(terms.begin(), terms.end());
  double sum = 0.0;
  for (const double value : terms) { sum += std::exp(value - high); }
  const double log_proposal = high + std::log(sum);
  auto channel = std::make_shared<ExcitationChannel>();
  channel->born = vacuum ? std::exp(-0.5 * log_proposal) : 0.0;
  channel->radius = radius_.back();
  channel->momentum = momentum_;
  channel->impact = radius_;
  channel->amplitude.resize(radius_.size());
  std::vector<double> bare(radius_.size()), screened(radius_.size());
  for (const auto& i : indices(bare)) {
    channel->amplitude[i] = std::exp(0.5 * (log_probability[i] - log_proposal));
    bare[i] = channel->born - channel->amplitude[i];
    screened[i] = channel->born - survival_[i] * channel->amplitude[i];
  }
  // The photon response vanishes beyond its support, leaving only the vacuum channel
  bare.back() = screened.back() = 0.0;
  channel->bare = kernel_ * bare;
  channel->screened = kernel_ * screened;
  state.channel = std::move(channel);
  return state;
}

}  // namespace gra::nuclear
