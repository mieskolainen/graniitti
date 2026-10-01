// Integrated spin densities and quantum correlation metrics
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Spin/MQMetrics.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <ostream>
#include <stdexcept>

#include "Graniitti/Math/MStatistics.h"
#include "Graniitti/Tech/MAux.h"

using gra::aux::indices;

namespace gra::spin {
namespace {

// Compute observables from one integrated spin density
std::array<double, 6> Metrics(const SpinSpec &spec, const MQMetrics::Density &rho) {
  if (spec.type == SpinType::Production) {
    const auto m = DensityMetrics(rho);
    return {m.purity, m.entropy, 0.0, 0.0, 0.0, 0.0};
  }
  const std::size_t first = spec.type == SpinType::ProtonCentral ? 4 : spec.helicity[0].size();
  const auto m = EntanglementMetrics(rho, first, spec.Dimension() / first);
  return {m.pair.purity, m.pair.entropy, m.entropy1, m.entropy2, m.negativity, m.log_negativity};
}

// Encode a complex matrix in row-major real and imaginary pairs
nlohmann::json Encode(const MQMetrics::Density &matrix) {
  auto out = nlohmann::json::array();
  for (const auto value : matrix.Elements()) { out.push_back({value.real(), value.imag()}); }
  return out;
}

// Decode a finite square matrix with its configured physical dimension
MQMetrics::Density Decode(const nlohmann::json &input, const std::size_t dimension) {
  if (!input.is_array() || input.size() != dimension * dimension) {
    throw std::invalid_argument("MQMetrics: invalid serialized matrix dimension");
  }
  MQMetrics::Density out(dimension, dimension, 0.0);
  for (const auto i : indices(input)) {
    if (!input[i].is_array() || input[i].size() != 2) {
      throw std::invalid_argument("MQMetrics: invalid serialized complex entry");
    }
    out.Elements()[i] = {input[i][0].get<double>(), input[i][1].get<double>()};
  }
  if (!out.IsFinite()) { throw std::invalid_argument("MQMetrics: non-finite serialized matrix"); }
  return out;
}

// Encode the physical meaning of a density independently of the amplitude implementation
nlohmann::json EncodeSpec(const SpinSpec &spec) {
  return {{"name", spec.name},   {"type", spec.type == SpinType::Production ? "production" :
                                        spec.type == SpinType::FinalPair ? "final_pair" :
                                        spec.type == SpinType::IntermediatePair ? "intermediate_pair" :
                                        spec.type == SpinType::ProtonCentral ? "proton_central" : "proton_pair"},
          {"pdg", spec.pdg},     {"helicity_x2", spec.helicity},
          {"frame", spec.frame}, {"measure", spec.measure}};
}

}  // namespace

// Validate physical spaces once before constructing worker-local accumulators
void MQMetrics::Configure(std::vector<SpinSpec> specs, const std::size_t groups, const std::size_t max_dimension,
                             std::string note, const std::array<int, 2> protons) {
  if (groups < 2 || max_dimension == 0) { throw std::invalid_argument("MQMetrics: invalid numerical controls"); }
  const std::size_t central_count = specs.size();
  std::vector<std::optional<std::size_t>> proton(central_count);
  if (protons[0] != 0 && protons[1] != 0) {
    for (std::size_t i = 0; i < central_count; ++i) {
      SpinSpec joint = specs[i];
      joint.name += "+pp";
      joint.type = SpinType::ProtonCentral;
      joint.pdg.insert(joint.pdg.begin(), protons.begin(), protons.end());
      joint.helicity.insert(joint.helicity.begin(), 2, std::vector<int>{-1, 1});
      joint.frame = "pp collider helicities; " + joint.frame;
      joint.measure += "; partition (p1,p2) | central, unpolarized incoming beams";
      SpinSpec pair;
      pair.name = specs[i].name + ":pp";
      pair.type = SpinType::ProtonPair;
      pair.pdg = {protons[0], protons[1]};
      pair.helicity = {{-1, 1}, {-1, 1}};
      pair.frame = "pp collider helicities";
      pair.measure = specs[i].measure + "; central spins traced, partition p1 | p2";
      proton[i] = specs.size();
      specs.push_back(std::move(joint));
      specs.push_back(std::move(pair));
    }
  }
  for (const auto i : indices(specs)) {
    const auto &s = specs[i];
    if (s.name.empty() || s.frame.empty() || s.measure.empty() || s.pdg.size() != s.helicity.size() ||
        s.helicity.size() != (s.type == SpinType::ProtonCentral ? 4 : 2) ||
        (s.type == SpinType::Production && s.helicity[1].size() != 1)) {
      throw std::invalid_argument("MQMetrics: invalid spin space " + s.name);
    }
    std::size_t dimension = 1;
    for (const auto &h : s.helicity) {
      if (h.empty() || dimension > max_dimension / h.size()) {
        throw std::invalid_argument("MQMetrics: spin space " + s.name + " exceeds metric_max_dimension or is empty");
      }
      dimension *= h.size();
    }
    for (const auto &h : s.helicity) {
      if (!std::is_sorted(h.begin(), h.end()) || std::adjacent_find(h.begin(), h.end()) != h.end()) {
        throw std::invalid_argument("MQMetrics: helicities must be distinct and ordered");
      }
    }
    for (std::size_t j = 0; j < i; ++j) {
      if (specs[j].name == s.name) { throw std::invalid_argument("MQMetrics: duplicate spin space " + s.name); }
    }
  }
  std::optional<std::size_t> pair;
  for (const auto i : indices(specs)) {
    if (specs[i].type == SpinType::FinalPair) {
      if (pair) { throw std::invalid_argument("MQMetrics: multiple full final spin spaces"); }
      pair = i;
    }
  }
  pair_          = pair;
  specs_         = std::move(specs);
  proton_        = std::move(proton);
  groups_        = groups;
  max_dimension_ = max_dimension;
  note_          = std::move(note);
  Reset();
}

// Find the prepared index of a named physical spin space
std::optional<std::size_t> MQMetrics::Find(const std::string &name) const {
  for (const auto i : indices(specs_)) {
    if (specs_[i].name == name) { return i; }
  }
  return std::nullopt;
}

// Clear sums and event scratch while retaining the immutable definitions
void MQMetrics::Reset() {
  event_.clear();
  for (const auto &spec : specs_) { event_.emplace_back(spec.Dimension(), spec.Dimension(), 0.0); }
  ready_.assign(specs_.size(), false);
  event_failed_ = ready_;
  failures_.assign(specs_.size(), 0);
  sums_.assign(groups_, Group{});
  for (auto &group : sums_) { group.rest = group.last = event_; }
  active_ = false;
}

// Activate density construction only for explicitly selected integration trials
void MQMetrics::Begin(const bool active) {
  active_ = active && !specs_.empty();
  if (active_) { ClearEvent(); }
}

// Discard Born or previous-event densities before collecting the physical result
void MQMetrics::ClearEvent() {
  for (auto &matrix : event_) { std::fill(matrix.Elements().begin(), matrix.Elements().end(), 0.0); }
  std::fill(ready_.begin(), ready_.end(), false);
  std::fill(event_failed_.begin(), event_failed_.end(), false);
}

// Trace unobserved indices after coherent amplitudes have been combined
// [REFERENCE: arXiv:2209.01405, Eq. (3)]
void MQMetrics::Add(const std::size_t index, const Density &amplitude, const double normalization,
                       const std::size_t proton_rows) {
  if (!Active()) { return; }
  try {
    if (amplitude.isEmpty()) { throw std::domain_error("empty spin amplitude"); }
    const std::size_t n    = specs_.at(index).Dimension();
    const auto        rows = amplitude.FlatTensorTranspose(amplitude.size_row() / n, n, amplitude.size_col());
    // rho_ij = sum_a A_ai A_aj^*, preserving the physical complex phase
    const auto rho = (rows.Transpose() * rows.Conj()) * normalization;
    if (!rho.IsFinite() || !(normalization > 0.0)) { throw std::domain_error("invalid spin amplitude"); }
    event_[index] += rho;
    ready_[index] = true;
    if (index < proton_.size() && proton_[index]) {
      if ((proton_rows != 4 && proton_rows != 16) || amplitude.size_row() % (proton_rows * n) != 0) {
        throw std::domain_error("invalid proton helicity layout");
      }
      const std::size_t spectators = amplitude.size_row() / (proton_rows * n);
      const std::size_t unobserved = amplitude.size_col();
      Density joint(4 * spectators * unobserved, 4 * n, 0.0);
      for (std::size_t h = 0; h < proton_rows; ++h) {
        const std::size_t initial = proton_rows == 4 ? h : h / 4;
        const std::size_t final = proton_rows == 4 ? h : h % 4;
        for (std::size_t a = 0; a < spectators; ++a) {
          for (std::size_t m = 0; m < n; ++m) {
            for (std::size_t u = 0; u < unobserved; ++u) {
              joint[(initial * spectators + a) * unobserved + u][final * n + m] =
                  amplitude[(h * spectators + a) * n + m][u];
            }
          }
        }
      }
      const auto density = (joint.Transpose() * joint.Conj()) * normalization;
      const auto j = *proton_[index];
      event_[j] += density;
      event_[j + 1] += density.PartialTrace(4, n, TensorFactor::Second);
      ready_[j] = ready_[j + 1] = true;
    }
  } catch (const std::exception &) {
    if (index >= event_failed_.size()) { throw; }
    event_failed_[index] = true;
    if (index < proton_.size() && proton_[index]) {
      const auto j = *proton_[index];
      event_failed_[j] = event_failed_[j + 1] = true;
    }
  }
}

// Commit one trial using the same proposal correction as the cross section
void MQMetrics::Observe(const double measure, const double log_inverse_density, const std::uint64_t sample) {
  if (!active_) { return; }
  auto      &group = sums_[sample % groups_];
  const bool last  = group.count == 0 || sample > group.last_sample;
  ++group.count;
  if (last) { group.last_sample = sample; }
  for (const auto i : indices(specs_)) {
    auto &density = event_[i];
    bool  failed  = !std::isfinite(measure) || measure < 0.0 || (measure > 0.0 && (!ready_[i] || event_failed_[i]));
    if (measure > 0.0 && !failed) {
      const double log_weight = std::log(measure) + log_inverse_density;
      for (auto &value : density.Elements()) {
        value = {statistics::ImportanceWeight(value.real(), log_weight),
                 statistics::ImportanceWeight(value.imag(), log_weight)};
      }
      failed = !density.IsFinite();
    }
    if (failed) { ++failures_[i]; }
    if (failed || !(measure > 0.0)) { std::fill(density.Elements().begin(), density.Elements().end(), 0.0); }
    if (last) {
      group.rest[i] += group.last[i];
      group.last[i] = density;
    } else {
      group.rest[i] += density;
    }
  }
  active_ = false;
}

// Merge independent workers without transferring active event scratch
void MQMetrics::Merge(const MQMetrics &other) {
  if (specs_ != other.specs_ || proton_ != other.proton_ || groups_ != other.groups_ || max_dimension_ != other.max_dimension_ ||
      note_ != other.note_) {
    throw std::invalid_argument("MQMetrics::Merge: incompatible spin definitions");
  }
  for (const auto g : indices(sums_)) {
    auto       &a    = sums_[g];
    const auto &b    = other.sums_[g];
    const bool  last = b.count > 0 && (a.count == 0 || b.last_sample > a.last_sample);
    for (const auto i : indices(specs_)) {
      a.rest[i] += b.rest[i];
      a.rest[i] += last ? a.last[i] : b.last[i];
    }
    if (last) {
      a.last_sample = b.last_sample;
      a.last        = b.last;
    }
    a.count += b.count;
  }
  for (const auto i : indices(failures_)) { failures_[i] += other.failures_[i]; }
}

// Normalize after integration and estimate metric errors from equal-size groups
std::vector<SpinResult> MQMetrics::Results() const {
  std::vector<SpinResult> out;
  if (!Configured()) { return out; }
  const auto bounds =
      std::minmax_element(sums_.begin(), sums_.end(), [](const auto &a, const auto &b) { return a.count < b.count; });
  const auto min_count = bounds.first->count;
  const auto max_count = bounds.second->count;
  for (const auto i : indices(specs_)) {
    SpinResult result;
    result.spec     = specs_[i];
    result.failures = failures_[i];
    Density              total(specs_[i].Dimension(), specs_[i].Dimension(), 0.0);
    std::vector<Density> blocks;
    for (const auto &group : sums_) {
      result.samples += group.count;
      total += group.rest[i];
      total += group.last[i];
      blocks.push_back(group.rest[i]);
      if (group.count == min_count) { blocks.back() += group.last[i]; }
    }
    result.status      = result.failures > 0 ? SpinStatus::Incomplete : SpinStatus::Empty;
    const double trace = total.Trace().real();
    if (!total.IsFinite() || !std::isfinite(trace) || trace < 0.0) { result.status = SpinStatus::Incomplete; }
    if (std::isfinite(trace) && trace > 0.0 && result.samples > 0 && result.status != SpinStatus::Incomplete) {
      try {
        result.rho      = total / trace;
        result.integral = trace / static_cast<double>(result.samples);
        result.value    = Metrics(specs_[i], result.rho);
        result.status   = SpinStatus::Ready;
        if (min_count > 0 && max_count - min_count <= 1) {
          std::vector<std::vector<double>> replicas;
          for (const auto g : indices(blocks)) {
            Density rho(total.size_row(), total.size_col(), 0.0);
            for (const auto h : indices(blocks)) {
              if (h != g) { rho += blocks[h]; }
            }
            if (!(rho.Trace().real() > 0.0)) { break; }
            const auto values = Metrics(specs_[i], rho);
            replicas.emplace_back(values.begin(), values.end());
          }
          if (replicas.size() == groups_) {
            const auto covariance = statistics::JackknifeCovariance(replicas);
            for (const auto j : indices(result.error)) { result.error[j] = std::sqrt(std::max(0.0, covariance[j][j])); }
            result.errors = true;
          }
        }
      } catch (const std::exception &) { result.status = SpinStatus::Incomplete; }
    }
    out.push_back(std::move(result));
  }
  return out;
}

// Save sufficient statistics together with the exact spin-space definitions
nlohmann::json MQMetrics::Serialize() const {
  nlohmann::json out = {
      {"version", 3}, {"groups", groups_}, {"max_dimension", max_dimension_}, {"note", note_}, {"failures", failures_}};
  out["specs"]   = nlohmann::json::array();
  out["sums"]    = nlohmann::json::array();
  out["results"] = nlohmann::json::array();
  for (const auto &spec : specs_) { out["specs"].push_back(EncodeSpec(spec)); }
  for (const auto &group : sums_) {
    nlohmann::json row = {{"count", group.count},
                          {"last_sample", group.last_sample},
                          {"rest", nlohmann::json::array()},
                          {"last", nlohmann::json::array()}};
    for (const auto &matrix : group.rest) { row["rest"].push_back(Encode(matrix)); }
    for (const auto &matrix : group.last) { row["last"].push_back(Encode(matrix)); }
    out["sums"].push_back(std::move(row));
  }
  for (const auto &r : Results()) {
    out["results"].push_back({{"name", r.spec.name},
                              {"status", r.status == SpinStatus::Ready   ? "ready"
                                         : r.status == SpinStatus::Empty ? "empty"
                                                                         : "incomplete"},
                              {"rho", Encode(r.rho)},
                              {"integral", r.integral},
                              {"metrics", r.value},
                              {"errors", r.errors ? nlohmann::json(r.error) : nlohmann::json(nullptr)},
                              {"samples", r.samples},
                              {"failures", r.failures}});
  }
  return out;
}

// Restore into a temporary accumulator so invalid input cannot alter current sums
void MQMetrics::Deserialize(const nlohmann::json &input) {
  nlohmann::json specs = nlohmann::json::array();
  for (const auto &spec : specs_) { specs.push_back(EncodeSpec(spec)); }
  if (input.at("version").get<int>() != 3 || input.at("specs") != specs ||
      input.at("groups").get<std::size_t>() != groups_ ||
      input.at("max_dimension").get<std::size_t>() != max_dimension_ || input.at("note").get<std::string>() != note_ ||
      input.at("sums").size() != groups_ || input.at("failures").size() != specs_.size()) {
    throw std::invalid_argument("MQMetrics::Deserialize: incompatible spin definitions");
  }
  MQMetrics restored = *this;
  restored.Reset();
  std::uint64_t samples = 0;
  for (const auto g : indices(sums_)) {
    const auto &row   = input.at("sums").at(g);
    auto       &group = restored.sums_[g];
    group.count       = row.at("count").get<std::uint64_t>();
    group.last_sample = row.at("last_sample").get<std::uint64_t>();
    if (!row.at("count").is_number_unsigned() || !row.at("last_sample").is_number_unsigned() ||
        (group.count > 0 && (group.last_sample % groups_ != g || group.count - 1 > group.last_sample / groups_)) ||
        group.count > std::numeric_limits<std::uint64_t>::max() - samples) {
      throw std::invalid_argument("MQMetrics::Deserialize: invalid sample counts");
    }
    samples += group.count;
    if (row.at("rest").size() != specs_.size() || row.at("last").size() != specs_.size()) {
      throw std::invalid_argument("MQMetrics::Deserialize: invalid matrix count");
    }
    for (const auto i : indices(specs_)) {
      group.rest[i] = Decode(row.at("rest").at(i), specs_[i].Dimension());
      group.last[i] = Decode(row.at("last").at(i), specs_[i].Dimension());
      for (const auto *matrix : {&group.rest[i], &group.last[i]}) {
        if (matrix->FrobNorm() > 0.0) { (void)DensityMetrics(*matrix); }
        if (matrix->FrobNorm() > 0.0 && (group.count == 0 || (matrix == &group.rest[i] && group.count == 1))) {
          throw std::invalid_argument("MQMetrics::Deserialize: nonzero density without trials");
        }
      }
    }
  }
  for (const auto i : indices(specs_)) {
    const auto &count = input.at("failures").at(i);
    if (!count.is_number_unsigned() || count.get<std::uint64_t>() > samples) {
      throw std::invalid_argument("MQMetrics::Deserialize: invalid failure count");
    }
    restored.failures_[i] = count.get<std::uint64_t>();
  }
  *this = std::move(restored);
}

// Print a compact spin table followed by definitions and any missing results
void MQMetrics::Print(std::ostream &output) const {
  if (!Configured()) { return; }
  output << '\n' << rang::style::bold << "Integrated spin quantum metrics" << rang::style::reset << "\n\n";
  std::vector<std::vector<std::string>> rows;
  std::vector<std::string>              notes;
  for (const auto &r : Results()) {
    std::string frame = r.spec.frame;
    const std::string proton_frame = "pp collider helicities; ";
    if (frame.starts_with(proton_frame)) { frame = frame.substr(proton_frame.size()); }
    if (frame == "Pair CM helicities") { frame = "Pair CM"; }
    if (r.spec.type == SpinType::ProtonPair) { frame = "Collider"; }
    std::vector<std::string> row = {r.spec.name, frame};
    for (const auto i : indices(r.value)) {
      std::string value = "-";
      if (r.spec.type != SpinType::Production || i < 2) {
        value = r.status == SpinStatus::Ready ? aux::ToString(r.value[i], 3) : "n/a";
        if (r.status == SpinStatus::Ready && r.errors) { value += " +- " + aux::ToString(r.error[i], 3); }
      }
      row.push_back(std::move(value));
    }
    rows.push_back(std::move(row));
    if (r.status != SpinStatus::Ready) {
      notes.push_back(r.spec.name +
                      (r.status == SpinStatus::Incomplete ? ": incomplete" : ": no samples or zero rate") +
                      ", trials=" + std::to_string(r.samples) + ", failures=" + std::to_string(r.failures));
    } else if (!r.errors) {
      notes.push_back(r.spec.name + ": MC error unavailable");
    }
    if (r.spec.type == SpinType::IntermediatePair) { notes.push_back(r.spec.name + ": " + r.spec.measure); }
  }
  if (!rows.empty()) {
    aux::PrintTable(
        {"State", "Frame", "Purity", "S [bits]", "S1 [bits]", "S2 [bits]", "Negativity", "Log negativity [bits]"}, rows,
        output);
    output << "\nNegativity > 0: entangled (zero is inconclusive). S1/S2 measure entanglement only for pure pairs\n";
  }
  if (!note_.empty()) { output << note_ << '\n'; }
  for (const auto &note : notes) { output << note << '\n'; }
}

}  // namespace gra::spin
