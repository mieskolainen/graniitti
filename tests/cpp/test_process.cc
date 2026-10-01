// Process and event record tests
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "support/models_test_support.hh"
#include <future>
#include "Graniitti/Process/MEventRecord.h"
#include "HepMC3/GenVertex.h"

// Check that the collinear unit amplitude survives common helicity normalization
TEST_CASE("DZ flux publishes a normalized unit helicity amplitude",
          "[gra::MProcess][flux][normalization]") {
  gra::PROC_031_QED_YY_DZ_FLUX process;
  gra::LORENTZSCALAR lts;
  lts.beam1.mass = lts.beam2.mass = gra::PDG::mp;
  lts.s = 10000.0;
  for (const double x : {0.01, 0.1, 0.3}) {
    lts.x1 = x;
    lts.x2 = 0.2;
    lts.s_hat = lts.s * lts.x1 * lts.x2;
    // Clear the prior event, then also exercise replacement of an existing bank
    for (const bool previous : {false, true}) {
      lts.hamp.clear();
      if (previous) { lts.hamp.assign(3, std::complex<double>(2.0, -1.0)); }
      const double scalar = process.Amp2(lts);
      REQUIRE(std::isfinite(scalar));
      REQUIRE(scalar > 0.0);
      REQUIRE(lts.hamp.size() == 1);
      CHECK(process.NormalizeAmplitude(lts.hamp).total_weight == Approx(scalar).epsilon(1e-12));
    }
  }
}

namespace {

using Complex = std::complex<double>;
using gra::aux::indices;

// Construct B0 -> D0 pi0, D0 -> K_S pi0 and K_S -> pi+ pi- with physical masses and lifetimes
gra::MDecayBranch MakeFlightTree() {
  const auto pdg = LoadedPDGTable();
  gra::MDecayBranch root;
  root.p = pdg.FindByPDG(511);
  root.p4.SetPxPyPzM(0.2, 0.4, 0.6, root.p.mass);
  root.legs.resize(2);
  auto &child = root.legs[0];
  child.p = pdg.FindByPDG(421);
  root.legs[1].p = pdg.FindByPDG(111);
  child.legs.resize(2);
  auto &grandchild = child.legs[0];
  grandchild.p = pdg.FindByPDG(310);
  child.legs[1].p = pdg.FindByPDG(111);
  grandchild.legs.resize(2);
  grandchild.legs[0].p = pdg.FindByPDG(211);
  grandchild.legs[1].p = pdg.FindByPDG(-211);
  gra::MRandom random;
  random.SetSeed(173);
  for (auto *branch : {&root, &child, &grandchild}) {
    std::vector<gra::M4Vec> daughters;
    const auto weight = gra::kinematics::TwoBodyPhaseSpace(
        branch->p4, branch->p.mass, {branch->legs[0].p.mass, branch->legs[1].p.mass}, daughters, random);
    REQUIRE(weight.Integral() > 0.0);
    REQUIRE(branch->p.tau > 0.0);
    for (const auto &i : indices(branch->legs)) { branch->legs[i].p4 = daughters[i]; }
  }
  return root;
}

// Retain only the defined m=0 GP rows when testing direct ladder kernels
void RetainMZeroGPLadderRows(nlohmann::json &card) {
  for (auto &[exchange, pairs] : card.items()) {
    (void)exchange;
    if (!pairs.is_object()) { continue; }
    for (auto &[pair, block] : pairs.items()) {
      (void)pair;
      if (!block.is_object()) { continue; }
      for (const std::string sector : {"same", "opposite", "self"}) {
        if (!block.contains(sector) || !block.at(sector).is_object()) { continue; }
        for (const std::string field : {"g_ls", "helicity"}) {
          if (!block.at(sector).contains(field)) { continue; }
          auto &rows = block.at(sector)[field];
          if (!rows.is_array() || rows.empty() || !rows.front().is_array() || rows.front().size() != 5) { continue; }
          rows.erase(std::remove_if(rows.begin(), rows.end(),
                                    [](const auto &row) {
                                      return !row.is_array() || row.size() != 5 || row.at(2).template get<int>() != 0;
                                    }),
                     rows.end());
        }
      }
    }
  }
}

// Configure primary and secondary pion exchange pairs for switch tests
void SetSecondaryPionChannels(nlohmann::json &card) {
  card["PARAM_REGGE"]["PARAM_CON"]["MP"]["[211,-211]"] = {{995, 995}, {995, 9915}};
  card["PARAM_REGGE"]["PARAM_CON"]["XP"]["[211,-211]"] = {{995, 995}, {995, 9915}};
  card["PARAM_REGGE"]["PARAM_CON"]["GP"]["[211,-211]"] = {{990, 990}, {990, 9910}};
}

// Compute the independent canonical azimuthal harmonics in row-major spin order
constexpr std::array<int, 16> ReferenceSpinHarmonics() { return {0, -1, 1, 0, 1, 0, 2, 1, -1, -2, 0, -1, 0, -1, 1, 0}; }

// Fill one dense pair operator with a distinct finite complex pattern
void FillReferencePairMatrix(gra::MMatrix<Complex> &matrix, const double scale) {
  for (std::size_t row = 0; row < matrix.size_row(); ++row) {
    for (std::size_t column = 0; column < matrix.size_col(); ++column) {
      matrix(row, column) =
          scale * Complex(0.031 * static_cast<double>(row + 1) + 0.007 * static_cast<double>(column + 1),
                          -0.019 * static_cast<double>(row + 1) + 0.011 * static_cast<double>(column + 1));
    }
  }
}

// Compute the squared distance between two helicity amplitude vectors
double VectorDistance2(const std::vector<Complex> &left, const std::vector<Complex> &right) {
  REQUIRE(left.size() == right.size());
  double out = 0.0;
  for (const auto &i : indices(left)) { out += std::norm(left[i] - right[i]); }
  return out;
}

// Construct sixteen independent nonzero radial pair spin matrices
gra::MEikonalMatrix::PairSpinBank ReferencePairSpinBank(const std::size_t pair_dimension) {
  gra::MEikonalMatrix::PairSpinBank bank;
  for (std::size_t entry = 0; entry < bank.size(); ++entry) {
    bank[entry]       = gra::MMatrix<Complex>(pair_dimension, pair_dimension, 0.0);
    const double sign = entry % 2 == 0 ? 1.0 : -1.0;
    FillReferencePairMatrix(bank[entry], sign * (0.07 + 0.013 * static_cast<double>(entry)));
  }
  return bank;
}

// Construct one full 16-row hard pair source with no spectator indices
gra::ProtonGoodWalkerAmplitude ReferencePairAmplitude(const gra::SoftModelPtr &model, const double scale,
                                                      const bool proton_identity = false) {
  gra::ScreeningMetadata metadata;
  metadata.spin_basis =
      proton_identity ? gra::ScreeningSpinBasis::ProtonIdentity : gra::ScreeningSpinBasis::ProtonHelicity;
  metadata.proton_mode             = gra::ProtonScreeningMode::Elastic;
  metadata.amplitude_normalization = proton_identity ? 1.0 : 0.25;
  metadata.spin_rows               = proton_identity ? 1 : 16;
  metadata.forward_noflip          = proton_identity;
  metadata.amplitude_type          = gra::ScreeningAmplitudeType::GoodWalker;
  metadata.PrepareSpinTransitions();

  const std::size_t              pair_dimension = model->GoodWalker().PairDimension();
  gra::ProtonGoodWalkerComponent component;
  component.coherence_group = 0;
  component.upper_sector    = gra::ProtonGoodWalkerSector::Elastic;
  component.lower_sector    = gra::ProtonGoodWalkerSector::Elastic;
  component.source          = gra::MMatrix<Complex>(metadata.spin_rows, pair_dimension, 0.0);
  for (std::size_t row = 0; row < component.source.size_row(); ++row) {
    for (std::size_t pair = 0; pair < pair_dimension; ++pair) {
      component.source(row, pair) =
          scale * Complex(0.13 * static_cast<double>(row + 1) + 0.017 * static_cast<double>(pair + 1),
                          -0.047 * static_cast<double>(row + 1) + 0.009 * static_cast<double>(pair + 1));
    }
  }
  return {model, model->GoodWalker().ChannelCount(), {std::move(component)}};
}

// Construct one nontrivial azimuthal quadrature node for direct convolution
gra::MEikonal::LoopConst ReferencePairLoop(const std::size_t pair_dimension) {
  gra::MEikonal::LoopConst loop;
  loop.initialized       = true;
  loop.kt2               = {0.25};
  loop.kt_x              = gra::MMatrix<double>(1, 1, 0.0);
  loop.kt_y              = gra::MMatrix<double>(1, 1, 0.0);
  loop.node_weight       = gra::MMatrix<Complex>(1, 1, 0.0);
  loop.kt_x(0, 0)        = 0.3;
  loop.kt_y(0, 0)        = 0.4;
  loop.node_weight(0, 0) = Complex(0.035, -0.012);
  // Cache w exp(i m phi) for the five proton-helicity harmonics m=-2,...,+2
  const Complex z = std::exp(Complex(0.0, std::atan2(loop.kt_y(0, 0), loop.kt_x(0, 0))));
  loop.pair_screening_harmonic_weight.push_back({loop.node_weight(0, 0) * std::conj(z * z),
                                                 loop.node_weight(0, 0) * std::conj(z), loop.node_weight(0, 0),
                                                 loop.node_weight(0, 0) * z, loop.node_weight(0, 0) * z * z});
  loop.pair_screening_spin.push_back(ReferencePairSpinBank(pair_dimension));
  return loop;
}

// Construct one valid quadrature node with a vanishing loop operator
gra::MEikonal::LoopConst ZeroPairLoop(const std::size_t pair_dimension) {
  auto loop = ReferencePairLoop(pair_dimension);
  for (auto &entry : loop.pair_screening_spin.front()) { entry *= 0.0; }
  return loop;
}

// Integrate a specified pair source with the production Good Walker kernel
std::vector<Complex> ConvolvePair(const gra::MEikonal &eikonal, const gra::ProtonGoodWalkerAmplitude &born,
                                  const gra::ProtonGoodWalkerAmplitude &shifted,
                                  const gra::MEikonal::LoopConst &loop, bool identity = false) {
  gra::ScreeningMetadata metadata;
  metadata.spin_basis = identity || born.components.front().source.size_row() == 1
                            ? gra::ScreeningSpinBasis::ProtonIdentity : gra::ScreeningSpinBasis::ProtonHelicity;
  metadata.proton_mode = gra::ProtonScreeningMode::Elastic;
  metadata.amplitude_normalization = metadata.spin_basis == gra::ScreeningSpinBasis::ProtonIdentity ? 1.0 : 0.25;
  metadata.spin_rows = metadata.spin_basis == gra::ScreeningSpinBasis::ProtonIdentity ? 1 : 16;
  metadata.forward_noflip = metadata.spin_rows == 1;
  metadata.amplitude_type = gra::ScreeningAmplitudeType::GoodWalker;
  metadata.PrepareSpinTransitions();
  gra::eikonal::MProtonGoodWalkerScreen screen(born, metadata, loop, eikonal.SoftModelHandle(), eikonal.GetChannelCount());
  for (const auto &i : indices(loop.kt2)) {
    for (std::size_t j = 0; j < loop.node_weight.size_col(); ++j) { screen.Add(i, j, shifted); }
  }
  return screen.Result();
}

// Contract Born and loop sources in four explicit Good Walker indices
std::vector<Complex> DirectPairConvolutionReference(const gra::ProtonGoodWalkerAmplitude &born,
                                                    const gra::ProtonGoodWalkerAmplitude &shifted,
                                                    const gra::MEikonal::LoopConst       &loop) {
  REQUIRE(shifted.components.size() == born.components.size());
  const std::size_t channels  = born.channel_count;
  const auto       &proton    = born.model->GoodWalker().ProtonVector();
  constexpr auto    harmonics = ReferenceSpinHarmonics();

  using CoherentKey = std::tuple<std::size_t, gra::ProtonGoodWalkerSector, gra::ProtonGoodWalkerSector>;
  std::map<CoherentKey, std::vector<Complex>> coherent;
  for (const auto &component : indices(born.components)) {
    const auto &born_component    = born.components[component];
    const auto &shifted_component = shifted.components[component];
    std::vector<Complex> projected(16, 0.0);
    for (std::size_t initial = 0; initial < 4; ++initial) {
      for (std::size_t final = 0; final < 4; ++final) {
        std::vector<Complex> pair(channels * channels, 0.0);
        for (std::size_t a = 0; a < channels; ++a) {
          for (std::size_t b = 0; b < channels; ++b) {
            const std::size_t output_pair = a * channels + b;
            pair[output_pair] =
                born_component.source(gra::spin::PairHelicityTransitionIndex(initial, final), output_pair);
            for (const auto &i : indices(loop.kt2)) {
              for (std::size_t j = 0; j < loop.node_weight.size_col(); ++j) {
                const double azimuth = std::atan2(loop.kt_y(i, j), loop.kt_x(i, j));
                for (std::size_t intermediate = 0; intermediate < 4; ++intermediate) {
                  const std::size_t transition = gra::spin::PairHelicityMatrixIndex(final, intermediate);
                  const auto &matrix = loop.pair_screening_spin[i][transition];
                  const Complex factor = loop.node_weight(i, j) * std::exp(Complex(0.0, harmonics[transition] * azimuth));
                  for (std::size_t c = 0; c < channels; ++c) {
                    for (std::size_t d = 0; d < channels; ++d) {
                      const std::size_t input_pair = c * channels + d;
                      pair[output_pair] += factor * matrix(output_pair, input_pair) * shifted_component.source(
                          gra::spin::PairHelicityTransitionIndex(initial, intermediate), input_pair);
                    }
                  }
                }
              }
            }
          }
        }
        for (std::size_t a = 0; a < channels; ++a) {
          for (std::size_t b = 0; b < channels; ++b) {
            projected[gra::spin::PairHelicityTransitionIndex(initial, final)] +=
                proton[a] * proton[b] * pair[a * channels + b];
          }
        }
      }
    }
    const CoherentKey key = {born_component.coherence_group, born_component.upper_sector, born_component.lower_sector};
    auto              found = coherent.try_emplace(key, projected.size(), Complex(0.0, 0.0)).first;
    gra::AddScaled(found->second, projected, Complex(1.0, 0.0));
  }
  std::vector<Complex> output;
  for (auto &[key, amplitude] : coherent) {
    (void)key;
    output.insert(output.end(), amplitude.begin(), amplitude.end());
  }
  return output;
}

// Prepare the exact local spin kernels needed by one toy pion ladder
void PrepareMultiReggeSpinCacheForTest(gra::LORENTZSCALAR &lts, const gra::ReggeProductionModel model) {
  lts.process.CONT_LADDER_POLE.clear();
  lts.process.REGGE_MODEL = model;
  const int exchange_pdg =
      model == gra::ReggeProductionModel::MP ? 995 : (model == gra::ReggeProductionModel::XP ? 995 : 990);

  for (const auto &first : indices(lts.decaytree)) {
    for (const auto &second : indices(lts.decaytree)) {
      if (lts.decaytree[first].p.pdg != -lts.decaytree[second].p.pdg) { continue; }
      gra::ReggeContinuumPole                         cache;
      const std::array<std::array<std::size_t, 2>, 2> ordering = {std::array<std::size_t, 2>{first, second},
                                                                  std::array<std::size_t, 2>{second, first}};
      for (const auto &leg : indices(ordering)) {
        if (model == gra::ReggeProductionModel::GP) {
          auto &vertex = cache.gp_vertex[leg];
          gra::gpom::InitCrossed(vertex, lts.decaytree[ordering[leg][0]].p.spinX2 / 2.0,
                                 lts.decaytree[ordering[leg][1]].p.spinX2 / 2.0, 0,
                                 "PrepareMultiReggeSpinCacheForTest");
          vertex.coupling_basis = gra::CouplingBasis::Helicity;
          vertex.T[0][0]        = 1.0 + 0.2 * static_cast<double>(leg);
          vertex.T_set[0][0]    = true;
        } else {
          const auto &exchange     = lts.PDG.FindByPDG(exchange_pdg);
          cache.pole_operator[leg] = gra::spin::PoleResidue(gra::spin::PreparePoleLS(
              exchange, lts.decaytree[ordering[leg][0]].p, lts.decaytree[ordering[leg][1]].p, {{2, 0, 1.0}}, 1.0, true,
              false, true, gra::spin::VertexContext::SubTUChannelExchange));
        }
      }
      const gra::ReggeContinuumPoleKey key = {lts.decaytree[first].p.pdg, lts.decaytree[second].p.pdg, exchange_pdg,
                                              exchange_pdg};
      lts.process.CONT_LADDER_POLE.insert_or_assign(key, std::move(cache));
    }
  }
}

// Build one physical alternating pion ladder with all multi-Regge invariants
gra::LORENTZSCALAR MultiReggeLadderLTSForTest(const std::size_t central_multiplicity) {
  REQUIRE((central_multiplicity == 4 || central_multiplicity == 6));
  gra::LORENTZSCALAR lts          = DirectCentralPairLTSForTest(211, -211);
  const gra::M4Vec   central      = lts.pfinal[0];
  const double       central_mass = central.M();
  const double       pion_mass    = lts.PDG.FindByPDG(211).mass;
  const double       energy       = central_mass / central_multiplicity;
  const double       momentum     = std::sqrt(gra::math::pow2(energy) - gra::math::pow2(pion_mass));
  REQUIRE(momentum > 0.0);

  const std::array<std::array<double, 3>, 6> directions = {{{{1.0, 0.0, 0.0}},
                                                            {{-1.0, 0.0, 0.0}},
                                                            {{0.0, 1.0, 0.0}},
                                                            {{0.0, -1.0, 0.0}},
                                                            {{0.0, 0.0, 1.0}},
                                                            {{0.0, 0.0, -1.0}}}};
  lts.decaytree.resize(central_multiplicity);
  for (std::size_t i = 0; i < central_multiplicity; ++i) {
    const int pdg      = (i % 2 == 0) ? 211 : -211;
    lts.decaytree[i].p = lts.PDG.FindByPDG(pdg);
    gra::M4Vec momentum_i(momentum * directions[i][0], momentum * directions[i][1], momentum * directions[i][2],
                          energy);
    gra::kinematics::LorentzBoost(central, central_mass, momentum_i, +1);
    lts.decaytree[i].p4 = momentum_i;
  }
  UpdateToyDurhamDerivedKinematics(lts);

  constexpr std::size_t offset = 3;
  for (std::size_t i = 0; i < central_multiplicity; ++i) {
    const std::size_t a = offset + i;
    const auto       &p = lts.decaytree[i].p4;
    lts.ss[1][a] = lts.ss[a][1] = (lts.pfinal[1] + p).M2();
    lts.ss[2][a] = lts.ss[a][2] = (lts.pfinal[2] + p).M2();
    lts.tt_1[a]                 = (lts.q1 - p).M2();
    lts.tt_2[a]                 = (lts.q2 - p).M2();
    for (std::size_t j = 0; j < central_multiplicity; ++j) {
      const std::size_t b     = offset + j;
      const auto       &other = lts.decaytree[j].p4;
      lts.ss[a][b]            = (p + other).M2();
      lts.tt_xy[a][b]         = (lts.q1 - p - other).M2();
    }
  }
  lts.process.CONT_LADDER_PERMUTATIONS.resize(1);
  for (std::size_t i = 0; i < central_multiplicity; ++i) {
    lts.process.CONT_LADDER_PERMUTATIONS[0].push_back(static_cast<int>(offset + i));
  }

  gra::ScreeningMetadata metadata;
  metadata.spin_basis              = gra::ScreeningSpinBasis::ProtonHelicity;
  metadata.proton_mode             = gra::ProtonScreeningMode::Elastic;
  metadata.amplitude_normalization = 0.25;
  metadata.spin_rows               = 16;
  metadata.forward_noflip          = false;
  metadata.amplitude_type          = gra::ScreeningAmplitudeType::GoodWalker;
  metadata.PrepareSpinTransitions();
  lts.hamp.Configure(metadata);
  PrepareMultiReggeSpinCacheForTest(lts, gra::ReggeProductionModel::MP);
  return lts;
}

// Build one physical multi-Regge event with the real initialized spin caches
gra::LORENTZSCALAR InitializedMultiReggeEventForTest(const std::string &family, const std::size_t central_multiplicity,
                                                     const gra::MModelTunePtr &model_tune) {
  const std::string  pair        = "pi+ pi-";
  const std::string  final_state = central_multiplicity == 4 ? pair + " " + pair : pair + " " + pair + " " + pair;
  ToyHelicityProcess process;
  ConfigureToyProductionProcess(process, family, "CON", final_state);
  process.SetModelTune(model_tune);
  process.SetHelicityConfig(model_tune);
  REQUIRE_NOTHROW(process.InitializeProcessAmplitude());

  gra::LORENTZSCALAR lts = MultiReggeLadderLTSForTest(central_multiplicity);
  lts.process            = process.state.lts.process;
  lts.hamp.Configure(process.state.lts.hamp.metadata);
  return lts;
}

// Build the small physical eikonal profile used by multi-Regge amplitude tests
gra::MEikonal BuildMultiReggeEikonalForTest(const gra::MModelTunePtr &model_tune, const double s) {
  std::filesystem::create_directories(gra::aux::GetBasePath(2) + "/eikonal");
  gra::MEikonal eikonal(model_tune);
  eikonal.Numerics.SetLoopDiscretization(2 - eikonal.Numerics.LOOP.radial_intervals,
                                         2 - eikonal.Numerics.LOOP.azimuth_nodes);
  eikonal.S3Constructor(s, ProtonInitialState(), false, 4, 4);
  return eikonal;
}

// Project an elastic pair source directly onto two physical proton vectors
std::vector<Complex> ElasticCNIBornProjection(const gra::ProtonGoodWalkerAmplitude &state) {
  const auto       &proton   = state.model->GoodWalker().ProtonVector();
  const std::size_t channels = state.channel_count;
  using CoherentKey          = std::tuple<std::size_t, gra::ProtonGoodWalkerSector, gra::ProtonGoodWalkerSector>;
  std::map<CoherentKey, std::vector<Complex>> coherent;
  for (const auto &component : state.components) {
    REQUIRE(component.upper_sector == gra::ProtonGoodWalkerSector::Elastic);
    REQUIRE(component.lower_sector == gra::ProtonGoodWalkerSector::Elastic);
    const CoherentKey key   = {component.coherence_group, component.upper_sector, component.lower_sector};
    auto              found = coherent.try_emplace(key, component.source.size_row(), Complex(0.0, 0.0)).first;
    REQUIRE(found->second.size() == component.source.size_row());
    for (std::size_t row = 0; row < component.source.size_row(); ++row) {
      for (std::size_t upper = 0; upper < channels; ++upper) {
        for (std::size_t lower = 0; lower < channels; ++lower) {
          found->second[row] += proton[upper] * proton[lower] * component.source(row, upper * channels + lower);
        }
      }
    }
  }
  std::vector<Complex> out;
  for (auto &[key, amplitude] : coherent) {
    (void)key;
    out.insert(out.end(), amplitude.begin(), amplitude.end());
  }
  return out;
}

}  // namespace

TEST_CASE("Typed event caches reject stale and cross-scope data", "[gra::MEventCache][cache]") {
  gra::MEventCacheScope              scope;
  gra::MEventCache<int, std::string> cache;
  CHECK_FALSE(scope.Active());
  CHECK_THROWS_AS(cache.Get(scope, 7, "event cache test"), std::logic_error);

  scope.Begin();
  const std::string &stored = cache.Store(scope, 7, "born decay", "event cache test");
  CHECK(stored == "born decay");
  CHECK(cache.Get(scope, 7, "event cache test") == "born decay");
  CHECK_THROWS_AS(cache.Store(scope, 7, "duplicate", "event cache test"), std::logic_error);

  gra::MEventCacheScope copied_scope = scope;
  CHECK_FALSE(copied_scope.Active());
  gra::MEventCache<int, std::string> copied = cache;
  CHECK_THROWS_AS(copied.Get(scope, 7, "event cache test"), std::logic_error);

  gra::MEventCacheScope other;
  other.Begin();
  CHECK_THROWS_AS(cache.Get(other, 7, "event cache test"), std::logic_error);

  scope.Begin();
  CHECK_THROWS_AS(cache.Get(scope, 7, "event cache test"), std::logic_error);
  CHECK(cache.Store(scope, 7, "next event", "event cache test") == "next event");
  scope.Abort();
  CHECK_THROWS_AS(cache.Get(scope, 7, "event cache test"), std::logic_error);
}

TEST_CASE("Regge pair screening matches a direct four-index contraction", "[gra::MProcess][screening][GoodWalker]") {
  const std::vector<std::vector<double>> angle_sets = {{}, {0.21}, {0.21, -0.13, 0.17}};
  for (std::size_t fixture = 0; fixture < angle_sets.size(); ++fixture) {
    const std::size_t channels = fixture + 1;
    CAPTURE(channels);
    const auto eikonal =
        BuildTestEikonal(angle_sets[fixture], "direct_pair_convolution_n" + std::to_string(channels), 4, 4, 0.0);
    const auto model = eikonal.SoftModelHandle();
    REQUIRE(model != nullptr);
    REQUIRE(model->GoodWalker().ChannelCount() == channels);

    const auto born      = ReferencePairAmplitude(model, 1.0);
    const auto shifted   = ReferencePairAmplitude(model, 0.73);
    const auto loop      = ReferencePairLoop(channels * channels);
    const auto reference = DirectPairConvolutionReference(born, shifted, loop);

    const auto screened = ConvolvePair(eikonal, born, shifted, loop);
    const double amp2 = 0.25 * gra::SquaredNorm(screened);
    REQUIRE(screened.size() == reference.size());
    double reference_amp2 = 0.0;
    for (std::size_t row = 0; row < reference.size(); ++row) {
      CAPTURE(row);
      const double tolerance = 3.0e-12 * std::max(1.0, std::abs(reference[row]));
      CHECK(std::abs(screened[row] - reference[row]) <= tolerance);
      reference_amp2 += std::norm(reference[row]);
    }
    reference_amp2 *= 0.25;
    CHECK(amp2 == Approx(reference_amp2).epsilon(2.0e-12));

  }
}

TEST_CASE("Scalar diagonal pair screening matches the dense fallback",
          "[gra::MProcess][screening][GoodWalker][fast_path]") {
  const std::vector<std::vector<double>> angle_sets = {
      {}, {0.21}, {0.21, -0.13, 0.17}, {0.21, -0.13, 0.17, 0.09, -0.07, 0.11}};
  for (std::size_t fixture = 0; fixture < angle_sets.size(); ++fixture) {
    const std::size_t channels = fixture + 1;
    CAPTURE(channels);
    const auto eikonal =
        BuildTestEikonal(angle_sets[fixture], "scalar_diagonal_pair_screening_n" + std::to_string(channels), 4, 4, 0.0);
    const auto model = eikonal.SoftModelHandle();
    REQUIRE(model != nullptr);
    const auto &fast_loop = eikonal.GetLoopConst(eikonal.InitializedMandelstamS());
    REQUIRE(fast_loop.pair_screening_cache_ready);
    REQUIRE(fast_loop.pair_screening_spin_scalar);
    REQUIRE(fast_loop.pair_screening_pair_diagonal);
    REQUIRE(fast_loop.pair_screening_diagonal.size() == fast_loop.kt2.size());
    REQUIRE(fast_loop.pair_screening_harmonic_weight.size() == fast_loop.kt2.size() * fast_loop.node_weight.size_col());

    auto dense_loop                       = fast_loop;
    dense_loop.pair_screening_pair_diagonal = false;
    for (const bool proton_identity : {false, true}) {
      CAPTURE(proton_identity);
      const auto                    born    = ReferencePairAmplitude(model, 1.0, proton_identity);
      const auto                    shifted = ReferencePairAmplitude(model, 0.73, proton_identity);
      const auto fast = ConvolvePair(eikonal, born, shifted, fast_loop, proton_identity);
      const auto dense = ConvolvePair(eikonal, born, shifted, dense_loop, proton_identity);
      const double fast_amp2 = gra::SquaredNorm(fast);
      const double dense_amp2 = gra::SquaredNorm(dense);
      CHECK(fast_amp2 == Approx(dense_amp2).epsilon(3.0e-12));
      RequireVectorNear(fast, dense, 3.0e-12);
    }
  }
}

TEST_CASE("Cached pair harmonics match direct spin rotations", "[gra::MProcess][screening][GoodWalker][harmonics]") {
  const auto eikonal = BuildTestEikonal({0.21}, "cached_pair_harmonic_screening", 4, 4, 0.37);
  const auto model   = eikonal.SoftModelHandle();
  REQUIRE(model != nullptr);
  const auto &cached_loop = eikonal.GetLoopConst(eikonal.InitializedMandelstamS());
  REQUIRE(cached_loop.pair_screening_cache_ready);
  REQUIRE_FALSE(cached_loop.pair_screening_spin_scalar);

  const auto                    born     = ReferencePairAmplitude(model, 1.0);
  const auto                    shifted  = ReferencePairAmplitude(model, 0.73);
  const auto cached = ConvolvePair(eikonal, born, shifted, cached_loop);
  const auto direct = DirectPairConvolutionReference(born, shifted, cached_loop);
  const double cached_amp2 = gra::SquaredNorm(cached);
  const double direct_amp2 = gra::SquaredNorm(direct);
  CHECK(cached_amp2 == Approx(direct_amp2).epsilon(3.0e-12));
  RequireVectorNear(cached, direct, 3.0e-12);
}

TEST_CASE("Regge forward sectors map to explicit Good Walker final bases",
          "[gra::MProcess][gra::MRegge][GoodWalker][validation]") {
  CHECK(gra::SectorFinalBasis(gra::ProtonGoodWalkerSector::Elastic) == gra::GoodWalkerFinalBasis::Proton);
  CHECK(gra::SectorFinalBasis(gra::ProtonGoodWalkerSector::TripleResolved) == gra::GoodWalkerFinalBasis::Excited);
  CHECK(gra::SectorFinalBasis(gra::ProtonGoodWalkerSector::TripleInclusive) == gra::GoodWalkerFinalBasis::Complete);
  CHECK(gra::SectorFinalBasis(gra::ProtonGoodWalkerSector::InelasticEPA) == gra::GoodWalkerFinalBasis::Complete);
  CHECK_THROWS_AS(gra::SectorFinalBasis(static_cast<gra::ProtonGoodWalkerSector>(99)), std::invalid_argument);
}

TEST_CASE("One-channel soft SD and DD preserve the Born norm with a zero loop",
          "[gra::MProcess][gra::MRegge][screening][N1]") {
  const auto eikonal = BuildTestEikonal({}, "soft_me2_zero_loop_n1", 4, 4, 0.0);
  const auto tune    = eikonal.ModelTuneHandle();
  REQUIRE(tune != nullptr);
  const auto model = tune->Soft();
  REQUIRE(model != nullptr);
  REQUIRE(model->GoodWalker().ChannelCount() == 1);

  for (const auto mode : {gra::MReggeInclusive::SD, gra::MReggeInclusive::DD}) {
    CAPTURE(mode);
    gra::LORENTZSCALAR lts = MakeToyCoherentPhotonLTS();
    lts.t                  = -0.17;
    lts.ss[1][1]           = 4.0;
    lts.ss[2][2]           = 9.0;
    lts.excite1            = true;
    lts.excite2            = mode == gra::MReggeInclusive::DD;
    gra::ScreeningMetadata metadata;
    metadata.spin_basis              = gra::ScreeningSpinBasis::ProtonIdentity;
    metadata.proton_mode             = gra::ProtonScreeningMode::TripleRegge;
    metadata.amplitude_normalization = 1.0;
    metadata.spin_rows               = 1;
    metadata.forward_noflip          = false;
    metadata.amplitude_type          = gra::ScreeningAmplitudeType::GoodWalker;
    metadata.PrepareSpinTransitions();
    lts.hamp.Configure(metadata);

    gra::MRegge   regge(lts, tune,
                        gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Soft,
                                                        mode == gra::MReggeInclusive::SD ? "soft_sd" : "soft_dd"));
    const Complex born_scalar = regge.ME2(lts, mode);
    REQUIRE(lts.proton_good_walker.has_value());
    const auto born_pair = *lts.proton_good_walker;
    REQUIRE(born_pair.components.size() == 1);
    REQUIRE(born_pair.components[0].source.size_row() == 1);
    REQUIRE(born_pair.components[0].source.size_col() == 1);

    const auto screened_amp = ConvolvePair(eikonal, born_pair, born_pair, ZeroPairLoop(model->GoodWalker().PairDimension()));
    const double screened = metadata.amplitude_normalization * gra::SquaredNorm(screened_amp);
    const double unscreened = std::norm(born_scalar);
    CHECK(screened == Approx(unscreened).epsilon(2.0e-12));
    REQUIRE(screened_amp.size() == 16);
    std::size_t nonzero_rows = 0;
    for (const auto amplitude : screened_amp) { nonzero_rows += std::abs(amplitude) > 0.0 ? 1 : 0; }
    CHECK(nonzero_rows == 4);
  }
}

TEST_CASE("Shifted Regge nodes retain only canonical pair sources", "[gra::MRegge][screening][source_only]") {
  const auto model = gra::MModelTune::Load(modelfile);

  SECTION("central continuum") {
    gra::LORENTZSCALAR     lts = MakeToyReggeLTSAsymmetric(0.3, 4.2, -4.6);
    gra::ScreeningMetadata metadata;
    metadata.spin_basis              = gra::ScreeningSpinBasis::ProtonHelicity;
    metadata.proton_mode             = gra::ProtonScreeningMode::Elastic;
    metadata.amplitude_normalization = 0.25;
    metadata.spin_rows               = 4;
    metadata.forward_noflip          = true;
    metadata.amplitude_type          = gra::ScreeningAmplitudeType::GoodWalker;
    metadata.PrepareSpinTransitions();
    lts.hamp.Configure(metadata);
    lts.hamp.push_back(Complex(7.0, -3.0));
    lts.screening.active = true;

    gra::MRegge  regge(lts, model,
                       gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::ContinuumTwoBody, "shifted_continuum"));
    const double placeholder = regge.Amp2(lts, gra::ReggeProductionModel::MP, gra::MReggeMode::ContinuumTwoBody);

    REQUIRE(std::isfinite(placeholder));
    CHECK(placeholder == Approx(0.0));
    CHECK(lts.hamp.empty());
    REQUIRE(lts.hamp.metadata.spin_basis == metadata.spin_basis);
    REQUIRE(lts.hamp.metadata.spin_rows == metadata.spin_rows);
    REQUIRE(lts.hamp.metadata.forward_noflip == metadata.forward_noflip);
    REQUIRE(lts.hamp.metadata.amplitude_type == metadata.amplitude_type);
    REQUIRE(lts.hamp.metadata.spin_transition_count == metadata.spin_transition_count);
    REQUIRE(lts.proton_good_walker.has_value());
    CHECK_FALSE(lts.proton_good_walker->components.empty());
  }

  SECTION("single and double diffraction") {
    for (const auto mode : {gra::MReggeInclusive::SD, gra::MReggeInclusive::DD}) {
      CAPTURE(mode);
      gra::LORENTZSCALAR lts = MakeToyCoherentPhotonLTS();
      lts.s                  = 100.0;
      lts.t                  = -0.17;
      lts.ss[1][1]           = 4.0;
      lts.ss[2][2]           = 9.0;
      lts.excite1            = true;
      lts.excite2            = mode == gra::MReggeInclusive::DD;
      gra::ScreeningMetadata metadata;
      metadata.spin_basis              = gra::ScreeningSpinBasis::ProtonIdentity;
      metadata.proton_mode             = gra::ProtonScreeningMode::TripleRegge;
      metadata.amplitude_normalization = 1.0;
      metadata.spin_rows               = 1;
      metadata.forward_noflip          = false;
      metadata.amplitude_type          = gra::ScreeningAmplitudeType::GoodWalker;
      metadata.PrepareSpinTransitions();
      lts.hamp.Configure(metadata);
      lts.hamp.push_back(Complex(-2.0, 5.0));
      lts.screening.active = true;

      gra::MRegge   regge(lts, model,
                          gra::MRegge::ProcessDefinitionFor(
                              gra::MReggeMode::Soft, mode == gra::MReggeInclusive::SD ? "shifted_sd" : "shifted_dd"));
      const Complex placeholder = regge.ME2(lts, mode);

      REQUIRE(std::isfinite(placeholder.real()));
      REQUIRE(std::isfinite(placeholder.imag()));
      CHECK(placeholder == Complex(0.0, 0.0));
      CHECK(lts.hamp.empty());
      REQUIRE(lts.hamp.metadata.spin_basis == metadata.spin_basis);
      REQUIRE(lts.hamp.metadata.spin_rows == metadata.spin_rows);
      REQUIRE(lts.hamp.metadata.amplitude_type == metadata.amplitude_type);
      REQUIRE(lts.hamp.metadata.spin_transition_count == metadata.spin_transition_count);
      REQUIRE(lts.proton_good_walker.has_value());
      CHECK_FALSE(lts.proton_good_walker->components.empty());
    }
  }

  SECTION("direct elastic remains a physical Born projection") {
    gra::LORENTZSCALAR lts = MakeToyCoherentPhotonLTS();
    lts.s                  = 100.0;
    lts.t                  = -0.17;
    gra::ScreeningMetadata metadata;
    metadata.spin_basis              = gra::ScreeningSpinBasis::ProtonHelicity;
    metadata.proton_mode             = gra::ProtonScreeningMode::Elastic;
    metadata.amplitude_normalization = 0.25;
    metadata.spin_rows               = 16;
    metadata.forward_noflip          = false;
    metadata.amplitude_type          = gra::ScreeningAmplitudeType::ElasticCNI;
    metadata.PrepareSpinTransitions();
    lts.hamp.Configure(metadata);
    lts.screening.active = true;

    gra::MRegge   regge(lts, model, gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Soft, "shifted_elastic"));
    const Complex amplitude = regge.ME2(lts, gra::MReggeInclusive::EL);

    REQUIRE(std::isfinite(amplitude.real()));
    REQUIRE(std::isfinite(amplitude.imag()));
    CHECK(std::abs(amplitude) > 0.0);
    CHECK_FALSE(lts.hamp.empty());
    CHECK_FALSE(lts.proton_good_walker.has_value());
  }
}

TEST_CASE("MRegge ME6 and ME8 keep the canonical pair source through screening",
          "[gra::MProcess][gra::MRegge][screening][multiregge]") {
  const gra::LORENTZSCALAR               event      = MultiReggeLadderLTSForTest(4);
  const std::vector<std::vector<double>> angle_sets = {
      {}, {0.21}, {0.21, -0.13, 0.17}, {0.21, -0.13, 0.17, 0.09, -0.07, 0.11}};
  for (std::size_t fixture = 0; fixture < angle_sets.size(); ++fixture) {
    const std::size_t channels = fixture + 1;
    CAPTURE(channels);
    const auto eikonal =
        BuildTestEikonalWithInitialState(angle_sets[fixture], "multiregge_pair_screening_n" + std::to_string(channels),
                                         ProtonInitialState(), 0.25, 1, 4, 4, event.s, 0.0, 0.0, 2, 2);
    const auto tune = eikonal.ModelTuneHandle();
    REQUIRE(tune != nullptr);
    const auto model = tune->Soft();
    REQUIRE(model != nullptr);
    REQUIRE(model->GoodWalker().ChannelCount() == channels);

    for (const std::size_t multiplicity : {4U, 6U}) {
      CAPTURE(multiplicity);
      gra::LORENTZSCALAR lts = MultiReggeLadderLTSForTest(multiplicity);
      gra::MRegge        regge(lts, tune,
                               gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::ContinuumTwoFourSixBody,
                                                          multiplicity == 4 ? "ME6_test" : "ME8_test"));
      const Complex      born_scalar = multiplicity == 4 ? regge.EvalContinuum4(lts, gra::ReggeProductionModel::MP)
                                                         : regge.EvalContinuum6(lts, gra::ReggeProductionModel::MP);
      REQUIRE(std::isfinite(born_scalar.real()));
      REQUIRE(std::isfinite(born_scalar.imag()));
      REQUIRE(std::abs(born_scalar) > 0.0);
      REQUIRE(lts.proton_good_walker.has_value());
      const auto born_pair = *lts.proton_good_walker;
      REQUIRE(born_pair.model == model);
      REQUIRE(born_pair.channel_count == channels);
      REQUIRE(lts.hamp.metadata.spin_basis == gra::ScreeningSpinBasis::ProtonHelicity);
      REQUIRE(lts.hamp.metadata.spin_rows == 16);
      REQUIRE(lts.hamp.metadata.amplitude_type == gra::ScreeningAmplitudeType::GoodWalker);
      for (const auto &component : born_pair.components) {
        REQUIRE(component.source.size_col() == model->GoodWalker().PairDimension());
      }

      const auto born_reference = ElasticCNIBornProjection(born_pair);
      RequireVectorNear(lts.hamp, born_reference, 3.0e-12);
      double born_norm = 0.0;
      for (const auto amplitude : born_reference) { born_norm += std::norm(amplitude); }
      CHECK(std::norm(born_scalar) == Approx(0.25 * born_norm).epsilon(3.0e-12));

      gra::LORENTZSCALAR shifted_lts = MultiReggeLadderLTSForTest(multiplicity);
      shifted_lts.hamp.push_back(Complex(4.0, -6.0));
      shifted_lts.screening.active = true;
      const double placeholder =
          regge.Amp2(shifted_lts, gra::ReggeProductionModel::MP, gra::MReggeMode::ContinuumTwoFourSixBody);
      REQUIRE(std::isfinite(placeholder));
      CHECK(placeholder == Approx(0.0));
      CHECK(shifted_lts.hamp.empty());
      REQUIRE(shifted_lts.hamp.metadata.spin_basis == gra::ScreeningSpinBasis::ProtonHelicity);
      REQUIRE(shifted_lts.hamp.metadata.spin_rows == 16);
      REQUIRE(shifted_lts.hamp.metadata.amplitude_type == gra::ScreeningAmplitudeType::GoodWalker);
      REQUIRE(shifted_lts.proton_good_walker.has_value());
      REQUIRE(shifted_lts.proton_good_walker->channel_count == channels);
      CHECK_FALSE(shifted_lts.proton_good_walker->components.empty());
      for (const auto &component : shifted_lts.proton_good_walker->components) {
        REQUIRE(component.source.size_col() == model->GoodWalker().PairDimension());
      }

      const auto                    loop               = ReferencePairLoop(channels * channels);
      const auto                    screened_reference = DirectPairConvolutionReference(born_pair, born_pair, loop);
      const auto screened_amp = ConvolvePair(eikonal, born_pair, born_pair, loop);
      const double screened = 0.25 * gra::SquaredNorm(screened_amp);
      RequireVectorNear(screened_amp, screened_reference, 4.0e-12);
      double screened_norm = 0.0;
      for (const auto amplitude : screened_reference) { screened_norm += std::norm(amplitude); }
      CHECK(screened == Approx(0.25 * screened_norm).epsilon(4.0e-12));
    }
  }
}

// Compare cached serial ladders with a direct coherent diagram sum
TEST_CASE("Cached serial pion ladders preserve the direct permutation sum", "[gra::MRegge][multiregge][physics]") {
  const gra::LORENTZSCALAR event = MultiReggeLadderLTSForTest(4);
  const auto eikonal = BuildTestEikonalWithInitialState({}, "multiregge_serial_cache_equivalence", ProtonInitialState(),
                                                        0.25, 1, 4, 4, event.s, 0.0, 0.0, 2, 2);
  const auto tune    = eikonal.ModelTuneHandle();
  REQUIRE(tune != nullptr);
  const auto model = tune->Soft();
  REQUIRE(model != nullptr);

  for (const gra::ReggeProductionModel spin :
       {gra::ReggeProductionModel::MP, gra::ReggeProductionModel::XP, gra::ReggeProductionModel::GP}) {
    for (const std::size_t multiplicity : {4U, 6U}) {
      CAPTURE(spin, multiplicity);
      gra::LORENTZSCALAR base = MultiReggeLadderLTSForTest(multiplicity);
      PrepareMultiReggeSpinCacheForTest(base, spin);
      base.process.CONT_LADDER_PERMUTATIONS.clear();
      base.process.MULTIREGGE_TOPOLOGIES = {{static_cast<int>(multiplicity)}};
      gra::MRegge probe(
          base, tune,
          gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::ContinuumTwoFourSixBody,
                                            (multiplicity == 4 ? "ME6" : "ME8") + std::string("_serial_cache_") +
                                                std::to_string(static_cast<int>(spin))));
      const auto param        = ReggeParametersForTest(probe, base);
      const auto permutations = gra::regge::LadderPermutations(base, *param, multiplicity, spin);
      REQUIRE(permutations.size() > 1);

      // Evaluate one explicit permutation bank through the physical serial API
      const auto evaluate = [&](const std::vector<std::vector<int>> &permutation_bank) {
        gra::LORENTZSCALAR lts               = base;
        lts.process.CONT_LADDER_PERMUTATIONS = permutation_bank;
        gra::MRegge regge(lts, tune, gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::ContinuumTwoFourSixBody, "serial_cache_component"));
        if (multiplicity == 4) {
          regge.EvalContinuum4(lts, spin);
        } else {
          regge.EvalContinuum6(lts, spin);
        }
        REQUIRE(lts.proton_good_walker.has_value());
        REQUIRE(lts.proton_good_walker->components.size() == 1);
        return lts;
      };

      const auto            optimized           = evaluate(permutations);
      const auto           &optimized_component = optimized.proton_good_walker->components.front();
      gra::MMatrix<Complex> direct_source(optimized_component.source.size_row(), optimized_component.source.size_col(),
                                          Complex(0.0, 0.0));
      std::vector<Complex>  direct_hamp(optimized.hamp.size(), Complex(0.0, 0.0));
      for (const auto &permutation : permutations) {
        const auto  direct = evaluate({permutation});
        const auto &source = direct.proton_good_walker->components.front().source;
        REQUIRE(source.size_row() == direct_source.size_row());
        REQUIRE(source.size_col() == direct_source.size_col());
        REQUIRE(direct.hamp.size() == direct_hamp.size());
        for (std::size_t row = 0; row < direct_source.size_row(); ++row) {
          for (std::size_t column = 0; column < direct_source.size_col(); ++column) {
            direct_source(row, column) += source(row, column);
          }
        }
        for (const auto &row : indices(direct_hamp)) { direct_hamp[row] += direct.hamp[row]; }
      }
      RequireMatrixNear(optimized_component.source, direct_source, 8.0e-11);
      RequireVectorNear(optimized.hamp, direct_hamp, 8.0e-11);
    }
  }
}

// Preserve numerical failures instead of dropping terms from a coherent topology sum
TEST_CASE("Parallel pion ladders propagate amplitude failures", "[gra::MRegge][multiregge][regression]") {
  const auto event   = MultiReggeLadderLTSForTest(4);
  const auto eikonal = BuildTestEikonalWithInitialState({}, "multiregge_failure", ProtonInitialState(), 0.25, 1, 4, 4,
                                                        event.s, 0.0, 0.0, 2, 2);
  for (const auto spin :
       {gra::ReggeProductionModel::MP, gra::ReggeProductionModel::XP, gra::ReggeProductionModel::GP}) {
    for (const bool bad_attachment : {false, true}) {
      CAPTURE(spin, bad_attachment);
      auto lts = event;
      PrepareMultiReggeSpinCacheForTest(lts, spin);
      lts.process.MULTIREGGE_TOPOLOGIES = {{2, 2}};
      if (!bad_attachment) { lts.process.MULTIREGGE_TOPOLOGIES.push_back({4}); }
      gra::MRegge regge(
          lts, eikonal.ModelTuneHandle(),
          gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::ContinuumTwoFourSixBody, "multiregge_failure"));
      if (bad_attachment) {
        lts.q1.SetPx(std::numeric_limits<double>::quiet_NaN());
      } else {
        lts.q2 = -lts.q1;
      }
      CHECK_THROWS_AS(regge.EvalContinuum4(lts, spin), gra::AmplitudeFailure);
    }
  }
}

TEST_CASE("Selected serial and parallel pion topologies add as one amplitude", "[gra::MRegge][multiregge][coherence]") {
  const gra::LORENTZSCALAR event = MultiReggeLadderLTSForTest(4);
  const auto eikonal = BuildTestEikonalWithInitialState({}, "multiregge_topology_coherence", ProtonInitialState(), 0.25,
                                                        1, 4, 4, event.s, 0.0, 0.0, 2, 2);
  const auto tune    = eikonal.ModelTuneHandle();
  REQUIRE(tune != nullptr);
  const auto model = tune->Soft();
  REQUIRE(model != nullptr);

  for (const gra::ReggeProductionModel spin :
       {gra::ReggeProductionModel::MP, gra::ReggeProductionModel::XP, gra::ReggeProductionModel::GP}) {
    CAPTURE(spin);
    gra::LORENTZSCALAR spin_event = event;
    PrepareMultiReggeSpinCacheForTest(spin_event, spin);

    // Evaluate one explicit topology bank through the physical ME6 API
    const auto evaluate = [&](const std::vector<std::vector<int>> &topologies, const std::string &label) {
      gra::LORENTZSCALAR lts            = spin_event;
      lts.process.MULTIREGGE_TOPOLOGIES = topologies;
      gra::MRegge regge(lts, tune, gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::ContinuumTwoFourSixBody, label));
      regge.EvalContinuum4(lts, spin);
      return lts;
    };

    const std::string suffix   = std::to_string(static_cast<int>(spin));
    const auto        serial   = evaluate({{4}}, "ME6_serial_coherence_" + suffix);
    const auto        parallel = evaluate({{2, 2}}, "ME6_parallel_coherence_" + suffix);
    const auto        combined = evaluate({{4}, {2, 2}}, "ME6_combined_coherence_" + suffix);

    REQUIRE(serial.proton_good_walker.has_value());
    REQUIRE(parallel.proton_good_walker.has_value());
    REQUIRE(combined.proton_good_walker.has_value());
    REQUIRE(serial.proton_good_walker->components.size() == 1);
    REQUIRE(parallel.proton_good_walker->components.size() == 1);
    REQUIRE(combined.proton_good_walker->components.size() == 1);
    const auto &serial_source   = serial.proton_good_walker->components.front().source;
    const auto &parallel_source = parallel.proton_good_walker->components.front().source;
    RequireMatrixNear(combined.proton_good_walker->components.front().source, serial_source + parallel_source, 2.0e-11);

    const std::size_t minus_minus        = gra::spin::BinaryPairHelicityIndexX2(-1, -1);
    const std::size_t minus_minus_noflip = gra::spin::PairHelicityTransitionIndex(minus_minus, minus_minus);
    REQUIRE(minus_minus == 0);
    REQUIRE(minus_minus_noflip == 0);
    REQUIRE(gra::spin::BinaryPairHelicityIndexX2(-1, +1) == 1);
    REQUIRE(gra::spin::BinaryPairHelicityIndexX2(+1, -1) == 2);
    REQUIRE(gra::spin::BinaryPairHelicityIndexX2(+1, +1) == 3);
    for (std::size_t pair = 0; pair < 4; ++pair) {
      const std::size_t noflip               = gra::spin::PairHelicityTransitionIndex(pair, pair);
      double            parallel_noflip_norm = 0.0;
      for (std::size_t column = 0; column < parallel_source.size_col(); ++column) {
        parallel_noflip_norm += std::norm(parallel_source(noflip, column));
      }
      REQUIRE(parallel_noflip_norm > 0.0);
    }

    auto projected_sum = serial.hamp;
    REQUIRE(projected_sum.size() == parallel.hamp.size());
    for (const auto &row : indices(projected_sum)) { projected_sum[row] += parallel.hamp[row]; }
    RequireVectorNear(combined.hamp, projected_sum, 2.0e-11);

    Complex overlap = 0.0;
    for (const auto &row : indices(serial.hamp)) { overlap += std::conj(serial.hamp[row]) * parallel.hamp[row]; }
    const double average       = combined.hamp.metadata.amplitude_normalization;
    const double serial_norm   = average * gra::SquaredNorm(serial.hamp);
    const double parallel_norm = average * gra::SquaredNorm(parallel.hamp);
    const double coherent_norm = average * gra::SquaredNorm(combined.hamp);
    const double interference  = 2.0 * average * std::real(overlap);
    CHECK(coherent_norm == Approx(serial_norm + parallel_norm + interference).epsilon(3.0e-11));
    CHECK(std::abs(interference) > 1.0e-12 * std::max(serial_norm, parallel_norm));
  }
}

// Check that all three six-pion topology orders remain in one coherent source
TEST_CASE("Selected six-pion topology orders add as one amplitude", "[gra::MRegge][multiregge][coherence]") {
  const gra::LORENTZSCALAR event = MultiReggeLadderLTSForTest(6);
  const auto eikonal = BuildTestEikonalWithInitialState({}, "multiregge_six_topology_coherence", ProtonInitialState(),
                                                        0.25, 1, 4, 4, event.s, 0.0, 0.0, 2, 2);
  const auto tune    = eikonal.ModelTuneHandle();
  REQUIRE(tune != nullptr);
  const auto model = tune->Soft();
  REQUIRE(model != nullptr);

  for (const gra::ReggeProductionModel spin :
       {gra::ReggeProductionModel::MP, gra::ReggeProductionModel::XP, gra::ReggeProductionModel::GP}) {
    CAPTURE(spin);
    gra::LORENTZSCALAR spin_event = event;
    PrepareMultiReggeSpinCacheForTest(spin_event, spin);

    // Evaluate one explicit six-pion topology bank through the physical ME8 API
    const auto evaluate = [&](const std::vector<std::vector<int>> &topologies, const std::string &label) {
      gra::LORENTZSCALAR lts            = spin_event;
      lts.process.MULTIREGGE_TOPOLOGIES = topologies;
      gra::MRegge regge(lts, tune, gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::ContinuumTwoFourSixBody, label));
      regge.EvalContinuum6(lts, spin);
      return lts;
    };

    const std::vector<std::vector<std::vector<int>>> topology_bank = {{{6}}, {{4, 2}}, {{2, 2, 2}}};
    std::vector<gra::LORENTZSCALAR>                  isolated;
    isolated.reserve(topology_bank.size());
    for (const auto &order : indices(topology_bank)) {
      CAPTURE(order);
      isolated.push_back(evaluate(topology_bank[order], "ME8_order_" + std::to_string(order + 1)));
    }
    const auto combined =
        evaluate({{6}, {4, 2}, {2, 2, 2}}, "ME8_combined_coherence_" + std::to_string(static_cast<int>(spin)));

    REQUIRE(combined.proton_good_walker.has_value());
    REQUIRE(combined.proton_good_walker->components.size() == 1);
    auto source_sum = isolated.front().proton_good_walker->components.front().source;
    auto projected_sum = isolated.front().hamp;
    REQUIRE_FALSE(projected_sum.empty());
    std::size_t active_orders = 0;
    for (const auto &order : indices(isolated)) {
      CAPTURE(order);
      REQUIRE(isolated[order].proton_good_walker.has_value());
      REQUIRE(isolated[order].proton_good_walker->components.size() == 1);
      const auto &isolated_component = isolated[order].proton_good_walker->components.front();
      REQUIRE(isolated_component.source.IsFinite());
      REQUIRE(gra::AllFinite(isolated[order].hamp));
      active_orders += isolated_component.source.FrobNorm() > 0.0 ? 1 : 0;

      if (order > 0) {
        source_sum += isolated_component.source;
        REQUIRE(projected_sum.size() == isolated[order].hamp.size());
        for (const auto &row : indices(projected_sum)) { projected_sum[row] += isolated[order].hamp[row]; }
      }
    }
    RequireMatrixNear(combined.proton_good_walker->components.front().source, source_sum, 3.0e-11);
    REQUIRE(active_orders > 0);
    REQUIRE(gra::SquaredNorm(combined.hamp) > 0.0);
    RequireVectorNear(combined.hamp, projected_sum, 3.0e-11);
  }
}

TEST_CASE("MP XP and GP evaluate direct four and six pion ladders", "[gra::MProcess][gra::MRegge][multiregge]") {
  const gra::LORENTZSCALAR event = MultiReggeLadderLTSForTest(4);
  const auto eikonal = BuildTestEikonalWithInitialState({}, "multiregge_spin_models", ProtonInitialState(), 0.25, 1, 4,
                                                        4, event.s, 0.0, 0.0, 2, 2);
  const auto tune    = eikonal.ModelTuneHandle();
  REQUIRE(tune != nullptr);
  const auto model = tune->Soft();
  REQUIRE(model != nullptr);

  for (const gra::ReggeProductionModel spin :
       {gra::ReggeProductionModel::MP, gra::ReggeProductionModel::XP, gra::ReggeProductionModel::GP}) {
    for (const std::size_t multiplicity : {4U, 6U}) {
      CAPTURE(spin, multiplicity);
      gra::LORENTZSCALAR lts = MultiReggeLadderLTSForTest(multiplicity);
      PrepareMultiReggeSpinCacheForTest(lts, spin);
      gra::MRegge   regge(lts, tune,
                          gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::ContinuumTwoFourSixBody,
                                                          multiplicity == 4 ? "ME6_spin_test" : "ME8_spin_test"));
      const Complex amplitude = multiplicity == 4 ? regge.EvalContinuum4(lts, spin) : regge.EvalContinuum6(lts, spin);
      CHECK(std::isfinite(amplitude.real()));
      CHECK(std::isfinite(amplitude.imag()));
      CHECK(std::abs(amplitude) > 0.0);

      gra::LORENTZSCALAR changed = MultiReggeLadderLTSForTest(multiplicity);
      PrepareMultiReggeSpinCacheForTest(changed, spin);
      for (auto &[key, cache] : changed.process.CONT_LADDER_POLE) {
        (void)key;
        if (spin == gra::ReggeProductionModel::GP) {
          cache.gp_vertex[0].T[0][0] *= 1.25;
        } else {
          auto vertex = cache.pole_operator[0].Pole();
          vertex.terms.front().coefficient *= 1.25;
          cache.pole_operator[0] = gra::spin::PoleResidue(vertex);
        }
      }
      const Complex modified = multiplicity == 4 ? regge.EvalContinuum4(changed, spin) : regge.EvalContinuum6(changed, spin);
      CHECK(std::abs(modified - amplitude) > 1.0e-10 * std::abs(amplitude));
    }
  }
}

TEST_CASE("GP four and six pion ladders are independent of resonance MMAX",
          "[gra::MProcess][gra::MRegge][multiregge][GP]") {
  const gra::LORENTZSCALAR event = MultiReggeLadderLTSForTest(4);
  const auto eikonal = BuildTestEikonalWithInitialState({}, "multiregge_gp_mmax", ProtonInitialState(), 0.25, 1, 4, 4,
                                                        event.s, 0.0, 0.0, 2, 2);
  const auto tune    = eikonal.ModelTuneHandle();
  REQUIRE(tune != nullptr);
  const auto model = tune->Soft();
  REQUIRE(model != nullptr);

  for (const std::size_t multiplicity : {4U, 6U}) {
    std::vector<Complex> amplitudes;
    for (const int mmax : {1, 4}) {
      gra::LORENTZSCALAR lts = MultiReggeLadderLTSForTest(multiplicity);
      lts.process.MMAX       = mmax;
      PrepareMultiReggeSpinCacheForTest(lts, gra::ReggeProductionModel::GP);
      gra::MRegge regge(lts, tune,
                        gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::ContinuumTwoFourSixBody,
                                                          multiplicity == 4 ? "ME6_gp_mmax" : "ME8_gp_mmax"));
      amplitudes.push_back(multiplicity == 4 ? regge.EvalContinuum4(lts, gra::ReggeProductionModel::GP)
                                             : regge.EvalContinuum6(lts, gra::ReggeProductionModel::GP));
    }
    CAPTURE(multiplicity, amplitudes);
    RequireComplexNear(amplitudes[1], amplitudes[0], 2.0e-12);
  }
}

// Check unsupported transverse photon transport fails before event sampling
TEST_CASE("GP multiparticle photon ladders are rejected during setup",
          "[gra::MProcess][gra::MRegge][multiregge][GP][photon]") {
  ModelParamRestoreGuard restore;
  const auto             tune = WriteModifiedPhotoVMTune("multiregge_gp_photon", [](auto &j) {
    j.at("PARAM_REGGE").at("PARAM_CON").at("GP").at("[211,-211]") = {{22, 22}};
  });
  gra::MODELPARAM             = tune.first;

  ToyHelicityProcess process;
  ConfigureToyProductionProcess(process, "GP", "CON", "pi+ pi- pi+ pi-");
  process.SetModelTune(gra::MModelTune::Load(tune.second));
  REQUIRE_THROWS_AS(process.InitializeProcessAmplitude(), std::invalid_argument);
}

// Check that every multibody topology uses the common continuum pair factor
TEST_CASE("Four and six body ladders use the continuum transfer form factor",
          "[gra::MProcess][gra::MRegge][multiregge][form-factor][physics]") {
  ModelParamRestoreGuard restore;
  const auto narrow_tune = WriteModifiedPhotoVMTune("multiregge_transfer_narrow", {}, {}, [](auto &card) {
    RetainMZeroGPLadderRows(card);
    SetContinuumField(card, "[211,211]", "FF_transfer", {{"type", "power"}, {"norm", "zero"}, {"Lambda2", 0.7}, {"n", 1.0}});
  });
  const auto broad_tune = WriteModifiedPhotoVMTune("multiregge_transfer_broad", {}, {}, [](auto &card) {
    RetainMZeroGPLadderRows(card);
    SetContinuumField(card, "[211,211]", "FF_transfer", {{"type", "power"}, {"norm", "zero"}, {"Lambda2", 2.3}, {"n", 1.0}});
  });
  const auto   narrow_model = gra::MModelTune::Load(narrow_tune.second);
  const auto   broad_model  = gra::MModelTune::Load(broad_tune.second);
  const double s            = MultiReggeLadderLTSForTest(4).s;

  gra::MODELPARAM           = narrow_tune.first;
  const auto narrow_eikonal = BuildMultiReggeEikonalForTest(narrow_model, s);
  gra::MODELPARAM           = broad_tune.first;
  const auto broad_eikonal  = BuildMultiReggeEikonalForTest(broad_model, s);

  const auto evaluate = [](gra::LORENTZSCALAR lts, const gra::MEikonal &eikonal, const gra::ReggeProductionModel model,
                           const gra::regge::Topology &topology) {
    lts.process.MULTIREGGE_TOPOLOGIES = {topology};
    gra::MRegge regge(lts, eikonal.ModelTuneHandle(),
                      gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::ContinuumTwoFourSixBody,
                                                        "continuum_transfer_" + std::to_string(topology.size())));
    if (lts.decaytree.size() == 4) {
      regge.EvalContinuum4(lts, model);
    } else {
      regge.EvalContinuum6(lts, model);
    }
    REQUIRE(gra::AllFinite(lts.hamp));
    REQUIRE(gra::SquaredNorm(lts.hamp) > 0.0);
    return lts.hamp;
  };

  const std::array<std::pair<const char *, gra::ReggeProductionModel>, 3> models = {
      {{"MP", gra::ReggeProductionModel::MP},
       {"XP", gra::ReggeProductionModel::XP},
       {"GP", gra::ReggeProductionModel::GP}}};
  for (const auto &[family, model] : models) {
    CAPTURE(family);
    gra::MODELPARAM        = narrow_tune.first;
    const auto narrow_four = InitializedMultiReggeEventForTest(family, 4, narrow_model);
    const auto narrow_six  = InitializedMultiReggeEventForTest(family, 6, narrow_model);
    gra::MODELPARAM        = broad_tune.first;
    const auto broad_four  = InitializedMultiReggeEventForTest(family, 4, broad_model);
    const auto broad_six   = InitializedMultiReggeEventForTest(family, 6, broad_model);

    for (const auto &[narrow_event, broad_event, topology] :
         {std::tuple{narrow_four, broad_four, gra::regge::Topology{4}},
          std::tuple{narrow_four, broad_four, gra::regge::Topology{2, 2}},
          std::tuple{narrow_six, broad_six, gra::regge::Topology{6}},
          std::tuple{narrow_six, broad_six, gra::regge::Topology{4, 2}},
          std::tuple{narrow_six, broad_six, gra::regge::Topology{2, 2, 2}}}) {
      CAPTURE(topology);
      const auto   narrow     = evaluate(narrow_event, narrow_eikonal, model, topology);
      const auto   broad      = evaluate(broad_event, broad_eikonal, model, topology);
      const double difference = VectorDistance2(narrow, broad);
      const double scale      = std::max(gra::SquaredNorm(narrow), gra::SquaredNorm(broad));
      CHECK(difference > 1.0e-16 * scale);
    }
  }
}

// Check the secondary-Reggeon switch in every serial and parallel amplitude
TEST_CASE("Multi-Regge secondary exchanges enter only connected subladders",
          "[gra::MProcess][gra::MRegge][multiregge][Reggeon][physics]") {
  ModelParamRestoreGuard restore;
  const auto             off_tune = WriteModifiedPhotoVMTune(
                  "multiregge_secondary_off",
                  [](auto &j) {
        j.at("PARAM_REGGE").at("PARAM_CON").at("MULTI").at("secondary_exchanges") = false;
        SetSecondaryPionChannels(j);
      },
                  {}, RetainMZeroGPLadderRows);
  const auto on_tune = WriteModifiedPhotoVMTune(
      "multiregge_secondary_on",
      [](auto &j) {
        j.at("PARAM_REGGE").at("PARAM_CON").at("MULTI").at("secondary_exchanges") = true;
        SetSecondaryPionChannels(j);
      },
      {}, RetainMZeroGPLadderRows);
  const auto   off_model = gra::MModelTune::Load(off_tune.second);
  const auto   on_model  = gra::MModelTune::Load(on_tune.second);
  const double s         = MultiReggeLadderLTSForTest(4).s;

  gra::MODELPARAM        = off_tune.first;
  const auto off_eikonal = BuildMultiReggeEikonalForTest(off_model, s);
  gra::MODELPARAM        = on_tune.first;
  const auto on_eikonal  = BuildMultiReggeEikonalForTest(on_model, s);

  const auto evaluate = [](gra::LORENTZSCALAR lts, const gra::MEikonal &eikonal, const gra::ReggeProductionModel spin,
                           const gra::regge::Topology &topology) {
    lts.process.MULTIREGGE_TOPOLOGIES = {topology};
    gra::MRegge regge(lts, eikonal.ModelTuneHandle(),
                      gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::ContinuumTwoFourSixBody,
                                                        "secondary_exchange_" + std::to_string(topology.size())));
    if (lts.decaytree.size() == 4) {
      regge.EvalContinuum4(lts, spin);
    } else {
      regge.EvalContinuum6(lts, spin);
    }
    REQUIRE(lts.proton_good_walker.has_value());
    REQUIRE(gra::SquaredNorm(lts.hamp) > 0.0);
    return lts;
  };

  const std::array<std::pair<const char *, gra::ReggeProductionModel>, 3> models = {
      {{"MP", gra::ReggeProductionModel::MP},
       {"XP", gra::ReggeProductionModel::XP},
       {"GP", gra::ReggeProductionModel::GP}}};
  for (const auto &[family, spin] : models) {
    CAPTURE(family);
    gra::MODELPARAM     = off_tune.first;
    const auto off_four = InitializedMultiReggeEventForTest(family, 4, off_model);
    const auto off_six  = InitializedMultiReggeEventForTest(family, 6, off_model);
    gra::MODELPARAM     = on_tune.first;
    const auto on_four  = InitializedMultiReggeEventForTest(family, 4, on_model);
    const auto on_six   = InitializedMultiReggeEventForTest(family, 6, on_model);
    REQUIRE(on_four.process.CONT_LADDER_POLE.size() > off_four.process.CONT_LADDER_POLE.size());
    REQUIRE(on_six.process.CONT_LADDER_POLE.size() > off_six.process.CONT_LADDER_POLE.size());

    for (const auto &[off_event, on_event, topology, has_secondary] :
         {std::tuple{off_four, on_four, gra::regge::Topology{4}, true},
          std::tuple{off_six, on_six, gra::regge::Topology{6}, true},
          std::tuple{off_six, on_six, gra::regge::Topology{4, 2}, true},
          std::tuple{off_four, on_four, gra::regge::Topology{2, 2}, false},
          std::tuple{off_six, on_six, gra::regge::Topology{2, 2, 2}, false}}) {
      CAPTURE(topology, has_secondary);
      const auto   off        = evaluate(off_event, off_eikonal, spin, topology);
      const auto   on         = evaluate(on_event, on_eikonal, spin, topology);
      const double difference = VectorDistance2(off.hamp, on.hamp);
      const double scale      = std::max(gra::SquaredNorm(off.hamp), gra::SquaredNorm(on.hamp));
      if (has_secondary) {
        CHECK(difference > 1.0e-16 * scale);
      } else {
        CHECK(difference < 2.0e-22 * scale);
      }
    }
  }
}

TEST_CASE("MSubProc copies lazily rebuild activated process state", "[gra::MSubProc][threading]") {
  ToyHelicityProcess master;
  master.ProcPtr = gra::MSubProc({"X"}, "Q");
  master.ProcPtr.Initialize("X", "ND");
  master.InitializeProcessAmplitude();
  gra::MSubProc &proc = master.ProcPtr;

  gra::LORENTZSCALAR lts;
  REQUIRE(proc.GetBareAmplitude2(lts) == Approx(1.0));

  gra::MSubProc copy(proc);
  REQUIRE(copy.GetBareAmplitude2(lts) == Approx(1.0));

  gra::MSubProc assigned;
  assigned = proc;
  REQUIRE(assigned.GetBareAmplitude2(lts) == Approx(1.0));

  gra::MSubProc invalid({"X"}, "Q");
  invalid.Initialize("X", "missing");
  REQUIRE_THROWS_AS(invalid.Processes(), std::invalid_argument);
}

TEST_CASE("Physical process classes initialize screening metadata once", "[gra::MSubProc][screening]") {
  ToyHelicityProcess elastic;
  elastic.SetInitialState({"p+", "p+"}, {2.5, 2.5});
  elastic.SetScreening(true);
  elastic.state.gcuts.q_t_abs_min          = 1.0e-4;
  elastic.state.gcuts.q_t_abs_max          = 1.0;
  elastic.state.lts.process.FORWARD_NOFLIP = true;
  elastic.ProcPtr                          = gra::MSubProc({"X"}, "Q");
  elastic.ProcPtr.Initialize("X", "EL");
  elastic.InitializeProcessAmplitude();
  const auto &elastic_metadata = elastic.state.lts.hamp.metadata;
  REQUIRE(elastic_metadata.spin_basis == gra::ScreeningSpinBasis::ProtonHelicity);
  REQUIRE(elastic_metadata.proton_mode == gra::ProtonScreeningMode::Elastic);
  REQUIRE(elastic_metadata.amplitude_normalization == Approx(0.25));
  REQUIRE_FALSE(elastic_metadata.forward_noflip);
  REQUIRE(elastic_metadata.spin_rows == 16);
  REQUIRE(elastic_metadata.spin_transition_count == 16);

  ToyHelicityProcess diffraction;
  diffraction.ProcPtr = gra::MSubProc({"X"}, "Q");
  diffraction.ProcPtr.Initialize("X", "SD");
  diffraction.InitializeProcessAmplitude();
  const auto &diffraction_metadata = diffraction.state.lts.hamp.metadata;
  REQUIRE(diffraction_metadata.spin_basis == gra::ScreeningSpinBasis::ProtonIdentity);
  REQUIRE(diffraction_metadata.proton_mode == gra::ProtonScreeningMode::TripleRegge);
  REQUIRE(diffraction_metadata.spin_transition_count == 4);

  ToyHelicityProcess nondiffractive;
  nondiffractive.ProcPtr = gra::MSubProc({"X"}, "Q");
  nondiffractive.ProcPtr.Initialize("X", "ND");
  nondiffractive.InitializeProcessAmplitude();
  REQUIRE(nondiffractive.state.lts.hamp.metadata.spin_basis == gra::ScreeningSpinBasis::Scalar);
  REQUIRE(nondiffractive.state.lts.hamp.metadata.spin_transition_count == 0);
}

TEST_CASE("External eikonals require the exact configured beam state",
          "[gra::MProcess][gra::MEikonal][initial-state]") {
  const auto eikonal = BuildTestEikonal({}, "external_eikonal_compatibility", 4, 4, 0.0);
  const auto tune    = eikonal.ModelTuneHandle();
  REQUIRE(tune != nullptr);
  const auto model = tune->Soft();
  REQUIRE(model != nullptr);

  SECTION("matching proton beams and energy") {
    gra::MFactorized process;
    process.state.lts.PDG = LoadedPDGTable();
    process.SetModelTune(tune);
    process.SetInitialState({"p+", "p+"}, {2.5, 2.5});
    REQUIRE_NOTHROW(process.SetEikonal(eikonal));
  }

  SECTION("Mandelstam s") {
    gra::MFactorized process;
    process.state.lts.PDG = LoadedPDGTable();
    process.SetModelTune(tune);
    process.SetInitialState({"p+", "p+"}, {2.6, 2.5});
    REQUIRE_THROWS_AS(process.SetEikonal(eikonal), std::invalid_argument);
  }

  SECTION("ordered signed beam PDG") {
    gra::MEikonal crossed(tune);
    crossed.S3Constructor(25.0, ProtonAntiprotonInitialState(), false, 4, 4);
    gra::MFactorized process;
    process.state.lts.PDG = LoadedPDGTable();
    process.SetModelTune(tune);
    process.SetInitialState({"p-", "p+"}, {2.5, 2.5});
    REQUIRE_THROWS_AS(process.SetEikonal(crossed), std::invalid_argument);
  }

  SECTION("beam mass") {
    auto altered_state = ProtonInitialState();
    altered_state[0].mass += 1.0e-6;
    gra::MEikonal altered(tune);
    altered.S3Constructor(25.0, altered_state, false, 4, 4);
    gra::MFactorized process;
    process.state.lts.PDG = LoadedPDGTable();
    process.SetModelTune(tune);
    process.SetInitialState({"p+", "p+"}, {2.5, 2.5});
    REQUIRE_THROWS_AS(process.SetEikonal(altered), std::invalid_argument);
  }

  SECTION("beam charge") {
    auto altered_state        = ProtonInitialState();
    altered_state[0].chargeX3 = 0;
    gra::MEikonal altered(tune);
    altered.S3Constructor(25.0, altered_state, false, 4, 4);
    gra::MFactorized process;
    process.state.lts.PDG = LoadedPDGTable();
    process.SetModelTune(tune);
    process.SetInitialState({"p+", "p+"}, {2.5, 2.5});
    REQUIRE_THROWS_AS(process.SetEikonal(altered), std::invalid_argument);
  }
}

TEST_CASE("MSubProc copied wrappers activate independently under parallel calls", "[gra::MSubProc][threading]") {
  ToyHelicityProcess master;
  master.ProcPtr = gra::MSubProc({"X"}, "Q");
  master.ProcPtr.Initialize("X", "ND");
  master.InitializeProcessAmplitude();
  gra::MSubProc &base = master.ProcPtr;

  gra::LORENTZSCALAR lts;
  REQUIRE(base.GetBareAmplitude2(lts) == Approx(1.0));

  constexpr std::size_t nthreads = 8;
  constexpr std::size_t repeats  = 32;

  std::atomic<unsigned int> failures{0};
  std::vector<double>       final_values(nthreads, 0.0);
  std::vector<std::thread>  workers;
  workers.reserve(nthreads);

  for (std::size_t i = 0; i < nthreads; ++i) {
    workers.emplace_back([&, i]() {
      try {
        gra::MSubProc local(base);
        for (std::size_t r = 0; r < repeats; ++r) {
          gra::LORENTZSCALAR local_lts;
          const double       amp2 = local.GetBareAmplitude2(local_lts);
          if (std::abs(amp2 - 1.0) > 1e-12) { ++failures; }
          final_values[i] = amp2;
        }
      } catch (...) { ++failures; }
    });
  }

  for (auto &worker : workers) { worker.join(); }

  REQUIRE(failures.load() == 0);
  for (double value : final_values) { REQUIRE(value == Approx(1.0)); }
}

TEST_CASE("MGraniitti reads the quasielastic t sampling range", "[MGraniitti][GENCUTS][QuasiElastic]") {
  nlohmann::json card = {
      {"GENCUTS", {{"<Q>", {{"Xi", {0.0, 0.1}}, {"t", {0.024, 0.08}}}}}},
  };
  const auto range = ParsedQuasiElasticAbsTRangeForTest(card);
  CHECK(range[0] == Approx(0.024));
  CHECK(range[1] == Approx(0.08));

  card["GENCUTS"]["<Q>"]["t"] = nullptr;
  CHECK_THROWS_AS(ParsedQuasiElasticAbsTRangeForTest(card), std::invalid_argument);

  card["GENCUTS"]["<Q>"].erase("t");
  CHECK_THROWS_AS(ParsedQuasiElasticAbsTRangeForTest(card), std::invalid_argument);

  CHECK_THROWS_AS(ParsedQuasiElasticAbsTRangeForTest(card, "X", "SD", false), std::invalid_argument);
  CHECK_THROWS_AS(ParsedQuasiElasticAbsTRangeForTest(card, "X", "DD", true), std::invalid_argument);

  card["GENCUTS"]["<Q>"]["t"] = {0.0, 0.08};
  CHECK_THROWS_AS(ParsedQuasiElasticAbsTRangeForTest(card), std::invalid_argument);
  CHECK_NOTHROW(ParsedQuasiElasticAbsTRangeForTest(card, "X", "SD", false));
  CHECK_NOTHROW(ParsedQuasiElasticAbsTRangeForTest(card, "X", "DD", true));

  card["GENCUTS"]["<Q>"]["t"] = {0.024, 0.08};
  CHECK_THROWS_AS(ParsedQuasiElasticAbsTRangeForTest(card, "X", "EL", false), std::invalid_argument);

  card["GENCUTS"]["<Q>"]["t"] = {0.024, std::numeric_limits<double>::infinity()};
  CHECK_THROWS_AS(ParsedQuasiElasticAbsTRangeForTest(card), std::invalid_argument);

  card["GENCUTS"]["<Q>"]["t"] = 0.08;
  CHECK_THROWS_AS(ParsedQuasiElasticAbsTRangeForTest(card), std::invalid_argument);

  card["GENCUTS"]["<Q>"]["t"] = {0.08, 0.024};
  CHECK_THROWS_AS(ParsedQuasiElasticAbsTRangeForTest(card), std::invalid_argument);

  nlohmann::json elastic_without_xi = {
      {"GENCUTS", {{"<Q>", {{"t", {0.024, 0.08}}}}}},
  };
  CHECK_NOTHROW(ParsedQuasiElasticAbsTRangeForTest(elastic_without_xi));
  CHECK_THROWS_AS(ParsedQuasiElasticAbsTRangeForTest(elastic_without_xi, "X", "SD", false), std::invalid_argument);
  CHECK_THROWS_AS(ParsedQuasiElasticAbsTRangeForTest(elastic_without_xi, "X", "DD", true), std::invalid_argument);
}

TEST_CASE("MGraniitti reads lowercase VETOCUTS selectors", "[MGraniitti][VETOCUTS]") {
  const nlohmann::json card = {{"VETOCUTS",
                                {{"active", true},
                                 {"0",
                                  {{"Eta", {-5.0, -3.0}},
                                   {"Pt", {0.1, 10.0}},
                                   {"sources", std::vector<std::string>{"forward"}},
                                   {"charge", "charged"}}},
                                 {"1",
                                  {{"Eta", {3.0, 5.0}},
                                   {"Pt", {0.2, 20.0}},
                                   {"sources", std::vector<std::string>{"central"}},
                                   {"charge", "neutral"}}},
                                 {"2",
                                  {{"Eta", {-2.5, 2.5}},
                                   {"Pt", {0.3, 30.0}},
                                   {"sources", std::vector<std::string>{"forward", "central"}},
                                   {"charge", "any"}}},
                                 {"3", {{"Eta", {-1.0, 1.0}}, {"Pt", {0.4, 40.0}}}}}}};

  const gra::VETOCUT veto = ParsedVetoCutsForTest(card);
  CHECK(veto.active);
  REQUIRE(veto.cuts.size() == 4);
  CHECK(veto.cuts[0].source_forward);
  CHECK_FALSE(veto.cuts[0].source_central);
  CHECK(veto.cuts[0].charge == gra::VetoCharge::Charged);
  CHECK_FALSE(veto.cuts[1].source_forward);
  CHECK(veto.cuts[1].source_central);
  CHECK(veto.cuts[1].charge == gra::VetoCharge::Neutral);
  CHECK(veto.cuts[2].source_forward);
  CHECK(veto.cuts[2].source_central);
  CHECK(veto.cuts[2].charge == gra::VetoCharge::Any);
  CHECK(veto.cuts[3].source_forward);
  CHECK(veto.cuts[3].source_central);
  CHECK(veto.cuts[3].charge == gra::VetoCharge::Any);

  auto uppercase_sources_key = card;
  uppercase_sources_key["VETOCUTS"]["0"].erase("sources");
  uppercase_sources_key["VETOCUTS"]["0"]["SOURCES"] = {"forward"};
  CHECK_THROWS_AS(ParsedVetoCutsForTest(uppercase_sources_key), std::invalid_argument);

  auto uppercase_charge_key = card;
  uppercase_charge_key["VETOCUTS"]["0"].erase("charge");
  uppercase_charge_key["VETOCUTS"]["0"]["CHARGE"] = "charged";
  CHECK_THROWS_AS(ParsedVetoCutsForTest(uppercase_charge_key), std::invalid_argument);

  auto uppercase_source_value                        = card;
  uppercase_source_value["VETOCUTS"]["0"]["sources"] = {"FORWARD"};
  CHECK_THROWS_AS(ParsedVetoCutsForTest(uppercase_source_value), std::invalid_argument);

  auto uppercase_charge_value                       = card;
  uppercase_charge_value["VETOCUTS"]["0"]["charge"] = "CHARGED";
  CHECK_THROWS_AS(ParsedVetoCutsForTest(uppercase_charge_value), std::invalid_argument);
}

TEST_CASE("MProcess enforces parsed generic and system fiducial cuts", "[MGraniitti][MProcess][FIDCUTS]") {
  const nlohmann::json card = {{"FIDCUTS",
                                {{"active", true},
                                 {"CENTRAL",
                                  {{"*", {{"Eta", {-3.0, 3.0}}, {"Pt", {15.0, 100000.0}}}},
                                   {"SYSTEM", {{"M", {30.0, 250.0}}, {"Rap", {-3.0, 3.0}}, {"Pt", {0.0, 5.0}}}}}}}}};
  const gra::FIDCUT    cuts = ParsedFiducialCutsForTest(card);
  REQUIRE(cuts.active);
  REQUIRE(cuts.particle_pt_active);
  REQUIRE(cuts.particle_eta_active);
  REQUIRE(cuts.system_M_active);
  REQUIRE(cuts.system_Rap_active);
  REQUIRE(cuts.system_Pt_active);

  const std::vector<gra::MDecayBranch> passing_tree = {MasslessFiducialLeafForTest(21, 15.0, 2.0),
                                                       MasslessFiducialLeafForTest(21, 25.0, -2.0)};
  ToyHelicityProcess                   proc;
  CHECK(proc.CommonCutsForTest(cuts, passing_tree, 100.0, 0.0, 1.0));

  auto low_pt_tree = passing_tree;
  low_pt_tree[0]   = MasslessFiducialLeafForTest(21, 14.9, 2.0);
  CHECK_FALSE(proc.CommonCutsForTest(cuts, low_pt_tree, 100.0, 0.0, 1.0));

  auto forward_tree = passing_tree;
  forward_tree[0]   = MasslessFiducialLeafForTest(21, 20.0, 3.1);
  CHECK_FALSE(proc.CommonCutsForTest(cuts, forward_tree, 100.0, 0.0, 1.0));

  gra::MDecayBranch cascade;
  cascade.legs = passing_tree;
  CHECK(proc.CommonCutsForTest(cuts, {cascade}, 100.0, 0.0, 1.0));
  cascade.legs[1] = MasslessFiducialLeafForTest(21, 14.9, -2.0);
  CHECK_FALSE(proc.CommonCutsForTest(cuts, {cascade}, 100.0, 0.0, 1.0));

  CHECK_FALSE(proc.CommonCutsForTest(cuts, passing_tree, 29.9, 0.0, 1.0));
  CHECK_FALSE(proc.CommonCutsForTest(cuts, passing_tree, 250.1, 0.0, 1.0));
  CHECK_FALSE(proc.CommonCutsForTest(cuts, passing_tree, 100.0, 3.1, 1.0));
  CHECK_FALSE(proc.CommonCutsForTest(cuts, passing_tree, 100.0, 0.0, 5.1));

  gra::FIDCUT zero_mass_cut;
  zero_mass_cut.system_M_active = true;
  zero_mass_cut.M_min           = 0.0;
  zero_mass_cut.M_max           = 1.0;
  gra::LORENTZSCALAR invalid_system;
  invalid_system.m2 = std::numeric_limits<double>::quiet_NaN();
  CHECK_FALSE(zero_mass_cut.PassCentralSystem(invalid_system));

  auto inactive   = cuts;
  inactive.active = false;
  CHECK(proc.CommonCutsForTest(inactive, low_pt_tree, 10.0, 10.0, 10.0));

  gra::FIDCUT rapidity_cut;
  rapidity_cut.active              = true;
  rapidity_cut.particle_rap_active = true;
  rapidity_cut.rap_min             = -1.0;
  rapidity_cut.rap_max             = 1.0;
  CHECK_FALSE(proc.CommonCutsForTest(rapidity_cut, passing_tree));

  gra::M4Vec massive_p4;
  massive_p4.SetPxPyPzM(1.0, 0.0, 0.0, 0.5);
  gra::FIDCUT transverse_energy_cut;
  transverse_energy_cut.active             = true;
  transverse_energy_cut.particle_Et_active = true;
  transverse_energy_cut.Et_min             = 1.2;
  transverse_energy_cut.Et_max             = 2.0;
  CHECK_FALSE(proc.CommonCutsForTest(transverse_energy_cut, {FiducialLeafForTest(211, massive_p4)}));
}

// Apply the jet threshold to each resolved parton while leaving leptons unrestricted
TEST_CASE("Generic jet fiducial cuts accept each physical flavour", "[MGraniitti][MProcess][FIDCUTS]") {
  const nlohmann::json card = {{"FIDCUTS", {{"active", true}, {"CENTRAL", {{"j", {{"Pt", {5.0, 1000.0}}}}}}}}};
  const auto cuts = ParsedFiducialCutsForTest(card);
  const auto muon = FiducialLeafForTest(13, gra::M4Vec(1.0, 0.0, 0.0, 1.1));
  for (int pdg : {1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 21}) {
    const auto hard = FiducialLeafForTest(pdg, gra::M4Vec(10.0, 0.0, 0.0, 10.0));
    const auto soft = FiducialLeafForTest(pdg, gra::M4Vec(1.0, 0.0, 0.0, 1.0));
    CHECK(cuts.PassSelectedParticles({muon, hard}));
    CHECK_FALSE(cuts.PassSelectedParticles({muon, soft}));
    CHECK_FALSE(cuts.PassSelectedParticles({hard, soft}));
  }
  CHECK_FALSE(cuts.PassSelectedParticles({muon}));
  CHECK_FALSE(cuts.PassSelectedParticles({FiducialLeafForTest(22, gra::M4Vec(10.0, 0.0, 0.0, 10.0))}));
}

TEST_CASE("MProcess PDG fiducial cuts select particles and systems", "[MProcess][FIDCUTS]") {
  ToyHelicityProcess proc;

  gra::M4Vec mu_plus_p4;
  mu_plus_p4.SetPxPyPzM(0.1, 0.0, 1.55, 0.105658);
  gra::M4Vec mu_minus_p4;
  mu_minus_p4.SetPxPyPzM(-0.1, 0.0, -1.55, 0.105658);
  const gra::M4Vec gamma_p4(0.30, 0.0, 0.20, std::sqrt(0.30 * 0.30 + 0.20 * 0.20));

  const std::vector<gra::MDecayBranch> tree = {FiducialLeafForTest(22, gamma_p4), FiducialLeafForTest(-13, mu_plus_p4),
                                               FiducialLeafForTest(13, mu_minus_p4)};
  const std::vector<gra::MDecayBranch> permuted_tree = {
      FiducialLeafForTest(13, mu_minus_p4), FiducialLeafForTest(22, gamma_p4), FiducialLeafForTest(-13, mu_plus_p4)};

  gra::FIDCUT    cuts       = BroadFiducialCutsForTest();
  gra::FIDPDGCUT gamma_cut  = PDGFiducialCutForTest({22}, {false});
  gamma_cut.Et              = FiducialRangeForTest(0.20, 100.0);
  gra::FIDPDGCUT dimuon_cut = PDGFiducialCutForTest({13, -13}, {false, false});
  dimuon_cut.M              = FiducialRangeForTest(3.0, 3.2);
  dimuon_cut.Rap            = FiducialRangeForTest(-0.1, 0.1);
  cuts.pdg_cuts             = {gamma_cut, dimuon_cut};

  CHECK(proc.CommonCutsForTest(cuts, tree));
  CHECK(proc.CommonCutsForTest(cuts, permuted_tree));

  auto reversed_selector                = cuts;
  reversed_selector.pdg_cuts[1].pdg     = {-13, 13};
  reversed_selector.pdg_cuts[1].pdg_abs = {false, false};
  CHECK(proc.CommonCutsForTest(reversed_selector, tree));

  auto gamma_fail           = cuts;
  gamma_fail.pdg_cuts[0].Et = FiducialRangeForTest(0.40, 100.0);
  CHECK_FALSE(proc.CommonCutsForTest(gamma_fail, tree));

  auto dimuon_fail          = cuts;
  dimuon_fail.pdg_cuts[1].M = FiducialRangeForTest(3.2, 4.0);
  CHECK_FALSE(proc.CommonCutsForTest(dimuon_fail, tree));

  auto missing_pdg                = cuts;
  missing_pdg.pdg_cuts[0].pdg     = {111};
  missing_pdg.pdg_cuts[0].pdg_abs = {false};
  CHECK_FALSE(proc.CommonCutsForTest(missing_pdg, tree));

  gra::M4Vec                           soft_gamma_p4(0.05, 0.0, 0.0, 0.05);
  const std::vector<gra::MDecayBranch> two_gamma_tree = {
      FiducialLeafForTest(22, gamma_p4), FiducialLeafForTest(22, soft_gamma_p4), FiducialLeafForTest(-13, mu_plus_p4),
      FiducialLeafForTest(13, mu_minus_p4)};
  auto all_gammas_fail           = BroadFiducialCutsForTest();
  all_gammas_fail.pdg_cuts       = {gamma_cut};
  all_gammas_fail.pdg_cuts[0].Et = FiducialRangeForTest(0.20, 100.0);
  CHECK_FALSE(proc.CommonCutsForTest(all_gammas_fail, two_gamma_tree));

  gra::M4Vec pip0;
  pip0.SetPxPyPzM(0.3, 0.0, 0.2, 0.13957);
  gra::M4Vec pim0;
  pim0.SetPxPyPzM(-0.2, 0.0, 0.1, 0.13957);
  gra::M4Vec pip1;
  pip1.SetPxPyPzM(0.1, 0.0, -0.3, 0.13957);
  gra::M4Vec pim1;
  pim1.SetPxPyPzM(-0.1, 0.0, -0.2, 0.13957);
  const std::vector<gra::MDecayBranch> four_pion_tree = {
      FiducialLeafForTest(211, pip0), FiducialLeafForTest(-211, pim0), FiducialLeafForTest(211, pip1),
      FiducialLeafForTest(-211, pim1)};

  auto           pion_pair_selector = BroadFiducialCutsForTest();
  gra::FIDPDGCUT pion_pair_cut      = PDGFiducialCutForTest({211, -211}, {false, false});
  pion_pair_cut.M                   = FiducialRangeForTest(0.0, 0.01);
  pion_pair_selector.pdg_cuts       = {pion_pair_cut};
  CHECK_FALSE(proc.CommonCutsForTest(pion_pair_selector, four_pion_tree));

  auto           incomplete_system_selector = BroadFiducialCutsForTest();
  gra::FIDPDGCUT incomplete_system_cut      = PDGFiducialCutForTest({13, -13, 111}, {false, false, false});
  incomplete_system_cut.M                   = FiducialRangeForTest(0.0, 100.0);
  incomplete_system_selector.pdg_cuts       = {incomplete_system_cut};
  CHECK_FALSE(proc.CommonCutsForTest(incomplete_system_selector, tree));

  auto           four_pion_selector = BroadFiducialCutsForTest();
  gra::FIDPDGCUT four_pion_cut      = PDGFiducialCutForTest({211, -211, 211, -211}, {false, false, false, false});
  four_pion_cut.M                   = FiducialRangeForTest(0.0, 0.01);
  four_pion_selector.pdg_cuts       = {four_pion_cut};
  CHECK_FALSE(proc.CommonCutsForTest(four_pion_selector, four_pion_tree));

  auto           abs_dimuon_selector = BroadFiducialCutsForTest();
  gra::FIDPDGCUT abs_dimuon_cut      = PDGFiducialCutForTest({13, 13}, {true, true});
  abs_dimuon_cut.M                   = FiducialRangeForTest(3.2, 4.0);
  abs_dimuon_selector.pdg_cuts       = {abs_dimuon_cut};
  CHECK_FALSE(proc.CommonCutsForTest(abs_dimuon_selector, tree));

  abs_dimuon_selector.pdg_cuts[0].M = FiducialRangeForTest(0.0, 100.0);
  CHECK(proc.CommonCutsForTest(abs_dimuon_selector, tree));
  auto excess_muon_tree = tree;
  excess_muon_tree.push_back(FiducialLeafForTest(13, mu_minus_p4));
  CHECK_FALSE(proc.CommonCutsForTest(abs_dimuon_selector, excess_muon_tree));

  auto           abs_muon_particle_selector = BroadFiducialCutsForTest();
  gra::FIDPDGCUT abs_muon_particle_cut      = PDGFiducialCutForTest({13}, {true});
  abs_muon_particle_cut.Pt                  = FiducialRangeForTest(0.11, 100.0);
  abs_muon_particle_selector.pdg_cuts       = {abs_muon_particle_cut};
  CHECK_FALSE(proc.CommonCutsForTest(abs_muon_particle_selector, tree));
}

TEST_CASE("MGraniitti accepts sparse FIDCUTS blocks and observables", "[MGraniitti][FIDCUTS]") {
  const nlohmann::json sparse_card = {{"FIDCUTS",
                                       {
                                           {"active", true},
                                           {"CENTRAL",
                                            {{"*", nlohmann::json::object()},
                                             {"SYSTEM", {{"M", {3.0, 4.0}}}},
                                             {"[22]", {{"Et", {0.2, 100000.0}}}},
                                             {"[13,-13]", {{"M", {0.0, 100.0}}, {"Rap", {-5.0, 5.0}}}}}},
                                           {"FORWARD", nlohmann::json::object()},
                                       }}};

  const gra::FIDCUT sparse = ParsedFiducialCutsForTest(sparse_card);
  CHECK(sparse.active);
  CHECK_FALSE(sparse.particle_eta_active);
  CHECK(sparse.eta_min == Approx(-30.0));
  CHECK(sparse.eta_max == Approx(30.0));
  CHECK(sparse.system_M_active);
  CHECK(sparse.M_min == Approx(3.0));
  CHECK(sparse.M_max == Approx(4.0));
  CHECK_FALSE(sparse.system_Rap_active);
  CHECK(sparse.Y_min == Approx(-30.0));
  CHECK(sparse.Y_max == Approx(30.0));
  CHECK_FALSE(sparse.forward_t_active);
  CHECK(sparse.forward_t_min == Approx(0.0));
  CHECK(sparse.forward_t_max == Approx(1000000.0));
  CHECK_FALSE(sparse.forward_xi_active);
  CHECK(sparse.forward_xi_min == Approx(0.0));
  CHECK(sparse.forward_xi_max == Approx(1.0));
  REQUIRE(sparse.pdg_cuts.size() == 2);
  const auto find_cut = [](const gra::FIDCUT &cuts, const std::vector<int> &pdg, const std::vector<bool> &pdg_abs) {
    for (const auto &cut : cuts.pdg_cuts) {
      if (cut.pdg == pdg && cut.pdg_abs == pdg_abs) { return &cut; }
    }
    return static_cast<const gra::FIDPDGCUT *>(nullptr);
  };
  const gra::FIDPDGCUT *gamma_cut  = find_cut(sparse, {22}, {false});
  const gra::FIDPDGCUT *dimuon_cut = find_cut(sparse, {13, -13}, {false, false});
  REQUIRE(gamma_cut != nullptr);
  CHECK(gamma_cut->Et.active);
  CHECK_FALSE(gamma_cut->Eta.active);
  REQUIRE(dimuon_cut != nullptr);
  CHECK(dimuon_cut->M.active);
  CHECK(dimuon_cut->Rap.active);
  CHECK_FALSE(dimuon_cut->Pt.active);

  const nlohmann::json pdg_only_card = {
      {"FIDCUTS", {{"active", true}, {"CENTRAL", {{"[22]", {{"Pt", {0.1, 10.0}}}}}}}}};
  const gra::FIDCUT pdg_only = ParsedFiducialCutsForTest(pdg_only_card);
  CHECK(pdg_only.active);
  CHECK_FALSE(pdg_only.particle_pt_active);
  CHECK(pdg_only.pt_min == Approx(0.0));
  CHECK(pdg_only.pt_max == Approx(1000000.0));
  REQUIRE(pdg_only.pdg_cuts.size() == 1);
  CHECK(pdg_only.pdg_cuts[0].pdg == std::vector<int>{22});
  CHECK(pdg_only.pdg_cuts[0].pdg_abs == std::vector<bool>{false});
  CHECK(pdg_only.pdg_cuts[0].Pt.active);

  const nlohmann::json pdg_scalar_card = {
      {"FIDCUTS", {{"active", true}, {"CENTRAL", {{"22", {{"Pt", {0.1, 10.0}}}}}}}}};
  const gra::FIDCUT pdg_scalar = ParsedFiducialCutsForTest(pdg_scalar_card);
  REQUIRE(pdg_scalar.pdg_cuts.size() == 1);
  CHECK(pdg_scalar.pdg_cuts[0].pdg == std::vector<int>{22});
  CHECK(pdg_scalar.pdg_cuts[0].pdg_abs == std::vector<bool>{false});

  const nlohmann::json pdg_abs_card = {
      {"FIDCUTS", {{"active", true}, {"CENTRAL", {{"ABS(13)", {{"Pt", {0.1, 10.0}}}}}}}}};
  const gra::FIDCUT pdg_abs = ParsedFiducialCutsForTest(pdg_abs_card);
  REQUIRE(pdg_abs.pdg_cuts.size() == 1);
  CHECK(pdg_abs.pdg_cuts[0].pdg == std::vector<int>{13});
  CHECK(pdg_abs.pdg_cuts[0].pdg_abs == std::vector<bool>{true});

  const nlohmann::json no_pdg_card = {{"FIDCUTS", {{"active", true}, {"CENTRAL", {{"*", {{"Pt", {0.2, 10.0}}}}}}}}};
  const gra::FIDCUT    no_pdg      = ParsedFiducialCutsForTest(no_pdg_card);
  CHECK(no_pdg.active);
  CHECK(no_pdg.particle_pt_active);
  CHECK(no_pdg.pt_min == Approx(0.2));
  CHECK(no_pdg.pt_max == Approx(10.0));
  CHECK(no_pdg.pdg_cuts.empty());

  gra::M4Vec high_eta_p4;
  high_eta_p4.SetPxPyPzM(1.0, 0.0, 100.0, 0.13957);
  const std::vector<gra::MDecayBranch> high_eta_tree = {FiducialLeafForTest(211, high_eta_p4)};
  ToyHelicityProcess                   proc;
  CHECK(proc.CommonCutsForTest(no_pdg, high_eta_tree));

  auto eta_restricted                = no_pdg;
  eta_restricted.particle_eta_active = true;
  eta_restricted.eta_min             = -1.0;
  eta_restricted.eta_max             = 1.0;
  CHECK_FALSE(proc.CommonCutsForTest(eta_restricted, high_eta_tree));

  const nlohmann::json system_pt_card = {
      {"FIDCUTS", {{"active", true}, {"CENTRAL", {{"SYSTEM", {{"Pt", {0.0, 1.0}}}}}}}}};
  const gra::FIDCUT system_pt = ParsedFiducialCutsForTest(system_pt_card);
  CHECK(system_pt.system_Pt_active);
  CHECK_FALSE(system_pt.system_M_active);
  CHECK(proc.CommonCutsForTest(system_pt, high_eta_tree));

  auto mass_restricted            = system_pt;
  mass_restricted.system_M_active = true;
  mass_restricted.M_min           = 0.0;
  mass_restricted.M_max           = 1.0;
  CHECK_FALSE(proc.CommonCutsForTest(mass_restricted, high_eta_tree));

  const nlohmann::json empty_pdg_card = {{"FIDCUTS", {{"active", true}, {"CENTRAL", nlohmann::json::object()}}}};
  const gra::FIDCUT    empty_pdg      = ParsedFiducialCutsForTest(empty_pdg_card);
  CHECK(empty_pdg.active);
  CHECK(empty_pdg.pdg_cuts.empty());

  const nlohmann::json empty_pdg_entry_card = {
      {"FIDCUTS", {{"active", true}, {"CENTRAL", {{"[22]", nlohmann::json::object()}}}}}};
  const gra::FIDCUT empty_pdg_entry = ParsedFiducialCutsForTest(empty_pdg_entry_card);
  CHECK(empty_pdg_entry.active);
  CHECK(empty_pdg_entry.pdg_cuts.empty());

  const nlohmann::json negative_pt_card = {
      {"FIDCUTS", {{"active", true}, {"CENTRAL", {{"*", {{"Pt", {-1.0, 10.0}}}}}}}}};
  CHECK_THROWS_AS(ParsedFiducialCutsForTest(negative_pt_card), std::invalid_argument);

  const nlohmann::json negative_system_mass_card = {
      {"FIDCUTS", {{"active", true}, {"CENTRAL", {{"SYSTEM", {{"M", {-1.0, 10.0}}}}}}}}};
  CHECK_THROWS_AS(ParsedFiducialCutsForTest(negative_system_mass_card), std::invalid_argument);

  const nlohmann::json negative_pdg_et_card = {
      {"FIDCUTS", {{"active", true}, {"CENTRAL", {{"22", {{"Et", {-1.0, 10.0}}}}}}}}};
  CHECK_THROWS_AS(ParsedFiducialCutsForTest(negative_pdg_et_card), std::invalid_argument);

  const nlohmann::json negative_forward_t_card = {{"FIDCUTS", {{"active", true}, {"FORWARD", {{"t", {-1.0, 10.0}}}}}}};
  CHECK_THROWS_AS(ParsedFiducialCutsForTest(negative_forward_t_card), std::invalid_argument);

  const nlohmann::json invalid_forward_xi_card = {
      {"FIDCUTS", {{"active", true}, {"FORWARD", {{"Xi", {-0.01, 1.01}}}}}}};
  CHECK_THROWS_AS(ParsedFiducialCutsForTest(invalid_forward_xi_card), std::invalid_argument);

  const double         nan                = std::numeric_limits<double>::quiet_NaN();
  const nlohmann::json nonfinite_eta_card = {
      {"FIDCUTS", {{"active", true}, {"CENTRAL", {{"*", {{"Eta", {nan, 10.0}}}}}}}}};
  CHECK_THROWS_AS(ParsedFiducialCutsForTest(nonfinite_eta_card), std::invalid_argument);

  const nlohmann::json non_integer_pdg_card = {
      {"FIDCUTS", {{"active", true}, {"CENTRAL", {{"[22.5]", {{"Pt", {0.1, 10.0}}}}}}}}};
  CHECK_THROWS_AS(ParsedFiducialCutsForTest(non_integer_pdg_card), std::invalid_argument);

  const nlohmann::json unknown_particle_observable = {
      {"FIDCUTS", {{"active", true}, {"CENTRAL", {{"*", {{"Pseudorapidity", {-2.5, 2.5}}}}}}}}};
  CHECK_THROWS_AS(ParsedFiducialCutsForTest(unknown_particle_observable), std::invalid_argument);

  const nlohmann::json unknown_pdg_observable = {
      {"FIDCUTS", {{"active", true}, {"CENTRAL", {{"[13,-13]", {{"Mass", {3.0, 4.0}}}}}}}}};
  CHECK_THROWS_AS(ParsedFiducialCutsForTest(unknown_pdg_observable), std::invalid_argument);

  const nlohmann::json unknown_forward_observable = {
      {"FIDCUTS", {{"active", true}, {"FORWARD", {{"phi", {0.0, 180.0}}}}}}};
  CHECK_THROWS_AS(ParsedFiducialCutsForTest(unknown_forward_observable), std::invalid_argument);

  const nlohmann::json invalid_abs_pdg_card = {
      {"FIDCUTS", {{"active", true}, {"CENTRAL", {{"ABS(-2147483648)", {{"Pt", {0.1, 10.0}}}}}}}}};
  CHECK_THROWS_AS(ParsedFiducialCutsForTest(invalid_abs_pdg_card), std::invalid_argument);
}

TEST_CASE("MProcess rejects unsupported fiducial configurations", "[MProcess][FIDCUTS]") {
  ToyHelicityProcess proc;

  gra::FIDCUT central = BroadFiducialCutsForTest();
  CHECK(gra::IsKnownUserCut(0));
  CHECK(gra::IsKnownUserCut(-3));
  CHECK(gra::IsKnownUserCut(1792394010));
  CHECK_NOTHROW(proc.ValidateFiducialCuts(central, 1792394010));
  CHECK_THROWS_AS(proc.ValidateFiducialCuts(central, 3000000000LL), std::invalid_argument);
  proc.ConfigureFiducialValidationForTest("Q", "EL", 0);
  CHECK_THROWS_AS(proc.ValidateFiducialCuts(central, 0), std::invalid_argument);

  gra::FIDCUT forward;
  forward.active           = true;
  forward.forward_t_active = true;
  proc.ConfigureFiducialValidationForTest("P", "CON", 0);
  CHECK_THROWS_AS(proc.ValidateFiducialCuts(forward, 0), std::invalid_argument);

  gra::FIDCUT forward_mass;
  forward_mass.active           = true;
  forward_mass.forward_M_active = true;
  proc.ConfigureFiducialValidationForTest("Q", "EL", 0);
  CHECK_THROWS_AS(proc.ValidateFiducialCuts(forward_mass, 0), std::invalid_argument);
  proc.ConfigureFiducialValidationForTest("F", "CON", 0);
  CHECK_THROWS_AS(proc.ValidateFiducialCuts(forward_mass, 0), std::invalid_argument);
  proc.ConfigureFiducialValidationForTest("F", "CON", 1);
  CHECK_NOTHROW(proc.ValidateFiducialCuts(forward_mass, 0));

  gra::FIDCUT inactive;
  CHECK_THROWS_AS(proc.ValidateFiducialCuts(inactive, 7), std::invalid_argument);
}

TEST_CASE("FORWARD dPhi cuts use degrees for elastic and dissociative systems", "[MGraniitti][MProcess][FIDCUTS]") {
  const nlohmann::json below_card = {{"FIDCUTS", {{"active", true}, {"FORWARD", {{"dPhi", {0.0, 90.0}}}}}}};
  const gra::FIDCUT    below      = ParsedFiducialCutsForTest(below_card);
  CHECK(below.forward_dPhi_active);
  CHECK(below.forward_dPhi_min == Approx(0.0));
  CHECK(below.forward_dPhi_max == Approx(90.0));
  CHECK(below.HasForwardCuts());

  ToyHelicityProcess proc;
  CHECK(proc.CommonForwardCutsForTest(below, MakeToyProductionLTSDphi(0.25 * gra::math::PI)));
  CHECK(proc.CommonForwardCutsForTest(below, MakeToyProductionLTSDphi(0.5 * gra::math::PI)));
  CHECK_FALSE(proc.CommonForwardCutsForTest(below, MakeToyProductionLTSDphi(0.75 * gra::math::PI)));

  gra::LORENTZSCALAR single_dissociation = MakeToyProductionLTSDphi(0.25 * gra::math::PI);
  single_dissociation.excite1            = true;
  single_dissociation.pfinal[1].SetE(single_dissociation.pfinal[1].E() + 1.0);
  CHECK(proc.CommonForwardCutsForTest(below, single_dissociation));

  gra::LORENTZSCALAR double_dissociation = MakeToyProductionLTSDphi(0.75 * gra::math::PI);
  double_dissociation.excite1            = true;
  double_dissociation.excite2            = true;
  double_dissociation.pfinal[1].SetE(double_dissociation.pfinal[1].E() + 1.0);
  double_dissociation.pfinal[2].SetE(double_dissociation.pfinal[2].E() + 1.0);
  CHECK_FALSE(proc.CommonForwardCutsForTest(below, double_dissociation));

  gra::LORENTZSCALAR incomplete;
  incomplete.pfinal.resize(2);
  CHECK_FALSE(proc.CommonForwardCutsForTest(below, incomplete));

  const gra::LORENTZSCALAR transfer_event = MakeToyProductionLTSDphi(0.25 * gra::math::PI);
  const double             maximum_abs_t  = std::max(std::abs(transfer_event.t1), std::abs(transfer_event.t2));
  gra::FIDCUT              transfer_cut;
  transfer_cut.active           = true;
  transfer_cut.forward_t_active = true;
  transfer_cut.forward_t_min    = 0.0;
  transfer_cut.forward_t_max    = maximum_abs_t;
  CHECK(proc.CommonForwardCutsForTest(transfer_cut, transfer_event));
  transfer_cut.forward_t_max = 0.5 * maximum_abs_t;
  CHECK_FALSE(proc.CommonForwardCutsForTest(transfer_cut, transfer_event));
  transfer_cut.forward_t_max            = maximum_abs_t;
  gra::LORENTZSCALAR nonfinite_transfer = transfer_event;
  nonfinite_transfer.t1                 = std::numeric_limits<double>::quiet_NaN();
  CHECK_FALSE(proc.CommonForwardCutsForTest(transfer_cut, nonfinite_transfer));

  gra::LORENTZSCALAR mass_event = MakeToyProductionLTSDphi(0.25 * gra::math::PI);
  mass_event.excite1            = true;
  mass_event.pfinal[1].SetE(mass_event.pfinal[1].E() + 1.0);
  const double excited_mass = mass_event.pfinal[1].M();
  gra::FIDCUT  mass_cut;
  mass_cut.active           = true;
  mass_cut.forward_M_active = true;
  mass_cut.forward_M_min    = excited_mass - 0.1;
  mass_cut.forward_M_max    = excited_mass + 0.1;
  CHECK(proc.CommonForwardCutsForTest(mass_cut, mass_event));
  mass_cut.forward_M_min = excited_mass + 0.1;
  mass_cut.forward_M_max = excited_mass + 1.0;
  CHECK_FALSE(proc.CommonForwardCutsForTest(mass_cut, mass_event));
  mass_cut.forward_M_min            = excited_mass - 0.1;
  mass_cut.forward_M_max            = excited_mass + 0.1;
  gra::LORENTZSCALAR nonfinite_mass = mass_event;
  nonfinite_mass.pfinal[1].SetE(std::numeric_limits<double>::quiet_NaN());
  CHECK_FALSE(proc.CommonForwardCutsForTest(mass_cut, nonfinite_mass));
  CHECK(proc.CommonForwardCutsForTest(mass_cut, MakeToyProductionLTSDphi(0.25 * gra::math::PI)));

  const nlohmann::json above_card = {{"FIDCUTS", {{"active", true}, {"FORWARD", {{"dPhi", {90.0, 180.0}}}}}}};
  const gra::FIDCUT    above      = ParsedFiducialCutsForTest(above_card);
  CHECK(proc.CommonForwardCutsForTest(above, MakeToyProductionLTSDphi(0.75 * gra::math::PI)));
  CHECK_FALSE(proc.CommonForwardCutsForTest(above, MakeToyProductionLTSDphi(0.25 * gra::math::PI)));

  const nlohmann::json full_card = {{"FIDCUTS", {{"active", true}, {"FORWARD", {{"dPhi", {0.0, 180.0}}}}}}};
  const gra::FIDCUT    full      = ParsedFiducialCutsForTest(full_card);
  CHECK(proc.CommonForwardCutsForTest(full, MakeToyProductionLTSDphi(0.0)));
  CHECK(proc.CommonForwardCutsForTest(full, MakeToyProductionLTSDphi(gra::math::PI)));

  const nlohmann::json wrapped_card = {{"FIDCUTS", {{"active", true}, {"FORWARD", {{"dPhi", {1.5, 2.5}}}}}}};
  const gra::FIDCUT    wrapped      = ParsedFiducialCutsForTest(wrapped_card);
  CHECK(proc.CommonForwardCutsForTest(wrapped,
                                      MakeToyProductionLTSDphi(gra::math::Deg2Rad(2.0), gra::math::Deg2Rad(179.0))));

  const nlohmann::json negative_card = {{"FIDCUTS", {{"active", true}, {"FORWARD", {{"dPhi", {-1.0, 90.0}}}}}}};
  CHECK_THROWS_AS(ParsedFiducialCutsForTest(negative_card), std::invalid_argument);

  const nlohmann::json overflow_card = {{"FIDCUTS", {{"active", true}, {"FORWARD", {{"dPhi", {90.0, 181.0}}}}}}};
  CHECK_THROWS_AS(ParsedFiducialCutsForTest(overflow_card), std::invalid_argument);
}

// Check photon virtuality and hadron transfer selections independently through the real parser
TEST_CASE("Independent forward transfers preserve the other beam phase space", "[MGraniitti][MProcess][FIDCUTS]") {
  const nlohmann::json card = {{"FIDCUTS", {{"active", true}, {"FORWARD", {{"t1", {0.0, 4.0}}}}}}};
  const gra::FIDCUT cuts = ParsedFiducialCutsForTest(card);
  CHECK(cuts.HasForwardCuts());
  gra::LORENTZSCALAR event;
  event.t1 = -4.0;
  event.t2 = -23.0;
  CHECK(cuts.PassForward(event));
  event.t1 = -4.1;
  CHECK_FALSE(cuts.PassForward(event));
  event.t1 = 0.0;
  CHECK(cuts.PassForward(event));
  auto reversed = card;
  reversed["FIDCUTS"]["FORWARD"] = {{"t2", {0.0, 4.0}}};
  const auto other = ParsedFiducialCutsForTest(reversed);
  std::swap(event.t1, event.t2);
  CHECK(other.PassForward(event));
  event.t2 = -4.1;
  CHECK_FALSE(other.PassForward(event));
  auto invalid = card;
  invalid["FIDCUTS"]["FORWARD"]["t1"] = {-1.0, 4.0};
  REQUIRE_THROWS_AS(ParsedFiducialCutsForTest(invalid), std::invalid_argument);
}

TEST_CASE("FORWARD Xi cuts use tagged losses independently of hard x", "[MGraniitti][MProcess][FIDCUTS]") {
  const nlohmann::json card = {{"FIDCUTS", {{"active", true}, {"FORWARD", {{"Xi", {0.02, 0.12}}}}}}};
  const gra::FIDCUT    cuts = ParsedFiducialCutsForTest(card);
  CHECK(cuts.forward_xi_active);
  CHECK(cuts.forward_xi_min == Approx(0.02));
  CHECK(cuts.forward_xi_max == Approx(0.12));
  CHECK(cuts.HasForwardCuts());

  ToyHelicityProcess proc;
  gra::LORENTZSCALAR event = MakeToyProductionLTSDphi(0.5 * gra::math::PI);
  event.x1                 = 0.004;
  event.x2                 = 0.006;
  event.xi1                = 0.02;
  event.xi2                = 0.12;
  event.has_xi1            = true;
  event.has_xi2            = true;
  CHECK(proc.CommonForwardCutsForTest(cuts, event));

  event.xi1 = 0.019;
  CHECK_FALSE(proc.CommonForwardCutsForTest(cuts, event));
  event.xi1 = 0.02;
  event.xi2 = 0.121;
  CHECK_FALSE(proc.CommonForwardCutsForTest(cuts, event));

  event.xi2     = 0.12;
  event.has_xi1 = false;
  event.has_xi2 = false;
  CHECK_FALSE(proc.CommonForwardCutsForTest(cuts, event));

  event.xi1     = 0.0;
  event.has_xi2 = true;
  CHECK(proc.CommonForwardCutsForTest(cuts, event));
}

TEST_CASE("MQuasiElastic applies inclusive forward fiducial ranges", "[MQuasiElastic][FIDCUTS]") {
  gra::LORENTZSCALAR event;
  event.t = 0.2;
  event.pfinal.resize(3);
  event.pfinal[1].SetPxPyPzM(0.1, 0.0, 0.0, 2.0);
  event.pfinal[2].SetPxPyPzM(-0.1, 0.0, 0.0, 2.5);

  ToyQuasiElasticProcess proc;
  gra::FIDCUT            transfer_cut;
  transfer_cut.active           = true;
  transfer_cut.forward_t_active = true;
  transfer_cut.forward_t_min    = 0.2;
  transfer_cut.forward_t_max    = 0.4;
  CHECK(proc.FiducialCutsForTest(transfer_cut, event, "EL"));
  transfer_cut.forward_t_min = 0.0;
  transfer_cut.forward_t_max = 0.2;
  CHECK(proc.FiducialCutsForTest(transfer_cut, event, "EL"));
  transfer_cut.forward_t_min = 0.21;
  transfer_cut.forward_t_max = 0.4;
  CHECK_FALSE(proc.FiducialCutsForTest(transfer_cut, event, "EL"));

  gra::FIDCUT mass_cut;
  mass_cut.active           = true;
  mass_cut.forward_M_active = true;
  mass_cut.forward_M_min    = 2.0;
  mass_cut.forward_M_max    = 2.5;
  CHECK(proc.FiducialCutsForTest(mass_cut, event, "SD"));
  CHECK(proc.FiducialCutsForTest(mass_cut, event, "DD"));
  mass_cut.forward_M_max = 2.4;
  CHECK_FALSE(proc.FiducialCutsForTest(mass_cut, event, "DD"));
  CHECK(proc.FiducialCutsForTest(mass_cut, event, "EL"));

  gra::FIDCUT xi_cut;
  xi_cut.active            = true;
  xi_cut.forward_xi_active = true;
  xi_cut.forward_xi_min    = 0.02;
  xi_cut.forward_xi_max    = 0.03;
  event.xi1                = 0.02;
  event.xi2                = 0.03;
  event.has_xi1            = true;
  event.has_xi2            = true;
  CHECK(proc.FiducialCutsForTest(xi_cut, event, "EL"));
  event.xi2 = 0.031;
  CHECK_FALSE(proc.FiducialCutsForTest(xi_cut, event, "EL"));
}

TEST_CASE("MGraniitti validates USERCUTS identifiers and event inputs", "[MGraniitti][USERCUTS]") {
  constexpr std::int64_t large_positive = 3000000000LL;
  constexpr std::int64_t large_negative = -3000000000LL;
  const nlohmann::json   positive_card  = {{"FIDCUTS", {{"active", true}, {"USERCUTS", large_positive}}}};
  CHECK_THROWS_AS(ParsedUserCutIDForTest(positive_card), std::invalid_argument);

  const nlohmann::json negative_card = {{"FIDCUTS", {{"active", true}, {"USERCUTS", large_negative}}}};
  CHECK_THROWS_AS(ParsedUserCutIDForTest(negative_card), std::invalid_argument);
  CHECK_FALSE(gra::IsKnownUserCut(large_positive));
  CHECK_FALSE(gra::IsKnownUserCut(large_negative));
  CHECK_FALSE(gra::UserCut(large_positive, gra::LORENTZSCALAR()));
  CHECK_FALSE(gra::UserCut(-3, gra::LORENTZSCALAR()));
  CHECK_FALSE(gra::UserCut(170804053, gra::LORENTZSCALAR()));

  gra::LORENTZSCALAR forward_protons;
  forward_protons.pfinal[1] = gra::M4Vec(0.1, 0.0, 1.0, 2.0);
  forward_protons.pfinal[2] = gra::M4Vec(-0.1, 0.0, -1.0, 2.0);
  CHECK(gra::UserCut(-3, forward_protons));
  forward_protons.pfinal[1] = gra::M4Vec(0.2, 0.0, 1.0, 2.0);
  forward_protons.pfinal[2] = gra::M4Vec(-0.2, 0.0, -1.0, 2.0);
  CHECK_FALSE(gra::UserCut(-3, forward_protons));

  gra::LORENTZSCALAR atlas_muons;
  atlas_muons.decaytree.resize(2);
  atlas_muons.decaytree[0].p4 = gra::M4Vec(7.0, 0.0, 0.0, 7.0);
  atlas_muons.decaytree[1].p4 = gra::M4Vec(-7.0, 0.0, 0.0, 7.0);
  atlas_muons.m2              = 400.0;
  CHECK(gra::UserCut(170804053, atlas_muons));
  atlas_muons.m2 = 100.0;
  CHECK_FALSE(gra::UserCut(170804053, atlas_muons));

  const nlohmann::json disabled_card = {{"FIDCUTS", {{"active", true}, {"USERCUTS", false}}}};
  CHECK(ParsedUserCutIDForTest(disabled_card) == 0);

  const nlohmann::json true_card = {{"FIDCUTS", {{"active", true}, {"USERCUTS", true}}}};
  CHECK_THROWS_AS(ParsedUserCutIDForTest(true_card), std::invalid_argument);

  const nlohmann::json fractional_card = {{"FIDCUTS", {{"active", true}, {"USERCUTS", 1.5}}}};
  CHECK_THROWS_AS(ParsedUserCutIDForTest(fractional_card), std::invalid_argument);

  const std::uint64_t  overflow      = static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max()) + 1U;
  const nlohmann::json overflow_card = {{"FIDCUTS", {{"active", true}, {"USERCUTS", overflow}}}};
  CHECK_THROWS_AS(ParsedUserCutIDForTest(overflow_card), std::invalid_argument);
}

// Check the slope normalization through the process amplitude API
TEST_CASE("Flat matrix elements apply B at cross-section level", "[gra::MProcess][slope]") {
  ProcessProbe process;
  process.SetModelTune(gra::MModelTune::Load(modelfile));
  process.SetFLATAMP(1);

  gra::LORENTZSCALAR lts;
  lts.s                 = 25.0;
  lts.t1                = -0.2;
  lts.t2                = -0.3;
  const double expected = gra::math::pow2(lts.s) * std::exp(process.GetModelTune()->Flat().B * (lts.t1 + lts.t2));

  CHECK(process.GetFlatAmp2(lts) == Approx(expected).margin(1e-15));
  CHECK_NOTHROW(process.SetFLATAMP(0));
  CHECK_NOTHROW(process.SetFLATAMP(4));
  CHECK_THROWS_AS(process.SetFLATAMP(-1), std::invalid_argument);
  CHECK_THROWS_AS(process.SetFLATAMP(5), std::invalid_argument);
}

// Check separate process tunes keep independent amplitude slopes
TEST_CASE("Flat process weights retain their immutable tune", "[gra::MProcess][slope][threading]") {
  const std::filesystem::path dir = std::filesystem::path("tmp") / "test_flat_tunes";
  std::filesystem::create_directories(dir);
  for (const auto &filename : {"NUMERICS.json", "CON_MP.json", "CON_XP.json", "CON_GP.json", "CON_TP.json"}) {
    std::filesystem::copy_file(gra::ResolveModelDataFile("TUNE0", filename), dir / filename,
                               std::filesystem::copy_options::overwrite_existing);
  }
  // Reload the same file to exercise independence from later tune changes
  const auto load = [&](double B) {
    auto card = nlohmann::json::parse(gra::aux::GetInputData(modelfile));
    card["PARAM_FLAT"]["B"] = B;
    const auto path = dir / "GENERAL.json";
    std::ofstream output(path);
    output << card;
    output.close();
    return gra::MModelTune::Load(path.string());
  };
  ProcessProbe first;
  first.SetModelTune(load(2.0));
  first.SetFLATAMP(1);
  ProcessProbe second;
  second.SetModelTune(load(7.0));
  second.SetFLATAMP(1);
  gra::LORENTZSCALAR lts;
  lts.s = 25.0;
  lts.s_hat = 4.0;
  lts.t1 = -0.2;
  lts.t2 = -0.3;
  for (int mode : {1, 2, 3}) {
    first.SetFLATAMP(mode);
    second.SetFLATAMP(mode);
    const double divisor = mode == 1 ? 1.0 : mode == 2 ? 2.0 : 4.0;
    auto a = std::async(std::launch::async, [&] { return first.GetFlatAmp2(lts); });
    auto b = std::async(std::launch::async, [&] { return second.GetFlatAmp2(lts); });
    CHECK(a.get() == Approx(625.0 * std::exp(-1.0) / divisor));
    CHECK(b.get() == Approx(625.0 * std::exp(-3.5) / divisor));
  }
}

// Check named decay syntax preserves baryon number and charge conjugation
TEST_CASE("PDG names and multiplets preserve particle-antiparticle identity", "[gra::MPDG]") {
  const auto &pdg = LoadedPDGTable();
  CHECK(pdg.FindByPDGName("Sigma+").pdg == 3222);
  CHECK(pdg.FindByPDGName("Sigma-").pdg == 3112);
  CHECK(pdg.FindByPDGName("Sigma-~").pdg == -3222);
  CHECK(pdg.FindByPDGName("Sigma+~").pdg == -3112);
  CHECK(pdg.FindByPDGName("p-").pdg == -2212);
  std::vector<gra::MDecayBranch> tree;
  pdg.TokenizeProcess("Sigma+ Sigma-~", 0, tree);
  REQUIRE(tree.size() == 2);
  CHECK(tree[0].p.pdg == -tree[1].p.pdg);
  CHECK(tree[0].p.chargeX3 == -tree[1].p.chargeX3);
  for (const auto &[id, particle] : pdg.PDG_table) {
    if (id <= 0 || particle.spinX2 < 0 || particle.spinX2 % 2 == 0) { continue; }
    CAPTURE(id);
    const auto &anti = pdg.FindByPDG(-id);
    CHECK(anti.P == -particle.P);
    CHECK(anti.chargeX3 == -particle.chargeX3);
    CHECK(anti.mass == Approx(particle.mass));
    CHECK(anti.width == Approx(particle.width));
    CHECK(pdg.FindByPDGName(particle.name).pdg == id);
    CHECK(pdg.FindByPDGName(anti.name).pdg == -id);
  }
}

// Check PDG spin and parity through allowed and forbidden physical decay waves
TEST_CASE("PDG quantum numbers give the physical meson and baryon decay waves", "[gra::MPDG][gra::spin]") {
  const auto &pdg = LoadedPDGTable();
  // [REFERENCE: PDG 2025, Mesons and Baryons Summary Tables]
  for (const auto &decay : std::vector<std::array<int, 4>>{
           {117, 211, -211, 3}, {9010225, 211, -211, 2}, {9030221, 211, -211, 0},
           {3124, 2212, -321, 2}, {22112, 2212, -211, 0}, {2222, 2212, 211, 0},
           {13224, 3122, 211, 2}, {104122, 4112, 211, 2},
           {23114, 3122, -211, 2}, {23214, 3122, 111, 2}, {23224, 3122, 211, 2},
           {104314, 4312, 111, 0}, {104324, 4322, 111, 0},
           {104312, 4312, 111, 2}, {104322, 4322, 111, 2}}) {
    CAPTURE(decay);
    auto mother = pdg.FindByPDG(decay[0]);
    const auto &a = pdg.FindByPDG(decay[1]);
    const auto &b = pdg.FindByPDG(decay[2]);
    gra::HELMatrix hel;
    hel.BR = 1.0;
    hel.P_symmetry = true;
    hel.C_symmetry = mother.spinX2 % 2 == 0;
    hel.alpha_ls.Set(decay[3], a.spinX2 + b.spinX2, 1.0);
    auto exchanged = hel;
    auto forbidden = hel;
    REQUIRE_NOTHROW(gra::spin::InitTMatrix(hel, mother, a, b, false, "PDG physical decay", false, false));
    REQUIRE_NOTHROW(gra::spin::InitTMatrix(exchanged, mother, b, a, false, "PDG exchanged daughters", false, false));
    CHECK(hel.T.FrobNorm2() == Approx(exchanged.T.FrobNorm2()));
    mother.P = -mother.P;
    REQUIRE_THROWS_AS(gra::spin::InitTMatrix(forbidden, mother, a, b, false, "PDG forbidden parity", false, false),
                      std::invalid_argument);
  }
  for (int id : {130, 310}) {
    CAPTURE(id);
    const auto &kaon = pdg.FindByPDG(id);
    CHECK(kaon.spinX2 == 0);
    CHECK(kaon.P == -1);
    CHECK_FALSE(pdg.PDG_table.contains(-id));
    gra::HELMatrix hel;
    hel.BR = 1.0;
    hel.P_symmetry = false;
    hel.C_symmetry = false;
    hel.alpha_ls.Set(0, 0, 1.0);
    REQUIRE_NOTHROW(gra::spin::InitTMatrix(hel, kaon, pdg.FindByPDG(211), pdg.FindByPDG(-211), false,
                                         "neutral kaon weak decay", false, false));
  }
}

// Check colored intermediates remain visible through recursively parsed decays
TEST_CASE("PDG color representations reach cascade validation", "[gra::MPDG][color]") {
  const auto &pdg = LoadedPDGTable();
  for (int id : {1, 2, 3, 4, 5, 6}) {
    CHECK(pdg.FindByPDG(id).color == 3);
    CHECK(pdg.FindByPDG(-id).color == -3);
  }
  CHECK(pdg.FindByPDG(21).color == 8);
  std::vector<gra::MDecayBranch> top, gluon, colorless;
  pdg.TokenizeProcess("t > {b W+}", 0, top);
  pdg.TokenizeProcess("Z > {g > {u u~} gamma}", 0, gluon);
  pdg.TokenizeProcess("Z > {b b~}", 0, colorless);
  CHECK(gra::DecayTreeHasColoredIntermediate(top));
  CHECK(gra::DecayTreeHasColoredIntermediate(gluon));
  CHECK_FALSE(gra::DecayTreeHasColoredIntermediate(colorless));
}

// Reject malformed fixed-column numbers without changing the PDG reference file
TEST_CASE("PDG fixed-column input rejects failed and partial numeric conversions", "[gra::MPDG][validation]") {
  ModelParamRestoreGuard restore;
  gra::MODELPARAM = "TUNE0";
  const std::filesystem::path dir = std::filesystem::path("tmp") / "test_pdg_numbers";
  std::filesystem::create_directories(dir);
  std::ifstream input(pdgfile);
  std::string proton;
  std::string line;
  while (std::getline(input, line)) {
    if (line.size() < 108 || line.front() == '*') { continue; }
    std::istringstream ids(line.substr(0, 32));
    int id = 0;
    if ((ids >> id) && id == 2212) { proton = line; break; }
  }
  REQUIRE_FALSE(proton.empty());
  const auto read = [&](std::string row) {
    const auto path = dir / "particle.mcd";
    std::ofstream output(path);
    output << row << '\n';
    output.close();
    gra::MPDG pdg;
    pdg.ReadParticleData(path.string());
    return pdg;
  };
  REQUIRE_NOTHROW(read(proton));
  for (std::size_t offset : {33U, 70U}) {
    for (const std::string number : {"BAD", "1.0junk", "NaN", "Inf", "1E999", "-1.0"}) {
      CAPTURE(offset, number);
      auto invalid = proton;
      invalid.replace(offset, 18, number + std::string(18 - number.size(), ' '));
      REQUIRE_THROWS_AS(read(invalid), std::invalid_argument);
    }
  }
  auto missing_mass = proton;
  missing_mass.replace(33, 18, 18, ' ');
  REQUIRE_THROWS_AS(read(missing_mass), std::invalid_argument);
  auto missing_width = proton;
  missing_width.replace(70, 36, 36, ' ');
  REQUIRE_NOTHROW(read(missing_width));
  auto invalid_id = proton;
  invalid_id.replace(0, 8, "  2212x ");
  REQUIRE_THROWS_AS(read(invalid_id), std::invalid_argument);
  REQUIRE_THROWS_AS(read(proton + ","), std::invalid_argument);
  REQUIRE_THROWS_AS(read(proton + ",  "), std::invalid_argument);

  // A failed reload must preserve both particles and their charge conjugates
  auto retained = read(proton);
  const auto size = retained.PDG_table.size();
  const auto path = dir / "failed_reload.mcd";
  {
    std::ofstream output(path);
    output << proton << '\n' << invalid_id << '\n';
  }
  REQUIRE_THROWS_AS(retained.ReadParticleData(path.string()), std::invalid_argument);
  CHECK(retained.PDG_table.size() == size);
  CHECK(retained.FindByPDGName("p+").pdg == 2212);
  CHECK(retained.FindByPDGName("p-").pdg == -2212);
}

// Check structured extra particles and fixed-column blank-line handling
TEST_CASE("MPDG reads structured extra particle table", "[gra::MPDG]") {
  const gra::MPDG            &pdg        = LoadedPDGTable();
  const std::filesystem::path extra_path = std::filesystem::path(pdgfile).parent_path() / "TUNE0" / "PDG_EXTRA.json";
  const auto                  extra      = nlohmann::json::parse(gra::aux::GetInputData(extra_path.string()));

  REQUIRE(extra.at("PARAM_PDG").is_object());
  REQUIRE_FALSE(extra.at("PARAM_PDG").empty());

  for (const auto &[label, entry] : extra.at("PARAM_PDG").items()) {
    CAPTURE(label);
    CHECK_FALSE(label.empty());
    CHECK_FALSE(std::all_of(label.begin(), label.end(), [](unsigned char c) { return std::isdigit(c); }));
    const int             pdg_id   = entry.at("PDG").get<int>();
    const gra::MParticle &particle = pdg.FindByPDG(pdg_id);

    CAPTURE(pdg_id);
    CHECK(particle.name == entry.at("name").get<std::string>());
    CHECK(particle.mass == Approx(entry.at("mass").get<double>()));
    CHECK(particle.width == Approx(entry.at("width").get<double>()));
    CHECK(particle.chargeX3 == entry.at("chargeX3").get<int>());
    CHECK(particle.spinX2 == entry.at("spinX2").get<int>());
    CHECK(particle.isospinX2 == entry.at("isospinX2").get<int>());
    CHECK(particle.P == entry.at("P").get<int>());
    CHECK(particle.C == entry.at("C").get<int>());
    CHECK(particle.G == entry.at("G").get<int>());
    CHECK(particle.L == static_cast<unsigned int>(entry.at("L").get<int>()));
    CHECK(particle.color == entry.at("color").get<int>());
  }

  CHECK(pdg.FindByPDG(211).C == 0);
  CHECK(pdg.FindByPDG(-211).C == 0);
  CHECK(pdg.FindByPDG(321).C == 0);
  CHECK(pdg.FindByPDG(-321).C == 0);
  CHECK(pdg.FindByPDG(311).C == 0);
  CHECK(pdg.FindByPDG(-311).C == 0);
  CHECK(pdg.FindByPDG(311).P == pdg.FindByPDG(-311).P);
  CHECK(pdg.FindByPDG(311).mass == Approx(pdg.FindByPDG(-311).mass));
  CHECK(pdg.FindByPDG(421).C == 0);
  CHECK(pdg.FindByPDG(-421).C == 0);
  for (const int neutrino : {12, 14, 16}) {
    CHECK(pdg.FindByPDG(neutrino).mass > 0.0);
    CHECK(pdg.FindByPDG(-neutrino).mass > 0.0);
  }
  CHECK(pdg.PDG_table.count(100223) == 1);
  CHECK(pdg.PDG_table.count(1000223) == 0);
  CHECK(pdg.PDG_table.count(9000113) == 0);
  CHECK(pdg.PDG_table.count(9000213) == 0);
  CHECK(pdg.FindByPDG(991).spinX2 == 0);
  CHECK(pdg.FindByPDG(993).spinX2 == 2);
  CHECK(pdg.FindByPDG(995).spinX2 == 4);
  CHECK(pdg.FindByPDG(9000221).P == 1);
  CHECK(pdg.FindByPDG(9000221).C == 1);
  CHECK(pdg.FindByPDG(9991).C == -1);
  CHECK(pdg.FindByPDG(9915).C == 1);
  CHECK(pdg.FindByPDG(9933).C == -1);
  CHECK(pdg.FindByPDG(990).P == 1);
  CHECK(pdg.FindByPDG(990).C == 1);
  CHECK(pdg.FindByPDG(9990).P == -1);
  CHECK(pdg.FindByPDG(9990).C == -1);
  CHECK(pdg.FindByPDG(9910).P == 1);
  CHECK(pdg.FindByPDG(9910).C == 1);
  CHECK(pdg.FindByPDG(9930).P == -1);
  CHECK(pdg.FindByPDG(9930).C == -1);
  CHECK(pdg.FindByPDG(990).isospinX2 == 0);
  CHECK(pdg.FindByPDG(990).G == 1);
  CHECK(pdg.FindByPDG(9990).isospinX2 == 0);
  CHECK(pdg.FindByPDG(9990).G == -1);
  CHECK(pdg.FindByPDG(9910).isospinX2 == 0);
  CHECK(pdg.FindByPDG(9910).G == 1);
  CHECK(pdg.FindByPDG(9920).isospinX2 == 2);
  CHECK(pdg.FindByPDG(9920).G == -1);
  CHECK(pdg.FindByPDG(9930).isospinX2 == 2);
  CHECK(pdg.FindByPDG(9930).G == 1);
  CHECK(pdg.FindByPDG(9940).isospinX2 == 0);
  CHECK(pdg.FindByPDG(9940).G == -1);

  const gra::MParticle &monopolium = pdg.FindByPDG(gra::PDG::PDG_monopolium);
  CHECK(monopolium.name == "monopolium(0)");
  CHECK(monopolium.spinX2 == 0);
  CHECK(monopolium.P == 1);
  CHECK(monopolium.C == 1);
  CHECK(pdg.PDG_table.count(-gra::PDG::PDG_monopolium) == 0);

  const gra::MParticle &monopole     = pdg.FindByPDG(gra::PDG::PDG_monopole);
  const gra::MParticle &antimonopole = pdg.FindByPDG(-gra::PDG::PDG_monopole);
  CHECK(monopole.name == "dirac-monopole");
  CHECK(monopole.spinX2 == 1);
  CHECK(antimonopole.spinX2 == 1);
  CHECK(monopole.pdg == -antimonopole.pdg);

  const std::filesystem::path blank_line_table = std::filesystem::path("tmp") / "testbench_pdg_blank_line.mcd";
  std::filesystem::create_directories(blank_line_table.parent_path());
  {
    std::ifstream input(pdgfile);
    std::ofstream output(blank_line_table);
    REQUIRE(input.is_open());
    REQUIRE(output.is_open());
    output << '\n' << input.rdbuf();
  }
  gra::MPDG blank_line_pdg;
  CHECK_NOTHROW(blank_line_pdg.ReadParticleData(blank_line_table.string()));
  CHECK(blank_line_pdg.FindByPDG(211).mass == Approx(pdg.FindByPDG(211).mass));

  gra::PROC_002_QED_YY_MONOPOLIUM0 process;
  gra::LORENTZSCALAR               lts;
  CHECK(process.RootResonancePDG(lts) == gra::PDG::PDG_monopolium);
}

// Check strict isospin and G-parity parsing for extra particles
//
TEST_CASE("MPDG rejects invalid extra-particle isospin and G-parity", "[gra::MPDG][validation]") {
  const std::string old_modelparam = gra::MODELPARAM;
  struct RestoreModelParam {
    std::string value;
    // Restore the tune selection after each test section
    ~RestoreModelParam() { gra::MODELPARAM = value; }
  } restore{old_modelparam};

  const std::filesystem::path tune_dir = std::filesystem::path("tmp") / "testbench_pdg_isospin_g";
  std::filesystem::create_directories(tune_dir);
  const auto valid =
      nlohmann::json::parse(gra::aux::GetInputData(gra::ResolveModelDataFile("TUNE0", "PDG_EXTRA.json")));

  const auto read_extra = [&](const nlohmann::json &extra) {
    std::ofstream out(tune_dir / "PDG_EXTRA.json");
    out << extra.dump();
    out.close();
    gra::MODELPARAM = tune_dir.string();
    gra::MPDG particle_data;
    particle_data.ReadParticleData(pdgfile);
  };

  SECTION("missing G-parity") {
    auto invalid = valid;
    invalid.at("PARAM_PDG").at("pomeron_analytic").erase("G");
    REQUIRE_THROWS_AS(read_extra(invalid), std::invalid_argument);
  }

  SECTION("half-integer isospin with defined G-parity") {
    auto invalid                                                   = valid;
    invalid.at("PARAM_PDG").at("pomeron_analytic").at("isospinX2") = 1;
    REQUIRE_THROWS_AS(read_extra(invalid), std::invalid_argument);
  }

  SECTION("neutral-state G-parity sign") {
    auto invalid                                              = valid;
    invalid.at("PARAM_PDG").at("f2_reggeon_analytic").at("G") = -1;
    REQUIRE_THROWS_AS(read_extra(invalid), std::invalid_argument);
  }
  SECTION("charged particles have undefined C-parity") {
    auto charged = valid;
    auto &entry = charged["PARAM_PDG"]["charged_scalar"];
    entry = {{"PDG", 8000001}, {"name", "charged_scalar"}, {"mass", 1.0}, {"width", 0.1},
             {"chargeX3", 3}, {"spinX2", 0}, {"P", 1}, {"C", 0}, {"isospinX2", 2}, {"G", 0},
             {"L", 0}, {"color", 0}};
    REQUIRE_NOTHROW(read_extra(charged));
    for (const int parity : {-1, 1}) {
      entry["C"] = parity;
      REQUIRE_THROWS_AS(read_extra(charged), std::invalid_argument);
    }
  }
  SECTION("oversized integer labels cannot narrow to valid quantum numbers") {
    for (const std::string field : {"PDG", "spinX2", "chargeX3", "P", "C", "isospinX2", "G", "L", "color"}) {
      CAPTURE(field);
      auto invalid = valid;
      auto &entry = invalid["PARAM_PDG"]["pomeron_analytic"];
      entry[field] = 4294967296ULL + entry[field].get<int>();
      REQUIRE_THROWS_AS(read_extra(invalid), std::invalid_argument);
    }
  }
}

// Test initializing Tensor Pomeron amplitudes
//
TEST_CASE("MProcess keeps weak-width stable leaves on shell", "[MProcess]") {
  ProcessProbe process;
  gra::MDecayBranch   neutron;
  neutron.p.mass  = 0.93956542052;
  neutron.p.width = 7.49e-28;

  const double mass = process.SampleBranchMass(neutron);

  CHECK(mass == Approx(neutron.p.mass));
  CHECK(neutron.mass_proposal_norm == Approx(1.0));
  CHECK_FALSE((neutron.mass_proposal == gra::MassProposal::BreitWigner));
  CHECK_FALSE((neutron.mass_proposal == gra::MassProposal::Uniform));
  CHECK((neutron.mass_proposal == gra::MassProposal::Fixed));
  CHECK(neutron.mass_proposal_min2 == Approx(gra::math::pow2(neutron.p.mass)));
  CHECK(neutron.mass_proposal_max2 == Approx(gra::math::pow2(neutron.p.mass)));
}

TEST_CASE("MProcess samples resolvable-width stable leaves off shell", "[MProcess]") {
  ProcessProbe process;
  gra::MDecayBranch   particle;
  particle.p.mass  = 1.0;
  particle.p.width = 0.1;

  const double mass = process.SampleBranchMass(particle);

  CHECK(mass >= 0.5);
  CHECK(mass <= 1.5);
  CHECK(particle.mass_proposal_norm > 1.0);
  CHECK((particle.mass_proposal == gra::MassProposal::BreitWigner));
  CHECK_FALSE((particle.mass_proposal == gra::MassProposal::Uniform));
  CHECK_FALSE((particle.mass_proposal == gra::MassProposal::Fixed));
  CHECK(particle.mass_proposal_min2 == Approx(0.25));
  CHECK(particle.mass_proposal_max2 == Approx(2.25));
  CHECK(pow2(mass) >= particle.mass_proposal_min2);
  CHECK(pow2(mass) <= particle.mass_proposal_max2);
}

// Keep a below-threshold pole's allowed Breit-Wigner tail in the generated decay mass
TEST_CASE("MProcess samples a nonempty resonance tail above the daughter threshold",
          "[MProcess][proposal][tails]") {
  ProcessProbe process;
  process.SetWIDTHMIN(0.0);
  process.SetOFFSHELL(3e20);
  gra::MDecayBranch branch;
  branch.p.mass = 1.0;
  branch.p.width = 1e-20;
  branch.legs.resize(2);
  for (auto &leg : branch.legs) { leg.p.mass = 1.0; }
  for (unsigned int i = 0; i < 256; ++i) {
    const double mass = process.SampleBranchMass(branch);
    REQUIRE(mass >= 2.0);
    REQUIRE(mass <= 4.0);
    REQUIRE((branch.mass_proposal == gra::MassProposal::BreitWigner));
    const double lower = branch.mass_proposal_min2;
    const double upper = branch.mass_proposal_max2;
    const double expected = (upper - lower) / ((lower - 1.0) * (upper - 1.0));
    CHECK(branch.mass_proposal_norm == Approx(expected).epsilon(1e-12));
  }
}

TEST_CASE("MProcess classifies fixed and continuous mass proposals with WIDTH_MIN", "[MProcess][proposal]") {
  ProcessProbe process;
  process.SetWIDTHMIN(1e-3);

  gra::MDecayBranch particle;
  particle.p.mass         = 1.0;
  particle.p.width        = 5e-4;
  const double fixed_mass = process.SampleBranchMass(particle);
  REQUIRE(fixed_mass == Approx(particle.p.mass));
  REQUIRE((particle.mass_proposal == gra::MassProposal::Fixed));
  REQUIRE_FALSE((particle.mass_proposal == gra::MassProposal::BreitWigner));
  REQUIRE_FALSE((particle.mass_proposal == gra::MassProposal::Uniform));

  particle.p.width          = 2e-3;
  const double sampled_mass = process.SampleBranchMass(particle);
  REQUIRE(sampled_mass >= 0.99);
  REQUIRE(sampled_mass <= 1.01);
  REQUIRE_FALSE((particle.mass_proposal == gra::MassProposal::Fixed));
  REQUIRE((particle.mass_proposal == gra::MassProposal::BreitWigner));
  REQUIRE_FALSE((particle.mass_proposal == gra::MassProposal::Uniform));

  process.SetWIDTHMIN(0.0);
  particle.p.width             = 0.0;
  const double zero_width_mass = process.SampleBranchMass(particle);
  REQUIRE(zero_width_mass == Approx(particle.p.mass));
  REQUIRE((particle.mass_proposal == gra::MassProposal::Fixed));
  REQUIRE_FALSE((particle.mass_proposal == gra::MassProposal::BreitWigner));
  REQUIRE_FALSE((particle.mass_proposal == gra::MassProposal::Uniform));
}

TEST_CASE("MProcess fixes floating-point collapsed mass windows at the pole", "[MProcess][proposal]") {
  ProcessProbe process;
  process.SetWIDTHMIN(0.0);

  gra::MDecayBranch particle;
  particle.p.mass  = 1.0;
  particle.p.width = 1.0e-30;

  for (const bool flat_mass_squared : {false, true}) {
    process.SetFLATMASS2(flat_mass_squared);
    const double mass = process.SampleBranchMass(particle);
    REQUIRE(mass == Approx(particle.p.mass));
    REQUIRE(particle.mass_proposal_norm == Approx(1.0));
    REQUIRE(particle.mass_proposal_min2 == Approx(1.0));
    REQUIRE(particle.mass_proposal_max2 == Approx(1.0));
    REQUIRE((particle.mass_proposal == gra::MassProposal::Fixed));
    REQUIRE_FALSE((particle.mass_proposal == gra::MassProposal::BreitWigner));
    REQUIRE_FALSE((particle.mass_proposal == gra::MassProposal::Uniform));
  }
}

TEST_CASE("MProcess zero OFFSHELL fixes resolvable widths at the pole mass", "[MProcess][proposal]") {
  ProcessProbe process;
  process.SetOFFSHELL(0.0);

  gra::MDecayBranch particle;
  particle.p.mass   = 1.0;
  particle.p.width  = 0.1;
  const double mass = process.SampleBranchMass(particle);

  REQUIRE(mass == Approx(particle.p.mass));
  REQUIRE((particle.mass_proposal == gra::MassProposal::Fixed));
  REQUIRE_FALSE((particle.mass_proposal == gra::MassProposal::BreitWigner));
  REQUIRE_FALSE((particle.mass_proposal == gra::MassProposal::Uniform));
  REQUIRE(particle.mass_proposal_min2 == Approx(1.0));
  REQUIRE(particle.mass_proposal_max2 == Approx(1.0));
  REQUIRE_THROWS_AS(process.SetOFFSHELL(-1.0), std::invalid_argument);
}

// Verify independent proper lifetimes, inherited positions and event length units
TEST_CASE("HepMC cascades propagate every unstable descendant from its production vertex",
          "[gra::MProcess][HepMC][lifetime]") {
  for (const auto units : {HepMC3::Units::MM, HepMC3::Units::CM}) {
    for (const bool prompt_root : {false, true}) {
      for (const bool prompt_child : {false, true}) {
        CAPTURE(units, prompt_root, prompt_child);
        auto branch = MakeFlightTree();
        if (prompt_child) { branch.legs[0].p.tau = 0.0; }
        HepMC3::GenEvent evt(HepMC3::Units::GEV, units);
        auto production = std::make_shared<HepMC3::GenVertex>(HepMC3::FourVector(1.0, 2.0, 3.0, 4.0));
        auto mother = std::make_shared<HepMC3::GenParticle>(
            gra::aux::M4Vec2HepMC3(branch.p4), branch.p.pdg, gra::PDG::PDG_DECAY);
        production->add_particle_out(mother);
        evt.add_vertex(production);
        gra::MRandom random;
        random.SetSeed(98513);
        gra::MRandom reference = random;
        gra::record::WriteBranch(branch, mother, evt, random, !prompt_root);

        const double scale = units == HepMC3::Units::MM ? 1e3 : 1e2;
        gra::M4Vec expected(1.0, 2.0, 3.0, 4.0);
        const std::array<const gra::MDecayBranch *, 3> chain = {
            &branch, &branch.legs[0], &branch.legs[0].legs[0]};
        auto particle = mother;
        for (const auto &i : indices(chain)) {
          const auto &decay = *chain[i];
          const bool displaced = !(i == 0 && prompt_root) && decay.p.tau > 0.0;
          if (displaced) {
            const double tau = reference.ExpRandom(1.0 / decay.p.tau);
            // x^mu = x_production^mu + p^mu c tau / m
            expected += decay.p4 * (gra::PDG::c * scale * tau / decay.p4.M());
          }
          REQUIRE(particle->end_vertex() != nullptr);
          const auto position = particle->end_vertex()->position();
          CHECK(position.x() == Approx(expected.X()).margin(1e-10));
          CHECK(position.y() == Approx(expected.Y()).margin(1e-10));
          CHECK(position.z() == Approx(expected.Z()).margin(1e-10));
          CHECK(position.t() == Approx(expected.E()).margin(1e-10));
          if (displaced) {
            CHECK(position.t() > particle->production_vertex()->position().t());
          } else {
            CHECK(position.t() == Approx(particle->production_vertex()->position().t()));
          }
          gra::M4Vec daughters;
          for (const auto &daughter : particle->end_vertex()->particles_out()) {
            daughters += gra::aux::HepMC2M4Vec(daughter->momentum());
          }
          CHECK(gra::math::CheckEMC(decay.p4 - daughters));
          particle = particle->end_vertex()->particles_out().front();
        }
      }
    }
  }
}

TEST_CASE("Common HepMC record writes the resolved root resonance identity", "[gra::MProcess][HepMC][root-resonance]") {
  gra::LORENTZSCALAR lts         = MakeToyPhotoZFFbar(13);
  lts.process.root_resonance_pdg = gra::PDG::PDG_monopolium;

  ToyPhotoZRecordProcess recorder;
  recorder.SetModelTune(gra::MModelTune::Load(modelfile));
  recorder.state.lts = lts;
  HepMC3::GenEvent evt(HepMC3::Units::GEV, HepMC3::Units::MM);
  REQUIRE(recorder.Record(evt));

  int monopolium_count = 0;
  int system_count     = 0;
  for (const auto &particle : evt.particles()) {
    if (particle->pid() == gra::PDG::PDG_monopolium) { ++monopolium_count; }
    if (particle->pid() == gra::PDG::PDG_system) { ++system_count; }
  }
  REQUIRE(monopolium_count == 1);
  REQUIRE(system_count == 0);
}

TEST_CASE("Direct integration reports reduced chi2 across batches", "[integration][chi2]") {
  gra::Stats             statistics;
  gra::MEventWeightState accepted;
  constexpr std::size_t  batch_events = 10;
  constexpr double       half_width   = 0.3;
  for (const double mean : {0.9, 1.1, 1.0}) {
    statistics.ResetIntegrationBatch();
    for (std::size_t i = 0; i < batch_events; ++i) {
      const double weight = mean + (i < batch_events / 2 ? -half_width : half_width);
      statistics.ObserveLogSample(accepted, weight, 0.0, gra::SamplingStage::Integration);
    }
    statistics.UpdateIntegrationChi2();
  }
  CHECK(statistics.chi2 == Approx(1.0).margin(1e-12));
}
