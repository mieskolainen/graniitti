// Regge model tests
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <limits>
#include <numeric>

// Own
#include "Graniitti/Kinematics/MEventKinematics.h"
#include "Graniitti/Kinematics/MFactorized.h"
#include "Graniitti/Math/MCombinatorics.h"
#include "Graniitti/Regge/MReggeMPXP.h"
#include "Graniitti/Regge/MReggeMP.h"
#include "Graniitti/Regge/MReggeMulti.h"
#include "Graniitti/Regge/MReggeNumerics.h"
#include "Graniitti/Regge/MReggeSig.h"
#include "Graniitti/Regge/MReggeSource.h"
#include "Graniitti/Regge/MReggeInit.h"
#include "Graniitti/Regge/MReggeMPInit.h"
#include "Graniitti/Regge/MReggeXPInit.h"
#include "Graniitti/Regge/MReggeGPInit.h"
#include "Graniitti/Regge/MReggeXP.h"
#include "Graniitti/Tech/MException.h"
#include "support/models_test_support.hh"

namespace {

// Select explicit fusion inputs for tests of MP pole and derivative operators
void SetMPFusion(gra::PARAM_RES &res) {
  res.spin_basis = "none";
  for (auto &channel : res.MP.channels) {
    channel.basis = gra::ReggeVertexBasis::AutoMinL;
    channel.Lambda = 1.0;
  }
}

// Build an exact central pair above the proton threshold with unequal beam transfers
gra::LORENTZSCALAR AsymmetricContinuumPairForTest(int first, int second) {
  auto lts = DirectCentralPairLTSForTest(first, second);
  lts.pfinal[1].SetPxPyPzM(0.16 * std::cos(0.21), 0.16 * std::sin(0.21), 97.7, lts.beam1.mass);
  lts.pfinal[2].SetPxPyPzM(0.60 * std::cos(2.61), 0.60 * std::sin(2.61), -98.1, lts.beam2.mass);
  lts.pfinal[0] = lts.pbeam1 + lts.pbeam2 - lts.pfinal[1] - lts.pfinal[2];
  const auto pair = TwoBodyRestKinematics(lts.pfinal[0].M(), lts.decaytree[0].p.mass,
                                         lts.decaytree[1].p.mass, 0.91, -0.38);
  for (const auto i : indices(lts.decaytree)) {
    lts.decaytree[i].p4 = BoostFromRestFrame(pair[i], lts.pfinal[0]);
  }
  RefreshToyDerivedKinematicsPreserveDecay(lts);
  return lts;
}

// Evaluate fixed-spin t/u amplitudes through freshly built poles and full exchange tensors
std::pair<MMatrix<std::complex<double>>, MMatrix<std::complex<double>>> FullPoleContinuum(
    const gra::LORENTZSCALAR &lts, gra::ReggeProductionModel model, double s0, std::size_t channel) {
  const auto                     &tree  = lts.process.CONT_PRODUCTIONTREE[channel];
  const auto                     &pole  = lts.process.CONTINUUM_POLE[channel];
  const std::array<gra::M4Vec, 2> final = {
      gra::kinematics::BoostToRestFrame(lts.decaytree[0].p4, lts.pfinal[0], "full continuum upper"),
      gra::kinematics::BoostToRestFrame(lts.decaytree[1].p4, lts.pfinal[0], "full continuum lower")};
  auto lower_axis = lts.q2_in_X;
  lower_axis.Flip3();
  const gra::spin::ForwardSpec forward{lts.process.FORWARD_VERTEX, s0, gra::ExchangeBasisType::ReducedRegge};
  const auto                   upper =
      gra::spin::Forward(lts, tree[0], lts.pbeam1, lts.pfinal[1], false,
                         gra::spin::Rows(tree[0], lts.process.FORWARD_NOFLIP), lts.process.PHOTON_VERTEX, forward);
  const auto lower =
      gra::spin::Forward(lts, tree[1], lts.pbeam2, lts.pfinal[2], true,
                         gra::spin::Rows(tree[1], lts.process.FORWARD_NOFLIP), lts.process.PHOTON_VERTEX, forward);
  // Build both orderings with the same physical internal numerator used by the model
  const auto central = [&](std::size_t ordering) {
    const auto &up     = pole[2 * ordering].Pole();
    const auto &dn     = pole[2 * ordering + 1].Pole();
    const auto  first  = up.helicity.exchange_basis == gra::ExchangeBasisType::ReducedRegge
                             ? gra::spin::EvaluateCrossedPole(up)
                             : gra::spin::EvaluatePoleSubvertex(up, final[ordering], lts.q1_in_X, false);
    const auto  second = dn.helicity.exchange_basis == gra::ExchangeBasisType::ReducedRegge
                             ? gra::spin::EvaluateCrossedPole(dn)
                             : gra::spin::EvaluatePoleSubvertex(dn, final[1 - ordering], lower_axis, true);
    if (model == gra::ReggeProductionModel::XP) {
      const auto numerator = gra::xpom::ReducedNumerator(lts.decaytree[ordering].p, lts.q1_in_X - final[ordering]);
      const gra::spin::InternalHelicityMetric metric{numerator.helicities, numerator.matrix};
      return gra::spin::Subchannel(first, second, lts.decaytree[0].p, lts.decaytree[1].p, ordering != 0, false,
                                   &metric);
    }
    return gra::spin::Subchannel(first, second, lts.decaytree[0].p, lts.decaytree[1].p, ordering != 0);
  };
  return std::pair{gra::spin::Contract(upper, lower, central(0)), gra::spin::Contract(upper, lower, central(1))};
}

// Evaluate a hadronic GP continuum through its full exchange-pair tensors
std::pair<MMatrix<std::complex<double>>, MMatrix<std::complex<double>>> FullGPContinuum(
    const gra::LORENTZSCALAR &lts, const gra::gpom::AmpCache &cache, std::size_t channel) {
  const auto &exchange = lts.process.CONT_PRODUCTION[channel];
  const auto &vertices = lts.process.CONTINUUM_GP[channel];
  const auto  source   = std::find_if(cache.sources.cbegin(), cache.sources.cend(), [&](const auto &entry) {
    return entry.key.up_pdg == exchange[0] && entry.key.dn_pdg == exchange[1] &&
           entry.key.MMAX == vertices[0].analytic_MMAX;
  });
  REQUIRE(source != cache.sources.cend());
  const auto upper      = gra::OuterProduct(source->upper_row_factor, source->upper_column_factor);
  const auto lower      = gra::OuterProduct(source->lower_row_factor, source->lower_column_factor);
  auto       lower_axis = lts.q2;
  lower_axis.Flip3();
  // Retain every exchange-pair row until the final forward contraction
  const auto central = [&](std::size_t ordering) {
    const auto &up = vertices[2 * ordering];
    const auto &dn = vertices[2 * ordering + 1];
    return gra::spin::Subchannel(gra::gpom::Crossed(up, source->basis.alpha1, lts.q1, false),
                                 gra::gpom::Crossed(dn, source->basis.alpha2, lower_axis, true), lts.decaytree[0].p,
                                 lts.decaytree[1].p, ordering != 0, true);
  };
  return std::pair{gra::spin::Contract(upper, lower, central(0)), gra::spin::Contract(upper, lower, central(1))};
}

// Configure primary and secondary exchange pairs for continuum channel tests
void SetSecondaryContinuumChannels(nlohmann::json &card) {
  constexpr std::array<const char *, 7> pairs = {"[211,-211]", "[321,-321]", "[2212,-2212]", "[111,111]",
                                                 "[311,-311]", "[113,113]",  "[333,333]"};
  for (const char *pair : pairs) {
    card["PARAM_REGGE"]["PARAM_CON"]["MP"][pair] = {{995, 995}, {995, 9915}};
    card["PARAM_REGGE"]["PARAM_CON"]["XP"][pair] = {{995, 995}, {995, 9915}};
    card["PARAM_REGGE"]["PARAM_CON"]["GP"][pair] = {{990, 990}, {990, 9910}};
  }
}

// Invert one local soft trajectory branch for an analytic spin value
double ReggeTransferForAlpha(const gra::regge::Param &param, const int pdg, const double target) {
  const std::size_t trajectory = gra::regge::TrajectoryIndex(param, pdg);
  const auto        exchange   = param.exchanges.at(trajectory).soft_exchange;
  const auto       &soft       = param.soft_model->Exchange(exchange);
  if (!std::isfinite(target) || soft.AlphaPrime() <= 0.0) {
    throw std::invalid_argument("ReggeTransferForAlpha requires a finite target and positive slope");
  }
  const auto residual = [&](const double transfer) { return param.soft_model->Alpha(exchange, transfer) - target; };

  const double center = (target - soft.Alpha0()) / soft.AlphaPrime();
  if (std::abs(residual(center)) <= 1.0e-14) { return center; }
  double span       = std::max(1.0, std::abs(center));
  double low        = center - span;
  double high       = center + span;
  double low_value  = residual(low);
  double high_value = residual(high);
  for (std::size_t step = 0; step < 64 && low_value * high_value > 0.0; ++step) {
    span *= 2.0;
    low        = center - span;
    high       = center + span;
    low_value  = residual(low);
    high_value = residual(high);
  }
  if (!std::isfinite(low_value) || !std::isfinite(high_value) || low_value * high_value > 0.0) {
    throw std::logic_error("ReggeTransferForAlpha could not bracket target");
  }

  for (std::size_t step = 0; step < 160; ++step) {
    const double middle       = 0.5 * (low + high);
    const double middle_value = residual(middle);
    if (std::abs(middle_value) <= 1.0e-14) { return middle; }
    if ((low_value <= 0.0 && middle_value >= 0.0) || (low_value >= 0.0 && middle_value <= 0.0)) {
      high = middle;
    } else {
      low       = middle;
      low_value = middle_value;
    }
  }
  return 0.5 * (low + high);
}

// Compute one populated analytic exchange-spin cache block
const gra::gpom::ReggeCGBlock &GPSpinBlock(const gra::gpom::ForwardSource &source, std::size_t two_s) {
  const auto block = std::find_if(source.regge_cg.cbegin(), source.regge_cg.cend(),
                                  [two_s](const auto &candidate) { return candidate.two_s == two_s; });
  REQUIRE(block != source.regge_cg.cend());
  return *block;
}

// Compute one populated analytic exchange-spin coefficient
const std::complex<double> &GPSpinValue(const gra::gpom::ForwardSource &source, const gra::gpom::ReggeCGBlock &block,
                                        int m1, int m2) {
  const std::size_t i1    = gra::gpom::AnalyticMIndex(m1, source.basis.MMAX, "GP cached upper spin");
  const std::size_t i2    = gra::gpom::AnalyticMIndex(m2, source.basis.MMAX, "GP cached lower spin");
  const auto       &value = block.coefficient[i1 * source.basis.nm + i2];
  REQUIRE(value.has_value());
  return *value;
}

// Build one sparse scalar crossed vertex in the finite Regge m basis
gra::HELMatrix ScalarFiniteMCrossed(const int mmax, const std::vector<std::pair<int, std::complex<double>>> &residue) {
  gra::HELMatrix hel;
  gra::gpom::InitCrossed(hel, 0.0, 0.0, mmax, "scalar finite-m crossed test");
  hel.C_symmetry     = true;
  hel.P_symmetry     = true;
  hel.coupling_basis = gra::CouplingBasis::Helicity;
  for (const auto &[m, value] : residue) {
    const std::size_t column = gra::gpom::AnalyticMIndex(m, mmax, "scalar finite-m crossed test");
    hel.T[0][column]         = value;
    hel.T_set[0][column]     = true;
  }
  gra::gpom::CheckHelicity(hel, "scalar finite-m crossed test", false);
  return hel;
}

// Replace scalar toy subchannels by a common finite Regge m vertex
void UseScalarFiniteMContinuum(gra::LORENTZSCALAR                                      &lts,
                               const std::vector<std::pair<int, std::complex<double>>> &residue) {
  REQUIRE(lts.decaytree.size() == 2);
  REQUIRE(lts.decaytree[0].p.spinX2 == 0);
  REQUIRE(lts.decaytree[1].p.spinX2 == 0);
  lts.process.CONTINUUM_POLE.clear();
  lts.process.CONTINUUM_GP.clear();
  for (const auto &tree : lts.process.CONT_PRODUCTIONTREE) {
    REQUIRE(tree.size() == 2);
    REQUIRE(tree[0].p.pdg != gra::PDG::PDG_gamma);
    REQUIRE(tree[1].p.pdg != gra::PDG::PDG_gamma);
    std::vector<gra::HELMatrix> vertices;
    vertices.reserve(4);
    for (std::size_t vertex = 0; vertex < 4; ++vertex) {
      vertices.push_back(ScalarFiniteMCrossed(lts.process.MMAX, residue));
    }
    lts.process.CONTINUUM_GP.push_back(std::move(vertices));
  }
}

// Apply the Q_m = (-1)^m section carried by a reflected GP continuum leg
gra::LORENTZSCALAR GPOddMSectionForTest(gra::LORENTZSCALAR lts) {
  for (auto &channel : lts.process.CONTINUUM_GP) {
    for (gra::HELMatrix &vertex : channel) {
      if (vertex.exchange_basis != gra::ExchangeBasisType::ReggeHelicity) { continue; }
      for (int m = -vertex.analytic_MMAX; m <= vertex.analytic_MMAX; ++m) {
        if (std::abs(m) % 2 == 0) { continue; }
        const std::size_t column = gra::gpom::AnalyticMIndex(m, vertex.analytic_MMAX, "GP reflected odd-m section");
        if (vertex.UsesLSCouplings()) {
          vertex.m_ls[column].Scale(-1.0);
          continue;
        }
        for (std::size_t row = 0; row < vertex.T.size_row(); ++row) {
          if (vertex.T_set[row][column]) { vertex.T[row][column] *= -1.0; }
        }
      }
    }
  }
  return lts;
}

// Keep physical photon vertices and replace the Regge leg by a finite m vertex
void UseScalarMixedFiniteMContinuum(gra::LORENTZSCALAR                                      &lts,
                                    const std::vector<std::pair<int, std::complex<double>>> &residue) {
  UseToyGPSubchannelHelicity(lts);
  REQUIRE(lts.process.CONTINUUM_GP.size() == 1);
  REQUIRE(lts.process.CONT_PRODUCTIONTREE.size() == 1);
  const auto &tree   = lts.process.CONT_PRODUCTIONTREE.front();
  auto       &vertex = lts.process.CONTINUUM_GP.front();
  REQUIRE(tree.size() == 2);
  REQUIRE(vertex.size() == 4);
  for (const auto &i : indices(vertex)) {
    const bool photon = tree[i % 2].p.pdg == gra::PDG::PDG_gamma;
    if (photon) {
      REQUIRE(vertex[i].exchange_basis == gra::ExchangeBasisType::HelicityTransport);
      REQUIRE(vertex[i].analytic_MMAX == 1);
      continue;
    }
    vertex[i] = ScalarFiniteMCrossed(lts.process.MMAX, residue);
  }
}

// Replace scalar toy subchannels by one explicit m zero vertex
void UseScalarZeroMContinuum(gra::LORENTZSCALAR &lts, const std::complex<double> residue) {
  REQUIRE(lts.decaytree.size() == 2);
  REQUIRE(lts.decaytree[0].p.spinX2 == 0);
  REQUIRE(lts.decaytree[1].p.spinX2 == 0);
  lts.process.CONTINUUM_POLE.clear();
  lts.process.CONTINUUM_GP.clear();
  for (const auto &tree : lts.process.CONT_PRODUCTIONTREE) {
    REQUIRE(tree.size() == 2);
    std::vector<gra::HELMatrix> vertices;
    vertices.reserve(4);
    for (std::size_t vertex = 0; vertex < 4; ++vertex) {
      gra::HELMatrix hel;
      gra::gpom::InitCrossed(hel, 0.0, 0.0, 0, "scalar m zero crossed test");
      hel.C_symmetry     = true;
      hel.P_symmetry     = true;
      hel.coupling_basis = gra::CouplingBasis::Helicity;
      hel.T[0][0]        = residue;
      hel.T_set[0][0]    = true;
      gra::gpom::CheckHelicity(hel, "scalar m zero crossed test", false);
      vertices.push_back(std::move(hel));
    }
    lts.process.CONTINUUM_GP.push_back(std::move(vertices));
  }
}

// Convert one complete pole tensor to independent direct GP representatives
gra::RES_PRODUCTION_CHANNEL DirectGPFromPole(const std::array<int, 2>                 &pair,
                                             const gra::MMatrix<std::complex<double>> &pole) {
  REQUIRE(pole.size_row() % 2 == 1);
  REQUIRE(pole.size_col() % 2 == 1);
  const int  m1max     = static_cast<int>((pole.size_row() - 1) / 2);
  const int  m2max     = static_cast<int>((pole.size_col() - 1) / 2);
  const bool identical = pair[0] == pair[1];
  REQUIRE_FALSE((identical && m1max != m2max));
  gra::RES_PRODUCTION_CHANNEL model;
  model.exchange   = pair;
  model.basis      = gra::ReggeVertexBasis::Helicity;
  model.C_symmetry = true;
  model.P_symmetry = true;
  for (int m1 = -m1max; m1 <= m1max; ++m1) {
    for (int m2 = -m2max; m2 <= m2max; ++m2) {
      const std::pair<int, int> coordinate = {m1, m2};
      std::pair<int, int>       orbit      = std::min(coordinate, std::make_pair(-m1, -m2));
      if (identical) { orbit = std::min({orbit, std::make_pair(m2, m1), std::make_pair(-m2, -m1)}); }
      if (coordinate != orbit) { continue; }
      const std::size_t i1 = static_cast<std::size_t>(m1 + m1max);
      const std::size_t i2 = static_cast<std::size_t>(m2 + m2max);
      if (std::abs(pole[i1][i2]) <= 1.0e-13) { continue; }
      model.helicity.push_back({static_cast<double>(m1), static_cast<double>(m2)});
      model.g_helicity.push_back(pole[i1][i2]);
    }
  }
  return model;
}

// Set one typed LS channel from compact steering rows
void SetChannelLS(gra::RES_PRODUCTION_CHANNEL &channel, const nlohmann::json &rows) {
  channel.basis = gra::ReggeVertexBasis::LS;
  channel.g_ls.Clear();
  channel.helicity.clear();
  channel.g_helicity.clear();
  for (const auto &row : rows) {
    channel.g_ls.Set(row[0].get<std::size_t>(), static_cast<std::size_t>(std::llround(2.0 * row[1].get<double>())),
                     std::polar(row[2].get<double>(), row[3].get<double>()));
  }
}

// Set one typed helicity channel from compact steering rows
void SetChannelHelicity(gra::RES_PRODUCTION_CHANNEL &channel, const nlohmann::json &rows) {
  channel.basis = gra::ReggeVertexBasis::Helicity;
  channel.g_ls.Clear();
  channel.helicity.clear();
  channel.g_helicity.clear();
  for (const auto &row : rows) {
    channel.helicity.push_back({row[0].get<double>(), row[1].get<double>()});
    channel.g_helicity.push_back(std::polar(row[2].get<double>(), row[3].get<double>()));
  }
}

}  // namespace

TEST_CASE("Forward reference density accepts fermion helicities", "[gra::MRegge][normalization][fermion]") {
  const gra::MMatrix<std::complex<double>> tensor = {
      std::vector<std::complex<double>>{1.0, 0.0}, std::vector<std::complex<double>>{0.0, 2.0},
      std::vector<std::complex<double>>{3.0, 0.0}, std::vector<std::complex<double>>{0.0, 4.0}};
  const gra::MMatrix<double> helicity = {std::vector<double>{-0.5, -0.5}, std::vector<double>{-0.5, 0.5},
                                         std::vector<double>{0.5, -0.5}, std::vector<double>{0.5, 0.5}};
  const double density = gra::spin::ForwardHelicityDensity(tensor, helicity, gra::spin::ForwardLegType::Hadron,
                                                           gra::spin::ForwardLegType::Hadron);
  CHECK(density == Approx(7.5));
}

TEST_CASE("Continuum metrics cover physical poles and reduced Regge residues",
          "[gra::MRegge][continuum][spin-metric][GP]") {
  for (int spin = 0; spin <= 4; ++spin) {
    const auto projections = gra::spin::SpinProjections(spin);
    {
      const auto fixed = gra::spin::ExchangeMetric(gra::ExchangeBasisType::HelicityTransport, projections, spin);
      for (const auto &row : gra::aux::indices(projections)) {
        for (const auto &column : gra::aux::indices(projections)) {
          const bool                 paired = std::abs(projections[row] + projections[column]) < 1.0e-12;
          const std::complex<double> expected =
              paired ? gra::spin::JacobWickSecondLegReversalPhase(spin, projections[row]) : 0.0;
          RequireComplexNear(fixed[row][column], expected, 1.0e-12);
        }
      }
    }
  }

  {
    const auto reduced = gra::spin::ExchangeMetric(gra::ExchangeBasisType::ReducedRegge, {0.0}, 2.0);
    REQUIRE(reduced.size_row() == 1);
    REQUIRE(reduced.size_col() == 1);
    RequireComplexNear(reduced[0][0], 1.0, 1.0e-12);
  }
  const auto gp_reduced = gra::spin::ExchangeMetric(gra::ExchangeBasisType::ReggeHelicity, {0.0}, 2.0);
  REQUIRE(gp_reduced.size_row() == 1);
  REQUIRE(gp_reduced.size_col() == 1);
  RequireComplexNear(gp_reduced[0][0], 1.0, 1.0e-12);
  REQUIRE_THROWS_AS(gra::spin::ExchangeMetric(gra::ExchangeBasisType::ReggeHelicity, {-1.0, 0.0, 1.0}, 1.0),
                    std::invalid_argument);
}

// Check a model, exchange and pair transfer override wins over the common form
TEST_CASE("Continuum transfer form selects the exact production-card override",
          "[gra::MRegge][continuum][form-factor]") {
  const std::string tune          = WriteModifiedContinuumTune("transfer_override_precedence", "GP", [](auto &card) {
    card.at("990").at("[211,211]")["FF_transfer"] = {
        {"type", "power"}, {"norm", "zero"}, {"Lambda2", 0.25}, {"n", 1.0}};
  });
  const std::string general       = tune + "/GENERAL.json";
  const auto        tune_snapshot = gra::MModelTune::Load(general);
  const auto        param         = gra::regge::ReadParam({211, -211}, LoadedPDGTable(), *tune_snapshot);

  const std::vector<int> pion_pair = {211, -211};
  const auto             &gp_pair  = gra::regge::Pair(param, pion_pair, gra::ReggeProductionModel::GP);
  const auto channel = std::find_if(gp_pair.channels.cbegin(), gp_pair.channels.cend(), [](const auto &value) {
    return value.first == 990 && value.second == 990;
  });
  REQUIRE(channel != gp_pair.channels.cend());
  const auto &exact = channel->transfer[0];
  REQUIRE(exact.type == gra::regge::FFType::Power);
  REQUIRE(exact.norm == gra::regge::FFNorm::Zero);
  REQUIRE(exact.param.size() == 2);
  CHECK(exact.param[0] == Approx(0.25));
  CHECK(exact.param[1] == Approx(1.0));

  const auto &mp_pair = gra::regge::Pair(param, pion_pair, gra::ReggeProductionModel::MP);
  const auto baseline = gra::regge::ReadParam(pion_pair, LoadedPDGTable(), *gra::MModelTune::Load(modelfile));
  const auto &mp_baseline = gra::regge::Pair(baseline, pion_pair, gra::ReggeProductionModel::MP);
  REQUIRE_FALSE(mp_pair.channels.empty());
  CHECK(mp_pair.channels.front().transfer[0].type == mp_baseline.channels.front().transfer[0].type);
  const std::vector<int> proton_pdgs = {2212, -2212};
  const auto            &proton_pair = gra::regge::Pair(param, proton_pdgs, gra::ReggeProductionModel::GP);
  const auto& proton_baseline = gra::regge::Pair(baseline, proton_pdgs, gra::ReggeProductionModel::GP);
  REQUIRE(proton_pair.channels.size() == proton_baseline.channels.size());
  for (const auto& i : gra::aux::indices(proton_pair.channels)) {
    CHECK(proton_pair.channels[i].transfer[0].type == proton_baseline.channels[i].transfer[0].type);
  }
}

TEST_CASE("GP nonsense factors have anchor-free natural residue zeros", "[gra::MRegge][continuum][spin-metric][GP]") {
  CHECK(gra::gpom::NonsenseZero(1.37, 0) == Approx(1.0));
  CHECK(gra::gpom::NonsenseZero(3.0, 3) == Approx(6.0));
  CHECK(gra::gpom::NonsenseZero(1.5, 3) == Approx(-0.375));
  CHECK(gra::gpom::NonsenseZero(1.5, -3) == Approx(-0.375));
  for (int zero = 0; zero < 3; ++zero) {
    CHECK(gra::gpom::NonsenseZero(static_cast<double>(zero), 3) == Approx(0.0).margin(1.0e-15));
  }
  CHECK_THROWS_AS(gra::gpom::NonsenseZero(std::numeric_limits<double>::quiet_NaN(), 1), gra::AmplitudeFailure);
}

TEST_CASE("Analytic trajectory aliases resolve canonical physical poles", "[gra::MRegge][GP][normalization][params]") {
  const auto pdg   = LoadedPDGTable();
  const auto tune  = gra::MModelTune::Load(modelfile);
  const auto param = gra::regge::ReadParam({211, -211}, pdg, *tune);
  CHECK(gra::regge::PoleRepresentative(param, pdg, 990).pdg == 995);
  CHECK(gra::regge::PoleRepresentative(param, pdg, 9910).pdg == 9915);
  CHECK(gra::regge::PoleRepresentative(param, pdg, 9930).pdg == 9933);
  CHECK(gra::regge::PoleRepresentative(param, pdg, 9990).pdg == 9993);
  CHECK(gra::gpom::PoleSpinX2(param, pdg, 990) == 4);
  CHECK(gra::gpom::PoleSpinX2(param, pdg, 9990) == 2);
}

TEST_CASE("GP raw coefficients do not depend on MMAX", "[gra::MRegge][GP][normalization]") {
  const std::complex<double> coupling = std::polar(0.73, -0.41);
  for (const int mmax : {1, 4}) {
    gra::HELMatrix hel;
    gra::gpom::InitHelicity(hel, 0.0, 0.0, mmax, "MMAX test");
    hel.coupling_basis      = gra::CouplingBasis::Helicity;
    const std::size_t zero  = gra::gpom::AnalyticMIndex(0, mmax, "MMAX test");
    const std::size_t empty = gra::gpom::AnalyticMIndex(-mmax, mmax, "MMAX zero row test");
    const std::size_t small = gra::gpom::AnalyticMIndex(mmax, mmax, "MMAX small row test");
    hel.T[0][zero]          = coupling;
    hel.T_set[0][zero]      = true;
    hel.T[0][empty]         = 0.0;
    hel.T_set[0][empty]     = true;
    hel.T[0][small]         = std::numeric_limits<double>::min();
    hel.T_set[0][small]     = true;
    gra::gpom::CheckHelicity(hel, "MMAX test", false);
    RequireComplexNear(hel.T[0][zero], coupling, 1.0e-14);
    REQUIRE_FALSE(hel.T_set[0][empty]);
    REQUIRE(hel.T_set[0][small]);
    REQUIRE(hel.T_active == std::vector<std::pair<std::size_t, std::size_t>>{{0, zero}, {0, small}});
  }
  CHECK_THROWS_AS(gra::gpom::AnalyticMIndex(0, (std::numeric_limits<int>::max() - 1) / 2 + 1, "oversized MMAX test"),
                  std::invalid_argument);
  CHECK_THROWS_AS(gra::gpom::AnalyticMIndex(std::numeric_limits<int>::min(), 2, "minimum m test"),
                  std::invalid_argument);
}

TEST_CASE("Scalar continuum ladders are one-dimensional matrix chains", "[gra::MRegge][continuum][spin-chain]") {
  std::vector<gra::MMatrix<std::complex<double>>> kernels;
  kernels.emplace_back(1, 1, 2.0);
  kernels.emplace_back(1, 1, 3.0);
  std::vector<gra::MMatrix<std::complex<double>>> metrics;
  metrics.emplace_back(1, 1, 5.0);
  const std::complex<double> amplitude = gra::ladder::Contract({7.0}, kernels, metrics, {11.0});
  RequireComplexNear(amplitude, 2310.0, 1.0e-12);
}

TEST_CASE("Regge parameters construct safely from one tune across threads", "[gra::MRegge][params][threading]") {
  const auto       model_tune = gra::MModelTune::Load(modelfile);
  const auto       pdg_table  = LoadedPDGTable();
  gra::MModelCache cache(model_tune);

  constexpr std::size_t                                 nthreads = 8;
  std::vector<std::shared_ptr<const gra::regge::Param>> handles(nthreads);
  std::vector<std::thread>                              workers;
  workers.reserve(nthreads);
  std::atomic<std::size_t> ready = 0;
  std::atomic<bool>        start = false;

  // Start every worker on the same genuinely cold cache key
  for (std::size_t i = 0; i < nthreads; ++i) {
    workers.emplace_back([i, &handles, &pdg_table, &ready, &cache, &start] {
      ready.fetch_add(1, std::memory_order_release);
      while (!start.load(std::memory_order_acquire)) { std::this_thread::yield(); }
      handles[i] = gra::regge::GetParam(cache, {211, -211}, pdg_table);
    });
  }
  while (ready.load(std::memory_order_acquire) != nthreads) { std::this_thread::yield(); }
  start.store(true, std::memory_order_release);
  for (auto &worker : workers) { worker.join(); }
  const auto first = handles.front();
  for (const auto &handle : handles) { REQUIRE(handle == first); }
  REQUIRE(first->con.form_pdgs == std::vector<int>{211, -211});
  REQUIRE(first->soft_model == model_tune->Soft());
  REQUIRE_FALSE(first->exchanges.empty());
  REQUIRE(first->pomeron_trajectory < first->exchanges.size());
  const auto kaons = gra::regge::GetParam(cache, {321, -321}, pdg_table);
  REQUIRE(kaons != first);
  REQUIRE(kaons->con.form_pdgs == std::vector<int>{321, -321});
  REQUIRE(gra::GetReggeNumerics(cache) == gra::GetReggeNumerics(cache));
}

TEST_CASE("MRegge master setup initializes parameters before worker creation", "[gra::MRegge][params][lifecycle]") {
  const std::array<std::string, 3> families = {"MP", "XP", "GP"};
  MRandom                          rng;

  for (std::size_t i = 0; i < families.size(); ++i) {
    const std::string &family = families[i];
    CAPTURE(family);
    const double s0         = 0.83 + 0.07 * static_cast<double>(i);
    const auto   tune       = WriteModifiedPhotoVMTune("regge_master_parameters_" + family,
                                                       [s0](auto &card) { card["PARAM_REGGE"]["s0"] = s0; });
    const auto   soft_model = gra::MModelTune::Load(tune.second);

    ToyHelicityProcess master;
    ConfigureToyProductionProcess(master, family, "RES", "pi+ pi-");
    master.SetModelTune(soft_model);
    const auto f0 = gra::resonance::Read("RES/f0_500.json", rng, gra::ParseReggeProductionModel(family));
    master.SetResonances({{"f0_500", f0}});

    REQUIRE_NOTHROW(master.InitializeProcessAmplitude());
    REQUIRE_NOTHROW(gra::MRegge(master.state.lts, soft_model,
                                gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Resonance, family + "_RES")));
  }
}

TEST_CASE("Regge parameters use immutable tune and PDG snapshots", "[gra::MRegge][params][snapshot]") {
  const auto first_tune =
      WriteModifiedPhotoVMTune("regge_snapshot_reload", [](auto &card) { card["PARAM_REGGE"]["s0"] = 0.91; });
  const auto first_model = gra::MModelTune::Load(first_tune.second);
  const auto first_pdg   = LoadedPDGTable();
  const auto first       = gra::regge::ReadParamPtr({211, -211}, first_pdg, *first_model);
  REQUIRE(first->s0 == Approx(0.91));

  const auto second_tune =
      WriteModifiedPhotoVMTune("regge_snapshot_reload", [](auto &card) { card["PARAM_REGGE"]["s0"] = 1.19; });
  const auto second_model = gra::MModelTune::Load(second_tune.second);
  REQUIRE(second_model->GeneralFile() == first_model->GeneralFile());
  REQUIRE(second_model != first_model);
  const auto second = gra::regge::ReadParamPtr({211, -211}, first_pdg, *second_model);
  REQUIRE(second != first);
  CHECK(first->s0 == Approx(0.91));
  CHECK(second->s0 == Approx(1.19));

  auto changed_pdg = first_pdg;
  changed_pdg.PDG_table.at(211).mass += 1.0e-3;
  const auto changed = gra::regge::ReadParamPtr({211, -211}, changed_pdg, *second_model);
  CHECK(changed != second);
  CHECK(changed->s0 == Approx(1.19));
}

TEST_CASE("MRegge process validates continuum event caches", "[gra::MRegge][process][continuum]") {
  SECTION("multi-body ladder cache") {
    constexpr auto     mode = gra::MReggeMode::ContinuumTwoFourSixBody;
    gra::LORENTZSCALAR lts;
    lts.decaytree.resize(4);
    for (std::size_t i = 0; i < lts.decaytree.size(); ++i) {
      lts.decaytree[i].p.pdg = i % 2 == 0 ? gra::PDG::PDG_pip : gra::PDG::PDG_pim;
    }

    REQUIRE_NOTHROW(gra::MRegge::ValidateInitializedProcess(mode, lts, "MP_CON"));

    lts.process.CONT_LADDER_PERMUTATIONS = {{3, 4, 5}};
    CHECK_THROWS_AS(gra::MRegge::ValidateInitializedProcess(mode, lts, "MP_CON"), std::invalid_argument);
    lts.process.CONT_LADDER_PERMUTATIONS = {{3, 4, 5, 5}};
    CHECK_THROWS_AS(gra::MRegge::ValidateInitializedProcess(mode, lts, "MP_CON"), std::invalid_argument);
    lts.process.CONT_LADDER_PERMUTATIONS = {{3, 4, 5, 6}};
    REQUIRE_NOTHROW(gra::MRegge::ValidateInitializedProcess(mode, lts, "MP_CON"));
    lts.process.CONT_LADDER_PERMUTATIONS = {{3, 4, 5, 6}, {3, 4, 5, 6}};
    CHECK_THROWS_AS(gra::MRegge::ValidateInitializedProcess(mode, lts, "MP_CON"), std::invalid_argument);

    lts.process.CONT_LADDER_PERMUTATIONS = {{3, 4, 5, 6}};
    lts.process.MULTIREGGE_TOPOLOGIES    = {{4}, {2, 2}};
    REQUIRE_NOTHROW(gra::MRegge::ValidateInitializedProcess(mode, lts, "MP_CON"));
    for (const std::vector<gra::regge::Topology> &malformed :
         {std::vector<gra::regge::Topology>{}, std::vector<gra::regge::Topology>{{}},
          std::vector<gra::regge::Topology>{{2, 4}}, std::vector<gra::regge::Topology>{{2}},
          std::vector<gra::regge::Topology>{{3, 1}}, std::vector<gra::regge::Topology>{{4}, {4}}}) {
      CAPTURE(malformed);
      lts.process.MULTIREGGE_TOPOLOGIES = malformed;
      if (malformed.empty()) {
        REQUIRE_NOTHROW(gra::MRegge::ValidateInitializedProcess(mode, lts, "MP_CON"));
      } else {
        CHECK_THROWS_AS(gra::MRegge::ValidateInitializedProcess(mode, lts, "MP_CON"), std::invalid_argument);
      }
    }
  }

  SECTION("two-body production and helicity caches") {
    constexpr auto     mode = gra::MReggeMode::ContinuumTwoBody;
    gra::LORENTZSCALAR lts;
    lts.decaytree.resize(2);
    lts.process.CONT_PRODUCTION     = {{991, 991}};
    lts.process.CONT_PRODUCTIONTREE = {std::vector<gra::MDecayBranch>(2)};
    lts.process.CONTINUUM_POLE      = {std::vector<gra::spin::PoleResidue>(4)};
    CHECK_THROWS_AS(gra::MRegge::ValidateInitializedProcess(mode, lts, "MP_CON"), std::invalid_argument);
    lts.process.CONTINUUM_POLE = MakeToyScalarContinuumLTSAsymmetric(0.3, 4.2, -4.6).process.CONTINUUM_POLE;
    REQUIRE_NOTHROW(gra::MRegge::ValidateInitializedProcess(mode, lts, "MP_CON"));

    lts.process.CONT_PRODUCTIONTREE[0].pop_back();
    CHECK_THROWS_AS(gra::MRegge::ValidateInitializedProcess(mode, lts, "MP_CON"), std::invalid_argument);
    lts.process.CONT_PRODUCTIONTREE[0].resize(2);
    lts.process.CONTINUUM_POLE[0].pop_back();
    CHECK_THROWS_AS(gra::MRegge::ValidateInitializedProcess(mode, lts, "MP_CON"), std::invalid_argument);
  }
}

TEST_CASE("Direct gamma-Pomeron helicity basis enforces P and C constraints", "[gra::spin][helicity]") {
  ToyHelicityProcess proc;

  const auto pomeron0 = ToyParticle(991, 0, 1, 1, "P0");
  const auto photon   = ToyParticle(22, 2, -1, -1, "gamma");
  const auto rho0     = ToyParticle(113, 2, -1, -1, "rho0");

  auto hel = DirectRhoGammaPomeronMatrix();
  REQUIRE_NOTHROW(
      gra::spin::ValidateDirectTMatrix(hel, rho0, photon, pomeron0, true, "direct helicity test", true, false));
  REQUIRE(hel.UsesHelicityCouplings());
  REQUIRE(hel.T.size_row() == 3);
  REQUIRE(hel.T.size_col() == 1);

  const double schc = 1.0 / std::sqrt(2.0);
  CHECK(std::abs(hel.T[0][0] - schc) < 1e-12);
  CHECK(std::abs(hel.T[2][0] - schc) < 1e-12);
  CHECK(hel.T_set[0][0]);
  CHECK(hel.T_set[2][0]);

  auto parity_breaking = hel;
  parity_breaking.T[2][0] *= -1.0;
  REQUIRE_THROWS(
      gra::spin::ValidateDirectTMatrix(parity_breaking, rho0, photon, pomeron0, true, "CON_MP.json", true, false));

  auto c_forbidden = rho0;
  c_forbidden.C    = 1;
  auto c_bad       = hel;
  REQUIRE_THROWS(
      gra::spin::ValidateDirectTMatrix(c_bad, c_forbidden, photon, pomeron0, true, "CON_MP.json", true, false));
}

TEST_CASE(
    "MP cards select generated central couplings per resonance and "
    "continuum vertex",
    "[gra::spin][MP][RES][CON]") {
  gra::MODELPARAM = "TUNE0";
  MRandom                        rng;
  const std::vector<std::string> modes = {"auto_min_L", "auto_min_S", "auto_equal_ls", "auto_equal_helicity"};

  SECTION("MP central fusion") {
    for (const auto &mode : modes) {
      CAPTURE(mode);
      ToyHelicityProcess proc;
      ConfigureToyProductionProcess(proc, "MP", "RES", "pi+ pi-");

      auto  f0         = gra::resonance::Read("RES/f0_500.json", rng, gra::ReggeProductionModel::MP);
      SetMPFusion(f0);
      auto &channel    = f0.MP.channels.front();
      channel.exchange = {993, 993};
      channel.basis    = gra::ParseReggeVertexBasis(mode, gra::ReggeVertexRole::Resonance);
      const std::complex<double> expected_coefficient = channel.g;
      proc.SetResonances({{"f0_500", f0}});

      REQUIRE_NOTHROW(proc.InitializeProcessAmplitude());
      const auto resonances = proc.GetResonances();
      REQUIRE(resonances.count("f0_500") == 1);
      const auto &configured = resonances.at("f0_500");
      REQUIRE(configured.production.size() == 1);
      const auto &hel    = configured.production.front().hel;
      const auto &vertex = configured.production.front().pole.value();
      RequireCanonicalPoleVertex(vertex);
      if (mode == "auto_min_L" || mode == "auto_min_S") {
        REQUIRE(vertex.terms.size() == 1);
        CHECK(vertex.terms.front().l == 0);
        CHECK(vertex.terms.front().two_s == 0);
      } else if (mode == "auto_equal_ls") {
        REQUIRE(vertex.terms.size() > 1);
        for (const auto &term : vertex.terms) { RequireComplexNear(term.coefficient, expected_coefficient, 1.0e-14); }
      }
      CHECK_FALSE(hel.UsesHelicityCouplings());
    }
  }

  SECTION("MP explicit central fusion bases") {
    for (const std::string basis : {"g_ls", "helicity"}) {
      CAPTURE(basis);
      ToyHelicityProcess proc;
      ConfigureToyProductionProcess(proc, "MP", "RES", "pi+ pi-");

      auto  f0         = gra::resonance::Read("RES/f0_500.json", rng, gra::ReggeProductionModel::MP);
      SetMPFusion(f0);
      auto &channel    = f0.MP.channels.front();
      channel.exchange = {993, 993};
      if (basis == "g_ls") {
        const auto pole = proc.state.lts.PDG.FindByPDG(993);
        const auto operators =
            gra::spin::CanonicalPoleOperators(f0.p, pole, pole, true, channel.C_symmetry, channel.P_symmetry);
        REQUIRE_FALSE(operators.empty());
        channel.basis = gra::ReggeVertexBasis::LS;
        channel.g_ls.Clear();
        for (const auto &i : indices(operators)) {
          const auto &term = operators[i].coupling;
          channel.g_ls.Set(term.l, term.two_s, i == 0 ? 1.0 : 0.0);
        }
      } else {
        SetChannelHelicity(channel, {{-1, -1, 1.0, 0.0}, {0, 0, 1.0, gra::math::PI}});
      }
      proc.SetResonances({{"f0_500", f0}});

      REQUIRE_NOTHROW(proc.InitializeProcessAmplitude());
      const auto &configured = proc.GetResonances().at("f0_500");
      REQUIRE(configured.production.size() == 1);
      const auto &hel = configured.production.front().hel;
      RequireCanonicalPoleVertex(configured.production.front().pole.value());
      CHECK_FALSE(hel.UsesHelicityCouplings());
    }
  }

  SECTION("MP tensor central fusion normalization") {
    for (const auto &mode : modes) {
      CAPTURE(mode);
      ToyHelicityProcess proc;
      ConfigureToyProductionProcess(proc, "MP", "RES", "pi+ pi-");

      auto  f2         = gra::resonance::Read("RES/f2_1270.json", rng, gra::ReggeProductionModel::MP);
      f2.spin_basis = "none";
      f2.a_Jz.clear();
      auto &channel    = f2.MP.channels.front();
      channel.exchange = {993, 993};
      channel.basis    = gra::ParseReggeVertexBasis(mode, gra::ReggeVertexRole::Resonance);
      const std::complex<double> expected_coefficient = channel.g;
      proc.SetResonances({{"f2_1270", f2}});

      REQUIRE_NOTHROW(proc.InitializeProcessAmplitude());
      const auto resonances = proc.GetResonances();
      REQUIRE(resonances.count("f2_1270") == 1);
      const auto &configured = resonances.at("f2_1270");
      REQUIRE(configured.production.size() == 1);
      const auto &hel = configured.production.front().hel;
      REQUIRE(configured.production.size() == 1);
      const auto &vertex = configured.production.front().pole.value();
      RequireCanonicalPoleVertex(vertex);
      CHECK_FALSE(hel.UsesHelicityCouplings());
      if (mode == "auto_min_L") {
        REQUIRE(vertex.terms.size() == 1);
        CHECK(vertex.terms.front().l == 0);
        CHECK(vertex.terms.front().two_s == 4);
      } else if (mode == "auto_min_S") {
        REQUIRE(vertex.terms.size() == 1);
        CHECK(vertex.terms.front().l == 2);
        CHECK(vertex.terms.front().two_s == 0);
      } else if (mode == "auto_equal_ls") {
        REQUIRE(vertex.terms.size() > 1);
        for (const auto &term : vertex.terms) { REQUIRE(term.coefficient == expected_coefficient); }
      }

      CHECK(configured.UsesUnrestrictedSpinBasis());
      CHECK(configured.a_Jz.empty());
      CHECK(configured.MP.filter.isEmpty());
    }
  }

  SECTION("MP t/u subchannel vertices") {
    for (const auto &mode : modes) {
      CAPTURE(mode);
      const std::string  tune = WriteModifiedContinuumTune("mp_generated_" + mode, "MP", [&mode](auto &card) {
        for (auto &[exchange, pairs] : card.items()) {
          (void)exchange;
          for (auto &[pair_key, block] : pairs.items()) {
            (void)pair_key;
            for (const std::string sector : {"same", "opposite", "self"}) {
              if (block.contains(sector)) { block.at(sector)["basis"] = "crossed_" + mode; }
            }
          }
        }
      });
      ToyHelicityProcess proc;
      ConfigureToyProductionProcess(proc, "MP", "CON", "pi+ pi-");
      proc.SetTuneForTest(tune);

      REQUIRE_NOTHROW(proc.InitializeProcessAmplitude());
      REQUIRE(proc.state.lts.process.CONT_TU_SIGN.size() == proc.state.lts.process.CONT_PRODUCTION.size());
      for (const double sign : proc.state.lts.process.CONT_TU_SIGN) { REQUIRE(sign == Approx(1.0).margin(1.0e-15)); }
      REQUIRE(proc.state.lts.process.CONTINUUM_POLE.size() == proc.state.lts.process.CONT_PRODUCTION.size());
      for (const auto &channel : proc.state.lts.process.CONTINUUM_POLE) {
        REQUIRE(channel.size() == 4);
        for (const auto &vertex : channel) { RequireCanonicalPoleVertex(vertex.Pole()); }
      }
    }
  }

  SECTION("MP identical scalar pairs use the matched spin-two Pomeron") {
    const auto         tune = WriteModifiedPhotoVMTune("mp_identical_scalar_pair", [](auto &card) {
      card.at("PARAM_REGGE").at("PARAM_CON").at("MP")["[9000221,9000221]"] = {{995, 995}};
    });
    ToyHelicityProcess proc;
    ConfigureToyProductionProcess(proc, "MP", "CON", "f(0)(500)0 f(0)(500)0");
    proc.SetTuneForTest(tune.first);

    REQUIRE_NOTHROW(proc.InitializeProcessAmplitude());
    REQUIRE(proc.state.lts.process.CONT_PRODUCTION == std::vector<std::vector<int>>{{995, 995}});
    REQUIRE(proc.state.lts.process.CONTINUUM_POLE.size() == 1);
    for (const auto &residue : proc.state.lts.process.CONTINUUM_POLE[0]) {
      const auto &vertex = residue.Pole();
      REQUIRE(vertex.terms.size() == 1);
      CHECK(vertex.terms.front().l == 2);
      CHECK(vertex.terms.front().two_s == 0);
    }
  }

  SECTION("fermion-antifermion subvertices retain the physical opposite sector") {
    for (const std::string model : {"MP", "XP", "GP"}) {
      CAPTURE(model);
      ToyHelicityProcess proc;
      ConfigureToyProductionProcess(proc, model, "CON", "p+ p-");

      REQUIRE_NOTHROW(proc.InitializeProcessAmplitude());
      if (model == "GP") {
        const auto &channels = proc.state.lts.process.CONTINUUM_GP;
        REQUIRE_FALSE(channels.empty());
        for (const auto &vertices : channels) {
          REQUIRE(vertices.size() == 4);
          for (const auto &vertex : vertices) {
            CHECK(vertex.UsesReggeDomain());
            CHECK(vertex.exchange_basis == gra::ExchangeBasisType::ReggeHelicity);
            const auto pole = gra::gpom::Crossed(vertex, vertex.J, gra::M4Vec{}, false).helicity;
            CHECK(pole.T.IsFinite());
            CHECK(pole.T.MaskedSquaredNorm(pole.T_set) > 0.0);
            const std::size_t zero = gra::gpom::AnalyticMIndex(0, vertex.analytic_MMAX, "physical p pbar parity test");
            std::size_t       active = 0;
            for (std::size_t row = 0; row < pole.T_set.size_row(); ++row) {
              if (!pole.T_set[row][zero]) { continue; }
              ++active;
              const double      lambda1 = vertex.lambda_values[row][0];
              const double      lambda2 = vertex.lambda_values[row][1];
              const std::size_t partner = gra::spin::DirectHelicityCoordinateIndex(
                  -lambda1, -lambda2, vertex.s1, vertex.s2, "physical p pbar parity test");
              REQUIRE(pole.T_set[partner][zero]);
              RequireComplexNear(pole.T[partner][zero], pole.T[row][zero], 1.0e-12);
            }
            CHECK(active == 4);
          }
        }
        continue;
      }
      REQUIRE_FALSE(proc.state.lts.process.CONTINUUM_POLE.empty());
      for (const auto &channel : indices(proc.state.lts.process.CONTINUUM_POLE)) {
        const auto &vertices = proc.state.lts.process.CONTINUUM_POLE[channel];
        const auto &exchange = proc.state.lts.process.CONT_PRODUCTION[channel];
        REQUIRE(vertices.size() == 4);
        for (const auto &vertex_index : indices(vertices)) {
          const auto &vertex       = vertices[vertex_index].Pole();
          const int   exchange_pdg = exchange[vertex_index % 2];
          const bool  c_odd        = exchange_pdg == 9933 || exchange_pdg == 9930;
          bool        active       = false;
          for (const auto &term : vertex.terms) {
            if (std::abs(term.coefficient) <= 0.0) { continue; }
            active = true;
            CHECK(term.l == (c_odd ? 0 : 1));
            CHECK(term.two_s == 2);
          }
          CHECK(active);
        }
      }
    }
  }

  SECTION("invalid typed selector is rejected") {
    ToyHelicityProcess proc;
    ConfigureToyProductionProcess(proc, "MP", "RES", "pi+ pi-");
    auto f0                      = gra::resonance::Read("RES/f0_500.json", rng, gra::ReggeProductionModel::MP);
    f0.MP.channels.front().basis = static_cast<gra::ReggeVertexBasis>(999);
    proc.SetResonances({{"f0_500", f0}});
    REQUIRE_THROWS_AS(proc.InitializeProcessAmplitude(), std::invalid_argument);
  }
}

TEST_CASE("XP accepts and GP rejects generated resonance couplings", "[gra::spin][XP][GP][RES]") {
  gra::MODELPARAM                      = "TUNE0";
  const std::vector<std::string> modes = {"auto_min_L", "auto_min_S", "auto_equal_ls", "auto_equal_helicity"};

  for (const std::string model : {"XP", "GP"}) {
    for (const auto &mode : modes) {
      CAPTURE(model, mode);
      ToyHelicityProcess proc;
      ConfigureToyProductionProcess(proc, model, "RES", "pi+ pi-");
      auto f2 = gra::resonance::Read("RES/f2_1270.json", proc.state.random, gra::ParseReggeProductionModel(model));

      if (model == "XP") {
        auto &channel = f2.XP.channels.front();
        channel.basis = gra::ParseReggeVertexBasis(mode, gra::ReggeVertexRole::Resonance);
        channel.g     = 0.2879;
        channel.g_ls.Clear();
      } else {
        auto &channel = f2.GP.channels.front();
        channel.basis = gra::ParseReggeVertexBasis(mode, gra::ReggeVertexRole::Resonance);
        channel.g     = 0.2879;
        channel.g_ls.Clear();
      }
      proc.SetResonances({{"f2_1270", f2}});

      if (model == "GP") {
        REQUIRE_THROWS_AS(proc.InitializeProcessAmplitude(), std::invalid_argument);
        continue;
      }
      REQUIRE_NOTHROW(proc.InitializeProcessAmplitude());
      const auto &configured = proc.GetResonances().at("f2_1270");
      REQUIRE(configured.production.size() == 1);
      REQUIRE(configured.production.size() == 1);
      const auto &vertex = configured.production.front().pole.value();
      RequireCanonicalPoleVertex(vertex);
      if (mode == "auto_min_L" || mode == "auto_min_S") {
        const std::size_t l     = mode == "auto_min_L" ? 0 : 2;
        const std::size_t two_s = mode == "auto_min_L" ? 4 : 0;
        REQUIRE(vertex.terms.size() == 1);
        CHECK(vertex.terms.front().l == l);
        CHECK(vertex.terms.front().two_s == two_s);
        CHECK(std::abs(vertex.terms.front().coefficient) > 0.0);
      }
    }
  }

  for (const std::string model : {"XP", "GP"}) {
    for (const auto &mode : modes) {
      CAPTURE(model, mode);
      const std::string tune =
          WriteModifiedContinuumTune("generated_" + model + "_" + mode, model, [&mode](auto &card) {
            for (auto &[exchange, pairs] : card.items()) {
              (void)exchange;
              for (auto &[pair_key, block] : pairs.items()) {
                (void)pair_key;
                for (const std::string sector : {"same", "opposite", "self"}) {
                  if (!block.contains(sector)) { continue; }
                  auto &vertex    = block.at(sector);
                  vertex["basis"] = "crossed_" + mode;
                  vertex.erase("g_ls");
                  vertex.erase("helicity");
                  vertex["g"] = {1.0, 0.0};
                }
              }
            }
          });
      ToyHelicityProcess proc;
      ConfigureToyProductionProcess(proc, model, "CON", "pi+ pi-");
      proc.SetTuneForTest(tune);

      if (model == "GP") {
        REQUIRE_THROWS_AS(proc.InitializeProcessAmplitude(), std::invalid_argument);
        continue;
      } else {
        REQUIRE_NOTHROW(proc.InitializeProcessAmplitude());
      }
      REQUIRE_FALSE(proc.state.lts.process.CONTINUUM_POLE.empty());
      for (const auto &channel : proc.state.lts.process.CONTINUUM_POLE) {
        for (const auto &vertex : channel) { RequireCanonicalPoleVertex(vertex.Pole()); }
      }
    }
  }
}

// Exercise sign steering through initialization and complex continuum amplitudes
TEST_CASE("MRegge resolves continuum t/u signs during initialization",
          "[gra::MRegge][process][continuum][validation][tu-sign]") {
  ModelParamRestoreGuard restore;
  const std::string family = GENERATE("MP", "XP", "GP");
  const std::string final = GENERATE("pi+ pi-", "p+ p-", "rho(770)0 rho(770)0");
  const std::array<std::string, 3> modes = {"auto", "positive", "negative"};
  std::array<std::vector<std::complex<double>>, 3> amplitude;
  std::size_t physical = 0;
  for (const auto i : indices(modes)) {
    CAPTURE(family, final, modes[i]);
    const auto tune = WriteModifiedPhotoVMTune("tu_sign_" + family + "_" + modes[i], [&](auto &card) {
      card.at("PARAM_REGGE").at("TU_SIGN").at(family) = modes[i];
    });
    ToyHelicityProcess process;
    ConfigureToyProductionProcess(process, family, "CON", final);
    process.SetTuneForTest(tune.first);
    REQUIRE_NOTHROW(process.InitializeProcessAmplitude());
    const auto &state = process.state.lts;
    physical = state.decaytree[0].p.spinX2 % 2 == 0 ? 1 : 2;
    const double sign = (i == 0 ? physical : i) == 1 ? 1.0 : -1.0;
    REQUIRE_FALSE(state.process.CONT_TU_SIGN.empty());
    for (const auto value : state.process.CONT_TU_SIGN) { CHECK(value == Approx(sign)); }

    auto lts = AsymmetricContinuumPairForTest(state.decaytree[0].p.pdg, state.decaytree[1].p.pdg);
    lts.process = state.process;
    lts.hamp.Configure(state.hamp.metadata);
    gra::MRegge regge(lts, process.state.model_tune,
                      gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::ContinuumTwoBody, family + "_tu_sign"));
    REQUIRE(TestReggeCon(regge, lts, gra::ParseReggeProductionModel(family)) > 0.0);
    amplitude[i].assign(lts.hamp.begin(), lts.hamp.end());
  }
  RequireVectorNear(amplitude[0], amplitude[physical], 1.0e-11);
  const double scale = std::max(gra::SquaredNorm(amplitude[1]), gra::SquaredNorm(amplitude[2]));
  // Check that both t and u diagrams survive the override
  for (const double sign : {-1.0, 1.0}) {
    CHECK(gra::SquaredNorm(LinearAmplitudeSectionForTest(amplitude[1], amplitude[2], 1.0, sign)) > 1.0e-12 * scale);
  }
}

TEST_CASE("GP setup expands the selected effective DL exchange bank", "[gra::MProcess][GP][physics]") {
  ToyHelicityProcess proc;
  ConfigureToyProductionProcess(proc, "GP", "CON", "pi+ pi-");
  REQUIRE_NOTHROW(proc.InitializeProcessAmplitude());

  const auto card     = nlohmann::json::parse(gra::aux::GetInputData(modelfile));
  const auto selected = card.at("PARAM_REGGE").at("PARAM_CON").at("GP").at("[211,-211]").get<std::vector<std::vector<int>>>();
  std::vector<std::vector<int>> expected;
  for (const auto &pair : selected) {
    REQUIRE(pair.size() == 2);
    expected.push_back(pair);
    if (pair[0] != pair[1]) { expected.push_back({pair[1], pair[0]}); }
  }
  REQUIRE(proc.state.lts.process.CONT_PRODUCTION == expected);
  REQUIRE(proc.state.lts.process.CONT_PRODUCTIONTREE.size() == expected.size());
  REQUIRE(proc.state.lts.process.CONTINUUM_GP.size() == expected.size());
  for (std::size_t i = 0; i < expected.size(); ++i) {
    REQUIRE(proc.state.lts.process.CONT_PRODUCTIONTREE[i].size() == 2);
    CHECK(proc.state.lts.process.CONT_PRODUCTIONTREE[i][0].p.pdg == expected[i][0]);
    CHECK(proc.state.lts.process.CONT_PRODUCTIONTREE[i][1].p.pdg == expected[i][1]);
    REQUIRE(proc.state.lts.process.CONTINUUM_GP[i].size() == 4);
  }
}

TEST_CASE("MP and XP initialize direct four and six pion ladders", "[gra::MProcess][continuum][multiregge]") {
  for (const std::string family : {"MP", "XP"}) {
    for (const std::string final_state : {"pi+ pi- pi+ pi-", "pi+ pi- pi+ pi- pi+ pi-"}) {
      CAPTURE(family, final_state);
      ToyHelicityProcess proc;
      ConfigureToyProductionProcess(proc, family, "CON", final_state);
      REQUIRE_NOTHROW(proc.InitializeProcessAmplitude());
      CHECK_FALSE(proc.state.lts.process.CONT_LADDER_PERMUTATIONS.empty());
      CHECK(proc.state.lts.hamp.metadata.spin_basis == gra::ScreeningSpinBasis::ProtonIdentity);
      CHECK(proc.state.lts.hamp.metadata.spin_rows == 1);
      CHECK(proc.state.lts.hamp.metadata.amplitude_normalization == Approx(1.0));
    }
  }
}

// Check that direct ladders reject unsupported external spin before sampling
TEST_CASE("MP XP and GP direct ladders require spin zero final states",
          "[gra::MProcess][continuum][multiregge][validation]") {
  for (const std::string family : {"MP", "XP", "GP"}) {
    for (const std::string final_state : {"p+ p- p+ p-", "p+ p- p+ p- p+ p-"}) {
      CAPTURE(family, final_state);
      ToyHelicityProcess proc;
      ConfigureToyProductionProcess(proc, family, "CON", final_state);
      REQUIRE_THROWS_WITH(proc.InitializeProcessAmplitude(),
                          Catch::Matchers::Contains("continuum ladders require spin-zero final states"));
    }
  }
}

// Reject finite m in direct GP ladders while retaining the supported m zero channel
TEST_CASE("GP direct ladders require m zero pair residues", "[gra::MProcess][GP][continuum][multiregge][validation]") {
  SECTION("finite m pion tune is rejected") {
    const auto tune = WriteModifiedContinuumTune("gp_ladder_finite_m", "GP", [](auto &card) {
      for (const std::string sector : {"same", "opposite"}) {
        card.at("990").at("[211,211]").at(sector).at("helicity") = {{0, 0, 0, 1.0, 0.0}, {0, 0, -1, 0.2, 0.0}};
      }
    });
    ToyHelicityProcess proc;
    ConfigureToyProductionProcess(proc, "GP", "CON", "pi+ pi- pi+ pi-");
    proc.SetTuneForTest(tune);
    REQUIRE_THROWS_WITH(proc.InitializeProcessAmplitude(),
                        Catch::Matchers::Contains("finite nonzero-m GP ladder transport is not defined"));
  }

  SECTION("explicit m zero kaon tune initializes four and six particle ladders") {
    const std::string tune = WriteModifiedContinuumTune("gp_ladder_m_zero", "GP", [](auto &card) {
      auto &pair = card.at("990").at("[321,321]");
      for (const std::string sector : {"same", "opposite"}) {
        auto &rows = pair.at(sector).at("helicity");
        rows.erase(std::remove_if(rows.begin(), rows.end(),
                                  [](const auto &row) { return row.at(2).template get<int>() != 0; }),
                   rows.end());
      }
    });
    for (const std::string final_state : {"K+ K- K+ K-", "K+ K- K+ K- K+ K-"}) {
      CAPTURE(final_state);
      ToyHelicityProcess proc;
      ConfigureToyProductionProcess(proc, "GP", "CON", final_state);
      proc.SetTuneForTest(tune);
      REQUIRE_NOTHROW(proc.InitializeProcessAmplitude());
      CHECK_FALSE(proc.state.lts.process.CONT_LADDER_PERMUTATIONS.empty());
    }
  }
}

// Reject invalid XP pole fields during initialization, including spin-disabled production
TEST_CASE("XP validates massive continuum poles before sampling", "[gra::MProcess][XP][continuum][params]") {
  for (const std::string decay : {"p+ p-", "rho(770)0 rho(770)0"}) {
    for (const bool spin : {false, true}) {
      for (const double mass : {0.0, -1.0, std::numeric_limits<double>::quiet_NaN()}) {
        CAPTURE(decay, spin, mass);
        ToyHelicityProcess process;
        ConfigureToyProductionProcess(process, "XP", "CON", decay);
        process.SetSPINGEN(spin);
        process.state.lts.decaytree[0].p.mass = mass;
        REQUIRE_THROWS_WITH(process.InitializeProcessAmplitude(),
                            Catch::Matchers::Contains("massive pole requires a positive finite mass"));
      }
    }
  }
}

TEST_CASE("Scalar hadronic continuum selects its physical proton spin basis", "[gra::MProcess][continuum][screening]") {
  for (const std::string family : {"MP", "XP"}) {
    CAPTURE(family);
    ToyHelicityProcess physical;
    ConfigureToyProductionProcess(physical, family, "CON", "K+ K-");
    REQUIRE_NOTHROW(physical.InitializeProcessAmplitude());
    CHECK(physical.state.lts.hamp.metadata.spin_basis == gra::ScreeningSpinBasis::ProtonHelicity);
    CHECK(physical.state.lts.hamp.metadata.spin_rows == 4);
    CHECK(physical.state.lts.hamp.metadata.amplitude_normalization == Approx(0.25));
    CHECK(physical.state.lts.hamp.metadata.spin_transition_count == 4);

    ToyHelicityProcess blind;
    ConfigureToyProductionProcess(blind, family, "CON", "K+ K-");
    blind.SetSPINGEN(false);
    REQUIRE_NOTHROW(blind.InitializeProcessAmplitude());
    CHECK(blind.state.lts.hamp.metadata.spin_basis == gra::ScreeningSpinBasis::ProtonIdentity);
    CHECK(blind.state.lts.hamp.metadata.spin_rows == 1);
    CHECK(blind.state.lts.hamp.metadata.amplitude_normalization == Approx(1.0));
    CHECK(blind.state.lts.hamp.metadata.spin_transition_count == 4);
  }
}

TEST_CASE("GP analytic subchannels require the full CP declaration", "[gra::MProcess][GP][validation]") {
  ModelParamRestoreGuard restore_model;

  // Configure one process and return its analytic Pomeron and charged-pion
  // particles
  auto make_process_inputs = [](ToyHelicityProcess &proc) {
    ConfigureToyProductionProcess(proc, "GP", "CON", "pi+ pi-");
    proc.state.lts.process.MMAX = 2;
    const auto trajectory       = proc.state.lts.PDG.FindByPDG(990);
    const auto pion             = proc.state.lts.PDG.FindByPDG(211);
    const auto antipion         = proc.state.lts.PDG.FindByPDG(-211);
    return std::make_tuple(trajectory, pion, antipion);
  };

  SECTION("C declaration is required") {
    const std::string tune = WriteModifiedContinuumTune(
        "gp_c_declaration", "GP", [](auto &j) { j.at("990").at("[211,211]").at("same").at("CP").at(0) = false; });
    ToyHelicityProcess proc;
    make_process_inputs(proc);
    gra::MODELPARAM = tune;
    REQUIRE_THROWS_AS(proc.SetTuneForTest(tune), std::invalid_argument);
  }

  SECTION("parity declaration is required") {
    const std::string tune = WriteModifiedContinuumTune(
        "gp_parity_declaration", "GP", [](auto &j) { j.at("990").at("[211,211]").at("same").at("CP").at(1) = false; });
    ToyHelicityProcess proc;
    make_process_inputs(proc);
    gra::MODELPARAM = tune;
    REQUIRE_THROWS_AS(proc.SetTuneForTest(tune), std::invalid_argument);
  }

  SECTION("all pair sectors validate their CP array") {
    const std::string  tune = WriteModifiedContinuumTune("gp_unselected_cp_shape", "GP", [](auto &j) {
      j.at("990").at("[211,211]").at("opposite").at("CP") = nlohmann::json::array({true});
    });
    ToyHelicityProcess proc;
    make_process_inputs(proc);
    gra::MODELPARAM = tune;
    REQUIRE_THROWS_AS(proc.SetTuneForTest(tune), std::invalid_argument);
  }

  SECTION("self cannot coexist with a charge sector") {
    const std::string  tune = WriteModifiedContinuumTune("gp_ambiguous_self_sector", "GP", [](auto &j) {
      auto &pair   = j.at("990").at("[111,111]");
      pair["same"] = pair.at("self");
    });
    ToyHelicityProcess proc;
    make_process_inputs(proc);
    gra::MODELPARAM = tune;
    REQUIRE_THROWS_AS(proc.SetTuneForTest(tune), std::invalid_argument);
  }

  SECTION("nonphoton rows retain sparse finite m labels") {
    const std::string  tune = WriteModifiedContinuumTune("gp_raw_continuum", "GP", [](auto &j) {
      j.at("990").at("[211,211]").at("opposite") = {
          {"basis", "crossed_helicity"},
          {"CP", {true, true}},
          {"helicity", {{0, 0, 2, 0.7, 0.2}, {0, 0, 1, 0.4, -0.3}, {0, 0, 0, 0.9, 0.5}}}};
    });
    ToyHelicityProcess proc;
    auto [trajectory, pion, antipion] = make_process_inputs(proc);
    proc.SetMMAX(3);
    gra::MODELPARAM = tune;
    proc.SetTuneForTest(tune);
    const auto hel = proc.ProcessHelicityStructure(trajectory, {pion, antipion}, true, true, "", false,
                                                   gra::spin::VertexContext::SubTUChannelExchange);
    REQUIRE(hel.UsesReggeDomain());
    REQUIRE(hel.UsesHelicityCouplings());
    REQUIRE(hel.exchange_basis == gra::ExchangeBasisType::ReggeHelicity);
    REQUIRE(hel.analytic_MMAX == 3);
    REQUIRE(hel.T.size_col() == 7);
    std::size_t set_rows = 0;
    for (std::size_t row = 0; row < hel.T_set.size_row(); ++row) {
      for (std::size_t column = 0; column < hel.T_set.size_col(); ++column) {
        set_rows += hel.T_set[row][column] ? 1U : 0U;
      }
    }
    CHECK(set_rows == 5);
    const std::array<std::pair<int, std::complex<double>>, 5> expected = {{
        {-2, std::polar(0.7, 0.2)},
        {-1, std::polar(0.4, -0.3)},
        {0, std::polar(0.9, 0.5)},
        {1, std::polar(0.4, -0.3)},
        {2, std::polar(0.7, 0.2)},
    }};
    for (const auto &[m, value] : expected) {
      const std::size_t column = gra::gpom::AnalyticMIndex(m, hel.analytic_MMAX, "parsed finite-m crossed residue");
      REQUIRE(hel.T_set[0][column]);
      RequireComplexNear(hel.T[0][column], value, 1.0e-12);
    }
    for (const int m : {-3, 3}) {
      const std::size_t column = gra::gpom::AnalyticMIndex(m, hel.analytic_MMAX, "parsed finite-m sparse padding");
      CHECK_FALSE(hel.T_set[0][column]);
    }

    const gra::M4Vec first_axis(0.6, -0.4, 0.8, 1.3);
    const gra::M4Vec second_axis(0.6, -0.4, -2.1, 2.7);
    const double     alpha  = 1.37;
    const auto       first  = gra::gpom::Crossed(hel, alpha, first_axis, false);
    const auto       second = gra::gpom::Crossed(hel, alpha, second_axis, false);
    RequireMatrixNear(first.frame, second.frame, 2.0e-13);
    for (const auto &[m, value] : expected) {
      const std::size_t column        = gra::gpom::AnalyticMIndex(m, hel.analytic_MMAX, "finite-m crossed section");
      const auto        upper_section = std::exp(-gra::math::zi * static_cast<double>(m) * first_axis.Phi());
      const double spherical = m == 1 ? -1.0 : 1.0;
      RequireComplexNear(first.frame[0][column], spherical * value * upper_section, 2.0e-13);
    }

    const auto lower               = gra::gpom::Crossed(hel, alpha, first_axis, true);
    const auto lower_polar_changed = gra::gpom::Crossed(hel, alpha, second_axis, true);
    RequireMatrixNear(lower.frame, lower_polar_changed.frame, 2.0e-13);
    for (const auto &[m, value] : expected) {
      const std::size_t column = gra::gpom::AnalyticMIndex(m, hel.analytic_MMAX, "finite-m lower crossed section");
      const auto        lower_section = std::exp(gra::math::zi * static_cast<double>(m) * first_axis.Phi());
      const double spherical = m == -1 ? -1.0 : 1.0;
      RequireComplexNear(lower.frame[0][column], spherical * value * lower_section, 2.0e-13);
    }
  }

  SECTION("crossed LS keeps independent sparse m coefficient blocks") {
    const std::string  tune = WriteModifiedContinuumTune("gp_finite_m_ls", "GP", [](auto &j) {
      j.at("990").at("[211,211]").at("opposite") = {
          {"basis", "crossed_ls"}, {"CP", {true, true}}, {"g_ls", {{2, 0, 2, 0.6, 0.3}, {2, 0, 0, 0.85, -0.4}}}};
    });
    ToyHelicityProcess proc;
    auto [trajectory, pion, antipion] = make_process_inputs(proc);
    proc.SetMMAX(3);
    gra::MODELPARAM = tune;
    proc.SetTuneForTest(tune);
    auto hel = proc.ProcessHelicityStructure(trajectory, {pion, antipion}, true, true, "", false,
                                             gra::spin::VertexContext::SubTUChannelExchange);
    REQUIRE(hel.UsesReggeDomain());
    REQUIRE(hel.UsesLSCouplings());
    REQUIRE_FALSE(hel.UsesHelicityCouplings());
    REQUIRE(hel.exchange_basis == gra::ExchangeBasisType::ReggeHelicity);
    REQUIRE(hel.analytic_MMAX == 3);
    REQUIRE(hel.m_ls.size() == 7);

    const auto g2 = std::polar(0.6, 0.3);
    const auto g0 = std::polar(0.85, -0.4);
    for (const int m : {-2, 2}) {
      const std::size_t column = gra::gpom::AnalyticMIndex(m, hel.analytic_MMAX, "parsed finite-m crossed LS");
      REQUIRE(hel.m_ls[column].Size() == 1);
      RequireComplexNear(hel.m_ls[column].At(2, 0), g2, 1.0e-12);
    }
    const std::size_t zero = gra::gpom::AnalyticMIndex(0, hel.analytic_MMAX, "parsed finite-m crossed LS zero");
    REQUIRE(hel.m_ls[zero].Size() == 1);
    RequireComplexNear(hel.m_ls[zero].At(2, 0), g0, 1.0e-12);
    for (const int m : {-3, -1, 1, 3}) {
      const std::size_t column = gra::gpom::AnalyticMIndex(m, hel.analytic_MMAX, "parsed finite-m crossed LS padding");
      CHECK(hel.m_ls[column].Empty());
    }

    gra::gpom::InitCrossedLS(hel, 2);
    const double     alpha = 1.37;
    const gra::M4Vec first_axis(0.5, 0.3, 0.9, 1.4);
    const gra::M4Vec second_axis(0.5, 0.3, -1.8, 2.3);
    const auto       finite        = gra::gpom::Crossed(hel, alpha, first_axis, false);
    const auto       polar_changed = gra::gpom::Crossed(hel, alpha, second_axis, false);
    RequireMatrixNear(finite.frame, polar_changed.frame, 2.0e-13);
    const auto lower               = gra::gpom::Crossed(hel, alpha, first_axis, true);
    const auto lower_polar_changed = gra::gpom::Crossed(hel, alpha, second_axis, true);
    RequireMatrixNear(lower.frame, lower_polar_changed.frame, 2.0e-13);
    const std::size_t positive =
        gra::gpom::AnalyticMIndex(2, hel.analytic_MMAX, "evaluated finite-m crossed LS positive");
    const std::size_t negative =
        gra::gpom::AnalyticMIndex(-2, hel.analytic_MMAX, "evaluated finite-m crossed LS negative");
    RequireComplexNear(finite.helicity.T[0][negative], finite.helicity.T[0][positive], 2.0e-12);

    gra::HELMatrix reduced;
    gra::gpom::InitCrossed(reduced, 0.0, 0.0, 0, "reduced crossed LS m=0 reference");
    reduced.C_symmetry     = true;
    reduced.P_symmetry     = true;
    reduced.coupling_basis = gra::CouplingBasis::LS;
    reduced.m_ls.resize(1);
    reduced.m_ls[0].Set(2, 0, g0);
    gra::gpom::InitCrossedLS(reduced, 2);
    const auto reference = gra::gpom::Crossed(reduced, alpha, gra::M4Vec{}, false);
    RequireComplexNear(finite.helicity.T[0][zero], reference.helicity.T[0][0], 2.0e-12);
    RequireComplexNear(finite.frame[0][zero], reference.frame[0][0], 2.0e-12);
  }

  SECTION("repeated symmetry representatives are rejected") {
    const std::string  tune = WriteModifiedContinuumTune("gp_raw_continuum_parity", "GP", [](auto &j) {
      j.at("990").at("[211,211]").at("opposite") = {{"basis", "crossed_helicity"},
                                                    {"CP", {true, true}},
                                                    {"helicity", {{0, 0, 2, 0.7, 0.2}, {0, 0, -2, 0.7, 0.2}}}};
    });
    ToyHelicityProcess proc;
    auto [trajectory, pion, antipion] = make_process_inputs(proc);
    proc.SetMMAX(2);
    gra::MODELPARAM = tune;
    proc.SetTuneForTest(tune);
    REQUIRE_THROWS(proc.ProcessHelicityStructure(trajectory, {pion, antipion}, true, true, "", false,
                                                 gra::spin::VertexContext::SubTUChannelExchange));
  }
}

TEST_CASE("GP continuum rows require exact ordered legs", "[gra::MProcess][GP][helicity][validation]") {
  ModelParamRestoreGuard restore_model;
  const std::string      tune = WriteModifiedContinuumTune("gp_ordered_analytic_rows", "GP", [](auto &j) {
    auto pair = j.at("990").at("[211,211]");
    j.at("990").erase("[211,211]");
    j.at("990")["[321,211]"] = std::move(pair);
  });

  ToyHelicityProcess proc;
  ConfigureToyProductionProcess(proc, "GP", "CON", "pi+ pi-");
  gra::MODELPARAM = tune;
  REQUIRE_THROWS_AS(proc.SetTuneForTest(tune), std::invalid_argument);
}

TEST_CASE("GP setup rejects duplicate order-equivalent final-state rows", "[gra::MProcess][GP][validation]") {
  ModelParamRestoreGuard restore_model;
  const auto             tune = WriteModifiedPhotoVMTune("gp_duplicate_runtime_fallback", [](auto &j) {
    auto &gp         = j.at("PARAM_REGGE").at("PARAM_CON").at("GP");
    gp["[-211,211]"] = gp.at("[211,-211]");
  });

  ToyHelicityProcess proc;
  ConfigureToyProductionProcess(proc, "GP", "CON", "pi+ pi-");
  proc.SetModelTune(gra::MModelTune::Load(tune.second));
  gra::MODELPARAM = tune.first;
  REQUIRE_THROWS(proc.InitializeProcessAmplitude());
}

// Check complete fixed-spin secondary-Reggeon data through process setup
TEST_CASE("MP and XP initialize every fixed secondary-Reggeon representative", "[gra::MProcess][MP][XP][Reggeon]") {
  ModelParamRestoreGuard restore_model;
  gra::MODELPARAM = "TUNE0";
  MRandom rng;

  for (const std::string model : {"MP", "XP"}) {
    for (const int exchange : {9925, 9943}) {
      CAPTURE(model, exchange);
      ToyHelicityProcess proc;
      ConfigureToyProductionProcess(proc, model, "RES", "pi+ pi-");
      auto f0 = gra::resonance::Read("RES/f0_500.json", rng, gra::ParseReggeProductionModel(model));
      if (model == "MP") {
        SetMPFusion(f0);
        f0.MP.channels.front().exchange = {exchange, exchange};
      } else {
        auto &channel    = f0.XP.channels.front();
        channel.exchange = {exchange, exchange};
        const auto pole  = proc.state.lts.PDG.FindByPDG(exchange);
        const auto operators =
            gra::spin::CanonicalPoleOperators(f0.p, pole, pole, true, channel.C_symmetry, channel.P_symmetry);
        REQUIRE_FALSE(operators.empty());
        const auto ls_reference = operators.front().coupling;
        channel.g_ls.Clear();
        for (const auto &operator_term : operators) {
          const bool selected =
              operator_term.coupling.l == ls_reference.l && operator_term.coupling.two_s == ls_reference.two_s;
          channel.g_ls.Set(operator_term.coupling.l, operator_term.coupling.two_s, selected ? 1.0 : 0.0);
        }
      }
      proc.SetResonances({{"f0_500", f0}});

      REQUIRE_NOTHROW(proc.InitializeProcessAmplitude());
      const auto configured = proc.GetResonances().at("f0_500");
      REQUIRE(configured.production.size() == 1);
      REQUIRE(configured.production.front().tree.size() == 2);
      CHECK(configured.production.front().tree[0].p.pdg == exchange);
      CHECK(configured.production.front().tree[1].p.pdg == exchange);
      CHECK(configured.production.front().tree[0].hel.lambda_values.size_row() == 4);
      CHECK(configured.production.front().tree[1].hel.lambda_values.size_row() == 4);
      CHECK(configured.production.front().tree[0].hel.Jz_values.size() ==
            static_cast<std::size_t>(configured.production.front().tree[0].p.spinX2 + 1));
      CHECK(configured.production.front().tree[1].hel.Jz_values.size() ==
            static_cast<std::size_t>(configured.production.front().tree[1].p.spinX2 + 1));
      if (model == "XP") { REQUIRE(configured.production.size() == 1); }
    }
  }
}

// Check production symmetries and physical diphoton coupling normalization
TEST_CASE("MP production C checks cover fixed-spin PP RP and OP exchange IDs", "[gra::spin][MP]") {
  gra::MODELPARAM = "TUNE0";
  MRandom rng;
  auto    f0  = gra::resonance::Read("RES/f0_500.json", rng, gra::ReggeProductionModel::MP);
  auto    rho = gra::resonance::Read("RES/rho_770.json", rng, gra::ReggeProductionModel::MP);

  SECTION("C-even f0 accepts fixed Pomeron-Pomeron") {
    ToyHelicityProcess proc;
    ConfigureToyProductionProcess(proc, "MP", "RES", "pi+ pi-");
    f0.MP.channels.front().exchange = {991, 991};
    proc.SetResonances({{"f0_500", f0}});
    REQUIRE_NOTHROW(proc.InitializeProcessAmplitude());

    const auto configured = proc.GetResonances().at("f0_500");
    REQUIRE(configured.production.size() == 1);
    REQUIRE(configured.production.size() == 1);
    CHECK(configured.production[0].tree[0].p.pdg == 991);
    CHECK(configured.production[0].tree[1].p.pdg == 991);
  }

  SECTION("C-even f0 accepts fixed f2-Reggeon-Pomeron") {
    ToyHelicityProcess proc;
    ConfigureToyProductionProcess(proc, "MP", "RES", "pi+ pi-");
    f0.MP.channels.front().exchange = {9915, 991};
    proc.SetResonances({{"f0_500", f0}});
    REQUIRE_NOTHROW(proc.InitializeProcessAmplitude());
  }

  SECTION("single excitation accepts mirrored Pomeron-Reggeon orderings") {
    ToyHelicityProcess proc;
    ConfigureToyProductionProcess(proc, "MP", "RES", "pi+ pi-");
    proc.SetExcitation(1);
    f0.MP.channels.front().exchange = {9915, 991};
    proc.SetResonances({{"f0_500", f0}});
    REQUIRE_NOTHROW(proc.InitializeProcessAmplitude());
  }

  SECTION("double excitation rejects Pomeron-Reggeon orderings") {
    ToyHelicityProcess proc;
    ConfigureToyProductionProcess(proc, "MP", "RES", "pi+ pi-");
    proc.SetExcitation(2);
    f0.MP.channels.front().exchange = {9915, 991};
    proc.SetResonances({{"f0_500", f0}});
    REQUIRE_THROWS(proc.InitializeProcessAmplitude());
  }

  SECTION("single excitation rejects Reggeon-Reggeon orderings") {
    ToyHelicityProcess proc;
    ConfigureToyProductionProcess(proc, "MP", "RES", "pi+ pi-");
    proc.SetExcitation(1);
    f0.MP.channels.front().exchange = {9915, 9915};
    proc.SetResonances({{"f0_500", f0}});
    REQUIRE_THROWS(proc.InitializeProcessAmplitude());
  }

  SECTION("double excitation accepts Pomeron-Pomeron ordering") {
    ToyHelicityProcess proc;
    ConfigureToyProductionProcess(proc, "MP", "RES", "pi+ pi-");
    proc.SetExcitation(2);
    f0.MP.channels.front().exchange = {991, 991};
    proc.SetResonances({{"f0_500", f0}});
    REQUIRE_NOTHROW(proc.InitializeProcessAmplitude());
  }

  SECTION(
      "C-odd rho expands fixed vector-Odderon-Pomeron into both spin "
      "orderings") {
    ToyHelicityProcess proc;
    ConfigureToyProductionProcess(proc, "MP", "RES", "pi+ pi-");
    SetMPFusion(rho);
    auto &channel    = rho.MP.channels.front();
    channel.exchange = {9993, 991};
    channel.basis    = gra::ReggeVertexBasis::AutoMinL;
    channel.g        = std::polar(0.7, 0.2);
    proc.SetResonances({{"rho_770", rho}});
    REQUIRE_NOTHROW(proc.InitializeProcessAmplitude());

    const auto configured = proc.GetResonances().at("rho_770");
    REQUIRE(configured.production.size() == 2);
    REQUIRE(configured.production.size() == 2);
    CHECK(configured.production[0].tree[0].p.pdg == 9993);
    CHECK(configured.production[0].tree[1].p.pdg == 991);
    CHECK(configured.production[1].tree[0].p.pdg == 991);
    CHECK(configured.production[1].tree[1].p.pdg == 9993);
    REQUIRE(configured.production.size() == 2);
    RequireCanonicalPoleVertex(configured.production[0].pole.value());
    RequireCanonicalPoleVertex(configured.production[1].pole.value());
  }

  SECTION("MP gamma-gamma coupling matches the physical transverse tensor") {
    ToyHelicityProcess proc;
    ConfigureToyProductionProcess(proc, "MP", "RES", "pi+ pi-");
    auto f2 = gra::resonance::Read("RES/f2_1270_yy.json", rng, gra::ReggeProductionModel::MP);
    proc.SetResonances({{"f2_1270_yy", f2}});
    REQUIRE_NOTHROW(proc.InitializeProcessAmplitude());

    const gra::PARAM_RES configured = proc.GetResonances().at("f2_1270_yy");
    REQUIRE(configured.production.size() == 1);
    REQUIRE(configured.production.size() == 1);
    CHECK(configured.production[0].tree[0].p.pdg == 22);
    CHECK(configured.production[0].tree[1].p.pdg == 22);
    REQUIRE(configured.production.size() == 1);
    const double norm2 =
        gra::spin::PoleLSReduced(configured.production.front().pole.value(), 0.5 * configured.p.mass).FrobNorm2();
    CHECK(norm2 == Approx(gra::math::pow2(gra::resonance::GammaGammaResonanceCoupling(configured.p))));
  }

  SECTION("GP gamma-gamma LS and helicity have the same width-normalized amplitude") {
    ModelParamRestoreGuard restore;
    for (const bool derivative : {false, true}) {
      const auto files = WriteModifiedPhotoVMTune("gp_diphoton_" + std::to_string(derivative), [&](auto &card) {
        card.at("PARAM_REGGE").at("DERIVATIVE_FACTOR").at("GP") = derivative;
      });
      const auto tune = gra::MModelTune::Load(files.second);
      for (const double scale : {0.7, 1.4}) {
        CAPTURE(derivative, scale);
        std::vector<gra::HelAmp> amplitudes;
        for (const auto basis : {gra::ReggeVertexBasis::LS, gra::ReggeVertexBasis::Helicity}) {
          ToyHelicityProcess proc;
          ConfigureToyProductionProcess(proc, "GP", "RES", "pi+ pi-");
          proc.SetModelTune(tune);
          proc.SetHelicityConfig(tune);
          auto f2 = gra::resonance::Read("RES/f2_1270_yy.json", rng, gra::ReggeProductionModel::GP);
          auto &channel = f2.GP.channels.front();
          if (basis == gra::ReggeVertexBasis::LS) {
            channel.basis = basis;
            channel.Lambda = scale;
            channel.helicity.clear();
            channel.g_helicity.clear();
            // The raw STF tensors cancel H0 for g_22 = 6 g_20
            channel.g_ls.Set(0, 4, 0.0);
            channel.g_ls.Set(2, 0, 1.0);
            channel.g_ls.Set(2, 4, 6.0);
            channel.g_ls.Set(4, 4, 0.0);
          }
          proc.SetResonances({{"f2_1270_yy", f2}});
          REQUIRE_NOTHROW(proc.InitializeProcessAmplitude());
          const auto configured = proc.GetResonances().at("f2_1270_yy");
          REQUIRE(configured.production.size() == 1);
          if (basis == gra::ReggeVertexBasis::Helicity) {
            const auto &hel = configured.production.front().hel;
            CHECK(hel.T.MaskedSquaredNorm(hel.T_set) == Approx(gra::math::pow2(
                gra::resonance::GammaGammaResonanceCoupling(configured.p))));
          }
          auto lts = MakeToyCoherentPhotonLTS();
          lts.process = proc.state.lts.process;
          const double q = 0.5 * configured.p.mass;
          lts.q1_in_X = gra::M4Vec(0.0, 0.0, q, q);
          lts.q2_in_X = gra::M4Vec(0.0, 0.0, -q, q);
          gra::MRegge regge(lts, tune, gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "diphoton width"));
          const auto param = ReggeParametersForTest(regge, lts);
          const auto production = gra::gpom::Resonance(lts, *param, configured, nullptr);
          REQUIRE(production.size() == 1);
          REQUIRE(production.front().FrobNorm2() > 0.0);
          amplitudes.push_back(production.front());
        }
        CHECK(amplitudes[0].FrobNorm2() == Approx(amplitudes[1].FrobNorm2()).epsilon(2.0e-12));
        RequireMatrixNear(amplitudes[0], amplitudes[1], 2.0e-12);
      }
    }
  }

  SECTION("XP gamma-gamma helicity input retains the diphoton width") {
    ToyHelicityProcess proc;
    ConfigureToyProductionProcess(proc, "XP", "RES", "pi+ pi-");
    auto f2 = gra::resonance::Read("RES/f2_1270_yy.json", rng, gra::ReggeProductionModel::XP);
    proc.SetResonances({{"f2_1270_yy", f2}});
    REQUIRE_NOTHROW(proc.InitializeProcessAmplitude());

    const gra::PARAM_RES configured = proc.GetResonances().at("f2_1270_yy");
    REQUIRE(configured.production.size() == 1);
    const auto &vertex = configured.production.front().pole.value();
    REQUIRE_FALSE(vertex.terms.empty());
    CHECK(vertex.leg1_transverse_pole);
    CHECK(vertex.leg2_transverse_pole);
    const double norm2 = gra::math::pow2(gra::resonance::GammaGammaResonanceCoupling(configured.p));
    const auto pole = gra::spin::PoleLSReduced(vertex, 0.5 * configured.p.mass);
    CHECK(pole.FrobNorm2() == Approx(norm2));
    RequireComplexNear(pole[0][2], pole[2][0], 1.0e-12);
    RequireComplexNear(pole[0][0], pole[2][2], 1.0e-12);
    // Direct helicities specify the LS tensor at q/Lambda = 1
    const auto reference = gra::spin::PoleLSReduced(vertex, vertex.Lambda);
    REQUIRE(reference.size_row() == 3);
    REQUIRE(reference.size_col() == 3);
    CHECK(reference.FrobNorm2() == Approx(2.0 * std::norm(reference[0][2])).epsilon(1.0e-12));
    CHECK(std::abs(reference[0][0]) == Approx(0.0).margin(1.0e-12));
    CHECK(std::abs(reference[2][2]) == Approx(0.0).margin(1.0e-12));
  }

  SECTION("XP rejects a central vertex with no active operator") {
    ToyHelicityProcess proc;
    ConfigureToyProductionProcess(proc, "XP", "RES", "pi+ pi-");
    auto f2                              = gra::resonance::Read("RES/f2_1270_yy.json", rng, gra::ReggeProductionModel::XP);
    f2.XP.channels.front().width_derived = false;
    for (auto &coupling : f2.XP.channels.front().g_helicity) { coupling = 0.0; }
    proc.SetResonances({{"f2_1270_yy", f2}});
    REQUIRE_THROWS(proc.InitializeProcessAmplitude());
  }

  SECTION("XP retains absolute operators in spin-disabled production") {
    ToyHelicityProcess proc;
    ConfigureToyProductionProcess(proc, "XP", "RES", "pi+ pi-");
    proc.SetSPINGEN(false);
    const auto f0 = gra::resonance::Read("RES/f0_980.json", rng, gra::ReggeProductionModel::XP);
    proc.SetResonances({{"f0_980", f0}});
    REQUIRE_NOTHROW(proc.InitializeProcessAmplitude());
    const auto &configured = proc.GetResonances().at("f0_980");
    REQUIRE(configured.production.size() == configured.production.size());
  }

  SECTION(
      "C-odd rho photoproduction uses automatic MP central "
      "tensor before spin steering") {
    ToyHelicityProcess proc;
    ConfigureToyProductionProcess(proc, "MP", "RES", "pi+ pi-");
    SetMPFusion(rho);
    auto &channel    = rho.MP.channels.front();
    channel.exchange = {22, 991};
    channel.basis    = gra::ReggeVertexBasis::AutoEqualLS;
    channel.g        = 1.0;
    proc.SetResonances({{"rho_770", rho}});
    REQUIRE_NOTHROW(proc.InitializeProcessAmplitude());

    const auto resonances = proc.GetResonances();
    REQUIRE(resonances.count("rho_770") == 1);
    const auto &configured = resonances.at("rho_770");
    REQUIRE(configured.production.size() == 2);
    const auto &hel = configured.production;
    REQUIRE(hel.size() == 2);
    CHECK(configured.production[0].tree[0].p.pdg == 22);
    CHECK(configured.production[0].tree[1].p.pdg == 991);
    CHECK(configured.production[1].tree[0].p.pdg == 991);
    CHECK(configured.production[1].tree[1].p.pdg == 22);

    REQUIRE(configured.production.size() == 2);
    RequireCanonicalPoleVertex(configured.production[0].pole.value());
    RequireCanonicalPoleVertex(configured.production[1].pole.value());
    REQUIRE(configured.production[0].pole.value().terms.size() == 2);
    REQUIRE(configured.production[1].pole.value().terms.size() == 2);
  }

  SECTION("C-even f0 rejects fixed Odderon-Pomeron") {
    ToyHelicityProcess proc;
    ConfigureToyProductionProcess(proc, "MP", "RES", "pi+ pi-");
    f0.MP.channels.front().exchange = {9991, 991};
    proc.SetResonances({{"f0_500", f0}});
    REQUIRE_THROWS(proc.InitializeProcessAmplitude());
  }

  SECTION("C-odd rho rejects fixed Pomeron-Pomeron") {
    ToyHelicityProcess proc;
    ConfigureToyProductionProcess(proc, "MP", "RES", "pi+ pi-");
    rho.MP.channels.front().exchange = {991, 991};
    proc.SetResonances({{"rho_770", rho}});
    REQUIRE_THROWS(proc.InitializeProcessAmplitude());
  }
}

// Read zero model couplings without discarding their quantum numbers or spin steering
TEST_CASE("Zero resonance model cards retain their channel definitions", "[gra::MProcess][resonance][zero]") {
  auto card = nlohmann::json::parse(gra::aux::GetInputData(gra::ResolveModelDataFile("TUNE0", "RES/f2_1270.json")));
  for (const std::string model : {"MP", "XP", "GP"}) {
    for (auto &[key, block] : card.at("PARAM_RES").at("MODELS").at(model).items()) {
      if (!key.starts_with("[")) { continue; }
      if (block.contains("g")) { block["g"] = {0.0, 0.0}; }
      for (const std::string field : {"g_ls", "helicity"}) {
        if (!block.contains(field)) { continue; }
        for (auto &row : block[field]) { row[2] = 0.0; row[3] = 0.0; }
      }
    }
  }
  const std::filesystem::path directory = "tmp/graniitti_zero_resonance_models";
  std::filesystem::create_directories(directory / "RES");
  std::ofstream output(directory / "RES/f2_1270.json");
  REQUIRE(output.good());
  output << card.dump(2);
  output.close();
  gra::MRandom rng;
  const auto zero = gra::resonance::Read("RES/f2_1270.json", rng, gra::ReggeProductionModel::GP, directory.string());
  for (const auto *channels : {&zero.MP.channels, &zero.XP.channels, &zero.GP.channels}) {
    REQUIRE_FALSE(channels->empty());
    for (const auto &channel : *channels) { CHECK_FALSE(channel.Active(0.0)); }
  }
}

// Validate zero production blocks before omitting them from the resonance plan
TEST_CASE("Zero resonance production is omitted after initialization", "[gra::MProcess][resonance][zero]") {
  for (const std::string model : {"MP", "XP", "GP"}) {
    for (const std::string name : {"f2_1270", "rho_770"}) {
      CAPTURE(model, name);
      ToyHelicityProcess proc;
      ConfigureToyProductionProcess(proc, model, "RES", "pi+ pi-");
      auto zero = gra::resonance::Read("RES/" + name + ".json", proc.state.random, gra::ParseReggeProductionModel(model));
      const auto active = gra::resonance::Read("RES/f0_980.json", proc.state.random, gra::ParseReggeProductionModel(model));
      auto &channels = model == "MP" ? zero.MP.channels : model == "XP" ? zero.XP.channels : zero.GP.channels;
      for (auto &channel : channels) {
        channel.g = 0.0;
        for (auto &term : channel.g_ls) { term.coefficient = 0.0; }
        for (auto &coupling : channel.g_helicity) { coupling = 0.0; }
      }
      proc.SetResonances({{"f0_980", active}, {name, zero}});
      REQUIRE_NOTHROW(proc.InitializeProcessAmplitude());
      REQUIRE(proc.GetResonances().size() == 1);
      CHECK(proc.GetResonances().contains("f0_980"));
      CHECK_FALSE(proc.GetResonances().at("f0_980").production.empty());
    }
  }
}

TEST_CASE("MP spin steering preserves the physical state normalization", "[gra::spin][resonance][MP][normalization]") {
  MMatrix<std::complex<double>> central(4, 3, 0.0);
  central[0][0] = {0.4, 0.1};
  central[0][2] = {-0.2, 0.3};
  central[1][1] = {0.7, -0.2};
  central[2][0] = {-0.1, 0.5};
  central[2][1] = {0.2, 0.4};
  central[3][2] = {0.6, -0.3};

  const std::vector<std::complex<double>> pure_state = {{0.3, 0.1}, {-0.2, 0.4}, {0.1, -0.15}};

  const auto combined = gra::spin::Steer(central, pure_state);
  CHECK(combined.FrobNorm2() == Approx(central.RightDiagonalFrobNorm2(pure_state)).epsilon(1e-12));

  const std::vector<std::complex<double>> unsupported = {0.0, 0.0, 1.0};
  MMatrix<std::complex<double>>           no_third_column(2, 3, 0.0);
  no_third_column[0][0] = 1.0;
  no_third_column[1][1] = 1.0;
  CHECK(gra::spin::Steer(no_third_column, unsupported).FrobNorm2() == Approx(0.0).margin(1e-24));
  REQUIRE_THROWS_AS(gra::spin::Steer(no_third_column, {1.0}), std::invalid_argument);
  no_third_column[0][2] = std::numeric_limits<double>::infinity();
  REQUIRE_FALSE(gra::spin::Steer(no_third_column, unsupported).IsFinite());
}

// Retain angular momentum zeros in the dynamical fusion prescription
TEST_CASE("MP fusion retains its physical production zeros", "[gra::MRegge][MP][spin][regression][spin_zero]") {
  auto lts = ScalarPolePhasePointForTest(0.0, 0.91, -0.38, 1.3);
  lts.process.MP_FRAME = "CM";
  auto res = MakeToyScalarMPResonance();
  res.p.spinX2 = 4;
  PrepareToyPoleOperators(res, gra::ReggeProductionModel::MP, {{{2, 0, 1.0}}});
  const auto central = gra::mpom::Fusion(lts, *res.production.front().pole);
  CHECK(central.RightDiagonalFrobNorm2(gra::HelVec{0.0, 0.0, 1.0, 0.0, 0.0}) > 0.0);
  CHECK(central.RightDiagonalFrobNorm2(gra::HelVec{1.0, 0.0, 0.0, 0.0, 0.0}) == Approx(0.0).margin(1e-24));
}

TEST_CASE("XP production constraints cover fixed-spin PP RP and OP exchange IDs", "[gra::spin][XP]") {
  gra::MODELPARAM        = "TUNE0";
  const auto &pdg        = LoadedPDGTable();
  const auto  f0         = pdg.FindByPDG(9000221);
  const auto  rho0       = pdg.FindByPDG(113);
  const auto  pomeron0   = pdg.FindByPDG(991);
  const auto  f2_reggeon = pdg.FindByPDG(9915);
  MRandom     rng;

  SECTION("C-even f0 accepts configured fixed Pomeron-Pomeron rows") {
    ToyHelicityProcess proc;
    ConfigureToyProductionProcess(proc, "XP", "RES", "pi+ pi-");
    const auto f0_card = gra::resonance::Read("RES/f0_500.json", rng, gra::ReggeProductionModel::XP);
    proc.SetResonances({{"f0_500", f0_card}});
    REQUIRE_NOTHROW(proc.InitializeProcessAmplitude());
    const auto resonances = proc.GetResonances();
    REQUIRE(resonances.count("f0_500") == 1);
    REQUIRE(resonances.at("f0_500").production.size() == 1);
  }

  SECTION("C-odd rho accepts configured fixed Odderon-Pomeron rows") {
    ToyHelicityProcess proc;
    ConfigureToyProductionProcess(proc, "XP", "RES", "pi+ pi-");
    const auto rho_card = gra::resonance::Read("RES/rho_770_odd.json", rng, gra::ReggeProductionModel::XP);
    proc.SetResonances({{"rho_770_odd", rho_card}});
    REQUIRE_NOTHROW(proc.InitializeProcessAmplitude());
    const auto resonances = proc.GetResonances();
    REQUIRE(resonances.count("rho_770_odd") == 1);
    REQUIRE(resonances.at("rho_770_odd").production.size() == 2);
    const auto &tree = resonances.at("rho_770_odd").production;
    REQUIRE(tree.size() == 2);
    CHECK(tree[0].tree[1].p.pdg == 9993);
    CHECK(tree[1].tree[0].p.pdg == 9993);
  }

  SECTION("fixed f2-Reggeon-Pomeron direct H obeys the generic JW subspace") {
    auto rp_alpha = SingleLSDefinition(2, 4);
    REQUIRE_NOTHROW(gra::spin::InitTMatrix(rp_alpha, f0, f2_reggeon, pomeron0, true, "CON_MP.json", false, false));
    REQUIRE(rp_alpha.T.FrobNorm2() > 1e-12);

    auto direct = DirectDefinitionFromComputedT(rp_alpha);
    REQUIRE_NOTHROW(
        gra::spin::ValidateDirectTMatrix(direct, f0, f2_reggeon, pomeron0, true, "CON_MP.json", true, false));

    auto c_bad = direct;
    REQUIRE_THROWS(
        gra::spin::ValidateDirectTMatrix(c_bad, rho0, f2_reggeon, pomeron0, true, "CON_MP.json", true, false));
  }
}

TEST_CASE("PARAM_CON MULTI auto fully permutes neutral particles", "[gra::MRegge][continuum][symmetry][physics]") {
  std::vector<gra::MDecayBranch> pi0(4);
  for (auto &branch : pi0) {
    branch.p          = ToyParticle("pi0", 111, 0, 0.13498);
    branch.p.chargeX3 = 0;
  }
  CHECK(gra::regge::PermCount(pi0, pi0.size(), gra::regge::PermType::Auto) == 1);
  CHECK(gra::regge::PermCount(pi0, pi0.size(), gra::regge::PermType::Charged) == 0);
  CHECK(gra::regge::PermCount(pi0, pi0.size(), gra::regge::PermType::All) == 1);

  std::vector<gra::MDecayBranch> charged(4);
  const std::array<int, 4>       pdgs = {211, -211, 321, -321};
  for (std::size_t i = 0; i < charged.size(); ++i) {
    charged[i].p          = ToyParticle("charged", pdgs[i], 0, 0.2);
    charged[i].p.chargeX3 = (pdgs[i] > 0) ? 3 : -3;
  }
  CHECK(gra::regge::PermCount(charged, charged.size(), gra::regge::PermType::Auto) == 0);

  std::vector<gra::MDecayBranch> distinct_neutral(4);
  const std::array<int, 4>       neutral_pdgs = {111, 221, 113, 223};
  for (std::size_t i = 0; i < distinct_neutral.size(); ++i) {
    distinct_neutral[i].p          = ToyParticle("neutral", neutral_pdgs[i], 0, 0.2);
    distinct_neutral[i].p.chargeX3 = 0;
  }
  CHECK(gra::regge::PermCount(distinct_neutral, distinct_neutral.size(), gra::regge::PermType::Auto) == 1);
}

TEST_CASE(
    "GP production C and naturality checks cover analytic PP RP and OP "
    "exchange IDs",
    "[gra::spin][GP]") {
  gra::MODELPARAM = "TUNE0";
  MRandom rng;
  auto    f0  = gra::resonance::Read("RES/f0_500.json", rng, gra::ReggeProductionModel::GP);
  auto    rho = gra::resonance::Read("RES/rho_770.json", rng, gra::ReggeProductionModel::GP);

  SECTION("C-even f0 accepts analytic Pomeron-Pomeron") {
    ToyHelicityProcess proc;
    ConfigureToyProductionProcess(proc, "GP", "RES", "pi+ pi-");
    f0.GP.channels.front().exchange = {990, 990};
    proc.SetResonances({{"f0_500", f0}});
    REQUIRE_NOTHROW(proc.InitializeProcessAmplitude());
  }

  SECTION("C-even f0 accepts analytic f2-Reggeon-Pomeron") {
    ToyHelicityProcess proc;
    ConfigureToyProductionProcess(proc, "GP", "RES", "pi+ pi-");
    f0.GP.channels.front().exchange = {9910, 990};
    proc.SetResonances({{"f0_500", f0}});

    REQUIRE_NOTHROW(proc.InitializeProcessAmplitude());

    const auto configured = proc.GetResonances().at("f0_500");
    REQUIRE(configured.production.size() == 2);
    REQUIRE(configured.production.size() == 2);
    CHECK(configured.production[0].tree[0].p.pdg == 9910);
    CHECK(configured.production[0].tree[1].p.pdg == 990);
    CHECK(configured.production[1].tree[0].p.pdg == 990);
    CHECK(configured.production[1].tree[1].p.pdg == 9910);
  }

  SECTION("C-odd rho accepts analytic Odderon-Pomeron") {
    ToyHelicityProcess proc;
    ConfigureToyProductionProcess(proc, "GP", "RES", "pi+ pi-");
    auto rho_odd = gra::resonance::Read("RES/rho_770_odd.json", proc.state.random, gra::ReggeProductionModel::GP);
    proc.SetResonances({{"rho_770_odd", rho_odd}});
    REQUIRE_NOTHROW(proc.InitializeProcessAmplitude());
  }

  SECTION("C-even f0 rejects analytic Odderon-Pomeron") {
    ToyHelicityProcess proc;
    ConfigureToyProductionProcess(proc, "GP", "RES", "pi+ pi-");
    f0.GP.channels.front().exchange = {9990, 990};
    proc.SetResonances({{"f0_500", f0}});
    REQUIRE_THROWS(proc.InitializeProcessAmplitude());
  }

  SECTION("C-odd rho rejects analytic Pomeron-Pomeron") {
    ToyHelicityProcess proc;
    ConfigureToyProductionProcess(proc, "GP", "RES", "pi+ pi-");
    rho.GP.channels.front().exchange = {990, 990};
    proc.SetResonances({{"rho_770", rho}});
    REQUIRE_THROWS(proc.InitializeProcessAmplitude());
  }

  SECTION("analytic f2-Reggeon-Pomeron rejects wrong naturality LS rows") {
    ToyHelicityProcess proc;
    ConfigureToyProductionProcess(proc, "GP", "RES", "pi+ pi-");
    auto &channel    = f0.GP.channels.front();
    channel.exchange = {9910, 990};
    channel.basis    = gra::ReggeVertexBasis::LS;
    channel.g_ls.Clear();
    channel.g_ls.Set(1, 2, 1.0);
    proc.SetResonances({{"f0_500", f0}});
    REQUIRE_THROWS(proc.InitializeProcessAmplitude());
  }
}

TEST_CASE("GP rejects old resonance basis names", "[gra::spin][GP][schema]") {
  ToyHelicityProcess process;
  ConfigureToyProductionProcess(process, "GP", "RES", "pi+ pi-");
  auto rho                      = gra::resonance::Read("RES/rho_770.json", process.state.random, gra::ReggeProductionModel::GP);
  rho.GP.channels.front().basis = static_cast<gra::ReggeVertexBasis>(999);
  process.SetResonances({{"rho_770", rho}});
  REQUIRE_THROWS(process.InitializeProcessAmplitude());
}

TEST_CASE("GP fusion helicity parity orbit is symmetry complete", "[gra::spin][GP][parity]") {
  gra::MODELPARAM = "TUNE0";

  SECTION("one representative generates its reflected row") {
    ToyHelicityProcess proc;
    ConfigureToyProductionProcess(proc, "GP", "RES", "pi+ pi-");
    auto rho                           = gra::resonance::Read("RES/rho_770.json", proc.state.random, gra::ReggeProductionModel::GP);
    rho.GP.channels.front().helicity   = {{-1.0, 0.0}};
    rho.GP.channels.front().g_helicity = {1.0};
    proc.SetResonances({{"rho_770", rho}});
    REQUIRE_NOTHROW(proc.InitializeProcessAmplitude());
    const auto resonances = proc.GetResonances();
    REQUIRE(resonances.count("rho_770") == 1);
    const auto &vertices = resonances.at("rho_770").production;
    REQUIRE(vertices.size() == 2);
    for (const auto &production : vertices) {
      const auto &vertex = production.hel;
      std::size_t active = 0;
      for (std::size_t row = 0; row < vertex.T_set.size_row(); ++row) {
        for (std::size_t column = 0; column < vertex.T_set.size_col(); ++column) {
          active += vertex.T_set[row][column] ? 1U : 0U;
          if (!vertex.T_set[row][column]) { continue; }
          const std::size_t parity_row    = vertex.T_set.size_row() - 1 - row;
          const std::size_t parity_column = vertex.T_set.size_col() - 1 - column;
          REQUIRE(vertex.T_set[parity_row][parity_column]);
          RequireComplexNear(vertex.T[parity_row][parity_column], vertex.T[row][column], 1.0e-12);
        }
      }
      CHECK(active == 2);
    }
  }

  SECTION("a second representative from the same orbit is rejected") {
    ToyHelicityProcess proc;
    ConfigureToyProductionProcess(proc, "GP", "RES", "pi+ pi-");
    auto rho                           = gra::resonance::Read("RES/rho_770.json", proc.state.random, gra::ReggeProductionModel::GP);
    rho.GP.channels.front().helicity   = {{-1.0, 0.0}, {1.0, 0.0}};
    rho.GP.channels.front().g_helicity = {1.0, std::polar(1.0, 0.4)};
    proc.SetResonances({{"rho_770", rho}});
    REQUIRE_THROWS(proc.InitializeProcessAmplitude());
  }

  SECTION("an odd Jacob-Wick phase generates the negative partner") {
    gra::LORENTZSCALAR lts = DirectCentralPairLTSForTest(211, -211);
    gra::MRegge        regge(lts, gra::MModelTune::Load(modelfile),
                             gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));
    const auto         param = ReggeParametersForTest(regge, lts);

    gra::MParticle mother                  = lts.PDG.FindByPDG(113);
    mother.P                               = 1;
    const std::vector<gra::MParticle> legs = {lts.PDG.FindByPDG(gra::PDG::PDG_gamma), lts.PDG.FindByPDG(990)};
    gra::RES_PRODUCTION_CHANNEL       model;
    model.exchange   = {gra::PDG::PDG_gamma, 990};
    model.basis      = gra::ReggeVertexBasis::Helicity;
    model.C_symmetry = false;
    model.P_symmetry = true;
    model.helicity   = {{-1.0, 0.0}};
    model.g_helicity = {1.0};

    const auto        vertex   = gra::gpom::PrepareResonance(mother, legs, model, *param, lts.PDG, 2, 0.0);
    const std::size_t negative = gra::gpom::AnalyticMIndex(-1, 2, "odd GP parity negative");
    const std::size_t zero     = gra::gpom::AnalyticMIndex(0, 2, "odd GP parity zero");
    const std::size_t positive = gra::gpom::AnalyticMIndex(1, 2, "odd GP parity positive");
    REQUIRE(vertex.T_set[negative][zero]);
    REQUIRE(vertex.T_set[positive][zero]);
    RequireComplexNear(vertex.T[positive][zero], -vertex.T[negative][zero], 1.0e-12);
  }
}

// Check direct GP tensors against the same physical pole LS operators
TEST_CASE("GP fusion LS and direct helicity are pole duals for all exchange pairs",
          "[gra::spin][GP][exchange][normalization][physics]") {
  gra::LORENTZSCALAR           lts = DirectCentralPairLTSForTest(211, -211);
  gra::MRegge                  regge(lts, gra::MModelTune::Load(modelfile),
                                     gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));
  const auto                   param     = ReggeParametersForTest(regge, lts);
  constexpr std::array<int, 4> exchanges = {990, 9910, 9930, 9990};

  for (const int pdg1 : exchanges) {
    for (const int pdg2 : exchanges) {
      CAPTURE(pdg1, pdg2);
      const auto leg1   = lts.PDG.FindByPDG(pdg1);
      const auto leg2   = lts.PDG.FindByPDG(pdg2);
      const auto pole1  = gra::regge::PoleRepresentative(*param, lts.PDG, pdg1);
      const auto pole2  = gra::regge::PoleRepresentative(*param, lts.PDG, pdg2);
      const int  pair_c = pole1.C * pole2.C;
      auto       mother = lts.PDG.FindByPDG(pair_c < 0 ? 113 : 225);
      mother.C          = pair_c;

      const auto operators =
          gra::spin::CanonicalPoleOperators(mother, pole1, pole2, true, true, true, gra::spin::VertexContext::Auto);
      REQUIRE_FALSE(operators.empty());
      std::vector<gra::spin::LSTerm> terms;
      terms.reserve(operators.size());
      const std::complex<double> common_phase = std::polar(1.0, 0.37);
      for (const auto &i : indices(operators)) {
        const auto                &row = operators[i].coupling;
        const std::complex<double> coefficient =
            common_phase *
            (i == 0 ? std::complex<double>(0.71, 0.0)
                    : std::polar(0.19 + 0.03 * static_cast<double>(i), -0.41 + 0.09 * static_cast<double>(i)));
        terms.push_back({row.l, row.two_s, coefficient});
      }
      const auto     pole_vertex = gra::spin::PreparePoleLS(mother, pole1, pole2, terms, 1.0, true, true, true,
                                                            gra::spin::VertexContext::Auto, 0.0, false);
      const auto     expected    = gra::spin::PoleLSReduced(pole_vertex, 1.0);
      gra::HELMatrix pole_direct;
      pole_direct.C_symmetry     = true;
      pole_direct.P_symmetry     = true;
      pole_direct.coupling_basis = gra::CouplingBasis::Helicity;
      pole_direct.T              = expected;
      pole_direct.T_set          = MMatrix<bool>(expected.size_row(), expected.size_col(), true);
      const auto expansion       = gra::spin::DirectHelicityToLSCoefficients(pole_direct, mother, pole1, pole2, true,
                                                                             gra::spin::VertexContext::Auto);
      REQUIRE(expansion.rows.size() == terms.size());
      for (const auto &i : indices(terms)) {
        CAPTURE(i, terms[i].l, terms[i].two_s, terms[i].coefficient, expansion.rows[i].l, expansion.rows[i].two_s,
                expansion.coefficients[i]);
        REQUIRE(expansion.rows[i].l == terms[i].l);
        REQUIRE(expansion.rows[i].two_s == terms[i].two_s);
        RequireComplexNear(
            expansion.coefficients[i],
            terms[i].coefficient * gra::spin::RawPoleLSNormalization(pole1.spinX2, pole2.spinX2, mother.spinX2,
                                                                     terms[i].l, static_cast<int>(terms[i].two_s)),
            2.0e-11);
      }
      const auto model = DirectGPFromPole({pdg1, pdg2}, expected);
      REQUIRE_FALSE(model.helicity.empty());
      if (pdg1 == pdg2) {
        MMatrix<std::complex<double>> reconstructed(expected.size_row(), expected.size_col(), 0.0);
        for (const auto &i : indices(model.helicity)) {
          const int                  m1    = static_cast<int>(model.helicity[i][0]);
          const int                  m2    = static_cast<int>(model.helicity[i][1]);
          const std::complex<double> value = model.g_helicity[i];
          for (const auto &[a, b] : {std::pair<int, int>{m1, m2}, {-m1, -m2}, {m2, m1}, {-m2, -m1}}) {
            reconstructed[static_cast<std::size_t>(a + pole1.spinX2 / 2)]
                         [static_cast<std::size_t>(b + pole2.spinX2 / 2)] = value;
          }
        }
        RequireMatrixNear(reconstructed, expected, 2.0e-11);
      }

      const auto actual = gra::gpom::PrepareResonance(mother, {leg1, leg2}, model, *param, lts.PDG, 2, 0.0);
      const int  j1     = pole1.spinX2 / 2;
      const int  j2     = pole2.spinX2 / 2;
      for (int m1 = -j1; m1 <= j1; ++m1) {
        for (int m2 = -j2; m2 <= j2; ++m2) {
          const std::size_t source1 = static_cast<std::size_t>(m1 + j1);
          const std::size_t source2 = static_cast<std::size_t>(m2 + j2);
          const std::size_t target1 = gra::gpom::AnalyticMIndex(m1, 2, "GP pole dual target 1");
          const std::size_t target2 = gra::gpom::AnalyticMIndex(m2, 2, "GP pole dual target 2");
          RequireComplexNear(actual.T[target1][target2], expected[source1][source2], 2.0e-11);
        }
      }
    }
  }
}

// Check raw GP rows remain independent when parity conservation is disabled
TEST_CASE("GP fusion helicity keeps independent rows without parity symmetry",
          "[gra::MRegge][GP][spin][parity][regression]") {
  gra::LORENTZSCALAR lts = DirectCentralPairLTSForTest(211, -211);
  gra::MRegge        regge(lts, gra::MModelTune::Load(modelfile),
                           gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));
  const auto         param = ReggeParametersForTest(regge, lts);

  gra::MParticle mother                  = lts.PDG.FindByPDG(113);
  mother.P                               = 0;
  const std::vector<gra::MParticle> legs = {lts.PDG.FindByPDG(gra::PDG::PDG_gamma), lts.PDG.FindByPDG(990)};
  gra::RES_PRODUCTION_CHANNEL       model;
  model.exchange   = {gra::PDG::PDG_gamma, 990};
  model.basis      = gra::ReggeVertexBasis::Helicity;
  model.C_symmetry = false;
  model.P_symmetry = false;
  model.helicity   = {{-1.0, 0.0}, {1.0, 0.0}};
  model.g_helicity = {{0.73, -0.21}, {-0.34, 0.58}};

  const auto        vertex   = gra::gpom::PrepareResonance(mother, legs, model, *param, lts.PDG, 2, 0.0);
  const std::size_t negative = gra::gpom::AnalyticMIndex(-1, 2, "unconstrained GP negative");
  const std::size_t zero     = gra::gpom::AnalyticMIndex(0, 2, "unconstrained GP zero");
  const std::size_t positive = gra::gpom::AnalyticMIndex(1, 2, "unconstrained GP positive");
  REQUIRE(vertex.T_set[negative][zero]);
  REQUIRE(vertex.T_set[positive][zero]);
  RequireComplexNear(vertex.T[negative][zero], model.g_helicity[0], 1.0e-12);
  RequireComplexNear(vertex.T[positive][zero], model.g_helicity[1], 1.0e-12);
}

TEST_CASE("GP rejects resonance rows above MMAX", "[gra::spin][GP][schema]") {
  ToyHelicityProcess proc;
  ConfigureToyProductionProcess(proc, "GP", "RES", "pi+ pi-");
  proc.state.lts.process.MMAX = 2;
  auto  rho                   = gra::resonance::Read("RES/rho_770.json", proc.state.random, gra::ReggeProductionModel::GP);
  auto &channel               = rho.GP.channels.front();
  channel.basis               = gra::ReggeVertexBasis::Helicity;
  channel.helicity            = {{-1.0, 3.0}};
  channel.g_helicity          = {1.0};
  proc.SetResonances({{"rho_770", rho}});
  REQUIRE_THROWS(proc.InitializeProcessAmplitude());
}

TEST_CASE("GP resonance keeps Regge LS coefficients without a pole cache", "[gra::spin][GP]") {
  ToyHelicityProcess proc;
  gra::MODELPARAM = "TUNE0";
  proc.SetProcessForTest("GP", "RES");
  proc.state.lts.PDG   = LoadedPDGTable();
  proc.state.lts.beam1 = proc.state.lts.PDG.FindByPDG(2212);
  proc.state.lts.beam2 = proc.state.lts.PDG.FindByPDG(2212);
  proc.SetDecayMode("pi+ pi-");

  gra::PARAM_RES f2        = gra::resonance::Read("RES/f2_2150.json", proc.state.random, gra::ReggeProductionModel::GP);
  auto          &couplings = f2.GP.channels.front().g_ls;
  for (auto &term : couplings) { term.coefficient = 0.1; }
  REQUIRE_FALSE(couplings.Empty());
  couplings.begin()->coefficient = 0.5 * proc.GetModelTune()->Global().coupling_min;
  proc.SetResonances({{"f2_2150", f2}});

  REQUIRE_NOTHROW(proc.InitializeProcessAmplitude());

  const auto resonances = proc.GetResonances();
  REQUIRE(resonances.count("f2_2150") == 1);
  const auto &vertices = resonances.at("f2_2150").production;
  REQUIRE(vertices.size() == 1);
  CHECK(vertices.front().hel.UsesLSCouplings());
  CHECK(vertices.front().hel.alpha_ls.IsFinite());
  REQUIRE(vertices.front().hel.alpha_ls.Size() + 1 == couplings.Size());
  for (const auto &term : vertices.front().hel.alpha_ls) {
    CHECK(std::abs(term.coefficient) > proc.GetModelTune()->Global().coupling_min);
  }
}



TEST_CASE("GP raw resonance retains an m equals three amplitude", "[gra::MRegge][spin][GP][RES]") {
  gra::LORENTZSCALAR lts = MakeToyReggeLTSAsymmetric(0.31, 4.2, -4.6);
  lts.process.MMAX       = 3;
  gra::PARAM_RES res     = MakeToyScalarGPResonance(lts.process.MMAX);
  auto          &hel     = res.production.front().hel;
  hel.coupling_basis     = gra::CouplingBasis::Helicity;
  hel.analytic_MMAX      = lts.process.MMAX;
  hel.T                  = MMatrix<std::complex<double>>(7, 7, 0.0);
  hel.T_set              = MMatrix<bool>(7, 7, false);
  const std::size_t m3   = gra::gpom::AnalyticMIndex(3, lts.process.MMAX, "raw m=3 test");
  hel.T[m3][m3]          = std::polar(0.73, 0.21);
  hel.T_set[m3][m3]      = true;
  gra::PruneHelicityCouplings(hel, 0.0, "GP raw m=3 test");
  REQUIRE(hel.T_active.size() == 1);
  CHECK(hel.T_active.front() == std::make_pair(m3, m3));

  std::size_t steered = 0;
  for (std::size_t row = 0; row < hel.T_set.size_row(); ++row) {
    for (std::size_t col = 0; col < hel.T_set.size_col(); ++col) { steered += hel.T_set[row][col] ? 1U : 0U; }
  }
  CHECK(steered == 1);

  gra::MRegge  regge(lts, gra::MModelTune::Load(modelfile),
                     gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Resonance, "GP_raw_m3"));
  const double amp2 = TestReggeRes(regge, lts, res, gra::ReggeProductionModel::GP);
  CHECK(std::isfinite(amp2));
  CHECK(amp2 > 0.0);
}

// Check the dedicated multi-Regge steering block and local form factors
TEST_CASE("PARAM_CON MULTI controls the continuum ladder structure", "[gra::MRegge][continuum][physics]") {
  const auto        tune = WriteModifiedPhotoVMTune("multiregge_enabled", [](auto &j) {
    auto &multi                                                  = j.at("PARAM_REGGE").at("PARAM_CON").at("MULTI");
    multi.at("FF_transfer_ext")                                  = false;
    multi.at("FF_transfer_int")                                  = true;
    multi.at("secondary_exchanges")                              = true;
    auto &channels                                               = j.at("PARAM_REGGE").at("PARAM_CON").at("MP").at("[211,-211]");
    channels.clear();
    const std::array<int, 3> exchanges = {995, 9915, 9933};
    for (const auto &upper : indices(exchanges)) {
      for (std::size_t lower = upper; lower < exchanges.size(); ++lower) {
        channels.push_back({exchanges[upper], exchanges[lower]});
      }
    }
  }, {}, [](auto &card) { SetContinuumField(card, "[211,211]", "reggeize", {{"active", false}, {"freeze_scale2", 1.0}}); });
  const std::string data = gra::aux::GetInputData(tune.second);
  const auto        card = nlohmann::json::parse(data);
  REQUIRE(card.at("PARAM_REGGE").at("PARAM_CON").contains("MULTI"));
  const auto param = gra::regge::ReadParam({211, -211}, LoadedPDGTable(), *gra::MModelTune::Load(tune.second));
  CHECK_FALSE(param.con.multiregge_transfer_ext);
  CHECK(param.con.multiregge_transfer_int);
  CHECK(param.con.multiregge_secondary_exchanges);
  CHECK(param.con.permutations == gra::regge::PermType::Auto);
  REQUIRE(param.con.multiregge_topologies.size() == 2);
  CHECK(param.con.multiregge_topologies.at(4) == std::vector<gra::regge::Topology>{{4}, {2, 2}});
  CHECK(param.con.multiregge_topologies.at(6) == std::vector<gra::regge::Topology>{{6}, {4, 2}, {2, 2, 2}});

  const auto *entry = gra::regge::FindPair(param, {211, -211}, gra::ReggeProductionModel::MP);
  REQUIRE(entry != nullptr);

  // Check every connected four-pion Pomeron and Reggeon ladder
  std::set<std::vector<std::size_t>> four_pion_chains;
  for (const auto &first : entry->channels) {
    for (const auto &second : entry->channels) {
      const std::size_t first_top = gra::regge::TrajectoryIndex(param, first.first);
      const std::size_t middle    = gra::regge::TrajectoryIndex(param, first.second);
      if (middle != gra::regge::TrajectoryIndex(param, second.first)) { continue; }
      four_pion_chains.insert({first_top, middle, gra::regge::TrajectoryIndex(param, second.second)});
    }
  }
  std::set<std::vector<std::size_t>> expected_four_pion_chains;
  for (std::size_t upper = 0; upper < 3; ++upper) {
    for (std::size_t middle = 0; middle < 3; ++middle) {
      for (std::size_t lower = 0; lower < 3; ++lower) { expected_four_pion_chains.insert({upper, middle, lower}); }
    }
  }
  CHECK(four_pion_chains == expected_four_pion_chains);

  // Check that the same connected-channel rule extends to the six-pion ladder
  std::set<std::vector<std::size_t>> six_pion_chains;
  for (const auto &first : entry->channels) {
    for (const auto &second : entry->channels) {
      const std::size_t first_top     = gra::regge::TrajectoryIndex(param, first.first);
      const std::size_t first_bottom  = gra::regge::TrajectoryIndex(param, first.second);
      const std::size_t second_top    = gra::regge::TrajectoryIndex(param, second.first);
      const std::size_t second_bottom = gra::regge::TrajectoryIndex(param, second.second);
      if (first_bottom != second_top) { continue; }
      for (const auto &third : entry->channels) {
        const std::size_t third_top = gra::regge::TrajectoryIndex(param, third.first);
        if (second_bottom != third_top) { continue; }
        six_pion_chains.insert(
            {first_top, first_bottom, second_bottom, gra::regge::TrajectoryIndex(param, third.second)});
      }
    }
  }
  std::set<std::vector<std::size_t>> expected_six_pion_chains;
  for (std::size_t upper = 0; upper < 3; ++upper) {
    for (std::size_t first_middle = 0; first_middle < 3; ++first_middle) {
      for (std::size_t second_middle = 0; second_middle < 3; ++second_middle) {
        for (std::size_t lower = 0; lower < 3; ++lower) {
          expected_six_pion_chains.insert({upper, first_middle, second_middle, lower});
        }
      }
    }
  }
  CHECK(six_pion_chains == expected_six_pion_chains);

  gra::LORENTZSCALAR lts;
  lts.PDG = LoadedPDGTable();
  gra::MDecayBranch left;
  gra::MDecayBranch right;
  left.p   = lts.PDG.FindByPDG(211);
  right.p  = lts.PDG.FindByPDG(-211);
  left.p4  = gra::M4Vec(0.2, 0.0, 0.6, 1.1);
  right.p4 = gra::M4Vec(-0.1, 0.0, -0.3, 0.9);

  const double               meson_t  = -2.3;
  const double               pion_m2  = gra::math::pow2(left.p.mass);
  const auto                &vertex   = entry->channels.front();
  const double               upper_offshell = gra::regge::FormFactor(meson_t, pion_m2, vertex.offshell[0]);
  const double               lower_offshell = gra::regge::FormFactor(meson_t, pion_m2, vertex.offshell[1]);
  const std::complex<double> meson_prop =
      gra::regge::OffshellProp(meson_t, pion_m2, entry->reggeize, gra::regge::MesonTraj{}, left.p4, right.p4);
  const std::complex<double> meson_exchange =
      gra::regge::MesonExchange(param, meson_t, pion_m2, *entry, vertex, left, right);
  CHECK(std::abs(meson_exchange - upper_offshell * meson_prop * lower_offshell) < 1.0e-14);
  const double upper_t = -0.31;
  const double lower_t = -0.47;
  const double expected_vertex =
      gra::regge::VertexSign(vertex, left, right, param) *
      gra::regge::TransferFF(upper_t, vertex.transfer[0]) *
      gra::regge::TransferFF(lower_t, vertex.transfer[1]);
  const double vertex_sign = gra::regge::VertexSign(vertex, left, right, param);
  const double upper_ff    = gra::regge::TransferFF(upper_t, vertex.transfer[0]);
  const double lower_ff    = gra::regge::TransferFF(lower_t, vertex.transfer[1]);
  CHECK(gra::regge::PairVertexFactor(param, vertex, left, right, upper_t, lower_t, true, true) ==
        Approx(expected_vertex).margin(1.0e-14));
  CHECK(gra::regge::PairVertexFactor(param, vertex, left, right, upper_t, lower_t, true, false) ==
        Approx(vertex_sign * upper_ff).margin(1.0e-14));
  CHECK(gra::regge::PairVertexFactor(param, vertex, left, right, upper_t, lower_t, false, true) ==
        Approx(vertex_sign * lower_ff).margin(1.0e-14));
  CHECK(gra::regge::PairVertexFactor(param, vertex, left, right, upper_t, lower_t, false, false) ==
        Approx(vertex_sign).margin(1.0e-14));
  const std::complex<double> pair_exchange = gra::regge::PairExchange(
      param, *entry, vertex, left, right, meson_t, pion_m2, upper_t, lower_t);
  CHECK(std::abs(pair_exchange - meson_exchange * expected_vertex) < 1.0e-14);
  CHECK_THROWS_AS(gra::regge::OffshellProp(pion_m2, pion_m2, {false, 1.0}, gra::regge::MesonTraj{}, left.p4, right.p4),
                  gra::AmplitudeFailure);
}

// Check that parallel production owns its transverse integration rule
TEST_CASE("Parallel Regge quadrature is independent of screening numerics", "[gra::MRegge][continuum][numerics]") {
  auto document = nlohmann::json::parse(gra::aux::GetInputData(gra::ResolveModelDataFile("TUNE0", "NUMERICS.json")));
  document["NUMERICS_EIKONAL"]["LOOP_INTEGRAL"]["NumberLoopKT"]  = 3;
  document["NUMERICS_EIKONAL"]["LOOP_INTEGRAL"]["NumberLoopPHI"] = 4;
  auto &parallel                                                 = document["NUMERICS_REGGE"]["LOOP_INTEGRAL"];
  parallel["kT_integrator"]                                      = "GL";
  parallel["phi_integrator"]                                     = "Trap";
  parallel["log_kT"]                                             = false;
  parallel["MinKT"]                                              = 0.0;
  parallel["MaxKT"]                                              = 2.0;
  parallel["NumberKT"]                                           = 5;
  parallel["NumberPHI"]                                          = 7;
  parallel["MaxBT"]                                              = 3.0;
  parallel["NumberBT"]                                           = 11;

  gra::MReggeNumerics numerics;
  numerics.ConfigureFromJson("independent Regge quadrature", document.dump());
  const auto &rule = numerics.ParallelRule();
  REQUIRE(rule.kt.size() == 5);
  REQUIRE(rule.radial_weight.size() == 5);
  REQUIRE(rule.measure_weight.size_row() == 5);
  REQUIRE(rule.measure_weight.size_col() == 7);
  REQUIRE(rule.kt_x.size_row() == 5);
  REQUIRE(rule.kt_x.size_col() == 7);
  REQUIRE(numerics.ParallelImpactCount() == 11);
  CHECK(numerics.ParallelImpactMax() == Approx(3.0 / gra::PDG::GeV2fm));

  double disk_measure = 0.0;
  for (const auto &radial : gra::aux::indices(rule.kt)) {
    for (std::size_t azimuth = 0; azimuth < 7; ++azimuth) { disk_measure += rule.measure_weight[radial][azimuth]; }
  }
  CHECK(disk_measure == Approx(4.0 * gra::math::PI).epsilon(1.0e-12));

  for (const std::string field : {"NumberKT", "NumberPHI", "NumberBT"}) {
    for (const nlohmann::json &value : {nlohmann::json(2.5), nlohmann::json(0), nlohmann::json(-1),
                                     nlohmann::json(static_cast<unsigned long long>(std::numeric_limits<unsigned int>::max()) + 4ULL)}) {
      CAPTURE(field, value);
      auto invalid_document = document;
      invalid_document["NUMERICS_REGGE"]["LOOP_INTEGRAL"][field] = value;
      gra::MReggeNumerics invalid_count;
      CHECK_THROWS_AS(invalid_count.ConfigureFromJson("invalid Regge count", invalid_document.dump()), std::invalid_argument);
    }
  }

  parallel["kT_integrator"] = "Boole";
  gra::MReggeNumerics invalid;
  CHECK_THROWS_AS(invalid.ConfigureFromJson("invalid Regge quadrature", document.dump()), std::invalid_argument);
}

// Check that malformed multi-Regge parameters fail during tune initialization
TEST_CASE("PARAM_CON MULTI rejects invalid steering", "[gra::MRegge][params][validation]") {
  SECTION("removed suppression") {
    const auto tune = WriteModifiedPhotoVMTune("multiregge_removed_suppression", [](auto &j) {
      j.at("PARAM_REGGE").at("PARAM_CON").at("MULTI")["suppression"] = {{"W0", 1.0}, {"a", 0.1}};
    });
    REQUIRE_THROWS(gra::regge::ReadParam({211, -211}, LoadedPDGTable(), *gra::MModelTune::Load(tune.second)));
  }

  SECTION("permutation mode") {
    const auto tune = WriteModifiedPhotoVMTune(
        "multiregge_bad_permutations", [](auto &j) { j.at("PARAM_REGGE").at("PARAM_CON").at("MULTI").at("permutations") = "neutral"; });
    REQUIRE_THROWS(gra::regge::ReadParam({211, -211}, LoadedPDGTable(), *gra::MModelTune::Load(tune.second)));
  }

  SECTION("secondary exchange switch") {
    const auto tune = WriteModifiedPhotoVMTune("multiregge_bad_secondary_exchanges", [](auto &j) {
      j.at("PARAM_REGGE").at("PARAM_CON").at("MULTI").at("secondary_exchanges") = 1;
    });
    REQUIRE_THROWS(gra::regge::ReadParam({211, -211}, LoadedPDGTable(), *gra::MModelTune::Load(tune.second)));
  }

  SECTION("external transfer switch") {
    const auto tune = WriteModifiedPhotoVMTune(
        "multiregge_bad_transfer_ext", [](auto &j) { j.at("PARAM_REGGE").at("PARAM_CON").at("MULTI").at("FF_transfer_ext") = 1; });
    REQUIRE_THROWS(gra::regge::ReadParam({211, -211}, LoadedPDGTable(), *gra::MModelTune::Load(tune.second)));
  }

  SECTION("internal transfer switch") {
    const auto tune = WriteModifiedPhotoVMTune(
        "multiregge_bad_transfer_int", [](auto &j) { j.at("PARAM_REGGE").at("PARAM_CON").at("MULTI").at("FF_transfer_int") = 1; });
    REQUIRE_THROWS(gra::regge::ReadParam({211, -211}, LoadedPDGTable(), *gra::MModelTune::Load(tune.second)));
  }

  SECTION("empty topology bank") {
    const auto tune = WriteModifiedPhotoVMTune("multiregge_empty_topology", [](auto &j) {
      j.at("PARAM_REGGE").at("PARAM_CON").at("MULTI").at("partitions").at("4") = nlohmann::json::array();
    });
    REQUIRE_THROWS(gra::regge::ReadParam({211, -211}, LoadedPDGTable(), *gra::MModelTune::Load(tune.second)));
  }

  SECTION("empty topology") {
    const auto tune = WriteModifiedPhotoVMTune("multiregge_empty_partition", [](auto &j) {
      j.at("PARAM_REGGE").at("PARAM_CON").at("MULTI").at("partitions").at("4") = nlohmann::json::array({nlohmann::json::array()});
    });
    REQUIRE_THROWS(gra::regge::ReadParam({211, -211}, LoadedPDGTable(), *gra::MModelTune::Load(tune.second)));
  }

  SECTION("unsupported topology block") {
    const auto tune = WriteModifiedPhotoVMTune("multiregge_bad_partition_block", [](auto &j) {
      j.at("PARAM_REGGE").at("PARAM_CON").at("MULTI").at("partitions").at("4") = nlohmann::json::array({nlohmann::json::array({3, 1})});
    });
    REQUIRE_THROWS(gra::regge::ReadParam({211, -211}, LoadedPDGTable(), *gra::MModelTune::Load(tune.second)));
  }

  SECTION("ascending topology") {
    const auto tune = WriteModifiedPhotoVMTune("multiregge_ascending_partition", [](auto &j) {
      j.at("PARAM_REGGE").at("PARAM_CON").at("MULTI").at("partitions").at("6") = nlohmann::json::array({nlohmann::json::array({2, 4})});
    });
    REQUIRE_THROWS(gra::regge::ReadParam({211, -211}, LoadedPDGTable(), *gra::MModelTune::Load(tune.second)));
  }

  SECTION("wrong topology sum") {
    const auto tune = WriteModifiedPhotoVMTune("multiregge_bad_partition_sum", [](auto &j) {
      j.at("PARAM_REGGE").at("PARAM_CON").at("MULTI").at("partitions").at("6") = nlohmann::json::array({nlohmann::json::array({4})});
    });
    REQUIRE_THROWS(gra::regge::ReadParam({211, -211}, LoadedPDGTable(), *gra::MModelTune::Load(tune.second)));
  }

  SECTION("duplicate topology") {
    const auto tune = WriteModifiedPhotoVMTune("multiregge_duplicate_partition", [](auto &j) {
      j.at("PARAM_REGGE").at("PARAM_CON").at("MULTI").at("partitions").at("4") =
          nlohmann::json::array({nlohmann::json::array({4}), nlohmann::json::array({4})});
    });
    REQUIRE_THROWS(gra::regge::ReadParam({211, -211}, LoadedPDGTable(), *gra::MModelTune::Load(tune.second)));
  }
}

// Check that the multi-Regge switch acts only while building direct ladders
TEST_CASE("PARAM_CON MULTI preserves every two-body model channel", "[gra::MRegge][continuum][physics]") {
  const auto tune = WriteModifiedPhotoVMTune("multiregge_no_secondary_exchanges", [](auto &j) {
    j.at("PARAM_REGGE").at("PARAM_CON").at("MULTI").at("secondary_exchanges") = false;
    SetSecondaryContinuumChannels(j);
  });
  const auto param =
      gra::regge::ReadParam({211, -211, 211, -211}, LoadedPDGTable(), *gra::MModelTune::Load(tune.second));
  CHECK_FALSE(param.con.multiregge_secondary_exchanges);

  for (const auto model :
       {gra::ReggeProductionModel::MP, gra::ReggeProductionModel::XP, gra::ReggeProductionModel::GP}) {
    const auto *entry = gra::regge::FindPair(param, {211, -211}, model);
    REQUIRE(entry != nullptr);
    REQUIRE(entry->channels.size() >= 2);

    bool found_pomeron   = false;
    bool found_secondary = false;
    for (const auto &channel : entry->channels) {
      const bool allowed   = gra::regge::VertexAllowed(param, channel);
      const bool secondary = gra::regge::IsSecondaryReggeonTrajectory(param, channel.first) ||
                             gra::regge::IsSecondaryReggeonTrajectory(param, channel.second);
      found_pomeron   = found_pomeron || (!secondary && allowed);
      found_secondary = found_secondary || secondary;
      if (secondary) { CHECK_FALSE(allowed); }
    }
    CHECK(found_pomeron);
    CHECK(found_secondary);
  }
}

// Check exact fixed-spin aliases and trajectory-level GP ladder connections
TEST_CASE("Multi-Regge internal exchanges follow each model spin basis",
          "[gra::MRegge][continuum][multiregge][exchange][physics]") {
  const auto param =
      gra::regge::ReadParam({211, -211, 211, -211, 211, -211}, LoadedPDGTable(), *gra::MModelTune::Load(modelfile));

  gra::regge::PairParam upper;
  gra::regge::PairParam middle;
  gra::regge::PairParam lower;
  upper.channels                                         = {{995, 991}};
  middle.channels                                        = {{993, 995}};
  lower.channels                                         = {{991, 995}};
  const std::vector<const gra::regge::PairParam *> two   = {&upper, &middle};
  const std::vector<const gra::regge::PairParam *> three = {&upper, &middle, &lower};

  CHECK_FALSE(gra::regge::ContinuumLadderEntriesConnect(param, two, gra::ReggeProductionModel::MP));
  CHECK_FALSE(gra::regge::ContinuumLadderEntriesConnect(param, two, gra::ReggeProductionModel::XP));
  CHECK(gra::regge::ContinuumLadderEntriesConnect(param, two, gra::ReggeProductionModel::GP));
  CHECK(gra::regge::ContinuumLadderEntriesConnect(param, three, gra::ReggeProductionModel::GP));

  middle.channels = {{991, 995}};
  CHECK(gra::regge::ContinuumLadderEntriesConnect(param, two, gra::ReggeProductionModel::MP));
  CHECK(gra::regge::ContinuumLadderEntriesConnect(param, two, gra::ReggeProductionModel::XP));

  middle.channels = {{9933, 995}};
  CHECK_FALSE(gra::regge::ContinuumLadderEntriesConnect(param, two, gra::ReggeProductionModel::GP));
}

TEST_CASE("PARAM_REGGE continuum lookup requires explicit final-state pairs", "[gra::MRegge][continuum]") {
  const auto param = gra::regge::ReadParam({211, -211, 321, -321}, LoadedPDGTable(), *gra::MModelTune::Load(modelfile));

  const auto *pipi = gra::regge::FindPair(param, {211, -211}, gra::ReggeProductionModel::MP);
  const auto *kk   = gra::regge::FindPair(param, {321, -321}, gra::ReggeProductionModel::MP);
  REQUIRE(pipi != nullptr);
  REQUIRE(kk != nullptr);
  CHECK(pipi->pdg == std::vector<int>{211, -211});
  CHECK(kk->pdg == std::vector<int>{321, -321});
  const auto pipi_ff       = pipi->forms.at(995).offshell;
  const auto kk_ff         = kk->forms.at(995).offshell;
  const int  pipi_reggeize = pipi->reggeize.active;
  const int  kk_reggeize   = kk->reggeize.active;

  CHECK(gra::regge::FindPair(param, {321, -211}, gra::ReggeProductionModel::MP) == nullptr);
  CHECK(gra::regge::FindPair(param, {-321, 211}, gra::ReggeProductionModel::MP) == nullptr);

  const auto param_reordered =
      gra::regge::ReadParam({321, 211, -321, -211}, LoadedPDGTable(), *gra::MModelTune::Load(modelfile));
  const auto *pipi_reordered = gra::regge::FindPair(param_reordered, {211, -211}, gra::ReggeProductionModel::MP);
  const auto *kk_reordered   = gra::regge::FindPair(param_reordered, {321, -321}, gra::ReggeProductionModel::MP);
  REQUIRE(pipi_reordered != nullptr);
  REQUIRE(kk_reordered != nullptr);
  CHECK(pipi_reordered->forms.at(995).offshell.type == pipi_ff.type);
  CHECK(kk_reordered->forms.at(995).offshell.type == kk_ff.type);
  CHECK(pipi_reordered->reggeize.active == pipi_reggeize);
  CHECK(kk_reordered->reggeize.active == kk_reggeize);
  REQUIRE(pipi_reordered->forms.at(995).offshell.param.size() == pipi_ff.param.size());
  REQUIRE(kk_reordered->forms.at(995).offshell.param.size() == kk_ff.param.size());
  for (std::size_t i = 0; i < pipi_ff.param.size(); ++i) {
    CHECK(pipi_reordered->forms.at(995).offshell.param[i] == Approx(pipi_ff.param[i]));
  }
  for (std::size_t i = 0; i < kk_ff.param.size(); ++i) {
    CHECK(kk_reordered->forms.at(995).offshell.param[i] == Approx(kk_ff.param[i]));
  }

  REQUIRE_THROWS(gra::regge::Pair(param, {321, -211}, gra::ReggeProductionModel::MP));
}

TEST_CASE("PARAM_REGGE rejects duplicate continuum final-state rows", "[gra::MRegge][continuum][validation]") {
  SECTION("obsolete model identifier") {
    const auto tune = WriteModifiedPhotoVMTune("continuum_obsolete_model_identifier",
                                               [](auto &j) { j.at("PARAM_REGGE").at("PARAM_CON")["MCON"] = j.at("PARAM_REGGE").at("PARAM_CON").at("MP"); });
    REQUIRE_THROWS(gra::regge::ReadParam({211, -211}, LoadedPDGTable(), *gra::MModelTune::Load(tune.second)));
  }

  SECTION("reversed pair key") {
    const auto tune = WriteModifiedPhotoVMTune("continuum_duplicate_active_pair", [](auto &j) {
      auto &gp         = j.at("PARAM_REGGE").at("PARAM_CON").at("GP");
      gp["[-211,211]"] = gp.at("[211,-211]");
    });
    REQUIRE_THROWS(gra::regge::ReadParam({211, -211}, LoadedPDGTable(), *gra::MModelTune::Load(tune.second)));
  }

  SECTION("obsolete default pair") {
    const auto tune = WriteModifiedPhotoVMTune("continuum_default_pair", [](auto &j) {
      j.at("PARAM_REGGE").at("PARAM_CON").at("GP")["default"] = {{990, 990}};
    });
    REQUIRE_THROWS(gra::regge::ReadParam({211, -211}, LoadedPDGTable(), *gra::MModelTune::Load(tune.second)));
  }
}

TEST_CASE("Continuum cards validate vertex and meson-line settings", "[gra::MRegge][continuum][validation]") {
  const auto mutation = [&](auto &card) {
    auto &pairs = card.begin().value();
    auto &block = pairs.at("[211,211]");
    SECTION("compact pair key") {
      pairs["[211, 211]"] = block;
      pairs.erase("[211,211]");
    }
    SECTION("selected form family") {
      block["FF_offshell"] = {{"type", "power"}, {"norm", "pole"}, {"Lambda2", 0.8}, {"n", 0.0}};
    }
    SECTION("veto scale") { block["pveto"]["M0"] = 0.0; }
    SECTION("invalid freezing virtuality") { block["reggeize"]["freeze_scale2"] = 0.0; }
    SECTION("missing freezing virtuality") { block["reggeize"].erase("freeze_scale2"); }
    SECTION("boolean reggeization") { block["reggeize"] = true; }
    SECTION("inconsistent freezing virtuality") { block["reggeize"]["freeze_scale2"] = 2.0; }
    SECTION("inconsistent shared propagator") { block["reggeize"]["active"] = !block.at("reggeize").at("active").template get<bool>(); }
  };
  const auto tune = WriteModifiedPhotoVMTune("continuum_invalid_vertex", {}, {}, mutation);
  REQUIRE_THROWS(gra::regge::ReadParam({211, -211}, LoadedPDGTable(), *gra::MModelTune::Load(tune.second)));
}

TEST_CASE("PARAM_REGGE continuum Reggeon vertices obey two-body quantum numbers", "[gra::MRegge][continuum][physics]") {
  const gra::MPDG &pdg   = LoadedPDGTable();
  const auto       param = gra::regge::ReadParam({211, -211}, pdg, *gra::MModelTune::Load(modelfile));

  const std::map<int, int> pion_pair_exchanges = {{9915, 9915}, {9925, 9925}, {9933, 9933}, {9943, 9943},
                                                  {9910, 9915}, {9920, 9925}, {9930, 9933}, {9940, 9943}};
  for (const auto &[exchange, representative] : pion_pair_exchanges) {
    const auto check = gra::regge::CheckVertex(param, pdg, exchange, {211, -211});
    CAPTURE(exchange);
    CHECK(check.applies);
    CHECK(check.allowed);
    CHECK(check.representative_pdg == representative);
  }

  for (const int exchange : {991, 993, 995, 9915, 9925, 990, 9910, 9920}) {
    CAPTURE(exchange);
    CHECK(gra::regge::CParity(param, exchange) == 1);
    CHECK(gra::math::IsExactEqual(gra::regge::AntiparticleSign(param, exchange, -321), 1.0));
  }
  for (const int exchange : {9933, 9943, 9991, 9993, 9997, 9930, 9940, 9990}) {
    CAPTURE(exchange);
    CHECK(gra::regge::CParity(param, exchange) == -1);
    CHECK(gra::math::IsExactEqual(gra::regge::AntiparticleSign(param, exchange, -321), -1.0));
  }
  CHECK(gra::regge::PairCParity(param, 995, 9915) == 1);
  CHECK(gra::regge::PairCParity(param, 995, 9933) == -1);
  CHECK(gra::regge::PairCParity(param, 9933, 9943) == 1);

  const auto charge_forbidden = gra::regge::CheckVertex(param, pdg, 9915, {211, 211});
  CHECK(charge_forbidden.applies);
  CHECK_FALSE(charge_forbidden.allowed);
  CHECK(charge_forbidden.reason.find("electric charge") != std::string::npos);

  for (const int exchange : {9915, 9910}) {
    const auto check = gra::regge::CheckVertex(param, pdg, exchange, {113, 113});
    CAPTURE(exchange);
    CHECK(check.applies);
    CHECK(check.allowed);
  }

  for (const int exchange : {9933, 9930}) {
    const auto check = gra::regge::CheckVertex(param, pdg, exchange, {113, 113});
    CAPTURE(exchange);
    CHECK(check.applies);
    CHECK_FALSE(check.allowed);
    CHECK(check.reason.find("J, P, C") != std::string::npos);
  }

  // Confirm that an invalid active channel is rejected while reading the tune
  const auto invalid_tune = WriteModifiedPhotoVMTune("forbidden_rho_rho_vertex", [](auto &j) {
    j["PARAM_REGGE"]["PARAM_CON"]["MP"]["[113,113]"].push_back({993, 9933});
  });
  REQUIRE_THROWS(gra::regge::ReadParam({113, 113}, pdg, *gra::MModelTune::Load(invalid_tune.second)));
}

TEST_CASE(
    "PARAM_CON MULTI uses library permutations with explicit "
    "continuum vertex filtering",
    "[gra::MRegge][continuum]") {
  const auto param = gra::regge::ReadParam({211, -211, 321, -321}, LoadedPDGTable(), *gra::MModelTune::Load(modelfile));

  struct Counts {
    std::size_t candidates = 0;
    std::size_t accepted   = 0;
    std::size_t rejected   = 0;
  };

  auto charge_class_order = [](const std::vector<int> &final_pdgs) {
    std::map<int, bool> has_antiparticle;
    for (const int pdg : final_pdgs) {
      if (pdg < 0) { has_antiparticle[std::abs(pdg)] = true; }
    }

    std::map<int, std::size_t> neutral_seen;
    std::vector<int>           plus;
    std::vector<int>           minus;
    for (std::size_t i = 0; i < final_pdgs.size(); ++i) {
      const int pdg     = final_pdgs[i];
      const int abs_pdg = std::abs(pdg);
      int       sign    = 0;
      if (has_antiparticle[abs_pdg]) { sign = (pdg >= 0) ? 1 : -1; }
      if (sign == 0) {
        const std::size_t seen = neutral_seen[pdg]++;
        sign                   = (seen % 2 == 0) ? 1 : -1;
      }
      if (sign > 0) {
        plus.push_back(static_cast<int>(i) + 3);
      } else {
        minus.push_back(static_cast<int>(i) + 3);
      }
    }

    REQUIRE(plus.size() == minus.size());
    std::vector<int> order;
    for (std::size_t i = 0; i < plus.size(); ++i) {
      order.push_back(plus[i]);
      order.push_back(minus[i]);
    }
    return order;
  };

  auto count = [&charge_class_order, &param](const std::vector<int> &final_pdgs, int mode) {
    Counts out;
    auto   permutations = gra::math::GetAmpPerm(4, mode);
    if (mode == 0) {
      const auto order = charge_class_order(final_pdgs);
      for (auto &permutation : permutations) {
        for (auto &index : permutation) { index = order.at(static_cast<std::size_t>(index - 3)); }
      }
    }
    out.candidates = permutations.size();
    for (const auto &permutation : permutations) {
      REQUIRE(permutation.size() == 4);
      const std::vector<int> first_pair  = {final_pdgs.at(static_cast<std::size_t>(permutation[0] - 3)),
                                            final_pdgs.at(static_cast<std::size_t>(permutation[1] - 3))};
      const std::vector<int> second_pair = {final_pdgs.at(static_cast<std::size_t>(permutation[2] - 3)),
                                            final_pdgs.at(static_cast<std::size_t>(permutation[3] - 3))};
      const bool             explicit_vertices =
          gra::regge::FindPair(param, first_pair, gra::ReggeProductionModel::MP) != nullptr &&
          gra::regge::FindPair(param, second_pair, gra::ReggeProductionModel::MP) != nullptr;
      if (explicit_vertices) {
        ++out.accepted;
      } else {
        ++out.rejected;
      }
    }
    return out;
  };

  const std::vector<std::vector<int>> final_state_orders = {
      {211, -211, 321, -321},
      {321, 211, -321, -211},
      {2212, 211, -211, -2212},
  };
  for (const auto &final_pdgs : final_state_orders) {
    CAPTURE(final_pdgs);
    const Counts charge_class     = count(final_pdgs, 0);
    const Counts all_combinations = count(final_pdgs, 1);

    CHECK(charge_class.candidates == 16);
    CHECK(all_combinations.candidates == 24);
    CHECK(charge_class.accepted == 8);
    CHECK(all_combinations.accepted == 8);
    CHECK(charge_class.rejected == 8);
    CHECK(all_combinations.rejected == 16);
  }
}

TEST_CASE("MRegge rotating eta mode", "[gra::MRegge]") {
  const auto         tune           = WriteModifiedPhotoVMTune("eta_mode", [](auto &j) {
    const std::string model                                      = j["PARAM_SOFT"]["active_model"];
    j["PARAM_SOFT"]["MODEL"][model]["EXCHANGE"]["P"]["eta_mode"] = "rotating";
  });
  const std::string &temp_modelfile = tune.second;

  const auto  model_tune = gra::MModelTune::Load(temp_modelfile);
  const auto  param      = gra::regge::ReadParam({211, -211}, LoadedPDGTable(), *model_tune);
  const auto  pomeron_id = param.exchanges.at(param.pomeron_trajectory).soft_exchange;
  const auto &pomeron    = model_tune->Soft()->Exchange(pomeron_id);
  REQUIRE(pomeron.Eta() == gra::EtaMode::Rotating);

  gra::LORENTZSCALAR lts;
  lts.PDG = LoadedPDGTable();
  MRegge regge(lts, model_tune, gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));

  const double               s     = 25.0;
  const double               t     = -0.1;
  const double               alpha = model_tune->Soft()->Alpha(pomeron_id, t);
  const std::complex<double> expected =
      gra::regge::EtaPhase(alpha, pomeron.Signature()) * std::pow(s / param.s0, alpha);
  const std::complex<double> actual = regge.PomeronKernel(s, t);

  REQUIRE(std::real(actual) == Approx(std::real(expected)).epsilon(1e-12));
  REQUIRE(std::imag(actual) == Approx(std::imag(expected)).epsilon(1e-12));
}

TEST_CASE("MRegge raw eta mode", "[gra::MRegge][ReggeSW]") {
  const auto         tune           = WriteModifiedPhotoVMTune("eta_mode_raw", [](auto &j) {
    const std::string model                                      = j["PARAM_SOFT"]["active_model"];
    j["PARAM_SOFT"]["MODEL"][model]["EXCHANGE"]["P"]["eta_mode"] = "raw";
  });
  const std::string &temp_modelfile = tune.second;

  const auto  model_tune = gra::MModelTune::Load(temp_modelfile);
  const auto  param      = gra::regge::ReadParam({211, -211}, LoadedPDGTable(), *model_tune);
  const auto  pomeron_id = param.exchanges.at(param.pomeron_trajectory).soft_exchange;
  const auto &pomeron    = model_tune->Soft()->Exchange(pomeron_id);
  REQUIRE(pomeron.Eta() == gra::EtaMode::Raw);

  gra::LORENTZSCALAR lts;
  lts.PDG = LoadedPDGTable();
  MRegge regge(lts, model_tune, gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));

  const double               s     = 25.0;
  const double               t     = -0.1;
  const double               alpha = model_tune->Soft()->Alpha(pomeron_id, t);
  const std::complex<double> expected =
      gra::regge::Pole(s, param.s0, {alpha, 0.0}, pomeron.Signature(), gra::regge::Rim::Lower);
  const std::complex<double> actual = regge.PomeronKernel(s, t);

  REQUIRE(std::real(actual) == Approx(std::real(expected)).epsilon(1e-12));
  REQUIRE(std::imag(actual) == Approx(std::imag(expected)).epsilon(1e-12));
}

// Check GP resonance and continuum amplitudes use both configured outer kernels
TEST_CASE("GP resonance and continuum use raw SW outer kernels", "[gra::MRegge][GP][ReggeSW][physics]") {
  ModelParamRestoreGuard restore;
  // Write one otherwise identical tune with the selected Pomeron eta mode
  const auto write_tune = [](const std::string &suffix, const std::string &mode) {
    return WriteModifiedPhotoVMTune(suffix, [&mode](auto &card) {
      auto             &soft                                           = card.at("PARAM_SOFT");
      const std::string model                                          = soft.at("active_model");
      soft.at("MODEL").at(model).at("EXCHANGE").at("P").at("eta_mode") = mode;
    });
  };
  const auto raw_tune       = write_tune("gp_outer_sw_raw", "raw");
  const auto rotating_tune  = write_tune("gp_outer_sw_rotating", "rotating");
  const auto raw_model      = gra::MModelTune::Load(raw_tune.second);
  const auto rotating_model = gra::MModelTune::Load(rotating_tune.second);

  gra::LORENTZSCALAR base     = MakeToyReggeLTSAsymmetric(0.31, 4.2, -4.6);
  base.process.MMAX           = 2;
  base.process.FORWARD_NOFLIP = true;
  SetToyContinuumExchangePair(base, 990, 990);
  UseToyGPSubchannelHelicity(base);

  struct Amplitudes {
    std::vector<std::complex<double>> resonance;
    std::vector<std::complex<double>> continuum_t;
    std::vector<std::complex<double>> continuum_u;
  };

  // Evaluate the physical GP amplitudes and resolve the two Born graphs
  const auto evaluate = [&base](const gra::MModelTunePtr &model_tune) {
    gra::LORENTZSCALAR event = base;
    event.model_cache.reset();
    gra::MRegge regge(event, model_tune, gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "gp_outer_sw"));
    gra::LORENTZSCALAR resonance_lts = event;
    gra::PARAM_RES     resonance     = MakeToyScalarGPResonance(resonance_lts.process.MMAX);
    TestReggeRes(regge, resonance_lts, resonance, gra::ReggeProductionModel::GP);

    const auto continuum = [&event, &regge](const double sign) {
      gra::LORENTZSCALAR lts   = event;
      lts.process.CONT_TU_SIGN = {sign};
      TestReggeCon(regge, lts, gra::ReggeProductionModel::GP);
      return lts.hamp;
    };
    const auto plus  = continuum(1.0);
    const auto minus = continuum(-1.0);
    REQUIRE(plus.size() == minus.size());
    Amplitudes out{resonance_lts.hamp, plus, plus};
    for (const auto &i : gra::aux::indices(plus)) {
      out.continuum_t[i] = 0.5 * (plus[i] + minus[i]);
      out.continuum_u[i] = 0.5 * (plus[i] - minus[i]);
    }
    return out;
  };

  const Amplitudes raw      = evaluate(raw_model);
  const Amplitudes rotating = evaluate(rotating_model);
  REQUIRE(gra::SquaredNorm(raw.resonance) > 0.0);
  REQUIRE(gra::SquaredNorm(raw.continuum_t) > 0.0);
  REQUIRE(gra::SquaredNorm(raw.continuum_u) > 0.0);

  gra::MRegge       raw_regge(base, raw_model,
                              gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "gp_outer_sw_reference"));
  const auto        param         = ReggeParametersForTest(raw_regge, base);
  const std::size_t trajectory    = gra::regge::TrajectoryIndex(*param, 990);
  const auto        exchange      = param->exchanges.at(trajectory).soft_exchange;
  const auto       &raw_soft      = raw_model->Soft()->Exchange(exchange);
  const auto       &rotating_soft = rotating_model->Soft()->Exchange(exchange);
  REQUIRE(raw_soft.Eta() == gra::EtaMode::Raw);
  REQUIRE(rotating_soft.Eta() == gra::EtaMode::Rotating);

  // Compute the exact raw over rotating kernel ratio of one outer leg
  const auto leg_ratio = [&](const double s, const double t) {
    const double raw_alpha      = raw_model->Soft()->Alpha(exchange, t);
    const double rotating_alpha = rotating_model->Soft()->Alpha(exchange, t);
    REQUIRE(raw_alpha == Approx(rotating_alpha).epsilon(2.0e-14));
    const std::complex<double> raw_kernel =
        gra::regge::Pole(s, param->s0, {raw_alpha, 0.0}, raw_soft.Signature(), gra::regge::Rim::Lower);
    const std::complex<double> rotating_kernel =
        gra::regge::EtaPhase(rotating_alpha, rotating_soft.Signature()) * std::pow(s / param->s0, rotating_alpha);
    return raw_kernel / rotating_kernel;
  };

  // Check one complete amplitude vector against its two outer kernels
  const auto require_outer_ratio = [](const auto &raw_amplitude, const auto &rotating_amplitude,
                                      const std::complex<double> ratio) {
    REQUIRE(raw_amplitude.size() == rotating_amplitude.size());
    auto expected = rotating_amplitude;
    for (auto &value : expected) { value *= ratio; }
    RequireVectorNear(raw_amplitude, expected, 3.0e-11);
  };

  require_outer_ratio(raw.resonance, rotating.resonance, leg_ratio(base.s1, base.t1) * leg_ratio(base.s2, base.t2));
  require_outer_ratio(raw.continuum_t, rotating.continuum_t,
                      leg_ratio(base.ss[1][3], base.t1) * leg_ratio(base.ss[2][4], base.t2));
  require_outer_ratio(raw.continuum_u, rotating.continuum_u,
                      leg_ratio(base.ss[1][4], base.t1) * leg_ratio(base.ss[2][3], base.t2));
}

TEST_CASE("MRegge keeps GP continuum entries separate from beam residues", "[gra::MRegge][params][physics]") {
  const auto tune = WriteModifiedPhotoVMTune("regge_proton_residue_ownership", [](auto &j) { (void)j; });

  const auto param = gra::regge::ReadParam({2212, -2212}, LoadedPDGTable(), *gra::MModelTune::Load(tune.second));

  const std::vector<int> central_ppbar = {2212, -2212};
  const auto             entry         = gra::regge::Pair(param, central_ppbar, gra::ReggeProductionModel::GP);
  const auto channel = std::find_if(entry.channels.begin(), entry.channels.end(), [](const auto &candidate) {
    return candidate.first == 990 && candidate.second == 990;
  });
  REQUIRE(channel != entry.channels.end());

  gra::LORENTZSCALAR        lts = MakeToyCoherentPhotonLTS();
  gra::MRegge               regge(lts, gra::MModelTune::Load(tune.second),
                                  gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));
  const gra::SoftExchangeId pomeron  = param.exchanges.at(param.pomeron_trajectory).soft_exchange;
  const auto               &soft     = *regge.SoftModelHandle();
  const auto                baseline = gra::SoftModel::LoadFromJson(modelfile, gra::aux::GetInputData(modelfile));
  const double              expected = baseline->PhysicalResidue(baseline->ExchangeId("P"), lts.t1);
  CHECK(soft.PhysicalResidue(pomeron, lts.t1) == Approx(expected).epsilon(1.0e-12));

  const auto  &photo                    = gra::regge::Photo(param, 443);
  const double photo_s                  = gra::math::pow2(photo.W0);
  const double physical_proton_coupling = soft.PhysicalCoupling(pomeron);
  CHECK(std::abs(soft.PhysicalResidue(pomeron, 0.0) * regge.PhotoKernel(photo_s, 0.0, 443)) ==
        Approx(photo_s * physical_proton_coupling).epsilon(1.0e-12));
}

TEST_CASE("MRegge central excitation keeps the full forward mass dependence", "[gra::MRegge][GoodWalker][excitation]") {
  gra::LORENTZSCALAR base           = MakeToyReggeLTSAsymmetric(0.3, 4.2, -4.6);
  const double       beam_energy    = 50.0;
  const double       beam_pz        = std::sqrt(gra::math::pow2(beam_energy) - gra::math::pow2(gra::PDG::mp));
  const double       forward_energy = 25.0;
  base.pbeam1                       = gra::M4Vec(0.0, 0.0, beam_pz, beam_energy);
  base.pbeam2                       = gra::M4Vec(0.0, 0.0, -beam_pz, beam_energy);
  const auto physical_forward       = [forward_energy](const gra::M4Vec &p, const double sign) {
    const double pz = std::sqrt(gra::math::pow2(forward_energy) - p.Pt2() - gra::math::pow2(gra::PDG::mp));
    return gra::M4Vec(p.Px(), p.Py(), sign * pz, forward_energy);
  };
  base.pfinal[1] = physical_forward(base.pfinal[1], 1.0);
  base.pfinal[2] = physical_forward(base.pfinal[2], -1.0);

  const auto rebuild_central_state = [](gra::LORENTZSCALAR &lts) {
    const gra::M4Vec central = lts.pbeam1 + lts.pbeam2 - lts.pfinal[1] - lts.pfinal[2];
    REQUIRE(central.M() > lts.decaytree[0].p.mass + lts.decaytree[1].p.mass);
    const auto momenta =
        TwoBodyRestKinematics(central.M(), lts.decaytree[0].p.mass, lts.decaytree[1].p.mass, 0.61, -0.37);
    lts.decaytree[0].p4 = BoostFromRestFrame(momenta[0], central);
    lts.decaytree[1].p4 = BoostFromRestFrame(momenta[1], central);
    UpdateToyDurhamDerivedKinematics(lts);
    const gra::M4Vec imbalance = lts.pbeam1 + lts.pbeam2 - lts.pfinal[0] - lts.pfinal[1] - lts.pfinal[2];
    REQUIRE(std::abs(imbalance.Px()) < 2.0e-12);
    REQUIRE(std::abs(imbalance.Py()) < 2.0e-12);
    REQUIRE(std::abs(imbalance.Pz()) < 2.0e-12);
    REQUIRE(std::abs(imbalance.E()) < 2.0e-12);
  };
  rebuild_central_state(base);
  const auto  tune = WriteModifiedPhotoVMTune("regge_forward_excitation_n2",
                                              [](auto &j) { j.at("PARAM_SOFT").at("active_model") = "double"; });
  gra::MRegge regge(base, gra::MModelTune::Load(tune.second),
                    gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));

  const auto source_profile = [&regge, &base, &rebuild_central_state](const double mass) {
    gra::LORENTZSCALAR lts = base;
    SetToyPhotoForwardExcitation(lts, 1, mass);
    rebuild_central_state(lts);
    const auto upper_state = gra::ResolveForwardLegState(lts, gra::ForwardBeamLeg::Upper);
    REQUIRE(upper_state.mass2 == Approx(gra::math::pow2(mass)).epsilon(2.0e-12));
    REQUIRE(upper_state.xi == Approx(lts.xi1).epsilon(1.0e-14));
    REQUIRE(upper_state.t == Approx(lts.q1.M2()).epsilon(1.0e-14));

    auto resonance                = MakeToyScalarMPResonance();
    resonance.MP.ff_prod          = {};
    resonance.hel_decay.ff_decay          = {};
    TestReggeRes(regge, lts, resonance, gra::ReggeProductionModel::MP);
    REQUIRE(lts.proton_good_walker.has_value());
    const auto &state = *lts.proton_good_walker;
    REQUIRE(state.model == regge.SoftModelHandle());
    REQUIRE(state.channel_count == regge.SoftModelHandle()->GoodWalker().ChannelCount());
    long double                       strength      = 0.0L;
    bool                              has_resolved  = false;
    bool                              has_inclusive = false;
    std::vector<std::complex<double>> projected;
    const auto                       &space = state.model->GoodWalker();
    for (const auto &component : state.components) {
      has_resolved |= component.upper_sector == gra::ProtonGoodWalkerSector::TripleResolved;
      has_inclusive |= component.upper_sector == gra::ProtonGoodWalkerSector::TripleInclusive;
      REQUIRE(component.lower_sector == gra::ProtonGoodWalkerSector::Elastic);
      REQUIRE(component.source.size_col() == space.PairDimension());
      const auto upper_basis = component.upper_sector == gra::ProtonGoodWalkerSector::TripleResolved
                                   ? gra::GoodWalkerFinalBasis::Excited
                                   : gra::GoodWalkerFinalBasis::Complete;
      for (std::size_t row = 0; row < component.source.size_row(); ++row) {
        for (std::size_t col = 0; col < component.source.size_col(); ++col) {
          const auto value = component.source(row, col);
          strength += static_cast<long double>(value.real()) * value.real() +
                      static_cast<long double>(value.imag()) * value.imag();
        }
        const std::vector<std::complex<double>> pair_source(component.source[row],
                                                            component.source[row] + component.source.size_col());
        const auto values = space.ProjectPair(pair_source, upper_basis, gra::GoodWalkerFinalBasis::Proton);
        projected.insert(projected.end(), values.begin(), values.end());
      }
    }
    REQUIRE(has_resolved);
    REQUIRE(has_inclusive);
    RequireVectorNear(projected, lts.hamp, 1.0e-12);

    const auto                        param   = ReggeParametersForTest(regge, lts);
    const auto                        pomeron = param->exchanges.at(param->pomeron_trajectory).soft_exchange;
    const auto                       &soft    = *regge.SoftModelHandle();
    std::vector<std::complex<double>> proton(space.ProtonVector().begin(), space.ProtonVector().end());
    const auto                        triple_source  = regge.TriplePomeronCouplingRoot() * proton;
    const auto                        elastic_source = soft.ResidueMatrix(pomeron, lts.t2) * proton;
    long double                       triple_norm    = 0.0L;
    long double                       elastic_norm   = 0.0L;
    for (const auto value : triple_source) {
      triple_norm +=
          static_cast<long double>(value.real()) * value.real() + static_cast<long double>(value.imag()) * value.imag();
    }
    for (const auto value : elastic_source) {
      elastic_norm +=
          static_cast<long double>(value.real()) * value.real() + static_cast<long double>(value.imag()) * value.imag();
    }

    const auto production = gra::rspin::Resonance(
        lts, resonance,
        (resonance.production_model == gra::ReggeProductionModel::MP ? gra::mpom::Fusion : gra::xpom::Fusion),
        param->s0, resonance.UsesUnrestrictedSpinBasis() ? nullptr : &resonance.MP.filter);
    REQUIRE(production.size() == 1);
    gra::spin::DecayAmp(lts, resonance, lts.process.MP_FRAME);
    const std::complex<double> common =
        resonance.hel_decay.g_decay * gra::resonance::LineShape(lts.m2, resonance) *
        ResonanceFormFactorReference(lts, resonance, resonance.production[0].tree);
    const auto  central      = production[0].MultiplyScaled(resonance.decay_f, common);
    long double central_norm = 0.0L;
    for (std::size_t row = 0; row < central.size_row(); ++row) {
      for (std::size_t col = 0; col < central.size_col(); ++col) {
        const auto value = central(row, col);
        central_norm += static_cast<long double>(value.real()) * value.real() +
                        static_cast<long double>(value.imag()) * value.imag();
      }
    }
    const std::complex<double> common_regge = std::pow(param->s0 / lts.m2, param->omega.at(resonance.production_model)) *
                                              regge.PomeronKernel(lts.s1, lts.t1) * regge.PomeronKernel(lts.s2, lts.t2);
    const long double common_norm = static_cast<long double>(common_regge.real()) * common_regge.real() +
                                    static_cast<long double>(common_regge.imag()) * common_regge.imag();
    const long double without_profile = central_norm * common_norm * triple_norm * elastic_norm;
    REQUIRE(without_profile > 0.0L);
    const double profile = soft.ForwardExcitationFactor(pomeron, upper_state.t, upper_state.mass2);
    REQUIRE(profile > 0.0);
    CHECK(static_cast<double>(strength / without_profile) == Approx(gra::math::pow2(profile)).epsilon(2.0e-10));
    return profile;
  };

  const double                kinematic_mass = std::sqrt(gra::math::pow2(forward_energy) - base.pfinal[1].Pt2());
  const std::array<double, 6> masses = {gra::PDG::mp + 2.0 * 0.13957061, 2.0, 5.0, 10.0, 15.0, 0.70 * kinematic_mass};
  double                      previous_profile = std::numeric_limits<double>::quiet_NaN();
  for (const double mass : masses) {
    CAPTURE(mass, kinematic_mass);
    REQUIRE(mass < kinematic_mass);
    const double profile = source_profile(mass);
    REQUIRE(std::isfinite(profile));
    if (std::isfinite(previous_profile)) { REQUIRE(profile != Approx(previous_profile).epsilon(1.0e-8)); }
    previous_profile = profile;
  }
}

TEST_CASE("MRegge secondary trajectories use mapped SOFT exchanges", "[gra::MRegge][physics][Reggeon]") {
  gra::LORENTZSCALAR lts = MakeToyCoherentPhotonLTS();
  gra::MRegge        regge(lts, gra::MModelTune::Load(modelfile),
                           gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));
  const auto         param      = gra::regge::ReadParam({211, -211}, LoadedPDGTable(), *regge.ModelTuneHandle());
  const std::size_t  f2         = gra::regge::TrajectoryIndex(param, 9910);
  const std::size_t  rho        = gra::regge::TrajectoryIndex(param, 9930);
  const auto        &soft_model = *regge.SoftModelHandle();
  const auto         f2_id      = param.exchanges.at(f2).soft_exchange;
  const auto         rho_id     = param.exchanges.at(rho).soft_exchange;
  const auto        &soft_f2    = soft_model.Exchange(f2_id);
  const auto        &soft_rho   = soft_model.Exchange(rho_id);
  CHECK(soft_f2.Name() == "R_f2");
  CHECK(soft_rho.Name() == "R_rho");
  CHECK(soft_f2.Eta() == gra::EtaMode::Rotating);
  CHECK(soft_rho.Eta() == gra::EtaMode::Rotating);

  const std::complex<double> f2_eta           = regge.ExchangeKernel(1.0, 0.0, 9910);
  const std::complex<double> rho_eta          = regge.ExchangeKernel(1.0, 0.0, 9930);
  const std::complex<double> expected_f2_eta  = gra::regge::EtaPhase(soft_f2.Alpha0(), soft_f2.Signature());
  const std::complex<double> expected_rho_eta = gra::regge::EtaPhase(soft_rho.Alpha0(), soft_rho.Signature());
  CHECK(std::real(f2_eta) == Approx(std::real(expected_f2_eta)).epsilon(1e-12));
  CHECK(std::imag(f2_eta) == Approx(std::imag(expected_f2_eta)).epsilon(1e-12));
  CHECK(std::real(rho_eta) == Approx(std::real(expected_rho_eta)).epsilon(1e-12));
  CHECK(std::imag(rho_eta) == Approx(std::imag(expected_rho_eta)).epsilon(1e-12));

  const auto projected_residue = [&lts, &soft_model](const gra::SoftExchangeId exchange) {
    const auto  matrix = soft_model.ResidueMatrix(exchange, lts.t1);
    const auto &proton = soft_model.GoodWalker().ProtonVector();
    double      value  = 0.0;
    for (std::size_t i = 0; i < proton.size(); ++i) {
      for (std::size_t j = 0; j < proton.size(); ++j) { value += proton[i] * matrix(i, j) * proton[j]; }
    }
    return value;
  };
  const double expected_f2_vertex  = projected_residue(f2_id);
  const double expected_rho_vertex = projected_residue(rho_id);
  CHECK(regge.SoftModelHandle()->PhysicalResidue(f2_id, lts.t1) == Approx(expected_f2_vertex).epsilon(1e-12));
  CHECK(regge.SoftModelHandle()->PhysicalResidue(rho_id, lts.t1) == Approx(expected_rho_vertex).epsilon(1e-12));
}

// Keep fixed and analytic PDG aliases grouped under the same SOFT trajectory
TEST_CASE("PARAM_REGGE groups fixed and analytic aliases by trajectory", "[gra::MRegge][params][Reggeon]") {
  gra::MPDG  pdg        = LoadedPDGTable();
  const auto model_tune = gra::MModelTune::Load(modelfile);

  const auto        param    = gra::regge::ReadParam({211, -211}, pdg, *model_tune);
  const std::size_t fixed    = gra::regge::TrajectoryIndex(param, 9925);
  const std::size_t analytic = gra::regge::TrajectoryIndex(param, 9920);
  CHECK(pdg.FindByPDG(9925).spinX2 / 2 == 2);
  CHECK(pdg.FindByPDG(9920).spinX2 == gra::aux::kNullSpinX2);
  CHECK(param.exchanges.at(fixed).soft_exchange == param.exchanges.at(analytic).soft_exchange);
}

// Read separate phase switches and reject incomplete model tables
TEST_CASE("PARAM_REGGE reads model specific decay phase steering", "[gra::MRegge][params][phase]") {
  const auto enabled_tune = WriteModifiedPhotoVMTune("regge_with_zeta", [](auto &j) {
    j.at("PARAM_REGGE").at("use_zeta") = {{"MP", true}, {"XP", true}, {"GP", true}};
  });
  const auto enabled_param =
      gra::regge::ReadParam({211, -211}, LoadedPDGTable(), *gra::MModelTune::Load(enabled_tune.second));
  CHECK(enabled_param.use_zeta.at(gra::ReggeProductionModel::GP));
  CHECK(enabled_param.use_zeta.at(gra::ReggeProductionModel::MP));
  CHECK(enabled_param.use_zeta.at(gra::ReggeProductionModel::XP));

  const auto disabled_tune = WriteModifiedPhotoVMTune(
      "regge_without_zeta", [](auto &j) {
        j.at("PARAM_REGGE").at("use_zeta") = {{"MP", true}, {"XP", true}, {"GP", false}};
      });
  const auto disabled_param = gra::regge::ReadParam(
      {211, -211}, LoadedPDGTable(), *gra::MModelTune::Load(disabled_tune.second));
  CHECK_FALSE(disabled_param.use_zeta.at(gra::ReggeProductionModel::GP));
  CHECK(disabled_param.use_zeta.at(gra::ReggeProductionModel::MP));
  CHECK(disabled_param.use_zeta.at(gra::ReggeProductionModel::XP));

  const auto missing_tune = WriteModifiedPhotoVMTune(
      "regge_missing_use_zeta", [](auto &j) { j.at("PARAM_REGGE").erase("use_zeta"); });
  REQUIRE_THROWS(
      gra::regge::ReadParam({211, -211}, LoadedPDGTable(), *gra::MModelTune::Load(missing_tune.second)));
  for (const std::string key : {"use_zeta", "omega"}) {
    const auto incomplete = WriteModifiedPhotoVMTune("regge_incomplete_" + key, [&](auto &j) {
      j.at("PARAM_REGGE").at(key).erase("XP");
    });
    REQUIRE_THROWS(gra::regge::ReadParam({211, -211}, LoadedPDGTable(), *gra::MModelTune::Load(incomplete.second)));
    const auto scalar = WriteModifiedPhotoVMTune("regge_scalar_" + key, [&](auto &j) {
      j.at("PARAM_REGGE").at(key) = 0.0;
    });
    REQUIRE_THROWS(gra::regge::ReadParam({211, -211}, LoadedPDGTable(), *gra::MModelTune::Load(scalar.second)));
  }
}

TEST_CASE("PARAM_REGGE rejects invalid trajectory steering", "[gra::MRegge][params][validation]") {
  SECTION("typed exchange mapping is required") {
    const auto tune =
        WriteModifiedPhotoVMTune("regge_missing_soft_mapping", [](auto &j) { j.at("PARAM_REGGE").erase("EXCHANGES"); });
    REQUIRE_THROWS(gra::regge::ReadParam({211, -211}, LoadedPDGTable(), *gra::MModelTune::Load(tune.second)));
  }


  SECTION("typed exchange mapping must not be empty") {
    const auto tune = WriteModifiedPhotoVMTune(
        "regge_short_soft_mapping", [](auto &j) { j.at("PARAM_REGGE").at("EXCHANGES") = nlohmann::json::array(); });
    REQUIRE_THROWS(gra::regge::ReadParam({211, -211}, LoadedPDGTable(), *gra::MModelTune::Load(tune.second)));
  }

  SECTION("one SOFT exchange cannot be mapped by multiple rows") {
    const auto tune = WriteModifiedPhotoVMTune("regge_duplicate_soft_mapping", [](auto &j) {
      auto &rows                     = j.at("PARAM_REGGE").at("EXCHANGES");
      rows.at(1).at("soft_exchange") = rows.at(0).at("soft_exchange");
    });
    REQUIRE_THROWS(gra::regge::ReadParam({211, -211}, LoadedPDGTable(), *gra::MModelTune::Load(tune.second)));
  }

  SECTION("unknown soft exchange") {
    const auto tune = WriteModifiedPhotoVMTune("regge_unknown_soft_exchange", [](auto &j) {
      j.at("PARAM_REGGE").at("EXCHANGES").at(0).at("soft_exchange") = "missing";
    });
    REQUIRE_THROWS(gra::regge::ReadParam({211, -211}, LoadedPDGTable(), *gra::MModelTune::Load(tune.second)));
  }

  SECTION("declared role must match the SOFT exchange") {
    const auto tune = WriteModifiedPhotoVMTune(
        "regge_wrong_exchange_role", [](auto &j) { j.at("PARAM_REGGE").at("EXCHANGES").at(1).at("role") = "odderon"; });
    REQUIRE_THROWS(gra::regge::ReadParam({211, -211}, LoadedPDGTable(), *gra::MModelTune::Load(tune.second)));
  }

  SECTION("disabled soft exchange") {
    const auto tune = WriteModifiedPhotoVMTune("regge_disabled_soft_exchange", [](auto &j) {
      const std::string model = j.at("PARAM_SOFT").at("active_model");
      j.at("PARAM_SOFT").at("MODEL").at(model).at("EXCHANGE").at("R_f2").at("on") = false;
    });
    REQUIRE_THROWS(gra::regge::ReadParam({211, -211}, LoadedPDGTable(), *gra::MModelTune::Load(tune.second)));
  }

  SECTION("soft quantum numbers must match aliases") {
    const auto tune = WriteModifiedPhotoVMTune("regge_wrong_soft_crossing", [](auto &j) {
      j.at("PARAM_SOFT").at("EXCHANGE_DEF").at("R_rho").at("crossing") = 1;
    });
    REQUIRE_THROWS(gra::regge::ReadParam({211, -211}, LoadedPDGTable(), *gra::MModelTune::Load(tune.second)));
  }

  SECTION("GP ppbar rows are not beam residues") {
    const auto tune = WriteModifiedPhotoVMTune("regge_no_gp_beam_residue", [](auto &j) { (void)j; });
    REQUIRE_NOTHROW(gra::regge::ReadParam({211, -211}, LoadedPDGTable(), *gra::MModelTune::Load(tune.second)));
  }

  SECTION("trajectory scale") {
    const auto tune = WriteModifiedPhotoVMTune("regge_bad_s0", [](auto &j) { j.at("PARAM_REGGE").at("s0") = 0.0; });
    REQUIRE_THROWS(gra::regge::ReadParam({211, -211}, LoadedPDGTable(), *gra::MModelTune::Load(tune.second)));
  }

  SECTION("photoproduction accepts only rotating modes") {
    const auto tune = WriteModifiedPhotoVMTune("regge_invalid_photo_eta",
                                               [](auto &j) { j.at("PARAM_REGGE").at("photoprod_eta_mode") = "gamma"; });
    REQUIRE_THROWS(gra::regge::ReadParam({211, -211}, LoadedPDGTable(), *gra::MModelTune::Load(tune.second)));
  }
}

TEST_CASE("MRegge Pomeron slot uses the nonlinear SOFT trajectory", "[gra::MRegge][Form]") {
  const std::string data                                       = gra::aux::GetInputData(modelfile);
  auto              j                                          = nlohmann::json::parse(data);
  const std::string model                                      = j["PARAM_SOFT"]["active_model"];
  j["PARAM_SOFT"]["MODEL"][model]["EXCHANGE"]["P"]["eta_mode"] = "raw";

  const auto         tune = WriteModifiedPhotoVMTune("regge_linear_pomeron", [&j](auto &general) { general = j; });
  const std::string &temp_modelfile = tune.second;

  const auto model_tune = gra::MModelTune::Load(temp_modelfile);
  const auto param      = gra::regge::ReadParam({211, -211}, LoadedPDGTable(), *model_tune);

  gra::LORENTZSCALAR lts;
  lts.PDG = LoadedPDGTable();
  MRegge regge(lts, model_tune, gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));

  const double s          = 25.0;
  const double t          = -0.1;
  const auto   pomeron_id = param.exchanges.at(param.pomeron_trajectory).soft_exchange;
  const auto  &pomeron    = model_tune->Soft()->Exchange(pomeron_id);
  const double alpha      = model_tune->Soft()->Alpha(pomeron_id, t);
  REQUIRE(alpha != Approx(pomeron.Alpha0() + pomeron.AlphaPrime() * t));
  const std::complex<double> expected =
      gra::regge::Pole(s, param.s0, {alpha, 0.0}, pomeron.Signature(), gra::regge::Rim::Lower);
  const std::complex<double> actual = regge.ExchangeKernel(s, t, 991);

  REQUIRE(std::isfinite(std::real(actual)));
  REQUIRE(std::isfinite(std::imag(actual)));
  REQUIRE(std::real(actual) == Approx(std::real(expected)).epsilon(1e-12));
  REQUIRE(std::imag(actual) == Approx(std::imag(expected)).epsilon(1e-12));
}

TEST_CASE("Physical Pomeron residues are invariant under eigenchannel relabeling",
          "[gra::MForm][gra::MRegge][physics]") {
  auto j                          = nlohmann::json::parse(gra::aux::GetInputData(modelfile));
  j["PARAM_SOFT"]["active_model"] = "double";
  auto &soft                      = j["PARAM_SOFT"]["MODEL"]["double"];
  soft["GW"]["theta"]             = {0.37};
  soft["EXCHANGE"]["P"]["g"]      = {{4.0, 0.0}, {0.0, 9.0}};
  soft["FF"]["P"]["param"]        = {{0.42, 1.31, 1.17}, {1.73, 2.24, 1.61}};
  soft["3P"]["g"]                 = 0.2;

  // Evaluate every physical vertex consumer for one eigenchannel labeling
  const auto evaluate = [](const nlohmann::json &card, const std::string &suffix) {
    const auto tune =
        WriteModifiedPhotoVMTune("soft_eigenchannel_permutation_" + suffix, [&card](auto &general) { general = card; });
    const std::string &path = tune.second;

    const double       t   = -0.23;
    gra::LORENTZSCALAR lts = MakeToyCoherentPhotonLTS();
    lts.t                  = t;
    lts.t1                 = t;
    lts.ss[1][1]           = 4.0;
    lts.ss[2][2]           = 9.0;
    gra::MRegge               regge(lts, gra::MModelTune::Load(path),
                                    gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));
    const auto                param   = gra::regge::ReadParam({211, -211}, LoadedPDGTable(), *regge.ModelTuneHandle());
    const auto               &photo   = gra::regge::Photo(param, 443);
    const double              photo_s = 1.25 * gra::math::pow2(photo.W0);
    const auto               &soft_model = *regge.SoftModelHandle();
    const gra::SoftExchangeId pomeron    = soft_model.ExchangeId("P");
    const auto photo_amplitude           = soft_model.PhysicalResidue(pomeron, t) * regge.PhotoKernel(photo_s, t, 443);
    lts.excite1                          = true;
    lts.excite2                          = false;
    const auto sd_amplitude              = regge.ME2(lts, gra::MReggeInclusive::SD);
    lts.excite2                          = true;
    const auto dd_amplitude              = regge.ME2(lts, gra::MReggeInclusive::DD);

    return std::array<double, 12>{soft_model.PhysicalCoupling(pomeron),
                                  soft_model.PhysicalResidue(pomeron, t),
                                  soft_model.PhysicalResidue(pomeron, t) / soft_model.PhysicalCoupling(pomeron),
                                  soft_model.Alpha(pomeron, t),
                                  soft_model.EffectiveTriplePomeronCoupling(pomeron),
                                  soft_model.ForwardExcitationFactor(pomeron, t, lts.ss[1][1]),
                                  std::real(photo_amplitude),
                                  std::imag(photo_amplitude),
                                  std::real(sd_amplitude),
                                  std::imag(sd_amplitude),
                                  std::real(dd_amplitude),
                                  std::imag(dd_amplitude)};
  };

  const auto reference = evaluate(j, "reference");

  const double theta     = soft["GW"]["theta"][0];
  soft["GW"]["theta"][0] = gra::math::PI / 2.0 - theta;
  for (auto it = soft["EXCHANGE"].begin(); it != soft["EXCHANGE"].end(); ++it) {
    auto      &coupling = it.value()["g"];
    const auto first    = coupling[0][0];
    coupling[0][0]      = coupling[1][1];
    coupling[1][1]      = first;
  }
  for (auto it = soft["FF"].begin(); it != soft["FF"].end(); ++it) {
    auto      &parameters = it.value()["param"];
    const auto first      = parameters[0];
    parameters[0]         = parameters[1];
    parameters[1]         = first;
  }

  const auto permuted = evaluate(j, "permuted");
  for (std::size_t i = 0; i < reference.size(); ++i) {
    CAPTURE(i, reference[i], permuted[i]);
    REQUIRE(permuted[i] == Approx(reference[i]).epsilon(1e-12).margin(1e-12));
  }
}

TEST_CASE("MRegge SD and DD use triple-Regge powers without S3FINEL", "[gra::MRegge][GoodWalker][triple-Regge]") {
  gra::LORENTZSCALAR lts = MakeToyCoherentPhotonLTS();
  lts.s                  = 100.0;
  lts.t                  = -0.17;
  lts.ss[1][1]           = 4.0;
  lts.ss[2][2]           = 9.0;
  gra::ScreeningMetadata metadata;
  metadata.spin_basis              = gra::ScreeningSpinBasis::ProtonIdentity;
  metadata.proton_mode             = gra::ProtonScreeningMode::TripleRegge;
  metadata.amplitude_normalization = 1.0;
  metadata.spin_rows               = 1;
  metadata.forward_noflip          = false;
  metadata.amplitude_type          = gra::ScreeningAmplitudeType::GoodWalker;
  metadata.PrepareSpinTransitions();
  lts.hamp.Configure(metadata);

  const auto  model = gra::MModelTune::Load(modelfile);
  gra::MRegge regge(lts, model, gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Soft, "triple_test"));
  const auto  parameters = ReggeParametersForTest(regge, lts);
  const auto  pomeron    = parameters->exchanges.at(parameters->pomeron_trajectory).soft_exchange;
  const auto &soft_model = *model->Soft();
  const auto &exchange   = soft_model.Exchange(pomeron);
  const auto &space      = soft_model.GoodWalker();
  std::vector<std::complex<double>> proton(space.ProtonVector().begin(), space.ProtonVector().end());
  const auto                        triple      = regge.TriplePomeronCouplingRoot() * proton;
  const auto                        elastic     = soft_model.ResidueMatrix(pomeron, lts.t) * proton;
  const auto                        vector_norm = [](const auto &source) {
    double out = 0.0;
    for (const auto value : source) { out += std::norm(value); }
    return out;
  };
  const double               triple_norm  = vector_norm(triple);
  const double               elastic_norm = vector_norm(elastic);
  const double               alpha_t      = soft_model.Alpha(pomeron, lts.t);
  const std::complex<double> eta =
      gra::regge::EtaFactor(alpha_t, exchange.Alpha0(), exchange.Signature(), soft_model.TriplePomeronEtaMode());

  for (const auto mode : {gra::MReggeInclusive::SD, gra::MReggeInclusive::DD}) {
    CAPTURE(mode);
    lts.excite1 = true;
    lts.excite2 = mode == gra::MReggeInclusive::DD;
    regge.ME2(lts, mode);
    REQUIRE(lts.proton_good_walker.has_value());
    double actual = 0.0;
    for (const auto &component : lts.proton_good_walker->components) {
      for (std::size_t row = 0; row < component.source.size_row(); ++row) {
        for (std::size_t col = 0; col < component.source.size_col(); ++col) {
          actual += std::norm(component.source(row, col));
        }
      }
    }

    const double mx2   = lts.ss[1][1];
    const double my2   = lts.ss[2][2];
    const double power = mode == gra::MReggeInclusive::SD
                             ? std::pow(lts.s / mx2, 2.0 * alpha_t) * std::pow(mx2 / parameters->s0, exchange.Alpha0())
                             : std::pow(lts.s * parameters->s0 / (mx2 * my2), 2.0 * alpha_t) *
                                   std::pow(mx2 / parameters->s0, exchange.Alpha0()) *
                                   std::pow(my2 / parameters->s0, exchange.Alpha0());
    const double source_norm =
        mode == gra::MReggeInclusive::SD ? triple_norm * elastic_norm : triple_norm * triple_norm;
    const double expected = std::norm(eta) * power * source_norm;
    CHECK(actual == Approx(expected).epsilon(2.0e-11));

    const double inelastic_profile =
        soft_model.ForwardExcitationFactor(pomeron, lts.t, mode == gra::MReggeInclusive::SD ? mx2 : my2);
    REQUIRE(std::abs(gra::math::pow2(inelastic_profile) - 1.0) > 1.0e-4);
    CHECK(actual != Approx(expected * gra::math::pow2(inelastic_profile)).epsilon(1e-6));
  }
}

TEST_CASE("MRegge SD and DD apply the mapped Pomeron sign exactly once",
          "[gra::MRegge][GoodWalker][triple-Regge][sign]") {
  const auto positive_tune  = WriteModifiedPhotoVMTune("positive_triple_regge_sign", [](auto &card) {
    const std::string active                                     = card["PARAM_SOFT"]["active_model"];
    card["PARAM_SOFT"]["MODEL"][active]["EXCHANGE"]["P"]["sign"] = 1;
  });
  const auto negative_tune  = WriteModifiedPhotoVMTune("negative_triple_regge_sign", [](auto &card) {
    const std::string active                                     = card["PARAM_SOFT"]["active_model"];
    card["PARAM_SOFT"]["MODEL"][active]["EXCHANGE"]["P"]["sign"] = -1;
  });
  const auto positive_model = gra::MModelTune::Load(positive_tune.second);
  const auto negative_model = gra::MModelTune::Load(negative_tune.second);
  REQUIRE(positive_model->Soft()->Exchange(positive_model->Soft()->ExchangeId("P")).ResidueSign() == 1);
  REQUIRE(negative_model->Soft()->Exchange(negative_model->Soft()->ExchangeId("P")).ResidueSign() == -1);

  const auto evaluate = [](const gra::MModelTunePtr &model, const gra::MReggeInclusive mode) {
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
    gra::MRegge                regge(lts, model, gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Soft, "sign_test"));
    const std::complex<double> amplitude = regge.ME2(lts, mode);
    REQUIRE(lts.proton_good_walker.has_value());
    return std::make_pair(amplitude, *lts.proton_good_walker);
  };

  for (const auto mode : {gra::MReggeInclusive::SD, gra::MReggeInclusive::DD}) {
    CAPTURE(mode);
    const auto positive = evaluate(positive_model, mode);
    const auto negative = evaluate(negative_model, mode);
    RequireComplexNear(negative.first, positive.first, 2.0e-12);
    REQUIRE(negative.second.components.size() == positive.second.components.size());
    for (std::size_t component = 0; component < positive.second.components.size(); ++component) {
      const auto &positive_source = positive.second.components[component].source;
      const auto &negative_source = negative.second.components[component].source;
      REQUIRE(negative_source.size_row() == positive_source.size_row());
      REQUIRE(negative_source.size_col() == positive_source.size_col());
      for (std::size_t row = 0; row < positive_source.size_row(); ++row) {
        for (std::size_t column = 0; column < positive_source.size_col(); ++column) {
          RequireComplexNear(negative_source(row, column), -positive_source(row, column), 2.0e-12);
        }
      }
    }
    CHECK(std::norm(negative.first) == Approx(std::norm(positive.first)).epsilon(2.0e-12));
  }
}

TEST_CASE(
    "gra::MRegge:: FORWARD_VERTEX modes produce distinct card-loaded "
    "helicity amplitudes",
    "[gra::MRegge][spin]") {
  ToyHelicityProcess proc;
  proc.SetProcessForTest("XP", "RES");
  proc.state.lts.PDG   = LoadedPDGTable();
  proc.state.lts.beam1 = proc.state.lts.PDG.FindByPDG(2212);
  proc.state.lts.beam2 = proc.state.lts.PDG.FindByPDG(2212);
  proc.SetDecayMode("pi+ pi-");

  gra::PARAM_RES f2 = gra::resonance::Read("RES/f2_1270.json", proc.state.random, gra::ReggeProductionModel::XP);
  proc.SetResonances({{"f2_1270", f2}});
  REQUIRE_NOTHROW(proc.InitializeProcessAmplitude());

  const auto resonances = proc.GetResonances();
  REQUIRE(resonances.count("f2_1270") == 1);
  gra::PARAM_RES regge = resonances.at("f2_1270");

  for (const bool FORWARD_NOFLIP : std::vector<bool>{true, false}) {
    CAPTURE(FORWARD_NOFLIP);
    gra::PARAM_RES pole = regge;

    gra::LORENTZSCALAR lts_regge     = MakeToyReggeLTSAsymmetric(0.3, 4.2, -4.6);
    lts_regge.process.MP_FRAME       = "HX";
    lts_regge.process.FORWARD_VERTEX = gra::ForwardVertexMode::HelicityResidue;
    lts_regge.process.FORWARD_NOFLIP = FORWARD_NOFLIP;
    lts_regge.PDG                    = LoadedPDGTable();
    gra::LORENTZSCALAR lts_pole      = lts_regge;
    lts_pole.process.FORWARD_VERTEX  = gra::ForwardVertexMode::UnitResidue;
    gra::MRegge regge_engine(lts_regge, gra::MModelTune::Load(modelfile),
                             gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));

    const std::complex<double> amp_regge = TestReggeRes(regge_engine, lts_regge, regge, gra::ReggeProductionModel::XP);
    const std::vector<std::complex<double>> hamp_regge = lts_regge.hamp;
    const std::complex<double> amp_pole = TestReggeRes(regge_engine, lts_pole, pole, gra::ReggeProductionModel::XP);
    const std::vector<std::complex<double>> hamp_pole = lts_pole.hamp;

    REQUIRE(hamp_regge.size() == hamp_pole.size());
    REQUIRE(regge.production.size() == 1);
    const auto        up_rows      = FORWARD_NOFLIP ? HelicityConservingRows(regge.production[0].tree[0].hel)
                                                    : AllHelicityRowsForTest(regge.production[0].tree[0].hel);
    const auto        dn_rows      = FORWARD_NOFLIP ? HelicityConservingRows(regge.production[0].tree[1].hel)
                                                    : AllHelicityRowsForTest(regge.production[0].tree[1].hel);
    const std::size_t initial_rows = up_rows.size() * dn_rows.size();
    REQUIRE(initial_rows > 0);
    REQUIRE(hamp_regge.size() % initial_rows == 0);
    const std::size_t decay_cols = hamp_regge.size() / initial_rows;

    double max_relative_diff = 0.0;
    for (std::size_t row = 0; row < initial_rows; ++row) {
      const std::size_t up_index = row / dn_rows.size();
      const std::size_t dn_index = row % dn_rows.size();
      REQUIRE(up_index < up_rows.size());
      REQUIRE(dn_index < dn_rows.size());

      for (std::size_t decay_col = 0; decay_col < decay_cols; ++decay_col) {
        const std::size_t index         = row * decay_cols + decay_col;
        const double      abs_diff      = std::abs(hamp_regge[index] - hamp_pole[index]);
        const double      scale         = std::max({std::abs(hamp_regge[index]), std::abs(hamp_pole[index]), 1e-300});
        const double      relative_diff = abs_diff / scale;
        max_relative_diff               = std::max(max_relative_diff, relative_diff);
      }
    }

    REQUIRE(std::abs(amp_regge - amp_pole) > 1e-12);
    REQUIRE(max_relative_diff > 1e-12);
  }
}

TEST_CASE("gra::spin:: vector-vector pseudoscalar XP has the forward dphi zero", "[gra::spin]") {
  for (const auto &frame : std::vector<std::string>{"HX", "CS", "CM"}) {
    CAPTURE(frame);
    gra::PARAM_RES res = MakeToyPseudoscalarVectorXP();

    gra::LORENTZSCALAR lts0      = MakeToyProductionLTSDphi(0.0);
    gra::LORENTZSCALAR lts90     = MakeToyProductionLTSDphi(0.5 * gra::math::PI);
    gra::LORENTZSCALAR lts180    = MakeToyProductionLTSDphi(gra::math::PI);
    gra::LORENTZSCALAR lts0_rot  = MakeToyProductionLTSDphi(0.0, 0.37);
    gra::LORENTZSCALAR lts90_rot = MakeToyProductionLTSDphi(0.5 * gra::math::PI, 0.37);
    lts0.process.MP_FRAME        = frame;
    lts90.process.MP_FRAME       = frame;
    lts180.process.MP_FRAME      = frame;
    lts0_rot.process.MP_FRAME    = frame;
    lts90_rot.process.MP_FRAME   = frame;

    const double amp2_0      = SteeredProductionAmp2(lts0, res);
    const double amp2_90     = SteeredProductionAmp2(lts90, res);
    const double amp2_180    = SteeredProductionAmp2(lts180, res);
    const double amp2_0_rot  = SteeredProductionAmp2(lts0_rot, res);
    const double amp2_90_rot = SteeredProductionAmp2(lts90_rot, res);

    CAPTURE(amp2_0, amp2_90, amp2_180, amp2_0_rot, amp2_90_rot);
    REQUIRE(amp2_90 > 1e-14);
    REQUIRE(amp2_0 < 1e-10 * amp2_90);
    REQUIRE(amp2_180 < 1e-10 * amp2_90);
    REQUIRE(amp2_90_rot == Approx(amp2_90).epsilon(1e-12));
    REQUIRE(amp2_0_rot < 1e-10 * amp2_90_rot);
  }
}

TEST_CASE("gra::spin:: vector-vector pseudoscalar MP keeps the forward dphi zero", "[gra::spin][MP]") {
  for (const auto &frame : std::vector<std::string>{"HX", "CS", "CM"}) {
    CAPTURE(frame);
    gra::PARAM_RES res = MakeToyPseudoscalarVectorRes();

    gra::LORENTZSCALAR lts0      = MakeToyProductionLTSDphi(0.0);
    gra::LORENTZSCALAR lts90     = MakeToyProductionLTSDphi(0.5 * gra::math::PI);
    gra::LORENTZSCALAR lts180    = MakeToyProductionLTSDphi(gra::math::PI);
    gra::LORENTZSCALAR lts0_rot  = MakeToyProductionLTSDphi(0.0, 0.37);
    gra::LORENTZSCALAR lts90_rot = MakeToyProductionLTSDphi(0.5 * gra::math::PI, 0.37);
    lts0.process.MP_FRAME        = frame;
    lts90.process.MP_FRAME       = frame;
    lts180.process.MP_FRAME      = frame;
    lts0_rot.process.MP_FRAME    = frame;
    lts90_rot.process.MP_FRAME   = frame;

    const double amp2_0      = SteeredProductionAmp2(lts0, res);
    const double amp2_90     = SteeredProductionAmp2(lts90, res);
    const double amp2_180    = SteeredProductionAmp2(lts180, res);
    const double amp2_0_rot  = SteeredProductionAmp2(lts0_rot, res);
    const double amp2_90_rot = SteeredProductionAmp2(lts90_rot, res);

    CAPTURE(amp2_0, amp2_90, amp2_180, amp2_0_rot, amp2_90_rot);
    REQUIRE(amp2_90 > 1e-14);
    REQUIRE(amp2_0 < 1e-10 * amp2_90);
    REQUIRE(amp2_180 < 1e-10 * amp2_90);
    REQUIRE(amp2_90_rot == Approx(amp2_90).epsilon(1e-12));
    REQUIRE(amp2_0_rot < 1e-10 * amp2_90_rot);
  }
}

TEST_CASE(
    "gra::spin:: vector-vector pseudoscalar XP keeps the full-basis "
    "dphi zero",
    "[gra::spin]") {
  for (const auto &frame : std::vector<std::string>{"HX", "CS", "CM"}) {
    CAPTURE(frame);
    gra::PARAM_RES res = MakeToyPseudoscalarVectorXP();

    gra::LORENTZSCALAR lts0       = MakeToyProductionLTSDphi(0.0);
    gra::LORENTZSCALAR lts90      = MakeToyProductionLTSDphi(0.5 * gra::math::PI);
    gra::LORENTZSCALAR lts180     = MakeToyProductionLTSDphi(gra::math::PI);
    lts0.process.MP_FRAME         = frame;
    lts0.process.FORWARD_NOFLIP   = false;
    lts90.process.MP_FRAME        = frame;
    lts90.process.FORWARD_NOFLIP  = false;
    lts180.process.MP_FRAME       = frame;
    lts180.process.FORWARD_NOFLIP = false;

    const double amp2_0   = SteeredProductionAmp2(lts0, res);
    const double amp2_90  = SteeredProductionAmp2(lts90, res);
    const double amp2_180 = SteeredProductionAmp2(lts180, res);

    CAPTURE(amp2_0, amp2_90, amp2_180);
    REQUIRE(amp2_90 > 1e-14);
    REQUIRE(amp2_0 < 1e-10 * amp2_90);
    REQUIRE(amp2_180 < 1e-10 * amp2_90);
  }
}

TEST_CASE("gra::spin:: tensor-tensor pseudoscalar XP keeps the forward dphi zero", "[gra::spin]") {
  for (const auto &frame : std::vector<std::string>{"HX", "CS", "CM"}) {
    CAPTURE(frame);
    gra::PARAM_RES res = MakeToyPseudoscalarTensorXP();

    gra::LORENTZSCALAR lts0      = MakeToyProductionLTSDphi(0.0);
    gra::LORENTZSCALAR lts90     = MakeToyProductionLTSDphi(0.5 * gra::math::PI);
    gra::LORENTZSCALAR lts180    = MakeToyProductionLTSDphi(gra::math::PI);
    gra::LORENTZSCALAR lts0_rot  = MakeToyProductionLTSDphi(0.0, 0.37);
    gra::LORENTZSCALAR lts90_rot = MakeToyProductionLTSDphi(0.5 * gra::math::PI, 0.37);
    lts0.process.MP_FRAME        = frame;
    lts90.process.MP_FRAME       = frame;
    lts180.process.MP_FRAME      = frame;
    lts0_rot.process.MP_FRAME    = frame;
    lts90_rot.process.MP_FRAME   = frame;

    const double amp2_0      = SteeredProductionAmp2(lts0, res);
    const double amp2_90     = SteeredProductionAmp2(lts90, res);
    const double amp2_180    = SteeredProductionAmp2(lts180, res);
    const double amp2_0_rot  = SteeredProductionAmp2(lts0_rot, res);
    const double amp2_90_rot = SteeredProductionAmp2(lts90_rot, res);

    CAPTURE(amp2_0, amp2_90, amp2_180, amp2_0_rot, amp2_90_rot);
    REQUIRE(amp2_90 > 1e-14);
    REQUIRE(amp2_0 < 1e-10 * amp2_90);
    REQUIRE(amp2_180 < 1e-10 * amp2_90);
    REQUIRE(amp2_90_rot == Approx(amp2_90).epsilon(1e-12));
    REQUIRE(amp2_0_rot < 1e-10 * amp2_90_rot);
  }
}

TEST_CASE(
    "gra::rspin:: MP Prod3 applies "
    "central MP tensor",
    "[gra::rspin]") {
  for (const auto &frame : std::vector<std::string>{"HX", "CS", "CM"}) {
    CAPTURE(frame);

    gra::LORENTZSCALAR lts = MakeToyProductionLTSAsymmetric(0.3, 4.2, -4.6);
    lts.process.MP_FRAME   = frame;
    gra::PARAM_RES res     = MakeToyMPResonance(false);

    const auto up_rows        = HelicityConservingRows(res.production[0].tree[0].hel);
    const auto dn_rows        = HelicityConservingRows(res.production[0].tree[1].hel);

    const auto f_up_first_full  = gra::spin::Forward(lts, res.production[0].tree[0], lts.pbeam1, lts.pfinal[1], false,
                                                     gra::spin::Rows(res.production[0].tree[0], false),
                                                     lts.process.PHOTON_VERTEX, gra::spin::ForwardSpec{});
    const auto f_dn_second_full = gra::spin::Forward(lts, res.production[0].tree[1], lts.pbeam2, lts.pfinal[2], true,
                                                     gra::spin::Rows(res.production[0].tree[1], false),
                                                     lts.process.PHOTON_VERTEX, gra::spin::ForwardSpec{});
    const auto f_in             = SelectRows(f_up_first_full, up_rows).Kronecker(SelectRows(f_dn_second_full, dn_rows));
    const auto f_x              = gra::mpom::Fusion(lts, res.production[0].pole.value());

    const auto expected = f_in * f_x;
    const auto actual   = gra::rspin::Resonance(
          lts, res,
          (res.production_model == gra::ReggeProductionModel::MP ? gra::mpom::Fusion : gra::xpom::Fusion), 1.0);

    REQUIRE(actual.size() == 1);
    REQUIRE(actual[0].size_row() == 4);
    REQUIRE(actual[0].size_col() == res.a_Jz.size());
    RequireMatrixNear(actual[0], expected);
  }
}

TEST_CASE("gra::gpom factorized analytic residues match the nested formulation", "[gra::MRegge][spin][equivalence]") {
  const auto   noflip_rows = gra::spin::SpinHalfTransitions(true);
  const auto   full_rows   = gra::spin::SpinHalfTransitions(false);
  const double s0          = 1.7;

  for (const bool second_exchange_daughter : {false, true}) {
    const gra::M4Vec q =
        second_exchange_daughter ? gra::M4Vec(-0.31, 0.27, -0.4, 0.9) : gra::M4Vec(0.37, -0.22, 0.5, 1.1);
    const double alpha_t = second_exchange_daughter ? 0.63 : 1.41;

    for (const bool drop_flip : {false, true}) {
      const auto &rows = drop_flip ? noflip_rows : full_rows;
      for (const bool use_exchange_helicity_barrier : {false, true}) {
        for (const int MMAX : {0, 1, 2, 3}) {
          CAPTURE(second_exchange_daughter, drop_flip, use_exchange_helicity_barrier, MMAX, alpha_t);

          std::vector<int>    m_values;
          std::vector<double> nonsense;
          for (int m = -MMAX; m <= MMAX; ++m) {
            m_values.push_back(m);
            double factor = 1.0;
            for (int k = 0; k < std::abs(m); ++k) { factor *= alpha_t - static_cast<double>(k); }
            nonsense.push_back(factor);
          }

          const std::size_t reference_column = static_cast<std::size_t>(MMAX);
          auto expected = NestedAnalyticReggeResidueForTest(rows, m_values, nonsense, q, s0, second_exchange_daughter,
                                                            use_exchange_helicity_barrier);
          const double collider_transfer_azimuth = (second_exchange_daughter ? q : -q).Phi();
          for (std::size_t row = 0; row < rows.size(); ++row) {
            expected.ScaleRow(
                row, gra::spin::SpinHalfForwardHelicitySectionFactor(
                         rows[row].first, rows[row].second, collider_transfer_azimuth, second_exchange_daughter));
          }
          const auto actual =
              gra::gpom::ResidueMatrix(rows, m_values, nonsense, q, s0, second_exchange_daughter,
                                       use_exchange_helicity_barrier, "test factorized analytic Regge residue");

          RequireMatrixNear(actual, expected, 1e-13);
          for (const auto &row : gra::aux::indices(rows)) {
            REQUIRE(std::abs(actual[row][reference_column]) == Approx(1.0).epsilon(1e-13));
          }
        }
      }
    }
  }
}

// Check the complete spin-half pair section including reverse flip signs
TEST_CASE("gra::gpom factorized residues obey collider helicity covariance",
          "[gra::MRegge][spin][helicity][covariance]") {
  const auto                rows     = gra::spin::SpinHalfTransitions(false);
  const std::vector<int>    m_values = {0};
  const std::vector<double> nonsense = {1.0};
  const gra::M4Vec          q1(-0.37, 0.0, 0.5, 1.1);
  const gra::M4Vec          q2(0.37, 0.0, -0.4, 0.9);
  const auto                destination_rows = gra::spin::CanonicalProtonPairSpinLayout::KroneckerDestinationRows(4, 4);

  // Build one scalar exchange pair at the common reference plane
  const auto pair_source = [&](gra::M4Vec upper_q, gra::M4Vec lower_q) {
    const auto upper = gra::gpom::ResidueMatrix(rows, m_values, nonsense, upper_q, 1.0, false, false,
                                                "test upper factorized collider section");
    const auto lower = gra::gpom::ResidueMatrix(rows, m_values, nonsense, lower_q, 1.0, true, false,
                                                "test lower factorized collider section");
    return upper.Kronecker(lower, destination_rows);
  };

  const auto     reference   = pair_source(q1, q2);
  constexpr auto helicity_x2 = gra::spin::BinaryHelicityLabelsX2();
  for (std::size_t i1 = 0; i1 < 2; ++i1) {
    for (std::size_t i2 = 0; i2 < 2; ++i2) {
      for (std::size_t f1 = 0; f1 < 2; ++f1) {
        for (std::size_t f2 = 0; f2 < 2; ++f2) {
          const std::size_t row           = gra::spin::CanonicalProtonPairSpinLayout::HardRow(i1, i2, f1, f2);
          const std::size_t transpose_row = gra::spin::CanonicalProtonPairSpinLayout::HardRow(f1, f2, i1, i2);
          const double      sign          = gra::spin::ColliderSpinHalfReciprocitySign(helicity_x2[i1], helicity_x2[i2],
                                                                                       helicity_x2[f1], helicity_x2[f2]);
          CAPTURE(i1, i2, f1, f2, row, reference[row][0], reference[transpose_row][0]);
          REQUIRE(std::abs(reference[row][0]) > 0.0);
          RequireComplexNear(reference[row][0], sign * reference[transpose_row][0], 1.0e-12);
        }
      }
    }
  }

  for (const double angle : {-1.21, 0.37, 2.04}) {
    gra::M4Vec rotated_q1 = q1;
    gra::M4Vec rotated_q2 = q2;
    rotated_q1.RotateZ(angle);
    rotated_q2.RotateZ(angle);
    const auto rotated = pair_source(rotated_q1, rotated_q2);
    for (std::size_t i1 = 0; i1 < 2; ++i1) {
      for (std::size_t i2 = 0; i2 < 2; ++i2) {
        for (std::size_t f1 = 0; f1 < 2; ++f1) {
          for (std::size_t f2 = 0; f2 < 2; ++f2) {
            const std::size_t row      = gra::spin::CanonicalProtonPairSpinLayout::HardRow(i1, i2, f1, f2);
            const int         harmonic = gra::spin::ColliderSpinHalfHelicityHarmonic(helicity_x2[i1], helicity_x2[i2],
                                                                                     helicity_x2[f1], helicity_x2[f2]);
            const std::complex<double> expected =
                std::exp(gra::math::zi * static_cast<double>(harmonic) * angle) * reference[row][0];
            CAPTURE(angle, i1, i2, f1, f2, harmonic, rotated[row][0], expected);
            RequireComplexNear(rotated[row][0], expected, 1.0e-12);
          }
        }
      }
    }
  }
}

TEST_CASE("gra::spin:: spin-disabled production uses an incoherent hidden-spin basis", "[gra::spin][MRegge]") {
  gra::LORENTZSCALAR lts = MakeToyReggeLTSAsymmetric(0.3, 4.2, -4.6);

  SECTION("X-mode flat production fallback") {
    gra::PARAM_RES res  = MakeToyCovariantXPonance();
    lts.process.SPINGEN = false;
    TestProd3Sum(lts, res, 1.0);
    REQUIRE(res.prod_f.size_row() == 12);
  }

  SECTION("MP flat production fallback") {
    gra::PARAM_RES res  = MakeToyMPResonance(false);
    lts.process.SPINGEN = false;
    TestProd3Sum(lts, res, 1.0);
    REQUIRE(res.prod_f.size_row() == 12);
  }

  SECTION("MP blind production does not reapply polarization") {
    gra::PARAM_RES res  = MakeToyMPResonance(false);
    lts.process.SPINGEN = false;
    const auto plain    = gra::rspin::Resonance(lts, res, gra::mpom::Fusion, 1.0, nullptr);
    const auto steered  = gra::rspin::Resonance(lts, res, gra::mpom::Fusion, 1.0, &res.MP.filter);
    REQUIRE(plain.size() == 1);
    REQUIRE(steered.size() == 1);
    RequireMatrixNear(plain.front(), steered.front());
  }

  SECTION(
      "MP flat production fallback keeps flip rows when "
      "requested") {
    gra::PARAM_RES res         = MakeToyMPResonance(false);
    lts.process.SPINGEN        = false;
    lts.process.FORWARD_NOFLIP = false;
    TestProd3Sum(lts, res, 1.0);
    REQUIRE(res.prod_f.size_row() == 48);
  }

  SECTION("continuum blind production keeps separate t and u residues") {
    gra::LORENTZSCALAR cont_lts = MakeToyContinuumLTSAsymmetric(0.3, 4.2, -4.6);
    cont_lts.process.SPINGEN    = false;
    auto &cache                 = cont_lts.process.CONTINUUM_POLE.front();
    auto  upper_t               = cache[0].Pole();
    auto  upper_u               = cache[2].Pole();
    upper_t.terms.front().coefficient *= 2.0;
    upper_u.terms.front().coefficient *= 3.0;
    cache[0]            = gra::spin::PoleResidue(upper_t);
    cache[2]            = gra::spin::PoleResidue(upper_u);
    const auto channels = gra::mpom::Continuum(cont_lts, 1.0);
    REQUIRE(channels.size() == 1);
    const auto       &tree           = cont_lts.process.CONT_PRODUCTIONTREE.front();
    const std::size_t initial_states = gra::spin::Rows(tree[0], cont_lts.process.FORWARD_NOFLIP).size() *
                                       gra::spin::Rows(tree[1], cont_lts.process.FORWARD_NOFLIP).size();
    const std::size_t final_states = gra::spin::FinalStateHelicityCount(cont_lts.decaytree[0].p, "test left final") *
                                     gra::spin::FinalStateHelicityCount(cont_lts.decaytree[1].p, "test right final");
    const double t_scale = std::sqrt(gra::spin::LeadingPoleDensity(cache[0].Pole(), cache[0].Pole().Lambda) *
                                     gra::spin::LeadingPoleDensity(cache[1].Pole(), cache[1].Pole().Lambda));
    const double u_scale = std::sqrt(gra::spin::LeadingPoleDensity(cache[2].Pole(), cache[2].Pole().Lambda) *
                                     gra::spin::LeadingPoleDensity(cache[3].Pole(), cache[3].Pole().Lambda));
    RequireMatrixNear(channels.front().first, gra::spin::Blind(initial_states, final_states, t_scale));
    RequireMatrixNear(channels.front().second, gra::spin::Blind(initial_states, final_states, u_scale));
  }

  SECTION("continuum flat production fallback keeps flip rows when requested") {
    gra::LORENTZSCALAR cont_lts     = MakeToyContinuumLTSAsymmetric(0.3, 4.2, -4.6);
    cont_lts.process.SPINGEN        = false;
    cont_lts.process.FORWARD_NOFLIP = false;
    const auto mats                 = TestProd4Sum(cont_lts, 1.0);
    REQUIRE(mats.first.size_row() == 144);
    REQUIRE(mats.second.size_row() == 144);
  }

  SECTION("spin-disabled covariant XP keeps its SOFT-reference scale") {
    gra::MRegge    regge(lts, gra::MModelTune::Load(modelfile),
                         gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));
    gra::PARAM_RES res  = MakeToyCovariantPhotoXPonance();
    lts.process.SPINGEN = false;

    const auto production = gra::rspin::Resonance(
        lts, res,
        (res.production_model == gra::ReggeProductionModel::MP ? gra::mpom::Fusion : gra::xpom::Fusion), 1.0);
    REQUIRE(production.size() == 1);
    const std::size_t initial_states =
        gra::spin::Rows(res.production.front().tree[0], lts.process.FORWARD_NOFLIP).size() *
        gra::spin::Rows(res.production.front().tree[1], lts.process.FORWARD_NOFLIP).size();
    const std::size_t spin_states = static_cast<std::size_t>(res.p.spinX2 + 1);
    const double      momentum    = gra::kinematics::PairRest(lts.q1_in_X, lts.q2_in_X).momentum;
    const double pole_scale = std::sqrt(gra::spin::LeadingPoleDensity(res.production.front().pole.value(), momentum));
    const std::complex<double> scale = pole_scale;
    RequireMatrixNear(production.front(), gra::spin::Blind(initial_states, spin_states, scale));
  }

  SECTION("GP reuses exact resonance spin coefficients within one amplitude") {
    gra::MRegge          regge(lts, gra::MModelTune::Load(modelfile),
                               gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));
    const auto           param    = ReggeParametersForTest(regge, lts);
    const gra::PARAM_RES res      = MakeToyScalarGPResonance(lts.process.MMAX);
    lts.process.DERIVATIVE_FACTOR = true;
    gra::gpom::AmpCache cache;
    const auto         &orbital = res.production.front().hel.gp_orbital;
    REQUIRE(orbital.nmu == 1);
    REQUIRE(orbital.spins == std::vector<std::size_t>{0});
    REQUIRE(orbital.terms.size() == 1);
    REQUIRE(orbital.su2.size() == 1);
    CHECK(orbital.su2.front() == Approx(std::sqrt(5.0)));

    const auto first = gra::gpom::Resonance(lts, *param, res, &cache);
    REQUIRE(first.size() == 1);
    REQUIRE(cache.sources.size() == 1);
    const std::size_t regge_count = cache.sources.front().regge_cg.size();
    REQUIRE(regge_count > 0);
    std::size_t value_count = 0;
    for (const auto &block : cache.sources.front().regge_cg) {
      value_count += static_cast<std::size_t>(std::count_if(block.coefficient.begin(), block.coefficient.end(),
                                                            [](const auto &value) { return value.has_value(); }));
    }
    REQUIRE(value_count > 0);

    const auto second = gra::gpom::Resonance(lts, *param, res, &cache);
    REQUIRE(second.size() == 1);
    RequireMatrixNear(second.front(), first.front());
    CHECK(cache.sources.front().regge_cg.size() == regge_count);
    std::size_t second_value_count = 0;
    for (const auto &block : cache.sources.front().regge_cg) {
      second_value_count += static_cast<std::size_t>(std::count_if(
          block.coefficient.begin(), block.coefficient.end(), [](const auto &value) { return value.has_value(); }));
    }
    CHECK(second_value_count == value_count);
  }

  SECTION("GP resonance derivative factor applies q^L to analytic LS rows") {
    gra::MRegge    regge(lts, gra::MModelTune::Load(modelfile),
                         gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));
    const auto     param = ReggeParametersForTest(regge, lts);
    gra::PARAM_RES res   = MakeToyTensorGPResonance(lts.process.MMAX);
    auto          &hel   = res.production.front().hel;
    hel.alpha_ls.Clear();
    hel.alpha_ls.Set(2, 4, {0.8, -0.2});
    hel.analytic_Lambda = 0.63;
    gra::gpom::InitResonanceLS(hel, res.p.spinX2 / 2, 4, 4);

    lts.process.DERIVATIVE_FACTOR = false;
    gra::gpom::AmpCache disabled_cache;
    const auto          disabled = gra::gpom::Resonance(lts, *param, res, &disabled_cache);
    REQUIRE(disabled.size() == 1);

    lts.process.DERIVATIVE_FACTOR = true;
    gra::gpom::AmpCache enabled_cache;
    const auto          enabled = gra::gpom::Resonance(lts, *param, res, &enabled_cache);
    REQUIRE(enabled.size() == 1);
    const double radial = pow2(lts.q1_in_X.P3mod() / hel.analytic_Lambda);
    CHECK(radial != Approx(1.0));
    RequireMatrixNear(enabled.front(), disabled.front() * radial);

    lts.process.SPINGEN           = false;
    lts.process.DERIVATIVE_FACTOR = false;
    gra::gpom::AmpCache blind_disabled_cache;
    const auto          blind_disabled = gra::gpom::Resonance(lts, *param, res, &blind_disabled_cache);

    lts.process.DERIVATIVE_FACTOR = true;
    gra::gpom::AmpCache blind_enabled_cache;
    const auto          blind_enabled = gra::gpom::Resonance(lts, *param, res, &blind_enabled_cache);
    RequireMatrixNear(blind_enabled.front(), blind_disabled.front() * radial);
  }

  SECTION("GP RES separates PP and RP forward sources") {
    lts.process.SPINGEN = true;
    lts.process.MMAX    = 1;
    gra::MRegge    regge(lts, gra::MModelTune::Load(modelfile),
                         gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));
    const auto     param   = ReggeParametersForTest(regge, lts);
    gra::PARAM_RES res     = MakeToyScalarGPResonance(lts.process.MMAX);
    auto           rp_tree = res.production.front().tree;
    rp_tree.front().p      = lts.PDG.FindByPDG(9910);
    res.production.push_back({std::move(rp_tree), res.production.front().hel});
    gra::gpom::AmpCache cache;

    const auto production = gra::gpom::Resonance(lts, *param, res, &cache);
    REQUIRE(production.size() == 2);
    REQUIRE(cache.sources.size() == 2);
    const auto pp = std::find_if(cache.sources.cbegin(), cache.sources.cend(), [](const auto &source) {
      return source.key.up_pdg == 990 && source.key.dn_pdg == 990;
    });
    const auto rp = std::find_if(cache.sources.cbegin(), cache.sources.cend(), [](const auto &source) {
      return source.key.up_pdg == 9910 && source.key.dn_pdg == 990;
    });
    CHECK(pp != cache.sources.cend());
    CHECK(rp != cache.sources.cend());
  }

  SECTION("GP RES rebuilds its orbital basis when initializing new LS couplings") {
    lts.process.MMAX = 1;
    lts.process.SPINGEN = true;
    gra::MRegge    regge(lts, gra::MModelTune::Load(modelfile),
                         gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));
    const auto     param = ReggeParametersForTest(regge, lts);
    gra::PARAM_RES res   = MakeToyScalarGPResonance(lts.process.MMAX);
    gra::gpom::AmpCache cache;
    const auto previous = gra::gpom::Resonance(lts, *param, res, &cache);
    REQUIRE(previous.size() == 1);

    auto &hel = res.production.front().hel;
    hel.alpha_ls.Clear();
    hel.alpha_ls.Set(2, 4, 1.0);
    // Initialize the changed couplings before evaluating another resonance
    gra::gpom::InitResonanceLS(hel, res.p.spinX2 / 2, 4, 4);
    const auto updated = gra::gpom::Resonance(lts, *param, res, &cache);
    const auto fresh = gra::gpom::Resonance(lts, *param, res, nullptr);
    REQUIRE(updated.size() == 1);
    REQUIRE(fresh.size() == 1);
    REQUIRE(fresh.front().FrobNorm2() > 0.0);
    RequireMatrixNear(updated.front(), fresh.front());
    REQUIRE((updated.front() - previous.front()).FrobNorm2() > 1.0e-8 * fresh.front().FrobNorm2());
  }

  SECTION("GP skips forbidden orbital projections before continuation") {
    lts.process.MMAX = 2;
    gra::MRegge    regge(lts, gra::MModelTune::Load(modelfile),
                         gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));
    const auto     param = ReggeParametersForTest(regge, lts);
    gra::PARAM_RES res   = MakeToyTensorGPResonance(lts.process.MMAX);
    res.production.front().hel.alpha_ls.Set(2, 0, 0.25);
    gra::gpom::InitResonanceLS(res.production.front().hel, res.p.spinX2 / 2, 4, 4);
    gra::gpom::AmpCache cache;

    const auto production = gra::gpom::Resonance(lts, *param, res, &cache);
    REQUIRE(production.size() == 1);
    REQUIRE(cache.sources.size() == 1);
    const auto &source = cache.sources.front();
    const auto  block  = std::find_if(source.regge_cg.cbegin(), source.regge_cg.cend(),
                                      [](const auto &candidate) { return candidate.two_s == 0; });
    REQUIRE(block != source.regge_cg.cend());
    const std::size_t i1          = gra::gpom::AnalyticMIndex(1, source.basis.MMAX, "test forbidden m1");
    const std::size_t i2          = gra::gpom::AnalyticMIndex(0, source.basis.MMAX, "test forbidden m2");
    const auto       &coefficient = block->coefficient[i1 * source.basis.nm + i2];
    CHECK_FALSE(coefficient.has_value());
  }

  SECTION("spin-disabled GP keeps the common physical pole scale") {
    gra::MRegge    regge(lts, gra::MModelTune::Load(modelfile),
                         gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));
    const auto     param = ReggeParametersForTest(regge, lts);
    gra::PARAM_RES res   = MakeToyScalarGPResonance(lts.process.MMAX);
    lts.process.SPINGEN  = false;
    gra::gpom::AmpCache cache;
    const auto          production = gra::gpom::Resonance(lts, *param, res, &cache);
    REQUIRE(production.size() == 1);
    CHECK(cache.sources.empty());
    const std::size_t proton_rows = gra::spin::SpinHalfTransitions(lts.process.FORWARD_NOFLIP).size();
    const std::size_t spin_states = static_cast<std::size_t>(res.p.spinX2 + 1);
    const auto       &hel         = res.production.front().hel;
    const double      scale       = std::sqrt(hel.alpha_ls.Norm2());
    RequireMatrixNear(production.front(), gra::spin::Blind(proton_rows * proton_rows, spin_states, scale));
  }
}

TEST_CASE("GP photon sources preserve the physical transverse density", "[gra::MRegge][spin][photon][normalization]") {
  for (const std::string mode : {"EPA", "QED"}) {
    DYNAMIC_SECTION(mode) {
      gra::LORENTZSCALAR lts     = MakeToyCoherentPhotonLTS();
      lts.process.FORWARD_NOFLIP = false;
      lts.process.MMAX           = 1;
      lts.process.PHOTON_VERTEX  = mode;
      gra::MRegge            regge(lts, gra::MModelTune::Load(modelfile),
                                   gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "photon density test"));
      const auto             param       = ReggeParametersForTest(regge, lts);
      const auto             transitions = gra::spin::SpinHalfTransitions(false);
      const std::vector<int> m_values    = {-1, 0, 1};
      const auto   upper_source          = gra::qed::PhotonSourceMatrixTransitions(lts, 1, transitions, m_values, mode);
      const auto   lower_source  = gra::qed::PhotonSourceMatrixTransitions(lts, 2, transitions, m_values, mode, true);
      const double upper_density = mode == "EPA" ? 1.0 : gra::qed::ElasticPhotonDensity(lts, 1);
      const double lower_density = mode == "EPA" ? 1.0 : gra::qed::ElasticPhotonDensity(lts, 2);
      CHECK(gra::spin::SourceSpinAveragedDensity(upper_source, 2, "test upper photon source") ==
            Approx(upper_density).epsilon(2.0e-12));
      CHECK(gra::spin::SourceSpinAveragedDensity(lower_source, 2, "test lower photon source") ==
            Approx(lower_density).epsilon(2.0e-12));

      gra::PARAM_RES      gamma_p = MakeToyPhotoGPResonance(lts.process.MMAX);
      gra::gpom::AmpCache gamma_p_cache;
      const auto          gamma_p_production = gra::gpom::Resonance(lts, *param, gamma_p, &gamma_p_cache);
      REQUIRE(gamma_p_production.size() == 1);
      REQUIRE(gamma_p_cache.sources.size() == 1);
      const auto &gamma_p_source = gamma_p_cache.sources.front();
      RequireMatrixNear(gamma_p_source.upper_residue, upper_source);
      CHECK(gamma_p_production.front().IsFinite());
      CHECK(gamma_p_production.front().FrobNorm2() > 0.0);

      gra::PARAM_RES       gamma_gamma         = MakeToyScalarGPResonance(lts.process.MMAX);
      const gra::MParticle photon              = lts.PDG.FindByPDG(gra::PDG::PDG_gamma);
      gamma_gamma.production.front().tree[0].p = photon;
      gamma_gamma.production.front().tree[1].p = photon;
      gra::gpom::AmpCache gamma_gamma_cache;
      const auto          gamma_gamma_production = gra::gpom::Resonance(lts, *param, gamma_gamma, &gamma_gamma_cache);
      REQUIRE(gamma_gamma_production.size() == 1);
      REQUIRE(gamma_gamma_cache.sources.size() == 1);
      const auto &gamma_gamma_source = gamma_gamma_cache.sources.front();
      RequireMatrixNear(gamma_gamma_source.upper_residue, upper_source);
      RequireMatrixNear(gamma_gamma_source.lower_residue, lower_source);
      CHECK(gamma_gamma_production.front().IsFinite());
      CHECK(gamma_gamma_production.front().FrobNorm2() > 0.0);
    }
  }
}

TEST_CASE("Spin-disabled MRegge QED rows retain the scalar EPA density", "[gra::MRegge][spin][photon][regression]") {
  const auto require_same_density = [](const std::string &label, const std::vector<std::complex<double>> &epa,
                                       const std::vector<std::complex<double>> &qed) {
    CAPTURE(label, epa.size(), qed.size());
    REQUIRE(qed.size() == epa.size());
    const double epa_density = gra::SquaredNorm(epa);
    const double qed_density = gra::SquaredNorm(qed);
    REQUIRE(std::isfinite(epa_density));
    REQUIRE(std::isfinite(qed_density));
    REQUIRE(epa_density > 0.0);
    CHECK(qed_density == Approx(epa_density).epsilon(2.0e-8));
  };

  const std::array<std::pair<gra::ReggeProductionModel, std::string>, 2> modes = {
      std::pair{gra::ReggeProductionModel::MP, std::string("M")},
      std::pair{gra::ReggeProductionModel::GP, std::string("G")}};

  for (const auto &[spin, label] : modes) {
    DYNAMIC_SECTION(label << " resonance") {
      gra::LORENTZSCALAR epa_lts    = MakeToyCoherentPhotonLTS();
      epa_lts.process.SPINGEN       = false;
      epa_lts.process.PHOTON_VERTEX = "EPA";
      epa_lts.process.MMAX          = 1;
      gra::LORENTZSCALAR qed_lts    = epa_lts;
      qed_lts.process.PHOTON_VERTEX = "QED";

      gra::PARAM_RES epa_res = spin == gra::ReggeProductionModel::GP
                                   ? MakeToyPhotoGPResonance(epa_lts.process.MMAX)
                                   : (spin == gra::ReggeProductionModel::XP ? MakeToyCovariantPhotoXPonance()
                                                                            : MakeToyPhotoMPResonance(false));
      gra::PARAM_RES qed_res = epa_res;
      gra::MRegge    epa_regge(epa_lts, gra::MModelTune::Load(modelfile),
                               gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "blind EPA resonance test"));
      gra::MRegge    qed_regge(qed_lts, gra::MModelTune::Load(modelfile),
                               gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "blind QED resonance test"));

      TestReggeRes(epa_regge, epa_lts, epa_res, spin);
      TestReggeRes(qed_regge, qed_lts, qed_res, spin);
      require_same_density(label + "RES", epa_lts.hamp, qed_lts.hamp);
    }

    DYNAMIC_SECTION(label << " continuum") {
      const auto tune = WriteModifiedPhotoVMTune("blind_photon_continuum_" + label, {}, {}, SetToyPhotonContinuum);
      gra::LORENTZSCALAR epa_lts    = MakeToyCoherentPhotonLTS();
      ConfigureToyPionPair(epa_lts);
      epa_lts.process.SPINGEN       = false;
      epa_lts.process.PHOTON_VERTEX = "EPA";
      epa_lts.process.MMAX          = 1;
      epa_lts.process.TU_SIGN       = "positive";
      SetToyContinuumExchangePair(epa_lts, 22, 22);
      epa_lts.process.CONT_PRODUCTIONTREE[0][0].p.pdg = 22;
      epa_lts.process.CONT_PRODUCTIONTREE[0][1].p.pdg = 22;
      if (spin == gra::ReggeProductionModel::GP) { UseToyGPSubchannelHelicity(epa_lts); }
      gra::LORENTZSCALAR qed_lts    = epa_lts;
      qed_lts.process.PHOTON_VERTEX = "QED";

      gra::MRegge epa_regge(epa_lts, gra::MModelTune::Load(tune.second),
                            gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "blind EPA continuum test"));
      gra::MRegge qed_regge(qed_lts, gra::MModelTune::Load(tune.second),
                            gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "blind QED continuum test"));

      TestReggeCon(epa_regge, epa_lts, spin);
      TestReggeCon(qed_regge, qed_lts, spin);
      require_same_density(label + "CON", epa_lts.hamp, qed_lts.hamp);
    }
  }
}

TEST_CASE("MRegge central production uses the configured Regge scale", "[gra::MRegge][spin][regression]") {
  const auto tune = WriteModifiedPhotoVMTune("regge_central_s0", [](auto &model) { model["PARAM_REGGE"]["s0"] = 4.0; });
  const std::string &model_path = tune.second;

  gra::LORENTZSCALAR lts = MakeToyReggeLTSAsymmetric(0.3, 4.2, -4.6);
  // Compare the scale in a common external helicity basis
  lts.process.MP_FRAME = "CM";
  gra::MRegge        regge(lts, gra::MModelTune::Load(model_path),
                           gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));
  gra::PARAM_RES     resonance = MakeToyMPResonance(false);
  TestReggeRes(regge, lts, resonance, gra::ReggeProductionModel::MP);
  const auto actual = lts.hamp;

  TestProd3Sum(lts, resonance, 4.0);
  gra::spin::DecayAmp(lts, resonance, lts.process.MP_FRAME);
  const auto production = gra::rspin::Resonance(
      lts, resonance,
      (resonance.production_model == gra::ReggeProductionModel::MP ? gra::mpom::Fusion : gra::xpom::Fusion),
      4.0, resonance.UsesUnrestrictedSpinBasis() ? nullptr : &resonance.MP.filter);
  REQUIRE(production.size() == 1);
  const auto expected = (production[0] * resonance.decay_f) * MPCommonFactor(regge, lts, resonance);
  RequireMatrixNear(ReshapeFlatMatrix(actual, 4, resonance.hel_decay.lambda_values.size_row()), expected);
}

// Check the central subenergy factor is absent from every photon production row
TEST_CASE("MRegge photon resonance production omits the central subenergy residue",
          "[gra::MRegge][photon][resonance][normalization]") {
  const auto flat_tune = WriteModifiedPhotoVMTune("photon_resonance_omega_zero", [](auto &model) {
    model["PARAM_REGGE"]["s0"]    = 4.0;
    model["PARAM_REGGE"]["omega"]["XP"] = 0.0;
  });
  const auto scaled_tune = WriteModifiedPhotoVMTune("photon_resonance_omega_one", [](auto &model) {
    model["PARAM_REGGE"]["s0"]    = 4.0;
    model["PARAM_REGGE"]["omega"]["XP"] = 1.0;
  });

  const auto evaluate = [](const std::string &tune, gra::PARAM_RES resonance) {
    gra::LORENTZSCALAR lts = MakeToyReggeLTSAsymmetric(0.3, 4.2, -4.6);
    gra::MRegge        regge(lts, gra::MModelTune::Load(tune),
                             gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));
    TestReggeRes(regge, lts, resonance, gra::ReggeProductionModel::XP);
    return std::pair{lts, lts.hamp};
  };

  gra::PARAM_RES photon_odderon = MakeToyCovariantPhotoXPonance();
  auto          &odderon        = photon_odderon.production.front().tree.back().p;
  odderon.pdg                    = 9993;
  odderon.P                      = -1;
  PrepareToyXPOperators(photon_odderon, {{{1, 2, 1.0}}});
  const auto photon_flat   = evaluate(flat_tune.second, photon_odderon);
  const auto photon_scaled = evaluate(scaled_tune.second, photon_odderon);
  RequireVectorNear(photon_scaled.second, photon_flat.second, 2.0e-12);

  const auto strong_flat   = evaluate(flat_tune.second, MakeToyCovariantXPonance());
  const auto strong_scaled = evaluate(scaled_tune.second, MakeToyCovariantXPonance());
  auto       strong_expected = strong_flat.second;
  for (auto &value : strong_expected) { value *= 4.0 / strong_flat.first.m2; }
  RequireVectorNear(strong_scaled.second, strong_expected, 2.0e-12);
}

TEST_CASE("MRegge amplitudes use the compact no-flip forward-proton basis", "[gra::MRegge][spin]") {
  gra::LORENTZSCALAR lts = MakeToyReggeLTSAsymmetric(0.3, 4.2, -4.6);
  // Compare row ordering in a common external helicity basis
  lts.process.MP_FRAME = "CM";
  lts.process.MMAX       = 2;
  gra::MRegge regge(lts, gra::MModelTune::Load(modelfile),
                    gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));

  const std::size_t continuum_dim = 9;

  TestReggeCon(regge, lts, gra::ReggeProductionModel::MP);
  const auto continuum_hamp = lts.hamp;
  REQUIRE(continuum_hamp.size() == 4 * continuum_dim);
  REQUIRE(std::isfinite(gra::SquaredNorm(continuum_hamp)));

  gra::PARAM_RES mp = MakeToyMPResonance(false);
  TestReggeRes(regge, lts, mp, gra::ReggeProductionModel::MP);
  const auto mp_hamp = lts.hamp;
  TestProd3Sum(lts, mp, 1.0);
  gra::spin::DecayAmp(lts, mp, lts.process.MP_FRAME);
  const auto mp_prod = gra::rspin::Resonance(
      lts, mp, (mp.production_model == gra::ReggeProductionModel::MP ? gra::mpom::Fusion : gra::xpom::Fusion),
      1.0, mp.UsesUnrestrictedSpinBasis() ? nullptr : &mp.MP.filter);
  REQUIRE(mp_prod.size() == 1);
  const auto mp_expected = (mp_prod[0] * mp.decay_f) * MPCommonFactor(regge, lts, mp);
  REQUIRE(mp_hamp.size() == 4 * mp.hel_decay.lambda_values.size_row());
  REQUIRE(std::isfinite(gra::SquaredNorm(mp_hamp)));
  REQUIRE(mp_hamp.size() == continuum_hamp.size());
  {
    INFO("MP coherent reference");
    RequireMatrixNear(ReshapeFlatMatrix(mp_hamp, 4, mp.hel_decay.lambda_values.size_row()), mp_expected);
  }

  gra::PARAM_RES rho_mp  = MakeToyMPResonance(false);
  rho_mp.spin_basis      = "rho";
  rho_mp.production.front().g = 1.0;
  rho_mp.MP.filter = gra::HelAmp::IdentityMatrix(3);
  TestReggeRes(regge, lts, rho_mp, gra::ReggeProductionModel::MP);
  REQUIRE(lts.hamp.size() == mp_hamp.size());
  REQUIRE(std::isfinite(gra::SquaredNorm(lts.hamp)));

  gra::PARAM_RES complex_rho  = MakeToyMPResonance(false);
  complex_rho.spin_basis      = "rho";
  complex_rho.production.front().g = 1.0;
  complex_rho.MP.filter = gra::RankOneProjector(gra::HelVec{{0.7, 0.2}, {-0.3, 0.5}, {0.4, -0.1}}).Transpose();
  TestReggeRes(regge, lts, complex_rho, gra::ReggeProductionModel::MP);
  gra::spin::DecayAmp(lts, complex_rho, lts.process.MP_FRAME);
  const auto reference_rho = gra::mpom::Resonance(lts, complex_rho, 1.0);
  const auto expected_rho = (reference_rho.front() * complex_rho.decay_f) * MPCommonFactor(regge, lts, complex_rho);
  {
    INFO("MP density reference");
    RequireMatrixNear(ReshapeFlatMatrix(lts.hamp, 4, complex_rho.hel_decay.lambda_values.size_row()), expected_rho,
                      2.0e-12);
  }

  gra::PARAM_RES xp_mode = MakeToyCovariantXPonance();
  TestReggeRes(regge, lts, xp_mode, gra::ReggeProductionModel::XP);
  const auto xp_mode_hamp = lts.hamp;
  TestProd3Sum(lts, xp_mode, 1.0);
  gra::spin::DecayAmp(lts, xp_mode, "CM");
  const auto xp_mode_expected = (xp_mode.prod_f * xp_mode.decay_f) * MPCommonFactor(regge, lts, xp_mode);
  REQUIRE(lts.hamp.size() == 4 * xp_mode.hel_decay.lambda_values.size_row());
  REQUIRE(std::isfinite(gra::SquaredNorm(lts.hamp)));
  REQUIRE(lts.hamp.size() == continuum_hamp.size());
  {
    INFO("XP hadronic reference");
    RequireMatrixNear(ReshapeFlatMatrix(xp_mode_hamp, 4, xp_mode.hel_decay.lambda_values.size_row()), xp_mode_expected);
  }

  gra::PARAM_RES photo = MakeToyPhotoMPResonance(false);
  TestReggeRes(regge, lts, photo, gra::ReggeProductionModel::MP);
  const auto photo_hamp = lts.hamp;
  gra::spin::DecayAmp(lts, photo, lts.process.MP_FRAME);
  const auto photo_prod = gra::rspin::Resonance(
      lts, photo,
      (photo.production_model == gra::ReggeProductionModel::MP ? gra::mpom::Fusion : gra::xpom::Fusion), 1.0,
      photo.UsesUnrestrictedSpinBasis() ? nullptr : &photo.MP.filter);
  const auto photo_expected = (photo_prod[0] * photo.decay_f) * PhotoCommonFactor(regge, lts, photo);
  REQUIRE(lts.hamp.size() == 4 * photo.hel_decay.lambda_values.size_row());
  REQUIRE(std::isfinite(gra::SquaredNorm(lts.hamp)));
  REQUIRE(lts.hamp.size() == continuum_hamp.size());
  REQUIRE(photo_prod.size() == 1);
  {
    INFO("MP photon reference");
    RequireMatrixNear(ReshapeFlatMatrix(photo_hamp, 4, photo.hel_decay.lambda_values.size_row()), photo_expected);
  }

  gra::PARAM_RES photox = MakeToyCovariantScalarPomeronPhotoXP();
  TestReggeRes(regge, lts, photox, gra::ReggeProductionModel::XP);
  const auto photox_hamp = lts.hamp;
  TestProd3Sum(lts, photox, 1.0);
  gra::spin::DecayAmp(lts, photox, "CM");
  const auto photox_prod = gra::rspin::Resonance(
      lts, photox,
      (photox.production_model == gra::ReggeProductionModel::MP ? gra::mpom::Fusion : gra::xpom::Fusion),
      1.0);
  const auto photox_expected = (photox_prod[0] * photox.decay_f) * PhotoCommonFactor(regge, lts, photox);
  REQUIRE(lts.hamp.size() == 4 * photox.hel_decay.lambda_values.size_row());
  REQUIRE(std::isfinite(gra::SquaredNorm(lts.hamp)));
  REQUIRE(photox_prod.size() == 1);
  {
    INFO("XP photon reference");
    RequireMatrixNear(ReshapeFlatMatrix(photox_hamp, 4, photox.hel_decay.lambda_values.size_row()), photox_expected);
  }

  gra::PARAM_RES photohel = MakeToyPhotoGPResonance(lts.process.MMAX);
  TestReggeRes(regge, lts, photohel, gra::ReggeProductionModel::GP);
  const auto photohel_uncached_hamp = lts.hamp;
  REQUIRE(lts.hamp.size() == 4 * photohel.hel_decay.lambda_values.size_row());
  REQUIRE(std::isfinite(gra::SquaredNorm(lts.hamp)));
  REQUIRE(FiniteHelicityAmp2(photohel_uncached_hamp) > 0.0);

  gra::PARAM_RES lower_photohel = MakeToyPhotoGPResonance(lts.process.MMAX);
  std::swap(lower_photohel.production[0].tree[0], lower_photohel.production[0].tree[1]);
  TestReggeRes(regge, lts, lower_photohel, gra::ReggeProductionModel::GP);
  const auto lower_photohel_uncached_hamp = lts.hamp;
  REQUIRE(lts.hamp.size() == 4 * lower_photohel.hel_decay.lambda_values.size_row());
  REQUIRE(std::isfinite(gra::SquaredNorm(lts.hamp)));
  REQUIRE(FiniteHelicityAmp2(lower_photohel_uncached_hamp) > 0.0);

  gra::LORENTZSCALAR ppbar_lts = lts;
  ppbar_lts.beam2.pdg          = -2212;
  ppbar_lts.beam2.chargeX3     = -3;
  const auto beam_param        = ReggeParametersForTest(regge, lts);
  REQUIRE(gra::regge::AntiparticleSign(*beam_param, 9993, ppbar_lts.beam2.pdg) == Approx(-1.0));
  REQUIRE(gra::regge::AntiparticleSign(*beam_param, 993, ppbar_lts.beam2.pdg) == Approx(1.0));

  gra::PARAM_RES upper_gamma_pp = MakeToyCovariantPhotoXPonance();
  TestReggeRes(regge, lts, upper_gamma_pp, gra::ReggeProductionModel::XP);
  const auto     upper_gamma_pp_hamp = lts.hamp;
  gra::PARAM_RES upper_gamma_ppbar   = MakeToyCovariantPhotoXPonance();
  TestReggeRes(regge, ppbar_lts, upper_gamma_ppbar, gra::ReggeProductionModel::XP);
  const auto upper_gamma_ppbar_hamp = ppbar_lts.hamp;
  RequireMatrixNear(ReshapeFlatMatrix(upper_gamma_ppbar_hamp, 4, upper_gamma_ppbar.hel_decay.lambda_values.size_row()),
                    ReshapeFlatMatrix(upper_gamma_pp_hamp, 4, upper_gamma_pp.hel_decay.lambda_values.size_row()));

  gra::PARAM_RES lower_gamma_pp = MakeToyCovariantPhotoXPonance();
  std::swap(lower_gamma_pp.production[0].tree[0], lower_gamma_pp.production[0].tree[1]);
  PrepareToyXPOperators(lower_gamma_pp, {{{0, 2, 1.0}}});
  TestReggeRes(regge, lts, lower_gamma_pp, gra::ReggeProductionModel::XP);
  const auto     lower_gamma_pp_hamp = lts.hamp;
  gra::PARAM_RES lower_gamma_ppbar   = lower_gamma_pp;
  TestReggeRes(regge, ppbar_lts, lower_gamma_ppbar, gra::ReggeProductionModel::XP);
  const auto lower_gamma_ppbar_hamp = ppbar_lts.hamp;
  RequireMatrixNear(ReshapeFlatMatrix(lower_gamma_ppbar_hamp, 4, lower_gamma_ppbar.hel_decay.lambda_values.size_row()),
                    ReshapeFlatMatrix(lower_gamma_pp_hamp, 4, lower_gamma_pp.hel_decay.lambda_values.size_row()) *
                        std::complex<double>(-1.0, 0.0));

  gra::PARAM_RES    mixed      = MakeToyCovariantPhotoXPonance();
  gra::MDecayBranch odderon_up = mixed.production[0].tree[1];
  odderon_up.p.pdg             = 9993;
  odderon_up.p.P               = -1;
  gra::MDecayBranch odderon_dn = odderon_up;
  gra::MDecayBranch pomeron_dn = mixed.production[0].tree[1];
  pomeron_dn.p.pdg             = 993;
  pomeron_dn.p.P               = -1;
  mixed.production.push_back({{odderon_up, pomeron_dn}});
  mixed.production.push_back({{mixed.production[0].tree[0], odderon_dn}});
  PrepareToyXPOperators(mixed, {{{0, 2, 1.0}}, {{1, 2, 1.0}}, {{1, 2, 1.0}}});
  TestReggeRes(regge, lts, mixed, gra::ReggeProductionModel::XP);
  const auto mixed_hamp = lts.hamp;
  gra::spin::DecayAmp(lts, mixed, "CM");
  const auto mixed_prod = gra::rspin::Resonance(
      lts, mixed,
      (mixed.production_model == gra::ReggeProductionModel::MP ? gra::mpom::Fusion : gra::xpom::Fusion), 1.0);
  const auto mixed_expected =
      (mixed_prod[0] * mixed.decay_f) * MixedProductionCommonFactor(regge, lts, mixed, mixed.production[0].tree) +
      (mixed_prod[1] * mixed.decay_f) * MixedProductionCommonFactor(regge, lts, mixed, mixed.production[1].tree) +
      (mixed_prod[2] * mixed.decay_f) * MixedProductionCommonFactor(regge, lts, mixed, mixed.production[2].tree);
  REQUIRE(mixed_prod.size() == 3);
  RequireMatrixNear(ReshapeFlatMatrix(mixed_hamp, 4, mixed.hel_decay.lambda_values.size_row()), mixed_expected);

  std::vector<std::complex<double>> coherent_sum = continuum_hamp;
  std::transform(mp_hamp.begin(), mp_hamp.end(), coherent_sum.begin(), coherent_sum.begin(), std::plus<>{});
  REQUIRE(std::isfinite(gra::SquaredNorm(coherent_sum)));
}

// Check photon production does not add the generic strong transfer factor
TEST_CASE("MRegge photon resonance routes disable strong transfer factors", "[gra::MRegge][photon][form-factor]") {
  gra::LORENTZSCALAR lts = MakeToyReggeLTSAsymmetric(0.3, 4.2, -4.6);
  lts.process.MMAX       = 2;
  gra::MRegge regge(lts, gra::MModelTune::Load(modelfile),
                    gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));

  for (const auto model :
       {gra::ReggeProductionModel::MP, gra::ReggeProductionModel::XP, gra::ReggeProductionModel::GP}) {
    CAPTURE(model);
    gra::PARAM_RES none            = model == gra::ReggeProductionModel::MP
                                         ? MakeToyPhotoMPResonance(false)
                                         : (model == gra::ReggeProductionModel::XP ? MakeToyCovariantScalarPomeronPhotoXP()
                                                                                   : MakeToyPhotoGPResonance(lts.process.MMAX));
    auto          &none_form       = ReggeResonanceFormForTest(none, model);
    none_form.ff_transfer          = {};
    TestReggeRes(regge, lts, none, model);
    const auto reference = lts.hamp;

    gra::PARAM_RES active      = none;
    auto          &active_form = ReggeResonanceFormForTest(active, model);
    active_form.ff_transfer    = {gra::regge::FFType::Power, gra::regge::FFNorm::Zero, {0.01, 4.0}};
    TestReggeRes(regge, lts, active, model);
    RequireVectorNear(lts.hamp, reference, 2.0e-12);

    // The production mass factor acts on every complex photon amplitude and is unity at the pole
    for (const bool at_pole : {false, true}) {
      if (at_pole) { none.p.mass = std::sqrt(lts.m2); }
      none_form.ff_prod = {};
      TestReggeRes(regge, lts, none, model);
      const auto bare = lts.hamp;
      REQUIRE(gra::SquaredNorm(bare) > 0.0);
      const double delta = lts.m2 - gra::math::pow2(none.p.mass);
      for (const double width : {0.8, 1.7}) {
        const double scale = width * (std::abs(delta) + 0.5);
        none_form.ff_prod = {gra::regge::FFType::Gaussian, gra::regge::FFNorm::Pole, {scale}};
        TestReggeRes(regge, lts, none, model);
        const double factor = std::exp(-gra::math::pow2(delta / scale));
        for (const auto &i : indices(bare)) { RequireComplexNear(lts.hamp[i], factor * bare[i], 2.0e-12); }
      }
    }
  }
}

TEST_CASE("MRegge photon vertex mode reaches resonance and continuum amplitudes", "[gra::MRegge][spin][photon]") {
  const auto tune = WriteModifiedPhotoVMTune("photon_vertex_mode", {}, {}, SetToyPhotonContinuum);
  gra::LORENTZSCALAR lts     = MakeToyCoherentPhotonLTS();
  lts.process.FORWARD_NOFLIP = false;
  lts.process.MMAX           = 1;
  gra::MRegge regge(lts, gra::MModelTune::Load(tune.second),
                    gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));

  auto require_nonzero = [](const std::string &label, const std::vector<std::complex<double>> &hamp) {
    const double amp2 = gra::SquaredNorm(hamp);
    CAPTURE(label, amp2);
    REQUIRE(std::isfinite(amp2));
    REQUIRE(amp2 > 0.0);
  };

  gra::PARAM_RES res = MakeToyPhotoMPResonance(false);
  TestReggeRes(regge, lts, res, gra::ReggeProductionModel::MP);
  const auto res_epa = lts.hamp;
  require_nonzero("MP EPA", res_epa);

  gra::PARAM_RES     res_qed        = MakeToyPhotoMPResonance(false);
  gra::LORENTZSCALAR res_qed_lts    = lts;
  res_qed_lts.process.PHOTON_VERTEX = "QED";
  TestReggeRes(regge, res_qed_lts, res_qed, gra::ReggeProductionModel::MP);
  const auto res_qed_hamp = res_qed_lts.hamp;
  require_nonzero("MP QED", res_qed_hamp);
  REQUIRE(MatrixDiffNorm2(ReshapeFlatMatrix(res_qed_hamp, 16, res_qed.hel_decay.lambda_values.size_row()),
                          ReshapeFlatMatrix(res_epa, 16, res.hel_decay.lambda_values.size_row())) > 1e-12);

  gra::PARAM_RES     xp_qed        = MakeToyCovariantPhotoXPonance();
  gra::LORENTZSCALAR xp_qed_lts    = lts;
  xp_qed_lts.process.PHOTON_VERTEX = "QED";
  TestReggeRes(regge, xp_qed_lts, xp_qed, gra::ReggeProductionModel::XP);
  require_nonzero("XP QED", xp_qed_lts.hamp);

  gra::PARAM_RES     gp_qed        = MakeToyPhotoGPResonance(lts.process.MMAX);
  gra::LORENTZSCALAR gp_qed_lts    = lts;
  gp_qed_lts.process.PHOTON_VERTEX = "QED";
  TestReggeRes(regge, gp_qed_lts, gp_qed, gra::ReggeProductionModel::GP);
  require_nonzero("GP QED", gp_qed_lts.hamp);

  gra::LORENTZSCALAR con_lts    = lts;
  ConfigureToyPionPair(con_lts);
  con_lts.process.PHOTON_VERTEX = "EPA";
  con_lts.process.TU_SIGN       = "positive";
  SetToyContinuumExchangePair(con_lts, 22, 22);
  con_lts.process.CONT_PRODUCTIONTREE[0][0].p.pdg = 22;
  con_lts.process.CONT_PRODUCTIONTREE[0][1].p.pdg = 22;
  gra::MRegge con_regge(con_lts, gra::MModelTune::Load(tune.second),
                        gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::ContinuumTwoBody, "photon continuum test"));
  const auto  con_epa_prod = gra::mpom::Continuum(con_lts, 1.0);
  REQUIRE(con_epa_prod.size() == 1);
  REQUIRE(MatrixNorm2(con_epa_prod[0].first) > 0.0);
  TestReggeCon(con_regge, con_lts, gra::ReggeProductionModel::MP);
  const auto con_epa = con_lts.hamp;
  REQUIRE(std::isfinite(gra::SquaredNorm(con_epa)));

  con_lts.process.PHOTON_VERTEX = "QED";
  const auto con_qed_prod       = gra::mpom::Continuum(con_lts, 1.0);
  REQUIRE(con_qed_prod.size() == 1);
  REQUIRE(MatrixNorm2(con_qed_prod[0].first) > 0.0);
  REQUIRE(MatrixDiffNorm2(con_qed_prod[0].first, con_epa_prod[0].first) > 1e-12);
  TestReggeCon(con_regge, con_lts, gra::ReggeProductionModel::MP);
  const auto con_qed = con_lts.hamp;
  REQUIRE(std::isfinite(gra::SquaredNorm(con_qed)));
  REQUIRE(con_qed.size() == con_epa.size());

  gra::LORENTZSCALAR gp_lts    = con_lts;
  gp_lts.process.PHOTON_VERTEX = "QED";
  SetToyContinuumExchangePair(gp_lts, 22, 22);
  gp_lts.process.CONT_PRODUCTIONTREE[0][0].p.pdg = 22;
  gp_lts.process.CONT_PRODUCTIONTREE[0][1].p.pdg = 22;
  UseToyGPSubchannelHelicity(gp_lts);
  TestReggeCon(con_regge, gp_lts, gra::ReggeProductionModel::GP);
  REQUIRE(std::isfinite(gra::SquaredNorm(gp_lts.hamp)));
}

// Match the complete physical photon sewing across every generic continuum
TEST_CASE("MP XP and GP gamma-gamma scalar continua share the pole matrix",
          "[gra::MRegge][spin][photon][continuum][normalization][regression]") {
  for (const std::string mode : {"EPA", "QED"}) {
    DYNAMIC_SECTION(mode) {
      gra::LORENTZSCALAR fixed = MakeToyScalarContinuumLTSAsymmetric(0.3, 4.2, -4.6);
      UpdateToyDerivedKinematics(fixed);
      fixed.model_cache            = std::make_shared<gra::MModelCache>(gra::MModelTune::Load(modelfile));
      fixed.process.FORWARD_NOFLIP = false;
      fixed.process.PHOTON_VERTEX  = mode;
      fixed.process.MMAX           = 1;
      SetToyContinuumExchangePair(fixed, gra::PDG::PDG_gamma, gra::PDG::PDG_gamma);

      gra::LORENTZSCALAR analytic = fixed;
      UseToyGPSubchannelHelicity(analytic);
      const std::array<std::complex<double>, 4> residue = {std::polar(0.83, 0.27), std::polar(1.11, -0.41),
                                                           std::polar(0.72, 0.63), std::polar(1.29, -0.19)};
      REQUIRE(fixed.process.CONTINUUM_POLE.size() == 1);
      REQUIRE(analytic.process.CONTINUUM_GP.size() == 1);
      for (const auto &vertex : indices(residue)) {
        auto  fixed_vertex    = fixed.process.CONTINUUM_POLE[0][vertex].Pole();
        auto &analytic_vertex = analytic.process.CONTINUUM_GP[0][vertex];
        REQUIRE(fixed_vertex.terms.size() == 1);
        REQUIRE(analytic_vertex.alpha_ls.Size() == 1);
        fixed_vertex.terms[0].coefficient             = residue[vertex];
        fixed.process.CONTINUUM_POLE[0][vertex]        = gra::spin::PoleResidue(fixed_vertex);
        analytic_vertex.alpha_ls.begin()->coefficient = residue[vertex];
      }

      gra::MRegge regge(
          fixed, gra::MModelTune::Load(modelfile),
          gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::ContinuumTwoBody, "matched gamma-gamma scalar continuum"));
      const auto param = ReggeParametersForTest(regge, fixed);
      const auto mp    = gra::mpom::Continuum(fixed, param->s0);
      const auto xp    = gra::xpom::Continuum(fixed, param->s0);
      const auto gp    = gra::gpom::Continuum(analytic, *param);

      REQUIRE(mp.size() == 1);
      REQUIRE(xp.size() == 1);
      REQUIRE(gp.size() == 1);
      REQUIRE(mp[0].first.FrobNorm2() > 0.0);
      REQUIRE(mp[0].second.FrobNorm2() > 0.0);
      RequireMatrixNear(xp[0].first, mp[0].first, 3.0e-12);
      RequireMatrixNear(xp[0].second, mp[0].second, 3.0e-12);
      RequireMatrixNear(gp[0].first, mp[0].first, 3.0e-12);
      RequireMatrixNear(gp[0].second, mp[0].second, 3.0e-12);

      // Spin averaging preserves the same physical photon coupling in every model
      auto blind_fixed = fixed;
      auto blind_analytic = analytic;
      blind_fixed.process.SPINGEN = false;
      blind_analytic.process.SPINGEN = false;
      const auto blind_gp = gra::gpom::Continuum(blind_analytic, *param);
      for (const auto model : {gra::ReggeProductionModel::MP, gra::ReggeProductionModel::XP}) {
        const auto blind = (model == gra::ReggeProductionModel::MP ? gra::mpom::Continuum : gra::xpom::Continuum)(
            blind_fixed, param->s0);
        RequireMatrixNear(blind_gp[0].first, blind[0].first, 3.0e-12);
        RequireMatrixNear(blind_gp[0].second, blind[0].second, 3.0e-12);
      }
      for (auto &hel : blind_analytic.process.CONTINUUM_GP[0]) {
        const auto reduced = gra::gpom::PhotonPoleHelicity(hel);
        hel.coupling_basis = gra::CouplingBasis::Helicity;
        hel.alpha_ls.Clear();
        for (std::size_t row = 0; row < hel.T.size_row(); ++row) {
          for (const std::size_t m : {0U, 2U}) {
            hel.T[row][m]     = reduced.T[hel.lambda_idx[row][0]][hel.lambda_idx[row][1]];
            hel.T_set[row][m] = true;
          }
        }
        gra::gpom::CheckHelicity(hel, "spin averaged photon pole", false);
      }
      const auto blind_direct = gra::gpom::Continuum(blind_analytic, *param);
      RequireMatrixNear(blind_direct[0].first, blind_gp[0].first, 3.0e-12);
      RequireMatrixNear(blind_direct[0].second, blind_gp[0].second, 3.0e-12);

      if (mode == "EPA") {
        const gra::M4Vec lower_final = gra::kinematics::BoostToRestFrame(fixed.decaytree[1].p4, fixed.pfinal[0],
                                                                         "test matched photon continuum lower final");
        gra::M4Vec       lower_axis  = fixed.q2_in_X;
        lower_axis.Flip3();
        const auto      &lower_vertex = fixed.process.CONTINUUM_POLE[0][1].Pole();
        const gra::M4Vec relative     = lower_final + lower_axis * 0.5;
        const auto       helicity     = gra::spin::PoleLSHelicity(lower_vertex, relative.P3mod());
        const auto       raw          = gra::spin::VirtualSubchannelFrame(helicity, relative, lower_axis, true);
        const auto       lower        = gra::spin::EvaluatePoleSubvertex(lower_vertex, lower_final, lower_axis, true);
        REQUIRE(raw.size_row() == 1);
        REQUIRE(raw.size_col() == 2);
        REQUIRE(std::abs(raw[0][0] - raw[0][1]) > 1.0e-6);
        RequireComplexNear(lower.frame[0][0], raw[0][1], 3.0e-13);
        RequireComplexNear(lower.frame[0][1], raw[0][0], 3.0e-13);
      }
    }
  }
}

TEST_CASE("Photon scalar continuum retains electromagnetic current tensors",
          "[gra::MRegge][spin][photon][continuum][physics][GP]") {
  gra::LORENTZSCALAR full     = MakeToyScalarContinuumLTSAsymmetric(0.3, 4.2, -4.6);
  full.model_cache            = std::make_shared<gra::MModelCache>(gra::MModelTune::Load(modelfile));
  full.process.FORWARD_NOFLIP = false;
  full.process.PHOTON_VERTEX  = "QED";
  full.process.MMAX           = 2;
  SetToyContinuumExchangePair(full, gra::PDG::PDG_gamma, gra::PDG::PDG_gamma);
  REQUIRE(full.process.SPINGEN);

  gra::LORENTZSCALAR blind = full;
  blind.process.SPINGEN    = false;
  for (const auto model : {gra::ReggeProductionModel::MP, gra::ReggeProductionModel::XP}) {
    CAPTURE(model);
    const auto full_prod =
        (model == gra::ReggeProductionModel::MP ? gra::mpom::Continuum : gra::xpom::Continuum)(full, 1.0);
    const auto blind_prod =
        (model == gra::ReggeProductionModel::MP ? gra::mpom::Continuum : gra::xpom::Continuum)(blind, 1.0);
    REQUIRE(full_prod.size() == 1);
    REQUIRE(blind_prod.size() == 1);
    REQUIRE(MatrixDiffNorm2(full_prod[0].first, blind_prod[0].first) > 1.0e-12);
  }

  UseToyGPSubchannelHelicity(full);
  blind                 = full;
  blind.process.SPINGEN = false;
  gra::MRegge regge(full, gra::MModelTune::Load(modelfile),
                    gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::ContinuumTwoBody, "photon_scalar_continuum"));
  const auto  param    = ReggeParametersForTest(regge, full);
  const auto  full_gp  = gra::gpom::Continuum(full, *param);
  const auto  blind_gp = gra::gpom::Continuum(blind, *param);
  REQUIRE(full_gp.size() == 1);
  REQUIRE(blind_gp.size() == 1);
  REQUIRE(MatrixDiffNorm2(full_gp[0].first, blind_gp[0].first) > 1.0e-12);

  for (const int mmax : {1, 4}) {
    gra::LORENTZSCALAR trial = full;
    trial.process.MMAX       = mmax;
    UseToyGPSubchannelHelicity(trial);
    for (const auto &vertex : trial.process.CONTINUUM_GP.front()) { CHECK(vertex.analytic_MMAX == 1); }
    const auto trial_production = gra::gpom::Continuum(trial, *param);
    RequireMatrixNear(trial_production[0].first, full_gp[0].first, 2.0e-13);
    RequireMatrixNear(trial_production[0].second, full_gp[0].second, 2.0e-13);
  }
}

TEST_CASE("MRegge GP[RES] plus GP[CON] matches the explicit coherent sum",
          "[gra::MRegge][spin][GP][continuum][resonance]") {
  gra::LORENTZSCALAR base = MakeToyReggeLTSAsymmetric(0.3, 4.2, -4.6);
  base.process.MMAX       = 2;
  base.process.TU_SIGN    = "positive";
  SetToyContinuumExchangePair(base, 990, 990);
  UseToyGPSubchannelHelicity(base);

  gra::PARAM_RES f0    = MakeToyGPResonance(base.process.MMAX);
  f0.hel_decay.g_decay = std::complex<double>(0.7, 0.1);
  gra::PARAM_RES f1    = MakeToyGPResonance(base.process.MMAX);
  f1.hel_decay.g_decay = std::complex<double>(0.25, -0.15);
  f1.production[0].hel.alpha_ls.Set(2, 2, std::complex<double>(0.35, 0.08));
  gra::gpom::InitResonanceLS(f1.production[0].hel, f1.p.spinX2 / 2, 4, 4);
  base.process.RESONANCES = {{"toy_f0", f0}, {"toy_f1", f1}};

  gra::LORENTZSCALAR explicit_lts = base;
  gra::MRegge        explicit_regge(explicit_lts, gra::MModelTune::Load(modelfile),
                                    gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));
  auto               explicit_resonances = std::move(explicit_lts.process.RESONANCES);
  TestReggeCon(explicit_regge, explicit_lts, gra::ReggeProductionModel::GP);
  std::vector<std::complex<double>> explicit_hamp = explicit_lts.hamp;
  for (auto &entry : explicit_resonances) {
    TestReggeRes(explicit_regge, explicit_lts, entry.second, gra::ReggeProductionModel::GP);
    REQUIRE(explicit_lts.hamp.size() == explicit_hamp.size());
    std::transform(explicit_lts.hamp.begin(), explicit_lts.hamp.end(), explicit_hamp.begin(), explicit_hamp.begin(),
                   std::plus<>{});
  }

  gra::LORENTZSCALAR cached_lts = base;
  gra::MRegge        cached_regge(cached_lts, gra::MModelTune::Load(modelfile),
                                  gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));
  cached_regge.Amp2(cached_lts, gra::ReggeProductionModel::GP, gra::MReggeMode::ResonanceContinuumTwoBody);

  REQUIRE(cached_lts.hamp.size() == explicit_hamp.size());
  RequireVectorNear(cached_lts.hamp, explicit_hamp, 1.0e-12);
}

TEST_CASE("MRegge composes finite-spin continuum and resonances", "[gra::MRegge][spin]") {
  gra::LORENTZSCALAR base = MakeToyReggeLTSAsymmetric(0.3, 4.2, -4.6);
  base.process.TU_SIGN    = "positive";

  for (const auto spin : {gra::ReggeProductionModel::MP, gra::ReggeProductionModel::XP}) {
    gra::LORENTZSCALAR explicit_lts = base;
    if (spin == gra::ReggeProductionModel::XP) { PrepareToyXPContinuumOperators(explicit_lts, {{1, 2, 1.0}}); }
    gra::MRegge    explicit_regge(explicit_lts, gra::MModelTune::Load(modelfile),
                                  gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));
    gra::PARAM_RES explicit_res =
        spin == gra::ReggeProductionModel::XP ? MakeToyCovariantXPonance() : MakeToyMPResonance(false);
    TestReggeCon(explicit_regge, explicit_lts, spin);
    auto expected = explicit_lts.hamp;
    TestReggeRes(explicit_regge, explicit_lts, explicit_res, spin);
    REQUIRE(explicit_lts.hamp.size() == expected.size());
    std::transform(explicit_lts.hamp.begin(), explicit_lts.hamp.end(), expected.begin(), expected.begin(),
                   std::plus<>{});

    gra::LORENTZSCALAR sum_lts = base;
    if (spin == gra::ReggeProductionModel::XP) { PrepareToyXPContinuumOperators(sum_lts, {{1, 2, 1.0}}); }
    gra::MRegge    sum_regge(sum_lts, gra::MModelTune::Load(modelfile),
                             gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));
    gra::PARAM_RES sum_res =
        spin == gra::ReggeProductionModel::XP ? MakeToyCovariantXPonance() : MakeToyMPResonance(false);
    sum_lts.process.RESONANCES = {{"toy", sum_res}};
    const double amp2          = sum_regge.Amp2(sum_lts, spin, gra::MReggeMode::ResonanceContinuumTwoBody);

    RequireVectorNear(sum_lts.hamp, expected, 1e-12);
    REQUIRE(amp2 == Approx(sum_lts.hamp.metadata.amplitude_normalization * gra::SquaredNorm(expected)));
  }

  gra::LORENTZSCALAR explicit_lts = base;
  gra::MRegge        explicit_regge(explicit_lts, gra::MModelTune::Load(modelfile),
                                    gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));
  gra::PARAM_RES     explicit_rho = MakeToyMPResonance(false);
  explicit_rho.spin_basis         = "rho";
  explicit_rho.production.front().g = 1.0;
  explicit_rho.MP.filter = gra::HelAmp::IdentityMatrix(3);
  TestReggeCon(explicit_regge, explicit_lts, gra::ReggeProductionModel::MP);
  auto expected = explicit_lts.hamp;
  TestReggeRes(explicit_regge, explicit_lts, explicit_rho, gra::ReggeProductionModel::MP);
  expected.insert(expected.end(), explicit_lts.hamp.begin(), explicit_lts.hamp.end());

  gra::LORENTZSCALAR sum_lts = base;
  gra::MRegge        sum_regge(sum_lts, gra::MModelTune::Load(modelfile),
                               gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));
  gra::PARAM_RES     sum_rho = MakeToyMPResonance(false);
  sum_rho.spin_basis         = "rho";
  sum_rho.production.front().g = 1.0;
  sum_rho.MP.filter = explicit_rho.MP.filter;
  sum_lts.process.RESONANCES = {{"toy", sum_rho}};
  sum_regge.Amp2(sum_lts, gra::ReggeProductionModel::MP, gra::MReggeMode::ResonanceContinuumTwoBody);
  RequireVectorNear(sum_lts.hamp, expected, 1e-12);
}

TEST_CASE("MRegge GP continuum uses GP-compatible proton helicity rows", "[gra::MRegge][spin]") {
  gra::LORENTZSCALAR lts = MakeToyReggeLTSAsymmetric(0.3, 4.2, -4.6);
  lts.process.MMAX       = 1;
  SetToyContinuumExchangePair(lts, 990, 990);
  UseToyGPSubchannelHelicity(lts);
  gra::MRegge regge(lts, gra::MModelTune::Load(modelfile),
                    gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));

  TestReggeCon(regge, lts, gra::ReggeProductionModel::GP);
  REQUIRE(lts.hamp.size() == 4 * 9);
  REQUIRE(std::isfinite(gra::SquaredNorm(lts.hamp)));

  lts.process.FORWARD_NOFLIP = false;
  TestReggeCon(regge, lts, gra::ReggeProductionModel::GP);
  REQUIRE(lts.hamp.size() == 16 * 9);
  REQUIRE(std::isfinite(gra::SquaredNorm(lts.hamp)));
}

TEST_CASE("MRegge continuum pair source equals the explicit t/u source sum", "[gra::MRegge][GoodWalker][continuum]") {
  gra::LORENTZSCALAR lts = MakeToyScalarContinuumLTSAsymmetric(0.3, 4.2, -4.6);
  UpdateToyDerivedKinematics(lts);
  lts.process.TU_SIGN = "positive";
  gra::MRegge regge(lts, gra::MModelTune::Load(modelfile),
                    gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));

  const auto production = gra::mpom::Continuum(lts, 1.0);
  const auto decay      = gra::spin::ContinuumDecayMatrix(lts, lts.process.MP_FRAME);
  REQUIRE(production.size() == 1);
  TestReggeCon(regge, lts, gra::ReggeProductionModel::MP);
  REQUIRE(lts.proton_good_walker.has_value());
  REQUIRE(lts.proton_good_walker->components.size() == 1);
  const auto &actual = lts.proton_good_walker->components.front();
  CHECK(actual.coherence_group == 0);
  CHECK(actual.upper_sector == gra::ProtonGoodWalkerSector::Elastic);
  CHECK(actual.lower_sector == gra::ProtonGoodWalkerSector::Elastic);

  const auto  param     = ReggeParametersForTest(regge, lts);
  const int   upper_pdg = lts.process.CONT_PRODUCTION.front()[0];
  const int   lower_pdg = lts.process.CONT_PRODUCTION.front()[1];
  const auto  upper_id  = param->exchanges.at(gra::regge::TrajectoryIndex(*param, upper_pdg)).soft_exchange;
  const auto  lower_id  = param->exchanges.at(gra::regge::TrajectoryIndex(*param, lower_pdg)).soft_exchange;
  const auto &space     = regge.SoftModelHandle()->GoodWalker();
  std::vector<std::complex<double>> proton(space.ProtonVector().begin(), space.ProtonVector().end());
  const auto pair_profile = gra::KroneckerProduct(regge.SoftModelHandle()->ResidueMatrix(upper_id, lts.t1) * proton,
                                                  regge.SoftModelHandle()->ResidueMatrix(lower_id, lts.t2) * proton);

  const double               mass2_t  = gra::math::pow2(lts.decaytree[0].p.mass);
  const double               mass2_u  = gra::math::pow2(lts.decaytree[1].p.mass);
  const std::vector<int>     pair_pdgs = {lts.decaytree[0].p.pdg, lts.decaytree[1].p.pdg};
  const auto                &entry     = gra::regge::Pair(*param, pair_pdgs, gra::ReggeProductionModel::MP);
  const auto vertex_it = std::find_if(entry.channels.cbegin(), entry.channels.cend(), [&](const auto &vertex) {
    return vertex.first == upper_pdg && vertex.second == lower_pdg;
  });
  // Reproduce the runtime fallback for a synthetic same-trajectory toy exchange
  gra::regge::VertexParam vertex{upper_pdg, lower_pdg};
  if (vertex_it != entry.channels.cend()) {
    vertex = *vertex_it;
  } else {
    vertex.transfer = {entry.forms.at(upper_pdg).transfer, entry.forms.at(lower_pdg).transfer};
    vertex.offshell = {entry.forms.at(upper_pdg).offshell, entry.forms.at(lower_pdg).offshell};
  }
  const std::complex<double> pair_t   = gra::regge::PairExchange(
      *param, entry, vertex, lts.decaytree[0], lts.decaytree[1], lts.t_hat, mass2_t, lts.t1, lts.t2);
  const std::complex<double> pair_u   = gra::regge::PairExchange(
      *param, entry, vertex, lts.decaytree[1], lts.decaytree[0], lts.u_hat, mass2_u, lts.t1, lts.t2);
  const std::complex<double> kernel_t = -regge.ExchangeKernel(lts.ss[1][3], lts.t1, upper_pdg) *
                                        regge.ExchangeKernel(lts.ss[2][4], lts.t2, lower_pdg) * pair_t;
  const std::complex<double> kernel_u = -regge.ExchangeKernel(lts.ss[1][4], lts.t1, upper_pdg) *
                                        regge.ExchangeKernel(lts.ss[2][3], lts.t2, lower_pdg) * pair_u;
  auto central = production[0].first * kernel_t;
  central.AddScaled(production[0].second, kernel_u);
  central = central.MultiplyScaled(decay, gra::regge::Veto(entry.veto, gra::math::msqrt(lts.s_hat)));
  const auto expected = gra::OuterProduct(central.Elements(), pair_profile);
  RequireMatrixNear(actual.source, expected, 1.0e-12);
}

TEST_CASE("MRegge zero fitted flip keeps full and compact Born norms equal", "[gra::MRegge][GoodWalker][helicity]") {
  for (const gra::ReggeProductionModel spin :
       {gra::ReggeProductionModel::MP, gra::ReggeProductionModel::XP, gra::ReggeProductionModel::GP}) {
    CAPTURE(static_cast<int>(spin));
    gra::LORENTZSCALAR compact = MakeToyScalarContinuumLTSAsymmetric(0.3, 4.2, -4.6);
    UpdateToyDerivedKinematics(compact);
    if (spin == gra::ReggeProductionModel::GP) {
      compact.process.MMAX = 1;
      SetToyContinuumExchangePair(compact, 990, 990);
      UseToyGPSubchannelHelicity(compact);
    } else {
      SetToyContinuumExchangePair(compact, 991, 991);
      if (spin == gra::ReggeProductionModel::XP) { PrepareToyXPContinuumOperators(compact, {{0, 0, 1.0}}); }
    }
    compact.process.FORWARD_NOFLIP = true;
    gra::LORENTZSCALAR full        = compact;
    full.process.FORWARD_NOFLIP    = false;

    const auto   model = gra::MModelTune::Load(modelfile);
    gra::MRegge  compact_regge(compact, model, gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));
    gra::MRegge  full_regge(full, model, gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));
    const double compact_amp2 = TestReggeCon(compact_regge, compact, spin);
    const double full_amp2    = TestReggeCon(full_regge, full, spin);

    REQUIRE(compact.hamp.size() == 4);
    REQUIRE(full.hamp.size() == 16);
    REQUIRE(full_amp2 == Approx(compact_amp2).epsilon(2.0e-12));
    for (const std::size_t initial : {0U, 1U, 2U, 3U}) {
      RequireComplexNear(full.hamp[gra::spin::PairHelicityTransitionIndex(initial, initial)], compact.hamp[initial],
                         2.0e-12);
    }
    for (const std::size_t initial : {0U, 1U, 2U, 3U}) {
      for (const std::size_t final : {0U, 1U, 2U, 3U}) {
        if (initial != final) {
          RequireComplexNear(full.hamp[gra::spin::PairHelicityTransitionIndex(initial, final)], 0.0, 2.0e-12);
        }
      }
    }
  }
}

TEST_CASE("MRegge fitted flip sources preserve Good Walker channel structure",
          "[gra::MRegge][GoodWalker][helicity][mapping]") {
  const auto tune  = WriteModifiedPhotoVMTune("regge_forward_flip_profiles", [](auto &card) {
    auto &model                                      = card["PARAM_SOFT"]["MODEL"]["double"];
    card["PARAM_SOFT"]["active_model"]               = "double";
    model["EIKONAL"]["helicity"]                     = true;
    model["EXCHANGE"]["P"]["helicity"]["kappa"]      = nlohmann::json{{0.12, -0.27}, {-0.27, 0.41}};
    model["EXCHANGE"]["P"]["helicity"]["B_kappa"]    = nlohmann::json{{0.31, 0.16}, {0.16, 0.57}};
    model["EXCHANGE"]["R_f2"]["helicity"]["kappa"]   = nlohmann::json{{-0.44, 0.19}, {0.19, 0.23}};
    model["EXCHANGE"]["R_f2"]["helicity"]["B_kappa"] = nlohmann::json{{0.28, 0.47}, {0.47, 0.11}};
  });
  const auto model = gra::MModelTune::Load(tune.second);

  const auto evaluate = [&](const int exchange_alias) {
    gra::LORENTZSCALAR lts = MakeToyScalarContinuumLTSAsymmetric(0.3, 4.2, -4.6);
    UpdateToyDerivedKinematics(lts);
    lts.process.MMAX = 1;
    SetToyContinuumExchangePair(lts, exchange_alias, exchange_alias);
    UseToyGPSubchannelHelicity(lts);
    lts.process.FORWARD_NOFLIP = false;
    gra::MRegge regge(lts, model, gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));
    TestReggeCon(regge, lts, gra::ReggeProductionModel::GP);
    REQUIRE(lts.proton_good_walker.has_value());
    REQUIRE(lts.proton_good_walker->components.size() == 1);
    return std::make_pair(lts, lts.proton_good_walker->components[0].source);
  };

  const auto [pomeron_lts, pomeron_source] = evaluate(990);
  const auto [reggeon_lts, reggeon_source] = evaluate(9910);
  REQUIRE(pomeron_source.size_row() == 16);
  REQUIRE(pomeron_source.size_col() == 4);
  REQUIRE(reggeon_source.size_row() == 16);
  REQUIRE(reggeon_source.size_col() == 4);

  const auto check_profile = [&](const gra::LORENTZSCALAR &lts, const gra::MMatrix<std::complex<double>> &source,
                                 const std::string &exchange_name) {
    std::vector<std::complex<double>> proton(model->Soft()->GoodWalker().ProtonVector().begin(),
                                             model->Soft()->GoodWalker().ProtonVector().end());
    const auto                        exchange      = model->Soft()->ExchangeId(exchange_name);
    const auto                        upper_nonflip = model->Soft()->ResidueMatrix(exchange, lts.t1) * proton;
    auto                              upper_flip = model->Soft()->HelicityFlipResidueMatrix(exchange, lts.t1) * proton;
    gra::Scale(upper_flip, std::sqrt(-lts.t1) / (2.0 * model->Soft()->Eikonal().HelicityMassScale()));
    const auto lower_nonflip = model->Soft()->ResidueMatrix(exchange, lts.t2) * proton;
    auto       lower_flip    = model->Soft()->HelicityFlipResidueMatrix(exchange, lts.t2) * proton;
    gra::Scale(lower_flip, std::sqrt(-lts.t2) / (2.0 * model->Soft()->Eikonal().HelicityMassScale()));
    const std::size_t    nonflip_row        = gra::spin::CanonicalProtonPairSpinLayout::HardRow(0, 0, 0, 0);
    const std::size_t    upper_flip_row     = gra::spin::CanonicalProtonPairSpinLayout::HardRow(0, 0, 1, 0);
    const std::size_t    lower_flip_row     = gra::spin::CanonicalProtonPairSpinLayout::HardRow(0, 0, 0, 1);
    const auto           nonflip_profile    = gra::KroneckerProduct(upper_nonflip, lower_nonflip);
    const auto           upper_flip_profile = gra::KroneckerProduct(upper_flip, lower_nonflip);
    const auto           lower_flip_profile = gra::KroneckerProduct(upper_nonflip, lower_flip);
    std::complex<double> nonflip_scale      = 0.0;
    std::complex<double> upper_flip_scale   = 0.0;
    std::complex<double> lower_flip_scale   = 0.0;
    for (const auto &col : gra::aux::indices(nonflip_profile)) {
      if (std::abs(nonflip_profile[col]) > 1.0e-14) {
        nonflip_scale = source[nonflip_row][col] / nonflip_profile[col];
        break;
      }
    }
    for (const auto &col : gra::aux::indices(upper_flip_profile)) {
      if (std::abs(upper_flip_profile[col]) > 1.0e-14) {
        upper_flip_scale = source[upper_flip_row][col] / upper_flip_profile[col];
        break;
      }
    }
    for (const auto &col : gra::aux::indices(lower_flip_profile)) {
      if (std::abs(lower_flip_profile[col]) > 1.0e-14) {
        lower_flip_scale = source[lower_flip_row][col] / lower_flip_profile[col];
        break;
      }
    }
    REQUIRE(std::abs(nonflip_scale) > 0.0);
    REQUIRE(std::abs(upper_flip_scale) > 0.0);
    REQUIRE(std::abs(lower_flip_scale) > 0.0);
    for (const auto &col : gra::aux::indices(nonflip_profile)) {
      RequireComplexNear(source[nonflip_row][col], nonflip_scale * nonflip_profile[col], 2.0e-12);
      RequireComplexNear(source[upper_flip_row][col], upper_flip_scale * upper_flip_profile[col], 2.0e-12);
      RequireComplexNear(source[lower_flip_row][col], lower_flip_scale * lower_flip_profile[col], 2.0e-12);
    }
    const std::complex<double> expected_upper_section =
        gra::spin::SpinHalfForwardHelicitySectionFactor(-0.5, 0.5, (lts.pfinal[1] - lts.pbeam1).Phi(), false);
    const std::complex<double> expected_lower_section =
        gra::spin::SpinHalfForwardHelicitySectionFactor(-0.5, 0.5, (lts.pbeam2 - lts.pfinal[2]).Phi(), true);
    RequireComplexNear(upper_flip_scale / nonflip_scale, expected_upper_section, 2.0e-12);
    RequireComplexNear(lower_flip_scale / nonflip_scale, expected_lower_section, 2.0e-12);
    std::complex<double> projection = 0.0;
    for (const auto &i : gra::aux::indices(upper_nonflip)) {
      projection += std::conj(upper_nonflip[i]) * upper_flip[i];
    }
    auto residual = upper_flip;
    gra::AddScaled(residual, upper_nonflip, -projection / gra::SquaredNorm(upper_nonflip));
    REQUIRE(gra::SquaredNorm(residual) > 1.0e-12);
  };

  check_profile(pomeron_lts, pomeron_source, "P");
  check_profile(reggeon_lts, reggeon_source, "R_f2");
  REQUIRE(MatrixDiffNorm2(pomeron_source, reggeon_source) > 1.0e-12);
}

TEST_CASE("MRegge rejects invalid generated forward flip transfers", "[gra::MRegge][helicity][failure]") {
  auto base = MakeToyScalarContinuumLTSAsymmetric(0.3, 4.2, -4.6);
  UpdateToyDerivedKinematics(base);
  base.process.MMAX           = 1;
  base.process.FORWARD_NOFLIP = false;
  SetToyContinuumExchangePair(base, 990, 990);
  UseToyGPSubchannelHelicity(base);
  const auto model = gra::MModelTune::Load(modelfile);

  const auto evaluate = [&](const double transfer) {
    auto lts = base;
    lts.t1   = transfer;
    gra::MRegge regge(lts, model, gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));
    return TestReggeCon(regge, lts, gra::ReggeProductionModel::GP);
  };

  REQUIRE_THROWS_AS(evaluate(std::numeric_limits<double>::quiet_NaN()), gra::AmplitudeFailure);
  REQUIRE_THROWS_AS(evaluate(0.5e-12), gra::AmplitudeFailure);
}

TEST_CASE("MRegge GP[RES] coherently sums both mixed Pomeron f2 beam orderings",
          "[gra::MRegge][spin][GP][RES][physics]") {
  gra::LORENTZSCALAR base = MakeToyReggeLTSAsymmetric(0.31, 4.2, -4.6);
  base.process.MMAX       = 1;

  // Build one beam ordered mixed exchange resonance channel
  auto ordered_resonance = [&base](int upper_pdg, int lower_pdg, std::complex<double> coupling) {
    gra::PARAM_RES res          = MakeToyScalarGPResonance(base.process.MMAX);
    res.production[0].tree[0].p = base.PDG.FindByPDG(upper_pdg);
    res.production[0].tree[1].p = base.PDG.FindByPDG(lower_pdg);
    res.production.front().hel.alpha_ls.Scale(coupling);
    res.production.front().hel.T *= coupling;
    return res;
  };

  const std::complex<double> pom_f2_coupling(0.8, 0.2);
  const std::complex<double> f2_pom_coupling(-0.35, 0.15);
  gra::PARAM_RES             pom_f2_res   = ordered_resonance(990, 9910, pom_f2_coupling);
  gra::PARAM_RES             f2_pom_res   = ordered_resonance(9910, 990, f2_pom_coupling);
  gra::PARAM_RES             coherent_res = pom_f2_res;
  coherent_res.production.push_back(f2_pom_res.production.front());

  gra::LORENTZSCALAR pom_f2   = base;
  gra::LORENTZSCALAR f2_pom   = base;
  gra::LORENTZSCALAR coherent = base;
  gra::MRegge        pom_f2_regge(pom_f2, gra::MModelTune::Load(modelfile),
                                  gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));
  gra::MRegge        f2_pom_regge(f2_pom, gra::MModelTune::Load(modelfile),
                                  gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));
  gra::MRegge        coherent_regge(coherent, gra::MModelTune::Load(modelfile),
                                    gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));
  TestReggeRes(pom_f2_regge, pom_f2, pom_f2_res, gra::ReggeProductionModel::GP);
  TestReggeRes(f2_pom_regge, f2_pom, f2_pom_res, gra::ReggeProductionModel::GP);
  TestReggeRes(coherent_regge, coherent, coherent_res, gra::ReggeProductionModel::GP);

  REQUIRE(coherent.hamp.size() == pom_f2.hamp.size());
  REQUIRE(coherent.hamp.size() == f2_pom.hamp.size());
  double ordering_difference = 0.0;
  for (std::size_t i = 0; i < coherent.hamp.size(); ++i) {
    const std::complex<double> expected = pom_f2.hamp[i] + f2_pom.hamp[i];
    const double               scale    = std::max({1.0, std::abs(expected), std::abs(coherent.hamp[i])});
    CHECK(std::abs(coherent.hamp[i] - expected) < 1e-11 * scale);
    ordering_difference += std::norm(pom_f2.hamp[i] - f2_pom.hamp[i]);
  }
  CHECK(FiniteHelicityAmp2(pom_f2.hamp) > 0.0);
  CHECK(FiniteHelicityAmp2(f2_pom.hamp) > 0.0);
  CHECK(ordering_difference > 1e-20);
}

TEST_CASE("MRegge GP[CON] coherently sums both mixed Pomeron f2 beam orderings",
          "[gra::MRegge][spin][GP][CON][physics]") {
  gra::LORENTZSCALAR base  = MakeToyScalarContinuumLTSAsymmetric(0.31, 4.2, -4.6);
  const auto         beams = ProtonInitialState();
  base.beam1               = beams[0];
  base.beam2               = beams[1];
  UpdateToyDerivedKinematics(base);
  base.process.MMAX    = 1;
  base.process.TU_SIGN = "auto";
  SetToyContinuumExchangePair(base, 990, 990);
  UseToyGPSubchannelHelicity(base);

  gra::LORENTZSCALAR pom_f2 = base;
  gra::LORENTZSCALAR f2_pom = base;
  SetToyContinuumExchangePair(pom_f2, 990, 9910);
  SetToyContinuumExchangePair(f2_pom, 9910, 990);
  UseToyGPSubchannelHelicity(pom_f2);
  UseToyGPSubchannelHelicity(f2_pom);

  gra::LORENTZSCALAR coherent = pom_f2;
  coherent.process.CONT_PRODUCTION.push_back(f2_pom.process.CONT_PRODUCTION.front());
  coherent.process.CONT_PRODUCTIONTREE.push_back(f2_pom.process.CONT_PRODUCTIONTREE.front());
  coherent.process.CONTINUUM_GP.push_back(f2_pom.process.CONTINUUM_GP.front());

  gra::MRegge pom_f2_regge(pom_f2, gra::MModelTune::Load(modelfile),
                           gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));
  gra::MRegge f2_pom_regge(f2_pom, gra::MModelTune::Load(modelfile),
                           gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));
  gra::MRegge coherent_regge(coherent, gra::MModelTune::Load(modelfile),
                             gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));
  TestReggeCon(pom_f2_regge, pom_f2, gra::ReggeProductionModel::GP);
  TestReggeCon(f2_pom_regge, f2_pom, gra::ReggeProductionModel::GP);
  TestReggeCon(coherent_regge, coherent, gra::ReggeProductionModel::GP);

  REQUIRE(coherent.hamp.size() == pom_f2.hamp.size());
  REQUIRE(coherent.hamp.size() == f2_pom.hamp.size());
  double ordering_difference = 0.0;
  for (std::size_t i = 0; i < coherent.hamp.size(); ++i) {
    const std::complex<double> expected = pom_f2.hamp[i] + f2_pom.hamp[i];
    const double               scale    = std::max({1.0, std::abs(expected), std::abs(coherent.hamp[i])});
    CHECK(std::abs(coherent.hamp[i] - expected) < 1e-11 * scale);
    ordering_difference += std::norm(pom_f2.hamp[i] - f2_pom.hamp[i]);
  }
  CHECK(FiniteHelicityAmp2(pom_f2.hamp) > 0.0);
  CHECK(FiniteHelicityAmp2(f2_pom.hamp) > 0.0);
  CHECK(ordering_difference > 1e-20);
}

// Check the complete configured exchange sum for every affected continuum
TEST_CASE("Call-local Regge source cache matches direct forward evaluation", "[gra::MRegge][cache][screening]") {
  gra::LORENTZSCALAR    lts = MakeToyContinuumLTSAsymmetric(0.3, 4.2, -4.6);
  gra::MRegge           regge(lts, gra::MModelTune::Load(modelfile),
                              gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));
  const auto            param     = ReggeParametersForTest(regge, lts);
  const auto            upper     = gra::ResolveForwardLegState(lts, gra::ForwardBeamLeg::Upper);
  const auto            lower     = gra::ResolveForwardLegState(lts, gra::ForwardBeamLeg::Lower);
  const int             up_pdg    = lts.process.CONT_PRODUCTION.front()[0];
  const int             dn_pdg    = lts.process.CONT_PRODUCTION.front()[1];
  constexpr double      s_forward = 100.0;
  gra::ReggeSourceCache cache(lts, regge, *param, upper, lower);

  const auto  direct_up = gra::LegProfile(lts, upper, regge, *param, up_pdg);
  const auto  direct_dn = gra::LegSources(lts, lower, regge, *param, dn_pdg, s_forward, false, 0);
  const auto &cached_up = cache.Profile(gra::ForwardBeamLeg::Upper, up_pdg);
  const auto &cached_dn = cache.Sources(gra::ForwardBeamLeg::Lower, dn_pdg, s_forward, false, 0);
  REQUIRE(cached_up.size() == direct_up.size());
  REQUIRE(cached_dn.size() == direct_dn.size());
  for (const auto &i : indices(direct_up)) {
    CHECK(cached_up[i].sector == direct_up[i].sector);
    RequireVectorNear(cached_up[i].nonflip, direct_up[i].nonflip);
    RequireVectorNear(cached_up[i].flip, direct_up[i].flip);
  }
  for (const auto &i : indices(direct_dn)) {
    CHECK(cached_dn[i].sector == direct_dn[i].sector);
    RequireVectorNear(cached_dn[i].nonflip, direct_dn[i].nonflip);
    RequireVectorNear(cached_dn[i].flip, direct_dn[i].flip);
  }

  const auto direct_kernel = gra::LegKernel(lts, upper, regge, *param, up_pdg, s_forward, false, 0);
  RequireComplexNear(cache.Kernel(gra::ForwardBeamLeg::Upper, up_pdg, s_forward, false, 0), direct_kernel);
  CHECK(&cached_up == &cache.Profile(gra::ForwardBeamLeg::Upper, up_pdg));
  CHECK(&cached_dn == &cache.Sources(gra::ForwardBeamLeg::Lower, dn_pdg, s_forward, false, 0));
}

TEST_CASE("MP XP and GP continuum exchange channels add coherently",
          "[gra::MRegge][continuum][exchange][coherence][physics]") {
  ModelParamRestoreGuard restore;
  const auto             tune = WriteModifiedPhotoVMTune("continuum_exchange_coherence", SetSecondaryContinuumChannels);
  const auto             model_tune = gra::MModelTune::Load(tune.second);
  gra::MODELPARAM                   = tune.first;
  struct ContinuumCase {
    const char               *family;
    gra::ReggeProductionModel model;
    const char               *final_state;
    std::array<int, 2>        pdg;
  };
  const std::array<ContinuumCase, 9> cases = {{
      {"MP", gra::ReggeProductionModel::MP, "pi+ pi-", {211, -211}},
      {"MP", gra::ReggeProductionModel::MP, "K+ K-", {321, -321}},
      {"MP", gra::ReggeProductionModel::MP, "p+ p-", {2212, -2212}},
      {"XP", gra::ReggeProductionModel::XP, "pi+ pi-", {211, -211}},
      {"XP", gra::ReggeProductionModel::XP, "K+ K-", {321, -321}},
      {"XP", gra::ReggeProductionModel::XP, "p+ p-", {2212, -2212}},
      {"GP", gra::ReggeProductionModel::GP, "pi+ pi-", {211, -211}},
      {"GP", gra::ReggeProductionModel::GP, "K+ K-", {321, -321}},
      {"GP", gra::ReggeProductionModel::GP, "p+ p-", {2212, -2212}},
  }};

  for (const auto &test : cases) {
    CAPTURE(test.family, test.final_state);
    ToyHelicityProcess process;
    ConfigureToyProductionProcess(process, test.family, "CON", test.final_state);
    process.SetModelTune(model_tune);
    REQUIRE_NOTHROW(process.InitializeProcessAmplitude());
    const auto configured = process.state.lts.process;
    REQUIRE(configured.CONT_PRODUCTION.size() > 1);
    REQUIRE(configured.CONT_PRODUCTION.size() == configured.CONT_PRODUCTIONTREE.size());
    const std::size_t cache_channels =
        test.model == gra::ReggeProductionModel::GP ? configured.CONTINUUM_GP.size() : configured.CONTINUUM_POLE.size();
    REQUIRE(configured.CONT_PRODUCTION.size() == cache_channels);
    REQUIRE(configured.CONT_PRODUCTION.size() == configured.CONT_TU_SIGN.size());

    const auto evaluate = [&](const std::vector<std::size_t> &selected) {
      gra::LORENTZSCALAR lts = DirectCentralPairLTSForTest(test.pdg[0], test.pdg[1]);
      RefreshToyDerivedKinematicsPreserveDecay(lts);
      lts.process = configured;
      lts.hamp.Configure(process.state.lts.hamp.metadata);
      lts.process.CONT_PRODUCTION.clear();
      lts.process.CONT_PRODUCTIONTREE.clear();
      lts.process.CONTINUUM_POLE.clear();
      lts.process.CONTINUUM_GP.clear();
      lts.process.CONT_TU_SIGN.clear();
      for (const std::size_t channel : selected) {
        lts.process.CONT_PRODUCTION.push_back(configured.CONT_PRODUCTION.at(channel));
        lts.process.CONT_PRODUCTIONTREE.push_back(configured.CONT_PRODUCTIONTREE.at(channel));
        if (test.model == gra::ReggeProductionModel::GP) {
          lts.process.CONTINUUM_GP.push_back(configured.CONTINUUM_GP.at(channel));
        } else {
          lts.process.CONTINUUM_POLE.push_back(configured.CONTINUUM_POLE.at(channel));
        }
        lts.process.CONT_TU_SIGN.push_back(configured.CONT_TU_SIGN.at(channel));
      }
      gra::MRegge regge(lts, process.state.model_tune,
                        gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::ContinuumTwoBody,
                                                          std::string(test.family) + "_exchange_coherence"));
      TestReggeCon(regge, lts, test.model);
      REQUIRE(lts.proton_good_walker.has_value());
      REQUIRE(lts.proton_good_walker->components.size() == 1);
      REQUIRE(gra::SquaredNorm(lts.hamp) > 0.0);
      return std::pair{std::vector<std::complex<double>>(lts.hamp.begin(), lts.hamp.end()),
                       lts.proton_good_walker->components.front().source};
    };

    std::vector<std::size_t> all(configured.CONT_PRODUCTION.size());
    std::iota(all.begin(), all.end(), 0U);
    const auto                         full = evaluate(all);
    std::vector<std::complex<double>>  amplitude_sum(full.first.size(), 0.0);
    gra::MMatrix<std::complex<double>> source_sum(full.second.size_row(), full.second.size_col(), 0.0);
    for (const std::size_t channel : all) {
      const auto isolated = evaluate({channel});
      REQUIRE(isolated.first.size() == amplitude_sum.size());
      REQUIRE(isolated.second.size_row() == source_sum.size_row());
      REQUIRE(isolated.second.size_col() == source_sum.size_col());
      gra::AddScaled(amplitude_sum, isolated.first, std::complex<double>(1.0, 0.0));
      source_sum += isolated.second;
    }
    RequireVectorNear(full.first, amplitude_sum, 3.0e-11);
    RequireMatrixNear(full.second, source_sum, 3.0e-11);
  }
}

TEST_CASE("MRegge GP proton residues generate expected dPhi harmonics", "[gra::MRegge][spin]") {
  auto harmonic = [](int n, double phi1, double phi2) {
    return gra::gpom::Residue(n, 0.0, 1.0, phi1, 1.0, false, false) *
               gra::gpom::Residue(n, 0.0, 1.0, phi2, 1.0, true, false) +
           gra::gpom::Residue(-n, 0.0, 1.0, phi1, 1.0, false, false) *
               gra::gpom::Residue(-n, 0.0, 1.0, phi2, 1.0, true, false);
  };

  for (const int n : {1, 2}) {
    const double phi1 = 0.4;
    const double phi2 = -0.7;
    RequireComplexNear(harmonic(n, phi1, phi2), std::complex<double>(2.0 * std::cos(n * (phi1 - phi2)), 0.0));
    RequireComplexNear(harmonic(n, 0.0, 0.0), std::complex<double>(2.0, 0.0));
    RequireComplexNear(harmonic(n, 0.0, gra::math::PI / n), std::complex<double>(-2.0, 0.0));
  }
}

TEST_CASE("GP finite m continuum has the exact reduced m zero limit",
          "[gra::MRegge][GP][CON][spin][normalization][regression]") {
  gra::LORENTZSCALAR base = MakeToyScalarContinuumLTSAsymmetric(0.31, 4.2, -4.6);
  UpdateToyDerivedKinematics(base);
  base.process.MMAX           = 2;
  base.process.FORWARD_NOFLIP = true;
  base.process.TU_SIGN        = "positive";
  SetToyContinuumExchangePair(base, 990, 990);
  const std::complex<double> residue(0.73, -0.28);

  gra::LORENTZSCALAR reduced = base;
  UseScalarZeroMContinuum(reduced, residue);
  gra::LORENTZSCALAR finite = base;
  UseScalarFiniteMContinuum(finite, {{0, residue}});

  gra::MRegge regge(base, gra::MModelTune::Load(modelfile),
                    gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));
  const auto  param        = ReggeParametersForTest(regge, base);
  const auto  reduced_prod = gra::gpom::Continuum(reduced, *param);
  const auto  finite_prod  = gra::gpom::Continuum(finite, *param);
  REQUIRE(reduced_prod.size() == 1);
  REQUIRE(finite_prod.size() == 1);
  RequireMatrixNear(finite_prod[0].first, reduced_prod[0].first, 2.0e-12);
  RequireMatrixNear(finite_prod[0].second, reduced_prod[0].second, 2.0e-12);
}

TEST_CASE("Finite m crossed LS and helicity bases share the spin two pole",
          "[gra::MRegge][MP][XP][GP][CON][spin][normalization]"
          "[unit_coupling][regression]") {
  constexpr int                                           mmax        = 2;
  const std::vector<std::pair<int, std::complex<double>>> coefficient = {
      {-2, {0.37, -0.14}}, {-1, {-0.21, 0.09}}, {0, {0.73, -0.28}}, {1, {-0.21, 0.09}}, {2, {0.37, -0.14}}};

  gra::HELMatrix crossed_ls;
  gra::gpom::InitCrossed(crossed_ls, 0.0, 0.0, mmax, "finite-m crossed LS pole normalization");
  crossed_ls.C_symmetry     = true;
  crossed_ls.P_symmetry     = true;
  crossed_ls.coupling_basis = gra::CouplingBasis::LS;
  crossed_ls.m_ls.assign(2 * mmax + 1, {});
  for (const auto &[m, value] : coefficient) {
    const std::size_t column = gra::gpom::AnalyticMIndex(m, mmax, "finite-m crossed LS pole normalization");
    crossed_ls.m_ls[column].Set(2, 0, value);
  }
  gra::gpom::InitCrossedLS(crossed_ls, 2);

  const double norm                 = std::sqrt(2.0 / 3.0);
  auto         helicity_coefficient = coefficient;
  for (auto &entry : helicity_coefficient) { entry.second *= norm; }
  const gra::HELMatrix crossed_helicity = ScalarFiniteMCrossed(mmax, helicity_coefficient);
  const gra::M4Vec     axis(0.7, 0.0, 0.9, 1.5);
  const auto           ls       = gra::gpom::Crossed(crossed_ls, 2.0, axis, false);
  const auto           helicity = gra::gpom::Crossed(crossed_helicity, 2.0, axis, false);
  for (const auto &[m, value] : coefficient) {
    const std::size_t column = gra::gpom::AnalyticMIndex(m, mmax, "finite-m crossed basis pole comparison");
    RequireComplexNear(ls.helicity.T[0][column], value * norm, 2.0e-12);
    RequireComplexNear(helicity.helicity.T[0][column], value * norm, 2.0e-12);
    RequireComplexNear(ls.frame[0][column], helicity.frame[0][column], 2.0e-12);
  }

  const auto                 pomeron  = LoadedPDGTable().FindByPDG(995);
  const auto                 pion     = LoadedPDGTable().FindByPDG(211);
  const auto                 antipion = LoadedPDGTable().FindByPDG(-211);
  const std::complex<double> g        = coefficient[2].second;
  const auto fixed_vertex = gra::spin::PreparePoleLS(pomeron, pion, antipion, {{2, 0, g}}, 1.0, true, true, true,
                                                     gra::spin::VertexContext::SubTUChannelExchange, 0.0, false);
  const auto fixed        = gra::spin::PoleLSHelicity(fixed_vertex, 1.0);
  REQUIRE(fixed.T.size_row() == 1);
  REQUIRE(fixed.T.size_col() == 1);
  const std::size_t zero = gra::gpom::AnalyticMIndex(0, mmax, "finite-m crossed fixed-pole comparison");
  RequireComplexNear(fixed.T[0][0], g * norm, 2.0e-12);
  RequireComplexNear(ls.helicity.T[0][zero], fixed.T[0][0], 2.0e-12);
  RequireComplexNear(helicity.helicity.T[0][zero], fixed.T[0][0], 2.0e-12);

  // The isotropic approximation uses physical pole coefficients in both input bases
  auto lts = MakeToyScalarContinuumLTSAsymmetric(0.3, 4.2, -4.6);
  UpdateToyDerivedKinematics(lts);
  SetToyContinuumExchangePair(lts, 990, 990);
  lts.process.SPINGEN = false;
  lts.process.MMAX = mmax;
  lts.process.CONTINUUM_GP = {std::vector<gra::HELMatrix>(4, crossed_ls)};
  gra::MRegge regge(lts, gra::MModelTune::Load(modelfile),
                    gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "crossed pole scale"));
  const auto param = ReggeParametersForTest(regge, lts);
  const auto blind_ls = gra::gpom::Continuum(lts, *param);
  lts.process.CONTINUUM_GP = {std::vector<gra::HELMatrix>(4, crossed_helicity)};
  const auto blind_helicity = gra::gpom::Continuum(lts, *param);
  RequireMatrixNear(blind_ls[0].first, blind_helicity[0].first, 2.0e-12);
  RequireMatrixNear(blind_ls[0].second, blind_helicity[0].second, 2.0e-12);

  // The single m zero vertex also matches MP and XP at the common pole
  auto zero_ls = crossed_ls;
  for (const auto &m : indices(zero_ls.m_ls)) {
    if (m != zero) { zero_ls.m_ls[m].Clear(); }
  }
  gra::gpom::InitCrossedLS(zero_ls, 2);
  lts.process.CONTINUUM_GP = {std::vector<gra::HELMatrix>(4, zero_ls)};
  const auto blind_zero = gra::gpom::Continuum(lts, *param);
  auto fixed_lts = lts;
  SetToyContinuumExchangePair(fixed_lts, 995, 995);
  fixed_lts.process.CONTINUUM_POLE = {std::vector<gra::spin::PoleResidue>(4, gra::spin::PoleResidue(fixed_vertex))};
  for (const auto model : {gra::ReggeProductionModel::MP, gra::ReggeProductionModel::XP}) {
    const auto blind =
        (model == gra::ReggeProductionModel::MP ? gra::mpom::Continuum : gra::xpom::Continuum)(fixed_lts, param->s0);
    RequireMatrixNear(blind_zero[0].first, blind[0].first, 2.0e-12);
    RequireMatrixNear(blind_zero[0].second, blind[0].second, 2.0e-12);
  }
}

// Compare prepared poles with full contractions at actual screening momentum shifts
TEST_CASE("MP XP GP continuum contractions preserve full amplitudes at screening nodes",
          "[gra::MRegge][MP][XP][GP][continuum][cache][screening]") {
  struct FinalPair {
    const char *name;
    int         first;
    int         second;
  };
  for (const auto model :
       {gra::ReggeProductionModel::MP, gra::ReggeProductionModel::XP, gra::ReggeProductionModel::GP}) {
    for (const FinalPair pair : {FinalPair{"pi+ pi-", 211, -211}, FinalPair{"p+ p-", 2212, -2212},
                                 FinalPair{"rho(770)0 rho(770)0", 113, 113}, FinalPair{"phi(1020)0 phi(1020)0", 333, 333}}) {
      CAPTURE(model, pair.name);
      ToyHelicityProcess process;
      ConfigureToyProductionProcess(process, gra::ReggeProductionModelName(model), "CON", pair.name);
      process.InitializeProcessAmplitude();
      auto base        = AsymmetricContinuumPairForTest(pair.first, pair.second);
      base.process     = process.state.lts.process;
      const auto param = gra::regge::GetParam(*process.state.lts.model_cache, {pair.first, pair.second}, base.PDG);
      for (const bool noflip : {false, true}) {
        for (const std::array<double, 2> shift : {std::array<double, 2>{0.0, 0.0}, {0.07, -0.03}, {-0.04, 0.11}}) {
          CAPTURE(noflip, shift[0], shift[1]);
          auto state                       = process.state;
          state.lts                        = base;
          state.lts.process.FORWARD_NOFLIP = noflip;
          PrepareScreeningPoint(state);
          state.lts.pfinal_orig = state.lts.pfinal;
          REQUIRE(gra::kinematics::RebuildScreeningKinematics(
              state.lts, {base.pfinal[1].Px() - shift[0], base.pfinal[1].Py() - shift[1]},
              {base.pfinal[2].Px() + shift[0], base.pfinal[2].Py() + shift[1]}, true));
          REQUIRE(gra::kinematics::SetLorentzScalars(state, 4, true));
          gra::gpom::AmpCache cache;
          const auto          actual =
              model == gra::ReggeProductionModel::GP
                           ? gra::gpom::Continuum(state.lts, *param, &cache)
                           : (model == gra::ReggeProductionModel::MP ? gra::mpom::Continuum : gra::xpom::Continuum)(state.lts,
                                                                                                           param->s0);
          for (const auto &channel : indices(actual)) {
            const auto full = model == gra::ReggeProductionModel::GP
                                  ? FullGPContinuum(state.lts, cache, channel)
                                  : FullPoleContinuum(state.lts, model, param->s0, channel);
            RequireMatrixNear(actual[channel].first, full.first, 3.0e-11);
            RequireMatrixNear(actual[channel].second, full.second, 3.0e-11);
            if (model == gra::ReggeProductionModel::XP) {
              // Identical vertices must give the same complex t/u amplitudes in MP and XP
              const auto mp = FullPoleContinuum(state.lts, gra::ReggeProductionModel::MP, param->s0, channel);
              RequireMatrixNear(actual[channel].first, mp.first, 3.0e-11);
              RequireMatrixNear(actual[channel].second, mp.second, 3.0e-11);
            }
          }
        }
      }
    }
  }
}

TEST_CASE("GP finite m continuum is padded covariantly and shares its cache",
          "[gra::MRegge][GP][CON][RES][spin][covariance][cache]"
          "[screening][regression]") {
  const std::vector<std::pair<int, std::complex<double>>> residue = {
      {-2, {0.41, -0.17}}, {-1, {-0.23, 0.11}}, {0, {0.79, 0.26}}, {1, {-0.23, 0.11}}, {2, {0.41, -0.17}}};
  const auto make_lts = [&](const int mmax) {
    gra::LORENTZSCALAR lts = MakeToyScalarContinuumLTSAsymmetric(0.31, 4.2, -4.6);
    UpdateToyDerivedKinematics(lts);
    lts.process.MMAX           = mmax;
    lts.process.FORWARD_NOFLIP = true;
    lts.process.TU_SIGN        = "positive";
    SetToyContinuumExchangePair(lts, 990, 990);
    UseScalarFiniteMContinuum(lts, residue);
    return lts;
  };

  gra::LORENTZSCALAR base = make_lts(2);
  gra::MRegge        regge(base, gra::MModelTune::Load(modelfile),
                           gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));
  const auto         param     = ReggeParametersForTest(regge, base);
  const auto         reference = gra::gpom::Continuum(base, *param);
  REQUIRE(reference.size() == 1);
  REQUIRE(reference[0].first.FrobNorm2() > 0.0);

  const std::complex<double> upper_section = std::polar(1.0, 0.37);
  const std::complex<double> lower_section = std::polar(1.0, -0.61);
  const std::complex<double> pair_section  = upper_section * lower_section;
  gra::LORENTZSCALAR         phased        = base;
  REQUIRE(phased.process.CONTINUUM_GP.size() == 1);
  REQUIRE(phased.process.CONTINUUM_GP[0].size() == 4);
  for (const std::size_t vertex : {0U, 2U}) { phased.process.CONTINUUM_GP[0][vertex].T *= upper_section; }
  for (const std::size_t vertex : {1U, 3U}) { phased.process.CONTINUUM_GP[0][vertex].T *= lower_section; }
  const auto phased_prod = gra::gpom::Continuum(phased, *param);
  REQUIRE(phased_prod.size() == 1);
  RequireMatrixNear(phased_prod[0].first, pair_section * reference[0].first, 2.0e-12);
  RequireMatrixNear(phased_prod[0].second, pair_section * reference[0].second, 2.0e-12);

  gra::LORENTZSCALAR padded      = make_lts(4);
  const auto         padded_prod = gra::gpom::Continuum(padded, *param);
  REQUIRE(padded_prod.size() == 1);
  RequireMatrixNear(padded_prod[0].first, reference[0].first, 2.0e-12);
  RequireMatrixNear(padded_prod[0].second, reference[0].second, 2.0e-12);

  for (const double angle : {-1.21, 0.37, 2.04}) {
    const gra::LORENTZSCALAR rotated      = RotateToyEventAroundZ(base, angle);
    const auto               rotated_prod = gra::gpom::Continuum(rotated, *param);
    REQUIRE(rotated_prod.size() == 1);
    CAPTURE(angle);
    RequireMatrixNear(rotated_prod[0].first, reference[0].first, 2.0e-11);
    RequireMatrixNear(rotated_prod[0].second, reference[0].second, 2.0e-11);
  }

  gra::gpom::AmpCache cache;
  const auto          cached_continuum = gra::gpom::Continuum(base, *param, &cache);
  REQUIRE(cache.sources.size() == 1);
  gra::PARAM_RES resonance        = MakeToyGPResonance(base.process.MMAX);
  const auto     cached_resonance = gra::gpom::Resonance(base, *param, resonance, &cache);
  REQUIRE(cache.sources.size() == 1);
  const auto direct_resonance = gra::gpom::Resonance(base, *param, resonance, nullptr);
  REQUIRE(cached_continuum.size() == reference.size());
  REQUIRE(cached_resonance.size() == 1);
  REQUIRE(direct_resonance.size() == 1);
  RequireMatrixNear(cached_continuum[0].first, reference[0].first, 2.0e-13);
  RequireMatrixNear(cached_continuum[0].second, reference[0].second, 2.0e-13);
  RequireMatrixNear(cached_resonance[0], direct_resonance[0], 2.0e-13);
}

TEST_CASE("GP mixed photon finite m continuum is padded covariantly",
          "[gra::MRegge][GP][CON][spin][photon][covariance][regression]") {
  const auto tune = WriteModifiedPhotoVMTune("gp_mixed_photon", {}, {}, SetToyPhotonContinuum);
  constexpr auto                                          helicity_x2 = gra::spin::BinaryHelicityLabelsX2();
  const std::vector<std::pair<int, std::complex<double>>> residue     = {
          {-2, {0.41, -0.17}}, {-1, {-0.23, 0.11}}, {0, {0.79, 0.26}}, {1, {-0.23, 0.11}}, {2, {0.41, -0.17}}};

  const auto make_lts = [&](const int mmax, const bool lower_photon, const std::string &photon_vertex) {
    gra::LORENTZSCALAR lts = MakeToyScalarContinuumLTSAsymmetric(0.31, 4.2, -4.6);
    UpdateToyDerivedKinematics(lts);
    lts.process.MMAX           = mmax;
    lts.process.FORWARD_NOFLIP = false;
    lts.process.FORWARD_VERTEX = gra::ForwardVertexMode::HelicityResidue;
    lts.process.PHOTON_VERTEX  = photon_vertex;
    lts.process.TU_SIGN        = "positive";
    const int upper            = lower_photon ? 990 : gra::PDG::PDG_gamma;
    const int lower            = lower_photon ? gra::PDG::PDG_gamma : 990;
    SetToyContinuumExchangePair(lts, upper, lower);
    UseScalarMixedFiniteMContinuum(lts, residue);
    return lts;
  };

  const auto evaluate = [&](gra::LORENTZSCALAR lts) {
    gra::MRegge regge(lts, gra::MModelTune::Load(tune.second),
                      gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "mixed photon finite-m test"));
    TestReggeCon(regge, lts, gra::ReggeProductionModel::GP);
    return lts.hamp;
  };

  for (const std::string photon_vertex : {"EPA", "QED"}) {
    for (const bool lower_photon : {false, true}) {
      CAPTURE(photon_vertex, lower_photon);
      const gra::LORENTZSCALAR base      = make_lts(2, lower_photon, photon_vertex);
      const auto               reference = evaluate(base);
      REQUIRE(reference.size() == 16);
      REQUIRE(std::isfinite(gra::SquaredNorm(reference)));
      REQUIRE(gra::SquaredNorm(reference) > 0.0);

      gra::LORENTZSCALAR zero_only = base;
      UseScalarMixedFiniteMContinuum(zero_only, {{0, residue[2].second}});
      const auto zero                = evaluate(zero_only);
      double     finite_m_difference = 0.0;
      for (const auto &i : indices(reference)) { finite_m_difference += std::norm(reference[i] - zero[i]); }
      REQUIRE(finite_m_difference > 1.0e-20);

      for (const int mmax : {3, 4}) {
        CAPTURE(mmax);
        const auto padded = evaluate(make_lts(mmax, lower_photon, photon_vertex));
        RequireVectorNear(padded, reference, 3.0e-11);
      }

      for (const double angle : {-1.21, 0.37, 2.04}) {
        const auto rotated = evaluate(RotateToyEventAroundZ(base, angle));
        REQUIRE(rotated.size() == reference.size());
        for (std::size_t in1 = 0; in1 < 2; ++in1) {
          for (std::size_t in2 = 0; in2 < 2; ++in2) {
            for (std::size_t out1 = 0; out1 < 2; ++out1) {
              for (std::size_t out2 = 0; out2 < 2; ++out2) {
                const std::size_t row = gra::spin::CanonicalProtonPairSpinLayout::HardRow(in1, in2, out1, out2);
                const int harmonic    = gra::spin::ColliderSpinHalfHelicityHarmonic(helicity_x2[in1], helicity_x2[in2],
                                                                                    helicity_x2[out1], helicity_x2[out2]);
                const std::complex<double> expected =
                    reference[row] * std::exp(gra::math::zi * static_cast<double>(harmonic) * angle);
                CAPTURE(angle, in1, in2, out1, out2, row, harmonic, rotated[row], expected);
                CHECK(std::abs(rotated[row] - expected) <= 5.0e-10 * std::max(1.0, std::abs(expected)));
              }
            }
          }
        }
      }
    }
  }
}

TEST_CASE("MRegge GP continuum is invariant under common azimuth rotations", "[gra::MRegge][spin]") {
  gra::LORENTZSCALAR base     = MakeToyReggeLTSAsymmetric(0.31, 4.2, -4.6);
  base.process.MMAX           = 2;
  base.process.FORWARD_NOFLIP = true;
  base.process.TU_SIGN        = "positive";
  SetToyContinuumExchangePair(base, 990, 990);
  UseToyGPSubchannelHelicity(base);

  auto evaluate = [](gra::LORENTZSCALAR lts) {
    gra::MRegge regge(lts, gra::MModelTune::Load(modelfile),
                      gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));
    TestReggeCon(regge, lts, gra::ReggeProductionModel::GP);
    return FiniteHelicityAmp2(lts.hamp);
  };

  const double reference = evaluate(base);
  REQUIRE(reference > 0.0);

  for (const double angle : {0.37, -1.21, 2.04}) {
    const gra::LORENTZSCALAR rotated = RotateToyEventAroundZ(base, angle);
    const double             actual  = evaluate(rotated);
    CAPTURE(angle, reference, actual);
    REQUIRE(actual / reference == Approx(1.0).epsilon(1e-10));
  }
}

TEST_CASE("MRegge GP continuum is independent of resonance MMAX", "[gra::MRegge][spin][GP][continuum]") {
  const gra::LORENTZSCALAR base = MakeToyReggeLTSAsymmetric(0.31, 4.2, -4.6);
  std::vector<double>      amp2;
  for (const int MMAX : {2, 3, 4, 5}) {
    gra::LORENTZSCALAR lts     = base;
    lts.process.MMAX           = MMAX;
    lts.process.FORWARD_NOFLIP = true;
    lts.process.TU_SIGN        = "positive";
    SetToyContinuumExchangePair(lts, 990, 990);
    UseToyGPSubchannelHelicity(lts);
    gra::MRegge regge(lts, gra::MModelTune::Load(modelfile),
                      gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));
    TestReggeCon(regge, lts, gra::ReggeProductionModel::GP);
    amp2.push_back(FiniteHelicityAmp2(lts.hamp));
    CAPTURE(MMAX, amp2.back());
    REQUIRE(std::isfinite(amp2.back()));
    REQUIRE(amp2.back() > 0.0);
  }

  for (std::size_t index = 1; index < amp2.size(); ++index) {
    CAPTURE(amp2, index);
    CHECK(amp2[index] / amp2.front() == Approx(1.0).epsilon(2.0e-13));
  }
}

TEST_CASE("MRegge crossed continuum is invariant under common azimuth rotations", "[gra::MRegge][spin]") {
  for (const gra::ForwardVertexMode forward_vertex :
       {gra::ForwardVertexMode::HelicityResidue, gra::ForwardVertexMode::UnitResidue}) {
    CAPTURE(static_cast<int>(forward_vertex));
    gra::LORENTZSCALAR base     = MakeToyReggeLTSAsymmetric(0.31, 4.2, -4.6);
    base.process.FORWARD_VERTEX = forward_vertex;
    base.process.FORWARD_NOFLIP = true;
    base.process.TU_SIGN        = "positive";
    SetToyContinuumExchangePair(base, 993, 993);

    auto evaluate = [](gra::LORENTZSCALAR lts) {
      gra::MRegge regge(lts, gra::MModelTune::Load(modelfile),
                        gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));
      TestReggeCon(regge, lts, gra::ReggeProductionModel::MP);
      return FiniteHelicityAmp2(lts.hamp);
    };

    const double reference = evaluate(base);
    REQUIRE(reference > 0.0);

    for (const double angle : {0.37, -1.21, 2.04}) {
      const gra::LORENTZSCALAR rotated = RotateToyEventAroundZ(base, angle);
      const double             actual  = evaluate(rotated);
      CAPTURE(angle, reference, actual);
      REQUIRE(actual / reference == Approx(1.0).epsilon(1e-10));
    }
  }
}

TEST_CASE("MRegge crossed continuum exposes every collider harmonic", "[gra::MRegge][spin][helicity][covariance]") {
  constexpr auto helicity_x2 = gra::spin::BinaryHelicityLabelsX2();

  for (const gra::ForwardVertexMode forward_vertex :
       {gra::ForwardVertexMode::HelicityResidue, gra::ForwardVertexMode::UnitResidue}) {
    for (const auto spin : {gra::ReggeProductionModel::MP, gra::ReggeProductionModel::XP}) {
      const bool xp_mode = spin == gra::ReggeProductionModel::XP;
      CAPTURE(static_cast<int>(forward_vertex), xp_mode);
      gra::LORENTZSCALAR base = MakeToyScalarContinuumLTSAsymmetric(0.31, 4.2, -4.6);
      UpdateToyDerivedKinematics(base);
      base.process.FORWARD_VERTEX = forward_vertex;
      base.process.FORWARD_NOFLIP = false;
      base.process.TU_SIGN        = "positive";
      if (spin == gra::ReggeProductionModel::XP) { PrepareToyXPContinuumOperators(base, {{0, 0, 1.0}}); }

      // Evaluate the public continuum amplitude in the full proton spin basis
      const auto evaluate = [spin](gra::LORENTZSCALAR lts) {
        gra::MRegge regge(lts, gra::MModelTune::Load(modelfile),
                          gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));
        TestReggeCon(regge, lts, spin);
        return lts.hamp;
      };

      const auto reference = evaluate(base);
      REQUIRE(reference.size() == 16);
      REQUIRE(gra::SquaredNorm(reference) > 0.0);

      for (const double angle : {-1.21, 0.37, 2.04}) {
        const auto rotated = evaluate(RotateToyEventAroundZ(base, angle));
        REQUIRE(rotated.size() == reference.size());
        for (std::size_t in1 = 0; in1 < 2; ++in1) {
          for (std::size_t in2 = 0; in2 < 2; ++in2) {
            for (std::size_t out1 = 0; out1 < 2; ++out1) {
              for (std::size_t out2 = 0; out2 < 2; ++out2) {
                const std::size_t row = gra::spin::CanonicalProtonPairSpinLayout::HardRow(in1, in2, out1, out2);
                const int harmonic    = gra::spin::ColliderSpinHalfHelicityHarmonic(helicity_x2[in1], helicity_x2[in2],
                                                                                    helicity_x2[out1], helicity_x2[out2]);
                const std::complex<double> expected =
                    reference[row] * std::exp(gra::math::zi * static_cast<double>(harmonic) * angle);
                CAPTURE(angle, in1, in2, out1, out2, row, harmonic, rotated[row], expected);
                CHECK(std::abs(rotated[row] - expected) <= 2.0e-10 * std::max(1.0, std::abs(expected)));
              }
            }
          }
        }
      }
    }
  }
}

TEST_CASE("MRegge fixed-spin resonances are invariant under common azimuth rotations", "[gra::MRegge][spin]") {
  for (const gra::ForwardVertexMode forward_vertex :
       {gra::ForwardVertexMode::HelicityResidue, gra::ForwardVertexMode::UnitResidue}) {
    CAPTURE(static_cast<int>(forward_vertex));
    gra::LORENTZSCALAR base     = MakeToyReggeLTSAsymmetric(0.31, 4.2, -4.6);
    base.process.FORWARD_VERTEX = forward_vertex;
    base.process.FORWARD_NOFLIP = true;

    auto evaluate_res = [](gra::LORENTZSCALAR lts, bool xp_mode) {
      gra::MRegge    regge(lts, gra::MModelTune::Load(modelfile),
                           gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));
      gra::PARAM_RES res = xp_mode ? MakeToyCovariantXPonance() : MakeToyMPResonance(false);
      if (xp_mode) {
        TestReggeRes(regge, lts, res, gra::ReggeProductionModel::XP);
      } else {
        TestReggeRes(regge, lts, res, gra::ReggeProductionModel::MP);
      }
      return FiniteHelicityAmp2(lts.hamp);
    };

    for (const bool xp_mode : {false, true}) {
      CAPTURE(xp_mode);
      const double reference = evaluate_res(base, xp_mode);
      REQUIRE(reference > 0.0);

      for (const double angle : {0.37, -1.21, 2.04}) {
        const gra::LORENTZSCALAR rotated = RotateToyEventAroundZ(base, angle);
        const double             actual  = evaluate_res(rotated, xp_mode);
        CAPTURE(angle, reference, actual);
        REQUIRE(actual / reference == Approx(1.0).epsilon(1e-10));
      }
    }
  }
}

// Check the complete Jacob-Wick phase before summing resonance spin components
TEST_CASE("MP and XP resonance components obey collider azimuth covariance",
          "[gra::MRegge][MP][XP][spin][helicity][resonance][covariance]"
          "[parity][regression]") {
  constexpr auto     helicity_x2 = gra::spin::BinaryHelicityLabelsX2();
  gra::LORENTZSCALAR base        = MakeToyProductionLTSAsymmetric(0.31, 4.2, -4.6);
  base.process.MP_FRAME          = "CM";
  base.process.SPINGEN           = true;
  base.process.FORWARD_NOFLIP    = false;
  base.process.FORWARD_VERTEX    = gra::ForwardVertexMode::HelicityResidue;
  const auto exchange            = base.PDG.FindByPDG(995);
  REQUIRE(exchange.spinX2 == 4);

  for (int spin = 0; spin <= 6; ++spin) {
    for (const int parity : {-1, 1}) {
      CAPTURE(spin, parity);
      gra::PARAM_RES seed = MakeToyMPResonance(false);
      seed.p.pdg          = 920000 + 10 * spin + parity;
      seed.p.spinX2       = 2 * spin;
      seed.p.P            = parity;
      seed.p.C            = 1;
      for (auto &branch : seed.production.front().tree) {
        branch.p   = exchange;
        branch.hel = ToyProtonLegHelicityMatrix(exchange.spinX2);
      }

      const auto operators = gra::spin::CanonicalPoleOperators(seed.p, exchange, exchange, true, true, true,
                                                               gra::spin::VertexContext::Auto);
      REQUIRE_FALSE(operators.empty());
      std::vector<gra::spin::LSTerm> terms;
      terms.reserve(operators.size());
      for (const auto &op : operators) {
        const double index = static_cast<double>(terms.size());
        terms.push_back({op.coupling.l, op.coupling.two_s, std::polar(0.61 + 0.013 * index, -0.37 + 0.071 * index)});
      }

      std::vector<MMatrix<std::complex<double>>> reference;
      for (const auto model : {gra::ReggeProductionModel::MP, gra::ReggeProductionModel::XP}) {
        CAPTURE(gra::ReggeProductionModelName(model));
        gra::PARAM_RES res = seed;
        PrepareToyPoleOperators(res, model, {terms}, true, true);
        const auto production = gra::rspin::Resonance(
            base, res, (model == gra::ReggeProductionModel::MP ? gra::mpom::Fusion : gra::xpom::Fusion), 1.0);
        REQUIRE(production.size() == 1);
        REQUIRE(production.front().size_row() == 16);
        REQUIRE(production.front().size_col() == static_cast<std::size_t>(2 * spin + 1));
        REQUIRE(production.front().FrobNorm2() > 0.0);
        reference.push_back(production.front());

        for (const double angle : {-1.21, 0.37, 2.04}) {
          const auto rotated_lts = RotateToyEventAroundZ(base, angle);
          const auto rotated     = gra::rspin::Resonance(
                  rotated_lts, res,
                  (model == gra::ReggeProductionModel::MP ? gra::mpom::Fusion : gra::xpom::Fusion), 1.0);
          REQUIRE(rotated.size() == 1);
          for (std::size_t in1 = 0; in1 < 2; ++in1) {
            for (std::size_t in2 = 0; in2 < 2; ++in2) {
              for (std::size_t out1 = 0; out1 < 2; ++out1) {
                for (std::size_t out2 = 0; out2 < 2; ++out2) {
                  const std::size_t row      = gra::spin::CanonicalProtonPairSpinLayout::HardRow(in1, in2, out1, out2);
                  const int         harmonic = gra::spin::ColliderSpinHalfHelicityHarmonic(
                              helicity_x2[in1], helicity_x2[in2], helicity_x2[out1], helicity_x2[out2]);
                  for (const auto &column : indices(res.production.front().pole.value().helicity.Jz_values)) {
                    const double projection = res.production.front().pole.value().helicity.Jz_values[column];
                    const std::complex<double> expected =
                        production.front()[row][column] *
                        std::exp(gra::math::zi * (static_cast<double>(harmonic) - projection) * angle);
                    CAPTURE(angle, in1, in2, out1, out2, row, harmonic, column, projection,
                            rotated.front()[row][column], expected);
                    CHECK(std::abs(rotated.front()[row][column] - expected) <=
                          5.0e-10 * std::max(1.0, std::abs(expected)));
                  }
                }
              }
            }
          }
        }
      }
      REQUIRE(reference.size() == 2);
      RequireMatrixNear(reference[0], reference[1], 5.0e-11);
    }
  }
}

// Check every physical photon RES component before decay-spin summation
TEST_CASE(
    "MP XP and GP photon resonance components obey collider azimuth "
    "covariance",
    "[gra::MRegge][MP][XP][GP][spin][photon][resonance][covariance]"
    "[regression]") {
  struct ModelCase {
    const char               *label;
    gra::ReggeProductionModel model;
  };
  const std::array<ModelCase, 3> models      = {ModelCase{"MP", gra::ReggeProductionModel::MP},
                                                ModelCase{"XP", gra::ReggeProductionModel::XP},
                                                ModelCase{"GP", gra::ReggeProductionModel::GP}};
  constexpr auto                 helicity_x2 = gra::spin::BinaryHelicityLabelsX2();
  const std::complex<double>     coupling    = std::polar(0.83, -0.41);

  gra::LORENTZSCALAR base     = MakeToyCoherentPhotonLTS();
  base.process.MP_FRAME       = "CM";
  base.process.SPINGEN        = true;
  base.process.FORWARD_NOFLIP = false;
  base.process.FORWARD_VERTEX = gra::ForwardVertexMode::HelicityResidue;
  base.process.MMAX           = 2;
  const auto  definition      = gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "photon RES covariance");
  gra::MRegge regge(base, gra::MModelTune::Load(modelfile), definition);
  const auto  param = ReggeParametersForTest(regge, base);

  for (const std::string mode : {"EPA", "QED"}) {
    for (const auto &test : models) {
      for (const bool lower_photon : {false, true}) {
        CAPTURE(mode, test.label, lower_photon);
        gra::PARAM_RES res;
        if (test.model == gra::ReggeProductionModel::MP) {
          res = MakeToyPhotoMPResonance(false);
        } else if (test.model == gra::ReggeProductionModel::XP) {
          res = MakeToyCovariantPhotoXPonance();
        } else {
          res = MakeToyPhotoGPResonance(base.process.MMAX);
        }

        if (lower_photon) {
          auto &tree = res.production.front().tree;
          std::swap(tree[0], tree[1]);
        }
        if (test.model == gra::ReggeProductionModel::GP) {
          res.production_model = gra::ReggeProductionModel::GP;
          auto &central        = res.production.front().hel;
          REQUIRE(central.alpha_ls.Size() == 1);
          central.alpha_ls.begin()->coefficient = coupling;
          const int j1x2                        = lower_photon ? 4 : 2;
          const int j2x2                        = lower_photon ? 2 : 4;
          gra::gpom::InitResonanceLS(central, res.p.spinX2 / 2, j1x2, j2x2);
        } else {
          PrepareToyPoleOperators(res, test.model, {{{0, 2, coupling}}}, false, true);
        }

        const auto evaluate = [&](const gra::LORENTZSCALAR &event) {
          std::vector<MMatrix<std::complex<double>>> production;
          if (test.model == gra::ReggeProductionModel::GP) {
            gra::gpom::AmpCache cache;
            production = gra::gpom::Resonance(event, *param, res, &cache);
          } else {
            production = gra::rspin::Resonance(
                event, res,
                (test.model == gra::ReggeProductionModel::MP ? gra::mpom::Fusion : gra::xpom::Fusion),
                param->s0);
          }
          REQUIRE(production.size() == 1);
          REQUIRE(production.front().size_row() == 16);
          const auto spin_states = static_cast<std::size_t>(res.p.spinX2 + 1);
          REQUIRE(production.front().size_col() == spin_states);
          REQUIRE(production.front().FrobNorm2() > 0.0);
          return production.front();
        };

        gra::LORENTZSCALAR reference_lts    = base;
        reference_lts.process.PHOTON_VERTEX = mode;
        const auto   reference              = evaluate(reference_lts);
        const double spin                   = 0.5 * static_cast<double>(res.p.spinX2);
        const auto   projections            = gra::spin::SpinProjections(spin);
        REQUIRE(projections.size() == reference.size_col());
        for (const double angle : {-1.21, 0.37, 2.04}) {
          const auto rotated_lts = RotateToyEventAroundZ(reference_lts, angle);
          const auto rotated     = evaluate(rotated_lts);
          for (std::size_t in1 = 0; in1 < 2; ++in1) {
            for (std::size_t in2 = 0; in2 < 2; ++in2) {
              for (std::size_t out1 = 0; out1 < 2; ++out1) {
                for (std::size_t out2 = 0; out2 < 2; ++out2) {
                  using PairLayout           = gra::spin::CanonicalProtonPairSpinLayout;
                  const std::size_t row      = PairLayout::HardRow(in1, in2, out1, out2);
                  const int         harmonic = gra::spin::ColliderSpinHalfHelicityHarmonic(
                              helicity_x2[in1], helicity_x2[in2], helicity_x2[out1], helicity_x2[out2]);
                  for (const auto &column : indices(projections)) {
                    const double               projection = projections[column];
                    const double               exponent   = (static_cast<double>(harmonic) - projection) * angle;
                    const std::complex<double> expected   = reference[row][column] * std::exp(gra::math::zi * exponent);
                    CAPTURE(angle, in1, in2, out1, out2, row, harmonic, column, projection, rotated[row][column],
                            expected);
                    const double error = std::abs(rotated[row][column] - expected);
                    const double scale = std::max(1.0, std::abs(expected));
                    CHECK(error <= 5.0e-10 * scale);
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

TEST_CASE(
    "MP and XP charged-pion continuum t/u terms are constructive at "
    "symmetric kinematics",
    "[gra::MRegge][spin][continuum][physics]") {
  for (const std::string family : {"MP", "XP"}) {
    CAPTURE(family);
    ToyHelicityProcess process;
    ConfigureToyProductionProcess(process, family, "CON", "pi+ pi-");
    REQUIRE_NOTHROW(process.InitializeProcessAmplitude());

    gra::LORENTZSCALAR lts         = process.state.lts;
    const double       beam_energy = 5.0;
    const double       proton_pz   = 4.2;
    const double       proton_pt   = 0.25;
    const double       proton_mass = lts.beam1.mass;
    const double       proton_energy =
        std::sqrt(gra::math::pow2(proton_mass) + gra::math::pow2(proton_pz) + gra::math::pow2(proton_pt));
    lts.pbeam1 = gra::M4Vec(0.0, 0.0, beam_energy, beam_energy);
    lts.pbeam2 = gra::M4Vec(0.0, 0.0, -beam_energy, beam_energy);
    lts.pfinal.resize(3);
    lts.pfinal[1]              = gra::M4Vec(proton_pt, 0.0, proton_pz, proton_energy);
    lts.pfinal[2]              = gra::M4Vec(-proton_pt, 0.0, -proton_pz, proton_energy);
    lts.pfinal[0]              = lts.pbeam1 + lts.pbeam2 - lts.pfinal[1] - lts.pfinal[2];
    const double pion_energy   = 0.5 * lts.pfinal[0].M();
    const double pion_momentum = std::sqrt(gra::math::pow2(pion_energy) - gra::math::pow2(lts.decaytree[0].p.mass));
    lts.decaytree[0].p4        = gra::M4Vec(0.0, pion_momentum, 0.0, pion_energy);
    lts.decaytree[1].p4        = gra::M4Vec(0.0, -pion_momentum, 0.0, pion_energy);
    RefreshToyDerivedKinematicsPreserveDecay(lts);

    const auto model    = family == "MP" ? gra::ReggeProductionModel::MP : gra::ReggeProductionModel::XP;
    const auto channels =
        (model == gra::ReggeProductionModel::MP ? gra::mpom::Continuum : gra::xpom::Continuum)(lts, 1.0);
    REQUIRE(channels.size() == lts.process.CONT_PRODUCTION.size());
    REQUIRE(lts.process.CONT_TU_SIGN.size() == channels.size());
    bool saw_pomeron_pair = false;
    for (const auto &channel : indices(channels)) {
      if (lts.process.CONT_PRODUCTION[channel] != std::vector<int>{995, 995}) { continue; }
      saw_pomeron_pair                              = true;
      const auto                &t                  = channels[channel].first;
      const auto                &u                  = channels[channel].second;
      const double               t_norm2            = t.FrobNorm2();
      const double               u_norm2            = u.FrobNorm2();
      const std::complex<double> overlap            = (t.Dagger() * u).Trace();
      const double               normalized_overlap = std::real(overlap) / std::sqrt(t_norm2 * u_norm2);
      CAPTURE(channel, lts.process.CONT_PRODUCTION[channel], t_norm2, u_norm2, overlap, normalized_overlap);
      REQUIRE(t_norm2 > 0.0);
      REQUIRE(u_norm2 == Approx(t_norm2).epsilon(1.0e-12));
      REQUIRE(normalized_overlap > 0.99);
    }
    REQUIRE(saw_pomeron_pair);
  }
}

// Check the extra Fermi crossing sign at symmetric baryon-pair kinematics
TEST_CASE(
    "MP and XP ppbar continuum t/u terms are destructive at symmetric "
    "kinematics",
    "[gra::MRegge][spin][continuum][baryon][C-parity][physics]"
    "[regression]") {
  for (const std::string family : {"MP", "XP"}) {
    CAPTURE(family);
    ToyHelicityProcess process;
    ConfigureToyProductionProcess(process, family, "CON", "p+ p-");
    REQUIRE_NOTHROW(process.InitializeProcessAmplitude());

    gra::LORENTZSCALAR lts         = process.state.lts;
    const double       beam_energy = 6.0;
    const double       proton_pz   = 4.2;
    const double       proton_pt   = 0.25;
    const double       proton_mass = lts.beam1.mass;
    const double       proton_energy =
        std::sqrt(gra::math::pow2(proton_mass) + gra::math::pow2(proton_pz) + gra::math::pow2(proton_pt));
    lts.pbeam1 = gra::M4Vec(0.0, 0.0, beam_energy, beam_energy);
    lts.pbeam2 = gra::M4Vec(0.0, 0.0, -beam_energy, beam_energy);
    lts.pfinal.resize(3);
    lts.pfinal[1]                = gra::M4Vec(proton_pt, 0.0, proton_pz, proton_energy);
    lts.pfinal[2]                = gra::M4Vec(-proton_pt, 0.0, -proton_pz, proton_energy);
    lts.pfinal[0]                = lts.pbeam1 + lts.pbeam2 - lts.pfinal[1] - lts.pfinal[2];
    const double baryon_energy   = 0.5 * lts.pfinal[0].M();
    const double baryon_momentum = std::sqrt(gra::math::pow2(baryon_energy) - gra::math::pow2(lts.decaytree[0].p.mass));
    lts.decaytree[0].p4          = gra::M4Vec(0.0, baryon_momentum, 0.0, baryon_energy);
    lts.decaytree[1].p4          = gra::M4Vec(0.0, -baryon_momentum, 0.0, baryon_energy);
    RefreshToyDerivedKinematicsPreserveDecay(lts);

    gra::MRegge regge(lts, process.state.model_tune,
                      gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::ContinuumTwoBody, family + "_ppbar_crossing"));
    const auto  param = ReggeParametersForTest(regge, lts);
    for (const auto &channel : indices(lts.process.CONT_PRODUCTION)) {
      const int    upper = lts.process.CONT_PRODUCTION[channel][0];
      const int    lower = lts.process.CONT_PRODUCTION[channel][1];
      const double c_t   = gra::regge::AntiparticleSign(*param, upper, lts.decaytree[0].p.pdg) *
                         gra::regge::AntiparticleSign(*param, lower, lts.decaytree[1].p.pdg);
      const double c_u = gra::regge::AntiparticleSign(*param, upper, lts.decaytree[1].p.pdg) *
                         gra::regge::AntiparticleSign(*param, lower, lts.decaytree[0].p.pdg);
      const int pair_c = gra::regge::PairCParity(*param, upper, lower);
      CAPTURE(channel, upper, lower, c_t, c_u, pair_c);
      REQUIRE(lts.process.CONT_TU_SIGN[channel] == Approx(-1.0));
      REQUIRE(lts.process.CONT_TU_SIGN[channel] * c_u / c_t == Approx(-pair_c));
    }

    const auto model    = family == "MP" ? gra::ReggeProductionModel::MP : gra::ReggeProductionModel::XP;
    const auto channels =
        (model == gra::ReggeProductionModel::MP ? gra::mpom::Continuum : gra::xpom::Continuum)(lts, 1.0);
    REQUIRE(channels.size() == lts.process.CONT_PRODUCTION.size());
    REQUIRE(lts.process.CONT_TU_SIGN.size() == channels.size());
    bool saw_pomeron_pair = false;
    for (const auto &channel : indices(channels)) {
      if (lts.process.CONT_PRODUCTION[channel] != std::vector<int>{995, 995}) { continue; }
      saw_pomeron_pair                              = true;
      const auto                &t                  = channels[channel].first;
      const auto                &u                  = channels[channel].second;
      const double               t_norm2            = t.FrobNorm2();
      const double               u_norm2            = u.FrobNorm2();
      const std::complex<double> overlap            = (t.Dagger() * u).Trace();
      const double               normalized_overlap = std::real(overlap) / std::sqrt(t_norm2 * u_norm2);
      CAPTURE(channel, lts.process.CONT_PRODUCTION[channel], t_norm2, u_norm2, overlap, normalized_overlap);
      REQUIRE(t_norm2 > 0.0);
      REQUIRE(u_norm2 == Approx(t_norm2).epsilon(1.0e-12));
      REQUIRE(normalized_overlap > 0.0);
      REQUIRE(lts.process.CONT_TU_SIGN[channel] == Approx(-1.0));
      REQUIRE(lts.process.CONT_TU_SIGN[channel] * normalized_overlap < 0.0);
    }
    REQUIRE(saw_pomeron_pair);
  }
}

// Check fixed and analytic continuum models retain their card coefficients
TEST_CASE("Continuum vertices keep direct steering coefficients",
          "[gra::MRegge][spin][continuum][normalization][physics]") {
  constexpr std::array             models       = {gra::ReggeProductionModel::MP, gra::ReggeProductionModel::XP,
                                                   gra::ReggeProductionModel::GP};
  const std::array<std::string, 5> final_states = {"p+ p-", "pi+ pi-", "K+ K-", "rho(770)0 rho(770)0",
                                                   "phi(1020)0 phi(1020)0"};

  for (const auto model : models) {
    const std::string model_name = gra::ReggeProductionModelName(model);
    for (const auto &final_state : final_states) {
      CAPTURE(model_name, final_state);
      ToyHelicityProcess process;
      ConfigureToyProductionProcess(process, model_name, "CON", final_state);
      REQUIRE_NOTHROW(process.InitializeProcessAmplitude());
      const auto &state = process.state.lts.process;
      REQUIRE_FALSE(state.CONT_PRODUCTION.empty());
      if (model == gra::ReggeProductionModel::GP) {
        REQUIRE(state.CONTINUUM_GP.size() == state.CONT_PRODUCTION.size());
      } else {
        REQUIRE(state.CONTINUUM_POLE.size() == state.CONT_PRODUCTION.size());
      }
      for (const auto &channel : indices(state.CONT_PRODUCTION)) {
        if (model == gra::ReggeProductionModel::GP) {
          REQUIRE(state.CONTINUUM_GP[channel].size() == 4);
          for (const auto &vertex : state.CONTINUUM_GP[channel]) {
            CHECK(vertex.UsesReggeDomain());
            CHECK(vertex.exchange_basis == gra::ExchangeBasisType::ReggeHelicity);
            CHECK(vertex.analytic_MMAX == state.MMAX);
            const auto pole = gra::gpom::Crossed(vertex, vertex.J, gra::M4Vec(0.3, 0.4, 0.7, 1.0), false);
            CHECK(pole.helicity.T.IsFinite());
            CHECK(pole.helicity.T.FrobNorm2() > 0.0);
          }
          continue;
        }
        REQUIRE(state.CONTINUUM_POLE[channel].size() == 4);
        for (std::size_t vertex = 0; vertex < 4; ++vertex) {
          const auto &pole = state.CONTINUUM_POLE[channel][vertex].Pole();
          RequireCanonicalPoleVertex(pole);
          std::size_t active = 0;
          for (const auto &term : pole.terms) {
            if (std::abs(term.coefficient) > 1.0e-14) {
              ++active;
              CHECK(std::isfinite(term.coefficient.real()));
              CHECK(std::isfinite(term.coefficient.imag()));
            }
          }
          CHECK(active == 1);
        }
      }
    }
  }
}

// Check the lowest DL tensor through every continuum model and input basis
TEST_CASE("MP XP GP crossed LS and helicity inputs share the lowest pole tensor",
          "[gra::MRegge][MP][XP][GP][CON][normalization][unit_coupling][parity]") {
  ModelParamRestoreGuard restore_model;
  const auto pdg = LoadedPDGTable();
  const auto pole = pdg.FindByPDG(995);
  const auto context = gra::spin::VertexContext::SubTUChannelExchange;
  const std::complex<double> coupling = std::polar(0.73, -0.28);
  for (const auto &pair : std::vector<std::pair<int, int>>{{211, -211}, {2212, -2212}, {2212, 2212}, {113, 113}, {333, 333}}) {
    const auto first = pdg.FindByPDG(pair.first);
    const auto second = pdg.FindByPDG(pair.second);
    const auto allowed = gra::spin::CanonicalPoleOperators(pole, first, second, true, true, true, context);
    const auto ls = allowed.front().coupling;
    const auto expected =
        gra::spin::PoleLSReduced(gra::spin::PreparePoleLS(pole, first, second, {{ls.l, ls.two_s, coupling}}, 1.0, true,
                                                          true, true, context, 0.0, false),
                                 1.0);
    for (const std::string family : {"MP", "XP", "GP"}) {
      const bool gp = family == "GP";
      const auto exchange = pdg.FindByPDG(gp ? 990 : 995);
      for (const std::string field : {"g_ls", "helicity"}) {
        CAPTURE(pair, family, field);
        auto rows = nlohmann::json::array();
        if (field == "g_ls") {
          for (const auto &op : allowed) {
            const bool active = op.coupling.l == ls.l && op.coupling.two_s == ls.two_s;
            auto row = nlohmann::json::array({op.coupling.l, 0.5 * op.coupling.two_s});
            if (gp) { row.push_back(0); }
            row.push_back(active ? std::abs(coupling) : 0.0);
            row.push_back(active ? std::arg(coupling) : 0.0);
            rows.push_back(row);
          }
        } else {
          std::set<std::pair<int, int>> covered;
          for (std::size_t i = 0; i < expected.size_row(); ++i) {
            for (std::size_t j = 0; j < expected.size_col(); ++j) {
              const int h1 = -first.spinX2 + 2 * static_cast<int>(i);
              const int h2 = -second.spinX2 + 2 * static_cast<int>(j);
              if (covered.contains({h1, h2})) { continue; }
              auto row = nlohmann::json::array({0.5 * h1, 0.5 * h2});
              if (gp) { row.push_back(0); }
              row.push_back(std::abs(expected[i][j]));
              row.push_back(std::arg(expected[i][j]));
              rows.push_back(row);
              covered.insert({h1, h2});
              covered.insert({-h1, -h2});
              if (gra::spin::PhysicalIdenticalPair(first, second, true, context)) {
                covered.insert({h2, h1});
                covered.insert({-h2, -h1});
              }
            }
          }
        }
        const auto pair_key = "[" + std::to_string(std::abs(pair.first)) + "," + std::to_string(std::abs(pair.second)) + "]";
        const std::string sector = first.C != 0 ? "self" : pair.second < 0 ? "opposite" : "same";
        const auto tune = WriteModifiedContinuumTune("dl_pole_" + family + field + std::to_string(pair.second), family,
            [&](auto &card) {
              card[std::to_string(exchange.pdg)][pair_key][sector] = {
                  {"basis", field == "g_ls" ? "crossed_ls" : "crossed_helicity"}, {"CP", {true, true}}, {field, rows}};
            });
        ToyHelicityProcess process;
        process.SetProcessForTest(family, "CON");
        process.SetMMAX(2);
        process.SetTuneForTest(tune);
        if (gp) {
          auto vertex = process.ProcessHelicityStructure(exchange, {first, second}, true, true, "", false, context);
          if (vertex.UsesLSCouplings()) { gra::gpom::InitCrossedLS(vertex, 2); }
          const auto actual = gra::gpom::Crossed(vertex, 2.0, gra::M4Vec(0.3, 0.4, 0.7, 1.0), false);
          const auto zero = gra::gpom::AnalyticMIndex(0, 2, "DL pole test");
          for (std::size_t i = 0; i < expected.size_row(); ++i) {
            for (std::size_t j = 0; j < expected.size_col(); ++j) {
              RequireComplexNear(actual.helicity.T[i * expected.size_col() + j][zero], expected[i][j], 2.0e-12);
            }
          }
        } else {
          const auto model = family == "MP" ? gra::ReggeProductionModel::MP : gra::ReggeProductionModel::XP;
          const auto vertex = process.ProcessPoleOperatorStructure(model, exchange, {first, second}, context);
          RequireMatrixNear(gra::spin::PoleLSReduced(vertex, 1.0), expected, 2.0e-12);
        }
      }
    }
  }
}

// Check that both continuum input bases preserve the selected C and parity constraints
TEST_CASE("MP and XP continuum helicity inputs preserve the configured symmetries",
          "[gra::MRegge][MP][XP][CON][normalization][parity][regression]") {
  ModelParamRestoreGuard restore_model;
  const auto pdg = LoadedPDGTable();
  const auto pole = pdg.FindByPDG(995);
  const auto first = pdg.FindByPDG(2212);
  const auto second = pdg.FindByPDG(-2212);
  const auto context = gra::spin::VertexContext::SubTUChannelExchange;
  for (const bool C : {false, true}) {
    for (const bool P : {false, true}) {
      const auto allowed = gra::spin::CanonicalPoleOperators(pole, first, second, true, C, P, context);
      REQUIRE_FALSE(allowed.empty());
      std::vector<gra::spin::LSTerm> terms;
      for (const auto &i : indices(allowed)) {
        const auto &ls = allowed[i].coupling;
        terms.push_back({ls.l, ls.two_s, std::polar(0.7 + 0.1 * i, 0.2 + 0.3 * i)});
      }
      const auto expected =
          gra::spin::PoleLSReduced(gra::spin::PreparePoleLS(pole, first, second, terms, 1.0, true, C, P, context), 1.0);
      auto rows = nlohmann::json::array();
      for (std::size_t i = 0; i < expected.size_row(); ++i) {
        for (std::size_t j = 0; j < expected.size_col(); ++j) {
          rows.push_back({-0.5 * first.spinX2 + i, -0.5 * second.spinX2 + j,
                          std::abs(expected[i][j]), std::arg(expected[i][j])});
        }
      }
      for (const std::string family : {"MP", "XP"}) {
        CAPTURE(family, C, P);
        const auto tune = WriteModifiedContinuumTune(
            "continuum_symmetry_" + family + std::to_string(C) + std::to_string(P), family, [&](auto &card) {
              card["995"]["[2212,2212]"]["opposite"] = {
                  {"basis", "crossed_helicity"}, {"CP", {C, P}}, {"helicity", rows}};
            });
        ToyHelicityProcess process;
        process.SetProcessForTest(family, "CON");
        process.SetTuneForTest(tune);
        const auto model = family == "MP" ? gra::ReggeProductionModel::MP : gra::ReggeProductionModel::XP;
        const auto vertex = process.ProcessPoleOperatorStructure(model, pole, {first, second}, context);
        RequireMatrixNear(gra::spin::PoleLSReduced(vertex, 1.0), expected, 2.0e-12);
      }
    }
  }
}

// Check the common DL pion couplings in the MP and XP pole bases
TEST_CASE("MP and XP charged-pion unit residues share the DL poles",
          "[gra::MRegge][spin][continuum][normalization][physics]"
          "[regression]") {
  struct ContinuumModel {
    const char *family;
    int         pomeron;
    int         reggeon;
  };
  constexpr std::array<ContinuumModel, 2> models = {{
      {"MP", 995, 9915},
      {"XP", 995, 9915},
  }};
  struct Trajectory {
    const char                *name;
    bool                       secondary;
    std::array<std::size_t, 4> vertices;
    std::size_t                count;
  };
  constexpr std::array<Trajectory, 2>                   trajectories = {{
                        {"Pomeron", false, {0, 1, 2, 3}, 4},
                        {"f2 Reggeon", true, {1, 3, 0, 0}, 2},
  }};
  std::array<std::complex<double>, trajectories.size()> reference{};
  std::array<bool, trajectories.size()>                 reference_set{};

  for (const auto &model : models) {
    CAPTURE(model.family);
    ToyHelicityProcess process;
    ConfigureToyProductionProcess(process, model.family, "CON", "pi+ pi-");
    const auto tune = WriteModifiedContinuumTune(std::string("unit_pion_") + model.family, model.family, [&](auto &card) {
      for (const int exchange : {model.pomeron, model.reggeon}) {
        for (const std::string sector : {"same", "opposite"}) {
          auto &vertex = card.at(std::to_string(exchange)).at("[211,211]").at(sector);
          vertex.erase("g");
          vertex.erase("helicity");
          vertex["basis"] = "crossed_ls";
          vertex["g_ls"] = {{2, 0, 1.0, 0.0}};
        }
      }
    });
    process.SetTuneForTest(tune);
    REQUIRE_NOTHROW(process.InitializeProcessAmplitude());

    const auto &state = process.state.lts.process;
    for (const auto &trajectory : indices(trajectories)) {
      const auto &test = trajectories[trajectory];
      CAPTURE(test.name);
      const std::vector<int> channel_id = test.secondary ? std::vector<int>{model.pomeron, model.reggeon}
                                                         : std::vector<int>{model.pomeron, model.pomeron};
      const auto channel = std::find(state.CONT_PRODUCTION.cbegin(), state.CONT_PRODUCTION.cend(), channel_id);
      REQUIRE(channel != state.CONT_PRODUCTION.cend());
      const std::size_t channel_index =
          static_cast<std::size_t>(std::distance(state.CONT_PRODUCTION.cbegin(), channel));
      for (std::size_t i = 0; i < test.count; ++i) {
        const std::size_t    vertex = test.vertices[i];
        const auto          &input   = state.CONTINUUM_POLE.at(channel_index).at(vertex).Pole();
        const auto           pole    = gra::spin::PoleLSHelicity(input, input.Lambda);
        std::complex<double> residue = pole.T[0][0];
        REQUIRE(pole.T.size_row() == 1);
        REQUIRE(pole.T.size_col() == 1);
        CAPTURE(vertex, residue);
        if (!reference_set[trajectory]) {
          reference[trajectory]     = residue;
          reference_set[trajectory] = true;
        }
        RequireComplexNear(residue, reference[trajectory], 2.0e-8);
      }
    }
  }
  CHECK(reference_set[0]);
  CHECK(reference_set[1]);
}

// Check scalar pair sewing without identifying the model dependent spin kernels
TEST_CASE(
    "TUNE0 MP and XP charged-pion DL pole residues sew with unit "
    "normalization",
    "[gra::MRegge][spin][continuum][normalization][phase][physics]"
    "[regression]") {
  struct ContinuumModel {
    const char               *family;
    gra::ReggeProductionModel model;
    int                       pomeron;
  };
  constexpr std::array<ContinuumModel, 2>        models        = {{
                    {"MP", gra::ReggeProductionModel::MP, 995},
                    {"XP", gra::ReggeProductionModel::XP, 995},
  }};
  constexpr double                               upper_phase   = 0.37;
  constexpr double                               lower_phase   = -0.61;
  const std::complex<double>                     upper_section = std::polar(1.0, upper_phase);
  const std::complex<double>                     lower_section = std::polar(1.0, lower_phase);
  const std::array<std::pair<double, double>, 3> angles        = {{{0.39, -0.62}, {0.91, 0.38}, {1.47, 1.13}}};

  for (const auto &test : models) {
    CAPTURE(test.family);
    ToyHelicityProcess process;
    ConfigureToyProductionProcess(process, test.family, "CON", "pi+ pi-");
    REQUIRE_NOTHROW(process.InitializeProcessAmplitude());

    const auto &state   = process.state.lts.process;
    const auto  channel = std::find(state.CONT_PRODUCTION.cbegin(), state.CONT_PRODUCTION.cend(),
                                    std::vector<int>{test.pomeron, test.pomeron});
    REQUIRE(channel != state.CONT_PRODUCTION.cend());
    const std::size_t channel_index = static_cast<std::size_t>(std::distance(state.CONT_PRODUCTION.cbegin(), channel));
    gra::ReggeContinuumPole pole;
    const auto             &vertices = state.CONTINUUM_POLE.at(channel_index);
    REQUIRE(vertices.size() == 4);
    pole.pole_operator = {vertices[0], vertices[1]};

    const auto reduced_residue = [&](const std::size_t leg) {
      const auto reduced =
          gra::spin::PoleLSHelicity(pole.pole_operator[leg].Pole(), pole.pole_operator[leg].Pole().Lambda);
      REQUIRE(reduced.T.size_row() == 1);
      REQUIRE(reduced.T.size_col() == 1);
      return reduced.T[0][0];
    };
    const double expected_norm2 = std::norm(reduced_residue(0)) * std::norm(reduced_residue(1));
    REQUIRE(expected_norm2 > 0.0);

    gra::ReggeContinuumPole    phased           = pole;
    const std::complex<double> expected_section = upper_section * lower_section;
    auto                       upper_pole       = phased.pole_operator[0].Pole();
    auto                       lower_pole       = phased.pole_operator[1].Pole();
    REQUIRE(upper_pole.terms.size() == 1);
    REQUIRE(lower_pole.terms.size() == 1);
    upper_pole.terms[0].coefficient *= upper_section;
    lower_pole.terms[0].coefficient *= lower_section;
    phased.pole_operator[0] = gra::spin::PoleResidue(upper_pole);
    phased.pole_operator[1] = gra::spin::PoleResidue(lower_pole);

    for (const auto &point : indices(angles)) {
      const auto [theta, phi] = angles[point];
      const auto lts          = ScalarPolePhasePointForTest(0.20, theta, phi, 1.30, 6500.0);
      const auto base    = gra::rspin::PairKernel(pole, lts.decaytree[0], lts.decaytree[1], lts.q1_in_X,
                                                  lts.q2_in_X);
      const auto rotated = gra::rspin::PairKernel(phased, lts.decaytree[0], lts.decaytree[1], lts.q1_in_X,
                                                  lts.q2_in_X);
      CAPTURE(point, theta, phi, expected_norm2, base.FrobNorm2());
      REQUIRE(base.IsFinite());
      REQUIRE(rotated.IsFinite());
      REQUIRE(base.size_row() == 1);
      REQUIRE(base.size_col() == 1);
      CHECK(base.FrobNorm2() == Approx(expected_norm2).epsilon(2.0e-10));
      REQUIRE(rotated.size_row() == base.size_row());
      REQUIRE(rotated.size_col() == base.size_col());
      const auto coherent = base + rotated;
      CHECK(coherent.FrobNorm2() == Approx(std::norm(1.0 + expected_section) * expected_norm2).epsilon(2.0e-10));
      std::size_t active = 0;
      for (std::size_t row = 0; row < base.size_row(); ++row) {
        for (std::size_t column = 0; column < base.size_col(); ++column) {
          if (std::abs(base[row][column]) <= 1.0e-14) {
            CHECK(std::abs(rotated[row][column]) <= 1.0e-12);
            continue;
          }
          ++active;
          RequireComplexNear(rotated[row][column], expected_section * base[row][column], 2.0e-10);
        }
      }
      CHECK(active > 0);
    }
  }
}

// Check the common DL baryon coupling in both fixed spin models
TEST_CASE("MP and XP ppbar Pomeron vertices share the physical pole",
          "[gra::MRegge][spin][continuum][normalization][physics]") {
  constexpr std::array models = {gra::ReggeProductionModel::MP, gra::ReggeProductionModel::XP};
  // Retain each model's exact canonical pole matrices for cross-model checks
  std::array<std::vector<gra::spin::PoleResidue>, models.size()> pole_vertices;
  for (const auto model : models) {
    const std::string family = gra::ReggeProductionModelName(model);
    CAPTURE(family);
    ToyHelicityProcess process;
    ConfigureToyProductionProcess(process, family, "CON", "p+ p-");
    REQUIRE_NOTHROW(process.InitializeProcessAmplitude());
    const int   pomeron    = 995;
    const auto &configured = process.state.lts.process;
    const auto  channel    = std::find(configured.CONT_PRODUCTION.begin(), configured.CONT_PRODUCTION.end(),
                                       std::vector<int>{pomeron, pomeron});
    REQUIRE(channel != configured.CONT_PRODUCTION.end());
    const std::size_t channel_index =
        static_cast<std::size_t>(std::distance(configured.CONT_PRODUCTION.begin(), channel));
    const auto &cache = configured.CONTINUUM_POLE.at(channel_index);
    REQUIRE(cache.size() == 4);
    pole_vertices.at(static_cast<std::size_t>(model) - 1) = cache;
    for (const auto &residue : cache) {
      const auto &vertex = residue.Pole();
      std::size_t active = 0;
      for (const auto &term : vertex.terms) {
        if (std::abs(term.coefficient) <= 1.0e-14) { continue; }
        ++active;
        CHECK(term.l == 1);
        CHECK(term.two_s == 2);
      }
      CHECK(active == 1);
    }
  }

  for (std::size_t model = 1; model < pole_vertices.size(); ++model) {
    for (std::size_t vertex = 0; vertex < 4; ++vertex) {
      CAPTURE(model, vertex);
      RequireMatrixNear(
          gra::spin::PoleLSReduced(pole_vertices[model][vertex].Pole(), pole_vertices[model][vertex].Pole().Lambda),
          gra::spin::PoleLSReduced(pole_vertices[0][vertex].Pole(), pole_vertices[0][vertex].Pole().Lambda), 2.0e-12);
    }
  }

  const auto tune       = WriteModifiedPhotoVMTune("regge_ppbar_physical_pole", [](auto &j) {
    j.at("PARAM_SOFT").at("EXCHANGE_DEF").at("P").at("trajectory_mode") = "linear";
    const std::string active = j.at("PARAM_SOFT").at("active_model").template get<std::string>();
    j.at("PARAM_SOFT").at("MODEL").at(active).at("EXCHANGE").at("P").at("alpha") = {2.0, 0.0};
  });
  const auto pole_model = gra::MModelTune::Load(tune.second);

  struct PoleEvaluation {
    std::vector<std::complex<double>>  amplitude;
    gra::MMatrix<std::complex<double>> t;
    gra::MMatrix<std::complex<double>> u;
    gra::MMatrix<std::complex<double>> decay;
    gra::MMatrix<std::complex<double>> source;
    double                             tu_sign = 0.0;
    std::vector<std::complex<double>>  upper_profile;
    std::vector<std::complex<double>>  lower_profile;
    gra::MMatrix<std::complex<double>> pair_profile;
    std::complex<double>               t_kernel   = 0.0;
    std::complex<double>               u_kernel   = 0.0;
    double                             t_crossing = 0.0;
    double                             u_crossing = 0.0;
  };
  const auto evaluate = [&](const gra::ReggeProductionModel model, const bool spingen) {
    const std::string  family = gra::ReggeProductionModelName(model);
    ToyHelicityProcess process;
    ConfigureToyProductionProcess(process, family, "CON", "p+ p-");
    REQUIRE_NOTHROW(process.InitializeProcessAmplitude());
    const auto &configured = process.state.lts.process;
    const int   pomeron    = 995;
    const auto  channel    = std::find(configured.CONT_PRODUCTION.begin(), configured.CONT_PRODUCTION.end(),
                                       std::vector<int>{pomeron, pomeron});
    REQUIRE(channel != configured.CONT_PRODUCTION.end());
    const std::size_t channel_index =
        static_cast<std::size_t>(std::distance(configured.CONT_PRODUCTION.begin(), channel));

    gra::LORENTZSCALAR lts = DirectCentralPairLTSForTest(2212, -2212);
    RefreshToyDerivedKinematicsPreserveDecay(lts);
    lts.process                     = configured;
    lts.process.SPINGEN             = spingen;
    lts.process.CONT_PRODUCTION     = {configured.CONT_PRODUCTION.at(channel_index)};
    lts.process.CONT_PRODUCTIONTREE = {configured.CONT_PRODUCTIONTREE.at(channel_index)};
    lts.process.CONTINUUM_POLE      = {configured.CONTINUUM_POLE.at(channel_index)};
    lts.process.CONT_TU_SIGN        = {configured.CONT_TU_SIGN.at(channel_index)};
    lts.hamp.Configure(process.state.lts.hamp.metadata);
    gra::MRegge regge(lts, pole_model,
                      gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::ContinuumTwoBody, family + "_ppbar_pole"));
    const auto  param      = ReggeParametersForTest(regge, lts);
    const auto  production =
        (model == gra::ReggeProductionModel::MP ? gra::mpom::Continuum : gra::xpom::Continuum)(lts, param->s0);
    REQUIRE(production.size() == 1);
    const auto decay =
        gra::spin::ContinuumDecayMatrix(lts, model == gra::ReggeProductionModel::MP ? lts.process.MP_FRAME : "CM");
    const auto post_decay_production =
        (model == gra::ReggeProductionModel::MP ? gra::mpom::Continuum : gra::xpom::Continuum)(lts, param->s0);
    REQUIRE(post_decay_production.size() == 1);
    RequireMatrixNear(post_decay_production.front().first, production.front().first, 2.0e-12);
    RequireMatrixNear(post_decay_production.front().second, production.front().second, 2.0e-12);
    const auto upper_state    = gra::ResolveForwardLegState(lts, gra::ForwardBeamLeg::Upper);
    const auto lower_state    = gra::ResolveForwardLegState(lts, gra::ForwardBeamLeg::Lower);
    const int  upper_exchange = lts.process.CONT_PRODUCTION.front()[0];
    const int  lower_exchange = lts.process.CONT_PRODUCTION.front()[1];
    const auto upper_profile  = gra::LegProfile(lts, upper_state, regge, *param, upper_exchange);
    const auto lower_profile  = gra::LegProfile(lts, lower_state, regge, *param, lower_exchange);
    REQUIRE(upper_profile.size() == 1);
    REQUIRE(lower_profile.size() == 1);
    const auto pair_profile = gra::PairSources(upper_profile, lower_profile, regge.SoftModelHandle()->GoodWalker(),
                                               gra::ProtonSpinLayout(lts), 1.0);
    REQUIRE(pair_profile.size() == 1);
    const auto t_kernel = gra::LegKernel(lts, upper_state, regge, *param, upper_exchange, lts.ss[1][3], false, 0) *
                          gra::LegKernel(lts, lower_state, regge, *param, lower_exchange, lts.ss[2][4], false, 0);
    const auto u_kernel = gra::LegKernel(lts, upper_state, regge, *param, upper_exchange, lts.ss[1][4], false, 0) *
                          gra::LegKernel(lts, lower_state, regge, *param, lower_exchange, lts.ss[2][3], false, 0);
    const double t_crossing = gra::regge::AntiparticleSign(*param, upper_exchange, lts.decaytree[0].p.pdg) *
                              gra::regge::AntiparticleSign(*param, lower_exchange, lts.decaytree[1].p.pdg);
    const double u_crossing = gra::regge::AntiparticleSign(*param, upper_exchange, lts.decaytree[1].p.pdg) *
                              gra::regge::AntiparticleSign(*param, lower_exchange, lts.decaytree[0].p.pdg);
    TestReggeCon(regge, lts, model);
    REQUIRE(gra::SquaredNorm(lts.hamp) > 0.0);
    REQUIRE(lts.proton_good_walker.has_value());
    REQUIRE(lts.proton_good_walker->components.size() == 1);
    return PoleEvaluation{std::vector<std::complex<double>>(lts.hamp.begin(), lts.hamp.end()),
                          production.front().first,
                          production.front().second,
                          decay,
                          lts.proton_good_walker->components.front().source,
                          lts.process.CONT_TU_SIGN.front(),
                          upper_profile.front().nonflip,
                          lower_profile.front().nonflip,
                          pair_profile.front().amplitude,
                          t_kernel,
                          u_kernel,
                          t_crossing,
                          u_crossing};
  };

  for (const bool spingen : {false, true}) {
    CAPTURE(spingen);
    const auto mp = evaluate(gra::ReggeProductionModel::MP, spingen);
    const auto xp = evaluate(gra::ReggeProductionModel::XP, spingen);
    INFO("XP physical-pole comparison");
    RequireMatrixNear(xp.t, mp.t, 2.0e-10);
    RequireMatrixNear(xp.u, mp.u, 2.0e-10);
    RequireMatrixNear(xp.decay, mp.decay, 2.0e-10);
    RequireMatrixNear(xp.source, mp.source, 2.0e-10);
    RequireVectorNear(xp.amplitude, mp.amplitude, 2.0e-10);
  }
}

// Check longitudinal boosts at fixed invariants through the public continuum amplitude
TEST_CASE("MP XP and GP continua preserve longitudinal boosts",
          "[gra::MRegge][spin][continuum][boost][physics][regression]") {
  for (const auto spin : {gra::ReggeProductionModel::MP, gra::ReggeProductionModel::XP,
                         gra::ReggeProductionModel::GP}) {
    const std::string family = gra::ReggeProductionModelName(spin);
    ToyHelicityProcess process;
    ConfigureToyProductionProcess(process, family, "CON", "pi+ pi-");
    REQUIRE_NOTHROW(process.InitializeProcessAmplitude());
    auto base = AsymmetricF2PhasePointForTest(0.91, -0.38);
    RefreshToyDerivedKinematicsPreserveDecay(base);
    base.process = process.state.lts.process;
    base.hamp.Configure(process.state.lts.hamp.metadata);
    const auto evaluate = [&](gra::LORENTZSCALAR lts) {
      gra::MRegge regge(lts, process.state.model_tune,
                       gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::ContinuumTwoBody, family + "_CON"));
      return TestReggeCon(regge, lts, spin);
    };
    const double reference = evaluate(base);
    REQUIRE(reference > 0.0);
    for (const double rapidity : {-1.6, -0.8, 0.8, 1.6}) {
      const auto shifted = BoostToyEventAlongZ(base, rapidity);
      CHECK(shifted.s == Approx(base.s).epsilon(2.0e-12));
      CHECK(shifted.t1 == Approx(base.t1).epsilon(2.0e-12));
      CHECK(shifted.t2 == Approx(base.t2).epsilon(2.0e-12));
      const double actual = evaluate(shifted);
      CAPTURE(family, rapidity, actual, reference);
      CHECK(actual == Approx(reference).epsilon(2.0e-9));
    }
  }
}

// Check the equal pion rapidity region through the public continuum amplitude
TEST_CASE("MP XP and GP charged-pion continuum has no equal rapidity node",
          "[gra::MRegge][spin][continuum][physics][regression]") {
  constexpr std::array<double, 5> polar_offsets = {-0.35, -0.18, 0.0, 0.18, 0.35};

  for (const auto spin : {gra::ReggeProductionModel::MP, gra::ReggeProductionModel::XP,
                         gra::ReggeProductionModel::GP}) {
    const std::string family = gra::ReggeProductionModelName(spin);
    CAPTURE(family);
    ToyHelicityProcess process;
    ConfigureToyProductionProcess(process, family, "CON", "pi+ pi-");
    REQUIRE_NOTHROW(process.InitializeProcessAmplitude());

    gra::LORENTZSCALAR base          = process.state.lts;
    const double       beam_pz       = 5.0;
    const double       proton_pz     = 4.2;
    const double       proton_mass   = base.beam1.mass;
    const double       beam_energy   = std::sqrt(gra::math::pow2(proton_mass) + gra::math::pow2(beam_pz));
    const double       proton_energy = std::sqrt(gra::math::pow2(proton_mass) + gra::math::pow2(proton_pz));
    base.pbeam1                      = gra::M4Vec(0.0, 0.0, beam_pz, beam_energy);
    base.pbeam2                      = gra::M4Vec(0.0, 0.0, -beam_pz, beam_energy);
    base.pfinal.resize(3);
    base.pfinal[1]             = gra::M4Vec(0.0, 0.0, proton_pz, proton_energy);
    base.pfinal[2]             = gra::M4Vec(0.0, 0.0, -proton_pz, proton_energy);
    base.pfinal[0]             = base.pbeam1 + base.pbeam2 - base.pfinal[1] - base.pfinal[2];
    const double pion_energy   = 0.5 * base.pfinal[0].M();
    const double pion_momentum = std::sqrt(gra::math::pow2(pion_energy) - gra::math::pow2(base.decaytree[0].p.mass));
    std::array<double, polar_offsets.size()> amp2{};
    std::array<double, polar_offsets.size()> rapidity_gap{};
    for (const auto &i : indices(polar_offsets)) {
      gra::LORENTZSCALAR lts   = base;
      const double       theta = 0.5 * gra::math::PI + polar_offsets[i];
      const double       py    = pion_momentum * std::sin(theta);
      const double       pz    = pion_momentum * std::cos(theta);
      lts.decaytree[0].p4      = gra::M4Vec(0.0, py, pz, pion_energy);
      lts.decaytree[1].p4      = gra::M4Vec(0.0, -py, -pz, pion_energy);
      RefreshToyDerivedKinematicsPreserveDecay(lts);

      gra::MRegge regge(lts, process.state.model_tune,
                        gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::ContinuumTwoBody, family + "_CON"));
      amp2[i]         = TestReggeCon(regge, lts, spin);
      rapidity_gap[i] = lts.decaytree[0].p4.Rap() - lts.decaytree[1].p4.Rap();
      CAPTURE(i, polar_offsets[i], rapidity_gap[i], amp2[i]);
      REQUIRE(std::isfinite(amp2[i]));
      REQUIRE(amp2[i] > 0.0);
    }

    const std::size_t center = polar_offsets.size() / 2;
    CHECK(rapidity_gap[center] == Approx(0.0).margin(1.0e-14));
    for (std::size_t i = 0; i < center; ++i) {
      const std::size_t mirror = polar_offsets.size() - 1 - i;
      CAPTURE(i, mirror, amp2[i], amp2[mirror]);
      CHECK(rapidity_gap[i] == Approx(-rapidity_gap[mirror]).epsilon(1.0e-12));
      CHECK(amp2[i] == Approx(amp2[mirror]).epsilon(1.0e-10));
    }
    const double adjacent = 0.5 * (amp2[center - 1] + amp2[center + 1]);
    CAPTURE(amp2[center], adjacent);
    CHECK(amp2[center] > 0.25 * adjacent);
  }
}

// Check card-loaded resonance production and the corresponding public amplitude
TEST_CASE("TUNE0 MP and XP resonance trees use the matched spin-two Pomeron",
          "[gra::MRegge][spin][resonance][physics][regression]") {
  struct ResonanceSpinCase {
    const char               *family;
    gra::ReggeProductionModel spin;
    int                       pomeron_spinX2;
  };
  constexpr std::array<ResonanceSpinCase, 2> cases = {
      {{"MP", gra::ReggeProductionModel::MP, 4}, {"XP", gra::ReggeProductionModel::XP, 4}}};

  for (const auto &test : cases) {
    CAPTURE(test.family, test.pomeron_spinX2);
    ToyHelicityProcess process;
    ConfigureToyProductionProcess(process, test.family, "RES", "pi+ pi-");
    const auto f0 = gra::resonance::Read("RES/f0_500.json", process.state.random, gra::ParseReggeProductionModel(test.family));
    process.SetResonances({{"f0_500", f0}});
    REQUIRE_NOTHROW(process.InitializeProcessAmplitude());
    gra::PARAM_RES configured = process.GetResonances().at("f0_500");
    REQUIRE_FALSE(configured.production.empty());
    for (const auto &channel : indices(configured.production)) {
      REQUIRE(configured.production[channel].tree.size() == 2);
      for (const auto &leg : configured.production[channel].tree) {
        REQUIRE(leg.p.pdg == 995);
        CHECK(leg.p.C == 1);
        CHECK(leg.p.spinX2 == test.pomeron_spinX2);
      }
      REQUIRE(configured.production[channel].pole.has_value());
      const auto reduced =
          gra::spin::PoleLSReduced(configured.production[channel].pole.value(), 0.5 * configured.p.mass);
      CHECK(reduced.FrobNorm2() > 0.0);
    }

    gra::LORENTZSCALAR lts     = ScalarPolePhasePointForTest(0.20, 0.91, -0.38, configured.p.mass);
    lts.process.MP_FRAME       = process.state.lts.process.MP_FRAME;
    lts.process.FORWARD_VERTEX = process.state.lts.process.FORWARD_VERTEX;
    lts.process.FORWARD_NOFLIP = process.state.lts.process.FORWARD_NOFLIP;
    gra::MRegge  regge(lts, process.state.model_tune,
                       gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Resonance, std::string(test.family) + "_RES"));
    const double amp2 = TestReggeRes(regge, lts, configured, test.spin);
    CAPTURE(amp2);
    CHECK(std::isfinite(amp2));
    CHECK(amp2 > 0.0);
  }
}

// Check that resonance derivative factors are selected per Regge model
TEST_CASE("PARAM_REGGE selects resonance derivative factors by model", "[gra::MRegge][spin][resonance][derivative]") {
  const bool derivative_factor = GENERATE(false, true);
  for (const std::string family : {"MP", "XP", "GP"}) {
    CAPTURE(family, derivative_factor);
    const auto tune = WriteModifiedPhotoVMTune("derivative_" + family + std::to_string(derivative_factor), [&](auto &card) {
      auto &factors = card.at("PARAM_REGGE").at("DERIVATIVE_FACTOR");
      for (auto &value : factors) { value = !derivative_factor; }
      factors.at(family) = derivative_factor;
    });
    ToyHelicityProcess process;
    ConfigureToyProductionProcess(process, family, "RES", "pi+ pi-");
    process.SetModelTune(gra::MModelTune::Load(tune.second));
    auto f0 = gra::resonance::Read("RES/f0_500.json", process.state.random, gra::ParseReggeProductionModel(family));
    if (family == "MP") { SetMPFusion(f0); }
    process.SetResonances({{"f0_500", f0}});
    REQUIRE_NOTHROW(process.InitializeProcessAmplitude());
    CHECK(process.state.lts.process.DERIVATIVE_FACTOR == derivative_factor);

    const auto &configured = process.GetResonances().at("f0_500");
    if (family != "GP") {
      REQUIRE_FALSE(configured.production.empty());
      for (const auto &production : configured.production) {
        REQUIRE(production.pole.has_value());
        CHECK(production.pole->derivative_factor == derivative_factor);
      }
    }
  }
}

// Check that an enabled MP derivative factor has a declared momentum scale
TEST_CASE("MP derivative factor requires a positive channel Lambda", "[gra::MRegge][spin][resonance][derivative]") {
  const auto tune = WriteModifiedPhotoVMTune(
      "mp_derivative_factor", [](auto &card) { card.at("PARAM_REGGE").at("DERIVATIVE_FACTOR").at("MP") = true; });
  const auto model = gra::MModelTune::Load(tune.second);

  ToyHelicityProcess missing;
  ConfigureToyProductionProcess(missing, "MP", "RES", "pi+ pi-");
  missing.SetModelTune(model);
  auto f0 = gra::resonance::Read("RES/f0_500.json", missing.state.random, gra::ReggeProductionModel::MP);
  SetMPFusion(f0);
  auto       no_scale                 = f0;
  no_scale.MP.channels.front().Lambda = 0.0;
  missing.SetResonances({{"f0_500", no_scale}});
  REQUIRE_THROWS_AS(missing.InitializeProcessAmplitude(), std::invalid_argument);

  ToyHelicityProcess declared;
  ConfigureToyProductionProcess(declared, "MP", "RES", "pi+ pi-");
  declared.SetModelTune(model);
  auto scaled                       = f0;
  scaled.MP.channels.front().Lambda = 0.73;
  declared.SetResonances({{"f0_500", scaled}});
  REQUIRE_NOTHROW(declared.InitializeProcessAmplitude());
  const auto &vertex = declared.GetResonances().at("f0_500").production.front().pole.value();
  CHECK(vertex.derivative_factor);
  CHECK(vertex.Lambda == Approx(0.73));
}

// Check every spin-parity class through the initialized production caches
TEST_CASE("MP XP and GP resonance classes have finite configured residues",
          "[gra::MRegge][spin][resonance][normalization][physics]") {
  struct ResonanceClass {
    const char *label;
    const char *card;
    const char *decay;
    int         spin_x2;
    int         parity;
    int         c_parity;
  };
  constexpr std::array<ResonanceClass, 10> resonances = {{
      {"eta", "RES/eta.json", "gamma gamma", 0, -1, 1},
      {"f0_500", "RES/f0_500.json", "pi+ pi-", 0, 1, 1},
      {"rho_770", "RES/rho_770.json", "pi+ pi-", 2, -1, -1},
      {"f1_1420", "RES/f1_1420.json", "K*(892)+ K-", 2, 1, 1},
      {"eta2_1645", "RES/eta2_1645.json", "a(2)(1320)0 pi0", 4, -1, 1},
      {"f2_1270", "RES/f2_1270.json", "pi+ pi-", 4, 1, 1},
      {"rho3_1690", "RES/rho3_1690.json", "pi+ pi-", 6, -1, -1},
      {"f4_2300", "RES/f4_2300.json", "pi+ pi-", 8, 1, 1},
      {"rho5_2350", "RES/rho5_2350.json", "pi+ pi-", 10, -1, -1},
      {"f6_2510", "RES/f6_2510.json", "pi+ pi-", 12, 1, 1},
  }};
  constexpr std::array                     models     = {gra::ReggeProductionModel::MP, gra::ReggeProductionModel::XP,
                                                         gra::ReggeProductionModel::GP};

  for (const auto &test : resonances) {
    for (const auto model : models) {
      const std::string model_name = gra::ReggeProductionModelName(model);
      CAPTURE(test.label, model_name, test.spin_x2, test.parity, test.c_parity);
      ToyHelicityProcess process;
      ConfigureToyProductionProcess(process, model_name, "RES", test.decay);
      const auto input = gra::resonance::Read(test.card, process.state.random, gra::ParseReggeProductionModel(model_name));
      process.SetResonances({{test.label, input}});
      REQUIRE_NOTHROW(process.InitializeProcessAmplitude());

      const auto &configured = process.GetResonances().at(test.label);
      CHECK(configured.p.spinX2 == test.spin_x2);
      CHECK(configured.p.P == test.parity);
      CHECK(configured.p.C == test.c_parity);
      REQUIRE_FALSE(configured.production.empty());
      REQUIRE(configured.hel_decay.T.IsFinite());
      CHECK(configured.hel_decay.J == Approx(0.5 * static_cast<double>(test.spin_x2)).margin(1.0e-12));
      CHECK(gra::spin::ReducedHelicityNorm2(configured.hel_decay) == Approx(1.0).epsilon(2.0e-12));

      for (const auto &channel : indices(configured.production)) {
        if (model == gra::ReggeProductionModel::GP) {
          REQUIRE(channel < configured.production.size());
          const auto  &hel      = configured.production[channel].hel;
          const double residue2 = hel.UsesLSCouplings() ? hel.alpha_ls.Norm2() : hel.T.MaskedSquaredNorm(hel.T_set);
          CAPTURE(channel, residue2);
          CHECK(std::isfinite(residue2));
          CHECK(residue2 > 0.0);
          continue;
        }
        if (model == gra::ReggeProductionModel::MP && !configured.UsesUnrestrictedSpinBasis()) {
          auto event = MakeToyProductionLTSAsymmetric(0.3, 4.2, -4.6);
          event.process = process.state.lts.process;
          const auto amplitude = gra::mpom::Resonance(event, configured, 1.0).at(channel);
          REQUIRE(amplitude.size_row() > 0);
          REQUIRE(amplitude.IsFinite());
          const auto reference = gra::rspin::Resonance(event, configured, gra::mpom::Fusion, 1.0).at(channel);
          RequireMatrixNear(amplitude, reference * configured.MP.filter, 2e-12);
          continue;
        }
        const double density = ConfiguredPoleReferenceDensity(configured, model, channel);
        CAPTURE(channel, density);
        CHECK(std::isfinite(density));
        CHECK(density > 0.0);
      }
    }
  }
}

// Check complete GP steering tables before sparse runtime cache construction
TEST_CASE("TUNE0 GP LS and helicity tables build sparse runtime terms",
          "[gra::MRegge][GP][spin][resonance][sparse][regression]") {
  struct SparseCase {
    const char *label;
    const char *card;
    const char *decay;
  };
  constexpr std::array<SparseCase, 20> resonances = {{
      {"chi_c0", "RES/chi_c0.json", "J/psi(1S)0 gamma"},
      {"chi_c1", "RES/chi_c1.json", "J/psi(1S)0 gamma"},
      {"chi_c2", "RES/chi_c2.json", "J/psi(1S)0 gamma"},
      {"eta", "RES/eta.json", "gamma gamma"},
      {"eta2_1645", "RES/eta2_1645.json", "a(2)(1320)0 pi0"},
      {"eta_prime", "RES/eta_prime.json", "gamma gamma"},
      {"f0_500", "RES/f0_500.json", "pi+ pi-"},
      {"f0_980", "RES/f0_980.json", "pi+ pi-"},
      {"f0_1500", "RES/f0_1500.json", "pi+ pi-"},
      {"f0_1710", "RES/f0_1710.json", "pi+ pi-"},
      {"f0_1710_neg_P", "RES/f0_1710_neg_P.json", "rho(770)0 rho(770)0"},
      {"f1_1420", "RES/f1_1420.json", "K*(892)+ K-"},
      {"f2_1270", "RES/f2_1270.json", "pi+ pi-"},
      {"f2_1525", "RES/f2_1525.json", "pi+ pi-"},
      {"f2_1950", "RES/f2_1950.json", "pi+ pi-"},
      {"f2_2150", "RES/f2_2150.json", "pi+ pi-"},
      {"f4_2300", "RES/f4_2300.json", "pi+ pi-"},
      {"f6_2510", "RES/f6_2510.json", "pi+ pi-"},
      {"phi_1020_odd", "RES/phi_1020_odd.json", "K+ K-"},
      {"rho_770_odd", "RES/rho_770_odd.json", "pi+ pi-"},
  }};

  for (const auto &test : resonances) {
    CAPTURE(test.label, test.card, test.decay);
    ToyHelicityProcess process;
    ConfigureToyProductionProcess(process, "GP", "RES", test.decay);
    const auto   input        = gra::resonance::Read(test.card, process.state.random, gra::ReggeProductionModel::GP);
    const double coupling_min = process.state.model_tune->Global().coupling_min;
    process.SetResonances({{test.label, input}});
    REQUIRE_NOTHROW(process.InitializeProcessAmplitude());
    const auto &configured = process.GetResonances().at(test.label);
    REQUIRE(configured.GP.channels.front().g_ls.Size() == input.GP.channels.front().g_ls.Size());
    REQUIRE_FALSE(configured.production.empty());
    for (const auto &production : configured.production) {
      const auto &runtime = production.hel;
      const auto &channel = configured.GP.channels.front();
      if (channel.basis == gra::ReggeVertexBasis::LS) {
        const std::size_t active = static_cast<std::size_t>(
            std::count_if(channel.g_ls.cbegin(), channel.g_ls.cend(),
                          [coupling_min](const auto &term) { return std::abs(term.coefficient) > coupling_min; }));
        REQUIRE(runtime.UsesLSCouplings());
        REQUIRE(runtime.alpha_ls.Size() == active);
        REQUIRE(runtime.gp_orbital.terms.size() == active);
        REQUIRE(runtime.gp_orbital.su2.size() == active * runtime.gp_orbital.nmu);
        for (const auto &term : channel.g_ls) {
          const bool retained = std::abs(term.coefficient) > coupling_min;
          CHECK(runtime.alpha_ls.Contains(term.l, term.two_s) == retained);
        }
      } else {
        REQUIRE(channel.basis == gra::ReggeVertexBasis::Helicity);
        CHECK_FALSE(runtime.UsesLSCouplings());
        CHECK(runtime.UsesHelicityCouplings());
        CHECK(runtime.T.IsFinite());
        CHECK_FALSE(runtime.T_active.empty());
      }
    }
  }
}

// Check complete fixed-spin tables before sparse pole-operator construction
TEST_CASE("MP and XP explicit LS tables validate before sparse runtime terms",
          "[gra::MRegge][MP][XP][spin][resonance][sparse][regression]") {
  for (const std::string family : {"MP", "XP"}) {
    CAPTURE(family);
    ToyHelicityProcess process;
    ConfigureToyProductionProcess(process, family, "RES", "pi+ pi-");
    auto       input = gra::resonance::Read("RES/f0_980.json", process.state.random, gra::ParseReggeProductionModel(family));
    const auto pole  = process.state.lts.PDG.FindByPDG(995);
    const auto allowed =
        gra::spin::CanonicalPoleOperators(input.p, pole, pole, true, true, true, gra::spin::VertexContext::Auto);
    REQUIRE(allowed.size() > 1);

    const auto configure = [&](gra::PARAM_RES &res, const bool omit) {
      if (family == "MP") { SetMPFusion(res); }
      auto &channel    = family == "MP" ? res.MP.channels.front() : res.XP.channels.front();
      channel.exchange = {995, 995};
      channel.basis    = gra::ReggeVertexBasis::LS;
      channel.Lambda   = 1.0;
      channel.g_ls.Clear();
      const std::size_t count = allowed.size() - (omit ? 1U : 0U);
      for (std::size_t i = 0; i < count; ++i) {
        channel.g_ls.Set(allowed[i].coupling.l, allowed[i].coupling.two_s, i == 0 ? 1.0 : 0.0);
      }
    };

    configure(input, false);
    process.SetResonances({{"f0_980", input}});
    REQUIRE_NOTHROW(process.InitializeProcessAmplitude());
    const auto &configured = process.GetResonances().at("f0_980");
    REQUIRE(configured.production.size() == 1);
    const auto &runtime = configured.production.front().pole.value();
    REQUIRE(runtime.terms.size() == 1);
    REQUIRE(runtime.reduced_basis.size() == runtime.terms.size());
    CHECK(std::abs(runtime.terms.front().coefficient) > process.state.model_tune->Global().coupling_min);
    const auto &channel = family == "MP" ? configured.MP.channels.front() : configured.XP.channels.front();
    CHECK(channel.g_ls.Size() == allowed.size());

    ToyHelicityProcess missing;
    ConfigureToyProductionProcess(missing, family, "RES", "pi+ pi-");
    auto omitted = gra::resonance::Read("RES/f0_980.json", missing.state.random, gra::ParseReggeProductionModel(family));
    configure(omitted, true);
    missing.SetResonances({{"f0_980", omitted}});
    REQUIRE_THROWS_AS(missing.InitializeProcessAmplitude(), std::invalid_argument);
  }
}

// Check unit scalar vertices at a matched spin 2 Pomeron pole
TEST_CASE(
    "MP XP GP and TP unit f0 resonance couplings share matched pole "
    "normalizations",
    "[gra::MRegge][MTensorPomeron][spin][resonance][normalization]"
    "[physics][unit_coupling]") {
  ModelParamRestoreGuard restore_model;
  gra::MODELPARAM = "TUNE0";

  const nlohmann::json scalar_ls       = {{0, 0, 1.0, 0.0}, {2, 2, 0.0, 0.0}, {4, 4, 0.0, 0.0}};
  const nlohmann::json scalar_helicity = {{-2, -2, 1.0, 0.0}, {-1, -1, 1.0, gra::math::PI}, {0, 0, 1.0, 0.0}};
  const std::array<std::string, 2>                  families = {"MP", "XP"};
  const std::array<std::string, 2>                  bases    = {"g_ls", "helicity"};
  std::array<gra::MMatrix<std::complex<double>>, 4> fixed;

  std::size_t cache = 0;
  for (const auto &family : families) {
    for (const auto &basis : bases) {
      CAPTURE(family, basis);
      ToyHelicityProcess process;
      ConfigureToyProductionProcess(process, family, "RES", "pi+ pi-");
      auto  f0         = gra::resonance::Read("RES/f0_980.json", process.state.random, gra::ParseReggeProductionModel(family));
      if (family == "MP") { SetMPFusion(f0); }
      auto &channel    = family == "MP" ? f0.MP.channels.front() : f0.XP.channels.front();
      channel.exchange = {995, 995};
      if (basis == "g_ls") {
        SetChannelLS(channel, scalar_ls);
      } else {
        SetChannelHelicity(channel, scalar_helicity);
      }
      process.SetResonances({{"f0_980", f0}});
      REQUIRE_NOTHROW(process.InitializeProcessAmplitude());
      const auto &configured = process.GetResonances().at("f0_980");
      REQUIRE(configured.production.size() == 1);
      fixed[cache++] = gra::spin::PoleLSReduced(configured.production.front().pole.value(), 0.0);
    }
  }
  REQUIRE(cache == fixed.size());
  for (std::size_t i = 1; i < fixed.size(); ++i) { RequireMatrixNear(fixed[i], fixed.front(), 2.0e-12); }

  ToyHelicityProcess gp_process;
  ConfigureToyProductionProcess(gp_process, "GP", "RES", "pi+ pi-");
  auto  gp_f0         = gra::resonance::Read("RES/f0_980.json", gp_process.state.random, gra::ReggeProductionModel::GP);
  auto &gp_channel    = gp_f0.GP.channels.front();
  gp_channel.exchange = {990, 990};
  SetChannelLS(gp_channel, scalar_ls);
  gp_process.SetResonances({{"f0_980", gp_f0}});
  REQUIRE_NOTHROW(gp_process.InitializeProcessAmplitude());
  const auto &gp_configured = gp_process.GetResonances().at("f0_980");
  REQUIRE(gp_configured.production.size() == 1);
  auto gp_pole = gp_configured.production.front().hel;
  REQUIRE(gp_pole.UsesLSCouplings());
  REQUIRE(gp_pole.alpha_ls.Size() == 1);
  RequireComplexNear(gp_pole.alpha_ls.At(0, 0), 1.0, 1.0e-14);
  REQUIRE(gp_pole.gp_orbital.su2.size() == 1);
  CHECK(gp_pole.gp_orbital.su2.front() == Approx(std::sqrt(5.0)));

  gra::LORENTZSCALAR tensor_lts;
  tensor_lts.PDG = LoadedPDGTable();
  gra::MTensorPomeron                tensor(tensor_lts, gra::MModelTune::Load(modelfile),
                                            gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
  const double                       exchange_mass = 0.9;
  const double                       momentum      = 1.0e-8 * exchange_mass;
  const double                       energy = std::sqrt(gra::math::pow2(exchange_mass) + gra::math::pow2(momentum));
  const gra::M4Vec                   q1(0.0, 0.0, momentum, energy);
  const gra::M4Vec                   q2(0.0, 0.0, -momentum, energy);
  const auto                         vertex = tensor.iG_PPS_0();
  gra::MMatrix<std::complex<double>> tp(5, 5, 0.0);
  for (int h1 = -2; h1 <= 2; ++h1) {
    const auto eps1 = tensor.EpsMassiveSpin2(q1, h1);
    for (int h2 = -2; h2 <= 2; ++h2) {
      const auto   eps2             = tensor.EpsMassiveSpin2(q2, h2);
      const double second_leg_phase = gra::spin::JacobWickSecondLegReversalPhase(2.0, h2);
      for (const auto &mu : tensor.LI) {
        for (const auto &nu : tensor.LI) {
          for (const auto &kappa : tensor.LI) {
            for (const auto &lambda : tensor.LI) {
              tp[static_cast<std::size_t>(h1 + 2)][static_cast<std::size_t>(h2 + 2)] +=
                  std::conj(eps1(mu, nu)) * (second_leg_phase * std::conj(eps2(kappa, lambda))) *
                  vertex(mu, nu, kappa, lambda);
            }
          }
        }
      }
    }
  }
  RequireMatrixNear(tp, fixed.front() * (2.0 * gra::math::zi), 2.0e-12);
}

TEST_CASE("MP and XP direct helicity preserve an arbitrary common phase",
          "[gra::MRegge][MP][XP][helicity][phase][physics]") {
  ModelParamRestoreGuard restore_model;
  gra::MODELPARAM = "TUNE0";

  for (const std::string family : {"MP", "XP"}) {
    CAPTURE(family);
    // Build one scalar pole coefficient with a selected common helicity phase
    const auto build = [&](const double phase) {
      ToyHelicityProcess process;
      ConfigureToyProductionProcess(process, family, "RES", "pi+ pi-");
      auto                 f0       = gra::resonance::Read("RES/f0_980.json", process.state.random, gra::ParseReggeProductionModel(family));
      if (family == "MP") { SetMPFusion(f0); }
      const nlohmann::json helicity = {{-2, -2, 1.0, phase}, {-1, -1, 1.0, -gra::math::PI + phase}, {0, 0, 1.0, phase}};
      auto                &channel  = family == "MP" ? f0.MP.channels.front() : f0.XP.channels.front();
      channel.exchange              = {995, 995};
      SetChannelHelicity(channel, helicity);
      process.SetResonances({{"f0_980", f0}});
      process.InitializeProcessAmplitude();
      const auto &configured = process.GetResonances().at("f0_980");
      REQUIRE(configured.production.size() == 1);
      return gra::spin::PoleLSReduced(configured.production.front().pole.value(), 0.0);
    };

    constexpr double phase     = 0.37;
    const auto       reference = build(0.0);
    auto             expected  = reference;
    expected *= std::polar(1.0, phase);
    RequireMatrixNear(build(phase), expected, 2.0e-12);
  }
}

// Check the analytic LS event evaluator at the canonical spin two pole
TEST_CASE("GP fusion LS runtime equals the canonical fixed pole",
          "[gra::MRegge][GP][spin][resonance][normalization][unit_coupling]") {
  gra::LORENTZSCALAR lts        = MakeToyReggeLTSAsymmetric(0.3, 4.2, -4.6);
  lts.process.MMAX              = 2;
  lts.process.SPINGEN           = true;
  lts.process.DERIVATIVE_FACTOR = false;
  gra::MRegge regge(lts, gra::MModelTune::Load(modelfile),
                    gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));
  const auto  param = ReggeParametersForTest(regge, lts);

  gra::PARAM_RES ls          = MakeToyScalarGPResonance(lts.process.MMAX);
  const auto     analytic    = lts.PDG.FindByPDG(9910);
  ls.production[0].tree[0].p = analytic;
  ls.production[0].tree[1].p = analytic;

  const std::size_t trajectory = gra::regge::TrajectoryIndex(*param, 9910);
  const auto       &soft       = param->soft_model->Exchange(param->exchanges.at(trajectory).soft_exchange);
  REQUIRE(soft.TrajectoryMode() == gra::SoftTrajectoryMode::Linear);
  REQUIRE(soft.AlphaPrime() > 0.0);
  const double pole_transfer = (2.0 - soft.Alpha0()) / soft.AlphaPrime();
  lts.t1                     = pole_transfer;
  lts.t2                     = pole_transfer;

  const auto &pole = gra::regge::PoleRepresentative(*param, lts.PDG, 9910);
  const auto     fixed  = gra::spin::PoleLSReduced(gra::spin::PreparePoleLS(ls.p, pole, pole, {{0, 0, 1.0}}, 1.0), 0.0);
  gra::PARAM_RES direct = ls;
  auto          &hel    = direct.production.front().hel;
  hel.coupling_basis    = gra::CouplingBasis::Helicity;
  hel.alpha_ls.Clear();
  hel.gp_orbital = {};
  hel.T          = fixed;
  hel.T_set      = fixed.Transform([](const auto &value) { return gra::math::abs2(value) > 1.0e-24; });
  gra::PruneHelicityCouplings(hel, 0.0, "GP direct pole comparison");

  gra::gpom::AmpCache ls_cache;
  gra::gpom::AmpCache direct_cache;
  const auto          ls_runtime     = gra::gpom::Resonance(lts, *param, ls, &ls_cache);
  const auto          direct_runtime = gra::gpom::Resonance(lts, *param, direct, &direct_cache);
  REQUIRE(ls_runtime.size() == 1);
  REQUIRE(direct_runtime.size() == 1);
  RequireMatrixNear(ls_runtime.front(), direct_runtime.front(), 2.0e-12);
}

// Check numerical steering through the real immutable model reader
TEST_CASE("GP trajectory floor is validated during initialization", "[gra::MRegge][GP][freeze][params]") {
  const auto tune = WriteModifiedPhotoVMTune("gp_freeze", {});
  const auto path = std::filesystem::path(tune.first) / "NUMERICS.json";
  auto numerics = nlohmann::json::parse(gra::aux::GetInputData(path.string()));
  for (const nlohmann::json& value : {nlohmann::json(0.0), nlohmann::json(0.25), nlohmann::json(-0.1),
                                    nlohmann::json(-0.5), nlohmann::json(-5.0), nlohmann::json(1.0), nlohmann::json(1.1),
                                    nlohmann::json(nullptr), nlohmann::json("zero")}) {
    CAPTURE(value);
    numerics["NUMERICS_REGGE"]["GP_alpha_min"] = value;
    std::ofstream(path) << numerics.dump();
    const auto model = gra::MModelTune::Load(tune.second);
    if (value.is_number() && value.get<double>() > -0.5 && value.get<double>() <= 1.0) {
      CHECK(gra::regge::ReadParam({211, -211}, LoadedPDGTable(), *model).gp_alpha_min == Approx(value.get<double>()));
    } else {
      CHECK_THROWS_AS(gra::regge::ReadParam({211, -211}, LoadedPDGTable(), *model), std::invalid_argument);
    }
  }
}

// Compare complex scalar and tensor vertices at and below the angular trajectory floor
TEST_CASE("GP angular trajectories freeze continuously with collider symmetries", "[gra::MRegge][GP][freeze][spin]") {
  auto lts = MakeToyReggeLTSAsymmetric(0.31, 4.2, -4.6);
  lts.process.MMAX = 2;
  lts.process.SPINGEN = true;
  lts.process.FORWARD_NOFLIP = true;
  lts.process.DERIVATIVE_FACTOR = false;
  auto param = gra::regge::ReadParam({211, -211}, lts.PDG, *gra::MModelTune::Load(modelfile));
  for (const double floor : {-0.25, 0.0, 0.25}) {
    param.gp_alpha_min = floor;
    for (auto res : {MakeToyScalarGPResonance(2), MakeToyTensorGPResonance(2)}) {
      CAPTURE(floor, res.p.spinX2);
      for (auto &branch : res.production.front().tree) { branch.p = lts.PDG.FindByPDG(9910); }
      // Use a complex LS coupling to retain sensitivity to phase changes
      res.production.front().hel.alpha_ls.At(0, res.p.spinX2) *= std::polar(1.0, 0.37);
      const auto evaluate = [&](const gra::LORENTZSCALAR &event) {
        return gra::gpom::Resonance(event, param, res, nullptr).front();
      };
      lts.t1 = ReggeTransferForAlpha(param, 9910, floor);
      lts.t2 = ReggeTransferForAlpha(param, 9910, 0.875);
      const auto reference = evaluate(lts);
      REQUIRE(reference.FrobNorm2() > 0.0);
      for (const double alpha : {floor - 1.0e-7, -3.25, -4336.125}) {
        CAPTURE(alpha);
        lts.t1 = ReggeTransferForAlpha(param, 9910, alpha);
        CHECK(gra::regge::Alpha(param, 9910, lts.t1) == Approx(alpha));
        RequireMatrixNear(evaluate(lts), reference, 2.0e-11);
      }
      auto rotated = RotateToyEventAroundZ(lts, 0.63);
      rotated.t1 = lts.t1;
      rotated.t2 = lts.t2;
      CHECK(evaluate(rotated).FrobNorm2() == Approx(reference.FrobNorm2()).epsilon(2.0e-10));
      auto reflected = ReflectToyEventInXZ(lts);
      reflected.t1 = lts.t1;
      reflected.t2 = lts.t2;
      CHECK(evaluate(reflected).FrobNorm2() == Approx(reference.FrobNorm2()).epsilon(2.0e-10));
      auto exchanged = BeamExchangeMirrorWithDecay(lts);
      exchanged.t1 = lts.t2;
      exchanged.t2 = lts.t1;
      CHECK(evaluate(exchanged).FrobNorm2() == Approx(reference.FrobNorm2()).epsilon(2.0e-10));
      if (res.p.spinX2 == 0) { RequireMatrixNear(evaluate(rotated), reference, 2.0e-11); }
      // Gamma square roots can give a cusp at integer trajectory spin
      double previous = std::numeric_limits<double>::infinity();
      for (const double offset : {1.0e-4, 1.0e-6, 1.0e-8}) {
        lts.t1 = ReggeTransferForAlpha(param, 9910, floor + offset);
        const double distance = std::sqrt((evaluate(lts) - reference).FrobNorm2() / reference.FrobNorm2());
        CAPTURE(offset, distance, previous);
        CHECK(distance < 0.2 * previous + 1.0e-10);
        previous = distance;
      }
      CHECK(previous < 1.0e-3);
    }
  }
}

// Exercise the GP continuum ladder entry with negative raw trajectory spins
TEST_CASE("GP ladder angular trajectories use the same floor", "[gra::MRegge][GP][freeze][continuum]") {
  auto event = MakeToyScalarContinuumLTSAsymmetric(0.3, 4.2, -4.6);
  auto param = gra::regge::ReadParam({211, -211}, event.PDG, *gra::MModelTune::Load(modelfile));
  gra::HELMatrix vertex;
  gra::gpom::InitCrossed(vertex, 0.0, 0.0, 0, "GP frozen crossed LS");
  vertex.coupling_basis = gra::CouplingBasis::LS;
  vertex.m_ls.resize(1);
  vertex.m_ls[0].Set(2, 0, {1.0, 0.4});
  gra::gpom::InitCrossedLS(vertex, 2);
  gra::ReggeContinuumPole pole;
  pole.gp_vertex = {vertex, vertex};
  for (const double floor : {-0.25, 0.0, 0.25}) {
    param.gp_alpha_min = floor;
    const auto reference = gra::gpom::PairKernel(pole, param, event.decaytree[0], event.decaytree[1], floor, 0.875);
    for (const double alpha : {-0.5, -3.25, -4336.125}) {
      const auto frozen = gra::gpom::PairKernel(pole, param, event.decaytree[0], event.decaytree[1], alpha, 0.875);
      RequireMatrixNear(frozen, reference, 2.0e-12);
    }
  }
}

// Check that the scalar block uses the common Gamma continuation off pole
TEST_CASE("GP scalar fusion uses the uniform Gamma continuation off pole",
          "[gra::MRegge][GP][spin][resonance][normalization][regression]") {
  gra::LORENTZSCALAR lts        = MakeToyReggeLTSAsymmetric(0.3, 4.2, -4.6);
  lts.process.MMAX              = 2;
  lts.process.SPINGEN           = true;
  lts.process.DERIVATIVE_FACTOR = false;
  gra::MRegge regge(lts, gra::MModelTune::Load(modelfile),
                    gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));
  const auto  param = ReggeParametersForTest(regge, lts);

  constexpr double alpha1 = 0.91;
  constexpr double alpha2 = 1.07;
  lts.t1                  = ReggeTransferForAlpha(*param, 990, alpha1);
  lts.t2                  = ReggeTransferForAlpha(*param, 990, alpha2);

  const gra::PARAM_RES res = MakeToyScalarGPResonance(lts.process.MMAX);
  gra::gpom::AmpCache  cache;
  const auto           production = gra::gpom::Resonance(lts, *param, res, &cache);
  REQUIRE(production.size() == 1);
  REQUIRE(cache.sources.size() == 1);
  const auto &source = cache.sources.front();
  const auto  block  = std::find_if(source.regge_cg.cbegin(), source.regge_cg.cend(),
                                    [](const auto &candidate) { return candidate.two_s == 0; });
  REQUIRE(block != source.regge_cg.cend());
  CHECK(block->reflection_phase == 1);

  const std::complex<double> section = std::exp(-2.0 * gra::math::zi * gra::math::PI * (alpha1 - alpha2));
  for (int m = -source.basis.MMAX; m <= source.basis.MMAX; ++m) {
    const std::size_t i     = gra::gpom::AnalyticMIndex(m, source.basis.MMAX, "off-pole GP scalar m");
    const auto       &value = block->coefficient[i * source.basis.nm + i];
    REQUIRE(value.has_value());
    const auto direct =
        gra::wigner::CGRegge(alpha1, alpha2, static_cast<double>(m), -static_cast<double>(m), 0.0, 0.0);
    const auto reflected =
        gra::wigner::CGRegge(alpha1, alpha2, -static_cast<double>(m), static_cast<double>(m), 0.0, 0.0);
    RequireComplexNear(*value, section * 0.5 * (direct + reflected), 2.0e-12);
  }

  const std::size_t zero  = gra::gpom::AnalyticMIndex(0, source.basis.MMAX, "off-pole GP scalar m zero");
  const auto       &value = block->coefficient[zero * source.basis.nm + zero];
  REQUIRE(value.has_value());
  const double               mean = 0.5 * (alpha1 + alpha2);
  const std::complex<double> scalar_mean =
      section * std::exp(gra::math::zi * gra::math::PI * mean) / std::sqrt(2.0 * mean + 1.0);
  CHECK(std::abs(*value - scalar_mean) > 1.0e-4);
}

// Check the completed Regge-spin block away from its integer pole
TEST_CASE("GP fusion LS preserves reflected Regge spin off pole",
          "[gra::MRegge][GP][spin][parity][resonance][regression]") {
  gra::LORENTZSCALAR lts        = MakeToyReggeLTSAsymmetric(0.3, 4.2, -4.6);
  lts.process.MMAX              = 2;
  lts.process.SPINGEN           = true;
  lts.process.DERIVATIVE_FACTOR = false;
  gra::MRegge regge(lts, gra::MModelTune::Load(modelfile),
                    gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));
  const auto  param = ReggeParametersForTest(regge, lts);

  constexpr double alpha    = 1.02;
  const double     transfer = ReggeTransferForAlpha(*param, 990, alpha);
  lts.t1                    = transfer;
  lts.t2                    = transfer;

  const gra::PARAM_RES res = MakeToyTensorGPResonance(lts.process.MMAX);
  gra::gpom::AmpCache  cache;
  const auto           production = gra::gpom::Resonance(lts, *param, res, &cache);
  REQUIRE(production.size() == 1);
  REQUIRE(cache.sources.size() == 1);
  const auto &source = cache.sources.front();
  CHECK(source.basis.alpha1 == Approx(alpha).epsilon(1.0e-12));
  CHECK(source.basis.alpha2 == Approx(alpha).epsilon(1.0e-12));
  const auto block = std::find_if(source.regge_cg.cbegin(), source.regge_cg.cend(),
                                  [](const auto &candidate) { return candidate.two_s == 4; });
  REQUIRE(block != source.regge_cg.cend());
  CHECK(block->reflection_phase == 1);

  std::size_t nonzero_pairs = 0;
  for (int m1 = -source.basis.MMAX; m1 <= source.basis.MMAX; ++m1) {
    for (int m2 = -source.basis.MMAX; m2 <= source.basis.MMAX; ++m2) {
      if (std::abs(m1 - m2) > 2) { continue; }
      const std::size_t i1    = gra::gpom::AnalyticMIndex(m1, source.basis.MMAX, "off-pole GP parity m1");
      const std::size_t i2    = gra::gpom::AnalyticMIndex(m2, source.basis.MMAX, "off-pole GP parity m2");
      const std::size_t p1    = gra::gpom::AnalyticMIndex(-m1, source.basis.MMAX, "off-pole GP parity reflected m1");
      const std::size_t p2    = gra::gpom::AnalyticMIndex(-m2, source.basis.MMAX, "off-pole GP parity reflected m2");
      const auto       &value = block->coefficient[i1 * source.basis.nm + i2];
      const auto       &reflected = block->coefficient[p1 * source.basis.nm + p2];
      REQUIRE(value.has_value());
      REQUIRE(reflected.has_value());
      CAPTURE(m1, m2, *value, *reflected);
      RequireComplexNear(*reflected, *value, 2.0e-12);
      nonzero_pairs += std::abs(*value) > 1.0e-10 ? 1U : 0U;
    }
  }
  CHECK(nonzero_pairs > 0);

  constexpr double   pole_alpha    = 2.0;
  gra::LORENTZSCALAR pole_lts      = lts;
  const double       pole_transfer = ReggeTransferForAlpha(*param, 990, pole_alpha);
  pole_lts.t1                      = pole_transfer;
  pole_lts.t2                      = pole_transfer;
  gra::gpom::AmpCache pole_cache;
  const auto          pole_production = gra::gpom::Resonance(pole_lts, *param, res, &pole_cache);
  REQUIRE(pole_production.size() == 1);
  REQUIRE(pole_cache.sources.size() == 1);
  const auto &pole_source = pole_cache.sources.front();
  CHECK(pole_source.basis.alpha1 == Approx(pole_alpha).epsilon(1.0e-12));
  CHECK(pole_source.basis.alpha2 == Approx(pole_alpha).epsilon(1.0e-12));
  const auto pole_block = std::find_if(pole_source.regge_cg.cbegin(), pole_source.regge_cg.cend(),
                                       [](const auto &candidate) { return candidate.two_s == 4; });
  REQUIRE(pole_block != pole_source.regge_cg.cend());
  const std::size_t probe_1    = gra::gpom::AnalyticMIndex(1, source.basis.MMAX, "off-pole GP alpha probe m1");
  const std::size_t probe_2    = gra::gpom::AnalyticMIndex(0, source.basis.MMAX, "off-pole GP alpha probe m2");
  const auto       &probe      = block->coefficient[probe_1 * source.basis.nm + probe_2];
  const auto       &pole_probe = pole_block->coefficient[probe_1 * pole_source.basis.nm + probe_2];
  REQUIRE(probe.has_value());
  REQUIRE(pole_probe.has_value());
  CHECK(std::abs(*pole_probe - *probe) > 1.0e-4);

  gra::PARAM_RES odd = MakeToyTensorGPResonance(lts.process.MMAX);
  odd.p.P            = -1;
  auto &odd_hel      = odd.production.front().hel;
  odd_hel.alpha_ls.Clear();
  odd_hel.alpha_ls.Set(1, 2, 1.0);
  gra::gpom::InitResonanceLS(odd_hel, odd.p.spinX2 / 2, 4, 4);
  gra::gpom::AmpCache odd_cache;
  const auto          odd_production = gra::gpom::Resonance(lts, *param, odd, &odd_cache);
  REQUIRE(odd_production.size() == 1);
  REQUIRE(odd_cache.sources.size() == 1);
  const auto &odd_source = odd_cache.sources.front();
  const auto  odd_block  = std::find_if(odd_source.regge_cg.cbegin(), odd_source.regge_cg.cend(),
                                        [](const auto &candidate) { return candidate.two_s == 2; });
  REQUIRE(odd_block != odd_source.regge_cg.cend());
  CHECK(odd_block->reflection_phase == -1);
  for (int m1 = -odd_source.basis.MMAX; m1 <= odd_source.basis.MMAX; ++m1) {
    for (int m2 = -odd_source.basis.MMAX; m2 <= odd_source.basis.MMAX; ++m2) {
      if (std::abs(m1 - m2) > 1) { continue; }
      const std::size_t i1        = gra::gpom::AnalyticMIndex(m1, odd_source.basis.MMAX, "odd GP parity m1");
      const std::size_t i2        = gra::gpom::AnalyticMIndex(m2, odd_source.basis.MMAX, "odd GP parity m2");
      const std::size_t p1        = gra::gpom::AnalyticMIndex(-m1, odd_source.basis.MMAX, "odd GP parity reflected m1");
      const std::size_t p2        = gra::gpom::AnalyticMIndex(-m2, odd_source.basis.MMAX, "odd GP parity reflected m2");
      const auto       &value     = odd_block->coefficient[i1 * odd_source.basis.nm + i2];
      const auto       &reflected = odd_block->coefficient[p1 * odd_source.basis.nm + p2];
      REQUIRE(value.has_value());
      REQUIRE(reflected.has_value());
      CAPTURE(m1, m2, *value, *reflected);
      RequireComplexNear(*reflected, -*value, 2.0e-12);
    }
  }

  gra::gpom::AmpCache odd_pole_cache;
  const auto          odd_pole_production = gra::gpom::Resonance(pole_lts, *param, odd, &odd_pole_cache);
  REQUIRE(odd_pole_production.size() == 1);
  REQUIRE(odd_pole_cache.sources.size() == 1);
  const auto &odd_pole_source = odd_pole_cache.sources.front();
  const auto  odd_pole_block  = std::find_if(odd_pole_source.regge_cg.cbegin(), odd_pole_source.regge_cg.cend(),
                                             [](const auto &candidate) { return candidate.two_s == 2; });
  REQUIRE(odd_pole_block != odd_pole_source.regge_cg.cend());
  const auto &odd_probe      = odd_block->coefficient[probe_1 * odd_source.basis.nm + probe_2];
  const auto &odd_pole_probe = odd_pole_block->coefficient[probe_1 * odd_pole_source.basis.nm + probe_2];
  REQUIRE(odd_probe.has_value());
  REQUIRE(odd_pole_probe.has_value());
  CHECK(std::abs(*odd_pole_probe - *odd_probe) > 1.0e-4);
}

// Check identical analytic legs after interchanging unequal trajectories
TEST_CASE("GP fusion LS preserves identical-leg exchange off pole",
          "[gra::MRegge][GP][spin][parity][beam-exchange][resonance]"
          "[regression]") {
  gra::LORENTZSCALAR lts        = MakeToyReggeLTSAsymmetric(0.3, 4.2, -4.6);
  lts.process.MMAX              = 2;
  lts.process.SPINGEN           = true;
  lts.process.DERIVATIVE_FACTOR = false;
  gra::MRegge      regge(lts, gra::MModelTune::Load(modelfile),
                         gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));
  const auto       param     = ReggeParametersForTest(regge, lts);
  constexpr double alpha1    = 1.02;
  constexpr double alpha2    = 0.82;
  const double     transfer1 = ReggeTransferForAlpha(*param, 990, alpha1);
  const double     transfer2 = ReggeTransferForAlpha(*param, 990, alpha2);
  lts.t1                     = transfer1;
  lts.t2                     = transfer2;
  const gra::PARAM_RES res   = MakeToyTensorGPResonance(lts.process.MMAX);
  gra::gpom::AmpCache  direct_cache;
  const auto           direct = gra::gpom::Resonance(lts, *param, res, &direct_cache);
  REQUIRE(direct.size() == 1);
  REQUIRE(direct_cache.sources.size() == 1);
  const auto &direct_source = direct_cache.sources.front();
  CHECK(direct_source.basis.alpha1 == Approx(alpha1).epsilon(1.0e-12));
  CHECK(direct_source.basis.alpha2 == Approx(alpha2).epsilon(1.0e-12));

  gra::LORENTZSCALAR exchanged_lts = lts;
  exchanged_lts.t1                 = transfer2;
  exchanged_lts.t2                 = transfer1;
  gra::gpom::AmpCache exchanged_cache;
  const auto          exchanged = gra::gpom::Resonance(exchanged_lts, *param, res, &exchanged_cache);
  REQUIRE(exchanged.size() == 1);
  REQUIRE(exchanged_cache.sources.size() == 1);
  const auto &exchanged_source = exchanged_cache.sources.front();
  CHECK(exchanged_source.basis.alpha1 == Approx(alpha2).epsilon(1.0e-12));
  CHECK(exchanged_source.basis.alpha2 == Approx(alpha1).epsilon(1.0e-12));

  const auto direct_block    = std::find_if(direct_source.regge_cg.cbegin(), direct_source.regge_cg.cend(),
                                            [](const auto &candidate) { return candidate.two_s == 4; });
  const auto exchanged_block = std::find_if(exchanged_source.regge_cg.cbegin(), exchanged_source.regge_cg.cend(),
                                            [](const auto &candidate) { return candidate.two_s == 4; });
  REQUIRE(direct_block != direct_source.regge_cg.cend());
  REQUIRE(exchanged_block != exchanged_source.regge_cg.cend());
  CHECK(direct_block->reflection_phase == 1);
  CHECK(exchanged_block->reflection_phase == 1);

  std::size_t nonzero_pairs = 0;
  for (int m1 = -direct_source.basis.MMAX; m1 <= direct_source.basis.MMAX; ++m1) {
    for (int m2 = -direct_source.basis.MMAX; m2 <= direct_source.basis.MMAX; ++m2) {
      if (std::abs(m1 - m2) > 2) { continue; }
      const std::size_t i1 = gra::gpom::AnalyticMIndex(m1, direct_source.basis.MMAX, "GP leg exchange direct m1");
      const std::size_t i2 = gra::gpom::AnalyticMIndex(m2, direct_source.basis.MMAX, "GP leg exchange direct m2");
      const std::size_t x1 = gra::gpom::AnalyticMIndex(m2, exchanged_source.basis.MMAX, "GP leg exchange swapped m1");
      const std::size_t x2 = gra::gpom::AnalyticMIndex(m1, exchanged_source.basis.MMAX, "GP leg exchange swapped m2");
      const auto       &value           = direct_block->coefficient[i1 * direct_source.basis.nm + i2];
      const auto       &exchanged_value = exchanged_block->coefficient[x1 * exchanged_source.basis.nm + x2];
      REQUIRE(value.has_value());
      REQUIRE(exchanged_value.has_value());
      CAPTURE(m1, m2, *value, *exchanged_value);
      RequireComplexNear(*exchanged_value, *value, 2.0e-12);
      nonzero_pairs += std::abs(*value) > 1.0e-10 ? 1U : 0U;
    }
  }
  CHECK(nonzero_pairs > 0);

  gra::LORENTZSCALAR pole_lts      = lts;
  const double       pole_transfer = ReggeTransferForAlpha(*param, 990, 2.0);
  pole_lts.t1                      = pole_transfer;
  pole_lts.t2                      = pole_transfer;
  gra::gpom::AmpCache pole_cache;
  const auto          pole = gra::gpom::Resonance(pole_lts, *param, res, &pole_cache);
  REQUIRE(pole.size() == 1);
  REQUIRE(pole_cache.sources.size() == 1);
  const auto &pole_source = pole_cache.sources.front();
  const auto  pole_block  = std::find_if(pole_source.regge_cg.cbegin(), pole_source.regge_cg.cend(),
                                         [](const auto &candidate) { return candidate.two_s == 4; });
  REQUIRE(pole_block != pole_source.regge_cg.cend());
  const std::size_t probe1 = gra::gpom::AnalyticMIndex(1, direct_source.basis.MMAX, "GP leg exchange alpha probe m1");
  const std::size_t probe2 = gra::gpom::AnalyticMIndex(0, direct_source.basis.MMAX, "GP leg exchange alpha probe m2");
  const auto       &probe  = direct_block->coefficient[probe1 * direct_source.basis.nm + probe2];
  const auto       &pole_probe = pole_block->coefficient[probe1 * pole_source.basis.nm + probe2];
  REQUIRE(probe.has_value());
  REQUIRE(pole_probe.has_value());
  REQUIRE(std::abs(*probe) > 1.0e-10);
  CHECK(std::abs(std::abs(*probe) - std::abs(*pole_probe)) > 1.0e-4);

  gra::PARAM_RES odd = MakeToyTensorGPResonance(lts.process.MMAX);
  odd.p.P            = -1;
  auto &odd_hel      = odd.production.front().hel;
  odd_hel.alpha_ls.Clear();
  odd_hel.alpha_ls.Set(1, 2, 1.0);
  gra::gpom::InitResonanceLS(odd_hel, odd.p.spinX2 / 2, 4, 4);

  gra::gpom::AmpCache odd_direct_cache;
  gra::gpom::AmpCache odd_exchanged_cache;
  const auto          odd_direct    = gra::gpom::Resonance(lts, *param, odd, &odd_direct_cache);
  const auto          odd_exchanged = gra::gpom::Resonance(exchanged_lts, *param, odd, &odd_exchanged_cache);
  REQUIRE(odd_direct.size() == 1);
  REQUIRE(odd_exchanged.size() == 1);
  REQUIRE(odd_direct_cache.sources.size() == 1);
  REQUIRE(odd_exchanged_cache.sources.size() == 1);
  const auto &odd_direct_source    = odd_direct_cache.sources.front();
  const auto &odd_exchanged_source = odd_exchanged_cache.sources.front();
  const auto  odd_direct_block = std::find_if(odd_direct_source.regge_cg.cbegin(), odd_direct_source.regge_cg.cend(),
                                              [](const auto &candidate) { return candidate.two_s == 2; });
  const auto  odd_exchanged_block =
      std::find_if(odd_exchanged_source.regge_cg.cbegin(), odd_exchanged_source.regge_cg.cend(),
                   [](const auto &candidate) { return candidate.two_s == 2; });
  REQUIRE(odd_direct_block != odd_direct_source.regge_cg.cend());
  REQUIRE(odd_exchanged_block != odd_exchanged_source.regge_cg.cend());
  CHECK(odd_direct_block->reflection_phase == -1);
  CHECK(odd_exchanged_block->reflection_phase == -1);

  std::size_t odd_nonzero_pairs = 0;
  for (int m1 = -odd_direct_source.basis.MMAX; m1 <= odd_direct_source.basis.MMAX; ++m1) {
    for (int m2 = -odd_direct_source.basis.MMAX; m2 <= odd_direct_source.basis.MMAX; ++m2) {
      if (std::abs(m1 - m2) > 1) { continue; }
      const std::size_t i1 =
          gra::gpom::AnalyticMIndex(m1, odd_direct_source.basis.MMAX, "odd GP leg exchange direct m1");
      const std::size_t i2 =
          gra::gpom::AnalyticMIndex(m2, odd_direct_source.basis.MMAX, "odd GP leg exchange direct m2");
      const std::size_t x1 =
          gra::gpom::AnalyticMIndex(m2, odd_exchanged_source.basis.MMAX, "odd GP leg exchange swapped m1");
      const std::size_t x2 =
          gra::gpom::AnalyticMIndex(m1, odd_exchanged_source.basis.MMAX, "odd GP leg exchange swapped m2");
      const auto &value           = odd_direct_block->coefficient[i1 * odd_direct_source.basis.nm + i2];
      const auto &exchanged_value = odd_exchanged_block->coefficient[x1 * odd_exchanged_source.basis.nm + x2];
      REQUIRE(value.has_value());
      REQUIRE(exchanged_value.has_value());
      CAPTURE(m1, m2, *value, *exchanged_value);
      RequireComplexNear(*exchanged_value, *value, 2.0e-12);
      odd_nonzero_pairs += std::abs(*value) > 1.0e-10 ? 1U : 0U;
    }
  }
  CHECK(odd_nonzero_pairs > 0);
}

// Check physical spin one and mixed exchange poles in every coupled-spin sector
TEST_CASE("GP fusion LS pole reflection covers spin one and mixed exchanges",
          "[gra::MRegge][GP][spin][parity][beam-exchange][resonance]"
          "[normalization][regression]") {
  struct PoleCase {
    const char *label;
    int         upper_pdg;
    int         lower_pdg;
    int         J;
    int         parity;
    int         c_parity;
    std::size_t l;
    std::size_t s;
    bool        identical;
  };
  constexpr std::array<PoleCase, 6> cases = {{
      {"P O S1", 990, 9990, 1, -1, -1, 0, 1, false},
      {"P O S2", 990, 9990, 1, -1, -1, 2, 2, false},
      {"P O S3", 990, 9990, 1, -1, -1, 2, 3, false},
      {"O O S0", 9990, 9990, 0, 1, 1, 0, 0, true},
      {"O O S1", 9990, 9990, 0, -1, 1, 1, 1, true},
      {"O O S2", 9990, 9990, 0, 1, 1, 2, 2, true},
  }};

  gra::LORENTZSCALAR base        = MakeToyReggeLTSAsymmetric(0.3, 4.2, -4.6);
  base.process.MMAX              = 2;
  base.process.SPINGEN           = true;
  base.process.DERIVATIVE_FACTOR = false;
  gra::MRegge regge(base, gra::MModelTune::Load(modelfile),
                    gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));
  const auto  param = ReggeParametersForTest(regge, base);

  for (const auto &test : cases) {
    CAPTURE(test.label, test.upper_pdg, test.lower_pdg, test.J, test.parity, test.l, test.s);
    // Only nonzero orbital CG factors require an exchange-spin continuation
    const auto active_mu = [&](const int mu) {
      return std::abs(mu) <= test.J &&
             !gra::math::IsZero(gra::wigner::CG(static_cast<double>(test.l), static_cast<double>(test.s), 0.0,
                                                static_cast<double>(mu), static_cast<double>(test.J),
                                                static_cast<double>(mu)));
    };
    const auto  upper      = base.PDG.FindByPDG(test.upper_pdg);
    const auto  lower      = base.PDG.FindByPDG(test.lower_pdg);
    const auto &upper_pole = gra::regge::PoleRepresentative(*param, base.PDG, test.upper_pdg);
    const auto &lower_pole = gra::regge::PoleRepresentative(*param, base.PDG, test.lower_pdg);
    const int   j1x2       = upper_pole.spinX2;
    const int   j2x2       = lower_pole.spinX2;
    REQUIRE(j1x2 >= 0);
    REQUIRE(j2x2 >= 0);

    gra::PARAM_RES ls =
        test.J == 0 ? MakeToyScalarGPResonance(base.process.MMAX) : MakeToyGPResonance(base.process.MMAX);
    ls.p.spinX2                     = 2 * test.J;
    ls.p.P                          = test.parity;
    ls.p.C                          = test.c_parity;
    ls.production.front().tree[0].p = upper;
    ls.production.front().tree[1].p = lower;

    gra::HELMatrix hel;
    gra::gpom::InitHelicity(hel, 0.5 * static_cast<double>(j1x2), 0.5 * static_cast<double>(j2x2), base.process.MMAX,
                            "GP exchange pole sweep");
    hel.BR             = 1.0;
    hel.P_symmetry     = true;
    hel.C_symmetry     = true;
    hel.coupling_basis = gra::CouplingBasis::LS;
    hel.alpha_ls.Set(test.l, 2 * test.s, 1.0);
    gra::gpom::InitResonanceLS(hel, test.J, j1x2, j2x2);
    ls.production.front().hel = hel;

    const auto vertex = gra::spin::PreparePoleLS(ls.p, upper_pole, lower_pole, {{test.l, 2 * test.s, 1.0}}, 1.0, true,
                                                 true, true, gra::spin::VertexContext::Auto, 0.0, false);
    REQUIRE(gra::spin::PoleLSReduced(vertex, 0.37).FrobNorm2() > 0.0);

    gra::LORENTZSCALAR pole_lts = base;
    pole_lts.t1                 = ReggeTransferForAlpha(*param, test.upper_pdg, 0.5 * static_cast<double>(j1x2));
    pole_lts.t2                 = ReggeTransferForAlpha(*param, test.lower_pdg, 0.5 * static_cast<double>(j2x2));
    gra::gpom::AmpCache pole_cache;
    const auto          pole = gra::gpom::Resonance(pole_lts, *param, ls, &pole_cache);
    REQUIRE(pole.size() == 1);
    REQUIRE(pole_cache.sources.size() == 1);
    const auto &pole_source      = pole_cache.sources.front();
    const auto &pole_block       = GPSpinBlock(pole_source, 2 * test.s);
    const int   reflection_phase = ((j1x2 + j2x2 - static_cast<int>(2 * test.s)) / 2) % 2 == 0 ? 1 : -1;
    CHECK(pole_block.reflection_phase == reflection_phase);

    std::size_t pole_nonzero = 0;
    const int   j1           = j1x2 / 2;
    const int   j2           = j2x2 / 2;
    for (int m1 = -j1; m1 <= j1; ++m1) {
      for (int m2 = -j2; m2 <= j2; ++m2) {
        const int mu = m1 - m2;
        if (!active_mu(mu)) { continue; }
        const auto &value = GPSpinValue(pole_source, pole_block, m1, m2);
        const auto  expected =
            gra::wigner::CG(static_cast<double>(j1), static_cast<double>(j2), static_cast<double>(m1),
                            -static_cast<double>(m2), static_cast<double>(test.s), static_cast<double>(mu));
        RequireComplexNear(value, expected, 2.0e-12);
        pole_nonzero += std::abs(value) > 1.0e-12 ? 1U : 0U;
      }
    }
    CHECK(pole_nonzero > 0);

    gra::LORENTZSCALAR off_lts = base;
    off_lts.t1                 = ReggeTransferForAlpha(*param, test.upper_pdg, 1.02);
    off_lts.t2                 = ReggeTransferForAlpha(*param, test.lower_pdg, 0.82);
    gra::gpom::AmpCache off_cache;
    const auto          off = gra::gpom::Resonance(off_lts, *param, ls, &off_cache);
    REQUIRE(off.size() == 1);
    REQUIRE(off_cache.sources.size() == 1);
    const auto &off_source  = off_cache.sources.front();
    const auto &off_block   = GPSpinBlock(off_source, 2 * test.s);
    std::size_t off_nonzero = 0;
    for (int m1 = -off_source.basis.MMAX; m1 <= off_source.basis.MMAX; ++m1) {
      for (int m2 = -off_source.basis.MMAX; m2 <= off_source.basis.MMAX; ++m2) {
        if (!active_mu(m1 - m2)) { continue; }
        const auto &value     = GPSpinValue(off_source, off_block, m1, m2);
        const auto &reflected = GPSpinValue(off_source, off_block, -m1, -m2);
        RequireComplexNear(reflected, static_cast<double>(reflection_phase) * value, 2.0e-12);
        off_nonzero += std::abs(value) > 1.0e-12 ? 1U : 0U;
      }
    }
    CHECK(off_nonzero > 0);

    if (test.identical) {
      gra::LORENTZSCALAR exchanged_lts = off_lts;
      std::swap(exchanged_lts.t1, exchanged_lts.t2);
      gra::gpom::AmpCache exchanged_cache;
      const auto          exchanged = gra::gpom::Resonance(exchanged_lts, *param, ls, &exchanged_cache);
      REQUIRE(exchanged.size() == 1);
      REQUIRE(exchanged_cache.sources.size() == 1);
      const auto &exchanged_source = exchanged_cache.sources.front();
      const auto &exchanged_block  = GPSpinBlock(exchanged_source, 2 * test.s);
      for (int m1 = -off_source.basis.MMAX; m1 <= off_source.basis.MMAX; ++m1) {
        for (int m2 = -off_source.basis.MMAX; m2 <= off_source.basis.MMAX; ++m2) {
          if (!active_mu(m1 - m2)) { continue; }
          const auto &value           = GPSpinValue(off_source, off_block, m1, m2);
          const auto &exchanged_value = GPSpinValue(exchanged_source, exchanged_block, m2, m1);
          RequireComplexNear(exchanged_value, value, 2.0e-12);
        }
      }
    }
  }
}

// Check coherent tensor GP LS fusion at and away from the fixed-spin pole
TEST_CASE("GP coherent complex fusion stays bilinear around the f2 pole",
          "[gra::MRegge][GP][spin][parity][resonance][analytic]"
          "[normalization][regression]") {
  gra::LORENTZSCALAR lts        = MakeToyReggeLTSAsymmetric(0.3, 4.2, -4.6);
  lts.process.MMAX              = 2;
  lts.process.SPINGEN           = true;
  lts.process.DERIVATIVE_FACTOR = false;
  gra::MRegge regge(lts, gra::MModelTune::Load(modelfile),
                    gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));
  const auto  param = ReggeParametersForTest(regge, lts);

  const double pole_transfer = ReggeTransferForAlpha(*param, 990, 2.0);
  lts.t1                     = pole_transfer;
  lts.t2                     = pole_transfer;

  const double exchange_mass   = 0.9;
  const double momentum        = 0.37;
  const double exchange_energy = std::sqrt(gra::math::pow2(exchange_mass) + gra::math::pow2(momentum));
  lts.q1_in_X                  = gra::M4Vec(0.0, 0.0, momentum, exchange_energy);
  lts.q2_in_X                  = gra::M4Vec(0.0, 0.0, -momentum, exchange_energy);

  const std::vector<gra::spin::LSTerm> terms = {{0, 4, std::complex<double>(0.63, -0.27)},
                                                {2, 0, std::complex<double>(-0.31, 0.44)}};

  // Build one coherent GP LS resonance with a common coupling phase
  const auto make_res = [&](const std::complex<double> phase) {
    gra::PARAM_RES res = MakeToyTensorGPResonance(lts.process.MMAX);
    auto          &hel = res.production.front().hel;
    hel.alpha_ls.Clear();
    for (const auto &term : terms) { hel.alpha_ls.Set(term.l, term.two_s, phase * term.coefficient); }
    gra::gpom::InitResonanceLS(hel, res.p.spinX2 / 2, 4, 4);
    return res;
  };

  const gra::PARAM_RES res    = make_res(1.0);
  const auto          &pole   = gra::regge::PoleRepresentative(*param, lts.PDG, 990);
  const auto           vertex = gra::spin::PreparePoleLS(res.p, pole, pole, terms, 1.0, true, true, true,
                                                         gra::spin::VertexContext::Auto, 0.0, false);
  const auto           fixed  = gra::xpom::Fusion(lts, vertex);
  REQUIRE(fixed.size_row() == 25);
  REQUIRE(fixed.size_col() == 5);
  REQUIRE(fixed.FrobNorm2() > 1.0e-20);

  gra::gpom::AmpCache cache;
  const auto          actual = gra::gpom::Resonance(lts, *param, res, &cache);
  REQUIRE(actual.size() == 1);
  REQUIRE(cache.sources.size() == 1);
  const auto &source = cache.sources.front();
  REQUIRE_FALSE(source.basis.upper_photon);
  REQUIRE_FALSE(source.basis.lower_photon);
  // Check whether one source retains a nontrivial azimuthal section
  const auto has_phase = [](const auto &factors) {
    return std::any_of(factors.cbegin(), factors.cend(),
                       [](const auto &value) { return std::abs(std::imag(value)) > 1.0e-6; });
  };
  REQUIRE(has_phase(source.upper_column_factor));
  REQUIRE(has_phase(source.lower_column_factor));
  const auto destination = gra::spin::CanonicalProtonPairSpinLayout::KroneckerDestinationRows(
      source.upper_row_factor.size(), source.lower_row_factor.size());
  const auto expected =
      fixed.FactorizedKroneckerMultiply(source.upper_row_factor, source.upper_column_factor, source.lower_row_factor,
                                        source.lower_column_factor, destination);
  RequireMatrixNear(actual.front(), expected, 2.0e-12);

  gra::LORENTZSCALAR off_lts = lts;
  off_lts.t1                 = ReggeTransferForAlpha(*param, 990, 1.37);
  off_lts.t2                 = ReggeTransferForAlpha(*param, 990, 1.61);
  gra::gpom::AmpCache  off_cache;
  const auto           off        = gra::gpom::Resonance(off_lts, *param, res, &off_cache);
  const gra::PARAM_RES phased     = make_res(gra::math::zi);
  const auto           off_phased = gra::gpom::Resonance(off_lts, *param, phased, &off_cache);
  REQUIRE(off.size() == 1);
  REQUIRE(off_phased.size() == 1);
  REQUIRE(off_cache.sources.size() == 1);
  REQUIRE(off.front().FrobNorm2() > 1.0e-20);
  RequireMatrixNear(off_phased.front(), off.front() * gra::math::zi, 2.0e-12);
}

// Check frame independence after sewing every exchange-helicity row to a physical scalar decay
TEST_CASE("MP pole production and decay agree in CM CS and HX",
          "[gra::MRegge][MP][spin][frame][normalization][regression]") {
  const auto                            asymmetric = AsymmetricF2PhasePointForTest(0.91, -0.38);
  const std::vector<gra::LORENTZSCALAR> events     = {asymmetric,
                                                      RotateToyEventAroundZ(asymmetric, 1.17),
                                                      BoostToyEventAlongZ(asymmetric, -asymmetric.pfinal[0].Rap()),
                                                      BoostToyEventAlongZ(asymmetric, 1.3),
                                                      BoostToyEventAlongZ(asymmetric, -1.1),
                                                      ScalarPolePhasePointForTest(0.20, 0.91, -0.38, 1.30)};
  const auto                            seed       = MakeRealisticTensorMPResonance();
  const auto                            exchange   = seed.production.front().tree.front().p;
  for (const int J : {0, 1, 2, 4, 6}) {
    gra::PARAM_RES res = J == 1 ? MakeToyPhotoMPResonance() : seed;
    res.p.spinX2       = 2 * J;
    res.p.C            = J == 1 ? -1 : 1;
    auto vertex        = res.production.front().pole.value();
    if (J != 1) {
      const auto operators = gra::spin::CanonicalPoleOperators(res.p, exchange, exchange);
      REQUIRE_FALSE(operators.empty());
      std::vector<gra::spin::LSTerm> terms;
      for (const auto i : indices(operators)) {
        terms.push_back(
            {operators[i].coupling.l, operators[i].coupling.two_s, std::polar(0.61 + 0.03 * i, -0.37 + 0.17 * i)});
      }
      vertex = gra::spin::PreparePoleLS(res.p, exchange, exchange, terms, 1.0);
    }
    const auto decay_vertex = gra::spin::PreparePoleLS(res.p, asymmetric.decaytree[0].p, asymmetric.decaytree[1].p,
                                                       {{static_cast<std::size_t>(J), 0, {0.73, -0.19}}}, 1.0, false);
    res.hel_decay           = gra::spin::PoleLSHelicity(decay_vertex, 1.0);
    for (const auto event : indices(events)) {
      auto lts             = events[event];
      lts.process.SPINGEN  = true;
      lts.process.SPINDEC  = true;
      lts.process.MP_FRAME = "CM";
      const auto cm        = gra::mpom::Fusion(lts, vertex);
      const auto reference = cm * gra::spin::ResonanceDecayMatrix(lts, res, "CM");
      REQUIRE(reference.FrobNorm2() > 1.0e-20);
      RequireMatrixNear(cm, gra::xpom::Fusion(lts, vertex), 2.0e-11);
      for (const std::string frame : {"CM", "CS", "HX"}) {
        CAPTURE(J, event, frame);
        lts.process.MP_FRAME = frame;
        const auto central   = gra::mpom::Fusion(lts, vertex);
        const auto amplitude = central * gra::spin::ResonanceDecayMatrix(lts, res, frame);
        CHECK(central.FrobNorm2() == Approx(cm.FrobNorm2()).epsilon(2.0e-11));
        RequireMatrixNear(amplitude, reference, 2.0e-10);
      }
    }
  }
}

// Check complex RES+CON amplitudes under frame steering and common azimuth rotations
TEST_CASE("Regge coherent amplitudes preserve their root spin basis",
          "[gra::MRegge][MP][XP][GP][spin][frame][interference][regression]") {
  for (const auto model :
       {gra::ReggeProductionModel::MP, gra::ReggeProductionModel::XP, gra::ReggeProductionModel::GP}) {
    for (const auto &final : std::array<std::pair<const char *, std::array<int, 2>>, 3>{
             {{"pi+ pi-", {211, -211}}, {"p+ p-", {2212, -2212}}, {"rho(770)0 rho(770)0", {113, 113}}}}) {
      CAPTURE(final.first);
      ToyHelicityProcess process;
      ConfigureToyProductionProcess(process, gra::ReggeProductionModelName(model), "RES+CON", final.first);
      const std::string name      = final.second[0] == 211 ? "f2_1270" : "f2_2150";
      auto              resonance = gra::resonance::Read("RES/" + name + ".json", process.state.random, model);
      resonance.spin_basis        = "none";
      if (model == gra::ReggeProductionModel::MP) {
        for (auto &channel : resonance.MP.channels) {
          channel.basis = gra::ReggeVertexBasis::LS;
          channel.g_ls.Clear();
          const auto operators = gra::spin::CanonicalPoleOperators(
              resonance.p, process.state.lts.PDG.FindByPDG(channel.exchange[0]),
              process.state.lts.PDG.FindByPDG(channel.exchange[1]), true, channel.C_symmetry, channel.P_symmetry);
          for (const auto i : indices(operators)) {
            const auto &term = operators[i].coupling;
            channel.g_ls.Set(term.l, term.two_s, std::polar(0.61 + 0.03 * i, -0.37 + 0.17 * i));
          }
        }
      }
      process.SetResonances({{name, resonance}});
      REQUIRE_NOTHROW(process.InitializeProcessAmplitude());
      auto base = DirectCentralPairLTSForTest(final.second[0], final.second[1]);
      RefreshToyDerivedKinematicsPreserveDecay(base);
      base.process                = process.state.lts.process;
      base.process.FORWARD_NOFLIP = true;
      base.hamp.Configure(process.state.lts.hamp.metadata);
      const auto evaluate = [&](const std::string &frame, double angle) {
        auto lts             = RotateToyEventAroundZ(base, angle);
        lts.process.MP_FRAME = frame;
        gra::MRegge regge(lts, process.state.model_tune,
                          gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::ResonanceContinuumTwoBody, "Regge_frame"));
        REQUIRE(regge.Amp2(lts, model, gra::MReggeMode::ResonanceContinuumTwoBody) > 0.0);
        return lts.hamp;
      };
      const auto reference = evaluate("CM", 0.0);
      for (const double angle : {0.0, 0.73, -1.21}) {
        for (const std::string frame : {"CM", "CS", "HX"}) {
          CAPTURE(model, angle, frame);
          // XP and GP always sew in CM, independently of MP frame steering
          RequireVectorNear(evaluate(frame, angle), reference, 2.0e-10);
        }
      }
    }
  }
}

// Check that a fixed M = 0 selection denotes different physical states in different axes
TEST_CASE("MP fixed polarization retains physical frame dependence",
          "[gra::MRegge][MP][spin][frame][polarization][regression]") {
  ToyHelicityProcess process;
  ConfigureToyProductionProcess(process, "MP", "RES", "pi+ pi-");
  auto resonance       = gra::resonance::Read("RES/f2_1270.json", process.state.random, gra::ReggeProductionModel::MP);
  resonance.spin_basis = "a_Jz";
  resonance.a_Jz       = {0.0, 0.0, 1.0, 0.0, 0.0};
  process.SetResonances({{"f2_1270", resonance}});
  process.InitializeProcessAmplitude();
  auto base    = AsymmetricF2PhasePointForTest(0.91, -0.38);
  base.process = process.state.lts.process;
  base.hamp.Configure(process.state.lts.hamp.metadata);
  const auto evaluate = [&](const std::string &frame) {
    auto lts             = base;
    lts.process.MP_FRAME = frame;
    gra::MRegge  regge(lts, process.state.model_tune,
                       gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Resonance, "MP_polarization"));
    const double norm = regge.Amp2(lts, gra::ReggeProductionModel::MP, gra::MReggeMode::Resonance);
    REQUIRE(norm > 0.0);
    return norm;
  };
  const double cm = evaluate("CM");
  for (const std::string frame : {"CS", "HX"}) {
    const double selected = evaluate(frame);
    CAPTURE(frame, cm, selected);
    CHECK(std::abs(selected / cm - 1.0) > 1.0e-7);
  }
}

// Check refreshed exchange momenta, retained decay matrices and restored Born amplitudes
TEST_CASE("Regge exchange rest momenta and amplitude caches follow screening shifts",
          "[gra::MRegge][MP][XP][GP][frame][cache][screening][regression]") {
  for (const auto model :
       {gra::ReggeProductionModel::MP, gra::ReggeProductionModel::XP, gra::ReggeProductionModel::GP}) {
    ToyHelicityProcess process;
    ConfigureToyProductionProcess(process, gra::ReggeProductionModelName(model), "RES+CON", "pi+ pi-");
    auto resonance       = gra::resonance::Read("RES/f2_1270.json", process.state.random, model);
    process.SetResonances({{"f2_1270", resonance}});
    process.InitializeProcessAmplitude();
    auto base    = AsymmetricF2PhasePointForTest(0.91, -0.38);
    base.process = process.state.lts.process;
    base.hamp.Configure(process.state.lts.hamp.metadata);
    process.state.lts = base;
    PrepareScreeningPoint(process.state);
    base = process.state.lts;
    for (const std::string frame : {"CM", "CS", "HX"}) {
      auto &lts            = process.state.lts;
      lts                  = base;
      lts.process.MP_FRAME = frame;
      lts.pfinal_orig      = lts.pfinal;
      lts.amplitude.BeginCentral();
      gra::MRegge regge(lts, process.state.model_tune,
                        gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::ResonanceContinuumTwoBody, "Regge_cache"));
      REQUIRE(regge.Amp2(lts, model, gra::MReggeMode::ResonanceContinuumTwoBody) > 0.0);
      const auto born      = lts.hamp;
      lts.screening.active = true;
      for (const std::array<double, 2> shift : {std::array<double, 2>{0.07, -0.03}, {-0.04, 0.11}, {0.07, -0.03}}) {
        CAPTURE(model, frame, shift[0], shift[1]);
        REQUIRE(process.LoopKinematics(
            {base.pfinal[1].Px() - shift[0], base.pfinal[1].Py() - shift[1]},
            {base.pfinal[2].Px() + shift[0], base.pfinal[2].Py() + shift[1]}));
        const auto q1 = gra::kinematics::BoostToRestFrame(lts.q1, lts.pfinal[0], "test upper exchange");
        const auto q2 = gra::kinematics::BoostToRestFrame(lts.q2, lts.pfinal[0], "test lower exchange");
        for (std::size_t component = 0; component < 4; ++component) {
          CHECK(lts.q1_in_X[component] == Approx(q1[component]).margin(2.0e-11));
          CHECK(lts.q2_in_X[component] == Approx(q2[component]).margin(2.0e-11));
        }
        CHECK((lts.q1_in_X + lts.q2_in_X).P3mod() < 2.0e-10);
        CHECK((lts.q1_in_X - base.q1_in_X).P3mod() > 1.0e-3);
        REQUIRE_NOTHROW(regge.Amp2(lts, model, gra::MReggeMode::ResonanceContinuumTwoBody));
        REQUIRE(lts.proton_good_walker.has_value());
        auto fresh = lts;
        fresh.amplitude.AbortCentral();
        fresh.screening.active = false;
        REQUIRE(regge.Amp2(fresh, model, gra::MReggeMode::ResonanceContinuumTwoBody) > 0.0);
        const auto& cached = lts.proton_good_walker->components;
        const auto& rebuilt = fresh.proton_good_walker->components;
        REQUIRE(cached.size() == rebuilt.size());
        for (const auto component : indices(cached)) {
          RequireMatrixNear(cached[component].source, rebuilt[component].source, 2.0e-10);
        }
        RequireVectorNear(gra::ProjectGoodWalker(*lts.proton_good_walker), fresh.hamp, 2.0e-10);
      }
      REQUIRE(process.RestoreBornKinematics());
      for (std::size_t component = 0; component < 4; ++component) {
        CHECK(lts.q1_in_X[component] == Approx(base.q1_in_X[component]).margin(2.0e-11));
        CHECK(lts.q2_in_X[component] == Approx(base.q2_in_X[component]).margin(2.0e-11));
      }
      REQUIRE_NOTHROW(regge.Amp2(lts, model, gra::MReggeMode::ResonanceContinuumTwoBody));
      REQUIRE(lts.proton_good_walker.has_value());
      RequireVectorNear(gra::ProjectGoodWalker(*lts.proton_good_walker), born, 2.0e-10);
      lts.screening.active = false;
      lts.amplitude.AbortCentral();
    }
  }
}

// Check every supported integer spin and both parities in all Regge models
TEST_CASE("MP XP and GP pole bases agree for spins zero through six",
          "[gra::MRegge][MP][XP][GP][spin][parity][resonance][normalization]"
          "[unit_coupling][regression]") {
  struct PoleCase {
    int         J;
    int         parity;
    std::size_t l;
    std::size_t s;
  };

  gra::LORENTZSCALAR lts        = MakeToyReggeLTSAsymmetric(0.3, 4.2, -4.6);
  lts.process.MMAX              = 2;
  lts.process.SPINGEN           = true;
  lts.process.DERIVATIVE_FACTOR = false;
  gra::MRegge  regge(lts, gra::MModelTune::Load(modelfile),
                     gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));
  const auto   param            = ReggeParametersForTest(regge, lts);
  const double pomeron_transfer = ReggeTransferForAlpha(*param, 990, 2.0);
  const double reggeon_transfer = ReggeTransferForAlpha(*param, 9910, 2.0);
  lts.t1                        = pomeron_transfer;
  lts.t2                        = pomeron_transfer;

  const double exchange_mass   = 0.9;
  const double momentum        = 0.37;
  const double exchange_energy = std::sqrt(gra::math::pow2(exchange_mass) + gra::math::pow2(momentum));
  lts.q1_in_X                  = gra::M4Vec(0.0, 0.0, momentum, exchange_energy);
  lts.q2_in_X                  = gra::M4Vec(0.0, 0.0, -momentum, exchange_energy);
  const auto &pomeron_pole     = gra::regge::PoleRepresentative(*param, lts.PDG, 990);
  const auto &reggeon_pole     = gra::regge::PoleRepresentative(*param, lts.PDG, 9910);
  const auto  pomeron          = lts.PDG.FindByPDG(990);
  const auto  reggeon          = lts.PDG.FindByPDG(9910);

  std::vector<PoleCase> cases;
  for (int J = 0; J <= 6; ++J) {
    for (const int parity : {-1, 1}) {
      auto mother   = MakeToyTensorGPResonance(lts.process.MMAX).p;
      mother.spinX2 = 2 * J;
      mother.P      = parity;
      const auto operators = gra::spin::CanonicalPoleOperators(mother, pomeron_pole, pomeron_pole, true, true, true,
                                                               gra::spin::VertexContext::Auto);
      REQUIRE_FALSE(operators.empty());
      for (const auto &op : operators) { cases.push_back({J, parity, op.coupling.l, op.coupling.two_s / 2}); }
    }
  }

  constexpr double theta = 0.83;
  constexpr double phi   = -0.47;
  const gra::M4Vec nonaxial_q1(momentum * std::sin(theta) * std::cos(phi), momentum * std::sin(theta) * std::sin(phi),
                               momentum * std::cos(theta), exchange_energy);
  const gra::M4Vec nonaxial_q2(-nonaxial_q1.Px(), -nonaxial_q1.Py(), -nonaxial_q1.Pz(), exchange_energy);

  for (const auto &test : cases) {
    CAPTURE(test.J, test.parity, test.l, test.s);
    gra::PARAM_RES ls = MakeToyTensorGPResonance(lts.process.MMAX);
    ls.p.spinX2       = 2 * test.J;
    ls.p.P            = test.parity;
    auto &ls_hel      = ls.production.front().hel;
    ls_hel.alpha_ls.Clear();
    ls_hel.alpha_ls.Set(test.l, 2 * test.s, 1.0);
    gra::gpom::InitResonanceLS(ls_hel, test.J, 4, 4);

    const auto vertex = gra::spin::PreparePoleLS(ls.p, pomeron_pole, pomeron_pole, {{test.l, 2 * test.s, 1.0}}, 1.0,
                                                 true, true, true, gra::spin::VertexContext::Auto, 0.0, false);
    const auto fixed  = gra::xpom::Fusion(lts, vertex);
    REQUIRE(fixed.size_row() == 25);
    REQUIRE(fixed.size_col() == static_cast<std::size_t>(2 * test.J + 1));
    REQUIRE(fixed.FrobNorm2() > 1.0e-20);

    gra::LORENTZSCALAR fixed_lts                                                        = lts;
    fixed_lts.process.MP_FRAME                                                          = "CM";
    fixed_lts.q1                                                                        = lts.q1_in_X;
    fixed_lts.q2                                                                        = lts.q2_in_X;
    fixed_lts.pfinal[0]                                                                 = fixed_lts.q1 + fixed_lts.q2;
    const std::array<std::pair<const char *, gra::ReggeProductionModel>, 2> pole_models = {{
        {"MP", gra::ReggeProductionModel::MP},
        {"XP", gra::ReggeProductionModel::XP},
    }};
    for (const auto &[label, model] : pole_models) {
      CAPTURE(label);
      const auto actual =
          (model == gra::ReggeProductionModel::MP ? gra::mpom::Fusion : gra::xpom::Fusion)(fixed_lts, vertex);
      RequireMatrixNear(actual, fixed, 2.0e-12);
    }

    gra::LORENTZSCALAR nonaxial_lts = fixed_lts;
    nonaxial_lts.q1                 = nonaxial_q1;
    nonaxial_lts.q2                 = nonaxial_q2;
    nonaxial_lts.q1_in_X            = nonaxial_q1;
    nonaxial_lts.q2_in_X            = nonaxial_q2;
    nonaxial_lts.pfinal[0]          = nonaxial_q1 + nonaxial_q2;
    const auto nonaxial_expected    = gra::xpom::Fusion(nonaxial_lts, vertex);
    const auto nonaxial_actual      = gra::mpom::Fusion(nonaxial_lts, vertex);
    RequireMatrixNear(nonaxial_actual, nonaxial_expected, 2.0e-12);

    gra::gpom::AmpCache ls_cache;
    const auto          ls_actual = gra::gpom::Resonance(lts, *param, ls, &ls_cache);
    REQUIRE(ls_actual.size() == 1);
    REQUIRE(ls_cache.sources.size() == 1);
    const auto &ls_source      = ls_cache.sources.front();
    const auto  ls_destination = gra::spin::CanonicalProtonPairSpinLayout::KroneckerDestinationRows(
         ls_source.upper_row_factor.size(), ls_source.lower_row_factor.size());
    const auto ls_expected =
        fixed.FactorizedKroneckerMultiply(ls_source.upper_row_factor, ls_source.upper_column_factor,
                                          ls_source.lower_row_factor, ls_source.lower_column_factor, ls_destination);
    RequireMatrixNear(ls_actual.front(), ls_expected, 2.0e-11);

    // Disabling production spin keeps the same pole coupling in every model
    auto blind_lts = fixed_lts;
    blind_lts.process.SPINGEN = false;
    const auto gp_blind = gra::gpom::Resonance(blind_lts, *param, ls, nullptr);
    auto fixed_res = ls;
    fixed_res.production.front().pole = vertex;
    for (const auto &[label, model] : pole_models) {
      CAPTURE(label);
      const auto fixed_blind = gra::rspin::Resonance(
          blind_lts, fixed_res,
          (model == gra::ReggeProductionModel::MP ? gra::mpom::Fusion : gra::xpom::Fusion), param->s0);
      RequireMatrixNear(gp_blind.front(), fixed_blind.front(), 2.0e-11);
    }

    const auto mixed_vertex =
        gra::spin::PreparePoleLS(ls.p, pomeron_pole, reggeon_pole, {{test.l, 2 * test.s, 1.0}}, 1.0, true, true, true,
                                 gra::spin::VertexContext::Auto, 0.0, false);
    const auto mixed_reduced = gra::spin::PoleLSReduced(mixed_vertex, momentum);
    REQUIRE(mixed_reduced.size_row() == 5);
    REQUIRE(mixed_reduced.size_col() == 5);

    gra::RES_PRODUCTION_CHANNEL direct_model;
    direct_model.exchange     = {990, 9910};
    direct_model.basis        = gra::ReggeVertexBasis::Helicity;
    direct_model.C_symmetry   = true;
    direct_model.P_symmetry   = true;
    std::size_t parity_orbits = 0;
    for (int m1 = -lts.process.MMAX; m1 <= lts.process.MMAX; ++m1) {
      for (int m2 = -lts.process.MMAX; m2 <= lts.process.MMAX; ++m2) {
        if (std::abs(m1 - m2) > test.J) { continue; }
        const auto projection = std::make_pair(m1, m2);
        const auto reflected  = std::make_pair(-m1, -m2);
        if (projection > reflected) { continue; }
        const std::size_t          i1 = gra::gpom::AnalyticMIndex(m1, lts.process.MMAX, "GP mixed direct upper spin");
        const std::size_t          i2 = gra::gpom::AnalyticMIndex(m2, lts.process.MMAX, "GP mixed direct lower spin");
        const std::complex<double> coupling = mixed_reduced[i1][i2];
        if (std::abs(coupling) <= 1.0e-13) { continue; }
        direct_model.helicity.push_back({static_cast<double>(m1), static_cast<double>(m2)});
        direct_model.g_helicity.push_back(coupling);
        parity_orbits += projection != reflected ? 1U : 0U;
      }
    }
    REQUIRE_FALSE(direct_model.helicity.empty());
    REQUIRE(parity_orbits > 0);
    const auto direct_hel =
        gra::gpom::PrepareResonance(ls.p, {pomeron, reggeon}, direct_model, *param, lts.PDG, lts.process.MMAX, 0.0);
    RequireMatrixNear(direct_hel.T, mixed_reduced, 2.0e-12);

    gra::PARAM_RES direct               = ls;
    direct.production.front().tree[0].p = pomeron;
    direct.production.front().tree[1].p = reggeon;
    direct.GP.channels                  = {direct_model};
    direct.production.front().hel       = direct_hel;
    gra::LORENTZSCALAR direct_lts       = lts;
    direct_lts.t1                       = pomeron_transfer;
    direct_lts.t2                       = reggeon_transfer;
    const auto          mixed_fixed     = gra::xpom::Fusion(direct_lts, mixed_vertex);
    gra::gpom::AmpCache direct_cache;
    const auto          direct_actual = gra::gpom::Resonance(direct_lts, *param, direct, &direct_cache);
    REQUIRE(direct_actual.size() == 1);
    REQUIRE(direct_cache.sources.size() == 1);
    const auto &direct_source = direct_cache.sources.front();
    CHECK(direct_source.basis.alpha1 == Approx(2.0).epsilon(1.0e-12));
    CHECK(direct_source.basis.alpha2 == Approx(2.0).epsilon(1.0e-12));
    const auto direct_destination = gra::spin::CanonicalProtonPairSpinLayout::KroneckerDestinationRows(
        direct_source.upper_row_factor.size(), direct_source.lower_row_factor.size());
    const auto direct_expected = mixed_fixed.FactorizedKroneckerMultiply(
        direct_source.upper_row_factor, direct_source.upper_column_factor, direct_source.lower_row_factor,
        direct_source.lower_column_factor, direct_destination);
    RequireMatrixNear(direct_actual.front(), direct_expected, 2.0e-11);
    direct_lts.process.SPINGEN = false;
    const auto direct_blind = gra::gpom::Resonance(direct_lts, *param, direct, nullptr);
    RequireMatrixNear(direct_blind.front(), gp_blind.front(), 2.0e-11);
  }
}

// Check physical photon projection and LS interference in spin-disabled GP production
TEST_CASE("GP spin-disabled photon poles keep LS and helicity normalization",
          "[gra::MRegge][GP][spin][photon][normalization][regression]") {
  auto lts = MakeToyReggeLTSAsymmetric(0.3, 4.2, -4.6);
  lts.process.MMAX = 2;
  lts.process.SPINGEN = false;
  lts.process.DERIVATIVE_FACTOR = false;
  gra::MRegge regge(lts, gra::MModelTune::Load(modelfile),
                    gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "photon_pole_scale"));
  const auto param = ReggeParametersForTest(regge, lts);
  const auto photon = lts.PDG.FindByPDG(gra::PDG::PDG_gamma);
  const auto pomeron = lts.PDG.FindByPDG(990);
  const auto pole = gra::regge::PoleRepresentative(*param, lts.PDG, 990);
  gra::PARAM_RES res = MakeToyTensorGPResonance(lts.process.MMAX);
  res.p = lts.PDG.FindByPDG(113);
  res.production.front().tree[0].p = photon;
  res.production.front().tree[1].p = pomeron;

  for (const auto &terms : std::vector<std::vector<gra::spin::LSTerm>>{
           {{0, 2, 1.0}}, {{2, 2, 1.0}}, {{2, 4, 1.0}}, {{2, 6, 1.0}}, {{4, 6, 1.0}},
           {{0, 2, {0.7, -0.2}}, {2, 2, {-0.3, 0.5}}, {2, 4, {0.4, 0.1}}}}) {
    const auto                  vertex = gra::spin::PreparePoleLS(res.p, photon, pole, terms, 1.0, true, true, true,
                                                                  gra::spin::VertexContext::Auto, 0.0, false);
    gra::RES_PRODUCTION_CHANNEL ls;
    ls.exchange = {22, 990};
    ls.basis = gra::ReggeVertexBasis::LS;
    ls.C_symmetry = true;
    ls.P_symmetry = true;
    for (const auto &op : gra::spin::CanonicalPoleOperators(res.p, photon, pole)) {
      ls.g_ls.Set(op.coupling.l, op.coupling.two_s, 0.0);
    }
    for (const auto &term : terms) { ls.g_ls.Set(term.l, term.two_s, term.coefficient); }
    res.production.front().hel = gra::gpom::PrepareResonance(res.p, {photon, pomeron}, ls, *param, lts.PDG, 2, 0.0);
    const auto ls_amplitude = gra::gpom::Resonance(lts, *param, res, nullptr).front();

    auto direct = ls;
    direct.basis = gra::ReggeVertexBasis::Helicity;
    direct.g_ls.Clear();
    const auto reduced = gra::spin::PoleLSReduced(vertex, 1.0);
    for (std::size_t i = 0; i < reduced.size_col(); ++i) {
      if (std::abs(reduced[0][i]) < 1.0e-13) { continue; }
      direct.helicity.push_back({-1.0, static_cast<double>(i) - 2.0});
      direct.g_helicity.push_back(reduced[0][i]);
    }
    res.production.front().hel = gra::gpom::PrepareResonance(res.p, {photon, pomeron}, direct, *param, lts.PDG, 2, 0.0);
    const auto helicity_amplitude = gra::gpom::Resonance(lts, *param, res, nullptr).front();
    RequireMatrixNear(ls_amplitude, helicity_amplitude, 2.0e-11);
    const std::size_t rows = gra::spin::SpinHalfTransitions(lts.process.FORWARD_NOFLIP).size();
    const auto expected    = gra::spin::Blind(rows * rows, 3, std::sqrt(gra::spin::LeadingPoleDensity(vertex, 1.0)));
    RequireMatrixNear(ls_amplitude, expected, 2.0e-11);
  }
}

// Check the Jacob Wick phase when one mixed GP helicity card is beam reversed
TEST_CASE("GP fusion helicity transports the mixed leg exchange phase",
          "[gra::MRegge][GP][spin][beam-exchange][phase][regression]") {
  gra::LORENTZSCALAR lts = MakeToyReggeLTSAsymmetric(0.3, 4.2, -4.6);
  lts.process.MMAX       = 2;
  gra::MRegge regge(lts, gra::MModelTune::Load(modelfile),
                    gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));
  const auto  param = ReggeParametersForTest(regge, lts);

  gra::MParticle mother     = lts.PDG.FindByPDG(113);
  mother.P                  = 1;
  mother.C                  = 1;
  const auto   pomeron      = lts.PDG.FindByPDG(990);
  const auto   reggeon      = lts.PDG.FindByPDG(9910);
  const auto  &pomeron_pole = gra::regge::PoleRepresentative(*param, lts.PDG, pomeron.pdg);
  const auto  &reggeon_pole = gra::regge::PoleRepresentative(*param, lts.PDG, reggeon.pdg);
  const double momentum     = 0.83;
  const auto   direct_vertex   = gra::spin::PreparePoleLS(mother, pomeron_pole, reggeon_pole, {{0, 2, 1.0}}, 1.0, true,
                                                          false, true, gra::spin::VertexContext::Auto, 0.0, false);
  const auto   reversed_vertex = gra::spin::PreparePoleLS(mother, reggeon_pole, pomeron_pole, {{0, 2, -1.0}}, 1.0, true,
                                                          false, true, gra::spin::VertexContext::Auto, 0.0, false);
  const auto   direct          = gra::spin::PoleLSReduced(direct_vertex, momentum);
  const auto   expected        = gra::spin::PoleLSReduced(reversed_vertex, momentum);

  gra::RES_PRODUCTION_CHANNEL ls_model;
  ls_model.exchange   = {pomeron.pdg, reggeon.pdg};
  ls_model.basis      = gra::ReggeVertexBasis::LS;
  ls_model.C_symmetry = false;
  ls_model.P_symmetry = true;
  ls_model.Lambda     = 1.0;
  for (const auto &pole : gra::spin::CanonicalPoleOperators(mother, pomeron_pole, reggeon_pole, true, false, true)) {
    const bool selected = pole.coupling.l == 0 && pole.coupling.two_s == 2;
    ls_model.g_ls.Set(pole.coupling.l, pole.coupling.two_s, selected ? 1.0 : 0.0);
  }
  const auto reversed_ls =
      gra::gpom::PrepareResonance(mother, {reggeon, pomeron}, ls_model, *param, lts.PDG, lts.process.MMAX, 0.0);
  REQUIRE(reversed_ls.alpha_ls.Size() == 1);
  RequireComplexNear(reversed_ls.alpha_ls.At(0, 2), -1.0, 2.0e-12);

  gra::RES_PRODUCTION_CHANNEL model;
  model.exchange   = {pomeron.pdg, reggeon.pdg};
  model.basis      = gra::ReggeVertexBasis::Helicity;
  model.C_symmetry = false;
  model.P_symmetry = true;
  for (int m1 = -lts.process.MMAX; m1 <= lts.process.MMAX; ++m1) {
    for (int m2 = -lts.process.MMAX; m2 <= lts.process.MMAX; ++m2) {
      const auto projection = std::make_pair(m1, m2);
      const auto reflected  = std::make_pair(-m1, -m2);
      if (projection > reflected || std::abs(m1 - m2) > mother.spinX2 / 2) { continue; }
      const std::size_t i1 = gra::gpom::AnalyticMIndex(m1, lts.process.MMAX, "GP direct exchange phase m1");
      const std::size_t i2 = gra::gpom::AnalyticMIndex(m2, lts.process.MMAX, "GP direct exchange phase m2");
      if (std::abs(direct[i1][i2]) <= 1.0e-13) { continue; }
      model.helicity.push_back({static_cast<double>(m1), static_cast<double>(m2)});
      model.g_helicity.push_back(direct[i1][i2]);
    }
  }
  REQUIRE_FALSE(model.helicity.empty());

  const auto reversed = gra::gpom::PrepareResonance(mother, {reggeon, pomeron}, model, *param, lts.PDG, lts.process.MMAX, 0.0);
  RequireMatrixNear(reversed.T, expected, 2.0e-12);
}

// Check negative orbital and spin reflection phases in a crossed LS vertex
TEST_CASE("GP crossed LS preserves reflected helicities off pole",
          "[gra::MRegge][GP][spin][parity][continuum][regression]") {
  gra::HELMatrix crossed;
  gra::gpom::InitCrossed(crossed, 1.0, 1.0, 0, "crossed GP parity test");
  crossed.coupling_basis = gra::CouplingBasis::LS;
  crossed.P_symmetry     = true;
  crossed.m_ls.resize(1);
  crossed.m_ls[0].Set(1, 2, 1.0);
  gra::gpom::InitCrossedLS(crossed, 1);

  const auto                 evaluated = gra::gpom::Crossed(crossed, 1.02, gra::M4Vec{}, false);
  std::optional<std::size_t> positive;
  std::optional<std::size_t> reflected;
  for (std::size_t row = 0; row < evaluated.helicity.lambda_values.size_row(); ++row) {
    const double lambda1 = evaluated.helicity.lambda_values[row][0];
    const double lambda2 = evaluated.helicity.lambda_values[row][1];
    if (std::abs(lambda1 - 1.0) < 1.0e-12 && std::abs(lambda2) < 1.0e-12) { positive = row; }
    if (std::abs(lambda1 + 1.0) < 1.0e-12 && std::abs(lambda2) < 1.0e-12) { reflected = row; }
  }
  REQUIRE(positive.has_value());
  REQUIRE(reflected.has_value());
  const auto value        = evaluated.helicity.T[*positive][0];
  const auto parity_value = evaluated.helicity.T[*reflected][0];
  REQUIRE(std::abs(value) > 1.0e-10);
  RequireComplexNear(parity_value, value, 2.0e-12);

  const auto pole       = gra::gpom::Crossed(crossed, 1.0, gra::M4Vec{}, false);
  const auto pole_value = pole.helicity.T[*positive][0];
  CHECK(std::abs(pole_value - value) > 1.0e-8);
}

// Check unit transverse photoproduction vertices at a matched spin 2 pole
TEST_CASE(
    "MP XP GP and TP unit rho resonance couplings share matched pole "
    "normalizations",
    "[gra::MRegge][MTensorPomeron][spin][resonance][normalization]"
    "[physics][unit_coupling]") {
  ModelParamRestoreGuard restore_model;
  gra::MODELPARAM = "TUNE0";

  const nlohmann::json rho_ls = {
      {0, 1, 1.0, 0.0}, {2, 1, 0.0, 0.0}, {2, 2, 0.0, 0.0}, {2, 3, 0.0, 0.0}, {4, 3, 0.0, 0.0}};
  const nlohmann::json                              helicity = {{-1, 0, 1.0, 0.0}};
  const std::array<std::string, 2>                  families = {"MP", "XP"};
  const std::array<std::string, 2>                  bases    = {"g_ls", "helicity"};
  std::array<gra::MMatrix<std::complex<double>>, 4> fixed;

  std::size_t cache = 0;
  for (const auto &family : families) {
    for (const auto &basis : bases) {
      CAPTURE(family, basis);
      ToyHelicityProcess process;
      ConfigureToyProductionProcess(process, family, "RES", "pi+ pi-");
      auto  rho              = gra::resonance::Read("RES/rho_770.json", process.state.random, gra::ParseReggeProductionModel(family));
      if (family == "MP") { SetMPFusion(rho); }
      auto &input_channel    = family == "MP" ? rho.MP.channels.front() : rho.XP.channels.front();
      input_channel.exchange = {22, 995};
      if (basis == "g_ls") {
        SetChannelLS(input_channel, rho_ls);
      } else {
        SetChannelHelicity(input_channel, helicity);
      }
      process.SetResonances({{"rho_770", rho}});
      REQUIRE_NOTHROW(process.InitializeProcessAmplitude());
      const auto &configured = process.GetResonances().at("rho_770");
      const auto  topology =
          std::find_if(configured.production.cbegin(), configured.production.cend(), [](const auto &production) {
            const auto &tree = production.tree;
            return tree.size() == 2 && tree[0].p.pdg == 22 && tree[1].p.pdg == 995;
          });
      REQUIRE(topology != configured.production.cend());
      const std::size_t channel_index =
          static_cast<std::size_t>(std::distance(configured.production.cbegin(), topology));
      REQUIRE(channel_index < configured.production.size());
      const auto &vertex = configured.production[channel_index].pole.value();
      fixed[cache++]     = gra::spin::PoleLSReduced(vertex, vertex.Lambda);
    }
  }
  REQUIRE(cache == fixed.size());
  RequireMatrixNear(fixed[2], fixed[0], 2.0e-12);
  RequireMatrixNear(fixed[3], fixed[1], 2.0e-12);

  for (const auto &basis : bases) {
    CAPTURE(basis);
    ToyHelicityProcess gp_process;
    ConfigureToyProductionProcess(gp_process, "GP", "RES", "pi+ pi-");
    auto  gp_rho      = gra::resonance::Read("RES/rho_770.json", gp_process.state.random, gra::ReggeProductionModel::GP);
    auto &gp_input    = gp_rho.GP.channels.front();
    gp_input.exchange = {22, 990};
    if (basis == "g_ls") {
      SetChannelLS(gp_input, rho_ls);
      gp_input.Lambda = 1.0;
    } else {
      SetChannelHelicity(gp_input, helicity);
    }
    gp_process.SetResonances({{"rho_770", gp_rho}});
    REQUIRE_NOTHROW(gp_process.InitializeProcessAmplitude());
    const auto &gp_configured = gp_process.GetResonances().at("rho_770");
    const auto  gp_channel =
        std::find_if(gp_configured.production.cbegin(), gp_configured.production.cend(), [](const auto &production) {
          const auto &tree = production.tree;
          return tree.size() == 2 && tree[0].p.pdg == 22 && tree[1].p.pdg == 990;
        });
    REQUIRE(gp_channel != gp_configured.production.cend());
    const std::size_t gp_index = static_cast<std::size_t>(std::distance(gp_configured.production.cbegin(), gp_channel));
    REQUIRE(gp_index < gp_configured.production.size());
    auto gp = gp_configured.production[gp_index].hel;

    if (basis == "g_ls") {
      REQUIRE(gp.UsesLSCouplings());
      RequireComplexNear(gp.alpha_ls.At(0, 2), 1.0, 1.0e-14);
      REQUIRE(gp.gp_orbital.terms.size() == gp.alpha_ls.Size());
      REQUIRE(gp.gp_orbital.su2.size() == gp.gp_orbital.terms.size() * gp.gp_orbital.nmu);
      continue;
    }

    REQUIRE(gp.UsesHelicityCouplings());
    const int mmax = gp.analytic_MMAX;
    REQUIRE(mmax >= 2);
    for (const int photon_helicity : {-1, 1}) {
      for (int pomeron_helicity = -2; pomeron_helicity <= 2; ++pomeron_helicity) {
        const std::size_t fixed_photon  = static_cast<std::size_t>(photon_helicity + 1);
        const std::size_t fixed_pomeron = static_cast<std::size_t>(pomeron_helicity + 2);
        const std::size_t gp_photon     = gra::gpom::AnalyticMIndex(photon_helicity, mmax, "rho unit pole benchmark");
        const std::size_t gp_pomeron    = gra::gpom::AnalyticMIndex(pomeron_helicity, mmax, "rho unit pole benchmark");
        RequireComplexNear(gp.T[gp_photon][gp_pomeron], fixed[1][fixed_photon][fixed_pomeron], 2.0e-12);
      }
    }
  }

  gra::LORENTZSCALAR tensor_lts;
  tensor_lts.PDG = LoadedPDGTable();
  gra::MTensorPomeron tensor(tensor_lts, gra::MModelTune::Load(modelfile),
                             gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
  const double        rho_mass        = 1.7;
  const double        exchange_mass   = 0.9;
  const double        photon_momentum = (gra::math::pow2(rho_mass) - gra::math::pow2(exchange_mass)) / (2.0 * rho_mass);
  const double        exchange_energy = std::sqrt(gra::math::pow2(exchange_mass) + gra::math::pow2(photon_momentum));
  const gra::M4Vec    photon(0.0, 0.0, photon_momentum, photon_momentum);
  const gra::M4Vec    tensor_pomeron(0.0, 0.0, -photon_momentum, exchange_energy);
  const gra::M4Vec    rho = photon + tensor_pomeron;
  const gra::regge::FFParam no_form_factor;
  const double tp_scale = std::sqrt(6.0) * rho_mass * photon_momentum / (3.0 * gra::math::pow2(exchange_mass));
  const std::array<std::complex<double>, 2> tp_unit = {
      -gra::math::zi * 4.0 * rho_mass * gra::math::pow2(photon_momentum) * (photon_momentum + exchange_energy) *
          tp_scale,
      -gra::math::zi *
          (2.0 * gra::math::pow2(photon_momentum) + 2.0 * photon_momentum * exchange_energy +
           gra::math::pow2(exchange_mass)) *
          tp_scale};
  for (const std::array<double, 2> coupling : {std::array<double, 2>{1.0, 0.0}, std::array<double, 2>{0.0, 1.0}}) {
    const auto                          vertex = tensor.iG_Pvv(rho, photon, coupling[0], coupling[1], no_form_factor);
    std::array<std::complex<double>, 2> tp{};
    for (const auto &i : indices(tp)) {
      const int  h           = i == 0 ? -1 : 1;
      const auto rho_eps     = tensor.EpsMassiveSpin1(rho, h);
      const auto photon_eps  = tensor.EpsSpin1(photon, h);
      const auto pomeron_eps = tensor.EpsMassiveSpin2(tensor_pomeron, 0);
      for (const auto &mu : tensor.LI) {
        for (const auto &nu : tensor.LI) {
          for (const auto &kappa : tensor.LI) {
            for (const auto &lambda : tensor.LI) {
              tp[i] +=
                  std::conj(rho_eps(mu)) * photon_eps(nu) * pomeron_eps(kappa, lambda) * vertex(mu, nu, kappa, lambda);
            }
          }
        }
      }
    }
    CAPTURE(coupling, tp);
    REQUIRE(std::abs(tp[0]) > 1.0e-12);
    RequireComplexNear(tp[0], coupling[0] * tp_unit[0] + coupling[1] * tp_unit[1], 2.0e-12);
    RequireComplexNear(tp[1], tp[0], 2.0e-12);
    for (const auto &i : indices(tp)) {
      RequireComplexNear(tp[i] / tp[0], fixed[1][2 * i][2] / fixed[1][0][2], 2.0e-12);
    }
  }
}

// Check arbitrary configured XP resonance spins through the public amplitude
TEST_CASE("TUNE0 XP resonance operators support arbitrary even spin", "[gra::MRegge][spin][XP][physics][regression]") {
  struct XPSpinCase {
    const char *label;
    const char *card;
    int         spinX2;
  };
  constexpr std::array<XPSpinCase, 3> cases = {{
      {"f2_1270", "RES/f2_1270.json", 4},
      {"f4_2300", "RES/f4_2300.json", 8},
      {"f6_2510", "RES/f6_2510.json", 12},
  }};

  for (const auto &test : cases) {
    CAPTURE(test.label, test.spinX2);
    ToyHelicityProcess process;
    ConfigureToyProductionProcess(process, "XP", "RES", "pi+ pi-");
    const auto resonance = gra::resonance::Read(test.card, process.state.random, gra::ReggeProductionModel::XP);
    process.SetResonances({{test.label, resonance}});
    REQUIRE_NOTHROW(process.InitializeProcessAmplitude());

    auto configured = process.GetResonances().at(test.label);
    REQUIRE(configured.p.spinX2 == test.spinX2);
    for (const auto &production : configured.production) {
      REQUIRE(production.pole.has_value());
      CHECK(production.pole->helicity.Jz_values.size() == static_cast<std::size_t>(test.spinX2 + 1));
    }

    gra::LORENTZSCALAR lts     = ScalarPolePhasePointForTest(0.20, 0.91, -0.38, configured.p.mass);
    lts.process.MP_FRAME       = process.state.lts.process.MP_FRAME;
    lts.process.FORWARD_VERTEX = process.state.lts.process.FORWARD_VERTEX;
    lts.process.FORWARD_NOFLIP = process.state.lts.process.FORWARD_NOFLIP;
    gra::MRegge  regge(lts, process.state.model_tune,
                       gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Resonance, "XP_RES_arbitrary_spin"));
    const double amp2 = TestReggeRes(regge, lts, configured, gra::ReggeProductionModel::XP);
    CAPTURE(amp2);
    CHECK(std::isfinite(amp2));
    CHECK(amp2 > 0.0);
  }
}

// Check XP photoproduction uses the declared scalar fixed-spin Pomeron
TEST_CASE("TUNE0 XP photoproduction uses the scalar Pomeron",
          "[gra::MRegge][spin][XP][photoproduction][physics][regression]") {
  ToyHelicityProcess process;
  ConfigureToyProductionProcess(process, "XP", "RES", "pi+ pi-");
  const auto rho = gra::resonance::Read("RES/rho_770.json", process.state.random, gra::ReggeProductionModel::XP);
  process.SetResonances({{"rho_770", rho}});
  REQUIRE_NOTHROW(process.InitializeProcessAmplitude());

  gra::PARAM_RES configured = process.GetResonances().at("rho_770");
  REQUIRE(configured.production.size() == 2);
  for (const auto &channel : indices(configured.production)) {
    REQUIRE(configured.production[channel].tree.size() == 2);
    bool saw_photon  = false;
    bool saw_pomeron = false;
    for (const auto &leg : configured.production[channel].tree) {
      if (leg.p.pdg == gra::PDG::PDG_gamma) {
        CHECK(leg.hel.Jz_values == std::vector<double>{-1.0, 1.0});
        saw_photon = true;
        continue;
      }
      CHECK(leg.p.C == 1);
      CHECK(leg.p.pdg == 991);
      CHECK(leg.p.spinX2 == 0);
      saw_pomeron = true;
    }
    CHECK(saw_photon);
    CHECK(saw_pomeron);
    CHECK(configured.production[channel].pole.value().helicity.lambda_values.size_row() ==
          configured.production[channel].tree[0].hel.Jz_values.size() *
              configured.production[channel].tree[1].hel.Jz_values.size());
    const auto reduced = gra::spin::PoleLSReduced(configured.production[channel].pole.value(),
                                                  configured.production[channel].pole.value().Lambda);
    CHECK(reduced.FrobNorm2() > 0.0);
    const double reference_density = gra::spin::LeadingPoleDensity(configured.production[channel].pole.value(),
                                                                   configured.production[channel].pole.value().Lambda);
    CHECK(std::isfinite(reference_density));
    CHECK(reference_density > 0.0);
  }

  gra::LORENTZSCALAR lts     = MakeToyCoherentPhotonLTS();
  lts.process.MP_FRAME       = process.state.lts.process.MP_FRAME;
  lts.process.FORWARD_VERTEX = process.state.lts.process.FORWARD_VERTEX;
  lts.process.FORWARD_NOFLIP = process.state.lts.process.FORWARD_NOFLIP;
  gra::MRegge  regge(lts, process.state.model_tune,
                     gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Resonance, "XP_RES_photoproduction"));
  const double amp2 = TestReggeRes(regge, lts, configured, gra::ReggeProductionModel::XP);
  CAPTURE(amp2);
  CHECK(std::isfinite(amp2));
  CHECK(amp2 > 0.0);
}

TEST_CASE("MRegge GP resonances are invariant under common azimuth rotations", "[gra::MRegge][spin][GP]") {
  gra::LORENTZSCALAR base     = MakeToyReggeLTSAsymmetric(0.31, 4.2, -4.6);
  base.process.MMAX           = 2;
  base.process.FORWARD_NOFLIP = true;

  auto evaluate = [](gra::LORENTZSCALAR lts) {
    gra::MRegge    regge(lts, gra::MModelTune::Load(modelfile),
                         gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));
    gra::PARAM_RES res = MakeToyTensorGPResonance(lts.process.MMAX);
    TestReggeRes(regge, lts, res, gra::ReggeProductionModel::GP);
    return FiniteHelicityAmp2(lts.hamp);
  };

  const double reference = evaluate(base);
  REQUIRE(reference > 0.0);

  for (const double angle : {0.37, -1.21, 2.04}) {
    const gra::LORENTZSCALAR rotated = RotateToyEventAroundZ(base, angle);
    const double             actual  = evaluate(rotated);
    CAPTURE(angle, reference, actual);
    REQUIRE(actual / reference == Approx(1.0).epsilon(1e-10));
  }
}

// Check the strict forward pseudoscalar density from the (L,S)=(1,1) state
TEST_CASE("MRegge forward pseudoscalar density follows the sine squared law",
          "[gra::MRegge][spin][pseudoscalar][physics][regression]") {
  struct ModelCase {
    const char               *family;
    gra::ReggeProductionModel model;
  };
  constexpr std::array<ModelCase, 3> models = {{{"MP", gra::ReggeProductionModel::MP},
                                                {"XP", gra::ReggeProductionModel::XP},
                                                {"GP", gra::ReggeProductionModel::GP}}};

  for (const auto &test : models) {
    ToyHelicityProcess process;
    ConfigureToyProductionProcess(process, test.family, "RES", "gamma gamma");
    const auto eta = gra::resonance::Read("RES/eta.json", process.state.random, gra::ParseReggeProductionModel(test.family));
    process.SetResonances({{"eta", eta}});
    REQUIRE_NOTHROW(process.InitializeProcessAmplitude());

    gra::PARAM_RES pseudoscalar = process.state.lts.process.RESONANCES.at("eta");
    if (test.model != gra::ReggeProductionModel::GP) {
      PrepareToyPoleOperators(pseudoscalar, test.model, {{{1, 2, 1.0}}}, true, true);
      for (auto &production : pseudoscalar.production) {
        REQUIRE(production.pole.has_value());
        production.pole->derivative_factor = test.model == gra::ReggeProductionModel::XP;
      }
    }

    // Evaluate before the common line shape with higher fixed pole helicities
    // power suppressed by the strict forward transfer
    const auto evaluate = [&](const double dphi) {
      gra::LORENTZSCALAR lts         = MakeToyProductionLTSDphi(dphi);
      constexpr double   beam_pz     = 1000.0;
      constexpr double   proton_pz   = 998.0;
      constexpr double   proton_pt   = 0.002;
      const double       beam_energy = std::sqrt(gra::math::pow2(beam_pz) + gra::math::pow2(gra::PDG::mp));
      const double       proton_energy =
          std::sqrt(gra::math::pow2(proton_pz) + gra::math::pow2(proton_pt) + gra::math::pow2(gra::PDG::mp));
      lts.pbeam1    = gra::M4Vec(0.0, 0.0, beam_pz, beam_energy);
      lts.pbeam2    = gra::M4Vec(0.0, 0.0, -beam_pz, beam_energy);
      lts.pfinal[1] = gra::M4Vec(proton_pt, 0.0, proton_pz, proton_energy);
      lts.pfinal[2] = gra::M4Vec(proton_pt * std::cos(dphi), proton_pt * std::sin(dphi), -proton_pz, proton_energy);
      UpdateToyDerivedKinematics(lts);
      lts.process                = process.state.lts.process;
      lts.process.FORWARD_NOFLIP = true;
      lts.process.MP_FRAME       = "CM";
      lts.decaytree              = process.state.lts.decaytree;
      gra::MRegge                                     regge(lts, process.state.model_tune,
                                                            gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Resonance,
                                                                                              std::string(test.family) + "_pseudoscalar_dphi"));
      const auto                                      param = ReggeParametersForTest(regge, lts);
      std::vector<gra::MMatrix<std::complex<double>>> production;
      if (test.model == gra::ReggeProductionModel::GP) {
        production = gra::gpom::Resonance(lts, *param, pseudoscalar, nullptr);
      } else {
        production = gra::rspin::Resonance(
            lts, pseudoscalar,
            (test.model == gra::ReggeProductionModel::MP ? gra::mpom::Fusion : gra::xpom::Fusion), param->s0);
      }
      REQUIRE(production.size() == 1);
      return MatrixNorm2(production.front());
    };

    const double reference = evaluate(0.5 * gra::math::PI);
    REQUIRE(reference > 0.0);
    for (const double dphi : {0.0, gra::math::PI / 6.0, gra::math::PI / 3.0, 2.0 * gra::math::PI / 3.0,
                              5.0 * gra::math::PI / 6.0, gra::math::PI}) {
      const double ratio    = evaluate(dphi) / reference;
      const double expected = gra::math::pow2(std::sin(dphi));
      CAPTURE(test.family, dphi, ratio, expected);
      CHECK(ratio == Approx(expected).margin(1.0e-4));
    }
  }
}

// Check the finite-transfer spin-two harmonics and the bilinear RES fusion
TEST_CASE("MP and XP pseudoscalar poles retain the finite transfer node",
          "[gra::MRegge][MP][XP][spin][pseudoscalar][resonance][physics]"
          "[regression]") {
  constexpr double s0        = 1.0;
  constexpr double z         = 0.5;
  const double     proton_pt = std::sqrt(z * s0);
  constexpr double phi0      = 0.41;

  const auto event = [proton_pt, phi0](const double dphi) {
    gra::LORENTZSCALAR lts         = MakeToyProductionLTSDphi(dphi, phi0);
    constexpr double   beam_pz     = 100.0;
    constexpr double   proton_pz   = 98.0;
    const double       beam_energy = std::sqrt(gra::math::pow2(beam_pz) + gra::math::pow2(gra::PDG::mp));
    const double       proton_energy =
        std::sqrt(gra::math::pow2(proton_pz) + gra::math::pow2(proton_pt) + gra::math::pow2(gra::PDG::mp));
    lts.pbeam1    = gra::M4Vec(0.0, 0.0, beam_pz, beam_energy);
    lts.pbeam2    = gra::M4Vec(0.0, 0.0, -beam_pz, beam_energy);
    lts.pfinal[1] = gra::M4Vec(proton_pt * std::cos(phi0), proton_pt * std::sin(phi0), proton_pz, proton_energy);
    lts.pfinal[2] =
        gra::M4Vec(proton_pt * std::cos(phi0 + dphi), proton_pt * std::sin(phi0 + dphi), -proton_pz, proton_energy);
    UpdateToyDerivedKinematics(lts);
    lts.process.MP_FRAME       = "CM";
    lts.process.SPINGEN        = true;
    lts.process.FORWARD_NOFLIP = true;
    lts.process.FORWARD_VERTEX = gra::ForwardVertexMode::HelicityResidue;
    return lts;
  };

  for (const auto model : {gra::ReggeProductionModel::MP, gra::ReggeProductionModel::XP}) {
    CAPTURE(gra::ReggeProductionModelName(model));
    gra::PARAM_RES eta = MakeToyPseudoscalarTensorXP();
    PrepareToyPoleOperators(eta, model, {{{1, 2, 1.0}}}, true, true);
    eta.production.front().pole.value().derivative_factor = false;

    const auto amplitude = [&](const double dphi) {
      const auto production = gra::rspin::Resonance(
          event(dphi), eta, (model == gra::ReggeProductionModel::MP ? gra::mpom::Fusion : gra::xpom::Fusion),
          s0);
      REQUIRE(production.size() == 1);
      REQUIRE(production.front().size_col() == 1);
      REQUIRE(production.front().size_row() == 4);
      return production.front()[0][0];
    };
    const std::complex<double> reference = amplitude(gra::math::PI / 2.0);
    REQUIRE(std::abs(reference) > 1.0e-12);
    const double node = std::acos(1.0 / (4.0 * z));
    for (const double dphi : {0.0, gra::math::PI / 6.0, node, gra::math::PI / 3.0, gra::math::PI / 2.0,
                              2.0 * gra::math::PI / 3.0, 5.0 * gra::math::PI / 6.0, gra::math::PI}) {
      const double shape = (z * std::sin(dphi) - 2.0 * z * z * std::sin(2.0 * dphi)) / z;
      CAPTURE(dphi, shape);
      RequireComplexNear(amplitude(dphi) / reference, shape, 2.0e-10);
    }

    const gra::LORENTZSCALAR probe      = event(0.91);
    const auto              &tree       = eta.production.front().tree;
    const auto               up_rows    = gra::spin::Rows(tree[0], true);
    const auto               dn_rows    = gra::spin::Rows(tree[1], true);
    gra::M4Vec               lower_axis = probe.q2_in_X;
    lower_axis.Flip3();
    const gra::spin::ForwardSpec forward{gra::ForwardVertexMode::HelicityResidue, s0};
    const auto upper   = gra::spin::Forward(probe, tree[0], probe.pbeam1, probe.pfinal[1], false, up_rows,
                                            probe.process.PHOTON_VERTEX, forward);
    const auto lower   = gra::spin::Forward(probe, tree[1], probe.pbeam2, probe.pfinal[2], true, dn_rows,
                                            probe.process.PHOTON_VERTEX, forward);
    const auto central = (model == gra::ReggeProductionModel::MP ? gra::mpom::Fusion : gra::xpom::Fusion)(
        probe, eta.production.front().pole.value());
    const auto bilinear   = gra::spin::Contract(upper, lower, central);
    const auto wrong_dual = gra::spin::Contract(upper, lower.Conj(), central);
    const auto production = gra::rspin::Resonance(
        probe, eta, (model == gra::ReggeProductionModel::MP ? gra::mpom::Fusion : gra::xpom::Fusion), s0);
    REQUIRE(production.size() == 1);
    RequireMatrixNear(production.front(), bilinear, 2.0e-12);
    REQUIRE((production.front() - wrong_dual).FrobNorm2() > 1.0e-4 * production.front().FrobNorm2());
  }
}

// Check that analytic GP resonance fusion returns the matched spin-two pole
TEST_CASE("GP pseudoscalar fusion reaches the full fixed-spin pole matrix",
          "[gra::MRegge][MP][GP][spin][pseudoscalar][resonance][analytic]"
          "[regression]") {
  constexpr int      mmax        = 2;
  gra::LORENTZSCALAR base        = MakeToyProductionLTSDphi(0.31, 0.41);
  base.process.MMAX              = mmax;
  base.process.MP_FRAME          = "CM";
  base.process.SPINGEN           = true;
  base.process.FORWARD_NOFLIP    = true;
  base.process.FORWARD_VERTEX    = gra::ForwardVertexMode::HelicityResidue;
  base.process.DERIVATIVE_FACTOR = false;
  gra::MRegge  regge(base, gra::MModelTune::Load(modelfile),
                     gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));
  const auto   param         = ReggeParametersForTest(regge, base);
  const double pole_transfer = ReggeTransferForAlpha(*param, 990, 2.0);

  gra::PARAM_RES fixed = MakeToyPseudoscalarTensorXP();
  PrepareToyPoleOperators(fixed, gra::ReggeProductionModel::MP, {{{1, 2, 1.0}}}, true, true);
  fixed.production.front().pole.value().derivative_factor = false;

  gra::PARAM_RES analytic = MakeToyGPResonance(mmax);
  analytic.p              = fixed.p;
  gra::HELMatrix central;
  gra::gpom::InitHelicity(central, 2.0, 2.0, mmax, "GP pseudoscalar fixed pole test");
  central.BR             = 1.0;
  central.P_symmetry     = true;
  central.C_symmetry     = true;
  central.coupling_basis = gra::CouplingBasis::LS;
  central.alpha_ls.Set(1, 2, 1.0);
  gra::gpom::InitResonanceLS(central, 0, 4, 4);
  analytic.production.front().hel = std::move(central);

  for (const double dphi : {0.31, 0.91, 1.71, 2.47}) {
    gra::LORENTZSCALAR lts = MakeToyProductionLTSDphi(dphi, 0.41);
    lts.process            = base.process;
    UpdateToyDerivedKinematics(lts);
    lts.t1                     = pole_transfer;
    lts.t2                     = pole_transfer;
    const auto fixed_matrix    = gra::rspin::Resonance(lts, fixed, gra::mpom::Fusion, param->s0);
    const auto analytic_matrix = gra::gpom::Resonance(lts, *param, analytic, nullptr);
    REQUIRE(fixed_matrix.size() == 1);
    REQUIRE(analytic_matrix.size() == 1);
    REQUIRE(fixed_matrix.front().FrobNorm2() > 0.0);
    CAPTURE(dphi);
    RequireMatrixNear(analytic_matrix.front(), std::complex<double>(4.0, 0.0) * fixed_matrix.front(), 2.0e-11);
  }
}

// Check the real spherical section of either ordered exchange source
TEST_CASE("Forward exchange residues obey spherical reality",
          "[gra::MRegge][MP][XP][GP][spin][parity][source][regression]") {
  constexpr double qt            = 0.37;
  constexpr double phi           = 0.61;
  constexpr double exchange_mass = 0.83;
  constexpr double flip_mass     = 0.94;
  for (const bool use_exchange_barrier : {false, true}) {
    for (const bool reverse_phase : {false, true}) {
      for (int m = 1; m <= 6; ++m) {
        const auto positive =
            gra::spin::HelicityFactor(m, 0.0, qt, phi, exchange_mass, flip_mass, reverse_phase, use_exchange_barrier);
        const auto negative =
            gra::spin::HelicityFactor(-m, 0.0, qt, phi, exchange_mass, flip_mass, reverse_phase, use_exchange_barrier);
        const double section = m % 2 == 0 ? 1.0 : -1.0;
        CAPTURE(use_exchange_barrier, reverse_phase, m);
        RequireComplexNear(negative, section * std::conj(positive), 2.0e-14);
      }
    }
  }
}

// Check the exact forward axial-vector density for both canonical waves
TEST_CASE("MRegge forward axial-vector density follows its noflip law",
          "[gra::MRegge][spin][axial][physics][regression]") {
  struct ModelCase {
    const char               *family;
    gra::ReggeProductionModel model;
  };
  struct WaveCase {
    const char                    *label;
    std::vector<gra::spin::LSTerm> terms;
  };
  constexpr std::array<ModelCase, 3> models = {{{"MP", gra::ReggeProductionModel::MP},
                                                {"XP", gra::ReggeProductionModel::XP},
                                                {"GP", gra::ReggeProductionModel::GP}}};
  const std::array<WaveCase, 3>      waves  = {{
            {"(2,2)", {{2, 4, 1.0}}},
            {"(4,4)", {{4, 8, 1.0}}},
            {"complex (2,2)+(4,4)", {{2, 4, std::complex<double>(0.81, -0.27)}, {4, 8, std::complex<double>(-0.34, 0.52)}}},
  }};

  for (const auto &test : models) {
    ToyHelicityProcess process;
    ConfigureToyProductionProcess(process, test.family, "RES", "K*(892)+ K-");
    const auto f1 = gra::resonance::Read("RES/f1_1420.json", process.state.random, gra::ParseReggeProductionModel(test.family));
    process.SetResonances({{"f1_1420", f1}});
    REQUIRE_NOTHROW(process.InitializeProcessAmplitude());

    const gra::PARAM_RES configured = process.state.lts.process.RESONANCES.at("f1_1420");
    for (const auto &wave : waves) {
      CAPTURE(test.family, wave.label);
      gra::PARAM_RES axial = configured;
      if (test.model == gra::ReggeProductionModel::GP) {
        gra::HELMatrix central;
        gra::gpom::InitHelicity(central, 2.0, 2.0, process.state.lts.process.MMAX, "forward axial-vector density");
        central.BR             = 1.0;
        central.P_symmetry     = true;
        central.C_symmetry     = true;
        central.coupling_basis = gra::CouplingBasis::LS;
        for (const auto &term : wave.terms) { central.alpha_ls.Set(term.l, term.two_s, term.coefficient); }
        gra::gpom::InitResonanceLS(central, axial.p.spinX2 / 2, 4, 4);
        axial.production.front().hel = std::move(central);
      } else {
        PrepareToyPoleOperators(axial, test.model, {wave.terms}, true, true);
        for (auto &production : axial.production) {
          REQUIRE(production.pole.has_value());
          production.pole->derivative_factor = test.model == gra::ReggeProductionModel::XP;
        }
      }

      // Evaluate the forward production density before the common line shape
      const auto evaluate = [&](const double dphi, const double phi0 = 0.0, const double upper_pt = 0.002,
                                const double lower_pt = 0.002) {
        gra::LORENTZSCALAR lts         = MakeToyProductionLTSDphi(dphi, phi0);
        constexpr double   beam_pz     = 1000.0;
        constexpr double   proton_pz   = 998.0;
        const double       beam_energy = std::sqrt(gra::math::pow2(beam_pz) + gra::math::pow2(gra::PDG::mp));
        const double       upper_energy =
            std::sqrt(gra::math::pow2(proton_pz) + gra::math::pow2(upper_pt) + gra::math::pow2(gra::PDG::mp));
        const double lower_energy =
            std::sqrt(gra::math::pow2(proton_pz) + gra::math::pow2(lower_pt) + gra::math::pow2(gra::PDG::mp));
        lts.pbeam1    = gra::M4Vec(0.0, 0.0, beam_pz, beam_energy);
        lts.pbeam2    = gra::M4Vec(0.0, 0.0, -beam_pz, beam_energy);
        lts.pfinal[1] = gra::M4Vec(upper_pt * std::cos(phi0), upper_pt * std::sin(phi0), proton_pz, upper_energy);
        lts.pfinal[2] =
            gra::M4Vec(lower_pt * std::cos(phi0 + dphi), lower_pt * std::sin(phi0 + dphi), -proton_pz, lower_energy);
        UpdateToyDerivedKinematics(lts);
        lts.process                = process.state.lts.process;
        lts.process.FORWARD_NOFLIP = true;
        lts.process.MP_FRAME       = "CM";
        lts.decaytree              = process.state.lts.decaytree;
        gra::MRegge regge(
            lts, process.state.model_tune,
            gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Resonance, std::string(test.family) + "_axial_dphi"));
        const auto                                      param = ReggeParametersForTest(regge, lts);
        std::vector<gra::MMatrix<std::complex<double>>> production;
        if (test.model == gra::ReggeProductionModel::GP) {
          production = gra::gpom::Resonance(lts, *param, axial, nullptr);
        } else {
          production = gra::rspin::Resonance(
              lts, axial,
              (test.model == gra::ReggeProductionModel::MP ? gra::mpom::Fusion : gra::xpom::Fusion),
              param->s0);
        }
        REQUIRE(production.size() == 1);
        return production.front();
      };

      const double reference = MatrixNorm2(evaluate(gra::math::PI / 2.0));
      REQUIRE(reference > 0.0);
      for (const double dphi : {0.0, gra::math::PI / 6.0, gra::math::PI / 3.0, gra::math::PI / 2.0,
                                2.0 * gra::math::PI / 3.0, 5.0 * gra::math::PI / 6.0, gra::math::PI}) {
        const double ratio    = MatrixNorm2(evaluate(dphi)) / reference;
        const double expected = 1.0 - std::cos(dphi);
        CAPTURE(test.family, wave.label, dphi, ratio, expected);
        CHECK(ratio == Approx(expected).margin(1.0e-4));
      }

      // At either planar transfer the axial state is normal to the plane
      constexpr double           phi0     = 0.37;
      const std::complex<double> rotation = std::exp(-2.0 * gra::math::zi * phi0);
      const std::size_t          minus    = gra::spin::SpinProjectionIndex(-1.0, 1.0, "axial planar Jz=-1");
      const std::size_t          zero     = gra::spin::SpinProjectionIndex(0.0, 1.0, "axial planar Jz=0");
      const std::size_t          plus     = gra::spin::SpinProjectionIndex(1.0, 1.0, "axial planar Jz=+1");
      for (const double dphi : {0.0, gra::math::PI}) {
        const auto normal = evaluate(dphi, phi0, 0.0015, 0.0025);
        REQUIRE(normal.size_col() == 3);
        for (std::size_t row = 0; row < normal.size_row(); ++row) {
          const double transverse = std::max(std::abs(normal[row][minus]), std::abs(normal[row][plus]));
          REQUIRE(transverse > 1.0e-14);
          CAPTURE(test.family, wave.label, dphi, row, transverse);
          RequireComplexNear(normal[row][plus] / transverse, rotation * normal[row][minus] / transverse, 2.0e-8);
          CHECK(std::abs(normal[row][zero]) < 2.0e-8 * transverse);
        }
      }
    }
  }
}

// Check full parity mirrors with unchanged production and decay couplings
TEST_CASE("MP XP GP continua and resonances obey parity", "[gra::MRegge][spin][parity][physics][regression]") {
  struct ModelCase {
    const char               *family;
    gra::ReggeProductionModel model;
  };
  struct ResonanceCase {
    const char *label;
    const char *card;
    const char *decay;
    int         first_pdg;
    int         second_pdg;
  };
  constexpr std::array<ModelCase, 3>     models     = {{{"MP", gra::ReggeProductionModel::MP},
                                                        {"XP", gra::ReggeProductionModel::XP},
                                                        {"GP", gra::ReggeProductionModel::GP}}};
  constexpr std::array<ResonanceCase, 10> resonances = {{
      {"f0_500", "RES/f0_500.json", "pi+ pi-", 211, -211},
      {"eta", "RES/eta.json", "gamma gamma", 22, 22},
      {"rho_770", "RES/rho_770.json", "pi+ pi-", 211, -211},
      {"f1_1420", "RES/f1_1420.json", "K*(892)+ K-", 323, -321},
      {"f2_1270", "RES/f2_1270.json", "pi+ pi-", 211, -211},
      {"eta2_1645", "RES/eta2_1645.json", "a(2)(1320)0 pi0", 115, 111},
      {"rho3_1690", "RES/rho3_1690.json", "pi+ pi-", 211, -211},
      {"f4_2300", "RES/f4_2300.json", "pi+ pi-", 211, -211},
      {"rho5_2350", "RES/rho5_2350.json", "pi+ pi-", 211, -211},
      {"f6_2510", "RES/f6_2510.json", "pi+ pi-", 211, -211},
  }};

  for (const auto &test : models) {
    CAPTURE(test.family);

    for (const auto &[decay, pdgs] : std::array<std::pair<const char *, std::array<int, 2>>, 4>{{
             {"pi+ pi-", {211, -211}}, {"rho(770)0 rho(770)0", {113, 113}},
             {"p+ p-", {2212, -2212}}, {"phi(1020)0 phi(1020)0", {333, 333}}}}) {
      DYNAMIC_SECTION(test.family << " " << decay << " continuum") {
        ToyHelicityProcess process;
        ConfigureToyProductionProcess(process, test.family, "CON", decay);
        REQUIRE_NOTHROW(process.InitializeProcessAmplitude());

        gra::LORENTZSCALAR base = AsymmetricContinuumPairForTest(pdgs[0], pdgs[1]);
        base.process = process.state.lts.process;
        base.hamp.Configure(process.state.lts.hamp.metadata);
        const auto evaluate = [&](gra::LORENTZSCALAR lts) {
          gra::MRegge regge(lts, process.state.model_tune,
                            gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::ContinuumTwoBody,
                                                              std::string(test.family) + "_parity"));
          return TestReggeCon(regge, lts, test.model);
        };

        const double reference = evaluate(base);
        const double mirrored  = evaluate(ReflectToyEventInXZ(base));
        const double exchanged = evaluate(BeamExchangeMirrorWithDecay(base));
        CAPTURE(reference, mirrored, exchanged);
        REQUIRE(reference > 0.0);
        CHECK(mirrored == Approx(reference).epsilon(2.0e-10));
        CHECK(exchanged == Approx(reference).epsilon(2.0e-10));
        CHECK(evaluate(RotateToyEventAroundZ(base, 0.71)) == Approx(reference).epsilon(2.0e-10));
        CHECK(evaluate(BoostToyEventAlongZ(base, 0.43)) == Approx(reference).epsilon(2.0e-10));
      }
    }

    for (const auto &resonance : resonances) {
      DYNAMIC_SECTION(test.family << " " << resonance.label) {
        ToyHelicityProcess process;
        ConfigureToyProductionProcess(process, test.family, "RES", resonance.decay);
        const auto input = gra::resonance::Read(resonance.card, process.state.random, gra::ParseReggeProductionModel(test.family));
        process.SetResonances({{resonance.label, input}});
        REQUIRE_NOTHROW(process.InitializeProcessAmplitude());

        gra::LORENTZSCALAR base = DirectCentralPairLTSForTest(resonance.first_pdg, resonance.second_pdg);
        RefreshToyDerivedKinematicsPreserveDecay(base);
        base.process = process.state.lts.process;
        base.hamp.Configure(process.state.lts.hamp.metadata);
        const auto evaluate = [&](gra::LORENTZSCALAR lts) {
          gra::MRegge regge(
              lts, process.state.model_tune,
              gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Resonance, std::string(test.family) + "_parity"));
          return regge.Amp2(lts, test.model, gra::MReggeMode::Resonance);
        };

        const double reference = evaluate(base);
        const double mirrored  = evaluate(ReflectToyEventInXZ(base));
        CAPTURE(reference, mirrored);
        REQUIRE(std::isfinite(reference));
        REQUIRE(std::isfinite(mirrored));
        if (reference < 1.0e-24) {
          CHECK(std::abs(mirrored) < 1.0e-24);
        } else {
          CHECK(mirrored == Approx(reference).epsilon(2.0e-10));
        }
      }
    }
  }
}

// Check beam exchange with the same physical couplings in every model
TEST_CASE("MP XP GP amplitudes obey beam exchange",
          "[gra::MRegge][spin][beam-exchange][physics][regression]") {
  struct ModelCase {
    const char               *family;
    gra::ReggeProductionModel model;
  };
  constexpr std::array<ModelCase, 3> models = {{{"MP", gra::ReggeProductionModel::MP},
                                                {"XP", gra::ReggeProductionModel::XP},
                                                {"GP", gra::ReggeProductionModel::GP}}};

  for (const auto &test : models) {
    const std::vector<std::string> channels = {"CON", "RES", "RES+CON"};
    for (const std::string &channel : channels) {
      DYNAMIC_SECTION(test.family << " " << channel) {
        ToyHelicityProcess process;
        ConfigureToyProductionProcess(process, test.family, channel, "pi+ pi-");
        if (channel != "CON") {
          const auto f2 = gra::resonance::Read("RES/f2_1270.json", process.state.random, gra::ParseReggeProductionModel(test.family));
          process.SetResonances({{"f2_1270", f2}});
        }
        REQUIRE_NOTHROW(process.InitializeProcessAmplitude());

        const std::vector<std::string> frames = test.model == gra::ReggeProductionModel::MP
                                                    ? std::vector<std::string>{"CS", "HX", "CM"}
                                                    : std::vector<std::string>{"CM"};
        for (const auto &frame : frames) {
          CAPTURE(frame);
          gra::LORENTZSCALAR base = AsymmetricF2PhasePointForTest(0.91, -0.38);
          REQUIRE(std::abs(base.t1 - base.t2) > 1.0e-3);
          REQUIRE(base.pfinal[0].M() > 1.20);
          REQUIRE(base.pfinal[0].M() < 1.40);
          RefreshToyDerivedKinematicsPreserveDecay(base);
          base.process          = process.state.lts.process;
          base.process.MP_FRAME = frame;
          base.hamp.Configure(process.state.lts.hamp.metadata);
          const auto evaluate_mode = [&](gra::LORENTZSCALAR lts, gra::MReggeMode mode) {
            gra::MRegge  regge(lts, process.state.model_tune,
                               gra::MRegge::ProcessDefinitionFor(mode, std::string(test.family) + "_beam_exchange"));
            const double amp2 = regge.Amp2(lts, test.model, mode);
            return std::make_pair(amp2, lts.hamp);
          };
          const gra::MReggeMode mode = channel == "CON"   ? gra::MReggeMode::ContinuumTwoBody
                                       : channel == "RES" ? gra::MReggeMode::Resonance
                                                          : gra::MReggeMode::ResonanceContinuumTwoBody;

          const auto   reference_state = evaluate_mode(base, mode);
          const auto   mirrored_state  = evaluate_mode(BeamExchangeMirrorWithDecay(base), mode);
          const double reference       = reference_state.first;
          const double mirrored        = mirrored_state.first;
          if (mode == gra::MReggeMode::ResonanceContinuumTwoBody) {
            const auto               continuum_state     = evaluate_mode(base, gra::MReggeMode::ContinuumTwoBody);
            const auto               resonance_state     = evaluate_mode(base, gra::MReggeMode::Resonance);
            const gra::LORENTZSCALAR mirror              = BeamExchangeMirrorWithDecay(base);
            const auto               continuum_mirror    = evaluate_mode(mirror, gra::MReggeMode::ContinuumTwoBody);
            const auto               resonance_mirror    = evaluate_mode(mirror, gra::MReggeMode::Resonance);
            const double             continuum           = continuum_state.first;
            const double             resonance           = resonance_state.first;
            const double             incoherent          = continuum + resonance;
            const double             mirrored_continuum  = continuum_mirror.first;
            const double             mirrored_resonance  = resonance_mirror.first;
            const double             mirrored_incoherent = mirrored_continuum + mirrored_resonance;
            const double             average             = continuum_state.second.metadata.amplitude_normalization;
            REQUIRE(average > 0.0);
            REQUIRE(resonance_state.second.metadata.amplitude_normalization == Approx(average));
            REQUIRE(continuum_mirror.second.metadata.amplitude_normalization == Approx(average));
            REQUIRE(resonance_mirror.second.metadata.amplitude_normalization == Approx(average));
            REQUIRE(reference_state.second.metadata.amplitude_normalization == Approx(average));
            REQUIRE(mirrored_state.second.metadata.amplitude_normalization == Approx(average));
            std::complex<double> interference          = 0.0;
            std::complex<double> mirrored_interference = 0.0;
            REQUIRE(continuum_state.second.size() == resonance_state.second.size());
            REQUIRE(continuum_mirror.second.size() == resonance_mirror.second.size());
            for (const auto &i : indices(continuum_state.second)) {
              interference += continuum_state.second[i] * std::conj(resonance_state.second[i]);
              mirrored_interference += continuum_mirror.second[i] * std::conj(resonance_mirror.second[i]);
            }
            const double coherent          = incoherent + 2.0 * average * std::real(interference);
            const double mirrored_coherent = mirrored_incoherent + 2.0 * average * std::real(mirrored_interference);
            const double coherent_scale    = std::max({std::abs(reference), std::abs(coherent), 1.0e-30});
            const double mirrored_coherent_scale = std::max({std::abs(mirrored), std::abs(mirrored_coherent), 1.0e-30});
            const double overlap_scale = std::max({std::abs(interference), std::abs(mirrored_interference), 1.0e-30});
            const double interference_scale = std::max({std::abs(coherent), std::abs(incoherent), 1.0e-30});
            CAPTURE(reference, mirrored, continuum, resonance, mirrored_continuum, mirrored_resonance, incoherent,
                    mirrored_incoherent, coherent, mirrored_coherent, interference, mirrored_interference, average);
            CHECK(std::abs(reference - coherent) <= 2.0e-10 * coherent_scale);
            CHECK(std::abs(mirrored - mirrored_coherent) <= 2.0e-10 * mirrored_coherent_scale);
            CHECK(std::abs(mirrored_interference - interference) <= 2.0e-10 * overlap_scale);
            CHECK(std::abs(coherent - incoherent) > 1.0e-12 * interference_scale);
          }
          CAPTURE(reference, mirrored);
          REQUIRE(reference > 0.0);
          CHECK(mirrored == Approx(reference).epsilon(2.0e-10));
        }
      }
    }
  }
}

// Check that odd m continuum components preserve parity without changing couplings
TEST_CASE(
    "GP continuum preserves odd m interference under reflection and beam "
    "exchange",
    "[gra::MRegge][GP][CON][spin][parity][beam-exchange][regression]") {
  for (const std::string &channel : {std::string("CON"), std::string("RES+CON")}) {
    DYNAMIC_SECTION(channel) {
      ToyHelicityProcess process;
      ConfigureToyProductionProcess(process, "GP", channel, "pi+ pi-");
      if (channel != "CON") {
        const auto f2 = gra::resonance::Read("RES/f2_1270.json", process.state.random, gra::ReggeProductionModel::GP);
        process.SetResonances({{"f2_1270", f2}});
      }
      const auto tune = WriteModifiedContinuumTune("gp_odd_m_reflection", "GP", [](auto &card) {
        for (const std::string sector : {"same", "opposite"}) {
          card.at("990").at("[211,211]").at(sector).at("helicity") = {{0, 0, 0, 1.0, 0.0}, {0, 0, -1, 0.2, 0.0}};
        }
      });
      process.SetTuneForTest(tune);
      REQUIRE_NOTHROW(process.InitializeProcessAmplitude());

      gra::LORENTZSCALAR base = AsymmetricF2PhasePointForTest(0.91, -0.38);
      REQUIRE(std::abs(base.t1 - base.t2) > 1.0e-3);
      RefreshToyDerivedKinematicsPreserveDecay(base);
      base.process          = process.state.lts.process;
      base.process.MP_FRAME = "CM";
      base.hamp.Configure(process.state.lts.hamp.metadata);
      const gra::MReggeMode mode =
          channel == "CON" ? gra::MReggeMode::ContinuumTwoBody : gra::MReggeMode::ResonanceContinuumTwoBody;
      const auto evaluate = [&](gra::LORENTZSCALAR lts) {
        gra::MRegge regge(lts, process.state.model_tune,
                          gra::MRegge::ProcessDefinitionFor(mode, "GP_finite_m_section"));
        return regge.Amp2(lts, gra::ReggeProductionModel::GP, mode);
      };

      const double reference = evaluate(base);
      const double section   = evaluate(GPOddMSectionForTest(base));
      const double reflected = evaluate(ReflectToyEventInXZ(base));
      const double exchanged = evaluate(BeamExchangeMirrorWithDecay(base));
      CAPTURE(reference, section, reflected, exchanged);
      REQUIRE(reference > 0.0);
      REQUIRE(section > 0.0);
      CHECK(reflected == Approx(reference).epsilon(2.0e-10));
      CHECK(exchanged == Approx(reference).epsilon(2.0e-10));
      CHECK(std::abs(reference - section) > 1.0e-6 * reference);
    }
  }
}

// Keep fixed-C channel selection while allowing explicit interference signs
TEST_CASE("MRegge continuum t/u overrides retain C selection", "[gra::MRegge][spin][tu-sign]") {
  gra::LORENTZSCALAR pion_base = MakeToyScalarContinuumLTSAsymmetric(0.3, 4.2, -4.6);
  UpdateToyDerivedKinematics(pion_base);
  gra::MRegge regge(pion_base, gra::MModelTune::Load(modelfile),
                    gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));

  auto eval_pion = [&](std::string mode, int upper_exchange, int lower_exchange) {
    gra::LORENTZSCALAR lts = pion_base;
    lts.process.TU_SIGN    = std::move(mode);
    SetToyContinuumExchangePair(lts, upper_exchange, lower_exchange);
    TestReggeCon(regge, lts, gra::ReggeProductionModel::MP);
    return lts.hamp;
  };

  RequireVectorNear(eval_pion("auto", 9993, 993), eval_pion("positive", 9993, 993), 1e-11);
  REQUIRE_NOTHROW(eval_pion("negative", 9993, 993));
  REQUIRE_NOTHROW(eval_pion("negative", 993, 993));
  REQUIRE_THROWS_AS(eval_pion("invalid", 993, 993), std::invalid_argument);

  gra::LORENTZSCALAR rho_base = MakeToyReggeLTSAsymmetric(0.3, 4.2, -4.6);
  gra::MRegge rho_regge(rho_base, gra::MModelTune::Load(modelfile),
                       gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));
  auto               eval_rho = [&](std::string mode, int upper_exchange, int lower_exchange) {
    gra::LORENTZSCALAR lts = rho_base;
    lts.process.TU_SIGN    = std::move(mode);
    SetToyContinuumExchangePair(lts, upper_exchange, lower_exchange);
    TestReggeCon(rho_regge, lts, gra::ReggeProductionModel::MP);
    return lts.hamp;
  };

  RequireVectorNear(eval_rho("auto", 993, 993), eval_rho("positive", 993, 993), 1e-11);
  REQUIRE_NOTHROW(eval_rho("negative", 993, 993));
  for (const std::string mode : {"auto", "positive", "negative"}) {
    REQUIRE_THROWS_AS(eval_rho(mode, 9993, 993), std::invalid_argument);
  }
}

TEST_CASE(
    "MRegge GP continuum t/u projector rejects C-incompatible "
    "identical bosons",
    "[gra::MRegge][spin]") {
  gra::LORENTZSCALAR base = MakeToyReggeLTSAsymmetric(0.3, 4.2, -4.6);
  base.process.MMAX       = 1;
  SetToyContinuumExchangePair(base, 990, 990);
  UseToyGPSubchannelHelicity(base);
  gra::MRegge regge(base, gra::MModelTune::Load(modelfile),
                    gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));

  auto eval = [&](std::string mode, int upper_exchange, int lower_exchange) {
    gra::LORENTZSCALAR lts = base;
    lts.process.TU_SIGN    = std::move(mode);
    SetToyContinuumExchangePair(lts, upper_exchange, lower_exchange);
    UseToyGPSubchannelHelicity(lts);
    TestReggeCon(regge, lts, gra::ReggeProductionModel::GP);
    return lts.hamp;
  };

  RequireVectorNear(eval("auto", 990, 990), eval("positive", 990, 990), 1e-11);
  REQUIRE_THROWS(eval("auto", 9930, 990));
}

TEST_CASE(
    "MRegge phase algebra and LOOPSCREEN convolution agree for every spin "
    "route",
    "[gra::MRegge][MProcess][screening][phase][regression]") {
  const double rotation = 0.71;
  struct ReggePhaseCase {
    const char               *label;
    gra::ReggeProductionModel spin;
    gra::MReggeMode           mode;
  };
  const std::array<ReggePhaseCase, 6> cases = {{
      {"MP_RES", gra::ReggeProductionModel::MP, gra::MReggeMode::Resonance},
      {"XP_RES", gra::ReggeProductionModel::XP, gra::MReggeMode::Resonance},
      {"GP_RES", gra::ReggeProductionModel::GP, gra::MReggeMode::Resonance},
      {"MP_CON", gra::ReggeProductionModel::MP, gra::MReggeMode::ContinuumTwoBody},
      {"XP_CON", gra::ReggeProductionModel::XP, gra::MReggeMode::ContinuumTwoBody},
      {"GP_CON", gra::ReggeProductionModel::GP, gra::MReggeMode::ContinuumTwoBody},
  }};

  const double screening_s = MakeToyCoherentPhotonLTS().s;
  for (const auto &test : cases) {
    CAPTURE(test.label);
    const auto path = WriteSoftModelFile({}, "regge_phase_" + std::string(test.label), 0.25, 1, 0.0, 0.0);
    for (const std::string model : {"MP", "XP", "GP"}) {
      const auto filename = std::filesystem::path(path).parent_path() / ("CON_" + model + ".json");
      auto card = nlohmann::json::parse(gra::aux::GetInputData(filename.string()));
      SetToyPhotonContinuum(card);
      std::ofstream output(filename);
      REQUIRE(output.good());
      output << card.dump(2);
    }
    const auto model_tune = gra::MModelTune::Load(path);
    gra::MEikonal eikonal(model_tune);
    eikonal.S3Constructor(screening_s, ProtonInitialState(), false, 4, 4);

    const std::vector<std::string>       frames = test.spin == gra::ReggeProductionModel::MP
                                                      ? std::vector<std::string>{"CS", "HX"}
                                                      : std::vector<std::string>{"CM"};
    const std::vector<ToyPhotonTopology> topologies =
        test.mode == gra::MReggeMode::Resonance
            ? std::vector<ToyPhotonTopology>{ToyPhotonTopology::Upper, ToyPhotonTopology::Lower,
                                             ToyPhotonTopology::Coherent}
            : std::vector<ToyPhotonTopology>{ToyPhotonTopology::Hadronic, ToyPhotonTopology::Upper,
                                             ToyPhotonTopology::Lower, ToyPhotonTopology::Coherent};
    for (const auto &frame : frames) {
      for (const auto topology : topologies) {
        CAPTURE(frame, topology);
        ToyReggePhaseScreeningProcess born(model_tune, test.spin, test.mode, frame, 0.0, topology);
        ToyReggePhaseScreeningProcess rotated(model_tune, test.spin, test.mode, frame, rotation, topology);
        const double                  born_amp2    = born.BornAmp2();
        const double                  rotated_amp2 = rotated.BornAmp2();
        REQUIRE_FALSE(born.state.lts.hamp.empty());
        RequireAzimuthalCovariance(rotated.state.lts.hamp, born.state.lts.hamp,
                                   std::vector<int>(born.state.lts.hamp.size(), 0), rotation, 2.0e-10);
        REQUIRE(rotated_amp2 == Approx(born_amp2).epsilon(2.0e-10).margin(2.0e-10));
      }
    }

    const std::string             screening_frame    = test.spin == gra::ReggeProductionModel::MP ? "CS" : "CM";
    const ToyPhotonTopology       screening_topology = ToyPhotonTopology::Coherent;
    ToyReggePhaseScreeningProcess unscreened(model_tune, test.spin, test.mode, screening_frame, 0.0,
                                             screening_topology);
    const double                  unscreened_amp2 = unscreened.BornAmp2();
    ToyReggePhaseScreeningProcess screened(model_tune, test.spin, test.mode, screening_frame, 0.0, screening_topology);
    screened.eikonal                                  = eikonal;
    screened.eikonal.Numerics.LOOP.radial_integrator  = "GL";
    screened.eikonal.Numerics.LOOP.azimuth_integrator = "Trap";
    screened.eikonal.Numerics.LOOP.radial_map         = gra::math::RadialMap::Linear;
    screened.eikonal.Numerics.LOOP.r_min              = 0.0;
    screened.eikonal.Numerics.LOOP.r_max              = 0.03;
    screened.eikonal.Numerics.LOOP.radial_intervals   = 2;
    screened.eikonal.Numerics.LOOP.azimuth_nodes      = 4;
    screened.eikonal.InitLoopWeightMatrix();

    const auto  &loop          = screened.eikonal.GetLoopConst(screening_s);
    const auto  &initial_state = screened.eikonal.InitialState();
    const double beta          = gra::kinematics::beta12(screening_s, initial_state[0].mass, initial_state[1].mass);
    RequireComplexNear(loop.norm, gra::math::zi / (8.0 * gra::math::PIPI * screening_s * beta), 1.0e-14);
    const double screened_amp2 = screened.ScreenedAmp2();
    const auto   reference     = ManualReggeScreeningFromTrace(screened.PairTrace(), screened.eikonal);
    CAPTURE(screened_amp2, reference.born, reference.interference, reference.loop_squared);
    RequireVectorNear(screened.state.lts.hamp, reference.amplitude, 2.0e-10);
    const double expanded_amp2 = reference.born + reference.interference + reference.loop_squared;
    REQUIRE(reference.born == Approx(unscreened_amp2).epsilon(2.0e-10).margin(2.0e-10));
    REQUIRE(screened_amp2 == Approx(expanded_amp2).epsilon(2.0e-10).margin(2.0e-10));
    REQUIRE(reference.interference < 0.0);
    REQUIRE(screened_amp2 < reference.born);

    // A rotation by one angular node spacing permutes the complete loop rule
    const double                  loop_rotation = 0.5 * gra::math::PI;
    ToyReggePhaseScreeningProcess rotated_screened(model_tune, test.spin, test.mode, screening_frame, loop_rotation,
                                                   screening_topology);
    rotated_screened.eikonal                                  = eikonal;
    rotated_screened.eikonal.Numerics.LOOP.radial_integrator  = "GL";
    rotated_screened.eikonal.Numerics.LOOP.azimuth_integrator = "Trap";
    rotated_screened.eikonal.Numerics.LOOP.radial_map         = gra::math::RadialMap::Linear;
    rotated_screened.eikonal.Numerics.LOOP.r_min              = 0.0;
    rotated_screened.eikonal.Numerics.LOOP.r_max              = 0.03;
    rotated_screened.eikonal.Numerics.LOOP.radial_intervals   = 2;
    rotated_screened.eikonal.Numerics.LOOP.azimuth_nodes      = 4;
    rotated_screened.eikonal.InitLoopWeightMatrix();
    const double rotated_screened_amp2 = rotated_screened.ScreenedAmp2();
    RequireAzimuthalCovariance(rotated_screened.state.lts.hamp, screened.state.lts.hamp,
                               std::vector<int>(screened.state.lts.hamp.size(), 0), loop_rotation, 3.0e-10);
    REQUIRE(rotated_screened_amp2 == Approx(screened_amp2).epsilon(3.0e-10).margin(3.0e-10));

    // A common resonance coupling phase commutes with the linear screening map
    const std::complex<double> phasor        = std::polar(1.0, 0.613);
    auto                       rotated_trace = screened.PairTrace();
    for (auto &entry : rotated_trace) {
      for (auto &component : entry.components) { component.source *= phasor; }
    }
    const auto rotated = ManualReggeScreeningFromTrace(rotated_trace, screened.eikonal);
    REQUIRE(rotated.born == Approx(reference.born).epsilon(2.0e-11));
    REQUIRE(rotated.interference == Approx(reference.interference).epsilon(2.0e-11));
    REQUIRE(rotated.loop_squared == Approx(reference.loop_squared).epsilon(2.0e-11));
    for (const auto &i : indices(reference.amplitude)) {
      RequireComplexNear(rotated.amplitude[i], phasor * reference.amplitude[i], 2.0e-10);
    }

    // A production mass factor commutes with the full screening integral at fixed central mass
    if (test.mode == gra::MReggeMode::Resonance) {
      REQUIRE(screened.state.lts.process.RESONANCES.size() == 1);
      auto &res = screened.state.lts.process.RESONANCES.begin()->second;
      auto &form = ReggeResonanceFormForTest(res, test.spin);
      const double delta = screened.state.lts.m2 - gra::math::pow2(res.p.mass);
      const double previous = form.ff_prod.param[0];
      form.ff_prod.param[0] = 1.3 * previous;
      const double factor = std::exp(gra::math::pow2(delta / previous) -
                                     gra::math::pow2(delta / form.ff_prod.param[0]));
      CHECK(screened.ScreenedAmp2() == Approx(gra::math::pow2(factor) * screened_amp2).epsilon(3.0e-10));
      for (const auto &i : indices(reference.amplitude)) {
        RequireComplexNear(screened.state.lts.hamp[i], factor * reference.amplitude[i], 3.0e-10);
      }
    }
  }
}

TEST_CASE(
    "MRegge Pomeron Reggeon and Odderon screening equals an independent "
    "N2 contraction",
    "[gra::MRegge][MProcess][screening][GoodWalker][multichannel]") {
  const double screening_s = MakeToyCoherentPhotonLTS().s;
  auto eikonal = BuildTestEikonalWithInitialState({0.31}, "regge_n2_pro_screening", ProtonInitialState(), 0.25, 1, 4, 4,
                                                  screening_s, 0.13);
  REQUIRE(eikonal.GetChannelCount() == 2);

  const std::vector<std::array<int, 2>> hard_channels = {{993, 993}, {9915, 993}, {9993, 993}};
  ToyReggePhaseScreeningProcess         process(eikonal.ModelTuneHandle(), gra::ReggeProductionModel::MP,
                                                gra::MReggeMode::ContinuumTwoBody, "CS", 0.0, ToyPhotonTopology::Hadronic,
                                                hard_channels);
  const double                          born_amp2 = process.BornAmp2();
  REQUIRE(std::isfinite(born_amp2));
  REQUIRE(born_amp2 > 0.0);
  REQUIRE(process.PairTrace().size() == 1);
  const auto &coherent_source = process.PairTrace().front().components.front().source;
  REQUIRE(coherent_source.size_col() == 4);

  gra::MMatrix<std::complex<double>> source_sum(coherent_source.size_row(), coherent_source.size_col(), 0.0);
  for (const auto &channel : hard_channels) {
    ToyReggePhaseScreeningProcess isolated(eikonal.ModelTuneHandle(), gra::ReggeProductionModel::MP,
                                           gra::MReggeMode::ContinuumTwoBody, "CS", 0.0, ToyPhotonTopology::Hadronic,
                                           {channel});
    REQUIRE(isolated.BornAmp2() > 0.0);
    REQUIRE(isolated.PairTrace().size() == 1);
    const auto &source = isolated.PairTrace().front().components.front().source;
    REQUIRE(source.size_row() == source_sum.size_row());
    REQUIRE(source.size_col() == source_sum.size_col());
    source_sum += source;
  }
  RequireMatrixNear(coherent_source, source_sum, 2.0e-10);

  process.eikonal                                  = eikonal;
  process.eikonal.Numerics.LOOP.radial_integrator  = "GL";
  process.eikonal.Numerics.LOOP.azimuth_integrator = "Trap";
  process.eikonal.Numerics.LOOP.radial_map         = gra::math::RadialMap::Linear;
  process.eikonal.Numerics.LOOP.r_min              = 0.0;
  process.eikonal.Numerics.LOOP.r_max              = 0.03;
  process.eikonal.Numerics.LOOP.radial_intervals   = 2;
  process.eikonal.Numerics.LOOP.azimuth_nodes      = 4;
  process.eikonal.InitLoopWeightMatrix();

  const double screened_amp2 = process.ScreenedAmp2();
  const auto   reference     = ManualReggeScreeningFromTrace(process.PairTrace(), process.eikonal);
  RequireVectorNear(process.state.lts.hamp, reference.amplitude, 3.0e-10);
  const double expanded = reference.born + reference.interference + reference.loop_squared;
  REQUIRE(reference.born == Approx(born_amp2).epsilon(3.0e-10).margin(1.0e-12));
  REQUIRE(screened_amp2 == Approx(expanded).epsilon(3.0e-10).margin(1.0e-12));
  REQUIRE(std::isfinite(screened_amp2));
  REQUIRE(std::abs(screened_amp2 - born_amp2) > 1.0e-12 * std::max(1.0, born_amp2));
}

TEST_CASE("MP XP GP and Tensor Pomeron share Breit-Wigner phase motion",
          "[gra::MRegge][MTensorPomeron][phase][interference][physics]") {
  const auto                  tensor_tune = WriteModifiedPhotoVMTune("tensor_pomeron_pole_phase", [](auto &j) {
    for (auto &[key, pairs] : j.at("PARAM_TENSORPOM").at("PARAM_CON").at("TP").items()) {
      (void)key;
      pairs = {{995, 995}};
    }
  });
  const double                pole_mass   = 1.2;
  const double                width       = 0.002;
  const std::array<double, 3> mass2       = {gra::math::pow2(pole_mass) - pole_mass * width, gra::math::pow2(pole_mass),
                                             gra::math::pow2(pole_mass) + pole_mass * width};

  // Build one physical pion-pair point at a selected central invariant mass
  const auto make_lts = [](double central_mass) {
    gra::LORENTZSCALAR lts           = MakeToyScalarContinuumLTSAsymmetric(0.2, 4.2, -4.2);
    const double       beam_pz       = 5.0;
    const double       beam_energy   = std::sqrt(gra::math::pow2(beam_pz) + gra::math::pow2(gra::PDG::mp));
    const double       proton_energy = beam_energy - 0.5 * central_mass;
    const double       proton_px     = 0.16;
    const double       proton_py     = 0.07;
    const double       proton_pz     = std::sqrt(gra::math::pow2(proton_energy) - gra::math::pow2(gra::PDG::mp) -
                                                 gra::math::pow2(proton_px) - gra::math::pow2(proton_py));
    lts.pbeam1                       = gra::M4Vec(0.0, 0.0, beam_pz, beam_energy);
    lts.pbeam2                       = gra::M4Vec(0.0, 0.0, -beam_pz, beam_energy);
    lts.pfinal[1]                    = gra::M4Vec(proton_px, proton_py, proton_pz, proton_energy);
    lts.pfinal[2]                    = gra::M4Vec(-proton_px, -proton_py, -proton_pz, proton_energy);
    lts.pfinal[0]                    = lts.pbeam1 + lts.pbeam2 - lts.pfinal[1] - lts.pfinal[2];
    const auto pair =
        TwoBodyRestKinematics(central_mass, lts.decaytree[0].p.mass, lts.decaytree[1].p.mass, 0.91, -0.38);
    lts.decaytree[0].p4        = pair[0];
    lts.decaytree[1].p4        = pair[1];
    lts.process.MP_FRAME       = "CM";
    lts.process.FORWARD_NOFLIP = true;
    lts.process.TU_SIGN        = "positive";
    RefreshToyDerivedKinematicsPreserveDecay(lts);
    return lts;
  };

  // Compute the basis-independent resonance-continuum interference overlap
  const auto overlap = [](const auto &resonance, const auto &continuum) {
    REQUIRE(resonance.size() == continuum.size());
    std::complex<double> out = 0.0;
    for (const auto &i : indices(resonance)) { out += resonance[i] * std::conj(continuum[i]); }
    REQUIRE(std::abs(out) > 0.0);
    return out;
  };

  // Evaluate the physical MP, XP, or GP resonance against its continuum
  const auto regge_motion = [&](gra::ReggeProductionModel model, const double phi = 0.0) {
    std::array<std::complex<double>, 3> interference{};
    for (const auto &i : indices(mass2)) {
      gra::LORENTZSCALAR lts       = make_lts(std::sqrt(mass2[i]));
      gra::PARAM_RES     resonance = MakeToyScalarMPResonance();
      if (model == gra::ReggeProductionModel::XP) {
        PrepareToyXPOperators(resonance, {{{0, 0, 1.0}}});
        PrepareToyXPContinuumOperators(lts, {{0, 0, 1.0}});
      } else if (model == gra::ReggeProductionModel::GP) {
        resonance = MakeToyScalarGPResonance(lts.process.MMAX);
        SetToyContinuumExchangePair(lts, 990, 990);
        UseToyGPSubchannelHelicity(lts);
      }
      resonance.p.mass            = pole_mass;
      resonance.p.width           = width;
      auto &form                  = ReggeResonanceFormForTest(resonance, model);
      form.ff_prod                = {};
      resonance.hel_decay.g_decay = 1.0;
      if (model == gra::ReggeProductionModel::MP) {
        resonance.MP.phi = phi;
      } else if (model == gra::ReggeProductionModel::XP) {
        resonance.XP.phi = phi;
      } else {
        resonance.GP.phi = phi;
      }

      gra::MRegge regge(lts, gra::MModelTune::Load(modelfile),
                        gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "pole_phase_scan"));
      TestReggeCon(regge, lts, model);
      const auto continuum = lts.hamp;
      TestReggeRes(regge, lts, resonance, model);
      interference[i] = overlap(lts.hamp, continuum);
    }
    return interference;
  };

  // Evaluate the same overlap with complete Tensor Pomeron graphs
  const auto tensor_motion = [&](const double phi = 0.0) {
    std::array<std::complex<double>, 3> interference{};
    const auto                          soft_model = gra::MModelTune::Load(tensor_tune.second);
    for (const auto &i : indices(mass2)) {
      gra::LORENTZSCALAR lts = make_lts(std::sqrt(mass2[i]));
      gra::PARAM_RES     resonance;
      resonance.p       = ToyParticle("f0_phase", 9001710, 0, pole_mass);
      resonance.p.P     = 1;
      resonance.p.C     = 1;
      resonance.p.width = width;
      SetToyTensorChannel(resonance, {0.37, -0.21});
      resonance.TP.phi                   = phi;
      resonance.hel_decay.g_decay_TP = {0.29};
      lts.process.RESONANCES             = {{"f0_phase", resonance}};

      gra::MTensorPomeron tensor(lts, soft_model,
                                 gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
      tensor.ME3(lts);
      const auto resonance_amplitude = lts.hamp;
      tensor.ME4(lts, gra::TensorContinuumMode::TensorPomeron);
      interference[i] = overlap(resonance_amplitude, lts.hamp);
    }
    return interference;
  };

  // Remove each model's constant coupling and continuum phase at the resonance
  const auto require_bw_motion = [](const std::string &label, const auto &interference) {
    const double below = std::arg(interference[0] / interference[1]);
    const double pole  = std::arg(interference[1] / interference[1]);
    const double above = std::arg(interference[2] / interference[1]);
    CAPTURE(label, below, pole, above);
    CHECK(below == Approx(-0.25 * gra::math::PI).margin(0.08));
    CHECK(pole == Approx(0.0).margin(1.0e-14));
    CHECK(above == Approx(0.25 * gra::math::PI).margin(0.08));
  };

  const auto mp = regge_motion(gra::ReggeProductionModel::MP);
  const auto xp = regge_motion(gra::ReggeProductionModel::XP);
  const auto gp = regge_motion(gra::ReggeProductionModel::GP);
  const auto tp = tensor_motion();
  require_bw_motion("MP", mp);
  require_bw_motion("XP", xp);
  require_bw_motion("GP", gp);
  require_bw_motion("TP", tp);

  // A model coupling phase rotates every mass point without changing intensity
  const double               phi              = 0.731;
  const std::complex<double> phasor           = std::polar(1.0, phi);
  const auto                 require_rotation = [&](const auto &reference, const auto &rotated) {
    for (const auto &i : indices(reference)) {
      RequireComplexNear(rotated[i], phasor * reference[i], 2.0e-11);
      CHECK(std::norm(rotated[i]) == Approx(std::norm(reference[i])).epsilon(2.0e-11));
    }
  };
  require_rotation(mp, regge_motion(gra::ReggeProductionModel::MP, phi));
  require_rotation(xp, regge_motion(gra::ReggeProductionModel::XP, phi));
  require_rotation(gp, regge_motion(gra::ReggeProductionModel::GP, phi));
  require_rotation(tp, tensor_motion(phi));
}

TEST_CASE("TP MP and XP resolve the same resonance and continuum phases",
          "[gra::MRegge][MTensorPomeron][continuum][resonance][phase]"
          "[physics]") {
  const auto tensor_tune  = WriteModifiedPhotoVMTune("tensor_pomeron_resonance_continuum_phase", [](auto &j) {
    for (auto &[key, pairs] : j.at("PARAM_TENSORPOM").at("PARAM_CON").at("TP").items()) {
      (void)key;
      pairs = {{995, 995}};
    }
  });
  const auto model_tune   = gra::MModelTune::Load(tensor_tune.second);
  const auto tensor_param = gra::ReadTensorPomeronParam(*model_tune, LoadedPDGTable());

  struct ContinuumCase {
    const char               *family;
    gra::ReggeProductionModel model;
    int                       pomeron_pdg;
  };
  const std::array<ContinuumCase, 2> cases = {{
      {"MP", gra::ReggeProductionModel::MP, 995},
      {"XP", gra::ReggeProductionModel::XP, 995},
  }};

  // Initialize one common card-loaded scalar pole in every model
  ToyHelicityProcess tensor_process;
  ConfigureToyProductionProcess(tensor_process, "TP", "RES+CON", "pi+ pi-");
  const auto tensor_input = gra::resonance::Read("RES/f0_500.json", tensor_process.state.random, gra::ReggeProductionModel::TP);
  tensor_process.SetResonances({{"f0_500", tensor_input}});
  REQUIRE_NOTHROW(tensor_process.InitializeProcessAmplitude());
  const auto                              tensor_resonance = tensor_process.GetResonances().at("f0_500");
  const double                            pole_mass        = tensor_resonance.p.mass;
  const std::array<gra::LORENTZSCALAR, 3> points           = {ScalarPolePhasePointForTest(0.12, 0.61, -0.72, pole_mass),
                                                              ScalarPolePhasePointForTest(0.20, 1.13, 0.38, pole_mass),
                                                              ScalarPolePhasePointForTest(0.31, 2.17, 1.29, pole_mass)};

  struct TensorPhasePoint {
    std::vector<std::complex<double>> resonance;
    std::vector<std::complex<double>> t;
    std::vector<std::complex<double>> u;
    std::complex<double>              forward_t;
    std::complex<double>              forward_u;
  };
  std::array<TensorPhasePoint, 3> tensor_points;
  for (const auto &point_index : indices(points)) {
    gra::LORENTZSCALAR lts = points[point_index];
    lts.process            = tensor_process.state.lts.process;
    lts.process.RESONANCES = {{"f0_500", tensor_resonance}};
    gra::MTensorPomeron tensor(lts, model_tune,
                               gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::ResonanceContinuum));
    tensor.ME3(lts);
    tensor_points[point_index].resonance = lts.hamp;
    tensor.ME4(lts, gra::TensorContinuumMode::TensorPomeron);
    const auto  full = lts.hamp;
    const auto  tu   = TensorPionContinuumTUForTest(tensor, *tensor_param, lts);
    const auto &p1   = lts.pfinal[1];
    const auto &p2   = lts.pfinal[2];
    const auto &p3   = lts.decaytree[0].p4;
    const auto &p4   = lts.decaytree[1].p4;
    tensor_points[point_index].forward_t =
        tensor.PomeronPropagatorFactor((p1 + p3).M2(), lts.t1) * tensor.PomeronPropagatorFactor((p2 + p4).M2(), lts.t2);
    tensor_points[point_index].forward_u =
        tensor.PomeronPropagatorFactor((p1 + p4).M2(), lts.t1) * tensor.PomeronPropagatorFactor((p2 + p3).M2(), lts.t2);
    REQUIRE(std::abs(tu[0] + tu[1]) > 0.0);
    const auto active = std::max_element(
        full.begin(), full.end(), [](const auto &left, const auto &right) { return std::abs(left) < std::abs(right); });
    REQUIRE(active != full.end());
    RequireComplexNear(*active, tu[0] + tu[1], 2.0e-11);
    tensor_points[point_index].t = full;
    tensor_points[point_index].u = full;
    for (const auto &i : indices(full)) {
      tensor_points[point_index].t[i] *= tu[0] / (tu[0] + tu[1]);
      tensor_points[point_index].u[i] *= tu[1] / (tu[0] + tu[1]);
    }
    RequireVectorNear(
        LinearAmplitudeSectionForTest(tensor_points[point_index].t, tensor_points[point_index].u, 1.0, 1.0), full,
        2.0e-11);
    REQUIRE(tensor_points[point_index].resonance.size() == full.size());
    const std::complex<double> tu_overlap =
        ComplexSectionOverlapForTest(tensor_points[point_index].t, tensor_points[point_index].u);
    const std::complex<double> rt_overlap =
        ComplexSectionOverlapForTest(tensor_points[point_index].resonance, tensor_points[point_index].t);
    const std::complex<double> ru_overlap =
        ComplexSectionOverlapForTest(tensor_points[point_index].resonance, tensor_points[point_index].u);
    CAPTURE(point_index, tu_overlap, rt_overlap, ru_overlap);
    CHECK(std::abs(tu_overlap.imag()) <= 2.0e-10 * std::abs(tu_overlap));
    CHECK(tu_overlap.real() > 0.0);
    CHECK(std::abs(rt_overlap.real()) <= 2.0e-10 * std::abs(rt_overlap));
    CHECK(std::abs(ru_overlap.real()) <= 2.0e-10 * std::abs(ru_overlap));
    CHECK(rt_overlap.imag() > 0.0);
    CHECK(ru_overlap.imag() > 0.0);
  }

  for (const auto &test : cases) {
    CAPTURE(test.family, test.pomeron_pdg);
    ToyHelicityProcess process;
    ConfigureToyProductionProcess(process, test.family, "RES+CON", "pi+ pi-");
    auto input = gra::resonance::Read("RES/f0_500.json", process.state.random, gra::ParseReggeProductionModel(test.family));
    // Compare production phases at the same physical pole in every model
    input.p = tensor_resonance.p;
    input.MP.phi = 0.0;
    input.XP.phi = 0.0;
    for (auto &channel : input.MP.channels) { channel.g = 1.0; }
    for (auto &channel : input.XP.channels) {
      for (auto &term : channel.g_ls) { term.coefficient = term.l == 0 ? 1.0 : 0.0; }
    }
    process.SetResonances({{"f0_500", input}});
    REQUIRE_NOTHROW(process.InitializeProcessAmplitude());
    const auto &configured = process.state.lts.process;
    const auto  channel    = std::find(configured.CONT_PRODUCTION.begin(), configured.CONT_PRODUCTION.end(),
                                       std::vector<int>{test.pomeron_pdg, test.pomeron_pdg});
    REQUIRE(channel != configured.CONT_PRODUCTION.end());
    const std::size_t channel_index =
        static_cast<std::size_t>(std::distance(configured.CONT_PRODUCTION.begin(), channel));

    for (const auto &point_index : indices(points)) {
      const auto evaluate_continuum = [&](double sign) {
        gra::LORENTZSCALAR lts          = points[point_index];
        lts.process                     = configured;
        lts.process.FORWARD_NOFLIP      = true;
        lts.process.TU_SIGN             = "positive";
        lts.process.CONT_PRODUCTION     = {configured.CONT_PRODUCTION[channel_index]};
        lts.process.CONT_PRODUCTIONTREE = {configured.CONT_PRODUCTIONTREE[channel_index]};
        lts.process.CONTINUUM_GP.clear();
        lts.process.CONTINUUM_POLE = {configured.CONTINUUM_POLE[channel_index]};
        lts.process.CONT_TU_SIGN   = {sign};
        gra::MRegge regge(lts, process.state.model_tune,
                          gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::ContinuumTwoBody,
                                                            std::string(test.family) + "_resolved_CON"));
        TestReggeCon(regge, lts, test.model);
        return lts.hamp;
      };
      const auto plus  = evaluate_continuum(1.0);
      const auto minus = evaluate_continuum(-1.0);
      const auto t     = LinearAmplitudeSectionForTest(plus, minus, 0.5, 0.5);
      const auto u     = LinearAmplitudeSectionForTest(plus, minus, 0.5, -0.5);

      gra::LORENTZSCALAR lts   = points[point_index];
      lts.process              = configured;
      gra::PARAM_RES resonance = process.GetResonances().at("f0_500");
      gra::MRegge    regge(
             lts, process.state.model_tune,
             gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Resonance, std::string(test.family) + "_resolved_RES"));
      TestReggeRes(regge, lts, resonance, test.model);
      const auto                &tensor          = tensor_points[point_index];
      const std::complex<double> regge_tu        = ComplexSectionOverlapForTest(t, u);
      const std::complex<double> tensor_tu       = ComplexSectionOverlapForTest(tensor.t, tensor.u);
      const std::complex<double> regge_rt        = ComplexSectionOverlapForTest(lts.hamp, t);
      const std::complex<double> tensor_rt       = ComplexSectionOverlapForTest(tensor.resonance, tensor.t);
      const std::complex<double> regge_ru        = ComplexSectionOverlapForTest(lts.hamp, u);
      const std::complex<double> tensor_ru       = ComplexSectionOverlapForTest(tensor.resonance, tensor.u);
      const std::complex<double> regge_forward_t = -regge.ExchangeKernel(lts.ss[1][3], lts.t1, test.pomeron_pdg) *
                                                   regge.ExchangeKernel(lts.ss[2][4], lts.t2, test.pomeron_pdg);
      const std::complex<double> regge_forward_u = -regge.ExchangeKernel(lts.ss[1][4], lts.t1, test.pomeron_pdg) *
                                                   regge.ExchangeKernel(lts.ss[2][3], lts.t2, test.pomeron_pdg);
      const std::complex<double> absolute_t =
          ComplexSectionOverlapForTest(t, tensor.t) / (regge_forward_t * std::conj(tensor.forward_t));
      const std::complex<double> absolute_u =
          ComplexSectionOverlapForTest(u, tensor.u) / (regge_forward_u * std::conj(tensor.forward_u));
      const std::complex<double> tu_ratio = regge_tu / tensor_tu;
      const std::complex<double> rt_ratio = regge_rt / tensor_rt;
      const std::complex<double> ru_ratio = regge_ru / tensor_ru;
      CAPTURE(point_index, points[point_index].t1, points[point_index].t2, regge_tu, tensor_tu, regge_rt, tensor_rt,
              regge_ru, tensor_ru, absolute_t, absolute_u, tu_ratio, rt_ratio, ru_ratio, std::arg(absolute_t),
              std::arg(absolute_u), std::arg(tu_ratio), std::arg(rt_ratio), std::arg(ru_ratio));
      CHECK(std::abs(absolute_t.imag()) <= 2.0e-10 * std::abs(absolute_t));
      CHECK(std::abs(absolute_u.imag()) <= 2.0e-10 * std::abs(absolute_u));
      CHECK(std::abs(tu_ratio.imag()) <= 2.0e-10 * std::abs(tu_ratio));
      CHECK(std::abs(rt_ratio.imag()) <= 2.0e-10 * std::abs(rt_ratio));
      CHECK(std::abs(ru_ratio.imag()) <= 2.0e-10 * std::abs(ru_ratio));
      CHECK(tu_ratio.real() > 0.0);
      CHECK(rt_ratio.real() > 0.0);
      CHECK(ru_ratio.real() > 0.0);
      CHECK(absolute_t.real() > 0.0);
      CHECK(absolute_u.real() > 0.0);
    }
  }
}

// Compare every generic Regge continuum against the covariant TP reference
TEST_CASE("MP XP and GP track the Tensor Pomeron pion continuum shape",
          "[gra::MRegge][MTensorPomeron][continuum][physics][regression]") {
  constexpr std::array<double, 6> masses         = {0.32, 0.40, 0.55, 0.80, 1.20, 1.60};
  constexpr std::size_t           reference_mass = 4;
  constexpr std::array<double, 5> cos_theta      = {-0.9061798459, -0.5384693101, 0.0, 0.5384693101, 0.9061798459};
  constexpr std::array<double, 5> angular_weight = {0.2369268851, 0.4786286705, 0.5688888889, 0.4786286705,
                                                    0.2369268851};
  constexpr std::array<double, 2> proton_pt      = {0.12, 0.28};
  const auto                      norm2          = [](const auto &amplitude) {
    double out = 0.0;
    for (const auto &value : amplitude) { out += std::norm(value); }
    return out;
  };

  const auto tune = WriteModifiedPhotoVMTune(
      "tensor_pomeron_pion_shape",
      [](auto &card) {
        for (auto &[key, pairs] : card.at("PARAM_TENSORPOM").at("PARAM_CON").at("TP").items()) {
          (void)key;
          pairs = {{995, 995}};
        }
      },
      [](auto &card) {
        for (auto &[exchange, pairs] : card.items()) {
          (void)exchange;
          for (auto &[pair, block] : pairs.items()) {
            (void)pair;
            block["FF_transfer"] = {{"type", "none"}};
            block["FF_offshell"] = {{"type", "none"}};
          }
        }
      },
      [](auto &card) {
        for (auto &[exchange, pairs] : card.items()) {
          (void)exchange;
          for (auto &[pair, block] : pairs.items()) {
            (void)pair;
            block["reggeize"]["active"] = false;
            block["pveto"]["active"] = false;
            block["FF_transfer"] = {{"type", "none"}};
            block["FF_offshell"] = {{"type", "none"}};
          }
        }
      });
  const auto         soft_model = gra::MModelTune::Load(tune.second);
  ToyHelicityProcess tp_process;
  ConfigureToyProductionProcess(tp_process, "TP", "CON", "pi+ pi-");
  tp_process.SetTuneForTest(tune.first);
  REQUIRE_NOTHROW(tp_process.InitializeProcessAmplitude());
  const auto tp_metadata = tp_process.state.lts.hamp.metadata;

  std::array<double, masses.size()>                               tp_profile{};
  std::array<std::vector<std::complex<double>>, cos_theta.size()> tp_reference;
  std::array<double, cos_theta.size()>                            tp_reference_norm{};
  for (const auto &mass_index : indices(masses)) {
    for (const auto &pt_index : indices(proton_pt)) {
      for (const auto &angle : indices(cos_theta)) {
        auto lts    = ScalarPolePhasePointForTest(proton_pt[pt_index], std::acos(cos_theta[angle]), 0.38,
                                                  masses[mass_index], 6500.0);
        lts.process = tp_process.state.lts.process;
        lts.hamp.Configure(tp_metadata);
        gra::MTensorPomeron tensor(lts, soft_model,
                                   gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
        tensor.ME4(lts, gra::TensorContinuumMode::TensorPomeron);
        const double intensity = norm2(lts.hamp);
        REQUIRE(std::isfinite(intensity));
        REQUIRE(intensity > 0.0);
        tp_profile[mass_index] += angular_weight[angle] * intensity;
        if (mass_index == reference_mass && pt_index == 0) {
          tp_reference[angle].assign(lts.hamp.begin(), lts.hamp.end());
          tp_reference_norm[angle] = intensity;
        }
      }
    }
  }
  REQUIRE(tp_profile[reference_mass] > 0.0);

  struct ContinuumModel {
    const char               *family;
    gra::ReggeProductionModel model;
    int                       pomeron;
  };
  constexpr std::array<ContinuumModel, 3> families = {{
      {"MP", gra::ReggeProductionModel::MP, 995},
      {"XP", gra::ReggeProductionModel::XP, 995},
      {"GP", gra::ReggeProductionModel::GP, 990},
  }};
  for (const auto &test : families) {
    const auto         name  = test.family;
    const auto         model = test.model;
    ToyHelicityProcess process;
    ConfigureToyProductionProcess(process, name, "CON", "pi+ pi-");
    process.SetTuneForTest(tune.first);
    REQUIRE_NOTHROW(process.InitializeProcessAmplitude());
    REQUIRE(process.state.lts.process.SPINGEN);
    REQUIRE(process.state.lts.process.FORWARD_VERTEX == gra::ForwardVertexMode::HelicityResidue);
    REQUIRE(process.state.lts.process.FORWARD_NOFLIP);
    const auto metadata = process.state.lts.hamp.metadata;
    REQUIRE(metadata.spin_basis == gra::ScreeningSpinBasis::ProtonHelicity);
    REQUIRE(metadata.spin_rows == 4);
    REQUIRE(metadata.forward_noflip);

    const auto channel =
        std::find(process.state.lts.process.CONT_PRODUCTION.begin(), process.state.lts.process.CONT_PRODUCTION.end(),
                  std::vector<int>{test.pomeron, test.pomeron});
    REQUIRE(channel != process.state.lts.process.CONT_PRODUCTION.end());
    const std::size_t channel_index =
        static_cast<std::size_t>(std::distance(process.state.lts.process.CONT_PRODUCTION.begin(), channel));
    const auto   selected_production = process.state.lts.process.CONT_PRODUCTION.at(channel_index);
    const auto   selected_tree       = process.state.lts.process.CONT_PRODUCTIONTREE.at(channel_index);
    const double selected_tu_sign    = process.state.lts.process.CONT_TU_SIGN.at(channel_index);

    const auto evaluate = [&](double mass, double pt, double theta) {
      auto lts    = ScalarPolePhasePointForTest(pt, theta, 0.38, mass, 6500.0);
      lts.process = process.state.lts.process;
      lts.hamp.Configure(metadata);
      lts.process.CONT_PRODUCTION     = {selected_production};
      lts.process.CONT_PRODUCTIONTREE = {selected_tree};
      if (model == gra::ReggeProductionModel::GP) {
        const auto selected      = lts.process.CONTINUUM_GP.at(channel_index);
        lts.process.CONTINUUM_GP = {selected};
        lts.process.CONTINUUM_POLE.clear();
      } else {
        const auto selected        = lts.process.CONTINUUM_POLE.at(channel_index);
        lts.process.CONTINUUM_POLE = {selected};
        lts.process.CONTINUUM_GP.clear();
      }
      lts.process.CONT_TU_SIGN = {selected_tu_sign};
      gra::MRegge regge(
          lts, process.state.model_tune,
          gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::ContinuumTwoBody, std::string(name) + "_tensor_shape"));
      TestReggeCon(regge, lts, model);
      REQUIRE(lts.hamp.size() == metadata.spin_rows);
      return std::vector<std::complex<double>>(lts.hamp.begin(), lts.hamp.end());
    };

    std::array<double, masses.size()>                               profile{};
    std::array<std::vector<std::complex<double>>, cos_theta.size()> reference;
    for (const auto &mass_index : indices(masses)) {
      for (const auto &pt_index : indices(proton_pt)) {
        for (const auto &angle : indices(cos_theta)) {
          const auto   amplitude = evaluate(masses[mass_index], proton_pt[pt_index], std::acos(cos_theta[angle]));
          const double intensity = norm2(amplitude);
          REQUIRE(std::isfinite(intensity));
          REQUIRE(intensity > 0.0);
          profile[mass_index] += angular_weight[angle] * intensity;
          if (mass_index == reference_mass && pt_index == 0) { reference[angle] = amplitude; }
        }
      }
    }
    REQUIRE(profile[reference_mass] > 0.0);
    for (const auto &mass_index : indices(masses)) {
      const double ratio =
          (profile[mass_index] / profile[reference_mass]) / (tp_profile[mass_index] / tp_profile[reference_mass]);
      CAPTURE(name, masses[mass_index], profile[mass_index], tp_profile[mass_index], ratio);
      CHECK(ratio == Approx(1.0).margin(model == gra::ReggeProductionModel::GP ? 0.30 : 0.25));
    }

    constexpr std::size_t      central_angle = 2;
    const std::complex<double> center =
        ComplexSectionOverlapForTest(reference[central_angle], tp_reference[central_angle]) /
        tp_reference_norm[central_angle];
    REQUIRE(std::abs(center) > 0.0);
    for (const auto &angle : indices(cos_theta)) {
      const std::complex<double> amplitude_shape =
          (ComplexSectionOverlapForTest(reference[angle], tp_reference[angle]) / tp_reference_norm[angle]) / center;
      CAPTURE(name, cos_theta[angle], amplitude_shape);
      CHECK(amplitude_shape.real() == Approx(1.0).margin(0.08));
      CHECK(std::abs(amplitude_shape.imag()) <= 0.08 * std::abs(amplitude_shape));
    }
  }
}

// Test initializing Durham model amplitudes
//

// Check the same complex rho production vertex through LS and helicity steering
TEST_CASE("MP and XP mixed rho production preserves LS leg exchange phases",
          "[gra::MRegge][MP][XP][resonance][helicity][normalization][beam-exchange]") {
  ModelParamRestoreGuard restore_model;
  gra::MODELPARAM = "TUNE0";
  gra::LORENTZSCALAR lts;
  lts.PDG = LoadedPDGTable();
  double symmetry = 1.0;
  gra::MProcessSetup setup(lts, gra::MModelTune::Load(modelfile), "MP", "RES",
                          gra::ReggeProcessInfo(gra::ReggeProductionModel::MP, gra::MReggeMode::Resonance),
                          false, false, 0, false, false, false, symmetry, {}, {}, {});
  gra::PARAM_RES res;
  res.p = lts.PDG.FindByPDG(113);
  const auto photon = lts.PDG.FindByPDG(22);
  const auto pomeron = lts.PDG.FindByPDG(995);
  for (const auto model : {gra::ReggeProductionModel::MP, gra::ReggeProductionModel::XP}) {
    gra::RES_PRODUCTION_CHANNEL channel;
    channel.exchange = {22, 995};
    channel.basis = gra::ReggeVertexBasis::LS;
    channel.Lambda = 1.0;
    const std::complex<double> coupling(0.7, -0.4);
    channel.g_ls.Set(0, 2, 0.0);
    channel.g_ls.Set(2, 2, 0.0);
    channel.g_ls.Set(2, 4, coupling);
    channel.g_ls.Set(2, 6, 0.0);
    channel.g_ls.Set(4, 6, 0.0);
    const auto direct = (model == gra::ReggeProductionModel::MP ? gra::mpom::PrepareResonance : gra::xpom::PrepareResonance)(setup, res, channel, {photon, pomeron});
    const auto reversed = (model == gra::ReggeProductionModel::MP ? gra::mpom::PrepareResonance : gra::xpom::PrepareResonance)(setup, res, channel, {pomeron, photon});
    REQUIRE(direct.terms.size() == 1);
    REQUIRE(reversed.terms.size() == 1);
    RequireComplexNear(direct.terms.front().coefficient, coupling, 1.0e-12);
    RequireComplexNear(reversed.terms.front().coefficient, -coupling, 1.0e-12);

    const auto h       = gra::spin::PoleLSReduced(direct, 1.0);
    channel.basis = gra::ReggeVertexBasis::Helicity;
    channel.P_symmetry = false;
    for (std::size_t i = 0; i < h.size_row(); ++i) {
      for (std::size_t j = 0; j < h.size_col(); ++j) {
        if (std::abs(h[i][j]) < 1.0e-12) { continue; }
        channel.helicity.push_back({static_cast<double>(i) - 1.0, static_cast<double>(j) - 2.0});
        channel.g_helicity.push_back(h[i][j]);
      }
    }
    const auto direct_h = (model == gra::ReggeProductionModel::MP ? gra::mpom::PrepareResonance : gra::xpom::PrepareResonance)(setup, res, channel, {photon, pomeron});
    const auto reversed_h = (model == gra::ReggeProductionModel::MP ? gra::mpom::PrepareResonance : gra::xpom::PrepareResonance)(setup, res, channel, {pomeron, photon});
    RequireMatrixNear(gra::spin::PoleLSReduced(direct_h, 1.0), h, 1.0e-12);
    RequireMatrixNear(gra::spin::PoleLSReduced(reversed_h, 1.0), gra::spin::PoleLSReduced(reversed, 1.0), 1.0e-12);
    for (const auto basis : {gra::ReggeVertexBasis::AutoMinL, gra::ReggeVertexBasis::AutoMinS,
                             gra::ReggeVertexBasis::AutoEqualLS, gra::ReggeVertexBasis::AutoEqualHelicity}) {
      channel.basis = basis;
      channel.P_symmetry = true;
      channel.g = coupling;
      const auto prepare = model == gra::ReggeProductionModel::MP ? gra::mpom::PrepareResonance : gra::xpom::PrepareResonance;
      const auto generated = prepare(setup, res, channel, {photon, pomeron});
      const auto reverse = prepare(setup, res, channel, {pomeron, photon});
      REQUIRE(generated.terms.size() == reverse.terms.size());
      auto explicit_ls = channel;
      explicit_ls.basis = gra::ReggeVertexBasis::LS;
      for (auto& term : explicit_ls.g_ls) { term.coefficient = 0.0; }
      for (const auto i : indices(generated.terms)) {
        const auto& term = generated.terms[i];
        const int exponent = term.l + (photon.spinX2 + pomeron.spinX2 - term.two_s) / 2;
        RequireComplexNear(reverse.terms[i].coefficient, (exponent % 2 == 0 ? 1.0 : -1.0) * term.coefficient, 1e-12);
        explicit_ls.g_ls.Set(term.l, term.two_s, term.coefficient);
      }
      const auto reference = prepare(setup, res, explicit_ls, {pomeron, photon});
      for (const double momentum : {0.7, 2.1}) {
        RequireMatrixNear(gra::spin::PoleLSReduced(reverse, momentum), gra::spin::PoleLSReduced(reference, momentum), 1e-12);
      }
    }
  }
}

// Read photon channels through the real card reader without a SOFT photon trajectory
TEST_CASE("Regge continuum cards accept photon and mixed exchange pairs",
          "[gra::MRegge][continuum][photon][params]") {
  const auto tune = WriteModifiedPhotoVMTune("regge_photon_channels", [](auto &card) {
    for (const std::string model : {"MP", "XP", "GP"}) {
      const int pomeron = model == "GP" ? 990 : 995;
      card["PARAM_REGGE"]["PARAM_CON"][model]["[211,-211]"] = {{22, 22}, {22, pomeron}};
    }
    card["PARAM_REGGE"]["PARAM_CON"]["MULTI"]["secondary_exchanges"] = false;
  }, {}, SetToyPhotonContinuum);
  const auto pdg = LoadedPDGTable();
  const std::vector<int> pions{211, -211};
  const auto param = gra::regge::ReadParam(pions, pdg, *gra::MModelTune::Load(tune.second));
  REQUIRE_THROWS_AS(gra::regge::TrajectoryIndex(param, 22), std::invalid_argument);
  CHECK_FALSE(gra::regge::IsSecondaryReggeonTrajectory(param, 22));
  CHECK_FALSE(gra::regge::CheckVertex(param, pdg, 22, {211, -211}).applies);
  for (const auto model : {gra::ReggeProductionModel::MP, gra::ReggeProductionModel::XP, gra::ReggeProductionModel::GP}) {
    const auto &pair = gra::regge::Pair(param, pions, model);
    REQUIRE(pair.channels.size() == 3);
    CHECK(pair.channels[0].first == 22);
    CHECK(pair.channels[0].second == 22);
    CHECK(pair.channels[1].first == pair.channels[2].second);
    CHECK(pair.channels[1].second == pair.channels[2].first);
    CHECK_FALSE(gra::regge::ContinuumLadderEntriesConnect(param, {&pair, &pair}, model));
    for (const auto &vertex : pair.channels) { CHECK_FALSE(gra::regge::VertexAllowed(param, vertex)); }
    auto hadronic = pair;
    const int pomeron = model == gra::ReggeProductionModel::GP ? 990 : 995;
    hadronic.channels.push_back({pomeron, pomeron});
    CHECK(gra::regge::ContinuumLadderEntriesConnect(param, {&hadronic, &hadronic}, model));
    auto secondary = param;
    secondary.con.multiregge_secondary_exchanges = true;
    CHECK_FALSE(gra::regge::ContinuumLadderEntriesConnect(secondary, {&pair, &pair}, model));
  }
}

// Read photon transfer forms without inventing a Regge trajectory for the photon
TEST_CASE("Regge continuum photon vertex forms use the production card",
          "[gra::MRegge][continuum][photon][form-factor]") {
  const auto tune = WriteModifiedContinuumTune("regge_photon_form", "MP", [](auto &card) {
    SetToyPhotonContinuum(card);
    card["22"]["[211,211]"]["FF_transfer"] = {
        {"type", "power"}, {"norm", "zero"}, {"Lambda2", 0.5}, {"n", 1.0}};
  });
  const std::string general_path = tune + "/GENERAL.json";
  auto general = nlohmann::json::parse(gra::aux::GetInputData(general_path));
  general["PARAM_REGGE"]["PARAM_CON"]["MP"]["[211,-211]"] = {{22, 22}};
  std::ofstream output(general_path);
  REQUIRE(output.good());
  output << general.dump(2);
  output.close();
  const std::vector<int> pions{211, -211};
  const auto param = gra::regge::ReadParam(pions, LoadedPDGTable(), *gra::MModelTune::Load(general_path));
  const auto &pair = gra::regge::Pair(param, pions, gra::ReggeProductionModel::MP);
  REQUIRE(pair.channels.size() == 1);
  for (const auto &form : pair.channels.front().transfer) {
    CHECK(form.type == gra::regge::FFType::Power);
    REQUIRE(form.param.size() == 2);
    CHECK(form.param[0] == Approx(0.5));
    CHECK(form.param[1] == Approx(1.0));
  }
}

// Validate trajectory coverage when reggeization is enabled, before event sampling
TEST_CASE("Regge continuum input requires every enabled meson trajectory",
          "[gra::MRegge][continuum][params][reggeization]") {
  for (const bool reggeize : {false, true}) {
    const auto tune = WriteModifiedPhotoVMTune(reggeize ? "missing_active_meson_trajectory" : "missing_unused_meson_trajectory",
                                               [reggeize](auto &card) {
      auto &regge = card.at("PARAM_REGGE");
      auto &rows = regge.at("meson_trajectories");
      for (auto row = rows.begin(); row != rows.end(); ++row) {
        if ((*row)[0].template get<int>() == 211) {
          rows.erase(row);
          break;
        }
      }
    }, {}, [reggeize](auto &card) { SetContinuumField(card, "[211,211]", "reggeize", {{"active", reggeize}, {"freeze_scale2", 1.0}}); });
    const auto model = gra::MModelTune::Load(tune.second);
    if (reggeize) {
      REQUIRE_THROWS_AS(gra::regge::ReadParam({211, -211}, LoadedPDGTable(), *model), std::invalid_argument);
    } else {
      REQUIRE_NOTHROW(gra::regge::ReadParam({211, -211}, LoadedPDGTable(), *model));
    }
  }
}

// Require each trajectory to meet the spin of its declared physical pole
TEST_CASE("Regge meson trajectories validate the physical pole spin",
          "[gra::MRegge][continuum][params][reggeization]") {
  const auto tune = WriteModifiedDurhamTune("wrong_meson_pole_spin", [](auto &card) {
    for (auto &row : card.at("PARAM_REGGE").at("meson_trajectories")) {
      if (row[0].template get<int>() == 211) { row[1] = 1.0; }
    }
  });
  REQUIRE_THROWS_AS(gra::regge::ReadParam({211, -211}, LoadedPDGTable(), *gra::MModelTune::Load(tune.second)),
                    std::invalid_argument);
}

// Trace hidden spin states in separate sectors while retaining equal-spin interference
TEST_CASE("MRegge isotropic spin sums combine continuum and resonance sectors", "[gra::MRegge][spin][physics]") {
  const auto tune = gra::MModelTune::Load(modelfile);
  for (const auto model : {gra::ReggeProductionModel::MP, gra::ReggeProductionModel::XP, gra::ReggeProductionModel::GP}) {
    for (const bool spingen : {false, true}) {
      CAPTURE(model, spingen);
      auto base = MakeToyReggeLTSAsymmetric(0.3, 4.2, -4.6);
      base.process.MMAX = 2;
      base.process.SPINGEN = spingen;
      base.process.SPINDEC = false;
      base.process.TU_SIGN = "positive";
      if (model == gra::ReggeProductionModel::GP) {
        SetToyContinuumExchangePair(base, 990, 990);
        UseToyGPSubchannelHelicity(base);
      } else if (model == gra::ReggeProductionModel::XP) {
        PrepareToyXPContinuumOperators(base, {{1, 2, 1.0}});
      }
      auto first = model == gra::ReggeProductionModel::MP ? MakeToyMPResonance(false)
                   : model == gra::ReggeProductionModel::XP ? MakeToyCovariantXPonance() : MakeToyGPResonance(2);
      auto second = first;
      second.hel_decay.g_decay *= std::complex<double>(0.5, 0.25);
      // Evaluate the full public process for each physical contribution
      const auto evaluate = [&](const auto &resonances, gra::MReggeMode mode) {
        auto lts = base;
        lts.process.RESONANCES = resonances;
        gra::MRegge regge(lts, tune, gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));
        return regge.Amp2(lts, model, mode);
      };
      using Resonances = decltype(base.process.RESONANCES);
      const double continuum = evaluate(Resonances{}, gra::MReggeMode::ContinuumTwoBody);
      const double resonance = evaluate(Resonances{{"first", first}}, gra::MReggeMode::Resonance);
      const double coherent = std::norm(std::complex<double>(1.5, 0.25)) * resonance;
      CHECK(evaluate(Resonances{{"first", first}, {"second", second}}, gra::MReggeMode::ResonanceContinuumTwoBody)
            == Approx(continuum + coherent).epsilon(1e-10));
      if (model == gra::ReggeProductionModel::GP) {
        const auto tensor = MakeToyTensorGPResonance(2);
        const auto scalar = MakeToyScalarGPResonance(2);
        const double tensor_norm = evaluate(Resonances{{"tensor", tensor}}, gra::MReggeMode::Resonance);
        const double scalar_norm = evaluate(Resonances{{"scalar", scalar}}, gra::MReggeMode::Resonance);
        CHECK(evaluate(Resonances{{"first", first}, {"second", second}, {"tensor", tensor}, {"scalar", scalar}}, gra::MReggeMode::ResonanceContinuumTwoBody)
              == Approx(continuum + coherent + tensor_norm + scalar_norm).epsilon(1e-10));
      }
    }
  }
}

// Check the complex continued LS sum and restrict MP/XP matching to the integer pole
TEST_CASE("GP photo vertex follows the sliding LS amplitude", "[gra::MRegge][GP][photoproduction][physics]") {
  const auto  tune     = gra::MModelTune::Load(modelfile);
  const auto  pdg      = LoadedPDGTable();
  const auto  param    = gra::regge::ReadParam({211, -211}, pdg, *tune);
  const auto& mother   = pdg.FindByPDG(113);
  const auto& photon   = pdg.FindByPDG(22);
  const auto& exchange = pdg.FindByPDG(990);
  const auto& pole     = gra::regge::PoleRepresentative(param, pdg, 990);
  for (const bool reverse : {false, true}) {
    const std::vector<gra::MParticle> legs = reverse ? std::vector{exchange, photon} : std::vector{photon, exchange};
    gra::RES_PRODUCTION_CHANNEL       channel;
    channel.exchange = {22, 990};
    channel.basis    = gra::ReggeVertexBasis::LS;
    for (const auto& entry :
         gra::spin::CanonicalPoleOperators(mother, photon, pole, true, channel.C_symmetry, channel.P_symmetry)) {
      channel.g_ls.Set(entry.coupling.l, entry.coupling.two_s, 0.0);
    }
    channel.g_ls.Set(0, 2, 1.0);
    channel.g_ls.Set(2, 2, {0.5, 0.25});
    gra::RES_PRODUCTION gp;
    gp.tree.resize(2);
    gp.tree[0].p = legs[0];
    gp.tree[1].p = legs[1];
    gp.hel       = gra::gpom::PrepareResonance(mother, legs, channel, param, pdg, 2, 0.0);
    for (const bool derivative : {false, true}) {
      for (const double momentum : {0.6, 1.0, 1.7}) {
        for (const double alpha : {1.08, 1.5, 2.0}) {
          CAPTURE(reverse, derivative, momentum, alpha);
          const double         j1       = reverse ? alpha : 1.0;
          const double         j2       = reverse ? 1.0 : alpha;
          const double         m1       = reverse ? 0.0 : 1.0;
          const double         m2       = reverse ? 1.0 : 0.0;
          const double         mu       = m1 - m2;
          std::complex<double> expected = 0.0;
          for (const auto& term : gp.hel.alpha_ls) {
            const double S = 0.5 * term.two_s;
            const double raw =
                gra::spin::RawPoleLSNormalization(2, pole.spinX2, 2, term.l, static_cast<int>(term.two_s));
            const double radial  = derivative ? std::pow(momentum / channel.Lambda, term.l) : 1.0;
            const double orbital = std::sqrt((2.0 * term.l + 1.0) / 3.0) * gra::wigner::CG(term.l, S, 0.0, mu, 1.0, mu);
            const auto   cg =
                0.5 * (gra::wigner::CGRegge(j1, j2, m1, -m2, S, mu) + gra::wigner::CGRegge(j1, j2, -m1, m2, S, -mu));
            expected += term.coefficient * radial * raw * orbital * cg;
          }
          const auto actual = gra::gpom::PhotoCoupling(gp, alpha, momentum, derivative);
          CHECK(std::abs(actual - expected) < 1e-12);
        }
      }
    }
    auto frozen = param;
    frozen.gp_alpha_min = 0.25;
    auto event = MakeToyCoherentPhotonLTS();
    event.process.DERIVATIVE_FACTOR = false;
    for (const double alpha : {-3.25, 0.75}) {
      const double t = ReggeTransferForAlpha(frozen, 990, alpha);
      RequireComplexNear(gra::gpom::PhotoCoupling(event, frozen, gp, t),
                         gra::gpom::PhotoCoupling(gp, std::max(alpha, frozen.gp_alpha_min), 0.0, false), 2.0e-12);
    }
    gra::RES_PRODUCTION reference;
    reference.pole =
        gra::spin::PreparePoleLS(mother, photon, pole, {{0, 2, 1.0}, {2, 2, {0.5, 0.25}}}, 1.0, true, false, true);
    const double expected = gra::rspin::PhotoCoupling(reference);
    CHECK(std::abs(gra::gpom::PhotoCoupling(gp, pole.spinX2 / 2.0, 1.0, true)) == Approx(expected).epsilon(1e-12));
    CHECK(std::abs(gra::gpom::PhotoCoupling(gp, 1.08, 1.0, true)) < 0.5 * expected);
    const auto helicity = gra::spin::PoleLSHelicity(*reference.pole, 1.0);
    channel.basis       = gra::ReggeVertexBasis::Helicity;
    for (int target = -2; target <= 0; ++target) {
      channel.helicity.push_back({-1.0, static_cast<double>(target)});
      channel.g_helicity.push_back(helicity.T[0][target + 2]);
    }
    gp.hel = gra::gpom::PrepareResonance(mother, legs, channel, param, pdg, 2, 0.0);
    CHECK(std::abs(gra::gpom::PhotoCoupling(gp, 2.0, 1.0, true)) == Approx(expected).epsilon(1e-12));
  }
}

// A transfer-suppressed helicity row must never become a nonzero forward reference
TEST_CASE("GP photo vertex keeps genuine forward zeros", "[gra::MRegge][GP][photoproduction][physics]") {
  const auto pdg   = LoadedPDGTable();
  const auto tune  = gra::MModelTune::Load(modelfile);
  const auto param = gra::regge::ReadParam({211, -211}, pdg, *tune);
  for (const bool reverse : {false, true}) {
    gra::RES_PRODUCTION gp;
    gp.tree.resize(2);
    gp.tree[0].p = pdg.FindByPDG(reverse ? 990 : 22);
    gp.tree[1].p = pdg.FindByPDG(reverse ? 22 : 990);
    gra::RES_PRODUCTION_CHANNEL channel;
    channel.exchange = {22, 990};
    channel.basis    = gra::ReggeVertexBasis::Helicity;
    channel.helicity = {{-1, 0}, {-1, -1}};
    for (const double nonflip : {0.0, 1.0}) {
      channel.g_helicity = {{nonflip, 0.0}, {2.0, 1.0}};
      gp.hel = gra::gpom::PrepareResonance(pdg.FindByPDG(113), {gp.tree[0].p, gp.tree[1].p}, channel, param, pdg, 2, 0.0);
      CHECK(std::abs(gra::gpom::PhotoCoupling(gp, 1.08, 1.0, true) - nonflip) < 1e-12);
    }
  }
}
