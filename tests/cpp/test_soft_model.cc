// Immutable SOFT model and arbitrary-channel Good Walker tests
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <string>
#include <utility>
#include <vector>

// Own
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Particle/MPDG.h"
#include "Graniitti/Regge/MSoftModel.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MException.h"

// Libraries
#include "json.hpp"
#include <catch.hpp>

using gra::aux::indices;

namespace {

using json = nlohmann::json;

// Load the checked-in TUNE0 document through the standard input reader
json LoadTuneDocument() {
  const std::string path =
      gra::aux::ResolveProjectPath("modeldata/TUNE0/GENERAL.json");
  return json::parse(gra::aux::GetInputData(path));
}

// Construct one immutable model snapshot without writing a temporary card
gra::SoftModelPtr LoadModel(const json &document, const std::string &label) {
  return gra::SoftModel::LoadFromJson(label, document.dump());
}

// Generate a complete four-channel SOFT fixture from the one-channel card
json FourChannelDocument() {
  json document = LoadTuneDocument();
  json model = document.at("PARAM_SOFT").at("MODEL").at("single");
  constexpr std::size_t channels = 4;

  model["GW"]["theta"] = {0.20, -0.10, 0.30, 0.15, -0.25, 0.18};
  model["GW"]["a_c"] = {1.0, 2.0, 3.0};

  for (auto item : model["FF"].items()) {
    auto &bank = item.value();
    const json row = bank.at("param").at(0);
    bank["param"] = json::array();
    for (std::size_t i = 0; i < channels; ++i) {
      json channel_row = row;
      channel_row.at(0) =
          row.at(0).get<double>() * (1.0 + 0.05 * static_cast<double>(i));
      bank["param"].push_back(std::move(channel_row));
    }
  }

  for (auto item : model["EXCHANGE"].items()) {
    auto &exchange = item.value();
    const double base = exchange.at("g").at(0).at(0).get<double>();
    exchange["g"] = json::array();
    for (std::size_t i = 0; i < channels; ++i) {
      json row = json::array();
      for (std::size_t j = 0; j < channels; ++j) {
        row.push_back(i == j ? base * (1.0 + 0.1 * i) : 0.0);
      }
      exchange["g"].push_back(std::move(row));
    }
  }
  model["EXCHANGE"]["P"]["g"] = {{4.0, 0.5, 0.2, 0.1},
                                 {0.5, 3.0, 0.3, 0.2},
                                 {0.2, 0.3, 2.0, 0.4},
                                 {0.1, 0.2, 0.4, 1.5}};

  document["PARAM_SOFT"]["active_model"] = "generated_four";
  document["PARAM_SOFT"]["MODEL"]["generated_four"] = std::move(model);
  return document;
}

// Multiply two matrices using an explicit reference loop
gra::MMatrix<double> Multiply(const gra::MMatrix<double> &first,
                              const gra::MMatrix<double> &second) {
  REQUIRE(first.size_col() == second.size_row());
  gra::MMatrix<double> product(first.size_row(), second.size_col(), 0.0);
  for (std::size_t i = 0; i < first.size_row(); ++i) {
    for (std::size_t j = 0; j < second.size_col(); ++j) {
      for (std::size_t k = 0; k < first.size_col(); ++k) {
        product(i, j) += first(i, k) * second(k, j);
      }
    }
  }
  return product;
}

// Check two real matrices entry by entry
void CheckMatrix(const gra::MMatrix<double> &value,
                 const gra::MMatrix<double> &reference,
                 const double tolerance = 1.0e-11) {
  REQUIRE(value.size_row() == reference.size_row());
  REQUIRE(value.size_col() == reference.size_col());
  for (std::size_t i = 0; i < value.size_row(); ++i) {
    for (std::size_t j = 0; j < value.size_col(); ++j) {
      CAPTURE(i, j);
      CHECK(value(i, j) == Approx(reference(i, j)).margin(tolerance));
    }
  }
}

// Check one matrix is a symmetric orthogonal projector
void CheckProjector(const gra::MMatrix<double> &projector) {
  REQUIRE(projector.size_row() == projector.size_col());
  CheckMatrix(projector, projector.Transpose());
  CheckMatrix(Multiply(projector, projector), projector);
}

// Build an identity matrix without using the projector implementation
gra::MMatrix<double> Identity(const std::size_t size) {
  gra::MMatrix<double> identity(size, size, 0.0);
  for (std::size_t i = 0; i < size; ++i) {
    identity(i, i) = 1.0;
  }
  return identity;
}

// Reconstruct the triple Pomeron target matrix from its principal root
void CheckTriplePomeronRoot(const gra::SoftModel &model,
                            const gra::SoftExchangeId exchange) {
  const auto root = model.TriplePomeronCouplingRoot(exchange);
  const auto reconstructed = Multiply(root.Transpose(), root);
  const auto target = model.Exchange(exchange).CouplingMatrix() *
                      model.EffectiveTriplePomeronCoupling(exchange);
  CheckMatrix(root, root.Transpose(), 2.0e-10);
  CheckMatrix(reconstructed, target, 2.0e-9);
}

// Relabel one square channel matrix with a new to old eigenstate permutation
template <typename T>
gra::MMatrix<T> PermuteChannelMatrix(const gra::MMatrix<T> &matrix,
                                     const std::vector<std::size_t> &order) {
  REQUIRE(matrix.size_row() == order.size());
  REQUIRE(matrix.size_col() == order.size());
  gra::MMatrix<T> permuted(order.size(), order.size(), T{});
  for (std::size_t i = 0; i < order.size(); ++i) {
    REQUIRE(order[i] < order.size());
    for (std::size_t j = 0; j < order.size(); ++j) {
      permuted(i, j) = matrix(order[i], order[j]);
    }
  }
  return permuted;
}

// Relabel one channel vector with a new to old eigenstate permutation
template <typename T>
std::vector<T> PermuteChannelVector(const std::vector<T> &source,
                                    const std::vector<std::size_t> &order) {
  REQUIRE(source.size() == order.size());
  std::vector<T> permuted(order.size());
  for (std::size_t i = 0; i < order.size(); ++i) {
    REQUIRE(order[i] < order.size());
    permuted[i] = source[order[i]];
  }
  return permuted;
}

} // namespace

TEST_CASE("Generic inner products support real and complex vectors",
          "[gra::InnerProduct]") {
  CHECK(gra::InnerProduct(std::vector<double>{1.0, 2.0, 3.0},
                          std::vector<double>{4.0, 5.0, 6.0}) == Approx(32.0));

  const std::vector<std::complex<double>> first{{1.0, 2.0}, {3.0, -1.0}};
  const std::vector<std::complex<double>> second{{2.0, -1.0}, {-1.0, 4.0}};
  CHECK(gra::InnerProduct(first, second) == std::complex<double>{-7.0, 6.0});

  const double inverse_sqrt_two = 1.0 / std::sqrt(2.0);
  const std::vector<std::complex<double>> state{{inverse_sqrt_two, 0.0},
                                                {0.0, inverse_sqrt_two}};
  const auto projector = gra::RankOneProjector(state);
  CHECK(projector(0, 0).real() == Approx(0.5));
  CHECK(projector(0, 1).imag() == Approx(-0.5));
  CHECK(projector(1, 0).imag() == Approx(0.5));
  CHECK(projector(1, 1).real() == Approx(0.5));

  CHECK_THROWS_AS(
      gra::InnerProduct(std::vector<float>{1.0F}, std::vector<float>{}),
      std::invalid_argument);
}

TEST_CASE("Forward excitation selection is independent of eikonal exchanges",
          "[gra::SoftModel][mapping][validation]") {
  json document = LoadTuneDocument();
  auto &soft = document["PARAM_SOFT"];
  const std::string active_model = soft["active_model"].get<std::string>();
  auto &model = soft["MODEL"][active_model];
  soft["EXCHANGE_DEF"]["P_aux"] = soft["EXCHANGE_DEF"]["P"];
  model["EXCHANGE"]["P_aux"] = model["EXCHANGE"]["P"];
  model["EIKONAL"]["screening_exchanges"] = {"P", "P_aux"};

  const auto primary = LoadModel(document, "two screening Pomerons fixture");
  REQUIRE(primary->Eikonal().ScreeningExchanges().size() == 2);
  CHECK(primary->ForwardExcitationExchange() == primary->ExchangeId("P"));

  SECTION("wildcard selects every eligible excitation Pomeron") {
    json wildcard = document;
    wildcard["PARAM_SOFT"]["MODEL"][active_model]["EIKONAL"]
            ["excitation_exchanges"] = {"*"};
    wildcard["PARAM_SOFT"]["MODEL"][active_model]["EIKONAL"]
            ["screening_exchanges"] = {"P"};
    wildcard["PARAM_SOFT"]["MODEL"][active_model]["EXCHANGE"]["P_aux"]["on"] =
        false;
    const auto selected =
        LoadModel(wildcard, "wildcard forward Pomeron fixture");
    CHECK(selected->ForwardExcitationExchange() == selected->ExchangeId("P"));
  }

  document["PARAM_SOFT"]["MODEL"][active_model]["EIKONAL"]
          ["excitation_exchanges"] = {"P_aux"};
  const auto auxiliary =
      LoadModel(document, "auxiliary forward Pomeron fixture");
  CHECK(auxiliary->ForwardExcitationExchange() ==
        auxiliary->ExchangeId("P_aux"));

  SECTION("unselected disabled Pomeron cannot supply an excitation factor") {
    json invalid = document;
    invalid["PARAM_SOFT"]["MODEL"][active_model]["EIKONAL"]
           ["excitation_exchanges"] = {"P"};
    invalid["PARAM_SOFT"]["MODEL"][active_model]["EIKONAL"]
           ["screening_exchanges"] = {"P"};
    invalid["PARAM_SOFT"]["MODEL"][active_model]["EXCHANGE"]["P_aux"]["on"] =
        false;
    const auto disabled =
        LoadModel(invalid, "unselected disabled forward Pomeron fixture");
    CHECK_THROWS(disabled->ForwardExcitationFactor(
        disabled->ExchangeId("P_aux"), -0.2, 25.0));
  }

  SECTION("missing selector") {
    json invalid = document;
    invalid["PARAM_SOFT"]["MODEL"][active_model]["EIKONAL"].erase(
        "excitation_exchanges");
    CHECK_THROWS(LoadModel(invalid, "missing forward selector fixture"));
  }
  SECTION("unknown exchange") {
    json invalid = document;
    invalid["PARAM_SOFT"]["MODEL"][active_model]["EIKONAL"]
           ["excitation_exchanges"] = {"missing"};
    CHECK_THROWS(LoadModel(invalid, "unknown forward exchange fixture"));
  }
  SECTION("disabled exchange") {
    json invalid = document;
    invalid["PARAM_SOFT"]["MODEL"][active_model]["EXCHANGE"]["P_aux"]["on"] =
        false;
    CHECK_THROWS(LoadModel(invalid, "disabled forward exchange fixture"));
  }
  SECTION("non-Pomeron exchange") {
    json invalid = document;
    invalid["PARAM_SOFT"]["MODEL"][active_model]["EIKONAL"]
           ["excitation_exchanges"] = {"R_f2"};
    CHECK_THROWS(LoadModel(invalid, "non-Pomeron forward fixture"));
  }
  SECTION("multiple excitation exchanges") {
    json invalid = document;
    invalid["PARAM_SOFT"]["MODEL"][active_model]["EIKONAL"]
           ["excitation_exchanges"] = {"P", "P_aux"};
    CHECK_THROWS(LoadModel(invalid, "multiple forward selector fixture"));
  }
}

TEST_CASE("SOFT exchange enable flags are required booleans",
          "[gra::SoftModel][mapping][validation]") {
  const json valid = LoadTuneDocument();
  const std::string active_model =
      valid.at("PARAM_SOFT").at("active_model").get<std::string>();

  SECTION("missing enable flag") {
    json invalid = valid;
    invalid["PARAM_SOFT"]["MODEL"][active_model]["EXCHANGE"]["P"].erase("on");
    CHECK_THROWS(LoadModel(invalid, "missing exchange on fixture"));
  }

  SECTION("nonboolean enable flag") {
    json invalid = valid;
    invalid["PARAM_SOFT"]["MODEL"][active_model]["EXCHANGE"]["P"]["on"] = 1;
    CHECK_THROWS(LoadModel(invalid, "nonboolean exchange on fixture"));
  }
}

TEST_CASE("SOFT exchange definitions reject fractional discrete values",
          "[gra::SoftModel][mapping][validation]") {
  const json valid = LoadTuneDocument();
  const std::string active_model =
      valid.at("PARAM_SOFT").at("active_model").get<std::string>();

  SECTION("fractional crossing parity") {
    json invalid = valid;
    invalid["PARAM_SOFT"]["EXCHANGE_DEF"]["P"]["crossing"] = 0.5;
    CHECK_THROWS(LoadModel(invalid, "fractional crossing fixture"));
  }
  SECTION("fractional residue sign") {
    json invalid = valid;
    invalid["PARAM_SOFT"]["MODEL"][active_model]["EXCHANGE"]["P"]["sign"] = 0.5;
    CHECK_THROWS(LoadModel(invalid, "fractional residue sign fixture"));
  }
}

// Check that failed generated transfers use the amplitude failure bookkeeping
TEST_CASE("SOFT generated kinematics use amplitude failures",
          "[gra::SoftModel][validation][failure]") {
  auto document = LoadTuneDocument();
  const auto model = LoadModel(document, "generated SOFT kinematics");
  const auto pomeron = model->ExchangeId("P");
  const double nan = std::numeric_limits<double>::quiet_NaN();
  CHECK_THROWS_AS(model->Alpha(pomeron, nan), gra::AmplitudeFailure);
  CHECK_THROWS_AS(model->FormFactor(pomeron, nan, 0), gra::AmplitudeFailure);
  CHECK_THROWS_AS(model->ForwardExcitationFactor(pomeron, -0.1, nan), gra::AmplitudeFailure);

  document["PARAM_SOFT"]["active_model"] = "single";
  auto &bank = document["PARAM_SOFT"]["MODEL"]["single"]["FF"]["P"];
  bank = {{"type", "GKERNEL"}, {"param", {{1.0, 1.5, 0.5, 0.1}}}};
  const auto kernel = LoadModel(document, "SOFT kernel domain");
  CHECK_THROWS_AS(kernel->FormFactor(kernel->ExchangeId("P"), 0.2, 0), gra::AmplitudeFailure);
  bank["param"][0][0] = -1.0;
  CHECK_THROWS_AS(LoadModel(document, "invalid SOFT kernel input"), std::invalid_argument);
}

TEST_CASE("SOFT helicity transitions accept scalar or symmetric matrices",
          "[gra::SoftModel][helicity][N2][validation]") {
  json valid = LoadTuneDocument();
  valid["PARAM_SOFT"]["active_model"] = "double";
  const std::string active_model = "double";

  SECTION("scalar controls broadcast over multiple channels") {
    json scalar = valid;
    scalar["PARAM_SOFT"]["MODEL"][active_model]["EIKONAL"]["helicity"] = true;
    auto &helicity = scalar["PARAM_SOFT"]["MODEL"][active_model]["EXCHANGE"]
                           ["P"]["helicity"];
    helicity["kappa"] = -0.21;
    helicity["B_kappa"] = 0.37;
    const auto model = LoadModel(scalar, "scalar helicity fixture");
    const auto &exchange = model->Exchange(model->ExchangeId("P"));
    for (std::size_t i = 0; i < 2; ++i) {
      for (std::size_t j = 0; j < 2; ++j) {
        CHECK(exchange.HelicityFlipCoupling(i, j) == Approx(-0.21));
        CHECK(exchange.HelicityFlipSlope(i, j) == Approx(0.37));
      }
    }
  }

  SECTION("symmetric channel matrices retain each entry") {
    json matrix = valid;
    matrix["PARAM_SOFT"]["MODEL"][active_model]["EIKONAL"]["helicity"] = true;
    auto &helicity = matrix["PARAM_SOFT"]["MODEL"][active_model]["EXCHANGE"]
                           ["P"]["helicity"];
    helicity["kappa"] = {{0.11, -0.23}, {-0.23, 0.31}};
    helicity["B_kappa"] = {{0.41, 0.17}, {0.17, 0.53}};
    const auto model = LoadModel(matrix, "matrix helicity fixture");
    const auto &exchange = model->Exchange(model->ExchangeId("P"));
    CHECK(exchange.HelicityFlipCoupling(0, 0) == Approx(0.11));
    CHECK(exchange.HelicityFlipCoupling(0, 1) == Approx(-0.23));
    CHECK(exchange.HelicityFlipCoupling(1, 0) == Approx(-0.23));
    CHECK(exchange.HelicityFlipCoupling(1, 1) == Approx(0.31));
    CHECK(exchange.HelicityFlipSlope(0, 0) == Approx(0.41));
    CHECK(exchange.HelicityFlipSlope(0, 1) == Approx(0.17));
    CHECK(exchange.HelicityFlipSlope(1, 0) == Approx(0.17));
    CHECK(exchange.HelicityFlipSlope(1, 1) == Approx(0.53));

    const double t = -0.37;
    const auto nonflip = model->ResidueMatrix(model->ExchangeId("P"), t);
    const auto flip =
        model->HelicityFlipResidueMatrix(model->ExchangeId("P"), t);
    const std::array<std::array<double, 2>, 2> kappa = {
        {{{0.11, -0.23}}, {{-0.23, 0.31}}}};
    const std::array<std::array<double, 2>, 2> slope = {
        {{{0.41, 0.17}}, {{0.17, 0.53}}}};
    for (const auto &i : gra::aux::indices(kappa)) {
      for (const auto &j : gra::aux::indices(kappa[i])) {
        CHECK(flip[i][j] == Approx(nonflip[i][j] * kappa[i][j] *
                                   std::exp(0.5 * slope[i][j] * t)));
      }
    }
  }

  SECTION("disabled helicity gives a zero flip residue") {
    const auto model = LoadModel(valid, "disabled helicity fixture");
    const auto flip =
        model->HelicityFlipResidueMatrix(model->ExchangeId("P"), -0.31);
    CHECK(flip.FrobNorm2() == Approx(0.0));
  }

  SECTION("asymmetric matrix is rejected") {
    json asymmetric = valid;
    asymmetric["PARAM_SOFT"]["MODEL"][active_model]["EIKONAL"]["helicity"] =
        true;
    asymmetric["PARAM_SOFT"]["MODEL"][active_model]["EXCHANGE"]["P"]["helicity"]
              ["kappa"] = {{0.1, 0.2}, {0.3, 0.4}};
    CHECK_THROWS(LoadModel(asymmetric, "asymmetric helicity fixture"));
  }

  SECTION("wrong matrix dimension is rejected") {
    json wrong_dimension = valid;
    wrong_dimension["PARAM_SOFT"]["MODEL"][active_model]["EIKONAL"]
                   ["helicity"] = true;
    wrong_dimension["PARAM_SOFT"]["MODEL"][active_model]["EXCHANGE"]["P"]
                   ["helicity"]["B_kappa"] = {{0.1}};
    CHECK_THROWS(
        LoadModel(wrong_dimension, "wrong dimension helicity fixture"));
  }
}

TEST_CASE("SOFT eikonal selections and kernel domains are validated on load",
          "[gra::SoftModel][eikonal][validation]") {
  const json valid = LoadTuneDocument();
  const std::string active_model =
      valid.at("PARAM_SOFT").at("active_model").get<std::string>();

  SECTION("scalar screening selection") {
    json invalid = valid;
    invalid["PARAM_SOFT"]["MODEL"][active_model]["EIKONAL"]
           ["screening_exchanges"] = "P";
    CHECK_THROWS(LoadModel(invalid, "scalar screening fixture"));
  }
  SECTION("empty screening selection") {
    json invalid = valid;
    invalid["PARAM_SOFT"]["MODEL"][active_model]["EIKONAL"]
           ["screening_exchanges"] = json::array();
    CHECK_THROWS(LoadModel(invalid, "empty screening fixture"));
  }
  SECTION("screening wildcard selects every enabled exchange") {
    json wildcard = valid;
    wildcard["PARAM_SOFT"]["MODEL"][active_model]["EIKONAL"]
            ["screening_exchanges"] = {"*"};
    const auto model = LoadModel(wildcard, "screening wildcard fixture");
    std::size_t enabled = 0;
    for (const auto &exchange : model->Exchanges()) {
      enabled += exchange.Enabled() ? 1U : 0U;
    }
    CHECK(model->Eikonal().ScreeningExchanges().size() == enabled);
  }
  SECTION("wildcard cannot be combined with explicit exchanges") {
    json invalid = valid;
    invalid["PARAM_SOFT"]["MODEL"][active_model]["EIKONAL"]
           ["screening_exchanges"] = {"*", "P"};
    CHECK_THROWS(LoadModel(invalid, "mixed screening wildcard fixture"));
  }
  SECTION("helicity switch must be boolean") {
    json invalid = valid;
    invalid["PARAM_SOFT"]["MODEL"][active_model]["EIKONAL"]["helicity"] =
        "false";
    CHECK_THROWS(LoadModel(invalid, "nonboolean helicity fixture"));
  }
  SECTION("nonpositive exponential form factor slope") {
    json invalid = valid;
    invalid["PARAM_SOFT"]["MODEL"][active_model]["FF"]["R"]["param"][0][0] =
        0.0;
    CHECK_THROWS(LoadModel(invalid, "invalid EXP slope fixture"));
  }

  for (const auto &[coordinate, value] :
       std::array<std::pair<std::size_t, double>, 4>{
           std::pair{0, -1.0}, std::pair{1, 0.0}, std::pair{2, -1.0},
           std::pair{3, -1.0}}) {
    DYNAMIC_SECTION("invalid GKERNEL coordinate " << coordinate) {
      json invalid = valid;
      auto &model = invalid["PARAM_SOFT"]["MODEL"][active_model];
      model["FF"]["UNUSED"] = {
          {"type", "GKERNEL"},
          {"param", {{1.0, 1.0, 0.0, 0.0}, {1.0, 1.0, 0.0, 0.0}}}};
      model["FF"]["UNUSED"]["param"][0][coordinate] = value;
      CHECK_THROWS(LoadModel(invalid, "invalid GKERNEL fixture"));
    }
  }
  SECTION("GKERNEL requires a finite forward slope") {
    json invalid = valid;
    auto &model = invalid["PARAM_SOFT"]["MODEL"][active_model];
    model["FF"]["UNUSED"] = {
        {"type", "GKERNEL"},
        {"param", {{1.0, 0.5, 0.0, 0.0}, {1.0, 1.0, 0.0, 0.0}}}};
    CHECK_THROWS(LoadModel(invalid, "nonfinite GKERNEL slope fixture"));
  }
  SECTION("GKERNEL optional prefactor must be nonnegative") {
    json invalid = valid;
    auto &model = invalid["PARAM_SOFT"]["MODEL"][active_model];
    model["FF"]["UNUSED"] = {
        {"type", "GKERNEL"},
        {"param", {{-0.1, 1.0, 1.0, 0.0, 0.0}, {0.0, 1.0, 1.0, 0.0, 0.0}}}};
    CHECK_THROWS(LoadModel(invalid, "negative GKERNEL prefactor fixture"));
  }
}

// Require each soft exchange definition to match an active model row
TEST_CASE("SOFT exchange definitions and model rows must agree",
          "[gra::SoftModel][mapping][validation]") {
  for (const bool missing_definition : {false, true}) {
    CAPTURE(missing_definition);
    json invalid = LoadTuneDocument();
    auto &soft = invalid["PARAM_SOFT"];
    const auto active = soft.at("active_model").get<std::string>();
    auto &rows = missing_definition ? soft["EXCHANGE_DEF"] : soft["MODEL"][active]["EXCHANGE"];
    rows.erase("R_f2");
    CHECK_THROWS_WITH(LoadModel(invalid, "inconsistent exchange rows"),
                      Catch::Matchers::Contains(missing_definition ? "unknown exchange R_f2" : "missing R_f2"));
  }
}

TEST_CASE("Checked-in SOFT channel models load through one immutable parser",
          "[gra::SoftModel][mapping]") {
  for (const auto &[name, channels] :
       std::vector<std::pair<std::string, std::size_t>>{
           {"single", 1}, {"double", 2}, {"triple", 3}}) {
    CAPTURE(name, channels);
    json document = LoadTuneDocument();
    document["PARAM_SOFT"]["active_model"] = name;
    const auto model = LoadModel(document, "TUNE0 " + name);
    REQUIRE(model->ActiveModel() == name);
    REQUIRE(model->GoodWalker().ChannelCount() == channels);
    REQUIRE(model->GoodWalker().PairDimension() == channels * channels);
    CheckTriplePomeronRoot(*model, model->ExchangeId("P"));
  }
}

TEST_CASE("Generated four-channel Good Walker algebra is complete",
          "[gra::GoodWalker][N4]") {
  const auto model =
      LoadModel(FourChannelDocument(), "generated four-channel fixture");
  const auto &space = model->GoodWalker();
  REQUIRE(space.ChannelCount() == 4);
  REQUIRE(space.PairDimension() == 16);
  REQUIRE(space.MixingAngles().size() == 6);
  REQUIRE(space.ResolvedCoefficients().size() == 3);
  const double inverse_norm = 1.0 / std::sqrt(14.0);
  CHECK(space.ResolvedCoefficients()[0] == Approx(inverse_norm));
  CHECK(space.ResolvedCoefficients()[1] == Approx(2.0 * inverse_norm));
  CHECK(space.ResolvedCoefficients()[2] == Approx(3.0 * inverse_norm));

  for (std::size_t i = 0; i < space.ChannelCount(); ++i) {
    for (std::size_t k = 0; k < space.ChannelCount(); ++k) {
      CHECK(space.PairIndex(i, k) == 4 * i + k);
    }
  }
  CHECK_THROWS_AS(space.PairIndex(4, 0), std::out_of_range);
  CHECK_THROWS_AS(space.PairIndex(0, 4), std::out_of_range);

  const auto pomeron = model->ExchangeId("P");
  const double transition_t = -0.27;
  const double geometric =
      std::sqrt(model->FormFactor(pomeron, transition_t, 0) *
                model->FormFactor(pomeron, transition_t, 1));
  CHECK(model->TransitionFormFactor(pomeron, transition_t, 0, 1) ==
        Approx(geometric).margin(1.0e-14));
  CHECK(model->ResidueMatrix(pomeron, transition_t)(0, 1) ==
        Approx(model->Exchange(pomeron).Coupling(0, 1) * geometric)
            .margin(1.0e-14));

  CHECK(gra::InnerProduct(space.ProtonVector(), space.ProtonVector()) ==
        Approx(1.0).margin(1.0e-13));
  CHECK(gra::InnerProduct(space.ResolvedVector(), space.ResolvedVector()) ==
        Approx(1.0).margin(1.0e-13));
  CHECK(gra::InnerProduct(space.ProtonVector(), space.ResolvedVector()) ==
        Approx(0.0).margin(1.0e-13));

  const auto identity = Identity(4);
  CheckProjector(space.ProtonProjector());
  CheckProjector(space.ExcitedProjector());
  CheckProjector(space.ResolvedProjector());
  CheckProjector(space.InclusiveProjector());
  CheckMatrix(space.ProtonProjector() + space.ExcitedProjector(), identity);
  CheckMatrix(space.ResolvedProjector() + space.InclusiveProjector(), identity);
  CheckMatrix(Multiply(space.ExcitedProjector(), space.ResolvedProjector()),
              space.ResolvedProjector());

  std::vector<std::complex<double>> pair_source(space.PairDimension(), 0.0);
  for (std::size_t i = 0; i < pair_source.size(); ++i) {
    pair_source[i] = {static_cast<double>(i + 1), -0.25 * i};
  }
  const auto complete =
      space.ProjectPair(pair_source, gra::GoodWalkerFinalBasis::Complete,
                        gra::GoodWalkerFinalBasis::Complete);
  CHECK(complete == pair_source);
  REQUIRE(space
              .ProjectPair(pair_source, gra::GoodWalkerFinalBasis::Proton,
                           gra::GoodWalkerFinalBasis::Proton)
              .size() == 1);
  REQUIRE(space
              .ProjectPair(pair_source, gra::GoodWalkerFinalBasis::Excited,
                           gra::GoodWalkerFinalBasis::Excited)
              .size() == 9);

  CheckTriplePomeronRoot(*model, pomeron);
}

TEST_CASE("Good Walker sources are covariant under eigenstate permutations",
          "[gra::GoodWalker][N4][permutation]") {
  const auto model =
      LoadModel(FourChannelDocument(), "permuted four-channel fixture");
  const auto pomeron = model->ExchangeId("P");
  const auto &space = model->GoodWalker();
  const std::vector<std::size_t> order = {2, 0, 3, 1};

  const auto coupling = model->Exchange(pomeron).CouplingMatrix();
  const auto coupling_permuted = PermuteChannelMatrix(coupling, order);
  const auto root = model->TriplePomeronCouplingRoot(pomeron);
  const auto root_permuted = PermuteChannelMatrix(root, order);
  const double triple_coupling = model->EffectiveTriplePomeronCoupling(pomeron);
  CheckMatrix(Multiply(root_permuted.Transpose(), root_permuted),
              coupling_permuted * triple_coupling, 2.0e-9);

  const auto proton_permuted =
      PermuteChannelVector(space.ProtonVector(), order);
  const auto resolved_permuted =
      PermuteChannelVector(space.ResolvedVector(), order);
  CheckMatrix(PermuteChannelMatrix(space.ProtonProjector(), order),
              gra::RankOneProjector(proton_permuted));
  CheckMatrix(PermuteChannelMatrix(space.ResolvedProjector(), order),
              gra::RankOneProjector(resolved_permuted));

  const auto residue = model->ResidueMatrix(pomeron, -0.27);
  const auto residue_permuted = PermuteChannelMatrix(residue, order);
  double original_projection = 0.0;
  double permuted_projection = 0.0;
  for (std::size_t i = 0; i < order.size(); ++i) {
    for (std::size_t j = 0; j < order.size(); ++j) {
      original_projection +=
          space.ProtonVector()[i] * residue(i, j) * space.ProtonVector()[j];
      permuted_projection +=
          proton_permuted[i] * residue_permuted(i, j) * proton_permuted[j];
    }
  }
  CHECK(permuted_projection == Approx(original_projection).margin(2.0e-13));
}

TEST_CASE("Normalized beam residue is distinct from the bare coupling",
          "[gra::SoftModel][residue][3P]") {
  json document = LoadTuneDocument();
  document["PARAM_SOFT"]["active_model"] = "double";
  auto &model_row = document["PARAM_SOFT"]["MODEL"]["double"];
  model_row["GW"]["theta"] = {0.41};
  model_row["EXCHANGE"]["P"]["g"] = {{4.0, 1.25}, {1.25, 2.0}};
  model_row["EXCHANGE"]["P"]["transition_ff"] = "diagonal";

  const auto model = LoadModel(document, "diagonal transition fixture");
  const auto pomeron = model->ExchangeId("P");
  const auto &proton = model->GoodWalker().ProtonVector();
  const auto &coupling = model->Exchange(pomeron).CouplingMatrix();

  double bare_reference = 0.0;
  double residue_reference = 0.0;
  for (std::size_t i = 0; i < proton.size(); ++i) {
    for (std::size_t j = 0; j < proton.size(); ++j) {
      bare_reference += proton[i] * coupling(i, j) * proton[j];
      if (i == j) {
        residue_reference += proton[i] * coupling(i, i) * proton[i];
      }
    }
  }
  REQUIRE(std::abs(bare_reference - residue_reference) > 1.0e-3);
  CHECK(model->PhysicalCoupling(pomeron) ==
        Approx(bare_reference).margin(1.0e-14));
  CHECK(model->PhysicalResidue(pomeron, 0.0) ==
        Approx(residue_reference).margin(1.0e-14));
  CHECK(model->NormalizedPhysicalResidue(pomeron, 0.0) ==
        Approx(1.0).margin(1.0e-14));

  const double t = -0.31;
  double residue_t = 0.0;
  for (std::size_t i = 0; i < proton.size(); ++i) {
    residue_t += proton[i] * coupling(i, i) * proton[i] *
                 model->FormFactor(pomeron, t, i);
  }
  CHECK(model->NormalizedPhysicalResidue(pomeron, t) ==
        Approx(residue_t / residue_reference).margin(1.0e-14));
  CHECK(model->EffectiveTriplePomeronCoupling(pomeron) ==
        Approx(model->TriplePomeronRatio() * bare_reference).margin(1.0e-14));
  CheckTriplePomeronRoot(*model, pomeron);
}

TEST_CASE("Resolved Good Walker direction is scale invariant",
          "[gra::GoodWalker][validation]") {
  const std::vector<double> angles = {0.17, -0.11, 0.23};
  const gra::GoodWalkerSpace first(3, angles, {0.6, -0.8});
  for (const double scale : {1e-300, 1.0, 1e300}) {
    CAPTURE(scale);
    const gra::GoodWalkerSpace scaled(3, angles, {3.0 * scale, -4.0 * scale});
    for (const auto &i : indices(first.ResolvedCoefficients())) {
      CHECK(first.ResolvedCoefficients()[i] ==
            Approx(scaled.ResolvedCoefficients()[i]).margin(1.0e-15));
    }
    for (const auto &i : indices(first.ResolvedVector())) {
      CHECK(first.ResolvedVector()[i] ==
            Approx(scaled.ResolvedVector()[i]).margin(1.0e-15));
    }
    CheckMatrix(first.ResolvedProjector(), scaled.ResolvedProjector());
    CheckMatrix(first.InclusiveProjector(), scaled.InclusiveProjector());
  }

  CHECK_THROWS(gra::GoodWalkerSpace(3, angles, {0.0, 0.0}));
  CHECK_THROWS(gra::GoodWalkerSpace(
      3, angles, {1.0, std::numeric_limits<double>::infinity()}));
}

TEST_CASE("One-channel Good Walker limit is wholly inclusive",
          "[gra::GoodWalker][N1]") {
  const gra::GoodWalkerSpace space(1, {}, {});
  REQUIRE(space.ChannelCount() == 1);
  REQUIRE(space.PairDimension() == 1);
  CHECK(space.PairIndex(0, 0) == 0);
  REQUIRE(space.ProtonVector() == std::vector<double>{1.0});
  REQUIRE(space.ResolvedVector() == std::vector<double>{0.0});
  CHECK(space.ProtonProjector()(0, 0) == Approx(1.0));
  CHECK(space.ExcitedProjector()(0, 0) == Approx(0.0));
  CHECK(space.ResolvedProjector()(0, 0) == Approx(0.0));
  CHECK(space.InclusiveProjector()(0, 0) == Approx(1.0));
  CHECK(space.FinalBasis(gra::GoodWalkerFinalBasis::Excited).size_col() == 0);

  const std::vector<std::complex<double>> source = {{2.0, -3.0}};
  CHECK(space.ProjectPair(source, gra::GoodWalkerFinalBasis::Complete,
                          gra::GoodWalkerFinalBasis::Complete) == source);
  CHECK(space
            .ProjectPair(source, gra::GoodWalkerFinalBasis::Excited,
                         gra::GoodWalkerFinalBasis::Complete)
            .empty());
  CHECK_THROWS(
      gra::GoodWalkerSpace(std::numeric_limits<std::size_t>::max(), {}, {}));
}

TEST_CASE("Triple Pomeron principal root rejects a non-PSD coupling",
          "[gra::SoftModel][3P][validation]") {
  json document = LoadTuneDocument();
  document["PARAM_SOFT"]["active_model"] = "double";
  document["PARAM_SOFT"]["MODEL"]["double"]["EXCHANGE"]["P"]["g"] = {
      {1.0, 2.0}, {2.0, 1.0}};

  const auto model = LoadModel(document, "indefinite Pomeron fixture");
  CHECK_THROWS(model->TriplePomeronCouplingRoot(model->ExchangeId("P")));
}

TEST_CASE("Triple Pomeron root handles zero and numerical null directions",
          "[gra::SoftModel][3P][validation]") {
  SECTION("zero coupling matrix") {
    json document = LoadTuneDocument();
    document["PARAM_SOFT"]["active_model"] = "single";
    document["PARAM_SOFT"]["MODEL"]["single"]["EXCHANGE"]["P"]["g"] = {{0.0}};
    const auto model = LoadModel(document, "zero Pomeron coupling fixture");
    const auto pomeron = model->ExchangeId("P");
    REQUIRE(model->EffectiveTriplePomeronCoupling(pomeron) == Approx(0.0));
    const auto root = model->TriplePomeronCouplingRoot(pomeron);
    REQUIRE(root.size_row() == 1);
    REQUIRE(root.size_col() == 1);
    REQUIRE(root(0, 0) == Approx(0.0));
  }

  SECTION("roundoff-scale negative eigenvalue") {
    json document = LoadTuneDocument();
    document["PARAM_SOFT"]["active_model"] = "double";
    document["PARAM_SOFT"]["MODEL"]["double"]["EXCHANGE"]["P"]["g"] = {
        {1.0, 1.0}, {1.0, 1.0 - 1.0e-15}};
    const auto model =
        LoadModel(document, "numerical Pomeron null direction fixture");
    const auto pomeron = model->ExchangeId("P");
    const double scale = model->EffectiveTriplePomeronCoupling(pomeron);
    const auto target = model->Exchange(pomeron).CouplingMatrix() * scale;
    REQUIRE(target(0, 0) * target(1, 1) - target(0, 1) * target(1, 0) < 0.0);
    const auto root = model->TriplePomeronCouplingRoot(pomeron);
    const auto reconstructed = Multiply(root.Transpose(), root);
    CheckMatrix(root, root.Transpose(), 1.0e-13);
    CheckMatrix(reconstructed, target, 1.0e-13);
  }
}
