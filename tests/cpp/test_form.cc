// Form factor and model parameter tests
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "support/models_test_support.hh"

#include "Graniitti/Regge/MReggeSig.h"
#include "Graniitti/Tech/MException.h"

// Compare the approximate spin propagators with the published MadWidth table
TEST_CASE("MadWidth propagator estimates use the published spin numerators", "[gra::resonance][physics]") {
  // [REFERENCE: Alwall et al., arXiv:1402.1178v2, Table 1]
  const double mass = 4.0;
  const double width = 0.5;
  const std::array<std::array<double, 6>, 3> reference = {{{2.0, 1.0, 2.0, 0.75, 1.0, 0.875},
                                                           {3.0, 1.0, 3.0, 0.4375, 0.875, 241.0 / 384.0},
                                                           {6.0, 1.0, 6.0, -1.25, -5.0, 37.0 / 24.0}}};
  for (const auto &row : reference) {
    const double energy2 = row[0] * row[0];
    const std::complex<double> denominator(energy2 - mass * mass, mass * width);
    for (int spinX2 = 0; spinX2 <= 4; ++spinX2) {
      CAPTURE(row[0], spinX2);
      RequireComplexNear(gra::resonance::JacksonLineShape(energy2, mass, width, 0.5 * spinX2),
                         row[spinX2 + 1] / denominator);
    }
  }
}

// Keep invalid line shapes visible to the central amplitude failure bookkeeping
TEST_CASE("Resonance square roots preserve invalid kinematics", "[gra::resonance][physics]") {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  CHECK_FALSE(std::isfinite(gra::resonance::BreitWignerAmplitude(nan, 1.0, 0.1)));
  CHECK_FALSE(std::isfinite(gra::resonance::BreitWignerAmplitude(1.0, 1.0, 0.0)));
  for (const double spin : {0.5, 1.5}) {
    CHECK_FALSE(std::isfinite(std::abs(gra::resonance::JacksonLineShape(-1.0, 1.0, 0.1, spin))));
  }
}

// Check physical partial widths through both orderings of the decay daughters
TEST_CASE("Electronic partial widths accept exchanged decay daughters", "[gra::form][gra::resonance]") {
  const auto dir = std::filesystem::path("tmp") / "test_partial_width_order";
  std::filesystem::create_directories(dir);
  gra::MParticle resonance;
  resonance.pdg = 443;
  resonance.width = 0.002;
  const auto write = [&](const nlohmann::json &table) {
    std::ofstream output(dir / "DECAYS.json");
    REQUIRE(output.is_open());
    output << table.dump();
  };
  for (const std::string key : {"[11,-11]", "[-11,11]"}) {
    CAPTURE(key);
    nlohmann::json table;
    table["443"][key]["BR"] = 0.04;
    write(table);
    CHECK(gra::resonance::ElectronicPartialWidth(resonance, dir.string()) == Approx(0.00008));
    for (const auto &value : std::vector<nlohmann::json>{0.0, -0.1, 1.1, true, "0.04"}) {
      table["443"][key]["BR"] = value;
      write(table);
      REQUIRE_THROWS_AS(gra::resonance::ElectronicPartialWidth(resonance, dir.string()), std::invalid_argument);
    }
  }
  nlohmann::json missing;
  missing["443"]["[13,-13]"]["BR"] = 0.04;
  write(missing);
  REQUIRE_THROWS_AS(gra::resonance::ElectronicPartialWidth(resonance, dir.string()), std::invalid_argument);
}

// Check shared exponential helpers connect amplitude and cross-section factors
TEST_CASE("Exponential slope helpers use one cross-section convention",
          "[gra::form][slope]") {
  const double B = 4.7;
  const double delta_t = -0.35;
  const double amplitude = gra::form::ExpSlopeAmplitude(B, delta_t);
  const double weight = gra::form::ExpSlopeWeight(B, delta_t);

  CHECK(amplitude == Approx(std::exp(0.5 * B * delta_t)).margin(1e-15));
  CHECK(weight == Approx(std::exp(B * delta_t)).margin(1e-15));
  CHECK(gra::math::pow2(amplitude) == Approx(weight).margin(1e-15));
}

// Check the flat process applies its configured slope directly to the weight
TEST_CASE("Model-card loaders reject malformed physics parameters",
          "[gra::form][params][validation]") {
  nlohmann::json card =
      nlohmann::json::parse(gra::aux::GetInputData(modelfile));
  card.at("PARAM_FLAT").erase("B");
  REQUIRE_THROWS_AS(
      gra::form::FlatParam::Read("missing B", card.dump()),
      std::invalid_argument);
  for (const auto &value : {nlohmann::json(-1.0), nlohmann::json(true), nlohmann::json("4.0"), nlohmann::json(nullptr)}) {
    card["PARAM_FLAT"]["B"] = value;
    REQUIRE_THROWS_AS(gra::form::FlatParam::Read("invalid B", card.dump()), std::invalid_argument);
  }

  nlohmann::json invalid_nstar =
      nlohmann::json::parse(gra::aux::GetInputData(modelfile));
  gra::MNstarParam nstar;
  for (const double p : {-0.1, 0.0, 1.0, 1.1}) {
    CAPTURE(p);
    invalid_nstar.at("PARAM_NSTAR").at("single_side_prob") = p;
    REQUIRE_THROWS_AS(nstar.Configure(invalid_nstar.at("PARAM_NSTAR"),
                                     "invalid N-star probability"), std::invalid_argument);
  }
  for (const double p : {0.2, 0.8}) {
    invalid_nstar.at("PARAM_NSTAR").at("single_side_prob") = p;
    REQUIRE_NOTHROW(nstar.Configure(invalid_nstar.at("PARAM_NSTAR"), "asymmetric N-star sampling"));
  }

  nlohmann::json invalid_structure =
      nlohmann::json::parse(gra::aux::GetInputData(modelfile));
  invalid_structure.at("PARAM_STRUCTURE").at("F2") = "unknown";
  REQUIRE_THROWS(gra::form::ParamStore::Read("invalid structure model",
                                             invalid_structure.dump()));
}

// Check log squared normalization, curvature, asymptotic bounds and input domains
TEST_CASE("LogExp form factors match exponential curvature", "[gra::form][logexp]") {
  using gra::regge::ReadFF;
  using gra::regge::FormFactor;
  for (const double b : {0.2, 0.9, 2.0}) {
    for (const double scale2 : {0.3, 1.0, 3.0}) {
      const auto ff = ReadFF({{"type", "logexp"}, {"norm", "pole"}, {"b", b}, {"Lambda2", scale2}}, "test");
      CHECK(FormFactor(1.0, 1.0, ff) == Approx(1.0));
      const double h = 1.0e-3 * scale2;
      const double logf = std::log(FormFactor(1.0 - h, 1.0, ff));
      CHECK(-logf / h == Approx(b).epsilon(2.0e-7));
      CHECK((logf + b * h) / (b * h * h * h / (6.0 * scale2 * scale2)) == Approx(1.0).epsilon(2.0e-3));
      double previous = 1.0;
      for (int i = 1; i <= 10000; ++i) {
        const double x = 20.0 * scale2 * i / 10000.0;
        const double value = FormFactor(1.0 - x, 1.0, ff);
        CHECK(value >= 0.0);
        CHECK(value <= previous);
        CHECK(value + 1.0e-14 >= std::exp(-b * x));
        CHECK(value <= std::exp(-b * scale2 * std::log1p(x / scale2)) + 1.0e-14);
        previous = value;
      }
    }
  }
  const nlohmann::json card = {{"type", "logexp"}, {"norm", "zero"}, {"b", 1.0}, {"Lambda2", 1.0}};
  const auto ff = ReadFF(card, "test");
  CHECK(gra::regge::TransferFF(0.0, ff) == Approx(1.0));
  CHECK(gra::regge::TransferFF(-1.0e200, ff) == Approx(0.0).margin(1.0e-300));
  CHECK(gra::regge::TransferFF(0.5, ff) == Approx(std::exp(std::log(2.0) - 0.5 * std::pow(std::log(2.0), 2))));
  CHECK_THROWS_AS(gra::regge::TransferFF(1.0, ff), gra::AmplitudeFailure);
  CHECK_THROWS_AS(gra::regge::TransferFF(2.0, ff), gra::AmplitudeFailure);
  CHECK_THROWS_AS(gra::regge::TransferFF(std::numeric_limits<double>::quiet_NaN(), ff), gra::AmplitudeFailure);
  auto changed = card;
  changed["b"] = 0.0;
  CHECK(gra::regge::TransferFF(-10.0, ReadFF(changed, "test")) == Approx(1.0));
  for (const std::string field : {"b", "Lambda2"}) {
    for (const double invalid : {-1.0, std::numeric_limits<double>::infinity()}) {
      changed = card;
      changed[field] = invalid;
      CHECK_THROWS_AS(ReadFF(changed, "test"), std::invalid_argument);
    }
    changed = card;
    changed.erase(field);
    CHECK_THROWS_AS(ReadFF(changed, "test"), std::invalid_argument);
  }
  changed = card;
  changed["Lambda2"] = 0.0;
  CHECK_THROWS_AS(ReadFF(changed, "test"), std::invalid_argument);
}

// Check the normalized named form factor families
TEST_CASE("FormFactor uses normalized named families", "[gra::form]") {
  const double t_hat = -0.7;
  const double M2 = 1.5;
  const double x = M2 - t_hat;
  const double B = 0.4;

  const gra::regge::FFParam none;
  CHECK(gra::regge::FormFactor(t_hat, M2, none) == 1.0);
  CHECK(gra::regge::TransferFF(t_hat, none) == 1.0);
  CHECK(gra::regge::MassFF(M2, M2, none) == 1.0);

  const auto exponential = gra::regge::FFParam{
      gra::regge::FFType::Exponential, gra::regge::FFNorm::Pole, {B}};
  CHECK(gra::regge::FormFactor(t_hat, M2, exponential) ==
        Approx(std::exp(-B * x)).margin(1e-15));

  const double power_scale2 = 2.5;
  const auto power_one = gra::regge::FFParam{
      gra::regge::FFType::Power, gra::regge::FFNorm::Pole, {power_scale2, 1.0}};
  CHECK(gra::regge::FormFactor(t_hat, M2, power_one) ==
        Approx(1.0 / (1.0 + x / power_scale2)).margin(1e-15));

  const double b_or = 0.7;
  const double a_or = 0.25;
  const auto orear = gra::regge::FFParam{
      gra::regge::FFType::Orear, gra::regge::FFNorm::Pole, {b_or, a_or}};
  CHECK(gra::regge::FormFactor(t_hat, M2, orear) ==
        Approx(std::exp(-b_or * (std::sqrt(x + a_or * a_or) - a_or)))
            .margin(1e-15));

  const double power_scale2_n = 1.7;
  const double n = 2.2;
  const auto power = gra::regge::FFParam{
      gra::regge::FFType::Power, gra::regge::FFNorm::Pole, {power_scale2_n, n}};
  CHECK(gra::regge::FormFactor(t_hat, M2, power) ==
        Approx(std::pow(1.0 + x / power_scale2_n, -n)).margin(1e-15));

  const double vector_scale2 = 4.0;
  const double vector_n = 0.5;
  const auto vector = gra::regge::FFParam{gra::regge::FFType::Vector,
                                          gra::regge::FFNorm::Pole,
                                          {vector_scale2, vector_n}};
  const double vector_base =
      1.0 + t_hat * (t_hat - M2) / gra::math::pow2(vector_scale2);
  CHECK(gra::regge::FormFactor(t_hat, M2, vector) ==
        Approx(std::pow(vector_base, -vector_n)).margin(1e-15));
  CHECK(gra::regge::FormFactor(0.0, M2, vector) == Approx(1.0));
  CHECK(gra::regge::FormFactor(M2, M2, vector) == Approx(1.0));

  const auto kernel = [](std::vector<double> param) {
    return gra::regge::FFParam{gra::regge::FFType::GKernel,
                               gra::regge::FFNorm::Pole, std::move(param)};
  };
  const auto named_kernel = gra::regge::ReadFF(
      {{"type", "gkernel"}, {"norm", "pole"},
       {"terms", {{{"a", b_or}, {"p", 0.5}, {"nu", 0.0}, {"mu2", a_or * a_or}}}}}, "gkernel");
  CHECK(gra::regge::FormFactor(t_hat, M2, named_kernel) ==
        Approx(gra::regge::FormFactor(t_hat, M2, orear)).margin(1e-15));
  CHECK(gra::regge::FormFactor(t_hat, M2, kernel({B, 1.0, 0.0, 0.0})) ==
        Approx(gra::regge::FormFactor(t_hat, M2, exponential)).margin(1e-15));
  CHECK(gra::regge::FormFactor(t_hat, M2,
                               kernel({1.0 / power_scale2, 1.0, 1.0, 0.0})) ==
        Approx(gra::regge::FormFactor(t_hat, M2, power_one)).margin(1e-15));
  CHECK(gra::regge::FormFactor(
            t_hat, M2, kernel({n / power_scale2_n, 1.0, 1.0 / n, 0.0})) ==
        Approx(gra::regge::FormFactor(t_hat, M2, power)).margin(1e-15));
  CHECK(gra::regge::FormFactor(t_hat, M2,
                               kernel({b_or, 0.5, 0.0, a_or * a_or})) ==
        Approx(gra::regge::FormFactor(t_hat, M2, orear)).margin(1e-15));
  CHECK(gra::regge::FormFactor(t_hat, M2,
                               kernel({B, 1.0, 0.0, 0.0, n / power_scale2_n,
                                       1.0, 1.0 / n, 0.0})) ==
        Approx(gra::regge::FormFactor(t_hat, M2, exponential) *
               gra::regge::FormFactor(t_hat, M2, power))
            .margin(1e-15));

  const double hard_x = 1e6;
  const double hard_t = M2 - hard_x;
  const double nu = 2.5;
  const double a = 0.7;
  const double scaled_tail =
      gra::regge::FormFactor(hard_t, M2, kernel({a, 1.0, nu, 0.0})) *
      std::pow(hard_x, 1.0 / nu);
  CHECK(scaled_tail == Approx(std::pow(nu * a, -1.0 / nu)).epsilon(2e-6));

  for (const char *removed :
       {"exponential", "QEXP", "KAPPA", "STEXP", "SPECMIX"}) {
    CAPTURE(removed);
    CHECK_THROWS_AS(gra::regge::ParseFFType(removed), std::invalid_argument);
  }

  CHECK_THROWS(gra::regge::FormFactor(2.0, 1.0, kernel({1.0, 1.0, 1.0, 0.0})));
  CHECK_THROWS_AS(
      gra::regge::FormFactor(std::numeric_limits<double>::quiet_NaN(), M2,
                             kernel({1.0, 1.0, 1.0, 0.0})),
      gra::AmplitudeFailure);
}

// Check the gauged meson factors at coincident virtualities, poles and large spacelike momenta
TEST_CASE("Meson form factor differences retain their analytic limit", "[gra::form][physics]") {
  const double pole = 0.8;
  for (const std::string type : {"none", "exp", "gaussian"}) {
    nlohmann::json block = {{"type", type}};
    if (type != "none") {
      block["norm"] = "pole";
      block[type == "exp" ? "b" : "Lambda2"] = 1.4;
    }
    const auto ff = gra::regge::ReadFF(block, "test divided meson form");
    for (const double q1 : {-100.0, -0.4, pole, 1.1}) {
      for (const double step : {0.0, 1.0e-13, 0.2}) {
        const double q2 = q1 + step;
        const auto f = gra::regge::MassFFPair(q1, q2, pole, ff);
        const auto crossed = gra::regge::MassFFPair(q2, q1, pole, ff);
        CAPTURE(type, q1, q2);
        CHECK(f[0] == Approx(gra::regge::MassFF(q1, pole, ff)).margin(1.0e-14));
        CHECK(f[1] == Approx(gra::regge::MassFF(q2, pole, ff)).margin(1.0e-14));
        CHECK(f[2] == Approx(crossed[2]).margin(1.0e-14));
        const double derivative = type == "none" ? 0.0 : type == "exp" ? ff.param[0] : -2.0 * (q1 - pole) / pow2(ff.param[0]);
        const double reference = step < 1.0e-10 ? derivative * f[0] : (f[0] - f[1]) / (q1 - q2);
        CHECK(f[2] == Approx(reference).margin(1.0e-12));
        CHECK(f[0] - f[1] == Approx((q1 - q2) * f[2]).margin(1.0e-14));
      }
    }
  }
}

TEST_CASE("Soft form factor GKERNEL covers product profiles", "[gra::form]") {
  const auto source = nlohmann::json::parse(gra::aux::GetInputData(modelfile));
  const double t = -0.7;
  const double x = -t;
  const auto evaluate = [&source, t](const std::vector<double> &parameters) {
    auto card = source;
    const std::string active = card.at("PARAM_SOFT").at("active_model");
    auto &profile =
        card.at("PARAM_SOFT").at("MODEL").at(active).at("FF").at("P");
    const std::size_t channels = profile.at("param").size();
    profile.at("type") = "GKERNEL";
    profile.at("param") = nlohmann::json::array();
    for (std::size_t channel = 0; channel < channels; ++channel) {
      profile.at("param").push_back(parameters);
    }
    const auto model =
        gra::SoftModel::LoadFromJson("in-memory GKERNEL fixture", card.dump());
    return model->FormFactor(model->ExchangeId("P"), t, 0);
  };

  const double b = 0.777;
  const double c = 1.216;
  const double d = 1.425;
  CHECK(evaluate({d / b, 1.0, 1.0 / d, 0.0, d / c, 1.0, 1.0 / d, 0.0}) ==
        Approx(std::pow((1.0 / (1.0 - t / b)) * (1.0 / (1.0 - t / c)), d))
            .margin(1e-15));

  const double bexp = 8.501;
  const double cexp = 0.25;
  const double dexp = 0.45;
  CHECK(evaluate({std::pow(bexp, dexp), dexp, 0.0, cexp}) ==
        Approx(std::exp(-std::pow(bexp * (cexp - t), dexp) +
                        std::pow(bexp * cexp, dexp)))
            .margin(1e-15));

  const double B = 6.0;
  const double r = 2.3;
  const double m2 = 5.0;
  CHECK(evaluate({r, 0.5 * B, 1.0, 0.0, 0.0, 2.0 / m2, 1.0, 0.5, 0.0}) ==
        Approx(std::exp(-0.5 * B * x) * (1.0 + r * x) /
               gra::math::pow2(1.0 + x / m2))
            .margin(1e-15));
}

TEST_CASE("Inelastic proton profile reads its structure parameters",
          "[gra::form][gra::form::ParamStore]") {
  const auto source = nlohmann::json::parse(gra::aux::GetInputData(modelfile));
  const auto load_model = [&](double s0, double a, const std::string &suffix) {
    auto parameters = source;
    parameters["PARAM_SOFT"]["FORWARD_EXCITATION"]["s0"] = s0;
    parameters["PARAM_SOFT"]["FORWARD_EXCITATION"]["a"] = a;
    return gra::SoftModel::LoadFromJson("FORWARD_EXCITATION " + suffix,
                                        parameters.dump());
  };

  SECTION("configured profile") {
    const double s0 = 2.5;
    const double a = 0.75;
    const double t = -0.2;
    const double M2 = 3.0;
    const auto model = load_model(s0, a, "configured");
    const auto pomeron = model->ExchangeId("P");
    REQUIRE(model->ForwardExcitationProfile().S0() == Approx(s0));
    REQUIRE(model->ForwardExcitationProfile().A() == Approx(a));
    const double pomeron_intercept = model->Exchange(pomeron).Alpha0();
    const double expected = std::pow(
        (s0 * std::abs(t)) / (M2 * (std::abs(t) + a)), 0.5 * pomeron_intercept);
    CHECK(model->ForwardExcitationFactor(pomeron, t, M2) ==
          Approx(expected).epsilon(1e-14));
    CHECK(model->ForwardExcitationFactor(pomeron, 0.0, M2) == Approx(0.0));
    CHECK_THROWS_AS(model->ForwardExcitationFactor(pomeron, t, 0.0),
                    gra::AmplitudeFailure);
  }

  SECTION("nonpositive scales") {
    REQUIRE_THROWS(load_model(0.0, 0.75, "invalid_s0"));

    REQUIRE_THROWS(load_model(2.5, -0.75, "invalid_a"));
  }
}

// Check the complex propagator with a smooth rapidity gap and configurable freezing virtuality
TEST_CASE("OffshellProp reggeizes with supplied subchannel rapidity gap",
          "[gra::form]") {
  const double t_hat = -0.6;
  const double M2 = 0.25;
  const gra::M4Vec left(0.1, 0.0, 0.5, 1.0);
  const gra::M4Vec right(-0.1, 0.0, -0.2, 0.8);
  const double base = 1.0 / (t_hat - M2);
  const double alpha = 0.7 * (std::max(-1.0, t_hat) - M2);
  const double dY = std::abs(left.Rap() - right.Rap());
  const double gap = 2.0 * std::log(std::cosh(0.5 * dY));
  const gra::regge::MesonTraj pion_trajectory{0.0, 0.7};

  const std::complex<double> nonreg =
      gra::regge::OffshellProp(t_hat, M2, {false, 1.0}, pion_trajectory, left, right);
  CHECK(nonreg.real() == Approx(base).margin(1e-15));
  CHECK(nonreg.imag() == Approx(0.0).margin(1e-15));

  const std::complex<double> reg =
      gra::regge::OffshellProp(t_hat, M2, {true, 1.0}, pion_trajectory, left, right);
  CHECK(reg.real() == Approx(base * std::exp(alpha * gap)).margin(1e-15));
  CHECK(reg.imag() == Approx(0.0).margin(1e-15));
  for (const double scale2 : {0.3, 1.0, 2.0}) {
    for (const double t : {-0.2, -0.8, -3.0}) {
      const auto value = gra::regge::OffshellProp(t, M2, {true, scale2}, pion_trajectory, left, right);
      const double expected = std::exp(0.7 * (std::max(-scale2, t) - M2) * gap) / (t - M2);
      CHECK(std::abs(value - std::complex<double>(expected, 0.0)) < 1e-14);
      const auto equal_y = gra::regge::OffshellProp(t, M2, {true, scale2}, pion_trajectory, left, left);
      CHECK(std::abs(equal_y - std::complex<double>(1.0 / (t - M2), 0.0)) < 1e-14);
    }
  }
}

// Check the exchanged meson trajectory in the complex propagator
TEST_CASE("OffshellProp uses the exchanged meson trajectory",
          "[gra::form][physics]") {
  const auto configured = gra::regge::ReadParam(
      {211, -211}, LoadedPDGTable(), *gra::MModelTune::Load(modelfile));
  auto param = configured;

  const double t_hat = -0.4;
  const gra::M4Vec left(0.2, 0.0, 0.6, 1.1);
  const gra::M4Vec right(-0.1, 0.0, -0.3, 0.9);
  const double dY = std::abs(left.Rap() - right.Rap());
  const double gap = 2.0 * std::log(std::cosh(0.5 * dY));
  const double pion_m2 = gra::math::pow2(0.13957);
  const double rho_m2 = gra::math::pow2(0.77526);

  const auto pion =
      gra::regge::OffshellProp(t_hat, pion_m2, {true, 1.0}, gra::regge::MesonTrajectory(param, 211), left, right);
  const auto rho =
      gra::regge::OffshellProp(t_hat, rho_m2, {true, 1.0}, gra::regge::MesonTrajectory(param, 113), left, right);
  const double pion_expected =
      std::exp(0.70 * (t_hat - pion_m2) * gap) / (t_hat - pion_m2);
  const double rho_expected =
      std::exp(0.90 * (t_hat - rho_m2) * gap) / (t_hat - rho_m2);

  CHECK(pion.real() == Approx(pion_expected).margin(1e-14));
  CHECK(rho.real() == Approx(rho_expected).margin(1e-14));
  CHECK(std::abs(pion - rho) > 1e-3);
  REQUIRE_THROWS(
      gra::regge::MesonTrajectory(param, 999999));
}

// Check smooth gap limits, the pole coefficient and exchange symmetry through the complex propagator
TEST_CASE("OffshellProp smooth gap removes the linear rapidity cusp", "[gra::form][physics]") {
  const gra::regge::MesonTraj trajectory{0.0, 0.7};
  const gra::M4Vec rest(0.0, 0.0, 0.0, 1.0);
  for (const double dy : {0.0, 1.0e-5, 0.1, 0.99, 1.01, 4.0, 14.0}) {
    const gra::M4Vec left(0.0, 0.0, std::sinh(0.5 * dy), std::cosh(0.5 * dy));
    const gra::M4Vec right(0.0, 0.0, -left.Pz(), left.E());
    const long double d = left.Rap() - right.Rap();
    const double gap = static_cast<double>(2.0L * std::log(std::cosh(0.5L * d)));
    const auto value = gra::regge::OffshellProp(-0.6, 0.25, {true, 1.0}, trajectory, left, right);
    RequireComplexNear(value, -std::exp(-0.7 * 0.85 * gap) / 0.85, 1.0e-13);
    RequireComplexNear(value, gra::regge::OffshellProp(-0.6, 0.25, {true, 1.0}, trajectory, right, left), 1.0e-13);
    if (dy > 10.0) { CHECK(gap == Approx(static_cast<double>(d) - 2.0 * std::log(2.0)).margin(2.0e-6)); }
    if (dy > 0.0 && dy < 1.0e-4) {
      CHECK(-std::log(-0.85 * value.real()) / (dy * dy) == Approx(0.7 * 0.85 / 4.0).epsilon(2.0e-5));
    }
    const double offset = 1.0e-7;
    const auto pole = gra::regge::OffshellProp(0.25 - offset, 0.25, {true, 1.0}, trajectory, left, right);
    RequireComplexNear((0.25 - offset - 0.25) * pole, 1.0, 1.0e-6);
  }
  const gra::M4Vec invalid(0.0, 0.0, 2.0, 1.0);
  CHECK_THROWS_AS(gra::regge::OffshellProp(-0.6, 0.25, {true, 1.0}, trajectory, invalid, rest), gra::AmplitudeFailure);
}

// Check compact coherent spin rows complete their parity reflection
TEST_CASE("resonance::Read infers parity-conjugate coherent a_Jz rows",
          "[gra::resonance]") {
  const std::string old_modelparam = gra::MODELPARAM;
  struct RestoreModelParam {
    std::string value;
    ~RestoreModelParam() { gra::MODELPARAM = value; }
  } restore{old_modelparam};

  const std::filesystem::path tune_dir =
      std::filesystem::path("tmp") / "testbench4_ajz_parity";
  std::filesystem::create_directories(tune_dir / "RES");
  const std::filesystem::path card = tune_dir / "RES" / "compact_f2.json";

  std::filesystem::copy_file(
      gra::ResolveModelDataFile("TUNE0", "PDG_EXTRA.json"),
      tune_dir / "PDG_EXTRA.json",
      std::filesystem::copy_options::overwrite_existing);

  std::ofstream out(card);
  out << R"json(
{
  "PARAM_RES": {
    "name": "RES_compact_f2",
    "PDG": 225,
    "spinX2": 4,
    "P": 1,
    "C": 1,
    "MODELS": {
      "GP": {
        "mass": 1.2754,
        "width": 0.1866,
        "BW": "kinematic-width",
        "phi": 0.0,
        "FF_transfer": {"type": "none"},
        "FF_prod": {"type": "none"},
        "[990,990]": {
          "basis": "g_ls",
          "Lambda": 1.0,
          "g_ls": [[0, 2, 0.35, 0.0]],
          "CP": [true, true]
        }
      },
      "XP": {
        "mass": 1.2754,
        "width": 0.1866,
        "BW": "kinematic-width",
        "phi": 0.0,
        "FF_transfer": {"type": "none"},
        "FF_prod": {"type": "none"},
        "[991,991]": {
          "basis": "g_ls",
          "Lambda": 1.0,
          "g_ls": [[2, 0, 0.35, 0.0]],
          "CP": [true, true]
        }
      },
      "MP": {
        "mass": 1.2754,
        "width": 0.1866,
        "BW": "kinematic-width",
        "phi": 0.0,
        "FF_transfer": {"type": "none"},
        "FF_prod": {"type": "none"},
        "[991,991]": {
          "basis": "auto_min_L",
          "Lambda": 1.0000,
          "g": [0.35, 0.0],
          "polarization": {
            "mode": "a_Jz",
            "a_Jz": [[-2, 0.4, 0.15],
                     [-1, 0.3, -0.2],
                     [0, 0.7071067811865476, 0.0]],
            "rho_mag": [[0.2, 0, 0, 0, 0],
                        [0, 0.2, 0, 0, 0],
                        [0, 0, 0.2, 0, 0],
                        [0, 0, 0, 0.2, 0],
                        [0, 0, 0, 0, 0.2]],
            "rho_phase": [[0, 0, 0, 0, 0],
                          [0, 0, 0, 0, 0],
                          [0, 0, 0, 0, 0],
                          [0, 0, 0, 0, 0],
                          [0, 0, 0, 0, 0]],
            "random_rho": false
          },
          "CP": [true, true]
        }
      },
      "TP": {
        "mass": 1.2754,
        "width": 0.1866,
        "phi": 0.0,
        "FF_transfer": {"type": "none"},
        "FF_prod": {"type": "none"},
        "[995,995]": {
          "g_tensor": [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        }
      }
    }
  }
}
)json";
  out.close();

  gra::MODELPARAM = tune_dir.string();
  MRandom rng;
  const gra::PARAM_RES f2 = gra::resonance::Read("RES/compact_f2.json", rng, gra::ReggeProductionModel::GP);

  CHECK(f2.spin_basis == "a_Jz");
  CHECK(f2.UsesCoherentSpinBasis());
  CHECK_FALSE(f2.UsesDensitySpinBasis());
  REQUIRE(f2.a_Jz.size() == 5);
  CHECK(std::abs(f2.a_Jz[4] - f2.a_Jz[0]) < 1e-12);
  CHECK(std::abs(f2.a_Jz[3] + f2.a_Jz[1]) < 1e-12);
  CHECK(std::abs(f2.a_Jz[2]) == Approx(std::sqrt(0.5)).margin(1e-12));

  const auto jz1_steered = gra::resonance::CoherentAJzFromDiagonalWeights(
      f2, {0.0, 0.5, 0.0, 0.5, 0.0});
  REQUIRE(jz1_steered.size() == 5);
  CHECK(std::abs(jz1_steered[1]) == Approx(std::sqrt(0.5)).margin(1e-12));
  CHECK(std::abs(jz1_steered[1] / std::abs(jz1_steered[1]) -
                 f2.a_Jz[1] / std::abs(f2.a_Jz[1])) < 1e-12);
  CHECK(std::abs(jz1_steered[3] + jz1_steered[1]) < 1e-12);

  const auto jz2_steered = gra::resonance::CoherentAJzFromDiagonalWeights(
      f2, {0.5, 0.0, 0.0, 0.0, 0.5});
  REQUIRE(jz2_steered.size() == 5);
  CHECK(std::abs(jz2_steered[0]) == Approx(std::sqrt(0.5)).margin(1e-12));
  CHECK(std::abs(jz2_steered[4] - jz2_steered[0]) < 1e-12);

  double norm2 = 0.0;
  for (const auto &amp : f2.a_Jz) {
    norm2 += std::norm(amp);
  }
  CHECK(norm2 == Approx(1.0).margin(1e-12));

  const std::string card_data = gra::aux::GetInputData(card.string());

  // A common input scale cannot change the normalized spin density
  auto near_unit = nlohmann::json::parse(card_data);
  auto &near_unit_rows = near_unit["PARAM_RES"]["MODELS"]["MP"]["[991,991]"]["polarization"]["a_Jz"];
  for (auto &row : near_unit_rows) { row[1] = row[1].get<double>() * std::sqrt(1.0 + 4.0e-7); }
  {
    std::ofstream output(tune_dir / "RES/near_unit_coherent.json");
    output << near_unit.dump();
  }
  const auto near_unit_res = gra::resonance::Read("RES/near_unit_coherent.json", rng, gra::ReggeProductionModel::GP);
  CHECK(gra::SquaredNorm(near_unit_res.a_Jz) == Approx(1.0).epsilon(1e-13));
  CHECK((near_unit_res.rho - f2.rho).FrobNorm2() < 1e-24);

  // Write one malformed derived card and require strict parser rejection
  const auto require_invalid = [&](const nlohmann::json &data,
                                   const std::string &filename) {
    CAPTURE(filename);
    std::ofstream invalid_out(tune_dir / "RES" / filename);
    invalid_out << data.dump();
    invalid_out.close();
    REQUIRE_THROWS_AS(gra::resonance::Read("RES/" + filename, rng, gra::ReggeProductionModel::GP),
                      std::invalid_argument);
  };

  // Require every model to declare both production form factors
  for (const std::string model : {"MP", "XP", "GP", "TP"}) {
    for (const std::string field : {"FF_transfer", "FF_prod"}) {
      auto missing = nlohmann::json::parse(card_data);
      missing["PARAM_RES"]["MODELS"][model].erase(field);
      require_invalid(missing, "missing_" + model + "_" + field + ".json");
    }
  }

  // Normalizing tiny coherent amplitudes must preserve the parity constraint
  auto scaled = nlohmann::json::parse(card_data);
  auto &scaled_rows = scaled["PARAM_RES"]["MODELS"]["MP"]["[991,991]"]["polarization"]["a_Jz"];
  for (auto &row : scaled_rows) { row[1] = row[1].get<double>() * 1e-12; }
  {
    std::ofstream output(tune_dir / "RES/scaled_coherent.json");
    output << scaled.dump();
  }
  const auto scaled_res = gra::resonance::Read("RES/scaled_coherent.json", rng, gra::ReggeProductionModel::GP);
  CHECK((scaled_res.rho - f2.rho).FrobNorm2() < 1e-20);
  auto tiny_odd = scaled;
  tiny_odd["PARAM_RES"]["P"] = -1;
  {
    std::ofstream output(tune_dir / "RES/tiny_odd_jz0.json");
    output << tiny_odd.dump();
  }
  const auto tiny_odd_res = gra::resonance::Read("RES/tiny_odd_jz0.json", rng, gra::ReggeProductionModel::GP);
  CHECK((tiny_odd_res.rho - f2.rho).FrobNorm2() < 1e-20);
  scaled_rows.push_back({2, 0.4e-12, 0.15 + gra::math::PI});
  require_invalid(scaled, "tiny_wrong_parity_pair.json");

  // Reject fractional and oversized quantum numbers before integer conversion
  for (const std::string field : {"PDG", "spinX2", "P", "C"}) {
    for (const auto &value : {nlohmann::json(1.5), nlohmann::json(true),
                             nlohmann::json(4294967300ULL), nlohmann::json(-4294967295LL)}) {
      CAPTURE(field, value);
      auto invalid = nlohmann::json::parse(card_data);
      invalid["PARAM_RES"][field] = value;
      require_invalid(invalid, "invalid_quantum_number.json");
    }
  }
  for (const std::string model : {"GP", "MP", "XP", "TP"}) {
    for (const std::string field : {"mass", "width"}) {
      for (const auto &value : {nlohmann::json(true), nlohmann::json(-1.0), nlohmann::json("invalid")}) {
        auto invalid = nlohmann::json::parse(card_data);
        invalid["PARAM_RES"]["MODELS"][model][field] = value;
        require_invalid(invalid, "invalid_mass_width.json");
      }
    }
  }
  for (const auto &value : {nlohmann::json(4294967296ULL), nlohmann::json(2147483647), nlohmann::json(-2147483648LL)}) {
    auto invalid = nlohmann::json::parse(card_data);
    invalid["PARAM_RES"]["MODELS"]["MP"]["[991,991]"]["polarization"]["a_Jz"][0][0] = value;
    require_invalid(invalid, "invalid_jz.json");
  }

  // Reject missing scales independently of spin steering
  for (const std::string mode : {"a_Jz", "none", "rho"}) {
    auto missing = nlohmann::json::parse(card_data);
    auto& block = missing["PARAM_RES"]["MODELS"]["MP"]["[991,991]"];
    block["polarization"]["mode"] = mode;
    block.erase("Lambda");
    require_invalid(missing, "missing_mp_lambda.json");
  }
  // Reject removed basis names and nested production dynamics
  for (const std::string basis : {"polarization", "fusion_ls", "fusion_helicity", "fusion_auto_min_L", "ls"}) {
    auto invalid = nlohmann::json::parse(card_data);
    invalid["PARAM_RES"]["MODELS"]["MP"]["[991,991]"]["basis"] = basis;
    require_invalid(invalid, "invalid_mp_basis.json");
  }
  auto nested = nlohmann::json::parse(card_data);
  nested["PARAM_RES"]["MODELS"]["MP"]["[991,991]"]["polarization"]["dynamics"] = "auto_min_L";
  require_invalid(nested, "nested_mp_dynamics.json");

  // Cache one canonical model phase for every resonance implementation
  nlohmann::json phase_json = nlohmann::json::parse(card_data);
  phase_json["PARAM_RES"]["MODELS"]["MP"]["phi"] = 0.21;
  phase_json["PARAM_RES"]["MODELS"]["XP"]["phi"] = -0.32;
  phase_json["PARAM_RES"]["MODELS"]["GP"]["phi"] = 0.43;
  phase_json["PARAM_RES"]["MODELS"]["TP"]["phi"] = -0.54;
  const std::filesystem::path phase_card =
      tune_dir / "RES" / "compact_f2_phase.json";
  std::ofstream phase_out(phase_card);
  phase_out << phase_json.dump();
  phase_out.close();
  const gra::PARAM_RES phased =
      gra::resonance::Read("RES/compact_f2_phase.json", rng, gra::ReggeProductionModel::GP);
  CHECK(phased.MP.phi == Approx(0.21));
  CHECK(phased.XP.phi == Approx(-0.32));
  CHECK(phased.GP.phi == Approx(0.43));
  CHECK(phased.TP.phi == Approx(-0.54));

  // Read independent typed channel vectors for every Regge resonance model
  nlohmann::json multichannel = nlohmann::json::parse(card_data);
  auto &models = multichannel["PARAM_RES"]["MODELS"];

  models["MP"]["[991,993]"] = models["MP"]["[991,991]"];
  models["XP"]["[991,993]"] = models["XP"]["[991,991]"];
  models["GP"]["[990,9910]"] = models["GP"]["[990,990]"];
  const std::filesystem::path multichannel_card = tune_dir / "RES" / "compact_f2_multichannel.json";
  std::ofstream multichannel_out(multichannel_card);
  multichannel_out << multichannel.dump();
  multichannel_out.close();
  const gra::PARAM_RES channel_res = gra::resonance::Read("RES/compact_f2_multichannel.json", rng, gra::ReggeProductionModel::GP);
  REQUIRE(channel_res.MP.channels.size() == 2);
  REQUIRE(channel_res.XP.channels.size() == 2);
  REQUIRE(channel_res.GP.channels.size() == 2);
  CHECK(channel_res.MP.channels[1].exchange == std::array<int, 2>{991, 993});
  CHECK(channel_res.XP.channels[1].exchange == std::array<int, 2>{991, 993});
  CHECK(channel_res.GP.channels[1].exchange == std::array<int, 2>{990, 9910});

  // Accept the parameter free disabled model override
  nlohmann::json disabled_mass = nlohmann::json::parse(card_data);
  disabled_mass["PARAM_RES"]["MODELS"]["MP"]["FF_prod"] = {{"type", "none"}};
  const std::filesystem::path disabled_mass_card =
      tune_dir / "RES" / "compact_f2_disabled_mass.json";
  std::ofstream disabled_mass_out(disabled_mass_card);
  disabled_mass_out << disabled_mass.dump();
  disabled_mass_out.close();
  CHECK(gra::resonance::Read("RES/compact_f2_disabled_mass.json", rng, gra::ReggeProductionModel::GP)
            .MP.ff_prod.type == gra::regge::FFType::None);

  // Preserve zero GP rows for validation before runtime pruning
  nlohmann::json zero_gp_ls = nlohmann::json::parse(card_data);
  zero_gp_ls["PARAM_RES"]["MODELS"]["GP"]["[990,990]"]["g_ls"].push_back(
      {2, 2, 0.0, 0.0});
  const std::filesystem::path zero_gp_ls_card =
      tune_dir / "RES" / "compact_f2_zero_gp_ls.json";
  std::ofstream zero_gp_ls_out(zero_gp_ls_card);
  zero_gp_ls_out << zero_gp_ls.dump();
  zero_gp_ls_out.close();
  const gra::PARAM_RES zero_gp_ls_res =
      gra::resonance::Read("RES/compact_f2_zero_gp_ls.json", rng, gra::ReggeProductionModel::GP);
  REQUIRE(zero_gp_ls_res.GP.channels.front().g_ls.Size() == 2);
  REQUIRE(zero_gp_ls_res.GP.channels.front().g_ls.Contains(2, 4));
  CHECK(std::abs(zero_gp_ls_res.GP.channels.front().g_ls.At(2, 4)) ==
        Approx(0.0));

  // Preserve a relative minus sign at the phase boundary
  nlohmann::json signed_gp_ls = nlohmann::json::parse(card_data);
  signed_gp_ls["PARAM_RES"]["MODELS"]["GP"]["[990,990]"]["g_ls"].push_back(
      {2, 2, 0.2, -gra::math::PI});
  const std::filesystem::path signed_gp_ls_card =
      tune_dir / "RES" / "compact_f2_signed_gp_ls.json";
  std::ofstream signed_gp_ls_out(signed_gp_ls_card);
  signed_gp_ls_out << signed_gp_ls.dump();
  signed_gp_ls_out.close();
  const gra::PARAM_RES signed_gp_ls_res =
      gra::resonance::Read("RES/compact_f2_signed_gp_ls.json", rng, gra::ReggeProductionModel::GP);
  const gra::spin::LSTerm *signed_term =
      signed_gp_ls_res.GP.channels.front().g_ls.Find(2, 4);
  REQUIRE(signed_term != nullptr);
  CHECK(signed_term->coefficient.real() == Approx(-0.2).margin(1.0e-15));
  CHECK(signed_term->coefficient.imag() == Approx(0.0).margin(1.0e-15));

  nlohmann::json noncanonical_gp_ls = signed_gp_ls;
  noncanonical_gp_ls["PARAM_RES"]["MODELS"]["GP"]["[990,990]"]["g_ls"][1][3] =
      gra::math::PI;
  require_invalid(noncanonical_gp_ls,
                  "compact_f2_noncanonical_gp_ls.json");

  nlohmann::json zero_gp_helicity = nlohmann::json::parse(card_data);
  auto &zero_gp_helicity_block =
      zero_gp_helicity["PARAM_RES"]["MODELS"]["GP"]["[990,990]"];
  zero_gp_helicity_block["basis"] = "helicity";
  zero_gp_helicity_block.erase("g_ls");
  zero_gp_helicity_block.erase("Lambda");
  zero_gp_helicity_block["helicity"] = {{0, 0, 0.35, 0.0}, {1, 1, 0.0, 0.0}};
  const std::filesystem::path zero_gp_helicity_card =
      tune_dir / "RES" / "compact_f2_zero_gp_helicity.json";
  std::ofstream zero_gp_helicity_out(zero_gp_helicity_card);
  zero_gp_helicity_out << zero_gp_helicity.dump();
  zero_gp_helicity_out.close();
  const gra::PARAM_RES zero_gp_helicity_res =
      gra::resonance::Read("RES/compact_f2_zero_gp_helicity.json", rng, gra::ReggeProductionModel::GP);
  const auto &zero_helicity = zero_gp_helicity_res.GP.channels.front();
  REQUIRE(zero_helicity.g_helicity.size() == 2);
  REQUIRE(zero_helicity.helicity.size() == 2);
  CHECK(std::abs(zero_helicity.g_helicity[1]) == Approx(0.0));

  // Reject missing, nonnumeric, and noncanonical model phases
  for (const std::string model : {"MP", "XP", "GP", "TP"}) {
    nlohmann::json missing = nlohmann::json::parse(card_data);
    missing["PARAM_RES"]["MODELS"][model].erase("phi");
    require_invalid(missing, "compact_f2_missing_" + model + "_phi.json");
  }
  std::size_t bad_phase_index = 0;
  for (const nlohmann::json &value :
       {nlohmann::json(gra::math::PI), nlohmann::json(-gra::math::PI - 1.0e-9),
        nlohmann::json("zero"), nlohmann::json(nullptr)}) {
    nlohmann::json invalid = nlohmann::json::parse(card_data);
    invalid["PARAM_RES"]["MODELS"]["GP"]["phi"] = value;
    require_invalid(invalid, "compact_f2_bad_phi_" +
                                 std::to_string(bad_phase_index++) + ".json");
  }

  // Preserve an arbitrary common LS phase and reject a zero-row phase
  nlohmann::json row_phase = nlohmann::json::parse(card_data);
  row_phase["PARAM_RES"]["MODELS"]["GP"]["[990,990]"]["g_ls"][0][3] = 0.2;
  const std::filesystem::path row_phase_card =
      tune_dir / "RES" / "compact_f2_free_reference_phase.json";
  std::ofstream row_phase_out(row_phase_card);
  row_phase_out << row_phase.dump();
  row_phase_out.close();
  const auto free_phase =
      gra::resonance::Read("RES/compact_f2_free_reference_phase.json", rng, gra::ReggeProductionModel::GP);
  const auto *reference = free_phase.GP.channels.front().g_ls.Find(0, 4);
  REQUIRE(reference != nullptr);
  RequireComplexNear(reference->coefficient, std::polar(0.35, 0.2), 1.0e-14);
  nlohmann::json zero_phase = nlohmann::json::parse(card_data);
  zero_phase["PARAM_RES"]["MODELS"]["GP"]["[990,990]"]["g_ls"].push_back(
      {2, 2, 0.0, 0.2});
  require_invalid(zero_phase, "compact_f2_bad_zero_phase.json");
  nlohmann::json tensor_sign = nlohmann::json::parse(card_data);
  tensor_sign["PARAM_RES"]["MODELS"]["TP"]["[995,995]"]["g_tensor"][0] = -1.0;
  require_invalid(tensor_sign, "compact_f2_bad_tensor_reference.json");
  nlohmann::json zero_model = nlohmann::json::parse(card_data);
  zero_model["PARAM_RES"]["MODELS"]["MP"]["[991,991]"]["g"] = {0.0, 0.0};
  // A zero scalar coupling is valid input and removes this channel from sampling
  {
    std::ofstream output(tune_dir / "RES/compact_f2_zero_mp_model.json");
    output << zero_model;
  }
  const auto inactive = gra::resonance::Read("RES/compact_f2_zero_mp_model.json", rng, gra::ReggeProductionModel::GP);
  CHECK_FALSE(inactive.MP.channels.front().Active(0.0));

  // Require strict string routing for each resonance denominator mode
  const std::array<std::pair<const char *, gra::BreitWigner>, 3> bw_modes = {{
      {"fixed-width", gra::BreitWigner::FixedWidth},
      {"kinematic-width", gra::BreitWigner::KinematicWidth},
      {"running-width", gra::BreitWigner::RunningWidth},
  }};
  for (const auto &[mode, expected] : bw_modes) {
    nlohmann::json mode_json = nlohmann::json::parse(card_data);
    mode_json["PARAM_RES"]["MODELS"]["GP"]["BW"] = mode;
    const std::string filename = std::string("compact_f2_") + mode + ".json";
    std::ofstream mode_out(tune_dir / "RES" / filename);
    mode_out << mode_json;
    mode_out.close();
    CHECK(gra::resonance::Read("RES/" + filename, rng, gra::ReggeProductionModel::GP).BW == expected);
  }

  // Reject obsolete numeric and unknown resonance denominator modes
  for (const nlohmann::json &mode :
       {nlohmann::json(2), nlohmann::json("constant")}) {
    nlohmann::json invalid_bw = nlohmann::json::parse(card_data);
    invalid_bw["PARAM_RES"]["MODELS"]["GP"]["BW"] = mode;
    const std::string filename = mode.is_number()
                                     ? "compact_f2_numeric_bw.json"
                                     : "compact_f2_unknown_bw.json";
    std::ofstream invalid_out(tune_dir / "RES" / filename);
    invalid_out << invalid_bw;
    invalid_out.close();
    REQUIRE_THROWS_AS(gra::resonance::Read("RES/" + filename, rng, gra::ReggeProductionModel::GP),
                      std::invalid_argument);
  }

  // Reject BW steering in the tensor model
  auto tensor_bw = nlohmann::json::parse(card_data);
  tensor_bw["PARAM_RES"]["MODELS"]["TP"]["BW"] = "fixed-width";
  require_invalid(tensor_bw, "tensor_bw.json");

  // Reject automatic GP couplings because analytic vertices must be explicit
  nlohmann::json automatic_gp_json = nlohmann::json::parse(card_data);
  auto &automatic_gp =
      automatic_gp_json["PARAM_RES"]["MODELS"]["GP"]["[990,990]"];
  automatic_gp["basis"] = "auto_min_L";
  automatic_gp.erase("g_ls");
  automatic_gp["g"] = {0.35, 0.0};
  const std::filesystem::path automatic_gp_card =
      tune_dir / "RES" / "compact_f2_automatic_gp.json";
  std::ofstream automatic_gp_out(automatic_gp_card);
  automatic_gp_out << automatic_gp_json.dump();
  automatic_gp_out.close();
  REQUIRE_THROWS_AS(
      gra::resonance::Read("RES/compact_f2_automatic_gp.json", rng, gra::ReggeProductionModel::GP),
      std::invalid_argument);

  // Fixed-spin XP vertices retain automatic coupling generation
  nlohmann::json automatic_xp_json = nlohmann::json::parse(card_data);
  auto &automatic_xp =
      automatic_xp_json["PARAM_RES"]["MODELS"]["XP"]["[991,991]"];
  automatic_xp["basis"] = "auto_min_S";
  automatic_xp.erase("g_ls");
  automatic_xp["g"] = {0.35, 0.0};
  const std::filesystem::path automatic_xp_card =
      tune_dir / "RES" / "compact_f2_automatic_xp.json";
  std::ofstream automatic_xp_out(automatic_xp_card);
  automatic_xp_out << automatic_xp_json.dump();
  automatic_xp_out.close();
  const gra::PARAM_RES automatic =
      gra::resonance::Read("RES/compact_f2_automatic_xp.json", rng, gra::ReggeProductionModel::GP);
  CHECK(automatic.GP.channels.front().basis == gra::ReggeVertexBasis::LS);
  REQUIRE(automatic.XP.channels.size() == 1);
  CHECK(automatic.XP.channels.front().basis == gra::ReggeVertexBasis::AutoMinS);

  // Accept valid dormant spin steering when mode none is active
  nlohmann::json none_json = nlohmann::json::parse(card_data);
  none_json["PARAM_RES"]["MODELS"]["MP"]["[991,991]"]["polarization"]["mode"] =
      "none";
  none_json["PARAM_RES"]["MODELS"]["MP"]["[991,991]"]["basis"] = "auto_min_L";
  none_json["PARAM_RES"]["MODELS"]["MP"]["[991,991]"]["Lambda"] = 1.0;
  const std::filesystem::path none_card =
      tune_dir / "RES" / "compact_f2_none_with_steering.json";
  std::ofstream none_out(none_card);
  none_out << none_json.dump();
  none_out.close();
  const gra::PARAM_RES f2_none =
      gra::resonance::Read("RES/compact_f2_none_with_steering.json", rng, gra::ReggeProductionModel::GP);
  CHECK(f2_none.spin_basis == "none");
  CHECK(f2_none.a_Jz.empty());
  CHECK(f2_none.rho.isEmpty());

  // Reject obsolete resonance model identifiers instead of accepting aliases
  nlohmann::json obsolete_model_json = nlohmann::json::parse(card_data);
  obsolete_model_json["PARAM_RES"]["MODELS"]["MRES"] =
      obsolete_model_json["PARAM_RES"]["MODELS"]["MP"];
  const std::filesystem::path obsolete_model_card =
      tune_dir / "RES" / "compact_f2_obsolete_model.json";
  std::ofstream obsolete_model_out(obsolete_model_card);
  obsolete_model_out << obsolete_model_json.dump();
  obsolete_model_out.close();
  REQUIRE_THROWS(
      gra::resonance::Read("RES/compact_f2_obsolete_model.json", rng, gra::ReggeProductionModel::GP));

  // Reject obsolete or misspelled XP fields instead of silently ignoring them
  nlohmann::json xp_unknown_json = nlohmann::json::parse(card_data);
  xp_unknown_json["PARAM_RES"]["MODELS"]["XP"]["g_channel"] = 1.0;
  const std::filesystem::path xp_unknown_card =
      tune_dir / "RES" / "compact_f2_bad_xp_field.json";
  std::ofstream xp_unknown_out(xp_unknown_card);
  xp_unknown_out << xp_unknown_json.dump();
  xp_unknown_out.close();
  REQUIRE_THROWS(gra::resonance::Read("RES/compact_f2_bad_xp_field.json", rng, gra::ReggeProductionModel::GP));

  nlohmann::json xp_block_unknown_json = nlohmann::json::parse(card_data);
  xp_block_unknown_json["PARAM_RES"]["MODELS"]["XP"]["[991,991]"]["alpha_ls"] =
      {{2, 0, 1.0, 0.0}};
  const std::filesystem::path xp_block_unknown_card =
      tune_dir / "RES" / "compact_f2_bad_xp_block_field.json";
  std::ofstream xp_block_unknown_out(xp_block_unknown_card);
  xp_block_unknown_out << xp_block_unknown_json.dump();
  xp_block_unknown_out.close();
  REQUIRE_THROWS(
      gra::resonance::Read("RES/compact_f2_bad_xp_block_field.json", rng, gra::ReggeProductionModel::GP));

  // Reject unused MP and GP fields so every normalization parameter has one
  // owner
  nlohmann::json mp_unknown_json = nlohmann::json::parse(card_data);
  mp_unknown_json["PARAM_RES"]["MODELS"]["MP"]["g_channel"] = 1.0;
  const std::filesystem::path mp_unknown_card =
      tune_dir / "RES" / "compact_f2_bad_mp_field.json";
  std::ofstream mp_unknown_out(mp_unknown_card);
  mp_unknown_out << mp_unknown_json.dump();
  mp_unknown_out.close();
  REQUIRE_THROWS(gra::resonance::Read("RES/compact_f2_bad_mp_field.json", rng, gra::ReggeProductionModel::GP));

  nlohmann::json gp_block_unknown_json = nlohmann::json::parse(card_data);
  gp_block_unknown_json["PARAM_RES"]["MODELS"]["GP"]["[990,990]"]
                       ["g_channel"] = 1.0;
  const std::filesystem::path gp_block_unknown_card =
      tune_dir / "RES" / "compact_f2_bad_gp_block_field.json";
  std::ofstream gp_block_unknown_out(gp_block_unknown_card);
  gp_block_unknown_out << gp_block_unknown_json.dump();
  gp_block_unknown_out.close();
  REQUIRE_THROWS(
      gra::resonance::Read("RES/compact_f2_bad_gp_block_field.json", rng, gra::ReggeProductionModel::GP));

  nlohmann::json gp_basis_mismatch_json = nlohmann::json::parse(card_data);
  gp_basis_mismatch_json["PARAM_RES"]["MODELS"]["GP"]["[990,990]"]
                        ["helicity"] = nlohmann::json::array(
                            {nlohmann::json::array({0, 0, 1.0, 0.0})});
  const std::filesystem::path gp_basis_mismatch_card =
      tune_dir / "RES" / "compact_f2_bad_gp_basis_fields.json";
  std::ofstream gp_basis_mismatch_out(gp_basis_mismatch_card);
  gp_basis_mismatch_out << gp_basis_mismatch_json.dump();
  gp_basis_mismatch_out.close();
  REQUIRE_THROWS(
      gra::resonance::Read("RES/compact_f2_bad_gp_basis_fields.json", rng, gra::ReggeProductionModel::GP));

  nlohmann::json parity_bad_json = nlohmann::json::parse(card_data);
  parity_bad_json["PARAM_RES"]["MODELS"]["MP"]["[991,991]"]["polarization"]["a_Jz"] = nlohmann::json::array(
      {nlohmann::json::array({-2, 0.4, 0.15}), nlohmann::json::array({-1, 0.3, -0.2}),
       nlohmann::json::array({0, 0.7071067811865476, 0.0}), nlohmann::json::array({1, 0.3, -0.2}),
       nlohmann::json::array({2, 0.4, 0.15})});

  const std::filesystem::path parity_bad_card =
      tune_dir / "RES" / "compact_f2_bad_ajz_parity.json";
  std::ofstream parity_bad_out(parity_bad_card);
  parity_bad_out << parity_bad_json.dump();
  parity_bad_out.close();

  REQUIRE_THROWS(
      gra::resonance::Read("RES/compact_f2_bad_ajz_parity.json", rng, gra::ReggeProductionModel::GP));

  nlohmann::json eta_steered_json = nlohmann::json::parse(card_data);
  eta_steered_json["PARAM_RES"]["MODELS"]["MP"]["[991,991]"]["eta_P"] = 1;

  const std::filesystem::path eta_steered_card =
      tune_dir / "RES" / "compact_f2_steered_eta_p.json";
  std::ofstream eta_steered_out(eta_steered_card);
  eta_steered_out << eta_steered_json.dump();
  eta_steered_out.close();

  REQUIRE_THROWS(
      gra::resonance::Read("RES/compact_f2_steered_eta_p.json", rng, gra::ReggeProductionModel::GP));

  nlohmann::json negative_parity_json = nlohmann::json::parse(card_data);
  negative_parity_json["PARAM_RES"]["P"] = -1;

  const std::filesystem::path negative_parity_card =
      tune_dir / "RES" / "compact_f2_negative_parity.json";
  std::ofstream negative_parity_out(negative_parity_card);
  negative_parity_out << negative_parity_json.dump();
  negative_parity_out.close();

  const gra::PARAM_RES f2_negative_parity =
      gra::resonance::Read("RES/compact_f2_negative_parity.json", rng, gra::ReggeProductionModel::GP);
  REQUIRE(f2_negative_parity.a_Jz.size() == 5);
  CHECK(std::abs(f2_negative_parity.a_Jz[4] - f2_negative_parity.a_Jz[0]) < 1e-12);
  CHECK(std::abs(f2_negative_parity.a_Jz[3] + f2_negative_parity.a_Jz[1]) <
        1e-12);
  CHECK(std::abs(f2_negative_parity.a_Jz[2]) == Approx(std::sqrt(0.5)).margin(1e-12));
  CHECK_NOTHROW(gra::resonance::CoherentAJzFromDiagonalWeights(
      f2_negative_parity, {0.25, 0.25, 0.0, 0.25, 0.25}));
  CHECK_NOTHROW(gra::resonance::CoherentAJzFromDiagonalWeights(f2_negative_parity, {0.125, 0.125, 0.5, 0.125, 0.125}));
  auto wrong_spin_dimension = f2_negative_parity;
  wrong_spin_dimension.a_Jz.clear();
  REQUIRE_THROWS_AS(gra::resonance::CoherentAJzFromDiagonalWeights(wrong_spin_dimension, {}),
                    std::invalid_argument);

  nlohmann::json parity_free_json = nlohmann::json::parse(card_data);
  parity_free_json["PARAM_RES"]["MODELS"]["MP"]["[991,991]"]["CP"] =
      nlohmann::json::array({true, false});

  const std::filesystem::path parity_free_card =
      tune_dir / "RES" / "compact_f2_parity_free.json";
  std::ofstream parity_free_out(parity_free_card);
  parity_free_out << parity_free_json.dump();
  parity_free_out.close();

  const gra::PARAM_RES f2_parity_free =
      gra::resonance::Read("RES/compact_f2_parity_free.json", rng, gra::ReggeProductionModel::GP);
  CHECK_FALSE(f2_parity_free.MP.channels.front().P_symmetry);
  CHECK(std::abs(f2_parity_free.a_Jz[4]) < 1e-12);
  CHECK(std::abs(f2_parity_free.a_Jz[3]) < 1e-12);

  nlohmann::json rho_json = nlohmann::json::parse(card_data);
  rho_json["PARAM_RES"]["MODELS"]["MP"]["[991,991]"]["polarization"]["mode"] =
      "rho";

  const std::filesystem::path rho_card =
      tune_dir / "RES" / "compact_f2_rho.json";
  std::ofstream rho_out(rho_card);
  rho_out << rho_json.dump();
  rho_out.close();

  const gra::PARAM_RES f2_rho =
      gra::resonance::Read("RES/compact_f2_rho.json", rng, gra::ReggeProductionModel::GP);
  CHECK(f2_rho.spin_basis == "rho");
  CHECK_FALSE(f2_rho.UsesCoherentSpinBasis());
  CHECK(f2_rho.UsesDensitySpinBasis());

  // A coherent projector remains valid as a density with nonzero complex coherences
  auto coherent_rho = rho_json;
  auto &coherent_pol = coherent_rho["PARAM_RES"]["MODELS"]["MP"]["[991,991]"]["polarization"];
  for (std::size_t i = 0; i < f2.rho.size_row(); ++i) {
    for (std::size_t j = 0; j < f2.rho.size_col(); ++j) {
      coherent_pol["rho_mag"][i][j] = std::abs(f2.rho[i][j]);
      coherent_pol["rho_phase"][i][j] = std::arg(f2.rho[i][j]);
    }
  }
  {
    std::ofstream out(tune_dir / "RES" / "coherent_rho.json");
    out << coherent_rho;
  }
  const auto coherent_density = gra::resonance::Read("RES/coherent_rho.json", rng, gra::ReggeProductionModel::GP);
  CHECK((coherent_density.rho - f2.rho).FrobNorm2() < 1.0e-20);

  // A polarized density must obey the same reflection as random_rho when P is active
  auto asymmetric_rho = rho_json;
  auto &asymmetric_pol = asymmetric_rho["PARAM_RES"]["MODELS"]["MP"]["[991,991]"]["polarization"];
  asymmetric_pol["rho_mag"] = {{1,0,0,0,0}, {0,0,0,0,0}, {0,0,0,0,0}, {0,0,0,0,0}, {0,0,0,0,0}};
  require_invalid(asymmetric_rho, "asymmetric_rho.json");
  auto wrong_reflection = rho_json;
  wrong_reflection["PARAM_RES"]["MODELS"]["MP"]["[991,991]"]["polarization"]["rho_mag"] =
      {{0.25,0.25,0,0.25,0.25}, {0.25,0.25,0,0.25,0.25}, {0,0,0,0,0},
       {0.25,0.25,0,0.25,0.25}, {0.25,0.25,0,0.25,0.25}};
  wrong_reflection["PARAM_RES"]["MODELS"]["MP"]["[991,991]"]["polarization"]["rho_phase"] =
      {{0,0,0,0,0}, {0,0,0,0,0}, {0,0,0,0,0}, {0,0,0,0,0}, {0,0,0,0,0}};
  require_invalid(wrong_reflection, "wrong_rho_reflection.json");
  asymmetric_rho["PARAM_RES"]["MODELS"]["MP"]["[991,991]"]["CP"][1] = false;
  {
    std::ofstream out(tune_dir / "RES" / "asymmetric_rho_allowed.json");
    out << asymmetric_rho;
  }
  const auto polarized = gra::resonance::Read("RES/asymmetric_rho_allowed.json", rng, gra::ReggeProductionModel::GP);
  CHECK(polarized.rho[0][0].real() == Approx(1.0));
  CHECK(std::abs(polarized.rho[4][4]) < 1.0e-12);

  // Reject extra density rows and columns instead of silently truncating them
  for (const std::string field : {"rho_mag", "rho_phase"}) {
    auto nonnumeric = rho_json;
    nonnumeric["PARAM_RES"]["MODELS"]["MP"]["[991,991]"]["polarization"][field][0][0] = true;
    require_invalid(nonnumeric, "nonnumeric_rho.json");
    auto oversized = rho_json;
    auto &rows = oversized["PARAM_RES"]["MODELS"]["MP"]["[991,991]"]["polarization"][field];
    for (auto &row : rows) { row.push_back(0.0); }
    rows.push_back({0,0,0,0,0,0});
    require_invalid(oversized, "oversized_rho.json");
  }

  nlohmann::json random_json = rho_json;
  random_json["PARAM_RES"]["MODELS"]["MP"]["[991,991]"]["polarization"]
             ["random_rho"] = true;
  random_json["PARAM_RES"]["MODELS"]["MP"]["[991,991]"]["CP"] =
      nlohmann::json::array({true, false});

  const std::filesystem::path random_card =
      tune_dir / "RES" / "compact_f2_random_rho.json";
  std::ofstream random_out(random_card);
  random_out << random_json.dump();
  random_out.close();

  MRandom parser_rng;
  MRandom direct_rng;
  parser_rng.SetSeed(17);
  direct_rng.SetSeed(17);
  const gra::PARAM_RES f2_random =
      gra::resonance::Read("RES/compact_f2_random_rho.json", parser_rng, gra::ReggeProductionModel::GP);
  const auto expected_rho = gra::spin::RandomRho(4, false, direct_rng);
  CHECK_FALSE(f2_random.MP.channels.front().P_symmetry);
  REQUIRE(f2_random.rho.size_row() == expected_rho.size_row());
  REQUIRE(f2_random.rho.size_col() == expected_rho.size_col());
  for (std::size_t i = 0; i < f2_random.rho.size_row(); ++i) {
    for (std::size_t j = 0; j < f2_random.rho.size_col(); ++j) {
      CHECK(std::abs(f2_random.rho[i][j] - expected_rho[i][j]) < 1e-12);
    }
  }
}

TEST_CASE("proton electromagnetic form factors preserve Sachs identities",
          "[gra::form][electromagnetic]") {
  const double infinity = std::numeric_limits<double>::infinity();
  REQUIRE(gra::math::IsZero(form::G_E(infinity)));
  REQUIRE_FALSE(std::isfinite(form::G_M_KELLY(infinity)));
  REQUIRE(gra::math::IsZero(form::F2xQ2(0.1, infinity)));

  gra::LORENTZSCALAR lts;
  lts.PDG = LoadedPDGTable();
  const auto model_tune = gra::MModelTune::Load(modelfile);
  MTensorPomeron tensor(lts, model_tune,
                        gra::MTensorPomeron::ProcessDefinitionFor(
                            gra::MTensorPomeronMode::Generic));

  for (const std::string mode : {"DIPOLE", "KELLY"}) {
    CAPTURE(mode);
    gra::form::ParamStore structure;
    structure.EM = mode;
    REQUIRE(form::G_E(0.0, structure) == Approx(1.0).margin(1e-15));
    REQUIRE(form::G_M(0.0, structure) ==
            Approx(form::mu_ratio()).margin(1e-15));
    REQUIRE(form::F1(0.0, structure) == Approx(1.0).margin(1e-15));
    REQUIRE(form::F2(0.0, structure) ==
            Approx(form::mu_ratio() - 1.0).margin(1e-15));

    for (const double q2 : {0.05, 0.5, 1.25, 5.0}) {
      CAPTURE(q2);
      const double tau = q2 / (4.0 * gra::math::pow2(gra::PDG::mp));
      const double f1 = form::F1(q2, structure);
      const double f2 = form::F2(q2, structure);

      REQUIRE(f1 - tau * f2 ==
              Approx(form::G_E(q2, structure)).margin(2e-15));
      REQUIRE(f1 + f2 ==
              Approx(form::G_M(q2, structure)).margin(2e-15));
      REQUIRE(form::F1(-q2, structure) == Approx(f1).margin(1e-15));
      REQUIRE(form::F2(-q2, structure) == Approx(f2).margin(1e-15));
    }
  }

  const auto &structure = model_tune->Structure();
  for (const double q2 : {0.0, 0.05, 0.5, 1.25, 5.0}) {
    REQUIRE(form::F1(q2, structure) ==
            Approx(tensor.F1_(-q2)).margin(2e-15));
    REQUIRE(form::F2(q2, structure) ==
            Approx(tensor.F2_(-q2)).margin(2e-15));
  }

  const double q2 = 1.25;
  const double dipole_reference = 1.0 / gra::math::pow2(1.0 + q2 / 0.71);
  REQUIRE(form::G_E_DIPOLE(q2) == Approx(dipole_reference).margin(1e-15));
  REQUIRE(form::G_E_DIPOLE(-q2) == Approx(dipole_reference).margin(1e-15));
  REQUIRE(form::G_E_DIPOLE(q2) != Approx(form::G_E_KELLY(q2)).epsilon(1e-6));
}

TEST_CASE("Four-point SW signature factor follows the GRANIITTI convention",
          "[ReggeSW]") {
  using gra::regge::EtaRaw;
  using gra::regge::Rim;
  using gra::regge::Signature;

  const std::array<std::complex<double>, 3> spins = {
      std::complex<double>(0.37, 0.23), std::complex<double>(-1.21, -0.17),
      std::complex<double>(2.44, 0.09)};
  for (const auto J : spins) {
    for (const Signature signature :
         {Signature::Positive, Signature::Negative}) {
      const double tau = gra::regge::Tau(signature);
      for (const Rim rim : {Rim::Lower, Rim::Upper}) {
        CAPTURE(J.real(), J.imag(), tau, rim);
        const double rim_sign = rim == Rim::Lower ? -1.0 : 1.0;
        const std::complex<double> expected =
            -(1.0 +
              tau * std::exp(rim_sign * gra::math::zi * gra::math::PI * J)) /
            std::sin(gra::math::PI * J);
        RequireComplexNear(EtaRaw(J, signature, rim), expected, 3.0e-14);
      }
    }
  }
}

TEST_CASE("Four-point SW rims obey reflection and periodicity", "[ReggeSW]") {
  using gra::regge::EtaRaw;
  using gra::regge::Rim;
  using gra::regge::Signature;

  for (const std::complex<double> J :
       {std::complex<double>(0.31, 0.27), std::complex<double>(-1.17, -0.19)}) {
    for (const Signature signature :
         {Signature::Positive, Signature::Negative}) {
      const int tau = gra::regge::Tau(signature);
      CAPTURE(J.real(), J.imag(), tau);
      RequireComplexNear(EtaRaw(std::conj(J), signature, Rim::Upper),
                         std::conj(EtaRaw(J, signature, Rim::Lower)), 3.0e-14);
      RequireComplexNear(EtaRaw(J + 2.0, signature, Rim::Lower),
                         EtaRaw(J, signature, Rim::Lower), 3.0e-14);
      RequireComplexNear(EtaRaw(J - 2.0, signature, Rim::Upper),
                         EtaRaw(J, signature, Rim::Upper), 3.0e-14);
    }
  }

  for (const double J : {-1.37, -0.42, 0.23, 1.61}) {
    for (const Signature signature :
         {Signature::Positive, Signature::Negative}) {
      const int tau = gra::regge::Tau(signature);
      CAPTURE(J, tau);
      RequireComplexNear(EtaRaw(J, signature, Rim::Upper),
                         std::conj(EtaRaw(J, signature, Rim::Lower)), 1.0e-15);
      REQUIRE(std::imag(EtaRaw(J, signature, Rim::Lower)) ==
              Approx(tau).margin(1.0e-15));
      REQUIRE(std::imag(EtaRaw(J, signature, Rim::Upper)) ==
              Approx(-tau).margin(1.0e-15));
    }
  }
}

TEST_CASE("Four-point SW contours remain finite at large imaginary spin",
          "[ReggeSW]") {
  using gra::regge::EtaRaw;
  using gra::regge::Rim;
  using gra::regge::Signature;

  for (const Signature signature : {Signature::Positive, Signature::Negative}) {
    const double tau = gra::regge::Tau(signature);
    CAPTURE(tau);
    const std::complex<double> upper_half(0.37, 400.0);
    const std::complex<double> lower_half(0.37, -400.0);
    const auto lower_rim_upper_half = EtaRaw(upper_half, signature, Rim::Lower);
    const auto upper_rim_upper_half = EtaRaw(upper_half, signature, Rim::Upper);
    const auto lower_rim_lower_half = EtaRaw(lower_half, signature, Rim::Lower);
    const auto upper_rim_lower_half = EtaRaw(lower_half, signature, Rim::Upper);

    REQUIRE(std::isfinite(std::abs(lower_rim_upper_half)));
    REQUIRE(std::isfinite(std::abs(upper_rim_upper_half)));
    REQUIRE(std::isfinite(std::abs(lower_rim_lower_half)));
    REQUIRE(std::isfinite(std::abs(upper_rim_lower_half)));
    RequireComplexNear(lower_rim_upper_half, 2.0 * gra::math::zi * tau,
                       1.0e-15);
    RequireComplexNear(upper_rim_upper_half, {0.0, 0.0}, 1.0e-15);
    RequireComplexNear(lower_rim_lower_half, {0.0, 0.0}, 1.0e-15);
    RequireComplexNear(upper_rim_lower_half, -2.0 * gra::math::zi * tau,
                       1.0e-15);
  }
}

TEST_CASE("Four-point SW integer signatures have the correct pole algebra",
          "[ReggeSW]") {
  using gra::regge::Classify;
  using gra::regge::EtaRaw;
  using gra::regge::IntegerType;
  using gra::regge::Rim;
  using gra::regge::Signature;

  constexpr double delta = 1.0e-5;
  for (int J = -4; J <= 4; ++J) {
    const Signature signature =
        J % 2 == 0 ? Signature::Positive : Signature::Negative;
    const Signature opposite = signature == Signature::Positive
                                   ? Signature::Negative
                                   : Signature::Positive;
    const double tau = gra::regge::Tau(signature);
    CAPTURE(J, tau);
    REQUIRE(Classify(J, signature) == IntegerType::Allowed);
    REQUIRE(Classify(J, opposite) == IntegerType::Opposite);

    for (const Rim rim : {Rim::Lower, Rim::Upper}) {
      CAPTURE(rim);
      const std::complex<double> residue =
          0.5 * delta *
          (EtaRaw(static_cast<double>(J) + delta, signature, rim) -
           EtaRaw(static_cast<double>(J) - delta, signature, rim));
      RequireComplexNear(residue, {-2.0 * tau / gra::math::PI, 0.0}, 2.0e-10);

      const double opposite_tau = gra::regge::Tau(opposite);
      const std::complex<double> expected = rim == Rim::Lower
                                                ? gra::math::zi * opposite_tau
                                                : -gra::math::zi * opposite_tau;
      RequireComplexNear(EtaRaw(static_cast<double>(J), opposite, rim),
                         expected, 1.0e-15);
      REQUIRE(std::isfinite(
          std::abs(EtaRaw(static_cast<double>(J) + 1.0e-12, opposite, rim))));
    }
  }
}

TEST_CASE("Four-point SW pole uses the positive-s power", "[ReggeSW]") {
  using gra::regge::EtaRaw;
  using gra::regge::Pole;
  using gra::regge::Power;
  using gra::regge::Rim;
  using gra::regge::Signature;

  const std::complex<double> J(1.23, -0.31);
  const double s = 17.0;
  const double s0 = 2.5;
  const std::complex<double> expected = std::exp(J * std::log(s / s0));
  RequireComplexNear(Power(9.0, 1.0, {0.5, 0.0}), {3.0, 0.0}, 1.0e-15);
  RequireComplexNear(Power(s, s0, J), expected, 2.0e-15);
  RequireComplexNear(Pole(s, s0, J, Signature::Positive, Rim::Lower),
                     EtaRaw(J, Signature::Positive, Rim::Lower) * expected,
                     3.0e-14);
  RequireComplexNear(
      std::exp(-gra::math::zi * gra::math::PI * J) * Power(s, s0, J),
      std::exp(J * std::complex<double>(std::log(s / s0), -gra::math::PI)),
      3.0e-14);
  RequireComplexNear(
      std::exp(gra::math::zi * gra::math::PI * J) * Power(s, s0, J),
      std::exp(J * std::complex<double>(std::log(s / s0), gra::math::PI)),
      3.0e-14);
  const auto half_pole =
      Pole(9.0, 1.0, {0.5, 0.0}, Signature::Positive, Rim::Lower);
  RequireComplexNear(
      half_pole, 3.0 * EtaRaw(0.5, Signature::Positive, Rim::Lower), 1.0e-15);
}

TEST_CASE("Four-point SW rejects invalid contour inputs", "[ReggeSW]") {
  using gra::regge::Classify;
  using gra::regge::EtaRaw;
  using gra::regge::Pole;
  using gra::regge::Power;
  using gra::regge::Rim;
  using gra::regge::Signature;
  using gra::regge::Tau;

  const double infinity = std::numeric_limits<double>::infinity();
  const double quiet_nan = std::numeric_limits<double>::quiet_NaN();
  const Signature invalid_signature = static_cast<Signature>(17);

  RequireComplexNear(
      EtaRaw(0.37, Signature::Negative, Rim::Upper),
      EtaRaw(std::complex<double>(0.37, 0.0), Signature::Negative, Rim::Upper),
      1.0e-15);
  RequireComplexNear(
      Pole(17.0, 2.5, {1.23, -0.31}, Signature::Negative, Rim::Upper),
      EtaRaw({1.23, -0.31}, Signature::Negative, Rim::Upper) *
          Power(17.0, 2.5, {1.23, -0.31}),
      3.0e-14);

  REQUIRE_FALSE(std::isfinite(std::norm(EtaRaw(infinity, Signature::Positive, Rim::Lower))));
  REQUIRE_FALSE(std::isfinite(std::norm(EtaRaw(std::complex<double>(0.5, quiet_nan),
                                              Signature::Positive, Rim::Lower))));
  REQUIRE_THROWS_AS(Tau(invalid_signature), std::invalid_argument);
  REQUIRE_THROWS_AS(Classify(2, invalid_signature), std::invalid_argument);
  REQUIRE_THROWS_AS(Power(0.0, 1.0, {0.5, 0.0}), gra::AmplitudeFailure);
  REQUIRE_THROWS_AS(Power(-1.0, 1.0, {0.5, 0.0}), gra::AmplitudeFailure);
  REQUIRE_FALSE(std::isfinite(std::norm(Power(1.0, 1.0, {infinity, 0.0}))));
  REQUIRE_THROWS_AS(Pole(0.0, 1.0, {0.5, 0.0}, Signature::Positive, Rim::Lower),
                    gra::AmplitudeFailure);
}

TEST_CASE("gra::regge::Eta current signature factor", "[ReggeEta]") {
  using gra::regge::Signature;
  constexpr double POLE_EPSILON = 1e-8;

  SECTION("analytic values away from poles") {
    constexpr double IDENTITY_EPSILON = 1e-12;
    struct SignaturePoint {
      double alpha;
      Signature signature;
      double real;
      double imag;
    };
    const std::array<SignaturePoint, 4> points = {
        SignaturePoint{0.5, Signature::Positive, -1.0, 1.0},
        SignaturePoint{0.5, Signature::Negative, -1.0, -1.0},
        SignaturePoint{1.5, Signature::Positive, 1.0, 1.0},
        SignaturePoint{1.5, Signature::Negative, 1.0, -1.0}};

    for (const SignaturePoint &point : points) {
      CAPTURE(point.alpha, gra::regge::Tau(point.signature));
      const std::complex<double> eta =
          gra::regge::Eta(point.alpha, point.signature, IDENTITY_EPSILON);
      REQUIRE(std::real(eta) == Approx(point.real).margin(2e-12));
      REQUIRE(std::imag(eta) == Approx(point.imag).margin(2e-12));
    }
  }

  SECTION("opposite-signature poles have finite analytic limits") {
    const std::complex<double> even_at_odd =
        gra::regge::Eta(1.0, Signature::Positive, POLE_EPSILON);
    const std::complex<double> odd_at_even =
        gra::regge::Eta(0.0, Signature::Negative, POLE_EPSILON);
    const std::complex<double> strongly_regulated_even_at_odd =
        gra::regge::Eta(1.0, Signature::Positive, 1.0);
    const std::complex<double> strongly_regulated_odd_at_even =
        gra::regge::Eta(0.0, Signature::Negative, 1.0);

    REQUIRE(std::real(even_at_odd) == Approx(0.0).margin(1e-14));
    REQUIRE(std::imag(even_at_odd) == Approx(1.0).margin(1e-14));
    REQUIRE(std::real(odd_at_even) == Approx(0.0).margin(1e-14));
    REQUIRE(std::imag(odd_at_even) == Approx(-1.0).margin(1e-14));
    REQUIRE(std::imag(strongly_regulated_even_at_odd) ==
            Approx(1.0).margin(1e-14));
    REQUIRE(std::imag(strongly_regulated_odd_at_even) ==
            Approx(-1.0).margin(1e-14));
  }

  SECTION("signature-allowed poles retain a finite regulator") {
    const std::complex<double> even_pole =
        gra::regge::Eta(0.0, Signature::Positive, POLE_EPSILON);
    const std::complex<double> odd_pole =
        gra::regge::Eta(1.0, Signature::Negative, POLE_EPSILON);
    const std::complex<double> wider_even_pole =
        gra::regge::Eta(0.0, Signature::Positive, 1e-4);

    REQUIRE(std::isfinite(std::real(even_pole)));
    REQUIRE(std::isfinite(std::imag(even_pole)));
    REQUIRE(std::isfinite(std::real(odd_pole)));
    REQUIRE(std::isfinite(std::imag(odd_pole)));
    REQUIRE(std::imag(even_pole) > 0.0);
    REQUIRE(std::imag(odd_pole) < 0.0);
    REQUIRE(std::abs(even_pole) > std::abs(wider_even_pole));
  }

  SECTION("invalid inputs are rejected") {
    const double infinity = std::numeric_limits<double>::infinity();
    const double quiet_nan = std::numeric_limits<double>::quiet_NaN();

    REQUIRE_THROWS_AS(gra::regge::Eta(infinity, Signature::Positive),
                      std::invalid_argument);
    REQUIRE_THROWS_AS(gra::regge::Eta(1.0, static_cast<Signature>(17)),
                      std::invalid_argument);
    REQUIRE_THROWS_AS(gra::regge::Eta(1.0, Signature::Positive, 0.0),
                      std::invalid_argument);
    REQUIRE_THROWS_AS(gra::regge::Eta(1.0, Signature::Positive, -1e-8),
                      std::invalid_argument);
    REQUIRE_THROWS_AS(gra::regge::Eta(1.0, Signature::Positive, infinity),
                      std::invalid_argument);
    REQUIRE_THROWS_AS(gra::regge::Eta(1.0, Signature::Positive, quiet_nan),
                      std::invalid_argument);
  }
}

TEST_CASE("gra::regge::EtaRaw unregulated signature factor", "[ReggeEta]") {
  using gra::regge::Rim;
  using gra::regge::Signature;
  SECTION("raw values retain the full moving signature factor") {
    struct SignaturePoint {
      double alpha;
      Signature signature;
      double real;
      double imag;
    };
    const std::array<SignaturePoint, 4> points = {
        SignaturePoint{0.5, Signature::Positive, -1.0, 1.0},
        SignaturePoint{0.5, Signature::Negative, -1.0, -1.0},
        SignaturePoint{1.5, Signature::Positive, 1.0, 1.0},
        SignaturePoint{1.5, Signature::Negative, 1.0, -1.0}};

    for (const SignaturePoint &point : points) {
      CAPTURE(point.alpha, gra::regge::Tau(point.signature));
      const std::complex<double> eta =
          gra::regge::EtaRaw(point.alpha, point.signature, Rim::Lower);
      REQUIRE(std::real(eta) == Approx(point.real).margin(2e-14));
      REQUIRE(std::imag(eta) == Approx(point.imag).margin(2e-14));
    }
  }

  SECTION("raw physical poles are unregulated") {
    const std::complex<double> even_pole =
        gra::regge::EtaRaw(0.0, Signature::Positive, Rim::Lower);
    const std::complex<double> odd_pole =
        gra::regge::EtaRaw(1.0, Signature::Negative, Rim::Lower);
    const std::complex<double> near_even_pole =
        gra::regge::EtaRaw(1.0e-12, Signature::Positive, Rim::Lower);

    REQUIRE_FALSE((std::isfinite(std::real(even_pole)) &&
                   std::isfinite(std::imag(even_pole))));
    REQUIRE_FALSE((std::isfinite(std::real(odd_pole)) &&
                   std::isfinite(std::imag(odd_pole))));
    REQUIRE(std::abs(near_even_pole) > 1.0e11);
  }
}

TEST_CASE("gra::regge::EtaPhase reduced signature factor", "[ReggeEta]") {
  using gra::regge::Signature;
  SECTION("rotating phases are pole-stripped and have unit modulus") {
    for (const double alpha : {0.0, 0.5, 1.0, 1.5}) {
      for (const Signature signature :
           {Signature::Positive, Signature::Negative}) {
        CAPTURE(alpha, gra::regge::Tau(signature));
        const std::complex<double> eta = gra::regge::EtaPhase(alpha, signature);
        const std::complex<double> shifted =
            gra::regge::EtaPhase(alpha + 2.0, signature);
        REQUIRE(std::abs(eta) == Approx(1.0).margin(1e-14));
        REQUIRE(std::real(shifted) == Approx(-std::real(eta)).margin(1e-14));
        REQUIRE(std::imag(shifted) == Approx(-std::imag(eta)).margin(1e-14));
      }
    }
  }

  SECTION("rotating phases have fixed integer anchors") {
    const std::complex<double> even_rotating =
        gra::regge::EtaPhase(0.0, Signature::Positive);
    const std::complex<double> odd_rotating =
        gra::regge::EtaPhase(0.0, Signature::Negative);

    REQUIRE(std::real(even_rotating) == Approx(-1.0).margin(1e-14));
    REQUIRE(std::imag(even_rotating) == Approx(0.0).margin(1e-14));
    REQUIRE(std::real(odd_rotating) == Approx(0.0).margin(1e-14));
    REQUIRE(std::imag(odd_rotating) == Approx(-1.0).margin(1e-14));
  }

  SECTION(
      "reduced phase inputs require a finite trajectory and exact signature") {
    const double infinity = std::numeric_limits<double>::infinity();
    REQUIRE_FALSE(std::isfinite(std::norm(gra::regge::EtaPhase(infinity, Signature::Positive))));
    REQUIRE_THROWS_AS(gra::regge::CheckEta(1.0, 0.25, static_cast<Signature>(17),
                                         gra::EtaMode::Rotating, "eta"), std::invalid_argument);
  }
}

TEST_CASE("gra::regge::EtaFactor complete prescriptions", "[ReggeEta]") {
  using gra::regge::Rim;
  using gra::regge::Signature;
  constexpr double alpha0 = 1.08;
  constexpr double alpha = 0.91;
  constexpr Signature signature = Signature::Positive;

  SECTION("fixed and moving modes dispatch exactly") {
    CHECK(gra::regge::EtaFactor(alpha, alpha0, signature, gra::EtaMode::Raw) ==
          gra::regge::EtaRaw(alpha, signature, Rim::Lower));
    CHECK(gra::regge::EtaFactor(alpha, alpha0, signature,
                                gra::EtaMode::RotatingT0) ==
          gra::regge::EtaPhase(alpha0, signature));
    CHECK(gra::regge::EtaFactor(alpha, alpha0, signature,
                                gra::EtaMode::Rotating) ==
          gra::regge::EtaPhase(alpha, signature));
  }

  SECTION("mode names and trajectory constraints are exact") {
    CHECK(gra::regge::ParseEta("raw", "eta") == gra::EtaMode::Raw);
    CHECK(gra::regge::EtaName(gra::EtaMode::Raw) == "raw");
    CHECK(gra::regge::EtaName(gra::EtaMode::RotatingT0) == "rotating_t0");
    REQUIRE_THROWS(gra::regge::ParseEta("gamma", "eta"));
    REQUIRE_THROWS(gra::regge::ParseEta("full_t0", "eta"));
    REQUIRE_THROWS(gra::regge::ParseEta("unknown", "eta"));
    CHECK(gra::regge::ParseSignature(1, "eta.tau") == Signature::Positive);
    CHECK(gra::regge::ParseSignature(-1, "eta.tau") == Signature::Negative);
    REQUIRE_THROWS(gra::regge::ParseSignature(0, "eta.tau"));
  }

  SECTION("sampling failures use amplitude bookkeeping") {
    for (const double invalid : {std::numeric_limits<double>::infinity(),
                                 std::numeric_limits<double>::quiet_NaN()}) {
      for (const auto mode : {gra::EtaMode::Raw, gra::EtaMode::RotatingT0, gra::EtaMode::Rotating}) {
        REQUIRE_THROWS_AS(gra::regge::EtaFactor(invalid, alpha0, signature, mode), gra::AmplitudeFailure);
      }
    }
    REQUIRE_FALSE(std::isfinite(std::norm(gra::regge::EtaFactor(0.0, alpha0, Signature::Positive, gra::EtaMode::Raw))));
    REQUIRE_FALSE(std::isfinite(std::norm(gra::regge::EtaFactor(1.0, alpha0, Signature::Negative, gra::EtaMode::Raw))));
    REQUIRE_NOTHROW(gra::regge::EtaFactor(1.0, alpha0, Signature::Positive, gra::EtaMode::Raw));
    REQUIRE_NOTHROW(gra::regge::EtaFactor(0.0, alpha0, Signature::Negative, gra::EtaMode::Raw));
  }

  SECTION("invalid steering remains fatal") {
    REQUIRE_THROWS_AS(gra::regge::EtaFactor(alpha, alpha0, signature, static_cast<gra::EtaMode>(17)),
                      std::invalid_argument);
    REQUIRE_THROWS_AS(gra::regge::CheckEta(alpha0, 0.25, static_cast<Signature>(17), gra::EtaMode::Raw, "eta"),
                      std::invalid_argument);
  }
}
