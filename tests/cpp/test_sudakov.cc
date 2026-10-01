// Sudakov model tests
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "support/models_test_support.hh"

#include "Graniitti/Math/MFloat.h"

#include <bit>
#include <cstdint>

TEST_CASE("MSudakovStore shares initialized tables across threads",
          "[gra::MSudakov][threading]") {
  gra::MODELPARAM = "TUNE0";

  constexpr double sqrts = 13000.0;
  const std::string pdfset = "MMHT2014lo68cl";
  const auto model_tune = gra::MModelTune::Load(
      gra::ResolveModelDataFile(gra::MODELPARAM, "GENERAL.json"));
  const auto soft_model = model_tune->Soft();
  gra::MSudakovStore store;
  const std::shared_ptr<const gra::MSudakov> first =
      store.GetSudakov(sqrts, pdfset, soft_model);
  const std::shared_ptr<const gra::MSudakov> second =
      store.GetSudakov(sqrts, pdfset, soft_model);
  REQUIRE(first == second);
  REQUIRE(first->initialized);

  constexpr std::size_t nthreads = 8;
  std::vector<std::shared_ptr<const gra::MSudakov>> handles(nthreads);
  std::vector<double> alphas(nthreads, 0.0);
  std::vector<double> fluxes(nthreads, 0.0);
  std::vector<std::thread> workers;
  workers.reserve(nthreads);

  for (std::size_t i = 0; i < nthreads; ++i) {
    workers.emplace_back(
        [i, &store, &handles, &alphas, &fluxes, &soft_model, pdfset, sqrts] {
          handles[i] = store.GetSudakov(sqrts, pdfset, soft_model);
          alphas[i] = handles[i]->AlphaS_Q2(10.0);
          fluxes[i] = handles[i]->fg_xQ2Mu(1e-3, 10.0, 10.0);
        });
  }
  for (auto &worker : workers) {
    worker.join();
  }

  for (std::size_t i = 0; i < nthreads; ++i) {
    REQUIRE(handles[i] == first);
    REQUIRE(alphas[i] == Approx(first->AlphaS_Q2(10.0)));
    REQUIRE(fluxes[i] == Approx(first->fg_xQ2Mu(1e-3, 10.0, 10.0)));
  }
}

// Require numerical stores to use a document captured with the model snapshot
TEST_CASE("MSudakovStore rejects an incomplete model snapshot",
          "[gra::MSudakov][snapshot]") {
  gra::MODELPARAM = "TUNE0";
  const std::string model_file =
      gra::ResolveModelDataFile("TUNE0", "GENERAL.json");
  const auto soft_model = gra::SoftModel::LoadFromJson(
      model_file, gra::aux::GetInputData(model_file));
  gra::MSudakovStore store;
  REQUIRE_THROWS(store.GetSudakov(13000.0, "MMHT2014lo68cl", soft_model));
}

TEST_CASE("MSudakovNumerics rejects malformed grid steering",
          "[gra::MSudakov][params]") {
  gra::MSudakovNumerics valid;
  valid.q2_MAX = 125.0;
  valid.SudakovIntegralN = 16;
  valid.ShuvaevIntegralN = 16;
  valid.SUDA_log_ON = {true, true};
  valid.SHUV_log_ON = {true, true};
  valid.SUDA_N = {8, 8};
  valid.SHUV_N = {8, 8};
  REQUIRE_NOTHROW(valid.Validate());

  SECTION("short log grid vector") {
    auto invalid = valid;
    invalid.SUDA_log_ON = {true};
    REQUIRE_THROWS_AS(invalid.Validate(), std::invalid_argument);
  }
  SECTION("extra grid-size entry") {
    auto invalid = valid;
    invalid.SHUV_N = {8, 8, 8};
    REQUIRE_THROWS_AS(invalid.Validate(), std::invalid_argument);
  }
  SECTION("degenerate grid size") {
    auto invalid = valid;
    invalid.SUDA_N = {1, 8};
    REQUIRE_THROWS_AS(invalid.Validate(), std::invalid_argument);
  }
  SECTION("grid size outside signed interpolation indices") {
    auto invalid = valid;
    invalid.SUDA_N[0] = std::numeric_limits<unsigned int>::max();
    REQUIRE_THROWS_AS(invalid.Validate(), std::invalid_argument);
    gra::IArray2D grid;
    REQUIRE_THROWS_AS(grid.Set(0, "q2", 1.0, 10.0, invalid.SUDA_N[0], true),
                      std::invalid_argument);
  }
  SECTION("grid counts are validated before JSON conversion") {
    const std::string file = gra::ResolveModelDataFile("TUNE0", "NUMERICS.json");
    const auto original = nlohmann::json::parse(gra::aux::GetInputData(file));
    const auto invalid_counts = nlohmann::json::parse(R"([-1,2.5,4294967304,2147483647,true,"8"])");
    for (const std::string axis : {"SUDA", "SHUV"}) {
      for (const auto &count : invalid_counts) {
        auto card = original;
        card["NUMERICS_SUDAKOV"][axis]["N"][0] = count;
        gra::MSudakovNumerics numerics;
        REQUIRE_THROWS_AS(numerics.ConfigureFromJson(file, card.dump()), std::invalid_argument);
      }
    }
  }
  SECTION("integral intervals must be even integers") {
    const std::string file     = gra::ResolveModelDataFile("TUNE0", "NUMERICS.json");
    const auto        original = nlohmann::json::parse(gra::aux::GetInputData(file));
    for (const std::string axis : {"SUDA", "SHUV"}) {
      for (const auto& count : nlohmann::json::parse(R"([0,1,3,-2,2.5,4294967304,true,"8"])")) {
        auto card                               = original;
        card["NUMERICS_SUDAKOV"][axis]["N_int"] = count;
        gra::MSudakovNumerics numerics;
        REQUIRE_THROWS_AS(numerics.ConfigureFromJson(file, card.dump()), std::invalid_argument);
      }
    }
  }
  SECTION("non-positive upper q2 boundary") {
    auto invalid = valid;
    invalid.q2_MAX = 0.0;
    REQUIRE_THROWS_AS(invalid.Validate(), std::invalid_argument);
  }
}

TEST_CASE("MSudakovModel rejects inadmissible infrared models",
          "[gra::MSudakov][params]") {
  gra::MSudakovModel valid;
  valid.mode_name = "IR_RKHS";
  valid.Q0 = 1.0;
  valid.ir_logq2_length = 1.0;
  valid.rkhs_logq2_length = 1.0;
  valid.rkhs_logx_length = 2.0;
  valid.rkhs_weight_exponent = 0.5;
  valid.rkhs_radius = 1.0;
  valid.member_s_anchors = {1.0};
  valid.member_x_anchors = {1e-3};
  valid.member_coefficients = {0.1};
  REQUIRE_NOTHROW(valid.Validate());

  SECTION("unknown mode") {
    auto invalid = valid;
    invalid.mode_name = "AD_HOC";
    REQUIRE_THROWS(invalid.Validate());
  }
  SECTION("non-positive matching scale") {
    auto invalid = valid;
    invalid.Q0 = 0.0;
    REQUIRE_THROWS(invalid.Validate());
  }
  SECTION("non-positive baseline correlation length") {
    auto invalid = valid;
    invalid.ir_logq2_length = 0.0;
    REQUIRE_THROWS(invalid.Validate());
  }
  SECTION("non-localized functional space") {
    auto invalid = valid;
    invalid.rkhs_weight_exponent = 0.0;
    REQUIRE_THROWS(invalid.Validate());
  }
  SECTION("non-positive logarithmic x correlation") {
    auto invalid = valid;
    invalid.rkhs_logx_length = 0.0;
    REQUIRE_THROWS(invalid.Validate());
  }
  SECTION("ill-conditioned functional anchors") {
    auto invalid = valid;
    invalid.member_s_anchors = {1.0, std::nextafter(1.0, 2.0)};
    invalid.member_x_anchors = {1e-3, std::nextafter(1e-3, 2e-3)};
    invalid.member_coefficients = {0.1, 0.1};
    REQUIRE_THROWS(invalid.Validate());
  }
  SECTION("mismatched functional vectors") {
    auto invalid = valid;
    invalid.member_coefficients.clear();
    REQUIRE_THROWS(invalid.Validate());
    REQUIRE_THROWS_AS(invalid.MemberNorm2(), std::invalid_argument);
  }
  SECTION("direct member norm rejects missing anchors") {
    auto invalid = valid;
    invalid.member_s_anchors.clear();
    REQUIRE_THROWS_AS(invalid.MemberNorm2(), std::invalid_argument);
    invalid = valid;
    invalid.member_x_anchors.clear();
    REQUIRE_THROWS_AS(invalid.MemberNorm2(), std::invalid_argument);
  }
  SECTION("member outside weighted-RKHS ball") {
    auto invalid = valid;
    invalid.rkhs_radius = 1e-3;
    REQUIRE_THROWS(invalid.Validate());
  }
  SECTION("functional member in central-only mode") {
    auto invalid = valid;
    invalid.mode_name = "IR_BASELINE";
    REQUIRE_THROWS(invalid.Validate());
  }
  SECTION("functional radius in central-only mode") {
    auto invalid = valid;
    invalid.mode_name = "IR_BASELINE";
    invalid.member_coefficients = {0.0};
    REQUIRE_THROWS(invalid.Validate());
  }
}


TEST_CASE("MSudakov derives matching and flavour thresholds from LHAPDF",
          "[gra::MSudakov][params]") {
  ModelParamRestoreGuard restore;
  gra::MODELPARAM = "TUNE0";
  gra::MLHAPDFStore pdf_store;
  const auto pdf = pdf_store.GetPDF("MMHT2014lo68cl", 0);
  const double qmin = pdf->info().get_entry_as<double>("QMin");
  const double charm = pdf->info().get_entry_as<double>("MCharm");
  const double bottom = pdf->info().get_entry_as<double>("MBottom");

  const std::string general_file =
      gra::ResolveModelDataFile(gra::MODELPARAM, "GENERAL.json");
  const std::string numerics_file =
      gra::ResolveModelDataFile(gra::MODELPARAM, "NUMERICS.json");
  const auto soft_model = gra::MModelTune::Load(general_file)->Soft();
  gra::MSudakov sudakov;
  REQUIRE_NOTHROW(sudakov.Init(20.0, "MMHT2014lo68cl", soft_model,
                               numerics_file,
                               gra::aux::GetInputData(numerics_file), false));
  gra::MSudakovModel configured_model;
  configured_model.ConfigureFromJson(general_file, soft_model->SourceJson());
  REQUIRE(sudakov.GetQMin() == Approx(configured_model.Q0));
  REQUIRE(sudakov.GetQ2Min() == Approx(gra::math::pow2(configured_model.Q0)));
  REQUIRE(sudakov.NumFlavor(std::nextafter(gra::math::pow2(charm), 0.0)) ==
          Approx(3.0));
  REQUIRE(sudakov.NumFlavor(gra::math::pow2(charm)) == Approx(4.0));
  REQUIRE(sudakov.NumFlavor(std::nextafter(gra::math::pow2(bottom), 0.0)) ==
          Approx(4.0));
  REQUIRE(sudakov.NumFlavor(gra::math::pow2(bottom)) == Approx(5.0));

  const double physical_q0 = qmin + 0.5;
  const std::string q0_tune = WriteModifiedSudakovModelTune(
      "physical_q0",
      [&](auto &j) { j["PARAM_SKEWED_UGD"]["Q0"] = physical_q0; },
      [](auto &) {});
  const std::string q0_general_file =
      gra::ResolveModelDataFile(q0_tune, "GENERAL.json");
  const std::string q0_numerics_file =
      gra::ResolveModelDataFile(q0_tune, "NUMERICS.json");
  gra::MSudakov q0_sudakov;
  REQUIRE_NOTHROW(q0_sudakov.Init(
      20.0, "MMHT2014lo68cl",
      gra::MModelTune::Load(q0_general_file)->Soft(),
      q0_numerics_file, gra::aux::GetInputData(q0_numerics_file), false));
  REQUIRE(q0_sudakov.GetQMin() == Approx(physical_q0));

  const std::string invalid_q0_tune = WriteModifiedSudakovModelTune(
      "invalid_q0", [&](auto &j) { j["PARAM_SKEWED_UGD"]["Q0"] = 0.9 * qmin; },
      [](auto &) {});
  const std::string invalid_general_file =
      gra::ResolveModelDataFile(invalid_q0_tune, "GENERAL.json");
  const std::string invalid_numerics_file =
      gra::ResolveModelDataFile(invalid_q0_tune, "NUMERICS.json");
  gra::MSudakov invalid_q0;
  REQUIRE_THROWS(
      invalid_q0.Init(20.0, "MMHT2014lo68cl",
                      gra::MModelTune::Load(invalid_general_file)->Soft(),
                      invalid_numerics_file,
                      gra::aux::GetInputData(invalid_numerics_file), false));
}

TEST_CASE("IArray2D throws outside its interpolation domain",
          "[gra::MSudakov][IArray2D]") {
  gra::IArray2D grid;
  grid.Set(0, "q2", 1.0, 3.0, 2, false);
  grid.Set(1, "mu", 2.0, 4.0, 2, false);
  grid.InitArray();
  for (std::size_t i = 0; i < grid.F.size(); ++i) {
    for (std::size_t j = 0; j < grid.F[i].size(); ++j) {
      const double q2 = grid.MIN[0] + i * grid.STEP[0];
      const double mu = grid.MIN[1] + j * grid.STEP[1];
      grid.F[i][j] = {q2, mu, q2 + 2.0 * mu, 3.0 * q2 - mu};
    }
  }

  const auto value = grid.Interpolate2D(2.0, 3.0);
  REQUIRE(value.first == Approx(8.0));
  REQUIRE(value.second == Approx(3.0));
  REQUIRE_NOTHROW(grid.Interpolate2D(std::nextafter(1.0, 0.0), 2.0));
  REQUIRE_THROWS_AS(grid.Interpolate2D(0.999, 3.0), std::out_of_range);
  REQUIRE_THROWS_AS(grid.Interpolate2D(3.001, 3.0), std::out_of_range);
  REQUIRE_THROWS_AS(grid.Interpolate2D(2.0, 1.999), std::out_of_range);
  REQUIRE_THROWS_AS(grid.Interpolate2D(2.0, 4.001), std::out_of_range);
}

// Check interpolation below the former small-x cutoff and reject invalid bounds
TEST_CASE("IArray2D supports positive small-x logarithmic grids",
          "[gra::MSudakov][IArray2D]") {
  gra::IArray2D grid;
  const double xmin = gra::math::pow2(1.7 / 100000.0);
  grid.Set(0, "q2", 1.0, 4.0, 2, true);
  REQUIRE_NOTHROW(grid.Set(1, "x", xmin, 0.1, 4, true));
  grid.InitArray();
  for (const auto &i : indices(grid.F)) {
    for (const auto &j : indices(grid.F[i])) {
      const double u = grid.MIN[0] + i * grid.STEP[0];
      const double v = grid.MIN[1] + j * grid.STEP[1];
      grid.F[i][j] = {u, v, u + 2.0 * v, std::exp(-u)};
    }
  }
  for (const double x : {xmin, 2.0 * xmin, 1e-9}) {
    const auto value = grid.Interpolate2D(2.0, x);
    REQUIRE(value.first == Approx(std::log(2.0) + 2.0 * std::log(x)).epsilon(1e-12));
    REQUIRE(value.second == Approx(0.5).epsilon(1e-12));
  }
  for (const double bound : {0.0, -1.0, std::numeric_limits<double>::quiet_NaN(),
                             std::numeric_limits<double>::infinity()}) {
    REQUIRE_THROWS_AS(grid.Set(1, "x", bound, 1.0, 4, true), std::invalid_argument);
  }
  REQUIRE_THROWS_AS(grid.Set(1, "x", xmin, std::numeric_limits<double>::infinity(), 4, true),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(grid.Set(1, "x", 1e-100, std::nextafter(1e-100, 1.0), 4, true),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(grid.Set(0, "q2", 1.0, std::nextafter(1.0, 2.0), 4, false),
                    std::invalid_argument);
}

TEST_CASE("IArray2D returns the derivative of its Hermite interpolant",
          "[gra::MSudakov][IArray2D]") {
  gra::IArray2D grid;
  grid.Set(0, "q2", 1.0, 9.0, 2, true);
  grid.Set(1, "x", 0.1, 0.3, 2, false);
  grid.InitArray();
  for (std::size_t i = 0; i < grid.F.size(); ++i) {
    for (std::size_t j = 0; j < grid.F[i].size(); ++j) {
      const double log_q2 = grid.MIN[0] + i * grid.STEP[0];
      const double x = grid.MIN[1] + j * grid.STEP[1];
      const double q2 = std::exp(log_q2);
      grid.F[i][j] = {log_q2, x, q2 * q2 + x, 2.0 * q2};
    }
  }

  const double q2 = 2.3;
  const double x = 0.17;
  const double step = 1e-5 * q2;
  double       second    = 0.0;
  const auto   value     = grid.Interpolate2D(q2, x, &second);
  const double numerical = (grid.Interpolate2D(q2 + step, x).first -
                            grid.Interpolate2D(q2 - step, x).first) /
                           (2.0 * step);
  REQUIRE(value.second == Approx(numerical).epsilon(1e-9));
  const double curvature =
      (grid.Interpolate2D(q2 + step, x).second - grid.Interpolate2D(q2 - step, x).second) / (2.0 * step);
  REQUIRE(second == Approx(curvature).epsilon(1e-9));
}

TEST_CASE("IArray2D atomically publishes exact cache shapes",
          "[gra::MSudakov][IArray2D][threading]") {
  auto make_grid = [](double offset) {
    gra::IArray2D grid;
    grid.Set(0, "q2", 1.0, 3.0, 2, true);
    grid.Set(1, "x", 0.1, 0.3, 2, false);
    grid.InitArray();
    for (std::size_t i = 0; i < grid.F.size(); ++i) {
      for (std::size_t j = 0; j < grid.F[i].size(); ++j) {
        grid.F[i][j] = {grid.MIN[0] + i * grid.STEP[0],
                        grid.MIN[1] + j * grid.STEP[1], offset + 10.0 * i + j,
                        offset - 10.0 * i - j};
      }
    }
    return grid;
  };

  const auto first = make_grid(0.0);
  const auto second = make_grid(1000.0);
  const std::string filename = "tmp/sudakov_nested/cache/atomic.json";
  std::atomic<unsigned int> failures{0};
  std::thread writer_a([&] {
    try {
      first.WriteArray(filename, true);
    } catch (...) {
      ++failures;
    }
  });
  std::thread writer_b([&] {
    try {
      second.WriteArray(filename, true);
    } catch (...) {
      ++failures;
    }
  });
  writer_a.join();
  writer_b.join();
  REQUIRE(failures.load() == 0);

  auto loaded = make_grid(-1.0);
  REQUIRE(loaded.ReadArray(filename));
  const double offset = loaded.F[0][0][2];
  REQUIRE((offset == Approx(0.0) || offset == Approx(1000.0)));
  for (std::size_t i = 0; i < loaded.F.size(); ++i) {
    for (std::size_t j = 0; j < loaded.F[i].size(); ++j) {
      REQUIRE(std::bit_cast<std::uint64_t>(loaded.F[i][j][0]) ==
              std::bit_cast<std::uint64_t>(first.F[i][j][0]));
      REQUIRE(loaded.F[i][j][2] == Approx(offset + 10.0 * i + j));
      REQUIRE(loaded.F[i][j][3] == Approx(offset - 10.0 * i - j));
    }
  }

  auto mismatched = loaded;
  mismatched.islog[0] = false;
  REQUIRE_FALSE(mismatched.ReadArray(filename));
  REQUIRE(loaded.ReadArray(filename));

  std::ofstream corrupted(filename, std::ios::out | std::ios::trunc);
  corrupted << "0,0,0,0\n";
  corrupted.close();
  REQUIRE_FALSE(loaded.ReadArray(filename));
}

TEST_CASE("MSudakov infrared modes obey smooth matching and coherent-member "
          "constraints",
          "[gra::MSudakov][infrared]") {
  ModelParamRestoreGuard restore;
  auto make_tune = [](const std::string &suffix, const std::string &mode,
                      double radius, double coefficient, double weight = 0.5) {
    return WriteModifiedSudakovModelTune(
        suffix,
        [&](auto &j) {
          auto &model = j["PARAM_SKEWED_UGD"];
          model["mode"] = mode;
          model["Q0"] = 1.0;
          model["ir_logq2_length"] = 1.0;
          model["rkhs_logq2_length"] = 1.0;
          model["rkhs_logx_length"] = 2.0;
          model["rkhs_weight_exponent"] = weight;
          model["rkhs_radius"] = radius;
          model["member_s_anchors"] = {1.0};
          model["member_x_anchors"] = {1e-2};
          model["member_coefficients"] = {coefficient};
        },
        [](auto &j) {
          auto &s = j["NUMERICS_SUDAKOV"];
          s["q2_MAX"] = 5.0;
          s["SUDA"]["N"] = {8, 6};
          s["SHUV"]["N"] = {8, 6};
        });
  };
  auto initialize = [](gra::MSudakov &sudakov, const std::string &tune) {
    const std::string general_file =
        gra::ResolveModelDataFile(tune, "GENERAL.json");
    const std::string numerics_file =
        gra::ResolveModelDataFile(tune, "NUMERICS.json");
    sudakov.Init(20.0, "MMHT2014lo68cl",
                 gra::MModelTune::Load(general_file)->Soft(), numerics_file,
                 gra::aux::GetInputData(numerics_file), true);
  };

  const std::string perturbative_tune =
      make_tune("perturbative", "PERTURBATIVE_ONLY", 0.0, 0.0);
  const std::string minimum_tune =
      make_tune("minimum", "IR_BASELINE", 0.0, 0.0);
  const std::string functional_tune =
      make_tune("functional", "IR_RKHS", 0.5, 0.1);
  const std::string tail_tune =
      make_tune("functional_tail", "IR_RKHS", 0.5, 1e5, 12.0);
  gra::MSudakov perturbative;
  gra::MSudakov minimum;
  gra::MSudakov functional;
  gra::MSudakov tail;
  REQUIRE_NOTHROW(initialize(perturbative, perturbative_tune));
  REQUIRE_NOTHROW(initialize(minimum, minimum_tune));
  REQUIRE_NOTHROW(initialize(functional, functional_tune));
  REQUIRE_NOTHROW(initialize(tail, tail_tune));

  const double q02 = minimum.GetQ2Min();
  const double x = 0.01;
  const double mu = 3.0;
  const double q2 = 0.4 * q02;

  // Compare the event-local cached evaluator with the full matched flux path
  auto check_prepared_flux = [&](gra::MSudakov &sudakov) {
    auto prepared = sudakov.PrepareFlux(x, mu);
    REQUIRE(gra::math::IsZero(
        prepared.OverQ2(std::numeric_limits<double>::quiet_NaN())));
    for (const double invalid : {std::numeric_limits<double>::quiet_NaN(),
                                 std::numeric_limits<double>::infinity(), -1.0}) {
      REQUIRE(gra::math::IsZero(sudakov.fg_xQ2MuOverQ2(x, invalid, mu)));
    }
    const std::array<double, 7> points = {
        0.0,       q02 * 1e-10, q2, q02, std::nextafter(q02, 2.0 * q02),
        2.0 * q02, 4.5};
    for (const double point : points) {
      CAPTURE(point);
      REQUIRE(prepared.OverQ2(point) ==
              Approx(sudakov.fg_xQ2MuOverQ2(x, point, mu))
                  .epsilon(2e-12)
                  .margin(1e-13));
    }
    // Check the finite quotient without first forming the underflowing gluon flux
    const double limit = prepared.OverQ2(0.0);
    for (const double point : {std::numeric_limits<double>::denorm_min(),
                               std::numeric_limits<double>::min(), 1e-250}) {
      REQUIRE(prepared.OverQ2(point) == Approx(limit).epsilon(1e-10));
      REQUIRE(sudakov.fg_xQ2MuOverQ2(x, point, mu) == Approx(limit).epsilon(1e-10));
    }
    // Compare with the separately evaluated flux where division remains well resolved
    REQUIRE(prepared.OverQ2(q2) ==
            Approx(sudakov.fg_xQ2Mu(x, q2, mu) / q2).epsilon(2e-12));
  };
  check_prepared_flux(perturbative);
  check_prepared_flux(minimum);
  check_prepared_flux(functional);
  check_prepared_flux(tail);

  REQUIRE(q02 == Approx(gra::math::pow2(minimum.GetQMin())));
  REQUIRE(perturbative.fg_xQ2Mu(x, q2, mu) == Approx(0.0));
  REQUIRE_THROWS_AS(perturbative.IntegratedFlux_xQ2Mu(x, q2, mu),
                    std::domain_error);
  REQUIRE(minimum.fg_xQ2Mu(x, q2, mu) > 0.0);
  REQUIRE(minimum.IntegratedFlux_xQ2Mu(x, q2, mu) > 0.0);

  const double epsilon = 1e-5;
  const double left =
      minimum.IntegratedFlux_xQ2Mu(x, q02 * std::exp(-epsilon), mu);
  const double right =
      minimum.IntegratedFlux_xQ2Mu(x, q02 * std::exp(epsilon), mu);
  const double boundary_flux = minimum.fg_xQ2Mu(x, q02, mu);
  const double numerical_flux = (right - left) / (2.0 * epsilon);
  REQUIRE(numerical_flux == Approx(boundary_flux).epsilon(5e-3));
  const double flux_left = minimum.fg_xQ2Mu(x, q02 * std::exp(-epsilon), mu);
  const double flux_center = minimum.fg_xQ2Mu(x, q02, mu);
  const double flux_right = minimum.fg_xQ2Mu(x, q02 * std::exp(epsilon), mu);
  const double slope_left = (flux_center - flux_left) / epsilon;
  const double slope_right = (flux_right - flux_center) / epsilon;
  REQUIRE(slope_left == Approx(slope_right).epsilon(2e-2).margin(1e-8));

  const unsigned int intervals = 1000;
  const double q2_low = q02 * 1e-4;
  const double log_step = std::log(q02 / q2_low) / intervals;
  double flux_integral = 0.0;
  for (unsigned int i = 0; i <= intervals; ++i) {
    const double point = q2_low * std::exp(i * log_step);
    const double weight = (i == 0 || i == intervals) ? 0.5 : 1.0;
    flux_integral += weight * minimum.fg_xQ2Mu(x, point, mu);
  }
  flux_integral *= log_step;
  const double integrated_difference =
      minimum.IntegratedFlux_xQ2Mu(x, q02, mu) -
      minimum.IntegratedFlux_xQ2Mu(x, q2_low, mu);
  REQUIRE(flux_integral == Approx(integrated_difference).epsilon(2e-4));

  const double origin_limit = minimum.fg_xQ2MuOverQ2(x, 0.0, mu);
  const double near_origin = minimum.fg_xQ2MuOverQ2(x, q02 * 1e-10, mu);
  REQUIRE(origin_limit > 0.0);
  REQUIRE(near_origin == Approx(origin_limit).epsilon(1e-6));
  REQUIRE_THROWS_AS(minimum.AlphaS_Q2(0.9 * q02), std::out_of_range);
  REQUIRE_THROWS_AS(minimum.PerturbativeShuvaevGluon_xQ2(x, 0.9 * q02),
                    std::out_of_range);

  REQUIRE(functional.GetIRMode() == gra::MSudakovIRMode::IR_RKHS);
  REQUIRE(functional.GetMemberNorm2() <= 0.25 + 1e-12);
  const double kernel_base = std::exp(-1.0);
  const double kernel_c0 = std::exp(-1.0);
  const double kernel_c1 = 0.5 * kernel_c0;
  const double kernel_projection =
      kernel_c0 * (1.25 * kernel_c0 + 0.5 * kernel_c1) +
      kernel_c1 * (0.5 * kernel_c0 + kernel_c1);
  REQUIRE(functional.GetMemberNorm2() ==
          Approx(0.01 * (kernel_base - kernel_projection)).epsilon(1e-13));
  const double selected = functional.fg_xQ2Mu(x, q2, mu);
  REQUIRE(functional.AlphaSFlux_xQ2Mu(x, q2, mu, q2) ==
          Approx(functional.AlphaS_Q2(q02) * selected).epsilon(1e-13));
  REQUIRE(selected != Approx(minimum.fg_xQ2Mu(x, q2, mu)));
  REQUIRE(selected >= 0.0);

  const auto opposite_member = functional.WithMemberCoefficients({-0.1});
  REQUIRE(opposite_member->GetIRMode() == gra::MSudakovIRMode::IR_RKHS);
  REQUIRE(opposite_member->fg_xQ2Mu(x, q2, mu) != Approx(selected));
  REQUIRE(opposite_member->fg_xQ2Mu(x, q2, mu) >= 0.0);

  // Require nonzero members to preserve f_g = dG/dln(Q2) inside the infrared domain
  for (const gra::MSudakov *member :
       std::array<const gra::MSudakov *, 3>{&functional, opposite_member.get(), &tail}) {
    const double step = 1e-6;
    for (const double fraction : {0.4, 0.8, 0.99}) {
      const double point = fraction * q02;
      const double numerical =
          (member->IntegratedFlux_xQ2Mu(x, point * std::exp(step), mu) -
           member->IntegratedFlux_xQ2Mu(x, point * std::exp(-step), mu)) /
          (2.0 * step);
      REQUIRE(numerical == Approx(member->fg_xQ2Mu(x, point, mu)).epsilon(2e-6));
    }
    const double matched_derivative =
        (member->IntegratedFlux_xQ2Mu(x, q02 * std::exp(step), mu) -
         member->IntegratedFlux_xQ2Mu(x, q02 * std::exp(-step), mu)) /
        (2.0 * step);
    REQUIRE(matched_derivative == Approx(member->fg_xQ2Mu(x, q02, mu)).epsilon(1e-4));
  }
  auto positive_prepared = functional.PrepareFlux(x, mu);
  auto negative_prepared = opposite_member->PrepareFlux(x, mu);
  for (unsigned int i = 0; i <= 40; ++i) {
    const double point = q02 * std::exp(-0.5 * i);
    CAPTURE(point);
    REQUIRE(positive_prepared.OverQ2(point) >= 0.0);
    REQUIRE(negative_prepared.OverQ2(point) >= 0.0);
  }
  REQUIRE_THROWS_AS(functional.WithMemberCoefficients({10.0}),
                    std::invalid_argument);

  const std::string nonpositive_ball_tune =
      make_tune("nonpositive_ball", "IR_RKHS", 100.0, 0.0);
  gra::MSudakov nonpositive_ball;
  REQUIRE_NOTHROW(initialize(nonpositive_ball, nonpositive_ball_tune));
  auto rejected_ir = nonpositive_ball.PrepareFlux(x, mu);
  REQUIRE(std::isfinite(rejected_ir.OverQ2(2.0 * q02)));
  REQUIRE_THROWS(rejected_ir.OverQ2(q2));
}

// Reject shell syntax and path components before LHAPDF lookup or download
TEST_CASE("LHAPDF access rejects invalid set identifiers",
          "[gra::MLHAPDF][params]") {
  gra::MLHAPDFStore store;
  for (const std::string name : {"", "null", "../CT10nlo", "/CT10nlo", "-CT10nlo",
                                 "bad;:", "bad$(true)", "bad`true`", "bad'name",
                                 "bad name", "bad\nname"}) {
    CAPTURE(name);
    REQUIRE_THROWS_AS(store.GetPDF(name), std::invalid_argument);
    REQUIRE_THROWS_AS(gra::aux::AutoDownloadLHAPDF(name), std::invalid_argument);
  }
  for (const std::string name : {"CT10nlo", "NNPDF31_lo_as_0118", "PDF4LHC21_40_pdfas",
                                 "NNPDF23_nlo_as_0118_qed", "set-1.0"}) {
    REQUIRE_NOTHROW(gra::aux::ValidateLHAPDFName(name));
  }
}

// Check that invalid member requests leave an installed set usable
TEST_CASE("LHAPDF rejects members outside the installed set", "[gra::MLHAPDF][params]") {
  gra::MLHAPDFStore store;
  const auto central = store.GetPDF("CT10nlo", 0);
  const LHAPDF::PDFSet set("CT10nlo");
  REQUIRE(set.size() < static_cast<std::size_t>(std::numeric_limits<int>::max()));
  REQUIRE_THROWS_AS(store.GetPDF("CT10nlo", static_cast<int>(set.size())), std::invalid_argument);
  REQUIRE_THROWS_AS(store.GetPDF("CT10nlo", std::numeric_limits<int>::max()), std::invalid_argument);
  REQUIRE_THROWS_AS(store.GetPDF("CT10nlo", -1), std::invalid_argument);
  REQUIRE(store.GetPDF("CT10nlo", 0) == central);
}

// Check shared PDF evaluation from several threads
TEST_CASE("MLHAPDFStore shares PDF objects across threads",
          "[gra::MLHAPDF][threading]") {
  gra::MLHAPDFStore store;
  const auto first = store.GetPDF("CT10nlo", 0);
  const auto second = store.GetPDF("CT10nlo", 0);
  REQUIRE(first == second);

  constexpr std::size_t nthreads = 8;
  std::vector<std::shared_ptr<const LHAPDF::PDF>> handles(nthreads);
  std::vector<std::thread> workers;
  workers.reserve(nthreads);

  for (std::size_t i = 0; i < nthreads; ++i) {
    workers.emplace_back([i, &handles, &store] {
      handles[i] = store.GetPDF("CT10nlo", 0);
    });
  }
  for (auto &worker : workers) {
    worker.join();
  }

  for (const auto &handle : handles) {
    REQUIRE(handle == first);
    REQUIRE(handle->xfxQ2(21, 1e-3, 10.0) ==
            Approx(first->xfxQ2(21, 1e-3, 10.0)));
  }
}

TEST_CASE("MSudakov preserves the signed Durham flux derivative",
          "[gra::MSudakov][flux]") {
  gra::MSudakov sudakov;

  const double central_mass = 10.0;
  const double mu_sud       = central_mass;
  for (const double kt2 : std::array<double, 5>{0.25, 1.0, 4.0, 9.0, 16.0}) {
    REQUIRE(sudakov.SudakovDelta(kt2, mu_sud) ==
            Approx(std::sqrt(kt2) / central_mass).epsilon(1e-14));
  }
  REQUIRE_THROWS_AS(sudakov.SudakovDelta(-1.0, mu_sud), std::invalid_argument);
  REQUIRE_THROWS_AS(sudakov.SudakovDelta(1.0, 0.0), std::invalid_argument);
  const double positive = sudakov.FluxDerivative(4.0, 2.0, 0.30, 0.25, 0.02);
  const double negative = sudakov.FluxDerivative(4.0, 2.0, -0.30, 0.25, 0.02);
  REQUIRE(positive == Approx(0.76).epsilon(1e-14));
  REQUIRE(negative == Approx(-0.44).epsilon(1e-14));
  REQUIRE(sudakov.FluxDerivative(4.0, 2.0, -0.30, 0.0, 0.02) == Approx(0.0));
}

// Check the empty radiation interval through the direct and prepared gluon APIs
TEST_CASE("Sudakov veto is unity above the hard scale between grid nodes", "[gra::MSudakov][flux][veto]") {
  for (const bool aligned : {false, true}) {
    const std::string tune = WriteModifiedSudakovModelTune(
        aligned ? "veto_aligned" : "veto_unaligned", [](auto &) {}, [aligned](auto &j) {
          auto &grid = j["NUMERICS_SUDAKOV"];
          grid["SUDA"]["N"] = {12, 10};
          grid["SHUV"]["N"] = aligned ? nlohmann::json{12, 10} : nlohmann::json{15, 9};
        });
    const auto model = gra::MModelTune::Load(gra::ResolveModelDataFile(tune, "GENERAL.json"));
    gra::MSudakov sudakov;
    sudakov.Init(20.0, "MMHT2014lo68cl", model->Soft(), model->NumericsFile(),
                 gra::aux::GetInputData(model->NumericsFile()));
    const double x = 0.03;
    for (const double scale : {0.5, 1.0, 1.000001, 1.001, 1.37, 2.13, 3.27}) {
      const double mu = scale * sudakov.GetMuMin();
      auto prepared = sudakov.PrepareFlux(x, mu);
      for (const double ratio : {1.0, 1.05, 1.5}) {
        const double q2 = std::max(sudakov.GetQ2Min(), ratio * mu * mu);
        CAPTURE(aligned, mu, q2);
        const double gluon = sudakov.PerturbativeShuvaevGluon_xQ2(x, q2);
        const double derivative = sudakov.PerturbativeShuvaevGluonDerivative_xQ2(x, q2);
        CHECK(sudakov.PerturbativeSqrtSudakov_Q2Mu(q2, mu) == Approx(1.0).epsilon(1.0e-13));
        CHECK(sudakov.IntegratedFlux_xQ2Mu(x, q2, mu) == Approx(gluon).epsilon(1.0e-12));
        CHECK(sudakov.fg_xQ2Mu(x, q2, mu) == Approx(q2 * derivative).epsilon(1.0e-12));
        CHECK(prepared.OverQ2(q2) == Approx(derivative).epsilon(1.0e-12));
      }
      if (0.8 * mu * mu > sudakov.GetQ2Min()) {
        const double q2 = 0.8 * mu * mu;
        const double h = 1.0e-5 * q2;
        const double derivative = (sudakov.IntegratedFlux_xQ2Mu(x, q2 + h, mu) -
                                   sudakov.IntegratedFlux_xQ2Mu(x, q2 - h, mu)) / (2.0 * h);
        const double sqrt_t = sudakov.PerturbativeSqrtSudakov_Q2Mu(q2, mu);
        CHECK(sqrt_t > 0.0);
        CHECK(sqrt_t < 1.0);
        CHECK(sudakov.fg_xQ2Mu(x, q2, mu) == Approx(q2 * derivative).epsilon(1.0e-7));
        CHECK(prepared.OverQ2(q2) == Approx(derivative).epsilon(1.0e-7));
      }
    }
  }
}

// Cover the x and hard-scale range of the Durham multijet validation
// [REFERENCE: Harland-Lang et al., arXiv:1508.02718, Figure 3]
TEST_CASE(
    "MSudakov infrared baseline remains non-negative for Durham multijets",
    "[gra::MSudakov][infrared][regression]") {
  ModelParamRestoreGuard restore;
  gra::MODELPARAM = "TUNE0";
  const std::string general_file =
      gra::ResolveModelDataFile(gra::MODELPARAM, "GENERAL.json");
  const std::string numerics_file =
      gra::ResolveModelDataFile(gra::MODELPARAM, "NUMERICS.json");
  gra::MSudakov sudakov;
  REQUIRE_NOTHROW(sudakov.Init(
      13000.0, "MMHT2014lo68cl",
      gra::MModelTune::Load(general_file)->Soft(),
      numerics_file, gra::aux::GetInputData(numerics_file), true));

  const std::array<double, 10> x_values = {2e-4, 5e-4, 1e-3, 2e-3, 5e-3,
                                           1e-2, 2e-2, 5e-2, 1e-1, 2e-1};
  const std::array<double, 6> mu_values = {20.0, 30.0, 40.0, 60.0, 80.0, 100.0};
  const std::array<double, 5> q2_fractions = {0.0, 1e-6, 0.1, 0.4, 0.9};
  for (const double x : x_values) {
    for (const double mu : mu_values) {
      CAPTURE(x, mu);
      auto flux = sudakov.PrepareFlux(x, mu);
      for (const double fraction : q2_fractions) {
        const double value = flux.OverQ2(fraction * sudakov.GetQ2Min());
        REQUIRE(std::isfinite(value));
        REQUIRE(value >= 0.0);
      }
    }
  }
}

// Check value and derivative continuity at the moving radiation boundary
TEST_CASE("Sudakov matching remains finite across the hard scale", "[gra::MSudakov][veto][regression]") {
  const auto         model = gra::MModelTune::Load(gra::ResolveModelDataFile("TUNE0", "GENERAL.json"));
  gra::MSudakovStore store;
  const auto         sud = store.GetSudakov(13000.0, "MMHT2014lo68cl", model->Soft());
  const double       q0  = sud->GetQMin();
  for (const double ratio : {0.9, 1.0, 1.000001, 1.0001, 1.001, 1.01, 1.37, 2.81}) {
    const double mu       = ratio * q0;
    auto         prepared = sud->PrepareFlux(0.001, mu);
    REQUIRE(std::isfinite(prepared.OverQ2(0.4 * q0 * q0)));
    if (ratio <= 1.0) { continue; }
    const double q2    = mu * mu;
    const double below = q2 * (1.0 - 1e-8);
    REQUIRE(sud->PerturbativeSqrtSudakov_Q2Mu(below, mu) == Approx(1.0).epsilon(1e-12));
    REQUIRE(prepared.OverQ2(below) == Approx(prepared.OverQ2(q2)).epsilon(1e-6));
  }
  for (const double q2 : {0.4, 1.0}) {
    const double lower = sud->PrepareFlux(0.001, (1.0 - 1e-7) * q0).OverQ2(q2);
    const double upper = sud->PrepareFlux(0.001, (1.0 + 1e-7) * q0).OverQ2(q2);
    REQUIRE(lower == Approx(upper).epsilon(1e-5));
  }
  for (const double mu : {2.37, 5.2, 30.0, 100.0}) {
    for (const double q2 : {q0 * q0, 4.3, 10.7, 81.0}) {
      CAPTURE(mu, q2);
      const double table = std::pow(sud->PerturbativeSqrtSudakov_Q2Mu(q2, mu), 2);
      REQUIRE(table == Approx(sud->Sudakov_T(q2, mu).first).epsilon(2e-3));
    }
  }
}

// Check radiation integration across PDF flavour thresholds before differentiating the table
TEST_CASE("Sudakov radiation converges across flavour thresholds", "[gra::MSudakov][veto][convergence]") {
  const auto model = gra::MModelTune::Load(gra::ResolveModelDataFile("TUNE0", "GENERAL.json"));
  gra::MLHAPDFStore pdf_store;
  const auto pdf = pdf_store.GetPDF("MMHT2014lo68cl", 0);
  const double bottom2 = gra::math::pow2(pdf->info().get_entry_as<double>("MBottom"));
  std::array<gra::MSudakov, 2> sudakov;
  for (std::size_t i = 0; i < sudakov.size(); ++i) {
    auto numerics = model->Numerics();
    numerics["NUMERICS_SUDAKOV"]["SUDA"]["N_int"] = i == 0 ? 64 : 4096;
    sudakov[i].Init(13000.0, "MMHT2014lo68cl", model->Soft(), model->NumericsFile(), numerics.dump(), false, pdf);
  }
  for (const double mu : {30.0, 100.0}) {
    for (const double ratio : {0.51, 0.97, 1.03, 2.91}) {
      const double q2 = bottom2 * ratio;
      const auto coarse = sudakov[0].Sudakov_T(q2, mu);
      const auto fine   = sudakov[1].Sudakov_T(q2, mu);
      CAPTURE(mu, q2);
      REQUIRE(coarse.first == Approx(fine.first).epsilon(2e-5));
      REQUIRE(coarse.second == Approx(fine.second).epsilon(2e-5));
      const double h = 1e-4 * q2;
      const double derivative = (sudakov[1].Sudakov_T(q2 + h, mu).first -
                                   sudakov[1].Sudakov_T(q2 - h, mu).first) / (2.0 * h);
      REQUIRE(derivative == Approx(fine.second).epsilon(2e-5));
    }
  }
}
