// Tensor Pomeron model tests
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Kinematics/MFactorized.h"
#include "Graniitti/Regge/MReggeSource.h"
#include "Graniitti/Tensor/MTensorPhoto.h"
#include "Graniitti/Tech/MException.h"
#include "support/models_test_support.hh"

// Check rotations, parity, beam exchange and longitudinal boosts of physical amplitudes
template <typename Evaluate>
void RequireTensorSymmetries(const gra::LORENTZSCALAR &lts, const Evaluate &evaluate) {
  auto original = lts;
  original.amplitude.BeginCentral();
  const double amp2 = evaluate(original);
  REQUIRE(amp2 > 0.0);
  auto transformed = std::array{RotateToyEventAroundZ(lts, 0.83), ReflectToyEventInXZ(lts),
                                 BeamExchangeMirrorWithDecay(lts), BoostToyEventAlongZ(lts, 0.37)};
  for (const auto i : indices(transformed)) {
    CAPTURE(i);
    transformed[i].amplitude.BeginCentral();
    CHECK(evaluate(transformed[i]) == Approx(amp2).epsilon(3.0e-9));
  }
}

// Check exchange and baryon species select the shared transfer at every strong vertex
TEST_CASE("Tensor baryon transfer steering multiplies each complex vertex", "[MTensorPomeron][form-factor][dirac]") {
  const bool noflip = GENERATE(true, false);
  const std::string em = GENERATE("DIPOLE", "KELLY");
  for (const int baryon : {2212, 3122}) {
    for (const int exchange : {995, 9915, 9925, 9993, 9933, 9943}) {
      if (baryon == 3122 && exchange != 995) { continue; }
      CAPTURE(noflip, em, baryon, exchange);
      const std::string pair = "[" + std::to_string(baryon) + "," + std::to_string(baryon) + "]";
      const auto general = [&](auto &j) {
        j.at("PARAM_STRUCTURE").at("EM") = em == "DIPOLE" ? "KELLY" : "DIPOLE";
        auto &tp = j.at("PARAM_TENSORPOM");
        tp.at("FORWARD_NOFLIP") = noflip;
        tp.at("PARAM_CON").at("TP")["[" + std::to_string(baryon) + ",-" + std::to_string(baryon) + "]"] =
            {{exchange, exchange}};
      };
      const auto evaluate = [&](bool active) {
        const auto tune = WriteModifiedPhotoVMTune("tensor_baryon_dirac_" + std::to_string(active), general, [&](auto &j) {
          SetContinuumField(j, "[2212,2212]", "FF_transfer", {{"type", "none"}});
          SetContinuumField(j, pair, "FF_transfer", {{"type", "none"}});
          if (active) {
            j.at(std::to_string(exchange)).at(pair).at("FF_transfer") =
                {{"type", "dirac"}, {"norm", "zero"}, {"EM", em}};
          }
        });
        auto lts = DirectCentralPairLTSForTest(-baryon, baryon);
        gra::MTensorPomeron tensor(lts, gra::MModelTune::Load(tune.second),
            gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
        REQUIRE(tensor.ME4(lts, gra::TensorContinuumMode::TensorPomeron) > 0.0);
        return std::vector<std::complex<double>>(lts.hamp.begin(), lts.hamp.end());
      };
      const auto without_ff = evaluate(false);
      const auto with_ff = evaluate(true);
      const auto point = DirectCentralPairLTSForTest(-baryon, baryon);
      gra::form::ParamStore structure;
      structure.EM = em;
      double factor = gra::form::F1(point.t1, structure) * gra::form::F1(point.t2, structure);
      if (baryon == 2212) { factor *= factor; }
      REQUIRE(with_ff.size() == without_ff.size());
      for (const auto &i : indices(with_ff)) { RequireComplexNear(with_ff[i], factor * without_ff[i], 2.0e-11); }
    }
  }
}

// Check strong PHOTO vertices use the selected exchange form factor independently of photon vertices
TEST_CASE("Tensor PHOTO uses baryon transfer steering", "[MTensorPomeron][photo][form-factor][dirac]") {
  const bool noflip = GENERATE(true, false);
  const int exchange = GENERATE(995, 9915, 9925, 9933, 9943);
  const auto point = AsymmetricCentralPairLTSForTest(321, -321, 0.91, -0.38);
  for (const bool upper : {true, false}) {
    const auto evaluate = [&](bool active) {
      const auto tune = WriteModifiedPhotoVMTune("tensor_photo_dirac_" + std::to_string(active), [&](auto &j) {
        j.at("PARAM_STRUCTURE").at("EM") = "KELLY";
        auto &tp = j.at("PARAM_TENSORPOM");
        tp.at("FORWARD_NOFLIP") = noflip;
        tp.at("PHOTO").at("photon_exchange") = false;
        tp.at("PHOTO").at("exchanges") = {exchange};
      }, [&](auto &j) {
        SetContinuumField(j, "[2212,2212]", "FF_transfer", {{"type", "none"}});
        if (active) {
          j.at(std::to_string(exchange)).at("[2212,2212]").at("FF_transfer") =
              {{"type", "dirac"}, {"norm", "zero"}, {"EM", "DIPOLE"}};
        }
      });
      auto lts = point;
      const auto model = gra::MModelTune::Load(tune.second);
      gra::MTensorPomeron tensor(lts, model,
          gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Photo), gra::MTensorPomeronMode::Photo);
      const gra::MTensorPhoto photo(tensor, gra::ReadTensorPomeronParam(*model, lts.PDG));
      return photo.DSCurrent(lts, upper, 0, 0);
    };
    const auto without_ff = evaluate(false);
    const auto with_ff = evaluate(true);
    gra::form::ParamStore structure;
    structure.EM = "DIPOLE";
    const double factor = gra::form::F1(upper ? point.t2 : point.t1, structure);
    REQUIRE(gra::SquaredNorm(without_ff) > 0.0);
    for (const auto &mu : indices(with_ff)) { RequireComplexNear(with_ff[mu], factor * without_ff[mu], 2.0e-11); }
  }
}

// Compare the PHOTO pole terms to direct unequal-subenergy Regge powers
// [REFERENCE: Lebiedowicz, Nachtmann and Szczurek, arXiv:2508.06334v2, Eqs. (2.21)-(2.23), (2.47)-(2.49)]
// [REFERENCE: Lebiedowicz et al., arXiv:2609.01285, Eqs. (2.15)-(2.18)]
TEST_CASE("Tensor PHOTO retains finite unequal subenergies", "[MTensorPomeron][photo][physics][numerics]") {
  for (const int point : {0, 1, 2}) {
    const bool asymmetric = point == 2;
    auto lts = asymmetric ? ScalarPolePhasePointForTest(0.02, 1.0e-5, 0.3, 1.0e4, 1.0e6)
                          : ScalarPolePhasePointForTest(0.2, point == 1 ? gra::math::PI / 2.0 : 0.8,
                                                       point == 1 ? gra::math::PI / 2.0 : 0.3, 1.3);
    const auto &q = lts.q1;
    const auto p = lts.pbeam2 + lts.pfinal[2];
    const auto &kp = lts.decaytree[0].p4;
    const auto &km = lts.decaytree[1].p4;
    const double t = (lts.pfinal[2] - lts.pbeam2).M2();
    const auto extended_p = p.Contravariant<long double>();
    const long double nu1 = 0.5L * gra::MinkowskiProduct(extended_p, kp.Contravariant<long double>());
    const long double nu2 = 0.5L * gra::MinkowskiProduct(extended_p, km.Contravariant<long double>());
    const long double nubar2 = 0.5L * (nu1 * nu1 + nu2 * nu2);
    // Contract perpendicular to the target current while retaining the off-shell contact term
    const gra::M4Vec n(p.Py(), -p.Px(), 0.0, 0.0);
    for (const int pdg : {995, 9933}) {
      for (const double exponent : {0.45, 9.0e-7, -9.0e-7, 0.0}) {
        CAPTURE(point, pdg, exponent);
        const auto tune = WriteModifiedPhotoVMTune("tensor_photo_subenergy_" + std::to_string(point) + "_" +
            std::to_string(pdg) + "_" + std::to_string(exponent), [&](auto &j) {
          auto &tp = j.at("PARAM_TENSORPOM");
          tp.at("FORWARD_NOFLIP") = true;
          tp.at("PHOTO").at("photon_exchange") = false;
          tp.at("PHOTO").at("exchanges") = {pdg};
          if (asymmetric) {
            tp.at("PHOTO").at("offshell").at("211") = {{"FF_transfer", {{"type", "none"}}}, {"FF_offshell", {{"type", "none"}}}};
          }
          auto &exchange = tp.at("EXCHANGES").at(std::to_string(pdg));
          const double slope = exchange.at("ap").template get<double>();
          exchange.at("delta") = (pdg == 995 ? 1.0 : 0.0) - 2.0 * exponent - slope * t;
        });
        const auto model = gra::MModelTune::Load(tune.second);
        auto event = lts;
        gra::MTensorPomeron tensor(event, model,
            gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Photo), gra::MTensorPomeronMode::Photo);
        const auto parameters = gra::ReadTensorPomeronParam(*model, lts.PDG);
        const gra::MTensorPhoto photo(tensor, parameters);
        const auto current = photo.DSCurrent(event, true, 0, 0);
        const auto &exchange = parameters->exchange.FindExchange(pdg);
        const double alpha = 1.0 + exchange.delta + exchange.ap * t;
        const long double lambda = 0.5 * ((pdg == 995 ? 2.0 : 1.0) - alpha);
        const long double wa = std::pow(nu2 * nu2 / nubar2, -lambda);
        const long double wb = std::pow(nu1 * nu1 / nubar2, -lambda);
        const double fm = parameters->photo.m0_2 / (parameters->photo.m0_2 - q.M2());
        const long double ha = fm * (n * (kp * 2.0 - q)) / (q.M2() - 2.0 * (kp * q)) +
                               (n * q) / (parameters->photo.m0_2 - q.M2());
        const long double hb = fm * (n * (km * 2.0 - q)) / (q.M2() - 2.0 * (km * q)) +
                               (n * q) / (parameters->photo.m0_2 - q.M2());
        const auto rp = km - kp + q;
        const auto rm = km - kp - q;
        const long double a = gra::MinkowskiProduct(extended_p, rp.Contravariant<long double>());
        const long double b = gra::MinkowskiProduct(extended_p, rm.Contravariant<long double>());
        const double pion_g = parameters->exchange.FindVertex(pdg, 211).g_tensor.front();
        const double proton_g = parameters->exchange.FindVertex(pdg, 2212).g_tensor.front();
        const auto propagator = parameters->exchange.PropagatorFactor(pdg, 2.0 * std::sqrt(static_cast<double>(nubar2)), t);
        const double form = gra::regge::TransferFF(t, parameters->exchange.FindVertex(pdg, 2212).ff_transfer);
        const auto &[transfer, offshell] = parameters->photo.offshell.at(211);
        const long double ff = gra::regge::TransferFF(t, transfer);
        const long double xa = q.M2() - 2.0L * (q * kp), xb = q.M2() - 2.0L * (q * km);
        const bool exponential = offshell.type == gra::regge::FFType::Exponential;
        const bool gaussian = offshell.type == gra::regge::FFType::Gaussian;
        const long double slope = exponential ? offshell.param[0] : 0.0;
        const long double aoff = gaussian ? -1.0L / pow2(offshell.param[0]) : 0.0;
        const long double ga = ff * std::exp(slope * xa + aoff * xa * xa);
        const long double gb = ff * std::exp(slope * xb + aoff * xb * xb);
        const long double h = std::abs(xa - xb) < 1.0e-12L ? ga * (slope + 2.0L * aoff * xa) : (ga - gb) / (xa - xb);
        std::complex<double> residue;
        long double expected;
        if (pdg == 995) {
          residue = gra::math::zi * 6.0 * pion_g * proton_g * form * propagator;
          const long double sa = wa * (2.0L * a * a - 2.0L * pow2(gra::PDG::mp) * rp.M2());
          const long double sb = wb * (2.0L * b * b - 2.0L * pow2(gra::PDG::mp) * rm.M2());
          expected = ga * ha * sa - gb * hb * sb -
                     4.0L * (ga + gb) * pow2(gra::PDG::mp) * (n * (km - kp)) +
                     (n * (km - kp)) * h * (sa + sb);
        } else {
          residue = -gra::math::zi * 0.5 * pion_g * proton_g * form * propagator;
          expected = ga * ha * wa * a - gb * hb * wb * b + (n * (km - kp)) * h * (wa * a + wb * b);
        }
        std::complex<double> observed = 0.0;
        for (const auto &mu : indices(current)) { observed += n[mu] * current[mu]; }
        observed /= gra::qed::e_QED() * residue;
        CAPTURE(expected, observed, n[0], n[1], n[2], n[3], current);
        double scale = std::abs(static_cast<double>(expected));
        for (const auto &mu : indices(current)) { scale += std::abs(current[mu] / (gra::qed::e_QED() * residue)); }
        CHECK(std::abs(observed - static_cast<double>(expected)) < 2.0e-9 * scale);
        if (asymmetric) { continue; }

        // Compare every complex current component using direct Regge differences for the contact term
        const auto d = km - kp;
        const long double dp = gra::MinkowskiProduct(extended_p, d.Contravariant<long double>());
        const long double qp = gra::MinkowskiProduct(extended_p, q.Contravariant<long double>());
        const long double sa = pdg == 995 ? 2.0L * a * a - 2.0L * pow2(gra::PDG::mp) * rp.M2() : a;
        const long double sb = pdg == 995 ? 2.0L * b * b - 2.0L * pow2(gra::PDG::mp) * rm.M2() : b;
        for (const auto &mu : indices(current)) {
          const long double pole_a = fm * ((kp * 2.0 - q) % mu) / (q.M2() - 2.0 * (kp * q)) +
                                     (q % mu) / (parameters->photo.m0_2 - q.M2());
          const long double pole_b = fm * ((km * 2.0 - q) % mu) / (q.M2() - 2.0 * (km * q)) +
                                     (q % mu) / (parameters->photo.m0_2 - q.M2());
          const long double contact = (pdg == 995 ? 8.0L * (dp * (p % mu) - pow2(gra::PDG::mp) * (d % mu))
                                                        : 2.0L * (p % mu)) +
                                      (p % mu) / qp * ((wa - 1.0L) * sa - (wb - 1.0L) * sb);
          const long double direct = ga * pole_a * wa * sa - gb * pole_b * wb * sb +
                                     (d % mu) * h * (wa * sa + wb * sb) + 0.5L * (ga + gb) * contact;
          const auto reference = gra::qed::e_QED() * residue * static_cast<double>(direct);
          CAPTURE(mu, reference, current[mu]);
          CHECK(std::abs(current[mu] - reference) < 2.0e-9 * std::abs(gra::qed::e_QED() * residue) * scale);
        }
      }
    }
  }
}

// Check Ward identities and charge conjugation for each isolated PHOTO exchange
TEST_CASE("Tensor PHOTO exchanges obey Ward identities and C parity", "[MTensorPomeron][photo][physics][gauge]") {
  const int meson = GENERATE(211, 321);
  const std::string form = GENERATE("none", "exp", "gaussian");
  CAPTURE(meson, form);
  for (const bool noflip : {false, true}) {
    for (const int pdg : {995, 9915, 9933, 22}) {
      const auto tune = WriteModifiedPhotoVMTune("tensor_photo_ward_" + std::to_string(noflip) + "_" + std::to_string(pdg),
          [&](auto &j) {
            auto &tp = j.at("PARAM_TENSORPOM");
            tp.at("FORWARD_NOFLIP") = noflip;
            tp.at("PHOTO").at("photon_exchange") = pdg == 22;
            tp.at("PHOTO").at("exchanges") = {pdg == 22 ? 995 : pdg};
            auto &ff = tp.at("PHOTO").at("offshell").at(std::to_string(meson)).at("FF_offshell");
            ff = {{"type", form}};
            if (form != "none") {
              ff["norm"] = "pole";
              ff[form == "exp" ? "b" : "Lambda2"] = 1.3;
            }
          }, [&](auto &j) {
            if (pdg == 22) { j.at("995").at("[" + std::to_string(meson) + "," + std::to_string(meson) + "]").at("g_tensor") = {0.0}; }
          });
      auto lts = AsymmetricCentralPairLTSForTest(meson, -meson, 0.91, -0.38);
      const auto model = gra::MModelTune::Load(tune.second);
      gra::MTensorPomeron tensor(lts, model,
          gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Photo), gra::MTensorPomeronMode::Photo);
      const gra::MTensorPhoto photo(tensor, gra::ReadTensorPomeronParam(*model, lts.PDG));
      for (const bool upper : {true, false}) {
        if (pdg == 22 && !upper) { continue; }
        const auto &q = upper ? lts.q1 : lts.q2;
        for (const auto initial : gra::spin::BinaryHelicityIndices()) {
          for (const auto final : gra::spin::BinaryHelicityIndices()) {
            if (noflip && initial != final) { continue; }
            const auto current = photo.DSCurrent(lts, upper, initial, final);
            std::complex<double> ward = 0.0;
            double scale = 0.0;
            for (const auto &mu : indices(current)) {
              ward += q[mu] * current[mu];
              scale += std::abs(q[mu] * current[mu]);
            }
            REQUIRE(scale > 0.0);
            CHECK(std::abs(ward) < 2.0e-9 * scale);
            std::swap(lts.decaytree[0].p4, lts.decaytree[1].p4);
            const auto crossed = photo.DSCurrent(lts, upper, initial, final);
            std::swap(lts.decaytree[0].p4, lts.decaytree[1].p4);
            const double parity = pdg == 995 || pdg == 9915 ? -1.0 : 1.0;
            for (const auto &mu : indices(current)) { RequireComplexNear(crossed[mu], parity * current[mu], 2.0e-9); }
          }
        }
      }
    }
  }
}

// Check event-dependent Regge failures use the sampling exception type
TEST_CASE("Tensor Regge propagators report invalid kinematics", "[MTensorPomeron][exchange][failure]") {
  const auto parameters = gra::ReadTensorPomeronParam(*gra::MModelTune::Load(modelfile), LoadedPDGTable());
  for (const int pdg : {995, 9933, 9993}) {
    for (const double s : {0.0, -1.0, std::numeric_limits<double>::infinity(),
                           std::numeric_limits<double>::quiet_NaN()}) {
      CHECK_THROWS_AS(parameters->exchange.PropagatorFactor(pdg, s, -0.1), gra::AmplitudeFailure);
    }
    CHECK_THROWS_AS(parameters->exchange.PropagatorFactor(pdg, 100.0,
        std::numeric_limits<double>::quiet_NaN()), gra::AmplitudeFailure);
    CHECK(std::isfinite(std::abs(parameters->exchange.PropagatorFactor(pdg, 100.0, -0.1))));
  }
  CHECK_THROWS_AS(parameters->exchange.PropagatorFactor(123456, 100.0, -0.1), std::invalid_argument);
}

// Check vector propagators distinguish invalid event momenta from invalid pole inputs
TEST_CASE("Tensor vector propagators report invalid kinematics", "[MTensorPomeron][vector][failure]") {
  auto lts = TensorRhoCascadeLTSForTest();
  const auto model = gra::MModelTune::Load(modelfile);
  gra::MTensorPomeron tensor(lts, model,
      gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
  const auto parameters = gra::ReadTensorPomeronParam(*model, lts.PDG);
  const auto &rho = parameters->FindVector(113);
  for (const double value : {std::numeric_limits<double>::infinity(),
                             std::numeric_limits<double>::quiet_NaN()}) {
    const gra::M4Vec invalid(0.0, 0.0, 0.0, value);
    CHECK_THROWS_AS(tensor.iD_V(invalid, rho.mass, 10.0, 113), gra::AmplitudeFailure);
    CHECK_THROWS_AS(tensor.iD_V(lts.q1, rho.mass, value, 113), gra::AmplitudeFailure);
    CHECK_THROWS_AS(tensor.iD_VMES(invalid, rho.mass, rho.width, 113, true, true), gra::AmplitudeFailure);
  }
  CHECK_THROWS_AS(tensor.iD_V(lts.q1, rho.mass, 0.0, 113), gra::AmplitudeFailure);
}

// Check PHOTO initialization and evaluation permit an active resonance without continuum
TEST_CASE("Tensor PHOTO accepts resonance production with zero continuum", "[MTensorPomeron][photo][process]") {
  const auto tune = WriteModifiedPhotoVMTune("tensor_photo_resonance_only", [](auto &j) {
    j["PARAM_TENSORPOM"]["PHOTO"]["photon_exchange"] = false;
  }, [](auto &j) {
    for (auto &[key, block] : j.items()) {
      (void)key;
      if (block.contains("[211,211]")) { block["[211,211]"]["g_tensor"] = {0.0}; }
    }
  });
  ToyHelicityProcess reader;
  reader.SetTuneForTest(tune.first);
  reader.SetProcessForTest("TP", "PHOTO");
  auto lts = AsymmetricCentralPairLTSForTest(211, -211, 0.91, -0.38);
  // Tensor PHOTO has no MP production frame
  lts.process.MP_FRAME = "null";
  gra::PARAM_RES res;
  res.p = lts.PDG.FindByPDG(113);
  SetToyTensorChannel(res, {0.43, -0.18});
  lts.process.RESONANCES = {{"rho", res}};
  const auto model = gra::MModelTune::Load(tune.second);
  double symmetry = 1.0;
  gra::MProcessSetup setup(lts, model, "TP", "PHOTO", gra::PROC_403_TP_PHOTO().Info(), false, false, 0, false, false, false, symmetry,
      [&](const auto &particle, const auto &daughters, bool production, bool strict,
          const auto &label, bool verbose, gra::spin::VertexContext context) {
        return reader.ProcessHelicityStructure(particle, daughters, production, strict, label, verbose, context);
      }, {}, {});
  gra::MTensorPomeron::InitializeParameters(setup);
  REQUIRE_NOTHROW(gra::MTensorPomeron::InitializeBranching(setup, gra::MTensorPomeronMode::Photo));
  gra::MTensorPomeron tensor(lts, model,
      gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Photo), gra::MTensorPomeronMode::Photo);
  REQUIRE(tensor.MEPhoto(lts) > 0.0);
  const auto amplitude = lts.hamp;
  tensor.ME3(lts, gra::MTensorPomeronMode::Photo);
  RequireVectorNear(lts.hamp, amplitude, 2.0e-10);

  const auto disabled = WriteModifiedPhotoVMTune("tensor_photo_disabled_resonance_source", [](auto &j) {
    auto &soft = j.at("PARAM_SOFT");
    soft.at("MODEL").at(soft.at("active_model").template get<std::string>()).at("EXCHANGE").at("R_f2").at("on") = false;
  });
  setup.soft_model = gra::MModelTune::Load(disabled.second)->Soft();
  lts.process.RESONANCES.at("rho").TP.channels.front().exchange = {22, 9915};
  CHECK_THROWS_AS(gra::MTensorPomeron::InitializeBranching(setup, gra::MTensorPomeronMode::Photo), std::invalid_argument);
  setup.soft_model = model->Soft();
  lts.process.RESONANCES.clear();
  CHECK_THROWS_AS(gra::MTensorPomeron::InitializeBranching(setup, gra::MTensorPomeronMode::Photo), std::invalid_argument);
}

// Check the full pion continuum crossing phase and charge-order independence
TEST_CASE("Tensor pion continuum obeys exchange C parity", "[MTensorPomeron][continuum][physics][crossing]") {
  for (const bool noflip : {true, false}) {
    for (const int second : {995, 9933}) {
      const auto tune = WriteModifiedPhotoVMTune("tensor_pion_crossing_" + std::to_string(noflip) + "_" + std::to_string(second),
          [=](auto &j) {
            j["PARAM_TENSORPOM"]["FORWARD_NOFLIP"] = noflip;
            j["PARAM_TENSORPOM"]["PARAM_CON"]["TP"]["[211,-211]"] = {{995, second}};
          });
      auto lts = AsymmetricCentralPairLTSForTest(211, -211, 0.91, -0.38);
      gra::MTensorPomeron tensor(lts, gra::MModelTune::Load(tune.second),
          gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
      REQUIRE(tensor.ME4(lts, gra::TensorContinuumMode::TensorPomeron) > 0.0);
      const auto original = lts.hamp;
      std::swap(lts.decaytree[0], lts.decaytree[1]);
      tensor.ME4(lts, gra::TensorContinuumMode::TensorPomeron);
      RequireVectorNear(lts.hamp, original, 2.0e-10);
      std::swap(lts.decaytree[0].p4, lts.decaytree[1].p4);
      tensor.ME4(lts, gra::TensorContinuumMode::TensorPomeron);
      const double parity = second == 995 ? 1.0 : -1.0;
      for (const auto &i : indices(original)) { RequireComplexNear(lts.hamp[i], parity * original[i], 2.0e-10); }
    }
  }
}

// Check the coherent PHOTO weight in a common collider and pion charge convention
TEST_CASE("Tensor PHOTO preserves rotations and pion branch order", "[MTensorPomeron][photo][physics][phase]") {
  for (const bool noflip : {true, false}) {
    for (const bool pauli : {true, false}) {
      const auto tune = WriteModifiedPhotoVMTune("tensor_photo_invariance_" + std::to_string(noflip) + "_" + std::to_string(pauli),
          [=](auto &j) {
            j["PARAM_TENSORPOM"]["FORWARD_NOFLIP"] = noflip;
            j["PARAM_TENSORPOM"]["PHOTO"]["proton_pauli"] = pauli;
            j["PARAM_TENSORPOM"]["PHOTO"]["photon_exchange"] = true;
          });
      auto lts = AsymmetricCentralPairLTSForTest(211, -211, 0.91, -0.38);
      gra::PARAM_RES res;
      res.p = lts.PDG.FindByPDG(113);
      SetToyTensorChannel(res, {0.43, -0.18});
      res.hel_decay.g_decay_TP = {0.71};
      lts.process.RESONANCES = {{"rho", res}};
      gra::MTensorPomeron tensor(lts, gra::MModelTune::Load(tune.second),
          gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Photo), gra::MTensorPomeronMode::Photo);
      const double original = tensor.MEPhoto(lts);
      REQUIRE(original > 0.0);
      const auto amplitude = lts.hamp;
      auto rotated = RotateToyEventAroundZ(lts, 0.83);
      CHECK(tensor.MEPhoto(rotated) == Approx(original).epsilon(2.0e-10));
      std::swap(lts.decaytree[0], lts.decaytree[1]);
      CHECK(tensor.MEPhoto(lts) == Approx(original).epsilon(2.0e-10));
      RequireVectorNear(lts.hamp, amplitude, 2.0e-10);
    }
  }
}

// Check a single gamma-gamma contribution with symmetric virtuality cuts
TEST_CASE("Tensor PHOTO gamma-gamma cuts preserve beam exchange", "[MTensorPomeron][photo][physics][beams]") {
  const auto tune = WriteModifiedPhotoVMTune("tensor_photo_two_photon_cut", [](auto &j) {
    j["PARAM_TENSORPOM"]["FORWARD_NOFLIP"] = true;
    j["PARAM_TENSORPOM"]["PHOTO"]["photon_exchange"] = true;
  }, [](auto &j) {
    for (auto &[key, block] : j.items()) {
      (void)key;
      if (block.contains("[211,211]")) { block["[211,211]"]["g_tensor"] = {0.0}; }
    }
  });
  auto lts = AsymmetricCentralPairLTSForTest(211, -211, 0.91, -0.38);
  const auto model = gra::MModelTune::Load(tune.second);
  gra::MTensorPomeron tensor(lts, model,
      gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Photo), gra::MTensorPomeronMode::Photo);
  auto parameters = std::make_shared<gra::MTensorPomeronParam>(*gra::ReadTensorPomeronParam(*model, lts.PDG));
  auto exchanged = lts;
  const auto rotate = [](const gra::M4Vec &p) { return gra::M4Vec(p.Px(), -p.Py(), -p.Pz(), p.E()); };
  exchanged.pbeam1 = rotate(lts.pbeam2);
  exchanged.pbeam2 = rotate(lts.pbeam1);
  for (auto &p : exchanged.pfinal) { p = rotate(p); }
  std::swap(exchanged.pfinal[1], exchanged.pfinal[2]);
  for (auto &branch : exchanged.decaytree) { branch.p4 = rotate(branch.p4); }
  RefreshToyDerivedKinematicsPreserveDecay(exchanged);
  for (const bool noflip : {true, false}) {
    for (const bool pauli : {true, false}) {
      parameters->FORWARD_NOFLIP = noflip;
      parameters->photo.proton_pauli = pauli;
      for (const bool restricted : {false, true}) {
        CAPTURE(noflip, pauli, restricted);
        parameters->photo.q2_max = restricted ? -0.5 * (lts.t1 + lts.t2) : 2.0 * std::max(-lts.t1, -lts.t2);
        const gra::MTensorPhoto photo(tensor, parameters);
        const double original = photo.Amp2(lts);
        const double swapped = photo.Amp2(exchanged);
        if (restricted) {
          CHECK(std::abs(original) < 1.0e-20);
          CHECK(std::abs(swapped) < 1.0e-20);
        } else {
          REQUIRE(original > 0.0);
          CHECK(swapped == Approx(original).epsilon(2.0e-10));
        }
      }
    }
  }
}

// Check event-dependent nonfinite PHOTO currents reach amplitude bookkeeping
TEST_CASE("Tensor PHOTO reports nonfinite amplitudes", "[MTensorPomeron][photo][failure]") {
  auto lts = AsymmetricCentralPairLTSForTest(211, -211, 0.91, -0.38);
  const auto model = gra::MModelTune::Load(modelfile);
  gra::MTensorPomeron tensor(lts, model,
      gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Photo), gra::MTensorPomeronMode::Photo);
  const gra::MTensorPhoto photo(tensor, gra::ReadTensorPomeronParam(*model, lts.PDG));
  lts.decaytree[0].p4 = gra::M4Vec(0.0, 0.0, 0.0, std::numeric_limits<double>::quiet_NaN());
  CHECK_THROWS_AS(photo.DSCurrent(lts, true, 0, 0), gra::AmplitudeFailure);
}

// Check overflow of finite Regge couplings is recorded as an amplitude failure
TEST_CASE("Tensor PHOTO reports Regge residue overflow", "[MTensorPomeron][photo][failure]") {
  const auto tune = WriteModifiedPhotoVMTune("tensor_photo_residue_overflow", [](auto &) {}, [](auto &j) {
    j["995"]["[211,211]"]["g_tensor"] = {1.0e200};
    j["995"]["[2212,2212]"]["g_tensor"] = {1.0e200};
  });
  auto lts = AsymmetricCentralPairLTSForTest(211, -211, 0.91, -0.38);
  const auto model = gra::MModelTune::Load(tune.second);
  gra::MTensorPomeron tensor(lts, model,
      gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Photo), gra::MTensorPomeronMode::Photo);
  const gra::MTensorPhoto photo(tensor, gra::ReadTensorPomeronParam(*model, lts.PDG));
  CHECK_FALSE(std::isfinite(photo.Amp2(lts)));

  // Exercise the shared amplitude failure boundary with the real Tensor process
  class TensorWeight : public gra::MFactorized {
   public:
    using gra::MFactorized::MFactorized;
    using gra::MProcess::GetAmp2;
  };
  TensorWeight process("TP[PHOTO]<F>", {}, model);
  process.state.lts = lts;
  process.state.lts.process.MP_FRAME = "null";
  process.SetModelTune(model);
  process.SetHelicityConfig(model);
  process.SetScreening(false);
  process.InitializeProcessAmplitude();
  gra::MEventWeightState aux;
  CHECK(gra::math::IsZero(process.GetAmp2(false, aux)));
  CHECK_FALSE(aux.amplitude_ok);
  CHECK(aux.technical_failure);
}

// Check daughter form factors enter the coherent vector cascade amplitude
TEST_CASE("Tensor cascades apply branch decay form factors", "[MTensorPomeron][cascade][form-factor]") {
  auto lts = TensorRhoCascadeLTSForTest();
  gra::MTensorPomeron tensor(lts, gra::MModelTune::Load(modelfile),
      gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
  const double original = tensor.ME6(lts);
  REQUIRE(original > 0.0);
  for (auto &branch : lts.decaytree) {
    branch.hel.ff_decay = gra::regge::ReadFF(
        nlohmann::json{{"type", "gaussian"}, {"norm", "pole"}, {"Lambda2", 0.1}}, "tensor cascade decay test");
  }
  const double dressed = tensor.ME6(lts);
  REQUIRE(dressed > 0.0);
  CHECK(std::abs(dressed - original) > 1.0e-3 * original);
  const auto amplitude = lts.hamp;
  std::swap(lts.decaytree[0].legs[0], lts.decaytree[0].legs[1]);
  tensor.ME6(lts);
  RequireVectorNear(lts.hamp, amplitude, 2.0e-10);
}

namespace {

// Compute the canonical Tensor hard-spin layout used by pair screening
gra::ScreeningMetadata TensorPairMetadata(const bool noflip) {
  gra::ScreeningMetadata metadata;
  metadata.spin_basis            = gra::ScreeningSpinBasis::ProtonHelicity;
  metadata.proton_mode           = gra::ProtonScreeningMode::ForwardExcitation;
  metadata.amplitude_normalization = 0.25;
  metadata.spin_rows             = noflip ? 4 : 16;
  metadata.forward_noflip        = noflip;
  metadata.amplitude_type       = gra::ScreeningAmplitudeType::GoodWalker;
  return metadata;
}

// Require the canonical pair source to reproduce the physical Tensor Born row
void RequireTensorBornProjection(const gra::LORENTZSCALAR &lts, const std::size_t channels) {
  REQUIRE(lts.proton_good_walker.has_value());
  const auto &pair = *lts.proton_good_walker;
  REQUIRE(pair.channel_count == channels);
  REQUIRE(pair.components.size() == 1);
  RequireVectorNear(gra::ProjectGoodWalker(pair), lts.hamp, 2.0e-10);
}

// Require one equal Tensor pair to carry normalized SOFT vertex directions
void RequireTensorResiduePair(const gra::LORENTZSCALAR &lts, const gra::SoftModelPtr &model,
                              const std::string &exchange_name) {
  REQUIRE(lts.proton_good_walker.has_value());
  REQUIRE(lts.proton_good_walker->components.size() == 1);
  const auto                       &source = lts.proton_good_walker->components.front().source;
  std::vector<std::complex<double>> proton(model->GoodWalker().ProtonVector().begin(),
                                           model->GoodWalker().ProtonVector().end());
  const auto                        exchange = model->ExchangeId(exchange_name);

  // Normalize one public SOFT vertex against the physical proton
  const auto leg_source = [&](const double t) {
    auto                       leg      = model->ResidueMatrix(exchange, t) * proton;
    const std::complex<double> physical = gra::InnerProduct(proton, leg);
    REQUIRE(std::abs(physical) > 1.0e-14);
    gra::Scale(leg, 1.0 / physical);
    return leg;
  };
  const auto expected = gra::KroneckerProduct(leg_source(lts.t1), leg_source(lts.t2));
  REQUIRE(source.size_row() == lts.hamp.size());
  REQUIRE(source.size_col() == expected.size());
  for (const auto &row : indices(lts.hamp)) {
    for (const auto &col : indices(expected)) {
      RequireComplexNear(source(row, col), lts.hamp[row] * expected[col], 2.0e-10);
    }
  }
}

// Run one real Tensor continuum through the production screening interface
class ToyTensorScreeningProcess : public gra::MFactorized {
 public:
  // Construct one elastic Tensor process in the canonical pair route
  explicit ToyTensorScreeningProcess(const gra::MModelTunePtr &model_tune) {
    state.screening = true;
    ProcPtr.ISTATE  = "TP";
    ProcPtr.CHANNEL = "CON";
    state.lts       = DirectCentralPairLTSForTest(211, -211);
    state.lts.hamp.Configure(TensorPairMetadata(true));
    SetModelTune(model_tune);
    PrepareScreeningPoint(state);
    tensor = std::make_unique<gra::MTensorPomeron>(
        state.lts, model_tune, gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
  }

  // Evaluate only the physical Tensor Born amplitude
  double BornAmp2() {
    pair_trace.clear();
    state.lts.proton_good_walker.reset();
    return ScreenedAmplitudeSquared(false);
  }

  // Evaluate the Tensor amplitude with its coupled channel convolution
  double ScreenedAmp2() {
    pair_trace.clear();
    state.lts.proton_good_walker.reset();
    return ScreenedAmplitudeSquared(true);
  }

  // Compute the Born source followed by every shifted loop source
  const std::vector<gra::ProtonGoodWalkerAmplitude> &PairTrace() const { return pair_trace; }

 protected:
  // Evaluate and retain one channel-resolved Tensor source
  double EvaluateBareAmplitude() override {
    const double amp2 = tensor->ME4(state.lts, gra::TensorContinuumMode::TensorPomeron);
    REQUIRE(state.lts.proton_good_walker.has_value());
    pair_trace.push_back(*state.lts.proton_good_walker);
    return amp2;
  }

 private:
  std::unique_ptr<gra::MTensorPomeron> tensor;
  std::vector<gra::ProtonGoodWalkerAmplitude> pair_trace;
};

// Contract a traced Tensor source with the public pair screening operators
std::vector<std::complex<double>> TensorScreeningReference(const std::vector<gra::ProtonGoodWalkerAmplitude> &trace,
                                                           const gra::MEikonal                        &eikonal) {
  REQUIRE_FALSE(trace.empty());
  const auto &born = trace.front();
  REQUIRE(born.model == eikonal.SoftModelHandle());
  REQUIRE(born.channel_count == eikonal.GetChannelCount());
  REQUIRE(born.components.size() == 1);
  gra::ScreeningMetadata metadata;
  metadata.spin_basis     = gra::ScreeningSpinBasis::ProtonHelicity;
  metadata.spin_rows      = 4;
  metadata.forward_noflip = true;
  metadata.PrepareSpinTransitions();
  REQUIRE(metadata.spin_basis == gra::ScreeningSpinBasis::ProtonHelicity);
  REQUIRE(metadata.forward_noflip);
  REQUIRE(metadata.spin_rows == 4);
  const std::size_t d           = born.model->GoodWalker().PairDimension();
  const auto       &born_source = born.components.front().source;
  REQUIRE(born_source.size_col() == d);
  REQUIRE(born_source.size_row() % metadata.spin_rows == 0);
  const std::size_t spectators = born_source.size_row() / metadata.spin_rows;
  const auto       &loop       = eikonal.GetLoopConst(eikonal.InitializedMandelstamS());
  const std::size_t n_phi      = loop.node_weight.size_col();
  const std::size_t nodes      = loop.kt2.size() * n_phi;
  REQUIRE(trace.size() == nodes + 1);

  std::vector<std::complex<double>> pair_amplitude(16 * spectators * d, 0.0);
  for (std::size_t index = 0; index < metadata.spin_transition_count; ++index) {
    const auto &transition = metadata.spin_transition[index];
    for (std::size_t spectator = 0; spectator < spectators; ++spectator) {
      const std::size_t source = transition.source_row * spectators + spectator;
      const std::size_t output =
          (gra::spin::PairHelicityTransitionIndex(transition.initial, transition.intermediate) * spectators +
           spectator) *
          d;
      for (std::size_t a = 0; a < d; ++a) { pair_amplitude[output + a] += born_source[source][a]; }
    }
  }

  constexpr auto spin_transition = gra::CanonicalProtonHelicityTransitions();
  for (std::size_t node = 0; node < nodes; ++node) {
    const std::size_t radial  = node / n_phi;
    const std::size_t azimuth = node % n_phi;
    const auto       &shifted = trace[node + 1];
    REQUIRE(shifted.model == born.model);
    REQUIRE(shifted.channel_count == born.channel_count);
    REQUIRE(shifted.components.size() == 1);
    const auto &source = shifted.components.front().source;
    REQUIRE(source.size_row() == born_source.size_row());
    REQUIRE(source.size_col() == d);
    const auto  &soft = loop.pair_screening_spin[radial];
    const double phi  = std::atan2(loop.kt_y(radial, azimuth), loop.kt_x(radial, azimuth));
    for (std::size_t index = 0; index < metadata.spin_transition_count; ++index) {
      const auto &transition = metadata.spin_transition[index];
      for (std::size_t final = 0; final < 4; ++final) {
        const std::size_t          spin = gra::spin::PairHelicityMatrixIndex(final, transition.intermediate);
        const std::complex<double> weight =
            loop.node_weight(radial, azimuth) * std::polar(1.0, spin_transition[spin].azimuth_harmonic * phi);
        for (std::size_t spectator = 0; spectator < spectators; ++spectator) {
          const std::size_t source_row = transition.source_row * spectators + spectator;
          const std::size_t output =
              (gra::spin::PairHelicityTransitionIndex(transition.initial, final) * spectators + spectator) * d;
          for (std::size_t a = 0; a < d; ++a) {
            for (std::size_t b = 0; b < d; ++b) {
              pair_amplitude[output + a] += weight * soft[spin][a][b] * source[source_row][b];
            }
          }
        }
      }
    }
  }

  std::vector<std::complex<double>> out;
  for (std::size_t row = 0; row < 16 * spectators; ++row) {
    const std::span<const std::complex<double>> source(pair_amplitude.data() + row * d, d);
    const auto projected = born.model->GoodWalker().ProjectPair(source, gra::GoodWalkerFinalBasis::Proton,
                                                                gra::GoodWalkerFinalBasis::Proton);
    REQUIRE(projected.size() == 1);
    out.push_back(projected.front());
  }
  return out;
}

// Evaluate the unreduced direct Tensor Pomeron continuum contractions
std::vector<std::complex<double>> TensorContinuumReferenceAmplitudes(gra::MTensorPomeron            &tensor,
                                                                     const gra::LORENTZSCALAR       &lts,
                                                                     const gra::MTensorPomeronParam &parameters) {
  using FTensor::Tensor2;
  using FTensor::Tensor4;

  const int   spin_x2        = lts.decaytree[0].p.spinX2;
  const int   pdg_left       = lts.decaytree[0].p.pdg;
  const int   abs_pdg        = std::abs(pdg_left);
  std::size_t anti_index     = 0;
  std::size_t particle_index = 1;
  if (spin_x2 == 1 && pdg_left > 0) {
    anti_index     = 1;
    particle_index = 0;
  }

  const gra::M4Vec pa             = lts.pbeam1;
  const gra::M4Vec pb             = lts.pbeam2;
  const gra::M4Vec p1             = lts.pfinal[1];
  const gra::M4Vec p2             = lts.pfinal[2];
  const gra::M4Vec p3             = lts.decaytree[anti_index].p4;
  const gra::M4Vec p4             = lts.decaytree[particle_index].p4;
  const gra::M4Vec pt             = pa - p1 - p3;
  const gra::M4Vec pu             = p4 - pa + p1;
  const double     mass           = lts.decaytree[particle_index].p.mass;
  const auto       upper_state    = gra::ResolveForwardLegState(lts, gra::ForwardBeamLeg::Upper);
  const auto       lower_state    = gra::ResolveForwardLegState(lts, gra::ForwardBeamLeg::Lower);
  const bool       forward_noflip = tensor.ForwardNoFlip(gra::MTensorPomeronMode::Continuum, lts);

  std::array<gra::MDirac::Spinor, 2> u_a{};
  std::array<gra::MDirac::Spinor, 2> u_b{};
  std::array<gra::MDirac::Spinor, 2> ubar_1{};
  std::array<gra::MDirac::Spinor, 2> ubar_2{};
  if (!forward_noflip) {
    u_a    = tensor.SpinorStates(pa, "u");
    u_b    = tensor.SpinorStates(pb, "u");
    ubar_1 = tensor.SpinorStates(p1, "ubar");
    ubar_2 = tensor.SpinorStates(p2, "ubar");
  }

  const auto iDP_13 = tensor.iD_P((p1 + p3).M2(), lts.t1);
  const auto iDP_24 = tensor.iD_P((p2 + p4).M2(), lts.t2);
  const auto iDP_14 = tensor.iD_P((p1 + p4).M2(), lts.t1);
  const auto iDP_23 = tensor.iD_P((p2 + p3).M2(), lts.t2);
  const auto iG_1a  = tensor.iG_PForwardHE(upper_state);
  const auto iG_2b  = tensor.iG_PForwardHE(lower_state);

  FTensor::Index<'a', 4> mu1;
  FTensor::Index<'b', 4> nu1;
  FTensor::Index<'c', 4> rho1;
  FTensor::Index<'d', 4> rho2;
  FTensor::Index<'e', 4> rho3;
  FTensor::Index<'f', 4> rho4;
  FTensor::Index<'g', 4> alpha1;
  FTensor::Index<'h', 4> beta1;
  FTensor::Index<'i', 4> alpha2;
  FTensor::Index<'j', 4> beta2;
  FTensor::Index<'k', 4> mu2;
  FTensor::Index<'l', 4> nu2;

  std::vector<std::complex<double>> out;
  for (const auto ha : gra::spin::BinaryHelicityIndices()) {
    for (const auto hb : gra::spin::BinaryHelicityIndices()) {
      for (const auto h1 : gra::spin::BinaryHelicityIndices()) {
        for (const auto h2 : gra::spin::BinaryHelicityIndices()) {
          if (forward_noflip && (ha != h1 || hb != h2)) { continue; }
          const auto proton_1 = forward_noflip ? iG_1a : tensor.iG_PForward(upper_state, ubar_1[h1], u_a[ha], ha, h1);
          const auto proton_2 = forward_noflip ? iG_2b : tensor.iG_PForward(lower_state, ubar_2[h2], u_b[hb], hb, h2);

          if (spin_x2 == 0) {
            const auto          &meson = parameters.FindPseudoscalar(abs_pdg);
            const auto           iG_ta = tensor.iG_Ppsps(pt, -p3, meson.gPPS);
            const auto           iG_tb = tensor.iG_Ppsps(p4, pt, meson.gPPS);
            const auto           iG_ua = tensor.iG_Ppsps(p4, pu, meson.gPPS);
            const auto           iG_ub = tensor.iG_Ppsps(pu, -p3, meson.gPPS);
            std::complex<double> M_t   = proton_1(mu1, nu1) * iDP_13(mu1, nu1, alpha1, beta1) * iG_ta(alpha1, beta1) *
                                       iG_tb(alpha2, beta2) * iDP_24(alpha2, beta2, mu2, nu2) * proton_2(mu2, nu2);
            std::complex<double> M_u = proton_1(mu1, nu1) * iDP_14(mu1, nu1, alpha1, beta1) * iG_ua(alpha1, beta1) *
                                       iG_ub(alpha2, beta2) * iDP_23(alpha2, beta2, mu2, nu2) * proton_2(mu2, nu2);
            M_t *= tensor.iD_MES0(pt, mass) *
                   gra::math::pow2(gra::regge::FormFactor(pt.M2(), gra::math::pow2(mass), meson.ff_offshell));
            M_u *= tensor.iD_MES0(pu, mass) *
                   gra::math::pow2(gra::regge::FormFactor(pu.M2(), gra::math::pow2(mass), meson.ff_offshell));
            out.push_back(-gra::math::zi * (M_t + M_u));
          } else if (spin_x2 == 1) {
            const auto &baryon = parameters.FindBaryon(abs_pdg);
            const auto &ff = parameters.exchange.FindVertex(995, abs_pdg).ff_transfer;
            const auto  iSF_t  = tensor.iD_F(pt, mass);
            const auto  iSF_u  = tensor.iD_F(pu, mass);
            const auto  v_3    = tensor.SpinorStates(p3, "v");
            const auto  ubar_4 = tensor.SpinorStates(p4, "ubar");
            for (const auto h3 : gra::spin::BinaryHelicityIndices()) {
              for (const auto h4 : gra::spin::BinaryHelicityIndices()) {
                const auto           iG_t = tensor.iG_PppbarP(p4, ubar_4[h4], pt, iSF_t, v_3[h3], -p3, baryon.gPBB, ff);
                const auto           iG_u = tensor.iG_PppbarP(p4, ubar_4[h4], pu, iSF_u, v_3[h3], -p3, baryon.gPBB, ff);
                std::complex<double> M_t  = proton_1(mu1, nu1) * iDP_13(mu1, nu1, alpha1, beta1) *
                                           iG_t(alpha2, beta2, alpha1, beta1) * iDP_24(alpha2, beta2, mu2, nu2) *
                                           proton_2(mu2, nu2);
                std::complex<double> M_u = proton_1(mu1, nu1) * iDP_14(mu1, nu1, alpha1, beta1) *
                                           iG_u(alpha1, beta1, alpha2, beta2) * iDP_23(alpha2, beta2, mu2, nu2) *
                                           proton_2(mu2, nu2);
                M_t *= gra::math::pow2(gra::regge::FormFactor(pt.M2(), gra::math::pow2(mass), baryon.ff_offshell));
                M_u *= gra::math::pow2(gra::regge::FormFactor(pu.M2(), gra::math::pow2(mass), baryon.ff_offshell));
                out.push_back(-gra::math::zi * (M_t + M_u));
              }
            }
          } else {
            const int   pdg    = lts.decaytree[0].p.pdg;
            const auto &vector = parameters.FindVector(pdg);
            const auto  iG_tA  = tensor.iG_Pvv(pt, -p3, vector.gPvv[0], vector.gPvv[1], vector.ff_transfer);
            const auto  iG_tB  = tensor.iG_Pvv(p4, pt, vector.gPvv[0], vector.gPvv[1], vector.ff_transfer);
            const auto  iG_uA  = tensor.iG_Pvv(p4, pu, vector.gPvv[0], vector.gPvv[1], vector.ff_transfer);
            const auto  iG_uB  = tensor.iG_Pvv(pu, -p3, vector.gPvv[0], vector.gPvv[1], vector.ff_transfer);
            const auto  iDV_t  = tensor.iD_V(pt, mass, lts.pfinal[0].M2(), pdg);
            const auto  iDV_u  = tensor.iD_V(pu, mass, lts.pfinal[0].M2(), pdg);
            Tensor2<std::complex<double>, 4, 4> M_t;
            Tensor2<std::complex<double>, 4, 4> M_u;
            Tensor2<std::complex<double>, 4, 4> A;
            Tensor2<std::complex<double>, 4, 4> B;
            A(rho1, rho3) = proton_1(mu1, nu1) * iDP_13(mu1, nu1, alpha1, beta1) * iG_tA(rho1, rho3, alpha1, beta1);
            B(rho4, rho2) = proton_2(mu2, nu2) * iDP_24(alpha2, beta2, mu2, nu2) * iG_tB(rho4, rho2, alpha2, beta2);
            M_t(rho3, rho4) =
                A(rho1, rho3) * iDV_t(rho1, rho2) * B(rho4, rho2) *
                gra::math::pow2(gra::regge::FormFactor(pt.M2(), gra::math::pow2(mass), vector.ff_offshell));
            A(rho4, rho1) = proton_1(mu1, nu1) * iDP_14(mu1, nu1, alpha1, beta1) * iG_uA(rho4, rho1, alpha1, beta1);
            B(rho2, rho3) = proton_2(mu2, nu2) * iDP_23(alpha2, beta2, mu2, nu2) * iG_uB(rho2, rho3, alpha2, beta2);
            M_u(rho3, rho4) =
                A(rho4, rho1) * iDV_u(rho1, rho2) * B(rho2, rho3) *
                gra::math::pow2(gra::regge::FormFactor(pu.M2(), gra::math::pow2(mass), vector.ff_offshell));
            Tensor2<std::complex<double>, 4, 4> total;
            for (const auto mu : tensor.LI) {
              for (const auto nu : tensor.LI) { total(mu, nu) = -gra::math::zi * (M_t(mu, nu) + M_u(mu, nu)); }
            }
            const auto helicities = tensor.MassiveSpin1PolSum(total, p3, p4);
            out.insert(out.end(), helicities.begin(), helicities.end());
          }
        }
      }
    }
  }
  return out;
}

// Contract one forward current with the unreduced rank-four Pomeron propagator
FTensor::Tensor2<std::complex<double>, 4, 4> DensePomeronLeg(
    const gra::MTensorPomeron &tensor, const FTensor::Tensor2<std::complex<double>, 4, 4> &current, double s, double t,
    bool upper) {
  FTensor::Tensor2<std::complex<double>, 4, 4> out;
  const auto                                   propagator = tensor.iD_P(s, t);
  for (const auto alpha : tensor.LI) {
    for (const auto beta : tensor.LI) {
      out(alpha, beta) = 0.0;
      for (const auto mu : tensor.LI) {
        for (const auto nu : tensor.LI) {
          out(alpha, beta) += upper ? current(mu, nu) * propagator(mu, nu, alpha, beta)
                                    : propagator(alpha, beta, mu, nu) * current(mu, nu);
        }
      }
    }
  }
  return out;
}

// Build every exact or high-energy Pomeron leg with the unreduced propagator
std::array<FTensor::Tensor2<std::complex<double>, 4, 4>, 4> DenseResonancePomeronLegs(const gra::MTensorPomeron &tensor,
                                                                                      const gra::LORENTZSCALAR  &lts,
                                                                                      gra::ForwardBeamLeg beam_leg,
                                                                                      bool                noflip) {
  const bool                                                  upper    = beam_leg == gra::ForwardBeamLeg::Upper;
  const auto                                                  state    = gra::ResolveForwardLegState(lts, beam_leg);
  const gra::M4Vec                                            incoming = upper ? lts.pbeam1 : lts.pbeam2;
  const gra::M4Vec                                            outgoing = upper ? lts.pfinal[1] : lts.pfinal[2];
  const gra::M4Vec                                            central  = lts.pfinal[0];
  const double s = upper ? (lts.pbeam1 + lts.pfinal[1]) * central
                         : (lts.pbeam2 + lts.pfinal[2]) * central;
  const double                                                t        = upper ? lts.t1 : lts.t2;
  std::array<FTensor::Tensor2<std::complex<double>, 4, 4>, 4> out;

  if (noflip) {
    const auto leg = DensePomeronLeg(tensor, tensor.iG_PForwardHE(state), s, t, upper);
    out[0]         = leg;
    out[3]         = leg;
    return out;
  }

  const auto initial = tensor.SpinorStates(incoming, "u");
  const auto final   = tensor.SpinorStates(outgoing, "ubar");
  for (const auto hi : gra::spin::BinaryHelicityIndices()) {
    for (const auto hf : gra::spin::BinaryHelicityIndices()) {
      const std::size_t index   = gra::spin::BinaryPairHelicityIndex(hi, hf);
      const auto        current = tensor.iG_PForward(state, final[hf], initial[hi], hi, hf);
      out[index]                = DensePomeronLeg(tensor, current, s, t, upper);
    }
  }
  return out;
}

// Evaluate one unreduced rank-two by rank-four by rank-two central scalar
std::complex<double> DenseRank4Central(const gra::MTensorPomeron                                &tensor,
                                       const FTensor::Tensor2<std::complex<double>, 4, 4>       &left,
                                       const FTensor::Tensor4<std::complex<double>, 4, 4, 4, 4> &vertex,
                                       const FTensor::Tensor2<std::complex<double>, 4, 4>       &right) {
  std::complex<double> out = 0.0;
  for (const auto mu : tensor.LI) {
    for (const auto nu : tensor.LI) {
      for (const auto kappa : tensor.LI) {
        for (const auto lambda : tensor.LI) {
          out += left(mu, nu) * vertex(mu, nu, kappa, lambda) * right(kappa, lambda);
        }
      }
    }
  }
  return out;
}

// Evaluate the unreduced PP axial current before spin-frame projection
FTensor::Tensor1<std::complex<double>, 4> DenseAxialCurrent(const gra::MTensorPomeron                          &tensor,
                                                            const FTensor::Tensor2<std::complex<double>, 4, 4> &left,
                                                            const gra::MTensor<std::complex<double>>           &vertex,
                                                            const FTensor::Tensor2<std::complex<double>, 4, 4> &right) {
  FTensor::Tensor1<std::complex<double>, 4> out;
  for (const auto alpha : tensor.LI) {
    out(alpha) = 0.0;
    for (const auto mu : tensor.LI) {
      for (const auto nu : tensor.LI) {
        for (const auto kappa : tensor.LI) {
          for (const auto lambda : tensor.LI) {
            out(alpha) += left(mu, nu) * vertex({mu, nu, kappa, lambda, alpha}) * right(kappa, lambda);
          }
        }
      }
    }
  }
  return out;
}

// Transform one lower-index complex current into the central CM frame
std::array<std::complex<double>, 4> AxialCMCurrentForTest(const gra::LORENTZSCALAR                        &lts,
                                                          const FTensor::Tensor1<std::complex<double>, 4> &current) {
  std::vector<gra::M4Vec> real = {
      gra::M4Vec(-std::real(current(1)), -std::real(current(2)), -std::real(current(3)), std::real(current(0)))};
  std::vector<gra::M4Vec> imag = {
      gra::M4Vec(-std::imag(current(1)), -std::imag(current(2)), -std::imag(current(3)), std::imag(current(0)))};
  gra::kinematics::CMframe(real, lts.pfinal[0]);
  gra::kinematics::CMframe(imag, lts.pfinal[0]);
  return {{{real[0].E(), imag[0].E()},
           {-real[0].Px(), -imag[0].Px()},
           {-real[0].Py(), -imag[0].Py()},
           {-real[0].Pz(), -imag[0].Pz()}}};
}

// Evaluate the full gamma-vector-Pomeron resonance chain without staging
std::complex<double> DenseVectorResonanceChain(const gra::MTensorPomeron                          &tensor,
                                               const FTensor::Tensor1<std::complex<double>, 4>    &photon,
                                               const FTensor::Tensor2<std::complex<double>, 4, 4> &gamma_vector,
                                               const FTensor::Tensor2<std::complex<double>, 4, 4> &vector_transfer,
                                               const FTensor::Tensor4<std::complex<double>, 4, 4, 4, 4> &pomeron_vector,
                                               const FTensor::Tensor2<std::complex<double>, 4, 4>       &resonance,
                                               const FTensor::Tensor1<std::complex<double>, 4>          &decay,
                                               const FTensor::Tensor2<std::complex<double>, 4, 4>       &pomeron) {
  std::complex<double> out = 0.0;
  for (const auto mu : tensor.LI) {
    for (const auto nu : tensor.LI) {
      for (const auto rho : tensor.LI) {
        for (const auto rho2 : tensor.LI) {
          for (const auto alpha : tensor.LI) {
            for (const auto beta : tensor.LI) {
              for (const auto nu2 : tensor.LI) {
                out += photon(mu) * gamma_vector(mu, nu) * vector_transfer(nu, rho) *
                       pomeron_vector(rho2, rho, alpha, beta) * resonance(rho2, nu2) * decay(nu2) *
                       pomeron(alpha, beta);
              }
            }
          }
        }
      }
    }
  }
  return out;
}

}  // namespace

TEST_CASE(
    "Tensor-Pomeron parameters construct safely from one tune across "
    "threads",
    "[gra::MTensorPomeron][params][threading]") {
  const auto       model_tune = gra::MModelTune::Load(modelfile);
  auto             pdg_table  = LoadedPDGTable();
  gra::MModelCache cache(model_tune);
  const auto       first  = gra::GetTensorParam(cache, pdg_table, {});
  const auto       second = gra::GetTensorParam(cache, pdg_table, {});
  REQUIRE(first == second);
  REQUIRE(first->initialized);
  REQUIRE_FALSE(first->VMD.empty());

  constexpr std::size_t                                        nthreads = 8;
  std::vector<std::shared_ptr<const gra::MTensorPomeronParam>> handles(nthreads);
  std::vector<std::thread>                                     workers;
  workers.reserve(nthreads);

  for (std::size_t i = 0; i < nthreads; ++i) {
    workers.emplace_back([i, &handles, &cache, &pdg_table] { handles[i] = gra::GetTensorParam(cache, pdg_table, {}); });
  }
  for (auto &worker : workers) { worker.join(); }

  for (const auto &handle : handles) { REQUIRE(handle == first); }
}

TEST_CASE("Tensor-Pomeron parameters use the immutable PDG mass snapshot", "[gra::MTensorPomeron][params][snapshot]") {
  const auto   model_tune       = gra::MModelTune::Load(modelfile);
  const auto   nominal_pdg      = LoadedPDGTable();
  const auto   nominal          = gra::ReadTensorPomeronParam(*model_tune, nominal_pdg);
  const double nominal_rho_mass = nominal->FindVMD(113).mass;

  auto shifted_pdg = nominal_pdg;
  shifted_pdg.PDG_table.at(113).mass += 0.0123;
  const auto shifted = gra::ReadTensorPomeronParam(*model_tune, shifted_pdg);
  REQUIRE(shifted != nominal);
  CHECK(nominal->FindVMD(113).mass == Approx(nominal_rho_mass));
  CHECK(shifted->FindVMD(113).mass == Approx(nominal_rho_mass + 0.0123));
}

// Check two continuum vertices retain the loaded tune and scale the physical amplitude
TEST_CASE("Tensor continuum amplitudes use the immutable coupling card",
          "[MTensorPomeron][continuum][normalization][snapshot]") {
  const auto tune = WriteModifiedPhotoVMTune("tensor_continuum_snapshot", [](auto &general) {
    general["PARAM_TENSORPOM"]["PARAM_CON"]["TP"]["[211,-211]"] = {{995, 995}};
  });
  const auto model = gra::MModelTune::Load(tune.second);
  const auto evaluate = [](const gra::MModelTunePtr &loaded) {
    auto lts = DirectCentralPairLTSForTest(211, -211);
    gra::MTensorPomeron tensor(lts, loaded,
        gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
    return tensor.ME4(lts, gra::TensorContinuumMode::TensorPomeron);
  };
  const double original = evaluate(model);
  REQUIRE(std::isfinite(original));
  REQUIRE(original > 0.0);
  const std::filesystem::path path = std::filesystem::path(tune.first) / "CON_TP.json";
  auto card = nlohmann::json::parse(gra::aux::GetInputData(path.string()));
  constexpr double scale = 1.7;
  auto &coupling = card["995"]["[211,211]"]["g_tensor"][0];
  coupling = scale * coupling.get<double>();
  std::ofstream output(path);
  REQUIRE(output.good());
  output << card.dump(2);
  output.close();
  CHECK(evaluate(model) == Approx(original).epsilon(1.0e-12));
  CHECK(evaluate(gra::MModelTune::Load(tune.second)) == Approx(std::pow(scale, 4) * original).epsilon(1.0e-12));
}

TEST_CASE("MTensorPomeron rejects beams without implemented forward currents",
          "[gra::MTensorPomeron][initial-state][validation]") {
  const auto soft_model = gra::MModelTune::Load(modelfile);
  const auto p_pbar     = ProtonAntiprotonInitialState();

  for (const bool antiproton_is_upper : {false, true}) {
    CAPTURE(antiproton_is_upper);
    gra::LORENTZSCALAR lts = DirectCentralPairLTSForTest(211, -211);
    lts.beam1              = antiproton_is_upper ? p_pbar[1] : p_pbar[0];
    lts.beam2              = antiproton_is_upper ? p_pbar[0] : p_pbar[1];
    REQUIRE_THROWS(gra::MTensorPomeron(lts, soft_model,
                                       gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic)));
  }
}

// Check that beam support follows the PDG identity, not a duplicated mass
TEST_CASE("MTensorPomeron identifies proton beams only from their PDG IDs",
          "[gra::MTensorPomeron][initial-state][validation]") {
  const auto         soft_model = gra::MModelTune::Load(modelfile);
  gra::LORENTZSCALAR lts        = DirectCentralPairLTSForTest(211, -211);
  REQUIRE(lts.beam1.pdg == gra::PDG::PDG_p);
  REQUIRE(lts.beam2.pdg == gra::PDG::PDG_p);
  lts.beam1.mass += 0.0123;
  lts.beam2.mass += 0.0456;

  REQUIRE_NOTHROW(gra::MTensorPomeron(lts, soft_model,
                                      gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic)));
}

TEST_CASE("MTensorPomeron rejects incomplete tensor-specific steering", "[gra::MTensorPomeron][params]") {
  const auto tune = WriteModifiedPhotoVMTune("tensor_missing_vector_trajectory",
                                             [](auto &j) { j.at("PARAM_TENSORPOM").at("VECTOR").erase("a0"); });
  REQUIRE_THROWS(gra::ReadTensorPomeronParam(*gra::MModelTune::Load(tune.second), LoadedPDGTable()));
}

TEST_CASE("MTensorPomeron requires its independent forward-spin steering", "[gra::MTensorPomeron][params][helicity]") {
  const auto tune = WriteModifiedPhotoVMTune("tensor_missing_forward_noflip",
                                             [](auto &j) { j.at("PARAM_TENSORPOM").erase("FORWARD_NOFLIP"); });
  REQUIRE_THROWS(gra::ReadTensorPomeronParam(*gra::MModelTune::Load(tune.second), LoadedPDGTable()));
}

// Require explicit phase steering independently of the current tune value
TEST_CASE("MTensorPomeron requires explicit decay phase steering", "[gra::MTensorPomeron][params][phase]") {
  const auto missing_tune = WriteModifiedPhotoVMTune(
      "tensor_missing_use_zeta", [](auto &j) { j.at("PARAM_TENSORPOM").erase("use_zeta"); });
  REQUIRE_THROWS(
      gra::ReadTensorPomeronParam(*gra::MModelTune::Load(missing_tune.second), LoadedPDGTable()));
}

TEST_CASE("MTensorPomeron accepts parameter free disabled form factors", "[gra::MTensorPomeron][params][form-factor]") {
  const auto tune = WriteModifiedPhotoVMTune(
      "tensor_disabled_form_factors",
      [](auto &) {},
      [](auto &tensor) {
        tensor.at("995").at("[111,111]")["FF_transfer"] = {{"type", "none"}};
      });
  const auto parameters = gra::ReadTensorPomeronParam(*gra::MModelTune::Load(tune.second), LoadedPDGTable());
  CHECK(parameters->exchange.FindVertex(995, 111).ff_transfer.type == gra::regge::FFType::None);
}

TEST_CASE("MTensorPomeron rejects unsupported vector prescriptions", "[gra::MTensorPomeron][params]") {
  SECTION("unsupported vector decay currents") {
    for (const int daughter : {13, 111, 130, 310}) {
      const auto tune = WriteModifiedPhotoVMTune("tensor_vector_daughter_" + std::to_string(daughter),
          [&](auto &j) { j["PARAM_TENSORPOM"]["VECTOR"]["dPDG"][1] = daughter; });
      CHECK_THROWS_AS(gra::ReadTensorPomeronParam(*gra::MModelTune::Load(tune.second), LoadedPDGTable()),
                      std::invalid_argument);
    }
  }

  SECTION("ambiguous vector and VMD species") {
    for (const std::string block : {"VECTOR", "VMD"}) {
      for (const bool duplicate : {false, true}) {
        CAPTURE(block, duplicate);
        const auto tune = WriteModifiedPhotoVMTune("tensor_vector_species_" + block + std::to_string(duplicate),
            [&](auto &j) {
              auto &pdgs = j["PARAM_TENSORPOM"][block]["PDG"];
              pdgs[1] = duplicate ? pdgs[0].template get<int>() : 211;
            });
        CHECK_THROWS_AS(gra::ReadTensorPomeronParam(*gra::MModelTune::Load(tune.second), LoadedPDGTable()),
                        std::invalid_argument);
      }
    }
  }

  SECTION("nonfinite vector pole data") {
    const auto model = gra::MModelTune::Load(modelfile);
    for (const double value : {std::numeric_limits<double>::infinity(),
                               std::numeric_limits<double>::quiet_NaN()}) {
      for (const bool width : {false, true}) {
        auto pdg = LoadedPDGTable();
        auto &rho = pdg.PDG_table.at(113);
        (width ? rho.width : rho.mass) = value;
        CHECK_THROWS_AS(gra::ReadTensorPomeronParam(*model, pdg), std::invalid_argument);
      }
    }
  }

  SECTION("P-wave pole at or below its decay threshold") {
    const auto model = gra::MModelTune::Load(modelfile);
    for (const double fraction : {0.9, 1.0}) {
      auto pdg = LoadedPDGTable();
      pdg.PDG_table.at(113).mass = fraction * 2.0 * pdg.FindByPDG(211).mass;
      CHECK_THROWS_AS(gra::ReadTensorPomeronParam(*model, pdg), std::invalid_argument);
    }
  }

  SECTION("Regge phase") {
    const auto tune = WriteModifiedPhotoVMTune("tensor_bad_vector_phase", [](auto &j) {
      j.at("PARAM_TENSORPOM").at("VECTOR").at("phase").at(0) = "SIGNATURE";
    });
    REQUIRE_THROWS(gra::ReadTensorPomeronParam(*gra::MModelTune::Load(tune.second), LoadedPDGTable()));
  }

  SECTION("Width model") {
    const auto tune = WriteModifiedPhotoVMTune(
        "tensor_bad_vector_width", [](auto &j) { j.at("PARAM_TENSORPOM").at("VECTOR").at("Wmode").at(1) = "RUNNING"; });
    REQUIRE_THROWS(gra::ReadTensorPomeronParam(*gra::MModelTune::Load(tune.second), LoadedPDGTable()));
  }
}

TEST_CASE("MTensorPomeron validates common trajectory parameters", "[gra::MTensorPomeron][params][validation]") {
  SECTION("legacy exchange field") {
    const auto tune = WriteModifiedPhotoVMTune("tensor_legacy_exchange_kind", [](auto &j) {
      auto &exchange   = j.at("PARAM_TENSORPOM").at("EXCHANGES").at("995");
      exchange["kind"] = exchange.at("type");
      exchange.erase("type");
    });
    REQUIRE_THROWS(gra::ReadTensorPomeronParam(*gra::MModelTune::Load(tune.second), LoadedPDGTable()));
  }

  SECTION("secondary slope") {
    const auto tune = WriteModifiedPhotoVMTune("tensor_bad_secondary_slope", [](auto &j) {
      j.at("PARAM_TENSORPOM").at("EXCHANGES").at("9933").at("ap") = 0.0;
    });
    REQUIRE_THROWS(gra::ReadTensorPomeronParam(*gra::MModelTune::Load(tune.second), LoadedPDGTable()));
  }

  SECTION("Odderon discrete phase sign") {
    const auto tune = WriteModifiedPhotoVMTune("tensor_bad_odderon_signature", [](auto &j) {
      j.at("PARAM_TENSORPOM").at("EXCHANGES").at("9993").at("eta") = 0;
    });
    REQUIRE_THROWS(gra::ReadTensorPomeronParam(*gra::MModelTune::Load(tune.second), LoadedPDGTable()));
  }

  SECTION("rank-one exchange parity") {
    auto pdg                 = LoadedPDGTable();
    pdg.PDG_table.at(9993).P = 1;
    REQUIRE_THROWS(gra::ReadTensorPomeronParam(*gra::MModelTune::Load(modelfile), pdg));
  }
}

// Check strict validation of coherent Tensor continuum exchange selections
TEST_CASE("MTensorPomeron rejects malformed continuum exchange banks",
          "[gra::MTensorPomeron][params][continuum][exchange][validation]") {
  const auto require_invalid = [](const std::string &label, const auto &mutate) {
    const auto tune = WriteModifiedPhotoVMTune(label, mutate);
    REQUIRE_THROWS(gra::ReadTensorPomeronParam(*gra::MModelTune::Load(tune.second), LoadedPDGTable()));
  };

  SECTION("empty exchange bank") {
    require_invalid("tensor_empty_exchange_bank",
                    [](auto &j) { j.at("PARAM_TENSORPOM").at("PARAM_CON").at("TP").at("[321,-321]") = nlohmann::json::array(); });
  }

  SECTION("unknown exchange") {
    require_invalid("tensor_unknown_continuum_exchange", [](auto &j) {
      j.at("PARAM_TENSORPOM").at("PARAM_CON").at("TP").at("[321,-321]") = nlohmann::json::array({nlohmann::json::array({995, 123456})});
    });
  }

  SECTION("unordered exchange duplicate") {
    require_invalid("tensor_duplicate_exchange_pair", [](auto &j) {
      j.at("PARAM_TENSORPOM").at("PARAM_CON").at("TP").at("[321,-321]") =
          nlohmann::json::array({nlohmann::json::array({995, 9915}), nlohmann::json::array({9915, 995})});
    });
  }

  SECTION("order equivalent final-state duplicate") {
    require_invalid("tensor_duplicate_final_pair", [](auto &j) {
      j.at("PARAM_TENSORPOM").at("PARAM_CON").at("TP")["[-321,321]"] = nlohmann::json::array({nlohmann::json::array({995, 995})});
    });
  }

  SECTION("missing hadron vertex") {
    require_invalid("tensor_missing_continuum_vertex", [](auto &j) {
      j.at("PARAM_TENSORPOM").at("PARAM_CON").at("TP").at("[211,-211]") = nlohmann::json::array({nlohmann::json::array({995, 9925})});
    });
  }
}

// Reject unsupported currents and forbidden quantum numbers even for zero couplings
TEST_CASE("Tensor continuum validates every crossed vertex before sampling",
          "[MTensorPomeron][continuum][physics][validation]") {
  for (const auto &[exchange, hadron, count] : std::vector<std::tuple<int, int, std::size_t>>{
           {9993, 321, 1}, {9933, 113, 2}, {9933, 111, 1}, {9925, 211, 1}, {9943, 211, 1}}) {
    for (const double coupling : {0.0, 0.7}) {
      CAPTURE(exchange, hadron, coupling);
      const auto tune = WriteModifiedPhotoVMTune(
          "tensor_forbidden_vertex_" + std::to_string(exchange) + "_" + std::to_string(hadron),
          [](auto &) {}, [&](auto &card) {
            const std::string pair = "[" + std::to_string(hadron) + "," + std::to_string(hadron) + "]";
            card[std::to_string(exchange)][pair] = {{"g_tensor", std::vector<double>(count, coupling)}};
          });
      CHECK_THROWS_AS(gra::ReadTensorPomeronParam(*gra::MModelTune::Load(tune.second), LoadedPDGTable()),
                      std::invalid_argument);
    }
  }
}

// Check the revised Drell-Soding current against its electromagnetic Ward
// identity
TEST_CASE("Tensor pion photoproduction current is gauge invariant", "[gra::MTensorPomeron][photo][gauge][physics]") {
  auto                lts       = DirectCentralPairLTSForTest(gra::PDG::PDG_pip, gra::PDG::PDG_pim);
  const auto          model     = gra::MModelTune::Load(modelfile);
  const auto          parameter = gra::ReadTensorPomeronParam(*model, lts.PDG);
  gra::MTensorPomeron tensor(lts, model, gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Photo),
                             gra::MTensorPomeronMode::Photo);
  const gra::MTensorPhoto photo(tensor, parameter);

  for (const bool photon_from_upper : {true, false}) {
    const auto           current = photo.DSCurrent(lts, photon_from_upper, 0, 0);
    const gra::M4Vec    &q       = photon_from_upper ? lts.q1 : lts.q2;
    std::complex<double> ward    = 0.0;
    double               scale   = 0.0;
    for (const auto &mu : indices(current)) {
      ward += q[mu] * current[mu];
      scale += std::abs(q[mu] * current[mu]);
    }
    REQUIRE(scale > 0.0);
    CHECK(std::abs(ward) < 2.0e-10 * scale);
  }
}

// Check PHOTO uses the shared Tensor forward proton spin setting
TEST_CASE("Tensor pion photoproduction uses shared forward-spin steering", "[gra::MTensorPomeron][photo][helicity]") {
  const auto evaluate = [](const bool noflip) {
    const auto tune  = WriteModifiedPhotoVMTune("tensor_photo_forward_" + std::to_string(noflip), [noflip](auto &j) {
      j.at("PARAM_TENSORPOM").at("FORWARD_NOFLIP") = noflip;
    });
    auto       lts   = DirectCentralPairLTSForTest(gra::PDG::PDG_pip, gra::PDG::PDG_pim);
    const auto model = gra::MModelTune::Load(tune.second);
    const auto parameter = gra::ReadTensorPomeronParam(*model, lts.PDG);
    gra::MTensorPomeron tensor(lts, model, gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Photo),
                               gra::MTensorPomeronMode::Photo);
    const gra::MTensorPhoto photo(tensor, parameter);
    const double                amp2 = photo.Amp2(lts);
    REQUIRE(std::isfinite(amp2));
    return lts.hamp.size();
  };

  REQUIRE(evaluate(true) == 4);
  REQUIRE(evaluate(false) == 16);
}

// Require distinct Minimal, analytic and Tensor Pomeron selector families
TEST_CASE("TP[CON] rejects direct four-stable topology", "[MTensorPomeron][tensor]") {
  ToyHelicityProcess proc;
  REQUIRE_THROWS_AS([&] {
    ConfigureToyProductionProcess(proc, "TP", "CON", "pi+ pi- pi+ pi-");
    proc.InitializeProcessAmplitude();
  }(), std::invalid_argument);
}

TEST_CASE(
    "Generated cascade stable-leaf proposal activates for coherent "
    "cascade channels",
    "[MProcess][continuum][MTensorPomeron][symmetry]") {
  ToyHelicityProcess proc;
  const auto         tree = TensorRhoCascadeLTSForTest().decaytree;

  CHECK_FALSE(proc.StableLeafProposalActiveForTest(tree, false, "MP", "RES"));
  CHECK(proc.StableLeafProposalActiveForTest(tree, true, "MP", "RES"));
  CHECK(proc.StableLeafProposalActiveForTest(tree, true, "XP", "RES"));
  CHECK(proc.StableLeafProposalActiveForTest(tree, true, "GP", "RES"));
  for (const auto &[family, channel] : std::vector<std::pair<std::string, std::string>>{
           {"MP", "CON"}, {"XP", "CON"}, {"GP", "CON"}, {"MP", "RES+CON"}, {"XP", "RES+CON"}, {"GP", "RES+CON"}}) {
    INFO("family=" << family << ", channel=" << channel);
    CHECK_FALSE(proc.StableLeafProposalActiveForTest(tree, false, family, channel));
    CHECK(proc.StableLeafProposalActiveForTest(tree, true, family, channel));
    CHECK_FALSE(proc.StableLeafProposalActiveForTest(tree, false, family, channel, false));
  }
  for (const auto channel : {"RES", "CON"}) {
    INFO(channel);
    CHECK(proc.StableLeafProposalActiveForTest(tree, false, "TP", channel));
    CHECK(proc.StableLeafProposalActiveForTest(tree, true, "TP", channel));
    CHECK_THROWS_AS(proc.StableLeafProposalActiveForTest(tree, false, "TP", channel, false), std::invalid_argument);
  }
  CHECK_THROWS_AS(proc.StableLeafProposalActiveForTest(tree, false, "TP", "RES+CON"), std::invalid_argument);
  CHECK_THROWS_AS(proc.StableLeafProposalActiveForTest(tree, true, "TP", "RES+CON"), std::invalid_argument);
  CHECK_THROWS_AS(proc.StableLeafProposalActiveForTest(tree, false, "TP", "RES+CON", false), std::invalid_argument);
}

// Check explicit TP couplings preserve the BR normalization used by other models
TEST_CASE("Tensor decay overrides preserve independent BR couplings", "[gra::spin][MTensorPomeron][physics]") {
  const std::string family = GENERATE("MP", "XP", "GP", "TP");
  const int pdg = GENERATE(333, 9010221);
  const auto evaluate = [&](bool explicit_tp) {
    const auto tune = WriteModifiedPhotoVMTune("tensor_decay_override_" + std::to_string(explicit_tp), {});
    const std::string path = tune.first + "/DECAYS.json";
    auto card = nlohmann::json::parse(gra::aux::GetInputData(path));
    auto &row = card.at(std::to_string(pdg)).at("[321,-321]");
    row.at("BR") = 0.23;
    row.erase("g_decay_TP");
    if (explicit_tp) { row["g_decay_TP"] = {0.73}; }
    std::ofstream output(path);
    output << card.dump();
    output.close();
    ToyHelicityProcess proc;
    proc.SetTuneForTest(tune.first);
    proc.SetProcessForTest(family, "RES");
    proc.state.lts.PDG = LoadedPDGTable();
    auto parent = proc.state.lts.PDG.FindByPDG(pdg);
    parent.mass = pdg == 333 ? 1.1 : 0.9;
    parent.width = 0.08;
    return proc.ProcessHelicityStructure(parent,
        {proc.state.lts.PDG.FindByPDG(321), proc.state.lts.PDG.FindByPDG(-321)}, false, true, "", false);
  };
  const auto derived = evaluate(false);
  const auto explicit_tp = evaluate(true);
  REQUIRE(std::abs(derived.g_decay) > 0.0);
  REQUIRE(derived.g_decay_TP.size() == 1);
  REQUIRE(derived.g_decay_TP[0] > 0.0);
  REQUIRE(explicit_tp.g_decay_TP.size() == 1);
  CHECK(explicit_tp.g_decay_TP[0] == Approx(0.73));
  RequireComplexNear(explicit_tp.g_decay, derived.g_decay, 2.0e-13);
}

TEST_CASE("MTensorPomeron resonance production uses common fast conventions", "[MTensorPomeron][physics][vertex]") {
  gra::LORENTZSCALAR lts;
  lts.PDG = LoadedPDGTable();
  const auto tune = gra::MModelTune::Load(modelfile);
  gra::MTensorPomeron tensor(lts, tune,
                             gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
  gra::MRandom random;
  const auto res = gra::resonance::Read("RES/f0_980.json", random, gra::ReggeProductionModel::TP);
  const auto &channel = ToyTensorChannel(res);
  const gra::M4Vec    q1(0.21, -0.11, 0.31, 0.09);
  const gra::M4Vec    q2(-0.15, 0.20, -0.25, 0.12);
  const double        mass = 1.3;

  auto norm4 = [&tensor](const auto &vertex) {
    double out = 0.0;
    for (const auto mu : tensor.LI) {
      for (const auto nu : tensor.LI) {
        for (const auto kappa : tensor.LI) {
          for (const auto lambda : tensor.LI) { out += std::norm(vertex(mu, nu, kappa, lambda)); }
        }
      }
    }
    return out;
  };

  SECTION("enum routing preserves the scalar and pseudoscalar structures") {
    const std::vector<double> coupling    = {0.7, -0.2};
    const double              form_factor = gra::regge::TransferFF(q1.M2(), channel.ff_transfer) *
                               gra::regge::TransferFF(q2.M2(), channel.ff_transfer) *
                               gra::regge::MassFF((q1 + q2).M2(), gra::math::pow2(mass), channel.ff_prod);
    const auto scalar        = tensor.iG_PPS_total(q1, q2, mass, gra::TensorResonanceType::Scalar, coupling, channel);
    const auto scalar0       = tensor.iG_PPS_0();
    const auto scalar1       = tensor.iG_PPS_1(q1, q2, coupling[1]);
    const auto pseudoscalar  = tensor.iG_PPS_total(q1, q2, mass, gra::TensorResonanceType::Pseudoscalar, coupling, channel);
    const auto pseudoscalar0 = tensor.iG_PPPS_0(q1, q2, coupling[0]);
    const auto pseudoscalar1 = tensor.iG_PPPS_1(q1, q2, coupling[1]);

    for (const auto mu : tensor.LI) {
      for (const auto nu : tensor.LI) {
        for (const auto kappa : tensor.LI) {
          for (const auto lambda : tensor.LI) {
            RequireComplexNear(
                scalar(mu, nu, kappa, lambda),
                form_factor * (coupling[0] * scalar0(mu, nu, kappa, lambda) + scalar1(mu, nu, kappa, lambda)), 1e-13);
            RequireComplexNear(
                pseudoscalar(mu, nu, kappa, lambda),
                form_factor * (pseudoscalar0(mu, nu, kappa, lambda) + pseudoscalar1(mu, nu, kappa, lambda)), 1e-13);
          }
        }
      }
    }
    REQUIRE_THROWS(tensor.iG_PPS_total(q1, q2, mass, gra::TensorResonanceType::Tensor, coupling, channel));
  }

  SECTION("all five production models share the configured coupling cutoff") {
    const double            inactive = tune->Global().coupling_min;
    const double            active   = 1.01 * inactive;
    gra::RES_TENSOR_CHANNEL axial_inactive;
    axial_inactive.g_tensor = {inactive, -inactive};
    auto axial_active       = axial_inactive;
    axial_active.g_tensor   = {active, 0.0};
    const std::vector<double> tensor_inactive(7, inactive);
    std::vector<double>       tensor_active(7, 0.0);
    tensor_active[0]                   = active;
    const gra::regge::FFParam transfer = {gra::regge::FFType::Power, gra::regge::FFNorm::Zero, {0.5, 1.0}};

    REQUIRE(norm4(tensor.iG_PPS_total(q1, q2, mass, gra::TensorResonanceType::Scalar, {inactive, -inactive}, channel)) ==
            Approx(0.0));
    REQUIRE(norm4(tensor.iG_PPS_total(q1, q2, mass, gra::TensorResonanceType::Pseudoscalar, {inactive, -inactive}, channel)) ==
            Approx(0.0));
    REQUIRE(norm4(tensor.iG_Pvv(q1, q2, inactive, -inactive, transfer)) == Approx(0.0));
    REQUIRE(DynamicTensorNorm2(tensor.iG_PPA_total(q1, q2, mass, axial_inactive)) == Approx(0.0));
    REQUIRE(DynamicTensorNorm2(tensor.iG_PPT_total(q1, q2, mass, tensor_inactive, channel)) == Approx(0.0));

    REQUIRE(norm4(tensor.iG_PPS_total(q1, q2, mass, gra::TensorResonanceType::Scalar, {active, 0.0}, channel)) > 0.0);
    REQUIRE(norm4(tensor.iG_PPS_total(q1, q2, mass, gra::TensorResonanceType::Pseudoscalar, {active, 0.0}, channel)) > 0.0);
    REQUIRE(norm4(tensor.iG_Pvv(q1, q2, active, 0.0, transfer)) > 0.0);
    REQUIRE(DynamicTensorNorm2(tensor.iG_PPA_total(q1, q2, mass, axial_active)) > 0.0);
    REQUIRE(DynamicTensorNorm2(tensor.iG_PPT_total(q1, q2, mass, tensor_active, channel)) > 0.0);
  }
}

TEST_CASE("Reduced PP tensor contractions equal the full six-index vertex",
          "[MTensorPomeron][physics][contraction][exact]") {
  using FTensor::Tensor2;

  std::vector<std::vector<double>> couplings;
  for (std::size_t mode = 0; mode < 7; ++mode) {
    std::vector<double> one(7, 0.0);
    one[mode] = 0.19 + 0.03 * mode;
    couplings.push_back(one);
  }
  couplings.push_back({0.19, -0.22, 0.25, -0.28, 0.31, -0.34, 0.37});

  for (const bool noflip : {true, false}) {
    CAPTURE(noflip);
    const auto tune       = WriteModifiedPhotoVMTune("tensor_exact_ppt_" + std::to_string(noflip), [noflip](auto &j) {
      j.at("PARAM_TENSORPOM").at("FORWARD_NOFLIP") = noflip;
    });
    const auto soft_model = gra::MModelTune::Load(tune.second);

    for (const auto &coupling : couplings) {
      gra::LORENTZSCALAR lts = DirectCentralPairLTSForTest(211, -211);
      gra::PARAM_RES     res;
      res.p       = ToyParticle("f2", 9000225, 4, lts.pfinal[0].M());
      res.p.P     = 1;
      res.p.C     = 1;
      res.p.width = 0.17;
      SetToyTensorChannel(res, coupling);
      res.hel_decay.g_decay_TP = {0.31};
      lts.process.RESONANCES       = {{"f2_exact", res}};

      gra::MTensorPomeron tensor(lts, soft_model,
                                 gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
      const auto          left  = DenseResonancePomeronLegs(tensor, lts, gra::ForwardBeamLeg::Upper, noflip);
      const auto          right = DenseResonancePomeronLegs(tensor, lts, gra::ForwardBeamLeg::Lower, noflip);
      const auto          full  = tensor.iG_PPT_total(lts.q1, lts.q2, res.p.mass, ToyTensorChannel(res).g_tensor, ToyTensorChannel(res));

      const auto decay_vertex =
          tensor.iG_f2psps(lts.decaytree[0].p4, lts.decaytree[1].p4, res.p.mass, res.hel_decay.g_decay_TP[0], res.hel_decay.ff_decay);
      const auto resonance_propagator = tensor.iD_TMES(lts.pfinal[0], res.p.mass, res.p.width, true);
      Tensor2<std::complex<double>, 4, 4> decay;
      for (const auto rho : tensor.LI) {
        for (const auto sigma : tensor.LI) {
          decay(rho, sigma) = 0.0;
          for (const auto alpha : tensor.LI) {
            for (const auto beta : tensor.LI) {
              decay(rho, sigma) += resonance_propagator(rho, sigma, alpha, beta) * decay_vertex(alpha, beta);
            }
          }
        }
      }

      std::vector<std::complex<double>> expected;
      for (const auto ha : gra::spin::BinaryHelicityIndices()) {
        for (const auto hb : gra::spin::BinaryHelicityIndices()) {
          for (const auto h1 : gra::spin::BinaryHelicityIndices()) {
            for (const auto h2 : gra::spin::BinaryHelicityIndices()) {
              if (noflip && (ha != h1 || hb != h2)) { continue; }
              const std::size_t    upper     = gra::spin::BinaryPairHelicityIndex(ha, h1);
              const std::size_t    lower     = gra::spin::BinaryPairHelicityIndex(hb, h2);
              std::complex<double> amplitude = 0.0;
              for (const auto mu : tensor.LI) {
                for (const auto nu : tensor.LI) {
                  for (const auto kappa : tensor.LI) {
                    for (const auto lambda : tensor.LI) {
                      for (const auto rho : tensor.LI) {
                        for (const auto sigma : tensor.LI) {
                          amplitude += left[upper](mu, nu) * full({mu, nu, kappa, lambda, rho, sigma}) *
                                       right[lower](kappa, lambda) * decay(rho, sigma);
                        }
                      }
                    }
                  }
                }
              }
              expected.push_back((-gra::math::zi) * amplitude);
            }
          }
        }
      }

      const double amp2 = tensor.ME3(lts);
      RequireTensorSymmetries(lts, [&](auto &event) { return tensor.ME3(event); });
      CAPTURE(coupling);
      RequireVectorNear(lts.hamp, expected, 3.0e-10);
      REQUIRE(amp2 == Approx(HelicityAverageAmp2ForTest(expected)).epsilon(3.0e-10));
    }
  }
}

// Check scalar and pseudoscalar ME3 production against unreduced rank-four sums
TEST_CASE("Staged spin-zero resonance production equals dense contractions",
          "[MTensorPomeron][physics][contraction][exact]") {
  const auto tune       = WriteModifiedPhotoVMTune("tensor_exact_spin0",
                                                   [](auto &j) { j.at("PARAM_TENSORPOM").at("FORWARD_NOFLIP") = false; });
  const auto soft_model = gra::MModelTune::Load(tune.second);

  for (const bool pseudoscalar : {false, true}) {
    CAPTURE(pseudoscalar);
    gra::LORENTZSCALAR lts =
        DirectCentralPairLTSForTest(pseudoscalar ? 22 : gra::PDG::PDG_pip, pseudoscalar ? 22 : gra::PDG::PDG_pim);
    gra::PARAM_RES res;
    res.p       = ToyParticle(pseudoscalar ? "eta" : "f0", pseudoscalar ? 9000111 : 9001710, 0, lts.pfinal[0].M());
    res.p.P     = pseudoscalar ? -1 : 1;
    res.p.C     = 1;
    res.p.width = 0.14;
    SetToyTensorChannel(res, {0.37, -0.21});
    res.hel_decay.g_decay_TP = {0.29};
    lts.process.RESONANCES       = {{pseudoscalar ? "eta_exact" : "f0_exact", res}};

    gra::MTensorPomeron tensor(lts, soft_model,
                               gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
    const auto          left  = DenseResonancePomeronLegs(tensor, lts, gra::ForwardBeamLeg::Upper, false);
    const auto          right = DenseResonancePomeronLegs(tensor, lts, gra::ForwardBeamLeg::Lower, false);
    const auto          type = pseudoscalar ? gra::TensorResonanceType::Pseudoscalar : gra::TensorResonanceType::Scalar;
    const auto          vertex = tensor.iG_PPS_total(lts.q1, lts.q2, res.p.mass, type, ToyTensorChannel(res).g_tensor, ToyTensorChannel(res));
    const std::complex<double> propagator = tensor.iD_MES(lts.pfinal[0], res.p.mass, res.p.width);

    std::vector<std::complex<double>> decay;
    if (pseudoscalar) {
      const auto decay_vertex =
          tensor.iG_psvv(lts.decaytree[0].p4, lts.decaytree[1].p4, res.p.mass, res.hel_decay.g_decay_TP[0], res.hel_decay.ff_decay);
      decay = tensor.MasslessSpin1PolSum(decay_vertex, lts.decaytree[0].p4, lts.decaytree[1].p4);
    } else {
      decay.push_back(
          tensor.iG_f0ss(lts.decaytree[0].p4, lts.decaytree[1].p4, res.p.mass, res.hel_decay.g_decay_TP[0], res.hel_decay.ff_decay));
    }

    std::vector<std::complex<double>> expected;
    for (const auto ha : gra::spin::BinaryHelicityIndices()) {
      for (const auto hb : gra::spin::BinaryHelicityIndices()) {
        for (const auto h1 : gra::spin::BinaryHelicityIndices()) {
          for (const auto h2 : gra::spin::BinaryHelicityIndices()) {
            const std::size_t          upper = gra::spin::BinaryPairHelicityIndex(ha, h1);
            const std::size_t          lower = gra::spin::BinaryPairHelicityIndex(hb, h2);
            const std::complex<double> production =
                (-gra::math::zi) * propagator * DenseRank4Central(tensor, left[upper], vertex, right[lower]);
            for (const auto &decay_amplitude : decay) { expected.push_back(production * decay_amplitude); }
          }
        }
      }
    }

    const double amp2 = tensor.ME3(lts);
    RequireTensorSymmetries(lts, [&](auto &event) { return tensor.ME3(event); });
    RequireVectorNear(lts.hamp, expected, 3.0e-11);
    REQUIRE(amp2 == Approx(HelicityAverageAmp2ForTest(expected)).epsilon(3.0e-11));
  }
}

// Construct a conserving K+ K- pi0 state for the inclusive axial decay amplitude
gra::LORENTZSCALAR AxialThreeBodyLTS() {
  auto lts = DirectCentralPairLTSForTest(gra::PDG::PDG_Kp, gra::PDG::PDG_Km);
  lts.decaytree.emplace_back().p = lts.PDG.FindByPDG(111);
  std::vector<double> masses;
  for (const auto &branch : lts.decaytree) { masses.push_back(branch.p.mass); }
  gra::MRandom random;
  random.SetSeed(3017);
  std::vector<gra::M4Vec> daughters;
  const auto weight = gra::kinematics::NBodyPhaseSpace(lts.pfinal[0], lts.pfinal[0].M(), masses, daughters, false, random);
  REQUIRE(weight.GetW() > 0.0);
  gra::M4Vec sum;
  for (const auto &i : indices(lts.decaytree)) {
    lts.decaytree[i].p4 = daughters[i];
    sum += daughters[i];
    REQUIRE(daughters[i].M2() == Approx(gra::math::pow2(masses[i])).epsilon(1e-10));
  }
  REQUIRE(gra::math::CheckEMC(lts.pfinal[0] - sum));
  return lts;
}

// Check axial ME3 production against the unreduced rank-five vertex
TEST_CASE("Staged axial resonance production equals dense contractions",
          "[MTensorPomeron][physics][contraction][axial][exact]") {
  const bool use_zeta = GENERATE(0, 1) != 0;
  CAPTURE(use_zeta);
  const auto tune = WriteModifiedPhotoVMTune("tensor_exact_axial_" + std::to_string(use_zeta), [&](auto &j) {
    j.at("PARAM_TENSORPOM").at("FORWARD_NOFLIP") = false;
    j.at("PARAM_TENSORPOM").at("use_zeta") = use_zeta;
  });
  gra::LORENTZSCALAR lts = AxialThreeBodyLTS();

  gra::PARAM_RES res;
  res.p       = ToyParticle("f1", 9002023, 2, lts.pfinal[0].M());
  res.p.P     = 1;
  res.p.C     = 1;
  res.p.width = 0.12;
  SetToyTensorChannel(res, {0.8, -0.25});
  res.hel_decay.g_decay           = {0.6, 0.1};
  res.hel_decay.ff_decay = {gra::regge::FFType::Gaussian, gra::regge::FFNorm::Pole, {0.73}};
  lts.process.RESONANCES = {{"f1_exact", res}};

  gra::MTensorPomeron        tensor(lts, gra::MModelTune::Load(tune.second),
                                    gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
  const auto                 left       = DenseResonancePomeronLegs(tensor, lts, gra::ForwardBeamLeg::Upper, false);
  const auto                 right      = DenseResonancePomeronLegs(tensor, lts, gra::ForwardBeamLeg::Lower, false);
  const auto                 vertex     = tensor.iG_PPA_total(lts.q1, lts.q2, res.p.mass, ToyTensorChannel(res));
  const std::complex<double> propagator = tensor.iD_MES(lts.pfinal[0], res.p.mass, res.p.width);

  gra::PARAM_RES decay_res = res;
  gra::spin::DecayAmp(lts, decay_res, "CM");
  REQUIRE(decay_res.decay_f.size_row() == 3);
  decay_res.decay_f = decay_res.decay_f *
                      gra::regge::MassFF(lts.m2, gra::math::pow2(res.p.mass), res.hel_decay.ff_decay);
  const gra::M4Vec rest(0.0, 0.0, 0.0, lts.pfinal[0].M());
  const auto       eps = tensor.MassiveSpin1States(rest, "conj", true);

  std::vector<std::complex<double>> expected;
  for (const auto ha : gra::spin::BinaryHelicityIndices()) {
    for (const auto hb : gra::spin::BinaryHelicityIndices()) {
      for (const auto h1 : gra::spin::BinaryHelicityIndices()) {
        for (const auto h2 : gra::spin::BinaryHelicityIndices()) {
          const std::size_t upper = gra::spin::BinaryPairHelicityIndex(ha, h1);
          const std::size_t lower = gra::spin::BinaryPairHelicityIndex(hb, h2);
          const auto current = AxialCMCurrentForTest(lts, DenseAxialCurrent(tensor, left[upper], vertex, right[lower]));
          gra::MMatrix<std::complex<double>> production(1, 3, 0.0);
          for (const auto h : indices(eps)) {
            for (const auto alpha : tensor.LI) { production[0][h] += current[alpha] * eps[h](alpha); }
            production[0][h] *= -gra::math::zi;
          }
          const std::complex<double> coupling = use_zeta ? decay_res.hel_decay.g_decay
                                                         : std::complex<double>(std::abs(decay_res.hel_decay.g_decay), 0.0);
          const gra::MMatrix<std::complex<double>> amplitude = (production * decay_res.decay_f) * (propagator * coupling);
          for (std::size_t column = 0; column < amplitude.size_col(); ++column) {
            expected.push_back(amplitude[0][column]);
          }
        }
      }
    }
  }

  const double amp2 = tensor.ME3(lts);
  RequireTensorSymmetries(lts, [&](auto &event) { return tensor.ME3(event); });
  RequireVectorNear(lts.hamp, expected, 3.0e-10);
  REQUIRE(amp2 == Approx(HelicityAverageAmp2ForTest(expected)).epsilon(3.0e-10));
}

// Check vector photoproduction against the unreduced gamma and Pomeron chain
TEST_CASE("Staged vector resonance production equals dense contractions",
          "[MTensorPomeron][physics][contraction][vector][exact]") {
  const bool mixing = GENERATE(false, true);
  const auto tune = WriteModifiedPhotoVMTune("tensor_exact_vector", [](auto &j) {
    j.at("PARAM_TENSORPOM").at("FORWARD_NOFLIP") = false;
    j.at("PARAM_TENSORPOM").at("VECTOR").at("Wmode").at(0) = "P_WAVE_2BODY";
  });
  const auto         soft_model = gra::MModelTune::Load(tune.second);
  gra::LORENTZSCALAR lts        = DirectCentralPairLTSForTest(gra::PDG::PDG_pip, gra::PDG::PDG_pim);

  gra::PARAM_RES res;
  res.p                        = ToyParticle("rho_exact", 113, 2, 0.81);
  res.p.P                      = -1;
  res.p.C                      = -1;
  res.p.width                  = 0.12;
  auto &channel                = SetToyTensorChannel(res, {0.43, -0.18});
  channel.ff_transfer          = TensorVectorTransferForTest(res.p.pdg);
  channel.ff_prod              = gra::regge::ReadFF(
      {{"type", "vector"}, {"norm", "pole"}, {"Lambda2", 3.4}, {"n", 0.7}}, "vector test");
  if (mixing) {
    channel.VMD_MIXING = {{223, {0.43, -0.18}, {}}};
    channel.g_tensor = {0.0, 0.0};
  }
  res.hel_decay.g_decay_TP = {0.71};
  lts.process.RESONANCES       = {{"rho_exact", res}};

  gra::MTensorPomeron tensor(lts, soft_model,
                             gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
  const auto          ff_transfer      = channel.ff_transfer;
  const auto          left             = DenseResonancePomeronLegs(tensor, lts, gra::ForwardBeamLeg::Upper, false);
  const auto          right            = DenseResonancePomeronLegs(tensor, lts, gra::ForwardBeamLeg::Lower, false);
  const auto          gamma_vector     = tensor.iG_yV(0.0, mixing ? 223 : res.p.pdg);
  const auto          vector_1         = mixing ? tensor.iD_VMD(lts.q1, 223) :
      tensor.iD_VMES(lts.q1, res.p.mass, res.p.width, res.p.pdg, true, true);
  const auto          vector_2         = mixing ? tensor.iD_VMD(lts.q2, 223) :
      tensor.iD_VMES(lts.q2, res.p.mass, res.p.width, res.p.pdg, true, true);
  const auto          resonance        = tensor.iD_VMES(lts.pfinal[0], res.p.mass, res.p.width, res.p.pdg, true, true);
  const auto          pomeron_vector_1 = tensor.iG_Pvv(lts.pfinal[0], lts.q1, 0.43, -0.18, ff_transfer);
  const auto          pomeron_vector_2 = tensor.iG_Pvv(lts.pfinal[0], lts.q2, 0.43, -0.18, ff_transfer);
  // Dress both vector virtualities independently of the decay form factor
  // [REFERENCE: Lebiedowicz et al., arXiv:2508.06334v2, Eqs. (2.26)-(2.29)]
  const double incoming_mass = mixing ? lts.PDG.FindByPDG(223).mass : res.p.mass;
  const auto vector_ff = [](double q2, double mass) {
    return std::pow(1.0 + q2 * (q2 - pow2(mass)) / pow2(3.4), -0.7);
  };
  const double ff_out = vector_ff(lts.m2, res.p.mass);
  const double ff_1 = ff_out * vector_ff(lts.q1.M2(), incoming_mass);
  const double ff_2 = ff_out * vector_ff(lts.q2.M2(), incoming_mass);
  const auto          decay =
      tensor.iG_vpsps(lts.decaytree[0].p4, lts.decaytree[1].p4, res.p.mass, res.hel_decay.g_decay_TP[0], res.hel_decay.ff_decay);

  const auto upper_state   = gra::ResolveForwardLegState(lts, gra::ForwardBeamLeg::Upper);
  const auto lower_state   = gra::ResolveForwardLegState(lts, gra::ForwardBeamLeg::Lower);
  const auto upper_initial = tensor.SpinorStates(lts.pbeam1, "u");
  const auto lower_initial = tensor.SpinorStates(lts.pbeam2, "u");
  const auto upper_final   = tensor.SpinorStates(lts.pfinal[1], "ubar");
  const auto lower_final   = tensor.SpinorStates(lts.pfinal[2], "ubar");

  std::vector<std::complex<double>> expected;
  for (const auto ha : gra::spin::BinaryHelicityIndices()) {
    for (const auto hb : gra::spin::BinaryHelicityIndices()) {
      for (const auto h1 : gra::spin::BinaryHelicityIndices()) {
        for (const auto h2 : gra::spin::BinaryHelicityIndices()) {
          const std::size_t upper = gra::spin::BinaryPairHelicityIndex(ha, h1);
          const std::size_t lower = gra::spin::BinaryPairHelicityIndex(hb, h2);
          const auto gamma_1 = tensor.iG_yForwardSources(lts, upper_state, upper_final[h1], upper_initial[ha], ha, h1);
          const auto gamma_2 = tensor.iG_yForwardSources(lts, lower_state, lower_final[h2], lower_initial[hb], hb, h2);
          REQUIRE(gamma_1.size() == 1);
          REQUIRE(gamma_2.size() == 1);

          const std::complex<double> gamma_pomeron = DenseVectorResonanceChain(
              tensor, gamma_1[0], gamma_vector, vector_1, pomeron_vector_1, resonance, decay, right[lower]);
          const std::complex<double> pomeron_gamma = DenseVectorResonanceChain(
              tensor, gamma_2[0], gamma_vector, vector_2, pomeron_vector_2, resonance, decay, left[upper]);
          expected.push_back((-gra::math::zi) * (ff_1 * gamma_pomeron + ff_2 * pomeron_gamma));
        }
      }
    }
  }

  const double amp2 = tensor.ME3(lts);
  RequireTensorSymmetries(lts, [&](auto &event) { return tensor.ME3(event); });
  RequireVectorNear(lts.hamp, expected, 3.0e-10);
  REQUIRE(amp2 == Approx(HelicityAverageAmp2ForTest(expected)).epsilon(3.0e-10));
}

TEST_CASE("Staged vector polarization contractions equal dense sums", "[MTensorPomeron][physics][contraction][exact]") {
  using FTensor::Tensor2;
  using FTensor::Tensor4;

  gra::LORENTZSCALAR lts;
  lts.PDG = LoadedPDGTable();
  gra::MTensorPomeron tensor(lts, gra::MModelTune::Load(modelfile),
                             gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));

  Tensor2<std::complex<double>, 4, 4>       rank2;
  Tensor4<std::complex<double>, 4, 4, 4, 4> rank4;
  for (const auto mu : tensor.LI) {
    for (const auto nu : tensor.LI) {
      rank2(mu, nu) = {0.03 * (1.0 + mu + 2.0 * nu), -0.02 * (1.0 + 3.0 * mu + nu)};
      for (const auto alpha : tensor.LI) {
        for (const auto beta : tensor.LI) {
          rank4(mu, nu, alpha, beta) = {0.01 * (1.0 + mu + 2.0 * nu + alpha + 3.0 * beta),
                                        -0.007 * (1.0 + 2.0 * mu + nu + 3.0 * alpha + beta)};
        }
      }
    }
  }

  const auto check = [&](const auto &eps3, const auto &eps4, const auto &scalar, const auto &matrix) {
    std::size_t state = 0;
    for (const auto &h3 : gra::aux::indices(eps3)) {
      for (const auto &h4 : gra::aux::indices(eps4)) {
        std::complex<double> expected_scalar = 0.0;
        for (const auto mu : tensor.LI) {
          for (const auto nu : tensor.LI) { expected_scalar += eps3[h3](mu) * eps4[h4](nu) * rank2(mu, nu); }
        }
        RequireComplexNear(scalar[state], expected_scalar, 2.0e-13);

        for (const auto alpha : tensor.LI) {
          for (const auto beta : tensor.LI) {
            std::complex<double> expected_matrix = 0.0;
            for (const auto mu : tensor.LI) {
              for (const auto nu : tensor.LI) {
                expected_matrix += eps3[h3](mu) * eps4[h4](nu) * rank4(mu, nu, alpha, beta);
              }
            }
            RequireComplexNear(matrix[state](alpha, beta), expected_matrix, 2.0e-13);
          }
        }
        ++state;
      }
    }
  };

  const gra::M4Vec p3_massless(0.31, -0.12, 0.45, std::sqrt(0.31 * 0.31 + 0.12 * 0.12 + 0.45 * 0.45));
  const gra::M4Vec p4_massless(-0.20, 0.14, -0.35, std::sqrt(0.20 * 0.20 + 0.14 * 0.14 + 0.35 * 0.35));
  const auto       massless3 = tensor.MasslessSpin1States(p3_massless, "conj", true);
  const auto       massless4 = tensor.MasslessSpin1States(p4_massless, "conj", true);
  check(massless3, massless4, tensor.MasslessSpin1PolSum(rank2, p3_massless, p4_massless),
        tensor.MasslessSpin1PolSum(rank4, p3_massless, p4_massless));

  const double     vector_mass = 0.77;
  const gra::M4Vec p3_massive(0.31, -0.12, 0.45,
                              std::sqrt(vector_mass * vector_mass + 0.31 * 0.31 + 0.12 * 0.12 + 0.45 * 0.45));
  const gra::M4Vec p4_massive(-0.20, 0.14, -0.35,
                              std::sqrt(vector_mass * vector_mass + 0.20 * 0.20 + 0.14 * 0.14 + 0.35 * 0.35));
  const auto       massive3 = tensor.MassiveSpin1States(p3_massive, "conj", true);
  const auto       massive4 = tensor.MassiveSpin1States(p4_massive, "conj", true);
  check(massive3, massive4, tensor.MassiveSpin1PolSum(rank2, p3_massive, p4_massive),
        tensor.MassiveSpin1PolSum(rank4, p3_massive, p4_massive));
}

// Recover partial widths from the actual scalar, vector and tensor decay vertices
TEST_CASE("Tensor decay vertices reproduce their partial widths", "[MTensorPomeron][decay][physics][normalization]") {
  gra::LORENTZSCALAR lts;
  lts.PDG = LoadedPDGTable();
  gra::MTensorPomeron tensor(lts, gra::MModelTune::Load(modelfile),
      gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
  const double mass = 1.4;
  const double width = 0.08;
  const double daughter = 0.2;
  const double br = 0.35;
  const gra::M4Vec parent(0.0, 0.0, 0.0, mass);
  for (const double theta : {0.31, 1.27, 2.35}) {
    const auto pair = TwoBodyRestKinematics(mass, daughter, daughter, theta, 0.47);
    for (const int spin : {0, 1, 2}) {
      for (const double symmetry : {1.0, 2.0}) {
        // Two identical spin-zero bosons cannot carry odd orbital angular momentum
        if (spin == 1 && symmetry > 1.5) { continue; }
        const double coupling = gra::MTensorPomeron::GDecay(spin, mass, width, daughter, br, symmetry);
        double summed = 0.0;
        if (spin == 0) {
          summed = std::norm(tensor.iG_f0ss(pair[0], pair[1], mass, coupling, gra::regge::ReadFF({{"type", "none"}}, "test decay")));
        } else if (spin == 1) {
          const auto vertex = tensor.iG_vpsps(pair[0], pair[1], mass, coupling, gra::regge::ReadFF({{"type", "none"}}, "test decay"));
          for (int h = -1; h <= 1; ++h) {
            const auto eps = tensor.EpsMassiveSpin1(parent, h);
            std::complex<double> amplitude = 0.0;
            for (const auto mu : tensor.LI) { amplitude += eps(mu) * vertex(mu); }
            summed += std::norm(amplitude);
          }
        } else {
          const auto vertex = tensor.iG_f2psps(pair[0], pair[1], mass, coupling, gra::regge::ReadFF({{"type", "none"}}, "test decay"));
          for (int h = -2; h <= 2; ++h) {
            const auto eps = tensor.EpsMassiveSpin2(parent, h);
            std::complex<double> amplitude = 0.0;
            for (const auto mu : tensor.LI) {
              for (const auto nu : tensor.LI) { amplitude += eps(mu, nu) * vertex(mu, nu); }
            }
            summed += std::norm(amplitude);
          }
        }
        const double recovered = gra::kinematics::PDW2body(pow2(mass), pow2(daughter), pow2(daughter),
                                                          summed / (2.0 * spin + 1.0), symmetry);
        CHECK(recovered == Approx(width * br).epsilon(2.0e-12));
      }
    }
  }
  const auto photons = TwoBodyRestKinematics(mass, 0.0, 0.0, 0.71, -0.42);
  const double coupling = gra::MTensorPomeron::GDecayPseudoscalarGammaGamma(mass, width, br);
  const auto amplitudes = tensor.MasslessSpin1PolSum(tensor.iG_psvv(photons[0], photons[1], mass, coupling, gra::regge::ReadFF({{"type", "none"}}, "test decay")),
                                                    photons[0], photons[1]);
  CHECK(gra::kinematics::PDW2body(pow2(mass), 0.0, 0.0, gra::SquaredNorm(amplitudes), 2.0) ==
        Approx(width * br).epsilon(2.0e-12));
}

// Check both photon helicity couplings against the published spin-two partial width
// [REFERENCE: Ewerz et al., arXiv:1309.3478, Eqs. (5.25)-(5.28)]
TEST_CASE("Tensor two photon couplings reproduce the helicity partial widths", "[MTensorPomeron][decay][physics][normalization]") {
  gra::LORENTZSCALAR lts;
  lts.PDG = LoadedPDGTable();
  gra::MTensorPomeron tensor(lts, gra::MModelTune::Load(modelfile),
      gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
  const double mass = 1.4;
  const gra::M4Vec parent(0.0, 0.0, 0.0, mass);
  for (const double theta : {0.31, 1.27, 2.35}) {
    const auto pair = TwoBodyRestKinematics(mass, 0.0, 0.0, theta, 0.47);
    for (const auto &coupling : std::array{std::array{0.31, 0.0}, std::array{0.0, 0.57}, std::array{0.31, 0.57}}) {
      const auto vertex = tensor.iG_f2yy(pair[0], pair[1], mass, coupling[0], coupling[1], {});
      double summed = 0.0;
      for (int h = -2; h <= 2; ++h) {
        const auto eps = tensor.EpsMassiveSpin2(parent, h);
        FTensor::Tensor2<std::complex<double>, 4, 4> current;
        for (const auto mu : tensor.LI) {
          for (const auto nu : tensor.LI) {
            current(mu, nu) = 0.0;
            for (const auto kappa : tensor.LI) {
              for (const auto lambda : tensor.LI) { current(mu, nu) += vertex(mu, nu, kappa, lambda) * eps(kappa, lambda); }
            }
          }
        }
        summed += gra::SquaredNorm(tensor.MasslessSpin1PolSum(current, pair[0], pair[1]));
      }
      const double expected = mass / (80.0 * gra::math::PI) *
          (std::pow(mass, 6) * pow2(coupling[0]) / 6.0 + pow2(mass * coupling[1]));
      CHECK(gra::kinematics::PDW2body(pow2(mass), 0.0, 0.0, summed / 5.0, 2.0) == Approx(expected).epsilon(3.0e-12));
    }
  }
}

// Check the massive spin-two projector and its five physical polarization states
TEST_CASE("Tensor resonance propagator has a transverse traceless spin sum", "[MTensorPomeron][physics][projector]") {
  gra::LORENTZSCALAR lts;
  lts.PDG = LoadedPDGTable();
  gra::MTensorPomeron tensor(lts, gra::MModelTune::Load(modelfile),
      gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
  const gra::M4Vec p(0.3, -0.2, 0.4, 1.6);
  const double mass = 1.3;
  const double width = 0.1;
  const auto up = tensor.iD_TMES(p, mass, width, true);
  const auto lo = tensor.iD_TMES(p, mass, width, false);
  const auto factor = gra::math::zi / std::complex<double>(p.M2() - pow2(mass), mass * width);
  std::complex<double> trace = 0.0;
  for (const auto mu : tensor.LI) {
    for (const auto nu : tensor.LI) {
      trace += tensor.g(mu, mu) * tensor.g(nu, nu) * up(mu, nu, mu, nu) / factor;
      for (const auto rho : tensor.LI) {
        std::complex<double> ward = 0.0;
        std::complex<double> pair_trace = 0.0;
        for (const auto sigma : tensor.LI) {
          ward += (p % sigma) * up(sigma, mu, nu, rho);
          pair_trace += tensor.g(sigma, sigma) * up(sigma, sigma, nu, rho);
          std::complex<double> spin_sum = 0.0;
          for (int h = -2; h <= 2; ++h) {
            const auto eps = tensor.EpsMassiveSpin2(p, h);
            spin_sum += eps(mu, nu) * std::conj(eps(rho, sigma));
          }
          RequireComplexNear(up(mu, nu, rho, sigma) / factor, spin_sum, 2.0e-12);
          RequireComplexNear(lo(mu, nu, rho, sigma), tensor.g(mu, mu) * tensor.g(nu, nu) *
              tensor.g(rho, rho) * tensor.g(sigma, sigma) * up(mu, nu, rho, sigma), 2.0e-12);
        }
        CHECK(std::abs(ward) < 2.0e-12 * std::abs(factor));
        CHECK(std::abs(pair_trace) < 2.0e-12 * std::abs(factor));
      }
    }
  }
  RequireComplexNear(trace, 5.0, 2.0e-12);
  CHECK_THROWS_AS(tensor.iD_TMES(gra::M4Vec(0.0, 0.0, 1.0, 1.0), mass, width, true), gra::AmplitudeFailure);
  CHECK_THROWS_AS(tensor.iD_TMES(gra::M4Vec(0.0, 0.0, 0.0, std::numeric_limits<double>::quiet_NaN()),
                                mass, width, true), gra::AmplitudeFailure);
}

// Check decay-coupling normalization with identical particles
TEST_CASE("MTensorPomeron decay couplings include identical-particle factors", "[MTensorPomeron][physics]") {
  const double mass          = 1.2;
  const double width         = 0.08;
  const double daughter_mass = 0.2;
  const double br            = 0.35;
  const double distinct      = gra::MTensorPomeron::GDecay(0, mass, width, daughter_mass, br, 1.0);
  const double identical     = gra::MTensorPomeron::GDecay(0, mass, width, daughter_mass, br, 2.0);
  REQUIRE(identical / distinct == Approx(std::sqrt(2.0)).epsilon(1e-13));

  const double gamma_coupling  = gra::MTensorPomeron::GDecayPseudoscalarGammaGamma(mass, width, br);
  const double recovered_width = pow2(gamma_coupling) * gra::math::pow3(mass) / (256.0 * gra::math::PI);
  REQUIRE(recovered_width == Approx(width * br).epsilon(1e-13));
}

TEST_CASE("MTensorPomeron resonance input rejects a missing decay tree", "[MTensorPomeron][topology]") {
  const auto definition = gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Resonance);
  REQUIRE_FALSE(definition->MatchProcess({}).has_value());
}

TEST_CASE(
    "MTensorPomeron vector vertices and propagators keep their "
    "literature conventions",
    "[MTensorPomeron][physics][vector]") {
  gra::LORENTZSCALAR lts;
  lts.PDG                        = LoadedPDGTable();
  const auto          soft_model = gra::MModelTune::Load(modelfile);
  gra::MTensorPomeron tensor(lts, soft_model,
                             gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
  const auto          parameters = gra::ReadTensorPomeronParam(*soft_model, LoadedPDGTable());

  SECTION("rho and phi Pomeron vertices use their configured form factors") {
    const gra::M4Vec prime(0.20, -0.10, 0.50, 1.20);
    const gra::M4Vec p(-0.10, 0.15, 0.20, 1.00);
    const double     t              = (prime - p).M2();
    const auto       rho_ff         = parameters->FindVector(113).ff_transfer;
    const auto       phi_ff         = parameters->FindVector(333).ff_transfer;
    const auto       rho_vertex     = tensor.iG_Pvv(prime, p, 0.49, 4.27, rho_ff);
    const auto       phi_vertex     = tensor.iG_Pvv(prime, p, 0.49, 4.27, phi_ff);
    const double     expected_ratio = gra::regge::TransferFF(t, phi_ff) / gra::regge::TransferFF(t, rho_ff);

    double norm = 0.0;
    for (const auto mu : tensor.LI) {
      for (const auto nu : tensor.LI) {
        for (const auto kappa : tensor.LI) {
          for (const auto lambda : tensor.LI) {
            norm += std::norm(rho_vertex(mu, nu, kappa, lambda));
            REQUIRE(std::abs(phi_vertex(mu, nu, kappa, lambda) - expected_ratio * rho_vertex(mu, nu, kappa, lambda)) ==
                    Approx(0.0).margin(1e-12));
          }
        }
      }
    }
    REQUIRE(norm > 0.0);
  }

  SECTION("Reggeized vector propagator retains the Feynman i phase") {
    const gra::M4Vec p(0.20, 0.0, 0.0, 0.0);
    const double     mass        = 0.8;
    const auto       prop        = tensor.iD_V(p, mass, 4.0 * mass * mass, 113);
    const double     denominator = p.M2() - mass * mass;

    for (const auto mu : tensor.LI) {
      for (const auto nu : tensor.LI) {
        const std::complex<double> expected = -gra::math::zi * tensor.g(mu, nu) / denominator;
        REQUIRE(std::abs(prop(mu, nu) - expected) == Approx(0.0).margin(1e-12));
      }
    }
  }

  SECTION("rho exchange uses the real Regge power of Eq. 2.22") {
    const gra::M4Vec           p(0.30, 0.0, 0.0, 0.0);
    const auto                &rho       = parameters->FindVector(113);
    const double               s34       = 8.0 * pow2(rho.mass);
    const double               threshold = 4.0 * pow2(rho.mass);
    const double               alpha     = rho.trajectory_intercept + rho.trajectory_slope * p.M2();
    const double               regge     = std::pow(s34 / threshold, alpha - 1.0);
    const std::complex<double> expected  = -gra::math::zi / (p.M2() - pow2(rho.mass)) * regge;
    const auto                 prop      = tensor.iD_V(p, rho.mass, s34, 113);

    REQUIRE(std::abs(prop(0, 0) - expected) == Approx(0.0).margin(1e-12));
    REQUIRE(std::real(prop(0, 0)) == Approx(0.0).margin(1e-14));
  }

  SECTION("phi exchange uses the published alpha_phi trajectory") {
    const gra::M4Vec p(0.30, 0.0, 0.0, 0.0);
    const double     mass      = parameters->FindVector(333).mass;
    const double     s34       = 8.0 * mass * mass;
    const double     threshold = 4.0 * mass * mass;
    const double     phase     = gra::math::PI / 2.0 * std::exp((threshold - s34) / threshold) - gra::math::PI / 2.0;
    const auto      &phi_parameters  = parameters->FindVector(333);
    const double     alpha           = phi_parameters.trajectory_intercept + phi_parameters.trajectory_slope * p.M2();
    const std::complex<double> regge = std::pow(std::exp(gra::math::zi * phase) * s34 / threshold, alpha - 1.0);
    const std::complex<double> expected = -gra::math::zi / (p.M2() - mass * mass) * regge;
    const auto                 prop     = tensor.iD_V(p, mass, s34, 333);
    REQUIRE(std::abs(prop(0, 0) - expected) == Approx(0.0).margin(1e-12));
  }

  SECTION("Conserved-current vector propagator is finite at null momentum") {
    const gra::M4Vec           null_p(0.0, 0.0, 1.0, 1.0);
    const double               mass  = 0.775;
    const double               width = 0.149;
    const auto                 prop  = tensor.iD_VMES(null_p, mass, width, 113, true, true);
    const std::complex<double> delta = 1.0 / (-mass * mass);

    for (const auto mu : tensor.LI) {
      for (const auto nu : tensor.LI) {
        const std::complex<double> expected = -gra::math::zi * tensor.g(mu, nu) * delta;
        REQUIRE(std::isfinite(std::real(prop(mu, nu))));
        REQUIRE(std::isfinite(std::imag(prop(mu, nu))));
        REQUIRE(std::abs(prop(mu, nu) - expected) == Approx(0.0).margin(1e-12));
      }
    }
    REQUIRE_THROWS_AS(tensor.iD_VMES(null_p, mass, width, 113, true, false), gra::AmplitudeFailure);
  }

  SECTION("phi running width follows the P-wave kaon threshold form") {
    const auto      &phi = parameters->FindVector(333);
    const double     s   = pow2(phi.mass + 0.02);
    const gra::M4Vec p(0.0, 0.0, 0.0, std::sqrt(s));
    const double     threshold        = 4.0 * pow2(phi.decay_daughter_mass);
    const double     pole_phase_space = pow2(phi.mass) - threshold;
    const double     imaginary_self_energy =
        phi.width * pow2(phi.mass) / std::sqrt(s) * std::pow((s - threshold) / pole_phase_space, 1.5);
    const std::complex<double> delta = 1.0 / (s - pow2(phi.mass) + gra::math::zi * imaginary_self_energy);
    const auto                 prop  = tensor.iD_VMES(p, phi.mass, phi.width, 333, true, true);

    for (const auto mu : tensor.LI) {
      for (const auto nu : tensor.LI) {
        const std::complex<double> expected = -gra::math::zi * tensor.g(mu, nu) * delta;
        REQUIRE(std::abs(prop(mu, nu) - expected) == Approx(0.0).margin(1e-12));
      }
    }
  }

  SECTION("phi propagator is real below the timelike cut") {
    const auto                &phi = parameters->FindVector(333);
    const gra::M4Vec           p(0.20, 0.0, 0.0, 0.0);
    const std::complex<double> delta = 1.0 / (p.M2() - pow2(phi.mass));
    const auto                 prop  = tensor.iD_VMES(p, phi.mass, phi.width, 333, true, true);

    for (const auto mu : tensor.LI) {
      for (const auto nu : tensor.LI) {
        const std::complex<double> expected = -gra::math::zi * tensor.g(mu, nu) * delta;
        REQUIRE(std::abs(prop(mu, nu) - expected) == Approx(0.0).margin(1e-12));
      }
    }
  }

  SECTION("rho running width follows the published P-wave threshold form") {
    const auto tune = WriteModifiedPhotoVMTune("tensor_rho_pwave", [](auto &j) {
      j.at("PARAM_TENSORPOM").at("VECTOR").at("Wmode").at(0) = "P_WAVE_2BODY";
    });
    auto event = lts;
    event.model_cache.reset();
    gra::MTensorPomeron running(event, gra::MModelTune::Load(tune.second),
        gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
    const auto      &rho = parameters->FindVector(113);
    const double     s   = 0.90 * 0.90;
    const gra::M4Vec p(0.0, 0.0, 0.0, std::sqrt(s));
    const double     threshold        = 4.0 * pow2(rho.decay_daughter_mass);
    const double     pole_phase_space = pow2(rho.mass) - threshold;
    const double     imaginary_self_energy =
        rho.width * pow2(rho.mass) / std::sqrt(s) * std::pow((s - threshold) / pole_phase_space, 1.5);
    const std::complex<double> delta = 1.0 / (s - pow2(rho.mass) + gra::math::zi * imaginary_self_energy);
    const auto                 prop  = running.iD_VMES(p, rho.mass, rho.width, 113, true, true);

    for (const auto mu : tensor.LI) {
      for (const auto nu : tensor.LI) {
        const std::complex<double> expected = -gra::math::zi * tensor.g(mu, nu) * delta;
        REQUIRE(std::abs(prop(mu, nu) - expected) == Approx(0.0).margin(1e-12));
      }
    }
  }

  SECTION("Dropping the longitudinal projector preserves conserved contractions") {
    const gra::M4Vec            p(0.20, -0.10, 0.30, 1.00);
    const auto                  full   = tensor.iD_VMES(p, 0.775, 0.149, 113, true, false);
    const auto                  stable = tensor.iD_VMES(p, 0.775, 0.149, 113, true, true);
    const std::array<double, 4> left   = {0.20, -1.0, 0.0, 0.0};
    const std::array<double, 4> right  = {-0.10, 0.0, -1.0, 0.0};

    std::complex<double> full_contraction   = 0.0;
    std::complex<double> stable_contraction = 0.0;
    for (const auto mu : tensor.LI) {
      for (const auto nu : tensor.LI) {
        full_contraction += left[mu] * full(mu, nu) * right[nu];
        stable_contraction += left[mu] * stable(mu, nu) * right[nu];
      }
    }
    REQUIRE(std::abs(full_contraction - stable_contraction) == Approx(0.0).margin(1e-12));
  }

  SECTION("resonance to vector vertices have no timelike daughter monopole poles") {
    const double     virtuality = parameters->meson_ff_transfer.param.at(0);
    const gra::M4Vec k1(0.0, 0.0, 0.0, std::sqrt(virtuality));
    const gra::M4Vec k2(0.0, 0.0, 0.0, std::sqrt(virtuality));
    const auto       scalar = tensor.iG_f0vv(k1, k2, 1.5, 0.7, -0.2, gra::regge::ReadFF({{"type", "none"}}, "test decay"));
    const auto       spin2  = tensor.iG_f2vv(k1, k2, 1.5, 0.7, -0.2, gra::regge::ReadFF({{"type", "none"}}, "test decay"));
    for (const auto mu : tensor.LI) {
      for (const auto nu : tensor.LI) {
        REQUIRE(std::isfinite(std::real(scalar(mu, nu))));
        REQUIRE(std::isfinite(std::imag(scalar(mu, nu))));
        for (const auto kappa : tensor.LI) {
          for (const auto lambda : tensor.LI) {
            REQUIRE(std::isfinite(std::real(spin2(mu, nu, kappa, lambda))));
            REQUIRE(std::isfinite(std::imag(spin2(mu, nu, kappa, lambda))));
          }
        }
      }
    }
  }

  SECTION("V to pseudoscalar pairs is the bare literature vertex on shell") {
    const double     mass          = 0.775;
    const double     daughter_mass = 0.13957;
    const double     coupling      = 11.51;
    const auto       daughters     = TwoBodyRestKinematics(mass, daughter_mass, daughter_mass, 0.73, -0.41);
    const auto       vertex        = tensor.iG_vpsps(daughters[0], daughters[1], mass, coupling, gra::regge::ReadFF({{"type", "none"}}, "test decay"));
    const gra::M4Vec difference    = daughters[0] - daughters[1];
    for (const auto mu : tensor.LI) {
      const std::complex<double> expected = -0.5 * gra::math::zi * coupling * (difference % mu);
      REQUIRE(std::abs(vertex(mu) - expected) == Approx(0.0).margin(1e-12));
    }
  }
}

TEST_CASE("MTensorPomeron ME4 routes QED and hadronic final states explicitly", "[MTensorPomeron][QED][continuum]") {
  gra::LORENTZSCALAR  muon_lts = DirectCentralPairLTSForTest(-13, 13);
  gra::MTensorPomeron tensor(muon_lts, gra::MModelTune::Load(modelfile),
                             gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
  const double        qed_forward = tensor.ME4(muon_lts, gra::TensorContinuumMode::QED);
  REQUIRE(std::isfinite(qed_forward));
  REQUIRE(qed_forward > 0.0);
  REQUIRE_FALSE(muon_lts.proton_good_walker.has_value());
  REQUIRE_FALSE(gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Continuum)
                    ->MatchProcess(muon_lts.decaytree).has_value());

  std::swap(muon_lts.decaytree[0], muon_lts.decaytree[1]);
  const double qed_reversed = tensor.ME4(muon_lts, gra::TensorContinuumMode::QED);
  REQUIRE(qed_reversed == Approx(qed_forward).epsilon(1e-10));
  REQUIRE_FALSE(muon_lts.proton_good_walker.has_value());

  gra::LORENTZSCALAR  pion_lts = DirectCentralPairLTSForTest(211, -211);
  gra::MTensorPomeron pion_tensor(pion_lts, gra::MModelTune::Load(modelfile),
                                  gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
  const double        tensor_pions = pion_tensor.ME4(pion_lts, gra::TensorContinuumMode::TensorPomeron);
  REQUIRE(std::isfinite(tensor_pions));
  REQUIRE(tensor_pions > 0.0);
  REQUIRE_FALSE(pion_lts.proton_good_walker.has_value());
  REQUIRE_FALSE(gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::QED)
                    ->MatchProcess(pion_lts.decaytree).has_value());
}

TEST_CASE("Tensor forward Pomeron excitation keeps the HE current normalization",
          "[MTensorPomeron][forward][dissociation]") {
  const auto                model_tune        = gra::MModelTune::Load(modelfile);
  const auto                tensor_parameters = gra::ReadTensorPomeronParam(*model_tune, LoadedPDGTable());
  const double              g_pnn             = tensor_parameters->FindBaryon(gra::PDG::PDG_p).gPBB;
  const gra::MDirac::Spinor zero_spinor{};

  for (const auto leg : {gra::ForwardBeamLeg::Upper, gra::ForwardBeamLeg::Lower}) {
    const int index = leg == gra::ForwardBeamLeg::Upper ? 1 : 2;
    CAPTURE(index);
    gra::LORENTZSCALAR lts = DirectCentralPairLTSForTest(211, -211);
    SetToyPhotoForwardExcitation(lts, index, 2.0);
    const auto          state = gra::ResolveForwardLegState(lts, leg);
    gra::MTensorPomeron tensor(lts, model_tune,
                               gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));

    const auto   current = tensor.iG_PForwardHE(state);
    const double profile = model_tune->Soft()->ForwardExcitationFactor(model_tune->Soft()->ForwardExcitationExchange(),
                                                                       state.t, state.mass2);
    const gra::M4Vec momentum_sum = state.outgoing + state.incoming;
    for (const auto mu : tensor.LI) {
      for (const auto nu : tensor.LI) {
        const std::complex<double> expected =
            -gra::math::zi * 3.0 * g_pnn * profile * (momentum_sum % mu) * (momentum_sum % nu);
        RequireComplexNear(current(mu, nu), expected, 2.0e-13);
      }
    }

    const auto reggeon      = tensor.iG_TForwardHE(state, 9915);
    const auto odderon      = tensor.iG_VForwardHE(state, 9993);
    double     reggeon_norm = 0.0;
    double     odderon_norm = 0.0;
    for (const auto mu : tensor.LI) {
      odderon_norm += std::norm(odderon(mu));
      for (const auto nu : tensor.LI) { reggeon_norm += std::norm(reggeon(mu, nu)); }
    }
    REQUIRE(reggeon_norm == Approx(0.0).margin(1.0e-30));
    REQUIRE(odderon_norm == Approx(0.0).margin(1.0e-30));

    const auto diagonal = tensor.iG_PForward(state, zero_spinor, zero_spinor, 0, 0);
    const auto flip     = tensor.iG_PForward(state, zero_spinor, zero_spinor, 0, 1);
    for (const auto mu : tensor.LI) {
      for (const auto nu : tensor.LI) {
        RequireComplexNear(diagonal(mu, nu), current(mu, nu), 2.0e-13);
        RequireComplexNear(flip(mu, nu), 0.0, 1.0e-15);
      }
    }
    REQUIRE_FALSE(lts.proton_good_walker.has_value());
  }
}

TEST_CASE("Tensor exact forward currents use the common collider section",
          "[MTensorPomeron][forward][helicity][phase]") {
  const auto                                 soft_model    = gra::MModelTune::Load(modelfile);
  gra::LORENTZSCALAR                         reference_lts = DirectCentralPairLTSForTest(211, -211);
  gra::MTensorPomeron                        tensor(reference_lts, soft_model,
                                                    gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
  constexpr double                           angle       = 0.47;
  const gra::LORENTZSCALAR                   rotated_lts = RotateToyEventAroundZ(reference_lts, angle);
  const double                               cosine      = std::cos(angle);
  const double                               sine        = std::sin(angle);
  const std::array<std::array<double, 4>, 4> rotation    = {{
         {{1.0, 0.0, 0.0, 0.0}},
         {{0.0, cosine, -sine, 0.0}},
         {{0.0, sine, cosine, 0.0}},
         {{0.0, 0.0, 0.0, 1.0}},
  }};

  for (const auto leg : {gra::ForwardBeamLeg::Upper, gra::ForwardBeamLeg::Lower}) {
    const bool lower           = leg == gra::ForwardBeamLeg::Lower;
    const auto reference_state = gra::ResolveForwardLegState(reference_lts, leg);
    const auto rotated_state   = gra::ResolveForwardLegState(rotated_lts, leg);
    const auto reference_in    = tensor.SpinorStates(reference_state.incoming, "u");
    const auto reference_out   = tensor.SpinorStates(reference_state.outgoing, "ubar");
    const auto rotated_in      = tensor.SpinorStates(rotated_state.incoming, "u");
    const auto rotated_out     = tensor.SpinorStates(rotated_state.outgoing, "ubar");

    for (std::size_t initial = 0; initial < 2; ++initial) {
      for (std::size_t final = 0; final < 2; ++final) {
        const double lambda_in  = initial == 0 ? -0.5 : 0.5;
        const double lambda_out = final == 0 ? -0.5 : 0.5;
        const auto   reference =
            tensor.iG_PForward(reference_state, reference_out[final], reference_in[initial], initial, final);
        const auto rotated = tensor.iG_PForward(rotated_state, rotated_out[final], rotated_in[initial], initial, final);
        const std::complex<double> section =
            gra::spin::ForwardHelicitySectionPhase(lambda_in, lambda_out, angle, lower);
        for (std::size_t mu = 0; mu < 4; ++mu) {
          for (std::size_t nu = 0; nu < 4; ++nu) {
            std::complex<double> expected = 0.0;
            for (std::size_t alpha = 0; alpha < 4; ++alpha) {
              for (std::size_t beta = 0; beta < 4; ++beta) {
                expected += rotation[mu][alpha] * rotation[nu][beta] * reference(alpha, beta);
              }
            }
            expected *= section;
            CAPTURE(lower, initial, final, mu, nu);
            CHECK(std::abs(rotated(mu, nu) - expected) <= 2.0e-10 * std::max(1.0, std::abs(expected)));
          }
        }

        const auto reference_photon = tensor.iG_yForwardSources(reference_lts, reference_state, reference_out[final],
                                                                reference_in[initial], initial, final);
        const auto rotated_photon   = tensor.iG_yForwardSources(rotated_lts, rotated_state, rotated_out[final],
                                                                rotated_in[initial], initial, final);
        REQUIRE(reference_photon.size() == 1);
        REQUIRE(rotated_photon.size() == 1);
        for (std::size_t mu = 0; mu < 4; ++mu) {
          std::complex<double> expected = 0.0;
          for (std::size_t alpha = 0; alpha < 4; ++alpha) {
            expected += rotation[mu][alpha] * reference_photon[0](alpha);
          }
          expected *= section;
          CAPTURE(lower, initial, final, mu);
          CHECK(std::abs(rotated_photon[0](mu) - expected) <= 2.0e-10 * std::max(1.0, std::abs(expected)));
        }
      }
    }
  }
}

TEST_CASE("Reduced Tensor Pomeron currents match the full covariant tensors",
          "[MTensorPomeron][continuum][contraction]") {
  gra::LORENTZSCALAR  lts        = DirectCentralPairLTSForTest(211, -211);
  const auto          soft_model = gra::MModelTune::Load(modelfile);
  gra::MTensorPomeron tensor(lts, soft_model,
                             gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
  const auto          upper_state = gra::ResolveForwardLegState(lts, gra::ForwardBeamLeg::Upper);
  const auto          incoming    = tensor.SpinorStates(upper_state.incoming, "u");
  const auto          outgoing    = tensor.SpinorStates(upper_state.outgoing, "ubar");
  const auto          current     = tensor.iG_PForward(upper_state, outgoing[1], incoming[0], 0, 1);
  const double        subenergy   = (lts.pfinal[1] + lts.decaytree[0].p4).M2();
  const auto          propagator  = tensor.iD_P(subenergy, lts.t1);
  const auto          reduced     = tensor.PomeronPropagatorCurrent(current, subenergy, lts.t1);

  for (const auto alpha : tensor.LI) {
    for (const auto beta : tensor.LI) {
      std::complex<double> expected = 0.0;
      for (const auto mu : tensor.LI) {
        for (const auto nu : tensor.LI) { expected += current(mu, nu) * propagator(mu, nu, alpha, beta); }
      }
      RequireComplexNear(reduced(alpha, beta), expected, 2.0e-13);
    }
  }

  const gra::M4Vec prime = upper_state.outgoing;
  const gra::M4Vec p     = upper_state.incoming;
  const auto       baryon_matrix =
      tensor.TensorBaryonCurrent(reduced, prime, p, 995, gra::PDG::PDG_p);
  const auto           full_vertex       = tensor.iG_Ppp(prime, p, outgoing[0], incoming[1]);
  std::complex<double> expected_bilinear = 0.0;
  for (const auto mu : tensor.LI) {
    for (const auto nu : tensor.LI) { expected_bilinear += reduced(mu, nu) * full_vertex(mu, nu); }
  }
  RequireComplexNear(baryon_matrix.BilinearForm(outgoing[0], incoming[1]), expected_bilinear, 2.0e-12);
}

TEST_CASE("Tensor elastic Pomeron exchange obeys collider pair conventions",
          "[MTensorPomeron][forward][helicity][covariance][reciprocity]") {
  constexpr double momentum = 18.0;
  constexpr double abs_t    = 0.07;
  const double     energy   = std::sqrt(momentum * momentum + gra::PDG::mp * gra::PDG::mp);

  // Construct one ordered elastic proton pair at a transfer azimuth
  const auto elastic_event = [energy, momentum, abs_t](const double azimuth) {
    const double       cosine = 1.0 - abs_t / (2.0 * momentum * momentum);
    const double       sine   = std::sqrt(1.0 - cosine * cosine);
    gra::LORENTZSCALAR lts;
    lts.PDG          = LoadedPDGTable();
    const auto beams = ProtonInitialState();
    lts.beam1        = beams[0];
    lts.beam2        = beams[1];
    lts.pbeam1       = gra::M4Vec(0.0, 0.0, momentum, energy);
    lts.pbeam2       = gra::M4Vec(0.0, 0.0, -momentum, energy);
    lts.pfinal.resize(3);
    lts.pfinal[1] =
        gra::M4Vec(momentum * sine * std::cos(azimuth), momentum * sine * std::sin(azimuth), momentum * cosine, energy);
    lts.pfinal[2] = lts.pbeam1 + lts.pbeam2 - lts.pfinal[1];
    lts.q1        = lts.pbeam1 - lts.pfinal[1];
    lts.q2        = lts.pbeam2 - lts.pfinal[2];
    lts.t1        = lts.q1.M2();
    lts.t2        = lts.q2.M2();
    lts.qt1       = lts.q1.Pt();
    lts.qt2       = lts.q2.Pt();
    lts.s         = (lts.pbeam1 + lts.pbeam2).M2();
    lts.sqrt_s    = std::sqrt(lts.s);
    return lts;
  };

  gra::LORENTZSCALAR  constructor_lts = elastic_event(0.0);
  gra::MTensorPomeron tensor(constructor_lts, gra::MModelTune::Load(modelfile),
                             gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
  constexpr auto      helicity_x2  = gra::spin::BinaryHelicityLabelsX2();
  constexpr auto      spinor_index = gra::spin::BinaryHelicityIndices();

  // Contract both exact beam currents through one Tensor Pomeron propagator
  const auto pair_matrix = [&tensor, &spinor_index](const gra::LORENTZSCALAR &lts) {
    const auto                           upper_state = gra::ResolveForwardLegState(lts, gra::ForwardBeamLeg::Upper);
    const auto                           lower_state = gra::ResolveForwardLegState(lts, gra::ForwardBeamLeg::Lower);
    const auto                           upper_in    = tensor.SpinorStates(upper_state.incoming, "u");
    const auto                           upper_out   = tensor.SpinorStates(upper_state.outgoing, "ubar");
    const auto                           lower_in    = tensor.SpinorStates(lower_state.incoming, "u");
    const auto                           lower_out   = tensor.SpinorStates(lower_state.outgoing, "ubar");
    const auto                           propagator  = tensor.iD_P(lts.s, lts.t1);
    std::array<std::complex<double>, 16> out{};

    for (std::size_t out1 = 0; out1 < 2; ++out1) {
      for (std::size_t out2 = 0; out2 < 2; ++out2) {
        const std::size_t row = gra::spin::BinaryPairHelicityIndex(out1, out2);
        for (std::size_t in1 = 0; in1 < 2; ++in1) {
          for (std::size_t in2 = 0; in2 < 2; ++in2) {
            const std::size_t col      = gra::spin::BinaryPairHelicityIndex(in1, in2);
            const std::size_t raw_in1  = spinor_index[in1];
            const std::size_t raw_in2  = spinor_index[in2];
            const std::size_t raw_out1 = spinor_index[out1];
            const std::size_t raw_out2 = spinor_index[out2];
            const auto        upper =
                tensor.iG_PForward(upper_state, upper_out[raw_out1], upper_in[raw_in1], raw_in1, raw_out1);
            const auto lower =
                tensor.iG_PForward(lower_state, lower_out[raw_out2], lower_in[raw_in2], raw_in2, raw_out2);
            for (const auto mu : tensor.LI) {
              for (const auto nu : tensor.LI) {
                for (const auto kappa : tensor.LI) {
                  for (const auto lambda : tensor.LI) {
                    out[gra::spin::PairHelicityMatrixIndex(row, col)] +=
                        upper(mu, nu) * propagator(mu, nu, kappa, lambda) * lower(kappa, lambda);
                  }
                }
              }
            }
          }
        }
      }
    }
    return out;
  };

  const auto reference = pair_matrix(constructor_lts);
  for (std::size_t row = 0; row < 4; ++row) {
    const std::size_t out1 = row / 2;
    const std::size_t out2 = row % 2;
    for (std::size_t col = 0; col < 4; ++col) {
      const std::size_t          in1   = col / 2;
      const std::size_t          in2   = col % 2;
      const double               sign  = gra::spin::ColliderSpinHalfReciprocitySign(helicity_x2[in1], helicity_x2[in2],
                                                                                    helicity_x2[out1], helicity_x2[out2]);
      const std::size_t          index = gra::spin::PairHelicityMatrixIndex(row, col);
      const std::size_t          transpose = gra::spin::PairHelicityMatrixIndex(col, row);
      const std::complex<double> expected  = sign * reference[transpose];
      CAPTURE(row, col, reference[index], expected);
      CHECK(std::abs(reference[index] - expected) <= 2.0e-10 * std::max(1.0, std::abs(expected)));
    }
  }

  for (const double azimuth : {-1.13, 0.41, 2.07}) {
    const gra::LORENTZSCALAR rotated_lts = elastic_event(azimuth);
    const auto               rotated     = pair_matrix(rotated_lts);
    for (std::size_t row = 0; row < 4; ++row) {
      const std::size_t out1 = row / 2;
      const std::size_t out2 = row % 2;
      for (std::size_t col = 0; col < 4; ++col) {
        const std::size_t in1      = col / 2;
        const std::size_t in2      = col % 2;
        const int         harmonic = gra::spin::ColliderSpinHalfHelicityHarmonic(helicity_x2[in1], helicity_x2[in2],
                                                                                 helicity_x2[out1], helicity_x2[out2]);
        const std::size_t index    = gra::spin::PairHelicityMatrixIndex(row, col);
        const std::complex<double> expected =
            reference[index] * std::exp(gra::math::zi * static_cast<double>(harmonic) * azimuth);
        CAPTURE(row, col, azimuth, harmonic, rotated[index], expected);
        CHECK(std::abs(rotated[index] - expected) <= 2.0e-10 * std::max(1.0, std::abs(expected)));
      }
    }
  }
}

TEST_CASE("Tensor excited photon sources reconstruct the incoherent EPA density",
          "[MTensorPomeron][EPA][dissociation]") {
  const auto                soft_model = gra::MModelTune::Load(modelfile);
  const gra::MDirac::Spinor zero_spinor{};
  for (const auto leg : {gra::ForwardBeamLeg::Upper, gra::ForwardBeamLeg::Lower}) {
    const int index = leg == gra::ForwardBeamLeg::Upper ? 1 : 2;
    CAPTURE(index);
    gra::LORENTZSCALAR lts = DirectCentralPairLTSForTest(-13, 13);
    SetToyPhotoForwardExcitation(lts, index, 2.0);
    const auto          state = gra::ResolveForwardLegState(lts, leg);
    gra::MTensorPomeron tensor(lts, soft_model,
                               gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
    const auto          sources = tensor.iG_yForwardSources(lts, state, zero_spinor, zero_spinor, 0, 0);
    REQUIRE(sources.size() == 2);
    const auto            density = gra::flux::ForwardPhotonFluxTransverse(state);
    std::array<double, 2> norm    = {0.0, 0.0};
    std::complex<double>  overlap = 0.0;
    for (const auto mu : tensor.LI) {
      norm[0] += std::norm(sources[0](mu));
      norm[1] += std::norm(sources[1](mu));
      overlap += std::conj(sources[0](mu)) * sources[1](mu);
    }
    CHECK(norm[0] == Approx(density.parallel / state.xi).epsilon(2.0e-12));
    CHECK(norm[1] == Approx(density.perpendicular / state.xi).epsilon(2.0e-12));
    CHECK(std::abs(overlap) <= 2.0e-12 * std::sqrt(std::max(0.0, norm[0] * norm[1])));

    const double qt = state.transfer.Pt();
    REQUIRE(qt > 0.0);
    const double               qx                  = state.transfer.Px() / qt;
    const double               qy                  = state.transfer.Py() / qt;
    const std::complex<double> parallel_along      = qx * sources[0](1) + qy * sources[0](2);
    const std::complex<double> parallel_cross      = -qy * sources[0](1) + qx * sources[0](2);
    const std::complex<double> perpendicular_along = qx * sources[1](1) + qy * sources[1](2);
    const std::complex<double> perpendicular_cross = -qy * sources[1](1) + qx * sources[1](2);
    CHECK(std::abs(parallel_cross) <= 2.0e-12 * std::max(1.0, std::abs(parallel_along)));
    CHECK(std::abs(perpendicular_along) <= 2.0e-12 * std::max(1.0, std::abs(perpendicular_cross)));

    const auto flip_sources = tensor.iG_yForwardSources(lts, state, zero_spinor, zero_spinor, 0, 1);
    REQUIRE(flip_sources.size() == 2);
    for (const auto &source : flip_sources) {
      for (const auto mu : tensor.LI) { RequireComplexNear(source(mu), 0.0, 1.0e-15); }
    }
    REQUIRE_FALSE(lts.proton_good_walker.has_value());
  }
}

// Build a minimal Tensor forward state for source-bank algebra tests
gra::ForwardLegState TensorForwardStateForTest(const gra::ForwardBeamLeg leg, const bool excited) {
  gra::ForwardLegState state;
  state.leg         = leg;
  state.final_state = excited ? gra::ForwardFinalState::InclusiveExcitation : gra::ForwardFinalState::Elastic;
  return state;
}

TEST_CASE("Tensor forward source banks keep ordered mechanisms incoherent", "[MTensorPomeron][forward][sources]") {
  const auto require_label = [](const gra::TensorForwardSourceBank &bank, const std::size_t component,
                                const gra::TensorForwardSource upper, const gra::TensorForwardSource lower) {
    REQUIRE(bank.Label(component) == std::array<gra::TensorForwardSource, 2>{upper, lower});
  };

  for (const bool upper_excited : {false, true}) {
    for (const bool lower_excited : {false, true}) {
      CAPTURE(upper_excited, lower_excited);
      const auto        upper    = TensorForwardStateForTest(gra::ForwardBeamLeg::Upper, upper_excited);
      const auto        lower    = TensorForwardStateForTest(gra::ForwardBeamLeg::Lower, lower_excited);
      const auto        hadronic = gra::TensorForwardSourceBank::Hadronic(upper, lower);
      const std::size_t expected_hadronic =
          !upper_excited && !lower_excited ? 1 : (upper_excited && lower_excited ? 5 : 3);
      REQUIRE(hadronic.Size() == expected_hadronic);

      const auto &pp            = hadronic.Components(gra::TensorForwardMechanism::PomeronPomeron);
      const auto &gamma_pomeron = hadronic.Components(gra::TensorForwardMechanism::GammaPomeron);
      const auto &pomeron_gamma = hadronic.Components(gra::TensorForwardMechanism::PomeronGamma);
      REQUIRE(pp.size() == 1);
      REQUIRE(gamma_pomeron.size() == (upper_excited ? 2 : 1));
      REQUIRE(pomeron_gamma.size() == (lower_excited ? 2 : 1));
      require_label(hadronic, pp[0],
                    upper_excited ? gra::TensorForwardSource::InclusivePomeron : gra::TensorForwardSource::Elastic,
                    lower_excited ? gra::TensorForwardSource::InclusivePomeron : gra::TensorForwardSource::Elastic);
      for (std::size_t i = 0; i < gamma_pomeron.size(); ++i) {
        require_label(hadronic, gamma_pomeron[i],
                      upper_excited ? (i == 0 ? gra::TensorForwardSource::PhotonParallel
                                              : gra::TensorForwardSource::PhotonPerpendicular)
                                    : gra::TensorForwardSource::Elastic,
                      lower_excited ? gra::TensorForwardSource::InclusivePomeron : gra::TensorForwardSource::Elastic);
      }
      for (std::size_t i = 0; i < pomeron_gamma.size(); ++i) {
        require_label(hadronic, pomeron_gamma[i],
                      upper_excited ? gra::TensorForwardSource::InclusivePomeron : gra::TensorForwardSource::Elastic,
                      lower_excited ? (i == 0 ? gra::TensorForwardSource::PhotonParallel
                                              : gra::TensorForwardSource::PhotonPerpendicular)
                                    : gra::TensorForwardSource::Elastic);
      }

      const auto        photon          = gra::TensorForwardSourceBank::PhotonFusion(upper, lower);
      const std::size_t expected_photon = (upper_excited ? 2 : 1) * (lower_excited ? 2 : 1);
      REQUIRE(photon.Size() == expected_photon);
      const auto &gamma_gamma = photon.Components(gra::TensorForwardMechanism::GammaGamma);
      REQUIRE(gamma_gamma.size() == expected_photon);
      std::size_t component = 0;
      for (std::size_t i = 0; i < (upper_excited ? 2U : 1U); ++i) {
        for (std::size_t j = 0; j < (lower_excited ? 2U : 1U); ++j) {
          require_label(
              photon, gamma_gamma[component++],
              upper_excited
                  ? (i == 0 ? gra::TensorForwardSource::PhotonParallel : gra::TensorForwardSource::PhotonPerpendicular)
                  : gra::TensorForwardSource::Elastic,
              lower_excited
                  ? (j == 0 ? gra::TensorForwardSource::PhotonParallel : gra::TensorForwardSource::PhotonPerpendicular)
                  : gra::TensorForwardSource::Elastic);
        }
      }
      REQUIRE_THROWS_AS(photon.Label(photon.Size()), std::out_of_range);
      REQUIRE_THROWS_AS(photon.Components(static_cast<gra::TensorForwardMechanism>(99)), std::invalid_argument);
    }
  }
}

TEST_CASE("Tensor PP and QED excitation use their own orthogonal source banks",
          "[MTensorPomeron][forward][dissociation][amplitude]") {
  const auto evaluate = [](const bool qed, const bool upper_excited, const bool lower_excited) {
    gra::LORENTZSCALAR lts = DirectCentralPairLTSForTest(qed ? -13 : 211, qed ? 13 : -211);
    if (upper_excited) { SetToyPhotoForwardExcitation(lts, 1, 2.0); }
    if (lower_excited) { SetToyPhotoForwardExcitation(lts, 2, 2.5); }
    lts.process.FORWARD_NOFLIP = true;
    gra::MTensorPomeron tensor(lts, gra::MModelTune::Load(modelfile),
                               gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
    const auto          mode = qed ? gra::TensorContinuumMode::QED : gra::TensorContinuumMode::TensorPomeron;
    const double        amp2 = tensor.ME4(lts, mode);
    const std::size_t   source_count =
        qed ? (upper_excited ? 2U : 1U) * (lower_excited ? 2U : 1U)
              : (!upper_excited && !lower_excited ? 1U : (upper_excited && lower_excited ? 5U : 3U));
    const std::size_t physical_rows = qed ? 16U : 4U;
    REQUIRE(lts.hamp.size() == physical_rows * source_count);
    REQUIRE_FALSE(lts.proton_good_walker.has_value());
    REQUIRE(std::isfinite(amp2));
    REQUIRE(amp2 > 0.0);

    double direct_norm = 0.0;
    for (const auto &amplitude : lts.hamp) {
      REQUIRE(std::isfinite(amplitude.real()));
      REQUIRE(std::isfinite(amplitude.imag()));
      direct_norm += std::norm(amplitude);
    }
    REQUIRE(amp2 == Approx(0.25 * direct_norm).epsilon(2.0e-12));
    if (!qed) {
      for (std::size_t row = 0; row < physical_rows; ++row) {
        std::size_t populated = 0;
        for (std::size_t source = 0; source < source_count; ++source) {
          populated += std::abs(lts.hamp[row * source_count + source]) > 0.0;
        }
        REQUIRE(populated == 1);
      }
    }
  };

  for (const bool qed : {false, true}) {
    evaluate(qed, false, false);
    evaluate(qed, true, false);
    evaluate(qed, false, true);
    evaluate(qed, true, true);
  }
}

TEST_CASE("Tensor gamma-P excitation separates crossed production directions",
          "[MTensorPomeron][forward][photoproduction][dissociation]") {
  const auto evaluate = [](const bool upper_excited, const bool lower_excited) {
    gra::LORENTZSCALAR lts = AsymmetricCentralPairLTSForTest(211, -211, 0.91, -0.38);
    if (upper_excited) { SetToyPhotoForwardExcitation(lts, 1, 2.0); }
    if (lower_excited) { SetToyPhotoForwardExcitation(lts, 2, 2.5); }

    gra::PARAM_RES rho;
    rho.p = lts.PDG.FindByPDG(113);
    SetToyTensorChannel(rho, {0.82, -0.21});
    rho.hel_decay.g_decay_TP = {11.5};
    lts.process.RESONANCES       = {{"rho_test", rho}};

    gra::MTensorPomeron tensor(lts, gra::MModelTune::Load(modelfile),
                               gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
    const double        amp2       = tensor.ME3(lts);
    const std::size_t source_count = !upper_excited && !lower_excited ? 1U : (upper_excited && lower_excited ? 5U : 3U);
    REQUIRE(lts.hamp.size() == 4 * source_count);
    REQUIRE_FALSE(lts.proton_good_walker.has_value());
    REQUIRE(std::isfinite(amp2));
    REQUIRE(amp2 > 0.0);

    double direct_norm = 0.0;
    for (std::size_t row = 0; row < 4; ++row) {
      std::size_t populated = 0;
      for (std::size_t source = 0; source < source_count; ++source) {
        const auto amplitude = lts.hamp[row * source_count + source];
        REQUIRE(std::isfinite(amplitude.real()));
        REQUIRE(std::isfinite(amplitude.imag()));
        direct_norm += std::norm(amplitude);
        populated += std::abs(amplitude) > 0.0;
      }
      const std::size_t expected_populated = upper_excited && lower_excited ? 4U : source_count;
      REQUIRE(populated == expected_populated);
    }
    REQUIRE(amp2 == Approx(0.25 * direct_norm).epsilon(2.0e-12));
  };

  evaluate(false, false);
  evaluate(true, false);
  evaluate(false, true);
  evaluate(true, true);
}

TEST_CASE("MTensorPomeron ME4 separates tensor and QED forward-spin steering",
          "[MTensorPomeron][QED][continuum][helicity]") {
  for (const bool tensor_noflip : {true, false}) {
    CAPTURE(tensor_noflip);
    const auto tune = WriteModifiedPhotoVMTune(
        "tensor_me4_forward_" + std::to_string(tensor_noflip),
        [tensor_noflip](auto &j) { j.at("PARAM_TENSORPOM").at("FORWARD_NOFLIP") = tensor_noflip; });

    // Tensor continuum must ignore the generic PARAM_SPIN-derived runtime flag
    auto tensor_amplitudes = [&tune](bool generic_noflip) {
      gra::LORENTZSCALAR lts     = DirectCentralPairLTSForTest(211, -211);
      lts.process.FORWARD_NOFLIP = generic_noflip;
      gra::MTensorPomeron tensor(lts, gra::MModelTune::Load(tune.second),
                                 gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
      tensor.ME4(lts, gra::TensorContinuumMode::TensorPomeron);
      return lts.hamp;
    };

    const auto        generic_full   = tensor_amplitudes(false);
    const auto        generic_noflip = tensor_amplitudes(true);
    const std::size_t expected_rows  = tensor_noflip ? 4 : 16;
    REQUIRE(generic_full.size() == expected_rows);
    REQUIRE(generic_noflip.size() == expected_rows);
    RequireVectorNear(generic_full, generic_noflip, 1e-13);
  }

  const auto tune = WriteModifiedPhotoVMTune("tensor_me4_qed_independent",
                                             [](auto &j) { j.at("PARAM_TENSORPOM").at("FORWARD_NOFLIP") = true; });
  // QED continuum remains controlled by the generic PARAM_SPIN-derived runtime
  // flag
  auto qed_rows = [&tune](bool generic_noflip) {
    gra::LORENTZSCALAR lts     = DirectCentralPairLTSForTest(-13, 13);
    lts.process.FORWARD_NOFLIP = generic_noflip;
    gra::MTensorPomeron tensor(lts, gra::MModelTune::Load(tune.second),
                               gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
    tensor.ME4(lts, gra::TensorContinuumMode::QED);
    return lts.hamp.size();
  };
  REQUIRE(qed_rows(false) == 64);
  REQUIRE(qed_rows(true) == 16);
}

TEST_CASE("Tensor continuum exposes collider phases in every proton spin row",
          "[MTensorPomeron][continuum][helicity][covariance]") {
  const auto tune       = WriteModifiedPhotoVMTune("tensor_me4_collider_phase",
                                                   [](auto &j) { j.at("PARAM_TENSORPOM").at("FORWARD_NOFLIP") = false; });
  const auto soft_model = gra::MModelTune::Load(tune.second);

  // Evaluate the physical Tensor continuum amplitude in its complete spin basis
  const auto evaluate = [&soft_model](gra::LORENTZSCALAR lts) {
    gra::MTensorPomeron tensor(lts, soft_model,
                               gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
    tensor.ME4(lts, gra::TensorContinuumMode::TensorPomeron);
    return lts.hamp;
  };

  gra::LORENTZSCALAR reference_lts = DirectCentralPairLTSForTest(211, -211);
  const auto         reference     = evaluate(reference_lts);
  REQUIRE(reference.size() == 16);
  REQUIRE(gra::SquaredNorm(reference) > 0.0);
  constexpr auto helicity_x2 = gra::spin::BinaryHelicityLabelsX2();

  for (const double angle : {-1.21, 0.37, 2.04}) {
    const auto rotated = evaluate(RotateToyEventAroundZ(reference_lts, angle));
    REQUIRE(rotated.size() == reference.size());
    for (std::size_t in1 = 0; in1 < 2; ++in1) {
      for (std::size_t in2 = 0; in2 < 2; ++in2) {
        for (std::size_t out1 = 0; out1 < 2; ++out1) {
          for (std::size_t out2 = 0; out2 < 2; ++out2) {
            const std::size_t row      = gra::spin::CanonicalProtonPairSpinLayout::HardRow(in1, in2, out1, out2);
            const int         harmonic = gra::spin::ColliderSpinHalfHelicityHarmonic(helicity_x2[in1], helicity_x2[in2],
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

TEST_CASE("Tensor process metadata follows its mode-specific forward steering",
          "[MTensorPomeron][QED][process][screening][helicity]") {
  const auto require_layout = [](const gra::ScreeningMetadata &metadata, const bool forward_noflip) {
    REQUIRE(metadata.spin_basis == gra::ScreeningSpinBasis::ProtonHelicity);
    REQUIRE(metadata.proton_mode == gra::ProtonScreeningMode::ForwardExcitation);
    CHECK(metadata.forward_noflip == forward_noflip);
    CHECK(metadata.spin_rows == (forward_noflip ? 4 : 16));
    CHECK(metadata.spin_transition_count == (forward_noflip ? 4 : 16));
  };

  for (const auto &[spin_noflip, tensor_noflip] :
       std::array<std::pair<bool, bool>, 2>{std::pair{true, false}, std::pair{false, true}}) {
    CAPTURE(spin_noflip, tensor_noflip);
    const auto             tune = WriteModifiedPhotoVMTune("tensor_process_layout_" + std::to_string(spin_noflip),
                                                           [spin_noflip, tensor_noflip](auto &j) {
                                                 j.at("PARAM_SPIN").at("FORWARD_NOFLIP")      = spin_noflip;
                                                 j.at("PARAM_TENSORPOM").at("FORWARD_NOFLIP") = tensor_noflip;
                                               });
    ModelParamRestoreGuard restore;

    ToyHelicityProcess qed;
    ConfigureToyProductionProcess(qed, "yy", "QED", "mu+ mu-");
    qed.SetExcitation(1);
    qed.SetModelTune(gra::MModelTune::Load(tune.second));
    gra::MODELPARAM = tune.first;
    REQUIRE_NOTHROW(qed.InitializeProcessAmplitude());
    require_layout(qed.state.lts.hamp.metadata, spin_noflip);

    ToyHelicityProcess tensor;
    ConfigureToyProductionProcess(tensor, "TP", "CON", "pi+ pi-");
    tensor.SetExcitation(1);
    tensor.SetModelTune(gra::MModelTune::Load(tune.second));
    gra::MODELPARAM = tune.first;
    REQUIRE_NOTHROW(tensor.InitializeProcessAmplitude());
    require_layout(tensor.state.lts.hamp.metadata, tensor_noflip);
  }
}

TEST_CASE(
    "Tensor elastic amplitudes retain exchange-resolved Good Walker "
    "sources",
    "[MTensorPomeron][screening][GoodWalker][Reggeon][Odderon]") {
  struct PairRun {
    std::vector<std::complex<double>> hamp;
    gra::ProtonGoodWalkerAmplitude           pair;
  };

  const auto evaluate = [](const std::string &label, const std::array<int, 2> final_state,
                           const std::vector<std::array<int, 2>> &exchanges, const std::string &soft_exchange = "") {
    const auto tune        = WriteModifiedPhotoVMTune("tensor_pair_source_" + label, [final_state, exchanges](auto &j) {
      j.at("PARAM_SOFT").at("active_model")        = "double";
      j.at("PARAM_TENSORPOM").at("FORWARD_NOFLIP") = true;
      const std::string key =
          "[" + std::to_string(std::abs(final_state[0])) + ",-" + std::to_string(std::abs(final_state[1])) + "]";
      j.at("PARAM_TENSORPOM").at("PARAM_CON").at("TP").at(key) = exchanges;
    });
    const auto model       = gra::MModelTune::Load(tune.second);
    gra::LORENTZSCALAR lts = DirectCentralPairLTSForTest(final_state[0], final_state[1]);
    lts.hamp.Configure(TensorPairMetadata(true));
    gra::MTensorPomeron tensor(lts, model, gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
    const double        amp2 = tensor.ME4(lts, gra::TensorContinuumMode::TensorPomeron);
    REQUIRE(std::isfinite(amp2));
    REQUIRE(amp2 > 0.0);
    RequireTensorBornProjection(lts, 2);
    if (!soft_exchange.empty()) { RequireTensorResiduePair(lts, model->Soft(), soft_exchange); }
    return PairRun{lts.hamp, *lts.proton_good_walker};
  };

  const PairRun pp = evaluate("pp", {211, -211}, {{{995, 995}}}, "P");
  evaluate("rr", {211, -211}, {{{9915, 9915}}}, "R_f2");
  evaluate("oo", {-2212, 2212}, {{{9993, 9993}}}, "O");
  const PairRun rp    = evaluate("rp", {211, -211}, {{{9915, 995}}});
  const PairRun pp_rp = evaluate("pp_rp", {211, -211}, {{{995, 995}}, {{9915, 995}}});
  REQUIRE(pp.hamp.size() == rp.hamp.size());
  REQUIRE(pp.hamp.size() == pp_rp.hamp.size());
  const auto &pp_source  = pp.pair.components.front().source;
  const auto &rp_source  = rp.pair.components.front().source;
  const auto &sum_source = pp_rp.pair.components.front().source;
  REQUIRE(pp_source.size_row() == rp_source.size_row());
  REQUIRE(pp_source.size_col() == rp_source.size_col());
  REQUIRE(sum_source.size_row() == pp_source.size_row());
  REQUIRE(sum_source.size_col() == pp_source.size_col());
  for (const auto &row : indices(pp.hamp)) {
    RequireComplexNear(pp_rp.hamp[row], pp.hamp[row] + rp.hamp[row], 2.0e-10);
    for (std::size_t col = 0; col < pp_source.size_col(); ++col) {
      RequireComplexNear(sum_source(row, col), pp_source(row, col) + rp_source(row, col), 2.0e-10);
    }
  }

  const PairRun po                = evaluate("po", {-2212, 2212}, {{{995, 9993}}});
  const auto    normalized_source = [](const PairRun &run) {
    const auto &source = run.pair.components.front().source;
    std::size_t row    = 0;
    while (row < run.hamp.size() && std::abs(run.hamp[row]) <= 1.0e-13) { ++row; }
    REQUIRE(row < run.hamp.size());
    std::vector<std::complex<double>> out(source.size_col(), 0.0);
    for (const auto &col : indices(out)) { out[col] = source(row, col) / run.hamp[row]; }
    return out;
  };
  const auto pp_profile = normalized_source(pp);
  const auto rp_profile = normalized_source(rp);
  const auto po_profile = normalized_source(po);
  REQUIRE(pp_profile.size() == 4);
  REQUIRE(rp_profile.size() == pp_profile.size());
  REQUIRE(po_profile.size() == pp_profile.size());
  double pp_rp_distance = 0.0;
  double pp_po_distance = 0.0;
  for (const auto &i : indices(pp_profile)) {
    pp_rp_distance += std::norm(pp_profile[i] - rp_profile[i]);
    pp_po_distance += std::norm(pp_profile[i] - po_profile[i]);
  }
  REQUIRE(pp_rp_distance > 1.0e-8);
  REQUIRE(pp_po_distance > 1.0e-8);
}

TEST_CASE("Tensor exchange aliases map to compatible SOFT residues",
          "[MTensorPomeron][screening][exchange][validation]") {
  const auto  model      = gra::MModelTune::Load(modelfile);
  const auto  parameters = gra::ReadTensorPomeronParam(*model, LoadedPDGTable());
  const auto &exchange   = parameters->exchange;
  CHECK(exchange.SoftId(995, *model->Soft()) == model->Soft()->ExchangeId("P"));
  CHECK(exchange.SoftId(9915, *model->Soft()) == model->Soft()->ExchangeId("R_f2"));
  CHECK(exchange.SoftId(9925, *model->Soft()) == model->Soft()->ExchangeId("R_f2"));
  CHECK(exchange.SoftId(9993, *model->Soft()) == model->Soft()->ExchangeId("O"));
  CHECK(exchange.SoftId(9933, *model->Soft()) == model->Soft()->ExchangeId("R_rho"));
  CHECK(exchange.SoftId(9943, *model->Soft()) == model->Soft()->ExchangeId("R_rho"));
}

// Check the selected exchange content and common off-shell convention
TEST_CASE("Tensor continuum selections resolve the selected DL exchanges",
          "[MTensorPomeron][continuum][exchange][Reggeon][physics]") {
  const auto  model      = gra::MModelTune::Load(modelfile);
  const auto  parameters = gra::ReadTensorPomeronParam(*model, LoadedPDGTable());
  const auto &exchange   = parameters->exchange;

  const auto                            card        = nlohmann::json::parse(gra::aux::GetInputData(modelfile));
  const auto                           &selected    = card.at("PARAM_TENSORPOM").at("PARAM_CON").at("TP");
  const auto                            kaon_baryon = selected.at("[321,-321]").get<std::vector<std::array<int, 2>>>();
  const auto                            pion        = selected.at("[211,-211]").get<std::vector<std::array<int, 2>>>();
  const auto                            proton = selected.at("[2212,-2212]").get<std::vector<std::array<int, 2>>>();
  const std::vector<std::array<int, 2>> equal_order = {{995, 995}};
  const std::vector<std::array<int, 2>> mixed_order = {{995, 9925}, {9925, 995}};
  CHECK(exchange.FindContinuumPairs(321, -321) == kaon_baryon);
  CHECK(exchange.FindContinuumPairs(2212, -2212) == proton);
  CHECK(exchange.FindContinuumPairs(211, -211) == pion);
  for (const auto &pairs : {kaon_baryon, pion, proton}) {
    REQUIRE_FALSE(pairs.empty());
    for (const auto &pair : pairs) {
      CHECK_NOTHROW(exchange.SoftId(pair[0], *model->Soft()));
      CHECK_NOTHROW(exchange.SoftId(pair[1], *model->Soft()));
    }
  }

  CHECK(gra::MTensorExchangeModel::OrderedPairs({995, 995}) == equal_order);
  CHECK(gra::MTensorExchangeModel::OrderedPairs({995, 9925}) == mixed_order);

  const auto &phi_odderon = exchange.FindVertex(995, 333, 9993);
  REQUIRE(phi_odderon.g_tensor.size() == 2);
  CHECK(std::fpclassify(phi_odderon.g_tensor[0]) == FP_ZERO);
  CHECK(phi_odderon.active_g_tensor == std::vector<std::size_t>{1});
  const auto phi_transfers = exchange.FindActiveTransfers(995, 995, 333, 333);
  CHECK(std::find(phi_transfers.cbegin(), phi_transfers.cend(), 333) != phi_transfers.cend());
  CHECK(std::find(phi_transfers.cbegin(), phi_transfers.cend(), 9993) != phi_transfers.cend());

  const auto expected = gra::regge::ReadFF(
      model->Continuum("TP").at("995").at("[321,321]").at("FF_offshell"),
      "test Tensor continuum kaon off-shell form factor");
  for (const int exchange_pdg : {995, 9915, 9925, 9933, 9943}) {
    const auto &vertex = exchange.FindVertex(exchange_pdg, 321);
    CHECK(vertex.ff_offshell.type == expected.type);
    REQUIRE(vertex.ff_offshell.param.size() == expected.param.size());
    for (const auto &i : indices(expected.param)) {
      CHECK(vertex.ff_offshell.param[i] == Approx(expected.param[i]));
    }
  }
}

// Check charge conjugation signs of rank-one and rank-two meson vertices
TEST_CASE("Tensor continuum vertices apply the exchange C parity sign",
          "[MTensorPomeron][continuum][exchange][C-parity][physics]") {
  gra::LORENTZSCALAR lts;
  lts.PDG = LoadedPDGTable();
  gra::MTensorPomeron tensor(lts, gra::MModelTune::Load(modelfile),
                             gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
  const gra::M4Vec    prime(0.17, -0.11, 0.31, 0.83);
  const gra::M4Vec    momentum(-0.09, 0.13, -0.22, 0.67);

  const auto vector_particle     = tensor.iG_Vpsps(prime, momentum, 9933, 321);
  const auto vector_antiparticle = tensor.iG_Vpsps(prime, momentum, 9933, -321);
  for (const auto &mu : tensor.LI) { RequireComplexNear(vector_antiparticle(mu), -vector_particle(mu), 1.0e-13); }

  const auto tensor_particle     = tensor.iG_Tpsps(prime, momentum, 9915, 321);
  const auto tensor_antiparticle = tensor.iG_Tpsps(prime, momentum, 9915, -321);
  for (const auto &mu : tensor.LI) {
    for (const auto &nu : tensor.LI) {
      RequireComplexNear(tensor_antiparticle(mu, nu), tensor_particle(mu, nu), 1.0e-13);
    }
  }
}

TEST_CASE("Tensor forward exchanges require compatible SOFT sources",
          "[MTensorPomeron][screening][exchange][validation]") {
  const auto         tune = WriteModifiedPhotoVMTune("tensor_disabled_soft_odderon", [](auto &j) {
    auto             &soft                                     = j.at("PARAM_SOFT");
    const std::string model                                    = soft.at("active_model");
    soft.at("MODEL").at(model).at("EXCHANGE").at("O").at("on") = false;
    j.at("PARAM_TENSORPOM").at("PARAM_CON").at("TP").at("[211,-211]")                = {{{995, 9993}}};
  });
  ToyHelicityProcess process;
  ConfigureToyProductionProcess(process, "TP", "CON", "pi+ pi-");
  process.SetModelTune(gra::MModelTune::Load(tune.second));
  REQUIRE_THROWS_AS(process.InitializeProcessAmplitude(), std::invalid_argument);
}

TEST_CASE("Tensor resonance production uses the global coupling threshold", "[MTensorPomeron][production][numerics]") {
  for (const double scale : {0.5, 2.0}) {
    CAPTURE(scale);
    ToyHelicityProcess process;
    ConfigureToyProductionProcess(process, "TP", "RES", "pi+ pi-");
    auto         resonance = gra::resonance::Read("RES/f0_500.json", process.state.random, gra::ReggeProductionModel::TP);
    const double coupling  = scale * process.GetModelTune()->Global().coupling_min;
    SetToyTensorChannel(resonance, {coupling, 0.0});
    process.SetResonances({{"f0_500", resonance}});
    if (scale < 1.0) {
      REQUIRE_THROWS_AS(process.InitializeProcessAmplitude(), std::invalid_argument);
    } else {
      REQUIRE_NOTHROW(process.InitializeProcessAmplitude());
    }
  }
}

TEST_CASE(
    "Tensor continuum screening equals an independent N2 pair "
    "contraction",
    "[MTensorPomeron][screening][GoodWalker][multichannel][physics]") {
  ModelParamRestoreGuard restore;
  const double           screening_s = DirectCentralPairLTSForTest(211, -211).s;
  const auto             tune        = WriteModifiedPhotoVMTune("tensor_continuum_screening_n2",
                                                                [](auto &j) { j.at("PARAM_SOFT").at("active_model") = "double"; });
  gra::MODELPARAM                    = tune.first;
  std::filesystem::create_directories(gra::aux::GetBasePath(2) + "/eikonal");
  gra::MEikonal eikonal(gra::MModelTune::Load(tune.second));
  eikonal.S3Constructor(screening_s, ProtonInitialState(), false, 4, 4);
  REQUIRE(eikonal.GetChannelCount() == 2);

  ToyTensorScreeningProcess process(eikonal.ModelTuneHandle());
  const double              born_amp2 = process.BornAmp2();
  REQUIRE(std::isfinite(born_amp2));
  REQUIRE(born_amp2 > 0.0);

  process.eikonal                                  = eikonal;
  process.eikonal.Numerics.LOOP.radial_integrator  = "GL";
  process.eikonal.Numerics.LOOP.azimuth_integrator = "Trap";
  process.eikonal.Numerics.LOOP.radial_map         = gra::math::RadialMap::Linear;
  process.eikonal.Numerics.LOOP.r_min              = 0.0;
  process.eikonal.Numerics.LOOP.r_max              = 0.04;
  process.eikonal.Numerics.LOOP.radial_intervals   = 2;
  process.eikonal.Numerics.LOOP.azimuth_nodes      = 4;
  process.eikonal.InitLoopWeightMatrix();

  const double screened_amp2 = process.ScreenedAmp2();
  const auto   reference     = TensorScreeningReference(process.PairTrace(), process.eikonal);
  RequireVectorNear(process.state.lts.hamp, reference, 3.0e-10);
  const double reference_amp2 = 0.25 * gra::SquaredNorm(reference);
  REQUIRE(screened_amp2 == Approx(reference_amp2).epsilon(3.0e-10).margin(1.0e-12));
  REQUIRE(std::isfinite(screened_amp2));
  REQUIRE(std::abs(screened_amp2 - born_amp2) > 1.0e-12 * std::max(1.0, born_amp2));
}

TEST_CASE("Tensor and full QED processes initialize forward excitation modes",
          "[MTensorPomeron][QED][process][excitation]") {
  ModelParamRestoreGuard restore;
  for (const std::string channel : {"RES", "CON", "RES+CON"}) {
    CAPTURE(channel);
    ToyHelicityProcess tensor;
    ConfigureToyProductionProcess(tensor, "TP", channel, "pi+ pi-");
    tensor.SetExcitation(0);
    if (channel != "CON") {
      tensor.SetResonances({{"f0_980", gra::resonance::Read("RES/f0_980.json", tensor.state.random, gra::ReggeProductionModel::TP)}});
    }
    REQUIRE_NOTHROW(tensor.InitializeProcessAmplitude());
    CHECK(tensor.state.lts.hamp.metadata.amplitude_type == gra::ScreeningAmplitudeType::GoodWalker);
  }

  for (const int excitation : {1, 2}) {
    CAPTURE(excitation);
    for (const std::string channel : {"RES", "CON", "RES+CON"}) {
      CAPTURE(channel);
      ToyHelicityProcess tensor;
      ConfigureToyProductionProcess(tensor, "TP", channel, "pi+ pi-");
      tensor.SetExcitation(excitation);
      if (channel != "CON") {
        tensor.SetResonances({{"f0_980", gra::resonance::Read("RES/f0_980.json", tensor.state.random, gra::ReggeProductionModel::TP)}});
      }
      REQUIRE_NOTHROW(tensor.InitializeProcessAmplitude());
      CHECK(tensor.state.lts.hamp.metadata.proton_mode == gra::ProtonScreeningMode::ForwardExcitation);
      CHECK(tensor.state.lts.hamp.metadata.amplitude_type == gra::ScreeningAmplitudeType::Physical);
      CHECK_FALSE(tensor.state.lts.proton_good_walker.has_value());
    }

    ToyHelicityProcess qed;
    ConfigureToyProductionProcess(qed, "yy", "QED", "mu+ mu-");
    qed.SetExcitation(excitation);
    REQUIRE_NOTHROW(qed.InitializeProcessAmplitude());
    CHECK(qed.state.lts.hamp.metadata.proton_mode == gra::ProtonScreeningMode::ForwardExcitation);
    CHECK(qed.state.lts.hamp.metadata.amplitude_type == gra::ScreeningAmplitudeType::Physical);
    CHECK_FALSE(qed.state.lts.proton_good_walker.has_value());
  }
}

TEST_CASE(
    "Tensor multichannel screening rejects unresolved forward "
    "excitation",
    "[MTensorPomeron][process][screening][excitation][validation]") {
  ModelParamRestoreGuard restore;
  const auto             tune = WriteModifiedPhotoVMTune("tensor_excitation_screening_n2",
                                                         [](auto &j) { j.at("PARAM_SOFT").at("active_model") = "double"; });
  gra::MODELPARAM             = tune.first;
  const auto    model         = gra::MModelTune::Load(tune.second);
  gra::MEikonal eikonal(model);
  eikonal.S3Constructor(25.0, ProtonInitialState(), false, 4, 4);
  REQUIRE(eikonal.GetChannelCount() == 2);

  const auto configure = [&model](ToyHelicityProcess &process) {
    ConfigureToyProductionProcess(process, "TP", "CON", "pi+ pi-");
    process.SetModelTune(model);
    process.SetInitialState({"p+", "p+"}, {2.5, 2.5});
    process.SetExcitation(1);
    process.SetScreening(true);
  };

  SECTION("external eikonal") {
    ToyHelicityProcess process;
    configure(process);
    REQUIRE_THROWS_AS(process.PrepareRun(eikonal), std::invalid_argument);
  }
  SECTION("internal eikonal is rejected before construction") {
    ToyHelicityProcess process;
    configure(process);
    REQUIRE_FALSE(process.GetEikonal().IsInitialized());
    REQUIRE_THROWS_AS(process.PrepareRun(), std::invalid_argument);
    REQUIRE_FALSE(process.GetEikonal().IsInitialized());
  }
}

TEST_CASE("Tensor RES+CON adds coherent Regge pair sources", "[MTensorPomeron][process][screening][GoodWalker]") {
  gra::LORENTZSCALAR base = DirectCentralPairLTSForTest(gra::PDG::PDG_pip, gra::PDG::PDG_pim);
  gra::PARAM_RES     resonance;
  resonance.p                        = ToyParticle("f0_pair_sum", 9000021, 0, base.pfinal[0].M());
  resonance.p.P                      = 1;
  resonance.p.C                      = 1;
  resonance.p.width                  = 0.12;
  resonance.hel_decay.g_decay_TP = {0.31};
  SetToyTensorChannel(resonance, {0.42, -0.17});
  base.process.RESONANCES = {{"f0_pair_sum", resonance}};

  const auto model    = gra::MModelTune::Load(modelfile);
  const auto evaluate = [&](const std::string &channel, const gra::MReggeMode mode) {
    gra::LORENTZSCALAR     lts = base;
    gra::MReggeCentralProc process("TP", channel, {"Tensor pair-source sum", "T-Pomeron", "", "", 401},
                                   gra::ReggeProductionModel::TP, mode);
    process.BindModelTune(model);
    lts.hamp.Configure(TensorPairMetadata(lts.process.FORWARD_NOFLIP));
    process.InitializeProcess(lts);
    const double amp2 = process.Amp2(lts);
    REQUIRE(lts.proton_good_walker.has_value());
    REQUIRE(lts.proton_good_walker->components.size() == 1);
    return std::make_pair(std::move(lts), amp2);
  };

  const auto  resonance_result = evaluate("RES", gra::MReggeMode::Resonance);
  const auto  continuum_result = evaluate("CON", gra::MReggeMode::ContinuumTwoBody);
  const auto  coherent_result  = evaluate("RES+CON", gra::MReggeMode::ResonanceContinuumTwoBody);
  const auto &resonance_lts    = resonance_result.first;
  const auto &continuum_lts    = continuum_result.first;
  const auto &coherent_lts     = coherent_result.first;

  REQUIRE(coherent_lts.hamp.size() == resonance_lts.hamp.size());
  REQUIRE(coherent_lts.hamp.size() == continuum_lts.hamp.size());
  for (const auto &i : gra::aux::indices(coherent_lts.hamp)) {
    RequireComplexNear(coherent_lts.hamp[i], resonance_lts.hamp[i] + continuum_lts.hamp[i], 1.0e-11);
  }

  const auto &resonance_source = resonance_lts.proton_good_walker->components.front().source;
  const auto &continuum_source = continuum_lts.proton_good_walker->components.front().source;
  const auto &coherent_source  = coherent_lts.proton_good_walker->components.front().source;
  REQUIRE(coherent_source.size_row() == resonance_source.size_row());
  REQUIRE(coherent_source.size_col() == resonance_source.size_col());
  REQUIRE(coherent_source.size_row() == continuum_source.size_row());
  REQUIRE(coherent_source.size_col() == continuum_source.size_col());
  for (std::size_t row = 0; row < coherent_source.size_row(); ++row) {
    for (std::size_t col = 0; col < coherent_source.size_col(); ++col) {
      RequireComplexNear(coherent_source[row][col], resonance_source[row][col] + continuum_source[row][col], 1.0e-11);
    }
  }

  const auto projected = gra::ProjectGoodWalker(*coherent_lts.proton_good_walker);
  RequireVectorNear(projected, coherent_lts.hamp, 1.0e-11);
  CHECK(coherent_result.second == Approx(0.25 * gra::SquaredNorm(coherent_lts.hamp)).epsilon(1.0e-11));
}

TEST_CASE("Tensor resonance decay cache matches shifted recomputation", "[MTensorPomeron][RES][cache][screening]") {
  gra::LORENTZSCALAR base = DirectCentralPairLTSForTest(gra::PDG::PDG_pip, gra::PDG::PDG_pim);
  gra::PARAM_RES     resonance;
  resonance.p                        = ToyParticle("f0_decay_cache", 9000022, 0, base.pfinal[0].M());
  resonance.p.P                      = 1;
  resonance.p.C                      = 1;
  resonance.p.width                  = 0.12;
  resonance.hel_decay.g_decay_TP = {0.31};
  SetToyTensorChannel(resonance, {0.42, -0.17});
  base.process.RESONANCES = {{"f0_decay_cache", resonance}};
  const auto model        = gra::MModelTune::Load(modelfile);

  const auto shift_forward = [](gra::LORENTZSCALAR &lts) {
    constexpr double dx = 0.013;
    constexpr double dy = -0.009;
    const auto      &p1 = lts.pfinal[1];
    const auto      &p2 = lts.pfinal[2];
    lts.pfinal[1]       = gra::M4Vec(p1.Px() - dx, p1.Py() - dy, p1.Pz(), p1.E());
    lts.pfinal[2]       = gra::M4Vec(p2.Px() + dx, p2.Py() + dy, p2.Pz(), p2.E());
    RefreshToyDerivedKinematicsPreserveDecay(lts);
  };
  gra::LORENTZSCALAR     cached_lts = base;
  gra::MReggeCentralProc cached_process("TP", "RES", {"Tensor decay cache", "T-Pomeron", "", "", 402},
                                        gra::ReggeProductionModel::TP, gra::MReggeMode::Resonance);
  cached_process.BindModelTune(model);
  cached_lts.hamp.Configure(TensorPairMetadata(cached_lts.process.FORWARD_NOFLIP));
  cached_process.InitializeProcess(cached_lts);
  cached_lts.amplitude.BeginCentral();
  REQUIRE(cached_process.Amp2(cached_lts) > 0.0);
  shift_forward(cached_lts);
  cached_lts.screening.active = true;
  const double cached_amp2      = cached_process.Amp2(cached_lts);
  REQUIRE(std::isfinite(cached_amp2));
  REQUIRE(cached_amp2 > 0.0);
  REQUIRE(cached_lts.proton_good_walker.has_value());
  const auto cached_pair = *cached_lts.proton_good_walker;

  gra::LORENTZSCALAR     direct_lts = base;
  gra::MReggeCentralProc direct_process("TP", "RES", {"Tensor decay cache", "T-Pomeron", "", "", 402},
                                        gra::ReggeProductionModel::TP, gra::MReggeMode::Resonance);
  direct_process.BindModelTune(model);
  direct_lts.hamp.Configure(TensorPairMetadata(direct_lts.process.FORWARD_NOFLIP));
  direct_process.InitializeProcess(direct_lts);
  shift_forward(direct_lts);
  REQUIRE(direct_process.Amp2(direct_lts) > 0.0);
  REQUIRE(direct_lts.proton_good_walker.has_value());
  const auto &direct_pair = *direct_lts.proton_good_walker;

  REQUIRE(cached_pair.components.size() == direct_pair.components.size());
  for (const auto &i : gra::aux::indices(direct_pair.components)) {
    CHECK(cached_pair.components[i].upper_sector == direct_pair.components[i].upper_sector);
    CHECK(cached_pair.components[i].lower_sector == direct_pair.components[i].lower_sector);
    RequireMatrixNear(cached_pair.components[i].source, direct_pair.components[i].source, 2.0e-12);
  }
}

TEST_CASE("MTensorPomeron stores physical QED proton rows negative first",
          "[MTensorPomeron][QED][continuum][helicity][pair-layout][exact]") {
  gra::LORENTZSCALAR full_lts     = DirectCentralPairLTSForTest(-13, 13);
  full_lts.process.FORWARD_NOFLIP = false;
  gra::MTensorPomeron tensor(full_lts, gra::MModelTune::Load(modelfile),
                             gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
  tensor.ME4(full_lts, gra::TensorContinuumMode::QED);
  REQUIRE(full_lts.hamp.size() == 64);

  const auto       upper        = gra::ResolveForwardLegState(full_lts, gra::ForwardBeamLeg::Upper);
  const auto       lower        = gra::ResolveForwardLegState(full_lts, gra::ForwardBeamLeg::Lower);
  const auto       u_a          = tensor.SpinorStates(full_lts.pbeam1, "u");
  const auto       u_b          = tensor.SpinorStates(full_lts.pbeam2, "u");
  const auto       ubar_1       = tensor.SpinorStates(full_lts.pfinal[1], "ubar");
  const auto       ubar_2       = tensor.SpinorStates(full_lts.pfinal[2], "ubar");
  const gra::M4Vec p3           = full_lts.decaytree[0].p4;
  const gra::M4Vec p4           = full_lts.decaytree[1].p4;
  const gra::M4Vec pt           = full_lts.pbeam1 - full_lts.pfinal[1] - p3;
  const gra::M4Vec pu           = p4 - full_lts.pbeam1 + full_lts.pfinal[1];
  const double     mass         = full_lts.decaytree[1].p.mass;
  const auto       propagator_t = tensor.iD_F(pt, mass);
  const auto       propagator_u = tensor.iD_F(pu, mass);
  const auto       v_3          = tensor.SpinorStates(p3, "v");
  const auto       ubar_4       = tensor.SpinorStates(p4, "ubar");

  std::vector<std::complex<double>> expected(64, 0.0);
  for (std::size_t ha = 0; ha < 2; ++ha) {
    for (std::size_t hb = 0; hb < 2; ++hb) {
      for (std::size_t h1 = 0; h1 < 2; ++h1) {
        for (std::size_t h2 = 0; h2 < 2; ++h2) {
          const auto upper_sources = tensor.iG_yForwardSources(full_lts, upper, ubar_1[h1], u_a[ha], ha, h1);
          const auto lower_sources = tensor.iG_yForwardSources(full_lts, lower, ubar_2[h2], u_b[hb], hb, h2);
          REQUIRE(upper_sources.size() == 1);
          REQUIRE(lower_sources.size() == 1);

          const std::size_t row = gra::spin::CanonicalProtonPairSpinLayout::HardRow(ha, hb, h1, h2);
          for (std::size_t h3 = 0; h3 < 2; ++h3) {
            for (std::size_t h4 = 0; h4 < 2; ++h4) {
              const auto           current_t = tensor.iG_yeebary(ubar_4[h4], propagator_t, v_3[h3]);
              const auto           current_u = tensor.iG_yeebary(ubar_4[h4], propagator_u, v_3[h3]);
              std::complex<double> amplitude = 0.0;
              for (std::size_t nu1 = 0; nu1 < 4; ++nu1) {
                for (std::size_t nu2 = 0; nu2 < 4; ++nu2) {
                  amplitude += (-gra::math::zi) * upper_sources[0](nu1) * (current_t(nu2, nu1) + current_u(nu1, nu2)) *
                               lower_sources[0](nu2);
                }
              }
              expected[4 * row + gra::spin::BinaryPairHelicityIndex(h3, h4)] = amplitude;
            }
          }
        }
      }
    }
  }
  RequireVectorNear(full_lts.hamp, expected, 2.0e-12);

  gra::LORENTZSCALAR compact_lts     = full_lts;
  compact_lts.process.FORWARD_NOFLIP = true;
  tensor.ME4(compact_lts, gra::TensorContinuumMode::QED);
  REQUIRE(compact_lts.hamp.size() == 16);
  std::array<std::size_t, 4> full_diagonal_rows{};
  for (std::size_t pair = 0; pair < full_diagonal_rows.size(); ++pair) {
    full_diagonal_rows[pair] = gra::spin::PairHelicityTransitionIndex(pair, pair);
  }
  for (std::size_t row = 0; row < full_diagonal_rows.size(); ++row) {
    for (std::size_t central = 0; central < 4; ++central) {
      CAPTURE(row, central);
      const auto compact = compact_lts.hamp[4 * row + central];
      const auto full    = full_lts.hamp[4 * full_diagonal_rows[row] + central];
      CHECK(std::real(compact) == Approx(std::real(full)).margin(2.0e-12));
      CHECK(std::imag(compact) == Approx(std::imag(full)).margin(2.0e-12));
    }
  }
}

// Reject unsupported resonance spin, parity, charge conjugation and charge
TEST_CASE("MTensorPomeron rejects every unsupported resonance JPC assignment",
          "[MTensorPomeron][physics][quantum-numbers]") {
  struct QuantumNumbers {
    int spinX2;
    int parity;
    int charge_conjugation;
    int chargeX3;
  };
  const std::array<QuantumNumbers, 5> unsupported = {QuantumNumbers{0, 1, -1, 0}, QuantumNumbers{2, 1, -1, 0},
                                                     QuantumNumbers{4, -1, 1, 0}, QuantumNumbers{10, -1, -1, 0},
                                                     QuantumNumbers{0, 1, 1, 3}};

  for (const auto &quantum_numbers : unsupported) {
    CAPTURE(quantum_numbers.spinX2, quantum_numbers.parity, quantum_numbers.charge_conjugation,
            quantum_numbers.chargeX3);
    gra::LORENTZSCALAR lts = TensorRhoCascadeLTSForTest();
    gra::PARAM_RES     res;
    res.p          = ToyParticle("unsupported", 9000999, quantum_numbers.spinX2, lts.pfinal[0].M());
    res.p.P        = quantum_numbers.parity;
    res.p.C        = quantum_numbers.charge_conjugation;
    res.p.chargeX3 = quantum_numbers.chargeX3;
    res.p.width    = 0.1;
    SetToyTensorChannel(res, std::vector<double>(7, 0.1));
    res.hel_decay.g_decay_TP = {0.2, 0.3};
    lts.process.RESONANCES       = {{"unsupported", res}};
    gra::MTensorPomeron tensor(lts, gra::MModelTune::Load(modelfile),
                               gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
    REQUIRE_THROWS(tensor.ME3(lts));
  }
}

TEST_CASE("MTensorPomeron scalar ME3 controls the helicity decay phase", "[MTensorPomeron][physics][phase]") {
  auto make_lts = [](double zeta) {
    gra::LORENTZSCALAR lts = TensorRhoCascadeLTSForTest();
    gra::PARAM_RES     res;
    res.p       = ToyParticle("f0", 9001712, 0, lts.pfinal[0].M());
    res.p.P     = 1;
    res.p.C     = 1;
    res.p.width = 0.15;
    SetToyTensorChannel(res, {0.42, -0.17});
    res.hel_decay.g_decay_TP = {0.31, -0.22};
    res.hel_decay.zeta           = zeta;
    lts.process.RESONANCES       = {{"f0_test", res}};
    return lts;
  };

  const auto disabled_tune = WriteModifiedPhotoVMTune(
      "tensor_scalar_no_zeta", [](auto &j) { j.at("PARAM_TENSORPOM").at("use_zeta") = false; });
  const auto model_tune = gra::MModelTune::Load(disabled_tune.second);
  auto                reference_lts = make_lts(0.0);
  auto                phased_lts    = make_lts(1.234);
  gra::MTensorPomeron reference(reference_lts, model_tune,
                                gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
  gra::MTensorPomeron phased(phased_lts, model_tune,
                             gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
  const double        reference_amp2 = reference.ME3(reference_lts);
  const double        phased_amp2    = phased.ME3(phased_lts);
  RequireVectorNear(phased_lts.hamp, reference_lts.hamp, 1e-12);
  REQUIRE(phased_amp2 == Approx(reference_amp2).epsilon(1e-13));

  const auto zeta_tune = WriteModifiedPhotoVMTune(
      "tensor_scalar_use_zeta", [](auto &j) { j.at("PARAM_TENSORPOM").at("use_zeta") = true; });
  auto enabled_lts = make_lts(1.234);
  gra::MTensorPomeron enabled(enabled_lts, gra::MModelTune::Load(zeta_tune.second),
                              gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
  const double enabled_amp2 = enabled.ME3(enabled_lts);
  const auto   zeta_phase   = std::polar(1.0, 1.234);
  REQUIRE(enabled_lts.hamp.size() == reference_lts.hamp.size());
  for (const auto &i : indices(enabled_lts.hamp)) {
    RequireComplexNear(enabled_lts.hamp[i], zeta_phase * reference_lts.hamp[i], 1e-12);
  }
  REQUIRE(enabled_amp2 == Approx(reference_amp2).epsilon(1e-13));

  auto barrier_lts                     = make_lts(0.0);
  barrier_lts.process.DECAY_BARRIER    = true;
  auto no_barrier_lts                  = barrier_lts;
  no_barrier_lts.process.DECAY_BARRIER = false;
  gra::MTensorPomeron barrier(barrier_lts, gra::MModelTune::Load(modelfile),
                              gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
  gra::MTensorPomeron no_barrier(no_barrier_lts, gra::MModelTune::Load(modelfile),
                                 gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
  const double        barrier_amp2    = barrier.ME3(barrier_lts);
  const double        no_barrier_amp2 = no_barrier.ME3(no_barrier_lts);
  RequireVectorNear(barrier_lts.hamp, no_barrier_lts.hamp, 1e-12);
  REQUIRE(barrier_amp2 == Approx(no_barrier_amp2).epsilon(1e-13));
}

TEST_CASE(
    "MTensorPomeron axial ME3 supports its direct multi-body spin-blind "
    "fallback",
    "[MTensorPomeron][physics][axial]") {
  {
    gra::LORENTZSCALAR lts = AxialThreeBodyLTS();

    gra::PARAM_RES res;
    res.p       = ToyParticle("f1", 9002023, 2, lts.pfinal[0].M());
    res.p.P     = 1;
    res.p.C     = 1;
    res.p.width = 0.12;
    SetToyTensorChannel(res, {0.8, -0.25});
    res.hel_decay.g_decay  = {0.6, 0.1};
    lts.process.RESONANCES = {{"f1_test", res}};

    gra::MTensorPomeron tensor(lts, gra::MModelTune::Load(modelfile),
                               gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
    const double        amp2 = tensor.ME3(lts);
    REQUIRE(std::isfinite(amp2));
    REQUIRE(amp2 > 0.0);
    REQUIRE(!lts.hamp.empty());
  }
}

TEST_CASE("MTensorPomeron processes reject unsupported amplitude topologies", "[MTensorPomeron][topology][process]") {
  const auto qed       = gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::QED);
  const auto continuum = gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Continuum);
  const auto resonance = gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Resonance);
  const auto resonance_continuum =
      gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::ResonanceContinuum);

  const auto muons                = DirectCentralPairLTSForTest(-13, 13).decaytree;
  const auto pions                = DirectCentralPairLTSForTest(gra::PDG::PDG_pip, gra::PDG::PDG_pim).decaytree;
  const auto rho_cascades         = TensorRhoCascadeLTSForTest().decaytree;
  auto       tune_vector_cascades = rho_cascades;
  for (auto &branch : tune_vector_cascades) {
    branch.p.pdg         = 9000113;
    branch.legs[0].p.pdg = 9000211;
    branch.legs[1].p.pdg = -9000211;
  }
  REQUIRE(qed->MatchProcess(muons).has_value());
  REQUIRE_FALSE(qed->MatchProcess(pions).has_value());
  REQUIRE(continuum->MatchProcess(pions).has_value());
  REQUIRE(continuum->MatchProcess(rho_cascades).has_value());
  REQUIRE(continuum->MatchProcess(tune_vector_cascades).has_value());
  REQUIRE_FALSE(continuum->MatchProcess(muons).has_value());
  REQUIRE(resonance->MatchProcess(rho_cascades).has_value());
  REQUIRE_FALSE(resonance->MatchProcess({rho_cascades.front()}).has_value());
  REQUIRE(resonance_continuum->MatchProcess(pions).has_value());
  REQUIRE_FALSE(resonance_continuum->MatchProcess(rho_cascades).has_value());

  auto nested_daughter = rho_cascades;
  nested_daughter[0].legs[0].legs.push_back(nested_daughter[0].legs[1]);
  REQUIRE_FALSE(continuum->MatchProcess(nested_daughter).has_value());
}

// Check scalar and tensor resonance cascades use branch-local decay amplitudes
TEST_CASE("MTensorPomeron resonance cascades require configured vector propagators",
          "[MTensorPomeron][decay][process][resonance]") {
  const auto definition = gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Resonance);
  for (const int resonance_spin_x2 : {0, 4}) {
    CAPTURE(resonance_spin_x2);
    gra::LORENTZSCALAR lts     = TensorRhoCascadeLTSForTest();

    gra::PARAM_RES resonance;
    resonance.p       = ToyParticle("f", 9000021 + resonance_spin_x2, resonance_spin_x2, lts.pfinal[0].M());
    resonance.p.P     = 1;
    resonance.p.C     = 1;
    resonance.p.width = 0.12;
    resonance.hel_decay.g_decay_TP = {0.7, -0.2};
    SetToyTensorChannel(resonance, resonance_spin_x2 == 0 ? std::vector<double>{0.7, 0.2}
                                                         : std::vector<double>{0.7, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
    lts.process.RESONANCES             = {{"f_test", resonance}};
    auto unconfigured = lts;
    unconfigured.decaytree[1].p.pdg = 9000113;
    REQUIRE_THROWS_AS(definition->ResolveProcess(unconfigured), std::invalid_argument);
    auto unresolved = lts;
    unresolved.process.TENSOR_MODEL_READY = false;
    REQUIRE_THROWS_AS(definition->ResolveProcess(unresolved), std::invalid_argument);

    const auto resolved = definition->ResolveProcess(lts);
    REQUIRE(resolved.has_value());
    REQUIRE(resolved->decay_structure.type == gra::DecayType::Full);

    gra::MTensorPomeron tensor(lts, gra::MModelTune::Load(modelfile), definition);
    const double amp2 = tensor.ME3(lts);
    REQUIRE(std::isfinite(amp2));
    REQUIRE(amp2 > 0.0);

    auto missing_coupling = lts;
    missing_coupling.decaytree[0].hel.g_decay_TP.clear();
    CHECK_THROWS_AS(definition->ResolveProcess(missing_coupling), std::invalid_argument);

    for (const bool width : {false, true}) {
      auto invalid = lts;
      auto &vector = invalid.decaytree[0].p;
      (width ? vector.width : vector.mass) = std::numeric_limits<double>::infinity();
      CHECK_THROWS_AS(definition->ResolveProcess(invalid), std::invalid_argument);
    }

    auto zero_daughter_momenta = lts;
    for (auto &branch : zero_daughter_momenta.decaytree) {
      for (auto &daughter : branch.legs) { daughter.p4 = gra::M4Vec(0.0, 0.0, 0.0, 0.0); }
    }
    REQUIRE(definition->ResolveProcess(zero_daughter_momenta).has_value());

    auto stale_daughter_momenta                    = lts;
    stale_daughter_momenta.decaytree[0].legs[1].p4 = gra::M4Vec(0.0, 0.0, 0.0, 0.5);
    REQUIRE(definition->ResolveProcess(stale_daughter_momenta).has_value());

    auto unequal_particle_masses = lts;
    unequal_particle_masses.decaytree[0].legs[1].p.mass += 0.05;
    REQUIRE_THROWS(definition->ResolveProcess(unequal_particle_masses));
  }
}

// Check continuum cascade validation uses tune channels but not current momenta
TEST_CASE("MTensorPomeron continuum process validates particle masses and spins",
          "[MTensorPomeron][decay][process][continuum]") {
  gra::LORENTZSCALAR lts               = TensorRhoCascadeLTSForTest();
  lts.process.TENSOR_MODEL_READY       = true;
  lts.process.TENSOR_VECTOR_DECAY_PDGS = {{std::abs(lts.decaytree[0].p.pdg), std::abs(lts.decaytree[0].legs[0].p.pdg)}};

  const auto definition = gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Continuum);
  REQUIRE(definition->ResolveProcess(lts).has_value());

  for (const bool width : {false, true}) {
    auto invalid = lts;
    auto &vector = invalid.decaytree[0].p;
    (width ? vector.width : vector.mass) = std::numeric_limits<double>::infinity();
    CHECK_THROWS_AS(definition->ResolveProcess(invalid), std::invalid_argument);
  }

  auto zero_daughter_momenta = lts;
  for (auto &branch : zero_daughter_momenta.decaytree) {
    for (auto &daughter : branch.legs) { daughter.p4 = gra::M4Vec(0.0, 0.0, 0.0, 0.0); }
  }
  REQUIRE(definition->ResolveProcess(zero_daughter_momenta).has_value());

  auto stale_daughter_momenta                    = lts;
  stale_daughter_momenta.decaytree[0].legs[1].p4 = gra::M4Vec(0.0, 0.0, 0.0, 0.5);
  REQUIRE(definition->ResolveProcess(stale_daughter_momenta).has_value());

  auto missing_channel = lts;
  missing_channel.process.TENSOR_VECTOR_DECAY_PDGS.clear();
  CHECK_THROWS_AS(definition->ResolveProcess(missing_channel), std::invalid_argument);

  auto unequal_particle_masses = lts;
  unequal_particle_masses.decaytree[0].legs[1].p.mass += 0.05;
  REQUIRE_THROWS(definition->ResolveProcess(unequal_particle_masses));
}

// Check vector resonance production retains its resolved tune-vector channel
TEST_CASE("MTensorPomeron vector resonance process requires tune data", "[MTensorPomeron][decay][process][resonance]") {
  gra::LORENTZSCALAR lts = DirectCentralPairLTSForTest(gra::PDG::PDG_pip, gra::PDG::PDG_pim);
  gra::PARAM_RES     vector;
  vector.p                             = ToyParticle("V", 9000113, 2, lts.pfinal[0].M());
  vector.p.P                           = -1;
  vector.p.C                           = -1;
  vector.p.width                       = 0.12;
  vector.hel_decay.g_decay_TP      = {0.7};
  lts.process.RESONANCES               = {{"V_test", vector}};
  lts.process.TENSOR_MODEL_READY       = true;
  lts.process.TENSOR_VECTOR_DECAY_PDGS = {{std::abs(vector.p.pdg), std::abs(lts.decaytree[0].p.pdg)}};

  const auto definition = gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Resonance);
  REQUIRE(definition->ResolveProcess(lts).has_value());

  auto missing_vector = lts;
  missing_vector.process.TENSOR_VECTOR_DECAY_PDGS.clear();
  CHECK_THROWS_AS(definition->ResolveProcess(missing_vector), std::invalid_argument);
}

TEST_CASE(
    "MTensorPomeron decay structure follows the supported decay "
    "topology",
    "[MTensorPomeron][decay]") {
  const auto definition = gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic);
  const gra::LORENTZSCALAR direct = DirectCentralPairLTSForTest(gra::PDG::PDG_Kp, gra::PDG::PDG_Km);
  const auto direct_decay = definition->DecayStructureFor(direct);
  REQUIRE(direct_decay.type == gra::DecayType::Full);

  gra::LORENTZSCALAR covariant_cascade = TensorRhoCascadeLTSForTest();
  const auto         covariant_decay   = definition->DecayStructureFor(covariant_cascade);
  REQUIRE(covariant_decay.type == gra::DecayType::Full);
  gra::MTensorPomeron                      tensor(covariant_cascade, gra::MModelTune::Load(modelfile),
                                                  gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
  const gra::amplitude::ProcessDefinition &tensor_definition = tensor;
  REQUIRE(tensor_definition.DecayStructureFor(covariant_cascade) == covariant_decay);

}

TEST_CASE(
    "MTensorPomeron axial Jacob-Wick cascade retains proposal-independent complex poles",
    "[MTensorPomeron][decay][axial]") {
  gra::LORENTZSCALAR lts             = DirectCentralPairLTSForTest(gra::PDG::PDG_Kp, gra::PDG::PDG_Km);
  lts.amplitude.DECAY_SYM            = false;
  lts.decay_symmetry_proposal_active = false;

  const gra::M4Vec vector_p4 = lts.decaytree[0].p4;
  lts.decaytree[0] = ContinuumVectorBranchForTest("R", 900001, vector_p4, 900011, 900012, 0.11, 0.13, 0.74, 0.51,
                                                  std::complex<double>(0.72, -0.31), 0.46, 0.08, 1.7);

  gra::PARAM_RES axial;
  axial.p       = ToyParticle("f1", 9002023, 2, lts.pfinal[0].M());
  axial.p.P     = 1;
  axial.p.C     = 1;
  axial.p.width = 0.12;
  SetToyTensorChannel(axial, {0.8, -0.25});
  axial.hel_decay = SpinOneToSpinOneScalarHelicityMatrix(
      {std::complex<double>(0.43, 0.17), std::complex<double>(-0.28, 0.39), std::complex<double>(0.61, -0.22)});

  const std::complex<double> bw = gra::spin::CascadeBWProduct(lts.decaytree);
  REQUIRE(std::abs(bw) > 0.0);

  gra::LORENTZSCALAR flat_lts                 = lts;
  flat_lts.decaytree[0].mass_proposal = gra::MassProposal::Uniform;

  lts.process.RESONANCES      = {{"f1_test", axial}};
  flat_lts.process.RESONANCES = {{"f1_test", axial}};
  gra::MTensorPomeron bw_tensor(lts, gra::MModelTune::Load(modelfile),
                                gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
  gra::MTensorPomeron flat_tensor(flat_lts, gra::MModelTune::Load(modelfile),
                                  gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
  lts.decay_structure = bw_tensor.DecayStructureFor(lts);
  flat_lts.decay_structure = flat_tensor.DecayStructureFor(flat_lts);
  REQUIRE(lts.decay_structure.type == gra::DecayType::JacobWickIncoherent);
  const double        bw_amp2   = bw_tensor.ME3(lts);
  const double        flat_amp2 = flat_tensor.ME3(flat_lts);
  REQUIRE(lts.decay_structure.type == gra::DecayType::JacobWickIncoherent);
  REQUIRE(flat_lts.decay_structure == lts.decay_structure);

  REQUIRE(bw_amp2 > 0.0);
  REQUIRE(flat_amp2 == Approx(bw_amp2).epsilon(1e-11));
  RequireVectorNear(flat_lts.hamp, lts.hamp, 1e-10);
  for (const bool flat : {false, true}) {
    auto mixed = flat ? flat_lts : lts;
    mixed.decay_symmetry_proposal_active = true;
    REQUIRE(bw_tensor.ME3(mixed) == Approx(bw_amp2).epsilon(1e-11));
    RequireVectorNear(mixed.hamp, lts.hamp, 1e-10);
  }
}

// Check TP axial history sums, proposal selection and normalization through the shared declaration
TEST_CASE("TP axial cascades use coherent and incoherent Jacob-Wick declarations",
          "[MTensorPomeron][decay][axial][cascade][proposal]") {
  auto physical = TensorRhoCascadeLTSForTest();
  for (auto &branch : physical.decaytree) {
    branch.hel.alpha_ls.Set(1, 0, 1.0);
    gra::spin::InitTMatrix(branch.hel, branch.p, branch.legs[0].p, branch.legs[1].p,
                          false, "rho decay", false, false);
    branch.hel.g_decay = {0.8, 0.2};
  }
  gra::PARAM_RES axial;
  axial.p = ToyParticle("f1", 9002023, 2, physical.pfinal[0].M());
  axial.p.P = axial.p.C = 1;
  axial.p.width = 0.12;
  SetToyTensorChannel(axial, {0.8, -0.25});
  axial.hel_decay.alpha_ls.Set(2, 4, 1.0);
  gra::spin::InitTMatrix(axial.hel_decay, axial.p, physical.decaytree[0].p, physical.decaytree[1].p,
                        false, "axial vector pair", false, false);
  axial.hel_decay.g_decay = {0.6, 0.1};
  physical.process.RESONANCES = {{"f1_test", axial}};
  const auto definition = gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Resonance);
  gra::MTensorPomeron tensor(physical, gra::MModelTune::Load(modelfile), definition);
  ToyHelicityProcess sampler;

  for (const bool coherent : {false, true}) {
    auto lts = physical;
    lts.amplitude.DECAY_SYM = coherent;
    lts.decay_structure = definition->DecayStructureFor(lts);
    const auto expected_type = coherent ? gra::DecayType::JacobWickCoherent : gra::DecayType::JacobWickIncoherent;
    REQUIRE(lts.decay_structure.type == expected_type);
    const double amp2 = tensor.ME3(lts);
    REQUIRE(amp2 > 0.0);
    REQUIRE(lts.decay_structure.type == expected_type);
    RequireTensorSymmetries(lts, [&](auto &event) { return tensor.ME3(event); });

    sampler.state.lts = lts;
    sampler.ProcPtr.Initialize("TP", "RES");
    sampler.PrepareDecaySymmetryProposal();
    REQUIRE(sampler.state.lts.decay_symmetry_proposal_active == coherent);
    REQUIRE(sampler.DecaySymmetryCompensationFactor() == Approx(coherent ? 1.0 : 2.0));

    std::vector<std::complex<double>> expected(lts.hamp.size(), 0.0);
    const auto terms = gra::spin::StableLeafAmplitudeTrees(lts, lts.decaytree);
    REQUIRE(terms.size() == 2);
    for (std::size_t i = 0; i < (coherent ? terms.size() : 1); ++i) {
      auto term = lts;
      term.decaytree = terms[i].tree;
      term.amplitude.DECAY_SYM = false;
      term.amplitude.BeginCentral();
      REQUIRE(tensor.ME3(term) > 0.0);
      for (const auto &j : indices(expected)) { expected[j] += term.hamp[j] * terms[i].statistics_sign; }
    }
    RequireVectorNear(lts.hamp, expected, 2e-10);

    auto uniform = lts;
    for (auto &branch : uniform.decaytree) { branch.mass_proposal = gra::MassProposal::Uniform; }
    uniform.amplitude.BeginCentral();
    REQUIRE(tensor.ME3(uniform) == Approx(amp2).epsilon(2e-11));
    RequireVectorNear(uniform.hamp, lts.hamp, 2e-10);
  }

  // A third primary branch selects the same incoherent fallback as Regge resonances
  auto multibody = AxialThreeBodyLTS();
  for (std::size_t i = 0; i < 2; ++i) {
    auto branch = physical.decaytree[i];
    const auto next = multibody.decaytree[i].p4;
    const auto daughters = TwoBodyRestKinematics(next.M(), branch.legs[0].p.mass, branch.legs[1].p.mass, 0.74, 0.51);
    branch.p4 = next;
    for (const auto &j : indices(branch.legs)) { branch.legs[j].p4 = BoostFromRestFrame(daughters[j], next); }
    multibody.decaytree[i] = branch;
  }
  multibody.process.RESONANCES = physical.process.RESONANCES;
  multibody.amplitude.DECAY_SYM = true;
  sampler.state.lts = multibody;
  sampler.ProcPtr.Initialize("TP", "RES");
  sampler.PrepareDecaySymmetryProposal();
  REQUIRE(sampler.state.lts.decay_structure.type == gra::DecayType::JacobWickIncoherent);
  REQUIRE_FALSE(sampler.state.lts.decay_symmetry_proposal_active);
  REQUIRE(sampler.DecaySymmetryCompensationFactor() == Approx(2.0));
  const double multibody_amp2 = tensor.ME3(multibody);
  REQUIRE(multibody_amp2 > 0.0);
  const auto multibody_hamp = multibody.hamp;
  multibody.amplitude.DECAY_SYM = false;
  REQUIRE(tensor.ME3(multibody) == Approx(multibody_amp2).epsilon(1e-12));
  RequireVectorNear(multibody.hamp, multibody_hamp, 1e-12);
}

TEST_CASE("MTensorPomeron isolated axial cascades remain production only", "[MTensorPomeron][decay][axial][isolated]") {
  auto make_isolated_lts = [](double left_mass, double right_mass) {
    gra::LORENTZSCALAR lts          = TensorRhoCascadeLTSForTest();
    const auto         vectors_in_X = TwoBodyRestKinematics(lts.pfinal[0].M(), left_mass, right_mass, 0.91, -0.42);
    const gra::M4Vec   left         = BoostFromRestFrame(vectors_in_X[0], lts.pfinal[0]);
    const gra::M4Vec   right        = BoostFromRestFrame(vectors_in_X[1], lts.pfinal[0]);

    const double rho_pole  = 0.77526;
    const double rho_width = 0.1491;
    const double pion_mass = 0.13957061;
    lts.decaytree = {TensorVectorBranchForTest(113, rho_pole, rho_width, left, gra::PDG::PDG_pip, gra::PDG::PDG_pim,
                                               pion_mass, 0.74, 0.51, 11.95),
                     TensorVectorBranchForTest(113, rho_pole, rho_width, right, gra::PDG::PDG_pip, gra::PDG::PDG_pim,
                                               pion_mass, 1.38, -0.81, 11.95)};
    lts.process.root_decay_mode = gra::RootDecayMode::Isolated;
    lts.PS_active               = false;
    lts.amplitude.DECAY_SYM     = false;

    gra::PARAM_RES axial;
    axial.p       = ToyParticle("f1", 9002023, 2, lts.pfinal[0].M());
    axial.p.P     = 1;
    axial.p.C     = 1;
    axial.p.width = 0.12;
    SetToyTensorChannel(axial, {0.8, -0.25});
    lts.process.RESONANCES = {{"f1_test", axial}};
    return lts;
  };

  gra::LORENTZSCALAR first           = make_isolated_lts(0.73, 0.94);
  gra::LORENTZSCALAR second          = make_isolated_lts(0.86, 0.81);
  const auto definition = gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Resonance);
  const auto decay_structure = definition->DecayStructureFor(first);
  REQUIRE(decay_structure.type == gra::DecayType::None);
  const auto resolved   = definition->ResolveProcess(first);
  REQUIRE(resolved.has_value());
  REQUIRE(resolved->decay_structure == decay_structure);

  gra::MTensorPomeron first_tensor(first, gra::MModelTune::Load(modelfile),
                                   gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
  gra::MTensorPomeron second_tensor(second, gra::MModelTune::Load(modelfile),
                                    gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
  const double        first_amp2  = first_tensor.ME3(first);
  const double        second_amp2 = second_tensor.ME3(second);

  REQUIRE(first_amp2 > 0.0);
  REQUIRE(second_amp2 == Approx(first_amp2).epsilon(1e-12));
  RequireVectorNear(second.hamp, first.hamp, 1e-11);
}

TEST_CASE("MTensorPomeron exposes the coherent raw full cascade amplitude", "[MTensorPomeron][tensor][physics]") {
  gra::LORENTZSCALAR  lts  = TensorRhoCascadeLTSForTest();
  const auto          tune = WriteModifiedPhotoVMTune("tensor_raw_cascade", [](auto &j) {
    j.at("PARAM_TENSORPOM").at("FORWARD_NOFLIP") = true;
    j.at("PARAM_TENSORPOM").at("PARAM_CON").at("TP").at("[113,113]")   = {{995, 995}};
  });
  gra::MTensorPomeron tensor(lts, gra::MModelTune::Load(tune.second),
                             gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));

  const std::vector<std::complex<double>> raw  = CoherentRawTensorCascadeAmplitudesForTest(lts, tensor, lts.decaytree);
  const double                            amp2 = tensor.ME6(lts);
  REQUIRE(lts.hamp.size() == raw.size());
  for (std::size_t i = 0; i < raw.size(); ++i) { RequireComplexNear(lts.hamp[i], raw[i], 1e-10); }
  REQUIRE(amp2 == Approx(HelicityAverageAmp2ForTest(lts.hamp)).epsilon(1e-12));
}

TEST_CASE("MTensorPomeron cascade amplitudes do not depend on DECAY_SYM switch", "[MTensorPomeron][tensor][symmetry]") {
  gra::LORENTZSCALAR lts_off  = TensorRhoCascadeLTSForTest();
  gra::LORENTZSCALAR lts_on   = TensorRhoCascadeLTSForTest();
  lts_off.amplitude.DECAY_SYM = false;
  lts_on.amplitude.DECAY_SYM  = true;

  gra::MTensorPomeron tensor(lts_off, gra::MModelTune::Load(modelfile),
                             gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
  const double        amp2_off = tensor.ME6(lts_off);
  const double        amp2_on  = tensor.ME6(lts_on);

  RequireVectorNear(lts_off.hamp, lts_on.hamp, 1e-10);
  REQUIRE(amp2_off == Approx(amp2_on).epsilon(1e-12));
}

TEST_CASE("MTensorPomeron rejects nonconserved vector decay currents", "[MTensorPomeron][tensor][physics][current]") {
  gra::LORENTZSCALAR lts      = TensorRhoCascadeLTSForTest();
  lts.decaytree[0].legs[1].p4 = lts.decaytree[0].legs[1].p4 + gra::M4Vec(0.0, 0.0, 0.0, 0.05);
  gra::MTensorPomeron tensor(lts, gra::MModelTune::Load(modelfile),
                             gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
  REQUIRE_THROWS_AS(tensor.ME6(lts), gra::AmplitudeFailure);
  lts.decaytree[0].legs[1].p4 = gra::M4Vec(0.0, 0.0, 0.0, std::numeric_limits<double>::quiet_NaN());
  REQUIRE_THROWS_AS(tensor.ME6(lts), gra::PhaseSpaceFailure);
}

TEST_CASE("MTensorPomeron cascade daughters follow tensor-vector steering", "[MTensorPomeron][tensor][params]") {
  const auto tune = WriteModifiedPhotoVMTune(
      "tensor_vector_daughters", [](auto &j) {
        auto &vector = j.at("PARAM_TENSORPOM").at("VECTOR");
        vector.at("dPDG").at(0) = 321;
        // A closed pole channel cannot normalize a two-body running width
        vector.at("Wmode").at(0) = "CONSTANT";
      });
  ToyHelicityProcess process;
  process.state.lts = TensorRhoCascadeLTSForTest();
  process.state.lts.process.MP_FRAME = "null";
  process.SetModelTune(gra::MModelTune::Load(tune.second));
  process.SetHelicityConfig(process.GetModelTune());
  process.ProcPtr.Initialize("TP", "CON");
  REQUIRE_THROWS_AS(process.InitializeProcessAmplitude(), std::invalid_argument);
}

// Check the selected Tensor Pomeron process reverses every mixed proposal exactly once
TEST_CASE("Tensor Pomeron process declarations give unbiased full cascade weights",
          "[MTensorPomeron][MProcess][cascade][proposal][physics]") {
  const auto tune = WriteModifiedPhotoVMTune("tensor_mixed_mass_proposals", [](auto &j) {
    j.at("PARAM_TENSORPOM").at("FORWARD_NOFLIP") = true;
    j.at("PARAM_TENSORPOM").at("PARAM_CON").at("TP").at("[113,113]") = {{995, 995}};
  });
  auto point = TensorRhoCascadeLTSForTest();
  point.amplitude.DECAY_SYM = true;
  gra::MTensorPomeron tensor(point, gra::MModelTune::Load(tune.second),
      gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
  ToyHelicityProcess sampler;
  sampler.ProcPtr.Initialize("TP", "CON");
  const double amp2 = tensor.ME6(point);
  REQUIRE(amp2 > 0.0);
  const auto raw = point.hamp;
  for (const bool mixture : {false, true}) {
    for (unsigned int selection = 0; selection < 4; ++selection) {
      CAPTURE(mixture, selection);
      auto lts = point;
      lts.decay_symmetry_proposal_active = mixture;
      unsigned int mask = selection;
      SetMixedMassProposalsForTest(lts.decaytree, mask);
      const auto structure = sampler.ProcPtr.DecayStructureFor(lts);
      REQUIRE(structure.type == gra::DecayType::Full);
      const double actual_amp2 = tensor.ME6(lts);
      REQUIRE(actual_amp2 == Approx(amp2).epsilon(1e-11));
      RequireVectorNear(lts.hamp, raw, 1e-10);
      const double density = mixture ? MixedHistoryDensityForTest(lts) : MixedMassDensityForTest(lts.decaytree);
      const double phase = sampler.CascadePhaseSpaceForTest(lts, structure);
      REQUIRE(density * phase * actual_amp2 ==
              Approx(InternalCascadePhaseSpaceForTest(lts.decaytree) * amp2).epsilon(1e-10));
    }
  }
}

// Exercise real branching initialization before evaluating either mass proposal
TEST_CASE("Tensor cascade initialization accepts both complete-amplitude mass proposals",
          "[MTensorPomeron][MProcess][cascade][proposal][initialization]") {
  double reference = 0.0;
  std::vector<std::complex<double>> reference_helicities;
  for (const bool flat : {false, true}) {
    CAPTURE(flat);
    ToyHelicityProcess sampler;
    sampler.state.lts = TensorRhoCascadeLTSForTest();
    sampler.state.lts.process.MP_FRAME = "null";
    sampler.SetModelTune(sampler.GetModelTune());
    sampler.ProcPtr.Initialize("TP", "CON");
    sampler.state.flat_mass2 = flat;
    auto setup = sampler.CreateProcessSetup();
    REQUIRE_NOTHROW(gra::MTensorPomeron::InitializeBranching(setup, gra::MTensorPomeronMode::Continuum));
    auto &lts = sampler.state.lts;
    const auto structure = sampler.ProcPtr.DecayStructureFor(lts);
    REQUIRE(structure.type == gra::DecayType::Full);
    unsigned int mask = flat ? 3U : 0U;
    SetMixedMassProposalsForTest(lts.decaytree, mask);
    gra::MTensorPomeron tensor(lts, sampler.GetModelTune(),
        gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Continuum), gra::MTensorPomeronMode::Continuum);
    lts.amplitude.DECAY_SYM = true;
    const double value = tensor.ME6(lts);
    REQUIRE(std::isfinite(value));
    REQUIRE(value > 0.0);
    if (!flat) {
      reference = value;
      reference_helicities.assign(lts.hamp.begin(), lts.hamp.end());
    }
    REQUIRE(value == Approx(reference).epsilon(1e-11));
    RequireVectorNear(lts.hamp, reference_helicities, 1e-10);
  }
}

TEST_CASE("MTensorPomeron stable-leaf proposal compensates coherent tensor amplitude",
          "[MTensorPomeron][tensor][symmetry]") {
  gra::LORENTZSCALAR lts             = TensorRhoCascadeLTSForTest();
  lts.decay_symmetry_proposal_active = true;
  lts.decay_structure                = {gra::DecayType::Full};
  const auto          tune           = WriteModifiedPhotoVMTune("tensor_raw_stable_leaf", [](auto &j) {
    j.at("PARAM_TENSORPOM").at("FORWARD_NOFLIP") = true;
    j.at("PARAM_TENSORPOM").at("PARAM_CON").at("TP").at("[113,113]")   = {{995, 995}};
  });
  gra::MTensorPomeron tensor(lts, gra::MModelTune::Load(tune.second),
                             gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));

  const double                            density = gra::decay::MixtureDensity(lts, lts.decaytree);
  const std::vector<std::complex<double>> raw = CoherentRawTensorCascadeAmplitudesForTest(lts, tensor, lts.decaytree);

  const double amp2 = tensor.ME6(lts);
  REQUIRE(lts.hamp.size() == raw.size());
  for (std::size_t i = 0; i < raw.size(); ++i) { RequireComplexNear(lts.hamp[i], raw[i], 1e-10); }
  REQUIRE(amp2 == Approx(HelicityAverageAmp2ForTest(lts.hamp)).epsilon(1e-12));

  ToyHelicityProcess sampler;
  const double       cascade_phase_space = sampler.CascadePhaseSpaceForTest(lts, lts.decay_structure);
  REQUIRE(cascade_phase_space == Approx(InternalCascadePhaseSpaceForTest(lts.decaytree) / density).epsilon(1e-11));
}

// Check that the common sampler alone compensates a flat stable-leaf mixture
TEST_CASE("MProcess flat stable-leaf mixture is compensated exactly once",
          "[MProcess][MTensorPomeron][proposal][flat][symmetry]") {
  gra::LORENTZSCALAR lts             = TensorRhoCascadeLTSForTest();
  lts.decay_symmetry_proposal_active = true;
  lts.decay_structure                = {gra::DecayType::Full};
  for (auto &branch : lts.decaytree) {
    branch.mass_proposal = gra::MassProposal::Uniform;
  }

  const auto terms = gra::spin::StableLeafAmplitudeTrees(lts, lts.decaytree);
  REQUIRE(terms.size() > 1);
  const double reference_phase_space = StableLeafDensityPhaseSpaceForTest(lts, lts.decaytree);
  double       expected_density      = 0.0;
  for (const auto &term : terms) {
    const double proposal_norm    = CascadeMassProposalNormForTest(term.tree);
    const double term_phase_space = StableLeafDensityPhaseSpaceForTest(lts, term.tree);
    REQUIRE(proposal_norm > 0.0);
    REQUIRE(term_phase_space > 0.0);
    expected_density += reference_phase_space / (proposal_norm * term_phase_space);
  }
  expected_density /= static_cast<double>(terms.size());

  const double density = gra::decay::MixtureDensity(lts, lts.decaytree);
  REQUIRE(density == Approx(expected_density).epsilon(1e-12));
  REQUIRE(density != Approx(1.0).epsilon(1e-6));

  const auto                              tune = WriteModifiedPhotoVMTune("tensor_flat_stable_leaf", [](auto &j) {
    j.at("PARAM_TENSORPOM").at("FORWARD_NOFLIP") = true;
    j.at("PARAM_TENSORPOM").at("PARAM_CON").at("TP").at("[113,113]")   = {{995, 995}};
  });
  gra::MTensorPomeron                     tensor(lts, gra::MModelTune::Load(tune.second),
                                                 gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
  const std::vector<std::complex<double>> raw  = CoherentRawTensorCascadeAmplitudesForTest(lts, tensor, lts.decaytree);
  const double                            amp2 = tensor.ME6(lts);
  RequireVectorNear(lts.hamp, raw, 1e-10);
  REQUIRE(amp2 == Approx(HelicityAverageAmp2ForTest(raw)).epsilon(1e-12));

  ToyHelicityProcess sampler;
  const double       base_phase_space    = InternalCascadePhaseSpaceForTest(lts.decaytree);
  const double       cascade_phase_space = sampler.CascadePhaseSpaceForTest(lts, lts.decay_structure);
  REQUIRE(cascade_phase_space == Approx(base_phase_space / density).epsilon(1e-11));
  REQUIRE(cascade_phase_space * amp2 ==
          Approx(base_phase_space * HelicityAverageAmp2ForTest(raw) / density).epsilon(1e-11));
}

TEST_CASE(
    "MTensorPomeron cascade proposal gives the direct 4-body transformed "
    "weight",
    "[MTensorPomeron][tensor][physics]") {
  gra::LORENTZSCALAR  lts  = TensorRhoCascadeLTSForTest();
  const auto          tune = WriteModifiedPhotoVMTune("tensor_raw_direct_weight", [](auto &j) {
    j.at("PARAM_TENSORPOM").at("FORWARD_NOFLIP") = true;
    j.at("PARAM_TENSORPOM").at("PARAM_CON").at("TP").at("[113,113]")   = {{995, 995}};
  });
  gra::MTensorPomeron tensor(lts, gra::MModelTune::Load(tune.second),
                             gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));

  const std::vector<std::complex<double>> raw = CoherentRawTensorCascadeAmplitudesForTest(lts, tensor, lts.decaytree);
  const double                            raw_amp2 = HelicityAverageAmp2ForTest(raw);
  const std::complex<double>              bw_ref   = gra::spin::CascadeBWProduct(lts.decaytree);
  REQUIRE(std::abs(bw_ref) > 0.0);

  const double full_amp2 = tensor.ME6(lts);
  ToyHelicityProcess sampler;
  const double cascade_weight = lts.DW.Integral() / (2.0 * gra::math::PI) *
                                sampler.CascadePhaseSpaceForTest(lts, {gra::DecayType::Full}) * full_amp2;

  const double direct_four_body_density = TensorDirectFourBodyPhaseSpaceForTest(lts) * raw_amp2;
  const double bw_proposal_transform =
      lts.decaytree[0].mass_proposal_norm * lts.decaytree[1].mass_proposal_norm / gra::math::abs2(bw_ref);
  const double transformed_direct_weight = direct_four_body_density * bw_proposal_transform;

  REQUIRE(cascade_weight == Approx(transformed_direct_weight).epsilon(1e-10));
}

TEST_CASE("MTensorPomeron scalar VV cascade uses coherent vector decay vertices", "[MTensorPomeron][tensor][physics]") {
  gra::LORENTZSCALAR lts = TensorRhoCascadeLTSForTest();

  gra::PARAM_RES res;
  res.p       = ToyParticle("f0", 9001710, 0, lts.pfinal[0].M());
  res.p.P     = 1;
  res.p.C     = 1;
  res.p.width = 0.15;
  SetToyTensorChannel(res, {0.42, -0.17});
  res.hel_decay.g_decay_TP = {0.31, -0.22};
  lts.process.RESONANCES       = {{"f0_test", res}};

  const auto          tune = WriteModifiedPhotoVMTune("tensor_raw_scalar_cascade",
                                                      [](auto &j) { j.at("PARAM_TENSORPOM").at("FORWARD_NOFLIP") = true; });
  gra::MTensorPomeron tensor(lts, gra::MModelTune::Load(tune.second),
                             gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
  const auto          terms = gra::spin::StableLeafAmplitudeTrees(lts, lts.decaytree);
  REQUIRE(terms.size() == 2);
  std::vector<std::complex<double>> raw;
  for (const auto &term : terms) {
    gra::LORENTZSCALAR term_lts = lts;
    term_lts.decaytree          = term.tree;
    const auto term_raw         = ScalarRawCascadeAmplitudesForTest(tensor, term_lts, res);
    if (raw.empty()) { raw.assign(term_raw.size(), 0.0); }
    REQUIRE(raw.size() == term_raw.size());
    for (std::size_t i = 0; i < raw.size(); ++i) { raw[i] += term.statistics_sign * term_raw[i]; }
  }
  const double amp2 = tensor.ME3(lts);
  REQUIRE(amp2 > 0.0);
  REQUIRE(lts.hamp.size() == raw.size());
  for (std::size_t i = 0; i < raw.size(); ++i) { RequireComplexNear(lts.hamp[i], raw[i], 1e-10); }
  REQUIRE(amp2 == Approx(HelicityAverageAmp2ForTest(lts.hamp)).epsilon(1e-12));

  ToyHelicityProcess sampler;
  sampler.ProcPtr.Initialize("TP", "RES");
  for (const bool mixture : {false, true}) {
    for (unsigned int selection = 0; selection < 4; ++selection) {
      CAPTURE(mixture, selection);
      auto sampled = lts;
      sampled.decay_symmetry_proposal_active = mixture;
      unsigned int mask = selection;
      SetMixedMassProposalsForTest(sampled.decaytree, mask);
      const auto structure = sampler.ProcPtr.DecayStructureFor(sampled);
      REQUIRE(structure.type == gra::DecayType::Full);
      const double sampled_amp2 = tensor.ME3(sampled);
      REQUIRE(sampled_amp2 == Approx(amp2).epsilon(1e-11));
      RequireVectorNear(sampled.hamp, raw, 1e-10);
      const double density = mixture ? MixedHistoryDensityForTest(sampled) : MixedMassDensityForTest(sampled.decaytree);
      const double phase = sampler.CascadePhaseSpaceForTest(sampled, structure);
      REQUIRE(density * phase * sampled_amp2 ==
              Approx(InternalCascadePhaseSpaceForTest(sampled.decaytree) * amp2).epsilon(1e-10));
    }
  }
}

TEST_CASE("MTensorPomeron ME3 uses only its tensor-local forward-spin switch", "[MTensorPomeron][tensor][helicity]") {
  for (const bool tensor_noflip : {true, false}) {
    CAPTURE(tensor_noflip);
    const auto tune = WriteModifiedPhotoVMTune(
        "tensor_me3_forward_" + std::to_string(tensor_noflip),
        [tensor_noflip](auto &j) { j.at("PARAM_TENSORPOM").at("FORWARD_NOFLIP") = tensor_noflip; });

    // Evaluate identical tensor physics with opposite generic Regge spin
    // switches
    auto evaluate = [&tune](bool generic_noflip) {
      gra::LORENTZSCALAR lts     = TensorRhoCascadeLTSForTest();
      lts.process.FORWARD_NOFLIP = generic_noflip;

      gra::PARAM_RES res;
      res.p       = ToyParticle("f0", 9001711, 0, lts.pfinal[0].M());
      res.p.P     = 1;
      res.p.C     = 1;
      res.p.width = 0.15;
      SetToyTensorChannel(res, {0.42, -0.17});
      res.hel_decay.g_decay_TP = {0.31, -0.22};
      lts.process.RESONANCES       = {{"f0_test", res}};

      gra::MTensorPomeron tensor(lts, gra::MModelTune::Load(tune.second),
                                 gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
      const double        amp2 = tensor.ME3(lts);
      return std::make_pair(lts.hamp, amp2);
    };

    const auto        generic_full   = evaluate(false);
    const auto        generic_noflip = evaluate(true);
    const std::size_t expected_rows  = tensor_noflip ? 4 : 16;
    REQUIRE(generic_full.first.size() == expected_rows);
    REQUIRE(generic_noflip.first.size() == expected_rows);
    RequireVectorNear(generic_full.first, generic_noflip.first, 1e-13);
    REQUIRE(generic_full.second == Approx(generic_noflip.second).epsilon(1e-13));
    REQUIRE(generic_full.second == Approx(HelicityAverageAmp2ForTest(generic_full.first)).epsilon(1e-12));
  }
}

TEST_CASE("MTensorPomeron ME6 uses only its tensor-local forward-spin switch",
          "[MTensorPomeron][tensor][cascade][helicity]") {
  for (const bool tensor_noflip : {true, false}) {
    CAPTURE(tensor_noflip);
    const auto tune = WriteModifiedPhotoVMTune(
        "tensor_me6_forward_" + std::to_string(tensor_noflip),
        [tensor_noflip](auto &j) { j.at("PARAM_TENSORPOM").at("FORWARD_NOFLIP") = tensor_noflip; });

    // Evaluate the same cascade with opposite generic Regge spin switches
    auto evaluate = [&tune](bool generic_noflip) {
      gra::LORENTZSCALAR lts     = TensorRhoCascadeLTSForTest();
      lts.process.FORWARD_NOFLIP = generic_noflip;
      gra::MTensorPomeron tensor(lts, gra::MModelTune::Load(tune.second),
                                 gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
      const double        amp2 = tensor.ME6(lts);
      return std::make_pair(lts.hamp, amp2);
    };

    const auto        generic_full   = evaluate(false);
    const auto        generic_noflip = evaluate(true);
    const std::size_t expected_rows  = tensor_noflip ? 4 : 16;
    REQUIRE(generic_full.first.size() == expected_rows);
    REQUIRE(generic_noflip.first.size() == expected_rows);
    RequireVectorNear(generic_full.first, generic_noflip.first, 1e-13);
    REQUIRE(generic_full.second == Approx(generic_noflip.second).epsilon(1e-13));
    REQUIRE(generic_full.second == Approx(HelicityAverageAmp2ForTest(generic_full.first)).epsilon(1e-12));
  }
}



// Compare pseudoscalar fusion to the covariant Levi-Civita contraction
TEST_CASE(
    "gra::spin:: pseudoscalar fusion matches covariant epsilon helicity "
    "structure",
    "[gra::spin][MTensorPomeron]") {
  gra::LORENTZSCALAR tensor_lts;
  tensor_lts.PDG = LoadedPDGTable();
  gra::MTensorPomeron tensor(tensor_lts, gra::MModelTune::Load(modelfile),
                             gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));

  const double     mass = 0.9;
  const gra::M4Vec q1(0.31, -0.17, 1.13, std::sqrt(mass * mass + 0.31 * 0.31 + 0.17 * 0.17 + 1.13 * 1.13));
  const gra::M4Vec q2(-q1.Px(), -q1.Py(), -q1.Pz(), q1.E());

  SECTION("spin-1 spin-1 uses the vector-vector pseudoscalar epsilon tensor") {
    const auto hel       = PseudoscalarFusionHelicityMatrix(993, 2, -1, 1, 2);
    const auto jw        = MHelicityScatterCentralAmplitudes(hel, q1, q2);
    const auto covariant = ProjectVectorPseudoscalarVertex(tensor, q1, q2);
    RequireProjectivelyEqual(jw, covariant);
  }

  SECTION(
      "spin-2 spin-2 tensor-pomeron epsilon structures stay in the "
      "pseudoscalar LS subspace") {
    const auto hel_11      = PseudoscalarFusionHelicityMatrix(995, 4, 1, 1, 2);
    const auto hel_33      = PseudoscalarFusionHelicityMatrix(995, 4, 1, 3, 6);
    const auto jw_11       = MHelicityScatterCentralAmplitudes(hel_11, q1, q2);
    const auto jw_33       = MHelicityScatterCentralAmplitudes(hel_33, q1, q2);
    const auto covariant_0 = ProjectTensorPseudoscalarVertex(tensor, q1, q2, 0);
    const auto covariant_1 = ProjectTensorPseudoscalarVertex(tensor, q1, q2, 1);
    RequireInTwoVectorSpan(covariant_0, jw_11, jw_33);
    RequireInTwoVectorSpan(covariant_1, jw_11, jw_33);

    const auto [c00, c01]                          = TwoVectorSpanCoefficients(covariant_0, jw_11, jw_33);
    const auto [c10, c11]                          = TwoVectorSpanCoefficients(covariant_1, jw_11, jw_33);
    const std::complex<double> change_of_basis_det = c00 * c11 - c01 * c10;

    CAPTURE(c00, c01, c10, c11, change_of_basis_det);
    REQUIRE(std::abs(change_of_basis_det) > 1e-10);
  }
}

TEST_CASE("Tensor-Pomeron axial vertices match the published reduced amplitudes",
          "[gra::spin][MTensorPomeron][axial]") {
  gra::LORENTZSCALAR tensor_lts;
  tensor_lts.PDG = LoadedPDGTable();
  gra::MTensorPomeron tensor(tensor_lts, gra::MModelTune::Load(modelfile),
                             gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));

  const double     exchange_mass = 0.9;
  const gra::M4Vec q1(0.0, 0.0, 1.13, std::sqrt(pow2(exchange_mass) + 1.13 * 1.13));
  const gra::M4Vec q2(-q1.Px(), -q1.Py(), -q1.Pz(), q1.E());
  const auto       g22         = tensor.iG_PPA_22(q1, q2, 1.0);
  const auto       g44         = tensor.iG_PPA_44(q1, q2, 1.0);
  const auto       g22_swapped = tensor.iG_PPA_22(q2, q1, 1.0);
  const auto       g44_swapped = tensor.iG_PPA_44(q2, q1, 1.0);
  RequireAxialVertexIdentities(g22, q1 + q2);
  RequireAxialVertexIdentities(g44, q1 + q2);

  for (std::size_t kappa = 0; kappa < 4; ++kappa) {
    for (std::size_t lambda = 0; lambda < 4; ++lambda) {
      for (std::size_t rho = 0; rho < 4; ++rho) {
        for (std::size_t sigma = 0; sigma < 4; ++sigma) {
          for (std::size_t alpha = 0; alpha < 4; ++alpha) {
            const std::vector<std::size_t> direct  = {kappa, lambda, rho, sigma, alpha};
            const std::vector<std::size_t> swapped = {rho, sigma, kappa, lambda, alpha};
            RequireComplexNear(g22(direct), g22_swapped(swapped), 1e-12);
            RequireComplexNear(g44(direct), g44_swapped(swapped), 1e-12);
          }
        }
      }
    }
  }

  const auto covariant22 = ProjectTensorAxialVertex(tensor, q1, q2, 0);
  const auto covariant44 = ProjectTensorAxialVertex(tensor, q1, q2, 1);
  const auto paper22     = PaperTensorAxialReducedTable(exchange_mass, q1.P3mod(), 0);
  const auto paper44     = PaperTensorAxialReducedTable(exchange_mass, q1.P3mod(), 1);
  RequireProjectivelyEqual(covariant22, paper22, 1e-9);
  RequireProjectivelyEqual(covariant44, paper44, 1e-9);

  const gra::M4Vec equal_q(0.2, -0.1, 0.3, 1.2);
  const auto       equal22 = tensor.iG_PPA_22(equal_q, equal_q, 1.0);
  const auto       equal44 = tensor.iG_PPA_44(equal_q, equal_q, 1.0);
  for (std::size_t kappa = 0; kappa < 4; ++kappa) {
    for (std::size_t lambda = 0; lambda < 4; ++lambda) {
      for (std::size_t rho = 0; rho < 4; ++rho) {
        for (std::size_t sigma = 0; sigma < 4; ++sigma) {
          for (std::size_t alpha = 0; alpha < 4; ++alpha) {
            const std::vector<std::size_t> index = {kappa, lambda, rho, sigma, alpha};
            CHECK(std::abs(equal22(index)) < 1e-14);
            CHECK(std::abs(equal44(index)) < 1e-14);
          }
        }
      }
    }
  }
}

// Validate the two alternative WA102 f1(1420) coupling benchmarks and form factor
// [REFERENCE: Lebiedowicz et al., Phys. Rev. D 102, 114003 (2020), Eqs. (2.15), (3.15), (3.16)]
TEST_CASE("TUNE0 f1(1420) uses a published tensor-Pomeron benchmark", "[gra::spin][MTensorPomeron][axial][card]") {
  ModelParamRestoreGuard restore;
  gra::MODELPARAM = "TUNE0";
  MRandom    rng;
  const auto f1 = gra::resonance::Read("RES/f1_1420.json", rng, gra::ReggeProductionModel::TP);

  REQUIRE(f1.TP.channels.size() == 1);
  const auto &f1_channel = ToyTensorChannel(f1);
  REQUIRE(f1_channel.g_tensor.size() == 2);

  gra::LORENTZSCALAR tensor_lts;
  tensor_lts.PDG = LoadedPDGTable();
  gra::MTensorPomeron tensor(tensor_lts, gra::MModelTune::Load(modelfile),
                             gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
  const gra::M4Vec    q1(0.31, -0.13, 0.42, 0.11);
  const gra::M4Vec    q2(-0.19, 0.07, -0.28, 0.09);
  REQUIRE(q1.M2() < 0.0);
  REQUIRE(q2.M2() < 0.0);

  constexpr double lambda_e          = 0.7;
  const double     paper_form_factor = std::exp((q1.M2() + q2.M2()) / pow2(lambda_e));
  const auto       bare22            = tensor.iG_PPA_22(q1, q2, 1.0);
  const auto       bare44            = tensor.iG_PPA_44(q1, q2, 1.0);
  auto             fitted22_channel  = f1_channel;
  fitted22_channel.ff_prod           = {};
  const auto fitted22                = tensor.iG_PPA_total(q1, q2, f1.p.mass, fitted22_channel);

  auto fitted44_channel     = fitted22_channel;
  fitted44_channel.g_tensor = {0.0, 4.20};
  const auto fitted44       = tensor.iG_PPA_total(q1, q2, f1.p.mass, fitted44_channel);

  for (std::size_t kappa = 0; kappa < 4; ++kappa) {
    for (std::size_t lambda = 0; lambda < 4; ++lambda) {
      for (std::size_t rho = 0; rho < 4; ++rho) {
        for (std::size_t sigma = 0; sigma < 4; ++sigma) {
          for (std::size_t alpha = 0; alpha < 4; ++alpha) {
            const std::vector<std::size_t> index = {kappa, lambda, rho, sigma, alpha};
            RequireComplexNear(fitted22(index), 2.39 * paper_form_factor * bare22(index), 1e-12);
            RequireComplexNear(fitted44(index), 4.20 * paper_form_factor * bare44(index), 1e-12);
          }
        }
      }
    }
  }
}

TEST_CASE("Tensor continuum reduced contractions reproduce covariant amplitudes",
          "[MTensorPomeron][continuum][contraction][physics][exact]") {
  const std::array<std::array<int, 2>, 6> channels = {{
      {211, -211},
      {321, -321},
      {-2212, 2212},
      {-3122, 3122},
      {113, 113},
      {333, 333},
  }};

  for (const bool forward_noflip : {true, false}) {
    const auto tune       = WriteModifiedPhotoVMTune("tensor_continuum_contraction_" + std::to_string(forward_noflip),
                                                     [forward_noflip](auto &j) {
                                                 j.at("PARAM_TENSORPOM").at("FORWARD_NOFLIP") = forward_noflip;
                                                 for (auto &[key, pairs] : j.at("PARAM_TENSORPOM").at("PARAM_CON").at("TP").items()) {
                                                   (void)key;
                                                   pairs = {{995, 995}};
                                                 }
                                               }, [](auto &j) {
                                                 j.at("995").at("[333,9993]").at("g_tensor") = {0.0, 0.0};
                                               });
    const auto soft_model = gra::MModelTune::Load(tune.second);

    for (const auto &channel : channels) {
      CAPTURE(forward_noflip, channel[0], channel[1]);
      gra::LORENTZSCALAR  lts = DirectCentralPairLTSForTest(channel[0], channel[1]);
      gra::MTensorPomeron tensor(lts, soft_model,
                                 gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
      const auto          parameters = gra::ReadTensorPomeronParam(*soft_model, lts.PDG);
      const auto          reference  = TensorContinuumReferenceAmplitudes(tensor, lts, *parameters);
      const double        amp2       = tensor.ME4(lts, gra::TensorContinuumMode::TensorPomeron);
      RequireTensorSymmetries(lts, [&](auto &event) { return tensor.ME4(event, gra::TensorContinuumMode::TensorPomeron); });

      REQUIRE(lts.hamp.size() == reference.size());
      for (const auto &i : gra::aux::indices(reference)) { RequireComplexNear(lts.hamp[i], reference[i], 1.0e-10); }
      REQUIRE(amp2 == Approx(0.25 * gra::SquaredNorm(reference)).epsilon(1.0e-10).margin(1.0e-12));
    }
  }
}

// Check coherent resonance summation across unequal rank-two exchanges
TEST_CASE("Tensor resonance channels add as coherent covariant amplitudes",
          "[MTensorPomeron][resonance][exchange][physics]") {
  auto evaluate = [](const std::vector<gra::RES_TENSOR_CHANNEL> &channels) {
    gra::LORENTZSCALAR lts = DirectCentralPairLTSForTest(211, -211);
    gra::PARAM_RES     resonance;
    resonance.p                        = ToyParticle("f0_exchange", 9001710, 0, lts.pfinal[0].M());
    resonance.p.P                      = 1;
    resonance.p.C                      = 1;
    resonance.p.width                  = 0.15;
    resonance.TP.channels              = channels;
    resonance.hel_decay.g_decay_TP = {0.31};
    lts.process.RESONANCES             = {{"f0_exchange", resonance}};
    gra::MTensorPomeron tensor(lts, gra::MModelTune::Load(modelfile),
                               gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
    tensor.ME3(lts);
    return lts.hamp;
  };

  gra::RES_TENSOR_CHANNEL pp;
  pp.exchange                 = {995, 995};
  pp.g_tensor                 = {0.37, -0.21};
  gra::RES_TENSOR_CHANNEL f2p = pp;
  f2p.exchange                = {9915, 995};
  f2p.g_tensor                = {-0.14, 0.09};
  const auto pp_amp           = evaluate({pp});
  const auto f2p_amp          = evaluate({f2p});
  const auto coherent         = evaluate({pp, f2p});
  REQUIRE(coherent.size() == pp_amp.size());
  REQUIRE(coherent.size() == f2p_amp.size());
  for (const auto &i : indices(coherent)) { RequireComplexNear(coherent[i], pp_amp[i] + f2p_amp[i], 1.0e-10); }
}

// Check every published secondary exchange enters the baryon continuum
TEST_CASE("Tensor continuum evaluates every configured secondary exchange",
          "[MTensorPomeron][continuum][Reggeon][Odderon][physics]") {
  const std::array<std::array<int, 2>, 5> pairs = {std::array<int, 2>{995, 9915}, std::array<int, 2>{995, 9925},
                                                   std::array<int, 2>{995, 9993}, std::array<int, 2>{995, 9933},
                                                   std::array<int, 2>{995, 9943}};
  for (const auto &pair : pairs) {
    CAPTURE(pair[0], pair[1]);
    const auto          tune = WriteModifiedPhotoVMTune("tensor_secondary_" + std::to_string(pair[1]), [pair](auto &j) {
      j.at("PARAM_TENSORPOM").at("FORWARD_NOFLIP")  = true;
      j.at("PARAM_TENSORPOM").at("PARAM_CON").at("TP").at("[2212,-2212]") = {pair};
    });
    gra::LORENTZSCALAR  lts  = DirectCentralPairLTSForTest(-2212, 2212);
    gra::MTensorPomeron tensor(lts, gra::MModelTune::Load(tune.second),
                               gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
    const double        amp2 = tensor.ME4(lts, gra::TensorContinuumMode::TensorPomeron);
    REQUIRE(std::isfinite(amp2));
    REQUIRE(amp2 > 0.0);
  }
}

// Check the published baryon pair charge conjugation exchange relations
TEST_CASE("Tensor baryon continuum has the published C exchange symmetry",
          "[MTensorPomeron][continuum][baryon][C-parity][physics]") {
  const auto evaluate = [](const std::array<int, 2> &exchange, const bool swap_momenta, const std::string &label) {
    const auto         tune = WriteModifiedPhotoVMTune(label, [exchange](auto &j) {
      j.at("PARAM_TENSORPOM").at("FORWARD_NOFLIP")  = true;
      j.at("PARAM_TENSORPOM").at("PARAM_CON").at("TP").at("[2212,-2212]") = {exchange};
    });
    gra::LORENTZSCALAR lts  = DirectCentralPairLTSForTest(-2212, 2212);
    if (swap_momenta) { std::swap(lts.decaytree[0].p4, lts.decaytree[1].p4); }
    gra::MTensorPomeron tensor(lts, gra::MModelTune::Load(tune.second),
                               gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
    REQUIRE(tensor.ME4(lts, gra::TensorContinuumMode::TensorPomeron) > 0.0);
    return std::vector<std::complex<double>>(lts.hamp.begin(), lts.hamp.end());
  };

  const std::array<int, 2> even_exchange = {995, 995};
  const auto               even          = evaluate(even_exchange, false, "tensor_baryon_c_even");
  const auto               even_swapped  = evaluate(even_exchange, true, "tensor_baryon_c_even_swapped");
  REQUIRE(even.size() == 16);
  REQUIRE(even_swapped.size() == even.size());

  struct ExchangeCase {
    std::array<int, 2> exchange;
    int                C;
  };
  const std::array<ExchangeCase, 5> cases = {{
      {{995, 9915}, 1},
      {{9933, 9933}, 1},
      {{995, 9993}, -1},
      {{995, 9933}, -1},
      {{995, 9943}, -1},
  }};
  for (const auto &test : cases) {
    CAPTURE(test.exchange[0], test.exchange[1], test.C);
    const std::string stem =
        "tensor_baryon_c_" + std::to_string(test.exchange[0]) + "_" + std::to_string(test.exchange[1]);
    const auto amplitude = evaluate(test.exchange, false, stem);
    const auto swapped   = evaluate(test.exchange, true, stem + "_swapped");
    REQUIRE(amplitude.size() == even.size());
    REQUIRE(swapped.size() == even.size());

    for (std::size_t forward = 0; forward < 4; ++forward) {
      double               even_norm            = 0.0;
      double               even_swapped_norm    = 0.0;
      double               norm                 = 0.0;
      double               swapped_norm         = 0.0;
      std::complex<double> interference         = 0.0;
      std::complex<double> swapped_interference = 0.0;
      for (const auto h3 : gra::spin::BinaryHelicityIndices()) {
        for (const auto h4 : gra::spin::BinaryHelicityIndices()) {
          const auto direct_index  = 4 * forward + gra::spin::BinaryPairHelicityIndex(h3, h4);
          const auto swapped_index = 4 * forward + gra::spin::BinaryPairHelicityIndex(h4, h3);
          even_norm += std::norm(even[direct_index]);
          even_swapped_norm += std::norm(even_swapped[swapped_index]);
          norm += std::norm(amplitude[direct_index]);
          swapped_norm += std::norm(swapped[swapped_index]);
          interference += even[direct_index] * std::conj(amplitude[direct_index]);
          swapped_interference += even_swapped[swapped_index] * std::conj(swapped[swapped_index]);
        }
      }
      REQUIRE(even_swapped_norm == Approx(even_norm).epsilon(2.0e-10));
      REQUIRE(swapped_norm == Approx(norm).epsilon(2.0e-10));
      RequireComplexNear(swapped_interference, static_cast<double>(test.C) * interference, 2.0e-10);
    }
  }
}

// Check exact coherent sums across every supported outer exchange pair
TEST_CASE("Tensor continuum outer exchange pairs add coherently",
          "[MTensorPomeron][continuum][exchange][coherence][physics]") {
  struct TensorContinuumCase {
    const char                     *label;
    const char                     *key;
    std::array<int, 2>              final_state;
    std::vector<std::array<int, 2>> exchange;
  };
  const std::array<TensorContinuumCase, 3> cases = {{
      {"kaon",
       "[321,-321]",
       {321, -321},
       {{995, 995},
        {995, 9915},
        {995, 9925},
        {995, 9933},
        {995, 9943},
        {9915, 9915},
        {9915, 9925},
        {9915, 9933},
        {9915, 9943},
        {9925, 9925},
        {9925, 9933},
        {9925, 9943},
        {9933, 9933},
        {9933, 9943},
        {9943, 9943}}},
      {"baryon",
       "[2212,-2212]",
       {-2212, 2212},
       {{995, 995},
        {995, 9915},
        {995, 9925},
        {995, 9993},
        {995, 9933},
        {995, 9943},
        {9915, 9915},
        {9915, 9925},
        {9915, 9933},
        {9915, 9943},
        {9925, 9925},
        {9925, 9933},
        {9925, 9943},
        {9933, 9933},
        {9933, 9943},
        {9943, 9943}}},
      {"vector", "[113,113]", {113, 113}, {{995, 995}, {995, 9915}}},
  }};

  struct TensorRun {
    std::vector<std::complex<double>>  amplitude;
    gra::MMatrix<std::complex<double>> source;
  };
  const auto evaluate = [](const TensorContinuumCase &test, const std::vector<std::array<int, 2>> &exchange,
                           const std::string &suffix) {
    const auto         tune       = WriteModifiedPhotoVMTune("tensor_exchange_sum_" + suffix, [&](auto &j) {
      j.at("PARAM_TENSORPOM").at("FORWARD_NOFLIP") = true;
      j.at("PARAM_TENSORPOM").at("PARAM_CON").at("TP").at(test.key)      = exchange;
    });
    const auto         soft_model = gra::MModelTune::Load(tune.second);
    gra::LORENTZSCALAR lts        = DirectCentralPairLTSForTest(test.final_state[0], test.final_state[1]);
    lts.hamp.Configure(TensorPairMetadata(true));
    gra::MTensorPomeron tensor(lts, soft_model,
                               gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
    const double        amp2 = tensor.ME4(lts, gra::TensorContinuumMode::TensorPomeron);
    REQUIRE(std::isfinite(amp2));
    REQUIRE(amp2 > 0.0);
    REQUIRE(lts.proton_good_walker.has_value());
    REQUIRE(lts.proton_good_walker->components.size() == 1);
    return TensorRun{std::vector<std::complex<double>>(lts.hamp.begin(), lts.hamp.end()),
                     lts.proton_good_walker->components.front().source};
  };

  for (const auto &test : cases) {
    CAPTURE(test.label);
    const auto                         full = evaluate(test, test.exchange, std::string(test.label) + "_full");
    std::vector<std::complex<double>>  amplitude_sum(full.amplitude.size(), 0.0);
    gra::MMatrix<std::complex<double>> source_sum(full.source.size_row(), full.source.size_col(), 0.0);
    for (const auto &channel : indices(test.exchange)) {
      const auto isolated =
          evaluate(test, {test.exchange[channel]}, std::string(test.label) + "_" + std::to_string(channel));
      REQUIRE(isolated.amplitude.size() == amplitude_sum.size());
      REQUIRE(isolated.source.size_row() == source_sum.size_row());
      REQUIRE(isolated.source.size_col() == source_sum.size_col());
      gra::AddScaled(amplitude_sum, isolated.amplitude, std::complex<double>(1.0, 0.0));
      source_sum += isolated.source;
    }
    RequireVectorNear(full.amplitude, amplitude_sum, 3.0e-10);
    RequireMatrixNear(full.source, source_sum, 3.0e-10);
  }

  const auto forward = evaluate(cases[0], {{995, 9925}}, "kaon_forward");
  const auto reverse = evaluate(cases[0], {{9925, 995}}, "kaon_reverse");
  RequireVectorNear(forward.amplitude, reverse.amplitude, 3.0e-11);
  RequireMatrixNear(forward.source, reverse.source, 3.0e-11);
}

// Check vector Odderon and tensor exchange fusion for a vector resonance
TEST_CASE("Vector Odderon tensor fusion produces the complete phi amplitude",
          "[MTensorPomeron][resonance][Odderon][physics]") {
  gra::LORENTZSCALAR lts = AsymmetricCentralPairLTSForTest(321, -321, 0.91, -0.38);
  gra::PARAM_RES     phi;
  phi.p = lts.PDG.FindByPDG(333);
  gra::RES_TENSOR_CHANNEL channel;
  channel.exchange             = {995, 9993};
  channel.g_tensor             = {-0.8, 1.6};
  channel.ff_transfer          = gra::regge::FFParam{gra::regge::FFType::Power, gra::regge::FFNorm::Zero, {0.5, 1.0}};
  phi.TP.channels              = {channel};
  phi.hel_decay.g_decay_TP = {4.48};
  lts.process.RESONANCES       = {{"phi_odderon", phi}};
  gra::MTensorPomeron tensor(lts, gra::MModelTune::Load(modelfile),
                             gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
  const double        amp2 = tensor.ME3(lts);
  REQUIRE(std::isfinite(amp2));
  REQUIRE(amp2 > 0.0);
  REQUIRE(gra::SquaredNorm(lts.hamp) > 0.0);
}

// Check coherent meson and vector Odderon transfer in stable and decaying phi pairs
TEST_CASE("Tensor phi-pair production sums meson and Odderon transfers",
          "[MTensorPomeron][continuum][Odderon][cascade][physics]") {
  const bool cascade = GENERATE(false, true);
  const bool noflip = GENERATE(false, true);
  auto evaluate = [&](const std::string &label, const bool meson, const bool odderon) {
    const auto tune = WriteModifiedPhotoVMTune(
        label, [&](auto &j) { j.at("PARAM_TENSORPOM").at("FORWARD_NOFLIP") = noflip; },
        [meson, odderon](auto &j) {
          if (!meson) { j.at("995").at("[333,333]").at("g_tensor") = {0.0, 0.0}; }
          if (!odderon) { j.at("995").at("[333,9993]").at("g_tensor") = {0.0, 0.0}; }
        });
    gra::LORENTZSCALAR  lts = TensorPhiCascadeLTSForTest();
    if (!cascade) {
      for (auto &branch : lts.decaytree) { branch.legs.clear(); }
    }
    const auto model = gra::MModelTune::Load(tune.second);
    gra::MTensorPomeron tensor(lts, model,
                               gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
    lts.hamp.Configure(TensorPairMetadata(noflip));
    const double        amp2 = cascade ? tensor.ME6(lts) : tensor.ME4(lts, gra::TensorContinuumMode::TensorPomeron);
    REQUIRE(std::isfinite(amp2));
    REQUIRE(amp2 > 0.0);
    RequireTensorBornProjection(lts, model->Soft()->GoodWalker().ChannelCount());
    RequireTensorSymmetries(lts, [&](auto &event) {
      return cascade ? tensor.ME6(event) : tensor.ME4(event, gra::TensorContinuumMode::TensorPomeron);
    });
    return lts.hamp;
  };

  const auto meson    = evaluate("tensor_phi_transfer", true, false);
  const auto odderon  = evaluate("tensor_odderon_transfer", false, true);
  const auto coherent = evaluate("tensor_phi_odderon_transfer", true, true);
  REQUIRE(coherent.size() == meson.size());
  REQUIRE(coherent.size() == odderon.size());
  for (const auto &i : indices(coherent)) { RequireComplexNear(coherent[i], meson[i] + odderon[i], 2.0e-10); }
}

// Check that the published phi pair threshold factor multiplies the amplitude
// once
TEST_CASE("Tensor phi-pair Odderon threshold follows Eq. 3.51",
          "[MTensorPomeron][continuum][Odderon][cascade][threshold][physics]") {
  const bool cascade = GENERATE(false, true);
  const auto evaluate = [cascade](const bool threshold, const std::string &label) {
    const auto tune = WriteModifiedPhotoVMTune(
        label, [](auto &j) { j.at("PARAM_TENSORPOM").at("FORWARD_NOFLIP") = true; },
        [threshold](auto &j) {
          j.at("995").at("[333,333]").at("g_tensor")   = {0.0, 0.0};
          j.at("995").at("[333,9993]").at("threshold") = threshold;
        });
    gra::LORENTZSCALAR  lts = TensorPhiCascadeLTSForTest();
    if (!cascade) {
      for (auto &branch : lts.decaytree) { branch.legs.clear(); }
    }
    gra::MTensorPomeron tensor(lts, gra::MModelTune::Load(tune.second),
                               gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
    REQUIRE((cascade ? tensor.ME6(lts) : tensor.ME4(lts, gra::TensorContinuumMode::TensorPomeron)) > 0.0);
    return std::pair{lts, std::vector<std::complex<double>>(lts.hamp.begin(), lts.hamp.end())};
  };

  const auto   unsuppressed = evaluate(false, "tensor_phi_threshold_off");
  const auto   suppressed   = evaluate(true, "tensor_phi_threshold_on");
  const double phi_mass     = unsuppressed.first.PDG.FindByPDG(333).mass;
  const double threshold    = 4.0 * pow2(phi_mass);
  const double factor       = 1.0 - std::exp((threshold - unsuppressed.first.pfinal[0].M2()) / threshold);
  REQUIRE(factor > 0.0);
  REQUIRE(factor < 1.0);
  REQUIRE(suppressed.second.size() == unsuppressed.second.size());
  for (const auto &i : indices(unsuppressed.second)) {
    RequireComplexNear(suppressed.second[i], factor * unsuppressed.second[i], 2.0e-10);
  }
}

// Check both vector legs and the symmetric traceless Pomeron indices off shell
// [REFERENCE: Ewerz et al., arXiv:1309.3478, Eqs. (3.18)-(3.22)]
TEST_CASE("Tensor vector vertices obey Ward and Bose identities", "[MTensorPomeron][vector][physics][gauge]") {
  auto lts = DirectCentralPairLTSForTest(333, 333);
  gra::MTensorPomeron tensor(lts, gra::MModelTune::Load(modelfile),
      gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
  const gra::M4Vec k1(0.23, -0.17, 0.51, 0.81), k2(-0.31, 0.22, -0.42, 0.14);
  for (const int mode : {0, 2}) {
    const auto vertex = mode == 0 ? tensor.Gamma0(k1, k2) : tensor.Gamma2(k1, k2);
    const auto crossed = mode == 0 ? tensor.Gamma0(k2, k1) : tensor.Gamma2(k2, k1);
    for (const auto mu : tensor.LI) {
      for (const auto nu : tensor.LI) {
        std::complex<double> trace = 0.0;
        for (const auto a : tensor.LI) {
          trace += tensor.g[a][a] * vertex(mu, nu, a, a);
          for (const auto b : tensor.LI) {
            RequireComplexNear(vertex(mu, nu, a, b), vertex(mu, nu, b, a), 1.0e-12);
            RequireComplexNear(vertex(mu, nu, a, b), crossed(nu, mu, a, b), 1.0e-12);
          }
          std::complex<double> ward1 = 0.0, ward2 = 0.0;
          for (const auto rho : tensor.LI) {
            ward1 += k1[rho] * vertex(rho, mu, nu, a);
            ward2 += k2[rho] * vertex(mu, rho, nu, a);
          }
          RequireComplexNear(ward1, 0.0, 1.0e-12);
          RequireComplexNear(ward2, 0.0, 1.0e-12);
        }
        RequireComplexNear(trace, 0.0, 1.0e-12);
      }
    }
  }
}

// Check the transverse vector vertex against the pseudoscalar normalization
// [REFERENCE: Ewerz et al., arXiv:1309.3478, Eqs. (7.25)-(7.28)]
TEST_CASE("Tensor vector and pseudoscalar forward residues agree", "[MTensorPomeron][physics][normalization]") {
  auto lts = DirectCentralPairLTSForTest(333, 333);
  const auto model = gra::MModelTune::Load(modelfile);
  gra::MTensorPomeron tensor(lts, model,
      gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
  const auto parameters = gra::ReadTensorPomeronParam(*model, lts.PDG);
  const double beta = parameters->exchange.FindVertex(995, 211).g_tensor.front();
  const double mass = 0.91, a = 0.37, b = 0.61;
  for (const double angle : {0.0, 0.63}) {
    gra::M4Vec k(0.0, 0.0, -0.71, std::sqrt(pow2(0.71) + pow2(mass))), n(0.0, 0.0, 1.0, 1.0);
    k.RotateY(angle);
    n.RotateY(angle);
    const auto eps = tensor.MassiveSpin1States(k, "conj", true);
    const auto vector = tensor.iG_Pvv(k, k, a, b, {});
    const auto scalar = tensor.iG_Tpsps(k, k, 995, 211);
    for (const auto h : {0, 2}) {
      std::complex<double> vector_residue = 0.0, scalar_residue = 0.0;
      for (const auto alpha : tensor.LI) {
        for (const auto beta_index : tensor.LI) {
          scalar_residue += n[alpha] * n[beta_index] * scalar(alpha, beta_index);
          for (const auto mu : tensor.LI) {
            for (const auto nu : tensor.LI) {
              vector_residue += n[alpha] * n[beta_index] * eps[h](mu) * std::conj(eps[h](nu)) *
                                vector(mu, nu, alpha, beta_index);
            }
          }
        }
      }
      RequireComplexNear(vector_residue, scalar_residue * (2.0 * pow2(mass) * a + b) / (4.0 * beta), 1.0e-12);
    }
  }
}

// Accept an Odderon-only vector continuum and reject a fully inactive model during setup
TEST_CASE("Tensor vector continuum initializes every active transfer", "[MTensorPomeron][continuum][process][Odderon]") {
  const bool active = GENERATE(false, true);
  const auto tune = WriteModifiedPhotoVMTune("tensor_vector_active_" + std::to_string(active), {}, [&](auto &j) {
    j.at("995").at("[333,333]").at("g_tensor") = {0.0, 0.0};
    j.at("995").at("[333,9993]").at("g_tensor") = {0.0, active ? 0.71 : 0.0};
  });
  ToyHelicityProcess process;
  ConfigureToyProductionProcess(process, "TP", "CON", "333 333");
  process.SetModelTune(gra::MModelTune::Load(tune.second));
  if (active) {
    REQUIRE_NOTHROW(process.InitializeProcessAmplitude());
  } else {
    REQUIRE_THROWS_AS(process.InitializeProcessAmplitude(), std::invalid_argument);
  }
}

// Reconstruct the complex four-kaon amplitude from physical phi helicities
// [REFERENCE: Lebiedowicz et al., arXiv:1901.11490, Eqs. (3.2)-(3.5)]
TEST_CASE("Tensor vector helicities reconstruct the coherent cascade", "[MTensorPomeron][cascade][physics][phase]") {
  const bool noflip = GENERATE(false, true);
  const auto tune = WriteModifiedPhotoVMTune("tensor_vector_decay_" + std::to_string(noflip), [&](auto &j) {
    j.at("PARAM_TENSORPOM").at("FORWARD_NOFLIP") = noflip;
  });
  auto lts = TensorPhiCascadeLTSForTest();
  gra::MTensorPomeron tensor(lts, gra::MModelTune::Load(tune.second),
      gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
  REQUIRE(tensor.ME6(lts) > 0.0);
  std::vector<std::complex<double>> expected(lts.hamp.size(), 0.0);
  for (const auto &term : gra::spin::StableLeafAmplitudeTrees(lts, lts.decaytree)) {
    auto stable = lts;
    stable.decaytree = term.tree;
    for (auto &branch : stable.decaytree) { branch.legs.clear(); }
    RefreshToyDerivedKinematicsPreserveDecay(stable);
    REQUIRE(tensor.ME4(stable, gra::TensorContinuumMode::TensorPomeron) > 0.0);
    REQUIRE(stable.hamp.size() == expected.size() * 9);
    std::array<std::array<std::complex<double>, 3>, 2> decay{};
    for (const auto i : indices(decay)) {
      const auto &branch = term.tree[i];
      const auto eps = tensor.MassiveSpin1States(branch.p4, "conj", true);
      const auto propagator = tensor.iD_VMES(branch.p4, branch.p.mass, branch.p.width, branch.p.pdg, true, true);
      const auto vertex = tensor.iG_vpsps(branch.legs[0].p4, branch.legs[1].p4, branch.p.mass,
                                          branch.hel.g_decay_TP[0], branch.hel.ff_decay);
      for (const auto h : indices(eps)) {
        for (const auto mu : tensor.LI) {
          for (const auto nu : tensor.LI) {
            decay[i][h] -= std::conj(eps[h](mu)) * tensor.g[mu][mu] * propagator(mu, nu) * vertex(nu);
          }
        }
      }
    }
    for (const auto row : indices(expected)) {
      for (const auto h3 : indices(decay[0])) {
        for (const auto h4 : indices(decay[1])) {
          expected[row] += term.statistics_sign * stable.hamp[row * 9 + h3 * 3 + h4] * decay[0][h3] * decay[1][h4];
        }
      }
    }
  }
  RequireVectorNear(lts.hamp, expected, 3.0e-10);
}
