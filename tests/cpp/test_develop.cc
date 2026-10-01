// Check generated GP coupling cards against the C++ pole and helicity API
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include <array>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>

#include "Graniitti/Regge/MReggeInit.h"
#include "Graniitti/Regge/MReggeGPInit.h"
#include "support/models_test_support.hh"

namespace {

// Read coupling JSON from the repository tools in the normal test environment
nlohmann::json ToolCouplings(const std::string& tool, const std::string& arguments) {
  return nlohmann::json::parse(gra::aux::ExecCommand("python develop/tools/" + tool + ".py " + arguments));
}

// Read a generated channel through the normal resonance card parser
gra::RES_PRODUCTION_CHANNEL CouplingChannel(const nlohmann::json& block, const std::array<int, 2>& exchange,
                                            const gra::MParticle& mother) {
  auto  card = nlohmann::json::parse(gra::aux::GetInputData(gra::ResolveModelDataFile("TUNE0", "RES/f0_500.json")));
  auto& resonance           = card.at("PARAM_RES");
  resonance["spinX2"]       = mother.spinX2;
  resonance["P"]            = mother.P;
  resonance["C"]            = mother.C;
  const std::string key     = "[" + std::to_string(exchange[0]) + "," + std::to_string(exchange[1]) + "]";
  auto &gp = resonance["MODELS"]["GP"];
  gp = {{"mass", gp.at("mass")}, {"width", gp.at("width")}, {"BW", gp.at("BW")},
        {"phi", 0.0}, {"FF_transfer", {{"type", "none"}}}, {"FF_prod", {{"type", "none"}}}, {key, block}};
  for (auto& [name, channel] : resonance["MODELS"]["MP"].items()) {
    if (!name.empty() && name.front() == '[') {
      channel["basis"] = "auto_min_L";
      channel["polarization"] = {{"mode", "none"}};
    }
  }
  for (auto& [name, channel] : resonance["MODELS"]["TP"].items()) {
    if (!name.empty() && name.front() == '[') {
      std::vector<double> tensor(mother.spinX2 == 4 ? 7 : 2, 0.0);
      tensor.front()      = 1.0;
      channel["g_tensor"] = tensor;
    }
  }
  std::filesystem::create_directories("tmp");
  std::string directory = "tmp/develop_cards_XXXXXX";
  REQUIRE(mkdtemp(directory.data()) != nullptr);
  std::filesystem::copy_file(gra::ResolveModelDataFile("TUNE0", "PDG_EXTRA.json"), directory + "/PDG_EXTRA.json");
  std::ofstream output(directory + "/RES.json");
  output << card.dump(2);
  output.close();
  REQUIRE(output.good());
  gra::MRandom rng;
  const auto   parsed = gra::resonance::Read("RES.json", rng, gra::ReggeProductionModel::GP, directory);
  REQUIRE(parsed.GP.channels.size() == 1);
  return parsed.GP.channels.front();
}

}  // namespace

// Check DL pole strengths with the real C++ continuum readers and complex tensors
TEST_CASE("DL continuum residues satisfy the stated spin averaged pole strength", "[develop][DL][normalization]") {
  for (const std::string shape : {"preserve", "lowest"}) {
    const auto report = ToolCouplings("DL_couplings", "--format json --shape " + shape);
    std::filesystem::create_directories("tmp");
    std::string directory = "tmp/dl_poles_XXXXXX";
    REQUIRE(mkdtemp(directory.data()) != nullptr);
    std::filesystem::copy(std::filesystem::path(modelfile).parent_path(), directory,
                          std::filesystem::copy_options::recursive | std::filesystem::copy_options::overwrite_existing);
    for (const std::string model : {"MP", "XP", "GP"}) {
      std::ofstream(directory + "/CON_" + model + ".json") << report.at("CON").at(model).dump(2);
    }
    auto general = nlohmann::json::parse(gra::aux::GetInputData(directory + "/GENERAL.json"));
    for (const std::string model : {"MP", "XP", "GP"}) {
      general["PARAM_REGGE"]["PARAM_CON"][model]["[311,-311]"] = model == "GP"
          ? nlohmann::json{{990, 990}, {9910, 9910}, {9930, 9930}}
          : nlohmann::json{{995, 995}, {9915, 9915}, {9933, 9933}};
    }
    std::ofstream(directory + "/GENERAL.json") << general.dump(2);
    for (const auto& final : std::array<std::pair<int, std::string>, 5>{{
             {211, "pi+ pi-"}, {311, "K0 K0~"}, {2212, "p+ p-"}, {113, "rho(770)0 rho(770)0"}, {333, "phi(1020)0 phi(1020)0"}}}) {
      std::map<std::string, gra::HELMatrix> reference;
      for (const std::string model : {"MP", "XP", "GP"}) {
        ToyHelicityProcess process;
        ConfigureToyProductionProcess(process, model, "CON", final.second);
        process.SetTuneForTest(directory);
        REQUIRE_NOTHROW(process.InitializeProcessAmplitude());
        const auto& configured = process.state.lts.process;
        for (const auto& channel : indices(configured.CONT_PRODUCTION)) {
          const auto& pair = configured.CONT_PRODUCTION[channel];
          if (pair[0] != pair[1]) { continue; }
          const auto scalar = std::find_if(report.at("normalization_checks").begin(), report.at("normalization_checks").end(),
              [&](const auto& row) {
                return row.at("model") == model && row.at("exchange") == std::to_string(pair[0]) &&
                       row.at("pair")[0] == final.first && row.at("sector") == (final.first == 113 || final.first == 333 ? "self" : "opposite");
              });
          REQUIRE(scalar != report.at("normalization_checks").end());
          const double residue = scalar->at("scalar_residue").template get<double>();
          gra::HELMatrix pole;
          if (model == "GP") {
            const auto& input = configured.CONTINUUM_GP[channel][0];
            pole = gra::gpom::Crossed(input, input.J, gra::M4Vec{}, false).helicity;
          } else {
            pole = gra::spin::PoleLSHelicity(configured.CONTINUUM_POLE[channel][0].Pole(), 1.0);
          }
          const double norm2 = model == "GP" ? gra::SquaredNorm(pole.T.Column(gra::gpom::AnalyticMIndex(0, pole.analytic_MMAX, "DL pole")))
                                             : gra::spin::ReducedHelicityNorm2(pole);
          CAPTURE(shape, model, final.first, pair[0], norm2, residue);
          CHECK(norm2 / (pole.s1 * 2.0 + 1.0) == Approx(residue * residue).epsilon(2e-8));
          const auto& exchanges = general.at("PARAM_REGGE").at("EXCHANGES");
          const auto mapping = std::find_if(exchanges.begin(), exchanges.end(), [&](const auto& row) {
            const auto& aliases = row.at("pdg");
            return std::find(aliases.begin(), aliases.end(), pair[0]) != aliases.end();
          });
          REQUIRE(mapping != exchanges.end());
          const std::string exchange = mapping->at("soft_exchange");
          if (model == "MP") { reference[exchange] = pole; }
          if (shape != "lowest" || model == "MP") { continue; }
          const auto& expected = reference.at(exchange);
          for (std::size_t row = 0; row < pole.lambda_values.size_row(); ++row) {
            const auto i = pole.lambda_idx[row][0], j = pole.lambda_idx[row][1];
            const auto value = model == "GP" ? pole.T[row][gra::gpom::AnalyticMIndex(0, pole.analytic_MMAX, "DL pole")] : pole.T[i][j];
            CHECK(std::abs(value - expected.T[i][j]) <= 2e-8 * std::max(1.0, residue));
          }
        }
      }
    }
  }
}

// Check complete pole tables and the normalization of the selected active coupling
TEST_CASE("GP card builder supplies complete physical pole LS tables", "[develop][GP][LS][normalization]") {
  auto        lts = DirectCentralPairLTSForTest(211, -211);
  gra::MRegge regge(lts, gra::MModelTune::Load(modelfile),
                    gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));
  const auto  param = ReggeParametersForTest(regge, lts);
  struct Channel {
    int         first;
    int         second;
    int         spin2;
    int         parity;
    int         cparity;
    std::size_t s;
  };
  const std::array<Channel, 4> channels = {
      {{990, 990, 0, 1, 1, 0}, {22, 990, 2, -1, -1, 1}, {990, 9990, 2, -1, -1, 1}, {22, 22, 0, 1, 1, 0}}};
  for (const auto& input : channels) {
    const double magnitude = input.first == 22 && input.second == 22 ? 1.0 : 0.42;
    const auto   output    = ToolCouplings(
        "resonance_card_builder", "--model GP --format card --basis g_ls --fuse " + std::to_string(input.first) +
                                      " " + std::to_string(input.second) + " --spinX2 " + std::to_string(input.spin2) +
                                      " --P " + std::to_string(input.parity) + " --C " + std::to_string(input.cparity) +
                                      " --active-ls 0 " + std::to_string(input.s) + " --channel-mag " +
                                      std::to_string(magnitude) + " --channel-phase 0.3");
    const std::string key = "[" + std::to_string(input.first) + "," + std::to_string(input.second) + "]";
    gra::MParticle    mother;
    mother.spinX2       = input.spin2;
    mother.P            = input.parity;
    mother.C            = input.cparity;
    const auto  channel = CouplingChannel(output.at("GP").at(key), {input.first, input.second}, mother);
    const auto& first   = lts.PDG.FindByPDG(input.first);
    const auto& second  = lts.PDG.FindByPDG(input.second);
    const auto& pole1   = input.first == 22 ? first : gra::regge::PoleRepresentative(*param, lts.PDG, input.first);
    const auto& pole2   = input.second == 22 ? second : gra::regge::PoleRepresentative(*param, lts.PDG, input.second);
    const auto  allowed = gra::spin::CanonicalPoleOperators(mother, pole1, pole2, true, true, true);
    REQUIRE(channel.g_ls.Size() == allowed.size());
    const int  mmax   = std::max(pole1.spinX2, pole2.spinX2) / 2;
    const auto vertex = gra::gpom::PrepareResonance(mother, {first, second}, channel, *param, lts.PDG, mmax, 0.0);
    REQUIRE(vertex.alpha_ls.Size() == 1);
    RequireComplexNear(vertex.alpha_ls.At(0, 2 * input.s), std::polar(magnitude, 0.3), 2.0e-12);
  }
}

// Check emitted independent orbits reconstruct the analytic Jz sectors in C++
TEST_CASE("Analytic Jz cards reproduce GP parity and exchange sectors", "[develop][GP][helicity][parity]") {
  auto              lts = DirectCentralPairLTSForTest(211, -211);
  gra::MRegge       regge(lts, gra::MModelTune::Load(modelfile),
                          gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));
  const auto        param     = ReggeParametersForTest(regge, lts);
  const auto&       leg       = lts.PDG.FindByPDG(990);
  const int         pole_spin = gra::regge::PoleRepresentative(*param, lts.PDG, 990).spinX2 / 2;
  const int         mmax      = std::min(2, pole_spin);
  const std::string arguments = "--J 0,1,2 --parity 1 --pole-spin " + std::to_string(pole_spin) +
                                " --mmax " + std::to_string(mmax);
  const auto        dense     = ToolCouplings("analytic_Jz", arguments + " --format json");
  const auto        cards     = ToolCouplings("analytic_Jz", arguments + " --format cards");
  for (const auto& model : dense.at("models")) {
    gra::MParticle mother;
    mother.spinX2 = 2 * model.at("J").get<int>();
    mother.P      = model.at("P").get<int>();
    mother.C      = 1;
    for (const auto& sector : model.at("sectors")) {
      if (!sector.at("found").get<bool>()) { continue; }
      const std::string key = "abs_Jz_" + std::to_string(sector.at("abs_Jz").get<int>());
      const auto channel    = CouplingChannel(cards.at(model.at("J_P").get<std::string>()).at(key), {990, 990}, mother);
      const auto vertex     = gra::gpom::PrepareResonance(mother, {leg, leg}, channel, *param, lts.PDG, mmax, 0.0);
      REQUIRE(vertex.T.FrobNorm2() == Approx(1.0).epsilon(2.0e-12));
      for (const auto& row : sector.at("helicity")) {
        const std::size_t i = gra::gpom::AnalyticMIndex(row[0].get<int>(), mmax, "Jz upper");
        const std::size_t j = gra::gpom::AnalyticMIndex(row[1].get<int>(), mmax, "Jz lower");
        RequireComplexNear(vertex.T[i][j], std::polar(row[2].get<double>(), row[3].get<double>()), 2.0e-12);
      }
    }
  }
}
