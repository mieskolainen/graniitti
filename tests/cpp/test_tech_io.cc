// Regression tests for steering input and Les Houches event conversion
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "catch.hpp"

#include <cmath>
#include <array>
#include <filesystem>
#include <fstream>
#include <limits>
#include <memory>
#include <string>

#include "Graniitti/Particle/MPDG.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MJsonOverride.h"
#include "Graniitti/Tech/MLHE.h"
#include "HepMC3/Attribute.h"
#include "HepMC3/GenCrossSection.h"
#include "HepMC3/GenPdfInfo.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/WriterAscii.h"
#include "json.hpp"

using gra::aux::indices;

namespace {

// Keep test outputs under the project temporary directory
std::string OutputPath(const std::string &name) {
  const auto dir = std::filesystem::path(gra::aux::ResolveProjectPath("tmp")) / "test_tech_io";
  std::filesystem::create_directories(dir);
  return (dir / name).string();
}

// Construct a momentum-conserving proton collision with a lepton pair
HepMC3::GenEvent LeptonEvent() {
  HepMC3::GenEvent event(HepMC3::Units::GEV, HepMC3::Units::MM);
  auto vertex = std::make_shared<HepMC3::GenVertex>();
  const double beam_pz = std::sqrt(100.0 - 0.9383 * 0.9383);
  const double muon_pz = std::sqrt(100.0 - 9.0 - 0.1057 * 0.1057);
  vertex->add_particle_in(std::make_shared<HepMC3::GenParticle>(HepMC3::FourVector(0, 0, beam_pz, 10), 2212, 4));
  vertex->add_particle_in(std::make_shared<HepMC3::GenParticle>(HepMC3::FourVector(0, 0, -beam_pz, 10), 2212, 4));
  vertex->add_particle_out(std::make_shared<HepMC3::GenParticle>(HepMC3::FourVector(3, 0, muon_pz, 10), 13, 1));
  vertex->add_particle_out(std::make_shared<HepMC3::GenParticle>(HepMC3::FourVector(-3, 0, -muon_pz, 10), -13, 1));
  event.add_vertex(vertex);
  auto xs = std::make_shared<HepMC3::GenCrossSection>();
  xs->set_cross_section(2.0, 0.2);
  event.set_cross_section(xs);
  event.weights() = {1.0};
  return event;
}

// Construct two color singlets in a momentum-conserving photon collision
HepMC3::GenEvent FourQuarkEvent() {
  HepMC3::GenEvent event(HepMC3::Units::GEV, HepMC3::Units::MM);
  auto vertex = std::make_shared<HepMC3::GenVertex>();
  for (const double sign : {-1.0, 1.0}) {
    vertex->add_particle_in(std::make_shared<HepMC3::GenParticle>(HepMC3::FourVector(0, 0, sign * 10, 10), 22, 4));
  }
  const std::array<int, 4> pdg{2, -2, 1, -1};
  const std::array<HepMC3::FourVector, 4> momenta{{{3, 0, 4, 5}, {-3, 0, -4, 5},
                                                  {0, 3, 4, 5}, {0, -3, -4, 5}}};
  for (const auto i : indices(pdg)) {
    vertex->add_particle_out(std::make_shared<HepMC3::GenParticle>(momenta[i], pdg[i], 1));
  }
  event.add_vertex(vertex);
  for (const auto i : indices(pdg)) {
    vertex->particles_out()[i]->add_attribute(pdg[i] > 0 ? "flow1" : "flow2",
        std::make_shared<HepMC3::IntAttribute>(701 + i / 2));
  }
  auto pdf = std::make_shared<HepMC3::GenPdfInfo>();
  pdf->set(22, 22, 1.0, 1.0, 10.0, 1.0, 1.0);
  event.set_pdf_info(pdf);
  return event;
}

// Construct spacelike photon fusion into a named resonance with exact proton recoils
HepMC3::GenEvent FusionEvent() {
  HepMC3::GenEvent event(HepMC3::Units::GEV, HepMC3::Units::MM);
  auto fusion = std::make_shared<HepMC3::GenVertex>();
  const double mass = 0.9383;
  const double beam_pz = std::sqrt(100.0 - mass * mass);
  const double recoil_pz = std::sqrt(36.0 - mass * mass - 1.0);
  for (const double sign : {1.0, -1.0}) {
    auto vertex = std::make_shared<HepMC3::GenVertex>();
    vertex->add_particle_in(std::make_shared<HepMC3::GenParticle>(HepMC3::FourVector(0, 0, sign * beam_pz, 10), 2212, 4));
    vertex->add_particle_out(std::make_shared<HepMC3::GenParticle>(HepMC3::FourVector(-sign, 0, sign * recoil_pz, 6), 2212, 1));
    auto photon = std::make_shared<HepMC3::GenParticle>(HepMC3::FourVector(sign, 0, sign * (beam_pz - recoil_pz), 4),
                                                      gra::PDG::PDG_propagator, gra::PDG::PDG_INTERMEDIATE);
    vertex->add_particle_out(photon);
    fusion->add_particle_in(photon);
    event.add_vertex(vertex);
  }
  auto resonance = std::make_shared<HepMC3::GenParticle>(HepMC3::FourVector(0, 0, 0, 8), 23, 2);
  fusion->add_particle_out(resonance);
  event.add_vertex(fusion);
  auto decay = std::make_shared<HepMC3::GenVertex>();
  decay->add_particle_in(resonance);
  const double momentum = std::sqrt(16.0 - mass * mass);
  for (const int sign : {1, -1}) {
    decay->add_particle_out(std::make_shared<HepMC3::GenParticle>(HepMC3::FourVector(sign * momentum, 0, 0, 4), sign * 2212, 1));
  }
  event.add_vertex(decay);
  auto pdf = std::make_shared<HepMC3::GenPdfInfo>();
  pdf->set(22, 22, 0.4, 0.4, 8.0, 1.0, 1.0);
  event.set_pdf_info(pdf);
  return event;
}

}  // namespace

// Check lexical comment handling without changing quoted input
TEST_CASE("Steering comments preserve JSON strings and token boundaries", "[MAux][input]") {
  const auto path = OutputPath("comments.json");
  {
    std::ofstream out(path);
    out << R"json({
      /* a comment
         spanning lines */
      "url": "https://example.org/a//b",
      "label": "two  spaces /* literal */ ,] ,}",
      "escape": "quote\" // literal and backslash\\",
      "values": [1, 2, /* trailing comment */],
    })json" << "// final comment without newline";
  }
  const auto card = nlohmann::json::parse(gra::aux::GetInputData(path));
  CHECK(card.at("url") == "https://example.org/a//b");
  CHECK(card.at("label") == "two  spaces /* literal */ ,] ,}");
  CHECK(card.at("escape") == "quote\" // literal and backslash\\");
  CHECK(card.at("values") == nlohmann::json::array({1, 2}));

  SECTION("comments do not concatenate numeric tokens") {
    { std::ofstream out(path); out << "{\"value\":1/* separator */2}"; }
    CHECK_THROWS([&] { return nlohmann::json::parse(gra::aux::GetInputData(path)); }());
  }
  SECTION("unterminated comments are fatal input errors") {
    { std::ofstream out(path); out << "{\"value\":1} /* unclosed"; }
    CHECK_THROWS_AS(gra::aux::GetInputData(path), std::invalid_argument);
  }
}

// Check data tables can be read without write access and missing inputs are fatal
TEST_CASE("CSV input reads data without write access", "[MAux][input]") {
  const auto path = OutputPath("readonly.csv");
  { std::ofstream out(path); out << "1.0,2.0\n3.0,4.0\n"; }
  std::filesystem::permissions(path, std::filesystem::perms::owner_read);
  std::vector<std::vector<std::string>> rows;
  CHECK_NOTHROW(gra::aux::ReadCSV(path, rows));
  std::filesystem::permissions(path, std::filesystem::perms::owner_read | std::filesystem::perms::owner_write);
  CHECK(rows == std::vector<std::vector<std::string>>{{"1.0", "2.0"}, {"3.0", "4.0"}});
  CHECK_THROWS_AS(gra::aux::ReadCSV(OutputPath("missing/file.csv"), rows), std::invalid_argument);
}

// Check physical four-momenta and unspecified LHA metadata through the writer
TEST_CASE("LHE converts MeV momenta and preserves GeV PDF scales", "[MLHE][units]") {
  auto event = LeptonEvent();
  const auto reference = gra::BuildLHERows(event);
  event.set_units(HepMC3::Units::MEV, HepMC3::Units::MM);
  const auto converted = gra::BuildLHERows(event);
  REQUIRE(converted.size() == reference.size());
  for (const auto i : indices(converted)) {
    CHECK(converted[i].pid == reference[i].pid);
    for (const auto j : indices(converted[i].momentum)) {
      CHECK(converted[i].momentum[j] == Approx(reference[i].momentum[j]));
    }
  }
  CHECK(event.momentum_unit() == HepMC3::Units::MEV);

  const auto path = OutputPath("units.lhe");
  SECTION("unknown scale and coupling use LHA sentinel values") {
    gra::MLHEWriter writer(path);
    writer.WriteEvent(event);
    writer.Close();
    CHECK_NOTHROW(writer.Close());
    CHECK_THROWS_AS(writer.WriteEvent(event), std::runtime_error);
    LHEF::Reader reader(path);
    CHECK(reader.heprup.EBMUP.first == Approx(10.0));
    CHECK(reader.heprup.EBMUP.second == Approx(10.0));
    CHECK(reader.heprup.PDFGUP.first == -1);
    CHECK(reader.heprup.PDFSUP.second == -1);
    REQUIRE(reader.readEvent());
    CHECK(reader.hepeup.SCALUP == Approx(-1.0));
    CHECK(reader.hepeup.AQEDUP == Approx(-1.0));
    CHECK(reader.hepeup.AQCDUP == Approx(-1.0));
    for (const double spin : reader.hepeup.SPINUP) { CHECK(spin == Approx(9.0)); }
    for (const auto i : indices(reference)) {
      for (const auto j : indices(reference[i].momentum)) {
        CHECK(reader.hepeup.PUP[i][j] == Approx(reference[i].momentum[j]));
      }
    }
  }
  SECTION("PDF scale remains in GeV for a MeV event") {
    auto pdf = std::make_shared<HepMC3::GenPdfInfo>();
    pdf->set(22, 22, 0.1, 0.1, 5.0, 1.0, 1.0);
    event.set_pdf_info(pdf);
    event.add_attribute("graniitti_alpha_qed", std::make_shared<HepMC3::DoubleAttribute>(0.0073));
    event.add_attribute("graniitti_alpha_qcd", std::make_shared<HepMC3::DoubleAttribute>(0.12));
    gra::MLHEWriter writer(path);
    writer.WriteEvent(event);
    writer.Close();
    LHEF::Reader reader(path);
    REQUIRE(reader.readEvent());
    CHECK(reader.hepeup.SCALUP == Approx(5.0));
    CHECK(reader.hepeup.scales.muf == Approx(5.0));
    CHECK(reader.hepeup.AQEDUP == Approx(0.0073));
    CHECK(reader.hepeup.AQCDUP == Approx(0.12));
  }
}

// Check normalization against the final running cross section for the whole file
TEST_CASE("HepMC conversion uses the final run cross section", "[MLHE][normalization]") {
  auto event = LeptonEvent();
  const bool weighted = GENERATE(false, true);
  const auto input = OutputPath(weighted ? "weighted.hepmc3" : "unit.hepmc3");
  const auto output = OutputPath(weighted ? "weighted.lhe" : "unit.lhe");
  {
    HepMC3::WriterAscii writer(input);
    event.cross_section()->set_cross_section(1.0, 0.1);
    writer.write_event(event);
    event.set_event_number(1);
    event.cross_section()->set_cross_section(2.0, 0.2);
    event.weights() = {weighted ? 3.0 : 1.0};
    writer.write_event(event);
    writer.close();
    REQUIRE_FALSE(writer.failed());
  }
  const auto stats = gra::ConvertHepMC3ToLHE(input, output, false);
  CHECK(stats.events == 2);
  LHEF::Reader reader(output);
  CHECK(reader.heprup.XSECUP.at(0) == Approx(2.0));
  CHECK(reader.heprup.XERRUP.at(0) == Approx(0.2));
  CHECK(reader.heprup.IDWTUP == (weighted ? 4 : 3));
  REQUIRE(reader.readEvent());
  CHECK(reader.hepeup.XWGTUP == Approx(1.0));
  REQUIRE(reader.readEvent());
  CHECK(reader.hepeup.XWGTUP == Approx(weighted ? 3.0 : 1.0));
  CHECK_FALSE(reader.readEvent());
}

// Check fatal output errors and fixed run validation through the public API
TEST_CASE("LHE reports output failures and incompatible streaming runs", "[MLHE][io]") {
  const auto directory = OutputPath("directory");
  std::filesystem::create_directories(directory);
  CHECK_THROWS_AS(gra::MLHEWriter(directory), std::runtime_error);
  auto event = LeptonEvent();
  gra::MLHEWriter writer(OutputPath("fixed_run.lhe"));
  writer.WriteEvent(event);
  event.cross_section()->set_cross_section(3.0, 0.3);
  CHECK_THROWS_AS(writer.WriteEvent(event), std::invalid_argument);
  writer.Close();
}

// Check a comma without a preceding value remains a fatal syntax error
TEST_CASE("Steering trailing commas cannot create missing values", "[MAux][input]") {
  const auto path = OutputPath("empty_comma.json");
  for (const auto &input : {"{,}", "[, ]", "[/* comment */,]", "[1,,]"}) {
    { std::ofstream out(path); out << input; }
    CHECK_THROWS([&] { return nlohmann::json::parse(gra::aux::GetInputData(path)); }());
  }
}

// Check conversion never destroys its input through an identical path or hard link
TEST_CASE("LHE conversion rejects input and output file aliases", "[MLHE][io]") {
  const auto input = OutputPath("same_file.hepmc3");
  {
    HepMC3::WriterAscii writer(input);
    writer.write_event(LeptonEvent());
    writer.close();
  }
  std::ifstream before(input);
  const std::string original((std::istreambuf_iterator<char>(before)), std::istreambuf_iterator<char>());
  CHECK_THROWS_AS(gra::ConvertHepMC3ToLHE(input, input, false), std::invalid_argument);
  const auto alias = OutputPath("same_file_alias.hepmc3");
  if (!std::filesystem::exists(alias)) { std::filesystem::create_hard_link(input, alias); }
  CHECK_THROWS_AS(gra::ConvertHepMC3ToLHE(input, alias, false), std::invalid_argument);
  std::ifstream after(input);
  const std::string preserved((std::istreambuf_iterator<char>(after)), std::istreambuf_iterator<char>());
  CHECK(preserved == original);
}

// Check a corrupt event after valid events is not converted as a successful prefix
TEST_CASE("LHE conversion rejects malformed HepMC input", "[MLHE][io]") {
  const auto input = OutputPath("malformed.hepmc3");
  const auto output = OutputPath("malformed.lhe");
  {
    auto event = LeptonEvent();
    HepMC3::WriterAscii writer(input);
    writer.write_event(event);
    event.set_event_number(1);
    writer.write_event(event);
    writer.close();
  }
  std::string contents;
  {
    std::ifstream stream(input);
    contents.assign(std::istreambuf_iterator<char>(stream), std::istreambuf_iterator<char>());
  }
  const auto start = contents.find("\nE 1 ");
  REQUIRE(start != std::string::npos);
  const auto end = contents.find('\n', start + 1);
  contents.replace(start + 1, end - start - 1, "E 1 1 9");
  { std::ofstream stream(input); stream << contents; }
  { std::ofstream stream(output); stream << "existing output"; }
  CHECK_THROWS_AS(gra::ConvertHepMC3ToLHE(input, output, false), std::invalid_argument);
  std::ifstream stream(output);
  std::string preserved;
  std::getline(stream, preserved);
  CHECK(preserved == "existing output");
}

// Check optional per-event cross-section attributes do not erase the latest estimate
TEST_CASE("LHE conversion retains the last available cross section", "[MLHE][normalization]") {
  const auto input = OutputPath("optional_xs.hepmc3");
  const auto output = OutputPath("optional_xs.lhe");
  {
    auto event = LeptonEvent();
    HepMC3::WriterAscii writer(input);
    writer.write_event(event);
    event.set_event_number(1);
    event.remove_attribute("GenCrossSection");
    writer.write_event(event);
    writer.close();
  }
  CHECK(gra::ConvertHepMC3ToLHE(input, output, false).events == 2);
  LHEF::Reader reader(output);
  CHECK(reader.heprup.XSECUP.at(0) == Approx(2.0));
  CHECK(reader.heprup.XERRUP.at(0) == Approx(0.2));
  REQUIRE(reader.readEvent());
  REQUIRE(reader.readEvent());
  CHECK_FALSE(reader.readEvent());
}

// Check numeric steering consumes a complete integer and never truncates a PDG id
TEST_CASE("Integer input rejects incomplete and trailing numeric text", "[MAux][MPDG][input]") {
  for (const auto input : {"", "+", "-", "--211", "211-1", "2.5", "2e3", "211junk"}) {
    CHECK_FALSE(gra::aux::IsIntegerDigits(input));
    CHECK_THROWS_AS(gra::aux::SplitStr2Int(input), std::invalid_argument);
  }
  CHECK(gra::aux::IsIntegerDigits("-211"));
  CHECK(gra::aux::IsIntegerDigits("+211"));
  CHECK(gra::aux::SplitStr2Int("-211, +211, 22") == std::vector<int>{-211, 211, 22});
  gra::MPDG pdg;
  pdg.ReadParticleData();
  CHECK(pdg.FindByPDGName("211").pdg == 211);
  CHECK_THROWS_AS(pdg.FindByPDGName("211-1"), std::invalid_argument);
}

// Check malformed command brackets and trailing text remain fatal input errors
TEST_CASE("Process commands require complete ordered syntax", "[MAux][input]") {
  for (const auto input : {"@PDG]211[{M:0.14}", "@PDG[[211]]{M:0.14}",
                           "@PDG[211]}M:0.14{", "@PDG[211]{M:0.14}ignored",
                           "@PDG[211]{M:{0.14}}", "@PDG{M:[211]0.14}",
                           "@PDG{M:0.14}[211]", "@", "@{M:0.14}"}) {
    CHECK_THROWS_AS(gra::aux::SplitCommands(input), std::invalid_argument);
  }
  const auto commands = gra::aux::SplitCommands("@PDG[211]{M:0.14,W:0.0} @FLATAMP:1");
  REQUIRE(commands.size() == 2);
  CHECK(commands[0].id == "PDG");
  CHECK(commands[0].target == std::vector<std::string>{"211"});
  CHECK(commands[0].arg.at("M") == "0.14");
}

// Check invalid momenta and parton fractions cannot enter an LHE particle record
TEST_CASE("LHE rejects non-finite momenta and invalid parton fractions", "[MLHE][validation]") {
  for (const double value : {std::numeric_limits<double>::quiet_NaN(),
                              std::numeric_limits<double>::infinity()}) {
    auto event = LeptonEvent();
    event.particles().back()->set_momentum(HepMC3::FourVector(value, 0.0, 0.0, 1.0));
    CHECK_THROWS_AS(gra::BuildLHERows(event), std::invalid_argument);
  }
  for (const double x : {0.0, -0.1, 1.1, std::numeric_limits<double>::quiet_NaN(),
                          std::numeric_limits<double>::infinity()}) {
    auto event = LeptonEvent();
    auto pdf = std::make_shared<HepMC3::GenPdfInfo>();
    pdf->set(22, 22, x, 0.1, 5.0, 1.0, 1.0);
    event.set_pdf_info(pdf);
    CHECK_THROWS_AS(gra::BuildLHERows(event), std::invalid_argument);
  }
}

// Check a rejected first event cannot fix the output run cross section
TEST_CASE("LHE run identity is set by the first accepted event", "[MLHE][io]") {
  const auto path = OutputPath("first_accepted.lhe");
  auto event = LeptonEvent();
  gra::MLHEWriter writer(path);
  event.weights() = {2.0};
  CHECK_THROWS_AS(writer.WriteEvent(event), std::invalid_argument);
  event.weights() = {1.0};
  event.cross_section()->set_cross_section(3.0, 0.3);
  REQUIRE_NOTHROW(writer.WriteEvent(event));
  writer.Close();
  LHEF::Reader reader(path);
  CHECK(reader.heprup.XSECUP.at(0) == Approx(3.0));
  REQUIRE(reader.readEvent());
  CHECK_FALSE(reader.readEvent());
}

// Check shared numeric readers preserve whitespace and reject suffixes and non-finite values
TEST_CASE("Numerical steering requires complete finite values", "[MAux][input]") {
  CHECK(gra::aux::ParseInt(" \t-211\r\n", "test") == -211);
  CHECK(gra::aux::ParseDouble(" \t1.4e-1\r\n", "test") == Approx(0.14));
  CHECK(gra::aux::SplitStr2Int("\t211, -211\n") == std::vector<int>{211, -211});
  for (const auto value : {"1junk", "1.5", "1e2", "999999999999999999999"}) {
    CHECK_THROWS_AS(gra::aux::ParseInt(value, "test"), std::invalid_argument);
  }
  for (const auto value : {"0.14GeV", "1 2", "1e", "nan", "inf", "-inf", "1e999"}) {
    CHECK_THROWS_AS(gra::aux::ParseDouble(value, "test"), std::invalid_argument);
    CHECK_THROWS_AS(gra::aux::SplitStr(value, 0.0), std::invalid_argument);
  }
  CHECK(gra::aux::SplitStr("-211, 211", 0) == std::vector<int>{-211, 211});
  CHECK(gra::aux::SplitStr("0.25, 0.5", 0.0) == std::vector<double>{0.25, 0.5});
  CHECK_THROWS_AS(gra::aux::SplitStr("211junk", 0), std::invalid_argument);
}

// Check escaped steering keys change the requested parameter without adding another key
TEST_CASE("JSON overrides decode quoted keys and reject missing separators", "[json-override][input]") {
  nlohmann::json card = {{"a\"b", 1}, {"c'd", 1}, {"e\\f", 1}};
  const auto specs = gra::json_override::ParseSpecs({R"(["a\"b"]=2)", R"(['c\'d']=3)", R"(["e\\f"]=4)"});
  gra::json_override::ApplyInputOverrides(card, specs);
  CHECK(card.size() == 3);
  CHECK(card.at("a\"b").get<int>() == 2);
  CHECK(card.at("c'd").get<int>() == 3);
  CHECK(card.at("e\\f").get<int>() == 4);
  CHECK_THROWS_AS(gra::json_override::ParseSpec("A[0]B=1"), std::invalid_argument);
  CHECK_THROWS_AS(gra::json_override::ApplyOverride(card, {}), std::invalid_argument);
}

// Check supplied color connections survive conversion without changing singlet pairings
TEST_CASE("LHE preserves supplied hard color flow and rejects ambiguous pairings", "[MLHE][color]") {
  auto event = FourQuarkEvent();
  const auto rows = gra::BuildLHERows(event);
  REQUIRE(rows.size() == 6);
  for (const auto &row : rows) {
    if (row.status != 1) { continue; }
    REQUIRE(row.particle != nullptr);
    const auto flow = row.particle->attribute<HepMC3::IntAttribute>(row.pid > 0 ? "flow1" : "flow2");
    REQUIRE(flow != nullptr);
    CHECK((row.pid > 0 ? row.colors.first : row.colors.second) == flow->value());
  }
  for (auto &particle : event.particles()) {
    particle->remove_attribute("flow1");
    particle->remove_attribute("flow2");
  }
  CHECK_THROWS_AS(gra::BuildLHERows(event), std::invalid_argument);
}

// Check the MG5 sextet convention survives the real HepMC to LHE conversion
TEST_CASE("LHE preserves signed sextet color tags", "[MLHE][color][MG2GRA]") {
  auto event = FourQuarkEvent();
  const std::array<int, 4> pdgs{6000001, -6000001, 22, 22};
  const std::array<std::pair<int, int>, 4> colors{{{701, -702}, {-701, 702}, {0, 0}, {0, 0}}};
  const auto outgoing = event.vertices().front()->particles_out();
  for (const auto i : indices(pdgs)) {
    outgoing[i]->set_pid(pdgs[i]);
    outgoing[i]->add_attribute("flow1", std::make_shared<HepMC3::IntAttribute>(colors[i].first));
    outgoing[i]->add_attribute("flow2", std::make_shared<HepMC3::IntAttribute>(colors[i].second));
  }
  const auto rows = gra::BuildLHERows(event);
  REQUIRE(rows.size() == 6);
  for (const auto i : indices(colors)) { CHECK(rows[i + 2].colors == colors[i]); }
}

// Check virtualities, central baryons and four-momentum conservation in the hard LHE view
TEST_CASE("LHE preserves exact resonance fusion kinematics and central baryons", "[MLHE][kinematics]") {
  auto event = FusionEvent();
  for (const auto unit : {HepMC3::Units::GEV, HepMC3::Units::MEV}) {
    event.set_units(unit, HepMC3::Units::MM);
    const auto rows = gra::BuildLHERows(event);
    REQUIRE(rows.size() == 5);
    CHECK(rows[0].momentum[0] == Approx(1.0));
    CHECK(rows[1].momentum[0] == Approx(-1.0));
    CHECK(rows[0].momentum[4] < 0.0);
    CHECK(rows[1].momentum[4] == Approx(rows[0].momentum[4]));
    std::array<double, 4> balance{};
    int baryons = 0;
    for (const auto &row : rows) {
      const auto &p = row.momentum;
      CHECK(p[3] * p[3] - p[0] * p[0] - p[1] * p[1] - p[2] * p[2] == Approx(p[4] * std::abs(p[4])));
      if (row.status == 1) { ++baryons; CHECK(std::abs(row.pid) == 2212); }
      if (row.status != 1 && row.status != -1) { continue; }
      for (const auto i : indices(balance)) { balance[i] += row.status * p[i]; }
    }
    CHECK(baryons == 2);
    for (const double value : balance) { CHECK(value == Approx(0.0).margin(1e-12)); }
  }
}

// Check LHE four-vectors transform covariantly under azimuthal rotation and beam exchange
TEST_CASE("LHE fusion rows retain rotations and beam exchange", "[MLHE][kinematics][symmetry]") {
  const bool exchange = GENERATE(false, true);
  auto event = FusionEvent();
  const auto reference = gra::BuildLHERows(event);
  for (auto &particle : event.particles()) {
    auto momentum = gra::aux::HepMC2M4Vec(particle->momentum());
    momentum.RotateZ(0.37);
    if (exchange) { momentum.RotateY(std::acos(-1.0)); }
    particle->set_momentum(gra::aux::M4Vec2HepMC3(momentum));
  }
  const auto rows = gra::BuildLHERows(event);
  REQUIRE(rows.size() == reference.size());
  for (const auto i : indices(rows)) {
    const auto j = exchange && i < 2 ? 1 - i : i;
    const auto &p = reference[j].momentum;
    gra::M4Vec momentum(p[0], p[1], p[2], p[3]);
    momentum.RotateZ(0.37);
    if (exchange) { momentum.RotateY(std::acos(-1.0)); }
    const std::array<double, 4> expected{momentum.Px(), momentum.Py(), momentum.Pz(), momentum.E()};
    for (const auto k : indices(expected)) { CHECK(rows[i].momentum[k] == Approx(expected[k]).margin(1e-12)); }
    CHECK(rows[i].momentum[4] == Approx(p[4]));
  }
}

// Check supplied masses survive ultrarelativistic four-vector cancellation and unit changes
TEST_CASE("LHE retains generated masses at large boosts", "[MLHE][units]") {
  auto event = LeptonEvent();
  for (auto &particle : event.particles()) {
    const double sign = particle->momentum().pz() > 0.0 ? 1.0 : -1.0;
    particle->set_momentum(HepMC3::FourVector(0, 0, sign * 1e9, 1e9));
    particle->set_generated_mass(particle->status() == 4 ? 0.9383 : 0.1057);
  }
  for (const auto unit : {HepMC3::Units::GEV, HepMC3::Units::MEV}) {
    event.set_units(unit, HepMC3::Units::MM);
    for (const auto &row : gra::BuildLHERows(event)) {
      CHECK(row.momentum[4] == Approx(row.status == -1 ? 0.9383 : 0.1057));
    }
  }
}

// Check a physical t-channel photon keeps its spacelike LHA status and mass
TEST_CASE("LHE marks spacelike propagators separately from resonances", "[MLHE][kinematics]") {
  HepMC3::GenEvent event(HepMC3::Units::GEV, HepMC3::Units::MM);
  const double mass = 0.000511;
  const double beam_pz = std::sqrt(100.0 - mass * mass);
  const double final_pz = std::sqrt(99.0 - mass * mass);
  auto photon = std::make_shared<HepMC3::GenParticle>(HepMC3::FourVector(-1, 0, beam_pz - final_pz, 0),
                                                     22, gra::PDG::PDG_INTERMEDIATE);
  for (const int sign : {1, -1}) {
    auto vertex = std::make_shared<HepMC3::GenVertex>();
    vertex->add_particle_in(std::make_shared<HepMC3::GenParticle>(HepMC3::FourVector(0, 0, sign * beam_pz, 10), 11, 4));
    vertex->add_particle_out(std::make_shared<HepMC3::GenParticle>(HepMC3::FourVector(sign, 0, sign * final_pz, 10), 11, 1));
    if (sign > 0) { vertex->add_particle_out(photon); }
    else { vertex->add_particle_in(photon); }
    event.add_vertex(vertex);
  }
  const auto rows = gra::BuildLHERows(event);
  REQUIRE(rows.size() == 5);
  for (const auto &row : rows) {
    if (row.pid != 22) { continue; }
    CHECK(row.status == -2);
    CHECK(row.momentum[4] == Approx(photon->momentum().m()));
    CHECK(row.momentum[4] < 0.0);
  }
}

// Check the sampled proper lifetime and decay ancestry through the real LHE writer
TEST_CASE("LHE retains displaced decay lengths in millimetres", "[MLHE][lifetime]") {
  HepMC3::GenEvent event(HepMC3::Units::GEV, HepMC3::Units::CM);
  const double mass = 0.1350;
  auto production = std::make_shared<HepMC3::GenVertex>(HepMC3::FourVector(1, 2, 3, 4));
  for (const int sign : {-1, 1}) {
    production->add_particle_in(std::make_shared<HepMC3::GenParticle>(HepMC3::FourVector(0, 0, sign * 1.25 * mass, 1.25 * mass), 22, 4));
  }
  auto pion = std::make_shared<HepMC3::GenParticle>(HepMC3::FourVector(0, 0, 0.75 * mass, 1.25 * mass), 111, 2);
  pion->set_generated_mass(mass);
  production->add_particle_out(pion);
  production->add_particle_out(std::make_shared<HepMC3::GenParticle>(HepMC3::FourVector(0, 0, -0.75 * mass, 1.25 * mass), 111, 1));
  event.add_vertex(production);
  auto decay = std::make_shared<HepMC3::GenVertex>(HepMC3::FourVector(1, 2, 3.3, 4.5));
  decay->add_particle_in(pion);
  for (const int sign : {-1, 1}) {
    decay->add_particle_out(std::make_shared<HepMC3::GenParticle>(HepMC3::FourVector(sign * 0.5 * mass, 0, 0.375 * mass, 0.625 * mass), 22, 1));
  }
  event.add_vertex(decay);
  auto xs = std::make_shared<HepMC3::GenCrossSection>();
  xs->set_cross_section(1.0, 0.0);
  event.set_cross_section(xs);
  for (const auto unit : {HepMC3::Units::CM, HepMC3::Units::MM}) {
    event.set_units(HepMC3::Units::GEV, unit);
    const auto path = OutputPath("lifetime_" + HepMC3::Units::name(unit) + ".lhe");
    gra::MLHEWriter writer(path);
    writer.WriteEvent(event);
    writer.Close();
    LHEF::Reader reader(path);
    REQUIRE(reader.readEvent());
    int mother = 0;
    for (const auto i : indices(reader.hepeup.IDUP)) {
      if (reader.hepeup.IDUP[i] == 111 && reader.hepeup.ISTUP[i] == 2) {
        mother = i + 1;
        CHECK(reader.hepeup.VTIMUP[i] == Approx(4.0));
      }
      if (reader.hepeup.IDUP[i] == 22 && reader.hepeup.ISTUP[i] == 1) {
        REQUIRE(mother > 0);
        CHECK(reader.hepeup.MOTHUP[i].first == mother);
      }
    }
    CHECK(mother > 0);
  }
}
