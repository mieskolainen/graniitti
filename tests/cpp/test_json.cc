// Generic JSON references and command line edits of shared parameters
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "catch.hpp"

#include <filesystem>
#include <fstream>
#include <future>

#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MJson.h"
#include "Graniitti/Tech/MJsonOverride.h"

namespace {

// Locate the reference fixtures shared with the Python tests
std::string Fixture(const std::string &name) {
  return gra::aux::ResolveProjectPath("tests/data/json_refs/" + name);
}

}  // namespace

TEST_CASE("JSON references expand arbitrary values and paths", "[json][references]") {
  const auto expected = nlohmann::json::parse(gra::aux::ReadFile(Fixture("expected.json")));
  auto value = nlohmann::json::parse(gra::aux::GetInputData(Fixture("valid.json")));
  REQUIRE(value == expected);
  value["object"]["rows"][1][0] = 99;
  CHECK(value["array"] == expected["array"]);
  std::vector<std::future<nlohmann::json>> reads;
  for (int i = 0; i < 4; ++i) {
    reads.push_back(std::async(std::launch::async, [] { return gra::MJson{}.Read(Fixture("valid.json")); }));
  }
  for (auto &read : reads) { CHECK(read.get() == expected); }
}

TEST_CASE("JSON references reject invalid targets and cycles", "[json][references]") {
  const auto cases = nlohmann::json::parse(gra::aux::ReadFile(Fixture("invalid.json")));
  for (const auto &test : cases) {
    CAPTURE(test.at("name"));
    REQUIRE_THROWS_AS(gra::MJson{}.Resolve(test.at("input"), Fixture("input.json")), std::invalid_argument);
  }
}

TEST_CASE("JSON model card overrides follow references across files", "[json][references]") {
  const auto dir = std::filesystem::path(gra::aux::ResolveProjectPath("tmp/test_json"));
  std::filesystem::create_directories(dir);
  const auto card = (dir / "card.json").string();
  const auto source = (dir / "source.json").string();
  std::ofstream(source) << R"({"values": {"x": 1, "array": [2,3]}})";
  std::ofstream(card) << R"({"a": {"$ref":"source.json#/values"}, "b": {"$ref":"#/a"}})";
  using namespace gra::json_override;
  ClearCardOverrides();
  RegisterCardOverrides({ParseSpec("card.json:a.x=7"), ParseSpec("card.json:b.array[0]=8")});
  const auto value = nlohmann::json::parse(gra::aux::GetInputData(card));
  CHECK(value["a"]["x"] == 7);
  CHECK(value["a"]["array"][0] == 8);
  CHECK(value["a"] == value["b"]);
  CHECK(nlohmann::json::parse(gra::aux::GetInputData(source))["values"] == value["a"]);
  RequireAllCardOverridesApplied();
  ClearCardOverrides();
  CHECK(gra::MJson{}.Read(card)["a"]["x"] == 1);
}
