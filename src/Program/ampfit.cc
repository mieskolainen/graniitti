// Prepare coherent helicity amplitudes on common HepMC3 events for amplitude fits
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "HepMC3/ReaderAscii.h"
#include "Graniitti/MGraniitti.h"
#include "Graniitti/Program/ampfit.h"
#include "Graniitti/Tech/MAux.h"
#include "cxxopts.hpp"

using gra::aux::indices;

namespace {

// Buffer complex rows without assuming the binary layout of std::complex
class Amplitudes {
 public:
  std::ofstream output;
  std::vector<double> values;
  std::size_t batch;
  std::size_t rows = 0;

  // Open the requested output with a positive event batch size
  Amplitudes(const std::string &path, std::size_t batch_events) : output(path, std::ios::binary), batch(batch_events) {
    if (!output) { throw std::runtime_error("Cannot open amplitude output: " + path); }
  }

  // Write the completed batch and retain its allocated storage
  void Flush() {
    output.write(reinterpret_cast<const char *>(values.data()), static_cast<std::streamsize>(values.size() * sizeof(double)));
    if (!output) { throw std::runtime_error("Cannot write amplitude output"); }
    values.clear();
    rows = 0;
  }

  // Append one physical external helicity row
  void Add(const std::vector<std::complex<double>> &amplitude) {
    for (const auto &value : amplitude) {
      values.push_back(value.real());
      values.push_back(value.imag());
    }
    if (++rows == batch) { Flush(); }
  }
};

// Evaluate one initialized component through the common generator failure bookkeeping
nlohmann::json Evaluate(const std::string &card, const std::string &events, Amplitudes &output,
                        std::ofstream &intensities, std::ofstream &kinematics, std::size_t &helicities,
                        bool source) {
  gra::MGraniitti generator;
  generator.ReadInput(nlohmann::json::parse(gra::aux::GetInputData(card)));
  generator.proc->PrepareRun();
  const auto *factorized = dynamic_cast<const gra::MFactorized *>(generator.proc);
  if (factorized == nullptr) { throw std::invalid_argument("ampfit requires factorized central kinematics"); }
  gra::program::AmplitudeProcess process(*factorized);
  nlohmann::json daughters = nlohmann::json::array();
  for (const auto &leg : process.state.lts.decaytree) {
    daughters.push_back({{"pdg", leg.p.pdg}, {"mass", leg.p.mass}});
  }
  HepMC3::ReaderAscii input(events);
  std::size_t count = 0, empty = 0, kinematic_failures = 0, amplitude_failures = 0;
  double closure = 0.0;
  while (!input.failed()) {
    HepMC3::GenEvent event;
    input.read_event(event);
    if (input.failed()) { break; }
    event.set_units(HepMC3::Units::GEV, HepMC3::Units::MM);
    gra::MEventWeightState weight;
    double intensity = 0.0;
    std::vector<std::complex<double>> amplitude;
    try {
      amplitude = process.Evaluate(event, weight, intensity);
    } catch (const gra::PhaseSpaceFailure &) {
      weight.kinematics_ok = false;
    }
    if (!weight.Valid()) {
      if (!weight.kinematics_ok) { ++kinematic_failures; } else { ++amplitude_failures; }
      intensity = 0.0;
    }
    intensities.write(reinterpret_cast<const char *>(&intensity), sizeof(double));
    if (source) {
      const auto &lts = process.state.lts;
      const std::array<double, 5> invariants = weight.Valid()
          ? std::array<double, 5>{lts.m2, lts.t1, lts.t2, lts.t_hat, lts.u_hat} : std::array<double, 5>{};
      kinematics.write(reinterpret_cast<const char *>(invariants.data()), sizeof(invariants));
    }
    ++count;
    if (helicities == 0 && !weight.Valid()) { ++empty; continue; }
    if (helicities == 0) {
      helicities = amplitude.size();
      for (std::size_t index = 0; index < empty; ++index) {
        output.Add(std::vector<std::complex<double>>(helicities, 0.0));
      }
    }
    if (!weight.Valid()) { amplitude.assign(helicities, 0.0); }
    if (amplitude.size() != helicities) { throw std::invalid_argument("ampfit cards have different external helicities"); }
    const double norm = gra::SquaredNorm(amplitude);
    const double scale = std::max(norm, intensity);
    if (scale > 0.0) { closure = std::max(closure, std::abs(norm - intensity) / scale); }
    output.Add(amplitude);
  }
  output.Flush();
  return {{"events", count}, {"kinematic_failures", kinematic_failures}, {"amplitude_failures", amplitude_failures},
          {"closure", closure}, {"daughters", daughters}};
}

// Validate CLI preparation controls before opening inputs or outputs
void Validate(const std::vector<std::string> &cards, const std::string &events, std::size_t batch, double tolerance) {
  if (cards.empty()) { throw std::invalid_argument("ampfit requires at least one --input card"); }
  if (batch == 0 || batch > std::numeric_limits<std::size_t>::max() / sizeof(double)) {
    throw std::invalid_argument("--batch-events must be a positive event count");
  }
  if (!std::isfinite(tolerance) || tolerance <= 0.0) {
    throw std::invalid_argument("--closure-rtol must be finite and positive");
  }
  for (const auto &path : cards) {
    if (!std::ifstream(path)) { throw std::invalid_argument("Cannot open input card: " + path); }
  }
  if (!std::ifstream(events)) { throw std::invalid_argument("Cannot open HepMC3 events: " + events); }
}

}  // namespace

// Prepare ordered basis amplitudes using ordinary generator cards and explicit CLI controls
int main(int argc, char **argv) {
  try {
    cxxopts::Options options(argv[0], "Prepare complex amplitudes with the normal GRANIITTI process");
    options.add_options()
        ("i,input", "Generator card, repeat in basis order with the source card first", cxxopts::value<std::string>())
        ("e,events", "Common HepMC3 event sample", cxxopts::value<std::string>())
        ("o,output", "Binary amplitude output", cxxopts::value<std::string>())
        ("reference", "Prepared source amplitude metadata", cxxopts::value<std::string>())
        ("batch-events", "Buffered events per output batch", cxxopts::value<std::size_t>())
        ("closure-rtol", "Relative tolerance for the physical helicity norm", cxxopts::value<double>())
        ("h,help", "Print usage");
    const auto args = options.parse(argc, argv);
    if (args.count("help") || argc == 1) { std::cout << options.help() << '\n'; return 0; }
    std::vector<std::string> cards;
    for (const auto &arg : args.arguments()) {
      if (arg.key() == "input") { cards.push_back(arg.value()); }
    }
    const auto events = args["events"].as<std::string>();
    const auto path = args["output"].as<std::string>();
    const auto batch = args["batch-events"].as<std::size_t>();
    const auto tolerance = args["closure-rtol"].as<double>();
    Validate(cards, events, batch, tolerance);
    const auto reference = args.count("reference")
        ? nlohmann::json::parse(gra::aux::GetInputData(args["reference"].as<std::string>())) : nlohmann::json();
    std::size_t count = reference.is_null() ? 0 : reference.at("shape").at(1).get<std::size_t>();
    std::size_t helicities = reference.is_null() ? 0 : reference.at("shape").at(2).get<std::size_t>();
    if (!reference.is_null() && (count == 0 || helicities == 0)) {
      throw std::invalid_argument("ampfit reference has no events or helicities");
    }
    gra::aux::CreateDirectory(std::filesystem::path(path).parent_path().string());
    Amplitudes output(path, batch);
    std::ofstream intensities(path + ".intensity", std::ios::binary);
    std::ofstream kinematics(path + ".kinematics", std::ios::binary);
    if (!intensities || !kinematics) { throw std::runtime_error("Cannot open amplitude normalization outputs"); }
    nlohmann::json components = nlohmann::json::array(), failures = nlohmann::json::array();
    for (const auto index : indices(cards)) {
      const auto result = Evaluate(cards[index], events, output, intensities, kinematics, helicities, index == 0);
      const auto size = result.at("events").get<std::size_t>();
      if (index == 0 && reference.is_null()) { count = size; }
      if (count == 0 || size != count || helicities == 0) {
        throw std::invalid_argument("ampfit event counts differ or no physical amplitudes were found");
      }
      if (result.at("closure").get<double>() > tolerance) {
        throw std::runtime_error("ampfit helicity norm disagrees with the generator intensity");
      }
      if (!reference.is_null() && result.at("daughters") != reference.at("daughters")) {
        throw std::invalid_argument("ampfit reference has different final particles");
      }
      failures.push_back(result.at("kinematic_failures").get<std::size_t>() + result.at("amplitude_failures").get<std::size_t>());
      components.push_back(result);
      std::cout << "ampfit component " << index + 1 << '/' << cards.size() << ", events " << count
                << ", failures " << failures.back() << '\n';
    }
    output.output.close();
    intensities.close();
    kinematics.close();
    if (!output.output || !intensities || !kinematics) { throw std::runtime_error("Cannot finish amplitude outputs"); }
    const nlohmann::json metadata = {{"shape", {cards.size(), count, helicities}}, {"failures", failures},
        {"components", components}, {"daughters", components.front().at("daughters")},
        {"layout", "component,event,helicity,re_im"}};
    std::ofstream metadata_output(path + ".json");
    metadata_output << metadata.dump(2) << '\n';
    metadata_output.close();
    if (!metadata_output) { throw std::runtime_error("Cannot finish amplitude metadata"); }
    return 0;
  } catch (const std::exception &error) {
    std::cerr << "ampfit: " << error.what() << '\n';
    return 1;
  }
}

