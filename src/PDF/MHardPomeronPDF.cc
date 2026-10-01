// Hard diffractive Pomeron PDF access
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

// Own
#include "Graniitti/PDF/MHardPomeronPDF.h"
#include "Graniitti/PDF/MLHAPDF.h"
#include "Graniitti/Particle/MForm.h"
#include "Graniitti/Tech/MAux.h"

// Other
#include "json.hpp"

using gra::aux::indices;

namespace gra {

// Validate model-card ranges before using LHAPDF
void MHardPomeronPDFParam::Validate() const {
  if (!std::isfinite(remnant_mass) || !(remnant_mass > 0.0)) {
    throw std::invalid_argument("MHardPomeronPDFParam::Validate: invalid remnant_mass");
  }
  if (DPDF_SET.empty() || DPDF_SET == "null") {
    throw std::invalid_argument(
        "MHardPomeronPDFParam::Validate: invalid DPDF_SET");
  }
  if (DPDF_MEMBER < 0) {
    throw std::invalid_argument(
        "MHardPomeronPDFParam::Validate: negative DPDF_MEMBER");
  }
  if (xi_range.size() != 2 || !std::isfinite(xi_range[0]) ||
      !std::isfinite(xi_range[1]) || !(xi_range[0] > 0.0) ||
      !(xi_range[0] < xi_range[1]) || !(xi_range[1] < 1.0)) {
    throw std::invalid_argument(
        "MHardPomeronPDFParam::Validate: invalid xi_range");
  }
  if (beta_range.size() != 2 || !std::isfinite(beta_range[0]) ||
      !std::isfinite(beta_range[1]) || !(beta_range[0] > 0.0) ||
      !(beta_range[0] < beta_range[1]) || beta_range[1] > 1.0) {
    throw std::invalid_argument(
        "MHardPomeronPDFParam::Validate: invalid beta_range");
  }
  if (t_range.size() != 2 || !std::isfinite(t_range[0]) ||
      !std::isfinite(t_range[1]) || !(t_range[0] < t_range[1]) ||
      t_range[1] > 0.0) {
    throw std::invalid_argument(
        "MHardPomeronPDFParam::Validate: invalid t_range");
  }
  if (!(Q2_min > 0.0) || !std::isfinite(Q2_min)) {
    throw std::invalid_argument(
        "MHardPomeronPDFParam::Validate: invalid Q2_min");
  }
  if (!std::isfinite(alpha0) || !std::isfinite(alpha_prime) ||
      !std::isfinite(B_flux) || !std::isfinite(norm)) {
    throw std::invalid_argument(
        "MHardPomeronPDFParam::Validate: non-finite parameter");
  }
  if (!(norm > 0.0)) {
    throw std::invalid_argument(
        "MHardPomeronPDFParam::Validate: non-positive norm");
  }
  if (parton_flavours.empty()) {
    throw std::invalid_argument(
        "MHardPomeronPDFParam::Validate: empty parton_flavours");
  }
  std::vector<int> unique_flavours = parton_flavours;
  std::sort(unique_flavours.begin(), unique_flavours.end());
  if (std::adjacent_find(unique_flavours.begin(), unique_flavours.end()) !=
      unique_flavours.end()) {
    throw std::invalid_argument(
        "MHardPomeronPDFParam::Validate: duplicate parton_flavours");
  }
  for (const int pid : parton_flavours) {
    if (pid != 21 && (std::abs(pid) < 1 || std::abs(pid) > 5)) {
      throw std::invalid_argument(
          "MHardPomeronPDFParam::Validate: invalid parton_flavours entry");
    }
  }
}

// Construct with the default model-card parameters
MHardPomeronPDF::MHardPomeronPDF() { InitPDF(); }

// Construct from a GRANIITTI GENERAL.json model card
MHardPomeronPDF::MHardPomeronPDF(const std::string &modelfile) {
  ReadParameters(modelfile);
}

// Construct from one already captured immutable model snapshot
MHardPomeronPDF::MHardPomeronPDF(const SoftModelPtr &soft_model) {
  if (soft_model == nullptr) {
    throw std::invalid_argument("MHardPomeronPDF: missing SOFT model");
  }
  ConfigureFromJson(soft_model->SourceFile(), soft_model->SourceJson());
  InitPDF();
}

// Construct from one model snapshot and a run owned LHAPDF store
MHardPomeronPDF::MHardPomeronPDF(const SoftModelPtr &soft_model,
                                 MLHAPDFStore &pdf_store) {
  if (soft_model == nullptr) {
    throw std::invalid_argument("MHardPomeronPDF: missing SOFT model");
  }
  ConfigureFromJson(soft_model->SourceFile(), soft_model->SourceJson());
  InitPDF(pdf_store);
}

// Read parameters and load a changed PDF before replacing the current state
void MHardPomeronPDF::ReadParameters(const std::string &modelfile) {
  std::string json_text;
  try {
    json_text = gra::aux::GetInputData(modelfile);
  } catch (const std::exception &e) {
    throw std::invalid_argument(
        "MHardPomeronPDF::ReadParameters: Error reading " + modelfile + ": " +
        e.what());
  }
  MHardPomeronPDF selected(*this);
  selected.ConfigureFromJson(modelfile, json_text);
  if (selected.dpdf == nullptr || selected.param.DPDF_SET != param.DPDF_SET ||
      selected.param.DPDF_MEMBER != param.DPDF_MEMBER) {
    selected.InitPDF();
  }
  param = std::move(selected.param);
  dpdf = std::move(selected.dpdf);
}

// Read and validate parameters from one exact GENERAL JSON text
void MHardPomeronPDF::ConfigureFromJson(const std::string &source_file,
                                        const std::string &json_text) {
  MHardPomeronPDFParam selected = param;
  try {
    const nlohmann::json j = nlohmann::json::parse(json_text);

    if (j.contains("PARAM_HARDPOMERON")) {
      const auto &block = j.at("PARAM_HARDPOMERON");
      selected.remnant_mass = block.at("remnant_mass");
      if (block.contains("DPDF_SET")) {
        selected.DPDF_SET = block.at("DPDF_SET");
      }
      if (block.contains("DPDF_MEMBER")) {
        selected.DPDF_MEMBER = block.at("DPDF_MEMBER");
      }
      if (block.contains("alpha0")) {
        selected.alpha0 = block.at("alpha0");
      }
      if (block.contains("alpha_prime")) {
        selected.alpha_prime = block.at("alpha_prime");
      }
      if (block.contains("B_flux")) {
        selected.B_flux = block.at("B_flux");
      }
      if (block.contains("norm")) {
        selected.norm = block.at("norm");
      }
      if (block.contains("parton_flavours")) {
        selected.parton_flavours =
            block.at("parton_flavours").get<std::vector<int>>();
      }
      if (block.contains("xi_range")) {
        selected.xi_range = block.at("xi_range").get<std::vector<double>>();
      }
      if (block.contains("beta_range")) {
        selected.beta_range = block.at("beta_range").get<std::vector<double>>();
      }
      if (block.contains("t_range")) {
        selected.t_range = block.at("t_range").get<std::vector<double>>();
      }
      if (block.contains("Q2_min")) {
        selected.Q2_min = block.at("Q2_min");
      }
    }
  } catch (const nlohmann::json::exception &e) {
    throw std::invalid_argument(
        "MHardPomeronPDF::ReadParameters: Error parsing " + source_file + ": " +
        e.what());
  }

  selected.Validate();
  param = std::move(selected);
}

// Compute the ordinary proton remnant mass in GeV
double MHardPomeronPDF::RemnantMass() const { return param.remnant_mass; }

// Read only the configured hard-process parton flavours without loading LHAPDF
std::vector<int>
MHardPomeronPDF::ReadPartonFlavours(const std::string &modelfile) {
  MHardPomeronPDFParam selected;
  try {
    const nlohmann::json j =
        nlohmann::json::parse(gra::aux::GetInputData(modelfile));
    if (j.contains("PARAM_HARDPOMERON") &&
        j.at("PARAM_HARDPOMERON").contains("parton_flavours")) {
      selected.parton_flavours = j.at("PARAM_HARDPOMERON")
                                     .at("parton_flavours")
                                     .get<std::vector<int>>();
    }
  } catch (const std::exception &e) {
    throw std::invalid_argument(
        "MHardPomeronPDF::ReadPartonFlavours: Error reading " + modelfile +
        ": " + e.what());
  }
  selected.Validate();
  return selected.parton_flavours;
}

// Compute the hard Pomeron trajectory
// alpha_P(t) = alpha_P(0) + alpha'_P t
double MHardPomeronPDF::AlphaP(double t) const {
  return param.alpha0 + param.alpha_prime * t;
}

// Compute the Pomeron flux in the proton
// f_IP/p(xi,t) = N exp(B t) xi^[1-2 alpha_P(t)]
double MHardPomeronPDF::Flux(double xi, double t) const {
  if (!InRange(xi, param.xi_range) || !InRange(t, param.t_range)) {
    return 0.0;
  }

  const double alpha = AlphaP(t);
  const double flux = param.norm * gra::form::ExpSlopeWeight(param.B_flux, t) /
                      std::pow(xi, 2.0 * alpha - 1.0);
  return std::isfinite(flux) ? flux : 0.0;
}

// Compute the parton density in the Pomeron from the LHAPDF xfx value
// f_i/IP(beta,Q2) = [beta f_i/IP(beta,Q2)]_LHAPDF / beta
double MHardPomeronPDF::PartonDensity(int pid, double beta, double Q2) const {
  if (!InRange(beta, param.beta_range) || !std::isfinite(Q2) ||
      !(Q2 >= param.Q2_min)) {
    return 0.0;
  }
  if (dpdf == nullptr) {
    throw std::invalid_argument(
        "MHardPomeronPDF::PartonDensity: DPDF is not initialized");
  }
  try {
    if (!dpdf->inRangeXQ2(beta, Q2)) {
      return 0.0;
    }
    const double xfx = dpdf->xfxQ2(pid, beta, Q2);
    const double f = xfx / beta;
    return std::isfinite(f) ? f : 0.0;
  } catch (...) {
    return 0.0;
  }
}

// Compute the factorized diffractive parton density
// f_i/p^D = f_IP/p(xi,t) f_i/IP(beta,Q2)
double MHardPomeronPDF::DiffractiveDensity(int pid, double xi, double beta,
                                           double t, double Q2) const {
  return Flux(xi, t) * PartonDensity(pid, beta, Q2);
}

// Compute x_hard = xi beta for event-record bookkeeping
double MHardPomeronPDF::HardX(double xi, double beta) const {
  return xi * beta;
}

// Compute max(hard_scale2, Q2_min) for DPDF evaluation
double MHardPomeronPDF::FactorizationQ2(double hard_scale2) const {
  return std::max(hard_scale2, param.Q2_min);
}

// Compute alpha_s from the same DPDF member used for diffractive densities
double MHardPomeronPDF::AlphaS(double Q2) const {
  if (dpdf == nullptr) {
    throw std::invalid_argument(
        "MHardPomeronPDF::AlphaS: DPDF is not initialized");
  }
  if (!(Q2 > 0.0) || !std::isfinite(Q2)) {
    throw std::invalid_argument("MHardPomeronPDF::AlphaS: invalid Q2");
  }
  if (!dpdf->hasAlphaS() || !dpdf->inRangeQ2(Q2)) {
    throw std::invalid_argument(
        "MHardPomeronPDF::AlphaS: Q2 is outside the DPDF alpha_s grid");
  }
  const double alpha_s = dpdf->alphasQ2(Q2);
  if (!(alpha_s > 0.0) || !std::isfinite(alpha_s)) {
    throw std::invalid_argument(
        "MHardPomeronPDF::AlphaS: DPDF returned invalid alpha_s");
  }
  return alpha_s;
}

// Compute whether one proton PDF uses the same alpha_s evolution
bool MHardPomeronPDF::MatchAlphaS(const LHAPDF::PDF &pdf, double Q2) const {
  if (dpdf == nullptr || !dpdf->hasAlphaS() || !pdf.hasAlphaS() ||
      !(Q2 > 0.0) || !std::isfinite(Q2) || !dpdf->inRangeQ2(Q2) ||
      !pdf.inRangeQ2(Q2)) {
    return false;
  }
  const auto &dpdf_info = dpdf->info();
  const auto &pdf_info = pdf.info();
  if (dpdf_info.get_entry_as<int>("OrderQCD", -1) !=
          pdf_info.get_entry_as<int>("OrderQCD", -2) ||
      dpdf_info.get_entry_as<int>("AlphaS_OrderQCD", -1) !=
          pdf_info.get_entry_as<int>("AlphaS_OrderQCD", -2)) {
    return false;
  }
  try {
    const double dpdf_alpha = dpdf->alphasQ2(Q2);
    const double pdf_alpha = pdf.alphasQ2(Q2);
    const double scale = std::max(std::abs(dpdf_alpha), std::abs(pdf_alpha));
    return std::isfinite(dpdf_alpha) && std::isfinite(pdf_alpha) &&
           dpdf_alpha > 0.0 && pdf_alpha > 0.0 && scale > 0.0 &&
           std::abs(dpdf_alpha - pdf_alpha) <= 1.0e-3 * scale;
  } catch (...) {
    return false;
  }
}

// Validate PDF support and report mixed coupling evolution once before sampling
void MHardPomeronPDF::ValidateAlphaS(const LHAPDF::PDF &pdf, double Q2_min, double Q2_max) const {
  if (!(Q2_min > 0.0) || !(Q2_max >= Q2_min) || !std::isfinite(Q2_max)) {
    throw std::invalid_argument("MHardPomeronPDF: invalid hard-scale interval");
  }
  if (dpdf == nullptr || !dpdf->hasAlphaS() || !pdf.hasAlphaS() ||
      !dpdf->inRangeQ2(Q2_min) || !dpdf->inRangeQ2(Q2_max) ||
      !pdf.inRangeQ2(Q2_min) || !pdf.inRangeQ2(Q2_max) ||
      dpdf->info().get_entry_as<int>("OrderQCD", -1) != pdf.info().get_entry_as<int>("OrderQCD", -2)) {
    throw std::invalid_argument("IPp requires PDF and DPDF scale support and matching QCD order");
  }
  std::vector<double> scale = {Q2_min, Q2_max};
  for (const auto *member : {dpdf.get(), &pdf}) {
    const auto &info = member->info();
    auto knots = info.get_entry_as<std::vector<double>>("AlphaS_Qs", {});
    for (const auto *quark : {"MCharm", "MBottom", "MTop"}) {
      knots.push_back(info.get_entry_as<double>(quark, 0.0));
    }
    for (const double Q : knots) {
      const double Q2 = Q * Q;
      if (Q2 > Q2_min && Q2 < Q2_max) { scale.push_back(Q2); }
    }
  }
  std::sort(scale.begin(), scale.end());
  for (const auto &i : indices(scale)) {
    const double midpoint = i > 0 ? std::sqrt(scale[i - 1]) * std::sqrt(scale[i]) : scale[i];
    if (!MatchAlphaS(pdf, scale[i]) || !MatchAlphaS(pdf, midpoint)) {
      const double Q2 = !MatchAlphaS(pdf, scale[i]) ? scale[i] : midpoint;
      std::ostringstream message;
      message << std::fixed << std::setprecision(4)
              << "IPp: mixed PDF and DPDF alpha_s evolution at Q = " << std::sqrt(Q2) << " GeV\n"
              << "  Proton PDF: " << pdf.set().name() << " / " << pdf.memberID()
              << ", alpha_s(Q) = " << pdf.alphasQ2(Q2) << '\n'
              << "  DPDF:       " << param.DPDF_SET << " / " << param.DPDF_MEMBER
              << ", alpha_s(Q) = " << dpdf->alphasQ2(Q2) << '\n'
              << "  The hard amplitude retains DPDF alpha_s. This is a mixed PDF model choice.";
      aux::PrintWarning();
      std::cout << message.str() << std::endl;
      return;
    }
  }
}

// Map one unit random variable to the configured xi range
double MHardPomeronPDF::MapXi(double unit) const {
  return param.xi_range[0] + unit * (param.xi_range[1] - param.xi_range[0]);
}

// Map one unit random variable to the configured beta range
double MHardPomeronPDF::MapBeta(double unit) const {
  return param.beta_range[0] +
         unit * (param.beta_range[1] - param.beta_range[0]);
}

// Map one unit random variable to the configured t range
double MHardPomeronPDF::MapT(double unit) const {
  return param.t_range[0] + unit * (param.t_range[1] - param.t_range[0]);
}

// Compute the exponential momentum-transfer slope of the flux at fixed xi
// d ln(f_IP/p)/dt = B - 2 alpha'_P ln(xi)
double MHardPomeronPDF::FluxTSlope(double xi) const {
  if (!InRange(xi, param.xi_range)) {
    return 0.0;
  }
  return param.B_flux - 2.0 * param.alpha_prime * std::log(xi);
}

// Access the configured xi support
const std::vector<double> &MHardPomeronPDF::XiRange() const {
  return param.xi_range;
}

// Access the configured beta support
const std::vector<double> &MHardPomeronPDF::BetaRange() const {
  return param.beta_range;
}

// Access the configured t support
const std::vector<double> &MHardPomeronPDF::TRange() const {
  return param.t_range;
}

// Compute the integration volume for one diffractive proton leg
// V = Delta xi Delta beta Delta t
double MHardPomeronPDF::DomainVolume() const {
  return (param.xi_range[1] - param.xi_range[0]) *
         (param.beta_range[1] - param.beta_range[0]) *
         (param.t_range[1] - param.t_range[0]);
}

// Access the configured hard-process parton flavours
const std::vector<int> &MHardPomeronPDF::PartonFlavours() const {
  return param.parton_flavours;
}

// Print the configured DPDF domain and set metadata
void MHardPomeronPDF::PrintSummary(std::ostream &os) const {
  os << "- DPDF set: " << param.DPDF_SET << std::endl;
  os << "- xi range: [" << param.xi_range[0] << ", " << param.xi_range[1] << "]"
     << std::endl;
  os << "- beta range: [" << param.beta_range[0] << ", " << param.beta_range[1]
     << "]" << std::endl;
  os << "- parton flavours:";
  for (const int pid : param.parton_flavours) {
    os << " " << pid;
  }
  os << std::endl;
}

// Initialize the shared read-only LHAPDF member
void MHardPomeronPDF::InitPDF() {
  MLHAPDFStore pdf_store;
  InitPDF(pdf_store);
}

// Initialize the PDF member through one run owned store
void MHardPomeronPDF::InitPDF(MLHAPDFStore &pdf_store) {
  param.Validate();
  dpdf = pdf_store.GetPDF(param.DPDF_SET, param.DPDF_MEMBER);
}

// Construct a standalone wrapper store with one private PDF store
MHardPomeronPDFStore::MHardPomeronPDFStore()
    : owned_pdf(std::make_unique<MLHAPDFStore>()),
      pdf_store(owned_pdf.get()) {}

// Construct a wrapper store sharing the run owned PDF store
MHardPomeronPDFStore::MHardPomeronPDFStore(MLHAPDFStore &pdf_store_in)
    : pdf_store(&pdf_store_in) {}

// Test whether a value is inside an inclusive two-element range
bool MHardPomeronPDF::InRange(double x,
                              const std::vector<double> &range) const {
  return std::isfinite(x) && range.size() == 2 && x >= range[0] &&
         x <= range[1];
}

// Provide strict weak ordering for hard-Pomeron store keys
bool MHardPomeronPDFStore::Key::operator<(const Key &other) const {
  if (modelfile != other.modelfile) {
    return modelfile < other.modelfile;
  }
  return model_json < other.model_json;
}

// Compute one shared read-only hard-Pomeron PDF wrapper
MHardPomeronPDFStore::HardPomeronPtr
MHardPomeronPDFStore::GetHardPomeronPDF(const std::string &modelfile) {
  if (modelfile.empty()) {
    throw std::invalid_argument(
        "MHardPomeronPDFStore::GetHardPomeronPDF: empty model file path");
  }
  return GetHardPomeronPDF(SoftModel::LoadFromJson(modelfile, aux::GetInputData(modelfile)));
}

// Compute a wrapper from one already captured immutable model snapshot
MHardPomeronPDFStore::HardPomeronPtr
MHardPomeronPDFStore::GetHardPomeronPDF(const SoftModelPtr &soft_model) {
  if (soft_model == nullptr) {
    throw std::invalid_argument(
        "MHardPomeronPDFStore::GetHardPomeronPDF: missing SOFT model");
  }

  const Key key = MakeKey(soft_model);
  return store.GetOrLoad(
      key, [this, &soft_model] { return LoadHardPomeronPDF(soft_model); });
}

// Build the full cache key from one immutable model snapshot
MHardPomeronPDFStore::Key
MHardPomeronPDFStore::MakeKey(const SoftModelPtr &soft_model) const {
  return {soft_model->SourceFile(), soft_model->SourceJson()};
}

// Construct one initialized wrapper from an immutable model snapshot
MHardPomeronPDFStore::HardPomeronPtr
MHardPomeronPDFStore::LoadHardPomeronPDF(const SoftModelPtr &soft_model) const {
  return std::make_shared<MHardPomeronPDF>(soft_model, *pdf_store);
}

} // namespace gra
