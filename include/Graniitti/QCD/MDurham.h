// Durham QCD processes and amplitudes
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MDURHAM_H
#define MDURHAM_H

// C++
#include <array>
#include <cmath>
#include <complex>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <optional>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Own
#include "Graniitti/Amplitude/MG5/Durham/AMP_MG5_DurhamRegistry.h"
#include "Graniitti/Amplitude/MG5/Runtime/AMP_MG5_ProcessRegistry.h"
#include "Graniitti/Kinematics/M4Vec.h"
#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Math/MPolarCut.h"
#include "Graniitti/PDF/MSudakov.h"
#include "Graniitti/Process/MAmpMatch.h"
#include "Graniitti/Process/MProcessSetup.h"
#include "Graniitti/Process/MProcessState.h"
#include "Graniitti/Regge/MSoftModel.h"
#include "Graniitti/Sampling/MRandom.h"
#include "Graniitti/Tech/MAux.h"

namespace gra {

// Select the Durham amplitude family exposed through the process definition
enum class MDurhamMode { Generic, Resonance, Parton, PhotonPair, MesonPair, Flux };

// Compute true for one direct meson pair supported by the Durham hard kernel
bool DurhamMesonPairSupported(const std::vector<MDecayBranch> &tree);

// Durham model parameters
struct MDurhamParam {
  std::string PDF_scale               = "MIN";      // Scheme
  double      muF                     = 1.0;        // Hard factorization scale mu_F = value M_X
  double      muR                     = 1.0;        // Hard alpha_s scale mu_R = value M_X
  double      loop_q2_cut             = 0.4;        // Common Durham virtuality cutoff in GeV^2
  double      MAXCOS                  = 0.9;        // Meson pair amplitude |cos(theta*)| < MAXCOS
  double      MESON_pt_min            = 2.0;        // Per-meson hard-frame pT cutoff in GeV
  double      f_pi                    = 0.1300;     // Pion decay constant in GeV
  double      eta_theta8_deg          = -21.2;      // eta octet decay-constant mixing angle
  double      eta_theta1_deg          = -9.2;       // eta singlet decay-constant mixing angle
  double      f_eta8_over_fpi         = 1.26;       // eta8 decay constant ratio
  double      f_eta0_over_fpi         = 1.17;       // singlet f1 decay constant ratio
  double      chic0_gg_width_fraction = 2.0 / 3.0;  // chi_c0 two-gluon branching fraction
  std::string JET_ALGO                = "anti-kt";  // Partonic MadGraph jet algorithm or none
  double      JET_R                   = 0.4;        // Partonic MadGraph jet radius
  double      JET_pt_min              = 1.0;        // Partonic MadGraph jet minimum pT in GeV
  double      JET_rap_max             = 10.0;       // Partonic MadGraph jet maximum lab-frame |y|

  bool initialized = false;

  // Validate Durham model parameter card values
  void Validate() const {
    if (!std::isfinite(muF) || muF <= 0.0) {
      throw std::invalid_argument("MDurhamParam::Validate: muF must be positive");
    }
    if (!std::isfinite(muR) || muR <= 0.0) {
      throw std::invalid_argument("MDurhamParam::Validate: muR must be positive");
    }
    if (!std::isfinite(loop_q2_cut) || loop_q2_cut < 0.0) {
      throw std::invalid_argument("MDurhamParam::Validate: loop_q2_cut must be non-negative");
    }
    if (!std::isfinite(MAXCOS) || MAXCOS < 0.0 || MAXCOS >= 1.0) {
      throw std::invalid_argument("MDurhamParam::Validate: MAXCOS must be in [0,1)");
    }
    if (!std::isfinite(MESON_pt_min) || MESON_pt_min <= 0.0) {
      throw std::invalid_argument("MDurhamParam::Validate: MESON_pt_min must be positive");
    }
    if (!std::isfinite(f_pi) || f_pi <= 0.0) {
      throw std::invalid_argument("MDurhamParam::Validate: f_pi must be positive");
    }
    if (!std::isfinite(eta_theta8_deg) || !std::isfinite(eta_theta1_deg)) {
      throw std::invalid_argument(
          "MDurhamParam::Validate: eta_theta8_deg and "
          "eta_theta1_deg must be finite");
    }
    if (!std::isfinite(f_eta8_over_fpi) || !std::isfinite(f_eta0_over_fpi) || f_eta8_over_fpi <= 0.0 ||
        f_eta0_over_fpi <= 0.0) {
      throw std::invalid_argument("MDurhamParam::Validate: eta decay-constant ratios must be positive");
    }
    if (!std::isfinite(chic0_gg_width_fraction) || chic0_gg_width_fraction <= 0.0 || chic0_gg_width_fraction > 1.0) {
      throw std::invalid_argument("MDurhamParam::Validate: chic0_gg_width_fraction must be in (0,1]");
    }
    if (JET_ALGO != "none" && JET_ALGO != "anti-kt" && JET_ALGO != "kt" && JET_ALGO != "CA") {
      throw std::invalid_argument("MDurhamParam::Validate: unknown JET_ALGO option " + JET_ALGO);
    }
    if (JET_ALGO != "none") {
      if (!std::isfinite(JET_R) || JET_R <= 0.0) {
        throw std::invalid_argument("MDurhamParam::Validate: JET_R must be positive");
      }
      if (!std::isfinite(JET_pt_min) || JET_pt_min < 0.0) {
        throw std::invalid_argument("MDurhamParam::Validate: JET_pt_min must be non-negative");
      }
      if (!std::isfinite(JET_rap_max) || JET_rap_max <= 0.0) {
        throw std::invalid_argument("MDurhamParam::Validate: JET_rap_max must be positive");
      }
    }
    if (PDF_scale != "MIN" && PDF_scale != "MAX" && PDF_scale != "IN" && PDF_scale != "EX" && PDF_scale != "AVG") {
      throw std::invalid_argument("MDurhamParam::Validate: unknown PDF_scale option " + PDF_scale);
    }
  }

  // Read parameters from one immutable GENERAL JSON text
  void ConfigureFromJson(const std::string &source_file, const std::string &json_text) {
    using json = nlohmann::json;
    json j;

    try {
      j = json::parse(json_text);

      // JSON block identifier
      const std::string XID = "PARAM_DURHAM";
      PDF_scale               = j.at(XID).at("PDF_scale");
      muF                     = j.at(XID).at("muF");
      muR                     = j.at(XID).at("muR");
      loop_q2_cut             = j.at(XID).at("loop_q2_cut");
      MAXCOS                  = j.at(XID).at("MAXCOS");
      MESON_pt_min            = j.at(XID).at("MESON_pt_min");
      f_pi                    = j.at(XID).at("f_pi");
      eta_theta8_deg          = j.at(XID).at("eta_theta8_deg");
      eta_theta1_deg          = j.at(XID).at("eta_theta1_deg");
      f_eta8_over_fpi         = j.at(XID).at("f_eta8_over_fpi");
      f_eta0_over_fpi         = j.at(XID).at("f_eta0_over_fpi");
      chic0_gg_width_fraction = j.at(XID).at("chic0_gg_width_fraction");
      JET_ALGO                = j.at(XID).at("JET_ALGO");
      JET_R                   = j.at(XID).at("JET_R");
      JET_pt_min              = j.at(XID).at("JET_pt_min");
      JET_rap_max             = j.at(XID).at("JET_rap_max");
    } catch (const json::exception &error) {
      throw std::invalid_argument("MDurhamParam::ConfigureFromJson: Error reading " + source_file + ": " +
                                  error.what());
    }

    Validate();

    const std::string XID = "PARAM_DURHAM";
    std::cout << "MDurham::ReadParameters: [PARAM_DURHAM]" << std::endl;
    std::cout << j.at(XID) << std::endl;
    std::cout << std::endl;

    initialized = true;
  }
};

// Physical Durham scales derived from one central-system mass
struct MDurhamScales {
  double muF   = 0.0;
  double muR   = 0.0;
};

// Compute the factorization and renormalization scales
// (muF,muR) = M_X(cF,cR)
inline MDurhamScales DurhamCentralScales(double central_mass, const MDurhamParam &param) {
  return {param.muF * central_mass, param.muR * central_mass};
}

// Durham gluon loop numerical integration parameters
struct MDurhamNumerics {
  math::PolarParam loop;
  std::string      helicity_projector = "";

  double qt2_MAX = 0.0;
  unsigned int N_x     = 0;

  bool initialized = false;

  // Validate Durham loop integration card values
  void Validate(double sudakov_q2_max) const {
    unsigned int multiplier = 0;
    try {
      multiplier = math::PolarIntervalMultiple(loop);
    } catch (const std::invalid_argument &) {
      throw std::invalid_argument("MDurhamNumerics::Validate: Unknown qt_integrator = " + loop.radial_integrator);
    }
    if (loop.azimuth_integrator != "Trap" && loop.azimuth_integrator != "GL") {
      throw std::invalid_argument("MDurhamNumerics::Validate: Unknown phi_integrator = " + loop.azimuth_integrator);
    }
    if (helicity_projector != "transverse" && helicity_projector != "jzp") {
      throw std::invalid_argument("MDurhamNumerics::Validate: Unknown HELICITY_PROJECTOR = " + helicity_projector);
    }
    if (loop.radial_intervals == 0 || loop.azimuth_nodes == 0) {
      throw std::invalid_argument("MDurhamNumerics::Validate: N_qt and N_phi must be positive");
    }
    if (loop.radial_integrator != "GL" && loop.radial_intervals % multiplier != 0) {
      throw std::invalid_argument("MDurhamNumerics::Validate: N_qt = " + std::to_string(loop.radial_intervals) +
                                  " is not compatible with qt_integrator " + loop.radial_integrator);
    }
    if (!std::isfinite(qt2_MAX) || qt2_MAX <= 0.0) {
      throw std::invalid_argument("MDurhamNumerics::Validate: qt2_MAX must be positive");
    }
    if (!std::isfinite(sudakov_q2_max) || sudakov_q2_max <= 0.0) {
      throw std::invalid_argument("MDurhamNumerics::Validate: invalid NUMERICS_SUDAKOV q2_MAX");
    }
    if (qt2_MAX > sudakov_q2_max * (1.0 + 1e-12)) {
      throw std::invalid_argument(
          "MDurhamNumerics::Validate: qt2_MAX must not "
          "exceed NUMERICS_SUDAKOV q2_MAX");
    }
  }

  // Read Durham loop integration parameters from one NUMERICS JSON text
  void ConfigureFromJson(const std::string &source_file, const std::string &json_text) {
    using json = nlohmann::json;
    json j;

    try {
      j = json::parse(json_text);

      const std::string XID           = "NUMERICS_DURHAM";
      qt2_MAX                         = j.at(XID).at("qt2_MAX");
      const std::string qt_integrator = j.at(XID).at("qt_integrator");
      loop.radial_integrator          = qt_integrator == "GL_IR" ? "GL" : qt_integrator;
      loop.radial_map                 = qt_integrator == "GL_IR" ? math::RadialMap::Square : math::RadialMap::Linear;
      if (j.at(XID).at("log_qt").get<bool>()) { loop.radial_map = math::RadialMap::Log; }
      loop.azimuth_integrator         = j.at(XID).at("phi_integrator");
      helicity_projector              = j.at(XID).at("HELICITY_PROJECTOR");

      // Require integer counts before JSON conversion can truncate or wrap
      const auto count = [&](const std::string& key) {
        const auto& value = j.at(XID).at(key);
        if (!value.is_number_integer() || value < 1 || value >= std::numeric_limits<int>::max()) {
          throw std::invalid_argument("MDurhamNumerics: " + key + " must be a positive integer");
        }
        return value.get<unsigned int>();
      };
      loop.radial_intervals = count("N_qt");
      loop.azimuth_nodes    = count("N_phi");
      N_x                   = count("N_x");

      loop.r_min = 0.0;
      loop.r_max = math::msqrt(qt2_MAX);

      const double sudakov_q2_max = j.at("NUMERICS_SUDAKOV").at("q2_MAX");
      Validate(sudakov_q2_max);

      std::cout << "MDurhamNumerics::ReadParameters: [" << XID << "]" << std::endl;
      std::cout << j.at(XID) << std::endl;
      std::cout << std::endl;
    } catch (const json::exception &error) {
      throw std::invalid_argument("MDurhamNumerics::ConfigureFromJson: Error reading " + source_file + ": " +
                                  error.what());
    }

    initialized = true;
  }
};

class MDurham : public amplitude::ProcessFamily {
 public:
  using DurhamInitialHelicity    = std::array<std::complex<double>, 4>;
  using DurhamProjectedAmp       = std::vector<DurhamInitialHelicity>;
  using DurhamTransverseMomentum = std::array<double, 2>;
  using DurhamLoopAmp            = std::function<std::vector<std::complex<double>>(
      const DurhamTransverseMomentum &, const DurhamTransverseMomentum &, const std::vector<std::complex<double>> &)>;

  MDurham(gra::LORENTZSCALAR &lts, MModelTunePtr model_tune, gra::MRandom &rng,
          std::shared_ptr<const amplitude::ProcessDefinition> definition);
  ~MDurham() {}

  // Build one immutable Durham process definition for an analytic channel
  static std::shared_ptr<const amplitude::ProcessDefinition> ProcessDefinitionFor(MDurhamMode        mode,
                                                                                  const std::string &process_name);

  // Initialize Durham and Sudakov parameters before worker process copies
  static void InitializeParameters(MProcessSetup &setup);

  // Build the tune-independent union defined by the process registry
  static std::shared_ptr<const amplitude::ProcessDefinition> ContinuumProcesses();

  double DurhamQCD(gra::LORENTZSCALAR &lts, const std::string &process);
  // Compute the status of the latest generated hard-kernel evaluation
  mg5helas::EvaluationStatus EvaluationStatus() const { return evaluation_status; }
  // Evaluate the reduced Durham loop and return pp-level channel amplitudes
  std::vector<std::complex<double>> DQtloopAmplitudes(gra::LORENTZSCALAR &lts, const DurhamProjectedAmp &Amp,
                                                      const std::vector<DurhamProjectedAmp> *ColorAmp = nullptr,
                                                      const DurhamLoopAmp                   *LoopAmp  = nullptr);
  // Evaluate the reduced Durham loop, fill lts.hamp, and return sum_h |M_h|^2
  double      DQtloop(gra::LORENTZSCALAR &lts, const DurhamProjectedAmp &Amp,
                      const std::vector<DurhamProjectedAmp> *ColorAmp = nullptr, const DurhamLoopAmp *LoopAmp = nullptr);
  bool        SampleColorFlow(gra::LORENTZSCALAR &lts);
  // Compute the two gluon PDF scales in the configured scheme
  void DScaleChoise(double qt2, double q1_2, double q2_2, double &Q1_2_scale, double &Q2_2_scale) const;

  DurhamLoopAmp Dgg2chic0(const gra::LORENTZSCALAR &lts) const;
  DurhamLoopAmp Dgg2chic1(const gra::LORENTZSCALAR &lts) const;
  DurhamLoopAmp Dgg2chic2(const gra::LORENTZSCALAR &lts) const;

  // Compute the four Durham transverse spin projectors
  void DHelicity(const DurhamTransverseMomentum &q1, const DurhamTransverseMomentum &q2,
                 std::vector<std::complex<double>> &JzP) const;

  // Contract initial helicity amplitudes with the Durham spin projectors
  std::complex<double> DHelProj(const std::vector<std::complex<double>> &A,
                                const std::vector<std::complex<double>> &JzP) const;
  // Contract a fixed initial helicity array with the Durham spin projectors
  std::complex<double> DHelProj(const DurhamInitialHelicity &A, const std::vector<std::complex<double>> &JzP) const;

  // Fill any generated finite-Nc Durham hard process through the common
  // registry
  void Dgg2Generated(gra::LORENTZSCALAR &lts, DurhamMG5Process &matrix_element, DurhamProjectedAmp &Amp,
                     std::vector<DurhamProjectedAmp> *ColorAmp = nullptr);

  void                Dgg2MMbar(const gra::LORENTZSCALAR &lts, DurhamProjectedAmp &Amp);
  double              phi_CZ(double x, double fM) const;
  std::vector<double> EvalPhi(const std::vector<double> &xval, int pdg) const;

 private:
  MModelTunePtr model_tune;
  SoftModelPtr  soft_model;
  using MesonWaveCacheKey = std::pair<int, std::vector<double>>;

 public:
  using LoopConst = math::PolarCutRule;

  // Shared immutable Durham parameters and loop numerics
  struct DurhamConfig {
    MDurhamParam                    param;
    MDurhamNumerics                 numerics;
    std::optional<LoopConst>        loop_const;
  };

 private:
  gra::MRandom                                            &rng;  // Process-local random number generator
  std::shared_ptr<const DurhamConfig>                      config;
  const MDurhamParam                                      &param;
  const MDurhamNumerics                                   &numerics;
  const LoopConst                                         &loop_const;
  const std::array<double, 2>                              q2_flavour;
  mutable std::mutex                                       meson_wave_cache_mutex;
  mutable std::map<MesonWaveCacheKey, std::vector<double>> meson_wave_cache;
  amplitude::ProcessRegistry                               generated_durham_registry;
  std::unique_ptr<DurhamMG5Process>                        generated_durham_process;
  std::optional<amplitude::Process>                        generated_durham_match;
  DurhamMG5Evaluation                                      generated_durham_evaluation;
  mg5helas::EvaluationStatus                               evaluation_status = mg5helas::EvaluationStatus::Success;

  // Compute generalized-kt jet power for the selected MadGraph jet algorithm
  int MadGraphJetPower() const;
  // Cluster partons with the selected generalized-kt MadGraph jet algorithm
  std::vector<gra::M4Vec> ClusterMadGraphJets(const std::vector<gra::M4Vec> &partons) const;
  // Apply resolved partonic jet cuts before evaluating MadGraph amplitudes
  bool PassMadGraphJetCuts(const gra::LORENTZSCALAR &lts,
                           const std::vector<int>   *final_color_representations = nullptr) const;
  // Reuse event-local jet-cut decisions throughout the screening convolution
  bool PassCachedMadGraphJetCuts(gra::LORENTZSCALAR     &lts,
                                 const std::vector<int> *final_color_representations = nullptr) const;
  // Reset projected hard amplitudes after a failed hard-process cut
  void ZeroProjectedAmplitudes(DurhamProjectedAmp &Amp, std::vector<DurhamProjectedAmp> *ColorAmp) const;
  // Resolve and cache the generated matrix_element matching the ordered final
  // state
  DurhamMG5Process &ResolveGeneratedDurhamProcess(const gra::LORENTZSCALAR &lts);

  void                       ClearCentralColorFlow(gra::LORENTZSCALAR &lts) const;
  const std::vector<double> &CachedMesonWaveFunction(const std::vector<double> &xval, int pdg) const;
};

// Construct one immutable Durham model and loop configuration
std::shared_ptr<const MDurham::DurhamConfig> ReadDurhamConfig(const MModelTune &tune);

// Compute the run owned immutable Durham configuration
std::shared_ptr<const MDurham::DurhamConfig> GetDurhamConfig(MModelCache &cache);

}  // namespace gra

#endif
