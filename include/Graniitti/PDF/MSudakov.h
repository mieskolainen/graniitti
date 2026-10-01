// Shuvaev PDF and Sudakov suppression class
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MSUDAKOV_H
#define MSUDAKOV_H

// C++
#include <array>
#include <complex>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

// LHAPDF
#include "LHAPDF/LHAPDF.h"

// Own
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MFixedStore.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/PDF/MLHAPDF.h"
#include "Graniitti/Regge/MSoftModel.h"
#include "Graniitti/Tech/MTimer.h"

// Libraries
#include "json.hpp"

using gra::aux::indices;

namespace gra {

extern std::string MODELPARAM;
std::string ResolveModelDataFile(const std::string& modelparam, const std::string& filename);

// Select the treatment of the unresolved skewed gluon below the PDF starting scale
enum class MSudakovIRMode {
  PERTURBATIVE_ONLY,
  IR_BASELINE,
  IR_RKHS
};

// Physics parameters for the matched skewed-gluon boundary condition
struct MSudakovModel {
  std::string mode_name            = "IR_BASELINE";
  double      Q0                   = 1.0;
  double      ir_logq2_length      = 1.0;
  double      rkhs_logq2_length    = 1.0;
  double      rkhs_logx_length     = 2.0;
  double      rkhs_weight_exponent = 0.5;
  double      rkhs_radius          = 0.0;
  std::vector<double> member_s_anchors;
  std::vector<double> member_x_anchors;
  std::vector<double> member_coefficients;

  // Compute the parsed infrared mode
  MSudakovIRMode Mode() const;

  // Compute the squared weighted-RKHS norm of the selected functional member
  double MemberNorm2() const;

  // Validate the matched-boundary physics parameters
  void Validate() const;

  // Read matched-boundary physics from one immutable GENERAL JSON text
  void ConfigureFromJson(const std::string &source_file,
                                  const std::string &json_text);
};

struct MSudakovNumerics {
  double       q2_MIN = 0.0;
  double       q2_MAX = 0.0;
  double       mu_MIN = 0.0;
  double       mu_MAX = 0.0;
  double       x_MIN  = 0.0;
  const double x_MAX  = 1.0 - 1E-9;

  // Numerical integral intervals
  unsigned int SudakovIntegralN = 0;
  unsigned int ShuvaevIntegralN = 0;

  // Logarithmic stepping true/false
  std::vector<bool> SUDA_log_ON;
  std::vector<bool> SHUV_log_ON;

  // Number of discrete intervals
  std::vector<unsigned int> SUDA_N;
  std::vector<unsigned int> SHUV_N;

  unsigned long config_hash              = 0;

  // CMS energy
  double sqrts = 0.0;

  bool DEBUG = false;

  // Validate Sudakov/Shuvaev numerical steering values
  void Validate() const {
    if (SudakovIntegralN < 2 || ShuvaevIntegralN < 2 || SudakovIntegralN % 2 || ShuvaevIntegralN % 2) {
      throw std::invalid_argument("MSudakovNumerics::Validate: integral intervals must be positive and even");
    }
    if (!std::isfinite(q2_MAX) || q2_MAX <= 0.0) {
      throw std::invalid_argument("MSudakovNumerics::Validate: q2_MAX must be positive");
    }
    if (SUDA_log_ON.size() != 2 || SHUV_log_ON.size() != 2 || SUDA_N.size() != 2 ||
        SHUV_N.size() != 2) {
      throw std::invalid_argument(
          "MSudakovNumerics::Validate: SUDA and SHUV log_ON/N arrays must have two entries");
    }
    if (SUDA_N[0] < 2 || SUDA_N[1] < 2 || SHUV_N[0] < 2 || SHUV_N[1] < 2) {
      throw std::invalid_argument(
          "MSudakovNumerics::Validate: SUDA and SHUV grid sizes must be at least two");
    }
    for (const auto &grid : {SUDA_N, SHUV_N}) {
      for (const unsigned int n : grid) {
        if (n >= static_cast<unsigned int>(std::numeric_limits<int>::max())) {
          throw std::invalid_argument("MSudakovNumerics::Validate: grid size exceeds index range");
        }
      }
    }
  }

  // Read Sudakov numerical steering from one immutable NUMERICS JSON text
  void ConfigureFromJson(const std::string &source_file,
                                  const std::string &json_text) {
    using json = nlohmann::json;
    json j;

    try {
      j = json::parse(json_text);

      // JSON block identifier
      const std::string XID = "NUMERICS_SUDAKOV";

      q2_MAX = j.at(XID).at("q2_MAX");
      // mu_MAX = j.at(XID).at("mu_MAX"); // Post-Setup
      // x_MIN = j.at(XID).at("x_MIN"); // Post-Setup
      // x_MAX = j.at(XID).at("x_MAX"); // Set above

      // Logarithmic stepping true/false [assign needed for <cast>]
      SUDA_log_ON = j.at(XID).at("SUDA").at("log_ON").get<std::vector<bool>>();
      SHUV_log_ON = j.at(XID).at("SHUV").at("log_ON").get<std::vector<bool>>();

      // Validate integer grid counts before conversion can wrap or truncate them
      const auto read_grid = [](const json &grid) {
        if (!grid.is_array() || grid.size() != 2) {
          throw std::invalid_argument("MSudakovNumerics: grid N must contain two integers");
        }
        std::vector<unsigned int> counts;
        counts.reserve(2);
        for (const auto &n : grid) {
          if (!n.is_number_integer() || n < 2 || n >= std::numeric_limits<int>::max()) {
            throw std::invalid_argument("MSudakovNumerics: grid N is outside the integer index range");
          }
          counts.push_back(n.get<unsigned int>());
        }
        return counts;
      };
      SUDA_N = read_grid(j.at(XID).at("SUDA").at("N"));
      SHUV_N = read_grid(j.at(XID).at("SHUV").at("N"));

      config_hash             = gra::aux::djb2hash(j.at(XID).dump());

      const auto intervals =
          read_grid(json::array({j.at(XID).at("SUDA").at("N_int"), j.at(XID).at("SHUV").at("N_int")}));
      SudakovIntegralN = intervals[0];
      ShuvaevIntegralN = intervals[1];

      DEBUG = j.at(XID).at("DEBUG");
      Validate();

      std::cout << "MSudakovNumerics::ReadParameters: [NUMERICS_SUDAKOV]" << std::endl;
      std::cout << j.at(XID) << std::endl;
      std::cout << std::endl;

    } catch (const nlohmann::json::exception &error) {
      throw std::invalid_argument(
          "MSudakovNumerics::ConfigureFromJson: Error reading " +
          source_file + ": " + error.what());
    }
  }
};


// Interpolation container
class IArray2D {
 public:
  // Hold one validated fixed coordinate on the second interpolation axis
  struct PreparedSecondCoordinate {
    const IArray2D *owner = nullptr;
    int             index = 0;
    double          fraction = 0.0;
  };

  IArray2D(){};
  ~IArray2D(){};

  std::string  name[2];
  double       MIN[2]   = {0};
  double       MAX[2]   = {0};
  unsigned int N[2]     = {0};
  double       STEP[2]  = {0};
  bool         islog[2] = {false};

  // Setup discretization
  void Set(unsigned int VAR, std::string _name, double _min, double _max,
           unsigned int _N, bool _logarithmic) {
    // Out of index
    if (VAR > 1) {
      throw std::invalid_argument("IArray2D::Set: Error: VAR = " + std::to_string(VAR) + " > 1");
    }

    // Require at least two intervals with representable signed interpolation indices
    if (_N < 2 || _N >= static_cast<unsigned int>(std::numeric_limits<int>::max())) {
      throw std::invalid_argument("IArray2D::Set: grid N is outside the integer index range");
    }

    // Require finite, ordered interpolation bounds
    if (!std::isfinite(_min) || !std::isfinite(_max) || _min >= _max) {
      throw std::invalid_argument("IArray2D::Set: Error: Variable " + _name + " MIN = " +
                                  std::to_string(_min) + " >= MAX = " + std::to_string(_max));
    }

    // Require a positive lower bound for logarithmic spacing
    if (_logarithmic && !(_min > 0.0)) {
      throw std::invalid_argument(
          "IArray2D::Set: Error: Variable " + _name +
          " is using logarithmic stepping with boundary MIN = " + std::to_string(_min));
    }
    // Check that transformed grid coordinates retain a finite resolved spacing
    const double lower = _logarithmic ? std::log(_min) : _min;
    const double upper = _logarithmic ? std::log(_max) : _max;
    const double step = (upper - lower) / _N;
    if (!std::isfinite(step) || !(step > 0.0) ||
        !(lower + step > lower) || !(upper - step < upper)) {
      throw std::invalid_argument("IArray2D::Set: unresolved grid spacing for " + _name);
    }
    islog[VAR] = _logarithmic;
    MIN[VAR] = lower;
    MAX[VAR] = upper;
    N[VAR]   = _N;

    STEP[VAR] = step;
    name[VAR] = _name;
  }

  // Call this last
  void InitArray() {
    // Note N+1 !
    F = std::vector<std::vector<std::array<double, 4>>>(
        N[0] + 1, std::vector<std::array<double, 4>>(N[1] + 1));
    first_axis_physical.resize(N[0] + 1);
    for (const auto &i : indices(first_axis_physical)) {
      const double coordinate = MIN[0] + i * STEP[0];
      first_axis_physical[i] = islog[0] ? std::exp(coordinate) : coordinate;
    }
  }

  // N+1 per dimension!
  std::vector<std::vector<std::array<double, 4>>> F;
  std::vector<double> first_axis_physical;

  // CMS energy and PDF used to generate the stored values
  double sqrts = 0.0;
  std::string pdf;
  std::string quantity;

  std::string GetHashString() const {
    std::string str = std::to_string(islog[0]) + std::to_string(islog[1]) + std::to_string(MIN[0]) +
                      std::to_string(MIN[1]) + std::to_string(MAX[0]) + std::to_string(MAX[1]) +
                      std::to_string(N[0]) + std::to_string(N[1]) + std::to_string(sqrts);
    return str;
  }

  // Describe the stored coordinates without assuming a particular grid or PDF
  nlohmann::json Metadata() const {
    return {{"axes", name}, {"log", islog}, {"sqrts", sqrts}, {"pdf", pdf}, {"quantity", quantity}};
  }

  bool                      WriteArray(const std::string &filename, bool overwrite) const;
  bool                      ReadArray(const std::string &filename);
  std::pair<double, double> Interpolate2D(double A, double B, double* second = nullptr) const;
  // Prepare a repeatedly used physical coordinate on the second axis
  PreparedSecondCoordinate PrepareSecondCoordinate(double B) const;
  // Interpolate with a pretransformed first coordinate and prepared second coordinate
  std::pair<double, double> InterpolatePrepared(double input_a, double coordinate_a,
                                                const PreparedSecondCoordinate& prepared_b,
                                                double*                         second = nullptr) const;
  // Interpolate prepared coordinates already known to lie inside the grid
  std::pair<double, double> InterpolatePreparedInDomain(double input_a, double coordinate_a,
                                                        const PreparedSecondCoordinate& prepared_b,
                                                        double*                         second = nullptr) const;
};

// Sudakov suppression and skewed pdf
//
// Take care when copying this class - it caches interpolation arrays
class MSudakov {
 public:
  // Hold event-local fixed-x and fixed-mu data for fast Durham flux evaluation
  struct PreparedFlux {
    // Compute f_g/Q2 with the same matched perturbative and infrared model
    double OverQ2(double q2);

   private:
    struct InterpolationField {
      double value_1 = 0.0;
      double value_2 = 0.0;
      double slope_1 = 0.0;
      double slope_2 = 0.0;
    };

    struct PerturbativeCell {
      int index = -1;
      double coordinate_1 = 0.0;
      double inverse_step = 0.0;
      InterpolationField shuvaev;
      InterpolationField sudakov;
    };

    friend class MSudakov;
    // Prepare the infrared continuation parameters and analytic zero limit
    void PrepareInfrared();
    // Compute the infrared continuation contribution to f_g/Q2
    double InfraredOverQ2(double q2) const;
    // Populate one direct-mapped interpolation cell for both perturbative fields
    PerturbativeCell &PreparePerturbativeCell(int index);
    // Evaluate one cached cubic-Hermite field and its physical derivative
    std::pair<double, double> EvaluatePerturbativeCell(
        const PerturbativeCell &cell, const InterpolationField &field,
        double input_q2, double coordinate_q2) const;
    const MSudakov *owner = nullptr;
    IArray2D::PreparedSecondCoordinate shuvaev_x;
    IArray2D::PreparedSecondCoordinate sudakov_mu;
    double x = 0.0;
    double mu = 0.0;
    MSudakovIRMode infrared_mode = MSudakovIRMode::PERTURBATIVE_ONLY;
    double q2_min = 0.0;
    double q2_max = 0.0;
    double infrared_lambda = 0.0;
    double infrared_c = 0.0;
    double infrared_d = 0.0;
    double boundary_integrated = 0.0;
    double boundary_flux = 0.0;
    double boundary_over_q2_min = 0.0;
    double zero_limit = 0.0;
    bool aligned_first_axis = false;
    bool valid = false;
    bool                               infrared_ready       = false;
    std::array<PerturbativeCell, 32> perturbative_cache = {};
  };

  MSudakov();
  ~MSudakov();

  void Init(double sqrts_in, const std::string &PDFSET,
            const SoftModelPtr &soft_model,
            const std::string &numerics_source,
            const std::string &numerics_json,
            bool init_arrays = true,
            std::shared_ptr<const LHAPDF::PDF> pdf = nullptr);
  void InitLHAPDF(const std::string &PDFSET,
                  std::shared_ptr<const LHAPDF::PDF> pdf = nullptr);

  double fg_xQ2Mu(double x, double q2, double mu) const;
  // Prepare one event-local fixed-x and fixed-mu Durham flux evaluator
  PreparedFlux PrepareFlux(double x, double mu) const;
  // Compute the matched integrated amplitude-level gluon H_g sqrt(T)
  double IntegratedFlux_xQ2Mu(double x, double q2, double mu) const;
  // Compute f_g / Q2 with the analytic color-neutral limit at Q2 = 0
  double fg_xQ2MuOverQ2(double x, double q2, double mu) const;
  // Compute alpha_s times f_g with alpha_s frozen at Q0 outside the RKHS variation
  double AlphaSFlux_xQ2Mu(double x, double q2, double mu, double alpha_scale2) const;
  // Compose the signed logarithmic derivative entering the Durham gluon flux
  double FluxDerivative(double q2, double Hg, double dHg, double Tg, double dTg) const;
  double AlphaS_Q2(double q2) const;
  double NumFlavor(double q2) const;
  double xg_xQ2(double x, double Q2) const;
  double PerturbativeShuvaevGluon_xQ2(double x, double q2) const;
  double PerturbativeShuvaevGluonDerivative_xQ2(double x, double q2) const;
  double PerturbativeSqrtSudakov_Q2Mu(double q2, double mu) const;
  // Compute the Sudakov soft-emission cutoff Delta = kt / mu
  double SudakovDelta(double kt2, double mu) const;
  // Compute the direct radiation veto for numerical convergence checks
  std::pair<double, double> Sudakov_T(double qt2, double mu) const;
  // Compute whether the hard scale has a defined radiation interval
  bool SupportsMu(double mu) const { return std::isfinite(mu) && mu > 0.0 && mu <= Numerics.mu_MAX; }
  // Compute the physical perturbative matching scale squared
  double GetQ2Min() const { return Numerics.q2_MIN; }
  // Access the immutable SOFT model that owns the UGD physics snapshot
  const SoftModelPtr &SoftModelHandle() const noexcept { return soft_model; }
  // Compute the physical perturbative matching scale
  double GetQMin() const { return Numerics.mu_MIN; }
  // Compute the upper q2 interpolation boundary
  double GetQ2Max() const { return Numerics.q2_MAX; }
  // Compute the lower factorization-scale interpolation boundary
  double GetMuMin() const { return Numerics.mu_MIN; }
  // Compute the upper factorization-scale interpolation boundary
  double GetMuMax() const { return Numerics.mu_MAX; }
  // Compute the configured infrared treatment
  MSudakovIRMode GetIRMode() const { return Model.Mode(); }
  // Compute true when q2 belongs to the resolved perturbative domain
  bool IsPerturbativeQ2(double q2) const { return q2 >= Numerics.q2_MIN; }
  // Compute the selected functional member weighted-RKHS norm squared
  double GetMemberNorm2() const { return Model.MemberNorm2(); }
  // Compute the admissible weighted-RKHS radius
  double GetRKHSRadius() const { return Model.rkhs_radius; }
  // Compute the representer anchors in logarithmic infrared distance
  const std::vector<double> &GetMemberSAnchors() const {
    return Model.member_s_anchors;
  }
  // Compute the representer anchors in Bjorken x
  const std::vector<double> &GetMemberXAnchors() const {
    return Model.member_x_anchors;
  }
  // Compute an immutable clone carrying one globally correlated functional member
  std::shared_ptr<const MSudakov> WithMemberCoefficients(
      const std::vector<double> &coefficients) const;
  void   TestPDF() const;

  bool initialized = false;

 private:
  struct MatchedFlux {
    double integrated = 0.0;
    double flux       = 0.0;
  };

  void   InitArrays();
  std::pair<double, double> Shuvaev_H(double q2, double x);
  std::pair<double, double> ShuvaevTransform(double q2, double x);
  std::pair<double, double> SudakovRadiator(double qt2, double mu);
  double xg_xQ2_raw(double x, double q2) const;
  double diff_xg_xQ2_wrt_Q2(double x, double q2) const;
  MatchedFlux PerturbativeFlux(double x, double q2, double mu) const;
  double PerturbativeFluxLogDerivative(double x, double q2, double mu) const;
  MatchedFlux InfraredFlux(double x, double q2, double mu) const;
  double FunctionalKernel(double s, double anchor, double x,
                          double x_anchor) const;
  double FunctionalKernelIntegral(double s, double anchor, double x,
                                  double x_anchor) const;
  double AP_gg(double delta) const;
  double AP_qg(double delta, double qt2) const;
  void   CalculateArray(IArray2D &arr, std::pair<double, double> (MSudakov::*f)(double, double));
  void   LoadOrBuildArray(IArray2D &arr, const std::string &filename,
                          std::pair<double, double> (MSudakov::*f)(double, double));

  std::string  PDFSETNAME;
  SoftModelPtr soft_model;
  std::shared_ptr<const LHAPDF::PDF> PdfPtr = nullptr;
  double pdf_q2_MIN = 0.0;
  double charm_mass = 0.0;
  double bottom_mass = 0.0;

  IArray2D veto;  // Sudakov veto
  IArray2D spdf;  // Shuvaev pdf

  MSudakovModel    Model;     // Matched-boundary physics
  MSudakovNumerics Numerics;  // Numerics
};

// Process-wide cache for initialized read-only Sudakov tables
class MSudakovStore {
 public:
  using SudakovPtr = std::shared_ptr<const MSudakov>;

  // Construct a standalone store with its own LHAPDF store
  MSudakovStore();

  // Construct a run owned store sharing one LHAPDF store
  explicit MSudakovStore(MLHAPDFStore &pdf_store);

  // Compute one shared read-only Sudakov/Shuvaev table set
  SudakovPtr GetSudakov(double sqrts, const std::string &pdfset,
                        const SoftModelPtr &soft_model);

 private:
  struct Key {
    SoftModelPtr soft_model;
    double sqrts = 0.0;
    std::string pdf_identity;

    // Provide strict weak ordering for cache lookup
    bool operator<(const Key &other) const;
  };

  // Build the full cache key from runtime energy and current numerics
  Key MakeKey(double sqrts, const std::string &pdf_identity,
              const SoftModelPtr &soft_model) const;

  // Construct one fully initialized Sudakov object
  SudakovPtr LoadSudakov(double sqrts, const std::string &pdfset,
                         const SoftModelPtr &soft_model,
                         const MLHAPDFStore::PDFPtr &pdf) const;

  MFixedStore<Key, MSudakov> store;
  std::unique_ptr<MLHAPDFStore> owned_pdf;
  MLHAPDFStore *pdf_store = nullptr;
};

}  // namespace gra

#endif
