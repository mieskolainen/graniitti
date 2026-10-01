// Coupled channel proton eikonal amplitudes and diffraction
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MEIKONAL_H
#define MEIKONAL_H

// C++
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Own
#include "Graniitti/Eikonal/MEikonalMatrix.h"
#include "Graniitti/Eikonal/MElasticCNI.h"
#include "Graniitti/MModelTune.h"
#include "Graniitti/Kinematics/M4Vec.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Math/MPolarQuadrature.h"
#include "Graniitti/Particle/MForm.h"
#include "Graniitti/Particle/MPDG.h"
#include "Graniitti/Regge/MSoftModel.h"
#include "Graniitti/Tech/MAux.h"

// External
#include "json.hpp"
#include "rang.hpp"

namespace gra {

// Inclusive physical Good-Walker final-state classes
enum class GWFinalClass {
  Elastic,
  SingleDissociation1,
  SingleDissociation2,
  DoubleDissociation
};

// One exclusive physical Good-Walker final-state pair
struct GWExclusiveState {
  std::size_t f1 = 0;
  std::size_t f2 = 0;
};

// Event-by-event loop 2D-integral
using LoopNumerics = math::PolarParam;

// Numerical integration parameters
// See NUMERICS.json
struct MEikonalNumerics {

  // S-matrix contractivity validation
  bool strict_unitarity = false;
  double unitarity_tolerance = 0.0;

  // Pre-computation arrays

  // kt^2-space lookup table for the screening amplitude
  // (MinKT2 is a numerical grid floor, not a physical cutoff)
  double MinKT2 = 0.0;
  double MaxKT2 = 0.0;
  unsigned int NumberKT2 = 0;
  bool logKT2 = false;

  // b_t-space
  double MinBT = 0.0;
  double MaxBT = 0.0;
  unsigned int NumberBT = 0;
  bool logBT = false;

  // Fourier-Bessel integral
  double FBIntegralMinKT = 0.0;
  double FBIntegralMaxKT = 0.0;
  unsigned int FBIntegralN = 0;

  ElasticCNINumerics CNI;

  LoopNumerics LOOP;

  int extra_NumberLoopKT = 0;
  int extra_NumberLoopPHI = 0;

  // Set independent signed offsets relative to the tuned loop discretization
  void SetLoopDiscretization(int extra_kt, int extra_phi) {
    extra_NumberLoopKT = extra_kt;
    extra_NumberLoopPHI = extra_phi;
  }

  // Unique hash
  std::string GetHashString() {
    // Event-by-event loop integral parameters not included in the hash
    // log(x) taken for numerical representation reasons
    std::string str =
        std::to_string(strict_unitarity) + std::to_string(unitarity_tolerance) +
        std::to_string(std::log10(MinKT2 + 1E-32)) + std::to_string(MaxKT2) +
        std::to_string(NumberKT2) + std::to_string(logKT2) +
        std::to_string(std::log10(MinBT + 1E-32)) + std::to_string(MaxBT) +
        std::to_string(NumberBT) + std::to_string(logBT) +
        std::to_string(std::log10(FBIntegralMinKT + 1E-32)) +
        std::to_string(FBIntegralMaxKT) + std::to_string(FBIntegralN);
    return str;
  }

  // Read parameters from one numerical steering file
  void ReadParameters(const std::string &inputfile) {
    ReadParameters(inputfile, gra::aux::GetInputData(inputfile));
  }

  // Read parameters from one immutable numerical steering snapshot
  void ReadParameters(const std::string &source, const std::string &data) {

    using json = nlohmann::json;

    json j;

    try {
      j = json::parse(data);

      // JSON block identifier
      const std::string XID = "NUMERICS_EIKONAL";

      strict_unitarity = j.at(XID).at("strict_unitarity");
      unitarity_tolerance = j.at(XID).at("unitarity_tolerance");

      MinKT2 = j.at(XID).at("MinKT2");
      MaxKT2 = j.at(XID).at("MaxKT2");
      NumberKT2 = j.at(XID).at("NumberKT2");
      logKT2 = j.at(XID).at("logKT2");

      MinBT = j.at(XID).at("MinBT");
      MinBT /= PDG::GeV2fm;
      MaxBT = j.at(XID).at("MaxBT");
      MaxBT /= PDG::GeV2fm;

      NumberBT = j.at(XID).at("NumberBT");
      logBT = j.at(XID).at("logBT");

      FBIntegralMinKT = j.at(XID).at("FBIntegralMinKT");
      FBIntegralMaxKT = j.at(XID).at("FBIntegralMaxKT");
      FBIntegralN = j.at(XID).at("FBIntegralN");

      CNI.NumberKT2 = j.at(XID).at("CNI").at("NumberKT2");
      CNI.logKT2 = j.at(XID).at("CNI").at("logKT2");
      CNI.FBIntegralMaxKT = j.at(XID).at("CNI").at("FBIntegralMaxKT");
      CNI.FBIntegralN = j.at(XID).at("CNI").at("FBIntegralN");
      CNI.tail_max_qb = j.at(XID).at("CNI").at("tail_max_qb");
      CNI.TailIntegralN = j.at(XID).at("CNI").at("TailIntegralN");
      CNI.tail_order = j.at(XID).at("CNI").at("tail_order");
      CNI.tail_split_qb = j.at(XID).at("CNI").at("tail_split_qb");
      CNI.interp_rel_tol = j.at(XID).at("CNI").at("interp_rel_tol");

      // Event-by-event loop integral parameters
      LOOP.radial_integrator =
          j.at(XID).at("LOOP_INTEGRAL").at("kT_integrator");
      LOOP.azimuth_integrator =
          j.at(XID).at("LOOP_INTEGRAL").at("phi_integrator");
      LOOP.radial_map = j.at(XID).at("LOOP_INTEGRAL").at("log_kT")
                            ? math::RadialMap::Log
                            : math::RadialMap::Linear;
      LOOP.r_min = j.at(XID).at("LOOP_INTEGRAL").at("MinLoopKT");
      LOOP.r_max = j.at(XID).at("LOOP_INTEGRAL").at("MaxLoopKT");
      LOOP.radial_intervals =
          j.at(XID).at("LOOP_INTEGRAL").at("NumberLoopKT");
      LOOP.azimuth_nodes =
          j.at(XID).at("LOOP_INTEGRAL").at("NumberLoopPHI");

      std::cout << "MEikonalNumerics::ReadParameters: [" << XID << "]"
                << std::endl;
      std::cout << j.at(XID) << std::endl;
      std::cout << std::endl;

    } catch (const std::exception &error) {
      throw std::invalid_argument(
          "MEikonalNumerics::ReadParameters: Error reading " + source + ": " +
          error.what());
    }

    unsigned int multiplier = 0;
    try {
      multiplier = math::PolarIntervalMultiple(LOOP);
    } catch (const std::invalid_argument &) {
      throw std::invalid_argument(
          "MEikonalNumerics::ReadParameters: Unknown LOOP.kT_integrator = " +
          LOOP.radial_integrator);
    }
    if (LOOP.azimuth_integrator != "Trap") {
      throw std::invalid_argument(
          "MEikonalNumerics::ReadParameters: Unknown LOOP.phi_integrator = " +
          LOOP.azimuth_integrator);
    }

    // User re-adjust
    auto AdjustLoopCount = [&](unsigned int base, int extra, int step,
                               const std::string &name) {
      const long long adjusted =
          static_cast<long long>(base) + static_cast<long long>(step) * extra;
      if (adjusted <= 0 ||
          adjusted > static_cast<long long>(
                         std::numeric_limits<unsigned int>::max())) {
        throw std::invalid_argument(
            "MEikonalNumerics::ReadParameters: " + name +
            " adjustment gives invalid loop count " + std::to_string(adjusted));
      }
      return static_cast<unsigned int>(adjusted);
    };
    LOOP.radial_intervals =
        AdjustLoopCount(LOOP.radial_intervals, extra_NumberLoopKT, multiplier,
                        "NumberLoopKT");
    LOOP.azimuth_nodes = AdjustLoopCount(
        LOOP.azimuth_nodes, extra_NumberLoopPHI, 1, "NumberLoopPHI");

    auto ValidateLoopCount = [&](unsigned int count, const std::string &name) {
      if (count % multiplier != 0) {
        throw std::invalid_argument("MEikonalNumerics::ReadParameters: LOOP." +
                                    name + " = " + std::to_string(count) +
                                    " is not compatible with kT_integrator " +
                                    LOOP.radial_integrator);
      }
    };
    if (LOOP.radial_integrator != "GL") {
      ValidateLoopCount(LOOP.radial_intervals, "NumberLoopKT");
    }

    if (LOOP.r_min < 0.0) {
      throw std::invalid_argument(
          "MEikonalNumerics::ReadParameters: LOOP.MinLoopKT must be >= 0");
    }
    if (LOOP.r_max <= LOOP.r_min) {
      throw std::invalid_argument("MEikonalNumerics::ReadParameters: "
                                  "LOOP.MinLoopKT < LOOP.MaxLoopKT required");
    }
    if (MaxKT2 <= 0.0) {
      throw std::invalid_argument(
          "MEikonalNumerics::ReadParameters: MaxKT2 must be > 0");
    }
    if (NumberBT < 2 || NumberBT % 2 != 0) {
      throw std::invalid_argument(
          "MEikonalNumerics::ReadParameters: NumberBT must be even and at "
          "least two");
    }
    if (!std::isfinite(unitarity_tolerance) || unitarity_tolerance <= 0.0) {
      throw std::invalid_argument(
          "MEikonalNumerics::ReadParameters: unitarity_tolerance must be "
          "positive");
    }
    const double loop_max_kt2 = gra::math::pow2(LOOP.r_max);
    if (loop_max_kt2 > MaxKT2 * (1.0 + 1e-12)) {
      throw std::invalid_argument("MEikonalNumerics::ReadParameters: "
                                  "LOOP.MaxLoopKT^2 must be <= MaxKT2");
    }
    if (LOOP.radial_map == math::RadialMap::Log && LOOP.r_min <= 0.0) {
      throw std::invalid_argument(
          "MEikonalNumerics::ReadParameters: LOOP.log_kT = true requires "
          "LOOP.MinLoopKT > 0");
    }
    if (CNI.NumberKT2 < 3 || !std::isfinite(CNI.FBIntegralMaxKT) ||
        CNI.FBIntegralMaxKT <= 0.0 || CNI.FBIntegralN < 2 ||
        !std::isfinite(CNI.tail_max_qb) || CNI.tail_max_qb <= 0.0 ||
        CNI.TailIntegralN < 2 || CNI.TailIntegralN % 2 != 0 ||
        CNI.tail_order < 1 || static_cast<unsigned int>(CNI.tail_order) > CNI.TailIntegralN ||
        !std::isfinite(CNI.tail_split_qb) || CNI.tail_split_qb <= 0.0 ||
        !std::isfinite(CNI.interp_rel_tol) || CNI.interp_rel_tol <= 0.0) {
      throw std::invalid_argument(
          "MEikonalNumerics::ReadParameters: invalid CNI configuration");
    }

    std::cout << "After steering card override:" << std::endl;
    std::cout << "- LOOP.kT_integrator  = " << LOOP.radial_integrator
              << std::endl;
    std::cout << "- LOOP.phi_integrator = " << LOOP.azimuth_integrator
              << std::endl;
    std::cout << "- LOOP.log_kT       = "
              << (LOOP.radial_map == math::RadialMap::Log) << std::endl;
    std::cout << "- LOOP.NumberLoopKT  = " << LOOP.radial_intervals
              << std::endl;
    std::cout << "- LOOP.NumberLoopPHI = " << LOOP.azimuth_nodes
              << std::endl;
    std::cout << std::endl;
  }
};

class MEikonal {
public:
  // One exclusive Good Walker output with its precontracted loop operators
  struct LoopGoodWalkerChannel {
    GWExclusiveState final_state;
    double born_coefficient = 0.0;
    bool spin_scalar = true;
    std::vector<std::complex<double>> scalar_screening_weight;
    std::vector<ProtonHelicityMatrix> screening_weight;
  };

  // Event-by-event screening loop constants evaluated once per eikonal setup
  struct LoopConst {
    // Retain the quadrature definition used to construct the cached nodes
    math::PolarParam quadrature;
    bool initialized = false;
    bool log_kT = false;
    bool physical_screening_is_scalar = true;
    bool pair_screening_cache_ready = false;
    bool pair_screening_spin_scalar = false;
    bool pair_screening_pair_diagonal = false;

    double StepKT = 0.0;
    double MinLogKT = 0.0;
    double StepPhi = 0.0;
    double quadrature_scale = 0.0;
    double s = 0.0;

    std::complex<double> norm = 0.0;

    std::vector<double> kt;
    std::vector<double> kt2;
    std::vector<double> jac;

    MMatrix<double> W;
    MMatrix<double> measure_weight;
    MMatrix<double> kt_x;
    MMatrix<double> kt_y;
    MMatrix<std::complex<double>> node_weight;
    // Physical proton screening amplitude for each radial node
    std::vector<std::complex<double>> physical_screening_amplitude;
    // Complete oriented pair-space screening operators for each radial node
    std::vector<MEikonalMatrix::PairSpinBank> pair_screening_spin;
    // Cached scalar pair-space diagonal for each radial node
    std::vector<std::vector<std::complex<double>>> pair_screening_diagonal;
    // Complete azimuthal harmonic weights for each integration node
    std::vector<std::array<std::complex<double>, 5>>
        pair_screening_harmonic_weight;
    // Complete physical proton screening weight for each integration node
    MMatrix<std::complex<double>> physical_screening_weight;
    // Complete pair-helicity screening operator for each integration node
    std::vector<ProtonHelicityMatrix> physical_screening_helicity_weight;
    // Physical low-mass Good Walker classes
    std::array<std::vector<LoopGoodWalkerChannel>, 4> good_walker_channels;
  };

  MEikonal();
  explicit MEikonal(MModelTunePtr tune);
  ~MEikonal();

  // Access the immutable SOFT model used by this eikonal object
  const SoftModelPtr &SoftModelHandle() const noexcept { return soft_model; }

  // Access the immutable complete tune used by this eikonal object
  const MModelTunePtr &ModelTuneHandle() const noexcept { return model_tune; }

  // Build the Good Walker eikonal and its b and qT space lookup tables
  void S3Constructor(double s_in,
                     const std::vector<gra::MParticle> &initialstate_in,
                     bool onlydensity = false, int NumberBT = 0,
                     int NumberKT2 = 0, double max_kt2 = 0.0);

  // Initialization already done
  bool IsInitialized() const {
    if (S3INIT == true) {
      return true;
    }
    return false;
  }

  // Compute the Mandelstam s value of the current eikonal initialization
  double InitializedMandelstamS() const noexcept { return s; }

  // Compute the ordered beam metadata of the current eikonal initialization
  const std::vector<MParticle> &InitialState() const noexcept {
    return INITIALSTATE;
  }

  // Compute the optical theorem total, elastic and inelastic cross sections
  void GetTotXS(double &tot, double &el, double &in) const;
  const MMatrix<double> &GetInclusiveDiffXS() const { return sigma_diff; }
  const MMatrix<double> &GetExclusiveDiffXS() const {
    return sigma_diff_exclusive;
  }

  // Compute the minimum effective cut-Pomeron opacity from initialization
  double MinOpacityEigenvalue() const noexcept {
    return min_opacity_eigenvalue;
  }

  const MMatrix<double> &GetMixingMatrix() const { return U; }
  std::size_t GetChannelCount() const { return U.size_row(); }

  // Extract the six independent proton helicity screening operators
  MEikonalMatrix::PairHelicityBank PairScreeningHelicityBank(double kt2) const;
  // Compute all spin transitions of the radial Good Walker screening operator
  MEikonalMatrix::PairSpinBank PairScreeningSpinBank(double kt2) const;

  // Compute the coupled-channel matrix tables for physical amplitude projections
  const MEikonalMatrix &GetMatrixRuntime() const {
    if (!matrix_runtime) {
      throw std::logic_error(
          "MEikonal::GetMatrixRuntime: matrix eikonal is not initialized");
    }
    return *matrix_runtime;
  }
  // Build coherent Coulomb nuclear interference in the elastic proton channel
  void InitializeElasticCNI(double min_abs_t, double max_abs_t);
  // Compute whether the physical elastic CNI runtime is initialized
  bool HasElasticCNI() const noexcept {
    return static_cast<bool>(elastic_cni_runtime);
  }
  // Compute the coherent strong plus electromagnetic proton helicity matrix
  ProtonHelicityMatrix PhysicalElasticHelicityMatrix(const M4Vec &p1_in,
                                                     const M4Vec &p2_in,
                                                     const M4Vec &p1_out,
                                                     const M4Vec &p2_out) const;
  // Rebuild the event-by-event screening loop quadrature constants
  void InitLoopWeightMatrix();
  // Compute the event-by-event screening loop quadrature weights
  const MMatrix<double> &GetLoopWeightMatrix() const;
  // Compute the event-by-event screening loop quadrature constants
  const LoopConst &GetLoopConst(double current_s) const;
  // Compute the event-by-event screening loop radial quadrature coefficient
  double GetLoopQuadratureCoefficient() const;

  // Compute one Regge pole Born amplitude between diffractive eigenstates
  static std::complex<double> SingleAmpElastic(const SoftModelPtr &model,
                                               double s, double t,
                                               SoftExchangeId exchange,
                                               std::size_t _i, std::size_t _k);
  // Compute crossing even and odd eigenchannel profiles Xi_ik(s,b)
  std::pair<std::complex<double>, std::complex<double>>
  S3DensityProfiles(double bt, std::size_t _i, std::size_t _k) const;
  // Form the beam eikonal with the physical C-odd crossing sign
  std::complex<double> S3Density(double bt, std::size_t _i,
                                 std::size_t _k) const;
  // Transform one unitarized eigenchannel profile from b to qT space
  std::complex<double> S3Screening(double kt2, std::size_t _i,
                                   std::size_t _k) const;
  // Project the spin averaged screening amplitude onto physical protons
  std::complex<double> S3PhysicalScreeningAmp(double kt2) const;
  // Project the full amplitude onto one exclusive Good Walker final state
  std::complex<double> S3ExclusiveAmp(double kt2, std::size_t f1,
                                      std::size_t f2) const;
  // Project selected screening exchanges between Good Walker pair states
  std::complex<double> S3ExclusiveScreeningAmp(double kt2, std::size_t f1,
                                               std::size_t f2,
                                               std::size_t i1 = 0,
                                               std::size_t i2 = 0) const;
  // Enumerate all exclusive states in one Good Walker diffraction class
  std::vector<GWExclusiveState>
  ExclusiveFinalStates(GWFinalClass final_class) const;
  // Sum the incoherent spin averaged norm of one diffraction class
  double S3DiffractiveAmpSquared(double kt2, GWFinalClass final_class) const;

  // Get random bt and the number of cut Pomerons
  template <typename T>
  void S3GetRandomCutsBt(unsigned int &m, double &bt, T &rng) {
    if (!S3INIT || cut_b.size() < 2) {
      throw std::invalid_argument(
          "MEikonal::S3GetRandomCutsBt: cut-Pomeron grid is not initialized");
    }
    if (P_array.size() < MCUT || !(MAX_P_array > 0.0) ||
        !std::isfinite(MAX_P_array)) {
      throw std::invalid_argument(
          "MEikonal::S3GetRandomCutsBt: invalid cut-Pomeron probability table");
    }
    for (std::size_t cut = 1; cut < MCUT; ++cut) {
      if (P_array[cut].size() != cut_b.size()) {
        throw std::invalid_argument(
            "MEikonal::S3GetRandomCutsBt: cut-Pomeron grid size mismatch");
      }
    }
    // Construct distributions from the current initialized object
    std::uniform_real_distribution<double> flat(0, 1);
    std::uniform_int_distribution<std::size_t> randbt(0, cut_b.size() - 1);
    std::uniform_int_distribution<unsigned int> randm(1, MCUT - 1); // 1,2,...

    // Acceptance-Rejection (this could be made more efficient)
    std::size_t bt_i = 0;
    while (true) {

      // Draw random impact parameter bt (index) from Uniform
      bt_i = randbt(rng);

      const int m_cut = randm(rng);
      if (flat(rng) < P_array[m_cut][bt_i] / MAX_P_array) {

        // Number of cut Pomerons
        m = m_cut;
        break;
      }
    }
    // Use the immutable nodes belonging to the sampled probability table
    bt = cut_b[bt_i];
  }

  // Get random number of cut Pomerons (bt-space already integrated over)
  template <typename T> unsigned int S3GetRandomCuts(T &rng) {
    if (!S3INIT || P_cut.size() < 2 || !(MAX_P_cut > 0.0) ||
        !std::isfinite(MAX_P_cut)) {
      throw std::invalid_argument("MEikonal::S3GetRandomCuts: cut-Pomeron "
                                  "probabilities are not initialized");
    }

    // Construct distributions from the current probability table
    std::uniform_real_distribution<double> flat(0, 1);
    std::uniform_int_distribution<unsigned int> RANDI(
        1, static_cast<unsigned int>(P_cut.size() - 1));

    // Acceptance-Rejection
    while (true) {
      const unsigned int m = RANDI(rng);
      if (flat(rng) < P_cut[m] / MAX_P_cut) {
        return m; // m cut Pomerons
      }
    }
  }

  // Parameters
  MEikonalNumerics Numerics;

private:
  static const unsigned int MCUT = 25; // Maximum number of cut Pomerons

  // S3 soft survival initialized
  bool S3INIT = false;
  // Momentum-space screening amplitudes initialized
  bool screening_amplitude_initialized = false;

  // Coupled-channel runtime tables owned by this eikonal initialization
  std::shared_ptr<const MEikonalMatrix> matrix_runtime;

  // Physical elastic full-spin electromagnetic interference tables
  std::shared_ptr<const MElasticCNI> elastic_cni_runtime;

  // Immutable SOFT model shared with every table built by this object
  MModelTunePtr model_tune;
  SoftModelPtr soft_model;

  void S3CalcEikonalXS(bool amplitude_diagnostics);
  void S3InitCutPomerons();
  // void S3CalcEikonalXSTest();
  // Compute the C-odd beam crossing sign relative to pp
  double OddProfileBeamSign() const;

  // Eigenbasis unitary (orthogonal) rotation matrix
  MMatrix<double> U;

  // Mandelstam s
  double s = 0.0;

  // Initial state
  std::vector<gra::MParticle> INITIALSTATE;

  // Integrated eikonal based total cross sections
  double sigma_tot = 0.0;
  double sigma_el = 0.0;
  double sigma_inel = 0.0;

  // Cut Pomeron probabilities with Simpson weighted joint impact parameter nodes
  std::vector<double> P_cut;
  std::vector<std::vector<double>> P_array;
  std::vector<double> cut_b;
  double MAX_P_array = 0.0;
  double MAX_P_cut = 0.0;
  double min_opacity_eigenvalue = 0.0;

  // Multi-Channel eikonal based cross sections
  MMatrix<double> sigma_diff = {{0.0, 0.0}, {0.0, 0.0}};
  MMatrix<double> sigma_diff_exclusive;

  // Event-by-event screening loop quadrature constants
  LoopConst loop_const;

  // Compute number of radial loop nodes for the current radial rule
  unsigned int GetLoopKTNodeCount() const;
  // Check that the cached quadrature matches the current numerical controls
  bool LoopNumericsMatch() const;
  // Cache the physical proton screening amplitude on the loop grid
  void InitLoopPhysicalScreeningCache();
  // Clear every runtime value owned by one eikonal initialization
  void ResetRuntime();

  friend class MEikonalMatrix;

  static const unsigned int X = 0;
  static const unsigned int Y = 1;
};

} // namespace gra

#endif
