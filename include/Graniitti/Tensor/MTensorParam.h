// Tensor Pomeron parameters
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MTENSORPARAM_H
#define MTENSORPARAM_H

#include <array>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "Graniitti/MModelCache.h"
#include "Graniitti/MModelTune.h"
#include "Graniitti/Particle/MPDG.h"
#include "Graniitti/Particle/MResonance.h"
#include "Graniitti/Photon/MPhotoDiss.h"
#include "Graniitti/Regge/MReggeForm.h"
#include "Graniitti/Regge/MSoftModel.h"
#include "Graniitti/Tensor/MTensorExchange.h"
#include "Graniitti/Tensor/MTensorVector.h"

// Libraries
#include "json.hpp"

namespace gra {

// Select the literature phase prescription for an exchanged vector trajectory
enum class TensorVectorReggePhase { None, Exponential };

// Select the transverse vector-meson line-shape prescription
enum class TensorVectorWidthModel { PWaveTwoBody, Constant, RhoOmega };

// Store the complete immutable Tensor Pomeron model parameters
struct MTensorPomeronParam {
  // Select the high-energy helicity-conserving forward-proton approximation
  bool FORWARD_NOFLIP = false;

  // Select the generic helicity decay phase for TP resonance decays
  bool use_zeta = false;

  // Transfer step used to derive the elementary optical profile
  double photo_dt = 0.0;
  std::map<int, flux::PhotoDissParam> photo_diss;

  MTensorExchangeModel exchange;
  regge::FFParam meson_ff_transfer;

  // Store gauge-restored low-Q2 charged-meson photoproduction steering
  struct PhotoParam {
    bool proton_pauli = false;
    bool photon_exchange = true;
    double m0_2 = 0.0;
    std::map<int, double> lambda;
    std::map<int, std::pair<regge::FFParam, regge::FFParam>> offshell;
    double q2_max = 0.0;
    std::vector<int> exchanges;
  };
  PhotoParam photo;

  // Store one Pomeron-pseudoscalar continuum channel
  struct PseudoscalarParam {
    int pdg = 0;
    double gPPS = 0;
    regge::FFParam ff_offshell;
  };
  std::vector<PseudoscalarParam> pseudoscalars;

  // Store one Pomeron-baryon continuum channel
  struct BaryonParam {
    int pdg = 0;
    double gPBB = 0;
    regge::FFParam ff_offshell;
  };
  std::vector<BaryonParam> baryons;

  // Store one exchanged-vector trajectory and line shape
  struct VectorParam {
    int pdg = 0;
    std::array<double, 2> gPvv = {0, 0};
    regge::FFParam ff_transfer;
    double trajectory_intercept = 0;
    double trajectory_slope = 0;
    regge::FFParam ff_offshell;
    TensorVectorReggePhase regge_phase = TensorVectorReggePhase::None;
    TensorVectorWidthModel width_model = TensorVectorWidthModel::Constant;
    int decay_daughter_pdg = 0;
    double decay_daughter_mass = 0;
    double mass = 0;
    double width = 0;
  };
  std::vector<VectorParam> vectors;
  tensor::RhoOmega rho_omega;

  // Store one vector-meson-dominance transition
  struct VMDParam {
    int pdg = 0;
    double gammaV2 = 0;
    int gammaV_sign = 0;
    double mass = 0;
  };
  std::vector<VMDParam> VMD;

  bool initialized = false;

  // Parse one exchanged-vector Regge phase prescription
  static TensorVectorReggePhase ParseVectorReggePhase(const std::string &mode);

  // Parse one transverse vector-meson width prescription
  static TensorVectorWidthModel ParseVectorWidthModel(const std::string &mode);

  // Find pseudoscalar parameters by absolute PDG id
  const PseudoscalarParam &FindPseudoscalar(int pdg) const;

  // Find baryon parameters by absolute PDG id
  const BaryonParam &FindBaryon(int pdg) const;

  // Find vector parameters by absolute PDG id
  const VectorParam &FindVector(int pdg) const;

  // Find VMD parameters by PDG id
  const VMDParam &FindVMD(int pdg) const;

  // Read parameters from immutable model cards and one PDG snapshot
  void Configure(const nlohmann::json &general, const nlohmann::json &continuum, const std::string &source,
                 const MPDG &pdg_table, double coupling_min);
};

using MTensorPomeronParamPtr = std::shared_ptr<const MTensorPomeronParam>;

// Read one immutable Tensor Pomeron parameter block from a tune
MTensorPomeronParamPtr ReadTensorPomeronParam(const MModelTune &tune,
                                             const MPDG &pdg_table);

// Compute the run owned immutable Tensor Pomeron parameter block
MTensorPomeronParamPtr GetTensorParam(MModelCache &cache,
                                     const MPDG &pdg_table, const std::map<std::string, PARAM_RES> &resonances);

} // namespace gra

#endif
