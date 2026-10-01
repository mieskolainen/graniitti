// Tensor Pomeron amplitudes
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MTENSORPOMERON_H
#define MTENSORPOMERON_H

// C++
#include <complex>
#include <memory>
#include <vector>

// Tensor algebra
#include "FTensor.hpp"

// Own
#include "Graniitti/Kinematics/M4Vec.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Math/MTensor.h"
#include "Graniitti/Process/MAmpMatch.h"
#include "Graniitti/Process/MProcessSetup.h"
#include "Graniitti/Process/MProcessState.h"
#include "Graniitti/Regge/MReggeForm.h"
#include "Graniitti/Tensor/MTensorForward.h"
#include "Graniitti/Photon/MFlux.h"
#include "Graniitti/Tensor/MTensorParam.h"
#include "Graniitti/Tensor/MTensorProcess.h"
#include "Graniitti/Spin/MDirac.h"

namespace gra {

// Select the physical interaction used by the shared covariant 2 to 4
// implementation
enum class TensorContinuumMode { TensorPomeron, QED };

// Select the elastic photon-fermion vertex used by one Tensor amplitude
enum class TensorPhotonVertex { DiracPauli, Dirac };

// Matrix element dimension: " GeV^" << -(2*external_legs - 8)
class MTensorPomeron : public MDirac, public amplitude::ProcessFamily {
 public:
  MTensorPomeron(gra::LORENTZSCALAR &lts, MModelTunePtr tune,
                 std::shared_ptr<const amplitude::ProcessDefinition> definition,
                 MTensorPomeronMode                                  mode = MTensorPomeronMode::Generic);
  ~MTensorPomeron() {}

  // Compute the forward proton spin layout owned by one Tensor amplitude mode
  bool ForwardNoFlip(MTensorPomeronMode mode, const gra::LORENTZSCALAR &lts) const;

  // Build one immutable Tensor Pomeron process definition for a selected mode
  static std::shared_ptr<const amplitude::ProcessDefinition> ProcessDefinitionFor(MTensorPomeronMode mode);

  // Initialize immutable Tensor Pomeron parameters before worker process copies
  static void InitializeParameters(const MProcessSetup &setup);

  // Initialize tensor Pomeron resonance and continuum branching structures
  static void InitializeBranching(MProcessSetup &setup, MTensorPomeronMode mode);

  // Compute the direct stable-final-state tensor-amplitude decay structure
  static constexpr DecayStructure DirectDecayStructure() { return {DecayType::Full}; }

  // Decay coupling (static so we can call it independently)
  static double GDecay(int J, double M, double Gamma, double mf, double BR, double symmetry = 1.0);

  // Compute the dimensionless pseudoscalar to two-photon coupling
  static double GDecayPseudoscalarGammaGamma(double M, double Gamma, double BR);

  // Amplitude squared
  double ME3(gra::LORENTZSCALAR &lts, MTensorPomeronMode mode = MTensorPomeronMode::Resonance) const;
  double ME4(gra::LORENTZSCALAR &lts, TensorContinuumMode mode) const;
  double ME6(gra::LORENTZSCALAR &lts) const;

  // Evaluate gauge-restored charged-meson photo and electroproduction
  double MEPhoto(gra::LORENTZSCALAR &lts) const;

  // Evaluate Tensor photoproduction with physical lepton and nuclear beams
  double PhotoProduction(gra::LORENTZSCALAR &lts, bool continuum) const;

  // Compute the HERA target current with the elastic forward normalization at W0
  FTensor::Tensor2<std::complex<double>, 4, 4> PhotoDissCurrent(
      const ForwardLegState &state, int exchange_pdg, int vector_pdg, double nu, double w2) const;

  // Derive the transverse vector-nucleon optical profile from Tensor production couplings
  flux::PhotoTargetProfile PhotoProfile(const PARAM_RES &res, double w2) const;

  // Scalar, Pseudoscalar, Tensor coupling structures
  FTensor::Tensor4<std::complex<double>, 4, 4, 4, 4> iG_PPS_0() const;
  FTensor::Tensor4<std::complex<double>, 4, 4, 4, 4> iG_PPS_1(const M4Vec &q1, const M4Vec &q2, double g_PPS) const;
  FTensor::Tensor4<std::complex<double>, 4, 4, 4, 4> iG_PPPS_0(const M4Vec &q1, const M4Vec &q2, double g_PPPS) const;
  FTensor::Tensor4<std::complex<double>, 4, 4, 4, 4> iG_PPPS_1(const M4Vec &q1, const M4Vec &q2, double g_PPPS) const;
  MTensor<std::complex<double>> iG_PPT_12(const M4Vec &q1, const M4Vec &q2, double g_PPT, int mode) const;
  MTensor<std::complex<double>> iG_PPT_03(const M4Vec &q1, const M4Vec &q2, double g_PPT) const;
  MTensor<std::complex<double>> iG_PPT_04(const M4Vec &q1, const M4Vec &q2, double g_PPT) const;
  MTensor<std::complex<double>> iG_PPT_05(const M4Vec &q1, const M4Vec &q2, double g_PPT) const;
  MTensor<std::complex<double>> iG_PPT_06(const M4Vec &q1, const M4Vec &q2, double g_PPT) const;
  MTensor<std::complex<double>> iG_PPA_22(const M4Vec &q1, const M4Vec &q2, double g_PPA) const;
  MTensor<std::complex<double>> iG_PPA_44(const M4Vec &q1, const M4Vec &q2, double g_PPA) const;

  // Vertex functions
  FTensor::Tensor2<std::complex<double>, 4, 4> iG_vv2psps(const std::vector<M4Vec> &p, int PDG,
                                                          double decay_coupling, const regge::FFParam &ff_decay) const;

  FTensor::Tensor4<std::complex<double>, 4, 4, 4, 4> iG_PPS_total(const M4Vec &q1, const M4Vec &q2, double M0,
                                                                  TensorResonanceType        type,
                                                                  const std::vector<double> &g_PPS,
                                                                  const RES_TENSOR_CHANNEL &channel) const;
  MTensor<std::complex<double>>                      iG_PPT_total(const M4Vec &q1, const M4Vec &q2, double M0,
                                                                  const std::vector<double> &g_PPT,
                                                                  const RES_TENSOR_CHANNEL &channel) const;
  MTensor<std::complex<double>>                      iG_PPA_total(const M4Vec &q1, const M4Vec &q2, double M0,
                                                                  const RES_TENSOR_CHANNEL &channel) const;

  FTensor::Tensor2<std::complex<double>, 4, 4> iG_PppHE(const M4Vec &prime, const M4Vec p) const;

  // Compute the elastic or inclusive high-energy Tensor Pomeron beam current
  FTensor::Tensor2<std::complex<double>, 4, 4> iG_PForwardHE(const ForwardLegState &state) const;

  // Compute one elastic exact or inclusive identity Tensor Pomeron beam current
  FTensor::Tensor2<std::complex<double>, 4, 4> iG_PForward(const ForwardLegState &state, const MDirac::Spinor &ubar,
                                                           const MDirac::Spinor &u, std::size_t initial_helicity,
                                                           std::size_t final_helicity) const;

  // Compute the high-energy current for one configured rank-two exchange
  FTensor::Tensor2<std::complex<double>, 4, 4> iG_TForwardHE(const ForwardLegState &state, int exchange_pdg) const;

  // Compute one exact elastic current for a configured rank-two exchange
  FTensor::Tensor2<std::complex<double>, 4, 4> iG_TForward(const ForwardLegState &state, int exchange_pdg,
                                                           const MDirac::Spinor &ubar, const MDirac::Spinor &u,
                                                           std::size_t initial_helicity,
                                                           std::size_t final_helicity) const;

  // Compute the high-energy current for one configured vector exchange
  FTensor::Tensor1<std::complex<double>, 4> iG_VForwardHE(const ForwardLegState &state, int exchange_pdg) const;

  // Compute one exact elastic current for a configured vector exchange
  FTensor::Tensor1<std::complex<double>, 4> iG_VForward(const ForwardLegState &state, int exchange_pdg,
                                                        const MDirac::Spinor &ubar, const MDirac::Spinor &u,
                                                        std::size_t initial_helicity, std::size_t final_helicity) const;

  // Compute orthogonal elastic or inclusive EPA sources after photon transfer
  std::vector<FTensor::Tensor1<std::complex<double>, 4>> iG_yForwardSources(
      const LORENTZSCALAR &lts, const ForwardLegState &state, const MDirac::Spinor &ubar, const MDirac::Spinor &u,
      std::size_t initial_helicity, std::size_t final_helicity,
      TensorPhotonVertex vertex = TensorPhotonVertex::DiracPauli) const;

  FTensor::Tensor1<std::complex<double>, 4> iG_yee(const M4Vec &prime, const M4Vec &p, const MDirac::Spinor &ubar,
                                                   const MDirac::Spinor &u) const;

  FTensor::Tensor2<std::complex<double>, 4, 4> iG_yeebary(const MDirac::Spinor                &ubar,
                                                          const MMatrix<std::complex<double>> &iSF,
                                                          const MDirac::Spinor                &v) const;

  FTensor::Tensor1<std::complex<double>, 4> iG_ypp(const M4Vec &prime, const M4Vec &p, const MDirac::Spinor &ubar,
                                                   const MDirac::Spinor &u) const;

  FTensor::Tensor2<std::complex<double>, 4, 4> iG_Ppp(const M4Vec &prime, const M4Vec &p, const MDirac::Spinor &ubar,
                                                      const MDirac::Spinor &u) const;

  FTensor::Tensor4<std::complex<double>, 4, 4, 4, 4> iG_PppbarP(const M4Vec &prime, const MDirac::Spinor &ubar,
                                                                const M4Vec                         &pt,
                                                                const MMatrix<std::complex<double>> &iSF,
                                                                const MDirac::Spinor &v, const M4Vec &p,
                                                                double gPBB, const regge::FFParam &ff_transfer) const;

  FTensor::Tensor2<std::complex<double>, 4, 4>       iG_Ppsps(const M4Vec &prime, const M4Vec &p, double g1) const;
  FTensor::Tensor4<std::complex<double>, 4, 4, 4, 4> iG_Pvv(const M4Vec &prime, const M4Vec &p, double g1, double g2,
                                                            const regge::FFParam &ff_transfer, bool forward = false) const;

  // Compute one exchange-specific pseudoscalar continuum vertex
  FTensor::Tensor2<std::complex<double>, 4, 4> iG_Tpsps(const M4Vec &prime, const M4Vec &p, int exchange_pdg,
                                                        int hadron_pdg) const;

  // Compute one exchange-specific vector continuum vertex
  FTensor::Tensor4<std::complex<double>, 4, 4, 4, 4> iG_Tvv(const M4Vec &prime, const M4Vec &p, int exchange_pdg,
                                                            int hadron_pdg) const;

  // Compute one vector-Reggeon pseudoscalar continuum vertex
  FTensor::Tensor1<std::complex<double>, 4> iG_Vpsps(const M4Vec &prime, const M4Vec &p, int exchange_pdg,
                                                     int hadron_pdg) const;

  std::complex<double>                         iG_f0ss(const M4Vec &p3, const M4Vec &p4, double M0, double g1,
                                                       const regge::FFParam &ff_decay) const;
  FTensor::Tensor2<std::complex<double>, 4, 4> iG_f0vv(const M4Vec &p3, const M4Vec &p4, double M0, double g1,
                                                       double g2, const regge::FFParam &ff_decay) const;

  FTensor::Tensor2<std::complex<double>, 4, 4> iG_psvv(const M4Vec &p3, const M4Vec &p4, double M0, double g1,
                                                       const regge::FFParam &ff_decay) const;
  FTensor::Tensor1<std::complex<double>, 4>    iG_vpsps(const M4Vec &k1, const M4Vec &k2, double M0, double g1,
                                                        const regge::FFParam &ff_decay) const;

  FTensor::Tensor2<std::complex<double>, 4, 4>       iG_f2psps(const M4Vec &k1, const M4Vec &k2, double M0, double g1,
                                                               const regge::FFParam &ff_decay) const;
  FTensor::Tensor4<std::complex<double>, 4, 4, 4, 4> iG_f2vv(const M4Vec &k1, const M4Vec &k2, double M0, double g1,
                                                             double g2, const regge::FFParam &ff_decay) const;
  FTensor::Tensor4<std::complex<double>, 4, 4, 4, 4> iG_f2yy(const M4Vec &k1, const M4Vec &k2, double M0, double g1,
                                                             double g2, const regge::FFParam &ff_decay) const;

  FTensor::Tensor2<std::complex<double>, 4, 4> iG_yV(double q2, int pdg) const;

  // Polarization sums
  std::vector<FTensor::Tensor2<std::complex<double>, 4, 4>> MassiveSpin1PolSum(
      const FTensor::Tensor4<std::complex<double>, 4, 4, 4, 4> &M, const M4Vec &p3, const M4Vec &p4) const;
  std::vector<FTensor::Tensor2<std::complex<double>, 4, 4>> MasslessSpin1PolSum(
      const FTensor::Tensor4<std::complex<double>, 4, 4, 4, 4> &M, const M4Vec &p3, const M4Vec &p4) const;

  std::vector<std::complex<double>> MasslessSpin1PolSum(const FTensor::Tensor2<std::complex<double>, 4, 4> &M,
                                                        const M4Vec &p3, const M4Vec &p4) const;
  std::vector<std::complex<double>> MassiveSpin1PolSum(const FTensor::Tensor2<std::complex<double>, 4, 4> &M,
                                                       const M4Vec &p3, const M4Vec &p4) const;

  // Propagators
  // Contract one forward tensor current with the spin-2 Pomeron propagator
  FTensor::Tensor2<std::complex<double>, 4, 4> PomeronPropagatorCurrent(
      const FTensor::Tensor2<std::complex<double>, 4, 4> &current, double s, double t) const;

  // Project one tensor current with a precomputed Pomeron Regge factor
  FTensor::Tensor2<std::complex<double>, 4, 4> PomeronPropagatorCurrent(
      const FTensor::Tensor2<std::complex<double>, 4, 4> &current, const std::complex<double> &factor) const;

  // Compute the scalar Regge factor of the spin-2 Pomeron propagator
  std::complex<double> PomeronPropagatorFactor(double s, double t) const;

  // Compute the scalar Regge factor of any configured rank-two exchange
  std::complex<double> TensorPropagatorFactor(int exchange_pdg, double s, double t) const;

  // Contract one vector current with a configured vector propagator
  FTensor::Tensor1<std::complex<double>, 4> VectorPropagatorCurrent(
      const FTensor::Tensor1<std::complex<double>, 4> &current, int exchange_pdg, double s, double t) const;

  // Contract one Pomeron current into a baryon vertex Dirac matrix
  MMatrix<std::complex<double>> PomeronBaryonCurrent(const FTensor::Tensor2<std::complex<double>, 4, 4> &current,
                                                     const M4Vec &prime, const M4Vec &p, double gPBB,
                                                     const regge::FFParam &ff_transfer) const;

  // Contract one rank-two exchange current into a continuum baryon vertex
  MMatrix<std::complex<double>> TensorBaryonCurrent(const FTensor::Tensor2<std::complex<double>, 4, 4> &current,
                                                    const M4Vec &prime, const M4Vec &p, int exchange_pdg,
                                                    int baryon_pdg) const;

  // Contract one vector exchange current into a continuum baryon vertex
  MMatrix<std::complex<double>> VectorBaryonCurrent(const FTensor::Tensor1<std::complex<double>, 4> &current,
                                                    const M4Vec &prime, const M4Vec &p, int exchange_pdg,
                                                    int baryon_pdg) const;

  FTensor::Tensor4<std::complex<double>, 4, 4, 4, 4> iD_P(double s, double t) const;
  FTensor::Tensor4<std::complex<double>, 4, 4, 4, 4> iD_2R(double s, double t) const;
  FTensor::Tensor2<std::complex<double>, 4, 4>       iD_VExchange(int exchange_pdg, double s, double t) const;
  FTensor::Tensor2<std::complex<double>, 4, 4>       iD_O(double s, double t) const;
  FTensor::Tensor2<std::complex<double>, 4, 4>       iD_1R(double s, double t) const;

  FTensor::Tensor2<std::complex<double>, 4, 4> iD_V(const M4Vec &p, double M0, double s34, int pdg) const;
  FTensor::Tensor2<std::complex<double>, 4, 4> iD_VMD(const M4Vec &p, int pdg) const;
  std::complex<double>                         iD_MES0(const M4Vec &p, double M0) const;
  std::complex<double>                         iD_MES(const M4Vec &p, double M0, double Gamma) const;
  FTensor::Tensor2<std::complex<double>, 4, 4> iD_VMES(const M4Vec &p, double M0, double Gamma, int pdg, bool INDEX_UP,
                                                       bool CONSERVED_CURRENT) const;
  FTensor::Tensor4<std::complex<double>, 4, 4, 4, 4> iD_TMES(const M4Vec &p, double M0, double Gamma,
                                                             bool INDEX_UP) const;

  // Compute the transverse vector spectral function with the selected line shape
  std::complex<double> VectorSpectrum(double s, double mass, double width, int pdg) const;

  // Compute a propagated vector decay current, including final rho-omega mixing
  FTensor::Tensor1<std::complex<double>, 4> VectorDecay(const M4Vec &positive, const M4Vec &negative,
      double mass, double width, int pdg, double coupling, const regge::FFParam &form) const;

  // Tensor functions
  FTensor::Tensor4<std::complex<double>, 4, 4, 4, 4> Gamma0(const M4Vec &k1, const M4Vec &k2) const;
  FTensor::Tensor4<std::complex<double>, 4, 4, 4, 4> Gamma2(const M4Vec &k1, const M4Vec &k2) const;
  void                                               CalcRTensor();

  // Trajectories
  double alpha_P(double t) const;
  double alpha_O(double t) const;
  double alpha_1R(double t) const;
  double alpha_2R(double t) const;

  // Form factors
  double F1_(double t) const;
  double F2_(double t) const;

  double GD(double t) const;

 private:
  // Store photon-tensor and vector-tensor kernels in both beam orders
  struct VectorKernels {
    std::map<int, FTensor::Tensor3<std::complex<double>, 4, 4, 4>> yT, Ty;
    std::map<std::pair<int, int>, FTensor::Tensor3<std::complex<double>, 4, 4, 4>> VT, TV;
  };

  // Build the shared covariant vector production kernels including the physical decay
  VectorKernels VectorProduction(const LORENTZSCALAR &lts, const PARAM_RES &res,
                                 const FTensor::Tensor1<std::complex<double>, 4> &decay_current) const;


  // Build the effective gauge-invariant spin-three photon-tensor kernel
  VectorKernels Spin3Production(const LORENTZSCALAR &lts, const PARAM_RES &res, const TensorDecayState &decay) const;

  // Compute whether one production coupling is large enough to evaluate
  bool ActiveProductionCoupling(double coupling) const;

  // Compute the common factorized PP-resonance form factor
  double PPResonanceFormFactor(const M4Vec &q1, const M4Vec &q2, double M0,
                               const RES_TENSOR_CHANNEL &channel) const;

  // Contract the PP-tensor-resonance vertex directly with both Pomeron currents
  FTensor::Tensor2<std::complex<double>, 4, 4> iG_PPT_contract(
      const FTensor::Tensor2<std::complex<double>, 4, 4> &left,
      const FTensor::Tensor2<std::complex<double>, 4, 4> &right, const M4Vec &q1, const M4Vec &q2, double M0,
      const std::vector<double> &g_PPT, const RES_TENSOR_CHANNEL &channel) const;

  // Compute one unsymmetrized Gamma8 element for the axial (2,2) vertex
  double GammaPPA22Base(std::size_t kappa, std::size_t lambda, std::size_t rho, std::size_t sigma, std::size_t mu,
                        std::size_t nu, std::size_t alpha, std::size_t beta) const;

  // Minkowski metric tensor
  FTensor::Tensor2<double, 4, 4> gT;

  // Epsilon tensors
  FTensor::Tensor4<double, 4, 4, 4, 4> eps_lo;
  FTensor::Tensor4<double, 4, 4, 4, 4> eps_hi;

  // Aux tensors
  FTensor::Tensor4<double, 4, 4, 4, 4> R_DDDD;
  FTensor::Tensor4<double, 4, 4, 4, 4> R_DDDU;
  FTensor::Tensor4<double, 4, 4, 4, 4> R_DDUU;
  FTensor::Tensor4<double, 4, 4, 4, 4> R_UUDD;

  MTensor<double> T1;
  MTensor<double> T2;
  MTensor<double> T3;

  // Pre-calculated kinematics-independent resonance tensors
  MTensor<std::complex<double>> PPT_00;
  MTensor<double>               PPA_22;

  MModelTunePtr                              model_tune;
  SoftModelPtr                               soft_model;
  std::shared_ptr<const MTensorPomeronParam> tensor_param_handle;
  const MTensorPomeronParam                 &tensor_param;
};

}  // namespace gra

#endif
