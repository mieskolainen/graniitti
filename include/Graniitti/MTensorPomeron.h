// Tensor Pomeron amplitudes
//
// (c) 2017-2020 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MTENSORPOMERON_H
#define MTENSORPOMERON_H

// C++
#include <complex>
#include <random>
#include <vector>

// Tensor algebra
#include "FTensor.hpp"

// Own
#include "Graniitti/MGlobals.h"
#include "Graniitti/M4Vec.h"
#include "Graniitti/MDirac.h"
#include "Graniitti/MKinematics.h"
#include "Graniitti/MMath.h"
#include "Graniitti/MTensor.h"


namespace gra {

// Tensor Pomeron model parameters
struct MTensorPomeronParam {

  // Trajectory parameters
  double delta_P  = 0;
  double delta_O  = 0;
  double delta_1R = 0;
  double delta_2R = 0;

  double ap_P  = 0;      // GeV^{-2}
  double ap_O  = 0;      // GeV^{-2}
  double ap_1R = 0;      // GeV^{-2}
  double ap_2R = 0;      // GeV^{-2}

  double eta_O = 1;      // +- 1
  
  // Scales
  double M_O  = 0;       // GeV
  double M_1R = 0;       // GeV
  
  // Couplings
  double gPNN   = 0;  // GeV^{-1}, Pomeron-Proton-Proton
  double gPpipi = 0;  // GeV^{-1}, Pomeron-Pion-Pion
  double gPKK   = 0;  // GeV^{-1}, Pomeron-Kaon-Kaon

  std::vector<double> gPrhorho; // [GeV^{-3}, GeV^{-1}], Pomeron-Rho-Rho
  std::vector<double> gPphiphi; // [GeV^{-3}, GeV^{-1}], Pomeron-Phi-Phi

  // QED
  double mu_ratio = 0;   // Proton magnetic moment in magneton units (mu_p / mu_N)
  double alpha_QED = 0;  // alpha_QED at q^2 ~ 0
  
  bool   initialized = false;

  // Read parameters from file
  void ReadParameters(const std::string& modelfile) {

    using json = nlohmann::json;
    const std::string data = gra::aux::GetInputData(modelfile);
    json j;

    try {
      j = json::parse(data);

      // JSON block identifier
      const std::string XID = "PARAM_TENSOR_POMERON";

      delta_P   = j.at(XID).at("delta_P");
      delta_O   = j.at(XID).at("delta_O");
      delta_1R  = j.at(XID).at("delta_1R");
      delta_2R  = j.at(XID).at("delta_2R");

      ap_P      = j.at(XID).at("ap_P");
      ap_O      = j.at(XID).at("ap_O");
      ap_1R     = j.at(XID).at("ap_1R");
      ap_2R     = j.at(XID).at("ap_2R");

      eta_O     = j.at(XID).at("eta_O");
      M_O       = j.at(XID).at("M_O");
      M_1R      = j.at(XID).at("M_1R");

      gPNN   = j.at(XID).at("gPNN");
      gPpipi = j.at(XID).at("gPpipi");
      gPKK   = j.at(XID).at("gPKK");
      
      std::vector<double> a = j.at(XID).at("gPrhorho");
      gPrhorho  = a;

      std::vector<double> b = j.at(XID).at("gPphiphi");
      gPphiphi  = b;
      
      alpha_QED = j.at(XID).at("alpha_QED");
      mu_ratio  = j.at(XID).at("mu_ratio");

      initialized = true;
    } catch (...) {
      std::string str = "MTensorPomeronParam::ReadParameters: Error parsing " + modelfile +
                        " (Check for extra/missing commas)";
      throw std::invalid_argument(str);
    }
  }
};


// Matrix element dimension: " GeV^" << -(2*external_legs - 8)
class MTensorPomeron : public MDirac {
 public:
  MTensorPomeron(gra::LORENTZSCALAR& lts, const std::string& modelfile);
  ~MTensorPomeron() {}

  // Decay coupling (static so we can call it independently)
  static double GDecay(int J, double M, double Gamma, double mf, double BR);

  // Amplitude squared
  double ME3(gra::LORENTZSCALAR &lts) const;
  double ME4(gra::LORENTZSCALAR &lts) const;
  double ME6(gra::LORENTZSCALAR &lts) const;

  // Scalar, Pseudoscalar, Tensor coupling structures
  FTensor::Tensor4<std::complex<double>, 4, 4, 4, 4> iG_PPS_0() const;
  FTensor::Tensor4<std::complex<double>, 4, 4, 4, 4> iG_PPS_1(const M4Vec &q1, const M4Vec &q2,
                                                              double g_PPS) const;
  FTensor::Tensor4<std::complex<double>, 4, 4, 4, 4> iG_PPPS_0(const M4Vec &q1, const M4Vec &q2,
                                                               double g_PPPS) const;
  FTensor::Tensor4<std::complex<double>, 4, 4, 4, 4> iG_PPPS_1(const M4Vec &q1, const M4Vec &q2,
                                                               double g_PPPS) const;
  MTensor<std::complex<double>>                      iG_PPT_00() const;
  MTensor<std::complex<double>> iG_PPT_12(const M4Vec &q1, const M4Vec &q2, double g_PPT,
                                          int mode) const;
  MTensor<std::complex<double>> iG_PPT_03(const M4Vec &q1, const M4Vec &q2, double g_PPT) const;
  MTensor<std::complex<double>> iG_PPT_04(const M4Vec &q1, const M4Vec &q2, double g_PPT) const;
  MTensor<std::complex<double>> iG_PPT_05(const M4Vec &q1, const M4Vec &q2, double g_PPT) const;
  MTensor<std::complex<double>> iG_PPT_06(const M4Vec &q1, const M4Vec &q2, double g_PPT) const;

  // Vertex functions
  FTensor::Tensor2<std::complex<double>, 4, 4> iG_vv2psps(const std::vector<M4Vec> &p,
                                                          int                       PDG) const;

  FTensor::Tensor4<std::complex<double>, 4, 4, 4, 4> iG_PPS_total(
      const M4Vec &q1, const M4Vec &q2, double M0, const std::string &mode,
      const std::vector<double> &g_PPS) const;
  MTensor<std::complex<double>> iG_PPT_total(const M4Vec &q1, const M4Vec &q2, double M0,
                                             const std::vector<double> &g_PPT) const;

  FTensor::Tensor2<std::complex<double>, 4, 4> iG_PppHE(const M4Vec &prime, const M4Vec p) const;

  FTensor::Tensor1<std::complex<double>, 4> iG_yee(
      const M4Vec &prime, const M4Vec &p, const std::vector<std::complex<double>> &ubar,
      const std::vector<std::complex<double>> &u) const;

  FTensor::Tensor2<std::complex<double>, 4, 4> iG_yeebary(
      const std::vector<std::complex<double>> &ubar, const MMatrix<std::complex<double>> &iSF,
      const std::vector<std::complex<double>> &v) const;

  FTensor::Tensor1<std::complex<double>, 4> iG_ypp(
      const M4Vec &prime, const M4Vec &p, const std::vector<std::complex<double>> &ubar,
      const std::vector<std::complex<double>> &u) const;

  FTensor::Tensor2<std::complex<double>, 4, 4> iG_Ppp(
      const M4Vec &prime, const M4Vec &p, const std::vector<std::complex<double>> &ubar,
      const std::vector<std::complex<double>> &u) const;

  FTensor::Tensor4<std::complex<double>, 4, 4, 4, 4> iG_PppbarP(
      const M4Vec &prime, const std::vector<std::complex<double>> &ubar, const M4Vec &pt,
      const MMatrix<std::complex<double>> &iSF, const std::vector<std::complex<double>> &v,
      const M4Vec &p) const;

  FTensor::Tensor2<std::complex<double>, 4, 4>       iG_Ppsps(const M4Vec &prime, const M4Vec &p,
                                                              double g1) const;
  FTensor::Tensor4<std::complex<double>, 4, 4, 4, 4> iG_Pvv(const M4Vec &prime, const M4Vec &p,
                                                            double g1, double g2) const;

  std::complex<double> iG_f0ss(const M4Vec &p3, const M4Vec &p4, double M0, double g1) const;
  FTensor::Tensor2<std::complex<double>, 4, 4> iG_f0vv(const M4Vec &p3, const M4Vec &p4, double M0,
                                                       double g1, double g2) const;

  FTensor::Tensor2<std::complex<double>, 4, 4> iG_psvv(const M4Vec &p3, const M4Vec &p4, double M0,
                                                       double g1) const;
  FTensor::Tensor1<std::complex<double>, 4>    iG_vpsps(const M4Vec &k1, const M4Vec &k2, double M0,
                                                        double g1) const;

  FTensor::Tensor2<std::complex<double>, 4, 4>       iG_f2psps(const M4Vec &k1, const M4Vec &k2,
                                                               double M0, double g1) const;
  FTensor::Tensor4<std::complex<double>, 4, 4, 4, 4> iG_f2vv(const M4Vec &k1, const M4Vec &k2,
                                                             double M0, double g1, double g2) const;
  FTensor::Tensor4<std::complex<double>, 4, 4, 4, 4> iG_f2yy(const M4Vec &k1, const M4Vec &k2,
                                                             double M0, double g1, double g2) const;

  FTensor::Tensor2<std::complex<double>, 4, 4> iG_yV(double q2, int pdg) const;

  // Polarization sums
  std::vector<FTensor::Tensor2<std::complex<double>, 4, 4>> MassiveSpin1PolSum(
      const FTensor::Tensor4<std::complex<double>, 4, 4, 4, 4> &M, const M4Vec &p3,
      const M4Vec &p4) const;
  std::vector<FTensor::Tensor2<std::complex<double>, 4, 4>> MasslessSpin1PolSum(
      const FTensor::Tensor4<std::complex<double>, 4, 4, 4, 4> &M, const M4Vec &p3,
      const M4Vec &p4) const;

  std::vector<std::complex<double>> MasslessSpin1PolSum(
      const FTensor::Tensor2<std::complex<double>, 4, 4> &M, const M4Vec &p3,
      const M4Vec &p4) const;
  std::vector<std::complex<double>> MassiveSpin1PolSum(
      const FTensor::Tensor2<std::complex<double>, 4, 4> &M, const M4Vec &p3,
      const M4Vec &p4) const;

  // Propagators
  FTensor::Tensor4<std::complex<double>, 4, 4, 4, 4> iD_P(double s, double t) const;
  FTensor::Tensor4<std::complex<double>, 4, 4, 4, 4> iD_2R(double s, double t) const;
  FTensor::Tensor2<std::complex<double>, 4, 4>       iD_O(double s, double t) const;
  FTensor::Tensor2<std::complex<double>, 4, 4>       iD_1R(double s, double t) const;

  FTensor::Tensor2<std::complex<double>, 4, 4> iD_V(const M4Vec &p, double M0, double s34) const;
  std::complex<double>                         iD_MES0(const M4Vec &p, double M0) const;
  std::complex<double> iD_MES(const M4Vec &p, double M0, double Gamma) const;
  FTensor::Tensor2<std::complex<double>, 4, 4> iD_VMES(const M4Vec &p, double M0, double Gamma,
                                                       bool INDEX_UP) const;
  FTensor::Tensor4<std::complex<double>, 4, 4, 4, 4> iD_TMES(const M4Vec &p, double M0,
                                                             double Gamma, bool INDEX_UP) const;

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
  double F1(double t) const;
  double F2(double t) const;

  double F1_(double t) const;
  double F2_(double t) const;

  double GD(double t) const;
  double FM(double q2, double LAMBDA) const;
  double F(double q2, double M0, double LAMBDA) const;

 private:
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

  // Parameters
  MTensorPomeronParam Param;
};

}  // namespace gra

#endif
