// "Durham QCD" Processes and Amplitudes
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MDURHAM_H
#define MDURHAM_H

// C++
#include <complex>
#include <random>
#include <vector>

// MADGRAPH
#include "Graniitti/Amplitude/MAmpMG5_gg_gg.h"
#include "Graniitti/Amplitude/MAmpMG5_gg_ggg.h"
#include "Graniitti/Amplitude/MAmpMG5_gg_qqbar.h"

// Own
#include "Graniitti/M4Vec.h"
#include "Graniitti/MAux.h"
#include "Graniitti/MGlobals.h"
#include "Graniitti/MKinematics.h"
#include "Graniitti/MSudakov.h"

namespace gra {

// Durham loop integral discretization technicals
struct MDurhamParam {
  unsigned int N_qt  = 0;  // (> 30)
  unsigned int N_phi = 0;  // (> 10)

  double qt2_MIN = 0;  // Loop momentum qt^2 minimum (GeV^2)
  double qt2_MAX = 0;  // Loop momentum qt^2 maximum (GeV^2)

  std::string PDF_scale    = "MIN";  // Scheme
  double      alphas_scale = 4.0;    // PDF factorization scale
  double      MAXCOS       = 0.9;    // Meson amplitude |cos(theta*)| < MAXCOS

  // THESE ARE CALCULATED FROM ABOVE
  double qt_MIN      = 0;
  double qt_MAX      = 0;
  double qt_STEP     = 0;
  double phi_STEP    = 0;
  bool   initialized = false;

  // Read parameters from file
  void ReadParameters() {
    using json = nlohmann::json;

    const std::string inputfile =
        gra::aux::GetBasePath(2) + "/modeldata/" + gra::MODELPARAM + "/GENERAL.json";
    const std::string data = gra::aux::GetInputData(inputfile);
    json              j;

    try {
      j = json::parse(data);

      // JSON block identifier
      const std::string XID = "PARAM_DURHAM_QCD";
      N_qt                  = j.at(XID).at("N_qt");
      N_phi                 = j.at(XID).at("N_phi");
      qt2_MIN               = j.at(XID).at("qt2_MIN");
      qt2_MAX               = j.at(XID).at("qt2_MAX");
      PDF_scale             = j.at(XID).at("PDF_scale");
      alphas_scale          = j.at(XID).at("alphas_scale");
      MAXCOS                = j.at(XID).at("MAXCOS");

      // Now calculate rest
      qt_MIN = math::msqrt(qt2_MIN);
      qt_MAX = math::msqrt(qt2_MAX);

      qt_STEP  = (qt_MIN - qt_MAX) / N_qt;
      phi_STEP = (2.0 * math::PI) / N_phi;

      initialized = true;
    } catch (...) {
      std::string str = "Durham::ReadParameters: Error parsing " + inputfile +
                        " (Check for extra/missing commas)";
      throw std::invalid_argument(str);
    }
  }
};

// Matrix element dimension: " GeV^" << -(2*external_legs - 8)
class MDurham {
 public:
  MDurham() {}
  ~MDurham() {}

  double      DurhamQCD(gra::LORENTZSCALAR &lts, const std::string &process);
  double      DQtloop(gra::LORENTZSCALAR &lts, std::vector<std::vector<std::complex<double>>> Amp);
  inline void DScaleChoise(double qt2, double q1_2, double q2_2, double &Q1_2_scale,
                           double &Q2_2_scale) const;

  inline void Dgg2chic0(const gra::LORENTZSCALAR &                      lts,
                        std::vector<std::vector<std::complex<double>>> &Amp,
                        const std::vector<double> &qt1, const std::vector<double> &qt2) const;

  inline void DHelicity(const std::vector<double> &q1, const std::vector<double> &q2,
                        std::vector<std::complex<double>> &JzP) const;

  inline std::complex<double> DHelProj(const std::vector<std::complex<double>> &A,
                                       const std::vector<std::complex<double>> &JzP) const;

  void Dgg2gg(const gra::LORENTZSCALAR &lts, std::vector<std::vector<std::complex<double>>> &Amp);

  void Dgg2qqbar(const gra::LORENTZSCALAR &                      lts,
                 std::vector<std::vector<std::complex<double>>> &Amp);

  void                Dgg2MMbar(const gra::LORENTZSCALAR &                      lts,
                                std::vector<std::vector<std::complex<double>>> &Amp);
  double              phi_CZ(double x, double fM) const;
  std::vector<double> EvalPhi(int N, int pdg) const;

  double Asum = 0.0;
  double Nsum = 0.0;

  // ------------------------------------------------------------------
  // MadGraph amplitudes here
  MAmpMG5_gg_gg    AmpMG5_gg_gg;
  MAmpMG5_gg_ggg   AmpMG5_gg_ggg;
  MAmpMG5_gg_qqbar AmpMG5_gg_qqbar;


 private:
  // Parameters
  MDurhamParam Param;
};

}  // namespace gra

#endif