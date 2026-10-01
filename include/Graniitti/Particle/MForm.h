// Form factors and proton structure functions
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MFORM_H
#define MFORM_H

// C++
#include <string>

namespace gra {

namespace form {

// Store one immutable flat-amplitude slope
struct FlatParam {
  double B = 0.0;

  // Read and validate one PARAM_FLAT block
  static FlatParam Read(const std::string &source_file, const std::string &json_text);
};

// Store one immutable proton-structure parameter selection
struct ParamStore {
  std::string F2 = "CKMT";
  std::string EM = "DIPOLE";
  std::string QED_alpha = "ZERO";

  // Read and validate one PARAM_STRUCTURE block
  static ParamStore Read(const std::string &source_file,
                         const std::string &json_text);
};

// Evaluate an exponential cross-section slope as an amplitude factor
double ExpSlopeAmplitude(double B, double delta_t);

// Evaluate an exponential cross-section slope as a weight factor
double ExpSlopeWeight(double B, double delta_t);

// Proton structure functions
double F2xQ2(double xi, double Q2, const ParamStore &param = {});
double RxQ2(double xi, double Q2);
double FLxQ2(double xi, double Q2, const ParamStore &param = {});
double F1xQ2(double xi, double Q2, const ParamStore &param = {});

// Proton electromagnetic form factors
double F1(double Q2, const ParamStore &param = {});
double F2(double Q2, const ParamStore &param = {});
double G_M(double Q2, const ParamStore &param = {});
double G_E(double Q2, const ParamStore &param = {});

double G_M_DIPOLE(double Q2);
double G_E_DIPOLE(double Q2);

double G_M_KELLY(double Q2);
double G_E_KELLY(double Q2);

double mu_ratio();

} // namespace form
} // namespace gra

#endif
