// Regge form factors and exchanged-meson kernels
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MREGGEFORM_H
#define MREGGEFORM_H

#include <array>
#include <complex>
#include <string>
#include <vector>

#include "Graniitti/Regge/MReggeParam.h"

namespace gra::regge {

// Parse one steering card form factor type
FFType ParseFFType(const std::string &type);

// Parse one steering card form factor normalization point
FFNorm ParseFFNorm(const std::string &norm);

// Read one minimal named form factor object
FFParam ReadFF(const nlohmann::json &block, const std::string &context);

// Validate one positive generated invariant mass squared
void CheckMass2(double mass2, const std::string &context);

// Validate one parsed form factor
void CheckFF(const FFParam &ff, const std::string &context);

// Evaluate one form factor at q^2 with its pole mass squared
double FormFactor(double q2, double pole2, const FFParam &ff);

// Evaluate the optional meson-side Reggeon vertex form factor
double TransferFF(double t, const FFParam &ff);

// Evaluate the invariant-mass resonance form factor
double MassFF(double mass2, double pole2, const FFParam &ff);

// Compute F(q1), F(q2) and their divided difference for none, exp or gaussian profiles
std::array<double, 3> MassFFPair(double q1, double q2, double pole2, const FFParam &ff);

// Evaluate one exchanged-meson propagator
std::complex<double> OffshellProp(double t_hat, double mass2, const ReggeizeParam &reggeize, const MesonTraj &trajectory, const M4Vec &left, const M4Vec &right);

// Evaluate one pair-specific exchanged-meson propagator
std::complex<double> MesonPropagator(const Param &param, double t_hat, double mass2, const PairParam &entry, const MDecayBranch &left, const MDecayBranch &right);

// Evaluate the two vertex form factors of one exchanged-meson line
double MesonVertexFactor(double t_hat, double mass2, const VertexParam &vertex);

// Evaluate one multi-Regge exchanged-meson line
std::complex<double> MesonExchange(const Param &param, double t_hat, double mass2, const PairParam &entry, const VertexParam &vertex, const MDecayBranch &left, const MDecayBranch &right);

// Evaluate the charge and transfer factor of one ordered continuum pair vertex
double PairVertexFactor(const Param &param, const VertexParam &vertex, const MDecayBranch &left, const MDecayBranch &right, double upper_t, double lower_t, bool upper_ff, bool lower_ff);

// Evaluate one complete ordered continuum pair scalar factor
std::complex<double> PairExchange(const Param &param, const PairParam &entry, const VertexParam &vertex, const MDecayBranch &left, const MDecayBranch &right, double t_hat, double mass2, double upper_t,
                                  double lower_t);

// Evaluate the zero-secondary-production amplitude factor
double Veto(const VetoParam &veto, double central_mass);

}  // namespace gra::regge

#endif
