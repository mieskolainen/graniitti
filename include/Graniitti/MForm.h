// Form factors, structure functions, Regge trajectories etc. parametrizations
//
// (c) 2017-2020 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MFORM_H
#define MFORM_H

// C++
#include <complex>
#include <map>
#include <vector>

// Own
#include "Graniitti/M4Vec.h"
#include "Graniitti/MAux.h"
#include "Graniitti/MKinematics.h"
#include "Graniitti/MParticle.h"
#include "Graniitti/MRandom.h"
#include "Graniitti/MResonance.h"

// Libraries
#include "json.hpp"

namespace gra {

// Model parameters
namespace PARAM_SOFT {
// Pomeron trajectory
extern double DELTA_P;
extern double ALPHA_P;

// Couplings
extern double gN_P;
extern double gN_O;

extern double g3P;
extern double gamma;

// Proton form factor
extern double fc1;
extern double fc2;
extern double fc3;

extern bool ODDERON_ON;

void        PrintParam();
std::string GetHashString();

void        ReadParameters(const std::string& modelfile);
extern bool initialized;
}  // namespace PARAM_SOFT

namespace PARAM_STRUCTURE {
extern std::string F2;
extern std::string EM;
extern std::string QED_alpha;

void        ReadParameters(const std::string& modelfile);
extern bool initialized;
}  // namespace PARAM_STRUCTURE

// Flat (DEBUG) amplitude parameters
namespace PARAM_FLAT {
extern double b;

void        ReadParameters(const std::string& modelfile);
extern bool initialized;
}  // namespace PARAM_FLAT

// Forward proton excitation
namespace PARAM_NSTAR {
extern std::string         fragment;
extern std::vector<double> rc;

void        ReadParameters(const std::string& modelfile);
extern bool initialized;
}  // namespace PARAM_NSTAR

namespace form {

// Read resonance
gra::PARAM_RES ReadResonance(const std::string& resparam_str, MRandom& rng);

// Regge signature
std::complex<double> ReggeEta(double alpha_t, double sigma);
std::complex<double> ReggeEtaLinear(double t, double alpha_t0, double ap, double sigma);

// Proton form factor, elastic and inelastic
double S3F(double t);
double S3FINEL(double t, double M2);

// Proton structure functions
double F2xQ2(double xi, double Q2);
double F1xQ2(double xi, double Q2);

// Pomeron trajectory
double S3PomAlpha(double t);

// Pion loop insert
double S3HPL(double tau, double t);

// t-integrated collinear EPA flux
double DZFlux(double x);

// Coherent gamma flux
double F1(double Q2);
double F2(double Q2);
double G_M(double Q2);
double G_E(double Q2);

double G_M_DIPOLE(double Q2);
double G_E_DIPOLE(double Q2);

double G_M_KELLY(double Q2);
double G_E_KELLY(double Q2);

double CohFlux(double x, double t, double pt);
double IncohFlux(double x, double t, double pt, double M2);

double mu_ratio();

// Breit-Wigner functions
double deltaBWxsec(double shat, double M0, double Gamma);
double deltaBWamp(double shat, double M0, double Gamma);

std::complex<double> CBW(const gra::LORENTZSCALAR& lts, const gra::PARAM_RES& resonance);
std::complex<double> CBW_FW(double m2, double M0, double Gamma);
std::complex<double> CBW_RW(double m2, double M0, double Gamma);
std::complex<double> CBW_BF(double m2, double M0, double Gamma, int J, double mA, double mB);
std::complex<double> CBW_JR(double m2, double M0, double Gamma, double J);

}  // namespace form
}  // namespace gra

#endif