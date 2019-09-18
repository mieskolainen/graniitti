// Form factors, structure functions, Regge trajectories etc. parametrizations
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MFORM_H
#define MFORM_H

// C++
#include <complex>
#include <map>
#include <vector>

// JSON library
#include "json.hpp"

// Keep this here to avoid backward/forward declaration problem
#include "Graniitti/M4Vec.h"
#include "Graniitti/MAux.h"
#include "Graniitti/MKinematics.h"
#include "Graniitti/MRandom.h"

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

void        PrintParam();
std::string GetHashString();
}

namespace PARAM_STRUCTURE {
extern std::string F2;
extern std::string EM;
}

// Flat (DEBUG) amplitude parameters
namespace PARAM_FLAT {
extern double b;
}

// Monopole wave functions and parameters
namespace PARAM_MONOPOLE {
extern int         En;        // Bound state energy level
extern double      M0;        // Monopole mass
extern double      Gamma0;    // Monopolium width
extern std::string coupling;  // Coupling scenarios
extern int         gn;        // Monopole charge n = 1,2,3,...

double EnergyMP(double n);
double GammaMP(double n, double alpha_g);
double PsiMP(double n);

void PrintParam(double sqrts, bool forceprint = false);

extern int _printcalls;
}

// N* excitation
namespace PARAM_NSTAR {
extern std::vector<double> rc;
}

namespace form {
// Read resonance
gra::PARAM_RES ReadResonance(const std::string &resparam_str, MRandom &rng);

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
double alpha_EM(double Q2);
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

// QED coupling alpha_EM at Q^2 = 0 is 1/137.035999679,
// gives 1/alpha_EM^2 = 0.085.
// Electric charge in natural units thus ~ 0.3
double e_EM(double Q2);
double e_EM();

// Breit-Wigner functions
double deltaBWxsec(double shat, double M0, double Gamma);
double deltaBWamp(double shat, double M0, double Gamma);

std::complex<double> CBW(const gra::LORENTZSCALAR &lts, const gra::PARAM_RES &resonance);
std::complex<double> CBW_FW(double m2, double M0, double Gamma);
std::complex<double> CBW_RW(double m2, double M0, double Gamma);
std::complex<double> CBW_BF(double m2, double M0, double Gamma, int J, double mA, double mB);
std::complex<double> CBW_JR(double m2, double M0, double Gamma, double J);

}  // form namespace
}  // gra namespace

#endif