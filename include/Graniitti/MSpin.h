// Spin polarization and correlation functions
//
// (c) 2017-2022 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MSPIN_H
#define MSPIN_H

// C++
#include <complex>
#include <random>
#include <vector>

// Own
#include "Graniitti/MKinematics.h"
#include "Graniitti/MRandom.h"

namespace gra {
namespace spin {

// ---------
// Recursive

void postOrder(const MDecayBranch& root, std::vector<MDecayBranch>& stack, const std::string dir = "");
int  tree_height(const MDecayBranch& node);
void get_ijk_amplitude(const MDecayBranch& i, const MDecayBranch& j, const MDecayBranch& k,
    MMatrix<std::complex<double>>& amplitude, std::string& name);
void backward_traverse_amplitude(std::vector<MDecayBranch>& stack);
int  count_cmark(const std::string& s, const char marker);
std::vector<std::size_t> normal_order(const std::vector<MDecayBranch>& stack, std::size_t i, std::size_t j);

void TreeRecursion(MDecayBranch &branch, const gra::LORENTZSCALAR &lts,
    const std::string& mother, const std::string& dir);
MMatrix<std::complex<double>> CalculateFMatrix(const MDecayBranch &branch, const gra::LORENTZSCALAR &lts);
// ---------

std::complex<double> ProdAmp(const gra::LORENTZSCALAR &lts, gra::PARAM_RES &res);
std::complex<double> DecayAmp(gra::LORENTZSCALAR &lts, gra::PARAM_RES &res);
void GetRhoRotation(const gra::LORENTZSCALAR &lts, const std::string &FRAME, double &theta_R,
                    double &phi_R);

MMatrix<std::complex<double>> DMatrix(double J, double theta_mother, double phi_mother);
MMatrix<std::complex<double>> fMatrix(const MMatrix<std::complex<double>> &T, double J, double s1,
                                      double s2, double theta, double phi);

std::vector<double> SpinProjections(double J);
void InitTMatrix(gra::HELMatrix &hc, const gra::MParticle &p, const gra::MParticle &p1,
                 const gra::MParticle &p2);

// Spin-Statistics
bool BoseSymmetry(int l, int s);
bool FermiSymmetry(int l, int s);

// Spin algebra functions
double               ClebschGordan(double j1, double j2, double m1, double m2, double j, double m);
double               CGrules(double j1, double j2, double m1, double m2, double j, double m);
bool                 chalfint(double x);
bool                 cint(double x);
bool                 isequal(double x, double y);
double               W3j(double j1, double j2, double j3, double m1, double m2, double m3);
bool                 W3jrules(double j1, double j2, double j3, double m1, double m2, double m3);
double               TriangleCoeff(double j1, double j2, double j3);
std::complex<double> WignerD(double theta, double phi, double m, double mp, double J);
double               WignerSmalld(double theta, double m, double mp, double J);

// Density matrix functions
bool                          Positivity(const MMatrix<std::complex<double>> &rho, unsigned int J);
MMatrix<std::complex<double>> RandomRho(unsigned int J, bool parity, MRandom &rng);
double                        VonNeumannEntropy(const MMatrix<std::complex<double>> &rho);

}  // namespace spin
}  // namespace gra

#endif