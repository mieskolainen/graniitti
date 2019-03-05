// Gamma-Gamma -> Lepton Pair

//==========================================================================
// This file has been automatically generated for C++ Standalone
// MadGraph5_aMC@NLO v. 2.6.0, 2017-08-16
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef MG5_sm_yy_ww_H
#define MG5_sm_yy_ww_H

#include <cmath>
#include <complex>

#include "Graniitti/Amplitude/Parameters_sm.h"

using namespace std;

namespace MG5_sm_yy_ww {
void ixxxxx(double p[4], double fmass, int nhel, int nsf,
            std::complex<double> fi[6]);

double Sgn(double e, double f);

void txxxxx(double p[4], double tmass, int nhel, int nst,
            std::complex<double> fi[18]);

void sxxxxx(double p[4], int nss, std::complex<double> sc[3]);

void oxxxxx(double p[4], double fmass, int nhel, int nsf,
            std::complex<double> fo[6]);

void vxxxxx(double p[4], double vmass, int nhel, int nsv,
            std::complex<double> v[6]);

void VVV1_0(std::complex<double> V1[], std::complex<double> V2[],
            std::complex<double> V3[], std::complex<double> COUP,
            std::complex<double>& vertex);

void VVV1_2(std::complex<double> V1[], std::complex<double> V3[],
            std::complex<double> COUP, double M2, double W2,
            std::complex<double> V2[]);

void VVVV2_0(std::complex<double> V1[], std::complex<double> V2[],
             std::complex<double> V3[], std::complex<double> V4[],
             std::complex<double> COUP, std::complex<double>& vertex);

void VVV1_3(std::complex<double> V1[], std::complex<double> V2[],
            std::complex<double> COUP, double M3, double W3,
            std::complex<double> V3[]);
}

#endif
