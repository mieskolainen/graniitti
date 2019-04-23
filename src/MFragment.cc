// KISS fragmentation class
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <complex>
#include <random>
#include <valarray>
#include <vector>

// Own
#include "Graniitti/M4Vec.h"
#include "Graniitti/MAux.h"
#include "Graniitti/MForm.h"
#include "Graniitti/MFragment.h"
#include "Graniitti/MKinematics.h"
#include "Graniitti/MMath.h"
#include "Graniitti/MMatrix.h"
#include "Graniitti/MPDG.h"
#include "Graniitti/MRandom.h"

// Libraries
#include "rang.hpp"

using gra::aux::indices;
using gra::math::pow2;
using gra::math::msqrt;

using gra::PDG::mp;
using gra::PDG::mpi;

using gra::PDG::PDG_p;
using gra::PDG::PDG_pip;
using gra::PDG::PDG_pim;
using gra::PDG::PDG_pi0;
using gra::PDG::PDG_gamma;
using gra::PDG::PDG_n;

namespace gra {
// Return decay status
//
// Should update this function to draw exponential decay lengths,
// and based on a user defined ctau say 1 cm, decide if stable or not
//
//
void MFragment::GetDecayStatus(const std::vector<int> &pdgcode, std::vector<bool> &isstable) {
  isstable.resize(pdgcode.size(), true);
  for (const auto &i : indices(isstable)) {
    isstable[i] = (pdgcode[i] == 111 || pdgcode[i] == 311) ? false : true;
  }
}

// Sample from Hagedorn type parametrization
//
// This function is useful to generate distribution with
// exponential at low pt, powerlaw at high pt.
// Does not generate jet like topologies,
// so this can be used as a strict null model (H0) for soft like topologies.
//
// [slow function, performance should be upgraded]
//
// Parameters <=>
//  p0 = T/(q-1)
//   n = 1/(q-1)
//
void MFragment::ExpPowRND(double q, double T, double maxpt, const std::vector<double> &mass,
                          std::vector<double> &x, MRandom &rng) {
  const unsigned int MAXTRIAL = 1e4;  // Safety break
  const unsigned int Nbins    = 1e4;
  const double       ptstep   = maxpt / Nbins;

  // Calculate pt-bin values
  std::vector<double> ptval(Nbins);
  for (std::size_t i = 0; i < Nbins; ++i) { ptval[i] = i * ptstep; }

  // Random integer from [0,NBins-1]
  // C++11, thread_local is also static
  thread_local std::uniform_int_distribution<int> RANDI(0, Nbins - 1);

  // Calculate for pion, kaon, proton masses
  const std::vector<double>        fixmass = {PDG::mpi0, PDG::mpi, PDG::mK, PDG::mp};
  std::vector<std::vector<double>> dsdpt(fixmass.size(), std::vector<double>(Nbins, 0.0));
  std::vector<double>              maxval(fixmass.size(), 0.0);

  // Pre-evaluate dsigma/dpt |_y=0 (mid rapidity)
  for (const auto &k : indices(fixmass)) {
    const double m2 = pow2(fixmass[k]);
    for (std::size_t i = 0; i < Nbins; ++i) {
      const double mtval = msqrt(pow2(ptval[i]) + m2);

      // PDF
      dsdpt[k][i] = ptval[i] * mtval * std::pow(1.0 + (q - 1.0) * mtval / T, q / (1.0 - q));
      // Save maximum
      maxval[k] = (dsdpt[k][i] > maxval[k]) ? dsdpt[k][i] : maxval[k];
    }
  }

  // For each particle, generate pt
  x.resize(mass.size(), 0.0);
  for (const auto &p : indices(mass)) {
    // Find the best distribution
    unsigned int best    = 0;
    double       mindist = 1e32;
    for (std::size_t k = 0; k < fixmass.size(); ++k) {
      const double dist = std::abs(fixmass[k] - mass[p]);
      if (dist < mindist) {
        mindist = dist;
        best    = k;
      }
    }

    // Acceptance-Rejection
    unsigned int trials = 0;
    while (true) {
      const int BIN = RANDI(rng.rng);
      if (rng.U(0.0, 1.0) * maxval[best] < dsdpt[best][BIN]) {
        x[p] = ptval[BIN];
        break;
      }
      ++trials;
      if (trials > MAXTRIAL) {
        break;  // failed
      }
    }
  }
}

// "Retro" N-body fragmentation with tube (cylinder) phase space
// (approximate target distribution is flat over rapidity, pt from a given
// distribution)
//
// [REFERENCE: Jadach, Computer Physics Communications, 9 (1975) 297-304]
// [REFERENCE: UA5, http://cds.cern.ch/record/173907/files/198701320.pdf]
//
// This function contains a mixture of dynamics and kinematics, i.e.,
// is not a pure phase space and is suitable for soft fragmentation studies.
//
// Return: 1.0 for a valid fragmentation and -1.0 for a kinematically impossible

double MFragment::TubeFragment(const M4Vec &mother, double M0, const std::vector<double> &m,
                               std::vector<M4Vec> &p, double q, double T, double maxpt,
                               MRandom &rng) {
  // Number of re-trials in a case of failing kinematics
  const unsigned int MAXTRIAL = 30;
  unsigned int       trials   = 0;

  const unsigned int    N = m.size();
  std::valarray<double> m2(N);
  for (const auto &i : indices(m)) { m2[i] = pow2(m[i]); }
  // transverse momentum px, py
  std::valarray<double> px(N);
  std::valarray<double> py(N);
  std::valarray<double> pt2(N);

  // rapidity and transverse mass
  std::valarray<double> y(N);
  std::valarray<double> mt(N);

  while (true) {
    // Random sample px,py and initial rapidity
    std::vector<double> pt(N);
    ExpPowRND(q, T, maxpt, m, pt, rng);

    for (const auto &i : indices(pt)) {
      // Sample angle, px, py
      const double phi = rng.U(0, 2.0 * gra::math::PI);
      px[i]            = pt[i] * std::cos(phi);
      py[i]            = pt[i] * std::sin(phi);

      // Rapidity
      y[i] = rng.U(0.0, 1.0);
    }
    std::sort(begin(y), end(y));  // Rapidity ordering

    // ---------------------------------------------------------------
    // Scale rapidities	to span the full range [0,1] (min=0, max=N-1)
    y = (y - y[0]) / (y[N - 1] - y[0]);

    // Set transverse momentum with zero sum
    px  = px - px.sum() / N;
    py  = py - py.sum() / N;
    pt2 = px * px + py * py;

    // Transverse mass for all particles
    mt = sqrt(m2 + pt2);

    // ---------------------------------------------------------------
    // Solve rapidity scale factor 'alpha'
    double alpha = 0.0;
    if (!SolveAlpha(alpha, M0, m, mt, y)) {
      ++trials;
      if (trials > MAXTRIAL) {
        return -1.0;
      } else {
        continue;
      }
    }
    // Scale all rapidities
    y = alpha * y;

    // Set longitudinal momentum and energy using 'alpha'
    // get boost factor to +- 0
    const double Q                 = (mt * exp(y)).sum();
    y                              = y + std::log(M0 / Q);
    const std::valarray<double> pz = mt * sinh(y);
    const std::valarray<double> E  = mt * cosh(y);

    // ---------------------------------------------------------------
    // Set all particles 4-momenta and boost to the original frame
    p.resize(N);
    const int sign = 1;  // positive
    for (const auto &i : indices(p)) {
      p[i] = M4Vec(px[i], py[i], pz[i], E[i]);
      gra::kinematics::LorentzBoost(mother, M0, p[i], sign);
    }

    // Check EM-conservation
    M4Vec p_sum(0, 0, 0, 0);
    for (const auto &n : p) { p_sum += n; }
    bool             valid = gra::math::CheckEMC(p_sum - mother);

    // Check rapidity (floating points can fail after boost in forward)
    for (const auto &i : indices(p)) {
      if (std::isnan(p[i].Rap()) || std::isinf(p[i].Rap())) {
        valid = false;
        break;
      }
    }
    if (!valid) {
      ++trials;
      if (trials > MAXTRIAL) {
        return -1.0;
      } else {
        continue;
      }
    }
    break;  // Fragmentation succesful!
  }
  return 1.0;
}

// ---------------------------------------------------------------
// Solve rapidity scale factor via Newton iteration
// (alternative strategies are viable, too, this does not converge
// always with non-gaussian pt-distributions)
//
// a_1 = a_0 - f(a_0)/f'(a_1), find root a such that f(a) = 0
//   iteratively via
// a_{n+1} = a_n - f(a_n) / f'(a_n)
//
bool MFragment::SolveAlpha(double &alpha, double M0, const std::vector<double> &m,
                           const std::valarray<double> &mt, const std::valarray<double> &y) {
  const double       STOP_EPS = 1e-8;
  const unsigned int MAXITER  = 20;

  // Starting value
  const int    N        = m.size();
  const double C        = std::log(M0 * M0);
  alpha                 = C - std::log(m[0] * m[N - 1]);
  std::vector<double> E = {0, 0, 0, 0};

  unsigned int iter = 0;
  while (true) {
    const std::valarray<double> x = exp(alpha * y);

    E[0] = (mt * x).sum();
    E[1] = (mt / x).sum();
    E[2] = (mt * y * x).sum();  // Derivative d/dy
    E[3] = (mt * y / x).sum();  // Derivative d/dy

    // Iterate the solution, d/dx ln(x) = 1/x
    const double DY = E[0] * E[1] * (C - std::log(E[0] * E[1])) / (E[0] * E[3] - E[1] * E[2]);
    alpha -= DY;

    if (std::abs(N * DY / alpha) < STOP_EPS) { return true; }
    if (std::isnan(alpha) || (++iter) > MAXITER) { return false; }
  }
}

// N* decay table [set manually according to experimental data]
// M0 is the N* mass
void MFragment::NstarDecayTable(double M0, std::vector<int> &pdgcode, MRandom &rng) {
  int decaymode = 0;

  if (M0 < (PDG::mp + 2 * PDG::mpi)) {
    // Only 2-body decay possible, mass below 3-body threshold
    decaymode = 0;
  } else {  // 2- or 3-body decay possible

    // C++11, thread_local is also static
    thread_local std::discrete_distribution<> d(
        {0.60, 0.40});       // 2->body / 3-body branching ratios from PDG
    decaymode = d(rng.rng);  // Draw random
  }

  // ----------------------------------------------------------------------
  // N*(J^P = 1/2^+) Decay Parameters from PDG
  // Only major decays implemented

  std::vector<int> decayids;

  // 2-Body channel
  if (decaymode == 0) {
    // Subchannels
    // C++11, thread_local is also static
    thread_local std::discrete_distribution<> subd(
        {2.0 / 3.0, 1.0 / 3.0});  // Branching ratios from Clebsch-Gordan
    const int channel = subd(rng.rng);

    // PDG-ID of subchannels
    const std::vector<std::vector<int>> ID = {{PDG_n, PDG_pip},   // neutron & pi+
                                              {PDG_p, PDG_pi0}};  // proton  & pi0

    // Choose the decay channel
    pdgcode = ID[channel];
  }

  // 3-body channel
  if (decaymode == 1) {
    // Subchannels
    // C++11, thread_local is also static
    thread_local std::discrete_distribution<> subd(
        {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0});  // Braching ratios ansatz
    const int channel = subd(rng.rng);

    // PDG-ID of subchannels
    const std::vector<std::vector<int>> ID = {
        {PDG_n, PDG_pip, PDG_pi0}, {PDG_p, PDG_pip, PDG_pim}, {PDG_p, PDG_pi0, PDG_pi0}};

    // Perhaps to add
    //{PDG_delta0, PDG_pip, PDG_pi0},     // delta0 & pi+ & pi0
    //{PDG_deltap, PDG_pip, PDG_pim},     // delta+ & pi+ & pi-
    //{PDG_deltap, PDG_pi0, PDG_pi0}};    // delta+ & pi0 & pi0

    // Draw the decay channel
    pdgcode = ID[channel];
  }
}

// Return excited forward proton masses
void MFragment::GetForwardMass(double &mass1, double &mass2, bool &excite1, bool &excite2,
                               unsigned int excite, MRandom &random) {
  if (excite == 0) {  // Fully elastic
    mass1   = PDG::mp;
    mass2   = PDG::mp;
    excite1 = false;
    excite2 = false;
  } else if (excite == 1) {  // Single excitation
    if (random.U(0, 1) < 0.5) {
      mass1 = PDG::mp;
      GetSingleForwardMass(mass2, random);
      excite1 = false;
      excite2 = true;
    } else {
      GetSingleForwardMass(mass1, random);
      mass2   = PDG::mp;
      excite1 = true;
      excite2 = false;
    }
  } else if (excite == 2) {  // Double excitation
    GetSingleForwardMass(mass1, random);
    GetSingleForwardMass(mass2, random);
    excite1 = true;
    excite2 = true;
  }
}

// Return excited forward system mass
void MFragment::GetSingleForwardMass(double &mass, MRandom &random) {
  // Good-Walker resonance excitation probabilities
  // C++11, thread_local is also static
  thread_local std::discrete_distribution<> d(
      {PARAM_NSTAR::rc[0], PARAM_NSTAR::rc[1], PARAM_NSTAR::rc[2]});

  const int state = d(random.rng);  // Draw random

  // Excited proton states
  double M0    = 0.0;
  double width = 0.0;

  if (state == 0) {
    M0    = 1.440;  // N*
    width = 0.325;
  } else if (state == 1) {
    M0    = 1.680;  // N**
    width = 0.140;
  } else if (state == 2) {
    M0    = 2.190;  // N***
    width = 0.450;
  }

  // Draw the excited state mass
  const double safe_margin = 0.01;
  double       mX          = 0;
  while (true) {
    mX = random.RelativisticBWRandom(M0, width, 1e6);
    if (mX > (PDG::mp + PDG::mpi + safe_margin)) {  // && mX < lts.sqrt_s) {
      mass = mX;
      return;
    }
  }
}

// Simple statistical toy particle pick-up, nothing more
//
int MFragment::PickParticles(double M, unsigned int N, int B, int S, int Q,
                             std::vector<double> &mass, std::vector<int> &pdgcode, const MPDG &PDG,
                             MRandom &rng) {
  const unsigned int MAXTRIAL = 1e5;

  // C++11, thread_local is also static
  thread_local std::uniform_int_distribution<int> RANDI3(0, 2);

  // Pion, Kaon, Proton ratios [MEASURED]
  const double              denom3 = (1 + 0.1 + 0.06);
  const std::vector<double> ratio3 = {1.0 / denom3, 0.1 / denom3, 0.06 / denom3};

  const std::vector<int> charged_pdg = {211, 321, 2212};  // pi+-, K+-, (anti)proton
  const std::vector<int> charged_B   = {0, 0, 1};
  const std::vector<int> charged_S   = {0, 1, 0};

  const std::vector<int> neutral_pdg = {111, 311, 2112};  // pi0, K0, neutron
  const std::vector<int> neutral_B   = {0, 0, 1};
  const std::vector<int> neutral_S   = {0, 1, 0};

  // -----------------------------------------------------------------

  // Now resize!
  mass.resize(N, 0);
  pdgcode.resize(N, 0);
  std::vector<int> Qcharges(N, 0);
  std::vector<int> Bcharges(N, 0);
  std::vector<int> Scharges(N, 0);

  unsigned int trials = 0;
  const double QProb  = 2.0 / 3.0;
  const double DOUBLE = 0.9;

  while (true) {
    unsigned int i = 0;

    do {  // Hadron picking

      // Charged
      if (rng.U(0, 1) < QProb) {
        // Charged pair
        if (rng.U(0, 1) < DOUBLE && i <= N - 2) {
          // Sample particle flavour
          while (true) {
            const int bin = RANDI3(rng.rng);

            if (rng.U(0, 1) < ratio3[bin]) {
              Qcharges[i] = 1;
              Bcharges[i] = charged_B[bin];
              Scharges[i] = charged_S[bin];
              pdgcode[i]  = charged_pdg[bin];
              ++i;  // Next slot

              Qcharges[i] = -1;
              Bcharges[i] = (-1) * charged_B[bin];
              Scharges[i] = (-1) * charged_S[bin];
              pdgcode[i]  = (-1) * charged_pdg[bin];
              ++i;
              break;
            }
          }

          // Single
        } else {
          while (true) {
            const int bin = RANDI3(rng.rng);

            if (rng.U(0, 1) < ratio3[bin]) {
              int sign = 0;
              sign     = (rng.U(0, 1) < 0.5) ? 1 : -1;

              Qcharges[i] = sign;
              Bcharges[i] = sign * charged_B[bin];
              Scharges[i] = sign * charged_S[bin];
              pdgcode[i]  = sign * charged_pdg[bin];
              ++i;
              break;
            }
          }
        }

        // Neutral
      } else {
        while (true) {
          // Sample particle flavour
          const int bin = RANDI3(rng.rng);
          if (rng.U(0, 1) < ratio3[bin]) {
            Qcharges[i] = 0;
            Bcharges[i] = neutral_B[bin];
            Scharges[i] = neutral_S[bin];
            pdgcode[i]  = neutral_pdg[bin];
            ++i;
            break;
          }
        }
      }
    } while (i < N);  // Picking loop

    // Sum charges
    const int Q_sum = std::accumulate(Qcharges.begin(), Qcharges.end(), 0);
    const int B_sum = std::accumulate(Bcharges.begin(), Bcharges.end(), 0);
    const int S_sum = std::accumulate(Scharges.begin(), Scharges.end(), 0);

    // Get corresponding masses
    for (std::size_t i = 0; i < N; ++i) {
      // printf("PDGCODE[i=%d] = %d [N=%d] \n", i, pdgcode[i], N);
      mass[i] = PDG.FindByPDG(pdgcode[i]).mass;
    }
    const double M_sum = std::accumulate(mass.begin(), mass.end(), 0.0);

    // Check charge and mass threshold
    bool Q_check = (Q_sum == Q) ? true : false;
    bool B_check = (B_sum == B) ? true : false;
    bool S_check = (S_sum == S) ? true : false;
    bool M_check = (M_sum < M) ? true : false;

    if (Q_check && B_check && S_check && M_check) {
      break;  // we are ok!
    } else {
      // printf("[PICK] Q_sum = %d, B_sum = %d, M_sum = %0.1f, M = %0.1f, N = %d
      // \n", Q_sum, B_sum, M_sum, M, N);
      ++trials;
      if (trials > MAXTRIAL) { return 1; }
    }
  }  // Out loop

  return 0;
}

}  // gra namespace ends
